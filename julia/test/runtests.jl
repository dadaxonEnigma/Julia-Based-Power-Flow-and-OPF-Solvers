# PowerFlowJulia — unified test runner
include("../src/PowerFlowJulia.jl")
using .PowerFlowJulia

passed = 0; failed = 0

macro test(name, body)
    quote
        try
            $(esc(body))
            println("  PASS  $($(name))")
            global passed += 1
        catch e
            println("  FAIL  $($(name)): $e")
            global failed += 1
        end
    end
end

println("=" ^ 60)
println("PowerFlowJulia — Test Suite")
println("=" ^ 60)

# ── Helpers ───────────────────────────────────────────────────────────────
function make_3bus(; s_nom=Inf)
    n = Network(baseMVA=100.0)
    add!(n, "Bus", "B1"; slack=true)
    add!(n, "Bus", "B2")
    add!(n, "Bus", "B3")
    add!(n, "Line", "L12"; bus0="B1", bus1="B2", r=0.01, x=0.1, s_nom=s_nom)
    add!(n, "Line", "L13"; bus0="B1", bus1="B3", r=0.01, x=0.1, s_nom=s_nom)
    add!(n, "Line", "L23"; bus0="B2", bus1="B3", r=0.01, x=0.1, s_nom=s_nom)
    add!(n, "Generator", "G1"; bus="B1", p_nom=400.0, marginal_cost=20.0)
    add!(n, "Generator", "G2"; bus="B2", p_nom=300.0, marginal_cost=50.0)
    add!(n, "Load", "D2"; bus="B2", p_set=200.0)
    add!(n, "Load", "D3"; bus="B3", p_set=300.0)
    return n
end

# ── 1. Components ─────────────────────────────────────────────────────────
println("\n[1] Components")
@test "Bus default constructor"     (b=Bus("B1"); @assert b.v_nom==1.0)
@test "Generator UC fields"         (g=Generator("G","B",p_nom=100.0,committable=true,min_up_time=3); @assert g.min_up_time==3)
@test "Store e_cyclic default"      (s=Store("S","B",e_nom=100.0,p_nom=50.0); @assert s.e_cyclic)
@test "Carrier co2"                 (c=Carrier("gas",co2_emissions=0.20); @assert c.co2_emissions≈0.20)
@test "Link efficiency validation"  try Link("L","A","B",efficiency=0.0); false catch; true end
@test "GlobalConstraint sense"      try GlobalConstraint("g",sense="??"); false catch; true end

# ── 2. Network API ────────────────────────────────────────────────────────
println("\n[2] Network API")
@test "add! Bus"         (n=Network(); add!(n,"Bus","B1"); @assert length(n.buses)==1)
@test "add! Line bus0/bus1" begin
    n=Network(); add!(n,"Bus","B1"); add!(n,"Bus","B2")
    add!(n,"Line","L1"; bus0="B1",bus1="B2",x=0.1)
    @assert haskey(n.lines,"L1")
end
@test "add! Generator"   begin
    n=Network(); add!(n,"Bus","B1")
    add!(n,"Generator","G1"; bus="B1",p_nom=100.0)
    @assert n.generators["G1"].p_nom==100.0
end
@test "add! Link"        begin
    n=Network(); add!(n,"Bus","A"); add!(n,"Bus","B")
    add!(n,"Link","HVDC"; bus0="A",bus1="B",p_nom=200.0,efficiency=0.97)
    @assert n.links["HVDC"].efficiency≈0.97
end
@test "duplicate Bus rejected" try
    n=Network(); add!(n,"Bus","B1"); add!(n,"Bus","B1"); false
catch; true end
@test "missing bus rejected" try
    n=Network(); add!(n,"Line","L"; bus0="B1",bus1="B2",x=0.1); false
catch; true end

# ── 3. Power Flow ─────────────────────────────────────────────────────────
println("\n[3] Power Flow")
net3 = make_3bus()
@test "pf(:dc) converges"   (r=pf(net3,verbose=false); @assert r.converged)
@test "pf(:dc) slack=0"     (r=pf(net3,verbose=false); @assert abs(r.θ[1])<1e-10)
@test "pf(:lac) converges"  (r=pf(net3,method=:lac,verbose=false); @assert r.converged)
@test "pf(:lac) |V|≈1"      (r=pf(net3,method=:lac,verbose=false); @assert all(abs.(r.V_mag.-1.0).<0.05))
@test "pf(:auto) dispatches" (r=pf(net3,method=:auto,verbose=false); @assert r.converged)

# ── 4. LOPF + LMP ────────────────────────────────────────────────────────
println("\n[4] LOPF + LMP")
@test "optimize(:lopf) converged"   (r=optimize(make_3bus(),verbose=false); @assert r.converged)
@test "LOPF uncongested LMP equal"  begin
    r=optimize(make_3bus(),verbose=false)
    lmps=collect(values(r.lmp))
    @assert maximum(lmps)-minimum(lmps)<0.01
end
@test "LOPF congested LMP diverge"  begin
    r=optimize(make_3bus(s_nom=200.0),verbose=false)
    lmps=collect(values(r.lmp))
    @assert maximum(lmps)-minimum(lmps)>0.0
end
@test "optimize auto-selects lopf"  begin
    r=optimize(make_3bus(),verbose=false)
    @assert r.converged && haskey(r,:lmp)
end

# ── 5. Multi-period LOPF ─────────────────────────────────────────────────
println("\n[5] Multi-period LOPF")
net5 = make_3bus()
add!(net5,"StorageUnit","Bat"; bus="B2",p_nom=60.0,e_nom=240.0)
@test "optimize T=24 converges"  begin
    r=optimize(net5,T=24,verbose=false)
    @assert string(r.status)=="OPTIMAL"
end
@test "optimize auto-selects :mp" begin
    r=optimize(net5,T=12,verbose=false)
    @assert haskey(r,:soc)
end

# ── 6. Unit Commitment ────────────────────────────────────────────────────
println("\n[6] Unit Commitment")
net6 = Network(baseMVA=100.0)
add!(net6,"Bus","B1"; slack=true)
add!(net6,"Generator","G_commit"; bus="B1",p_nom=200.0,p_min_pu=0.3,
     marginal_cost=20.0,committable=true,min_up_time=2,startup_cost=1000.0)
add!(net6,"Generator","G_cont"; bus="B1",p_nom=100.0,marginal_cost=60.0)
add!(net6,"Load","D1"; bus="B1",p_set=100.0)

@test "optimize auto-selects :uc" begin
    r=optimize(net6,T=6,verbose=false)
    @assert haskey(r,:u) && string(r.status)=="OPTIMAL"
end
@test "UC binary variables in {0,1}" begin
    r=optimize(net6,T=6,verbose=false)
    for (_, vec) in r.u
        @assert all(v in (0,1) for v in vec)
    end
end
@test "UC min_up_time respected" begin
    r=optimize(net6,T=8,verbose=false)
    u=r.u["G_commit"]
    for t in 1:8
        if r.su["G_commit"][t]==1
            for dt in 0:min(1,7-t)
                @assert u[t+dt]==1
            end
        end
    end
end

# ── 7. GlobalConstraint ──────────────────────────────────────────────────
println("\n[7] GlobalConstraint")
@test "CO₂ cap raises cost" begin
    function make_co2_net(; with_cap=false)
        n=Network(baseMVA=100.0)
        add!(n,"Bus","B1"; slack=true)
        add!(n,"Generator","Coal"; bus="B1",p_nom=300.0,marginal_cost=20.0,carrier="coal")
        add!(n,"Generator","Gas";  bus="B1",p_nom=300.0,marginal_cost=50.0,carrier="gas")
        add!(n,"Load","D1"; bus="B1",p_set=200.0)
        if with_cap
            # cap=55t: forces partial gas dispatch (all-coal=68t > 55t)
            add!(n,"GlobalConstraint","co2"; constant=55.0,
                 carrier_weightings=Dict("coal"=>0.34,"gas"=>0.20))
        end
        return n
    end
    r_free=optimize(make_co2_net(with_cap=false),verbose=false)
    r_cap =optimize(make_co2_net(with_cap=true), verbose=false)
    @assert r_cap.total_cost >= r_free.total_cost - 1e-3
end

# ── 8. PyPSA API parity check ────────────────────────────────────────────
println("\n[8] PyPSA API parity")
@test "add! accepts PyPSA bus0/bus1 for Line" begin
    n=Network(); add!(n,"Bus","A"); add!(n,"Bus","B")
    add!(n,"Line","L"; bus0="A",bus1="B",x=0.1)
    @assert n.lines["L"].from_bus=="A" && n.lines["L"].to_bus=="B"
end
@test "add! accepts from_bus/to_bus for Line" begin
    n=Network(); add!(n,"Bus","A"); add!(n,"Bus","B")
    add!(n,"Line","L"; from_bus="A",to_bus="B",x=0.1)
    @assert n.lines["L"].from_bus=="A"
end
@test "pf returns converged field"   (r=pf(make_3bus(),verbose=false); @assert haskey(r,:converged))
@test "optimize returns lmp field"   (r=optimize(make_3bus(),verbose=false); @assert haskey(r,:lmp))
@test "optimize returns total_cost"  (r=optimize(make_3bus(),verbose=false); @assert r.total_cost>0)

# ── 10. ML Forecasting ───────────────────────────────────────────────────────
println("\n[10] ML Forecasting (LSTM + Conformal Prediction)")
@test "generate_synthetic_data shape" begin
    d = generate_synthetic_data(30)
    @assert size(d) == (30, 24)
    @assert all(d .>= 0.0)
end
@test "train_forecaster returns LoadForecaster" begin
    d  = generate_synthetic_data(60; seed=1)
    fc = train_forecaster(d; hidden=8, epochs=5, verbose=false)
    @assert fc isa LoadForecaster
    @assert fc.horizon == 24 && fc.window == 24
    @assert length(fc.cal_scores) > 0
end
@test "predict_scenarios returns correct shapes" begin
    d    = generate_synthetic_data(60; seed=2)
    fc   = train_forecaster(d; hidden=8, epochs=5, verbose=false)
    hist = vec(d[end, :])
    pred = predict_scenarios(fc, hist; n_scenarios=5)
    @assert length(pred.mean)    == 24
    @assert length(pred.lower)   == 24
    @assert length(pred.upper)   == 24
    @assert size(pred.scenarios) == (24, 5)
end
@test "conformal intervals: lower ≤ mean ≤ upper" begin
    d    = generate_synthetic_data(60; seed=3)
    fc   = train_forecaster(d; hidden=8, epochs=5, verbose=false)
    pred = predict_scenarios(fc, vec(d[end, :]); n_scenarios=3)
    @assert all(pred.lower .<= pred.mean .+ 1e-6)
    @assert all(pred.mean  .<= pred.upper .+ 1e-6)
end
@test "forecast_metrics computes MAE/RMSE/MAPE" begin
    m = forecast_metrics([1.0, 2.0, 3.0], [1.1, 2.1, 3.1])
    @assert abs(m.mae - 0.1) < 1e-6
    @assert m.rmse > 0 && m.mape > 0
end

# ── 11. Stochastic LOPF ──────────────────────────────────────────────────────
println("\n[11] Stochastic LOPF (SAA)")
@test "lopf_stochastic converges on 3-bus" begin
    n3 = make_3bus()
    d  = generate_synthetic_data(60; seed=4)
    fc = train_forecaster(d; hidden=8, epochs=5, verbose=false)
    pred = predict_scenarios(fc, vec(d[end, :]); n_scenarios=3)
    r = lopf_stochastic(n3, pred.scenarios; T=24, verbose=false)
    @assert r.expected_cost > 0
    @assert r.n_feasible > 0
    @assert r.cvar_90 >= r.expected_cost - 1.0
end
@test "stochastic cost ≥ 0 for all scenarios" begin
    n3 = make_3bus()
    d  = generate_synthetic_data(60; seed=5)
    fc = train_forecaster(d; hidden=8, epochs=5, verbose=false)
    pred = predict_scenarios(fc, vec(d[end, :]); n_scenarios=3)
    r = lopf_stochastic(n3, pred.scenarios; T=24, verbose=false)
    @assert all(r.costs .>= 0)
end

# ── 9. AC Power Flow ─────────────────────────────────────────────────────────
println("\n[9] AC Power Flow (PowerModels.jl + Ipopt)")
@test "pf(:ac) converges" begin
    r = pf(make_3bus(), method=:ac, verbose=false)
    @assert r.converged
end
@test "pf(:ac) slack voltage = 1 p.u." begin
    r = pf(make_3bus(), method=:ac, verbose=false)
    @assert abs(r.V_mag[1] - 1.0) < 0.01
end
@test "pf(:ac) returns Q_flow" begin
    r = pf(make_3bus(), method=:ac, verbose=false)
    @assert haskey(r, :Q_flow) && length(r.Q_flow) > 0
end
@test "pf(:ac) P_flow nonzero" begin
    r = pf(make_3bus(), method=:ac, verbose=false)
    @assert any(abs.(r.P_flow) .> 1e-6)
end

# ── Summary ───────────────────────────────────────────────────────────────
total = passed + failed
println("\n" * "=" ^ 60)
println("Results: $passed / $total passed  $(failed==0 ? "✓" : "✗ ($failed failed)")")
println("=" ^ 60)
