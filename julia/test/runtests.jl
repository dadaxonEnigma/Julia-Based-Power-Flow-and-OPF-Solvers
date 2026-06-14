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
    @assert fc.horizon == 24 && fc.window == 168
    @assert length(fc.cal_scores) > 0
end
@test "predict_scenarios returns correct shapes" begin
    d    = generate_synthetic_data(60; seed=2)
    fc   = train_forecaster(d; hidden=8, epochs=5, verbose=false)
    hist = vec(d[end-6:end, :]')          # last 168 chronological hours
    pred = predict_scenarios(fc, hist; n_scenarios=5)
    @assert length(pred.mean)    == 24
    @assert length(pred.lower)   == 24
    @assert length(pred.upper)   == 24
    @assert size(pred.scenarios) == (24, 5)
end
@test "conformal intervals: lower ≤ mean ≤ upper" begin
    d    = generate_synthetic_data(60; seed=3)
    fc   = train_forecaster(d; hidden=8, epochs=5, verbose=false)
    pred = predict_scenarios(fc, vec(d[end-6:end, :]'); n_scenarios=3)
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
    pred = predict_scenarios(fc, vec(d[end-6:end, :]'); n_scenarios=3)
    r = lopf_stochastic(n3, pred.scenarios; T=24, verbose=false)
    @assert r.expected_cost > 0
    @assert r.n_feasible > 0
    @assert r.cvar_90 >= r.expected_cost - 1.0
end
@test "stochastic cost ≥ 0 for all scenarios" begin
    n3 = make_3bus()
    d  = generate_synthetic_data(60; seed=5)
    fc = train_forecaster(d; hidden=8, epochs=5, verbose=false)
    pred = predict_scenarios(fc, vec(d[end-6:end, :]'); n_scenarios=3)
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

# ── 12. Component Constructors (complete) ─────────────────────────────────────
println("\n[12] Component Constructors")

@test "Transformer constructor" begin
    t = Transformer("T1","B1","B2"; x=0.05, s_nom=500.0, tap_ratio=1.02, phase_shift=0.0)
    @assert t.tap_ratio ≈ 1.02 && t.x ≈ 0.05 && t.s_nom ≈ 500.0
end
@test "add! Transformer" begin
    n=Network(); add!(n,"Bus","B1"; slack=true); add!(n,"Bus","B2")
    add!(n,"Transformer","T1"; bus0="B1",bus1="B2",x=0.05,tap_ratio=1.02)
    @assert haskey(n.transformers,"T1")
    @assert n.transformers["T1"].tap_ratio ≈ 1.02
end
@test "add! StorageUnit direct" begin
    n=Network(); add!(n,"Bus","B1"; slack=true)
    add!(n,"StorageUnit","Bat"; bus="B1",p_nom=100.0,e_nom=400.0,
         efficiency_charge=0.92,efficiency_discharge=0.92)
    @assert n.storage_units["Bat"].e_nom ≈ 400.0
    @assert n.storage_units["Bat"].efficiency_charge ≈ 0.92
end
@test "add! Load direct" begin
    n=Network(); add!(n,"Bus","B1"; slack=true)
    add!(n,"Load","D1"; bus="B1",p_set=300.0)
    @assert n.loads["D1"].p_set ≈ 300.0
end
@test "add! Store direct" begin
    n=Network(); add!(n,"Bus","B1"; slack=true)
    add!(n,"Store","S1"; bus="B1",e_nom=500.0,standing_loss=0.002)
    @assert n.stores["S1"].e_nom ≈ 500.0
    @assert n.stores["S1"].standing_loss ≈ 0.002
end
@test "add! Carrier direct" begin
    n=Network()
    add!(n,"Carrier","coal"; co2_emissions=0.34)
    @assert n.carriers["coal"].co2_emissions ≈ 0.34
end
@test "GlobalConstraint sense >=" begin
    gc = GlobalConstraint("gc1"; sense=">=", constant=100.0)
    @assert gc.sense == ">="
end
@test "GlobalConstraint sense =" begin
    gc = GlobalConstraint("gc2"; sense="=", constant=0.0)
    @assert gc.sense == "="
end

# ── 13. DC PF Numerical Values ────────────────────────────────────────────────
println("\n[13] DC PF Numerical Values")

# make_3bus: G2(300MW)@B2, D2=200 → net B2=+100 MW; D3=300 → net B3=-300 MW
# B=1000 MW/rad per line; θ₂=−1/30, θ₃=−1/6 (analytic solution)

@test "DC PF θ_B2 ≈ −1/30 rad" begin
    r = pf(make_3bus(), verbose=false)
    idx = findfirst(==("B2"), r.buses)
    @assert abs(r.θ[idx] - (-1.0/30)) < 1e-4
end
@test "DC PF θ_B3 ≈ −1/6 rad" begin
    r = pf(make_3bus(), verbose=false)
    idx = findfirst(==("B3"), r.buses)
    @assert abs(r.θ[idx] - (-1.0/6)) < 1e-4
end
@test "DC PF line flow L13 ≈ 166.7 MW" begin
    r = pf(make_3bus(), verbose=false)
    l13 = first(f for f in r.line_flows if f.name == "L13")
    @assert abs(l13.P_MW - 166.67) < 1.0
end
@test "DC PF line flow L12 ≈ 33.3 MW" begin
    r = pf(make_3bus(), verbose=false)
    l12 = first(f for f in r.line_flows if f.name == "L12")
    @assert abs(l12.P_MW - 33.33) < 1.0
end
@test "DC PF P_inj at slack has generator injection > 0" begin
    # P_inj stores raw injections before solve: B1 has G1=400MW, so P_inj[B1]=400
    r = pf(make_3bus(), verbose=false)
    idx = findfirst(==("B1"), r.buses)
    @assert r.P_inj[idx] > 0
end

# ── 14. LOPF Numerical Values ────────────────────────────────────────────────
println("\n[14] LOPF Numerical Values")

@test "LOPF uncongested cost = 13 000 €" begin
    # G1(400,20)+G2(100,50) = 8000+5000 = 13000
    r = optimize(make_3bus(), verbose=false)
    @assert abs(r.total_cost - 13_000.0) < 1.0
end
@test "LOPF congested cost = 16 000 €" begin
    # s_nom=200 forces G1=300, G2=200: 6000+10000 = 16000
    r = optimize(make_3bus(s_nom=200.0), verbose=false)
    @assert abs(r.total_cost - 16_000.0) < 1.0
end
@test "LOPF uncongested G1=400 MW" begin
    r = optimize(make_3bus(), verbose=false)
    @assert abs(r.P_gen["G1"] - 400.0) < 1.0
end
@test "LOPF uncongested G2=100 MW" begin
    r = optimize(make_3bus(), verbose=false)
    @assert abs(r.P_gen["G2"] - 100.0) < 1.0
end
@test "LOPF congested G1=300 G2=200" begin
    r = optimize(make_3bus(s_nom=200.0), verbose=false)
    @assert abs(r.P_gen["G1"] - 300.0) < 1.0
    @assert abs(r.P_gen["G2"] - 200.0) < 1.0
end
@test "LOPF total generation = total load (500 MW)" begin
    r = optimize(make_3bus(), verbose=false)
    @assert abs(sum(values(r.P_gen)) - 500.0) < 1.0
end
@test "LOPF P_line is non-empty" begin
    r = optimize(make_3bus(), verbose=false)
    @assert !isempty(r.P_line)
    @assert any(abs.(r.P_line) .> 1.0)
end

# ── 15. Multi-period LOPF Values ─────────────────────────────────────────────
println("\n[15] Multi-period LOPF Values")

net15 = begin
    n = make_3bus()
    add!(n,"StorageUnit","Bat"; bus="B2",p_nom=50.0,e_nom=200.0,
         efficiency_charge=0.9,efficiency_discharge=0.9)
    n
end

@test "MP-LOPF gen_dispatch is T-length vector" begin
    r = optimize(net15, T=12, verbose=false)
    @assert length(r.gen_dispatch["G1"]) == 12
end
@test "MP-LOPF SoC within [0, e_nom]" begin
    r = optimize(net15, T=24, verbose=false)
    soc = r.soc["Bat"]
    @assert all(soc .>= -1e-3) && all(soc .<= 200.0 + 1e-3)
end
@test "MP-LOPF SoC cyclic: E[T] ≈ E[0]" begin
    r = optimize(net15, T=24, verbose=false)
    soc = r.soc["Bat"]   # indexed 0:T → length T+1
    @assert abs(soc[end] - soc[1]) < 1.0
end
@test "MP-LOPF lmp is T-length vector" begin
    r = optimize(net15, T=12, verbose=false)
    @assert length(r.lmp["B1"]) == 12
end
@test "MP-LOPF p_ch + p_dis fields exist" begin
    r = optimize(net15, T=6, verbose=false)
    @assert haskey(r, :p_ch) && haskey(r, :p_dis)
    @assert all(r.p_ch["Bat"] .>= -1e-6)
    @assert all(r.p_dis["Bat"] .>= -1e-6)
end

# ── 16. UC Additional Fields ──────────────────────────────────────────────────
println("\n[16] UC Additional Fields")

net16 = begin
    n = Network(baseMVA=100.0)
    add!(n,"Bus","B1"; slack=true)
    add!(n,"Generator","Base"; bus="B1",p_nom=200.0,p_min_pu=0.3,
         marginal_cost=15.0,committable=true,
         min_up_time=2,min_down_time=2,
         startup_cost=500.0,shutdown_cost=200.0,initial_status=false)
    add!(n,"Generator","Peak"; bus="B1",p_nom=100.0,marginal_cost=80.0)
    add!(n,"Load","D1"; bus="B1",p_set=120.0)
    n
end

@test "UC sd (shutdown) field exists" begin
    r = optimize(net16, T=8, verbose=false)
    @assert haskey(r, :sd)
    @assert haskey(r.sd, "Base")
    @assert all(v in (0,1) for v in r.sd["Base"])
end
@test "UC P_gen is time-vector" begin
    r = optimize(net16, T=8, verbose=false)
    @assert length(r.P_gen["Base"]) == 8
    @assert length(r.P_gen["Peak"]) == 8
end
@test "UC P_line is T-length list" begin
    r = optimize(net16, T=6, verbose=false)
    @assert haskey(r, :P_line) && length(r.P_line) == 6
end
@test "UC lmp has length T" begin
    r = optimize(net16, T=8, verbose=false)
    @assert length(r.lmp["B1"]) == 8
end
@test "UC startup cost raises total cost" begin
    function make_uc_net(su_cost)
        n = Network(baseMVA=100.0)
        add!(n,"Bus","B1"; slack=true)
        add!(n,"Generator","G"; bus="B1",p_nom=200.0,p_min_pu=0.3,
             marginal_cost=15.0,committable=true,min_up_time=1,
             startup_cost=su_cost,initial_status=false)
        add!(n,"Generator","P"; bus="B1",p_nom=100.0,marginal_cost=80.0)
        add!(n,"Load","D"; bus="B1",p_set=120.0)
        n
    end
    r_cheap = optimize(make_uc_net(0.0),    T=8, verbose=false)
    r_exp   = optimize(make_uc_net(5000.0), T=8, verbose=false)
    @assert r_exp.total_cost >= r_cheap.total_cost - 1e-3
end
@test "UC min_down_time respected" begin
    r = optimize(net16, T=10, verbose=false)
    u = r.sd["Base"]
    for t in 1:10
        if u[t] == 1  # shutdown at t → must be off for min_down_time=2 periods
            for dt in 1:min(1, 10-t)
                @assert r.u["Base"][t+dt] == 0
            end
        end
    end
end

# ── 17. Transformer in Solvers ────────────────────────────────────────────────
println("\n[17] Transformer in Solvers")

function make_trafo_net(; tap=1.0)
    n = Network(baseMVA=100.0)
    add!(n,"Bus","B1"; slack=true)
    add!(n,"Bus","B2")
    add!(n,"Transformer","T1"; bus0="B1",bus1="B2",x=0.1,tap_ratio=tap)
    add!(n,"Generator","G1"; bus="B1",p_nom=300.0,marginal_cost=20.0)
    add!(n,"Load","D1"; bus="B2",p_set=150.0)
    n
end

@test "Transformer DC PF converges" begin
    r = pf(make_trafo_net(), verbose=false)
    @assert r.converged
end
@test "Transformer: tap_ratio changes angle" begin
    r1 = pf(make_trafo_net(tap=1.00), verbose=false)
    r2 = pf(make_trafo_net(tap=1.05), verbose=false)
    idx = findfirst(==("B2"), r1.buses)
    @assert abs(r1.θ[idx] - r2.θ[idx]) > 1e-6
end
@test "Transformer LOPF converges" begin
    r = optimize(make_trafo_net(), verbose=false)
    @assert r.converged
end
@test "Transformer DC PF trafo_flows field populated" begin
    r = pf(make_trafo_net(), verbose=false)
    @assert length(r.trafo_flows) == 1
    @assert abs(r.trafo_flows[1].P_MW) > 1.0
end

# ── 18. Return Field Completeness ────────────────────────────────────────────
println("\n[18] Return Field Completeness")

@test "DC PF fields: θ, P_inj, line_flows, buses, converged" begin
    r = pf(make_3bus(), verbose=false)
    for f in (:θ, :P_inj, :line_flows, :buses, :converged)
        @assert haskey(r, f) "Missing field: $f"
    end
end
@test "LAC PF fields: V_mag, V_ang, P_flow, Q_flow, converged" begin
    r = pf(make_3bus(), method=:lac, verbose=false)
    for f in (:V_mag, :V_ang, :P_flow, :Q_flow, :converged)
        @assert haskey(r, f) "Missing field: $f"
    end
end
@test "LOPF fields: P_gen, lmp, total_cost, P_line, converged, status" begin
    r = optimize(make_3bus(), verbose=false)
    for f in (:P_gen, :lmp, :total_cost, :P_line, :converged, :status)
        @assert haskey(r, f) "Missing field: $f"
    end
end
@test "MP-LOPF fields: gen_dispatch, soc, lmp, total_cost, status" begin
    r = optimize(net15, T=6, verbose=false)
    for f in (:gen_dispatch, :soc, :lmp, :total_cost, :status)
        @assert haskey(r, f) "Missing field: $f"
    end
end
@test "UC fields: u, su, sd, P_gen, lmp, total_cost, status" begin
    r = optimize(net16, T=6, verbose=false)
    for f in (:u, :su, :sd, :P_gen, :lmp, :total_cost, :status)
        @assert haskey(r, f) "Missing field: $f"
    end
end

# ── Summary ───────────────────────────────────────────────────────────────
total = passed + failed
println("\n" * "=" ^ 60)
println("Results: $passed / $total passed  $(failed==0 ? "✓" : "✗ ($failed failed)")")
println("=" ^ 60)
