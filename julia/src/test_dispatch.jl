include("dispatch.jl")

# ─── Helper: build base 3-bus network (no wind, no storage) ─────────────────
function make_3bus(name; s_nom=Inf)
    n = Network(name=name, baseMVA=100.0)
    add_bus!(n, "Bus1", v_nom=380.0, slack=true,  bus_type=3)
    add_bus!(n, "Bus2", v_nom=380.0, slack=false, bus_type=2)
    add_bus!(n, "Bus3", v_nom=380.0, slack=false, bus_type=1)
    add_line!(n, "L1-2", "Bus1", "Bus2", r=0.01, x=0.1, s_nom=s_nom)
    add_line!(n, "L1-3", "Bus1", "Bus3", r=0.01, x=0.1, s_nom=s_nom)
    add_line!(n, "L2-3", "Bus2", "Bus3", r=0.01, x=0.1, s_nom=s_nom)
    add_generator!(n, "G1", "Bus1", p_nom=400.0, marginal_cost=20.0)
    add_generator!(n, "G2", "Bus2", p_nom=300.0, marginal_cost=50.0)
    add_load!(n, "D2", "Bus2", p_set=200.0, q_set=0.0)
    add_load!(n, "D3", "Bus3", p_set=300.0, q_set=0.0)
    return n
end

println("Network summary:")
println(make_3bus("demo"))
println()

# ─── TEST 1: DC Power Flow ───────────────────────────────────────────────────
println("─"^65)
println("TEST 1: DC Power Flow")
println("─"^65)

# Net injection: G1=400, G2=300, D2=200, D3=300 → Bus1: +400, Bus2: +100, Bus3: -300
dcpf = dc_pf(make_3bus("3-bus DC PF"))
@assert dcpf.converged
println("  PASS: DC PF converged")

# ─── TEST 2: Linearized AC Power Flow ───────────────────────────────────────
println("\n" * "─"^65)
println("TEST 2: Linearized AC Power Flow")
println("─"^65)

lacpf = linear_ac_pf(make_3bus("3-bus LACPF"))
@assert lacpf.converged
@assert all(lacpf.V_mag .> 0.99)      "voltage magnitude should be near 1 p.u."
@assert abs(lacpf.V_ang[1]) < 1e-10   "slack bus angle must be 0"
println("  PASS: LACPF converged, slack θ = $(lacpf.V_ang[1]) rad")

# ─── TEST 3a: LOPF — no congestion ──────────────────────────────────────────
println("\n" * "─"^65)
println("TEST 3a: LOPF — unconstrained dispatch")
println("─"^65)
# Load = 500 MW; G1 (cheap) should carry 400, G2 = 100
res_A = lopf(make_3bus("3-bus LOPF-A"))
@assert res_A.converged
@assert abs(res_A.P_gen["G1"] - 400.0) < 0.5  "G1 should dispatch 400 MW"
@assert abs(res_A.P_gen["G2"] - 100.0) < 0.5  "G2 should dispatch 100 MW"
@assert abs(res_A.total_cost - 13_000.0) < 1.0
println("  G1 = $(round(res_A.P_gen["G1"],digits=1)) MW  (expected 400.0) ✓")
println("  G2 = $(round(res_A.P_gen["G2"],digits=1)) MW  (expected 100.0) ✓")
println("  Cost = $(round(res_A.total_cost,digits=0)) €/h  (expected 13000) ✓")

# ─── TEST 3b: LOPF — line congestion ────────────────────────────────────────
println("\n" * "─"^65)
println("TEST 3b: LOPF — line capacity 200 MW (congestion)")
println("─"^65)
res_B = lopf(make_3bus("3-bus LOPF-B", s_nom=200.0))
@assert res_B.converged
@assert abs(res_B.P_gen["G1"] - 300.0) < 0.5  "Congested: G1 should dispatch 300 MW"
@assert abs(res_B.P_gen["G2"] - 200.0) < 0.5  "Congested: G2 should dispatch 200 MW"
@assert abs(res_B.total_cost - 16_000.0) < 1.0
pct = (res_B.total_cost - res_A.total_cost) / res_A.total_cost * 100
println("  G1 = $(round(res_B.P_gen["G1"],digits=1)) MW  (expected 300.0) ✓")
println("  G2 = $(round(res_B.P_gen["G2"],digits=1)) MW  (expected 200.0) ✓")
println("  Cost increase: $(round(pct,digits=1))%  (expected +23.1%) ✓")

# ─── TEST 4: Multi-period LOPF with storage and wind ────────────────────────
println("\n" * "─"^65)
println("TEST 4: Multi-period LOPF (24h, storage + wind)")
println("─"^65)

net_mp = Network(name="3-bus multi-period", baseMVA=100.0)
add_bus!(net_mp, "Bus1", v_nom=380.0, slack=true,  bus_type=3)
add_bus!(net_mp, "Bus2", v_nom=380.0, slack=false, bus_type=2)
add_bus!(net_mp, "Bus3", v_nom=380.0, slack=false, bus_type=1)
add_line!(net_mp, "L1-2", "Bus1", "Bus2", r=0.01, x=0.1)
add_line!(net_mp, "L1-3", "Bus1", "Bus3", r=0.01, x=0.1)
add_line!(net_mp, "L2-3", "Bus2", "Bus3", r=0.01, x=0.1)
add_generator!(net_mp, "G1",    "Bus1", p_nom=270.0, marginal_cost=20.0, carrier="gas")
add_generator!(net_mp, "G2",    "Bus2", p_nom=100.0, marginal_cost=50.0, carrier="coal")
add_generator!(net_mp, "Wind3", "Bus3", p_nom=150.0, marginal_cost=0.0,  carrier="wind")
add_load!(net_mp, "D2", "Bus2", p_set=250.0)
add_load!(net_mp, "D3", "Bus3", p_set=175.0)
add_storage_unit!(net_mp, "Bat2", "Bus2",
    p_nom=100.0, e_nom=300.0,
    efficiency_charge=0.95, efficiency_discharge=0.95,
    cyclic_state_of_charge=true)

res_mp = lopf_multiperiod(net_mp, T=24)

@assert string(res_mp.status) == "OPTIMAL"
@assert haskey(res_mp.gen_dispatch, "G1")
@assert haskey(res_mp.soc, "Bat2")
@assert res_mp.total_cost > 0
# Cyclic SoC: E[end] == E[0]
soc_vec = res_mp.soc["Bat2"]
@assert abs(soc_vec[end] - soc_vec[1]) < 1.0  "Cyclic SoC violated"
println("  Status: $(res_mp.status) ✓")
println("  Total cost: $(round(res_mp.total_cost, digits=0)) €")
println("  Bat2 SoC: $(round(soc_vec[1],digits=1)) → $(round(soc_vec[end],digits=1)) MWh (cyclic) ✓")

# ─── Summary ─────────────────────────────────────────────────────────────────
println()
println("="^65)
println("All 4 dispatch tests passed.")
println("="^65)
