include("dispatch.jl")
mean(v) = sum(v) / length(v)

println("="^65)
println("LMP VALIDATION TESTS")
println("="^65)

# ── TEST 1: Uncongested single-bus ────────────────────────────────────────────
# One generator at 20 €/MWh, one load. LMP should equal marginal cost = 20.
println("\n[1] Single-bus uncongested — LMP == marginal cost")
net1 = Network(name="1-bus", baseMVA=100.0)
add_bus!(net1, "B1", slack=true)
add_generator!(net1, "G1", "B1", p_nom=400.0, marginal_cost=20.0)
add_load!(net1, "D1", "B1", p_set=200.0)

r1 = lopf(net1, verbose=false)
@assert r1.converged
lmp_b1 = r1.lmp["B1"]
@printf("  LMP(B1) = %.4f €/MWh  (expected 20.0)\n", lmp_b1)
@assert abs(lmp_b1 - 20.0) < 0.01  "LMP should equal marginal cost"
println("  PASS ✓")

# ── TEST 2: Uncongested 3-bus — all LMPs equal cheapest dispatched generator ─
# G1 cheap (20), G2 expensive (50). Load = 500 MW. G1 dispatched fully, G2 at margin.
# LMP at all buses should equal marginal cost of marginal generator.
println("\n[2] 3-bus uncongested — all LMPs equal marginal generator cost")
net2 = Network(name="3-bus uncongested", baseMVA=100.0)
add_bus!(net2, "B1", slack=true)
add_bus!(net2, "B2"); add_bus!(net2, "B3")
add_line!(net2, "L12", "B1", "B2", x=0.1)
add_line!(net2, "L13", "B1", "B3", x=0.1)
add_line!(net2, "L23", "B2", "B3", x=0.1)
add_generator!(net2, "G1", "B1", p_nom=400.0, marginal_cost=20.0)
add_generator!(net2, "G2", "B2", p_nom=300.0, marginal_cost=50.0)
add_load!(net2, "D2", "B2", p_set=200.0)
add_load!(net2, "D3", "B3", p_set=300.0)

r2 = lopf(net2, verbose=false)
@assert r2.converged
# G1 fully dispatched (400 MW), G2 = 100 MW → G2 is marginal → LMP = 50
@printf("  G1 = %.1f MW, G2 = %.1f MW\n", r2.P_gen["G1"], r2.P_gen["G2"])
@printf("  LMP: B1=%.4f  B2=%.4f  B3=%.4f  (expected ≈ 50.0)\n",
        r2.lmp["B1"], r2.lmp["B2"], r2.lmp["B3"])
# All LMPs should be equal (no congestion) and equal to G2's marginal cost
lmp_vals = collect(values(r2.lmp))
@assert maximum(lmp_vals) - minimum(lmp_vals) < 0.01  "No congestion → equal LMPs"
@assert abs(mean(lmp_vals) - 50.0) < 0.1              "LMP = marginal gen cost"
println("  PASS ✓")

# ── TEST 3: Congested 2-bus — LMP diverge due to line limit ──────────────────
# Classic congestion: cheap gen at Bus1, load at Bus2, line limited to 80 MW.
# Bus2 needs expensive local gen → LMP(B2) > LMP(B1).
println("\n[3] 2-bus congested — LMP divergence (congestion rent)")
net3 = Network(name="2-bus congested", baseMVA=100.0)
add_bus!(net3, "B1", slack=true)
add_bus!(net3, "B2")
add_line!(net3, "L12", "B1", "B2", x=0.1, s_nom=80.0)   # tight limit
add_generator!(net3, "G1", "B1", p_nom=500.0, marginal_cost=10.0)  # cheap
add_generator!(net3, "G2", "B2", p_nom=200.0, marginal_cost=60.0)  # expensive
add_load!(net3, "D2", "B2", p_set=150.0)

r3 = lopf(net3, verbose=true)
@assert r3.converged
lmp1 = r3.lmp["B1"]; lmp2 = r3.lmp["B2"]
congestion_rent = lmp2 - lmp1
@printf("  LMP(B1)=%.4f, LMP(B2)=%.4f, rent=%.4f €/MWh\n",
        lmp1, lmp2, congestion_rent)
@assert lmp2 > lmp1 + 0.1  "Congestion: LMP at load bus > LMP at gen bus"
@assert congestion_rent > 0
println("  PASS ✓ — congestion rent = $(round(congestion_rent, digits=4)) €/MWh")

# ── TEST 4: Multi-period LMP ──────────────────────────────────────────────────
println("\n[4] Multi-period LOPF — time-varying LMP")
net4 = Network(name="multiperiod LMP", baseMVA=100.0)
add_bus!(net4, "B1", slack=true); add_bus!(net4, "B2")
add_line!(net4, "L12", "B1", "B2", x=0.1, s_nom=100.0)
add_generator!(net4, "G1", "B1", p_nom=300.0, marginal_cost=20.0)
add_generator!(net4, "G2", "B2", p_nom=200.0, marginal_cost=80.0)
add_load!(net4, "D2", "B2", p_set=200.0)

r4 = lopf_multiperiod(net4, T=6, verbose=false)
@assert string(r4.status) == "OPTIMAL"
@assert haskey(r4.lmp, "B1") && haskey(r4.lmp, "B2")
@assert length(r4.lmp["B1"]) == 6

# LMP at B2 should be ≥ LMP at B1 (or equal when uncongested)
for t in 1:6
    @assert r4.lmp["B2"][t] >= r4.lmp["B1"][t] - 1e-4
end
@printf("  B1 LMP: avg=%.2f  min=%.2f  max=%.2f\n",
        mean(r4.lmp["B1"]), minimum(r4.lmp["B1"]), maximum(r4.lmp["B1"]))
@printf("  B2 LMP: avg=%.2f  min=%.2f  max=%.2f\n",
        mean(r4.lmp["B2"]), minimum(r4.lmp["B2"]), maximum(r4.lmp["B2"]))
println("  PASS ✓")

# ── TEST 5: LMP with GlobalConstraint (CO₂ cap) ──────────────────────────────
# CO₂ cap forces more expensive generator → LMP rises above unconstrained level
println("\n[5] LMP with CO₂ GlobalConstraint")
net5 = Network(name="co2 LMP", baseMVA=100.0)
add_carrier!(net5, "coal", co2_emissions=0.34)
add_carrier!(net5, "gas",  co2_emissions=0.20)
add_bus!(net5, "B1", slack=true)
add_generator!(net5, "Coal", "B1", p_nom=300.0, marginal_cost=20.0, carrier="coal")
add_generator!(net5, "Gas",  "B1", p_nom=300.0, marginal_cost=50.0, carrier="gas")
add_load!(net5, "D1", "B1", p_set=200.0)

# Without cap: only coal dispatched → LMP = 20
r5a = lopf(net5, verbose=false)
# With CO₂ cap = 55 t:
#   All-coal = 200×0.34 = 68 t  (too high, violates cap)
#   All-gas  = 200×0.20 = 40 t  (feasible)
#   Constraint forces coal ≤ 107 MW, rest from gas → gas is marginal → LMP = 50
add_global_constraint!(net5, "co2",
    constant=55.0, carrier_weightings=Dict("coal"=>0.34, "gas"=>0.20))
r5b = lopf(net5, verbose=false)

@printf("  Without cap: LMP=%.2f €/MWh  (coal marginal cost=20)\n",  r5a.lmp["B1"])
@printf("  With CO₂ cap: LMP=%.2f €/MWh  (gas forced, marginal cost=50)\n", r5b.lmp["B1"])
@assert r5a.lmp["B1"] < r5b.lmp["B1"]  "CO₂ cap should raise LMP"
println("  PASS ✓ — CO₂ cap raises LMP by $(round(r5b.lmp["B1"]-r5a.lmp["B1"],digits=2)) €/MWh")

println("\n" * "="^65)
println("All LMP tests passed ✓")
println("="^65)
