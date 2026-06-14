include("dispatch.jl")

# ─── Carriers ────────────────────────────────────────────────────────────────
println("─"^60)
println("TEST 1: Carrier")
println("─"^60)

net = Network(name="component test", baseMVA=100.0)
add_carrier!(net, "gas",  co2_emissions=0.20, color="#ff8c00")
add_carrier!(net, "coal", co2_emissions=0.34, color="#333333")
add_carrier!(net, "wind", co2_emissions=0.0,  color="#4caf50")

@assert length(net.carriers) == 3
@assert net.carriers["gas"].co2_emissions ≈ 0.20
println("  PASS: 3 carriers added, co2 values correct")

# ─── Store ───────────────────────────────────────────────────────────────────
println("\n" * "─"^60)
println("TEST 2: Store in single-period LOPF")
println("─"^60)

net2 = Network(name="Store test", baseMVA=100.0)
add_bus!(net2, "B1", slack=true)
add_bus!(net2, "B2")
add_line!(net2, "L12", "B1", "B2", x=0.1)
add_generator!(net2, "G1", "B1", p_nom=200.0, marginal_cost=30.0)
add_load!(net2, "D2", "B2", p_set=150.0)
# Store at B2: 100 MWh, 50 MW power — can help meet load
add_store!(net2, "S1", "B2", e_nom=100.0, p_nom=50.0, marginal_cost=5.0)

res2 = lopf(net2, verbose=false)
@assert res2.converged
@assert haskey(res2.P_store, "S1")
println("  Status: $(res2.status)")
println("  G1 = $(round(res2.P_gen["G1"], digits=1)) MW")
println("  Store S1 = $(round(res2.P_store["S1"], digits=1)) MW  (+ = inject)")
println("  Cost = $(round(res2.total_cost, digits=1)) €/h")
println("  PASS")

# ─── Link ─────────────────────────────────────────────────────────────────────
println("\n" * "─"^60)
println("TEST 3: Link (HVDC) in LOPF")
println("─"^60)

# Two isolated areas connected only by a Link (HVDC cable, η=0.97)
# Area A: cheap generator (10 €/MWh), Area B: expensive generator (80 €/MWh) + load
# Optimal: export from A to B via HVDC
net3 = Network(name="HVDC test", baseMVA=100.0)
add_bus!(net3, "A", slack=true)
add_bus!(net3, "B")
add_generator!(net3, "GA", "A", p_nom=500.0, marginal_cost=10.0)
add_generator!(net3, "GB", "B", p_nom=200.0, marginal_cost=80.0)
add_load!(net3, "DB", "B", p_set=300.0)
add_link!(net3, "HVDC", "A", "B", p_nom=400.0, efficiency=0.97, marginal_cost=2.0)

res3 = lopf(net3, verbose=true)
@assert res3.converged
@assert haskey(res3.P_link, "HVDC")
p_hvdc = res3.P_link["HVDC"]
# HVDC should carry power (cheap GA exports to expensive area B)
@assert p_hvdc > 0  "HVDC should flow A→B"
println("  HVDC flow: $(round(p_hvdc, digits=1)) MW (A→B)")
println("  Delivered to B: $(round(p_hvdc * 0.97, digits=1)) MW")
println("  PASS")

# ─── GlobalConstraint (CO₂ cap) ──────────────────────────────────────────────
println("\n" * "─"^60)
println("TEST 4: GlobalConstraint — CO₂ cap in multi-period LOPF")
println("─"^60)

net4 = Network(name="CO2 cap test", baseMVA=100.0)
add_carrier!(net4, "gas",  co2_emissions=0.20)
add_carrier!(net4, "coal", co2_emissions=0.34)

add_bus!(net4, "B1", slack=true)
add_bus!(net4, "B2")
add_line!(net4, "L12", "B1", "B2", x=0.1)
add_generator!(net4, "Gas1",  "B1", p_nom=300.0, marginal_cost=40.0, carrier="gas")
add_generator!(net4, "Coal1", "B1", p_nom=300.0, marginal_cost=20.0, carrier="coal")
add_load!(net4, "D2", "B2", p_set=200.0)

# Without CO₂ cap: coal (cheap) should dominate
res_nocap = lopf_multiperiod(net4, T=6, verbose=false)
coal_total = sum(res_nocap.gen_dispatch["Coal1"])
gas_total  = sum(res_nocap.gen_dispatch["Gas1"])
println("  Without cap: Coal=$(round(coal_total,digits=0)) MWh, Gas=$(round(gas_total,digits=0)) MWh")

# With CO₂ cap: shift dispatch from coal to gas
# Total load ≈ 682 MWh (200 MW × 6h × avg_profile ≈ 0.568)
# All-coal: 682 × 0.34 ≈ 232 t | All-gas: 682 × 0.20 ≈ 136 t
# Cap = 180 t → forces max coal ≈ 311 MWh (rest from gas)
add_global_constraint!(net4, "co2_cap",
    type               = "co2_limit",
    sense              = "<=",
    constant           = 180.0,
    carrier_weightings = Dict("coal" => 0.34, "gas" => 0.20))

res_cap = lopf_multiperiod(net4, T=6, verbose=false)
coal_capped = sum(res_cap.gen_dispatch["Coal1"])
gas_capped  = sum(res_cap.gen_dispatch["Gas1"])

@assert coal_capped < coal_total  "CO₂ cap should reduce coal dispatch"
@assert res_cap.total_cost > res_nocap.total_cost  "CO₂ cap should increase cost"
println("  With cap:    Coal=$(round(coal_capped,digits=0)) MWh, Gas=$(round(gas_capped,digits=0)) MWh")
println("  Cost increase: $(round(res_cap.total_cost - res_nocap.total_cost, digits=0)) €")
println("  PASS: CO₂ cap correctly constrains coal dispatch")

# ─── Summary ─────────────────────────────────────────────────────────────────
println()
println("="^60)
println("All new-component tests passed.")
println("Carrier, Store, Link, GlobalConstraint — fully integrated.")
println("="^60)
