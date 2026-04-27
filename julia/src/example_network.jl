include("network.jl")

# ── IEEE 3-bus example (same topology used across all solvers) ──────────────

net = Network(name="IEEE 3-bus example", baseMVA=100.0)

add_bus!(net, "Bus1", v_nom=380.0, slack=true,  bus_type=3)
add_bus!(net, "Bus2", v_nom=380.0, slack=false, bus_type=2)
add_bus!(net, "Bus3", v_nom=380.0, slack=false, bus_type=1)

add_line!(net, "L12", "Bus1", "Bus2", r=0.01, x=0.1, s_nom=200.0)
add_line!(net, "L13", "Bus1", "Bus3", r=0.01, x=0.1, s_nom=200.0)
add_line!(net, "L23", "Bus2", "Bus3", r=0.01, x=0.1, s_nom=200.0)

add_generator!(net, "G1", "Bus1", p_nom=400.0, marginal_cost=20.0, carrier="gas")
add_generator!(net, "G2", "Bus2", p_nom=300.0, marginal_cost=50.0, carrier="coal")

add_load!(net, "D2", "Bus2", p_set=200.0)
add_load!(net, "D3", "Bus3", p_set=300.0)

add_storage_unit!(net, "Bat1", "Bus3",
    p_nom=50.0, e_nom=200.0,
    efficiency_charge=0.95, efficiency_discharge=0.95,
    cyclic_state_of_charge=true)

println(net)
println()
println("Slack buses  : ", slack_buses(net))
println("Gen capacity : ", total_p_nom(net), " MW")
println("Peak load    : ", total_load(net), " MW")
println("Bus index    : ", bus_index(net))
println("Branches@Bus1: ", length(branches_at(net, "Bus1")), " branches")

# ── Transformer example (IEEE 14-bus style) ──────────────────────────────────

net14 = Network(name="IEEE 14-bus (partial)", baseMVA=100.0)
add_bus!(net14, "B1",  v_nom=132.0, slack=true, bus_type=3)
add_bus!(net14, "B4",  v_nom=132.0, bus_type=1)
add_bus!(net14, "B5",  v_nom=132.0, bus_type=1)
add_bus!(net14, "B7",  v_nom=69.0,  bus_type=1)
add_bus!(net14, "B9",  v_nom=69.0,  bus_type=1)

add_transformer!(net14, "T4-7",  "B4", "B7",  x=0.20912, s_nom=100.0, tap_ratio=0.978)
add_transformer!(net14, "T4-9",  "B4", "B9",  x=0.55618, s_nom=100.0, tap_ratio=0.969)
add_transformer!(net14, "T5-6",  "B5", "B7",  x=0.25202, s_nom=100.0, tap_ratio=0.932)

println()
println(net14)
