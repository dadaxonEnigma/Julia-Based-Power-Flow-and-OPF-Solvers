# Unit Commitment demo — 3-bus, 24h
include("../src/PowerFlowJulia.jl")
using .PowerFlowJulia
using Printf

net = Network(name="UC demo", baseMVA=100.0)

add!(net, "Carrier", "gas";  co2_emissions=0.20)
add!(net, "Carrier", "coal"; co2_emissions=0.34)
add!(net, "Carrier", "wind"; co2_emissions=0.00)

add!(net, "Bus", "B1"; v_nom=380.0, slack=true)
add!(net, "Bus", "B2"; v_nom=380.0)
add!(net, "Bus", "B3"; v_nom=380.0)
add!(net, "Line", "L12"; bus0="B1", bus1="B2", x=0.1)
add!(net, "Line", "L13"; bus0="B1", bus1="B3", x=0.1)
add!(net, "Line", "L23"; bus0="B2", bus1="B3", x=0.1)

# Baseload: cheap, committable, slow start
add!(net, "Generator", "BaseLoad";
     bus="B1", p_nom=300.0, p_min_pu=0.4, marginal_cost=15.0,
     carrier="coal", committable=true,
     min_up_time=4, min_down_time=4,
     startup_cost=5000.0, shutdown_cost=1000.0,
     initial_status=true)

# Peaker: expensive, fast, no min time
add!(net, "Generator", "Peaker";
     bus="B2", p_nom=150.0, p_min_pu=0.0, marginal_cost=70.0,
     carrier="gas", committable=true,
     min_up_time=1, min_down_time=1,
     startup_cost=300.0, initial_status=false)

# Wind: non-committable, zero cost
add!(net, "Generator", "Wind";
     bus="B3", p_nom=120.0, marginal_cost=0.0,
     carrier="wind", committable=false)

add!(net, "Load", "D2"; bus="B2", p_set=150.0)
add!(net, "Load", "D3"; bus="B3", p_set=100.0)

# CO₂ cap
add!(net, "GlobalConstraint", "co2_budget";
     type="co2_limit", sense="<=", constant=1500.0,
     carrier_weightings=Dict("gas"=>0.20, "coal"=>0.34))

println(net)

# optimize! auto-selects Unit Commitment because committable=true generators exist
result = optimize(net, method=:uc, T=24)

println("\nUnit Commitment result:")
@printf("Total cost: %.0f €\n", result.total_cost)
println("Status: ", result.status)
