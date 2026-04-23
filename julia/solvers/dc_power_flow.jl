using LinearAlgebra
using SparseArrays
using Printf

struct PowerNetwork
    n_buses::Int
    buses::Vector{String}
    
    lines::Vector{Tuple{Int, Int, Float64}}
    line_names::Vector{String}
    
    generators::Dict{Int, Float64}
    
    loads::Dict{Int, Float64}
    
    slack_bus::Int
end

function build_susceptance_matrix(network::PowerNetwork, v_base=380.0)
    n = network.n_buses
    B = zeros(n, n)
    

    v_factor = v_base^2  # 380^2 = 144400
    
    for (from_bus, to_bus, x) in network.lines
        b = v_factor / x  
        
        B[from_bus, from_bus] += b
        B[to_bus, to_bus] += b
        
        B[from_bus, to_bus] -= b
        B[to_bus, from_bus] -= b
    end
    
    return B
end


function calculate_power_injections(network::PowerNetwork)
    P = zeros(network.n_buses)
    
    for (bus, power) in network.generators
        P[bus] += power
    end
    
    for (bus, power) in network.loads
        P[bus] -= power
    end
    
    return P
end

function solve_dc_power_flow(network::PowerNetwork)
    println("="^60)
    println("DC POWER FLOW - JULIA IMPLEMENTATION")
    println("="^60)
    

    println("\nStep 1: Building susceptance matrix B (with base voltage)...")
    v_base = 380.0  # кВ
    B = build_susceptance_matrix(network, v_base)
    println("Base voltage: $v_base kV")
    println("Matrix B ($(size(B))):")
    display(B)
    
    println("\n\nStep 2: Calculating power injections...")
    P = calculate_power_injections(network)
    println("Power injections P:")
    for i in 1:network.n_buses
        @printf("  %s: %.1f MW\n", network.buses[i], P[i])
    end
    
    println("\nStep 3: Removing slack bus ($(network.buses[network.slack_bus]))...")
    slack = network.slack_bus
    
    non_slack_indices = [i for i in 1:network.n_buses if i != slack]
    
    B_reduced = B[non_slack_indices, non_slack_indices]
    P_reduced = P[non_slack_indices]
    
    println("Reduced B matrix ($(size(B_reduced))):")
    display(B_reduced)
    println("\n\nReduced P vector:")
    display(P_reduced)
    
    println("\n\nStep 4: Solving B*θ = P...")
    θ_reduced = B_reduced \ P_reduced
    
    θ = zeros(network.n_buses)
    θ[non_slack_indices] = θ_reduced
    
    println("Voltage angles θ (radians):")
    for i in 1:network.n_buses
        @printf("  %s: %.6f rad\n", network.buses[i], θ[i])
    end
    
    println("\n\nStep 5: Calculating line flows...")
    line_flows = Float64[]
    
    for (idx, (from_bus, to_bus, x)) in enumerate(network.lines)
        b = v_base^2 / x 
        flow = (θ[from_bus] - θ[to_bus]) * b
        push!(line_flows, flow)
        
        @printf("  %s: %.2f MW\n", network.line_names[idx], flow)
    end
    
    println("\n" * "="^60)
    println("DC POWER FLOW COMPLETED!")
    println("="^60)
    
    return (θ=θ, line_flows=line_flows, B=B, P=P)
end


println("\n\nCreating test network (same as PyPSA example)...")

buses = ["Bus 0", "Bus 1", "Bus 2"]

lines = [
    (1, 2, 0.1),  # Line 0-1
    (1, 3, 0.1),  # Line 0-2
    (2, 3, 0.1)   # Line 1-2
]
line_names = ["Line 0-1", "Line 0-2", "Line 1-2"]

generators = Dict(1 => 500.0)  # Gen 0 на Bus 0

loads = Dict(
    2 => 300.0,  # Load 1 на Bus 1
    3 => 200.0   # Load 2 на Bus 2
)

network = PowerNetwork(
    3,           # n_buses
    buses,
    lines,
    line_names,
    generators,
    loads,
    1            # slack_bus (Bus 0 = index 1)
)

results = solve_dc_power_flow(network)

println("\n\n" * "="^60)
println("COMPARISON WITH PyPSA RESULTS")
println("="^60)

println("\nExpected (from PyPSA):")
println("  Bus 0 angle: 0.0 rad")
println("  Bus 1 angle: -0.000185 rad")
println("  Bus 2 angle: -0.000162 rad")
println("  Line 0-1 flow: 266.67 MW")
println("  Line 0-2 flow: 233.33 MW")
println("  Line 1-2 flow: -33.33 MW")

println("\nActual (Julia):")
@printf("  Bus 0 angle: %.6f rad\n", results.θ[1])
@printf("  Bus 1 angle: %.6f rad\n", results.θ[2])
@printf("  Bus 2 angle: %.6f rad\n", results.θ[3])
@printf("  Line 0-1 flow: %.2f MW\n", results.line_flows[1])
@printf("  Line 0-2 flow: %.2f MW\n", results.line_flows[2])
@printf("  Line 1-2 flow: %.2f MW\n", results.line_flows[3])