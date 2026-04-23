using JuMP
using HiGHS
using Printf
using LinearAlgebra


function solve_lopf(buses, lines, generators, loads, line_names;
                    line_capacity = Inf,
                    baseMVA       = 100.0,
                    verbose       = true)

    n_buses = length(buses)
    n_lines = length(lines)

    B = zeros(n_buses, n_buses)
    susceptances = Float64[]
    for (from, to, r, x) in lines
        b = baseMVA / x          # [МВт/рад], x в per-unit (0.1 pu → b=1000 МВт/рад)
        push!(susceptances, b)
        B[from, from] += b
        B[to,   to  ] += b
        B[from, to  ] -= b
        B[to,   from] -= b
    end

    P_load = zeros(n_buses)
    for (bus, P) in loads
        P_load[bus] = P
    end

    gen_buses = sort(collect(keys(generators)))
    P_max_gen = Dict(bus => generators[bus][1] for bus in gen_buses)
    costs     = Dict(bus => generators[bus][2] for bus in gen_buses)
    ref_bus   = gen_buses[1]

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, θ[1:n_buses])                                         # углы [рад]
    @variable(model, P_gen[bus in gen_buses],
              lower_bound = 0.0,
              upper_bound = P_max_gen[bus])                                 # мощность [МВт]

    @constraint(model, ref_angle, θ[ref_bus] == 0.0)

    for k in 1:n_buses
        P_inj_gen = (k in gen_buses) ? P_gen[k] : 0.0
        @constraint(model,
            sum(B[k, m] * θ[m] for m in 1:n_buses) == P_inj_gen - P_load[k])
    end

    if isfinite(line_capacity)
        for (i, (from, to, r, x)) in enumerate(lines)
            b = susceptances[i]
            @constraint(model, b * (θ[from] - θ[to]) <=  line_capacity)
            @constraint(model, b * (θ[from] - θ[to]) >= -line_capacity)
        end
    end

    @objective(model, Min, sum(costs[bus] * P_gen[bus] for bus in gen_buses))

    optimize!(model)

    status = termination_status(model)
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        println("⚠️  LOPF status: $status")
        return (converged = false, status = status)
    end

    θ_val     = value.(θ)
    P_gen_val = Dict(bus => value(P_gen[bus]) for bus in gen_buses)
    P_line    = [susceptances[i] * (θ_val[from] - θ_val[to])
                 for (i, (from, to, r, x)) in enumerate(lines)]
    total_cost = objective_value(model)

    if verbose
        println("="^62)
        println("LOPF RESULTS")
        println("="^62)

        println("\n1. GENERATOR DISPATCH:")
        @printf("%-20s %10s %12s %14s\n", "Generator", "P (MW)", "P_max (MW)", "Cost (€/MWh)")
        println("-"^58)
        for bus in gen_buses
            @printf("G%-3d (%s)  %10.2f %12.2f %14.2f\n",
                    bus, buses[bus], P_gen_val[bus], P_max_gen[bus], costs[bus])
        end

        cap_str = isfinite(line_capacity) ? @sprintf("%.0f", line_capacity) : "∞"
        println("\n2. LINE FLOWS (capacity = $cap_str MW):")
        @printf("%-15s %10s %12s %12s\n", "Line", "P (MW)", "Limit (MW)", "Loading (%)")
        println("-"^52)
        for (i, name) in enumerate(line_names)
            if isfinite(line_capacity)
                loading = abs(P_line[i]) / line_capacity * 100
                @printf("%-15s %10.2f %12.1f %11.1f%%\n",
                        name, P_line[i], line_capacity, loading)
            else
                @printf("%-15s %10.2f %12s %12s\n", name, P_line[i], "∞", "—")
            end
        end

        println("\n3. VOLTAGE ANGLES:")
        for k in 1:n_buses
            @printf("  %-10s θ = %+.4f rad  (%+.3f°)\n",
                    buses[k], θ_val[k], θ_val[k] * 180/π)
        end

        @printf("\n💰 Total generation cost: %.2f €/h\n", total_cost)
        println("="^62)
    end

    return (
        θ          = θ_val,
        P_gen      = P_gen_val,
        P_line     = P_line,
        total_cost = total_cost,
        converged  = true,
        status     = status
    )
end



println("="^62)
println("LINEAR OPTIMAL POWER FLOW (LOPF)")
println("Julia  ·  JuMP + HiGHS")
println("="^62)

buses = ["Bus 1", "Bus 2", "Bus 3"]

lines      = [(1, 2, 0.01, 0.1), (1, 3, 0.01, 0.1), (2, 3, 0.01, 0.1)]
line_names = ["Line 1-2", "Line 1-3", "Line 2-3"]

generators = Dict(
    1 => (400.0, 20.0),
    2 => (300.0, 50.0),
)

loads = Dict(2 => 200.0, 3 => 300.0)

println("\n📌 SCENARIO A — No line limits (unconstrained dispatch)")
res_A = solve_lopf(buses, lines, generators, loads, line_names,
                   line_capacity = Inf)

println("\n📌 SCENARIO B — Line capacity = 200 MW (congestion case)")
res_B = solve_lopf(buses, lines, generators, loads, line_names,
                   line_capacity = 200.0)

if res_A.converged && res_B.converged
    println("\n" * "="^62)
    println("SCENARIO COMPARISON  (Unconstrained vs Congested grid)")
    println("="^62)
    @printf("\n%-28s %14s %14s\n", "Metric", "Scenario A", "Scenario B")
    println("-"^58)
    @printf("%-28s %14.1f %14.1f\n", "G1 dispatch (MW)", res_A.P_gen[1], res_B.P_gen[1])
    @printf("%-28s %14.1f %14.1f\n", "G2 dispatch (MW)", res_A.P_gen[2], res_B.P_gen[2])
    @printf("%-28s %14.2f %14.2f\n", "Max line flow (MW)",
            maximum(abs.(res_A.P_line)), maximum(abs.(res_B.P_line)))
    @printf("%-28s %14.2f %14.2f\n", "Total cost (€/h)",
            res_A.total_cost, res_B.total_cost)
    pct = (res_B.total_cost - res_A.total_cost) / res_A.total_cost * 100
    @printf("%-28s %14s %13.1f%%\n", "Cost increase", "—", pct)
    println("\n→ Grid congestion forces expensive G2 online, raising cost by $(round(pct,digits=1))%")

    println("\n" * "="^62)
    println("VALIDATION vs PyPSA")
    println("="^62)
    println("\nScenario A:")
    @printf("  G1: Julia=%.1f MW,  PyPSA expected=400.0 MW  → Δ=%.2f MW\n",
            res_A.P_gen[1], abs(res_A.P_gen[1] - 400.0))
    @printf("  G2: Julia=%.1f MW,  PyPSA expected=100.0 MW  → Δ=%.2f MW\n",
            res_A.P_gen[2], abs(res_A.P_gen[2] - 100.0))
    println("Scenario B:")
    @printf("  G1: Julia=%.1f MW,  PyPSA expected=300.0 MW  → Δ=%.2f MW\n",
            res_B.P_gen[1], abs(res_B.P_gen[1] - 300.0))
    @printf("  G2: Julia=%.1f MW,  PyPSA expected=200.0 MW  → Δ=%.2f MW\n",
            res_B.P_gen[2], abs(res_B.P_gen[2] - 200.0))
end
