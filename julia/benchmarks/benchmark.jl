using LinearAlgebra
using SparseArrays
using JuMP
using HiGHS
using Random
using Printf
using Statistics

function generate_network(n_buses; seed=42)
    rng = MersenneTwister(seed)


    lines = Tuple{Int,Int,Float64}[]
    for i in 1:n_buses-1
        x = 0.05 + rand(rng) * 0.45         
        push!(lines, (i, i+1, x))
    end
    for _ in 1:max(1, n_buses ÷ 3)          
        u = rand(rng, 1:n_buses-1)
        v = rand(rng, u+1:n_buses)
        x = 0.05 + rand(rng) * 0.45
        push!(lines, (u, v, x))
    end

    loads      = Dict{Int,Float64}()
    total_load = 0.0
    for bus in 2:n_buses
        if rand(rng) > 0.3
            p = 50.0 + rand(rng) * 450.0
            loads[bus] = p
            total_load += p
        end
    end
    isempty(loads) && (loads[2] = 200.0; total_load = 200.0)

    generators = Dict{Int,Float64}(1 => total_load * 1.1)
    for bus in 2:4:n_buses
        generators[bus] = get(loads, bus, 0.0) * 0.5 + 50.0
    end

    return n_buses, lines, generators, loads
end

function dc_pf_solve(n_buses, lines, generators, loads; v_base=380.0)
    B = zeros(n_buses, n_buses)
    for (from, to, x) in lines
        b = v_base^2 / x
        B[from, from] += b;  B[to, to]   += b
        B[from, to]   -= b;  B[to, from] -= b
    end

    P = zeros(n_buses)
    for (bus, p) in generators;  P[bus] += p; end
    for (bus, p) in loads;       P[bus] -= p; end

    idx   = 2:n_buses
    θ_red = B[idx, idx] \ P[idx]

    θ = zeros(n_buses)
    θ[idx] = θ_red
    return θ
end

function lopf_solve(n_buses, lines, generators, loads; baseMVA=100.0)
    gen_buses = sort(collect(keys(generators)))
    P_load    = zeros(n_buses)
    for (bus, p) in loads; P_load[bus] = p; end

    B = zeros(n_buses, n_buses)
    sus = Float64[]
    for (from, to, r, x) in lines        
        b = baseMVA / x
        push!(sus, b)
        B[from, from] += b;  B[to, to]   += b
        B[from, to]   -= b;  B[to, from] -= b
    end

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, θ[1:n_buses])
    @variable(model, P_gen[bus in gen_buses],
              lower_bound = 0.0,
              upper_bound = generators[bus])

    @constraint(model, θ[gen_buses[1]] == 0.0)
    for k in 1:n_buses
        P_inj = (k in gen_buses) ? P_gen[k] : 0.0
        @constraint(model, sum(B[k,m]*θ[m] for m in 1:n_buses) == P_inj - P_load[k])
    end

    @objective(model, Min, sum(20.0 * P_gen[bus] for bus in gen_buses))
    optimize!(model)
    return objective_value(model)
end


function time_median(f, n_runs)
    times = Float64[]
    for _ in 1:n_runs
        t = @elapsed f()
        push!(times, t)
    end
    return median(times), minimum(times)
end


println("="^70)
println("BENCHMARK: Julia — DC Power Flow & LOPF")
println("="^70)

DC_SIZES   = [3, 10, 50, 100, 500, 1000, 2000]
LOPF_SIZES = [3, 10, 50, 100, 500]

println("\n📊 DC POWER FLOW BENCHMARK")
println("-"^70)
@printf("%-10s %12s %12s %10s\n", "Buses", "Median (ms)", "Min (ms)", "Lines")
println("-"^70)

dc_results = Dict{Int, Float64}()

for n in DC_SIZES
    n_buses, lines, generators, loads = generate_network(n, seed=42)

    lines_dc = [(f, t, x) for (f, t, x) in lines]

    dc_pf_solve(n_buses, lines_dc, generators, loads)

    n_runs = n <= 100 ? 200 : (n <= 500 ? 50 : 10)
    med, mn = time_median(n_runs) do
        dc_pf_solve(n_buses, lines_dc, generators, loads)
    end

    dc_results[n] = med * 1000   # в мс
    @printf("%-10d %12.4f %12.4f %10d\n", n, med*1000, mn*1000, length(lines))
end

println("\n📊 LOPF BENCHMARK  (LP: JuMP + HiGHS)")
println("-"^70)
@printf("%-10s %12s %12s %10s\n", "Buses", "Median (ms)", "Min (ms)", "Lines")
println("-"^70)

lopf_results = Dict{Int, Float64}()

for n in LOPF_SIZES
    n_buses, lines_raw, generators, loads = generate_network(n, seed=42)

    lines_lopf = [(f, t, 0.01, x) for (f, t, x) in lines_raw]

    lopf_solve(n_buses, lines_lopf, generators, loads)

    n_runs = n <= 50 ? 20 : (n <= 100 ? 10 : 5)
    med, mn = time_median(n_runs) do
        lopf_solve(n_buses, lines_lopf, generators, loads)
    end

    lopf_results[n] = med * 1000
    @printf("%-10d %12.3f %12.3f %10d\n", n, med*1000, mn*1000, length(lines_lopf))
end

open("results/julia_benchmark.csv", "w") do io
    println(io, "module,n_buses,time_ms")
    for n in DC_SIZES
        println(io, "DC_PF,$n,$(dc_results[n])")
    end
    for n in LOPF_SIZES
        println(io, "LOPF,$n,$(lopf_results[n])")
    end
end

println("\n✓ Results saved to results/julia_benchmark.csv")
println("\nRun python/benchmark.py to get Python times for comparison.")
