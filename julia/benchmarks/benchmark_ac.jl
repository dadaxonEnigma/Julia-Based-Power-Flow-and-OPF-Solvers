using PowerModels
using Ipopt
using JuMP
using Random
using Printf
using Statistics

include("../solvers/ac_power_flow.jl")

# ----------------------------------------------------------------
#  Silent solver (reused across all benchmark calls)
# ----------------------------------------------------------------
const SILENT_IPOPT = optimizer_with_attributes(
    Ipopt.Optimizer, "print_level" => 0, "sb" => "yes"
)

# ----------------------------------------------------------------
#  Network generator — AC format (same random seed as DC/LOPF)
# ----------------------------------------------------------------
function generate_network_ac(n_buses; seed=42)
    rng = MersenneTwister(seed)

    buses = ["Bus$i" for i in 1:n_buses]

    lines      = Tuple{Int,Int,Float64,Float64}[]
    line_names = String[]

    # Spanning tree — guarantees connectivity
    for i in 1:n_buses-1
        x = 0.05 + rand(rng) * 0.45
        push!(lines,      (i, i+1, 0.01, x))   # (from, to, r, x)
        push!(line_names, "L$i-$(i+1)")
    end
    # Additional mesh edges (~n/3)
    for _ in 1:max(1, n_buses ÷ 3)
        u = rand(rng, 1:n_buses-1)
        v = rand(rng, u+1:n_buses)
        x = 0.05 + rand(rng) * 0.45
        push!(lines,      (u, v, 0.01, x))
        push!(line_names, "L$u-$v")
    end

    # Loads at ~70% of non-slack buses
    loads      = Dict{Int, Tuple{Float64,Float64}}()
    total_load = 0.0
    for bus in 2:n_buses
        if rand(rng) > 0.3
            p = 50.0 + rand(rng) * 450.0
            loads[bus] = (p, 0.0)           # (P_MW, Q_MVAr=0)
            total_load += p
        end
    end
    isempty(loads) && (loads[2] = (200.0, 0.0); total_load = 200.0)

    # Generators: slack at bus 1, PV at every 4th bus
    generators = Dict{Int, Tuple{Float64,Float64}}()
    generators[1] = (total_load * 1.1, 1.0)    # (P_MW, V_set p.u.)
    for bus in 2:4:n_buses
        p_gen = get(loads, bus, (0.0, 0.0))[1] * 0.5 + 50.0
        generators[bus] = (p_gen, 1.0)
    end

    return buses, lines, line_names, generators, loads
end

# ----------------------------------------------------------------
#  Benchmark kernel: data construction + AC PF solve
# ----------------------------------------------------------------
function ac_pf_bench(buses, lines, generators, loads)
    pm_data = network_to_powermodels(buses, lines, generators, loads)
    result  = solve_ac_pf(pm_data, SILENT_IPOPT)
    return result["termination_status"]
end

function time_median(f, n_runs)
    times = Float64[]
    for _ in 1:n_runs
        push!(times, @elapsed f())
    end
    return median(times), minimum(times)
end

# ================================================================
#  BENCHMARK
# ================================================================
println("="^70)
println("BENCHMARK: Julia — AC Power Flow (PowerModels.jl + Ipopt)")
println("="^70)
println("Note: times include pm_data construction + Ipopt solve (post-JIT).")

AC_SIZES = [3, 10, 50, 100]

println("\n[AC POWER FLOW BENCHMARK]")
println("-"^70)
@printf("%-10s %12s %12s %10s\n", "Buses", "Median (ms)", "Min (ms)", "Lines")
println("-"^70)

ac_results = Dict{Int,Float64}()

for n in AC_SIZES
    buses, lines, line_names, generators, loads = generate_network_ac(n)

    # JIT warmup — compile all PowerModels / Ipopt code paths
    ac_pf_bench(buses, lines, generators, loads)

    n_runs = n <= 10 ? 20 : (n <= 50 ? 10 : 5)
    med, mn = time_median(n_runs) do
        ac_pf_bench(buses, lines, generators, loads)
    end

    ac_results[n] = med * 1000
    @printf("%-10d %12.3f %12.3f %10d\n", n, med*1000, mn*1000, length(lines))
end

# ----------------------------------------------------------------
#  Save CSV (appends AC_PF rows; run alongside DC/LOPF benchmarks)
# ----------------------------------------------------------------
open(joinpath(@__DIR__, "..", "..", "results", "julia_ac_benchmark.csv"), "w") do io
    println(io, "module,n_buses,time_ms")
    for n in AC_SIZES
        println(io, "AC_PF,$n,$(ac_results[n])")
    end
end

println("\n[OK] Saved to results/julia_ac_benchmark.csv")
println("Run python/benchmark_ac.py to get Python times for comparison.")
