"""
Benchmark: Julia — full AC Power Flow, via the PowerFlowJulia LIBRARY API.

Benchmarks the exported `pf(net; method=:ac)` entry point (PowerModels.jl + Ipopt),
NOT a standalone kernel. Same seed-based topology generator as the DC/LOPF
benchmark, matching python/benchmark_ac.py.

Statistics: SEEDS random topologies × runs_per_seed timed solves; mean ± std + min.
"""

using Random
using Printf
using Statistics
using LinearAlgebra

BLAS.set_num_threads(1)   # fair single-thread linear algebra

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia

const SEEDS = [1, 2, 3]

function generate_network_ac(n_buses; seed=42)
    rng = MersenneTwister(seed)

    lines = Tuple{Int,Int,Float64,Float64}[]
    for i in 1:n_buses-1
        x = 0.05 + rand(rng) * 0.45
        push!(lines, (i, i+1, 0.01, x))
    end
    for _ in 1:max(1, n_buses ÷ 3)
        u = rand(rng, 1:n_buses-1)
        v = rand(rng, u+1:n_buses)
        x = 0.05 + rand(rng) * 0.45
        push!(lines, (u, v, 0.01, x))
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

function build_network_ac(n_buses, lines, generators, loads)
    net = Network(baseMVA=100.0)
    for b in 1:n_buses
        add!(net, "Bus", "Bus$b"; v_nom=380.0, slack=(b == 1))
    end
    for (i, (f, t, r, x)) in enumerate(lines)
        add!(net, "Line", "L$i"; bus0="Bus$f", bus1="Bus$t", r=r, x=x, s_nom=1e6)
    end
    for (bus, p_nom) in generators
        add!(net, "Generator", "G$bus"; bus="Bus$bus", p_nom=p_nom, marginal_cost=20.0)
    end
    for (bus, p) in loads
        add!(net, "Load", "D$bus"; bus="Bus$bus", p_set=p)
    end
    return net
end

function bench_stats(build, solve, runs_per_seed)
    samples = Float64[]
    for s in SEEDS
        net = build(s)
        solve(net)                              # JIT / warm-up
        for _ in 1:runs_per_seed
            push!(samples, @elapsed solve(net))
        end
    end
    n = length(samples)
    return (mean = mean(samples) * 1000,
            std  = (n > 1 ? std(samples) : 0.0) * 1000,
            min  = minimum(samples) * 1000,
            n    = n)
end

println("="^74)
println("BENCHMARK: Julia (PowerFlowJulia API) — AC Power Flow")
println("pf(net; method=:ac) → PowerModels.jl + Ipopt  |  seeds $SEEDS")
println("="^74)

AC_SIZES = [3, 10, 50, 100]

println("\n[AC POWER FLOW]")
println("-"^74)
@printf("%-8s %12s %12s %12s %6s\n", "Buses", "Mean (ms)", "Std (ms)", "Min (ms)", "N")
println("-"^74)

ac = Dict{Int,NamedTuple}()
for n in AC_SIZES
    runs = n <= 10 ? 10 : (n <= 50 ? 8 : 5)
    st = bench_stats(s -> build_network_ac(generate_network_ac(n, seed=s)...),
                     net -> pf(net; method=:ac, verbose=false), runs)
    ac[n] = st
    @printf("%-8d %12.3f %12.3f %12.3f %6d\n", n, st.mean, st.std, st.min, st.n)
end

open(joinpath(@__DIR__, "..", "..", "results", "benchmarks", "julia_ac_benchmark.csv"), "w") do io
    println(io, "module,n_buses,time_ms,std_ms,min_ms,n_samples")
    for n in AC_SIZES; s=ac[n]; println(io, "AC_PF,$n,$(s.mean),$(s.std),$(s.min),$(s.n)"); end
end

println("\n[OK] results/benchmarks/julia_ac_benchmark.csv  |  run python/benchmark_ac.py for PyPSA")
