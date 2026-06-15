"""
Benchmark: Julia — DC Power Flow & LOPF, via the PowerFlowJulia LIBRARY API.

Benchmarks the *exported* `pf(net)` / `optimize(net)` entry points — exactly the
API the thesis describes — NOT a hand-written kernel, so the Julia-vs-PyPSA
comparison is apples-to-apples (both build a network through their public API
and call their public solver, both LP via HiGHS).

Statistics: each size is measured over SEEDS random topologies (identical to
python/benchmark.py), with `runs_per_seed` timed solves each after a JIT
warm-up. Reported: mean ± std and min over the pooled samples. The CSV keeps
`time_ms` = mean for backward compatibility, plus `std_ms`, `min_ms`.
"""

using Random
using Printf
using Statistics
using LinearAlgebra

BLAS.set_num_threads(1)   # fair single-thread linear algebra (DC PF solve)

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia

const SEEDS = [1, 2, 3]

# ── Network generator (identical to python/benchmark.py) ─────────────────────
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

# ── Build a PowerFlowJulia Network from the generated topology ───────────────
function build_network(n_buses, lines, generators, loads)
    net = Network(baseMVA=100.0)
    for b in 1:n_buses
        add!(net, "Bus", "Bus$b"; v_nom=380.0, slack=(b == 1))
    end
    for (i, (f, t, x)) in enumerate(lines)
        add!(net, "Line", "L$i"; bus0="Bus$f", bus1="Bus$t", x=x, r=0.01, s_nom=1e6)
    end
    for (bus, p_nom) in generators
        add!(net, "Generator", "G$bus"; bus="Bus$bus", p_nom=p_nom, marginal_cost=20.0)
    end
    for (bus, p) in loads
        add!(net, "Load", "D$bus"; bus="Bus$bus", p_set=p)
    end
    return net
end

# Collect timing samples over SEEDS × runs_per_seed; return mean/std/min (ms).
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
println("BENCHMARK: Julia (PowerFlowJulia API) — DC Power Flow & LOPF")
println("Seeds: $SEEDS  |  mean ± std over pooled samples")
println("="^74)

DC_SIZES   = [3, 10, 50, 100, 300]
LOPF_SIZES = [3, 10, 50, 100, 300]

println("\n[DC POWER FLOW]  pf(net; method=:dc)")
println("-"^74)
@printf("%-8s %12s %12s %12s %6s\n", "Buses", "Mean (ms)", "Std (ms)", "Min (ms)", "N")
println("-"^74)

dc = Dict{Int,NamedTuple}()
for n in DC_SIZES
    runs = n <= 100 ? 10 : (n <= 500 ? 8 : (n <= 1000 ? 5 : 3))
    st = bench_stats(s -> build_network(generate_network(n, seed=s)...),
                     net -> pf(net; method=:dc, verbose=false), runs)
    dc[n] = st
    @printf("%-8d %12.4f %12.4f %12.4f %6d\n", n, st.mean, st.std, st.min, st.n)
end

println("\n[LOPF]  optimize(net; method=:lopf)   (JuMP + HiGHS)")
println("-"^74)
@printf("%-8s %12s %12s %12s %6s\n", "Buses", "Mean (ms)", "Std (ms)", "Min (ms)", "N")
println("-"^74)

lp = Dict{Int,NamedTuple}()
for n in LOPF_SIZES
    runs = n <= 50 ? 10 : (n <= 100 ? 8 : 4)
    st = bench_stats(s -> build_network(generate_network(n, seed=s)...),
                     net -> optimize(net; method=:lopf, verbose=false), runs)
    lp[n] = st
    @printf("%-8d %12.4f %12.4f %12.4f %6d\n", n, st.mean, st.std, st.min, st.n)
end

open(joinpath(@__DIR__, "..", "..", "results", "benchmarks", "julia_benchmark.csv"), "w") do io
    println(io, "module,n_buses,time_ms,std_ms,min_ms,n_samples")
    for n in DC_SIZES;   s=dc[n]; println(io, "DC_PF,$n,$(s.mean),$(s.std),$(s.min),$(s.n)"); end
    for n in LOPF_SIZES; s=lp[n]; println(io, "LOPF,$n,$(s.mean),$(s.std),$(s.min),$(s.n)"); end
end

println("\n[OK] results/benchmarks/julia_benchmark.csv  |  run python/benchmark.py for PyPSA")
