"""
Benchmark: Unit Commitment scaling — Julia (PowerFlowJulia API, JuMP + HiGHS).

Exported `unit_commitment(net; ...)` on add!-built networks at 3, 14, 30, 40 buses
(T=24). Matches python/benchmark_uc_scale.py.

Statistics: SEEDS random topologies × runs_per_seed timed solves; mean ± std + min.
"""

using Statistics
using Printf
using Random
using LinearAlgebra

BLAS.set_num_threads(1)   # fair single-thread linear algebra

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia

const SEEDS = [1, 2, 3]

function build_uc_network(n_buses; seed=42)
    rng = MersenneTwister(seed)
    net = Network(baseMVA=100.0)

    for i in 1:n_buses
        add!(net, "Bus", "B$i"; v_nom=380.0, slack=(i == 1))
    end

    added = Set{Tuple{Int,Int}}()
    for i in 1:n_buses-1
        x = 0.05 + rand(rng) * 0.45
        add!(net, "Line", "L$(i)_$(i+1)"; bus0="B$i", bus1="B$(i+1)", x=x, r=0.01)
        push!(added, (i, i+1))
    end
    for _ in 1:max(1, n_buses ÷ 3)
        u = rand(rng, 1:n_buses-1)
        v = rand(rng, u+1:n_buses)
        (u, v) in added && continue
        x = 0.05 + rand(rng) * 0.45
        add!(net, "Line", "Lx$(u)_$(v)"; bus0="B$u", bus1="B$v", x=x, r=0.01)
        push!(added, (u, v))
    end

    total_load = 0.0
    for bus in 2:n_buses
        rand(rng) > 0.3 || continue
        p = 50.0 + rand(rng) * 200.0
        add!(net, "Load", "D$bus"; bus="B$bus", p_set=p)
        total_load += p
    end
    if total_load == 0.0
        add!(net, "Load", "D2"; bus="B2", p_set=200.0)
        total_load = 200.0
    end

    add!(net, "Generator", "G_base"; bus="B1", p_nom=total_load * 1.3, marginal_cost=20.0)
    for bus in 2:4:n_buses
        p_cap = 80.0 + rand(rng) * 120.0
        mc    = 50.0 + rand(rng) * 40.0
        if rand(rng) < 0.4
            add!(net, "Generator", "G_pk$bus"; bus="B$bus", p_nom=p_cap, marginal_cost=mc,
                 committable=true, min_up_time=2, min_down_time=1,
                 startup_cost=500.0, p_min_pu=0.3)
        else
            add!(net, "Generator", "G_ct$bus"; bus="B$bus", p_nom=p_cap, marginal_cost=mc)
        end
    end
    return net
end

const LP = PowerFlowJulia.DEFAULT_LOAD_PROFILE
const WP = PowerFlowJulia.DEFAULT_WIND_PROFILE

solve_uc(net) = unit_commitment(net; T=24, load_profile=LP, wind_profile=WP,
                                compute_lmp=true, verbose=false)

function bench_stats(n_buses, runs_per_seed)
    samples = Float64[]
    for s in SEEDS
        net = build_uc_network(n_buses; seed=s)
        solve_uc(net)                           # JIT / warm-up
        for _ in 1:runs_per_seed
            push!(samples, @elapsed solve_uc(net))
        end
    end
    n = length(samples)
    return (mean = mean(samples) * 1000,
            std  = (n > 1 ? std(samples) : 0.0) * 1000,
            min  = minimum(samples) * 1000,
            n    = n)
end

println("="^74)
println("BENCHMARK: Unit Commitment Scaling  (PowerFlowJulia API, T=24)")
println("seeds $SEEDS  |  mean ± std")
println("="^74)
@printf("%-8s %12s %12s %12s %6s\n", "Buses", "Mean (ms)", "Std (ms)", "Min (ms)", "N")
println("-"^74)

SIZES = [3, 14, 30]
res = Dict{Int,NamedTuple}()
for n in SIZES
    runs = n <= 14 ? 8 : 5
    st = bench_stats(n, runs)
    res[n] = st
    @printf("%-8d %12.3f %12.3f %12.3f %6d\n", n, st.mean, st.std, st.min, st.n)
end

open(joinpath(@__DIR__, "..", "..", "results", "julia_uc_scale.csv"), "w") do io
    println(io, "module,n_buses,T,time_ms,std_ms,min_ms,n_samples")
    for n in SIZES; s=res[n]; println(io, "UC,$n,24,$(s.mean),$(s.std),$(s.min),$(s.n)"); end
end

println("\n[OK] results/julia_uc_scale.csv  |  run python/benchmark_uc_scale.py for PyPSA")
