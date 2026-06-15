"""
Benchmark: Unit Commitment + LMP — PowerFlowJulia API (JuMP + HiGHS).

Exported `unit_commitment(net; ...)` on an add!-built 3-bus network, varying T.
Same network as python/benchmark_uc.py.

Network: 3-bus | G_base=270MW | G_peak=150MW committable | Wind=150MW.
Fixed network ⇒ no seed sweep; mean ± std + min reported over `n_runs` solves.
"""

using Statistics
using Printf

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia

tile(base, T) = repeat(base, (T ÷ length(base)) + 1)[1:T]

function build_network()
    net = Network(baseMVA=100.0)
    add!(net, "Bus", "Bus1"; v_nom=380.0, slack=true)
    add!(net, "Bus", "Bus2"; v_nom=380.0)
    add!(net, "Bus", "Bus3"; v_nom=380.0)
    add!(net, "Line", "L12"; bus0="Bus1", bus1="Bus2", r=0.01, x=0.1, s_nom=1e6)
    add!(net, "Line", "L13"; bus0="Bus1", bus1="Bus3", r=0.01, x=0.1, s_nom=1e6)
    add!(net, "Line", "L23"; bus0="Bus2", bus1="Bus3", r=0.01, x=0.1, s_nom=1e6)
    add!(net, "Generator", "G_base"; bus="Bus1", p_nom=270.0, marginal_cost=20.0)
    add!(net, "Generator", "G_peak"; bus="Bus1", p_nom=150.0, marginal_cost=80.0,
         committable=true, min_up_time=2, min_down_time=1,
         startup_cost=500.0, shutdown_cost=0.0, p_min_pu=0.3)
    add!(net, "Generator", "Wind3"; bus="Bus3", p_nom=150.0, marginal_cost=0.0, carrier="wind")
    add!(net, "Load", "Load2"; bus="Bus2", p_set=250.0)
    add!(net, "Load", "Load3"; bus="Bus3", p_set=175.0)
    return net
end

function bench_stats(solve, n_runs)
    solve()                                     # JIT / warm-up
    samples = [@elapsed solve() for _ in 1:n_runs]
    return (mean = mean(samples) * 1000,
            std  = (n_runs > 1 ? std(samples) : 0.0) * 1000,
            min  = minimum(samples) * 1000,
            n    = n_runs)
end

println("="^60)
println("BENCHMARK: Unit Commitment + LMP  (PowerFlowJulia API)")
println("3-bus | G_base=270 | G_peak=150 committable | Wind=150")
println("="^60)
@printf("%-8s %12s %12s %12s %6s\n", "T (h)", "Mean (ms)", "Std (ms)", "Min (ms)", "N")
println("-"^60)

HORIZONS = [6, 12, 24]
net = build_network()
res = Dict{Int,NamedTuple}()
for T in HORIZONS
    lp = tile(DEFAULT_LOAD_PROFILE, T)
    wp = tile(DEFAULT_WIND_PROFILE, T)
    n_runs = T <= 24 ? 12 : 8
    st = bench_stats(() -> unit_commitment(net; T=T, load_profile=lp, wind_profile=wp,
                                           compute_lmp=true, verbose=false), n_runs)
    res[T] = st
    @printf("%-8d %12.3f %12.3f %12.3f %6d\n", T, st.mean, st.std, st.min, st.n)
end

open(joinpath(@__DIR__, "..", "..", "results", "benchmarks", "julia_uc_benchmark.csv"), "w") do io
    println(io, "module,T,time_ms,std_ms,min_ms,n_samples")
    for T in HORIZONS; s=res[T]; println(io, "UC,$T,$(s.mean),$(s.std),$(s.min),$(s.n)"); end
end

println("\n[OK] results/benchmarks/julia_uc_benchmark.csv  |  run python/benchmark_uc.py for PyPSA")
