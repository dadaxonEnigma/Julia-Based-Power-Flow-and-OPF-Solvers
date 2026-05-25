"""
Benchmark: Unit Commitment + LMP — PowerFlowJulia (JuMP + HiGHS)
Same network as python/benchmark_uc.py for direct comparison.

Network: 3-bus, G_base (continuous 270 MW), G_peak (committable 150 MW),
         Wind 150 MW (zero cost), Load2=250 MW·profile, Load3=175 MW·profile
"""

using Statistics
using Printf

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia

const LOAD_PROFILE_24 = [
    0.60, 0.57, 0.55, 0.54, 0.55, 0.60,
    0.70, 0.80, 0.88, 0.90, 0.92, 0.91,
    0.90, 0.89, 0.88, 0.87, 0.89, 0.95,
    1.00, 0.98, 0.93, 0.85, 0.75, 0.65,
]

const WIND_PROFILE_24 = [
    0.80, 0.82, 0.85, 0.83, 0.78, 0.70,
    0.60, 0.55, 0.50, 0.45, 0.42, 0.40,
    0.38, 0.37, 0.40, 0.43, 0.50, 0.58,
    0.65, 0.70, 0.74, 0.76, 0.78, 0.80,
]

function tile_profile(base, T)
    reps = (T ÷ length(base)) + 1
    return repeat(base, reps)[1:T]
end

function build_network()
    net = Network(baseMVA=100.0)
    add!(net, "Bus", "Bus1"; v_nom=380.0, slack=true)
    add!(net, "Bus", "Bus2"; v_nom=380.0)
    add!(net, "Bus", "Bus3"; v_nom=380.0)
    add!(net, "Line", "L12"; bus0="Bus1", bus1="Bus2", r=0.01, x=0.1, s_nom=1e6)
    add!(net, "Line", "L13"; bus0="Bus1", bus1="Bus3", r=0.01, x=0.1, s_nom=1e6)
    add!(net, "Line", "L23"; bus0="Bus2", bus1="Bus3", r=0.01, x=0.1, s_nom=1e6)
    # Continuous base-load generator
    add!(net, "Generator", "G_base";
         bus="Bus1", p_nom=270.0, marginal_cost=20.0)
    # Committable peaker: binary on/off, min_up=2h, startup=500€, p_min=30%
    add!(net, "Generator", "G_peak";
         bus="Bus1", p_nom=150.0, marginal_cost=80.0,
         committable=true, min_up_time=2, min_down_time=1,
         startup_cost=500.0, shutdown_cost=0.0, p_min_pu=0.3)
    # Zero-cost wind (variable output via wind_profile)
    add!(net, "Generator", "Wind3";
         bus="Bus3", p_nom=150.0, marginal_cost=0.0, carrier="wind")
    add!(net, "Load", "Load2"; bus="Bus2", p_set=250.0)
    add!(net, "Load", "Load3"; bus="Bus3", p_set=175.0)
    return net
end

function time_median(f, n_runs)
    times = Float64[]
    for _ in 1:n_runs
        push!(times, @elapsed f())
    end
    return median(times), minimum(times)
end

# ── JIT warmup ────────────────────────────────────────────────────────────────
let net = build_network(),
    lp  = tile_profile(LOAD_PROFILE_24, 6),
    wp  = tile_profile(WIND_PROFILE_24, 6)
    unit_commitment(net; T=6, load_profile=lp, wind_profile=wp,
                    compute_lmp=true, verbose=false)
end

println("=" ^ 65)
println("BENCHMARK: Unit Commitment + LMP  (Julia — JuMP + HiGHS)")
println("=" ^ 65)
println("Network: 3-bus | G_base=270MW | G_peak=150MW committable | Wind=150MW")
println()
@printf("%-8s %12s %12s\n", "T (h)", "Median (ms)", "Min (ms)")
println("-" ^ 38)

HORIZONS = [6, 12, 24, 48]
uc_results = Dict{Int,Float64}()
net = build_network()   # immutable — safe to reuse across horizons

for T in HORIZONS
    lp = tile_profile(LOAD_PROFILE_24, T)
    wp = tile_profile(WIND_PROFILE_24, T)
    n_runs = T <= 24 ? 10 : 5

    med, mn = time_median(n_runs) do
        unit_commitment(net; T=T, load_profile=lp, wind_profile=wp,
                        compute_lmp=true, verbose=false)
    end

    uc_results[T] = med * 1000
    @printf("%-8d %12.3f %12.3f\n", T, med * 1000, mn * 1000)
end

# ── LMP snapshot for T=24 ─────────────────────────────────────────────────────
println("\n--- LMP snapshot (T=24) ---")
lp24 = tile_profile(LOAD_PROFILE_24, 24)
wp24 = tile_profile(WIND_PROFILE_24, 24)
r24  = unit_commitment(net; T=24, load_profile=lp24, wind_profile=wp24,
                       compute_lmp=true, verbose=false)
for b in ["Bus1", "Bus2", "Bus3"]
    v = r24.lmp[b]
    @printf("  %-5s  avg=%7.3f  min=%7.3f  max=%7.3f  €/MWh\n",
            b, sum(v)/24, minimum(v), maximum(v))
end

# ── Save CSV ──────────────────────────────────────────────────────────────────
csv_path = joinpath(@__DIR__, "..", "..", "results", "julia_uc_benchmark.csv")
open(csv_path, "w") do io
    println(io, "module,T,time_ms")
    for T in HORIZONS
        println(io, "UC,$T,$(uc_results[T])")
    end
end

println("\n[OK] Saved to results/julia_uc_benchmark.csv")
println("Run python/benchmark_uc.py for Python comparison.")
