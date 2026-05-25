"""
Benchmark: Unit Commitment scaling — Julia (JuMP + HiGHS)
Tests UC + LMP at 3, 14, and 30 buses, T=24.
Network generated deterministically (seed=42) for reproducibility.

Compare with python/benchmark_uc_scale.py.
"""

using Statistics
using Printf
using Random

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia

# ── Network generator ────────────────────────────────────────────────────────
"""
Build a random UC-ready network of n_buses buses.
Seed is fixed for reproducibility. Roughly 30% of non-base generators
are committable peakers.
"""
function build_uc_network(n_buses; seed=42)
    rng = MersenneTwister(seed)
    net = Network(baseMVA=100.0)

    # Buses
    for i in 1:n_buses
        add!(net, "Bus", "B$i"; v_nom=380.0, slack=(i == 1))
    end

    # Spanning tree
    added = Set{Tuple{Int,Int}}()
    for i in 1:n_buses-1
        x = 0.05 + rand(rng) * 0.45
        add!(net, "Line", "L$(i)_$(i+1)"; bus0="B$i", bus1="B$(i+1)", x=x, r=0.01)
        push!(added, (i, i+1))
    end
    # Extra cross-edges (≈ n/3)
    for _ in 1:max(1, n_buses ÷ 3)
        u = rand(rng, 1:n_buses-1)
        v = rand(rng, u+1:n_buses)
        (u, v) in added && continue
        x = 0.05 + rand(rng) * 0.45
        add!(net, "Line", "Lx$(u)_$(v)"; bus0="B$u", bus1="B$v", x=x, r=0.01)
        push!(added, (u, v))
    end

    # Loads: random buses 2..n with p in [50, 250] MW
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

    # Base generator at B1 (always on, cheap)
    add!(net, "Generator", "G_base"; bus="B1",
         p_nom=total_load * 1.3, marginal_cost=20.0)

    # Additional generators: every 4th bus, ~40% committable
    for bus in 2:4:n_buses
        p_cap = 80.0 + rand(rng) * 120.0
        mc    = 50.0 + rand(rng) * 40.0
        if rand(rng) < 0.4
            add!(net, "Generator", "G_pk$bus"; bus="B$bus",
                 p_nom=p_cap, marginal_cost=mc,
                 committable=true, min_up_time=2, min_down_time=1,
                 startup_cost=500.0, p_min_pu=0.3)
        else
            add!(net, "Generator", "G_ct$bus"; bus="B$bus",
                 p_nom=p_cap, marginal_cost=mc)
        end
    end

    return net
end

function time_median(f, n_runs)
    times = Float64[]
    for _ in 1:n_runs
        push!(times, @elapsed f())
    end
    return median(times), minimum(times)
end

# ── Shared profile (24h) ─────────────────────────────────────────────────────
const LP = PowerFlowJulia.DEFAULT_LOAD_PROFILE
const WP = PowerFlowJulia.DEFAULT_WIND_PROFILE   # reuse as dummy wind

# ── JIT warmup (3-bus) ───────────────────────────────────────────────────────
let net = build_uc_network(3)
    unit_commitment(net; T=6, load_profile=LP[1:6], wind_profile=WP[1:6],
                    compute_lmp=true, verbose=false)
end

println("=" ^ 70)
println("BENCHMARK: Unit Commitment Scaling  (Julia — JuMP + HiGHS, T=24)")
println("=" ^ 70)
println("n_buses  committable  continuous  T   Median(ms)   Min(ms)")
println("-" ^ 70)

SIZES   = [3, 14, 30]
T_FIXED = 24
results = Dict{Int,Float64}()

for n in SIZES
    net = build_uc_network(n)

    n_commit = count(g.committable for g in values(net.generators) if g.carrier != "wind")
    n_cont   = count(!g.committable for g in values(net.generators) if g.carrier != "wind")

    # Warmup for this size
    unit_commitment(net; T=T_FIXED, load_profile=LP, wind_profile=WP,
                    compute_lmp=true, verbose=false)

    n_runs = n <= 10 ? 10 : 5
    med, mn = time_median(n_runs) do
        unit_commitment(net; T=T_FIXED, load_profile=LP, wind_profile=WP,
                        compute_lmp=true, verbose=false)
    end

    results[n] = med * 1000
    @printf("%-8d  %-11d  %-10d  %-3d  %10.3f  %10.3f\n",
            n, n_commit, n_cont, T_FIXED, med * 1000, mn * 1000)
end

# ── Save CSV ─────────────────────────────────────────────────────────────────
csv_path = joinpath(@__DIR__, "..", "..", "results", "julia_uc_scale.csv")
open(csv_path, "w") do io
    println(io, "module,n_buses,T,time_ms")
    for n in SIZES
        println(io, "UC,$n,$T_FIXED,$(results[n])")
    end
end

println("\n[OK] Saved to results/julia_uc_scale.csv")
println("Run python/benchmark_uc_scale.py for Python comparison.")
