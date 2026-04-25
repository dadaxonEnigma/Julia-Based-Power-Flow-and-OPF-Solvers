using JuMP
using HiGHS
using Random
using Printf
using Statistics

include("../solvers/lopf_multiperiod.jl")

# ----------------------------------------------------------------
#  Fixed 3-bus network (same as lopf_multiperiod.jl main script)
# ----------------------------------------------------------------
const BUSES_MP = ["Bus 1", "Bus 2", "Bus 3"]
const LINES_MP = [(1, 2, 0.01, 0.1), (1, 3, 0.01, 0.1), (2, 3, 0.01, 0.1)]
const GENS_MP  = Dict(1 => (270.0, 20.0), 2 => (100.0, 50.0))
const LOADS_MP = Dict(2 => 250.0, 3 => 175.0)
const STOR_MP  = Dict(2 => (100.0, 300.0, 0.95, 0.95))
const WIND_MP  = Dict(3 => 150.0)

function time_median(f, n_runs)
    times = Float64[]
    for _ in 1:n_runs
        push!(times, @elapsed f())
    end
    return median(times), minimum(times)
end

println("="^65)
println("BENCHMARK: Multi-Period LOPF  (Julia — JuMP + HiGHS)")
println("="^65)
println("Network: 3 buses, storage, wind; varying horizon T\n")

HORIZONS = [6, 12, 24, 48, 96]

@printf("%-8s %12s %12s %8s\n", "T (h)", "Median (ms)", "Min (ms)", "Vars")
println("-"^45)

mp_results = Dict{Int,Float64}()

for T in HORIZONS
    # JIT warmup
    solve_lopf_multiperiod(BUSES_MP, LINES_MP, GENS_MP, LOADS_MP,
                           STOR_MP, WIND_MP; T=T)

    n_runs = T <= 24 ? 20 : (T <= 48 ? 10 : 5)
    med, mn = time_median(n_runs) do
        solve_lopf_multiperiod(BUSES_MP, LINES_MP, GENS_MP, LOADS_MP,
                               STOR_MP, WIND_MP; T=T)
    end

    # Approximate variable count: T*(n_buses + n_gens + 3*n_stor)
    n_vars = T * (3 + 2 + 3)
    mp_results[T] = med * 1000
    @printf("%-8d %12.3f %12.3f %8d\n", T, med*1000, mn*1000, n_vars)
end

open(joinpath(@__DIR__, "..", "..", "results", "julia_mp_benchmark.csv"), "w") do io
    println(io, "module,T,time_ms")
    for T in HORIZONS
        println(io, "MLOPF,$T,$(mp_results[T])")
    end
end

println("\n[OK] Saved to results/julia_mp_benchmark.csv")
println("Run python/benchmark_multiperiod.py for Python comparison.")
