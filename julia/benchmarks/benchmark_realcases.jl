#=
benchmark_realcases.jl — performance benchmark on STANDARD published networks
(IEEE and PEGASE MATPOWER cases), not random topologies.

Loads each MATPOWER `.m` case into a PowerFlowJulia `Network` and times the
**exported API** (`pf` for DC power flow, `optimize` for LOPF) — the same entry
points a user calls — after a JIT warm-up. Results are written to
`results/julia_realcases_benchmark.csv`.

Cases (place the .m files in data/):
  case118.m        IEEE 118-bus
  case300.m        IEEE 300-bus
  case1354pegase.m PEGASE 1354-bus (part of the European transmission grid)
  case2869pegase.m PEGASE 2869-bus (European grid)

Run:
  julia --project=. benchmarks/benchmark_realcases.jl
=#
include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia
using PowerModels
using BenchmarkTools
using Statistics, Printf
using LinearAlgebra

BLAS.set_num_threads(1)   # fair single-thread linear algebra (DC PF solve)
PowerModels.silence()

DATA = joinpath(@__DIR__, "..", "..", "data")
CASES = ["case118.m", "case300.m", "case1354pegase.m", "case2869pegase.m"]

# ── MATPOWER .m → PowerFlowJulia Network (via exported add! API) ─────────────
# Branches with tap ≠ 1 (transformers) are folded into an equivalent line
# reactance x·tap; this preserves problem size and DC behaviour for the timing
# study. Per-unit quantities from PowerModels are rescaled to MW by baseMVA.
function load_matpower(path::String)
    data = PowerModels.parse_file(path)
    base = data["baseMVA"]
    net  = Network(name = basename(path), baseMVA = base)

    for (_, b) in data["bus"]
        i  = b["index"]
        kv = get(b, "base_kv", 0.0); kv <= 0 && (kv = 380.0)
        add!(net, "Bus", "b$(i)"; v_nom = kv, slack = (b["bus_type"] == 3))
    end
    for (k, br) in data["branch"]
        br["br_status"] == 0 && continue
        tap  = get(br, "tap", 1.0); tap == 0 && (tap = 1.0)
        rate = get(br, "rate_a", 0.0)
        snom = rate > 0 ? rate * base : 1.0e6
        add!(net, "Line", "l$(k)";
             bus0 = "b$(br["f_bus"])", bus1 = "b$(br["t_bus"])",
             r = get(br, "br_r", 0.0), x = br["br_x"] * tap, s_nom = snom)
    end
    for (k, g) in data["gen"]
        g["gen_status"] == 0 && continue
        pmax = g["pmax"] * base
        pmax <= 0 && continue
        cost = (haskey(g, "cost") && length(g["cost"]) >= 2) ? g["cost"][end-1] : 10.0
        cost <= 0 && (cost = 10.0)
        add!(net, "Generator", "g$(k)"; bus = "b$(g["gen_bus"])",
             p_nom = pmax, marginal_cost = float(cost))
    end
    for (k, ld) in data["load"]
        add!(net, "Load", "d$(k)"; bus = "b$(ld["load_bus"])", p_set = ld["pd"] * base)
    end
    return net
end

# Time a zero-argument call with BenchmarkTools: auto-tuned repeated samples,
# GC-aware, reporting the minimum (noise only adds time) alongside median/std.
# evals=1 so each sample is a single full call (model build + solve).
function timed(f)
    f()                                   # warm-up — exclude JIT compilation
    b = @benchmark ($f)() samples = 50 seconds = 25 evals = 1
    return median(b).time / 1e6, std(b.times) / 1e6, minimum(b).time / 1e6
end

rows = Tuple[]
println("="^70)
println("  Real-case benchmark — standard IEEE/PEGASE networks (exported API)")
println("="^70)
for c in CASES
    path = joinpath(DATA, c)
    if !isfile(path)
        @printf("  %-18s  SKIP (not found)\n", c); continue
    end
    net = load_matpower(path)
    nb  = length(net.buses)
    @printf("\n[%s]  %d buses, %d lines, %d gens, %d loads\n",
            c, nb, length(net.lines), length(net.generators), length(net.loads))

    for (label, f) in (("DC_PF", () -> pf(net; verbose = false)),
                       ("LOPF",  () -> optimize(net; verbose = false)))
        try
            m, s, mn = timed(f)
            push!(rows, (c, nb, label, m, s, mn))
            @printf("    %-6s  median=%.2f ms  std=%.2f  min=%.2f\n", label, m, s, mn)
        catch e
            @printf("    %-6s  FAILED: %s\n", label, sprint(showerror, e))
        end
    end
end

out = joinpath(@__DIR__, "..", "..", "results", "julia_realcases_benchmark.csv")
open(out, "w") do io
    println(io, "case,n_buses,module,time_ms,std_ms,min_ms")
    for r in rows
        @printf(io, "%s,%d,%s,%.6f,%.6f,%.6f\n", r...)
    end
end
println("\n[OK] $(length(rows)) rows → $out")
