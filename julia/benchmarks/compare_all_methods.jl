"""
compare_all_methods.jl — Julia vs PyPSA: full benchmark comparison.

Reads all results CSVs and produces:
  1. Six-panel figure (results/benchmark_all_methods.png)
  2. Console speedup table for thesis Chapter 3

Methods compared:
  DC_PF      — by network size (buses)
  LOPF       — by network size (buses)
  AC_PF      — by network size (buses)
  UC (time)  — by horizon T (hours)
  UC (scale) — by network size (buses), T=24 fixed
  MLOPF      — by horizon T (hours)

Run from repo root:
  julia julia/benchmarks/compare_all_methods.jl
"""

using Plots
using CSV
using DataFrames
using Printf
using Statistics

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT = joinpath(@__DIR__, "..", "..")
RES  = joinpath(ROOT, "results")

jl_dc   = CSV.read(joinpath(RES, "julia_benchmark.csv"),    DataFrame)
py_dc   = CSV.read(joinpath(RES, "python_benchmark.csv"),   DataFrame)
jl_ac   = CSV.read(joinpath(RES, "julia_ac_benchmark.csv"), DataFrame)
py_ac   = CSV.read(joinpath(RES, "python_ac_benchmark.csv"),DataFrame)
jl_uc   = CSV.read(joinpath(RES, "julia_uc_benchmark.csv"), DataFrame)
py_uc   = CSV.read(joinpath(RES, "python_uc_benchmark.csv"),DataFrame)
jl_ucs  = CSV.read(joinpath(RES, "julia_uc_scale.csv"),     DataFrame)
py_ucs  = CSV.read(joinpath(RES, "python_uc_scale.csv"),    DataFrame)
jl_mp   = CSV.read(joinpath(RES, "julia_mp_benchmark.csv"), DataFrame)
py_mp   = CSV.read(joinpath(RES, "python_mp_benchmark.csv"),DataFrame)

# Filter by method label
filter_mod(df, m) = df[df.module .== m, :]

dc_jl   = filter_mod(jl_dc,  "DC_PF")
dc_py   = filter_mod(py_dc,  "DC_PF")
lp_jl   = filter_mod(jl_dc,  "LOPF")
lp_py   = filter_mod(py_dc,  "LOPF")
ac_jl   = jl_ac
ac_py   = py_ac
uc_jl   = jl_uc
uc_py   = py_uc
ucs_jl  = jl_ucs
ucs_py  = py_ucs
mp_jl   = jl_mp
mp_py   = py_mp

# ── Color palette ────────────────────────────────────────────────────────────
C_PY = RGB(0.80, 0.15, 0.15)   # crimson  — PyPSA
C_JL = RGB(0.13, 0.47, 0.71)   # steel blue — Julia
GRID_A = 0.25

BASE = (
    tickfontsize  = 9,
    guidefontsize = 10,
    titlefontsize = 11,
    legendfontsize= 8,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1,1,1,0.85),
    framestyle = :box,
    grid       = true,
    gridalpha  = GRID_A,
    dpi        = 300,
)

# ── Helper: build one log-log time panel ─────────────────────────────────────
function time_panel(x_jl, y_jl, x_py, y_py;
                    title="", xlabel="Buses", ylabel="Time (ms, log)")
    p = plot(; BASE..., title=title, xlabel=xlabel, ylabel=ylabel,
               yscale=:log10, xscale=:log10, legend=:topleft)
    plot!(p, x_py, y_py; label="PyPSA", color=C_PY,
          linewidth=2, marker=:circle, markersize=5, markerstrokewidth=0)
    plot!(p, x_jl, y_jl; label="Julia", color=C_JL,
          linewidth=2, marker=:circle, markersize=5, markerstrokewidth=0)
    return p
end

# ── Six panels ───────────────────────────────────────────────────────────────
p1 = time_panel(dc_jl.n_buses, dc_jl.time_ms, dc_py.n_buses, dc_py.time_ms;
                title="(a) DC Power Flow", xlabel="Buses")

p2 = time_panel(lp_jl.n_buses, lp_jl.time_ms, lp_py.n_buses, lp_py.time_ms;
                title="(b) LOPF", xlabel="Buses")

p3 = time_panel(ac_jl.n_buses, ac_jl.time_ms, ac_py.n_buses, ac_py.time_ms;
                title="(c) AC Power Flow", xlabel="Buses")

p4 = time_panel(uc_jl.T, uc_jl.time_ms, uc_py.T, uc_py.time_ms;
                title="(d) Unit Commitment", xlabel="Horizon T (hours)")

p5 = time_panel(ucs_jl.n_buses, ucs_jl.time_ms, ucs_py.n_buses, ucs_py.time_ms;
                title="(e) UC Scaling (T=24)", xlabel="Buses")

p6 = time_panel(mp_jl.T, mp_jl.time_ms, mp_py.T, mp_py.time_ms;
                title="(f) Multi-Period LOPF", xlabel="Horizon T (hours)")

fig = plot(p1, p2, p3, p4, p5, p6;
           layout        = (2, 3),
           size          = (1500, 800),
           plot_title    = "Julia vs PyPSA — Execution Time Comparison (all methods)",
           plot_titlefontsize = 13,
           left_margin   = 8Plots.mm,
           right_margin  = 4Plots.mm,
           bottom_margin = 8Plots.mm,
           top_margin    = 6Plots.mm,
           dpi           = 300)

out = joinpath(RES, "benchmark_all_methods.png")
savefig(fig, out)
println("Saved: $out")

# ── Speedup table ─────────────────────────────────────────────────────────────
function speedup_table(label, x_col, jl_df, py_df)
    println("\n  $label")
    println("  " * "─"^58)
    @printf("  %-14s  %12s  %10s  %10s\n", x_col, "PyPSA (ms)", "Julia (ms)", "Speedup")
    println("  " * "─"^58)
    for row in eachrow(jl_df)
        x   = row[Symbol(x_col)]
        py_row = filter(r -> r[Symbol(x_col)] == x, py_df)
        isempty(py_row) && continue
        jt  = row.time_ms
        pt  = py_row.time_ms[1]
        sp  = pt / jt
        @printf("  %-14s  %12.2f  %10.3f  %8.1f×\n", string(x), pt, jt, sp)
    end
end

println("\n" * "="^62)
println("  SPEEDUP SUMMARY  Julia vs PyPSA  (median solve time)")
println("="^62)
speedup_table("DC Power Flow",      "n_buses", dc_jl,  dc_py)
speedup_table("LOPF",               "n_buses", lp_jl,  lp_py)
speedup_table("AC Power Flow",      "n_buses", ac_jl,  ac_py)
speedup_table("Unit Commitment",    "T",       uc_jl,  uc_py)
speedup_table("UC Scaling (T=24)",  "n_buses", ucs_jl, ucs_py)
speedup_table("Multi-Period LOPF",  "T",       mp_jl,  mp_py)

println("\n  Six-panel figure → results/benchmark_all_methods.png")
println("="^62)
