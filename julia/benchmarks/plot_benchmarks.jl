using Plots
using StatsPlots
using CSV
using DataFrames
using Printf

# ─── Load data ──────────────────────────────────────────────────────────────
python_data = CSV.read("../results/python_benchmark.csv", DataFrame)
julia_data  = CSV.read("../results/julia_benchmark.csv",  DataFrame)

dc_py   = python_data[python_data.module .== "DC_PF", :]
dc_jl   = julia_data[julia_data.module   .== "DC_PF", :]
lopf_py = python_data[python_data.module .== "LOPF",  :]
lopf_jl = julia_data[julia_data.module   .== "LOPF",  :]

# ─── Compute speedups ────────────────────────────────────────────────────────
speedup_dc   = Float64[]
speedup_lopf = Float64[]
sizes_dc     = Int[]
sizes_lopf   = Int[]

for (py_mod, jl_mod, sp_vec, sz_vec) in [
        (dc_py,   dc_jl,   speedup_dc,   sizes_dc),
        (lopf_py, lopf_jl, speedup_lopf, sizes_lopf)]
    for i in 1:nrow(py_mod)
        n      = py_mod.n_buses[i]
        jl_idx = findfirst(jl_mod.n_buses .== n)
        if jl_idx !== nothing
            push!(sp_vec, py_mod.time_ms[i] / jl_mod.time_ms[jl_idx])
            push!(sz_vec, n)
        end
    end
end

# ─── Color palette (academic) ────────────────────────────────────────────────
C_PY_DC   = RGB(0.80, 0.15, 0.15)   # crimson  — Python DC
C_JL_DC   = RGB(0.13, 0.47, 0.71)   # steel blue — Julia DC
C_PY_LOPF = RGB(0.93, 0.53, 0.18)   # orange   — Python LOPF
C_JL_LOPF = RGB(0.17, 0.63, 0.17)   # forest green — Julia LOPF
GRID_ALPHA = 0.25

# Shared style keyword sets
BASE = (
    tickfontsize   = 10,
    guidefontsize  = 12,
    titlefontsize  = 14,
    legendfontsize = 10,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.85),
    framestyle = :box,
    grid       = true,
    gridalpha  = GRID_ALPHA,
    dpi  = 300,
)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 1 — Execution Time (two-panel side by side)                    ║
# ╚══════════════════════════════════════════════════════════════════════════╝

pa = plot(;
    BASE...,
    title   = "DC Power Flow",
    xlabel  = "Network size (buses)",
    ylabel  = "Execution time (ms, log scale)",
    yscale  = :log10,
    legend  = :topleft,
    size    = (600, 480),
    xscale  = :log10,
)

plot!(pa, dc_py.n_buses, dc_py.time_ms;
      label="PyPSA (Python)", color=C_PY_DC,
      linewidth=2.5, marker=:circle, markersize=6,
      markerstrokewidth=0)
plot!(pa, dc_jl.n_buses, dc_jl.time_ms;
      label="Julia", color=C_JL_DC,
      linewidth=2.5, marker=:circle, markersize=6,
      markerstrokewidth=0)

# Shade the gap between curves
fill_x  = vcat(dc_py.n_buses, reverse(dc_jl.n_buses))
fill_y  = vcat(dc_py.time_ms, reverse(dc_jl.time_ms))
plot!(pa, fill_x, fill_y;
      seriestype=:shape, fillalpha=0.10,
      color=C_PY_DC, linewidth=0, label="")

pb = plot(;
    BASE...,
    title   = "LOPF (Linear Optimal Power Flow)",
    xlabel  = "Network size (buses)",
    ylabel  = "",
    yscale  = :log10,
    legend  = :topleft,
    size    = (600, 480),
    xscale  = :log10,
)

plot!(pb, lopf_py.n_buses, lopf_py.time_ms;
      label="PyPSA (Python)", color=C_PY_LOPF,
      linewidth=2.5, marker=:square, markersize=6,
      markerstrokewidth=0)
plot!(pb, lopf_jl.n_buses, lopf_jl.time_ms;
      label="Julia", color=C_JL_LOPF,
      linewidth=2.5, marker=:square, markersize=6,
      markerstrokewidth=0)

fill_x2 = vcat(lopf_py.n_buses, reverse(lopf_jl.n_buses))
fill_y2 = vcat(lopf_py.time_ms, reverse(lopf_jl.time_ms))
plot!(pb, fill_x2, fill_y2;
      seriestype=:shape, fillalpha=0.10,
      color=C_PY_LOPF, linewidth=0, label="")

p1 = plot(pa, pb;
          layout    = (1, 2),
          size      = (1200, 520),
          plot_title = "Execution Time: Julia vs Python/PyPSA",
          plot_titlefontsize = 15,
          left_margin  = 5Plots.mm,
          right_margin = 5Plots.mm,
          bottom_margin = 6Plots.mm,
          dpi = 300)

savefig(p1, "../results/benchmark_time.png")
println("Saved: results/benchmark_time.png")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 2 — Speedup factor                                             ║
# ╚══════════════════════════════════════════════════════════════════════════╝

p2 = plot(;
    BASE...,
    title   = "Speedup Factor: Julia vs Python/PyPSA",
    xlabel  = "Network size (buses)",
    ylabel  = "Speedup  (×,  log scale)",
    yscale  = :log10,
    xscale  = :log10,
    legend  = :bottomleft,
    size    = (900, 580),
)

plot!(p2, sizes_dc, speedup_dc;
      label="DC Power Flow", color=C_JL_DC,
      linewidth=3, marker=:circle, markersize=8,
      markerstrokewidth=1, markerstrokecolor=:white)
plot!(p2, sizes_lopf, speedup_lopf;
      label="LOPF", color=C_JL_LOPF,
      linewidth=3, marker=:square, markersize=8,
      markerstrokewidth=1, markerstrokecolor=:white)
hline!(p2, [1.0];
       color=:gray40, linestyle=:dash, linewidth=1.5,
       label="1× (no speedup)")

# Annotate speedup values next to markers
for (n, sp) in zip(sizes_dc, speedup_dc)
    lbl = sp > 999 ? @sprintf("%dk×", round(Int, sp/1000)) : @sprintf("%d×", round(Int, sp))
    annotate!(p2, n * 1.25, sp * 1.8,
              Plots.text(lbl, 8, C_JL_DC, :left))
end
for (n, sp) in zip(sizes_lopf, speedup_lopf)
    lbl = @sprintf("%d×", round(Int, sp))
    annotate!(p2, n * 1.25, sp / 2.5,
              Plots.text(lbl, 8, C_JL_LOPF, :left))
end

savefig(p2, "../results/benchmark_speedup.png")
println("Saved: results/benchmark_speedup.png")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 3 — Grouped bar chart with speedup annotations                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝

selected = [10, 100, 500]
sp_dc_sel   = [speedup_dc[findfirst(sizes_dc .== s)]
               for s in selected if s in sizes_dc]
sp_lopf_sel = [speedup_lopf[findfirst(sizes_lopf .== s)]
               for s in selected if s in sizes_lopf]

bar_data   = [sp_dc_sel sp_lopf_sel]   # n_groups × n_series
group_lbls = string.(selected) .* " buses"

p3 = groupedbar(bar_data;
    bar_position = :dodge,
    bar_width    = 0.6,
    label        = ["DC Power Flow" "LOPF"],
    color        = [C_JL_DC C_JL_LOPF],
    alpha        = 0.85,
    title        = "Speedup by Network Size  (Julia vs Python/PyPSA)",
    xlabel       = "Network size",
    ylabel       = "Speedup  (×)",
    xticks       = (1:length(selected), group_lbls),
    legend       = :topright,
    BASE...,
    size         = (860, 560),
    yscale       = :log10,
    ylims        = (10, 300_000),
    dpi          = 300,
)

# Add value labels above each bar
n_grp = length(selected)
bar_w = 0.6
offsets = [-bar_w/4, bar_w/4]   # two series: left, right

for (s, series_vals) in enumerate([sp_dc_sel, sp_lopf_sel])
    for (g, val) in enumerate(series_vals)
        x_pos = g + offsets[s]
        lbl   = val > 999 ?
                @sprintf("%.0fk×", val / 1000) :
                @sprintf("%d×", round(Int, val))
        annotate!(p3, x_pos, val * 2.2,
                  Plots.text(lbl, 9, :center,
                             s == 1 ? C_JL_DC : C_JL_LOPF))
    end
end

savefig(p3, "../results/benchmark_speedup_bars.png")
println("Saved: results/benchmark_speedup_bars.png")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 4 — Combined 2×2 "hero" figure for thesis                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# Panel A: DC time
pA = plot(;
    BASE..., title="(a) DC Power Flow — Time",
    xlabel="Buses", ylabel="Time (ms, log)",
    yscale=:log10, xscale=:log10, legend=:topleft,
)
plot!(pA, dc_py.n_buses, dc_py.time_ms; label="PyPSA",
      color=C_PY_DC,  linewidth=2, marker=:circle, markersize=5)
plot!(pA, dc_jl.n_buses, dc_jl.time_ms; label="Julia",
      color=C_JL_DC,  linewidth=2, marker=:circle, markersize=5)

# Panel B: LOPF time
pB = plot(;
    BASE..., title="(b) LOPF — Time",
    xlabel="Buses", ylabel="Time (ms, log)",
    yscale=:log10, xscale=:log10, legend=:topleft,
)
plot!(pB, lopf_py.n_buses, lopf_py.time_ms; label="PyPSA",
      color=C_PY_LOPF, linewidth=2, marker=:square, markersize=5)
plot!(pB, lopf_jl.n_buses, lopf_jl.time_ms; label="Julia",
      color=C_JL_LOPF, linewidth=2, marker=:square, markersize=5)

# Panel C: DC speedup
pC = plot(;
    BASE..., title="(c) DC Power Flow — Speedup",
    xlabel="Buses", ylabel="Speedup (×, log)", yscale=:log10, xscale=:log10,
    legend=false,
)
plot!(pC, sizes_dc, speedup_dc;
      color=C_JL_DC, linewidth=2.5, marker=:circle, markersize=6)
hline!(pC, [1.0]; color=:gray50, linestyle=:dash, linewidth=1)

# Panel D: LOPF speedup
pD = plot(;
    BASE..., title="(d) LOPF — Speedup",
    xlabel="Buses", ylabel="Speedup (×)", yscale=:log10, xscale=:log10,
    legend=false,
)
plot!(pD, sizes_lopf, speedup_lopf;
      color=C_JL_LOPF, linewidth=2.5, marker=:square, markersize=6)
hline!(pD, [1.0]; color=:gray50, linestyle=:dash, linewidth=1)

p4 = plot(pA, pB, pC, pD;
          layout        = (2, 2),
          size          = (1200, 900),
          plot_title    = "Julia vs Python/PyPSA  —  Power System Benchmark Results",
          plot_titlefontsize = 14,
          left_margin   = 6Plots.mm,
          right_margin  = 4Plots.mm,
          bottom_margin = 6Plots.mm,
          top_margin    = 4Plots.mm,
          dpi           = 300)

savefig(p4, "../results/benchmark_combined.png")
println("Saved: results/benchmark_combined.png")

# ─── Summary table ───────────────────────────────────────────────────────────
println("\n" * "="^62)
println("  BENCHMARK SUMMARY  (for thesis Chapter 3)")
println("="^62)

println("\n  DC POWER FLOW")
@printf("  %-12s  %12s  %10s  %12s\n", "Size (buses)", "Python (ms)", "Julia (ms)", "Speedup")
println("  " * "-"^52)
for i in 1:length(sizes_dc)
    @printf("  %-12d  %12.2f  %10.4f  %10.1f×\n",
            sizes_dc[i], dc_py.time_ms[i], dc_jl.time_ms[i], speedup_dc[i])
end

println("\n  LOPF")
@printf("  %-12s  %12s  %10s  %12s\n", "Size (buses)", "Python (ms)", "Julia (ms)", "Speedup")
println("  " * "-"^52)
for i in 1:length(sizes_lopf)
    @printf("  %-12d  %12.2f  %10.4f  %10.1f×\n",
            sizes_lopf[i], lopf_py.time_ms[i], lopf_jl.time_ms[i], speedup_lopf[i])
end

println("\n  4 figures saved to results/")
