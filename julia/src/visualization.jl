"""
    visualization.jl

Publication-quality plots for PowerFlowJulia results.
Requires Plots.jl (add via `Pkg.add("Plots")`).

Exported functions
  plot_dispatch      — generator dispatch (single-period or multi-period stacked area)
  plot_lmp           — LMP bar (single-period) or heatmap (multi-period)
  plot_soc           — StorageUnit state-of-charge profile
  plot_uc_schedule   — Unit commitment Gantt chart
  plot_network       — Schematic bus–line topology diagram
"""

using Plots
gr()

# ── Colour palette ────────────────────────────────────────────────────────────

const _C = ["#1565C0","#2E7D32","#F57F17","#B71C1C",
            "#6A1B9A","#00838F","#4E342E","#37474F"]
_col(i) = _C[mod1(i, length(_C))]

# ── Internal helpers ──────────────────────────────────────────────────────────

# Sort generator names by marginal cost (cheapest first = bottom of stack).
function _sorted_gens(data, net)
    names = sort(collect(keys(data)))
    net === nothing && return names
    mc = [haskey(net.generators, n) ? net.generators[n].marginal_cost : 999.0
          for n in names]
    return names[sortperm(mc)]
end

# ── plot_dispatch ─────────────────────────────────────────────────────────────

"""
    plot_dispatch(result; net=nothing, title="Generator Dispatch", savepath=nothing)

Stacked area chart (multi-period) or horizontal bar (single-period).
Pass `net` to sort generators by marginal cost (cheapest at bottom).
"""
function plot_dispatch(result; net=nothing,
                       title="Generator Dispatch", savepath=nothing)
    data = hasproperty(result, :gen_dispatch) ? result.gen_dispatch :
           hasproperty(result, :P_gen)        ? result.P_gen        :
           error("plot_dispatch: result has no P_gen / gen_dispatch field")

    names = _sorted_gens(data, net)
    fv    = first(values(data))

    plt = fv isa AbstractVector ?
          _dispatch_multi(data, names, length(fv), title) :
          _dispatch_single(data, names, title)

    savepath !== nothing && savefig(plt, savepath)
    return plt
end

function _dispatch_single(data, names, title)
    vals = [max(0.0, get(data, n, 0.0)) for n in names]
    bar(names, vals;
        title=title, xlabel="Generator", ylabel="Power (MW)",
        color=reshape([_col(i) for i in 1:length(names)], 1, :),
        legend=false, bar_width=0.6, linewidth=0,
        size=(600, 380), titlefontsize=12, guidefontsize=10)
end

function _dispatch_multi(data, names, T, title)
    plt = plot(title=title, xlabel="Period", ylabel="Power (MW)",
               legend=:outertopright, size=(820, 400),
               titlefontsize=12, guidefontsize=10)
    bottom = zeros(T)
    for (i, n) in enumerate(names)
        top = bottom .+ max.(0.0, get(data, n, zeros(T))[1:T])
        plot!(plt, 1:T, top; fillrange=bottom, fillalpha=0.82,
              label=n, color=_col(i), linewidth=0.8, linecolor=_col(i))
        bottom = top
    end
    return plt
end

# ── plot_lmp ──────────────────────────────────────────────────────────────────

"""
    plot_lmp(result; title="Locational Marginal Prices", savepath=nothing)

Single-period: bar chart per bus.
Multi-period: heatmap (bus × time period).
"""
function plot_lmp(result; title="Locational Marginal Prices (€/MWh)",
                  savepath=nothing)
    lmp   = result.lmp
    buses = sort(collect(keys(lmp)))
    fv    = first(values(lmp))

    if fv isa AbstractVector
        T   = length(fv)
        mat = [lmp[b][t] for b in buses, t in 1:T]
        plt = heatmap(1:T, buses, mat;
                      title=title, xlabel="Period", ylabel="Bus",
                      color=:RdYlGn_r, colorbar_title="€/MWh",
                      size=(820, max(300, 55 * length(buses))),
                      titlefontsize=12, guidefontsize=10)
    else
        vals = [lmp[b] for b in buses]
        plt  = bar(buses, vals;
                   title=title, xlabel="Bus", ylabel="€/MWh",
                   color="#1565C0", legend=false, linewidth=0,
                   size=(600, 380), titlefontsize=12, guidefontsize=10)
    end

    savepath !== nothing && savefig(plt, savepath)
    return plt
end

# ── plot_soc ──────────────────────────────────────────────────────────────────

"""
    plot_soc(result; title="Storage State of Charge", savepath=nothing)

Line plot of StorageUnit SoC over T periods.
`result.soc[unit]` has T+1 values: index 1 = initial, 2..T+1 = end-of-period.
"""
function plot_soc(result; title="Storage State of Charge (MWh)", savepath=nothing)
    hasproperty(result, :soc) || error("plot_soc: result has no :soc field")
    soc   = result.soc
    isempty(soc) && error("plot_soc: no storage units in result")

    units = sort(collect(keys(soc)))
    T     = length(first(values(soc))) - 1   # soc has T+1 entries

    plt = plot(title=title, xlabel="Period", ylabel="Energy (MWh)",
               legend=:outertopright, size=(700, 360),
               titlefontsize=12, guidefontsize=10)
    for (i, u) in enumerate(units)
        plot!(plt, 1:T, soc[u][2:end];   # drop initial value
              label=u, color=_col(i + 3),
              linewidth=2.5, marker=:circle, markersize=5)
    end

    savepath !== nothing && savefig(plt, savepath)
    return plt
end

# ── plot_uc_schedule ──────────────────────────────────────────────────────────

"""
    plot_uc_schedule(result; title="Unit Commitment Schedule", savepath=nothing)

Gantt-style heatmap: green = ON, grey = OFF.
"""
function plot_uc_schedule(result; title="Unit Commitment Schedule", savepath=nothing)
    hasproperty(result, :u) || error("plot_uc_schedule: result has no :u field")
    u_dict = result.u
    isempty(u_dict) && error("plot_uc_schedule: no committable generators")

    gens = sort(collect(keys(u_dict)))
    T    = length(first(values(u_dict)))
    mat  = Float64[Float64(u_dict[g][t]) for g in gens, t in 1:T]

    plt = heatmap(1:T, gens, mat;
                  title=title, xlabel="Period", ylabel="Generator",
                  color=cgrad([:lightgray, "#43A047"]),
                  clims=(0, 1), colorbar=false,
                  linewidth=0.8, linecolor=:white,
                  size=(max(620, 28 * T), max(280, 58 * length(gens))),
                  titlefontsize=12, guidefontsize=10)

    savepath !== nothing && savefig(plt, savepath)
    return plt
end

# ── plot_network ──────────────────────────────────────────────────────────────

"""
    plot_network(net; result=nothing, title="Network Topology", savepath=nothing)

Schematic topology: buses on a circle, lines as edges.
Bus colour: gold = slack, blue = has generator, grey = load-only.
If `result` is provided with P_line data, edges are coloured by loading.
"""
function plot_network(net; result=nothing, title="Network Topology", savepath=nothing)
    buses = sort(collect(keys(net.buses)))
    n     = length(buses)
    n == 0 && error("plot_network: network has no buses")

    # Circular layout
    θ   = range(π/2, π/2 + 2π * (n-1)/n, length=n)
    pos = Dict(buses[i] => (cos(θ[i]), sin(θ[i])) for i in 1:n)

    plt = plot(title=title, axis=false, grid=false, legend=false,
               aspect_ratio=:equal, size=(620, 620), titlefontsize=12)

    # Edges — coloured by line loading only if result.P_line is a Dict{String,...}
    p_line_data = (result !== nothing &&
                   hasproperty(result, :P_line) &&
                   result.P_line isa AbstractDict) ? result.P_line : nothing

    for (lname, l) in net.lines
        x0, y0 = pos[l.from_bus]
        x1, y1 = pos[l.to_bus]
        loading = 0.0
        if p_line_data !== nothing && haskey(p_line_data, lname)
            flow     = p_line_data[lname]
            flow_val = flow isa AbstractVector ? maximum(abs.(flow)) : abs(flow)
            loading  = l.s_nom > 0 ? min(1.0, flow_val / l.s_nom) : 0.0
        end
        lc = loading > 0.8 ? "#C62828" : (loading > 0.5 ? "#F57F17" : "#1565C0")
        plot!(plt, [x0, x1], [y0, y1]; color=lc, linewidth=2 + 2*loading, alpha=0.75)
        mx, my = (x0+x1)/2, (y0+y1)/2
        annotate!(plt, mx, my, text(lname, :gray50, :center, 6))
    end

    # Nodes
    gens_at  = Dict{String,Bool}()
    for g in values(net.generators)
        gens_at[g.bus] = true
    end
    loads_at = Dict{String,Float64}()
    for l in values(net.loads)
        loads_at[l.bus] = get(loads_at, l.bus, 0.0) + l.p_set
    end

    for b in buses
        x, y   = pos[b]
        is_slack = net.buses[b].slack
        has_gen  = get(gens_at, b, false)
        nc = is_slack ? "#FFC107" : (has_gen ? "#1565C0" : "#78909C")
        scatter!(plt, [x], [y]; markersize=18, color=nc,
                 markerstrokecolor=:white, markerstrokewidth=2)
        annotate!(plt, x, y, text(b, :white, :center, 8, :bold))
        if haskey(loads_at, b)
            annotate!(plt, x, y - 0.18,
                      text("$(round(Int, loads_at[b])) MW", :gray40, :center, 7))
        end
    end

    savepath !== nothing && savefig(plt, savepath)
    return plt
end
