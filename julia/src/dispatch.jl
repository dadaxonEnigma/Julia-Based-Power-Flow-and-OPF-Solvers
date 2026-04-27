"""
    dispatch.jl

Network-aware solver dispatch via Julia multiple dispatch.
Each function accepts a `Network` struct, extracts component data,
and runs the corresponding power flow / optimisation algorithm.

Solvers provided:
  dc_pf(net)              — DC Power Flow
  linear_ac_pf(net)       — Linearized AC Power Flow
  lopf(net)               — Single-period Linear OPF
  lopf_multiperiod(net)   — 24-hour LOPF with storage and wind
"""

include("network.jl")

using LinearAlgebra
using JuMP
using HiGHS
using Printf

# ────────────────────────────────────────────────────────────────────────────
#  Default time-series profiles (24-hour, normalised)
# ────────────────────────────────────────────────────────────────────────────

const DEFAULT_LOAD_PROFILE = [
    0.60, 0.57, 0.55, 0.54, 0.55, 0.60,   # 00-05  night
    0.70, 0.80, 0.88, 0.90, 0.92, 0.91,   # 06-11  morning ramp
    0.90, 0.89, 0.88, 0.87, 0.89, 0.95,   # 12-17  midday
    1.00, 0.98, 0.93, 0.85, 0.75, 0.65    # 18-23  evening peak
]

const DEFAULT_WIND_PROFILE = [
    0.80, 0.82, 0.85, 0.83, 0.78, 0.70,   # 00-05  strong night wind
    0.60, 0.55, 0.50, 0.45, 0.42, 0.40,   # 06-11  morning drop
    0.38, 0.37, 0.40, 0.43, 0.50, 0.58,   # 12-17  afternoon recovery
    0.65, 0.70, 0.74, 0.76, 0.78, 0.80    # 18-23  evening wind
]

# ────────────────────────────────────────────────────────────────────────────
#  Shared internal helpers
# ────────────────────────────────────────────────────────────────────────────

# Build susceptance matrix B and per-branch susceptance vector.
# Transformer model (simplified DC, consistent with PyPSA / MATPOWER):
#   b_eff = (baseMVA / x) / tap_ratio
#   Both diagonal entries receive b_eff (symmetric approximation).
#   Phase shift is noted in the Transformer struct but not applied in DC PF.
function _build_B(net::Network, bidx::Dict{String,Int})
    n   = length(bidx)
    B   = zeros(n, n)
    sus = Float64[]
    branch_names = String[]

    for l in values(net.lines)
        f = bidx[l.from_bus]; t = bidx[l.to_bus]
        b = net.baseMVA / l.x
        push!(sus, b); push!(branch_names, l.name)
        B[f,f]+=b; B[t,t]+=b; B[f,t]-=b; B[t,f]-=b
    end
    for tr in values(net.transformers)
        f = bidx[tr.from_bus]; t = bidx[tr.to_bus]
        b_eff = net.baseMVA / (tr.x * tr.tap_ratio)   # PyPSA-style tap correction
        push!(sus, b_eff); push!(branch_names, tr.name)
        B[f,f]+=b_eff; B[t,t]+=b_eff; B[f,t]-=b_eff; B[t,f]-=b_eff
    end
    return B, sus, branch_names
end

# Build full nodal admittance matrix Y = G + jB for LACPF.
# Transformer model (simplified): y_eff = y / tap_ratio (both sides symmetric).
# Phase shift is not yet applied (conservative linearisation).
function _build_Y(net::Network, bidx::Dict{String,Int})
    n = length(bidx)
    G = zeros(n, n); B = zeros(n, n)
    for l in values(net.lines)
        f = bidx[l.from_bus]; t = bidx[l.to_bus]
        d    = l.r^2 + l.x^2
        g_ij = l.r / d; b_ij = l.x / d
        G[f,f]+=g_ij; G[t,t]+=g_ij; G[f,t]-=g_ij; G[t,f]-=g_ij
        B[f,f]+=b_ij; B[t,t]+=b_ij; B[f,t]-=b_ij; B[t,f]-=b_ij
    end
    for tr in values(net.transformers)
        f = bidx[tr.from_bus]; t = bidx[tr.to_bus]
        a    = tr.tap_ratio
        d    = tr.r^2 + tr.x^2
        g_ij = tr.r / d / a    # g_eff = g / a
        b_ij = tr.x / d / a    # b_eff = b / a
        G[f,f]+=g_ij; G[t,t]+=g_ij; G[f,t]-=g_ij; G[t,f]-=g_ij
        B[f,f]+=b_ij; B[t,t]+=b_ij; B[f,t]-=b_ij; B[t,f]-=b_ij
    end
    return G .* net.baseMVA, B .* net.baseMVA  # scale to [MW/rad], [MW/p.u.-V]
end

# Phase-shift correction vector for DC PF and LOPF.
# A transformer with phase shift φ [deg] contributes a fixed equivalent injection:
#   P_shift[from] += b_eff · φ_rad
#   P_shift[to]   -= b_eff · φ_rad
# The nodal equation becomes:  B·θ = P_inj + P_shift
# Branch flow:  P_ft = b_eff · (θ_f − θ_t − φ_rad)
function _phase_shift_injection(net::Network, bidx::Dict{String,Int})
    n       = length(bidx)
    P_shift = zeros(n)
    for tr in values(net.transformers)
        tr.phase_shift == 0.0 && continue
        f     = bidx[tr.from_bus]; t = bidx[tr.to_bus]
        φ     = tr.phase_shift * π / 180
        b_eff = net.baseMVA / (tr.x * tr.tap_ratio)
        P_shift[f] += b_eff * φ
        P_shift[t] -= b_eff * φ
    end
    return P_shift
end

function _slack_index(net::Network, bidx::Dict{String,Int})
    sl = slack_buses(net)
    isempty(sl) && error("Network \"$(net.name)\" has no slack bus")
    return bidx[sl[1]]
end

# ────────────────────────────────────────────────────────────────────────────
#  DC Power Flow
# ────────────────────────────────────────────────────────────────────────────

"""
    dc_pf(net::Network; verbose=true) → NamedTuple

DC Power Flow:  B·θ = P_inj,  solved via sparse LU factorisation.

Generator active injections are taken as `p_nom × p_max_pu` (rated output).
For a given dispatch, set `p_max_pu` accordingly before calling.

Returns:
  (θ, buses, line_flows, B, P_inj, converged)
  where `line_flows` is a Vector of (name, from, to, P_MW).

PyPSA equivalent: `network.lpf()`
"""
function dc_pf(net::Network; verbose=true)
    isempty(net.buses) && error("Network has no buses")

    bnames = bus_names(net)
    bidx   = bus_index(net)
    n      = length(bnames)
    slack  = _slack_index(net, bidx)

    B, sus, _ = _build_B(net, bidx)
    P_shift   = _phase_shift_injection(net, bidx)

    # Net active-power injections [MW]
    P = zeros(n)
    for g in values(net.generators); P[bidx[g.bus]] += g.p_nom * g.p_max_pu; end
    for l in values(net.loads);      P[bidx[l.bus]] -= l.p_set;              end

    # Reduce system: remove slack row/column
    # RHS = P_inj + P_shift  (phase-shift transformers add equivalent injections)
    ns    = [i for i in 1:n if i != slack]
    θ_red = B[ns, ns] \ (P[ns] .+ P_shift[ns])
    θ     = zeros(n); θ[ns] = θ_red

    # Active line flows  P_km = b_km · (θ_k − θ_m)
    line_list  = collect(values(net.lines))
    line_flows = [(name = l.name,
                   from = l.from_bus,
                   to   = l.to_bus,
                   P_MW = (net.baseMVA / l.x) * (θ[bidx[l.from_bus]] - θ[bidx[l.to_bus]]),
                   kind = :line)
                  for l in line_list]

    # Transformer flows  P_km = b_eff · (θ_k − θ_m − φ)
    trafo_list  = collect(values(net.transformers))
    trafo_flows = [(name = tr.name,
                    from = tr.from_bus,
                    to   = tr.to_bus,
                    P_MW = (net.baseMVA / (tr.x * tr.tap_ratio)) *
                           (θ[bidx[tr.from_bus]] - θ[bidx[tr.to_bus]]
                            - tr.phase_shift * π / 180),
                    kind = :transformer)
                   for tr in trafo_list]

    all_flows = vcat(line_flows, trafo_flows)

    if verbose
        println("="^55)
        println("DC POWER FLOW — $(net.name)")
        println("="^55)
        @printf("\n  %-14s  %12s\n", "Bus", "θ (rad)")
        println("  " * "─"^29)
        for (i, b) in enumerate(bnames)
            @printf("  %-14s  %12.6f\n", b, θ[i])
        end
        println("\n  Branch active power flows:")
        @printf("  %-16s  %-12s  %10s\n", "Branch", "Type", "P (MW)")
        println("  " * "─"^42)
        for f in all_flows
            @printf("  %-16s  %-12s  %10.2f\n", f.name, string(f.kind), f.P_MW)
        end
        println("="^55)
    end

    return (θ=θ, buses=bnames, line_flows=line_flows, trafo_flows=trafo_flows,
            all_flows=all_flows, B=B, P_inj=P, converged=true)
end

# ────────────────────────────────────────────────────────────────────────────
#  Linearized AC Power Flow
# ────────────────────────────────────────────────────────────────────────────

"""
    linear_ac_pf(net::Network; verbose=true) → NamedTuple

Linearized AC Power Flow (LACPF).

Solves the decoupled system linearised around flat start (|V|=1, θ=0):
    [ B'  -G'] [Δθ ]   [P_inj / S_base]
    [-G'  -B'] [Δ|V|] = [Q_inj / S_base]

Generator reactive injection assumed zero (Q_gen = 0).
Reactive loads taken from `load.q_set`.

Returns:
  (V_mag, V_ang, P_flow, Q_flow, buses, converged)

PyPSA equivalent: `network.pf()` (linearised approximation)
"""
function linear_ac_pf(net::Network; verbose=true)
    bnames = bus_names(net)
    bidx   = bus_index(net)
    n      = length(bnames)
    slack  = _slack_index(net, bidx)

    G_mw, B_mw = _build_Y(net, bidx)
    P_shift    = _phase_shift_injection(net, bidx)

    P = zeros(n); Q = zeros(n)
    for g in values(net.generators); P[bidx[g.bus]] += g.p_nom * g.p_max_pu; end
    for l in values(net.loads)
        P[bidx[l.bus]] -= l.p_set
        Q[bidx[l.bus]] -= l.q_set
    end
    P .+= P_shift   # phase-shift equivalent injections (active power only)

    ns = [i for i in 1:n if i != slack]
    nm = n - 1

    B_r = B_mw[ns, ns]; G_r = G_mw[ns, ns]
    A   = [B_r  -G_r;
           -G_r  -B_r]
    rhs = [P[ns]; Q[ns]]

    x_sol = A \ rhs
    Δθ_r  = x_sol[1:nm];   Δv_r = x_sol[nm+1:end]

    θ   = zeros(n); θ[ns]   = Δθ_r
    Δv  = zeros(n); Δv[ns]  = Δv_r
    V_mag = 1.0 .+ Δv

    line_list = collect(values(net.lines))
    P_flow = Float64[]; Q_flow = Float64[]
    for l in line_list
        f = bidx[l.from_bus]; t = bidx[l.to_bus]
        d    = l.r^2 + l.x^2
        g_ij = l.r / d; b_ij = l.x / d
        Δθ_ft = θ[f] - θ[t]; ΔV_ft = Δv[f] - Δv[t]
        push!(P_flow, (b_ij*Δθ_ft + g_ij*ΔV_ft) * net.baseMVA)
        push!(Q_flow, (b_ij*ΔV_ft - g_ij*Δθ_ft) * net.baseMVA)
    end

    if verbose
        println("="^60)
        println("LINEARIZED AC POWER FLOW — $(net.name)")
        println("="^60)
        @printf("\n  %-14s  %10s  %10s\n", "Bus", "|V| (p.u.)", "θ (rad)")
        println("  " * "─"^38)
        for (i, b) in enumerate(bnames)
            @printf("  %-14s  %10.6f  %10.6f\n", b, V_mag[i], θ[i])
        end
        println("\n  Line flows:")
        @printf("  %-14s  %10s  %10s\n", "Line", "P (MW)", "Q (MVAr)")
        println("  " * "─"^38)
        for (k, l) in enumerate(line_list)
            @printf("  %-14s  %10.4f  %10.4f\n", l.name, P_flow[k], Q_flow[k])
        end
        println("="^60)
    end

    return (V_mag=V_mag, V_ang=θ, P_flow=P_flow, Q_flow=Q_flow,
            buses=bnames, converged=true)
end

# ────────────────────────────────────────────────────────────────────────────
#  Single-period Linear OPF
# ────────────────────────────────────────────────────────────────────────────

"""
    lopf(net::Network; line_capacity=Inf, verbose=true) → NamedTuple

Linear Optimal Power Flow (LOPF):
    min  Σᵢ cᵢ · Pᵢ
    s.t. B·θ = P_inj          (nodal power balance)
         |P_km| ≤ s_nom        (thermal limits, per-line or global)
         p_min_pu·p_nom ≤ Pᵢ ≤ p_max_pu·p_nom

Generators with carrier == "wind" are treated as fixed injection (zero cost).
`line_capacity` is the global fallback when a line has s_nom = Inf.

Returns:
  (θ, P_gen, P_line, total_cost, converged, status)

PyPSA equivalent: `network.optimize()`
"""
function lopf(net::Network; line_capacity=Inf, verbose=true)
    bnames = bus_names(net)
    bidx   = bus_index(net)
    n      = length(bnames)
    ref    = _slack_index(net, bidx)

    B, _, _ = _build_B(net, bidx)
    P_shift   = _phase_shift_injection(net, bidx)

    P_load = zeros(n)
    for l in values(net.loads); P_load[bidx[l.bus]] += l.p_set; end

    disp_gens  = [g for g in values(net.generators) if g.carrier != "wind"]
    wind_inj   = zeros(n)
    for g in values(net.generators)
        g.carrier == "wind" && (wind_inj[bidx[g.bus]] += g.p_nom * g.p_max_pu)
    end
    isempty(disp_gens) && error("No dispatchable generators in network")

    link_units  = collect(values(net.links))
    store_units = collect(values(net.stores))

    model = Model(HiGHS.Optimizer); set_silent(model)

    @variable(model, θ[1:n])
    @variable(model, P_gen[g in disp_gens],
              lower_bound = g.p_nom * g.p_min_pu,
              upper_bound = g.p_nom * g.p_max_pu)

    # Link variables: p_link > 0 → power flows bus0 → bus1
    if !isempty(link_units)
        @variable(model, p_link[l in link_units],
                  lower_bound = l.p_nom * l.p_min_pu,
                  upper_bound = l.p_nom * l.p_max_pu)
    end

    # Store: single net power (p > 0 = inject to grid). Single period → no SoC.
    if !isempty(store_units)
        @variable(model, p_store[s in store_units],
                  lower_bound = -s.p_nom,
                  upper_bound =  s.p_nom)
    end

    @constraint(model, θ[ref] == 0.0)

    for k in 1:n
        gen_inj = sum(P_gen[g] for g in disp_gens if bidx[g.bus] == k;
                      init = AffExpr(0.0))
        link_inj = isempty(link_units) ? 0.0 :
            sum( (l.bus1 == bnames[k] ?  l.efficiency : 0.0) * p_link[l] +
                 (l.bus0 == bnames[k] ? -1.0          : 0.0) * p_link[l]
                 for l in link_units; init = AffExpr(0.0))
        store_inj = isempty(store_units) ? AffExpr(0.0) :
            sum(p_store[s] for s in store_units if bidx[s.bus] == k;
                init = AffExpr(0.0))
        @constraint(model,
            sum(B[k,m]*θ[m] for m in 1:n) ==
            gen_inj + wind_inj[k] + link_inj + store_inj - P_load[k] + P_shift[k])
    end

    line_list = collect(values(net.lines))
    for l in line_list
        lim = isfinite(l.s_nom) ? l.s_nom : line_capacity
        isfinite(lim) || continue
        f = bidx[l.from_bus]; t = bidx[l.to_bus]
        b = net.baseMVA / l.x
        @constraint(model, b*(θ[f]-θ[t]) <=  lim)
        @constraint(model, b*(θ[f]-θ[t]) >= -lim)
    end

    for tr in values(net.transformers)
        lim = isfinite(tr.s_nom) ? tr.s_nom : line_capacity
        isfinite(lim) || continue
        f     = bidx[tr.from_bus]; t = bidx[tr.to_bus]
        b_eff = net.baseMVA / (tr.x * tr.tap_ratio)
        φ     = tr.phase_shift * π / 180
        @constraint(model, b_eff*(θ[f]-θ[t]-φ) <=  lim)
        @constraint(model, b_eff*(θ[f]-θ[t]-φ) >= -lim)
    end

    # GlobalConstraints (e.g. CO₂ cap)
    for gc in values(net.global_constraints)
        isfinite(gc.constant) || continue
        expr = sum(P_gen[g] * get(gc.carrier_weightings, g.carrier, 0.0)
                   for g in disp_gens; init = AffExpr(0.0))
        if gc.sense == "<="
            @constraint(model, expr <= gc.constant)
        elseif gc.sense == ">="
            @constraint(model, expr >= gc.constant)
        else
            @constraint(model, expr == gc.constant)
        end
    end

    obj = sum(g.marginal_cost * P_gen[g] for g in disp_gens)
    if !isempty(link_units)
        obj += sum(l.marginal_cost * p_link[l]
                   for l in link_units if l.marginal_cost > 0; init = AffExpr(0.0))
    end
    if !isempty(store_units)
        obj += sum(s.marginal_cost * p_store[s]
                   for s in store_units if s.marginal_cost > 0; init = AffExpr(0.0))
    end
    @objective(model, Min, obj)

    optimize!(model)
    status    = termination_status(model)
    converged = status ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
    !converged && @warn "LOPF status: $status"

    θ_val      = value.(θ)
    P_gen_val  = Dict(g.name => value(P_gen[g]) for g in disp_gens)
    P_link_val = isempty(link_units) ? Dict{String,Float64}() :
                 Dict(l.name => value(p_link[l]) for l in link_units)
    P_store_val= isempty(store_units) ? Dict{String,Float64}() :
                 Dict(s.name => value(p_store[s]) for s in store_units)
    P_line     = [(net.baseMVA / l.x) *
                  (θ_val[bidx[l.from_bus]] - θ_val[bidx[l.to_bus]])
                  for l in line_list]
    total_cost = objective_value(model)

    if verbose && converged
        println("="^62)
        println("LOPF — $(net.name)")
        println("="^62)
        println("\n  Generator dispatch:")
        @printf("  %-16s  %8s  %10s  %12s\n",
                "Generator", "P (MW)", "P_max (MW)", "Cost (€/MWh)")
        println("  " * "─"^52)
        for g in disp_gens
            @printf("  %-16s  %8.2f  %10.2f  %12.2f\n",
                    g.name, P_gen_val[g.name], g.p_nom*g.p_max_pu, g.marginal_cost)
        end
        if !isempty(link_units)
            println("\n  Link flows:")
            @printf("  %-14s  %8s  %8s\n", "Link", "P0 (MW)", "P1 (MW)")
            println("  " * "─"^34)
            for l in link_units
                p0 = P_link_val[l.name]
                @printf("  %-14s  %8.2f  %8.2f\n", l.name, p0, l.efficiency*p0)
            end
        end
        println("\n  Line flows:")
        @printf("  %-14s  %8s  %10s\n", "Line", "P (MW)", "Loading")
        println("  " * "─"^36)
        for (i, l) in enumerate(line_list)
            lim = isfinite(l.s_nom) ? l.s_nom : (isfinite(line_capacity) ? line_capacity : Inf)
            loading = isfinite(lim) ? @sprintf("%7.1f%%", abs(P_line[i])/lim*100) : "  —"
            @printf("  %-14s  %8.2f  %10s\n", l.name, P_line[i], loading)
        end
        @printf("\n  Total cost: %.2f €/h\n", total_cost)
        println("="^62)
    end

    return (θ=θ_val, P_gen=P_gen_val, P_link=P_link_val, P_store=P_store_val,
            P_line=P_line, total_cost=total_cost, converged=converged, status=status)
end

# ────────────────────────────────────────────────────────────────────────────
#  Multi-period LOPF with storage and wind
# ────────────────────────────────────────────────────────────────────────────

"""
    lopf_multiperiod(net; T, load_profile, wind_profile, line_capacity, verbose)

Multi-period LOPF over T time steps.

Storage SoC dynamics per unit s at time t:
    E(t) = E(t-1) + η_ch·P_ch(t) − P_dis(t)/η_dis

Wind generators (carrier == "wind") output is bounded by:
    P_wind(t) ≤ p_nom × wind_profile[t]   (zero marginal cost)

Returns:
  (gen_dispatch, soc, p_ch, p_dis, total_cost, status)
  All dispatch dicts are indexed by component name.

PyPSA equivalent: `network.optimize()` with multi-period snapshots
"""
function lopf_multiperiod(net::Network;
                           T             = 24,
                           load_profile  = DEFAULT_LOAD_PROFILE,
                           wind_profile  = DEFAULT_WIND_PROFILE,
                           line_capacity = Inf,
                           verbose       = true)
    bnames = bus_names(net)
    bidx   = bus_index(net)
    n      = length(bnames)
    ref    = _slack_index(net, bidx)

    B, _, _ = _build_B(net, bidx)
    P_shift   = _phase_shift_injection(net, bidx)   # time-invariant

    P_load_base = zeros(n)
    for l in values(net.loads); P_load_base[bidx[l.bus]] += l.p_set; end

    disp_gens   = [g for g in values(net.generators) if g.carrier != "wind"]
    wind_gens   = [g for g in values(net.generators) if g.carrier == "wind"]
    stor_units  = collect(values(net.storage_units))
    store_units = collect(values(net.stores))
    link_units  = collect(values(net.links))

    model = Model(HiGHS.Optimizer); set_silent(model)

    # ── Variables ──────────────────────────────────────────────────
    @variable(model, θ[1:n, 1:T])
    @variable(model, P_gen[g in disp_gens, 1:T],
              lower_bound = g.p_nom * g.p_min_pu,
              upper_bound = g.p_nom * g.p_max_pu)

    # Link variables
    if !isempty(link_units)
        @variable(model, p_link[l in link_units, 1:T],
                  lower_bound = l.p_nom * l.p_min_pu,
                  upper_bound = l.p_nom * l.p_max_pu)
    end

    # Store variables: single net power + SoC
    if !isempty(store_units)
        @variable(model, p_store[s in store_units, 1:T],
                  lower_bound = -s.p_nom,
                  upper_bound =  s.p_nom)
        @variable(model, e_store[s in store_units, 0:T],
                  lower_bound = s.e_nom * s.e_min_pu,
                  upper_bound = s.e_nom * s.e_max_pu)
        for s in store_units
            e0 = s.e_initial > 0 ? s.e_initial : s.e_nom * 0.5
            @constraint(model, e_store[s,0] == e0)
            for t in 1:T
                # E(t) = (1-sl)·E(t-1) - p(t)  [p>0 discharges]
                @constraint(model,
                    e_store[s,t] == (1-s.standing_loss)*e_store[s,t-1] - p_store[s,t])
            end
            s.e_cyclic && @constraint(model, e_store[s,T] == e_store[s,0])
        end
    end

    if !isempty(stor_units)
        @variable(model, p_ch[s in stor_units, 1:T],
                  lower_bound=0.0, upper_bound=s.p_nom)
        @variable(model, p_dis[s in stor_units, 1:T],
                  lower_bound=0.0, upper_bound=s.p_nom)
        @variable(model, soc_var[s in stor_units, 0:T],
                  lower_bound=0.0, upper_bound=s.e_nom)

        for s in stor_units
            e0 = s.e_initial > 0 ? s.e_initial : s.e_nom * 0.5
            @constraint(model, soc_var[s, 0] == e0)
            for t in 1:T
                @constraint(model,
                    soc_var[s,t] == soc_var[s,t-1]
                                  + s.efficiency_charge    * p_ch[s,t]
                                  - p_dis[s,t] / s.efficiency_discharge)
            end
            s.cyclic_state_of_charge &&
                @constraint(model, soc_var[s,T] == soc_var[s,0])
        end
    end

    # ── Per-period constraints ──────────────────────────────────────
    line_list = collect(values(net.lines))
    for t in 1:T
        prof_t = load_profile[mod1(t, length(load_profile))]
        wprf_t = wind_profile[mod1(t, length(wind_profile))]

        @constraint(model, θ[ref, t] == 0.0)

        for k in 1:n
            gen_inj = sum(P_gen[g,t] for g in disp_gens if bidx[g.bus]==k;
                          init=AffExpr(0.0))
            wind_inj = sum(g.p_nom * wprf_t
                           for g in wind_gens if bidx[g.bus]==k;
                           init=0.0)
            stor_net = isempty(stor_units) ? AffExpr(0.0) :
                       sum(p_dis[s,t] - p_ch[s,t]
                           for s in stor_units if bidx[s.bus]==k;
                           init=AffExpr(0.0))

            link_inj = isempty(link_units) ? 0.0 :
                sum( (l.bus1 == bnames[k] ?  l.efficiency : 0.0) * p_link[l,t] +
                     (l.bus0 == bnames[k] ? -1.0          : 0.0) * p_link[l,t]
                     for l in link_units; init = AffExpr(0.0))
            store_inj = isempty(store_units) ? AffExpr(0.0) :
                sum(p_store[s,t] for s in store_units if bidx[s.bus]==k;
                    init = AffExpr(0.0))

            @constraint(model,
                sum(B[k,m]*θ[m,t] for m in 1:n) ==
                gen_inj + wind_inj + stor_net + link_inj + store_inj
                - P_load_base[k]*prof_t + P_shift[k])
        end

        for l in line_list
            lim = isfinite(l.s_nom) ? l.s_nom : line_capacity
            isfinite(lim) || continue
            f = bidx[l.from_bus]; tt = bidx[l.to_bus]
            b = net.baseMVA / l.x
            @constraint(model, b*(θ[f,t]-θ[tt,t]) <=  lim)
            @constraint(model, b*(θ[f,t]-θ[tt,t]) >= -lim)
        end

        for tr in values(net.transformers)
            lim = isfinite(tr.s_nom) ? tr.s_nom : line_capacity
            isfinite(lim) || continue
            f     = bidx[tr.from_bus]; tt = bidx[tr.to_bus]
            b_eff = net.baseMVA / (tr.x * tr.tap_ratio)
            φ     = tr.phase_shift * π / 180
            @constraint(model, b_eff*(θ[f,t]-θ[tt,t]-φ) <=  lim)
            @constraint(model, b_eff*(θ[f,t]-θ[tt,t]-φ) >= -lim)
        end
    end

    # ── GlobalConstraints (CO₂ cap etc.) ───────────────────────────
    for gc in values(net.global_constraints)
        isfinite(gc.constant) || continue
        expr = sum(P_gen[g,t] * get(gc.carrier_weightings, g.carrier, 0.0)
                   for g in disp_gens, t in 1:T; init = AffExpr(0.0))
        if gc.sense == "<="
            @constraint(model, expr <= gc.constant)
        elseif gc.sense == ">="
            @constraint(model, expr >= gc.constant)
        else
            @constraint(model, expr == gc.constant)
        end
    end

    # ── Objective ──────────────────────────────────────────────────
    obj = sum(g.marginal_cost * P_gen[g,t] for g in disp_gens, t in 1:T)
    if !isempty(link_units)
        obj += sum(l.marginal_cost * p_link[l,t]
                   for l in link_units, t in 1:T if l.marginal_cost > 0;
                   init = AffExpr(0.0))
    end
    if !isempty(store_units)
        obj += sum(s.marginal_cost * p_store[s,t]
                   for s in store_units, t in 1:T if s.marginal_cost > 0;
                   init = AffExpr(0.0))
    end
    @objective(model, Min, obj)

    optimize!(model)
    status = termination_status(model)
    status != MOI.OPTIMAL && @warn "Multi-period LOPF status: $status"

    gen_dispatch = Dict(g.name => [value(P_gen[g,t]) for t in 1:T]
                        for g in disp_gens)

    soc_out  = isempty(stor_units) ? Dict{String,Vector{Float64}}() :
               Dict(s.name => [value(soc_var[s,t]) for t in 0:T] for s in stor_units)
    pch_out  = isempty(stor_units) ? Dict{String,Vector{Float64}}() :
               Dict(s.name => [value(p_ch[s,t])    for t in 1:T] for s in stor_units)
    pdis_out = isempty(stor_units) ? Dict{String,Vector{Float64}}() :
               Dict(s.name => [value(p_dis[s,t])   for t in 1:T] for s in stor_units)

    store_e_out = isempty(store_units) ? Dict{String,Vector{Float64}}() :
                  Dict(s.name => [value(e_store[s,t]) for t in 0:T] for s in store_units)
    store_p_out = isempty(store_units) ? Dict{String,Vector{Float64}}() :
                  Dict(s.name => [value(p_store[s,t]) for t in 1:T] for s in store_units)
    link_p_out  = isempty(link_units) ? Dict{String,Vector{Float64}}() :
                  Dict(l.name => [value(p_link[l,t]) for t in 1:T] for l in link_units)

    total_cost = objective_value(model)

    if verbose
        println("="^65)
        println("MULTI-PERIOD LOPF ($(T)h) — $(net.name)")
        println("="^65)
        @printf("Total generation cost: %.2f €\n", total_cost)
        println("\n  Generator summary:")
        @printf("  %-16s  %8s  %8s  %8s\n", "Generator", "Avg MW", "Max MW", "Min MW")
        println("  " * "─"^45)
        for (name, vec) in gen_dispatch
            @printf("  %-16s  %8.1f  %8.1f  %8.1f\n",
                    name, sum(vec)/T, maximum(vec), minimum(vec))
        end
        if !isempty(soc_out)
            println("\n  StorageUnit SoC range [MWh]:")
            for (name, vec) in soc_out
                @printf("  %-16s  min=%6.1f  max=%6.1f\n", name, minimum(vec), maximum(vec))
            end
        end
        if !isempty(store_e_out)
            println("\n  Store energy range [MWh]:")
            for (name, vec) in store_e_out
                @printf("  %-16s  min=%6.1f  max=%6.1f\n", name, minimum(vec), maximum(vec))
            end
        end
        if !isempty(link_p_out)
            println("\n  Link avg flow [MW]:")
            for (name, vec) in link_p_out
                @printf("  %-16s  avg=%6.1f\n", name, sum(vec)/T)
            end
        end
        println("="^65)
    end

    return (gen_dispatch=gen_dispatch,
            soc=soc_out, p_ch=pch_out, p_dis=pdis_out,
            store_e=store_e_out, store_p=store_p_out,
            link_p=link_p_out,
            total_cost=total_cost, status=status)
end
