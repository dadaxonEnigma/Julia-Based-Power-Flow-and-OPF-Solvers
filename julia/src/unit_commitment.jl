"""
    unit_commitment.jl

Unit Commitment MILP solver. Include this file after dispatch.jl.
"""

"""
    unit_commitment(net; T, load_profile, wind_profile, line_capacity,
                    compute_lmp, verbose) -> NamedTuple

Unit Commitment: Mixed-Integer Linear Program over T time steps.

For generators with `committable = true`:
  - Binary commitment:  u[g,t]  in {0,1}
  - Startup:            su[g,t] in {0,1}  (1 if starts at t)
  - Shutdown:           sd[g,t] in {0,1}  (1 if stops at t)
  - Dispatch bounds:    p_min*u <= P <= p_max*u
  - State transition:   u[t] - u[t-1] = su[t] - sd[t]
  - Min up time:   sum(su[t'] for t' in (t-T_up+1):t) <= u[t]
  - Min down time: sum(sd[t'] for t' in (t-T_dn+1):t) <= 1 - u[t]

Objective: dispatch cost + startup cost + shutdown cost

LMP (price-consistent): after MILP, binaries are fixed and LP is re-solved.

Returns: (u, su, sd, P_gen, P_line, lmp, total_cost, status)

PyPSA equivalent: network.optimize() with committable=True generators
"""
function unit_commitment(net::Network;
                         T             = 24,
                         load_profile  = DEFAULT_LOAD_PROFILE,
                         wind_profile  = DEFAULT_WIND_PROFILE,
                         line_capacity = Inf,
                         compute_lmp   = true,
                         verbose       = true)

    bnames = bus_names(net)
    bidx   = bus_index(net)
    n      = length(bnames)
    ref    = _slack_index(net, bidx)

    B, _, _ = _build_B(net, bidx)
    P_shift = _phase_shift_injection(net, bidx)

    P_load_base = zeros(n)
    for l in values(net.loads); P_load_base[bidx[l.bus]] += l.p_set; end

    commit_gens = [g for g in values(net.generators)
                   if g.committable && g.carrier != "wind"]
    cont_gens   = [g for g in values(net.generators)
                   if !g.committable && g.carrier != "wind"]
    wind_gens   = [g for g in values(net.generators) if g.carrier == "wind"]
    stor_units  = collect(values(net.storage_units))
    link_units  = collect(values(net.links))

    # ── MILP model ───────────────────────────────────────────────────────────
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, P_cont[g in cont_gens, 1:T],
              lower_bound = g.p_nom * g.p_min_pu,
              upper_bound = g.p_nom * g.p_max_pu)
    @variable(model, P_com[g in commit_gens, 1:T],
              lower_bound = 0.0,
              upper_bound = g.p_nom * g.p_max_pu)
    @variable(model, u[g  in commit_gens, 1:T], Bin)
    @variable(model, su[g in commit_gens, 1:T], Bin)
    @variable(model, sd[g in commit_gens, 1:T], Bin)
    @variable(model, theta[1:n, 1:T])

    # Storage
    if !isempty(stor_units)
        @variable(model, p_ch[s in stor_units, 1:T],
                  lower_bound = 0.0, upper_bound = s.p_nom)
        @variable(model, p_dis[s in stor_units, 1:T],
                  lower_bound = 0.0, upper_bound = s.p_nom)
        @variable(model, soc_v[s in stor_units, 0:T],
                  lower_bound = 0.0, upper_bound = s.e_nom)
        for s in stor_units
            e0 = s.e_initial > 0 ? s.e_initial : s.e_nom * 0.5
            @constraint(model, soc_v[s, 0] == e0)
            for t in 1:T
                @constraint(model,
                    soc_v[s, t] == soc_v[s, t-1]
                        + s.efficiency_charge   * p_ch[s, t]
                        - p_dis[s, t] / s.efficiency_discharge)
            end
            s.cyclic_state_of_charge &&
                @constraint(model, soc_v[s, T] == soc_v[s, 0])
        end
    end

    # Links
    if !isempty(link_units)
        @variable(model, p_link[l in link_units, 1:T],
                  lower_bound = l.p_nom * l.p_min_pu,
                  upper_bound = l.p_nom * l.p_max_pu)
    end

    # ── UC constraints ───────────────────────────────────────────────────────
    for g in commit_gens
        u0 = g.initial_status ? 1 : 0
        for t in 1:T
            # Dispatch bounded by commitment status
            @constraint(model, P_com[g, t] <= g.p_nom * g.p_max_pu * u[g, t])
            @constraint(model, P_com[g, t] >= g.p_nom * g.p_min_pu * u[g, t])
            # State transition
            u_prev = (t == 1) ? u0 : u[g, t-1]
            @constraint(model, u[g, t] - u_prev == su[g, t] - sd[g, t])
            @constraint(model, su[g, t] + sd[g, t] <= 1)
            # Minimum up time
            if g.min_up_time > 1
                t0 = max(1, t - g.min_up_time + 1)
                @constraint(model, sum(su[g, tau] for tau in t0:t) <= u[g, t])
            end
            # Minimum down time
            if g.min_down_time > 1
                t0 = max(1, t - g.min_down_time + 1)
                @constraint(model, sum(sd[g, tau] for tau in t0:t) <= 1 - u[g, t])
            end
        end
    end

    # ── Power balance per period ──────────────────────────────────────────────
    line_list   = collect(values(net.lines))
    balance_con = Matrix{Any}(undef, n, T)

    for t in 1:T
        prof_t = load_profile[mod1(t, length(load_profile))]
        wprf_t = wind_profile[mod1(t, length(wind_profile))]

        @constraint(model, theta[ref, t] == 0.0)

        for k in 1:n
            cont_inj = sum(P_cont[g, t] for g in cont_gens if bidx[g.bus] == k;
                           init = AffExpr(0.0))
            com_inj  = sum(P_com[g, t]  for g in commit_gens if bidx[g.bus] == k;
                           init = AffExpr(0.0))
            wind_inj = sum(g.p_nom * wprf_t
                           for g in wind_gens if bidx[g.bus] == k; init = 0.0)
            stor_net = isempty(stor_units) ? AffExpr(0.0) :
                       sum(p_dis[s, t] - p_ch[s, t]
                           for s in stor_units if bidx[s.bus] == k;
                           init = AffExpr(0.0))
            link_inj = isempty(link_units) ? 0.0 :
                       sum((l.bus1 == bnames[k] ?  l.efficiency : 0.0) * p_link[l, t] +
                           (l.bus0 == bnames[k] ? -1.0          : 0.0) * p_link[l, t]
                           for l in link_units; init = AffExpr(0.0))

            balance_con[k, t] = @constraint(model,
                sum(B[k, m] * theta[m, t] for m in 1:n) ==
                cont_inj + com_inj + wind_inj + stor_net + link_inj
                - P_load_base[k] * prof_t + P_shift[k])
        end

        for l in line_list
            lim = isfinite(l.s_nom) ? l.s_nom : line_capacity
            isfinite(lim) || continue
            f  = bidx[l.from_bus]; tt = bidx[l.to_bus]
            b  = net.baseMVA / l.x
            @constraint(model, b * (theta[f, t] - theta[tt, t]) <=  lim)
            @constraint(model, b * (theta[f, t] - theta[tt, t]) >= -lim)
        end
    end

    # GlobalConstraints
    for gc in values(net.global_constraints)
        isfinite(gc.constant) || continue
        all_g = vcat(cont_gens, commit_gens)
        expr  = sum(
            (g in cont_gens ? P_cont[g, t] : P_com[g, t]) *
            get(gc.carrier_weightings, g.carrier, 0.0)
            for g in all_g, t in 1:T; init = AffExpr(0.0))
        if     gc.sense == "<=" ; @constraint(model, expr <= gc.constant)
        elseif gc.sense == ">=" ; @constraint(model, expr >= gc.constant)
        else                    ; @constraint(model, expr == gc.constant)
        end
    end

    # ── Objective ─────────────────────────────────────────────────────────────
    @objective(model, Min,
        sum(g.marginal_cost * P_cont[g, t] for g in cont_gens,   t in 1:T) +
        sum(g.marginal_cost * P_com[g, t]  for g in commit_gens, t in 1:T) +
        sum(g.startup_cost  * su[g, t]     for g in commit_gens, t in 1:T) +
        sum(g.shutdown_cost * sd[g, t]     for g in commit_gens, t in 1:T))

    optimize!(model)
    status = termination_status(model)
    status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) && @warn "UC status: $status"

    # ── Extract MILP solution ─────────────────────────────────────────────────
    u_val  = Dict(g.name => [round(Int, value(u[g, t]))  for t in 1:T]
                  for g in commit_gens)
    su_val = Dict(g.name => [round(Int, value(su[g, t])) for t in 1:T]
                  for g in commit_gens)
    sd_val = Dict(g.name => [round(Int, value(sd[g, t])) for t in 1:T]
                  for g in commit_gens)

    P_gen_val = merge(
        Dict(g.name => [value(P_cont[g, t]) for t in 1:T] for g in cont_gens),
        Dict(g.name => [value(P_com[g, t])  for t in 1:T] for g in commit_gens))

    P_line_val = [[net.baseMVA / l.x *
                   (value(theta[bidx[l.from_bus], t]) - value(theta[bidx[l.to_bus], t]))
                   for l in line_list] for t in 1:T]

    total_cost = objective_value(model)

    # ── LMP: fix binaries, re-solve LP ───────────────────────────────────────
    lmp_t = Dict{String, Vector{Float64}}()
    if compute_lmp && status ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        lp = Model(HiGHS.Optimizer)
        set_silent(lp)

        @variable(lp, Pc[g in cont_gens, 1:T],
                  lower_bound = g.p_nom * g.p_min_pu,
                  upper_bound = g.p_nom * g.p_max_pu)
        @variable(lp, Pm[g in commit_gens, 1:T],
                  lower_bound = 0.0,
                  upper_bound = g.p_nom * g.p_max_pu)
        @variable(lp, th[1:n, 1:T])

        if !isempty(stor_units)
            @variable(lp, pc2[s in stor_units, 1:T],
                      lower_bound = 0.0, upper_bound = s.p_nom)
            @variable(lp, pd2[s in stor_units, 1:T],
                      lower_bound = 0.0, upper_bound = s.p_nom)
            @variable(lp, sc2[s in stor_units, 0:T],
                      lower_bound = 0.0, upper_bound = s.e_nom)
            for s in stor_units
                e0 = s.e_initial > 0 ? s.e_initial : s.e_nom * 0.5
                @constraint(lp, sc2[s, 0] == e0)
                for t in 1:T
                    @constraint(lp,
                        sc2[s, t] == sc2[s, t-1]
                            + s.efficiency_charge   * pc2[s, t]
                            - pd2[s, t] / s.efficiency_discharge)
                end
                s.cyclic_state_of_charge &&
                    @constraint(lp, sc2[s, T] == sc2[s, 0])
            end
        end

        lp_bc = Matrix{Any}(undef, n, T)
        for t in 1:T
            pf = load_profile[mod1(t, length(load_profile))]
            wf = wind_profile[mod1(t, length(wind_profile))]
            @constraint(lp, th[ref, t] == 0.0)

            # Fix binaries: bound Pm by fixed u_val
            for g in commit_gens
                uf = u_val[g.name][t]
                @constraint(lp, Pm[g, t] <= g.p_nom * g.p_max_pu * uf)
                @constraint(lp, Pm[g, t] >= g.p_nom * g.p_min_pu * uf)
            end

            for k in 1:n
                ci = sum(Pc[g, t] for g in cont_gens   if bidx[g.bus] == k; init=AffExpr(0.0))
                mi = sum(Pm[g, t] for g in commit_gens if bidx[g.bus] == k; init=AffExpr(0.0))
                wi = sum(g.p_nom * wf for g in wind_gens if bidx[g.bus] == k; init=0.0)
                si = isempty(stor_units) ? AffExpr(0.0) :
                     sum(pd2[s, t] - pc2[s, t]
                         for s in stor_units if bidx[s.bus] == k; init=AffExpr(0.0))
                lp_bc[k, t] = @constraint(lp,
                    sum(B[k, m] * th[m, t] for m in 1:n) ==
                    ci + mi + wi + si - P_load_base[k] * pf + P_shift[k])
            end

            for l in line_list
                lim = isfinite(l.s_nom) ? l.s_nom : line_capacity
                isfinite(lim) || continue
                f  = bidx[l.from_bus]; tt = bidx[l.to_bus]; b = net.baseMVA / l.x
                @constraint(lp, b * (th[f, t] - th[tt, t]) <=  lim)
                @constraint(lp, b * (th[f, t] - th[tt, t]) >= -lim)
            end
        end

        @objective(lp, Min,
            sum(g.marginal_cost * Pc[g, t] for g in cont_gens,   t in 1:T) +
            sum(g.marginal_cost * Pm[g, t] for g in commit_gens, t in 1:T))
        optimize!(lp)

        if termination_status(lp) ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
            lmp_t = Dict(bnames[k] => [-dual(lp_bc[k, t]) for t in 1:T]
                         for k in 1:n)
        end
    end

    # ── Verbose output ────────────────────────────────────────────────────────
    if verbose && status ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        println("=" ^ 68)
        println("UNIT COMMITMENT ($(T)h) — $(net.name)")
        println("=" ^ 68)
        @printf("Total cost: %.2f €\n\n", total_cost)

        if !isempty(commit_gens)
            println("  Commitment schedule  (1=ON  0=OFF):")
            cols = min(T, 24)
            @printf("  %-14s  %s\n", "Generator",
                    join([@sprintf("%3d", t) for t in 1:cols], ""))
            println("  " * "─" ^ (16 + 3 * cols))
            for g in commit_gens
                sch = join([u_val[g.name][t] == 1 ? "  1" : "  ·"
                            for t in 1:cols], "")
                @printf("  %-14s %s   starts=%d\n",
                        g.name, sch, sum(su_val[g.name]))
            end
        end

        println("\n  Dispatch summary:")
        @printf("  %-14s  %8s  %8s  %8s\n", "Generator", "Avg MW", "Max MW", "Min MW")
        println("  " * "─" ^ 44)
        for name in sort(collect(keys(P_gen_val)))
            v = P_gen_val[name]
            @printf("  %-14s  %8.1f  %8.1f  %8.1f\n",
                    name, sum(v)/T, maximum(v), minimum(v))
        end

        if !isempty(lmp_t)
            println("\n  LMP summary (€/MWh):")
            @printf("  %-14s  %8s  %8s  %8s\n", "Bus", "Avg", "Min", "Max")
            println("  " * "─" ^ 44)
            for b in bnames
                v = lmp_t[b]
                @printf("  %-14s  %8.3f  %8.3f  %8.3f\n",
                        b, sum(v)/T, minimum(v), maximum(v))
            end
        end
        println("=" ^ 68)
    end

    return (u = u_val, su = su_val, sd = sd_val,
            P_gen = P_gen_val, P_line = P_line_val,
            lmp = lmp_t, total_cost = total_cost, status = status)
end
