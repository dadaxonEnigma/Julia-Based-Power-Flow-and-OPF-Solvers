using JuMP
using HiGHS
using Printf
using Statistics

# ----------------------------------------------------------------
#  24-hour load profile (normalised, peak = 1.0)
# ----------------------------------------------------------------
const LOAD_PROFILE = [
    0.60, 0.57, 0.55, 0.54, 0.55, 0.60,   # 00-05  night
    0.70, 0.80, 0.88, 0.90, 0.92, 0.91,   # 06-11  morning ramp
    0.90, 0.89, 0.88, 0.87, 0.89, 0.95,   # 12-17  midday + pre-peak
    1.00, 0.98, 0.93, 0.85, 0.75, 0.65    # 18-23  evening peak + decay
]

# 24-hour wind capacity factor (fraction of installed capacity)
const WIND_PROFILE = [
    0.80, 0.82, 0.85, 0.83, 0.78, 0.70,   # 00-05  strong night wind
    0.60, 0.55, 0.50, 0.45, 0.42, 0.40,   # 06-11  morning drop
    0.38, 0.37, 0.40, 0.43, 0.50, 0.58,   # 12-17  afternoon recovery
    0.65, 0.70, 0.74, 0.76, 0.78, 0.80    # 18-23  evening wind
]

# ================================================================
#  Core solver
# ================================================================
"""
    solve_lopf_multiperiod(buses, lines, generators, loads,
                           storage, wind_buses; kwargs)

Multi-period Linear OPF with StorageUnit and wind generation.

# Arguments
- `buses`      : vector of bus names
- `lines`      : vector of (from, to, r, x) tuples [x in p.u.]
- `generators` : Dict(bus => (P_max_MW, cost_€_per_MWh))
- `loads`      : Dict(bus => P_base_MW)   (scaled by LOAD_PROFILE)
- `storage`    : Dict(bus => (P_nom_MW, E_nom_MWh, η_ch, η_dis))
- `wind_buses` : Dict(bus => P_nom_MW)    (scaled by WIND_PROFILE)
- `T`          : number of time periods (default 24)
- `line_capacity` : thermal limit MW (default Inf)
- `baseMVA`    : system base (default 100.0)
"""
function solve_lopf_multiperiod(
        buses, lines, generators, loads,
        storage  = Dict{Int,NTuple{4,Float64}}(),
        wind_buses = Dict{Int,Float64}();
        T            = 24,
        line_capacity = Inf,
        baseMVA      = 100.0,
        verbose      = true)

    n   = length(buses)
    T_h = T                     # number of hours

    # ── Susceptance matrix ──────────────────────────────────────
    B   = zeros(n, n)
    sus = Float64[]
    for (from, to, r, x) in lines
        b = baseMVA / x
        push!(sus, b)
        B[from,from] += b;  B[to,to]   += b
        B[from,to]   -= b;  B[to,from] -= b
    end

    # ── Base load vector ────────────────────────────────────────
    P_load_base = zeros(n)
    for (bus, p) in loads;  P_load_base[bus] = p;  end

    # ── Generator buses & parameters ───────────────────────────
    gen_buses = sort(collect(keys(generators)))
    P_max_gen = Dict(bus => generators[bus][1] for bus in gen_buses)
    costs     = Dict(bus => generators[bus][2] for bus in gen_buses)
    ref_bus   = gen_buses[1]

    # ── Wind generation (deterministic) ─────────────────────────
    # Tile 24-h profile cyclically for T > 24
    P_wind = Dict{Int,Vector{Float64}}()
    for (bus, p_nom) in wind_buses
        P_wind[bus] = [p_nom * WIND_PROFILE[mod1(t, 24)] for t in 1:T_h]
    end

    # ── JuMP model ───────────────────────────────────────────────
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    # Decision variables
    @variable(model, θ[1:n, 1:T_h])
    @variable(model, P_gen[bus in gen_buses, 1:T_h],
              lower_bound = 0.0,
              upper_bound = P_max_gen[bus])

    # Storage variables
    stor_buses = sort(collect(keys(storage)))
    if !isempty(stor_buses)
        stor_pnom  = Dict(bus => storage[bus][1] for bus in stor_buses)
        stor_enom  = Dict(bus => storage[bus][2] for bus in stor_buses)
        stor_η_ch  = Dict(bus => storage[bus][3] for bus in stor_buses)
        stor_η_dis = Dict(bus => storage[bus][4] for bus in stor_buses)

        @variable(model, p_ch[bus in stor_buses, 1:T_h],
                  lower_bound = 0.0,
                  upper_bound = stor_pnom[bus])
        @variable(model, p_dis[bus in stor_buses, 1:T_h],
                  lower_bound = 0.0,
                  upper_bound = stor_pnom[bus])
        @variable(model, soc[bus in stor_buses, 0:T_h],
                  lower_bound = 0.0,
                  upper_bound = stor_enom[bus])

        # Initial SOC = 50% capacity
        for bus in stor_buses
            @constraint(model, soc[bus, 0] == stor_enom[bus] * 0.5)
        end

        # SOC dynamics:  soc[t] = soc[t-1] + η_ch*p_ch[t] - p_dis[t]/η_dis
        for bus in stor_buses, t in 1:T_h
            @constraint(model,
                soc[bus, t] == soc[bus, t-1]
                             + stor_η_ch[bus]  * p_ch[bus, t]
                             - p_dis[bus, t]   / stor_η_dis[bus])
        end

        # Cyclicity: SOC at end = SOC at start
        for bus in stor_buses
            @constraint(model, soc[bus, T_h] == soc[bus, 0])
        end
    end

    # ── Constraints ──────────────────────────────────────────────
    for t in 1:T_h
        # Reference angle
        @constraint(model, θ[ref_bus, t] == 0.0)

        # Power balance at each bus (tile profile cyclically for T > 24)
        P_load_t = P_load_base .* LOAD_PROFILE[mod1(t, 24)]

        for k in 1:n
            # Net injection: generators + wind - load - storage
            P_inj  = (k in gen_buses)  ? P_gen[k, t] : 0.0
            P_wind_k = haskey(P_wind, k) ? P_wind[k][t] : 0.0

            P_stor_net = 0.0
            if k in stor_buses
                P_stor_net = p_dis[k, t] - p_ch[k, t]
            end

            @constraint(model,
                sum(B[k, m] * θ[m, t] for m in 1:n) ==
                P_inj + P_wind_k + P_stor_net - P_load_t[k])
        end

        # Line thermal limits
        if isfinite(line_capacity)
            for (i, (from, to, r, x)) in enumerate(lines)
                f = sus[i] * (θ[from, t] - θ[to, t])
                @constraint(model, -line_capacity <= f <= line_capacity)
            end
        end
    end

    # ── Objective: minimise total generation cost ────────────────
    @objective(model, Min,
        sum(costs[bus] * P_gen[bus, t]
            for bus in gen_buses, t in 1:T_h))

    optimize!(model)

    status = termination_status(model)
    status != OPTIMAL && @warn "LOPF status: $status"

    # ── Extract results ──────────────────────────────────────────
    gen_dispatch = Dict(
        bus => [value(P_gen[bus, t]) for t in 1:T_h]
        for bus in gen_buses)

    stor_soc = isempty(stor_buses) ? Dict() :
        Dict(bus => [value(soc[bus, t]) for t in 0:T_h]
             for bus in stor_buses)

    stor_ch  = isempty(stor_buses) ? Dict() :
        Dict(bus => [value(p_ch[bus, t]) for t in 1:T_h]
             for bus in stor_buses)

    stor_dis = isempty(stor_buses) ? Dict() :
        Dict(bus => [value(p_dis[bus, t]) for t in 1:T_h]
             for bus in stor_buses)

    total_cost = objective_value(model)

    return (
        status       = status,
        gen_dispatch = gen_dispatch,
        soc          = stor_soc,
        p_ch         = stor_ch,
        p_dis        = stor_dis,
        total_cost   = total_cost,
        wind_gen     = P_wind,
        load_profile = [P_load_base .* LOAD_PROFILE[mod1(t, 24)] for t in 1:T_h]
    )
end

# ================================================================
#  Main: 3-bus test network with storage + wind
# ================================================================
if abspath(PROGRAM_FILE) == @__FILE__

println("="^65)
println("MULTI-PERIOD LOPF  (24h, Storage + Wind)")
println("="^65)

buses = ["Bus 1", "Bus 2", "Bus 3"]
lines = [(1, 2, 0.01, 0.1),
         (1, 3, 0.01, 0.1),
         (2, 3, 0.01, 0.1)]

# G1 deliberately limited to 270 MW so peak demand forces G2 / storage
generators = Dict(
    1 => (270.0, 20.0),   # 270 MW max, 20 €/MWh — cheap baseload
    2 => (100.0, 50.0)    # 100 MW max, 50 €/MWh — expensive backup
)

# Base loads (MW) — scaled hourly by LOAD_PROFILE; total peak = 425 MW
loads = Dict(2 => 250.0, 3 => 175.0)

# Storage at Bus 2: 100 MW / 300 MWh, η=0.95 charge & discharge
storage = Dict(2 => (100.0, 300.0, 0.95, 0.95))

# Wind at Bus 3: 150 MW installed capacity, zero cost
wind_buses = Dict(3 => 150.0)

println("\nNetwork: 3 buses, 3 lines")
println("Generators: G1 (270 MW, 20 €/MWh), G2 (100 MW, 50 €/MWh)")
println("Storage:    Bus 2 — 100 MW / 300 MWh, η=0.95")
println("Wind:       Bus 3 — 150 MW installed")
println("Peak load:  425 MW  →  G1 constrained, forces storage/G2 dispatch\n")

res = solve_lopf_multiperiod(buses, lines, generators, loads,
                              storage, wind_buses)

println("Termination status: $(res.status)")
@printf("Total generation cost (24h): %.2f €\n", res.total_cost)
@printf("Average cost per MWh:        %.2f €/MWh\n",
        res.total_cost / sum(sum(res.load_profile[t]) for t in 1:24))

println("\nHour |  Load (MW) | Wind (MW) | G1 (MW) | G2 (MW) | Stor (MW) | SOC (MWh)")
println("-"^80)
for t in 1:24
    load_t  = sum(res.load_profile[t])
    wind_t  = res.wind_gen[3][t]
    g1_t    = res.gen_dispatch[1][t]
    g2_t    = res.gen_dispatch[2][t]
    ch_t    = res.p_ch[2][t]
    dis_t   = res.p_dis[2][t]
    net_s   = dis_t - ch_t
    soc_t   = res.soc[2][t+1]   # SOC at end of hour t (after dispatch)
    @printf("%3d  |  %8.1f  |  %7.1f  | %7.1f | %7.1f | %9.1f | %8.1f\n",
            t-1, load_t, wind_t, g1_t, g2_t, net_s, soc_t)
end

println("\n" * "="^65)
println("Done. Run python/lopf_multiperiod.py for PyPSA comparison.")
println("="^65)

end # if abspath
