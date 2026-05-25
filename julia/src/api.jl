"""
    api.jl

PyPSA-compatible high-level API.
Provides add!, pf, and optimize as single entry points,
mirroring the PyPSA workflow as closely as possible.
"""

# ─────────────────────────────────────────────────────────────────────────────
#  add!  —  PyPSA-style generic component addition
# ─────────────────────────────────────────────────────────────────────────────

"""
    add!(net, component_type, name; kwargs...)

Generic component addition, mirroring PyPSA's `network.add()`.

Supported component types (String or Symbol):
  "Bus", "Line", "Transformer", "Generator", "Load",
  "StorageUnit", "Store", "Carrier", "Link", "GlobalConstraint"

PyPSA-style aliases accepted:
  - Line:      `bus0`/`bus1`  →  `from_bus`/`to_bus`
  - Link:      `bus0`/`bus1`  →  `bus0`/`bus1` (already matching)
  - Generator: `bus`          →  `bus` (matching)

Examples:
    add!(net, "Bus",       "B1"; v_nom=380.0, slack=true)
    add!(net, "Generator", "G1"; bus="B1", p_nom=400.0, marginal_cost=20.0)
    add!(net, "Line",      "L1"; bus0="B1", bus1="B2", x=0.1, s_nom=200.0)
    add!(net, "Load",      "D1"; bus="B1", p_set=300.0)
    add!(net, "Link",      "HVDC1"; bus0="B1", bus1="B2", p_nom=500.0, efficiency=0.97)
    add!(net, "Carrier",   "gas"; co2_emissions=0.20)
    add!(net, "GlobalConstraint", "co2_cap";
         constant=500.0, carrier_weightings=Dict("gas"=>0.20))

PyPSA equivalent: `network.add("ComponentType", "name", **params)`
"""
function add!(net::Network, ctype::Union{String,Symbol}, name::String; kwargs...)
    d = Dict{Symbol,Any}(kwargs)
    ct = string(ctype)

    if ct == "Bus"
        add_bus!(net, name; d...)

    elseif ct == "Line"
        # Accept bus0/bus1 (PyPSA) or from_bus/to_bus (our native names)
        fb = pop!(d, :bus0, get(d, :from_bus, nothing))
        tb = pop!(d, :bus1, get(d, :to_bus,  nothing))
        delete!(d, :from_bus); delete!(d, :to_bus)
        (fb === nothing || tb === nothing) &&
            error("Line '$name': provide bus0+bus1 or from_bus+to_bus")
        add_line!(net, name, string(fb), string(tb); d...)

    elseif ct == "Transformer"
        fb = pop!(d, :bus0, get(d, :from_bus, nothing))
        tb = pop!(d, :bus1, get(d, :to_bus,  nothing))
        delete!(d, :from_bus); delete!(d, :to_bus)
        (fb === nothing || tb === nothing) &&
            error("Transformer '$name': provide bus0+bus1 or from_bus+to_bus")
        add_transformer!(net, name, string(fb), string(tb); d...)

    elseif ct == "Generator"
        bus = pop!(d, :bus, nothing)
        bus === nothing && error("Generator '$name': provide bus=...")
        add_generator!(net, name, string(bus); d...)

    elseif ct == "Load"
        bus = pop!(d, :bus, nothing)
        bus === nothing && error("Load '$name': provide bus=...")
        add_load!(net, name, string(bus); d...)

    elseif ct == "StorageUnit"
        bus = pop!(d, :bus, nothing)
        bus === nothing && error("StorageUnit '$name': provide bus=...")
        add_storage_unit!(net, name, string(bus); d...)

    elseif ct == "Store"
        bus = pop!(d, :bus, nothing)
        bus === nothing && error("Store '$name': provide bus=...")
        add_store!(net, name, string(bus); d...)

    elseif ct == "Carrier"
        add_carrier!(net, name; d...)

    elseif ct == "Link"
        b0 = pop!(d, :bus0, nothing)
        b1 = pop!(d, :bus1, nothing)
        (b0 === nothing || b1 === nothing) &&
            error("Link '$name': provide bus0 and bus1")
        add_link!(net, name, string(b0), string(b1); d...)

    elseif ct == "GlobalConstraint"
        add_global_constraint!(net, name; d...)

    else
        error("Unknown component type: '$ct'. " *
              "Valid: Bus, Line, Transformer, Generator, Load, " *
              "StorageUnit, Store, Carrier, Link, GlobalConstraint")
    end
    return net
end

# ─────────────────────────────────────────────────────────────────────────────
#  pf  —  Power Flow entry point
# ─────────────────────────────────────────────────────────────────────────────

"""
    pf(net; method=:dc, verbose=true) → NamedTuple

Power flow entry point. Mirrors PyPSA's `network.pf()`.

`method` options:
  :dc   — DC Power Flow (default, fast, lossless)
  :lac  — Linearized AC Power Flow (recovers Q flows and |V| deviations)
  :ac   — Full nonlinear AC Power Flow via PowerModels.jl + Ipopt (Newton-Raphson)
  :auto — :lac if any line has r > 0, otherwise :dc

Returns the same NamedTuple as the underlying solver.

PyPSA equivalent:
  `network.lpf()`  →  pf(net, method=:dc)
  `network.pf()`   →  pf(net, method=:ac)
"""
function pf(net::Network; method=:dc, verbose=true, kwargs...)
    if method == :auto
        has_r = any(l.r > 0 for l in values(net.lines))
        method = has_r ? :lac : :dc
    end
    if method == :dc
        return dc_pf(net; verbose=verbose, kwargs...)
    elseif method == :lac
        return linear_ac_pf(net; verbose=verbose, kwargs...)
    elseif method == :ac
        return ac_pf(net; verbose=verbose, kwargs...)
    else
        error("pf: unknown method '$method'. Use :dc, :lac, :ac, or :auto")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
#  optimize  —  Optimisation entry point
# ─────────────────────────────────────────────────────────────────────────────

"""
    optimize(net; method=:auto, T=1, kwargs...) → NamedTuple

Optimisation entry point. Mirrors PyPSA's `network.optimize()`.

`method` options:
  :auto  — selects based on network:
              T=1  + no committable gens → :lopf
              T>1  + no committable gens → :mp
              any committable gen        → :uc
  :lopf  — single-period Linear OPF
  :mp    — multi-period Linear OPF (24-hour with storage and wind)
  :uc    — Unit Commitment MILP

`T` — number of time periods (default 1; only used for :mp and :uc)

Examples:
    result = optimize(net)                       # single-period LOPF
    result = optimize(net, T=24)                 # 24-hour multi-period LOPF
    result = optimize(net, method=:uc, T=24)     # unit commitment

PyPSA equivalent: `network.optimize()`
"""
function optimize(net::Network; method=:auto, T=1, kwargs...)
    if method == :auto
        has_committable = any(g.committable for g in values(net.generators)
                              if g.carrier != "wind")
        if has_committable
            method = :uc
        elseif T > 1
            method = :mp
        else
            method = :lopf
        end
    end

    if method == :lopf
        return lopf(net; kwargs...)
    elseif method == :mp
        return lopf_multiperiod(net; T=T, kwargs...)
    elseif method == :uc
        return unit_commitment(net; T=T, kwargs...)
    else
        error("optimize: unknown method '$method'. Use :lopf, :mp, :uc, or :auto")
    end
end
