"""
    network.jl

Network container — analogue of `pypsa.Network`.
Aggregates all component dictionaries and exposes `add_*!` mutating functions.
"""

include("components.jl")

# ============================================================
#  Network
# ============================================================
"""
    Network(; name, baseMVA)

Top-level container for a power network.
Components are stored in `Dict{String, T}` indexed by component name.

Usage:
    net = Network(name="IEEE 3-bus", baseMVA=100.0)
    add_bus!(net, "Bus1", v_nom=380.0, slack=true)
    add_line!(net, "L1", "Bus1", "Bus2", x=0.1, s_nom=200.0)
    add_generator!(net, "G1", "Bus1", p_nom=400.0, marginal_cost=20.0)
    add_load!(net, "D1", "Bus3", p_set=300.0)
"""
mutable struct Network
    name::String
    baseMVA::Float64
    buses::Dict{String, Bus}
    lines::Dict{String, Line}
    transformers::Dict{String, Transformer}
    generators::Dict{String, Generator}
    loads::Dict{String, Load}
    storage_units::Dict{String, StorageUnit}
end

function Network(; name = "Network", baseMVA = 100.0)
    Network(
        name,
        Float64(baseMVA),
        Dict{String, Bus}(),
        Dict{String, Line}(),
        Dict{String, Transformer}(),
        Dict{String, Generator}(),
        Dict{String, Load}(),
        Dict{String, StorageUnit}(),
    )
end

# ============================================================
#  add_*! functions
# ============================================================

function add_bus!(net::Network, name::String; kwargs...)
    haskey(net.buses, name) && error("Bus '$name' already exists")
    net.buses[name] = Bus(name; kwargs...)
    return net
end

function add_line!(net::Network, name::String,
                   from_bus::String, to_bus::String; kwargs...)
    haskey(net.lines, name)       && error("Line '$name' already exists")
    haskey(net.buses, from_bus)   || error("from_bus '$from_bus' not in network")
    haskey(net.buses, to_bus)     || error("to_bus '$to_bus' not in network")
    net.lines[name] = Line(name, from_bus, to_bus; kwargs...)
    return net
end

function add_transformer!(net::Network, name::String,
                           from_bus::String, to_bus::String; kwargs...)
    haskey(net.transformers, name)  && error("Transformer '$name' already exists")
    haskey(net.buses, from_bus)     || error("from_bus '$from_bus' not in network")
    haskey(net.buses, to_bus)       || error("to_bus '$to_bus' not in network")
    net.transformers[name] = Transformer(name, from_bus, to_bus; kwargs...)
    return net
end

function add_generator!(net::Network, name::String, bus::String; kwargs...)
    haskey(net.generators, name) && error("Generator '$name' already exists")
    haskey(net.buses, bus)       || error("bus '$bus' not in network")
    net.generators[name] = Generator(name, bus; kwargs...)
    return net
end

function add_load!(net::Network, name::String, bus::String; kwargs...)
    haskey(net.loads, name) && error("Load '$name' already exists")
    haskey(net.buses, bus)  || error("bus '$bus' not in network")
    net.loads[name] = Load(name, bus; kwargs...)
    return net
end

function add_storage_unit!(net::Network, name::String, bus::String; kwargs...)
    haskey(net.storage_units, name) && error("StorageUnit '$name' already exists")
    haskey(net.buses, bus)          || error("bus '$bus' not in network")
    net.storage_units[name] = StorageUnit(name, bus; kwargs...)
    return net
end

# ============================================================
#  Query helpers
# ============================================================

"""Return names of all slack buses."""
slack_buses(net::Network) = [b.name for b in values(net.buses) if b.slack]

"""Return all generators attached to a given bus."""
generators_at(net::Network, bus::String) =
    [g for g in values(net.generators) if g.bus == bus]

"""Return all loads attached to a given bus."""
loads_at(net::Network, bus::String) =
    [l for l in values(net.loads) if l.bus == bus]

"""Return all branches (lines + transformers) incident to a bus."""
function branches_at(net::Network, bus::String)
    lines = [l for l in values(net.lines)
             if l.from_bus == bus || l.to_bus == bus]
    trafos = [t for t in values(net.transformers)
              if t.from_bus == bus || t.to_bus == bus]
    return vcat(lines, trafos)
end

"""Total installed generation capacity [MW]."""
total_p_nom(net::Network) = sum(g.p_nom for g in values(net.generators); init=0.0)

"""Total peak load [MW]."""
total_load(net::Network) = sum(l.p_set for l in values(net.loads); init=0.0)

# ============================================================
#  Bus/component index helpers (used by solvers)
# ============================================================

"""Ordered list of bus names (sorted for reproducibility)."""
bus_names(net::Network) = sort(collect(keys(net.buses)))

"""Dict mapping bus name → integer index (1-based)."""
bus_index(net::Network) = Dict(b => i for (i, b) in enumerate(bus_names(net)))

# ============================================================
#  show
# ============================================================

function Base.show(io::IO, net::Network)
    println(io, "Network: \"$(net.name)\"  (baseMVA = $(net.baseMVA) MVA)")
    println(io, "  Buses        : $(length(net.buses))")
    println(io, "  Lines        : $(length(net.lines))")
    println(io, "  Transformers : $(length(net.transformers))")
    println(io, "  Generators   : $(length(net.generators))")
    println(io, "  Loads        : $(length(net.loads))")
    print(io,   "  Storage units: $(length(net.storage_units))")
end
