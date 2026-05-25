"""
    PowerFlowJulia

Julia implementation of power flow and optimal power flow solvers,
providing a PyPSA-compatible API with Julia-native performance.

Quick start:
    include("src/PowerFlowJulia.jl")
    using .PowerFlowJulia

    net = Network(name="My Grid", baseMVA=100.0)
    add!(net, "Bus",       "B1"; v_nom=380.0, slack=true)
    add!(net, "Bus",       "B2"; v_nom=380.0)
    add!(net, "Line",      "L1"; bus0="B1", bus1="B2", x=0.1)
    add!(net, "Generator", "G1"; bus="B1", p_nom=400.0, marginal_cost=20.0)
    add!(net, "Load",      "D1"; bus="B2", p_set=300.0)

    result = pf(net)         # DC power flow
    result = optimize(net)   # LOPF

PyPSA docs reference: https://pypsa.readthedocs.io
"""
module PowerFlowJulia

using Printf   # re-export so users get @printf without extra import

# ── Core source files ─────────────────────────────────────────────────────
include("components.jl")
include("network.jl")
include("dispatch.jl")
include("unit_commitment.jl")
include("ac_pf.jl")
include("api.jl")
include("visualization.jl")

# ── Component types ───────────────────────────────────────────────────────
export Network
export Bus, Line, Transformer, Generator, Load, StorageUnit
export Store, Carrier, Link, GlobalConstraint

# ── Network construction API ──────────────────────────────────────────────
export add!                          # generic PyPSA-style: add!(net, "Bus", ...)
export add_bus!, add_line!, add_transformer!
export add_generator!, add_load!, add_storage_unit!
export add_store!, add_carrier!, add_link!, add_global_constraint!

# ── Query helpers ─────────────────────────────────────────────────────────
export bus_names, bus_index, slack_buses
export generators_at, loads_at, branches_at, links_at
export total_p_nom, total_load

# ── High-level solver API (PyPSA-style) ───────────────────────────────────
export pf          # power flow:  pf(net; method=:dc/:lac/:auto)
export optimize    # optimisation: optimize(net; method=:auto, T=1)

# ── Low-level solvers ─────────────────────────────────────────────────────
export dc_pf, linear_ac_pf, ac_pf
export lopf, lopf_multiperiod
export unit_commitment

# ── Visualisation ─────────────────────────────────────────────────────────
export plot_dispatch, plot_lmp, plot_soc, plot_uc_schedule, plot_network

# ── Default profiles ──────────────────────────────────────────────────────
export DEFAULT_LOAD_PROFILE, DEFAULT_WIND_PROFILE

end # module PowerFlowJulia
