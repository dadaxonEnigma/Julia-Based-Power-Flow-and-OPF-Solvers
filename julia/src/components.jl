"""
    components.jl

Type-stable structs for power network components.
Attribute naming follows PyPSA conventions where applicable.
"""

# ============================================================
#  Bus
# ============================================================
"""
    Bus(name; v_nom, v_mag_pu, v_ang, slack, bus_type)

Network bus (node). Voltage reference is per-unit with base `v_nom`.

PyPSA equivalent: `network.buses`
| Julia field    | PyPSA attribute  | Unit  |
|----------------|------------------|-------|
| v_nom          | v_nom            | kV    |
| v_mag_pu       | v_mag_pu_set     | p.u.  |
| v_ang          | —                | rad   |
| bus_type       | 1=PQ 2=PV 3=slack| —     |
"""
struct Bus
    name::String
    v_nom::Float64       # nominal voltage [kV]
    v_mag_pu::Float64    # initial voltage magnitude [p.u.]
    v_ang::Float64       # initial voltage angle [rad]
    slack::Bool
    bus_type::Int        # 1=PQ, 2=PV, 3=slack (MATPOWER convention)
end

function Bus(name::String;
             v_nom    = 1.0,
             v_mag_pu = 1.0,
             v_ang    = 0.0,
             slack    = false,
             bus_type = slack ? 3 : 1)
    Bus(name, Float64(v_nom), Float64(v_mag_pu),
        Float64(v_ang), slack, bus_type)
end

# ============================================================
#  Line
# ============================================================
"""
    Line(name, from_bus, to_bus; r, x, b, s_nom, length)

AC transmission line. Impedance parameters in per-unit on system base.

PyPSA equivalent: `network.lines`
| Julia field | PyPSA attribute | Unit  |
|-------------|-----------------|-------|
| r           | r               | p.u.  |
| x           | x               | p.u.  |
| b           | b               | p.u.  |
| s_nom       | s_nom           | MVA   |
| length      | length          | km    |
"""
struct Line
    name::String
    from_bus::String
    to_bus::String
    r::Float64           # resistance [p.u.]
    x::Float64           # reactance [p.u.]
    b::Float64           # total shunt susceptance [p.u.]
    s_nom::Float64       # thermal limit [MVA]; Inf = unconstrained
    length::Float64      # geographic length [km]
end

function Line(name::String, from_bus::String, to_bus::String;
              r      = 0.0,
              x      = 0.1,
              b      = 0.0,
              s_nom  = Inf,
              length = 1.0)
    x == 0.0 && error("Line $name: reactance x must be non-zero")
    Line(name, from_bus, to_bus,
         Float64(r), Float64(x), Float64(b),
         Float64(s_nom), Float64(length))
end

# ============================================================
#  Transformer
# ============================================================
"""
    Transformer(name, from_bus, to_bus; r, x, s_nom, tap_ratio, phase_shift)

Two-winding transformer with off-nominal tap and phase shift.
Modelled as a π-equivalent branch with complex turns ratio:
    a = tap_ratio · exp(j · phase_shift · π/180)

PyPSA equivalent: `network.transformers`
| Julia field  | PyPSA attribute | Unit    |
|--------------|-----------------|---------|
| r            | r               | p.u.    |
| x            | x               | p.u.    |
| s_nom        | s_nom           | MVA     |
| tap_ratio    | tap_ratio       | p.u.    |
| phase_shift  | phase_shift     | degrees |
"""
struct Transformer
    name::String
    from_bus::String
    to_bus::String
    r::Float64
    x::Float64
    s_nom::Float64
    tap_ratio::Float64   # off-nominal turns ratio [p.u.], 1.0 = nominal
    phase_shift::Float64 # phase shift [degrees]
end

function Transformer(name::String, from_bus::String, to_bus::String;
                     r           = 0.0,
                     x           = 0.1,
                     s_nom       = Inf,
                     tap_ratio   = 1.0,
                     phase_shift = 0.0)
    Transformer(name, from_bus, to_bus,
                Float64(r), Float64(x), Float64(s_nom),
                Float64(tap_ratio), Float64(phase_shift))
end

# ============================================================
#  Generator
# ============================================================
"""
    Generator(name, bus; p_nom, p_min_pu, p_max_pu, marginal_cost, carrier, committable)

Synchronous or non-synchronous generator.
Dispatch bounds: p_min_pu·p_nom ≤ P ≤ p_max_pu·p_nom.

PyPSA equivalent: `network.generators`
| Julia field    | PyPSA attribute  | Unit    |
|----------------|------------------|---------|
| p_nom          | p_nom            | MW      |
| p_min_pu       | p_min_pu         | p.u.    |
| p_max_pu       | p_max_pu         | p.u.    |
| marginal_cost  | marginal_cost    | €/MWh   |
| carrier        | carrier          | string  |
| committable    | committable      | bool    |
"""
struct Generator
    name::String
    bus::String
    p_nom::Float64
    p_min_pu::Float64    # fraction of p_nom (lower bound)
    p_max_pu::Float64    # fraction of p_nom (upper bound)
    marginal_cost::Float64  # €/MWh
    carrier::String      # "gas", "coal", "wind", "solar", "hydro", ...
    committable::Bool    # true → UC binary variable required
end

function Generator(name::String, bus::String;
                   p_nom         = 0.0,
                   p_min_pu      = 0.0,
                   p_max_pu      = 1.0,
                   marginal_cost = 0.0,
                   carrier       = "gas",
                   committable   = false)
    Generator(name, bus,
              Float64(p_nom), Float64(p_min_pu), Float64(p_max_pu),
              Float64(marginal_cost), carrier, committable)
end

# ============================================================
#  Load
# ============================================================
"""
    Load(name, bus; p_set, q_set)

Fixed demand. In multi-period problems p_set is scaled by an external profile.

PyPSA equivalent: `network.loads`
| Julia field | PyPSA attribute | Unit  |
|-------------|-----------------|-------|
| p_set       | p_set           | MW    |
| q_set       | q_set           | MVAr  |
"""
struct Load
    name::String
    bus::String
    p_set::Float64       # active power demand [MW]
    q_set::Float64       # reactive power demand [MVAr]
end

function Load(name::String, bus::String;
              p_set = 0.0,
              q_set = 0.0)
    Load(name, bus, Float64(p_set), Float64(q_set))
end

# ============================================================
#  StorageUnit
# ============================================================
"""
    StorageUnit(name, bus; p_nom, e_nom, efficiency_charge, efficiency_discharge,
                standing_loss, e_initial, cyclic_state_of_charge)

Battery or pumped-hydro storage with explicit SoC dynamics:
    E(t+1) = (1 − standing_loss)·E(t) + η_ch·P_ch(t) − P_dis(t)/η_dis

PyPSA equivalent: `network.storage_units`
| Julia field              | PyPSA attribute              | Unit  |
|--------------------------|------------------------------|-------|
| p_nom                    | p_nom                        | MW    |
| e_nom                    | e_nom (= p_nom·max_hours)    | MWh   |
| efficiency_charge        | efficiency_store             | p.u.  |
| efficiency_discharge     | efficiency_dispatch          | p.u.  |
| standing_loss            | standing_loss                | 1/h   |
| e_initial                | state_of_charge_initial      | MWh   |
| cyclic_state_of_charge   | cyclic_state_of_charge       | bool  |
"""
struct StorageUnit
    name::String
    bus::String
    p_nom::Float64              # rated charge/discharge power [MW]
    e_nom::Float64              # energy capacity [MWh]
    efficiency_charge::Float64
    efficiency_discharge::Float64
    standing_loss::Float64      # fractional energy loss per hour
    e_initial::Float64          # initial state of charge [MWh]
    cyclic_state_of_charge::Bool
end

function StorageUnit(name::String, bus::String;
                     p_nom                  = 0.0,
                     e_nom                  = 0.0,
                     efficiency_charge      = 0.9,
                     efficiency_discharge   = 0.9,
                     standing_loss          = 0.0,
                     e_initial              = 0.0,
                     cyclic_state_of_charge = true)
    StorageUnit(name, bus,
                Float64(p_nom), Float64(e_nom),
                Float64(efficiency_charge), Float64(efficiency_discharge),
                Float64(standing_loss), Float64(e_initial),
                cyclic_state_of_charge)
end
