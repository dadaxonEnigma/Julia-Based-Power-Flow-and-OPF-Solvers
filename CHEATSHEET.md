# PowerFlowJulia — Cheat Sheet
## PyPSA vs Julia API comparison

---

## 1. Network Creation

| PyPSA (Python) | PowerFlowJulia (Julia) |
|---|---|
| `import pypsa` | `include("julia/src/PowerFlowJulia.jl")` |
| `n = pypsa.Network()` | `using .PowerFlowJulia` |
| `n = pypsa.Network(name="G", override_components=...)` | `net = Network(name="My Grid", baseMVA=100.0)` |

```python
# PyPSA
import pypsa
n = pypsa.Network()
n.set_snapshots(range(24))
```

```julia
# Julia
net = Network(name="My Grid", baseMVA=100.0)
```

---

## 2. Components

### Bus

| Field | PyPSA | Julia | Default | Unit |
|---|---|---|---|---|
| name | positional | positional | — | — |
| v_nom | `v_nom` | `v_nom` | 1.0 | kV |
| slack bus | `n.buses.loc["B1","control"]="Slack"` | `slack=true` | false | — |
| bus type | auto-assigned | `bus_type` | 1=PQ, 3=slack | — |

```python
# PyPSA
n.add("Bus", "B1", v_nom=380)
n.add("Bus", "B2", v_nom=380, control="Slack")
```

```julia
# Julia
add!(net, "Bus", "B1"; v_nom=380.0)
add!(net, "Bus", "B2"; v_nom=380.0, slack=true)
```

---

### Line

| Field | PyPSA | Julia | Default | Unit |
|---|---|---|---|---|
| from/to | `bus0`, `bus1` | `bus0`, `bus1` | — | — |
| resistance | `r` | `r` | 0.0 | p.u. |
| reactance | `x` | `x` | 0.1 | p.u. |
| shunt susceptance | `b` | `b` | 0.0 | p.u. |
| thermal limit | `s_nom` | `s_nom` | Inf | MVA |
| length | `length` | `length` | 1.0 | km |

```python
# PyPSA
n.add("Line", "L12", bus0="B1", bus1="B2", x=0.1, r=0.01, s_nom=200)
```

```julia
# Julia  (bus0/bus1 are PyPSA aliases — both work)
add!(net, "Line", "L12"; bus0="B1", bus1="B2", x=0.1, r=0.01, s_nom=200.0)
```

---

### Transformer

| Field | PyPSA | Julia | Default | Unit |
|---|---|---|---|---|
| from/to | `bus0`, `bus1` | `bus0`, `bus1` | — | — |
| reactance | `x` | `x` | 0.1 | p.u. |
| rated power | `s_nom` | `s_nom` | Inf | MVA |
| tap ratio | `tap_ratio` | `tap_ratio` | 1.0 | p.u. |
| phase shift | `phase_shift` | `phase_shift` | 0.0 | degrees |

```python
# PyPSA
n.add("Transformer", "T1", bus0="B1", bus1="B2", x=0.05, s_nom=500, tap_ratio=1.02)
```

```julia
# Julia
add!(net, "Transformer", "T1"; bus0="B1", bus1="B2", x=0.05, s_nom=500.0, tap_ratio=1.02)
```

---

### Generator

| Field | PyPSA | Julia | Default | Unit |
|---|---|---|---|---|
| bus | `bus` | `bus` | — | — |
| capacity | `p_nom` | `p_nom` | 0.0 | MW |
| min output | `p_min_pu` | `p_min_pu` | 0.0 | p.u. |
| max output | `p_max_pu` | `p_max_pu` | 1.0 | p.u. |
| fuel cost | `marginal_cost` | `marginal_cost` | 0.0 | €/MWh |
| fuel type | `carrier` | `carrier` | "gas" | — |
| UC mode | `committable=True` | `committable=true` | false | — |
| min on time | `min_up_time` | `min_up_time` | 1 | h |
| min off time | `min_down_time` | `min_down_time` | 1 | h |
| startup cost | `start_up_cost` | `startup_cost` | 0.0 | € |
| shutdown cost | `shut_down_cost` | `shutdown_cost` | 0.0 | € |

```python
# PyPSA
n.add("Generator", "G1", bus="B1", p_nom=400, marginal_cost=30, carrier="coal")
n.add("Generator", "G2", bus="B1", p_nom=200, marginal_cost=5,  carrier="wind",
      p_max_pu=[0.8, 0.6, ...])  # time-varying
n.add("Generator", "G3", bus="B2", p_nom=300, marginal_cost=65,
      committable=True, min_up_time=3, start_up_cost=5000)
```

```julia
# Julia
add!(net, "Generator", "G1"; bus="B1", p_nom=400.0, marginal_cost=30.0, carrier="coal")
add!(net, "Generator", "G2"; bus="B1", p_nom=200.0, marginal_cost=5.0,  carrier="wind")
add!(net, "Generator", "G3"; bus="B2", p_nom=300.0, marginal_cost=65.0,
     committable=true, min_up_time=3, startup_cost=5000.0)
# Wind profile is passed at solve time: optimize(net; wind_profile=[0.8,0.6,...])
```

---

### Load

| Field | PyPSA | Julia | Default | Unit |
|---|---|---|---|---|
| bus | `bus` | `bus` | — | — |
| active demand | `p_set` | `p_set` | 0.0 | MW |
| reactive demand | `q_set` | `q_set` | 0.0 | MVAr |

```python
# PyPSA
n.add("Load", "D1", bus="B2", p_set=300)
n.add("Load", "D2", bus="B2", p_set=[300, 320, ...])  # time-varying
```

```julia
# Julia (scalar; time-varying via load_profile at solve time)
add!(net, "Load", "D1"; bus="B2", p_set=300.0)
# Multi-period: optimize(net; T=24, load_profile=[0.6, 0.7, ...])
```

---

### StorageUnit

| Field | PyPSA | Julia | Default | Unit |
|---|---|---|---|---|
| bus | `bus` | `bus` | — | — |
| rated power | `p_nom` | `p_nom` | 0.0 | MW |
| energy capacity | `p_nom * max_hours` | `e_nom` | 0.0 | MWh |
| charge efficiency | `efficiency_store` | `efficiency_charge` | 0.9 | p.u. |
| discharge efficiency | `efficiency_dispatch` | `efficiency_discharge` | 0.9 | p.u. |
| standby loss | `standing_loss` | `standing_loss` | 0.0 | 1/h |
| initial SoC | `state_of_charge_initial` | `e_initial` | 0.0 | MWh |
| cyclic SoC | `cyclic_state_of_charge` | `cyclic_state_of_charge` | true | — |

```python
# PyPSA
n.add("StorageUnit", "Bat", bus="B3",
      p_nom=100, max_hours=4,
      efficiency_store=0.92, efficiency_dispatch=0.92,
      standing_loss=0.001, cyclic_state_of_charge=True)
```

```julia
# Julia
add!(net, "StorageUnit", "Bat"; bus="B3",
     p_nom=100.0, e_nom=400.0,          # e_nom = p_nom * max_hours
     efficiency_charge=0.92, efficiency_discharge=0.92,
     standing_loss=0.001)
```

---

### Store

Simplified storage — **no separate charge/discharge efficiency**, single net power variable.

| Field | PyPSA | Julia | Default | Unit |
|---|---|---|---|---|
| energy capacity | `e_nom` | `e_nom` | 0.0 | MWh |
| power limit | derived | `p_nom` | 0.0 | MW |
| min SoC | `e_min_pu` | `e_min_pu` | 0.0 | p.u. |
| max SoC | `e_max_pu` | `e_max_pu` | 1.0 | p.u. |
| cyclic | `e_cyclic` | `e_cyclic` | true | — |

```python
# PyPSA
n.add("Store", "S1", bus="B1", e_nom=500, e_cyclic=True, standing_loss=0.002)
```

```julia
# Julia
add!(net, "Store", "S1"; bus="B1", e_nom=500.0, e_cyclic=true, standing_loss=0.002)
```

---

### Carrier

Energy type metadata for CO₂ accounting.

| Field | PyPSA | Julia | Default | Unit |
|---|---|---|---|---|
| co2 factor | `co2_emissions` | `co2_emissions` | 0.0 | tCO₂/MWh_th |
| color | `color` | `color` | "#888888" | hex |

```python
# PyPSA
n.add("Carrier", "gas",  co2_emissions=0.20)
n.add("Carrier", "coal", co2_emissions=0.34)
```

```julia
# Julia
add!(net, "Carrier", "gas";  co2_emissions=0.20)
add!(net, "Carrier", "coal"; co2_emissions=0.34)
```

---

### Link

Controllable converter between two buses (HVDC, pump, heat pump).

| Field | PyPSA | Julia | Default | Unit |
|---|---|---|---|---|
| from/to | `bus0`, `bus1` | `bus0`, `bus1` | — | — |
| capacity | `p_nom` | `p_nom` | 0.0 | MW |
| min flow | `p_min_pu` | `p_min_pu` | 0.0 | p.u. |
| max flow | `p_max_pu` | `p_max_pu` | 1.0 | p.u. |
| conversion | `efficiency` | `efficiency` | 1.0 | p.u. |
| bidirectional | `p_min_pu=-1` | `p_min_pu=-1.0` | — | — |

```python
# PyPSA
n.add("Link", "HVDC", bus0="B1", bus1="B2",
      p_nom=500, efficiency=0.97, p_min_pu=-1)  # bidirectional
```

```julia
# Julia
add!(net, "Link", "HVDC"; bus0="B1", bus1="B2",
     p_nom=500.0, efficiency=0.97, p_min_pu=-1.0)
```

---

### GlobalConstraint

System-wide constraint (e.g., CO₂ emission cap).

| Field | PyPSA | Julia | Default |
|---|---|---|---|
| type | `type` | `type` | "co2_limit" |
| direction | `sense` | `sense` | "<=" |
| right-hand side | `constant` | `constant` | Inf |
| carrier weights | `carrier_attribute` | `carrier_weightings` | Dict() |

```python
# PyPSA
n.add("GlobalConstraint", "co2_limit",
      type="primary_energy", carrier_attribute="co2_emissions",
      sense="<=", constant=500.0)
```

```julia
# Julia
add!(net, "GlobalConstraint", "co2_cap";
     constant=500.0,
     carrier_weightings=Dict("gas"=>0.20, "coal"=>0.34))
```

---

## 3. Power Flow Solvers — `pf()`

| Method | Julia call | PyPSA equivalent | Notes |
|---|---|---|---|
| DC PF | `pf(net)` or `pf(net, method=:dc)` | `n.lpf()` | Default; lossless, no Q |
| Linearized AC PF | `pf(net, method=:lac)` | — | Recovers \|V\| deviation and Q |
| Full AC PF | `pf(net, method=:ac)` | `n.pf()` | Nonlinear, Newton-Raphson via Ipopt |
| Auto | `pf(net, method=:auto)` | — | :lac if any r>0, else :dc |

### Return fields — `pf()`

Field names differ between DC and AC — they are listed exactly as the solver
returns them (note `θ` is a Unicode field: access as `r.θ`).

```julia
# DC power flow — pf(net) / pf(net, method=:dc)
r = pf(net)
r.converged          # Bool
r.θ                  # Dict{String,Float64}  — voltage angles [rad]
r.line_flows         # Dict{String,Float64}  — line active flows [MW]
r.trafo_flows        # Dict{String,Float64}  — transformer flows [MW]
r.all_flows          # Dict{String,Float64}  — lines + transformers combined
r.P_inj              # Dict{String,Float64}  — net injection per bus [MW]
r.B                  # Matrix                — susceptance matrix
r.buses              # Vector{String}        — bus name order

# AC power flow — pf(net, method=:ac); linearized AC is method=:lac (same fields)
r = pf(net, method=:ac)
r.converged          # Bool
r.V_mag              # Dict{String,Float64}  — |V| [p.u.]
r.V_ang              # Dict{String,Float64}  — voltage angles [rad]  (NOT r.θ)
r.P_flow             # Dict{String,Float64}  — branch active flows [MW]
r.Q_flow             # Dict{String,Float64}  — reactive flows [MVAr]
r.buses              # Vector{String}
```

```python
# PyPSA equivalent fields:
n.buses_t.v_ang["B1"]        # → r.θ["B1"] (DC) / r.V_ang["B1"] (AC)
n.lines_t.p0["L12"]          # → r.line_flows["L12"] (DC) / r.P_flow["L12"] (AC)
n.buses_t.v_mag_pu["B1"]     # → r.V_mag["B1"]  (AC only)
```

---

## 4. Optimisation Solvers — `optimize()`

| Method | Julia call | PyPSA equivalent | Notes |
|---|---|---|---|
| Single LOPF | `optimize(net)` | `n.optimize()` | T=1, LP via HiGHS |
| Multi-period LOPF | `optimize(net, T=24)` | `n.optimize(snapshots=...)` | 24h + storage + wind |
| Unit Commitment | `optimize(net, method=:uc, T=24)` | `n.optimize()` with committable | MILP, binary on/off |
| Stochastic LOPF | `optimize(net, method=:stochastic, T=24, load_scenarios=S)` | — | SAA over S scenarios |
| Auto | `optimize(net, T=24)` | `n.optimize()` | Selects :lopf / :mp / :uc automatically |

### Return fields — single-period LOPF — `optimize(net)`

Scalars per generator/bus/line (no time dimension). NB: line flows are `P_line`,
not `P_flow`; there is **no** `soc`/`gen_dispatch` here (those are multi-period).

```julia
r = optimize(net)

r.status             # MOI status
r.converged          # Bool
r.total_cost         # Float64  — objective value [€]
r.P_gen              # Dict{String,Float64}  — generator dispatch [MW]
r.lmp                # Dict{String,Float64}  — LMP per bus [€/MWh]
r.P_line             # Dict{String,Float64}  — line flows [MW]
r.P_link             # Dict{String,Float64}  — link flows [MW]
r.P_store            # Dict{String,Float64}  — store net power [MW]
r.θ                  # Dict{String,Float64}  — voltage angles [rad]
```

### Return fields — multi-period LOPF — `optimize(net, T=24)`

Time series per component (`Vector` of length T). NB: generator dispatch is
`gen_dispatch` here (**not** `P_gen`); there is no line-flow field in this return.

```julia
r = optimize(net, T=24)

r.status             # MOI status
r.total_cost         # Float64  — objective value [€]
r.gen_dispatch       # Dict{String,Vector{Float64}}  — dispatch per hour [MW]
r.lmp                # Dict{String,Vector{Float64}}  — LMP per bus per hour [€/MWh]
r.soc                # Dict{String,Vector{Float64}}  — storage SoC [MWh]
r.p_ch / r.p_dis     # Dict{String,Vector{Float64}}  — storage charge / discharge [MW]
r.store_e / r.store_p# Dict{String,Vector{Float64}}  — Store energy / power
r.link_p             # Dict{String,Vector{Float64}}  — link flows [MW]
```

```python
# PyPSA equivalent:
n.generators_t.p["G1"]           # → r.P_gen["G1"] (1p) / r.gen_dispatch["G1"] (mp)
n.buses_t.marginal_price["B1"]   # → r.lmp["B1"]
n.lines_t.p0["L12"]              # → r.P_line["L12"] (1p only)
n.storage_units_t.state_of_charge["Bat"]  # → r.soc["Bat"] (mp only)
```

### Return fields — `optimize(method=:uc)`

NB: commitment/startup/shutdown are `u` / `su` / `sd` (short names).

```julia
r = optimize(net, method=:uc, T=24)

r.status             # MOI status
r.total_cost         # Float64
r.P_gen              # Dict{String,Vector{Float64}}  — dispatch per hour [MW]
r.P_line             # Dict{String,Vector{Float64}}  — line flows [MW]
r.u                  # Dict{String,Vector{Float64}}  — commitment u(t) ∈ {0,1}
r.su                 # Dict{String,Vector{Float64}}  — startup events
r.sd                 # Dict{String,Vector{Float64}}  — shutdown events
r.lmp                # Dict{String,Vector{Float64}}  — from LP relaxation
```

### Return fields — `optimize(method=:stochastic)`

```julia
scenarios = predict_scenarios(fc, history; n_scenarios=7)
r = optimize(net, method=:stochastic, T=24,
             load_scenarios=scenarios.scenarios)

r.expected_cost      # Float64  — E[cost] = weighted mean [€]
r.cvar_90            # Float64  — CVaR at 90%: expected cost in worst 10% [€]
r.worst_cost         # Float64  — max scenario cost [€]
r.best_cost          # Float64  — min scenario cost [€]
r.std_cost           # Float64  — std deviation [€]
r.costs              # Vector{Float64}  — per-scenario cost (length S)
r.mean_lmp           # Dict{String,Float64}  — time+scenario average LMP
r.n_feasible         # Int  — feasible scenario count
r.status             # "OK" or "PARTIAL (k/S feasible)"
```

---

## 5. AI Component — `forecasting.jl`

```julia
# Step 1 — generate / load historical data
data = generate_synthetic_data(365; noise_std=0.05, seed=42)
# returns Matrix{Float64}(365 × 24) — daily load profiles [p.u.]

# Step 2 — train LSTM with conformal calibration
fc = train_forecaster(data;
     hidden=32, epochs=100, lr=1e-3,
     val_frac=0.15, cal_frac=0.10, verbose=true)
# returns LoadForecaster struct

# Step 3 — predict next 24h with uncertainty
pred = predict_scenarios(fc, last_24h_observations;
       n_scenarios=7, α=0.10)
# pred.mean      → point forecast [p.u.]  (length 24)
# pred.lower     → conformal lower bound (90% coverage)
# pred.upper     → conformal upper bound
# pred.scenarios → Matrix{Float64}(24 × 7) — 7 load scenarios

# Step 4 — evaluate forecast quality
m = forecast_metrics(actual_values, pred.mean)
# m.mae, m.rmse, m.mape
```

**PyPSA equivalent**: no built-in ML — you would pre-generate profiles externally and pass as `p_set` time series.

---

## 6. Visualization

```julia
# All plots return a Plots.Plot object; pass savefig=true to save automatically.

plot_dispatch(result, net)              # stacked area: generation dispatch
plot_lmp(result, net)                  # heatmap: LMP[bus × hour]
plot_soc(result, net)                  # line: storage SoC over time
plot_uc_schedule(result, net)          # Gantt chart: on/off commitment
plot_network(net)                      # graph: bus/line topology
```

---

## 7. Quick Example — Full Pipeline

```julia
include("julia/src/PowerFlowJulia.jl")
using .PowerFlowJulia

net = Network(baseMVA=100.0)
add!(net, "Bus",       "B1"; v_nom=380.0, slack=true)
add!(net, "Bus",       "B2"; v_nom=380.0)
add!(net, "Line",      "L12"; bus0="B1", bus1="B2", x=0.1, s_nom=200.0)
add!(net, "Generator", "Coal"; bus="B1", p_nom=300.0, marginal_cost=30.0)
add!(net, "Generator", "Gas";  bus="B1", p_nom=200.0, marginal_cost=65.0)
add!(net, "Generator", "Wind"; bus="B2", p_nom=150.0, marginal_cost=5.0, carrier="wind")
add!(net, "Load",      "D1";  bus="B2", p_set=350.0)
add!(net, "StorageUnit","Bat"; bus="B2", p_nom=80.0, e_nom=320.0,
     efficiency_charge=0.92, efficiency_discharge=0.92)

# Power flow
r_dc = pf(net)                            # DC power flow
r_ac = pf(net, method=:ac)               # Full AC power flow

# Optimisation
r1 = optimize(net)                        # single-period LOPF
r2 = optimize(net, T=24)                  # 24h multi-period LOPF
r3 = optimize(net, method=:uc, T=24)      # unit commitment

# AI pipeline
data = generate_synthetic_data(365)
fc   = train_forecaster(data; verbose=false)
pred = predict_scenarios(fc, data[end,:]; n_scenarios=7)
r4   = optimize(net, method=:stochastic, T=24, load_scenarios=pred.scenarios)

println("E[cost] = ", r4.expected_cost, " €")
println("CVaR90  = ", r4.cvar_90, " €")
```

---

## 8. What is NOT implemented

| Feature | PyPSA | PowerFlowJulia | Priority |
|---|---|---|---|
| Capacity expansion (LCOE) | `n.optimize()` with `p_nom_extendable=True` | ❌ not implemented | low |
| N-1 security (SCOPF) | contingency analysis tools | ❌ not implemented | low |
| ShuntImpedance | `network.shunt_impedances` | ❌ not implemented | low |
| DCLine | `network.lines` with DC flag | ❌ not implemented | low |
| Multi-node Links (bus2, bus3…) | `network.links` multi-bus | ❌ single output bus only | low |
| Time-varying `p_max_pu` per generator | `generators_t.p_max_pu` | via `wind_profile` kwarg (global) | partial |
