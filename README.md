# PowerFlowJulia

**AI-Assisted Simulation and Optimization of Power Market Projects**
Master's Thesis — Innopolis University, 2026
Author: Dadakhon Turgunboev · Supervisor: Leonard Johard

---

## Overview

PowerFlowJulia is a Julia implementation of core power system simulation and optimization algorithms, reimplementing the key functionality of [PyPSA](https://pypsa.org/) with a PyPSA-compatible API. The library achieves **7×–750× speedups** over PyPSA across all solver methods and integrates an **LSTM-based load forecaster** with conformal prediction intervals for uncertainty-aware stochastic dispatch.

---

## Key Results

| Method | Network size | Julia (ms) | PyPSA (ms) | Speedup |
|---|---|---|---|---|
| DC Power Flow | 500 buses | 10.6 | 811.8 | **77×** |
| LOPF | 500 buses | 21.0 | 15,857 | **753×** |
| AC Power Flow | 100 buses | 47.6 | 341.0 | **7×** |
| Unit Commitment | 14 buses, T=24 | 63.3 | 2,066 | **33×** |
| Multi-Period LOPF | T=24 | 6.4 | 1,485 | **188×** |

**AI component:** LSTM load forecaster (MAPE 13.5%) reduces dispatch cost by **7.9%** vs naive profile. Stochastic LOPF provides CVaR₉₀ = 349,879 € for risk quantification.

---

## Implemented Methods

| Module | Description | Solver |
|---|---|---|
| DC Power Flow | Linearized lossless PF: **B·θ = P** | LAPACK |
| Linearized AC PF | First-order Taylor expansion, recovers \|V\| and Q | LAPACK |
| AC Power Flow | Full nonlinear Newton-Raphson via PowerModels.jl | Ipopt |
| Single-period LOPF | Economic dispatch + DC network constraints | HiGHS LP |
| Multi-period LOPF | 24h dispatch with storage SoC dynamics and wind | HiGHS LP |
| Unit Commitment | MILP with binary on/off, min up/down time, startup costs | HiGHS MILP |
| **LSTM Forecaster** | Sequence-to-sequence load forecast, window=24h → horizon=24h | Flux.jl |
| **Conformal Prediction** | Distribution-free 90% coverage intervals | — |
| **Stochastic LOPF** | Sample Average Approximation over S demand scenarios | HiGHS LP × S |

---

## Project Structure

```
thesis/
├── julia/
│   ├── src/
│   │   ├── PowerFlowJulia.jl     ← main module entry point
│   │   ├── components.jl         ← Bus, Line, Generator, Load, StorageUnit, ...
│   │   ├── network.jl            ← Network container + add_*! functions
│   │   ├── dispatch.jl           ← dc_pf, linear_ac_pf, lopf, lopf_multiperiod
│   │   ├── unit_commitment.jl    ← UC (MILP)
│   │   ├── ac_pf.jl              ← AC PF via PowerModels.jl + Ipopt
│   │   ├── forecasting.jl        ← LSTM forecaster + conformal prediction
│   │   ├── stochastic_lopf.jl    ← scenario-based SAA stochastic LOPF
│   │   ├── api.jl                ← add!(), pf(), optimize() — PyPSA-compatible
│   │   └── visualization.jl      ← plot_dispatch, plot_lmp, plot_soc, ...
│   ├── examples/
│   │   ├── quickstart.jl         ← 3-bus tutorial (mirrors PyPSA quickstart)
│   │   ├── ml_forecast_demo.jl   ← full AI pipeline: LSTM → stochastic LOPF
│   │   ├── unit_commitment_demo.jl
│   │   └── visualization_demo.jl
│   ├── benchmarks/
│   │   ├── benchmark.jl          ← DC PF + LOPF scaling
│   │   ├── benchmark_ac.jl       ← AC PF scaling
│   │   ├── benchmark_uc.jl       ← UC by horizon T
│   │   ├── benchmark_uc_scale.jl ← UC by network size
│   │   ├── benchmark_multiperiod.jl
│   │   └── compare_all_methods.jl ← 6-panel Julia vs PyPSA figure
│   ├── test/
│   │   └── runtests.jl           ← 43 tests (36 core + 7 ML/stochastic)
│   ├── validation/               ← IEEE 14-bus vs MATPOWER
│   └── Project.toml
├── python/                       ← PyPSA reference implementations + benchmarks
├── data/
│   └── case14.m                  ← IEEE 14-bus MATPOWER case
├── results/                      ← benchmark CSVs and PNG figures
├── latex/                        ← Master's thesis (LaTeX)
└── CHEATSHEET.md                 ← PyPSA vs Julia API quick reference
```

---

## Quick Start

```julia
# Install dependencies
# cd julia/ && julia -e 'using Pkg; Pkg.instantiate()'

include("julia/src/PowerFlowJulia.jl")
using .PowerFlowJulia

# Build network
net = Network(baseMVA=100.0)
add!(net, "Bus",       "B1"; v_nom=380.0, slack=true)
add!(net, "Bus",       "B2"; v_nom=380.0)
add!(net, "Line",      "L12"; bus0="B1", bus1="B2", x=0.1, s_nom=200.0)
add!(net, "Generator", "G1"; bus="B1", p_nom=400.0, marginal_cost=20.0)
add!(net, "Generator", "G2"; bus="B2", p_nom=200.0, marginal_cost=50.0)
add!(net, "Load",      "D1"; bus="B2", p_set=300.0)

# Power flow
r = pf(net)                          # DC power flow
r = pf(net, method=:ac)              # Full AC power flow

# Optimisation
r = optimize(net)                    # single-period LOPF
r = optimize(net, T=24)              # 24h multi-period LOPF
r = optimize(net, method=:uc, T=24)  # unit commitment

println("Cost: ", r.total_cost, " €")
println("LMP B1: ", r.lmp["B1"], " €/MWh")
```

---

## AI Pipeline

```julia
# 1. Generate / load historical load data
data = generate_synthetic_data(365; noise_std=0.05, seed=42)

# 2. Train LSTM forecaster with conformal calibration
fc = train_forecaster(data; hidden=32, epochs=100, verbose=true)

# 3. Predict next 24h: point forecast + 90% intervals + 7 scenarios
pred = predict_scenarios(fc, data[end,:]; n_scenarios=7, α=0.10)
# pred.mean      — LSTM point forecast
# pred.lower/upper — conformal bounds (≥90% coverage guaranteed)
# pred.scenarios — 24×7 matrix of sampled profiles

# 4. Stochastic LOPF over scenarios
r = optimize(net, method=:stochastic, T=24, load_scenarios=pred.scenarios)
println("E[cost] = ", r.expected_cost, " €")
println("CVaR90  = ", r.cvar_90, " €")
```

---

## Running Tests

```bash
julia julia/test/runtests.jl
# 43 / 43 passed  (~90 seconds, includes LSTM training)
```

## Running Benchmarks

```bash
# All 6 methods → results/benchmark_all_methods.png
julia julia/benchmarks/compare_all_methods.jl

# Full ML demo → results/plots/ml_*.png
julia julia/examples/ml_forecast_demo.jl
```

---

## Requirements

**Julia ≥ 1.10**

```bash
cd julia
julia -e 'using Pkg; Pkg.instantiate()'
```

Packages: `JuMP`, `HiGHS`, `Ipopt`, `PowerModels`, `Flux`, `Plots`, `StatsPlots`, `CSV`, `DataFrames`

**Python ≥ 3.9** (reference benchmarks only)

```bash
pip install pypsa highspy pandas numpy
```

---

## Validation

| Test case | Method | Julia vs reference | Max error |
|---|---|---|---|
| 3-bus | DC PF | vs PyPSA | < 10⁻¹⁰ rad |
| 3-bus | LOPF | vs PyPSA | 0 € (exact) |
| 3-bus | AC PF | vs PyPSA | < 5×10⁻⁷ p.u. |
| IEEE 14-bus | AC PF | vs MATPOWER | 1.33×10⁻³ p.u. |
| IEEE 14-bus | DC PF | vs PyPSA | identical |
| IEEE 14-bus | LOPF | vs PyPSA | 0 € (exact) |

---

## PyPSA Compatibility

All component and solver argument names match PyPSA where applicable:

| PyPSA | PowerFlowJulia |
|---|---|
| `n = pypsa.Network()` | `net = Network()` |
| `n.add("Bus", "B1", v_nom=380)` | `add!(net, "Bus", "B1"; v_nom=380.0)` |
| `n.add("Generator", ..., committable=True)` | `add!(net, "Generator", ...; committable=true)` |
| `n.pf()` | `pf(net)` |
| `n.optimize()` | `optimize(net)` |
| `n.generators_t.p["G1"]` | `result.P_gen["G1"]` |
| `n.buses_t.marginal_price["B1"]` | `result.lmp["B1"]` |

See [CHEATSHEET.md](CHEATSHEET.md) for the full API comparison.

---

## Not Implemented

- Capacity expansion planning (`p_nom_extendable`)
- Security-constrained OPF (N-1 contingency)
- Sector coupling (heat, hydrogen, transport)
- PyPSA netCDF/HDF5 data format reader

---

## License

Academic use. Innopolis University, 2026.
