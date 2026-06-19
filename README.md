# PowerFlowJulia

**AI-Assisted Simulation and Optimization of Power Market Projects**
Master's Thesis — Innopolis University, 2026
Author: Dadakhon Turgunboev · Supervisor: Leonard Johard

---

## Overview

PowerFlowJulia is a Julia implementation of core power system simulation and optimization algorithms, reimplementing the key functionality of [PyPSA](https://pypsa.org/) with a PyPSA-compatible API. The library delivers large speedups over PyPSA's model-build layer — e.g. **118–283× on synthetic LOPF** and **58× on the 2869-bus PEGASE European grid** — using the *same* HiGHS/Ipopt solvers on both sides (AC power flow is the honest exception: slower than PyPSA's Newton–Raphson beyond ~50 buses). It also integrates an **LSTM-based load forecaster** with conformal prediction intervals for uncertainty-aware stochastic dispatch.

---

## Key Results

| Method | Network size | Julia (ms) | PyPSA (ms) | Speedup |
|---|---|---|---|---|
| DC Power Flow | 300 buses (synthetic) | 2.14 | 360.8 | **169×** |
| LOPF | 300 buses (synthetic) | 17.1 | 2,022 | **118×** |
| AC Power Flow | 3 buses (synthetic) | 7.2 | 177.5 | **24.6×** |
| AC Power Flow | 100 buses (synthetic) | 549.0 | 264.2 | **0.5× (slower)** |
| Unit Commitment | 14 buses, T=24 | 42.6 | 730.0 | **17×** |
| Multi-Period LOPF | T=24 | 8.7 | 1,141 | **131×** |
| LOPF | PEGASE 2869 (real EU grid) | 491 | 28,427 | **58×** |
| DC Power Flow | PEGASE 2869 (real EU grid) | 423 | 2,005 | **4.7×** |

Speedups come from Julia's leaner model-build layer (no Python/pandas/linopy overhead), **not** a faster solver: LP/MILP use the same HiGHS on both sides, AC PF uses Ipopt. Small-network ratios are dominated by PyPSA's fixed Python overhead; the compute-bound figures (large real networks) are the honest ones. AC PF is the one method where the port does **not** win: Julia (PowerModels.jl + Ipopt) is faster only on small networks (24.6× at 3 buses), reaches parity near 50 buses (1.5×), and is about **2× slower** than PyPSA's specialised Newton–Raphson at 100 buses (0.5×) — an algorithmic, not linguistic, gap. All figures are the minimum over multiple seeds on pinned PyPSA 1.2.2 / Julia 1.12.5.

**AI component:** LSTM day-ahead load forecaster trained on **real** German load (OPSD/ENTSO-E). On a strict out-of-time test it beats the persistence baseline (MAPE **4.03%** vs 7.73%) but **not** the seasonal-naive baseline (2.77%) — aggregate national load has a near-deterministic weekly cycle. The contribution is therefore **not** point accuracy but the calibrated **conformal prediction interval** (≥90% coverage, distribution-free), which feeds a scenario-based **stochastic LOPF** yielding E[cost] and CVaR₉₀ for risk-aware dispatch — a risk profile no deterministic forecast provides.

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
| **LSTM Forecaster** | Day-ahead load forecast, weekly window=168h + calendar/temperature features → horizon=24h | Flux.jl |
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
# 1. Load real German load data (OPSD/ENTSO-E); synthetic generator is a fallback
r = load_real_data(; normalize=:peak)          # (n_days × 24) pu matrix
data = r.data
# data = generate_synthetic_data(365; noise_std=0.05, seed=42)   # fallback only

# 2. Train LSTM forecaster with conformal calibration
fc = train_forecaster(data; hidden=64, epochs=80, window=168, verbose=true)

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

**Python 3.12** (reference benchmarks only; versions pinned for reproducibility)

```bash
pip install -r python/requirements.txt   # PyPSA 1.2.2, linopy 0.8.0, highspy 1.14.0, ...
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
