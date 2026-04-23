# PyPSA → Julia: Power Flow Porting

**Master's thesis — Innopolis University, 2026**

> "Porting Key PyPSA Modules to Julia for Improved Power System Simulation Performance"

---

## Overview

This project ports three core power flow modules from [PyPSA](https://pypsa.org/) (Python) to Julia, comparing correctness and performance of both implementations.

| Module | Julia | Python (PyPSA) | Status |
|---|---|---|---|
| DC Power Flow | `julia/dc_power_flow.jl` | `python/dc_power_flow.py` | ✅ Done |
| AC Power Flow | `julia/ac_power_flow_pm.jl` | `python/ac_power_flow.py` | ✅ Done |
| Linear OPF (LOPF) | `julia/lopf.jl` | `python/lopf.py` | ✅ Done |

---

## Methods

### DC Power Flow
Linearised (lossless) power flow. Solves **B·θ = P** via sparse LU factorisation.

### AC Power Flow
Full nonlinear power flow using Newton–Raphson.
Julia implementation uses [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) + Ipopt solver.

### Linear Optimal Power Flow (LOPF)
Economic dispatch with DC network constraints:

```
min  Σ cᵢ · Pᵢ
s.t. B·θ = Pᵢₙⱼ        (power balance)
     |bₖₘ·(θₖ−θₘ)| ≤ Pₘₐₓ   (line limits)
     0 ≤ Pᵢ ≤ Pᵢₘₐₓ          (generator limits)
```

Julia: [JuMP.jl](https://jump.dev/) + [HiGHS](https://highs.dev/) solver.
Python: `network.optimize()` in PyPSA (also uses HiGHS internally).

---

## Test Network

3-bus system used for validation across all three modules:

```
G1 (cheap, 20 €/MWh)       G2 (expensive, 50 €/MWh)
     [Bus 1] ──── Line 1-2 ──── [Bus 2]
         \                          /
       Line 1-3              Line 2-3
               \              /
              [Bus 3] (load 300 MW)
```

- **Lines:** x = 0.1 pu, r = 0.01 pu
- **Loads:** Bus 2 = 200 MW, Bus 3 = 300 MW
- **Generators:** G1 at Bus 1 (400 MW max), G2 at Bus 2 (300 MW max)

### Validation Results (Julia vs PyPSA)

**DC Power Flow:**
| Bus | Julia θ (rad) | PyPSA θ (rad) | Δ |
|---|---|---|---|
| Bus 0 | 0.000000 | 0.000000 | 0 |
| Bus 1 | −0.000185 | −0.000185 | < 1e-6 |
| Bus 2 | −0.000162 | −0.000162 | < 1e-6 |

**AC Power Flow:**
| Metric | Julia | PyPSA | Max Δ |
|---|---|---|---|
| V_mag | 1.0 / 0.999982 / 0.999984 | same | 4.8×10⁻⁷ |
| V_ang (rad) | 0.0 / −0.000185 / −0.000162 | same | 4.1×10⁻⁷ |
| Line P (MW) | 266.67 / 233.34 / −33.33 | same | 3.5×10⁻³ |

**LOPF — Scenario A (no line limits):**
| | Julia | PyPSA |
|---|---|---|
| G1 | 400.0 MW | 400.0 MW |
| G2 | 100.0 MW | 100.0 MW |
| Cost | 13 000 €/h | 13 000 €/h |

**LOPF — Scenario B (line limit 200 MW → congestion):**
| | Julia | PyPSA |
|---|---|---|
| G1 | 300.0 MW | 300.0 MW |
| G2 | 200.0 MW | 200.0 MW |
| Cost | 16 000 €/h (+23.1%) | 16 000 €/h (+23.1%) |

---

## Project Structure

```
thesis/
├── julia/
│   ├── solvers/                    # Core algorithm implementations
│   │   ├── dc_power_flow.jl        # DC Power Flow (B·θ = P)
│   │   ├── ac_power_flow.jl        # AC Power Flow (PowerModels.jl + Ipopt)
│   │   └── lopf.jl                 # LOPF (JuMP + HiGHS)
│   ├── benchmarks/                 # Performance measurement
│   │   ├── benchmark.jl            # Benchmark runner (DC PF + LOPF)
│   │   └── plot_benchmarks.jl      # Benchmark result visualisation
│   └── visualization/              # Network and result plots
│       ├── visualize_dcpf.jl       # DC PF network visualisation
│       └── visualize_results.jl    # General result visualisation
├── python/
│   ├── dc_power_flow.py            # DC Power Flow (PyPSA lpf)
│   ├── ac_power_flow.py            # AC Power Flow (PyPSA pf)
│   ├── lopf.py                     # LOPF (PyPSA optimize)
│   └── explore.py                  # PyPSA API exploration
├── latex/
│   ├── thesis.tex                  # Main LaTeX document
│   ├── introduction.tex            # Introduction + Related Work
│   ├── ch2_methodology.tex         # Chapter 2: Design and Methodology
│   ├── ch3_implementation.tex      # Chapter 3: Implementation and Experiments
│   ├── ch4_conclusions.tex         # Chapter 4: Conclusions and Future Work
│   └── literature.bib              # Bibliography
├── results/                        # Generated plots (PNG)
│   ├── benchmark_time.png
│   ├── benchmark_combined.png
│   └── results_validation.png
└── tests/
    └── test_julia.jl
```

---

## Requirements

**Julia** (≥ 1.9):
```julia
using Pkg
Pkg.add(["PowerModels", "Ipopt", "JuMP", "HiGHS", "Plots"])
```

**Python** (≥ 3.9):
```bash
pip install pypsa highspy
```

---

## Running

```bash
# DC Power Flow
julia julia/dc_power_flow.jl
python python/dc_power_flow.py

# AC Power Flow
julia julia/ac_power_flow_pm.jl
python python/ac_power_flow.py

# LOPF
julia julia/lopf.jl
python python/lopf.py
```
