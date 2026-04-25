# Julia-Based Power Flow and Optimal Power Flow Solvers

**Master's Thesis — Innopolis University, 2026**

> *Porting Key PyPSA Modules to Julia for Improved Power System Simulation Performance*

---

## Abstract

This project implements and benchmarks five power system analysis modules in Julia, each paired with an equivalent PyPSA (Python) reference implementation. The work evaluates Julia's suitability as a high-performance replacement for Python-based power system toolchains, targeting the formulations used in PyPSA: DC Power Flow, AC Power Flow, Linearized AC Power Flow, Linear Optimal Power Flow (LOPF), and Multi-Period LOPF with storage and renewable generation.

---

## Implemented Methods

| Method | Julia | Python (PyPSA) | Network |
|---|---|---|---|
| DC Power Flow (DCPF) | `julia/solvers/dc_power_flow.jl` | `python/dc_power_flow.py` | 3-bus |
| AC Power Flow (ACPF) | `julia/solvers/ac_power_flow.jl` | `python/ac_power_flow.py` | 3-bus |
| Linearized AC Power Flow (LACPF) | `julia/solvers/linear_ac_pf.jl` | `python/linear_ac_pf.py` | 3-bus |
| Linear OPF (LOPF) | `julia/solvers/lopf.jl` | `python/lopf.py` | 3-bus |
| Multi-Period LOPF | `julia/solvers/lopf_multiperiod.jl` | `python/lopf_multiperiod.py` | 3-bus |

Validation against MATPOWER reference data is performed on the **IEEE 14-bus** test case (`julia/validation/ieee14_validation.jl`).

---

## Methods

### DC Power Flow (DCPF)

The linearised (lossless) DC approximation neglects resistances and reactive power, reducing the power flow equations to a linear system:

```
B · θ = P_inj
```

where **B** is the susceptance matrix and **P_inj** is the vector of net active power injections. Solved via sparse LU factorisation.

### AC Power Flow (ACPF)

The full nonlinear AC power flow is formulated as a system of 2(n−1) equations in voltage magnitudes |V| and angles θ, solved by Newton–Raphson iteration. The Julia implementation uses [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) with the Ipopt interior-point solver.

### Linearized AC Power Flow (LACPF)

LACPF linearises the full AC power flow equations around the flat-start operating point (|V| = 1 p.u., θ = 0 rad) using a first-order Taylor expansion. Unlike DCPF, it simultaneously solves for both active and reactive power, recovering approximate voltage magnitudes:

```
[ B'  -G' ] [ Δθ  ]   [ P_inj / S_base ]
[-G'  -B' ] [ Δ|V|] = [ Q_inj / S_base ]
```

where G' + jB' is the reduced nodal admittance matrix (slack bus eliminated).

### Linear Optimal Power Flow (LOPF)

Single-period economic dispatch subject to DC network constraints:

```
min  Σᵢ cᵢ · Pᵢ
s.t. B · θ = P_inj           (nodal power balance)
     |b_km · (θ_k − θ_m)| ≤ P_max   (line thermal limits)
     0 ≤ Pᵢ ≤ Pᵢ_max               (generator capacity)
```

Julia: [JuMP.jl](https://jump.dev/) + [HiGHS](https://highs.dev/) LP solver.
Python: `network.optimize()` via PyPSA (HiGHS backend).

### Multi-Period LOPF with Storage and Wind

Extends LOPF to a 24-hour horizon with time-coupled storage state-of-charge constraints, curtailable wind generation, and an hourly load profile:

```
min  Σ_t Σᵢ cᵢ · Pᵢ(t)
s.t.  Power balance per bus per period
      Generator capacity bounds
      Line thermal limits
      Storage SoC dynamics:  E(t+1) = E(t) + η_ch·P_ch(t) − P_dis(t)/η_dis
      SoC bounds:  0 ≤ E(t) ≤ E_nom
      Wind output:  P_wind(t) ≤ CF(t) · P_nom
```

---

## Validation

### 3-Bus Test Network

```
G1 (cheap, 20 €/MWh)       G2 (expensive, 50 €/MWh)
     [Bus 1] ──── Line 1-2 ──── [Bus 2]
         \                          /
       Line 1-3              Line 2-3
               \              /
              [Bus 3] (load 300 MW)
```

- Lines: x = 0.1 pu, r = 0.01 pu
- Loads: Bus 2 = 200 MW, Bus 3 = 300 MW
- Generators: G1 at Bus 1 (400 MW max), G2 at Bus 2 (300 MW max)

**DC Power Flow — Julia vs PyPSA:**
| Bus | Julia θ (rad) | PyPSA θ (rad) | |Δ| |
|---|---|---|---|
| Bus 0 | 0.000000 | 0.000000 | 0 |
| Bus 1 | −0.000185 | −0.000185 | < 1×10⁻⁶ |
| Bus 2 | −0.000162 | −0.000162 | < 1×10⁻⁶ |

**AC Power Flow — Julia vs PyPSA:**
| Metric | Julia | PyPSA | Max |Δ| |
|---|---|---|---|
| V_mag (p.u.) | 1.0 / 0.999982 / 0.999984 | same | 4.8×10⁻⁷ |
| V_ang (rad) | 0.0 / −0.000185 / −0.000162 | same | 4.1×10⁻⁷ |
| P_line (MW) | 266.67 / 233.34 / −33.33 | same | 3.5×10⁻³ |

**LOPF — Scenario A (no congestion):**
| | Julia | PyPSA |
|---|---|---|
| G1 dispatch | 400.0 MW | 400.0 MW |
| G2 dispatch | 100.0 MW | 100.0 MW |
| Total cost | 13 000 €/h | 13 000 €/h |

**LOPF — Scenario B (line limit 200 MW, congestion):**
| | Julia | PyPSA |
|---|---|---|
| G1 dispatch | 300.0 MW | 300.0 MW |
| G2 dispatch | 200.0 MW | 200.0 MW |
| Total cost | 16 000 €/h (+23.1%) | 16 000 €/h (+23.1%) |

### IEEE 14-Bus Validation

The AC power flow and LOPF solvers are validated against the MATPOWER case14 reference solution (`data/case14.m`). Voltage magnitudes and angles are compared across all 14 buses; deviations are reported in `julia/validation/ieee14_validation.jl`.

---

## Project Structure

```
thesis/
├── julia/
│   ├── solvers/
│   │   ├── dc_power_flow.jl        # DC Power Flow (B·θ = P)
│   │   ├── ac_power_flow.jl        # AC Power Flow (PowerModels.jl + Ipopt)
│   │   ├── linear_ac_pf.jl         # Linearized AC PF (LACPF)
│   │   ├── lopf.jl                 # Single-period LOPF (JuMP + HiGHS)
│   │   └── lopf_multiperiod.jl     # 24-hour LOPF with storage and wind
│   ├── benchmarks/
│   │   ├── benchmark.jl            # DCPF + LOPF benchmark runner
│   │   ├── benchmark_ac.jl         # ACPF benchmark runner
│   │   ├── benchmark_multiperiod.jl# Multi-period LOPF benchmark runner
│   │   └── plot_benchmarks.jl      # Benchmark result visualisation
│   ├── validation/
│   │   └── ieee14_validation.jl    # IEEE 14-bus validation vs MATPOWER
│   └── visualization/
│       ├── visualize_dcpf.jl
│       └── visualize_results.jl
├── python/
│   ├── dc_power_flow.py            # DC PF (PyPSA lpf)
│   ├── ac_power_flow.py            # AC PF (PyPSA pf)
│   ├── linear_ac_pf.py             # Linearized AC PF (PyPSA)
│   ├── lopf.py                     # Single-period LOPF (PyPSA optimize)
│   ├── lopf_multiperiod.py         # Multi-period LOPF (PyPSA)
│   ├── ieee14_validation.py        # IEEE 14-bus reference comparison
│   ├── benchmark.py                # DCPF + LOPF benchmarks
│   ├── benchmark_ac.py             # ACPF benchmarks
│   ├── benchmark_multiperiod.py    # Multi-period LOPF benchmarks
│   └── compare_benchmarks.py       # Cross-language performance comparison
├── data/
│   └── case14.m                    # IEEE 14-bus MATPOWER case
├── results/                        # Generated benchmark outputs (CSV + PNG)
└── tests/
    └── test_julia.jl
```

---

## Requirements

**Julia** (≥ 1.10):
```julia
using Pkg
Pkg.add(["PowerModels", "Ipopt", "JuMP", "HiGHS", "Plots", "LinearAlgebra"])
```

**Python** (≥ 3.9):
```bash
pip install pypsa highspy
```

---

## Running

```bash
# DC Power Flow
julia julia/solvers/dc_power_flow.jl
python python/dc_power_flow.py

# AC Power Flow
julia julia/solvers/ac_power_flow.jl
python python/ac_power_flow.py

# Linearized AC Power Flow
julia julia/solvers/linear_ac_pf.jl
python python/linear_ac_pf.py

# Single-period LOPF
julia julia/solvers/lopf.jl
python python/lopf.py

# Multi-period LOPF (24-hour, with storage and wind)
julia julia/solvers/lopf_multiperiod.jl
python python/lopf_multiperiod.py

# IEEE 14-bus validation
julia julia/validation/ieee14_validation.jl
python python/ieee14_validation.py

# Benchmarks
julia julia/benchmarks/benchmark.jl
julia julia/benchmarks/benchmark_ac.jl
julia julia/benchmarks/benchmark_multiperiod.jl
python python/compare_benchmarks.py
```

---

## Performance Notes

Julia's performance advantages over Python/PyPSA in this context stem from three complementary factors. First, JIT compilation via LLVM produces native machine code on the first call, eliminating Python's interpreter overhead for subsequent runs. Second, type-stable Julia code enables aggressive compiler optimisations that are unavailable in dynamically typed Python. Third, multiple dispatch allows the JuMP modelling layer to specialise constraint and objective code paths at compile time, reducing solver interface overhead. These properties are particularly significant for iterative solvers (Newton–Raphson in ACPF) and large-scale LP problems (multi-period LOPF), where inner-loop costs dominate total runtime.
