# Defense Q&A + Study Roadmap — PowerFlowJulia

Master's thesis: *AI-Assisted Simulation and Optimization of Power Market Projects*
Innopolis University, 2026.

This is a defense-preparation document: anticipated committee questions with honest
answers, plus a study roadmap (what to master deeply vs. what to understand at a
high level). Every number here must come from a re-runnable script — see
`results/` CSVs. Where a number is currently inconsistent across sources, it is
flagged in Part 0.

> Russian version: [DEFENSE_QA_RU.md](DEFENSE_QA_RU.md).

---

## Part 0 — Open issues to FIX before defense

1. **AC PF numbers — RESOLVED 2026-06-19.** The honest CSV figures are:

   | Buses | Julia min (ms) | PyPSA min (ms) | Speedup (min) |
   |---|---|---|---|
   | 3   | 7.22  | 177.5 | 24.6× |
   | 10  | 29.9  | 188.3 | 6.3×  |
   | 50  | 151.1 | 225.2 | 1.49× |
   | 100 | 549.0 | 264.2 | **0.48× (Julia ~2× slower)** |

   Status: the **LaTeX thesis is already correct** (`tab:ac_benchmark` shows
   24.6×→0.5×, analysis paragraph and figure caption state Julia is slower at 100
   buses). The stale 50.6×→7.2× numbers were only in the abandoned `thesis_docx/`
   markdown (gitignored) and in README — README is now fixed. No thesis action
   needed; just don't reintroduce the docx numbers. Conclusion to defend: Julia AC
   PF wins only on small networks, parity ~50 buses, slower than PyPSA's
   Newton–Raphson at 100 — algorithmic (Ipopt vs NR), not linguistic.

2. **`real_data_forecast.csv` artifact — DONE 2026-06-19.** Re-ran both configs;
   committed `results/ml/real_data_forecast_temp.csv` + `_notemp.csv`. The old
   3.67%/5.83% LSTM numbers did NOT reproduce — authoritative is now
   **2.77 / 4.03 / 5.42 / 7.73** (seasonal / LSTM+temp / LSTM-calendar / persistence).
   Thesis `tab:forecast_metrics`, text, README and `fig_ml_accuracy.png` all updated;
   the figure now reads the CSVs directly so it cannot drift. Conclusion unchanged.

3. **Audit thesis table/figure → source CSV — DONE 2026-06-19.** Every benchmark
   table (DC, LOPF, AC, realcases, MP, UC-by-T, UC-scale) and the stochastic table
   match their CSVs exactly (checked by min-time, speedups recompute). Two text fixes:
   (a) "seven load scenarios" → "thirty" (code/CSV use 30); (b) README Key-Results
   table was stale (500-bus rows absent from current CSVs, UC/MP mismatched) — now
   synced to the CSVs, headline reframed (the old "750×" had no source).

---

## Part 1 — Anticipated questions and honest answers

### A. Motivation, scope, contribution

**Q1. What is the contribution? Is this a new algorithm?**
No. It is an implementation + performance study + an AI add-on. The numerical
solving is delegated to the same class of tools PyPSA uses (HiGHS for LP/MILP,
Ipopt for AC). The contribution is: (a) a PyPSA-compatible power-flow/OPF library
in Julia with a leaner model-build layer, (b) a fair, reproducible benchmark of
where that layer wins and loses, and (c) an AI layer (LSTM + conformal +
stochastic LOPF) that adds calibrated uncertainty quantification.

**Q2. Why Julia and not just optimizing PyPSA?**
PyPSA's per-solve cost is dominated by Python/pandas/linopy model assembly that
scales with problem size. Julia assembles the same optimization model in compiled,
type-stable code via JuMP, eliminating interpreter overhead — without changing the
solver. So the comparison isolates the language/model-build layer.

**Q3. Why "AI-assisted" — isn't the AI bolted on?**
The thesis title requires an AI component, and we are honest that the AI is a
distinct layer, not woven into the solver. Its justified role is uncertainty
quantification: the LSTM's conformal intervals feed a scenario-based stochastic
LOPF that produces a risk profile (E[cost], CVaR₉₀) no deterministic forecast can.

### B. Power flow theory

**Q4. Derive the DC power flow. Why is it "DC"?**
Assumptions: flat voltage (|V|≈1 p.u.), small angle differences (sin θ ≈ θ),
negligible line resistance (r ≪ x), lossless. The active power balance linearizes
to **B·θ = P**, where B is the network susceptance (Laplacian-like) matrix, θ the
bus voltage angles, P the net injections. The slack bus fixes the angle reference.
"DC" is by analogy — the linear form resembles a resistive DC circuit, not because
it models DC current. Line flow: p_km = b_km (θ_k − θ_m).

**Q5. How is AC PF solved here vs. in PyPSA?**
We cast AC PF as a nonlinear feasibility/optimization problem and solve it with
PowerModels.jl + **Ipopt** (interior-point). PyPSA applies **Newton–Raphson**
directly to the power-flow mismatch equations. Different algorithms — which is
exactly why the AC benchmark measures algorithm choice, not just language overhead
(see Part 0 #1: NR wins at scale).

**Q6. What is the Linearized AC PF (LACPF) and when is it used?**
A first-order Taylor expansion around the operating point that recovers |V|
deviations and reactive flows Q, which plain DC PF discards. Solved as a single
coupled linear system (2(n−1)×2(n−1)) via LAPACK backslash. `pf(net, method=:lac)`;
auto-selected when any line has r > 0.

**Q7. The 3-bus DC angles differ from PyPSA by a factor of 1444 — bug?**
No, a documented per-unit convention difference. PyPSA treats `Line.x` as ohmic and
divides by V²/S_base = 380²/100 = 1444; we follow MATPOWER and treat x as already
per-unit (b = S_base/x). The angle scale cancels in the flow expression, so flows,
dispatch and prices are identical. On IEEE-14 (physically consistent per-unit data)
angles agree to < 0.001°.

### C. Optimization, markets, storage

**Q8. Write the single-period LOPF (DC-OPF).**
Minimize Σ c_g·p_g subject to: nodal balance B·θ = P_inj(p_g − d) (the balance
constraint carries the LMP dual), generator limits p_min ≤ p_g ≤ p_max, line limits
|p_km| ≤ s_nom. LP solved by HiGHS.

**Q9. How do you get LMPs? Why the minus sign?**
LMP at bus k = **−dual(balance_con[k])**. Load enters the balance RHS as −d, so
∂(objective)/∂d = −dual; the minus follows JuMP's minimization sign convention. The
LMP is the marginal cost of serving one more MW at that bus; differences across
buses signal congestion (congestion rent).

**Q10. How is multi-period LOPF different? How is storage modeled?**
Couples T snapshots through a StorageUnit state-of-charge recursion:
SoC_t = (1−loss)·SoC_{t−1} + η_ch·p_ch,t − p_dis,t/η_dis, with cyclic boundary
SoC_0 = SoC_T. This enables time-arbitrage: charge when cheap/abundant (wind),
discharge at peak. Wind enters as a time-varying p_max_pu profile.

**Q11. What is Unit Commitment and why is it harder?**
UC adds binary on/off variables u_g,t ∈ {0,1} with min-up/min-down-time and
startup/shutdown costs → a MILP, solved by HiGHS branch-and-bound. Harder because
the feasible set is non-convex (combinatorial). LMPs are recovered by fixing the
binaries and re-solving the LP relaxation to extract duals.

**Q12. Your UC cost is lower than PyPSA's (34,100 vs 36,800) — is that a bug?**
No — verified it is a PyPSA min-up-time boundary artifact at the horizon edge; our
solution is the true optimum. (Documented in the validation notes.)

### D. Performance & benchmarking — the honesty landmines

**Q13. "You only beat PyPSA's fixed Python overhead, not real compute."**
Partly true and stated openly. Small-network ratios (10⁴–10⁵×) are PyPSA's ~130 ms
fixed overhead divided by a near-instant Julia solve — not compute superiority. The
honest, compute-bound figures are the large real networks: PEGASE-2869 LOPF 58×,
DC PF 4.7×. Same HiGHS / same sparse solve on both sides — the gap is model-build.

**Q14. Then why does LOPF speedup stay high (58×) while DC PF collapses (4.7×)?**
DC PF is a single sparse linear solve — once the network is large, both sides spend
their time in the same O(n^~2.4) solve, so the ratio → language overhead only. LOPF
builds a whole LP (O(n) vars, O(n+|E|) constraints); PyPSA's linopy assembles it via
Python object loops, Julia's JuMP in compiled code — that assembly gap grows with
problem size and dominates PyPSA's runtime (its solve time is the same HiGHS).

**Q15. Is AC PF faster?** No — see Part 0 #1. Faster only at small N, slower than
NR at 100 buses. We report it as a limitation; a native Julia NR (no Ipopt wrapper)
is future work.

**Q16. Are the benchmarks fair / reproducible?**
Same synthetic topology generator (seed=42) on both sides; pinned PyPSA 1.2.2 /
Julia 1.12.5; HiGHS and BLAS pinned to a single thread; network built once, only
the per-solve (model-build + solve) timed; JIT warm-up excluded; 3 seeds, report
**minimum** (noise-robust — OS jitter and GC only add time). Also validated on
standard IEEE/PEGASE cases, not only synthetic.

**Q17. Why report minimum, not mean?**
Minimum is the standard estimator of intrinsic compute time: garbage collection and
OS scheduling can only add latency. We also report mean ± std; on large Julia DC PF
the std is large (GC pauses), which is why min is the cleaner figure.

### E. Julia & implementation

**Q18. What specifically makes Julia's model-build faster?**
JIT compilation to native code, type stability (immutable structs, concrete types →
no boxing), and multiple dispatch (the compiler specializes assembly code per type).
JuMP builds the constraint matrix directly in compiled code vs. Python object loops.

**Q19. Isn't JIT warm-up a hidden cost you're hiding?**
It is a one-time 2–5 s time-to-first-execution, excluded by a warm-up call — fair
because the target use case (energy-system studies) runs thousands of solves per
run, amortizing it to negligible. We state this explicitly.

**Q20. How is the API PyPSA-compatible?**
`add!(net, "Generator", "G1"; bus=..., p_nom=...)` mirrors `n.add(...)`; attribute
names match (`p_nom`, `marginal_cost`, `x`, `s_nom`, `efficiency_charge`, ...);
`pf()`↔`n.pf()`/`n.lpf()`, `optimize()`↔`n.optimize()`. Migration table in CLAUDE.md
and CHEATSHEET.md.

### F. AI component

**Q21. Describe the LSTM forecaster.**
Flux.jl LSTM (hidden=64), weekly input window (168 h) plus calendar features
(sin/cos hour, sin/cos day-of-week, weekend) and optional temperature; predicts the
24 h day-ahead residual to a seasonal-naive target; trained batched over full
sequences on real German load (OPSD/ENTSO-E), with a held-out conformal calibration
set.

**Q22. Honest result: does the LSTM beat the baselines?**
It beats persistence (MAPE 3.67% vs 7.73%, +61.5% RMSE skill) but **not**
seasonal-naive (2.77% MAPE). Aggregate national load has a near-deterministic weekly
cycle that seasonal-naive copies for free; a learned model only adds value through an
exogenous signal (temperature closes the gap from −84% to −23%). We do not spin this.

**Q23. Then what is the AI worth?**
Calibrated uncertainty, not point accuracy. Conformal prediction gives a
distribution-free ≥90% coverage interval (empirical coverage 100% on the 60-day test
— conservative because the non-conformity score is the per-day max error). That
interval drives the stochastic LOPF.

**Q24. Explain conformal prediction's guarantee.**
Given exchangeable data, for miscoverage α the interval [ŷ − q, ŷ + q], where q is
the ⌈(1−α)(n+1)⌉/n empirical quantile of calibration non-conformity scores, contains
the true value with probability ≥ 1−α — distribution-free, no model assumptions. We
use α = 0.10.

**Q25. What is the stochastic LOPF and CVaR₉₀?**
Sample Average Approximation: draw S demand scenarios from the conformal band, solve
an independent LOPF per scenario, aggregate. Outputs E[cost], spread σ, and
**CVaR₉₀** = mean cost over the worst 10% of scenarios (the tail-risk reserve). On
the 3-bus T=24 demo: E[cost] 236,080 €, CVaR₉₀ 261,283 €, σ 13,895 € — a
+10.7% reserve a risk-averse operator must budget. No deterministic forecast
(including the more accurate seasonal-naive) yields this.

**Q26. Why not a better model (XGBoost, weather model) to beat seasonal-naive?**
Out of scope and deliberately cut (XGBoost was tried then removed to keep the AI
layer minimal and honest). Beating seasonal-naive needs a real weather-forecast
feature or a residual model — named as future work. The thesis claim is the
uncertainty layer, not SOTA accuracy.

### G. Validation & reproducibility

**Q27. How is correctness validated?**
3-bus vs PyPSA (DC exact to <10⁻¹⁰ rad, LOPF exact 0 €, AC < 5×10⁻⁷ p.u.);
IEEE-14 vs MATPOWER (AC |V| error 1.33×10⁻³ p.u.); IEEE-14 DC vs PyPSA identical.
Plus the standard IEEE/PEGASE benchmark cases.

**Q28. Can the committee reproduce your numbers?**
Yes: pinned versions (`python/requirements.txt`, `Project.toml`), fixed seeds,
scripts in `benchmarks/` and `examples/` regenerate every CSV in `results/`.

### H. Limitations & future work

**Q29. What is NOT implemented?**
Capacity expansion planning, SCOPF (N-1), ShuntImpedance, DCLine, multi-node Links,
sector coupling. All low-priority and listed; the AC PF speedup ceiling and the
LSTM accuracy ceiling are the substantive limitations.

**Q30. If you had more time?**
Native Julia Newton–Raphson AC PF (recover the AC speedup), weather-informed
residual forecaster (beat seasonal-naive), capacity expansion, SCOPF.

---

## Part 2 — Study roadmap

### Master deeply (must defend cold, derive on the whiteboard)

1. **DC power flow derivation** — assumptions → B·θ = P → flow p_km = b(θ_k−θ_m).
   The single most likely derivation question.
2. **LP duality → LMP** — why LMP = −dual(balance), congestion rent, KKT intuition.
   This is the market core of the thesis.
3. **DC-OPF / LOPF formulation** — objective, all constraints, what each dual means.
4. **Multi-period storage SoC recursion** + cyclic constraint (time-arbitrage).
5. **Unit Commitment MILP** — binaries, min up/down, why non-convex, how LMP is
   recovered via LP relaxation.
6. **Why the speedup exists** — model-build vs solver; why DC ratio falls and LOPF
   ratio holds; why AC is slower (the honesty story). Know the real numbers.
7. **Conformal prediction guarantee** — the one ML proof that is your contribution;
   be able to state the quantile construction and the exchangeability assumption.
8. **CVaR definition** and why it (not the mean) is the risk metric.
9. **Benchmark methodology** — seeds, min vs mean, warm-up, thread pinning, fairness
   (same solver both sides).

### Understand the work (explain at a high level; delegated to libraries)

1. **Ipopt interior-point internals** — barrier method idea; you delegate to it.
   Know NR vs IPM trade-off (Q5), not the line-search math.
2. **HiGHS simplex / branch-and-bound internals** — know what B&B does, not its
   pivoting rules.
3. **LSTM gate equations** — conceptual (memory cell, gates handle long sequences);
   Flux + Zygote autodiff compute gradients; you need not derive backprop-through-time.
4. **JuMP/Zygote/Flux implementation details** — know what they provide, not their
   internals.
5. **PowerModels.jl data format** — know it needs full per-unit branch data; the
   field-completeness pain is an implementation note, not theory.
6. **Newton–Raphson for AC PF** — know the mismatch-equation + Jacobian iteration
   conceptually (relevant to the AC limitation and future work).

### Background to refresh (assumed, low risk but askable)

- Per-unit system and why it matters (the 1444 factor, Q7).
- Power system components (bus/line/transformer/generator/load) and their parameters.
- Basic LP/MILP definitions; convex vs non-convex.
- Time-series error metrics (MAPE, RMSE) and naive baselines.

### Suggested study order

1. Power-flow + DCOPF + LMP/duality (deep block 1–3) — the spine of Q4–Q9.
2. Storage + UC (deep 4–5) — Q10–Q12.
3. Performance story + benchmark methodology (deep 6, 9) — Q13–Q19, the landmines.
4. AI: conformal + CVaR + stochastic LOPF (deep 7–8) — Q21–Q26.
5. Skim the "understand" list so you can gracefully say "delegated to <library>,
   here's the trade-off" instead of bluffing internals.
