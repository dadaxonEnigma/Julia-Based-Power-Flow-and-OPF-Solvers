# Conclusions and Future Work

This chapter draws together the findings of the thesis. The first section summarises the key quantitative results established in Chapter 3. The second section states the principal conclusions regarding Julia's viability as a high-performance alternative to Python for power system computation, addressing in turn each research question posed in the Introduction. The final section identifies a structured programme of future work that extends the implemented modules toward a broader Julia-based power system analysis framework.

## Summary of Results

This thesis implemented and evaluated eight computational contributions: three core power system solvers ported from PyPSA to Julia (DC Power Flow, AC Power Flow, and single-period LOPF), an extension to multi-period LOPF with storage and wind, Unit Commitment (MILP), an LSTM load forecaster with conformal prediction intervals, a scenario-based stochastic LOPF, and a systematic validation campaign on standard IEEE test cases. The key findings are summarised below.

**Numerical correctness on 3-bus network.**
For the LP-based methods --- single-period LOPF, multi-period LOPF, locational marginal prices, and emission-constrained dispatch --- the Julia and PyPSA results are numerically identical (zero difference in generator dispatch, total cost, and nodal prices). For the power flow modules, active and reactive power flows agree to better than $10^{-3}$, while two quantities differ by a constant scale factor: the DC power flow voltage angles and the transformer branch flow. Both differences trace to a documented per-unit convention difference --- PyPSA normalises branch susceptance by $V_\text{nom}^2$ and the transformer reactance by its own rating $s_\text{nom}$, whereas the present implementation follows the MATPOWER convention of expressing susceptance in MW/rad on the system base. Crucially, these convention differences do not affect any economically meaningful quantity (dispatch, cost, or price), confirming the correctness of the porting.

**IEEE 14-bus standard test case validation.**
The implementations were additionally validated on the IEEE 14-bus benchmark network [@zimmerman2011matpower]. The AC power flow solution matches the MATPOWER reference to within $1.33 \times 10^{-3}$ p.u. in voltage magnitude and $0.017^{\circ}$ in angle. The DC power flow and LOPF solutions from Julia and PyPSA are identical to numerical precision. This confirms that the implementations generalise correctly beyond the construction network.

**DC Power Flow performance.**
Julia achieves speedups from $220{,}000\times$ (3 buses) to $27\times$ (2000 buses) over PyPSA. The decreasing speedup reflects PyPSA's fixed $\approx 130\,\text{ms}$ Python overhead dominating at small sizes, while Julia's efficient linear algebra narrows the gap at large sizes.

**AC Power Flow performance.**
Benchmarks over 3--100 buses yield speedups of $50\times$ (3 buses) to $7\times$ (100 buses). The smaller speedup compared to DC power flow reflects an algorithmic difference: Julia uses Ipopt (interior-point NLP) while PyPSA uses Newton--Raphson, a specialised iterative method for power flow. Despite this disadvantage, Julia remains faster at all tested sizes due to PyPSA's fixed Python overhead.

**Single-period LOPF performance.**
The Julia LOPF (JuMP + HiGHS) achieves speedups from $246\times$ (3 buses) to $753\times$ (500 buses), with the speedup *growing* with network size. Since both implementations use HiGHS as the LP solver, the entire speedup comes from JuMP's efficient model construction versus PyPSA/linopy's Python-based constraint loop. This trend is directly applicable to large-scale energy system studies such as PyPSA-Eur.

**Multi-period LOPF with storage and wind.**
The 24-hour multi-period LOPF with a StorageUnit and variable wind generation was successfully implemented and validated. Both Julia and PyPSA find the same optimal cost of $126{,}790.79\,\text{\euro}$. The LP correctly exploits time-arbitrage: the storage charges at night ($\approx 22\,\text{\euro/MWh}$ round-trip) and discharges at peak to displace expensive backup generation ($50\,\text{\euro/MWh}$). Performance benchmarks over horizons $T \in \{6, 12, 24, 48, 96\}$ hours show Julia achieving $87$--$371\times$ speedups. The Python solve time is nearly constant ($\sim 1.4\,\text{s}$) across all horizons, confirming that PyPSA's fixed Python overhead dominates for small networks. Julia's solve time scales linearly with $T$, from 3.5 ms at $T = 6$ to 17 ms at $T = 96$.

**Unit Commitment.**
The MILP Unit Commitment solver was implemented and benchmarked. Julia achieves $19$--$77\times$ speedups over PyPSA across horizons $T \in \{6, 12, 24, 48\}$ and network sizes of 3--40 buses. The speedup is lower than for continuous LOPF because the MILP branch-and-bound search --- executed by HiGHS in both implementations --- constitutes the dominant cost. LMPs are extracted via LP relaxation with binaries fixed, following the standard electricity market methodology [@oneill2005].

**AI component: LSTM forecasting and stochastic dispatch.**
An LSTM load forecaster (Flux.jl, hidden size 32, window 24 h, horizon 24 h) was trained on 365 days of synthetic load data. Early stopping at epoch 19 achieved a validation RMSE of 0.827 p.u. On a held-out test day, the model achieves MAE $= 0.103$ p.u. and MAPE $= 13.5\%$. Conformal calibration on a held-out set of 871 samples produces 90%-coverage prediction intervals with mean width 0.689 p.u.

Three dispatch strategies were compared on a 3-bus network over 24 hours:
(1) the naive default profile cost $320{,}328\,\text{\euro}$;
(2) the ML point-forecast deterministic LOPF cost $294{,}935\,\text{\euro}$ ($-7.9\%$);
(3) the stochastic SAA expected cost $309{,}711\,\text{\euro}$ with CVaR$_{90} = 349{,}879\,\text{\euro}$.
The 7.9% cost reduction relative to the naive baseline, and the CVaR$_{90}$--$\mathbb{E}[\text{cost}]$ gap of $40{,}168\,\text{\euro}$ (+13%), demonstrate the end-to-end forecasting-to-dispatch pipeline.

This AI component is presented as a methodological proof of concept rather than a production forecaster. Two limitations qualify its results and motivate the future work in Section [1.3](#sec:future){reference-type="ref" reference="sec:future"}. First, the LSTM is trained on *synthetic* load data; its MAPE of 13.5% and the wide conformal band (mean width 0.689 p.u.) reflect both the limited and noisy training signal and the deliberately compact model (hidden size 32). On real historical load series, short-term forecasters typically achieve 1--3% MAPE [@hong2016probabilistic], so the present figures should be read as an upper bound on error for this pipeline, not as the achievable accuracy. Second, the 7.9% improvement is measured against a static naive profile; a stronger baseline (persistence or seasonal ARIMA) would provide a more rigorous benchmark. The value of the component therefore lies in demonstrating that conformal-calibrated forecasts can be propagated into a stochastic LOPF entirely within Julia --- not in the specific accuracy numbers, which are bounded by the synthetic data.

## Conclusions

The experimental results demonstrate that porting PyPSA's core algorithms to Julia is both feasible and highly beneficial. The main conclusions are:

1.  **Julia is a viable replacement for Python in power system computation.**
    The mathematical formulations translate directly from PyPSA's Python code, and the Julia ecosystem (JuMP, HiGHS, PowerModels.jl) provides mature, well-maintained tools for every required component. Validation on both a custom 3-bus network and the standard IEEE 14-bus case confirms numerical equivalence.

2.  **The performance advantage is largest for LP-based optimisation.**
    The single-period LOPF speedup of $250$--$750\times$ grows with problem size, making Julia particularly valuable for iterative large-scale studies. A country-level PyPSA-Eur model that currently requires hours per run could, in principle, be reduced to minutes.

3.  **Multi-period optimisation with storage scales naturally.**
    Extending the single-period LOPF to 24 time periods with a StorageUnit required only incremental additions to the JuMP model --- SOC dynamics, charge/discharge variables, and a cyclicity constraint. The Julia implementation produces results identical to PyPSA, confirming that the time-coupling constraints are correctly formulated.

4.  **Algorithm choice matters more than language for AC power flow.**
    The AC power flow speedup ($7$--$50\times$) is lower than for DC/LOPF because Julia uses Ipopt (a general NLP solver) while PyPSA uses Newton--Raphson (a specialised iterative method). A native Newton--Raphson implementation in Julia would eliminate this algorithmic disadvantage and achieve much higher speedups, consistent with the pattern seen for DC power flow.

5.  **JIT warmup is not a practical barrier.**
    For realistic use cases --- time-series analysis, scenario sweeps, iterative optimisation --- functions are invoked thousands of times, making the one-time JIT compilation cost ($2$--$5\,\text{s}$) negligible. The benchmarks in this thesis measure steady-state (post-warmup) performance, which is the relevant metric for production use.

6.  **The existing Julia power system landscape has a gap.**
    As documented in the Related Work section, no actively maintained, self-contained Julia reimplementation of PyPSA exists. PSA.jl and EnergyModels.jl were abandoned years ago; PowerModels.jl and PowerSimulations.jl follow independent architectures. This thesis provides a validated starting point for filling that gap.

7.  **AI-assisted dispatch reduces cost and quantifies risk.**
    The LSTM forecaster reduces dispatch cost by 7.9% relative to a naive profile by providing a more accurate 24-hour load shape to the LOPF. The conformal prediction framework provides distribution-free 90% coverage intervals without any assumption on the residual distribution. The stochastic SAA provides operators with a CVaR$_{90}$ estimate that is inaccessible to any deterministic formulation, enabling explicit tail-risk management. These results demonstrate that data-driven forecasting is a practical, deployable complement to physics-based power flow optimisation.

## Future Work {#sec:future}

This thesis provides a foundation for a broader Julia-based power system analysis framework. Several concrete directions for future work are identified.

**Native Newton--Raphson AC power flow.**
The current AC power flow implementation delegates to Ipopt, which imposes an algorithmic penalty relative to PyPSA's Newton--Raphson. Implementing a direct Newton--Raphson solver in Julia --- building and factorising the Jacobian matrix natively --- would achieve speedups consistent with the DC power flow results ($27$--$220{,}000\times$) and constitute the most impactful single improvement.

**Parallel stochastic LOPF.**
The SAA implemented in Section [\[sec:ai_component\]](#sec:ai_component){reference-type="ref" reference="sec:ai_component"} solves $S$ independent LOPFs sequentially. Each scenario solve is fully independent, making the SAA embarrassingly parallel. Exploiting `Threads.@threads` or `Distributed.pmap` in Julia would reduce stochastic LOPF wall-clock time by a factor of $S$ on a multi-core machine, enabling real-time stochastic dispatch.

**Larger standard test cases.**
Validation was performed on the 3-bus construction network and the IEEE 14-bus benchmark. Extending to the IEEE 118-bus and 300-bus cases, and eventually to a subset of the PyPSA-Eur European transmission network [@horsch2018pypsa_eur], would demonstrate scalability to real-world problem sizes and provide a stronger empirical foundation for the performance claims.

**Real historical load data.**
The LSTM forecaster was trained on synthetic data. Substituting real smart-meter or TSO load datasets (e.g., ENTSO-E transparency platform) would provide a more rigorous evaluation of forecast accuracy and allow calibration of conformal quantiles on genuinely out-of-distribution future observations.

**Security-constrained OPF.**
Adding N-1 contingency constraints for transmission line or generator outages would produce a security-constrained LOPF (SCLOPF), which is a standard requirement for transmission system operator studies.

**Capacity expansion planning.**
Extending the multi-period LOPF with investment variables (`p_nom_extendable`) would enable greenfield or retrofit generation planning studies, equivalent to PyPSA's capacity expansion mode.

**PyPSA data format interoperability.**
Developing a Julia reader for PyPSA's native netCDF/HDF5 network format would enable seamless use of existing PyPSA-Eur model files in the Julia solvers, without requiring manual data conversion.

The results of this thesis demonstrate that a Julia-based power system analysis framework can deliver order-of-magnitude performance improvements over Python while maintaining numerical equivalence with PyPSA, and that integrating a data-driven forecasting layer enables uncertainty-aware, risk-quantified dispatch decisions that are inaccessible to deterministic formulations.
