"""
Benchmark: Python/PyPSA — DC Power Flow & LOPF.
Run alongside julia/benchmarks/benchmark.jl for comparison.

Statistics: each size is measured over SEEDS random topologies (identical seeds
and generator logic as the Julia side), runs_per_seed timed solves each after a
warm-up. CSV: time_ms = mean, plus std_ms, min_ms, n_samples.
"""
import os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"          # fair single-thread linear algebra (set before numpy)
import numpy as np
import pypsa
import time
import statistics
import csv
import warnings
import logging

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

SEEDS = [1, 2, 3]


def generate_network(n_buses, seed=42):
    rng = np.random.default_rng(seed)

    lines = []
    for i in range(1, n_buses):
        x = 0.05 + rng.random() * 0.45
        lines.append((i, i + 1, x))
    for _ in range(max(1, n_buses // 3)):
        u = int(rng.integers(1, n_buses))
        v = int(rng.integers(u + 1, n_buses + 1))
        x = 0.05 + rng.random() * 0.45
        lines.append((u, v, x))

    loads = {}
    total_load = 0.0
    for bus in range(2, n_buses + 1):
        if rng.random() > 0.3:
            p = 50.0 + rng.random() * 450.0
            loads[bus] = p
            total_load += p
    if not loads:
        loads[2] = 200.0
        total_load = 200.0

    generators = {1: total_load * 1.1}
    for bus in range(2, n_buses + 1, 4):
        generators[bus] = loads.get(bus, 0.0) * 0.5 + 50.0

    return n_buses, lines, generators, loads


def build_pypsa_network(n_buses, lines, generators, loads, line_capacity=1e6):
    n = pypsa.Network()
    for bus in range(1, n_buses + 1):
        n.add("Bus", f"Bus{bus}", v_nom=380.0)
    for idx, (f, t, x) in enumerate(lines):
        n.add("Line", f"L{idx}", bus0=f"Bus{f}", bus1=f"Bus{t}",
              x=x, r=0.01, s_nom=line_capacity)
    for bus, p_max in generators.items():
        ctrl = "Slack" if bus == 1 else "PQ"
        n.add("Generator", f"G{bus}", bus=f"Bus{bus}",
              p_nom=p_max, marginal_cost=20.0, control=ctrl)
    for bus, p in loads.items():
        n.add("Load", f"Load{bus}", bus=f"Bus{bus}", p_set=p)
    return n


def bench_stats(run_factory, runs_per_seed):
    """run_factory(seed) -> zero-arg callable. Pool samples over SEEDS × runs."""
    samples = []
    for s in SEEDS:
        run = run_factory(s)
        run()                                   # warm-up
        for _ in range(runs_per_seed):
            t0 = time.perf_counter()
            run()
            samples.append(time.perf_counter() - t0)
    mean = statistics.mean(samples) * 1000
    std = (statistics.stdev(samples) if len(samples) > 1 else 0.0) * 1000
    return mean, std, min(samples) * 1000, len(samples)


print("=" * 74)
print("BENCHMARK: Python/PyPSA — DC Power Flow & LOPF")
print(f"Seeds: {SEEDS}  |  mean +/- std over pooled samples")
print("=" * 74)

DC_SIZES   = [3, 10, 50, 100, 300]
LOPF_SIZES = [3, 10, 50, 100, 300]

print("\n[DC POWER FLOW]  net.lpf()")
print("-" * 74)
print(f"{'Buses':<8} {'Mean (ms)':>12} {'Std (ms)':>12} {'Min (ms)':>12} {'N':>6}")
print("-" * 74)

dc_results = {}
for n in DC_SIZES:
    runs = 10 if n <= 100 else (8 if n <= 500 else (5 if n <= 1000 else 3))

    def factory(s, n=n):
        net = build_pypsa_network(*generate_network(n, seed=s))
        return lambda: net.lpf()

    mean, sd, mn, ns = bench_stats(factory, runs)
    dc_results[n] = (mean, sd, mn, ns)
    print(f"{n:<8} {mean:>12.4f} {sd:>12.4f} {mn:>12.4f} {ns:>6}")

print("\n[LOPF]  net.optimize(solver='highs')")
print("-" * 74)
print(f"{'Buses':<8} {'Mean (ms)':>12} {'Std (ms)':>12} {'Min (ms)':>12} {'N':>6}")
print("-" * 74)

lopf_results = {}
for n in LOPF_SIZES:
    runs = 10 if n <= 50 else (8 if n <= 100 else 4)

    def factory(s, n=n):
        # Build the network ONCE (outside the timed call), matching the Julia
        # side, so the timing measures only optimisation-model assembly + solve,
        # not PyPSA network construction. optimize() rebuilds the linopy model
        # internally each call, just as Julia's solve() rebuilds the JuMP model.
        net = build_pypsa_network(*generate_network(n, seed=s))
        return lambda: net.optimize(solver_name="highs", solver_options={"threads": 1})

    mean, sd, mn, ns = bench_stats(factory, runs)
    lopf_results[n] = (mean, sd, mn, ns)
    print(f"{n:<8} {mean:>12.4f} {sd:>12.4f} {mn:>12.4f} {ns:>6}")

with open("results/python_benchmark.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["module", "n_buses", "time_ms", "std_ms", "min_ms", "n_samples"])
    for n in DC_SIZES:
        writer.writerow(["DC_PF", n, *dc_results[n]])
    for n in LOPF_SIZES:
        writer.writerow(["LOPF", n, *lopf_results[n]])

print("\n[OK] results/python_benchmark.csv")
