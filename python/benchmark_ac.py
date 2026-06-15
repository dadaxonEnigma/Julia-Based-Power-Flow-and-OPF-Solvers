"""
Benchmark: Python/PyPSA — AC Power Flow (Newton-Raphson).
Run alongside julia/benchmarks/benchmark_ac.jl for comparison.

Statistics: SEEDS random topologies × runs_per_seed timed solves; mean ± std + min.
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


def generate_network_ac(n_buses, seed=42):
    rng = np.random.default_rng(seed)

    lines = []
    for i in range(1, n_buses):
        x = 0.05 + rng.random() * 0.45
        lines.append((i, i + 1, 0.01, x))
    for _ in range(max(1, n_buses // 3)):
        u = int(rng.integers(1, n_buses))
        v = int(rng.integers(u + 1, n_buses + 1))
        x = 0.05 + rng.random() * 0.45
        lines.append((u, v, 0.01, x))

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


def build_pypsa_network_ac(n_buses, lines, generators, loads):
    net = pypsa.Network()
    net.set_snapshots([0])
    for bus in range(1, n_buses + 1):
        net.add("Bus", f"Bus{bus}", v_nom=380.0)
    for idx, (f, t, r, x) in enumerate(lines):
        net.add("Line", f"L{idx}", bus0=f"Bus{f}", bus1=f"Bus{t}", x=x, r=r, s_nom=1e6)
    for bus, p_max in generators.items():
        ctrl = "Slack" if bus == 1 else "PV"
        net.add("Generator", f"G{bus}", bus=f"Bus{bus}",
                p_nom=p_max, p_set=p_max * 0.5, control=ctrl)
    for bus, p in loads.items():
        net.add("Load", f"Load{bus}", bus=f"Bus{bus}", p_set=p)
    return net


def bench_stats(run_factory, runs_per_seed):
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
print("BENCHMARK: Python/PyPSA — AC Power Flow (Newton-Raphson)")
print(f"Seeds: {SEEDS}")
print("=" * 74)

AC_SIZES = [3, 10, 50, 100]

print(f"\n{'Buses':<8} {'Mean (ms)':>12} {'Std (ms)':>12} {'Min (ms)':>12} {'N':>6}")
print("-" * 74)

ac_results = {}
for n in AC_SIZES:
    runs = 10 if n <= 10 else (8 if n <= 50 else 5)

    def factory(s, n=n):
        net = build_pypsa_network_ac(*generate_network_ac(n, seed=s))
        return lambda: net.pf()

    mean, sd, mn, ns = bench_stats(factory, runs)
    ac_results[n] = (mean, sd, mn, ns)
    print(f"{n:<8} {mean:>12.3f} {sd:>12.3f} {mn:>12.3f} {ns:>6}")

with open("results/benchmarks/python_ac_benchmark.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["module", "n_buses", "time_ms", "std_ms", "min_ms", "n_samples"])
    for n in AC_SIZES:
        writer.writerow(["AC_PF", n, *ac_results[n]])

print("\n[OK] results/benchmarks/python_ac_benchmark.csv")
