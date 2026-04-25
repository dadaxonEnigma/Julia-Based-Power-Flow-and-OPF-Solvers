"""
Benchmark: Python/PyPSA — AC Power Flow (Newton-Raphson)
Run alongside julia/benchmarks/benchmark_ac.jl for comparison.
"""
import numpy as np
import pypsa
import time
import statistics
import csv
import warnings
import logging

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)


# ----------------------------------------------------------------
#  Network generator — same logic and seed as julia/benchmark_ac.jl
# ----------------------------------------------------------------
def generate_network_ac(n_buses, seed=42):
    rng = np.random.default_rng(seed)

    lines = []
    # Spanning tree
    for i in range(1, n_buses):
        x = 0.05 + rng.random() * 0.45
        lines.append((i, i + 1, 0.01, x))
    # Additional mesh edges
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
        net.add("Line", f"L{idx}",
                bus0=f"Bus{f}", bus1=f"Bus{t}",
                x=x, r=r, s_nom=1e6)

    for bus, p_max in generators.items():
        ctrl = "Slack" if bus == 1 else "PV"
        net.add("Generator", f"G{bus}", bus=f"Bus{bus}",
                p_nom=p_max, p_set=p_max * 0.5,
                control=ctrl)

    for bus, p in loads.items():
        net.add("Load", f"Load{bus}", bus=f"Bus{bus}", p_set=p)

    return net


def time_median(func, n_runs):
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        func()
        times.append(time.perf_counter() - t0)
    return statistics.median(times), min(times)


# ================================================================
#  BENCHMARK
# ================================================================
print("=" * 70)
print("BENCHMARK: Python/PyPSA — AC Power Flow (Newton-Raphson)")
print("=" * 70)

AC_SIZES = [3, 10, 50, 100]

print("\n[AC POWER FLOW BENCHMARK]")
print("-" * 70)
print(f"{'Buses':<10} {'Median (ms)':>12} {'Min (ms)':>12} {'Lines':>10}")
print("-" * 70)

ac_results = {}

for n in AC_SIZES:
    n_buses, lines, generators, loads = generate_network_ac(n, seed=42)
    net = build_pypsa_network_ac(n_buses, lines, generators, loads)

    n_runs = 20 if n <= 10 else (10 if n <= 50 else 5)
    med, mn = time_median(lambda: net.pf(), n_runs)

    ac_results[n] = med * 1000
    print(f"{n:<10} {med*1000:>12.3f} {mn*1000:>12.3f} {len(lines):>10}")

with open("../results/python_ac_benchmark.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["module", "n_buses", "time_ms"])
    for n in AC_SIZES:
        writer.writerow(["AC_PF", n, ac_results[n]])

print("\n[OK] Saved to results/python_ac_benchmark.csv")
print("Run julia/benchmarks/benchmark_ac.jl to get Julia times for comparison.")
