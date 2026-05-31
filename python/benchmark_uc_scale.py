"""
Benchmark: Unit Commitment scaling — PyPSA (linopy + HiGHS)
Tests UC + LMP at 3, 14, and 30 buses, T=24.
Same deterministic network generator as julia/benchmarks/benchmark_uc_scale.jl.

Compare with julia/benchmarks/benchmark_uc_scale.jl.
"""
import csv
import time
import statistics
import warnings
import logging

import numpy as np
import pypsa

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

# ── Shared 24-hour profile ────────────────────────────────────────────────────
LOAD_PROFILE = np.array([
    0.60, 0.57, 0.55, 0.54, 0.55, 0.60,
    0.70, 0.80, 0.88, 0.90, 0.92, 0.91,
    0.90, 0.89, 0.88, 0.87, 0.89, 0.95,
    1.00, 0.98, 0.93, 0.85, 0.75, 0.65,
])

# ── Deterministic network generator (matches Julia build_uc_network) ──────────
def build_uc_network(n_buses, seed=42):
    """
    Deterministic random network matching Julia's build_uc_network.
    Uses numpy default_rng with the same sequence logic.
    """
    rng = np.random.default_rng(seed)

    lines   = []
    added   = set()
    # Spanning tree
    for i in range(1, n_buses):
        x = 0.05 + rng.random() * 0.45
        lines.append((i, i + 1, x))
        added.add((i, i + 1))
    # Extra cross-edges
    for _ in range(max(1, n_buses // 3)):
        u = int(rng.integers(1, n_buses))
        v = int(rng.integers(u + 1, n_buses + 1))
        if (u, v) not in added:
            x = 0.05 + rng.random() * 0.45
            lines.append((u, v, x))
            added.add((u, v))

    # Loads
    total_load = 0.0
    loads = {}
    for bus in range(2, n_buses + 1):
        if rng.random() > 0.3:
            p = 50.0 + rng.random() * 200.0
            loads[bus] = p
            total_load += p
    if not loads:
        loads[2] = 200.0
        total_load = 200.0

    # Generator assignment
    gen_params = {}
    # Base at bus 1
    gen_params[1] = dict(p_nom=total_load * 1.3, mc=20.0, commit=False)
    for bus in range(2, n_buses + 1, 4):
        p_cap = 80.0 + rng.random() * 120.0
        mc    = 50.0 + rng.random() * 40.0
        commit = rng.random() < 0.4
        gen_params[bus] = dict(p_nom=p_cap, mc=mc, commit=commit)

    return lines, loads, gen_params


def build_pypsa(n_buses, T, seed=42):
    lines, loads, gen_params = build_uc_network(n_buses, seed)
    net = pypsa.Network()
    net.set_snapshots(range(T))

    for i in range(1, n_buses + 1):
        net.add("Bus", f"B{i}", v_nom=380.0)

    for idx, (f, t, x) in enumerate(lines):
        net.add("Line", f"L{idx}", bus0=f"B{f}", bus1=f"B{t}",
                x=x, r=0.01, s_nom=1e6)

    for bus, gp in gen_params.items():
        control = "Slack" if bus == 1 else "PQ"
        if gp["commit"]:
            net.add("Generator", f"G_pk{bus}", bus=f"B{bus}",
                    p_nom=gp["p_nom"], marginal_cost=gp["mc"],
                    control=control,
                    committable=True,
                    min_up_time=2, min_down_time=1,
                    start_up_cost=500.0, shut_down_cost=0.0,
                    p_min_pu=0.3, p_max_pu=1.0)
        else:
            net.add("Generator", f"G_ct{bus}", bus=f"B{bus}",
                    p_nom=gp["p_nom"], marginal_cost=gp["mc"],
                    control=control)

    for bus, p_base in loads.items():
        net.add("Load", f"D{bus}", bus=f"B{bus}",
                p_set=p_base * LOAD_PROFILE[:T])

    return net


def time_median(func, n_runs):
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        func()
        times.append(time.perf_counter() - t0)
    return statistics.median(times), min(times)


# ── Benchmark ─────────────────────────────────────────────────────────────────
T_FIXED = 24
SIZES   = [3, 14, 30, 40]

print("=" * 70)
print("BENCHMARK: Unit Commitment Scaling  (Python/PyPSA — linopy + HiGHS, T=24)")
print("=" * 70)
print(f"{'n_buses':<8}  {'T':<4}  {'Median (ms)':>12}  {'Min (ms)':>12}")
print("-" * 70)

uc_results = {}

for n in SIZES:
    def run_uc(n=n):
        net = build_pypsa(n, T_FIXED)
        net.optimize(solver_name="highs")

    n_runs = 10 if n <= 10 else 5
    med, mn = time_median(run_uc, n_runs)

    uc_results[n] = med * 1000
    print(f"{n:<8}  {T_FIXED:<4}  {med*1000:>12.3f}  {mn*1000:>12.3f}")

# ── Save CSV ───────────────────────────────────────────────────────────────────
with open("../results/python_uc_scale.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["module", "n_buses", "T", "time_ms"])
    for n in SIZES:
        writer.writerow(["UC", n, T_FIXED, uc_results[n]])

print("\n[OK] Saved to results/python_uc_scale.csv")
print("Run julia/benchmarks/benchmark_uc_scale.jl for Julia comparison.")
