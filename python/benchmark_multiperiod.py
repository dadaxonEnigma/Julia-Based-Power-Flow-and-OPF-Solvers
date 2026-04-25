"""
Benchmark: Multi-Period LOPF — PyPSA (linopy + HiGHS)
Same network and parameters as julia/benchmarks/benchmark_multiperiod.jl
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

LOAD_PROFILE_FULL = np.array([
    0.60, 0.57, 0.55, 0.54, 0.55, 0.60,
    0.70, 0.80, 0.88, 0.90, 0.92, 0.91,
    0.90, 0.89, 0.88, 0.87, 0.89, 0.95,
    1.00, 0.98, 0.93, 0.85, 0.75, 0.65,
])

WIND_PROFILE_FULL = np.array([
    0.80, 0.82, 0.85, 0.83, 0.78, 0.70,
    0.60, 0.55, 0.50, 0.45, 0.42, 0.40,
    0.38, 0.37, 0.40, 0.43, 0.50, 0.58,
    0.65, 0.70, 0.74, 0.76, 0.78, 0.80,
])


def get_profile(T):
    """Tile the 24-h profile to length T."""
    reps = (T // 24) + 1
    return (np.tile(LOAD_PROFILE_FULL, reps)[:T],
            np.tile(WIND_PROFILE_FULL, reps)[:T])


def build_network(T):
    net = pypsa.Network()
    load_p, wind_cf = get_profile(T)
    net.set_snapshots(range(T))

    for i in [1, 2, 3]:
        net.add("Bus", f"Bus{i}", v_nom=380.0)

    net.add("Line", "L12", bus0="Bus1", bus1="Bus2", x=0.1, r=0.01, s_nom=1e6)
    net.add("Line", "L13", bus0="Bus1", bus1="Bus3", x=0.1, r=0.01, s_nom=1e6)
    net.add("Line", "L23", bus0="Bus2", bus1="Bus3", x=0.1, r=0.01, s_nom=1e6)

    net.add("Generator", "G1", bus="Bus1", p_nom=270.0,
            marginal_cost=20.0, control="Slack")
    net.add("Generator", "G2", bus="Bus2", p_nom=100.0, marginal_cost=50.0)

    net.add("Load", "Load2", bus="Bus2", p_set=250.0 * load_p)
    net.add("Load", "Load3", bus="Bus3", p_set=175.0 * load_p)

    net.add("Generator", "Wind3", bus="Bus3",
            p_nom=150.0, marginal_cost=0.0, p_max_pu=wind_cf)

    net.add("StorageUnit", "Stor2", bus="Bus2",
            p_nom=100.0, max_hours=3.0,
            efficiency_store=0.95, efficiency_dispatch=0.95,
            cyclic_state_of_charge=True,
            state_of_charge_initial=150.0)
    return net


def time_median(func, n_runs):
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        func()
        times.append(time.perf_counter() - t0)
    return statistics.median(times), min(times)


print("=" * 65)
print("BENCHMARK: Multi-Period LOPF  (Python/PyPSA — linopy + HiGHS)")
print("=" * 65)
print(f"{'T (h)':<8} {'Median (ms)':>12} {'Min (ms)':>12} {'Vars':>8}")
print("-" * 45)

HORIZONS = [6, 12, 24, 48, 96]
mp_results = {}

for T in HORIZONS:
    net = build_network(T)
    n_vars = T * 8  # approximate

    n_runs = 20 if T <= 24 else (10 if T <= 48 else 5)
    med, mn = time_median(lambda: net.optimize(solver_name="highs"), n_runs)

    mp_results[T] = med * 1000
    print(f"{T:<8} {med*1000:>12.3f} {mn*1000:>12.3f} {n_vars:>8}")

with open("../results/python_mp_benchmark.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["module", "T", "time_ms"])
    for T in HORIZONS:
        writer.writerow(["MLOPF", T, mp_results[T]])

print("\n[OK] Saved to results/python_mp_benchmark.csv")
