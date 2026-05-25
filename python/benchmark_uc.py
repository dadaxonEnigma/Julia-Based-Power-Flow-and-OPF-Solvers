"""
Benchmark: Unit Commitment + LMP — PyPSA (linopy + HiGHS)
Same network as julia/benchmarks/benchmark_uc.jl for direct comparison.

Network: 3-bus, G_base (continuous 270 MW), G_peak (committable 150 MW),
         Wind 150 MW (zero cost), Load2=250 MW·profile, Load3=175 MW·profile
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

LOAD_PROFILE_24 = np.array([
    0.60, 0.57, 0.55, 0.54, 0.55, 0.60,
    0.70, 0.80, 0.88, 0.90, 0.92, 0.91,
    0.90, 0.89, 0.88, 0.87, 0.89, 0.95,
    1.00, 0.98, 0.93, 0.85, 0.75, 0.65,
])

WIND_PROFILE_24 = np.array([
    0.80, 0.82, 0.85, 0.83, 0.78, 0.70,
    0.60, 0.55, 0.50, 0.45, 0.42, 0.40,
    0.38, 0.37, 0.40, 0.43, 0.50, 0.58,
    0.65, 0.70, 0.74, 0.76, 0.78, 0.80,
])


def tile_profile(base, T):
    reps = (T // len(base)) + 1
    return np.tile(base, reps)[:T]


def build_network(T):
    load_p = tile_profile(LOAD_PROFILE_24, T)
    wind_cf = tile_profile(WIND_PROFILE_24, T)

    net = pypsa.Network()
    net.set_snapshots(range(T))

    for i in [1, 2, 3]:
        net.add("Bus", f"Bus{i}", v_nom=380.0)

    net.add("Line", "L12", bus0="Bus1", bus1="Bus2", x=0.1, r=0.01, s_nom=1e6)
    net.add("Line", "L13", bus0="Bus1", bus1="Bus3", x=0.1, r=0.01, s_nom=1e6)
    net.add("Line", "L23", bus0="Bus2", bus1="Bus3", x=0.1, r=0.01, s_nom=1e6)

    # Continuous base-load generator (always dispatchable)
    net.add("Generator", "G_base",
            bus="Bus1", p_nom=270.0, marginal_cost=20.0, control="Slack")

    # Committable peaker: binary on/off, min up-time 2 h, startup cost 500 €
    net.add("Generator", "G_peak",
            bus="Bus1", p_nom=150.0, marginal_cost=80.0,
            committable=True,
            min_up_time=2, min_down_time=1,
            start_up_cost=500.0, shut_down_cost=0.0,
            p_min_pu=0.3, p_max_pu=1.0)

    net.add("Load", "Load2", bus="Bus2", p_set=250.0 * load_p)
    net.add("Load", "Load3", bus="Bus3", p_set=175.0 * load_p)

    net.add("Generator", "Wind3",
            bus="Bus3", p_nom=150.0, marginal_cost=0.0,
            p_max_pu=wind_cf)

    return net


def time_median(func, n_runs):
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        func()
        times.append(time.perf_counter() - t0)
    return statistics.median(times), min(times)


print("=" * 65)
print("BENCHMARK: Unit Commitment + LMP  (Python/PyPSA — linopy + HiGHS)")
print("=" * 65)
print(f"Network: 3-bus | G_base=270MW | G_peak=150MW committable | Wind=150MW")
print(f"\n{'T (h)':<8} {'Median (ms)':>12} {'Min (ms)':>12}")
print("-" * 38)

HORIZONS = [6, 12, 24, 48]
uc_results = {}

for T in HORIZONS:
    def run_uc(T=T):
        net = build_network(T)
        net.optimize(solver_name="highs")

    n_runs = 10 if T <= 24 else 5
    med, mn = time_median(run_uc, n_runs)

    uc_results[T] = med * 1000
    print(f"{T:<8} {med*1000:>12.3f} {mn*1000:>12.3f}")

# ── LMP snapshot for T=24 ────────────────────────────────────────────────────
print("\n--- LMP snapshot (T=24) ---")
net24 = build_network(24)
net24.optimize(solver_name="highs")
try:
    lmp = net24.buses_t.marginal_price
    for bus in ["Bus1", "Bus2", "Bus3"]:
        v = lmp[bus].values
        print(f"  {bus}: avg={v.mean():.3f}  min={v.min():.3f}  max={v.max():.3f}  €/MWh")
except Exception as e:
    print(f"  LMP unavailable: {e}")

# ── Save CSV ─────────────────────────────────────────────────────────────────
with open("../results/python_uc_benchmark.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["module", "T", "time_ms"])
    for T in HORIZONS:
        writer.writerow(["UC", T, uc_results[T]])

print("\n[OK] Saved to results/python_uc_benchmark.csv")
print("Run julia/benchmarks/benchmark_uc.jl for Julia comparison.")
