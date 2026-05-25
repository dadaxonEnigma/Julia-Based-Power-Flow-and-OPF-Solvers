"""
Compare Julia vs PyPSA for UC scaling (3, 14, 30 buses, T=24).
Reads results/julia_uc_scale.csv and results/python_uc_scale.csv.
"""
import csv
import os


def read_scale_csv(path):
    results = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            results[int(row["n_buses"])] = float(row["time_ms"])
    return results


julia_path  = "../results/julia_uc_scale.csv"
python_path = "../results/python_uc_scale.csv"

for p in [julia_path, python_path]:
    if not os.path.exists(p):
        print(f"Missing: {p}")
        print("Run benchmark_uc_scale.jl and benchmark_uc_scale.py first.")
        raise SystemExit(1)

julia  = read_scale_csv(julia_path)
python = read_scale_csv(python_path)

print()
print("=" * 66)
print("  UC Scaling (T=24)  —  Julia (JuMP+HiGHS) vs Python (PyPSA+HiGHS)")
print("=" * 66)
print(f"{'Buses':<8} {'Julia (ms)':>12} {'Python (ms)':>13} {'Speedup':>10}")
print("-" * 66)

buses = sorted(set(julia) | set(python))
for n in buses:
    j = julia.get(n)
    p = python.get(n)
    if j and p:
        speedup = p / j
        bar = "#" * min(30, int(speedup))
        print(f"{n:<8} {j:>12.3f} {p:>13.3f} {speedup:>9.1f}x  {bar}")
    elif j:
        print(f"{n:<8} {j:>12.3f} {'—':>13}")
    else:
        print(f"{n:<8} {'—':>12} {p:>13.3f}")

print()
j_vals = [julia[n] for n in buses if n in julia]
p_vals = [python[n] for n in buses if n in python]
if len(j_vals) == len(p_vals) == len(buses):
    avg_speedup = sum(p / j for j, p in zip(j_vals, p_vals)) / len(buses)
    print(f"  Average speedup across all sizes: {avg_speedup:.1f}x")
