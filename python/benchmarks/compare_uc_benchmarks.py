"""
Compare Julia vs PyPSA for Unit Commitment + LMP.
Reads results/benchmarks/julia_uc_benchmark.csv and results/benchmarks/python_uc_benchmark.csv.
"""
import csv


def read_uc_csv(path):
    results = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            results[int(row["T"])] = float(row["time_ms"])
    return results


try:
    julia  = read_uc_csv("../results/benchmarks/julia_uc_benchmark.csv")
    python = read_uc_csv("../results/benchmarks/python_uc_benchmark.csv")
except FileNotFoundError as e:
    print(f"Missing file: {e}")
    print("Run julia/benchmarks/benchmark_uc.jl and python/benchmark_uc.py first.")
    raise SystemExit(1)

print()
print("=" * 62)
print("  Unit Commitment + LMP  —  Julia vs Python/PyPSA")
print("=" * 62)
print(f"{'T (h)':<8} {'Julia (ms)':>12} {'Python (ms)':>13} {'Speedup':>10}")
print("-" * 62)

for T in sorted(set(julia) | set(python)):
    j = julia.get(T)
    p = python.get(T)
    if j and p:
        speedup = p / j
        bar = "#" * min(35, int(speedup))
        print(f"{T:<8} {j:>12.3f} {p:>13.3f} {speedup:>9.1f}x  {bar}")
    elif j:
        print(f"{T:<8} {j:>12.3f} {'—':>13}")
    else:
        print(f"{T:<8} {'—':>12} {p:>13.3f}")
