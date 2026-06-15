"""
Читает results/benchmarks/julia_benchmark.csv и results/benchmarks/python_benchmark.csv,
выводит сводную таблицу со speedup Julia vs Python.
"""
import csv

def read_csv(path):
    results = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            key = (row["module"], int(row["n_buses"]))
            results[key] = float(row["time_ms"])
    return results

julia  = read_csv("results/benchmarks/julia_benchmark.csv")
python = read_csv("results/benchmarks/python_benchmark.csv")

for module in ["DC_PF", "LOPF"]:
    sizes = sorted({n for (m, n) in julia if m == module})
    print(f"\n{'='*68}")
    print(f"  {module}  —  Julia vs Python/PyPSA")
    print(f"{'='*68}")
    print(f"{'Buses':<8} {'Julia (ms)':>12} {'Python (ms)':>13} {'Speedup':>10}")
    print("-"*68)
    for n in sizes:
        j = julia.get((module, n))
        p = python.get((module, n))
        if j and p:
            speedup = p / j
            bar = "#" * min(40, int(speedup))
            print(f"{n:<8} {j:>12.3f} {p:>13.3f} {speedup:>9.1f}x  {bar}")
        elif j:
            print(f"{n:<8} {j:>12.3f} {'—':>13}")
