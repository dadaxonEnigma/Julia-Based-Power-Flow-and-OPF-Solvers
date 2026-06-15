"""
compare_validation.py — Diff table: PyPSA vs Julia for all tests.

Reads:
  results/validation/pypsa_validation_full.csv   (from python/validation_full.py)
  results/validation/julia_validation_full.csv   (from julia/validation/validation_full.jl)

Saves:
  results/validation/validation_comparison.csv   — full diff table
  results/validation/validation_summary.csv      — one row per test, pass/warn/fail counts

Usage:
  python python/compare_validation.py
"""

import csv, math

# ── Tolerances per variable type ──────────────────────────────────────────────
TOLS = {
    "theta":      {"warn": 1e-4, "fail": 1e-2},   # bus angles (rad)
    "flow":       {"warn": 0.1,  "fail": 1.0},    # line flows (MW)
    "P_":         {"warn": 0.1,  "fail": 1.0},    # generator output (MW)
    "p_":         {"warn": 0.1,  "fail": 1.0},    # link/storage power (MW)
    "LMP":        {"warn": 0.01, "fail": 0.5},    # locational prices (€/MWh)
    "SoC":        {"warn": 0.5,  "fail": 5.0},    # state of charge (MWh)
    "total_cost": {"warn": 1.0,  "fail": 10.0},   # total cost (€)
    "u_":         {"warn": 0.01, "fail": 0.5},    # commitment binary
    "default":    {"warn": 0.1,  "fail": 1.0},
}

def get_tol(variable):
    for prefix, tol in TOLS.items():
        if variable.startswith(prefix):
            return tol
    return TOLS["default"]

def status(abs_diff, variable):
    tol = get_tol(variable)
    if abs_diff > tol["fail"]:  return "FAIL"
    if abs_diff > tol["warn"]:  return "WARN"
    return "PASS"

# ── Load CSVs ─────────────────────────────────────────────────────────────────
def load_csv(path):
    data = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            key = (row["test_id"], row["variable"], int(row["t"]))
            data[key] = float(row["value"])
    return data

try:
    pypsa_data = load_csv("../results/validation/pypsa_validation_full.csv")
except FileNotFoundError:
    print("[ERROR] pypsa_validation_full.csv not found.")
    print("        Run: python python/validation_full.py")
    exit(1)

try:
    julia_data = load_csv("../results/validation/julia_validation_full.csv")
except FileNotFoundError:
    print("[ERROR] julia_validation_full.csv not found.")
    print("        Run: julia julia/validation/validation_full.jl")
    exit(1)

# ── Build comparison ──────────────────────────────────────────────────────────
all_keys = sorted(set(pypsa_data.keys()) | set(julia_data.keys()))

comparison_rows = []
summary = {}   # test_id → {pass, warn, fail}

for key in all_keys:
    test_id, variable, t = key
    py_val  = pypsa_data.get(key, None)
    jl_val  = julia_data.get(key, None)

    if py_val is None:
        abs_diff = rel_diff = float("nan")
        st = "ONLY_JULIA"
    elif jl_val is None:
        abs_diff = rel_diff = float("nan")
        st = "ONLY_PYPSA"
    else:
        abs_diff = abs(jl_val - py_val)
        denom    = max(abs(py_val), 1e-6)
        rel_diff = abs_diff / denom
        st       = status(abs_diff, variable)

    comparison_rows.append({
        "test_id":  test_id,
        "variable": variable,
        "t":        t,
        "pypsa":    f"{py_val:.6f}" if py_val is not None else "—",
        "julia":    f"{jl_val:.6f}" if jl_val is not None else "—",
        "abs_diff": f"{abs_diff:.2e}" if not math.isnan(abs_diff) else "—",
        "rel_diff": f"{rel_diff:.2e}" if not math.isnan(rel_diff) else "—",
        "status":   st,
    })

    if test_id not in summary:
        summary[test_id] = {"PASS": 0, "WARN": 0, "FAIL": 0, "OTHER": 0}
    bucket = st if st in ("PASS","WARN","FAIL") else "OTHER"
    summary[test_id][bucket] += 1

# ── Save comparison CSV ───────────────────────────────────────────────────────
out_comp = "../results/validation/validation_comparison.csv"
with open(out_comp, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["test_id","variable","t",
                                      "pypsa","julia","abs_diff","rel_diff","status"])
    w.writeheader()
    w.writerows(comparison_rows)

# ── Save summary CSV ──────────────────────────────────────────────────────────
out_sum = "../results/validation/validation_summary.csv"
with open(out_sum, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["test_id","PASS","WARN","FAIL","OTHER","overall"])
    w.writeheader()
    for tid, counts in sorted(summary.items()):
        overall = "FAIL" if counts["FAIL"] > 0 else ("WARN" if counts["WARN"] > 0 else "PASS")
        w.writerow({"test_id": tid, **counts, "overall": overall})

# ── Print console report ──────────────────────────────────────────────────────
W = 80
print("=" * W)
print("  VALIDATION COMPARISON  —  Julia vs PyPSA")
print("=" * W)

prev_test = None
fail_count = 0

for row in comparison_rows:
    if row["test_id"] != prev_test:
        if prev_test is not None:
            print()
        s = summary[row["test_id"]]
        overall = "FAIL" if s["FAIL"] > 0 else ("WARN" if s["WARN"] > 0 else "PASS")
        marker = "✓" if overall == "PASS" else ("?" if overall == "WARN" else "✗")
        print(f"\n  {marker}  {row['test_id']}  "
              f"[PASS={s['PASS']}  WARN={s['WARN']}  FAIL={s['FAIL']}]")
        print(f"  {'Variable':<22} {'t':>3}  {'PyPSA':>14}  {'Julia':>14}  "
              f"{'|Δ|':>10}  Status")
        print("  " + "─" * 72)
        prev_test = row["test_id"]

    # Only print non-PASS rows plus first few PASS per test
    if row["status"] != "PASS":
        marker = "  FAIL !" if row["status"] == "FAIL" else "  WARN ?"
        fail_count += 1
    else:
        marker = "       "

    print(f"  {row['variable']:<22} {row['t']:>3}  {row['pypsa']:>14}  "
          f"{row['julia']:>14}  {row['abs_diff']:>10}  {row['status']}{marker}")

print("\n" + "=" * W)
print(f"  Total variables: {len(comparison_rows)}")
total_fail = sum(1 for r in comparison_rows if r["status"] == "FAIL")
total_warn = sum(1 for r in comparison_rows if r["status"] == "WARN")
total_pass = sum(1 for r in comparison_rows if r["status"] == "PASS")
print(f"  PASS: {total_pass}   WARN: {total_warn}   FAIL: {total_fail}")

print(f"\n  Saved: {out_comp}")
print(f"  Saved: {out_sum}")

if total_fail == 0 and total_warn == 0:
    print("\n  ALL CHECKS PASSED ✓")
elif total_fail == 0:
    print(f"\n  {total_warn} warnings — review WARN rows above")
else:
    print(f"\n  {total_fail} FAILURES — review FAIL rows above")
print("=" * W)
