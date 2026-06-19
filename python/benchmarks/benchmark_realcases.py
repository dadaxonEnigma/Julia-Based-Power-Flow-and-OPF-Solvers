"""
benchmark_realcases.py — PyPSA counterpart of julia/benchmarks/benchmark_realcases.jl.

Loads the same standard MATPOWER cases (IEEE 118/300, PEGASE 1354/2869) into a
PyPSA network and times the same operations as the Julia side:
  DC_PF -> net.lpf()                      (linear power flow)
  LOPF  -> net.optimize(solver='highs')

Same modeling choices as the Julia loader (branches as lines with x*tap, loads
from bus Pd, generators with Pmax and the gencost linear term), so the timing
comparison is apples-to-apples on identical problem dimensions.

Run:
  pypsa_env/Scripts/python.exe python/benchmark_realcases.py
"""
import os
# Fair single-thread linear algebra (must be set before numpy is imported).
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import re, time, csv, statistics, warnings, logging
import pypsa

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

DATA = os.path.join(os.path.dirname(__file__), "..", "..", "data")
CASES = ["case118.m", "case300.m", "case1354pegase.m", "case2869pegase.m"]


def parse_m(path):
    txt = open(path).read()

    def matrix(name):
        m = re.search(r"mpc\." + name + r"\s*=\s*\[(.*?)\]", txt, re.S)
        if not m:
            return []
        rows = []
        for line in m.group(1).strip().split("\n"):
            line = line.split("%")[0].strip().rstrip(";").strip()
            if not line:
                continue
            rows.append([float(p) for p in re.split(r"\s+", line)])
        return rows

    base = float(re.search(r"mpc\.baseMVA\s*=\s*([\d.]+)", txt).group(1))
    return base, matrix("bus"), matrix("gen"), matrix("branch"), matrix("gencost")


def build_network(path):
    base, bus, gen, branch, gencost = parse_m(path)
    n = pypsa.Network()
    ref = next((int(r[0]) for r in bus if int(r[1]) == 3), int(bus[0][0]))

    for r in bus:
        n.add("Bus", f"b{int(r[0])}", v_nom=380.0)
    for r in bus:                                   # loads from bus Pd (col 3)
        if abs(r[2]) > 1e-9:
            n.add("Load", f"d{int(r[0])}", bus=f"b{int(r[0])}", p_set=r[2])
    for k, r in enumerate(branch):                  # branches -> lines
        if len(r) >= 11 and r[10] == 0:
            continue
        tap = r[8] if len(r) > 8 and r[8] != 0 else 1.0
        rate = r[5] if len(r) > 5 else 0.0
        n.add("Line", f"l{k}", bus0=f"b{int(r[0])}", bus1=f"b{int(r[1])}",
              x=r[3] * tap, r=r[2], s_nom=(rate if rate > 0 else 1e6))
    for k, r in enumerate(gen):                     # generators
        if len(r) > 7 and r[7] == 0:
            continue
        pmax = r[8]
        if pmax <= 0:
            continue
        cost = 20.0
        if k < len(gencost):
            gc = gencost[k]
            if len(gc) >= 4 and int(gc[0]) == 2:
                coeffs = gc[4:4 + int(gc[3])]
                if len(coeffs) >= 2 and coeffs[-2] > 0:
                    cost = coeffs[-2]
        ctrl = "Slack" if int(r[0]) == ref else "PQ"
        n.add("Generator", f"g{k}", bus=f"b{int(r[0])}",
              p_nom=pmax, marginal_cost=cost, control=ctrl)
    return n


def timed(fn, n=7):
    fn()                                            # warm-up
    ts = []
    for _ in range(n):
        t0 = time.perf_counter()
        fn()
        ts.append((time.perf_counter() - t0) * 1000)
    return statistics.mean(ts), (statistics.stdev(ts) if len(ts) > 1 else 0.0), min(ts)


rows = []
print("=" * 70)
print("  Real-case benchmark — PyPSA on standard IEEE/PEGASE networks")
print("=" * 70)
for c in CASES:
    path = os.path.join(DATA, c)
    if not os.path.isfile(path):
        print(f"  {c:18s} SKIP (not found)"); continue
    net = build_network(path)
    nb = len(net.buses)
    print(f"\n[{c}]  {nb} buses, {len(net.lines)} lines, "
          f"{len(net.generators)} gens, {len(net.loads)} loads")
    for label, fn in (("DC_PF", lambda net=net: net.lpf()),
                      ("LOPF",  lambda net=net: net.optimize(
                          solver_name="highs", solver_options={"threads": 1}))):
        try:
            m, s, mn = timed(fn)
            rows.append((c, nb, label, m, s, mn))
            print(f"    {label:6s} mean={m:.2f} ms  std={s:.2f}  min={mn:.2f}")
        except Exception as e:
            print(f"    {label:6s} FAILED: {e}")

out = os.path.join(os.path.dirname(__file__), "..", "..", "results", "benchmarks", "pypsa_realcases_benchmark.csv")
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["case", "n_buses", "module", "time_ms", "std_ms", "min_ms"])
    for r in rows:
        w.writerow([r[0], r[1], r[2], f"{r[3]:.6f}", f"{r[4]:.6f}", f"{r[5]:.6f}"])
print(f"\n[OK] {len(rows)} rows -> {out}")
