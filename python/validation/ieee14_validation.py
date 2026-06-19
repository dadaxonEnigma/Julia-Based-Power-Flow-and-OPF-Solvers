"""
IEEE 14-Bus validation: PyPSA vs MATPOWER reference.
Reads data/case14.m directly — no pandapower required.
"""
import re
import numpy as np
import pypsa
import warnings
import logging

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

# ----------------------------------------------------------------
#  Parse MATPOWER .m file
# ----------------------------------------------------------------
def parse_matpower(path):
    with open(path) as f:
        text = f.read()

    def extract_matrix(name):
        pat = rf"mpc\.{name}\s*=\s*\[(.*?)\]"
        m = re.search(pat, text, re.DOTALL)
        if not m:
            return np.array([])
        rows = []
        for line in m.group(1).strip().split('\n'):
            line = line.strip().rstrip(';').strip()
            if line and not line.startswith('%'):
                rows.append([float(x) for x in line.split()])
        return np.array(rows)

    bus    = extract_matrix("bus")
    gen    = extract_matrix("gen")
    branch = extract_matrix("branch")
    gencost = extract_matrix("gencost")
    return bus, gen, branch, gencost


# ----------------------------------------------------------------
#  Build PyPSA network from parsed MATPOWER data
# ----------------------------------------------------------------
def build_pypsa_ieee(bus, gen, branch, gencost, base_mva=100.0):
    net = pypsa.Network()
    net.set_snapshots([0])

    # Buses
    bus_type_map = {1: "PQ", 2: "PV", 3: "Slack"}
    for row in bus:
        i     = int(row[0])
        btype = int(row[1])
        net.add("Bus", f"Bus{i}",
                v_nom=1.0,          # work in per-unit
                control=bus_type_map.get(btype, "PQ"))

    # All branches modelled as Lines in per-unit.
    # Transformer tap ratios are ignored (DC PF / linearised AC approximation).
    # This matches the DC PF assumption and gives consistent results with Julia.
    for k, row in enumerate(branch):
        f, t = int(row[0]), int(row[1])
        r, x, b = row[2], row[3], row[4]
        s_nom = row[5] if row[5] > 0 else 9999.0
        net.add("Line", f"L{k+1}",
                bus0=f"Bus{f}", bus1=f"Bus{t}",
                r=r, x=max(x, 1e-6), b=b, s_nom=s_nom)

    # Loads
    load_idx = 1
    for row in bus:
        i  = int(row[0])
        pd = row[2]   # MW
        qd = row[3]   # MVAr
        if pd != 0 or qd != 0:
            net.add("Load", f"Load{load_idx}",
                    bus=f"Bus{i}",
                    p_set=pd / base_mva,
                    q_set=qd / base_mva)
            load_idx += 1

    # Set voltage magnitude setpoints for PV buses
    for row in bus:
        i     = int(row[0])
        btype = int(row[1])
        vm    = row[7]   # Vm setpoint from case
        if btype in (2, 3):
            net.buses.loc[f"Bus{i}", "v_mag_pu_set"] = vm

    # Generators
    for k, row in enumerate(gen):
        bus_i = int(row[0])
        pg    = row[1] / base_mva
        pmax  = row[8] / base_mva
        pmin  = row[9] / base_mva
        ctrl  = "Slack" if bus_i == 1 else "PV"

        mc = float(gencost[k, -2]) if k < len(gencost) else 1.0

        net.add("Generator", f"G{k+1}",
                bus=f"Bus{bus_i}",
                p_nom=pmax,
                p_set=pg,
                p_min_pu=pmin/pmax if pmax > 0 else 0,
                marginal_cost=mc,
                control=ctrl)

    return net


# ----------------------------------------------------------------
#  Reference solution from MATPOWER case14.m (stored solved state)
# ----------------------------------------------------------------
REF_VM = [1.060, 1.045, 1.010, 1.019, 1.020,
          1.070, 1.062, 1.090, 1.056, 1.051,
          1.057, 1.055, 1.050, 1.036]

REF_VA_DEG = [0.00, -4.98, -12.72, -10.33, -8.78,
              -14.22, -13.37, -13.36, -14.94, -15.10,
              -14.79, -15.07, -15.16, -16.04]

# ================================================================
CASE14 = "../data/case14.m"
bus, gen, branch, gencost = parse_matpower(CASE14)
net = build_pypsa_ieee(bus, gen, branch, gencost)

SEP = "=" * 60

# ----------------------------------------------------------------
#  1. AC Power Flow
# ----------------------------------------------------------------
print(SEP)
print("IEEE 14-Bus: AC Power Flow  (PyPSA Newton-Raphson)")
print(SEP)

net.pf()

buses = [f"Bus{i}" for i in range(1, 15)]
vm = net.buses_t.v_mag_pu.loc[0, buses].values
va = np.degrees(net.buses_t.v_ang.loc[0, buses].values)

print(f"\n{'Bus':>4}  {'Vm PyPSA':>9}  {'Vm Ref':>8}  {'|ΔVm|':>8}  "
      f"{'Va PyPSA':>9}  {'Va Ref':>8}  {'|ΔVa|':>8}")
print("-" * 70)
for i in range(14):
    dvm = abs(vm[i] - REF_VM[i])
    dva = abs(va[i] - REF_VA_DEG[i])
    print(f"{i+1:>4}  {vm[i]:>9.4f}  {REF_VM[i]:>8.4f}  {dvm:>8.2e}  "
          f"{va[i]:>9.4f}  {REF_VA_DEG[i]:>8.4f}  {dva:>8.2e}")

max_dvm = max(abs(vm[i] - REF_VM[i]) for i in range(14))
max_dva = max(abs(va[i] - REF_VA_DEG[i]) for i in range(14))
print(f"\nMax |ΔVm| = {max_dvm:.2e} p.u.   Max |ΔVa| = {max_dva:.2e} deg")

# ----------------------------------------------------------------
#  2. DC Power Flow
# ----------------------------------------------------------------
print(f"\n{SEP}")
print("IEEE 14-Bus: DC Power Flow  (PyPSA lpf)")
print(SEP)

net.lpf()
va_dc = np.degrees(net.buses_t.v_ang.loc[0, buses].values)

print(f"\n{'Bus':>4}  {'θ PyPSA (deg)':>13}  {'θ Ref (deg)':>11}  {'|Δθ|':>8}")
print("-" * 45)
for i in range(14):
    dθ = abs(va_dc[i] - REF_VA_DEG[i])
    print(f"{i+1:>4}  {va_dc[i]:>13.4f}  {REF_VA_DEG[i]:>11.4f}  {dθ:>8.4f}")

max_dθ = max(abs(va_dc[i] - REF_VA_DEG[i]) for i in range(14))
print(f"\nMax |Δθ| = {max_dθ:.4f} deg  (DC vs full-AC reference — expected ~0.1–0.5°)")

# ----------------------------------------------------------------
#  3. LOPF
# ----------------------------------------------------------------
print(f"\n{SEP}")
print("IEEE 14-Bus: LOPF  (PyPSA + HiGHS)")
print(SEP)

net.optimize(solver_name="highs")

print(f"\n{'Generator':>11}  {'P (MW)':>10}  {'Cost (€/MWh)':>13}")
print("-" * 38)
total_cost = 0.0
for g in net.generators.index:
    p  = net.generators_t.p.loc[0, g] * 100   # p.u. → MW
    mc = net.generators.loc[g, "marginal_cost"]
    total_cost += mc * p
    print(f"{g:>11}  {p:>10.2f}  {mc:>13.4f}")

total_load = net.loads.p_set.sum() * 100
print(f"\nTotal cost: {total_cost:.2f} €/h")
print(f"Total load: {total_load:.2f} MW")

print(f"\n{SEP}")
print("IEEE 14-Bus validation complete.")
print(SEP)
