"""
check_dc_pf.py — DC Power Flow reference values (PyPSA).
Run side-by-side with julia/validation/check_dc_pf.jl

Network:  3-bus (Bus1=slack, Bus2, Bus3)
          G1=400MW at Bus1, D2=200MW at Bus2, D3=300MW at Bus3
          L12: x=0.1, L13: x=0.1, L23: x=0.15  (all r=0.01, s_nom=inf)
"""
import pypsa
import warnings
warnings.filterwarnings("ignore")

def make_net():
    net = pypsa.Network()
    net.set_snapshots([0])
    for b in ["Bus1", "Bus2", "Bus3"]:
        net.add("Bus", b, v_nom=380.0)
    net.add("Line", "L12", bus0="Bus1", bus1="Bus2", r=0.01, x=0.1,  s_nom=1e6)
    net.add("Line", "L13", bus0="Bus1", bus1="Bus3", r=0.01, x=0.1,  s_nom=1e6)
    net.add("Line", "L23", bus0="Bus2", bus1="Bus3", r=0.01, x=0.15, s_nom=1e6)
    net.add("Generator", "G1", bus="Bus1", p_nom=400.0, marginal_cost=20.0, control="Slack")
    net.add("Generator", "G2", bus="Bus2", p_nom=300.0, marginal_cost=50.0)
    net.add("Load", "D2", bus="Bus2", p_set=200.0)
    net.add("Load", "D3", bus="Bus3", p_set=300.0)
    return net

W = 52

# ── Test 1: DC Power Flow ──────────────────────────────────────────────────────
print("=" * W)
print("  TEST 1 — DC Power Flow (lpf)")
print("=" * W)
net = make_net()
net.lpf()

angles = net.buses_t.v_ang.iloc[0]  # radians
print(f"  θ Bus1 (slack) : {angles['Bus1']:>12.6f} rad")
print(f"  θ Bus2         : {angles['Bus2']:>12.6f} rad")
print(f"  θ Bus3         : {angles['Bus3']:>12.6f} rad")
print()

for line in ["L12", "L13", "L23"]:
    flow = net.lines_t.p0[line].iloc[0]
    print(f"  Flow {line}       : {flow:>12.4f} MW")
print()
print(f"  G1 dispatch    : {net.generators_t.p['G1'].iloc[0]:>12.4f} MW  (= total load = 500 MW)")

# ── Test 2: LOPF uncongested ───────────────────────────────────────────────────
print("\n" + "=" * W)
print("  TEST 2 — LOPF (no congestion, s_nom=∞)")
print("=" * W)
net = make_net()
net.optimize(solver_name="highs")

print(f"  P_G1           : {net.generators_t.p['G1'].iloc[0]:>12.4f} MW")
print(f"  P_G2           : {net.generators_t.p['G2'].iloc[0]:>12.4f} MW")
print(f"  Total cost     : {net.objective:>12.4f} €")
print()
for b in ["Bus1", "Bus2", "Bus3"]:
    lmp = net.buses_t.marginal_price[b].iloc[0]
    print(f"  LMP {b}     : {lmp:>12.4f} €/MWh")

# ── Test 3: LOPF congested ────────────────────────────────────────────────────
print("\n" + "=" * W)
print("  TEST 3 — LOPF (s_nom=200 MW on all lines)")
print("=" * W)
net2 = pypsa.Network()
net2.set_snapshots([0])
for b in ["Bus1", "Bus2", "Bus3"]:
    net2.add("Bus", b, v_nom=380.0)
net2.add("Line", "L12", bus0="Bus1", bus1="Bus2", r=0.01, x=0.1,  s_nom=200.0)
net2.add("Line", "L13", bus0="Bus1", bus1="Bus3", r=0.01, x=0.1,  s_nom=200.0)
net2.add("Line", "L23", bus0="Bus2", bus1="Bus3", r=0.01, x=0.15, s_nom=200.0)
net2.add("Generator", "G1", bus="Bus1", p_nom=400.0, marginal_cost=20.0, control="Slack")
net2.add("Generator", "G2", bus="Bus2", p_nom=300.0, marginal_cost=50.0)
net2.add("Load",      "D2", bus="Bus2", p_set=200.0)
net2.add("Load",      "D3", bus="Bus3", p_set=300.0)
net2.optimize(solver_name="highs")

print(f"  P_G1           : {net2.generators_t.p['G1'].iloc[0]:>12.4f} MW")
print(f"  P_G2           : {net2.generators_t.p['G2'].iloc[0]:>12.4f} MW")
print(f"  Total cost     : {net2.objective:>12.4f} €")
print()
for b in ["Bus1", "Bus2", "Bus3"]:
    lmp = net2.buses_t.marginal_price[b].iloc[0]
    print(f"  LMP {b}     : {lmp:>12.4f} €/MWh")

# ── Test 4: Multi-period LOPF T=6 ─────────────────────────────────────────────
print("\n" + "=" * W)
print("  TEST 4 — Multi-Period LOPF (T=6, StorageUnit)")
print("=" * W)
import numpy as np
LP6 = np.array([0.60, 0.57, 0.55, 0.54, 0.55, 0.60])
net3 = pypsa.Network()
net3.set_snapshots(range(6))
for b in ["Bus1", "Bus2", "Bus3"]:
    net3.add("Bus", b, v_nom=380.0)
net3.add("Line", "L12", bus0="Bus1", bus1="Bus2", r=0.01, x=0.1,  s_nom=1e6)
net3.add("Line", "L13", bus0="Bus1", bus1="Bus3", r=0.01, x=0.1,  s_nom=1e6)
net3.add("Line", "L23", bus0="Bus2", bus1="Bus3", r=0.01, x=0.15, s_nom=1e6)
net3.add("Generator", "G1", bus="Bus1", p_nom=400.0, marginal_cost=20.0, control="Slack")
net3.add("Generator", "G2", bus="Bus2", p_nom=300.0, marginal_cost=50.0)
net3.add("Load",      "D2", bus="Bus2", p_set=200.0 * LP6)
net3.add("Load",      "D3", bus="Bus3", p_set=300.0 * LP6)
net3.add("StorageUnit", "Bat", bus="Bus2", p_nom=60.0,
         max_hours=4.0, efficiency_store=0.9, efficiency_dispatch=0.9,
         cyclic_state_of_charge=False, state_of_charge_initial=120.0)
net3.optimize(solver_name="highs")

print(f"  Total cost     : {net3.objective:>12.4f} €")
print(f"\n  {'t':<4} {'P_G1':>10} {'P_G2':>10} {'SoC_Bat':>10} {'LMP_Bus1':>10}")
print("  " + "─"*46)
for t in range(6):
    g1  = net3.generators_t.p["G1"].iloc[t]
    g2  = net3.generators_t.p["G2"].iloc[t]
    soc = net3.storage_units_t.state_of_charge["Bat"].iloc[t]
    lmp = net3.buses_t.marginal_price["Bus1"].iloc[t]
    print(f"  {t:<4} {g1:>10.3f} {g2:>10.3f} {soc:>10.3f} {lmp:>10.4f}")

# ── Test 5: Unit Commitment T=6 ───────────────────────────────────────────────
print("\n" + "=" * W)
print("  TEST 5 — Unit Commitment (T=6, G_peak committable)")
print("=" * W)
net4 = pypsa.Network()
net4.set_snapshots(range(6))
for b in ["Bus1", "Bus2", "Bus3"]:
    net4.add("Bus", b, v_nom=380.0)
net4.add("Line", "L12", bus0="Bus1", bus1="Bus2", r=0.01, x=0.1,  s_nom=1e6)
net4.add("Line", "L13", bus0="Bus1", bus1="Bus3", r=0.01, x=0.1,  s_nom=1e6)
net4.add("Line", "L23", bus0="Bus2", bus1="Bus3", r=0.01, x=0.15, s_nom=1e6)
net4.add("Generator", "G1", bus="Bus1", p_nom=400.0, marginal_cost=20.0, control="Slack")
net4.add("Generator", "G_peak", bus="Bus2", p_nom=300.0, marginal_cost=50.0,
         committable=True, min_up_time=2, min_down_time=1,
         start_up_cost=1000.0, p_min_pu=0.3, p_max_pu=1.0, initial_status=0)
net4.add("Load", "D2", bus="Bus2", p_set=200.0 * LP6)
net4.add("Load", "D3", bus="Bus3", p_set=300.0 * LP6)
net4.optimize(solver_name="highs")

print(f"  Total cost     : {net4.objective:>12.4f} €")
print(f"\n  {'t':<4} {'P_G1':>10} {'P_G_peak':>10} {'u_G_peak':>10}")
print("  " + "─"*38)
for t in range(6):
    g1   = net4.generators_t.p["G1"].iloc[t]
    gpk  = net4.generators_t.p["G_peak"].iloc[t]
    try:
        u = net4.generators_t["status"]["G_peak"].iloc[t]
    except:
        u = float("nan")
    print(f"  {t:<4} {g1:>10.3f} {gpk:>10.3f} {u:>10.0f}")

print("\n" + "=" * W)
print("  Done. Compare with:  julia julia/validation/check_dc_pf.jl")
print("=" * W)
