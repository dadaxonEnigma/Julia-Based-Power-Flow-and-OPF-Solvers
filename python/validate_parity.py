"""
Numerical parity validation: PyPSA reference results.
Saves results/pypsa_parity.csv — read by julia/validation/validate_parity.jl.

Cases
-----
  LOPF_uncon   : 3-bus LOPF, no line limits
  LOPF_con     : 3-bus LOPF, s_nom=200 MW on all lines
  MP_t6        : 3-bus multi-period LOPF T=6 with StorageUnit (no wind)
  UC_t6        : 3-bus UC T=6 with one committable generator
"""

import csv
import numpy as np
import pypsa
import warnings
import logging

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

# ── Shared profiles ────────────────────────────────────────────────────────────
LOAD_PROFILE_6 = np.array([0.60, 0.57, 0.55, 0.54, 0.55, 0.60])


# ── Network builders ──────────────────────────────────────────────────────────
def base_3bus(s_nom=1e6, snapshots=None):
    """3-bus network matching Julia test suite make_3bus()."""
    net = pypsa.Network()
    if snapshots is not None:
        net.set_snapshots(snapshots)
    else:
        net.set_snapshots([0])

    for b in ["Bus1", "Bus2", "Bus3"]:
        net.add("Bus", b, v_nom=380.0)

    for f, t, name in [("Bus1", "Bus2", "L12"),
                        ("Bus1", "Bus3", "L13"),
                        ("Bus2", "Bus3", "L23")]:
        net.add("Line", name, bus0=f, bus1=t, r=0.01, x=0.1, s_nom=s_nom)

    net.add("Generator", "G1", bus="Bus1", p_nom=400.0,
            marginal_cost=20.0, control="Slack")
    net.add("Generator", "G2", bus="Bus2", p_nom=300.0, marginal_cost=50.0)
    return net


# ── Row helpers ───────────────────────────────────────────────────────────────
rows = []

def save(case, variable, t, value):
    rows.append({"case": case, "variable": variable, "t": t, "value": round(float(value), 6)})


# ═══════════════════════════════════════════════════════════════════════════════
# Case 1 — LOPF uncongested
# ═══════════════════════════════════════════════════════════════════════════════
print("Case 1: LOPF uncongested ...", end=" ")
net = base_3bus(s_nom=1e6)
net.add("Load", "D2", bus="Bus2", p_set=200.0)
net.add("Load", "D3", bus="Bus3", p_set=300.0)
net.optimize(solver_name="highs")

for g in ["G1", "G2"]:
    save("LOPF_uncon", f"P_{g}", 0, net.generators_t.p[g].iloc[0])
for b in ["Bus1", "Bus2", "Bus3"]:
    save("LOPF_uncon", f"LMP_{b}", 0, net.buses_t.marginal_price[b].iloc[0])
save("LOPF_uncon", "cost", 0, net.objective)
print("done")


# ═══════════════════════════════════════════════════════════════════════════════
# Case 2 — LOPF congested  (s_nom = 200 MW)
# ═══════════════════════════════════════════════════════════════════════════════
print("Case 2: LOPF congested (s_nom=200) ...", end=" ")
net = base_3bus(s_nom=200.0)
net.add("Load", "D2", bus="Bus2", p_set=200.0)
net.add("Load", "D3", bus="Bus3", p_set=300.0)
net.optimize(solver_name="highs")

for g in ["G1", "G2"]:
    save("LOPF_con", f"P_{g}", 0, net.generators_t.p[g].iloc[0])
for b in ["Bus1", "Bus2", "Bus3"]:
    save("LOPF_con", f"LMP_{b}", 0, net.buses_t.marginal_price[b].iloc[0])
save("LOPF_con", "cost", 0, net.objective)
print("done")


# ═══════════════════════════════════════════════════════════════════════════════
# Case 3 — Multi-period LOPF  T=6, StorageUnit, no wind
# ═══════════════════════════════════════════════════════════════════════════════
print("Case 3: Multi-period LOPF T=6 ...", end=" ")
net = base_3bus(s_nom=1e6, snapshots=range(6))

p_load2 = 200.0 * LOAD_PROFILE_6
p_load3 = 300.0 * LOAD_PROFILE_6
net.add("Load", "D2", bus="Bus2", p_set=p_load2)
net.add("Load", "D3", bus="Bus3", p_set=p_load3)
net.add("StorageUnit", "Bat",
        bus="Bus2", p_nom=60.0,
        max_hours=4.0,               # e_nom = 60 * 4 = 240 MWh
        efficiency_store=0.9,
        efficiency_dispatch=0.9,
        cyclic_state_of_charge=False,   # Julia fixes e_initial; match that
        state_of_charge_initial=120.0)
net.optimize(solver_name="highs")

for t in range(6):
    for g in ["G1", "G2"]:
        save("MP_t6", f"P_{g}", t, net.generators_t.p[g].iloc[t])
    save("MP_t6", "SoC_Bat", t, net.storage_units_t.state_of_charge["Bat"].iloc[t])
    for b in ["Bus1", "Bus2", "Bus3"]:
        save("MP_t6", f"LMP_{b}", t, net.buses_t.marginal_price[b].iloc[t])
save("MP_t6", "cost", 0, net.objective)
print("done")


# ═══════════════════════════════════════════════════════════════════════════════
# Case 4 — Unit Commitment  T=6
# ═══════════════════════════════════════════════════════════════════════════════
print("Case 4: Unit Commitment T=6 ...", end=" ")
net = base_3bus(s_nom=1e6, snapshots=range(6))

p_load2 = 200.0 * LOAD_PROFILE_6
p_load3 = 300.0 * LOAD_PROFILE_6
net.add("Load", "D2", bus="Bus2", p_set=p_load2)
net.add("Load", "D3", bus="Bus3", p_set=p_load3)
# Replace G2 with a committable peaker
net.remove("Generator", "G2")
net.add("Generator", "G_peak",
        bus="Bus2", p_nom=300.0, marginal_cost=50.0,
        committable=True,
        min_up_time=2, min_down_time=1,
        start_up_cost=1000.0, shut_down_cost=0.0,
        p_min_pu=0.3, p_max_pu=1.0,
        initial_status=0)   # match Julia default (initial_status=false)
net.optimize(solver_name="highs")

for t in range(6):
    save("UC_t6", "P_G1",     t, net.generators_t.p["G1"].iloc[t])
    save("UC_t6", "P_G_peak", t, net.generators_t.p["G_peak"].iloc[t])
    save("UC_t6", "u_G_peak", t,
         net.generators_t["status"]["G_peak"].iloc[t]
         if "status" in net.generators_t else 1)  # fallback if not extracted
save("UC_t6", "cost", 0, net.objective)
print("done")


# ═══════════════════════════════════════════════════════════════════════════════
# Save CSV
# ═══════════════════════════════════════════════════════════════════════════════
out = "../results/pypsa_parity.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["case", "variable", "t", "value"])
    w.writeheader()
    w.writerows(rows)

print(f"\n[OK] Saved {len(rows)} rows to {out}")
print("Now run: julia julia/validation/validate_parity.jl")
