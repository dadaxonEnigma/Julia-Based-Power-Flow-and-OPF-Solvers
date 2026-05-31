"""
validation_full.py — PyPSA reference values for all components and solvers.

Saves: results/pypsa_validation_full.csv
Compare with: julia julia/validation/validation_full.jl
Then diff: python python/compare_validation.py

Tests covered
─────────────
  01_dc_pf        DC Power Flow (3-bus, lpf)
  02_lopf_free    LOPF uncongested (3-bus)
  03_lopf_cong    LOPF congested s_nom=200 MW
  04_lopf_trafo   LOPF with Transformer (tap=1)
  05_lopf_link    LOPF with HVDC Link
  06_lopf_co2     LOPF with CO2 GlobalConstraint
  07_mp_storage   Multi-period LOPF T=6 + StorageUnit
  08_uc           Unit Commitment T=6
"""

import csv, warnings, logging
import numpy as np
import pypsa

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

rows = []

def save(test_id, variable, t, value):
    rows.append({
        "test_id": test_id,
        "variable": variable,
        "t": int(t),
        "value": round(float(value), 8),
    })

LP6 = np.array([0.60, 0.57, 0.55, 0.54, 0.55, 0.60])

# ─────────────────────────────────────────────────────────────────────────────
#  Helpers
# ─────────────────────────────────────────────────────────────────────────────
def base_3bus(s_nom=1e9, snapshots=None, x_l23=0.15):
    net = pypsa.Network()
    if snapshots is not None:
        net.set_snapshots(snapshots)
    else:
        net.set_snapshots([0])
    net.add("Bus", "Bus1", v_nom=380.0)
    net.add("Bus", "Bus2", v_nom=380.0)
    net.add("Bus", "Bus3", v_nom=380.0)
    net.add("Line", "L12", bus0="Bus1", bus1="Bus2", r=0.01, x=0.1,   s_nom=s_nom)
    net.add("Line", "L13", bus0="Bus1", bus1="Bus3", r=0.01, x=0.1,   s_nom=s_nom)
    net.add("Line", "L23", bus0="Bus2", bus1="Bus3", r=0.01, x=x_l23, s_nom=s_nom)
    return net

def add_gens_loads(net, p_g1=400, p_g2=300, p_d2=200, p_d3=300):
    net.add("Generator", "G1", bus="Bus1", p_nom=p_g1,
            marginal_cost=20.0, control="Slack")
    net.add("Generator", "G2", bus="Bus2", p_nom=p_g2, marginal_cost=50.0)
    net.add("Load",      "D2", bus="Bus2", p_set=p_d2)
    net.add("Load",      "D3", bus="Bus3", p_set=p_d3)

# ─────────────────────────────────────────────────────────────────────────────
#  01  DC Power Flow
# ─────────────────────────────────────────────────────────────────────────────
print("01 DC Power Flow ...", end=" ", flush=True)
net = base_3bus()
# Only G1 at slack (p_nom=600 so it covers all load=500 MW)
net.add("Generator", "G1", bus="Bus1", p_nom=600.0,
        p_max_pu=1.0, marginal_cost=20.0, control="Slack")
net.add("Load", "D2", bus="Bus2", p_set=200.0)
net.add("Load", "D3", bus="Bus3", p_set=300.0)
net.lpf()

for b in ["Bus1", "Bus2", "Bus3"]:
    save("01_dc_pf", f"theta_{b}", 0, net.buses_t.v_ang[b].iloc[0])
for l in ["L12", "L13", "L23"]:
    save("01_dc_pf", f"flow_{l}", 0, net.lines_t.p0[l].iloc[0])
print("done")

# ─────────────────────────────────────────────────────────────────────────────
#  02  LOPF uncongested
# ─────────────────────────────────────────────────────────────────────────────
print("02 LOPF uncongested ...", end=" ", flush=True)
net = base_3bus()
add_gens_loads(net)
net.optimize(solver_name="highs")

save("02_lopf_free", "P_G1",     0, net.generators_t.p["G1"].iloc[0])
save("02_lopf_free", "P_G2",     0, net.generators_t.p["G2"].iloc[0])
save("02_lopf_free", "total_cost", 0, net.objective)
for b in ["Bus1", "Bus2", "Bus3"]:
    save("02_lopf_free", f"LMP_{b}", 0, net.buses_t.marginal_price[b].iloc[0])
print("done")

# ─────────────────────────────────────────────────────────────────────────────
#  03  LOPF congested (s_nom = 200 MW)
# ─────────────────────────────────────────────────────────────────────────────
print("03 LOPF congested ...", end=" ", flush=True)
net = base_3bus(s_nom=200.0)
add_gens_loads(net)
net.optimize(solver_name="highs")

save("03_lopf_cong", "P_G1",     0, net.generators_t.p["G1"].iloc[0])
save("03_lopf_cong", "P_G2",     0, net.generators_t.p["G2"].iloc[0])
save("03_lopf_cong", "total_cost", 0, net.objective)
for b in ["Bus1", "Bus2", "Bus3"]:
    save("03_lopf_cong", f"LMP_{b}", 0, net.buses_t.marginal_price[b].iloc[0])
print("done")

# ─────────────────────────────────────────────────────────────────────────────
#  04  LOPF with Transformer (tap_ratio=1)
# ─────────────────────────────────────────────────────────────────────────────
print("04 LOPF with Transformer ...", end=" ", flush=True)
net = pypsa.Network()
net.set_snapshots([0])
net.add("Bus", "Bus1", v_nom=380.0)
net.add("Bus", "Bus2", v_nom=380.0)
net.add("Bus", "Bus3", v_nom=380.0)
# Replace L12 with a transformer
net.add("Transformer", "T12",
        bus0="Bus1", bus1="Bus2",
        x=0.1, s_nom=500.0,
        tap_ratio=1.0, phase_shift=0.0)
net.add("Line", "L13", bus0="Bus1", bus1="Bus3", r=0.01, x=0.1,  s_nom=1e9)
net.add("Line", "L23", bus0="Bus2", bus1="Bus3", r=0.01, x=0.15, s_nom=1e9)
net.add("Generator", "G1", bus="Bus1", p_nom=400.0, marginal_cost=20.0, control="Slack")
net.add("Generator", "G2", bus="Bus2", p_nom=300.0, marginal_cost=50.0)
net.add("Load",      "D2", bus="Bus2", p_set=200.0)
net.add("Load",      "D3", bus="Bus3", p_set=300.0)
net.optimize(solver_name="highs")

save("04_lopf_trafo", "P_G1",     0, net.generators_t.p["G1"].iloc[0])
save("04_lopf_trafo", "P_G2",     0, net.generators_t.p["G2"].iloc[0])
save("04_lopf_trafo", "total_cost", 0, net.objective)
save("04_lopf_trafo", "flow_T12", 0, net.transformers_t.p0["T12"].iloc[0])
for b in ["Bus1", "Bus2", "Bus3"]:
    save("04_lopf_trafo", f"LMP_{b}", 0, net.buses_t.marginal_price[b].iloc[0])
print("done")

# ─────────────────────────────────────────────────────────────────────────────
#  05  LOPF with HVDC Link
# ─────────────────────────────────────────────────────────────────────────────
print("05 LOPF with HVDC Link ...", end=" ", flush=True)
net = pypsa.Network()
net.set_snapshots([0])
net.add("Bus", "Bus1", v_nom=380.0)
net.add("Bus", "Bus2", v_nom=380.0)
# AC line Bus1-Bus2
net.add("Line", "L12", bus0="Bus1", bus1="Bus2", r=0.01, x=0.2, s_nom=150.0)
# HVDC Link Bus1-Bus2 (parallel to AC line, efficiency=0.97)
net.add("Link", "HVDC",
        bus0="Bus1", bus1="Bus2",
        p_nom=200.0, efficiency=0.97, marginal_cost=0.0)
net.add("Generator", "G1", bus="Bus1", p_nom=400.0, marginal_cost=20.0, control="Slack")
net.add("Generator", "G2", bus="Bus2", p_nom=200.0, marginal_cost=80.0)
net.add("Load", "D2", bus="Bus2", p_set=350.0)
net.optimize(solver_name="highs")

save("05_lopf_link", "P_G1",     0, net.generators_t.p["G1"].iloc[0])
save("05_lopf_link", "P_G2",     0, net.generators_t.p["G2"].iloc[0])
save("05_lopf_link", "p_HVDC",   0, net.links_t.p0["HVDC"].iloc[0])
save("05_lopf_link", "total_cost", 0, net.objective)
for b in ["Bus1", "Bus2"]:
    save("05_lopf_link", f"LMP_{b}", 0, net.buses_t.marginal_price[b].iloc[0])
print("done")

# ─────────────────────────────────────────────────────────────────────────────
#  06  LOPF with CO2 GlobalConstraint
# ─────────────────────────────────────────────────────────────────────────────
print("06 LOPF with CO2 cap ...", end=" ", flush=True)
net = pypsa.Network()
net.set_snapshots([0])
net.add("Bus", "B1", v_nom=380.0)
net.add("Carrier", "coal", co2_emissions=0.34)
net.add("Carrier", "gas",  co2_emissions=0.20)
# Load=300 MW; Gas p_nom=300 → Coal can be 0 (all-gas CO2=60≤80).
# Cap=80 forces: Coal*0.14 ≤ 20 → Coal≤142.86 MW → partial gas dispatch required.
net.add("Generator", "Coal", bus="B1", p_nom=300.0, marginal_cost=20.0, carrier="coal")
net.add("Generator", "Gas",  bus="B1", p_nom=300.0, marginal_cost=50.0, carrier="gas")
net.add("Load", "D1", bus="B1", p_set=300.0)
net.add("GlobalConstraint", "co2_cap",
        carrier_attribute="co2_emissions",
        sense="<=",
        constant=80.0)          # all-coal=102t > 80, all-gas=60t ≤ 80 → mixed dispatch
net.optimize(solver_name="highs")

save("06_lopf_co2", "P_Coal",    0, net.generators_t.p["Coal"].iloc[0])
save("06_lopf_co2", "P_Gas",     0, net.generators_t.p["Gas"].iloc[0])
save("06_lopf_co2", "total_cost", 0, net.objective)
save("06_lopf_co2", "LMP_B1",   0, net.buses_t.marginal_price["B1"].iloc[0])
print("done")

# ─────────────────────────────────────────────────────────────────────────────
#  07  Multi-period LOPF T=6, StorageUnit
# ─────────────────────────────────────────────────────────────────────────────
print("07 Multi-period LOPF T=6 + StorageUnit ...", end=" ", flush=True)
net = base_3bus(snapshots=range(6))
net.add("Generator", "G1", bus="Bus1", p_nom=400.0, marginal_cost=20.0, control="Slack")
net.add("Generator", "G2", bus="Bus2", p_nom=300.0, marginal_cost=50.0)
net.add("Load", "D2", bus="Bus2", p_set=200.0 * LP6)
net.add("Load", "D3", bus="Bus3", p_set=300.0 * LP6)
net.add("StorageUnit", "Bat",
        bus="Bus2", p_nom=60.0,
        max_hours=4.0,                    # e_nom = 60*4 = 240 MWh
        efficiency_store=0.9,
        efficiency_dispatch=0.9,
        cyclic_state_of_charge=False,
        state_of_charge_initial=120.0)
net.optimize(solver_name="highs")

save("07_mp_storage", "total_cost", 0, net.objective)
for t in range(6):
    save("07_mp_storage", "P_G1",    t, net.generators_t.p["G1"].iloc[t])
    save("07_mp_storage", "P_G2",    t, net.generators_t.p["G2"].iloc[t])
    save("07_mp_storage", "SoC_Bat", t, net.storage_units_t.state_of_charge["Bat"].iloc[t])
    save("07_mp_storage", "p_ch",    t, net.storage_units_t.p_store["Bat"].iloc[t]
                                        if "p_store" in net.storage_units_t else 0.0)
    for b in ["Bus1", "Bus2", "Bus3"]:
        save("07_mp_storage", f"LMP_{b}", t, net.buses_t.marginal_price[b].iloc[t])
print("done")

# ─────────────────────────────────────────────────────────────────────────────
#  08  Unit Commitment T=6
# ─────────────────────────────────────────────────────────────────────────────
print("08 Unit Commitment T=6 ...", end=" ", flush=True)
net = base_3bus(snapshots=range(6))
net.add("Generator", "G1", bus="Bus1", p_nom=400.0, marginal_cost=20.0, control="Slack")
net.add("Generator", "G_peak",
        bus="Bus2", p_nom=300.0, marginal_cost=50.0,
        committable=True,
        min_up_time=2, min_down_time=1,
        start_up_cost=1000.0, shut_down_cost=0.0,
        p_min_pu=0.3, p_max_pu=1.0,
        initial_status=0)
net.add("Load", "D2", bus="Bus2", p_set=200.0 * LP6)
net.add("Load", "D3", bus="Bus3", p_set=300.0 * LP6)
net.optimize(solver_name="highs")

save("08_uc", "total_cost", 0, net.objective)
for t in range(6):
    save("08_uc", "P_G1",     t, net.generators_t.p["G1"].iloc[t])
    save("08_uc", "P_G_peak", t, net.generators_t.p["G_peak"].iloc[t])
    try:
        u_val = net.generators_t["status"]["G_peak"].iloc[t]
    except Exception:
        u_val = 0.0
    save("08_uc", "u_G_peak", t, u_val)
    for b in ["Bus1", "Bus2", "Bus3"]:
        save("08_uc", f"LMP_{b}", t, net.buses_t.marginal_price[b].iloc[t])
print("done")

# ─────────────────────────────────────────────────────────────────────────────
#  Save CSV
# ─────────────────────────────────────────────────────────────────────────────
out = "../results/pypsa_validation_full.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["test_id", "variable", "t", "value"])
    w.writeheader()
    w.writerows(rows)

print(f"\n[OK] {len(rows)} rows → {out}")
print("Next: julia julia/validation/validation_full.jl")
print("Then: python python/compare_validation.py")
