"""
Multi-period LOPF with StorageUnit and Wind — PyPSA reference.
Same network, profiles, and parameters as julia/solvers/lopf_multiperiod.jl.
"""
import numpy as np
import pypsa
import warnings
import logging

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

# ----------------------------------------------------------------
#  24-hour profiles (identical to lopf_multiperiod.jl)
# ----------------------------------------------------------------
LOAD_PROFILE = np.array([
    0.60, 0.57, 0.55, 0.54, 0.55, 0.60,
    0.70, 0.80, 0.88, 0.90, 0.92, 0.91,
    0.90, 0.89, 0.88, 0.87, 0.89, 0.95,
    1.00, 0.98, 0.93, 0.85, 0.75, 0.65
])

WIND_PROFILE = np.array([
    0.80, 0.82, 0.85, 0.83, 0.78, 0.70,
    0.60, 0.55, 0.50, 0.45, 0.42, 0.40,
    0.38, 0.37, 0.40, 0.43, 0.50, 0.58,
    0.65, 0.70, 0.74, 0.76, 0.78, 0.80
])

# ----------------------------------------------------------------
#  Build PyPSA network
# ----------------------------------------------------------------
net = pypsa.Network()
hours = list(range(24))
net.set_snapshots(hours)

# Buses
for i in [1, 2, 3]:
    net.add("Bus", f"Bus{i}", v_nom=380.0)

# Lines
net.add("Line", "L12", bus0="Bus1", bus1="Bus2", x=0.1, r=0.01, s_nom=1e6)
net.add("Line", "L13", bus0="Bus1", bus1="Bus3", x=0.1, r=0.01, s_nom=1e6)
net.add("Line", "L23", bus0="Bus2", bus1="Bus3", x=0.1, r=0.01, s_nom=1e6)

# G1 limited to 270 MW so peak demand forces G2 / storage activation
net.add("Generator", "G1",
        bus="Bus1", p_nom=270.0, marginal_cost=20.0, control="Slack")
net.add("Generator", "G2",
        bus="Bus2", p_nom=100.0, marginal_cost=50.0)

# Loads (time-varying); total peak = 425 MW
p_base_bus2 = 250.0
p_base_bus3 = 175.0
net.add("Load", "Load2", bus="Bus2",
        p_set=p_base_bus2 * LOAD_PROFILE)   # pandas Series over snapshots
net.add("Load", "Load3", bus="Bus3",
        p_set=p_base_bus3 * LOAD_PROFILE)

# Wind generator — zero marginal cost, variable p_max_pu
net.add("Generator", "Wind3",
        bus="Bus3",
        p_nom=150.0,
        marginal_cost=0.0,
        p_max_pu=WIND_PROFILE)

# Storage unit at Bus 2
# PyPSA StorageUnit: p_nom [MW], max_hours [h] → e_nom = p_nom * max_hours
net.add("StorageUnit", "Stor2",
        bus="Bus2",
        p_nom=100.0,
        max_hours=3.0,          # E_nom = 100 * 3 = 300 MWh
        efficiency_store=0.95,
        efficiency_dispatch=0.95,
        cyclic_state_of_charge=True,
        state_of_charge_initial=150.0)  # 50% of 300 MWh

# ----------------------------------------------------------------
#  Solve
# ----------------------------------------------------------------
print("=" * 65)
print("MULTI-PERIOD LOPF  (24h, Storage + Wind)  — PyPSA")
print("=" * 65)

net.optimize(solver_name="highs")

# ----------------------------------------------------------------
#  Results
# ----------------------------------------------------------------
total_cost = (
    net.generators_t.p["G1"].values * 20.0 +
    net.generators_t.p["G2"].values * 50.0
).sum()

total_load = (net.loads_t.p_set["Load2"] + net.loads_t.p_set["Load3"]).sum()
avg_cost   = total_cost / total_load

print(f"\nTotal generation cost (24h): {total_cost:.2f} €")
print(f"Average cost per MWh:        {avg_cost:.2f} €/MWh")

print(f"\n{'Hour':>4} | {'Load (MW)':>9} | {'Wind (MW)':>9} | "
      f"{'G1 (MW)':>7} | {'G2 (MW)':>7} | {'Stor (MW)':>9} | {'SOC (MWh)':>9}")
print("-" * 80)

soc = net.storage_units_t.state_of_charge["Stor2"].values
p_stor = net.storage_units_t.p["Stor2"].values   # positive = discharging

for t in range(24):
    load_t = (net.loads_t.p_set["Load2"].iloc[t] +
              net.loads_t.p_set["Load3"].iloc[t])
    wind_t = net.generators_t.p["Wind3"].iloc[t]
    g1_t   = net.generators_t.p["G1"].iloc[t]
    g2_t   = net.generators_t.p["G2"].iloc[t]
    print(f"{t:>4} | {load_t:>9.1f} | {wind_t:>9.1f} | "
          f"{g1_t:>7.1f} | {g2_t:>7.1f} | {p_stor[t]:>9.1f} | {soc[t]:>9.1f}")

print("\n" + "=" * 65)
print("Done. Compare with julia/solvers/lopf_multiperiod.jl")
print("=" * 65)
