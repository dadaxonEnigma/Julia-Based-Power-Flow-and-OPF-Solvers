"""Fresh re-run of the 08_uc validation case to settle the Julia/PyPSA cost
discrepancy. Tests with and without control='Slack' on G1."""
import warnings, logging
import numpy as np
import pypsa
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)

LP6 = np.array([0.60, 0.57, 0.55, 0.54, 0.55, 0.60])

def build(slack):
    net = pypsa.Network(); net.set_snapshots(range(6))
    for b in ("Bus1", "Bus2", "Bus3"):
        net.add("Bus", b, v_nom=380.0)
    net.add("Line", "L12", bus0="Bus1", bus1="Bus2", r=0.01, x=0.10, s_nom=1e9)
    net.add("Line", "L13", bus0="Bus1", bus1="Bus3", r=0.01, x=0.10, s_nom=1e9)
    net.add("Line", "L23", bus0="Bus2", bus1="Bus3", r=0.01, x=0.15, s_nom=1e9)
    if slack:
        net.add("Generator", "G1", bus="Bus1", p_nom=400.0, marginal_cost=20.0, control="Slack")
    else:
        net.add("Generator", "G1", bus="Bus1", p_nom=400.0, marginal_cost=20.0)
    net.add("Generator", "G_peak", bus="Bus2", p_nom=300.0, marginal_cost=50.0,
            committable=True, min_up_time=2, min_down_time=1,
            start_up_cost=1000.0, shut_down_cost=0.0, p_min_pu=0.3, p_max_pu=1.0,
            initial_status=0)
    net.add("Load", "D2", bus="Bus2", p_set=200.0 * LP6)
    net.add("Load", "D3", bus="Bus3", p_set=300.0 * LP6)
    return net

for slack in (True, False):
    net = build(slack)
    net.optimize(solver_name="highs")
    gp = net.generators_t.p["G_peak"].values
    g1 = net.generators_t.p["G1"].values
    print(f"control={'Slack' if slack else 'None'}:  objective={net.objective:.1f}  "
          f"G1[0]={g1[0]:.1f}  G_peak={np.round(gp,1).tolist()}")
print("Demand per hour:", np.round(500 * LP6, 1).tolist(), " (G1 alone, max 400, suffices)")
print("Expected optimum = 1705 MWh * 20 =", 1705 * 20)

print("\n=== diagnostic variants (control=None) ===")
def variant(**ch):
    net = pypsa.Network(); net.set_snapshots(range(6))
    for b in ("Bus1","Bus2","Bus3"): net.add("Bus", b, v_nom=380.0)
    net.add("Line","L12",bus0="Bus1",bus1="Bus2",r=0.01,x=0.10,s_nom=1e9)
    net.add("Line","L13",bus0="Bus1",bus1="Bus3",r=0.01,x=0.10,s_nom=1e9)
    net.add("Line","L23",bus0="Bus2",bus1="Bus3",r=0.01,x=0.15,s_nom=1e9)
    net.add("Generator","G1",bus="Bus1",p_nom=400.0,marginal_cost=20.0)
    gp=dict(p_nom=300.0,marginal_cost=50.0,committable=True,min_up_time=2,
            min_down_time=1,start_up_cost=1000.0,shut_down_cost=0.0,
            p_min_pu=0.3,p_max_pu=1.0,initial_status=0)
    gp.update(ch)
    net.add("Generator","G_peak",bus="Bus2",**gp)
    net.add("Load","D2",bus="Bus2",p_set=200.0*LP6)
    net.add("Load","D3",bus="Bus3",p_set=300.0*LP6)
    net.optimize(solver_name="highs")
    return net.objective, net.generators_t.p["G_peak"].round(1).tolist()

for label,ch in [("min_up_time=1",dict(min_up_time=1)),
                 ("p_min_pu=0.0",dict(p_min_pu=0.0)),
                 ("initial_status=1",dict(initial_status=1)),
                 ("min_up_time=1,p_min_pu=0",dict(min_up_time=1,p_min_pu=0.0))]:
    obj,gp=variant(**ch)
    print(f"  {label:28s} objective={obj:.1f}  G_peak={gp}")
