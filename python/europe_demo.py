"""
Europe 10-bus demo — PyPSA reference with interactive geographic map.

10 major European cities as buses, realistic interconnections,
base + peak generators, demand proportional to country size.

Outputs:
  results/europe_map.html      — interactive Plotly map (open in browser)
  results/europe_lopf.csv      — LOPF results for Julia comparison
"""

import warnings, logging
warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)

import numpy as np
import pypsa
import pandas as pd
import csv, time

# ── Bus definitions (city, lat, lon, peak_load_MW) ────────────────────────────
BUSES = {
    "London":    {"y": 51.51, "x": -0.13, "load": 5600},
    "Paris":     {"y": 48.85, "x":  2.35, "load": 4750},
    "Berlin":    {"y": 52.52, "x": 13.40, "load": 4000},
    "Madrid":    {"y": 40.42, "x": -3.70, "load": 3500},
    "Rome":      {"y": 41.90, "x": 12.50, "load": 3250},
    "Amsterdam": {"y": 52.37, "x":  4.90, "load": 2250},
    "Warsaw":    {"y": 52.23, "x": 21.01, "load": 2500},
    "Vienna":    {"y": 48.21, "x": 16.37, "load": 1875},
    "Stockholm": {"y": 59.33, "x": 18.07, "load": 2125},
    "Zurich":    {"y": 47.38, "x":  8.54, "load": 1500},
}

# ── Interconnections (from, to, s_nom_MW, x_pu) ──────────────────────────────
LINES = [
    ("London",    "Amsterdam", 2500, 0.08),
    ("London",    "Paris",     3000, 0.06),
    ("Paris",     "Berlin",    2000, 0.10),
    ("Paris",     "Madrid",    1500, 0.14),
    ("Paris",     "Zurich",    2000, 0.07),
    ("Berlin",    "Amsterdam", 1800, 0.07),
    ("Berlin",    "Warsaw",    1500, 0.09),
    ("Berlin",    "Vienna",    2000, 0.09),
    ("Berlin",    "Stockholm", 1200, 0.13),
    ("Zurich",    "Rome",      1800, 0.09),
    ("Zurich",    "Vienna",    2000, 0.06),
    ("Vienna",    "Warsaw",    1200, 0.10),
    ("Vienna",    "Rome",      1500, 0.10),
    ("Amsterdam", "London",    2500, 0.08),   # duplicate direction intentional
]
# Deduplicate by sorted tuple
seen = set()
LINES_CLEAN = []
for f, t, s, x in LINES:
    k = tuple(sorted([f, t]))
    if k not in seen:
        LINES_CLEAN.append((f, t, s, x))
        seen.add(k)

# ── Generator mix per bus (marginal_cost €/MWh, p_nom_MW) ────────────────────
GENERATORS = {
    "London":    [("nuclear",  20.0, 12000), ("gas",    65.0, 8000)],
    "Paris":     [("nuclear",  18.0, 20000), ("hydro",  10.0, 5000)],
    "Berlin":    [("coal",     40.0, 10000), ("wind",    5.0, 15000), ("gas", 70.0, 5000)],
    "Madrid":    [("gas",      55.0,  8000), ("solar",   3.0, 12000)],
    "Rome":      [("gas",      60.0,  7000), ("hydro",  12.0,  4000)],
    "Amsterdam": [("gas",      55.0,  5000), ("wind",    5.0,  6000)],
    "Warsaw":    [("coal",     45.0,  8000), ("gas",    75.0,  3000)],
    "Vienna":    [("hydro",    10.0,  4000), ("gas",    65.0,  3000)],
    "Stockholm": [("hydro",     8.0,  8000), ("nuclear",20.0,  5000)],
    "Zurich":    [("hydro",     9.0,  3000), ("nuclear",22.0,  2000)],
}

# ── 24-hour load profile (normalised) ────────────────────────────────────────
LOAD_PROFILE = np.array([
    0.60, 0.57, 0.55, 0.54, 0.55, 0.60,
    0.70, 0.80, 0.88, 0.90, 0.92, 0.91,
    0.90, 0.89, 0.88, 0.87, 0.89, 0.95,
    1.00, 0.98, 0.93, 0.85, 0.75, 0.65,
])

# ── Build PyPSA network ───────────────────────────────────────────────────────
T = 24
net = pypsa.Network()
net.set_snapshots(range(T))

for city, props in BUSES.items():
    net.add("Bus", city, x=props["x"], y=props["y"], v_nom=380.0)

for f, t, s_nom, x in LINES_CLEAN:
    net.add("Line", f"{f[:3]}-{t[:3]}",
            bus0=f, bus1=t, x=x, r=0.01, s_nom=s_nom)

for city, gens in GENERATORS.items():
    for i, (fuel, mc, pnom) in enumerate(gens):
        ctrl = "Slack" if (city == "Paris" and fuel == "nuclear" and i == 0) else "PQ"
        net.add("Generator", f"{city[:3]}_{fuel}",
                bus=city, p_nom=pnom, marginal_cost=mc, control=ctrl)

for city, props in BUSES.items():
    p_set = props["load"] * LOAD_PROFILE
    net.add("Load", f"D_{city[:3]}", bus=city, p_set=p_set)

# ── Run LOPF ──────────────────────────────────────────────────────────────────
print("Running PyPSA LOPF (T=24) ...")
t0 = time.perf_counter()
net.optimize(solver_name="highs")
elapsed = (time.perf_counter() - t0) * 1000
print(f"  Solved in {elapsed:.1f} ms  |  Total cost: {net.objective/1e6:.2f} M€")

# ── LMP summary ───────────────────────────────────────────────────────────────
lmp = net.buses_t.marginal_price
print("\nLMP summary (€/MWh):")
print(f"  {'Bus':<12}  {'Mean':>8}  {'Min':>8}  {'Max':>8}")
for city in BUSES:
    v = lmp[city]
    print(f"  {city:<12}  {v.mean():>8.2f}  {v.min():>8.2f}  {v.max():>8.2f}")

# ── Save results CSV for Julia comparison ────────────────────────────────────
rows = []
for city in BUSES:
    for t in range(T):
        for gen in net.generators.index:
            if net.generators.loc[gen, "bus"] == city:
                rows.append({"bus": city, "generator": gen, "t": t,
                              "p_mw": net.generators_t.p[gen].iloc[t],
                              "lmp":  lmp[city].iloc[t]})

with open("../results/europe_lopf.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["bus", "generator", "t", "p_mw", "lmp"])
    w.writeheader()
    w.writerows(rows)
print(f"\n[OK] Saved {len(rows)} rows to results/europe_lopf.csv")

# ── Interactive map ───────────────────────────────────────────────────────────
print("\nGenerating interactive map ...")
try:
    fig = net.iplot(
        bus_sizes=net.loads_t.p.mean() / 3e4,
        line_widths=0.5,
        title="European 10-Bus Network — LOPF Results",
        return_fig=True,
    )
    fig.write_html("../results/europe_map.html")
    print("[OK] Map saved to results/europe_map.html  (open in browser)")
except Exception as e:
    # Fallback: plain plotly scatter map
    import plotly.graph_objects as go

    fig = go.Figure()

    # Lines
    for f, t, *_ in LINES_CLEAN:
        fig.add_trace(go.Scattergeo(
            lon=[BUSES[f]["x"], BUSES[t]["x"]],
            lat=[BUSES[f]["y"], BUSES[t]["y"]],
            mode="lines",
            line=dict(width=1.5, color="#1565C0"),
            showlegend=False,
        ))

    # Bus LMP (mean over 24h) as colour
    lmp_mean = [lmp[c].mean() for c in BUSES]
    fig.add_trace(go.Scattergeo(
        lon=[BUSES[c]["x"] for c in BUSES],
        lat=[BUSES[c]["y"] for c in BUSES],
        text=[f"{c}<br>LMP: {lmp[c].mean():.1f} €/MWh" for c in BUSES],
        mode="markers+text",
        textposition="top center",
        marker=dict(
            size=12,
            color=lmp_mean,
            colorscale="RdYlGn_r",
            colorbar=dict(title="LMP (€/MWh)"),
            showscale=True,
        ),
        showlegend=False,
    ))

    fig.update_layout(
        title="European 10-Bus Network — Mean LMP (€/MWh)",
        geo=dict(
            scope="europe",
            projection_type="natural earth",
            showland=True, landcolor="#F5F5F5",
            showocean=True, oceancolor="#E3F2FD",
            showcountries=True, countrycolor="#BDBDBD",
            center=dict(lat=50, lon=10),
        ),
        height=600,
    )
    fig.write_html("../results/europe_map.html")
    print("[OK] Map saved to results/europe_map.html  (open in browser)")
