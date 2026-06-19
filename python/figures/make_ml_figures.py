"""
make_ml_figures.py — the four AI figures, from the honest real-data exports
produced by julia/examples/ml_figures_data.jl.

Figures (results/plots/):
  fig_ml_accuracy.png        MAPE bar chart: LSTM (+/- temperature) vs naive baselines
  fig_ml_forecast.png        sample-day forecast with 90% conformal band vs actual
  fig_ml_coverage.png        per-hour empirical coverage vs nominal 90%
  fig_ml_scenario_costs.png  stochastic-LOPF cost distribution with E[cost] and CVaR90

Run (after the Julia export):
  pypsa_env/Scripts/python.exe python/make_ml_figures.py
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif", "font.size": 12, "axes.labelsize": 12,
    "legend.fontsize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10,
    "axes.grid": True, "grid.alpha": 0.3, "grid.linestyle": ":",
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
})
RES = os.path.join(os.path.dirname(__file__), "..", "..", "results")
PLT = os.path.join(RES, "plots")
os.makedirs(PLT, exist_ok=True)


def save(fig, name):
    fig.savefig(os.path.join(PLT, name)); plt.close(fig)
    print("  saved plots/" + name)


# ── 1. Forecast accuracy (verified multi-seed MAPE, real out-of-time test) ───
# Source: real_data_forecast.jl, 3 seeds, 60-day out-of-time test. Read straight
# from the committed CSVs so this figure can never drift from tab:forecast_metrics.
_temp   = pd.read_csv(os.path.join(RES, "ml", "real_data_forecast_temp.csv")).set_index("model")
_notemp = pd.read_csv(os.path.join(RES, "ml", "real_data_forecast_notemp.csv")).set_index("model")
methods = ["Seasonal-naive", "LSTM\n(+temperature)", "LSTM\n(calendar)", "Persistence"]
mape    = [float(_temp.loc["Seasonal-naive", "mape_mean"]),
           float(_temp.loc["LSTM", "mape_mean"]),
           float(_notemp.loc["LSTM", "mape_mean"]),
           float(_temp.loc["Persistence", "mape_mean"])]
colors  = ["tab:green", "tab:blue", "tab:cyan", "tab:gray"]
fig, ax = plt.subplots(figsize=(6.5, 4))
bars = ax.bar(methods, mape, color=colors, edgecolor="black", linewidth=0.6)
for b, v in zip(bars, mape):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.1, f"{v:.2f}", ha="center", fontsize=10)
ax.set_ylabel("MAPE (%)  — lower is better")
ax.set_ylim(0, max(mape) * 1.18)
ax.axhline(mape[0], color="tab:green", linestyle="--", linewidth=1, alpha=0.7)
save(fig, "fig_ml_accuracy.png")

# ── 2. Sample-day forecast with 90% conformal band ──────────────────────────
sd = pd.read_csv(os.path.join(RES, "ml", "ml_sample_day.csv"))
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.fill_between(sd.hour, sd.lower, sd.upper, color="tab:blue", alpha=0.20,
                label="90% conformal interval")
ax.plot(sd.hour, sd["mean"], color="tab:blue", marker="o", markersize=4, label="LSTM forecast")
ax.plot(sd.hour, sd.actual, color="black", linestyle="--", marker="s", markersize=4, label="Actual load")
ax.set_xlabel("Hour of day")
ax.set_ylabel("Load (per unit)")
ax.set_xlim(0, 23)
ax.legend(frameon=True, loc="best")
save(fig, "fig_ml_forecast.png")

# ── 3. Stochastic LOPF cost distribution with E[cost] and CVaR90 ─────────────
sc = pd.read_csv(os.path.join(RES, "ml", "ml_scenario_costs.csv"))
costs = sc[sc.key.str.startswith("scenario_")]["value"].astype(float).values
ecost = float(sc.loc[sc.key == "expected", "value"].iloc[0])
cvar  = float(sc.loc[sc.key == "cvar90", "value"].iloc[0])
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.hist(costs, bins=12, color="tab:blue", alpha=0.75, edgecolor="black", linewidth=0.5)
ax.axvline(ecost, color="tab:green", linestyle="-", linewidth=2, label=f"E[cost] = {ecost:,.0f}")
ax.axvline(cvar, color="tab:red", linestyle="--", linewidth=2, label=f"CVaR$_{{90}}$ = {cvar:,.0f}")
ax.set_xlabel("Total dispatch cost (EUR)")
ax.set_ylabel("Number of scenarios")
ax.legend(frameon=True, loc="upper left")
save(fig, "fig_ml_scenario_costs.png")

print("done.")
