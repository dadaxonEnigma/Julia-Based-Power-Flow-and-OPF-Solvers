"""
make_benchmark_figures.py — generate publication-quality benchmark figures and a
consolidated Julia-vs-PyPSA speedup table from the results/*.csv files.

Academic conventions applied:
  - no in-figure title (the description goes in the LaTeX caption)
  - axis labels carry units in parentheses
  - distinct markers AND line styles so the figures read in grayscale
  - serif font, 300 dpi, legends with frames
  - the reported metric is the minimum solve time (min_ms), matching the thesis

Outputs -> results/plots/*.png  and  results/benchmarks/speedup_summary.csv

Run:
  pypsa_env/Scripts/python.exe python/make_benchmark_figures.py
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "lines.linewidth": 1.8,
    "lines.markersize": 6,
})
# Julia = blue circle solid; PyPSA = red square dashed (distinct in grayscale).
JL = dict(color="tab:blue", marker="o", linestyle="-",  label="Julia (PowerFlowJulia)")
PY = dict(color="tab:red",  marker="s", linestyle="--", label="Python (PyPSA)")
RES = os.path.join(os.path.dirname(__file__), "..", "..", "results")
PLT = os.path.join(RES, "plots")
os.makedirs(PLT, exist_ok=True)


def load(name):
    p = os.path.join(RES, "benchmarks", name)
    return pd.read_csv(p) if os.path.isfile(p) else None


def save(fig, name):
    path = os.path.join(PLT, name)
    fig.savefig(path)
    plt.close(fig)
    print("  saved", os.path.relpath(path, RES))


def panel_label(ax, text):
    # Left-aligned title placed ABOVE the axes so it never overlaps the data.
    ax.set_title(text, loc="left", fontweight="bold", fontsize=11)


def add_headroom(ax, factor=3.0):
    # Expand the top of a log y-axis so a near-flat top curve is not glued to
    # the frame (and does not collide with the legend).
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi * factor)


speedup_rows = []


def add_speedup(method, scenario, jl, py):
    if jl and py and jl > 0:
        speedup_rows.append((method, scenario, jl, py, py / jl))


# ── 1. DC PF & LOPF scaling vs network size (synthetic) ──────────────────────
jb, pb = load("julia_benchmark.csv"), load("python_benchmark.csv")
if jb is not None and pb is not None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, mod, lbl in zip(axes, ["DC_PF", "LOPF"], ["(a) DC power flow", "(b) LOPF"]):
        j, p = jb[jb.module == mod], pb[pb.module == mod]
        ax.plot(j.n_buses, j.min_ms, **JL)
        ax.plot(p.n_buses, p.min_ms, **PY)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Number of buses")
        ax.set_ylabel("Minimum solve time (ms)")
        ax.legend(frameon=True)
        panel_label(ax, lbl)
        for _, r in j.merge(p, on="n_buses", suffixes=("_jl", "_py")).iterrows():
            add_speedup(mod, f"{int(r.n_buses)} buses", r.min_ms_jl, r.min_ms_py)
    fig.tight_layout()
    save(fig, "fig_scaling_dcpf_lopf.png")

    # Speedup vs size
    fig, ax = plt.subplots(figsize=(6, 4))
    styles = {"DC_PF": dict(color="tab:blue", marker="o", linestyle="-", label="DC power flow"),
              "LOPF":  dict(color="tab:green", marker="^", linestyle="-.", label="LOPF")}
    for mod, st in styles.items():
        m = jb[jb.module == mod].merge(pb[pb.module == mod], on="n_buses",
                                       suffixes=("_jl", "_py"))
        ax.plot(m.n_buses, m.min_ms_py / m.min_ms_jl, **st)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Number of buses")
    ax.set_ylabel("Speedup (PyPSA time / Julia time)")
    ax.axhline(1, color="gray", linestyle="--", linewidth=1)
    ax.legend(frameon=True)
    fig.tight_layout()
    save(fig, "fig_speedup_size.png")

# ── 2. AC power flow scaling (honest crossover) ──────────────────────────────
ja, pa = load("julia_ac_benchmark.csv"), load("python_ac_benchmark.csv")
if ja is not None and pa is not None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(ja.n_buses, ja.min_ms, **dict(JL, label="Julia (PowerModels.jl + Ipopt)"))
    ax.plot(pa.n_buses, pa.min_ms, **dict(PY, label="Python (PyPSA, Newton-Raphson)"))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Number of buses")
    ax.set_ylabel("Minimum solve time (ms)")
    ax.legend(frameon=True)
    fig.tight_layout()
    save(fig, "fig_scaling_acpf.png")
    for _, r in ja.merge(pa, on="n_buses", suffixes=("_jl", "_py")).iterrows():
        add_speedup("AC_PF", f"{int(r.n_buses)} buses", r.min_ms_jl, r.min_ms_py)

# ── 3. Multi-period LOPF & UC vs horizon T ───────────────────────────────────
jm, pm = load("julia_mp_benchmark.csv"), load("python_mp_benchmark.csv")
ju, pu = load("julia_uc_benchmark.csv"), load("python_uc_benchmark.csv")
panels = [(jm, pm, "(a) Multi-period LOPF"), (ju, pu, "(b) Unit commitment")]
if all(x is not None for x in (jm, pm, ju, pu)):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, (jd, pd_, lbl) in zip(axes, panels):
        ax.plot(jd["T"], jd.min_ms, **JL)
        ax.plot(pd_["T"], pd_.min_ms, **PY)
        ax.set_yscale("log")
        ax.set_xlabel("Planning horizon T (hours)")
        ax.set_ylabel("Minimum solve time (ms)")
        panel_label(ax, lbl)
        add_headroom(ax, 4.0)
        ax.legend(frameon=True, loc="center right")
        mod = "MP_LOPF" if "LOPF" in lbl else "UC_T"
        for _, r in jd.merge(pd_, on="T", suffixes=("_jl", "_py")).iterrows():
            add_speedup(mod, f"T={int(r['T'])}", r.min_ms_jl, r.min_ms_py)
    fig.tight_layout()
    save(fig, "fig_scaling_mp_uc.png")

jus, pus = load("julia_uc_scale.csv"), load("python_uc_scale.csv")
if jus is not None and pus is not None:
    for _, r in jus.merge(pus, on="n_buses", suffixes=("_jl", "_py")).iterrows():
        add_speedup("UC_buses", f"{int(r.n_buses)} buses", r.min_ms_jl, r.min_ms_py)

# ── 4. Standard real networks (IEEE / PEGASE) ────────────────────────────────
jr, pr = load("julia_realcases_benchmark.csv"), load("pypsa_realcases_benchmark.csv")
if jr is not None and pr is not None:
    import numpy as np
    order = jr.drop_duplicates("case")[["case", "n_buses"]].sort_values("n_buses")
    labels = [c.replace(".m", "").replace("case", "").replace("pegase", " PEGASE")
              for c in order.case]
    labels = [f"{l}\n({int(n)} buses)" for l, n in zip(labels, order.n_buses)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, mod, lbl in zip(axes, ["DC_PF", "LOPF"], ["(a) DC power flow", "(b) LOPF"]):
        jm_ = order.merge(jr[jr.module == mod][["case", "min_ms"]], on="case").rename(columns={"min_ms": "jl"})
        pm_ = jm_.merge(pr[pr.module == mod][["case", "min_ms"]], on="case").rename(columns={"min_ms": "py"})
        x = np.arange(len(pm_))
        ax.bar(x - 0.2, pm_.jl, 0.4, color="tab:blue", label="Julia", edgecolor="black", linewidth=0.5)
        ax.bar(x + 0.2, pm_.py, 0.4, color="tab:red", label="PyPSA", hatch="//", edgecolor="black", linewidth=0.5)
        ax.set_yscale("log")
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("Minimum solve time (ms)")
        ax.legend(frameon=True)
        panel_label(ax, lbl)
        for _, r in pm_.iterrows():
            add_speedup(f"real_{mod}", r.case.replace(".m", ""), r.jl, r.py)
    fig.tight_layout()
    save(fig, "fig_realcases.png")

# ── Consolidated speedup table ───────────────────────────────────────────────
df = pd.DataFrame(speedup_rows, columns=["method", "scenario", "julia_ms", "pypsa_ms", "speedup"])
df.to_csv(os.path.join(RES, "benchmarks", "speedup_summary.csv"), index=False, float_format="%.3f")
print("\nspeedup_summary.csv  ({} rows)".format(len(df)))
if len(df):
    print("  speedup range: {:.1f} .. {:.1f}  (median {:.1f})".format(
        df.speedup.min(), df.speedup.max(), df.speedup.median()))
