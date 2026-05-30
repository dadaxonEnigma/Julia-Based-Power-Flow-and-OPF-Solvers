"""
Europe 10-bus demo — PowerFlowJulia vs PyPSA on identical network.
Runs multi-period LOPF (T=24), prints LMP summary, saves PNG plots.
"""

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia
using Printf

outdir = joinpath(@__DIR__, "..", "..", "results", "plots")
isdir(outdir) || mkpath(outdir)

# ── Bus definitions (same as python/europe_demo.py) ──────────────────────────
const BUSES = [
    ("London",    51.51, -0.13, 5600.0),
    ("Paris",     48.85,  2.35, 4750.0),
    ("Berlin",    52.52, 13.40, 4000.0),
    ("Madrid",    40.42, -3.70, 3500.0),
    ("Rome",      41.90, 12.50, 3250.0),
    ("Amsterdam", 52.37,  4.90, 2250.0),
    ("Warsaw",    52.23, 21.01, 2500.0),
    ("Vienna",    48.21, 16.37, 1875.0),
    ("Stockholm", 59.33, 18.07, 2125.0),
    ("Zurich",    47.38,  8.54, 1500.0),
]

const LINES = [
    ("London",    "Amsterdam", 2500.0, 0.08),
    ("London",    "Paris",     3000.0, 0.06),
    ("Paris",     "Berlin",    2000.0, 0.10),
    ("Paris",     "Madrid",    1500.0, 0.14),
    ("Paris",     "Zurich",    2000.0, 0.07),
    ("Berlin",    "Amsterdam", 1800.0, 0.07),
    ("Berlin",    "Warsaw",    1500.0, 0.09),
    ("Berlin",    "Vienna",    2000.0, 0.09),
    ("Berlin",    "Stockholm", 1200.0, 0.13),
    ("Zurich",    "Rome",      1800.0, 0.09),
    ("Zurich",    "Vienna",    2000.0, 0.06),
    ("Vienna",    "Warsaw",    1200.0, 0.10),
    ("Vienna",    "Rome",      1500.0, 0.10),
]

# (city, fuel, marginal_cost, p_nom_MW)
const GENERATORS = [
    ("London",    "nuclear",  20.0, 12000.0),
    ("London",    "gas",      65.0,  8000.0),
    ("Paris",     "nuclear",  18.0, 20000.0),   # slack bus
    ("Paris",     "hydro",    10.0,  5000.0),
    ("Berlin",    "coal",     40.0, 10000.0),
    ("Berlin",    "wind",      5.0, 15000.0),
    ("Berlin",    "gas",      70.0,  5000.0),
    ("Madrid",    "gas",      55.0,  8000.0),
    ("Madrid",    "solar",     3.0, 12000.0),
    ("Rome",      "gas",      60.0,  7000.0),
    ("Rome",      "hydro",    12.0,  4000.0),
    ("Amsterdam", "gas",      55.0,  5000.0),
    ("Amsterdam", "wind",      5.0,  6000.0),
    ("Warsaw",    "coal",     45.0,  8000.0),
    ("Warsaw",    "gas",      75.0,  3000.0),
    ("Vienna",    "hydro",    10.0,  4000.0),
    ("Vienna",    "gas",      65.0,  3000.0),
    ("Stockholm", "hydro",     8.0,  8000.0),
    ("Stockholm", "nuclear",  20.0,  5000.0),
    ("Zurich",    "hydro",     9.0,  3000.0),
    ("Zurich",    "nuclear",  22.0,  2000.0),
]

const LOAD_PROFILE_24 = DEFAULT_LOAD_PROFILE

# ── Build network ─────────────────────────────────────────────────────────────
net = Network(baseMVA=1000.0)   # baseMVA=1000 MVA for large European system

for (city, lat, lon, _) in BUSES
    is_slack = city == "Paris"
    add!(net, "Bus", city; v_nom=380.0, slack=is_slack)
end

for (f, t, s_nom, x) in LINES
    lname = "$(f[1:3])-$(t[1:3])"
    add!(net, "Line", lname; bus0=f, bus1=t, x=x, r=0.01, s_nom=s_nom)
end

for (city, fuel, mc, pnom) in GENERATORS
    add!(net, "Generator", "$(city[1:3])_$fuel";
         bus=city, p_nom=pnom, marginal_cost=mc)
end

for (city, _, _, peak_load) in BUSES
    add!(net, "Load", "D_$(city[1:3])"; bus=city, p_set=peak_load)
end

# ── Run multi-period LOPF ─────────────────────────────────────────────────────
println("Running Julia LOPF (T=24, 10-bus European network) ...")
t0 = time()
r = lopf_multiperiod(net; T=24, load_profile=LOAD_PROFILE_24,
                      wind_profile=LOAD_PROFILE_24, verbose=false)
elapsed = (time() - t0) * 1000

@printf("\n  Solved in %.1f ms  |  Total cost: %.2f M€\n",
        elapsed, r.total_cost / 1e6)

# ── LMP summary ───────────────────────────────────────────────────────────────
println("\nLMP summary (€/MWh):")
@printf("  %-12s  %8s  %8s  %8s\n", "Bus", "Mean", "Min", "Max")
for (city, _, _, _) in BUSES
    if haskey(r.lmp, city)
        v = r.lmp[city]
        @printf("  %-12s  %8.2f  %8.2f  %8.2f\n",
                city, sum(v)/length(v), minimum(v), maximum(v))
    end
end

# ── Dispatch summary ──────────────────────────────────────────────────────────
println("\nTop generators by total energy (MWh/day):")
totals = [(n, sum(v)) for (n, v) in r.gen_dispatch]
sort!(totals, by=x->x[2], rev=true)
for (n, e) in totals[1:min(8, length(totals))]
    @printf("  %-20s  %10.0f MWh\n", n, e)
end

# ── Plots ─────────────────────────────────────────────────────────────────────
println("\nGenerating plots ...")

p1 = plot_dispatch(r; net=net,
                   title="European Grid — Generator Dispatch (T=24h)",
                   savepath=joinpath(outdir, "europe_dispatch.png"))
println("  saved europe_dispatch.png")

p2 = plot_lmp(r;
              title="European Grid — LMP Heatmap (€/MWh, T=24h)",
              savepath=joinpath(outdir, "europe_lmp.png"))
println("  saved europe_lmp.png")

p3 = plot_network(net;
                  title="European 10-Bus Network Topology",
                  savepath=joinpath(outdir, "europe_network.png"))
println("  saved europe_network.png")

println("\n[OK] Done. Open results/europe_map.html for interactive map (run python/europe_demo.py first).")
