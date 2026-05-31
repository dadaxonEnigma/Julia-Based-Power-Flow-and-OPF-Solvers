"""
ml_forecast_demo.jl — AI-Assisted Power Dispatch Demo

Demonstrates the full ML forecasting → stochastic LOPF pipeline:

  1. Generate synthetic historical load data (365 days)
  2. Train LSTM forecaster with conformal prediction
  3. Predict next 24 hours: point forecast + 90% intervals + 7 scenarios
  4. Run deterministic LOPF with naive (average) profile
  5. Run deterministic LOPF with ML-predicted profile
  6. Run stochastic LOPF across 7 scenarios (SAA)
  7. Compare results and plot

Thesis framing:
  "AI-assisted" = LSTM forecasting feeds uncertainty-aware scenarios into LOPF,
  producing a stochastic dispatch that is robust to demand forecast errors.
"""

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia
using Plots
using StatsPlots
using Printf
using Statistics

outdir = joinpath(@__DIR__, "..", "..", "results", "plots")
isdir(outdir) || mkpath(outdir)

println("=" ^ 65)
println("  AI-Assisted Power Dispatch  —  PowerFlowJulia ML Demo")
println("=" ^ 65)

# ── 1. Build 3-bus test network ──────────────────────────────────────────────
println("\n[1] Building 3-bus network with storage and wind...")

net = Network(baseMVA=100.0)
add!(net, "Bus", "B1"; v_nom=380.0, slack=true)
add!(net, "Bus", "B2"; v_nom=380.0)
add!(net, "Bus", "B3"; v_nom=380.0)
add!(net, "Line", "L12"; bus0="B1", bus1="B2", x=0.1, r=0.01, s_nom=400.0)
add!(net, "Line", "L13"; bus0="B1", bus1="B3", x=0.1, r=0.01, s_nom=400.0)
add!(net, "Line", "L23"; bus0="B2", bus1="B3", x=0.15, r=0.01, s_nom=300.0)

# Generators
add!(net, "Generator", "Gas_B1";  bus="B1", p_nom=500.0, marginal_cost=65.0)
add!(net, "Generator", "Coal_B1"; bus="B1", p_nom=300.0, marginal_cost=30.0)
add!(net, "Generator", "Wind_B2"; bus="B2", p_nom=200.0, marginal_cost=5.0,  carrier="wind")

# Loads (nominal peak values; scaled by load_profile each hour)
add!(net, "Load", "D2"; bus="B2", p_set=350.0)
add!(net, "Load", "D3"; bus="B3", p_set=250.0)

# Battery storage at B3
add!(net, "StorageUnit", "Bat_B3";
     bus="B3", p_nom=80.0, e_nom=320.0,
     efficiency_charge=0.92, efficiency_discharge=0.92,
     standing_loss=0.001)

println("  Network: 3 buses, 3 lines, 3 generators, 2 loads, 1 storage")

# ── 2. Generate synthetic historical data ────────────────────────────────────
println("\n[2] Generating synthetic historical load data (365 days × 24h)...")
hist_data = generate_synthetic_data(365; noise_std=0.05, seed=42)
println(@sprintf("  Data: %.0f days  |  mean=%.3f pu  std=%.3f pu",
                 size(hist_data, 1), mean(hist_data), std(hist_data)))

# ── 3. Train LSTM forecaster ─────────────────────────────────────────────────
println("\n[3] Training LSTM forecaster (window=24h → horizon=24h)...")
fc = train_forecaster(hist_data;
     hidden=32, epochs=100, lr=1e-3,
     val_frac=0.15, cal_frac=0.10,
     seed=42, verbose=true)

println(@sprintf("  Conformal calibration: %d scores  |  q90=%.4f pu",
                 length(fc.cal_scores), fc.cal_scores[ceil(Int, 0.9*length(fc.cal_scores))]))

# ── 4. Predict next 24 hours with uncertainty ────────────────────────────────
println("\n[4] Predicting next 24h with conformal intervals (α=0.10)...")

# Use the LAST day of historical data as recent observation
last_day = vec(hist_data[end, :])
pred = predict_scenarios(fc, last_day; n_scenarios=7, α=0.10, seed=1)

println(@sprintf("  Point forecast : mean=%.3f  min=%.3f  max=%.3f (pu)",
                 mean(pred.mean), minimum(pred.mean), maximum(pred.mean)))
println(@sprintf("  Interval width : mean=%.3f pu  (90%% coverage guarantee)",
                 mean(pred.upper .- pred.lower)))

# ── 5. Run three LOPF variants ───────────────────────────────────────────────
T = 24
println("\n[5] Running LOPF variants (T=$T hours)...")

# Reference: deterministic with naive (average) profile
println("  5a. Deterministic LOPF — naive profile (DEFAULT_LOAD_PROFILE)...")
r_naive = optimize(net; T=T, method=:mp,
                   load_profile=DEFAULT_LOAD_PROFILE,
                   wind_profile=DEFAULT_WIND_PROFILE,
                   verbose=false)
println(@sprintf("      cost = %.2f €", r_naive.total_cost))

# ML-informed: deterministic with LSTM point forecast
println("  5b. Deterministic LOPF — ML point forecast...")
r_ml = optimize(net; T=T, method=:mp,
                load_profile=pred.mean,
                wind_profile=DEFAULT_WIND_PROFILE,
                verbose=false)
println(@sprintf("      cost = %.2f €", r_ml.total_cost))

# Stochastic: 7 scenarios from conformal intervals
println("  5c. Stochastic LOPF — 7 scenarios (SAA)...")
r_stoch = lopf_stochastic(net, pred.scenarios;
                          T=T, verbose=false)
println(@sprintf("      E[cost]=%.2f €  CVaR90=%.2f €  σ=%.2f €",
                 r_stoch.expected_cost, r_stoch.cvar_90, r_stoch.std_cost))

# ── 6. Comparison table ───────────────────────────────────────────────────────
println()
compare_deterministic_stochastic(r_ml, r_stoch)

Δ_naive_ml  = r_naive.total_cost - r_ml.total_cost
Δ_ml_stoch  = r_stoch.expected_cost - r_ml.total_cost
println(@sprintf("\n  Cost improvement (naive → ML forecast):  %.2f € (%+.1f%%)",
                 Δ_naive_ml, 100*Δ_naive_ml/max(r_naive.total_cost,1)))
println(@sprintf("  Stochastic premium (ML → SAA):           %.2f € (%+.1f%%)",
                 Δ_ml_stoch, 100*Δ_ml_stoch/max(r_ml.total_cost,1)))

# Forecast quality metrics
metrics = forecast_metrics(last_day, pred.mean)
println(@sprintf("\n  Forecast metrics (on last observed day):")         )
println(@sprintf("    MAE  = %.4f pu   RMSE = %.4f pu   MAPE = %.1f%%",
                 metrics.mae, metrics.rmse, metrics.mape))

# ── 7. Plots ─────────────────────────────────────────────────────────────────
println("\n[6] Generating plots...")

hours = 1:24

# (a) Forecast bands + scenarios
p1 = plot(hours, pred.mean;
    ribbon=(pred.mean .- pred.lower, pred.upper .- pred.mean),
    fillalpha=0.25, fillcolor=:steelblue,
    label="LSTM forecast (90% CI)", linewidth=2.5, color=:steelblue,
    xlabel="Hour", ylabel="Load (pu)",
    title="LSTM Load Forecast with Conformal Intervals",
    legend=:topright, dpi=200)
for s in 1:size(pred.scenarios, 2)
    plot!(p1, hours, pred.scenarios[:, s];
          label=(s==1 ? "Scenarios (n=$(size(pred.scenarios,2)))" : ""),
          color=:gray60, alpha=0.4, linewidth=1)
end
plot!(p1, hours, DEFAULT_LOAD_PROFILE;
      label="Naive (default)", color=:red, linestyle=:dash, linewidth=1.5)
savefig(p1, joinpath(outdir, "ml_forecast_bands.png"))
println("  Saved: results/plots/ml_forecast_bands.png")

# (b) Scenario cost distribution
p2 = histogram(r_stoch.costs;
    bins=max(3, length(r_stoch.costs)),
    label="Scenario costs",
    color=:steelblue, alpha=0.7,
    xlabel="Total cost (€)", ylabel="Count",
    title="Cost Distribution across $(length(r_stoch.costs)) Scenarios",
    dpi=200)
vline!(p2, [r_ml.total_cost];
       label="ML deterministic", color=:green, linewidth=2, linestyle=:dash)
vline!(p2, [r_stoch.expected_cost];
       label="E[cost] stochastic", color=:red, linewidth=2)
savefig(p2, joinpath(outdir, "ml_scenario_costs.png"))
println("  Saved: results/plots/ml_scenario_costs.png")

# (c) LMP comparison: naive vs ML-informed
bus_list = sort(collect(keys(r_naive.lmp)))
lmp_naive_mean = [mean(r_naive.lmp[b]) for b in bus_list]
lmp_ml_mean    = [mean(r_ml.lmp[b])    for b in bus_list]
lmp_stoch_mean = [get(r_stoch.mean_lmp, b, 0.0) for b in bus_list]

p3 = groupedbar(hcat(lmp_naive_mean, lmp_ml_mean, lmp_stoch_mean);
    bar_width=0.6,
    label=["Naive" "ML forecast" "Stochastic"],
    color=[:gray :steelblue :crimson], alpha=0.85,
    xticks=(1:length(bus_list), bus_list),
    xlabel="Bus", ylabel="Mean LMP (€/MWh)",
    title="Locational Marginal Prices: Naive vs ML vs Stochastic",
    legend=:topright, dpi=200)
savefig(p3, joinpath(outdir, "ml_lmp_comparison.png"))
println("  Saved: results/plots/ml_lmp_comparison.png")

# (d) ML dispatch vs naive dispatch over 24h (G_Gas)
if haskey(r_naive.gen_dispatch, "Gas_B1") && haskey(r_ml.gen_dispatch, "Gas_B1")
    p4 = plot(hours, r_naive.gen_dispatch["Gas_B1"];
        label="Gas dispatch — naive", color=:gray, linestyle=:dash, linewidth=2,
        xlabel="Hour", ylabel="Power (MW)",
        title="Gas Generator Dispatch: Naive vs ML-Informed",
        legend=:topright, dpi=200)
    plot!(p4, hours, r_ml.gen_dispatch["Gas_B1"];
        label="Gas dispatch — ML", color=:crimson, linewidth=2)
    savefig(p4, joinpath(outdir, "ml_dispatch_comparison.png"))
    println("  Saved: results/plots/ml_dispatch_comparison.png")
end

println("\n" * "=" ^ 65)
println("  ML Demo complete.")
println(@sprintf("  Naive → ML improvement  : %.2f €", Δ_naive_ml))
println(@sprintf("  Stochastic E[cost]      : %.2f €  (CVaR90=%.2f €)",
                 r_stoch.expected_cost, r_stoch.cvar_90))
println("=" ^ 65)
