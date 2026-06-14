#=
ml_figures_data.jl — export the data behind the four AI figures, computed
honestly on REAL OPSD German load with a strict out-of-time test.

Produces (in results/):
  ml_sample_day.csv     hour, actual, mean, lower, upper   (one test day, weather-informed)
  ml_coverage.csv       hour, coverage   + a final overall row   (conformal coverage on the test block)
  ml_scenario_costs.csv scenario, cost   + E[cost]/CVaR90 rows   (stochastic LOPF on the sample day)

The accuracy bar chart uses the verified multi-seed numbers from
real_data_forecast.jl and is assembled directly in the Python plotting script.

Run:
  julia --project=. examples/ml_figures_data.jl
=#
include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia
using Statistics, Printf, Dates

dow0(date) = dayofweek(date) - 1          # Julia Mon=1..Sun=7 -> 0..6
RES = joinpath(@__DIR__, "..", "..", "results")

TRAIN_DAYS = 1095
TEST_DAYS  = 60
α          = 0.10

println("[1] Loading real OPSD load + temperature ...")
r = load_real_data(; normalize=:peak, max_days = TRAIN_DAYS + TEST_DAYS,
                   temp_path = DEFAULT_TEMP_CSV)
full, full_temp, first_day = r.data, r.temp, r.first_day
n = size(full, 1)
split = train_test_split_days(full; test_days=TEST_DAYS)
train      = split.train
train_temp = full_temp[1:n-TEST_DAYS, :]
test_start = n - TEST_DAYS + 1

println("[2] Training weather-informed LSTM (hidden=64, real data) ...")
fc = train_forecaster(train; hidden=64, epochs=80, window=168,
                      start_dow=dow0(first_day), temp=train_temp,
                      val_frac=0.15, cal_frac=0.10, seed=1, verbose=false)

# Predict one test day d (returns mean/lower/upper/scenarios), weather-informed.
function predict_day(d; n_scen=1)
    window  = vec(permutedims(full[d-7:d-1, :]))
    win_dow = dow0(first_day + Day(d - 7 - 1))
    th = vec(permutedims(full_temp[d-7:d-1, :]))
    tf = full_temp[d, :]
    predict_scenarios(fc, window; n_scenarios=n_scen, α=α, hist_start_dow=win_dow,
                      temp_history=th, temp_future=tf, seed=0)
end

# ── Sample day (middle of the test block) ────────────────────────────────────
println("[3] Sample-day forecast + conformal band ...")
d_sample = test_start + TEST_DAYS ÷ 2
p = predict_day(d_sample; n_scen=1)
actual = full[d_sample, :]
open(joinpath(RES, "ml_sample_day.csv"), "w") do io
    println(io, "hour,actual,mean,lower,upper")
    for h in 1:24
        @printf(io, "%d,%.5f,%.5f,%.5f,%.5f\n", h-1, actual[h], p.mean[h], p.lower[h], p.upper[h])
    end
end

# ── Conformal coverage over the whole test block ─────────────────────────────
println("[4] Conformal coverage on the test block ...")
covered = zeros(Int, 24); total = 0
for d in test_start:n
    pp = predict_day(d; n_scen=1)
    a  = full[d, :]
    for h in 1:24
        (a[h] >= pp.lower[h] - 1e-9 && a[h] <= pp.upper[h] + 1e-9) && (covered[h] += 1)
    end
    global total += 1
end
open(joinpath(RES, "ml_coverage.csv"), "w") do io
    println(io, "hour,coverage")
    for h in 1:24
        @printf(io, "%d,%.4f\n", h-1, covered[h] / total)
    end
    @printf(io, "overall,%.4f\n", sum(covered) / (24 * total))
end
@printf("    overall coverage = %.1f%%  (nominal %.0f%%)\n", 100*sum(covered)/(24*total), 100*(1-α))

# ── Stochastic LOPF on the sample day: cost distribution + CVaR ──────────────
println("[5] Stochastic LOPF on sample-day scenarios ...")
net = Network(baseMVA=100.0)
add!(net, "Bus", "B1"; v_nom=380.0, slack=true)
add!(net, "Bus", "B2"; v_nom=380.0)
add!(net, "Bus", "B3"; v_nom=380.0)
add!(net, "Line", "L12"; bus0="B1", bus1="B2", x=0.1, r=0.01, s_nom=400.0)
add!(net, "Line", "L13"; bus0="B1", bus1="B3", x=0.1, r=0.01, s_nom=400.0)
add!(net, "Line", "L23"; bus0="B2", bus1="B3", x=0.15, r=0.01, s_nom=300.0)
add!(net, "Generator", "Gas";  bus="B1", p_nom=500.0, marginal_cost=65.0)
add!(net, "Generator", "Coal"; bus="B1", p_nom=300.0, marginal_cost=30.0)
add!(net, "Generator", "Wind"; bus="B2", p_nom=200.0, marginal_cost=5.0, carrier="wind")
add!(net, "Load", "D2"; bus="B2", p_set=350.0)
add!(net, "Load", "D3"; bus="B3", p_set=250.0)
add!(net, "StorageUnit", "Bat"; bus="B3", p_nom=80.0, e_nom=320.0,
     efficiency_charge=0.92, efficiency_discharge=0.92)

ps   = predict_day(d_sample; n_scen=30)
res  = lopf_stochastic(net, ps.scenarios; T=24, verbose=false)
open(joinpath(RES, "ml_scenario_costs.csv"), "w") do io
    println(io, "key,value")
    for (i, c) in enumerate(res.costs)
        @printf(io, "scenario_%d,%.2f\n", i, c)
    end
    @printf(io, "expected,%.2f\n", res.expected_cost)
    @printf(io, "cvar90,%.2f\n",   res.cvar_90)
end
@printf("    E[cost]=%.0f  CVaR90=%.0f  spread=%.0f\n",
        res.expected_cost, res.cvar_90, res.std_cost)

# Deterministic baselines on the SAME network/day, for a consistent cost table:
#   static  = a flat default profile;  point = the LSTM mean forecast.
det_static = optimize(net; T=24, load_profile=DEFAULT_LOAD_PROFILE, verbose=false)
det_point  = optimize(net; T=24, load_profile=ps.mean,             verbose=false)
open(joinpath(RES, "ml_cost_table.csv"), "w") do io
    println(io, "strategy,cost")
    @printf(io, "static,%.2f\n",   det_static.total_cost)
    @printf(io, "point,%.2f\n",    det_point.total_cost)
    @printf(io, "expected,%.2f\n", res.expected_cost)
    @printf(io, "cvar90,%.2f\n",   res.cvar_90)
    @printf(io, "best,%.2f\n",     res.best_cost)
    @printf(io, "worst,%.2f\n",    res.worst_cost)
    @printf(io, "std,%.2f\n",      res.std_cost)
end
@printf("    static=%.0f  point=%.0f\n", det_static.total_cost, det_point.total_cost)
println("[OK] wrote ml_sample_day.csv, ml_coverage.csv, ml_scenario_costs.csv, ml_cost_table.csv")
