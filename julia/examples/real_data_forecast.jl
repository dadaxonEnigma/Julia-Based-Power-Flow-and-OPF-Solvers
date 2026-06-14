"""
real_data_forecast.jl — honest LSTM forecasting on REAL load data (OPSD/ENTSO-E)

Unlike `ml_forecast_demo.jl` (synthetic, self-generated data), this script
trains and evaluates on real German hourly load and reports metrics on a
strictly out-of-time hold-out, compared against two naive baselines:

  • Persistence      : tomorrow = today
  • Seasonal-naive   : tomorrow = same weekday last week (day − 7)

The LSTM is only scientifically justified if it beats these baselines on the
held-out test period. This is the experiment to cite in Chapter 4 (Results).

Run:
  julia julia/scripts/download_opsd.jl     # once, to fetch data
  julia julia/examples/real_data_forecast.jl
"""

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia
using Statistics
using Printf
using Dates

# ── Config (override via env vars) ──────────────────────────────────────────
TRAIN_HORIZON_DAYS = parse(Int, get(ENV, "PF_TRAIN_DAYS", "1095"))  # ~3 years
TEST_DAYS          = parse(Int, get(ENV, "PF_TEST_DAYS",  "60"))
HIDDEN             = parse(Int, get(ENV, "PF_HIDDEN",     "64"))
EPOCHS             = parse(Int, get(ENV, "PF_EPOCHS",     "80"))
USE_TEMP           = get(ENV, "PF_TEMP", "1") == "1"   # weather-informed mode

dow0(date) = dayofweek(date) - 1     # Julia: Mon=1..Sun=7  →  0..6

println("="^65)
println("  Real-Data LSTM Load Forecast  —  OPSD (Germany, ENTSO-E)")
println("  Temperature feature: $(USE_TEMP ? "ON (weather-informed)" : "OFF")")
println("="^65)

# ── 1. Load real data ───────────────────────────────────────────────────────
println("\n[1] Loading real hourly load data...")
r = load_real_data(; normalize=:peak, max_days = TRAIN_HORIZON_DAYS + TEST_DAYS,
                   temp_path = USE_TEMP ? DEFAULT_TEMP_CSV : nothing)
full      = r.data
full_temp = r.temp
n         = size(full, 1)
first_day = r.first_day

tt_split = train_test_split_days(full; test_days=TEST_DAYS)
train      = tt_split.train
train_temp = USE_TEMP ? full_temp[1:n-TEST_DAYS, :] : nothing
@printf("    train: %d days   test (out-of-time): %d days   (start %s)\n",
        size(train,1), TEST_DAYS, string(first_day))

# Seeds for the multi-seed evaluation (neural-net init is random — a single run
# is noisy, so we report mean ± std across seeds).
SEEDS = parse.(Int, split(get(ENV, "PF_SEEDS", "1,2,3,4,5"), ","))
test_start = n - TEST_DAYS + 1
rmse(e) = sqrt(mean(e .^ 2))

# ── 2. Seed-independent naive baselines (computed once) ──────────────────────
seas_err = Float64[]; seas_ape = Float64[]; pers_err = Float64[]; pers_ape = Float64[]
for d in test_start:n
    a = full[d, :]
    append!(seas_err, abs.(full[d-7, :] .- a)); append!(seas_ape, abs.(full[d-7, :] .- a) ./ a)
    append!(pers_err, abs.(full[d-1, :] .- a)); append!(pers_ape, abs.(full[d-1, :] .- a) ./ a)
end
seas_rmse = rmse(seas_err); seas_mape = 100 * mean(seas_ape)
pers_rmse = rmse(pers_err); pers_mape = 100 * mean(pers_ape)

# Roll one trained forecaster over the test block → (RMSE, MAPE).
function eval_forecaster(fc)
    err = Float64[]; ape = Float64[]
    for d in test_start:n
        a       = full[d, :]
        window  = vec(full[d-7:d-1, :]')
        win_dow = dow0(first_day + Day(d - 7 - 1))
        # Weather-informed: actual target-day temperature proxies a perfect
        # weather forecast (an idealisation that bounds the gain — stated honestly).
        th = USE_TEMP ? vec(full_temp[d-7:d-1, :]') : nothing
        tf = USE_TEMP ? full_temp[d, :]             : nothing
        p  = predict_scenarios(fc, window; n_scenarios=1, α=0.10, hist_start_dow=win_dow,
                               temp_history=th, temp_future=tf, seed=0)
        append!(err, abs.(p.mean .- a)); append!(ape, abs.(p.mean .- a) ./ a)
    end
    return rmse(err), 100 * mean(ape)
end

# ── 3. Train + evaluate across seeds ─────────────────────────────────────────
println("\n[2] Training LSTM over $(length(SEEDS)) seeds (weekly window, calendar" *
        (USE_TEMP ? " + temperature" : "") * "; test never seen)...")
lstm_rmses = Float64[]; lstm_mapes = Float64[]; skills = Float64[]
for sd in SEEDS
    fc = train_forecaster(train; hidden=HIDDEN, epochs=EPOCHS, window=168,
                          start_dow=dow0(first_day), temp=train_temp,
                          val_frac=0.15, cal_frac=0.10, seed=sd, verbose=false)
    rm, mp = eval_forecaster(fc)
    push!(lstm_rmses, rm); push!(lstm_mapes, mp); push!(skills, 100 * (1 - rm / seas_rmse))
    @printf("    seed %-3d  RMSE=%.4f  MAPE=%.2f%%  skill_vs_seasonal=%+.1f%%\n",
            sd, rm, mp, skills[end])
end

# ── 4. Report (mean ± std across seeds) ──────────────────────────────────────
println("\n[3] Out-of-time results on $(TEST_DAYS) days (pu), $(length(SEEDS)) seeds:")
println("    " * "-"^54)
@printf("    %-16s  %14s  %16s\n", "Model", "RMSE", "MAPE %")
println("    " * "-"^54)
@printf("    %-16s  %6.4f ± %.4f  %7.2f ± %.2f\n", "LSTM",
        mean(lstm_rmses), std(lstm_rmses), mean(lstm_mapes), std(lstm_mapes))
@printf("    %-16s  %6.4f          %7.2f\n", "Seasonal-naive", seas_rmse, seas_mape)
@printf("    %-16s  %6.4f          %7.2f\n", "Persistence",    pers_rmse, pers_mape)
println("    " * "-"^54)
@printf("    LSTM skill vs seasonal-naive (RMSE): %+.1f ± %.1f %%\n", mean(skills), std(skills))
@printf("    LSTM skill vs persistence    (RMSE): %+.1f %%\n", 100 * (1 - mean(lstm_rmses) / pers_rmse))
println("    (positive ⇒ LSTM beats the baseline)")

println("\n" * "="^65)
println("  Done. $(length(SEEDS)) seeds × $(TEST_DAYS) out-of-time days.")
println("="^65)
