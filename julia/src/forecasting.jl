"""
    forecasting.jl

LSTM-based load/wind forecaster with conformal prediction intervals.
Used as the AI component of the PowerFlowJulia thesis project.

Pipeline:
  data = generate_synthetic_data(365)       # synthetic 365-day load matrix
  fc   = train_forecaster(data)             # train LSTM, calibrate intervals
  pred = predict_scenarios(fc, history)     # point forecast + bounds + scenarios

The conformal prediction intervals provide distribution-free coverage guarantees:
with α=0.10, the interval [lower, upper] contains the true value with ≥90%
probability on data drawn from the same distribution as the calibration set.
"""

using Flux
using Statistics
using Random
using Printf

# ─────────────────────────────────────────────────────────────────────────────
#  Forecaster struct
# ─────────────────────────────────────────────────────────────────────────────

"""
    LoadForecaster

Trained LSTM forecaster with conformal calibration scores.
Created by `train_forecaster`; used via `predict_scenarios`.

Fields:
  model      — Flux Chain (Float32 LSTM weights)
  μ, σ       — z-score normalisation constants (Float32)
  cal_scores — sorted non-conformity scores from calibration set
  horizon    — number of hours predicted (default 24)
  window     — input window length in hours (default 24)
"""
struct LoadForecaster
    model                        # Flux Chain
    μ :: Float32
    σ :: Float32
    cal_scores :: Vector{Float32}
    horizon    :: Int
    window     :: Int
end

# ─────────────────────────────────────────────────────────────────────────────
#  Synthetic data generation
# ─────────────────────────────────────────────────────────────────────────────

"""
    generate_synthetic_data(n_days; profile, noise_std, seed) → Matrix{Float64}

Generates an (n_days × 24) matrix of synthetic hourly load data (pu).

Each row is one day. The daily shape follows `profile` (24-hour normalized curve)
with superimposed:
  - Gaussian noise (noise_std in pu)
  - Seasonal variation ±15% (winter peak in December)
  - Weekend demand dip −8% (Saturday/Sunday)

Used to create training data for `train_forecaster` when real historical
load data is unavailable.
"""
function generate_synthetic_data(n_days::Int = 365;
        profile   :: Vector{Float64} = DEFAULT_LOAD_PROFILE,
        noise_std :: Float64 = 0.05,
        seed      :: Int     = 42)

    rng  = MersenneTwister(seed)
    data = Matrix{Float64}(undef, n_days, 24)
    for d in 1:n_days
        seasonal = 1.0 + 0.15 * sin(2π * (d - 355) / 365)  # winter peak ~Jan
        weekend  = (d % 7 >= 6) ? 0.92 : 1.0               # Sat/Sun −8%
        noise    = randn(rng, 24) .* noise_std
        data[d, :] = clamp.(profile .* seasonal .* weekend .+ noise, 0.20, 1.50)
    end
    return data
end

# ─────────────────────────────────────────────────────────────────────────────
#  Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

function _znorm(x::AbstractVector{<:Real}, μ::Float32, σ::Float32)
    Float32.((x .- μ) ./ σ)
end

function _zdenorm(z::AbstractVector, μ::Float32, σ::Float32)
    Float64.(z .* σ .+ μ)
end

# LSTM + Dense forecasting head
function _build_model(hidden::Int, horizon::Int)
    Chain(
        LSTM(1 => hidden),
        Dense(hidden => horizon),
    )
end

# One forward pass over a sequence of scalars; returns the horizon-length output
# at the LAST timestep. Must be called with model already reset.
# Flux 0.14+ LSTM requires 2D input (features=1, batch=1).
function _run_sequence(m, x_norm::AbstractVector{Float32})
    outputs = [m(reshape(Float32[v], 1, 1)) for v in x_norm]
    return vec(outputs[end])           # shape (horizon,)
end

# ─────────────────────────────────────────────────────────────────────────────
#  Training
# ─────────────────────────────────────────────────────────────────────────────

"""
    train_forecaster(data; hidden, epochs, lr, window, horizon,
                     val_frac, cal_frac, seed, verbose) → LoadForecaster

Trains an LSTM on (n_days × 24) load data with:
  - Adam optimiser, early stopping (patience = 10 epochs)
  - Held-out validation set for model selection
  - Conformal calibration on a separate calibration set

The conformal non-conformity score is the maximum pointwise absolute error
over the 24-hour horizon, which gives a uniform prediction band.

Returns a `LoadForecaster` ready for `predict_scenarios`.
"""
function train_forecaster(data::Matrix{Float64};
        hidden   :: Int     = 32,
        epochs   :: Int     = 100,
        lr       :: Float64 = 1e-3,
        window   :: Int     = 24,
        horizon  :: Int     = 24,
        val_frac :: Float64 = 0.15,
        cal_frac :: Float64 = 0.10,
        seed     :: Int     = 42,
        verbose  :: Bool    = true)

    # 1 — Flatten and normalise ------------------------------------------------
    series = Float64.(vec(data'))          # column-major flatten → 1D time series
    μ = Float32(mean(series))
    σ = Float32(std(series));  σ < 1f-6 && (σ = 1f0)
    z = _znorm(series, μ, σ)

    # 2 — Sliding-window dataset -----------------------------------------------
    n_total = length(z) - window - horizon + 1
    n_total < 10 && error("Too little data: need ≥ $(window + horizon + 10) hours.")

    n_val = max(2, round(Int, n_total * val_frac))
    n_cal = max(2, round(Int, n_total * cal_frac))
    n_trn = n_total - n_val - n_cal

    X = [z[i           : i+window-1  ] for i in 1:n_total]
    Y = [z[i+window    : i+window+horizon-1] for i in 1:n_total]

    X_trn = X[1:n_trn];               Y_trn = Y[1:n_trn]
    X_val = X[n_trn+1:n_trn+n_val];   Y_val = Y[n_trn+1:n_trn+n_val]
    X_cal = X[n_trn+n_val+1:end];     Y_cal = Y[n_trn+n_val+1:end]

    # 3 — Model and optimiser --------------------------------------------------
    rng       = MersenneTwister(seed)
    model     = _build_model(hidden, horizon)
    opt_state = Flux.setup(Adam(lr), model)

    best_val   = Inf32
    best_state = Flux.state(model)
    patience   = 0

    for ep in 1:epochs
        # Stochastic gradient descent (one sample at a time, shuffled)
        for i in randperm(rng, n_trn)
            x, y = X_trn[i], Y_trn[i]
            Flux.reset!(model)
            _, grads = Flux.withgradient(model) do m
                ŷ = _run_sequence(m, x)
                mean((ŷ .- y).^2)
            end
            Flux.update!(opt_state, model, grads[1])
        end

        # Validation loss
        val_loss = Float32(mean(begin
            Flux.reset!(model)
            ŷ = _run_sequence(model, x)
            mean((ŷ .- y).^2)
        end for (x, y) in zip(X_val, Y_val)))

        if val_loss < best_val - 1f-5
            best_val   = val_loss
            best_state = deepcopy(Flux.state(model))
            patience   = 0
        else
            patience += 1
            if patience >= 10
                verbose && @printf("[Forecaster] early stop epoch %d  val_RMSE=%.4f\n",
                                   ep, sqrt(best_val))
                break
            end
        end

        if verbose && (ep == 1 || ep % 25 == 0)
            @printf("[Forecaster] epoch %3d/%d  val_RMSE=%.4f\n",
                    ep, epochs, sqrt(val_loss))
        end
    end

    Flux.loadmodel!(model, best_state)
    verbose && @printf("[Forecaster] done  best_val_RMSE=%.4f\n", sqrt(best_val))

    # 4 — Conformal calibration ------------------------------------------------
    cal_scores = Float32[]
    for (x, y) in zip(X_cal, Y_cal)
        Flux.reset!(model)
        ŷ = _run_sequence(model, x)
        push!(cal_scores, maximum(abs.(ŷ .- y)))   # max abs error over horizon
    end
    sort!(cal_scores)

    return LoadForecaster(model, μ, σ, cal_scores, horizon, window)
end

# ─────────────────────────────────────────────────────────────────────────────
#  Prediction with conformal intervals and scenario sampling
# ─────────────────────────────────────────────────────────────────────────────

"""
    predict_scenarios(fc, history; n_scenarios, α, seed) → NamedTuple

Given at least `fc.window` recent hourly observations (`history`, original units),
returns forecasts for the next `fc.horizon` hours.

Returns a NamedTuple with:
  mean      :: Vector{Float64} (horizon)    — LSTM point forecast
  lower     :: Vector{Float64} (horizon)    — conformal lower bound (1-α coverage)
  upper     :: Vector{Float64} (horizon)    — conformal upper bound
  scenarios :: Matrix{Float64} (horizon × n_scenarios) — sampled load profiles
  rmse_cal  :: Float64                      — calibration set RMSE

Conformal intervals guarantee ≥ (1-α) marginal coverage on i.i.d. data.
Scenarios are Gaussian samples clipped to [lower, upper].

PyPSA equivalent: passing different `load_profile` vectors to network.optimize()
"""
function predict_scenarios(fc::LoadForecaster, history::AbstractVector{Float64};
        n_scenarios :: Int     = 7,
        α           :: Float64 = 0.10,
        seed        :: Int     = 0)

    length(history) < fc.window &&
        error("history must have ≥ $(fc.window) values; got $(length(history))")

    rng   = MersenneTwister(seed)
    x_win = _znorm(collect(Float64, history[end-fc.window+1:end]), fc.μ, fc.σ)

    Flux.reset!(fc.model)
    ẑ = _run_sequence(fc.model, x_win)       # z-score forecast (Float32, length horizon)

    # Conformal quantile q at level ceil((1-α)(n+1)) / n
    n_cal = length(fc.cal_scores)
    q_idx = clamp(ceil(Int, (1 - α) * (n_cal + 1)), 1, n_cal)
    q     = fc.cal_scores[q_idx]             # Float32

    # Back to original scale
    ŷ     = _zdenorm(ẑ,     fc.μ, fc.σ)
    lower = _zdenorm(ẑ .- q, fc.μ, fc.σ)
    upper = _zdenorm(ẑ .+ q, fc.μ, fc.σ)
    lower = max.(lower, 0.0)
    upper = max.(upper, 0.0)

    # Sample n_scenarios profiles from N(ŷ, ((upper-lower)/4)²) clipped to [lower,upper]
    half_band = (upper .- lower) ./ 4.0
    scen = Matrix{Float64}(undef, fc.horizon, n_scenarios)
    for s in 1:n_scenarios
        scen[:, s] = clamp.(ŷ .+ randn(rng, fc.horizon) .* half_band,
                            lower, upper)
    end

    rmse_cal = Float64(sqrt(mean(fc.cal_scores .^ 2)))
    return (mean=ŷ, lower=lower, upper=upper, scenarios=scen, rmse_cal=rmse_cal)
end

# ─────────────────────────────────────────────────────────────────────────────
#  Forecast metrics
# ─────────────────────────────────────────────────────────────────────────────

"""
    forecast_metrics(actual, predicted) → NamedTuple

Computes MAE, RMSE, and MAPE between actual and predicted vectors.
"""
function forecast_metrics(actual::AbstractVector{Float64},
                          predicted::AbstractVector{Float64})
    n    = length(actual)
    err  = predicted .- actual
    mae  = mean(abs.(err))
    rmse = sqrt(mean(err .^ 2))
    mape = 100 * mean(abs.(err) ./ max.(abs.(actual), 1e-6))
    return (mae=mae, rmse=rmse, mape=mape, n=n)
end
