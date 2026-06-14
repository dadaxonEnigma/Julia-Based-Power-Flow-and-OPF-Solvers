"""
    forecasting.jl

LSTM-based load/wind forecaster with conformal prediction intervals.
AI component of the PowerFlowJulia thesis project.

Pipeline:
  data = load_real_data().data              # real (n_days × 24) load matrix
  fc   = train_forecaster(data)             # train LSTM, calibrate intervals
  pred = predict_scenarios(fc, history)     # point forecast + bounds + scenarios

Architecture (v2 — sequence-to-vector, multi-feature):
  • Input window = 168 h (one week) so the LSTM sees the full weekly cycle that
    dominates electricity demand (the weekday/weekend pattern a naive
    seasonal model exploits).
  • Per-hour features = [load(z), sin/cos hour-of-day, sin/cos day-of-week,
    weekend flag] — calendar context lets the model know *where* in the week it
    is without having to infer it from the load shape alone.
  • Training is **batched over full sequences** via the 3-D recurrent call
    `enc(X)::(hidden, batch, time)`, taking the last timestep. This is both
    correct (the LSTM actually carries state across the window) and fast enough
    to train on multi-year data — unlike a per-timestep loop, which under
    Flux's stateless layers degenerates to a memoryless model.

The conformal prediction intervals provide distribution-free coverage: with
α=0.10 the band [lower, upper] contains the true value with ≥90% probability on
data drawn from the calibration distribution.
"""

using Flux
using Statistics
using Random
using Printf

# ─────────────────────────────────────────────────────────────────────────────
#  Sequence-to-vector model wrapper
# ─────────────────────────────────────────────────────────────────────────────

# enc: LSTM(F => hidden) processes (F, B, W) → (hidden, B, W); we keep the last
# timestep and map it to the horizon with a small MLP head.
struct Seq2Vec{E, H}
    enc  :: E
    head :: H
end
Flux.@layer Seq2Vec

# Sequence only (no exogenous future input).
(m::Seq2Vec)(x::AbstractArray{<:Real,3}) = m.head(m.enc(x)[:, :, end])

# Sequence + known future exogenous (e.g. target-day temperature): concatenate
# the future vector to the LSTM's last hidden state before the head.
(m::Seq2Vec)(x::AbstractArray{<:Real,3}, fut::AbstractMatrix) =
    m.head(vcat(m.enc(x)[:, :, end], fut))

# ─────────────────────────────────────────────────────────────────────────────
#  Forecaster struct
# ─────────────────────────────────────────────────────────────────────────────

"""
    LoadForecaster

Trained LSTM forecaster with conformal calibration scores.
Created by `train_forecaster`; used via `predict_scenarios`.

Fields:
  model      — Seq2Vec (LSTM encoder + MLP head, Float32)
  μ, σ       — z-score normalisation constants for load (Float32)
  cal_scores — sorted non-conformity scores from the calibration set (z-space)
  horizon    — number of hours predicted (default 24)
  window     — input window length in hours (default 168 = one week)
  n_features — input features per timestep (load + calendar)
  start_dow  — day-of-week of the first training day (0=Mon … 6=Sun)
"""
struct LoadForecaster
    model
    μ          :: Float32
    σ          :: Float32
    cal_scores :: Vector{Float32}
    horizon    :: Int
    window     :: Int
    n_features :: Int
    start_dow  :: Int
    residual   :: Bool        # if true, model predicts the correction to the
                              # seasonal-naive baseline (load one week earlier)
    has_temp   :: Bool        # temperature feature + future-temp head input
    μt         :: Float32     # temperature normalisation (°C)
    σt         :: Float32
end

# Seasonal lag used by the residual baseline: one week, capped by the window.
_seasonal_lag(window::Int) = min(168, window)

# ─────────────────────────────────────────────────────────────────────────────
#  Synthetic data generation (fallback when real data is unavailable)
# ─────────────────────────────────────────────────────────────────────────────

"""
    generate_synthetic_data(n_days; profile, noise_std, seed) → Matrix{Float64}

Generates an (n_days × 24) matrix of synthetic hourly load data (pu).

Each row is one day, following `profile` (24-hour normalized curve) with
Gaussian noise, ±15% seasonal variation (winter peak) and a weekend dip.
Used as a fallback for `train_forecaster` when real historical load data
(see `load_real_data`) is unavailable.
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

_znorm(x::AbstractVector{<:Real}, μ::Float32, σ::Float32) = Float32.((x .- μ) ./ σ)
_zdenorm(z::AbstractVector, μ::Float32, σ::Float32)        = Float64.(z .* σ .+ μ)

# Calendar features for `npos` consecutive hours starting at (start_hour, start_dow).
# Rows: [sin h, cos h, sin dow, cos dow, weekend]   (5 × npos), all Float32.
function _calendar_matrix(npos::Int, start_hour::Int, start_dow::Int)
    C = Matrix{Float32}(undef, 5, npos)
    for p in 0:npos-1
        tot = start_hour + p
        h   = tot % 24
        dow = (start_dow + tot ÷ 24) % 7
        C[1, p+1] = sin(2f0π * h   / 24f0)
        C[2, p+1] = cos(2f0π * h   / 24f0)
        C[3, p+1] = sin(2f0π * dow / 7f0)
        C[4, p+1] = cos(2f0π * dow / 7f0)
        C[5, p+1] = dow >= 5 ? 1f0 : 0f0     # Sat(5)/Sun(6)
    end
    return C
end

# Full feature matrix (F × L) over a z-normalised load series. F = 6 (load +
# 5 calendar), or 7 when a z-normalised temperature series `temp_z` is supplied.
function _feature_matrix(z::AbstractVector{Float32}, start_dow::Int;
                         temp_z::Union{Nothing,AbstractVector{Float32}} = nothing)
    L   = length(z)
    cal = _calendar_matrix(L, 0, start_dow)
    base = vcat(reshape(z, 1, L), cal)        # (6 × L)
    return temp_z === nothing ? base : vcat(base, reshape(temp_z, 1, L))  # (7 × L)
end

# Build the batched sliding-window dataset: X (F × n × W), Y (horizon × n).
function _make_dataset(feat::Matrix{Float32}, z::Vector{Float32},
                       window::Int, horizon::Int)
    F = size(feat, 1)
    L = size(feat, 2)
    n = L - window - horizon + 1
    X = Array{Float32,3}(undef, F, n, window)
    Y = Array{Float32,2}(undef, horizon, n)
    for i in 1:n
        @views X[:, i, :] .= feat[:, i:i+window-1]
        @views Y[:, i]    .= z[i+window : i+window+horizon-1]
    end
    return X, Y
end

# LSTM encoder + 2-layer MLP head. `head_extra` widens the head input to accept
# known future exogenous values (e.g. target-day temperature) appended to the
# LSTM's last hidden state.
function _build_model(n_features::Int, hidden::Int, horizon::Int; head_extra::Int = 0)
    Seq2Vec(
        LSTM(n_features => hidden),
        Chain(Dense(hidden + head_extra => hidden, relu), Dense(hidden => horizon)),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
#  Training
# ─────────────────────────────────────────────────────────────────────────────

"""
    train_forecaster(data; hidden, epochs, lr, window, horizon, batch_size,
                     start_dow, val_frac, cal_frac, seed, verbose) → LoadForecaster

Trains an LSTM on (n_days × 24) load data with:
  - Adam, mini-batched full-sequence training, early stopping (patience = 12)
  - Held-out validation set (chronological) for model selection
  - Conformal calibration on a separate, later calibration set

`start_dow` is the day-of-week (0=Mon … 6=Sun) of `data[1, :]`; pass it from the
real-data loader so the calendar features are aligned to actual dates.

The conformal non-conformity score is the maximum pointwise absolute error
over the horizon (z-space), giving a uniform prediction band.
"""
function train_forecaster(data::Matrix{Float64};
        hidden     :: Int     = 64,
        epochs     :: Int     = 80,
        lr         :: Float64 = 1e-3,
        window     :: Int     = 168,
        horizon    :: Int     = 24,
        batch_size :: Int     = 64,
        start_dow  :: Int     = 0,
        residual   :: Bool    = true,
        temp       :: Union{Nothing,Matrix{Float64}} = nothing,
        val_frac   :: Float64 = 0.15,
        cal_frac   :: Float64 = 0.10,
        seed       :: Int     = 42,
        verbose    :: Bool    = true)

    # 1 — Flatten chronologically and normalise (load only) --------------------
    series = Float64.(vec(data'))          # row-major flatten → 1-D hourly series
    μ = Float32(mean(series))
    σ = Float32(std(series));  σ < 1f-6 && (σ = 1f0)
    z = _znorm(series, μ, σ)

    n_total = length(z) - window - horizon + 1
    n_total < 10 && error("Too little data: need ≥ $(window + horizon + 10) hours.")

    # Optional exogenous temperature (same n_days × 24 layout as load).
    has_temp = temp !== nothing
    μt = σt = 1f0
    temp_z = nothing
    if has_temp
        size(temp) == size(data) || error("temp must have the same shape as data")
        tser = Float64.(vec(temp'))
        μt = Float32(mean(tser)); σt = Float32(std(tser)); σt < 1f-6 && (σt = 1f0)
        temp_z = _znorm(tser, μt, σt)
    end

    feat = _feature_matrix(z, start_dow; temp_z=temp_z)
    X, Y = _make_dataset(feat, z, window, horizon)

    # Known future (target-day) temperature for the head, when available.
    FUT = Array{Float32,2}(undef, has_temp ? horizon : 0, size(X, 2))
    if has_temp
        for i in 1:size(X, 2)
            @views FUT[:, i] .= temp_z[i+window : i+window+horizon-1]
        end
    end

    # Seasonal-naive baseline per sample: load one week (L hours) before the
    # target. With window=168 this is the first 24 h of the input window.
    L  = _seasonal_lag(window)
    SN = Array{Float32,2}(undef, horizon, size(X, 2))
    for i in 1:size(X, 2)
        @views SN[:, i] .= z[i+window-L : i+window-L+horizon-1]
    end
    # Training target: residual (actual − seasonal-naive) or the level itself.
    Yt = residual ? (Y .- SN) : Y

    # 2 — Chronological train / val / cal split --------------------------------
    n_val = max(2, round(Int, n_total * val_frac))
    n_cal = max(2, round(Int, n_total * cal_frac))
    n_trn = n_total - n_val - n_cal
    n_trn < 1 && error("Not enough windows after split; provide more data.")

    tr = 1:n_trn
    va = n_trn+1 : n_trn+n_val
    ca = n_trn+n_val+1 : n_total

    Xva, Yva, Fva = X[:, va, :], Yt[:, va], FUT[:, va]

    # 3 — Model and optimiser --------------------------------------------------
    rng       = MersenneTwister(seed)
    model     = _build_model(size(feat, 1), hidden, horizon;
                             head_extra = has_temp ? horizon : 0)
    opt_state = Flux.setup(Adam(lr), model)

    # Forward that branches on whether future temperature is supplied.
    fwd(m, xb, fb) = has_temp ? m(xb, fb) : m(xb)

    best_val   = Inf32
    best_state = deepcopy(Flux.state(model))
    patience   = 0

    for ep in 1:epochs
        idx = randperm(rng, n_trn)
        for s in 1:batch_size:n_trn
            b   = idx[s : min(s+batch_size-1, n_trn)]
            Xb  = X[:, tr[b], :]
            Yb  = Yt[:, tr[b]]
            Fb  = FUT[:, tr[b]]
            gs  = Flux.gradient(m -> mean(abs2, fwd(m, Xb, Fb) .- Yb), model)
            Flux.update!(opt_state, model, gs[1])
        end

        val_loss = mean(abs2, fwd(model, Xva, Fva) .- Yva)
        if val_loss < best_val - 1f-5
            best_val   = val_loss
            best_state = deepcopy(Flux.state(model))
            patience   = 0
        else
            patience += 1
            if patience >= 12
                verbose && @printf("[Forecaster] early stop epoch %d  val_RMSE=%.4f\n",
                                   ep, sqrt(best_val))
                break
            end
        end
        if verbose && (ep == 1 || ep % 20 == 0)
            @printf("[Forecaster] epoch %3d/%d  val_RMSE=%.4f\n", ep, epochs, sqrt(val_loss))
        end
    end

    Flux.loadmodel!(model, best_state)
    verbose && @printf("[Forecaster] done  best_val_RMSE=%.4f\n", sqrt(best_val))

    # 4 — Conformal calibration (z-space max abs error on the FINAL forecast) --
    Ŷca = fwd(model, X[:, ca, :], FUT[:, ca])   # residual or level
    residual && (Ŷca = Ŷca .+ SN[:, ca])        # reconstruct level before scoring
    Yca = Y[:, ca]
    cal_scores = vec(maximum(abs.(Ŷca .- Yca); dims=1))
    sort!(cal_scores)

    return LoadForecaster(model, μ, σ, cal_scores, horizon, window,
                          size(feat, 1), start_dow, residual, has_temp, μt, σt)
end

# ─────────────────────────────────────────────────────────────────────────────
#  Prediction with conformal intervals and scenario sampling
# ─────────────────────────────────────────────────────────────────────────────

"""
    predict_scenarios(fc, history; n_scenarios, α, hist_start_dow, seed) → NamedTuple

Given at least `fc.window` recent hourly observations (`history`, original units),
returns forecasts for the next `fc.horizon` hours.

`hist_start_dow` is the day-of-week (0=Mon … 6=Sun) of the first hour of the
window actually used (`history[end-window+1]`); pass it for real-date alignment.

Returns a NamedTuple with:
  mean      :: Vector{Float64} (horizon)    — LSTM point forecast
  lower     :: Vector{Float64} (horizon)    — conformal lower bound (1-α coverage)
  upper     :: Vector{Float64} (horizon)    — conformal upper bound
  scenarios :: Matrix{Float64} (horizon × n_scenarios) — sampled load profiles
  rmse_cal  :: Float64                      — calibration set RMSE (z-space)

Conformal intervals guarantee ≥ (1-α) marginal coverage on i.i.d. data.
Scenarios are Gaussian samples clipped to [lower, upper].
"""
function predict_scenarios(fc::LoadForecaster, history::AbstractVector{Float64};
        n_scenarios    :: Int     = 7,
        α              :: Float64 = 0.10,
        hist_start_dow :: Int     = 0,
        temp_history   :: Union{Nothing,AbstractVector{Float64}} = nothing,
        temp_future    :: Union{Nothing,AbstractVector{Float64}} = nothing,
        seed           :: Int     = 0)

    length(history) < fc.window &&
        error("history must have ≥ $(fc.window) values; got $(length(history))")
    if fc.has_temp
        (temp_history === nothing || temp_future === nothing) &&
            error("this forecaster needs temp_history (≥window) and temp_future (horizon)")
        length(temp_history) < fc.window && error("temp_history must have ≥ $(fc.window) values")
        length(temp_future) != fc.horizon && error("temp_future must have $(fc.horizon) values")
    end

    rng  = MersenneTwister(seed)
    hwin = collect(Float64, history[end-fc.window+1:end])
    z    = _znorm(hwin, fc.μ, fc.σ)
    cal  = _calendar_matrix(fc.window, 0, hist_start_dow)

    X = Array{Float32,3}(undef, fc.n_features, 1, fc.window)
    @views X[1, 1, :]   .= z
    @views X[2:6, 1, :] .= cal
    if fc.has_temp
        twin = _znorm(collect(Float64, temp_history[end-fc.window+1:end]), fc.μt, fc.σt)
        @views X[7, 1, :] .= twin
    end

    if fc.has_temp
        futz = reshape(_znorm(collect(Float64, temp_future), fc.μt, fc.σt), :, 1)
        ẑ = vec(fc.model(X, futz))            # residual or level (z-space)
    else
        ẑ = vec(fc.model(X))
    end

    # Reconstruct the level: add back the seasonal-naive baseline (load one week
    # earlier = first `horizon` hours of the z-normalised window).
    if fc.residual
        L = _seasonal_lag(fc.window)
        @views ẑ = ẑ .+ z[fc.window-L+1 : fc.window-L+fc.horizon]
    end

    # Conformal quantile at level ceil((1-α)(n+1)) / n
    n_cal = length(fc.cal_scores)
    q_idx = clamp(ceil(Int, (1 - α) * (n_cal + 1)), 1, n_cal)
    q     = fc.cal_scores[q_idx]

    ŷ     = _zdenorm(ẑ,        fc.μ, fc.σ)
    lower = max.(_zdenorm(ẑ .- q, fc.μ, fc.σ), 0.0)
    upper = max.(_zdenorm(ẑ .+ q, fc.μ, fc.σ), 0.0)

    half_band = (upper .- lower) ./ 4.0
    scen = Matrix{Float64}(undef, fc.horizon, n_scenarios)
    for s in 1:n_scenarios
        scen[:, s] = clamp.(ŷ .+ randn(rng, fc.horizon) .* half_band, lower, upper)
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
