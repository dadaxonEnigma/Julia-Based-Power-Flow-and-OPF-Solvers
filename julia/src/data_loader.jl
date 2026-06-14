"""
    data_loader.jl

Loader for **real** historical load data, replacing `generate_synthetic_data`
as the training source for the LSTM forecaster.

Data source: Open Power System Data (OPSD), `time_series` package, which
re-publishes hourly actual-load series from the ENTSO-E Transparency Platform.
Free, no API key, citable.

  Müller, J. et al. "Open Power System Data — Time series."
  https://doi.org/10.25832/time_series/2020-10-06

Pipeline (drop-in replacement for the synthetic generator):

  data = load_real_data()                    # (n_days × 24) matrix, pu
  fc   = train_forecaster(data)              # same call as before
  pred = predict_scenarios(fc, history)

The returned matrix has the same `(n_days × 24)` layout and `pu` scale
(peak-normalised, so max ≈ 1.0) as `generate_synthetic_data`, so every
downstream function (`train_forecaster`, `predict_scenarios`, `lopf_stochastic`)
works unchanged.
"""

using CSV
using DataFrames
using Dates
using Statistics

# Default location of the cleaned single-column CSV (timestamp, load_MW).
# Produced by extracting one country column from the full OPSD download.
const DEFAULT_LOAD_CSV = normpath(joinpath(@__DIR__, "..", "..", "data", "opsd_de_load.csv"))
const DEFAULT_TEMP_CSV = normpath(joinpath(@__DIR__, "..", "..", "data", "opsd_de_temp.csv"))

# ─────────────────────────────────────────────────────────────────────────────
#  Helpers
# ─────────────────────────────────────────────────────────────────────────────

# Parse OPSD UTC timestamps, accepting either pre-parsed DateTime (CSV.jl may
# auto-detect) or ISO-8601 strings like "2015-01-01T00:00:00Z".
_to_datetime(x::DateTime) = x
_to_datetime(x::AbstractString) = DateTime(replace(String(x), "Z" => ""), dateformat"yyyy-mm-ddTHH:MM:SS")

# Group an hourly (utc_timestamp, value) CSV into Date → 24-element vector
# (missing for absent hours). Shared by the load and temperature loaders.
function _group_by_day(path::AbstractString, col::AbstractString)
    df = CSV.read(path, DataFrame; select=["utc_timestamp", col])
    rename!(df, col => :v)
    dropmissing!(df)
    dts  = _to_datetime.(df.utc_timestamp)
    days = Date.(dts); hrs = hour.(dts); vals = Float64.(df.v)
    by = Dict{Date, Vector{Union{Missing,Float64}}}()
    for i in eachindex(vals)
        slot = get!(by, days[i], Vector{Union{Missing,Float64}}(missing, 24))
        slot[hrs[i] + 1] = vals[i]
    end
    return by
end

_complete_days(by) = sort!([d for (d, v) in by if !any(ismissing, v)])

# ─────────────────────────────────────────────────────────────────────────────
#  Main loader
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_real_data(path = DEFAULT_LOAD_CSV; column, normalize, max_days, verbose)
                                                        → NamedTuple

Read a cleaned hourly-load CSV and reshape it into a `(n_days × 24)` matrix in
per-unit, ready for `train_forecaster`.

Arguments:
  path      :: String   — CSV with a timestamp column and a numeric load column

Keyword arguments:
  column    :: String   — name of the load column (default "load_MW")
  normalize :: Symbol   — `:peak` (÷ max → max≈1.0, matches DEFAULT_LOAD_PROFILE),
                          `:mean` (÷ mean → mean≈1.0), or `:none` (raw MW)
  max_days  :: Int      — keep only the most recent `max_days` complete days
                          (0 = keep all; useful to cap LSTM training time)
  verbose   :: Bool     — print a short data summary

Only **complete** days (24 non-missing hours) are kept, so the hour-of-day
alignment is exact. Returns a NamedTuple:

  data      :: Matrix{Float64}  — (n_days × 24), per-unit (unless :none)
  n_days    :: Int
  peak_MW   :: Float64          — scaling constant used for :peak
  mean_MW   :: Float64          — scaling constant used for :mean
  normalize :: Symbol
  first_day :: Date
  last_day  :: Date
  source    :: String

PyPSA equivalent: building the `loads_t.p_set` time series from a CSV before
calling `network.optimize()`.
"""
function load_real_data(path::AbstractString = DEFAULT_LOAD_CSV;
        column    :: AbstractString = "load_MW",
        normalize :: Symbol  = :peak,
        max_days  :: Int     = 0,
        temp_path :: Union{Nothing,AbstractString} = nothing,
        verbose   :: Bool    = true)

    isfile(path) || error("Load data not found at $path.\n" *
        "Run `julia julia/scripts/download_opsd.jl` first to fetch and extract it.")
    normalize in (:peak, :mean, :none) ||
        error("normalize must be :peak, :mean or :none (got :$normalize)")

    by_load  = _group_by_day(path, column)
    complete = _complete_days(by_load)
    isempty(complete) && error("No complete 24-hour days found in $path.")

    # Optional temperature: keep only days complete in BOTH series ------------
    has_temp = temp_path !== nothing
    by_temp  = nothing
    if has_temp
        isfile(temp_path) || error("Temperature data not found at $temp_path.\n" *
            "Run `python python/download_temperature.py` first.")
        by_temp  = _group_by_day(temp_path, "temp_C")
        complete = sort!(intersect(complete, _complete_days(by_temp)))
        isempty(complete) && error("No days with both load and temperature.")
    end

    if max_days > 0 && length(complete) > max_days
        complete = complete[end-max_days+1:end]   # keep most recent
    end

    n_days = length(complete)
    data   = Matrix{Float64}(undef, n_days, 24)
    for (r, d) in enumerate(complete)
        data[r, :] = Float64.(by_load[d])
    end

    temp = nothing
    if has_temp
        temp = Matrix{Float64}(undef, n_days, 24)
        for (r, d) in enumerate(complete)
            temp[r, :] = Float64.(by_temp[d])
        end
    end

    peak_MW = maximum(data)
    mean_MW = mean(data)
    if normalize === :peak
        data ./= peak_MW
    elseif normalize === :mean
        data ./= mean_MW
    end

    if verbose
        println("[data_loader] $(basename(path)):  $n_days complete days  " *
                "($(complete[1]) → $(complete[end]))")
        println("[data_loader] raw load  peak=$(round(peak_MW,digits=1)) MW  " *
                "mean=$(round(mean_MW,digits=1)) MW  |  normalize=:$normalize" *
                (has_temp ? "  |  +temperature" : ""))
    end

    return (data=data, temp=temp, n_days=n_days, peak_MW=peak_MW, mean_MW=mean_MW,
            normalize=normalize, first_day=complete[1], last_day=complete[end],
            source=basename(path))
end

# ─────────────────────────────────────────────────────────────────────────────
#  Honest train / test split
# ─────────────────────────────────────────────────────────────────────────────

"""
    train_test_split_days(data; test_days) → (train, test)

Chronological split of a `(n_days × 24)` matrix into an earlier training block
and a held-out, **later** test block of `test_days` rows.

Unlike the random split used internally for validation, this is a strict
out-of-time hold-out: the LSTM never sees the test period during training,
which is the honest way to report forecast metrics in the thesis.
"""
function train_test_split_days(data::Matrix{Float64}; test_days::Int = 30)
    n = size(data, 1)
    test_days = clamp(test_days, 1, n - 1)
    train = data[1:n-test_days, :]
    test  = data[n-test_days+1:end, :]
    return (train=train, test=test)
end
