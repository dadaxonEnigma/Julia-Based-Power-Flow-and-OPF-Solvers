"""
download_opsd.jl — fetch real load data for the LSTM forecaster.

Downloads the Open Power System Data (OPSD) hourly `time_series` package and
extracts a single country's actual-load column into a compact CSV that
`load_real_data` reads.

  Output:  data/opsd_de_load.csv   (timestamp, load_MW)   — committed (~1.4 MB)
  Cache :  data/raw/opsd_time_series.csv  (~125 MB)       — git-ignored

Usage:
  julia julia/scripts/download_opsd.jl            # Germany (DE), default
  julia julia/scripts/download_opsd.jl FR         # France, etc.

Citation (use in the thesis, IEEE style):
  Open Power System Data, "Time series" (2020-10-06).
  https://doi.org/10.25832/time_series/2020-10-06
"""

using CSV
using DataFrames
using Downloads

const OPSD_URL = "https://data.open-power-system-data.org/time_series/2020-10-06/time_series_60min_singleindex.csv"

country = isempty(ARGS) ? "DE" : uppercase(ARGS[1])
col     = "$(country)_load_actual_entsoe_transparency"

data_dir = normpath(joinpath(@__DIR__, "..", "..", "data"))
raw_dir  = joinpath(data_dir, "raw")
raw_path = joinpath(raw_dir, "opsd_time_series.csv")
out_path = joinpath(data_dir, "opsd_$(lowercase(country))_load.csv")
isdir(raw_dir) || mkpath(raw_dir)

# 1 — Download (skip if cached) ------------------------------------------------
if isfile(raw_path) && filesize(raw_path) > 1_000_000
    println("[download_opsd] using cached $raw_path ($(round(filesize(raw_path)/1e6,digits=1)) MB)")
else
    println("[download_opsd] downloading OPSD time_series (~125 MB)...")
    Downloads.download(OPSD_URL, raw_path)
    println("[download_opsd] saved $raw_path")
end

# 2 — Extract single country column -------------------------------------------
println("[download_opsd] extracting column $col ...")
df = CSV.read(raw_path, DataFrame; select=["utc_timestamp", col])
rename!(df, col => "load_MW")
n0 = nrow(df)
dropmissing!(df)
CSV.write(out_path, df)

println("[download_opsd] wrote $out_path")
println("[download_opsd]   rows: $(nrow(df)) (dropped $(n0 - nrow(df)) missing)")
println("[download_opsd]   range: $(df.utc_timestamp[1]) → $(df.utc_timestamp[end])")
