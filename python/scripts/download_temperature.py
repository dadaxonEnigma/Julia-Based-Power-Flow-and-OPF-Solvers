"""
download_temperature.py — population-weighted hourly temperature for Germany.

Fetches 2-m air temperature from the Open-Meteo ERA5 archive (free, no API key)
for the major German metro areas, weights them by population, and writes a single
hourly series aligned to the OPSD load period (2015-01-01 .. 2020-09-30).

  Output: data/opsd_de_temp.csv  (utc_timestamp, temp_C)  — committed (~1 MB)

Population-weighted national temperature is the standard exogenous driver for
load forecasting (heating/cooling demand). Source:
  Hersbach et al., ERA5 (Copernicus C3S); served via https://open-meteo.com
"""
import csv
import json
import os
import urllib.parse
import urllib.request

START, END = "2015-01-01", "2020-09-30"

# (name, lat, lon, metro population [millions]) — weights for the average.
CITIES = [
    ("Berlin",    52.52, 13.40, 3.7),
    ("Hamburg",   53.55,  9.99, 1.9),
    ("Munich",    48.14, 11.58, 1.5),
    ("Cologne",   50.94,  6.96, 1.1),
    ("Frankfurt", 50.11,  8.68, 0.76),
    ("Stuttgart", 48.78,  9.18, 0.63),
]

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "opsd_de_temp.csv")


def fetch(lat, lon):
    q = urllib.parse.urlencode({
        "latitude": lat, "longitude": lon,
        "start_date": START, "end_date": END,
        "hourly": "temperature_2m", "timezone": "UTC",
    })
    url = f"https://archive-api.open-meteo.com/v1/archive?{q}"
    with urllib.request.urlopen(url, timeout=120) as r:
        h = json.load(r)["hourly"]
    return h["time"], h["temperature_2m"]


print(f"[temp] fetching {len(CITIES)} cities, {START}..{END} ...")
times = None
weighted_sum = None
wsum = 0.0
for name, lat, lon, pop in CITIES:
    t, temp = fetch(lat, lon)
    if times is None:
        times = t
        weighted_sum = [0.0] * len(t)
    assert len(t) == len(times), f"{name}: length mismatch"
    for i, v in enumerate(temp):
        if v is not None:
            weighted_sum[i] += v * pop
    wsum += pop
    print(f"[temp]   {name:<10} {len(t)} hours  (pop weight {pop})")

avg = [round(s / wsum, 2) for s in weighted_sum]

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["utc_timestamp", "temp_C"])
    for ts, v in zip(times, avg):
        w.writerow([ts, v])

print(f"[temp] wrote {OUT}")
print(f"[temp]   {len(times)} hours  ({times[0]} .. {times[-1]})")
print(f"[temp]   temp range: min={min(avg):.1f}  max={max(avg):.1f}  C")
