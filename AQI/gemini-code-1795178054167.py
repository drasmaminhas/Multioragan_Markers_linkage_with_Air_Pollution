cat << 'EOF' > cross_validate_openaq.py
import os
import requests
import pandas as pd
import numpy as np
import sqlite3

API_KEY = os.getenv("OPENAQ_API_KEY")
if not API_KEY:
    raise SystemExit("Error: Set OPENAQ_API_KEY first -> export OPENAQ_API_KEY='your_key'")

HEADERS = {"X-API-Key": API_KEY}
BASE_URL = "https://api.openaq.org/v3"

OKARA_COORDS = "30.8081,73.4534"
START_DATE = "2025-01-01"
END_DATE = "2026-11-20"

def pm25_to_us_aqi(pm):
    if pd.isna(pm) or pm < 0:
        return np.nan
    bps = [
        (0.0, 9.0, 0, 50),
        (9.1, 35.4, 51, 100),
        (35.5, 55.4, 101, 150),
        (55.5, 125.4, 151, 200),
        (125.5, 225.4, 201, 300),
        (225.5, 500.4, 301, 500),
    ]
    pm_r = round(float(pm), 1)
    for c_low, c_high, i_low, i_high in bps:
        if c_low <= pm_r <= c_high:
            return round(((i_high - i_low) / (c_high - c_low)) * (pm_r - c_low) + i_low)
    return 500 if pm_r > 500.4 else np.nan

# 1. Query only Okara's local OpenAQ ground locations (25 km radius)
print("Fetching OpenAQ v3 ground locations for Okara...")
loc_resp = requests.get(
    f"{BASE_URL}/locations",
    headers=HEADERS,
    params={"coordinates": OKARA_COORDS, "radius": 25000, "limit": 100},
    timeout=30,
)
loc_resp.raise_for_status()
locations = loc_resp.json().get("results", [])
print(f"Found {len(locations)} Okara stations:")
for loc in locations:
    print(f"  - {loc.get('name')} (ID: {loc.get('id')})")

records = []
for loc in locations:
    loc_name = loc.get("name", "Unknown")
    for sensor in loc.get("sensors", []):
        param_name = sensor.get("parameter", {}).get("name")
        if param_name not in ["pm25", "pm10", "no2", "so2", "o3", "co"]:
            continue
        sid = sensor["id"]
        day_resp = requests.get(
            f"{BASE_URL}/sensors/{sid}/days",
            headers=HEADERS,
            params={"date_from": START_DATE, "date_to": END_DATE, "limit": 250},
            timeout=30,
        )
        if day_resp.status_code != 200:
            continue
        for item in day_resp.json().get("results", []):
            dt_str = item.get("period", {}).get("datetimeFrom", {}).get("local", "")[:10]
            val = item.get("value")
            if dt_str and val is not None:
                records.append({
                    "date": pd.to_datetime(dt_str).date(),
                    "station": loc_name,
                    "parameter": param_name,
                    "value": float(val),
                })

if records:
    raw_df = pd.DataFrame(records)
    openaq_df = (
        raw_df.groupby(["date", "parameter"])["value"]
        .mean()
        .unstack("parameter")
        .reset_index()
    )
    rename_map = {
        "pm25": "openaq_okara_pm2_5",
        "pm10": "openaq_okara_pm10",
        "no2": "openaq_okara_no2",
        "so2": "openaq_okara_so2",
        "o3": "openaq_okara_o3",
        "co": "openaq_okara_co",
    }
    openaq_df.rename(columns=rename_map, inplace=True)
    if "openaq_okara_pm2_5" in openaq_df.columns:
        openaq_df["openaq_okara_us_aqi"] = openaq_df["openaq_okara_pm2_5"].apply(pm25_to_us_aqi)
else:
    openaq_df = pd.DataFrame(columns=["date", "openaq_okara_pm2_5", "openaq_okara_us_aqi"])

# 2. Load Source 1 (Open-Meteo CAMS) and keep only Okara columns
cams_df = pd.read_csv("okara_and_regional_aqi_2026.csv")
cams_df["date"] = pd.to_datetime(cams_df["date"]).dt.date
okara_cols = ["date"] + [c for c in cams_df.columns if c.startswith("okara_")]
cams_okara = cams_df[okara_cols].copy()

# 3. Merge Source 1 (Open-Meteo CAMS) + Source 2 (OpenAQ Okara Ground Sensors)
merged = cams_okara.merge(openaq_df, on="date", how="left")

# Consensus mean across both Okara sources
merged["consensus_okara_pm2_5"] = merged[["okara_pm2_5", "openaq_okara_pm2_5"]].mean(axis=1)
merged["consensus_okara_us_aqi"] = merged[["okara_us_aqi", "openaq_okara_us_aqi"]].mean(axis=1)

# 4. Cross-validation summary
valid_both = merged.dropna(subset=["okara_pm2_5", "openaq_okara_pm2_5"])
print("\n=== OKARA TWO-SOURCE CONFIRMATION (Jan–Jul 2026) ===")
print(f"Total days in window:                 {len(merged)}")
print(f"Days with OpenAQ ground confirmation: {len(valid_both)}")
if len(valid_both) >= 3:
    corr = valid_both["okara_pm2_5"].corr(valid_both["openaq_okara_pm2_5"])
    mae = (valid_both["okara_pm2_5"] - valid_both["openaq_okara_pm2_5"]).abs().mean()
    print(f"Source 1 (Open-Meteo CAMS) PM2.5:     {valid_both['okara_pm2_5'].mean():.2f} µg/m³")
    print(f"Source 2 (OpenAQ Ground) PM2.5:       {valid_both['openaq_okara_pm2_5'].mean():.2f} µg/m³")
    print(f"Two-Source Consensus PM2.5:           {valid_both['consensus_okara_pm2_5'].mean():.2f} µg/m³")
    print(f"Pearson r (Source 1 vs Source 2):     {corr:.3f} (MAE = {mae:.2f} µg/m³)")
elif len(valid_both) == 0:
    print("Note: The 3 OpenAQ stations in Okara had 0 active daily records between Jan–Jul 2026.")

merged.to_csv("okara_consensus_aqi_2026.csv", index=False)
conn = sqlite3.connect("gynae_env_exposure.db")
merged.to_sql("consensus_air_quality", conn, if_exists="replace", index=False)
conn.close()
print("\nSaved okara_consensus_aqi_2026.csv and SQLite table 'consensus_air_quality'.")
EOF
python cross_validate_openaq.py