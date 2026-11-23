# =====================================================================================
# 01_merge_exposure.R  -  attach daily PM2.5 / AQI / weather exposure to every lab report
# Input : claude_cleaned_tables/output/analysis_ready_data.csv + AQI/okara_consensus_aqi_2026.csv (+ weather)
# Output: analysis/results/data/analysis_exposure.csv, analysis/results/data/aqi_daily_processed.csv
#         Fig 2 (timeline), Fig S1 (seasonality by year), table of monthly air-quality summary
# =====================================================================================
source("00_setup.R")

a   <- read_csv(DATA_CSV, show_col_types = FALSE, guess_max = 20000) |> mutate(Est_Date = as.Date(Est_Date))
aqi <- read_csv(AQI_CSV,  show_col_types = FALSE) |> mutate(date = as.Date(date)) |> arrange(date)
stopifnot(all(c("okara_pm2_5", "okara_us_aqi", "consensus_okara_pm2_5") %in% names(aqi)))

# ---- complete daily calendar, short-gap interpolation (<= 3 days) -------------------
day <- tibble(date = seq(min(aqi$date), max(aqi$date), by = "day")) |> left_join(aqi, by = "date")
interp <- function(x, maxgap = 3) {
  if (all(is.na(x))) return(x)
  idx <- seq_along(x); ok <- !is.na(x); y <- approx(idx[ok], x[ok], xout = idx, rule = 1)$y
  r <- rle(is.na(x)); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  for (i in which(r$values & r$lengths > maxgap)) y[starts[i]:ends[i]] <- NA
  y
}
n_missing_before <- sum(is.na(day$okara_pm2_5))
day <- day |> mutate(across(c(okara_pm2_5, okara_pm10, okara_nitrogen_dioxide, okara_sulphur_dioxide,
                              okara_ozone, okara_carbon_monoxide, okara_us_aqi, consensus_okara_pm2_5), interp))
message("AQI days: ", nrow(day), " | PM2.5 missing before/after interpolation: ", n_missing_before, " / ", sum(is.na(day$okara_pm2_5)))

roll_right <- function(x, k) as.numeric(stats::filter(x, rep(1 / k, k), sides = 1))   # mean of day t-k+1 ... t

day <- day |> mutate(
  pm25 = okara_pm2_5, pm10 = okara_pm10, no2 = okara_nitrogen_dioxide, so2 = okara_sulphur_dioxide,
  o3 = okara_ozone, co = okara_carbon_monoxide, aqi = okara_us_aqi, cons_pm25 = consensus_okara_pm2_5,
  # PM2.5 rolling windows & lags
  pm25_L0 = pm25, pm25_L1 = lag(pm25, 1),
  pm25_ma3 = roll_right(pm25, 3), pm25_ma7 = roll_right(pm25, 7),
  pm25_ma14 = roll_right(pm25, 14), pm25_ma28 = roll_right(pm25, 28),
  # Co-pollutant rolling windows (lag 0, lag 1, 3d, 7d, 14d, 28d)
  pm10_L0 = pm10, pm10_L1 = lag(pm10, 1), pm10_ma3 = roll_right(pm10, 3), pm10_ma7 = roll_right(pm10, 7), pm10_ma14 = roll_right(pm10, 14), pm10_ma28 = roll_right(pm10, 28),
  no2_L0  = no2,  no2_L1  = lag(no2, 1),  no2_ma3  = roll_right(no2, 3),  no2_ma7  = roll_right(no2, 7),  no2_ma14  = roll_right(no2, 14),  no2_ma28  = roll_right(no2, 28),
  so2_L0  = so2,  so2_L1  = lag(so2, 1),  so2_ma3  = roll_right(so2, 3),  so2_ma7  = roll_right(so2, 7),  so2_ma14  = roll_right(so2, 14),  so2_ma28  = roll_right(so2, 28),
  o3_L0   = o3,   o3_L1   = lag(o3, 1),   o3_ma3   = roll_right(o3, 3),   o3_ma7   = roll_right(o3, 7),   o3_ma14   = roll_right(o3, 14),   o3_ma28   = roll_right(o3, 28),
  co_L0   = co,   co_L1   = lag(co, 1),   co_ma3   = roll_right(co, 3),   co_ma7   = roll_right(co, 7),   co_ma14   = roll_right(co, 14),   co_ma28   = roll_right(co, 28),
  # AQI & consensus PM2.5
  aqi_ma7 = roll_right(aqi, 7), cons_pm25_ma7 = roll_right(cons_pm25, 7), cons_pm25_ma14 = roll_right(cons_pm25, 14),
  days_aqi200_14 = roll_right(as.numeric(aqi > 200), 14) * 14,
  pm25_lead7 = lead(pm25_ma7, 7)          # mean of the 7 days AFTER the test: negative-control exposure
)
for (k in 0:LAG_MAX) day[[paste0("pm25_L", k)]] <- lag(day$pm25, k)

# ---- weather (loads existing file OR fetches Open-Meteo archive directly into DAT_DIR) ----
wx_local <- file.path(DAT_DIR, "okara_weather_daily.csv")
wx_source <- if (file.exists(WEATHER_CSV)) WEATHER_CSV else if (file.exists(wx_local)) wx_local else NA_character_

if (is.na(wx_source)) {
  message("Fetching Open-Meteo Archive weather for Okara (30.8081N, 73.4534E)...")
  wx_url <- sprintf(
    "https://archive-api.open-meteo.com/v1/archive?latitude=30.8081&longitude=73.4534&start_date=%s&end_date=%s&daily=temperature_2m_mean,temperature_2m_max,temperature_2m_min,relative_humidity_2m_mean&timezone=Asia%%2FKarachi&format=csv",
    min(day$date), max(day$date)
  )
  tryCatch({
    wx_raw <- read_csv(wx_url, skip = 3, show_col_types = FALSE) |>
      transmute(
        date = as.Date(time),
        temp_mean = `temperature_2m_mean (°C)`,
        temp_max  = `temperature_2m_max (°C)`,
        temp_min  = `temperature_2m_min (°C)`,
        rh_mean   = `relative_humidity_2m_mean (%)`
      )
    write_csv(wx_raw, wx_local)
    wx_source <- wx_local
    message("Saved fetched weather to: ", wx_local)
  }, error = function(e) warning("Open-Meteo fetch failed: ", conditionMessage(e)))
}

if (!is.na(wx_source) && file.exists(wx_source)) {
  wx <- read_csv(wx_source, show_col_types = FALSE) |> mutate(date = as.Date(date)) |> arrange(date)
  if (!"temp_mean" %in% names(wx) && "temperature_2m_mean" %in% names(wx)) wx$temp_mean <- wx$temperature_2m_mean
  if (!"rh_mean"   %in% names(wx) && "relative_humidity_2m_mean" %in% names(wx)) wx$rh_mean <- wx$relative_humidity_2m_mean
  wx <- wx |> mutate(
    temp_L0 = temp_mean, temp_L1 = lag(temp_mean, 1),
    temp_ma3 = roll_right(temp_mean, 3), temp_ma7 = roll_right(temp_mean, 7),
    temp_ma14 = roll_right(temp_mean, 14), temp_ma28 = roll_right(temp_mean, 28),
    rh_L0 = rh_mean, rh_L1 = lag(rh_mean, 1),
    rh_ma3 = roll_right(rh_mean, 3), rh_ma7 = roll_right(rh_mean, 7),
    rh_ma14 = roll_right(rh_mean, 14), rh_ma28 = roll_right(rh_mean, 28)
  )
  day <- day |> left_join(wx |> select(date, starts_with("temp_L"), starts_with("temp_ma"), starts_with("rh_L"), starts_with("rh_ma")), by = "date")
  message("Weather merged (temperature + humidity across all lags) from: ", wx_source)
} else {
  warning("No weather file available -> temperature/humidity-adjusted models will be SKIPPED.")
}
write_csv(day, file.path(DAT_DIR, "aqi_daily_processed.csv"))

# ---- merge onto lab reports ----------------------------------------------------------
m <- a |> left_join(day |> select(-any_of(c("okara_us_aqi", "okara_european_aqi", "okara_pm2_5", "okara_pm10",
        "okara_nitrogen_dioxide", "okara_sulphur_dioxide", "okara_ozone", "okara_carbon_monoxide",
        "openaq_okara_pm2_5", "openaq_okara_us_aqi", "consensus_okara_pm2_5", "consensus_okara_us_aqi"))),
        by = c("Est_Date" = "date"))
miss <- sum(is.na(m$pm25_ma7)); message("Reports without 7-day PM2.5 exposure: ", miss, " of ", nrow(m))
write_csv(m, file.path(DAT_DIR, "analysis_exposure.csv"))

# ---- monthly air-quality summary table (for Methods / Results text) ----------------------
mon <- day |> mutate(month = floor_date(date, "month")) |> group_by(month) |>
  summarise(days = n(), pm25_mean = mean(pm25, na.rm = TRUE), pm25_median = median(pm25, na.rm = TRUE),
            pm25_max = max(pm25, na.rm = TRUE), pct_days_gt15 = 100 * mean(pm25 > 15, na.rm = TRUE),
            pct_days_gt35 = 100 * mean(pm25 > 35, na.rm = TRUE), pct_days_gt75 = 100 * mean(pm25 > 75, na.rm = TRUE),
            pct_days_gt150 = 100 * mean(pm25 > 150, na.rm = TRUE), aqi_mean = mean(aqi, na.rm = TRUE), .groups = "drop") |>
  mutate(across(where(is.numeric), ~ round(.x, 1)))
write_tab(mon, "aqi_monthly_summary")

# ---- Figure 2: PM2.5 timeline + reports per month -------------------------------------
study <- range(a$Est_Date, na.rm = TRUE)
smog_B <- as.Date(c("2025-11-15", "2026-02-15"))
rep_m <- a |> mutate(mo = floor_date(Est_Date, "month") + 14) |> count(mo)
p_top <- ggplot(day, aes(date)) +
  annotate("rect", xmin = study[1], xmax = study[2], ymin = -Inf, ymax = Inf, fill = "grey90") +
  annotate("rect", xmin = smog_B[1], xmax = smog_B[2], ymin = -Inf, ymax = Inf, fill = COL_SMOG, alpha = 0.12) +
  geom_line(aes(y = pm25), colour = "grey60", linewidth = 0.3) +
  geom_line(aes(y = pm25_ma7), colour = "black", linewidth = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 15, linetype = "dashed", colour = COL_WARM) +
  annotate("text", x = min(day$date), y = 15, label = "WHO 24-h guideline (15)", hjust = 0, vjust = -0.6,
           size = 2.8, colour = COL_WARM) +
  annotate("text", x = mean(smog_B), y = max(day$pm25, na.rm = TRUE) * 0.97, label = "Smog window\n(15 Nov - 15 Feb)",
           size = 2.8, colour = COL_SMOG, fontface = "bold") +
  annotate("text", x = study[2], y = max(day$pm25, na.rm = TRUE) * 0.97, label = "Lab-report window (grey)",
           size = 2.8, hjust = 1, colour = "grey30") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y", expand = expansion(mult = 0.01)) +
  labs(x = NULL, y = expression(PM[2.5]~"("*mu*g/m^3*")"), title = "Daily PM2.5 in Okara (grey = daily, black = 7-day mean)")
p_bot <- ggplot(rep_m, aes(mo, n)) + geom_col(width = 25, fill = "grey40") +
  geom_text(aes(label = n), vjust = -0.3, size = 2.8) +
  scale_x_date(limits = range(day$date), date_breaks = "2 months", date_labels = "%b\n%Y", expand = expansion(mult = 0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + labs(x = NULL, y = "Lab reports")
save_fig(p_top / p_bot + plot_layout(heights = c(3, 1)), "fig2_pm25_timeline", 180, 120, tiff = TRUE)

# ---- Figure S1: seasonal pattern by calendar year -----------------------------------
sea <- day |> mutate(year = factor(year(date)), mon = month(date, label = TRUE)) |>
  group_by(year, mon) |> summarise(pm25 = mean(pm25, na.rm = TRUE), n = n(), .groups = "drop")
pS1 <- ggplot(sea, aes(mon, pm25, colour = year, group = year)) + geom_line(linewidth = 0.8) + geom_point() +
  scale_colour_manual(values = c("2025" = COL_WARM, "2026" = COL_SMOG)) +
  geom_hline(yintercept = 15, linetype = "dashed") +
  labs(x = NULL, y = expression("Monthly mean PM"[2.5]~"("*mu*g/m^3*")"), colour = "Year",
       title = "Seasonality of PM2.5 in Okara (November 2026 is partial)")
save_fig(pS1, "figS1_pm25_seasonality", 120, 80)
message("01 done.")