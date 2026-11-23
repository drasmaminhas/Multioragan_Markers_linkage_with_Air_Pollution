# =====================================================================================
# 07_sensitivity.R  -  does the conclusion survive reasonable alternative choices?
#   1. Scenario forest : exclusions + 3 smog definitions + season-boundary grid (smog vs warm, age/sex adjusted)
#   2. Date-jitter     : estimated test dates are only +/- a few days accurate -> re-run 200x with random date error
#   3. Exposure source : CAMS vs two-source consensus PM2.5, 7-day vs 14-day, PM10
#   4. Negative controls: (a) hepatitis/HIV positivity (b) FUTURE 7-day PM2.5 as a fake exposure
#   5. Multiple imputation (mice) for the CBC outcomes
#   6. E-values for the smog-vs-warm estimates
# =====================================================================================
source("00_setup.R")
d   <- load_analysis()
day <- read_csv(file.path(DAT_DIR, "aqi_daily_processed.csv"), show_col_types = FALSE) |> mutate(date = as.Date(date))
covs <- c("Age_Group", "Sex")
SENS <- MARKERS |> filter(var %in% c("NLR", "SII", "CBC_WBC_TLC_10e3_uL", "CBC_Platelet_Count_10e3_uL", "CBC_Hemoglobin_HB_g_dL",
                                     "eGFR_CKDEPI2021", "Blood_Glucose_Random_mg_dL"))
smog_by_dates <- function(dt, start, end, warm = as.Date("2026-03-01")) case_when(dt >= warm ~ 0L, dt >= start & dt <= end ~ 1L, TRUE ~ NA_integer_)
run_scn <- function(dat, x, scen) pmap_dfr(SENS, function(var, label, system, tr, primary) {
  r <- fit_cont(dat, var, tr, x, covs); if (is.null(r)) NULL else mutate(r, label = label, scenario = scen)
})

# ------------------------------------------------------------------ 1. scenario forest
d$smog_primary <- d$smog
scn <- list(
  "Primary (15 Nov-15 Feb vs warm)" = list(dat = d, x = "smog_primary"),
  "Definition A: whole months Nov-Feb" = list(dat = d |> mutate(s = Smog_A_WholeMonths_NovFeb), x = "s"),
  "Definition C: 15 Nov-31 Jan" = list(dat = d |> mutate(s = Smog_C_Nov15_Jan31), x = "s"),
  "Exclude CBC-suspect rows" = list(dat = d |> filter(CBC_Suspect %in% c(0, NA)), x = "smog_primary")
)
if ("Possible_Repeat_Patient" %in% names(d)) {
  scn[["Exclude possible repeat patients"]] <- list(dat = d |> filter(Possible_Repeat_Patient == 0), x = "smog_primary")
}
if ("Est_Date_Confidence" %in% names(d)) {
  scn[["High-confidence dates only (excl. Mar & Jul)"]] <- list(dat = d |> filter(Est_Date_Confidence == "high"), x = "smog_primary")
}
scn[["Exclude age < 40"]] <- list(dat = d |> filter(Age_Years >= 40), x = "smog_primary")

grid <- expand_grid(start = as.Date(c("2025-11-01", "2025-11-29")), end = as.Date(c("2026-01-31", "2026-02-28")))
for (i in seq_len(nrow(grid)))
  scn[[sprintf("Boundary: %s to %s", format(grid$start[i], "%d %b"), format(grid$end[i], "%d %b"))]] <-
    list(dat = d |> mutate(s = smog_by_dates(Est_Date, grid$start[i], grid$end[i])), x = "s")
sres <- imap_dfr(scn, function(z, nm) run_scn(z$dat, z$x, nm))
write_tab(sres |> select(scenario, label, n, est, lo, hi, std, std_lo, std_hi, p), "tab7a_scenarios")
sres$scenario <- factor(sres$scenario, levels = rev(names(scn))); sres$label <- factor(sres$label, levels = SENS$label)
p_s <- ggplot(sres, aes(std, scenario, colour = scenario == names(scn)[1])) + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = std_lo, xmax = std_hi), size = 0.3) + scale_colour_manual(values = c(`TRUE` = COL_SMOG, `FALSE` = "grey30"), guide = "none") +
  facet_wrap(~ label, ncol = 4) + labs(x = "Smog minus warm (SD of outcome, 95% CI; age group + sex adjusted)", y = NULL,
  title = "Robustness of the smog-vs-warm difference to analytic choices") + theme(axis.text.y = element_text(size = 6.5), axis.text.x = element_text(size = 7))
save_fig(p_s, "fig15_sensitivity_scenarios", 190, 130, tiff = TRUE)

# ------------------------------------------------------------------ 2. date-jitter simulation (PM2.5 within Nov-Feb)
pm7 <- setNames(day$pm25_ma7, as.character(day$date)); dnf <- d |> filter(Smog_A_WholeMonths_NovFeb == 1)
half <- if ("Est_Date_Confidence" %in% names(dnf)) if_else(dnf$Est_Date_Confidence == "high", 2L, 4L) else rep(3L, nrow(dnf))
base <- run_scn(dnf, "pm25_ma7_50", "observed") |> select(label, obs = std)
set.seed(SEED)
jit <- map_dfr(seq_len(N_SIM_JITTER), function(i) {
  dj <- dnf; shift <- sapply(half, function(h) sample(-h:h, 1)); nd <- dnf$Est_Date + shift
  nd <- pmin(pmax(nd, min(day$date) + 28), max(day$date)); dj$pm25_ma7_50 <- unname(pm7[as.character(nd)]) / PM_UNIT
  run_scn(dj, "pm25_ma7_50", "jitter") |> mutate(sim = i)
})
jsum <- jit |> group_by(label) |> summarise(jit_mean = mean(std), jit_lo = quantile(std, .025), jit_hi = quantile(std, .975), .groups = "drop") |> left_join(base, by = "label")
write_tab(jsum, "tab7b_date_jitter")
jit$label <- factor(jit$label, levels = SENS$label); base$label <- factor(base$label, levels = SENS$label)
p_j <- ggplot(jit, aes(std)) + geom_density(fill = "grey80", colour = "grey40") + geom_vline(data = base, aes(xintercept = obs), colour = COL_SMOG, linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") + facet_wrap(~ label, scales = "free", ncol = 4) +
  labs(x = "Effect of +50 ug/m3 7-day PM2.5 (SD of outcome), Nov-Feb only", y = "Density",
       title = sprintf("Date-error simulation (%d runs): grey = with random date error, red = as observed", N_SIM_JITTER))
save_fig(p_j, "figS5_date_jitter", 190, 100)

# ------------------------------------------------------------------ 3. exposure source / window (Nov-Feb)
ex <- c(pm25_ma7_50 = "PM2.5 CAMS 7-d", cons_pm25_ma7_50 = "PM2.5 consensus 7-d", pm25_ma14_50 = "PM2.5 CAMS 14-d", pm10_ma7_50 = "PM10 7-d")
ex <- ex[names(ex) %in% names(dnf)]
eres <- imap_dfr(ex, function(nm, x) run_scn(dnf, x, nm)); write_tab(eres |> select(scenario, label, n, std, std_lo, std_hi, p), "tab7c_exposure_sources")
eres$label <- factor(eres$label, levels = SENS$label)
p_e <- ggplot(eres, aes(std, scenario)) + geom_vline(xintercept = 0, linetype = "dashed") + geom_pointrange(aes(xmin = std_lo, xmax = std_hi), size = 0.3) +
  facet_wrap(~ label, ncol = 4) + labs(x = "Effect per +50 ug/m3 (SD of outcome), Nov-Feb only", y = NULL, title = "Exposure source and metric")
save_fig(p_e, "figS6_exposure_sources", 190, 90)

# ------------------------------------------------------------------ 4. negative controls
with_col <- function(r, ...) if (is.null(r)) NULL else mutate(r, ...)
nc_vars <- NEG_CONTROL[NEG_CONTROL %in% names(d)]
nco <- imap_dfr(nc_vars, function(lb, v) bind_rows(
  with_col(fit_bin(d |> filter(!is.na(smog)), v, "smog", covs), exposure = "Smog vs warm", outcome = lb),
  with_col(fit_bin(dnf, v, "pm25_ma7_50", covs), exposure = "+50 ug/m3 PM2.5 (Nov-Feb)", outcome = lb)))
write_tab(nco, "tab7d_negative_control_outcomes")
nce <- bind_rows(run_scn(dnf, "pm25_ma7_50", "Real exposure: preceding 7 days"), run_scn(dnf, "pm25_lead7_50", "Negative control: FOLLOWING 7 days")) |>
  mutate(label = factor(label, levels = SENS$label)); write_tab(nce |> select(scenario, label, n, std, std_lo, std_hi, p), "tab7e_negative_control_exposure")
p_nc <- ggplot(nce, aes(std, scenario, colour = scenario)) + geom_vline(xintercept = 0, linetype = "dashed") + geom_pointrange(aes(xmin = std_lo, xmax = std_hi), size = 0.3) +
  scale_colour_manual(values = c("grey40", COL_SMOG), guide = "none") + facet_wrap(~ label, ncol = 4) +
  labs(x = "Effect per +50 ug/m3 (SD of outcome), Nov-Feb", y = NULL, title = "Negative-control exposure: air quality AFTER the test should not 'cause' the result",
       caption = "Note: PM2.5 is autocorrelated, so a following-week effect can still appear if the preceding-week effect is real; compare magnitudes.")
save_fig(p_nc, "figS7_negative_control_exposure", 190, 90)

# ------------------------------------------------------------------ 5. multiple imputation (CBC outcomes)
if (requireNamespace("mice", quietly = TRUE)) {
  md <- d |> filter(!is.na(smog)) |> transmute(smog, Age_Years, Sex_Female, lNLR = log(NLR), Hb = CBC_Hemoglobin_HB_g_dL,
                                               Plt = CBC_Platelet_Count_10e3_uL, WBC = CBC_WBC_TLC_10e3_uL, lSII = log(SII))
  md <- md |> filter(is.finite(lNLR) | is.na(lNLR)); md$lNLR[!is.finite(md$lNLR)] <- NA; md$lSII[!is.finite(md$lSII)] <- NA
  set.seed(SEED); imp <- suppressWarnings(mice::mice(md, m = 20, maxit = 10, method = "pmm", printFlag = FALSE))
  mi <- map_dfr(c("lNLR", "lSII", "WBC", "Plt", "Hb"), function(y) {
    fits <- with(imp, lm(as.formula(paste(y, "~ smog + Age_Years + Sex_Female"))))
    s <- summary(mice::pool(fits), conf.int = TRUE); r <- s[s$term == "smog", ]
    sdy <- sd(md[[y]], na.rm = TRUE); tibble(outcome = y, est = r$estimate, lo = r$`2.5 %`, hi = r$`97.5 %`, p = r$p.value, std = r$estimate / sdy, std_lo = r$`2.5 %` / sdy, std_hi = r$`97.5 %` / sdy, n_imputed = nrow(md))
  })
  write_tab(mi, "tab7f_multiple_imputation")
} else message("  {mice} not installed: multiple-imputation sensitivity skipped.")

# ------------------------------------------------------------------ 6. E-values (smog vs warm, primary model)
ev <- run_scn(d, "smog_primary", "Primary") |> mutate(evalue_point = evalue_from_d(std),
        evalue_ci = if_else(std_lo > 0, evalue_from_d(std_lo), if_else(std_hi < 0, evalue_from_d(std_hi), 1)))
write_tab(ev |> select(label, n, std, std_lo, std_hi, p, evalue_point, evalue_ci), "tab7g_evalues")
message("07 done.")