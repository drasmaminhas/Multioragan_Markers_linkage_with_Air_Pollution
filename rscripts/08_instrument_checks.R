# =====================================================================================
# 08_instrument_checks.R  -  is the smog-vs-warm difference real biology, or analyser / rounding / case-mix?
#   A. Red-cell indices over time (MCHC, MCV, MCH, Hb, HCT, RBC): step change (= analyser/batch) or smooth drift?
#   B. Differential-count heaping (values rounded to 5s, default 3/2 monocyte-eosinophil pairs) by month
#      + NLR/SII re-estimated ONLY on reports whose differential does not look rounded
#   C. Case-mix by month (age, sex, tests ordered, anaemia prevalence by sex)
#   D. Weekly medians of key markers drawn under the daily PM2.5 curve (do markers FOLLOW pollution in time?)
# =====================================================================================
source("00_setup.R")
d   <- load_analysis()
day <- read_csv(file.path(DAT_DIR, "aqi_daily_processed.csv"), show_col_types = FALSE) |> mutate(date = as.Date(date))
mcol <- c("Smog months (Nov-Feb)" = COL_SMOG, "Warm months (Mar-Jul)" = COL_WARM)

# ------------------------------------------------------------------ A. analyser indices by month
IDX <- c(CBC_MCHC_g_dL = "MCHC (g/dL)", CBC_MCV_fL = "MCV (fL)", CBC_MCH_pg = "MCH (pg)",
         CBC_Hemoglobin_HB_g_dL = "Haemoglobin (g/dL)", CBC_HCT_Percent = "Haematocrit (%)", CBC_Total_RBCs_10e6_uL = "RBC (10^6/uL)")
idx_long <- d |> select(Month_Label, Month_Season, Sex, all_of(names(IDX))) |>
  pivot_longer(all_of(names(IDX)), names_to = "var", values_to = "v") |> filter(!is.na(v)) |>
  mutate(label = factor(IDX[var], levels = IDX))
mon_tab <- idx_long |> group_by(label, Month_Label) |> summarise(n = n(), median = median(v), q1 = quantile(v, .25), q3 = quantile(v, .75), .groups = "drop")
write_tab(mon_tab, "tab8a_analyser_indices_by_month")
kw <- idx_long |> group_by(label) |> summarise(kw_p = kruskal.test(v ~ Month_Label)$p.value, .groups = "drop"); write_tab(kw, "tab8a_kruskal_by_month")
pA <- idx_long |> group_by(label, Month_Label) |> filter(n() >= 5) |> ungroup() |>
  ggplot(aes(Month_Label, v, fill = Month_Season, colour = Month_Season)) + geom_boxplot(alpha = 0.35, outlier.shape = NA, linewidth = 0.4) +
  ggbeeswarm::geom_quasirandom(size = 0.4, alpha = 0.4, width = 0.2) + scale_fill_manual(values = mcol, name = NULL) + scale_colour_manual(values = mcol, name = NULL) +
  facet_wrap(~ label, scales = "free_y", ncol = 3) + coord_cartesian() +
  labs(x = NULL, y = NULL, title = "Red-cell indices by month", subtitle = "A step between two adjacent months = analyser/batch; a smooth drift = biology or case-mix") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_fig(pA, "fig16_analyser_indices_by_month", 190, 120, tiff = TRUE)

# ------------------------------------------------------------------ B. differential heaping
d <- d |> mutate(
  neut_m5 = CBC_Neutrophils_Percent %% 5 == 0, lym_m5 = CBC_Lymphocytes_Percent %% 5 == 0,
  mono_eos_32 = CBC_Monocytes_Percent == 3 & CBC_Eosinophils_Percent == 2,
  mono_eos_55 = CBC_Monocytes_Percent == 5 & CBC_Eosinophils_Percent == 5,
  diff_rounded = neut_m5 & lym_m5)
heap <- d |> filter(!is.na(CBC_Neutrophils_Percent)) |> group_by(Month_Label) |>
  summarise(n = n(), pct_neut_mult5 = 100 * mean(neut_m5), pct_lym_mult5 = 100 * mean(lym_m5), pct_both_mult5 = 100 * mean(diff_rounded),
            pct_mono3_eos2 = 100 * mean(mono_eos_32, na.rm = TRUE), pct_mono5_eos5 = 100 * mean(mono_eos_55, na.rm = TRUE), .groups = "drop")
write_tab(heap, "tab8b_differential_heaping_by_month")
pB <- heap |> pivot_longer(starts_with("pct"), names_to = "metric", values_to = "pct") |>
  mutate(metric = recode(metric, pct_neut_mult5 = "Neutrophil % multiple of 5", pct_lym_mult5 = "Lymphocyte % multiple of 5",
                         pct_both_mult5 = "Both multiples of 5", pct_mono3_eos2 = "Monocyte 3 & eosinophil 2", pct_mono5_eos5 = "Monocyte 5 & eosinophil 5")) |>
  ggplot(aes(Month_Label, pct, group = metric, colour = metric)) + geom_line(linewidth = 0.7) + geom_point(size = 2) +
  scale_colour_manual(values = OKABE[c(2, 3, 5, 6, 7)], name = NULL) + guides(colour = guide_legend(nrow = 2)) +
  labs(x = NULL, y = "% of reports", title = "Rounding / default patterns in the differential count",
       subtitle = "~20% would be multiples of 5 by chance; much higher = manual/rounded differential")
save_fig(pB, "figS8_differential_heaping", 150, 90)

heap_sens <- map_dfr(list(`All reports` = d, `Differential NOT rounded` = filter(d, !diff_rounded %in% TRUE),
                          `Differential rounded` = filter(d, diff_rounded %in% TRUE)), .id = "subset", function(dd) {
  MARKERS |> filter(var %in% c("NLR", "SII", "CBC_WBC_TLC_10e3_uL", "Eos_Abs_10e3_uL", "CBC_Hemoglobin_HB_g_dL", "CBC_Platelet_Count_10e3_uL")) |>
    pmap_dfr(function(var, label, system, tr, primary) { r <- fit_cont(dd, var, tr, "smog", c("Age_Group", "Sex")); if (is.null(r)) NULL else mutate(r, label = label) }) })
write_tab(heap_sens |> select(subset, label, n, std, std_lo, std_hi, p), "tab8b_smog_effect_by_differential_quality")
pB2 <- heap_sens |> mutate(label = factor(label, levels = unique(label))) |>
  ggplot(aes(std, subset, colour = subset)) + geom_vline(xintercept = 0, linetype = "dashed") + geom_pointrange(aes(xmin = std_lo, xmax = std_hi), size = 0.35) +
  scale_colour_manual(values = c("grey30", COL_SMOG, COL_WARM), guide = "none") + facet_wrap(~ label, ncol = 3) +
  labs(x = "Smog minus warm (SD of outcome, age group + sex adjusted)", y = NULL, title = "Does the smog effect depend on how the differential was reported?")
save_fig(pB2, "figS9_effect_by_differential_quality", 190, 100)

# ------------------------------------------------------------------ C. case-mix by month
cm <- d |> group_by(Month_Label) |> summarise(
  n = n(), age_median = median(Age_Years, na.rm = TRUE), pct_female = 100 * mean(Sex_Female == 1, na.rm = TRUE),
  pct_age60plus = 100 * mean(Age_Years >= 60, na.rm = TRUE), pct_has_cbc = 100 * mean(Has_CBC == 1), pct_has_renal = 100 * mean(Has_RFT == 1),
  pct_has_trop = 100 * mean(Has_Troponin == 1), pct_has_glucose = 100 * mean(Has_Glucose == 1),
  anaemia_pct_all = 100 * mean(Anaemia == 1, na.rm = TRUE),
  anaemia_pct_female = 100 * mean(Anaemia[Sex_Female == 1] == 1, na.rm = TRUE), anaemia_pct_male = 100 * mean(Anaemia[Sex_Female == 0] == 1, na.rm = TRUE),
  hb_median_female = median(CBC_Hemoglobin_HB_g_dL[Sex_Female == 1], na.rm = TRUE), hb_median_male = median(CBC_Hemoglobin_HB_g_dL[Sex_Female == 0], na.rm = TRUE),
  .groups = "drop") |> mutate(across(where(is.numeric), ~ round(.x, 1)))
write_tab(cm, "tab8c_casemix_by_month")

# ------------------------------------------------------------------ D. weekly medians under the PM2.5 curve
TS <- c("NLR", "SII", "CBC_WBC_TLC_10e3_uL", "CBC_Platelet_Count_10e3_uL", "CBC_Hemoglobin_HB_g_dL", "CBC_MCHC_g_dL")
tsl <- d |> select(Est_Date, all_of(TS)) |> pivot_longer(all_of(TS), names_to = "var", values_to = "v") |> filter(!is.na(v)) |>
  mutate(week = floor_date(Est_Date, "week", week_start = 1), label = factor(c(setNames(MARKERS$label, MARKERS$var), CBC_MCHC_g_dL = "MCHC")[var],
         levels = c("NLR", "SII", "WBC", "Platelets", "Haemoglobin", "MCHC")))
wk <- tsl |> group_by(label, week) |> summarise(n = n(), med = median(v), q1 = quantile(v, .25), q3 = quantile(v, .75), .groups = "drop") |> filter(n >= 8)
write_tab(wk, "tab8d_weekly_medians")
xl <- range(d$Est_Date, na.rm = TRUE) + c(-3, 3)
pTop <- ggplot(filter(day, date >= xl[1], date <= xl[2]), aes(date)) + geom_line(aes(y = pm25), colour = "grey65", linewidth = 0.3) +
  geom_line(aes(y = pm25_ma7), linewidth = 0.8) + scale_x_date(limits = xl, date_breaks = "1 month", date_labels = "%b") +
  labs(x = NULL, y = expression(PM[2.5]), title = "PM2.5 (top) and weekly marker medians (weeks with n >= 8)") + theme(axis.text.x = element_blank())
pBot <- ggplot(wk, aes(week, med)) + geom_ribbon(aes(ymin = q1, ymax = q3), fill = "grey85") + geom_line(colour = "black") +
  geom_point(aes(size = n), colour = COL_SMOG, alpha = 0.8) + scale_size_continuous(range = c(1, 3.5), name = "Reports/week") +
  scale_x_date(limits = xl, date_breaks = "1 month", date_labels = "%b") + facet_wrap(~ label, scales = "free_y", ncol = 1, strip.position = "right") +
  labs(x = NULL, y = "Weekly median (IQR band)") + theme(strip.text.y = element_text(angle = 0))
save_fig(pTop / pBot + plot_layout(heights = c(1, 5)), "fig17_markers_vs_pm25_over_time", 180, 220, tiff = TRUE)
message("08 done.")
