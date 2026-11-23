# =====================================================================================
# 09_within_month.R  -  the cleanest test available in this dataset
#   Question: on days with higher PM2.5, do patients seen in the SAME month (or same half-month) have different markers?
#   Month / half-month fixed effects remove everything that is constant inside a period:
#   analyser calibration, reagent lots, reporting habits (rounded differentials), department / referral mix, season, temperature trend.
#   What is left is day-to-day PM2.5 variation inside the period.
#   Also shown: the same model with FOLLOWING-week PM2.5 (negative-control exposure; should be ~0).
# =====================================================================================
source("00_setup.R")
d <- load_analysis()
d$Period <- factor(paste0(format(d$Est_Date, "%Y-%m"), ifelse(as.integer(format(d$Est_Date, "%d")) <= 15, "a", "b")))

OUT <- MARKERS |> filter(var %in% c("NLR", "SII", "CBC_WBC_TLC_10e3_uL", "CBC_Platelet_Count_10e3_uL", "CBC_Hemoglobin_HB_g_dL",
          "CBC_HCT_Percent", "CBC_Total_RBCs_10e6_uL", "CBC_MCV_fL", "eGFR_CKDEPI2021", "Blood_Glucose_Random_mg_dL"))
ANALYSER_AFFECTED <- c("CBC_Hemoglobin_HB_g_dL")      # Hb/MCH/MCHC drift across months (script 08); shown but flagged
MIN_CELL <- 15                                       # drop months / half-months with fewer reports than this (per outcome)
EXPS <- c(pm25_ma7_50 = "Preceding 7 days (real exposure)", pm25_ma14_50 = "Preceding 14 days", pm25_L0_50 = "Same day",
          pm25_lead7_50 = "FOLLOWING 7 days (negative control)")

within_sd <- d |> filter(Smog_A_WholeMonths_NovFeb == 1 | TRUE) |> group_by(Month_Label) |>
  summarise(n = n(), pm25_ma7_mean = mean(pm25_ma7), pm25_ma7_sd = sd(pm25_ma7), pm25_ma7_min = min(pm25_ma7), pm25_ma7_max = max(pm25_ma7), .groups = "drop")
write_tab(within_sd, "tab9a_within_month_pm25_variation")

run_fe <- function(cell, label_fe) pmap_dfr(OUT, function(var, label, system, tr, primary) {
  imap_dfr(EXPS, function(elab, x) {
    dd <- d[complete.cases(d[, c(var, x, "Age_Group", "Sex", cell)]), ]
    keep <- names(which(table(dd[[cell]]) >= MIN_CELL)); dd <- dd[dd[[cell]] %in% keep, ]; dd[[cell]] <- droplevels(factor(dd[[cell]]))
    if (nrow(dd) < 60 || nlevels(dd[[cell]]) < 2) return(NULL)
    r <- fit_cont(dd, var, tr, x, c(cell, "Age_Group", "Sex")); if (is.null(r)) return(NULL)
    mutate(r, label = label, exposure = elab, fe = label_fe, n_cells = nlevels(dd[[cell]]),
           analyser_flag = var %in% ANALYSER_AFFECTED)
  })
})
res <- bind_rows(run_fe("Month_Label", "Month fixed effects"), run_fe("Period", "Half-month fixed effects"))
write_tab(res |> select(fe, exposure, label, n, n_cells, est, std, std_lo, std_hi, p, analyser_flag), "tab9b_within_month_effects")

pl <- res |> filter(exposure %in% c("Preceding 7 days (real exposure)", "FOLLOWING 7 days (negative control)")) |>
  mutate(label = factor(ifelse(analyser_flag, paste0(label, " *"), label), levels = unique(ifelse(OUT$label %in% lab_of(ANALYSER_AFFECTED), paste0(OUT$label, " *"), OUT$label))))
p <- ggplot(pl, aes(std, label, colour = exposure)) + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = std_lo, xmax = std_hi), position = position_dodge(0.6), size = 0.35) +
  scale_colour_manual(values = c("Preceding 7 days (real exposure)" = COL_SMOG, "FOLLOWING 7 days (negative control)" = "grey50"), name = NULL) +
  facet_wrap(~ fe) + guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Effect of +50 ug/m3 PM2.5 (SD of outcome, 95% CI) after removing month / half-month differences", y = NULL,
       title = "Within-period PM2.5 effects", caption = "* haemoglobin shows analyser drift between months (script 08). Adjusted for age group and sex.")
save_fig(p, "fig18_within_month_effects", 190, 110, tiff = TRUE)
message("09 done.")
