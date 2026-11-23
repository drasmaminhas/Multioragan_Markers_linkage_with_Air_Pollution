# =====================================================================================
# 04_adjusted_models.R  -  adjusted effects of smog / PM2.5 on every marker
#   Exposures : E1 smog vs warm | E2 7-day mean PM2.5 (+50) | E3 14-day mean PM2.5 (+50)
#               E4 7-day mean PM2.5 (+50) restricted to Nov-Feb (season held fixed)
#   Models    : M0 unadjusted | M1 + age group + sex | M2 + temperature & humidity (if weather file)
#               M3 + eGFR (non-renal outcomes only; smaller n)
#   Estimators: linear regression with HC3 robust SE (log-transformed outcomes where skewed);
#               Firth logistic regression for flags (needs {logistf}; falls back to ordinary logistic)
# =====================================================================================
source("00_setup.R")
d <- load_analysis()

EXPOS <- tribble(
  ~id,  ~x,             ~desc,                                          ~subset,
  "E1", "smog",         "Smog vs warm (difference smog - warm)",        "all",
  "E2", "pm25_ma7_50",  "PM2.5 7-day mean, per +50 ug/m3 (all months)", "all",
  "E3", "pm25_ma14_50", "PM2.5 14-day mean, per +50 ug/m3 (all months)","all",
  "E4", "pm25_ma7_50",  "PM2.5 7-day mean, per +50 ug/m3 (Nov-Feb only)","novfeb"
)
MODELS <- c("M0", "M1", "M2", "M3")
subset_data <- function(d, s) if (s == "novfeb") filter(d, Smog_A_WholeMonths_NovFeb == 1) else d

# ------------------------------------------------------------------ continuous outcomes
cont_res <- pmap_dfr(EXPOS, function(id, x, desc, subset) {
  dd <- subset_data(d, subset)
  if (id == "E1") dd <- dd |> filter(!is.na(smog))
  map_dfr(MODELS, function(m) {
    pmap_dfr(MARKERS |> filter(var %in% names(dd)), function(var, label, system, tr, primary) {
      cv <- covars_for(m, dd, var); if (is.null(cv)) return(NULL)
      r <- fit_cont(dd, var, tr, x, cv); if (is.null(r)) return(NULL)
      mutate(r, exposure = id, exposure_desc = desc, model = m, label = label, system = system, tr = tr, primary = primary)
    })
  })
})
cont_res <- cont_res |> group_by(exposure, model) |> mutate(p_fdr = p.adjust(p, "BH")) |> ungroup()
write_tab(cont_res, "tab4_adjusted_continuous_all_models")

# ------------------------------------------------------------------ binary outcomes (Firth logistic)
bin_res <- pmap_dfr(EXPOS, function(id, x, desc, subset) {
  dd <- subset_data(d, subset)
  if (id == "E1") dd <- dd |> filter(!is.na(smog))
  map_dfr(c("M0", "M1", "M2"), function(m) {          # M3 skipped for binary outcomes (too few events)
    pmap_dfr(BINARY |> filter(var %in% names(dd)), function(var, label, system, primary) {
      cv <- covars_for(m, dd, var); if (is.null(cv)) return(NULL)
      r <- fit_bin(dd, var, x, cv); if (is.null(r)) return(NULL)
      mutate(r, exposure = id, exposure_desc = desc, model = m, label = label, system = system, primary = primary)
    })
  })
})
if (nrow(bin_res)) bin_res <- bin_res |> group_by(exposure, model) |> mutate(p_fdr = p.adjust(p, "BH")) |> ungroup()
write_tab(bin_res, "tab4b_adjusted_binary_all_models")

# ------------------------------------------------------------------ Table 3 (paper table): M1 and M2 per exposure
within_model <- if (has_weather(d)) "M2" else "M1"
message("Within-season adjusted model = ", within_model, " | Between-season primary model = M1 (plus M2 sensitivity)")
t3 <- cont_res |> filter(model %in% c("M1", within_model), exposure %in% c("E1", "E2", "E4")) |>
  mutate(effect = if_else(tr == "none", sprintf("%.2f (%.2f, %.2f)", est, lo, hi),
                          sprintf("%.1f%% (%.1f, %.1f)", pct, pct_lo, pct_hi)),
         std_effect = sprintf("%.2f (%.2f, %.2f)", std, std_lo, std_hi)) |>
  select(exposure, model, label, system, n, effect, std_effect, p, p_fdr) |> arrange(exposure, model, system, label)
write_tab(t3, "tab4c_main_adjusted_effects")

# ------------------------------------------------------------------ Figure 5: forest of standardised adjusted effects
fplot <- function(ex, mod, title) {
  f <- cont_res |> filter(exposure == ex, model == mod) |> mutate(sig = if_else(p_fdr < 0.05, "FDR < 0.05", "n.s."),
         system = factor(system, levels = unique(MARKERS$system))) |> arrange(system, std) |> mutate(label = factor(label, levels = unique(label)))
  ggplot(f, aes(std, label, colour = sig)) + geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = std_lo, xmax = std_hi), height = 0) + geom_point(size = 2) +
    scale_colour_manual(values = c("FDR < 0.05" = COL_SMOG, "n.s." = "grey40"), name = NULL) +
    facet_grid(system ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "Standardised effect (SD of outcome, 95% CI)", y = NULL, title = title) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1))
}
p5 <- fplot("E1", "M1", "A. Smog vs warm (M1: age + sex)") +
      fplot("E4", within_model, paste0("B. +50 ug/m3 PM2.5, Nov-Feb only (", within_model, ")")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(p5, "fig5_adjusted_forest", 190, 160, tiff = TRUE)

# ------------------------------------------------------------------ Figure 5b: binary adjusted odds ratios
if (nrow(bin_res)) {
  fb <- bin_res |> filter(model == "M1", exposure %in% c("E1", "E4"), is.finite(or_lo), is.finite(or_hi)) |>
    mutate(ex = if_else(exposure == "E1", "Smog vs warm", "+50 ug/m3 PM2.5 (Nov-Feb)"))
  p5b <- ggplot(fb, aes(or, reorder(label, or), colour = ex)) + geom_vline(xintercept = 1, linetype = "dashed") +
    geom_errorbarh(aes(xmin = or_lo, xmax = or_hi), height = 0, position = position_dodge(0.5)) +
    geom_point(position = position_dodge(0.5), size = 2) + scale_x_log10() +
    scale_colour_manual(values = c("Smog vs warm" = COL_SMOG, "+50 ug/m3 PM2.5 (Nov-Feb)" = "black"), name = NULL) +
    labs(x = "Adjusted odds ratio (age group + sex; Firth where {logistf} installed)", y = NULL, title = "Adjusted odds of abnormality")
  save_fig(p5b, "fig5b_adjusted_binary_forest", 150, 110)
}

# ------------------------------------------------------------------ Figure 6: GAM exposure-response curves (7-day PM2.5)
gam_curve <- function(var, tr, label) {
  dd <- d[complete.cases(d[, c(var, "pm25_ma7", "Age_Group", "Sex")]), ]; dd$.y <- tf(dd[[var]], tr); dd <- dd[is.finite(dd$.y), ]
  if (nrow(dd) < 60) return(NULL)
  g <- gam(.y ~ s(pm25_ma7, k = 4) + Age_Group + Sex, data = dd, method = "REML")
  grid <- tibble(pm25_ma7 = seq(quantile(dd$pm25_ma7, .02), quantile(dd$pm25_ma7, .98), length.out = 100),
                 Age_Group = factor("60-69", levels = levels(dd$Age_Group)), Sex = factor("Male", levels = levels(dd$Sex)))
  pr <- predict(g, grid, se.fit = TRUE); bt <- if (tr %in% c("log", "log1p")) exp else identity
  st <- summary(g)$s.table
  list(curve = grid |> mutate(fit = bt(pr$fit), lo = bt(pr$fit - 1.96 * pr$se.fit), hi = bt(pr$fit + 1.96 * pr$se.fit), label = label),
       rug = dd |> transmute(pm25_ma7, label = label), ann = tibble(label = label, edf = st[1, "edf"], p = st[1, "p-value"], n = nrow(dd)))
}
gl <- MARKERS |> filter(var %in% KEY6) |> pmap(function(var, label, system, tr, primary) gam_curve(var, tr, label)) |> compact()
gc <- map_dfr(gl, "curve"); gr <- map_dfr(gl, "rug"); ga <- map_dfr(gl, "ann") |>
  mutate(txt = sprintf("n = %d\nedf = %.1f, %s", n, edf, ifelse(p < 0.001, "p < 0.001", paste0("p = ", fmt_p(p)))))
gc$label <- factor(gc$label, levels = lab_of(KEY6)); gr$label <- factor(gr$label, levels = lab_of(KEY6)); ga$label <- factor(ga$label, levels = lab_of(KEY6))
p6 <- ggplot(gc, aes(pm25_ma7, fit)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = COL_SMOG, alpha = 0.2) +
  geom_line(colour = COL_SMOG, linewidth = 1) + geom_rug(data = gr, aes(x = pm25_ma7), inherit.aes = FALSE, alpha = 0.15, sides = "b") +
  geom_text(data = ga, aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.2, size = 2.8, inherit.aes = FALSE) +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(x = expression("7-day mean PM"[2.5]~"("*mu*g/m^3*")"), y = "Predicted value",
       title = "Exposure-response curves (GAM)",
       caption = "Adjusted for age group and sex; reference = male, 60-69 y. Log-scaled markers shown as geometric means. Shaded band = 95% CI.")
save_fig(p6, "fig6_gam_exposure_response", 190, 120, tiff = TRUE)

# ------------------------------------------------------------------ Troponin positivity by PM2.5 tertile
tv <- if ("Trop_Elevated_Any" %in% names(d) && sum(d$Trop_Elevated_Any, na.rm = TRUE) >= 5) "Trop_Elevated_Any" else "Trop_Elevated_QL"
dt <- d |> filter(!is.na(.data[[tv]]))
ev <- sum(dt[[tv]] == 1)
message("Troponin outcome used: ", tv, " | tested n = ", nrow(dt), " | positives = ", ev)
if (ev >= 5) {
  dt$tert <- cut(dt$pm25_ma7, quantile(d$pm25_ma7, c(0, 1/3, 2/3, 1), na.rm = TRUE), include.lowest = TRUE, labels = c("T1 (low)", "T2", "T3 (high)"))
  tt <- dt |> group_by(tert) |> summarise(n = n(), x = sum(.data[[tv]] == 1), .groups = "drop") |>
    mutate(pct = 100 * x / n, lo = 100 * map2_dbl(x, n, ~ wilson(.x, .y)[1]), hi = 100 * map2_dbl(x, n, ~ wilson(.x, .y)[2]))
  ptr <- prop.trend.test(tt$x, tt$n)$p.value; write_tab(tt, "tab4d_troponin_by_pm25_tertile")
  p9 <- ggplot(tt, aes(tert, pct)) + geom_col(fill = COL_SMOG, alpha = 0.8, width = 0.6) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
    geom_text(aes(label = sprintf("%d/%d", x, n)), vjust = -0.5, hjust = -0.2, size = 3) +
    labs(x = "Tertile of 7-day mean PM2.5", y = "Troponin-positive (%)", title = "Troponin positivity by PM2.5 tertile",
         subtitle = sprintf("Cochran-Armitage trend p = %s; outcome = %s", fmt_p(ptr), tv))
  save_fig(p9, "fig9_troponin_by_pm25_tertile", 110, 90)
} else message("  Too few troponin-positive patients (<5): troponin-by-tertile figure skipped. ",
               "Set the hs-troponin cut-off in script 02 of data prep (--hs-cutoff) to add more events.")
message("04 done.")