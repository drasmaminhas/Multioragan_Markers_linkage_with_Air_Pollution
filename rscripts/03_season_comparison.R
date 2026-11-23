# =====================================================================================
# 03_season_comparison.R  -  smog vs warm comparison of every marker (unadjusted)
#   Tests   : Mann-Whitney (Hodges-Lehmann difference + 95% CI), rank-biserial r, Cohen's d (transformed scale),
#             Fisher exact for flags; BH-FDR across all markers, Holm across the primary markers
#   Months  : Kruskal-Wallis, pairwise Wilcoxon (Holm), Spearman trend across Nov-Feb
#   Figures : forest plots (continuous, binary), monthly box plots + PM2.5 bars, prevalence by month
# =====================================================================================
source("00_setup.R")
d   <- load_analysis()
day <- read_csv(file.path(DAT_DIR, "aqi_daily_processed.csv"), show_col_types = FALSE)

# ------------------------------------------------------------------ comparison functions
cont_compare <- function(d, var, tr, smog_var) {
  dd <- d[!is.na(d[[smog_var]]) & !is.na(d[[var]]), ]
  x1 <- dd[[var]][dd[[smog_var]] == 1]; x0 <- dd[[var]][dd[[smog_var]] == 0]
  if (length(x1) < 5 || length(x0) < 5) return(NULL)
  w  <- suppressWarnings(wilcox.test(x1, x0, conf.int = TRUE, exact = FALSE))
  y1 <- tf(x1, tr); y0 <- tf(x0, tr); y1 <- y1[is.finite(y1)]; y0 <- y0[is.finite(y0)]; cd <- cohens_d(y1, y0)
  dig <- if (var == "Troponin_I_hs_pg_mL") 2 else 1
  tibble(var = var, n_smog = length(x1), n_warm = length(x0), smog_median_iqr = med_iqr(x1, dig), warm_median_iqr = med_iqr(x0, dig),
         hl_diff = unname(w$estimate), hl_lo = w$conf.int[1], hl_hi = w$conf.int[2], p = w$p.value,
         rank_biserial = rank_biserial(x1, x0), d = cd["d"], d_lo = cd["lo"], d_hi = cd["hi"],
         geo_mean_ratio = if (tr %in% c("log", "log1p")) exp(mean(y1) - mean(y0)) else NA_real_)
}
bin_compare <- function(d, var, smog_var) {
  dd <- d[!is.na(d[[smog_var]]) & !is.na(d[[var]]), ]
  a <- sum(dd[[var]] == 1 & dd[[smog_var]] == 1); n1 <- sum(dd[[smog_var]] == 1)
  b <- sum(dd[[var]] == 1 & dd[[smog_var]] == 0); n0 <- sum(dd[[smog_var]] == 0)
  if (n1 < 5 || n0 < 5 || (a + b) < 2) return(NULL)
  ft <- fisher.test(matrix(c(a, n1 - a, b, n0 - b), 2, byrow = TRUE))
  tibble(var = var, n_smog = n1, n_warm = n0, ev_smog = a, ev_warm = b, prev_smog = 100 * a / n1, prev_warm = 100 * b / n0,
         risk_diff = 100 * (a / n1 - b / n0), or = unname(ft$estimate), or_lo = ft$conf.int[1], or_hi = ft$conf.int[2], p = ft$p.value)
}
add_meta <- function(r, ...) if (is.null(r)) NULL else mutate(r, ...)
run_definition <- function(smog_var) {
  cc <- MARKERS |> filter(var %in% names(d)) |> pmap_dfr(function(var, label, system, tr, primary)
          add_meta(cont_compare(d, var, tr, smog_var), label = label, system = system, primary = primary))
  bb <- BINARY  |> filter(var %in% names(d)) |> pmap_dfr(function(var, label, system, primary)
          add_meta(bin_compare(d, var, smog_var), label = label, system = system, primary = primary))
  allp <- c(cc$p, bb$p); fdr <- p.adjust(allp, "BH")
  cc$p_fdr <- fdr[seq_len(nrow(cc))]; bb$p_fdr <- fdr[-seq_len(nrow(cc))]
  cc$p_holm_primary <- NA_real_; bb$p_holm_primary <- NA_real_
  pr <- c(cc$p[cc$primary], bb$p[bb$primary]); ph <- p.adjust(pr, "holm")
  cc$p_holm_primary[cc$primary] <- ph[seq_len(sum(cc$primary))]; bb$p_holm_primary[bb$primary] <- ph[-seq_len(sum(cc$primary))]
  list(cont = cc, bin = bb)
}
res <- map(c(A = "Smog_A_WholeMonths_NovFeb", B = "Smog_B_Nov15_Feb15", C = "Smog_C_Nov15_Jan31"), run_definition)
for (k in names(res)) { write_tab(res[[k]]$cont, paste0("tab3_smog_vs_warm_continuous_def", k))
                        write_tab(res[[k]]$bin,  paste0("tab3_smog_vs_warm_binary_def", k)) }
prim <- if (grepl("_A_", PRIMARY_SMOG)) "A" else if (grepl("_C_", PRIMARY_SMOG)) "C" else "B"
cc <- res[[prim]]$cont; bb <- res[[prim]]$bin
cat("\nPrimary definition:", PRIMARY_SMOG, "| markers significant after BH-FDR (continuous):", sum(cc$p_fdr < 0.05), "of", nrow(cc), "\n")

# ------------------------------------------------------------------ Figure 4A: forest of standardised differences
sys_levels <- unique(MARKERS$system)
fc <- cc |> mutate(system = factor(system, levels = sys_levels), sig = if_else(p_fdr < 0.05, "FDR < 0.05", "n.s.")) |>
  arrange(system, d) |> mutate(label = factor(label, levels = unique(label)))
p4a <- ggplot(fc, aes(d, label, colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed") + geom_errorbarh(aes(xmin = d_lo, xmax = d_hi), height = 0, linewidth = 0.5) +
  geom_point(size = 2) + scale_colour_manual(values = c("FDR < 0.05" = COL_SMOG, "n.s." = "grey40"), name = NULL) +
  facet_grid(system ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Standardised difference, smog minus warm (Cohen's d, 95% CI)", y = NULL, title = "A. Continuous markers (log-transformed where skewed)") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1))
fb <- bb |> filter(is.finite(or_lo), is.finite(or_hi), or_lo > 0) |> mutate(system = factor(system, levels = unique(BINARY$system)),
        sig = if_else(p_fdr < 0.05, "FDR < 0.05", "n.s.")) |> arrange(system, or) |> mutate(label = factor(label, levels = unique(label)))
p4b <- ggplot(fb, aes(or, label, colour = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed") + geom_errorbarh(aes(xmin = or_lo, xmax = or_hi), height = 0) + geom_point(size = 2) +
  scale_x_log10() + scale_colour_manual(values = c("FDR < 0.05" = COL_SMOG, "n.s." = "grey40"), name = NULL) +
  labs(x = "Odds ratio, smog vs warm (Fisher, log scale)", y = NULL, title = "B. Abnormality flags")
save_fig(p4a + p4b + plot_layout(widths = c(1.2, 1), guides = "collect") & theme(legend.position = "bottom"),
         "fig4_forest_smog_vs_warm", 190, 160, tiff = TRUE)

# ------------------------------------------------------------------ Figure 3: monthly distributions + PM2.5 bars
mon_pm <- day |> mutate(Month_Label = format(date, "%b-%Y")) |> filter(Month_Label %in% levels(d$Month_Label)) |>
  group_by(Month_Label) |> summarise(pm25 = mean(pm25, na.rm = TRUE), .groups = "drop") |>
  mutate(Month_Label = factor(Month_Label, levels = levels(d$Month_Label)),
         Month_Season = if_else(str_detect(as.character(Month_Label), "^(Nov|Dec|Jan|Feb)"), "Smog months (Nov-Feb)", "Warm months (Mar-Jul)"))
mcol <- c("Smog months (Nov-Feb)" = COL_SMOG, "Warm months (Mar-Jul)" = COL_WARM)
p_pm <- ggplot(mon_pm, aes(Month_Label, pm25, fill = Month_Season)) + geom_col(width = 0.7) +
  geom_text(aes(label = round(pm25)), vjust = -0.3, size = 3) + scale_fill_manual(values = mcol, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + labs(x = NULL, y = expression(PM[2.5]~"("*mu*g/m^3*")"),
  title = "A. Monthly mean PM2.5 (all days)") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
key <- MARKERS |> filter(var %in% KEY6) |> mutate(var = factor(var, levels = KEY6)) |> arrange(var)
long <- d |> select(Month_Label, Month_Season, all_of(KEY6)) |> pivot_longer(all_of(KEY6), names_to = "var", values_to = "value") |>
  filter(!is.na(value)) |> mutate(label = factor(lab_of(var), levels = key$label))
long <- long |> group_by(label, Month_Label) |> filter(n() >= 5) |> ungroup()
logsc <- key$label[key$tr == "log"]
pts_layer <- if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
  ggbeeswarm::geom_quasirandom(size = 0.5, alpha = 0.5, width = 0.2)
} else {
  geom_jitter(size = 0.5, alpha = 0.5, width = 0.15)
}
mk_box <- function(lbl) {
  dd <- long |> filter(label == lbl)
  g <- ggplot(dd, aes(Month_Label, value, fill = Month_Season, colour = Month_Season)) +
    geom_boxplot(alpha = 0.35, outlier.shape = NA, linewidth = 0.4) + pts_layer +
    scale_fill_manual(values = mcol, guide = "none") + scale_colour_manual(values = mcol, guide = "none") + labs(x = NULL, y = NULL, title = lbl) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), plot.title = element_text(size = 9))
  if (lbl %in% logsc) g <- g + scale_y_log10()
  g
}
p_box <- wrap_plots(map(key$label, mk_box), ncol = 3) + plot_annotation(title = "B. Monthly distributions (points = patients; months with n<5 omitted)",
         theme = theme(plot.title = element_text(face = "bold")))
save_fig(wrap_elements(p_pm) / wrap_elements(p_box) + plot_layout(heights = c(1, 3.2)) , "fig3_monthly_markers", 190, 200, tiff = TRUE)

# ------------------------------------------------------------------ Monthly tests: KW, pairwise Wilcoxon, trend Nov-Feb
mt <- MARKERS |> filter(var %in% names(d)) |> pmap_dfr(function(var, label, system, tr, primary) {
  dd <- d[!is.na(d[[var]]), ]; kw <- tryCatch(kruskal.test(dd[[var]], dd$Month_Label)$p.value, error = function(e) NA_real_)
  nf <- dd[dd$Lab_Month %in% c(11, 12, 1, 2), ]
  sp <- tryCatch(suppressWarnings(cor.test(nf[[var]], nf$Month_Index, method = "spearman", exact = FALSE)), error = function(e) NULL)
  tibble(var = var, label = label, kw_p_months = kw, trend_rho_NovFeb = if (!is.null(sp)) unname(sp$estimate) else NA,
         trend_p_NovFeb = if (!is.null(sp)) sp$p.value else NA)
}) |> mutate(kw_p_fdr = p.adjust(kw_p_months, "BH"), trend_p_fdr = p.adjust(trend_p_NovFeb, "BH"))
write_tab(mt, "tab3b_monthly_tests")
pw <- map_dfr(KEY6, function(v) {
  dd <- d[!is.na(d[[v]]), ]; r <- pairwise.wilcox.test(dd[[v]], dd$Month_Label, p.adjust.method = "holm", exact = FALSE)$p.value
  as_tibble(r, rownames = "month_a") |> pivot_longer(-month_a, names_to = "month_b", values_to = "p_holm") |> filter(!is.na(p_holm)) |> mutate(var = v)
})
write_tab(pw, "tab3c_pairwise_months_key6")

# ------------------------------------------------------------------ Figure: prevalence of key flags by month (Wilson CI)
fl <- c("Trop_Elevated_QL", "Anaemia", "Leukocytosis", "Thrombocytosis", "Neutrophilia_ANC", "Eosinophilia_AEC", "eGFR_lt60", "Hyperglycaemia_Random_ge200")
pv <- map_dfr(fl, function(v) d |> filter(!is.na(.data[[v]])) |> group_by(Month_Label, Month_Season) |>
  summarise(n = n(), x = sum(.data[[v]] == 1), .groups = "drop") |> filter(n >= 10) |>
  mutate(var = v, pct = 100 * x / n, lo = 100 * map2_dbl(x, n, ~ wilson(.x, .y)[1]), hi = 100 * map2_dbl(x, n, ~ wilson(.x, .y)[2])))
pv$label <- factor(lab_of(pv$var), levels = lab_of(fl))
pPrev <- ggplot(pv, aes(Month_Label, pct, colour = Month_Season)) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) + geom_point(size = 2) +
  scale_colour_manual(values = mcol, name = NULL) + facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Prevalence (%) with 95% Wilson CI", title = "Prevalence of abnormalities by month (months with n < 10 tested omitted)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_fig(pPrev, "fig3b_prevalence_by_month", 190, 140)
message("03 done.")