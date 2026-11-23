# =====================================================================================
# 05_lag_subgroup.R  -  WHEN does PM2.5 matter, and FOR WHOM?
#   A. Distributed-lag model (natural-spline over lags 0-14 d; linear in PM2.5), per +50 ug/m3
#   B. Which exposure window fits best (same-day, 3-, 7-, 14-, 28-day means)
#   C. Effect modification: age (<60 / 60+), sex, kidney function (eGFR <60 / >=60)
#   Each analysis is run in (i) all months and (ii) Nov-Feb only (season held fixed).
# =====================================================================================
source("00_setup.R")
d <- load_analysis()
OUTC <- MARKERS |> filter(var %in% c("NLR", "SII", "CBC_WBC_TLC_10e3_uL", "CBC_Platelet_Count_10e3_uL",
                                     "CBC_Hemoglobin_HB_g_dL", "eGFR_CKDEPI2021", "Blood_Glucose_Random_mg_dL", "Eos_Abs_10e3_uL"))
datasets <- list(`All months` = d, `Nov-Feb only` = filter(d, Smog_A_WholeMonths_NovFeb == 1))
DSCOL <- c(`All months` = "black", `Nov-Feb only` = COL_SMOG)
covs <- c("Age_Group", "Sex")

# ------------------------------------------------------------------ A. distributed-lag model
B <- dlm_basis(3)
dlm_fit <- function(dat, var, tr) {
  dd <- dat[complete.cases(dat[, c(var, LAG_COLS, covs)]), ]; dd$.y <- tf(dd[[var]], tr); dd <- dd[is.finite(dd$.y), ]
  if (nrow(dd) < 80) return(NULL)
  Z <- (as.matrix(dd[, LAG_COLS]) / PM_UNIT) %*% B; colnames(Z) <- paste0("cb", seq_len(ncol(Z)))
  dd <- bind_cols(dd, as_tibble(Z))
  fit <- lm(reformulate(c(colnames(Z), covs), response = ".y"), data = dd)
  if (anyNA(coef(fit)[colnames(Z)])) return(NULL)
  V <- sandwich::vcovHC(fit, type = "HC3")[colnames(Z), colnames(Z)]; b <- coef(fit)[colnames(Z)]; sdy <- sd(dd$.y)
  lag_est <- as.numeric(B %*% b); lag_se <- sqrt(diag(B %*% V %*% t(B)))
  w <- colSums(B); cum <- sum(w * b); cum_se <- sqrt(as.numeric(t(w) %*% V %*% w))
  list(lag = tibble(lag = 0:LAG_MAX, est = lag_est / sdy, lo = (lag_est - 1.96 * lag_se) / sdy, hi = (lag_est + 1.96 * lag_se) / sdy),
       cum = tibble(est = cum / sdy, lo = (cum - 1.96 * cum_se) / sdy, hi = (cum + 1.96 * cum_se) / sdy,
                    p = 2 * pnorm(-abs(cum / cum_se)), n = nrow(dd)))
}
dl <- imap_dfr(datasets, function(dat, nm) pmap_dfr(OUTC, function(var, label, system, tr, primary) {
  r <- dlm_fit(dat, var, tr); if (is.null(r)) return(NULL)
  bind_rows(mutate(r$lag, type = "lag"), mutate(r$cum, lag = NA_integer_, type = "cumulative")) |> mutate(label = label, dataset = nm, var = var)
}))
write_tab(dl |> filter(type == "cumulative") |> select(dataset, label, n, est, lo, hi, p), "tab5a_dlm_cumulative_effects")
write_tab(dl |> filter(type == "lag"), "tab5a_dlm_lag_specific")

pl <- dl |> filter(type == "lag") |> mutate(label = factor(label, levels = OUTC$label))
p7 <- ggplot(pl, aes(lag, est, colour = dataset, fill = dataset)) + geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 0.9) +
  scale_colour_manual(values = DSCOL, name = NULL) + scale_fill_manual(values = DSCOL, name = NULL) +
  facet_wrap(~ label, ncol = 4) +
  labs(x = "Lag (days before the test)", y = "Effect per +50 ug/m3 PM2.5 (SD of outcome)", title = "Distributed-lag association, lags 0-14 d (age group + sex adjusted)")
save_fig(p7, "fig7_distributed_lag", 190, 110, tiff = TRUE)

# ------------------------------------------------------------------ B. exposure-window comparison
WIN <- c(pm25_L0_50 = "Same day", pm25_ma3_50 = "3-day mean", pm25_ma7_50 = "7-day mean", pm25_ma14_50 = "14-day mean", pm25_ma28_50 = "28-day mean")
win <- imap_dfr(datasets, function(dat, nm) pmap_dfr(OUTC, function(var, label, system, tr, primary) {
  base <- dat[complete.cases(dat[, c(var, names(WIN), covs)]), ]      # identical patients for every window
  map_dfr(names(WIN), function(x) {
    r <- fit_cont(base, var, tr, x, covs); if (is.null(r)) return(NULL)
    dd <- base; dd$.y <- tf(dd[[var]], tr); aic <- AIC(lm(reformulate(c(x, covs), ".y"), data = dd[is.finite(dd$.y), ]))
    mutate(r, window = WIN[[x]], label = label, dataset = nm, aic = aic)
  })
}))
win <- win |> group_by(dataset, label) |> mutate(delta_aic = aic - min(aic)) |> ungroup()
write_tab(win, "tab5b_exposure_window_comparison")
win$window <- factor(win$window, levels = WIN); win$label <- factor(win$label, levels = OUTC$label)
p8 <- ggplot(win, aes(window, std, colour = dataset)) + geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = std_lo, ymax = std_hi), position = position_dodge(0.5), size = 0.35) +
  scale_colour_manual(values = DSCOL, name = NULL) + facet_wrap(~ label, ncol = 4) +
  labs(x = NULL, y = "Effect per +50 ug/m3 (SD of outcome)", title = "Which exposure window fits best? (identical patients across windows)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p8, "fig8_exposure_windows", 190, 110)

# ------------------------------------------------------------------ C. effect modification
SUB <- list(Age = "Age60", Sex = "Sex", Kidney = "CKD_grp")
d$CKD_grp <- factor(if_else(d$eGFR_CKDEPI2021 < 60, "eGFR <60", "eGFR >=60"), levels = c("eGFR >=60", "eGFR <60"))
datasets <- list(`All months` = d, `Nov-Feb only` = filter(d, Smog_A_WholeMonths_NovFeb == 1))
sub_res <- imap_dfr(datasets, function(dat, nm) pmap_dfr(OUTC, function(var, label, system, tr, primary) {
  imap_dfr(SUB, function(sv, sname) {
    if (var == "eGFR_CKDEPI2021" && sv == "CKD_grp") return(NULL)
    dd <- dat[complete.cases(dat[, c(var, "pm25_ma7_50", sv, covs)]), ]; dd$.y <- tf(dd[[var]], tr); dd <- dd[is.finite(dd$.y), ]
    if (nrow(dd) < 80 || nlevels(droplevels(factor(dd[[sv]]))) < 2) return(NULL)
    cv <- setdiff(covs, if (sv == "Sex") "Sex" else character(0))
    fi <- lm(as.formula(paste(".y ~ pm25_ma7_50 *", sv, "+", paste(cv, collapse = " + "))), data = dd)
    nm_int <- grep("^pm25_ma7_50:", names(coef(fi)), value = TRUE)[1]; V <- sandwich::vcovHC(fi, type = "HC3")
    p_int <- 2 * pt(-abs(coef(fi)[nm_int] / sqrt(V[nm_int, nm_int])), fi$df.residual)
    strat <- map_dfr(levels(droplevels(factor(dd[[sv]]))), function(lv) {
      r <- fit_cont(dd[dd[[sv]] == lv, ], var, tr, "pm25_ma7_50", cv); if (is.null(r)) return(NULL)
      mutate(r, level = lv)
    })
    if (!nrow(strat)) return(NULL)
    strat |> mutate(subgroup = sname, p_interaction = unname(p_int), label = label, dataset = nm)
  })
}))
write_tab(sub_res, "tab5c_effect_modification")
sub_res <- sub_res |> mutate(label = factor(label, levels = OUTC$label), level = factor(level, levels = unique(level)))
pint <- sub_res |> group_by(dataset, label, subgroup) |> summarise(p = first(p_interaction), .groups = "drop") |>
  mutate(sig = p < 0.05)
p9 <- ggplot(sub_res, aes(std, level, colour = dataset)) + geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = std_lo, xmax = std_hi), position = position_dodge(0.6), size = 0.35) +
  scale_colour_manual(values = DSCOL, name = NULL) +
  facet_grid(subgroup ~ label, scales = "free_y", space = "free_y") +
  labs(x = "Effect per +50 ug/m3 7-day PM2.5 (SD of outcome)", y = NULL,
       title = "Effect modification by age, sex and kidney function",
       caption = "Interaction p-values are in tab5c_effect_modification.csv; interpret as exploratory.") +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(size = 7))
save_fig(p9, "fig10_effect_modification", 190, 110)
message("05 done.")