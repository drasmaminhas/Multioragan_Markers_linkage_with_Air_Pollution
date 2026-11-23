# =====================================================================================
# 06_exploratory.R  -  hypothesis-generating analyses (label as EXPLORATORY in the paper)
#   A. Spearman correlation: air-pollution / weather windows  x  markers   (BH-FDR)
#   B. Marker-marker correlation matrix (hierarchically ordered)
#   C. PCA of the blood-count profile; do scores track season / PM2.5?
#   D. k-means phenotypes on PCs; are winter patients over-represented in a cluster?
#   E. Random forest: can the blood profile tell smog from warm?  (5-fold CV AUC + label-permutation null)
# =====================================================================================
source("00_setup.R")
for (pkg in c("ggrepel", "cluster", "ranger", "pROC")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, quiet = TRUE)
}
d <- load_analysis()
EXP <- c(pm25_L0 = "PM2.5 same day", pm25_ma3 = "PM2.5 3-d", pm25_ma7 = "PM2.5 7-d", pm25_ma14 = "PM2.5 14-d", pm25_ma28 = "PM2.5 28-d",
         pm10_ma7 = "PM10 7-d", no2_ma7 = "NO2 7-d", so2_ma7 = "SO2 7-d", o3_ma7 = "O3 7-d", co_ma7 = "CO 7-d", aqi_ma7 = "US AQI 7-d",
         temp_ma7 = "Temperature 7-d", rh_ma7 = "Humidity 7-d")
EXP <- EXP[names(EXP) %in% names(d)]
mk  <- MARKERS |> filter(var %in% names(d))

# ------------------------------------------------------------------ A. exposure x marker correlation
corr_tab <- function(dat, tag) {
  expand_grid(e = names(EXP), m = mk$var) |> pmap_dfr(function(e, m) {
    ok <- complete.cases(dat[[e]], dat[[m]]); if (sum(ok) < 30) return(NULL)
    ct <- suppressWarnings(cor.test(dat[[e]][ok], dat[[m]][ok], method = "spearman", exact = FALSE))
    tibble(exposure = EXP[[e]], var = m, rho = unname(ct$estimate), p = ct$p.value, n = sum(ok), dataset = tag)
  }) |> mutate(p_fdr = p.adjust(p, "BH"), label = lab_of(var))
}
cr <- bind_rows(corr_tab(d, "All months"), corr_tab(filter(d, Smog_A_WholeMonths_NovFeb == 1), "Nov-Feb only"))
write_tab(cr, "tab6a_exposure_marker_spearman")
cr <- cr |> mutate(exposure = factor(exposure, levels = EXP), label = factor(label, levels = rev(mk$label)),
                   star = case_when(p_fdr < 0.001 ~ "***", p_fdr < 0.01 ~ "**", p_fdr < 0.05 ~ "*", TRUE ~ ""))
p_a <- ggplot(cr, aes(exposure, label, fill = rho)) + geom_tile(colour = "white") +
  geom_text(aes(label = paste0(sprintf("%.2f", rho), star)), size = 2.2) +
  scale_fill_gradient2(low = COL_WARM, mid = "white", high = COL_SMOG, limits = c(-0.5, 0.5), oob = squish, name = "Spearman rho") +
  facet_wrap(~ dataset) + labs(x = NULL, y = NULL, title = "Exposure-marker correlations (* BH-FDR < 0.05, ** < 0.01, *** < 0.001)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_fig(p_a, "fig11_exposure_marker_heatmap", 190, 140, tiff = TRUE)

# ------------------------------------------------------------------ B. marker-marker matrix
cm_vars <- mk$var[mk$var %in% names(d)]
M <- cor(d[, cm_vars], use = "pairwise.complete.obs", method = "spearman"); M[is.na(M)] <- 0
ord <- hclust(as.dist(1 - abs(M)))$order; lv <- mk$label[match(cm_vars, mk$var)][ord]
cmg <- as_tibble(M, rownames = "a") |> pivot_longer(-a, names_to = "b", values_to = "rho") |>
  mutate(a = factor(lab_of(a), levels = lv), b = factor(lab_of(b), levels = lv))
p_b <- ggplot(cmg, aes(a, b, fill = rho)) + geom_tile(colour = "white") +
  scale_fill_gradient2(low = COL_WARM, mid = "white", high = COL_SMOG, limits = c(-1, 1), name = "Spearman rho") +
  labs(x = NULL, y = NULL, title = "Marker-marker correlations (hierarchically ordered)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank())
save_fig(p_b, "figS4_marker_correlation_matrix", 170, 150)

# ------------------------------------------------------------------ C. PCA
pc_vars <- c("CBC_Hemoglobin_HB_g_dL", "CBC_HCT_Percent", "CBC_Total_RBCs_10e6_uL", "CBC_MCV_fL", "CBC_WBC_TLC_10e3_uL",
             "Neut_Abs_10e3_uL", "Lymph_Abs_10e3_uL", "CBC_Platelet_Count_10e3_uL", "NLR", "Eos_Abs_10e3_uL")
pc_tr   <- c(NLR = "log", Eos_Abs_10e3_uL = "log1p", CBC_Platelet_Count_10e3_uL = "log", Neut_Abs_10e3_uL = "log", Lymph_Abs_10e3_uL = "log")
pcd <- d |> filter(CBC_Suspect %in% c(0, NA) | is.na(CBC_Suspect)) |> filter(complete.cases(across(all_of(pc_vars))))
X <- pcd |> select(all_of(pc_vars)) |> mutate(across(all_of(names(pc_tr)), ~ tf(.x, pc_tr[[cur_column()]])))
X <- X[, apply(X, 2, function(z) all(is.finite(z)))]; keep <- complete.cases(X); X <- X[keep, ]; pcd <- pcd[keep, ]
pca <- prcomp(X, scale. = TRUE); ve <- summary(pca)$importance[2, ] * 100
sc  <- bind_cols(pcd |> select(Study_ID, Season, smog, pm25_ma7, Month_Label), as_tibble(pca$x[, 1:4]))
ld  <- as_tibble(pca$rotation[, 1:2], rownames = "var") |> mutate(label = lab_of(var))
ld$label <- ifelse(is.na(ld$label), str_remove(ld$var, "CBC_|_.*$"), ld$label)
pc_tests <- tibble(PC = paste0("PC", 1:4),
  var_explained_pct = ve[1:4],
  rho_pm25_ma7 = map_dbl(1:4, ~ suppressWarnings(cor(sc[[paste0("PC", .x)]], sc$pm25_ma7, method = "spearman"))),
  p_pm25 = map_dbl(1:4, ~ suppressWarnings(cor.test(sc[[paste0("PC", .x)]], sc$pm25_ma7, method = "spearman", exact = FALSE)$p.value)),
  p_smog_vs_warm = map_dbl(1:4, ~ { s <- sc[!is.na(sc$smog), ]; suppressWarnings(wilcox.test(s[[paste0("PC", .x)]] ~ s$smog)$p.value) })) |>
  mutate(p_pm25_fdr = p.adjust(p_pm25, "BH"), p_season_fdr = p.adjust(p_smog_vs_warm, "BH"))
write_tab(pc_tests, "tab6b_pca_tests"); write_tab(as_tibble(pca$rotation, rownames = "var"), "tab6b_pca_loadings")
sc_plot <- sc |> filter(Season != "Shoulder")
mult <- min(max(abs(sc$PC1)) / max(abs(ld$PC1)), max(abs(sc$PC2)) / max(abs(ld$PC2))) * 0.8
p_c1 <- ggplot(sc_plot, aes(PC1, PC2, colour = Season)) + geom_point(alpha = 0.5, size = 1.2) + stat_ellipse(level = 0.68, linewidth = 0.8) +
  geom_segment(data = ld, aes(x = 0, y = 0, xend = PC1 * mult, yend = PC2 * mult), inherit.aes = FALSE,
               arrow = arrow(length = unit(1.5, "mm")), colour = "grey30") +
  ggrepel::geom_text_repel(data = ld, aes(PC1 * mult, PC2 * mult, label = label), inherit.aes = FALSE, size = 2.6, colour = "grey20") +
  scale_colour_manual(values = SEASON_COLS, name = NULL) +
  labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]), title = "A. PCA of the blood-count profile")
p_c2 <- ggplot(sc, aes(PC1, PC2, colour = pm25_ma7)) + geom_point(alpha = 0.7, size = 1.2) +
  scale_colour_viridis_c(option = "magma", direction = -1, name = expression(PM[2.5]~"7-d")) +
  labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]), title = "B. Same plot coloured by 7-day PM2.5")
save_fig(p_c1 + p_c2, "fig12_pca", 190, 90, tiff = TRUE)

# ------------------------------------------------------------------ D. clustering
npc <- which(cumsum(ve) >= 80)[1]; Z <- pca$x[, 1:npc, drop = FALSE]; set.seed(SEED)
sil <- map_dfr(2:5, function(kk) { km <- kmeans(Z, kk, nstart = 25); s <- cluster::silhouette(km$cluster, dist(Z)); tibble(k = kk, sil = mean(s[, 3])) })
kbest <- sil$k[which.max(sil$sil)]; km <- kmeans(Z, kbest, nstart = 50); sc$cluster <- factor(km$cluster)
write_tab(sil, "tab6c_cluster_silhouette")
Xz <- as_tibble(scale(X)) |> mutate(cluster = sc$cluster) |> group_by(cluster) |> summarise(across(everything(), mean), .groups = "drop") |>
  pivot_longer(-cluster, names_to = "var", values_to = "z") |> mutate(label = lab_of(var)); Xz$label[is.na(Xz$label)] <- str_remove(Xz$var[is.na(Xz$label)], "CBC_|_.*$")
comp <- sc |> filter(Season != "Shoulder") |> count(cluster, Season) |> group_by(cluster) |> mutate(pct = 100 * n / sum(n)) |> ungroup()
chi <- suppressWarnings(chisq.test(table(sc$cluster[sc$Season != "Shoulder"], sc$Season[sc$Season != "Shoulder"]))$p.value)
write_tab(comp, "tab6c_cluster_by_season")
p_d1 <- ggplot(Xz, aes(cluster, label, fill = z)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.1f", z)), size = 2.6) +
  scale_fill_gradient2(low = COL_WARM, mid = "white", high = COL_SMOG, name = "Mean z") + labs(x = "Cluster", y = NULL, title = "A. Cluster profiles (z-scores)")
p_d2 <- ggplot(comp, aes(cluster, pct, fill = Season)) + geom_col(width = 0.6) + geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3, colour = "white") +
  scale_fill_manual(values = SEASON_COLS) + labs(x = "Cluster", y = "% of cluster", title = "B. Season composition", subtitle = sprintf("k = %d (best silhouette); chi-square p = %s", kbest, fmt_p(chi)))
save_fig(p_d1 + p_d2, "fig13_clusters", 190, 90)

# ------------------------------------------------------------------ E. random forest: smog vs warm from the blood profile
rf_vars <- c(pc_vars, "CBC_MCH_pg", "CBC_MCHC_g_dL", "CBC_Monocytes_Percent", "Age_Years", "Sex_Female")
rd <- d |> filter(!is.na(smog)) |> select(smog, any_of(rf_vars)) |> filter(complete.cases(across(everything()))) |> mutate(smog = factor(smog, c(0, 1), c("Warm", "Smog")))
message("Random forest: n = ", nrow(rd), " (smog ", sum(rd$smog == "Smog"), ", warm ", sum(rd$smog == "Warm"), ")")
cv_auc <- function(dat, ntree = 300, folds = 5, seed = 1) {
  set.seed(seed); f <- sample(rep_len(1:folds, nrow(dat))); pr <- numeric(nrow(dat))
  for (i in 1:folds) { m <- ranger::ranger(smog ~ ., dat[f != i, ], num.trees = ntree, probability = TRUE); pr[f == i] <- predict(m, dat[f == i, ])$predictions[, "Smog"] }
  as.numeric(pROC::auc(pROC::roc(dat$smog, pr, quiet = TRUE, levels = c("Warm", "Smog"), direction = "<")))
}
auc_obs <- mean(map_dbl(1:3, ~ cv_auc(rd, seed = .x)))
set.seed(SEED); auc_null <- map_dbl(seq_len(N_PERM_RF), function(i) { r2 <- rd; r2$smog <- sample(r2$smog); cv_auc(r2, ntree = 150, seed = i) })
p_perm <- (sum(auc_null >= auc_obs) + 1) / (N_PERM_RF + 1)
rf <- ranger::ranger(smog ~ ., rd, num.trees = 1000, importance = "permutation", probability = TRUE)
extra_labs <- c(Age_Years = "Age", Sex_Female = "Female sex", CBC_MCH_pg = "MCH", CBC_MCHC_g_dL = "MCHC", CBC_Monocytes_Percent = "Monocytes %")
imp <- tibble(var = names(rf$variable.importance), imp = rf$variable.importance) |>
  mutate(label = coalesce(lab_of(var), extra_labs[var], var))
write_tab(imp |> arrange(desc(imp)), "tab6d_rf_importance"); write_tab(tibble(cv_auc = auc_obs, perm_null_mean = mean(auc_null), perm_p = p_perm, n = nrow(rd)), "tab6d_rf_auc")
p_e1 <- ggplot(imp, aes(imp, reorder(label, imp))) + geom_col(fill = COL_SMOG, width = 0.7) + labs(x = "Permutation importance", y = NULL, title = "A. What separates smog from warm?")
p_e2 <- ggplot(tibble(auc = auc_null), aes(auc)) + geom_histogram(bins = 20, fill = "grey70", colour = "white") +
  geom_vline(xintercept = auc_obs, colour = COL_SMOG, linewidth = 1) +
  labs(x = "Cross-validated AUC", y = "Permutations", title = "B. Observed AUC vs label-shuffled null",
       subtitle = sprintf("Observed AUC = %.2f; permutation p = %s", auc_obs, fmt_p(p_perm)))
save_fig(p_e1 + p_e2, "fig14_random_forest", 190, 90)
message("06 done.")