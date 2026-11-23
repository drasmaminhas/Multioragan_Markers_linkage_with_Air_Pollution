# =====================================================================================
# 02_descriptives.R  -  STROBE-style flow, Table 1 (by season and by PM2.5 tertile),
#                       sampling-fraction and test-availability figures
# =====================================================================================
source("00_setup.R")
d_all <- load_analysis(adults_only = FALSE)
d     <- load_analysis()

# ------------------------------------------------------------------ Figure 1: flow diagram
n_raw   <- sum(d_all$Records_Merged); n_uniq <- nrow(d_all); n_child <- n_uniq - nrow(d)
s_n     <- table(d$Season)
panel_n <- c(CBC = sum(d$Has_CBC == 1), `Renal` = sum(d$Has_RFT == 1), Electrolytes = sum(d$Has_Electrolytes == 1),
             LFT = sum(d$Has_LFT == 1), Glucose = sum(d$Has_Glucose == 1), Troponin = sum(d$Has_Troponin == 1),
             Virology = sum(d$Has_Virology == 1))
box <- function(x, y, txt, w = 62, h = 14, fill = "white")
  list(annotate("rect", xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill, colour = "black"),
       annotate("text", x = x, y = y, label = txt, size = 3))
arr <- function(x, y, xe, ye) annotate("segment", x = x, y = y, xend = xe, yend = ye, arrow = arrow(length = unit(2, "mm")))
pf <- ggplot() + xlim(0, 200) + ylim(0, 112) + theme_void() +
  box(70, 100, sprintf("Lab reports extracted from PDFs\nn = %d", n_raw)) +
  box(150, 86, sprintf("Duplicate extractions removed\nn = %d", n_raw - n_uniq), w = 56, fill = "grey95") +
  box(70, 72, sprintf("Unique lab reports\nn = %d", n_uniq)) +
  box(150, 58, sprintf("Excluded from adult analyses:\nage < 18, n = %d", n_child), w = 56, fill = "grey95") +
  box(70, 44, sprintf("Analysis set\nn = %d", nrow(d)), fill = "grey90") +
  box(25, 20, sprintf("Smog\n15 Nov - 15 Feb\nn = %d", s_n["Smog"]), w = 44, h = 18, fill = "#F6D5C0") +
  box(75, 20, sprintf("Shoulder\n1-14 Nov, 16-28 Feb\nn = %d", s_n["Shoulder"]), w = 44, h = 18, fill = "grey90") +
  box(125, 20, sprintf("Warm\n1 Mar - end Jul\nn = %d", s_n["Warm"]), w = 44, h = 18, fill = "#C9DFF0") +
  box(100, 3, paste(sprintf("%s %d", names(panel_n), panel_n), collapse = "  |  "), w = 196, h = 8) +
  arr(70, 93, 70, 79) + arr(70, 65, 70, 51) + arr(70, 86, 122, 86) + arr(70, 58, 122, 58) + arr(70, 37, 25, 29) +
  arr(70, 37, 75, 29) + arr(70, 37, 125, 29) +
  labs(title = "Study flow and available tests (patients with each panel)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_fig(pf, "fig1_flow", 180, 130)

# ------------------------------------------------------------------ Table 1 builder
cont_vars <- c("Age_Years", MARKERS$var[MARKERS$var %in% names(d)])
cont_lab  <- c("Age (years)", MARKERS$label[MARKERS$var %in% names(d)])
bin_vars  <- c("Sex_Female", BINARY$var[BINARY$var %in% names(d)])
bin_lab   <- c("Female sex", BINARY$label[BINARY$var %in% names(d)])

make_table1 <- function(data, g) {
  data <- data[!is.na(data[[g]]), ]; lv <- levels(factor(data[[g]]))
  hdr <- tibble(Variable = "N patients", p = NA_real_)
  for (l in lv) hdr[[l]] <- as.character(sum(data[[g]] == l))
  rows_c <- map2_dfr(cont_vars, cont_lab, function(v, lb) {
    dig <- if (v == "Troponin_I_hs_pg_mL") 2 else 1
    r <- tibble(Variable = paste0(lb, ", median (IQR)"))
    for (l in lv) { x <- data[[v]][data[[g]] == l]; r[[l]] <- sprintf("%s [n=%d]", med_iqr(x, dig), sum(!is.na(x))) }
    ok <- !is.na(data[[v]]); r$p <- if (length(unique(data[[g]][ok])) > 1) kruskal.test(data[[v]][ok], factor(data[[g]][ok]))$p.value else NA
    r })
  rows_b <- map2_dfr(bin_vars, bin_lab, function(v, lb) {
    r <- tibble(Variable = paste0(lb, ", n (%)"))
    for (l in lv) { x <- data[[v]][data[[g]] == l]; r[[l]] <- sprintf("%s [n=%d]", n_pct(x), sum(!is.na(x))) }
    ok <- !is.na(data[[v]]); r$p <- cat_p(data[[v]][ok], data[[g]][ok]); r })
  out <- bind_rows(hdr, rows_c, rows_b)
  out |> mutate(p_fdr = p.adjust(p, "BH"), p = fmt_p(p), p_fdr = fmt_p(p_fdr)) |> relocate(p, p_fdr, .after = last_col())
}
d$PM25_tertile <- cut(d$pm25_ma7, quantile(d$pm25_ma7, c(0, 1/3, 2/3, 1), na.rm = TRUE), include.lowest = TRUE,
                      labels = c("T1 (lowest PM2.5)", "T2", "T3 (highest PM2.5)"))
t1_season <- make_table1(d, "Season"); t1_pm <- make_table1(d, "PM25_tertile")
write_tab(t1_season, "table1_by_season"); write_tab(t1_pm, "table1_by_pm25_tertile")
tq <- d |> group_by(PM25_tertile) |> summarise(min = min(pm25_ma7), max = max(pm25_ma7), n = n()); write_tab(tq, "pm25_tertile_cutpoints")

if (requireNamespace("flextable", quietly = TRUE) && requireNamespace("officer", quietly = TRUE)) {
  doc <- officer::read_docx()
  for (nm in c("Table 1A. Characteristics by season (Kruskal-Wallis / chi-square; BH-adjusted p)",
               "Table 1B. Characteristics by tertile of 7-day mean PM2.5")) {
    tb <- if (grepl("1A", nm)) t1_season else t1_pm
    doc <- officer::body_add_par(doc, nm, style = "heading 2")
    ft <- flextable::autofit(flextable::fontsize(flextable::font(flextable::flextable(tb), fontname = "Arial"), size = 8, part = "all"))
    doc <- flextable::body_add_flextable(doc, ft)
  }
  print(doc, target = file.path(TAB_DIR, "table1.docx")); message("  table saved: table1.docx")
}

# ------------------------------------------------------------------ Figure S2: sampling fraction + test availability
mon <- d |> group_by(Month_Label) |> summarise(n = n(), reg = first(Month_Registrations_Est), .groups = "drop") |>
  mutate(frac = 100 * n / reg)
pA <- ggplot(mon, aes(Month_Label, n)) + geom_col(fill = "grey40", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%% of\nregistrations", frac)), vjust = -0.2, size = 2.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + labs(x = NULL, y = "Lab reports (n)", title = "A. Reports per month and sampling fraction")
av <- d |> group_by(Month_Label) |>
  summarise(across(c(Has_CBC, Has_RFT, Has_Electrolytes, Has_LFT, Has_Glucose, Has_Troponin, Has_Virology), ~ 100 * mean(.x == 1)), .groups = "drop") |>
  pivot_longer(-Month_Label, names_to = "Panel", values_to = "pct") |> mutate(Panel = str_remove(Panel, "Has_"))
pB <- ggplot(av, aes(Month_Label, Panel, fill = pct)) + geom_tile(colour = "white") +
  geom_text(aes(label = round(pct)), size = 2.8) + scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 100), name = "% of reports") +
  labs(x = NULL, y = NULL, title = "B. Percentage of reports with each test panel")
save_fig(pA / pB, "figS2_sampling_and_availability", 180, 150)

# ------------------------------------------------------------------ Figure S3: missingness map (rows = patients in time order)
mv <- c("Troponin_I_hs_pg_mL", "Blood_Glucose_Random_mg_dL", "RFT_Creatinine_mg_dL", "Electrolytes_Sodium_mmol_L",
        "LFT_SGPT_ALT_U_L", "CBC_Hemoglobin_HB_g_dL", "CBC_WBC_TLC_10e3_uL", "CBC_Platelet_Count_10e3_uL", "NLR", "HBsAg_Pos")
short <- c(Troponin_I_hs_pg_mL = "Troponin (hs)", Blood_Glucose_Random_mg_dL = "Glucose", RFT_Creatinine_mg_dL = "Creatinine",
           Electrolytes_Sodium_mmol_L = "Sodium", LFT_SGPT_ALT_U_L = "ALT", CBC_Hemoglobin_HB_g_dL = "Haemoglobin",
           CBC_WBC_TLC_10e3_uL = "WBC", CBC_Platelet_Count_10e3_uL = "Platelets", NLR = "NLR", HBsAg_Pos = "HBsAg")
mm <- d |> arrange(Est_Date) |> mutate(row = row_number()) |>
  transmute(row, across(any_of(mv), ~ !is.na(.x))) |>
  pivot_longer(-row, names_to = "var", values_to = "present") |> mutate(var = factor(short[var], levels = unname(short)))
pM <- ggplot(mm, aes(var, row, fill = present)) + geom_raster() + scale_y_reverse() +
  scale_fill_manual(values = c(`TRUE` = "grey25", `FALSE` = "grey90"), labels = c("missing", "present"), name = NULL) +
  labs(x = NULL, y = "Patients (chronological order)", title = "Missingness pattern") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(pM, "figS3_missingness_map", 150, 130)
message("02 done.")