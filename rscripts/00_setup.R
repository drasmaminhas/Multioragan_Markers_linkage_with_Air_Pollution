# =====================================================================================
# 00_setup.R  -  shared settings, marker definitions, helper functions, plot theme
# Every other script starts with source("00_setup.R").  Run scripts from THIS folder.
# =====================================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(purrr)
  library(stringr); library(forcats); library(tibble); library(lubridate)
  library(scales); library(patchwork); library(sandwich); library(mgcv)
})

# ------------------------------------------------------------------ PATHS (edit once)
PROJECT_DIR <- if (.Platform$OS.type == "windows") "C:/Users/Home/Desktop/Asma_Asghar/Mphil" else Sys.getenv("MPHIL_PROJECT", "/home/szd/win_desktop/Asma_Asghar/Mphil")

# Auto-locate cleaned clinical CSV inside claude_cleaned_tables/output/ (or root fallback)
DATA_CANDIDATES <- c(
  file.path(PROJECT_DIR, "claude_cleaned_tables", "output", "analysis_ready_data.csv"),
  file.path(PROJECT_DIR, "claude_cleaned_tables", "analysis_ready_data.csv")
)
DATA_CSV <- DATA_CANDIDATES[file.exists(DATA_CANDIDATES)][1]
if (is.na(DATA_CSV)) DATA_CSV <- file.path(PROJECT_DIR, "claude_cleaned_tables", "output", "analysis_ready_data.csv")

AQI_CSV <- file.path(PROJECT_DIR, "AQI", "okara_consensus_aqi_2026.csv")

# Auto-locate daily weather CSV across AQI/, analysis/results/data/, or gemini_pipeline_v1/data/
WEATHER_CANDIDATES <- c(
  file.path(PROJECT_DIR, "AQI", "okara_weather_daily.csv"),
  file.path(PROJECT_DIR, "AQI", "okara_daily_weather_2025_2026.csv"),
  file.path(PROJECT_DIR, "analysis", "results", "data", "okara_weather_daily.csv"),
  file.path(PROJECT_DIR, "gemini_pipeline_v1", "data", "okara_daily_weather_2025_2026.csv")
)
WEATHER_CSV <- WEATHER_CANDIDATES[file.exists(WEATHER_CANDIDATES)][1]
if (is.na(WEATHER_CSV)) WEATHER_CSV <- file.path(PROJECT_DIR, "analysis", "results", "data", "okara_weather_daily.csv")

# Isolated output subfolder inside analysis/
OUT_DIR <- file.path(PROJECT_DIR, "analysis", "results")

# ------------------------------------------------------------------ ANALYSIS SETTINGS
PRIMARY_SMOG <- "Smog_B_Nov15_Feb15"  # "Smog_A_WholeMonths_NovFeb" | "Smog_B_Nov15_Feb15" | "Smog_C_Nov15_Jan31"
PM_UNIT      <- 50                    # effects are reported per +50 ug/m3 PM2.5
LAG_MAX      <- 14                    # distributed-lag window (days)
TROP_HS_UNIT_VERIFIED <- TRUE         # verified 0.00 - 0.06 pg/mL (Abbott Architect Plus i1000SR)
N_SIM_JITTER <- 200                   # date-jitter simulations (script 07)
N_PERM_RF    <- 100                   # label permutations for the random-forest check (script 06)
FONT         <- "sans"                # e.g. "Arial" if installed
SEED         <- 20260101
set.seed(SEED)

TAB_DIR <- file.path(OUT_DIR, "tables"); FIG_DIR <- file.path(OUT_DIR, "figures"); DAT_DIR <- file.path(OUT_DIR, "data")
for (p in c(TAB_DIR, FIG_DIR, DAT_DIR)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------ COLOURS / THEME
COL_SMOG <- "#D55E00"; COL_WARM <- "#0072B2"; COL_SHOULDER <- "#999999"
OKABE <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
SEASON_COLS <- c(Smog = COL_SMOG, Shoulder = COL_SHOULDER, Warm = COL_WARM)

theme_paper <- function(base = 10) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
          legend.position = "bottom", plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
}
theme_set(theme_paper())

save_fig <- function(p, name, w = 180, h = 120, tiff = FALSE) {
  f <- file.path(FIG_DIR, name)
  pdf_dev <- tryCatch({ grDevices::cairo_pdf; grDevices::cairo_pdf }, error = function(e) grDevices::pdf)
  ggsave(paste0(f, ".pdf"), p, width = w, height = h, units = "mm", device = pdf_dev)
  ggsave(paste0(f, ".png"), p, width = w, height = h, units = "mm", dpi = 300, device = ragg::agg_png)
  if (tiff) ggsave(paste0(f, ".tiff"), p, width = w, height = h, units = "mm", dpi = 600,
                   device = ragg::agg_tiff, compression = "lzw")
  message("  figure saved: ", name)
}
write_tab <- function(df, name) { write_csv(df, file.path(TAB_DIR, paste0(name, ".csv"))); message("  table saved: ", name) }

# ------------------------------------------------------------------ MARKER DEFINITIONS
# tr = transformation used for modelling: none | log | log1p
MARKERS <- tribble(
  ~var,                          ~label,                 ~system,           ~tr,     ~primary,
  "NLR",                         "NLR",                  "Inflammation",    "log",   TRUE,
  "SII",                         "SII",                  "Inflammation",    "log",   FALSE,
  "PLR",                         "PLR",                  "Inflammation",    "log",   FALSE,
  "CBC_WBC_TLC_10e3_uL",         "WBC",                  "Inflammation",    "none",  TRUE,
  "Neut_Abs_10e3_uL",            "Neutrophils (abs)",    "Inflammation",    "none",  FALSE,
  "Lymph_Abs_10e3_uL",           "Lymphocytes (abs)",    "Inflammation",    "none",  FALSE,
  "CBC_Platelet_Count_10e3_uL",  "Platelets",            "Inflammation",    "none",  TRUE,
  "Eos_Abs_10e3_uL",             "Eosinophils (abs)",    "Airway / allergy","log1p", FALSE,
  "CBC_Hemoglobin_HB_g_dL",      "Haemoglobin",          "Oxygen / blood",  "none",  TRUE,
  "CBC_HCT_Percent",             "Haematocrit",          "Oxygen / blood",  "none",  FALSE,
  "CBC_Total_RBCs_10e6_uL",      "RBC count",            "Oxygen / blood",  "none",  FALSE,
  "CBC_MCV_fL",                  "MCV",                  "Oxygen / blood",  "none",  FALSE,
  "eGFR_CKDEPI2021",             "eGFR",                 "Renal",           "none",  TRUE,
  "RFT_Creatinine_mg_dL",        "Creatinine",           "Renal",           "log",   FALSE,
  "RFT_Urea_mg_dL",              "Urea",                 "Renal",           "log",   FALSE,
  "Electrolytes_Sodium_mmol_L",  "Sodium",               "Electrolytes",    "none",  FALSE,
  "Electrolytes_Potassium_mmol_L","Potassium",           "Electrolytes",    "none",  FALSE,
  "Electrolytes_Chloride_mmol_L","Chloride",             "Electrolytes",    "none",  FALSE,
  "Blood_Glucose_Random_mg_dL",  "Random glucose",       "Metabolic",       "log",   FALSE,
  "LFT_SGPT_ALT_U_L",            "ALT",                  "Liver",           "log",   FALSE,
  "LFT_SGOT_AST_U_L",            "AST",                  "Liver",           "log",   FALSE,
  "LFT_ALP_U_L",                 "ALP",                  "Liver",           "log",   FALSE,
  "LFT_Bilirubin_Total_mg_dL",   "Total bilirubin",      "Liver",           "log",   FALSE
)
if (TROP_HS_UNIT_VERIFIED)
  MARKERS <- bind_rows(tibble(var = "Troponin_I_hs_pg_mL", label = "hs-Troponin I", system = "Cardiac",
                              tr = "log", primary = TRUE), MARKERS)

BINARY <- tribble(
  ~var,                          ~label,                       ~system,           ~primary,
  "Trop_Elevated_QL",            "Troponin positive (qual.)",  "Cardiac",         TRUE,
  "Anaemia",                     "Anaemia",                    "Oxygen / blood",  FALSE,
  "Erythrocytosis",              "Erythrocytosis",             "Oxygen / blood",  FALSE,
  "Leukocytosis",                "Leukocytosis",               "Inflammation",    FALSE,
  "Neutrophilia_ANC",            "Neutrophilia",               "Inflammation",    FALSE,
  "Lymphopenia_ALC",             "Lymphopenia",                "Inflammation",    FALSE,
  "Thrombocytosis",              "Thrombocytosis",             "Inflammation",    FALSE,
  "Eosinophilia_AEC",            "Eosinophilia",               "Airway / allergy",FALSE,
  "eGFR_lt60",                   "eGFR < 60",                  "Renal",           FALSE,
  "Hyponatraemia",               "Hyponatraemia",              "Electrolytes",    FALSE,
  "Hypokalaemia",                "Hypokalaemia",               "Electrolytes",    FALSE,
  "Hyperkalaemia",               "Hyperkalaemia",              "Electrolytes",    FALSE,
  "Hyperglycaemia_Random_ge200", "Random glucose >= 200",      "Metabolic",       FALSE,
  "ALT_High",                    "ALT high",                   "Liver",           FALSE
)
NEG_CONTROL <- c(HBsAg_Pos = "HBsAg positive", AntiHCV_Pos = "Anti-HCV positive", AntiHIV_Pos = "Anti-HIV positive")
KEY6 <- c("NLR", "CBC_WBC_TLC_10e3_uL", "CBC_Hemoglobin_HB_g_dL", "CBC_Platelet_Count_10e3_uL",
          "eGFR_CKDEPI2021", "Blood_Glucose_Random_mg_dL")
RENAL_VARS <- c("eGFR_CKDEPI2021", "RFT_Creatinine_mg_dL", "RFT_Urea_mg_dL")

lab_of <- function(v) { m <- bind_rows(MARKERS, BINARY); m$label[match(v, m$var)] }
tf     <- function(x, tr) switch(tr, none = x, log = log(x), log1p = log1p(x))
tr_of  <- function(v) MARKERS$tr[match(v, MARKERS$var)]

# ------------------------------------------------------------------ LOAD ANALYSIS FRAME
load_analysis <- function(path = file.path(DAT_DIR, "analysis_exposure.csv"), adults_only = TRUE) {
  if (!file.exists(path)) stop("Run 01_merge_exposure.R first (missing ", path, ")")
  d <- read_csv(path, show_col_types = FALSE, guess_max = 20000)
  lv <- d |> distinct(Month_Label, Month_Index) |> arrange(Month_Index) |> pull(Month_Label)
  d <- d |> mutate(
    Est_Date    = as.Date(Est_Date),
    Season      = factor(Season_B_Label, levels = c("Smog", "Shoulder", "Warm")),
    Month_Label = factor(Month_Label, levels = lv),
    Month_Season= factor(if_else(Lab_Month %in% c(11, 12, 1, 2), "Smog months (Nov-Feb)", "Warm months (Mar-Jul)")),
    Sex         = factor(if_else(Sex_Female == 1, "Female", "Male"), levels = c("Male", "Female")),
    Age_Group   = factor(Age_Group, levels = c("<40", "40-49", "50-59", "60-69", "70-79", "80+")),
    Age60       = factor(if_else(Age_Years >= 60, "60+", "<60"), levels = c("<60", "60+")),
    smog        = .data[[PRIMARY_SMOG]],
    pm25_ma7_50 = pm25_ma7 / PM_UNIT, pm25_ma14_50 = pm25_ma14 / PM_UNIT,
    pm25_ma3_50 = pm25_ma3 / PM_UNIT, pm25_ma28_50 = pm25_ma28 / PM_UNIT,
    pm25_L0_50  = pm25_L0 / PM_UNIT, pm25_lead7_50 = pm25_lead7 / PM_UNIT,
    cons_pm25_ma7_50 = cons_pm25_ma7 / PM_UNIT, pm10_ma7_50 = pm10_ma7 / PM_UNIT
  )
  if (adults_only) d <- d |> filter(is.na(Adult_Flag) | Adult_Flag == 1)
  d
}
has_weather <- function(d) all(c("temp_ma7", "rh_ma7") %in% names(d)) && sum(!is.na(d$temp_ma7)) > 0
covars_for <- function(model, d, outcome = NULL) {
  base <- c("Age_Group", "Sex")
  switch(model,
    M0 = character(0),
    M1 = base,
    M2 = if (has_weather(d)) c(base, "temp_ma7", "rh_ma7") else NULL,
    M3 = if (!is.null(outcome) && outcome %in% RENAL_VARS) NULL
         else c(base, if (has_weather(d)) c("temp_ma7", "rh_ma7"), "eGFR_CKDEPI2021"))
}

# ------------------------------------------------------------------ STAT HELPERS
fmt_p   <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
med_iqr <- function(x, dig = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  q <- quantile(x, c(.25, .5, .75)); sprintf(paste0("%.", dig, "f (%.", dig, "f-%.", dig, "f)"), q[2], q[1], q[3])
}
n_pct <- function(x) { x <- x[!is.na(x)]; if (!length(x)) "" else sprintf("%d (%.1f%%)", sum(x == 1), 100 * mean(x == 1)) }
cohens_d <- function(y1, y0) {
  y1 <- y1[!is.na(y1)]; y0 <- y0[!is.na(y0)]; n1 <- length(y1); n0 <- length(y0)
  sp <- sqrt(((n1 - 1) * var(y1) + (n0 - 1) * var(y0)) / (n1 + n0 - 2)); dd <- (mean(y1) - mean(y0)) / sp
  se <- sqrt((n1 + n0) / (n1 * n0) + dd^2 / (2 * (n1 + n0)))
  c(d = dd, lo = dd - 1.96 * se, hi = dd + 1.96 * se)
}
rank_biserial <- function(x1, x0) {
  w <- suppressWarnings(wilcox.test(x1, x0, exact = FALSE))$statistic
  unname(2 * w / (length(x1) * length(x0)) - 1)
}
cat_p <- function(x, g) {
  tab <- table(x, g); tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  exp_ok <- all(suppressWarnings(chisq.test(tab)$expected) >= 5)
  if (exp_ok) chisq.test(tab)$p.value else fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value
}
wilson <- function(x, n) { if (n == 0) return(c(NA_real_, NA_real_)); prop.test(x, n, correct = FALSE)$conf.int }

# Robust (HC3) linear model: effect of one term, with standardised effect (per SD of outcome)
fit_cont <- function(d, var, tr, x, covars, label_extra = NULL) {
  cols <- c(var, x, covars); dd <- d[complete.cases(d[, cols]), ]
  if (nrow(dd) < 30 || length(unique(dd[[x]])) < 2) return(NULL)
  dd$.y <- tf(dd[[var]], tr); dd <- dd[is.finite(dd$.y), ]
  fit <- tryCatch(lm(reformulate(c(x, covars), response = ".y"), data = dd), error = function(e) NULL)
  if (is.null(fit) || is.na(coef(fit)[x])) return(NULL)
  V <- sandwich::vcovHC(fit, type = "HC3"); b <- unname(coef(fit)[x]); se <- sqrt(V[x, x]); dfr <- fit$df.residual
  sdy <- sd(dd$.y)
  tibble(var = var, n = nrow(dd), est = b, lo = b - qt(.975, dfr) * se, hi = b + qt(.975, dfr) * se,
         p = 2 * pt(-abs(b / se), dfr), sd_y = sdy, std = b / sdy, std_lo = (b - qt(.975, dfr) * se) / sdy,
         std_hi = (b + qt(.975, dfr) * se) / sdy,
         pct = if (tr %in% c("log", "log1p")) (exp(b) - 1) * 100 else NA_real_,
         pct_lo = if (tr %in% c("log", "log1p")) (exp(b - qt(.975, dfr) * se) - 1) * 100 else NA_real_,
         pct_hi = if (tr %in% c("log", "log1p")) (exp(b + qt(.975, dfr) * se) - 1) * 100 else NA_real_)
}
# Binary outcome: Firth logistic if {logistf} installed, otherwise ordinary logistic (Wald)
fit_bin <- function(d, var, x, covars) {
  cols <- c(var, x, covars); dd <- d[complete.cases(d[, cols]), ]
  ev <- sum(dd[[var]] == 1); n <- nrow(dd)
  if (n < 30 || ev < 3 || ev > n - 3 || length(unique(dd[[x]])) < 2) return(NULL)
  f <- reformulate(c(x, covars), response = var)
  out <- tryCatch({
    if (requireNamespace("logistf", quietly = TRUE)) {
      fit <- logistf::logistf(f, data = dd)
      tibble(b = unname(fit$coefficients[x]), lo = unname(fit$ci.lower[x]), hi = unname(fit$ci.upper[x]),
             p = unname(fit$prob[x]), method = "Firth")
    } else {
      fit <- glm(f, data = dd, family = binomial); s <- summary(fit)$coefficients[x, ]
      tibble(b = s[1], lo = s[1] - 1.96 * s[2], hi = s[1] + 1.96 * s[2], p = s[4], method = "Wald")
    }
  }, error = function(e) NULL)
  if (is.null(out)) return(NULL)
  out |> mutate(var = var, n = n, events = ev, or = exp(b), or_lo = exp(lo), or_hi = exp(hi))
}
# E-value for a standardised mean difference (VanderWeele & Ding) via RR ~ exp(0.91 d)
evalue_from_d <- function(d) { rr <- exp(0.91 * abs(d)); rr + sqrt(rr * (rr - 1)) }

# Distributed-lag basis (natural spline over lags 0..LAG_MAX; linear in exposure)
dlm_basis <- function(df_lag = 3) splines::ns(0:LAG_MAX, df = df_lag, intercept = TRUE)   # intercept=TRUE: lag-0 effect is estimated, not forced to 0
LAG_COLS  <- paste0("pm25_L", 0:LAG_MAX)

message("00_setup.R loaded. Project: ", PROJECT_DIR)
message("  Clinical CSV: ", DATA_CSV, " (exists: ", file.exists(DATA_CSV), ")")
message("  Outputs isolated to: ", OUT_DIR)
