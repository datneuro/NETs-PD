# =============================================================================
# Baseline_Clinical_Correlations.R
# DaT scan + UPDRS clinical correlations (Supplementary Fig. S3)
# Robust (percentage-bend) + longitudinal LMM + Bayesian (brms) estimation.
# =============================================================================
# Data:
#   df_normed_filtered_annotated.RData  -> NAME.log2.cpm.filtered.norm.df
#   metaDataIR3.csv
#   PPMI_Curated_Data_Cut_July_2024.xlsx -> con_caudate, con_putamen, con_striatum
#   PAR.UPDRS.df.RData -> PAR.UPDRS (PATNO, EVENT_ID, ENROLL_AGE, HY_Stage, Part1-3, TOTAL_UPDRS)
# Usage:  Rscript Baseline_Clinical_Correlations.R [DATA_DIR] [OUT_DIR] [PPMI_DIR]
# Output: <OUT_DIR>/{data,pdf,png}/Iplus_*
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(sandwich)
  library(lmtest)
  library(WRS2)        # pbcor robust correlation
  library(boot)
  library(brms)
  library(bayestestR)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(corrplot)
})

# ── CLI arguments ──────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[1] else
  "data"
PPMI_DIR <- if (length(args) >= 3) args[3] else
  "data"
OUT_DIR  <- if (length(args) >= 2) args[2] else
  "results"

dir.create(file.path(OUT_DIR, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "pdf"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "png"),  recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "Iplus_correlations_log.md")
log_md   <- function(...) cat(paste0(..., "\n"), file = LOG_FILE, append = TRUE)
cat("", file = LOG_FILE)
log_md("# DaT/UPDRS Correlation Log — Baseline_Clinical_Correlations.R")
log_md("Date: ", format(Sys.time()))
log_md("")

set.seed(42)

TRANSCRIPTS <- c(
  "PADI4-201" = "ENST00000375448.4",
  "PADI4-202" = "ENST00000375453.5",
  "MPO-201"   = "ENST00000225275.3"
)

BRMS_CHAINS <- 4
BRMS_ITER   <- 4000
BRMS_WARMUP <- 2000
BRMS_SEED   <- 42
RUN_BRMS    <- TRUE

custom_colors <- c("HC" = "#648FFF", "Prodromal" = "#785EF0", "PD" = "#DC267F")

# =============================================================================
# PART 1: Load and harmonize all data sources
# =============================================================================
message("Loading expression data...")
load(file.path(DATA_DIR, "Copy of df_normed_filtered_annotated.RData"))
NAME.log2.cpm.filtered.norm.df <- NAME.log2.cpm.filtered.norm.df %>%
  dplyr::select(-last_col())

samples <- read_csv(file.path(DATA_DIR, "metaDataIR3.csv"), show_col_types = FALSE) %>%
  filter(DIAGNOSIS %in% c("PD", "Prodromal", "Control")) %>%
  mutate(sample_id = paste0(PATNO, "_", CLINICAL_EVENT),
         DIAGNOSIS = recode(DIAGNOSIS, Control = "HC"),
         DIAGNOSIS = factor(DIAGNOSIS, levels = c("HC", "Prodromal", "PD")),
         timepoint_num = case_when(
           CLINICAL_EVENT == "BL"  ~ 0,
           CLINICAL_EVENT == "V02" ~ 6,
           CLINICAL_EVENT == "V04" ~ 12,
           CLINICAL_EVENT == "V06" ~ 24,
           TRUE                    ~ 36
         ),
         sex = factor(GENDER))

net_long <- NAME.log2.cpm.filtered.norm.df %>%
  filter(geneID %in% TRANSCRIPTS) %>%
  pivot_longer(cols = `3174_V08`:last_col(),
               names_to  = "sample_id",
               values_to = "log2cpm") %>%
  mutate(isoform = names(TRANSCRIPTS)[match(geneID, TRANSCRIPTS)])

expr_data <- samples %>%
  dplyr::select(sample_id, PATNO, CLINICAL_EVENT, DIAGNOSIS, sex,
                `RIN Value`, Plate, timepoint_num) %>%
  inner_join(net_long, by = "sample_id") %>%
  pivot_wider(names_from = isoform, values_from = log2cpm)

# ── UPDRS ──────────────────────────────────────────────────────────────────────
message("Loading UPDRS data...")
updrs_path <- file.path(DATA_DIR, "PAR.UPDRS.df.RData")
if (file.exists(updrs_path)) {
  load(updrs_path)
  PAR.UPDRS <- PAR.UPDRS %>%
    mutate(
      sample_id     = paste0(PATNO, "_", EVENT_ID),
      CLINICAL_EVENT = EVENT_ID,
      Duration_Days  = as.numeric(INFODT - PD_Diagnosis_Date),
      Duration_Years = Duration_Days / 365.25
    ) %>%
    dplyr::select(PATNO, sample_id, CLINICAL_EVENT,
                  ENROLL_AGE, HY_Stage, Part1, Part2, Part3, TOTAL_UPDRS,
                  Duration_Years)
} else {
  log_md("WARNING: PAR.UPDRS.df.RData not found — UPDRS analysis skipped")
  PAR.UPDRS <- NULL
}

# ── DaT scan ───────────────────────────────────────────────────────────────────
message("Loading DaT scan data...")
curated_path <- file.path(PPMI_DIR,
                           "PPMI_Curated_Data_Cut_July_2024.xlsx")
if (file.exists(curated_path)) {
  curated <- readxl::read_xlsx(curated_path)
  curated$sample_id <- paste0(curated$PATNO, "_", curated$EVENT_ID)
  dat_data <- curated %>%
    dplyr::select(PATNO, sample_id,
                  any_of(c("con_caudate", "con_putamen", "con_striatum",
                            "COHORT", "subgroup")))
} else {
  log_md("WARNING: PPMI Curated data not found — DaT analysis skipped")
  dat_data <- NULL
}

# ── Merged dataset ─────────────────────────────────────────────────────────────
combined_data <- expr_data
if (!is.null(PAR.UPDRS)) {
  combined_data <- combined_data %>%
    left_join(dplyr::select(PAR.UPDRS, -CLINICAL_EVENT), by = c("sample_id", "PATNO"))
}
if (!is.null(dat_data)) {
  combined_data <- combined_data %>%
    left_join(dplyr::select(dat_data, -any_of("sample_id")),
              by = "PATNO")
  # If DaT has per-visit rows, join by sample_id
  if ("sample_id" %in% names(dat_data)) {
    dat_visit <- dat_data %>% dplyr::select(sample_id, starts_with("con_"))
    combined_data <- combined_data %>%
      left_join(dat_visit, by = "sample_id", suffix = c("", ".visit"))
    for (col in c("con_caudate", "con_putamen", "con_striatum")) {
      visit_col <- paste0(col, ".visit")
      if (visit_col %in% names(combined_data)) {
        combined_data[[col]] <- coalesce(combined_data[[paste0(col, ".visit")]],
                                          combined_data[[col]])
        combined_data[[visit_col]] <- NULL
      }
    }
  }
}

log_md("## Data summary")
log_md("- Combined dataset rows: ", nrow(combined_data))
log_md("- Subjects: ", n_distinct(combined_data$PATNO))
log_md("")

# =============================================================================
# PART 2: Cross-sectional robust correlations (baseline)
# =============================================================================
message("Cross-sectional robust correlations...")
log_md("## Cross-sectional robust correlations (baseline, PD only)")

bl_pd <- combined_data %>%
  filter(CLINICAL_EVENT == "BL", DIAGNOSIS == "PD")

corr_pairs <- list(
  list(x = "PADI4-201", y = "TOTAL_UPDRS", label = "PADI4-201 vs Total UPDRS"),
  list(x = "PADI4-202", y = "TOTAL_UPDRS", label = "PADI4-202 vs Total UPDRS"),
  list(x = "MPO-201",   y = "TOTAL_UPDRS", label = "MPO-201 vs Total UPDRS"),
  list(x = "PADI4-201", y = "Part3",        label = "PADI4-201 vs MDS-UPDRS Part 3"),
  list(x = "PADI4-202", y = "Part3",        label = "PADI4-202 vs MDS-UPDRS Part 3"),
  list(x = "MPO-201",   y = "Part3",        label = "MPO-201 vs MDS-UPDRS Part 3"),
  list(x = "PADI4-201", y = "con_striatum", label = "PADI4-201 vs DaT striatum"),
  list(x = "PADI4-202", y = "con_striatum", label = "PADI4-202 vs DaT striatum"),
  list(x = "MPO-201",   y = "con_striatum", label = "MPO-201 vs DaT striatum"),
  list(x = "PADI4-201", y = "con_caudate",  label = "PADI4-201 vs DaT caudate"),
  list(x = "PADI4-201", y = "con_putamen",  label = "PADI4-201 vs DaT putamen")
)

cross_results <- list()

for (pair in corr_pairs) {
  x_var <- pair$x; y_var <- pair$y; lbl <- pair$label

  df_pair <- bl_pd %>%
    dplyr::select(all_of(c(x_var, y_var, "PATNO"))) %>%
    rename(x = all_of(x_var), y = all_of(y_var)) %>%
    filter(!is.na(x), !is.na(y))

  if (nrow(df_pair) < 10) {
    log_md("### ", lbl, ": skipped (n = ", nrow(df_pair), ")")
    next
  }

  # WRS2 percentage-bend robust correlation
  pb_res <- tryCatch(
    WRS2::pbcor(df_pair$x, df_pair$y),
    error = function(e) NULL
  )

  # Spearman as sensitivity
  sp_res <- cor.test(df_pair$x, df_pair$y, method = "spearman", exact = FALSE)

  # Bootstrap CI for Pearson r
  boot_fun <- function(d, i) cor(d$x[i], d$y[i], method = "pearson")
  boot_ci <- tryCatch({
    b <- boot::boot(df_pair, boot_fun, R = 5000)
    boot::boot.ci(b, type = "perc")$percent[4:5]
  }, error = function(e) c(NA, NA))

  result <- tibble(
    label       = lbl,
    x_var       = x_var,
    y_var       = y_var,
    n           = nrow(df_pair),
    pbcor_r     = if (!is.null(pb_res)) round(pb_res$cor, 4)    else NA,
    pbcor_p     = if (!is.null(pb_res)) round(pb_res$p.value, 4) else NA,
    spearman_r  = round(sp_res$estimate, 4),
    spearman_p  = round(sp_res$p.value, 4),
    pearson_ci_lo = round(boot_ci[1], 4),
    pearson_ci_hi = round(boot_ci[2], 4)
  )

  log_md("### ", lbl)
  log_md("- n = ", result$n)
  log_md("- pbcor: r = ", result$pbcor_r, ", p = ", result$pbcor_p)
  log_md("- Spearman: r = ", result$spearman_r, ", p = ", result$spearman_p)
  log_md("- Pearson bootstrap 95% CI: [", result$pearson_ci_lo,
         ", ", result$pearson_ci_hi, "]")
  log_md("")

  cross_results[[lbl]] <- result
}

cross_table <- bind_rows(cross_results)
write_csv(cross_table, file.path(OUT_DIR, "data", "Iplus_crosssectional_corr.csv"))

# =============================================================================
# PART 3: Scatter plots with robust smooth
# =============================================================================
make_scatter <- function(df_pair, x_label, y_label, title_text, corr_row) {
  subtitle <- if (!is.na(corr_row$pbcor_r)) {
    sprintf("pbcor r=%.3f, p=%.4f | Spearman r=%.3f, p=%.4f | n=%d",
            corr_row$pbcor_r, corr_row$pbcor_p,
            corr_row$spearman_r, corr_row$spearman_p,
            corr_row$n)
  } else {
    sprintf("Spearman r=%.3f, p=%.4f | n=%d",
            corr_row$spearman_r, corr_row$spearman_p, corr_row$n)
  }

  ggplot(df_pair, aes(x = x, y = y)) +
    geom_point(color = "#DC267F", alpha = 0.5, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "#DC267F",
                fill = "#DC267F", alpha = 0.12) +
    geom_smooth(method = "loess", se = FALSE, color = "grey40",
                linetype = "dashed", linewidth = 0.8) +
    labs(title    = title_text,
         subtitle = subtitle,
         x = x_label, y = y_label) +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8, face = "italic"))
}

# Key scatter plots
scatter_plots <- list()
key_pairs <- list(
  list("PADI4-201", "TOTAL_UPDRS", "PADI4-201 mRNA (log2 CPM)", "Total MDS-UPDRS"),
  list("PADI4-201", "con_striatum", "PADI4-201 mRNA (log2 CPM)", "DaT SBR (contralateral striatum)")
)

for (kp in key_pairs) {
  x_var  <- kp[[1]]; y_var <- kp[[2]]
  x_lbl  <- kp[[3]]; y_lbl <- kp[[4]]
  title  <- paste0(x_var, " vs ", y_var)
  key    <- paste0(x_var, " vs ", y_var)

  df_pair <- bl_pd %>%
    dplyr::select(all_of(c(x_var, y_var))) %>%
    rename(x = all_of(x_var), y = all_of(y_var)) %>%
    filter(!is.na(x), !is.na(y))

  if (nrow(df_pair) >= 5) {
    corr_row <- filter(cross_table, x_var == !!x_var, y_var == !!y_var)
    if (nrow(corr_row) == 0) {
      sp  <- cor.test(df_pair$x, df_pair$y, method = "spearman", exact = FALSE)
      corr_row <- tibble(pbcor_r = NA, pbcor_p = NA,
                          spearman_r = round(sp$estimate, 4),
                          spearman_p = round(sp$p.value, 4),
                          n = nrow(df_pair))
    }
    scatter_plots[[key]] <- make_scatter(df_pair, x_lbl, y_lbl, title, corr_row[1, ])
  }
}

if (length(scatter_plots) >= 2) {
  fig_scatter <- scatter_plots[[1]] + scatter_plots[[2]] +
    plot_annotation(tag_levels = "A")
  ggsave(file.path(OUT_DIR, "pdf", "Iplus_scatter_key.pdf"),
         fig_scatter, width = 12, height = 5)
  ggsave(file.path(OUT_DIR, "png", "Iplus_scatter_key.png"),
         fig_scatter, width = 12, height = 5, dpi = 300)
}

# =============================================================================
# PART 4: Longitudinal LMM — DaT and UPDRS trajectories
# =============================================================================
message("Longitudinal LMMs...")
log_md("## Longitudinal LMM")

longit_results <- list()

for (outcome in c("TOTAL_UPDRS", "con_striatum", "con_caudate")) {
  for (predictor in names(TRANSCRIPTS)) {
    lm_data <- combined_data %>%
      filter(DIAGNOSIS == "PD") %>%
      dplyr::select(PATNO, timepoint_num, DIAGNOSIS, sex,
                    ENROLL_AGE, all_of(c(predictor, outcome))) %>%
      rename(predictor_expr = all_of(predictor),
             outcome_val    = all_of(outcome)) %>%
      filter(!is.na(predictor_expr), !is.na(outcome_val))

    if (nrow(lm_data) < 20 || n_distinct(lm_data$PATNO) < 5) {
      log_md("### ", predictor, " → ", outcome, ": skipped (n = ", nrow(lm_data), ")")
      next
    }

    key_lm <- paste0(predictor, "_", outcome)
    log_md("### ", predictor, " → ", outcome,
           " (n = ", nrow(lm_data), ", n_subj = ", n_distinct(lm_data$PATNO), ")")

    lmm_fit <- tryCatch(
      lmer(outcome_val ~ predictor_expr + timepoint_num + ENROLL_AGE + sex +
             (1 | PATNO),
           data    = lm_data,
           REML    = TRUE,
           control = lmerControl(optimizer = "bobyqa",
                                  optCtrl   = list(maxfun = 2e5))),
      error = function(e) NULL
    )

    if (!is.null(lmm_fit)) {
      coef_tbl <- as.data.frame(coef(summary(lmm_fit)))

      # HC3 robust SE
      robust_lmm <- tryCatch(
        coeftest(lmm_fit, vcov. = vcovHC(lmm_fit, type = "HC3")),
        error = function(e) NULL
      )

      log_md("#### LMM fixed effects:")
      log_md("```"); log_md(capture.output(print(coef_tbl))); log_md("```")

      if (!is.null(robust_lmm)) {
        log_md("#### HC3 robust SE:")
        log_md("```"); log_md(capture.output(print(robust_lmm))); log_md("```")
      }

      longit_results[[key_lm]] <- list(
        model = lmm_fit,
        coef  = coef_tbl,
        robust = robust_lmm,
        n      = nrow(lm_data),
        n_subj = n_distinct(lm_data$PATNO)
      )
    }
    log_md("")
  }
}

# =============================================================================
# PART 5: Bayesian correlation models
# =============================================================================
if (RUN_BRMS) {
  message("Bayesian models with brms...")
  log_md("## Bayesian models (brms)")

  brms_corr_results <- list()

  key_outcomes <- intersect(
    c("TOTAL_UPDRS", "con_striatum"),
    names(combined_data)
  )

  for (outcome in key_outcomes) {
    for (predictor in c("PADI4-201", "PADI4-202")) {
      bdata <- combined_data %>%
        filter(DIAGNOSIS == "PD") %>%
        dplyr::select(PATNO, ENROLL_AGE, sex,
                      all_of(c(predictor, outcome))) %>%
        rename(predictor_expr = all_of(predictor),
               outcome_val    = all_of(outcome)) %>%
        filter(!is.na(predictor_expr), !is.na(outcome_val))

      if (nrow(bdata) < 15 || n_distinct(bdata$PATNO) < 5) next

      key_b <- paste0(predictor, "_", outcome, "_bayes")
      message("  brms: ", predictor, " → ", outcome,
              " (n = ", nrow(bdata), ")")
      log_md("### Bayesian: ", predictor, " → ", outcome)

      # Center predictors for stable sampling
      bdata$pred_z   <- scale(bdata$predictor_expr)[, 1]
      bdata$age_z    <- scale(bdata$ENROLL_AGE)[, 1]
      y_sd           <- sd(bdata$outcome_val, na.rm = TRUE)

      brms_fit <- tryCatch(
        brm(
          outcome_val ~ pred_z + age_z + sex + (1 | PATNO),
          data   = bdata,
          family = gaussian(),
          prior  = c(
            prior(normal(0, 1),           class = "b"),
            prior(normal(0, 5),           class = "Intercept"),
            prior(student_t(3, 0, 1),     class = "sd"),
            prior(student_t(3, 0, y_sd),  class = "sigma")
          ),
          chains  = BRMS_CHAINS,
          iter    = BRMS_ITER,
          warmup  = BRMS_WARMUP,
          seed    = BRMS_SEED,
          cores   = 4,
          refresh = 0,
          silent  = 2
        ),
        error = function(e) { message("  brms error: ", conditionMessage(e)); NULL }
      )

      if (!is.null(brms_fit)) {
        summ <- summary(brms_fit)
        pred_row <- summ$fixed["pred_z", ]

        log_md("- beta (standardized) = ", round(pred_row["Estimate"], 4),
               " [", round(pred_row["l-95% CI"], 4),
               ", ", round(pred_row["u-95% CI"], 4), "]")

        post <- as_draws_df(brms_fit)
        pp   <- mean(post$b_pred_z > 0)
        log_md("- P(effect > 0) = ", round(pp, 4))
        log_md("  (P < 0.5 = negative association is more probable than positive)")
        log_md("")

        saveRDS(brms_fit,
                file.path(OUT_DIR, "data", paste0("Iplus_brms_", key_b, ".rds")))
        brms_corr_results[[key_b]] <- list(
          fit        = brms_fit,
          beta       = pred_row["Estimate"],
          ci_lo      = pred_row["l-95% CI"],
          ci_hi      = pred_row["u-95% CI"],
          post_prob  = pp
        )
      }
    }
  }
}

# =============================================================================
# PART 6: Correlation matrix figure (baseline PD)
# =============================================================================
message("Generating correlation matrix figure...")

corr_vars <- intersect(
  c("PADI4-201", "PADI4-202", "MPO-201",
    "TOTAL_UPDRS", "Part1", "Part2", "Part3", "HY_Stage",
    "con_caudate", "con_putamen", "con_striatum",
    "Duration_Years", "ENROLL_AGE"),
  names(combined_data)
)

corr_bl_pd <- combined_data %>%
  filter(CLINICAL_EVENT == "BL", DIAGNOSIS == "PD") %>%
  dplyr::select(all_of(corr_vars)) %>%
  filter(rowSums(is.na(.)) < length(corr_vars) * 0.5)

if (nrow(corr_bl_pd) >= 10 && ncol(corr_bl_pd) >= 3) {
  corr_mat <- cor(corr_bl_pd, use = "pairwise.complete.obs", method = "spearman")

  pdf(file.path(OUT_DIR, "pdf", "Iplus_corrplot_PD_baseline.pdf"),
      width = 10, height = 10)
  corrplot(corr_mat,
           method   = "color",
           type     = "upper",
           tl.cex   = 0.9,
           tl.col   = "black",
           addCoef.col = "black",
           number.cex  = 0.65,
           col       = colorRampPalette(c("#648FFF", "white", "#DC267F"))(200),
           title    = "Spearman correlations: PD baseline",
           mar      = c(0, 0, 2, 0))
  dev.off()

  write_csv(as.data.frame(corr_mat) %>% rownames_to_column("var"),
            file.path(OUT_DIR, "data", "Iplus_spearman_matrix_PD_bl.csv"))
}

# =============================================================================
# Summary table
# =============================================================================
longit_table <- bind_rows(lapply(names(longit_results), function(k) {
  res <- longit_results[[k]]
  coef_tbl <- as.data.frame(res$coef) %>%
    rownames_to_column("term") %>%
    filter(grepl("predictor_expr", term))
  coef_tbl$analysis <- k
  coef_tbl$n        <- res$n
  coef_tbl$n_subj   <- res$n_subj
  coef_tbl
}))

if (nrow(longit_table) > 0) {
  write_csv(longit_table, file.path(OUT_DIR, "data", "Iplus_longit_lmm.csv"))
}

log_md("")
log_md("## Output files")
log_md("- data/Iplus_crosssectional_corr.csv")
log_md("- data/Iplus_longit_lmm.csv")
log_md("- data/Iplus_spearman_matrix_PD_bl.csv")
log_md("- data/Iplus_brms_*.rds")
log_md("- pdf/Iplus_scatter_key.pdf")
log_md("- pdf/Iplus_corrplot_PD_baseline.pdf")

message("Done — Baseline_Clinical_Correlations.R")
