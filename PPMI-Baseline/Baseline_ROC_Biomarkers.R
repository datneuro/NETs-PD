# =============================================================================
# Baseline_ROC_Biomarkers.R
# Diagnostic ROC curves for baseline NET transcripts (Supplementary Fig. S1)
# Bootstrap (5000 iter) + DeLong CIs + a glmnet elastic-net combined panel.
# =============================================================================
# Inputs (place in ./data or pass as the first command-line argument):
#   df_normed_filtered_annotated.RData  -> NAME.log2.cpm.filtered.norm.df
#   metaDataIR3.csv
# Transcripts: PADI4-201, PADI4-202, MPO-201 | Comparison: HC vs PD at baseline
# Usage:  Rscript Baseline_ROC_Biomarkers.R [DATA_DIR] [OUT_DIR]
# Output: <OUT_DIR>/{data,pdf,png}/roc_*
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ── CLI arguments ──────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[1] else
  "data"
OUT_DIR  <- if (length(args) >= 2) args[2] else
  "results"

dir.create(file.path(OUT_DIR, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "pdf"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "png"),  recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "I_roc_log.md")
log_md   <- function(...) cat(paste0(..., "\n"), file = LOG_FILE, append = TRUE)
cat("", file = LOG_FILE)
log_md("# ROC Panel Log — Baseline_ROC_Biomarkers.R")
log_md("Date: ", format(Sys.time()))
log_md("")

set.seed(42)

# ── Constants ──────────────────────────────────────────────────────────────────
TRANSCRIPTS <- c(
  "PADI4-201" = "ENST00000375448.4",
  "PADI4-202" = "ENST00000375453.5",
  "MPO-201"   = "ENST00000225275.3"
)

ROC_COLORS <- c(
  "PADI4-201"    = "#DC267F",
  "PADI4-202"    = "#785EF0",
  "MPO-201"      = "#648FFF",
  "Combined panel" = "#FE6100"
)

BOOT_N      <- 5000
BOOT_STRAT  <- TRUE
NFOLDS_GLMNET <- 10

# ── Load data ──────────────────────────────────────────────────────────────────
message("Loading data...")
load(file.path(DATA_DIR, "Copy of df_normed_filtered_annotated.RData"))
NAME.log2.cpm.filtered.norm.df <- NAME.log2.cpm.filtered.norm.df %>%
  dplyr::select(-last_col())

samples <- read_csv(file.path(DATA_DIR, "metaDataIR3.csv"), show_col_types = FALSE) %>%
  filter(DIAGNOSIS %in% c("PD", "Control")) %>%
  mutate(sample_id = paste0(PATNO, "_", CLINICAL_EVENT),
         DIAGNOSIS = recode(DIAGNOSIS, Control = "HC"),
         DIAGNOSIS = factor(DIAGNOSIS, levels = c("HC", "PD")))

net_long <- NAME.log2.cpm.filtered.norm.df %>%
  filter(geneID %in% TRANSCRIPTS) %>%
  pivot_longer(cols = `3174_V08`:last_col(),
               names_to  = "sample_id",
               values_to = "log2cpm") %>%
  mutate(isoform = names(TRANSCRIPTS)[match(geneID, TRANSCRIPTS)])

# Baseline HC vs PD
bl_data <- samples %>%
  filter(CLINICAL_EVENT == "BL") %>%
  dplyr::select(sample_id, PATNO, DIAGNOSIS) %>%
  inner_join(net_long, by = "sample_id") %>%
  dplyr::select(sample_id, PATNO, DIAGNOSIS, isoform, log2cpm)

bl_wide <- bl_data %>%
  pivot_wider(names_from = isoform, values_from = log2cpm) %>%
  filter(!is.na(`PADI4-201`) & !is.na(`PADI4-202`) & !is.na(`MPO-201`))

log_md("## Data summary")
log_md("- Baseline samples: ", nrow(bl_wide))
log_md("- HC: ", sum(bl_wide$DIAGNOSIS == "HC"), " | PD: ", sum(bl_wide$DIAGNOSIS == "PD"))
log_md("")

# =============================================================================
# INDIVIDUAL ROC curves with bootstrap + DeLong
# =============================================================================
message("Computing individual ROC curves...")

roc_results <- list()

for (iso in names(TRANSCRIPTS)) {
  message("  ROC: ", iso)

  roc_obj <- roc(bl_wide$DIAGNOSIS, bl_wide[[iso]],
                 levels = c("HC", "PD"), direction = "<",
                 quiet = TRUE)

  # DeLong CI
  ci_delong <- ci.auc(roc_obj, method = "delong")

  # Bootstrap CI
  ci_boot <- tryCatch(
    ci.auc(roc_obj, method = "bootstrap",
            boot.n            = BOOT_N,
            boot.stratified   = BOOT_STRAT,
            progress          = "none"),
    error = function(e) { message("  bootstrap error: ", conditionMessage(e)); ci_delong }
  )

  # Youden's optimal threshold
  coords_best <- coords(roc_obj, "best", best.method = "youden",
                         ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"))

  result <- tibble(
    isoform     = iso,
    AUC         = round(auc(roc_obj), 4),
    CI_lo_boot  = round(ci_boot[1], 4),
    CI_hi_boot  = round(ci_boot[3], 4),
    CI_lo_delong = round(ci_delong[1], 4),
    CI_hi_delong = round(ci_delong[3], 4),
    threshold   = round(coords_best$threshold, 4),
    sensitivity = round(coords_best$sensitivity, 4),
    specificity = round(coords_best$specificity, 4),
    n_HC        = sum(bl_wide$DIAGNOSIS == "HC"),
    n_PD        = sum(bl_wide$DIAGNOSIS == "PD"),
    boot_n      = BOOT_N
  )

  log_md("### ", iso)
  log_md("- AUC = ", result$AUC,
         " [DeLong: ", result$CI_lo_delong, "–", result$CI_hi_delong, "]",
         " [Boot: ", result$CI_lo_boot, "–", result$CI_hi_boot, "]")
  log_md("- Youden threshold = ", result$threshold,
         " (Sens = ", result$sensitivity,
         ", Spec = ", result$specificity, ")")
  log_md("")

  roc_results[[iso]] <- list(roc = roc_obj, summary = result)
}

# =============================================================================
# COMBINED PANEL: glmnet 10-fold CV
# =============================================================================
message("Computing glmnet combined panel...")

X <- as.matrix(bl_wide[, names(TRANSCRIPTS)])
y <- ifelse(bl_wide$DIAGNOSIS == "PD", 1, 0)

cv_fit <- cv.glmnet(X, y,
                     family   = "binomial",
                     nfolds   = NFOLDS_GLMNET,
                     type.measure = "auc",
                     alpha    = 0.5,   # elastic net (LASSO + ridge)
                     standardize = TRUE)

# Cross-validated predictions at lambda.1se
cv_pred <- predict(cv_fit, newx = X, s = "lambda.1se", type = "response")

roc_panel <- roc(y, as.numeric(cv_pred),
                  levels = c(0, 1), direction = "<", quiet = TRUE)

ci_panel_delong <- ci.auc(roc_panel, method = "delong")
ci_panel_boot   <- tryCatch(
  ci.auc(roc_panel, method = "bootstrap",
          boot.n = BOOT_N, boot.stratified = BOOT_STRAT, progress = "none"),
  error = function(e) ci_panel_delong
)

panel_result <- tibble(
  isoform      = "Combined panel",
  AUC          = round(auc(roc_panel), 4),
  CI_lo_boot   = round(ci_panel_boot[1], 4),
  CI_hi_boot   = round(ci_panel_boot[3], 4),
  CI_lo_delong = round(ci_panel_delong[1], 4),
  CI_hi_delong = round(ci_panel_delong[3], 4),
  lambda       = round(cv_fit$lambda.1se, 6),
  n_HC         = sum(bl_wide$DIAGNOSIS == "HC"),
  n_PD         = sum(bl_wide$DIAGNOSIS == "PD"),
  boot_n       = BOOT_N
)

log_md("### Combined panel (PADI4-201 + PADI4-202 + MPO-201)")
log_md("- AUC = ", panel_result$AUC,
       " [DeLong: ", panel_result$CI_lo_delong, "–", panel_result$CI_hi_delong, "]")
log_md("- Lambda.1se = ", panel_result$lambda)
log_md("- Coefficients:")
coef_panel <- coef(cv_fit, s = "lambda.1se")
log_md("```")
log_md(capture.output(print(coef_panel)))
log_md("```")
log_md("")

# DeLong test
best_single <- names(which.max(sapply(roc_results, function(r) auc(r$roc))))
delong_compare <- roc.test(roc_results[[best_single]]$roc, roc_panel, method = "delong")
log_md("### DeLong test: Combined vs ", best_single)
log_md("- p = ", round(delong_compare$p.value, 4))
log_md("")

roc_results[["Combined panel"]] <- list(roc = roc_panel, summary = panel_result)

# =============================================================================
# FIGURE: ROC plot
# =============================================================================
# Build ROC curve data frames
roc_df_list <- lapply(names(roc_results), function(nm) {
  r <- roc_results[[nm]]$roc
  tibble(
    FPR      = 1 - r$specificities,
    TPR      = r$sensitivities,
    isoform  = nm,
    AUC      = round(auc(r), 3),
    is_panel = nm == "Combined panel"
  )
})
roc_df_all <- bind_rows(roc_df_list)

# AUC labels
auc_labels <- sapply(names(roc_results), function(nm) {
  res <- roc_results[[nm]]$summary
  if (nm == "Combined panel") {
    sprintf("%s\nAUC=%.3f [%.3f–%.3f] (DeLong)",
            nm, res$AUC, res$CI_lo_delong, res$CI_hi_delong)
  } else {
    sprintf("%s\nAUC=%.3f [%.3f–%.3f]",
            nm, res$AUC, res$CI_lo_boot, res$CI_hi_boot)
  }
})

fig_roc <- ggplot(roc_df_all, aes(x = FPR, y = TPR,
                                   color = isoform, linetype = isoform)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey70", linewidth = 0.6) +
  geom_path(linewidth = 1.1) +
  scale_color_manual(
    values = ROC_COLORS,
    labels = auc_labels,
    name   = NULL
  ) +
  scale_linetype_manual(
    values = c("PADI4-201" = "solid", "PADI4-202" = "solid",
               "MPO-201"   = "solid", "Combined panel" = "longdash"),
    labels = auc_labels,
    name   = NULL
  ) +
  scale_x_continuous(labels = percent_format(), limits = c(0, 1),
                     name = "1 − Specificity (False Positive Rate)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1),
                     name = "Sensitivity (True Positive Rate)") +
  labs(
    title    = "ROC curves: NET transcripts for PD vs HC (baseline)",
    subtitle = paste0("Bootstrap CI: ", BOOT_N, " iterations, stratified | ",
                      "Combined panel: glmnet elastic net (α=0.5), 10-fold CV\n",
                      "DeLong combined vs best single p = ",
                      round(delong_compare$p.value, 4))
  ) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(size = 9, face = "italic", lineheight = 1.2),
    legend.position = c(0.65, 0.25),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.text     = element_text(size = 9),
    axis.title      = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUT_DIR, "pdf", "roc_panel.pdf"),
       fig_roc, width = 8, height = 8)
ggsave(file.path(OUT_DIR, "png", "roc_panel.png"),
       fig_roc, width = 8, height = 8, dpi = 300)

# Save AUC
auc_table <- bind_rows(lapply(names(roc_results), function(nm) {
  roc_results[[nm]]$summary
}))
write_csv(auc_table, file.path(OUT_DIR, "data", "I_auc_summary.csv"))

# Save glmnet model
saveRDS(cv_fit, file.path(OUT_DIR, "data", "I_glmnet_panel.rds"))

log_md("## Output files")
log_md("- data/I_auc_summary.csv")
log_md("- data/I_glmnet_panel.rds")
log_md("- pdf/roc_panel.pdf")
log_md("- png/roc_panel.png  (300 dpi)")

message("Done — Baseline_ROC_Biomarkers.R")
