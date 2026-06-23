# =============================================================================
# Baseline_Neutrophil_Deconvolution.R
# Neutrophil deconvolution + per-cell PADI4 activation (Supplementary Fig. S2)
# Tests whether the systemic PADI4 signal reflects neutrophil abundance or a
# per-cell change in expression (addresses NET-signature specificity).
# =============================================================================
# Inputs (place in ./data or pass as the first command-line argument):
#   CIBERSORT_ALL_Result_Dec30.RData -> `result` (BL samples x 22 LM22 cell types)
#   PPMI_Blood_Chemistry_Hematology.csv -> clinical neutrophil % (HMT15, EVENT_ID SC)
#   df_normed_filtered_annotated.RData, metaDataIR3.csv
#
# CONFIG: NEUT_SOURCE switches the neutrophil source used for adjustment
#   "deconvolution" — CIBERSORTx LM22 neutrophil fraction
#   "lab_blood"     — PPMI clinical hematology neutrophil %
#   "both"          — run both, side-by-side (default)
# Usage:  Rscript Baseline_Neutrophil_Deconvolution.R [DATA_DIR] [OUT_DIR]
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(limma)
  library(edgeR)
  library(emmeans)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
})


# ── CONFIG ─────────────────────────────────────────────────────────────────────
# Set to "deconvolution", "lab_blood", or "both"
NEUT_SOURCE <- "both"

# ── CLI arguments ──────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[1] else
  "data"
OUT_DIR  <- if (length(args) >= 2) args[2] else
  "results"

# Path to archived CIBERSORT results (already run Dec 2024)
CIBER_RDATA <- file.path("data", "CIBERSORT_ALL_Result_Dec30.RData")

# Path to PPMI blood chemistry hematology data
BLOOD_LAB_CSV <- file.path("data", "PPMI_Blood_Chemistry_Hematology.csv")

dir.create(file.path(OUT_DIR, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "pdf"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "png"),  recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "DA_deconvolution_log.md")
log_md   <- function(...) cat(paste0(..., "\n"), file = LOG_FILE, append = TRUE)
cat("", file = LOG_FILE)
log_md("# Deconvolution Analysis Log — Baseline_Neutrophil_Deconvolution.R")
log_md("Date: ", format(Sys.time()))
log_md("NEUT_SOURCE: ", NEUT_SOURCE)
log_md("")

# ── Transcript IDs of interest ─────────────────────────────────────────────────
TRANSCRIPTS <- c(
  "PADI4-201" = "ENST00000375448.4",
  "PADI4-202" = "ENST00000375453.5",
  "MPO-201"   = "ENST00000225275.3"
)
PADI4_GENE <- "PADI4"   # for gene-level deconvolution

custom_colors <- c("HC" = "#648FFF", "Prodromal" = "#785EF0", "PD" = "#DC267F")

# =============================================================================
# PART 1: Load existing CIBERSORTx results (already run Dec 2024)
# =============================================================================
message("Part 1: Loading pre-computed CIBERSORTx results...")
log_md("## Part 1: CIBERSORTx results (pre-computed Dec 2024)")

if (!file.exists(CIBER_RDATA)) {
  stop("CIBERSORT RData not found at: ", CIBER_RDATA,
       "\nThis file contains previously-run CIBERSORTx LM22 results (582 BL samples).")
}

# Load: creates object `result` — matrix 582 × 25 (22 cell types + P-value, Correlation, RMSE)
load(CIBER_RDATA)
ciber_matrix <- result
rm(result)

log_md("- Loaded: ", nrow(ciber_matrix), " samples × ", ncol(ciber_matrix), " columns")
log_md("- Source: Copy of CIBERSORT_ALL_Result_Dec30.RData")
log_md("- Cell types: ", paste(colnames(ciber_matrix)[1:22], collapse = ", "))

# Convert to data frame with PATNO
ciber_df <- as.data.frame(ciber_matrix) %>%
  rownames_to_column("PATNO") %>%
  mutate(PATNO = as.numeric(PATNO))

# Extract neutrophil fraction (CIBERSORTx LM22 has single "Neutrophils" column)
ciber_neut <- ciber_df %>%
  dplyr::select(PATNO, neutrophil_frac_deconv = Neutrophils)

message("  CIBERSORTx loaded: ", nrow(ciber_neut), " BL samples")
log_md("- Neutrophil fraction range: [",
       round(min(ciber_neut$neutrophil_frac_deconv), 4), ", ",
       round(max(ciber_neut$neutrophil_frac_deconv), 4), "]")
log_md("")

# =============================================================================
# PART 1b: Load lab blood neutrophil % (PPMI hematology)
# =============================================================================
message("Part 1b: Loading lab blood neutrophil data...")
log_md("## Part 1b: Lab blood neutrophil % (PPMI hematology)")

if (!file.exists(BLOOD_LAB_CSV)) {
  warning("Blood Chemistry CSV not found at: ", BLOOD_LAB_CSV)
  log_md("- WARNING: Blood lab CSV not found — lab_blood mode unavailable")
  lab_neut <- NULL
} else {
  blood_raw <- read_csv(BLOOD_LAB_CSV, show_col_types = FALSE)

  # HMT15 = Neutrophils (%), at Screening visit (SC)
  lab_neut <- blood_raw %>%
    filter(LTSTCODE == "HMT15", EVENT_ID == "SC") %>%
    dplyr::select(PATNO, LSIRES) %>%
    mutate(
      PATNO = as.numeric(PATNO),
      neu_blood_pct = as.numeric(LSIRES)
    ) %>%
    filter(!is.na(neu_blood_pct)) %>%
    distinct(PATNO, .keep_all = TRUE) %>%
    dplyr::select(PATNO, neu_blood_pct)

  message("  Lab neutrophil data loaded: ", nrow(lab_neut), " subjects")
  log_md("- Loaded: ", nrow(lab_neut), " subjects with Neutrophils (%) at SC")
  log_md("- Source: Blood Chemistry Hematology Nov 15 2024 (1).csv")
  log_md("- LTSTCODE = HMT15, EVENT_ID = SC")
  log_md("- Neutrophils (%) range: [",
         round(min(lab_neut$neu_blood_pct, na.rm = TRUE), 1), ", ",
         round(max(lab_neut$neu_blood_pct, na.rm = TRUE), 1), "]")

  # Overlap with CIBERSORT
  overlap_n <- sum(ciber_neut$PATNO %in% lab_neut$PATNO)
  log_md("- Overlap with CIBERSORTx PATNOs: ", overlap_n, " / ", nrow(ciber_neut))
  log_md("")
}

# =============================================================================
# PART 1c: Deconvolution vs Lab-blood validation scatter
# =============================================================================
if (!is.null(lab_neut)) {
  message("Part 1c: Deconvolution vs Lab-blood validation...")
  log_md("## Part 1c: Deconvolution vs Lab-blood validation")

  validation_df <- inner_join(ciber_neut, lab_neut, by = "PATNO")
  cor_test <- cor.test(validation_df$neutrophil_frac_deconv,
                       validation_df$neu_blood_pct,
                       method = "pearson")
  spearman_test <- cor.test(validation_df$neutrophil_frac_deconv,
                            validation_df$neu_blood_pct,
                            method = "spearman", exact = FALSE)

  log_md("- n matched: ", nrow(validation_df))
  log_md("- Pearson r = ", round(cor_test$estimate, 3),
         ", p = ", format.pval(cor_test$p.value, digits = 3))
  log_md("- Spearman rho = ", round(spearman_test$estimate, 3),
         ", p = ", format.pval(spearman_test$p.value, digits = 3))

  fig_validation <- ggplot(validation_df,
                           aes(x = neutrophil_frac_deconv, y = neu_blood_pct)) +
    geom_point(alpha = 0.4, size = 2, color = "#4A90D9") +
    geom_smooth(method = "lm", se = TRUE, color = "#DC267F", alpha = 0.15) +
    labs(
      title    = "CIBERSORTx Deconvolution vs Blood Test Neutrophils",
      subtitle = paste0("Pearson r = ", round(cor_test$estimate, 3),
                        ", p = ", format.pval(cor_test$p.value, digits = 3),
                        " | Spearman ρ = ", round(spearman_test$estimate, 3)),
      x = "Neutrophil fraction (CIBERSORTx LM22 deconvolution)",
      y = "Neutrophils % (Blood test, SC visit)"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(OUT_DIR, "pdf", "DA_validation_deconv_vs_lab.pdf"),
         fig_validation, width = 7, height = 6)
  ggsave(file.path(OUT_DIR, "png", "DA_validation_deconv_vs_lab.png"),
         fig_validation, width = 7, height = 6, dpi = 300)
  write_csv(validation_df, file.path(OUT_DIR, "data", "DA_validation_deconv_vs_lab.csv"))
  log_md("")
}

# =============================================================================
# PART 1d: Build metadata + neutrophil data for downstream
# =============================================================================
message("Part 1d: Merging metadata with neutrophil sources...")

samples <- read_csv(file.path(DATA_DIR, "metaDataIR3.csv"), show_col_types = FALSE) %>%
  filter(DIAGNOSIS %in% c("PD", "Prodromal", "Control")) %>%
  mutate(
    sample_id     = paste0(PATNO, "_", CLINICAL_EVENT),
    DIAGNOSIS     = recode(DIAGNOSIS, Control = "HC"),
    DIAGNOSIS     = factor(DIAGNOSIS, levels = c("HC", "Prodromal", "PD")),
    timepoint_num = case_when(
      CLINICAL_EVENT == "BL"  ~ 0,
      CLINICAL_EVENT == "V02" ~ 6,
      CLINICAL_EVENT == "V04" ~ 12,
      CLINICAL_EVENT == "V06" ~ 24,
      TRUE                    ~ 36
    )
  )

# Join neutrophil data — keep both sources if available
ciber_meta <- samples %>%
  filter(CLINICAL_EVENT == "BL") %>%
  left_join(ciber_neut, by = "PATNO")

if (!is.null(lab_neut)) {
  ciber_meta <- ciber_meta %>%
    left_join(lab_neut, by = "PATNO")
}

log_md("## Merged metadata (BL only)")
log_md("- n BL samples with CIBERSORTx: ", sum(!is.na(ciber_meta$neutrophil_frac_deconv)))
if (!is.null(lab_neut)) {
  log_md("- n BL samples with lab blood:  ", sum(!is.na(ciber_meta$neu_blood_pct)))
}

# Create the working neutrophil_fraction column based on NEUT_SOURCE
# For "both" mode, we run analyses on each source separately (handled in loops below)
if (NEUT_SOURCE == "lab_blood" && !is.null(lab_neut)) {
  ciber_meta$neutrophil_fraction <- ciber_meta$neu_blood_pct / 100  # convert % to fraction
} else {
  ciber_meta$neutrophil_fraction <- ciber_meta$neutrophil_frac_deconv
}

# Neutrophil fraction by diagnosis (for primary source)
neut_summ <- ciber_meta %>%
  filter(!is.na(neutrophil_fraction)) %>%
  group_by(DIAGNOSIS) %>%
  summarise(
    median_neut = round(median(neutrophil_fraction), 4),
    IQR_neut    = round(IQR(neutrophil_fraction), 4),
    n           = n(),
    .groups = "drop"
  )
log_md("### Neutrophil fractions by diagnosis (primary source):")
log_md("```"); log_md(capture.output(print(neut_summ))); log_md("```")

# Kruskal-Wallis
kw_neut <- kruskal.test(neutrophil_fraction ~ DIAGNOSIS,
                        data = ciber_meta %>% filter(!is.na(neutrophil_fraction)))
log_md("Kruskal-Wallis p = ", round(kw_neut$p.value, 4))

# Figure: neutrophil fraction violin
fig_neut <- ggplot(ciber_meta %>% filter(!is.na(neutrophil_fraction)),
                   aes(x = DIAGNOSIS, y = neutrophil_fraction,
                       fill = DIAGNOSIS, color = DIAGNOSIS)) +
  geom_violin(alpha = 0.3, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.07, alpha = 0.25, size = 1.5) +
  scale_fill_manual(values  = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title    = "Estimated neutrophil fraction (CIBERSORTx LM22)",
       subtitle = paste0("Kruskal-Wallis p = ", round(kw_neut$p.value, 4)),
       x = "", y = "Neutrophil fraction") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

ggsave(file.path(OUT_DIR, "pdf", "DA_neutrophil_violin_deconv.pdf"),
       fig_neut, width = 6, height = 5)
ggsave(file.path(OUT_DIR, "png", "DA_neutrophil_violin_deconv.png"),
       fig_neut, width = 6, height = 5, dpi = 300)

# If lab data available, also make lab version
if (!is.null(lab_neut)) {
  ciber_meta_lab <- ciber_meta %>% filter(!is.na(neu_blood_pct))
  kw_lab <- kruskal.test(neu_blood_pct ~ DIAGNOSIS, data = ciber_meta_lab)

  fig_neut_lab <- ggplot(ciber_meta_lab,
                         aes(x = DIAGNOSIS, y = neu_blood_pct,
                             fill = DIAGNOSIS, color = DIAGNOSIS)) +
    geom_violin(alpha = 0.3, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.07, alpha = 0.25, size = 1.5) +
    scale_fill_manual(values  = custom_colors) +
    scale_color_manual(values = custom_colors) +
    labs(title    = "Neutrophils % (Blood test, SC visit)",
         subtitle = paste0("Kruskal-Wallis p = ", round(kw_lab$p.value, 4)),
         x = "", y = "Neutrophils (%)") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")

  ggsave(file.path(OUT_DIR, "pdf", "DA_neutrophil_violin_labblood.pdf"),
         fig_neut_lab, width = 6, height = 5)
  ggsave(file.path(OUT_DIR, "png", "DA_neutrophil_violin_labblood.png"),
         fig_neut_lab, width = 6, height = 5, dpi = 300)

  log_md("### Lab blood neutrophils by diagnosis:")
  lab_summ <- ciber_meta_lab %>%
    group_by(DIAGNOSIS) %>%
    summarise(median_pct = round(median(neu_blood_pct), 1),
              IQR_pct = round(IQR(neu_blood_pct), 1),
              n = n(), .groups = "drop")
  log_md("```"); log_md(capture.output(print(lab_summ))); log_md("```")
  log_md("Kruskal-Wallis p = ", round(kw_lab$p.value, 4))
}

# Save merged data
write_csv(ciber_meta, file.path(OUT_DIR, "data", "DA_neutrophil_fractions.csv"))
# Save full 22-cell-type deconvolution results with metadata
ciber_full <- ciber_df %>%
  inner_join(samples %>% filter(CLINICAL_EVENT == "BL") %>%
               dplyr::select(PATNO, DIAGNOSIS, GENDER), by = "PATNO")
write_csv(ciber_full, file.path(OUT_DIR, "data", "DA_CIBERSORTx_full_22celltypes.csv"))
log_md("")

# =============================================================================
# HELPER: Run per-cell analysis for a given neutrophil source
# =============================================================================
# This function is called once or twice depending on NEUT_SOURCE
run_percell_analysis <- function(data_with_neut, neut_col_name, source_label,
                                 source_tag, net_long_data, out_dir) {

  message("  Running per-cell PADI4 analysis [", source_label, "]...")

  # Baseline samples only
  bl_expr <- net_long_data %>%
    inner_join(data_with_neut %>%
                 dplyr::select(sample_id, PATNO, CLINICAL_EVENT, DIAGNOSIS,
                               GENDER, `RIN Value`, Plate,
                               all_of(neut_col_name)),
               by = "sample_id") %>%
    filter(CLINICAL_EVENT == "BL")

  # Rename working neutrophil column
  bl_expr$neutrophil_fraction <- bl_expr[[neut_col_name]]
  bl_expr <- bl_expr %>% filter(!is.na(neutrophil_fraction))

  log_md("### Per-cell PADI4 analysis [", source_label, "]")
  log_md("- n samples: ", n_distinct(bl_expr$sample_id))

  per_cell_results <- list()

  for (iso in names(TRANSCRIPTS)) {
    iso_pc <- bl_expr %>% filter(isoform == iso)

    # Scatter: expression vs neutrophil fraction per group
    fig_scatter <- ggplot(iso_pc, aes(x = neutrophil_fraction, y = log2cpm,
                                      color = DIAGNOSIS)) +
      geom_point(alpha = 0.5, size = 2) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
      scale_color_manual(values = custom_colors) +
      labs(title    = paste0(iso, " vs Neutrophil Fraction [", source_label, "]"),
           subtitle = "Baseline samples; slope indicates per-cell association",
           x        = paste0("Neutrophil fraction (", source_label, ")"),
           y        = bquote(.(iso) ~ "mRNA (log"["2"] * "CPM)")) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    # Residual analysis: expression ~ neutrophil_fraction; test residuals HC vs PD
    iso_pc_hcpd <- iso_pc %>% filter(DIAGNOSIS %in% c("HC", "PD"))

    if (nrow(iso_pc_hcpd) < 10) {
      log_md("#### ", iso, " — SKIPPED (n < 10 HC+PD)")
      next
    }

    resid_model <- lm(log2cpm ~ neutrophil_fraction + GENDER + `RIN Value`,
                      data = iso_pc_hcpd)
    iso_pc_hcpd <- iso_pc_hcpd %>%
      mutate(resid_PADI4 = residuals(resid_model))

    wilcox_resid <- wilcox.test(resid_PADI4 ~ DIAGNOSIS,
                                data = iso_pc_hcpd, exact = FALSE)

    log_md("#### ", iso)
    log_md("- Residual Wilcoxon (HC vs PD) p = ", round(wilcox_resid$p.value, 4))

    # Save residuals CSV
    write_csv(
      iso_pc_hcpd %>%
        dplyr::select(sample_id, PATNO, DIAGNOSIS, neutrophil_fraction,
                      log2cpm, resid_PADI4) %>%
        mutate(isoform = iso, source = source_tag),
      file.path(out_dir, "data", paste0("DA_residuals_", iso, "_", source_tag, ".csv"))
    )
    log_md("- Residuals CSV: DA_residuals_", iso, "_", source_tag, ".csv")

    # Robust SE on full model
    full_model <- lm(log2cpm ~ DIAGNOSIS + neutrophil_fraction + GENDER + `RIN Value`,
                     data = iso_pc_hcpd)
    robust_full <- tryCatch(
      coeftest(full_model, vcov. = vcovHC(full_model, type = "HC3")),
      error = function(e) NULL
    )
    if (!is.null(robust_full)) {
      log_md("##### Full model (DIAGNOSIS + neutrophil_fraction) HC3 robust SE:")
      log_md("```")
      log_md(capture.output(print(robust_full)))
      log_md("```")
    }

    # Figure: residual violin per group
    fig_resid <- ggplot(iso_pc_hcpd, aes(x = DIAGNOSIS, y = resid_PADI4,
                                          fill = DIAGNOSIS, color = DIAGNOSIS)) +
      geom_violin(alpha = 0.3, trim = FALSE) +
      geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.07, alpha = 0.3, size = 1.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      scale_fill_manual(values  = c("HC" = "#648FFF", "PD" = "#DC267F")) +
      scale_color_manual(values = c("HC" = "#648FFF", "PD" = "#DC267F")) +
      labs(title    = paste0(iso, " residual expression [", source_label, "]"),
           subtitle = paste0("Wilcoxon p = ", round(wilcox_resid$p.value, 4),
                             " — residual elevation = per-cell activation"),
           y = "Residual log2 CPM", x = "") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), legend.position = "none")

    combined <- fig_scatter + fig_resid + plot_annotation(tag_levels = "A")
    ggsave(file.path(out_dir, "pdf", paste0("DA_percell_", iso, "_", source_tag, ".pdf")),
           combined, width = 12, height = 5)
    ggsave(file.path(out_dir, "png", paste0("DA_percell_", iso, "_", source_tag, ".png")),
           combined, width = 12, height = 5, dpi = 300)

    per_cell_results[[iso]] <- list(
      wilcox_p    = wilcox_resid$p.value,
      full_model  = full_model,
      robust_coef = robust_full,
      source      = source_tag
    )
  }

  per_cell_results
}

# =============================================================================
# HELPER: Isoform-level sensitivity — with vs without neutrophil adjustment
# =============================================================================
# Uses isoform-level log2CPM (already loaded from df_normed_filtered_annotated.RData)
# instead of gene-level raw counts (Txi_gene, 2.4 GB — not needed).
# For each target isoform, compares DIAGNOSIS coefficient ± neutrophil covariate.
run_isoform_sensitivity <- function(meta_bl, neut_col_name, source_label,
                                     source_tag, net_long_data, out_dir) {

  message("  Running isoform-level sensitivity [", source_label, "]...")
  log_md("### Isoform sensitivity: ± neutrophil adjustment [", source_label, "]")

  bl_samples <- meta_bl %>%
    filter(!is.na(.data[[neut_col_name]])) %>%
    dplyr::select(sample_id, PATNO, DIAGNOSIS, GENDER, `RIN Value`, Plate,
                  all_of(neut_col_name))
  bl_samples$neutrophil_fraction <- bl_samples[[neut_col_name]]

  # Join expression data (BL only, HC + PD only)
  bl_expr <- net_long_data %>%
    inner_join(bl_samples, by = "sample_id") %>%
    filter(DIAGNOSIS %in% c("HC", "PD"))

  if (nrow(bl_expr) < 20) {
    log_md("- SKIPPED: < 20 BL HC+PD samples with neutrophil data")
    return(NULL)
  }

  log_md("- n BL samples (HC+PD): ", n_distinct(bl_expr$sample_id))

  sensitivity_rows <- list()

  for (iso in names(TRANSCRIPTS)) {
    iso_data <- bl_expr %>% filter(isoform == iso)
    if (nrow(iso_data) < 20) next

    # Model 1: Primary (without neutrophil adjustment)
    fit1 <- lm(log2cpm ~ DIAGNOSIS + GENDER + `RIN Value` + Plate, data = iso_data)
    robust1 <- tryCatch(
      coeftest(fit1, vcov. = vcovHC(fit1, type = "HC3")),
      error = function(e) NULL
    )

    # Model 2: Sensitivity (with neutrophil adjustment)
    fit2 <- lm(log2cpm ~ DIAGNOSIS + GENDER + `RIN Value` + Plate + neutrophil_fraction,
               data = iso_data)
    robust2 <- tryCatch(
      coeftest(fit2, vcov. = vcovHC(fit2, type = "HC3")),
      error = function(e) NULL
    )

    # Extract DIAGNOSISPD coefficient from both
    diag_coef_name <- "DIAGNOSISPD"
    if (diag_coef_name %in% rownames(coef(summary(fit1)))) {
      s1 <- coef(summary(fit1))[diag_coef_name, ]
      s2 <- coef(summary(fit2))[diag_coef_name, ]

      # HC3 robust p-values if available
      p_robust1 <- if (!is.null(robust1) && diag_coef_name %in% rownames(robust1))
        robust1[diag_coef_name, "Pr(>|t|)"] else NA
      p_robust2 <- if (!is.null(robust2) && diag_coef_name %in% rownames(robust2))
        robust2[diag_coef_name, "Pr(>|t|)"] else NA

      row1 <- tibble(
        isoform       = iso,
        model         = "Primary (no neutrophil adj)",
        logFC         = s1["Estimate"],
        SE            = s1["Std. Error"],
        P.Value       = s1["Pr(>|t|)"],
        P.Value_HC3   = p_robust1,
        source        = source_tag
      )
      row2 <- tibble(
        isoform       = iso,
        model         = paste0("Sensitivity (+", source_label, " adj)"),
        logFC         = s2["Estimate"],
        SE            = s2["Std. Error"],
        P.Value       = s2["Pr(>|t|)"],
        P.Value_HC3   = p_robust2,
        source        = source_tag
      )
      sensitivity_rows <- c(sensitivity_rows, list(row1, row2))

      # Log the comparison
      log_md("#### ", iso)
      log_md("- **Without** neutrophil adj: logFC = ", round(s1["Estimate"], 4),
             ", p = ", format.pval(s1["Pr(>|t|)"], digits = 3),
             if (!is.na(p_robust1)) paste0(", p_HC3 = ", format.pval(p_robust1, digits = 3)) else "")
      log_md("- **With** ", source_label, " adj: logFC = ", round(s2["Estimate"], 4),
             ", p = ", format.pval(s2["Pr(>|t|)"], digits = 3),
             if (!is.na(p_robust2)) paste0(", p_HC3 = ", format.pval(p_robust2, digits = 3)) else "")

      # Change in logFC
      pct_change <- round((s2["Estimate"] - s1["Estimate"]) / abs(s1["Estimate"]) * 100, 1)
      log_md("- logFC change: ", round(s2["Estimate"] - s1["Estimate"], 4),
             " (", pct_change, "%)")
      if (abs(pct_change) < 20) {
        log_md("  → **logFC stable** (< 20% change) → per-cell activation, not cell count artifact")
      } else {
        log_md("  → logFC shifted > 20% → cell abundance partially contributes")
      }
    }
  }

  if (length(sensitivity_rows) > 0) {
    sensitivity_df <- bind_rows(sensitivity_rows)
    write_csv(sensitivity_df,
              file.path(out_dir, "data", paste0("DA_sensitivity_", source_tag, ".csv")))
    log_md("")
    log_md("Interpretation: if logFC stays similar after neutrophil adjustment →")
    log_md("PADI4 upregulation is per-cell activation, not just cell count increase")
    log_md("")
    return(sensitivity_df)
  }
  NULL
}

# =============================================================================
# PART 3: Load expression data + run per-cell + limma analyses
# =============================================================================
message("Part 3: Loading transcript-level expression...")

load(file.path(DATA_DIR, "Copy of df_normed_filtered_annotated.RData"))
NAME.log2.cpm.filtered.norm.df <- NAME.log2.cpm.filtered.norm.df %>%
  dplyr::select(-last_col())

# Also load gene annotation for limma (used inside helper)
gene_annot <- NAME.log2.cpm.filtered.norm.df %>%
  dplyr::select(geneID, external_gene_name) %>%
  distinct() %>%
  filter(!is.na(external_gene_name), external_gene_name != "")

net_long <- NAME.log2.cpm.filtered.norm.df %>%
  filter(geneID %in% TRANSCRIPTS) %>%
  pivot_longer(cols = `3174_V08`:last_col(),
               names_to  = "sample_id",
               values_to = "log2cpm") %>%
  mutate(isoform = names(TRANSCRIPTS)[match(geneID, TRANSCRIPTS)])

# Determine which sources to run
sources_to_run <- switch(
  NEUT_SOURCE,
  "deconvolution" = list(
    list(col = "neutrophil_frac_deconv", label = "CIBERSORTx Deconvolution", tag = "deconv")
  ),
  "lab_blood" = list(
    list(col = "neu_blood_pct", label = "Blood Test Neutrophils (%)", tag = "labblood")
  ),
  "both" = {
    sources <- list(
      list(col = "neutrophil_frac_deconv", label = "CIBERSORTx Deconvolution", tag = "deconv")
    )
    if (!is.null(lab_neut)) {
      sources <- c(sources, list(
        list(col = "neu_blood_pct", label = "Blood Test Neutrophils (%)", tag = "labblood")
      ))
    }
    sources
  }
)

# ── PART 4: Per-cell PADI4 analysis ──────────────────────────────────────────
message("Part 4: Per-cell PADI4 activation analysis...")
log_md("## Part 4: Per-cell PADI4 analysis")

all_percell_results <- list()

for (src in sources_to_run) {
  # For lab_blood, convert % to fraction for consistency in per-cell analysis
  meta_for_percell <- ciber_meta
  if (src$tag == "labblood") {
    meta_for_percell[[src$col]] <- meta_for_percell[[src$col]]  # keep as % for scatter labeling
  }

  results <- run_percell_analysis(
    data_with_neut = meta_for_percell,
    neut_col_name  = src$col,
    source_label   = src$label,
    source_tag     = src$tag,
    net_long_data  = net_long,
    out_dir        = OUT_DIR
  )
  all_percell_results[[src$tag]] <- results
}

# ── PART 5: Isoform-level sensitivity ± neutrophil adjustment ─────────────────
message("Part 5: Isoform sensitivity analysis...")
log_md("## Part 5: Isoform sensitivity (± neutrophil adjustment)")

all_sensitivity_results <- list()

for (src in sources_to_run) {
  meta_for_sens <- ciber_meta
  # For sensitivity, need fraction (0-1 scale) for lab blood
  if (src$tag == "labblood") {
    meta_for_sens$neutrophil_fraction_sens <- meta_for_sens[[src$col]] / 100
    sens_col <- "neutrophil_fraction_sens"
  } else {
    sens_col <- src$col
  }

  sens_res <- run_isoform_sensitivity(
    meta_bl        = meta_for_sens,
    neut_col_name  = sens_col,
    source_label   = src$label,
    source_tag     = src$tag,
    net_long_data  = net_long,
    out_dir        = OUT_DIR
  )
  all_sensitivity_results[[src$tag]] <- sens_res
}

# =============================================================================
# PART 6: Summary tables
# =============================================================================
message("Part 6: Generating summary tables...")
log_md("## Part 6: Summary")

# Per-cell summary
summary_rows <- list()
for (src_tag in names(all_percell_results)) {
  for (iso in names(all_percell_results[[src_tag]])) {
    summary_rows <- c(summary_rows, list(tibble(
      isoform  = iso,
      source   = src_tag,
      wilcox_p = all_percell_results[[src_tag]][[iso]]$wilcox_p
    )))
  }
}
if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  write_csv(summary_df, file.path(OUT_DIR, "data", "DA_percell_summary.csv"))
  log_md("### Per-cell Wilcoxon p-values:")
  log_md("```"); log_md(capture.output(print(summary_df))); log_md("```")
}

# Sensitivity comparison (if both sources)
if (NEUT_SOURCE == "both" && length(all_sensitivity_results) == 2) {
  combined_sens <- bind_rows(all_sensitivity_results)
  write_csv(combined_sens, file.path(OUT_DIR, "data", "DA_sensitivity_combined.csv"))
  log_md("### Combined sensitivity comparison (both sources) saved")
}

# =============================================================================
# Output manifest
# =============================================================================
log_md("")
log_md("## Output files")
log_md("### Data:")
log_md("- data/DA_neutrophil_fractions.csv         [merged BL metadata + both neutrophil sources]")
log_md("- data/DA_CIBERSORTx_full_22celltypes.csv  [full 22 cell-type deconvolution + metadata]")
log_md("- data/DA_validation_deconv_vs_lab.csv      [deconvolution vs lab-blood scatter data]")
log_md("- data/DA_percell_summary.csv               [Wilcoxon p-values per isoform per source]")
log_md("- data/DA_residuals_*_deconv.csv            [per-sample residuals, deconvolution source]")
log_md("- data/DA_residuals_*_labblood.csv          [per-sample residuals, lab-blood source]")
log_md("- data/DA_sensitivity_deconv.csv            [logFC with/without deconvolution covariate]")
log_md("- data/DA_sensitivity_labblood.csv          [logFC with/without lab-blood covariate]")
log_md("- data/DA_sensitivity_combined.csv          [combined comparison, if NEUT_SOURCE='both']")
log_md("")
log_md("### Figures:")
log_md("- pdf/DA_validation_deconv_vs_lab.pdf       [validation scatter: deconv vs lab neutrophils]")
log_md("- pdf/DA_neutrophil_violin_deconv.pdf       [neutrophil fraction violin by diagnosis]")
log_md("- pdf/DA_neutrophil_violin_labblood.pdf     [lab blood neutrophil violin by diagnosis]")
log_md("- pdf/DA_percell_*_deconv.pdf               [scatter + residual violin, deconv source]")
log_md("- pdf/DA_percell_*_labblood.pdf             [scatter + residual violin, lab-blood source]")

message("Done — Baseline_Neutrophil_Deconvolution.R")
