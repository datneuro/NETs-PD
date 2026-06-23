# =============================================================================
# Serum_Biomarker_ROC.R
# Diagnostic ROC for serum MPO-DNA and CitH3-DNA blood biomarkers (Figure 3)
# Merged ROC curves, bootstrap 95% CI bands, Youden cut-off / sensitivity / specificity.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(pROC)
  library(ggplot2)
})

# ── Directory Configurations ──────────────────────────────────────────────────
MPO_FILE   <- file.path("data", "mpo_final.xlsx")
CITH3_FILE <- file.path("data", "CitH3_final.xlsx")
OUT_DIR    <- "results"

if (!file.exists(MPO_FILE) || !file.exists(CITH3_FILE)) {
  stop("Critical Error: One or both source Excel data matrices are missing. Verify paths.")
}

# ── 1. Load and Compute MPO-DNA ROC Matrix ─────────────────────────────────────
message("Extracting MPO-DNA Matrix & Executing Bootstrap Resampling...")
df_mpo <- readxl::read_excel(MPO_FILE) %>%
  filter(!is.na(lg_MPO_DNA), !is.na(diag)) %>%
  mutate(response = ifelse(diag == "PD", 1, 0))

roc_mpo    <- pROC::roc(df_mpo$response, df_mpo$lg_MPO_DNA, direction = ">", quiet = TRUE)
ci_auc_mpo <- pROC::ci.auc(roc_mpo)

# Restored: Extract the mathematically optimal Youden coordinates
coords_mpo <- pROC::coords(roc_mpo, "best", ret = c("threshold", "specificity", "sensitivity"))

set.seed(42) # Lock seed for reproducible bootstrap tracks
ci_se_mpo  <- pROC::ci.se(roc_mpo, specificities = seq(0, 1, 0.01), method = "bootstrap", boot.n = 500, boot.stratified = TRUE)

dat_mpo <- data.frame(
  Specificity = as.numeric(rownames(ci_se_mpo)),
  Lower       = ci_se_mpo[, 1],
  Sensitivity = ci_se_mpo[, 2],
  Upper       = ci_se_mpo[, 3],
  Biomarker   = "MPO-DNA"
)

# ── 2. Load and Compute CitH3-DNA ROC Data ───────────────────────────────────
message("Extracting CitH3-DNA Matrix & Executing Bootstrap Resampling...")
df_cith3 <- readxl::read_excel(CITH3_FILE) %>%
  filter(!is.na(lg10CitH3), !is.na(diag)) %>%
  mutate(response = ifelse(diag == "PD", 1, 0))

roc_cith3    <- pROC::roc(df_cith3$response, df_cith3$lg10CitH3, direction = ">", quiet = TRUE)
ci_auc_cith3 <- pROC::ci.auc(roc_cith3)

# Restored: Extract the mathematically optimal Youden coordinates
coords_cith3 <- pROC::coords(roc_cith3, "best", ret = c("threshold", "specificity", "sensitivity"))

set.seed(42)
ci_se_cith3  <- pROC::ci.se(roc_cith3, specificities = seq(0, 1, 0.01), method = "bootstrap", boot.n = 500, boot.stratified = TRUE)

dat_cith3 <- data.frame(
  Specificity = as.numeric(rownames(ci_se_cith3)),
  Lower       = ci_se_cith3[, 1],
  Sensitivity = ci_se_cith3[, 2],
  Upper       = ci_se_cith3[, 3],
  Biomarker   = "CitH3-DNA"
)

# ── 3. Merge Datasets for Consolidated Plotting ──────────────────────────────
df_merged_roc <- bind_rows(dat_mpo, dat_cith3) %>%
  mutate(Biomarker = factor(Biomarker, levels = c("MPO-DNA", "CitH3-DNA")))

# ── 4. Setup Precision Nomenclature Subtitles (FULL INTEGRATION) ─────────────
# Re-assembled with perfect synchronization of AUC, CI, Cut-offs, Sen, and Spec metrics
stat_subtitle <- paste0(
  sprintf("MPO-DNA AUC: %.3f (95%% CI: [%.3f, %.3f]) | Optimal Cut-off: %.3f | Sen: %.1f%% | Spec: %.1f%%\n", 
          roc_mpo$auc, ci_auc_mpo[1], ci_auc_mpo[3], coords_mpo$threshold, coords_mpo$sensitivity*100, coords_mpo$specificity*100),
  sprintf("CitH3-DNA AUC: %.3f (95%% CI: [%.3f, %.3f]) | Optimal Cut-off: %.3f | Sen: %.1f%% | Spec: %.1f%%\n", 
          roc_cith3$auc, ci_auc_cith3[1], ci_auc_cith3[3], coords_cith3$threshold, coords_cith3$sensitivity*100, coords_cith3$specificity*100),
  "95% Confidence Interval bands computed via Stratified Bootstrap Resampling (R = 500 iterations)"
)

# Corporate Color Palette Mapping
biomarker_colors <- c("MPO-DNA" = "#648FFF", "CitH3-DNA" = "#DC267F")

# ─────────────────────────────────────────────────────────────────────────────
# 5. RENDER CONSOLIDATED PUBLICATION-READY ROC PLOT
# ─────────────────────────────────────────────────────────────────────────────
message("Generating official unified ROC asset...")

p_roc <- ggplot(df_merged_roc, aes(x = 1 - Specificity, y = Sensitivity, color = Biomarker)) +
  # Layer 1: Plot Reference Diagonal Line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  # Layer 2: Render Bootstrap 95% Confidence Interval Bands (Ribbons)
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Biomarker), color = NA, alpha = 0.12) +
  # Layer 3: Render Exact Solid ROC Curves
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = biomarker_colors) +
  scale_fill_manual(values = biomarker_colors) +
  xlim(0, 1) + ylim(0, 1) +
  labs(
    title = "a", 
    subtitle = stat_subtitle, 
    x = "1 - Specificity (False Positive Rate)", 
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title          = element_text(face = "bold", size = 13, hjust = 0), 
    plot.title.position = "plot", 
    plot.subtitle       = element_text(size = 7.5, face = "bold", color = "grey30", lineheight = 1.25),
    axis.text.x         = element_text(color = "black"),
    axis.text.y         = element_text(color = "black"),
    axis.title.x        = element_text(face = "bold", size = 9.5),
    axis.title.y        = element_text(face = "bold", size = 9.5),
    legend.title        = element_blank(),
    legend.text         = element_text(face = "bold", size = 9),
    legend.position     = c(0.80, 0.15), # Locked neatly in the empty bottom-right quadrant
    legend.background   = element_rect(fill = "white", color = "grey90", linewidth = 0.3),
    panel.grid.minor    = element_blank(),
    plot.margin         = margin(t = 10, r = 15, b = 10, l = 15)
  )

# ── 6. Export Asset (Micro-expanded width to comfortably hold clinical text lines) ──
ggsave(file.path(OUT_DIR, "pdf", "Supplementary_Figure_ROC_Biomarkers.pdf"), p_roc, width = 6.0, height = 4.8, device = "pdf")
ggsave(file.path(OUT_DIR, "png", "Supplementary_Figure_ROC_Biomarkers.png"), p_roc, width = 6.0, height = 4.8, dpi = 300)

message("Success: ROC figure saved to ", file.path(OUT_DIR, "pdf"))