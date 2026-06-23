# =============================================================================
# Longitudinal_Trajectories.R
# Longitudinal whole-blood trajectories of NET transcripts (Figure 2)
# Random-slope linear mixed models + Bayesian (brms) sensitivity analysis.
# =============================================================================
# Inputs (place in ./data or pass as the first command-line argument):
#   df_normed_filtered_annotated.RData  -> NAME.log2.cpm.filtered.norm.df
#   metaDataIR3.csv
# Transcripts: PADI4-201, PADI4-202, MPO-201
# Usage:  Rscript Longitudinal_Trajectories.R [DATA_DIR] [OUT_DIR]
# Output: <OUT_DIR>/{data,pdf,png}/B_longitudinal_*
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(sandwich)
  library(lmtest)
  library(brms)
  library(bayestestR)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
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

LOG_FILE <- file.path(OUT_DIR, "B_longitudinal_log.md")
log_md   <- function(...) cat(paste0(..., "\n"), file = LOG_FILE, append = TRUE)
cat("", file = LOG_FILE)
log_md("# Longitudinal Analysis Log — Longitudinal_Trajectories.R")
log_md("Date: ", format(Sys.time()))
log_md("")

# ── Constants ──────────────────────────────────────────────────────────────────
TRANSCRIPTS <- c(
  "PADI4-201" = "ENST00000375448.4",
  "PADI4-202" = "ENST00000375453.5",
  "MPO-201"   = "ENST00000225275.3"
)

custom_colors <- c("HC" = "#648FFF", "Prodromal" = "#785EF0", "PD" = "#DC267F")

BRMS_CHAINS <- 4
BRMS_ITER   <- 4000
BRMS_WARMUP <- 2000
BRMS_SEED   <- 42
RUN_BRMS    <- TRUE

# ── Adaptive fitter ────────────────────────────────────────────────────────────
fit_lmm <- function(formula, data, label) {
  m <- tryCatch(
    lmer(formula, data = data, REML = TRUE,
         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
    error = function(e) NULL
  )

  # Only accept the model if it successfully converged AND is NOT singular
  if (!is.null(m) && !lme4::isSingular(m)) {
    return(list(model = m, type = "lmer_bobyqa", label = label))
  }

  # Fallback: random intercept only if the primary model failed OR was singular
  fi <- as.character(formula)
  fi_simple <- reformulate(
    sub("\\(timepoint_num\\|PATNO\\)", "(1|PATNO)", fi[3]),
    response = fi[2]
  )
  m2 <- tryCatch(
    lmer(fi_simple, data = data, REML = TRUE),
    error = function(e) NULL
  )
  if (!is.null(m2)) {
    type_str <- if (is.null(m)) "lmer_intercept_only [fallback error]" else "lmer_intercept_only [fallback singular]"
    return(list(model = m2, type = type_str, label = paste(label, "[fallback RI]")))
  }
  return(NULL)
}

# ── Load data ──────────────────────────────────────────────────────────────────
message("Loading expression data...")
load(file.path(DATA_DIR, "Copy of df_normed_filtered_annotated.RData"))
NAME.log2.cpm.filtered.norm.df <- NAME.log2.cpm.filtered.norm.df %>%
  dplyr::select(-last_col())

samples <- read_csv(file.path(DATA_DIR, "metaDataIR3.csv"), show_col_types = FALSE) %>%
  filter(DIAGNOSIS %in% c("PD", "Prodromal", "Control")) %>%
  mutate(sample_id = paste0(PATNO, "_", CLINICAL_EVENT))

net_long <- NAME.log2.cpm.filtered.norm.df %>%
  filter(geneID %in% TRANSCRIPTS) %>%
  pivot_longer(cols = `3174_V08`:last_col(),
               names_to  = "sample_id",
               values_to = "log2cpm")

long_data <- samples %>%
  dplyr::select(sample_id, PATNO, CLINICAL_EVENT, DIAGNOSIS, GENDER,
                `RIN Value`, Plate) %>%
  inner_join(net_long, by = "sample_id") %>%
  mutate(
    DIAGNOSIS = recode(DIAGNOSIS, Control = "HC"),
    DIAGNOSIS = factor(DIAGNOSIS, levels = c("HC", "Prodromal", "PD")),
    isoform   = names(TRANSCRIPTS)[match(geneID, TRANSCRIPTS)],
    # timepoint_num (numeric 0,6,12,24,36): used ONLY in the random-slope term (timepoint_num|PATNO)
    # so lme4 can estimate a continuous rate-of-change per subject across time.
    # MONTH_FOLLOW (ordered factor): used ONLY in fixed effects so emmeans can generate
    # labeled pairwise contrasts per timepoint. Using a factor in fixed effects + a numeric
    # in the random slope is the standard mixed-model approach when you need both categorical
    # contrasts AND subject-level slope heterogeneity. The two variables carry identical
    # information encoded differently — this is intentional design, not duplication.
    timepoint_num = case_when(
      CLINICAL_EVENT == "BL"  ~ 0.0,
      CLINICAL_EVENT == "V02" ~ 0.5,
      CLINICAL_EVENT == "V04" ~ 1.0,
      CLINICAL_EVENT == "V06" ~ 2.0,
      TRUE                    ~ 3.0
    ),
    MONTH_FOLLOW = factor(
      case_when(
        CLINICAL_EVENT == "BL"  ~ "0 month",
        CLINICAL_EVENT == "V02" ~ "6 months",
        CLINICAL_EVENT == "V04" ~ "12 months",
        CLINICAL_EVENT == "V06" ~ "24 months",
        TRUE                    ~ "36 months"
      ),
      levels = c("0 month", "6 months", "12 months", "24 months", "36 months")
    ),
    sex   = factor(GENDER),
    plate = factor(Plate),
    rin   = as.numeric(`RIN Value`)
  )

log_md("## Data summary")
log_md("- n_subjects: ", n_distinct(long_data$PATNO))
log_md("- n_samples:  ", nrow(long_data))
log_md("- **Variable design:** MONTH_FOLLOW (factor) used in fixed effects for emmeans categorical")
log_md("  contrasts; timepoint_num (numeric 0/6/12/24/36) used in (timepoint_num|PATNO) random")
log_md("  slope for continuous per-subject rate-of-change. Both encode the same timepoints in")
log_md("  different representations required by their respective model components.")
log_md("")

# =============================================================================
# FREQUENTIST: Random slopes LMM (HC vs PD, PD vs Prodromal)
# =============================================================================

all_freq_results <- list()
all_freq_emm     <- list()

for (iso in names(TRANSCRIPTS)) {
  message("Frequentist LMM: ", iso)
  log_md("## ", iso)

  iso_data <- long_data %>% filter(isoform == iso)

  for (comparison in list(c("HC", "PD"), c("HC", "Prodromal", "PD"))) {
    comp_label <- paste(comparison, collapse = "_vs_")
    cdata      <- iso_data %>%
      filter(DIAGNOSIS %in% comparison) %>%
      mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = comparison))

    # Primary model: random slopes
    # MONTH_FOLLOW (factor) → fixed categorical effects + emmeans contrasts
    # (timepoint_num|PATNO) → continuous random slope per subject
    fit <- fit_lmm(
      log2cpm ~ DIAGNOSIS * MONTH_FOLLOW + sex + rin + (timepoint_num | PATNO),
      data  = cdata,
      label = paste0(iso, " [", comp_label, "]")
    )

    if (is.null(fit)) {
      log_md("### ", comp_label, ": model failed")
      next
    }

    log_md("### ", comp_label, " — Model type: ", fit$type)
    log_md("- n = ", nrow(cdata), " | n_subjects = ", n_distinct(cdata$PATNO))

    anova_res <- anova(fit$model)
    log_md("```")
    log_md(capture.output(print(anova_res)))
    log_md("```")

    # emmeans pairwise with Holm
    emm <- emmeans(fit$model,
                   specs = pairwise ~ DIAGNOSIS | MONTH_FOLLOW,
                   adjust = "holm")
    emm_df <- as.data.frame(emm$contrasts)

    log_md("#### emmeans contrasts (Holm-adjusted):")
    log_md("```")
    log_md(capture.output(print(emm_df)))
    log_md("```")

    # HC3 robust SE
    robust_coef <- tryCatch(
      coeftest(fit$model, vcov. = vcovHC(fit$model, type = "HC3")),
      error = function(e) NULL
    )
    if (!is.null(robust_coef)) {
      log_md("#### HC3 Robust SE:")
      log_md("```")
      log_md(capture.output(print(robust_coef)))
      log_md("```")
    }

    all_freq_results[[paste0(iso, "_", comp_label)]] <- list(
      anova    = anova_res,
      emmeans  = emm_df,
      model    = fit$model,
      type     = fit$type,
      n        = nrow(cdata),
      n_subj   = n_distinct(cdata$PATNO)
    )
    all_freq_emm[[paste0(iso, "_", comp_label)]] <- emm
  }
  log_md("")
}

# =============================================================================
# PREDICTED TRAJECTORIES PLOT (emmeans means ± 95% CI)
# =============================================================================

plot_predicted_trajectories <- function(iso, comparison = c("HC", "PD")) {
  cdata <- long_data %>%
    filter(isoform == iso, DIAGNOSIS %in% comparison) %>%
    mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = comparison))

  # MONTH_FOLLOW (factor) → categorical fixed effects; (timepoint_num|PATNO) → continuous random slope
  fit <- fit_lmm(
    log2cpm ~ DIAGNOSIS * MONTH_FOLLOW + sex + rin + (timepoint_num | PATNO),
    data = cdata, label = iso
  )
  if (is.null(fit)) return(NULL)

  pred_means <- emmeans(fit$model, ~ MONTH_FOLLOW * DIAGNOSIS) %>%
    as.data.frame()

  n_table <- cdata %>%
    group_by(MONTH_FOLLOW) %>%
    summarise(n = n(), .groups = "drop")

  x_labels <- setNames(
    paste0(n_table$MONTH_FOLLOW, "\n(n=", n_table$n, ")"),
    as.character(n_table$MONTH_FOLLOW)
  )

  anova_res <- anova(fit$model)
  row_names <- rownames(anova_res)

  t_idx <- grep("MONTH_FOLLOW$", row_names)
  d_idx <- grep("DIAGNOSIS$", row_names)
  i_idx <- grep(":", row_names)

  t_res <- if (length(t_idx) > 0) anova_res[t_idx[1], ] else NULL
  d_res <- if (length(d_idx) > 0) anova_res[d_idx[1], ] else NULL
  i_res <- if (length(i_idx) > 0) anova_res[i_idx[1], ] else NULL

  fmt_p <- function(p) {
    if (is.na(p) || is.nan(p)) {
      return("NA")
    }
    if (p < 0.0001) sprintf("%.2e", p) else sprintf("%.4f", p)
  }

  subtitle <- sprintf(
    "Time: F(%.1f,%.1f)=%.2f, p=%s | Diagnosis: F(%.1f,%.1f)=%.2f, p=%s\nInteraction: F(%.1f,%.1f)=%.2f, p=%s | Model: %s",
    t_res$NumDF, t_res$DenDF, t_res$`F value`, fmt_p(t_res$`Pr(>F)`),
    d_res$NumDF, d_res$DenDF, d_res$`F value`, fmt_p(d_res$`Pr(>F)`),
    i_res$NumDF, i_res$DenDF, i_res$`F value`, fmt_p(i_res$`Pr(>F)`),
    fit$type
  )

  ggplot(pred_means, aes(x = MONTH_FOLLOW, y = emmean,
                          color = DIAGNOSIS, group = DIAGNOSIS)) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = DIAGNOSIS),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.8) +
    scale_color_manual(values = custom_colors[comparison]) +
    scale_fill_manual(values  = custom_colors[comparison]) +
    scale_x_discrete(labels  = x_labels) +
    labs(
      title    = iso,
      subtitle = subtitle,
      y        = bquote(bold(.(iso) ~ "mRNA (log"["2"] * "CPM)")),
      x        = "",
      color    = "Diagnosis", fill = "Diagnosis"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 8, face = "italic", lineheight = 1.1),
      axis.title    = element_text(face = "bold"),
      axis.text.x   = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# Generate and save trajectory plots
traj_plots <- list()
for (iso in names(TRANSCRIPTS)) {
  p_hc_pd <- plot_predicted_trajectories(iso, c("HC", "PD"))
  p_all   <- plot_predicted_trajectories(iso, c("HC", "Prodromal", "PD"))

  if (!is.null(p_hc_pd) && !is.null(p_all)) {
    combined <- p_hc_pd + p_all + plot_annotation(tag_levels = "A")
    ggsave(file.path(OUT_DIR, "pdf", paste0("B_traj_", iso, ".pdf")),
           combined, width = 14, height = 6)
    ggsave(file.path(OUT_DIR, "png", paste0("B_traj_", iso, ".png")),
           combined, width = 14, height = 6, dpi = 300)
    traj_plots[[iso]] <- list(hc_pd = p_hc_pd, all = p_all)
  }
}

# Individual spaghetti + boxplot per diagnosis (HC and PD separately)
plot_individual_trajectories <- function(iso, diag_group) {
  cdata <- long_data %>%
    filter(isoform == iso, DIAGNOSIS %in% c("HC", "PD")) %>%
    mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("HC", "PD"))) %>%
    filter(DIAGNOSIS == diag_group)

  # MONTH_FOLLOW (factor) → categorical; (timepoint_num|PATNO) → continuous slope
  fit <- fit_lmm(
    log2cpm ~ MONTH_FOLLOW + sex + rin + (timepoint_num | PATNO),
    data = cdata, label = paste0(iso, "_", diag_group)
  )

  posthoc <- if (!is.null(fit)) {
    tryCatch(
      emmeans(fit$model, specs = trt.vs.ctrl ~ MONTH_FOLLOW,
              ref = 1, adjust = "dunnett"),
      error = function(e) NULL
    )
  } else NULL

  max_y    <- max(cdata$log2cpm, na.rm = TRUE)
  n_table  <- cdata %>% group_by(MONTH_FOLLOW) %>% summarise(n = n(), .groups = "drop")
  x_labels <- paste0(n_table$MONTH_FOLLOW, "\n(n=", n_table$n, ")")

  stats_df <- data.frame()
  if (!is.null(posthoc)) {
    stats_df <- as.data.frame(posthoc$contrasts) %>%
      mutate(
        group1 = "0 month",
        group2 = gsub(" - 0 month", "", contrast),
        p_star = case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01  ~ "**",
          p.value < 0.05  ~ "*",
          TRUE            ~ "ns"
        )
      ) %>%
      filter(p_star != "ns") %>%
      mutate(y.position = seq(max_y + 0.05 * max_y, by = 0.35, length.out = n()))
  }

  diag_color <- custom_colors[diag_group]

  p <- ggplot(cdata, aes(x = MONTH_FOLLOW, y = log2cpm)) +
    geom_line(aes(group = PATNO), color = "grey80", alpha = 0.4) +
    geom_boxplot(fill = diag_color, outlier.shape = NA, alpha = 0.45,
                 show.legend = FALSE, width = 0.5) +
    geom_jitter(color = diag_color, width = 0.12, alpha = 0.35,
                size = 1.8, show.legend = FALSE) +
    scale_x_discrete(labels = x_labels) +
    labs(
      title = paste0(iso, " — ", diag_group),
      y     = bquote(bold(.(iso) ~ "mRNA (log"["2"] * "CPM)")),
      x     = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 13),
      axis.text.x      = element_text(size = 9),
      panel.grid.minor = element_blank()
    )

  if (nrow(stats_df) > 0) {
    p <- p + stat_pvalue_manual(
      stats_df, label = "p_star", tip.length = 0.01,
      size = 6, bracket.size = 0.6, inherit.aes = FALSE
    )
  }
  p
}

for (iso in names(TRANSCRIPTS)) {
  p_hc <- plot_individual_trajectories(iso, "HC")
  p_pd <- plot_individual_trajectories(iso, "PD")
  combined <- p_hc + p_pd + plot_annotation(tag_levels = "A")
  ggsave(file.path(OUT_DIR, "pdf", paste0("B_individual_", iso, ".pdf")),
         combined, width = 14, height = 6)
  ggsave(file.path(OUT_DIR, "png", paste0("B_individual_", iso, ".png")),
         combined, width = 14, height = 6, dpi = 300)
}

# =============================================================================
# BAYESIAN: brms with random slopes + LKJ(2) prior on correlation
# =============================================================================
if (RUN_BRMS) {
  message("Bayesian random slopes with brms...")
  log_md("## Bayesian Sensitivity (brms)")

  brms_results <- list()

  for (iso in names(TRANSCRIPTS)) {
    cdata <- long_data %>%
      filter(isoform == iso, DIAGNOSIS %in% c("HC", "PD")) %>%
      mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("HC", "PD")))

    message("  brms: ", iso, " (n=", nrow(cdata), ")")

    brms_fit <- tryCatch(
      brm(
        # MONTH_FOLLOW (factor) → categorical fixed effects; (timepoint_num|PATNO) → continuous random slope
        log2cpm ~ DIAGNOSIS * MONTH_FOLLOW + sex + rin + (timepoint_num | PATNO),
        data   = cdata,
        family = gaussian(),
        prior  = c(
          prior(normal(0, 1),        class = "b"),
          prior(normal(0, 5),        class = "Intercept"),
          prior(student_t(3, 0, 1),  class = "sd"),
          prior(lkj(2),              class = "cor"),
          prior(student_t(3, 0, 1),  class = "sigma")
        ),
        chains  = BRMS_CHAINS,
        iter    = BRMS_ITER,
        warmup  = BRMS_WARMUP,
        seed    = BRMS_SEED,
        cores   = 4,
        refresh = 0,
        silent  = 2,
        control = list(adapt_delta = 0.99, max_treedepth = 15)
      ),
      error = function(e) { message("  brms error: ", conditionMessage(e)); NULL }
    )

    if (!is.null(brms_fit)) {
      summ <- summary(brms_fit)
      log_md("### ", iso)
      log_md("```")
      log_md(capture.output(print(summ$fixed)))
      log_md("```")

      # Posterior probability that PD effect > 0 at each timepoint
      post_draws <- as_draws_df(brms_fit)
      diag_coefs <- grep("^b_DIAGNOSISPD", names(post_draws), value = TRUE)
      pp_df <- tibble(
        parameter   = diag_coefs,
        post_prob_gt0 = sapply(diag_coefs, function(x) mean(post_draws[[x]] > 0))
      )
      log_md("#### Posterior P(PD effect > 0):")
      log_md("```")
      log_md(capture.output(print(pp_df)))
      log_md("```")

      # ROPE analysis
      rope_res <- tryCatch(
        bayestestR::rope(brms_fit, range = c(-0.1, 0.1), ci = 0.95),
        error = function(e) NULL
      )
      if (!is.null(rope_res)) {
        log_md("#### ROPE (|effect| < 0.1 as negligible threshold):")
        log_md("```")
        log_md(capture.output(print(rope_res)))
        log_md("```")
      }

      saveRDS(brms_fit, file.path(OUT_DIR, "data", paste0("B_brms_", iso, ".rds")))
      brms_results[[iso]] <- list(fit = brms_fit, pp = pp_df)
    }
    log_md("")
  }
}

# =============================================================================
# Save frequency model table
# =============================================================================
model_summary <- bind_rows(lapply(names(all_freq_results), function(k) {
  res <- all_freq_results[[k]]
  tibble(
    analysis   = k,
    Model_Type = res$type,
    n          = res$n,
    n_subjects = res$n_subj
  )
}))
write_csv(model_summary, file.path(OUT_DIR, "data", "B_model_summary.csv"))

emmeans_all <- bind_rows(lapply(names(all_freq_results), function(k) {
  all_freq_results[[k]]$emmeans %>% mutate(analysis = k)
}))
write_csv(emmeans_all, file.path(OUT_DIR, "data", "B_emmeans_contrasts.csv"))

log_md("")
log_md("## Output files")
log_md("- data/B_model_summary.csv")
log_md("- data/B_emmeans_contrasts.csv")
log_md("- data/B_brms_*.rds (if brms ran)")
log_md("- pdf/B_traj_*.pdf, B_individual_*.pdf")

message("Done — Longitudinal_Trajectories.R")
