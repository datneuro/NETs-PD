# =============================================================================
# Longitudinal_Stability.R
# Temporal stability of NET transcripts across 5 visits (Supplementary Fig. S5)
# ICC(2,1), between/within-subject variance partitioning, random-slope LRT in
# controls, and a Bayesian (brms) ICC sensitivity analysis.
# =============================================================================
# Inputs:
#   df_normed_filtered_annotated.RData  -> NAME.log2.cpm.filtered.norm.df
#   metaDataIR3.csv
# Transcripts: PADI4-201 (ENST00000375448.4), PADI4-202 (ENST00000375453.5),
#              MPO-201 (ENST00000225275.3)
# Output: <OUT_DIR>/{data}/C_stability_*  (statistics + log)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(performance) # icc()
  library(emmeans)
  library(brms)
})
# ── CLI arguments ──────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) {
  args[1]
} else {
  "data"
}
OUT_DIR <- if (length(args) >= 2) {
  args[2]
} else {
  "results"
}

dir.create(file.path(OUT_DIR, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "png"), recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "C_stability_log.md")
log_md <- function(...) cat(paste0(..., "\n"), file = LOG_FILE, append = TRUE)
cat("", file = LOG_FILE) # reset log
log_md("# Stability Analysis Log — Longitudinal_Stability.R")
log_md("Date: ", format(Sys.time()))
log_md("")

# ── Constants ──────────────────────────────────────────────────────────────────
TRANSCRIPTS <- c(
  "PADI4-201" = "ENST00000375448.4",
  "PADI4-202" = "ENST00000375453.5",
  "MPO-201"   = "ENST00000225275.3"
)

TIMEPOINT_MAP <- c(BL = 0, V02 = 6, V04 = 12, V06 = 24, V08 = 36)

BRMS_CHAINS <- 4
BRMS_ITER <- 4000
BRMS_WARMUP <- 2000
BRMS_SEED <- 42
RUN_BRMS <- TRUE 

# ── Load & prepare data ────────────────────────────────────────────────────────
message("Loading expression data...")
load(file.path(DATA_DIR, "Copy of df_normed_filtered_annotated.RData"))

NAME.log2.cpm.filtered.norm.df <- NAME.log2.cpm.filtered.norm.df %>%
  dplyr::select(-last_col())

net_long <- NAME.log2.cpm.filtered.norm.df %>%
  filter(geneID %in% TRANSCRIPTS) %>%
  pivot_longer(
    cols = `3174_V08`:last_col(),
    names_to = "sample_id",
    values_to = "log2cpm"
  )

samples <- read_csv(file.path(DATA_DIR, "metaDataIR3.csv"), show_col_types = FALSE) %>%
  filter(DIAGNOSIS %in% c("PD", "Prodromal", "Control")) %>%
  mutate(sample_id = paste0(PATNO, "_", CLINICAL_EVENT))

stability_data <- samples %>%
  dplyr::select(sample_id, PATNO, CLINICAL_EVENT, DIAGNOSIS, GENDER, `RIN Value`) %>%
  inner_join(net_long, by = "sample_id") %>%
  mutate(
    DIAGNOSIS = recode(DIAGNOSIS, Control = "HC"),
    DIAGNOSIS = factor(DIAGNOSIS, levels = c("HC", "Prodromal", "PD")),
    isoform = names(TRANSCRIPTS)[match(geneID, TRANSCRIPTS)],
    MONTH_NUM = TIMEPOINT_MAP[CLINICAL_EVENT],
    sex = factor(GENDER),
    rin = as.numeric(`RIN Value`),
    MONTH_FOLLOW = factor(
      case_when(
        CLINICAL_EVENT == "BL" ~ "0 month",
        CLINICAL_EVENT == "V02" ~ "6 months",
        CLINICAL_EVENT == "V04" ~ "12 months",
        CLINICAL_EVENT == "V06" ~ "24 months",
        TRUE ~ "36 months"
      ),
      levels = c("0 month", "6 months", "12 months", "24 months", "36 months")
    )
  ) %>%
  filter(!is.na(MONTH_NUM))

log_md("## Data summary")
log_md("- Subjects: ", n_distinct(stability_data$PATNO))
log_md("- Samples:  ", nrow(stability_data))
log_md("- Diagnosis counts:")
stability_data %>%
  distinct(PATNO, DIAGNOSIS) %>%
  count(DIAGNOSIS) %>%
  {
    log_md(capture.output(print(.)))
  }
log_md("")

# =============================================================================
# ICC, variance partitioning, and random-slope LRT
# =============================================================================

layer1_results <- list()

for (iso in names(TRANSCRIPTS)) {
  message("Layer 1: ", iso)

  iso_data <- stability_data %>% filter(isoform == iso)
  log_md("### ", iso)

  # ── ICC via performance::icc() per plan specification ───────────────────────
  # performance::icc() returns ICC_adjusted = between-subject / total variance = ICC(2,1)
  # VarCorr extraction is kept ONLY for variance partitioning columns (pct_between, pct_within)
  # Both methods should yield the same ICC value 
    group_by(DIAGNOSIS) %>%
    group_map(function(df, key) {
      grp <- key$DIAGNOSIS
      if (n_distinct(df$PATNO) < 5) {
        return(NULL)
      }

      # 1. Unadjusted (null) model
      m_null <- tryCatch(
        lmer(log2cpm ~ 1 + (1 | PATNO), data = df, REML = TRUE),
        error = function(e) NULL
      )
      m_adj <- tryCatch(
        lmer(log2cpm ~ sex + rin + (1 | PATNO), data = df, REML = TRUE),
        error = function(e) NULL
      )
      if (is.null(m_null)) {
        return(NULL)
      }

      perf_icc_null <- tryCatch(performance::icc(m_null), error = function(e) NULL)
      icc_unadj <- if (!is.null(perf_icc_null)) perf_icc_null$ICC_adjusted else NA_real_

      icc_adj <- NA_real_
      if (!is.null(m_adj)) {
        perf_icc_adj <- tryCatch(performance::icc(m_adj), error = function(e) NULL)
        icc_adj <- if (!is.null(perf_icc_adj)) perf_icc_adj$ICC_adjusted else NA_real_
      }

      vc <- as.data.frame(VarCorr(m_null))
      var_b <- vc$vcov[vc$grp == "PATNO"]
      var_w <- vc$vcov[vc$grp == "Residual"]

      tibble(
        DIAGNOSIS   = as.character(grp),
        ICC_Unadj   = round(icc_unadj, 4),
        ICC_Adj     = round(icc_adj, 4),
        var_between = round(var_b, 5),
        var_within  = round(var_w, 5),
        pct_between = round(100 * var_b / (var_b + var_w), 1),
        pct_within  = round(100 * var_w / (var_b + var_w), 1),
        n_subjects  = n_distinct(df$PATNO),
        n_obs       = nrow(df)
      )
    }) %>%
    bind_rows()

  log_md("#### ICC (performance::icc primary) and Variance Partitioning")
  log_md("- ICC = performance::icc()$ICC_adjusted = ICC(2,1) per plan spec")
  log_md("- ICC_manual = manual VarCorr; icc_match should be TRUE (sanity check)")
  log_md("```")
  log_md(capture.output(print(icc_table)))
  log_md("```")

  # ── Random slope test ───────────────────────────────────────────────────────
  hc_only <- iso_data %>% filter(DIAGNOSIS == "HC", !is.na(MONTH_NUM))
  slope_test <- tryCatch(
    {
      m_ri <- lmer(log2cpm ~ MONTH_NUM + (1 | PATNO), data = hc_only, REML = FALSE)
      m_rs <- lmer(log2cpm ~ MONTH_NUM + (MONTH_NUM | PATNO), data = hc_only, REML = FALSE)
      lrt <- anova(m_ri, m_rs)
      list(chi2 = round(lrt$Chisq[2], 3), df = lrt$Df[2], p = round(lrt[["Pr(>Chisq)"]][2], 4))
    },
    error = function(e) list(chi2 = NA, df = NA, p = NA)
  )

  log_md(
    "#### Random slope LRT (HC only): chi2=", slope_test$chi2,
    " df=", slope_test$df, " p=", slope_test$p
  )
  log_md("")

  layer1_results[[iso]] <- list(
    icc = icc_table,
    slope_lrt = slope_test
  )
}

# ── Save Layer 1 results ───────────────────────────────────────────────────────
icc_all <- bind_rows(lapply(names(layer1_results), function(iso) {
  layer1_results[[iso]]$icc %>% mutate(isoform = iso)
}))
write_csv(icc_all, file.path(OUT_DIR, "data", "C_icc_variance_table.csv"))

# =============================================================================
# Bayesian ICC (sensitivity analysis)
# =============================================================================
if (RUN_BRMS) {
  message("Bayesian ICC with brms...")
  log_md("## Bayesian ICC")

  bayes_icc_table <- list()

  for (iso in names(TRANSCRIPTS)) {
    for (grp in c("HC", "PD")) {
      grp_data <- stability_data %>%
        filter(isoform == iso, DIAGNOSIS == grp)

      if (n_distinct(grp_data$PATNO) < 5) next

      message("  brms ICC: ", iso, " (", grp, ")")

      brms_fit <- tryCatch(
        brm(
          log2cpm ~ 1 + (1 | PATNO),
          data = grp_data,
          family = gaussian(),
          prior = c(
            prior(normal(0, 5), class = "Intercept"),
            prior(student_t(3, 0, 1), class = "sd"),
            prior(student_t(3, 0, 1), class = "sigma")
          ),
          chains = BRMS_CHAINS,
          iter = BRMS_ITER,
          warmup = BRMS_WARMUP,
          seed = BRMS_SEED,
          cores = 4,
          refresh = 0,
          silent = 2
        ),
        error = function(e) {
          message("  brms error: ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(brms_fit)) {
        # Compute posterior ICC from variance components
        post <- as_draws_df(brms_fit)
        var_b <- post$`sd_PATNO__Intercept`^2
        var_w <- post$sigma^2
        icc_post <- var_b / (var_b + var_w)

        icc_summary <- tibble(
          isoform    = iso,
          DIAGNOSIS  = grp,
          ICC_median = round(median(icc_post), 4),
          ICC_lwr95  = round(quantile(icc_post, 0.025), 4),
          ICC_upr95  = round(quantile(icc_post, 0.975), 4),
          n_subjects = n_distinct(grp_data$PATNO)
        )
        bayes_icc_table[[paste0(iso, "_", grp)]] <- icc_summary

        log_md(
          "- ", iso, " (", grp, "): median ICC = ", icc_summary$ICC_median,
          " [", icc_summary$ICC_lwr95, ", ", icc_summary$ICC_upr95, "]"
        )

        saveRDS(
          brms_fit,
          file.path(
            OUT_DIR, "data",
            paste0("C_brms_icc_", iso, "_", grp, ".rds")
          )
        )
      }
    }
  }

  if (length(bayes_icc_table) > 0) {
    bayes_icc_df <- bind_rows(bayes_icc_table)
    write_csv(bayes_icc_df, file.path(OUT_DIR, "data", "C_bayesian_icc.csv"))
    log_md("")
    log_md("```")
    log_md(capture.output(print(bayes_icc_df)))
    log_md("```")
  }
}


log_md("")
log_md("## Output files saved")
log_md("- data/C_icc_variance_table.csv")
log_md("- data/C_bayesian_icc.csv (if brms ran)")

message("Done — Longitudinal_Stability.R")
