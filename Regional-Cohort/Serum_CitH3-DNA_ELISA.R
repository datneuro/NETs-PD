# =============================================================================
# Serum_CitH3-DNA_ELISA.R
# Serum CitH3-DNA complex ELISA (NET-specific), HC vs PD + severity (Figure 3b)
# LM + emmeans + HC3 robust SE + stratified bootstrap + Bayesian (brms).
# Usage:  Rscript Serum_CitH3-DNA_ELISA.R [DATA_FILE] [OUT_DIR] [OUTCOME_VAR] [RUN_BRMS]
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(emmeans)
  library(sandwich)
  library(lmtest)
  library(boot)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(brms)
})

# -- CLI args -----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
DATA_FILE <- if (length(args) >= 1) {
  args[1]
} else {
  file.path("data", "CitH3_final.xlsx")
}

OUT_DIR <- if (length(args) >= 2) {
  args[2]
} else {
  "results"
}

OUTCOME_VAR <- if (length(args) >= 3) args[3] else "lg10CitH3"

RUN_BRMS <- if (length(args) >= 4) {
  !tolower(args[4]) %in% c("0", "false", "no", "off", "skip")
} else {
  TRUE
}
BRMS_CHAINS <- 4
BRMS_ITER <- 4000
BRMS_WARMUP <- 2000
BRMS_SEED <- 42

set.seed(42)

dir.create(file.path(OUT_DIR, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "png"), recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "CitH3_E_serum_log.md")
log_md <- function(...) cat(paste0(..., "\n"), file = LOG_FILE, append = TRUE)
cat("", file = LOG_FILE)
log_md("# Serum ELISA Analysis Log - Serum_CitH3-DNA_ELISA.R")
log_md("Date: ", format(Sys.time()))
log_md("Data: ", DATA_FILE)
log_md("Outcome: ", OUTCOME_VAR)
log_md("RUN_BRMS: ", RUN_BRMS)
log_md("")

custom_colors <- c("HC" = "#648FFF", "Mild PD" = "#785EF0", "Moderate PD" = "#DC267F", "PD" = "#DC267F")

fmt_p <- function(p) {
  if (is.na(p)) {
    return("NA")
  }
  if (p < 1e-4) {
    return(format(p, scientific = TRUE, digits = 2))
  }
  sprintf("%.4f", p)
}

bootstrap_lm_term <- function(data, formula_obj, term_name, R = 5000, strata = NULL) {
  vars_needed <- all.vars(formula_obj)
  d <- data %>% dplyr::select(any_of(vars_needed))
  cc <- complete.cases(d)
  d <- d[cc, , drop = FALSE]
  if (nrow(d) < 20) {
    return(tibble(term = term_name, estimate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  }

  if (is.null(strata)) {
    strata_d <- factor(rep("all", nrow(d)))
  } else {
    if (length(strata) == nrow(data)) {
      strata_d <- as.factor(strata[cc])
    } else if (length(strata) == nrow(d)) {
      strata_d <- as.factor(strata)
    } else {
      strata_d <- factor(rep("all", nrow(d)))
    }
  }

  idx_by_group <- split(seq_len(nrow(d)), strata_d)

  boot_est <- replicate(R, {
    boot_idx <- unlist(lapply(idx_by_group, function(idx) sample(idx, size = length(idx), replace = TRUE)))
    m <- tryCatch(lm(formula_obj, data = d[boot_idx, , drop = FALSE]), error = function(e) NULL)
    if (is.null(m)) {
      return(NA_real_)
    }
    co <- coef(m)
    if (!term_name %in% names(co)) {
      return(NA_real_)
    }
    as.numeric(co[[term_name]])
  })

  boot_t <- as.numeric(boot_est)
  boot_t <- boot_t[is.finite(boot_t)]

  if (length(boot_t) == 0) {
    return(tibble(term = term_name, estimate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  }

  ci <- as.numeric(quantile(boot_t, probs = c(0.025, 0.975), na.rm = TRUE))

  tibble(
    term = term_name,
    estimate = mean(boot_t, na.rm = TRUE),
    ci_lo = ci[1],
    ci_hi = ci[2]
  )
}

resolve_term_name <- function(fit, wanted = "diagPD") {
  term_names <- names(coef(fit))
  if (wanted %in% term_names) {
    return(wanted)
  }
  diag_terms <- grep("^diag", term_names, value = TRUE)
  if (length(diag_terms) > 0) {
    return(diag_terms[1])
  }
  NA_character_
}

safe_pbcor <- function(x, y) {
  if (!requireNamespace("WRS2", quietly = TRUE)) {
    return(list(r = NA_real_, p = NA_real_, status = "WRS2 not installed"))
  }
  out <- tryCatch(WRS2::pbcor(x, y), error = function(e) NULL)
  if (is.null(out)) {
    return(list(r = NA_real_, p = NA_real_, status = "pbcor failed"))
  }
  list(r = as.numeric(out$cor), p = as.numeric(out$p.value), status = "ok")
}

# =============================================================================
# Load and preprocess data
# =============================================================================
message("Loading serum ELISA data...")

if (!file.exists(DATA_FILE)) {
  stop("Input file not found: ", DATA_FILE)
}

serum_raw <- readxl::read_excel(DATA_FILE)

required_cols <- c("diag", "age", "sex", "subDX", OUTCOME_VAR)
missing_cols <- setdiff(required_cols, names(serum_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

serum <- serum_raw %>%
  mutate(
    diag = factor(diag, levels = c("HC", "PD")),
    subDX = factor(subDX, levels = c("HC", "Mild PD", "Moderate PD")),
    sex = factor(sex),
    age = as.numeric(age),
    outcome = as.numeric(.data[[OUTCOME_VAR]])
  )

log_md("## Data summary")
log_md("- n total: ", nrow(serum))
log_md("- n HC/PD with complete outcome: ",
       sum(!is.na(serum$outcome) & !is.na(serum$diag)))
log_md("- Outcome range: [", round(min(serum$outcome, na.rm = TRUE), 4), ", ",
       round(max(serum$outcome, na.rm = TRUE), 4), "]")
log_md("- Group counts (diag):")
log_md("```")
log_md(capture.output(print(table(serum$diag, useNA = "ifany"))))
log_md("```")
log_md("- Group counts (subDX):")
log_md("```")
log_md(capture.output(print(table(serum$subDX, useNA = "ifany"))))
log_md("```")
log_md("")

# =============================================================================
# 1) Between-group HC vs PD
# =============================================================================
message("Running HC vs PD model...")

df_group <- serum %>%
  filter(!is.na(outcome), !is.na(diag), !is.na(age), !is.na(sex))

fit_group <- lm(outcome ~ diag + age + sex, data = df_group)

anova_group <- anova(fit_group)
coef_group <- summary(fit_group)$coefficients
robust_group <- lmtest::coeftest(fit_group, vcov. = sandwich::vcovHC(fit_group, type = "HC3"))
emm_group <- emmeans(fit_group, pairwise ~ diag, adjust = "holm")
emm_group_df <- as.data.frame(emm_group$contrasts)

term_diag <- resolve_term_name(fit_group, wanted = "diagPD")

if (is.na(term_diag)) {
  stop("Cannot resolve diagnosis term in HC vs PD model coefficients.")
}
boot_group <- bootstrap_lm_term(
  data = df_group,
  formula_obj = outcome ~ diag + age + sex,
  term_name = term_diag,
  R = 5000,
  strata = df_group$diag
)

group_table <- broom::tidy(fit_group) %>%
  left_join(
    tibble(
      term = rownames(robust_group),
      estimate_hc3 = robust_group[, 1],
      se_hc3 = robust_group[, 2],
      stat_hc3 = robust_group[, 3],
      p_hc3 = robust_group[, 4]
    ),
    by = "term"
  )

write_csv(group_table, file.path(OUT_DIR, "data", "CitH3_E_lm_HC_vs_PD_coefficients.csv"))
write_csv(emm_group_df, file.path(OUT_DIR, "data", "CitH3_E_lm_HC_vs_PD_emmeans.csv"))
write_csv(boot_group, file.path(OUT_DIR, "data", "CitH3_E_lm_HC_vs_PD_bootstrap.csv"))

if (any(!is.finite(unlist(boot_group[, c("estimate", "ci_lo", "ci_hi")]))) ) {
  warning("Bootstrap HC vs PD produced NA/Inf. Check factor coding and resampling settings.")
}

log_md("## HC vs PD: LM + emmeans + HC3 + bootstrap")
log_md("- Model: outcome ~ diag + age + sex")
log_md("- n used: ", nrow(df_group))
log_md("### ANOVA")
log_md("```")
log_md(capture.output(print(anova_group)))
log_md("```")
log_md("### Coefficients (classical)")
log_md("```")
log_md(capture.output(print(coef_group)))
log_md("```")
log_md("### Coefficients (HC3 robust)")
log_md("```")
log_md(capture.output(print(robust_group)))
log_md("```")
log_md("### emmeans pairwise (Holm)")
log_md("```")
log_md(capture.output(print(emm_group_df)))
log_md("```")
log_md("### Bootstrap CI (term: ", term_diag, ")")
log_md("```")
log_md(capture.output(print(boot_group)))
log_md("```")
log_md("")

# =============================================================================
# 2) Severity model: HC vs Mild PD vs Moderate PD + polynomial trend
# =============================================================================
message("Running severity model...")

df_sev <- serum %>%
  filter(!is.na(outcome), !is.na(subDX), !is.na(age), !is.na(sex))

fit_sev <- lm(outcome ~ subDX + age + sex, data = df_sev)
coef_sev <- summary(fit_sev)$coefficients
robust_sev <- lmtest::coeftest(fit_sev, vcov. = sandwich::vcovHC(fit_sev, type = "HC3"))

emm_sev <- emmeans(fit_sev, pairwise ~ subDX, adjust = "holm")
emm_sev_df <- as.data.frame(emm_sev$contrasts)

poly_trend <- tryCatch(
  contrast(emmeans(fit_sev, ~ subDX), method = "poly"),
  error = function(e) NULL
)
poly_trend_df <- if (!is.null(poly_trend)) as.data.frame(poly_trend) else tibble()

write_csv(broom::tidy(fit_sev), file.path(OUT_DIR, "data", "CitH3_E_lm_severity_coefficients.csv"))
write_csv(emm_sev_df, file.path(OUT_DIR, "data", "CitH3_E_lm_severity_emmeans.csv"))
write_csv(poly_trend_df, file.path(OUT_DIR, "data", "CitH3_E_lm_severity_polytrend.csv"))

log_md("## Severity: LM + emmeans + trend")
log_md("- Model: outcome ~ subDX + age + sex")
log_md("- n used: ", nrow(df_sev))
log_md("### Coefficients (classical)")
log_md("```")
log_md(capture.output(print(coef_sev)))
log_md("```")
log_md("### Coefficients (HC3 robust)")
log_md("```")
log_md(capture.output(print(robust_sev)))
log_md("```")
log_md("### emmeans pairwise (Holm)")
log_md("```")
log_md(capture.output(print(emm_sev_df)))
log_md("```")
if (nrow(poly_trend_df) > 0) {
  log_md("### Polynomial trend (linear/quadratic)")
  log_md("```")
  log_md(capture.output(print(poly_trend_df)))
  log_md("```")
}
log_md("")

# =============================================================================
# 3) Clinical correlations (PD only)
# =============================================================================
message("Running clinical correlations...")

clinical_vars <- c("HY", "updrs1", "updrs2", "updrs3", "updrst")
clinical_vars <- clinical_vars[clinical_vars %in% names(serum)]

pd_only <- serum %>% filter(diag == "PD")

boot_pearson_ci <- function(x, y, R = 5000) {
  d <- tibble(x = x, y = y) %>% filter(!is.na(x), !is.na(y))
  if (nrow(d) < 10) {
    return(c(NA_real_, NA_real_))
  }
  stat_fun <- function(df, i) cor(df$x[i], df$y[i], method = "pearson")
  b <- boot::boot(d, stat_fun, R = R)
  ci <- tryCatch(boot::boot.ci(b, type = "perc")$percent[4:5], error = function(e) c(NA_real_, NA_real_))
  ci
}

corr_results <- lapply(clinical_vars, function(v) {
  d <- pd_only %>%
    transmute(outcome = outcome, clinical = as.numeric(.data[[v]])) %>%
    filter(!is.na(outcome), !is.na(clinical))

  if (nrow(d) < 10) {
    return(tibble(
      variable = v, n = nrow(d),
      pbcor_r = NA_real_, pbcor_p = NA_real_, pbcor_status = "n<10",
      spearman_r = NA_real_, spearman_p = NA_real_,
      pearson_r = NA_real_, pearson_ci_lo = NA_real_, pearson_ci_hi = NA_real_
    ))
  }

  pb <- safe_pbcor(d$outcome, d$clinical)
  sp <- suppressWarnings(cor.test(d$outcome, d$clinical, method = "spearman", exact = FALSE))
  pe <- suppressWarnings(cor.test(d$outcome, d$clinical, method = "pearson"))
  ci <- boot_pearson_ci(d$outcome, d$clinical, R = 5000)

  tibble(
    variable = v,
    n = nrow(d),
    pbcor_r = pb$r,
    pbcor_p = pb$p,
    pbcor_status = pb$status,
    spearman_r = as.numeric(sp$estimate),
    spearman_p = as.numeric(sp$p.value),
    pearson_r = as.numeric(pe$estimate),
    pearson_ci_lo = ci[1],
    pearson_ci_hi = ci[2]
  )
}) %>% bind_rows()

write_csv(corr_results, file.path(OUT_DIR, "data", "CitH3_E_correlation_PD_clinical.csv"))

log_md("## Correlations with clinical variables (PD only)")
log_md("```")
log_md(capture.output(print(corr_results)))
log_md("```")
log_md("")

# =============================================================================
# 4) Bayesian sensitivity
# =============================================================================
if (RUN_BRMS) {
  message("Running Bayesian models (brms)...")
  log_md("## Bayesian sensitivity")

  # 4a. HC vs PD model
  fit_group_brm <- tryCatch(
    brm(
      outcome ~ diag + age + sex,
      data = df_group,
      family = gaussian(),
      prior = c(
        prior(normal(0, 1), class = "b"),
        prior(student_t(3, 0, 2.5), class = "Intercept"),
        prior(student_t(3, 0, 2.5), class = "sigma")
      ),
      chains = BRMS_CHAINS,
      iter = BRMS_ITER,
      warmup = BRMS_WARMUP,
      seed = BRMS_SEED,
      cores = 4,
      refresh = 0,
      silent = 2
    ),
    error = function(e) NULL
  )

  if (!is.null(fit_group_brm)) {
    saveRDS(fit_group_brm, file.path(OUT_DIR, "data", "CitH3_E_brms_HC_vs_PD.rds"))
    post <- as_draws_df(fit_group_brm)
    diag_coef_name <- grep("^b_diag", names(post), value = TRUE)[1]
    p_pd_lt_hc <- mean(post[[diag_coef_name]] < 0)

    fixed_tbl <- as.data.frame(summary(fit_group_brm)$fixed) %>%
      rownames_to_column("term")
    write_csv(fixed_tbl, file.path(OUT_DIR, "data", "CitH3_E_brms_HC_vs_PD_fixed.csv"))

    log_md("### Bayesian HC vs PD")
    log_md("- coefficient term: ", diag_coef_name)
    log_md("- Posterior P(PD < HC): ", round(p_pd_lt_hc, 4))
    log_md("```")
    log_md(capture.output(print(fixed_tbl)))
    log_md("```")
  } else {
    log_md("### Bayesian HC vs PD: brms model failed")
  }

  # 4b. Bayesian clinical association in PD: outcome_clinical ~ outcome + age
  bayes_clinical <- list()
  for (v in clinical_vars) {
    d <- pd_only %>%
      transmute(clinical = as.numeric(.data[[v]]), outcome = outcome, age = age) %>%
      filter(!is.na(clinical), !is.na(outcome), !is.na(age))

    if (nrow(d) < 15) {
      next
    }

    d <- d %>% mutate(
      outcome_z = as.numeric(scale(outcome)),
      age_z = as.numeric(scale(age))
    )

    fit_b <- tryCatch(
      brm(
        clinical ~ outcome_z + age_z,
        data = d,
        family = gaussian(),
        prior = c(
          prior(normal(0, 1), class = "b"),
          prior(student_t(3, 0, 2.5), class = "Intercept"),
          prior(student_t(3, 0, 2.5), class = "sigma")
        ),
        chains = BRMS_CHAINS,
        iter = BRMS_ITER,
        warmup = BRMS_WARMUP,
        seed = BRMS_SEED,
        cores = 4,
        refresh = 0,
        silent = 2
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_b)) {
      saveRDS(fit_b, file.path(OUT_DIR, "data", paste0("CitH3_E_brms_PD_", v, ".rds")))
      post_b <- as_draws_df(fit_b)
      p_pos <- mean(post_b$b_outcome_z > 0)
      fsum <- as.data.frame(summary(fit_b)$fixed) %>% rownames_to_column("term")
      bet <- fsum %>% filter(term == "outcome_z")

      bayes_clinical[[v]] <- tibble(
        variable = v,
        beta = bet$Estimate,
        l95 = bet$`l-95% CI`,
        u95 = bet$`u-95% CI`,
        post_prob_gt0 = p_pos,
        n = nrow(d)
      )
    }
  }

  if (length(bayes_clinical) > 0) {
    bayes_clinical_df <- bind_rows(bayes_clinical)
    write_csv(bayes_clinical_df, file.path(OUT_DIR, "data", "CitH3_E_brms_PD_clinical_summary.csv"))
    log_md("### Bayesian clinical models (PD only)")
    log_md("```")
    log_md(capture.output(print(bayes_clinical_df)))
    log_md("```")
  }

  log_md("")
}

# =============================================================================
# 5) Figures
# =============================================================================
message("Saving figures...")

fig_group <- ggplot(serum %>% filter(!is.na(diag), !is.na(outcome)),
                    aes(x = diag, y = outcome, color = diag, fill = diag)) +
  geom_violin(alpha = 0.15, trim = FALSE, linewidth = 0.2, adjust = 1.5) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.35, linewidth = 0.4) +
  geom_jitter(width = 0.08, alpha = 0.65, size = 1.8) +
  stat_compare_means(method = "wilcox.test") +
  scale_color_manual(values = custom_colors[c("HC", "PD")]) +
  scale_fill_manual(values = custom_colors[c("HC", "PD")]) +
  labs(
    title = "Serum CitH3-DNA by diagnosis",
    subtitle = paste0("Outcome: ", OUTCOME_VAR),
    x = "Group",
    y = OUTCOME_VAR
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

fig_sev <- ggplot(serum %>% filter(!is.na(subDX), !is.na(outcome)),
                  aes(x = subDX, y = outcome, color = subDX, fill = subDX)) +
  geom_violin(alpha = 0.15, trim = FALSE, linewidth = 0.2, adjust = 1.5) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.35, linewidth = 0.4) +
  geom_jitter(width = 0.08, alpha = 0.65, size = 1.8) +
  stat_compare_means(comparisons = list(c("HC", "Mild PD"), c("Mild PD", "Moderate PD"), c("HC", "Moderate PD")), method = "wilcox.test") +
  scale_color_manual(values = custom_colors[c("HC", "Mild PD", "Moderate PD")]) +
  scale_fill_manual(values = custom_colors[c("HC", "Mild PD", "Moderate PD")]) +
  labs(
    title = "Serum CitH3-DNA by severity",
    subtitle = paste0("Outcome: ", OUTCOME_VAR),
    x = "Severity",
    y = OUTCOME_VAR
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

fig_combined <- fig_group + fig_sev + plot_annotation(tag_levels = "A")

ggsave(file.path(OUT_DIR, "pdf", "CitH3_E_serum_group_panels.pdf"), fig_combined, width = 12, height = 5)
ggsave(file.path(OUT_DIR, "png", "CitH3_E_serum_group_panels.png"), fig_combined, width = 12, height = 5, dpi = 300)

# Correlation panel for key clinical outcomes if available
key_clinical <- intersect(c("updrst", "HY"), clinical_vars)
if (length(key_clinical) > 0) {
  corr_figs <- lapply(key_clinical, function(v) {
    d <- pd_only %>%
      transmute(clinical = as.numeric(.data[[v]]), outcome = outcome) %>%
      filter(!is.na(clinical), !is.na(outcome))
    if (nrow(d) < 10) {
      return(NULL)
    }
    sp <- suppressWarnings(cor.test(d$outcome, d$clinical, method = "spearman", exact = FALSE))
    ggplot(d, aes(x = outcome, y = clinical)) +
      geom_point(color = "#DC267F", alpha = 0.65, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "#DC267F", fill = "#DC267F", alpha = 0.15) +
      labs(
        title = paste0(OUTCOME_VAR, " vs ", v, " (PD only)"),
        subtitle = paste0("Spearman r = ", round(as.numeric(sp$estimate), 3),
                          ", p = ", fmt_p(sp$p.value), ", n = ", nrow(d)),
        x = OUTCOME_VAR,
        y = v
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(size = 8, face = "italic"))
  })
  corr_figs <- Filter(Negate(is.null), corr_figs)
  if (length(corr_figs) >= 1) {
    p_corr <- wrap_plots(corr_figs)
    ggsave(file.path(OUT_DIR, "pdf", "CitH3_E_serum_corr_key.pdf"), p_corr, width = 10, height = 4)
    ggsave(file.path(OUT_DIR, "png", "CitH3_E_serum_corr_key.png"), p_corr, width = 10, height = 4, dpi = 300)
  }
}

log_md("## Output files")
log_md("- data/E_lm_HC_vs_PD_coefficients.csv")
log_md("- data/E_lm_HC_vs_PD_emmeans.csv")
log_md("- data/E_lm_HC_vs_PD_bootstrap.csv")
log_md("- data/E_lm_severity_coefficients.csv")
log_md("- data/E_lm_severity_emmeans.csv")
log_md("- data/E_lm_severity_polytrend.csv")
log_md("- data/E_correlation_PD_clinical.csv")
log_md("- data/E_brms_*.rds / E_brms_*_summary.csv (if brms ran)")
log_md("- pdf/E_serum_group_panels.pdf")
log_md("- png/E_serum_group_panels.png")
log_md("- pdf/png E_serum_corr_key.* (if enough data)")

message("Done - Serum_CitH3-DNA_ELISA.R")
