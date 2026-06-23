#!/usr/bin/env Rscript
# 06_NET_Statistics.R — Final statistics for the NET confocal pipeline.
# Input: results_stats/{fov_level.csv, net_level.csv, master_spatial_events.csv}
# Output: results_stats/06_results/*.csv  (summary tables used for figures)
# Approach: frequentist (lme4) + Bayesian primary (brms, weak priors) for small n.

suppressPackageStartupMessages({
  library(tidyverse); library(lme4); library(lmerTest)
})
has_brms <- requireNamespace("brms", quietly=TRUE); if(has_brms) suppressPackageStartupMessages(library(brms))

RS  <- "results_stats"
OUT <- file.path(RS, "06_results"); dir.create(OUT, showWarnings=FALSE)
SN_2X_EXCLUDE <- c("T243-10-2__scene01","T243-10-2__scene03","T243-10-2__scene08",
                   "T243-10__scene01","T243-10__scene03","T243-10__scene06","T243-10__scene08")

fov <- read_csv(file.path(RS,"fov_level.csv"), show_col_types=FALSE)
net <- read_csv(file.path(RS,"net_level.csv"), show_col_types=FALSE)
fov$condition <- factor(fov$condition, levels=c("HC","PD"))
net$condition <- factor(net$condition, levels=c("HC","PD"))

bayes_dir <- function(draws_col, fit, gt0=TRUE){
  d <- as_draws_df(fit); v <- d[[draws_col]]
  list(est=median(v), lo=quantile(v,.025), hi=quantile(v,.975),
       p_dir=if(gt0) mean(v>0) else mean(v<0))
}
res <- list()

## ───────────────────────── SN: HC vs PD ─────────────────────────
sn <- fov %>% filter(region=="SN")

# 1) ABUNDANCE: count/FOV 
cat("\n[1] SN abundance count/FOV\n")
g <- glmer(n_net_real ~ condition + (1|case_id), data=sn, family=poisson)
fe <- summary(g)$coefficients["conditionPD",]
row <- tibble(metric="count_per_FOV", model="glmer_poisson",
              estimate_logRR=fe[1], p_value=fe[4],
              mean_HC=mean(sn$n_net_real[sn$condition=="HC"]),
              mean_PD=mean(sn$n_net_real[sn$condition=="PD"]))
if(has_brms){
  b <- brm(n_net_real ~ condition + (1|case_id), data=sn, family=poisson(),
           prior=c(prior(normal(0,1),class="b"), prior(student_t(3,0,1),class="sd")),
           chains=4, iter=2000, refresh=0, seed=42,
           control=list(adapt_delta=0.99))
  bd <- bayes_dir("b_conditionPD", b)
  row <- row %>% mutate(brms_logRR=bd$est, CrI_lo=bd$lo, CrI_hi=bd$hi, P_PD_gt_HC=bd$p_dir)
}
res$count <- row; print(row)

# 2) SIZE: volume per-NET (log) — net-level, nested
cat("\n[2] SN volume per-NET\n")
sn_net <- net %>% filter(region=="SN") %>% mutate(logvol=log(volume_um3))
m <- lmer(logvol ~ condition + (1|case_id/image_id), data=sn_net)
fe <- summary(m)$coefficients["conditionPD",]
row <- tibble(metric="volume_perNET_log", model="lmer",
              estimate=fe[1], p_value=fe[5],
              median_HC=median(sn_net$volume_um3[sn_net$condition=="HC"]),
              median_PD=median(sn_net$volume_um3[sn_net$condition=="PD"]),
              n_HC=sum(sn_net$condition=="HC"), n_PD=sum(sn_net$condition=="PD"))
if(has_brms){
  b <- brm(logvol ~ condition + (1|case_id/image_id), data=sn_net, family=gaussian(),
           prior=c(prior(normal(0,1),class="b"), prior(student_t(3,0,1),class="sd")),
           chains=4, iter=2000, refresh=0, seed=42,
           control=list(adapt_delta=0.99))
  bd <- bayes_dir("b_conditionPD", b)
  row <- row %>% mutate(brms_est=bd$est, CrI_lo=bd$lo, CrI_hi=bd$hi, P_PD_gt_HC=bd$p_dir)
}
res$volume <- row; print(row)

# 3) SHAPE: solidity per-NET
cat("\n[3] SN solidity per-NET\n")
m <- lmer(solidity ~ condition + (1|case_id/image_id), data=sn_net)
fe <- summary(m)$coefficients["conditionPD",]
res$solidity <- tibble(metric="solidity_perNET", model="lmer", estimate=fe[1], p_value=fe[5],
              median_HC=median(sn_net$solidity[sn_net$condition=="HC"]),
              median_PD=median(sn_net$solidity[sn_net$condition=="PD"]))
print(res$solidity)

## ───────────────────────── SPATIAL (Bayesian) ─────────────────────────
cat("\n[4] Spatial NET<->Lewy (Bayesian)\n")
sp_path <- file.path(RS,"master_spatial_events.csv")
if(file.exists(sp_path)){
  sp <- read_csv(sp_path, show_col_types=FALSE) %>% filter(!image_id %in% SN_2X_EXCLUDE)
  fovm <- sp %>% group_by(image_id, case_id, event_type) %>%
    summarise(md=mean(dist_to_nearest_lewy_um), .groups="drop") %>%
    pivot_wider(names_from=event_type, values_from=md) %>%
    filter(!is.na(True_NET), !is.na(Random_Point)) %>%
    mutate(diff = True_NET - Random_Point)   # <0 = NET closer to Lewy than random
  row <- tibble(metric="NET_to_Lewy_dist_diff", n_FOV=nrow(fovm),
                mean_true=mean(fovm$True_NET), mean_random=mean(fovm$Random_Point),
                mean_diff=mean(fovm$diff))
  if(has_brms && nrow(fovm)>=3){
    b <- brm(diff ~ 1 + (1|case_id), data=fovm, family=gaussian(),
             prior=c(prior(student_t(3,0,10),class="Intercept"), prior(student_t(3,0,10),class="sd")),
             chains=4, iter=2000, refresh=0, seed=42,
           control=list(adapt_delta=0.99))
    bd <- bayes_dir("b_Intercept", b, gt0=FALSE)  # P(diff<0) = P(NET closer)
    row <- row %>% mutate(brms_diff=bd$est, CrI_lo=bd$lo, CrI_hi=bd$hi, P_NET_closer=bd$p_dir)
  }
  res$spatial <- row; print(row)
}

## ───────────────────────── CORTEX (descriptive) ─────────────────────────
cat("\n[5] Cortex (PD) descriptive\n")
ctx <- fov %>% filter(region=="cortex"); ctx_net <- net %>% filter(region=="cortex")
res$cortex <- tibble(metric="cortex_PD_descriptive", n_FOV=nrow(ctx),
  mean_count_FOV=mean(ctx$n_net_real), mean_per_mm2=mean(ctx$net_per_mm2),
  median_volume=median(ctx_net$volume_um3), median_solidity=median(ctx_net$solidity))
print(res$cortex)

## ───────────────────────── EXPORT ─────────────────────────
for(nm in names(res)) write_csv(res[[nm]], file.path(OUT, paste0("06_", nm, ".csv")))
# Combine into a multi-sheet Excel (stats + data points for figures)
if(requireNamespace("openxlsx", quietly=TRUE)){
  sheets <- res
  sheets[["DATA_fov_level"]]  <- fov
  sheets[["DATA_net_level"]]  <- net
  openxlsx::write.xlsx(sheets, file.path(OUT, "NET_confocal_results.xlsx"), overwrite=TRUE)
  cat("→ Excel: NET_confocal_results.xlsx\n")
}
cat("\n✅ DONE. Results ->", OUT, "\n")
