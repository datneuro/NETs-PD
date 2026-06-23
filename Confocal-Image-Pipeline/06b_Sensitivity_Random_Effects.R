#!/usr/bin/env Rscript
# 06b — Field-level (no case random effect) Bayesian + frequentist. Adds CrIs for plot annotation.
# Output: 03_results/03b_RE_sensitivity.csv  (est, CrI, P(dir), freq_p cho field-level)
suppressPackageStartupMessages({library(tidyverse); library(lme4); library(lmerTest)})
has_brms <- requireNamespace("brms", quietly=TRUE); if(has_brms) suppressPackageStartupMessages(library(brms))
RS  <- "results_stats"
OUT <- file.path(RS,"06_results")
SN_2X <- c("T243-10-2__scene01","T243-10-2__scene03","T243-10-2__scene08",
           "T243-10__scene01","T243-10__scene03","T243-10__scene06","T243-10__scene08")
fov <- read_csv(file.path(RS,"fov_level.csv"),show_col_types=FALSE) %>% filter(region=="SN")
net <- read_csv(file.path(RS,"net_level.csv"),show_col_types=FALSE) %>% filter(region=="SN") %>% mutate(logvol=log(volume_um3))
fov$condition <- factor(fov$condition,levels=c("HC","PD")); net$condition <- factor(net$condition,levels=c("HC","PD"))

bstat <- function(fit,col,gt0=TRUE){v<-as_draws_df(fit)[[col]]
  tibble(brms_est=median(v), CrI_lo=as.numeric(quantile(v,.025)), CrI_hi=as.numeric(quantile(v,.975)),
         bayes_Pdir=if(gt0) mean(v>0) else mean(v<0))}
rows <- list()

## COUNT (field-level)
p0 <- summary(glm(n_net_real~condition,data=fov,family=poisson))$coef["conditionPD",4]
b <- brm(n_net_real~condition,data=fov,family=poisson(),prior=prior(normal(0,1),class="b"),
         chains=4,iter=2000,refresh=0,seed=1,control=list(adapt_delta=0.99))
rows[["count_field"]] <- tibble(metric="count/FOV", level="field", freq_p=p0) %>% bind_cols(bstat(b,"b_conditionPD"))
## VOLUME (lesion-level)
p0 <- summary(lm(logvol~condition,data=net))$coef["conditionPD",4]
b <- brm(logvol~condition,data=net,family=gaussian(),prior=prior(normal(0,1),class="b"),
         chains=4,iter=2000,refresh=0,seed=1,control=list(adapt_delta=0.99))
rows[["volume_field"]] <- tibble(metric="volume/NET", level="lesion", freq_p=p0) %>% bind_cols(bstat(b,"b_conditionPD"))
## SPATIAL (field-level)
sp_path <- file.path(RS,"master_spatial_events.csv")
if(file.exists(sp_path)){
  sp <- read_csv(sp_path,show_col_types=FALSE) %>% filter(!image_id %in% SN_2X)
  fovm <- sp %>% group_by(image_id,case_id,event_type) %>% summarise(md=mean(dist_to_nearest_lewy_um),.groups="drop") %>%
    pivot_wider(names_from=event_type,values_from=md) %>% filter(!is.na(True_NET),!is.na(Random_Point)) %>% mutate(diff=True_NET-Random_Point)
  p0 <- t.test(fovm$diff)$p.value
  b <- brm(diff~1,data=fovm,family=gaussian(),chains=4,iter=2000,refresh=0,seed=1,control=list(adapt_delta=0.99))
  rows[["spatial_field"]] <- tibble(metric="NET-Lewy dist", level="field", freq_p=p0) %>% bind_cols(bstat(b,"b_Intercept",gt0=FALSE))
}
out <- bind_rows(rows); print(out); write_csv(out, file.path(OUT,"03b_RE_sensitivity.csv"))
cat("\nDONE 03b\n")
