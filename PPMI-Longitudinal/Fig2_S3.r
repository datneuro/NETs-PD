#Library ----
library(tidyverse)
library(edgeR)
library(tidyplots)
library(multcomp)
library(ggplot2)
library(stringr)
library(ggpubr) 
library(GGally)
library(psych)
library(lme4)

custom_colors2 <- c(
  "HC"  = "#648FFF",   
  "Prodromal" = "#785EF0", 
  "PD"   = "#DC267F"   
)

#I. Load DATA of baseline ---- 

load("df_normed_filtered_annotated.RData")
NAME.log2.cpm.filtered.norm.df <- NAME.log2.cpm.filtered.norm.df %>% dplyr::select (-last_col())
all_data_log2 <- NAME.log2.cpm.filtered.norm.df
net_core <- all_data_log2 %>% dplyr::filter(external_gene_name %in% c("MPO", "ELANE", "PADI4")) %>% 
  pivot_longer(cols = '3174_V08':last_col(),
               names_to = "sample_id",
               values_to = "expression")
samples<- read_csv("metaDataIR3.csv") %>%
  dplyr::filter(DIAGNOSIS %in% c("PD","Prodromal", "Control")) 

samples$sample_id <- paste0(samples$PATNO, "_", samples$CLINICAL_EVENT)
net_core2 <- samples %>% 
  dplyr::select (sample_id, PATNO, CLINICAL_EVENT, DIAGNOSIS, GENDER) %>% 
  right_join(net_core, by = "sample_id")
net_core2 <- net_core2 %>% rename("exprss" = "expression")
net_core2$DIAGNOSIS[net_core2$DIAGNOSIS=="Control"] <- "HC"
net_core3 <- net_core2 %>%
  mutate(isoform_id = case_when(
    geneID == "ENST00000225275.3" ~ "MPO-201",
    geneID == "ENST00000375448.4" ~ "PADI4-201",
    geneID == "ENST00000375453.5" ~ "PADI4-202",
    TRUE ~ NA_character_ 
  ))


# II. Function definition for LMM analysis (Individual detais) plots ----
analyze_interaction_separate_plots <- function(target_gene_id, data) {
  gene_data_combined <- data %>%
    filter(geneID == target_gene_id, DIAGNOSIS %in% c("PD", "Control")) %>%
    mutate(MONTH_FOLLOW = case_when(
      CLINICAL_EVENT == "BL"  ~ "0 month",
      CLINICAL_EVENT == "V02" ~ "6 months",
      CLINICAL_EVENT == "V04" ~ "12 months",
      CLINICAL_EVENT == "V06" ~ "24 months",
      TRUE                    ~ "36 months"
    )) %>%
    mutate(
      MONTH_FOLLOW = factor(MONTH_FOLLOW, levels = c("0 month", "6 months", "12 months", "24 months", "36 months")),
      DIAGNOSIS = factor(DIAGNOSIS, levels = c("Control", "PD"))
    )
  
  if (nrow(gene_data_combined) < 10) return(NULL)
  common_gene_name <- gene_data_combined$isoform_id[1]
  interaction_model <- lmer(exprss ~ MONTH_FOLLOW * DIAGNOSIS + (1 | PATNO), data = gene_data_combined)
  anova_results <- anova(interaction_model)
  time_res <- anova_results["MONTH_FOLLOW",]
  diag_res <- anova_results["DIAGNOSIS",]
  int_res  <- anova_results["MONTH_FOLLOW:DIAGNOSIS",]
  
  fmt_p <- function(p) sprintf("%.2e", p)
  
  line1 <- sprintf(
    "Time effect: F(%.1f, %.1f)=%.2f, p=%s  |  Diagnosis: F(%.1f, %.1f)=%.2f, p=%s", 
    time_res$NumDF, time_res$DenDF, time_res$`F value`, fmt_p(time_res$`Pr(>F)`),
    diag_res$NumDF, diag_res$DenDF, diag_res$`F value`, fmt_p(diag_res$`Pr(>F)`)
  )
  
  line2 <- sprintf(
    "Time effect × Diagnosis: F(%.1f, %.1f)=%.2f, p=%s", 
    int_res$NumDF, int_res$DenDF, int_res$`F value`, fmt_p(int_res$`Pr(>F)`)
  )
  
  full_subtitle <- paste(line1, line2, sep = "\n")
  

  posthoc_results <- emmeans(interaction_model, specs = trt.vs.ctrl ~ MONTH_FOLLOW | DIAGNOSIS, ref = 1, adjust = "dunnett")
  

  create_single_plot <- function(subset_data, diagnosis_name, contrast_results, color_palette_func, show_stats_subtitle) {
    
    max_y <- max(subset_data$exprss, na.rm = TRUE)
    
    stats_df <- as.data.frame(contrast_results$contrasts) %>%
      filter(DIAGNOSIS == diagnosis_name) %>%
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
      filter(p_star != "ns") 
    
    if(nrow(stats_df) > 0) {
      stats_df <- stats_df %>%
        mutate(y.position = seq(from = max_y + (0.05 * max_y), by = 0.4, length.out = n()))
    }
    
    sample_sizes <- subset_data %>%
      group_by(MONTH_FOLLOW) %>%
      summarise(n = n(), .groups = 'drop')
    x_labels_n <- paste0(sample_sizes$MONTH_FOLLOW, "\n(n=", sample_sizes$n, ")")
    
    if (show_stats_subtitle) {
      final_subtitle <- full_subtitle
    } else {
      final_subtitle <- NULL 
    }
    
    p <- ggplot(subset_data, aes(x = MONTH_FOLLOW, y = exprss)) +
      geom_line(aes(group = PATNO), color = "grey85", alpha = 0.5) +
      geom_boxplot(aes(fill = MONTH_FOLLOW), outlier.shape = NA, alpha = 0.5, show.legend = FALSE) +
      geom_jitter(aes(color = MONTH_FOLLOW), width = 0.15, alpha = 0.4, size = 2, show.legend = FALSE) +
      
      color_palette_func() + 
      scale_color_bmj() + 
      
      scale_x_discrete(labels = x_labels_n) +
      labs(
        title = paste0(common_gene_name, " in ", diagnosis_name),
        subtitle = final_subtitle,
        y = bquote(bold(.(common_gene_name) ~ "mRNA (log"["2"] * "CPM)")),
        x = ""
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        # Adjusted size slightly down to fit the long subtitle
        plot.subtitle = element_text(size = 9, face = "italic", lineheight = 1.2),
        axis.text.x = element_text(size = 10),
        panel.grid.minor = element_blank()
      )
    
    if(nrow(stats_df) > 0) {
      p <- p + stat_pvalue_manual(
        stats_df, label = "p_star", tip.length = 0.01, size = 6, 
        bracket.size = 0.6, inherit.aes = FALSE
      )
    }
    return(p)
  }
  

  plot_control <- create_single_plot(
    subset_data = gene_data_combined %>% filter(DIAGNOSIS == "Control"),
    diagnosis_name = "Control",
    contrast_results = posthoc_results,
    color_palette_func = scale_fill_bmj,
    show_stats_subtitle = TRUE 
  )
  
  plot_pd <- create_single_plot(
    subset_data = gene_data_combined %>% filter(DIAGNOSIS == "PD"),
    diagnosis_name = "PD",
    contrast_results = posthoc_results,
    color_palette_func = scale_fill_bmj,
    show_stats_subtitle = TRUE
  )
  
  return(list(plot_PD = plot_pd, plot_Control = plot_control, anova_table = anova_results))
}

# EXECUTION ----
PADI201_res <- analyze_interaction_separate_plots("ENST00000375448.4", net_core3)
print(PADI201_res$plot_Control) 
print(PADI201_res$plot_PD) 

PADI202_res <- analyze_interaction_separate_plots("ENST00000375453.5", net_core3)
print(PADI202_res$plot_Control) 
print(PADI202_res$plot_PD)    

MPO_res <- analyze_interaction_separate_plots("ENST00000225275.3", net_core3)
print(MPO_res$plot_PD)
print(MPO_res$plot_Control)


# Function definition for plots of predicted means ---- 
plot_predicted_trajectories_final <- function(target_gene_id, data) {
  
  # 1. Prepare Data
  gene_data <- data %>%
    filter(geneID == target_gene_id, DIAGNOSIS %in% c("PD", "HC")) %>%
    mutate(MONTH_NUM = case_when(
      CLINICAL_EVENT == "BL"  ~ 0,
      CLINICAL_EVENT == "V02" ~ 6,
      CLINICAL_EVENT == "V04" ~ 12,
      CLINICAL_EVENT == "V06" ~ 24,
      TRUE                    ~ 36
    )) %>%
    # FIXED: Ensure strings and levels match exactly (all plural "months" except 0)
    mutate(
      MONTH_LABEL = ifelse(MONTH_NUM == 0, "0 month", paste0(MONTH_NUM, " months")),
      MONTH_FOLLOW = factor(MONTH_LABEL, 
                            levels = c("0 month", "6 months", "12 months", "24 months", "36 months")),
      DIAGNOSIS = factor(DIAGNOSIS, levels = c("PD", "HC"))
    )
  
  common_gene_name <- gene_data$isoform_id[1]
  
  # --- 2. CALCULATE SAMPLE SIZES (n=...) ---
  sample_sizes <- gene_data %>%
    group_by(MONTH_FOLLOW) %>%
    summarise(n = n(), .groups = 'drop')
  
  x_labels_with_n <- setNames(paste0(sample_sizes$MONTH_FOLLOW, "\n(n=", sample_sizes$n, ")"), 
                              sample_sizes$MONTH_FOLLOW)
  
  # 3. Fit the Linear Mixed Model
  model <- lmer(exprss ~ MONTH_FOLLOW * DIAGNOSIS + (1 | PATNO), data = gene_data)
  
  # 4. Extract Predicted Means and CI
  pred_means <- emmeans(model, ~ MONTH_FOLLOW * DIAGNOSIS) %>% as.data.frame()
  
  # 5. Extract Stats with Exact Formatting
  anova_res <- anova(model)
  fmt_p_exact <- function(p) {
    if (p < 0.0001) return(sprintf("%.2e", p))
    return(sprintf("%.4f", p))
  }
  
  t_res <- anova_res["MONTH_FOLLOW", ]
  d_res <- anova_res["DIAGNOSIS", ]
  i_res <- anova_res["MONTH_FOLLOW:DIAGNOSIS", ]
  
  line1 <- sprintf(
    "Time effect: F(%.1f, %.1f)=%.2f, p=%s | Diagnosis: F(%.1f, %.1f)=%.2f, p=%s",
    t_res$NumDF, t_res$DenDF, t_res$`F value`, fmt_p_exact(t_res$`Pr(>F)`),
    d_res$NumDF, d_res$DenDF, d_res$`F value`, fmt_p_exact(d_res$`Pr(>F)`)
  )
  line2 <- sprintf(
    "Interaction (Time effect × Diagnosis): F(%.1f, %.1f)=%.2f, p=%s",
    i_res$NumDF, i_res$DenDF, i_res$`F value`, fmt_p_exact(i_res$`Pr(>F)`)
  )
  full_subtitle <- paste(line1, line2, sep = "\n")
  
  # 6. Create the Plot
  p <- ggplot(pred_means, aes(x = MONTH_FOLLOW, y = emmean, color = DIAGNOSIS, group = DIAGNOSIS)) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = DIAGNOSIS), alpha = 0.12, color = NA) +
    geom_line(size = 1.1) +
    geom_point(size = 2.5) +
    
    # Professional Aesthetic (High Contrast for Publication)
    scale_color_manual(values = c("HC" = "#648FFF", "PD" = "#DC267F")) + 
    scale_fill_manual(values = c("HC" = "#648FFF", "PD" = "#DC267F")) +
    scale_x_discrete(labels = x_labels_with_n) +
    
    labs(
      #title = paste0("Longitudinal Trajectory: ", common_gene_name),
      subtitle = full_subtitle,
      y = bquote(bold(.(common_gene_name) ~ "mRNA (log"["2"] * "CPM)")),
      x = "",
      color = "Diagnosis", fill = "Diagnosis"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 8.5, face = "italic", lineheight = 1.1),
      legend.position = "right", 
      axis.title = element_text(face = "bold", size = 11),
      axis.text.x = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Execution - Plots for longitudinal lines ---- 
fig_PADI4_201 <- plot_predicted_trajectories_final("ENST00000375448.4", net_core3)
fig_PADI4_202 <- plot_predicted_trajectories_final("ENST00000375453.5", net_core3)
fig_MPO_201  <- plot_predicted_trajectories_final("ENST00000225275.3", net_core3)
print(fig_PADI4_201)
print(fig_PADI4_202)
print(fig_MPO_201)