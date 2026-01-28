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
#II. NETs isoforms dataframe ----

padi4.201.bl.df <- (filter(net_core2,geneID=="ENST00000375448.4",CLINICAL_EVENT == "BL")%>% 
                      mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("HC","Prodromal", "PD"))))


padi4.202.bl.df <- (filter(net_core2,geneID=="ENST00000375453.5",CLINICAL_EVENT == "BL")%>% 
                      mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("HC","Prodromal", "PD"))))

mpo.201.bl.df <- (filter(net_core2,geneID=="ENST00000225275.3",CLINICAL_EVENT == "BL")%>% 
                      mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("HC","Prodromal", "PD"))))


# III. Function Definition for Plots with ANOVA-Dunnett's test ----
analyze_and_plot_dunnett <- function(df, y_axis_label) {
  
  # --- 1. Data Prep ---
  df$DIAGNOSIS <- factor(df$DIAGNOSIS)
  if("HC" %in% levels(df$DIAGNOSIS)) {
    df$DIAGNOSIS <- relevel(df$DIAGNOSIS, ref = "HC")
  }
  
  # --- 2. Statistical Analysis ---
  model_aov <- aov(exprss ~ DIAGNOSIS, data = df)
  dunnett_test <- glht(model_aov, linfct = mcp(DIAGNOSIS = "Dunnett"))
  summ <- summary(dunnett_test)
  
  # --- 3. Create Clean Stats Table ---
  stats_table <- data.frame(
    comparison = names(summ$test$coefficients),
    estimate   = as.numeric(summ$test$coefficients),
    p.adj      = as.numeric(summ$test$pvalues)
  ) %>%
    mutate(
      group1 = "HC",
      group2 = str_remove(comparison, " - HC"),
      significance = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        p.adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    )
  
  # --- 4. Generate Plot ---
  plot_stats <- stats_table %>% filter(significance != "ns")
  
  if(nrow(plot_stats) > 0) {
    max_y <- max(df$exprss, na.rm = TRUE)
    plot_stats <- plot_stats %>%
      mutate(y.position = max_y + 1.5 + (row_number() - 1) * 0.5)
  }
  
  p <- df %>%
    tidyplot(x = DIAGNOSIS, y = exprss, color = DIAGNOSIS) %>%
    add_boxplot(alpha = 0.5, show_outliers = FALSE, show_whiskers = FALSE) %>%
    add_violin(alpha = 0.2) %>%
    adjust_y_axis_title("Total MPO expression", face = "bold") %>%
    adjust_size(width = 125, height = 100) %>%
    adjust_font(fontsize = 15, family = "Arial") %>%
    remove_caption() %>%
    adjust_colors(new_colors = custom_colors2) %>%
    remove_legend() %>%
    remove_x_axis_title() %>%
    add_data_points_jitter(size = 2.5, alpha = 0.5) +
    labs(y = y_axis_label)
  
  if(nrow(plot_stats) > 0) {
    p <- p + stat_pvalue_manual(
      plot_stats, label = "significance", tip.length = 0.01,
      size = 5.5, bracket.size = 0.6, fontface = "bold", inherit.aes = FALSE
    )
  }
  
  # --- 5. AUTO-PRINTING (Now includes ANOVA) ---
  cat("\n=========================================\n")
  cat("          ANOVA SUMMARY (Global)         \n")
  cat("=========================================\n")
  print(summary(model_aov))  # <--- PRINTS ANOVA HERE
  
  cat("\n\n=========================================\n")
  cat("      DUNNETT'S TEST (vs Control)        \n")
  cat("=========================================\n")
  print(stats_table)         # <--- PRINTS DUNNETT HERE
  
  print(p)                   # <--- SHOWS PLOT
  
  # Return list invisibly
  return(invisible(list(plot = p, dunnett = stats_table, anova = model_aov)))
}

#IV. Execute above function ---- 
label_201 <- expression(bold("PADI4-201 mRNA (log"["2"] * "CPM)"))
analyze_and_plot_dunnett(padi4.201.bl.df, label_201)

label_202 <- expression(bold("PADI4-202 mRNA (log"["2"] * "CPM)"))
analyze_and_plot_dunnett(padi4.202.bl.df, label_202)

label_mpo201 <- expression(bold("MPO-201 mRNA (log"["2"] * "CPM)"))
mpo201.bl <- analyze_and_plot_dunnett(mpo.201.bl.df, label_mpo201)



#V. Baseline clinical correlation (Supplementary Fig. S1) ---- 
net_core4_bl <- net_core3 %>%
  filter(CLINICAL_EVENT == "BL", isoform_id %in% c("PADI4-202","PADI4-201", "MPO-201")) %>% 
  dplyr::select(PATNO, sample_id, DIAGNOSIS, isoform_id, exprss) %>%
  pivot_wider(names_from = isoform_id, values_from = exprss)


PAR.UPDRS$Duration_Days <- as.numeric(PAR.UPDRS$INFODT - PAR.UPDRS$PD_Diagnosis_Date)
PAR.UPDRS$Duration_Years <- PAR.UPDRS$Duration_Days / 365.25

PAR.UPDRS2.BL <- PAR.UPDRS %>%
  filter(EVENT_ID == "BL") %>%
  dplyr::select(PATNO,Duration_Years, ENROLL_AGE, HY_Stage, Part1, Part2, Part3, TOTAL_UPDRS)



net_core5 <- net_core4_bl %>% 
  left_join(PAR.UPDRS2.BL, by = "PATNO") %>% 
  dplyr::select(PATNO, sample_id, DIAGNOSIS, ENROLL_AGE,Duration_Years,
                HY_Stage, Part1, Part2, Part3, TOTAL_UPDRS, 
                `PADI4-201`, `PADI4-202`, `MPO-201`)
net_core5$DIAGNOSIS <- factor(net_core5$DIAGNOSIS, levels = c("Control", "Prodromal", "PD"))



netCORE_CORR <- net_core5 %>% 
  dplyr::select(-PATNO, -sample_id) %>%
  rename(
    "Duration" = Duration_Years,
    "Age" = ENROLL_AGE,
    "H&Y Stage" = HY_Stage,
    "MDS-UPDRS Part 1" = Part1,
    "MDS-UPDRS Part 2" = Part2,
    "MDS-UPDRS Part 3" = Part3,
    "Total MDS-UPDRS" = TOTAL_UPDRS,
  )


ggpairs(
  netCORE_CORR,
  columns = 2:11, 
  aes(colour = `DIAGNOSIS`, alpha = 0.8),
  lower = list(continuous="smooth"),
  upper = list(continuous = GGally::wrap("cor", method="pearson", stars=TRUE, digits=4)),
  diag=list(continuous="densityDiag")) +
  scale_color_manual(values = unique(custom_colors2)) + 
  scale_fill_manual(values = unique(custom_colors2))