#1.Library ----
library(readxl)
library(dplyr)
library(tidyplots)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggstatsplot)
library(GGally)
library(finalfit)
library(EnvStats)
library(broom)

custom_colors2 <- c(
  "HC"  = "#648FFF",   
  "Mild PD" = "#785EF0", 
  "Moderate PD"   = "#DC267F"   
)

custom_colors3 <- c(
  "Control"  = "#2A6EBBFF",   
  "PD" = "#C50084FF", 
  "MSA"   = "#69BE28FF"
)

#2.Data import ----

mpo.dna.intergate <- read_excel("mpo_final_26Jan16.xlsx")
hc <- mpo.dna.intergate %>% filter (diag == "HC")
pd <- mpo.dna.intergate %>% filter (diag == "PD")

Varname <- "MPO-DNA"

#3. Serum NETs-ELISA (Figure 3a,b)---- 
fig3a<- mpo.dna.intergate |>
  tidyplot(x = diag, y = lg_MPO_DNA,color=diag) |>
  add_violin(alpha=0.2) |>
  add_boxplot(alpha=0.5) |>
  add_data_points_beeswarm(size = 2,alpha =0.8)|>
  add_test_asterisks(method = "wilcox_test",bracket.nudge.y = 0.5) |> 
  adjust_x_axis_title("Group")|>
  adjust_size (width = 100, height = 90) |> 
  adjust_colors(new_colors = custom_colors1)|> 
  adjust_legend_title("Group") |>
  adjust_font(fontsize=10, family = "Arial") + ggplot2::labs(y = bquote(.(Varname) ~ "(" * Log[10] ~ "Absorbance" * ")"))


fig3b <- mpo.dna.intergate |>
  tidyplot(x = subDX, y = lg_MPO_DNA,color=subDX) |>
  add_violin(alpha=0.2) |>
  add_boxplot(alpha=0.5) |>
  add_data_points_beeswarm(size = 2,alpha =0.8)|>
  add_test_asterisks(method = "wilcox_test", bracket.nudge.y = 0.5) |> 
  adjust_x_axis_title("Group")|>
  adjust_size (width = 100, height = 90) |> 
  adjust_colors(new_colors = custom_colors2)|> 
  adjust_legend_title("Group") |>
  adjust_font(fontsize=10, family = "Arial") + ggplot2::labs(y = bquote(.(Varname) ~ "(" * Log[10] ~ "Absorbance" * ")"))


#4. Correlation of NETs-ELISA and clinical (Fig 3c) ---- 
cor.var.1 <- mpo.dna.intergate %>% select (subDX, age, duration, HY, updrs1, 
                                           updrs2, updrs3, updrst, lg_MPO_DNA
                                           ) %>% 
  rename(`Age` = `age`,
         `Disease duration` = `duration`,
         `H&Y stage` = `HY`,
         `MDS-UPDRS part I` = `updrs1`,
         `MDS-UPDRS part II` = `updrs2`,
         `MDS-UPDRS part III` = `updrs3`,
         `MDS-UPDRS total score` = `updrst`,
         `MPO-DNA` = `lg_MPO_DNA`
         )


ggpairs(
  cor.var.1,
  columns = 2:9, 
  aes(colour = `subDX`, alpha = 0.8),
  lower = list(continuous="smooth"),
  upper = list(continuous = GGally::wrap("cor", method="pearson", stars=TRUE,
                                         na.rm = TRUE,digits=4)),
  diag=list(continuous="densityDiag")) +
  scale_color_manual(values = unique(custom_colors2)) + 
  scale_fill_manual(values = unique(custom_colors2))


#5. Western blot of MPO brain (Supplementary Figure S4) ----
mpo.wb <- read_excel("mpo_wb.xlsx")
mpo.wb$dx <- factor(mpo.wb$dx, levels = c("Control", "PD", "MSA"))


mpo.wb |> 
  tidyplot(x = dx, y = MPO, color = dx) |> 
  add_boxplot(alpha=0.5) |>
  add_violin(alpha=0.2) |> 
  add_test_asterisks(method = "t.test", 
                     bracket.nudge.y = 0.5,
                     label.size = 5.5) |> 
  adjust_y_axis_title("Total MPO expression", face="bold")|>
  adjust_size (width = 125, height = 100) |> 
  adjust_font(fontsize=15, family = "Arial")|> 
  remove_caption()|>
  remove_legend() |>
  remove_x_axis_title()|>
  adjust_colors(new_colors = custom_colors3)|> 
  add_data_points_jitter(size = 5,alpha =0.8) 



#6. ELISA data of NETs in brain (Fig.4) ----
net.brain <- read_excel("net.brain.xlsx")
net.brain$dx <- factor(net.brain$dx, levels = c("Control", "PD", "MSA"))


net.brain |> 
  tidyplot(x = dx, y = net, color = dx) |> 
  add_boxplot(alpha=0.5) |>
  add_violin(alpha=0.2) |>
  add_test_asterisks(method = "t.test", 
                     bracket.nudge.y = 0.5,
                     label.size = 5.5) |> 
  adjust_y_axis_title("NETs", face="bold")|>
  adjust_size (width = 100, height = 90) |> 
  adjust_font(fontsize=15, family = "Arial")|> 
  remove_caption()|>
  remove_legend() |>
  remove_x_axis_title()|>
  adjust_colors(new_colors = custom_colors3)|> 
  add_data_points_jitter(size = 5,alpha =0.8) 
