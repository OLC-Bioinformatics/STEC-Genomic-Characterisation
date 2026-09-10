#set working directory
setwd("/Users/noorshubair/Desktop/R_code_files/chapter2/")

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(readxl)



# KMA Only Comparison -----------------------------------------------------


kmadf <- read_excel("all_tools_subtype_data_forR.xlsx", sheet="kma_only_data")
kmadf

#make a new column that has all of our values concatenated
kmadf$Tool2 <- paste(kmadf$Tool,kmadf$ID,kmadf$Version, sep = "_")
kmadf

#remove the columns we dont need
kmadf2 <- kmadf[,-c(1:3)]
kmadf2


#get the metrics for KMA
#get the metrics for confidence intervals
kma_metrics <- kmadf2 %>%
  rowwise() %>%
  mutate(
    # Sensitivity - true positive rate
    Sens_Est = binom.test(TP, TP + FN)$estimate,
    Sens_Lwr = binom.test(TP, TP + FN)$conf.int[1],
    Sens_Upr = binom.test(TP, TP + FN)$conf.int[2],
    
    # Specificity - true negative rate
    Spec_Est = binom.test(TN, TN + FP)$estimate,
    Spec_Lwr = binom.test(TN, TN + FP)$conf.int[1],
    Spec_Upr = binom.test(TN, TN + FP)$conf.int[2],
    
    # PPV - precision
    PPV_Est  = binom.test(TP, TP + FP)$estimate,
    PPV_Lwr  = binom.test(TP, TP + FP)$conf.int[1],
    PPV_Upr  = binom.test(TP, TP + FP)$conf.int[2],
    
    # NPV 
    NPV_Est  = binom.test(TN, TN + FN)$estimate,
    NPV_Lwr  = binom.test(TN, TN + FN)$conf.int[1],
    NPV_Upr  = binom.test(TN, TN + FN)$conf.int[2],

    # Accuracy - proportion of all subtype predictions classified correctly
    Acc_Est = binom.test(TP + TN, TP + TN + FP + FN)$estimate,
    Acc_Lwr = binom.test(TP + TN, TP + TN + FP + FN)$conf.int[1],
    Acc_Upr = binom.test(TP + TN, TP + TN + FP + FN)$conf.int[2]
  ) %>%
  ungroup() %>%
  select(Tool2, contains("_"))

kma_metrics

#reshape to long format for ggplot
kmametrics_long <- kma_metrics %>%
  pivot_longer(
    cols = -Tool2,
    names_to = c("Metric", ".value"),
    names_pattern = "([A-Za-z]+)_(Est|Lwr|Upr)"
  ) %>%
  rename(Estimate = Est, Lower_CI = Lwr, Upper_CI = Upr) %>%
  mutate(Metric = case_when(
    Metric == "Sens" ~ "Sensitivity/Recall",
    Metric == "Spec" ~ "Specificity",
    Metric == "PPV" ~ "Precision (PPV)",
    Metric == "Acc" ~ "Accuracy",
    TRUE ~ Metric
  )) %>%
  mutate(Metric = factor(
    Metric,
    levels = c("NPV", "Sensitivity/Recall", "Precision (PPV)",
               "Specificity", "Accuracy")
  ))

kmametrics_long

#now we are going to split the Tool2 column back into the separate values..
kma_long <- separate_wider_delim(
  data = kmametrics_long,
  cols = Tool2,
  delim = "_",
  names = c("Tool","ID","Ver")
)

kma_long

# format ID legend labels to have 2 decimal places (instead of 1 for some)
kma_long <- kma_long %>%
  mutate(ID = sprintf("%.2f", as.numeric(ID)))

kma_long

#now we will plot it
kma_estimate_metrics <- ggplot(kma_long, 
                               aes(x = Metric, y = Estimate, group = interaction(Tool, ID), color = ID, fill = ID, shape = Ver)) +
  # Error bars show exact binomial 95% confidence intervals.
  geom_errorbar(
    aes(ymin = Lower_CI, ymax = Upper_CI, linetype = Ver),
    width = 0.5, linewidth = 0.8,
    position = position_dodge(0.69)
  ) +
  geom_point(position = position_dodge(0.7),size = 4, color = "black") +
  scale_shape_manual(values = c(22,21)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 15),
        legend.position = "top", legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), legend.key.width = unit(3, "line")
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "Estimated Value", color = "Identity Cutoff", fill = "Identity Cutoff", 
       shape = "ConClave Mode", linetype = "ConClave Mode") +
  guides(
    color = guide_legend(order = 1),
    fill = guide_legend(
      order = 1,
      override.aes = list(shape = 21)
    ),
    shape = guide_legend(order = 2, nrow = 2),
    linetype = guide_legend(order = 2, nrow = 2)
  )

kma_estimate_metrics

dpi = 300
jpeg("performance_KMA_only_plot_w_CI_allmetrics.jpg", width = 10*dpi, height = 8*dpi, res = dpi)
kma_estimate_metrics
dev.off()



# All Tools Comparison ----------------------------------------------------


alltoolsdf <- read_excel("all_tools_subtype_data_forR.xlsx", sheet="all_tools_data")
alltoolsdf


#get the metrics for confidence intervals
alltools_metrics <- alltoolsdf %>%
  rowwise() %>%
  mutate(
    # Sensitivity - true positive rate
    Sens_Est = binom.test(TP, TP + FN)$estimate,
    Sens_Lwr = binom.test(TP, TP + FN)$conf.int[1],
    Sens_Upr = binom.test(TP, TP + FN)$conf.int[2],
    
    # Specificity - true negative rate
    Spec_Est = binom.test(TN, TN + FP)$estimate,
    Spec_Lwr = binom.test(TN, TN + FP)$conf.int[1],
    Spec_Upr = binom.test(TN, TN + FP)$conf.int[2],
    
    # PPV - precision
    PPV_Est  = binom.test(TP, TP + FP)$estimate,
    PPV_Lwr  = binom.test(TP, TP + FP)$conf.int[1],
    PPV_Upr  = binom.test(TP, TP + FP)$conf.int[2],
    
    # NPV 
    NPV_Est  = binom.test(TN, TN + FN)$estimate,
    NPV_Lwr  = binom.test(TN, TN + FN)$conf.int[1],
    NPV_Upr  = binom.test(TN, TN + FN)$conf.int[2]
  ) %>%
  ungroup() %>%
  select(Tool, contains("_"))

alltools_metrics

#reshape to long format for ggplot
alltools_long <- alltools_metrics %>%
  pivot_longer(
  cols = -Tool,
  names_to = c("Metric", ".value"),
  names_pattern = "([A-Za-z]+)_(Est|Lwr|Upr)"
) %>%
  rename(Estimate = Est, Lower_CI = Lwr, Upper_CI = Upr) %>%
  mutate(Metric = case_when(
    Metric == "Sens" ~ "Sensitivity/Recall",
    Metric == "Spec" ~ "Specificity",
    Metric == "PPV" ~ "Precision (PPV)",
    TRUE ~ Metric
  ))

alltools_long

#plot it
performance_alltools_CI_plot <- ggplot(alltools_long, 
       aes(x = Metric, y = Estimate,group = Tool, color = Tool,fill = Tool, shape = Tool)) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.5, linewidth = 0.8,
                position = position_dodge(0.69)) +
  geom_point(position = position_dodge(0.7),size = 4, color = "black") +
  scale_shape_manual(values = c(22,21,24,25)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 15),
        legend.position = "top", legend.title = element_text(size = 12), legend.text = element_text(size = 10)
        ) +
  # Use manually selected, high-contrast colours.
  # KMA is orange so it is easier to distinguish from the white background.
  scale_color_manual(values = c(
    "GeneSeekr" = "#E69F00",
    "KMA 70% ConClave 1" = "#66C2A5",
    "StxTyper" = "#8DA0CB",
    "Sipprverse" = "#E76F51"
  )) +
  scale_fill_manual(values = c(
    "GeneSeekr" = "#E69F00",
    "KMA 70% ConClave 1" = "#66C2A5",
    "StxTyper" = "#8DA0CB",
    "Sipprverse" = "#E76F51"
  )) +
  labs(y = "Estimated Value")

performance_alltools_CI_plot


dpi = 300
jpeg("performance_all_tools_plot_w_CI.jpg", width = 10*dpi, height = 8*dpi, res = dpi)
performance_alltools_CI_plot
dev.off()

