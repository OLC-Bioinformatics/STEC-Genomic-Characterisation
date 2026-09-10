# Secondary KMA performance figure (same style as KMA-only figure)
# Parameter setting: 70% identity cutoff, ConClave mode 1
# Validation levels: nucleotide operon, variant, subtype

# -------------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# -------------------------------------------------------------------------
# Input
# -------------------------------------------------------------------------
# Set your working directory to the folder containing the Excel workbook,
# or replace input_file with the full path to the file.

setwd("/Users/noorshubair/Desktop/chapter2/raw_reads/kma/three_level_validation/")

input_file <- "kmaresults_70_conclave1_secondary_three_level_comparison_19.xlsx"

metrics <- read_excel(input_file, sheet = "Metrics") %>%
  filter(Parameter_Set == "kmaresults_70_conclave1")

# -------------------------------------------------------------------------
# Wilson 95% confidence interval for a binomial proportion
# -------------------------------------------------------------------------
wilson_ci <- function(x, n, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- (z / denom) * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2)))
  c(lower = centre - half, upper = centre + half)
}

# -------------------------------------------------------------------------
# Calculate the five metrics shown in the original KMA-only plot
# -------------------------------------------------------------------------
plot_data <- metrics %>%
  rowwise() %>%
  mutate(
    NPV = TN / (TN + FN),
    Sensitivity_Recall = TP / (TP + FN),
    Precision_PPV = TP / (TP + FP),
    Specificity = TN / (TN + FP),
    Accuracy = (TP + TN) / (TP + TN + FP + FN),

    NPV_low = wilson_ci(TN, TN + FN)["lower"],
    NPV_high = wilson_ci(TN, TN + FN)["upper"],

    Sensitivity_low = wilson_ci(TP, TP + FN)["lower"],
    Sensitivity_high = wilson_ci(TP, TP + FN)["upper"],

    Precision_low = wilson_ci(TP, TP + FP)["lower"],
    Precision_high = wilson_ci(TP, TP + FP)["upper"],

    Specificity_low = wilson_ci(TN, TN + FP)["lower"],
    Specificity_high = wilson_ci(TN, TN + FP)["upper"],

    Accuracy_low = wilson_ci(TP + TN, TP + TN + FP + FN)["lower"],
    Accuracy_high = wilson_ci(TP + TN, TP + TN + FP + FN)["upper"]
  ) %>%
  ungroup()

# Put estimates and CI limits into long format
plot_long <- bind_rows(
  plot_data %>% transmute(Level, Metric = "NPV", Estimate = NPV,
                          Lower = NPV_low, Upper = NPV_high),
  plot_data %>% transmute(Level, Metric = "Sensitivity/Recall", Estimate = Sensitivity_Recall,
                          Lower = Sensitivity_low, Upper = Sensitivity_high),
  plot_data %>% transmute(Level, Metric = "Precision (PPV)", Estimate = Precision_PPV,
                          Lower = Precision_low, Upper = Precision_high),
  plot_data %>% transmute(Level, Metric = "Specificity", Estimate = Specificity,
                          Lower = Specificity_low, Upper = Specificity_high),
  plot_data %>% transmute(Level, Metric = "Accuracy", Estimate = Accuracy,
                          Lower = Accuracy_low, Upper = Accuracy_high)
) %>%
  mutate(
    Level = recode(
      Level,
      Operon = "Nucleotide operon",
      Variant = "Variant",
      Subtype = "Subtype"
    ),
    Level = factor(Level,
                   levels = c("Nucleotide operon", "Variant", "Subtype")),
    Metric = factor(Metric,
                    levels = c("NPV", "Sensitivity/Recall", "Precision (PPV)",
                               "Specificity", "Accuracy"))
  )

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------
# dodge separates the three validation levels within each performance metric
pd <- position_dodge(width = 0.28)

secondary_kma_metrics_plot <- ggplot(
  plot_long,
  aes(x = Metric, y = Estimate, colour = Level, shape = Level)
) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    position = pd,
    width = 0.08,
    linewidth = 0.8
  ) +
  geom_point(
    position = pd,
    size = 3.5,
    stroke = 0.8
  ) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_y_continuous(
    limits = c(0.85, 1.005),
    breaks = seq(0.85, 1.00, by = 0.05),
    minor_breaks = seq(0.875, 0.975, by = 0.05)
  ) +
  labs(
    x = NULL,
    y = "Estimated Value",
    colour = "Validation level",
    shape = "Validation level"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(linewidth = 0.5),
    panel.grid.minor = element_line(linewidth = 0.35),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 13),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  guides(
    colour = guide_legend(title.position = "left", title.hjust = 0.5),
    shape = guide_legend(title.position = "left", title.hjust = 0.5)
  )

# Display
secondary_kma_metrics_plot

# -------------------------------------------------------------------------
# Print the plotted values and CIs for checking
# -------------------------------------------------------------------------
print(
  plot_long %>%
    arrange(Metric, Level) %>%
    mutate(across(c(Estimate, Lower, Upper), ~ round(.x, 6)))
)

# -------------------------------------------------------------------------
# Save figure
# -------------------------------------------------------------------------

dpi = 300
jpeg("/Users/noorshubair/Desktop/R_code_files/chapter2/secondary_KMA70_conclave1_allmetrics_CI.jpg", width = 10*dpi, height = 8*dpi, res = dpi)
secondary_kma_metrics_plot
dev.off()
