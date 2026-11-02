library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(readr)

setwd('/mnt/data4/表型研究2026/ccRCC/RHWU/6.HI和临床变量的相关性')
# ======================
# 1. Load required packages
# ======================
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(readr)

# ======================
# 2. Read data
# ======================
data <- read.csv("HI.csv", stringsAsFactors = FALSE)

# ======================
# 3. Data cleaning and preparation
# ======================

# Extract numeric pT stage (T1->1, T2->2, T3->3, T4->4)
data$pT_numeric <- as.numeric(gsub("T", "", data$pT_stage))
# Extract numeric ISUP grade (G1->1, G2->2, G3->3, G4->4, G5->5)
data$ISUP_numeric <- as.numeric(gsub("G", "", data$ISUP_grade))

# ======================
# 4. Correlation analysis
# ======================
# Correlation between HI and pT stage (Spearman)
cor_pT <- cor.test(data$Heterogeneity_Index, data$pT_numeric, method = "spearman")
cat("Spearman correlation between HI and pT stage:\n")
cat("rho =", round(cor_pT$estimate, 3), "\n")
cat("p =", format.pval(cor_pT$p.value, digits = 3), "\n\n")

# Correlation between HI and ISUP grade (Spearman)
cor_ISUP <- cor.test(data$Heterogeneity_Index, data$ISUP_numeric, method = "spearman")
cat("Spearman correlation between HI and ISUP grade:\n")
cat("rho =", round(cor_ISUP$estimate, 3), "\n")
cat("p =", format.pval(cor_ISUP$p.value, digits = 3), "\n\n")

# ======================
# 5. Draw correlation plots
# ======================
# Figure 1: scatter plot of HI vs pT stage
p1 <- ggplot(data, aes(x = pT_numeric, y = Heterogeneity_Index)) +
  geom_point(aes(color = pT_stage), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE, linetype = "dashed") +
  labs(title = paste("Heterogeneity Index vs pT Stage\nSpearman rho =",
                     round(cor_pT$estimate, 3),
                     ", p =", format.pval(cor_pT$p.value, digits = 3)),
       x = "pT Stage (numeric)",
       y = "Heterogeneity Index (HI)",
       color = "pT Stage") +
  scale_x_continuous(breaks = 1:4, labels = c("T1", "T2", "T3", "T4")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

# Figure 2: scatter plot of HI vs ISUP grade
p2 <- ggplot(data, aes(x = ISUP_numeric, y = Heterogeneity_Index)) +
  geom_point(aes(color = ISUP_grade), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "blue", se = TRUE, linetype = "dashed") +
  labs(title = paste("Heterogeneity Index vs ISUP Grade\nSpearman rho =",
                     round(cor_ISUP$estimate, 3),
                     ", p =", format.pval(cor_ISUP$p.value, digits = 3)),
       x = "ISUP Grade (numeric)",
       y = "Heterogeneity Index (HI)",
       color = "ISUP Grade") +
  scale_x_continuous(breaks = 2:5, labels = c("G2", "G3", "G4", "G5")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

# Figure 3: grouped boxplot - HI by pT stage
p3 <- ggplot(data, aes(x = pT_stage, y = Heterogeneity_Index, fill = pT_stage)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  labs(title = "Heterogeneity Index Distribution by pT Stage",
       x = "pT Stage",
       y = "Heterogeneity Index (HI)",
       fill = "pT Stage") +
  stat_compare_means(method = "kruskal.test",
                     label = "p.format",
                     label.y = max(data$Heterogeneity_Index) * 1.1) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

# Figure 4: grouped boxplot - HI by ISUP grade
p4 <- ggplot(data, aes(x = ISUP_grade, y = Heterogeneity_Index, fill = ISUP_grade)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  labs(title = "Heterogeneity Index Distribution by ISUP Grade",
       x = "ISUP Grade",
       y = "Heterogeneity Index (HI)",
       fill = "ISUP Grade") +
  stat_compare_means(method = "kruskal.test",
                     label = "p.format",
                     label.y = max(data$Heterogeneity_Index) * 1.1) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

# ======================
# 6. Combine plots and save
# ======================
# Combine all plots
combined_plot <- ggarrange(p1, p2, p3, p4,
                           ncol = 2, nrow = 2,
                           labels = c("A", "B", "C", "D"))

# Display figure
print(combined_plot)

# Save figure as PDF
ggsave("HI_correlation_plots.pdf", combined_plot,
       width = 14, height = 10, dpi = 300)

# Save figure as PNG
ggsave("HI_correlation_plots.png", combined_plot,
       width = 14, height = 10, dpi = 300)

cat("\nFigures saved as HI_correlation_plots.pdf and HI_correlation_plots.png\n")

# ======================
# 7. Output statistical results
# ======================
# Create results data frame
results_df <- data.frame(
  Comparison = c("HI vs pT Stage", "HI vs ISUP Grade"),
  Spearman_rho = c(round(cor_pT$estimate, 3), round(cor_ISUP$estimate, 3)),
  P_value = c(format.pval(cor_pT$p.value, digits = 3),
              format.pval(cor_ISUP$p.value, digits = 3)),
  Interpretation = c(
    ifelse(cor_pT$p.value < 0.05, "Significant", "Not significant"),
    ifelse(cor_ISUP$p.value < 0.05, "Significant", "Not significant")
  )
)

print(results_df)

# Save results as CSV
write.csv(results_df, "HI_correlation_results.csv", row.names = FALSE)
cat("\nStatistical results saved as HI_correlation_results.csv\n")
