# ... (Sections 1-11 identical to 0.Data_Cleaning.R) ...

# ==================== 13. PCA dimensionality reduction analysis (simplified version) ====================
cat("\n========== PCA dimensionality reduction analysis (simplified version) ==========\n")

# Perform PCA using significant features
pca_result <- prcomp(X_sig, center = TRUE, scale. = TRUE)

# Extract the first two principal components
pc_scores <- as.data.frame(pca_result$x[, 1:2])
colnames(pc_scores) <- c("PC1", "PC2")

# Standardize PC1 and PC2 (Z-score standardization) - correctly extract vectors
pc_scores$Standardized_PC1 <- as.vector(scale(pc_scores$PC1))
pc_scores$Standardized_PC2 <- as.vector(scale(pc_scores$PC2))

# Add sample information
pc_scores$Image <- data$Image
pc_scores$RITH_Phenotype <- data$RITH_Phenotype
pc_scores$Heterogeneity_Index <- data$Heterogeneity_Index

# Reorder columns
pc_scores <- pc_scores[, c("Image", "RITH_Phenotype", "Heterogeneity_Index",
                           "PC1", "PC2", "Standardized_PC1", "Standardized_PC2")]

# Save standardized PC scores
write_csv(pc_scores, "Standardized_PC_Scores.csv")
cat("Standardized PC scores saved to: Standardized_PC_Scores.csv\n")

# View standardized PC values of the first few samples
cat("\nStandardized PC values of the first 5 samples:\n")
print(head(pc_scores[, c("Image", "RITH_Phenotype", "Standardized_PC1", "Standardized_PC2")], 5))

# Calculate mean standardized PC values by phenotype
cat("\nMean standardized PC values by phenotype:\n")
pc_summary <- pc_scores %>%
  group_by(RITH_Phenotype) %>%
  summarise(
    Mean_PC1 = mean(Standardized_PC1, na.rm = TRUE),
    SD_PC1 = sd(Standardized_PC1, na.rm = TRUE),
    Mean_PC2 = mean(Standardized_PC2, na.rm = TRUE),
    SD_PC2 = sd(Standardized_PC2, na.rm = TRUE),
    n = n()
  )
print(pc_summary)

# ==================== 14. Visualize standardized PC1 and PC2 ====================
cat("\n========== Visualize standardized PC1 and PC2 ==========\n")

phenotype_colors <- c("Low RITH" = "#377EB8", "High RITH" = "#E41A1C")

# Figure 1: scatter plot of standardized PC1 and PC2 colored by phenotype
p1 <- ggplot(pc_scores, aes(x = Standardized_PC1, y = Standardized_PC2, color = RITH_Phenotype)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = phenotype_colors) +
  labs(
    title = "Standardized PCA: PC1 vs PC2",
    x = "Standardized PC1 (Z-score)",
    y = "Standardized PC2 (Z-score)",
    color = "RITH Phenotype"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "right",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  # Add zero reference lines on axes
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  # Add ellipse (95% confidence interval)
  stat_ellipse(level = 0.95, linetype = 2, size = 01)

# Figure 2: distribution of HI in standardized PCA space
p2 <- ggplot(pc_scores, aes(x = Standardized_PC1, y = Standardized_PC2, color = Heterogeneity_Index)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_gradient2(
    low = "#377EB8",
    mid = "white",
    high = "#E41A1C",
    midpoint = 0,
    name = "Heterogeneity Index"
  ) +
  labs(
    title = "Heterogeneity Index in Standardized PCA Space",
    x = "Standardized PC1 (Z-score)",
    y = "Standardized PC2 (Z-score)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "right",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", alpha = 0.5)

# Save figures
ggsave("Standardized_PCA_Phenotype.pdf", p1, width = 8, height = 6, dpi = 300)
ggsave("Standardized_PCA_HI.png", p2, width = 8, height = 6, dpi = 300)
ggsave("Standardized_PCA_Phenotype_HI_filtered.pdf", p1, width = 10, height = 8, device = "pdf")
ggsave("Standardized_PCA_HI_filtered.pdf", p2, width = 10, height = 8, device = "pdf")

# Create combined figure
library(gridExtra)
combined_pca <- grid.arrange(p1, p2, ncol = 2)
ggsave("Standardized_PCA_Combined.png", combined_pca, width = 16, height = 6, dpi = 300)

cat("✅ Standardized PCA visualization figures saved\n")

# ==================== 15. Simple statistical analysis ====================
cat("\n========== Statistical analysis of standardized PC values ==========\n")

# Test whether standardized PC1 can distinguish phenotypes
cat("Between-group comparison of standardized PC1:\n")
pc1_low <- pc_scores$Standardized_PC1[pc_scores$RITH_Phenotype == "Low RITH"]
pc1_high <- pc_scores$Standardized_PC1[pc_scores$RITH_Phenotype == "High RITH"]
pc1_test <- wilcox.test(pc1_low, pc1_high)
cat("  Wilcoxon test p value:", format.pval(pc1_test$p.value), "\n")
cat("  Mean in Low RITH group:", round(mean(pc1_low, na.rm = TRUE), 3), "\n")
cat("  Mean in High RITH group:", round(mean(pc1_high, na.rm = TRUE), 3), "\n")

# Test whether standardized PC2 can distinguish phenotypes
cat("\nBetween-group comparison of standardized PC2:\n")
pc2_low <- pc_scores$Standardized_PC2[pc_scores$RITH_Phenotype == "Low RITH"]
pc2_high <- pc_scores$Standardized_PC2[pc_scores$RITH_Phenotype == "High RITH"]
pc2_test <- wilcox.test(pc2_low, pc2_high)
cat("  Wilcoxon test p value:", format.pval(pc2_test$p.value), "\n")
cat("  Mean in Low RITH group:", round(mean(pc2_low, na.rm = TRUE), 3), "\n")
cat("  Mean in High RITH group:", round(mean(pc2_high, na.rm = TRUE), 3), "\n")

# ==================== 16. Save final standardized PC data ====================
cat("\n========== Save final data ==========\n")

# Create data frame containing all important information
final_data <- data.frame(
  Image = pc_scores$Image,
  RITH_Phenotype = pc_scores$RITH_Phenotype,
  Cluster_Label = data$Cluster_Label,
  Heterogeneity_Index = pc_scores$Heterogeneity_Index,
  Standardized_PC1 = pc_scores$Standardized_PC1,
  Standardized_PC2 = pc_scores$Standardized_PC2
)

# Sort by phenotype
final_data <- final_data[order(final_data$RITH_Phenotype), ]

write_csv(final_data, "Final_RITH_Results_with_PC.csv")
cat("Final results (including standardized PCs) saved to: Final_RITH_Results_with_PC.csv\n")

# Display data summary
cat("\nData summary:\n")
cat("Total number of samples:", nrow(final_data), "\n")
cat("Low RITH group:", sum(final_data$RITH_Phenotype == "Low RITH"), "samples\n")
cat("High RITH group:", sum(final_data$RITH_Phenotype == "High RITH"), "samples\n")
cat("Range of standardized PC1:", round(min(final_data$Standardized_PC1), 2), "to",
    round(max(final_data$Standardized_PC1), 2), "\n")
cat("Range of standardized PC2:", round(min(final_data$Standardized_PC2), 2), "to",
    round(max(final_data$Standardized_PC2), 2), "\n")

cat("\n✅ Standardized PCA analysis completed!\n")
