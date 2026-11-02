# ... (Sections 1-11 identical to 0.Data_Cleaning.R) ...

# ==================== 17. Violin plot and boxplot of HI ====================
cat("\n========== Visualization of HI violin plot ==========\n")

# Use specified colors: blue = Low RITH, red = High RITH
phenotype_colors <- c("Low RITH" = "#377EB8",   # Blue
                      "High RITH" = "#E41A1C")  # Red

# Figure 1: HI violin plot - increase font size
p_hi_violin <- ggplot(data, aes(x = RITH_Phenotype, y = Heterogeneity_Index, fill = RITH_Phenotype)) +
  geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = RITH_Phenotype), width = 0.01, size = 1, alpha = 0.5) +
  scale_fill_manual(values = phenotype_colors) +
  scale_color_manual(values = phenotype_colors) +
  labs(
    title = "",
    x = "RITH Phenotype",
    y = "Heterogeneity Index (HI)",
    fill = "Phenotype"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "none",
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# Manually add significance annotation
y_max <- max(data$Heterogeneity_Index, na.rm = TRUE)
y_min <- min(data$Heterogeneity_Index, na.rm = TRUE)
y_range <- y_max - y_min

signif_label <- ifelse(hi_p_value < 0.001, "***",
                       ifelse(hi_p_value < 0.01, "**",
                              ifelse(hi_p_value < 0.05, "*", "ns")))

# Add significance annotation (adjust size)
p_hi_violin <- p_hi_violin +
  geom_segment(x = 1, xend = 2,
               y = y_max + 0.1 * y_range,
               yend = y_max + 0.1 * y_range,
               color = "black", linewidth = 1) +
  geom_segment(x = 1, xend = 1,
               y = y_max + 0.08 * y_range,
               yend = y_max + 0.1 * y_range,
               color = "black", linewidth = 1) +
  geom_segment(x = 2, xend = 2,
               y = y_max + 0.08 * y_range,
               yend = y_max + 0.1 * y_range,
               color = "black", linewidth = 1) +
  annotate("text", x = 1.5, y = y_max + 0.13 * y_range,
           label = signif_label,
           size = 8, fontface = "bold", color = "black")

# Figure 2: HI boxplot - increase font size
p_hi_box <- ggplot(data, aes(x = RITH_Phenotype, y = Heterogeneity_Index, fill = RITH_Phenotype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 19, outlier.size = 3, outlier.alpha = 0.5) +
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.4, color = "black") +
  scale_fill_manual(values = phenotype_colors) +
  labs(
    title = "Heterogeneity Index (HI) Boxplot",
    x = "RITH Phenotype",
    y = "Heterogeneity Index (HI)",
    fill = "Phenotype"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "none",
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(..y.., 1)),
    vjust = -1.5,
    size = 5,
    fontface = "bold"
  ) +
  geom_segment(x = 1, xend = 2,
               y = y_max + 0.1 * y_range,
               yend = y_max + 0.1 * y_range,
               color = "black", linewidth = 1) +
  geom_segment(x = 1, xend = 1,
               y = y_max + 0.08 * y_range,
               yend = y_max + 0.1 * y_range,
               color = "black", linewidth = 1) +
  geom_segment(x = 2, xend = 2,
               y = y_max + 0.08 * y_range,
               yend = y_max + 0.1 * y_range,
               color = "black", linewidth = 1) +
  annotate("text", x = 1.5, y = y_max + 0.13 * y_range,
           label = signif_label,
           size = 8, fontface = "bold", color = "black")

# Figure 3: HI density plot - increase font size
p_hi_density <- ggplot(data, aes(x = Heterogeneity_Index, fill = RITH_Phenotype)) +
  geom_density(alpha = 0.6, color = NA) +
  scale_fill_manual(values = phenotype_colors) +
  labs(
    title = "Density Distribution of Heterogeneity Index",
    x = "Heterogeneity Index (HI)",
    y = "Density",
    fill = "RITH Phenotype"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  ) +
  geom_vline(
    data = data %>% group_by(RITH_Phenotype) %>% summarise(median_HI = median(Heterogeneity_Index)),
    aes(xintercept = median_HI, color = RITH_Phenotype),
    linetype = "dashed",
    linewidth = 1.2,
    show.legend = FALSE
  ) +
  scale_color_manual(values = phenotype_colors)

# Figure 4: HI histogram (facetted) - increase font size
p_hi_histogram <- ggplot(data, aes(x = Heterogeneity_Index, fill = RITH_Phenotype)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "white", linewidth = 0.3) +
  facet_wrap(~ RITH_Phenotype, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = phenotype_colors) +
  labs(
    title = "Histogram of Heterogeneity Index by Phenotype",
    x = "Heterogeneity Index (HI)",
    y = "Count",
    fill = "RITH Phenotype"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "none",
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# Save individual figures (PDF, high quality)
ggsave("HI_Violin_Plot.pdf", p_hi_violin, width = 8, height = 7, dpi = 300)
ggsave("HI_Box_Plot.pdf", p_hi_box, width = 8, height = 7, dpi = 300)
ggsave("HI_Density_Plot.pdf", p_hi_density, width = 10, height = 7, dpi = 300)
ggsave("HI_Histogram.pdf", p_hi_histogram, width = 8, height = 10, dpi = 300)

# Also save PNG versions for quick viewing
ggsave("HI_Violin_Plot.png", p_hi_violin, width = 8, height = 7, dpi = 300)
ggsave("HI_Box_Plot.png", p_hi_box, width = 8, height = 7, dpi = 300)
ggsave("HI_Density_Plot.png", p_hi_density, width = 10, height = 7, dpi = 300)
ggsave("HI_Histogram.png", p_hi_histogram, width = 8, height = 10, dpi = 300)

cat("✅ HI visualization figures saved\n")

# ==================== 18. HI statistical summary ====================
cat("\n========== HI statistical summary ==========\n")

# Calculate detailed statistical summary
hi_stats <- data %>%
  group_by(RITH_Phenotype) %>%
  summarise(
    n = n(),
    Mean = round(mean(Heterogeneity_Index, na.rm = TRUE), 2),
    SD = round(sd(Heterogeneity_Index, na.rm = TRUE), 2),
    Median = round(median(Heterogeneity_Index, na.rm = TRUE), 2),
    IQR = round(IQR(Heterogeneity_Index, na.rm = TRUE), 2),
    Min = round(min(Heterogeneity_Index, na.rm = TRUE), 2),
    Max = round(max(Heterogeneity_Index, na.rm = TRUE), 2),
    Q1 = round(quantile(Heterogeneity_Index, 0.25, na.rm = TRUE), 2),
    Q3 = round(quantile(Heterogeneity_Index, 0.75, na.rm = TRUE), 2)
  )

cat("Detailed statistical summary of HI:\n")
print(hi_stats)

# Save statistical summary
write_csv(hi_stats, "HI_Statistics_Summary.csv")
cat("HI statistical summary saved to: HI_Statistics_Summary.csv\n")

# Calculate effect size (Cohen's d)
calculate_cohens_d <- function(x, y) {
  mean_diff <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
  n1 <- length(x)
  n2 <- length(y)
  pooled_sd <- sqrt(((n1-1)*var(x, na.rm = TRUE) + (n2-1)*var(y, na.rm = TRUE)) / (n1 + n2 - 2))
  return(mean_diff / pooled_sd)
}

hi_low <- data$Heterogeneity_Index[data$RITH_Phenotype == "Low RITH"]
hi_high <- data$Heterogeneity_Index[data$RITH_Phenotype == "High RITH"]
cohens_d <- calculate_cohens_d(hi_high, hi_low)

cat("\nEffect size analysis:\n")
cat("Cohen's d =", round(cohens_d, 3), "\n")
if (abs(cohens_d) < 0.2) {
  cat("Effect size: very small\n")
} else if (abs(cohens_d) < 0.5) {
  cat("Effect size: small\n")
} else if (abs(cohens_d) < 0.8) {
  cat("Effect size: medium\n")
} else {
  cat("Effect size: large\n")
}

cat("\n✅ HI visualization analysis completed!\n")
