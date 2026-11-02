setwd('/mnt/data4/表型研究2026/ccRCC/RHWU/7.病理MITH')

# Clear workspace
rm(list = ls())

# Load required packages
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(tidyr)
library(dplyr)

# Load data
df <- read.csv("Feats_T.csv")

# Extract morphological feature columns (starting with Morph_)
morph_cols <- names(df)[grepl("^Morph_", names(df))]
cat("Found", length(morph_cols), "morphological features\n")

# ==================== Part 1: Create output directory ====================
output_dir <- "Morphological_Analysis_Results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Set figure dimensions (increased to accommodate larger font sizes)
plot_width <- 18
plot_height <- 14
dpi <- 300

# ==================== Part 2: Correlation plots between HI and each morphological feature ====================
cat("\n=== Draw correlation plots between HI and each morphological feature ===\n")

# Create data frame to store correlation results
correlation_results <- data.frame(
  Feature = character(),
  Correlation = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Calculate correlation with HI for each morphological feature and draw scatter plot
cor_plots <- list()

for (i in seq_along(morph_cols)) {
  feature <- morph_cols[i]

  # Calculate correlation
  cor_test <- cor.test(df[[feature]], df$Heterogeneity_Index,
                       use = "complete.obs", method = "pearson")

  # Save results
  correlation_results <- rbind(correlation_results, data.frame(
    Feature = feature,
    Correlation = cor_test$estimate,
    p_value = cor_test$p.value
  ))

  # Create scatter plot - increase font size
  p <- ggplot(df, aes_string(x = "Heterogeneity_Index", y = feature)) +
    geom_point(aes(color = RITH_Phenotype), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", size = 1) +
    theme_minimal(base_size = 14) +
    labs(title = feature,
         subtitle = sprintf("r = %.3f, p = %.3e",
                            cor_test$estimate,
                            cor_test$p.value),
         x = "Heterogeneity Index (HI)",
         y = "Feature Value") +
    scale_color_manual(values = c("High RITH" = "#E41A1C", "Low RITH" = "#377EB8")) +
    theme(legend.position = "none",
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 14))

  cor_plots[[feature]] <- p

  # Show progress every 10 features processed
  if (i %% 10 == 0) {
    cat(sprintf("Processed %d/%d features...\n", i, length(morph_cols)))
  }
}

# Combine and save all correlation scatter plots
cat("Combining all correlation scatter plots...\n")

# Calculate required rows and columns
n_features <- length(morph_cols)
n_cols <- 4
n_rows <- ceiling(n_features / n_cols)

# Create combined figure
all_cor_plots <- do.call(grid.arrange, c(cor_plots, ncol = n_cols))

# Save combined figure
ggsave(file.path(output_dir, "01_HI_vs_Morphological_Features_Correlation.pdf"),
       all_cor_plots,
       width = plot_width,
       height = n_rows * 4,
       dpi = dpi,
       limitsize = FALSE)

cat("Saved: 01_HI_vs_Morphological_Features_Correlation.pdf\n")

# ==================== Part 3: Violin distribution plots for all morphological features (with statistical analysis) ====================
cat("\n=== Draw violin distribution plots for all morphological features (with statistical analysis) ===\n")

# First convert wide data to long data
df_long <- df %>%
  select(id, RITH_Phenotype, all_of(morph_cols)) %>%
  pivot_longer(cols = all_of(morph_cols),
               names_to = "Feature",
               values_to = "Value")

# Create list of violin plots (with statistical analysis)
violin_plots <- list()

# Create violin plot for each feature
unique_features <- unique(df_long$Feature)

for (i in seq_along(unique_features)) {
  feature <- unique_features[i]

  # Extract data for this feature
  feature_data <- df_long[df_long$Feature == feature, ]

  # Check whether there is enough data
  high_data <- feature_data$Value[feature_data$RITH_Phenotype == "High RITH"]
  low_data <- feature_data$Value[feature_data$RITH_Phenotype == "Low RITH"]

  # Create violin plot - with statistical analysis (increase font size)
  p <- ggplot(feature_data, aes(x = RITH_Phenotype, y = Value, fill = RITH_Phenotype)) +
    geom_violin(alpha = 0.7, trim = TRUE, color = "black", size = 0.8) +
    geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA,
                 color = "black", size = 0.8) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1.2, color = "black") +

    # Add statistical comparison
    stat_compare_means(
      aes(label = ..p.signif..),
      method = "t.test",
      comparisons = list(c("High RITH", "Low RITH")),
      method.args = list(var.equal = FALSE),
      label.y = max(feature_data$Value, na.rm = TRUE) * 1.1,
      size = 10,
      bracket.size = 1.2,
      tip.length = 0.02,
      vjust = 0.5
    ) +

    labs(title = feature,
         x = "RITH Phenotype",
         y = "Feature Value") +

    scale_fill_manual(values = c("High RITH" = "#E41A1C", "Low RITH" = "#377EB8")) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      axis.line = element_line(color = "black", size = 0.8),
      axis.ticks = element_line(color = "black", size = 0.8),
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +

    # Adjust y-axis range to accommodate significance asterisks
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

  violin_plots[[feature]] <- p

  # Show progress every 10 features processed
  if (i %% 10 == 0) {
    cat(sprintf("Processed violin plots for %d/%d features...\n", i, length(unique_features)))
  }
}

# Combine and save all violin plots
cat("Combining all violin plots...\n")

# Create combined figure
all_violin_plots <- do.call(grid.arrange, c(violin_plots, ncol = n_cols))

# Save combined figure
ggsave(file.path(output_dir, "02_Morphological_Features_Violin_Plots_with_Stats.pdf"),
       all_violin_plots,
       width = plot_width,
       height = n_rows * 4,
       dpi = dpi,
       limitsize = FALSE)

cat("Saved: 02_Morphological_Features_Violin_Plots_with_Stats.pdf\n")

# ==================== Part 4: Save correlation results and statistical analysis ====================
cat("\n=== Save analysis results ===\n")

# Calculate more detailed statistics for violin plots
violin_stats <- df_long %>%
  group_by(Feature, RITH_Phenotype) %>%
  summarise(
    n = n(),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    q25 = quantile(Value, 0.25, na.rm = TRUE),
    q75 = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate statistical tests for between-group differences
ttest_results <- df_long %>%
  group_by(Feature) %>%
  summarise(
    t_statistic = tryCatch({
      test_result <- t.test(Value ~ RITH_Phenotype, na.action = na.omit)
      test_result$statistic
    }, error = function(e) NA),
    p_value = tryCatch({
      test_result <- t.test(Value ~ RITH_Phenotype, na.action = na.omit)
      test_result$p.value
    }, error = function(e) NA),
    p_signif = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    mean_diff = mean(Value[RITH_Phenotype == "High RITH"], na.rm = TRUE) -
      mean(Value[RITH_Phenotype == "Low RITH"], na.rm = TRUE)
  )

# Merge all results
final_results <- correlation_results %>%
  left_join(ttest_results, by = c("Feature" = "Feature")) %>%
  arrange(p_value.x)

# Rename columns
colnames(final_results) <- c("Feature", "Correlation_with_HI", "Correlation_p_value",
                             "T_statistic", "T_test_p_value", "Significance", "Mean_difference_High_vs_Low")

# Save results
write.csv(final_results,
          file.path(output_dir, "03_Statistical_Analysis_Results.csv"),
          row.names = FALSE)

write.csv(violin_stats,
          file.path(output_dir, "04_Violin_Plot_Statistics.csv"),
          row.names = FALSE)

cat("Saved: 03_Statistical_Analysis_Results.csv\n")
cat("Saved: 04_Violin_Plot_Statistics.csv\n")

# ==================== Part 5: Create individual plots for the most significant features (with more detailed statistical analysis) ====================
cat("\n=== Create individual plots for the most significant features ===\n")

# Identify features most strongly correlated with HI
top_cor_features <- final_results %>%
  arrange(desc(abs(Correlation_with_HI))) %>%
  head(10) %>%
  pull(Feature)

# Identify features with most significant between-group differences
top_diff_features <- final_results %>%
  arrange(T_test_p_value) %>%
  head(10) %>%
  pull(Feature)

# Create scatter plots for features most strongly correlated with HI (increase font size)
cat("Creating scatter plots for features most strongly correlated with HI...\n")
top_cor_plots <- list()

for (i in seq_along(top_cor_features)) {
  feature <- top_cor_features[i]
  cor_value <- final_results$Correlation_with_HI[final_results$Feature == feature]
  p_value_cor <- final_results$Correlation_p_value[final_results$Feature == feature]

  p <- ggplot(df, aes_string(x = "Heterogeneity_Index", y = feature)) +
    geom_point(aes(color = RITH_Phenotype), alpha = 0.7, size = 4) +
    geom_smooth(method = "lm", se = TRUE, color = "black", size = 1.2) +
    theme_minimal(base_size = 16) +
    labs(title = sprintf("%s\nCorrelation: %.3f (p = %.3e)", feature, cor_value, p_value_cor),
         x = "Heterogeneity Index",
         y = "Feature Value") +
    scale_color_manual(values = c("High RITH" = "#E41A1C", "Low RITH" = "#377EB8")) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 16))

  ggsave(file.path(output_dir, sprintf("Top_Correlation_%02d_%s.pdf", i, feature)),
         p, width = 8, height = 6, dpi = dpi)

  top_cor_plots[[feature]] <- p
}

# Combine and save
if (length(top_cor_plots) > 0) {
  top_cor_grid <- do.call(grid.arrange, c(top_cor_plots, ncol = 2))
  ggsave(file.path(output_dir, "05_Top_Correlation_Features.pdf"),
         top_cor_grid, width = 18, height = 7 * ceiling(length(top_cor_plots)/2),
         dpi = dpi, limitsize = FALSE)
  cat("Saved: 05_Top_Correlation_Features.pdf\n")
}

# Create violin plots for features with most significant between-group differences (with detailed statistical analysis, increase font size)
cat("Creating violin plots for features with most significant between-group differences...\n")
top_diff_plots <- list()

for (i in seq_along(top_diff_features)) {
  feature <- top_diff_features[i]
  p_value <- final_results$T_test_p_value[final_results$Feature == feature]
  significance <- final_results$Significance[final_results$Feature == feature]

  # Extract data
  plot_data <- df_long[df_long$Feature == feature, ]

  # Calculate maximum value to determine significance asterisk position
  y_max <- max(plot_data$Value, na.rm = TRUE)

  p <- ggplot(plot_data, aes(x = RITH_Phenotype, y = Value, fill = RITH_Phenotype)) +
    geom_violin(alpha = 0.7, color = "black", size = 0.8) +
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA,
                 color = "black", size = 0.8) +
    geom_jitter(width = 0.05, alpha = 0.5, size = 2, color = "black") +

    # Add significance comparison (with bracket and asterisks)
    stat_compare_means(
      aes(label = paste0(..p.signif.., "\n", "p=", ..p.format..)),
      method = "t.test",
      comparisons = list(c("High RITH", "Low RITH")),
      method.args = list(var.equal = FALSE),
      label.y = y_max * 1.15,
      size = 10,
      bracket.size = 1.2,
      tip.length = 0.03,
      vjust = 0.5
    ) +

    theme_classic(base_size = 16) +
    labs(title = feature,
         x = "RITH Phenotype",
         y = "Feature Value") +
    scale_fill_manual(values = c("High RITH" = "#E41A1C", "Low RITH" = "#377EB8")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 18),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      axis.line = element_line(color = "black", size = 0.8),
      axis.ticks = element_line(color = "black", size = 0.8),
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  ggsave(file.path(output_dir, sprintf("Top_Difference_%02d_%s.pdf", i, feature)),
         p, width = 8, height = 7, dpi = dpi)

  top_diff_plots[[feature]] <- p
}

# Combine and save
if (length(top_diff_plots) > 0) {
  top_diff_grid <- do.call(grid.arrange, c(top_diff_plots, ncol = 2))
  ggsave(file.path(output_dir, "06_Top_Difference_Features_with_Stats.pdf"),
         top_diff_grid, width = 18, height = 8 * ceiling(length(top_diff_plots)/2),
         dpi = dpi, limitsize = FALSE)
  cat("Saved: 06_Top_Difference_Features_with_Stats.pdf\n")
}

# ==================== Part 6: Create significance heatmap (optional) ====================
cat("\n=== Create significance heatmap ===\n")

# Prepare heatmap data
heatmap_data <- final_results %>%
  select(Feature, Correlation_with_HI, T_test_p_value) %>%
  mutate(
    Correlation_group = cut(abs(Correlation_with_HI),
                            breaks = c(0, 0.1, 0.3, 0.5, 1),
                            labels = c("Very Weak", "Weak", "Moderate", "Strong")),
    Pvalue_group = cut(T_test_p_value,
                       breaks = c(0, 0.001, 0.01, 0.05, 1),
                       labels = c("***", "**", "*", "ns"))
  )

# Save heatmap data
write.csv(heatmap_data,
          file.path(output_dir, "07_Significance_Heatmap_Data.csv"),
          row.names = FALSE)

cat("Saved: 07_Significance_Heatmap_Data.csv\n")

# ==================== Part 7: Generate analysis report ====================
cat("\n=== Generate analysis report ===\n")

sink(file.path(output_dir, "08_Analysis_Report.txt"))

cat("========================================\n")
cat("    Morphological feature analysis report (with statistical analysis)\n")
cat("========================================\n\n")

cat("Analysis date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("Analysis time:", format(Sys.time(), "%H:%M:%S"), "\n\n")

cat("=== Dataset information ===\n")
cat("Total number of samples:", nrow(df), "\n")
cat("Number of High RITH samples:", sum(df$RITH_Phenotype == "High RITH"), "\n")
cat("Number of Low RITH samples:", sum(df$RITH_Phenotype == "Low RITH"), "\n")
cat("Number of morphological features analyzed:", length(morph_cols), "\n\n")

cat("=== Explanation of statistical significance symbols ===\n")
cat("ns: p >= 0.05\n")
cat("*: p < 0.05\n")
cat("**: p < 0.01\n")
cat("***: p < 0.001\n\n")

cat("=== Most correlated features (with HI) ===\n")
cat("Top 10 features most correlated with Heterogeneity Index:\n")
print(head(final_results[order(abs(final_results$Correlation_with_HI), decreasing = TRUE),
                         c("Feature", "Correlation_with_HI", "Correlation_p_value", "Significance")], 10))

cat("\n=== Features with most significant between-group differences ===\n")
cat("Top 10 features with most significant differences between High and Low RITH groups:\n")
print(head(final_results[order(final_results$T_test_p_value),
                         c("Feature", "Mean_difference_High_vs_Low", "T_test_p_value", "Significance")], 10))

cat("\n=== Features that are significantly correlated and significantly different ===\n")
significant_both <- final_results %>%
  filter(Correlation_p_value < 0.05 & T_test_p_value < 0.05) %>%
  arrange(desc(abs(Correlation_with_HI)))

if (nrow(significant_both) > 0) {
  cat("Features significantly correlated with HI and significantly different between groups:\n")
  print(significant_both[, c("Feature", "Correlation_with_HI", "Correlation_p_value",
                             "T_test_p_value", "Significance")])
} else {
  cat("No features satisfying both conditions were found\n")
}

cat("\n=== List of generated files ===\n")
cat("01_HI_vs_Morphological_Features_Correlation.pdf - correlation scatter plots between HI and all morphological features\n")
cat("02_Morphological_Features_Violin_Plots_with_Stats.pdf - violin distribution plots for all morphological features (with statistical analysis)\n")
cat("03_Statistical_Analysis_Results.csv - statistical analysis results\n")
cat("04_Violin_Plot_Statistics.csv - detailed statistics for violin plots\n")
cat("05_Top_Correlation_Features.pdf - features most correlated with HI\n")
cat("06_Top_Difference_Features_with_Stats.pdf - features with most significant between-group differences (with statistical analysis)\n")
cat("07_Significance_Heatmap_Data.csv - significance heatmap data\n")
cat("08_Analysis_Report.txt - this analysis report\n")

cat("\n=== Statistical analysis notes ===\n")
cat("1. Correlation analysis: Pearson correlation coefficient\n")
cat("2. Between-group difference: two-sample t-test (unequal variance)\n")
cat("3. Significance annotation: standard significance symbols (*, **, ***, ns)\n")
cat("4. Comparison line in violin plots: shows significance of High vs Low RITH\n")
cat("5. All tests are two-sided\n")

cat("\n========================================\n")
cat("          Analysis completed\n")
cat("========================================\n")

sink()

cat("Saved: 08_Analysis_Report.txt\n")

# ==================== Part 8: Display completion message ====================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Analysis completed! All results saved to the '", output_dir, "' folder\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
