# ==================== Load required R packages ====================
library(readr)
library(dplyr)
library(SNFtool)
library(Spectrum)
library(ggplot2)
library(pheatmap)
library(factoextra)
library(cluster)
library(NbClust)
library(tidyr)
library(stringr)

setwd('/mnt/data4/表型研究2026/ccRCC')
# ==================== 1. Read data ====================
data <- read_csv("/mnt/data4/表型研究2026/ccRCC/RITH_Data_Cleaned.csv")
head(data)

# ==================== 2. Extract all feature names ====================
all_features <- colnames(data)
all_features <- all_features[all_features != "Image"]

cat("Total number of features:", length(all_features), "\n")
cat("Number of samples:", nrow(data), "\n")

# ==================== 3. Feature classification ====================
cat("\n========== Feature classification ==========\n")

identify_filter_type <- function(feature_name) {
  if (grepl("^original_", feature_name)) {
    return("original")
  } else if (grepl("^log-sigma", feature_name)) {
    sigma_match <- str_extract(feature_name, "log-sigma-([0-9.]+)-mm-3D")
    if (!is.na(sigma_match)) {
      sigma_value <- str_extract(sigma_match, "[0-9.]+")
      return(paste0("log-sigma-", sigma_value))
    }
    return("log-sigma")
  } else if (grepl("^wavelet-", feature_name)) {
    wavelet_match <- str_extract(feature_name, "wavelet-[A-Z]{3}")
    if (!is.na(wavelet_match)) {
      return(wavelet_match)
    }
    return("wavelet")
  }
  return("unknown")
}

identify_feature_class <- function(feature_name) {
  if (grepl("_firstorder_", feature_name)) {
    return("firstorder")
  } else if (grepl("_glcm_", feature_name)) {
    return("glcm")
  } else if (grepl("_glrlm_", feature_name)) {
    return("glrlm")
  } else if (grepl("_glszm_", feature_name)) {
    return("glszm")
  } else if (grepl("_gldm_", feature_name)) {
    return("gldm")
  } else if (grepl("_ngtdm_", feature_name)) {
    return("ngtdm")
  }
  return("unknown")
}

identify_feature_name <- function(feature_name) {
  parts <- str_split(feature_name, "_", simplify = TRUE)
  if (length(parts) >= 3) {
    return(parts[length(parts)])
  }
  return(feature_name)
}

feature_info <- data.frame(
  Feature = all_features,
  Filter_Type = sapply(all_features, identify_filter_type),
  Feature_Class = sapply(all_features, identify_feature_class),
  Feature_Name = sapply(all_features, identify_feature_name)
)

cat("Distribution of filter types:\n")
print(table(feature_info$Filter_Type))
cat("\nDistribution of feature classes:\n")
print(table(feature_info$Feature_Class))

# ==================== 4. MAD-based feature selection ====================
cat("\n========== MAD-based feature selection ==========\n")

X <- data[, all_features] %>% as.matrix()
rownames(X) <- data$Image

mad_values <- apply(X, 2, mad, na.rm = TRUE)
significant_features <- names(mad_values[mad_values > 1])

cat("Total number of features:", length(all_features), "\n")
cat("Number of significant features with MAD > 1:", length(significant_features), "\n")

significant_info <- feature_info[feature_info$Feature %in% significant_features, ]
cat("\nDistribution of significant features by filter type:\n")
print(table(significant_info$Filter_Type))
cat("\nDistribution of significant features by feature class:\n")
print(table(significant_info$Feature_Class))

X_sig <- X[, significant_features]

# ==================== 5. Determine feature direction based on positive correlation with ITH ====================
cat("\n========== Determine feature direction based on positive correlation with ITH ==========\n")

# Principle: most features are positively correlated with ITH.
# Therefore, samples with higher feature values should have higher heterogeneity.
# We first calculate a "heterogeneity tendency" score and then determine feature direction
# based on its correlation with this score.

# Step 1: calculate the "heterogeneity tendency score" for each sample.
# Since most features are positively correlated with ITH, we can simply use the mean feature value
# as the initial heterogeneity tendency.
feature_means <- rowMeans(X_sig, na.rm = TRUE)

# Step 2: calculate the correlation between each feature and the heterogeneity tendency score.
z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

X_z <- apply(X_sig, 2, z_score)
rownames(X_z) <- data$Image

# Calculate the Spearman correlation between each feature and the heterogeneity tendency.
correlations <- apply(X_z, 2, function(feat) {
  cor(feat, feature_means, method = "spearman", use = "complete.obs")
})

# Step 3: determine feature direction based on correlation.
# Positive correlation (>0): HHR feature (higher value indicates higher heterogeneity).
# Negative correlation (<0): LHR feature (lower value indicates higher heterogeneity).
# Note: according to your description "most features are positively correlated with ITH",
# the threshold is set to 0.
feature_direction <- ifelse(correlations > 0, "HHR_feature", "LHR_feature")
names(feature_direction) <- significant_features

# Summary statistics
HHR_features <- names(feature_direction[feature_direction == "HHR_feature"])
LHR_features <- names(feature_direction[feature_direction == "LHR_feature"])

cat("\nDistribution of feature directions based on correlation with heterogeneity tendency:\n")
print(table(feature_direction))

cat("\nHHR features (positively correlated with heterogeneity):", length(HHR_features), "features (",
    round(length(HHR_features)/length(significant_features)*100, 1), "%)\n")
cat("LHR features (negatively correlated with heterogeneity):", length(LHR_features), "features (",
    round(length(LHR_features)/length(significant_features)*100, 1), "%)\n")

# Save correlation analysis results
cor_df <- data.frame(
  Feature = significant_features,
  Feature_Name = sapply(significant_features, identify_feature_name),
  Feature_Class = sapply(significant_features, identify_feature_class),
  Filter_Type = sapply(significant_features, identify_filter_type),
  Correlation_with_Heterogeneity = correlations,
  Direction = feature_direction
)

# Sort by absolute correlation
cor_df$AbsCorrelation <- abs(cor_df$Correlation_with_Heterogeneity)
cor_df <- cor_df[order(-cor_df$AbsCorrelation), ]

write_csv(cor_df, "Feature_Correlation_with_Heterogeneity.csv")
cat("\nFeature-heterogeneity correlation analysis saved to: Feature_Correlation_with_Heterogeneity.csv\n")

cat("\nTop 10 features most strongly correlated with heterogeneity:\n")
print(head(cor_df, 10))

# ==================== 6. Save feature direction definitions ====================
cat("\n========== Save feature direction definitions ==========\n")

direction_definitions <- data.frame(
  Feature = significant_features,
  Filter_Type = sapply(significant_features, identify_filter_type),
  Feature_Class = sapply(significant_features, identify_feature_class),
  Feature_Name = sapply(significant_features, identify_feature_name),
  Correlation = correlations,
  Direction = feature_direction,
  Description = ifelse(feature_direction == "HHR_feature",
                       "Positively correlated with heterogeneity (higher value = higher heterogeneity)",
                       "Negatively correlated with heterogeneity (lower value = higher heterogeneity)")
)

write_csv(direction_definitions, "Feature_Direction_Definitions_CorrelationBased.csv")
cat("Correlation-based feature direction definitions saved to: Feature_Direction_Definitions_CorrelationBased.csv\n")

# ==================== 7. Perform SNF clustering analysis (2 clusters) ====================
cat("\n========== Perform SNF clustering analysis (2 clusters) ==========\n")

set.seed(123)
optimal_k <- 2
cat("Using number of clusters K =", optimal_k, "\n")

# Calculate similarity matrix
W1 <- affinityMatrix(dist2(X_sig, X_sig), K = 20, sigma = 0.5)

# Perform SNF clustering
cluster_labels <- spectralClustering(W1, K = optimal_k)

# Note: at this point we do not know which label corresponds to high ITH and which to low ITH.
data$Cluster_Label <- cluster_labels

cat("Initial cluster label distribution:\n")
print(table(cluster_labels))

# ==================== 8. Determine which label is the high-RITH subgroup ====================
cat("\n========== Determine high-RITH subgroup (red) and low-RITH subgroup (blue) ==========\n")

# Principle: based on "most features are positively correlated with ITH",
# samples with higher feature values should be high RITH.
# Calculate mean feature value for each cluster label.
cluster_means <- tapply(feature_means, cluster_labels, mean, na.rm = TRUE)
cat("\nMean feature value of each cluster label:\n")
print(cluster_means)

# The cluster with higher feature values should be the high-RITH subgroup (red).
if (cluster_means[1] > cluster_means[2]) {
  high_ith_label <- 1
  low_ith_label <- 2
} else {
  high_ith_label <- 2
  low_ith_label <- 1
}

cat("\nLabel", high_ith_label, "determined as high-RITH subgroup (red)\n")
cat("Label", low_ith_label, "determined as low-RITH subgroup (blue)\n")

# Create final phenotype labels
data$RITH_Phenotype <- ifelse(cluster_labels == high_ith_label, "High RITH", "Low RITH")
data$RITH_Phenotype <- factor(data$RITH_Phenotype, levels = c("Low RITH", "High RITH"))

cat("\nFinal phenotype distribution:\n")
print(table(data$RITH_Phenotype))

# Validate the phenotype definition: mean feature value in High RITH should be significantly higher than in Low RITH.
cat("\nValidate the phenotype definition:\n")
cat("Mean feature value in Low RITH group:", mean(feature_means[data$RITH_Phenotype == "Low RITH"], na.rm = TRUE), "\n")
cat("Mean feature value in High RITH group:", mean(feature_means[data$RITH_Phenotype == "High RITH"], na.rm = TRUE), "\n")

p_value <- wilcox.test(feature_means[data$RITH_Phenotype == "Low RITH"],
                       feature_means[data$RITH_Phenotype == "High RITH"])$p.value
cat("P value for difference in mean feature values between groups:", format.pval(p_value), "\n")

if (p_value < 0.05) {
  cat("✅ Phenotype definition is correct: feature values are significantly higher in High RITH than in Low RITH (p < 0.05)\n")
} else {
  cat("⚠️ Warning: the difference in mean feature values is not significant; rechecking may be needed\n")
}

# ==================== 9. Calculate heterogeneity index (HI) ====================
cat("\n========== Calculate heterogeneity index (HI) ==========\n")

# HI = sum of HHR features - sum of LHR features
HI <- rowSums(X_z[, HHR_features], na.rm = TRUE) -
  rowSums(X_z[, LHR_features], na.rm = TRUE)

data$Heterogeneity_Index <- HI

# Validate the association between HI and phenotype
cat("\nValidate the association between HI and phenotype:\n")
cat("Mean HI in Low RITH group:", mean(HI[data$RITH_Phenotype == "Low RITH"], na.rm = TRUE), "\n")
cat("Mean HI in High RITH group:", mean(HI[data$RITH_Phenotype == "High RITH"], na.rm = TRUE), "\n")

hi_p_value <- wilcox.test(HI[data$RITH_Phenotype == "Low RITH"],
                          HI[data$RITH_Phenotype == "High RITH"])$p.value
cat("P value for HI difference between groups:", format.pval(hi_p_value), "\n")

if (hi_p_value < 0.05) {
  cat("✅ HI significantly distinguishes High/Low RITH phenotypes (p < 0.05)\n")
} else {
  cat("⚠️ HI did not significantly distinguish High/Low RITH phenotypes\n")
}

# ==================== 10. Save phenotype and HI results ====================
cat("\n========== Save phenotype and HI results ==========\n")

results <- data.frame(
  Image = data$Image,
  Cluster_Label = data$Cluster_Label,
  RITH_Phenotype = data$RITH_Phenotype,
  Feature_Mean = feature_means,
  Heterogeneity_Index = data$Heterogeneity_Index
)

write_csv(results, "RITH_Phenotype_and_HI_Final.csv")
cat("Phenotype and HI results saved to: RITH_Phenotype_and_HI_Final.csv\n")

# ==================== 11. Optimize heatmap: separate vertically by direction ====================
cat("\n========== Draw heatmap separated vertically by direction ==========\n")

# Prepare heatmap data (Z-score standardized data)
heatmap_data <- X_z

# Sort samples by phenotype (Low RITH first, High RITH last)
phenotype_order <- order(data$RITH_Phenotype)
heatmap_data <- heatmap_data[phenotype_order, ]

# Ensure the features used exist in the data
HHR_features <- intersect(HHR_features, colnames(heatmap_data))
LHR_features <- intersect(LHR_features, colnames(heatmap_data))

cat("Number of HHR features used in heatmap:", length(HHR_features), "\n")
cat("Number of LHR features used in heatmap:", length(LHR_features), "\n")

# Reorder data: LHR features first, then HHR features
heatmap_data_ordered <- heatmap_data[, c(LHR_features, HHR_features)]

# Process HI values: truncate extreme values to between -200 and 200
HI_truncated <- data$Heterogeneity_Index
HI_truncated[HI_truncated > 200] <- 200
HI_truncated[HI_truncated < -200] <- -100

# Calculate split position
split_position <- length(LHR_features)

# Data processing
heatmap_data_scaled <- scale(heatmap_data_ordered)

heatmap_data_scaled <- t(heatmap_data_scaled)  # Transpose: rows are features, columns are samples

# Truncate extreme values
heatmap_data_scaled[heatmap_data_scaled > 2] <- 2
heatmap_data_scaled[heatmap_data_scaled < -2] <- -2

# Create sample annotation
annotation_row <- data.frame(
  RITH_Phenotype = factor(data$RITH_Phenotype[phenotype_order], levels = c("Low RITH", "High RITH")),
  Heterogeneity_Index = HI_truncated[phenotype_order]
)
rownames(annotation_row) <- rownames(heatmap_data_ordered)

# Create feature annotation
feature_direction_vector <- c(
  rep("LHR_feature", length(LHR_features)),
  rep("HHR_feature", length(HHR_features))
)

annotation_col <- data.frame(
  Direction = factor(feature_direction_vector, levels = c("LHR_feature", "HHR_feature")),
  row.names = rownames(heatmap_data_scaled)
)

# Color settings (as described: red = High RITH, blue = Low RITH)
phenotype_colors <- c("Low RITH" = "#377EB8",   # Blue
                      "High RITH" = "#E41A1C")  # Red

direction_colors <- c(
  "LHR_feature" = "#377EB8",  # Blue: lower value indicates higher heterogeneity
  "HHR_feature" = "#E41A1C"   # Red: higher value indicates higher heterogeneity
)

hi_colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(50)
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

heatmap_data_scaled <- scale(heatmap_data_scaled)
heatmap_data_scaled <- scale(heatmap_data_scaled)
# Draw heatmap
pheatmap(
  heatmap_data_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_col = annotation_row,
  annotation_row = annotation_col,
  annotation_colors = list(
    RITH_Phenotype = phenotype_colors,
    Direction = direction_colors,
    Heterogeneity_Index = hi_colors
  ),
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 50),
  scale = "none",
  fontsize_row = 8,
  fontsize_col = 8,
  gaps_row = split_position,
  annotation_legend = TRUE,
  legend = TRUE,
  border_color = NA,
  filename = "RITH_Direction_Separated_Heatmap.png",
  width = 10,
  height = 10
)

cat("✅ Direction-separated heatmap saved as: RITH_Direction_Separated_Heatmap.png\n")

# ==================== 12. Summary report ====================
cat("\n========== Analysis summary report ==========\n")
cat("1. Total number of samples:", nrow(data), "\n")
cat("2. Number of significant features (MAD > 1):", length(significant_features), "\n")
cat("3. Number of positively correlated features (HHR):", length(HHR_features),
    "(", round(length(HHR_features)/length(significant_features)*100, 1), "%)\n")
cat("4. Number of negatively correlated features (LHR):", length(LHR_features),
    "(", round(length(LHR_features)/length(significant_features)*100, 1), "%)\n")
cat("5. Low-RITH subgroup (blue):", sum(data$RITH_Phenotype == "Low RITH"), "samples\n")
cat("6. High-RITH subgroup (red):", sum(data$RITH_Phenotype == "High RITH"), "samples\n")
cat("7. Phenotype validation p value:", format.pval(p_value), "\n")
cat("8. HI discrimination p value:", format.pval(hi_p_value), "\n")
cat("\n✅ Analysis completed!\n")

# ==================== 13. Remove samples with HI < 0 in the High RITH group ====================
cat("\n========== Remove samples with HI < 0 in the High RITH group ==========\n")

# Identify samples with HI < 0 in the High RITH group
high_rith_samples <- data$Image[data$RITH_Phenotype == "High RITH"]
low_hi_samples <- data$Image[data$RITH_Phenotype == "High RITH" & data$Heterogeneity_Index < 0]

cat("Total number of samples in High RITH group:", length(high_rith_samples), "\n")
cat("Number of samples with HI < 0 in High RITH group:", length(low_hi_samples), "\n")
cat("Proportion of samples with HI < 0 in High RITH group:", round(length(low_hi_samples)/length(high_rith_samples)*100, 1), "%\n")

# Show detailed information for these samples
if (length(low_hi_samples) > 0) {
  cat("\nSamples with HI < 0 in High RITH group:\n")
  for (sample_id in low_hi_samples) {
    hi_value <- data$Heterogeneity_Index[data$Image == sample_id]
    cat(sample_id, sprintf(" (HI = %.2f)", hi_value), "\n")
  }
}

# Remove these samples
data_cleaned <- data[!(data$Image %in% low_hi_samples), ]

cat("\nTotal samples before removal:", nrow(data), "\n")
cat("Total samples after removal:", nrow(data_cleaned), "\n")
cat("Number of removed samples:", nrow(data) - nrow(data_cleaned), "\n")

# Recount phenotype distribution
cat("\nNew phenotype distribution:\n")
print(table(data_cleaned$RITH_Phenotype))

# Validate data after removal
cat("\nValidate data after removal:\n")
cat("Mean HI in Low RITH group:", mean(data_cleaned$Heterogeneity_Index[data_cleaned$RITH_Phenotype == "Low RITH"], na.rm = TRUE), "\n")
cat("Mean HI in High RITH group:", mean(data_cleaned$Heterogeneity_Index[data_cleaned$RITH_Phenotype == "High RITH"], na.rm = TRUE), "\n")

# Ensure all samples in High RITH group have HI >= 0
if (all(data_cleaned$Heterogeneity_Index[data_cleaned$RITH_Phenotype == "High RITH"] >= 0)) {
  cat("✅ Validation passed: all High RITH samples have HI >= 0\n")
} else {
  warning("⚠️ Some High RITH samples still have HI < 0!")
}

# ==================== 14. Save cleaned data ====================
cat("\n========== Save cleaned data ==========\n")

# Save complete cleaned data (including original feature values)
write_csv(data_cleaned, "RITH_Data_Cleaned.csv")
cat("Cleaned complete data saved to: RITH_Data_Cleaned.csv\n")
