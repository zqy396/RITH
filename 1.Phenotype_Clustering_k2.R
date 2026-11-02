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

setwd('/mnt/data4/表型研究2026/ccRCC/RHWU/1.表型聚类')
# ==================== 1. Read data ====================
data <- read_csv("/mnt/data4/表型研究2026/ccRCC/RHWU/RITH_Data_Cleaned.csv")
head(data)

# ... (Sections 2-6 are identical to 0.Data_Cleaning.R, use the English code above) ...

# ==================== 7. Perform SNF clustering analysis (2 clusters) ====================
cat("\n========== Perform SNF clustering analysis (2 clusters) ==========\n")

set.seed(1234)
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

# ... (Sections 8-10 are identical to 0.Data_Cleaning.R, use the English code above) ...

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
  "LHR_feature" = "#FF9800",  # Orange
  "HHR_feature" = "#E41A1C"   # Red
)

hi_colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(50)
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

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
  filename = "RITH_Direction_Separated_Heatmap.pdf",
  width = 10,
  height = 10
)

cat("✅ Direction-separated heatmap saved as: RITH_Direction_Separated_Heatmap.pdf\n")

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
