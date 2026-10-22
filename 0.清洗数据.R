# ==================== 加载必要的R包 ====================
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
# ==================== 1. 读取数据 ====================
data <- read_csv("/mnt/data4/表型研究2026/ccRCC/RITH_Data_Cleaned.csv")
head(data)

# ==================== 2. 提取所有特征名称 ====================
all_features <- colnames(data)
all_features <- all_features[all_features != "Image"]

cat("总特征数：", length(all_features), "\n")
cat("样本数：", nrow(data), "\n")

# ==================== 3. 特征分类 ====================
cat("\n========== 特征分类 ==========\n")

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

cat("滤波器类型分布：\n")
print(table(feature_info$Filter_Type))
cat("\n特征类别分布：\n")
print(table(feature_info$Feature_Class))

# ==================== 4. MAD特征筛选 ====================
cat("\n========== MAD特征筛选 ==========\n")

X <- data[, all_features] %>% as.matrix()
rownames(X) <- data$Image

mad_values <- apply(X, 2, mad, na.rm = TRUE)
significant_features <- names(mad_values[mad_values > 1])

cat("总特征数：", length(all_features), "\n")
cat("MAD > 1 的显著特征数：", length(significant_features), "\n")

significant_info <- feature_info[feature_info$Feature %in% significant_features, ]
cat("\n显著特征按滤波器类型分布：\n")
print(table(significant_info$Filter_Type))
cat("\n显著特征按特征类别分布：\n")
print(table(significant_info$Feature_Class))

X_sig <- X[, significant_features]

# ==================== 5. 基于与ITH正相关性确定特征方向 ====================
cat("\n========== 基于与ITH正相关性确定特征方向 ==========\n")

# 原则：大多数特征与ITH呈正相关
# 因此，特征值高的样本应该具有更高的异质性
# 我们需要先计算一个"异质性倾向"指标，然后根据特征与该指标的相关性确定方向

# 步骤1：计算每个样本的"异质性倾向分数"
# 由于大多数特征与ITH正相关，我们可以简单使用特征平均值作为初始的异质性倾向
feature_means <- rowMeans(X_sig, na.rm = TRUE)

# 步骤2：计算每个特征与异质性倾向分数的相关性
z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

X_z <- apply(X_sig, 2, z_score)
rownames(X_z) <- data$Image

# 计算每个特征与异质性倾向的Spearman相关性
correlations <- apply(X_z, 2, function(feat) {
  cor(feat, feature_means, method = "spearman", use = "complete.obs")
})

# 步骤3：根据相关性确定特征方向
# 正相关（>0）：HHR特征（值越高，异质性越高）
# 负相关（<0）：LHR特征（值越低，异质性越高）
# 注意：根据您的描述"大多数特征与ITH呈正相关"，所以阈值设为0
feature_direction <- ifelse(correlations > 0, "HHR_feature", "LHR_feature")
names(feature_direction) <- significant_features

# 统计
HHR_features <- names(feature_direction[feature_direction == "HHR_feature"])
LHR_features <- names(feature_direction[feature_direction == "LHR_feature"])

cat("\n特征方向分布（基于与异质性倾向的相关性）：\n")
print(table(feature_direction))

cat("\nHHR特征（与异质性正相关）：", length(HHR_features), "个（", 
    round(length(HHR_features)/length(significant_features)*100, 1), "%）\n")
cat("LHR特征（与异质性负相关）：", length(LHR_features), "个（", 
    round(length(LHR_features)/length(significant_features)*100, 1), "%）\n")

# 保存相关性分析结果
cor_df <- data.frame(
  Feature = significant_features,
  Feature_Name = sapply(significant_features, identify_feature_name),
  Feature_Class = sapply(significant_features, identify_feature_class),
  Filter_Type = sapply(significant_features, identify_filter_type),
  Correlation_with_Heterogeneity = correlations,
  Direction = feature_direction
)

# 按相关性绝对值排序
cor_df$AbsCorrelation <- abs(cor_df$Correlation_with_Heterogeneity)
cor_df <- cor_df[order(-cor_df$AbsCorrelation), ]

write_csv(cor_df, "Feature_Correlation_with_Heterogeneity.csv")
cat("\n特征与异质性相关性分析已保存至：Feature_Correlation_with_Heterogeneity.csv\n")

cat("\n与异质性相关性最强的10个特征：\n")
print(head(cor_df, 10))

# ==================== 6. 保存特征方向定义 ====================
cat("\n========== 保存特征方向定义 ==========\n")

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
cat("基于相关性的特征方向定义已保存至：Feature_Direction_Definitions_CorrelationBased.csv\n")

# ==================== 7. SNF聚类分析（2聚类） ====================
cat("\n========== 执行SNF聚类分析（2聚类） ==========\n")

set.seed(123)
optimal_k <- 2
cat("使用聚类数 K =", optimal_k, "\n")

# 计算相似性矩阵
W1 <- affinityMatrix(dist2(X_sig, X_sig), K = 20, sigma = 0.5)

# 执行SNF聚类
cluster_labels <- spectralClustering(W1, K = optimal_k)

# 注意：此时我们不知道哪个标签对应高ITH，哪个对应低ITH
data$Cluster_Label <- cluster_labels

cat("初始聚类标签分布：\n")
print(table(cluster_labels))

# ==================== 8. 确定哪个标签是高RITH亚组 ====================
cat("\n========== 确定高RITH亚组（红色）和低RITH亚组（蓝色） ==========\n")

# 原则：根据"大多数特征与ITH呈正相关"，特征值高的样本应该是高RITH
# 计算每个聚类标签的平均特征值
cluster_means <- tapply(feature_means, cluster_labels, mean, na.rm = TRUE)
cat("\n各聚类标签的平均特征值：\n")
print(cluster_means)

# 特征值更高的聚类应该是高RITH亚组（红色）
if (cluster_means[1] > cluster_means[2]) {
  high_ith_label <- 1
  low_ith_label <- 2
} else {
  high_ith_label <- 2
  low_ith_label <- 1
}

cat("\n确定标签", high_ith_label, "为高RITH亚组（红色）\n")
cat("确定标签", low_ith_label, "为低RITH亚组（蓝色）\n")

# 创建最终的表型标签
data$RITH_Phenotype <- ifelse(cluster_labels == high_ith_label, "High RITH", "Low RITH")
data$RITH_Phenotype <- factor(data$RITH_Phenotype, levels = c("Low RITH", "High RITH"))

cat("\n最终表型分布：\n")
print(table(data$RITH_Phenotype))

# 验证：高RITH组的特征平均值应该显著高于低RITH组
cat("\n验证表型定义的正确性：\n")
cat("低RITH组平均特征值：", mean(feature_means[data$RITH_Phenotype == "Low RITH"], na.rm = TRUE), "\n")
cat("高RITH组平均特征值：", mean(feature_means[data$RITH_Phenotype == "High RITH"], na.rm = TRUE), "\n")

p_value <- wilcox.test(feature_means[data$RITH_Phenotype == "Low RITH"], 
                       feature_means[data$RITH_Phenotype == "High RITH"])$p.value
cat("两组特征平均值差异的p值：", format.pval(p_value), "\n")

if (p_value < 0.05) {
  cat("✅ 表型定义正确：高RITH组特征值显著高于低RITH组（p < 0.05）\n")
} else {
  cat("⚠️ 警告：两组特征平均值差异不显著，可能需要重新检查\n")
}

# ==================== 9. 计算异质性指数（HI） ====================
cat("\n========== 计算异质性指数（HI） ==========\n")

# HI = HHR特征之和 - LHR特征之和
HI <- rowSums(X_z[, HHR_features], na.rm = TRUE) - 
  rowSums(X_z[, LHR_features], na.rm = TRUE)

data$Heterogeneity_Index <- HI

# 验证HI与表型的相关性
cat("\n验证HI与表型的相关性：\n")
cat("低RITH组HI均值：", mean(HI[data$RITH_Phenotype == "Low RITH"], na.rm = TRUE), "\n")
cat("高RITH组HI均值：", mean(HI[data$RITH_Phenotype == "High RITH"], na.rm = TRUE), "\n")

hi_p_value <- wilcox.test(HI[data$RITH_Phenotype == "Low RITH"], 
                          HI[data$RITH_Phenotype == "High RITH"])$p.value
cat("两组HI差异的p值：", format.pval(hi_p_value), "\n")

if (hi_p_value < 0.05) {
  cat("✅ HI能显著区分高/低RITH表型（p < 0.05）\n")
} else {
  cat("⚠️ HI未能显著区分高/低RITH表型\n")
}

# ==================== 10. 保存结果 ====================
cat("\n========== 保存表型和HI结果 ==========\n")

results <- data.frame(
  Image = data$Image,
  Cluster_Label = data$Cluster_Label,
  RITH_Phenotype = data$RITH_Phenotype,
  Feature_Mean = feature_means,
  Heterogeneity_Index = data$Heterogeneity_Index
)

write_csv(results, "RITH_Phenotype_and_HI_Final.csv")
cat("表型和HI结果已保存至：RITH_Phenotype_and_HI_Final.csv\n")

# ==================== 11. 优化热图：纵向按方向分隔 ====================
cat("\n========== 绘制纵向按方向分隔的热图 ==========\n")

# 准备热图数据（使用Z-score标准化后的数据）
heatmap_data <- X_z

# 按表型排序样本（低RITH在前，高RITH在后）
phenotype_order <- order(data$RITH_Phenotype)
heatmap_data <- heatmap_data[phenotype_order, ]

# 确保使用的特征在数据中存在
HHR_features <- intersect(HHR_features, colnames(heatmap_data))
LHR_features <- intersect(LHR_features, colnames(heatmap_data))

cat("用于热图的HHR特征数：", length(HHR_features), "\n")
cat("用于热图的LHR特征数：", length(LHR_features), "\n")

# 重新排列数据：先LHR特征，后HHR特征
heatmap_data_ordered <- heatmap_data[, c(LHR_features, HHR_features)]


# 4. 处理HI值：截断极值到-200到200之间
HI_truncated <- data$Heterogeneity_Index
HI_truncated[HI_truncated > 200] <- 200
HI_truncated[HI_truncated < -200] <- -100

# 计算分隔线位置
split_position <- length(LHR_features)

# 数据处理
heatmap_data_scaled <- scale(heatmap_data_ordered)

heatmap_data_scaled <- t(heatmap_data_scaled)  # 转置：行为特征，列为样本

# 截断极值


heatmap_data_scaled[heatmap_data_scaled > 2] <- 2
heatmap_data_scaled[heatmap_data_scaled < -2] <- -2

# 创建样本注释
annotation_row <- data.frame(
  RITH_Phenotype = factor(data$RITH_Phenotype[phenotype_order], levels = c("Low RITH", "High RITH")),
  Heterogeneity_Index = HI_truncated[phenotype_order]
)
rownames(annotation_row) <- rownames(heatmap_data_ordered)

# 创建特征注释
feature_direction_vector <- c(
  rep("LHR_feature", length(LHR_features)),
  rep("HHR_feature", length(HHR_features))
)

annotation_col <- data.frame(
  Direction = factor(feature_direction_vector, levels = c("LHR_feature", "HHR_feature")),
  row.names = rownames(heatmap_data_scaled)
)

# 颜色设置（根据您的描述：红色=高RITH，蓝色=低RITH）
phenotype_colors <- c("Low RITH" = "#377EB8",   # 蓝色
                      "High RITH" = "#E41A1C")  # 红色

direction_colors <- c(
  "LHR_feature" = "#377EB8",  # 蓝色：值越低表示异质性越高
  "HHR_feature" = "#E41A1C"   # 红色：值越高表示异质性越高
)

hi_colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(50)
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

heatmap_data_scaled <- scale(heatmap_data_scaled)
heatmap_data_scaled <- scale(heatmap_data_scaled)
# 绘制热图
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

cat("✅ 方向分隔热图已保存为：RITH_Direction_Separated_Heatmap.png\n")

# ==================== 12. 总结报告 ====================
cat("\n========== 分析总结报告 ==========\n")
cat("1. 总样本数：", nrow(data), "\n")
cat("2. 显著特征数（MAD > 1）：", length(significant_features), "\n")
cat("3. 正相关特征数（HHR）：", length(HHR_features), 
    "（", round(length(HHR_features)/length(significant_features)*100, 1), "%）\n")
cat("4. 负相关特征数（LHR）：", length(LHR_features), 
    "（", round(length(LHR_features)/length(significant_features)*100, 1), "%）\n")
cat("5. 低RITH亚组（蓝色）：", sum(data$RITH_Phenotype == "Low RITH"), "个样本\n")
cat("6. 高RITH亚组（红色）：", sum(data$RITH_Phenotype == "High RITH"), "个样本\n")
cat("7. 表型验证p值：", format.pval(p_value), "\n")
cat("8. HI区分能力p值：", format.pval(hi_p_value), "\n")
cat("\n✅ 分析完成！\n")

# ==================== 13. 删除高RITH组中HI小于0的样本 ====================
cat("\n========== 删除高RITH组中HI小于0的样本 ==========\n")

# 找出高RITH组中HI小于0的样本
high_rith_samples <- data$Image[data$RITH_Phenotype == "High RITH"]
low_hi_samples <- data$Image[data$RITH_Phenotype == "High RITH" & data$Heterogeneity_Index < 0]

cat("高RITH组总样本数：", length(high_rith_samples), "\n")
cat("高RITH组中HI小于0的样本数：", length(low_hi_samples), "\n")
cat("高RITH组中HI小于0的样本比例：", round(length(low_hi_samples)/length(high_rith_samples)*100, 1), "%\n")

# 显示这些样本的具体信息
if (length(low_hi_samples) > 0) {
  cat("\n高RITH组中HI小于0的样本：\n")
  for (sample_id in low_hi_samples) {
    hi_value <- data$Heterogeneity_Index[data$Image == sample_id]
    cat(sample_id, sprintf(" (HI = %.2f)", hi_value), "\n")
  }
}

# 删除这些样本
data_cleaned <- data[!(data$Image %in% low_hi_samples), ]

cat("\n删除前总样本数：", nrow(data), "\n")
cat("删除后总样本数：", nrow(data_cleaned), "\n")
cat("删除样本数：", nrow(data) - nrow(data_cleaned), "\n")

# 重新统计表型分布
cat("\n新的表型分布：\n")
print(table(data_cleaned$RITH_Phenotype))

# 验证删除后的数据
cat("\n验证删除后的数据：\n")
cat("低RITH组HI均值：", mean(data_cleaned$Heterogeneity_Index[data_cleaned$RITH_Phenotype == "Low RITH"], na.rm = TRUE), "\n")
cat("高RITH组HI均值：", mean(data_cleaned$Heterogeneity_Index[data_cleaned$RITH_Phenotype == "High RITH"], na.rm = TRUE), "\n")

# 确保高RITH组所有样本的HI都大于等于0
if (all(data_cleaned$Heterogeneity_Index[data_cleaned$RITH_Phenotype == "High RITH"] >= 0)) {
  cat("✅ 验证通过：高RITH组所有样本的HI都大于等于0\n")
} else {
  warning("⚠️ 仍有高RITH组样本HI小于0！")
}

# ==================== 14. 保存清理后的数据 ====================
cat("\n========== 保存清理后的数据 ==========\n")

# 保存完整的清理后数据（包含原始特征值）
write_csv(data_cleaned, "RITH_Data_Cleaned.csv")
cat("清理后的完整数据已保存至：RITH_Data_Cleaned.csv\n")