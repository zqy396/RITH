
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

setwd('/mnt/data4/表型研究2026/ccRCC/RHWU/3.HI箱线图')
# ==================== 1. 读取数据 ====================
data <- read_csv("/mnt/data4/表型研究2026/ccRCC/RHWU/RITH_Data_Cleaned.csv")
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
  "LHR_feature" = "#FF9800",  # 蓝色：值越低表示异质性越高
  "HHR_feature" = "#E41A1C"   # 红色：值越高表示异质性越高
)

hi_colors <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(50)
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

heatmap_data_scaled <- scale(heatmap_data_scaled)
heatmap_data_scaled <- scale(heatmap_data_scaled)


# ==================== 17. HI的小提琴图和箱线图 ====================
cat("\n========== HI的小提琴图可视化 ==========\n")

# 使用您指定的配色：蓝色=低RITH，红色=高RITH
phenotype_colors <- c("Low RITH" = "#377EB8",   # 蓝色
                      "High RITH" = "#E41A1C")  # 红色

# 图1：HI的小提琴图 - 调大字号
p_hi_violin <- ggplot(data, aes(x = RITH_Phenotype, y = Heterogeneity_Index, fill = RITH_Phenotype)) +
  geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = RITH_Phenotype), width = 0.01, size = 1, alpha = 0.5) +  # 点大小从0.5调大到1
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # 从14调大到20
    legend.position = "none",
    axis.title.x = element_text(size = 18, face = "bold"),  # X轴标题字号
    axis.title.y = element_text(size = 18, face = "bold"),  # Y轴标题字号
    axis.text = element_text(size = 16),  # 从11调大到16
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# 手动添加显著性标记
y_max <- max(data$Heterogeneity_Index, na.rm = TRUE)
y_min <- min(data$Heterogeneity_Index, na.rm = TRUE)
y_range <- y_max - y_min

signif_label <- ifelse(hi_p_value < 0.001, "***", 
                       ifelse(hi_p_value < 0.01, "**",
                              ifelse(hi_p_value < 0.05, "*", "ns")))

# 添加显著性标记（调整大小）
p_hi_violin <- p_hi_violin +
  geom_segment(x = 1, xend = 2, 
               y = y_max + 0.1 * y_range, 
               yend = y_max + 0.1 * y_range,
               color = "black", linewidth = 1) +  # linewidth从0.8调大到1
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
           size = 8, fontface = "bold", color = "black")  # size从5调大到8

# 图2：HI的箱线图 - 调大字号
p_hi_box <- ggplot(data, aes(x = RITH_Phenotype, y = Heterogeneity_Index, fill = RITH_Phenotype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 19, outlier.size = 3, outlier.alpha = 0.5) +  # outlier.size从2调大到3
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.4, color = "black") +  # size从1.5调大到2
  scale_fill_manual(values = phenotype_colors) +
  labs(
    title = "Heterogeneity Index (HI) Boxplot",
    x = "RITH Phenotype",
    y = "Heterogeneity Index (HI)",
    fill = "Phenotype"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # 从14调大到20
    legend.position = "none",
    axis.title.x = element_text(size = 18, face = "bold"),  # X轴标题字号
    axis.title.y = element_text(size = 18, face = "bold"),  # Y轴标题字号
    axis.text = element_text(size = 16),  # 从11调大到16
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(..y.., 1)),
    vjust = -1.5,
    size = 5,  # 从4调大到5
    fontface = "bold"
  ) +
  geom_segment(x = 1, xend = 2, 
               y = y_max + 0.1 * y_range, 
               yend = y_max + 0.1 * y_range,
               color = "black", linewidth = 1) +  # linewidth从0.8调大到1
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
           size = 8, fontface = "bold", color = "black")  # size从5调大到8

# 图3：HI的密度图 - 调大字号
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # 从14调大到20
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),  # 图例标题字号
    legend.text = element_text(size = 14),  # 图例文字字号
    axis.title.x = element_text(size = 18, face = "bold"),  # X轴标题字号
    axis.title.y = element_text(size = 18, face = "bold"),  # Y轴标题字号
    axis.text = element_text(size = 16),  # 从11调大到16
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  ) +
  geom_vline(
    data = data %>% group_by(RITH_Phenotype) %>% summarise(median_HI = median(Heterogeneity_Index)),
    aes(xintercept = median_HI, color = RITH_Phenotype),
    linetype = "dashed",
    linewidth = 1.2,  # 从1调大到1.2
    show.legend = FALSE
  ) +
  scale_color_manual(values = phenotype_colors)

# 图4：HI的直方图（分面显示）- 调大字号
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # 从14调大到20
    legend.position = "none",
    strip.text = element_text(size = 16, face = "bold"),  # 分面标题从12调大到16
    axis.title.x = element_text(size = 18, face = "bold"),  # X轴标题字号
    axis.title.y = element_text(size = 18, face = "bold"),  # Y轴标题字号
    axis.text = element_text(size = 16),  # 从11调大到16
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# 保存单个图（PDF格式，高质量）
ggsave("HI_Violin_Plot.pdf", p_hi_violin, width = 8, height = 7, dpi = 300)  # 调整尺寸
ggsave("HI_Box_Plot.pdf", p_hi_box, width = 8, height = 7, dpi = 300)
ggsave("HI_Density_Plot.pdf", p_hi_density, width = 10, height = 7, dpi = 300)
ggsave("HI_Histogram.pdf", p_hi_histogram, width = 8, height = 10, dpi = 300)

# 同时保存PNG版本以备快速查看
ggsave("HI_Violin_Plot.png", p_hi_violin, width = 8, height = 7, dpi = 300)
ggsave("HI_Box_Plot.png", p_hi_box, width = 8, height = 7, dpi = 300)
ggsave("HI_Density_Plot.png", p_hi_density, width = 10, height = 7, dpi = 300)
ggsave("HI_Histogram.png", p_hi_histogram, width = 8, height = 10, dpi = 300)


cat("✅ HI可视化图已保存\n")

# ==================== 18. HI统计摘要 ====================
cat("\n========== HI统计摘要 ==========\n")

# 计算详细的统计摘要
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

cat("HI的详细统计摘要：\n")
print(hi_stats)

# 保存统计摘要
write_csv(hi_stats, "HI_Statistics_Summary.csv")
cat("HI统计摘要已保存至：HI_Statistics_Summary.csv\n")

# 计算效应量（Cohen's d）
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

cat("\n效应量分析：\n")
cat("Cohen's d =", round(cohens_d, 3), "\n")
if (abs(cohens_d) < 0.2) {
  cat("效应量：很小\n")
} else if (abs(cohens_d) < 0.5) {
  cat("效应量：小\n")
} else if (abs(cohens_d) < 0.8) {
  cat("效应量：中等\n")
} else {
  cat("效应量：大\n")
}

cat("\n✅ HI可视化分析完成！\n")