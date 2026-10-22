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
setwd('/mnt/data4/表型研究2026/ccRCC/RHWU/4.特征和HI相关性')

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

# ==================== 2. 特征名称处理函数 ====================
identify_feature_name <- function(feature_name) {
  parts <- str_split(feature_name, "_", simplify = TRUE)
  if (length(parts) >= 3) {
    return(parts[length(parts)])
  }
  return(feature_name)
}

# ==================== 3. 计算特征与HI的相关性 ====================
cat("\n========== 计算特征与HI的相关性 ==========\n")

# 首先需要确保有HI列
if (!"Heterogeneity_Index" %in% colnames(data)) {
  stop("数据中没有找到Heterogeneity_Index列，请先计算HI")
}

# 提取所有特征（排除非特征列）
non_feature_cols <- c("Image", "RITH_Phenotype", "Heterogeneity_Index", "Cluster_Label")
all_features <- setdiff(colnames(data), non_feature_cols)

# 计算每个特征与HI的Spearman相关性
cor_results <- data.frame()

for (feat in all_features) {
  cor_val <- cor(data[[feat]], data$Heterogeneity_Index, 
                 method = "spearman", use = "complete.obs")
  
  cor_results <- rbind(cor_results, data.frame(
    Feature = feat,
    Feature_Name = identify_feature_name(feat),
    Correlation = cor_val
  ))
}

# 按相关性排序
cor_results <- cor_results[order(-cor_results$Correlation), ]
rownames(cor_results) <- NULL

# 保存相关性结果
write_csv(cor_results, "All_Features_HI_Correlation.csv")
cat("特征与HI相关性分析已保存至：All_Features_HI_Correlation.csv\n")

# ==================== 4. 选取正负前5的特征 ====================
cat("\n========== 选取正负相关系数前5的特征 ==========\n")

# 选取正相关前5的特征
top_pos_features <- cor_results %>%
  filter(Correlation > 0) %>%
  arrange(desc(Correlation)) %>%
  head(5)

# 选取负相关前5的特征（相关性最小，即最负的）
top_neg_features <- cor_results %>%
  filter(Correlation < 0) %>%
  arrange(Correlation) %>%  # 从最小（最负）开始
  head(5)

cat("\n正相关前5的特征：\n")
print(top_pos_features)

cat("\n负相关前5的特征：\n")
print(top_neg_features)

# 合并选中的特征
selected_features <- rbind(top_pos_features, top_neg_features)
selected_feature_names <- selected_features$Feature_Name
selected_feature_codes <- selected_features$Feature

# ==================== 5. 准备绘图数据（使用完整特征名） ====================
cat("\n========== 准备绘图数据（使用完整特征名） ==========\n")

# 为selected_features添加完整信息
selected_features <- selected_features %>%
  mutate(
    # 获取滤波器类型
    Filter_Type = sapply(Feature, identify_filter_type),
    # 获取特征类别
    Feature_Class = sapply(Feature, identify_feature_class),
    # 获取原始特征名
    Base_Feature_Name = sapply(Feature, identify_feature_name),
    # 创建完整显示名（包含滤波器信息）
    Full_Feature_Name = Feature,  # 直接用原始特征名
    # 创建更易读的显示名
    Display_Full = ifelse(
      Filter_Type == "original",
      paste0(Base_Feature_Name, " (Original)"),
      paste0(Base_Feature_Name, " (", Filter_Type, ")")
    ),
    # 如果还有重复，添加特征类别
    Display_Full = paste0(Display_Full, " [", Feature_Class, "]")
  )

# 检查唯一性
cat("检查完整特征名唯一性：\n")
dup_check <- selected_features %>%
  group_by(Display_Full) %>%
  summarise(count = n()) %>%
  filter(count > 1)

if (nrow(dup_check) > 0) {
  cat("发现重复的显示名，添加序号：\n")
  print(dup_check)
  
  selected_features <- selected_features %>%
    group_by(Display_Full) %>%
    mutate(
      dup_num = row_number(),
      Display_Full = ifelse(dup_num > 1,
                            paste0(Display_Full, " #", dup_num),
                            Display_Full)
    ) %>%
    ungroup()
}

# 创建绘图数据框
plot_data <- data.frame(
  HI = data$Heterogeneity_Index,
  RITH_Phenotype = data$RITH_Phenotype
)

# 添加选中的特征
for (i in 1:nrow(selected_features)) {
  feat_code <- selected_features$Feature[i]
  plot_data[[feat_code]] <- data[[feat_code]]
}

cat("绘图数据准备完成，包含", ncol(plot_data) - 2, "个特征\n")

# ==================== 6. 创建三种显示名 ====================
cat("\n========== 创建三种特征显示名 ==========\n")

selected_features <- selected_features %>%
  mutate(
    # 方法1：完整特征名 + 相关性
    Display_Name1 = paste0(
      Display_Full,
      "\nr = ", round(Correlation, 3), 
      ifelse(Correlation > 0, " (+) ", " (-) ")
    ),
    
    # 方法2：分组显示（正负相关分组）
    Display_Name2 = paste0(
      ifelse(Correlation > 0, "↗ ", "↘ "),
      Display_Full,
      "\n(r = ", round(Correlation, 3), ")"
    ),
    
    # 方法3：详细显示（包含所有信息）
    Display_Name3 = paste0(
      "Feature ", sprintf("%02d", row_number()), ": ",
      Base_Feature_Name, "\n",
      "Filter: ", Filter_Type, 
      " | Class: ", Feature_Class,
      "\nr = ", round(Correlation, 3),
      ifelse(Correlation > 0, " (Positive)", " (Negative)")
    )
  )

# 打印特征信息
cat("\n选择的特征信息：\n")
for (i in 1:nrow(selected_features)) {
  cat("\n", i, ". ", selected_features$Feature[i], "\n", sep = "")
  cat("   简名: ", selected_features$Base_Feature_Name[i], "\n", sep = "")
  cat("   滤波器: ", selected_features$Filter_Type[i], "\n", sep = "")
  cat("   类别: ", selected_features$Feature_Class[i], "\n", sep = "")
  cat("   相关性: ", round(selected_features$Correlation[i], 3), "\n", sep = "")
}

# ==================== 7. 将数据转换为长格式 ====================
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = all_of(selected_features$Feature),
    names_to = "Feature_Code",
    values_to = "Feature_Value"
  ) %>%
  left_join(selected_features, by = c("Feature_Code" = "Feature"))

# ==================== 8. 去除极端值函数 ====================
cat("\n========== 去除极端值 ==========\n")

remove_outliers <- function(x, lower_quantile = 0.025, upper_quantile = 0.975) {
  lower_bound <- quantile(x, lower_quantile, na.rm = TRUE)
  upper_bound <- quantile(x, upper_quantile, na.rm = TRUE)
  x[x < lower_bound] <- NA
  x[x > upper_bound] <- NA
  return(x)
}

plot_data_clean <- plot_data_long %>%
  group_by(Feature_Code) %>%
  mutate(
    Feature_Value_Clean = remove_outliers(Feature_Value)
  ) %>%
  ungroup() %>%
  filter(!is.na(Feature_Value_Clean))

# ==================== 9. 方法1：基本散点图（完整特征名） ====================
cat("\n========== 方法1：基本散点图 ==========\n")

# 按相关性排序
plot_data_clean$Display_Name1 <- factor(
  plot_data_clean$Display_Name1,
  levels = unique(plot_data_clean$Display_Name1[order(-plot_data_clean$Correlation)])
)

p_method1 <- ggplot(plot_data_clean, 
                    aes(x = Feature_Value_Clean, y = HI, color = RITH_Phenotype)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.7) +
  facet_wrap(~ Display_Name1, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c("Low RITH" = "#377EB8", "High RITH" = "#E41A1C")) +
  labs(
    title = "",
    subtitle = paste("Top", nrow(selected_features), "features correlated with Heterogeneity Index"),
    x = "Feature Value (95% central range)",
    y = "Heterogeneity Index (HI)",
    color = "RITH Phenotype"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b = 15)),
    strip.text = element_text(size = 9, face = "bold", lineheight = 1.1),
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )

# 保存图形
ggsave("HI_Scatter_Method1_FullNames.pdf", p_method1, 
       width = 16, height = 12, dpi = 300)
cat("方法1图形已保存：HI_Scatter_Method1_FullNames.png\n")



# ==================== 12. 保存特征信息表（修复版） ====================
cat("\n========== 保存特征信息 ==========\n")

# 先使用mutate创建新列，再select
feature_info_table <- selected_features %>%
  mutate(
    Direction = ifelse(Correlation > 0, "Positive", "Negative"),
    Correlation_Strength = case_when(
      abs(Correlation) >= 0.7 ~ "Very Strong",
      abs(Correlation) >= 0.5 ~ "Strong",
      abs(Correlation) >= 0.3 ~ "Moderate",
      TRUE ~ "Weak"
    )
  ) %>%
  select(
    Feature_Code = Feature,
    Base_Name = Base_Feature_Name,
    Filter_Type,
    Feature_Class,
    Correlation,
    Direction,
    Correlation_Strength
  ) %>%
  arrange(desc(Correlation))

write_csv(feature_info_table, "Selected_Features_Detailed_Info.csv")
cat("详细特征信息已保存：Selected_Features_Detailed_Info.csv\n")

# 显示特征信息
cat("\n选择的特征详细列表：\n")
print(feature_info_table)