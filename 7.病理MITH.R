setwd('F:\\桌面\\ITH表型\\ccRCC\\RHWU\\7.病理MITH')

# 清除工作空间
rm(list = ls())

# 加载必要的包
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(tidyr)
library(dplyr)

# 加载数据
df <- read.csv("Feats_T.csv")

# 提取形态特征列（以 Morph_ 开头）
morph_cols <- names(df)[grepl("^Morph_", names(df))]
cat("找到", length(morph_cols), "个形态特征\n")

# ==================== 第一部分：创建输出目录 ====================
output_dir <- "Morphological_Analysis_Results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# 设置图片尺寸（调大以适应更大的字号）
plot_width <- 16  # 从16调大到18
plot_height <- 12  # 从12调大到14
dpi <- 300

# ==================== 第二部分：HI和每个形态学特征的相关性图 ====================
cat("\n=== 绘制HI和每个形态学特征的相关性图 ===\n")

# 创建数据框存储相关性结果
correlation_results <- data.frame(
  Feature = character(),
  Correlation = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# 为每个形态特征计算与HI的相关性并绘制散点图
cor_plots <- list()

for (i in seq_along(morph_cols)) {
  feature <- morph_cols[i]
  
  # 计算相关性
  cor_test <- cor.test(df[[feature]], df$Heterogeneity_Index, 
                       use = "complete.obs", method = "pearson")
  
  # 保存结果
  correlation_results <- rbind(correlation_results, data.frame(
    Feature = feature,
    Correlation = cor_test$estimate,
    p_value = cor_test$p.value
  ))
  
  # 创建散点图 - 调大字号
  p <- ggplot(df, aes_string(x = "Heterogeneity_Index", y = feature)) +
    geom_point(aes(color = RITH_Phenotype), alpha = 0.6, size = 2.5) +  # 点大小从2调大到2.5
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", size = 1) +  # 线宽从0.8调大到1
    theme_minimal(base_size = 14) +  # 基础字号从12调大到14
    labs(title = feature,
         subtitle = sprintf("r = %.3f, p = %.3e", 
                            cor_test$estimate, 
                            cor_test$p.value),
         x = "Heterogeneity Index (HI)",
         y = "Feature Value") +
    scale_color_manual(values = c("High RITH" = "#E41A1C", "Low RITH" = "#377EB8")) +
    theme(legend.position = "none",
          plot.title = element_text(size = 16, face = "bold"),  # 从10调大到14
          plot.subtitle = element_text(size = 14),  # 从9调大到12
          axis.title = element_text(size = 16, face = "bold"),  # 新增坐标轴标题字号
          axis.text = element_text(size = 14))  # 新增坐标轴刻度字号
  
  cor_plots[[feature]] <- p
  
  # 每处理10个特征显示一次进度
  if (i %% 10 == 0) {
    cat(sprintf("已处理 %d/%d 个特征...\n", i, length(morph_cols)))
  }
}

# 将所有相关性散点图组合并保存
cat("组合所有相关性散点图...\n")

# 计算需要多少行多少列
n_features <- length(morph_cols)
n_cols <- 4
n_rows <- ceiling(n_features / n_cols)

# 创建组合图
all_cor_plots <- do.call(grid.arrange, c(cor_plots, ncol = n_cols))

# 保存组合图
ggsave(file.path(output_dir, "01_HI_vs_Morphological_Features_Correlation.pdf"),
       all_cor_plots,
       width = plot_width,
       height = n_rows * 4,  # 从3调大到4，适应更大的字号
       dpi = dpi,
       limitsize = FALSE)

cat("已保存: 01_HI_vs_Morphological_Features_Correlation.pdf\n")

# ==================== 第三部分：所有形态学特征的小提琴分布图（带统计学分析）====================
cat("\n=== 绘制所有形态学特征的小提琴分布图（带统计学分析）===\n")

# 首先需要将宽数据转换为长数据
df_long <- df %>%
  select(id, RITH_Phenotype, all_of(morph_cols)) %>%
  pivot_longer(cols = all_of(morph_cols),
               names_to = "Feature",
               values_to = "Value")

# 创建小提琴图列表（带统计学分析）
violin_plots <- list()

# 为每个特征创建小提琴图
unique_features <- unique(df_long$Feature)

for (i in seq_along(unique_features)) {
  feature <- unique_features[i]
  
  # 提取该特征的数据
  feature_data <- df_long[df_long$Feature == feature, ]
  
  # 检查是否有足够的数据
  high_data <- feature_data$Value[feature_data$RITH_Phenotype == "High RITH"]
  low_data <- feature_data$Value[feature_data$RITH_Phenotype == "Low RITH"]
  
  # 创建小提琴图 - 带统计学分析（调大字号）
  p <- ggplot(feature_data, aes(x = RITH_Phenotype, y = Value, fill = RITH_Phenotype)) +
    geom_violin(alpha = 0.7, trim = TRUE, color = "black", size = 0.8) +  # 线宽从0.5调大到0.8
    geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA, 
                 color = "black", size = 0.8) +  # 线宽从0.5调大到0.8
    geom_jitter(width = 0.1, alpha = 0.3, size = 1.2, color = "black") +  # 点大小从0.8调大到1.2
    
    # 添加统计学比较
    stat_compare_means(
      aes(label = ..p.signif..),
      method = "t.test",
      comparisons = list(c("High RITH", "Low RITH")),
      method.args = list(var.equal = FALSE),
      label.y = max(feature_data$Value, na.rm = TRUE) * 1.1,
      size = 10,  # 显著性标记字号从8调大到10
      bracket.size = 1.2,  # 括号线宽从1.0调大到1.2
      tip.length = 0.02,
      vjust = 0.5
    ) +
    
    labs(title = feature,
         x = "RITH Phenotype",
         y = "Feature Value") +
    
    scale_fill_manual(values = c("High RITH" = "#E41A1C", "Low RITH" = "#377EB8")) +
    theme_classic(base_size = 14) +  # 基础字号从12调大到14
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold"),  # 从10调大到14
      axis.title = element_text(size = 16, face = "bold"),  # 坐标轴标题字号
      axis.text = element_text(size = 14),  # 坐标轴刻度字号
      axis.line = element_line(color = "black", size = 0.8),  # 线宽从0.5调大到0.8
      axis.ticks = element_line(color = "black", size = 0.8),  # 线宽从0.5调大到0.8
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),  # 边框线宽从1调大到1.2
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    
    # 调整y轴范围以容纳显著性星号
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  violin_plots[[feature]] <- p
  
  # 每处理10个特征显示一次进度
  if (i %% 10 == 0) {
    cat(sprintf("已处理 %d/%d 个特征的小提琴图...\n", i, length(unique_features)))
  }
}

# 将所有小提琴图组合并保存
cat("组合所有小提琴图...\n")

# 创建组合图
all_violin_plots <- do.call(grid.arrange, c(violin_plots, ncol = n_cols))

# 保存组合图
ggsave(file.path(output_dir, "02_Morphological_Features_Violin_Plots_with_Stats.pdf"),
       all_violin_plots,
       width = plot_width,
       height = n_rows * 4,  # 从3.2调大到4，适应更大的字号
       dpi = dpi,
       limitsize = FALSE)

ggsave(file.path(output_dir, "02_Morphological_Features_Violin_Plots_with_Stats.png"),
       all_violin_plots,
       width = plot_width,
       height = n_rows * 4,
       dpi = dpi,
       limitsize = FALSE)

cat("已保存: 02_Morphological_Features_Violin_Plots_with_Stats.pdf\n")

# ==================== 第四部分：保存相关性结果和统计分析 ====================
cat("\n=== 保存分析结果 ===\n")

# 为小提琴图计算更详细的统计信息
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

# 计算组间差异的统计检验
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

# 合并所有结果
final_results <- correlation_results %>%
  left_join(ttest_results, by = c("Feature" = "Feature")) %>%
  arrange(p_value.x)

# 重命名列
colnames(final_results) <- c("Feature", "Correlation_with_HI", "Correlation_p_value",
                             "T_statistic", "T_test_p_value", "Significance", "Mean_difference_High_vs_Low")

# 保存结果
write.csv(final_results, 
          file.path(output_dir, "03_Statistical_Analysis_Results.csv"), 
          row.names = FALSE)

write.csv(violin_stats, 
          file.path(output_dir, "04_Violin_Plot_Statistics.csv"), 
          row.names = FALSE)

cat("已保存: 03_Statistical_Analysis_Results.csv\n")
cat("已保存: 04_Violin_Plot_Statistics.csv\n")

# ==================== 第五部分：创建最显著的特征的单独图表（带更详细的统计学分析）====================
cat("\n=== 创建最显著特征的单独图表 ===\n")

# 找出与HI相关性最强的特征
top_cor_features <- final_results %>%
  arrange(desc(abs(Correlation_with_HI))) %>%
  head(10) %>%
  pull(Feature)

# 找出组间差异最显著的特征
top_diff_features <- final_results %>%
  arrange(T_test_p_value) %>%
  head(10) %>%
  pull(Feature)

# 创建相关性最强的特征的散点图（调大字号）
cat("创建相关性最强的特征散点图...\n")
top_cor_plots <- list()

for (i in seq_along(top_cor_features)) {
  feature <- top_cor_features[i]
  cor_value <- final_results$Correlation_with_HI[final_results$Feature == feature]
  p_value_cor <- final_results$Correlation_p_value[final_results$Feature == feature]
  
  p <- ggplot(df, aes_string(x = "Heterogeneity_Index", y = feature)) +
    geom_point(aes(color = RITH_Phenotype), alpha = 0.7, size = 4) +  # 点大小从3调大到4
    geom_smooth(method = "lm", se = TRUE, color = "black", size = 1.2) +  # 线宽从1调大到1.2
    theme_minimal(base_size = 16) +  # 基础字号从14调大到16
    labs(title = sprintf("%s\nCorrelation: %.3f (p = %.3e)", feature, cor_value, p_value_cor),
         x = "Heterogeneity Index",
         y = "Feature Value") +
    scale_color_manual(values = c("High RITH" = "#E41A1C", "Low RITH" = "#377EB8")) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16, face = "bold"),  # 图例标题字号
          legend.text = element_text(size = 14),  # 图例文字字号
          plot.title = element_text(size = 16, face = "bold"),  # 标题字号
          axis.title = element_text(size = 16, face = "bold"),  # 坐标轴标题字号
          axis.text = element_text(size = 16))  # 坐标轴刻度字号
  
  ggsave(file.path(output_dir, sprintf("Top_Correlation_%02d_%s.pdf", i, feature)),
         p, width = 8, height = 6, dpi = dpi)  # 尺寸从8x6调大到9x7
  ggsave(file.path(output_dir, sprintf("Top_Correlation_%02d_%s.png", i, feature)),
         p, width = 8, height = 6, dpi = dpi)
  
  top_cor_plots[[feature]] <- p
}

# 组合并保存
if (length(top_cor_plots) > 0) {
  top_cor_grid <- do.call(grid.arrange, c(top_cor_plots, ncol = 2))
  ggsave(file.path(output_dir, "05_Top_Correlation_Features.pdf"),
         top_cor_grid, width = 18, height = 7 * ceiling(length(top_cor_plots)/2),  # 高度从6调大到7
         dpi = dpi, limitsize = FALSE)
  cat("已保存: 05_Top_Correlation_Features.pdf\n")
}

# 创建组间差异最显著的特征的小提琴图（带详细统计学分析，调大字号）
cat("创建组间差异最显著的特征小提琴图...\n")
top_diff_plots <- list()

for (i in seq_along(top_diff_features)) {
  feature <- top_diff_features[i]
  p_value <- final_results$T_test_p_value[final_results$Feature == feature]
  significance <- final_results$Significance[final_results$Feature == feature]
  
  # 提取数据
  plot_data <- df_long[df_long$Feature == feature, ]
  
  # 计算最大值以确定显著性星号位置
  y_max <- max(plot_data$Value, na.rm = TRUE)
  
  p <- ggplot(plot_data, aes(x = RITH_Phenotype, y = Value, fill = RITH_Phenotype)) +
    geom_violin(alpha = 0.7, color = "black", size = 0.8) +  # 线宽从0.5调大到0.8
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA, 
                 color = "black", size = 0.8) +  # 线宽从0.5调大到0.8
    geom_jitter(width = 0.05, alpha = 0.5, size = 2, color = "black") +  # 点大小从1.5调大到2
    
    # 添加显著性比较（带括号和星号）
    stat_compare_means(
      aes(label = paste0(..p.signif.., "\n", "p=", ..p.format..)),
      method = "t.test",
      comparisons = list(c("High RITH", "Low RITH")),
      method.args = list(var.equal = FALSE),
      label.y = y_max * 1.15,
      size = 10,  # 显著性标记字号从8调大到10
      bracket.size = 1.2,  # 括号线宽从1.0调大到1.2
      tip.length = 0.03,
      vjust = 0.5
    ) +
    
    theme_classic(base_size = 18) +  # 基础字号从14调大到16
    labs(title = feature,
         x = "RITH Phenotype",
         y = "Feature Value") +
    scale_fill_manual(values = c("High RITH" = "#E41A1C", "Low RITH" = "#377EB8")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 18),  # 从16调大到18
      axis.title = element_text(size = 18, face = "bold"),  # 坐标轴标题字号
      axis.text = element_text(size = 16),  # 坐标轴刻度字号
      axis.line = element_line(color = "black", size = 0.8),  # 线宽从0.5调大到0.8
      axis.ticks = element_line(color = "black", size = 0.8),  # 线宽从0.5调大到0.8
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),  # 边框线宽从1调大到1.2
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(output_dir, sprintf("Top_Difference_%02d_%s.pdf", i, feature)),
         p, width = 7, height = 7, dpi = dpi)  # 尺寸从8x7调大到9x8
  ggsave(file.path(output_dir, sprintf("Top_Difference_%02d_%s.png", i, feature)),
         p, width = 7, height = 7, dpi = dpi)  # 尺寸从8x7调大到9x8
  top_diff_plots[[feature]] <- p
}

# 组合并保存
if (length(top_diff_plots) > 0) {
  top_diff_grid <- do.call(grid.arrange, c(top_diff_plots, ncol = 2))
  ggsave(file.path(output_dir, "06_Top_Difference_Features_with_Stats.pdf"),
         top_diff_grid, width = 18, height = 8 * ceiling(length(top_diff_plots)/2),  # 高度从7调大到8
         dpi = dpi, limitsize = FALSE)
  cat("已保存: 06_Top_Difference_Features_with_Stats.pdf\n")
}
if (length(top_diff_plots) > 0) {
  ggsave(file.path(output_dir, "06_Top_Difference_Features_with_Stats.png"),
         top_diff_grid, width = 18, height = 8 * ceiling(length(top_diff_plots)/2),
         dpi = dpi, limitsize = FALSE)
  cat("已保存: 06_Top_Difference_Features_with_Stats.png\n")
}
# ==================== 第六部分：创建显著性热图（可选）====================
cat("\n=== 创建显著性热图 ===\n")

# 准备热图数据
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

# 保存热图数据
write.csv(heatmap_data, 
          file.path(output_dir, "07_Significance_Heatmap_Data.csv"), 
          row.names = FALSE)

cat("已保存: 07_Significance_Heatmap_Data.csv\n")

# ==================== 第七部分：生成分析报告 ====================
cat("\n=== 生成分析报告 ===\n")

sink(file.path(output_dir, "08_Analysis_Report.txt"))

cat("========================================\n")
cat("    形态学特征分析报告（带统计学分析）\n")
cat("========================================\n\n")

cat("分析日期:", format(Sys.Date(), "%Y年%m月%d日"), "\n")
cat("分析时间:", format(Sys.time(), "%H:%M:%S"), "\n\n")

cat("=== 数据集信息 ===\n")
cat("总样本数:", nrow(df), "\n")
cat("High RITH 样本数:", sum(df$RITH_Phenotype == "High RITH"), "\n")
cat("Low RITH 样本数:", sum(df$RITH_Phenotype == "Low RITH"), "\n")
cat("分析的形态特征数:", length(morph_cols), "\n\n")

cat("=== 统计学显著性符号说明 ===\n")
cat("ns: p ≥ 0.05\n")
cat("*: p < 0.05\n")
cat("**: p < 0.01\n")
cat("***: p < 0.001\n\n")

cat("=== 最相关的特征（与HI）===\n")
cat("前10个与Heterogeneity Index最相关的特征:\n")
print(head(final_results[order(abs(final_results$Correlation_with_HI), decreasing = TRUE), 
                         c("Feature", "Correlation_with_HI", "Correlation_p_value", "Significance")], 10))

cat("\n=== 组间差异最显著的特征 ===\n")
cat("前10个在High/Low RITH组间差异最显著的特征:\n")
print(head(final_results[order(final_results$T_test_p_value), 
                         c("Feature", "Mean_difference_High_vs_Low", "T_test_p_value", "Significance")], 10))

cat("\n=== 显著相关且显著差异的特征 ===\n")
significant_both <- final_results %>%
  filter(Correlation_p_value < 0.05 & T_test_p_value < 0.05) %>%
  arrange(desc(abs(Correlation_with_HI)))

if (nrow(significant_both) > 0) {
  cat("同时与HI显著相关且在组间有显著差异的特征:\n")
  print(significant_both[, c("Feature", "Correlation_with_HI", "Correlation_p_value", 
                             "T_test_p_value", "Significance")])
} else {
  cat("没有找到同时满足两个条件的特征\n")
}

cat("\n=== 生成的文件列表 ===\n")
cat("01_HI_vs_Morphological_Features_Correlation.pdf - HI与所有形态特征的相关性散点图\n")
cat("02_Morphological_Features_Violin_Plots_with_Stats.pdf - 所有形态特征的小提琴分布图（带统计学分析）\n")
cat("03_Statistical_Analysis_Results.csv - 统计分析结果\n")
cat("04_Violin_Plot_Statistics.csv - 小提琴图详细统计\n")
cat("05_Top_Correlation_Features.pdf - 与HI最相关的特征\n")
cat("06_Top_Difference_Features_with_Stats.pdf - 组间差异最显著的特征（带统计学分析）\n")
cat("07_Significance_Heatmap_Data.csv - 显著性热图数据\n")
cat("08_Analysis_Report.txt - 本分析报告\n")

cat("\n=== 统计学分析说明 ===\n")
cat("1. 相关性分析：Pearson相关系数\n")
cat("2. 组间差异：两样本t检验（不等方差）\n")
cat("3. 显著性标注：使用标准显著性符号（*, **, ***, ns）\n")
cat("4. 小提琴图中的比较线：显示High vs Low RITH的显著性\n")
cat("5. 所有检验均使用双侧检验\n")

cat("\n========================================\n")
cat("          分析完成\n")
cat("========================================\n")

sink()

cat("已保存: 08_Analysis_Report.txt\n")

