# 加载必要的包
library(limma)
library(edgeR)
library(ggplot2)
library(org.Hs.eg.db)

setwd('/mnt/data4/表型研究2026/ccRCC/RHWU/8.生物学')
# 1. 读取数据
# 假设文件名为 "RHWU_ccRCC_bulk_RNA_Seq_gene_counts.csv"
counts <- read.csv("RHWU_ccRCC_bulk_RNA_Seq_gene_counts.csv", row.names = 1, check.names = FALSE)

# 查看数据维度
dim(counts)  # 应包含约 6 万基因，8 个样本

# 2. 样本分组
# 定义对照组和实验组的样本名称（按列名）
control_samples <- c("王风云11010926")  #,"彭瑞娟11013049" "陶炎兰10965751",
treat_samples   <- c("任志海11029686","付喜明10959018","何兰珍11017067","徐子军10985575","张苹英10977361")

# 检查列名是否都存在
stopifnot(all(control_samples %in% colnames(counts)))
stopifnot(all(treat_samples %in% colnames(counts)))

# 按顺序排列列
samples <- c(control_samples, treat_samples)
counts <- counts[, samples]
group <- factor(c(rep("Control", length(control_samples)),
                  rep("Treatment", length(treat_samples))))

# ------------------------------ 4. 过滤低表达基因 ----------------------------------
keep <- rowSums(cpm(counts) > 1) >= 2
counts_filt <- counts[keep, ]
cat("过滤后基因数:", nrow(counts_filt), "\n")

# ------------------------------ 5. DGEList 与归一化 ----------------------------------
dge <- DGEList(counts = counts_filt, group = group)
dge <- calcNormFactors(dge)

# ------------------------------ 6. voom 转换与线性建模 ----------------------------------
design <- model.matrix(~ group)
v <- voom(dge, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

# ------------------------------ 7. 提取差异表达结果 ----------------------------------
res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
res$ensembl <- rownames(res)

# ------------------------------ 8. ENSG -> Gene Symbol 转换（使用 org.Hs.eg.db） ------------------------------
symbols <- mapIds(org.Hs.eg.db,
                  keys = res$ensembl,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")
res$gene_symbol <- symbols
# 未匹配的保留原始 ENSG
res$gene_symbol[is.na(res$gene_symbol)] <- res$ensembl[is.na(res$gene_symbol)]

# ------------------------------ 9. 标记显著性及上下调 ----------------------------------
logFC_cutoff <- 1
pval_cutoff  <- 0.05

res$regulation <- "Not significant"
res$regulation[res$logFC > logFC_cutoff & res$adj.P.Val < pval_cutoff] <- "Up"
res$regulation[res$logFC < -logFC_cutoff & res$adj.P.Val < pval_cutoff] <- "Down"

cat("上调基因数:", sum(res$regulation == "Up"), "\n")
cat("下调基因数:", sum(res$regulation == "Down"), "\n")

# ------------------------------ 10. 保存完整差异结果 ----------------------------------
write.csv(res, "DEG_results_limma.csv", row.names = FALSE)

# ------------------------------ 11. 火山图：基础版（含辅助线） ----------------------------------
p_base <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Not significant" = "grey")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black") +
  labs(title = "Differential Expression Analysis",
       subtitle = paste0("Treatment vs Control (|log2FC| > ", logFC_cutoff, ", adj.P.Val < ", pval_cutoff, ")"),
       x = "log2 Fold Change", y = "-log10(adjusted P-value)",
       color = "Regulation") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"))

ggsave("volcano_basic.pdf", p_base, width = 8, height = 7)
ggsave("volcano_basic.png", p_base, width = 8, height = 7, dpi = 300)

# ------------------------------ 12. 火山图：标注 top 基因 ----------------------------------
top_up <- res %>%
  filter(regulation == "Up") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

top_down <- res %>%
  filter(regulation == "Down") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 20)

top_genes <- bind_rows(top_up, top_down) %>%
  filter(!is.na(gene_symbol))

library(ggrepel)   # 重要：必须加载


p_labeled <- p_base +
  geom_label_repel(data = top_genes,
                   aes(x = logFC, y = -log10(adj.P.Val), label = gene_symbol),
                   size = 3, max.overlaps = 50, box.padding = 0.5,
                   point.padding = 0.3, segment.color = "grey50")

ggsave("volcano_labeled.pdf", p_labeled, width = 8, height = 8)
ggsave("volcano_labeled.png", p_labeled, width = 8, height = 8, dpi = 300)

# ------------------------------ 13. ECM 基因高亮火山图 ----------------------------------
# 用户提供的 ECM 基因列表
ecm_genes <- c(
  "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL4A3", "COL4A4",
  "COL4A5", "COL4A6", "COL5A1", "COL5A2", "COL5A3", "COL6A1", "COL6A2",
  "COL6A3", "COL7A1", "COL8A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3",
  "COL10A1", "COL11A1", "COL11A2", "COL12A1", "COL14A1", "COL15A1", "COL16A1",
  "COL17A1", "COL18A1", "COL19A1", "COL20A1", "COL21A1", "COL22A1", "COL23A1",
  "COL24A1", "COL25A1", "COL26A1", "COL27A1", "COL28A1",
  "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMB3",
  "LAMB4", "LAMC1", "LAMC2", "LAMC3",
  "ITGA1", "ITGA2", "ITGA2B", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGA7",
  "ITGA8", "ITGA9", "ITGA10", "ITGA11", "ITGAV", "ITGB1", "ITGB2", "ITGB3",
  "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8",
  "DCN", "BGN", "LUM", "VCAN", "ACAN", "NCAN", "BCAN", "HSPG2", "SDC1",
  "SDC2", "SDC3", "SDC4", "GPC1", "GPC2", "GPC3", "GPC4", "GPC5", "GPC6",
  "AGRN", "PRG2", "PRG3", "PRG4",
  "FN1", "FN2", "TNC", "TNN", "TNXB", "THBS1", "THBS2", "THBS3", "THBS4",
  "COMP", "SPARC", "SPARCL1", "SPP1", "BSP", "IBSP", "MEPE", "DMP1",
  "POSTN", "MGP", "MATN1", "MATN2", "MATN3", "MATN4", "VWF", "VTN",
  "MMP1", "MMP2", "MMP3", "MMP7", "MMP8", "MMP9", "MMP10", "MMP11", "MMP12",
  "MMP13", "MMP14", "MMP15", "MMP16", "MMP17", "MMP19", "MMP20", "MMP21",
  "MMP24", "MMP25", "MMP26", "MMP27", "MMP28",
  "TIMP1", "TIMP2", "TIMP3", "TIMP4",
  "ADAMTS1", "ADAMTS2", "ADAMTS3", "ADAMTS4", "ADAMTS5", "ADAMTS6",
  "ADAMTS7", "ADAMTS8", "ADAMTS9", "ADAMTS10", "ADAMTS12", "ADAMTS13",
  "ADAMTS14", "ADAMTS15", "ADAMTS16", "ADAMTS17", "ADAMTS18", "ADAMTS19",
  "ADAMTS20", "ADAM9", "ADAM10", "ADAM12", "ADAM15", "ADAM17",
  "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4",
  "TGFB1", "TGFB2", "TGFB3", "TGFBI", "CTGF", "CYR61", "NOV", "WISP1",
  "WISP2", "WISP3", "FBN1", "FBN2", "FBN3", "MFAP2", "MFAP4", "MFAP5",
  "EMILIN1", "EMILIN2", "EMILIN3", "LTBP1", "LTBP2", "LTBP3", "LTBP4",
  "FBLN1", "FBLN2", "FBLN3", "FBLN4", "FBLN5", "FBLN7", "EFEMP1", "EFEMP2"
)

res$is_ecm <- ifelse(res$gene_symbol %in% ecm_genes, "ECM", "Other")

p_ecm <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(data = subset(res, is_ecm == "Other"),
             aes(color = regulation), alpha = 0.5, size = 1) +
  geom_point(data = subset(res, is_ecm == "ECM"),
             aes(color = regulation), size = 3, shape = 18, alpha = 0.9) +
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Not significant" = "grey")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black") +
  geom_label_repel(data = subset(res, is_ecm == "ECM" & regulation != "Not significant"),
                   aes(label = gene_symbol), size = 3, max.overlaps = 50,
                   box.padding = 0.5, point.padding = 0.3) +
  labs(title = "Volcano Plot Highlighting ECM Genes",
       x = "log2 Fold Change", y = "-log10(adjusted P-value)",
       color = "Regulation") +
  theme_classic()

ggsave("volcano_ECM_highlight.pdf", p_ecm, width = 10, height = 8)
ggsave("volcano_ECM_highlight.png", p_ecm, width = 10, height = 8, dpi = 300)

# ------------------------------ 14. PCA 图（检查分组效果） ----------------------------------
pca_data <- t(v$E)   # 样本×基因
pca <- prcomp(pca_data, scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Group = group)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = "dashed") +
  theme_classic() +
  labs(title = "PCA of Samples",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"))

ggsave("PCA_plot.pdf", pca_plot, width = 6, height = 5)

# ------------------------------ 完成 ----------------------------------
cat("\n分析完成！结果文件已保存至当前目录。\n")
print(table(res$regulation))



# ============================================================================
# GO 和 KEGG 富集分析（基于 clusterProfiler）
# 参考 TCGA 分析代码，适配当前数据
# ============================================================================

# ------------------------------ 加载包 ----------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)      # 用于圆圈图（若不需要可注释）
library(ComplexHeatmap)

# ------------------------------ 设置参数 ----------------------------------
pvalueFilter <- 0.05       # P 值阈值
qvalueFilter <- 0.05       # 校正后 Q 值阈值
logFC_cutoff <- 1          # log2FC 阈值

# ------------------------------ 提取显著差异基因 ----------------------------------
# 根据 adj.P.Val 和 logFC 筛选显著基因
sig_genes <- res %>%
  filter(adj.P.Val < pvalueFilter & abs(logFC) > logFC_cutoff)

cat("显著上调基因数:", sum(sig_genes$logFC > 0), "\n")
cat("显著下调基因数:", sum(sig_genes$logFC < 0), "\n")

# 获取所有显著基因的 gene_symbol 列表
gene_symbols <- unique(sig_genes$gene_symbol)
gene_symbols <- gene_symbols[!is.na(gene_symbols)]  # 去除 NA

# 转换为 ENTREZID
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez <- entrez_ids$ENTREZID

cat("转换后获得 ENTREZID 数量:", length(gene_entrez), "\n")

# ------------------------------ GO 富集分析（BP, CC, MF 分别进行） ----------------------------------
go_all <- enrichGO(gene = gene_entrez,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",                # 同时分析 BP, CC, MF
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE)

# 保存 GO 结果
go_df <- as.data.frame(go_all)
write.csv(go_df, "GO_all.csv", row.names = FALSE)

# 按 ONTOLOGY 分类筛选显著条目
go_sig <- go_df %>% filter(pvalue < pvalueFilter & qvalue < qvalueFilter)

# ------------------------------ KEGG 富集分析 ----------------------------------
kegg_all <- enrichKEGG(gene = gene_entrez,
                       organism = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

kegg_df <- as.data.frame(kegg_all)
write.csv(kegg_df, "KEGG_all.csv", row.names = FALSE)

# ------------------------------ 可视化 1：条形图（GO） ----------------------------------
# 取每个 ONTOLOGY 前 10 个条目
go_top <- go_sig %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = pvalue, n = 10) %>%
  ungroup()

if (nrow(go_top) > 0) {
  p_bar <- ggplot(go_top, aes(x = Count, y = reorder(Description, Count), fill = ONTOLOGY)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = c("BP" = "#E64B35FF", "CC" = "#4DBBD5FF", "MF" = "#00A087FF")) +
    labs(title = "GO Enrichment (Top 10 per Ontology)",
         x = "Gene Count", y = "Term") +
    theme_classic() +
    theme(axis.text.y = element_text(size = 10),
          legend.position = "right")
  ggsave("GO_barplot.pdf", p_bar, width = 8, height = 10)
  print(p_bar)
}

# ------------------------------ 可视化 2：气泡图（GO） ----------------------------------
if (nrow(go_top) > 0) {
  p_dot <- dotplot(go_all, showCategory = 10, split = "ONTOLOGY", label_format = 60) +
    facet_grid(ONTOLOGY ~ ., scale = "free") +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = "GO Enrichment Dotplot") +
    theme(axis.text.y = element_text(size = 8))
  ggsave("GO_dotplot.pdf", p_dot, width = 10, height = 12)
  print(p_dot)
}

# ------------------------------ 可视化 3：条形图（KEGG） ----------------------------------
if (!is.null(kegg_all) && nrow(kegg_all) > 0) {
  kegg_top <- head(kegg_all, 15)
  p_kegg_bar <- ggplot(kegg_top, aes(x = Count, y = reorder(Description, Count))) +
    geom_bar(stat = "identity", fill = "#E64B35FF") +
    labs(title = "KEGG Enrichment (Top 15)",
         x = "Gene Count", y = "Pathway") +
    theme_classic()
  ggsave("KEGG_barplot.pdf", p_kegg_bar, width = 8, height = 6)
  print(p_kegg_bar)
  
  # 气泡图（KEGG）
  p_kegg_dot <- dotplot(kegg_all, showCategory = 15) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = "KEGG Enrichment Dotplot")
  ggsave("KEGG_dotplot.pdf", p_kegg_dot, width = 8, height = 6)
  print(p_kegg_dot)
}



cat("\n富集分析完成！所有结果已保存至当前目录。\n")









# 加载包
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# 设置工作目录
setwd("/mnt/data4/表型研究2026/ccRCC/RHWU/8.生物学")

# 从差异结果中准备排序基因列表（与之前一致）
geneList <- res$logFC
names(geneList) <- res$gene_symbol
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)

# 读取 KEGG GMT 文件（请确保文件存在）
gmt_file <- "c2.cp.kegg.Hs.symbols.gmt"
gmt <- read.gmt(gmt_file)

# 运行 GSEA（对所有通路）
gsea_res <- GSEA(geneList,
                 TERM2GENE = gmt,
                 pvalueCutoff = 0.05,
                 minGSSize = 15,
                 maxGSSize = 500,
                 verbose = FALSE)

# 查看结果
gsea_df <- as.data.frame(gsea_res)
head(gsea_df)

# 保存完整结果
write.csv(gsea_df, "GSEA_KEGG_all.csv", row.names = FALSE)


# 通路名称（KEGG 中通常为 "KEGG_ECM_RECEPTOR_INTERACTION"）
ecm_pathway <- "KEGG_DRUG_METABOLISM_CYTOCHROME_P450"

# 检查通路是否存在
if (ecm_pathway %in% gsea_res@result$ID) {
  
  # 绘制单通路图（带统计表）
  p <- gseaplot2(gsea_res,
                 geneSetID = ecm_pathway,
                 title = "ECM Receptor Interaction",
                 base_size = 12
                 ,pvalue_table = TRUE,  # 显示p值和NES
                 ES_geom = "line"
                 )
  
  # 保存 PDF
  pdf("GSEA_ECM_pathway.pdf", width = 7, height = 6)
  print(p)
  dev.off()
  
  # 显示该通路信息
  cat("\n=== ECM 受体相互作用通路富集结果 ===\n")
  print(gsea_res@result[gsea_res@result$ID == ecm_pathway, 
                        c("ID", "NES", "pvalue", "p.adjust", "setSize")])
  
} else {
  cat("通路", ecm_pathway, "不在 GSEA 结果中。\n")
  cat("可用的通路有（前10个）：\n")
  print(head(gsea_res@result$ID, 10))
}



