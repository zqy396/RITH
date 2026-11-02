# Load required packages
library(limma)
library(edgeR)
library(ggplot2)
library(org.Hs.eg.db)

setwd('/mnt/data4/表型研究2026/ccRCC/RHWU/8.生物学')
# 1. Read data
# Assume the file name is "RHWU_ccRCC_bulk_RNA_Seq_gene_counts.csv"
counts <- read.csv("RHWU_ccRCC_bulk_RNA_Seq_gene_counts.csv", row.names = 1, check.names = FALSE)

# Check data dimensions
dim(counts)  # Should contain about 60,000 genes and 8 samples

# 2. Sample grouping
# Define sample names for control and treatment groups (by column names)
control_samples <- c("王风云11010926")  #,"彭瑞娟11013049" "陶炎兰10965751",
treat_samples   <- c("任志海11029686","付喜明10959018","何兰珍11017067","徐子军10985575","张苹英10977361")

# Check whether all column names exist
stopifnot(all(control_samples %in% colnames(counts)))
stopifnot(all(treat_samples %in% colnames(counts)))

# Arrange columns in order
samples <- c(control_samples, treat_samples)
counts <- counts[, samples]
group <- factor(c(rep("Control", length(control_samples)),
                  rep("Treatment", length(treat_samples))))

# ------------------------------ 4. Filter low-expression genes ----------------------------------
keep <- rowSums(cpm(counts) > 1) >= 2
counts_filt <- counts[keep, ]
cat("Number of genes after filtering:", nrow(counts_filt), "\n")

# ------------------------------ 5. DGEList and normalization ----------------------------------
dge <- DGEList(counts = counts_filt, group = group)
dge <- calcNormFactors(dge)

# ------------------------------ 6. voom transformation and linear modeling ----------------------------------
design <- model.matrix(~ group)
v <- voom(dge, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

# ------------------------------ 7. Extract differential expression results ----------------------------------
res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
res$ensembl <- rownames(res)

# ------------------------------ 8. ENSG -> Gene Symbol conversion (using org.Hs.eg.db) ------------------------------
symbols <- mapIds(org.Hs.eg.db,
                  keys = res$ensembl,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")
res$gene_symbol <- symbols
# Keep original ENSG for unmatched entries
res$gene_symbol[is.na(res$gene_symbol)] <- res$ensembl[is.na(res$gene_symbol)]

# ------------------------------ 9. Mark significance and up/downregulation ----------------------------------
logFC_cutoff <- 1
pval_cutoff  <- 0.05

res$regulation <- "Not significant"
res$regulation[res$logFC > logFC_cutoff & res$adj.P.Val < pval_cutoff] <- "Up"
res$regulation[res$logFC < -logFC_cutoff & res$adj.P.Val < pval_cutoff] <- "Down"

cat("Number of upregulated genes:", sum(res$regulation == "Up"), "\n")
cat("Number of downregulated genes:", sum(res$regulation == "Down"), "\n")

# ------------------------------ 10. Save complete differential expression results ----------------------------------
write.csv(res, "DEG_results_limma.csv", row.names = FALSE)

# ------------------------------ 11. Volcano plot: basic version (with reference lines) ----------------------------------
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

# ------------------------------ 12. Volcano plot: label top genes ----------------------------------
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

library(ggrepel)   # Important: must load

p_labeled <- p_base +
  geom_label_repel(data = top_genes,
                   aes(x = logFC, y = -log10(adj.P.Val), label = gene_symbol),
                   size = 3, max.overlaps = 50, box.padding = 0.5,
                   point.padding = 0.3, segment.color = "grey50")

ggsave("volcano_labeled.pdf", p_labeled, width = 8, height = 8)
ggsave("volcano_labeled.png", p_labeled, width = 8, height = 8, dpi = 300)

# ------------------------------ 13. Volcano plot highlighting ECM genes ----------------------------------
# User-provided ECM gene list
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

# ------------------------------ 14. PCA plot (to check grouping) ----------------------------------
pca_data <- t(v$E)   # samples x genes
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

# ------------------------------ Completed ----------------------------------
cat("\nAnalysis completed! Result files saved to the current directory.\n")
print(table(res$regulation))

# ============================================================================
# GO and KEGG enrichment analysis (based on clusterProfiler)
# Adapted from TCGA analysis code for the current data
# ============================================================================

# ------------------------------ Load packages ----------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

# ------------------------------ Set parameters ----------------------------------
pvalueFilter <- 0.05
qvalueFilter <- 0.05
logFC_cutoff <- 1

# ------------------------------ Extract significant differentially expressed genes ----------------------------------
sig_genes <- res %>%
  filter(adj.P.Val < pvalueFilter & abs(logFC) > logFC_cutoff)

cat("Number of significantly upregulated genes:", sum(sig_genes$logFC > 0), "\n")
cat("Number of significantly downregulated genes:", sum(sig_genes$logFC < 0), "\n")

# Get gene_symbol list of all significant genes
gene_symbols <- unique(sig_genes$gene_symbol)
gene_symbols <- gene_symbols[!is.na(gene_symbols)]

# Convert to ENTREZID
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez <- entrez_ids$ENTREZID

cat("Number of ENTREZIDs obtained after conversion:", length(gene_entrez), "\n")

# ------------------------------ GO enrichment analysis (BP, CC, MF separately) ----------------------------------
go_all <- enrichGO(gene = gene_entrez,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE)

# Save GO results
go_df <- as.data.frame(go_all)
write.csv(go_df, "GO_all.csv", row.names = FALSE)

# Filter significant terms by ONTOLOGY
go_sig <- go_df %>% filter(pvalue < pvalueFilter & qvalue < qvalueFilter)

# ------------------------------ KEGG enrichment analysis ----------------------------------
kegg_all <- enrichKEGG(gene = gene_entrez,
                       organism = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

kegg_df <- as.data.frame(kegg_all)
write.csv(kegg_df, "KEGG_all.csv", row.names = FALSE)

# ------------------------------ Visualization 1: bar plot (GO) ----------------------------------
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

# ------------------------------ Visualization 2: dot plot (GO) ----------------------------------
if (nrow(go_top) > 0) {
  p_dot <- dotplot(go_all, showCategory = 10, split = "ONTOLOGY", label_format = 60) +
    facet_grid(ONTOLOGY ~ ., scale = "free") +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = "GO Enrichment Dotplot") +
    theme(axis.text.y = element_text(size = 8))
  ggsave("GO_dotplot.pdf", p_dot, width = 10, height = 12)
  print(p_dot)
}

# ------------------------------ Visualization 3: bar plot (KEGG) ----------------------------------
if (!is.null(kegg_all) && nrow(kegg_all) > 0) {
  kegg_top <- head(kegg_all, 15)
  p_kegg_bar <- ggplot(kegg_top, aes(x = Count, y = reorder(Description, Count))) +
    geom_bar(stat = "identity", fill = "#E64B35FF") +
    labs(title = "KEGG Enrichment (Top 15)",
         x = "Gene Count", y = "Pathway") +
    theme_classic()
  ggsave("KEGG_barplot.pdf", p_kegg_bar, width = 8, height = 6)
  print(p_kegg_bar)

  # Dot plot (KEGG)
  p_kegg_dot <- dotplot(kegg_all, showCategory = 15) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = "KEGG Enrichment Dotplot")
  ggsave("KEGG_dotplot.pdf", p_kegg_dot, width = 8, height = 6)
  print(p_kegg_dot)
}

cat("\nEnrichment analysis completed! All results saved to the current directory.\n")

# ------------------------------ GSEA ----------------------------------

# Load packages
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Set working directory
setwd("/mnt/data4/表型研究2026/ccRCC/RHWU/8.生物学")

# Prepare ranked gene list from differential results (same as before)
geneList <- res$logFC
names(geneList) <- res$gene_symbol
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)

# Read KEGG GMT file (ensure the file exists)
gmt_file <- "c2.cp.kegg.Hs.symbols.gmt"
gmt <- read.gmt(gmt_file)

# Run GSEA (for all pathways)
gsea_res <- GSEA(geneList,
                 TERM2GENE = gmt,
                 pvalueCutoff = 0.05,
                 minGSSize = 15,
                 maxGSSize = 500,
                 verbose = FALSE)

# View results
gsea_df <- as.data.frame(gsea_res)
head(gsea_df)

# Save complete results
write.csv(gsea_df, "GSEA_KEGG_all.csv", row.names = FALSE)

# Pathway name (usually "KEGG_ECM_RECEPTOR_INTERACTION" in KEGG)
ecm_pathway <- "KEGG_STARCH_AND_SUCROSE_METABOLISM"

# Check whether the pathway exists
if (ecm_pathway %in% gsea_res@result$ID) {

  # Draw single-pathway plot (with statistics table)
  p <- gseaplot2(gsea_res,
                 geneSetID = ecm_pathway,
                 title = "ECM Receptor Interaction",
                 base_size = 12
                 ,pvalue_table = TRUE,
                 ES_geom = "line"
                 )

  # Save PDF
  pdf("GSEA_ECM_pathway.pdf", width = 11, height = 6)
  print(p)
  dev.off()

  # Display information for this pathway
  cat("\n=== Enrichment results for ECM receptor interaction pathway ===\n")
  print(gsea_res@result[gsea_res@result$ID == ecm_pathway,
                        c("ID", "NES", "pvalue", "p.adjust", "setSize")])

} else {
  cat("Pathway", ecm_pathway, "is not in the GSEA results.\n")
  cat("Available pathways (top 10):\n")
  print(head(gsea_res@result$ID, 10))
}
