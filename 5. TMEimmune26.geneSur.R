# install.packages("survival")
# Remove problematic packages
# remove.packages(c("survminer", "xfun", "ggplot2", "ggpubr"))

# Restart R session (or RStudio)
# Then reinstall in order:
# install.packages("xfun")
# install.packages("ggplot2")
# install.packages("ggpubr")
# install.packages("survminer")

# Load required packages
library(survival)
library(survminer)

inputFile <- "survival.txt"
inputFile2 <- "survival.txt" # Input file name
gene <- "RITH_Phenotype"     # Gene/phenotype name
setwd("/mnt/data4/表型研究2026/ccRCC/RHWU/5.生存分析/Surv")    # Set working directory

# Read input files
rt <- read.table(inputFile, header = T, sep = ",", check.names = F)
rt2 <- read.table(inputFile2, header = T, sep = "\t", check.names = F)

# Group samples by the median value of the gene
# a <- ifelse(rt[, gene] <= 0.5, "low", "high")
# a <- ifelse(rt[, gene] <= median(rt[, gene]), "low", "high")
a <- ifelse(rt$RITH_Phenotype == "Low RITH", "Low RITH", "High RITH")

# Survival difference between high and low groups
diff <- survdiff(Surv(PFI.time, PFI) ~ a, data = rt)
pValue <- 1 - pchisq(diff$chisq, df = 1)

# Calculate number of patients at each time point
fit <- survfit(Surv(PFI.time, PFI) ~ a, data = rt)
summary(fit)

# Draw survival curve
if (pValue < 0.01) {
  pValue <- "p < 0.01"
} else {
  pValue <- paste0("p = ", sprintf("%.03f", pValue))
}

titleName <- gene

surPlot <- ggsurvplot(
  fit,
  data = rt,
  surv.median.line = "hv", # Median survival time
  # conf.int = TRUE,
  # pval = pValue,
  # pval.size = 6,
  # pval.coord = c(12, 3),
  risk.table = T,
  # ncensor.plot = TRUE,
  legend.labs = c("High RITH", "Low RITH"),
  legend.title = "",
  size = 1.6,
  xlab = "Time (month)",
  ylab = "PFI probability",
  font.x = 16,
  font.y = 16,
  break.time.by = 12,
  risk.table.title = "",
  legend = c(0.1, 0.13),
  censor.size = 6,
  main = "Survival curves",
  palette = c("maroon3", "royalblue3"),
  risk.table.height = .25
)

tiff(file = paste0("sur.", gene, ".tiff"), width = 2200, height = 1700, res = 300)
print(surPlot)
dev.off()

diff
