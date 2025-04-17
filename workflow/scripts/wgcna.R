# 导入数据
library(tidyverse)
library(WGCNA)
dataA <- read.csv("age.csv")
dataA <- column_to_rownames(dataA, var = "Name")

dataY <- read.csv("year.csv")
dataY <- column_to_rownames(dataY, var = "Name")

## 因为WGCNA针对的是基因进行聚类，这个时候需要转置
# 过滤掉平均表达量小于1的基因
datExprA <- t(dataA[apply(dataA, 1, mean) >= 1, ])


powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft <- pickSoftThreshold(datExprA, powerVector = powers, verbose = 0)

png("step2-beta-value.png", width = 800, height = 600)
# Plot the results:
## sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
  main = paste("Scale independence")
)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers, cex = cex1, col = "red"
)
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.90, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
  main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# 构建网络
net <- blockwiseModules(
  datExprA,
  power = 10,
  maxBlockSize = 100,
  TOMType = "unsigned", minModuleSize = 5,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F,
  verbose = 3
)

# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
# table(mergedColors)
moduleColors <- mergedColors
# Plot the dendrogram and the module colors underneath
png("step4-genes-modules.png", width = 800, height = 600)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()
## assign all of the gene to their corresponding module
## hclust for the genes.

dataAout <- as.data.frame(net$colors)

needed <- c("LK_T_001667_c00_g02_i02.p1", "LK_I_c16100_79332", "LK_I_c41088_26002")

dataAout$color <- moduleColors

dataAout[needed, ]

write.csv(dataAout[needed, ], file = "step4-genes-modules.csv")

# 提取模块3的基因

# Recalculate topological overlap
TOM <- TOMsimilarityFromExpr(datExprA, power = 10)
# Select module
module <- "turquoise"
# Select module probes
probes <- colnames(datExprA) ## 我们例子里面的probe就是基因
inModule <- (moduleColors == module)
modProbes <- probes[inModule]
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)
## 模块对应的基因关系矩
cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse = "-"), ".txt", sep = ""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse = "-"), ".txt", sep = ""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes,
  nodeAttr = moduleColors[inModule]
)
