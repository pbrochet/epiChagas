###############################################
#                   Library                   #
###############################################

library(ggplot2)
library(reshape2)
library(dplyr)

# Load our own R script
source("04.Scripts/volcanoPlot.R")
source("04.Scripts/PCA.R")
source("04.Scripts/Heatmap.R")


###############################################
#               Load parameters               #
###############################################

args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
path <- paste(dir, "/", sep = "")
normDataName <- args[2]
testDataName <- args[3]
phenoTableName <- args[4]
cutoffFC <- as.numeric(as.character(args[5]))
lncRNA2TargetName <- args[6]
lncTarDName <- args[7]
lncRNAdiseaseName <- args[8]

###############################################
#                  Load data                  #
###############################################

# Normalized data
normData <- read.table(normDataName, sep = "\t", header = TRUE)
# Test data
testData <- read.table(testDataName, sep = "\t", header = TRUE)
# Phenotype table
phenoTable <- read.table(phenoTableName, sep = ",", header=TRUE)
# Databases
lncRNA2Target <- read.table(lncRNA2TargetName, sep = "\t", header = TRUE, fill = TRUE)
lncTarD <- read.table(lncTarDName, sep = "\t", header = TRUE, fill = TRUE)
lncRNAdisease <- read.table(lncRNAdiseaseName, sep = "\t", header = TRUE, fill = TRUE)


###############################################
#            DE lncRNAs visualization         #
###############################################

# We keep only ncRNAs
testData <- testData[testData$Type %in% c("lincRNA", "miRNA", "snRNA", "snoRNA", "3prime_overlapping_ncrna"),]
DEGs <- testData[which(abs(testData$log2FoldChange) >= cutoffFC  & testData$padj <= 0.05),]

# VolcanoPlot
myData <- testData[,c("EnsemblID", "log2FoldChange", "padj")]
colnames(myData) <- c("ID", "FC", "padj")
myVolcanoPlot(data = myData, path = path, xlab = "log2(FC)", ylab = "-log10(adjusted pvalue)", cutoffFC = cutoffFC, name = "ncRNA")

# PCA
myData <- normData[as.character(DEGs$EnsemblID), ]
myPCA(data = myData, path = path, phenoTable = phenoTable, sex = TRUE, age = TRUE, name = "ncRNA")

# Heatmap
myData <- normData[as.character(DEGs$EnsemblID), ]
myHeatmap(data = myData, path = path, phenoTable = phenoTable, name = "ncRNA")


###############################################
#              ncRNAs enrichment              #
###############################################

# 1. lncRNA2Target
lncRNA2Target <- lncRNA2Target[which(lncRNA2Target$lncRNA_name_from_paper %in% as.character(DEGs$Name) | lncRNA2Target$LncRNA_official_symbol %in% as.character(DEGs$Name)),]
write.table(lncRNA2Target, paste(path, "lncRNA2Target_enrichment.tab", sep = ""), sep = "\t", row.names=FALSE)

# 2. lncTarD
lncTarD <- lncTarD[lncTarD$Regulator %in% as.character(DEGs$Name),]
write.table(lncTarD, paste(path, "lncTarD_enrichment.tab", sep = ""), sep = "\t", row.names=FALSE)

# 3. lncRNAdisease
lncRNAdisease <- lncRNAdisease[lncRNAdisease$ncRNA.Symbol %in% as.character(DEGs$Name),]
write.table(lncRNAdisease, paste(path, "lncRNAdisease_enrichment.tab", sep = ""), sep = "\t", row.names=FALSE)
