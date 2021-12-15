###############################################
#                   Library                   #
###############################################

# Test d'expression différentielle
library(DESeq2)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ADAPTS)
library(parallel)
library(ggpubr)


###############################################
#               Load parameters               #
###############################################

# Load arguments
args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
# Path for tables
countTableName <- args[2]
phenoTableName <- args[3]
conversionTablePath <- args[4]
test <- args[5]
cond1 <- strsplit(test, "_")[[1]][1]
cond2 <- strsplit(test, "_")[[1]][2]

# Create path for results
path <- paste(dir, "/", test, "/", sep = "")

# Extract conditions
cond1 <- strsplit(test, "_")[[1]][1]
cond2 <- strsplit(test, "_")[[1]][2]

###############################################
#                  Load data                  #
###############################################

# Phenotype table
phenoTable <- read.table(phenoTableName, sep = ",", header=TRUE)

# Conversion ensembl ID <-> Gene Name
conversionTable <- read.table(conversionTablePath, sep = "\t", header=TRUE)

# Counts
count <- read.table(countTableName, sep = "\t", header = TRUE)

# We use GeneID as row names
rownames(count) <- as.character(count[,1])

# We remove useless columns (for this analysis, like chromosome, start, end...)
count <- count[,-c(1,2,3,4,5,6)]

# We rename the columns
colNames <- colnames(count)
newColNames <- c()
for(i in colNames){
  myCol <- strsplit(i, "\\.")[[1]]
  newColNames <- c(newColNames, myCol[length(myCol) - 4])
}
colnames(count) <- newColNames

# We remove all genes whose expression is null
count <- count[which(rowSums(count)>0),]
count <- count[which(rownames(count) %in% as.character(conversionTable$EnsemblID)),]

# Finally, we select the desire samples, and sort the phenoTable and count table in the same way
phenoTable <- phenoTable[phenoTable$Phenotype %in% c(cond1, cond2),]
count <- count[,colnames(count) %in% as.character(phenoTable$ID)]
rownames(phenoTable) <- phenoTable$ID
phenoTable <- phenoTable[colnames(count),]
phenoTable$Phenotype <- as.character(phenoTable$Phenotype)


###############################################
#               Normalization                 #
###############################################

# Creation of DESeq2 object
deseq <- DESeqDataSetFromMatrix(countData = count, colData = phenoTable, design=~Phenotype)
    # Bien spécifier de ne pas suivre l'ordre alphabétique mais l'ordre des groupes :
deseq$Phenotype <- factor(deseq$Phenotype, levels = c(cond1, cond2))

# Normalization
deseq <- DESeq(deseq)
deseq_sizeF <- estimateSizeFactors(deseq)
deseq_norm <- counts(deseq_sizeF, normalized=TRUE)
write.table(deseq_norm, paste(path, "Normalized_data.tab", sep = ""), sep = "\t")

# Plots comparison before and after normalization
  # Before normalization
pseudoCount <- log2(count + 1)
pseudoCount <- melt(pseudoCount)
colnames(pseudoCount) <- c("ID", "Count")
pseudoCount <- full_join(pseudoCount, phenoTable)
pseudoCount$ID <- factor(pseudoCount$ID, levels = unique(pseudoCount$ID[order(pseudoCount$Phenotype)]))
tiff(paste(path, "Counts_distribution_before_normalization.tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
ggplot(pseudoCount, aes(x = pseudoCount$ID, y = Count, fill = Phenotype)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](Count + 1))) + ggtitle("Boxplot of pseudo-count before normalization")  +
  theme(axis.text.x = element_text(angle = 55, hjust = 1))
dev.off()
  # After normalization
pseudoCount <- log2(as.data.frame(deseq_norm) + 1)
pseudoCount <- melt(pseudoCount)
colnames(pseudoCount) <- c("ID", "Count")
pseudoCount <- full_join(pseudoCount, phenoTable)
pseudoCount$ID <- factor(pseudoCount$ID, levels = unique(pseudoCount$ID[order(pseudoCount$Phenotype)]))
tiff(paste(path, "Counts_distribution_after_normalization.tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
ggplot(pseudoCount, aes(x = pseudoCount$ID, y = Count, fill = Phenotype)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](Count + 1))) + ggtitle("Boxplot of pseudo-count after normalization and batch effect correction") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1))
dev.off()

print("Normalization : done")



#################################################
#          Perform gene type enrichment         #
#################################################

heartDB <- read.table(paste(dir, "/00.Database/heartDB.tab", sep = ""),sep = "\t", header=TRUE)
colnames(heartDB) <- c("Cardiomyocytes", "Endothelial cells", "Fibroblasts", "Macrophages", "Smooth muscle cells")

deseq_norm <- read.table(paste(path, "Normalized_data.tab", sep = ""), sep = "\t", header = TRUE)

mySamples <- colnames(deseq_norm)
deseq_norm$EnsemblID <- rownames(deseq_norm)


# Phenotype table
phenoTable <- read.table(phenoTableName, sep = ",", header=TRUE)
# Conversion ensembl ID <-> Gene Name
conversionTable <- read.table(conversionTablePath, sep = "\t", header=TRUE)



#######################################
# heartDB <- read.table("../00.Database/heartDB.tab", sep = "\t", header = TRUE)
# colnames(heartDB) <- c("Cardiomyocytes", "Endothelial cells", "Fibroblasts", "Macrophages", "Smooth muscle cells")
# deseq_norm <- read.table("Normalized_data.tab", sep = "\t", header = TRUE)
# mySamples <- colnames(deseq_norm)
# deseq_norm$EnsemblID <- rownames(deseq_norm)
# phenoTable <- read.table("../../../../04.Papier_RNAseq_meth/00.Data/ref/phenoTable.csv", sep = ",", header = TRUE)
# conversionTable <- read.table("../../../../04.Papier_RNAseq_meth/00.Data/ref/Conversion_ID_Name.tab", sep = "\t", header = TRUE)
#######################################


# We add gene name to normalized data
data <- left_join(deseq_norm, conversionTable)
data <- na.omit(data)
data <- data[,c("Name", mySamples)]


myAdapts <- function(data, database, name, path, phenoTable){
    # We select only genes available in the database and put them in row names
    myData <- data[data$Name %in% rownames(database),]
    rownames(myData) <- make.unique(as.character(myData$Name))
    myData$Name <- NULL
    # We estimate the percent of cell type by sample
    cellCounts <- estCellPercent(database, myData, method = "DCQ")
    # We remove cell types not present in our data
    res <- data.frame(cellCounts)
    res <- subset(res, rowSums(res) > 0)
    # We plot the result
    res$type <- rownames(res)
    resPlot <- melt(res)
    colnames(resPlot) <- c("CellType", "ID", "Value")
    resPlot <- left_join(resPlot, phenoTable)
    myCols <- as.character(resPlot$Color)
    names(myCols) <- as.character(resPlot$Phenotype)
    resPlot$Phenotype <- factor(resPlot$Phenotype, levels = sort(unique(c(resPlot$Phenotype))))
    plot <- ggplot(resPlot, aes(x = CellType, y = Value, fill = Phenotype)) +
        geom_boxplot() +
        theme_classic() +
        scale_fill_manual("Phenotype", values = myCols) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
        xlab("Cell type") +
        ylab("Percent of cells") 
        # stat_compare_means(aes(group = Phenotype), label = "p.signif", size = 5)

    tiff(paste(path, name, ".tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
    plot(plot)
    dev.off()
    return(cellCounts)
}


#################################################
#        Juste signature tissu cardiaque        #
#################################################

rownames(phenoTable) <- as.character(phenoTable$ID)

# Run analysis
enrichHeart <- myAdapts(data, heartDB, "heartDBonly", path, phenoTable)
enrichHeart <- data.frame(enrichHeart)
enrichHeart <- subset(enrichHeart, rowSums(enrichHeart) > 0)

phenoTable <- phenoTable[colnames(enrichHeart), ]

# Compute some statistics
pval <- c()
for(i in 1:nrow(enrichHeart)){
    test <- wilcox.test(as.numeric(enrichHeart[i,]) ~ as.character(phenoTable$Phenotype))
    pval <- c(pval, test$p.value)
}
enrichHeart$pvalue <- pval
enrichHeart$padj <- p.adjust(enrichHeart$pvalue, method = "BH")

write.table(enrichHeart, paste(path, "HeartDB_enrichment.tab", sep = ""), sep = "\t")

        
#################################################
#        Juste signature immune cells           #
#################################################

LM22 <- ADAPTS::LM22
colnames(LM22) <- c("B cells naive", "B cells memory", "Plasma cells", "T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated", "T cells follicular helper", "T cells regulatory Tregs", "T cells gamma delta", "NK cells resting", "NK cells activated", "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated", "Mast cells resting", "Mast.cells.activated", "Eosinophils", "Neutrophils")

enrichImmune <- myAdapts(data, LM22, "immuneDBonly",  path, phenoTable)
enrichImmune <- data.frame(enrichImmune)
enrichImmune <- subset(enrichImmune, rowSums(enrichImmune) > 0)

pval <- c()
for(i in 1:nrow(enrichImmune)){
    test <- wilcox.test(as.numeric(enrichImmune[i,]) ~ as.character(phenoTable$Phenotype))
    pval <- c(pval, test$p.value)
}

enrichImmune$pvalue <- pval
enrichImmune$padj <- p.adjust(enrichImmune$pvalue, method = "BH")

write.table(enrichImmune, paste(path, "ImmuneDB_enrichment.tab", sep = ""), sep = "\t")

res <- rbind(enrichHeart, enrichImmune)
write.table(res, paste(path, "All_enrichment.tab", sep = ""), sep = "\t")