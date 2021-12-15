###############################################
#                   Library                   #
###############################################

# Test d'expression différentielle
library( DESeq2)
library( ggplot2)
library( reshape2)
library( dplyr)

# Load our own R script
source( "04.Scripts/volcanoPlot.R")
source( "04.Scripts/PCA.R")
source( "04.Scripts/Heatmap.R")
source( "04.Scripts/GeneOntology.R")
source( "04.Scripts/Kegg_pathways.R")


###############################################
#               Load parameters               #
###############################################

args <- commandArgs( trailingOnly=TRUE)
dir <- args[1]
path <- paste( dir, "/04.Differential_expression/", sep = "")
countTableName <- args[2]
phenoTableName <- args[3]
conversionTablePath <- args[4]
cutoffFC <- as.numeric( as.character( args[5]))

test <- args[6]
cond1 <- strsplit( test, "_")[[1]][1]
cond2 <- strsplit( test, "_")[[1]][2]

path <- paste( path, test, "/", sep = "")

print( cutoffFC)
print( class( cutoffFC))

# path <- "./"
# cutoffFC <- 1.5

###############################################
#                  Load data                  #
###############################################

# Phenotype table
phenoTable <- read.table( phenoTableName, sep = ",", header=TRUE)

# Conversion ensembl ID <-> Gene Name
conversionTable <- read.table( conversionTablePath, sep = "\t", header=TRUE)

# Counts
count <- read.table( countTableName, sep = "\t", header = TRUE)

# We use GeneID as row names
rownames( count) <- as.character( count[,1])

# We remove useless columns ( for this analysis, like chromosome, start, end...)
count <- count[,-c( 1,2,3,4,5,6)]

# We rename the columns
colNames <- colnames( count)
newColNames <- c( )
for( i in colNames){
  myCol <- strsplit( i, "\\.")[[1]]
  newColNames <- c( newColNames, myCol[length( myCol) - 4])
}
colnames( count) <- newColNames

# We remove all genes whose expression is null
count <- count[which( rowSums( count)>0),]

# We remove noising genes ( ex. Ig)
remove <- c( "IG_V_gene", "IG_J_gene", "IG_C_gene", "IG_V_pseudogene", "IG_C_pseudogene", "IG_D_gene", "TR_C_gene", "TR_V_gene", "TR_J_pseudogene", "TR_J_gene", "TR_V_pseudogene", "TR_D_gene", "IG_J_pseudogene")
conversionTable <- conversionTable[!conversionTable$Type %in% remove,] # 62.831
count <- count[which( rownames( count) %in% as.character( conversionTable$EnsemblID)),] # 45.889

# Finally, we select the desire samples, and sort the phenoTable and count table in the same way
phenoTable <- phenoTable[phenoTable$Phenotype %in% c( cond1, cond2),]
count <- count[,colnames( count) %in% as.character( phenoTable$ID)]
rownames( phenoTable) <- phenoTable$ID
phenoTable <- phenoTable[colnames( count),]

phenoTable$Phenotype <- as.character( phenoTable$Phenotype)

print( phenoTable)

###############################################
#               Normalization                 #
###############################################

# Creation of DESeq2 object
deseq <- DESeqDataSetFromMatrix( countData = count, colData = phenoTable, design=~Phenotype)
    # Bien spécifier de ne pas suivre l'ordre alphabétique mais l'ordre des groupes :
deseq$Phenotype <- factor( deseq$Phenotype, levels = c( cond1, cond2))

# Normalization
deseq <- DESeq( deseq)
deseq_sizeF <- estimateSizeFactors( deseq)
deseq_norm <- counts( deseq_sizeF, normalized=TRUE)
write.table( deseq_norm, paste( path, "Normalized_data.tab", sep = ""), sep = "\t")

# Plots comparison before and after normalization
  # Before normalization
pseudoCount <- log2( count + 1)
pseudoCount <- melt( pseudoCount)
colnames( pseudoCount) <- c( "ID", "Count")
pseudoCount <- full_join( pseudoCount, phenoTable)
pseudoCount$ID <- factor( pseudoCount$ID, levels = unique( pseudoCount$ID[order( pseudoCount$Phenotype)]))
tiff( paste( path, "Counts_distribution_before_normalization.tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
ggplot( pseudoCount, aes( x = pseudoCount$ID, y = Count, fill = Phenotype)) + geom_boxplot( ) + xlab( "") +
  ylab( expression( log[2]( Count + 1))) + ggtitle( "Boxplot of pseudo-count before normalization")  +
  theme( axis.text.x = element_text( angle = 55, hjust = 1))
dev.off( )
  # After normalization
pseudoCount <- log2( as.data.frame( deseq_norm) + 1)
pseudoCount <- melt( pseudoCount)
colnames( pseudoCount) <- c( "ID", "Count")
pseudoCount <- full_join( pseudoCount, phenoTable)
pseudoCount$ID <- factor( pseudoCount$ID, levels = unique( pseudoCount$ID[order( pseudoCount$Phenotype)]))
tiff( paste( path, "Counts_distribution_after_normalization.tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
ggplot( pseudoCount, aes( x = pseudoCount$ID, y = Count, fill = Phenotype)) + geom_boxplot( ) + xlab( "") +
  ylab( expression( log[2]( Count + 1))) + ggtitle( "Boxplot of pseudo-count after normalization and batch effect correction") +
  theme( axis.text.x = element_text( angle = 55, hjust = 1))
dev.off( )

print( "Normalization : done")

###############################################
#       Differential expression test          #
###############################################

# Without shrinkage
testDEG <- results( deseq, contrast = c( "Phenotype", cond2, cond1))
tiff( paste( path, "Mean_normalized_counts_without_shrinkage.tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
DESeq2::plotMA( testDEG)
dev.off( )

# With shrinkage
testDEG <- lfcShrink( deseq, coef=paste( "Phenotype_", cond2, "_vs_", cond1, sep = ""), type = "apeglm")
tiff( paste( path, "Mean_normalized_counts_with_shrinkage.tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
DESeq2::plotMA( testDEG)
dev.off( )

# We add missing informations ( gene name, type...)
testDEG$EnsemblID <- rownames( testDEG)
testDEG <- as.data.frame( testDEG)
testDEG <- inner_join( testDEG, conversionTable)

# We keep the significative genes
DEGs <- testDEG[which( abs( testDEG$log2FoldChange) >= cutoffFC  & testDEG$padj <= 0.05),]
DEGs <- as.data.frame( DEGs)
write.table( testDEG, paste( path, "DEGs_test.tab", sep = ""), sep = "\t", row.names=FALSE)

# We summarize information about the number of genes by type
nGeneType <- data.frame( table( as.character( testDEG$Type)))
colnames( nGeneType) <- c( "Type", "nGene")
nGeneType$percentGene <- round( nGeneType$nGene * 100 / sum(  nGeneType$nGene),2)

nDEGType <- data.frame( table( as.character( DEGs$Type)))
colnames( nDEGType) <- c( "Type", "nDEG")
nDEGType$percentDEG <- round( nDEGType$nDEG * 100 / sum(  nDEGType$nDEG),2)

syntheseNGene <- full_join( nGeneType, nDEGType)
syntheseNGene[is.na( syntheseNGene)] <- 0

# And compute the enrichment in each gene type
chisqPvalue <- c()
chisqOR <- c()

for( i in 1:nrow( syntheseNGene)){
  contigencyTable <- data.frame( raw = c( syntheseNGene$nGene[i], sum( syntheseNGene$nGene)), deg = c( syntheseNGene$nDEG[i], sum( syntheseNGene$nDEG)))
  chisqPvalue <- c( chisqPvalue,chisq.test( contigencyTable)$p.value)
  chisqOR <- c( chisqOR,chisq.test( contigencyTable)$statistic)
}
syntheseNGene$chisqPvalue = signif( chisqPvalue,2)
syntheseNGene$chisqOR = signif( chisqOR,2)

write.table( syntheseNGene, paste( path, "GeneType_distribution.tab", sep = ""), sep = "\t", row.names=FALSE)


###############################################
#                 DEGs visualization          #
###############################################

# VolcanoPlot
myData <- testDEG[,c( "EnsemblID", "log2FoldChange", "padj")]
colnames( myData) <- c( "ID", "FC", "padj")
myVolcanoPlot( data = myData, path = path, xlab = "log2( FC)", ylab = "-log10( adjusted pvalue)", cutoffFC = cutoffFC, name = test)

# PCA
myData <- deseq_norm[as.character( DEGs$EnsemblID), ]
myPCA( data = myData, path = path, phenoTable = phenoTable, sex = TRUE, age = TRUE, name = test)

# Heatmap
myData <- deseq_norm[as.character( DEGs$EnsemblID), ]
myHeatmap( data = myData, path = path, phenoTable = phenoTable, name = test)


###############################################
#                 DEGs enrichment             #
###############################################

# Gene ontology
myData <- as.character( DEGs$Name)
myGO( data = myData, path = path, name = test)

# Kegg pathways
myData <- DEGs[,c( "EnsemblID", "log2FoldChange")]
colnames( myData)[2] <- "FC"
myKegg( myData, path, name = test)

