# This script do methylation analysis from Illumina EPIC kit. It allow analysis of DMP (differentially
# methylated position) and DMR (differentially methylated region) with CHAMP R package
# (which include DMRCate tool).

###############################################
#                   Library                   #
###############################################

# R library
library(ChAMP)
library(plyr)
library(ggplot2)
library(dplyr)
library(ReMapEnrich)
library(GenomicRanges)

# Load our own R script
source("04.Scripts/volcanoPlot.R")
source("04.Scripts/PCA.R")
source("04.Scripts/Heatmap.R")
source("04.Scripts/GeneOntology.R")
source("04.Scripts/Kegg_pathways.R")


###############################################
#               Load parameters               #
###############################################

set.seed(123)

args <- commandArgs(trailingOnly=TRUE)
pathData <- args[1]
dirResults <- args[2]
path <- paste(dirResults, "/", sep = "")
phenoTableName <- args[3]
conversionTablePath <- args[4]
cutoffFC <- as.numeric(as.character(args[5]))

DEG_test <- args[6]
DEG_cutoffFC <- as.numeric(as.character(args[7]))

ReMAP_db <- args[8]
genePromoter <- args[9]

# pathData <- "00.Data/data/Methylation/Tissue"
# dirResults <- "."
# phenoTableName <- "00.Data/ref/phenoTable.csv"
# conversionTablePath <- "00.Data/ref/Conversion_ID_Name.tab"
# cutoffFC <- 0.2



###############################################
#                  Load data                  #
###############################################
set.seed(123)
# Load raw data
# champ.load() function include some quality control
myLoad <- champ.load(directory = pathData, 
  filterBeads = TRUE, # probes with a beadcount less than 3 will be removed
  detPcut = 0.01, # The detection p-value threshhold. Probes about this cutoff will be filtered out
  filterDetP = TRUE, # probes above the detPcut will be filtered out
  filterMultiHit = TRUE,  # probes in which the probe aligns to multiple locations with bwa as defined in Nordlund et al are removed
  filterSNPs = TRUE, # probes in which the probed CpG falls near a SNP as defined in Nordlund et al are removed
  filterXY = TRUE, # probes from X and Y chromosomes are removed
  arraytype = "EPIC",
  method = "minfi",
  force = TRUE)

# Phenotype & conv table
phenoTable <- read.table(phenoTableName, sep = ",", header=TRUE)
conversionTable <- read.table(conversionTablePath, sep = "\t", header=TRUE)

# Subset phenotype table according to available samples
sampleOrder <- myLoad$pd$Sample_Name
rownames(phenoTable) <- phenoTable$ID
phenoTable <- phenoTable[sampleOrder,]

# Add covariable to analysis
myLoad$pd$Sample_Group <- phenoTable$Phenotype
myLoad$pd$Sex <- phenoTable$Sex
myLoad$pd$Sample_Group <- as.character(myLoad$pd$Sample_Group)


###############################################
#             Raw quality control             #
###############################################

# Quality control
champ.QC(resultsDir= paste(path, "01.QualityControl", sep = ""))


###############################################
#                Normalization                #
###############################################

# Normalization
myNorm <- champ.norm(arraytype="EPIC", method = "BMIQ", cores=4, plotBMIQ=TRUE, resultsDir= paste(path, "02.Normalization", sep = ""))
# Quality control
champ.QC(resultsDir=paste(path, "02.Normalization", sep = ""))

print("Normalization : done")

###############################################
#            Batch effect correction          #
###############################################

# Batch effect correction
myCombat <- champ.runCombat(batchname=c("Slide"))

# # Quality control
champ.QC(resultsDir=paste(path, "03.BatchEffect", sep = ""))
# # Save normalized and corrected table
write.table(myCombat, paste(path, "03.BatchEffect/betaValues.tab", sep = ""), sep = "\t")

print("Batch effect correction : done")

###############################################
#    Differential methylation position test   #
###############################################

myDMP <- champ.DMP(beta = myCombat, 
  pheno = myLoad$pd$Sample_Group, 
  adjPVal = 1, 
  adjust.method = "BH", 
  arraytype = "EPIC")

# Select results
CpG <- myDMP$CTRL_to_sevCCC
CpG$cpgID <- rownames(CpG)

# # Save results
system(paste("mkdir -p ", path, "04.DMP", sep = ""))
write.table(CpG, paste(path, "04.DMP/CpG_test.tab", sep = ""), sep = "\t")

DMPs <- CpG[which(abs(CpG$deltaBeta) >= cutoffFC & CpG$adj.P.Val <= 0.05),]
DMPid <- as.character(DMPs$cpgID)

print("Statistical test : done")


########################################################
#    Differential methylation position visualization   #
########################################################

# PCA
myData <- myCombat[DMPid,]
myPCA(data = myData, path = paste(path, "04.DMP/", sep = ""), phenoTable = phenoTable, sex = TRUE, age = TRUE, name = "CTRL_sevCCC")

# Heatmap
myData <- myCombat[DMPid, ]
myHeatmap(data = myData, path = paste(path, "04.DMP/", sep = ""), phenoTable = phenoTable, name = "CTRL_sevCCC")

print("DMPs visualization : done")

#####################################################
#    Differential methylation position enrichment   #
#####################################################

# Gene ontology
myData <- as.character(DMPs$gene)
myGO(data = myData, path = paste(path, "04.DMP/", sep = ""), name = "CTRL_sevCCC")

# Kegg pathways
myData <- DMPs[,c("gene", "deltaBeta")]
colnames(myData)[1] <- "Name"
myData <- left_join(myData, conversionTable)
myData <- myData[,c("EnsemblID", "deltaBeta")]
myData <- na.omit(myData)
colnames(myData)[2] <- "FC"
myKegg(myData, paste(path, "04.DMP/", sep = ""), name = "CTRL_sevCCC")

print("Gene with DMP enrichment : done")

########################################
#    Differential methylation region   #
########################################

DEG <- read.table(DEG_test, sep = "\t", header = TRUE)
DEG <- DEG[which(abs(DEG$log2FoldChange) >= DEG_cutoffFC & DEG$padj <= 0.05),]
DEG <- as.character(DEG$Name)

DMPs <- DMPs[DMPs$feature %in% c("1stExon", "5'UTR", "TSS200", "TSS1500"),] # 4209 DMPs
DMPs <- DMPs[DMPs$gene %in% DEG,] # 509 DMPs

# In a first time, we create un factor containing the phenotype of each sample
Pheno <- factor(as.character(phenoTable$Phenotype), levels = c("CTRL", "sevCCC"))
# Then, we build the model
design <- model.matrix(~Pheno, coef = 2)

myData <- myCombat[rownames(DMPs),]
dataAnnotation <- cpg.annotate(datatype = "array", myData, what = "Beta", arraytype = "EPIC", analysis.type = "differential", design=design, coef = 2)

myDMR <- dmrcate(dataAnnotation, lambda = 400, C = 2)


### Table with DMRs ###

# 1. We extract the DMRs position
chr <- c()
start <- c()
end <- c()
for(i in 1:length(myDMR@coord)){
  dmr_info <- strsplit(myDMR@coord[i], ":")[[1]]
  chr <- c(chr, dmr_info[1])
  start <- c(start, strsplit(dmr_info[2], "-")[[1]][1])
  end <- c(end, strsplit(dmr_info[2], "-")[[1]][2])
}

# 2. We create the final table, like a bed file
DMRs <- data.frame(chr = chr,
                    start = as.numeric(start),
                    end = as.numeric(end),
                    name = seq(1, length(myDMR@coord), 1),
                    score = rep(0, length(myDMR@coord)),
                    brin = rep(".", length(myDMR@coord)))



myDMP <- DMPs
myDMP$CHR <- paste("chr", myDMP$CHR, sep = "")
myDMP$feature <- as.character(myDMP$feature)
myDMP$feature[myDMP$feature == "TSS200"] <- "TSS"
myDMP$feature[myDMP$feature == "TSS1500"] <- "TSS"


listeGene <- c()
listeFeature <- c()
listeDMR <- c()
for(i in 1:nrow(DMRs)){
  # Récupère les infos par DMR (tous les dmps, leur gene, etc)
  myInfo <- myDMP[which(myDMP$CHR == DMRs$chr[i] & 
                  myDMP$MAPINFO <= DMRs$end[i] & 
                  myDMP$MAPINFO >= DMRs$start[i]),]
  # On recupere les genes uniques
  gene <- unique(as.character(myInfo$gene))
  # On parle avec ',' si différents genes
  if(length(gene) > 1){
    gene <- paste(gene, collapse = ",")
  }
  # Idem pour les features
  feature <- unique(as.character(myInfo$feature))
  if(length(feature) > 1){
    feature <- paste(feature, collapse = ",")
  }
  # On récupère les infos
  listeGene <- c(listeGene, gene)
  listeFeature <- c(listeFeature, feature)
  listeDMR <- c(listeDMR, i)
}

# On les ajoute
DMRs$feature <- listeFeature
DMRs$Gene <- listeGene

system(paste("mkdir -p ", path, "05.DMR", sep = ""))
write.table(DMRs, paste(path, "05.DMR/myDMR_DEG.tab", sep = ""), sep = "\t", row.names=FALSE)
write.table(DMRs[,1:7], paste(path, "05.DMR/myDMR_DEG.bed", sep = ""), sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)

print("DMR : done")


##########################################
# Creation of files for ologram analysis #
##########################################

system(paste("mkdir -p ", path, "06.TF", sep = ""))

promoter <- read.table(genePromoter, sep = "\t", header = FALSE)
colnames(promoter) <- c("chr", "start", "end", "geneName")

promoter$start <- ifelse(promoter$start < 0, 0, promoter$start)

promoterDMR <- promoter[promoter$geneName %in% unique(DMRs$Gene),]
promoterOther <- promoter

write.table(promoterDMR, paste(path, "06.TF/DifferentiallyMethylated_RegRegion.bed", sep = ""), sep ="\t", col.names=FALSE, row.names=FALSE, quote = FALSE)
write.table(promoterOther, paste(path, "06.TF/otherRegRegion.bed", sep = ""), sep ="\t", col.names=FALSE, row.names=FALSE, quote = FALSE)
