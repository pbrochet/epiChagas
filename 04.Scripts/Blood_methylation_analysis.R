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

# Load our own R script
source("04.Scripts/volcanoPlot.R")
source("04.Scripts/PCA.R")
source("04.Scripts/Heatmap.R")
source("04.Scripts/GeneOntology.R")
source("04.Scripts/Kegg_pathways.R")

set.seed(123)



###############################################
#               Load parameters               #
###############################################

args <- commandArgs(trailingOnly=TRUE)
pathData <- args[1]
dirResults <- args[2]
path <- paste(dirResults, "/", sep = "")
phenoTableName <- args[3]
conversionTablePath <- args[4]

# pathData <- "00.Data/data/Methylation/Blood"
# dirResults <- "."
# phenoTableName <- "00.Data/ref/phenoTable.csv"
# conversionTablePath <- "00.Data/ref/Conversion_ID_Name.tab"

###############################################
#                  Load data                  #
###############################################

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

# Phenotype table
phenoTable <- read.table(phenoTableName, sep = ",", header=TRUE)
conversionTable <- read.table(conversionTablePath, sep = "\t", header=TRUE)

# Subset phenotype table according to available samples
sampleOrder <- myLoad$pd$Sample_Name
rownames(phenoTable) <- phenoTable$ID
phenoTable <- phenoTable[sampleOrder,]

# Add covariable to analysis
myLoad$pd$Sample_Group <- as.character(phenoTable$Phenotype)
myLoad$pd$Sex <- phenoTable$Sex
myLoad$pd$Sample_Group <- as.character(myLoad$pd$Sample_Group)
myLoad$pd$Form <- ifelse(myLoad$pd$Sample_Group %in% "CTRL", "CTRL", "CCC")

for(i in names(myLoad$pd)){
  myLoad$pd[[i]] <- as.character(myLoad$pd[[i]])
}

###############################################
#             Raw quality control             #
###############################################

# Quality control
champ.QC(resultsDir= paste(path, "01.QualityControl", sep = ""))


###############################################
#                Normalization                #
###############################################

# Normalization
myNorm <- champ.norm(arraytype="EPIC", method = "BMIQ", cores=4, plotBMIQ=TRUE, resultsDir= paste(path, "/02.Normalization", sep = ""))
# Quality control
champ.QC(resultsDir=paste(path, "02.Normalization", sep = ""))


###############################################
#            Batch effect correction          #
###############################################

# Batch effect correction
myCombat <- champ.runCombat(beta = myNorm, batchname=c("Slide", "Array", "Sex"))

# Quality control
champ.QC(resultsDir=paste(path, "03.BatchEffect", sep = ""))
# Save normalized and corrected table
write.table(myCombat, paste(path, "03.BatchEffect/betaValues.tab", sep = ""), sep = "\t")


###############################################
#    Differential methylation position test   #
###############################################

system(paste("mkdir -p ", path, "04.DMP", sep = ""))
myPath <- paste(path, "04.DMP/", sep = "")

######################
### CTRL vs sevCCC ###
######################

# DMP test
myDMP <- champ.DMP(beta = myCombat, 
  pheno = myLoad$pd$Sample_Group, 
  adjPVal = 1, 
  adjust.method = "BH", 
  arraytype = "EPIC",
  compare.group = c("CTRL", "sevCCC"))

# Select results
CpG <- myDMP$CTRL_to_sevCCC
CpG$CpGID <- rownames(CpG)


# Save results
write.table(CpG, paste(myPath, "CTRL_sevCCC.tab", sep = ""), sep = "\t")
DMPs <- CpG[which(CpG$adj.P.Val <= 0.05),]
DMPid <- as.character(DMPs$CpGID)
DMPid_test1 <- DMPid

# Differential methylation position visualization

# PCA
phenoData <- phenoTable[phenoTable$Phenotype %in% c("CTRL", "sevCCC"),]
myData <- myCombat[DMPid,as.character(phenoData$ID)]
myPCA(data = myData, path = myPath, phenoTable = phenoData, sex = TRUE, age = FALSE, name = "CTRL_sevCCC")

# Heatmap
myHeatmap(data = myData, path = myPath, phenoTable = phenoData, name = "CTRL_sevCCC")


# Differential methylation position enrichment

# Gene ontology
myData <- as.character(DMPs$gene)
myGO(data = myData, path = myPath, name = "CTRL_sevCCC")

# Kegg pathways
myData <- DMPs[,c("gene", "deltaBeta")]
colnames(myData)[1] <- "Name"
myData <- left_join(myData, conversionTable)
myData <- myData[,c("EnsemblID", "deltaBeta")]
myData <- na.omit(myData)
colnames(myData)[2] <- "FC"
myKegg(myData, myPath, name = "CTRL_sevCCC")



########################
### modCCC vs sevCCC ###
########################

system(paste("mkdir -p ", path, "04.DMP", sep = ""))
myPath <- paste(path, "04.DMP/", sep = "")

# DMP test
myDMP <- champ.DMP(beta = myCombat, 
  pheno = myLoad$pd$Sample_Group, 
  adjPVal = 1, 
  adjust.method = "BH", 
  arraytype = "EPIC",
  compare.group = c("modCCC", "sevCCC"))

# Select results
CpG <- myDMP$modCCC_to_sevCCC
CpG$CpGID <- rownames(CpG)

# Save results
write.table(CpG, paste(myPath, "modCCC_sevCCC.tab", sep = ""), sep = "\t")
DMPs <- CpG[which(CpG$adj.P.Val <= 0.05),]
DMPid <- as.character(DMPs$CpGID)


# Differential methylation position visualization

# PCA
phenoData <- phenoTable[phenoTable$Phenotype %in% c("modCCC", "sevCCC"),]
myData <- myCombat[DMPid,as.character(phenoData$ID)]
myPCA(data = myData, path = myPath, phenoTable = phenoData, sex = TRUE, age = FALSE, name = "modCCC_sevCCC")

# Heatmap
myHeatmap(data = myData, path = myPath, phenoTable = phenoData, name = "modCCC_sevCCC")

# PCA 3 groups
myData <- myCombat[unique(c(DMPid, DMPid_test1)),]
myPCA(data = myData, path = myPath, phenoTable = phenoTable, sex = TRUE, age = FALSE, name = "CTRL_modCCC_sevCCC")


# Differential methylation position enrichment

# Gene ontology
myData <- as.character(DMPs$gene)
myGO(data = myData, path = myPath, name = "modCCC_sevCCC")

# Kegg pathways
myData <- DMPs[,c("gene", "deltaBeta")]
colnames(myData)[1] <- "Name"
myData <- left_join(myData, conversionTable)
myData <- myData[,c("EnsemblID", "deltaBeta")]
myData <- na.omit(myData)
colnames(myData)[2] <- "FC"
myKegg(myData, myPath, name = "modCCC_sevCCC")



# # 5. Biomarkers 
# #############################################

# pheno2 <- read.table("newPhenoTable.tab", sep = "\t", header = TRUE)

# myLoad$pd$Test1 <- pheno2$Test1
# myLoad$pd$Test1 <- as.character(myLoad$pd$Test1)
# myLoad$pd$Test2 <- pheno2$Test2
# myLoad$pd$Test2 <- as.character(myLoad$pd$Test2)

# # DMP test
# myDMP <- champ.DMP(beta = myCombat, 
#   pheno = myLoad$pd$Test1, 
#   adjPVal = 1, 
#   adjust.method = "BH", 
#   arraytype = "EPIC",
#   compare.group = c("CTRLtrain", "CCCtrain"))
# # Select results
# CpG <- myDMP$CCCtrain_to_CTRLtrain

# CpG$deltaBeta <- -CpG$deltaBeta
# # Create a column containing CpG ID
# CpG$CpGID <- rownames(CpG)


# # Save results
# write.table(CpG, "04.DMP/CTRL_CCC_train.tab", sep = "\t")

# # DMP test
# myDMP <- champ.DMP(beta = myCombat, 
#   pheno = myLoad$pd$Test2, 
#   adjPVal = 1, 
#   adjust.method = "BH", 
#   arraytype = "EPIC",
#   compare.group = c("modCCCtrain", "sevCCCtrain"))

# # Select results
# CpG <- myDMP$modCCCtrain_to_sevCCCtrain

# # Create a column containing CpG ID
# CpG$CpGID <- rownames(CpG)


# # Save results
# write.table(CpG, "04.DMP/modCCC_sevCCC_train.tab", sep = "\t")

# write.table(phenoTable, "newPhenoTable.tab", sep = "\t")



