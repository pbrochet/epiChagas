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

# Create path for results
path <- paste(dir, "/", sep = "")


#################################################
#           Create heart cell database          #
#################################################

set.seed(42)
doParallel::registerDoParallel(cores = parallel::detectCores())

# We read all the informations about the scRNAseq data analysis
db <- read.table(paste(dir, "/00.Database/GSE109816_normal_heart_umi_matrix.csv", sep = ""),sep = ",", header = TRUE)
info <- read.table(paste(dir, "/00.Database/GSE109816_normal_heart_cell_info.txt", sep = ""), sep = "\t", header = TRUE)
cluster <- read.table(paste(dir, "/00.Database/GSE109816_normal_heart_cell_cluster_info.txt", sep = ""), sep = "\t", header = TRUE)

# We pass the gene name in row names
rownames(db) <- db$X
db$X <- NULL

# We keep only cells from left ventricular
cells <- as.character(info$ID[info$Type %in% c("N_LV_CM", "N_LV_NCM")])
db <- db[,cells]

# We add the cell type information for each cell
db <- data.frame(t(db), check.names=FALSE)
db$ID <- rownames(db)
cluster <- cluster[,c("ID", "CellType")]
db <- left_join(db, cluster)

# We remove cells which cell type is not identified
db <- db[!is.na(db$CellType),]

# We put again gene name in row
myCellType <- as.character(db$CellType)
db$CellType <- db$ID <- NULL
db <- data.frame(t(db), check.names=FALSE)
for(i in 1:ncol(db)){
    colnames(db)[i] <- myCellType[i]
}

db <- db[which(rowSums(db) > 0),]


#################################################
#           Create matrice signature            #
#################################################

data <- log(db+1)

# Step 1: Randomly split the dataset into training and test set
trainTestSet<-splitSCdata(data,numSets = 2,randomize = TRUE)
trainSet<-as.matrix(trainTestSet[[1]])
testSet<-as.matrix(trainTestSet[[2]])
trainSet.30sam <- ADAPTS::scSample(RNAcounts = trainSet, groupSize = 30, randomize = TRUE)
trainSet.3sam <- ADAPTS::scSample(RNAcounts = trainSet, groupSize = 3, randomize = TRUE)

# Step 2: Build a deconvolution seed matrix using ranger forest and estimate the accuracy on pseudo bulk test set
seedMat<-ADAPTS::buildSeed(trainSet=trainSet, trainSet.3sam = trainSet.3sam, trainSet.30sam = trainSet.30sam, genesInSeed = 100, groupSize = 30, randomize = TRUE, num.trees = 1000, plotIt = TRUE)
pseudobulk.test <- data.frame(test=rowSums(testSet))
pseudobulk.test.counts<-table(colnames(testSet))
actFrac.test <- 100 * pseudobulk.test.counts / sum(pseudobulk.test.counts)
estimates.test <- as.data.frame(ADAPTS::estCellPercent.DCQ(seedMat, pseudobulk.test))
colnames(estimates.test)<-'seed'
estimates.test$actFrac<-round(actFrac.test[rownames(estimates.test)],2)
# Calculate correlation coefficient, pvalues etc betwen deconvolution predictions
seedAcc<-ADAPTS::calcAcc(estimates=estimates.test[,1], reference=estimates.test[,2])

# Step 3: Build a deconvolution matrix using all the genes and estimate the accuracy on pseudo bulk test set
allGeneSig <- apply(trainSet.3sam, 1, function(x){tapply(x, colnames(trainSet.3sam), mean, na.rm=TRUE)})
estimates.allGene <- as.data.frame(ADAPTS::estCellPercent.DCQ(t(allGeneSig), pseudobulk.test))
colnames(estimates.allGene)<-'all'
estimates.test<-cbind(estimates.allGene,estimates.test)
allAcc<-ADAPTS::calcAcc(estimates=estimates.test[,1], reference=estimates.test[,3])

# Step 4: Augment the seed matrix and estimate the accuray on pseudo bulk test set
gList <- ADAPTS::gListFromRF(trainSet=trainSet.30sam)
augTrain <- ADAPTS::AugmentSigMatrix(origMatrix = seedMat, fullData = trainSet.3sam, gList = gList, nGenes = 1:100, newData = trainSet.3sam, plotToPDF = FALSE, pdfDir = '.')
estimates.augment <- as.data.frame(ADAPTS::estCellPercent.DCQ(augTrain, pseudobulk.test))
colnames(estimates.augment) <- 'aug'
estimates.test <- cbind(estimates.augment, estimates.test)
augAcc<-ADAPTS::calcAcc(estimates=estimates.test[,1], reference=estimates.test[,4])

# Step 5: Shrink the augmented matrix to minimize the condition number and estimate the accuray on pseudo bulk test set
augTrain.shrink <- ADAPTS::shrinkSigMatrix(augTrain, numChunks=NULL, verbose=FALSE, plotIt=TRUE, sigGenesList=NULL, singleCore=FALSE, fastStop=FALSE)
estimates.shrink <- as.data.frame(ADAPTS::estCellPercent.DCQ(augTrain.shrink, pseudobulk.test))
colnames(estimates.shrink)<-'shrink'
estimates.test <- cbind(estimates.shrink, estimates.test)
titleStr <- paste('Shrunk Signature Matrix,', '# Genes:',nrow(augTrain),'->',nrow(augTrain.shrink))
# pheatmap(augTrain.shrink, main=titleStr)
augshrunkAcc<-ADAPTS::calcAcc(estimates=estimates.test[,1], reference=estimates.test[,5])

# Donc table finale :
heartDB <- augTrain.shrink

myInfo <- data.frame(CellID = c("EC", "MP", "SMC", "FB", "CM"),
                        Name = c("Endothelial cells", "Macrophages", "Smooth muscle cells", "Fibroblasts", "Cardiomyocytes"))

colnames(heartDB) <- c("Cardiomyocytes", "Endothelial cells", "Fibroblasts", "Macrophages", "Smooth muscle cells")

write.table(heartDB, paste(path, "00.Database/heartDB.tab", sep = ""),sep = "\t")
