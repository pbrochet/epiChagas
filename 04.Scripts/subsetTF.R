
###############################################
#               Load parameters               #
###############################################

args <- commandArgs(trailingOnly=TRUE)
TF_Genes <- args[1]
myTFs <- args[2]
myDMRs <- args[3]
path <- args[4]


# TF_Genes <- "06.TF/TFs_Genes.tab"
# myTFs <- "06.TF/TFs_list.txt"
# myDMRs <- "06.TF/DifferentiallyMethylated_RegRegion.bed"


###############################################
#                   Load data                 #
###############################################

data <- read.table(TF_Genes, sep = "\t", header = FALSE)
data <- data[,c("V4", "V8")]
data <- unique(data)
colnames(data) <- c("TF", "Gene")

myTFs <- read.table(myTFs, sep = "\t", header = FALSE)
myTFs <- as.character(myTFs$V1)

myDMRs <- read.table(myDMRs, sep = "\t", header = FALSE)
totGene <- unique(as.character(myDMRs$V4))


###############################################
#                Perform subset               #
###############################################

data <- data[data$TF %in% myTFs,]

nbCible <- c()
for(i in unique(as.character(data$TF))){
    nbCible <- c(nbCible, nrow(data[data$TF == i,]))
}

res <- data.frame(TF = unique(as.character(data$TF)), nbCible = nbCible)

totGene <- unique(as.character(myDMRs$V4))

res <- res[res$nbCible >= round(length(totGene) * 0.5),]

write.table(as.character(res$TF), paste(path, "/listeTF_subset.txt", sep = ""), sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)


# rm -r 05.Results/03.Tissue_methylation_analysis/Test3/06.TF/myTFs/
# rm -r 05.Results/03.Tissue_methylation_analysis/Test3/06.TF/Ologram/
# rm -r 05.Results/03.Tissue_methylation_analysis/Test3/06.TF/listeTF_subset.txt 