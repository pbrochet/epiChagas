###############################################
#               Load parameters               #
###############################################

args <- commandArgs(trailingOnly=TRUE)
ologramOutput <- args[1]
dirResults <- args[2]
path <- paste(dirResults, "/", sep = "")
complexSize <- as.numeric(as.character(args[3]))


#########################################
#               Load data               #
#########################################

data <- read.table(ologramOutput, sep = "\t", header = TRUE)
data$summed_bp_overlaps_padj <- p.adjust(data$summed_bp_overlaps_pvalue, method = "BH")
data <- data[which(data$summed_bp_overlaps_padj <= 0.05),]

# We check for unique transcription factor
tf <- c()
listeTF <- c()
nbTF <- c()

for(i in 1:nrow(data)){
    myTF <- as.character(data$feature_type[i])
    myTF <- strsplit(myTF, split = "[ + ]")[[1]]
    end <- length(myTF) - 2
    myTF <- myTF[2:end]
    myTF <- unique(myTF)
    myTF <- myTF[! myTF %in% ""]

    # We save unique TF
    tf <- c(tf, myTF)
    tf <- unique(tf)
    nbTF <- c(nbTF, length(myTF))
    # We create a better notation
    myTF <- paste(myTF, collapse=",")
    listeTF <- c(listeTF, myTF)
}

data$listeTF <- listeTF
data$nbTF <- nbTF

data <- data[order(data$nbTF, data$listeTF),]
data <- data[which(data$combination_order <= complexSize + 1),]

write.table(data, paste(path, "Ologram_subset.tsv", sep = ""), sep = "\t", row.names=FALSE, quote = FALSE)