# This script do methylation analysis from Illumina EPIC kit. It allow analysis of DMP (differentially
# methylated position) and DMR (differentially methylated region) with CHAMP R package
# (which include DMRCate tool).


###############################################
#               Load parameters               #
###############################################

set.seed(123)

args <- commandArgs(trailingOnly=TRUE)
dirResults <- args[1]
path <- paste(dirResults, "/", sep = "")
dataName <- args[2]


data <- read.table(dataName, sep = "\t", header = TRUE)
data$summed_bp_overlaps_padj <- p.adjust(data$summed_bp_overlaps_pvalue, method = "BH")
sign <- data[which(data$summed_bp_overlaps_padj <= 0.05 & data$summed_bp_overlaps_padj > 0),]

write.table(sign, paste(dirResults, "/Ologram_significative_TF.tab", sep = ""), sep = "\t", row.names=FALSE)