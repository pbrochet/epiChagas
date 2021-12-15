library(pathview)
library(gage)
library(gageData)
library(biomaRt)
library(dplyr)
library(httr)

myKegg <- function(data, path, name){
    # 1. Convert gene ID in entrez ID
    degID <- as.character(data$EnsemblID)
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    degEntrez <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_id"), values=degID, mart= mart, useCache = FALSE)
    colnames(degEntrez) <- c("EnsemblID", "EntrezId")
    data <- full_join(data, degEntrez)

    # 2. Create a vector containing the FC and names by the entrez id
    data <- na.omit(data)
    FC <- data$FC
    names(FC) <- data$EntrezId # 1063

    # 3. Perform the Kegg enrichment
    kegg <- kegg.gsets("hsa")
    kg <- kegg$kg.sets
    kres <- gage(FC, gsets=kg, same.dir=TRUE)
    tableKres <- data.frame(id = rownames(kres$greater), pval = kres$greater[,3], size = kres$greater[,5])
    tableK <- tableKres[order(tableKres$size, decreasing = TRUE),c("pval", "size")]

    # 4. Write the results of Kegg enrichment
    write.table(tableK, paste(path, "Kegg_enrichment_", name, ".tab", sep = ""), sep = "\t")
}