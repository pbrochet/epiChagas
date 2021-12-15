library(gplots)

# This function plot a heatmap according to 
# the data provided

myHeatmap <- function(data, path, phenoTable, name){
    rownames(phenoTable) <- phenoTable$ID
    phenoTable <- phenoTable[colnames(data),]
    tiff(paste(path, "Heatmap_", name, ".tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
    heatmap <- heatmap.2(as.matrix(data), 
            scale="row", 
            hclust=function(x) hclust(x,method="average"), 
            distfun=function(x) as.dist((1-cor(t(x), method = "spearman"))), 
            trace="none", 
            density="none", 
            labRow="",
            cexCol=1.5, 
            col = greenred,
            colCol = as.character(phenoTable$Color),
            srtCol=45)
    dev.off()
}