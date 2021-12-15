library(ggplot2)
library(gProfileR)

# This function perform an enrichment analysis on the 
# DEGs, plot the top results for the three GO category
# and save the results in a file

myGO <- function(data, path, name){
    # Perform enrichment with gProfileR
    gprofil <- gprofiler(data,organism = "hsapiens")
    gprofil <- gprofil[order(gprofil$p.value),]

    # Select the top five function in the three category
    BP <- gprofil[gprofil$domain == "BP",]
    BP <- BP[1:5, c("p.value", "domain", "term.name")]
    CC <- gprofil[gprofil$domain == "CC",]
    CC <- CC[1:5, c("p.value", "domain", "term.name")]
    MF <- gprofil[gprofil$domain == "MF",]
    MF <- MF[1:5, c("p.value", "domain", "term.name")]

    # Merge the result in a new table
    enrich <- rbind(BP, CC)
    enrich <- rbind(enrich, MF)
    enrich$term.name = with(enrich, factor(term.name, levels=term.name[order(ave(-p.value, domain, FUN=min),-p.value)]))
    enrich$log <- -log(enrich$p.value)

    # Plo the results
    p <- ggplot(enrich, aes(x = term.name, y = log, fill = domain)) +
        geom_bar(stat="identity", position = "dodge") +
        coord_flip() +
        xlab("Gene ontology") +
        ylab("-log10(p.value)")
    tiff(paste(path, "Top_GO_", name, ".tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
    plot(p)
    dev.off()

    # Save the results
    write.table(gprofil, paste(path, "GO_enrichment_", name, ".tab", sep = ""), sep = "\t", row.names=FALSE)

}



