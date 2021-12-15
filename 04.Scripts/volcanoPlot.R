library(ggplot2)

# This function create and save a volcano plot according to 
# the data provided

myVolcanoPlot <- function(data, path, xlab, ylab, cutoffFC, name){
    # We sort data by adjusted p value
    data <- data[order(data$padj),]

    # We create a new column containing a color depending on the parameters
    # (FC, adjusted pvalue)
    col <- c()
    for (i in 1:nrow(data)){
        if(data$FC[i] >= cutoffFC & is.na(data$padj[i]) == FALSE & data$padj[i] <= 0.05){
            col <- c(col, "red")
        }
        else if(data$FC[i] <= -cutoffFC & is.na(data$padj[i]) == FALSE & data$padj[i] <= 0.05){
            col <- c(col, "green")
        }
        else{
            col <- c(col, "black")
        }
    }
    data$Genes <- col

    # We log convert the adjusted pvalue
    data$logpval <- -log10(data$padj) 

    # We create the plot
    p <- ggplot(data=data, aes(FC, logpval, color=Genes))+
        geom_point() +     
        scale_color_manual(labels = c("No differentially expressed", "Down-regulated", "Up-regulated"), values=c("black", "green", "red"))+
        geom_hline(yintercept=1.30, linetype="dashed", color = "black") + 
        geom_vline(xintercept = -cutoffFC, linetype="dotted", color = "black") +
        geom_vline(xintercept = cutoffFC, linetype="dotted", color = "black") +
        labs(x = xlab, y = ylab) +
        theme_classic() +
        theme(axis.title.x = element_text(size = rel(1.8)),
                axis.title.y = element_text(size = rel(1.8)),
                axis.text.x = element_text(size = rel(1.8)),
                axis.text.y = element_text(size = rel(1.8)))

    # We save the plot
    tiff(paste(path, "VolcanoPlot_", name, ".tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
    plot(p)
    dev.off()
}