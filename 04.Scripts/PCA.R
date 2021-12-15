library(ggplot2)
library(factoextra)
library(FactoMineR)
library(dplyr)
library(ggrepel)

# This function perform a PCA analysis and plot the results according to 
# the data provided

myPCA <- function(data, path, phenoTable, sex, age, name){
    data <- data[which(apply(data, 1, var, na.rm=TRUE) != 0),]
    # data <- data[,which(apply(data, 2, var, na.rm=TRUE) != 0)]
    data <- data[apply(data, 1, function(x) !all(x==0)),]
    # Run the PCA analysis
    dataPCA <- as.data.frame(t(data))
    PCA <- prcomp(dataPCA, scale = TRUE)

    # Extract the eigen value and create the plot for the PCA, according to
    # the 1st and 2nd dimension
    contrib <- get_eigenvalue(PCA)
    PCAplot <- data.frame(PC1 = PCA$x[,"PC1"], PC2 = PCA$x[,"PC2"])
    PCAplot$ID <- rownames(PCAplot)
    PCAplot <- inner_join(phenoTable, PCAplot)

    myColorInfo <- unique(PCAplot[,c("Phenotype", "Color")])
    myCols <- as.character(myColorInfo$Color)
    names(myCols) <- as.character(myColorInfo$Phenotype)

    p <- ggplot(data = PCAplot, aes(x = PC1, y = PC2, color = Phenotype)) +
        geom_point(alpha = 0.8) +
        scale_color_manual("Phenotype", values = myCols) +
        xlab(paste("Principal component 1 (", round(contrib$variance.percent[1], 2), "%)", sep = "")) +
        ylab(paste("Principal component 2 (", round(contrib$variance.percent[2], 2), "%)", sep = "")) +
        theme_classic() +
        geom_text(aes(label=ID), size = 4,position=position_jitter(width=1,height=1)) +
        geom_hline(yintercept = 0, lty = 2) +
        geom_vline(xintercept = 0, lty = 2) +
        xlim(min(PCAplot$PC1) - 30*abs(min(PCAplot$PC1))/100, max(PCAplot$PC1) + 30*abs(max(PCAplot$PC1))/100) +
        theme(axis.title.x = element_text(size = rel(1.8)),
                axis.title.y = element_text(size = rel(1.8)),
                axis.text.x = element_text(size = rel(1.8)),
                axis.text.y = element_text(size = rel(1.8)))
        # geom_text_repel(aes(label=ID), size = 5, max.overlaps = Inf)


    tiff(paste(path, "PCA_", name, ".tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
    plot(p)
    dev.off()

    # PCA colors according to sex
    if(sex == TRUE){
        PCAplot <- data.frame(PC1 = PCA$x[,"PC1"], PC2 = PCA$x[,"PC2"])
        PCAplot$ID <- rownames(PCAplot)
        PCAplot <- inner_join(phenoTable, PCAplot)
        PCAplot <- PCAplot[!is.na(PCAplot$Sex),]
        PCAplot <- PCAplot[!PCAplot$Sex %in% "",]
        p <- ggplot(data = PCAplot, aes(x = PC1, y = PC2, color = Sex)) +
            geom_point(alpha = 0.8) +
            scale_color_manual("Sex", values=c("deeppink2", "blue")) +
            xlab(paste("Principal component 1 (", round(contrib$variance.percent[1], 2), "%)", sep = "")) +
            ylab(paste("Principal component 2 (", round(contrib$variance.percent[2], 2), "%)", sep = "")) +
            theme_classic() +
            geom_text(aes(label=ID), size = 4,position=position_jitter(width=1,height=1)) +
            geom_hline(yintercept = 0, lty = 2) +
            geom_vline(xintercept = 0, lty = 2) +
            xlim(min(PCAplot$PC1) - 30*abs(min(PCAplot$PC1))/100, max(PCAplot$PC1) + 30*abs(max(PCAplot$PC1))/100) +
            theme(axis.title.x = element_text(size = rel(1.8)),
                    axis.title.y = element_text(size = rel(1.8)),
                    axis.text.x = element_text(size = rel(1.8)),
                    axis.text.y = element_text(size = rel(1.8)))
            # geom_text_repel(aes(label=ID), size = 5, max.overlaps = Inf)
        tiff(paste(path, "PCA_sex_", name, ".tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
        plot(p)
        dev.off()
    }

    # PCA colors according to age
    if(age == TRUE){
        PCAplot <- data.frame(PC1 = PCA$x[,"PC1"], PC2 = PCA$x[,"PC2"])
        PCAplot$ID <- rownames(PCAplot)
        PCAplot <- inner_join(phenoTable, PCAplot)
        PCAplot <- PCAplot[!is.na(PCAplot$Age),]
        p <- ggplot(data = PCAplot, aes(x = PC1, y = PC2, color = Age)) +
            geom_point(alpha = 0.8) +
            scale_color_gradient2(low="grey", mid="skyblue",high="darkblue", midpoint=40) +
            xlab(paste("Principal component 1 (", round(contrib$variance.percent[1], 2), "%)", sep = "")) +
            ylab(paste("Principal component 2 (", round(contrib$variance.percent[2], 2), "%)", sep = "")) +
            theme_classic() +
            geom_text(aes(label=ID), size = 4,position=position_jitter(width=1,height=1)) +
            geom_hline(yintercept = 0, lty = 2) +
            geom_vline(xintercept = 0, lty = 2) +
            xlim(min(PCAplot$PC1) - 30*abs(min(PCAplot$PC1))/100, max(PCAplot$PC1) + 30*abs(max(PCAplot$PC1))/100) +
            theme(axis.title.x = element_text(size = rel(1.8)),
                    axis.title.y = element_text(size = rel(1.8)),
                    axis.text.x = element_text(size = rel(1.8)),
                    axis.text.y = element_text(size = rel(1.8)))
            # geom_text_repel(aes(label=ID), size = 5, max.overlaps = Inf)
        tiff(paste(path, "PCA_age_", name, ".tiff", sep = ""), width = 10, height = 10, units = 'in', res = 300)
        plot(p)
        dev.off()
    }
}













