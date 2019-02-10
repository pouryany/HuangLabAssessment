wgcna.reporter <- function(mRNA.fitlered, nameID = "mRNA"){
  
    powers <-  c(1:10)
    sft    <- pickSoftThreshold(t(mRNA.fitlered),dataIsExpr = TRUE,powerVector = powers,
                               corFnc = cor,corOptions = list(use = 'p'))
    
    
    sizeGrWindow(9, 5)
    par(mfrow = c(1,2))
    cex1 = 0.9
    pdf(paste0("images/",nameID,"powerlaw.pdf"))
    # Plot the results
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit, signed R^2",
         type="n", main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    
    # Red line corresponds to using an R^2 cut-off
    abline(h=0.80,col="red")
    
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
         ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()
    
 
    soft.pwr <- which( sft$fitIndices[,2] ==(max(sft$fitIndices[,2])))
    
    TOM      <- TOMsimilarityFromExpr(t(mRNA.fitlered), power = 4);
    library("flashClust")
    #hierarchical clustering of the genes based on the TOM dissimilarity measure
    colnames(TOM) =rownames(TOM) =rownames(mRNA.fitlered)
    dissTOM=1-TOM
    
    geneTree = flashClust(as.dist(dissTOM),method="average");
    
    #plot the resulting clustering tree (dendrogram)
    pdf(paste0("images/",nameID,"wgcnaDendo.pdf"))
    plot(geneTree, xlab="", sub="",cex=0.3);
    dev.off()
    # Set the minimum module size
    minModuleSize = 20;
    
    # Module identification using dynamic tree cut
    
    dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
    #dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
    
    #the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
    table(dynamicMods)
    dynamicColors = labels2colors(dynamicMods)
    table(dynamicColors)
    pdf(paste0("images/",nameID,"wgcnaDendoCluster.pdf"))
    plotDendroAndColors(geneTree, dynamicColors,
                        "Dynamic Tree Cut", 
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05, 
                        main = "Gene dendrogram and module colors")
    dev.off()
    
    diag(dissTOM) = NA;
    
    module_colors= setdiff(unique(dynamicColors), "grey")
    for (color in module_colors){
        module=rownames(mRNA.fitlered)[which(dynamicColors==color)]
        write.table(module, paste("Data/WGCNA/module_",color, ".txt",sep=""),
                    sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
    }
    
    return(dynamicColors)  
}
