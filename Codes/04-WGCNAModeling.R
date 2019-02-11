# Just doing some WGCNA analysis
rm(list= ls())
library(WGCNA)
library(CADIA)
library(heatmap3)
library(mclust)
library(RColorBrewer)


prot.filtered <- readRDS("Data/protfiltered.rds")
mRNA.filtered <- readRDS("Data/mRNAfiltered.rds")

mod.prot      <- readRDS("Data/modprot.rds")
mod.mRNA      <- readRDS("Data/modmRNA.rds")
imm.prots.all <- readRDS("Data/immprotsall.rds")


source("Codes/WGCNArunner.R")

### For protein data 
dynamicColors <- wgcna.reporter(prot.filtered,"prot")
coul          <- colorRampPalette(brewer.pal(10, "RdBu"))(50)



module.order <- unlist(tapply(1:ncol(t(prot.filtered)),as.factor(dynamicColors),I))
#m<-t(t(t(prot.filtered)[,module.order])/apply(t(prot.filtered)[,module.order],2,max))
m <- t(prot.filtered)[,module.order]
dynamic.col.color <- brewer.pal(4, "Set1")[mod.prot$classification]
col.order <- unlist(tapply(1:ncol((prot.filtered)),as.factor(dynamic.col.color),I))

m <- m[col.order,]


pdf("images/ProteinClustersWGCNA.pdf")
heatmap3(t(m),col=coul,Rowv=NA,Colv=NA,scale=c("none"),
         RowSideColors=dynamicColors[module.order],
         ColSideColors=dynamic.col.color[col.order],
         ColSideLabs = "Sample Clusters",
         RowSideLabs = "Proteins Clusters")
dev.off()





prot.res <- rownames(prot.filtered)
dynaCos  <- as.character(dynamicColors)
clust1   <- (prot.res[dynaCos == "grey" ])
clust2   <- (prot.res[dynaCos == "turquoise"])
clust3   <- (prot.res[dynaCos == "blue"])

gene.clust1     <- bitr(clust1, fromType = "SYMBOL",
                        toType = c("ENTREZID","ENSEMBL"),
                        OrgDb = org.Hs.eg.db)
gene.clust2     <- bitr(clust2, fromType = "SYMBOL",
                        toType = c("ENTREZID","ENSEMBL"),
                        OrgDb = org.Hs.eg.db)
gene.clust3     <- bitr(clust3, fromType = "SYMBOL",
                        toType = c("ENTREZID","ENSEMBL"),
                        OrgDb = org.Hs.eg.db)


set.seed(1)
cadia.res1  <- CADIA::causalDisturbance(gene.clust1$ENTREZID,imm.prots.all$ENTREZID,
                                        iter = 5000)
set.seed(1)
cadia.res2  <- CADIA::causalDisturbance(gene.clust2$ENTREZID,imm.prots.all$ENTREZID,
                                        iter = 5000)
set.seed(1)
cadia.res3  <- CADIA::causalDisturbance(gene.clust3$ENTREZID,imm.prots.all$ENTREZID,
                                        iter = 5000)


# The three genes cluster repsent different functionalities 

# Cluster 1 is associated with ECM Receptor interactions.

# Cluster 2's associations are weak


# Cluster 3 is associated with NOD-like receptor signaling,T cell receptor












### For mRNA data 
dynamicColors <- wgcna.reporter(mRNA.filtered,"mRNA")



module.order <- unlist(tapply(1:ncol(t(mRNA.filtered)),as.factor(dynamicColors),I))
#m<-t(t(t(mRNA.filtered)[,module.order])/apply(t(mRNA.filtered)[,module.order],2,max))
m <- t(mRNA.filtered)[,module.order]
dynamic.col.color <- brewer.pal(4, "Set1")[mod.mRNA$classification]
col.order <- unlist(tapply(1:ncol((mRNA.filtered)),as.factor(dynamic.col.color),I))

m <- m[col.order,]


pdf("images/mRNAClustersWGCNA.pdf")
heatmap3(t(m),col=coul,Rowv=NA,Colv=NA,scale=c("none"),
         RowSideColors=dynamicColors[module.order],
         ColSideColors=dynamic.col.color[col.order],
         ColSideLabs = "Sample Clusters",
         RowSideLabs = "mRNA Clusters")
dev.off()


mRNA.res <- rownames(mRNA.filtered)
dynaCos  <- as.character(dynamicColors)
clust1   <- (mRNA.res[dynaCos == "grey" ])
clust2   <- (mRNA.res[dynaCos == "turquoise"])
clust3   <- (mRNA.res[dynaCos == "blue"])
clust4   <- (mRNA.res[dynaCos == "brown"])
gene.clust1     <- bitr(clust1, fromType = "SYMBOL",
                        toType = c("ENTREZID","ENSEMBL"),
                        OrgDb = org.Hs.eg.db)
gene.clust2     <- bitr(clust2, fromType = "SYMBOL",
                        toType = c("ENTREZID","ENSEMBL"),
                        OrgDb = org.Hs.eg.db)
gene.clust3     <- bitr(clust3, fromType = "SYMBOL",
                        toType = c("ENTREZID","ENSEMBL"),
                        OrgDb = org.Hs.eg.db)
gene.clust4     <- bitr(clust4, fromType = "SYMBOL",
                        toType = c("ENTREZID","ENSEMBL"),
                        OrgDb = org.Hs.eg.db)



set.seed(1)
cadia.res1  <- CADIA::causalDisturbance(gene.clust1$ENTREZID,imm.prots.all$ENTREZID,
                                        iter = 5000)
set.seed(1)
cadia.res2  <- CADIA::causalDisturbance(gene.clust2$ENTREZID,imm.prots.all$ENTREZID,
                                        iter = 5000)
set.seed(1)
cadia.res3  <- CADIA::causalDisturbance(gene.clust3$ENTREZID,imm.prots.all$ENTREZID,
                                        iter = 5000)
set.seed(1)
cadia.res4  <- CADIA::causalDisturbance(gene.clust4$ENTREZID,imm.prots.all$ENTREZID,
                                        iter = 5000)

# The three genes cluster repsent different functionalities 

# Cluster 1 is associated with ECM Receptor interactions.

# Cluster 2's associations are weak


# Cluster 3 is associated with NOD-like receptor signaling,T cell receptor

# Cluster 4 ECM-Receptor interaction and focal adhesion.


