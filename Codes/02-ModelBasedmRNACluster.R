rm(list = ls())
library(CADIA)
library(heatmap3)
library(mclust)
library(RColorBrewer)


prot.filtered <- readRDS("Data/protfiltered.rds")
mRNA.filtered <- readRDS("Data/mRNAfiltered.rds")
mod.prot      <- readRDS("Data/modprot.rds")
mod.mRNA      <- readRDS("Data/modmRNA.rds")
imm.prots.all <- readRDS("Data/immprotsall.rds")


clust.order <- unlist(tapply(1:ncol((mRNA.filtered)),
                             as.factor(mod.mRNA$classification),I)
                      ,use.names = F)
mRNA.mat    <- data.matrix(mRNA.filtered)[,clust.order]



# Just some color setting for visualization

my_group    <- as.numeric(as.factor(mod.mRNA$classification))
my_col1     <- brewer.pal(8, "Set1")[my_group]
my_col1     <- my_col1[clust.order]

my_group    <- as.numeric(as.factor(mod.prot$classification))
my_col2     <- brewer.pal(8, "Set2")[my_group]
my_col2     <- my_col2[clust.order]

coul        <- colorRampPalette(brewer.pal(10, "RdBu"))(50)
my_col      <- cbind("mRNA"= my_col1, "Protein" = my_col2)



# Distance based clustering 

d  <- dist( data.matrix(mRNA.filtered ), method = "euclidian")
hc <- hclust(d,method = "complete")

plot(hc,labels=rownames(mRNA.filtered),cex=0.5)

# Three seems suitable so lets just do the rest of coloring and heatmap   
hclusters  <- cutree(hc, h=80)
row.order  <- unlist(tapply(1:nrow((mRNA.filtered)),
                            as.factor(hclusters),I),use.names = F)

mRNA.mat   <- mRNA.mat[row.order,]

my_group   <- as.numeric(as.factor(hclusters))
my_col1    <- brewer.pal(9, "Set1")[my_group]
my_row_col <- brewer.pal(9, "Set1")[my_group]
my_row_col <- my_row_col[row.order]

jpeg("images/Model-basedmRNACluster.jpg")
heatmap3(mRNA.mat,Colv = NA,Rowv = NA,showColDendro = F,showRowDendro = F,
         scale = "none", col = coul,RowSideColors =  my_row_col,
         ColSideColors=my_col, RowSideLabs = "mRNA Probs")
dev.off()

clust1 <- names(hclusters[hclusters ==1])
clust2 <- names(hclusters[hclusters ==2])
clust3 <- names(hclusters[hclusters ==3])

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
cadia.res1  <- CADIA::causalDisturbance(gene.clust1$ENTREZID,
                                        imm.prots.all$ENTREZID,
                                        iter = 5000)
set.seed(1)
cadia.res2  <- CADIA::causalDisturbance(gene.clust2$ENTREZID,
                                        imm.prots.all$ENTREZID,
                                        iter = 5000)
set.seed(1)
cadia.res3  <- CADIA::causalDisturbance(gene.clust3$ENTREZID,
                                        imm.prots.all$ENTREZID,
                                        iter = 5000)

# The three genes cluster repsent different functionalities 

# Cluster 1 is associated with T cell receptor pathway and
# Th17 cell differentiation.

# Cluster 2 is associated with several signaling pathways including Wnt,
# Calcium, Phospolipase D, cAMP,...


# Cluster 3 is associated with Focal adhesion
# Endocrine resistance, ECM-receptor interaction

