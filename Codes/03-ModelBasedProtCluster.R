rm(list = ls())
library(CADIA)
library(heatmap3)
library(mclust)
library(RColorBrewer)


prot.filtered <- readRDS("Data/protfiltered.rds")
mod.prot      <- readRDS("Data/modprot.rds")
mod.mRNA      <- readRDS("Data/modmRNA.rds")
imm.prots.all <- readRDS("Data/immprotsall.rds")


clust.order <- unlist(tapply(1:ncol((prot.filtered)),
                             as.factor(mod.prot$classification),I)
                      ,use.names = F)
prot.mat    <- data.matrix(prot.filtered)[,clust.order]



# Just some color setting for visualization

my_group    <- as.numeric(as.factor(mod.prot$classification))
my_col1     <- brewer.pal(8, "Set2")[my_group]
my_col1     <- my_col1[clust.order]

my_group    <- as.numeric(as.factor(mod.mRNA$classification))
my_col2     <- brewer.pal(8, "Set1")[my_group]
my_col2     <- my_col2[clust.order]

coul        <- colorRampPalette(brewer.pal(10, "RdBu"))(50)
my_col      <- cbind("prot"= my_col1, "mRNA" = my_col2)



# Distance based clustering 

d  <- dist( data.matrix(prot.filtered ), method = "euclidian")
hc <- hclust(d,method = "complete")

plot(hc,labels=rownames(prot.filtered),cex=0.5)

# Two seems suitable so lets just do the rest of coloring and heatmap   
hclusters  <- cutree(hc, h=45)
row.order  <- unlist(tapply(1:nrow((prot.filtered)),
                            as.factor(hclusters),I),use.names = F)

prot.mat   <- prot.mat[row.order,]

my_group   <- as.numeric(as.factor(hclusters))
my_col1    <- brewer.pal(9, "Set1")[my_group]
my_row_col <- brewer.pal(9, "Set1")[my_group]
my_row_col <- my_row_col[row.order]

pdf("images/Model-basedProteinCluster.pdf")
heatmap3(prot.mat,Colv = NA,Rowv = NA,showColDendro = F,showRowDendro = F,
         scale = "none", col = coul,RowSideColors =  my_row_col,
         ColSideColors=my_col)
dev.off()


clust1 <- names(hclusters[hclusters ==1])
clust2 <- names(hclusters[hclusters ==2])


gene.clust1     <- bitr(clust1, fromType = "SYMBOL",
                        toType = c("ENTREZID","ENSEMBL"),
                        OrgDb = org.Hs.eg.db)
gene.clust2     <- bitr(clust2, fromType = "SYMBOL",
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




