rm(list = ls())
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(CADIA)
library(dplyr)

# Parsing and inspecting the Mutation, Expression, and Protein Data

mutated    <- readRDS("Data/mutationcleaned.rds")
expression <- readRDS("Data/mRNAcleaned.rds")
proteins2  <- readRDS("Data/protcleaned.rds")




# Analysis of highly variable genes


vars        <- apply(expression, 1, var)
sd(vars)
expression1 <- expression[vars > 0,]
expression2 <- expression[vars > 2,]
geneList    <- (rownames(expression1))
deGenes     <- (rownames(expression2))

set.seed(1)
gene.df    <- bitr(geneList, fromType = "SYMBOL",
                   toType = c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)
deGenes.df <- bitr(deGenes, fromType = "SYMBOL",
                   toType = c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)



##  network based analysis on KEGG pathways using CADIA
set.seed(1)
cadia.res  <- CADIA::causalDisturbance(deGenes.df$ENTREZID,gene.df$ENTREZID,
                                       iter = 5000)



# Just keeping a list of enriched pathways and their genes for future.
cadia.res %<>% filter(.,cadia < 0.05) %>% select(., Name, KEGGID)
cadia.res

## CADIA uniquely identifies enrichment of a number of pathways including
# Wnt Signaling,and PI3K-Akt signaling pathway. 


# Moving on to proteins analysis.
# Loading cytokines and immune proteins. As downloaded form AmiGO

    # Immunoproteins
imm.tab   <- read.csv("Data/immunoproteins",header = F, sep = "\t")
imm.prot  <- as.character(imm.tab$V2)
imm.prot2 <- intersect(imm.prot, rownames(proteins2))

imm.prots.all     <- bitr(imm.prot, fromType = "SYMBOL",
                        toType = c("ENTREZID","ENSEMBL"),
                        OrgDb = org.Hs.eg.db)



    # Cytokines
# cyt.tab   <- read.csv("Data/cytokines",header = F, sep = "\t")
# cyt.prot  <- as.character(cyt.tab$V2)
# cyt.prot  <- intersect(cyt.prot, rownames(proteins2))
# cyt.imm   <- union(cyt.prot,imm.prot2)
# cyt.imm2  <- intersect(cyt.imm, rownames(proteins2))
# 


# Aligning the two pieces of data using the immunoproteins. 


mRNA.norm <- expression2[intersect(imm.prot2, rownames(expression2)),]
prot.norm <- proteins2[intersect(imm.prot2, rownames(expression2)),]

mRNA.norm <- mRNA.norm[,intersect(colnames(proteins2),
                                  colnames(expression2))]
prot.norm <- prot.norm[,intersect(colnames(proteins2),
                                  colnames(expression2))]


# Correlation Analysis

prot.mRNA.cors <- sapply(1:nrow(mRNA.norm), function(X){
    cor(t(mRNA.norm[X,]),
        t(prot.norm[X,]),
        use='complete.obs')})



prot.mRNA.pvals <- sapply(1:nrow(mRNA.norm), function(X){
    zz <- cor.test(t(mRNA.norm[X,]),
                   t(prot.norm[X,]))
    return(zz$p.value)})

plotDensities(prot.mRNA.cors)


nrow(mRNA.norm)
nrow(expression)
nrow(proteins2)
ncol(prot.norm)
ncol(mRNA.norm)


# We have 163 mRNA and proteins signals are variable and consistent
# between the two datasets. This is acroos 74 samples. 


prot.mRNA.fdr <- p.adjust(prot.mRNA.pvals)
plotDensities(prot.mRNA.cors)
plotDensities(prot.mRNA.cors[prot.mRNA.fdr < 0.005])
length(prot.mRNA.cors[prot.mRNA.fdr < 0.005])

## For 147 of the samples we have a consistent positive correlation
## Let's just work on these. 

mRNA.filtered <- mRNA.norm[prot.mRNA.fdr < 0.005 & prot.mRNA.cors > 0.4,]
prot.filtered <- prot.norm[prot.mRNA.fdr < 0.005 & prot.mRNA.cors > 0.4,]




### Let's dive into some clustering. The data is high dimentional


library(mclust)

BIC <- mclustBIC(t(mRNA.filtered))
plot(BIC)
summary(BIC)

BIC <- mclustBIC(t(prot.filtered))
plot(BIC)
summary(BIC)



mod.mRNA <- Mclust(t(mRNA.filtered), G = 4, modelName = "VEI")
mod.prot <- Mclust(t(prot.filtered), G = 4, modelName = "VEI")






saveRDS(object = prot.filtered , "Data/protfiltered.rds")
saveRDS(object = mRNA.filtered , "Data/mRNAfiltered.rds")
saveRDS(object = mod.prot      ,"Data/modprot.rds")
saveRDS(object = mod.mRNA      ,"Data/modmRNA.rds")
saveRDS(object = imm.prots.all ,"Data/immprotsall.rds")


#
#
#
#
#



# What if I could define a joint data matrix?

zzz  <- (data.matrix(mRNA.filtered)) %*% t(data.matrix(prot.filtered))
zzz  <- zzz%*% t(zzz)

# Clustering of this matrix may prove useful. I leave that for some later time
heatmap3(zzz, scale = "none", col = coul, showColDendro = F,showRowDendro = F)



## So let's do some dimentionality reduction. 



##  PCA on mRNA data
    mRNA.pca <- prcomp(t(mRNA.filtered), center = TRUE, scale. = TRUE)
    plot(mRNA.pca,type = "l")

# The plot tells us after the 4th PC we can't get much  variance explained

    mRNA.sum <- mRNA.pca$x[,1:4]
# Let's work on 4 PCs


        
    library(ggbiplot)
    set.seed(1)
    ggbiplot(mRNA.pca , obs.scale = 1, var.scale =1 , var.axes = F,
                  ellipse = T,groups = as.factor(mod.mRNA$classification))+
        geom_point(aes( color =  as.factor(mod.mRNA$classification)),size = 3) +
        theme_bw()
    ggsave("images/PCAmRNA.jpg")
    library("pca3d")
    
    
    pca3d(mRNA.pca,group = mod.mRNA$classification)
    
 
    
##  PCA on prot data
prot.pca <- prcomp(t(prot.filtered), center = TRUE, scale. = TRUE)
plot(prot.pca,type = "l")

# Similar to mRNA, the elbow happens arond 4th PC
# after the 4th PC we can't get much  variance explained

prot.sum <- prot.pca$x[,1:4]
# Let's work on 4 PCs

set.seed(1)
ggbiplot(prot.pca , obs.scale = 1, var.scale =1 , var.axes = F,
              ellipse = T,groups = as.factor(mod.prot$classification))+
    geom_point(aes( color =  as.factor(mod.prot$classification))
               ,size = 3) +
    theme_bw()

pca3d(prot.pca,group = mod.mRNA$classification)










