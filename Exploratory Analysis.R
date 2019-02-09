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



enriched   <- enrichGO(  gene          = deGenes.df$ENSEMBL,
                         universe      = gene.df$ENSEMBL,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.5,
                         readable      = TRUE)
#

#
 head(enriched,10)
 ego <- gofilter(enriched,4)
 head(ego)
 ego <- summary(ego)
 ego[ego$p.adjust < 0.01,]$Description
 # Too many terms the first three most significant enrichments might be 
 # Informative. Let's take a more focused approach.
 rm(list = c("ego","enriched"))

##  network based analysis on KEGG pathways using CADIA
set.seed(1)
cadia.res  <- CADIA::causalDisturbance(deGenes.df$ENTREZID,gene.df$ENTREZID,
                                       iter = 5000)



# Just keeping a list of enriched pathways and their genes for future.
cadia.res %<>% filter(.,cadia < 0.05) %>% select(., Name, KEGGID)

pathLists <- CADIA::geneReport(deGenes.df$ENTREZID,cadia.res$KEGGID)
pathLists <- as_tibble(pathLists)
pathLists <- mutate(pathLists,genes = as.character(genes))
thislist  <- data.frame()
for (i in 1:nrow(pathLists)){
    thislist1 <- unlist(strsplit(unlist(pathLists[i,2]), split = "[/]"))
    thislist1 <- mapIds(org.Hs.eg.db, keys=thislist1, column="SYMBOL",
                        keytype="ENTREZID", multiVals="first")
    thislist1 <- data.frame(as.character(pull(pathLists[i,1])),thislist1)
    thislist  <- rbind(thislist,thislist1)
}



## CADIA uniquely identifies enrichment of a number of pathways including
# Wnt Signaling,and PI3K-Akt signaling pathway. 


# Moving on to proteins analysis.
# Loading cytokines and immune proteins. As downloaded form AmiGO

    # Immunoproteins
imm.tab   <- read.csv("Data/immunoproteins",header = F, sep = "\t")
imm.prot  <- as.character(imm.tab$V2)
imm.prot2 <- intersect(imm.prot, rownames(proteins2))

    # Cytokines
cyt.tab   <- read.csv("Data/cytokines",header = F, sep = "\t")
cyt.prot  <- as.character(cyt.tab$V2)
cyt.imm   <- union(cyt.prot,imm.prot)
cyt.imm2  <- intersect(cyt.imm, rownames(proteins2))



# Aligning the two pieces of data. 


mRNA.norm <- expression2[intersect(cyt.imm2, rownames(expression2)),]
prot.norm <- proteins2[intersect(cyt.imm2, rownames(expression2)),]

mRNA.norm <- mRNA.norm[,intersect(colnames(proteins2),
                                  colnames(expression2))]
prot.norm <- prot.norm[,intersect(colnames(proteins2),
                                  colnames(expression2))]



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
nrow(proteins)
ncol(prot.norm)
ncol(mRNA.norm)


# We have 190 mRNA and proteins signals are variable and consistent
# between the two datasets. This is acroos 74 samples. 


prot.mRNA.fdr <- p.adjust(prot.mRNA.pvals)
plotDensities(prot.mRNA.cors[prot.mRNA.fdr < 0.05])
length(prot.mRNA.cors[prot.mRNA.fdr < 0.05])

## For 147 of the samples we have a consistent positive correlation
## Let's just work on these. 

mRNA.fitlered <- mRNA.norm[prot.mRNA.fdr < 0.05 & prot.mRNA.cors > 0.4,]
prot.fitlered <- prot.norm[prot.mRNA.fdr < 0.05 & prot.mRNA.cors > 0.4,]




### Let's dive into some clustering. The data is high dimentional


library(mclust)

BIC <- mclustBIC(t(mRNA.fitlered))
plot(BIC)
summary(BIC)



BIC <- mclustBIC(t(prot.fitlered))
plot(BIC)
summary(BIC)




mod.mRNA <- Mclust(t(mRNA.fitlered), G = 3, modelName = "VEI")
mod.prot <- Mclust(t(prot.fitlered), G = 3, modelName = "VEI")

mod.mRNA$classification


BIC <- mclustBIC(mRNA.sum)
plot(BIC)
summary(BIC)

BIC <- mclustBIC(prot.sum)
plot(BIC)
summary(BIC)

boot <- MclustBootstrap(mod, nboot = 999, type = "bs")
summary(boot, what = "se")


library(RColorBrewer)
library("gplots")

my_group=as.numeric(as.factor(mod.prot$classification))
my_col=brewer.pal(9, "Set1")[my_group]
coul = colorRampPalette(brewer.pal(10, "RdBu"))(25)

heatmap.2(data.matrix(mRNA.fitlered), scale = "none", col = coul,
          ColSideColors=my_col)











## So let's do some dimentionality reduction. 



##  PCA on mRNA data
mRNA.pca <- prcomp(t(mRNA.fitlered), center = TRUE, scale. = TRUE)
plot(mRNA.pca,type = "l")

# The plot tells us after the 4th PC we can't get much  variance explained

mRNA.sum <- mRNA.pca$x[,1:4]
# Let's work on 4 PCs

##  PCA on prot data
prot.pca <- prcomp(t(prot.fitlered), center = TRUE, scale. = TRUE)
plot(prot.pca,type = "l")

# Similar to mRNA, the elbow happens arond 4th PC
# after the 4th PC we can't get much  variance explained

prot.sum <- prot.pca$x[,1:4]
# Let's work on 4 PCs







library(RColorBrewer)



hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

# perform clustering on rows and columns
cl.row <- hclustfunc(distfunc(data.matrix(mRNA.fitlered)))
cl.col <- hclustfunc(distfunc(t(data.matrix(mRNA.fitlered))))


library("gplots")

heatmap.2(data.matrix(mRNA.fitlered), scale = "none", col = coul)

## Some PCA plot
seq.pca <- prcomp(t(mRNA.fitlered), center = TRUE, scale. = TRUE)
plot(seq.pca,type = "l")

# The plot tells us after the 5th PC we can't get much additional variance explained

scoreTot <- seq.pca$x[,1:4]
# Let's work on 5 PCs


library(mclust)

BIC <- mclustBIC(scoreTot)
plot(BIC)
summary(BIC)

ICL <- mclustICL(scoreTot)
summary(ICL)


mod <- Mclust(scoreTot, G = 2, modelName = "EEV")
summary(mod, parameters = TRUE)

plot(mod, what = "classification", main = FALSE)
boot <- MclustBootstrap(mod, nboot = 999, type = "bs")
summary(boot, what = "se")


modTot <- kmeans(scoreTot, 2)



library(ggbiplot)
set.seed(1)
g <- ggbiplot(seq.pca , obs.scale = 1, var.scale =1 , var.axes = F,
              ellipse = T,groups = as.factor(modTot$cluster))+
    geom_point(aes( color = as.factor(modTot$cluster)), size = 3) + theme_bw()





