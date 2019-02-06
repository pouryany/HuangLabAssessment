library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(CADIA)

# Parsing and inspecting the Mutation, Expression, and Protein Data

mutated    <- read.csv("BRCA_MC3_SOMATIC_formatted_amino_acid.txt",
                       stringsAsFactors = F, sep = "\t", row.names = 1)
expression <- read.csv("BRCA_mRNA_formatted_normalized.txt",
                       stringsAsFactors = F, sep = "\t", row.names = 1)
proteins   <- read.csv("BRCA_PRO_formatted_normalized.txt",
                       stringsAsFactors = F, sep = "\t", row.names = 1)



head(mutated)
head(proteins)
head(expression)



# Insepcting mutation data
# Ideally we would like to treat different mutation locations separately
# Let's make a simplifying assumption that the mutation causes loss of function
# Just some peaking into the table and making a Zero/One matrix
    dim(mutated)
    sum(is.na(mutated))
    
    mutated[mutated != "wt"] <- 1
    mutated[mutated == "wt"] <- 0
    
    mutated <- as.data.frame(mutated)
    mutated <- data.matrix(mutated)
    colSums(mutated)
    rowSums(mutated)
    
# TP53 the usual suspect ...
# Let's do some clustering on the mutation matrix and see what comes out
# Why doing clustering? maybe finding subgroups that we can further contrast
# I am not sure where mutations happen together or correlate, but maybe clusters
# together can bring some higher level descriptions.

    heatmap(mutated, Colv=F,scale ="none")
    d <- dist((mutated) )
    clusters <- hclust(d)
    plot(clusters,cex=0.5)
    # We can't get much out this heatmap the way it is

# For now let's just investigate the data.
# I will run a scenario by doing differential expression between
# mutated profiles later



# Parsing and inspecting the mRNA expression Data
# Just inspecting the data distribution. We know that the data is normalized
# So there shouldn't be any problems here. Distributions should align
   

    plotDensities(expression, legend = F)
    boxplot(expression)
    plotMD(expression,column = 1)
    abline(h=0,col="red")


    means  <- apply(expression, 1, mean)
    vars   <- apply(expression, 1, var)
    plot(means,vars)
    abline(h=10,col="red")
    expression2 <- expression[vars !=0,]
    dim(expression)
    dim(expression2)

    # We see some reduction becasue of zero variance genes
    
    plotDensities(vars)
    plot(ecdf(vars))

    # Somewhere in (2,3) the second derivative of the density plot/CDF of 
    #variance becomes positive.let's just pick that as the variance threshold.
    # Let's do some preliminary analysis of mRNA Expressions

    
    expression  <- expression[vars > 0,]
    
    expression2 <- expression[vars > 3,]

    
    
    geneList   <- (rownames(expression))
    deGenes    <- (rownames(expression2))
    gene.df    <- bitr(geneList, fromType = "SYMBOL",
                       toType = c("ENTREZID","ENSEMBL"),
                       OrgDb = org.Hs.eg.db)
    deGenes.df <- bitr(deGenes, fromType = "SYMBOL",
                       toType = c("ENTREZID","ENSEMBL"),
                       OrgDb = org.Hs.eg.db)
    
    length(unique(geneList))
    length((geneList))
    enriched   <- enrichGO( gene          = deGenes.df$ENSEMBL,
                            universe      = gene.df$ENSEMBL,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENSEMBL',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            readable      = TRUE)
    
    
    ego <- enrichKEGG(gene          = deGenes.df$ENTREZID,
                      organism      = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.1)
    
    head(ego)
    clusterProfiler::enrichMap(enriched)
    ego <- gofilter(enriched,3)
    bp2 <- simplify(ego, cutoff=0.5, by="p.adjust", select_fun=min)

    head(bp2)

    rm(list = c("ego","bp2","enriched"))


    set.seed(1)
    cadia.res <- CADIA::causalDisturbance(deGenes.df$ENTREZID,gene.df$ENTREZID,
                                          iter = 10000)
    
    head(proteins)




# Just a prefiltering step on proteins
    head(proteins)
    missings   <- apply(proteins, 1, function(X){
        sum(is.na(X))
    })

    sum(missings > 0)
# There are 3628 proteins with missing values
# We have to filter and potentially impute...

    # For now lets say no more than 30 proteins missing.
    proteins2  <- proteins[!(missings > 0),]



    prot.means <- apply(proteins2, 1, mean, na.rm =T)
    prot.vars   <- apply(proteins2, 1, var,  na.rm =T)

    sum(is.na(prot.means))
    sum(is.na(prot.vars))
    proteins <- proteins[!is.na(prot.means),]


# Just some data inspection. The normalization is described is the handout so
# we should be careful how to analyzed and filter

    plot(prot.means,prot.vars)
    plotDensities(proteins2)
    limma::plotMA(proteins2)
    boxplot(proteins2)

    plotMD(proteins2,column = 60)
    abline(h=0,col="red")

# This is so much different from standard MA plot of expression.
# IDK for now what are the implications. Definitely no two-sample testing.
# Let's just see what median polish does.

    prot.polish <- medpolish(proteins2, eps = 0.01, maxiter = 10,
                             trace.iter = TRUE,na.rm = T)


    plotMD(prot.polish$residuals,column = 60)
    abline(h=0,col="red")



# Let's see what's the relation between proteins and the mRNAs 
    
    mRNA.norm <- expression[intersect(rownames(proteins2), 
                                       rownames(expression)),]
    prot.norm <- proteins[intersect(rownames(proteins2), 
                                    rownames(expression)),]
    
    mRNA.norm <- mRNA.norm[,intersect(colnames(proteins2), 
                                      colnames(expression))]
    prot.norm <- prot.norm[,intersect(colnames(proteins2),
                                      colnames(expression))]
    
   
    
        
 
       


# We compare the variable genes and proteins. Aligning them, we dont have more
# than 300 genes left
nrow(mRNA.norm)
nrow(expression)
nrow(proteins)
ncol(prot.norm)
ncol(mRNA.norm)


prot.mRNA.cors <- sapply(1:nrow(mRNA.norm), function(X){
                         cor(t(mRNA.norm[X,]),
                             t(prot.norm[X,]),
                             use='complete.obs')})



prot.mRNA.pvals <- sapply(1:nrow(mRNA.norm), function(X){
                          zz <- cor.test(t(mRNA.norm[X,]),
                                         t(prot.norm[X,]))
                          return(zz$p.value)})

plotDensities(prot.mRNA.cors)


# The negative correlations are not statistically significant
prot.mRNA.fdr <- p.adjust(prot.mRNA.pvals)
hist(prot.mRNA.cors[prot.mRNA.fdr < 0.05])
plotDensities(prot.mRNA.cors[prot.mRNA.fdr < 0.05])
hist(prot.mRNA.fdr)
length(prot.mRNA.cors[prot.mRNA.fdr < 0.05])



mRNA.fitlered <- mRNA.norm[prot.mRNA.fdr < 0.05 & prot.mRNA.cors > 0.5,]


heatmap(data.matrix(mRNA.fitlered), scale = "none")

# This is not much informative. Lets do some filtering on the expressions


mRNA.norm <- expression2[intersect(rownames(proteins2), rownames(expression2)),]
prot.norm <- proteins[intersect(rownames(proteins2), rownames(expression2)),]

mRNA.norm <- mRNA.norm[,intersect(colnames(proteins2), colnames(expression2))]
prot.norm <- prot.norm[,intersect(colnames(proteins2), colnames(expression2))]



prot.mRNA.cors <- sapply(1:nrow(mRNA.norm), function(X){
    cor(t(mRNA.norm[X,]),
        t(prot.norm[X,]),
        use='complete.obs')})



prot.mRNA.pvals <- sapply(1:nrow(mRNA.norm), function(X){
    zz <- cor.test(t(mRNA.norm[X,]),
                   t(prot.norm[X,]))
    return(zz$p.value)})

plotDensities(prot.mRNA.cors)


# The negative correlations are not statistically significant
prot.mRNA.fdr <- p.adjust(prot.mRNA.pvals)
hist(prot.mRNA.cors[prot.mRNA.fdr < 0.05])
plotDensities(prot.mRNA.cors[prot.mRNA.fdr < 0.05])
hist(prot.mRNA.fdr)
length(prot.mRNA.cors[prot.mRNA.fdr < 0.05])



mRNA.fitlered <- mRNA.norm[prot.mRNA.fdr < 0.05 & prot.mRNA.cors > 0.5,]


heatmap(data.matrix(mRNA.fitlered), scale = "none")

# Now we can emerge four row clusters and four column clusters visually

d <- dist( data.matrix(mRNA.fitlered ))
hc <- hclust(d)
hc
plot(hc,labels=rownames(mRNA.fitlered),cex=0.5)
abline(h = 65 ,col = "red")
hclusters <- cutree(hc, h=65)
table(true=rownames(mRNA.fitlered), cluster=hclusters)

ex.names <- names(which(hclusters == 4))
write.csv(ex.names,"exampleCluster.csv")



## But let's also do some cluster analysis using pca 

## Some PCA plot
seq.pca <- prcomp(t(mRNA.fitlered), center = TRUE, scale. = TRUE)
plot(seq.pca,type = "l")

# The plot tells us after the 5th PC we can't get much additional variance explained

scoreTot <- seq.pca$x[,1:5]
# Let's work on 5 PCs


library(mclust)

BIC <- mclustBIC(scoreTot)
plot(BIC)
summary(BIC)

ICL <- mclustICL(scoreTot)
summary(ICL)


mod <- Mclust(scoreTot, G = 4, modelName = "VII")
summary(mod, parameters = TRUE)

plot(mod, what = "classification", main = FALSE)
boot <- MclustBootstrap(mod, nboot = 999, type = "bs")
summary(boot, what = "se")

modTot <- kmeans(scoreTot, 4)


#
#
#
## Listing immuno proteins from GO
im.tab <- read.csv("select_immunoprotein2", sep = "\t")
head(im.tab)



library(ggbiplot)
set.seed(1)
g <- ggbiplot(seq.pca , obs.scale = 1, var.scale =1 , var.axes = F,  ellipse = T)+
    geom_point(aes( ), size = 3) + theme_bw()




library("pca3d")


pca3d(seq.pca,group = modTot$cluster)
# Parsing and  the proto expression Data


### Just some summary of the types of analyses that I think I should do
#   Without mutation profiles.
#       1. Get the networks of changes.
#               a. Get highly variable genes.
#               b. Get highly variable proteins and contrast.
#               c. Trying to get the networks through different base data
#               d. Run the types of analysis such as WGCNA, SPACEJAM, etc.
#               e. Run enrichment analysis in the clusters of networks.
#                   It would be very nice to land CADIA now.
#               f. Comprate that with the ground truth networks, e.g. STRINGdb.
#       2.





### Some enrichment analysis to be provided


### Some network stuff

library(minet)
library(RBGL)
library(graph)
library(igraph)
mim <- build.mim(t(mRNA.fitlered))
arac.net <- round(minet::aracne(mim, eps=0.3), 3)

arac.net <- arac.net > 0.4

sum(arac.net)
plot(graph.adjacency(arac.net),layout= layout.circle(graph.adjacency(arac.net)))

