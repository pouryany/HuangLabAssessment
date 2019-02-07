# Test case analyzing between TP53 profiles
library(limma)

# Parsing and inspecting the Mutation Data

mutated <- read.csv("BRCA_MC3_SOMATIC_formatted_amino_acid.txt",
                    stringsAsFactors = F, sep = "\t", row.names = 1)


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
ncol(mutated)
# TP53 the usual suspect ...



p53.mutated <- names(which(mutated["TP53",] == 1))
no.53       <- names(which(mutated["TP53",] == 0))



expression <- read.csv("BRCA_mRNA_formatted_normalized.txt",
                       stringsAsFactors = F, sep = "\t", row.names = 1)

head(expression)
dim(expression)




expression <- expression[,intersect(colnames(expression),
                                     colnames(mutated))]


color.mute <-  rep("blue", ncol(expression))
color.mute[colnames(expression) %in%
               names(mutated["TP53",]
                     [mutated["TP53",] ==1])] <- "red"

plotMDS(expression, col = color.mute)

types <- rep("No_53",ncol(expression))
types[colnames(expression) %in%
          names(which(mutated["TP53",] == 1))] <- "TP53"

# Some filtering
means <- apply(expression, 1, mean)
vars   <- apply(expression, 1, var)
plot(means,vars)

expression2 <- expression[vars !=0,]
dim(expression2)



# Somewhere in (2,3) the second derivative of the density plot/CDF of variance
# becomes positive. let's just pick that as the variance threshold for now.

#expression2 <- expression[vars > 3,]




# Doing some differential expression analysis to see whether there's any 
# differences caused by TP53 mutations




library(limma)
library(KEGGgraph)


sml <- c()

# set up the data and proceed with analysis
sml     <- types    # set group names
fl      <- as.factor(sml)

design  <- model.matrix(~ fl + 0, expression2)
colnames(design) <- levels(fl)
fit     <- lmFit(expression2, design)
cont.matrix <- makeContrasts(TP53-No_53, levels=design)
fit2    <- contrasts.fit(fit, cont.matrix)
fit2    <- eBayes(fit2, 0.01)
tT      <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

head(tT)

tT.deGenes <- tT[tT$adj.P.Val < 0.005, ]
tT.deGenes <- tT.deGenes[abs(tT.deGenes$logFC) >1,]

nrow(tT.deGenes)


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.8")



geneList <- (rownames(tT))

deGenes  <- (rownames(tT.deGenes))

library(clusterProfiler)
library(org.Hs.eg.db)

gene.df    <- bitr(geneList, fromType = "SYMBOL",
                    toType = c("ENTREZID","ENSEMBL"),
                    OrgDb = org.Hs.eg.db)
deGenes.df <- bitr(deGenes, fromType = "SYMBOL",
                   toType = c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)

?enrichGO
ego <- enrichGO(gene          = deGenes.df$ENSEMBL,
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

head(ego,50)



set.seed(1)
cadia.res  <- CADIA::causalDisturbance(deGenes.df$ENTREZID,gene.df$ENTREZID,
                                       iter = 5000)






head(proteins)
missings   <- apply(proteins, 1, function(X){
    sum(is.na(X))
})

sum(missings > 0)
# There are 3628 proteins with missing values
# We have to filter and potentially impute...

# For now lets say no more than 30 proteins missing.
proteins2  <- proteins
prot.means <- apply(proteins2, 1, mean, na.rm =T)
prot.vars   <- apply(proteins2, 1, var,  na.rm =T)

sum(is.na(prot.means))
sum(is.na(prot.vars))


proteins2  <- proteins[!(missings > 0),]






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
                         trace.iter = TRUE,na.rm = F)


plotMD(prot.polish$residuals,column = 60)
boxplot(prot.polish$residuals)
abline(h=0,col="red")
proteins3 <- proteins2

# Loading cytokines and immune proteins

imm.tab <- read.csv("immunoproteins",header = F, sep = "\t")
head(imm.tab)

imm.prot  <- as.character(imm.tab$V2)

imm.prot2 <- intersect(imm.prot, rownames(proteins2)) 

cyt.tab <- read.csv("cytokines",header = F, sep = "\t")
head(cyt.tab)

cyt.prot <- as.character(cyt.tab$V2)

intersect(cyt.prot, rownames(proteins2))  

cyt.imm <- union(cyt.prot,imm.prot)    
sum(is.na(cyt.imm))

cyt.imm2 <- intersect(cyt.imm, rownames(proteins3))  




mRNA.norm <- expression2[intersect(cyt.imm2, rownames(tT.deGenes)),]
prot.norm <- proteins3[intersect(cyt.imm2, rownames(tT.deGenes)),]

mRNA.norm <- mRNA.norm[,intersect(colnames(proteins3), 
                                  colnames(expression2))]
prot.norm <- prot.norm[,intersect(colnames(proteins3),
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





prot.mRNA.fdr <- p.adjust(prot.mRNA.pvals)
plotDensities(prot.mRNA.cors[prot.mRNA.fdr < 0.05])
hist(prot.mRNA.fdr)
length(prot.mRNA.cors[prot.mRNA.fdr < 0.05])



mRNA.fitlered <- mRNA.norm[prot.mRNA.fdr < 0.05,]
prot.fitlered <- prot.norm[prot.mRNA.fdr < 0.05,]

library(RColorBrewer)
coul = colorRampPalette(brewer.pal(10, "RdBu"))(25)

my_group=as.numeric(as.factor(types))
my_col=brewer.pal(9, "Set1")[my_group]




hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")

# perform clustering on rows and columns
cl.row <- hclustfunc(distfunc(data.matrix(mRNA.fitlered)))
cl.col <- hclustfunc(distfunc(t(data.matrix(mRNA.fitlered))))

cl.col$labels %in%  names(mutated["TP53",]
                          [mutated["TP53",] ==1])

library("gplots")

heatmap.2(data.matrix(prot.fitlered), scale = "none", col = coul,
        ColSideColors=my_col)

## Some PCA plot
seq.pca <- prcomp(t(mRNA.fitlered), center = TRUE, scale. = TRUE)
plot(seq.pca,type = "l")

# The plot tells us after the 5th PC we can't get much additional variance explained

scoreTot <- seq.pca$x[,1:5]
# Let's work on 5 PCs


library(mclust)

BIC <- mclustBIC(t(mRNA.fitlered))
plot(BIC)
summary(BIC)

ICL <- mclustICL(t(mRNA.fitlered))
summary(ICL)


mod <- Mclust(scoreTot, G = 3, modelName = "VEI")
summary(mod, parameters = TRUE)

plot(mod, what = "classification", main = FALSE)
boot <- MclustBootstrap(mod, nboot = 999, type = "bs")
summary(boot, what = "se")


modTot <- kmeans(scoreTot, 3)



library(ggbiplot)
set.seed(1)
g <- ggbiplot(seq.pca , obs.scale = 1, var.scale =1 , var.axes = F,
              ellipse = T,groups = as.factor(modTot$cluster))+
    geom_point(aes(shape = types ), size = 3) + theme_bw()




library(GENIE3)

weightMat <- GENIE3(data.matrix(prot.fitlered))
hist(weightMat)
weightMat[weightMat < 0.1] <- 0
weightMat[weightMat > 0.1] <- 1


plot(graph.adjacency(weightMat),layout= layout.circle(graph.adjacency(arac.net)))

library(GGally)
apop.network <- network::as.network.matrix(weightMat)
set.seed(2)
ras.plot <- ggnet2(apop.network, node.label = colnames(weightMat), arrow.size = 5,label.size = 5,label.trim = T,arrow.gap = 0.01915,
                   mode = "fruchtermanreingold", layout.par = list(cell.jitter =0.001, niter = 1000 )) +
    theme(legend.position="none",
          legend.title=element_blank())

??igraph::plot

