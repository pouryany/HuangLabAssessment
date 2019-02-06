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




expression2 <- expression[,intersect(colnames(expression),
                                     colnames(mutated))]

color.mute <-  rep("blue", ncol(expression2))
color.mute[colnames(expression2) %in%
               names(mutated["TP53",]
                     [mutated["TP53",] ==1])] <- "red"

plotMDS(expression2, col = color.mute)

types <- rep("No_53",ncol(expression))
types[colnames(expression) %in%
          names(which(mutated["TP53",] == 1))] <- "TP53"

# Some filtering
means <- apply(expression2, 1, mean)
vars   <- apply(expression2, 1, var)
plot(means,vars)

expression2 <- expression2[vars !=0,]
dim(expression2)



# Somewhere in (2,3) the second derivative of the density plot/CDF of variance
# becomes positive. let's just pick that as the variance threshold for now.

expression3 <- expression2[vars > 3,]




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

tT.deGenes <- tT[tT$adj.P.Val < 0.05, ]
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

nrow(ego)
sampleGOdata <- new("topGOdata", description = "Simple session",
                     ontology = "BP", allGenes = geneList,
                     geneSel = topDiffGenes, nodeSize = 10,
                     annot = annFUN.db, affyLib = affyLib)

??topGO
# Parsing and  the proto expression Data

proteins <- read.csv("BRCA_PRO_formatted_normalized.txt",
                     stringsAsFactors = F, sep = "\t", row.names = 1)
head(proteins)




# Just a prefiltering step
missings   <- apply(proteins, 1, function(X){
    sum(is.na(X))
})

sum(missings > 0)
# There are 3828 proteins with missing values
# We have to filter and potentially impute...

# For now lets say no more than 30 proteins missing.
proteins  <- proteins[!(missings > 30),]



prot.means <- apply(proteins, 1, mean, na.rm =T)
prot.vars   <- apply(proteins, 1, var,  na.rm =T)

sum(is.na(prot.means))
sum(is.na(prot.vars))
proteins <- proteins[!is.na(prot.means),]


# Just some data inspection. The normalization is described is the handout so
# we should be careful how to analyzed and filter

plot(prot.means,prot.vars)
plotDensities(proteins)
plotMA(proteins)
boxplot(proteins)

plotMD(proteins,column = 60)
abline(h=0,col="grey")

# This is so much different from standard MA plot of expression.
# IDK for now what are the implications. Definitely no two-sample testing.
# Let's just see what median polish does.

prot.polish <- medpolish(proteins, eps = 0.01, maxiter = 10, trace.iter = TRUE,
                         na.rm = T)


plotMD(prot.polish$residuals,column = 60)
abline(h=0,col="grey")


#   proteins <- prot.polish$residuals
#  dim(proteins)

mRNA.norm <- expression2[intersect(rownames(proteins), rownames(expression2)),]
prot.norm <- proteins[intersect(rownames(proteins), rownames(expression2)),]

mRNA.norm <- mRNA.norm[,intersect(colnames(proteins), colnames(expression2))]
prot.norm <- prot.norm[,intersect(colnames(proteins), colnames(expression2))]



