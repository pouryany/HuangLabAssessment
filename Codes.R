# Parsing and inspecting the Mutation Data

mutated <- read.csv("BRCA_MC3_SOMATIC_formatted_amino_acid.txt",
                    stringsAsFactors = F, sep = "\t", row.names = 1)

dim(mutated)
sum(is.na(mutated))

mutated[mutated != "wt"] <- 1
mutated[mutated == "wt"] <- 0

mutated <- as.data.frame(mutated)
mutated <- as.matrix(mutated)


# Parsing and inspecting the mRNA expression Data
expression <- read.csv("BRCA_mRNA_formatted_normalized.txt",
                       stringsAsFactors = F, sep = "\t", row.names = 1)
head(expression)

means <- apply(expression, 1, mean)
vars   <- apply(expression, 1, var)
plot(means,vars)


expression <- expression[vars !=0,]

plot(1:length(sds),sds)
plotDensities(sds)
plot(ecdf(vars))

# Getting some plots for the expression data.
library(limma)
plotDensities(expression)
plotMA(expression)
boxplot(expression)


dim(expression)



# Parsing and  the proto expression Data

proteins <- read.csv("BRCA_PRO_formatted_normalized.txt",
                     stringsAsFactors = F, sep = "\t", row.names = 1)
head(proteins)

prot.means <- apply(proteins, 1, mean)
prot.vars   <- apply(proteins, 1, var)

sum(is.na(prot.means))
sum(is.na(prot.vars))


proteins <- proteins[!is.na(prot.means),]

plot(prot.means,prot.vars)


plotDensities(proteins)
plotMA(proteins)
boxplot(proteins)




mRNA.norm <- expression[intersect(rownames(proteins), rownames(expression)),]
prot.norm <- proteins[intersect(rownames(proteins), rownames(expression)),]

mRNA.norm <- mRNA.norm[,intersect(colnames(proteins), colnames(expression))]
prot.norm <- prot.norm[,intersect(colnames(proteins), colnames(expression))]


prot.mRNA.cors <- sapply(1:nrow(mRNA.norm), function(X){
    cor(t(mRNA.norm[X,]),t(prot.norm[X,]))
})

prot.mRNA.pvals <- sapply(1:nrow(mRNA.norm), function(X){
     zz <- cor.test(t(mRNA.norm[X,]),t(prot.norm[X,]))
     return(zz$p.value)
})
plotDensities(prot.mRNA.cors)



prot.mRNA.fdr <- p.adjust(prot.mRNA.pvals)
hist(prot.mRNA.cors[prot.mRNA.fdr < 0.005])
hist(prot.mRNA.fdr)


