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

# TP53 the usual suspect ...


# Let's do some clustering on the mutation matrix and see what comes out
# Why doing clustering? maybe finding subgroups that we can further contrast
# I am not sure where mutations happen together or correlate, but maybe clusters
# together can bring some higher level descriptions. 

heatmap(mutated, Colv=F,distfun = "manhattan",scale ="none")

?heatmap
d <- dist((mutated) )
clusters <- hclust(d)
plot(clusters,cex=0.5)


# For now let's just investigate the data.
# I will run a scenario by doing differential expression between 
# mutated profiles later

# Parsing and inspecting the mRNA expression Data
expression <- read.csv("BRCA_mRNA_formatted_normalized.txt",
                       stringsAsFactors = F, sep = "\t", row.names = 1)
head(expression)

# Just inspecting the data distribution. We know that the data is normalized 
# So there shouldn't be any problems here. Distributions should align
plotDensities(expression, legend = F)
plotMA(expression)
boxplot(expression)

# Some filtering
means <- apply(expression, 1, mean)
vars   <- apply(expression, 1, var)
plot(means,vars)

expression2 <- expression[vars !=0,]
dim(expression2)

    plot(1:length(vars),vars)
    plotDensities(vars)
    plot(ecdf(vars))

# Somewhere in (2,3) the second derivative of the density plot/CDF of variance
# becomes positive. let's just pick that as the variance threshold for now. 

    expression2 <- expression[vars > 3,]

    # Getting some plots for the filtered expression data.
    plotDensities(expression2, legend = F)
    plotMA(expression2)
    boxplot(expression2)





# Parsing and  the proto expression Data

proteins <- read.csv("BRCA_PRO_formatted_normalized.txt",
                     stringsAsFactors = F, sep = "\t", row.names = 1)
head(proteins)


# Just a prefiltering step
    prot.means <- apply(proteins, 1, mean)
    prot.vars   <- apply(proteins, 1, var)
    
    sum(is.na(prot.means))
    sum(is.na(prot.vars))
    proteins <- proteins[!is.na(prot.means),]

# Just some data inspection. The normalization is described is the handout so 
# we should be careful how to analyzed and filter

    plot(prot.means,prot.vars)
    plotDensities(proteins)
    plotMA(proteins)
    boxplot(proteins)




mRNA.norm <- expression2[intersect(rownames(proteins), rownames(expression2)),]
prot.norm <- proteins[intersect(rownames(proteins), rownames(expression2)),]

mRNA.norm <- mRNA.norm[,intersect(colnames(proteins), colnames(expression2))]
prot.norm <- prot.norm[,intersect(colnames(proteins), colnames(expression2))]

# We compare the variable genes and proteins. Aligning them, we dont have more
# than 300 genes left
nrow(mRNA.norm)

prot.mRNA.cors <- sapply(1:nrow(mRNA.norm), function(X){
    cor(t(mRNA.norm[X,]),t(prot.norm[X,]))
})

prot.mRNA.pvals <- sapply(1:nrow(mRNA.norm), function(X){
     zz <- cor.test(t(mRNA.norm[X,]),t(prot.norm[X,]))
     return(zz$p.value)
})
plotDensities(prot.mRNA.cors)


# The negative correlations are not statistically significant
prot.mRNA.fdr <- p.adjust(prot.mRNA.pvals)
hist(prot.mRNA.cors[prot.mRNA.fdr < 0.005])
hist(prot.mRNA.fdr)


