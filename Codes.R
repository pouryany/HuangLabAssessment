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


# Let's do some clustering on the mutation matrix and see what comes out
# Why doing clustering? maybe finding subgroups that we can further contrast
# I am not sure where mutations happen together or correlate, but maybe clusters
# together can bring some higher level descriptions.

heatmap(mutated, Colv=F,scale ="none")

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
dim(expression)
# Just inspecting the data distribution. We know that the data is normalized
# So there shouldn't be any problems here. Distributions should align
    plotDensities(expression, legend = F)
    plotMA(expression)
    boxplot(expression)
    plotMD(expression,column = 1)
    abline(h=0,col="grey")



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





## Some PCA plot
    seq.pca <- prcomp(t(expression2), center = TRUE, scale. = TRUE)
    plot(seq.pca,type = "l")


    library(ggbiplot)
    set.seed(1)
    g <- ggbiplot(seq.pca , obs.scale = 1, var.scale =1 , var.axes = F,  ellipse = T,
                  groups = types)+
        geom_point(aes( colour = types ), size = 3) + theme_bw()


library("pca3d")


    pca3d(seq.pca,group = types)
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

# We compare the variable genes and proteins. Aligning them, we dont have more
# than 300 genes left
nrow(mRNA.norm)
nrow(expression2)
nrow(proteins)
ncol(prot.norm)
ncol(mRNA.norm)


prot.mRNA.cors <- sapply(1:nrow(mRNA.norm), function(X){
    cor(t(mRNA.norm[X,]),t(prot.norm[X,]),use='complete.obs')
})


cor(t(mRNA.norm[1,]),(prot.norm[1,]))

prot.mRNA.pvals <- sapply(1:nrow(mRNA.norm), function(X){
     zz <- cor.test(t(mRNA.norm[X,]),t(prot.norm[X,]))
     return(zz$p.value)
})
plotDensities(prot.mRNA.cors)


# The negative correlations are not statistically significant
prot.mRNA.fdr <- p.adjust(prot.mRNA.pvals)
hist(prot.mRNA.cors[prot.mRNA.fdr < 0.005])
hist(prot.mRNA.fdr)
length(prot.mRNA.cors[prot.mRNA.fdr < 0.005])



mRNA.fitlered <- mRNA.norm[prot.mRNA.fdr < 0.05 & prot.mRNA.cors > 0.5,]


heatmap(data.matrix(mRNA.fitlered), scale = "none")

d <- dist( data.matrix(mRNA.fitlered ))
hc <- hclust(d)
hc
plot(hc,labels=rownames(mRNA.fitlered),cex=0.5)
hclusters <- cutree(hc, h=65)
table(true=rownames(mRNA.fitlered), cluster=hclusters)


ex.names <- names(which(hclusters == 4))
write.csv(ex.names,"exampleCluster.csv")

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

