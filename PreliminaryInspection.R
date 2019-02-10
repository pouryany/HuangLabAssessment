rm(list = ls())
library(limma)
# Parsing and inspecting the Mutation, Expression, and Protein Data

mutated    <- read.csv("Data/BRCA_MC3_SOMATIC_formatted_amino_acid.txt",
                       stringsAsFactors = F, sep = "\t", row.names = 1)
expression <- read.csv("Data/BRCA_mRNA_formatted_normalized.txt",
                       stringsAsFactors = F, sep = "\t", row.names = 1)
proteins   <- read.csv("Data/BRCA_PRO_formatted_normalized.txt",
                       stringsAsFactors = F, sep = "\t", row.names = 1)



head(mutated,5)
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

# TP53 the usual suspect ... PI3KCA  also frequent
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
sds    <- apply(expression, 1, sd)
plot(means,vars)
abline(h=10,col="red")
expression2 <- expression[vars !=0 & means !=0,]
dim(expression)
dim(expression2)

cof.var <- sds/means

# We see some reduction becasue of zero variance genes

plotDensities(vars)
plot(ecdf(vars))





# Moving on proteins    
# Just a prefiltering step on proteins
head(proteins)
missings   <- apply(proteins, 1, function(X){
    sum(is.na(X))
})

sum(missings > 0)
# There are 3628 proteins with missing values
# We have to filter and potentially impute...

# For now lets say no proteins missing from the rest of the analysis.
# Though this assumption can be extended to some level missing
proteins2   <- proteins
prot.means  <- apply(proteins2, 1, mean, na.rm =T)
prot.vars   <- apply(proteins2, 1, var,  na.rm =T)

sum(is.na(prot.means))
sum(is.na(prot.vars))
sum(prot.means == 0)
sum(prot.vars  == 0)

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
proteins3 <- prot.polish$residuals


# It seems that some sort further normalization such as median polish would 
# be useful for the proteins analysis and give a more clear picture.
# However, given my unfamiliarity with the data I just proceed the way it is.

saveRDS(mutated    ,"Data/mutationcleaned.rds")
saveRDS(expression2,"Data/mRNAcleaned.rds")
saveRDS(proteins2  ,"Data/protcleaned.rds")
