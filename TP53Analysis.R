# Test case analyzing between TP53 profiles
rm(list = ls())
library(limma)

# Parsing and inspecting the Mutation Data

mutated     <- readRDS("Data/mutationcleaned.rds")
p53.mutated <- names(which(mutated["TP53",] == 1))
no.53       <- names(which(mutated["TP53",] == 0))


proteins    <- readRDS("Data/protcleaned.rds")
expression  <- readRDS("Data/mRNAcleaned.rds")
expression  <- expression[,intersect(colnames(expression),
                                     colnames(mutated))]


color.mute <-  rep("blue", ncol(expression))
color.mute[colnames(expression) %in%
               names(mutated["TP53",]
                     [mutated["TP53",] ==1])] <- "red"

types <- rep("No_53",ncol(expression))
types[colnames(expression) %in%
          names(which(mutated["TP53",] == 1))] <- "TP53"

# Some filtering
means <- apply(expression, 1, mean)
vars   <- apply(expression, 1, var)
expression2 <- expression[vars !=0,]
dim(expression2)


# Differential expression analysis

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


# There are 308 DEG associated with TP53 mutation.



geneList <- (rownames(tT))
deGenes  <- (rownames(tT.deGenes))


gene.df    <- bitr(geneList, fromType = "SYMBOL",
                    toType = c("ENTREZID","ENSEMBL"),
                    OrgDb = org.Hs.eg.db)
deGenes.df <- bitr(deGenes, fromType = "SYMBOL",
                   toType = c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)


set.seed(1)
cadia.res  <- CADIA::causalDisturbance(deGenes.df$ENTREZID,gene.df$ENTREZID,
                                       iter = 5000)


# Just looking at KEGG pathway enrichments on TP53 Mutations.



# Again filtering on proteins
missings   <- apply(proteins, 1, function(X){
    sum(is.na(X))
})

# There are 3628 proteins with missing values

# For now lets say no more than 30 proteins missing.
proteins2   <- proteins
prot.means  <- apply(proteins2, 1, mean, na.rm =T)
prot.vars   <- apply(proteins2, 1, var,  na.rm =T)
proteins2   <- proteins[!(missings > 0),]





# Loading cytokines and immune proteins

imm.tab  <- read.csv("Data/immunoproteins",header = F, sep = "\t")
imm.prot  <- as.character(imm.tab$V2)
imm.prot2 <- intersect(imm.prot, rownames(proteins2)) 

cyt.tab <- read.csv("Data/cytokines",header = F, sep = "\t")
head(cyt.tab)

cyt.prot <- as.character(cyt.tab$V2)

intersect(cyt.prot, rownames(proteins2))  

cyt.imm <- union(cyt.prot,imm.prot)    
sum(is.na(cyt.imm))

cyt.imm2 <- intersect(cyt.imm, rownames(proteins))  




mRNA.norm <- expression2[intersect(cyt.imm2, rownames(tT.deGenes)),]
prot.norm <- proteins2[intersect(cyt.imm2, rownames(tT.deGenes)),]

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





prot.mRNA.fdr <- p.adjust(prot.mRNA.pvals)
plotDensities(prot.mRNA.cors[prot.mRNA.fdr < 0.05])
hist(prot.mRNA.fdr)
length(prot.mRNA.cors[prot.mRNA.fdr < 0.05])



mRNA.fitlered <- mRNA.norm[prot.mRNA.fdr < 0.05,]
prot.fitlered <- prot.norm[prot.mRNA.fdr < 0.05,]

library(RColorBrewer)




library(mclust)

BIC <- mclustBIC(t(mRNA.fitlered))
plot(BIC)
summary(BIC)

BIC <- mclustBIC(t(prot.fitlered))
plot(BIC)
summary(BIC)



mod.mRNA <- Mclust(t(mRNA.fitlered), G = 3, modelName = "VEI")
mod.prot <- Mclust(t(prot.fitlered), G = 3, modelName = "VEI")


library(RColorBrewer)


clust.order <- unlist(tapply(1:ncol((mRNA.fitlered)),
                             as.factor(mod.mRNA$classification),I)
                      ,use.names = F)
mRNA.mat    <- data.matrix(mRNA.fitlered)[,clust.order]



# Just some color setting for visualization

my_group    <- as.numeric(as.factor(mod.mRNA$classification))
my_col1     <- brewer.pal(8, "Set1")[my_group]
my_col1     <- my_col1[clust.order]

my_group    <- as.numeric(as.factor(mod.prot$classification))
my_col2     <- brewer.pal(8, "Set2")[my_group]
my_col2     <- my_col2[clust.order]

coul        <- colorRampPalette(brewer.pal(10, "RdBu"))(50)
my_col      <- cbind("mRNA"= my_col1, "Protein" = my_col2)



# Distance based clustering 

d  <- dist( data.matrix(mRNA.fitlered ), method = "euclidian")
hc <- hclust(d,method = "complete")

plot(hc,labels=rownames(mRNA.fitlered),cex=0.5)

# Three seems suitable so lets just do the rest of coloring and heatmap   
hclusters  <- cutree(hc, h=45)
row.order  <- unlist(tapply(1:nrow((mRNA.fitlered)),
                            as.factor(hclusters),I),use.names = F)

mRNA.mat   <- mRNA.mat[row.order,]

my_group   <- as.numeric(as.factor(hclusters))
my_col1    <- brewer.pal(9, "Set1")[my_group]
my_row_col <- brewer.pal(9, "Set1")[my_group]
my_row_col <- my_row_col[row.order]

heatmap3(mRNA.mat,Colv = NA,Rowv = NA,showColDendro = F,showRowDendro = F,
         scale = "none", col = coul,RowSideColors =  my_row_col,
         ColSideColors=my_col)


# In the three clusters of the samples, there is strong confirmation
 # between clusters 2 of proteins and 2 of mRNAs. 




## Some PCA plot
seq.pca <- prcomp(t(mRNA.fitlered), center = TRUE, scale. = TRUE)
plot(seq.pca,type = "l")

# The plot tells us after the 4th PC we can't get much additional variance explained

scoreTot <- seq.pca$x[,1:4]
# Let's work on 5 PCs


library(ggbiplot)
set.seed(1)
 ggbiplot(seq.pca , obs.scale = 1, var.scale =1 , var.axes = F,
              ellipse = T,groups = as.factor(mod.mRNA$classification))+
    geom_point(aes(shape = types ), size = 3) + theme_bw()


# PCA at this stage does not show the distinction between the subgroups of 
 # mutations 
 
 
 
# Now that we have small number of features and enough samples, 
 # Let's do state of the art link prediction

library(GENIE3)
library(GGally)
library(sna)
 
weightMat <- GENIE3(data.matrix(prot.fitlered))
plotDensities(weightMat)
weightMat[weightMat < 0.1] <- 0
weightMat[weightMat > 0.1] <- 1




apop.network1 <- network::as.network.matrix(weightMat)
set.seed(2)
ggnet2(apop.network1, node.label = colnames(weightMat), 
       arrow.size = 5,label.size = 5,label.trim = T,arrow.gap = 0.01915,
        mode = "fruchtermanreingold", 
       layout.par = list(cell.jitter =0.001, niter = 1000 )) +
 theme(legend.position="none",
          legend.title=element_blank())

ggsave("images/proteinNetworkGENIE3.pdf")

weightMat2 <- GENIE3(data.matrix(mRNA.fitlered))
hist(weightMat2)
weightMat2[weightMat2 < 0.1] <- 0
weightMat2[weightMat2 > 0.1] <- 1



apop.network2 <- network::as.network.matrix(weightMat2)
set.seed(2)
ggnet2(apop.network2, node.label = colnames(weightMat2), 
       arrow.size = 5,label.size = 5,label.trim = T,arrow.gap = 0.01915,
       mode = "fruchtermanreingold", layout.par = list(cell.jitter =0.001,
                                                       niter = 1000 )) +
theme(legend.position="none",
      legend.title=element_blank())

ggsave("images/mRNANetworkGENIE3.pdf")



weightMat3 <- weightMat2 + weightMat 

weightMat3[weightMat3 >0 ] <- 1

string.table <- read.csv("Data/string_interactions.tsv", stringsAsFactors = F, sep = "\t")
string.table <- string.table[,c(1,2)]

for ( i in 1:nrow(string.table)){
     from <- string.table[1,i]
     to   <- string.table[2,i]
     weightMat3[from,to] <- weightMat3[from,to] + 2
     weightMat3[to,from] <- weightMat3[to,from] + 2
}


apop.network3 <- network::as.network.matrix(weightMat3,
                                            names.eval = "weights",
                                            ignore.eval = FALSE)

set.seed(2)
network::set.edge.attribute(apop.network3, "lty", ifelse(apop.network3 %e% "weights" < 2, 1, 2))
ggnet2(apop.network3,   label = TRUE, edge.size = "weights",
       edge.color = "weights",edge.lty = "lty",arrow.gap = 0.01915)



 # Just adding later:  combined predicted networks. 
