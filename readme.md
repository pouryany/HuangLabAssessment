Huang Lab Assessment Project
================
Pourya Naderi Yeganeh
2/1/2019

Project overview
----------------

Initial inspection
------------------

We have three sets of files that contain mutations profiles, mRNA expressions, and relative protein expressions. The first phase to insepct and pre-process the data sources. For this, we run some standard procedures. But first, some library loading.

``` r
rm(list = ls())
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



head(mutated[,1:5],1)
```

    ##          A2.A0CM.01A A2.A0D2.01A A2.A0EQ.01A A2.A0EV.01A A2.A0EX.01A
    ## ARHGAP35          wt p.G169Afs*2          wt          wt          wt

``` r
head(proteins[,1:5],3)
```

    ##       C8.A12T.01A A2.A0D2.01A C8.A12U.01A AR.A1AS.01A A2.A0EV.01A
    ## A1BG    0.9087671  0.06445806  -0.1471702  -0.6461039   1.2123155
    ## A2M     0.1346704 -0.11119245   0.3026183  -1.1478470   0.8032429
    ## A2ML1  -3.8240789 -0.65919211  -1.5357392  -1.7804773          NA

``` r
head(expression[,1:5],1)
```

    ##      A2.A0CM.01A A2.A0D2.01A A2.A0EQ.01A A2.A0EV.01A A2.A0EX.01A
    ## A1BG    5.656971    4.703566    4.825658    6.081495    6.911484

Mutation data initial inspection
================================

On the mutations data, we have a charachter matrix with either entries that are 'wt' or specific mutations. Each entry might have a number of mutations and the mutations of a single gene might be happening in different sites. We make a simplifying assumption that a mutation affect the protein level activities in the same manor regardless of number of mutated sites. We changes the mutation matrix into a binary matrix where 1 indicates incidence of a mutation at a specific gene/sample pair and zero otherwise.

``` r
    dim(mutated)
```

    ## [1] 65 73

``` r
    sum(is.na(mutated))
```

    ## [1] 0

``` r
    mutated[mutated != "wt"] <- 1
    mutated[mutated == "wt"] <- 0
    
    mutated <- as.data.frame(mutated)
    mutated <- data.matrix(mutated)
```

TP53 and PI3KCA also the most frequent mutations across the samples. Just for one quick inspection, let's do some clustering on the mutation matrix and see what comes out. Why doing clustering? maybe finding subgroups that we can further contrast. I am not sure where mutations happen together or correlate, but maybe clusters together can bring some higher level descriptions.

``` r
    heatmap(mutated, Colv=F,scale ="none")
```

![](readme_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
    d <- dist((mutated) )
    clusters <- hclust(d)
    plot(clusters,cex=0.5)
```

![](readme_files/figure-markdown_github/unnamed-chunk-3-2.png)

Unsprisingly, we can't get much out this heatmap the way it is. In later analysis we will use subgroups of frequent mutations. For now, we will move on to investigate the expression data.

mRNA expression intitial inspection
===================================

Just inspecting the data distribution. We know that the data is normalized. So, there shouldn't be any problems here. Distributions should align. Just to verify.

``` r
    plotDensities(expression, legend = F)
```

![](readme_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
    boxplot(expression)
```

![](readme_files/figure-markdown_github/unnamed-chunk-4-2.png)

``` r
    plotMD(expression,column = 1)
    abline(h=0,col="red")
```

![](readme_files/figure-markdown_github/unnamed-chunk-4-3.png)

``` r
    means  <- apply(expression, 1, mean)
    vars   <- apply(expression, 1, var)
    sds    <- apply(expression, 1, sd)
    plot(means,vars)
    abline(h=10,col="red")
```

![](readme_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
    expression2 <- expression[vars !=0 & means !=0,]
    dim(expression)
```

    ## [1] 20502    74

``` r
    dim(expression2)
```

    ## [1] 19913    74

We see some reduction becasue of zero variance genes

``` r
   plotDensities(vars)
```

![](readme_files/figure-markdown_github/unnamed-chunk-6-1.png)

Somewhere in (2,3) the second derivative of the density plot/CDF of variance becomes positive.let's just pick that as the variance threshold. Let's do some preliminary analysis of mRNA Expressions
