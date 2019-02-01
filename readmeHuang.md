## Programming assignment for Huang lab ##
## immuno proteomics and genomics analyses of breast cancer ##

* Contact: Kuan-lin Huang @ MSSM
* kuan-lin.huang@mssm.edu

## Background ##
Diverse immune response and infiltration in tumor microenvironment contributes to varying tumor development across individuals, often having important implications in predicting response to immuno-therapy. Gene expression of immune gene markers have been utilized to de-convolute member profiles of the immune cell types in tumor or other tissues, using tools such as CYBERSORT (https://www.nature.com/articles/nmeth.3337) or xCell (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5688663/). However, the immuno-proteomes in tumors remain poorly characterized. 

Working with the CPTAC consortium, we have recently published comprehensive proteomic profiles of breast cancer (https://www.nature.com/articles/nature18003). In the 77 primary tumor samples passing quality control, we have identified more than 10K proteins. Such data provide us with a unique opportunity to better understand the immune variation through dissecting the tumor proteome.  

## Key questions for you to answer ##
# Please note some of these may be open-ended questions and you may run out of time to address them fully. We'd like to simply access how you think about these problems, break them down into pieces, tackle them with good approaches, and interpret the results putting them into biological contexts. If you do not have complete answers to all of them it is okay. Feel free to leverage existing tool/algorithms or develop your own as you see suitable. If you find fit, you may want to attack questions that you are interested in with higher depth than the others. #
0. Quick exploratory analyses of data quality. 
1. How are the major cytokines and immune proteins expressed differently across tumors?
2. How can we classify tumors based on their immune-related proteins?
[2.1 Optional: leveraging baseline immune cell proteome (https://www.nature.com/articles/ni.3693); how can we deconvolute the tumor proteomes into cell types. How does this behave differently from the immune profiles generated using the same methods using tumor transcriptomes?]
3. How does the immune tumor proteome and inferred cell type differ between mutation carriers and non-carriers in key driver genes?

4. Technical question: 1) How and why do you choose the specific algorithm or statistical model in your analyses? 2) How and why did you define the function or class the way you coded them? 
5. Presentation: imagine you are working on for a conference submission due in two days for a one-page summary of result and figure. Please highlight your key methods and findings in this two-page document, with additional methods and figures in Supplement (if applicable). 

## Description of files in directory:
# note given data availability not all samples have mutation and gene expression data, but the sample names (column names) in these matrixes should match up
1. BRCA_MC3_SOMATIC_formatted_amino_acid.txt Sample-gene matrix of somatic mutations listing out the amino acid change introduced by the mutation.
2. BRCA_mRNA_formatted_normalized.txt        Sample-gene matrix of mRNA. mRNA expression is measured in RSEM and normalized across samples at 75% quantile as downloaded from firehose, and further normalized by log2(RSEM+1).
3. BRCA_PRO_formatted.txt	Sample-gene matrix of protein expression, calculated as log2(sample/reference channel ion intensity) using iTRAQ quantitative proteomics data. Note that between different runs (files) the reference used was different and thus this data is inherently different from gene expression data. The expression ratio further normalized so all samples have mean of 0 and SD of 1 in their proteos. 
