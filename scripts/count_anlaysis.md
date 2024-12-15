---
title: "RNA-Seq Analysis"
output: html_document
author: Mohit Poudel
editor_options: 
  markdown: 
    wrap: 72
---
    R

### Install and Load Necessary Packages
    install.packages("BiocManager")
    BiocManager::install(c('limma','DESeq2','AnnotationDbi','org.At.tair.db','pathview','gage','gageData','GO.db','GOstats'))
    install.packages(c('dplyr','gplots','ggplot2','ggrepel'))

### Load necessary libraries for the analysis
    library(DESeq2) 
    library(dplyr) 
    library(ggplot2)

### Get the current working directory to confirm file path
    getwd()

### Load the count data from a CSV file, assuming space-separated values
    countData = read.csv("/home/labuser/MohitPoudel/extra/rna_seq_analysis/output/merged_counts.csv",sep=" ", header=TRUE)

### View the column names of the loaded data
    colnames(countData)

### Set the row names of countData to the sample ID column
    rownames(countData) = countData$X.sample_id

### Remove the 'X.sample_id' column from countData
    countData = countData[,-1]

### Rename the columns of the data frame to reflect sample names
    colnames(countData) <- c("tr720","tr721","tr722","con723","con725")

### View the first few rows of the count data
    head(countData)

### Remove the last column of the countData, usually to drop unwanted data
    countData <- countData[, -ncol(countData)]

### Check the first few rows of the updated data
    head(countData)

### Get a summary of the count data (basic statistics)
    summary(countData)

### Calculate the total counts for each sample
    colSums(countData)

### Plot the total counts per sample in millions
    par(mar=c(8,4,4,1)+0.1)
    barplot(colSums(countData)/1e6, las=3)

### Plot a histogram of counts for one sample (e.g., tr720)
    hist(countData$tr720, br=100)


## Log Transformation of Counts

### Perform log transformation on the count data (adding 1 to avoid log(0))
    logCount = log2(1 + countData)
    logCount

### Plot log-transformed counts between two samples (tr720 vs tr721)
    plot(logCount[,1], logCount[,2])

### Plot log-transformed counts between tr720 and con725
    plot(logCount[,1], logCount[,5])

## Prepare Metadata for DESeq2

### Define the sample types (treatment or control) for each sample
    sample_type = c("treatment","treatment","treatment","control","control")

### Create a metadata dataframe with sample names and sample types
    colData = as.data.frame(cbind(colnames(countData), sample_type))

### View the metadata
    colData

## DESeq2 Anlaysis

### Create DESeqDataSet object from the count data and metadata, with sample_type as the factor
    dds = DESeqDataSetFromMatrix(countData, colData, design = \~sample_type)

### Run the DESeq2 analysis to estimate size factors, dispersion, and differential expression
    dds0 = DESeq(dds)

### Check the number of rows (genes) after DESeq2 processing
    nrow(dds0)

### Filter out rows (genes) with low counts (only keep genes with counts \> 5 in at least one sample)
    dds = dds0[rowSums(counts(dds0)) \> 5,]

### Check how many rows (genes) are left after filtering
    dim(dds)[1]

### Check the size factors, which are used to normalize the count data
    sizeFactors(dds)

## PCA plot
### Perform regularized log transformation (rlog) for PCA
    rld = rlog(dds)

### Plot PCA with samples grouped by sample_type (treatment vs control)
    plotPCA(rld, intgroup = c("sample_type"))