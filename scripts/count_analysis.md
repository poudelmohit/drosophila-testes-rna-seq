### Install Packages


```R
install.packages("BiocManager")
BiocManager::install(c('limma','DESeq2','AnnotationDbi','org.At.tair.db','pathview','gage','gageData','GO.db','GOstats'))
install.packages(c("dplyr","gplots","ggplot2","ggrepel"))
```

### Load Libraries


```R
library(DESeq2) 
library(dplyr) 
library(ggplot2)
```


```R
getwd()
```


'/home/labuser/MohitPoudel/extra/rna_seq_analysis/scripts'



```R

```

### Load Count Data and Preprocess:


```R
count_file_path = file.path(dirname(getwd()),"output/counts/merged_counts.csv")
count_file_path
```


'/home/labuser/MohitPoudel/extra/rna_seq_analysis/output/counts/merged_counts.csv'



```R
countData = read.csv(count_file_path, sep=" ",header=TRUE)
```


```R
head(countData)
```


<table class="dataframe">
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>id</th><th scope=col>tr721</th><th scope=col>tr722</th><th scope=col>con723</th><th scope=col>con725</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>FBgn0085804</td><td>   0</td><td>  0</td><td>  0</td><td>  0</td></tr>
	<tr><th scope=row>2</th><td>FBgn0267431</td><td> 104</td><td>131</td><td> 77</td><td> 57</td></tr>
	<tr><th scope=row>3</th><td>FBgn0039987</td><td>   0</td><td>  0</td><td>  0</td><td>  0</td></tr>
	<tr><th scope=row>4</th><td>FBgn0058182</td><td>   0</td><td>  0</td><td>  0</td><td>  0</td></tr>
	<tr><th scope=row>5</th><td>FBgn0267430</td><td>1045</td><td>992</td><td>700</td><td>729</td></tr>
	<tr><th scope=row>6</th><td>FBgn0266747</td><td> 132</td><td>134</td><td> 82</td><td>126</td></tr>
</tbody>
</table>




```R
colnames(countData)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'id'</li><li>'tr721'</li><li>'tr722'</li><li>'con723'</li><li>'con725'</li></ol>




```R
rownames(countData) = countData$id
```


```R
countData=countData[,-1]
```


```R
summary(countData)
# difference is already observed between control and treatment samples
```


         tr721            tr722            con723             con725        
     Min.   :     0   Min.   :     0   Min.   :     0.0   Min.   :     0.0  
     1st Qu.:     8   1st Qu.:     8   1st Qu.:     6.0   1st Qu.:     8.0  
     Median :   183   Median :   206   Median :   141.0   Median :   176.0  
     Mean   :  2454   Mean   :  2836   Mean   :  2175.1   Mean   :  2110.7  
     3rd Qu.:  1160   3rd Qu.:  1408   3rd Qu.:   852.5   3rd Qu.:   978.5  
     Max.   :728263   Max.   :526754   Max.   :987462.0   Max.   :782962.0  



```R
colSums(countData)
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>tr721</dt><dd>43096987</dd><dt>tr722</dt><dd>49796402</dd><dt>con723</dt><dd>38193053</dd><dt>con725</dt><dd>37060965</dd></dl>




```R
# Open a png device in the desired directory (relative path)
png("../output/plots/feature_counts.png")

# Set margins
par(mar=c(8,4,4,1) + 0.5)

# Create the barplot
barplot(colSums(countData)/1e6, 
        col="skyblue",
        main="Feature Counts (in Millions)", 
        ylab="Read Counts (Millions)", 
        xlab="Samples",
        cex.axis=1.5,   
        cex.lab=1.5,    
        cex.names=1.5,  
        ylim=c(0, 50))   # Set y-axis limit to 50

# Close the device to finalize the Plot
dev.off()


```


<strong>pdf:</strong> 2



```R
# Histogram to see the value distribution for each samples

## Open a png device in the desired directory (relative path)
png("../output/plots/value_dis.png")

par(mfrow=c(2, 2))

## Create histograms for each sample
hist(countData$tr721, br=10, main="Sample tr721", xlab="Counts", xlim=c(0, 500000))
hist(countData$tr722, br=10, main="Sample tr722", xlab="Counts", xlim=c(0, 500000))
hist(countData$con723, br=10, main="Sample con723", xlab="Counts", xlim=c(0, 500000))
hist(countData$con725, bre=10, main="Sample con725", xlab="Counts", xlim=c(0, 500000))

# Reset the plotting layout to default
par(mfrow=c(1, 1))

# Close the device to finalize the Plot
dev.off()

```


<strong>pdf:</strong> 2



```R
logCount = log2(1 + countData)
# logCount

```


```R
# Histogram to see the value distribution for each samples

png("../output/plots/log_dis.png")
par(mfrow=c(2, 2))

# Create histograms for each sample
plot(logCount[,1], logCount[,2])
plot(logCount[,1], logCount[,3])
plot(logCount[,1], logCount[,4])
plot(logCount[,2], logCount[,3])

# Reset the plotting layout to default
par(mfrow=c(1, 1))

dev.off()
```


<strong>pdf:</strong> 2


### Loading in DEseq 


```R
sample_type = c("treatment","treatment","control","control")
```


```R
colData = as.data.frame(cbind(colnames(countData), sample_type))
```


```R
colData
```


<table class="dataframe">
<caption>A data.frame: 4 × 2</caption>
<thead>
	<tr><th scope=col>V1</th><th scope=col>sample_type</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>tr721 </td><td>treatment</td></tr>
	<tr><td>tr722 </td><td>treatment</td></tr>
	<tr><td>con723</td><td>control  </td></tr>
	<tr><td>con725</td><td>control  </td></tr>
</tbody>
</table>




```R
dds = DESeqDataSetFromMatrix(countData, colData, design = ~sample_type)
```

    Warning message in DESeqDataSet(se, design = design, ignoreRank):
    “some variables in design formula are characters, converting to factors”



```R
dds0 = DESeq(dds) # model fitting
```

    estimating size factors
    
    estimating dispersions
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    



```R
nrow(dds0)  # no. of genes in the resulting object
dds = dds0[rowSums(counts(dds0)) > 5,] # filter out genes with too-low count
```


17559



```R
dim(dds)[1] 
```


14725



```R
sizeFactors(dds)
# DESeq2 calculates size factors by using a method that is based on the median of ratios approach, which works as follows:
# For each sample, it computes the ratio of the read count of each gene to the geometric mean of that gene's count across all samples.
# These ratios are then used to calculate a size factor for each sample. This factor is the median of these ratios for each sample, \n 
# ensuring that the gene expression levels are scaled to a common baseline.
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>tr721</dt><dd>1.07188001150061</dd><dt>tr722</dt><dd>1.27907671458853</dd><dt>con723</dt><dd>0.821044902067994</dd><dt>con725</dt><dd>0.903500597477744</dd></dl>



## PCA Plot


```R
rld = rlog(dds) # obtaining regularized log counts
```


```R
pcaData <- plotPCA(rld, intgroup = c("sample_type"), returnData = TRUE)
pcaData
```

    using ntop=500 top features by variance
    



<table class="dataframe">
<caption>A data.frame: 4 × 5</caption>
<thead>
	<tr><th></th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>group</th><th scope=col>sample_type</th><th scope=col>name</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>tr721</th><td>-11.03007</td><td> 1.1572897</td><td>treatment</td><td>treatment</td><td>tr721 </td></tr>
	<tr><th scope=row>tr722</th><td>-20.04683</td><td> 0.0549765</td><td>treatment</td><td>treatment</td><td>tr722 </td></tr>
	<tr><th scope=row>con723</th><td> 12.15737</td><td>-5.4424311</td><td>control  </td><td>control  </td><td>con723</td></tr>
	<tr><th scope=row>con725</th><td> 18.91953</td><td> 4.2301649</td><td>control  </td><td>control  </td><td>con725</td></tr>
</tbody>
</table>




```R
# Extract PCA data
pcaData <- plotPCA(rld, intgroup = c("sample_type"), returnData = TRUE)

# Add sample names to the PCA data
pcaData$Sample <- rownames(pcaData)

# Plot PCA with publication-quality style and additional space for points
library(ggplot2)

png("../output/plots/pca.png")

ggplot(pcaData, aes(x = PC1, y = PC2, color = sample_type, label = Sample)) +
  geom_point(size = 4, alpha = 0.7, shape = 16) +  # Larger, semi-transparent points
  geom_text(vjust = -0.5, size = 5, fontface = "italic", color = "black") +  # Larger labels with italic style
  labs(title = "Principal Component Analysis", x = "PC1", y = "PC2") + 
  scale_color_brewer(palette = "Set1") +  # Better color palette for color-blind accessibility
  theme_minimal(base_size = 16) +  # Minimal theme with larger base font size
  theme(
    axis.title = element_text(size = 18, face = "bold"),   # Bold axis titles
    axis.text = element_text(size = 14),    # Larger axis text
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Bold and centered title
    legend.title = element_text(size = 14, face = "bold"), # Bold legend title
    legend.text = element_text(size = 12),   # Larger legend text
    legend.key.size = unit(1.2, "cm"),  # Adjust legend key size
    plot.margin = margin(1, 1, 1, 2, "cm")  # Increase right margin for extra space
  ) +
  coord_cartesian(clip = "off")  # Allow points outside the default plot area

dev.off()

```

    using ntop=500 top features by variance
    



<strong>pdf:</strong> 2


### Custom Functions:


```R
detectGroups=function(x){
    term=gsub("[0-9]*$","",x)
    # term=gsub("_$","",term)
    # term=gsub("_Rep$","",term)
    return (term)}
# The function returns the category name (control or treatment) from the sample id
```


```R
detectGroups("con721")

```


'con'



```R
dist2=function(x){
    as.dist(1-cor(t(x), method="pearson"))}

# calculates the distance matrix for the heatmap based on Pearson correlation
```


```R
hclust2 = function(x, method="average"){
    hclust(x, method=method)
}
# returns the hierarchical clustering result.

```

#### Heatmap Code:


```R
library(gplots)
```

    
    Attaching package: ‘gplots’
    
    
    The following object is masked from ‘package:IRanges’:
    
        space
    
    
    The following object is masked from ‘package:S4Vectors’:
    
        space
    
    
    The following object is masked from ‘package:stats’:
    
        lowess
    
    



```R

n = 100 # top 100 genes with the highest variance selected
x = assay(rld) # extracts the expression data from rlog-transformed gene expression data
if(n > dim(x)[1]) n = dim(x)[1] # ensureing n does not exceed the number of rows
x = x[order(apply(x, 1, sd), decreasing = TRUE),] # orders the genes (rows) based on their standard deviation

```

#### Preprocessing:


```R
x = as.matrix(x[1:n,]) - apply(x[1:n,], 1, mean) 
#  Normalizes the data by centering each gene (row) to have a mean of zero.
cutoff = median(unlist(x)) + 2 * sd(unlist(x)) 
# Calculates a cutoff value to limit extreme values (outliers). 
# This helps in controlling large values that may skew the color scale.
x[x > cutoff] = cutoff # Caps the values greater than the cutoff to the cutoff value.

```

#### Grouping and Color Mapping:


```R
groups = detectGroups(colnames(x))
groups.colors = rainbow(length(unique(groups)))
# creats color palette with distinct colors for each unique group extracted using detectGroups function
```

### Heatmap Creation


```R
png("../output/plots/heatmap.png")

lmat = rbind(c(5, 4), c(0, 1), c(3, 2)) # defines layout
lwid = c(1.5, 4)                        # width of plot
lhei = c(1, 0.2, 4)                       # height of plot

heatmap.2(x, distfun = dist2, hclustfun = hclust2,
          col = greenred(75), density.info = "none", trace = "none",
          scale = "none", keysize = 0.2, key = T,
          symkey = F, ColSideColors = groups.colors[as.factor(groups)],
          margins = c(12, 12),              # Increase margins for better label visibility
          cexRow = 1.5,                    # Increase font size of y-axis labels
          cexCol = 1.5,                    # Increase font size of x-axis labels
          srtCol = 45,                     # Rotate x-axis labels
          lmat = lmat, lwid = lwid, lhei = lhei)

dev.off()

```


<strong>pdf:</strong> 2


## DESeq2 SDEG Analysis:


```R
resultsNames(dds)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Intercept'</li><li>'sample_type_treatment_vs_control'</li></ol>




```R
res=results(dds)
head(res)
```


    log2 fold change (MLE): sample type treatment vs control 
    Wald test p-value: sample type treatment vs control 
    DataFrame with 6 rows and 6 columns
                  baseMean log2FoldChange     lfcSE      stat    pvalue      padj
                 <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
    FBgn0267431   89.07857      0.3484648  0.476756  0.730908 0.4648351  0.805881
    FBgn0267430  852.47894      0.0766696  0.248772  0.308192 0.7579362  0.935481
    FBgn0266747  116.81037     -0.0726535  0.428796 -0.169436 0.8654538  0.965876
    FBgn0086917    1.62885     -1.3409422  3.145707 -0.426277 0.6699061        NA
    FBgn0010247 1939.39505      0.5853793  0.266768  2.194341 0.0282109  0.199505
    FBgn0086378  675.35212      0.2385825  0.301530  0.791240 0.4288041  0.783478



```R
summary(res)
```

    
    out of 14725 with nonzero total read count
    adjusted p-value < 0.1
    LFC > 0 (up)       : 436, 3%
    LFC < 0 (down)     : 927, 6.3%
    outliers [1]       : 0, 0%
    low counts [2]     : 857, 5.8%
    (mean count < 4)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    


### MA plot


```R
png("../output/plots/ma_plot.png")

res=results(dds, lfcThreshold=0.01)
DESeq2::plotMA(res,ylim=c(-5,5))

dev.off()
```


<strong>pdf:</strong> 2


### Volcano Plot:


```R
png("../output/plots/volcano_plot.png")

# Filter the dataset to remove rows with NA values in log2FoldChange or padj
res1 <- as.data.frame(res) %>%
  filter(!is.na(log2FoldChange) & !is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj < 0.1 & abs(log2FoldChange) >= 1.0 ~ "FDR<0.1 & |Log2FC|≥1",
      TRUE ~ "Not Sig"
    )
  )

# Create the volcano plot
ggplot(res1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), size = 2, alpha = 0.8) +  # Use points with transparency
  scale_color_manual(
    values = c("FDR<0.1 & |Log2FC|≥1" = "red", "Not Sig" = "black"),
    name = "Significance"
  ) +
  theme_minimal(base_size = 14) +  # Clean theme with larger font size
  labs(
    title = "Volcano Plot",
    x = expression(Log[2] ~ Fold ~ Change),
    y = expression(-Log[10] ~ Adjusted ~ p ~ Value)
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    legend.position = "top",  # Move legend to top
    axis.text = element_text(size = 12),  # Larger axis text size
    axis.title = element_text(size = 14, face = "bold")  # Larger axis titles
  ) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "blue", linewidth = 0.7) +  # Significance threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", linewidth = 0.7)  # Fold-change threshold lines

dev.off()

```


<strong>pdf:</strong> 2


### Normalized Count Plot for the Top Gene


```R
png("../output/plots/top_gene_norm_count_plot.png")

res=res[order(abs(res$log2FoldChange),decreasing=TRUE),]
topGene=rownames(res)[1]
plotCounts(dds, gene=topGene, intgroup=c("sample_type"))

dev.off()
```


<strong>pdf:</strong> 2


## GO ANALYSIS:


```R
library(AnnotationDbi)
BiocManager::install("org.Dm.eg.db") # For Drosophila melanogaster
library(org.Dm.eg.db)  

```

    
    Attaching package: ‘AnnotationDbi’
    
    
    The following object is masked from ‘package:dplyr’:
    
        select
    
    
    'getOption("repos")' replaces Bioconductor standard repositories, see
    'help("repositories", package = "BiocManager")' for details.
    Replacement repositories:
        CRAN: https://cran.r-project.org
    
    Bioconductor version 3.20 (BiocManager 1.30.25), R 4.4.2 (2024-10-31)
    
    Warning message:
    “package(s) not installed when version(s) same as or greater than current; use
      `force = TRUE` to re-install: 'org.Dm.eg.db'”
    Old packages: 'BH', 'Biostrings', 'gower', 'openssl', 'parallelly', 'pillar',
      'survival', 'textshaping'
    
    
    



```R
library(GO.db)
library(GOstats)
library(gage)
```

    
    
    Loading required package: Category
    
    Loading required package: Matrix
    
    
    Attaching package: ‘Matrix’
    
    
    The following object is masked from ‘package:S4Vectors’:
    
        expand
    
    
    Loading required package: graph
    
    
    Attaching package: ‘GOstats’
    
    
    The following object is masked from ‘package:AnnotationDbi’:
    
        makeGOGraph
    
    



```R
columns(org.Dm.eg.db)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'ACCNUM'</li><li>'ALIAS'</li><li>'ENSEMBL'</li><li>'ENSEMBLPROT'</li><li>'ENSEMBLTRANS'</li><li>'ENTREZID'</li><li>'ENZYME'</li><li>'EVIDENCE'</li><li>'EVIDENCEALL'</li><li>'FLYBASE'</li><li>'FLYBASECG'</li><li>'FLYBASEPROT'</li><li>'GENENAME'</li><li>'GENETYPE'</li><li>'GO'</li><li>'GOALL'</li><li>'MAP'</li><li>'ONTOLOGY'</li><li>'ONTOLOGYALL'</li><li>'PATH'</li><li>'PMID'</li><li>'REFSEQ'</li><li>'SYMBOL'</li><li>'UNIPROT'</li></ol>




```R
rownames(head(res))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'FBgn0263407'</li><li>'FBgn0261399'</li><li>'FBgn0053300'</li><li>'FBgn0266070'</li><li>'FBgn0267663'</li><li>'FBgn0035845'</li></ol>




```R

```


```R
# Map FBgn IDs to Ensembl IDs
res$ENTREZ = mapIds(org.Dm.eg.db,
                     key = rownames(res),  # Assuming rownames(res) are FBgn IDs
                     column = "ENTREZID",   # Specify the Ensembl column
                     keytype = "FLYBASE",  # Specify the keytype as FLYBASE
                     multiVals = "first")  # Handle multiple matches by taking the first match

```

    'select()' returned 1:1 mapping between keys and columns
    



```R
res$SYMBOL = mapIds(org.Dm.eg.db,
                    key=rownames(res),
                    column="SYMBOL",
                    keytype="FLYBASE",
                    multiVals="first")
```

    'select()' returned 1:1 mapping between keys and columns
    



```R
head(res)
```


    log2 fold change (MLE): sample type treatment vs control 
    Wald test p-value: sample type treatment vs control 
    DataFrame with 6 rows and 8 columns
                 baseMean log2FoldChange     lfcSE      stat      pvalue
                <numeric>      <numeric> <numeric> <numeric>   <numeric>
    FBgn0263407  186.4330      -11.21204   1.50016  -7.47392 7.79388e-14
    FBgn0261399  235.7716       11.10903   1.48313   7.49027 6.88198e-14
    FBgn0053300  119.7278      -10.57336   1.50571  -7.02219 2.18661e-12
    FBgn0266070  107.0638        9.96861   1.51273   6.58982 4.40789e-11
    FBgn0267663   43.2547       -9.10615   1.58908  -5.73044 1.00238e-08
    FBgn0035845 1283.9707       -8.93338   1.15271  -7.74986 9.22022e-15
                       padj      ENTREZ         SYMBOL
                  <numeric> <character>    <character>
    FBgn0263407 9.82596e-12    12798501 lncRNA:CR43453
    FBgn0261399 8.82015e-12     5740113         Pp1-Y1
    FBgn0053300 2.19739e-10     2768931         Muc30E
    FBgn0266070 3.82053e-09    19834998 lncRNA:CR44823
    FBgn0267663 6.20578e-07    26067322 lncRNA:CR46001
    FBgn0035845 1.33194e-12       38907        CG13675


## SDEG hypergeometric test


```R
res_05=as.data.frame(subset(res,padj<0.05))
head(res_05)
```


<table class="dataframe">
<caption>A data.frame: 6 × 8</caption>
<thead>
	<tr><th></th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th><th scope=col>ENTREZ</th><th scope=col>SYMBOL</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>FBgn0263407</th><td> 186.43302</td><td>-11.212043</td><td>1.500156</td><td>-7.473920</td><td>7.793883e-14</td><td>9.825961e-12</td><td>12798501</td><td>lncRNA:CR43453</td></tr>
	<tr><th scope=row>FBgn0261399</th><td> 235.77158</td><td> 11.109027</td><td>1.483127</td><td> 7.490273</td><td>6.881975e-14</td><td>8.820147e-12</td><td>5740113 </td><td>Pp1-Y1        </td></tr>
	<tr><th scope=row>FBgn0053300</th><td> 119.72783</td><td>-10.573362</td><td>1.505707</td><td>-7.022189</td><td>2.186612e-12</td><td>2.197387e-10</td><td>2768931 </td><td>Muc30E        </td></tr>
	<tr><th scope=row>FBgn0266070</th><td> 107.06385</td><td>  9.968609</td><td>1.512729</td><td> 6.589820</td><td>4.407885e-11</td><td>3.820535e-09</td><td>19834998</td><td>lncRNA:CR44823</td></tr>
	<tr><th scope=row>FBgn0267663</th><td>  43.25472</td><td> -9.106148</td><td>1.589084</td><td>-5.730440</td><td>1.002376e-08</td><td>6.205783e-07</td><td>26067322</td><td>lncRNA:CR46001</td></tr>
	<tr><th scope=row>FBgn0035845</th><td>1283.97072</td><td> -8.933378</td><td>1.152714</td><td>-7.749864</td><td>9.220217e-15</td><td>1.331937e-12</td><td>38907   </td><td>CG13675       </td></tr>
</tbody>
</table>




```R
sig_lfc=1
selectGenesUp = unique(res_05[res_05$log2FoldChange > sig_lfc,"ENTREZ"])
selectGenesDn = unique(res_05[res_05$log2FoldChange < (-sig_lfc),"ENTREZ"])
universeGenes = unique(res_05$ENTREZ)


```


```R
dim(res_05)[1]
range(res_05$log2FoldChange)

```


1081



<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>-11.2120429505963</li><li>11.109026668946</li></ol>




```R
length(selectGenesUp)
length(selectGenesDn)
```


152



679


### GO Biological Process:


```R
upParams = new("GOHyperGParams",
               geneIds = selectGenesUp,
               universeGeneIds = universeGenes,
               annotation = "org.Dm.eg.db",
               ontology = "BP",
               pvalueCutoff = 0.01,
               conditional = FALSE,
               testDirection = "over")

dnParams = new("GOHyperGParams",
               geneIds = selectGenesDn,
               universeGeneIds = universeGenes,
               annotation = "org.Dm.eg.db",
               ontology = "BP",
               pvalueCutoff = 0.01,
               conditional = FALSE,
               testDirection = "over")

upBp = hyperGTest(upParams)
dnBp = hyperGTest(dnParams)
```


```R


```


```R
summary(upBp)[1:10,]
```


<table class="dataframe">
<caption>A data.frame: 10 × 7</caption>
<thead>
	<tr><th></th><th scope=col>GOBPID</th><th scope=col>Pvalue</th><th scope=col>OddsRatio</th><th scope=col>ExpCount</th><th scope=col>Count</th><th scope=col>Size</th><th scope=col>Term</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>GO:0090304</td><td>0.0001087850</td><td> 2.872810</td><td>13.058405</td><td>26</td><td>103</td><td>nucleic acid metabolic process                  </td></tr>
	<tr><th scope=row>2</th><td>GO:0043170</td><td>0.0002857317</td><td> 2.284441</td><td>28.018519</td><td>43</td><td>221</td><td>macromolecule metabolic process                 </td></tr>
	<tr><th scope=row>3</th><td>GO:0050794</td><td>0.0005117877</td><td> 2.273786</td><td>22.440171</td><td>36</td><td>177</td><td>regulation of cellular process                  </td></tr>
	<tr><th scope=row>4</th><td>GO:0006139</td><td>0.0008023183</td><td> 2.429820</td><td>14.579772</td><td>26</td><td>115</td><td>nucleobase-containing compound metabolic process</td></tr>
	<tr><th scope=row>5</th><td>GO:0051716</td><td>0.0008465448</td><td> 2.649127</td><td>10.776353</td><td>21</td><td> 85</td><td>cellular response to stimulus                   </td></tr>
	<tr><th scope=row>6</th><td>GO:0051093</td><td>0.0012089516</td><td>12.103175</td><td> 1.014245</td><td> 5</td><td>  8</td><td>negative regulation of developmental process    </td></tr>
	<tr><th scope=row>7</th><td>GO:0048522</td><td>0.0012106665</td><td> 2.843381</td><td> 8.113960</td><td>17</td><td> 64</td><td>positive regulation of cellular process         </td></tr>
	<tr><th scope=row>8</th><td>GO:0031323</td><td>0.0014026244</td><td> 2.516681</td><td>11.156695</td><td>21</td><td> 88</td><td>regulation of cellular metabolic process        </td></tr>
	<tr><th scope=row>9</th><td>GO:0006259</td><td>0.0017043505</td><td> 4.946502</td><td> 2.535613</td><td> 8</td><td> 20</td><td>DNA metabolic process                           </td></tr>
	<tr><th scope=row>10</th><td>GO:0007606</td><td>0.0018326653</td><td> 5.728997</td><td> 2.028490</td><td> 7</td><td> 16</td><td>sensory perception of chemical stimulus         </td></tr>
</tbody>
</table>




```R
summary(dnBp)[1:10,]
```


<table class="dataframe">
<caption>A data.frame: 10 × 7</caption>
<thead>
	<tr><th></th><th scope=col>GOBPID</th><th scope=col>Pvalue</th><th scope=col>OddsRatio</th><th scope=col>ExpCount</th><th scope=col>Count</th><th scope=col>Size</th><th scope=col>Term</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>GO:0045924</td><td>8.541596e-05</td><td>      Inf</td><td>  9.880342</td><td> 17</td><td> 17</td><td>regulation of female receptivity             </td></tr>
	<tr><th scope=row>2</th><td>GO:0060180</td><td>8.541596e-05</td><td>      Inf</td><td>  9.880342</td><td> 17</td><td> 17</td><td>female mating behavior                       </td></tr>
	<tr><th scope=row>3</th><td>GO:0046692</td><td>4.562821e-04</td><td>      Inf</td><td>  8.136752</td><td> 14</td><td> 14</td><td>sperm competition                            </td></tr>
	<tr><th scope=row>4</th><td>GO:0046008</td><td>7.958947e-04</td><td>      Inf</td><td>  7.555556</td><td> 13</td><td> 13</td><td>regulation of female receptivity, post-mating</td></tr>
	<tr><th scope=row>5</th><td>GO:0019953</td><td>9.709946e-04</td><td> 1.749826</td><td>111.589744</td><td>130</td><td>192</td><td>sexual reproduction                          </td></tr>
	<tr><th scope=row>6</th><td>GO:0044703</td><td>1.938649e-03</td><td>11.183206</td><td>  9.299145</td><td> 15</td><td> 16</td><td>multi-organism reproductive process          </td></tr>
	<tr><th scope=row>7</th><td>GO:0007320</td><td>3.186017e-03</td><td>10.411168</td><td>  8.717949</td><td> 14</td><td> 15</td><td>insemination                                 </td></tr>
	<tr><th scope=row>8</th><td>GO:0007620</td><td>3.186017e-03</td><td>10.411168</td><td>  8.717949</td><td> 14</td><td> 15</td><td>copulation                                   </td></tr>
	<tr><th scope=row>9</th><td>GO:0044706</td><td>3.186017e-03</td><td>10.411168</td><td>  8.717949</td><td> 14</td><td> 15</td><td>multi-multicellular organism process         </td></tr>
	<tr><th scope=row>10</th><td>GO:0006629</td><td>3.575362e-03</td><td> 2.786667</td><td> 24.410256</td><td> 33</td><td> 42</td><td>lipid metabolic process                      </td></tr>
</tbody>
</table>



### GO Cellular Component:


```R
upParams = new("GOHyperGParams",
               geneIds = selectGenesUp,
               universeGeneIds = universeGenes,
               annotation = "org.Dm.eg.db",
               ontology = "CC", # cellular components
               pvalueCutoff = 0.01,
               conditional = FALSE,
               testDirection = "over")

dnParams = new("GOHyperGParams",
               geneIds = selectGenesDn,
               universeGeneIds = universeGenes,
               annotation = "org.Dm.eg.db",
               ontology = "CC",
               pvalueCutoff = 0.01,
               conditional = FALSE,
               testDirection = "over")

upCC = hyperGTest(upParams)
dnCC = hyperGTest(dnParams)
```


```R
summary(upCC)[1:5,]

```


<table class="dataframe">
<caption>A data.frame: 5 × 7</caption>
<thead>
	<tr><th></th><th scope=col>GOCCID</th><th scope=col>Pvalue</th><th scope=col>OddsRatio</th><th scope=col>ExpCount</th><th scope=col>Count</th><th scope=col>Size</th><th scope=col>Term</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>GO:0032991</td><td>0.0001548667</td><td>2.703448</td><td>14.835905</td><td>28</td><td>123</td><td>protein-containing complex              </td></tr>
	<tr><th scope=row>2</th><td>GO:0005634</td><td>0.0021556880</td><td>2.172414</td><td>17.127630</td><td>28</td><td>142</td><td>nucleus                                 </td></tr>
	<tr><th scope=row>3</th><td>GO:0140535</td><td>0.0025066520</td><td>4.490842</td><td> 2.653576</td><td> 8</td><td> 22</td><td>intracellular protein-containing complex</td></tr>
	<tr><th scope=row>4</th><td>GO:1990234</td><td>0.0076404521</td><td>4.627500</td><td> 1.929874</td><td> 6</td><td> 16</td><td>transferase complex                     </td></tr>
	<tr><th scope=row>5</th><td>GO:0005886</td><td>0.0085012825</td><td>2.285714</td><td> 8.805049</td><td>16</td><td> 73</td><td>plasma membrane                         </td></tr>
</tbody>
</table>




```R
summary(dnCC)[1:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 7</caption>
<thead>
	<tr><th></th><th scope=col>GOCCID</th><th scope=col>Pvalue</th><th scope=col>OddsRatio</th><th scope=col>ExpCount</th><th scope=col>Count</th><th scope=col>Size</th><th scope=col>Term</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>GO:0005576</td><td>2.033933e-19</td><td>6.031503</td><td>111.27069</td><td>161</td><td>188</td><td>extracellular region</td></tr>
	<tr><th scope=row>2</th><td>GO:0005615</td><td>3.431644e-18</td><td>6.333551</td><td> 98.24965</td><td>144</td><td>166</td><td>extracellular space </td></tr>
	<tr><th scope=row>3</th><td>GO:0015629</td><td>4.083513e-03</td><td>6.065432</td><td> 11.24544</td><td> 17</td><td> 19</td><td>actin cytoskeleton  </td></tr>
	<tr><th scope=row>NA</th><td>NA        </td><td>          NA</td><td>      NA</td><td>       NA</td><td> NA</td><td> NA</td><td>NA                  </td></tr>
	<tr><th scope=row>NA.1</th><td>NA        </td><td>          NA</td><td>      NA</td><td>       NA</td><td> NA</td><td> NA</td><td>NA                  </td></tr>
</tbody>
</table>



### GO Molecular Functions


```R
upParams = new("GOHyperGParams",
               geneIds = selectGenesUp,
               universeGeneIds = universeGenes,
               annotation = "org.Dm.eg.db",
               ontology = "MF", # molecular functions
               pvalueCutoff = 0.01,
               conditional = FALSE,
               testDirection = "over")

dnParams = new("GOHyperGParams",
               geneIds = selectGenesDn,
               universeGeneIds = universeGenes,
               annotation = "org.Dm.eg.db",
               ontology = "MF",
               pvalueCutoff = 0.01,
               conditional = FALSE,
               testDirection = "over")

upMF = hyperGTest(upParams)
dnMF = hyperGTest(dnParams)
```


```R
summary(upMF)[1:5,]
summary(dnMF)[1:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 7</caption>
<thead>
	<tr><th></th><th scope=col>GOMFID</th><th scope=col>Pvalue</th><th scope=col>OddsRatio</th><th scope=col>ExpCount</th><th scope=col>Count</th><th scope=col>Size</th><th scope=col>Term</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>GO:0004888</td><td>0.0004995426</td><td> 7.846154</td><td> 1.6903409</td><td> 7</td><td> 14</td><td>transmembrane signaling receptor activity</td></tr>
	<tr><th scope=row>2</th><td>GO:0003676</td><td>0.0023963043</td><td> 2.380000</td><td>11.5909091</td><td>21</td><td> 96</td><td>nucleic acid binding                     </td></tr>
	<tr><th scope=row>3</th><td>GO:0022836</td><td>0.0053023253</td><td> 5.147679</td><td> 1.8110795</td><td> 6</td><td> 15</td><td>gated channel activity                   </td></tr>
	<tr><th scope=row>4</th><td>GO:0019904</td><td>0.0062243461</td><td>22.609756</td><td> 0.4829545</td><td> 3</td><td>  4</td><td>protein domain specific binding          </td></tr>
	<tr><th scope=row>5</th><td>GO:0005488</td><td>0.0086905151</td><td> 1.787013</td><td>39.2400568</td><td>50</td><td>325</td><td>binding                                  </td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 5 × 7</caption>
<thead>
	<tr><th></th><th scope=col>GOMFID</th><th scope=col>Pvalue</th><th scope=col>OddsRatio</th><th scope=col>ExpCount</th><th scope=col>Count</th><th scope=col>Size</th><th scope=col>Term</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>GO:0030246</td><td>7.786161e-06</td><td>11.07071</td><td>20.67045</td><td>32</td><td>34</td><td>carbohydrate binding            </td></tr>
	<tr><th scope=row>2</th><td>GO:0004866</td><td>1.414871e-05</td><td>     Inf</td><td>13.37500</td><td>22</td><td>22</td><td>endopeptidase inhibitor activity</td></tr>
	<tr><th scope=row>3</th><td>GO:0030414</td><td>1.414871e-05</td><td>     Inf</td><td>13.37500</td><td>22</td><td>22</td><td>peptidase inhibitor activity    </td></tr>
	<tr><th scope=row>4</th><td>GO:0048029</td><td>1.414871e-05</td><td>     Inf</td><td>13.37500</td><td>22</td><td>22</td><td>monosaccharide binding          </td></tr>
	<tr><th scope=row>5</th><td>GO:0061135</td><td>1.414871e-05</td><td>     Inf</td><td>13.37500</td><td>22</td><td>22</td><td>endopeptidase regulator activity</td></tr>
</tbody>
</table>




```R
summary(upMF)$Term
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'transmembrane signaling receptor activity'</li><li>'nucleic acid binding'</li><li>'gated channel activity'</li><li>'protein domain specific binding'</li><li>'binding'</li><li>'ligand-gated monoatomic ion channel activity'</li><li>'ligand-gated channel activity'</li></ol>




```R
summary(dnMF)$Term
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'carbohydrate binding'</li><li>'endopeptidase inhibitor activity'</li><li>'peptidase inhibitor activity'</li><li>'monosaccharide binding'</li><li>'endopeptidase regulator activity'</li><li>'serine-type endopeptidase inhibitor activity'</li><li>'peptidase regulator activity'</li><li>'procollagen-proline 4-dioxygenase activity'</li><li>'procollagen-proline dioxygenase activity'</li><li>'peptidyl-proline dioxygenase activity'</li><li>'peptidyl-proline 4-dioxygenase activity'</li><li>'oxidoreductase activity'</li><li>'vitamin binding'</li><li>'carboxylic acid binding'</li><li>'L-ascorbic acid binding'</li><li>'organic acid binding'</li><li>'oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen'</li><li>'enzyme inhibitor activity'</li><li>'molecular function inhibitor activity'</li><li>'enzyme regulator activity'</li><li>'iron ion binding'</li><li>'2-oxoglutarate-dependent dioxygenase activity'</li><li>'dioxygenase activity'</li><li>'molecular function regulator activity'</li><li>'carboxylic ester hydrolase activity'</li></ol>



### KEGG ANALYSIS:


```R
upKeggParams = new("KEGGHyperGParams",
                   geneIds=selectGenesUp,
                   universeGeneIds=universeGenes,
                   annotation="org.Dm.eg.db",
                   pvalueCutoff=1,
                   testDirection="over")

dnKeggParams = new("KEGGHyperGParams",
                   geneIds=selectGenesDn,
                   universeGeneIds=universeGenes,
                   annotation="org.Dm.eg.db",
                   pvalueCutoff=1,
                   testDirection="over")

upKEGG=hyperGTest(upKeggParams)
dnKEGG=hyperGTest(dnKeggParams)
```


```R
summary(upKEGG)[1:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 7</caption>
<thead>
	<tr><th></th><th scope=col>KEGGID</th><th scope=col>Pvalue</th><th scope=col>OddsRatio</th><th scope=col>ExpCount</th><th scope=col>Count</th><th scope=col>Size</th><th scope=col>Term</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>00563</td><td>0.09558824</td><td>     Inf</td><td>0.09558824</td><td>1</td><td>1</td><td>NA</td></tr>
	<tr><th scope=row>2</th><td>00670</td><td>0.09558824</td><td>     Inf</td><td>0.09558824</td><td>1</td><td>1</td><td>NA</td></tr>
	<tr><th scope=row>3</th><td>04140</td><td>0.09558824</td><td>     Inf</td><td>0.09558824</td><td>1</td><td>1</td><td>NA</td></tr>
	<tr><th scope=row>4</th><td>03450</td><td>0.09558824</td><td>     Inf</td><td>0.09558824</td><td>1</td><td>1</td><td>NA</td></tr>
	<tr><th scope=row>5</th><td>04080</td><td>0.18267974</td><td>10.16667</td><td>0.19117647</td><td>1</td><td>2</td><td>NA</td></tr>
</tbody>
</table>




```R
summary(dnKEGG)[1:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 7</caption>
<thead>
	<tr><th></th><th scope=col>KEGGID</th><th scope=col>Pvalue</th><th scope=col>OddsRatio</th><th scope=col>ExpCount</th><th scope=col>Count</th><th scope=col>Size</th><th scope=col>Term</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>01100</td><td>0.0001214777</td><td>3.929825</td><td>36.948529</td><td>48</td><td>67</td><td>NA</td></tr>
	<tr><th scope=row>2</th><td>00330</td><td>0.0005150782</td><td>     Inf</td><td> 6.617647</td><td>12</td><td>12</td><td>NA</td></tr>
	<tr><th scope=row>3</th><td>04142</td><td>0.0479516594</td><td>     Inf</td><td> 2.757353</td><td> 5</td><td> 5</td><td>NA</td></tr>
	<tr><th scope=row>4</th><td>00903</td><td>0.0585963330</td><td>4.022727</td><td> 6.066176</td><td> 9</td><td>11</td><td>NA</td></tr>
	<tr><th scope=row>5</th><td>00480</td><td>0.0891495640</td><td>     Inf</td><td> 2.205882</td><td> 4</td><td> 4</td><td>NA</td></tr>
</tbody>
</table>




```R

```
