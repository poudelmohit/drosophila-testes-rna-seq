## Introduction

Reproducing and validating research findings is a cornerstone of scientific progress. In this project, I focused on replicating the transcriptomic analysis presented in the paper *A transcriptomic (RNA-seq) analysis of Drosophila melanogaster adult testes overexpressing microRNA-2b-1* by Manickam et al. This study aimed to understand the role of the miR-2b-1 microRNA in testicular bulging and stem cell differentiation in *Drosophila melanogaster*. Leveraging RNA sequencing (RNA-seq) data, the authors identified key genes and pathways involved in this biological process, providing valuable insights into microRNA-mediated regulation in stem cells.

---

## Methods

The original study used Illumina HiSeq to generate RNA-seq data and analyzed differential gene expression with tools like RNA STAR and edgeR. Gene ontology (GO) enrichment was performed to classify biological processes, cellular components, and molecular functions.

For my reproduction, I followed a similar workflow but implemented some differences:

1. **Preprocessing**: Downloaded raw data from GEO and processed it using `fastq-dump` and `fastqc` to assess quality. Reference genome files were retrieved from Ensembl, and alignment was performed with STAR, producing sorted BAM files.
2. **Differential Expression Analysis**: I used DESeq2 for differential expression analysis instead of edgeR. Feature counts from BAM files were merged into a single CSV file for downstream analysis.
3. **Visualization**: I created heatmaps, PCA plots, MA plots, and volcano plots to visualize the results. These provide clear insights into sample clustering, gene expression distributions, and significantly differentially expressed genes (SDEGs).

---

## Results

The analysis highlighted several key findings:

1. **Heatmap**: The hierarchical clustering heatmap illustrates the separation between control and treatment groups based on top variable genes. Genes are clearly upregulated or downregulated in response to miR-2b-1 overexpression.
2. **PCA Plot**: Principal Component Analysis (PCA) showed distinct clustering of the treatment and control groups, reaffirming clear differences in gene expression profiles.
3. **MA Plot**: The MA plot visualized log-fold changes and revealed hundreds of SDEGs under an FDR of 0.1.
4. **Volcano Plot**: Highlighted significantly upregulated and downregulated genes, with stringent cutoffs for log2 fold change and adjusted p-values. This plot emphasizes the miR-2b-1 impact on gene expression.
5. **GO Enrichment**: Biological processes such as DNA repair, regulation of immune response, and cellular differentiation were significantly enriched. Cellular components and molecular functions aligned with previously reported findings on microRNA involvement in stem cell regulation.

---

## Discussion

This reproduction confirms key findings of the original paper, including the impact of miR-2b-1 overexpression on differential gene expression and enriched biological pathways. However, some differences in statistical analysis and toolsets may influence subtle variations in results. Moving forward:

- **Improved Data Integration**: Incorporating additional datasets could provide broader context.
- **Comparative Analysis**: Comparing edgeR and DESeq2 results may enhance understanding of analysis discrepancies.
- **Functional Validation**: Experimental validation of identified genes would further corroborate computational findings.

By bridging computational reproducibility with biological insights, this project reinforces the role of miRNAs in developmental biology and provides a robust framework for future studies.
