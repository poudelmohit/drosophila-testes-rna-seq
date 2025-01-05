
    mkdir -p raw_data working_data resources scripts output/plots metadata
    
    git init
    
    conda install -c bioconda sra-tools -y
    conda install bioconda::fastqc -y
    conda install bioconda::samtools -y
    conda install -c bioconda star -y
    conda install -c bioconda subread -y

    conda install conda-forge::r-base=4.4.2 

    conda install -c conda-forge r-essentials -y
    conda install -c conda-forge r-irkernel -y

    conda install -c bioconda bioconductor-deseq2 -y
    conda install -c bioconda seaborn -y
    conda install pandas -y
    conda install matplotlib -y








