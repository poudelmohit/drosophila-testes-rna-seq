
## Download raw files

    cd raw_data/

    cat SRR_Acc_List.txt | xargs -n1 prefetch
    ls -lh */*.sra
    
    for file in */*.sra;do
        out=$(basename $file .sra)
        echo $out
        fastq-dump --split-3 $file -v -O $out
    done

    ls -lh */*.fastq
    

## Quality Assessment of raw reads:

    mkdir -p ../working_data/initial_fastqc/
    fastqc */*.fastq -threads 20 --outdir ../working_data/initial_fastqc/

    ### Based on the fastqc reports, no quality control is required !!
    ### There is no presence of adapter sequences as well.


## Download Reference Transcriptome sequence:
        
    mkdir -p ../working_data/assembly && cd $_

    wget ftp://ftp.ensembl.org/pub/release-87/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz -v
    
    wget ftp://ftp.ensembl.org/pub/release-87/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.87.gtf.gz -v

    gunzip *.gz

    mkdir index/
    bowtie2-build *.fa.gz index/fly --verbose

### Alignment of RNA-seq reads to the reference transcriptome:

    STAR --runThreadN 20 \
         --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles Dro*.fa \
         --sjdbGTFfile Dro*.gtf \
         --sjdbOverhang 149
         --genomeSAindexNbases 12
    
    mkdir star_mapped/

    for read1 in ../../raw_data/*/*_1.fastq;do
        base=$(basename $read1 "_1.fastq")
        read2=$(echo $read1 | sed 's/_1/_2/')
        echo $base
        echo $read1
        echo $read2
        STAR --runThreadN 20 \
            --genomeDir star_index \
            --readFilesIn $read1 $read2 \
            --outFileNamePrefix star_mapped/${base}_ \
            --outSAMtype BAM SortedByCoordinate\
            --outFilterIntronMotifs RemoveNoncanonical \
            --outReadsUnmapped Fastx
    done

    ls -lh star_mapped/*.bam

    samtools quickcheck star_mapped/*.out.bam
    samtools view -H star_mapped/SRR21091725_Aligned.sortedByCoord.out.bam
    
## Gene Expression Quantificaction

    mkdir -p ../../output/counts/

    ls star_mapped/*.bam
    
    # iterating over all bam outputs to genarate feature count table:

    for bam in star_mapped/*.bam; do
        base=$(basename $bam '_Aligned.sortedByCoord.out.bam')
        echo $base
        featureCounts -p -t exon -g gene_id --verbose -a Dro*.gtf -o ../../output/counts/${base}_counts.txt $bam
    done

### Merge multiple feature count files into single csv file
    
    head ../../output/counts/*.txt.summary
    head -n6 ../../output/counts/*.txt


    for file in ../../output/*_counts.txt; do
        base=$(basename $file _counts.txt)
        awk 'NR > 2 {print $1, $7}' $file > ../../output/counts/${base}_cleaned.txt
    done

    
    printf "id tr721 tr722 con723 con725\n" > ../../output/counts/merged_counts.csv
    
    paste ../../output/counts/*_cleaned.txt | awk '{print $1, $2, $4, $6, $8, $10, $12}' | sed 's/[[:space:]]*$//' >> ../../output/counts/merged_counts.csv
    ## Here, I merged all files based on the values in first column (geneid). This approach works only if all samples have exactly same geneid present across all rows. For eg: FBgn0267431 if present in 2nd row of sample 1, it should be true for sample 2,3,4,etc. (random sample name used for example)
    ## This was checked manually in excel, but can be done with pandas or dplyr as well.
    
    head -n4 ../../output/counts/merged_counts.csv
    head -n4 merged_counts.csv

## Downstream Analysis with jupyter notebook
    
    cd ../../scripts/
    jupyter-notebook









        