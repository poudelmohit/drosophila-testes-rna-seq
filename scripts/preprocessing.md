
# Download raw files

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

# Based on the fastqc reports, no quality control is required !!
# There is no presence of adapter sequences as well.

## Adapter Trimming (SKIPPING FOR NOW AS FASTQC REPORT IS OKAY !!):

    trimmomatic PE -phred33 \
    input_R1.fastq.gz input_R2.fastq.gz \
    output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
    output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

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
    samtools view -H star_mapped/SRR21091720_Aligned.sortedByCoord.out.bam
    
## Gene Expression Quantificaction


    mkdir -p ../../output/

    ls star_mapped/*.bam
    
    for bam in star_mapped/*.bam; do
        base=$(basename $bam '_Aligned.sortedByCoord.out.bam')
        echo $base
        featureCounts -p -t exon -g gene_id --verbose -a Dro*.gtf -o ../../output/${base}_counts.txt $bam
    done

    # worked for 5 samples, failed for 1:
    
     featureCounts -p -t exon -g gene_id \
    --verbose -a Dro*.gtf -o ../../output/read_counts_24.txt \
    star_mapped/SRR21091724_Aligned.sortedByCoord.out.bam

    head ../../output/*.txt.summary
    head -n6 ../../output/*.txt


    for file in ../../output/*_counts.txt; do
        base=$(basename $file _counts.txt)
        awk 'NR > 2 {print $1, $7}' $file > ../../output/${base}_cleaned.txt
    done


    printf "#sample_id con720 con721 con722 tr723 tr725 \n" > ../../output/merged_counts.csv

    paste ../../output/*_cleaned.txt | awk '{print $1, $2, $4, $6, $8, $10, $12}' >> ../../output/merged_counts.csv
    head -n4 ../../output/merged_counts.csv







        