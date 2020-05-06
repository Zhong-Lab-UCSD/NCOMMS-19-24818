# Create a folder per each file
for i in {32877..32884}; do
    mkdir "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/LNA_rna_seq/Sample_"$i
done

PROJECT_PATH=/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/LNA_rna_seq
printf '%s\n' "Sample" "Input reads" "Average input read length" "Uniquely mapped reads" "% of uniquely mapped reads" "Number of reads mapped to multiple loci" "% of reads mapped to multiple loci" "Number of reads mapped to too many loci" "% of reads mapped to too many loci" "% of reads mapped to too many loci" "% of reads unmapped: too short" "% of reads unmapped: other" | paste -sd ',' >$PROJECT_PATH'/info_mapping.csv'

for i in {32877..32884}; do

    # Alignment job
    mkdir $PROJECT_PATH/Sample_$i/mapping
    STAR \
        --genomeDir /home/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/STARindex_withSJ \
        --readFilesIn $PROJECT_PATH/fastq/$i*.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --outSAMattributes All \
        --outFileNamePrefix $PROJECT_PATH/Sample_$i/mapping/ \
        --sjdbGTFfile /home/sysbio/Genomes/Homo_sapiens/Ensembl/GRCH38_hg38/Annotation/Genes/Homo_sapiens.GRCh38.84.chr.gtf \
        --runThreadN 32

    # Indexing the bam file
    samtools index $PROJECT_PATH/Sample_$i/mapping/Aligned.sortedByCoord.out.bam

    # Parse log file and write on csv file
    n1="Sample_"$i

    n2=$(grep 'Number of input reads |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    n3=$(grep 'Average input read length |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    n4=$(grep 'Uniquely mapped reads number |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    n5=$(grep 'Uniquely mapped reads % |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    n6=$(grep 'Number of reads mapped to multiple loci |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    n7=$(grep '% of reads mapped to multiple loci |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    n8=$(grep 'Number of reads mapped to too many loci |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    n9=$(grep '% of reads mapped to too many loci |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    n10=$(grep '% of reads mapped to too many loci |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    n11=$(grep '% of reads unmapped: too short |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    n12=$(grep '% of reads unmapped: other |' $PROJECT_PATH/Sample_$i/mapping/Log.final.out | sed 's/^.*|\t//')

    printf '%s\n' $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 $n9 $n10 $n11 $n12 | paste -sd ',' >>$PROJECT_PATH'/info_mapping.csv'

done

for i in {32877..32884}; do
    mkdir $PROJECT_PATH/Sample_$i/featureCounts

    featureCounts -T 32 \
        -t exon \
        -g gene_id \
        -s 2 \
        -a /home/sysbio/Genomes/Homo_sapiens/Ensembl/GRCH38_hg38/Annotation/Genes/Homo_sapiens.GRCh38.84.chr.gtf \
        -o $PROJECT_PATH/Sample_$i/featureCounts/counts.txt \
        $PROJECT_PATH/Sample_$i/mapping/Aligned.sortedByCoord.out.bam

done
