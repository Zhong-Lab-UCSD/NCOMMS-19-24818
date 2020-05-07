### HUVEC replicate 1
InputRawDataR1="/mnt/extraids/OceanStor-SysCmn-1/chenweizhong/Project____MARGI2/igm_storage1_ucsd/180613_K00180_0625_AHVCVTBBXX_PE100_Combo/HUVECcontrollibrary2_WW_S3_L005_R1_001.fastq.gz"
InputRawDataR2="/mnt/extraids/OceanStor-SysCmn-1/chenweizhong/Project____MARGI2/igm_storage1_ucsd/180613_K00180_0625_AHVCVTBBXX_PE100_Combo/HUVECcontrollibrary2_WW_S3_L005_R2_001.fastq.gz"
OutputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm/"

InputRawDataR1="/mnt/extraids/OceanStor-SysCmn-1/chenweizhong/Project____MARGI2/igm_storage1_ucsd/180918_K00180_0676_BHWCYNBBXX_PE100_Combo/HUVEC-H_T3d-P6_IGM_S1_L006_R1_001.fastq.gz"
InputRawDataR2="/mnt/extraids/OceanStor-SysCmn-1/chenweizhong/Project____MARGI2/igm_storage1_ucsd/180918_K00180_0676_BHWCYNBBXX_PE100_Combo/HUVEC-H_T3d-P6_IGM_S1_L006_R2_001.fastq.gz"
OutputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered"

InputRawDataR1="/mnt/extraids/OceanStor-SysCmn-1/chenweizhong/Project____MARGI2/igm_storage1_ucsd/180918_K00180_0676_BHWCYNBBXX_PE100_Combo/HUVEC-H_T7d-P7_IGM_S2_L007_R1_001.fastq.gz"
InputRawDataR2="/mnt/extraids/OceanStor-SysCmn-1/chenweizhong/Project____MARGI2/igm_storage1_ucsd/180918_K00180_0676_BHWCYNBBXX_PE100_Combo/HUVEC-H_T7d-P7_IGM_S2_L007_R2_001.fastq.gz"
OutputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered"


### HUVEC replicate 2
InputRawDataR1="/mnt/extraids/OceanStor-SysCmn-2/projects/margi/igm_data/190118_K00180_0743_BH2TV7BBXY_PE100_Combo/HUVEC_Control_S1_L006_R1_001.fastq.gz"
InputRawDataR2="/mnt/extraids/OceanStor-SysCmn-2/projects/margi/igm_data/190118_K00180_0743_BH2TV7BBXY_PE100_Combo/HUVEC_Control_S1_L006_R2_001.fastq.gz"
OutputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_control_igm____R2filtered/"

InputRawDataR1="/mnt/extraids/OceanStor-SysCmn-2/projects/margi/igm_data/190118_K00180_0743_BH2TV7BBXY_PE100_Combo/HUVEC_HT7d_S3_L008_R1_001.fastq.gz"
InputRawDataR2="/mnt/extraids/OceanStor-SysCmn-2/projects/margi/igm_data/190118_K00180_0743_BH2TV7BBXY_PE100_Combo/HUVEC_HT7d_S3_L008_R2_001.fastq.gz"
OutputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered/"

#### This code is to be executed for every sample above

CpuN=30

cd     ${OutputPath}
mkdir  ${OutputPath}1____raw_and_trimmed
mkdir  ${OutputPath}2____R1_qc
mkdir  ${OutputPath}2____R2_qc
mkdir  ${OutputPath}3____rmdup
mkdir  ${OutputPath}4____R1_mapping
mkdir  ${OutputPath}4____R2_mapping
mkdir  ${OutputPath}5____MARGI
mkdir  ${OutputPath}6____Biological



echo  R1    Trim first two bases
zcat "${InputRawDataR1}" | awk 'NR%2==0{print substr($0,3,length($0))} NR%2!=0{print $0}'    \
>>   ${OutputPath}1____raw_and_trimmed/R1_raw_afterTrim.fastq

echo  R1    Copy raw files
zcat "${InputRawDataR1}"    \
>>   ${OutputPath}1____raw_and_trimmed/R1_raw_withoutTrim.fastq

echo  R2    Copy raw files
zcat "${InputRawDataR2}"    \
>>   ${OutputPath}1____raw_and_trimmed/R2_raw.fastq






echo  R1    fastqc
fastqc   ${OutputPath}1____raw_and_trimmed/R1_raw_afterTrim.fastq    \
--outdir=${OutputPath}2____R1_qc    --extract

echo  R2    fastqc
fastqc   ${OutputPath}1____raw_and_trimmed/R2_raw.fastq    \
--outdir=${OutputPath}2____R2_qc    --extract






echo  Collapse duplicated pair-end reads
echo ${OutputPath}1____raw_and_trimmed/R1_raw_withoutTrim.fastq    \
>    ${OutputPath}3____rmdup/input_list.txt
echo ${OutputPath}1____raw_and_trimmed/R2_raw.fastq                \
>>   ${OutputPath}3____rmdup/input_list.txt

fastuniq    \
-i  ${OutputPath}3____rmdup/input_list.txt    \
-o  ${OutputPath}3____rmdup/R1_rmdup_withoutTrim.fa    \
-p  ${OutputPath}3____rmdup/R2_rmdup.fa                \
-t f    \
-c 0

echo  Removing first two bases from R1 reads
awk 'NR%2==0{print substr($0,3,length($0))} NR%2!=0{print $0}'    \
${OutputPath}3____rmdup/R1_rmdup_withoutTrim.fa >    \
${OutputPath}3____rmdup/R1_rmdup_afterTrim.fa

echo  Statistics R1 reads
python2.7    /mnt/extraids/OceanStor-SysCmn-1/rcalandrelli/Project____MARGI2/R1_RNA_fasta_reads_stats.py    -i ${OutputPath}3____rmdup/R1_rmdup_afterTrim.fa    -o ${OutputPath}3____rmdup/R1_rmdup_afterTrim.stats2    -n 2

echo  Filtering R2 reads
python2.7    /mnt/extraids/OceanStor-SysCmn-1/rcalandrelli/Project____MARGI2/R2_DNA_fasta_reads_filtering.py    -i ${OutputPath}3____rmdup/R2_rmdup.fa    -o ${OutputPath}3____rmdup/R2_rmdup_afterFilter.fa    -m CT    -n 2


gzip  -c  ${OutputPath}3____rmdup/R1_rmdup_withoutTrim.fa  > ${OutputPath}3____rmdup/R1_rmdup_withoutTrim.fa.gz
gzip  -c  ${OutputPath}3____rmdup/R1_rmdup_afterTrim.fa    > ${OutputPath}3____rmdup/R1_rmdup_afterTrim.fa.gz
gzip  -c  ${OutputPath}3____rmdup/R2_rmdup.fa              > ${OutputPath}3____rmdup/R2_rmdup.fa.gz
gzip  -c  ${OutputPath}3____rmdup/R2_rmdup_afterFilter.fa  > ${OutputPath}3____rmdup/R2_rmdup_afterFilter.fa.gz





echo  R1    mapping
STAR    \
--genomeDir    /home/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/STARindex_withSJ/    \
--readFilesIn          ${OutputPath}3____rmdup/R1_rmdup_afterTrim.fa.gz    \
--readFilesCommand zcat    \
--outSAMtype BAM SortedByCoordinate    \
--outSAMstrandField intronMotif    \
--outReadsUnmapped Fastx    \
--outSAMattributes All    \
--outFileNamePrefix    ${OutputPath}4____R1_mapping/R1    \
--quantMode GeneCounts    \
--sjdbGTFfile    /home/sysbio/Genomes/Homo_sapiens/Ensembl/GRCH38_hg38/Annotation/Genes/Homo_sapiens.GRCh38.84.chr.gtf    \
--alignEndsType Extend5pOfRead1    \
--seedSearchStartLmax 25    \
--runThreadN ${CpuN}    \
--outFilterScoreMin 21    \
--outFilterScoreMinOverLread 0.13    \
--outFilterMatchNmin 21    \
--outFilterMatchNminOverLread 0.13    \
--limitBAMsortRAM 46624240928
echo

echo    R2    mapping
STAR    \
--genomeDir    /home/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/STARindex_noSJ    \
--readFilesIn          ${OutputPath}3____rmdup/R2_rmdup_afterFilter.fa.gz    \
--readFilesCommand zcat    \
--outSAMtype BAM SortedByCoordinate    \
--outSAMstrandField intronMotif    \
--outReadsUnmapped Fastx    \
--outSAMattributes All    \
--outFileNamePrefix    ${OutputPath}4____R2_mapping/R2    \
--alignIntronMax 1    \
--alignEndsType Extend5pOfRead1    \
--seedSearchStartLmax 25    \
--runThreadN ${CpuN}    \
--outFilterScoreMin 21    \
--outFilterScoreMinOverLread 0.13    \
--outFilterMatchNmin 21    \
--outFilterMatchNminOverLread 0.13    \
--limitBAMsortRAM 35720943283
echo







echo    R1    sorting......
samtools sort -n -@ ${CpuN}    --output-fmt BAM    \
   ${OutputPath}4____R1_mapping/R1Aligned.sortedByCoord.out.bam    \
-o ${OutputPath}4____R1_mapping/R1Aligned.sortedByName.out.bam

samtools view    -h    --output-fmt SAM    \
-o ${OutputPath}4____R1_mapping/R1Aligned.sortedByCoord.out.sam    \
   ${OutputPath}4____R1_mapping/R1Aligned.sortedByCoord.out.bam

samtools view    -h    --output-fmt SAM    \
-o ${OutputPath}4____R1_mapping/R1Aligned.sortedByName.out.sam    \
   ${OutputPath}4____R1_mapping/R1Aligned.sortedByName.out.bam
echo

echo    R2    sorting......
samtools sort -n -@ ${CpuN}    --output-fmt BAM    \
   ${OutputPath}4____R2_mapping/R2Aligned.sortedByCoord.out.bam    \
-o ${OutputPath}4____R2_mapping/R2Aligned.sortedByName.out.bam

samtools view    -h    --output-fmt SAM    \
-o ${OutputPath}4____R2_mapping/R2Aligned.sortedByCoord.out.sam    \
   ${OutputPath}4____R2_mapping/R2Aligned.sortedByCoord.out.bam

samtools view    -h    --output-fmt SAM    \
-o ${OutputPath}4____R2_mapping/R2Aligned.sortedByName.out.sam    \
   ${OutputPath}4____R2_mapping/R2Aligned.sortedByName.out.bam
echo







echo Merging mates
samtools merge -nr -@ ${CpuN}    \
${OutputPath}5____MARGI/tmp.bam    \
${OutputPath}4____R1_mapping/R1Aligned.sortedByCoord.out.bam    \
${OutputPath}4____R2_mapping/R2Aligned.sortedByCoord.out.bam

echo Sort by name
samtools sort -n -@ ${CpuN}     --output-fmt BAM    \
   ${OutputPath}5____MARGI/tmp.bam    \
-o ${OutputPath}5____MARGI/R1R2.sortedByName.bam

echo Fixing mates and computing distances between mates
python2.7  /mnt/extraids/ExtSpace1/rcalandrelli/tools/bin/fixBAMmates    \
-b ${OutputPath}5____MARGI/R1R2.sortedByName.bam    \
-g R1Aligned.sortedByCoord.out,R2Aligned.sortedByCoord.out    \
-o ${OutputPath}5____MARGI/Fixed_mates.sortedByName.bam

echo Sorting by coordinate
samtools sort -@ ${CpuN}     --output-fmt BAM    \
   ${OutputPath}5____MARGI/Fixed_mates.sortedByName.bam    \
-o ${OutputPath}5____MARGI/Fixed_mates.sortedByCoord.bam

echo Indexing bam file
samtools index    ${OutputPath}5____MARGI/Fixed_mates.sortedByCoord.bam

echo Counting read types
chmod +w    ${OutputPath}5____MARGI/readTypes.txt
python2.7  /mnt/extraids/ExtSpace1/rcalandrelli/tools/bin/countReadTypes   \
-a ${OutputPath}5____MARGI/Fixed_mates.sortedByCoord.bam    \
-l 2000    \
-o ${OutputPath}5____MARGI/readTypes.txt
totalReads=$(samtools flagstat ${OutputPath}5____MARGI/Fixed_mates.sortedByCoord.bam | awk 'NR==5{print $1}')
echo -e totalReads"\t"$totalReads >> ${OutputPath}5____MARGI/readTypes.txt
echo


echo annotateBAM

python2.7  /mnt/extraids/OceanStor-SysCmn-1/rcalandrelli/Project____MARGI2/_MARGI_Programs_/annotateBAM_update   \
-b    ${OutputPath}5____MARGI/Fixed_mates.sortedByName.bam    \
-r    /home/sysbio/Genomes/Homo_sapiens/Ensembl/GRCH38_hg38/Annotation/Genes/Homo_sapiens.GRCh38.84.chr.gtf    \
-o    ${OutputPath}5____MARGI/annot_exonsIntrons.txt


# This is the final BEDPE file used for the analysis

echo "Convert from txt to BEDPE"

awk -v OFS='\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1, $2-=1, $3, $5, $6-=1, $7, $11, "1", $4, $8, $9, ".", $10, ".", abs($2-$7)-1}'  annot_exonsIntrons.txt > temp.bedpe

awk 'BEGIN {OFS=FS="\t"} $1!=$4 {$15="Inter-chromosome"}1' /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered/5____MARGI/temp.bedpe > /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered/5____MARGI/annot_exonsIntrons.bedpe



