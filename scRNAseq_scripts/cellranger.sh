######################## In vitro model ########################

### cellranger command
counts_dir="/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/"
cd $counts_dir

# R1 and R2 fastq files from each sample to be in separate folders under "fastq_dir"
fastq_dir="/dataOS/rcalandrelli/MARGI/RNAseq_02152019/fastq"

for i in $(find $fastq_dir -maxdepth 1 -mindepth 1 -type d -printf '%f\n'); do
	cellranger count --id="count_""$i" \
	--transcriptome=/dataOS/rcalandrelli/Software/refdata-cellranger-GRCh38-3.0.0 \
	--fastqs=/dataOS/rcalandrelli/MARGI/RNAseq_02152019/fastq/"$i" \
	--sample="$i" \
	--localcores=24
done

### Unzip matrices files
for i in $(find $fastq_dir -maxdepth 1 -mindepth 1 -type d -printf '%f\n'); do
	cd $counts_dir$i"/outs/filtered_feature_bc_matrix/"
	gunzip -k *.gz
done

######################## In vivo model ########################

### cellranger command
counts_dir="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/counts/"
cd $counts_dir

# R1 and R2 fastq files from each sample to be in separate folders under "fastq_dir"
fastq_dir="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/fastq"

for i in $(find $fastq_dir -maxdepth 1 -mindepth 1 -type d -printf '%f\n'); do
cellranger count --id="count_""$i" \
--transcriptome=/dataOS/rcalandrelli/Software/refdata-cellranger-GRCh38-3.0.0 \
--fastqs=/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/fastq/"$i" \
--sample="$i" \
--localcores=24
done

### Unzip matrices files
for i in $(find $fastq_dir -maxdepth 1 -mindepth 1 -type d -printf '%f\n'); do
	cd $counts_dir$i"/outs/filtered_feature_bc_matrix/"
	gunzip -k *.gz
done