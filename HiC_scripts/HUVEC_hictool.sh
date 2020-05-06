### Downsampling
less /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control/04-10-2019-Tri-HiC-HUVEC-control_S3_L003_R1_001.fastq | head -760000000 > /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled/04-10-2019-Tri-HiC-HUVEC-control_S3_L003_R1_001.fastq
less /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control/04-10-2019-Tri-HiC-HUVEC-control_S3_L003_R2_001.fastq | head -760000000 > /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled/04-10-2019-Tri-HiC-HUVEC-control_S3_L003_R2_001.fastq

less /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3/04-10-2019-Tri-HiC-HUVEC-d3_S4_L004_R1_001.fastq | head -760000000 > /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled/04-10-2019-Tri-HiC-HUVEC-d3_S4_L004_R1_001.fastq
less /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3/04-10-2019-Tri-HiC-HUVEC-d3_S4_L004_R2_001.fastq | head -760000000 > /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled/04-10-2019-Tri-HiC-HUVEC-d3_S4_L004_R2_001.fastq

less /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7/HiC_library_for_Seq_in_EB_HUVEC_d7_S3_L003_R1_001.fastq | head -760000000 > /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled/HiC_library_for_Seq_in_EB_HUVEC_d7_S3_L003_R1_001.fastq
less /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7/HiC_library_for_Seq_in_EB_HUVEC_d7_S3_L003_R2_001.fastq | head -760000000 > /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled/HiC_library_for_Seq_in_EB_HUVEC_d7_S3_L003_R2_001.fastq

## Pre-processing
chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_run_preprocessing.sh

## Day 0
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_run_preprocessing.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-o /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled/ \
-1 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled/04-10-2019-Tri-HiC-HUVEC-control_S3_L003_R1_001.fastq \
-2 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled/04-10-2019-Tri-HiC-HUVEC-control_S3_L003_R2_001.fastq \
-e [MboI,Hinfl] \
-q 30 \
-g /home/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome \
-p 32 \
-c 50000000

cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled
bedtools bamtobed -i HiCfile_pair1.bam | sort -k 4,4 > HiCfile_pair1.bed
bedtools bamtobed -i HiCfile_pair2.bam | sort -k 4,4 > HiCfile_pair2.bed
paste HiCfile_pair1.bed HiCfile_pair2.bed | awk -v OFS='\t' '{print $1, $2, $3, $7, $8, $9, $4, ".", $6, $12}' > HiCfile_paired.bedpe

## Day 3
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_run_preprocessing.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-o /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled/ \
-1 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled/04-10-2019-Tri-HiC-HUVEC-d3_S4_L004_R1_001.fastq \
-2 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled/04-10-2019-Tri-HiC-HUVEC-d3_S4_L004_R2_001.fastq \
-e [MboI,Hinfl] \
-q 30 \
-g /home/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome \
-p 32 \
-c 50000000

cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled
bedtools bamtobed -i HiCfile_pair1.bam | sort -k 4,4 > HiCfile_pair1.bed
bedtools bamtobed -i HiCfile_pair2.bam | sort -k 4,4 > HiCfile_pair2.bed
paste HiCfile_pair1.bed HiCfile_pair2.bed | awk -v OFS='\t' '{print $1, $2, $3, $7, $8, $9, $4, ".", $6, $12}' > HiCfile_paired.bedpe

## Day 7
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_run_preprocessing.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-o /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled/ \
-1 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled/HiC_library_for_Seq_in_EB_HUVEC_d7_S3_L003_R1_001.fastq \
-2 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled/HiC_library_for_Seq_in_EB_HUVEC_d7_S3_L003_R2_001.fastq \
-e [MboI,Hinfl] \
-q 30 \
-g /home/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome \
-p 32 \
-c 50000000

cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled
bedtools bamtobed -i HiCfile_pair1.bam | sort -k 4,4 > HiCfile_pair1.bed
bedtools bamtobed -i HiCfile_pair2.bam | sort -k 4,4 > HiCfile_pair2.bed
paste HiCfile_pair1.bed HiCfile_pair2.bed | awk -v OFS='\t' '{print $1, $2, $3, $7, $8, $9, $4, ".", $6, $12}' > HiCfile_paired.bedpe

### Data normalization

### Control
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
-f /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/fend_files/arima_hg38/arima_hg38_gc_map.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e [MboI,Hinfl] \
-m Hi-Corrector

# 1 mb resolution
chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh \
-i HiC_project_object.hdf5 \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-b 1000000 \
-s hg38 \
-p 24

chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_run_ic_mes.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-i ./observed_1000000/HiCtool_observed_global_1000000.txt \
-b 1000000 \
-q 100 \
-m 32000 \
-s hg38 \
-u 100

# 200 kb resolution
chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh \
-i HiC_project_object.hdf5 \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-b 200000 \
-s hg38 \
-p 24

chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-i ./observed_200000/HiCtool_observed_global_200000.txt \
-b 200000 \
-q 100 \
-m 32000 \
-s hg38 \
-u 100

# 40 kb resolution
chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh \
-i HiC_project_object.hdf5 \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-b 40000 \
-s hg38 \
-p 24

chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-i ./observed_40000/HiCtool_observed_global_40000.txt \
-b 40000 \
-q 100 \
-m 32000 \
-s hg38 \
-u 100

### Day 3
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
-f /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/fend_files/arima_hg38/arima_hg38_gc_map.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e [MboI,Hinfl] \
-m Hi-Corrector

# 1 mb resolution
chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh \
-i HiC_project_object.hdf5 \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-b 1000000 \
-s hg38 \
-p 24

chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-i ./observed_1000000/HiCtool_observed_global_1000000.txt \
-b 1000000 \
-q 100 \
-m 32000 \
-s hg38 \
-u 100

# 200 kb resolution
chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh \
-i HiC_project_object.hdf5 \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-b 200000 \
-s hg38 \
-p 24

chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-i ./observed_200000/HiCtool_observed_global_200000.txt \
-b 200000 \
-q 100 \
-m 32000 \
-s hg38 \
-u 100

# 40 kb resolution
chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh \
-i HiC_project_object.hdf5 \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-b 40000 \
-s hg38 \
-p 24

chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-i ./observed_40000/HiCtool_observed_global_40000.txt \
-b 40000 \
-q 100 \
-m 32000 \
-s hg38 \
-u 100

# Day 7
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
-f /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/fend_files/arima_hg38/arima_hg38_gc_map.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e [MboI,Hinfl] \
-m Hi-Corrector

# 1 mb
chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh \
-i HiC_project_object.hdf5 \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-b 1000000 \
-s hg38 \
-p 24

chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-i ./observed_1000000/HiCtool_observed_global_1000000.txt \
-b 1000000 \
-q 100 \
-m 32000 \
-s hg38 \
-u 100

# 200 kb resolution
chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh \
-i HiC_project_object.hdf5 \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-b 200000 \
-s hg38 \
-p 24

chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-i ./observed_200000/HiCtool_observed_global_200000.txt \
-b 200000 \
-q 100 \
-m 32000 \
-s hg38 \
-u 100

# 40 kb resolution
chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_observed.sh \
-i HiC_project_object.hdf5 \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-b 40000 \
-s hg38 \
-p 24

chmod u+x /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_normalize_global_matrix.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-i ./observed_40000/HiCtool_observed_global_40000.txt \
-b 40000 \
-q 100 \
-m 32000 \
-s hg38 \
-u 100

#### Plotting chr6 matrix (Supplementary Fig. 3a)

# To be executed once per each sample
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled

python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_analysis.py \
--action plot_map \
-i ./normalized_1000000/HiCtool_normalized_global_1000000.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--isGlobal 1 \
--tab_sep 1 \
--data_type normalized \
--chr_row 6 \
--chr_col 6 \
--my_colormap [white,red] \
--cutoff_type contact \
--cutoff 1 \
--max_color "#460000"

python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_analysis.py \
--action plot_map \
-i ./normalized_200000/chr6_chr6_200000.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 200000 \
-s hg38 \
--isGlobal 0 \
--tab_sep 1 \
--data_type normalized \
--chr_row 6 \
--chr_col 6 \
--chr_row_coord [80000000,120000000] \
--chr_col_coord [80000000,120000000] \
--my_colormap [white,red] \
--cutoff_type contact \
--cutoff 15 \
--max_color "#460000"

##### TAD analysis
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled

for i in "${chromosomes[@]}"; do
    python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_TAD_analysis.py \
	--action full_tad_analysis \
	-i "./normalized_40000/chr"$i"_chr"$i"_40000.txt" \
	-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
	-s hg38 \
	--isGlobal 0 \
	--tab_sep 1 \
	--chr $i \
	--data_type normalized
done


####### Yaffe-Tanay model to obtain observed/expected contact matrix and correlation matrix

# Day 0
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
-f /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/fend_files/arima_hg38/arima_hg38_gc_map.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e [MboI,Hinfl] \
-m Yaffe-Tanay

# Day 3
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
-f /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/fend_files/arima_hg38/arima_hg38_gc_map.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e [MboI,Hinfl] \
-m Yaffe-Tanay

# Day 7
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
-f /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/fend_files/arima_hg38/arima_hg38_gc_map.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e [MboI,Hinfl] \
-m Yaffe-Tanay

### Calculate and plot (Supplementary Fig. 3a) correlation matrices for chromosome 6

# Day 0
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_yaffe_tanay.py \
--action normalize_enrich \
-i HiC_norm_binning.hdf5 \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--save_obs 1 \
--save_expect 1

python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_correlation \
-i ./yaffe_tanay_1000000/chr6_1000000_correlation_matrix.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6

# Day 3
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_yaffe_tanay.py \
--action normalize_enrich \
-i HiC_norm_binning.hdf5 \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--save_obs 1 \
--save_expect 1

python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_correlation \
-i ./yaffe_tanay_1000000/chr6_1000000_correlation_matrix.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6

# Day 7
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_yaffe_tanay.py \
--action normalize_enrich \
-i HiC_norm_binning.hdf5 \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--save_obs 1 \
--save_expect 1

python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_correlation \
-i ./yaffe_tanay_1000000/chr6_1000000_correlation_matrix.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6

### A/B compartment analysis

# Day 0
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_compartment_analysis.py \
--action calculate_pc \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1

python2.7 /HiCtool-master/scripts/HiCtool_compartment_analysis.py \
--action plot_pc \
-i ./yaffe_tanay_1000000/chr6_1000000_PC1.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1 \
--plot_grid 0 \
--plot_axis 0

# Day 3
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_compartment_analysis.py \
--action calculate_pc \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1

python2.7 /HiCtool-master/scripts/HiCtool_compartment_analysis.py \
--action plot_pc \
-i ./yaffe_tanay_1000000/chr6_1000000_PC1.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1 \
--plot_grid 0 \
--plot_axis 0

# Day 7
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_compartment_analysis.py \
--action calculate_pc \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1

python2.7 /HiCtool-master/scripts/HiCtool_compartment_analysis.py \
--action plot_pc \
-i ./yaffe_tanay_1000000/chr6_1000000_PC1.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1 \
--plot_grid 0 \
--plot_axis 0