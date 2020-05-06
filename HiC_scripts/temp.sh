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
-o /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control/ \
-1 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control/04-10-2019-Tri-HiC-HUVEC-control_S3_L003_R1_001.fastq \
-2 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control/04-10-2019-Tri-HiC-HUVEC-control_S3_L003_R2_001.fastq \
-e [MboI,Hinfl] \
-g /home/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome \
-p 32 \
-m 50000000

cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control
bedtools bamtobed -i HiCfile_pair1.bam | sort -k 4,4 | awk -v OFS='\t' '{print $1, $2+=1, $3, $6}' > HiCfile_pair1.txt
bedtools bamtobed -i HiCfile_pair2.bam | sort -k 4,4 | awk -v OFS='\t' '{print $1, $2+=1, $3, $6}' > HiCfile_pair2.txt
paste HiCfile_pair1.txt HiCfile_pair2.txt > HiCfile_paired.txt

## Day 3
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_run_preprocessing.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-o /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3/ \
-1 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3/04-10-2019-Tri-HiC-HUVEC-d3_S4_L004_R1_001.fastq \
-2 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3/04-10-2019-Tri-HiC-HUVEC-d3_S4_L004_R2_001.fastq \
-e [MboI,Hinfl] \
-g /home/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome \
-p 32 \
-m 50000000

cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3
bedtools bamtobed -i HiCfile_pair1.bam | sort -k 4,4 | awk -v OFS='\t' '{print $1, $2+=1, $3, $6}' > HiCfile_pair1.txt
bedtools bamtobed -i HiCfile_pair2.bam | sort -k 4,4 | awk -v OFS='\t' '{print $1, $2+=1, $3, $6}' > HiCfile_pair2.txt
paste HiCfile_pair1.txt HiCfile_pair2.txt > HiCfile_paired.txt

## Day 7
/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_run_preprocessing.sh \
-h /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/ \
-o /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7/ \
-1 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7/HiC_library_for_Seq_in_EB_HUVEC_d7_S3_L003_R1_001.fastq \
-2 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7/HiC_library_for_Seq_in_EB_HUVEC_d7_S3_L003_R2_001.fastq \
-e [MboI,Hinfl] \
-g /home/sysbio/Genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome \
-p 32 \
-m 50000000

cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7
bedtools bamtobed -i HiCfile_pair1.bam | sort -k 4,4 | awk -v OFS='\t' '{print $1, $2+=1, $3, $6}' > HiCfile_pair1.txt
bedtools bamtobed -i HiCfile_pair2.bam | sort -k 4,4 | awk -v OFS='\t' '{print $1, $2+=1, $3, $6}' > HiCfile_pair2.txt
paste HiCfile_pair1.txt HiCfile_pair2.txt > HiCfile_paired.txt

### Data normalization

### Control
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control
python /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
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
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3
python /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
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
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7
python /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
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

### Plotting chr6 matrix (to be executed once per each sample)
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7

python /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_analysis.py \
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

python /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_global_map_analysis.py \
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
--cutoff 20 \
--max_color "#460000"

##### TAD analysis
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7

# To be executed per each sample above at a time
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


#### Plot PCC matrix
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_yaffe_tanay.py \
--action plot_correlation \
-i /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/GSM1551550/HiCtool_chr6_1mb_correlation_matrix.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6


### PCC analysis
### Control
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control
python /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
-f /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/fend_files/arima_hg38/arima_hg38_gc_map.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e [MboI,Hinfl] \
-m Yaffe-Tanay

import hifive
# Filtering HiC fends
hic = hifive.HiC('HiC_project_object.hdf5')
hic.filter_fends(mininteractions=1, mindistance=0, maxdistance=0)
# Finding HiC distance function
hic.find_distance_parameters(numbins=90, minsize=200, maxsize=0)
hic.save('HiC_project_object_with_distance_parameters.hdf5')
# Learning correction parameters using the binning algorithm
hic.find_binning_fend_corrections(max_iterations=1000,
                                              mindistance=500000,
                                              maxdistance=0,
                                              num_bins=[20,20],
                                              model=['len','distance'],
                                              parameters=['even','even'],
                                              usereads='cis',
                                              learning_threshold=1.0)
hic.save('HiC_norm_binning.hdf5')

python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_yaffe_tanay.py \
--action normalize_enrich \
-i HiC_norm_binning.hdf5 \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--save_obs 1 \
--save_expect 1

### Day 3
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3
python /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
-f /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/fend_files/arima_hg38/arima_hg38_gc_map.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e [MboI,Hinfl] \
-m Yaffe-Tanay


### Day 7
cd /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7
python /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_hifive.py \
-f /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/fend_files/arima_hg38/arima_hg38_gc_map.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e [MboI,Hinfl] \
-m Yaffe-Tanay




######## Compartment analysis
python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_compartment_analysis.py \
--action calculate_pc \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1

python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_compartment_analysis.py \
--action plot_pc \
-i HiCtool_chr6_1mb_PC1.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1