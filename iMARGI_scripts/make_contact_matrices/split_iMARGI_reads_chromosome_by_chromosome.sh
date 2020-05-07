### To split iMARGI reads chromosome by chromosome

ChrSize="/mnt/extraids/OceanStor-SysCmn-1/rcalandrelli/Project____InterChrom/Data/4DN_files_reference/4DNFI823LSII.hg38.mainonly.chrom.sizes"

### HUVEC replicate 1
inputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm/5____MARGI/"
outputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm/01_04____Split_iMARGI_ChrBrChr/"


inputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered/5____MARGI/"
outputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered/01_04____Split_iMARGI_ChrBrChr/"

inputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered/5____MARGI/"
outputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered/01_04____Split_iMARGI_ChrBrChr/"


### HUVEC replicate 2
inputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_control_igm____R2filtered/5____MARGI/"
outputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_control_igm____R2filtered/01_04____Split_iMARGI_ChrBrChr/"

inputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered/5____MARGI/"
outputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered/01_04____Split_iMARGI_ChrBrChr/"


# To be executed for every sample
python3.5  /mnt/extraids/OceanStor-SysCmn-1/rcalandrelli/Project____InterChrom/ProgramsV3/bin/ChrByChr_split_of_iMARGI_reads.PY3.py    -i $inputPath"annot_exonsIntrons.bedpe"    -o $outputPath    --name Single_HUVECcontrol2    --chrSize  ${ChrSize}









