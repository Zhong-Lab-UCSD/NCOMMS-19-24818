chrSizeFile="/mnt/extraids/OceanStor-SysCmn-1/rcalandrelli/Project____InterChrom/Data/4DN_files_reference/4DNFI823LSII.hg38.mainonly.chrom.sizes"

ChrList="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"


### Contact maps at 200 kb

# HUVEC replicate 1
InputFileBase="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm/01_04____Split_iMARGI_ChrBrChr/Single_HUVECcontrol2."
OutputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm/HUVECcontrol__Matrix__Resolution200000/"

InputFileBase="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered/01_04____Split_iMARGI_ChrBrChr/Single_HUVECcontrol2."
OutputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered/HUVECcontrol__Matrix__Resolution200000/"

InputFileBase="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered/01_04____Split_iMARGI_ChrBrChr/Single_HUVECcontrol2."
OutputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered/HUVECcontrol__Matrix__Resolution200000/"


# HUVEC replicate 2
InputFileBase="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_control_igm____R2filtered/01_04____Split_iMARGI_ChrBrChr/Single_HUVECcontrol2."
OutputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_control_igm____R2filtered/HUVECcontrol__Matrix__Resolution200000/"

InputFileBase="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered/01_04____Split_iMARGI_ChrBrChr/Single_HUVECcontrol2."
OutputPath="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered/HUVECcontrol__Matrix__Resolution200000/"


# To be executed per each sample
Resolution=200000

for chrA in ${ChrList}
do
    for chrB in ${ChrList}
    do
        python2.7  /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/RNAvsDNA_matrix_between_a_pair_of_chromosomes.PY27.py    -i ${InputFileBase}${chrA}_${chrB}.bedpe     -o ${OutputPath}    --chrSize ${chrSizeFile}    -c 200000    -r ${Resolution}    --rna_chr ${chrA}    --dna_chr ${chrB}
    done
done






