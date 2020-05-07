### Figure 2d
select_map="[75000000,125000000,75000000,125000000]"

input_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/Project____InterChrom/ResultsV3/05_01____InterChrom_matrix/HUVECcontrol__Matrix__Resolution200000/chr1_chr1.matrix"
se_chr1_200kb="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/index_overlap_SE_chr1_200kb.csv"
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm/HUVEC_control_chr1_75000000_125000000.pdf"
python2.7 plot_iMARGI.py -i $input_file -o $output_file --isGlobal 0 --chr_row 1 --chr_col 1 --bin_size 200000 --data_type Contact --select_part_map $select_map --plot_SE $se_chr1_200kb --cutoff_type contact --cutoff 15 --my_colormap "[#ffffff,#ff9999]" --max_color "#ff9999"


### Figure 3b
select_map="[195000000,235000000,80000000,120000000]"

# Day 0
input_file_1="/mnt/extraids/OceanStor-SysCmn-2/chenweizhong/Project____InterChrom/ResultsV3/05_01____InterChrom_matrix/HUVECcontrol__Matrix__Resolution200000/chr2_chr7.matrix"
input_file_2="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_control_igm____R2filtered/HUVECcontrol__Matrix__Resolution200000/chr2_chr7.matrix"
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/endoMT_analysis_NEW/chr2_chr7_200kb_control_MERGED.pdf"
python2.7 plot_iMARGI.py -i "["$input_file_1","$input_file_2"]" -o $output_file --isGlobal 0 --chr_row 2 --chr_col 7 --bin_size 200000 --data_type Contact --cutoff_type contact --cutoff 15 --select_part_map $select_map --scale_factor 1.275924

# Day 3
input_file="/mnt/extraids/OceanStor-SysCmn-2/chenweizhong/Project____InterChrom/ResultsV3/05_01____InterChrom_matrix/HUVEC_HT3dP6__Matrix__Resolution200000/chr2_chr7.matrix"
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/endoMT_analysis_NEW/chr2_chr7_200kb_T3.pdf"
python2.7 plot_iMARGI.py -i $input_file -o $output_file --isGlobal 0 --chr_row 2 --chr_col 7 --bin_size 200000 --data_type Contact --cutoff_type contact --cutoff 15 --select_part_map $select_map --scale_factor 1.87647

# Day 7
input_file_1="/mnt/extraids/OceanStor-SysCmn-2/chenweizhong/Project____InterChrom/ResultsV3/05_01____InterChrom_matrix/HUVEC_HT7dP7__Matrix__Resolution200000/chr2_chr7.matrix"
input_file_2="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered/HUVEC_T7d_P7__Matrix__Resolution200000/chr2_chr7.matrix"
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/endoMT_analysis_NEW/chr2_chr7_200kb_T7_MERGED.pdf"
python2.7 plot_iMARGI.py -i "["$input_file_1","$input_file_2"]" -o $output_file --isGlobal 0 --chr_row 2 --chr_col 7 --bin_size 200000 --data_type Contact --cutoff_type contact --cutoff 15 --select_part_map $select_map --scale_factor 1.517179



### Supplementary Figure 5b
select_map="[80000000,120000000,195000000,235000000]"

# Day 0
input_file_1="/mnt/extraids/OceanStor-SysCmn-2/chenweizhong/Project____InterChrom/ResultsV3/05_01____InterChrom_matrix/HUVECcontrol__Matrix__Resolution200000/chr7_chr2.matrix"
input_file_2="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_control_igm____R2filtered/HUVECcontrol__Matrix__Resolution200000/chr7_chr2.matrix"
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/endoMT_analysis_NEW/chr7_chr2_200kb_control_MERGED.pdf"
python2.7 plot_iMARGI.py -i "["$input_file_1","$input_file_2"]" -o $output_file --isGlobal 0 --chr_row 7 --chr_col 2 --bin_size 200000 --data_type Contact --cutoff_type contact --cutoff 15 --select_part_map $select_map --scale_factor 1.275924

# Day 3
input_file="/mnt/extraids/OceanStor-SysCmn-2/chenweizhong/Project____InterChrom/ResultsV3/05_01____InterChrom_matrix/HUVEC_HT3dP6__Matrix__Resolution200000/chr7_chr2.matrix"
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/endoMT_analysis_NEW/chr7_chr2_200kb_T3.pdf"
python2.7 plot_iMARGI.py -i $input_file -o $output_file --isGlobal 0 --chr_row 7 --chr_col 2 --bin_size 200000 --data_type Contact --cutoff_type contact --cutoff 15 --select_part_map $select_map --scale_factor 1.87647

# Day 7
input_file_1="/mnt/extraids/OceanStor-SysCmn-2/chenweizhong/Project____InterChrom/ResultsV3/05_01____InterChrom_matrix/HUVEC_HT7dP7__Matrix__Resolution200000/chr7_chr2.matrix"
input_file_2="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered/HUVEC_T7d_P7__Matrix__Resolution200000/chr7_chr2.matrix"
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/endoMT_analysis_NEW/chr7_chr2_200kb_T7_MERGED.pdf"
python2.7 plot_iMARGI.py -i "["$input_file_1","$input_file_2"]" -o $output_file --isGlobal 0 --chr_row 7 --chr_col 2 --bin_size 200000 --data_type Contact --cutoff_type contact --cutoff 15 --select_part_map $select_map --scale_factor 1.517179
