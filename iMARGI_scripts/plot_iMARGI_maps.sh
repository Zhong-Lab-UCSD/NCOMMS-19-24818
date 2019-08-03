############################################## Control
########## 1mb
input_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/4_2__WindowToWindow_matrix_of_RNAvsDNA/HUVECcontrollibrary2_igm_DistalInterChrom____resolution1000000.txt"
output_plot="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/normalized/control_MARGI_global_1mb.pdf"
python plot_iMARGI.py -i $input_file -o $output_plot --isGlobal 1 --bin_size 1000000 --data_type "Control iMARGI contact" --cutoff 20 --cutoff_type contact --scale_factor "1.517931"

############################################## T3d
########## 1mb
input_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/4_2__WindowToWindow_matrix_of_RNAvsDNA/HUVEC_H_T3d_P6_igm_DistalInterChrom____resolution1000000.txt"
output_plot="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/normalized/Day3_MARGI_global_1mb.pdf"
python plot_iMARGI.py -i $input_file -o $output_plot --isGlobal 1 --bin_size 1000000 --data_type "Day 3 iMARGI contact" --cutoff 20 --cutoff_type contact --scale_factor "1.87647"

############################################## T7d
########## 1mb
input_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/4_2__WindowToWindow_matrix_of_RNAvsDNA/HUVEC_H_T7d_P7_igm_DistalInterChrom____resolution1000000.txt"
output_plot="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/normalized/Day7_MARGI_global_1mb.pdf"
python plot_iMARGI.py -i $input_file -o $output_plot --isGlobal 1 --bin_size 1000000 --data_type "Day 7 iMARGI contact" --cutoff 20 --cutoff_type contact --scale_factor "1.998735"


################################### Super enhancer maps
control_SE='/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm/matrix_only_SE_sort.txt'
T3d_SE='/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered/matrix_only_SE_sort.txt'
T7d_SE='/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered/matrix_only_SE_sort.txt'
input_files="["$control_SE","$T3d_SE","$T7d_SE"]"

#### Plotting SE x SE contact map
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm/control_heatmap_SE_global.pdf"
python plot_iMARGI.py -i $control_SE -o $output_file --isGlobal 1 --bin_size 1 --data_type "control" --species "SE" --cutoff_type contact --cutoff 15
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered/T3d_heatmap_SE_global.pdf"
python plot_iMARGI.py -i $T3d_SE -o $output_file --isGlobal 1 --bin_size 1 --data_type "T3d_P6" --species "SE" --cutoff_type contact --cutoff 15
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered/T7d_heatmap_SE_global.pdf"
python plot_iMARGI.py -i $T7d_SE -o $output_file --isGlobal 1 --bin_size 1 --data_type "T7d_P7" --species "SE" --cutoff_type contact --cutoff 15

#### Plotting SE x SE contact map for chr2-chr7 and vice versa on a side-by-side view
output_file="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/endoMT_analysis_NEW/HUVEC_SE_chr2_chr7_time.pdf"
python plot_iMARGI_SE_time.py -i $input_files -o $output_file --isGlobal 1 --chr_row "[2,7]" --chr_col "[7,2]" --time_points "[Control,T3d_P6,T7d_P7]" --bin_size 1 --data_type "" --species SE --cutoff_type contact --cutoff 10
