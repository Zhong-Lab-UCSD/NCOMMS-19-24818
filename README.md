# Dynamic changes in RNA-chromatin interactome promote endothelial dysfunction (NCOMMS-19-24818)

Repository of codes used in the paper "Dynamic changes in RNA-chromatin interactome promote endothelial dysfunction" (NCOMMS-19-24818).

## iMARGI

Data were processed using the [iMARGI pipeline](https://sysbio.ucsd.edu/imargi_pipeline/). The data analysis has been done using the following scripts:

- [``HUVEC_iMARGI.r``]() is the main R script used for analyzing and parsing data in order to generate data structures suitable for network or heatmap plotting.
- [``plot_network.r``]() contains a custom function built in-house using exploiting the R package igraph in order to plot and customize super enhancer networks, such as node dimension, color, labels, edge width, etc. The function is called in the main script ``HUVEC_iMARGI.r``.
- [``plot_iMARGI.py``]() used to plot iMARGI heatmaps.
- [``plot_iMARGI_SE_time.py``]() used to plot iMARGI super enhancer heatmap (SE x SE contact map) with the three samples on a side-by-side view.
- [``plot_iMARGI_maps.sh``]() is a bash script used to run the two python scripts above for plotting.


## Hi-C

All the codes used for Hi-C data analysis and visualization (Figure 3A in the manuscript) are open source and available from the published software [HiCtool (v2.1)](https://github.com/Zhong-Lab-UCSD/HiCtool).


## scRNA-seq

Two scripts were used for single-cell RNA-seq analysis, a bash script to run cellranger from 10X Genomics, and an R script to run Seurat functions for data analysis and visualization (Figure 2 in the manuscript):

- [``cellranger.sh``](./scRNAseq_scripts/cellranger.sh)
- [``HUVEC_scRNAseq.r``](./scRNAseq_scripts/HUVEC_scRNAseq.r)
