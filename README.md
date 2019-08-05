# Dynamic changes in RNA-chromatin interactome promote endothelial dysfunction (NCOMMS-19-24818)

Repository of codes used in the paper "Dynamic changes in RNA-chromatin interactome promote endothelial dysfunction" (NCOMMS-19-24818).

## iMARGI

Data were processed using the [iMARGI pipeline](https://sysbio.ucsd.edu/imargi_pipeline/). Data analysis and visualization has been done using the following scripts:

- [``HUVEC_iMARGI.r``](./iMARGI_scripts/HUVEC_iMARGI.r) is the main R script used for analyzing and parsing data in order to generate data structures suitable for network or heatmap plotting, and generating data summary or reports.
- [``plot_network_function.r``](./iMARGI_scripts/plot_network_function.r) contains a custom function built in-house using the R package igraph in order to plot and customize super enhancer networks, such as node dimension, color, labels, edge width, etc. The function is called in the main script ``HUVEC_iMARGI.r``.
- [``plot_iMARGI.py``](./iMARGI_scripts/plot_iMARGI.py) is used to plot iMARGI contact heatmaps.
- [``plot_iMARGI_SE_time.py``](./iMARGI_scripts/plot_iMARGI_SE_time.py) is used to plot iMARGI super enhancer (SE) contact heatmap (SEs x SEs contact map) with the three samples on a side-by-side view.
- [``plot_iMARGI_maps.sh``](./iMARGI_scripts/plot_iMARGI_maps.sh) is a bash script used to run the two python scripts above for plotting iMARGI heatmaps.


## Hi-C

All the codes used for Hi-C data analysis and visualization are open source and available from the published software [HiCtool (v2.1)](https://github.com/Zhong-Lab-UCSD/HiCtool).


## scRNA-seq

Two scripts were used for single-cell RNA-seq analysis:

- [``cellranger.sh``](./scRNAseq_scripts/cellranger.sh) is a bash script used to run [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) from 10X Genomics.
- [``HUVEC_scRNAseq.r``](./scRNAseq_scripts/HUVEC_scRNAseq.r) is the main R script used for data analysis and visualization, based on functions from the package [Seurat (2.3.4)](https://satijalab.org/seurat/).
