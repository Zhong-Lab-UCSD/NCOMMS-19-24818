# Dynamic changes in RNA-chromatin interactome promote endothelial dysfunction (NCOMMS-19-24818)

Repository of in-house codes and useful files used in the paper "Dynamic changes in RNA-chromatin interactome promote endothelial dysfunction".

## iMARGI

Data were processed using the [iMARGI pipeline](https://sysbio.ucsd.edu/imargi_pipeline/). Data analysis and visualization has been done using the following scripts:

- [``HUVEC_iMARGI.r``](./iMARGI_scripts/HUVEC_iMARGI.r) is the main R script used for analyzing and parsing data in order to generate data structures suitable for network or heatmap plotting, and generating data summary or reports.
- [``plot_network_function.r``](./iMARGI_scripts/plot_network_function.r) contains a custom function built in-house using the R package igraph in order to plot and customize super enhancer networks, such as node dimension, color, labels, edge width, etc. The function is called in the main script ``HUVEC_iMARGI.r``.
- [``plot_iMARGI.py``](./iMARGI_scripts/plot_iMARGI.py) is used to plot iMARGI contact heatmaps.
- [``plot_iMARGI_SE_time.py``](./iMARGI_scripts/plot_iMARGI_SE_time.py) is used to plot iMARGI super enhancer (SE) contact heatmap (SEs x SEs contact map) with the three samples on a side-by-side view.
- [``plot_iMARGI_maps.sh``](./iMARGI_scripts/plot_iMARGI_maps.sh) is a bash script used to run the two python scripts above for plotting iMARGI heatmaps.


## Hi-C

Hi-C data analysis and visualization were mainly performed using the published software [HiCtool (v2.2)](https://github.com/Zhong-Lab-UCSD/HiCtool). [``HUVEC_hictool.sh``](./hic_scripts/HUVEC_hictool.sh) contains the HiCtool commands for: pre-processing raw Hi-C data, data normalization, contact heatmap and observed/expected heatmap visualization, TAD and A/B compartment analysis.

General Hi-C data statistics, Measure of Concordance (MoC) of TAD boundaries between samples, average interaction frequency by genomic distance curves, proportion of reads mapped within TADs, were calculated using the R script [``HUVEC_HiC.r``](./hic_scripts/HUVEC_HiC.r).


## RNA-seq

These scripts were used for RNA-seq analysis:

- [HUVEC_RNAseq.sh](./RNAseq_scripts/HUVEC_RNAseq.sh) is a bash script used for alignment job and ``featureCounts``, to obtain the input data for DESeq2 in R.
- [``HUVEC_RNAseq.r``](./RNAseq_scripts/HUVEC_RNAseq.r) contains the code from the package [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to perform RNA-seq analysis.


## scRNA-seq

These scripts were used for single-cell RNA-seq analysis:

- [``cellranger.sh``](./scRNAseq_scripts/cellranger.sh) is a bash script used to run [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) from 10X Genomics both for in vitro and in vivo models.
- [``HUVEC_scRNAseq.r``](./scRNAseq_scripts/HUVEC_scRNAseq.r) is the main R script used for data analysis and visualization of the **in VITRO model**, based on functions from the package [Seurat (2.3.4)](https://satijalab.org/seurat/).
- [``HUVEC_scRNAseq_human_vascular.r``](./scRNAseq_scripts/HUVEC_scRNAseq_human_vascular.r) is the main R script used for data analysis and visualization of the **in VIVO model**, based on functions from the package [Seurat (2.3.4)](https://satijalab.org/seurat/).



