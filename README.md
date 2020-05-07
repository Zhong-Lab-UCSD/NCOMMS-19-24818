# Dynamic changes in RNA-chromatin interactome promote endothelial dysfunction (NCOMMS-19-24818)

Repository of in-house codes and useful files used in the paper **Dynamic changes in RNA-chromatin interactome promote endothelial dysfunction**.

## iMARGI

### Data processing and computation of contact matrices

For this part, the following software are needed:
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [FastUniq](https://sourceforge.net/projects/fastuniq/)
- [STAR](https://github.com/alexdobin/STAR)
- [Samtools](http://www.htslib.org/)

These are the scripts used for this part:
- [``iMARGI_data_processing_MAIN.sh``](./iMARGI_scripts/raw_data_processing/iMARGI_data_processing_MAIN.sh) is the main script for raw data processing, from fastq files to the BEDPE file with the uniquely mapped read pairs.
- [``split_iMARGI_reads_chromosome_by_chromosome.sh``](./iMARGI_scripts/make_contact_matrices/split_iMARGI_reads_chromosome_by_chromosome.sh) is used to split the BEDPE file in 576 single BEDPE files, one per each chromosome pair.
- [``Make_WholeGenome_Matrix.sh``](./iMARGI_scripts/make_contact_matrices/Make_WholeGenome_Matrix.sh) is used to generate all the contact matrices from the BEDPE files, at a specified resolution and one per each chromosome pair.

### Data analysis and visualization

These are the scripts used for this part:

- [``HUVEC_iMARGI.r``](./iMARGI_scripts/HUVEC_iMARGI.r) is the main script used for analysis and parsing of the data in order to generate data structures suitable for network plotting, and generating data summary or reports.
- [``plot_network_function.r``](./iMARGI_scripts/plot_network_function.r) contains a custom function built in-house using the R package igraph in order to plot and customize super enhancer networks, such as node dimension, color, labels, edge width, etc. The function is called in the main script ``HUVEC_iMARGI.r``.
- [``HUVEC_iMARGI_replicate_2.r``](./iMARGI_scripts/HUVEC_iMARGI_replicate_2.r) is used for analysis and parsing of the data of biological replicate 2.
- [``HUVEC_iMARGI_summary.r``](./iMARGI_scripts/HUVEC_iMARGI_summary.r) is used for generating all the summary plots for iMARGI.
- [``plot_iMARGI_maps.sh``](./iMARGI_scripts/plot_iMARGI_maps.sh) is a bash script used to plot iMARGI contact matrices by running the python script ``plot_iMARGI.py``.
- [``plot_iMARGI.py``](./iMARGI_scripts/plot_iMARGI.py) contains the functions to plot iMARGI contact matrices.


## Hi-C

- The published software [HiCtool (v2.2)](https://github.com/Zhong-Lab-UCSD/HiCtool) was used for Hi-C data analysis and visualization, such as data pre-processing, data normalization, contact heatmap and correlation heatmap visualization, TAD and A/B compartment analyses.
- [``HUVEC_HiC.r``](./hic_scripts/HUVEC_HiC.r) serves to calculate general Hi-C data statistics, Measure of Concordance (MoC) of TAD boundaries between samples, average interaction frequency by genomic distance curves, and proportion of reads mapped within TADs.

## RNA-seq

RNA-seq analysis was mainly performed using the R package [DESeq2 (v1.24.0)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). These scripts were used for RNA-seq analysis:

- [``HUVEC_RNAseq.sh``](./RNAseq_scripts/HUVEC_RNAseq.sh) is a bash script used for alignment (performed using [STAR (v2.5.4b)](https://github.com/alexdobin/STAR)) and ``featureCounts`` from the package [Subread (v2.0.0)](http://subread.sourceforge.net/), to obtain the raw count data to input in DESeq2.
- [``HUVEC_RNAseq.r``](./RNAseq_scripts/HUVEC_RNAseq.r) contains the code from DESeq2 to perform the RNA-seq analysis.


## scRNA-seq

scRNA-seq analysis was mainly performed using the R package [Seurat (v2.3.4)](https://satijalab.org/seurat/). These scripts were used for single-cell RNA-seq analysis:

- [``cellranger.sh``](./scRNAseq_scripts/cellranger.sh) is a bash script used to run [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) from 10X Genomics both for in vitro and in vivo models.
- [``HUVEC_scRNAseq.r``](./scRNAseq_scripts/HUVEC_scRNAseq.r) is the main R script used for data analysis and visualization of the **in VITRO model**, based on functions from Seurat.
- [``HUVEC_scRNAseq_human_vascular.r``](./scRNAseq_scripts/HUVEC_scRNAseq_human_vascular.r) is the main R script used for data analysis and visualization of the **in VIVO model**, based on functions from Seurat.



