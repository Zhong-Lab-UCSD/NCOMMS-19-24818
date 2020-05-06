################ scRNA-seq analysis with Seurat

annotation <- read.table('/dataOS/rcalandrelli/MARGI/Homo_sapiens.GRCh38.84.chr.gtf_to_geneTable.tsv', stringsAsFactors = F)
colnames(annotation)=annotation[1,]
annotation=annotation[-1,]
annotation$start = as.numeric(annotation$start)
annotation$end = as.numeric(annotation$end)

### Load the datasets
library(Seurat) # version 2.3.4
library(ggbio)
library(dplyr)
library(cowplot)

min_cells = 0
min_genes = 0

########## normal_1

# normal_1_1
normal_EC_1.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/counts/count_normal_EC/outs/filtered_feature_bc_matrix")
normal_EC_1 <- CreateSeuratObject(raw.data = normal_EC_1.data, min.cells = min_cells, min.genes = min_genes, project = "normal_EC_1")
normal_EC_1@meta.data$sample <- "n.1"
normal_EC_1@meta.data$sample_type <- "normal"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = normal_EC_1@data), value = TRUE)
percent.mito <- Matrix::colSums(normal_EC_1@raw.data[mito.genes, ])/Matrix::colSums(normal_EC_1@raw.data)
normal_EC_1 <- AddMetaData(object = normal_EC_1, metadata = percent.mito, col.name = "percent.mito")

# normal_1_2
normal_SMC_1.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/counts/count_normal_SMC/outs/filtered_feature_bc_matrix")
normal_SMC_1 <- CreateSeuratObject(raw.data = normal_SMC_1.data, min.cells = min_cells, min.genes = min_genes, project = "normal_SMC_1")
normal_SMC_1@meta.data$sample <- "n.1"
normal_SMC_1@meta.data$sample_type <- "normal"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = normal_SMC_1@data), value = TRUE)
percent.mito <- Matrix::colSums(normal_SMC_1@raw.data[mito.genes, ])/Matrix::colSums(normal_SMC_1@raw.data)
normal_SMC_1 <- AddMetaData(object = normal_SMC_1, metadata = percent.mito, col.name = "percent.mito")

# Merge the data
normal_1 = MergeSeurat(normal_EC_1,normal_SMC_1, do.normalize = F, project = "normal_1")
normal_1@meta.data$sample <- "n.1"
normal_1@meta.data$sample_type <- "normal"

########## normal_2

# normal_2_1
normal_2_SA48770.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/counts/count_normal_2_SA48770/outs/filtered_feature_bc_matrix")
normal_2_SA48770 <- CreateSeuratObject(raw.data = normal_2_SA48770.data, min.cells = min_cells, min.genes = min_genes, project = "normal_2_SA48770")
normal_2_SA48770@meta.data$sample <- "n.2"
normal_2_SA48770@meta.data$sample_type <- "normal"
mito.genes <- grep(pattern = "^MT-", x = rownames(x = normal_2_SA48770@data), value = TRUE)
percent.mito <- Matrix::colSums(normal_2_SA48770@raw.data[mito.genes, ])/Matrix::colSums(normal_2_SA48770@raw.data)
normal_2_SA48770 <- AddMetaData(object = normal_2_SA48770, metadata = percent.mito, col.name = "percent.mito")

# normal_2_2
normal_2_SA48771.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/counts/count_normal_2_SA48771/outs/filtered_feature_bc_matrix")
normal_2_SA48771 <- CreateSeuratObject(raw.data = normal_2_SA48771.data, min.cells = min_cells, min.genes = min_genes, project = "normal_2_SA48771")
normal_2_SA48771@meta.data$sample <- "n.2"
normal_2_SA48771@meta.data$sample_type <- "normal"
mito.genes <- grep(pattern = "^MT-", x = rownames(x = normal_2_SA48771@data), value = TRUE)
percent.mito <- Matrix::colSums(normal_2_SA48771@raw.data[mito.genes, ])/Matrix::colSums(normal_2_SA48771@raw.data)
normal_2_SA48771 <- AddMetaData(object = normal_2_SA48771, metadata = percent.mito, col.name = "percent.mito")

# normal_2_3
normal_2_SA48772.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/counts/count_normal_2_SA48772/outs/filtered_feature_bc_matrix")
normal_2_SA48772 <- CreateSeuratObject(raw.data = normal_2_SA48772.data, min.cells = min_cells, min.genes = min_genes, project = "normal_2_SA48772")
normal_2_SA48772@meta.data$sample <- "n.2"
normal_2_SA48772@meta.data$sample_type <- "normal"
mito.genes <- grep(pattern = "^MT-", x = rownames(x = normal_2_SA48772@data), value = TRUE)
percent.mito <- Matrix::colSums(normal_2_SA48772@raw.data[mito.genes, ])/Matrix::colSums(normal_2_SA48772@raw.data)
normal_2_SA48772 <- AddMetaData(object = normal_2_SA48772, metadata = percent.mito, col.name = "percent.mito")

# normal_2_4
normal_2_SA48773.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/counts/count_normal_2_SA48773/outs/filtered_feature_bc_matrix")
normal_2_SA48773 <- CreateSeuratObject(raw.data = normal_2_SA48773.data, min.cells = min_cells, min.genes = min_genes, project = "normal_2_SA48773")
normal_2_SA48773@meta.data$sample <- "n.2"
normal_2_SA48773@meta.data$sample_type <- "normal"
mito.genes <- grep(pattern = "^MT-", x = rownames(x = normal_2_SA48773@data), value = TRUE)
percent.mito <- Matrix::colSums(normal_2_SA48773@raw.data[mito.genes, ])/Matrix::colSums(normal_2_SA48773@raw.data)
normal_2_SA48773 <- AddMetaData(object = normal_2_SA48773, metadata = percent.mito, col.name = "percent.mito")

# Merging data
normal_2 = MergeSeurat(normal_2_SA48770,normal_2_SA48771, do.normalize = F, project = "normal_2", add.cell.id1 = "normal.2.1", add.cell.id2 = "normal.2.2")
normal_2 = MergeSeurat(normal_2,normal_2_SA48772, do.normalize = F, project = "normal_2", add.cell.id2 = "normal.2.3")
normal_2 = MergeSeurat(normal_2,normal_2_SA48773, do.normalize = F, project = "normal_2", add.cell.id2 = "normal.2.4")
normal_2@meta.data$sample <- "n.2"
normal_2@meta.data$sample_type <- "normal"

########## T2D_1

# T2D_1_1
T2D_EC_1.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/counts/count_T2D_EC/outs/filtered_feature_bc_matrix")
T2D_EC_1 <- CreateSeuratObject(raw.data = T2D_EC_1.data, min.cells = min_cells, min.genes = min_genes, project = "T2D_EC_1")
T2D_EC_1@meta.data$sample <- "T2D.1"
T2D_EC_1@meta.data$sample_type <- "T2D"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T2D_EC_1@data), value = TRUE)
percent.mito <- Matrix::colSums(T2D_EC_1@raw.data[mito.genes, ])/Matrix::colSums(T2D_EC_1@raw.data)
T2D_EC_1 <- AddMetaData(object = T2D_EC_1, metadata = percent.mito, col.name = "percent.mito")

# T2D_1_2
T2D_SMC_1.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/counts/count_T2D_SMC/outs/filtered_feature_bc_matrix")
T2D_SMC_1 <- CreateSeuratObject(raw.data = T2D_SMC_1.data, min.cells = min_cells, min.genes = min_genes, project = "T2D_SMC_1")
T2D_SMC_1@meta.data$sample <- "T2D.1"
T2D_SMC_1@meta.data$sample_type <- "T2D"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T2D_SMC_1@data), value = TRUE)
percent.mito <- Matrix::colSums(T2D_SMC_1@raw.data[mito.genes, ])/Matrix::colSums(T2D_SMC_1@raw.data)
T2D_SMC_1 <- AddMetaData(object = T2D_SMC_1, metadata = percent.mito, col.name = "percent.mito")

# Merge the data
T2D_1 = MergeSeurat(T2D_EC_1,T2D_SMC_1, do.normalize = F, project = "T2D_1", add.cell.id1 = "T2D.1.1", add.cell.id2 = "T2D.1.2")
T2D_1@meta.data$sample <- "T2D.1"
T2D_1@meta.data$sample_type <- "T2D"

########## T2D_2

# T2D_2_1
T2D_EC_2.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/counts/count_T2D_EC_2/outs/filtered_feature_bc_matrix")
T2D_EC_2 <- CreateSeuratObject(raw.data = T2D_EC_2.data, min.cells = min_cells, min.genes = min_genes, project = "T2D_EC_2")
T2D_EC_2@meta.data$sample <- "T2D.2"
T2D_EC_2@meta.data$sample_type <- "T2D"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T2D_EC_2@data), value = TRUE)
percent.mito <- Matrix::colSums(T2D_EC_2@raw.data[mito.genes, ])/Matrix::colSums(T2D_EC_2@raw.data)
T2D_EC_2 <- AddMetaData(object = T2D_EC_2, metadata = percent.mito, col.name = "percent.mito")

# T2D_2_2
T2D_SMC_2.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/counts/count_T2D_SMC_2/outs/filtered_feature_bc_matrix")
T2D_SMC_2 <- CreateSeuratObject(raw.data = T2D_SMC_2.data, min.cells = min_cells, min.genes = min_genes, project = "T2D_SMC_2")
T2D_SMC_2@meta.data$sample <- "T2D.2"
T2D_SMC_2@meta.data$sample_type <- "T2D"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T2D_SMC_2@data), value = TRUE)
percent.mito <- Matrix::colSums(T2D_SMC_2@raw.data[mito.genes, ])/Matrix::colSums(T2D_SMC_2@raw.data)
T2D_SMC_2 <- AddMetaData(object = T2D_SMC_2, metadata = percent.mito, col.name = "percent.mito")

# Merge the data
T2D_2 = MergeSeurat(T2D_EC_2,T2D_SMC_2, do.normalize = F, project = "T2D_2", add.cell.id1 = "T2D.2.1", add.cell.id2 = "T2D.2.2")
T2D_2@meta.data$sample <- "T2D.2"
T2D_2@meta.data$sample_type <- "T2D"

### Filter cells and normalize

thres_cells = 0.98
thres_mito = 0.20
thres_min_genes = 300

thresh.normal_1.mito = thres_mito
thresh.normal_1.genes = round(quantile(normal_1@meta.data$nGene,thres_cells))
normal_1_filter <- FilterCells(object = normal_1, subset.names = c("nGene", "percent.mito"), low.thresholds = c(thres_min_genes, -Inf), high.thresholds = c(thresh.normal_1.genes, thresh.normal_1.mito))

thresh.normal_2.mito = thres_mito
thresh.normal_2.genes = round(quantile(normal_2@meta.data$nGene,thres_cells))
normal_2_filter <- FilterCells(object = normal_2, subset.names = c("nGene", "percent.mito"), low.thresholds = c(thres_min_genes, -Inf), high.thresholds = c(thresh.normal_2.genes, thresh.normal_2.mito))

thresh.T2D_1.mito = thres_mito
thresh.T2D_1.genes = round(quantile(T2D_1@meta.data$nGene,thres_cells))
T2D_1_filter <- FilterCells(object = T2D_1, subset.names = c("nGene", "percent.mito"), low.thresholds = c(thres_min_genes, -Inf), high.thresholds = c(thresh.T2D_1.genes, thresh.T2D_1.mito))

thresh.T2D_2.mito = thres_mito
thresh.T2D_2.genes = round(quantile(T2D_2@meta.data$nGene,thres_cells))
T2D_2_filter <- FilterCells(object = T2D_2, subset.names = c("nGene", "percent.mito"), low.thresholds = c(thres_min_genes, -Inf), high.thresholds = c(thresh.T2D_2.genes, thresh.T2D_2.mito))


### Merge the object together and normalize
sc_full_data = MergeSeurat(normal_1_filter, normal_2_filter, add.cell.id1 = "n.1", add.cell.id2 = "n.2")
sc_full_data = MergeSeurat(sc_full_data, T2D_1_filter, add.cell.id2 = "T2D.1")
sc_full_data = MergeSeurat(sc_full_data, T2D_2_filter, add.cell.id2 = "T2D.2")


### Detection of highly variable genes across the single cells
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/result/gene_dispersion_plot.png", width = 800, height = 500, units = "px")
sc_full_data <- FindVariableGenes(object = sc_full_data, mean.function = ExpMean, dispersion.function = LogVMR, 
                                  x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.5, 
                                  do.text = FALSE, do.contour = FALSE)
dev.off()

length(sc_full_data@var.genes) # number oh high variable genes detected according to the thresholds above
# Select the HVG as the top 1000 genes sorted by VMR
hv.genes = head(rownames(sc_full_data@hvg.info),1000)

### Scaling the data and removing unwanted sources of variation
sc_full_data <- ScaleData(object = sc_full_data, genes.use = hv.genes, 
                          vars.to.regress = c("nUMI","percent.mito"), do.par = TRUE, num.cores = 8)

### Perform linear dimensional reduction
sc_full_data <- RunPCA(object = sc_full_data, pc.genes = hv.genes, pcs.compute = 50,
                       do.print = FALSE, pcs.print = 1:5, genes.print = 5)

# Elbow plot
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/result/pc_elbow_plot.png", width = 600, height = 500, units = "px")
PCElbowPlot(object = sc_full_data, num.pc = 50)
dev.off()


### Clustering
top_PCs = 20
sc_full_data <- FindClusters(object = sc_full_data, reduction.type = "pca", dims.use = 1:top_PCs, resolution = 0.4, print.output = 0, save.SNN = TRUE)

# Save cell cluster data
cell_clusters = cbind(rownames(sc_full_data@meta.data), sc_full_data@meta.data[,"res.0.4"])
write.table(cell_clusters,"/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/result/cells_clusters.txt", row.names = F, col.names = T, sep ="\t", quote = F)

### t-SNE
sc_full_data <- RunTSNE(object = sc_full_data, dims.use = 1:top_PCs, reduction.use = "pca", do.fast = TRUE)

# Saving tSNE coordinates
tsne_coord = sc_full_data@dr$tsne@cell.embeddings
write.table(tsne_coord,"/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/result/tsne_coordinates.txt", row.names = T, col.names = T, sep = '\t', quote = F)

# Plots
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/result/Supplementary_Figure_11a-b.png", width = 2500, height = 2000, units = "px", res = 200)
p1 <- TSNEPlot(sc_full_data, do.return = T, pt.size = 0.5, group.by = "sample") + 
  theme(text = element_text(size=24),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24))
p2 <- TSNEPlot(sc_full_data, do.return = T, pt.size = 0.5, do.label = T, label.size = 6) + 
  theme(text = element_text(size=24),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24))
plot_grid(p1, p2, nrow=1, ncol=2)
dev.off()


# Plot CDH5 expression
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/result/Supplementary_Figure_11c.png", width = 4, height = 4, units = "in", res = 300)
FeaturePlot(object = sc_full_data, features.plot = "CDH5", nCol = 2,
            min.cutoff = "q10", max.cutoff = "q90", cols.use = c("white", "red"), pt.size = 0.5,
            dark.theme = TRUE, no.legend = FALSE)
dev.off()

heatmap_intensities = as.matrix(sc_full_data@data["CDH5",])
write.table(heatmap_intensities,"/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/result/gene_heatmap_intensities.txt", row.names = F, col.names = T, sep = '\t', quote = F)


### DE gene analysis 

my.data=FetchData(sc_full_data,c("ident","PC1","nGene","orig.ident"))

normal_cells_1 = rownames(my.data[which(my.data$orig.ident == "normal_EC_1" | my.data$orig.ident == "normal_SMC_1"),])
normal_cells_2 = rownames(my.data[which(my.data$orig.ident %in% c("normal_2_SA48770","normal_2_SA48771","normal_2_SA48772","normal_2_SA48773")),])
T2D_cells_1 = rownames(my.data[which(my.data$orig.ident == "T2D_EC_1" | my.data$orig.ident == "T2D_SMC_1"),])
T2D_cells_2 = rownames(my.data[which(my.data$orig.ident == "T2D_EC_2" | my.data$orig.ident == "T2D_SMC_2"),])

normal_cells = c(normal_cells_1,normal_cells_2)
T2D_cells = c(T2D_cells_1,T2D_cells_2)

### Extract only EC from clusters and gene expression of markers
EC_normal = rownames(sc_full_data@meta.data[which(sc_full_data@meta.data$res.0.4 %in% c(0,1,6) & 
                                                    sc_full_data@meta.data$sample_type == "normal"),])
EC_T2D = rownames(sc_full_data@meta.data[which(sc_full_data@meta.data$res.0.4 %in% c(0,1,6) & 
                                                 sc_full_data@meta.data$sample_type == "T2D"),])

# Change ident to cells
sc_full_data_orig = SetAllIdent(sc_full_data, "orig.ident")
sc_full_data_orig = SetIdent(sc_full_data_orig, cells.use = EC_normal, ident.use = "EC_normal")
sc_full_data_orig = SetIdent(sc_full_data_orig, cells.use = EC_T2D, ident.use = "EC_T2D")

# T2D over normal
de_genes_T2D_normal <- FindMarkers(sc_full_data_orig, ident.1 = "EC_T2D", ident.2 = "EC_normal", logfc.threshold = 0.25)
de_genes_T2D_normal = de_genes_T2D_normal[order(de_genes_T2D_normal$avg_logFC, decreasing = T),]
de_genes_T2D_normal = cbind(de_genes_T2D_normal,seq(1,nrow(de_genes_T2D_normal)))
colnames(de_genes_T2D_normal)[6] = "rank"

de_genes_T2D_normal[which(rownames(de_genes_T2D_normal)=="LINC00607"),]

temp = cbind(rownames(de_genes_T2D_normal),de_genes_T2D_normal)
colnames(temp)[1] = 'gene'
write.table(temp,"/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/result/de_genes_T2D_normal.txt",row.names = F,col.names = T,quote = F, sep = '\t')

de_genes_T2D_normal = read.table("/dataOS/rcalandrelli/MARGI/RNAseq_human_vascular/result/de_genes_T2D_normal.txt", stringsAsFactors = F)
colnames(de_genes_T2D_normal) = de_genes_T2D_normal[1,]
de_genes_T2D_normal = de_genes_T2D_normal[-1,]
rownames(de_genes_T2D_normal) = de_genes_T2D_normal[,1]
de_genes_T2D_normal$avg_logFC = as.numeric(de_genes_T2D_normal$avg_logFC)

# Plot DE genes
sc_full_data_plot = SetAllIdent(sc_full_data, "orig.ident")
sc_full_data_plot = SetIdent(sc_full_data_plot, cells.use = EC_normal, ident.use = "EC_normal")
sc_full_data_plot = SetIdent(sc_full_data_plot, cells.use = EC_T2D, ident.use = "EC_T2D")

sc_full_data_plot <- ScaleData(object = sc_full_data_plot, genes.use = rownames(de_genes_T2D_normal), 
                               vars.to.regress = c("nUMI","percent.mito"), do.par = TRUE, num.cores = 4)

# Heatmap plotted manually by binning single cells
make_binned_sc_heatmap<-function(data, # data to be binned
                                 cell_binned, # how many single cells in each bin
                                 ordering, # TRUE to order single cell before binning
                                 which_gene # gene whose expression levels are used to order single cells
                                 ){
  # Healty ECs
  data_normal = data[,EC_normal]
  if (ordering == T){
    data_normal = data_normal[,order(data_normal[which_gene,])]
  }
  data_normal_binned = matrix(nrow=nrow(data_normal),ncol=as.integer(ncol(data_normal)/cell_binned))
  k = 0
  i = 1
  while (k < ncol(data_normal)){
    ind_start = k + 1
    ind_end = k + cell_binned
    if (ind_end < ncol(data_normal)){
      data_normal_binned[,i] = rowMeans(data_normal[,ind_start:ind_end])
    }
    k = k + cell_binned
    i = i + 1
  }
  rownames(data_normal_binned) = rownames(data_normal)
  colnames(data_normal_binned) = paste0("HC.bin.",seq(1,ncol(data_normal_binned)))
  
  # T2D ECs
  data_T2D = data[,EC_T2D]
  if (ordering == T){
    data_T2D = data_T2D[,order(data_T2D[which_gene,])]
  }
  data_T2D_binned = matrix(nrow=nrow(data_T2D),ncol=as.integer(ncol(data_T2D)/cell_binned))
  k = 0
  i = 1
  while (k < ncol(data_T2D)){
    ind_start = k + 1
    ind_end = k + cell_binned
    if (ind_end < ncol(data_T2D)){
      data_T2D_binned[,i] = rowMeans(data_T2D[,ind_start:ind_end])
    }
    k = k + cell_binned
    i = i + 1
  }
  rownames(data_T2D_binned) = rownames(data_T2D)
  colnames(data_T2D_binned) = paste0("T2D.bin.",seq(1,ncol(data_T2D_binned)))
  
  ### Merge the data
  data_merged = cbind(data_normal_binned, NA, NA, data_T2D_binned)
  return(data_merged)
}

temp = make_binned_sc_heatmap(sc_full_data_plot@scale.data[temp_induced[,1],], 25, TRUE, "SERPINE1")

temp[which(temp < -2.5, arr.ind = T)] = -2.5
temp[which(temp > 2.5, arr.ind = T)] = 2.5

pdf("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/result/heatmap_binned_25.pdf", height=8, width=5)
pheatmap(temp,
         color = colorRampPalette(c("blue", "yellow"))(100),
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         show_colnames = F,
         na_col = "white")
dev.off()

### Odds ratio analysis

total_control = 65879147
total_T7d = 50031654

super_enhancers_sort_annotated_new = read.table("/dataOS/rcalandrelli/MARGI/super_enhancers_sort_annotated_new.txt", stringsAsFactors = F)
colnames(super_enhancers_sort_annotated_new) = super_enhancers_sort_annotated_new[1,]
super_enhancers_sort_annotated_new = super_enhancers_sort_annotated_new[-1,]

# Day 0
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm'
heatmap_control = read.table(paste0(directory,"/matrix_only_SE_sort.txt"), stringsAsFactors = F)
rownames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new

# Day 7
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered'
heatmap_T7d = read.table(paste0(directory,"/matrix_only_SE_sort.txt"), stringsAsFactors = F)
rownames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new

# LINC00607 analysis
LINC00607_signal_control = (heatmap_control[145,] + heatmap_control[,145] + 1) / total_control # added pseudocounts
LINC00607_signal_T7d = (heatmap_T7d[145,] + heatmap_T7d[,145] + 1) / total_T7d # added pseudocounts

se_with_genes = super_enhancers_sort_annotated_new[which(super_enhancers_sort_annotated_new$SE_genes != '.'),]

LINC00607_signal_control = LINC00607_signal_control[se_with_genes$SE_index_new]
LINC00607_signal_T7d = LINC00607_signal_T7d[se_with_genes$SE_index_new]

### Addind the number of iMARGI read pairs (not normalized) in Day 0 and Day 7
read_pairs_control = heatmap_control[145,] + heatmap_control[,145]
read_pairs_control = read_pairs_control[se_with_genes$SE_index_new]

read_pairs_T7d = heatmap_T7d[145,] + heatmap_T7d[,145]
read_pairs_T7d = read_pairs_T7d[se_with_genes$SE_index_new]

read_pairs_control_RNA_DNA = heatmap_control[145,]
read_pairs_control_RNA_DNA = read_pairs_control_RNA_DNA[se_with_genes$SE_index_new]

read_pairs_T7d_RNA_DNA = heatmap_T7d[145,]
read_pairs_T7d_RNA_DNA = read_pairs_T7d_RNA_DNA[se_with_genes$SE_index_new]

logfc_LINC00607_signal = log(LINC00607_signal_T7d/LINC00607_signal_control)

se_with_genes_LINC = cbind(se_with_genes,t(logfc_LINC00607_signal),t(read_pairs_control),t(read_pairs_T7d),t(read_pairs_control_RNA_DNA),t(read_pairs_T7d_RNA_DNA))
colnames(se_with_genes_LINC)[7:11] = c("logfc_LINC00607_signal","read_pairs_control","read_pairs_T7d","read_pairs_control_RNA_DNA","read_pairs_T7d_RNA_DNA")
se_with_genes_LINC = se_with_genes_LINC[order(se_with_genes_LINC$logfc_LINC00607_signal, decreasing = T),]
se_with_genes_LINC = cbind(se_with_genes_LINC,seq(1:nrow(se_with_genes_LINC)))
colnames(se_with_genes_LINC)[12] = "rank"

### Extract each gene and keep track of SE and rank
genes_ranked = matrix(nrow=0,ncol=8)
colnames(genes_ranked) = c("gene","SE_index_new","logfc_LINC00607_signal","read_pairs_control","read_pairs_T7d","read_pairs_control_RNA_DNA","read_pairs_T7d_RNA_DNA","rank")
for (i in 1:nrow(se_with_genes_LINC)){
  temp_genes = strsplit(se_with_genes_LINC[i,"SE_genes"],";")[[1]]
  temp_mat = cbind(temp_genes,
                   se_with_genes_LINC[i,"SE_index_new"],
                   se_with_genes_LINC[i,"logfc_LINC00607_signal"],
                   se_with_genes_LINC[i,"read_pairs_control"],
                   se_with_genes_LINC[i,"read_pairs_T7d"],
                   se_with_genes_LINC[i,"read_pairs_control_RNA_DNA"],
                   se_with_genes_LINC[i,"read_pairs_T7d_RNA_DNA"],
                   se_with_genes_LINC[i,"rank"])
  genes_ranked = rbind(genes_ranked,temp_mat)
}

# In genes_ranked each gene could be present more than once (with a different logfc_LINC00607_signal) if it overlaps multiple super enhancers
genes_ranked = data.frame(genes_ranked)
genes_ranked$gene = as.character(genes_ranked$gene)
genes_ranked$SE_index_new = as.numeric(as.character(genes_ranked$SE_index_new))
genes_ranked$logfc_LINC00607_signal = as.numeric(as.character(genes_ranked$logfc_LINC00607_signal))
genes_ranked$read_pairs_control = as.numeric(as.character(genes_ranked$read_pairs_control))
genes_ranked$read_pairs_T7d = as.numeric(as.character(genes_ranked$read_pairs_T7d))
genes_ranked$read_pairs_control_RNA_DNA = as.numeric(as.character(genes_ranked$read_pairs_control_RNA_DNA))
genes_ranked$read_pairs_T7d_RNA_DNA = as.numeric(as.character(genes_ranked$read_pairs_T7d_RNA_DNA))
genes_ranked$rank = as.numeric(as.character(genes_ranked$rank))

genes_ranked_sc = genes_ranked[which(genes_ranked$gene %in% rownames(sc_full_data@data)),]

logfc_categories_and_odds_ratio <- function(x,normal_cells,T2D_cells){
  my_gene = x[1]
  temp = sc_full_data@data[c("LINC00607", my_gene),normal_cells]
  temp_0_0_normal = sum(colSums(temp)==0) # 0,0
  temp_p_p_normal = sum(temp[1,]*temp[2,] != 0) # +,+
  temp_p_0_normal = sum(xor(temp[1,],temp[2,]) & temp[1,]>0) # +,0
  temp_0_p_normal = sum(xor(temp[1,],temp[2,]) & temp[2,]>0) # 0,+
  mat_normal = (matrix(c(temp_0_p_normal,temp_0_0_normal,temp_p_p_normal,temp_p_0_normal),2,2)+1)/length(normal_cells)
  
  temp = sc_full_data@data[c("LINC00607", my_gene),T2D_cells]
  temp_0_0_T2D = sum(colSums(temp)==0) # 0,0
  temp_p_p_T2D = sum(temp[1,]*temp[2,] != 0) # +,+
  temp_p_0_T2D = sum(xor(temp[1,],temp[2,]) & temp[1,]>0) # +,0
  temp_0_p_T2D = sum(xor(temp[1,],temp[2,]) & temp[2,]>0) # 0,+
  mat_T2D = (matrix(c(temp_0_p_T2D,temp_0_0_T2D,temp_p_p_T2D,temp_p_0_T2D),2,2)+1)/length(T2D_cells)
  
  mat = log(mat_T2D/mat_normal)
  df.mat = expand.grid(x = c("0", "+"), y = c("0", "+"))
  df.mat$logfc = c(mat[2,1],mat[2,2],mat[1,1],mat[1,2])
  df.mat$cells_normal = c(temp_0_0_normal,temp_p_0_normal,temp_0_p_normal,temp_p_p_normal)
  df.mat$cells_T2D = c(temp_0_0_T2D,temp_p_0_T2D,temp_0_p_T2D,temp_p_p_T2D)
  odds_ratio = (mat_T2D[1,2]*mat_normal[2,2])/(mat_normal[1,2]*mat_T2D[2,2])
  return(c(t(df.mat[,3:5])[1,],t(df.mat[,3:5])[2,],t(df.mat[,3:5])[3,],odds_ratio))
}

### Using all the cells
temp = data.frame(t(apply(genes_ranked_sc,1,logfc_categories_and_odds_ratio,normal_cells,T2D_cells)))

colnames(temp) = c("logfc_0_0","logfc_p_0","logfc_0_p","logfc_p_p",
                   "cells_normal_0_0","cells_normal_p_0","cells_normal_0_p","cells_normal_p_p",
                   "cells_T2D_0_0","cells_T2D_p_0","cells_T2D_0_p","cells_T2D_p_p",
                   "odds_ratio")
genes_ranked_sc_logfcCateg = cbind(genes_ranked_sc,temp)

# Remove one duplicate line of MALAT1
temp_ind = which(genes_ranked_sc_logfcCateg$gene=="MALAT1")[1]
genes_ranked_sc_logfcCateg = genes_ranked_sc_logfcCateg[-temp_ind,]

### Add information of the gene super enhancer to select inter-intra chromosomal respect to LINC00607 (chr2)
genes_ranked_sc_logfcCateg = merge(genes_ranked_sc_logfcCateg,super_enhancers_sort_annotated_new,by="SE_index_new")

### Plot only odds ratios in a violin plot
hubs_interacting_with_linc = c(as.numeric(indexes_T7d_hubs[which(indexes_T7d_hubs[,1]==145),2]),
  as.numeric(indexes_T7d_hubs[which(indexes_T7d_hubs[,2]==145),1]))

temp = genes_ranked_sc_logfcCateg[which(genes_ranked_sc_logfcCateg$SE_index_new %in% hubs_interacting_with_linc),]
temp = merge(temp, annotation[,c("gene_name","gene_biotype")], by.x = "gene", by.y = "gene_name", all.x = T)
temp = temp[which(temp$gene_biotype %in% c("protein_coding","lincRNA")),]
temp = temp[order(temp$odds_ratio, decreasing = T),]
write.table(temp,"/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/result/violin_plot_genes.txt", row.names = F, col.names = T, sep = '\t', quote = F)

temp$gene_label = "Other gene"
temp[which(temp$gene == "SERPINE1"),"gene_label"] = "SERPINE1"
temp[which(temp$gene == "RUNX1"),"gene_label"] = "RUNX1"
temp[which(temp$gene == "TRIO"),"gene_label"] = "TRIO"

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_human_vascular/result/Figure_5c.png", height=1000, width=1000, units = "px", res = 200)
ggplot(temp, aes(x=0, y=odds_ratio)) + 
  geom_violin() +
  geom_point(position = position_jitter(0.3), aes(shape = gene_label), size = 2) +
  scale_shape_manual(values = c(4,2,5,0)) +
  scale_y_continuous(position = "right", breaks = c(0,2,4,6), labels = c("0","2","4","6"), limits = c(0,6)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text=element_text(size=14),
        legend.title = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.y = element_text(angle = 0)) +
  labs(y = "Odds ratio (T2D/normal)") +
  labs(x = "")
dev.off()
