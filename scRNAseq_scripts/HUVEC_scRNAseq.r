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

# Day0_1
control1.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27138_Control1/outs/filtered_feature_bc_matrix")
control1 <- CreateSeuratObject(raw.data = control1.data, min.cells = min_cells, min.genes = min_genes, project = "C.1")
control1@meta.data$sample <- "C.1"
control1@meta.data$time <- "Control"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = control1@data), value = TRUE)
percent.mito <- Matrix::colSums(control1@raw.data[mito.genes, ])/Matrix::colSums(control1@raw.data)
control1 <- AddMetaData(object = control1, metadata = percent.mito, col.name = "percent.mito")

# Day0_2
control2.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27139_Control2/outs/filtered_feature_bc_matrix")
control2 <- CreateSeuratObject(raw.data = control2.data, min.cells = min_cells, min.genes = min_genes, project = "C.2")
control2@meta.data$sample <- "C.2"
control2@meta.data$time <- "Control"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = control2@data), value = TRUE)
percent.mito <- Matrix::colSums(control2@raw.data[mito.genes, ])/Matrix::colSums(control2@raw.data)
control2 <- AddMetaData(object = control2, metadata = percent.mito, col.name = "percent.mito")

# Day3_1
T3d_1.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27140_3day1/outs/filtered_feature_bc_matrix")
T3d_1 <- CreateSeuratObject(raw.data = T3d_1.data, min.cells = min_cells, min.genes = min_genes, project = "T3d.1")
T3d_1@meta.data$sample <- "T3d.1"
T3d_1@meta.data$time <- "T3d"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T3d_1@data), value = TRUE)
percent.mito <- Matrix::colSums(T3d_1@raw.data[mito.genes, ])/Matrix::colSums(T3d_1@raw.data)
T3d_1 <- AddMetaData(object = T3d_1, metadata = percent.mito, col.name = "percent.mito")

# Day3_2
T3d_2.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27141_3day2/outs/filtered_feature_bc_matrix")
T3d_2 <- CreateSeuratObject(raw.data = T3d_2.data, min.cells = min_cells, min.genes = min_genes, project = "T3d.2")
T3d_2@meta.data$sample <- "T3d.2"
T3d_2@meta.data$time <- "T3d"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T3d_2@data), value = TRUE)
percent.mito <- Matrix::colSums(T3d_2@raw.data[mito.genes, ])/Matrix::colSums(T3d_2@raw.data)
T3d_2 <- AddMetaData(object = T3d_2, metadata = percent.mito, col.name = "percent.mito")

# Day7_1
T7d_1.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27142_7day1/outs/filtered_feature_bc_matrix")
T7d_1 <- CreateSeuratObject(raw.data = T7d_1.data, min.cells = min_cells, min.genes = min_genes, project = "T7d.1")
T7d_1@meta.data$sample <- "T7d.1"
T7d_1@meta.data$time <- "T7d"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T7d_1@data), value = TRUE)
percent.mito <- Matrix::colSums(T7d_1@raw.data[mito.genes, ])/Matrix::colSums(T7d_1@raw.data)
T7d_1 <- AddMetaData(object = T7d_1, metadata = percent.mito, col.name = "percent.mito")

# Day7_2
T7d_2.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27143_7day2/outs/filtered_feature_bc_matrix")
T7d_2 <- CreateSeuratObject(raw.data = T7d_2.data, min.cells = min_cells, min.genes = min_genes, project = "T7d.2")
T7d_2@meta.data$sample <- "T7d.2"
T7d_2@meta.data$time <- "T7d"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T7d_2@data), value = TRUE)
percent.mito <- Matrix::colSums(T7d_2@raw.data[mito.genes, ])/Matrix::colSums(T7d_2@raw.data)
T7d_2 <- AddMetaData(object = T7d_2, metadata = percent.mito, col.name = "percent.mito")

### Filter cells and normalize

thres_cells = 0.99

# Day 0
thresh.control1.mito = quantile(control1@meta.data$percent.mito,thres_cells)
thresh.control1.genes = round(quantile(control1@meta.data$nGene,thres_cells))
thresh.control2.mito = quantile(control2@meta.data$percent.mito,thres_cells)
thresh.control2.genes = round(quantile(control2@meta.data$nGene,thres_cells))

control1 <- FilterCells(object = control1, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(thresh.control1.genes, thresh.control1.mito))
control2 <- FilterCells(object = control2, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(thresh.control2.genes, thresh.control2.mito))

# Day 3
thresh.T3d_1.mito = quantile(T3d_1@meta.data$percent.mito,thres_cells)
thresh.T3d_1.genes = round(quantile(T3d_1@meta.data$nGene,thres_cells))
thresh.T3d_2.mito = quantile(T3d_2@meta.data$percent.mito,thres_cells)
thresh.T3d_2.genes = round(quantile(T3d_2@meta.data$nGene,thres_cells))

T3d_1 <- FilterCells(object = T3d_1, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(thresh.T3d_1.genes, thresh.T3d_1.mito))
T3d_2 <- FilterCells(object = T3d_2, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(thresh.T3d_2.genes, thresh.T3d_2.mito))

# Day 7
thresh.T7d_1.mito = quantile(T7d_1@meta.data$percent.mito,thres_cells)
thresh.T7d_1.genes = round(quantile(T7d_1@meta.data$nGene,thres_cells))
thresh.T7d_2.mito = quantile(T7d_2@meta.data$percent.mito,thres_cells)
thresh.T7d_2.genes = round(quantile(T7d_2@meta.data$nGene,thres_cells))

T7d_1 <- FilterCells(object = T7d_1, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(thresh.T7d_1.genes, thresh.T7d_1.mito))
T7d_2 <- FilterCells(object = T7d_2, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(thresh.T7d_2.genes, thresh.T7d_2.mito))

### Merge the object together and normalize
sc_full_data = MergeSeurat(control1,control2, add.cell.id1 = "C.1", add.cell.id2 = "C.2")
sc_full_data = MergeSeurat(sc_full_data,T3d_1, add.cell.id2 = "T3d.1")
sc_full_data = MergeSeurat(sc_full_data,T3d_2, add.cell.id2 = "T3d.2")
sc_full_data = MergeSeurat(sc_full_data,T7d_1, add.cell.id2 = "T7d.1")
sc_full_data = MergeSeurat(sc_full_data,T7d_2, add.cell.id2 = "T7d.2")

n_cells_after_filter = c()
for (i in names(table(as.factor(sc_full_data@meta.data$sample)))){
  n_cells_after_filter = c(n_cells_after_filter,sum(sc_full_data@cell.names %like% paste0("%",i,"%")))
}

### Detection of highly variable genes across the single cells
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/gene_dispersion_plot.png", width = 800, height = 500, units = "px")
sc_full_data <- FindVariableGenes(object = sc_full_data, mean.function = ExpMean, dispersion.function = LogVMR, 
                              x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.5, 
                              do.text = FALSE, do.contour = FALSE)
dev.off()

length(sc_full_data@var.genes) # number oh high variable genes detected according to the thresholds above
# Select the HVG as the top 1000 genes sorted by VMR
hv.genes = head(rownames(sc_full_data@hvg.info),1000)

### Scaling the data and removing unwanted sources of variation
# You can perform gene scaling on only the HVG, dramatically improving speed and memory use. 
# Since dimensional reduction is run only on HVG, this will not affect downstream results.
sc_full_data <- ScaleData(object = sc_full_data, genes.use = hv.genes, 
                          vars.to.regress = c("nUMI","percent.mito"), do.par = TRUE, num.cores = 4)

### Perform linear dimensional reduction
sc_full_data <- RunPCA(object = sc_full_data, pc.genes = hv.genes, pcs.compute = 50,
                       do.print = FALSE, pcs.print = 1:5, genes.print = 5)

# Saving PCA coordinates
PCA_coord = sc_full_data@dr$pca@cell.embeddings[,1:2]
cell_labels = rownames(PCA_coord)
cell_labels = gsub("C.1_","Day0.1_",cell_labels)
cell_labels = gsub("C.2_","Day0.2_",cell_labels)
cell_labels = gsub("T3d.1_","Day3.1_",cell_labels)
cell_labels = gsub("T3d.2_","Day3.1_",cell_labels)
cell_labels = gsub("T7d.1_","Day7.1_",cell_labels)
cell_labels = gsub("T7d.2_","Day7.2_",cell_labels)
rownames(PCA_coord) = cell_labels
write.table(PCA_coord,"/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/PCA_coordinates.txt", row.names = T, col.names = T, sep = '\t', quote = F)

# Elbow plot
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/pc_elbow_plot.png", width = 600, height = 500, units = "px")
PCElbowPlot(object = sc_full_data, num.pc = 50)
dev.off()

### Clustering
top_PCs = 20
sc_full_data <- FindClusters(object = sc_full_data, reduction.type = "pca", dims.use = 1:top_PCs, resolution = 0.6, print.output = 0, save.SNN = TRUE)

### t-SNE
sc_full_data <- RunTSNE(object = sc_full_data, dims.use = 1:top_PCs, reduction.use = "pca", do.fast = TRUE)

# Saving tSNE coordinates
tsne_coord = sc_full_data@dr$tsne@cell.embeddings
cell_labels = rownames(tsne_coord)
cell_labels = gsub("C.1_","Day0.1_",cell_labels)
cell_labels = gsub("C.2_","Day0.2_",cell_labels)
cell_labels = gsub("T3d.1_","Day3.1_",cell_labels)
cell_labels = gsub("T3d.2_","Day3.1_",cell_labels)
cell_labels = gsub("T7d.1_","Day7.1_",cell_labels)
cell_labels = gsub("T7d.2_","Day7.2_",cell_labels)
rownames(tsne_coord) = cell_labels
write.table(tsne_coord,"/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/tsne_coordinates.txt", row.names = T, col.names = T, sep = '\t', quote = F)

### Plot with PCA and t-SNE
pdf("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/Figure_1b-c.pdf", width = 16, height = 6)
p1 <-TSNEPlot(sc_full_data, do.return = T, pt.size = 0.5, group.by = "sample") + 
  theme(text = element_text(size=32),
        axis.text.x = element_text(size=32),
        axis.text.y = element_text(size=32)) +
  scale_x_continuous(breaks=c(-30,-15,0,15,30))
p2<-PCAPlot(object = sc_full_data, dim.1 = 1, dim.2 = 2, group.by = "time", do.return=TRUE) + 
  theme(text = element_text(size=32),
        axis.text.x = element_text(size=32),
        axis.text.y = element_text(size=32))
plot_grid(p1, p2)
dev.off()

# Plot specific genes
features_to_plot = c("ICAM1", "SERPINE1", "FN1", "COL4A2", "SMAD3", "CTGF", "NOS3", "ACTA2", "LINC00607", "LINC01013","LINC02154","LINC01235" )

# Save data to file
heatmap_intensities = as.matrix(t(sc_full_data@data[features_to_plot,]))
write.table(heatmap_intensities,"/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/gene_heatmap_intensities.txt", row.names = F, col.names = T, sep = '\t', quote = F)

p<-FeaturePlot(object = sc_full_data, features.plot = features_to_plot, nCol = 4,
               min.cutoff = "q10", max.cutoff = "q90", cols.use = c("white", "red"), pt.size = 0.5,
               dark.theme = TRUE, no.legend = FALSE, do.return = T)

temp = lapply(p, function(x){x + theme(legend.title = element_blank(),
                                       plot.title = element_text(size=24),
                                       axis.text = element_text(size=16), 
                                       axis.title = element_text(size=16),
                                       legend.text = element_text(size=16))})
p1 = temp[[1]]
p2 = temp[[2]]
p3 = temp[[3]]
p4 = temp[[4]]
p5 = temp[[5]]
p6 = temp[[6]]
p7 = temp[[7]] + ggtitle("NOS3 (eNOS)")
p8 = temp[[8]] +  ggtitle(expression(paste("ACTA2 (", alpha, "-SMA)")))
p9 = temp[[9]]
p10 = temp[[10]]
p11 = temp[[11]]
p12 = temp[[12]]

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/Figure_1e.png", width = 15, height = 10, units = "in", res = 300)
print(p)
dev.off()


### Differentially expressed genes

my.data=FetchData(sc_full_data,c("ident","PC1","nGene","orig.ident"))

control_cells = rownames(my.data[which(my.data$orig.ident == "C.1" | my.data$orig.ident == "C.2"),])
T3d_cells = rownames(my.data[which(my.data$orig.ident == "T3d.1" | my.data$orig.ident == "T3d.2"),])
T7d_cells = rownames(my.data[which(my.data$orig.ident == "T7d.1" | my.data$orig.ident == "T7d.2"),])

# Change ident to cells
sc_full_data_orig = SetAllIdent(sc_full_data, "orig.ident")
sc_full_data_orig = SetIdent(sc_full_data_orig, cells.use = control_cells, ident.use = "control")
sc_full_data_orig = SetIdent(sc_full_data_orig, cells.use = T3d_cells, ident.use = "T3d")
sc_full_data_orig = SetIdent(sc_full_data_orig, cells.use = T7d_cells, ident.use = "T7d")

levels(sc_full_data@ident)
levels(sc_full_data_orig@ident)

# Day 3 over Day 0
de_genes_T3_control <- FindMarkers(sc_full_data_orig, ident.1 = "T3d", ident.2 = "control", logfc.threshold = 0.25)
de_genes_T3_control = cbind(de_genes_T3_control,seq(1,nrow(de_genes_T3_control)))
colnames(de_genes_T3_control)[6] = "rank"
de_genes_T3_control[which(rownames(de_genes_T3_control)=="NOS3"),]

temp = cbind(rownames(de_genes_T3_control),de_genes_T3_control)
colnames(temp)[1] = 'gene'
write.table(temp,"/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/de_genes_T3_control.txt",row.names = F,col.names = T,quote = F, sep = '\t')

de_genes_T3_control = read.table("/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/de_genes_T3_control.txt", stringsAsFactors = F)
colnames(de_genes_T3_control) = de_genes_T3_control[1,]
de_genes_T3_control = de_genes_T3_control[-1,]
rownames(de_genes_T3_control) = de_genes_T3_control[,1]

temp=annotation[which(annotation$gene_name %in% de_genes_T3_control$gene),c("gene_name","gene_id")]
write.table(temp[,2],"/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/de_genes_T3_control_IDs.txt",row.names = F,col.names = F,quote = F)

# Day 7 over Day 0
de_genes_T7_control <- FindMarkers(sc_full_data_orig, ident.1 = "T7d", ident.2 = "control", logfc.threshold = 0.25)
de_genes_T7_control = cbind(de_genes_T7_control,seq(1,nrow(de_genes_T7_control)))
colnames(de_genes_T7_control)[6] = "rank"
de_genes_T7_control[which(rownames(de_genes_T7_control)=="NOS3"),]

temp = cbind(rownames(de_genes_T7_control),de_genes_T7_control)
colnames(temp)[1] = 'gene'
write.table(de_genes_T7_control,"/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/de_genes_T7_control.txt",row.names = F,col.names = T,quote = F, sep = '\t')

de_genes_T7_control = read.table("/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/de_genes_T7_control.txt", stringsAsFactors = F)
colnames(de_genes_T7_control) = de_genes_T7_control[1,]
de_genes_T7_control = de_genes_T7_control[-1,]
rownames(de_genes_T7_control) = de_genes_T7_control[,1]

temp=annotation[which(annotation$gene_name %in% de_genes_T7_control[,'gene']),c("gene_name","gene_id")]
write.table(temp[,2],"/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/de_genes_T7_control_IDs.txt",row.names = F,col.names = F,quote = F)

gene_list_to_david = merge(temp,de_genes_T7_control,by.x="gene_name",by.y="gene",all.x=T)
gene_list_to_david = gene_list_to_david[,-3]
colnames(gene_list_to_david)[c(3,4,5,6,7)] = c("avg_logFC_d7d0","pct_d7","pct_d0","p_val_adj_d7d0","rank_d7_d0")
gene_list_to_david = merge(gene_list_to_david,de_genes_T3_control,by.x="gene_name",by.y='gene',all.x=T)
gene_list_to_david = gene_list_to_david[,c(-8,-11)]
colnames(gene_list_to_david)[c(8,9,10,11)] = c("avg_logFC_d3d0","pct_d3","p_val_adj_d3d0","rank_d3d0")
gene_list_to_david$rank_d7_d0 = as.numeric(gene_list_to_david$rank_d7_d0)
gene_list_to_david = gene_list_to_david[order(gene_list_to_david$rank_d7_d0),]
write.table(gene_list_to_david,"/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/gene_list_to_david.txt",row.names = F,col.names = T,sep='\t',quote=F)

# Re scale all genes in order to plot the heatmap. It may be that a DE gene wasn't before a HVG: since scaled data are plotted in the
# heatmap, those genes couldn't be plotted if not scaled before.
sc_full_data_orig <- ScaleData(object = sc_full_data_orig, genes.use = rownames(de_genes_T7_control), 
                               vars.to.regress = c("nUMI","percent.mito"), do.par = TRUE, num.cores = 4)
# saveRDS(sc_full_data_orig, file="/dataOS/rcalandrelli/MARGI/RNAseq_02152019/sc_full_data_2_99th_orig.rds")

pdf("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/de_heatmap_T7_control.pdf", height=8, width=8)
DoHeatmap(sc_full_data_orig, genes.use = rownames(head(de_genes_T7_control,60)), slim.col.label = TRUE, remove.key = FALSE,
          col.low = "blue", col.mid = "black", col.high = "red", use.scaled = TRUE)
dev.off()

### Plot significative genes coming out from DAVID analysis

temp = read.xlsx("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/david_genes.xlsx", sheetName = "Heatmap V4")
david_genes = as.character(temp$gene_name)
david_genes = c(david_genes[1:53],"NOS3",david_genes[54:92])

sc_full_data_orig <- ScaleData(object = sc_full_data_orig, genes.use = david_genes, 
                               vars.to.regress = c("nUMI","percent.mito"), do.par = TRUE, num.cores = 4)

# Heatmap plotted manually by binning single cells
make_binned_sc_heatmap<-function(data, # data to be binned
                                 cell_binned, # how many single cells in each bin
                                 ordering, # TRUE to order single cell before binning
                                 which_gene # gene whose expression levels are used to order single cells
                                 ){
  # Control cells
  data_control = data[,control_cells]
  if (ordering == T){
    data_control = data_control[,order(data_control[which_gene,])]
  }
  data_control_binned = matrix(nrow=nrow(data_control),ncol=as.integer(ncol(data_control)/cell_binned))
  k = 0
  i = 1
  while (k < ncol(data_control)){
    ind_start = k + 1
    ind_end = k + cell_binned
    if (ind_end < ncol(data_control)){
      data_control_binned[,i] = rowMeans(data_control[,ind_start:ind_end])
    }
    k = k + cell_binned
    i = i + 1
  }
  rownames(data_control_binned) = rownames(data_control)
  colnames(data_control_binned) = paste0("Day0.bin.",seq(1,ncol(data_control_binned)))
  
  # Day 3
  data_T3d = data[,T3d_cells]
  if (ordering == T){
    data_T3d = data_T3d[,order(data_T3d[which_gene,])]
  }
  data_T3d_binned = matrix(nrow=nrow(data_T3d),ncol=as.integer(ncol(data_T3d)/cell_binned))
  k = 0
  i = 1
  while (k < ncol(data_T3d)){
    ind_start = k + 1
    ind_end = k + cell_binned
    if (ind_end < ncol(data_T3d)){
      data_T3d_binned[,i] = rowMeans(data_T3d[,ind_start:ind_end])
    }
    k = k + cell_binned
    i = i + 1
  }
  rownames(data_T3d_binned) = rownames(data_T3d)
  colnames(data_T3d_binned) = paste0("Day3.bin.",seq(1,ncol(data_T3d_binned)))
  
  # Day 7
  data_T7d = data[,T7d_cells]
  if (ordering == T){
    data_T7d = data_T7d[,order(data_T7d[which_gene,])]
  }
  data_T7d_binned = matrix(nrow=nrow(data_T7d),ncol=as.integer(ncol(data_T7d)/cell_binned))
  k = 0
  i = 1
  while (k < ncol(data_T7d)){
    ind_start = k + 1
    ind_end = k + cell_binned
    if (ind_end < ncol(data_T7d)){
      data_T7d_binned[,i] = rowMeans(data_T7d[,ind_start:ind_end])
    }
    k = k + cell_binned
    i = i + 1
  }
  rownames(data_T7d_binned) = rownames(data_T7d)
  colnames(data_T7d_binned) = paste0("Day7.bin.",seq(1,ncol(data_T7d_binned)))
  
  ### Merge the data
  data_merged = cbind(data_control_binned, NA, NA, data_T3d_binned, NA, NA, data_T7d_binned)
  return(data_merged)
}

temp = make_binned_sc_heatmap(sc_full_data_orig@scale.data[david_genes,], 100, TRUE, "SERPINE1")

temp[which(temp < -2.5, arr.ind = T)] = -2.5
temp[which(temp > 2.5, arr.ind = T)] = 2.5

write.table(temp,"/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/heatmap_binned_100.txt", row.names = T, col.names = T, sep="\t", quote = F)

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/Figure_1d.png", height=15, width=10, units = "in", res=300)
pheatmap(temp,
         color = colorRampPalette(c("blue", "yellow"))(100),
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         show_colnames = F,
         na_col = "white")
dev.off()


#### Intersection between Day3 and Day 7
temp = merge(de_genes_T3_control,de_genes_T7_control,by='gene')
temp[,"avg_logFC.x"] = as.numeric(temp[,"avg_logFC.x"])
temp[,"avg_logFC.y"] = as.numeric(temp[,"avg_logFC.y"])

temp_induced = temp[which(temp$avg_logFC.x>0 & temp$avg_logFC.y>0),]
temp_induced = cbind(temp_induced, rowMeans(temp_induced[,c("avg_logFC.x", "avg_logFC.y")]))
colnames(temp_induced)[14] = "mean_avg_logFC"
temp_induced = temp_induced[order(temp_induced$mean_avg_logFC, decreasing = T),]

temp_suppressed = temp[which(temp$avg_logFC.x<0 & temp$avg_logFC.y<0),]
temp_suppressed = cbind(temp_suppressed, rowMeans(temp_suppressed[,c("avg_logFC.x", "avg_logFC.y")]))
colnames(temp_suppressed)[14] = "mean_avg_logFC"
temp_suppressed = temp_suppressed[order(temp_suppressed$mean_avg_logFC),]

# Venn diagram
library(VennDiagram)

temp_T3 = rownames(de_genes_T3_control)
temp_T7 = rownames(de_genes_T7_control)

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/Supplementary_Figure_1up.png", height=800, width=1000, units = "px", res = 300)
draw.pairwise.venn(525, 211, 181,
                   #category = c("Day7-Day0 (525)","Day3-Day0 (211)"),
                   scaled = TRUE,
                   fill = c("#6e9bf8","#53b74c"),
                   cex = 0,
                   lwd = c(0.8,0.8),
                   fontfamily = "arial",
                   cat.pos = c(230,150),
                   cat.cex = c(1.5,1.5),
                   cat.fontfamily = "arial",
                   cat.fontface = "bold",
                   cat.dist = c(0.06, 0.05),
                   margin = 0.15)
dev.off()

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/Supplementary_Figure_1down.png", height=800, width=1000, units = "px", res = 300)
draw.pairwise.venn(301, 246, 173,
                   #category = c("Day3-Day0 (246)", "Day7-Day0 (301)"),
                   scaled = TRUE,
                   fill = c("#6e9bf8","#53b74c"),
                   cex = 0,
                   lwd = c(0.8,0.8),
                   fontfamily = "arial",
                   cat.pos = c(230,150),
                   cat.cex = c(1.5,1.5),
                   cat.fontfamily = "arial",
                   cat.fontface = "bold",
                   cat.dist = c(0.06, 0.05),
                   margin = 0.15)
dev.off()

### Parse data from DAVID software
temp=read.table("/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/DAVID/GAD_DISEASE.txt", sep='\t', header = T, stringsAsFactors = F)
temp = temp[,c(1,2,3,4,5,12)]

strsplit(temp[15,2],"\t")

temp = temp[order(temp[,"Count"], decreasing = T),]


term_selected=c("GO:0007155~cell adhesion",
                "GO:0001525~angiogenesis",
                "GO:0006954~inflammatory response",
                "GO:0006955~immune response",
                "GO:0098609~cell-cell adhesion",
                "GO:0001666~response to hypoxia",
                "GO:0030198~extracellular matrix organization",
                "GO:0016477~cell migration",
                "GO:0050900~leukocyte migration",
                "GO:0007568~aging",
                "GO:0007179~transforming growth factor beta receptor signaling pathway",
                "GO:0033209~tumor necrosis factor-mediated signaling pathway",
                "GO:0006979~response to oxidative stress",
                "GO:0050776~regulation of immune response",
                "GO:0008286~insulin receptor signaling pathway",
                "GO:0071260~cellular response to mechanical stimulus",
                "GO:0034097~response to cytokine",
                "GO:0071356~cellular response to tumor necrosis factor",
                "GO:0007219~Notch signaling pathway",
                "GO:0009749~response to glucose",
                "GO:0050727~regulation of inflammatory response",
                "GO:0050729~positive regulation of inflammatory response",
                "GO:0032869~cellular response to insulin stimulus",
                "GO:0045429~positive regulation of nitric oxide biosynthetic process",
                "GO:0071560~cellular response to transforming growth factor beta stimulus",
                "GO:0007249~I-kappaB kinase/NF-kappaB signaling",
                "GO:0038061~NIK/NF-kappaB signaling",
                "GO:0006006~glucose metabolic process",
                "GO:0070098~chemokine-mediated signaling pathway",
                "GO:0090023~positive regulation of neutrophil chemotaxis",
                "GO:0001974~blood vessel remodeling",
                "GO:0051276~chromosome organization")

GOTERM_BP_DIRECT=read.table("/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/DAVID/GOTERM_BP_DIRECT.txt", sep='\t', header = T, stringsAsFactors = F)

GOTERM_BP_DIRECT_selected = GOTERM_BP_DIRECT[which(GOTERM_BP_DIRECT$Term %in% term_selected),]

out_selected = matrix(nrow=0,ncol=2)
for (i in 1:nrow(GOTERM_BP_DIRECT_selected)){
  temp_genes = gsub(" ","",strsplit(GOTERM_BP_DIRECT_selected[i,"Genes"],",")[[1]])
  out_selected = rbind(out_selected,cbind(GOTERM_BP_DIRECT_selected[i,"Term"],temp_genes))
}
colnames(out_selected) = c("Term","gene_id")

temp = merge(out_selected, annotation, by="gene_id", all.x=TRUE)
temp = temp[,c(2,1,7,8)]
temp = temp[order(temp$Term),]
write.table(temp,"/dataOS/rcalandrelli/MARGI/RNAseq_02152019/result/DAVID/GOTERM_BP_DIRECT_selected.txt",row.names = F,col.names = T,sep='\t',quote = F)

### Odds ratio analysis

super_enhancers_sort_annotated_new = read.table("/dataOS/rcalandrelli/MARGI/super_enhancers_sort_annotated_new.txt", stringsAsFactors = F)
colnames(super_enhancers_sort_annotated_new) = super_enhancers_sort_annotated_new[1,]
super_enhancers_sort_annotated_new = super_enhancers_sort_annotated_new[-1,]

total_control = 65879147
total_T7d = 50031654

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

logfc_categories_and_odds_ratio <- function(x){
  my_gene = x[1]
  temp = sc_full_data@data[c("LINC00607", my_gene),control_cells]
  temp_0_0_control = sum(colSums(temp)==0) # 0,0
  temp_p_p_control = sum(temp[1,]*temp[2,] != 0) # +,+
  temp_p_0_control = sum(xor(temp[1,],temp[2,]) & temp[1,]>0) # +,0
  temp_0_p_control = sum(xor(temp[1,],temp[2,]) & temp[2,]>0) # 0,+
  mat_control = (matrix(c(temp_0_p_control,temp_0_0_control,temp_p_p_control,temp_p_0_control),2,2)+1)/length(control_cells)
  
  temp = sc_full_data@data[c("LINC00607", my_gene),T7d_cells]
  temp_0_0_T7d = sum(colSums(temp)==0) # 0,0
  temp_p_p_T7d = sum(temp[1,]*temp[2,] != 0) # +,+
  temp_p_0_T7d = sum(xor(temp[1,],temp[2,]) & temp[1,]>0) # +,0
  temp_0_p_T7d = sum(xor(temp[1,],temp[2,]) & temp[2,]>0) # 0,+
  mat_T7d = (matrix(c(temp_0_p_T7d,temp_0_0_T7d,temp_p_p_T7d,temp_p_0_T7d),2,2)+1)/length(T7d_cells)
  
  mat = log(mat_T7d/mat_control)
  df.mat = expand.grid(x = c("0", "+"), y = c("0", "+"))
  df.mat$logfc = c(mat[2,1],mat[2,2],mat[1,1],mat[1,2])
  df.mat$cells_control = c(temp_0_0_control,temp_p_0_control,temp_0_p_control,temp_p_p_control)
  df.mat$cells_T7d = c(temp_0_0_T7d,temp_p_0_T7d,temp_0_p_T7d,temp_p_p_T7d)
  odds_ratio = (mat_T7d[1,2]*mat_control[2,2])/(mat_control[1,2]*mat_T7d[2,2])
  return(c(t(df.mat[,3:5])[1,],t(df.mat[,3:5])[2,],t(df.mat[,3:5])[3,],odds_ratio))
}

temp = data.frame(t(apply(genes_ranked_sc,1,logfc_categories_and_odds_ratio)))
colnames(temp) = c("logfc_0_0","logfc_p_0","logfc_0_p","logfc_p_p",
                   "cells_control_0_0","cells_control_p_0","cells_control_0_p","cells_control_p_p",
                   "cells_T7d_0_0","cells_T7d_p_0","cells_T7d_0_p","cells_T7d_p_p",
                   "odds_ratio")
genes_ranked_sc_logfcCateg = cbind(genes_ranked_sc,temp)

# Remove one duplicate line of MALAT1
temp_ind = which(genes_ranked_sc_logfcCateg$gene=="MALAT1")[1]
genes_ranked_sc_logfcCateg = genes_ranked_sc_logfcCateg[-temp_ind,]

### Add information of the gene super enhancer
genes_ranked_sc_logfcCateg = merge(genes_ranked_sc_logfcCateg,super_enhancers_sort_annotated_new,by="SE_index_new")

### Plot only odds ratios in a violin plot

# "indexes_T7d_hubs" in HUVEC_iMARGI.r"

hubs_interacting_with_linc = unique(c(as.numeric(indexes_T7d_hubs[which(indexes_T7d_hubs[,1]==145),2]),
                                          as.numeric(indexes_T7d_hubs[which(indexes_T7d_hubs[,2]==145),1])))

temp = genes_ranked_sc_logfcCateg[which(genes_ranked_sc_logfcCateg$SE_index_new %in% hubs_interacting_with_linc),]
temp = merge(temp, annotation[,c("gene_name","gene_biotype")], by.x = "gene", by.y = "gene_name", all.x = T)
temp = temp[which(temp$gene_biotype %in% c("protein_coding","lincRNA")),]
temp = temp[order(temp$odds_ratio, decreasing = T),]

temp$gene_label = "Other gene"
temp[which(temp$gene == "SERPINE1"),"gene_label"] = "SERPINE1"
temp[which(temp$gene == "RUNX1"),"gene_label"] = "RUNX1"
temp[which(temp$gene == "TRIO"),"gene_label"] = "TRIO"

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/Figure_5b.png", height=1000, width=1000, units = "px", res = 200)
ggplot(temp, aes(x=0, y=odds_ratio)) + 
  geom_violin() +
  geom_point(position = position_jitter(0.3), aes(shape = gene_label), size = 2) +
  scale_shape_manual(values = c(4,2,5,0)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text=element_text(size=14),
        legend.title = element_blank(),
        axis.ticks.length = unit(0.2, "cm")) +
  labs(y = "Odds ratio (Day 7/Day 0)") +
  labs(x = "")
dev.off()

