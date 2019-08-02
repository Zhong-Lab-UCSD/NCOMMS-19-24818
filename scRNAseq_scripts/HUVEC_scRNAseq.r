################ scRNA-seq analysis with Seurat
library(ggbio)

summary = read.csv("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/summary.csv", stringsAsFactors = F)

median(apply(T7d_1@data,2,function(x) sum(x>0))) # to see the median number of genes per cell

number.cells = c(11458, 15511,	4358,	13413,	8095,	6770)
median.genes = c(1226, 1018, 2332, 1196, 1706, 1928)

sample=rep(c("C.1", "C.2", "T3d.1", "T3d.2", "T7d.1","T7d.2"), 2)
feature=c(rep("Number of cells" , 6), rep("Median genes per cell" , 6))
value = c(number.cells, median.genes)
data=data.frame(sample,feature,value)

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/FigureS1.png", width = 14, height = 10, units = "in", res=200)
ggplot(data, aes(fill=feature, y=value, x=sample)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Sample", y = "") + 
  theme(legend.title = element_blank(),
        text = element_text(size=30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30))
dev.off()


### Load the datasets
library(Seurat) # version 2.3.4
library(dplyr)
library(cowplot)

min_cells = 0
min_genes = 0

# Control 1
control1.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27138_Control1/outs/filtered_feature_bc_matrix")
control1 <- CreateSeuratObject(raw.data = control1.data, min.cells = min_cells, min.genes = min_genes, project = "C.1")
control1@meta.data$sample <- "C.1"
control1@meta.data$time <- "Control"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = control1@data), value = TRUE)
percent.mito <- Matrix::colSums(control1@raw.data[mito.genes, ])/Matrix::colSums(control1@raw.data)
control1 <- AddMetaData(object = control1, metadata = percent.mito, col.name = "percent.mito")

# Control 2
control2.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27139_Control2/outs/filtered_feature_bc_matrix")
control2 <- CreateSeuratObject(raw.data = control2.data, min.cells = min_cells, min.genes = min_genes, project = "C.2")
control2@meta.data$sample <- "C.2"
control2@meta.data$time <- "Control"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = control2@data), value = TRUE)
percent.mito <- Matrix::colSums(control2@raw.data[mito.genes, ])/Matrix::colSums(control2@raw.data)
control2 <- AddMetaData(object = control2, metadata = percent.mito, col.name = "percent.mito")

# T3d_1
T3d_1.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27140_3day1/outs/filtered_feature_bc_matrix")
T3d_1 <- CreateSeuratObject(raw.data = T3d_1.data, min.cells = min_cells, min.genes = min_genes, project = "T3d.1")
T3d_1@meta.data$sample <- "T3d.1"
T3d_1@meta.data$time <- "T3d"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T3d_1@data), value = TRUE)
percent.mito <- Matrix::colSums(T3d_1@raw.data[mito.genes, ])/Matrix::colSums(T3d_1@raw.data)
T3d_1 <- AddMetaData(object = T3d_1, metadata = percent.mito, col.name = "percent.mito")

# T3d_2
T3d_2.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27141_3day2/outs/filtered_feature_bc_matrix")
T3d_2 <- CreateSeuratObject(raw.data = T3d_2.data, min.cells = min_cells, min.genes = min_genes, project = "T3d.2")
T3d_2@meta.data$sample <- "T3d.2"
T3d_2@meta.data$time <- "T3d"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T3d_2@data), value = TRUE)
percent.mito <- Matrix::colSums(T3d_2@raw.data[mito.genes, ])/Matrix::colSums(T3d_2@raw.data)
T3d_2 <- AddMetaData(object = T3d_2, metadata = percent.mito, col.name = "percent.mito")

# T7d_1
T7d_1.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27142_7day1/outs/filtered_feature_bc_matrix")
T7d_1 <- CreateSeuratObject(raw.data = T7d_1.data, min.cells = min_cells, min.genes = min_genes, project = "T7d.1")
T7d_1@meta.data$sample <- "T7d.1"
T7d_1@meta.data$time <- "T7d"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T7d_1@data), value = TRUE)
percent.mito <- Matrix::colSums(T7d_1@raw.data[mito.genes, ])/Matrix::colSums(T7d_1@raw.data)
T7d_1 <- AddMetaData(object = T7d_1, metadata = percent.mito, col.name = "percent.mito")

# T7d_2
T7d_2.data <- Read10X(data.dir = "/dataOS/rcalandrelli/MARGI/RNAseq_02152019/counts/count_27143_7day2/outs/filtered_feature_bc_matrix")
T7d_2 <- CreateSeuratObject(raw.data = T7d_2.data, min.cells = min_cells, min.genes = min_genes, project = "T7d.2")
T7d_2@meta.data$sample <- "T7d.2"
T7d_2@meta.data$time <- "T7d"

mito.genes <- grep(pattern = "^MT-", x = rownames(x = T7d_2@data), value = TRUE)
percent.mito <- Matrix::colSums(T7d_2@raw.data[mito.genes, ])/Matrix::colSums(T7d_2@raw.data)
T7d_2 <- AddMetaData(object = T7d_2, metadata = percent.mito, col.name = "percent.mito")

### Merge the object together just to plot samples together at this stage before any filtering
temp = MergeSeurat(control1,control2, add.cell.id1 = "C.1", add.cell.id2 = "C.2")
temp = MergeSeurat(temp,T3d_1, add.cell.id2 = "T3d.1")
temp = MergeSeurat(temp,T3d_2, add.cell.id2 = "T3d.2")
temp = MergeSeurat(temp,T7d_1, add.cell.id2 = "T7d.1")
temp = MergeSeurat(temp,T7d_2, add.cell.id2 = "T7d.2")

textsize = 30

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/violin_plot.png", width = 12, height = 10, units = "in", res = 200)
VlnPlot(object = temp, features.plot = c("nGene"), nCol = 3,
              size.x.use = 18, size.y.use = 18, size.title.use = 22, point.size.use = 0, group.by = "sample", do.return = TRUE) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
dev.off()

p2 <- VlnPlot(object = temp, features.plot = c("nUMI"), nCol = 3,
              size.x.use = 18, size.y.use = 18, size.title.use = 22, point.size.use = 0, group.by = "sample", do.return = TRUE) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))

p3 <- VlnPlot(object = temp, features.plot = c("percent.mito"), nCol = 3,
              size.x.use = 18, size.y.use = 18, size.title.use = 22, point.size.use = 0, group.by = "sample", do.return = TRUE) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
plot_grid(p1,p2,p3,nrow = 1,ncol = 3)
dev.off()

### Scatter plots

thres_cells = 0.99

# Control
thresh.control1.mito = quantile(control1@meta.data$percent.mito,thres_cells)
thresh.control1.genes = round(quantile(control1@meta.data$nGene,thres_cells))
thresh.control2.mito = quantile(control2@meta.data$percent.mito,thres_cells)
thresh.control2.genes = round(quantile(control2@meta.data$nGene,thres_cells))

# T3d
thresh.T3d_1.mito = quantile(T3d_1@meta.data$percent.mito,thres_cells)
thresh.T3d_1.genes = round(quantile(T3d_1@meta.data$nGene,thres_cells))
thresh.T3d_2.mito = quantile(T3d_2@meta.data$percent.mito,thres_cells)
thresh.T3d_2.genes = round(quantile(T3d_2@meta.data$nGene,thres_cells))

# T7d
thresh.T7d_1.mito = quantile(T7d_1@meta.data$percent.mito,thres_cells)
thresh.T7d_1.genes = round(quantile(T7d_1@meta.data$nGene,thres_cells))
thresh.T7d_2.mito = quantile(T7d_2@meta.data$percent.mito,thres_cells)
thresh.T7d_2.genes = round(quantile(T7d_2@meta.data$nGene,thres_cells))

textsize = 18
labelsize = 18

# Control
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/scatter_plot_control.png", width = 18, height = 5, units = "in", res = 200)
p1<-ggplot(control1@meta.data, aes(x=nUMI, y=percent.mito)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.control1.mito, color = "red") + ggtitle("Control.1") +
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p2<-ggplot(control1@meta.data, aes(x=nUMI, y=nGene)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.control1.genes, color = "red") + ggtitle("Control.1") +
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p3<-ggplot(control2@meta.data, aes(x=nUMI, y=percent.mito)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.control2.mito, color = "red") + ggtitle("Control.2") +
  theme(text = element_text(size=labelsize), 
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p4<-ggplot(control2@meta.data, aes(x=nUMI, y=nGene)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.control2.genes, color = "red") + ggtitle("Control.2") +
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
plot_grid(p1,p2,p3,p4,nrow = 1,ncol = 4)
dev.off()

# T3d
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/scatter_plot_T3d.png", width = 18, height = 5, units = "in", res = 200)
p1<-ggplot(T3d_1@meta.data, aes(x=nUMI, y=percent.mito)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.T3d_1.mito, color = "red") + ggtitle("T3d.1") +
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=textsize))
p2<-ggplot(T3d_1@meta.data, aes(x=nUMI, y=nGene)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.T3d_1.genes, color = "red") + ggtitle("T3d.1") +
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=textsize))
p3<-ggplot(T3d_2@meta.data, aes(x=nUMI, y=percent.mito)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.T3d_2.mito, color = "red") + ggtitle("T3d.2") +
  theme(text = element_text(size=labelsize), 
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p4<-ggplot(T3d_2@meta.data, aes(x=nUMI, y=nGene)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.T3d_2.genes, color = "red") + ggtitle("T3d.2") +
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
plot_grid(p1,p2,p3,p4,nrow = 1,ncol = 4)
dev.off()

# T7d
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/scatter_plot_T7d.png", width = 18, height = 5, units = "in", res = 200)
p1<-ggplot(T7d_1@meta.data, aes(x=nUMI, y=percent.mito)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.T7d_1.mito, color = "red") + ggtitle("T7d.1") +
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p2<-ggplot(T7d_1@meta.data, aes(x=nUMI, y=nGene)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.T7d_1.genes, color = "red") + ggtitle("T7d.1") +
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p3<-ggplot(T7d_2@meta.data, aes(x=nUMI, y=percent.mito)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.T7d_2.mito, color = "red") + ggtitle("T7d.2") +
  theme(text = element_text(size=labelsize), 
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p4<-ggplot(T7d_2@meta.data, aes(x=nUMI, y=nGene)) + geom_point(size=0.5) + geom_hline(yintercept=thresh.T7d_2.genes, color = "red") + ggtitle("T7d.2") +
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
plot_grid(p1,p2,p3,p4,nrow = 1,ncol = 4)
dev.off()


### Filter cells and normalize
control1 <- FilterCells(object = control1, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(thresh.control1.genes, thresh.control1.mito))
control2 <- FilterCells(object = control2, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(thresh.control2.genes, thresh.control2.mito))

T3d_1 <- FilterCells(object = T3d_1, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(thresh.T3d_1.genes, thresh.T3d_1.mito))
T3d_2 <- FilterCells(object = T3d_2, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(thresh.T3d_2.genes, thresh.T3d_2.mito))

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

textsize = 20
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/violin_plot_after_filter.png", width = 1500, height = 800, units = "px")
p1 <- VlnPlot(object = sc_full_data, features.plot = c("nGene"), nCol = 3,
              size.x.use = 18, size.y.use = 18, size.title.use = 22, point.size.use = 0, group.by = "sample", do.return = TRUE) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p2 <- VlnPlot(object = sc_full_data, features.plot = c("nUMI"), nCol = 3,
              size.x.use = 18, size.y.use = 18, size.title.use = 22, point.size.use = 0, group.by = "sample", do.return = TRUE) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p3 <- VlnPlot(object = sc_full_data, features.plot = c("percent.mito"), nCol = 3,
              size.x.use = 18, size.y.use = 18, size.title.use = 22, point.size.use = 0, group.by = "sample", do.return = TRUE) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
plot_grid(p1,p2,p3,nrow = 1,ncol = 3)
dev.off()

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

# Elbow plot
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/pc_elbow_plot.png", width = 600, height = 500, units = "px")
PCElbowPlot(object = sc_full_data, num.pc = 50)
dev.off()


### Clustering
top_PCs = 19
sc_full_data <- FindClusters(object = sc_full_data, reduction.type = "pca", dims.use = 1:top_PCs, resolution = 0.6, print.output = 0, save.SNN = TRUE)

# Plot clusters
textsize = 18
labelsize = 20

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/pc_scatter_plot_bysample.png", width = 1000, height = 400, units = "px")
p1<-PCAPlot(object = sc_full_data, dim.1 = 1, dim.2 = 2, do.return=TRUE) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p2<-PCAPlot(object = sc_full_data, dim.1 = 1, dim.2 = 2, group.by = "sample", do.return=TRUE) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
plot_grid(p1, p2)
dev.off()

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/pc_scatter_plot_bytime.png", width = 1000, height = 400, units = "px")
p1<-PCAPlot(object = sc_full_data, dim.1 = 1, dim.2 = 2, do.return=TRUE) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p2<-PCAPlot(object = sc_full_data, dim.1 = 1, dim.2 = 2, group.by = "time", do.return=TRUE) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
plot_grid(p1, p2)
dev.off()

### t-SNE
sc_full_data <- RunTSNE(object = sc_full_data, dims.use = 1:top_PCs, reduction.use = "pca", do.fast = TRUE)
# saveRDS(sc_full_data, file="/dataOS/rcalandrelli/MARGI/RNAseq_02152019/sc_full_data_2_99th.rds")

### USE THIS!
# sc_full_data = readRDS("/dataOS/rcalandrelli/MARGI/RNAseq_02152019/sc_full_data_2_99th.rds", refhook = NULL)

# Plot clusters
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/tsne_scatter_plot_bysample.png", width = 1000, height = 400, units = "px")
p1 <- TSNEPlot(sc_full_data, do.return = T, pt.size = 0.5, do.label = T, label.size = 8) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p2 <- TSNEPlot(sc_full_data, do.return = T, pt.size = 0.5, group.by = "sample") + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
plot_grid(p1, p2)
dev.off()

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/tsne_scatter_plot_bytime.png", width = 1000, height = 400, units = "px")
p1 <- TSNEPlot(sc_full_data, do.return = T, pt.size = 0.5, do.label = T, label.size = 8) + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
p2 <- TSNEPlot(sc_full_data, do.return = T, pt.size = 0.5, group.by = "time") + 
  theme(text = element_text(size=labelsize),
        axis.text.x = element_text(size=textsize),
        axis.text.y = element_text(size=textsize))
plot_grid(p1, p2)
dev.off()

### Plot paper with PCA and t-SNE
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/PCA-tsne_scatter_plot_bytime.png", width = 13, height = 5, units = "in", res = 300)
p1 <- TSNEPlot(sc_full_data, do.return = T, pt.size = 0.5, group.by = "time") + 
  theme(text = element_text(size=22),
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22)) +
  scale_fill_manual(name="", values=c("brown1","darkolivegreen4","burlywood3", labels=c("condition1", "condition2", "condition3")))
p2<-PCAPlot(object = sc_full_data, dim.1 = 1, dim.2 = 2, group.by = "time", do.return=TRUE) + 
  theme(text = element_text(size=22),
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22))
plot_grid(p2, p1)
dev.off()

# Plot specific genes
features_to_plot = c("ICAM1", "SERPINE1", "FN1", "COL4A2", "SMAD3", "CTGF", "NOS3", "ACTA2", "LINC00607", "LINC01013","LINC02154","LINC01235" )
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/tsne_genes_plot2.png", width = 15, height = 10, units = "in", res = 300)
FeaturePlot(object = sc_full_data, features.plot = features_to_plot, nCol = 4,
            min.cutoff = "q10", max.cutoff = "q90", cols.use = c("white", "red"), pt.size = 0.5,
            dark.theme = TRUE, no.legend = FALSE)
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

# T3 over control
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

# Re scale all genes in order to plot the heatmap. It may be that a DE gene wasn't before a HVG: since scaled data are plotted in the
# heatmap, those genes couldn't be plotted if not scaled before.
sc_full_data_orig <- ScaleData(object = sc_full_data_orig, genes.use = rownames(de_genes_T3_control), 
                               vars.to.regress = c("nUMI","percent.mito"), do.par = TRUE, num.cores = 4)
# saveRDS(sc_full_data_orig, file="/dataOS/rcalandrelli/MARGI/RNAseq_02152019/sc_full_data_2_99th_orig.rds")

pdf("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/de_heatmap_T3_control.pdf", height=8, width=8)
DoHeatmap(sc_full_data_orig, genes.use = rownames(head(de_genes_T3_control,60)), slim.col.label = TRUE, remove.key = FALSE,
          col.low = "blue", col.mid = "black", col.high = "red", use.scaled = TRUE)
dev.off()

### Top 20 DE LINC RNAs

temp = de_genes_T3_control[grepl("LINC", rownames(de_genes_T3_control)),]
temp = cbind(temp,matrix(0,nrow(temp),3))
for (i in 1:nrow(temp)){
  temp_linc = rownames(temp[i,])
  avg_express = c(mean(sc_full_data@data[temp_linc,control_cells]),mean(sc_full_data@data[temp_linc,T3d_cells]),mean(sc_full_data@data[temp_linc,T7d_cells]))
  temp[i,c(7,8,9)] = avg_express
}
colnames(temp)[c(7,8,9)] = c("avg_expr_control","avg_expr_T3d","avg_expr_T7d")
temp = cbind(rownames(temp),temp)
colnames(temp)[1] = 'gene'
write.table(temp,"/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/DE_LINC_T3d.txt",row.names = F,col.names = T,sep='\t',quote=F)

# T7 over control
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

# Plot significative genes coming out from DAVID analysis

temp = read.xlsx("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/david_genes.xlsx", sheetName = "Heatmap V4")
david_genes = as.character(temp$gene_name)

sc_full_data_orig <- ScaleData(object = sc_full_data_orig, genes.use = david_genes, 
                               vars.to.regress = c("nUMI","percent.mito"), do.par = TRUE, num.cores = 4)
my_colors = c(rep("#ff0000", 27), rep("#75015f",3), rep("#0208ff",17), rep("#05b703",6), rep("#ff7300",11), rep("#ff00ce",10), rep("#b7b701",14), rep("#97a8a7",4))
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/david_genes.png", height=15, width=10, units = "in", res=300)
p<-DoHeatmap(sc_full_data_orig, genes.use = david_genes, slim.col.label = TRUE, remove.key = FALSE,
          col.low = "blue", col.mid = "black", col.high = "red", use.scaled = TRUE)
print(p)
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

pdf("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/de_heatmap_common_induced.pdf", height=8, width=8)
DoHeatmap(sc_full_data_orig, genes.use = head(temp_induced[,1],60), slim.col.label = TRUE, remove.key = FALSE,
          col.low = "blue", col.mid = "black", col.high = "red", use.scaled = TRUE)
dev.off()

pdf("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/de_heatmap_common_suppressed.pdf", height=8, width=8)
DoHeatmap(sc_full_data_orig, genes.use = head(temp_suppressed[,1],60), slim.col.label = TRUE, remove.key = FALSE,
          col.low = "blue", col.mid = "black", col.high = "red", use.scaled = TRUE)
dev.off()

# Venn diagram
library(VennDiagram)

temp_T3 = rownames(de_genes_T3_control)
temp_T7 = rownames(de_genes_T7_control)

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/de_genes_venn.png", height=800, width=800, units = "px", res = 100)
draw.pairwise.venn(457, 826, 356,
                   #category = c("Day3-Day0 (457)", "Day7-Day0 (826)"),
                   scaled = TRUE,
                   fill = c("red","blue"),
                   cex = 2,
                   fontfamily = "arial",
                   cat.pos = c(230,150),
                   cat.cex = c(1.5,1.5),
                   cat.fontfamily = "arial",
                   cat.fontface = "bold",
                   cat.dist = c(0.06, 0.05),
                   margin = 0.08)
dev.off()

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/de_up_genes_venn.png", height=800, width=1000, units = "px", res = 300)
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

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/RNAseq_02152019/result/de_down_genes_venn.png", height=800, width=1000, units = "px", res = 300)
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






