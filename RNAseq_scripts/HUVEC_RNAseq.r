library("DESeq2")
library("ggplot2")
library(xlsx)
library(cowplot)

directory = "/dataOS/rcalandrelli/MARGI/LNA_rna_seq/"
annotation <- read.table('/dataOS/rcalandrelli/MARGI/Homo_sapiens.GRCh38.84.chr.gtf_to_geneTable.tsv', stringsAsFactors = F)
colnames(annotation) = annotation[1,]
annotation = annotation[-1,]

samples = c("Scr-NG-1","Scr-HT-1","LNA1-HT-1","LNA4-HT-1","Scr-NG-2","Scr-HT-2","LNA1-HT-2","LNA4-HT-2")
condition = c("Scr_NG","Scr_HT","LNA1_HT","LNA4_HT","Scr_NG","Scr_HT","LNA1_HT","LNA4_HT")

### featureCounts statistics
alignment_summary = read.xlsx(paste0(directory,"info_mapping.xlsx"), sheetIndex = 1)
alignment_summary = alignment_summary[,c(1,2,4,5)]
alignment_summary$assigned = 0
alignment_summary$percentage_assigned = 0

k = 0
for (i in 32877:32884){
  k = k + 1
  fc_summary = read.table(paste0(directory,"Sample_",i,"/featureCounts/counts.txt.summary"), skip = 1)
  alignment_summary[k,"assigned"] = fc_summary[1,2]
  alignment_summary[k,"percentage_assigned"] = round(fc_summary[1,2]/alignment_summary[k,"Uniquely.mapped.reads"], digits = 4)
}

colnames(alignment_summary) = c("Sample","Input reads","Uniquely mapped reads","% of uniquely mapped reads (UMR)","Assigned reads to features","% of assigned reads to features (vs UMR)")
write.table(alignment_summary,paste0(directory,"featureCounts_summary.txt"),row.names = F, col.names = T, sep = "\t", quote = F)

### Importing the data into DESeq2
coldata <- data.frame(condition)
rownames(coldata) <- samples

temp = read.table(paste0(directory,"Sample_32877/featureCounts/counts.txt"),sep="\t")
colnames(temp) = sapply(temp[1,],function(x){return(as.character(x))})
temp = temp[-1,]

cts = matrix(0, nrow = 60675, ncol = 8)
rownames(cts) = temp[,1]
colnames(cts) = samples

k = 0
for (i in c(32877:32884)){
  k = k + 1
  temp = read.table(paste0(directory,"Sample_",i,"/featureCounts/counts.txt"),sep="\t")
  colnames(temp) = sapply(temp[1,],function(x){return(as.character(x))})
  temp = temp[-1,]
  rownames(temp) = temp[,1]
  cts[rownames(temp),k] = as.numeric(as.character(temp[,7]))
}


#################### LNA1_HT as reference sample
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
# Updating rownames using gene names
temp = annotation[which(annotation[,5] %in% rownames(dds)),6]
rownames(dds) = temp

# Filtering genes with low read counts (less than 10 for all the samples)
hist(rowSums(counts(dds))[rowSums(counts(dds))<10])
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

### Normalizing using the "median of ratios method"
sizeFactors(estimateSizeFactors(dds)) # normalization factors
dds <- estimateSizeFactors(dds) # normalizing the data
counts(dds, normalized = T)["SERPINE1",]

### Differential expression analysis
dds$condition <- relevel(dds$condition, ref = "LNA1_HT") # selecting the reference condition
dds = DESeq(dds)

# 1 - Scr_NG vs LNA1_HT
# 2 - Scr_HT vs LNA1_HT
# 3 - LNA4_HT vs LNA1_HT

my_alpha = 0.05
res_1 <- results(dds, contrast=c("condition","Scr_NG","LNA1_HT"), alpha = my_alpha)
res_2 <- results(dds, contrast=c("condition","Scr_HT","LNA1_HT"), alpha = my_alpha)
res_3 <- results(dds, contrast=c("condition","LNA4_HT","LNA1_HT"), alpha = my_alpha)

# If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA. 
# These outlier counts are detected by Cook's distance.
# If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA.
res_1[which(is.na(res_1$pvalue)),]
res_2[which(is.na(res_2$pvalue)),]
res_3[which(is.na(res_3$pvalue)),]

# Summary of the results from DE analysis
summary(res_1) 
summary(res_2)
summary(res_3)

# Checking the DE genes based of a threshold of the FDR
sum(res_1$padj < my_alpha, na.rm=TRUE)
sum(res_2$padj < my_alpha, na.rm=TRUE)
sum(res_3$padj < my_alpha, na.rm=TRUE)

res_1_de = res_1[which(res_1$padj< my_alpha),]
res_2_de = res_2[which(res_2$padj< my_alpha),]
res_3_de = res_3[which(res_3$padj< my_alpha),]

# adding normalized count data and gene_id
temp_1 = merge(res_1_de, counts(dds, normalized = T), by="row.names")
temp_2 = merge(res_2_de, counts(dds, normalized = T), by="row.names")
temp_3 = merge(res_3_de, counts(dds, normalized = T), by="row.names")

colnames(temp_1)[1] = "gene_name"
colnames(temp_2)[1] = "gene_name"
colnames(temp_3)[1] = "gene_name"

temp_1 = merge(temp_1,annotation[,c("gene_name","gene_id")], by="gene_name", all.x=T)
temp_2 = merge(temp_2,annotation[,c("gene_name","gene_id")], by="gene_name", all.x=T)
temp_3 = merge(temp_3,annotation[,c("gene_name","gene_id")], by="gene_name", all.x=T)

# Up and downregulated genes
nrow(temp_1[which(temp_1$log2FoldChange>0),])
nrow(temp_2[which(temp_2$log2FoldChange>0),])
nrow(temp_3[which(temp_3$log2FoldChange>0),])

nrow(temp_1[which(temp_1$log2FoldChange<0),])
nrow(temp_2[which(temp_2$log2FoldChange<0),])
nrow(temp_3[which(temp_3$log2FoldChange<0),])

# Sorting results
res_1_Ordered <- temp_1[order(temp_1$pvalue),]
res_2_Ordered <- temp_2[order(temp_2$pvalue),]
res_3_Ordered <- temp_3[order(temp_3$pvalue),]

write.table(res_1_Ordered,paste0(directory,"result/Scr_NG_vs_LNA1_HT.txt"), sep='\t', row.names = T, col.names = T, quote = F)
write.table(res_2_Ordered,paste0(directory,"result/Scr_HT_vs_LNA1_HT.txt"), sep='\t', row.names = T, col.names = T, quote = F)
write.table(res_3_Ordered,paste0(directory,"result/LNA4_HT_vs_LNA1_HT.txt"), sep='\t', row.names = T, col.names = T, quote = F)

write.xlsx(res_1_Ordered,paste0(directory,"result/DE_analysis_vs_LNA1_HT.xlsx"), sheetName = "Scr_NG_vs_LNA1_HT", row.names = F)
write.xlsx(res_2_Ordered,paste0(directory,"result/DE_analysis_vs_LNA1_HT.xlsx"), sheetName = "Scr_HT_vs_LNA1_HT", row.names = F, append = T)
write.xlsx(res_3_Ordered,paste0(directory,"result/DE_analysis_vs_LNA1_HT.xlsx"), sheetName = "LNA4_HT_vs_LNA1_HT", row.names = F, append = T)


########## Data transformations and visualization

# Variance stabilizing transformations (VST)
vsd <- vst(dds, blind=FALSE)
assay(vsd)["SERPINE1",]

### Heatmap of the sample-to-sample distances
library("RColorBrewer")
library("pheatmap")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
png(paste0(directory,"result/distance_heatmap.png"), width = 10, height = 8, res = 200, units = "in")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         fontsize = 18)
dev.off()

### Principal component analysis
library(ggrepel)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
pcaData$nudge_custom = 1
pcaData[c(4,6),"nudge_custom"] = 0
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(paste0(directory,"result/PCA_plot.png"), width = 10, height = 10, res = 200, units = "in")
    ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=5) +
  geom_text(aes(label=name), nudge_x = ifelse(pcaData$nudge_custom==1,0.8,-0.8), nudge_y = ifelse(pcaData$nudge_custom==1,0.6,-0.6), size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() + 
  theme(axis.text=element_text(size=20), 
        axis.title=element_text(size=20),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

### Heatmap of selected genes DE between scr-HT and LNA1-HT and known pathways
library("pheatmap")

temp_table = c("SERPINE1",
"CXCL10",
"CXCL11",
"VCAM1",
"ICAM1",
"CCL2",
"LGALS9",
"CTSK",
"DHX58",
"EBI3",
"IRF7",
"BATF2",
"IL4I1",
"IFI27",
"CSF1",
"JAK3",
"RSAD2",
"MMP10",
"MMP2",
"FN1",
"COL8A1",
"COL5A1",
"COL5A2",
"COL4A2",
"COL4A1",
"ITGA2",
"ITGB1",
"SCX",
"TIMP2",
"THBS3",
"TGFBR1",
"SMAD3",
"TGFBR2",
"PDGFB",
"ACTA2")

# Check log2fc is positive for all the genes
temp=res_2_Ordered[which(res_2_Ordered$gene_name %in% temp_table),]

### Plot z-score
temp_avg = counts(dds, normalized = T)[as.character(temp_table[,1]),]
temp_avg[,"Scr-NG-1"] = rowMeans(temp_avg[,c("Scr-NG-1","Scr-NG-2")])
temp_avg[,"Scr-HT-1"] = rowMeans(temp_avg[,c("Scr-HT-1","Scr-HT-2")])
temp_avg[,"LNA1-HT-1"] = rowMeans(temp_avg[,c("LNA1-HT-1","LNA1-HT-2")])
temp_avg[,"LNA4-HT-1"] = rowMeans(temp_avg[,c("LNA4-HT-1","LNA4-HT-2")])
temp_avg = temp_avg[,c(-5:-8)]
colnames(temp_avg) = c("NM-scr","HT-scr","HT-LNA1","HT-LNA2")
temp_avg = log2(temp_avg + 1)

z_temp_avg = (temp_avg - rowMeans(temp_avg)) / apply(temp_avg,1,sd)

write.table(t(z_temp_avg),"/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/LNA_rna_seq/result/Scr_HT_vs_Scr_NG_gene_heatmap_avg_duplicates_Z_score.txt", row.names = T, col.names = T, sep = "\t", quote = F)

pdf("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/LNA_rna_seq/result/Figure_4d.pdf", height = 9, width = 3)
pheatmap(z_temp_avg[,c(4,3,2,1)],
         color = colorRampPalette(c("blue", "yellow"))(100),
         cluster_rows=FALSE,
         cluster_cols=FALSE)
dev.off()

### Plot z-score Table S5
genes_s5 = c("PRKAR1B",
"TRMT10B",
"CUBN",
"TRIO",
"TRIM56",
"SERPINE1",
"COL4A1",
"COL4A2",
"THSD4",
"SIPA1L3",
"ERAP1",
"CASC15",
"LPP",
"NCOR2")

temp_avg = counts(dds, normalized = T)[genes_s5,]
temp_avg[,"Scr-NG-1"] = rowMeans(temp_avg[,c("Scr-NG-1","Scr-NG-2")])
temp_avg[,"Scr-HT-1"] = rowMeans(temp_avg[,c("Scr-HT-1","Scr-HT-2")])
temp_avg[,"LNA1-HT-1"] = rowMeans(temp_avg[,c("LNA1-HT-1","LNA1-HT-2")])
temp_avg[,"LNA4-HT-1"] = rowMeans(temp_avg[,c("LNA4-HT-1","LNA4-HT-2")])
temp_avg = temp_avg[,c(-5:-8)]
colnames(temp_avg) = c("NM-scr","HT-scr","HT-LNA1","HT-LNA2")
temp_avg = log2(temp_avg + 1)

z_temp_avg = (temp_avg - rowMeans(temp_avg)) / apply(temp_avg,1,sd)

pdf("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/LNA_rna_seq/result/TableS5_heatmap.pdf", height = 5, width = 3)
pheatmap(z_temp_avg,
         color = colorRampPalette(c("blue", "yellow"))(100),
         cluster_rows=FALSE,
         cluster_cols=FALSE)
dev.off()

