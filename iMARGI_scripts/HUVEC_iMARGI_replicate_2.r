################# Main script for iMARGI analysis
library(GenomicRanges)
library(GenomicAlignments)
library(gdata)
library(dplyr)
library(ggbio)
library(reshape2)
library(DescTools)
library(igraph)
library(plyr)

hg38_chromosomes = c(paste0('chr',c(1:22)),c('chrX','chrY','chrM')) # UCSC
hg38_lengths = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569)
names(hg38_lengths) = hg38_chromosomes

Gr_hg38 <- GRanges(
  seqnames = Rle(names(hg38_lengths[1:23])),
  ranges = IRanges(rep(1,23), end = as.numeric(hg38_lengths[1:23]), names = c(1:length(hg38_lengths[1:23]))),
  strand = Rle(strand('*')))
seqlengths(Gr_hg38) <- hg38_lengths[names(seqlengths(Gr_hg38))]

### Loading super enhancers data
super_enhancers <- read.table('/dataOS/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm/SE_huvec_hg38.bed', stringsAsFactors = FALSE)
colnames(super_enhancers) = c('SE_chr','SE_start','SE_end','SE_id','SE_num')

Gr_super_enhancers <- GRanges(
  seqnames = Rle(super_enhancers[,1]),
  ranges = IRanges(as.numeric(super_enhancers[,2]), end = as.numeric(super_enhancers[,3]), names = c(1:nrow(super_enhancers))),
  strand = Rle(strand('*')),
  SE_id = as.character(super_enhancers[,4]))

super_enhancers_sort = super_enhancers
super_enhancers_sort$SE_chr <- as.factor(super_enhancers_sort$SE_chr)
super_enhancers_sort$SE_chr <- reorder.factor(super_enhancers_sort$SE_chr, levels=hg38_chromosomes)
super_enhancers_sort = arrange(super_enhancers_sort, SE_chr, SE_start, SE_end)

### To Annotate super enhancers
annotation <- read.table('/dataOS/rcalandrelli/MARGI/Homo_sapiens.GRCh38.84.chr.gtf_to_geneTable.tsv', stringsAsFactors = F)
colnames(annotation)=annotation[1,]
annotation=annotation[-1,]
annotation$start = as.numeric(annotation$start)
annotation$end = as.numeric(annotation$end)

Gr_annotation <- GRanges(
  seqnames = Rle(annotation[,1]),
  ranges = IRanges(as.numeric(annotation[,2]), end = as.numeric(annotation[,3]), names = c(1:nrow(annotation))),
  strand = Rle(strand(annotation[,4])),
  gene_id = annotation[,5],
  gene_name = annotation[,6],
  gene_biotype = annotation[,7],
  transcript_id = annotation[,8])

annotate_SE <- function(SE,only_name){
  Gr_SE <- GRanges( 
    seqnames = Rle(as.character(SE[1])),
    ranges = IRanges(as.numeric(SE[2]), end = as.numeric(SE[3]), names = "1"),
    strand = Rle('*'))
  
  temp <- countOverlaps(Gr_annotation, Gr_SE, ignore.strand=TRUE)
  mcols(Gr_annotation)['overlap'] = temp
  annotation <- Gr_annotation[mcols(Gr_annotation)[,"overlap"] == 1]
  if (length(annotation) == 1){
    gene_id=mcols(annotation)['gene_id'][1,1]
    gene_name=mcols(annotation)['gene_name'][1,1]
    gene_biotype=mcols(annotation)['gene_biotype'][1,1]
    if (only_name == TRUE){
      my_annotation <- gene_name
    }
    else {
      my_annotation <- paste0(gene_id,'|',gene_name,'|',gene_biotype)
    }
  } else if (length(annotation) == 0){
    my_annotation = '.'
  } else {
    gene_id=mcols(annotation[1,])['gene_id'][1,1]
    gene_name=mcols(annotation[1,])['gene_name'][1,1]
    gene_biotype=mcols(annotation[1,])['gene_biotype'][1,1]
    if (only_name == TRUE){
      my_annotation <- gene_name
    } else {
      my_annotation <- paste0(gene_id,'|',gene_name,'|',gene_biotype)
    }
    for (i in 2:length(annotation)){
      gene_id=mcols(annotation[i,])['gene_id'][1,1]
      gene_name=mcols(annotation[i,])['gene_name'][1,1]
      gene_biotype=mcols(annotation[i,])['gene_biotype'][1,1]
      if (only_name == TRUE){
        my_annotation <- paste0(my_annotation, ";", gene_name)
      } else {
        my_annotation <- paste0(my_annotation,';',gene_id,'|',gene_name,'|',gene_biotype)
      }
    }
  }
  return(my_annotation)
}

super_enhancers_sort_annotated = cbind(super_enhancers_sort,apply(super_enhancers_sort,1,annotate_SE,only_name=TRUE))
colnames(super_enhancers_sort_annotated)[6] = "genes"

# Add label for plotting network. For each super enhancers as C followed by the chromosome, dash (-) the order of that super enhancer for that chromosome.
labels = c()
for (i in 1:23){
  c = gsub("chr","",hg38_chromosomes[i])
  temp = super_enhancers_sort[which(super_enhancers_sort[,1] == paste0("chr",c)),]
  labels = c(labels,paste0("C",c,"-",seq(1,nrow(temp))))
}
SE_index = seq(1,912)

super_enhancers_sort = cbind(super_enhancers_sort,labels,SE_index)
super_enhancers_sort_annotated = cbind(super_enhancers_sort_annotated,labels,SE_index)

Gr_super_enhancers_sort_annotated <- GRanges(
  seqnames = Rle(super_enhancers_sort_annotated[,1]),
  ranges = IRanges(as.numeric(super_enhancers_sort_annotated[,2]), end = as.numeric(super_enhancers_sort_annotated[,3]), names = c(1:nrow(super_enhancers_sort_annotated))),
  strand = Rle(strand('*')),
  SE_id = as.character(super_enhancers_sort_annotated[,4]),
  SE_genes = as.character(super_enhancers_sort_annotated[,6]),
  SE_labels = as.character(super_enhancers_sort_annotated[,7]),
  SE_index = as.numeric(super_enhancers_sort_annotated[,8]))

# Super enhancers not overlapping any gene
nrow(super_enhancers_sort_annotated[which(super_enhancers_sort_annotated$genes=='.'),]) # 84

# Super enhancers fully embedded into one or multiple genes
overlaps = findOverlaps(Gr_super_enhancers_sort_annotated,Gr_annotation, type="within", ignore.strand=T)
temp = queryHits(overlaps)
se_embedded = unique(temp)
length(se_embedded) # 379

# New data frame to save the embedded super enhancers with updated coordinate
super_enhancers_sort_annotated_new = as.data.frame(Gr_super_enhancers_sort_annotated)
super_enhancers_sort_annotated_new = super_enhancers_sort_annotated_new[,c(-6,-8,-9)]

for (i in se_embedded){
  # Select genes associated with the super enhancer
  genes = strsplit(super_enhancers_sort_annotated_new[i,"SE_genes"],";")[[1]]
  Gr_genes_anno = Gr_annotation[mcols(Gr_annotation)[,"gene_name"] %in% genes]
  # Select the genes that fully contain the super enhancer
  overlaps = findOverlaps(Gr_super_enhancers_sort_annotated[i], Gr_genes_anno, type="within", ignore.strand=T)
  genes_anno_selected = as.data.frame(Gr_genes_anno[subjectHits(overlaps)])
  # Select the gene used to expand the coordinates (if multiple genes) based on the biggest width
  start_update = genes_anno_selected[which.min(genes_anno_selected$start),"start"]
  end_update = genes_anno_selected[which.max(genes_anno_selected$end),"end"]
  # Update the super enhancer coordinates
  super_enhancers_sort_annotated_new[i,"start"] = start_update
  super_enhancers_sort_annotated_new[i,"end"] = end_update
}

# Check that the algorithm worked
# sum(super_enhancers_sort_annotated_new[,"width"]<super_enhancers_sort_annotated_new[,"end"]-super_enhancers_sort_annotated_new[,"start"]+1)

# Remove super enhancers that now are duplicated and re-annotate super enhancers
super_enhancers_sort_annotated_new = super_enhancers_sort_annotated_new[!duplicated(super_enhancers_sort_annotated_new[,c("seqnames","start","end")]),]
super_enhancers_sort_annotated_new[,"SE_genes"] = apply(super_enhancers_sort_annotated_new,1,annotate_SE,only_name=TRUE)
super_enhancers_sort_annotated_new = super_enhancers_sort_annotated_new[,c(-4,-5)]

labels_new = c()
for (i in 1:23){
  c = gsub("chr","",hg38_chromosomes[i])
  temp = super_enhancers_sort_annotated_new[which(super_enhancers_sort_annotated_new[,1] == paste0("chr",c)),]
  labels_new = c(labels_new,paste0("C",c,"-",seq(1,nrow(temp))))
}

SE_index_new = seq(1,nrow(super_enhancers_sort_annotated_new))

super_enhancers_sort_annotated_new = cbind(super_enhancers_sort_annotated_new,labels_new,SE_index_new)
colnames(super_enhancers_sort_annotated_new)[1:3] = c("SE_chr","SE_start","SE_end")
rownames(super_enhancers_sort_annotated_new) = super_enhancers_sort_annotated_new$SE_index_new

Gr_super_enhancers_sort_annotated_new <- GRanges(
  seqnames = Rle(super_enhancers_sort_annotated_new[,1]),
  ranges = IRanges(as.numeric(super_enhancers_sort_annotated_new[,2]), end = as.numeric(super_enhancers_sort_annotated_new[,3]), names = c(1:nrow(super_enhancers_sort_annotated_new))),
  strand = Rle(strand('*')),
  SE_genes = as.character(super_enhancers_sort_annotated_new[,4]),
  labels_new = as.character(super_enhancers_sort_annotated_new[,5]),
  SE_index_new = as.numeric(super_enhancers_sort_annotated_new[,6]))


# Loading MARGI data
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_control_igm____R2filtered'
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered'

hg38_chromosomes = c(paste0('chr',c(1:22)),c('chrX','chrY','chrM')) # UCSC

data_file <- read.table(paste0(directory,'/5____MARGI/annot_exonsIntrons.txt'), stringsAsFactors = F)
data_file <- data_file[which(data_file[,1] %in% hg38_chromosomes & data_file[,5] %in% hg38_chromosomes),]
nrow(data_file[which(data_file[,1]!=data_file[,5]),]) # inter-chromosomal pairs
nrow(data_file[which(data_file[,1]==data_file[,5]),]) # intra-chromosomal pairs

Gr_data_file_RNA <- GRanges(
  seqnames = Rle(data_file[,1]),
  ranges = IRanges(data_file[,2], end = data_file[,3], names = c(1:nrow(data_file))),
  strand = Rle(strand(data_file[,4])),
  DNA_chr = data_file[,5],
  DNA_start = data_file[,6],
  DNA_stop = data_file[,7],
  DNA_strand = data_file[,8],
  annotation_RNA = data_file[,9],
  annotation_DNA = data_file[,10])

Gr_data_file_DNA <- GRanges(
  seqnames = Rle(data_file[,5]),
  ranges = IRanges(data_file[,6], end = data_file[,7], names = c(1:nrow(data_file))),
  strand = Rle(strand(data_file[,8])),
  RNA_chr = data_file[,1],
  RNA_start = data_file[,2],
  RNA_stop = data_file[,3],
  RNA_strand = data_file[,4],
  annotation_RNA = data_file[,9],
  annotation_DNA = data_file[,10])

### Enhancers
enhancers = read.table("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/super_enhancers_call/HUVEC/HUVEC_enhancers_hg38.bed", stringsAsFactors = F)
enhancers = enhancers[which(enhancers[,1] %in% hg38_chromosomes),]
Gr_enhancers = GRanges(
  seqnames = Rle(enhancers[,1]),
  ranges = IRanges(enhancers[,2], end = enhancers[,3], names = c(1:nrow(enhancers))),
  strand = Rle(strand('*')))

overlaps_RNA = countOverlaps(Gr_data_file_RNA,Gr_enhancers, ignore.strand=T)
mcols(Gr_data_file_RNA)["overlap_enhancer"] = overlaps_RNA
Gr_data_file_RNA_enhancer = Gr_data_file_RNA[mcols(Gr_data_file_RNA)[,"overlap_enhancer"] >= 1]
length(Gr_data_file_RNA_enhancer)

overlaps_DNA = countOverlaps(Gr_data_file_DNA,Gr_enhancers, ignore.strand=T)
mcols(Gr_data_file_DNA)["overlap_enhancer"] = overlaps_DNA
Gr_data_file_DNA_enhancer = Gr_data_file_DNA[mcols(Gr_data_file_DNA)[,"overlap_enhancer"] >= 1]
length(Gr_data_file_DNA_enhancer)

pairs_in_enhancers = intersect(names(Gr_data_file_RNA_enhancer),names(Gr_data_file_DNA_enhancer))
length(pairs_in_enhancers) # to extract the number of pairs over enhancers

### super enhancers (all 912)
super_enhancers <- read.table('/dataOS/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm/SE_huvec_hg38.bed', stringsAsFactors = FALSE)
colnames(super_enhancers) = c('SE_chr','SE_start','SE_end','SE_id','SE_num')

Gr_super_enhancers <- GRanges(
  seqnames = Rle(super_enhancers[,1]),
  ranges = IRanges(as.numeric(super_enhancers[,2]), end = as.numeric(super_enhancers[,3]), names = c(1:nrow(super_enhancers))),
  strand = Rle(strand('*')),
  SE_id = as.character(super_enhancers[,4]))

overlaps_RNA = countOverlaps(Gr_data_file_RNA,Gr_super_enhancers, ignore.strand=T)
mcols(Gr_data_file_RNA)["overlap_se_old"] = overlaps_RNA
Gr_data_file_RNA_se_old = Gr_data_file_RNA[mcols(Gr_data_file_RNA)[,"overlap_se_old"] >= 1]
length(Gr_data_file_RNA_se_old)

overlaps_DNA = countOverlaps(Gr_data_file_DNA,Gr_super_enhancers, ignore.strand=T)
mcols(Gr_data_file_DNA)["overlap_se_old"] = overlaps_DNA
Gr_data_file_DNA_se_old = Gr_data_file_DNA[mcols(Gr_data_file_DNA)[,"overlap_se_old"] >= 1]
length(Gr_data_file_DNA_se_old)

pairs_in_SE_old = intersect(names(Gr_data_file_RNA_se_old),names(Gr_data_file_DNA_se_old))
length(pairs_in_SE_old) # to extract the number of pairs over super enhancers

### Heatmap SEs x SEs where each entry is the number of RNA-DNA pairs falling withing the correspondend SEs

overlaps_RNA = countOverlaps(Gr_data_file_RNA,Gr_super_enhancers_sort_annotated_new, ignore.strand=T)
mcols(Gr_data_file_RNA)["overlap"] = overlaps_RNA
Gr_data_file_RNA_SE = Gr_data_file_RNA[mcols(Gr_data_file_RNA)[,"overlap"] >= 1]
length(Gr_data_file_RNA_SE)

overlaps_DNA = countOverlaps(Gr_data_file_DNA,Gr_super_enhancers_sort_annotated_new, ignore.strand=T)
mcols(Gr_data_file_DNA)["overlap"] = overlaps_DNA
Gr_data_file_DNA_SE = Gr_data_file_DNA[mcols(Gr_data_file_DNA)[,"overlap"] >= 1]
length(Gr_data_file_DNA_SE)

pairs_in_SE = intersect(names(Gr_data_file_DNA_SE),names(Gr_data_file_RNA_SE))
length(pairs_in_SE) # to extract the number of pairs over super enhancers

Gr_data_file_RNA_SE_intersect = Gr_data_file_RNA_SE[names(Gr_data_file_RNA_SE) %in% pairs_in_SE]
Gr_data_file_DNA_SE_intersect = Gr_data_file_DNA_SE[names(Gr_data_file_DNA_SE) %in% pairs_in_SE]

mcols(Gr_data_file_RNA_SE_intersect)["SE_num_RNA_end"] = rep("0",length(Gr_data_file_RNA_SE_intersect))
mcols(Gr_data_file_DNA_SE_intersect)["SE_num_DNA_end"] = rep("0",length(Gr_data_file_DNA_SE_intersect))

### RNA end
f_overlaps_RNA = findOverlaps(Gr_data_file_RNA_SE_intersect,Gr_super_enhancers_sort_annotated_new, ignore.strand=TRUE)
f_overlaps_RNA = as.data.frame(f_overlaps_RNA)

tab = table(as.factor(f_overlaps_RNA$queryHits))
# Reads overlapping one SE
temp = as.numeric(names(tab[tab==1]))
# Add information about which super enhancer the read overlaps for those reads
mcols(Gr_data_file_RNA_SE_intersect[temp])["SE_num_RNA_end"] = as.character(f_overlaps_RNA[which(f_overlaps_RNA$queryHits %in% temp),2])

# Work on reads that overlap more than one SE
temp = as.numeric(names(tab[tab>1]))
concatenate_SEs <- function(x){
  temp = f_overlaps_RNA[which(f_overlaps_RNA$queryHits == x),]
  return(paste(temp$subjectHits, collapse = ';'))
}
output = sapply(temp,concatenate_SEs)
mcols(Gr_data_file_RNA_SE_intersect[temp])["SE_num_RNA_end"] = output

### DNA end
f_overlaps_DNA = findOverlaps(Gr_data_file_DNA_SE_intersect,Gr_super_enhancers_sort_annotated_new, ignore.strand=TRUE)
f_overlaps_DNA = as.data.frame(f_overlaps_DNA)

tab = table(as.factor(f_overlaps_DNA$queryHits))
# Reads overlapping one SE
temp = as.numeric(names(tab[tab==1]))
# Add information about which super enhancer the read overlaps for those reads
mcols(Gr_data_file_DNA_SE_intersect[temp])["SE_num_DNA_end"] = as.character(f_overlaps_DNA[which(f_overlaps_DNA$queryHits %in% temp),2])

# Work on reads that overlap more than one SE
temp = as.numeric(names(tab[tab>1]))
concatenate_SEs <- function(x){
  temp = f_overlaps_DNA[which(f_overlaps_DNA$queryHits == x),]
  return(paste(temp$subjectHits, collapse = ';'))
}
output = sapply(temp,concatenate_SEs)
mcols(Gr_data_file_DNA_SE_intersect[temp])["SE_num_DNA_end"] = output

# Count number of interchromosomal pairs within super enhnacers
sum(as.character(seqnames(Gr_data_file_RNA_SE_intersect))!=mcols(Gr_data_file_RNA_SE_intersect)["DNA_chr"][,1])

# Count number of intrachromosomal pairs within super enhnacers (exluding those that are within the same SE)
intra_RNA = Gr_data_file_RNA_SE_intersect[mcols(Gr_data_file_RNA_SE_intersect)[,"DNA_chr"] == seqnames(Gr_data_file_RNA_SE_intersect)]
intra_DNA = Gr_data_file_DNA_SE_intersect[mcols(Gr_data_file_DNA_SE_intersect)[,"RNA_chr"] == seqnames(Gr_data_file_DNA_SE_intersect)]

temp_RNA = mcols(intra_RNA)[,"SE_num_RNA_end"]
temp_DNA = mcols(intra_DNA)[,"SE_num_DNA_end"]

select_not_pseudo <- function(x){
  rna = strsplit(temp_RNA[x],";")[[1]]
  dna = strsplit(temp_DNA[x],";")[[1]]
  return(length(intersect(rna,dna))==0)
}

output = sapply(seq(1,length(temp_RNA)),select_not_pseudo)
sum(output) # number of non pseudo intra pairs within SEs

### Generate matrix
df_data_file_RNA_SE_intersect = as.data.frame(Gr_data_file_RNA_SE_intersect)
df_data_file_DNA_SE_intersect = as.data.frame(Gr_data_file_DNA_SE_intersect)

heatmap_table_global_onlySE_sort = matrix(0,nrow(super_enhancers_sort_annotated_new),nrow(super_enhancers_sort_annotated_new))

for (i in 1:nrow(df_data_file_RNA_SE_intersect)){
  x = as.numeric(strsplit(df_data_file_RNA_SE_intersect[i,"SE_num_RNA_end"],";")[[1]])
  y = as.numeric(strsplit(df_data_file_DNA_SE_intersect[i,"SE_num_DNA_end"],";")[[1]])
  heatmap_table_global_onlySE_sort[x,y] = heatmap_table_global_onlySE_sort[x,y] + 1
}

write.table(heatmap_table_global_onlySE_sort,paste0(directory,"/matrix_only_SE_sort.txt"),sep='\t',row.names=F,col.names = F, quote=F)

# heatmap_table_global_onlySE_sort = read.table(paste0(directory,"/matrix_only_SE_sort.txt"), stringsAsFactors = F)
# rownames(heatmap_table_global_onlySE_sort) = super_enhancers_sort_annotated_new$SE_index_new
# colnames(heatmap_table_global_onlySE_sort) = super_enhancers_sort_annotated_new$SE_index_new

# Only interchromosomal pairs
heatmap_table_global_onlySE_sort_inter = heatmap_table_global_onlySE_sort
for (i in 1:nrow(super_enhancers_sort_annotated_new)){
  SE_row_chr = super_enhancers_sort_annotated_new[i, "SE_chr"]
  for (j in 1:nrow(super_enhancers_sort_annotated_new)){
    SE_col_chr = super_enhancers_sort_annotated_new[j, "SE_chr"]
    if (SE_row_chr == SE_col_chr){
      heatmap_table_global_onlySE_sort_inter[i,j] = 0
    }
  }
}
write.table(heatmap_table_global_onlySE_sort_inter,paste0(directory,"/matrix_only_SE_sort_inter.txt"),sep='\t',row.names=F,col.names = F, quote=F)

# Only intrachromosomal pairs
heatmap_table_global_onlySE_sort_intra = heatmap_table_global_onlySE_sort
for (i in 1:nrow(super_enhancers_sort_annotated_new)){
  SE_row_chr = super_enhancers_sort_annotated_new[i, "SE_chr"]
  for (j in 1:nrow(super_enhancers_sort_annotated_new)){
    SE_col_chr = super_enhancers_sort_annotated_new[j, "SE_chr"]
    if (SE_row_chr != SE_col_chr){
      heatmap_table_global_onlySE_sort_intra[i,j] = 0
    }
  }
}
write.table(heatmap_table_global_onlySE_sort_intra,paste0(directory,"/matrix_only_SE_sort_intra.txt"),sep='\t',row.names=F,col.names = F, quote=F)

### Indexes inter chromosomal interactions
is_inter_chr <- function(x){
  SE_row_chr = super_enhancers_sort_annotated_new[x[1], "SE_chr"]
  SE_col_chr = super_enhancers_sort_annotated_new[x[2], "SE_chr"]
  if (SE_row_chr != SE_col_chr){
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

half_indexes = t(combn(seq(1,875),2))
full_indexes = rbind(half_indexes,cbind(half_indexes[,2],half_indexes[,1]))
output = apply(full_indexes,1,is_inter_chr)
full_indexes = cbind(full_indexes,output)
colnames(full_indexes) = c("Row","Col","is_inter")

###################

total_control_2 = 12495439
total_T7d_2 = 15880128
thres = 2.5*10^-7

# This is to check how many super enhancer pairs there are
edges_number <- function(network_indexes){
  out = c()
  for (i in 1:nrow(network_indexes)){
    x = as.numeric(which(network_indexes[,1]==network_indexes[i,2] & network_indexes[,2]==network_indexes[i,1]))
    if (length(x) != 0){
      if (x > i){
        out = c(out,x)
      }
    }
  }
  if (length(out) > 0){
    out_indexes = network_indexes[-out,]
    return(nrow(out_indexes))
  } else {
    return(nrow(network_indexes))
  }
}

thres=2.5*10^-7
edges_number(which(heatmap_control>thres, arr.ind = T))
edges_number(which(heatmap_T7d>thres, arr.ind = T))

# Day 0
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_control_igm____R2filtered'
heatmap_control = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_control_intra_2 = heatmap_control/total_control_2
indexes_control = which(heatmap_control>thres, arr.ind = T)
length(unique(c(indexes_control[,1],indexes_control[,2]))) # number of nodes

directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_control_igm____R2filtered'
heatmap_control = read.table(paste0(directory,"/matrix_only_SE_sort_intra.txt"), stringsAsFactors = F)
diag(heatmap_control) = 0
rownames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_control_inter_2 = heatmap_control/total_control_2
indexes_control = which(heatmap_control>thres, arr.ind = T)
length(unique(c(indexes_control[,1],indexes_control[,2])))

heatmap_control_2 = heatmap_control_inter_2 + heatmap_control_intra_2
quantile(heatmap_control_2[which(heatmap_control_2>0, arr.ind = T)],0.95)

# Day 7
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered'
heatmap_T7d = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T7d = heatmap_T7d/total_T7d_2
indexes_T7d = which(heatmap_T7d>thres, arr.ind = T)
length(unique(c(indexes_T7d[,1],indexes_T7d[,2])))

directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/replicate_2/HUVEC_T7d_P7_igm____R2filtered'
heatmap_T7d = read.table(paste0(directory,"/matrix_only_SE_sort_intra.txt"), stringsAsFactors = F)
diag(heatmap_T7d) = 0
rownames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T7d = heatmap_T7d/total_T7d_2
indexes_T7d = which(heatmap_T7d>thres, arr.ind = T)
length(unique(c(indexes_T7d[,1],indexes_T7d[,2])))


