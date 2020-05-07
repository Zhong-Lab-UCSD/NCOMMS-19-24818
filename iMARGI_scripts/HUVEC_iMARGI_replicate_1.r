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

### Add field to show if a super enhancer contains only NC RNAs or not
non_coding_RNAs_category = c("snoRNA","snRNA","misc_RNA","antisense","miRNA","lincRNA","processed_transcript","pseudogene")

find_only_nc_RNAs_within_SE <- function(x){
  genes_within_se = strsplit(x["SE_genes"],";")[[1]]
  temp_annotation = annotation[which(annotation$gene_name %in% genes_within_se),"gene_biotype"]
  temp_annotation[grepl("pseudogene",temp_annotation)] = "pseudogene"
  if (sum(temp_annotation %in% non_coding_RNAs_category) == length(temp_annotation)){
    return(1)
  } else {
    return(0)
  }
}

# only_nc_rna=(0,1,2)
# 0 --> only_nc_rna = F
# 1 --> only_nc_rna = T
# 2 --> no genes at all
super_enhancers_sort_annotated_new$only_nc_rna = apply(super_enhancers_sort_annotated_new,1,find_only_nc_RNAs_within_SE)
super_enhancers_sort_annotated_new[which(super_enhancers_sort_annotated_new$SE_genes == "."),"only_nc_rna"] = 2

### Add field to show if a SE overlaps with at least one lincRNA
find_at_least_one_lincRNA_within_SE <- function(x){
  genes_within_se = strsplit(x["SE_genes"],";")[[1]]
  temp_annotation = annotation[which(annotation$gene_name %in% genes_within_se),"gene_biotype"]
  if ("lincRNA" %in% temp_annotation){
    return(1)
  } else {
    return(0)
  }
}
super_enhancers_sort_annotated_new$lincRNA = apply(super_enhancers_sort_annotated_new,1,find_at_least_one_lincRNA_within_SE)

Gr_super_enhancers_sort_annotated_new <- GRanges(
  seqnames = Rle(super_enhancers_sort_annotated_new[,1]),
  ranges = IRanges(as.numeric(super_enhancers_sort_annotated_new[,2]), end = as.numeric(super_enhancers_sort_annotated_new[,3]), names = c(1:nrow(super_enhancers_sort_annotated_new))),
  strand = Rle(strand('*')),
  SE_genes = as.character(super_enhancers_sort_annotated_new[,4]),
  labels_new = as.character(super_enhancers_sort_annotated_new[,5]),
  SE_index_new = as.numeric(super_enhancers_sort_annotated_new[,6]))


# Loading MARGI data
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm'
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered'
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered'

hg38_chromosomes = c(paste0('chr',c(1:22)),c('chrX','chrY','chrM')) # UCSC

data_file <- read.table(paste0(directory,'/annot_exonsIntrons.txt'), stringsAsFactors = F)
data_file <- data_file[which(data_file[,1] %in% hg38_chromosomes & data_file[,5] %in% hg38_chromosomes),]
#nrow(data_file[which(data_file[,1]!=data_file[,5]),]) # to check how many inter-chromosomal read pairs

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

### Enhancer analysis
enhancers = read.table("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/super_enhancers_call/HUVEC/HUVEC_enhancers_hg38.bed", stringsAsFactors = F)
enhancers = enhancers[which(enhancers[,1] %in% hg38_chromosomes),] # 120382426 total length enhancers
Gr_enhancers = GRanges(
  seqnames = Rle(enhancers[,1]),
  ranges = IRanges(enhancers[,2], end = enhancers[,3], names = c(1:nrow(enhancers))),
  strand = Rle(strand('*')))

overlaps_RNA = countOverlaps(Gr_data_file_RNA,Gr_enhancers, ignore.strand=T)
mcols(Gr_data_file_RNA)["overlap_enhancer"] = overlaps_RNA
Gr_data_file_RNA_enhancer = Gr_data_file_RNA[mcols(Gr_data_file_RNA)[,"overlap_enhancer"] >= 1]
length(Gr_data_file_RNA_enhancer) # to extract the number of pairs with RNA end over super enhancers

overlaps_DNA = countOverlaps(Gr_data_file_DNA,Gr_enhancers, ignore.strand=T)
mcols(Gr_data_file_DNA)["overlap_enhancer"] = overlaps_DNA
Gr_data_file_DNA_enhancer = Gr_data_file_DNA[mcols(Gr_data_file_DNA)[,"overlap_enhancer"] >= 1]
length(Gr_data_file_DNA_enhancer) # to extract the number of pairs with RNA end over super enhancers

pairs_in_enhancers = intersect(names(Gr_data_file_RNA_enhancer),names(Gr_data_file_DNA_enhancer))
length(pairs_in_enhancers) # to extract the number of pairs over enhancers

### Super enhancers (all 912)
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
out_9 = length(Gr_data_file_RNA_se_old)

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

overlaps_DNA = countOverlaps(Gr_data_file_DNA,Gr_super_enhancers_sort_annotated_new, ignore.strand=T)
mcols(Gr_data_file_DNA)["overlap"] = overlaps_DNA
Gr_data_file_DNA_SE = Gr_data_file_DNA[mcols(Gr_data_file_DNA)[,"overlap"] >= 1]

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

##### 1. Are there changes in the chromosomal conformation during endoMT?

total_control = 65879147
total_T3d = 53291562
total_T7d = 50031654
thres = 2*10^-7

# Control
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm'
heatmap_control = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_control = heatmap_control/total_control
indexes_control = which(heatmap_control>thres, arr.ind = T)
length(unique(c(indexes_control[,1],indexes_control[,2]))) # number of nodes
plot_network_NEW(indexes_control,paste0(directory,"/control_inter_pairs_above_",as.character(thres),"_network.png"), vertex_size = 1.5)

directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm'
heatmap_control = read.table(paste0(directory,"/matrix_only_SE_sort_intra.txt"), stringsAsFactors = F)
diag(heatmap_control) = 0
rownames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_control = heatmap_control/total_control
indexes_control = which(heatmap_control>thres, arr.ind = T)
length(unique(c(indexes_control[,1],indexes_control[,2])))
plot_network_NEW(indexes_control,paste0(directory,"/control_intra_pairs_above_",as.character(thres),"_network.png"), vertex_size = 1.5, my_layout = layout.auto, edge_width = 1,
                 plot_chr_label = TRUE)

# T3d_P6
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered'
heatmap_T3d = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T3d= heatmap_T3d/total_T3d
indexes_T3d = which(heatmap_T3d>thres, arr.ind = T)
length(unique(c(indexes_T3d[,1],indexes_T3d[,2])))
plot_network_NEW(indexes_T3d,paste0(directory,"/T3d_inter_pairs_above_",as.character(thres),"_network.png"), vertex_size = 1.5)

directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered'
heatmap_T3d = read.table(paste0(directory,"/matrix_only_SE_sort_intra.txt"), stringsAsFactors = F)
diag(heatmap_T3d) = 0
rownames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T3d= heatmap_T3d/total_T3d
indexes_T3d = which(heatmap_T3d>thres, arr.ind = T)
length(unique(c(indexes_T3d[,1],indexes_T3d[,2])))
plot_network_NEW(indexes_T3d,paste0(directory,"/T3d_intra_pairs_above_",as.character(thres),"_network.png"), vertex_size = 1.5, my_layout = layout.auto, edge_width = 1,
                 plot_chr_label = TRUE)

# T7d_P7
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered'
heatmap_T7d = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T7d = heatmap_T7d/total_T7d
indexes_T7d = which(heatmap_T7d>thres, arr.ind = T)
length(unique(c(indexes_T7d[,1],indexes_T7d[,2])))
plot_network_NEW(indexes_T7d,paste0(directory,"/T7d_inter_pairs_above_",as.character(thres),"_network.png"), vertex_size = 1.5)

directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered'
heatmap_T7d = read.table(paste0(directory,"/matrix_only_SE_sort_intra.txt"), stringsAsFactors = F)
diag(heatmap_T7d) = 0
rownames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T7d = heatmap_T7d/total_T7d
indexes_T7d = which(heatmap_T7d>thres, arr.ind = T)
length(unique(c(indexes_T7d[,1],indexes_T7d[,2])))
plot_network_NEW(indexes_T7d,paste0(directory,"/T7d_intra_pairs_above_",as.character(thres),"_network.png"), vertex_size = 1.5, my_layout = layout.auto, edge_width = 1,
                 plot_chr_label = TRUE)


# Boxplot log scale for pairs distribution of interacting super enhancers
# Before this load the indexes with thres=0

temp1 = heatmap_control[indexes_control]
temp_mat1 = cbind(temp1,1)
colnames(temp_mat1) = c("num","Sample")

temp2 = heatmap_T3d[indexes_T3d]
temp_mat2 = cbind(temp2,2)
colnames(temp_mat2) = c("num","Sample")

temp3 = heatmap_T7d[indexes_T7d]
temp_mat3 = cbind(temp3,3)
colnames(temp_mat3) = c("num","Sample")

temp_mat=rbind(temp_mat1,temp_mat2,temp_mat3)
temp_mat = as.data.frame(temp_mat)
temp_mat[which(temp_mat$Sample==1),"Sample"]='Control (Day 0)'
temp_mat[which(temp_mat$Sample==2),"Sample"]='Day 3'
temp_mat[which(temp_mat$Sample==3),"Sample"]='Day 7'
temp_mat[,1] = log10(temp_mat[,1])

x_lab = 0.55
y_lab = 0.05
lab_size = 7

p<-ggplot(temp_mat, aes(x=Sample, y=num)) + 
  geom_boxplot() +
  labs(x="", y=expression("log"["10"]*"(R_SE_total)")) +
  theme_bw() + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20)) + 
  geom_hline(yintercept=log10(8*10^-8), color = "red") +
  annotate("text", x = x_lab, y = log10(8*10^-8)+y_lab, label = "8e-08", color = "red", size = lab_size) +
  geom_hline(yintercept=log10(10^-7), color = "blue") + 
  annotate("text", x = x_lab, y = log10(10^-7)+y_lab, label = "1e-07", color = "blue", size = lab_size) +
  geom_hline(yintercept=log10(2*10^-7), color = "green") +
  annotate("text", x = x_lab, y = log10(2*10^-7)+y_lab, label = "2e-07", color = "green", size = lab_size)

png("/dataOS/rcalandrelli/MARGI/summary_data_NEW/boxplot_interchromosomal.png", height = 10, width = 10, units = "in", res = 200)
print(p)
dev.off()

### Degree of nodes per each network
thres = 2*10^-7
indexes_control = which(heatmap_control>thres, arr.ind = T)
indexes_T3d = which(heatmap_T3d>thres, arr.ind = T)
indexes_T7d = which(heatmap_T7d>thres, arr.ind = T)

indexes_T7d_full = cbind(indexes_T7d,0)
indexes_T7d_full = as.data.frame(indexes_T7d_full)
for (i in 1:nrow(indexes_T7d_full)){
  if (indexes_T7d_full[i,1] | indexes_T7d_full[i,2] %in% hubs_T7d){
    indexes_T7d_full[i,3] = 1
  }
}
colnames(indexes_T7d_full)[3] = "isHub"
indexes_T7d_full = arrange(indexes_T7d_full,row,col)

indexes_T7d_full = cbind(indexes_T7d_full[,1],super_enhancers_sort_annotated_new[indexes_T7d_full[,1],"SE_genes"],
                         indexes_T7d_full[,2],super_enhancers_sort_annotated_new[indexes_T7d_full[,2],"SE_genes"],
                         indexes_T7d_full[,3])
colnames(indexes_T7d_full) = c("SE_rna","SE_rna_genes","SE_dna","SE_dna_genes","isHub")
write.table(indexes_T7d_full,"/dataOS/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered/indexes_T7d_full.txt",sep='\t',row.names = F,col.names = T,quote = F)

thres_DON = 60 # DON = degree of node

# Control
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm'
temp_control = c(as.numeric(indexes_control[,1]),as.numeric(indexes_control[,2]))
tab_control = table(as.factor(temp_control))
tab_control[names(tab_control)==391]
df <- data.frame(
  se=names(tab_control),
  don=as.numeric(tab_control)
)

p1<-ggplot(df, aes(x=don)) + 
  geom_histogram(color="black", fill="darkblue") + 
  labs(x="Degree of nodes", y=expression("Number of nodes")) +
  theme_bw() + 
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(aes(xintercept=thres_DON),
             color="red", size=0.8)
png(paste0(directory,"/control_histogram_above_",as.character(thres),"_DON.png"))
print(p1)
dev.off()

hubs_control = as.numeric(names(tab_control[tab_control>thres_DON]))
indexes_control_hubs = rbind(indexes_control[which(indexes_control[,1] %in% hubs_control),],indexes_control[which(indexes_control[,2] %in% hubs_control),])
indexes_control_hubs = indexes_control_hubs[!duplicated(indexes_control_hubs),]
length(unique(c(indexes_control_hubs[,1],indexes_control_hubs[,2])))
plot_network_NEW(indexes_control_hubs,paste0(directory,"/control_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),".png"),
                 v_label_cex = 1.9, vertex_size = 1.5, vertex_to_increase = hubs_control, increase_v_size=10, color_increased = "#FF0000", mark_T7_hubs = TRUE)

df_hubs = df[which(df[,2]>thres_DON),]
df_hubs = df_hubs[order(df_hubs[,"don"], decreasing = T),]
df_hubs = cbind(super_enhancers_sort_annotated_new[as.character(df_hubs$se),],df_hubs[,"don"])
colnames(df_hubs)[7] = "Degree of node"
write.table(df_hubs,paste0(directory,"/control_hubs_above_",as.character(thres),".txt"),row.names = F, col.names = T, sep='\t', quote = F)

# Enhance the interaction between LINC and SERPINE
plot_network_NEW(indexes_control_hubs,paste0(directory,"/control_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),"_LINC.png"),
                 v_label_cex=1.9, vertex_size = 1.5, vertex_to_increase = hubs_control,increase_v_size=10, color_increased = "#FF0000",
                 edges_to_increase = c(391,145), increase_e_size = 8, label_LINC = TRUE, mark_T7_hubs = TRUE)

# T3d
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered'
temp_T3d = c(as.numeric(indexes_T3d[,1]),as.numeric(indexes_T3d[,2]))
tab_T3d = table(as.factor(temp_T3d))
tab_T3d[names(tab_T3d)==391]
df <- data.frame(
  se=names(tab_T3d),
  don=as.numeric(tab_T3d)
)

p2<-ggplot(df, aes(x=don)) + 
  geom_histogram(color="black", fill="darkblue") + 
  labs(x="Degree of nodes", y=expression("Number of nodes")) +
  theme_bw() + 
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(aes(xintercept=thres_DON),
             color="red", size=0.8)
png(paste0(directory,"/T3d_histogram_above_",as.character(thres),"_DON.png"))
print(p2)
dev.off()

hubs_T3d = as.numeric(names(tab_T3d[tab_T3d>thres_DON]))
indexes_T3d_hubs = rbind(indexes_T3d[which(indexes_T3d[,1] %in% hubs_T3d),],indexes_T3d[which(indexes_T3d[,2] %in% hubs_T3d),])
indexes_T3d_hubs = indexes_T3d_hubs[!duplicated(indexes_T3d_hubs),]
length(unique(c(indexes_T3d_hubs[,1],indexes_T3d_hubs[,2])))
indexes_T3d_hubs_no_MALAT = indexes_T3d_hubs[which(indexes_T3d_hubs[,1] != 557),]
indexes_T3d_hubs_no_MALAT = indexes_T3d_hubs_no_MALAT[which(indexes_T3d_hubs_no_MALAT[,2] != 557),]
plot_network_NEW(indexes_T3d_hubs,paste0(directory,"/T3d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),".png"),
                 v_label_cex=1.9, vertex_size = 1.5, vertex_to_increase = hubs_T3d, v_size_range=c(7,10), color_increased = "#FF0000", mark_T7_hubs = TRUE)

plot_network_NEW(indexes_T3d_hubs_no_MALAT,paste0(directory,"/T3d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),"_no_MALAT.png"),
                 v_label_cex=1.9, vertex_size = 1.5, vertex_to_increase = hubs_T3d[hubs_T3d!=557], v_size_range=c(7,11), color_increased = "#FF0000", mark_T7_hubs = TRUE)

df_hubs = df[which(df[,2]>thres_DON),]
df_hubs = df_hubs[order(df_hubs[,"don"], decreasing = T),]
df_hubs = cbind(super_enhancers_sort_annotated_new[as.character(df_hubs$se),],df_hubs[,"don"])
colnames(df_hubs)[7] = "Degree of node"
write.table(df_hubs,paste0(directory,"/T3d_hubs_above_",as.character(thres),".txt"),row.names = F, col.names = T, sep='\t', quote = F)

# Enhance the interaction between LINC and SERPINE
plot_network_NEW(indexes_T3d_hubs,paste0(directory,"/T3d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),"_LINC.png"),
                 v_label_cex=1.9, vertex_size = 1.5, vertex_to_increase = hubs_T3d, v_size_range=c(7,10), color_increased = "#FF0000", mark_T7_hubs = TRUE,
                 edges_to_increase = list(c(391,145)), increase_e_size = 8, label_LINC = TRUE)

plot_network_NEW(indexes_T3d_hubs_no_MALAT,paste0(directory,"/T3d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),"_no_MALAT_LINC.png"),
                 v_label_cex=1.9, vertex_size = 1.5, vertex_to_increase = hubs_T3d[hubs_T3d!=557], v_size_range=c(7,11), color_increased = "#FF0000", mark_T7_hubs = TRUE,
                 edges_to_increase = list(c(391,145)), increase_e_size = 8, label_LINC = TRUE)

# T7d
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered'
temp_T7d = c(as.numeric(indexes_T7d[,1]),as.numeric(indexes_T7d[,2]))
tab_T7d = table(as.factor(temp_T7d))
tab_T7d[names(tab_T7d)==391]
df <- data.frame(
  se=names(tab_T7d),
  don=as.numeric(tab_T7d)
)

p3<-ggplot(df, aes(x=don)) + 
  geom_histogram(color="black", fill="darkblue") + 
  labs(x="Degree of nodes", y=expression("Number of nodes")) +
  theme_bw() + 
  theme(axis.text=element_text(size=28), axis.title=element_text(size=28),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(aes(xintercept=thres_DON),
             color="red", size=0.8)
png(paste0(directory,"/T7d_histogram_above_",as.character(thres),"_DON.png"))
print(p3)
dev.off()

hubs_T7d = as.numeric(names(tab_T7d[tab_T7d>thres_DON]))
indexes_T7d_hubs = rbind(indexes_T7d[which(indexes_T7d[,1] %in% hubs_T7d),],indexes_T7d[which(indexes_T7d[,2] %in% hubs_T7d),])
indexes_T7d_hubs = indexes_T7d_hubs[!duplicated(indexes_T7d_hubs),]
length(unique(c(indexes_T7d_hubs[,1],indexes_T7d_hubs[,2])))
indexes_T7d_hubs_no_MALAT = indexes_T7d_hubs[which(indexes_T7d_hubs[,1] != 557),]
indexes_T7d_hubs_no_MALAT = indexes_T7d_hubs_no_MALAT[which(indexes_T7d_hubs_no_MALAT[,2] != 557),]
plot_network_NEW(indexes_T7d_hubs,paste0(directory,"/T7d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),".png"),
                 v_label_cex=1.9, vertex_size = 1.5, vertex_to_increase = hubs_T7d, v_size_range=c(7,10), color_increased = "#FF0000")

plot_network_NEW(indexes_T7d_hubs_no_MALAT,paste0(directory,"/T7d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),"_no_MALAT.png"),
                 v_label_cex=1.9, vertex_size = 1.5, vertex_to_increase = hubs_T7d[hubs_T7d!=557], v_size_range=c(7,11), color_increased = "#FF0000")

df_hubs = df[which(df[,2]>thres_DON),]
df_hubs = df_hubs[order(df_hubs[,"don"], decreasing = T),]
df_hubs = cbind(super_enhancers_sort_annotated_new[as.character(df_hubs$se),],df_hubs[,"don"])
colnames(df_hubs)[7] = "Degree of node"
write.table(df_hubs,paste0(directory,"/T7d_hubs_above_",as.character(thres),".txt"),row.names = F, col.names = T, sep='\t', quote = F)

### Plot DON barplots together
png("/dataOS/rcalandrelli/MARGI/summary_data_NEW/DON_barplots.png", width = 18, height = 6, units = "in", res = 100)
plot_grid(p1,p2,p3,nrow = 1,ncol = 3)
dev.off()


# Enhance the interaction between LINC and SERPINE
plot_network_NEW(indexes_T7d_hubs,paste0(directory,"/T7d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),"_LINC.png"),
                 v_label_cex=1.9, vertex_size = 1.5, vertex_to_increase = hubs_T7d, v_size_range=c(7,10), color_increased = "#FF0000",
                 edges_to_increase = matrix(c(391,145),nrow=1,ncol=2,byrow = T), increase_e_size = 8, label_LINC = TRUE)

plot_network_NEW(indexes_T7d_hubs,paste0(directory,"/T7d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),"_LINC_ALL.png"),
                 v_label_cex=1.9, vertex_size = 1.5, vertex_to_increase = hubs_T7d, v_size_range=c(7,10), color_increased = "#FF0000",
                 edges_to_increase = list(hubs_T7d,145), increase_e_size = 8, label_LINC = TRUE)

plot_network_NEW(indexes_T7d_hubs_no_MALAT,paste0(directory,"/T7d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),"_no_MALAT_LINC.png"),
                 v_label_cex=1.9, vertex_size = 1.5, vertex_to_increase = hubs_T7d[hubs_T7d!=557], v_size_range=c(7,11), color_increased = "#FF0000",
                 edges_to_increase = list(c(391,145)), increase_e_size = 8, label_LINC = TRUE)



### Summary for intrachromosomal super enhancer pairs (excluding diagonal)

# Boxplot log scale for pairs distribution of interacting super enhancers
# Before this load the indexes with thres=0

temp1 = heatmap_control[indexes_control]
temp_mat1 = cbind(temp1,1)
colnames(temp_mat1) = c("num","Sample")

temp2 = heatmap_T3d[indexes_T3d]
temp_mat2 = cbind(temp2,2)
colnames(temp_mat2) = c("num","Sample")

temp3 = heatmap_T7d[indexes_T7d]
temp_mat3 = cbind(temp3,3)
colnames(temp_mat3) = c("num","Sample")

temp_mat=rbind(temp_mat1,temp_mat2,temp_mat3)
temp_mat = as.data.frame(temp_mat)
temp_mat[which(temp_mat$Sample==1),"Sample"]='Control (Day 0)'
temp_mat[which(temp_mat$Sample==2),"Sample"]='Day 3'
temp_mat[which(temp_mat$Sample==3),"Sample"]='Day 7'
temp_mat[,1] = log10(temp_mat[,1])

x_lab = 0.55
y_lab = 0.07
lab_size = 7

p<-ggplot(temp_mat, aes(x=Sample, y=num)) + 
  geom_boxplot() +
  labs(x="", y=expression("log"["10"]*"(R_SE_total)")) +
  theme_bw() + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20)) + 
  geom_hline(yintercept=log10(2*10^-7), color = "red") +
  annotate("text", x = x_lab, y = log10(2*10^-7)+y_lab, label = "2e-07", color = "red", size = lab_size) +
  geom_hline(yintercept=log10(3*10^-7), color = "blue") + 
  annotate("text", x = x_lab, y = log10(3*10^-7)+y_lab, label = "3e-07", color = "blue", size = lab_size) +
  geom_hline(yintercept=log10(3*10^-6), color = "green") +
  annotate("text", x = x_lab, y = log10(3*10^-6)+y_lab, label = "3e-06", color = "green", size = lab_size)

png("/dataOS/rcalandrelli/MARGI/summary_data_NEW/boxplot_intrachromosomal.png", height = 10, width = 10, units = "in", res = 200)
print(p)
dev.off()



### Plot network of LINC00607 with two steps of connection
thres = 2*10^-7
# T3d_P6
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered'
heatmap_T3d = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T3d= heatmap_T3d/total_T3d
indexes_T3d = which(heatmap_T3d>thres, arr.ind = T)

# T7d_P7
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered'
heatmap_T7d = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T7d = heatmap_T7d/total_T7d
indexes_T7d = which(heatmap_T7d>thres, arr.ind = T)

### Control (no LINC00607 at all)
indexes_control[which(indexes_control[,1]==145),]
indexes_control[which(indexes_control[,2]==145),]

### T3d
directory = "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered"
temp1 = indexes_T3d[which(indexes_T3d[,1]==145),]
temp2 = indexes_T3d[which(indexes_T3d[,2]==145),]

# Do not care if LINC is on the RNA end or DNA end. Merge the data and remove duplicates to find all the nodes interacting with LINC
indexes_T3d_linc = rbind(temp1,temp2[,c(2,1)])
indexes_T3d_linc = indexes_T3d_linc[!duplicated(indexes_T3d_linc),]

indexes_T3d_linc_2step = indexes_T3d_linc

for (i in indexes_T3d_linc[,2]){
  temp1 = indexes_T3d[which(indexes_T3d[,1]==i & indexes_T3d[,2]!=145),]
  temp2 = indexes_T3d[which(indexes_T3d[,2]==i & indexes_T3d[,1]!=145),]
  if (length(temp1) != 0 | length(temp2) != 0){
    if (length(temp1) == 0){
      temp = temp2[c(2,1)]
      indexes_T3d_linc_2step = rbind(indexes_T3d_linc_2step,temp)
    }
    else if (length(temp2) == 0){
      temp = temp1
      indexes_T3d_linc_2step = rbind(indexes_T3d_linc_2step,temp)
    }
    else {
      if (is.null(dim(temp2))){
        temp = rbind(temp1,temp2[c(2,1)])
      } else {
        temp = rbind(temp1,temp2[,c(2,1)])
      }
      temp = temp[!duplicated(temp),]
      indexes_T3d_linc_2step = rbind(indexes_T3d_linc_2step,temp)
    }
  }
}

length(unique(c(indexes_T3d_linc_2step[,1],indexes_T3d_linc_2step[,2])))

plot_network_NEW(indexes_T3d_linc_2step,paste0(directory,"/T3d_inter_pairs_above_",as.character(thres),"_network_LINC_2step.png"),
                 v_label_cex=1.9, vertex_size = 1.5, 
                 vertex_to_increase = list(indexes_T3d_linc[,2],145), increase_v_size = list(c(7,10),7), increase_v_color = c("#FF0000","yellow"))

# Extract subnetwork where we plot only links that end up into MES markers
mes_markers_index = super_enhancers_sort_annotated_new[which(super_enhancers_sort_annotated_new$SE_genes %like% paste0("%",c("CTGF","VWF","SERPINE1","FN1"),"%")),"SE_index_new"]
se_step1_going_to_mes = c()
for (i in setdiff(mes_markers_index,indexes_T3d_linc_2step[,1])){ # loop over the se that are not directly linked to LINC00607
  se_step1_going_to_mes = rbind(se_step1_going_to_mes,indexes_T3d_linc_2step[which(indexes_T3d_linc_2step[,2]==i),])
}

plot_network_NEW(indexes_T3d_linc_2step,paste0(directory,"/T3d_inter_pairs_above_",as.character(thres),"_network_LINC_2step_MES.png"),
                 v_label_cex=1.9, vertex_size = 1.5,
                 # Vertex to increase:
                 # 1) SE interacting with LINC00607 and also with MES
                 # 2) mes interacting with 1)
                 # 3) LINC00607
                 vertex_to_increase = list(unique(se_step1_going_to_mes[,1]),setdiff(mes_markers_index,indexes_T3d_linc_2step[,1]),145),
                 increase_v_size = list(c(7,10),7,8), increase_v_color = c("#FF0000","green","yellow"),
                 # Edges to increase:
                 # 1) edges between LINC and 1) above
                 # 2) edges between 1) above and mes markers
                 # edges_to_increase = list(cbind(145,unique(se_step1_going_to_mes[,1])),se_step1_going_to_mes), 
                 edges_to_increase = list(se_step1_going_to_mes),
                 # increase_e_size=c(10,6), increase_e_color = c("red","black"))
                 increase_e_size=c(5), increase_e_color = c("black"))

### T7d

# SERPINE1 391
# CTGF 352
# VWF 583

directory = "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered"
temp1 = indexes_T7d[which(indexes_T7d[,1]==145),]
temp2 = indexes_T7d[which(indexes_T7d[,2]==145),]

# Do not care if LINC is on the RNA end or DNA end. Merge the data and remove duplicates to find all the nodes interacting with LINC
indexes_T7d_linc = rbind(temp1,temp2[,c(2,1)])
indexes_T7d_linc = indexes_T7d_linc[!duplicated(indexes_T7d_linc),]

indexes_T7d_linc_2step = indexes_T7d_linc

for (i in indexes_T7d_linc[,2]){
  temp1 = indexes_T7d[which(indexes_T7d[,1]==i & indexes_T7d[,2]!=145),]
  temp2 = indexes_T7d[which(indexes_T7d[,2]==i & indexes_T7d[,1]!=145),]
  if (length(temp1) != 0 | length(temp2) != 0){
    if (length(temp1) == 0){
      temp = temp2[c(2,1)]
      indexes_T7d_linc_2step = rbind(indexes_T7d_linc_2step,temp)
    }
    else if (length(temp2) == 0){
      temp = temp1
      indexes_T7d_linc_2step = rbind(indexes_T7d_linc_2step,temp)
    }
    else {
      if (is.null(dim(temp2))){
        temp = rbind(temp1,temp2[c(2,1)])
      } else {
        temp = rbind(temp1,temp2[,c(2,1)])
      }
      temp = temp[!duplicated(temp),]
      indexes_T7d_linc_2step = rbind(indexes_T7d_linc_2step,temp)
    }
  }
}

length(unique(c(indexes_T7d_linc_2step[,1],indexes_T7d_linc_2step[,2])))

plot_network_NEW(indexes_T7d_linc_2step,paste0(directory,"/T7d_inter_pairs_above_",as.character(thres),"_network_LINC_2step.png"),
                 v_label_cex=1.9, vertex_size = 1.5, 
                 vertex_to_increase = list(indexes_T7d_linc[,2],145), increase_v_size = list(c(7,10),7), increase_v_color = c("#FF0000","yellow"))

# Extract subnetwork where we plot only links that end up into MES markers
mes_markers_index = super_enhancers_sort_annotated_new[which(super_enhancers_sort_annotated_new$SE_genes %like% paste0("%",c("CTGF","VWF","SERPINE1","FN1"),"%")),"SE_index_new"]
se_step1_going_to_mes = c()
for (i in setdiff(mes_markers_index,indexes_T7d_linc_2step[,1])){ # loop over the se that are not directly linked to LINC00607
  se_step1_going_to_mes = rbind(se_step1_going_to_mes,indexes_T7d_linc_2step[which(indexes_T7d_linc_2step[,2]==i),])
}

plot_network_NEW(indexes_T7d_linc_2step,paste0(directory,"/T7d_inter_pairs_above_",as.character(thres),"_network_LINC_2step_MES.png"),
                 v_label_cex=1.9, vertex_size = 1.5,
                 # Vertex to increase:
                 # 1) SE interacting with LINC00607 and also with MES
                 # 2) mes interacting with 1)
                 # 3) LINC00607
                 vertex_to_increase = list(unique(se_step1_going_to_mes[,1]),setdiff(mes_markers_index,indexes_T7d_linc_2step[,1]),145),
                 increase_v_size = list(c(7,10),7,8), increase_v_color = c("#FF0000","green","yellow"),
                 # Edges to increase:
                 # 1) edges between LINC and 1) above
                 # 2) edges between 1) above and mes markers
                 # edges_to_increase = list(cbind(145,unique(se_step1_going_to_mes[,1])),se_step1_going_to_mes), 
                 edges_to_increase = list(se_step1_going_to_mes),
                 # increase_e_size=c(10,6), increase_e_color = c("red","black"))
                 increase_e_size=c(5), increase_e_color = c("black"))


temp_se_to_print = sort(unique(c(as.numeric(indexes_T3d_linc[,2]),as.numeric(indexes_T7d_linc[,2]),c(144,352,391,507,583))))
temp_se_to_print = cbind(temp_se_to_print,super_enhancers_sort_annotated_new[temp_se_to_print,"SE_genes"])
colnames(temp_se_to_print) = c("SE_index","SE_genes")
write.table(temp_se_to_print,"/dataOS/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered/LINC_step1_network_SE.txt",row.names = F,col.names = T,sep='\t',quote = F)



###################################
### Plot coverage track for interactions of LINC00607 RNA SE vs genome (1mb resolution)

library(biovizBase)
hg38IdeogramCyto <- getIdeogram("hg38", cytoband = TRUE)
seqlevels(hg38IdeogramCyto, pruning.mode="coarse") <- hg38_chromosomes[1:24] # keeping only valid chromosomes

PlotTagsProfile <- function(tags, genome_gr, window_size, amplifier, threshold, names, colors, y_lim, plot_chrom) {
  
  genome_window <- tileGenome(seqinfo(genome_gr), tilewidth = window_size, cut.last.tile.in.chrom = T)
  profile_list <- list()
  for(i in 1:length(tags)){
    profile_tmp  <- countOverlaps(genome_window, tags[[i]]) * amplifier[i]
    profile_tmp[profile_tmp > threshold[i]] <- threshold[i]
    profile_list <- c(profile_list, list(profile_tmp))
  }
  
  genome_window@elementMetadata <- setNames(DataFrame(profile_list), nm = names)
  
  if(plot_chrom==TRUE){
    p <- ggplot(genome_gr) + layout_karyogram(cytoband=TRUE) + guides(fill=F) + theme_gray() + xlab("")
    p <- autoplot(Gr_super_enhancers_sort_annotated_new[c(145,391)], layout = "karyogram")
  } else {
    p <- autoplot(Gr_super_enhancers_plot, layout = "karyogram")
  }
  
  for(i in 1:length(tags)){
    p <- p + layout_karyogram(genome_window, aes_string(x = "start", y = names[i]),
                              geom = "line", ylim = y_lim[[i]],
                              color = colors[i], lwd = 0.35)
  }
  return(p)
}


overlaps_RNA = countOverlaps(Gr_data_file_RNA,Gr_super_enhancers_sort_annotated_new[145], ignore.strand=T)
mcols(Gr_data_file_RNA)["overlap"] = overlaps_RNA

overlaps_DNA = countOverlaps(Gr_data_file_DNA,Gr_super_enhancers_sort_annotated_new[145], ignore.strand=T)
mcols(Gr_data_file_DNA)["overlap"] = overlaps_DNA

# Load each file above and run each line separately
Gr_data_file_RNA_LINC_control = Gr_data_file_RNA[mcols(Gr_data_file_RNA)[,"overlap"] >= 1]
Gr_data_file_RNA_LINC_T3d = Gr_data_file_RNA[mcols(Gr_data_file_RNA)[,"overlap"] >= 1]
Gr_data_file_RNA_LINC_T7d = Gr_data_file_RNA[mcols(Gr_data_file_RNA)[,"overlap"] >= 1]

Gr_data_file_DNA_LINC_control = Gr_data_file_DNA[mcols(Gr_data_file_DNA)[,"overlap"] >= 1]
Gr_data_file_DNA_LINC_T3d = Gr_data_file_DNA[mcols(Gr_data_file_DNA)[,"overlap"] >= 1]
Gr_data_file_DNA_LINC_T7d = Gr_data_file_DNA[mcols(Gr_data_file_DNA)[,"overlap"] >= 1]

### Extract interaction reads

### DNA reads interacting with LINC RNA
Gr_data_file_DNA_control <- GRanges(
  seqnames = Rle(mcols(Gr_data_file_RNA_LINC_control)[,"DNA_chr"]),
  ranges = IRanges(mcols(Gr_data_file_RNA_LINC_control)[,"DNA_start"], end = mcols(Gr_data_file_RNA_LINC_control)[,"DNA_stop"], names = c(1:length(Gr_data_file_RNA_LINC_control))),
  strand = Rle(strand(mcols(Gr_data_file_RNA_LINC_control)[,"DNA_strand"])))
Gr_data_file_DNA_T3d <- GRanges(
  seqnames = Rle(mcols(Gr_data_file_RNA_LINC_T3d)[,"DNA_chr"]),
  ranges = IRanges(mcols(Gr_data_file_RNA_LINC_T3d)[,"DNA_start"], end = mcols(Gr_data_file_RNA_LINC_T3d)[,"DNA_stop"], names = c(1:length(Gr_data_file_RNA_LINC_T3d))),
  strand = Rle(strand(mcols(Gr_data_file_RNA_LINC_T3d)[,"DNA_strand"])))
Gr_data_file_DNA_T7d <- GRanges(
  seqnames = Rle(mcols(Gr_data_file_RNA_LINC_T7d)[,"DNA_chr"]),
  ranges = IRanges(mcols(Gr_data_file_RNA_LINC_T7d)[,"DNA_start"], end = mcols(Gr_data_file_RNA_LINC_T7d)[,"DNA_stop"], names = c(1:length(Gr_data_file_RNA_LINC_T7d))),
  strand = Rle(strand(mcols(Gr_data_file_RNA_LINC_T7d)[,"DNA_strand"])))

colors <- c("#E41A1C","#0cad01","#001d68") # red, green, blue     yellow(#fffa00)
pdf("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/endoMT_analysis_NEW/DNA_reads__LINC_RNA.pdf", height = 18, width = 12)
PlotTagsProfile(tags = list(Gr_data_file_DNA_control,Gr_data_file_DNA_T3d,Gr_data_file_DNA_T7d),
                genome_gr = hg38IdeogramCyto, window_size = 1000000,
                amplifier = c(1,1,1), threshold = c(80,80,80), names = c("Control","T3d_P6","T7d_P7"),
                colors = colors, y_lim = list(c(20,100),c(20,100),c(20,100)), plot_chrom = TRUE)
dev.off()

### RNA reads interacting with LINC DNA
Gr_data_file_RNA_control <- GRanges(
  seqnames = Rle(mcols(Gr_data_file_DNA_LINC_control)[,"RNA_chr"]),
  ranges = IRanges(mcols(Gr_data_file_DNA_LINC_control)[,"RNA_start"], end = mcols(Gr_data_file_DNA_LINC_control)[,"RNA_stop"], names = c(1:length(Gr_data_file_DNA_LINC_control))),
  strand = Rle(strand(mcols(Gr_data_file_DNA_LINC_control)[,"RNA_strand"])))
Gr_data_file_RNA_T3d <- GRanges(
  seqnames = Rle(mcols(Gr_data_file_DNA_LINC_T3d)[,"RNA_chr"]),
  ranges = IRanges(mcols(Gr_data_file_DNA_LINC_T3d)[,"RNA_start"], end = mcols(Gr_data_file_DNA_LINC_T3d)[,"RNA_stop"], names = c(1:length(Gr_data_file_DNA_LINC_T3d))),
  strand = Rle(strand(mcols(Gr_data_file_DNA_LINC_T3d)[,"RNA_strand"])))
Gr_data_file_RNA_T7d <- GRanges(
  seqnames = Rle(mcols(Gr_data_file_DNA_LINC_T7d)[,"RNA_chr"]),
  ranges = IRanges(mcols(Gr_data_file_DNA_LINC_T7d)[,"RNA_start"], end = mcols(Gr_data_file_DNA_LINC_T7d)[,"RNA_stop"], names = c(1:length(Gr_data_file_DNA_LINC_T7d))),
  strand = Rle(strand(mcols(Gr_data_file_DNA_LINC_T7d)[,"RNA_strand"])))

colors <- c("#E41A1C","#0cad01","#001d68") # red, green, blue     yellow(#fffa00)
pdf("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/endoMT_analysis_NEW/RNA_reads__LINC_DNA.pdf", height = 18, width = 12)
PlotTagsProfile(tags = list(Gr_data_file_RNA_control,Gr_data_file_RNA_T3d,Gr_data_file_RNA_T7d),
                genome_gr = hg38IdeogramCyto, window_size = 1000000,
                amplifier = c(1,1,1), threshold = c(80,80,80), names = c("Control","T3d_P6","T7d_P7"),
                colors = colors, y_lim = list(c(20,100),c(20,100),c(20,100)), plot_chrom = TRUE)
dev.off()



######## Summary plot for number of reads to be included in the manuscript

### Overall intra- and inter-chromosomal read pairs
sample=rep(c("Ctrl", "Day 3", "Day 7"), 2)
feature=c(rep("Intra-chromosomal" , 3), rep("Inter-chromosomal" , 3))
value = c(0.684,	0.422,	0.297, 0.316,	0.578,	0.703)
data=data.frame(sample,feature,value)
data$feature <- factor(data$feature, levels = c("Intra-chromosomal", "Inter-chromosomal"))

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/summary_data_NEW/barplot_total_reads.png", width = 10, height = 10, units = "in", res = 100)
ggplot(data, aes(fill=feature, y=value, x=sample)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "", y = "Number of read pairs / total pairs") + 
  theme(legend.title = element_blank(),
        text = element_text(size=32),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32))
dev.off()

# Super enhancers coverage over MARGI
sample=rep(c("Ctrl", "Day 3", "Day 7"), 3)
feature=c(rep("RNA end over SEs" , 3), rep("DNA end over SEs" , 3), rep("Both ends over SEs" , 3))
value = c(10574927,	7691827,	7206673,
          6316504,	2967038,	2532498,
          4729619,	1506712, 931370)
data=data.frame(sample,feature,value)
data$feature <- factor(data$feature, levels = c("RNA end over SEs", "DNA end over SEs", "Both ends over SEs"))
perc_label = c("16%","14.4%","14.4%","9.6%","5.6%","5.1%","7.2%","2.8%","1.9%")

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/summary_data_NEW/barplot_SE_coverage.png", width = 15, height = 10, units = "in", res = 100)
ggplot(data, aes(fill=feature, y=value, x=sample)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "", y = "Number of read pairs") + 
  theme(legend.title = element_blank(),
        text = element_text(size=32),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32))
dev.off()

# Intra and inter over SEs
sample=rep(c("Ctrl", "Day 3", "Day 7"), 2)
feature=c(rep("Intra-chromosomal" , 3), rep("Inter-chromosomal" , 3))
value = c(1.5*10^-3,	0.73*10^-3,	0.55*10^-3, 2.7*10^-3,	4.6*10^-3,	5.2*10^-3)
data=data.frame(sample,feature,value)
data$feature <- factor(data$feature, levels = c("Intra-chromosomal", "Inter-chromosomal"))

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/summary_data_NEW/barplot_reads_over_SEs.png", width = 10, height = 10, units = "in", res = 100)
ggplot(data, aes(fill=feature, y=value, x=sample)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "", y = "Number of read pairs / total pairs") + 
  theme(legend.title = element_blank(),
        text = element_text(size=32),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32))
dev.off()
