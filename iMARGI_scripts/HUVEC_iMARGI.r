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

### Heatmap SEs x SEs where each entry is the number of RNA-DNA pairs falling within the corresponding SEs

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

############################################

### Check number of edges (super enhancer pairs): an edge between two super enhancers is counted only once independently if among them there are one- or two-ways read pairs.
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

edges_number(which(heatmap_control>thres, arr.ind = T))
edges_number(which(heatmap_T3d>thres, arr.ind = T))
edges_number(which(heatmap_T7d>thres, arr.ind = T))

### Day 0
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm'
heatmap_control = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_control_inter = heatmap_control/total_control
indexes_control = which(heatmap_control>thres, arr.ind = T)
length(unique(c(indexes_control[,1],indexes_control[,2]))) # number of nodes
# Supplementary Figure 4d
plot_network_NEW(indexes_control,paste0(directory,"/control_inter_pairs_above_",as.character(thres),"_network.png"), vertex_size = 2, edge_width = 0.2, resolution = 400)

directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm'
heatmap_control = read.table(paste0(directory,"/matrix_only_SE_sort_intra.txt"), stringsAsFactors = F)
diag(heatmap_control) = 0
rownames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_control_intra = heatmap_control/total_control
indexes_control = which(heatmap_control>thres, arr.ind = T)
length(unique(c(indexes_control[,1],indexes_control[,2]))) # number of nodes
# Supplementary Figure 4e
plot_network_NEW(indexes_control,paste0(directory,"/control_intra_pairs_above_",as.character(thres),"_network.png"), 
                 vertex_size = 2, my_layout = layout.auto, edge_width = 0.5,
                 plot_chr_label = T, v_label_cex = 0.8, resolution=400, pixels = 2000)

heatmap_control = heatmap_control_inter + heatmap_control_intra
indexes_control = which(heatmap_control>thres, arr.ind = T)
quantile(heatmap_control[which(heatmap_control>0, arr.ind = T)],0.95) # 95th percentile of the normalized counts among SEs

### Day 3
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered'
heatmap_T3d = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T3d= heatmap_T3d/total_T3d
indexes_T3d = which(heatmap_T3d>thres, arr.ind = T)
length(unique(c(indexes_T3d[,1],indexes_T3d[,2]))) # number of nodes
# Supplementary Figure 4d
plot_network_NEW(indexes_T3d,paste0(directory,"/T3d_inter_pairs_above_",as.character(thres),"_network.png"), vertex_size = 2, edge_width = 0.2, resolution = 400)

directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered'
heatmap_T3d = read.table(paste0(directory,"/matrix_only_SE_sort_intra.txt"), stringsAsFactors = F)
diag(heatmap_T3d) = 0
rownames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T3d= heatmap_T3d/total_T3d
indexes_T3d = which(heatmap_T3d>thres, arr.ind = T)
length(unique(c(indexes_T3d[,1],indexes_T3d[,2]))) # number of nodes
# Supplementary Figure 4e
plot_network_NEW(indexes_T3d,paste0(directory,"/T3d_intra_pairs_above_",as.character(thres),"_network.png"), 
                 vertex_size = 2, my_layout = layout.auto, edge_width = 0.5,
                 plot_chr_label = T, v_label_cex = 0.8, resolution = 400, pixels = 2000)

### Day 7
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered'
heatmap_T7d = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T7d = heatmap_T7d/total_T7d
indexes_T7d = which(heatmap_T7d>thres, arr.ind = T)
length(unique(c(indexes_T7d[,1],indexes_T7d[,2]))) # number of nodes
# Supplementary Figure 4d
plot_network_NEW(indexes_T7d,paste0(directory,"/T7d_inter_pairs_above_",as.character(thres),"_network.png"), vertex_size = 2, edge_width = 0.2, resolution = 400)

directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered'
heatmap_T7d = read.table(paste0(directory,"/matrix_only_SE_sort_intra.txt"), stringsAsFactors = F)
diag(heatmap_T7d) = 0
rownames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T7d = heatmap_T7d/total_T7d
indexes_T7d = which(heatmap_T7d>thres, arr.ind = T)
length(unique(c(indexes_T7d[,1],indexes_T7d[,2]))) # number of nodes
# Supplementary Figure 4e
plot_network_NEW(indexes_T7d,paste0(directory,"/T7d_intra_pairs_above_",as.character(thres),"_network.png"), 
                 vertex_size = 2, my_layout = layout.auto, edge_width = 0.5,
                 plot_chr_label = T, v_label_cex = 0.8, resolution = 400, pixels = 2000)


######### Super enhancer hub analysis
thres = 2*10^-7

### Day 0
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm'
heatmap_control = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_control) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_control = heatmap_control/total_control

### Day 3
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered'
heatmap_T3d = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T3d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T3d= heatmap_T3d/total_T3d

### Day 7
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered'
heatmap_T7d = read.table(paste0(directory,"/matrix_only_SE_sort_inter.txt"), stringsAsFactors = F)
rownames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
colnames(heatmap_T7d) = super_enhancers_sort_annotated_new$SE_index_new
heatmap_T7d = heatmap_T7d/total_T7d

indexes_control = which(heatmap_control>thres, arr.ind = T)
indexes_T3d = which(heatmap_T3d>thres, arr.ind = T)
indexes_T7d = which(heatmap_T7d>thres, arr.ind = T)

thres_DON = 60 # DON = degree of node

### Day 0
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_20180613_library2_control_igm'
temp_control = c(as.numeric(indexes_control[,1]),as.numeric(indexes_control[,2]))
tab_control = table(as.factor(temp_control))
tab_control[names(tab_control)==391]
df <- data.frame(
  se=names(tab_control),
  don=as.numeric(tab_control)
)

quantile(df$don, 0.95)

hubs_control = as.numeric(names(tab_control[tab_control>thres_DON]))
indexes_control_hubs = rbind(indexes_control[which(indexes_control[,1] %in% hubs_control),],indexes_control[which(indexes_control[,2] %in% hubs_control),])
indexes_control_hubs = indexes_control_hubs[!duplicated(indexes_control_hubs),]
temp = unique(c(indexes_control_hubs[,1],indexes_control_hubs[,2]))
length(unique(c(indexes_control_hubs[,1],indexes_control_hubs[,2]))) # number of nodes
edges_number(indexes_control_hubs) # number of edges
table(as.factor(super_enhancers_sort_annotated_new[unique(c(indexes_control_hubs[,1],indexes_control_hubs[,2])),"only_nc_rna"])) # number of nodes overlapping with only ncRNAs

temp = super_enhancers_sort_annotated_new[unique(c(indexes_control_hubs[,1],indexes_control_hubs[,2])),]
nrow(temp[which(temp$only_nc_rna==1 & temp$lincRNA==1),]) # number of nodes overlapping with only ncRNAs and with at least one lincRNA

# Figure 3a
plot_network_NEW(indexes_control_hubs,
                 paste0(directory,"/control_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),".png"),
                 v_label_cex = 4, vertex_size = 2, edge_width = 0.2,
                 vertex_to_increase = list(hubs_control,setdiff(hubs_T7d,hubs_control)[setdiff(hubs_T7d,hubs_control) %in% temp]), 
                 increase_v_size=list(11,6), increase_v_color = c("#FF0000","green"),
                 plot_vertex_label = F, resolution = 400)

### Day 3
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T3d_P6_igm____R2filtered'
temp_T3d = c(as.numeric(indexes_T3d[,1]),as.numeric(indexes_T3d[,2]))
tab_T3d = table(as.factor(temp_T3d))
tab_T3d[names(tab_T3d)==391]
df <- data.frame(
  se=names(tab_T3d),
  don=as.numeric(tab_T3d)
)

hubs_T3d = as.numeric(names(tab_T3d[tab_T3d>thres_DON]))
indexes_T3d_hubs = rbind(indexes_T3d[which(indexes_T3d[,1] %in% hubs_T3d),],indexes_T3d[which(indexes_T3d[,2] %in% hubs_T3d),])
indexes_T3d_hubs = indexes_T3d_hubs[!duplicated(indexes_T3d_hubs),]
temp1 = unique(c(indexes_T3d_hubs[,1],indexes_T3d_hubs[,2]))
length(unique(c(indexes_T3d_hubs[,1],indexes_T3d_hubs[,2]))) # number of nodes
edges_number(indexes_T3d_hubs) # number of edges
table(as.factor(super_enhancers_sort_annotated_new[unique(c(indexes_T3d_hubs[,1],indexes_T3d_hubs[,2])),"only_nc_rna"])) # number of nodes overlapping with only ncRNAs

temp = super_enhancers_sort_annotated_new[unique(c(indexes_T3d_hubs[,1],indexes_T3d_hubs[,2])),]
nrow(temp[which(temp$only_nc_rna==1 & temp$lincRNA==1),]) # number of nodes overlapping with only ncRNAs and with at least one lincRNA

# Figure 3a
plot_network_NEW(indexes_T3d_hubs,paste0(directory,"/T3d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),".png"),
                 v_label_cex = 4, vertex_size = 2, edge_width = 0.2,
                 vertex_to_increase = list(hubs_T3d,setdiff(hubs_T7d,hubs_T3d)[setdiff(hubs_T7d,hubs_T3d) %in% temp1],145), 
                 increase_v_size=list(c(6,10),6,6), increase_v_color = c("red","green","yellow"),
                 plot_vertex_label = F, resolution = 400)

### Day 7
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/HUVEC_H_T7d_P7_igm____R2filtered'
temp_T7d = c(as.numeric(indexes_T7d[,1]),as.numeric(indexes_T7d[,2]))
tab_T7d = table(as.factor(temp_T7d))
tab_T7d[names(tab_T7d)==391]
df <- data.frame(
  se=names(tab_T7d),
  don=as.numeric(tab_T7d)
)

hubs_T7d = as.numeric(names(tab_T7d[tab_T7d>thres_DON]))
indexes_T7d_hubs = rbind(indexes_T7d[which(indexes_T7d[,1] %in% hubs_T7d),],indexes_T7d[which(indexes_T7d[,2] %in% hubs_T7d),])
indexes_T7d_hubs = indexes_T7d_hubs[!duplicated(indexes_T7d_hubs),]
length(unique(c(indexes_T7d_hubs[,1],indexes_T7d_hubs[,2]))) # number of nodes
edges_number(indexes_T7d_hubs) # number of edges

table(as.factor(super_enhancers_sort_annotated_new[unique(c(indexes_T7d_hubs[,1],indexes_T7d_hubs[,2])),"only_nc_rna"])) # number of nodes overlapping with only ncRNAs
temp = super_enhancers_sort_annotated_new[unique(c(indexes_T7d_hubs[,1],indexes_T7d_hubs[,2])),]
nrow(temp[which(temp$only_nc_rna==1 & temp$lincRNA==1),]) # number of nodes overlapping with only ncRNAs and with at least one lincRNA

### Check how many lincRNAs are contained in the SEs extracted above
temp = temp[which(temp$only_nc_rna==1 & temp$lincRNA==1),]
temp_genes = c()
for (i in 1:nrow(temp)){
  temp_genes = c(temp_genes,strsplit(temp[i,"SE_genes"],";")[[1]])
}
temp_genes = unique(temp_genes)
temp_anno = c()
for (i in temp_genes){
  temp_anno = c(temp_anno,annotation[which(annotation$gene_name == i),"gene_biotype"][1])
}
table(as.factor(temp_anno))

temp_anno_full = annotation[which(annotation$gene_name %in% temp_genes),]
temp_anno_full_linc = temp_anno_full[which(temp_anno_full$gene_biotype == "lincRNA"),]
write.table(temp_anno_full_linc$gene_name,paste0(directory,"/lincRNAs_in_SEs_hub_network_day7.txt"), row.names = F, col.names = F, sep ="\t", quote = F)

# Figure 3a
plot_network_NEW(indexes_T7d_hubs,paste0(directory,"/T7d_inter_pairs_above_",as.character(thres),"_network_DON_above_",as.character(thres_DON),".png"),
                 v_label_cex=4, vertex_size = 2, edge_width = 0.2,
                 vertex_to_increase = list(hubs_T7d,145), increase_v_size=list(c(6,10),6), increase_v_color = c("red","yellow"),
                 plot_vertex_label = F, resolution = 400)


### Save Table S-HUBS
df_hubs = df[which(df[,2]>thres_DON),]
df_hubs = df_hubs[order(df_hubs[,"don"], decreasing = T),]
df_hubs = cbind(super_enhancers_sort_annotated_new[as.character(df_hubs$se),],df_hubs[,"don"])
colnames(df_hubs)[9] = "Degree of node"
df_hubs = df_hubs[order(df_hubs$SE_index_new),]
df_hubs$emerging_time_point = paste0("Day ", c(3,3,7,7,3,7,7,3,3,3,3,7,3,3,3,0,7,7,7,3,3,7,7,7,3))
df_hubs$emerging_time_point = factor(df_hubs$emerging_time_point, levels = c("Day 0", "Day 3", "Day 7"))
df_hubs = df_hubs[order(df_hubs$emerging_time_point, df_hubs$SE_index_new),]
write.table(df_hubs,paste0(directory,"/T7d_hubs_above_",as.character(thres),".txt"),row.names = F, col.names = T, sep='\t', quote = F)

select_only_mRNA_lincRNA <- function(x){
  temp_genes = strsplit(x["SE_genes"],";")[[1]]
  selected_genes = c()
  for (i in 1:length(temp_genes)){
    if (annotation[which(annotation$gene_name == temp_genes[i]),"gene_biotype"] %in% c("protein_coding","lincRNA")){
      selected_genes = c(selected_genes, temp_genes[i])
    }
  }
  return(paste0(selected_genes,collapse = "; "))
}

df_hubs$selected_genes = apply(df_hubs,1,select_only_mRNA_lincRNA)
df_hubs_to_save = df_hubs[,c(6,10,11)]
colnames(df_hubs_to_save) = c("Index","Emerging time point","Genes contained in this SE")
write.table(df_hubs_to_save,paste0(directory,"/Supplementary_Table_3.txt"),row.names = F, col.names = T, sep='\t', quote = F)


### Degree of nodes and proportions of SEs

# Day 0
df_control <- data.frame(
  se=names(tab_control),
  don=as.numeric(tab_control),
  sample="Day 0"
)

df_control_proportions <- data.frame(matrix(0,ncol = 3, nrow=length(table(as.factor(df_control$don)))))
x <- c("don", "proportion", "sample")
colnames(df_control_proportions) <- x
df_control_proportions$don <- as.numeric(names(table(as.factor(df_control$don))))
df_control_proportions$proportion = as.numeric(table(as.factor(df_control$don))) / sum(as.numeric(table(as.factor(df_control$don))))
df_control_proportions$sample = "Day 0"

# Day 3
df_day3 <- data.frame(
  se=names(tab_T3d),
  don=as.numeric(tab_T3d),
  sample="Day 3"
)

df_day3_proportions <- data.frame(matrix(0,ncol = 3, nrow=length(table(as.factor(df_day3$don)))))
x <- c("don", "proportion", "sample")
colnames(df_day3_proportions) <- x
df_day3_proportions$don <- as.numeric(names(table(as.factor(df_day3$don))))
df_day3_proportions$proportion = as.numeric(table(as.factor(df_day3$don))) / sum(as.numeric(table(as.factor(df_day3$don))))
df_day3_proportions$sample = "Day 3"

df_day7 <- data.frame(
  se=names(tab_T7d),
  don=as.numeric(tab_T7d),
  sample="Day 7"
)

# Day 7
df_day7_proportions <- data.frame(matrix(0,ncol = 3, nrow=length(table(as.factor(df_day7$don)))))
x <- c("don", "proportion", "sample")
colnames(df_day7_proportions) <- x
df_day7_proportions$don <- as.numeric(names(table(as.factor(df_day7$don))))
df_day7_proportions$proportion = as.numeric(table(as.factor(df_day7$don))) / sum(as.numeric(table(as.factor(df_day7$don))))
df_day7_proportions$sample = "Day 7"

df = rbind(df_control,df_day3,df_day7)
df_proportions = rbind(df_control_proportions,df_day3_proportions,df_day7_proportions)

write.table(df_proportions,"/dataOS/rcalandrelli/MARGI/summary_data_NEW/don_proportions_scatter.txt", row.names = F, col.names = T, sep = "\t", quote = F)

png("/dataOS/rcalandrelli/MARGI/summary_data_NEW/Supplementary_Figure_4c.png", width = 12, height = 10, units = "in", res = 200)
ggplot(df_proportions, aes(x=don, y=proportion, color=sample)) + 
  geom_point(aes(size=sample)) + 
  scale_size_manual(values=c(3,3,3)) +
  labs(x="Degree of nodes", y="Proportion of SEs (nodes)") +
  theme_bw() + 
  theme(axis.text = element_text(size=34), 
        axis.title = element_text(size=34),
        legend.text = element_text(size=34),
        legend.title = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm")) +
  scale_x_continuous(trans='log10', breaks = c(1,5,10,50,100,500)) +
  scale_y_continuous(trans='log10', limits = c(0.001,1), breaks = c(1,0.1,0.01,0.001), labels = c("100%","10%","1%","0.1%")) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
dev.off()


########## Plot reads over LINC00607 and SERPINE1 using Gviz (Figure 2c)
library(Gviz)

# Chromosome ideaogram track
itrack_chr2 <- IdeogramTrack(genome="hg38", chromosome="chr2")
itrack_chr7 <- IdeogramTrack(genome="hg38", chromosome="chr7")

# Genomic axis track
gtrack <- GenomeAxisTrack(col = "black")

# Gene track
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
GR <- transcripts(txdb_hg38)
GR.data = data.frame(GR)

createGeneTrack <- function(txdb_object, gene_id, gene_name){
  data(geneModels)
  txdb_temp = select(txdb_object, keys = gene_id, columns=columns(txdb_object), keytype="GENEID")
  txdb_out = txdb_temp[,c("EXONCHROM","EXONSTART","EXONEND","EXONSTRAND","GENEID","EXONID","TXID")]
  txdb_out$WIDTH = txdb_out$EXONEND - txdb_out$EXONSTART + 1
  txdb_out$GENENAME = "LINC00607"
  txdb_out$FEATURE = "lincRNA"
  txdb_out = txdb_out[,c(1:3,8,4,10,5:7,9)]
  for (i in c(1,5:10)){
    txdb_out[,i] = as.factor(txdb_out[,i])
  }
  colnames(txdb_out) = colnames(geneModels)
  
  Gr_txdb_out = GRanges(
    seqnames = Rle(txdb_out$chromosome),
    ranges = IRanges(txdb_out$start, end = txdb_out$end, names = 1:nrow(txdb_out)),
    strand = Rle(strand(txdb_out$strand)))
  Gr_txdb_out_reduced = reduce(Gr_txdb_out)
  txdb_out_reduced = data.frame(Gr_txdb_out_reduced)
  txdb_out_reduced$feature = "x"
  txdb_out_reduced$gene = "x"
  txdb_out_reduced$exon = c(1:length(Gr_txdb_out_reduced))
  txdb_out_reduced$trascript = "1"
  txdb_out_reduced$symbol = gene_name
  colnames(txdb_out_reduced) = colnames(geneModels)
  
  grtrack <- GeneRegionTrack(txdb_out_reduced, genome="hg38", chromosome=as.character(txdb_out_reduced[1,1]), name="Gene", transcriptAnnotation="symbol")
  
  return(grtrack)
}

grtrack_linc607 = createGeneTrack(txdb_hg38, "646324", "LINC00607")
grtrack_serpine1 = createGeneTrack(txdb_hg38, "5054", "SERPINE1")

### Reads track

# Extract reads of LINC00607 interacting with SERPINE1

extract_reads_to_target <- function(Gr, target, target_type){ # Gr is the object that contains both source and target, target_type = "RNA" or "DNA"
  if (target_type == "DNA"){
    Gr_paired = GRanges(
      seqnames = Rle(as.character(mcols(Gr)["DNA_chr"][,1])),
      ranges = IRanges(as.numeric(mcols(Gr)["DNA_start"][,1]), end = as.numeric(mcols(Gr)["DNA_stop"][,1]), names = names(Gr)),
      strand = Rle(strand(Gr)))
    overlaps = countOverlaps(Gr_paired, target, ignore.strand = T)
    mcols(Gr_paired)["overlap"] = overlaps
    Gr_target = Gr_paired[mcols(Gr_paired)[,"overlap"] >= 1]
    Gr_source = Gr[names(Gr_target)]
    mcols(Gr_source)["DNA_chr"] = NULL
    mcols(Gr_source)["DNA_start"] = NULL
    mcols(Gr_source)["DNA_stop"] = NULL
    mcols(Gr_source)["DNA_strand"] = NULL
    mcols(Gr_source)["annotation_RNA"] = NULL
    mcols(Gr_source)["annotation_DNA"] = NULL
    mcols(Gr_source)["overlap"] = NULL
    mcols(Gr_target)["overlap"] = NULL
  }
  else if (target_type == "RNA"){
    Gr_paired = GRanges(
      seqnames = Rle(as.character(mcols(Gr)["RNA_chr"][,1])),
      ranges = IRanges(as.numeric(mcols(Gr)["RNA_start"][,1]), end = as.numeric(mcols(Gr)["RNA_stop"][,1]), names = names(Gr)),
      strand = Rle(strand(Gr)))
    overlaps = countOverlaps(Gr_paired, target, ignore.strand = T)
    mcols(Gr_paired)["overlap"] = overlaps
    Gr_target = Gr_paired[mcols(Gr_paired)[,"overlap"] >= 1]
    Gr_source = Gr[names(Gr_target)]
    mcols(Gr_source)["RNA_chr"] = NULL
    mcols(Gr_source)["RNA_start"] = NULL
    mcols(Gr_source)["RNA_stop"] = NULL
    mcols(Gr_source)["RNA_strand"] = NULL
    mcols(Gr_source)["annotation_RNA"] = NULL
    mcols(Gr_source)["annotation_DNA"] = NULL
    mcols(Gr_source)["overlap"] = NULL
    mcols(Gr_target)["overlap"] = NULL
  }
  return(list(Gr_source,Gr_target))
}

# Reads LINC00607
reads_linc_control_rna = extract_reads_to_target(Gr_data_file_RNA_LINC_control,Gr_super_enhancers_sort_annotated_new[391],target_type = "DNA")[[1]]
reads_linc_day3_rna = extract_reads_to_target(Gr_data_file_RNA_LINC_T3d,Gr_super_enhancers_sort_annotated_new[391],target_type = "DNA")[[1]]
reads_linc_day7_rna = extract_reads_to_target(Gr_data_file_RNA_LINC_T7d,Gr_super_enhancers_sort_annotated_new[391],target_type = "DNA")[[1]]

aTrack.linc.control.rna <- AnnotationTrack(start=c(start(ranges(reads_linc_control_rna)),rep(250000000,length(reads_linc_control_rna))), # 250000000 is just there to make a read that goes beyond the plot boundaries
                                           width=c(width(ranges(reads_linc_control_rna)),rep(100,length(reads_linc_control_rna))),
                                           chromosome=as.character(seqnames(reads_linc_control_rna))[1],
                                           strand=rep(c("*"), length(reads_linc_control_rna)*2),
                                           group=rep(c(1:length(reads_linc_control_rna)),2),
                                           genome="hg38", name="Day 0")
aTrack.linc.day3.rna <- AnnotationTrack(start=c(start(ranges(reads_linc_day3_rna)),rep(250000000,length(reads_linc_day3_rna))), 
                                        width=c(width(ranges(reads_linc_day3_rna)),rep(100,length(reads_linc_day3_rna))),
                                        chromosome=as.character(seqnames(reads_linc_day3_rna))[1],
                                        strand=rep(c("*"), length(reads_linc_day3_rna)*2),
                                        group=rep(c(1:length(reads_linc_day3_rna)),2),
                                        genome="hg38", name="Day 3")
aTrack.linc.day7.rna <- AnnotationTrack(start=c(start(ranges(reads_linc_day7_rna)),rep(250000000,length(reads_linc_day7_rna))), 
                                        width=c(width(ranges(reads_linc_day7_rna)),rep(100,length(reads_linc_day7_rna))),
                                        chromosome=as.character(seqnames(reads_linc_day7_rna))[1],
                                        strand=rep(c("*"), length(reads_linc_day7_rna)*2),
                                        group=rep(c(1:length(reads_linc_day7_rna)),2),
                                        genome="hg38", name="Day 7")

reads_linc_control_dna = extract_reads_to_target(Gr_data_file_DNA_LINC_control,Gr_super_enhancers_sort_annotated_new[391],target_type = "RNA")[[1]]
reads_linc_day3_dna = extract_reads_to_target(Gr_data_file_DNA_LINC_T3d,Gr_super_enhancers_sort_annotated_new[391],target_type = "RNA")[[1]]
reads_linc_day7_dna = extract_reads_to_target(Gr_data_file_DNA_LINC_T7d,Gr_super_enhancers_sort_annotated_new[391],target_type = "RNA")[[1]]

aTrack.linc.control.dna <- AnnotationTrack(start=c(start(ranges(reads_linc_control_dna)),rep(250000000,length(reads_linc_control_dna))), # 250000000 is just there to make a read that goes beyond the plot boundaries
                                           width=c(width(ranges(reads_linc_control_dna)),rep(100,length(reads_linc_control_dna))),
                                           chromosome=as.character(seqnames(reads_linc_control_dna))[1],
                                           strand=rep(c("*"), length(reads_linc_control_dna)*2),
                                           group=rep(c(1:length(reads_linc_control_dna)),2),
                                           genome="hg38", name="Day 0")
aTrack.linc.day3.dna <- AnnotationTrack(start=c(start(ranges(reads_linc_day3_dna)),rep(250000000,length(reads_linc_day3_dna))), 
                                        width=c(width(ranges(reads_linc_day3_dna)),rep(100,length(reads_linc_day3_dna))),
                                        chromosome=as.character(seqnames(reads_linc_day3_dna))[1],
                                        strand=rep(c("*"), length(reads_linc_day3_dna)*2),
                                        group=rep(c(1:length(reads_linc_day3_dna)),2),
                                        genome="hg38", name="Day 3")
aTrack.linc.day7.dna <- AnnotationTrack(start=c(start(ranges(reads_linc_day7_dna)),rep(250000000,length(reads_linc_day7_dna))), 
                                        width=c(width(ranges(reads_linc_day7_dna)),rep(100,length(reads_linc_day7_dna))),
                                        chromosome=as.character(seqnames(reads_linc_day7_dna))[1],
                                        strand=rep(c("*"), length(reads_linc_day7_dna)*2),
                                        group=rep(c(1:length(reads_linc_day7_dna)),2),
                                        genome="hg38", name="Day 7")

reads_linc_control <- unlist(as(list(reads_linc_control_rna, reads_linc_control_dna), "GRangesList"))
reads_linc_day3 <- unlist(as(list(reads_linc_day3_rna, reads_linc_day3_dna), "GRangesList"))
reads_linc_day7 <- unlist(as(list(reads_linc_day7_rna, reads_linc_day7_dna), "GRangesList"))

aTrack.linc.control <- AnnotationTrack(start=c(start(ranges(reads_linc_control)),rep(250000000,length(reads_linc_control))), # 250000000 is just there to make a read that goes beyond the plot boundaries
                                       width=c(width(ranges(reads_linc_control)),rep(100,length(reads_linc_control))),
                                       chromosome=as.character(seqnames(reads_linc_control))[1],
                                       strand=rep(c("*"), length(reads_linc_control)*2),
                                       group=rep(c(1:length(reads_linc_control)),2),
                                       genome="hg38", name="Day 0")
aTrack.linc.day3 <- AnnotationTrack(start=c(start(ranges(reads_linc_day3)),rep(250000000,length(reads_linc_day3))), 
                                    width=c(width(ranges(reads_linc_day3)),rep(100,length(reads_linc_day3))),
                                    chromosome=as.character(seqnames(reads_linc_day3))[1],
                                    strand=rep(c("*"), length(reads_linc_day3)*2),
                                    group=rep(c(1:length(reads_linc_day3)),2),
                                    genome="hg38", name="Day 3")
aTrack.linc.day7 <- AnnotationTrack(start=c(start(ranges(reads_linc_day7)),rep(250000000,length(reads_linc_day7))), 
                                    width=c(width(ranges(reads_linc_day7)),rep(100,length(reads_linc_day7))),
                                    chromosome=as.character(seqnames(reads_linc_day7))[1],
                                    strand=rep(c("*"), length(reads_linc_day7)*2),
                                    group=rep(c(1:length(reads_linc_day7)),2),
                                    genome="hg38", name="Day 7")

# Reads SERPINE1
reads_serpine_control_dna = extract_reads_to_target(Gr_data_file_RNA_LINC_control,Gr_super_enhancers_sort_annotated_new[391],target_type = "DNA")[[2]]
reads_serpine_day3_dna = extract_reads_to_target(Gr_data_file_RNA_LINC_T3d,Gr_super_enhancers_sort_annotated_new[391],target_type = "DNA")[[2]]
reads_serpine_day7_dna = extract_reads_to_target(Gr_data_file_RNA_LINC_T7d,Gr_super_enhancers_sort_annotated_new[391],target_type = "DNA")[[2]]

aTrack.serpine.control.dna <- AnnotationTrack(start=c(start(ranges(reads_serpine_control_dna)),rep(90000000,length(reads_serpine_control_dna))), 
                                              width=c(width(ranges(reads_serpine_control_dna)),rep(100,length(reads_serpine_control_dna))),
                                              chromosome=as.character(seqnames(reads_serpine_control_dna))[1],
                                              strand=rep(c("*"), length(reads_serpine_control_dna)*2),
                                              group=rep(c(1:length(reads_serpine_control_dna)),2),
                                              genome="hg38", name="Day 0")
aTrack.serpine.day3.dna <- AnnotationTrack(start=c(start(ranges(reads_serpine_day3_dna)),rep(90000000,length(reads_serpine_day3_dna))), 
                                           width=c(width(ranges(reads_serpine_day3_dna)),rep(100,length(reads_serpine_day3_dna))),
                                           chromosome=as.character(seqnames(reads_serpine_day3_dna))[1],
                                           strand=rep(c("*"), length(reads_serpine_day3_dna)*2),
                                           group=rep(c(1:length(reads_serpine_day3_dna)),2),
                                           genome="hg38", name="Day 3")
aTrack.serpine.day7.dna <- AnnotationTrack(start=c(start(ranges(reads_serpine_day7_dna)),rep(90000000,length(reads_serpine_day7_dna))), 
                                           width=c(width(ranges(reads_serpine_day7_dna)),rep(100,length(reads_serpine_day7_dna))),
                                           chromosome=as.character(seqnames(reads_serpine_day7_dna))[1],
                                           strand=rep(c("*"), length(reads_serpine_day7_dna)*2),
                                           group=rep(c(1:length(reads_serpine_day7_dna)),2),
                                           genome="hg38", name="Day 7")

reads_serpine_control_rna = extract_reads_to_target(Gr_data_file_DNA_LINC_control,Gr_super_enhancers_sort_annotated_new[391],target_type = "RNA")[[2]]
reads_serpine_day3_rna = extract_reads_to_target(Gr_data_file_DNA_LINC_T3d,Gr_super_enhancers_sort_annotated_new[391],target_type = "RNA")[[2]]
reads_serpine_day7_rna = extract_reads_to_target(Gr_data_file_DNA_LINC_T7d,Gr_super_enhancers_sort_annotated_new[391],target_type = "RNA")[[2]]

aTrack.serpine.control.rna <- AnnotationTrack(start=c(start(ranges(reads_serpine_control_rna)),rep(90000000,length(reads_serpine_control_rna))), 
                                              width=c(width(ranges(reads_serpine_control_rna)),rep(100,length(reads_serpine_control_rna))),
                                              chromosome=as.character(seqnames(reads_serpine_control_rna))[1],
                                              strand=rep(c("*"), length(reads_serpine_control_rna)*2),
                                              group=rep(c(1:length(reads_serpine_control_rna)),2),
                                              genome="hg38", name="Day 0")
aTrack.serpine.day3.rna <- AnnotationTrack(start=c(start(ranges(reads_serpine_day3_rna)),rep(90000000,length(reads_serpine_day3_rna))), 
                                           width=c(width(ranges(reads_serpine_day3_rna)),rep(100,length(reads_serpine_day3_rna))),
                                           chromosome=as.character(seqnames(reads_serpine_day3_rna))[1],
                                           strand=rep(c("*"), length(reads_serpine_day3_rna)*2),
                                           group=rep(c(1:length(reads_serpine_day3_rna)),2),
                                           genome="hg38", name="Day 3")
aTrack.serpine.day7.rna <- AnnotationTrack(start=c(start(ranges(reads_serpine_day7_rna)),rep(90000000,length(reads_serpine_day7_rna))), 
                                           width=c(width(ranges(reads_serpine_day7_rna)),rep(100,length(reads_serpine_day7_rna))),
                                           chromosome=as.character(seqnames(reads_serpine_day7_rna))[1],
                                           strand=rep(c("*"), length(reads_serpine_day7_rna)*2),
                                           group=rep(c(1:length(reads_serpine_day7_rna)),2),
                                           genome="hg38", name="Day 7")

reads_serpine_control <- unlist(as(list(reads_serpine_control_dna, reads_serpine_control_rna), "GRangesList"))
reads_serpine_day3 <- unlist(as(list(reads_serpine_day3_dna, reads_serpine_day3_rna), "GRangesList"))
reads_serpine_day7 <- unlist(as(list(reads_serpine_day7_dna, reads_serpine_day7_rna), "GRangesList"))

aTrack.serpine.control <- AnnotationTrack(start=c(start(ranges(reads_serpine_control)),rep(90000000,length(reads_serpine_control))), 
                                          width=c(width(ranges(reads_serpine_control)),rep(100,length(reads_serpine_control))),
                                          chromosome=as.character(seqnames(reads_serpine_control))[1],
                                          strand=rep(c("*"), length(reads_serpine_control)*2),
                                          group=rep(c(1:length(reads_serpine_control)),2),
                                          genome="hg38", name="Day 0")
aTrack.serpine.day3 <- AnnotationTrack(start=c(start(ranges(reads_serpine_day3)),rep(90000000,length(reads_serpine_day3))), 
                                       width=c(width(ranges(reads_serpine_day3)),rep(100,length(reads_serpine_day3))),
                                       chromosome=as.character(seqnames(reads_serpine_day3))[1],
                                       strand=rep(c("*"), length(reads_serpine_day3)*2),
                                       group=rep(c(1:length(reads_serpine_day3)),2),
                                       genome="hg38", name="Day 3")
aTrack.serpine.day7 <- AnnotationTrack(start=c(start(ranges(reads_serpine_day7)),rep(90000000,length(reads_serpine_day7))), 
                                       width=c(width(ranges(reads_serpine_day7)),rep(100,length(reads_serpine_day7))),
                                       chromosome=as.character(seqnames(reads_serpine_day7))[1],
                                       strand=rep(c("*"), length(reads_serpine_day7)*2),
                                       group=rep(c(1:length(reads_serpine_day7)),2),
                                       genome="hg38", name="Day 7")


# Plot LINC00607 RNA and DNA end
size_panels = c(0.01,0.02,0.005,0.006,0.003,0.012,0.01667,0.03333)
set_params = list(background.panel="#ffffff", col=NULL, fontcolor.title="black", fontfamily.title="arial", fontsize = 14, lwd = 0.8, min.height = 0.5, min.width = 6)

displayPars(itrack_chr2) <- list(fontsize = 12)
displayPars(gtrack) <- list(fontsize = 12)
displayPars(grtrack_linc607) <- list(fill = "black", cex.group = 0.8, arrowHeadWidth=30)
displayPars(aTrack.linc.control.rna) <- list(fill = "#f8766d", col.line="black")
displayPars(aTrack.linc.day3.rna) <- list(fill = "#f8766d", col.line="black")
displayPars(aTrack.linc.day3.dna) <- list(fill = "#00bdc4", col.line="black")
displayPars(aTrack.linc.day7.rna) <- list(fill = "#f8766d", col.line="black")
displayPars(aTrack.linc.day7.dna) <- list(fill = "#00bdc4", col.line="black")

displayPars(grtrack_linc607) <- set_params
displayPars(aTrack.linc.control.rna) <- set_params
displayPars(aTrack.linc.day3.rna) <- set_params
displayPars(aTrack.linc.day3.dna) <- set_params
displayPars(aTrack.linc.day7.rna) <- set_params
displayPars(aTrack.linc.day7.dna) <- set_params

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/summary_data_NEW/LINC00607_reads_Gviz_RNA_DNA.png", width = 6, height = 4, units = "in", res = 300)
plotTracks(list(itrack_chr2, gtrack, grtrack_linc607, 
                aTrack.linc.control.rna, 
                aTrack.linc.day3.rna, aTrack.linc.day3.dna,
                aTrack.linc.day7.rna,aTrack.linc.day7.dna), 
           chromosome = "chr2", from = 214000000, to = 217000000, sizes=size_panels)
dev.off()

# Plot SERPINE1 RNA and DNA end
displayPars(itrack_chr7) <- list(fontsize = 12)
displayPars(gtrack) <- list(fontsize = 12)
displayPars(grtrack_serpine1) <- list(fill = "black", cex.group = 0.8)
displayPars(aTrack.serpine.control.dna) <- list(fill = "#00bdc4", col.line="black")
displayPars(aTrack.serpine.day3.dna) <- list(fill = "#00bdc4", col.line="black")
displayPars(aTrack.serpine.day3.rna) <- list(fill = "#f8766d", col.line="black")
displayPars(aTrack.serpine.day7.dna) <- list(fill = "#00bdc4", col.line="black")
displayPars(aTrack.serpine.day7.rna) <- list(fill = "#f8766d", col.line="black")

displayPars(grtrack_serpine1) <- set_params
displayPars(aTrack.serpine.control.dna) <- set_params
displayPars(aTrack.serpine.day3.dna) <- set_params
displayPars(aTrack.serpine.day3.rna) <- set_params
displayPars(aTrack.serpine.day7.dna) <- set_params
displayPars(aTrack.serpine.day7.rna) <- set_params

png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/summary_data_NEW/SERPINE1_reads_Gviz_RNA_DNA.png", width = 6, height = 4, units = "in", res = 300)
plotTracks(list(itrack_chr7, gtrack, grtrack_serpine1, 
                aTrack.serpine.control.dna, 
                aTrack.serpine.day3.dna, aTrack.serpine.day3.rna,
                aTrack.serpine.day7.dna, aTrack.serpine.day7.rna), 
           chromosome = "chr7", from = 101000000, to = 101300000, sizes=size_panels)
dev.off()


######## Read coverage tracks

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

### Adding replicate 2 for data merging (Gr_data_file_RNA and Gr_data_file_DNA are calculated in the script HUVEC_iMARGI_replicate_2.r)
overlaps_RNA = countOverlaps(Gr_data_file_RNA,Gr_super_enhancers_sort_annotated_new[145], ignore.strand=T)
mcols(Gr_data_file_RNA)["overlap"] = overlaps_RNA

overlaps_DNA = countOverlaps(Gr_data_file_DNA,Gr_super_enhancers_sort_annotated_new[145], ignore.strand=T)
mcols(Gr_data_file_DNA)["overlap"] = overlaps_DNA

# Load each file from replicate 2 and run each line separately
Gr_data_file_RNA_LINC_control_2 = Gr_data_file_RNA[mcols(Gr_data_file_RNA)[,"overlap"] >= 1]
Gr_data_file_RNA_LINC_T7d_2 = Gr_data_file_RNA[mcols(Gr_data_file_RNA)[,"overlap"] >= 1]

Gr_data_file_DNA_LINC_control_2 = Gr_data_file_DNA[mcols(Gr_data_file_DNA)[,"overlap"] >= 1]
Gr_data_file_DNA_LINC_T7d_2 = Gr_data_file_DNA[mcols(Gr_data_file_DNA)[,"overlap"] >= 1]


#################################### Coverage plots of LINC

# Replicate 1
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

# Replicate 2
Gr_data_file_DNA_control_2 <- GRanges(
  seqnames = Rle(mcols(Gr_data_file_RNA_LINC_control_2)[,"DNA_chr"]),
  ranges = IRanges(mcols(Gr_data_file_RNA_LINC_control_2)[,"DNA_start"], end = mcols(Gr_data_file_RNA_LINC_control_2)[,"DNA_stop"], names = c(1:length(Gr_data_file_RNA_LINC_control_2))),
  strand = Rle(strand(mcols(Gr_data_file_RNA_LINC_control_2)[,"DNA_strand"])))
Gr_data_file_DNA_T7d_2 <- GRanges(
  seqnames = Rle(mcols(Gr_data_file_RNA_LINC_T7d_2)[,"DNA_chr"]),
  ranges = IRanges(mcols(Gr_data_file_RNA_LINC_T7d_2)[,"DNA_start"], end = mcols(Gr_data_file_RNA_LINC_T7d_2)[,"DNA_stop"], names = c(1:length(Gr_data_file_RNA_LINC_T7d_2))),
  strand = Rle(strand(mcols(Gr_data_file_RNA_LINC_T7d_2)[,"DNA_strand"])))

# Merging replciates
Gr_data_file_DNA_control_merge = do.call("c",list(Gr_data_file_DNA_control,Gr_data_file_DNA_control_2))
Gr_data_file_DNA_T7d_merge = do.call("c",list(Gr_data_file_DNA_T7d,Gr_data_file_DNA_T7d_2))

total_control_merge = total_control + total_control_2
total_T7d_merge = total_T7d + total_T7d_2

### Saving bed files of LINC-chr7 in order to run homer for calling peaks and check the presence of a peak over SERPINE1
linc_rna_chr7_control = as.data.frame(Gr_data_file_DNA_control)
linc_rna_chr7_control = linc_rna_chr7_control[which(linc_rna_chr7_control$seqnames == "chr7"),c(1,2,3)]
linc_rna_chr7_control[,2] = linc_rna_chr7_control[,2] - 1
write.table(linc_rna_chr7_control,"/dataOS/rcalandrelli/MARGI/endoMT_analysis_NEW/linc_rna_chr7_control.bed", row.names = F, col.names = F, sep='\t', quote = F)

linc_rna_chr7_day3 = as.data.frame(Gr_data_file_DNA_T3d)
linc_rna_chr7_day3 = linc_rna_chr7_day3[which(linc_rna_chr7_day3$seqnames == "chr7"),c(1,2,3)]
linc_rna_chr7_day3[,2] = linc_rna_chr7_day3[,2] - 1
write.table(linc_rna_chr7_day3,"/dataOS/rcalandrelli/MARGI/endoMT_analysis_NEW/linc_rna_chr7_day3.bed", row.names = F, col.names = F, sep='\t', quote = F)

linc_rna_chr7_day7 = as.data.frame(Gr_data_file_DNA_T7d)
linc_rna_chr7_day7 = linc_rna_chr7_day7[which(linc_rna_chr7_day7$seqnames == "chr7"),c(1,2,3)]
linc_rna_chr7_day7[,2] = linc_rna_chr7_day7[,2] - 1
write.table(linc_rna_chr7_day7,"/dataOS/rcalandrelli/MARGI/endoMT_analysis_NEW/linc_rna_chr7_day7.bed", row.names = F, col.names = F, sep='\t', quote = F)

#### Plot
library("karyoploteR", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")

generate_karyoploteR_data <- function(tags,genome_gr, window_size, amplifier, threshold, names) {
  genome_window <- tileGenome(seqinfo(genome_gr), tilewidth = window_size, cut.last.tile.in.chrom = T)
  profile_list <- list()
  for(i in 1:length(tags)){
    profile_tmp  <- countOverlaps(genome_window, tags[[i]]) * amplifier[i]
    profile_tmp[profile_tmp > threshold[i]] <- threshold[i]
    profile_list <- c(profile_list, list(profile_tmp))
  }
  genome_window@elementMetadata <- setNames(DataFrame(profile_list), nm = names)
  return(data.frame(genome_window))
}

cov_data = generate_karyoploteR_data(tags = list(Gr_data_file_DNA_control_merge,Gr_data_file_DNA_T3d,Gr_data_file_DNA_T7d_merge),
                                     genome_gr = hg38IdeogramCyto, 
                                     window_size = 200000,
                                     amplifier = c((1/total_control_merge)*10^8,(1/total_T3d)*10^8,(1/total_T7d_merge)*10^8),
                                     threshold = c(40,40,40),
                                     names = c("Control","T3d_P6","T7d_P7"))
temp_merge_200kb = extract_coverage(cov_data, Gr_super_enhancers_sort_annotated_new[391], samples = c("Control","T3d_P6","T7d_P7"))

### Saving track values for the detail region below
cov_data_to_save = cov_data[which(cov_data$seqnames == "chr7" & cov_data$start >= 80000000 & cov_data$end <= 120000000),]
write.table(cov_data_to_save,"/dataOS/rcalandrelli/MARGI/endoMT_analysis_NEW/coverage_tracks_over_serpine1_200kb.txt", row.names = F, col.names = T, sep = "\t", quote = F)

detail.region <- toGRanges(data.frame("chr7", 80000000, 120000000))

r0_text = 0
r1_text = 0.05
y_text = 0.5

r0_se = r1_text + 0.01
r1_se = r0_se + 0.05

r0_signal = r1_se + 0.05
r1_signal = 0.8
ymin_signal = 0
ymax_signal = 40
line_width = 2

plot_params <- getDefaultPlotParams(plot.type=1)
plot_params$data1inmargin <- 5
plot_params$ideogramheight <- 10
kp <- plotKaryotype(genome="hg38", plot.type=1, plot.params = plot_params, zoom = detail.region, cex=cex_axis)
kpAddBaseNumbers(kp, tick.dist = 5000000, tick.len = 6, tick.col="#4d4d4d", cex=1.2,
                 minor.tick.dist = 1000000, minor.tick.len = 3, minor.tick.col = "#4d4d4d")
kpDataBackground(kp, data.panel = 1, r0=r0_signal, r1=r1_signal)
kpAxis(kp, ymin=ymin_signal, ymax=ymax_signal, r0=r0_signal, r1=r1_signal, col="gray50", cex=cex_axis, numticks = 3)
kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data$Control,
        col="#E41A1C", ymin=ymin_signal, ymax=ymax_signal, r0=r0_signal, r1=r1_signal, lwd=line_width)
kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data$T3d_P6,
        col="#0cad01", ymin=ymin_signal, ymax=ymax_signal, r0=r0_signal, r1=r1_signal, lwd=line_width)
kpLines(kp, chr=cov_data$seqnames, x=rowMeans(cov_data[,2:3]), y=cov_data$T7d_P7,
        col="#001d68", ymin=ymin_signal, ymax=ymax_signal, r0=r0_signal, r1=r1_signal, lwd=line_width)
kpAddLabels(kp, "Read coverage", label.margin=0.08,  side="left", pos=NULL, offset=1, r0=0.5, r1=1, data.panel=1, srt = 90, cex = cex_axis)

# Super enhancer track
kpDataBackground(kp, data.panel = 2, r0=r0_se, r1=r1_se)
kpRect(kp, chr=super_enhancers_sort_annotated_new$SE_chr, x0=super_enhancers_sort_annotated_new$SE_start, x1=super_enhancers_sort_annotated_new$SE_end, y0=0, y1=1, r0=r0_se, r1=r1_se, col=color_body, border=color_border)
kpRect(kp, chr="chr7", x0=super_enhancers_sort_annotated_new[391,"SE_start"], x1=super_enhancers_sort_annotated_new[391,"SE_end"], y0=0, y1=1, r0=r0_se, r1=r1_se, col=color_body_serpine, border=color_border_serpine)
kpText(kp, chr="chr7", x=mean(c(super_enhancers_sort_annotated_new[391,"SE_start"],super_enhancers_sort_annotated_new[391,"SE_end"])), 
       y=y_text, col="#ff00ff", r0=r0_text, r1=r1_text, labels="SERPINE1", cex=1)
kpText(kp, chr="chr7", x=mean(c(start(detail.region),end(detail.region))), 
       y=y_text, col="#000000", r0=r1_se + 0.005, r1=r1_se + 0.04, labels="Super Enhancer Track", cex=0.8)


