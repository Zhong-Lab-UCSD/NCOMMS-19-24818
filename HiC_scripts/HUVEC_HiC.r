### HUVEC Hi-C data

hg38_chromosomes = c(paste0('chr',c(1:22)),c('chrX','chrY','chrM')) # UCSC

### Read downsampled data and stats
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled'
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled'
directory = '/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled'

hg38_chromosomes = c(paste0('chr',c(1:22)),c('chrX','chrY','chrM')) # UCSC

data_file <- read.table(paste0(directory,'/HiCfile_paired.bedpe'), stringsAsFactors = F)
data_file <- data_file[which(data_file[,1] %in% hg38_chromosomes & data_file[,4] %in% hg38_chromosomes),]
nrow(data_file)
nrow(data_file[which(data_file[,1]==data_file[,4]),])

# Day 0: 64942544; intra = 59282163; inter = 5660381
# Day 3: 64710653 intra = 58262655; inter = 6447998
# Day 7: 56396995 intra = 52332121; inter = 4064874

##### Plot summary of coverage inter vs intra
sample=rep(c("0", "3", "7"), 2)
feature=c(rep("Intra" , 3), rep("Inter" , 3))
value=c(59282163/64942544, 58262655/64710653, 52332121/56396995, 
        5660381/64942544, 6447998/64710653, 4064874/56396995)
data_hic = data.frame(sample,feature,value)
data_hic$feature <- factor(data_hic$feature, levels = c("Intra", "Inter"))

directory = "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq"
png(paste0(directory,"/Figure_2a.png"), width = 8, height = 10, units = "in", res = 200)
ggplot(data_hic, aes(fill=feature, y=value, x=sample)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#ffeb00","blue")) +
  labs(x = "Day", y = "Number of read pairs / total pairs") + 
  scale_y_continuous(name="Number of read pairs / total pairs", breaks = c(0,0.5,1), labels = c("0","0.5","1")) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=34),
        axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none")
dev.off()

##### Calculate total number of TADs for MoC analysis

# Day 0
n_tads_control = 0
tads_control = matrix(nrow=0,ncol=3)
for (i in hg38_chromosomes[1:24]){
  temp_file = read.table(paste0("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled/tad_analysis/HiCtool_",i,"_topological_domains.txt"), stringsAsFactors = F)
  n_tads_control = n_tads_control+ nrow(temp_file)
  tads_control = rbind(tads_control,cbind(i,temp_file))
}
colnames(tads_control) = c("chr","start","end")
Gr_tads_control <- GRanges(
  seqnames = Rle(tads_control[,1]),
  ranges = IRanges(as.numeric(tads_control[,2]), end = as.numeric(tads_control[,3]), names = c(1:nrow(tads_control))),
  strand = Rle(strand("*")))

# Day 3
n_tads_day3 = 0
tads_day3 = matrix(nrow=0,ncol=3)
for (i in hg38_chromosomes[1:24]){
  temp_file = read.table(paste0("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled/tad_analysis/HiCtool_",i,"_topological_domains.txt"), stringsAsFactors = F)
  n_tads_day3 = n_tads_day3+ nrow(temp_file)
  tads_day3 = rbind(tads_day3,cbind(i,temp_file))
}
colnames(tads_day3) = c("chr","start","end")
Gr_tads_day3 <- GRanges(
  seqnames = Rle(tads_day3[,1]),
  ranges = IRanges(as.numeric(tads_day3[,2]), end = as.numeric(tads_day3[,3]), names = c(1:nrow(tads_day3))),
  strand = Rle(strand("*")))

# Day 7
n_tads_day7 = 0
tads_day7 = matrix(nrow=0,ncol=3)
for (i in hg38_chromosomes[1:24]){
  temp_file = read.table(paste0("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled/tad_analysis/HiCtool_",i,"_topological_domains.txt"), stringsAsFactors = F)
  n_tads_day7 = n_tads_day7+ nrow(temp_file)
  tads_day7 = rbind(tads_day7,cbind(i,temp_file))
}
colnames(tads_day7) = c("chr","start","end")
Gr_tads_day7 <- GRanges(
  seqnames = Rle(tads_day7[,1]),
  ranges = IRanges(as.numeric(tads_day7[,2]), end = as.numeric(tads_day7[,3]), names = c(1:nrow(tads_day7))),
  strand = Rle(strand("*")))


### Find degree of overlap using MoC (measure of concordance)

calculate_MoC <- function(Gr1,Gr2){
  hits <- findOverlaps(Gr1,Gr2)
  overlaps <- pintersect(Gr1[queryHits(hits)], Gr2[subjectHits(hits)])
  F_ij <- as.numeric(width(overlaps))
  P_i <- as.numeric(width(Gr1[queryHits(hits)]))
  Q_j <- as.numeric(width(Gr2[subjectHits(hits)]))
  
  MoC <- 1/(sqrt(length(Gr1) * length(Gr2)) - 1) * (sum(F_ij^2 / (P_i * Q_j)) - 1)
  return(MoC)
}

# Numbers in Supplementary Figure 3c
MoC_control_day3 <- calculate_MoC(Gr_tads_control,Gr_tads_day3)
MoC_control_day7 <- calculate_MoC(Gr_tads_control,Gr_tads_day7)
MoC_day3_day7 <- calculate_MoC(Gr_tads_day3,Gr_tads_day7)


#### Interaction frequency curve over genomic distance
extract_contact_frequency <- function(directory, resolution, max_distance, n_points){
  contact_frequency_list = list()
  start_bin_distance = 1
  end_bin_distance = max_distance / resolution
  step_bin_distance = end_bin_distance / n_points
  for (i in hg38_chromosomes[1:24]){
    contact_matrix = as.matrix(read.table(paste0(directory,i,"_",i,"_",as.character(resolution),".txt"), stringsAsFactors = F))
    contact_frequency = c()
    for (bin_distance in seq(start_bin_distance, end_bin_distance, step_bin_distance)){
      k = 1
      contact_frequency_bin_distance = 0
      while(k + bin_distance < nrow(contact_matrix)){
        contact_frequency_bin_distance = contact_frequency_bin_distance + contact_matrix[k,k+bin_distance]
        k = k + 1
      }
      contact_frequency = c(contact_frequency, contact_frequency_bin_distance / k)
    }
    contact_frequency_list[[i]] = contact_frequency
  }
  return(contact_frequency_list)
}

### Day 0
temp = extract_contact_frequency("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled/normalized_40000/",
                                 40000,
                                 20000000,
                                 10)
contact_frequency_control = colMeans(do.call("rbind", temp))

### Day 3
temp = extract_contact_frequency("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled/normalized_40000/",
                                 40000,
                                 20000000,
                                 10)
contact_frequency_day3 = colMeans(do.call("rbind", temp))

### Day 7
temp = extract_contact_frequency("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled/normalized_40000/",
                                 40000,
                                 20000000,
                                 10)
contact_frequency_day7 = colMeans(do.call("rbind", temp))


### Plot
temp1 = cbind(seq(from = 40000, to = 20000000, by = 40000*50),contact_frequency_control,1)
colnames(temp1) = c("distance","contact","Sample")

temp2 = cbind(seq(from = 40000, to = 20000000, by = 40000*50),contact_frequency_day3,2)
colnames(temp2) = c("distance","contact","Sample")

temp3 = cbind(seq(from = 40000, to = 20000000, by = 40000*50),contact_frequency_day7,3)
colnames(temp3) = c("distance","contact","Sample")

temp_mat = rbind(temp1,temp2,temp3)
temp_mat = as.data.frame(temp_mat)
temp_mat[which(temp_mat$Sample==1),"Sample"]='Day 0'
temp_mat[which(temp_mat$Sample==2),"Sample"]='Day 3'
temp_mat[which(temp_mat$Sample==3),"Sample"]='Day 7'

write.table(temp_mat,"/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/HiC_distance_interaction_plot_downsampled.txt", row.names = F, col.names = T, sep = "\t", quote = F)

library(scales)
png("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Supplementary_Figure_7a.png", height = 8, width = 8, units = "in", res = 200)
ggplot(temp_mat, aes(x=distance, y=contact, color=Sample)) + 
  geom_line(size=1) +
  scale_x_continuous(limits = c(10000,100000000), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(limits = c(0.01,100), trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x="Genomic Distance (bp)", y="Interaction Frequency") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks()
dev.off()


### Calculate proportion of read pairs mapped within TADs

calculate_read_pairs_mapped_within_TADs <- function(directory){
  read_pairs_within_tads = c()
  for (i in hg38_chromosomes[1:24]){
    print(paste0("Calculating read pairs mapped within TADs in ",i))
    input_matrix = read.table(paste0(directory,"/normalized_40000/",i,"_",i,"_40000.txt"), stringsAsFactors = F)
    input_boundaries = read.table(paste0(directory,"/tad_analysis/HiCtool_",i,"_topological_domains.txt"), stringsAsFactors = F)
    input_boundaries_bin = input_boundaries / 40000 + 1
    for (j in 1:nrow(input_boundaries_bin)){
      tad = input_boundaries_bin[j,]
      start_tad = as.numeric(tad[1])
      end_tad = as.numeric(tad[2])
      read_pairs_within_tads = c(read_pairs_within_tads, sum(input_matrix[start_tad:end_tad,start_tad:end_tad]))
    }
  }
  write.table(read_pairs_within_tads,paste0(directory,"/read_pairs_within_tads.txt"), sep="\t", row.names = F, col.names = F, quote = F)
  return(read_pairs_within_tads)
}

read_pairs_within_tads_day0 = calculate_read_pairs_mapped_within_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/control_downsampled")
read_pairs_within_tads_day3 = calculate_read_pairs_mapped_within_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day3_downsampled")
read_pairs_within_tads_day7 = calculate_read_pairs_mapped_within_TADs("/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/HUVEC/HiSeq/Day7_downsampled")

temp_mat1 = cbind(read_pairs_within_tads_day0/190000000,1)
colnames(temp_mat1) = c("num","Sample")

temp_mat2 = cbind(read_pairs_within_tads_day3/190000000,2)
colnames(temp_mat2) = c("num","Sample")

temp_mat3 = cbind(read_pairs_within_tads_day7/190000000,3)
colnames(temp_mat3) = c("num","Sample")

temp_mat=rbind(temp_mat1,temp_mat2,temp_mat3)
temp_mat = as.data.frame(temp_mat)
temp_mat[which(temp_mat$Sample==1),"Sample"]='Day 0'
temp_mat[which(temp_mat$Sample==2),"Sample"]='Day 3'
temp_mat[which(temp_mat$Sample==3),"Sample"]='Day 7'

write.table(temp_mat,"/dataOS/rcalandrelli/HiCtool/HUVEC/HiSeq/HiC_boxplot_proportion_of_reads_mapped_within_TADs.txt", row.names = F, col.names = T, sep = "\t", quote = F)

png("/dataOS/rcalandrelli/HiCtool/HUVEC/HiSeq/Supplementary_Figure_7b.png", height = 8, width = 8, units = "in", res = 200)
ggplot(temp_mat, aes(x=Sample, y=num)) + 
  geom_boxplot() +
  labs(x="", y="Proportion of reads mapped within TADs") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=26),
        axis.text.x = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.ticks.length = unit(0.2, "cm"))
dev.off()

