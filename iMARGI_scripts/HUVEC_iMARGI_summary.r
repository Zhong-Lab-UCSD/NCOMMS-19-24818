##### iMARGI statistics

total_control = 65879147
total_T3d = 53291562
total_T7d = 50031654

total_control_2 = 12495439
total_T7d_2 = 15880128

thres = 2*10^-7

######## Summary plot for number of reads

n = 2
directory = "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/iMARGI_merged_replicates"

######################## Figure 2b
sample=rep(c("0", "3", "7"), 2)
feature=c(rep("Intra" , 3), rep("Inter" , 3))

intra_control_merge = 0.5 * (45041974/total_control + 6140189/total_control_2)
error_intra_control_merge = sd(c(45041974/total_control, 6140189/total_control_2)) / sqrt(2)

intra_T3d_merge = 22872181 / total_T3d

intra_T7d_merge = 0.5 * (14858759/total_T7d + 6713214/total_T7d_2)
error_intra_T7d_merge =  sd(c(14858759/total_T7d , 6713214/total_T7d_2)) / sqrt(2)

inter_control_merge = 0.5 * (20837173/total_control + 6355250/total_control_2)
error_inter_control_merge = sd(c(20837173/total_control , 6355250/total_control_2)) / sqrt(2)

inter_T3d_merge = 30419381 / total_T3d

inter_T7d_merge = 0.5 * (35172895/total_T7d + 9166914/total_T7d_2)
error_inter_T7d_merge = sd(c(35172895/total_T7d , 9166914/total_T7d_2)) / sqrt(2)

value = c(intra_control_merge,	intra_T3d_merge,	intra_T7d_merge, inter_control_merge,	inter_T3d_merge,	inter_T7d_merge)
error = c(error_intra_control_merge, 0, error_intra_T7d_merge, error_inter_control_merge, 0, error_inter_T7d_merge)
data=data.frame(sample,feature,value,error)
data$feature <- factor(data$feature, levels = c("Intra", "Inter"))
data$error_pos = NA
data$error_pos[data$feature == "Intra"] = NA
data$error_pos[data$feature == "Inter"] = data$value[data$feature == "Inter"]

directory = "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/iMARGI_merged_replicates"
png(paste0(directory,"/Figure_2b.png"), width = 8, height = 10, units = "in", res = 200)
ggplot(data, aes(fill=feature, y=value, x=sample)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#ffeb00","blue")) +
  labs(x = "Day", y = "Number of read pairs / total pairs") + 
  geom_errorbar(aes(ymin=error_pos-error, ymax=error_pos+error), width=.2, size = 0.8) +
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

######################## Figure 2e
sample=rep(c("N+M", "H+T"), 2)
feature=c(rep("RNA end in SE" , 2), rep("DNA end in SE" , 2))

value_1 = c((5006772)/total_control,	(3819106 + 3572835)/(total_T3d+total_T7d),
            (2867551)/total_control,	(1176042 + 984091)/(total_T3d+total_T7d))

value_2 = c((929342)/total_control_2,	(1179473)/total_T7d_2,
            (422729)/total_control_2,	(453641)/total_T7d_2)

value_means = (value_1 + value_2) / 2

errors = c()
for (i in 1:4){
  errors = c(errors,sd(c(value_1[i],value_2[i]))/sqrt(2))
}

data=data.frame(sample,feature,value_means,errors)
data$feature <- factor(data$feature, levels = c("RNA end in SE", "DNA end in SE"))
data$sample <- factor(data$sample, levels = c("N+M","H+T"))

# Data update to plot first RNA and then DNA
data$sample = c("N+M","N+M","H+T","H+T")
data$sample <- factor(data$sample, levels = c("N+M","H+T"))
data$feature <- c("a","b","c","d")

directory = "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/iMARGI_merged_replicates"
png(paste0(directory,"/Figure_2e.png"), width = 8, height = 10, units = "in", res = 200)
ggplot(data, aes(fill=feature, y=value_means, x=sample)) +
  geom_bar(stat="identity", position=position_dodge(width=1), width=0.95) +
  scale_fill_manual(values=c("#f8766d","#f8766d","#00bdc4","#00bdc4")) +
  labs(x = "", y = "Number of read pairs / total pairs") + 
  geom_errorbar(aes(ymin=value_means-errors, ymax=value_means+errors), width=.2, size = 1, position=position_dodge(1)) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=34),
        axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none") +
  geom_hline(yintercept = 0.009,  color='black', linetype="dashed", size = 1.5)
dev.off()


######################## Supplementary Figure 4a
sample=rep(c("N+M", "H+T"), 2)
feature=c(rep("RNA end in E" , 2), rep("DNA end in E" , 2))

value_1 = c((6567581)/total_control,	(4681588 + 4287159)/(total_T3d+total_T7d),
            (5208882)/total_control,	(3365440 + 3116817)/(total_T3d+total_T7d))

value_2 = c((1136509)/total_control_2,	(1396885)/total_T7d_2,
            (912612)/total_control_2,	(1050455)/total_T7d_2)

value_means = (value_1 + value_2) / 2

errors = c()
for (i in 1:4){
  errors = c(errors,sd(c(value_1[i],value_2[i]))/sqrt(2))
}

data=data.frame(sample,feature,value_means,errors)
data$feature <- factor(data$feature, levels = c("RNA end in E", "DNA end in E"))
data$sample <- factor(data$sample, levels = c("N+M","H+T"))

# Data update to plot first RNA and then DNA
data$sample = c("N+M","N+M","H+T","H+T")
data$sample <- factor(data$sample, levels = c("N+M","H+T"))
data$feature <- c("a","b","c","d")

directory = "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/iMARGI_merged_replicates"
png(paste0(directory,"/Supplementary_Figure_4a.png"), width = 8, height = 10, units = "in", res = 200)
ggplot(data, aes(fill=feature, y=value_means, x=sample)) +
  geom_bar(stat="identity", position=position_dodge(width=1), width=0.95) +
  scale_fill_manual(values=c("#f8766d","#f8766d","#00bdc4","#00bdc4")) +
  labs(x = "", y = "Number of read pairs / total pairs") + 
  geom_errorbar(aes(ymin=value_means-errors, ymax=value_means+errors), width=.2, size = 1, position=position_dodge(1)) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=34),
        axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none") +
  geom_hline(yintercept = 0.009,  color='black', linetype="dashed", size = 1.5)
dev.off()

######################## Supplementary Figure 4b
sample=rep(c("Day0", "Day7"), 2)
feature=c(rep("Intra" , 2), rep("Inter" , 2))

value_1 = c(1281,	532,
            506,	3253)

value_2 = c(923,	910,
            1286,	2318)

value_means = (value_1 + value_2) / 2

errors = c()
for (i in 1:4){
  errors = c(errors,sd(c(value_1[i],value_2[i]))/sqrt(2))
}

data=data.frame(sample,feature,value_means,errors)

# Data update to plot first RNA and then DNA
data$sample = c("Day0","Day0","Day7","Day7")
data$sample <- factor(data$sample, levels = c("Day0","Day7"))
data$feature <- c("a","b","c","d")

directory = "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/iMARGI_merged_replicates"
png(paste0(directory,"/Supplementary_Figure_4b.png"), width = 8, height = 10, units = "in", res = 200)
ggplot(data, aes(fill=feature, y=value_means, x=sample)) +
  geom_bar(stat="identity", position=position_dodge(width=1), width=0.95) +
  scale_fill_manual(values=c("#f8766d","#f8766d","#00bdc4","#00bdc4")) +
  labs(x = "", y = "Number of SE pairs") + 
  geom_errorbar(aes(ymin=value_means-errors, ymax=value_means+errors), width=.2, size = 1, position=position_dodge(1)) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=34),
        axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none")
dev.off()


######################## Supplementary Figure 4f
sample=rep(c("Nodes"), 2)
feature=c(rep("Day0" , 1), rep("Day7" , 1))

value_1 = c(131,	456)

value_2 = c(303,	419)

value_means = (value_1 + value_2) / 2

errors = c()
for (i in 1:2){
  errors = c(errors,sd(c(value_1[i],value_2[i]))/sqrt(2))
}

data=data.frame(sample,feature,value_means,errors)

# Data update to plot first RNA and then DNA
data$feature <- c("a","b")

directory = "/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/MARGI/iMARGI_merged_replicates"
png(paste0(directory,"/Supplementary_Figure_4f.png"), width = 8, height = 10, units = "in", res = 200)
ggplot(data, aes(fill=feature, y=value_means, x=sample)) +
  geom_bar(stat="identity", position=position_dodge(width=1), width=0.95) +
  scale_fill_manual(values=c("#f8766d","#f8766d")) +
  labs(x = "", y = "Number of SE nodes") + 
  geom_errorbar(aes(ymin=value_means-errors, ymax=value_means+errors), width=.2, size = 1, position=position_dodge(1)) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=34),
        axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none")
dev.off()


