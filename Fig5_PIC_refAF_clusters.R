##############################################################################################
## Calculcate PIC values and allele frequencies along the genome between different clusters ##
##############################################################################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(zoo)


colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Read in SNP matrix published in Schreiber et al., 2023 as hapmap file. The genotypes need to be encoded as A, C, G, T

hap<-fread("Hapmap_haploid.hmp.txt")


# Read in the SNP information in numeric transposed format
num<-fread("SNP_numeric.txt")
num<-as.data.frame(num)

hapnames<-names(hap)[-c(1:11)]

# Read in a metadata file generated from Online Resource 1. 
metadata <- read.table("Metadata.txt", head=TRUE, sep="\t")


# Make sure that the accession names in the snpmatrix and the metadata are identical
SamplesKeep <- colnames(hap)[-c(1:11)]
metadata <- metadata[metadata$Line%in%SamplesKeep,]

metadata %<>% arrange(factor(Line, levels = hapnames))

names(hap)[-c(1:11)] == metadata$Line

metadata<-as.data.table(metadata)

num$Line<-names(hap)[-c(1:11)]
num%<>%relocate(Line)
num_t<-transpose_df_Col1(num, "snp")

snps_stat<-num_t

#PIC and MAF per cluster at k = 2

group_1 <- metadata[metadata$Cluster_k2 == 1]$Line
group_2 <- metadata[metadata$Cluster_k2 == 2]$Line

snps_1 <-snps_stat%>%dplyr::select(snp, all_of(group_1))
snps_2 <-snps_stat%>%dplyr::select(snp, all_of(group_2))

snps_1_1 <- snps_1
snps_2_1 <- snps_2

snps_1_1$AF1<-((rowSums(snps_1 == 0, na.rm = TRUE)*2)+(rowSums(snps_1 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_1[,2:ncol(snps_1)]))*2)
snps_2_1$AF1<-((rowSums(snps_2 == 0, na.rm = TRUE)*2)+(rowSums(snps_2 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_2[,2:ncol(snps_2)]))*2)

snps_1_1$AF2<-((rowSums(snps_1 == 2, na.rm = TRUE)*2)+(rowSums(snps_1 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_1[,2:ncol(snps_1)]))*2)
snps_2_1$AF2<-((rowSums(snps_2 == 2, na.rm = TRUE)*2)+(rowSums(snps_2 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_2[,2:ncol(snps_2)]))*2)

snps_1_1$MAF<-pmin(snps_1_1$AF1, snps_1_1$AF2)
snps_2_1$MAF<-pmin(snps_2_1$AF1, snps_2_1$AF2)

snps_1_1$PIC<-1-(snps_1_1$MAF^2+(1-snps_1_1$MAF)^2)-(2*(snps_1_1$MAF^2)*(1-snps_1_1$MAF)^2)
snps_2_1$PIC<-1-(snps_2_1$MAF^2+(1-snps_2_1$MAF)^2)-(2*(snps_2_1$MAF^2)*(1-snps_2_1$MAF)^2)

snps_1_1$Cluster <- 1
snps_2_1$Cluster <- 2

snps_stats<-rbind(snps_1_1[, c(1,(ncol(snps_1_1)-4):ncol(snps_1_1))], snps_2_1[, c(1,(ncol(snps_2_1)-4):ncol(snps_2_1))])  

snps_stats<-as.data.table(snps_stats)  
snps_stats$Cluster<-as.factor(snps_stats$Cluster)


# Calculate PIC and reference AF in rolling windows of 10000 SNPs
snppos<-hap[,1:4]
names(snppos)[1]<-"snp"

snps_stats_pos<-left_join(snps_stats, snppos, by="snp")

snps_stats_pos$PIC_roll<-rollmean(snps_stats_pos$PIC,10000,fill = list(NA, NULL, NA))
snps_stats_pos$AF_roll<-rollmean(snps_stats_pos$AF1, 10000,fill = list(NA, NULL, NA))


data_vline_p1<-data.frame(chrom = c(  "5H", "5H" ),
                          vline = c(  68.78, 320.04))


data_vline_centro<-data.frame(chrom = c("1H", "2H", "3H", "4H", "5H",  "6H", "7H" ),
                              vline = c(211.460194, 298.511764, 270.896959, 273.503980, 214.765154, 252.168825, 329.584237))

scaleFUN <- function(x) sprintf("%.2f", x)

p1<-ggplot(snps_stats_pos, aes(x=pos/1000000, y=PIC_roll, color=Cluster)) + geom_point(size=0.3) +     theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200), labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) +   scale_colour_manual(values = colorBlindGrey8)+  facet_grid(~chrom, scales='free_x', space="free_x")  + ylab("PIC") + labs(color="Release period")+ theme(legend.position = "none")+ theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+scale_y_continuous(labels=scaleFUN)+ geom_vline(data_vline_p1, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red")+geom_vline(data_vline_centro, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="lightgrey")

p2<-ggplot(snps_stats_pos, aes(x=pos/1000000, y=AF_roll, color=Cluster)) + geom_point(size=0.3) +    theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200), labels = function(x) format(x, scientific = FALSE))+
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) +  scale_colour_manual(values = colorBlindGrey8)+ facet_grid(~chrom, scales='free_x', space="free_x")+ xlab("Position (Mbp)")  + ylab("Ref. allele freq.") + labs(color="Cluster")+theme(legend.position = "none") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+scale_y_continuous(labels=scaleFUN)+ geom_vline(data_vline_p1, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red")+geom_vline(data_vline_centro, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="lightgrey")


#PIC and MAF per cluster at k = 3

snps_stats<-num_t

group_1 <- metadata[metadata$Cluster_k3 == 1]$Line
group_2 <- metadata[metadata$Cluster_k3 == 2]$Line
group_3 <- metadata[metadata$Cluster_k3 == 3]$Line

snps_1 <-snps_stat%>%select(snp, all_of(group_1))
snps_2 <-snps_stat%>%select(snp, all_of(group_2))
snps_3 <-snps_stat%>%select(snp, all_of(group_3))

snps_1_1 <- snps_1
snps_2_1 <- snps_2
snps_3_1 <- snps_3

snps_1_1$AF1<-((rowSums(snps_1 == 0, na.rm = TRUE)*2)+(rowSums(snps_1 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_1[,2:ncol(snps_1)]))*2)
snps_2_1$AF1<-((rowSums(snps_2 == 0, na.rm = TRUE)*2)+(rowSums(snps_2 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_2[,2:ncol(snps_2)]))*2)
snps_3_1$AF1<-((rowSums(snps_3 == 0, na.rm = TRUE)*2)+(rowSums(snps_3 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_3[,2:ncol(snps_3)]))*2)

snps_1_1$AF2<-((rowSums(snps_1 == 2, na.rm = TRUE)*2)+(rowSums(snps_1 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_1[,2:ncol(snps_1)]))*2)
snps_2_1$AF2<-((rowSums(snps_2 == 2, na.rm = TRUE)*2)+(rowSums(snps_2 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_2[,2:ncol(snps_2)]))*2)
snps_3_1$AF2<-((rowSums(snps_3 == 2, na.rm = TRUE)*2)+(rowSums(snps_3 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_3[,2:ncol(snps_3)]))*2)

snps_1_1$MAF<-pmin(snps_1_1$AF1, snps_1_1$AF2)
snps_2_1$MAF<-pmin(snps_2_1$AF1, snps_2_1$AF2)
snps_3_1$MAF<-pmin(snps_3_1$AF1, snps_3_1$AF2)

snps_1_1$PIC<-1-(snps_1_1$MAF^2+(1-snps_1_1$MAF)^2)-(2*(snps_1_1$MAF^2)*(1-snps_1_1$MAF)^2)
snps_2_1$PIC<-1-(snps_2_1$MAF^2+(1-snps_2_1$MAF)^2)-(2*(snps_2_1$MAF^2)*(1-snps_2_1$MAF)^2)
snps_3_1$PIC<-1-(snps_3_1$MAF^2+(1-snps_3_1$MAF)^2)-(2*(snps_3_1$MAF^2)*(1-snps_3_1$MAF)^2)

snps_1_1$Cluster <- 1
snps_2_1$Cluster <- 2
snps_3_1$Cluster <- 3

snps_stats<-rbind(snps_1_1[, c(1,(ncol(snps_1_1)-4):ncol(snps_1_1))], snps_2_1[, c(1,(ncol(snps_2_1)-4):ncol(snps_2_1))], snps_3_1[, c(1,(ncol(snps_3_1)-4):ncol(snps_3_1))])  

snps_stats<-as.data.table(snps_stats)  
snps_stats$Cluster<-as.factor(snps_stats$Cluster)
 


# Calculate PIC and reference AF in rolling windows of 10000 SNPs
snppos<-hap[,1:4]
names(snppos)[1]<-"snp"

snps_stats_pos<-left_join(snps_stats, snppos, by="snp")

snps_stats_pos$PIC_roll<-rollmean(snps_stats_pos$PIC,10000,fill = list(NA, NULL, NA))
snps_stats_pos$AF_roll<-rollmean(snps_stats_pos$AF1, 10000,fill = list(NA, NULL, NA))


##### Combined plot of PIC and refaf
data_vline_p1<-data.frame(chrom = c( "3H","3H", "5H", "5H" ),
                          vline = c(  89.7, 444, 68.78, 320.04))


scaleFUN <- function(x) sprintf("%.2f", x)

p3<-ggplot(snps_stats_pos, aes(x=pos/1000000, y=PIC_roll, color=Cluster)) + geom_point(size=0.3) +     theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200), labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) +   scale_colour_manual(values = colorBlindGrey8)+  facet_grid(~chrom, scales='free_x', space="free_x")  + ylab("PIC") + labs(color="Release period")+ theme(legend.position = "none")+ theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+scale_y_continuous(labels=scaleFUN)+ geom_vline(data_vline_p1, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red")+geom_vline(data_vline_centro, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="lightgrey")

p4<-ggplot(snps_stats_pos, aes(x=pos/1000000, y=AF_roll, color=Cluster)) + geom_point(size=0.3) +    theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200), labels = function(x) format(x, scientific = FALSE))+
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) +  scale_colour_manual(values = colorBlindGrey8)+ facet_grid(~chrom, scales='free_x', space="free_x")+ xlab("Position (Mbp)")  + ylab("Ref. allele freq.") + labs(color="Cluster")+ theme(legend.position = "none")+scale_y_continuous(labels=scaleFUN)+ geom_vline(data_vline_p1, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red")+geom_vline(data_vline_centro, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="lightgrey")


p5<-ggplot(snps_stats_pos, aes(x=pos/1000000, y=AF_roll, color=Cluster)) + geom_point(size=0.3) +    theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200), labels = function(x) format(x, scientific = FALSE))+
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) +  scale_colour_manual(values = colorBlindGrey8)+ facet_grid(~chrom, scales='free_x', space="free_x")+ xlab("Position (Mbp)")  + ylab("Ref. allele freq.") + labs(color="Cluster")+ theme(legend.position ="bottom",legend.title = element_text(size=9),legend.text = element_text(size=8), legend.spacing.y = unit(-0.1, 'cm'))+scale_y_continuous(labels=scaleFUN)+ geom_vline(data_vline_p1, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red")+geom_vline(data_vline_centro, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="lightgrey")+ guides(colour = guide_legend(override.aes = list(size=1)))+ labs(color = "k-means cluster")

legend<-get_legend(p5)
legendplot<-as_ggplot(legend)


tiff("Fig5.tiff", height=18, width=17.4, units = "cm", res=600)
ggarrange(p1, p2, p3, p4, legendplot,  nrow=5, heights=c(0.8, 0.8, 0.8, 1, 0.15), labels = c("a", "b", "c", "d"))
dev.off()
