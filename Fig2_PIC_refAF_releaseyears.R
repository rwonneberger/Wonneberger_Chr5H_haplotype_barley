#############################################################################################################################
## Calculcate PIC values and allele frequencies along the genome between groups of cultivars released in different periods ##
#############################################################################################################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(zoo)

colorBlind4   <- c("#ffe119", "#f58231", "#dcbeff", "#4363d8")


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
snps_stat<-transpose_df_Col1(num, "snp")


# PIC and MAF per 20 years

# Make subsets of the population according to release period
metadata$Range="NA"
metadata[metadata$Year <= 1959]$Range <- "1830 - 1959"
metadata[metadata$Year >= 1960 & metadata$Year <= 1979]$Range <- "1960 - 1979"
metadata[metadata$Year >= 1980 & metadata$Year <= 1999]$Range <- "1980 - 1999"
metadata[metadata$Year >= 2000]$Range <- "2000 - 2014"


group_1830 <- metadata[metadata$Range == "1830 - 1959"]$Line
group_1960 <- metadata[metadata$Range == "1960 - 1979"]$Line
group_1980 <- metadata[metadata$Range == "1980 - 1999"]$Line
group_2000 <- metadata[metadata$Range == "2000 - 2014"]$Line


snps_1830<-snps_stat%>%select(snp, all_of(group_1830))
snps_1960<-snps_stat%>%select(snp, all_of(group_1960))
snps_1980<-snps_stat%>%select(snp, all_of(group_1980))
snps_2000<-snps_stat%>%select(snp, all_of(group_2000))


snps_1830_1 <- snps_1830
snps_1960_1 <- snps_1960
snps_1980_1 <- snps_1980
snps_2000_1 <- snps_2000


# Calculate polymorphic information content and allele frequencies
snps_1830_1$AF1<-((rowSums(snps_1830 == 0, na.rm = TRUE)*2)+(rowSums(snps_1830 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_1830[,2:ncol(snps_1830)]))*2)
snps_1960_1$AF1<-((rowSums(snps_1960 == 0, na.rm = TRUE)*2)+(rowSums(snps_1960 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_1960[,2:ncol(snps_1960)]))*2)
snps_1980_1$AF1<-((rowSums(snps_1980 == 0, na.rm = TRUE)*2)+(rowSums(snps_1980 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_1980[,2:ncol(snps_1980)]))*2)
snps_2000_1$AF1<-((rowSums(snps_2000 == 0, na.rm = TRUE)*2)+(rowSums(snps_2000 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_2000[,2:ncol(snps_2000)]))*2)

snps_1830_1$AF2<-((rowSums(snps_1830 == 2, na.rm = TRUE)*2)+(rowSums(snps_1830 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_1830[,2:ncol(snps_1830)]))*2)
snps_1960_1$AF2<-((rowSums(snps_1960 == 2, na.rm = TRUE)*2)+(rowSums(snps_1960 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_1960[,2:ncol(snps_1960)]))*2)
snps_1980_1$AF2<-((rowSums(snps_1980 == 2, na.rm = TRUE)*2)+(rowSums(snps_1980 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_1980[,2:ncol(snps_1980)]))*2)
snps_2000_1$AF2<-((rowSums(snps_2000 == 2, na.rm = TRUE)*2)+(rowSums(snps_2000 == 1, na.rm = TRUE)))/(rowSums(!is.na(snps_2000[,2:ncol(snps_2000)]))*2)

snps_1830_1$MAF<-pmin(snps_1830_1$AF1, snps_1830_1$AF2)
snps_1960_1$MAF<-pmin(snps_1960_1$AF1, snps_1960_1$AF2)
snps_1980_1$MAF<-pmin(snps_1980_1$AF1, snps_1980_1$AF2)
snps_2000_1$MAF<-pmin(snps_2000_1$AF1, snps_2000_1$AF2)

snps_1830_1$PIC<-1-(snps_1830_1$MAF^2+(1-snps_1830_1$MAF)^2)-(2*(snps_1830_1$MAF^2)*(1-snps_1830_1$MAF)^2)
snps_1960_1$PIC<-1-(snps_1960_1$MAF^2+(1-snps_1960_1$MAF)^2)-(2*(snps_1960_1$MAF^2)*(1-snps_1960_1$MAF)^2)
snps_1980_1$PIC<-1-(snps_1980_1$MAF^2+(1-snps_1980_1$MAF)^2)-(2*(snps_1980_1$MAF^2)*(1-snps_1980_1$MAF)^2)
snps_2000_1$PIC<-1-(snps_2000_1$MAF^2+(1-snps_2000_1$MAF)^2)-(2*(snps_2000_1$MAF^2)*(1-snps_2000_1$MAF)^2)

snps_1830_1$Range <- "1830 - 1959"
snps_1960_1$Range <- "1960 - 1979"
snps_1980_1$Range <- "1980 - 1999"
snps_2000_1$Range <- "2000 - 2014"


snps_stats<-rbind(snps_1830_1[, c(1,(ncol(snps_1830_1)-4):ncol(snps_1830_1))], snps_1960_1[, c(1,(ncol(snps_1960_1)-4):ncol(snps_1960_1))], snps_1980_1[, c(1,(ncol(snps_1980_1)-4):ncol(snps_1980_1))], snps_2000_1[, c(1,(ncol(snps_2000_1)-4):ncol(snps_2000_1))])  

snps_stats<-as.data.table(snps_stats)  


# add physical positons
snppos<-hap[,1:4]
names(snppos)[1]<-"snp"


snps_stats_pos<-left_join(snps_stats, snppos, by="snp")

# Calculate PIC and reference AF in rolling windows oa 10000 SNPs
snps_stats_pos$PIC_roll<-rollmean(snps_stats_pos$PIC,10000,fill = list(NA, NULL, NA))

snps_stats_pos$AF_roll<-rollmean(snps_stats_pos$AF1, 10000,fill = list(NA, NULL, NA))



##### Combined plot of PIC and ref AF
# with vertical lines indicating the regions of interest
data_vline_p1<-data.frame(chrom = c("1H", "2H", "3H", "3H", "5H", "5H", "6H", "6H", "7H" ),
                          vline = c(330.42, 81.4, 283.7, 379.4, 68.78, 320.04, 45.6, 382, 410))


scaleFUN <- function(x) sprintf("%.2f", x)

p1<-ggplot(snps_stats_pos, aes(x=pos/1000000, y=PIC_roll, color=Range)) + geom_point(size=0.3) +     theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200), labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) +   scale_colour_manual(values = colorBlind4)+  facet_grid(~chrom, scales='free_x', space="free_x")  + ylab("PIC") + labs(color="Release period")+ theme(legend.position = "none")+ theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+scale_y_continuous(labels=scaleFUN) + geom_vline(data_vline_p1, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red")

p2<-ggplot(snps_stats_pos, aes(x=pos/1000000, y=AF_roll, color=Range)) + geom_point(size=0.3) +    theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200), labels = function(x) format(x, scientific = FALSE))+
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) +  scale_colour_manual(values = colorBlind4)+ facet_grid(~chrom, scales='free_x', space="free_x")+ xlab("Position (Mbp)")  + ylab("Ref. allele freq.") + labs(color="Release period")+ theme(legend.position = "none")+scale_y_continuous(labels=scaleFUN) + geom_vline(data_vline_p1, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red")


p3<-ggplot(snps_stats_pos, aes(x=pos/1000000, y=AF_roll, color=Range)) + geom_point(size=0.3) +    theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200), labels = function(x) format(x, scientific = FALSE))+
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) +  scale_colour_manual(values = colorBlind4)+ facet_grid(~chrom, scales='free_x', space="free_x")+ xlab("Position (Mbp)")  + ylab("Ref. allele freq.") + labs(color="Release period")+theme(legend.position = "bottom")+scale_y_continuous(labels=scaleFUN) + geom_vline(data_vline_p1, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red")+ guides(colour = guide_legend(override.aes = list(size=2)))

legend<-get_legend(p3)
p4<-as_ggplot(legend)

tiff("Fig2.tiff", height=12, width=17.4, units = "cm", res=600)
ggarrange(p1, p2, p4,  nrow=3,  heights = c(0.9, 1, 0.1), labels = c("a", "b"))

dev.off()



