######################################
## Calculate various PCAs and PCoAs ##
######################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")


library(factoextra)
library(FactoMineR)
library(ggplot2)
library(broom)
library(adegenet)
library(RColorBrewer)
library(ggfortify)
library(ggrepel)
library(ggpubr)



colorBlindGrey8   <- c( "#D55E00", "#E69F00","#F0E442","#009E73", 
                        "#0072B2","#56B4E9","#CC79A7", "#999999")

#########################
## PCoA using all SNPs ##
#########################

# Read in SNP matrix published in Schreiber et al., 2023 as hapmap file. The genotypes need to be encoded as A, C, G, T
hmp<-fread("Hapmap_haploid.hmp.txt")

# Read in the SNP information in numeric transposed format
num<-fread("SNP_numeric.txt")
num<-as.data.frame(num)

# Read in a metadata file generated from Online Resource 1. 
metadata <- read.table("Metadata.txt", head=TRUE, sep="\t")

# Make sure that the accession names in the snpmatrix and the metadata are identical
SamplesKeep <- colnames(hmp)[-c(1:11)]
metadata_reduced <- metadata[metadata$Line%in%SamplesKeep,]

metadata_reduced<-metadata_reduced[match(names(hmp)[-c(1:11)], metadata_reduced$Line),]


table(names(hmp)[-c(1:11)]  == metadata_reduced$Line)


# Make a distance matrix from the numeric SNP data
geno_dist <- dist(num)

#Perform multidimensional scaling (mds)
geno_dist_ind <- cmdscale(geno_dist,k=3,eig=TRUE)
mds_geno_ind_df <- data.frame(geno_dist_ind$points,check.names = TRUE)
colnames(mds_geno_ind_df) <- c("PCo1","PCo2", "PCo3")

# Add passport data to the df
mds_geno_ind_df$Genotype<-metadata_reduced$Line
mds_geno_ind_df$`Release period` <- metadata_reduced$Range
mds_geno_ind_df$Country <- metadata_reduced$Country
mds_geno_ind_df[mds_geno_ind_df$Country == "" , ]$Country<-"Unknown"

# Calculate % variation explained by the PCos
mds.var.per<-round(geno_dist_ind$eig/sum(geno_dist_ind$eig)*100, 1)



# Plot colored by release period
p1<-ggplot(mds_geno_ind_df, aes(x=PCo1,y=PCo2, color = `Release period`)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", mds.var.per[1], "%)") , y = paste0("PC2 (", mds.var.per[2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colorBlindGrey8)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7)) +   theme(legend.position = "none")


p2<-ggplot(mds_geno_ind_df, aes(x=PCo1,y=PCo3, color = `Release period`)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", mds.var.per[1], "%)") , y = paste0("PC3 (", mds.var.per[3], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colorBlindGrey8)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7)) +   theme(legend.position = "none")

p3<-ggplot(mds_geno_ind_df, aes(x=PCo2,y=PCo3, color = `Release period`)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC2 (", mds.var.per[2], "%)") , y = paste0("PC3 (", mds.var.per[3], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colorBlindGrey8)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7)) +   theme(legend.position = "none")

#make a new plot to be able to extract the legend from it9
p4<-ggplot(mds_geno_ind_df, aes(x=PCo2,y=PCo3, color = `Release period`)) +
  geom_point(size=0.5) + theme_classic()+
  labs(x = paste0("PC2 (", mds.var.per[2], "%)") , y = paste0("PC3 (", mds.var.per[3], "%)"))+
  scale_colour_manual(values = colorBlindGrey8) +  
  theme(legend.title = element_text(size=7),
        legend.text = element_text(size=6), 
        legend.spacing.x = unit(0, 'cm'))


legend<-get_legend(p4)
p5<-as_ggplot(legend)

tiff("Fig3.tiff", unit="cm", res=600, width=17.4, height=6)
figure<-ggarrange(p1, p2, p3, p5, ncol=4, widths=c(1,1,1,0.4))

annotate_figure(figure,
                top = text_grob("", size = 14))


dev.off()



#######################################
## Plot colored by country of origin ##
#######################################

countrycolor<-c('#9A6324', '#3cb44b', '#ffe119', '#4363d8', '#42d4f4','#f58231', '#911eb4',  '#f032e6', '#bfef45', '#fabed4', '#469990',  '#e6194B', '#a9a9a9', '#ffffff', '#000000')


p1<-ggplot(mds_geno_ind_df, aes(x=PCo1,y=PCo2, color = Country)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", mds.var.per[1], "%)") , y = paste0("PC2 (", mds.var.per[2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = countrycolor)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7)) +   theme(legend.position = "none")

p2<-ggplot(mds_geno_ind_df, aes(x=PCo1,y=PCo3, color = Country)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", mds.var.per[1], "%)") , y = paste0("PC3 (", mds.var.per[3], "%)"))+
  theme_classic()+
  scale_colour_manual(values = countrycolor)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7)) +   theme(legend.position = "none")

p3<-ggplot(mds_geno_ind_df, aes(x=PCo2,y=PCo3, color = Country)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC2 (", mds.var.per[2], "%)") , y = paste0("PC3 (", mds.var.per[3], "%)"))+
  theme_classic()+
  scale_colour_manual(values = countrycolor)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7)) +   theme(legend.position = "none")

#make a new plot to be able to extract the legend from it9
p4<-ggplot(mds_geno_ind_df, aes(x=PCo2,y=PCo3, color = Country)) +
  geom_point(size=0.5) + theme_classic()+
  labs(x = paste0("PC2 (", mds.var.per[2], "%)") , y = paste0("PC3 (", mds.var.per[3], "%)"))+
  scale_colour_manual(values = countrycolor) +
  theme(legend.title = element_text(size=7),
        legend.text = element_text(size=6)) +
  theme(legend.spacing.y = unit(-0.2, 'cm')) +
  guides(color = guide_legend(byrow = TRUE))


legend<-get_legend(p4)
p5<-as_ggplot(legend)

tiff("ESM5.tiff", unit="cm", res=600, width=17.4, height=6)
figure<-ggarrange(p1, p2, p3, p5, ncol=4, widths=c(1,1,1,0.45))

annotate_figure(figure,
                top = text_grob("", size = 14))


dev.off()




################################################################
## PCoA color-coded by cluster membership at DAPC k=2 and k=3 ##
################################################################

# Color by DAPC k=2


mds_geno_ind_df$Cluster<-as.factor(metadata_reduced$Cluster_k2)

colorBlindGrey7   <- c( "#E69F00", "#56B4E9", "#009E73", 
                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colorscale   <- c( "#800000",  "#4363d8", "#ffe119")

p1<-ggplot(mds_geno_ind_df, aes(x=PCo1,y=PCo2, color = Cluster, label=Genotype)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", mds.var.per[1], "%)") , y = paste0("PC2 (", mds.var.per[2], "%)"))+
  theme_classic()+
  
  scale_colour_manual(values = colorscale)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7))  + theme(legend.position = "none")


# Color by DAPC k=3

mds_geno_ind_df$Cluster<-as.factor(metadata_reduced$Cluster_k3)

p2<-ggplot(mds_geno_ind_df, aes(x=PCo1,y=PCo2, color = Cluster, label=Genotype)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", mds.var.per[1], "%)") , y = paste0("PC2 (", mds.var.per[2], "%)"))+
  theme_classic()+
  
  scale_colour_manual(values = colorscale)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7))  + theme(legend.position = "none")

p3<-ggplot(mds_geno_ind_df, aes(x=PCo1,y=PCo2, color = Cluster, label=Genotype)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", mds.var.per[1], "%)") , y = paste0("PC2 (", mds.var.per[2], "%)"))+
  theme_classic()+
  
  scale_colour_manual(values = colorscale)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7))  + labs(color="k-means cluster")+
  theme(legend.title = element_text(size=7),
        legend.text = element_text(size=6))


legend<-get_legend(p3)
legendplot<-as_ggplot(legend)


tiff("ESM7.tiff",  unit="cm", res=600, width=17.4, height=8)
ggarrange(p1, p2, legendplot,ncol=3, widths=c(1,1,0.3),  labels = c("a", "b")) 
dev.off()



#######################################################################################################################
## PCoA on the 5H SNPs in the region identified by FST between cluster 1 and cluster 2 at k=2 (68.78 and 320.04 Mbp) ##
#######################################################################################################################

start<-68780000
end<-320040000


# Read in SNP matrix published in Schreiber et al., 2023 as hapmap file. The genotypes need to be encoded as A, C, G, T

hmp<-fread("Hapmap_haploid.hmp.txt")


# Read in the SNP information in numeric transposed format
num<-fread("SNP_numeric.txt")
num<-as.data.frame(num)

# Read in a metadata file generated from Online Resource 1. 
metadata <- read.table("Metadata.txt", head=TRUE, sep="\t")


SamplesKeep <- colnames(hmp)[-c(1:11)]
metadata_reduced <- metadata[metadata$Line%in%SamplesKeep,]

#Sort metadata_reduced$Line by the order of accessions in the Hap file
metadata_reduced<-metadata_reduced[match(names(hmp)[-c(1:11)], metadata_reduced$Line),]


table(metadata_reduced$Line %in% names(hmp)[-c(1:11)])
table(names(hmp)[-c(1:11)] %in% metadata_reduced$Line) 

table(names(hmp)[-c(1:11)]  == metadata_reduced$Line)

# Select ony the SNPs in the region of interest
hmp_sub<-hmp%>%filter(chrom == "5H", pos >= start, pos <=end)

num_sub<-num[, hmp_sub$`rs#`]


geno_dist <- dist(num_sub)

#Perform multidimensional scaling (mds)
geno_dist_ind <- cmdscale(geno_dist,k=4,eig=TRUE)
mds_geno_ind_df <- data.frame(geno_dist_ind$points,check.names = TRUE)
colnames(mds_geno_ind_df) <- c("PCo1","PCo2", "PCo3", "PCo4")


# Add passport data
mds_geno_ind_df$`Release period` <- metadata_reduced$Range
mds_geno_ind_df$Genotype<-metadata_reduced$Line
mds_geno_ind_df$Cluster5H<-metadata_reduced$Haplotype

mds.var.per<-round(geno_dist_ind$eig/sum(geno_dist_ind$eig)*100, 1)



tiff("ESM15.tiff", unit="cm", res=600, width=8.4, height=7)
ggplot(mds_geno_ind_df, aes(x=PCo1,y=PCo2, color = Cluster5H)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", mds.var.per[1], "%)") , y = paste0("PC2 (", mds.var.per[2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = haplo_colors)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7))+labs(color="Chr 5H haplotype")+
  theme(legend.title = element_text(size=7),
        legend.text = element_text(size=6))
dev.off()
