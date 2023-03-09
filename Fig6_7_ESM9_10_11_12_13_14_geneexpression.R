###############################################################
## k-means clustering, PCoAs and heatmaps of gene expression ##
###############################################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library(ggplotify)
library(pheatmap)
library(factoextra)
library(FactoMineR)
library(broom)
library(adegenet)
library(RColorBrewer)
library(ggfortify)
library(ggrepel)
library(ggpubr)
library(cowplot)



colorBlindGrey8   <- c("#999999", "#D55E00", "#E69F00","#F0E442","#009E73", "#0072B2","#56B4E9","#CC79A7")
colorBlindGrey5   <- c("#D55E00" , "#E69F00","#56B4E9","#009E73",  "#999999")
haplo_colors=c('#E69F00', '#009E73', '#56B4E9', "#D55E00")
colors=c('#e6194B','#ffe119', '#3cb44b',  '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075' )


# Format gene expression file to the correct format for plotting a PCoA on cultivars
make_pca_5Hgenes_ind<-function(df){
  expmat_t<-transpose_df(as.data.frame(df))
  
  geno_dist_ind <- dist(expmat_t)
  
  #Perform multidimensional scaling (mds)
  mds_geno_ind <- cmdscale(geno_dist_ind,k=2,eig=TRUE)
  mds_geno_ind_df <- data.frame(mds_geno_ind$points,check.names = TRUE)
  colnames(mds_geno_ind_df) <- c("PCo1","PCo2")
  
  SamplesKeep <- names(df)[-1]
  metadata_reduced <- metadata[metadata$Line%in%SamplesKeep,]
  
  metadata_reduced<-metadata_reduced[match(names(df)[-1], metadata_reduced$Line),]
  
  print(table(metadata_reduced$Line %in% names(df)[-1]))
  print(table(names(df)[-1] %in% metadata_reduced$Line))
  print(table(names(df)[-1]  == metadata_reduced$Line))
  
  mds_geno_ind_df$`Release period` <- metadata_reduced$Range
  mds_geno_ind_df$Name<-metadata_reduced$Line
  mds_geno_ind_df$Haplotype<-metadata_reduced$Haplotype

  
  mds.var.per<-round(mds_geno_ind$eig/sum(mds_geno_ind$eig)*100, 1)
  
  output<-list(mds_geno_ind_df, mds.var.per)
  
  return(output)
}

# Format gene expression file to the correct format for plotting a PCoA on genes
make_pca_5Hgenes_gene<-function(df, goi){
  
  expmat_t<-transpose_df(as.data.frame(df))
  names(df)[1]<-"Gene"
  mds_geno_ind <- dist(df[, -1])
  
  #Perform multidimensional scaling (mds)
  mds_geno_ind <- cmdscale(mds_geno_ind,k=4,eig=TRUE)
  mds_geno_ind_df <- data.frame(mds_geno_ind$points,check.names = TRUE)
  colnames(mds_geno_ind_df) <- c("PCo1","PCo2", "PCo3", "PCo4")
  
  df$avg<-rowMeans(df[,-1])
  mds_geno_ind_df$Gene<-df$Gene
  
  
  mds_geno_ind_df$geneofinterest<-mds_geno_ind_df$Gene%in%goi
  mds_geno_ind_df<-left_join(mds_geno_ind_df, df, by="Gene")
  
  mds.var.per<-round(mds_geno_ind$eig/sum(mds_geno_ind$eig)*100, 1)
  
  output<-list(mds_geno_ind_df, mds.var.per)
  
  return(output)
}


##############################

start<-68780000
end<-320040000


# Read in SNP matrix published in Schreiber et al., 2023 as hapmap file. The genotypes need to be encoded as A, C, G, T

hmp<-fread("Hapmap_haploid.hmp.txt")

# Read in the SNP information in numeric transposed format
num<-fread("SNP_numeric.txt")
num<-as.data.frame(num)


# Read in a metadata file generated from Online Resource 1 and make sure that hmp and Metadata contain the same accessions
metadata <- read.table("Metadata.txt", head=TRUE, sep="\t")

SamplesKeep <- colnames(hmp)[-c(1:11)]
metadata_reduced <- metadata[metadata$Line%in%SamplesKeep,]

#Sort metadata_reduced$Line by the order of accessions in the Hap file
metadata_reduced<-metadata_reduced[match(names(hmp)[-c(1:11)], metadata_reduced$Line),]

table(metadata_reduced$Line %in% names(hmp)[-c(1:11)])
table(names(hmp)[-c(1:11)] %in% metadata_reduced$Line) 

table(names(hmp)[-c(1:11)]  == metadata_reduced$Line)

#######################################################################################################################################################
## PCoA on the cultivars using expression of the genes in the 5H region identified by FST between k1 and k3 (68.78 - 320.04 Mbp) and flanking 10 kbp ##
#######################################################################################################################################################

geneloc<-fread("genelocation_BaRT2v18.txt") # Read in a file containing the start and end positions and chromosome locations of each BaRT2V18 gene
names(geneloc)<-c("Gene", "Chr", "Start", "End")

geneloc_filt<-geneloc%>%filter(Chr == "5H", Start >= start-10000, End<= end+10000)

# The following files containing the log cpm read data for each tissue can be found in Scghreiber et al., 2023
# Read in cpm values - note: For this analysis we filtered the dataset to keep only genes with > 1 cpm in > 30% of accessions
crown_log_cpm<-fread("Log_cpm.Crown.0.3.1.txt")
root_log_cpm<-fread("Log_cpm.Root.0.3.1.txt")
devinf_log_cpm<-fread("Log_cpm.DevInf.0.3.1.txt")
ped_log_cpm<-fread("Log_cpm.Ped.0.3.1.txt")
ga_log_cpm<-fread("Log_cpm.GA.0.3.1.txt")
pa_log_cpm<-fread("Log_cpm.PA.0.3.1.txt")

# Keep only the genes in the region of interest
crown_sub<-crown_log_cpm[crown_log_cpm$V1%in%geneloc_filt$Gene,]
root_sub<-root_log_cpm[root_log_cpm$V1%in%geneloc_filt$Gene,]
devinf_sub<-devinf_log_cpm[devinf_log_cpm$V1%in%geneloc_filt$Gene,]
ped_sub<-ped_log_cpm[ped_log_cpm$V1%in%geneloc_filt$Gene,]
ga_sub<-ga_log_cpm[ga_log_cpm$V1%in%geneloc_filt$Gene,]
pa_sub<-pa_log_cpm[pa_log_cpm$V1%in%geneloc_filt$Gene,]


# Make a distance matrix and do multidimensional scaling on the gene expression data - per individual

crown_ind<-make_pca_5Hgenes_ind(crown_sub)
root_ind<-make_pca_5Hgenes_ind(root_sub)
devinf_ind<-make_pca_5Hgenes_ind(devinf_sub)
ped_ind<-make_pca_5Hgenes_ind(ped_sub)
ga_ind<-make_pca_5Hgenes_ind(ga_sub)
pa_ind<-make_pca_5Hgenes_ind(pa_sub)


crownplot<-ggplot(crown_ind[[1]], aes(x=PCo1,y=PCo2, color = Haplotype )) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", crown_ind[[2]][1], "%)") , y = paste0("PC2 (", crown_ind[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = haplo_colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8))+theme(legend.position = "none") 


rootplot<-ggplot(root_ind[[1]], aes(x=PCo1,y=PCo2, color = Haplotype)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", root_ind[[2]][1], "%)") , y = paste0("PC2 (", root_ind[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = haplo_colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8))+theme(legend.position = "none") 

devinfplot<-ggplot(devinf_ind[[1]], aes(x=PCo1,y=PCo2, color = Haplotype)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", devinf_ind[[2]][1], "%)") , y = paste0("PC2 (", devinf_ind[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = haplo_colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8))+theme(legend.position = "none") 

pedplot<-ggplot(ped_ind[[1]], aes(x=PCo1,y=PCo2, color = Haplotype)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", ped_ind[[2]][1], "%)") , y = paste0("PC2 (", ped_ind[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = haplo_colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8))+theme(legend.position = "none") 

gaplot<-ggplot(ga_ind[[1]], aes(x=PCo1,y=PCo2, color = Haplotype)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", ga_ind[[2]][1], "%)") , y = paste0("PC2 (", ga_ind[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = haplo_colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8))+theme(legend.position = "none") 


paplot<-ggplot(pa_ind[[1]], aes(x=PCo1,y=PCo2, color = Haplotype)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", pa_ind[[2]][1], "%)") , y = paste0("PC2 (", ga_ind[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = haplo_colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8))+theme(legend.position = "none") 


legend<-ggplot(pa_ind[[1]], aes(x=PCo1,y=PCo2, color = Haplotype)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", pa_ind[[2]][1], "%)") , y = paste0("PC2 (", ga_ind[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = haplo_colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8))+theme(legend.position="bottom", legend.box = "horizontal")+ labs(color = "Chr 5H haplotype")+
  guides(color=guide_legend( title.position="top"))+
  theme(legend.title = element_text(size=8),
        legend.text = element_text(size=7)) +
  theme(legend.spacing.x = unit(0.01, 'cm')) 


legend<-get_legend(legend)
legendplot<-as_ggplot(legend)

figure1<-annotate_figure(crownplot, top = text_grob("Crown", size = 9))
figure2<-annotate_figure(rootplot, top = text_grob("Root", size = 9))
figure3<-annotate_figure(devinfplot, top = text_grob("Dev. inflorescence", size = 9))
figure4<-annotate_figure(pedplot, top = text_grob("Peduncle", size = 9))
figure5<-annotate_figure(gaplot, top = text_grob("Spikelet", size = 9))
figure6<-annotate_figure(paplot, top = text_grob("Developing grain", size = 9))

tiff("Fig_6.tiff", res=600, units="cm", height=14, width=8.4)

ggdraw()+
  draw_plot(figure1, x = 0, y = 0.7, width = 0.5, height = 0.3)+
  draw_plot(figure2, x = 0.5, y = 0.7, width = 0.5, height = 0.3)+
  draw_plot(figure3, x = 0, y = 0.4, width = 0.5, height = 0.3) +
  draw_plot(figure4, x = 0.5, y = 0.4, width = 0.5, height = 0.3)+
  draw_plot(figure5, x = 0, y = 0.1, width = 0.5, height = 0.3)+
  draw_plot(figure6, x = 0.5, y = 0.1, width = 0.5, height = 0.3)+
  draw_plot(legend, x = 0, y = 0, width = 1, height = 0.1)+
  draw_plot_label(label = c("a", "b", "c", "d", "e", "f"), size = 9,
                  x = c(0, 0.5, 0, 0.5, 0, 0.5), y = c(1, 1, 0.7, 0.7,  0.4, 0.4 ))

dev.off()


###################################################################################################################################################
## PCoA on the genes using expression of the genes in the 5H region identified by FST between k1 and k3 (68.78 - 320.04 Mbp) and flanking 10 kbp ##
###################################################################################################################################################

# Do a k-means clustering on the gene expression

crown_sub1<-crown_sub[,-1]
root_sub1<-root_sub[,-1]
devinf_sub1<-devinf_sub[,-1]
ped_sub1<-ped_sub[,-1]
ga_sub1<-ga_sub[,-1]
pa_sub1<-pa_sub[,-1]


# You would run the outcommented part below on different numbers of potential clusters (k) and overlay the results on a PCoA to see if the results make sense. We did this and concluded that k=4 is suitable to identify the cluster of genes we are interested in

# grp_10_4<-find.clusters(crown_sub1, max.n.clust=50,  n.pca=10, n.clust=4)
# 
# crown_grpout<-as.data.frame(grp_10_4$grp)
# names(crown_grpout)[1]<-"clusters"
# crown_grpout$accession<-rownames(crown_grpout)
# crown_grpout<-cbind(crown_grpout, crown_sub)
# fwrite(as.data.frame(crown_grpout), "//filer/projekte/barn/PCAs_5HGenes/IPK_crown_grp_10_4.txt", sep="\t", row.names=F)
# 
# 
# 
# grp_10_4<-find.clusters(root_sub1, max.n.clust=50,  n.pca=10, n.clust=4)
# 
# root_grpout<-as.data.frame(grp_10_4$grp)
# names(root_grpout)[1]<-"clusters"
# root_grpout$accession<-rownames(root_grpout)
# root_grpout<-cbind(root_grpout, root_sub)
# fwrite(as.data.frame(root_grpout), "//filer/projekte/barn/PCAs_5HGenes/IPK_root_grp_10_4.txt", sep="\t", row.names=F)
# 
# 
# grp_10_4<-find.clusters(devinf_sub1, max.n.clust=50,  n.pca=10, n.clust=4)
# 
# devinf_grpout<-as.data.frame(grp_10_4$grp)
# names(devinf_grpout)[1]<-"clusters"
# devinf_grpout$accession<-rownames(devinf_grpout)
# devinf_grpout<-cbind(devinf_grpout, devinf_sub)
# fwrite(as.data.frame(devinf_grpout), "//filer/projekte/barn/PCAs_5HGenes/IPK_devinf_grp_10_4.txt", sep="\t", row.names=F)
# 
# 
# grp_10_4<-find.clusters(ped_sub1, max.n.clust=50,  n.pca=10, n.clust=4)
# 
# ped_grpout<-as.data.frame(grp_10_4$grp)
# names(ped_grpout)[1]<-"clusters"
# ped_grpout$accession<-rownames(ped_grpout)
# ped_grpout<-cbind(ped_grpout, ped_sub)
# fwrite(as.data.frame(ped_grpout), "//filer/projekte/barn/PCAs_5HGenes/IPK_ped_grp_10_4.txt", sep="\t", row.names=F)
# 
# 
# grp_10_4<-find.clusters(ga_sub1, max.n.clust=50,  n.pca=10, n.clust=4)
# 
# ga_grpout<-as.data.frame(grp_10_4$grp)
# names(ga_grpout)[1]<-"clusters"
# ga_grpout$accession<-rownames(ga_grpout)
# ga_grpout<-cbind(ga_grpout, ga_sub)
# fwrite(as.data.frame(ga_grpout), "//filer/projekte/barn/PCAs_5HGenes/IPK_ga_grp_10_4.txt", sep="\t", row.names=F)
# 
# 
# grp_10_4<-find.clusters(pa_sub1, max.n.clust=50,  n.pca=10, n.clust=4)
# 
# pa_grpout<-as.data.frame(grp_10_4$grp)
# names(pa_grpout)[1]<-"clusters"
# pa_grpout$accession<-rownames(pa_grpout)
# pa_grpout<-cbind(pa_grpout, pa_sub)
# fwrite(as.data.frame(pa_grpout), "//filer/projekte/barn/PCAs_5HGenes/IPK_pa_grp_10_4.txt", sep="\t", row.names=F)


crown_grpout<-fread("//filer/projekte/barn/PCAs_5HGenes/IPK_crown_grp_10_4_namesfixed.txt", sep="\t")
root_grpout<-fread("//filer/projekte/barn/PCAs_5HGenes/IPK_root_grp_10_4_namesfixed.txt", sep="\t")
devinf_grpout<-fread("//filer/projekte/barn/PCAs_5HGenes/IPK_devinf_grp_10_4_namesfixed.txt", sep="\t")
ped_grpout<-fread("//filer/projekte/barn/PCAs_5HGenes/IPK_ped_grp_10_4_namesfixed.txt", sep="\t")
ga_grpout<-fread("//filer/projekte/barn/PCAs_5HGenes/IPK_ga_grp_10_4_namesfixed.txt", sep="\t")
pa_grpout<-fread("//filer/projekte/barn/PCAs_5HGenes/IPK_pa_grp_10_4_namesfixed.txt", sep="\t")


crown_grpout$clusters<-as.character(crown_grpout$clusters)
root_grpout$clusters<-as.character(root_grpout$clusters)
devinf_grpout$clusters<-as.character(devinf_grpout$clusters)
ped_grpout$clusters<-as.character(ped_grpout$clusters)
ga_grpout$clusters<-as.character(ga_grpout$clusters)
pa_grpout$clusters<-as.character(pa_grpout$clusters)

# The cluster of interest will have a different number in the different tissues. To allow the use of the same color for this cluster in each tissue, we rename the cluster numbers below
crown_1 <- which(crown_grpout$clusters==1)
crown_2 <- which(crown_grpout$clusters==2)
crown_3 <- which(crown_grpout$clusters==3)
crown_4 <- which(crown_grpout$clusters==4)

root_1 <- which(root_grpout$clusters==1)
root_2 <- which(root_grpout$clusters==2)
root_3 <- which(root_grpout$clusters==3)
root_4 <- which(root_grpout$clusters==4)

devinf_1 <- which(devinf_grpout$clusters==1)
devinf_2 <- which(devinf_grpout$clusters==2)
devinf_3 <- which(devinf_grpout$clusters==3)
devinf_4 <- which(devinf_grpout$clusters==4)

ped_1 <- which(ped_grpout$clusters==1)
ped_2 <- which(ped_grpout$clusters==2)
ped_3 <- which(ped_grpout$clusters==3)
ped_4 <- which(ped_grpout$clusters==4)

ga_1 <- which(ga_grpout$clusters==1)
ga_2 <- which(ga_grpout$clusters==2)
ga_3 <- which(ga_grpout$clusters==3)
ga_4 <- which(ga_grpout$clusters==4)

pa_1 <- which(pa_grpout$clusters==1)
pa_2 <- which(pa_grpout$clusters==2)
pa_3 <- which(pa_grpout$clusters==3)
pa_4 <- which(pa_grpout$clusters==4)


crown_grpout$clusters[crown_1] <-4
crown_grpout$clusters[crown_2] <-3
crown_grpout$clusters[crown_3] <-2
crown_grpout$clusters[crown_4] <-1

root_grpout$clusters[root_1] <-3
root_grpout$clusters[root_2] <-2
root_grpout$clusters[root_3] <-4
root_grpout$clusters[root_4] <-1

devinf_grpout$clusters[devinf_1] <-4
devinf_grpout$clusters[devinf_2] <-1
devinf_grpout$clusters[devinf_3] <-3
devinf_grpout$clusters[devinf_4] <-2

ped_grpout$clusters[ped_1] <-1
ped_grpout$clusters[ped_2] <-2
ped_grpout$clusters[ped_3] <-4
ped_grpout$clusters[ped_4] <-3

ga_grpout$clusters[ga_1] <-3
ga_grpout$clusters[ga_2] <-2
ga_grpout$clusters[ga_3] <-1
ga_grpout$clusters[ga_4] <-4

pa_grpout$clusters[pa_1] <-4
pa_grpout$clusters[pa_2] <-1
pa_grpout$clusters[pa_3] <-3
pa_grpout$clusters[pa_4] <-2


crown_grpout$clusters<-as.factor(crown_grpout$clusters)
root_grpout$clusters<-as.factor(root_grpout$clusters)
devinf_grpout$clusters<-as.factor(devinf_grpout$clusters)
ped_grpout$clusters<-as.factor(ped_grpout$clusters)
ga_grpout$clusters<-as.factor(ga_grpout$clusters)
pa_grpout$clusters<-as.factor(pa_grpout$clusters)



# Make a distance matrix and do multidimensional scaling on the gene expression data - per gene
crown_mds_geno<-make_pca_5Hgenes_gene(crown_sub, geneloc_filt$Gene)
root_mds_geno<-make_pca_5Hgenes_gene(root_sub, geneloc_filt$Gene)
devinf_mds_geno<-make_pca_5Hgenes_gene(devinf_sub, geneloc_filt$Gene)
ped_mds_geno<-make_pca_5Hgenes_gene(ped_sub, geneloc_filt$Gene)
ga_mds_geno<-make_pca_5Hgenes_gene(ga_sub, geneloc_filt$Gene)
pa_mds_geno<-make_pca_5Hgenes_gene(pa_sub, geneloc_filt$Gene)


crown_1_2<-ggplot(crown_mds_geno[[1]], aes(x=PCo1,y=PCo2,  color=crown_grpout$clusters)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", crown_mds_geno[[2]][1], "%)") , y = paste0("PC2 (", crown_mds_geno[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8)) +theme(legend.position = "none") 


root_1_2<-ggplot(root_mds_geno[[1]], aes(x=PCo1,y=PCo2,  color=root_grpout$clusters)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", root_mds_geno[[2]][1], "%)") , y = paste0("PC2 (", root_mds_geno[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8)) +theme(legend.position = "none") 


devinf_1_2<-ggplot(devinf_mds_geno[[1]], aes(x=PCo1,y=PCo2,  color=devinf_grpout$clusters)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", devinf_mds_geno[[2]][1], "%)") , y = paste0("PC2 (", devinf_mds_geno[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8)) +theme(legend.position = "none") 


ped_1_2<-ggplot(ped_mds_geno[[1]], aes(x=PCo1,y=PCo2,  color=ped_grpout$clusters)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", ped_mds_geno[[2]][1], "%)") , y = paste0("PC2 (", ped_mds_geno[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8)) +theme(legend.position = "none") 


ga_1_2<-ggplot(ga_mds_geno[[1]], aes(x=PCo1,y=PCo2,  color=ga_grpout$clusters)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", ga_mds_geno[[2]][1], "%)") , y = paste0("PC2 (", ga_mds_geno[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8)) +theme(legend.position = "none") 


pa_1_2<-ggplot(pa_mds_geno[[1]], aes(x=PCo1,y=PCo2,  color=pa_grpout$clusters)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", pa_mds_geno[[2]][1], "%)") , y = paste0("PC2 (", pa_mds_geno[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colors)+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8)) +theme(legend.position="none") + labs(color = "k-means cluster ")



legend<-ggplot(pa_mds_geno[[1]], aes(x=PCo1,y=PCo2,  color=pa_grpout$clusters)) +
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", pa_mds_geno[[2]][1], "%)") , y = paste0("PC2 (", pa_mds_geno[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colors)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7))  +theme(legend.position="bottom") + labs(color = "k-means cluster")  + theme(legend.title = element_text(size=8), legend.text = element_text(size=7)) 

legend<-get_legend(legend)
legendplot<-as_ggplot(legend)

figure1<-annotate_figure(crown_1_2, top = text_grob("Crown", size = 9))
figure2<-annotate_figure(root_1_2, top = text_grob("Root", size = 9))
figure3<-annotate_figure(devinf_1_2, top = text_grob("Dev. inflorescence", size = 9))
figure4<-annotate_figure(ped_1_2, top = text_grob("Peduncle", size = 9))
figure5<-annotate_figure(ga_1_2, top = text_grob("Spikelet", size = 9))
figure6<-annotate_figure(pa_1_2, top = text_grob("Developing grain", size = 9))

tiff("Fig_7.tiff", res=600, units="cm", height=14, width=8.4)

ggdraw()+
  draw_plot(figure1, x = 0, y = 0.7, width = 0.5, height = 0.3)+
  draw_plot(figure2, x = 0.5, y = 0.7, width = 0.5, height = 0.3)+
  draw_plot(figure3, x = 0, y = 0.4, width = 0.5, height = 0.3) +
  draw_plot(figure4, x = 0.5, y = 0.4, width = 0.5, height = 0.3)+
  draw_plot(figure5, x = 0, y = 0.1, width = 0.5, height = 0.3)+
  draw_plot(figure6, x = 0.5, y = 0.1, width = 0.5, height = 0.3)+
  draw_plot(legend, x = 0, y = 0, width = 1, height = 0.1)+
  draw_plot_label(label = c("a", "b", "c", "d", "e", "f"), size = 9,
                  x = c(0, 0.5, 0, 0.5, 0, 0.5), y = c(1, 1, 0.7, 0.7,  0.4, 0.4 ))

dev.off()



########################################################################
## Heatmaps of the expression of the genes of interest in each tissue ##
########################################################################

## Manually select the genes in the cluster of interest - genes that are upregulated in the new haplotype
crown_goi<-crown_grpout[crown_grpout$clusters == 1, ]
root_goi<-root_grpout[root_grpout$clusters == 1, ]
devinf_goi<-devinf_grpout[devinf_grpout$clusters == 1, ]
ped_goi<-ped_grpout[ped_grpout$clusters == 1, ]
ga_goi<-ga_grpout[ga_grpout$clusters == 1, ]
pa_goi<-pa_grpout[pa_grpout$clusters == 1, ]


all_goi_up<-c(crown_goi$V1, root_goi$V1, devinf_goi$V1, ped_goi$V1,ga_goi$V1, pa_goi$V1 )
all_goi_up<-unique(all_goi_up)


## Manually select the genes in the cluster of interest - extract the two downregulated genes
crown_down<-crown_mds_geno[[1]]%>%filter(PCo2 < -30)%>%select(Gene)
root_down<-root_mds_geno[[1]]%>%filter(PCo2 < -50)%>%select(Gene)
devinf_down<-devinf_mds_geno[[1]]%>%filter(PCo2 > 30)%>%select(Gene)
ped_down<-ped_mds_geno[[1]]%>%filter(PCo2 < -40)%>%select(Gene)
ga_down<-ga_mds_geno[[1]]%>%filter(PCo2 < -35)%>%select(Gene)
pa_down<-pa_mds_geno[[1]]%>%filter(PCo2 < -23)%>%select(Gene)

all_goi_down<-c(crown_down$Gene, root_down$Gene, devinf_down$Gene, ped_down$Gene,ga_down$Gene, pa_down$Gene )
all_goi_down<-unique(all_goi_down)

# Combine up- and downregulated genes
all_goi<-c(all_goi_up, all_goi_down)


# Heatmaps showing the expression of the genes of interest in each tissue
names(metadata)[1]<-"Genotype"

names(crown_log_cpm)[1]<-"Genotype"
names(root_log_cpm)[1]<-"Genotype"
names(ga_log_cpm)[1]<-"Genotype"
names(pa_log_cpm)[1]<-"Genotype"
names(devinf_log_cpm)[1]<-"Genotype"
names(ped_log_cpm)[1]<-"Genotype"

crown_log_cpm$Tissue <- "Crown"
root_log_cpm$Tissue <- "Root"
ga_log_cpm$Tissue <- "Spikelet"
pa_log_cpm$Tissue <- "Grain"
devinf_log_cpm$Tissue <- "Inflorescence"
ped_log_cpm$Tissue <- "Peduncle"



allcpm<-plyr::rbind.fill(crown_log_cpm, root_log_cpm, ga_log_cpm, pa_log_cpm, devinf_log_cpm, ped_log_cpm)
dim(allcpm)
allcpm<-allcpm%>%dplyr::select(order(colnames(allcpm)))
allcpm%<>%relocate(Genotype, Tissue)
names(allcpm)[1]<-"Gene"
length(unique(allcpm$Gene))

allcpm$Tissue<-factor(allcpm$Tissue, levels = c("Crown", "Inflorescence", "Peduncle", "Spikelet", "Grain", "Root"))


genelist<-all_goi

Crown_sub<-crown_log_cpm[crown_log_cpm$Genotype%in%genelist, ]
Root_sub<-root_log_cpm[root_log_cpm$Genotype%in%genelist, ]
GA_sub<-ga_log_cpm[ga_log_cpm$Genotype%in%genelist, ]
PA_sub<-pa_log_cpm[pa_log_cpm$Genotype%in%genelist, ]
DevInf_sub<-devinf_log_cpm[devinf_log_cpm$Genotype%in%genelist, ]
Ped_sub<-ped_log_cpm[ped_log_cpm$Genotype%in%genelist, ]

Crown_sub$Tissue <- NULL
Root_sub$Tissue <- NULL
GA_sub$Tissue <- NULL
PA_sub$Tissue <- NULL
DevInf_sub$Tissue <- NULL
Ped_sub$Tissue <- NULL

Crown_sub_t<-transpose_df(Crown_sub)
Root_sub_t<-transpose_df(Root_sub)
GA_sub_t<-transpose_df(GA_sub)
PA_sub_t<-transpose_df(PA_sub)
DevInf_sub_t<-transpose_df(DevInf_sub)
Ped_sub_t<-transpose_df(Ped_sub)

Crown_sub_t$Genotype<-rownames(Crown_sub_t)
Root_sub_t$Genotype<-rownames(Root_sub_t)
GA_sub_t$Genotype<-rownames(GA_sub_t)
PA_sub_t$Genotype<-rownames(PA_sub_t)
DevInf_sub_t$Genotype<-rownames(DevInf_sub_t)
Ped_sub_t$Genotype<-rownames(Ped_sub_t)

Crown_sub_t%<>%relocate(Genotype)
Root_sub_t%<>%relocate(Genotype)
GA_sub_t%<>%relocate(Genotype)
PA_sub_t%<>%relocate(Genotype)
DevInf_sub_t%<>%relocate(Genotype)
Ped_sub_t%<>%relocate(Genotype)


Crown_sub_t<-left_join(Crown_sub_t, metadata, by="Genotype")
Root_sub_t<-left_join(Root_sub_t, metadata, by="Genotype")
GA_sub_t<-left_join(GA_sub_t, metadata, by="Genotype")
PA_sub_t<-left_join(PA_sub_t, metadata, by="Genotype")
DevInf_sub_t<-left_join(DevInf_sub_t, metadata, by="Genotype")
Ped_sub_t<-left_join(Ped_sub_t, metadata, by="Genotype")

Crown_sub_t<-Crown_sub_t%>%arrange(Haplotype, Year)
Root_sub_t<-Root_sub_t%>%arrange(Haplotype, Year)
GA_sub_t<-GA_sub_t%>%arrange(Haplotype, Year)
PA_sub_t<-PA_sub_t%>%arrange(Haplotype, Year)
DevInf_sub_t<-DevInf_sub_t%>%arrange(Haplotype, Year)
Ped_sub_t<-Ped_sub_t%>%arrange(Haplotype, Year)



Crown_sub_heatmap<-Crown_sub_t[, c(2:(ncol(Crown_sub_t)-14))]         
rownames(Crown_sub_heatmap)<-Crown_sub_t$Genotype
rownames(Crown_sub_heatmap)<-gsub("_", " ", rownames(Crown_sub_heatmap))
Crown_sub_anno<-data.frame("Cluster5H" = Crown_sub_t$Haplotype)
rownames(Crown_sub_anno) = rownames(Crown_sub_heatmap)
names(Crown_sub_anno)<-"Chr 5H haplotype"

Root_sub_heatmap<-Root_sub_t[, c(2:(ncol(Root_sub_t)-14))]         
rownames(Root_sub_heatmap)<-Root_sub_t$Genotype
rownames(Root_sub_heatmap)<-gsub("_", " ", rownames(Root_sub_heatmap))
Root_sub_anno<-data.frame("Cluster5H" = Root_sub_t$Haplotype)
rownames(Root_sub_anno) = rownames(Root_sub_heatmap)
names(Root_sub_anno)<-"Chr 5H haplotype"

DevInf_sub_heatmap<-DevInf_sub_t[, c(2:(ncol(DevInf_sub_t)-14))]         
rownames(DevInf_sub_heatmap)<-DevInf_sub_t$Genotype
rownames(DevInf_sub_heatmap)<-gsub("_", " ", rownames(DevInf_sub_heatmap))
DevInf_sub_anno<-data.frame("Cluster5H" = DevInf_sub_t$Haplotype)
rownames(DevInf_sub_anno) = rownames(DevInf_sub_heatmap)
names(DevInf_sub_anno)<-"Chr 5H haplotype"

Ped_sub_heatmap<-Ped_sub_t[, c(2:(ncol(Ped_sub_t)-14))]         
rownames(Ped_sub_heatmap)<-Ped_sub_t$Genotype
rownames(Ped_sub_heatmap)<-gsub("_", " ", rownames(Ped_sub_heatmap))
Ped_sub_anno<-data.frame("Cluster5H" = Ped_sub_t$Haplotype)
rownames(Ped_sub_anno) = rownames(Ped_sub_heatmap)
names(Ped_sub_anno)<-"Chr 5H haplotype"

GA_sub_heatmap<-GA_sub_t[, c(2:(ncol(GA_sub_t)-14))]         
rownames(GA_sub_heatmap)<-GA_sub_t$Genotype
rownames(GA_sub_heatmap)<-gsub("_", " ", rownames(GA_sub_heatmap))
GA_sub_anno<-data.frame("Cluster5H" = GA_sub_t$Haplotype)
rownames(GA_sub_anno) = rownames(GA_sub_heatmap)
names(GA_sub_anno)<-"Chr 5H haplotype"

PA_sub_heatmap<-PA_sub_t[, c(2:(ncol(PA_sub_t)-14))]         
rownames(PA_sub_heatmap)<-PA_sub_t$Genotype
rownames(PA_sub_heatmap)<-gsub("_", " ", rownames(PA_sub_heatmap))
PA_sub_anno<-data.frame("Cluster5H" = PA_sub_t$Haplotype)
rownames(PA_sub_anno) = rownames(PA_sub_heatmap)
names(PA_sub_anno)<-"Chr 5H haplotype"

Crown_sub_anno$`Chr 5H haplotype`<-as.factor(Crown_sub_anno$`Chr 5H haplotype`)


ann_colors = list(
  `Chr 5H haplotype` = c(`Haplotype 1` = haplo_colors[1], `Haplotype 2` = haplo_colors[2], `Haplotype 3` = haplo_colors[3], `Haplotype 4` = haplo_colors[4])
)

p1<-pheatmap(Crown_sub_heatmap,  cluster_cols = F, cluster_rows = F,  annotation_row = Crown_sub_anno, annotation_colors =ann_colors[1], annotation_legend=F, fontsize=5)
p2<-pheatmap(Root_sub_heatmap,  cluster_cols = F, cluster_rows = F, annotation_row = Root_sub_anno, annotation_colors =ann_colors[1], annotation_legend=F, fontsize=5)
p3<-pheatmap(DevInf_sub_heatmap,  cluster_cols = F, cluster_rows = F, annotation_row = DevInf_sub_anno, annotation_colors =ann_colors[1], annotation_legend=F, fontsize=5)
p4<-pheatmap(Ped_sub_heatmap,  cluster_cols = F, cluster_rows = F, annotation_row = Ped_sub_anno, annotation_colors =ann_colors[1], annotation_legend=F, fontsize=5)
p5<-pheatmap(GA_sub_heatmap,  cluster_cols = F, cluster_rows = F, annotation_row = GA_sub_anno, annotation_colors =ann_colors[1], annotation_legend=F, fontsize=5)
p6<-pheatmap(PA_sub_heatmap,  cluster_cols = F, cluster_rows = F, annotation_row = PA_sub_anno, annotation_colors =ann_colors[1], annotation_legend=F, fontsize=5)

p1gg <- as.ggplot(p1)
p2gg <- as.ggplot(p2)
p3gg <- as.ggplot(p3)
p4gg <- as.ggplot(p4)
p5gg <- as.ggplot(p5)
p6gg <- as.ggplot(p6)


legend<-ggplot(pa_ind[[1]], aes(x=PCo1,y=PCo2, color = Haplotype)) +
  geom_point( size = 2) +
  labs(x = paste0("PC1 (", pa_ind[[2]][1], "%)") , y = paste0("PC2 (", pa_ind[[2]][2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = haplo_colors)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7))+ labs(color = "Chr 5H haplotype  ") + theme(legend.box.background = element_rect(colour = "black") ,legend.title = element_text(size=7),
                                                                                                                                                  legend.text = element_text(size=6))

legend<-get_legend(legend)
legendplot<-as_ggplot(legend)


tiff("ESM9.tiff", height=28, width=11, units="cm", res=600)
ggdraw()+
  draw_plot(p1gg, x=, y=0, width=1, height=1)+
  draw_plot(legend, x = 0.9, y = 0.02, width = 0.01, height = 0.06)

dev.off()

tiff("ESM10.tiff",  height=28, width=11, units="cm", res=600)
ggdraw()+
  draw_plot(p2gg, x=, y=0, width=1, height=1)+
  draw_plot(legend, x = 0.9, y = 0.02, width = 0.01, height = 0.06)
dev.off()

tiff("ESM11.tiff",  height=28, width=11, units="cm", res=600)
ggdraw()+
  draw_plot(p3gg, x=, y=0, width=1, height=1)+
  draw_plot(legend, x = 0.9, y = 0.02, width = 0.01, height = 0.06)
dev.off()

tiff("ESM12.tiff",  height=28, width=11, units="cm", res=600)
ggdraw()+
  draw_plot(p4gg, x=, y=0, width=1, height=1)+
  draw_plot(legend, x = 0.9, y = 0.02, width = 0.01, height = 0.06)
dev.off()

tiff("ESM13.tiff",  height=28, width=11, units="cm", res=600)
ggdraw()+
  draw_plot(p5gg, x=, y=0, width=1, height=1)+
  draw_plot(legend, x = 0.9, y = 0.02, width = 0.01, height = 0.06)
dev.off()

tiff("ESM14.tiff",  height=28, width=11, units="cm", res=600)
ggdraw()+
  draw_plot(p6gg, x=, y=0, width=1, height=1)+
  draw_plot(legend, x = 0.9, y = 0.02, width = 0.01, height = 0.06)
dev.off()






