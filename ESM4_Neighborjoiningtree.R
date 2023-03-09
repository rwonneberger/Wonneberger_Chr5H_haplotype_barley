###########################
## Neighbor-joining tree ##
###########################

# Code taken from https://epirhandbook.com/en/phylogenetic-trees-1.html and modified

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library(ape)
library(ggtree)
library(ggnewscale)
library(adegenet)

start<-68780000
end<-320040000


# Read in numeric SNP data 
# Read in SNP matrix published in Schreiber et al., 2023 as hapmap file. The genotypes need to be encoded as A, C, G, T

num<-fread("SNP_numeric.txt")
num<-as.data.frame(num)

hmp<-fread("Hapmap_haploid.hmp.txt")

hapnames<-names(hmp)[-c(1:11)]

# Read in a metadata file generated from Online Resource 1. 
metadata <- read.table("Metadata.txt", head=TRUE, sep="\t")

# Only keep samples present in all files
SamplesKeep <- colnames(hmp)[-c(1:11)]
metadata_reduced <- metadata[metadata$Line%in%SamplesKeep,]

metadata_reduced %<>% arrange(factor(Line, levels = hapnames))

names(hmp)[-c(1:11)] == metadata_reduced$Line


colorBlindGrey8   <- c( "#D55E00", "#E69F00","#F0E442","#009E73", 
                        "#0072B2","#56B4E9","#CC79A7", "#999999")


# Make distance matrix to calculate the tree
geno_dist <- dist(num)
names(hmp)[-c(1:11)]<-gsub("_", " ", names(hmp)[-c(1:11)])
dimnames(geno_dist)<-names(hmp)[-c(1:11)]

# Calculate tree
tree <- nj(geno_dist)

metadata_reduced$Line<-gsub("_", " ", metadata_reduced$Line)

# Add tip labels
tree2 <- tree %>% 
  left_join(metadata_reduced, by=c("label" = "Line"))


# Get the dataframes for additional information to be shown in circles around the phylogenetic tree

Cluster_5H <- data.frame("Cluster 5H" = metadata_reduced$Haplotype)
rownames(Cluster_5H) <- metadata_reduced$Line


metadata_reduced$Cluster_k3<-paste0("Cluster ", metadata_reduced$Cluster_k3)

Clusterk3 <- data.frame("Cluster at k = 3" = metadata_reduced$Cluster_k3)
rownames(Clusterk3) <- metadata_reduced$Line


metadata_reduced$Cluster_k2<-paste0("Cluster ", metadata_reduced$Cluster_k2)

Clusterk2 <- data.frame("Cluster at k = 2" = metadata_reduced$Cluster_k2)
rownames(Clusterk2) <- metadata_reduced$Line


tiff("ESM4.tiff", width=40, height=40, res=600, units="cm")

# PLot the tree
p<-ggtree(tree2, layout='circular', branch.length='none') + geom_tiplab(aes(color=Range)) + 
  scale_colour_manual(values = colorBlindGrey8)  +
  labs( col="Release period") 

# Plot the haplotypes

h1 <-  gheatmap(p, Cluster_5H,                                 
                offset = 9,                               
                width = 0.05,                             
                color = NULL,                         
                colnames = FALSE) +                    
  scale_fill_manual(name = "Haplotype 5H",             
                    values = c('#E69F00', '#009E73', '#56B4E9', "#D55E00"),
                    breaks = c("Haplotype 1", "Haplotype 2", "Haplotype 3", "Haplotype 4"),
                    labels = c("Haplotype 1", "Haplotype 2", "Haplotype 3", "Haplotype 4")) +
  theme(legend.position = "bottom")

h2 <- h1 + new_scale_fill() 

# Plot the clustering at k = 2
h3 <-  gheatmap(h2, Clusterk2,                          
                offset = 13,                       
                width = 0.05,                 
                color = NULL,                  
                colnames = FALSE) +                         
  scale_fill_manual(name = "Cluster at k=2",                   
                    values = c("#800000",  "#4363d8"),
                    breaks = c("Cluster 1", "Cluster 2"),
                    labels = c("Cluster 1", "Cluster 2")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin())

h4 <- h3 + new_scale_fill() 

# Plot the clustering at k = 3
h5 <- gheatmap(h4, Clusterk3,                           
               offset = 11,                          
               width = 0.05,                           
               color = NULL,                            
               colnames = FALSE) +                              
  scale_fill_manual(name = "Cluster at k=3",                       
                    values = c("#800000",  "#4363d8", "#ffe119"),
                    breaks = c("Cluster 1", "Cluster 2", "Cluster 3"),
                    labels = c("Cluster 1", "Cluster 2", "Cluster 3")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin())
h5

dev.off()

