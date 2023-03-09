##########################################################################################################################
## k-means clustering to identify the most likely number of subpopulations and assign cultivars to these subpopulations ##
##########################################################################################################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")
library(adegenet)

# Read in SNP matrix published in Schreiber et al., 2023 as hapmap file. The genotypes need to be encoded as AA, CC, GG, TT

hmp<-fread("Hapmap_diploid.hmp.txt", na.strings = "NN")

# Make a subset of 1000000 SNPs or whatever your computational resources allow
hmp<-hmp[sample(nrow(hmp), 1000000), ]

hmp_sub<-as.data.frame(hmp[, -c(2:11)])

hmp_sub_t<-transpose_df(hmp_sub)
hmp_sub_t$Line<-colnames(hmp_sub)[-1]
hmp_sub_t%<>%relocate(Line)
hmp_sub_t%<>%arrange(Line)

# Read in a metadata file generated from Online Resource 1. 
metadata <- read.table("Metadata.txt", head=TRUE, check.names=FALSE, sep="\t", na.strings="")

metadata%<>%arrange(Line)

# Only keep samples present in all files
SamplesKeep <- names(hmp)[-c(1:11)]
metadata_reduced <- metadata[metadata$Line%in%SamplesKeep,]

hmp_sub_t_reduced<-hmp_sub_t[ hmp_sub_t$Line%in%metadata_reduced$Line,]

table(hmp_sub_t_reduced$Line == metadata_reduced$Line)

df_join<-cbind(metadata_reduced, hmp_sub_t_reduced)


locus<-df_join[,-c(1:8)]

ind <- as.character(df_join$Line) 

# Make genind object
genind_obj <- df2genind(locus, ploidy = 2,ind.names = ind, sep = "")
genind_obj


# Select the number of PCs that explain  ~90% of the variance

tiff("ESM6a.tiff", width=8.4, height=8.4, res=100, units = "cm")
grp<-find.clusters(genind_obj, max.n.clust=50) #110 PCs
dev.off()

tiff("ESM6b.tiff", width=8.4, height=8.4, res=100, units = "cm")
grp<-find.clusters(genind_obj, max.n.clust=50) #110 PCs
dev.off()


grp$grp # Contains the group membership of each cultivar at the  chosen number of k

