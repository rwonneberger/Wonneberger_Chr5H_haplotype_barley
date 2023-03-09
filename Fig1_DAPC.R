#############################################################
## DAPC analysis using release period as pre-defined group ##
#############################################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library(adegenet)
library(poppr)


colorBlindGrey8   <- c( "#D55E00", "#E69F00","#F0E442","#009E73", 
                        "#0072B2","#56B4E9","#CC79A7", "#999999")



# Read in SNP matrix published in Schreiber et al., 2023 as hapmap file. The genotypes need to be encoded as AA, CC, GG, TT

hmp<-fread("Hapmap_diploid.hmp.txt", na.strings = "NN")

# Read in a metadata file generated from Online Resource 1. 
metadata <- read.table("Metadata.txt", head=TRUE, check.names=FALSE, sep="\t", na.strings="")

# Make a random subset of 600000 markers
set.seed(5)
hmp<-hmp[sample(nrow(hmp), 600000), ] 

hmp_sub<-as.data.frame(hmp[, -c(2:11)])

hmp_sub_t<-transpose_df(hmp_sub)
hmp_sub_t$Line<-colnames(hmp_sub)[-1]
hmp_sub_t%<>%relocate(Line)%>%arrange(Line)


# Make sure that the accession names in the snpmatrix and the metadata are identical
metadata%<>%arrange(Line)

SamplesKeep <- names(hmp)[-c(1:11)]
metadata_reduced <- metadata[metadata$Line%in%SamplesKeep,]

hmp_sub_t_reduced<-hmp_sub_t[ hmp_sub_t$Line%in%metadata_reduced$Line,]

table(hmp_sub_t_reduced$Line == metadata_reduced$Line)

# Merge both datasets
df_join<-cbind(metadata_reduced, hmp_sub_t_reduced)
df_join<-mydata_join[, -c(15)] # Some of these lines need to be changed depending on how you make your metadata file and how many columns it contains. 

# Extract only the SNP data
locus<-df_join[,-c(1:15)] # Some of these lines need to be changed depending on how you make your metadata file and how many columns it contains. 

# Make a genind object
ind <- as.character(mydata_join$Line) 
population <- as.factor(mydata_join$Range)
genind_obj <- df2genind(locus, ploidy = 2,ind.names = ind,pop=population, sep = "")


# add strata
Range<-mydata_join$Range
stra<-data.frame(Range)
strata(genind_obj)<-stra
nameStrata(genind_obj)<- ~Range

# Make a genclone object
genclone_obj<-as.genclone(genind_obj)


# Run DAPC cross validation. Warning, this might run for a few days! Reduce the number of markers if your computational resources are limited
datax<-xvalDapc(tab(genclone_obj, NA.method="mean"), pop(MyGC), n.pca.max = 200, training.set = 0.9,
                result = "groupMean", 2 = TRUE, scale = FALSE,
                n.pca = NULL, n.rep = 100, xval.plot = TRUE, parallel="multicore")


# Optimal number of PCs is 40

# Run dapc  
dapc_out<-dapc(genind_obj, n.pca=40, n.da=6)


# Prepare data for plotting
plotdf<-as.data.frame(dapc_out$ind.coord)

plotdf$Line<-rownames(plotdf)

plotdf$Line == mydata_join$Line

plotdf$Range<-mydata_join$Range


mds.var.per<-round(dapc_out$eig/sum(dapc_out$eig)*100, 1)



tiff("Fig1.tiff", width = 8.4, height = 7, units = "cm", res = 600 )
ggplot(plotdf, aes(x=LD1,y=LD2, color = Range)) +
  stat_ellipse()+
  geom_point(size=0.5) +
  labs(x = paste0("PC1 (", mds.var.per[], "%)") , y = paste0("PC2 (", mds.var.per[2], "%)"))+
  theme_classic()+
  scale_colour_manual(values = colorBlindGrey8)+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7))+ labs(color="Release period") +  
  theme(legend.title = element_text(size=7),
        legend.text = element_text(size=6))


dev.off()