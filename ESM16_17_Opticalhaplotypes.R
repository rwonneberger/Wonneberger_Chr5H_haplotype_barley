#######################################################################################
## Plot the SNP alleles for each cultivar in the region of interest on chromosome 5H ##
#######################################################################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library(pheatmap)
library(scales)


colorBlindGrey8   <- c("#999999", "#D55E00", "#E69F00","#F0E442","#009E73", "#0072B2","#56B4E9","#CC79A7")

################################################
## Plot the entire haplotype region +- 10 kbp ##
################################################

start<-58780000
end<-330040000

# Read in numeric SNP data
num<-fread("SNP_numeric.txt")
num<-as.data.frame(num)

# Read in SNP matrix published in Schreiber et al., 2023 as hapmap file. The genotypes need to be encoded as A, C, G, T
hmp<-fread("Hapmap_haploid.hmp.txt")

hapnames<-names(hmp)[-c(1:11)]

# Read in a metadata file generated from Online Resource 1. 
metadata <- read.table("Metadata.txt", head=TRUE, sep="\t")

# Only keep samples present in both files
SamplesKeep <- colnames(hmp)[-c(1:11)]
metadata_reduced <- metadata[metadata$Line%in%SamplesKeep,]

metadata_reduced<-metadata_reduced[match(names(hmp)[-c(1:11)], metadata_reduced$Line),]

table(metadata_reduced$Line %in% names(hmp)[-c(1:11)])
table(names(hmp)[-c(1:11)] %in% metadata_reduced$Line) 

table(names(hmp)[-c(1:11)]  == metadata_reduced$Line)


# Filter SNPs located in the region of interest
hmp_sub<-hmp%>%filter(chrom == "5H", pos >= start, pos <=end)

num_sub<-num[, hmp_sub$`rs#`]

geno_dist <- dist(num_sub)

#Perform multidimensional scaling (mds)
mds_geno_ind <- cmdscale(geno_dist,k=2,eig=TRUE)
mds_geno_ind_df <- data.frame(mds_geno_ind$points,check.names = TRUE)
colnames(mds_geno_ind_df) <- c("PCo1","PCo2")

# Add passport data
mds_geno_ind_df$accession<-metadata_reduced$Line
mds_geno_ind_df$Haplotype<-metadata_reduced$Haplotype

loc<-transpose_df_Col1(hmp_sub[,-c(2:11)], "accession")

table(loc$accession == mds_geno_ind_df$accession)

df1<-cbind(num_sub, mds_geno_ind_df[, 4])
df1%<>%relocate(`mds_geno_ind_df[, 4]`)
names(df1)[1]<-"Haplotype"

df1%<>%arrange(Haplotype)

df1$number<-1:nrow(df1)
df1%<>%relocate(number)
df1$accession<-NULL

# Make a random subset of 10k markers for visualization
cols<-sort(sample(c(3:ncol(df1)), 10000) ) 
df2<-df1[, c(1:4, cols)]

df2<-reshape2::melt(df2, id.vars = c("Haplotype", "number"))

colors=c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075' )

haplo_colors=c('#E69F00', '#009E73', '#56B4E9', "#D55E00")

df2$value<-as.factor(df2$value)

df2$variable<-as.numeric(gsub("S5H_", "", df2$variable))


p1<-ggplot(df2, aes(variable/1000000, number, fill=value, color=value)) +
  geom_tile()  +
  theme_classic()+
  scale_colour_manual(values = colors)+
  scale_fill_manual(values = colors) + 
  theme(legend.position = "none") +
  xlab("Position (Mbp)")+
  ylab("Accession") +
  theme(axis.text.x = element_text(angle = 90, vjust=.5),
        axis.line = element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
        ) +
  scale_x_continuous(expand = c(0, 0))


p2<-ggplot(df2, aes( x=100, y=number, fill=Haplotype, color=Haplotype))+      geom_tile()  +
  theme_classic() +
  scale_fill_manual(name="Haplotype", values = haplo_colors)+ 
  scale_colour_manual(name="Haplotype", values = haplo_colors)+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust=.5, color = "white"),
        axis.title.x = element_text(color="white"),
        axis.ticks.x=element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
  scale_x_continuous( breaks =100, labels = 100, expand = c(0, 0))


tiff("ESM16.tiff", height=8, width=17.4, res=600, units="cm")
figure<-ggarrange(p1, p2, nrow=1, widths=c(1, 0.04))  
figure
dev.off()


#########################################
## Zoom into the start and end regions ##
#########################################

#########################
## Start of the region ##
#########################

#################################################################
## 30 Mbp upstream and 4 Mnp downstream of the haplotype start ##
#################################################################

start<-38780000
end<-72780000

hmp_sub<-hmp%>%filter(chrom == "5H", pos >= start, pos <=end)

num_sub<-num[, hmp_sub$`rs#`]

geno_dist <- dist(num_sub)

#Perform multidimensional scaling (mds)
mds_geno_ind<- cmdscale(geno_dist,k=2,eig=TRUE)
mds_geno_ind_df <- data.frame(mds_geno_ind$points,check.names = TRUE)
colnames(mds_geno_ind_df) <- c("PCo1","PCo2")


mds_geno_ind_df$accession<-metadata_reduced$Line
mds_geno_ind_df$Haplotype<-metadata_reduced$Haplotype

loc<-transpose_df_Col1(hmp_sub[,-c(2:11)], "accession")


table(loc$accession == mds_geno_ind_df$accession)

df1<-cbind(num_sub, mds_geno_ind_df[, 4])
df1%<>%relocate(`mds_geno_ind_df[, 4]`)
names(df1)[1]<-"Haplotype"

df1%<>%arrange(Haplotype)

df1$number<-1:nrow(df1)
df1%<>%relocate(number)
df1$accession<-NULL



df2<-df1

df2<-reshape2::melt(df2, id.vars = c("Haplotype", "number"))

colors=c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075' )

haplo_colors=c('#E69F00', '#009E73', '#56B4E9', "#D55E00")

df2$value<-as.factor(df2$value)


df2$variable<-as.numeric(gsub("S5H_", "", df2$variable))


p1<-ggplot(df2, aes(variable/1000000, number, fill=value, color=value)) +
  geom_tile() +
  theme_classic()+
  scale_colour_manual(values = colors)+
  scale_fill_manual(values = colors) + 
  theme(legend.position = "none") +
  xlab("Position (Mbp)")+
  ylab("Accession") +
  theme(axis.text.x = element_text(angle = 90, vjust=.5),
        axis.line = element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0))+ 
  geom_vline(xintercept = 68.78, size=0.7)


p2<-ggplot(df1, aes( x=10, y=number, fill=Haplotype, color=Haplotype))+      geom_tile()  +
  theme_classic() +
  scale_fill_manual(name="Haplotype", values = haplo_colors)+ 
  scale_colour_manual(name="Haplotype", values = haplo_colors)+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust=.5, color = "white"),
        axis.title.x = element_text(color="white"),
        axis.ticks.x=element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
  scale_x_continuous( breaks =10, labels = 10, expand = c(0, 0))



tiff("Rasterplot_start.tiff", height=8, width=17.4, res=100, units="cm")
figure<-ggarrange(p1, p2,ncol=2, widths=c(1,0.04))  
figure
dev.off()


# Save this for a combined plot
figure1<-figure  



#######################
## End of the region ##
#######################

###############################################################
## 4 Mbp upstream and 30 Mnp downstream of the haplotype end ##
###############################################################

start<-316040000
end<-350040000

hmp_sub<-hmp%>%filter(chrom == "5H", pos >= start, pos <=end)


num_sub<-num[, hmp_sub$`rs#`]

geno_dist <- dist(num_sub)

#Perform multidimensional scaling (mds)
mds_geno_ind <- cmdscale(geno_dist,k=2,eig=TRUE)
mds_geno_ind_df <- data.frame(mds_geno_ind$points,check.names = TRUE)
colnames(mds_geno_ind_df) <- c("PCo1","PCo2")


mds_geno_ind_df$accession<-metadata_reduced$Line
mds_geno_ind_df$Haplotype<-metadata_reduced$Haplotype

loc<-transpose_df_Col1(hmp_sub[,-c(2:11)], "accession")


table(loc$accession == mds_geno_ind_df$accession)

df1<-cbind(num_sub, mds_geno_ind_df[, 4])
df1%<>%relocate(`mds_geno_ind_df[, 4]`)
names(df1)[1]<-"Haplotype"

df1%<>%arrange(Haplotype)

df1$number<-1:nrow(df1)
df1%<>%relocate(number)
df1$accession<-NULL


df2<-df1

df2<-reshape2::melt(df2, id.vars = c("Haplotype", "number"))

colors=c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075' )

haplo_colors=c('#E69F00', '#009E73', '#56B4E9', "#D55E00")

df2$value<-as.factor(df2$value)


df3<-reshape2::melt(df3, id.vars = c("Haplotype", "number"))


p1<-ggplot(df2, aes(variable/1000000, number, fill=value, color=value)) +
  geom_tile()  +
  theme_classic()+
  scale_colour_manual(values = colors)+
  scale_fill_manual(values = colors) + 
  theme(legend.position = "none") +
  xlab("Position (Mbp)")+
  ylab("Accession") +
  theme(axis.text.x = element_text(angle = 90, vjust=.5),
        axis.line = element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0))+ 
  geom_vline(xintercept = 319.99, size=0.7)


p2<-ggplot(df1, aes( x=100, y=number, fill=Haplotype, color=Haplotype))+      geom_tile()  +
  theme_classic() +
  scale_fill_manual(name="Haplotype", values = haplo_colors)+ 
  scale_colour_manual(name="Haplotype", values = haplo_colors)+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust=.5, color = "white"),
        axis.title.x = element_text(color="white"),
        axis.ticks.x=element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
  scale_x_continuous(breaks = 100, labels=100, expand = c(0, 0))


tiff("Rasterplot_end.tiff", height=8, width=17.4, res=100, units="cm")
figure<-ggarrange(p1, p2,ncol=2, widths=c(1,0.04))  
figure
dev.off()


# Save this for a combined plot
figure2<-figure  




############################################
## Start of the region, zoomed in further ##
############################################

####################################################################
## 0.5 Mbp upstream and 0.5 Mnp downstream of the haplotype start ##
####################################################################

start<-68280000
end<-69280000

hmp_sub<-hmp%>%filter(chrom == "5H", pos >= start, pos <=end)

num_sub<-num[, hmp_sub$`rs#`]

geno_dist <- dist(num_sub)

#Perform multidimensional scaling (mds)
mds_geno_ind<- cmdscale(geno_dist,k=2,eig=TRUE)
mds_geno_ind_df <- data.frame(mds_geno_ind$points,check.names = TRUE)
colnames(mds_geno_ind_df) <- c("PCo1","PCo2")

mds_geno_ind_df$accession<-metadata_reduced$Line
mds_geno_ind_df$Haplotype<-metadata_reduced$Haplotype

loc<-transpose_df_Col1(hmp_sub[,-c(2:11)], "accession")

table(loc$accession == mds_geno_ind_df$accession)

df1<-cbind(num_sub, mds_geno_ind_df[, 4])
df1%<>%relocate(`mds_geno_ind_df[, 4]`)
names(df1)[1]<-"Haplotype"

df1%<>%arrange(Haplotype)

df1$number<-1:nrow(df1)
df1%<>%relocate(number)
df1$accession<-NULL

df2<-df1

df2<-reshape2::melt(df2, id.vars = c("Haplotype", "number"))

colors=c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075' )

haplo_colors=c('#E69F00', '#009E73', '#56B4E9', "#D55E00")

df2$value<-as.factor(df2$value)


df2$variable<-as.numeric(gsub("S5H_", "", df2$variable))


p1<-ggplot(df2, aes(variable/1000000, number, fill=value, color=value)) +
  geom_tile()  +
  theme_classic()+
  scale_colour_manual(values = colors)+
  scale_fill_manual(values = colors) + 
  theme(legend.position = "none") +
  xlab("Position (Mbp)")+
  ylab("Accession") +
  theme(axis.text.x = element_text(angle = 90, vjust=.5),
        axis.line = element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0))+ 
  geom_vline(xintercept = 68.78, size=0.7)

scaleFUN <- function(x) sprintf("%.2f", x)

p2<-ggplot(df1, aes( x=10, y=number, fill=Haplotype, color=Haplotype))+      geom_tile()  +
  theme_classic() +
  scale_fill_manual(name="Haplotype", values = haplo_colors)+ 
  scale_colour_manual(name="Haplotype", values = haplo_colors)+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust=.5, color = "white"),
        axis.title.x = element_text(color="white"),
        axis.ticks.x=element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
  scale_x_continuous(labels=scaleFUN, expand = c(0, 0))



tiff("Rasterplot_start_fine.tiff", height=8, width=17.4, res=100, units="cm")
figure<-ggarrange(p1, p2,ncol=2, widths=c(1,0.04))  
figure
dev.off()


# Save this for a combined plot
figure3<-figure  



##########################################
## End of the region, zoomed in further ##
##########################################

##################################################################
## 0.5 Mbp upstream and 0.5 Mnp downstream of the haplotype end ##
##################################################################

start<-319540000
end<-320540000

hmp_sub<-hmp%>%filter(chrom == "5H", pos >= start, pos <=end)

num_sub<-num[, hmp_sub$`rs#`]

geno_dist <- dist(num_sub)

#Perform multidimensional scaling (mds)
mds_geno_ind <- cmdscale(geno_dist,k=2,eig=TRUE)
mds_geno_ind_df <- data.frame(mds_geno_ind$points,check.names = TRUE)
colnames(mds_geno_ind_df) <- c("PCo1","PCo2")

mds_geno_ind_df$accession<-metadata_reduced$Line
mds_geno_ind_df$Haplotype<-metadata_reduced$Haplotype

loc<-transpose_df_Col1(hmp_sub[,-c(2:11)], "accession")


table(loc$accession == mds_geno_ind_df$accession)

df1<-cbind(num_sub, mds_geno_ind_df[, 4])
df1%<>%relocate(`mds_geno_ind_df[, 4]`)
names(df1)[1]<-"Haplotype"

df1%<>%arrange(Haplotype)

df1$number<-1:nrow(df1)
df1%<>%relocate(number)
df1$accession<-NULL


df2<-df1

df2<-reshape2::melt(df2, id.vars = c("Haplotype", "number"))

colors=c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075' )

haplo_colors=c('#E69F00', '#009E73', '#56B4E9', "#D55E00")

df2$value<-as.factor(df2$value)

df2$variable<-as.numeric(gsub("S5H_", "", df2$variable))

p1<-ggplot(df2, aes(variable/1000000, number, fill=value, color=value)) +
  geom_tile()  +
  theme_classic()+
  scale_colour_manual(values = colors)+
  scale_fill_manual(values = colors) + 
  theme(legend.position = "none") +
  xlab("Position (Mbp)")+
  ylab("Accession") +
  theme(axis.text.x = element_text(angle = 90, vjust=.5),
        axis.line = element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0))+ 
  geom_vline(xintercept = 319.99, size=0.7)

scaleFUN <- function(x) sprintf("%.2f", x)

p2<-ggplot(df1, aes( x=100, y=number, fill=Haplotype, color=Haplotype))+      geom_tile()  +
  theme_classic() +
  scale_fill_manual(name="Haplotype", values = haplo_colors)+ 
  scale_colour_manual(name="Haplotype", values = haplo_colors)+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust=.5, color = "white"),
        axis.title.x = element_text(color="white"),
        axis.ticks.x=element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  ) +
scale_x_continuous(labels=scaleFUN, expand = c(0, 0))

tiff("Rasterplot_end_fine.tiff", height=8, width=17.4, res=100, units="cm")
figure<-ggarrange(p1, p2,ncol=2, widths=c(1,0.04))  
figure
dev.off()


# Save this for a combined plot
figure4<-figure  


tiff("ESM17.tiff", width=17.4, height=24, units="cm", res=600)


ggarrange(figure1, figure2, figure3, figure4,  nrow=4, labels = c("a", "b", "c", "d"))

dev.off()
