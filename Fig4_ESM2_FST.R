#######################################################
## Plot FST values between clusters along the genome ##
#######################################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

window="100kb"
step="10kb"

###########################################################################
## Pairwise FST between k-means clusters in windows along the chromosome ##
###########################################################################

# This script requires the outputs from the FST_per_kmeanscluster.sh script

# Read in FST data
FST_k2<-fread(paste0("FST_DAPC_k2_C1_C2_", window, "_window_", step, "_step.windowed.weir.fst"))

FST_k3_C1_C2<-fread(paste0("FST_DAPC_k3_C1_C2_", window, "_window_", step, "_step.windowed.weir.fst"))
FST_k3_C1_C3<-fread(paste0("FST_DAPC_k3_C1_C3_", window, "_window_", step, "_step.windowed.weir.fst"))
FST_k3_C2_C3<-fread(paste0("FST_DAPC_k3_C2_C3_", window, "_window_", step, "_step.windowed.weir.fst"))

FST_k2$Range<-"k1 vs k2"

FST_k3_C1_C2$Range<-"C1 vs C2"
FST_k3_C1_C3$Range<-"C1 vs C3"
FST_k3_C2_C3$Range<-"C2 vs C3"

# Remove windows wiht less than 5 variants
FST_k2<-FST_k2[FST_k2$CHROM != "UN" & FST_k2$N_VARIANTS >= 5,]

FST_k3_C1_C2<-FST_k3_C1_C2[FST_k3_C1_C2$CHROM != "UN" & FST_k3_C1_C2$N_VARIANTS >= 5,]
FST_k3_C1_C3<-FST_k3_C1_C3[FST_k3_C1_C3$CHROM != "UN" & FST_k3_C1_C3$N_VARIANTS >= 5,]
FST_k3_C2_C3<-FST_k3_C2_C3[FST_k3_C2_C3$CHROM != "UN" & FST_k3_C2_C3$N_VARIANTS >= 5,]

# Merge files and remove chromosome "unmapped"
FST_k3<-rbind(FST_k3_C1_C2, FST_k3_C1_C3, FST_k3_C2_C3)
FST_k3<-FST_k3[FST_k3$CHROM != "UN",]

FST_k3$Range<-as.factor(FST_k3$Range)

#### FST k=2 and 3 together
scaleFUN <- function(x) sprintf("%.2f", x)

data_vline_p1<-data.frame(CHROM = c("3H", "4H", "5H", "5H" ),
                          vline = c(267.7, 34.4,  68.78, 320.04))

data_vline_centro<-data.frame(CHROM = c("1H", "2H", "3H", "4H", "5H",  "6H", "7H" ),
                              vline = c(211.460194, 298.511764, 270.896959, 273.503980, 214.765154, 252.168825, 329.584237))



p1<-ggplot(FST_k2,aes(x=BIN_START/1000000,y=WEIGHTED_FST))  + geom_point(size=0.4)+   theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
  scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200), labels = function(x) format(x, scientific = FALSE))+
  theme(axis.text.x = element_text(angle = 90, vjust=.5)) + 
  facet_grid(~CHROM, scales='free_x', space="free_x")+ xlab("Position (Mbp)")+ ylab("FST")  +scale_y_continuous(labels=scaleFUN)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+geom_vline(data_vline_centro, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="lightgrey") + geom_vline(data_vline_p1, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red")

data_vline_p2<-data.frame(CHROM = c("3H", "4H", "5H", "5H", "6H", "7H"),
                          vline = c(267.7, 34.4,  68.78, 320.04, 248.2, 316.6 ))

data_vline_p3<-data.frame(CHROM = c("3H", "4H", "5H", "5H" ),
                          vline = c(267.7, 34.4,  68.78, 320.04))
data_vline_p4<-data.frame(CHROM = c("3H", "3H",  "6H", "7H" ),
                          vline = c( 254.2, 381.1, 248.2, 316.6 ))


data_vline<-list(data_vline_p2, data_vline_p3, data_vline_p4)

scaleFUN <- function(x) sprintf("%.2f", x)

p=list()
count = 0

for (i in unique(FST_k3$Range)) {
  count = count +1
  nr1 <- substr(i, start = 2, stop = 2)
  nr2 <- substr(i, start = 8, stop = 8)
  if (count < 3) {
    sub<-FST_k3[FST_k3$Range == i, ]
    
    p[[i]]<-ggplot(sub,aes(x=BIN_START/1000000,y=WEIGHTED_FST)) + geom_point(size=0.4)+
      theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
      scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_y_continuous(labels=scaleFUN)+
      facet_grid(~CHROM, scales='free_x', space="free_x")+ xlab("Position (Mbp)")+ ylab("FST")+scale_y_continuous(labels=scaleFUN)+geom_vline(data_vline_centro, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="lightgrey")  + geom_vline(data_vline[[count]], mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red") } else  {
        sub<-FST_k3[FST_k3$Range == i, ]
        
        p[[i]]<-ggplot(sub,aes(x=BIN_START/1000000,y=WEIGHTED_FST))  + geom_point(size=0.4)+
          scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
          theme_classic() + theme(panel.grid.minor = element_line(colour = "grey90"), panel.grid.major = element_line(colour = "grey90"), panel.border = element_rect(fill=NA))+ 
          scale_x_continuous(breaks = seq(from = 0, to = 600, by = 200))+
          theme(axis.text.x = element_text(angle = 90, vjust=.5)) +
          scale_y_continuous(labels=scaleFUN)+
          facet_grid(~CHROM, scales='free_x', space="free_x")+ xlab("Position (Mbp)") + ylab("FST") +scale_y_continuous(labels=scaleFUN)+geom_vline(data_vline_centro, mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="lightgrey")+ geom_vline(data_vline[[count]], mapping=aes(xintercept = vline), linetype="dotted", size=1.1, color="red")
      }
  print(p[[i]])
}


tiff("Fig4.tiff", height=18, width=17.4, units="cm", res=600)
ggarrange(p1, p[[1]], p[[2]], p[[3]],  nrow=4, heights=c(1, 1, 1, 1.3), labels = c("a", "b", "c", "d"))

dev.off()



#################################################################################
## Calculate overall FST between groups of cultivars clustered by release year ##
#################################################################################

library(hierfstat)

# Read in SNP matrix published in Schreiber et al., 2023 as hapmap file. The genotypes need to be encoded as AA, CC, GG, TT

hmp<-fread("Hapmap_diploid.hmp.txt", na.strings = "NN")

# Make subset of 400000 markers
set.seed(5)
hmp<-hmp[sample(nrow(hmp), 400000), ]

hmp_sub<-as.data.frame(hmp[, -c(2:11)])

hmp_sub_t<-transpose_df(hmp_sub)
hmp_sub_t$Line<-colnames(hmp_sub)[-1]
hmp_sub_t%<>%relocate(Line)
hmp_sub_t%<>%arrange(Line)

# Read in metadata and define release year ranges
metadata <- read.table("Metadata.txt", head=TRUE, check.names=FALSE, sep="\t", na.strings="")

metadata$Range<-""
metadata[metadata$Year < 1960,]$Range<-"1830-1959"
metadata[metadata$Year >= 1960 & metadata$Year < 1980,]$Range<-"1960-1979"
metadata[metadata$Year >= 1980 & metadata$Year < 2000,]$Range<-"1980-1999"
metadata[metadata$Year >= 2000,]$Range<-"2000-2014"


metadata%<>%arrange(Line)

# Make sure that the accession names in the snpmatrix and the metadata are identical
SamplesKeep <- names(hmp)[-c(1:11)]
metadata_reduced <- metadata[metadata$Line%in%SamplesKeep,]

hmp_sub_t_reduced<-hmp_sub_t[ hmp_sub_t$Line%in%metadata_reduced$Line,]

table(hmp_sub_t_reduced$Line == metadata_reduced$Line)

df_join<-merge(metadata_reduced, hmp_sub_t_reduced)



df_join<-df_join[, -c(1:3, 5:7)]

df_join<-as.matrix(df_join)

df_join<-gsub("AA", "11", df_join)
df_join<-gsub("CC", "22", df_join)
df_join<-gsub("GG", "33", df_join)
df_join<-gsub("TT", "44", df_join)
df_join<-gsub("AC", "12", df_join)
df_join<-gsub("AG", "13", df_join)
df_join<-gsub("AT", "14", df_join)
df_join<-gsub("CA", "21", df_join)
df_join<-gsub("AG", "23", df_join)
df_join<-gsub("AT", "24", df_join)
df_join<-gsub("GA", "31", df_join)
df_join<-gsub("GC", "32", df_join)
df_join<-gsub("GT", "34", df_join)
df_join<-gsub("TA", "41", df_join)
df_join<-gsub("TC", "42", df_join)
df_join<-gsub("TG", "13", df_join)

df_join<-as.data.frame(df_join)

grp<-df_join[,1]
geno<-df_join[,-1]

geno[] <- sapply( geno, as.numeric )

alldata<-cbind(grp, geno)


# Calculate FST between groups
fst<-pairwise.WCfst(alldata)

fst<-as.data.frame(fst)


fst


