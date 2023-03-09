##################################################
## AMOVA of cultivars grouped by release period ##
##################################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library(adegenet)
library(poppr)

## AMOVA with release period as strata

# Read in SNP matrix published in Schreiber et al., 2023 as hapmap file. The genotypes need to be encoded as A, C, G, T

hmp<-fread("Hapmap_diploid.hmp.txt", na.strings = "NN")

# Make a subset of 600000 SNPs or whatever your computational resources allow
hmp<-hmp[sample(nrow(hmp), 600000), ]

hmp_sub<-as.data.frame(hmp[, -c(2:11)])

hmp_sub_t<-transpose_df(hmp_sub)
hmp_sub_t$Line<-colnames(hmp_sub)[-1]
hmp_sub_t%<>%relocate(Line)
hmp_sub_t%<>%arrange(Line)

# # Read in a metadata file generated from Online Resource 1 and define release periods
metadata <- read.table("Metadata.txt", head=TRUE, check.names=FALSE, sep="\t", na.strings="")

metadata$Range<-""
metadata[metadata$Year < 1960,]$Range<-"1830-1959"
metadata[metadata$Year >= 1960 & metadata$Year < 1980,]$Range<-"1960-1979"
metadata[metadata$Year >= 1980 & metadata$Year < 2000,]$Range<-"1980-1999"
metadata[metadata$Year >= 2000,]$Range<-"2000-2014"

metadata%<>%arrange(Line)

SamplesKeep <- names(hmp)[-c(1:11)]
metadata_reduced <- metadata[metadata$Line%in%SamplesKeep,]

hmp_sub_t_reduced<-hmp_sub_t[ hmp_sub_t$Line%in%metadata_reduced$Line,]

table(hmp_sub_t_reduced$Line == metadata_reduced$Line)

hmp_sub_t_join<-cbind(metadata_reduced, hmp_sub_t_reduced)


locus<-hmp_sub_t_join[,-c(1:8)]

# Make genind object
ind <- as.character(hmp_sub_t_join$Line) 
population <- as.character(hmp_sub_t_join$Range)
genind_obj <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep = "")


#Add strata
Range<-hmp_sub_t_join$Range
stra<-data.frame(Range)
strata(genind_obj)<-stra
nameStrata(genind_obj)<- ~Range

# Make genclone object
genclone_obj<-as.genclone(genind_obj)

# Run AMOVA using the release periods as groups
table(strata(genclone_obj, ~Range, combine = FALSE))
amova_Range <- poppr.amova(genind_obj, ~Range, threads=18, missing="mean", within = FALSE)

set.seed(1999)
rand_Range <- randtest(amova_Range, nrepet = 999)
