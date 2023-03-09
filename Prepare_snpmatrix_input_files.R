library(genoCN)
library(data.table)

# All analyses of variant data are based on the hapmap file published in Schreiber et al., 2023 and available at https://ics.hutton.ac.uk/germinate-barn/. The genotypes are encoded as A, C, G, T. For some of the analyses, the data needs to be encoded as AA, CC, GG, TT, and for other analyses a numeric file is needed. 
# This script converts the original file into the needed formats.

# Convert the hapmap file to diploid encoding (run on command line, requires Tassel):

# run_pipeline.pl -Xms10g -Xmx100g -fork1 -h Hapmap_haploid.hmp.txt -export Hapmap_diploid -exportType HapmapDiploid -runfork1


# Convert the hapmap file to numeric format and transpose:

# Path to Genotype hapmap file
hmp <- read.table("Hapmap_diploid.hmp.txt", header=T, stringsAsFactors=FALSE,na.strings = "NN",comment.char = "", check.names=FALSE)


# Transforming genotype variant file to numeric format

hmp_t <- t(hmp[,c(1,12:dim(hmp)[2])])
colnames(hmp_t) <- hmp_t[1,]
hmp_t <- hmp_t[-1,]
output <- matrix(NA, nrow=dim(hmp_t)[1], ncol=dim(hmp_t)[2])

for(i in 1:dim(hmp_t)[2]){
  
  output[,i] <- code.genotype(hmp_t[,i])
}

rownames(output) <- rownames(hmp_t)
colnames(output) <- colnames(hmp_t)

fwrite(output, "SNP_numeric.txt", sep="\t")




