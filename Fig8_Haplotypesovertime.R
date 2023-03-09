####################################################
## Development of haplotype frequencies over time ##
####################################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

# Read in a metadata file generated from Online Resource 1. 
df<-fread("Metadata.txt", na.strings = "")

df<-df[!is.na(df$Haplotype),]

# % of haplotypes per decade

df1<-table(df$Range, df$Haplotype)
df1<-as.data.table(df1)
names(df1)<-c("Range", "Haplotype", "Count")
df1<-reshape2::dcast(df1, Range ~ Haplotype)

df2<-df1
df2[,2]<-df1$`Haplotype 1`/rowSums(df1[,2:4])*100
df2[,3]<-df1$`Haplotype 2`/rowSums(df1[,2:4])*100
df2[,4]<-df1$`Haplotype 3`/rowSums(df1[,2:4])*100
df2[,5]<-df1$`Haplotype 4`/rowSums(df1[,2:4])*100

df3<-melt(df2)

tiff("Fig8.tiff", width=8.4, height=6, res=600, units = "cm")
ggplot(df3, aes(x=Range, y=value, group=variable)) + geom_line(aes(color=variable)) + theme_classic() +  
  theme(axis.text.x = element_text(angle = 45, vjust=.5)) + scale_color_manual(name="Haplotype", labels=c("Hap 1","Hap 2","Hap 3", "Hap 4"), values=c('#E69F00', '#009E73', '#56B4E9', "#D55E00")) + xlab("Release period")+ ylab("% of lines")+  
  theme(axis.text=element_text(size=6),axis.title=element_text(size=7))+ 
  theme(legend.title = element_text(size=7),
        legend.text = element_text(size=6))
dev.off()

