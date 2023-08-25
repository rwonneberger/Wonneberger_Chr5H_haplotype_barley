#####################################################################################################################
## Calculate the frequency of haplotype 2 in the IPK genebank barley collection for each country and plot on a map ##
#####################################################################################################################

source("https://raw.githubusercontent.com/rwonneberger/R_functions/main/Ronja_functions.R")

library(tibble)
library(maps)

# Make map plots

## Code from https://sarahpenir.github.io/r/making-maps/
plain <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  panel.background = element_rect(fill = "white"),
  plot.title = element_text(hjust = 0.5)
)

world <- map_data("world")


# This script requires a text file containing passport information of the gene bank accessions, as well as country information in a format that is compatible with the R package 'maps'.
# The file should contain the following columns: accession (containing the accession name, e.g. HORxxxx), status_genebank_information (e.g. landrace, cultivar), annual_growth_habit (winter, spring), row_type, and Region. Such a file can be found in Milner et al., 2019 

# You will also need a file containing the cluster/haplotype assignments based on the k-means clustering run on the genebank collection. Script ESM20_kmeans_genebank.R generates this dataset.

samples<-fread("Passport_data.txt", sep="\t") 

# Clean the dataset
samples%<>%filter(status_PCA != "outlier")

# Rename the samples so they match the accession names in the SNP file
table(samples$panel)
samples_Chinese<-samples%>%filter(panel == "Chinese_Genebank")
samples_Swiss<-samples%>%filter(panel == "Swiss_Genebank")
samples_IPK<-samples%>%filter(panel == "IPK_Genebank")

samples_Chinese$accession<-paste0("Chinese_Genebank_", samples_Chinese$accession)
samples_Swiss$accession<-paste0("Swiss_Genebank_", samples_Swiss$accession)
samples_IPK$accession<-paste0("BRIDGE_", samples_IPK$accession)
samples1<-rbind(samples_Chinese, samples_Swiss, samples_IPK )
#Leave out the Pourkheirandish_et_al_2015 lines because they're also included in the IPK collection


##################################################################################################################################
## Calculate the percentage of accessions carrying haplotype 2 per country, on subsets defined by e.g. row-type or growth habit ##
##################################################################################################################################

#################################################
## Landraces (LR) by row-type and growth habit ##
#################################################

samples_LR<-samples1%>%filter(status_genebank_information == "landrace")

## Make a df that contains the number of two-rowed, six-rowed-spring. winter etc LR per country
table(samples_LR$country_of_origin)
cluster_df1<-as.data.frame(table(samples_LR$country_of_origin))

table(samples_LR$country_of_origin, samples_LR$row_type)
cluster_df<-as.data.frame(table(samples_LR$country_of_origin, samples_LR$row_type))
cluster_df2<-dcast(cluster_df, Var1 ~Var2)

table(samples_LR$country_of_origin, samples_LR$annual_growth_habit)
cluster_df<-as.data.frame(table(samples_LR$country_of_origin, samples_LR$annual_growth_habit))
cluster_df3<-dcast(cluster_df, Var1 ~Var2)

clusters<-merge(cluster_df1, cluster_df2, by="Var1")
clusters<-merge(clusters, cluster_df3, by="Var1")
names(clusters)<-c("Country", "Total_lines_all", "unknown_rowtype_all", "2-rowed_all", "6-rowed_all", "deficiens_all", "intermedium_all", "labile_all", "unknown_growthhabit_all", "intermediate_all", "spring_all", "winter_all")

clusters$Country <- as.character(clusters$Country)


## Read in the clustering of the domesticated barleys
df<-fread( "PCA_all_domesticated_clusters.txt", sep="\t")

# Make a subset with only landraces having haplotype 2 and repeat the above calculations
df%<>%filter(Hap_cluster == "Hap 2", status_genebank_information == "landrace")

table(df$country_of_origin)
cluster_df1<-as.data.frame(table(df$country_of_origin))


table(df$country_of_origin, df$row_type)
cluster_df<-as.data.frame(table(df$country_of_origin, df$row_type))
cluster_df2<-dcast(cluster_df, Var1 ~Var2)

table(df$country_of_origin, df$annual_growth_habit)
cluster_df<-as.data.frame(table(df$country_of_origin, df$annual_growth_habit))
cluster_df3<-dcast(cluster_df, Var1 ~Var2)


clusters1<-merge(cluster_df1, cluster_df2, by="Var1")
clusters1<-merge(clusters1, cluster_df3, by="Var1")
clusters1<-add_column(clusters1, labile = NA , .after = 7) # add a dummy column "labile" to facilitate calculations

names(clusters1)<-c("Country", "Total_lines", "unknown_rowtype", "2-rowed", "6-rowed", "deficiens", "intermedium", "labile", "unknown_growthhabit", "intermediate", "spring", "winter")

clusters1$Country <- as.character(clusters1$Country)

# Fix some country names. It is problematic to replace Soviet Union by Russia but we don't have the exact origin of Soviet Union accessions but need to assign them to something
clusters$Country[clusters$Country == "SUN"] <- "RUS"
clusters$Country[clusters$Country == "DDR"] <- "DEU"
clusters$Country[clusters$Country == "GER"] <- "DEU"

clusters<-setDT(clusters)[, lapply(.SD, sum), .(Country)]

clusters1$Country[clusters1$Country == "SUN"] <- "RUS"
clusters1$Country[clusters1$Country == "DDR"] <- "DEU"
clusters1$Country[clusters1$Country == "GER"] <- "DEU"

clusters1<-setDT(clusters1)[, lapply(.SD, sum), .(Country)]

# Calculate percentages of accessions per group (defined by passport data) per country that carry haplotype 2
alldf<-left_join(clusters, clusters1, by="Country")

alldf<-as.data.frame(alldf)

alldf$perc_total          <-alldf[,13]/alldf[,2]*100
alldf$perc_unknown_rowtype<-alldf[,14]/alldf[,3]*100
alldf$perc_2rowed         <-alldf[,15]/alldf[,4]*100
alldf$perc_6rowed         <-alldf[,16]/alldf[,5]*100
alldf$perc_deficiens      <-alldf[,17]/alldf[,6]*100
alldf$perc_intermedium    <-alldf[,18]/alldf[,7]*100
alldf$perc_labile         <-alldf[,19]/alldf[,8]*100
alldf$perc_unknown_gh     <-alldf[,20]/alldf[,9]*100
alldf$perc_intermediate   <-alldf[,21]/alldf[,10]*100
alldf$perc_spring         <-alldf[,22]/alldf[,11]*100
alldf$perc_winter         <-alldf[,23]/alldf[,12]*100



as.data.frame(alldf)

alldf[is.na(alldf)] <- 0

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

alldf[is.nan(alldf)] <- 0


alldfperc<-alldf

regions<-fread("Maps_world_region.txt") # Helper file to add the full country name to the three-letter code used in the passport file. See the 'maps' documentation for details
alldfperc<-alldfperc[-1,]
finaldf<-left_join(regions, alldfperc, by="Country")


finaldf%<>%filter(Country != "unknown", Country != "Total", Country != "other" )

workdf<-finaldf

###################
## All landraces ##
###################

workdf[workdf$Total_lines_all < 15, ][,4:35] <- NA

names(workdf)[1]<-"region"

worldSubset <- inner_join(world, workdf, by = "region")

worldHap2_LR_all <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_total)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Hap 2")+
  plain + theme(legend.position="none")+ ggtitle("Landraces")+ theme(plot.title = element_text(size=9))+
  theme(plot.margin = unit(c(0.1,-1,0,-1), "cm"))


workdf<-finaldf

######################
## All two-rowed LR ##
######################

workdf[workdf$`2-rowed_all` < 15, ][,4:35] <- NA

names(workdf)[1]<-"region"

worldSubset <- inner_join(world, workdf, by = "region")

worldHap2_LR_2row <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_2rowed)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Hap 2")+
  plain + theme(legend.position="none") + ggtitle("Two-rowed landraces")+ theme(plot.title = element_text(size=9))+
  theme(plot.margin = unit(c(0.2,-1,0,-1), "cm"))


workdf<-finaldf

####################
## All 6-rowed LR ##
####################

workdf[workdf$`6-rowed_all` < 15, ][,4:35] <- NA

names(workdf)[1]<-"region"

worldSubset <- inner_join(world, workdf, by = "region")

worldHap2_LR_6row <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_6rowed)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Hap 2")+
  plain + theme(legend.position="none")+ ggtitle("Six-rowed landraces")+ theme(plot.title = element_text(size=9))+
  theme(plot.margin = unit(c(0.2,-1,0,-1), "cm"))


workdf<-finaldf

###################
## All spring LR ##
###################

workdf[workdf$spring_all < 15, ][,4:35] <- NA

names(workdf)[1]<-"region"

worldSubset <- inner_join(world, workdf, by = "region")

worldHap2_LR_spring <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_spring)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Hap 2")+
  plain + theme(legend.position="none")+ ggtitle("Spring landraces")+ theme(plot.title = element_text(size=9))+
  theme(plot.margin = unit(c(0.2,-1,0,-1), "cm"))


workdf<-finaldf

###################
## All winter LR ##
###################

workdf[workdf$winter_all < 15, ][,4:35] <- NA

names(workdf)[1]<-"region"

worldSubset <- inner_join(world, workdf, by = "region")

worldHap2_LR_winter <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_winter)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Hap 2")+
  plain + theme(legend.position="none")+ ggtitle("Winter landraces")+ theme(plot.title = element_text(size=9))+
  theme(plot.margin = unit(c(0.2,-1,0,-1), "cm"))




############################################
## Cultivars by row-type and growth habit ##
############################################

# Do the same calculations as above but this time for cultivars

###################
## All cultivars ##
###################

samples_cult<-samples1%>%filter(status_genebank_information == "cultivar")

## Make a df that contains the number of two-rowed, six-rowed-spring. winter etc LR per country
table(samples_cult$country_of_origin)
cluster_df1<-as.data.frame(table(samples_cult$country_of_origin))


table(samples_cult$country_of_origin, samples_cult$row_type)
cluster_df<-as.data.frame(table(samples_cult$country_of_origin, samples_cult$row_type))
cluster_df2<-dcast(cluster_df, Var1 ~Var2)

table(samples_cult$country_of_origin, samples_cult$annual_growth_habit)
cluster_df<-as.data.frame(table(samples_cult$country_of_origin, samples_cult$annual_growth_habit))
cluster_df3<-dcast(cluster_df, Var1 ~Var2)



clusters<-merge(cluster_df1, cluster_df2, by="Var1")
clusters<-merge(clusters, cluster_df3, by="Var1")
names(clusters)<-c("Country", "Total_lines_all", "unknown_rowtype_all", "2-rowed_all", "6-rowed_all", "deficiens_all", "intermedium_all", "labile_all", "unknown_growthhabit_all", "intermediate_all", "spring_all", "winter_all")

clusters$Country <- as.character(clusters$Country)

## Read in the clustering of the domesticated barleys
df<-fread( "PCA_all_domesticated_clusters.txt", sep="\t")

# Make a subset with only cultivars having haplotype 2 and repeat the above calculations
df%<>%filter(Hap_cluster == "Hap 2", status_genebank_information == "cultivar")

table(df$country_of_origin)
cluster_df1<-as.data.frame(table(df$country_of_origin))


table(df$country_of_origin, df$row_type)
cluster_df<-as.data.frame(table(df$country_of_origin, df$row_type))
cluster_df2<-dcast(cluster_df, Var1 ~Var2)

table(df$country_of_origin, df$annual_growth_habit)
cluster_df<-as.data.frame(table(df$country_of_origin, df$annual_growth_habit))
cluster_df3<-dcast(cluster_df, Var1 ~Var2)



clusters1<-merge(cluster_df1, cluster_df2, by="Var1")
clusters1<-merge(clusters1, cluster_df3, by="Var1")
clusters1<-add_column(clusters1, labile = NA , .after = 7) # add a dummy column "labile" to facilitate calculations

names(clusters1)<-c("Country", "Total_lines", "unknown_rowtype", "2-rowed", "6-rowed", "deficiens", "intermedium", "labile", "unknown_growthhabit", "intermediate", "spring", "winter")

clusters1$Country <- as.character(clusters1$Country)


# Fix some country names. It is problematic to replace Soviet Union by Russia but we don't have the exact origin of Soviet Union accessions but need to assign them to something
clusters$Country[clusters$Country == "SUN"] <- "RUS"
clusters$Country[clusters$Country == "DDR"] <- "DEU"
clusters$Country[clusters$Country == "GER"] <- "DEU"

clusters<-setDT(clusters)[, lapply(.SD, sum), .(Country)]

clusters1$Country[clusters1$Country == "SUN"] <- "RUS"
clusters1$Country[clusters1$Country == "DDR"] <- "DEU"
clusters1$Country[clusters1$Country == "GER"] <- "DEU"

clusters1<-setDT(clusters1)[, lapply(.SD, sum), .(Country)]

# Calculate percentages of accessions per group (defined by passport data) per country that carry haplotype 2

alldf<-left_join(clusters, clusters1, by="Country")

alldf<-as.data.frame(alldf)

alldf$perc_total          <-alldf[,13]/alldf[,2]*100
alldf$perc_unknown_rowtype<-alldf[,14]/alldf[,3]*100
alldf$perc_2rowed         <-alldf[,15]/alldf[,4]*100
alldf$perc_6rowed         <-alldf[,16]/alldf[,5]*100
alldf$perc_deficiens      <-alldf[,17]/alldf[,6]*100
alldf$perc_intermedium    <-alldf[,18]/alldf[,7]*100
alldf$perc_labile         <-alldf[,19]/alldf[,8]*100
alldf$perc_unknown_gh     <-alldf[,20]/alldf[,9]*100
alldf$perc_intermediate   <-alldf[,21]/alldf[,10]*100
alldf$perc_spring         <-alldf[,22]/alldf[,11]*100
alldf$perc_winter         <-alldf[,23]/alldf[,12]*100


as.data.frame(alldf)

alldf[is.na(alldf)] <- 0

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

alldf[is.nan(alldf)] <- 0


alldfperc<-alldf

regions<-fread("Maps_world_region.txt") # Helper file to add the full country name to the three-letter code used in the passport file. See the 'maps' documentation for details
alldfperc<-alldfperc[-1,]
finaldf<-left_join(regions, alldfperc, by="Country")



finaldf%<>%filter(Country != "unknown", Country != "Total", Country != "other" )

workdf<-finaldf

###################
## All cultivars ##
###################

workdf[workdf$Total_lines_all < 15, ][,4:35] <- NA

names(workdf)[1]<-"region"

worldSubset <- inner_join(world, workdf, by = "region")

worldHap2_cult_all <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_total)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Hap 2")+
  plain + theme(legend.position="none")+ ggtitle("Cultivars") + theme(plot.title = element_text(size=9))+
  theme(plot.margin = unit(c(0.1,-1,0,-1), "cm"))


workdf<-finaldf

###########################
## All 2-rowed cultivars ##
###########################

workdf[workdf$`2-rowed_all` < 15, ][,4:35] <- NA

names(workdf)[1]<-"region"

worldSubset <- inner_join(world, workdf, by = "region")

worldHap2_cult_2row <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_2rowed)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Hap 2")+
  plain + theme(legend.position="none") + ggtitle("Two-rowed cultivars")+ theme(plot.title = element_text(size=9))+
  theme(plot.margin = unit(c(0.2,-1,0,-1), "cm"))


workdf<-finaldf

###########################
## All 6-rowed cultivars ##
###########################

workdf[workdf$`6-rowed_all` < 15, ][,4:35] <- NA

names(workdf)[1]<-"region"

worldSubset <- inner_join(world, workdf, by = "region")

worldHap2_cult_6row <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_6rowed)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Hap 2")+
  plain + theme(legend.position="none")+ ggtitle("Six-rowed cultivars")+ theme(plot.title = element_text(size=9))+
  theme(plot.margin = unit(c(0.2,-1,0,-1), "cm"))


workdf<-finaldf

##########################
## All spring cultivars ##
##########################

workdf[workdf$spring_all < 15, ][,4:35] <- NA

names(workdf)[1]<-"region"

worldSubset <- inner_join(world, workdf, by = "region")

worldHap2_cult_spring <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_spring)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Hap 2")+
  plain + theme(legend.position="none")+ ggtitle("Spring cultivars")+ theme(plot.title = element_text(size=9))+
  theme(plot.margin = unit(c(0.2,-1,0,-1), "cm"))


workdf<-finaldf

##########################
## All winter cultivars ##
##########################

workdf[workdf$winter_all < 15, ][,4:35] <- NA

names(workdf)[1]<-"region"

worldSubset <- inner_join(world, workdf, by = "region")

worldHap2_cult_winter <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_winter)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Hap 2")+
  plain + theme(legend.position="none")+ ggtitle("Winter cultivars")+ theme(plot.title = element_text(size=9))+
  theme(plot.margin = unit(c(0.2,-1,0,-1), "cm"))


legend <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = perc_winter)) +
  scale_fill_distiller(palette ="YlOrRd", direction = 1, limits = c(0, 100)) + # or direction=1
  labs(fill = "% Haplotype 2")+
  plain + theme(legend.position="bottom")+ ggtitle("Winter cultivars") +
  theme(legend.key.size = unit(0.4, 'cm'))+  theme(legend.title = element_text(size=7),
                                                   legend.text = element_text(size=7))

legend<-get_legend(legend)
legendplot<-as_ggplot(legend)


#################################
## All landraces and cultivars ##
#################################

tiff("Fig10.tiff", width=8.4, height=13, units="cm", res=600)

ggarrange(worldHap2_LR_all, worldHap2_cult_all, legend,ncol=1, heights=c(1,1,0.2),  labels = c("a", "b"),  font.label = list(color = "black", size = 8))

dev.off()

#########################################
## Landraces and cultivars by category ##
#########################################

tiff("ESM24.tiff", width=17.4, height=23.4, units="cm", res=600)

legend_y = 0.05
height = (1-legend_y)/4
pos1 = legend_y
pos2 = pos1 + height
pos3 = pos2 + height
pos4 = pos3 + height

ggdraw()+
  draw_plot(worldHap2_LR_2row, x = 0, y = pos4, width = 0.5, height = height)+
  draw_plot(worldHap2_cult_2row, x = 0.5, y = pos4, width = 0.5, height = height)+
  draw_plot(worldHap2_LR_6row, x = 0, y = pos3, width = 0.5, height = height) +
  draw_plot(worldHap2_cult_6row, x = 0.5, y = pos3, width = 0.5, height = height)+
  draw_plot(worldHap2_LR_spring, x = 0, y = pos2, width = 0.5, height = height)+
  draw_plot(worldHap2_cult_spring, x = 0.5, y = pos2, width = 0.5, height = height)+
  draw_plot(worldHap2_LR_winter, x = 0, y = pos1, width = 0.5, height = height)+
  draw_plot(worldHap2_cult_winter, x = 0.5, y = pos1, width = 0.5, height = height)+
  draw_plot(legend, x = 0, y = 0, width = 1, height = legend_y)+
  draw_plot_label(label = c("a", "b", "c", "d", "e", "f", "g", "h"), size = 12,
                  x = c(0, 0.5, 0, 0.5, 0, 0.5, 0, 0.5), y = c(1, 1, pos4, pos4,  pos3, pos3, pos2, pos2 ))

dev.off()



