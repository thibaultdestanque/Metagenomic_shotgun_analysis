#!/usr/bin/Rscript



###########################
#######  Metagenomic ######
###########################

# Clean environnement
ls()
rm(list=ls())
ls()



#######################
# Installing packages #
#######################

install.packages(c("vegan", "metacoder", "taxa", "ggplot2", "dplyr", "readr", "stringr", "agricolae", "ape"),
                 repos = "http://cran.rstudio.com",
                 dependencies = TRUE)


library(BiocManager)
BiocManager::install("microbiome")
library(microbiome)  



####################
# Loading packages #
####################

library(Matrix)
library(purrr)
library("ps")
library("ade4")
library(readr)
library("igraph")
library("phyloseq")

library(microbiome)  
# Some function are masked... see link below to know what it mean :
# https://stackoverflow.com/questions/39137110/what-does-the-following-object-is-masked-from-packagexxx-mean


library(vegan)
library(metacoder)
library(taxa)
library(ggplot2)
library(dplyr)
library(readr) # pas fonctionn�
library(stringr) # pas fonctionn� ?
library(agricolae)
library(ape)

#install.packages("gplots")
library("gplots")

# to install packages from Bioconductor:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.12")

# don't work:
#biocLite("Heatplus")  # annHeatmap or annHeatmap2
#install.packages("Heatplus")
#library(Heatplus)

# load the vegan package for hierachical clustering if you want to use distance functions not specified in dist.
#install.packages("vegan")
library(vegan)

# load the RColorBrewer package for better colour options
#install.packages("RColorBrewer")
library(RColorBrewer)
library("plyr")
library("pheatmap")
library("RColorBrewer")
library("reshape2")
library("dplyr")
library(colorspace) # get nice colors
library(lattice)
library(tidyr) # Pas fonctionn�
library(knitr)
library(ggplot2)
library(ggfortify)
library(stringr)
library(gridExtra)
library(grid)
library("Rmisc")


# d = read.delim(file = "clipboard", stringsAsFactors = F, check.names = F, header = TRUE)
# rownames(d) = d$what
# d$what = NULL
# d
# grid.table(d)

###########################
## Set working directory ##
###########################

setwd("~/01_MetAmox_AntiSelfish_Project/03_results/01_3_Taxonomic_affiliation_reads")



#################
## Import data ##
#################

df_taxo_affi_reads_phylum   = read.table("1_taxo_affi_reads_phylum.tsv", sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE)
df_taxo_affi_reads_class    = read.table("2_taxo_affi_reads_class.tsv",  sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE)
df_taxo_affi_reads_order    = read.table("3_taxo_affi_reads_order.tsv",  sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE)
df_taxo_affi_reads_family   = read.table("4_taxo_affi_reads_family.tsv", sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE)
df_taxo_affi_reads_genus    = read.table("5_taxo_affi_reads_genus.tsv",  sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE)

# Eliminer avec un sed les characters speciaux qui pose pb avant de faire l'import:
# cat taxo_affi_reads_species.tsv | sed 's/'\''//g' > taxo_affi_reads_species_modif.tsv # to remoove "'" special characters
df_taxo_affi_reads_species  = read.table("6_taxo_affi_reads_species_modif.tsv",sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE)



#######################
### Import metadata ###
#######################

### WARNING ### 
# > get rid of special characters in metadata.txt
metadata    = data.frame(read.table("metadata_modif.txt",  sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE))
rownames(metadata) = metadata[,"names_metadata"] 

### TO DO ###
# replace Na by 0 in number of inject column and volume






                          #>    #########################################     <#
                          #>    ### Choose level of taxonomy to study ###     <#
                          #>    #########################################     <#


# Parameter to change:
Level_tax = "genus"
# Possible values: "phylum", "order", "family", "genus", "species"

# Affect level of taxo to study (do not touch):
if (Level_tax=="species"){
  data_used = df_taxo_affi_reads_species
} else if(Level_tax=="genus"){
  data_used = df_taxo_affi_reads_genus
} else if(Level_tax=="family"){
  data_used = df_taxo_affi_reads_family
} else if(Level_tax=="order"){
  data_used = df_taxo_affi_reads_order
} else if(Level_tax=="phylum"){
  data_used = df_taxo_affi_reads_phylum
}
                          


###############
# FORMAT DATA #
###############

# set first column as rownames
rownames(data_used) = data_used[,1]
data_used = data_used[,-1]

# Select specific columns "reads":
tmp = data_used %>% select(starts_with('reads')) 



##########################
# NB of reads per sample #
##########################

nb_reads_per_sample = tmp

nb_reads_tot_per_sample = as.data.frame(colSums(nb_reads_per_sample))
colnames(nb_reads_tot_per_sample)="Nb_of_reads"

rownames(nb_reads_tot_per_sample) = gsub("reads_", "", rownames(nb_reads_tot_per_sample)) # Retrieve only counts of reads.
nb_reads_tot_per_sample$Samples = rownames(nb_reads_tot_per_sample)

tmp_metadata = metadata # need to change colname of metadata for this figure
colnames(tmp_metadata)[1] = "Samples"

nb_reads_tot_per_sample_info = join(nb_reads_tot_per_sample, tmp_metadata) # join the 2 data frame 
sum(nb_reads_tot_per_sample_info$Nb_of_reads) # make the sum of column Nb_of_reads

# Plot number of reads per samples
ggplot(nb_reads_tot_per_sample_info, aes(fill=Prelevement, y=Nb_of_reads, x=Samples)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2, size = 7))

# Plot number of reads per Calf
ggplot(nb_reads_tot_per_sample_info, aes(fill=Prelevement, y=Nb_of_reads, x=Calf)) +
 geom_bar(stat="identity", position=position_dodge()) +
 theme(axis.text.x = element_text(angle=90))



                          #>    ##########################################     <#
                          #>    ### FORMAT DATA 2: choose specific run ###     <#
                          #>    ##########################################     <#

# Parameter to change:

# Filter of maximum relative abundance desire:
my_filter=0.01

Keepdata = "Metamox"
# Possible values: "MetamoxAntiselfi", "Metamox", "Antiselfi", "metamox_24-11-21_988385074", "metamox_24-11-21_995018042"

# select specific data base on colnames (do not touch):
if (Keepdata=="MetamoxAntiselfi"){
  # Nothing to do here, we keep all the data ;) 
} else if(Keepdata=="Metamox"){
  tmp = tmp %>% select(contains('Metamox')) 
} else if(Keepdata=="Antiselfi"){
  tmp = tmp %>% select(contains('Antiselfi')) 
} else if(Keepdata=="Metamox_24-11-21_988385074"){
  tmp = tmp %>% select(contains('Metamox_24-11-21_988385074')) 
} else if(Keepdata=="Metamox_24-11-21_995018042"){
  tmp = tmp %>% select(contains('Metamox_24-11-21_995018042')) 
}

# remove unclassified reads
which(rownames(tmp)=="unclassified")
tmp2 = tmp[-which(rownames(tmp)=="unclassified"),]
which(rownames(tmp2)=="unclassified")

# Remove "cannot be assigned to a (non-viral)" base on phylum, order, family, genus, species
if (Level_tax=="species"){
  which(rownames(tmp2)=="cannot be assigned to a (non-viral) species")
  tmp3 = tmp2[-which(rownames(tmp2)=="cannot be assigned to a (non-viral) species"),]
  which(rownames(tmp3)=="cannot be assigned to a (non-viral) species")
} else if(Level_tax=="genus"){
  which(rownames(tmp2)=="cannot be assigned to a (non-viral) genus")
  tmp3 = tmp2[-which(rownames(tmp2)=="cannot be assigned to a (non-viral) genus"),]
  which(rownames(tmp3)=="cannot be assigned to a (non-viral) genus")
} else if(Level_tax=="family"){
 which(rownames(tmp2)=="cannot be assigned to a (non-viral) family")
  tmp3 = tmp2[-which(rownames(tmp2)=="cannot be assigned to a (non-viral) family"),]
  which(rownames(tmp3)=="cannot be assigned to a (non-viral) family")
} else if(Level_tax=="order"){
  which(rownames(tmp2)=="cannot be assigned to a (non-viral) order")
  tmp3 = tmp2[-which(rownames(tmp2)=="cannot be assigned to a (non-viral) order"),]
  which(rownames(tmp3)=="cannot be assigned to a (non-viral) order")
} else if(Level_tax=="phylum"){
  which(rownames(tmp2)=="cannot be assigned to a (non-viral) phylum")
  tmp3 = tmp2[-which(rownames(tmp2)=="cannot be assigned to a (non-viral) phylum"),]
  which(rownames(tmp3)=="cannot be assigned to a (non-viral) phylum")
}

# Remove Viruses
which(rownames(tmp3)=="Viruses")
tmp4 = tmp3[-which(rownames(tmp3)=="Viruses"),]
which(rownames(tmp4)=="Viruses")

# Transpose data 
dataf_t = as.data.frame(t(tmp4))



###############################
## Nb of reads after filter: ##
###############################

nb_reads_per_sample = tmp4
nb_reads_tot_per_sample = as.data.frame(colSums(nb_reads_per_sample))
colnames(nb_reads_tot_per_sample)="Nb_of_reads"
rownames(nb_reads_tot_per_sample) = gsub("reads_", "", rownames(nb_reads_tot_per_sample)) # Retrieve only counts of reads.
nb_reads_tot_per_sample$Samples = rownames(nb_reads_tot_per_sample)
tmp_metadata = metadata
colnames(tmp_metadata)[1] = "Samples"
nb_reads_tot_per_sample_info = join(nb_reads_tot_per_sample, tmp_metadata)
sum(nb_reads_tot_per_sample_info$Nb_of_reads)

# Plot number of reads per samples
ggplot(nb_reads_tot_per_sample_info, aes(fill=Prelevement, y=Nb_of_reads, x=Samples)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2))

# Plot number of reads per Calf
ggplot(nb_reads_tot_per_sample_info, aes(fill=Prelevement, y=Nb_of_reads, x=Calf)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle=90))



###########################
### compare run metamox ###
###########################

# if (Keepdata=="Metamox"){
#   # Compare runs
#   df_74 <- dataf_t[which(rownames(dataf_t)=="Run_pool-Metamox_24-11-21_988385074"),]
#   df_42 <- dataf_t[which(rownames(dataf_t)=="Run_pool-Metamox_24-11-21_995018042"),]
#   
#   
#   tmp = data_used %>% select(starts_with('reads')) 
#   
#   
#   tmp <- rownames(dataf_t) %in% select(contains('Run_pool-Metamox_24-11-21_988385074'))
#   
#   
#   select(contains('Metamox'))
#   
#   
#   
# } else {
#   # Nothing to do
# }






######################
### Aggregate data ###
######################

# Encapsule cette partie sans modifier les autres (=> pas de modifications faite sur dataf_t)
my_data = dataf_t
rownames(my_data) = gsub("reads_", "", rownames(my_data)) # Retrieve only counts of reads.
my_data$info = rownames(my_data) # Copy rownames to a column "info"

# Visualize data
head(rownames(metadata))
head(rownames(my_data))
colnames(metadata)[1] = "info" # Rename first column as "info"
metadata$info
my_data$info

# Join information of metadata to my_data
my_data_info = join(my_data, metadata)

# PENSER A SEPARER LES ECHANTILLONS ANTIBIO / TEMOINS
# Pour antiselfi pas de témoins: seulement les échantillons 3,7,9,11,14,15
# A voir pour les échantillons MetAmox Metamox_24-11-21_988385074 : 2,6,11,15 et témoins: TN3, TN4, TN5
# A voir pour les échantillons MetAmox Metamox_24-11-21_995018042 : 2,6,11,15 et témoins: TN3, TN4, TN5

# S�pare les T�moins des echantillons

TN_data = my_data_info[which(str_detect(my_data_info$Calf, "TN")),]
Samples_data = my_data_info[which(str_detect(my_data_info$Calf, "V")),]

#### separe calfs and TN
 
Cond = c("TN", "V")


m<-list()
n<-list()
o<-list()
p<-list()


pdf("Barplot.pdf")
  for (i in Cond) {
    print(i)

    # Retrieve info of interest
    my_data_info_cond = my_data_info[which(str_detect(my_data_info$Calf, i)),]

    # Retrieve "Prelevement" needed in "info" column
    my_data_info_cond$info = my_data_info_cond$Prelevement 

    # Selectionne juste le tableau avec les infos utiles: valeurs et colonne "info"
    my_data_info_cond = my_data_info_cond[,1:which(colnames(my_data_info_cond)=="info")] 

    # Aggregate les counts pour chaque condition "prelevement"
    my_data_aggregate = aggregate(.~info, my_data_info_cond, sum) 
    rownames(my_data_aggregate) = my_data_aggregate$info # change les rownames par la colonne "info" (Prelevement)
    my_data_aggregate$info = NULL # Supprime la colonne info
    head(my_data_aggregate)
    
    # Normalisation
    # my_data_aggregate_t = data.frame(t(my_data_aggregate)) # transpose le data frame
    my_data_aggregate_prop <- my_data_aggregate/rowSums(my_data_aggregate)
    head(my_data_aggregate_prop)

    # Determine the maximum relative abundance for each column
    maxab <- apply(my_data_aggregate_prop, 2, max)
    head(maxab)

    # Remove the species with less than X% (variable my_filter) as their maximum relative abundance
    n1 <- names(which(maxab < my_filter))
    my_data_aggregate_prop_filter <- my_data_aggregate_prop[, -which(names(my_data_aggregate_prop) %in% n1)]

    # Melt data.frame for ggpplot 
    my_data_aggregate_prop_filter = as.matrix(my_data_aggregate_prop_filter)
    my_data_aggregate_prop_melt = melt(my_data_aggregate_prop_filter)
    colnames(my_data_aggregate_prop_melt) = c("Condition", "LevelTax", "Value")

    # Convert "factor"  to "character" : melt induce that by default...
    my_data_aggregate_prop_melt$LevelTax = as.character(my_data_aggregate_prop_melt$LevelTax)
    my_data_aggregate_prop_melt$Condition = as.character(my_data_aggregate_prop_melt$Condition)

    # Order by column LevelTax : alphabetical order:
    my_data_aggregate_prop_melt <- my_data_aggregate_prop_melt[order(my_data_aggregate_prop_melt$LevelTax),]

    ###############
    ### BARPLOT ###
    ###############

    if (i == "TN") {
    my_title = "Temoins negatifs"
    } else if (i == "V") {
    my_title = "Injected Calfs"
    }

    #4 p[[i]]<-
   m[[i]]<-ggplot(my_data_aggregate_prop_melt, aes(fill=Condition, y=Value, x=LevelTax)) +
            geom_bar(position="fill", stat="identity") +
            theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2)) +
            ggtitle(paste("Barplot:", my_title)) + 
            coord_flip()
  
   n[[i]]<-ggplot(my_data_aggregate_prop_melt, aes(fill=Condition, y=Value, x=LevelTax)) +
            geom_bar(stat="identity", position=position_dodge()) +
            theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2)) +
            ggtitle(paste("Barplot:", my_title))+ 
     coord_flip()

   o[[i]]<-ggplot(my_data_aggregate_prop_melt, aes(fill=Condition, y=Value, x=LevelTax)) +
            geom_bar(stat="identity", position=position_dodge()) +
            ylim(0,0.05) +
            theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2)) +
            ggtitle(paste("Barplot:", my_title, "(zoom 0.05)"))
  
   p[[i]]<-ggplot(my_data_aggregate_prop_melt, aes(fill=LevelTax, y=Value, x=Condition)) +
            geom_bar(position="fill", stat="identity") +
            ggtitle(paste("Barplot:", my_title, "(inverse)"))
   
   
   # Choose specific value on LevelTax column:
   levelTax_ofInterest = c("Bacteroides", "Enterococcus", "Escherichia", "Ruminococcus")
   my_data_aggregate_prop_melt_interest = my_data_aggregate_prop_melt[which(my_data_aggregate_prop_melt$LevelTax %in% levelTax_ofInterest),]
   
   ggplot(my_data_aggregate_prop_melt_interest, aes(fill=Condition, y=Value, x=LevelTax)) +
     geom_bar(stat="identity", position=position_dodge()) +
     theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2)) +
     ggtitle(paste("Barplot:", my_title))+ 
     coord_flip()
   
  }
dev.off()

layout_m = matrix(c(1,2), nrow = 2, byrow = TRUE)
layout_n = matrix(c(1,2), nrow = 2, byrow = TRUE)
layout_o = matrix(c(1,2), nrow = 2, byrow = TRUE)
layout_p = matrix(c(1,2), nrow = 2, byrow = TRUE)

multiplot(m[[1]],m[[2]], layout = layout_m)
multiplot(n[[1]],n[[2]], layout = layout_n)
multiplot(o[[1]],o[[2]], layout = layout_o)
multiplot(p[[1]],p[[2]], layout = layout_p)


#### MIX Calfs and TN:

# Retrieve "Prelevement" needed in "info" column
my_data_info$info = my_data_info$Prelevement 

# Selectionne juste le tableau avec les infos utiles: valeurs et colonne "info"
my_data_info = my_data_info[,1:which(colnames(my_data_info)=="info")] 

# Aggregate les counts pour chaque condition "prelevement"
my_data_aggregate = aggregate(.~info, my_data_info, sum) 
rownames(my_data_aggregate) = my_data_aggregate$info # change les rownames par la colonne "info" (Prelevement)
my_data_aggregate$info = NULL # Supprime la colonne info
head(my_data_aggregate)

# my_data_aggregate_t = data.frame(t(my_data_aggregate)) # transpose le data frame
my_data_aggregate_prop <- my_data_aggregate/rowSums(my_data_aggregate)

# Determine the maximum relative abundance for each column
maxab <- apply(my_data_aggregate_prop, 2, max)
head(maxab)

# Remove the species with less than 1% (< 0.01) as their maximum relative abundance
n1 <- names(which(maxab < 0.01))
my_data_aggregate_prop_filter <- my_data_aggregate_prop[, -which(names(my_data_aggregate_prop) %in% n1)]

# Melt data.frame for ggpplot 
my_data_aggregate_prop_filter = as.matrix(my_data_aggregate_prop_filter)
my_data_aggregate_prop_melt = melt(my_data_aggregate_prop_filter)
colnames(my_data_aggregate_prop_melt) = c("Condition", "LevelTax", "Value")

# Convert "factor"  to "character" : melt induce that by default...
my_data_aggregate_prop_melt$LevelTax = as.character(my_data_aggregate_prop_melt$LevelTax)
my_data_aggregate_prop_melt$Condition = as.character(my_data_aggregate_prop_melt$Condition)

# Order by column LevelTax : alphabetical order:
my_data_aggregate_prop_melt <- my_data_aggregate_prop_melt[order(my_data_aggregate_prop_melt$LevelTax),]



###############
### BARPLOT ###
###############

ggplot(my_data_aggregate_prop_melt, aes(fill=Condition, y=Value, x=LevelTax)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2))

ggplot(my_data_aggregate_prop_melt, aes(fill=Condition, y=Value, x=LevelTax)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2))

ggplot(my_data_aggregate_prop_melt, aes(fill=Condition, y=Value, x=LevelTax)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ylim(0,0.05) +
  theme(axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2))

ggplot(my_data_aggregate_prop_melt, aes(fill=LevelTax, y=Value, x=Condition)) +
  geom_bar(position="fill", stat="identity")



###########################
# Apply a bigger filter ? #
###########################

## 10 % to see what happen for majority bacteria

# Remove the species with less than 10% (< 0.1) as their maximum relative abundance
n1 <- names(which(maxab < 0.1))
my_data_aggregate_prop_filter <- my_data_aggregate_prop[, -which(names(my_data_aggregate_prop) %in% n1)]

# Melt data.frame for ggpplot 
my_data_aggregate_prop_filter = as.matrix(my_data_aggregate_prop_filter)
my_data_aggregate_prop_melt = melt(my_data_aggregate_prop_filter)
colnames(my_data_aggregate_prop_melt) = c("Condition", "LevelTax", "Value")

# Convert "factor"  to "character" : melt induce that by default...
my_data_aggregate_prop_melt$LevelTax = as.character(my_data_aggregate_prop_melt$LevelTax)
my_data_aggregate_prop_melt$Condition = as.character(my_data_aggregate_prop_melt$Condition)

# Order by column LevelTax : alphabetical order:
my_data_aggregate_prop_melt <- my_data_aggregate_prop_melt[order(my_data_aggregate_prop_melt$LevelTax),]

# Make a barplot
ggplot(my_data_aggregate_prop_melt, aes(fill=LevelTax, y=Value, x=Condition)) +
  geom_bar(position="fill", stat="identity")



#################
###  HEATMAP  ###
#################

# Convert nb of reads in proportion
data.prop <- dataf_t/rowSums(dataf_t)
data.prop[1:3, 1:3] # show first lines and columns

# https://www.molecularecologist.com/2013/08/20/making-heatmaps-with-r-for-microbiome-analysis/
set.seed(1) # to fix the random number used in heatmap and always keep the same cluster after a restart session.

# colorRampPalette is in the RColorBrewer package.  This creates a colour palette that shades from light yellow to red in RGB space with 100 unique colours
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)

# basic heatmap
heatmap(as.matrix(data.prop), Rowv = NA, Colv = NA, col = scaleyellowred)

# determine the maximum relative abundance for each column
maxab <- apply(data.prop, 2, max)
head(maxab)

# remove the species with less than 1% (< 0.01) as their maximum relative abundance
n1 <- names(which(maxab < 0.01))
data.prop.1 <- data.prop[, -which(names(data.prop) %in% n1)]

# Remove "reads_" from rownames(data.prop.1)
rownames(data.prop.1) = gsub("reads_", "", rownames(data.prop.1))
#data.prop.1[,"names_metadata"] = rownames(data.prop.1)

########################################
# join the 2 files ?
#met = join(data.prop.1,metadata, )
#rownames(met) = met[,"names_metadata"]
########################################

# the margins command sets the width of the white space around the plot. The first element is the bottom margin and the second is the right margin
heatmap(as.matrix(data.prop.1), Rowv = NA, Colv = NA, col = scaleyellowred, margins = c(10, 2))

## add dendrogram for the sample with Vegan package:
# calculate the Bray-Curtis dissimilarity matrix on the full dataset:
data.dist <- vegdist(data.prop, method = "bray")

# Do average linkage hierarchical clustering. Options are 'aver', complete' or 'single'. You'll need to choose the one that best fits the needs of your situation and your data.
row.clus <- hclust(data.dist, "complete")

library(dendextend)
# make the heatmap with Rowv = as.dendrogram(row.clus)
#heatmap(as.matrix(data.prop.1),
#        Rowv = as.dendrogram(row.clus),
#        Colv = NA,
#        col = scaleyellowred,
#        xlab = "Phylum", main = "Heatmap filter 1%",
#        margins = c(16, 14))


dev.off()
par(cex.main=1, cex.lab=0.7, cex.axis=0.7)
heatmap.2(as.matrix(data.prop.1),
        xlab = "Phylum", main = "Heatmap filter 1%",
        col = scaleyellowred,
        margins = c(12, 18),
        cexRow = 0.6,
        cexCol = 0.8,
        Rowv = as.dendrogram(row.clus))


# Make comparison between P0 & P2 and P0 & P1
data_prop_P0 = data.prop.1[grep("P0", rownames(data.prop.1)), ]
data_prop_P1 = data.prop.1[grep("P1", rownames(data.prop.1)), ]
data_prop_P2 = data.prop.1[grep("P2", rownames(data.prop.1)), ]

data_prop_P0vsP1 = rbind(data_prop_P0,data_prop_P1) 
data_prop_P0vsP2 = rbind(data_prop_P0,data_prop_P2) 

# P0 & P1
data.dist_P0vsP1 <- vegdist(data_prop_P0vsP1, method = "bray")
row.clus_P0vsP1  <- hclust(data.dist_P0vsP1, "complete")

heatmap.2(as.matrix(data_prop_P0vsP1),
          xlab = "Phylum", main = "Heatmap filter 1%",
          col = scaleyellowred,
          margins = c(12, 18),
          cexRow = 0.6,
          cexCol = 0.8,
          Rowv = as.dendrogram(row.clus_P0vsP1))

# P0 & P2
data.dist_P0vsP2 <- vegdist(data_prop_P0vsP2, method = "bray")
row.clus_P0vsP2  <- hclust(data.dist_P0vsP2, "complete")

# Add color:
#join(data_prop_P0vsP2,metadata$Prelevement)
metadata_colorPrelevement = cbind(rownames(metadata), metadata[,"Prelevement"])

metadata_colorPrelevement <- replace(metadata_colorPrelevement, which(metadata_colorPrelevement == "P0"), "green")
metadata_colorPrelevement <- replace(metadata_colorPrelevement, which(metadata_colorPrelevement == "P2"), "red")
metadata_colorPrelevement # contain info for P0, P1 and P2

# as the same as above we only need P0vsP1 and P0vsP2
metadata_colorPrelevement_P0 = metadata_colorPrelevement[grep("P0", metadata_colorPrelevement[,1]), ]
metadata_colorPrelevement_P1 = metadata_colorPrelevement[grep("P1", metadata_colorPrelevement[,1]), ]
metadata_colorPrelevement_P2 = metadata_colorPrelevement[grep("P2", metadata_colorPrelevement[,1]), ]

metadata_colorPrelevement_P0vsP1 = rbind(metadata_colorPrelevement_P0,metadata_colorPrelevement_P1) 
metadata_colorPrelevement_P0vsP2 = rbind(metadata_colorPrelevement_P0,metadata_colorPrelevement_P2) 


# heatmap
heatmap.2(as.matrix(data_prop_P0vsP2),
          xlab = "tax", main = "Heatmap filter 1% : P0 vs P2",
          col = scaleyellowred,
          margins = c(12, 18),
          cexRow = 0.6,
          cexCol = 0.8,
          Rowv = as.dendrogram(row.clus_P0vsP2))

#heatmap.2(as.matrix(data_prop_P0vsP2),
#          xlab = "tax", main = "Heatmap filter 1% : P0 vs P2",
#          col = scaleyellowred,
#          margins = c(12, 18),
#          cexRow = 0.6,
#          cexCol = 0.8,
#          Rowv = as.dendrogram(row.clus_P0vsP2),
#          RowSideColors = metadata_colorPrelevement_P0vsP2[,1])


# Rowv = as.dendrogram(row.clus),



###########
### ACP ###
###########

## retrieve data from n1 filter step and remove "reads_"
data.prop.1 <- data.prop[, -which(names(data.prop) %in% n1)]
rownames(data.prop.1) = gsub("reads_", "", rownames(data.prop.1))

# Pass rowanmes as column to match with metadata
data.prop.1$info = row.names(data.prop.1)

## make a new column info with PrelevementCalf_Run
metadata$info = str_c(metadata$Prelevement, '', metadata$Calf, '_', metadata$Run)

head(data.prop.1)
head(metadata)

# Add info from metadata file
df_info = join(data.prop.1, metadata)

# Replace number>=1 for this column by "yes" and "NA" by "No":
df_info[which(df_info$No_of_amoxicillin_injections >= 1),"No_of_amoxicillin_injections"] = "yes"
df_info[which(is.na(df_info$No_of_amoxicillin_injections)),"No_of_amoxicillin_injections"] = "No"

rownames(df_info) = df_info$info
head(df_info)

# Prepare data frame for prcomp function:
df <- df_info[,1:which(colnames(df_info)=="info")-1]


###### Tests
# Select specific species
# df = df %>% select(starts_with('Lactobacillus'))
######

#df <- iris[1:4]
pca_res <- prcomp(df, scale. = TRUE)
plot(pca_res)

# "Calf", "Farm"

# Prelevement
autoplot(pca_res, data = df_info, colour = 'Prelevement')
autoplot(pca_res, data = df_info, colour = 'Prelevement', label.size = 2,shape = FALSE)


# Run: les runs metamox se superposent:
autoplot(pca_res, data = df_info, colour = 'Farm')
autoplot(pca_res, data = df_info, colour = 'Farm',label = TRUE, label.size = 3,shape = FALSE)
autoplot(pca_res, data = df_info, colour = 'Farm',label = TRUE, label.size = 3,shape = FALSE, loadings = TRUE) # drive par pas mal de donn�es.

# try a different package for better representation of data: pb with overlapping name with autoplot
library(factoextra)
library(ggsignif)

fviz_eig(pca_res)

fviz_pca_ind(pca_res,
             axes = c(1, 2),
             col.ind = df_info$Prelevement, # Color by the quality of representation
             repel = TRUE,     # Avoid text overlapping
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Prelevements",
             labelsize = 2) 

fviz_pca_ind(pca_res,
             axes = c(1, 3),
             col.ind = df_info$Prelevement, # Color by the quality of representation
             repel = TRUE,     # Avoid text overlapping
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Prelevements",
             labelsize = 2)


fviz_pca_ind(pca_res,
             axes = c(2, 3),
             col.ind = df_info$Prelevement, # Color by the quality of representation
             repel = TRUE,     # Avoid text overlapping
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Prelevements",
             labelsize = 2)



# Graph of variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
###


#####################################################









########################################################################################################################
############################################## ZONE DE TEST ############################################################
########################################################################################################################















######################################################################################################
#############                  NORMALIZATION WITH DESEQ2 PACKAGE                         #############
######################################################################################################

# R 4.2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")




library("BiocManager")
BiocManager::install(pkgs = "DESeq2") # Oui A tout
BiocManager::install(pkgs = "XML") # NON a TOUT
BiocManager::install(pkgs = "locfit") # "a" puis oui


library("locfit")
library("XML")
library("DESeq2")


################
# Loading data #
################

coldata <- data.frame(read.table("metadata_modif.txt",  sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE))
coldata$info = str_c(coldata$Prelevement, '', coldata$Calf, '_', coldata$Run)

rownames(coldata) <- coldata$info
coldata <- coldata[order(rownames(coldata)),] 
coldata$Ind <- NULL
rownames(coldata)
coldata
    #             Condition
    #  J13-LU-TL1       umbo
    #  J13-LU-TL2       umbo
    #  J13-LU-TL3       umbo
    #  J22-LO-TL1a   oeillee
    #  J22-LO-TL1b   oeillee
    #  J22-LO-TL1c   oeillee
    #  J4-LD1              D
    #  J4-LD2              D
    #  J4-LD3              D
    #  J8-LV-TL1     velyger
    #  J8-LV-TL2     velyger
    #  J8-LV-TL3     velyger



cts <- read.table("6_taxo_affi_reads_species_modif.tsv",sep = "\t", header=T, stringsAsFactors = F, check.names=FALSE)

# set first column as rownames
rownames(cts) = cts[,1]
cts = cts[,-1]
cts = cts[,-1]

# Select specific columns "reads":
tmp = cts %>% select(starts_with('reads')) 

colnames(tmp) = gsub("reads_", "", colnames(tmp))

# remove unclassified reads
which(rownames(tmp)=="unclassified")
tmp2 = tmp[-which(rownames(tmp)=="unclassified"),]
which(rownames(tmp2)=="unclassified")
which(rownames(tmp2)=="cannot be assigned to a (non-viral) species")
tmp3 = tmp2[-which(rownames(tmp2)=="cannot be assigned to a (non-viral) species"),]
which(rownames(tmp3)=="cannot be assigned to a (non-viral) species")

cts = tmp3

colnames(cts)
cts <- cts[,order(colnames(cts))] 
cts <- cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(as.numeric(colnames(cts)))]
rownames(coldata)
colnames(cts)
head(cts)

###########
## CHECK ##
###########

# check : verifie si le nombre de colonnes de design est le meme que celui de join_devlarve.txt
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
coldata$Prelevement <- as.factor(coldata$Prelevement) 
str(cts)



##############
### DEseq2 ###
##############

# DESeqDataSet is a subclass of RangedSummarizedExperiment, used to store the input values, intermediate calculations and results of an analysis of differential expression
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Prelevement)
# Name_filtered <-rownames(dds)
colSums(cts)

# filtering
dds <- estimateSizeFactors(dds) # this function estimates the size factors using the "median ratio method" described by Equation 5 in Anders and Huber (2010)
dds$sizeFactor



#################################
### Normalisation des donnees ###
#################################

# Pour chaque colonne, v?rifie si le nombre de reads est >=10 + il faut que cela soit vrai au moins 3fois (par ligne).
idx <- rowSums(counts(dds,normalized=TRUE) >= 10 ) >= 3   
dds <- dds[idx,]                                          # Filtre sur les valeurs pr?c?dentes
dim(dds)

# Matrix log
norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 2)
#write.table(log.norm.counts,file="log_matrix.txt",quote=F)
norm.counts<-as.data.frame(norm.counts)



















##################################
# Define best number of clusters #
##################################

# Bayesian approach
#d_clust <- Mclust(as.matrix(logncsubset_minRowMeans), G=1:15, modelNames = mclust.options("emModelNames"))
#d_clust$BIC
#plot(d_clust)












