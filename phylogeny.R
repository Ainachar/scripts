library("phyloseq")
library("ape")
library("ggplot2")
library("vegan")


#The OTU table:
abund_table<-
  read.csv("~/ec34/biosin5410/microbiome/data/All_Good_P2_C03.csv",row.names=1,check.names=FALSE
  )
abund_table

#Transpose (flip) the data to have sample names on rows and OTU on column
abund_table<-t(abund_table)

#Now load the taxonomy and metainformation

OTU_taxonomy<-
  read.csv("~/ec34/biosin5410/microbiome/data/All_Good_P2_C03_Taxonomy.csv",row.names=1,check.names=FALSE)

meta_table<-
  read.csv("~/ec34/biosin5410/microbiome/data/ENV_pitlatrine.csv",row.names=1,check.names=FALSE)

#Get grouping information from the meta table. See that functions can be used inside functions.

grouping_info<-
  data.frame(row.names=rownames(meta_table),t(as.data.frame(strsplit(rownames(
    meta_table),"_"))))

#Name the columns:
colnames(grouping_info)<-c("Country","Latrine","Depth")
meta_table<-data.frame(meta_table,grouping_info)

#Filter out samples not present in meta_table. 
#NB! must be equally many rows for both meta data and OTU count

abund_table<-abund_table[rownames(abund_table) %in% rownames(meta_table),]

#Read in the phylogenetic tree file. need this. uses package ape to read this in phyloseq

OTU_tree <- read_tree("~/ec34/biosin5410/microbiome/data/All_Good_P2_C03.tre")

abund_table
OTU_taxonomy
meta_table

str(abund_table)
str(OTU_taxonomy)     
head(meta_table)
head(OTU_taxonomy)
head(abund_table)

dim(abund_table)
dim(OTU_taxonomy)

grouping_info
head(grouping_info)

rowMeans(abund_table)
colMeans(abund_table)
max.col(abund_table)
rowSums(abund_table)
max(rowSums(abund_table))
min(rowSums(abund_table))
rowMeans(abund_table)
max(colMeans(abund_table))

#what is C9 OTU, including taxa, abundance mean, min, max etc. C9 is in column 10 

mean(abund_table[,10])
max(abund_table[,10])
min(abund_table[,10])
head(OTU_taxonomy, 10:10)
OTU_taxonomy[10,]

#Convert the data to phyloseq format. This will become an phyloseq spesific object.

OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE) #only numbers is ok
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,OTU_tree)  #merge all three files into one object


#Phyloseq has functions for estimating richness or alpha diversity. estimate_richness will make a data frame of all the richness indexes. 
richness<-estimate_richness(physeq, split = TRUE, measures = NULL)

#Visualise the richness (alpha diversity) with the plot_richness function.
#The plot:
  p<-plot_richness(physeq, x = "Country", color = "Depth", measures=c("Observed","Shannon","InvSimpson"))

  p
#Make a nice background

  p<-p+theme_bw()
p

#Aim: To calculate and visualise community composition
#Calculate
distances_bray<-distance(physeq, method="bray")
distances_unifrac<-distance(physeq, method="unifrac")

str(distances_bray)
distances_unifrac

#Make an PCoA (https://en.wikipedia.org/wiki/Multidimensional_scaling#Types )
 #             ordination plot based on abundance based unifrac distances with the following commands
              ord <- ordinate(physeq, method="PCoA", distance="bray", weighted=TRUE)
              p <- plot_ordination(physeq, ord, color="Country",title="Phyloseq's Weighted Unifrac")
              p <- p + geom_point(size=5) + theme_bw()
p

ord <- ordinate(physeq, method="PCoA", distance="unifrac", weighted=TRUE)
p2 <- plot_ordination(physeq, ord, color="Country",title="Phyloseq's Weighted Unifrac")
p2 <- p2 + geom_point(size=5) + theme_bw()
p2
