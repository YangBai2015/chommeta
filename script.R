
# script to reproduce the analysis and figures from Hacquard et al., 2015
#
# if you use any of the following code, please cite:
#
# Stephane Hacquard, Ruben Garrido-Oter, Antonio González Peña, Stijn Spaepen,
# Gail Ackermann, Sarah Lebeis, Alice C. McHardy, Jeffrey L. Dangl, Rob Knight,
# Ruth Ley Paul Schulze-Lefert. Microbiota and Host Nutrition: Structure,
# Diversification and Functions across Plant and Animal Kingdoms,
# Cell Host and Microbe, 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# plot figures

system("mkdir -p ./figures")
system("rm -f ./figures/*pdf")

source("figure_2.R")
source("figure_4.R")
source("figure_5.R")
source("figure_s1.R")

# calculate numbers for the MS

source("sample_types.R")

# Proteobacteria abundances in plant root and fish and human gut

phylum_table <- aggregate(otu_table_norm, by=list(taxonomy$phylum), FUN=sum) * 100
rownames(phylum_table) <- levels(taxonomy$phylum)

print("perc. Proteobacteria in plant root samples")
print(mean(as.numeric(phylum_table[rownames(phylum_table)=="Proteobacteria",
                                   colnames(phylum_table) %in% plant.root.samples])))

print("perc. Proteobacteria in plant fish gut samples")
print(mean(as.numeric(phylum_table[rownames(phylum_table)=="Proteobacteria",
                      colnames(phylum_table) %in% fish.gut.samples])))

print("perc. Proteobacteria in plant human gut samples")
print(mean(as.numeric(phylum_table[rownames(phylum_table)=="Proteobacteria",
                      colnames(phylum_table) %in% human.gut.samples])))

# Bacteroidales and Clostridiales abundances in mammal gut and plant root

order_table <- aggregate(otu_table_norm, by=list(taxonomy$order), FUN=sum) * 100
rownames(order_table) <- levels(taxonomy$order)

print("perc. Bacteroidales in human and wild aninal samples")
print(mean(as.numeric(order_table[rownames(order_table)=="Bacteroidales",
                                  colnames(order_table) %in% mammal.gut.samples])))

print("perc. Clostridiales in human and wild animal samples")
print(mean(as.numeric(order_table[rownames(order_table)=="Clostridiales",
                                  colnames(order_table) %in% mammal.gut.samples])))

print("perc. Bacteroidales in plant root samples")
print(mean(as.numeric(order_table[rownames(order_table)=="Bacteroidales",
                                  colnames(order_table) %in% plant.root.samples])))

print("perc. Clostridiales in plant root samples")
print(mean(as.numeric(order_table[rownames(order_table)=="Clostridiales",
                                  colnames(order_table) %in% plant.root.samples])))

