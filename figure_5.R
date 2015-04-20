
# script to reproduce the tree figure of rel. abundances of 
# abundant community members (Fig. 5) from Hacquard et al., 2015
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

# plotting stuff

source("plotting.functions.R")

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    legend.position="top",
                    text=element_text(family="sans"))

colors <- data.frame(group=c("arabidopsis_and_rel", "barley", "fish",
                             "grapevine", "human", "hydra",
                             "maize", "rice", "wild_animal"),
                     color=c("#f8766d", "brown", "#00b9e3", 
                             "#00ba38", "#619cff", "#000000",
                             "#93aa00", "#00c19f", "#db72fb"))

# load metadata

mapping.file <- "./data/metadata.txt"
mapping <- read.table <- read.table(mapping.file, sep="\t", header=T)
mapping_subset <- mapping[mapping$include, ]

# load OTU table

otu_table.file <- "./data/otu_table.txt"

tab5rows <- read.table(otu_table.file, sep="\t", header=T, nrows=5, check.names=F)
classes <- sapply(tab5rows, class)
otu_table <- read.table(otu_table.file, sep="\t", header=T, colClasses=classes, check.names=F)

# parse taxonomy info

taxonomy <- unlist(sapply(otu_table$taxonomy, function(x) strsplit(as.character(x), "; ")))
taxonomy <- gsub("^.__", "", taxonomy)
taxonomy <- t(matrix(taxonomy, 7, dim(otu_table)[1]))
colnames(taxonomy) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxonomy <- data.frame(OTUID=rownames(otu_table), taxonomy)

# subset table and normalize

otu_table_subset <- otu_table[, colnames(otu_table) %in% mapping_subset$SampleID]
otu_table_norm <- apply(otu_table_subset, 2, function(x) x / sum(x))

# calculate average across group

mean_abundances <- data.frame(OTUID=rownames(otu_table_norm))

for (i in 1:length(levels(mapping_subset$host_description))){

    mean_abundances[, i] <- 0

    group <- levels(mapping_subset$host_description)[i]
    group_samples <- mapping_subset$SampleID[mapping_subset$host_description==group]
    group_samples <- group_samples[group_samples %in% colnames(otu_table_norm)]
    mean_abundances[, i] <- apply(otu_table_norm[, colnames(otu_table_norm) %in% group_samples], 1, mean)
    
    colnames(mean_abundances)[i] <- group 

}

mean_abundances <- mean_abundances[, !colnames(mean_abundances) %in% c("freshwater", "soil")]
rownames(mean_abundances) <- rownames(otu_table_norm)

# get abundant community members:
# OTUs with a mean higher than thresold in at least one group

threshold <- .1 / 100

ACM <- rownames(mean_abundances[rowSums(mean_abundances > threshold)!=0, ])

taxonomy_ACM <- taxonomy[taxonomy$OTUID %in% ACM, ]

ACM <- taxonomy_ACM$OTUID[taxonomy_ACM$kingdom=="Bacteria"]

write.table(ACM, "./data/ACM.txt",  quote=F, col.names=F, row.names=F)

# extract rep seqs from GreenGenes DB
# note: the aligned rep. seqs. and tree are provided in the data folder
# to recalculate uncomment the next lines
 
# system("rm -f ./data/ACM_aligned.fasta")
# system("./get_sequences.sh ./data/ACM.txt ./data/97_otus.fasta >> ./data/ACM_aligned.fasta")
# 
# # build tree
# 
# system("rm -f ./data/ACM_aligned.tree")
# system("FastTree -nt -gtr ACM_aligned.fasta >> ./data/ACM_aligned.tree")

# load tree and generate feature vectors for plotting

tree <- read.tree("./data/ACM_aligned.tree")

ntips <- length(tree$tip.label)

tips.cex <- rep(.3, ntips)

tips.shape <- rep(19, ntips)

phyla <- data.frame(group=c("Chloroflexi", "Acidobacteria", "Actinobacteria",
                            "Bacteroidetes", "Firmicutes", "Alphaproteobacteria",
                            "Betaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria",
                            "Other"),
                    color=c("black", c_red, c_yellow,
                            c_blue, c_orange, c_green,
                            c_sea_green, c_dark_green, c_very_dark_green,
                            "grey"))

idx <- taxonomy$phylum == "Proteobacteria"
taxonomy_group <- as.character(taxonomy$phylum)
taxonomy_group[idx] <- as.character(taxonomy$class[idx])
taxonomy_group[!taxonomy_group %in% phyla$group] <- "Other"
taxonomy$group <- factor(taxonomy_group, levels=c("Acidobacteria", "Actinobacteria", "Alphaproteobacteria",
                                                  "Bacteroidetes", "Betaproteobacteria", "Chloroflexi",
                                                  "Deltaproteobacteria", "Firmicutes", "Gammaproteobacteria",
                                                  "Other"))

taxonomy_ACM <- taxonomy[taxonomy$OTUID %in% ACM, ]
taxonomy_ACM <- taxonomy_ACM[match(tree$tip.label, taxonomy_ACM$OTUID), ]

tips.col <- as.character(phyla$color[match(taxonomy_ACM$group, phyla$group)])

mean_abundances <- mean_abundances[, c("barley", "arabidopsis_and_rel", "maize", "rice",
                                       "grapevine", "fish", "hydra", "wild_animal", "human")]

scaling.factor <- 1e4
barplots.height <- apply(mean_abundances, 2, function(x) log2(1 + x * scaling.factor))
barplots.height <- barplots.height[match(tree$tip.label, rownames(mean_abundances)), ]
barplots.height <- barplots.height / 400

barplots.col <- t(matrix(rep(colors$color[match(colnames(mean_abundances), colors$group)], dim(barplots.height)[1]), 
                         ncol=dim(barplots.height)[1], nrow=dim(barplots.height)[2]))

barplots.height <- cbind(.01, barplots.height)
barplots.col <- cbind(tips.col, barplots.col)

tree$tip.label <- gsub("\\]", "", gsub("\\[", "", taxonomy_ACM$order))

plot_btree(tree=tree, output.file="./figures/tree.pdf", type="phylogram", edge.width=.5, 
           color.labels="black", align.labels=T, label.offset=.01, cex.labels=.05, color.alg.lines="transparent",
           tips.col=tips.col, tips.shape=tips.shape, tips.cex=tips.cex,
           barplots.height=barplots.height, barplots.width=1, barplots.col=barplots.col)

