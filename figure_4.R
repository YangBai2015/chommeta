
# script to reproduce the stacked barplot of phylum-level relative abundances
# per sample (Fig. 4) from Hacquard et al., 2015
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

phyla <- data.frame(group=c("Chloroflexi", "Acidobacteria", "Actinobacteria",
                            "Bacteroidetes", "Firmicutes", "Alphaproteobacteria",
                            "Betaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria",
                            "Other"),
                    color=c("black", c_red, c_yellow,
                            c_blue, c_orange, c_green,
                            c_sea_green, c_dark_green, c_very_dark_green,
                            "white"))

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

idx <- taxonomy$phylum == "Proteobacteria"
taxonomy_group <- as.character(taxonomy$phylum)
taxonomy_group[idx] <- as.character(taxonomy$class[idx])
taxonomy_group[!taxonomy_group %in% phyla$group] <- "Other"
taxonomy$group <- factor(taxonomy_group,
                         levels=c("Acidobacteria", "Actinobacteria", "Alphaproteobacteria",
                                   "Bacteroidetes", "Betaproteobacteria", "Chloroflexi",
                                   "Deltaproteobacteria", "Firmicutes", "Gammaproteobacteria",
                                   "Other"))

# subset table and normalize

otu_table_subset <- otu_table[, colnames(otu_table) %in% mapping_subset$SampleID]
otu_table_norm <- apply(otu_table_subset, 2, function(x) x / sum(x))

### stacked barplot of phylum abundances (Fig. 4)
 
phylum_table <- aggregate(otu_table_norm, by=list(taxonomy$group), FUN=sum)
rownames(phylum_table) <- levels(taxonomy$group)

df <- melt(phylum_table)
colnames(df) <- c("group", "SampleID", "abundance")

df$host_description <- mapping$host_description[match(df$SampleID, mapping$SampleID)]
df$host_description <- factor(df$host_description,
                              levels=c("barley", "arabidopsis_and_rel", "maize",
                                       "rice", "grapevine", "fish", "hydra",
                                       "wild_animal", "human"))

# sort samples according to their type

source("sample_types.R")

order.samples <- order.samples[order.samples %in% df$SampleID]

idx <- match(order.samples, levels(df$SampleID))
df$SampleID <- factor(df$SampleID, levels=levels(df$SampleID)[idx])

cols <- as.character(phyla$color[match(levels(df$group), phyla$group)])

# plot phylum stacked barplots per sample

p <- ggplot(df, aes(x=SampleID, y=abundance, fill=group)) +
            geom_bar(stat='identity') +
            scale_fill_manual(values=cols) +
            main_theme +
            theme(axis.ticks.x=element_blank(),
                  axis.text.x=element_blank(),
                  legend.position="none") 

ggsave("./figures/phylum_barplot.pdf", p, width=7, height=5)

