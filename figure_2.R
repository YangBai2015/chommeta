
# script to reproduce the analyses of alpha- and beta-diverisity
# and corresponding figures (Fig. 2A and 2B) from Hacquard et al., 2015
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

library(ggplot2)

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

shapes <- data.frame(group=c("root", "feces", "hydra_bodily_fluid", "hydra_polyp",
                             "rhizoplane", "rhizosphere", "endosphere"),
                     shape=c(3, 16, 18, 18,
                             3, 17, 3))

# load metadata

mapping.file <- "./data/metadata.txt"
mapping <- read.table <- read.table(mapping.file, sep="\t", header=T)
mapping_subset <- mapping[mapping$include, ]

### beta diversity (Fig. 2A)

# load unweighted unifrac distance 

uw_unifrac.file <- "./data/unweighted_unifrac_dm.txt"
uw_unifrac <- read.table <- read.table(uw_unifrac.file, sep="\t", header=T, check.names=F)
idx <- colnames(uw_unifrac) %in% mapping$SampleID
uw_unifrac_subset <- uw_unifrac[idx, idx]

# PCoA

k <- 2
pcoa <- cmdscale(uw_unifrac_subset, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")
points <- cbind(points, mapping[match(rownames(points), mapping$SampleID), ])

points$compartment <- factor(points$compartment, levels=shapes$group)
points$host_description <- factor(points$host_description, levels=colors$group)
points$y <- -points$y

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=host_description, shape=compartment)) +
     geom_point(alpha=.8, size=1.5) +
     scale_colour_manual(values=as.character(colors$color)) +
     scale_shape_manual(values=shapes$shape) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
          y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
     main_theme +
     theme(legend.position="none")

ggsave("./figures/PCoA.pdf", p)

### alpha diversity (Fig. 2B)

# note: in the MS figure, plant samples are sepparated
# by compartment (root / rhizosphere)

PD_whole_tree.file <- "./data/PD_whole_tree.txt"

alpha <- read.table(PD_whole_tree.file, sep="\t", header=T)
alpha <- merge(alpha, mapping, by="SampleID")
alpha <- alpha[alpha$include, ]

order <- sapply(levels(alpha$host_description),
                function(i) median(alpha$PD_whole_tree[alpha$host_description==i],
                na.rm=T))
order <- order[!is.na(order)]
order <- names(sort(order))
alpha$host_description <- factor(alpha$host_description, levels=order)

cols <- as.character(colors$color[match(levels(alpha$host_description), colors$group)])

p <- ggplot(alpha, aes(x=host_description, y=PD_whole_tree, color=host_description)) +
            geom_boxplot(alpha=1, outlier.size=.5, size=1) +
            scale_colour_manual(values=cols) +
            main_theme +
            theme(legend.position="none")

ggsave("./figures/PD.pdf", p)

# alpha .v beta diversity comparison

beta <- mapping_subset[match(colnames(uw_unifrac), mapping_subset$SampleID), ]
beta$mean_dist_intergroup <- NA

# avg. distance from each sample to every other in the same group

for (i in 1:dim(uw_unifrac)[1]) {

    sample_subset <- beta$SampleID[beta$host_description==beta$host_description[i]]
    d <- uw_unifrac[i, ]
    beta$mean_dist_intergroup[i] <- mean(as.numeric(d[names(d) %in% sample_subset]))

}

beta <- beta[, c(1, 8)]
comparison <- merge(alpha, beta, by="SampleID")

order <- c("arabidopsis_and_rel", "barley", "fish", "grapevine",
           "human", "hydra", "maize", "rice", "wild_animal")
comparison$host_description <- factor(comparison$host_description, levels=order)
comparison$compartment <- factor(comparison$compartment, levels=shapes$group)

cols <- as.character(colors$color[match(levels(comparison$host_description), colors$group)])

# plot alpha-diversity v. inter-group uw. UniFrac distance

p <- ggplot(comparison, aes(x=PD_whole_tree, y=mean_dist_intergroup,
                            color=host_description, shape=compartment)) +
     geom_point(alpha=.8, size=1.5) +
     scale_colour_manual(values=cols) +
     scale_shape_manual(values=shapes$shape) +
     labs(x="alpha diversity (PD)",
          y="mean unweighted UniFrac distance (within group)") + 
     main_theme +
     theme(legend.position="none")

ggsave("./figures/alpha_beta_comparison.pdf", p)

