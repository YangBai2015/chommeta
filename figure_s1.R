
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

### analysis of dispersion per host / study (supplementary Fig. S1A)

source("sample_types.R")

# load unweighted unifrac distance 

uw_unifrac.file <- "./data/unweighted_unifrac_dm.txt"
uw_unifrac <- read.table <- read.table(uw_unifrac.file, sep="\t", header=T, check.names=F)
idx <- colnames(uw_unifrac) %in% mapping$SampleID
uw_unifrac_subset <- uw_unifrac[idx, idx]

study_variance <- NULL

for (s in levels(mapping_subset$study)) {
    
    samples <- mapping_subset$SampleID[mapping_subset$study==s &
                                       mapping_subset$SampleID %in% order.samples]
    idx <- colnames(uw_unifrac_subset) %in% samples
    v <- as.vector(as.dist(uw_unifrac_subset[idx, idx]))
    h <- as.character(mapping_subset$host_description[mapping_subset$study==s][1])
    
    study_variance <- rbind(study_variance, data.frame(study=s, variance=v, host_description=h))

}

host_variance <- NULL

for (h in levels(colors$group)) {

    samples <- mapping_subset$SampleID[mapping_subset$host_description==h &
                                       mapping_subset$SampleID %in% order.samples]
    idx <- colnames(uw_unifrac_subset) %in% samples
    v <- as.vector(as.dist(uw_unifrac_subset[idx, idx]))

    host_variance <- rbind(host_variance, data.frame(variance=v, host_description=h))

}

host_variance$study <- gsub("$", " (all studies)", host_variance$host_description)

# plotting

variance <- rbind(study_variance, host_variance)

lvls <- c("Franzenburg et al.",
          "arabidopsis_and_rel (all studies)", "Schlaeppi et al.", "Lundberg et al.",
          "Bulgarelli et al.",
          "Zarraonaindia et al.",
          "Peiffer et al.",
          "Edwards et al.",
          "fish (all studies)", "Roeselers et al.", "Bolnick et al.",
          "human (all studies)", "David et al.", "Koenig et al.", "Wu et al.", "Goodrich et al.", 
          "Muegge et al.")
variance$study <- factor(variance$study, levels=lvls)
variance <- variance[!is.na(variance$study), ]

cols <- as.character(colors$color[match(levels(variance$host_description), colors$group)])

p <- ggplot(variance, aes(x=study, y=variance, color=host_description)) +
            geom_boxplot(alpha=1, outlier.size=0, size=.7) +
            scale_colour_manual(values=cols) +
            labs(x="study / group", y="within-group UniFrac pairwise distance") +
            main_theme +
            theme(legend.position="none",
                  axis.text.x=element_text(angle=45, hjust=1))

ggsave("./figures/Fig_S1A.pdf", p)

### supplementary figure S1B

k <- 2
pcoa <- cmdscale(uw_unifrac_subset, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")
points <- cbind(points, mapping[match(rownames(points), mapping$SampleID), ])

points$compartment <- factor(points$compartment, levels=shapes$group)
points$host_description <- factor(points$host_description, levels=colors$group)
points$study <- factor(points$study, levels=lvls)
points$y <- -points$y

# plotting

cols <- rep("grey", length(levels(points$study)))
cols[levels(points$study)=="Lundberg et al."] <- "red"
cols[levels(points$study)=="Bulgarelli et al."] <- "brown"
cols[levels(points$study)=="Roeselers et al."] <- "yellow"
cols[levels(points$study)=="Bolnick et al."] <- "green"
cols[levels(points$study)=="David et al."] <- "blue"
cols[levels(points$study)=="Koenig et al."] <- "pink"
cols[levels(points$study)=="Wu et al."] <- "purple"
cols[levels(points$study)=="Goodrich et al."] <- "darkblue"

p <- ggplot(points, aes(x=x, y=y, color=study, shape=compartment)) +
     geom_point(alpha=.8, size=1.5) +
     scale_colour_manual(values=cols) +
     scale_shape_manual(values=shapes$shape) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
          y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
     main_theme +
     theme(legend.position="left")

ggsave("./figures/Fig_S1B.pdf", p, height=7, width=10)

### CPCoA analysis for study and host factors (for revision)
# note: CPCoA for such a large dataset is very slow if run serially (commented code bellow)

### library(vegan)
### 
### variability_table <- function(cca){
### 
###         chi <- c(cca$tot.chi,
###                        cca$CCA$tot.chi, cca$CA$tot.chi)
###         variability_table <- cbind(chi, chi/chi[1])
###         colnames(variability_table) <- c("inertia", "proportion")
###         rownames(variability_table) <- c("total", "constrained", "unconstrained")
###         return(variability_table)
### 
### }
### 
### mapping <- mapping[match(rownames(uw_unifrac), mapping$SampleID), ]
### rownames(mapping) <- mapping$SampleID
### 
### # CPCoA constrained by host and study
### 
### cpcoa_study <- capscale(uw_unifrac ~ study + Condition(host_description), data=mapping)
### cpcoa_host <- capscale(uw_unifrac ~ host_description, data=mapping)
### 
### # ANOVA-like permutation analysis
### # anova.cca is a vegan wrapper for CCA, which uses the function permutest
### 
### perm_anova_host <- anova.cca(cpcoa_host)
### print(perm_anova_host)
### 
### perm_anova_study <- anova.cca(cpcoa_study)
### print(perm_anova_study)
### 
### # generate variability tables and calculate confidence intervals for the variance
### 
### var_tbl_host <- variability_table(cpcoa_host)
### print(var_tbl_host)
### 
### var_tbl_study <- variability_table(cpcoa_study)
### print(var_tbl_study)
### 
### # extract the weighted average (sample) scores
### 
### wa_host <- cpcoa_host$CCA$wa
### wa_study <- cpcoa_study$CCA$wa
### 
