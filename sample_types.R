
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

# sample types for sorting the phylum-lvl. stacked barplot

# Arabidopsis samples
at_rel.rhizo.samples <- mapping_subset$SampleID[mapping_subset$host_description=="arabidopsis_and_rel" &
                                                      mapping_subset$compartment=="rhizosphere" &
                                                      mapping_subset$study!="Lundberg et al."]
at_rel.root.samples <- mapping_subset$SampleID[mapping_subset$host_description=="arabidopsis_and_rel" &
                                               mapping_subset$compartment!="rhizosphere" &
                                               mapping_subset$study!="Lundberg et al."]
arabidopsis.root.samples <- mapping_subset$SampleID[mapping_subset$host_description=="arabidopsis_and_rel" &
                                                    mapping_subset$compartment!="rhizosphere" &
                                                    mapping_subset$study=="Lundberg et al."]

# Barley samples
barley.rhizo.samples <- mapping_subset$SampleID[mapping_subset$host_description=="barley" &
                                                mapping_subset$compartment=="rhizosphere"]
barley.root.samples <- mapping_subset$SampleID[mapping_subset$host_description=="barley" &
                                               mapping_subset$compartment!="rhizosphere"]

# Maize samples
maize.rhizo.samples <- mapping_subset$SampleID[mapping_subset$host_description=="maize"]

# Rice samples
rice.rhizo.samples <- mapping_subset$SampleID[mapping_subset$host_description=="rice" &
                                              mapping_subset$compartment=="rhizosphere"]
rice.root.samples <- mapping_subset$SampleID[mapping_subset$host_description=="rice" &
                                             mapping_subset$compartment!="rhizosphere"]

# Grapevine samples
grapevine.rhizo.samples <- mapping_subset$SampleID[mapping_subset$host_description=="grapevine" &
                                                   mapping_subset$compartment=="rhizosphere"]
grapevine.root.samples <- mapping_subset$SampleID[mapping_subset$host_description=="grapevine" &
                                                  mapping_subset$compartment!="rhizosphere"]

# Fish samples
zebrafish.samples <- mapping_subset$SampleID[mapping_subset$host_common_name=="Zebrafish"]
stickleback.samples <- mapping_subset$SampleID[mapping_subset$host_common_name=="Stickleback"]

# Hydra samples
hydra.samples <- mapping_subset$SampleID[mapping_subset$host_description=="hydra"]

# Wild animal samples
wild_animal.samples <- mapping_subset$SampleID[mapping_subset$host_description=="wild_animal"]

#Human samples
human.infant.samples <- mapping_subset$SampleID[mapping_subset$host_description=="human" &
                                                mapping_subset$study=="Koenig et al."]
human.mit.samples <- mapping_subset$SampleID[mapping_subset$host_description=="human" &
                                             mapping_subset$study=="David et al."]
human.enterotypes.samples <- mapping_subset$SampleID[mapping_subset$host_description=="human" &
                                                     mapping_subset$study=="Wu et al."]
human.twins.samples <- mapping_subset$SampleID[mapping_subset$host_description=="human" &
                                               mapping_subset$study=="Goodrich et al."]

order.samples <- c(as.character(barley.rhizo.samples), as.character(barley.root.samples),
                   as.character(at_rel.rhizo.samples), as.character(at_rel.root.samples),
                   as.character(arabidopsis.root.samples),
                   as.character(maize.rhizo.samples),
                   as.character(rice.rhizo.samples), as.character(rice.root.samples),
                   as.character(grapevine.rhizo.samples), as.character(grapevine.root.samples),
                   as.character(zebrafish.samples), as.character(stickleback.samples),
                   as.character(hydra.samples),
                   as.character(wild_animal.samples),
                   as.character(human.infant.samples), as.character(human.mit.samples),
                   as.character(human.enterotypes.samples), as.character(human.twins.samples))

# general sample types (for MS numbers)

mammal.gut.samples <- mapping_subset$SampleID[mapping_subset$host_description=="human" |
                                              mapping_subset$host_description=="wild_animal"]

plant.root.samples <- mapping_subset$SampleID[mapping_subset$host_kingdom=="plant"]

fish.gut.samples <- mapping_subset$SampleID[mapping_subset$host_description=="fish"]

human.gut.samples <- mapping_subset$SampleID[mapping_subset$host_description=="human"]

