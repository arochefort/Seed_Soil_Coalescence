library(phyloseq)
library(ggplot2)
library("vegan")

dITS <- readRDS("PS_ITS1.rds")
dgyrB <- readRDS("PS_gyrB.rds")

plantgyrB <- subset_samples(dgyrB, Compartment != "Seed" & Compartment != "Soil")
plantgyrB <- prune_taxa(taxa_sums(plantgyrB) > 0, plantgyrB)
plantgyrB
plantITS <- subset_samples(dITS, Compartment != "Seed" & Compartment != "Soil")
plantITS <- prune_taxa(taxa_sums(plantITS) > 0, plantITS)
plantITS

#add new columns in sample_data with pH NO3 and NH4
sample_data_gyrB <- read.csv("sample_data_gyrB.csv", sep=";", dec=".", check.names=FALSE, row.names=1)
sample_data_ITS1 <- read.csv("sample_data_ITS1.csv", sep=";", dec=".", check.names=FALSE, row.names=1)

sample_data_gyrB <- sample_data(sample_data_gyrB)
plantgyrB = merge_phyloseq(plantgyrB, sample_data_gyrB)

sample_data_ITS1 <- sample_data(sample_data_ITS1)
plantITS = merge_phyloseq(plantITS, sample_data_ITS1)

#Subset by plant compartments
rootgyrB <- subset_samples(plantgyrB, Compartment== "Root")
rootgyrB <- filter_taxa(rootgyrB, function(x) sum(x) > 0, TRUE)
stemgyrB <- subset_samples(plantgyrB, Compartment== "Stem")
stemgyrB <- filter_taxa(stemgyrB, function(x) sum(x) > 0, TRUE)
rootITS <- subset_samples(plantITS, Compartment== "Root")
rootITS <- filter_taxa(rootITS, function(x) sum(x) > 0, TRUE)
stemITS <- subset_samples(plantITS, Compartment== "Stem")
stemITS <- filter_taxa(stemITS, function(x) sum(x) > 0, TRUE)

#log transform
root_log_gyrB <- transform_sample_counts(rootgyrB, function(x) log(1 + x))
stem_log_gyrB <- transform_sample_counts(stemgyrB, function(x) log(1 + x))
root_log_ITS1 <-  transform_sample_counts(rootITS, function(x) log(1 + x))
stem_log_ITS1 <-  transform_sample_counts(stemITS, function(x) log(1 + x))

#calculate wUF and BC distances
root.wUF.gyrB <- ordinate(root_log_gyrB, "PCoA", "wunifrac")
stem.wUF.gyrB <- ordinate(stem_log_gyrB, "PCoA", "wunifrac")
root.bc.ITS <- ordinate(root_log_ITS1, "PCoA", "bray")
stem.bc.ITS <- ordinate(stem_log_ITS1, "PCoA", "bray")

#Stat by habitat with richness
metadata.stem.gyrB <- as(sample_data(stem_log_gyrB), "data.frame") #Convert sample_data to data.frame
dist.wUF.stem.gyrB <- phyloseq::distance(stem_log_gyrB, method = "wunifrac")
adonis.wUF.stem.gyrB <- adonis2(dist.wUF.stem.gyrB ~ Fung_rich + pH + NO3 + NH4 + Combined + Stage + Bact_rich,
                                data = metadata.stem.gyrB, perm = 999) 
print(adonis.wUF.stem.gyrB)
metadata.root.gyrB <- as(sample_data(root_log_gyrB), "data.frame") #Convert sample_data to data.frame
dist.wUF.root.gyrB <- phyloseq::distance(root_log_gyrB, method = "wunifrac")
adonis.wUF.root.gyrB <- adonis2(dist.wUF.root.gyrB ~ Bact_rich + Fung_rich + pH + NO3 + NH4 + Combined + Stage,
                                data = metadata.root.gyrB, perm = 999) 
print(adonis.wUF.root.gyrB)

metadata.stem.ITS <- as(sample_data(stem_log_ITS1), "data.frame") #Convert sample_data to data.frame
dist.bc.stem.ITS <- phyloseq::distance(stem_log_ITS1, method = "bray")
adonis.BC.stem.ITS <- adonis2(dist.bc.stem.ITS ~ Bact_rich + Fung_rich + pH + NO3 + NH4 + Combined + Stage,
                              data = metadata.stem.ITS, perm = 999)
print(adonis.BC.stem.ITS)

metadata.root.ITS <- as(sample_data(root_log_ITS1), "data.frame") #Convert sample_data to data.frame
dist.bc.root.ITS <- phyloseq::distance(root_log_ITS1, method = "bray")
adonis.BC.root.ITS <- adonis2(dist.bc.root.ITS ~ Bact_rich + Fung_rich + pH + NO3 + NH4 + Combined + Stage,
                              data = metadata.root.ITS, perm = 999)
print(adonis.BC.root.ITS)


##with Faith's PD 
adonis.wUF.stem.gyrB_PD <- adonis2(dist.wUF.stem.gyrB ~ Bact_PD + Fung_PD + pH + NO3 + NH4 + Combined + Stage,
                                data = metadata.stem.gyrB, perm = 999) 
print(adonis.wUF.stem.gyrB_PD)

adonis.wUF.root.gyrB_PD <- adonis2(dist.wUF.root.gyrB ~ Bact_PD + Fung_PD + pH + NO3 + NH4 + Combined + Stage,
                                data = metadata.root.gyrB, perm = 999) 
print(adonis.wUF.root.gyrB_PD)

adonis.BC.stem.ITS_PD <- adonis2(dist.bc.stem.ITS ~ Bact_PD + Fung_PD + pH + NO3 + NH4 + Combined + Stage,
                              data = metadata.stem.ITS, perm = 999)
print(adonis.BC.stem.ITS_PD)

adonis.BC.root.ITS_PD <- adonis2(dist.bc.root.ITS ~ Bact_PD + Fung_PD + pH + NO3 + NH4 + Combined + Stage,
                              data = metadata.root.ITS, perm = 999)
print(adonis.BC.root.ITS_PD)
