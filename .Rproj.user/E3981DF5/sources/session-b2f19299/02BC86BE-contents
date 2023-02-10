# Tue Dec 13 12:30:44 2022 ------------------------------
# Birds complete

# Packages
library(picante)
library(RRphylo)
library(Rphylopars)
library(phytools)
library(ape)
library(letsR)
library(raster)
library(visdat)
library(tidyverse)
library(ggExtra)
library(cowplot)
library(geiger)
library(rgdal)

# Traits
setwd("C:/Users/Cliente/Desktop/vertebrates_evolution/traits/birds")
birds_traits <- as_tibble(readxl::read_xlsx("AVONET Supplementary dataset 1.xlsx", sheet=4))

short <- birds_traits %>% select(Species3, Mass)
short$Species3 <- gsub(" ", "_",short$Species3)
short

# Phylogeny
setwd("C:/Users/Cliente/Desktop/vertebrates_evolution/phylogeny")
birds_phy <- read.nexus("birds_consensus.tre")

# Transposition
birds_trans <- t(short)
colnames(birds_trans) <- birds_trans[1,]
View(birds_trans)++++++

#--- Match with phylogeny
birds_short_phy <- prune.sample(birds_trans, birds_phy)
birds_short_phy$edge.length

# Evolution rate
# Preparing the data
# Anura
mass_birds <- as.numeric(short$Mass)
names(mass_birds) <- short$Species3
mass_birds

# body size
rates.birds.mass <- RRphylo(tree= birds_short_phy, y= mass_birds)
memory.limit(size=number)

gc()
memory.limit (9999999999)
