---
title: "Bring diversity data together & modify files"
output: github_document
editor_options: 
  chunk_output_type: console
---
### Load packages

``` r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

``` r
library(parallel)
library(stringr)
library(reshape2)
library(phytools)
```

```
## Loading required package: ape
```

```
## 
## Attaching package: 'ape'
```

```
## The following object is masked from 'package:dplyr':
## 
##     where
```

```
## Loading required package: maps
```

``` r
library(tidyr)
```

```
## 
## Attaching package: 'tidyr'
```

```
## The following object is masked from 'package:reshape2':
## 
##     smiths
```

### Load in files

``` r
env_by_lake <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/environment/env_by_lake.csv")

# Fish present in the marine lakes and reference pool by lake code. fish_matched and fish_matched2
pres_abs_by_lake <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/fish_presence_matrix_by_lake.csv")

# Fish present in the marine lakes and reference pool by species. fish_matched and fish_matched2
pres_abs_by_species <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/fish_presence_matrix_by_species.csv")

# Fish present in marine lakes only
surveyed_sites_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/surveyed_site_presence_by_species.csv")

# Ocean site fish presence by species 
ocean_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/ocean_site_presence_by_species.csv")

# Holomictic fish presence by species 
holomictic_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/holomictic_presence_by_species.csv")

# Meromictic fish presence by species 
meromictic_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/meromictic_presence_by_species.csv")

# Mixed fish presence by species 
mixed_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/mixed_presence_by_species.csv")

# Stratified fish presence by species 
stratified_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/stratified_presence_by_species.csv")

# Surveyed sites fish presence by species each in separate files
BCM_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/BCM_fish.csv")
CLM_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/CLM_fish.csv")
FLK_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/FLK_fish.csv")
GLK_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/GLK_fish.csv")
HLM_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/HLM_fish.csv")
HLO_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/HLO_fish.csv")
IBK_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/IBK_fish.csv")
LLN_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/LLN_fish.csv")
LCN_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/LCN_fish.csv")
MLN_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/MLN_fish.csv")
NCN_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/NCN_fish.csv")
NLK_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/NLK_fish.csv")
NLN_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/NLN_fish.csv")
NLU_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/NLU_fish.csv")
OLO_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/OLO_fish.csv")
OCO_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/OCO_fish.csv")
OTM_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/OTM_fish.csv")
OCM_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/OCM_fish.csv")
RCA_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/RCA_fish.csv")
SLN_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/SLN_fish.csv")
TLN_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/TLN_fish.csv")
ULN_fish <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/ULN_fish.csv")

# Phylogenetic tree of species present in marine lakes, tree
tree <- read.tree("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/phylo/palau_fish_tree.tre")

# Traits for fish present in the marine lakes and reference pool, fish_traits_complete_matched and fish_traits_complete_matched2
species_functional_traits <- read.csv("/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/traits/final_traits_gnathostomata.csv")
```

## Check species names
- Check that names are consistent across files. 

``` r
## Check species names in incidence matrix, trait matrix, and phylogenetic tree
# Check trait data with occurrence data
# If differences are found then check meromictic, holomictic, and ocean site dataframes too
not_present_trait <- setdiff(species_functional_traits$Species, pres_abs_by_species$X)
not_present_trait
```

```
## character(0)
```

``` r
# Check tree data with occurrence data
# If differences are found then check meromictic, holomictic, and ocean site dataframes too
not_present_tree <- setdiff(tree$tip.label, pres_abs_by_species$X)
not_present_tree
```

```
## character(0)
```

``` r
# Check tree data with trait data
not_present <- setdiff(tree$tip.label, species_functional_traits$Species)
not_present
```

```
## character(0)
```

``` r
# Find what lakes pres/abs and env data
site_env_only <- setdiff(env_by_lake$X, pres_abs_by_lake$X)
site_env_only
```

```
## character(0)
```

``` r
# Delete lakes from above with no pres/abs data, LCM, LSM, SLM, ULU
# env_by_lake <- env_by_lake[-c(8, 11, 23, 27),]
```

## Modify environmental data

``` r
# Contains environmental data for lakes sampled from
env_by_lake$Location_type <- factor(env_by_lake$Location_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))

env_by_lake$logArea <- log(env_by_lake$surface_area_m2)

row.names(env_by_lake) <- env_by_lake$X
env <- env_by_lake
env$C <- env$temperature_median
env$S <- env$salinity_median
env$O <- env$oxygen_median
env$minD <- env$distance_to_ocean_min_m
env$M <- env$max_depth
env$LA <- env$logArea

scaled_env <- as.data.frame(scale(env[,c(2,6,8,10,22:32)]))

# Create a formula for aggregation
formula <- as.formula(paste(". ~ Location_type"))

# Use the aggregate function to calculate the mean for each variable by Location_type
Location_type_env <- aggregate(formula, data = env[, c(2,4,6,8,10,19,22:32)], FUN = mean)
```

## Modify incidence matrices

``` r
# Presence/absence by lake
row.names(pres_abs_by_lake) <- pres_abs_by_lake$X
pres_abs_by_lake[,c(2:1725)] <- sapply(pres_abs_by_lake[,c(2:1725)], as.numeric)
presabs_lake <- pres_abs_by_lake[,-1]

# Kept this here because HLM is holomictic but is mixed, so it differs between datasets #
com_presabs_lake <- presabs_lake
com_presabs_lake[24,] <- colSums(com_presabs_lake[c(7,9,11,16,18:19),])
rownames(com_presabs_lake)[rownames(com_presabs_lake) == "24"] <- "Ocean sites"
com_presabs_lake[25,] <- colSums(com_presabs_lake[c(3,5:6,8,10,13:15,23),])
rownames(com_presabs_lake)[rownames(com_presabs_lake) == "25"] <- "Holomictic lakes"
com_presabs_lake[26,] <- colSums(com_presabs_lake[c(1:2,4,12,17,21,22),])
rownames(com_presabs_lake)[rownames(com_presabs_lake) == "26"] <- "Meromictic lakes"

strat_presabs_lake <- presabs_lake
strat_presabs_lake[24,] <- colSums(strat_presabs_lake[c(7,9,11,16,18:19),])
rownames(strat_presabs_lake)[rownames(strat_presabs_lake) == "24"] <- "Ocean sites"
strat_presabs_lake[25,] <- colSums(strat_presabs_lake[c(3,6,8,10,13:15,23),])
rownames(strat_presabs_lake)[rownames(strat_presabs_lake) == "25"] <- "Mixed lakes"
strat_presabs_lake[26,] <- colSums(strat_presabs_lake[c(1:2,4,5,12,17,21,22),])
rownames(strat_presabs_lake)[rownames(strat_presabs_lake) == "26"] <- "Stratified lakes"

# Add total number of fish per site and how many environments each species is present in
pres_abs_lake_sums <- pres_abs_by_lake
pres_abs_lake_sums$row_sum <- rowSums(pres_abs_lake_sums[,-1])
pres_abs_lake_sums[24,-1] <- colSums(pres_abs_lake_sums[,-1])
site_sums <- pres_abs_lake_sums[, c(1,1726)]
species_sums <- pres_abs_lake_sums[c(24),-1]
rownames(species_sums)[rownames(species_sums) == "24"] <- "col_sum"
pres_abs_lake_sums <- pres_abs_lake_sums[,-1]

#Merge the species richness with env data
SR_env <- merge(site_sums, env, by = 'X', sort = F)
row.names(SR_env) <- SR_env$X
# data <- SR_env[,c(1:2,25)]
# row.names(data) <- data$X
# data <- data[,-1]
# colnames(data)[1] <- 's'
# colnames(data)[2] <- 'a'
# name <- "Palau fish data"
# SAR <- list(name, data)
# data(power);data(expo);data(negexpo);data(monod);data(ratio);data(logist);data(lomolino);data(weibull)
# mods <- c("power","expo","negexpo","monod","logist","ratio","lomolino","weibull")
# resAverage <- multiSAR(modelList = mods, SAR)


# Presence/absence by species
row.names(pres_abs_by_species) <- pres_abs_by_species$X
pres_abs_by_species[,c(2:24)] <- sapply(pres_abs_by_species[,c(2:24)], as.numeric)
presabs_species <- pres_abs_by_species[,-c(1)]

pres_abs_species_sums <- pres_abs_by_species
pres_abs_species_sums$row_sum <- rowSums(pres_abs_species_sums[,-1])
pres_abs_species_sums <- pres_abs_species_sums[,-1]


# Surveyed Locations
surveyed_sites_fish[,c(2:23)] <- sapply(surveyed_sites_fish[,c(2:23)], as.numeric)
row.names(surveyed_sites_fish) <- surveyed_sites_fish$X
surveyed_sites_species <- surveyed_sites_fish[,-1]
# Create dataframe with lakes as rows
surveyed_sites_lake <- as.data.frame(t(surveyed_sites_species))
# Reorder lake names alphabetically
surveyed_sites_lake <- surveyed_sites_lake[order(row.names(surveyed_sites_lake)),]

# Kept this here because HLM is holomictic but is mixed, so it differs between datasets #
surveyed_sites_com_abund <- surveyed_sites_lake
surveyed_sites_com_abund[25,] <- colSums(surveyed_sites_com_abund[c(7,9,11,16,18:19),])
rownames(surveyed_sites_com_abund)[rownames(surveyed_sites_com_abund) == "23"] <- "Ocean sites"
surveyed_sites_com_abund[24,] <- colSums(surveyed_sites_com_abund[c(3,5:6,8,10,13:15,22),])
rownames(surveyed_sites_com_abund)[rownames(surveyed_sites_com_abund) == "24"] <- "Holomictic lakes"
surveyed_sites_com_abund[23,] <- colSums(surveyed_sites_com_abund[c(1:2,4,12,17,20,21),])
rownames(surveyed_sites_com_abund)[rownames(surveyed_sites_com_abund) == "25"] <- "Meromictic lakes"
surveyed_sites_com <- surveyed_sites_com_abund
surveyed_sites_com[surveyed_sites_com > "1"] <- 1

surveyed_sites_strat_abund <- surveyed_sites_lake
surveyed_sites_strat_abund[25,] <- colSums(surveyed_sites_strat_abund[c(7,9,11,16,18:19),])
rownames(surveyed_sites_strat_abund)[rownames(surveyed_sites_strat_abund) == "23"] <- "Ocean sites"
surveyed_sites_strat_abund[24,] <- colSums(surveyed_sites_strat_abund[c(3,6,8,10,13:15,22),])
rownames(surveyed_sites_strat_abund)[rownames(surveyed_sites_strat_abund) == "24"] <- "Mixed lakes"
surveyed_sites_strat_abund[23,] <- colSums(surveyed_sites_strat_abund[c(1:2,4,5,12,17,20,21),])
rownames(surveyed_sites_strat_abund)[rownames(surveyed_sites_strat_abund) == "25"] <- "Stratified lakes"
surveyed_sites_strat <- surveyed_sites_strat_abund
surveyed_sites_strat[surveyed_sites_strat > "1"] <- 1

surveyed_fish <- surveyed_sites_fish[,1]
surveyed_fishe <- surveyed_sites_fish[,-c(8,9,15)]
surveyed_fishe <- surveyed_fishe[which(rowSums(surveyed_fishe[, -1]) > 0),]
surveyed_fishe <- surveyed_fishe[,1]
surveyed_fishb <- surveyed_sites_fish[,-c(10)]
surveyed_fishb <- surveyed_fishb[which(rowSums(surveyed_fishb[, -1]) > 0),]
surveyed_fishb <- surveyed_fishb[,1]


# Ocean sites surveyed
ocean_fish[,c(2:7)] <- sapply(ocean_fish[,c(2:7)], as.numeric)
row.names(ocean_fish) <- ocean_fish$X
ocean_species <- ocean_fish[,-1]
ocean_species <- ocean_species[order(row.names(ocean_species)),]
# Create dataframe with lakes as rows
ocean_site <- as.data.frame(t(ocean_species))
# Reorder lake names alphabetically
ocean_site <- ocean_site[order(row.names(ocean_site)),]


# Holomictic lake surveyed sites
holomictic_fish[,c(2:10)] <- sapply(holomictic_fish[,c(2:10)], as.numeric)
row.names(holomictic_fish) <- holomictic_fish$X
holomictic_species <- holomictic_fish[,-1]
holomictic_species <- holomictic_species[order(row.names(holomictic_species)),]
# Create dataframe with lakes as rows
holomictic_lake <- as.data.frame(t(holomictic_species))
# Reorder lake names alphabetically
holomictic_lake <- holomictic_lake[order(row.names(holomictic_lake)),]


# Meromictic lake surveyed sites
meromictic_fish[,c(2:7)] <- sapply(meromictic_fish[,c(2:7)], as.numeric)
row.names(meromictic_fish) <- meromictic_fish$X
meromictic_species <- meromictic_fish[,-1]
meromictic_species <- meromictic_species[order(row.names(meromictic_species)),]
# Create dataframe with lakes as rows
meromictic_lake <- as.data.frame(t(meromictic_species))
# Reorder lake names alphabetically
meromictic_lake <- meromictic_lake[order(row.names(meromictic_lake)),]


# Mixed lake surveyed sites
mixed_fish[,c(2:9)] <- sapply(mixed_fish[,c(2:9)], as.numeric)
row.names(mixed_fish) <- mixed_fish$X
mixed_species <- mixed_fish[,-1]
mixed_species <- mixed_species[order(row.names(mixed_species)),]
# Create dataframe with lakes as rows
mixed_lake <- as.data.frame(t(mixed_species))
# Reorder lake names alphabetically
mixed_lake <- mixed_lake[order(row.names(mixed_lake)),]


# Stratified lake surveyed sites
stratified_fish[,c(2:9)] <- sapply(stratified_fish[,c(2:9)], as.numeric)
row.names(stratified_fish) <- stratified_fish$X
stratified_species <- stratified_fish[,-1]
stratified_species <- stratified_species[order(row.names(stratified_species)),]
# Create dataframe with lakes as rows
stratified_lake <- as.data.frame(t(stratified_species))
# Reorder lake names alphabetically
stratified_lake <- stratified_lake[order(row.names(stratified_lake)),]


# Give sites pres abs matrices row names and make files with only species presence
row.names(BCM_fish) <- BCM_fish$X
BCM_species <- BCM_fish["BCM"]
row.names(CLM_fish) <- CLM_fish$X
CLM_species <- CLM_fish["CLM"]
row.names(FLK_fish) <- FLK_fish$X
FLK_species <- FLK_fish["FLK"]
row.names(GLK_fish) <- GLK_fish$X
GLK_species <- GLK_fish["GLK"]
row.names(HLM_fish) <- HLM_fish$X
HLM_species <- HLM_fish["HLM"]
row.names(HLO_fish) <- HLO_fish$X
HLO_species <- HLO_fish["HLO"]
row.names(IBK_fish) <- IBK_fish$X
IBK_species <- IBK_fish["IBK"]
row.names(LLN_fish) <- LLN_fish$X
LLN_species <- LLN_fish["LLN"]
row.names(LCN_fish) <- LCN_fish$X
LCN_species <- LCN_fish["LCN"]
row.names(MLN_fish) <- MLN_fish$X
MLN_species <- MLN_fish["MLN"]
row.names(NCN_fish) <- NCN_fish$X
NCN_species <- NCN_fish["NCN"]
row.names(NLK_fish) <- NLK_fish$X
NLK_species <- NLK_fish["NLK"]
row.names(NLN_fish) <- NLN_fish$X
NLN_species <- NLN_fish["NLN"]
row.names(NLU_fish) <- NLU_fish$X
NLU_species <- NLU_fish["NLU"]
row.names(OLO_fish) <- OLO_fish$X
OLO_species <- OLO_fish["OLO"]
row.names(OCO_fish) <- OCO_fish$X
OCO_species <- OCO_fish["OCO"]
row.names(OTM_fish) <- OTM_fish$X
OTM_species <- OTM_fish["OTM"]
row.names(OCM_fish) <- OCM_fish$X
OCM_species <- OCM_fish["OCM"]
row.names(RCA_fish) <- RCA_fish$X
RCA_species <- RCA_fish["RCA"]
row.names(SLN_fish) <- SLN_fish$X
SLN_species <- SLN_fish["SLN"]
row.names(TLN_fish) <- TLN_fish$X
TLN_species <- TLN_fish["TLN"]
row.names(ULN_fish) <- ULN_fish$X
ULN_species <- ULN_fish["ULN"]
```

## Modify phylogeny

``` r
# Make a tree of only surveyed site species and trees that include only sites with environmental or biogeographic data
stree <- drop.tip(tree, setdiff(tree$tip.label, surveyed_fish), trim.internal = TRUE, rooted = is.rooted(tree))

# Env tree
etree <- drop.tip(tree, setdiff(tree$tip.label, surveyed_fishe), trim.internal = TRUE, rooted = is.rooted(tree))

# Biogeo tree
btree <- drop.tip(tree, setdiff(tree$tip.label, surveyed_fishb), trim.internal = TRUE, rooted = is.rooted(tree))
```

## Modify trait data frames

``` r
# Salt, brack, fresh amalgamation
# Edit column trait information to simplify for analyses
species_functional_traits$WaterPref <- "all"
species_functional_traits$WaterPref[species_functional_traits$Fresh > 0 & species_functional_traits$Brack < 1 & species_functional_traits$Saltwater < 1] <- "fresh"
species_functional_traits$WaterPref[species_functional_traits$Fresh > 0 & species_functional_traits$Brack > 0 & species_functional_traits$Saltwater < 1] <- "fresh-brack"
species_functional_traits$WaterPref[species_functional_traits$Fresh < 1 & species_functional_traits$Brack > 0 & species_functional_traits$Saltwater < 1] <- "brack"
species_functional_traits$WaterPref[species_functional_traits$Fresh < 1 & species_functional_traits$Brack > 0 & species_functional_traits$Saltwater > 0] <- "brack-salt"
species_functional_traits$WaterPref[species_functional_traits$Fresh < 1 & species_functional_traits$Brack < 1 & species_functional_traits$Saltwater > 0] <- "salt"
species_functional_traits$WaterPref[species_functional_traits$Fresh > 0 & species_functional_traits$Brack < 1 & species_functional_traits$Saltwater > 0] <- "fresh-salt"

# Categorical variables
species_functional_traits <- species_functional_traits %>%
    mutate(BodyShapeI = dplyr::recode(BodyShapeI, 'elongated' = '4e', "short and / or deep" = "2s", 'fusiform / normal' =  '3f', 'eel-like' = '5l', 'other' = '1o'))
species_functional_traits <- species_functional_traits %>%
    mutate(DemersPelag = dplyr::recode(DemersPelag, 'bathydemersal' = '7bd', 'benthopelagic' = '6bp', 'demersal' = '5d', 'pelagic' =  '3p', 'pelagic-neritic' = '2pn', 'pelagic-oceanic' = '4po', 'reef-associated' = '1r'))
species_functional_traits <- species_functional_traits %>%
    mutate(FeedingPath = dplyr::recode(FeedingPath, 'benthic' = 'b', 'pelagic' = 'p'))
species_functional_traits <- species_functional_traits %>%
    mutate(RepGuild1 = dplyr::recode(RepGuild1, 'guarders' = '2g', 'bearers' = '1b', 'nonguarders' =  '3n'))
species_functional_traits <- species_functional_traits %>%
    mutate(RepGuild2 = dplyr::recode(RepGuild2, 'internal live bearers' = '1ib', 'Internal live bearers' = '1ib', 'open water/substratum egg scatterers' = '6s', 'nesters' =  '3n', 'Nesters' =  '3n', 'brood hiders' = '5h', 'clutch tenders' = '4t', 'external brooders' = '2eb'))
species_functional_traits <- species_functional_traits %>%
    mutate(ParentalCare = dplyr::recode(ParentalCare, 'none' = '4n', 'paternal' = '3p', 'maternal' = '2m', 'biparental' = '1b'))
species_functional_traits <- species_functional_traits %>%
    mutate(WaterPref = dplyr::recode(WaterPref, 'all' = '3a', 'salt' = '1s', 'brack-salt' = '2bs', 'brack' = '4b', 'fresh-brack' = '5fb', 'fresh' = '6f'))

# Operculum 
species_functional_traits$OperculumPresent[species_functional_traits$OperculumPresent > 0] <- "yes"
species_functional_traits$OperculumPresent[species_functional_traits$OperculumPresent < 1] <- "no"

# DorsalSpinesMean
species_functional_traits$DorsalSpinesMean <- rowMeans(species_functional_traits[c("DorsalSpinesMax", "DorsalSpinesMin")])

# Get rid of Salt, brack, fresh column and Dorsal Spines columns
drop <- c("Saltwater", "Brack", "Fresh", "DorsalSpinesMax", "DorsalSpinesMin")
species_functional_traits <- species_functional_traits[, -which(names(species_functional_traits) %in% drop)]

# Edit fields
factor <- c("BodyShapeI", "DemersPelag", "OperculumPresent", "FeedingPath", "RepGuild1", "RepGuild2", "ParentalCare", "WaterPref")
species_functional_traits[,factor] <- lapply(species_functional_traits[,factor], function(x) as.factor(tolower(trimws(x))))
# species_functional_traits$BodyShapeI <- as.factor(species_functional_traits$BodyShapeI)
# species_functional_traits$DemersPelag <- as.factor(species_functional_traits$DemersPelag)
# species_functional_traits$OperculumPresent <- as.factor(species_functional_traits$OperculumPresent)
# species_functional_traits$FeedingPath <- as.factor(species_functional_traits$FeedingPath)
# species_functional_traits$RepGuild1 <- as.factor(species_functional_traits$RepGuild1)
# species_functional_traits$RepGuild2 <- as.factor(species_functional_traits$RepGuild2)
# species_functional_traits$ParentalCare <- as.factor(species_functional_traits$ParentalCare)
# species_functional_traits$WaterPref <- as.factor(species_functional_traits$WaterPref)

str(species_functional_traits)
```

```
## 'data.frame':	1724 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_pyroferus" "Acanthurus_sp" "Acanthurus_bariene" "Acanthurus_blochii" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Acanthuriformes" "Acanthuriformes" "Acanthuriformes" ...
##  $ Family          : chr  "Acanthuridae" "Acanthuridae" "Acanthuridae" "Acanthuridae" ...
##  $ Genus           : chr  "Acanthurus" "Acanthurus" "Acanthurus" "Acanthurus" ...
##  $ Species         : chr  "Acanthurus_pyroferus" "Acanthurus_sp" "Acanthurus_bariene" "Acanthurus_blochii" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 2 2 2 2 2 2 2 2 2 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 2 2 2 2 1 2 2 2 ...
##  $ MaxLengthTL     : num  29 37.7 50 54.9 54 ...
##  $ Troph           : num  2 2.28 2 2 2 ...
##  $ DepthMin        : num  4 2.67 6 2 4 ...
##  $ DepthMax        : num  60 51.4 50 15 131 ...
##  $ TempPrefMin     : num  25 24.7 24.5 24.9 20.3 ...
##  $ TempPrefMax     : num  28.8 28.8 28.4 28.8 28.9 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 1 1 1 1 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 3 3 3 3 3 3 3 3 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 6 6 6 6 6 6 6 6 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 4 4 4 4 4 4 4 4 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ DorsalSpinesMean: num  8 8.82 9 9 9 ...
```

``` r
# Make Species column the row names and get rid of columns that are not traits
row.names(species_functional_traits) <- species_functional_traits$X
# Reorder species names alphabetically
species_functional_traits <- species_functional_traits[order(row.names(species_functional_traits)),]
# Trim off other taxonomic levels other than species name
traits <- species_functional_traits[,-c(2:6)]
rs <- pres_abs_by_species[,c(1,21)]
rt <- full_join(traits, rs, by = "X")
rt$LocationSums <- rt$REF
rt <- rt[,-17]
rt$Strat <- "R"
rt <- cbind(rt, env[20, c(2,4,6,8,10,19,22:32)])
```

```
## Warning in data.frame(..., check.names = FALSE): row names were found from a short
## variable and have been discarded
```

``` r
# Warning message:
# In data.frame(..., check.names = FALSE) :
#   row names were found from a short variable and have been discarded
rt <- rt %>% relocate(Location_type, .before = temperature_median)
traits <- traits[,-1]

# Surveyed
surveyed_sites_traits <- species_functional_traits[species_functional_traits$X %in% surveyed_sites_fish$X,]
str(surveyed_sites_traits)
```

```
## 'data.frame':	249 obs. of  21 variables:
##  $ X               : chr  "Abudefduf_lorenzi" "Abudefduf_septemfasciatus" "Abudefduf_sexfasciatus" "Acanthurus_lineatus" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Ovalentaria/misc" "Ovalentaria/misc" "Ovalentaria/misc" "Acanthuriformes" ...
##  $ Family          : chr  "Pomacentridae" "Pomacentridae" "Pomacentridae" "Acanthuridae" ...
##  $ Genus           : chr  "Abudefduf" "Abudefduf" "Abudefduf" "Acanthurus" ...
##  $ Species         : chr  "Abudefduf_lorenzi" "Abudefduf_septemfasciatus" "Abudefduf_sexfasciatus" "Acanthurus_lineatus" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 2 2 2 2 2 2 3 4 4 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 2 2 2 2 2 2 1 1 ...
##  $ MaxLengthTL     : num  18 23 19 38 50.3 ...
##  $ Troph           : num  2.68 3.01 2.7 2 2.99 ...
##  $ DepthMin        : num  1 0 1 0 1 ...
##  $ DepthMax        : num  6 3 20 15 30 ...
##  $ TempPrefMin     : num  26 25.6 24.7 24.7 25.3 ...
##  $ TempPrefMax     : num  29.4 29.3 29.3 29.3 29.3 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 1 1 1 1 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 2 2 3 3 3 3 3 2 2 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 3 3 6 6 6 6 6 3 3 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 3 3 4 4 4 4 4 3 3 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 2 1 1 1 3 2 ...
##  $ DorsalSpinesMean: num  13 13 13 9 9 ...
```

``` r
# Make Species column the row names and get rid of columns that are not traits
row.names(surveyed_sites_traits) <- surveyed_sites_traits$X
# Reorder species names alphabetically
surveyed_sites_traits <- surveyed_sites_traits[order(row.names(surveyed_sites_traits)),]
# Trim off other taxonomic levels other than species name
straits <- surveyed_sites_traits[,-c(2:6)]
surveyed_sites_fish$LocationSums <- rowSums(surveyed_sites_fish[,-1])
# # ss <- surveyed_sites_fish[,c(1,24)]
# # st <- full_join(straits, ss, by = "X")
# # st$Strat <- "S"
straits <- straits[,-1]
```

## Modify Location_type trait data frames

``` r
# Ocean
ocean_traits <- species_functional_traits[species_functional_traits$X %in% ocean_fish$X,]
str(ocean_traits)
```

```
## 'data.frame':	166 obs. of  21 variables:
##  $ X               : chr  "Abudefduf_lorenzi" "Abudefduf_septemfasciatus" "Abudefduf_sexfasciatus" "Acanthurus_lineatus" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Ovalentaria/misc" "Ovalentaria/misc" "Ovalentaria/misc" "Acanthuriformes" ...
##  $ Family          : chr  "Pomacentridae" "Pomacentridae" "Pomacentridae" "Acanthuridae" ...
##  $ Genus           : chr  "Abudefduf" "Abudefduf" "Abudefduf" "Acanthurus" ...
##  $ Species         : chr  "Abudefduf_lorenzi" "Abudefduf_septemfasciatus" "Abudefduf_sexfasciatus" "Acanthurus_lineatus" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 2 2 2 2 2 2 3 2 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 2 2 2 2 2 2 2 1 ...
##  $ MaxLengthTL     : num  18 23 19 38 50.3 ...
##  $ Troph           : num  2.68 3.01 2.7 2 2.99 ...
##  $ DepthMin        : num  1 0 1 0 1 ...
##  $ DepthMax        : num  6 3 20 15 30 ...
##  $ TempPrefMin     : num  26 25.6 24.7 24.7 25.3 ...
##  $ TempPrefMax     : num  29.4 29.3 29.3 29.3 29.3 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 1 1 1 1 2 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 2 2 3 3 3 3 3 2 2 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 3 3 6 6 6 6 6 3 3 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 3 3 4 4 4 4 4 3 3 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 2 1 1 1 1 1 ...
##  $ DorsalSpinesMean: num  13 13 13 9 9 ...
```

``` r
# Make Species column the row names and get rid of columns that are not traits
row.names(ocean_traits) <- ocean_traits$X
# Reorder species names alphabetically
ocean_traits <- ocean_traits[order(row.names(ocean_traits)),]
# Trim off other taxonomic levels other than species name
otraits <- ocean_traits[,-c(2:6)]
ocean_fish$LocationSums <- rowSums(ocean_fish[,-1])
os <- ocean_fish[,c(1,8)]
ot <- full_join(otraits, os, by = "X")
ot$Strat <- "O"
ot <- cbind(ot, Location_type_env[1,])
```

```
## Warning in data.frame(..., check.names = FALSE): row names were found from a short
## variable and have been discarded
```

``` r
# Warning message:
# In data.frame(..., check.names = FALSE) :
#   row names were found from a short variable and have been discarded
otraits <- otraits[,-1]

# Mixed
mixed_traits <- species_functional_traits[species_functional_traits$X %in% mixed_fish$X,]
str(mixed_traits)
```

```
## 'data.frame':	185 obs. of  21 variables:
##  $ X               : chr  "Abudefduf_lorenzi" "Abudefduf_sexfasciatus" "Acanthurus_nigricauda" "Acanthurus_xanthopterus" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Ovalentaria/misc" "Ovalentaria/misc" "Acanthuriformes" "Acanthuriformes" ...
##  $ Family          : chr  "Pomacentridae" "Pomacentridae" "Acanthuridae" "Acanthuridae" ...
##  $ Genus           : chr  "Abudefduf" "Abudefduf" "Acanthurus" "Acanthurus" ...
##  $ Species         : chr  "Abudefduf_lorenzi" "Abudefduf_sexfasciatus" "Acanthurus_nigricauda" "Acanthurus_xanthopterus" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 2 2 3 4 2 3 4 4 4 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 6 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 2 2 1 2 1 1 1 1 ...
##  $ MaxLengthTL     : num  18 19 50.3 70 14 ...
##  $ Troph           : num  2.68 2.7 2.99 2.87 3.38 ...
##  $ DepthMin        : num  1 1 1 1 5 ...
##  $ DepthMax        : num  6 20 30 100 35 40 15 10 52 35 ...
##  $ TempPrefMin     : num  26 24.7 25.3 23.3 25.2 ...
##  $ TempPrefMax     : num  29.4 29.3 29.3 29 29.1 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 1 2 1 1 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 2 3 3 2 2 2 2 2 2 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 3 6 6 3 3 3 3 3 3 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 3 4 4 3 3 3 3 3 3 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 2 1 2 1 1 1 1 1 ...
##  $ DorsalSpinesMean: num  13 13 9 8.5 6.5 13 6.5 6.5 7 7 ...
```

``` r
# Make Species column the row names and get rid of columns that are not traits
row.names(mixed_traits) <- mixed_traits$X
# Reorder species names alphabetically
mixed_traits <- mixed_traits[order(row.names(mixed_traits)),]
# Trim off other taxonomic levels other than species name
mxtraits <- mixed_traits[,-c(2:6)]
mixed_fish$LocationSums <- rowSums(mixed_fish[,-1])
mxs <- mixed_fish[,c(1,10)]
mxt <- full_join(mxtraits, mxs, by = "X")
mxt$Strat <- "M"
mxt <- cbind(mxt, Location_type_env[2,])
```

```
## Warning in data.frame(..., check.names = FALSE): row names were found from a short
## variable and have been discarded
```

``` r
# Warning message:
# In data.frame(..., check.names = FALSE) :
#   row names were found from a short variable and have been discarded
mxtraits <- mxtraits[,-1]


# Stratified
stratified_traits <- species_functional_traits[species_functional_traits$X %in% stratified_fish$X,]
str(stratified_traits)
```

```
## 'data.frame':	30 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_xanthopterus" "Acentrogobius_janthinopterus" "Arothron_reticularis" "Atherinomorus_endrachtensis" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Gobiiformes" "Tetraodontiformes" "Atheriniformes" ...
##  $ Family          : chr  "Acanthuridae" "Gobiidae" "Tetraodontidae" "Atherinidae" ...
##  $ Genus           : chr  "Acanthurus" "Acentrogobius" "Arothron" "Atherinomorus" ...
##  $ Species         : chr  "Acanthurus_xanthopterus" "Acentrogobius_janthinopterus" "Arothron_reticularis" "Atherinomorus_endrachtensis" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 3 4 2 4 4 3 3 3 2 4 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 1 1 2 2 2 2 2 2 1 ...
##  $ MaxLengthTL     : num  70 12.5 54.9 11 13.5 ...
##  $ Troph           : num  2.87 3.41 3.44 3.4 3.27 ...
##  $ DepthMin        : num  1 0 1 0 0.167 ...
##  $ DepthMax        : num  100 20 25 50 31.8 ...
##  $ TempPrefMin     : num  23.3 25.2 25.1 26.7 25.9 ...
##  $ TempPrefMax     : num  29 29.3 29.3 29.1 28.9 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 2 1 2 2 2 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 2 3 3 3 3 3 3 3 2 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 3 6 6 6 6 6 6 6 3 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 3 4 4 4 4 4 4 4 3 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 3 2 2 1 1 2 3 1 1 ...
##  $ DorsalSpinesMean: num  8.5 6.5 0 6 6.83 ...
```

``` r
# Make Species column the row names and get rid of columns that are not traits
row.names(stratified_traits) <- stratified_traits$X
# Reorder species names alphabetically
stratified_traits <- stratified_traits[order(row.names(stratified_traits)),]
# Trim off other taxonomic levels other than species name
sttraits <- stratified_traits[,-c(2:6)]
stratified_fish$LocationSums <- rowSums(stratified_fish[,-1])
sts <- stratified_fish[,c(1,10)]
stt <- full_join(sttraits, sts, by = "X")
stt$Strat <- "S"
stt <- cbind(stt, Location_type_env[3,])
```

```
## Warning in data.frame(..., check.names = FALSE): row names were found from a short
## variable and have been discarded
```

``` r
# Warning message:
# In data.frame(..., check.names = FALSE) :
#   row names were found from a short variable and have been discarded
sttraits <- sttraits[,-1]
```

### Modify Location_type trait data

``` r
# All
Location_type_at <- rbind(rt, ot, mxt, stt)

Location_type_at$Strat <- factor(Location_type_at$Strat, levels = c("R", "O", "M", "S"))

Location_type_at_weighted <- Location_type_at %>%
  group_by(Location_type) %>%
  uncount(weights = LocationSums)

Location_type_st <- rbind(ot, mxt, stt)

Location_type_st$Strat <- factor(Location_type_st$Strat, levels = c("O", "M", "S"))

Location_type_st_weighted <- Location_type_st %>%
  group_by(Location_type) %>%
  uncount(weights = LocationSums)
```

### Location_type trait data tests

``` r
# # Apply the normality test function to each column in your dataframe
# shapiro_results <- lapply(Location_type_at[,c(2:17)], test_normality)
# 
# ## For numerical traits
# # Create an empty vector to store numerical trait p-values
# Location_type_at_nt_p_values <- c()
# 
# # Loop through each trait (numerical variable)
# for (trait in colnames(Location_type_at[,c(4:10,12:13)])) {
#   # Use the Kruskal-Wallis test
#   kruskal_result <- kruskal.test(Location_type_at[[trait]] ~ Location_type, data = Location_type_at)
#   # Extract and store the p-value
#   Location_type_at_nt_p_values <- c(Location_type_at_nt_p_values, kruskal_result$p.value)
# }
# 
# # Display the p-values for each trait
# Location_type_at_nt_p_values
# 
# ## For categorical traits
# # Create an empty vector to store Location_type_ategorical trait p-values
# Location_type_at_ct_p_values <- numeric()
# 
# # Loop through each categorical variable (trait) in your data
# for (trait in colnames(Location_type_at[,c(2:3,11,14:17)])) {
#   # Create a contingency table for the current trait and Location_type
#   contingency_table <- table(Location_type_at$Location_type, Location_type_at[[trait]])
#   # Perform Fisher's Exact Test
#   fisher_test_result <- fisher.test(contingency_table, simulate.p.value = T)
#   # Extract and store the p-value
#   p_value <- fisher_test_result$p.value
#   Location_type_at_ct_p_values <- c(Location_type_at_ct_p_values, p_value)
# }
# 
# # Display the p-values for each trait
# Location_type_at_ct_p_values
# 
# # Example usage
# Location_type_at_traits <- colnames(st[,c(2:17)])
# Location_type_at_env_vars <- colnames(st[,c(20:36)])
# Location_type_at_test_results <- perform_tests(Location_type_at, Location_type_at_traits, Location_type_at_env_vars)
# 
# # Usage example
# Location_type_at_p_values <- extract_p_values(Location_type_at_test_results)
```

## Modify site trait data frames

``` r
# BCM
BCM_traits <- species_functional_traits[species_functional_traits$X %in% BCM_fish$X,]
str(BCM_traits)
```

```
## 'data.frame':	5 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_xanthopterus" "Acentrogobius_janthinopterus" "Atherinomorus_sp" "Eleotris_fusca" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Gobiiformes" "Atheriniformes" "Gobiiformes" ...
##  $ Family          : chr  "Acanthuridae" "Gobiidae" "Atherinidae" "Eleotridae" ...
##  $ Genus           : chr  "Acanthurus" "Acentrogobius" "Atherinomorus" "Eleotris" ...
##  $ Species         : chr  "Acanthurus_xanthopterus" "Acentrogobius_janthinopterus" "Atherinomorus_sp" "Eleotris_fusca" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 3 4 4 4 2
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 5 1
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 1 2 1 1
##  $ MaxLengthTL     : num  70 12.5 13.5 26 11
##  $ Troph           : num  2.87 3.41 3.27 3.79 3.98
##  $ DepthMin        : num  1 0 0.167 0 0
##  $ DepthMax        : num  100 20 31.8 5 3
##  $ TempPrefMin     : num  23.3 25.2 25.9 25.4 25.6
##  $ TempPrefMax     : num  29 29.3 28.9 29.3 29.3
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 2
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 2 3 2 1
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 3 6 4 2
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 3 4 3 3
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 3 1 3 3
##  $ DorsalSpinesMean: num  8.5 6.5 6.83 7 7
```

``` r
# Reorder species alphabetically
BCM_traits <- BCM_traits[order(row.names(BCM_traits)),]
# Trim off other taxonomic levels other than species name
BCMt <- BCM_traits[,-c(1:6)]
# BCMt$Location <- "BCM"
BCMt <- cbind(BCMt, env[1,])

# CLM
CLM_traits <- species_functional_traits[species_functional_traits$X %in% CLM_fish$X,]
str(CLM_traits)
```

```
## 'data.frame':	3 obs. of  21 variables:
##  $ X               : chr  "Acentrogobius_janthinopterus" "Fibramia_lateralis" "Sphaeramia_orbicularis"
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei"
##  $ Order           : chr  "Gobiiformes" "Kurtiformes" "Kurtiformes"
##  $ Family          : chr  "Gobiidae" "Apogonidae" "Apogonidae"
##  $ Genus           : chr  "Acentrogobius" "Fibramia" "Sphaeramia"
##  $ Species         : chr  "Acentrogobius_janthinopterus" "Fibramia_lateralis" "Sphaeramia_orbicularis"
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 4 2 2
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 1 1 1
##  $ MaxLengthTL     : num  12.5 11 10
##  $ Troph           : num  3.41 3.98 3.64
##  $ DepthMin        : num  0 0 0
##  $ DepthMax        : num  20 3 5
##  $ TempPrefMin     : num  25.2 25.6 26.3
##  $ TempPrefMax     : num  29.3 29.3 29.3
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 2 2
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 1 1
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 2 2
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 3 3
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 3 3 1
##  $ DorsalSpinesMean: num  6.5 7 8
```

``` r
# Reorder species alphabetically
CLM_traits <- CLM_traits[order(row.names(CLM_traits)),]
# Trim off other taxonomic levels other than species name
CLMt <- CLM_traits[,-c(1:6)]
# CLMt$Location <- "CLM"
CLMt <- cbind(CLMt, env[2,])

# FLK
FLK_traits <- species_functional_traits[species_functional_traits$X %in% FLK_fish$X,]
str(FLK_traits)
```

```
## 'data.frame':	32 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_xanthopterus" "Amblygobius_buanensis" "Amblygobius_linki" "Asterropteryx_sp" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Gobiiformes" "Gobiiformes" "Gobiiformes" ...
##  $ Family          : chr  "Acanthuridae" "Gobiidae" "Gobiidae" "Gobiidae" ...
##  $ Genus           : chr  "Acanthurus" "Amblygobius" "Amblygobius" "Asterropteryx" ...
##  $ Species         : chr  "Acanthurus_xanthopterus" "Amblygobius_buanensis" "Amblygobius_linki" "Asterropteryx_sp" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 3 3 4 4 4 3 2 2 2 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 1 1 1 2 2 2 2 2 2 ...
##  $ MaxLengthTL     : num  70 9.15 6.5 3.81 13.46 ...
##  $ Troph           : num  2.87 2.72 2.81 3.26 3.27 ...
##  $ DepthMin        : num  1 0 1 12.333 0.167 ...
##  $ DepthMax        : num  100 15 10 39.2 31.8 ...
##  $ TempPrefMin     : num  23.3 28.7 25.4 25.8 25.9 ...
##  $ TempPrefMax     : num  29 29.5 29.3 29.1 28.9 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 1 2 1 1 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 2 2 2 3 3 3 3 3 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 3 3 3 6 6 6 6 6 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 3 3 3 4 4 4 4 4 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 1 3 1 1 1 1 ...
##  $ DorsalSpinesMean: num  8.5 6.5 6.5 6.94 6.83 ...
```

``` r
# Reorder species alphabetically
FLK_traits <- FLK_traits[order(row.names(FLK_traits)),]
# Trim off other taxonomic levels other than species name
FLKt <- FLK_traits[,-c(1:6)]
# FLKt$Location <- "FLK"
FLKt <- cbind(FLKt, env[3,])

# GLK
GLK_traits <- species_functional_traits[species_functional_traits$X %in% GLK_fish$X,]
str(GLK_traits)
```

```
## 'data.frame':	5 obs. of  21 variables:
##  $ X               : chr  "Acentrogobius_janthinopterus" "Eleotris_fusca" "Exyrias_puntang" "Megalops_cyprinoides" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Gobiiformes" "Gobiiformes" "Gobiiformes" "Elopiformes" ...
##  $ Family          : chr  "Gobiidae" "Eleotridae" "Gobiidae" "Megalopidae" ...
##  $ Genus           : chr  "Acentrogobius" "Eleotris" "Exyrias" "Megalops" ...
##  $ Species         : chr  "Acentrogobius_janthinopterus" "Eleotris_fusca" "Exyrias_puntang" "Megalops_cyprinoides" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 4 4 3 4 4
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 5 1 6 1
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 1 1 2 2 2
##  $ MaxLengthTL     : num  12.5 26 16.2 150 5.49
##  $ Troph           : num  3.41 3.79 3.5 3.48 3.21
##  $ DepthMin        : num  0 0 0 50 0
##  $ DepthMax        : num  20 5 5 200 15
##  $ TempPrefMin     : num  25.2 25.4 26.2 19.8 26.1
##  $ TempPrefMax     : num  29.3 29.3 29.3 27.6 29.3
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 1
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 2 2 3 2
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 4 3 6 3
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 3 3 4 3
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 3 3 2 3 2
##  $ DorsalSpinesMean: num  6.5 7 7 0 6.5
```

``` r
# Reorder species alphabetically
GLK_traits <- GLK_traits[order(row.names(GLK_traits)),]
# Trim off other taxonomic levels other than species name
GLKt <- GLK_traits[,-c(1:6)]
# GLKt$Location <- "GLK"
GLKt <- cbind(GLKt, env[4,])

# HLM
HLM_traits <- species_functional_traits[species_functional_traits$X %in% HLM_fish$X,]
str(HLM_traits)
```

```
## 'data.frame':	18 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_xanthopterus" "Acentrogobius_janthinopterus" "Atherinomorus_endrachtensis" "Caesio_cuning" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Gobiiformes" "Atheriniformes" "Eupercaria/misc" ...
##  $ Family          : chr  "Acanthuridae" "Gobiidae" "Atherinidae" "Caesionidae" ...
##  $ Genus           : chr  "Acanthurus" "Acentrogobius" "Atherinomorus" "Caesio" ...
##  $ Species         : chr  "Acanthurus_xanthopterus" "Acentrogobius_janthinopterus" "Atherinomorus_endrachtensis" "Caesio_cuning" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 3 4 4 3 3 3 4 4 4 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 5 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 1 2 2 2 2 1 1 1 2 ...
##  $ MaxLengthTL     : num  70 12.5 11 60 129.9 ...
##  $ Troph           : num  2.87 3.41 3.4 3.4 4.49 ...
##  $ DepthMin        : num  1 0 0 1 0 0 0 2 0 0 ...
##  $ DepthMax        : num  100 20 50 60 190 146 30 18 5 5 ...
##  $ TempPrefMin     : num  23.3 25.2 26.7 26.1 23.2 ...
##  $ TempPrefMax     : num  29 29.3 29.1 29.1 29 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 2 2 2 2 1 1 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 2 3 3 3 3 2 2 2 2 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 3 6 6 6 6 3 3 4 3 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 3 4 4 4 4 3 3 3 3 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 3 2 1 2 3 1 1 3 2 ...
##  $ DorsalSpinesMean: num  8.5 6.5 6 10 9 9 6.5 6.5 7 7 ...
```

``` r
# Reorder species alphabetically
HLM_traits <- HLM_traits[order(row.names(HLM_traits)),]
# Trim off other taxonomic levels other than species name
HLMt <- HLM_traits[,-c(1:6)]
# HLMt$Location <- "HLM"
HLMt <- cbind(HLMt, env[5,])

# HLO
HLO_traits <- species_functional_traits[species_functional_traits$X %in% HLO_fish$X,]
str(HLO_traits)
```

```
## 'data.frame':	53 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_xanthopterus" "Amblyglyphidodon_curacao" "Amblygobius_buanensis" "Caesio_caerulaurea" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Ovalentaria/misc" "Gobiiformes" "Eupercaria/misc" ...
##  $ Family          : chr  "Acanthuridae" "Pomacentridae" "Gobiidae" "Caesionidae" ...
##  $ Genus           : chr  "Acanthurus" "Amblyglyphidodon" "Amblygobius" "Caesio" ...
##  $ Species         : chr  "Acanthurus_xanthopterus" "Amblyglyphidodon_curacao" "Amblygobius_buanensis" "Caesio_caerulaurea" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 3 2 3 3 3 3 2 2 2 2 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 1 2 2 2 2 2 2 1 ...
##  $ MaxLengthTL     : num  70 11 9.15 50.39 60 ...
##  $ Troph           : num  2.87 2.63 2.72 3.4 3.4 ...
##  $ DepthMin        : num  1 1 0 1 1 0 1 0 0 3 ...
##  $ DepthMax        : num  100 40 15 50 60 190 60 30 170 30 ...
##  $ TempPrefMin     : num  23.3 24.7 28.7 24.7 26.1 23.2 25.5 25 23.7 24.7 ...
##  $ TempPrefMax     : num  29 29 29.5 29 29.1 29 29 29.3 29 29.3 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 2 1 2 2 2 1 1 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 2 2 3 3 3 3 3 3 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 3 3 6 6 6 6 6 6 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 3 3 4 4 4 4 4 4 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 1 2 1 1 1 1 ...
##  $ DorsalSpinesMean: num  8.5 13 6.5 10 10 9 12.5 13 12 13.5 ...
```

``` r
# Reorder species alphabetically
HLO_traits <- HLO_traits[order(row.names(HLO_traits)),]
# Trim off other taxonomic levels other than species name
HLOt <- HLO_traits[,-c(1:6)]
# HLOt$Location <- "HLO"
HLOt <- cbind(HLOt, env[6,])

# IBK
IBK_traits <- species_functional_traits[species_functional_traits$X %in% IBK_fish$X,]
str(IBK_traits)
```

```
## 'data.frame':	87 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_nigricauda" "Acanthurus_xanthopterus" "Amblyglyphidodon_curacao" "Amblygobius_decussatus" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Acanthuriformes" "Ovalentaria/misc" "Gobiiformes" ...
##  $ Family          : chr  "Acanthuridae" "Acanthuridae" "Pomacentridae" "Gobiidae" ...
##  $ Genus           : chr  "Acanthurus" "Acanthurus" "Amblyglyphidodon" "Amblygobius" ...
##  $ Species         : chr  "Acanthurus_nigricauda" "Acanthurus_xanthopterus" "Amblyglyphidodon_curacao" "Amblygobius_decussatus" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 3 2 4 2 2 3 3 2 2 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 2 1 2 2 2 2 1 1 ...
##  $ MaxLengthTL     : num  50.3 70 11 11.6 18.3 ...
##  $ Troph           : num  2.99 2.87 2.63 2.72 2.87 ...
##  $ DepthMin        : num  1 1 1 3 1 1 1 1 6 1 ...
##  $ DepthMax        : num  30 100 40 25 60 50 40 60 70 15 ...
##  $ TempPrefMin     : num  25.3 23.3 24.7 26.3 25.3 25 25.3 26.1 24.3 25.1 ...
##  $ TempPrefMax     : num  29.3 29 29 29 29 29 29.1 29.1 28.9 29.3 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 2 1 1 1 1 2 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 3 2 2 2 2 3 3 2 2 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 6 3 3 3 3 6 6 3 3 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 4 3 3 3 1 4 4 4 2 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 2 1 1 2 1 1 1 1 1 1 ...
##  $ DorsalSpinesMean: num  9 8.5 13 7 10 3 9 10 2 0 ...
```

``` r
# Reorder species alphabetically
IBK_traits <- IBK_traits[order(row.names(IBK_traits)),]
# Trim off other taxonomic levels other than species name
IBKt <- IBK_traits[,-c(1:6)]
# IBKt$Location <- "IBK"
IBKt <- cbind(IBKt, env[7,])

# LLN
LLN_traits <- species_functional_traits[species_functional_traits$X %in% LLN_fish$X,]
str(LLN_traits)
```

```
## 'data.frame':	106 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_xanthopterus" "Amblyglyphidodon_curacao" "Amblygobius_buanensis" "Amblygobius_phalaena" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Ovalentaria/misc" "Gobiiformes" "Gobiiformes" ...
##  $ Family          : chr  "Acanthuridae" "Pomacentridae" "Gobiidae" "Gobiidae" ...
##  $ Genus           : chr  "Acanthurus" "Amblyglyphidodon" "Amblygobius" "Amblygobius" ...
##  $ Species         : chr  "Acanthurus_xanthopterus" "Amblyglyphidodon_curacao" "Amblygobius_buanensis" "Amblygobius_phalaena" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 3 2 3 4 4 4 4 2 3 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 1 1 1 2 2 2 2 2 ...
##  $ MaxLengthTL     : num  70 11 9.15 15 6.5 ...
##  $ Troph           : num  2.87 2.63 2.72 3.63 2.36 ...
##  $ DepthMin        : num  1 1 0 1.93 1 ...
##  $ DepthMax        : num  100 40 15 52 20 ...
##  $ TempPrefMin     : num  23.3 24.7 28.7 26.1 24.7 ...
##  $ TempPrefMax     : num  29 29 29.5 29.3 29.3 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 2 1 1 1 2 1 1 1 2 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 2 2 2 2 3 3 2 3 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 3 3 3 3 6 6 3 6 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 3 3 3 3 4 4 1 4 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 1 2 1 1 1 1 ...
##  $ DorsalSpinesMean: num  8.5 13 6.5 7 7 ...
```

``` r
# Reorder species alphabetically
LLN_traits <- LLN_traits[order(row.names(LLN_traits)),]
# Trim off other taxonomic levels other than species name
LLNt <- LLN_traits[,-c(1:6)]
# LLNt$Location <- "LLN"
LLNt <- cbind(LLNt, env[8,])

# LCN
LCN_traits <- species_functional_traits[species_functional_traits$X %in% LCN_fish$X,]
str(LCN_traits)
```

```
## 'data.frame':	27 obs. of  21 variables:
##  $ X               : chr  "Abudefduf_septemfasciatus" "Acanthurus_triostegus" "Amblygobius_buanensis" "Atherinomorus_sp" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Ovalentaria/misc" "Acanthuriformes" "Gobiiformes" "Atheriniformes" ...
##  $ Family          : chr  "Pomacentridae" "Acanthuridae" "Gobiidae" "Atherinidae" ...
##  $ Genus           : chr  "Abudefduf" "Acanthurus" "Amblygobius" "Atherinomorus" ...
##  $ Species         : chr  "Abudefduf_septemfasciatus" "Acanthurus_triostegus" "Amblygobius_buanensis" "Atherinomorus_sp" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 2 3 4 2 2 3 2 3 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 1 2 1 2 2 2 2 2 ...
##  $ MaxLengthTL     : num  23 27 9.15 13.46 12.1 ...
##  $ Troph           : num  3.01 2.78 2.72 3.27 3.03 ...
##  $ DepthMin        : num  0 0 0 0.167 10 ...
##  $ DepthMax        : num  3 90 15 31.8 36 ...
##  $ TempPrefMin     : num  25.6 23.4 28.7 25.9 26.1 ...
##  $ TempPrefMax     : num  29.3 28.9 29.5 28.9 29 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 1 1 1 1 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 3 2 3 2 3 3 2 2 2 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 6 3 6 3 6 6 3 3 3 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 4 3 4 2 4 4 3 3 3 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 1 1 1 1 1 2 ...
##  $ DorsalSpinesMean: num  13 9 6.5 6.83 0 ...
```

``` r
# Reorder species alphabetically
LCN_traits <- LCN_traits[order(row.names(LCN_traits)),]
# Trim off other taxonomic levels other than species name
LCNt <- LCN_traits[,-c(1:6)]
# LCNt$Location <- "LCN"
LCNt <- cbind(LCNt, env[9,])

# MLN
MLN_traits <- species_functional_traits[species_functional_traits$X %in% MLN_fish$X,]
str(MLN_traits)
```

```
## 'data.frame':	64 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_nigricauda" "Acanthurus_xanthopterus" "Amblyeleotris_gymnocephala" "Amblyglyphidodon_curacao" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Acanthuriformes" "Gobiiformes" "Ovalentaria/misc" ...
##  $ Family          : chr  "Acanthuridae" "Acanthuridae" "Gobiidae" "Pomacentridae" ...
##  $ Genus           : chr  "Acanthurus" "Acanthurus" "Amblyeleotris" "Amblyglyphidodon" ...
##  $ Species         : chr  "Acanthurus_nigricauda" "Acanthurus_xanthopterus" "Amblyeleotris_gymnocephala" "Amblyglyphidodon_curacao" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 3 4 2 3 4 2 3 3 2 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 1 2 1 1 2 2 1 2 ...
##  $ MaxLengthTL     : num  50.28 70 14 11 9.15 ...
##  $ Troph           : num  2.99 2.87 3.38 2.63 2.72 ...
##  $ DepthMin        : num  1 1 5 1 0 ...
##  $ DepthMax        : num  30 100 35 40 15 52 50 50 50 30 ...
##  $ TempPrefMin     : num  25.3 23.3 25.2 24.7 28.7 ...
##  $ TempPrefMax     : num  29.3 29 29.1 29 29.5 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 2 1 1 1 2 2 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 3 2 2 2 2 2 3 3 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 6 3 3 3 3 3 6 6 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 4 3 3 3 3 1 4 4 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 2 1 2 1 1 1 1 1 2 1 ...
##  $ DorsalSpinesMean: num  9 8.5 6.5 13 6.5 7 3 10 9 13 ...
```

``` r
# Reorder species alphabetically
MLN_traits <- MLN_traits[order(row.names(MLN_traits)),]
# Trim off other taxonomic levels other than species name
MLNt <- MLN_traits[,-c(1:6)]
# MLNt$Location <- "MLN"
MLNt <- cbind(MLNt, env[10,])

# NCN
NCN_traits <- species_functional_traits[species_functional_traits$X %in% NCN_fish$X,]
str(NCN_traits)
```

```
## 'data.frame':	52 obs. of  21 variables:
##  $ X               : chr  "Amblygobius_decussatus" "Atherinomorus_endrachtensis" "Bolbometopon_muricatum" "Caesio_cuning" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Gobiiformes" "Atheriniformes" "Eupercaria/misc" "Eupercaria/misc" ...
##  $ Family          : chr  "Gobiidae" "Atherinidae" "Scaridae" "Caesionidae" ...
##  $ Genus           : chr  "Amblygobius" "Atherinomorus" "Bolbometopon" "Caesio" ...
##  $ Species         : chr  "Amblygobius_decussatus" "Atherinomorus_endrachtensis" "Bolbometopon_muricatum" "Caesio_cuning" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 4 4 3 3 2 2 2 2 2 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 1 2 2 2 2 2 1 2 2 2 ...
##  $ MaxLengthTL     : num  11.6 11 130 60 23 ...
##  $ Troph           : num  2.72 3.4 2.67 3.4 3.69 ...
##  $ DepthMin        : num  3 0 1 1 1 0 3 2 5 1 ...
##  $ DepthMax        : num  25 50 40 60 60 30 30 30 30 100 ...
##  $ TempPrefMin     : num  26.3 26.7 25.3 26.1 25.5 25 24.7 25 24.7 24.9 ...
##  $ TempPrefMax     : num  29 29.1 29.1 29.1 29 29.3 29.3 29.3 29.3 28.8 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 2 1 2 1 1 1 1 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 3 3 3 3 3 3 3 3 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 6 6 6 6 6 6 6 6 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 4 4 4 4 4 4 4 4 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 2 2 1 1 1 1 1 1 1 1 ...
##  $ DorsalSpinesMean: num  7 6 9 10 12.5 13 13.5 12 13 9 ...
```

``` r
# Reorder species alphabetically
NCN_traits <- NCN_traits[order(row.names(NCN_traits)),]
# Trim off other taxonomic levels other than species name
NCNt <- NCN_traits[,-c(1:6)]
# NCNt$Location <- "NCN"
NCNt <- cbind(NCNt, env[11,])

# NLK
NLK_traits <- species_functional_traits[species_functional_traits$X %in% NLK_fish$X,]
str(NLK_traits)
```

```
## 'data.frame':	6 obs. of  21 variables:
##  $ X               : chr  "Acentrogobius_janthinopterus" "Exyrias_puntang" "Lutjanus_argentimaculatus" "Lutjanus_fulvus" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Gobiiformes" "Gobiiformes" "Eupercaria/misc" "Eupercaria/misc" ...
##  $ Family          : chr  "Gobiidae" "Gobiidae" "Lutjanidae" "Lutjanidae" ...
##  $ Genus           : chr  "Acentrogobius" "Exyrias" "Lutjanus" "Lutjanus" ...
##  $ Species         : chr  "Acentrogobius_janthinopterus" "Exyrias_puntang" "Lutjanus_argentimaculatus" "Lutjanus_fulvus" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 4 3 3 3 3 2
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 1 2 2 2 2 1
##  $ MaxLengthTL     : num  12.5 16.2 150 40 43 ...
##  $ Troph           : num  3.41 3.5 3.58 3.61 2 ...
##  $ DepthMin        : num  0 0 1 1 0 0
##  $ DepthMax        : num  20 5 120 75 25 5
##  $ TempPrefMin     : num  25.2 26.2 24.3 24.6 24.9 26.3
##  $ TempPrefMax     : num  29.3 29.3 29.1 29 29 29.3
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 1 2
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 2 3 3 3 1
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 3 6 6 6 2
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 3 4 4 4 3
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 3 2 3 3 1 1
##  $ DorsalSpinesMean: num  6.5 7 10 10 13 8
```

``` r
# Reorder species alphabetically
NLK_traits <- NLK_traits[order(row.names(NLK_traits)),]
# Trim off other taxonomic levels other than species name
NLKt <- NLK_traits[,-c(1:6)]
# NLKt$Location <- "NLK"
NLKt <- cbind(NLKt, env[12,])

# NLN
NLN_traits <- species_functional_traits[species_functional_traits$X %in% NLN_fish$X,]
str(NLN_traits)
```

```
## 'data.frame':	65 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_xanthopterus" "Amblygobius_buanensis" "Ancistrogobius_dipus" "Apogonichthyoides_melas" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Gobiiformes" "Gobiiformes" "Kurtiformes" ...
##  $ Family          : chr  "Acanthuridae" "Gobiidae" "Gobiidae" "Apogonidae" ...
##  $ Genus           : chr  "Acanthurus" "Amblygobius" "Ancistrogobius" "Apogonichthyoides" ...
##  $ Species         : chr  "Acanthurus_xanthopterus" "Amblygobius_buanensis" "Ancistrogobius_dipus" "Apogonichthyoides_melas" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 3 3 4 2 4 4 2 3 3 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 6 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 1 1 1 2 2 2 2 2 2 ...
##  $ MaxLengthTL     : num  70 9.15 5.17 12.2 10.98 ...
##  $ Troph           : num  2.87 2.72 3.18 3.5 3.4 ...
##  $ DepthMin        : num  1 0 15 1 0 ...
##  $ DepthMax        : num  100 15 35 15 50 ...
##  $ TempPrefMin     : num  23.3 28.7 25 27.1 26.7 ...
##  $ TempPrefMax     : num  29 29.5 28.9 29.3 29.1 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 2 1 1 2 2 2 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 2 2 1 3 3 2 3 3 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 3 3 2 6 6 3 6 6 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 3 3 3 4 4 1 4 4 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 2 1 1 1 1 2 ...
##  $ DorsalSpinesMean: num  8.5 6.5 7 8 6 ...
```

``` r
# Reorder species alphabetically
NLN_traits <- NLN_traits[order(row.names(NLN_traits)),]
# Trim off other taxonomic levels other than species name
NLNt <- NLN_traits[,-c(1:6)]
# NLNt$Location <- "NLN"
NLNt <- cbind(NLNt, env[13,])

# NLU
NLU_traits <- species_functional_traits[species_functional_traits$X %in% NLU_fish$X,]
str(NLU_traits)
```

```
## 'data.frame':	21 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_xanthopterus" "Amblygobius_buanensis" "Atherinomorus_sp" "Chaetodon_ephippium" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Gobiiformes" "Atheriniformes" "Acanthuriformes" ...
##  $ Family          : chr  "Acanthuridae" "Gobiidae" "Atherinidae" "Chaetodontidae" ...
##  $ Genus           : chr  "Acanthurus" "Amblygobius" "Atherinomorus" "Chaetodon" ...
##  $ Species         : chr  "Acanthurus_xanthopterus" "Amblygobius_buanensis" "Atherinomorus_sp" "Chaetodon_ephippium" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 3 3 4 2 4 3 4 2 3 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 1 2 2 1 2 2 1 2 2 ...
##  $ MaxLengthTL     : num  70 9.15 13.46 30 21.96 ...
##  $ Troph           : num  2.87 2.72 3.27 3.01 4.12 ...
##  $ DepthMin        : num  1 0 0.167 0 2 ...
##  $ DepthMax        : num  100 15 31.8 30 10 ...
##  $ TempPrefMin     : num  23.3 28.7 25.9 25 28 ...
##  $ TempPrefMax     : num  29 29.5 28.9 29.3 29.3 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 1 1 1 2 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 2 3 3 1 3 2 1 3 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 3 6 6 2 6 3 2 6 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 3 4 4 3 4 3 3 4 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 1 1 1 3 3 2 ...
##  $ DorsalSpinesMean: num  8.5 6.5 6.83 13 7 ...
```

``` r
# Reorder species alphabetically
NLU_traits <- NLU_traits[order(row.names(NLU_traits)),]
# Trim off other taxonomic levels other than species name
NLUt <- NLU_traits[,-c(1:6)]
# NLUt$Location <- "NLU"
NLUt <- cbind(NLUt, env[14,])

# OLO
OLO_traits <- species_functional_traits[species_functional_traits$X %in% OLO_fish$X,]
str(OLO_traits)
```

```
## 'data.frame':	19 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_xanthopterus" "Amblygobius_buanensis" "Atherinomorus_endrachtensis" "Atherinomorus_sp" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Gobiiformes" "Atheriniformes" "Atheriniformes" ...
##  $ Family          : chr  "Acanthuridae" "Gobiidae" "Atherinidae" "Atherinidae" ...
##  $ Genus           : chr  "Acanthurus" "Amblygobius" "Atherinomorus" "Atherinomorus" ...
##  $ Species         : chr  "Acanthurus_xanthopterus" "Amblygobius_buanensis" "Atherinomorus_endrachtensis" "Atherinomorus_sp" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 3 3 4 4 3 2 3 3 3 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 1 2 2 2 2 2 2 2 2 ...
##  $ MaxLengthTL     : num  70 9.15 10.98 13.46 120 ...
##  $ Troph           : num  2.87 2.72 3.4 3.27 4.5 ...
##  $ DepthMin        : num  1 0 0 0.167 0 ...
##  $ DepthMax        : num  100 15 50 31.8 146 ...
##  $ TempPrefMin     : num  23.3 28.7 26.7 25.9 22.3 ...
##  $ TempPrefMax     : num  29 29.5 29.1 28.9 28.9 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 2 1 2 1 1 1 2 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 2 3 3 3 3 3 2 3 2 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 3 6 6 6 6 6 3 6 3 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 3 4 4 4 4 4 3 4 3 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 2 1 3 1 1 1 1 2 ...
##  $ DorsalSpinesMean: num  8.5 6.5 6 6.83 9 ...
```

``` r
# Reorder species alphabetically
OLO_traits <- OLO_traits[order(row.names(OLO_traits)),]
# Trim off other taxonomic levels other than species name
OLOt <- OLO_traits[,-c(1:6)]
# OLOt$Location <- "OLO"
OLOt <- cbind(OLOt, env[15,])

# OCO
OCO_traits <- species_functional_traits[species_functional_traits$X %in% OCO_fish$X,]
str(OCO_traits)
```

```
## 'data.frame':	25 obs. of  21 variables:
##  $ X               : chr  "Abudefduf_lorenzi" "Abudefduf_sexfasciatus" "Amblyglyphidodon_curacao" "Atherinomorus_endrachtensis" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Ovalentaria/misc" "Ovalentaria/misc" "Ovalentaria/misc" "Atheriniformes" ...
##  $ Family          : chr  "Pomacentridae" "Pomacentridae" "Pomacentridae" "Atherinidae" ...
##  $ Genus           : chr  "Abudefduf" "Abudefduf" "Amblyglyphidodon" "Atherinomorus" ...
##  $ Species         : chr  "Abudefduf_lorenzi" "Abudefduf_sexfasciatus" "Amblyglyphidodon_curacao" "Atherinomorus_endrachtensis" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 2 2 4 4 3 3 3 2 2 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 2 2 2 2 2 2 2 2 ...
##  $ MaxLengthTL     : num  18 19 11 11 13.5 ...
##  $ Troph           : num  2.68 2.7 2.63 3.4 3.27 ...
##  $ DepthMin        : num  1 1 1 0 0.167 ...
##  $ DepthMax        : num  6 20 40 50 31.8 ...
##  $ TempPrefMin     : num  26 24.7 24.7 26.7 25.9 ...
##  $ TempPrefMax     : num  29.4 29.3 29 29.1 28.9 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 2 2 1 2 2 2 2 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 2 2 3 3 3 3 3 1 2 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 3 3 6 6 6 6 6 2 3 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 3 3 4 4 4 4 4 3 3 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 2 1 1 2 1 1 2 ...
##  $ DorsalSpinesMean: num  13 13 13 6 6.83 ...
```

``` r
# Reorder species alphabetically
OCO_traits <- OCO_traits[order(row.names(OCO_traits)),]
# Trim off other taxonomic levels other than species name
OCOt <- OCO_traits[,-c(1:6)]
# OCOt$Location <- "OCO"
OCOt <- cbind(OCOt, env[16,])

# OTM
OTM_traits <- species_functional_traits[species_functional_traits$X %in% OTM_fish$X,]
str(OTM_traits)
```

```
## 'data.frame':	3 obs. of  21 variables:
##  $ X               : chr  "Acentrogobius_janthinopterus" "Atherinomorus_endrachtensis" "Sphaeramia_orbicularis"
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei"
##  $ Order           : chr  "Gobiiformes" "Atheriniformes" "Kurtiformes"
##  $ Family          : chr  "Gobiidae" "Atherinidae" "Apogonidae"
##  $ Genus           : chr  "Acentrogobius" "Atherinomorus" "Sphaeramia"
##  $ Species         : chr  "Acentrogobius_janthinopterus" "Atherinomorus_endrachtensis" "Sphaeramia_orbicularis"
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 4 4 2
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 1 2 1
##  $ MaxLengthTL     : num  12.5 11 10
##  $ Troph           : num  3.41 3.4 3.64
##  $ DepthMin        : num  0 0 0
##  $ DepthMax        : num  20 50 5
##  $ TempPrefMin     : num  25.2 26.7 26.3
##  $ TempPrefMax     : num  29.3 29.1 29.3
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 2 2
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 3 1
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 6 2
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 4 3
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 3 2 1
##  $ DorsalSpinesMean: num  6.5 6 8
```

``` r
# Reorder species alphabetically
OTM_traits <- OTM_traits[order(row.names(OTM_traits)),]
# Trim off other taxonomic levels other than species name
OTMt <- OTM_traits[,-c(1:6)]
# OTMt$Location <- "OTM"
OTMt <- cbind(OTMt, env[17,])

# OCM
OCM_traits <- species_functional_traits[species_functional_traits$X %in% OCM_fish$X,]
str(OCM_traits)
```

```
## 'data.frame':	72 obs. of  21 variables:
##  $ X               : chr  "Abudefduf_lorenzi" "Abudefduf_sexfasciatus" "Acanthurus_lineatus" "Amblyglyphidodon_curacao" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Ovalentaria/misc" "Ovalentaria/misc" "Acanthuriformes" "Ovalentaria/misc" ...
##  $ Family          : chr  "Pomacentridae" "Pomacentridae" "Acanthuridae" "Pomacentridae" ...
##  $ Genus           : chr  "Abudefduf" "Abudefduf" "Acanthurus" "Amblyglyphidodon" ...
##  $ Species         : chr  "Abudefduf_lorenzi" "Abudefduf_sexfasciatus" "Acanthurus_lineatus" "Amblyglyphidodon_curacao" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 2 2 2 3 4 4 3 2 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 2 2 1 2 2 2 2 2 ...
##  $ MaxLengthTL     : num  18 19 38 11 9.15 ...
##  $ Troph           : num  2.68 2.7 2 2.63 2.72 ...
##  $ DepthMin        : num  1 1 0 1 0 ...
##  $ DepthMax        : num  6 20 15 40 15 ...
##  $ TempPrefMin     : num  26 24.7 24.7 24.7 28.7 ...
##  $ TempPrefMax     : num  29.4 29.3 29.3 29 29.5 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 2 1 2 1 2 1 2 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 2 3 2 2 3 3 3 2 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 3 6 3 3 6 6 6 3 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 3 4 3 3 4 4 4 1 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 1 2 1 2 1 1 ...
##  $ DorsalSpinesMean: num  13 13 9 13 6.5 ...
```

``` r
# Reorder species alphabetically
OCM_traits <- OCM_traits[order(row.names(OCM_traits)),]
# Trim off other taxonomic levels other than species name
OCMt <- OCM_traits[,-c(1:6)]
# OCMt$Location <- "OCM"
OCMt <- cbind(OCMt, env[18,])

# RCA
RCA_traits <- species_functional_traits[species_functional_traits$X %in% RCA_fish$X,]
str(RCA_traits)
```

```
## 'data.frame':	54 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_sp" "Acanthurus_xanthopterus" "Amblyglyphidodon_curacao" "Amblygobius_buanensis" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Acanthuriformes" "Ovalentaria/misc" "Gobiiformes" ...
##  $ Family          : chr  "Acanthuridae" "Acanthuridae" "Pomacentridae" "Gobiidae" ...
##  $ Genus           : chr  "Acanthurus" "Acanthurus" "Amblyglyphidodon" "Amblygobius" ...
##  $ Species         : chr  "Acanthurus_sp" "Acanthurus_xanthopterus" "Amblyglyphidodon_curacao" "Amblygobius_buanensis" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 3 2 3 4 4 3 3 2 2 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 2 1 1 1 2 2 2 2 ...
##  $ MaxLengthTL     : num  37.7 70 11 9.15 10.37 ...
##  $ Troph           : num  2.28 2.87 2.63 2.72 2.74 ...
##  $ DepthMin        : num  2.67 1 1 0 2 ...
##  $ DepthMax        : num  51.4 100 40 15 15 ...
##  $ TempPrefMin     : num  24.7 23.3 24.7 28.7 28.2 ...
##  $ TempPrefMax     : num  28.8 29 29 29.5 29.3 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 2 1 1 1 2 2 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 3 2 2 2 2 3 3 3 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 6 3 3 3 3 6 6 6 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 4 3 3 3 3 4 4 4 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ DorsalSpinesMean: num  8.82 8.5 13 6.5 6.5 ...
```

``` r
# Reorder species alphabetically
RCA_traits <- RCA_traits[order(row.names(RCA_traits)),]
# Trim off other taxonomic levels other than species name
RCAt <- RCA_traits[,-c(1:6)]
# RCAt$Location <- "RCA"
RCAt <- cbind(RCAt, env[19,])

# SLN
SLN_traits <- species_functional_traits[species_functional_traits$X %in% SLN_fish$X,]
str(SLN_traits)
```

```
## 'data.frame':	3 obs. of  21 variables:
##  $ X               : chr  "Acentrogobius_janthinopterus" "Exyrias_puntang" "Mugilogobius_cavifrons"
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei"
##  $ Order           : chr  "Gobiiformes" "Gobiiformes" "Gobiiformes"
##  $ Family          : chr  "Gobiidae" "Gobiidae" "Gobiidae"
##  $ Genus           : chr  "Acentrogobius" "Exyrias" "Mugilogobius"
##  $ Species         : chr  "Acentrogobius_janthinopterus" "Exyrias_puntang" "Mugilogobius_cavifrons"
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 4 3 4
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 5
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 1 2 1
##  $ MaxLengthTL     : num  12.5 16.2 4
##  $ Troph           : num  3.41 3.5 3.29
##  $ DepthMin        : num  0 0 0.333
##  $ DepthMax        : num  20 5 10.6
##  $ TempPrefMin     : num  25.2 26.2 23.7
##  $ TempPrefMax     : num  29.3 29.3 29.2
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 2 2
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 3 4
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 3 3
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 3 2 3
##  $ DorsalSpinesMean: num  6.5 7 6.93
```

``` r
# Reorder species alphabetically
SLN_traits <- SLN_traits[order(row.names(SLN_traits)),]
# Trim off other taxonomic levels other than species name
SLNt <- SLN_traits[,-c(1:6)]
# SLNt$Location <- "SLN"
SLNt <- cbind(SLNt, env[21,])

# TLN
TLN_traits <- species_functional_traits[species_functional_traits$X %in% TLN_fish$X,]
str(TLN_traits)
```

```
## 'data.frame':	17 obs. of  21 variables:
##  $ X               : chr  "Acanthurus_xanthopterus" "Arothron_reticularis" "Atherinomorus_endrachtensis" "Atherinomorus_sp" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Acanthuriformes" "Tetraodontiformes" "Atheriniformes" "Atheriniformes" ...
##  $ Family          : chr  "Acanthuridae" "Tetraodontidae" "Atherinidae" "Atherinidae" ...
##  $ Genus           : chr  "Acanthurus" "Arothron" "Atherinomorus" "Atherinomorus" ...
##  $ Species         : chr  "Acanthurus_xanthopterus" "Arothron_reticularis" "Atherinomorus_endrachtensis" "Atherinomorus_sp" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 3 2 4 4 2 4 3 2 4 3 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 1 2 2 2 1 2 1 1 2 ...
##  $ MaxLengthTL     : num  70 54.9 11 13.5 20 ...
##  $ Troph           : num  2.87 3.44 3.4 3.27 3.7 ...
##  $ DepthMin        : num  1 1 0 0.167 0 ...
##  $ DepthMax        : num  100 25 50 31.8 170 ...
##  $ TempPrefMin     : num  23.3 25.1 26.7 25.9 23.7 ...
##  $ TempPrefMax     : num  29 29.3 29.1 28.9 29 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 2 1 1 1 1 2 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 3 3 3 3 3 2 2 1 1 3 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 6 6 6 6 6 3 3 2 2 6 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 4 4 4 4 4 3 3 3 3 4 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 2 2 1 1 1 2 3 1 3 ...
##  $ DorsalSpinesMean: num  8.5 0 6 6.83 12 ...
```

``` r
# Reorder species alphabetically
TLN_traits <- TLN_traits[order(row.names(TLN_traits)),]
# Trim off other taxonomic levels other than species name
TLNt <- TLN_traits[,-c(1:6)]
# TLNt$Location <- "TLN"
TLNt <- cbind(TLNt, env[22,])

# ULN
ULN_traits <- species_functional_traits[species_functional_traits$X %in% ULN_fish$X,]
str(ULN_traits)
```

```
## 'data.frame':	71 obs. of  21 variables:
##  $ X               : chr  "Abudefduf_lorenzi" "Abudefduf_sexfasciatus" "Acanthurus_nigricauda" "Acanthurus_xanthopterus" ...
##  $ Class           : chr  "Teleostei" "Teleostei" "Teleostei" "Teleostei" ...
##  $ Order           : chr  "Ovalentaria/misc" "Ovalentaria/misc" "Acanthuriformes" "Acanthuriformes" ...
##  $ Family          : chr  "Pomacentridae" "Pomacentridae" "Acanthuridae" "Acanthuridae" ...
##  $ Genus           : chr  "Abudefduf" "Abudefduf" "Acanthurus" "Acanthurus" ...
##  $ Species         : chr  "Abudefduf_lorenzi" "Abudefduf_sexfasciatus" "Acanthurus_nigricauda" "Acanthurus_xanthopterus" ...
##  $ BodyShapeI      : Factor w/ 5 levels "1o","2s","3f",..: 2 2 2 3 3 4 4 4 4 2 ...
##  $ DemersPelag     : Factor w/ 7 levels "1r","2pn","3p",..: 1 1 1 1 1 5 1 1 1 1 ...
##  $ OperculumPresent: Factor w/ 2 levels "no","yes": 2 2 2 2 1 1 1 2 2 2 ...
##  $ MaxLengthTL     : num  18 19 50.28 70 9.15 ...
##  $ Troph           : num  2.68 2.7 2.99 2.87 2.72 ...
##  $ DepthMin        : num  1 1 1 1 0 ...
##  $ DepthMax        : num  6 20 30 100 15 ...
##  $ TempPrefMin     : num  26 24.7 25.3 23.3 28.7 ...
##  $ TempPrefMax     : num  29.4 29.3 29.3 29 29.5 ...
##  $ FeedingPath     : Factor w/ 2 levels "b","p": 1 1 1 1 1 1 1 2 1 1 ...
##  $ RepGuild1       : Factor w/ 3 levels "1b","2g","3n": 2 2 3 3 2 2 2 3 3 2 ...
##  $ RepGuild2       : Factor w/ 6 levels "1ib","2eb","3n",..: 3 3 6 6 3 3 3 6 6 3 ...
##  $ ParentalCare    : Factor w/ 4 levels "1b","2m","3p",..: 3 3 4 4 3 3 3 4 4 1 ...
##  $ WaterPref       : Factor w/ 6 levels "1s","2bs","3a",..: 1 1 2 1 1 2 1 2 1 1 ...
##  $ DorsalSpinesMean: num  13 13 9 8.5 6.5 ...
```

``` r
# Reorder species alphabetically
ULN_traits <- ULN_traits[order(row.names(ULN_traits)),]
# Trim off other taxonomic levels other than species name
ULNt <- ULN_traits[,-c(1:6)]
# ULNt$Location <- "ULN"
ULNt <- cbind(ULNt, env[23,])

# All sites
Locations_at <- rbind(BCMt, CLMt, FLKt, GLKt, HLMt, HLOt, IBKt, LLNt, LCNt, MLNt, NCNt, NLKt, NLNt, NLUt, OLOt, OCOt, OTMt, OCMt, RCAt, SLNt, TLNt, ULNt)
# Rename column
Locations_at <- Locations_at %>%
  rename(Location = X)

# Calculate median and mode for columns by "Location"
Locations_at_summary <- Locations_at[1:16] %>%
  group_by(Location) %>%
  summarise(across(.cols = where(is.numeric), 
                .fns = list(median = median), 
                .names = "{.col}_{.fn}", 
                na.rm = TRUE),
         across(.cols = where(is.factor), 
                .fns = ~ names(which.max(table(.))), 
                .names = "{.col}_mode"))
```

```
## Warning: There was 1 warning in `summarise()`.
## ℹ In argument: `across(...)`.
## ℹ In group 1: `Location = "BCM"`.
## Caused by warning:
## ! The `...` argument of `across()` is deprecated as of dplyr 1.1.0.
## Supply arguments directly to `.fns` through an anonymous function instead.
## 
##   # Previously
##   across(a:b, mean, na.rm = TRUE)
## 
##   # Now
##   across(a:b, \(x) mean(x, na.rm = TRUE))
```

``` r
Locations_at_summary <- cbind(Locations_at_summary,env[-20,-1])

row.names(Locations_at_summary) <- Locations_at_summary$Location

# Marine lakes
Lakes_at <- rbind(BCMt, CLMt, FLKt, GLKt, HLMt, HLOt, LLNt, MLNt, NLKt, NLNt, NLUt, OLOt, OTMt, SLNt, TLNt, ULNt)
# Rename column
Lakes_at <- Lakes_at %>%
  rename(Location = X)

Ocean_at <- rbind(IBKt, LCNt, NCNt, OCOt, OCMt, RCAt)
Mixed_at <- rbind(FLKt, HLOt, LLNt, MLNt, NLNt, NLUt, OLOt, ULNt)
Stratified_at <- rbind(BCMt, CLMt, GLKt, HLMt, NLKt, OTMt, SLNt, TLNt)
```

## Load out modified files

``` r
# Environment
write.csv(env,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/environment/env.csv")

# Fish presence by lake
write.csv(pres_abs_lake_sums,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/pres_abs_lake_sums.csv")

# Fish presence by species
write.csv(pres_abs_species_sums,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/pres_abs_species_sums.csv")

# Fish presence in surveyed sites by Location_type abundance
write.csv(surveyed_sites_strat_abund,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/surveyed_sites_strat_abund.csv")

# Fish presence in ocean sites
write.csv(ocean_site,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/pres_abs/ocean_site_presence_by_site.csv")

# Phylogenetic tree
write.tree(tree,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/phylo/fish_tree.tre")

# Phylogeny of surveyed sites
write.tree(stree,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/phylo/surveyed_fish_tree.tre")

# Phylogeny of sites with environmental data
write.tree(etree,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/phylo/env_fish_tree.tre")

# Phylogeny of sites with biogeographic data
write.tree(btree,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/phylo/biogeo_fish_tree.tre")

# Traits
write.csv(species_functional_traits,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/traits/traits.csv")

# Traits in surveyed sites
write.csv(surveyed_sites_traits,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/traits/surveyed_sites_traits.csv")

# Traits in ocean sites
write.csv(ocean_traits,"/Users/bailey/Documents/research/Fish_Biodiversity/data/collection/fish/traits/ocean_traits.csv")

sessionInfo()
```

```
## R version 4.4.3 (2025-02-28)
## Platform: aarch64-apple-darwin20
## Running under: macOS 26.0.1
## 
## Matrix products: default
## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/Los_Angeles
## tzcode source: internal
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] tidyr_1.3.1    phytools_2.4-4 maps_3.4.2.1   ape_5.8-1      reshape2_1.4.4
## [6] stringr_1.5.1  dplyr_1.1.4    knitr_1.50    
## 
## loaded via a namespace (and not attached):
##  [1] generics_0.1.3          stringi_1.8.7           lattice_0.22-7         
##  [4] digest_0.6.37           magrittr_2.0.3          evaluate_1.0.3         
##  [7] grid_4.4.3              iterators_1.0.14        DEoptim_2.2-8          
## [10] fastmap_1.2.0           foreach_1.5.2           doParallel_1.0.17      
## [13] plyr_1.8.9              Matrix_1.7-3            optimParallel_1.0-2    
## [16] combinat_0.0-8          purrr_1.0.4             codetools_0.2-20       
## [19] numDeriv_2016.8-1.1     mnormt_2.1.1            cli_3.6.4              
## [22] rlang_1.1.6             expm_1.0-0              scatterplot3d_0.3-44   
## [25] withr_3.0.2             yaml_2.3.10             tools_4.4.3            
## [28] coda_0.19-4.1           fastmatch_1.1-6         vctrs_0.6.5            
## [31] R6_2.6.1                lifecycle_1.0.4         MASS_7.3-65            
## [34] pkgconfig_2.0.3         pillar_1.10.2           glue_1.8.0             
## [37] phangorn_2.12.1         Rcpp_1.0.14             clusterGeneration_1.3.8
## [40] xfun_0.52               tibble_3.2.1            tidyselect_1.2.1       
## [43] rstudioapi_0.17.1       htmltools_0.5.8.1       nlme_3.1-168           
## [46] igraph_2.1.4            rmarkdown_2.29          compiler_4.4.3         
## [49] quadprog_1.5-8
```
