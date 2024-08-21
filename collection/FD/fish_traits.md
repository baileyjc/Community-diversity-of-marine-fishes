Fish traits
================

#### R Markdown

## Load packages

- These are the packages that are needed to extract data from online
  repositories and manipulate the data

``` r
library(rfishbase) # data manipulation
library(tidyr) # data manipulation
library(dplyr) # data manipulation
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(stringr) # data manipulation
```

## Possible trait cateogries

- We won’t use this but it is a helpful tool to check out to find
  categories of traits

``` r
# fb_tables(server = "fishbase", version = "latest")
```

## Find certain traits

- Can be useful to see the variety of traits available from FishBase

``` r
# # Use these to explore data tables for possible traits
# fb_tbl("spawnagg", server = "fishbase")
# 
# spawnagg <- rfishbase:::endpoint("spawnagg")
# spawnagg <- spawnagg(existent_species)
# 
# larvaepresence <- rfishbase:::endpoint("larvaepresence")
# larvaepresence <- larvaepresence(existent_species)
```

## Load in the species list

``` r
existing_species <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/species_lists/existing_species.csv")

existing_species <- existing_species[,-1]
```

## Extract traits for species

``` r
## Species Traits
# Download traits for fish and drop any traits you do not want
trait <- species(existing_species)
```

    ## Joining with `by = join_by(SpecCode)`

``` r
#trait <- trait[, c(2:3, 15, 20:23, 38)]

## Morphological traits 
# Download traits for fish and drop any traits you do not want
morphology <- morphology(existing_species)
```

    ## Joining with `by = join_by(SpecCode)`

``` r
## Estimated traits
# Download traits for fish and drop any traits you do not want
# Trophic levels
estimate <- estimate(existing_species)
```

    ## Joining with `by = join_by(SpecCode)`

``` r
## Ecological traits 
# Download traits for fish and drop any traits you do not want
ecology <- ecology(existing_species)
```

    ## Joining with `by = join_by(SpecCode)`

``` r
## Reproduction traits
# Download traits for fish and drop any traits you do not want
reproduction <- reproduction(existing_species)
```

    ## Joining with `by = join_by(SpecCode)`

``` r
# Bind the four trait dataframes by columns, duplicated column names will be deleted with the second line of code.
traits <- full_join(trait, morphology, by = "Species")
traits <- full_join(traits, estimate, by = "Species")
traits <- full_join(traits, ecology, by = "Species")
traits <- full_join(traits, reproduction, by = "Species")

# traits <- cbind(trait, morphology, estimate, ecology, reproduction)
traits <- traits[,!duplicated(colnames(traits))]
```

## Completeness Function

- You will run a function written by Gio Rapacciuolo that will search
  out traits that are a certain percentage complete for all species.
  \#From trait-beta-div-functions.R

``` r
##### -- Check for the completeness of vectors -- #####
completeness <- function(vec){
  if (is.list(vec)) return(apply(vec, 2, function(x) (sum(complete.cases(x))/length(x)) * 100))
  else return((sum(complete.cases(vec))/length(vec)) * 100)
}
```

## Filter for 90% complete traits

- Now isolate for completeness of traits because you need to have traits
  that are recorded for most species to look at the diversity. This will
  get rid of traits that are lacking sufficient data.

``` r
# Isolate the most complete traits (completeness > 75% - 90%)
# complete_traits_75 <- traits[, which(completeness(traits) > 75)]
# 
# complete_traits_80 <- traits[, which(completeness(traits) > 80)]
# 
# complete_traits_85 <- traits[, which(completeness(traits) > 85)]

complete_traits_90 <- traits[, which(completeness(traits) > 90)]

# check <- complete_traits_75[setdiff(names(complete_traits_75), names(complete_traits_90))]

# Isolate traits that include variation
# complete_traits_75 <- complete_traits_75[, which(apply(complete_traits_75, 2, function(x) length(unique(na.omit(x)))) > 1)]
# 
# complete_traits_80 <- complete_traits_80[, which(apply(complete_traits_80, 2, function(x) length(unique(na.omit(x)))) > 1)]
# 
# complete_traits_85 <- complete_traits_85[, which(apply(complete_traits_85, 2, function(x) length(unique(na.omit(x)))) > 1)]

complete_traits_90 <- complete_traits_90[, which(apply(complete_traits_90, 2, function(x) length(unique(na.omit(x)))) > 1)]
```

## Select specific traits

- Edit how traits are viewed by R. Also fishbase records for
  environmental data presence is -1 however, programs have trouble with
  that so we need to convert all those -1 to 1 so that they can be
  better dealt with in R.

``` r
complete_traits <- complete_traits_90

# # Calculate the estimated Weight with the following columns.
# complete_traits$Weight <- complete_traits$a * complete_traits$MaxLengthTL^complete_traits$b
# 
# # Calculate the estimated caudal fin length with the following columns.
# complete_traits$CaudalFinLength <- complete_traits$MaxLengthTL - complete_traits$MaxLengthSL

# Get rid of unnecessary or irrelevant traits
keep <- c("Species", "Genus", "BodyShapeI.x", "Fresh", "Brack", "Saltwater", "DemersPelag", "OperculumPresent", "DorsalSpinesMax", "MaxLengthTL", "Troph", "DepthMin", "DepthMax", "TempPrefMin", "TempPrefMax", "FeedingPath")
complete_traits <- complete_traits[,keep]
# Traits we are looking to keep:
# Species   Genus   BodyShapeI.x Fresh Brack Saltwater  DemersPelag OperculumPresent DorsalSpinesMax MaxLengthTL Troph  DepthMin DepthMax TempPrefMin TempPrefMax   FeedingPath

# complete_traits <- complete_traits[, c(1,3,10,13:16,21,55,57,59:62,68,70:71,73:76,83,84,87:91,93,95,97,98)]
# Traits we are looking to keep:
# Species   Family  Genus   BodyShapeI.x    Fresh   Brack   Saltwater   DemersPelag Vulnerability   OperculumPresent    DorsalSpinesMin DorsalSpinesMax DorsalSoftRaysMin   DorsalSoftRaysMax   Araymin Araymax MaxLengthTL Troph   seTroph a   sd_log10a   b   sd_b    DepthMin    DepthMax    PredPreyRatioMin    PredPreyRatioMax    TempPrefMin TempPrefMean    TempPrefMax FeedingPath MaxLengthSL Weight  CaudalFinLength RepGuild1   RepGuild2   ParentalCare

# Make -1 in 1 by multiplication
complete_traits$OperculumPresent <- complete_traits$OperculumPresent*-1

# RepGuild1 RepGuild2   ParentalCare
keep <- c("Species", "RepGuild1", "RepGuild2", "ParentalCare")
repro <- reproduction[,keep]

# DorsalSpinesMin
keep <- c("Species", "DorsalSpinesMin")
dorsal <- morphology[,keep]

complete_traits <- full_join(complete_traits, repro, by = "Species")
complete_traits <- full_join(complete_traits, dorsal, by = "Species")

# Make integer fields numeric
numerical <- c("DorsalSpinesMax", "MaxLengthTL", "Troph", "DepthMin", "DepthMax", "TempPrefMin", "TempPrefMax", "DorsalSpinesMin")
complete_traits[,numerical] <- sapply(complete_traits[,numerical], as.numeric)

# Make character fields into factors to use in trait analyses
factor <- c("BodyShapeI.x", "Fresh", "Brack", "Saltwater", "DemersPelag", "OperculumPresent", "FeedingPath", "RepGuild1", "RepGuild2", "ParentalCare")
# Convert to lowercase and then to factors
complete_traits[,factor] <- lapply(complete_traits[,factor], function(x) as.factor(tolower(trimws(x))))

# Identify which rows have NAs and will be filled in with Genus level data
fromGenus <- apply(complete_traits, 1, function(x){any(is.na(x))})
sum(fromGenus) 
```

    ## [1] 1083

``` r
complete_traits <-  cbind(complete_traits, fromGenus)
```

## Find and fill in missing data using genera averages

- Next we are going to replace an NA in species trait data information
  with a summary of trait from all species within that genus. That way
  we can make sure we are not missing information from species. We do
  much of what we did above but this is now information only on species
  that are lacking information. So we download all the data for the
  species with that Genus and then summarize it to give an estimate to
  the missing results.

``` r
# Run to check for duplicates
duplicates_check <- which(duplicated(complete_traits$Species))
duplicates_check
```

    ## integer(0)

``` r
# Add genera of species with incomplete data
# Last time I went through this, 1083 species, plus the 48 species only identified to genera had NAs and needed averaged data from the genus.
species_incomplete_genus <- complete_traits %>% dplyr::filter(!complete.cases(.)) %>% dplyr::select(Genus)
species_incomplete_genus <- unique(species_incomplete_genus$Genus)

## Download fishbase data from all species in those additional genera
# Generate the list of all species from those genera
species_incomplete_genus <- unlist(lapply(species_incomplete_genus, function(x) species_list(Genus = x)))

## Extract traits from fishbase
# Main traits
trait_incomplete_genus <- species(species_incomplete_genus)

## Morphology
morphology_incomplete_genus <- morphology(species_incomplete_genus)
# Used to determine duplicates of Species in a column
duplicates_check <- which(duplicated(morphology_incomplete_genus$Species))
duplicates_check 
```

    ## [1] 1585 1586

``` r
# Removed duplicates
morphology_incomplete_genus <- morphology_incomplete_genus %>%
  distinct(Species, .keep_all = TRUE)

# Estimates
estimate_incomplete_genus <- estimate(species_incomplete_genus)

## Ecology
ecology_incomplete_genus <- ecology(species_incomplete_genus)
# Used to determine duplicates of Species in a column
duplicates_check <- which(duplicated(ecology_incomplete_genus$Species))
duplicates_check 
```

    ## [1]   29 2323

``` r
# Removed duplicates
ecology_incomplete_genus <- ecology_incomplete_genus %>%
  distinct(Species, .keep_all = TRUE)

# Reproduction
reproduction_incomplete_genus <- reproduction(species_incomplete_genus)

# Merge the five data.frames by columns
traits_incomplete_genus <- full_join(trait_incomplete_genus, morphology_incomplete_genus, by = "Species")
traits_incomplete_genus <- full_join(traits_incomplete_genus, estimate_incomplete_genus, by = "Species")
traits_incomplete_genus <- full_join(traits_incomplete_genus, ecology_incomplete_genus, by = "Species")
traits_incomplete_genus <- full_join(traits_incomplete_genus, reproduction_incomplete_genus, by = "Species")

# Get rid of duplicated columns
traits_incomplete_genus <- traits_incomplete_genus[,!duplicated(colnames(traits_incomplete_genus))]

# # Calculate the estimated Weight with the following columns.
# traits_incomplete_genus$Weight <- traits_incomplete_genus$a * traits_incomplete_genus$MaxLengthTL^traits_incomplete_genus$b
# 
# # Calculate the estimated caudal fin length with the following columns.
# traits_incomplete_genus$CaudalFinLength <- traits_incomplete_genus$MaxLengthTL - traits_incomplete_genus$MaxLengthSL

# Make -1 in 1 by multiplication
traits_incomplete_genus$OperculumPresent <- traits_incomplete_genus$OperculumPresent*-1

# Filter the fields of interest
traits_incomplete_genus <- traits_incomplete_genus[names(complete_traits)[-ncol(complete_traits)]]

# Check columns to see data types
str(traits_incomplete_genus)
```

    ## tibble [5,450 × 20] (S3: tbl_df/tbl/data.frame)
    ##  $ Species         : chr [1:5450] "Abalistes filamentosus" "Abalistes stellatus" "Ablabys binotatus" "Ablabys gymnothorax" ...
    ##  $ Genus           : chr [1:5450] "Abalistes" "Abalistes" "Ablabys" "Ablabys" ...
    ##  $ BodyShapeI.x    : chr [1:5450] "short and / or deep" "short and / or deep" "short and / or deep" "short and / or deep" ...
    ##  $ Fresh           : int [1:5450] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Brack           : int [1:5450] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Saltwater       : int [1:5450] 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ DemersPelag     : chr [1:5450] "pelagic-neritic" "demersal" "benthopelagic" "benthopelagic" ...
    ##  $ OperculumPresent: num [1:5450] 0 1 0 0 0 0 0 0 0 0 ...
    ##  $ DorsalSpinesMax : int [1:5450] 3 3 15 16 16 17 18 2 2 2 ...
    ##  $ MaxLengthTL     : num [1:5450] 39.7 60 15 10.1 20 ...
    ##  $ Troph           : num [1:5450] 3.38 3.43 3.22 3.17 3.25 ...
    ##  $ DepthMin        : int [1:5450] 61 7 NA NA 1 44 1 3 1 3 ...
    ##  $ DepthMax        : int [1:5450] 180 350 NA 27 20 65 78 2000 80 900 ...
    ##  $ TempPrefMin     : num [1:5450] 18.8 22.9 NA NA 25.1 NA 23.5 19.2 22.8 20 ...
    ##  $ TempPrefMax     : num [1:5450] 26.9 28.3 NA NA 28.9 NA 29 27.9 29 28 ...
    ##  $ FeedingPath     : chr [1:5450] "benthic" "benthic" "benthic" "benthic" ...
    ##  $ RepGuild1       : chr [1:5450] NA "guarders" NA NA ...
    ##  $ RepGuild2       : chr [1:5450] NA "nesters" NA NA ...
    ##  $ ParentalCare    : chr [1:5450] NA NA NA NA ...
    ##  $ DorsalSpinesMin : int [1:5450] 3 3 15 16 15 16 17 2 2 2 ...

``` r
## Calculate genus-level summary traits
# Split species by genera
traits_incomplete_genus_list <- split(traits_incomplete_genus, as.factor(traits_incomplete_genus$Genus))

# Summarize values across all species in the genus. Take the most common character for character traits and take the mean value for numerical traits.
traits_incomplete_genus_summary <- lapply(traits_incomplete_genus_list, function(x){
  c(
  unlist(c(
    sapply(x[,factor], function(y){
      out_y <- table(y)
      ifelse(length(out_y) > 0, names(which.max(out_y)), NA)
    }))),    
  apply(x[,numerical], 2, mean, na.rm=T)
    )
})

# Bind rows into data.frame
complete_traits_incomplete_genus_summary <- data.frame(do.call("rbind", traits_incomplete_genus_summary))

## Check order of columns because they most likely shifted
str(complete_traits_incomplete_genus_summary)
```

    ## 'data.frame':    414 obs. of  18 variables:
    ##  $ BodyShapeI.x    : chr  "short and / or deep" "short and / or deep" "elongated" "elongated" ...
    ##  $ Fresh           : chr  "0" "0" "0" "0" ...
    ##  $ Brack           : chr  "0" "0" "1" "0" ...
    ##  $ Saltwater       : chr  "1" "1" "1" "1" ...
    ##  $ DemersPelag     : chr  "demersal" "benthopelagic" "demersal" "reef-associated" ...
    ##  $ OperculumPresent: chr  "0" "0" "0" "0" ...
    ##  $ FeedingPath     : chr  "benthic" "benthic" "benthic" "benthic" ...
    ##  $ RepGuild1       : chr  "guarders" NA NA "bearers" ...
    ##  $ RepGuild2       : chr  "nesters" NA NA "external brooders" ...
    ##  $ ParentalCare    : chr  NA NA NA "paternal" ...
    ##  $ DorsalSpinesMax : chr  "3" "16.4" "6.8" "0" ...
    ##  $ MaxLengthTL     : chr  "49.8250007629395" "13.2980000495911" "10.0370369663945" "5.65000009536743" ...
    ##  $ Troph           : chr  "3.40500009059906" "3.1960000038147" "3.36296295236658" "3.35666672388713" ...
    ##  $ DepthMin        : chr  "34" "15.3333333333333" "2.35714285714286" "1" ...
    ##  $ DepthMax        : chr  "265" "47.5" "22.625" "33.5" ...
    ##  $ TempPrefMin     : chr  "20.85" "24.3" "25.1636363636364" "24.8" ...
    ##  $ TempPrefMax     : chr  "27.6" "28.95" "29.0909090909091" "29.15" ...
    ##  $ DorsalSpinesMin : chr  "3" "15.8" "6.4" "0" ...

``` r
# Make integer fields numeric
complete_traits_incomplete_genus_summary[,numerical] <- apply(complete_traits_incomplete_genus_summary[,numerical], 2, as.numeric)

# Edit fields
complete_traits_incomplete_genus_summary[,factor] <- lapply(complete_traits_incomplete_genus_summary[,factor], function(x) as.factor(tolower(trimws(x))))

# Remake Genus column
complete_traits_incomplete_genus_summary$Genus <- row.names(complete_traits_incomplete_genus_summary)

# Make new file just in case code doesn't work so you don't have to go back
complete_imputed <- complete_traits

complete_imputed$BodyShapeI.x[is.na(complete_imputed$BodyShapeI.x)] <- complete_traits_incomplete_genus_summary$BodyShapeI.x[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$BodyShapeI.x))]

complete_imputed$Fresh[is.na(complete_imputed$Fresh)] <- complete_traits_incomplete_genus_summary$Fresh[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$Fresh))]

complete_imputed$Brack[is.na(complete_imputed$Brack)] <- complete_traits_incomplete_genus_summary$Brack[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$Brack))]

complete_imputed$Saltwater[is.na(complete_imputed$Saltwater)] <- complete_traits_incomplete_genus_summary$Saltwater[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$Saltwater))]

complete_imputed$DemersPelag[is.na(complete_imputed$DemersPelag)] <- complete_traits_incomplete_genus_summary$DemersPelag[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$DemersPelag))]

complete_imputed$OperculumPresent[is.na(complete_imputed$OperculumPresent)] <- complete_traits_incomplete_genus_summary$OperculumPresent[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$OperculumPresent))]

complete_imputed$DorsalSpinesMin[is.na(complete_imputed$DorsalSpinesMin)] <- complete_traits_incomplete_genus_summary$DorsalSpinesMin[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$DorsalSpinesMin))]

complete_imputed$DorsalSpinesMax[is.na(complete_imputed$DorsalSpinesMax)] <- complete_traits_incomplete_genus_summary$DorsalSpinesMax[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$DorsalSpinesMax))]

complete_imputed$MaxLengthTL[is.na(complete_imputed$MaxLengthTL)] <- complete_traits_incomplete_genus_summary$MaxLengthTL[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$MaxLengthTL))]

complete_imputed$Troph[is.na(complete_imputed$Troph)] <- complete_traits_incomplete_genus_summary$Troph[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$Troph))]

complete_imputed$DepthMin[is.na(complete_imputed$DepthMin)] <- complete_traits_incomplete_genus_summary$DepthMin[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$DepthMin))]

complete_imputed$DepthMax[is.na(complete_imputed$DepthMax)] <- complete_traits_incomplete_genus_summary$DepthMax[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$DepthMax))]

complete_imputed$TempPrefMin[is.na(complete_imputed$TempPrefMin)] <- complete_traits_incomplete_genus_summary$TempPrefMin[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$TempPrefMin))]

complete_imputed$TempPrefMax[is.na(complete_imputed$TempPrefMax)] <- complete_traits_incomplete_genus_summary$TempPrefMax[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$TempPrefMax))]

complete_imputed$FeedingPath[is.na(complete_imputed$FeedingPath)] <- complete_traits_incomplete_genus_summary$FeedingPath[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$FeedingPath))]

complete_imputed$RepGuild1[is.na(complete_imputed$RepGuild1)] <- complete_traits_incomplete_genus_summary$RepGuild1[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$RepGuild1))]

complete_imputed$RepGuild2[is.na(complete_imputed$RepGuild2)] <- complete_traits_incomplete_genus_summary$RepGuild2[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$RepGuild2))]

complete_imputed$ParentalCare[is.na(complete_imputed$ParentalCare)] <- complete_traits_incomplete_genus_summary$ParentalCare[match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)][which(is.na(complete_imputed$ParentalCare))]

# Replace "x" with NA in specified numeric columns
complete_imputed <- complete_imputed %>%
  mutate_at(vars(numerical), ~ ifelse(. == "NaN", NA, .))
```

    ## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
    ## ℹ Please use `all_of()` or `any_of()` instead.
    ##   # Was:
    ##   data %>% select(numerical)
    ## 
    ##   # Now:
    ##   data %>% select(all_of(numerical))
    ## 
    ## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
## Fill in NAs in trait data set "complete_imputed" with genus-level summaries
# complete_imputed <- complete_imputed %>%
#   mutate(across(everything(), ~ ifelse(is.na(.), complete_traits_incomplete_genus_summary[[cur_column()]][match(complete_imputed$Genus, complete_traits_incomplete_genus_summary$Genus)], .)))

str(complete_imputed)
```

    ## 'data.frame':    1675 obs. of  21 variables:
    ##  $ Species         : chr  "Abudefduf lorenzi" "Abudefduf notatus" "Abudefduf septemfasciatus" "Abudefduf sexfasciatus" ...
    ##  $ Genus           : chr  "Abudefduf" "Abudefduf" "Abudefduf" "Abudefduf" ...
    ##  $ BodyShapeI.x    : Factor w/ 5 levels "eel-like","elongated",..: 5 5 5 5 5 5 5 5 5 5 ...
    ##  $ Fresh           : Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Brack           : Factor w/ 2 levels "0","1": 1 1 1 1 2 1 1 1 1 1 ...
    ##  $ Saltwater       : Factor w/ 2 levels "0","1": 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ DemersPelag     : Factor w/ 7 levels "bathydemersal",..: 7 7 7 7 7 7 3 7 7 7 ...
    ##  $ OperculumPresent: Factor w/ 2 levels "0","1": 2 2 2 2 2 2 2 1 2 1 ...
    ##  $ DorsalSpinesMax : num  13 13 13 13 13 13 3 3 3 3 ...
    ##  $ MaxLengthTL     : num  18 17 23 19 24 ...
    ##  $ Troph           : num  2.68 2.73 3.01 2.7 2.69 ...
    ##  $ DepthMin        : num  1 1 0 1 0 1 7 2 10 0 ...
    ##  $ DepthMax        : num  6 12 3 20 3 15 350 21 75 10 ...
    ##  $ TempPrefMin     : num  26 22.9 25.6 24.7 24.7 21.9 22.9 27 24.2 26.3 ...
    ##  $ TempPrefMax     : num  29.4 29.3 29.3 29.3 29.3 29.3 28.3 29.4 29 29.3 ...
    ##  $ FeedingPath     : Factor w/ 2 levels "benthic","pelagic": 1 1 1 1 1 1 1 2 2 2 ...
    ##  $ RepGuild1       : Factor w/ 3 levels "bearers","guarders",..: 2 2 2 2 2 2 2 3 3 3 ...
    ##  $ RepGuild2       : Factor w/ 6 levels "brood hiders",..: 5 5 5 5 5 5 5 6 6 6 ...
    ##  $ ParentalCare    : Factor w/ 4 levels "biparental","maternal",..: 4 4 4 4 4 4 NA 3 3 3 ...
    ##  $ DorsalSpinesMin : num  13 13 13 13 13 13 3 3 3 3 ...
    ##  $ fromGenus       : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...

## Trait data for incomplete species

- For species that were only identified to the genus taxonomic level we
  will use the average trait data for the genus to fill in the missing
  trait data

``` r
#### Add genera from species not in fishbase
#### Use imcomplete_to_add file from fish_species_data.Rmd
### For species not included in fishbase, extract all species in the same genus
## Extract genus names
# Extract out the species that are not present on FishBase
# missing_species <- setdiff(existing_species, existent_species)
# missing_species
incomplete_to_add <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/species_lists/incomplete_to_add.csv")

missing_species <- incomplete_to_add[,-1]

# Delete anything after the space to get the genus level for each species
missing_genus <- str_extract(missing_species, "[^ ]+")
missing_genus
```

    ##  [1] "Acanthurus"     "Apogon"         "Asterropteryx"  "Atherinomorus" 
    ##  [5] "Bathygobius"    "Bryaninops"     "Callogobius"    "Callogobius"   
    ##  [9] "Cercamia"       "Chromis"        "Chromis"        "Chrysiptera"   
    ## [13] "Cirrhilabrus"   "Cryptocentrus"  "Drombus"        "Dussumieria"   
    ## [17] "Epinephelus"    "Eviota"         "Eviota"         "Eviota"        
    ## [21] "Fusigobius"     "Glossogobius"   "Gobiodon"       "Gobiodon"      
    ## [25] "Gorgasia"       "Liopropoma"     "Liopropoma"     "Lubbockichthys"
    ## [29] "Myripristis"    "Ophichthus"     "Opistognathus"  "Pandaka"       
    ## [33] "Parioglossus"   "Platax"         "Plectranthias"  "Pomacentrus"   
    ## [37] "Rhabdoblennius" "Schindleria"    "Sebastapistes"  "Silhouettea"   
    ## [41] "Soleichthys"    "Tetronarce"     "Tomiyamichthys" "Tosanoides"    
    ## [45] "Trimma"         "Trimma"         "Trimmatom"      "Xenisthmus"    
    ## [49] "Xenisthmus"

``` r
# Remove duplicated genera
missing_genus <- missing_genus[!duplicated(missing_genus)]

# Extract all species of the genera in the vector
species_missing_genus_added <- unlist(lapply(missing_genus, function(x) species_list(Genus = x)))

### Download fishbase data from all species in those additional genera
## Extract traits from fishbase
# Main traits
trait_missing_genus_added <- species(species_missing_genus_added)

# Morphology
morphology_missing_genus_added <- morphology(species_missing_genus_added)
# Check if there are duplicated morphologies
duplicates_check <- which(duplicated(morphology_missing_genus_added$Species))
duplicates_check
```

    ## [1] 169

``` r
# Deleted duplicates
morphology_missing_genus_added <- morphology_missing_genus_added %>%
  distinct(Species, .keep_all = TRUE)

# Estimates
estimate_missing_genus_added <- estimate(species_missing_genus_added)

## Ecology
ecology_missing_genus_added <- ecology(species_missing_genus_added)

# Reproduction
reproduction_missing_genus_added <- reproduction(species_missing_genus_added)

# Merge the five trait data.frames
traits_missing_genus_added <- full_join(trait_missing_genus_added, morphology_missing_genus_added, by = "Species")
traits_missing_genus_added <- full_join(traits_missing_genus_added, estimate_missing_genus_added, by = "Species")
traits_missing_genus_added <- full_join(traits_missing_genus_added, ecology_missing_genus_added, by = "Species")
traits_missing_genus_added <- full_join(traits_missing_genus_added, reproduction_missing_genus_added, by = "Species")

# Get rid of duplicated columns
traits_missing_genus_added <- traits_missing_genus_added[,!duplicated(colnames(traits_missing_genus_added))]

# # Calculate the estimated Weight with the following columns.
# traits_missing_genus_added$Weight <- traits_missing_genus_added$a * traits_missing_genus_added$MaxLengthTL^traits_missing_genus_added$b
# 
# # Calculate the estimated caudal fin length with the following columns.
# traits_missing_genus_added$CaudalFinLength <- traits_missing_genus_added$MaxLengthTL - traits_missing_genus_added$MaxLengthSL

# Make -1 in 1 by multiplication
traits_missing_genus_added$OperculumPresent <- traits_missing_genus_added$OperculumPresent*-1

## Filter the fields of interest
traits_missing_genus_added <- traits_missing_genus_added[names(complete_traits)[-ncol(complete_traits)]]

# Check columns to see data types
str(traits_missing_genus_added)
```

    ## tibble [1,308 × 20] (S3: tbl_df/tbl/data.frame)
    ##  $ Species         : chr [1:1308] "Acanthurus achilles" "Acanthurus albimento" "Acanthurus albipectoralis" "Acanthurus auranticavus" ...
    ##  $ Genus           : chr [1:1308] "Acanthurus" "Acanthurus" "Acanthurus" "Acanthurus" ...
    ##  $ BodyShapeI.x    : chr [1:1308] "short and / or deep" "short and / or deep" "fusiform / normal" "short and / or deep" ...
    ##  $ Fresh           : int [1:1308] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Brack           : int [1:1308] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Saltwater       : int [1:1308] 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ DemersPelag     : chr [1:1308] "reef-associated" "reef-associated" "reef-associated" "reef-associated" ...
    ##  $ OperculumPresent: num [1:1308] 1 0 0 1 NA 1 1 1 0 1 ...
    ##  $ DorsalSpinesMax : int [1:1308] 9 9 9 9 NA 9 9 9 8 9 ...
    ##  $ MaxLengthTL     : num [1:1308] 24 30.8 33 54.9 37.7 ...
    ##  $ Troph           : num [1:1308] 2 2.16 3.4 2 2.1 ...
    ##  $ DepthMin        : int [1:1308] 0 NA 5 1 2 6 2 2 0 2 ...
    ##  $ DepthMax        : int [1:1308] 10 NA 20 20 80 50 15 25 6 50 ...
    ##  $ TempPrefMin     : num [1:1308] 25.8 NA 24.7 25.6 NA 24.5 24.9 26.1 NA 23.7 ...
    ##  $ TempPrefMax     : num [1:1308] 28.3 NA 28.4 29.3 NA 28.4 28.8 28 NA 28.1 ...
    ##  $ FeedingPath     : chr [1:1308] "benthic" "benthic" "pelagic" "benthic" ...
    ##  $ RepGuild1       : chr [1:1308] "nonguarders" NA "nonguarders" "nonguarders" ...
    ##  $ RepGuild2       : chr [1:1308] "open water/substratum egg scatterers" NA "open water/substratum egg scatterers" "open water/substratum egg scatterers" ...
    ##  $ ParentalCare    : chr [1:1308] "none" NA "none" "none" ...
    ##  $ DorsalSpinesMin : int [1:1308] 9 9 8 9 NA 9 9 9 8 9 ...

``` r
#### Calculate genus-level summary traits
### Split species by genera
traits_missing_genus_added_list <- split(traits_missing_genus_added, as.factor(traits_missing_genus_added$Genus))

### Summarize values across all species in the genus
traits_missing_genus_added_summary <- lapply(traits_missing_genus_added_list, function(x){
  c(
  unlist(c(
    sapply(x[,factor], function(y){
      out_y <- table(y)
      ifelse(length(out_y) > 0, names(which.max(out_y)), NA)
    }))),    
  apply(x[,numerical], 2, mean, na.rm=T)
    )
})

### Bind rows into data.frame
complete_traits_missing_genus_added_summary <- data.frame(do.call("rbind", traits_missing_genus_added_summary))

# Make integer fields numeric
complete_traits_missing_genus_added_summary[,numerical] <- sapply(complete_traits_missing_genus_added_summary[,numerical], as.numeric)

# Edit fields
complete_traits_missing_genus_added_summary[,factor] <- lapply(complete_traits_missing_genus_added_summary[,factor], function(x) as.factor(tolower(trimws(x))))

## Check order of columns because they most likely shifted
str(complete_traits_missing_genus_added_summary)
```

    ## 'data.frame':    41 obs. of  18 variables:
    ##  $ BodyShapeI.x    : Factor w/ 5 levels "eel-like","elongated",..: 5 3 2 2 2 2 2 5 5 3 ...
    ##  $ Fresh           : Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Brack           : Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Saltwater       : Factor w/ 2 levels "0","1": 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ DemersPelag     : Factor w/ 4 levels "benthopelagic",..: 4 4 4 4 2 4 4 4 4 4 ...
    ##  $ OperculumPresent: Factor w/ 2 levels "0","1": 2 1 1 2 1 1 1 1 2 2 ...
    ##  $ FeedingPath     : Factor w/ 2 levels "benthic","pelagic": 1 1 1 1 1 2 1 1 2 1 ...
    ##  $ RepGuild1       : Factor w/ 3 levels "bearers","guarders",..: 3 1 2 NA 2 NA NA 1 2 2 ...
    ##  $ RepGuild2       : Factor w/ 5 levels "clutch tenders",..: 5 2 NA NA NA NA NA NA 4 4 ...
    ##  $ ParentalCare    : Factor w/ 2 levels "none","paternal": 1 2 NA NA NA NA NA NA 2 2 ...
    ##  $ DorsalSpinesMax : num  8.86 7.08 7 7.78 7 ...
    ##  $ MaxLengthTL     : num  37.7 8.6 3.81 13.46 9.9 ...
    ##  $ Troph           : num  2.28 3.45 3.26 3.27 3.45 ...
    ##  $ DepthMin        : num  2.667 5.86 12.333 0.167 0.105 ...
    ##  $ DepthMax        : num  51.4 46.5 39.2 31.8 15 ...
    ##  $ TempPrefMin     : num  24.7 24 25.8 25.9 24.2 ...
    ##  $ TempPrefMax     : num  28.8 28.3 29.1 28.9 28.9 ...
    ##  $ DorsalSpinesMin : num  8.78 7.08 6.89 5.89 6.89 ...

``` r
# Make Genus column
complete_traits_missing_genus_added_summary$Genus <- row.names(complete_traits_missing_genus_added_summary)
complete_traits_missing_genus_added_summary <- complete_traits_missing_genus_added_summary[!duplicated(complete_traits_missing_genus_added_summary[c('Genus')]), ]

# Duplicate the missing Genus rows that were deleted earlier
complete_traits_missing_genus_added_summary <- rbind(complete_traits_missing_genus_added_summary, complete_traits_missing_genus_added_summary[c("Callogobius", "Chromis", "Eviota", "Eviota", "Gobiodon", "Liopropoma", "Trimma", "Xenisthmus"), ])

# Re-order rows to organize duplicated rows
complete_traits_missing_genus_added_summary <- complete_traits_missing_genus_added_summary[order(row.names(complete_traits_missing_genus_added_summary)), ]

### Add species and genus names
complete_traits_missing_genus_added_summary$Species <- sort(missing_species)
complete_traits_missing_genus_added_summary$fromGenus <- TRUE

### Add final species to trait data.frame 
complete_imputed <- data.frame(rbind(complete_imputed, complete_traits_missing_genus_added_summary))

# Identify which rows have NAs and will be filled in with Genus level data
fromFamily <- apply(complete_imputed, 1, function(x){any(is.na(x))})
sum(fromFamily)
```

    ## [1] 723

``` r
complete_imputed_genus <-  cbind(complete_imputed, fromFamily)
```

## Higher taxonomic classifications

- We need higher taxonomic information for each species to fill in
  missing data using family trait averages.

``` r
# Added from Family field
#complete_imputed_genus$fromFamily <- 0
genera <- complete_imputed_genus$Genus

# Download taxonomic information from FishBase
genera <- rfishbase::load_taxa() %>% 
  filter(Genus %in% genera) %>%
  collect()
```

    ## Joining with `by = join_by(Subfamily, GenCode, FamCode)`
    ## Joining with `by = join_by(FamCode)`
    ## Joining with `by = join_by(Order, Ordnum, Class, ClassNum)`
    ## Joining with `by = join_by(Class, ClassNum)`

``` r
# Choose which classification levels you want to keep
keep <- c("Genus", "Family", "Order", "Class")
genera <- genera[,keep]

# Merge dataframes together
complete_imputed_genus <- merge(genera, complete_imputed_genus, by = "Genus", all.x = F, all.y = F)

# Remove duplicates in dataframe
complete_imputed_genus <- complete_imputed_genus[!duplicated(complete_imputed_genus[c('Species')]), ]

# Move family column to the first column of the dataframe
complete_imputed_genus <- complete_imputed_genus %>% relocate(Family, .before = Genus) %>% 
  relocate(Order, .before = Family) %>% 
  relocate(Class, .before = Order)
```

## Find and fill in missing data using family averages

- Next we are going to replace an NA in species trait data information
  with a summary of trait from all species within that family. That way
  we can make sure we are not missing information from species. We do
  much of what we did above but this is now information only on species
  that are lacking information after the Genus level summary. So we
  download all the data for the species with that Family and then
  summarize it to give an estimate to the missing results.

``` r
# Run to check for duplicates
duplicates_check <- which(duplicated(complete_imputed_genus$Species))
duplicates_check
```

    ## integer(0)

``` r
# Add genera of species with incomplete data
# Last time I went through this 724 species had NAs and needed data from the family level
species_incomplete_family <- complete_imputed_genus %>% dplyr::filter(!complete.cases(.)) %>% dplyr::select(Family)
species_incomplete_family <- unique(species_incomplete_family$Family)

## Download fishbase data from all species in those additional genera
# Generate the list of all species from those genera
species_incomplete_family <- unlist(lapply(species_incomplete_family, function(x) species_list(Family = x)))

# Get rid of everything after the space to get genera only
species_names_family <- str_extract(species_incomplete_family, "[^ ]+")

# Download taxonomic information from FishBase
species_names_family <- rfishbase::load_taxa() %>% 
  filter(Genus %in% species_names_family) %>%
  collect()

# Choose which classification levels you want to keep
keep <- c("Species", "Family")
species_names_family <- species_names_family[,keep]

## Extract traits from fishbase
# Main traits
trait_incomplete_family <- species(species_incomplete_family)

## Morphology
morphology_incomplete_family <- morphology(species_incomplete_family)
# Used to determine duplicates of Species in a column
duplicates_check <- which(duplicated(morphology_incomplete_family$Species))
duplicates_check 
```

    ## [1] 2465 2466

``` r
# Removed duplicates
morphology_incomplete_family <- morphology_incomplete_family %>%
  distinct(Species, .keep_all = TRUE)

# Estimates
estimate_incomplete_family <- estimate(species_incomplete_family)

## Ecology
ecology_incomplete_family <- ecology(species_incomplete_family)
# Used to determine duplicates of Species in a column
duplicates_check <- which(duplicated(ecology_incomplete_family$Species))
duplicates_check 
```

    ## [1] 3245

``` r
# Removed duplicates
ecology_incomplete_family <- ecology_incomplete_family %>%
  distinct(Species, .keep_all = TRUE)

# Reproduction
reproduction_incomplete_family <- reproduction(species_incomplete_family)

# Merge the five data.frames by columns
traits_incomplete_family <- full_join(trait_incomplete_family, morphology_incomplete_family, by = "Species")
traits_incomplete_family <- full_join(traits_incomplete_family, estimate_incomplete_family, by = "Species")
traits_incomplete_family <- full_join(traits_incomplete_family, ecology_incomplete_family, by = "Species")
traits_incomplete_family <- full_join(traits_incomplete_family, reproduction_incomplete_family, by = "Species")

# Get rid of duplicated columns
traits_incomplete_family <- traits_incomplete_family[,!duplicated(colnames(traits_incomplete_family))]

# # Calculate the estimated Weight with the following columns.
# traits_incomplete_family$Weight <- traits_incomplete_family$a * traits_incomplete_family$MaxLengthTL^traits_incomplete_family$b
# 
# # Calculate the estimated caudal fin length with the following columns.
# traits_incomplete_family$CaudalFinLength <- traits_incomplete_family$MaxLengthTL - traits_incomplete_family$MaxLengthSL

# Make -1 in 1 by multiplication
traits_incomplete_family$OperculumPresent <- traits_incomplete_family$OperculumPresent*-1

# Filter the fields of interest
traits_incomplete_family <- traits_incomplete_family[names(complete_traits)[-ncol(complete_traits)]]

# Check columns to see data types
str(traits_incomplete_family)
```

    ## tibble [8,861 × 20] (S3: tbl_df/tbl/data.frame)
    ##  $ Species         : chr [1:8861] "Abalistes filamentosus" "Abalistes stellatus" "Ablabys binotatus" "Ablabys gymnothorax" ...
    ##  $ Genus           : chr [1:8861] "Abalistes" "Abalistes" "Ablabys" "Ablabys" ...
    ##  $ BodyShapeI.x    : chr [1:8861] "short and / or deep" "short and / or deep" "short and / or deep" "short and / or deep" ...
    ##  $ Fresh           : int [1:8861] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Brack           : int [1:8861] 0 0 0 0 0 0 0 0 0 1 ...
    ##  $ Saltwater       : int [1:8861] 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ DemersPelag     : chr [1:8861] "pelagic-neritic" "demersal" "benthopelagic" "benthopelagic" ...
    ##  $ OperculumPresent: num [1:8861] 0 1 0 0 0 0 0 0 NA 0 ...
    ##  $ DorsalSpinesMax : int [1:8861] 3 3 15 16 16 17 18 NA NA 2 ...
    ##  $ MaxLengthTL     : num [1:8861] 39.7 60 15 10.1 20 ...
    ##  $ Troph           : num [1:8861] 3.38 3.43 3.22 3.17 3.25 ...
    ##  $ DepthMin        : int [1:8861] 61 7 NA NA 1 44 1 NA 1 1 ...
    ##  $ DepthMax        : int [1:8861] 180 350 NA 27 20 65 78 NA 50 50 ...
    ##  $ TempPrefMin     : num [1:8861] 18.8 22.9 NA NA 25.1 NA 23.5 NA 15.2 14.3 ...
    ##  $ TempPrefMax     : num [1:8861] 26.9 28.3 NA NA 28.9 NA 29 NA 18.4 21 ...
    ##  $ FeedingPath     : chr [1:8861] "benthic" "benthic" "benthic" "benthic" ...
    ##  $ RepGuild1       : chr [1:8861] NA "guarders" NA NA ...
    ##  $ RepGuild2       : chr [1:8861] NA "nesters" NA NA ...
    ##  $ ParentalCare    : chr [1:8861] NA NA NA NA ...
    ##  $ DorsalSpinesMin : int [1:8861] 3 3 15 16 15 16 17 NA NA 2 ...

``` r
# Add the family column into the dataframe
traits_incomplete_family <- merge(species_names_family, traits_incomplete_family, by = "Species", all.x = F)

## Calculate Family-level summary traits
# Split species by genera
traits_incomplete_family_list <- split(traits_incomplete_family, as.factor(traits_incomplete_family$Family))

# Summarize values across all species in the Family. Take the most common character for character traits and take the mean value for numerical traits.
traits_incomplete_family_summary <- lapply(traits_incomplete_family_list, function(x){
  c(
  unlist(c(
    sapply(x[,factor], function(y){
      out_y <- table(y)
      ifelse(length(out_y) > 0, names(which.max(out_y)), NA)
    }))),    
  apply(x[,numerical], 2, mean, na.rm=T)
    )
})

# Bind rows into data.frame
complete_traits_incomplete_family_summary <- data.frame(do.call("rbind", traits_incomplete_family_summary))

## Check order of columns because they most likely shifted
str(complete_traits_incomplete_family_summary)
```

    ## 'data.frame':    92 obs. of  18 variables:
    ##  $ BodyShapeI.x    : chr  "other" "fusiform / normal" "fusiform / normal" "Elongated" ...
    ##  $ Fresh           : chr  "0" "0" "0" "0" ...
    ##  $ Brack           : chr  "1" "0" "0" "0" ...
    ##  $ Saltwater       : chr  "1" "1" "1" "1" ...
    ##  $ DemersPelag     : chr  "benthopelagic" "demersal" "reef-associated" "demersal" ...
    ##  $ OperculumPresent: chr  "0" "0" "0" "0" ...
    ##  $ FeedingPath     : chr  "benthic" "pelagic" "pelagic" "pelagic" ...
    ##  $ RepGuild1       : chr  "bearers" NA "nonguarders" NA ...
    ##  $ RepGuild2       : chr  "internal live bearers" NA "open water/substratum egg scatterers" NA ...
    ##  $ ParentalCare    : chr  NA NA "none" NA ...
    ##  $ DorsalSpinesMax : chr  "0" "5.14285714285714" "10.1307692307692" "12.6875" ...
    ##  $ MaxLengthTL     : chr  "143.910003662109" "17.1455554962158" "17.2376595659459" "6.99795919048543" ...
    ##  $ Troph           : chr  "3.72999992370605" "3.45000002119276" "3.53719149447502" "3.31448981713276" ...
    ##  $ DepthMin        : chr  "0.666666666666667" "32.375" "56.3097826086956" "10.8085106382979" ...
    ##  $ DepthMax        : chr  "76.6666666666667" "235.333333333333" "160.490476190476" "56.8541666666667" ...
    ##  $ TempPrefMin     : chr  "23.7666666666667" "21.7333333333333" "18.6956204379562" "23.1911764705882" ...
    ##  $ TempPrefMax     : chr  "29.0666666666667" "27.8333333333333" "24.4510948905109" "27.6529411764706" ...
    ##  $ DorsalSpinesMin : chr  "0" "4.42857142857143" "10.0923076923077" "11.90625" ...

``` r
# Make integer fields numeric
complete_traits_incomplete_family_summary[,numerical] <- apply(complete_traits_incomplete_family_summary[,numerical], 2, as.numeric)

# Edit fields
complete_traits_incomplete_family_summary[,factor] <- lapply(complete_traits_incomplete_family_summary[,factor], function(x) as.factor(tolower(trimws(x))))

# Remake Family column
complete_traits_incomplete_family_summary$Family <- row.names(complete_traits_incomplete_family_summary)

# Make new file just in case code doesn't work so you don't have to go back
complete_imputed_family <- complete_imputed_genus

complete_imputed_family$BodyShapeI.x[is.na(complete_imputed_family$BodyShapeI.x)] <- complete_traits_incomplete_family_summary$BodyShapeI.x[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$BodyShapeI.x))]

complete_imputed_family$Fresh[is.na(complete_imputed_family$Fresh)] <- complete_traits_incomplete_family_summary$Fresh[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$Fresh))]

complete_imputed_family$Brack[is.na(complete_imputed_family$Brack)] <- complete_traits_incomplete_family_summary$Brack[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$Brack))]

complete_imputed_family$Saltwater[is.na(complete_imputed_family$Saltwater)] <- complete_traits_incomplete_family_summary$Saltwater[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$Saltwater))]

complete_imputed_family$DemersPelag[is.na(complete_imputed_family$DemersPelag)] <- complete_traits_incomplete_family_summary$DemersPelag[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$DemersPelag))]

complete_imputed_family$OperculumPresent[is.na(complete_imputed_family$OperculumPresent)] <- complete_traits_incomplete_family_summary$OperculumPresent[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$OperculumPresent))]

complete_imputed_family$DorsalSpinesMin[is.na(complete_imputed_family$DorsalSpinesMin)] <- complete_traits_incomplete_family_summary$DorsalSpinesMin[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$DorsalSpinesMin))]

complete_imputed_family$DorsalSpinesMax[is.na(complete_imputed_family$DorsalSpinesMax)] <- complete_traits_incomplete_family_summary$DorsalSpinesMax[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$DorsalSpinesMax))]

complete_imputed_family$MaxLengthTL[is.na(complete_imputed_family$MaxLengthTL)] <- complete_traits_incomplete_family_summary$MaxLengthTL[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$MaxLengthTL))]

complete_imputed_family$Troph[is.na(complete_imputed_family$Troph)] <- complete_traits_incomplete_family_summary$Troph[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$Troph))]

complete_imputed_family$DepthMin[is.na(complete_imputed_family$DepthMin)] <- complete_traits_incomplete_family_summary$DepthMin[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$DepthMin))]

complete_imputed_family$DepthMax[is.na(complete_imputed_family$DepthMax)] <- complete_traits_incomplete_family_summary$DepthMax[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$DepthMax))]

complete_imputed_family$TempPrefMin[is.na(complete_imputed_family$TempPrefMin)] <- complete_traits_incomplete_family_summary$TempPrefMin[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$TempPrefMin))]

complete_imputed_family$TempPrefMax[is.na(complete_imputed_family$TempPrefMax)] <- complete_traits_incomplete_family_summary$TempPrefMax[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$TempPrefMax))]

complete_imputed_family$FeedingPath[is.na(complete_imputed_family$FeedingPath)] <- complete_traits_incomplete_family_summary$FeedingPath[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$FeedingPath))]

complete_imputed_family$RepGuild1[is.na(complete_imputed_family$RepGuild1)] <- complete_traits_incomplete_family_summary$RepGuild1[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$RepGuild1))]

complete_imputed_family$RepGuild2[is.na(complete_imputed_family$RepGuild2)] <- complete_traits_incomplete_family_summary$RepGuild2[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$RepGuild2))]

complete_imputed_family$ParentalCare[is.na(complete_imputed_family$ParentalCare)] <- complete_traits_incomplete_family_summary$ParentalCare[match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Family)][which(is.na(complete_imputed_family$ParentalCare))]

## Fill in NAs in trait data set "complete_imputed" with genus-level summaries
# complete_imputed_family <- complete_imputed_family %>%
#   mutate(across(everything(), ~ ifelse(is.na(.), complete_traits_incomplete_family_summary[[cur_column()]][match(complete_imputed_family$Family, complete_traits_incomplete_family_summary$Genus)], .)))
```

## Last file checks

- Edit rows and fields to make sure that final complete_imputed file has
  the necessary format for analyses. Make sure certain variables are
  numeric and others are factors. You can run the str command at the end
  to be sure everything is in the correct format.

``` r
## Final editing steps
# Change species names that were altered by FishBase to the subspecies level because we only want species level information 
# complete$Species[complete$Species=="Platybelone argalus platyura"] <- "Platybelone argalus"
# complete$Species[complete$Species=="Tylosurus acus melanotus"] <- "Tylosurus melanotus"

complete <- complete_imputed_family

# Re-order rows
complete <- complete[order(complete$Species), ]

# Erase rownmaes since they need to be updated
rownames(complete) <- NULL

duplicates_check <- which(duplicated(complete$Species))
duplicates_check
```

    ## integer(0)

``` r
### Rename rows
complete$Species <- gsub(" ", "_", complete$Species)
row.names(complete) <- complete$Species
# complete <- complete[,-1]

# Isolate the most complete traits (completeness > 99.9%), you will lose traits here
complete_best_traits <- complete[, which(completeness(complete) > 99.9)]

### Select all traits with more than 98% coverage, you will lose species and traits here
complete_best_species <- complete[, which(completeness(complete) >= 98)] %>% dplyr::filter(complete.cases(.))

# Check your final species count with the species list you started with
complete_species <- complete[,5]
existing_species <- gsub(" ", "_", existing_species)
diff_check <- setdiff(existing_species, complete_species)
diff_check
```

    ## character(0)

## Load out files

``` r
## Load out files
# Traits for fish that are present in the marine lakes and ocean populations
# Has all or most species but loses traits. Species with 0 or very few NAs
write.csv(complete_best_traits,  "/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/traits/complete_best_all_species_gnathostomata.csv")

# Has all or most traits but loses species. Traits with 0 or very few NAs
write.csv(complete_best_species,  "/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/traits/complete_best_all_traits_gnathostomata.csv")

#Traits for fish that are present in Palau
# Traits with more NAs
# Edit species names, replace space with an underscore
write.csv(complete, "/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/traits/complete_traits_gnathostomata.csv")

sessionInfo()
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.6.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] stringr_1.5.1   dplyr_1.1.4     tidyr_1.3.0     rfishbase_4.1.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] jsonlite_1.8.8    compiler_4.3.1    crayon_1.5.2      tidyselect_1.2.0 
    ##  [5] progress_1.2.3    yaml_2.3.7        fastmap_1.1.1     readr_2.1.4      
    ##  [9] R6_2.5.1          generics_0.1.3    curl_5.2.0        knitr_1.45       
    ## [13] tibble_3.2.1      openssl_2.1.1     DBI_1.1.3         pillar_1.9.0     
    ## [17] tzdb_0.4.0        rlang_1.1.2       utf8_1.2.4        cachem_1.0.8     
    ## [21] stringi_1.8.2     xfun_0.41         fs_1.6.3          memoise_2.0.1    
    ## [25] cli_3.6.1         withr_2.5.2       magrittr_2.0.3    digest_0.6.33    
    ## [29] rstudioapi_0.15.0 dbplyr_2.4.0      askpass_1.2.0     hms_1.1.3        
    ## [33] lifecycle_1.0.4   prettyunits_1.2.0 vctrs_0.6.5       evaluate_0.23    
    ## [37] glue_1.6.2        duckdb_0.9.2-1    contentid_0.0.18  fansi_1.0.6      
    ## [41] httr_1.4.7        rmarkdown_2.25    purrr_1.0.2       tools_4.3.1      
    ## [45] pkgconfig_2.0.3   htmltools_0.5.7

## Post R processing

- Before proceeding to analyses with your trait file
  complete_traits_gnathostomata.csv you need to do some additional work
  with the dorsal spines and reproductive traits. Review the current
  traits NA values and use FishBase to see if the family page for the
  species contains any information you can use to fill in these NA
  values. Sometimes the R package rfishbase is not able to see some of
  the family information on the website. You find the family page for
  species by searching the species name (<https://www.fishbase.us/>) and
  then clicking the family name under the “Classification / Names”
  section on the page which should be the first section. After using
  FishBase to complete the dorsal spines and reproductive information
  your file will be named
  complete_traits_gnathostomata_with_FishBase_corrections.csv

- Then you use the book Fish Reproduction by Rocha et al. (2008) and use
  chapter 9, Reproductive Strategies of Fish by Patzner, R. A. to fill
  in some of the blanks of the RepGuild2 column. You will only use
  families where the RepGuild2 is the same for all species and fits with
  the categorization in our table already. Mainly the “open
  water/substratum egg scatterers”. You will then name your file
  complete_traits_gnathostomata_with_Patzner_2008_data.csv

- After that you can determine the traits that you would like to include
  in your analyses which may be all or a subset. For this project we
  chose a subset that have quality definitions according to FishBase and
  are not repetitive traits. We named our final file
  final_traits_gnathostomata.csv which we use in all downstream
  analyses.
