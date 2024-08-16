Stratification NMDS envfit analyses
================

#### R Markdown

## Load stratification NMDS envfit analyses packages

``` r
# Load the knitr package if not already loaded
library(knitr)

# Source the R Markdown files
knit("/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.Rmd", output = "/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md")
```

    ## 
    ## 
    ## processing file: /Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.Rmd

    ##   |                  |          |   0%  |                  |          |   3%                                                                          |                  |.         |   6% [Bringing everything together load modifying files packages]             |                  |.         |   9%                                                                          |                  |.         |  12% [Bringing everything together load in modifying files]                   |                  |..        |  16%                                                                          |                  |..        |  19% [Check species names across files]                                       |                  |..        |  22%                                                                          |                  |..        |  25% [Modify environment data]                                                |                  |...       |  28%                                                                          |                  |...       |  31% [Modify incidence matrices]                                              |                  |...       |  34%                                                                          |                  |....      |  38% [Modify phylogeny]                                                       |                  |....      |  41%                                                                          |                  |....      |  44% [Modify trait data]                                                      |                  |.....     |  47%                                                                          |                  |.....     |  50% [Modify community data frames]                                           |                  |.....     |  53%                                                                          |                  |......    |  56% [Modify community trait data]                                            |                  |......    |  59%                                                                          |                  |......    |  62% [Community trait data tests]                                             |                  |.......   |  66%                                                                          |                  |.......   |  69% [Modify stratification data frames]                                      |                  |.......   |  72%                                                                          |                  |........  |  75% [Modify stratification trait data]                                       |                  |........  |  78%                                                                          |                  |........  |  81% [Stratification trait data tests]                                        |                  |........  |  84%                                                                          |                  |......... |  88% [Modify site trait data frames]                                          |                  |......... |  91%                                                                          |                  |......... |  94% [Site trait data tests]                                                  |                  |..........|  97%                                                                          |                  |..........| 100% [Bringing everything together load out modified files and session info]

    ## output file: /Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md

    ## [1] "/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md"

``` r
knit("/Users/bailey/Documents/research/fish_biodiversity/src/analyses/dbFD_results.Rmd", output = "/Users/bailey/Documents/research/fish_biodiversity/src/analyses/dbFD_results.md")
```

    ## 
    ## 
    ## processing file: /Users/bailey/Documents/research/fish_biodiversity/src/analyses/dbFD_results.Rmd

    ##   |                           |                   |   0%  |                           |.                  |   6%                                                     |                           |..                 |  11% [FD analyses results load packages and files]       |                           |...                |  17%                                                     |                           |....               |  22% [Functional Diversity FRic, FEve, FDiv, FDis, CWM]  |                           |.....              |  28%                                                     |                           |......             |  33% [Edit FD files]                                     |                           |.......            |  39%                                                     |                           |........           |  44% [Biplot of FRic and FDis]                           |                           |..........         |  50%                                                     |                           |...........        |  56% [FD analyses results session info]                  |                           |............       |  61%                                                     |                           |.............      |  67% [Null models of functional diversity]               |                           |..............     |  72%                                                     |                           |...............    |  78% [CWM with null]                                     |                           |................   |  83%                                                     |                           |.................  |  89% [FRic with null]                                    |                           |.................. |  94%                                                     |                           |...................| 100% [FDis with null]                                  

    ## output file: /Users/bailey/Documents/research/fish_biodiversity/src/analyses/dbFD_results.md

    ## [1] "/Users/bailey/Documents/research/fish_biodiversity/src/analyses/dbFD_results.md"

``` r
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

#Analyses
library(vegan)
library(pairwiseAdonis)
```

    ## Loading required package: cluster

    ## 
    ## Attaching package: 'cluster'

    ## The following object is masked from 'package:maps':
    ## 
    ##     votes.repub

``` r
library(picante)
```

    ## Loading required package: nlme

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

``` r
library(ade4)
library(adespatial)
```

    ## Registered S3 methods overwritten by 'adegraphics':
    ##   method         from
    ##   biplot.dudi    ade4
    ##   kplot.foucart  ade4
    ##   kplot.mcoa     ade4
    ##   kplot.mfa      ade4
    ##   kplot.pta      ade4
    ##   kplot.sepan    ade4
    ##   kplot.statis   ade4
    ##   scatter.coa    ade4
    ##   scatter.dudi   ade4
    ##   scatter.nipals ade4
    ##   scatter.pco    ade4
    ##   score.acm      ade4
    ##   score.mix      ade4
    ##   score.pca      ade4
    ##   screeplot.dudi ade4

    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape

    ## Registered S3 methods overwritten by 'adespatial':
    ##   method             from       
    ##   plot.multispati    adegraphics
    ##   print.multispati   ade4       
    ##   summary.multispati ade4

    ## 
    ## Attaching package: 'adespatial'

    ## The following object is masked from 'package:ade4':
    ## 
    ##     multispati

``` r
library(car)
```

    ## Loading required package: carData

    ## 
    ## Attaching package: 'car'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode

``` r
library(rms)
```

    ## Warning: package 'rms' was built under R version 4.3.3

    ## Loading required package: Hmisc

    ## 
    ## Attaching package: 'Hmisc'

    ## The following object is masked from 'package:ape':
    ## 
    ##     zoom

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     src, summarize

    ## The following objects are masked from 'package:base':
    ## 
    ##     format.pval, units

    ## 
    ## Attaching package: 'rms'

    ## The following objects are masked from 'package:car':
    ## 
    ##     Predict, vif

    ## The following object is masked from 'package:vegan':
    ## 
    ##     calibrate

``` r
#Graphing
library(ggplot2)
library(viridis)
library(ggrepel)

stratification_group <- env[-c(20),19]

surveyed_sites <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "LCN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOO", "OTM", "OOM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_wo_LCN <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOO", "OTM", "OOM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_env <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OTM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_geo <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LCN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOO", "OTM", "OOM", "RCA", "SLN", "TLN", "ULN")

ocean_mixed_sites <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA", "FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")

ocean_mixed_sites_env <- c("IBK", "NCN", "RCA", "FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")

ocean_mixed_sites_geo <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA", "FLK", "HLO", "MLN", "NLN", "NLU", "OLO", "ULN")

ocean_stratified_sites <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA", "BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")

ocean_stratified_sites_env <- c("IBK", "NCN", "RCA", "BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")

mixed_stratified_lakes <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "LLN", "MLN", "NLK", "NLN", "NLU", "OLO", "OTM", "SLN", "TLN", "ULN")

mixed_stratified_lakes_geo <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "MLN", "NLK", "NLN", "NLU", "OLO", "OTM", "SLN", "TLN", "ULN")

stratified_lakes <- c("BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")

# Define your custom colors
custom_colors <- c("Reference" = "black", "Ocean" = "#EE6363", "Mixed" = "#87CEFA", "Stratified" = "#6E8B3D")

env_cont <- "#FFB90F"

trait_cont <- "#68228B"

trait_cat <- "#DEB887"
```

## Load stratification NMDS envfit analyses functions

- Need to run this function for downstream analyses in this Markdown
  file

``` r
# Function p.adjust.envfit
# Calculates adjusted P values for results stored in envfit object,
# created using envfit function from vegan, which fits supplementary variables
# onto axes of unconstrained ordination. 
# Arguments: 
# x - envfit object
# method - method for correction of multiple testing issue (default = 'bonferroni',
#          see ''?p.adjust' for more options)
# n - optional, number of tests for which to correct; if not given, the number is
#          taken as the number of tests conducted by envfit function (for both vectors and factors).
# Author: David Zeleny
p.adjust.envfit <- function (x, method = 'bonferroni', n)
{
  x.new <- x
  if (!is.null (x$vectors)) pval.vectors <- x$vectors$pvals else pval.vectors <- NULL
  if (!is.null (x$factors)) pval.factors <- x$factors$pvals else pval.factors <- NULL
  if (missing (n)) n <- length (pval.vectors) + length (pval.factors)
  if (!is.null (x$vectors)) x.new$vectors$pvals <- p.adjust (x$vectors$pvals, method = method, n = n)
  if (!is.null (x$factors)) x.new$factors$pvals <- p.adjust (x$factors$pvals, method = method, n = n)
  cat ('Adjustment of significance by', method, 'method')
  return (x.new)
}
```

## Incidence dissimilarity distances

``` r
# 1.BCM (S), 2.CLM (S), 3.FLK (M), 4.GLK (S), 5.HLM (S), 6.HLO (M), 7.IBK (O), 8.LLN (M), 9.LLNC (O), 10.MLN (M), 11.NCN (O), 12.NLK (S), 13.NLN (M), 14.NLU (M), 15.OLO (M), 16.OLOO (O), 17.OTM (S), 18.OTMO (O), 19.RCA (O), 20.REF (R), 21.SLN (S), 22.TLN (S), 23.ULN (M), 24.Ocean sites, 25.Mixed lakes, 26.Stratified lakes

### Regular
## Jaccard
# All sites
jpres_dist <- vegdist(presabs_lake[surveyed_sites,], method = "jaccard", binary = TRUE)
# Mixed and stratified lakes
jpresms_dist <- vegdist(presabs_lake[mixed_stratified_lakes,], method = "jaccard", binary = TRUE)
# Ocean sites and mixed lakes
jpresom_dist <- vegdist(presabs_lake[ocean_mixed_sites,], method = "jaccard", binary = TRUE)
# Stratified lakes and ocean sites
jpresso_dist <- vegdist(presabs_lake[ocean_stratified_sites,], method = "jaccard", binary = TRUE)
# Stratified lakes 
jpress_dist <- vegdist(presabs_lake[stratified_lakes,], method = "jaccard", binary = TRUE)

## Dice-Sørensen
# All sites
spres_dist <- vegdist(presabs_lake[surveyed_sites,], method = "bray", binary = TRUE)
# Mixed and stratified lakes
spresms_dist <- vegdist(presabs_lake[mixed_stratified_lakes,], method = "bray", binary = TRUE)
# Ocean sites and mixed lakes
spresom_dist <- vegdist(presabs_lake[ocean_mixed_sites,], method = "bray", binary = TRUE)
# Stratified lakes and ocean sites
spresso_dist <- vegdist(presabs_lake[ocean_stratified_sites,], method = "bray", binary = TRUE)


### Environmental
## Jaccard
# All sites
jprese_dist <- vegdist(presabs_lake[surveyed_sites_env,], method = "jaccard", binary = TRUE)
# Mixed and stratified lakes
jpresems_dist <- vegdist(presabs_lake[mixed_stratified_lakes,], method = "jaccard", binary = TRUE)
# Ocean sites and mixed lakes
jpreseom_dist <- vegdist(presabs_lake[ocean_mixed_sites_env,], method = "jaccard", binary = TRUE)
# Stratified lakes and ocean sites
jpreseso_dist <- vegdist(presabs_lake[ocean_stratified_sites_env,], method = "jaccard", binary = TRUE)

## Dice-Sørensen
# All sites
sprese_dist <- vegdist(presabs_lake[surveyed_sites_env,], method = "bray", binary = TRUE)
# Mixed and stratified lakes
spresems_dist <- vegdist(presabs_lake[mixed_stratified_lakes,], method = "bray", binary = TRUE)
# Ocean sites and mixed lakes
spreseom_dist <- vegdist(presabs_lake[ocean_mixed_sites_env,], method = "bray", binary = TRUE)
# Stratified lakes and ocean sites
spreseso_dist <- vegdist(presabs_lake[ocean_stratified_sites_env,], method = "bray", binary = TRUE)


### Geographic
## Jaccard
# All sites
jpresb_dist <- vegdist(presabs_lake[surveyed_sites_geo,], method = "jaccard", binary = TRUE)
# Mixed and stratified lakes
jpresbms_dist <- vegdist(presabs_lake[mixed_stratified_lakes_geo,], method = "jaccard", binary = TRUE)
# Ocean sites and mixed lakes 
jpresbom_dist <- vegdist(presabs_lake[ocean_mixed_sites_geo,], method = "jaccard", binary = TRUE)
# Stratified lakes and ocean sites
jpresbso_dist <- vegdist(presabs_lake[ocean_stratified_sites,], method = "jaccard", binary = TRUE)

## Dice-Sørensen
# All sites
spresb_dist <- vegdist(presabs_lake[surveyed_sites_geo,], method = "bray", binary = TRUE)
# Mixed and stratified lakes
spresbms_dist <- vegdist(presabs_lake[mixed_stratified_lakes_geo,], method = "bray", binary = TRUE)
# Ocean sites and mixed lakes
spresbom_dist <- vegdist(presabs_lake[ocean_mixed_sites_geo,], method = "bray", binary = TRUE)
# Stratified lakes and ocean sites
spresbso_dist <- vegdist(presabs_lake[ocean_stratified_sites,], method = "bray", binary = TRUE)
```

## Incidence NMDS

- To constrain dissimilarities we perform Nonmetric Multidimensional
  Scaling (NMDS), which tries to find a stable solution using the
  metaMDS package.

``` r
### Regular
## Jaccard
# All sites 
jpres_NMDS <- metaMDS(jpres_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
round(jpres_NMDS$stress, digits = 2)
# Mixed and stratified lakes
jpresms_NMDS <- metaMDS(jpresms_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
jpresom_NMDS <- metaMDS(jpresom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
jpresso_NMDS <- metaMDS(jpresso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes
jpress_NMDS <- metaMDS(jpress_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)

## Dice-Sørensen
# All sites 
spres_NMDS <- metaMDS(spres_dist, try = 1000, parallel = 4)
# Mixed and stratified lakes
spresms_NMDS <- metaMDS(spresms_dist, try = 1000, parallel = 4)
# Ocean sites and mixed lakes
spresom_NMDS <- metaMDS(spresom_dist, try = 1000, parallel = 4)
# Stratified lakes and ocean sites
spresso_NMDS <- metaMDS(spresso_dist, try = 1000, parallel = 4)


### Environmental
## Jaccard
# All sites 
jprese_NMDS <- metaMDS(jprese_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
jpresems_NMDS <- metaMDS(jpresems_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
jpreseom_NMDS <- metaMDS(jpreseom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
jpreseso_NMDS <- metaMDS(jpreseso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(jpreseso_dist, try = 1000, parallel = 4, previous.best, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
## Dice-Sørensen
# All sites 
sprese_NMDS <- metaMDS(sprese_dist, try = 1000, parallel = 4)
# Mixed and stratified lakes
spresems_NMDS <- metaMDS(spresems_dist, try = 1000, parallel = 4)
# Ocean sites and mixed lakes
spreseom_NMDS <- metaMDS(spreseom_dist, try = 1000, parallel = 4)
# Stratified lakes and ocean sites
spreseso_NMDS <- metaMDS(spreseso_dist, try = 1000, parallel = 4)
```

    ## Warning in metaMDS(spreseso_dist, try = 1000, parallel = 4): stress is (nearly)
    ## zero: you may have insufficient data

``` r
### Geographic
## Jaccard
# All sites 
jpresb_NMDS <- metaMDS(jpresb_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
jpresbms_NMDS <- metaMDS(jpresbms_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
jpresbom_NMDS <- metaMDS(jpresbom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
jpresbso_NMDS <- metaMDS(jpresbso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)

## Dice-Sørensen
# All sites 
spresb_NMDS <- metaMDS(spresb_dist, try = 1000, parallel = 4)
# Mixed and stratified lakes
spresbms_NMDS <- metaMDS(spresbms_dist, try = 1000, parallel = 4)
# Ocean sites and mixed lakes
spresbom_NMDS <- metaMDS(spresbom_dist, try = 1000, parallel = 4)
# Stratified lakes and ocean sites
spresbso_NMDS <- metaMDS(spresbso_dist, try = 1000, parallel = 4)
```

## Determine stratification homogeneity using incidences

- We use betadisper to determine homogeneity of the presence absence
  data based on the stratification category. Is the dispersion of
  incidences similar between the different stratifications.

``` r
## Stratification
# Jaccard
stratification_group_jpres_bd <- betadisper(jpres_dist, stratification_group)
anova(stratification_group_jpres_bd)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.001495 0.0007473  0.1275  0.881
    ## Residuals 19 0.111356 0.0058609

``` r
permutest(stratification_group_jpres_bd, pairwise = TRUE)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     2 0.001495 0.0007473 0.1275    999   0.89
    ## Residuals 19 0.111356 0.0058609                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              Mixed   Ocean Stratified
    ## Mixed              0.83300      0.432
    ## Ocean      0.83147              0.882
    ## Stratified 0.41321 0.85898

``` r
(stratification_group_jpres_bd.HSD <- TukeyHSD(stratification_group_jpres_bd))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff         lwr       upr     p adj
    ## Ocean-Mixed      0.010393222 -0.09464200 0.1154284 0.9658276
    ## Stratified-Mixed 0.019315395 -0.07792832 0.1165591 0.8699882
    ## Stratified-Ocean 0.008922173 -0.09611305 0.1139574 0.9746887

``` r
boxplot(stratification_group_jpres_bd)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Determine%20stratification%20homogeneity%20using%20incidence-1.png)<!-- -->

``` r
# Dice-Sørensen
stratification_group_spres_bd <- betadisper(spres_dist, stratification_group)
anova(stratification_group_spres_bd)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.004075 0.0020373  0.1745 0.8412
    ## Residuals 19 0.221851 0.0116764

``` r
permutest(stratification_group_spres_bd, pairwise = TRUE)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     2 0.004075 0.0020373 0.1745    999  0.829
    ## Residuals 19 0.221851 0.0116764                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              Mixed   Ocean Stratified
    ## Mixed              0.72700      0.378
    ## Ocean      0.71901              0.929
    ## Stratified 0.36777 0.93989

``` r
(stratification_group_spres_bd.HSD <- TukeyHSD(stratification_group_spres_bd))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff        lwr       upr     p adj
    ## Ocean-Mixed      0.024910614 -0.1233440 0.1731653 0.9049381
    ## Stratified-Mixed 0.030232187 -0.1070249 0.1674893 0.8428610
    ## Stratified-Ocean 0.005321574 -0.1429331 0.1535762 0.9954271

``` r
boxplot(stratification_group_spres_bd)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Determine%20stratification%20homogeneity%20using%20incidence-2.png)<!-- -->

## Determine stratification dissimilarity using incidences

- We use adonis to determine if dissimilarities of species incidences by
  stratification categories are significant and how much of the
  variation is explained by incidence dissimilarities.

``` r
### Stratification
## Jaccard
# All sites
jpres_pms <- adonis2(jpres_dist ~ env[surveyed_sites,19], permutations = 999)
jpres_pms
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = jpres_dist ~ env[surveyed_sites, 19], permutations = 999)
    ##                         Df SumOfSqs      R2      F Pr(>F)    
    ## env[surveyed_sites, 19]  2   2.1042 0.25523 3.2556  0.001 ***
    ## Residual                19   6.1401 0.74477                  
    ## Total                   21   8.2443 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratifcation
jpres_pms_pair <- pairwise.adonis(jpres_dist, env[surveyed_sites,19], p.adjust.m = "bonferroni", perm = 999)
jpres_pms_pair
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.3049256 3.903468 0.2454476   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.4999441 1.560433 0.1150725   0.041      0.123

``` r
## Dice-Sørensen
# All sites
spres_pms <- adonis2(spres_dist ~ env[surveyed_sites,19], permutations = 999)
spres_pms
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = spres_dist ~ env[surveyed_sites, 19], permutations = 999)
    ##                         Df SumOfSqs     R2      F Pr(>F)    
    ## env[surveyed_sites, 19]  2   2.5015 0.3572 5.2792  0.001 ***
    ## Residual                19   4.5016 0.6428                  
    ## Total                   21   7.0032 1.0000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratifcation
spres_pms_pair <- pairwise.adonis(spres_dist, env[surveyed_sites,19], p.adjust.m = "bonferroni", perm = 999)
spres_pms_pair
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.6177889 7.163751 0.3384916   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.6459069 6.509883 0.3516977   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.4361964 1.864338 0.1344700   0.068      0.204

## Envfit environmental influence on incidences

- We use envfit to determine significantly correlated environmental
  variables to our NMDS. We will use these results, if significant, in
  our figure.

``` r
### Environmental
## Jaccard
# All sites for figure
jprese_ef <- envfit(jprese_NMDS, env[surveyed_sites_env,c(36)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_env,19])
jprese_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##        NMDS1   NMDS2     r2   Pr(>r)   
    ## [1,] 0.94641 0.32296 0.8332 0.001998 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jprese_efp <- p.adjust.envfit(jprese_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jprese_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##        NMDS1   NMDS2     r2   Pr(>r)   
    ## [1,] 0.94641 0.32296 0.8332 0.001998 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
jpresea_ef <- envfit(jprese_NMDS, env[surveyed_sites_env,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_env,19])
jpresea_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.48785 -0.87293 0.0586 0.823177   
    ## salinity_median     0.99783 -0.06587 0.7388 0.028971 * 
    ## oxygen_median       0.79309  0.60911 0.5361 0.824176   
    ## pH_median           0.94641  0.32296 0.8332 0.002997 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpresea_efp <- p.adjust.envfit(jpresea_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpresea_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.48785 -0.87293 0.0586 1.00000  
    ## salinity_median     0.99783 -0.06587 0.7388 0.11588  
    ## oxygen_median       0.79309  0.60911 0.5361 1.00000  
    ## pH_median           0.94641  0.32296 0.8332 0.01199 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
jpresems_ef <- envfit(jpresems_NMDS, env[mixed_stratified_lakes,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jpresems_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.56523 -0.82494 0.0778 0.789211   
    ## salinity_median     0.99173  0.12834 0.7458 0.036963 * 
    ## oxygen_median       0.94269  0.33366 0.3552 0.895105   
    ## pH_median           0.96334  0.26828 0.8029 0.002997 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpresems_efp <- p.adjust.envfit(jpresems_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpresems_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.56523 -0.82494 0.0778 1.00000  
    ## salinity_median     0.99173  0.12834 0.7458 0.14785  
    ## oxygen_median       0.94269  0.33366 0.3552 1.00000  
    ## pH_median           0.96334  0.26828 0.8029 0.01199 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
jpreseom_ef <- envfit(jpreseom_NMDS, env[ocean_mixed_sites_env,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
jpreseom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.10614  0.99435 0.0705 0.69930  
    ## salinity_median    -0.78731 -0.61655 0.7842 0.01399 *
    ## oxygen_median      -0.96712  0.25434 0.7188 0.04595 *
    ## pH_median          -0.81133  0.58458 0.6420 0.02098 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpreseom_efp <- p.adjust.envfit(jpreseom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpreseom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.10614  0.99435 0.0705 1.00000  
    ## salinity_median    -0.78731 -0.61655 0.7842 0.05594 .
    ## oxygen_median      -0.96712  0.25434 0.7188 0.18382  
    ## pH_median          -0.81133  0.58458 0.6420 0.08392 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
jpreseso_ef <- envfit(jpreseso_NMDS, env[ocean_stratified_sites_env,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
jpreseso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median -0.00006239 -1.00000000 0.0951 0.4895
    ## salinity_median     0.00016684 -1.00000000 0.6308 0.1838
    ## oxygen_median       0.00031725  1.00000000 0.7720 0.2128
    ## pH_median           0.00074615 -1.00000000 0.6784 0.6643
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpreseso_efp <- p.adjust.envfit(jpreseso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpreseso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median -0.00006239 -1.00000000 0.0951 1.0000
    ## salinity_median     0.00016684 -1.00000000 0.6308 0.7353
    ## oxygen_median       0.00031725  1.00000000 0.7720 0.8511
    ## pH_median           0.00074615 -1.00000000 0.6784 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
## Dice-Sørensen
# All sites 
sprese_ef <- envfit(sprese_NMDS, env[surveyed_sites_env,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_env,19])
sprese_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.49140 -0.87093 0.0577 0.824176   
    ## salinity_median     0.99678 -0.08024 0.7389 0.028971 * 
    ## oxygen_median       0.76080  0.64898 0.5414 0.801199   
    ## pH_median           0.95618  0.29278 0.8277 0.001998 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
sprese_efp <- p.adjust.envfit(sprese_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
sprese_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.49140 -0.87093 0.0577 1.000000   
    ## salinity_median     0.99678 -0.08024 0.7389 0.115884   
    ## oxygen_median       0.76080  0.64898 0.5414 1.000000   
    ## pH_median           0.95618  0.29278 0.8277 0.007992 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
spresems_ef <- envfit(spresems_NMDS, env[mixed_stratified_lakes,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
spresems_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.56408 -0.82572 0.0779 0.830170   
    ## salinity_median     0.99171  0.12850 0.7458 0.036963 * 
    ## oxygen_median       0.94199  0.33563 0.3553 0.903097   
    ## pH_median           0.96324  0.26863 0.8030 0.004995 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
spresems_efp <- p.adjust.envfit(spresems_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
spresems_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.56408 -0.82572 0.0779 1.00000  
    ## salinity_median     0.99171  0.12850 0.7458 0.14785  
    ## oxygen_median       0.94199  0.33563 0.3553 1.00000  
    ## pH_median           0.96324  0.26863 0.8030 0.01998 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
spreseom_ef <- envfit(spreseom_NMDS, env[ocean_mixed_sites_env,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
spreseom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.09920  0.99507 0.0762 0.64835  
    ## salinity_median    -0.83313 -0.55308 0.7697 0.01898 *
    ## oxygen_median      -0.98317  0.18269 0.7138 0.04995 *
    ## pH_median          -0.83773  0.54608 0.6380 0.01499 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
spreseom_efp <- p.adjust.envfit(spreseom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
spreseom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.09920  0.99507 0.0762 1.00000  
    ## salinity_median    -0.83313 -0.55308 0.7697 0.07592 .
    ## oxygen_median      -0.98317  0.18269 0.7138 0.19980  
    ## pH_median          -0.83773  0.54608 0.6380 0.05994 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
spreseso_ef <- envfit(spreseso_NMDS, env[ocean_stratified_sites_env,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
spreseso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median  0.00008823  1.00000000 0.1476 0.2887
    ## salinity_median    -0.00109152 -1.00000000 0.4948 0.6613
    ## oxygen_median      -0.00075278  1.00000000 0.7500 0.3127
    ## pH_median          -0.00060409  1.00000000 0.7247 0.2318
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
spreseso_efp <- p.adjust.envfit(spreseso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
spreseso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median  0.00008823  1.00000000 0.1476 1.0000
    ## salinity_median    -0.00109152 -1.00000000 0.4948 1.0000
    ## oxygen_median      -0.00075278  1.00000000 0.7500 1.0000
    ## pH_median          -0.00060409  1.00000000 0.7247 0.9271
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes
jpress_ef <- envfit(jpress_NMDS, env[stratified_lakes,c(2,6,8,10,22:32)], permutations = 1000, na.rm = TRUE)
jpress_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median         -0.35881 -0.93341 0.1964 0.6344
    ## salinity_median            -0.70036  0.71379 0.5991 0.1119
    ## oxygen_median               0.98543 -0.17006 0.5185 0.1748
    ## pH_median                  -0.98992 -0.14163 0.5058 0.1668
    ## volume_m3_w_chemocline     -0.99942 -0.03419 0.1624 0.6793
    ## volume_m3                  -0.44002 -0.89799 0.1263 0.6603
    ## surface_area_m2            -0.35025 -0.93666 0.2598 0.4496
    ## distance_to_ocean_min_m     0.63373 -0.77355 0.2151 0.5275
    ## distance_to_ocean_mean_m    0.01326 -0.99991 0.1726 0.5994
    ## distance_to_ocean_median_m -0.00456 -0.99999 0.2084 0.5495
    ## tidal_lag_time_minutes     -0.95579 -0.29405 0.3022 0.3786
    ## tidal_efficiency            0.35277  0.93571 0.0640 0.8362
    ## perimeter_fromSat          -0.28273 -0.95920 0.3722 0.3097
    ## max_depth                  -0.46293 -0.88640 0.0126 0.9600
    ## logArea                    -0.26513 -0.96421 0.2144 0.5365
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpress_efp <- p.adjust.envfit(jpress_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpress_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median         -0.35881 -0.93341 0.1964      1
    ## salinity_median            -0.70036  0.71379 0.5991      1
    ## oxygen_median               0.98543 -0.17006 0.5185      1
    ## pH_median                  -0.98992 -0.14163 0.5058      1
    ## volume_m3_w_chemocline     -0.99942 -0.03419 0.1624      1
    ## volume_m3                  -0.44002 -0.89799 0.1263      1
    ## surface_area_m2            -0.35025 -0.93666 0.2598      1
    ## distance_to_ocean_min_m     0.63373 -0.77355 0.2151      1
    ## distance_to_ocean_mean_m    0.01326 -0.99991 0.1726      1
    ## distance_to_ocean_median_m -0.00456 -0.99999 0.2084      1
    ## tidal_lag_time_minutes     -0.95579 -0.29405 0.3022      1
    ## tidal_efficiency            0.35277  0.93571 0.0640      1
    ## perimeter_fromSat          -0.28273 -0.95920 0.3722      1
    ## max_depth                  -0.46293 -0.88640 0.0126      1
    ## logArea                    -0.26513 -0.96421 0.2144      1
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Geographic
## Jaccard
# All sites for figure
jpresb_ef <- envfit(jpresb_NMDS, env[surveyed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresb_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.68417 -0.72932 0.1623 0.43656  
    ## volume_m3                   0.64628 -0.76310 0.1557 0.40360  
    ## surface_area_m2             0.43194 -0.90190 0.2799 0.19381  
    ## distance_to_ocean_min_m    -0.98168 -0.19051 0.6466 0.10190  
    ## distance_to_ocean_mean_m   -0.95198 -0.30617 0.5726 0.53946  
    ## distance_to_ocean_median_m -0.94460 -0.32822 0.5670 0.53646  
    ## tidal_lag_time_minutes     -0.99902 -0.04417 0.5908 0.69530  
    ## tidal_efficiency            0.99868  0.05130 0.6114 0.40260  
    ## perimeter_fromSat           0.53003 -0.84798 0.2175 0.38661  
    ## max_depth                  -0.22622 -0.97408 0.1435 0.76923  
    ## logArea                     0.20933 -0.97785 0.3086 0.04795 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpresb_efp <- p.adjust.envfit(jpresb_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpresb_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.68417 -0.72932 0.1623 1.0000
    ## volume_m3                   0.64628 -0.76310 0.1557 1.0000
    ## surface_area_m2             0.43194 -0.90190 0.2799 1.0000
    ## distance_to_ocean_min_m    -0.98168 -0.19051 0.6466 1.0000
    ## distance_to_ocean_mean_m   -0.95198 -0.30617 0.5726 1.0000
    ## distance_to_ocean_median_m -0.94460 -0.32822 0.5670 1.0000
    ## tidal_lag_time_minutes     -0.99902 -0.04417 0.5908 1.0000
    ## tidal_efficiency            0.99868  0.05130 0.6114 1.0000
    ## perimeter_fromSat           0.53003 -0.84798 0.2175 1.0000
    ## max_depth                  -0.22622 -0.97408 0.1435 1.0000
    ## logArea                     0.20933 -0.97785 0.3086 0.5275
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
jpresba_ef <- envfit(jpresb_NMDS, env[surveyed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.68417 -0.72932 0.1623 0.42757  
    ## volume_m3                   0.64628 -0.76310 0.1557 0.40260  
    ## surface_area_m2             0.43194 -0.90190 0.2799 0.17582  
    ## distance_to_ocean_min_m    -0.98168 -0.19051 0.6466 0.11888  
    ## distance_to_ocean_mean_m   -0.95198 -0.30617 0.5726 0.55544  
    ## distance_to_ocean_median_m -0.94460 -0.32822 0.5670 0.54246  
    ## tidal_lag_time_minutes     -0.99902 -0.04417 0.5908 0.69530  
    ## tidal_efficiency            0.99868  0.05130 0.6114 0.41259  
    ## perimeter_fromSat           0.53003 -0.84798 0.2175 0.36364  
    ## max_depth                  -0.22622 -0.97408 0.1435 0.76523  
    ## logArea                     0.20933 -0.97785 0.3086 0.05994 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpresba_efp <- p.adjust.envfit(jpresba_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpresba_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.68417 -0.72932 0.1623 1.0000
    ## volume_m3                   0.64628 -0.76310 0.1557 1.0000
    ## surface_area_m2             0.43194 -0.90190 0.2799 1.0000
    ## distance_to_ocean_min_m    -0.98168 -0.19051 0.6466 1.0000
    ## distance_to_ocean_mean_m   -0.95198 -0.30617 0.5726 1.0000
    ## distance_to_ocean_median_m -0.94460 -0.32822 0.5670 1.0000
    ## tidal_lag_time_minutes     -0.99902 -0.04417 0.5908 1.0000
    ## tidal_efficiency            0.99868  0.05130 0.6114 1.0000
    ## perimeter_fromSat           0.53003 -0.84798 0.2175 1.0000
    ## max_depth                  -0.22622 -0.97408 0.1435 1.0000
    ## logArea                     0.20933 -0.97785 0.3086 0.6593
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
jpresbms_ef <- envfit(jpresbms_NMDS, env[mixed_stratified_lakes_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.35182  0.93607 0.0099 0.9790
    ## volume_m3                  -0.68238 -0.73100 0.1272 0.9031
    ## surface_area_m2            -0.93219  0.36197 0.0832 0.9471
    ## distance_to_ocean_min_m    -0.94998  0.31230 0.5323 0.1648
    ## distance_to_ocean_mean_m   -0.85542  0.51793 0.4286 0.4985
    ## distance_to_ocean_median_m -0.81025  0.58608 0.4367 0.4565
    ## tidal_lag_time_minutes     -0.74304  0.66925 0.4175 0.7073
    ## tidal_efficiency            0.58875 -0.80831 0.5328 0.1848
    ## perimeter_fromSat          -0.25659  0.96652 0.0917 0.7383
    ## max_depth                  -0.35005 -0.93673 0.1947 0.7123
    ## logArea                    -0.88146  0.47225 0.0201 0.9820
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpresbms_efp <- p.adjust.envfit(jpresbms_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpresbms_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.35182  0.93607 0.0099      1
    ## volume_m3                  -0.68238 -0.73100 0.1272      1
    ## surface_area_m2            -0.93219  0.36197 0.0832      1
    ## distance_to_ocean_min_m    -0.94998  0.31230 0.5323      1
    ## distance_to_ocean_mean_m   -0.85542  0.51793 0.4286      1
    ## distance_to_ocean_median_m -0.81025  0.58608 0.4367      1
    ## tidal_lag_time_minutes     -0.74304  0.66925 0.4175      1
    ## tidal_efficiency            0.58875 -0.80831 0.5328      1
    ## perimeter_fromSat          -0.25659  0.96652 0.0917      1
    ## max_depth                  -0.35005 -0.93673 0.1947      1
    ## logArea                    -0.88146  0.47225 0.0201      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
jpresbom_ef <- envfit(jpresbom_NMDS, env[ocean_mixed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline     -0.94020 -0.34063 0.1657 0.54446  
    ## volume_m3                  -0.93986 -0.34156 0.1652 0.54545  
    ## surface_area_m2            -0.73086  0.68253 0.2650 0.50450  
    ## distance_to_ocean_min_m     0.84332 -0.53741 0.4301 0.02997 *
    ## distance_to_ocean_mean_m    0.72901 -0.68450 0.2299 0.32967  
    ## distance_to_ocean_median_m  0.70885 -0.70536 0.2282 0.33866  
    ## tidal_lag_time_minutes      0.85515 -0.51838 0.2043 0.28771  
    ## tidal_efficiency           -0.76016  0.64973 0.2206 0.24575  
    ## perimeter_fromSat          -0.92430  0.38168 0.1845 0.57942  
    ## max_depth                  -0.66524 -0.74663 0.5244 0.03097 *
    ## logArea                    -0.77683  0.62971 0.3273 0.24076  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpresbom_efp <- p.adjust.envfit(jpresbom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpresbom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.94020 -0.34063 0.1657 1.0000
    ## volume_m3                  -0.93986 -0.34156 0.1652 1.0000
    ## surface_area_m2            -0.73086  0.68253 0.2650 1.0000
    ## distance_to_ocean_min_m     0.84332 -0.53741 0.4301 0.3297
    ## distance_to_ocean_mean_m    0.72901 -0.68450 0.2299 1.0000
    ## distance_to_ocean_median_m  0.70885 -0.70536 0.2282 1.0000
    ## tidal_lag_time_minutes      0.85515 -0.51838 0.2043 1.0000
    ## tidal_efficiency           -0.76016  0.64973 0.2206 1.0000
    ## perimeter_fromSat          -0.92430  0.38168 0.1845 1.0000
    ## max_depth                  -0.66524 -0.74663 0.5244 0.3407
    ## logArea                    -0.77683  0.62971 0.3273 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
jpresbso_ef <- envfit(jpresbso_NMDS, env[ocean_stratified_sites,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.94201 -0.33558 0.2252 0.7483
    ## volume_m3                  -0.95253 -0.30443 0.2135 0.7542
    ## surface_area_m2            -0.78575  0.61854 0.3178 0.6873
    ## distance_to_ocean_min_m     0.89003  0.45590 0.6701 0.1818
    ## distance_to_ocean_mean_m    0.89207  0.45190 0.6471 0.5335
    ## distance_to_ocean_median_m  0.87780  0.47902 0.6385 0.5345
    ## tidal_lag_time_minutes      0.98365  0.18009 0.7191 0.9411
    ## tidal_efficiency           -0.97292 -0.23114 0.7970 0.6004
    ## perimeter_fromSat          -0.99079  0.13543 0.2877 0.7143
    ## max_depth                   0.96977  0.24401 0.0579 0.9740
    ## logArea                    -0.54158  0.84065 0.2156 0.5325
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpresbso_efp <- p.adjust.envfit(jpresbso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpresbso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.94201 -0.33558 0.2252      1
    ## volume_m3                  -0.95253 -0.30443 0.2135      1
    ## surface_area_m2            -0.78575  0.61854 0.3178      1
    ## distance_to_ocean_min_m     0.89003  0.45590 0.6701      1
    ## distance_to_ocean_mean_m    0.89207  0.45190 0.6471      1
    ## distance_to_ocean_median_m  0.87780  0.47902 0.6385      1
    ## tidal_lag_time_minutes      0.98365  0.18009 0.7191      1
    ## tidal_efficiency           -0.97292 -0.23114 0.7970      1
    ## perimeter_fromSat          -0.99079  0.13543 0.2877      1
    ## max_depth                   0.96977  0.24401 0.0579      1
    ## logArea                    -0.54158  0.84065 0.2156      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
## Dice-Sørensen
# All sites 
spresb_ef <- envfit(spresb_NMDS, env[surveyed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
spresb_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.68421 -0.72928 0.1623 0.42957  
    ## volume_m3                   0.64632 -0.76306 0.1557 0.40060  
    ## surface_area_m2             0.43194 -0.90190 0.2799 0.18581  
    ## distance_to_ocean_min_m    -0.98169 -0.19049 0.6466 0.09291 .
    ## distance_to_ocean_mean_m   -0.95197 -0.30618 0.5726 0.51748  
    ## distance_to_ocean_median_m -0.94460 -0.32823 0.5670 0.50849  
    ## tidal_lag_time_minutes     -0.99902 -0.04418 0.5908 0.71029  
    ## tidal_efficiency            0.99868  0.05131 0.6114 0.42957  
    ## perimeter_fromSat           0.53005 -0.84797 0.2175 0.36464  
    ## max_depth                  -0.22623 -0.97407 0.1435 0.76424  
    ## logArea                     0.20932 -0.97785 0.3086 0.05794 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
spresb_efp <- p.adjust.envfit(spresb_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
spresb_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.68421 -0.72928 0.1623 1.0000
    ## volume_m3                   0.64632 -0.76306 0.1557 1.0000
    ## surface_area_m2             0.43194 -0.90190 0.2799 1.0000
    ## distance_to_ocean_min_m    -0.98169 -0.19049 0.6466 1.0000
    ## distance_to_ocean_mean_m   -0.95197 -0.30618 0.5726 1.0000
    ## distance_to_ocean_median_m -0.94460 -0.32823 0.5670 1.0000
    ## tidal_lag_time_minutes     -0.99902 -0.04418 0.5908 1.0000
    ## tidal_efficiency            0.99868  0.05131 0.6114 1.0000
    ## perimeter_fromSat           0.53005 -0.84797 0.2175 1.0000
    ## max_depth                  -0.22623 -0.97407 0.1435 1.0000
    ## logArea                     0.20932 -0.97785 0.3086 0.6374
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
spresbms_ef <- envfit(spresbms_NMDS, env[mixed_stratified_lakes_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
spresbms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.35176  0.93609 0.0099 0.9760
    ## volume_m3                  -0.68243 -0.73095 0.1272 0.9101
    ## surface_area_m2            -0.93211  0.36216 0.0832 0.9371
    ## distance_to_ocean_min_m    -0.94997  0.31235 0.5323 0.1568
    ## distance_to_ocean_mean_m   -0.85541  0.51795 0.4286 0.4905
    ## distance_to_ocean_median_m -0.81024  0.58610 0.4367 0.4446
    ## tidal_lag_time_minutes     -0.74305  0.66924 0.4175 0.7253
    ## tidal_efficiency            0.58876 -0.80831 0.5328 0.1828
    ## perimeter_fromSat          -0.25655  0.96653 0.0917 0.7383
    ## max_depth                  -0.35005 -0.93673 0.1947 0.7343
    ## logArea                    -0.88132  0.47253 0.0201 0.9750
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
spresbms_efp <- p.adjust.envfit(spresbms_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
spresbms_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.35176  0.93609 0.0099      1
    ## volume_m3                  -0.68243 -0.73095 0.1272      1
    ## surface_area_m2            -0.93211  0.36216 0.0832      1
    ## distance_to_ocean_min_m    -0.94997  0.31235 0.5323      1
    ## distance_to_ocean_mean_m   -0.85541  0.51795 0.4286      1
    ## distance_to_ocean_median_m -0.81024  0.58610 0.4367      1
    ## tidal_lag_time_minutes     -0.74305  0.66924 0.4175      1
    ## tidal_efficiency            0.58876 -0.80831 0.5328      1
    ## perimeter_fromSat          -0.25655  0.96653 0.0917      1
    ## max_depth                  -0.35005 -0.93673 0.1947      1
    ## logArea                    -0.88132  0.47253 0.0201      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
spresbom_ef <- envfit(spresbom_NMDS, env[ocean_mixed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
spresbom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline     -0.94026 -0.34047 0.1657 0.54246  
    ## volume_m3                  -0.93992 -0.34140 0.1652 0.54845  
    ## surface_area_m2            -0.73074  0.68266 0.2650 0.53247  
    ## distance_to_ocean_min_m     0.84326 -0.53751 0.4301 0.03397 *
    ## distance_to_ocean_mean_m    0.72885 -0.68467 0.2299 0.36963  
    ## distance_to_ocean_median_m  0.70869 -0.70552 0.2282 0.37263  
    ## tidal_lag_time_minutes      0.85498 -0.51866 0.2043 0.32468  
    ## tidal_efficiency           -0.76001  0.64991 0.2206 0.28671  
    ## perimeter_fromSat          -0.92421  0.38187 0.1846 0.58741  
    ## max_depth                  -0.66523 -0.74664 0.5244 0.03397 *
    ## logArea                    -0.77677  0.62979 0.3273 0.21179  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
spresbom_efp <- p.adjust.envfit(spresbom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
spresbom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.94026 -0.34047 0.1657 1.0000
    ## volume_m3                  -0.93992 -0.34140 0.1652 1.0000
    ## surface_area_m2            -0.73074  0.68266 0.2650 1.0000
    ## distance_to_ocean_min_m     0.84326 -0.53751 0.4301 0.3736
    ## distance_to_ocean_mean_m    0.72885 -0.68467 0.2299 1.0000
    ## distance_to_ocean_median_m  0.70869 -0.70552 0.2282 1.0000
    ## tidal_lag_time_minutes      0.85498 -0.51866 0.2043 1.0000
    ## tidal_efficiency           -0.76001  0.64991 0.2206 1.0000
    ## perimeter_fromSat          -0.92421  0.38187 0.1846 1.0000
    ## max_depth                  -0.66523 -0.74664 0.5244 0.3736
    ## logArea                    -0.77677  0.62979 0.3273 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
spresbso_ef <- envfit(spresbso_NMDS, env[ocean_stratified_sites,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
spresbso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.94237 -0.33456 0.2252 0.7233
    ## volume_m3                  -0.95288 -0.30334 0.2135 0.7303
    ## surface_area_m2            -0.78551  0.61884 0.3179 0.6603
    ## distance_to_ocean_min_m     0.88981  0.45632 0.6701 0.2008
    ## distance_to_ocean_mean_m    0.89201  0.45201 0.6471 0.5375
    ## distance_to_ocean_median_m  0.87774  0.47914 0.6385 0.5265
    ## tidal_lag_time_minutes      0.98379  0.17935 0.7191 0.9301
    ## tidal_efficiency           -0.97296 -0.23096 0.7970 0.6064
    ## perimeter_fromSat          -0.99069  0.13612 0.2878 0.7063
    ## max_depth                   0.96929  0.24591 0.0579 0.9720
    ## logArea                    -0.54148  0.84071 0.2156 0.5415
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
spresbso_efp <- p.adjust.envfit(spresbso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
spresbso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.94237 -0.33456 0.2252      1
    ## volume_m3                  -0.95288 -0.30334 0.2135      1
    ## surface_area_m2            -0.78551  0.61884 0.3179      1
    ## distance_to_ocean_min_m     0.88981  0.45632 0.6701      1
    ## distance_to_ocean_mean_m    0.89201  0.45201 0.6471      1
    ## distance_to_ocean_median_m  0.87774  0.47914 0.6385      1
    ## tidal_lag_time_minutes      0.98379  0.17935 0.7191      1
    ## tidal_efficiency           -0.97296 -0.23096 0.7970      1
    ## perimeter_fromSat          -0.99069  0.13612 0.2878      1
    ## max_depth                   0.96929  0.24591 0.0579      1
    ## logArea                    -0.54148  0.84071 0.2156      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Traits
# All sites for figure 
jprest_ef <- envfit(jpres_NMDS, FD_total_env[surveyed_sites,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites,19])
jprest_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.13931  0.99025 0.3396 0.002997 **
    ## Troph            -0.96689 -0.25518 0.7043 0.011988 * 
    ## DepthMin         -0.50885  0.86085 0.0092 0.902098   
    ## DepthMax          0.57886  0.81543 0.2272 0.105894   
    ## TempPrefMin       0.53547 -0.84456 0.3146 0.150849   
    ## TempPrefMax       0.13832 -0.99039 0.0173 0.717283   
    ## DorsalSpinesMean  0.99298 -0.11832 0.5996 0.018981 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.2929 -0.2733
    ## BodyShapeI3f         0.2803  0.1104
    ## BodyShapeI4e        -1.2708 -0.2089
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.1634 -0.3024
    ## OperculumPresentyes  0.2585  0.0672
    ## FeedingPathb         0.1050  0.0528
    ## FeedingPathp        -1.0503 -0.5276
    ## RepGuild11b         -1.1059 -0.4106
    ## RepGuild12g         -1.5644  0.0244
    ## RepGuild13n          0.2229  0.0190
    ## RepGuild22eb        -1.1059 -0.4106
    ## RepGuild23n         -1.5644  0.0244
    ## RepGuild26s          0.2229  0.0190
    ## ParentalCare3p      -1.0912 -0.1387
    ## ParentalCare4n       0.5092  0.0647
    ## WaterPref1s          0.3595  0.0225
    ## WaterPref3a         -1.2224 -0.0764
    ## 
    ## Goodness of fit:
    ##                      r2   Pr(>r)   
    ## BodyShapeI       0.4385 0.037962 * 
    ## DemersPelag      0.0000 1.000000   
    ## OperculumPresent 0.3654 0.091908 . 
    ## FeedingPath      0.1572 0.364635   
    ## RepGuild1        0.3744 0.004995 **
    ## RepGuild2        0.3744 0.004995 **
    ## ParentalCare     0.6427 0.102897   
    ## WaterPref        0.5021 0.046953 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jprest_efp <- p.adjust.envfit(jprest_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jprest_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL       0.13931  0.99025 0.3396 0.04496 *
    ## Troph            -0.96689 -0.25518 0.7043 0.17982  
    ## DepthMin         -0.50885  0.86085 0.0092 1.00000  
    ## DepthMax          0.57886  0.81543 0.2272 1.00000  
    ## TempPrefMin       0.53547 -0.84456 0.3146 1.00000  
    ## TempPrefMax       0.13832 -0.99039 0.0173 1.00000  
    ## DorsalSpinesMean  0.99298 -0.11832 0.5996 0.28472  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.2929 -0.2733
    ## BodyShapeI3f         0.2803  0.1104
    ## BodyShapeI4e        -1.2708 -0.2089
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.1634 -0.3024
    ## OperculumPresentyes  0.2585  0.0672
    ## FeedingPathb         0.1050  0.0528
    ## FeedingPathp        -1.0503 -0.5276
    ## RepGuild11b         -1.1059 -0.4106
    ## RepGuild12g         -1.5644  0.0244
    ## RepGuild13n          0.2229  0.0190
    ## RepGuild22eb        -1.1059 -0.4106
    ## RepGuild23n         -1.5644  0.0244
    ## RepGuild26s          0.2229  0.0190
    ## ParentalCare3p      -1.0912 -0.1387
    ## ParentalCare4n       0.5092  0.0647
    ## WaterPref1s          0.3595  0.0225
    ## WaterPref3a         -1.2224 -0.0764
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.4385 0.56943  
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.3654 1.00000  
    ## FeedingPath      0.1572 1.00000  
    ## RepGuild1        0.3744 0.07493 .
    ## RepGuild2        0.3744 0.07493 .
    ## ParentalCare     0.6427 1.00000  
    ## WaterPref        0.5021 0.70430  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
jpresta_ef <- envfit(jpres_NMDS, FD_total_env[surveyed_sites,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.13931  0.99025 0.3396 0.002997 **
    ## Troph            -0.96689 -0.25518 0.7043 0.014985 * 
    ## DepthMin         -0.50885  0.86085 0.0092 0.917083   
    ## DepthMax          0.57886  0.81543 0.2272 0.097902 . 
    ## TempPrefMin       0.53547 -0.84456 0.3146 0.170829   
    ## TempPrefMax       0.13832 -0.99039 0.0173 0.735265   
    ## DorsalSpinesMean  0.99298 -0.11832 0.5996 0.020979 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.2929 -0.2733
    ## BodyShapeI3f         0.2803  0.1104
    ## BodyShapeI4e        -1.2708 -0.2089
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.1634 -0.3024
    ## OperculumPresentyes  0.2585  0.0672
    ## FeedingPathb         0.1050  0.0528
    ## FeedingPathp        -1.0503 -0.5276
    ## RepGuild11b         -1.1059 -0.4106
    ## RepGuild12g         -1.5644  0.0244
    ## RepGuild13n          0.2229  0.0190
    ## RepGuild22eb        -1.1059 -0.4106
    ## RepGuild23n         -1.5644  0.0244
    ## RepGuild26s          0.2229  0.0190
    ## ParentalCare3p      -1.0912 -0.1387
    ## ParentalCare4n       0.5092  0.0647
    ## WaterPref1s          0.3595  0.0225
    ## WaterPref3a         -1.2224 -0.0764
    ## 
    ## Goodness of fit:
    ##                      r2   Pr(>r)   
    ## BodyShapeI       0.4385 0.037962 * 
    ## DemersPelag      0.0000 1.000000   
    ## OperculumPresent 0.3654 0.081918 . 
    ## FeedingPath      0.1572 0.352647   
    ## RepGuild1        0.3744 0.003996 **
    ## RepGuild2        0.3744 0.003996 **
    ## ParentalCare     0.6427 0.127872   
    ## WaterPref        0.5021 0.045954 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jpresta_efp <- p.adjust.envfit(jpresta_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jpresta_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL       0.13931  0.99025 0.3396 0.04496 *
    ## Troph            -0.96689 -0.25518 0.7043 0.22478  
    ## DepthMin         -0.50885  0.86085 0.0092 1.00000  
    ## DepthMax          0.57886  0.81543 0.2272 1.00000  
    ## TempPrefMin       0.53547 -0.84456 0.3146 1.00000  
    ## TempPrefMax       0.13832 -0.99039 0.0173 1.00000  
    ## DorsalSpinesMean  0.99298 -0.11832 0.5996 0.31469  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.2929 -0.2733
    ## BodyShapeI3f         0.2803  0.1104
    ## BodyShapeI4e        -1.2708 -0.2089
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.1634 -0.3024
    ## OperculumPresentyes  0.2585  0.0672
    ## FeedingPathb         0.1050  0.0528
    ## FeedingPathp        -1.0503 -0.5276
    ## RepGuild11b         -1.1059 -0.4106
    ## RepGuild12g         -1.5644  0.0244
    ## RepGuild13n          0.2229  0.0190
    ## RepGuild22eb        -1.1059 -0.4106
    ## RepGuild23n         -1.5644  0.0244
    ## RepGuild26s          0.2229  0.0190
    ## ParentalCare3p      -1.0912 -0.1387
    ## ParentalCare4n       0.5092  0.0647
    ## WaterPref1s          0.3595  0.0225
    ## WaterPref3a         -1.2224 -0.0764
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.4385 0.56943  
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.3654 1.00000  
    ## FeedingPath      0.1572 1.00000  
    ## RepGuild1        0.3744 0.05994 .
    ## RepGuild2        0.3744 0.05994 .
    ## ParentalCare     0.6427 1.00000  
    ## WaterPref        0.5021 0.68931  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
jprestms_ef <- envfit(jpresms_NMDS, FD_total_env[mixed_stratified_lakes,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)    
    ## MaxLengthTL       0.24333  0.96994 0.4596 0.022977 *  
    ## Troph            -0.99983 -0.01840 0.5841 0.074925 .  
    ## DepthMin         -0.11859  0.99294 0.1175 0.297702    
    ## DepthMax          0.41149  0.91141 0.4799 0.025974 *  
    ## TempPrefMin       0.14160 -0.98992 0.6620 0.000999 ***
    ## TempPrefMax       0.01259 -0.99992 0.4448 0.006993 ** 
    ## DorsalSpinesMean  0.98429 -0.17657 0.6164 0.020979 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.8457 -0.2645
    ## BodyShapeI3f         0.4293  0.0093
    ## BodyShapeI4e        -0.9691  0.0405
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -0.9042 -0.0953
    ## OperculumPresentyes  0.3014  0.0318
    ## FeedingPathb         0.1221  0.0497
    ## FeedingPathp        -0.8545 -0.3481
    ## RepGuild11b         -0.8457 -0.2645
    ## RepGuild12g         -1.1597  0.1376
    ## RepGuild13n          0.2435 -0.0008
    ## RepGuild22eb        -0.8457 -0.2645
    ## RepGuild23n         -1.1597  0.1376
    ## RepGuild26s          0.2435 -0.0008
    ## ParentalCare3p      -0.7717 -0.0099
    ## ParentalCare4n       0.6002  0.0077
    ## WaterPref1s          0.3939 -0.0081
    ## WaterPref3a         -0.8665  0.0179
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.6398 0.02398 *
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.4290 0.08891 .
    ## FeedingPath      0.1893 0.29071  
    ## RepGuild1        0.4168 0.01698 *
    ## RepGuild2        0.4168 0.01698 *
    ## ParentalCare     0.7213 0.12188  
    ## WaterPref        0.5315 0.07193 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jprestms_efp <- p.adjust.envfit(jprestms_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jprestms_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL       0.24333  0.96994 0.4596 0.34466  
    ## Troph            -0.99983 -0.01840 0.5841 1.00000  
    ## DepthMin         -0.11859  0.99294 0.1175 1.00000  
    ## DepthMax          0.41149  0.91141 0.4799 0.38961  
    ## TempPrefMin       0.14160 -0.98992 0.6620 0.01499 *
    ## TempPrefMax       0.01259 -0.99992 0.4448 0.10490  
    ## DorsalSpinesMean  0.98429 -0.17657 0.6164 0.31469  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.8457 -0.2645
    ## BodyShapeI3f         0.4293  0.0093
    ## BodyShapeI4e        -0.9691  0.0405
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -0.9042 -0.0953
    ## OperculumPresentyes  0.3014  0.0318
    ## FeedingPathb         0.1221  0.0497
    ## FeedingPathp        -0.8545 -0.3481
    ## RepGuild11b         -0.8457 -0.2645
    ## RepGuild12g         -1.1597  0.1376
    ## RepGuild13n          0.2435 -0.0008
    ## RepGuild22eb        -0.8457 -0.2645
    ## RepGuild23n         -1.1597  0.1376
    ## RepGuild26s          0.2435 -0.0008
    ## ParentalCare3p      -0.7717 -0.0099
    ## ParentalCare4n       0.6002  0.0077
    ## WaterPref1s          0.3939 -0.0081
    ## WaterPref3a         -0.8665  0.0179
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.6398 0.3596
    ## DemersPelag      0.0000 1.0000
    ## OperculumPresent 0.4290 1.0000
    ## FeedingPath      0.1893 1.0000
    ## RepGuild1        0.4168 0.2547
    ## RepGuild2        0.4168 0.2547
    ## ParentalCare     0.7213 1.0000
    ## WaterPref        0.5315 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
jprestom_ef <- envfit(jpresom_NMDS, FD_total_env[ocean_mixed_sites,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
jprestom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.99904 -0.04388 0.4748 0.018981 * 
    ## Troph             0.66398  0.74775 0.6802 0.001998 **
    ## DepthMin         -0.37268 -0.92796 0.2011 0.184815   
    ## DepthMax          0.02827  0.99960 0.0194 0.907093   
    ## TempPrefMin      -0.54984  0.83527 0.2882 0.120879   
    ## TempPrefMax      -0.08619  0.99628 0.0085 0.948052   
    ## DorsalSpinesMean -0.77597  0.63077 0.3617 0.078921 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.7194 -0.0735
    ## BodyShapeI3f         0.1199  0.0122
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentyes  0.0000  0.0000
    ## FeedingPathb         0.0000  0.0000
    ## RepGuild13n          0.0000  0.0000
    ## RepGuild26s          0.0000  0.0000
    ## ParentalCare4n       0.0000  0.0000
    ## WaterPref1s          0.0000  0.0000
    ## 
    ## Goodness of fit:
    ##                     r2 Pr(>r)
    ## BodyShapeI       0.117 0.6823
    ## DemersPelag      0.000 1.0000
    ## OperculumPresent 0.000 1.0000
    ## FeedingPath      0.000 1.0000
    ## RepGuild1        0.000 1.0000
    ## RepGuild2        0.000 1.0000
    ## ParentalCare     0.000 1.0000
    ## WaterPref        0.000 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jprestom_efp <- p.adjust.envfit(jprestom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jprestom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL       0.99904 -0.04388 0.4748 0.28472  
    ## Troph             0.66398  0.74775 0.6802 0.02997 *
    ## DepthMin         -0.37268 -0.92796 0.2011 1.00000  
    ## DepthMax          0.02827  0.99960 0.0194 1.00000  
    ## TempPrefMin      -0.54984  0.83527 0.2882 1.00000  
    ## TempPrefMax      -0.08619  0.99628 0.0085 1.00000  
    ## DorsalSpinesMean -0.77597  0.63077 0.3617 1.00000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.7194 -0.0735
    ## BodyShapeI3f         0.1199  0.0122
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentyes  0.0000  0.0000
    ## FeedingPathb         0.0000  0.0000
    ## RepGuild13n          0.0000  0.0000
    ## RepGuild26s          0.0000  0.0000
    ## ParentalCare4n       0.0000  0.0000
    ## WaterPref1s          0.0000  0.0000
    ## 
    ## Goodness of fit:
    ##                     r2 Pr(>r)
    ## BodyShapeI       0.117      1
    ## DemersPelag      0.000      1
    ## OperculumPresent 0.000      1
    ## FeedingPath      0.000      1
    ## RepGuild1        0.000      1
    ## RepGuild2        0.000      1
    ## ParentalCare     0.000      1
    ## WaterPref        0.000      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
jprestso_ef <- envfit(jpresso_NMDS, FD_total_env[ocean_stratified_sites,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL      -0.10571 -0.99440 0.4038 0.001998 **
    ## Troph             0.78314  0.62184 0.8719 0.007992 **
    ## DepthMin          0.19809 -0.98018 0.0304 0.696304   
    ## DepthMax         -0.48520 -0.87440 0.2052 0.074925 . 
    ## TempPrefMin      -0.38266  0.92389 0.3950 0.095904 . 
    ## TempPrefMax       0.04025  0.99919 0.0448 0.546454   
    ## DorsalSpinesMean -0.99785 -0.06555 0.6362 0.042957 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.4149  0.1729
    ## BodyShapeI3f        -0.2908 -0.1445
    ## BodyShapeI4e         0.8201  0.1232
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno   0.7469  0.2218
    ## OperculumPresentyes -0.2988 -0.0887
    ## FeedingPathb        -0.1148 -0.0594
    ## FeedingPathp         0.6890  0.3566
    ## RepGuild11b          0.7281  0.2733
    ## RepGuild12g          1.0431 -0.0618
    ## RepGuild13n         -0.2559 -0.0136
    ## RepGuild22eb         0.7281  0.2733
    ## RepGuild23n          1.0431 -0.0618
    ## RepGuild26s         -0.2559 -0.0136
    ## ParentalCare3p       0.6983  0.0440
    ## ParentalCare4n      -0.6983 -0.0440
    ## WaterPref1s         -0.4397 -0.0019
    ## WaterPref3a          0.7914  0.0035
    ## 
    ## Goodness of fit:
    ##                      r2   Pr(>r)   
    ## BodyShapeI       0.4042 0.068931 . 
    ## DemersPelag      0.0000 1.000000   
    ## OperculumPresent 0.3356 0.113886   
    ## FeedingPath      0.1386 0.377622   
    ## RepGuild1        0.3465 0.001998 **
    ## RepGuild2        0.3465 0.001998 **
    ## ParentalCare     0.6765 0.130869   
    ## WaterPref        0.4808 0.055944 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
jprestso_efp <- p.adjust.envfit(jprestso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
jprestso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL      -0.10571 -0.99440 0.4038 0.02997 *
    ## Troph             0.78314  0.62184 0.8719 0.11988  
    ## DepthMin          0.19809 -0.98018 0.0304 1.00000  
    ## DepthMax         -0.48520 -0.87440 0.2052 1.00000  
    ## TempPrefMin      -0.38266  0.92389 0.3950 1.00000  
    ## TempPrefMax       0.04025  0.99919 0.0448 1.00000  
    ## DorsalSpinesMean -0.99785 -0.06555 0.6362 0.64436  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.4149  0.1729
    ## BodyShapeI3f        -0.2908 -0.1445
    ## BodyShapeI4e         0.8201  0.1232
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno   0.7469  0.2218
    ## OperculumPresentyes -0.2988 -0.0887
    ## FeedingPathb        -0.1148 -0.0594
    ## FeedingPathp         0.6890  0.3566
    ## RepGuild11b          0.7281  0.2733
    ## RepGuild12g          1.0431 -0.0618
    ## RepGuild13n         -0.2559 -0.0136
    ## RepGuild22eb         0.7281  0.2733
    ## RepGuild23n          1.0431 -0.0618
    ## RepGuild26s         -0.2559 -0.0136
    ## ParentalCare3p       0.6983  0.0440
    ## ParentalCare4n      -0.6983 -0.0440
    ## WaterPref1s         -0.4397 -0.0019
    ## WaterPref3a          0.7914  0.0035
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.4042 1.00000  
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.3356 1.00000  
    ## FeedingPath      0.1386 1.00000  
    ## RepGuild1        0.3465 0.02997 *
    ## RepGuild2        0.3465 0.02997 *
    ## ParentalCare     0.6765 1.00000  
    ## WaterPref        0.4808 0.83916  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Traits
# All sites 
sprest_ef <- envfit(spres_NMDS, FD_total_env[surveyed_sites,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites,19])
sprest_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.13491  0.99086 0.3491 0.004995 **
    ## Troph            -0.95246 -0.30465 0.7103 0.005994 **
    ## DepthMin         -0.44727  0.89440 0.0106 0.910090   
    ## DepthMax          0.56463  0.82535 0.2272 0.115884   
    ## TempPrefMin       0.47615 -0.87937 0.3424 0.114885   
    ## TempPrefMax       0.11154 -0.99376 0.0256 0.623377   
    ## DorsalSpinesMean  0.99970 -0.02435 0.6034 0.031968 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.3319 -0.3129
    ## BodyShapeI3f         0.3161  0.1229
    ## BodyShapeI4e        -1.4344 -0.2261
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.3108 -0.3346
    ## OperculumPresentyes  0.2913  0.0743
    ## FeedingPathb         0.1188  0.0594
    ## FeedingPathp        -1.1878 -0.5941
    ## RepGuild11b         -1.2429 -0.4654
    ## RepGuild12g         -1.7642  0.0453
    ## RepGuild13n          0.2511  0.0197
    ## RepGuild22eb        -1.2429 -0.4654
    ## RepGuild23n         -1.7642  0.0453
    ## RepGuild26s          0.2511  0.0197
    ## ParentalCare3p      -1.2300 -0.1499
    ## ParentalCare4n       0.5740  0.0699
    ## WaterPref1s          0.4045  0.0243
    ## WaterPref3a         -1.3752 -0.0826
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)   
    ## BodyShapeI       0.4417 0.04296 * 
    ## DemersPelag      0.0000 1.00000   
    ## OperculumPresent 0.3665 0.10190   
    ## FeedingPath      0.1589 0.37063   
    ## RepGuild1        0.3767 0.00999 **
    ## RepGuild2        0.3767 0.00999 **
    ## ParentalCare     0.6457 0.12587   
    ## WaterPref        0.5030 0.05794 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
sprest_efp <- p.adjust.envfit(sprest_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
sprest_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL       0.13491  0.99086 0.3491 0.07493 .
    ## Troph            -0.95246 -0.30465 0.7103 0.08991 .
    ## DepthMin         -0.44727  0.89440 0.0106 1.00000  
    ## DepthMax          0.56463  0.82535 0.2272 1.00000  
    ## TempPrefMin       0.47615 -0.87937 0.3424 1.00000  
    ## TempPrefMax       0.11154 -0.99376 0.0256 1.00000  
    ## DorsalSpinesMean  0.99970 -0.02435 0.6034 0.47952  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.3319 -0.3129
    ## BodyShapeI3f         0.3161  0.1229
    ## BodyShapeI4e        -1.4344 -0.2261
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.3108 -0.3346
    ## OperculumPresentyes  0.2913  0.0743
    ## FeedingPathb         0.1188  0.0594
    ## FeedingPathp        -1.1878 -0.5941
    ## RepGuild11b         -1.2429 -0.4654
    ## RepGuild12g         -1.7642  0.0453
    ## RepGuild13n          0.2511  0.0197
    ## RepGuild22eb        -1.2429 -0.4654
    ## RepGuild23n         -1.7642  0.0453
    ## RepGuild26s          0.2511  0.0197
    ## ParentalCare3p      -1.2300 -0.1499
    ## ParentalCare4n       0.5740  0.0699
    ## WaterPref1s          0.4045  0.0243
    ## WaterPref3a         -1.3752 -0.0826
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.4417 0.6444
    ## DemersPelag      0.0000 1.0000
    ## OperculumPresent 0.3665 1.0000
    ## FeedingPath      0.1589 1.0000
    ## RepGuild1        0.3767 0.1499
    ## RepGuild2        0.3767 0.1499
    ## ParentalCare     0.6457 1.0000
    ## WaterPref        0.5030 0.8691
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
sprestms_ef <- envfit(spresms_NMDS, FD_total_env[mixed_stratified_lakes,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
sprestms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)    
    ## MaxLengthTL       0.24339  0.96993 0.4593 0.028971 *  
    ## Troph            -0.99980 -0.01987 0.5842 0.068931 .  
    ## DepthMin         -0.11843  0.99296 0.1177 0.299700    
    ## DepthMax          0.41165  0.91134 0.4797 0.030969 *  
    ## TempPrefMin       0.14158 -0.98993 0.6621 0.000999 ***
    ## TempPrefMax       0.01259 -0.99992 0.4447 0.005994 ** 
    ## DorsalSpinesMean  0.98430 -0.17648 0.6164 0.032967 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -1.0368 -0.3243
    ## BodyShapeI3f         0.5262  0.0114
    ## BodyShapeI4e        -1.1879  0.0498
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.1085 -0.1167
    ## OperculumPresentyes  0.3695  0.0389
    ## FeedingPathb         0.1497  0.0609
    ## FeedingPathp        -1.0476 -0.4266
    ## RepGuild11b         -1.0368 -0.3243
    ## RepGuild12g         -1.4216  0.1690
    ## RepGuild13n          0.2985 -0.0010
    ## RepGuild22eb        -1.0368 -0.3243
    ## RepGuild23n         -1.4216  0.1690
    ## RepGuild26s          0.2985 -0.0010
    ## ParentalCare3p      -0.9461 -0.0120
    ## ParentalCare4n       0.7358  0.0094
    ## WaterPref1s          0.4828 -0.0100
    ## WaterPref3a         -1.0622  0.0220
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.6398 0.01798 *
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.4290 0.08791 .
    ## FeedingPath      0.1893 0.28671  
    ## RepGuild1        0.4168 0.02098 *
    ## RepGuild2        0.4168 0.02098 *
    ## ParentalCare     0.7213 0.12288  
    ## WaterPref        0.5315 0.06993 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
sprestms_efp <- p.adjust.envfit(sprestms_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
sprestms_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL       0.24339  0.96993 0.4593 0.43457  
    ## Troph            -0.99980 -0.01987 0.5842 1.00000  
    ## DepthMin         -0.11843  0.99296 0.1177 1.00000  
    ## DepthMax          0.41165  0.91134 0.4797 0.46454  
    ## TempPrefMin       0.14158 -0.98993 0.6621 0.01499 *
    ## TempPrefMax       0.01259 -0.99992 0.4447 0.08991 .
    ## DorsalSpinesMean  0.98430 -0.17648 0.6164 0.49451  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -1.0368 -0.3243
    ## BodyShapeI3f         0.5262  0.0114
    ## BodyShapeI4e        -1.1879  0.0498
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.1085 -0.1167
    ## OperculumPresentyes  0.3695  0.0389
    ## FeedingPathb         0.1497  0.0609
    ## FeedingPathp        -1.0476 -0.4266
    ## RepGuild11b         -1.0368 -0.3243
    ## RepGuild12g         -1.4216  0.1690
    ## RepGuild13n          0.2985 -0.0010
    ## RepGuild22eb        -1.0368 -0.3243
    ## RepGuild23n         -1.4216  0.1690
    ## RepGuild26s          0.2985 -0.0010
    ## ParentalCare3p      -0.9461 -0.0120
    ## ParentalCare4n       0.7358  0.0094
    ## WaterPref1s          0.4828 -0.0100
    ## WaterPref3a         -1.0622  0.0220
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.6398 0.2697
    ## DemersPelag      0.0000 1.0000
    ## OperculumPresent 0.4290 1.0000
    ## FeedingPath      0.1893 1.0000
    ## RepGuild1        0.4168 0.3147
    ## RepGuild2        0.4168 0.3147
    ## ParentalCare     0.7213 1.0000
    ## WaterPref        0.5315 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
sprestom_ef <- envfit(spresom_NMDS, FD_total_env[ocean_mixed_sites,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
sprestom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.99903 -0.04411 0.4748 0.020979 * 
    ## Troph             0.66405  0.74769 0.6801 0.003996 **
    ## DepthMin         -0.37271 -0.92795 0.2010 0.176823   
    ## DepthMax          0.02830  0.99960 0.0194 0.917083   
    ## TempPrefMin      -0.54974  0.83533 0.2882 0.093906 . 
    ## TempPrefMax      -0.08616  0.99628 0.0085 0.925075   
    ## DorsalSpinesMean -0.77587  0.63089 0.3617 0.070929 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.5434 -0.0554
    ## BodyShapeI3f         0.0906  0.0092
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentyes  0.0000  0.0000
    ## FeedingPathb         0.0000  0.0000
    ## RepGuild13n          0.0000  0.0000
    ## RepGuild26s          0.0000  0.0000
    ## ParentalCare4n       0.0000  0.0000
    ## WaterPref1s          0.0000  0.0000
    ## 
    ## Goodness of fit:
    ##                     r2 Pr(>r)
    ## BodyShapeI       0.117 0.6683
    ## DemersPelag      0.000 1.0000
    ## OperculumPresent 0.000 1.0000
    ## FeedingPath      0.000 1.0000
    ## RepGuild1        0.000 1.0000
    ## RepGuild2        0.000 1.0000
    ## ParentalCare     0.000 1.0000
    ## WaterPref        0.000 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
sprestom_efp <- p.adjust.envfit(sprestom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
sprestom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL       0.99903 -0.04411 0.4748 0.31469  
    ## Troph             0.66405  0.74769 0.6801 0.05994 .
    ## DepthMin         -0.37271 -0.92795 0.2010 1.00000  
    ## DepthMax          0.02830  0.99960 0.0194 1.00000  
    ## TempPrefMin      -0.54974  0.83533 0.2882 1.00000  
    ## TempPrefMax      -0.08616  0.99628 0.0085 1.00000  
    ## DorsalSpinesMean -0.77587  0.63089 0.3617 1.00000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.5434 -0.0554
    ## BodyShapeI3f         0.0906  0.0092
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentyes  0.0000  0.0000
    ## FeedingPathb         0.0000  0.0000
    ## RepGuild13n          0.0000  0.0000
    ## RepGuild26s          0.0000  0.0000
    ## ParentalCare4n       0.0000  0.0000
    ## WaterPref1s          0.0000  0.0000
    ## 
    ## Goodness of fit:
    ##                     r2 Pr(>r)
    ## BodyShapeI       0.117      1
    ## DemersPelag      0.000      1
    ## OperculumPresent 0.000      1
    ## FeedingPath      0.000      1
    ## RepGuild1        0.000      1
    ## RepGuild2        0.000      1
    ## ParentalCare     0.000      1
    ## WaterPref        0.000      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
sprestso_ef <- envfit(spresso_NMDS, FD_total_env[ocean_stratified_sites,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
sprestso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL      -0.10568 -0.99440 0.4038 0.001998 **
    ## Troph             0.78323  0.62174 0.8718 0.012987 * 
    ## DepthMin          0.19846 -0.98011 0.0303 0.697303   
    ## DepthMax         -0.48496 -0.87454 0.2052 0.106893   
    ## TempPrefMin      -0.38273  0.92386 0.3950 0.086913 . 
    ## TempPrefMax       0.04024  0.99919 0.0448 0.551449   
    ## DorsalSpinesMean -0.99781 -0.06608 0.6362 0.053946 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.5874  0.2450
    ## BodyShapeI3f        -0.4117 -0.2048
    ## BodyShapeI4e         1.1611  0.1746
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno   1.0575  0.3141
    ## OperculumPresentyes -0.4230 -0.1257
    ## FeedingPathb        -0.1626 -0.0841
    ## FeedingPathp         0.9754  0.5048
    ## RepGuild11b          1.0308  0.3869
    ## RepGuild12g          1.4770 -0.0872
    ## RepGuild13n         -0.3622 -0.0193
    ## RepGuild22eb         1.0308  0.3869
    ## RepGuild23n          1.4770 -0.0872
    ## RepGuild26s         -0.3622 -0.0193
    ## ParentalCare3p       0.9887  0.0624
    ## ParentalCare4n      -0.9887 -0.0624
    ## WaterPref1s         -0.6225 -0.0028
    ## WaterPref3a          1.1205  0.0051
    ## 
    ## Goodness of fit:
    ##                      r2   Pr(>r)   
    ## BodyShapeI       0.4042 0.059940 . 
    ## DemersPelag      0.0000 1.000000   
    ## OperculumPresent 0.3355 0.115884   
    ## FeedingPath      0.1386 0.353646   
    ## RepGuild1        0.3465 0.004995 **
    ## RepGuild2        0.3465 0.004995 **
    ## ParentalCare     0.6765 0.133866   
    ## WaterPref        0.4808 0.041958 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
sprestso_efp <- p.adjust.envfit(sprestso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
sprestso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL      -0.10568 -0.99440 0.4038 0.02997 *
    ## Troph             0.78323  0.62174 0.8718 0.19481  
    ## DepthMin          0.19846 -0.98011 0.0303 1.00000  
    ## DepthMax         -0.48496 -0.87454 0.2052 1.00000  
    ## TempPrefMin      -0.38273  0.92386 0.3950 1.00000  
    ## TempPrefMax       0.04024  0.99919 0.0448 1.00000  
    ## DorsalSpinesMean -0.99781 -0.06608 0.6362 0.80919  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.5874  0.2450
    ## BodyShapeI3f        -0.4117 -0.2048
    ## BodyShapeI4e         1.1611  0.1746
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno   1.0575  0.3141
    ## OperculumPresentyes -0.4230 -0.1257
    ## FeedingPathb        -0.1626 -0.0841
    ## FeedingPathp         0.9754  0.5048
    ## RepGuild11b          1.0308  0.3869
    ## RepGuild12g          1.4770 -0.0872
    ## RepGuild13n         -0.3622 -0.0193
    ## RepGuild22eb         1.0308  0.3869
    ## RepGuild23n          1.4770 -0.0872
    ## RepGuild26s         -0.3622 -0.0193
    ## ParentalCare3p       0.9887  0.0624
    ## ParentalCare4n      -0.9887 -0.0624
    ## WaterPref1s         -0.6225 -0.0028
    ## WaterPref3a          1.1205  0.0051
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.4042 0.89910  
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.3355 1.00000  
    ## FeedingPath      0.1386 1.00000  
    ## RepGuild1        0.3465 0.07493 .
    ## RepGuild2        0.3465 0.07493 .
    ## ParentalCare     0.6765 1.00000  
    ## WaterPref        0.4808 0.62937  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

## Incidence mantel tests

- We used mantel tests to determine significance between the incidence
  and environmental distance matrices.

``` r
### Environmental
## Jaccard
# All sites
enve_dist_t <- dist(scaled_env[surveyed_sites_env,c(1)], method = "euclidean")
jpresea_mant_t <- mantel(jprese_dist, enve_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
jpresea_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jprese_dist, ydis = enve_dist_t, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.241 
    ##       Significance: 0.271 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.306 0.332 0.356 0.392 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_s <- dist(scaled_env[surveyed_sites_env,c(2)], method = "euclidean")
jpresea_mant_s <- mantel(jprese_dist, enve_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
jpresea_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jprese_dist, ydis = enve_dist_s, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6697 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.565 0.600 0.625 0.639 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_o <- dist(scaled_env[surveyed_sites_env,c(3)], method = "euclidean")
jpresea_mant_o <- mantel(jprese_dist, enve_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
jpresea_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jprese_dist, ydis = enve_dist_o, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r:  0.47 
    ##       Significance: 0.467 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.589 0.632 0.666 0.700 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_p <- dist(scaled_env[surveyed_sites_env,c(4)], method = "euclidean")
jpresea_mant_p <- mantel(jprese_dist, enve_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
jpresea_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jprese_dist, ydis = enve_dist_p, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7163 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.540 0.574 0.601 0.643 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresea_mant_pv <- rbind(jpresea_mant_t$signif, jpresea_mant_s$signif, jpresea_mant_o$signif, jpresea_mant_p$signif)
jpresea_mant_pv <- jpresea_mant_pv[,1]
jpresea_mant_pv <- p.adjust(jpresea_mant_pv, method = "bonferroni")
jpresea_mant_pv
```

    ## [1] 1.000 0.008 1.000 0.008

``` r
# Mixed and stratified lakes
envems_dist_t <- dist(scaled_env[mixed_stratified_lakes,c(1)], method = "euclidean")
jpresems_mant_t <- mantel(jpresems_dist, envems_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jpresems_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresems_dist, ydis = envems_dist_t, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2368 
    ##       Significance: 0.33 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.337 0.376 0.404 0.429 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_s <- dist(scaled_env[mixed_stratified_lakes,c(2)], method = "euclidean")
jpresems_mant_s <- mantel(jpresems_dist, envems_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jpresems_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresems_dist, ydis = envems_dist_s, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6676 
    ##       Significance: 0.013 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.569 0.603 0.631 0.677 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_o <- dist(scaled_env[mixed_stratified_lakes,c(3)], method = "euclidean")
jpresems_mant_o <- mantel(jpresems_dist, envems_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jpresems_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresems_dist, ydis = envems_dist_o, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2979 
    ##       Significance: 0.607 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.477 0.518 0.563 0.628 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_p <- dist(scaled_env[mixed_stratified_lakes,c(4)], method = "euclidean")
jpresems_mant_p <- mantel(jpresems_dist, envems_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jpresems_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresems_dist, ydis = envems_dist_p, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6357 
    ##       Significance: 0.003 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.431 0.470 0.508 0.544 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresems_mant_pv <- rbind(jpresems_mant_t$signif, jpresems_mant_s$signif, jpresems_mant_o$signif, jpresems_mant_p$signif)
jpresems_mant_pv <- jpresems_mant_pv[,1]
jpresems_mant_pv <- p.adjust(jpresems_mant_pv, method = "bonferroni")
jpresems_mant_pv
```

    ## [1] 1.000 0.052 1.000 0.012

``` r
# Ocean sites and mixed lakes
enveom_dist_t <- dist(scaled_env[ocean_mixed_sites_env,c(1)], method = "euclidean")
jpreseom_mant_t <- mantel(jpreseom_dist, enveom_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
jpreseom_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpreseom_dist, ydis = enveom_dist_t, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1214 
    ##       Significance: 0.754 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.124 0.170 0.206 0.378 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_dist_s <- dist(scaled_env[ocean_mixed_sites_env,c(2)], method = "euclidean")
jpreseom_mant_s <- mantel(jpreseom_dist, enveom_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
jpreseom_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpreseom_dist, ydis = enveom_dist_s, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4258 
    ##       Significance: 0.056 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.341 0.429 0.493 0.565 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_dist_o <- dist(scaled_env[ocean_mixed_sites_env,c(3)], method = "euclidean")
jpreseom_mant_o <- mantel(jpreseom_dist, enveom_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
jpreseom_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpreseom_dist, ydis = enveom_dist_o, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5319 
    ##       Significance: 0.032 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.326 0.456 0.552 0.659 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_dist_p <- dist(scaled_env[ocean_mixed_sites_env,c(4)], method = "euclidean")
jpreseom_mant_p <- mantel(jpreseom_dist, enveom_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
jpreseom_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpreseom_dist, ydis = enveom_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5793 
    ##       Significance: 0.016 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.316 0.403 0.476 0.609 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpreseom_mant_pv <- rbind(jpreseom_mant_t$signif, jpreseom_mant_s$signif, jpreseom_mant_o$signif, jpreseom_mant_p$signif)
jpreseom_mant_pv <- jpreseom_mant_pv[,1]
jpreseom_mant_pv <- p.adjust(jpreseom_mant_pv, method = "bonferroni")
jpreseom_mant_pv
```

    ## [1] 1.000 0.224 0.128 0.064

``` r
# Stratified lakes and ocean sites
enveso_dist_t <- dist(scaled_env[ocean_stratified_sites_env,c(1)], method = "euclidean")
jpreseso_mant_t <- mantel(jpreseso_dist, enveso_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
jpreseso_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpreseso_dist, ydis = enveso_dist_t, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.008188 
    ##       Significance: 0.209 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0221 0.0520 0.0676 0.0781 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveso_dist_s <- dist(scaled_env[ocean_stratified_sites_env,c(2)], method = "euclidean")
jpreseso_mant_s <- mantel(jpreseso_dist, enveso_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
jpreseso_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpreseso_dist, ydis = enveso_dist_s, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5528 
    ##       Significance: 0.015 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.453 0.498 0.525 0.558 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveso_dist_o <- dist(scaled_env[ocean_stratified_sites_env,c(3)], method = "euclidean")
jpreseso_mant_o <- mantel(jpreseso_dist, enveso_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
jpreseso_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpreseso_dist, ydis = enveso_dist_o, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6567 
    ##       Significance: 0.356 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.708 0.732 0.749 0.767 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveso_dist_p <- dist(scaled_env[ocean_stratified_sites_env,c(4)], method = "euclidean")
jpreseso_mant_p <- mantel(jpreseso_dist, enveso_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
jpreseso_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpreseso_dist, ydis = enveso_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6796 
    ##       Significance: 0.115 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.684 0.708 0.725 0.745 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpreseso_mant_pv <- rbind(jpreseso_mant_t$signif, jpreseso_mant_s$signif, jpreseso_mant_o$signif, jpreseso_mant_p$signif)
jpreseso_mant_pv <- jpreseso_mant_pv[,1]
jpreseso_mant_pv <- p.adjust(jpreseso_mant_pv, method = "bonferroni")
jpreseso_mant_pv
```

    ## [1] 0.836 0.060 1.000 0.460

``` r
## Dice-Sørensen
# All sites
enve_dist <- dist(scaled_env[surveyed_sites_env,c(2,3,4)], method = "euclidean")
sprese_mant <- mantel(sprese_dist, enve_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
sprese_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = sprese_dist, ydis = enve_dist, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.8057 
    ##       Significance: 0.003 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.703 0.727 0.747 0.768 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
envems_dist <- dist(scaled_env[mixed_stratified_lakes,c(2,3,4)], method = "euclidean")
spresems_mant <- mantel(spresems_dist, envems_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
spresems_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresems_dist, ydis = envems_dist, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7628 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.616 0.646 0.677 0.713 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
enveom_dist <- dist(scaled_env[ocean_mixed_sites_env,c(2,3,4)], method = "euclidean")
spreseom_mant <- mantel(spreseom_dist, enveom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
spreseom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spreseom_dist, ydis = enveom_dist, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6501 
    ##       Significance: 0.008 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.359 0.453 0.533 0.623 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
enveso_dist <- dist(scaled_env[ocean_stratified_sites_env,c(2,3,4)], method = "euclidean")
spreseso_mant <- mantel(spreseso_dist, enveso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
spreseso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spreseso_dist, ydis = enveso_dist, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.8145 
    ##       Significance: 0.026 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.778 0.798 0.814 0.834 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes
enves_dist_t <- dist(scaled_env[stratified_lakes,c(1)], method = "euclidean")
jpreses_mant_t <- mantel(jpress_dist, enves_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[stratified_lakes,19])
jpreses_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpress_dist, ydis = enves_dist_t, method = "spearman",      permutations = 999, strata = env[stratified_lakes, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1252 
    ##       Significance: 0.77 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.215 0.264 0.330 0.416 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enves_dist_s <- dist(scaled_env[stratified_lakes,c(2)], method = "euclidean")
jpreses_mant_s <- mantel(jpress_dist, enves_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[stratified_lakes,19])
jpreses_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpress_dist, ydis = enves_dist_s, method = "spearman",      permutations = 999, strata = env[stratified_lakes, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4001 
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.258 0.316 0.362 0.406 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enves_dist_o <- dist(scaled_env[stratified_lakes,c(3)], method = "euclidean")
jpreses_mant_o <- mantel(jpress_dist, enves_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[stratified_lakes,19])
jpreses_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpress_dist, ydis = enves_dist_o, method = "spearman",      permutations = 999, strata = env[stratified_lakes, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.307 
    ##       Significance: 0.054 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.244 0.316 0.398 0.429 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enves_dist_p <- dist(scaled_env[stratified_lakes,c(4)], method = "euclidean")
jpreses_mant_p <- mantel(jpress_dist, enves_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[stratified_lakes,19])
jpreses_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpress_dist, ydis = enves_dist_p, method = "spearman",      permutations = 999, strata = env[stratified_lakes, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.07194 
    ##       Significance: 0.344 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.218 0.280 0.339 0.389 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpreses_mant_pv <- rbind(jpreses_mant_t$signif, jpreses_mant_s$signif, jpreses_mant_o$signif, jpreses_mant_p$signif)
jpreses_mant_pv <- jpreses_mant_pv[,1]
jpreses_mant_pv <- p.adjust(jpreses_mant_pv, method = "bonferroni")
jpreses_mant_pv
```

    ## [1] 1.000 0.048 0.216 1.000

``` r
### Geographic
## Jaccard
# All sites
envb_dist_vc <- dist(scaled_env[surveyed_sites_geo,c(5)], method = "euclidean")
jpresba_mant_vc <- mantel(jpresb_dist, envb_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_vc, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.05716 
    ##       Significance: 0.637 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.145 0.165 0.179 0.203 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_v <- dist(scaled_env[surveyed_sites_geo,c(6)], method = "euclidean")
jpresba_mant_v <- mantel(jpresb_dist, envb_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_v, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1399 
    ##       Significance: 0.465 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.226 0.251 0.285 0.314 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_sa <- dist(scaled_env[surveyed_sites_geo,c(7)], method = "euclidean")
jpresba_mant_sa <- mantel(jpresb_dist, envb_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_sa, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1295 
    ##       Significance: 0.478 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.210 0.236 0.255 0.280 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_dmin <- dist(scaled_env[surveyed_sites_geo,c(8)], method = "euclidean")
jpresba_mant_dmin <- mantel(jpresb_dist, envb_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_dmin, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.517 
    ##       Significance: 0.065 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.492 0.528 0.551 0.577 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_dmean <- dist(scaled_env[surveyed_sites_geo,c(9)], method = "euclidean")
jpresba_mant_dmean <- mantel(jpresb_dist, envb_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_dmean, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4547 
    ##       Significance: 0.442 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.562 0.589 0.615 0.655 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_dmed <- dist(scaled_env[surveyed_sites_geo,c(10)], method = "euclidean")
jpresba_mant_dmed <- mantel(jpresb_dist, envb_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_dmed, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4347 
    ##       Significance: 0.43 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.538 0.582 0.612 0.645 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_tl <- dist(scaled_env[surveyed_sites_geo,c(11)], method = "euclidean")
jpresba_mant_tl <- mantel(jpresb_dist, envb_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_tl, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.449 
    ##       Significance: 0.614 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.566 0.595 0.619 0.645 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_te <- dist(scaled_env[surveyed_sites_geo,c(12)], method = "euclidean")
jpresba_mant_te <- mantel(jpresb_dist, envb_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_te, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4672 
    ##       Significance: 0.385 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.540 0.573 0.593 0.616 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_p <- dist(scaled_env[surveyed_sites_geo,c(13)], method = "euclidean")
jpresba_mant_p <- mantel(jpresb_dist, envb_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_p, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.07963 
    ##       Significance: 0.532 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.149 0.162 0.181 0.193 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_md <- dist(scaled_env[surveyed_sites_geo,c(14)], method = "euclidean")
jpresba_mant_md <- mantel(jpresb_dist, envb_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_md, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.03459 
    ##       Significance: 0.761 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.200 0.236 0.268 0.293 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_la <- dist(scaled_env[surveyed_sites_geo,c(15)], method = "euclidean")
jpresba_mant_la <- mantel(jpresb_dist, envb_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
jpresba_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist_la, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.07283 
    ##       Significance: 0.429 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.125 0.140 0.155 0.171 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresba__mant_pv <- rbind(jpresba_mant_vc$signif, jpresba_mant_v$signif, jpresba_mant_sa$signif, jpresba_mant_dmin$signif, jpresba_mant_dmean$signif, jpresba_mant_dmed$signif, jpresba_mant_tl$signif, jpresba_mant_te$signif, jpresba_mant_p$signif, jpresba_mant_md$signif, jpresba_mant_la$signif)
jpresba__mant_pv <- jpresba__mant_pv[,1]
jpresba__mant_pv <- p.adjust(jpresba__mant_pv, method = "bonferroni")
jpresba__mant_pv
```

    ##  [1] 1.000 1.000 1.000 0.715 1.000 1.000 1.000 1.000 1.000 1.000 1.000

``` r
# Mixed and stratified lakes
envbms_dist_vc <- dist(scaled_env[mixed_stratified_lakes_geo,c(5)], method = "euclidean")
jpresbms_mant_vc <- mantel(jpresbms_dist, envbms_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_vc, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.09114 
    ##       Significance: 0.822 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.120 0.168 0.193 0.230 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_v <- dist(scaled_env[mixed_stratified_lakes_geo,c(6)], method = "euclidean")
jpresbms_mant_v <- mantel(jpresbms_dist, envbms_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_v, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.09354 
    ##       Significance: 0.658 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.316 0.375 0.462 0.509 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_sa <- dist(scaled_env[mixed_stratified_lakes_geo,c(7)], method = "euclidean")
jpresbms_mant_sa <- mantel(jpresbms_dist, envbms_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_sa, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.06898 
    ##       Significance: 0.638 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.278 0.318 0.398 0.441 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_dmin <- dist(scaled_env[mixed_stratified_lakes_geo,c(8)], method = "euclidean")
jpresbms_mant_dmin <- mantel(jpresbms_dist, envbms_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_dmin, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.344 
    ##       Significance: 0.088 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.332 0.389 0.422 0.464 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_dmean <- dist(scaled_env[mixed_stratified_lakes_geo,c(9)], method = "euclidean")
jpresbms_mant_dmean <- mantel(jpresbms_dist, envbms_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_dmean, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2753 
    ##       Significance: 0.38 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.421 0.464 0.502 0.546 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_dmed <- dist(scaled_env[mixed_stratified_lakes_geo,c(10)], method = "euclidean")
jpresbms_mant_dmed <- mantel(jpresbms_dist, envbms_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_dmed, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2789 
    ##       Significance: 0.33 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.397 0.441 0.488 0.525 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_tl <- dist(scaled_env[mixed_stratified_lakes_geo,c(11)], method = "euclidean")
jpresbms_mant_tl <- mantel(jpresbms_dist, envbms_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_tl, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2164 
    ##       Significance: 0.511 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.375 0.420 0.454 0.508 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_te <- dist(scaled_env[mixed_stratified_lakes_geo,c(12)], method = "euclidean")
jpresbms_mant_te <- mantel(jpresbms_dist, envbms_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_te, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2201 
    ##       Significance: 0.438 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.362 0.400 0.429 0.448 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_p <- dist(scaled_env[mixed_stratified_lakes_geo,c(13)], method = "euclidean")
jpresbms_mant_p <- mantel(jpresbms_dist, envbms_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_p, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.04053 
    ##       Significance: 0.697 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.122 0.155 0.178 0.223 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_md <- dist(scaled_env[mixed_stratified_lakes_geo,c(14)], method = "euclidean")
jpresbms_mant_md <- mantel(jpresbms_dist, envbms_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_md, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.03234 
    ##       Significance: 0.861 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.230 0.279 0.316 0.345 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_la <- dist(scaled_env[mixed_stratified_lakes_geo,c(15)], method = "euclidean")
jpresbms_mant_la <- mantel(jpresbms_dist, envbms_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
jpresbms_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist_la, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02183 
    ##       Significance: 0.783 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.139 0.171 0.196 0.219 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresbms__mant_pv <- rbind(jpresbms_mant_vc$signif, jpresbms_mant_v$signif, jpresbms_mant_sa$signif, jpresbms_mant_dmin$signif, jpresbms_mant_dmean$signif, jpresbms_mant_dmed$signif, jpresbms_mant_tl$signif, jpresbms_mant_te$signif, jpresbms_mant_p$signif, jpresbms_mant_md$signif, jpresbms_mant_la$signif)
jpresbms__mant_pv <- jpresbms__mant_pv[,1]
jpresbms__mant_pv <- p.adjust(jpresbms__mant_pv, method = "bonferroni")
jpresbms__mant_pv
```

    ##  [1] 1.000 1.000 1.000 0.968 1.000 1.000 1.000 1.000 1.000 1.000 1.000

``` r
# Ocean sites and mixed lakes
envbom_dist_vc <- dist(scaled_env[ocean_mixed_sites_geo,c(5)], method = "euclidean")
jpresbom_mant_vc <- mantel(jpresbom_dist, envbom_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_vc, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.211 
    ##       Significance: 0.18 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.249 0.296 0.320 0.366 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_v <- dist(scaled_env[ocean_mixed_sites_geo,c(6)], method = "euclidean")
jpresbom_mant_v <- mantel(jpresbom_dist, envbom_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_v, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1989 
    ##       Significance: 0.219 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.257 0.301 0.343 0.393 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_sa <- dist(scaled_env[ocean_mixed_sites_geo,c(7)], method = "euclidean")
jpresbom_mant_sa <- mantel(jpresbom_dist, envbom_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_sa, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.143 
    ##       Significance: 0.339 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.238 0.267 0.290 0.315 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_dmin <- dist(scaled_env[ocean_mixed_sites_geo,c(8)], method = "euclidean")
jpresbom_mant_dmin <- mantel(jpresbom_dist, envbom_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_dmin, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2562 
    ##       Significance: 0.06 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.213 0.264 0.294 0.314 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_dmean <- dist(scaled_env[ocean_mixed_sites_geo,c(9)], method = "euclidean")
jpresbom_mant_dmean <- mantel(jpresbom_dist, envbom_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_dmean, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.08771 
    ##       Significance: 0.491 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.248 0.300 0.354 0.385 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_dmed <- dist(scaled_env[ocean_mixed_sites_geo,c(10)], method = "euclidean")
jpresbom_mant_dmed <- mantel(jpresbom_dist, envbom_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_dmed, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.05185 
    ##       Significance: 0.592 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.258 0.308 0.378 0.424 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_tl <- dist(scaled_env[ocean_mixed_sites_geo,c(11)], method = "euclidean")
jpresbom_mant_tl <- mantel(jpresbom_dist, envbom_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_tl, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.01695 
    ##       Significance: 0.361 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.165 0.217 0.259 0.307 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_te <- dist(scaled_env[ocean_mixed_sites_geo,c(12)], method = "euclidean")
jpresbom_mant_te <- mantel(jpresbom_dist, envbom_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_te, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.04564 
    ##       Significance: 0.322 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.151 0.201 0.261 0.279 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_p <- dist(scaled_env[ocean_mixed_sites_geo,c(13)], method = "euclidean")
jpresbom_mant_p <- mantel(jpresbom_dist, envbom_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1691 
    ##       Significance: 0.333 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.254 0.280 0.304 0.329 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_md <- dist(scaled_env[ocean_mixed_sites_geo,c(14)], method = "euclidean")
jpresbom_mant_md <- mantel(jpresbom_dist, envbom_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_md, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1423 
    ##       Significance: 0.111 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.153 0.239 0.303 0.383 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_la <- dist(scaled_env[ocean_mixed_sites_geo,c(15)], method = "euclidean")
jpresbom_mant_la <- mantel(jpresbom_dist, envbom_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
jpresbom_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist_la, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1925 
    ##       Significance: 0.13 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.208 0.240 0.282 0.320 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresbom__mant_pv <- rbind(jpresbom_mant_vc$signif, jpresbom_mant_v$signif, jpresbom_mant_sa$signif, jpresbom_mant_dmin$signif, jpresbom_mant_dmean$signif, jpresbom_mant_dmed$signif, jpresbom_mant_tl$signif, jpresbom_mant_te$signif, jpresbom_mant_p$signif, jpresbom_mant_md$signif, jpresbom_mant_la$signif)
jpresbom__mant_pv <- jpresbom__mant_pv[,1]
jpresbom__mant_pv <- p.adjust(jpresbom__mant_pv, method = "bonferroni")
jpresbom__mant_pv
```

    ##  [1] 1.00 1.00 1.00 0.66 1.00 1.00 1.00 1.00 1.00 1.00 1.00

``` r
# Stratified lakes and ocean sites
envbso_dist_vc <- dist(scaled_env[ocean_stratified_sites,c(5)], method = "euclidean")
jpresbso_mant_vc <- mantel(jpresbso_dist, envbso_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_vc, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1298 
    ##       Significance: 0.463 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.186 0.198 0.210 0.222 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_v <- dist(scaled_env[ocean_stratified_sites,c(6)], method = "euclidean")
jpresbso_mant_v <- mantel(jpresbso_dist, envbso_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_v, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1665 
    ##       Significance: 0.479 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.222 0.234 0.244 0.260 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_sa <- dist(scaled_env[ocean_stratified_sites,c(7)], method = "euclidean")
jpresbso_mant_sa <- mantel(jpresbso_dist, envbso_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_sa, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2151 
    ##       Significance: 0.282 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.248 0.266 0.279 0.294 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_dmin <- dist(scaled_env[ocean_stratified_sites,c(8)], method = "euclidean")
jpresbso_mant_dmin <- mantel(jpresbso_dist, envbso_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_dmin, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5463 
    ##       Significance: 0.127 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.561 0.591 0.606 0.644 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_dmean <- dist(scaled_env[ocean_stratified_sites,c(9)], method = "euclidean")
jpresbso_mant_dmean <- mantel(jpresbso_dist, envbso_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_dmean, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5794 
    ##       Significance: 0.476 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.673 0.701 0.724 0.756 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_dmed <- dist(scaled_env[ocean_stratified_sites,c(10)], method = "euclidean")
jpresbso_mant_dmed <- mantel(jpresbso_dist, envbso_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_dmed, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5612 
    ##       Significance: 0.49 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.659 0.690 0.705 0.730 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_tl <- dist(scaled_env[ocean_stratified_sites,c(11)], method = "euclidean")
jpresbso_mant_tl <- mantel(jpresbso_dist, envbso_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_tl, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6658 
    ##       Significance: 0.696 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.753 0.768 0.789 0.799 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_te <- dist(scaled_env[ocean_stratified_sites,c(12)], method = "euclidean")
jpresbso_mant_te <- mantel(jpresbso_dist, envbso_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_te, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7086 
    ##       Significance: 0.356 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.759 0.777 0.790 0.806 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_p <- dist(scaled_env[ocean_stratified_sites,c(13)], method = "euclidean")
jpresbso_mant_p <- mantel(jpresbso_dist, envbso_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2634 
    ##       Significance: 0.452 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.330 0.342 0.352 0.359 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_md <- dist(scaled_env[ocean_stratified_sites,c(14)], method = "euclidean")
jpresbso_mant_md <- mantel(jpresbso_dist, envbso_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_md, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.00291 
    ##       Significance: 0.753 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.182 0.221 0.248 0.289 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_la <- dist(scaled_env[ocean_stratified_sites,c(15)], method = "euclidean")
jpresbso_mant_la <- mantel(jpresbso_dist, envbso_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jpresbso_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist_la, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2524 
    ##       Significance: 0.217 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.286 0.301 0.319 0.336 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresbso__mant_pv <- rbind(jpresbso_mant_vc$signif, jpresbso_mant_v$signif, jpresbso_mant_sa$signif, jpresbso_mant_dmin$signif, jpresbso_mant_dmean$signif, jpresbso_mant_dmed$signif, jpresbso_mant_tl$signif, jpresbso_mant_te$signif, jpresbso_mant_p$signif, jpresbso_mant_md$signif, jpresbso_mant_la$signif)
jpresbso__mant_pv <- jpresbso__mant_pv[,1]
jpresbso__mant_pv <- p.adjust(jpresbso__mant_pv, method = "bonferroni")
jpresbso__mant_pv
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1

``` r
## Dice-Sørensen
# All sites
envb_dist <- dist(scaled_env[surveyed_sites_geo,c(5:15)], method = "euclidean")
spresb_mant <- mantel(spresb_dist, envb_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
spresb_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresb_dist, ydis = envb_dist, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3881 
    ##       Significance: 0.269 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.430 0.461 0.486 0.505 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
envbms_dist <- dist(scaled_env[mixed_stratified_lakes_geo,c(5:15)], method = "euclidean")
spresbms_mant <- mantel(spresbms_dist, envbms_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
spresbms_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresbms_dist, ydis = envbms_dist, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2892 
    ##       Significance: 0.356 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.409 0.446 0.500 0.555 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
envbom_dist <- dist(scaled_env[ocean_mixed_sites_geo,c(5:15)], method = "euclidean")
spresbom_mant <- mantel(spresbom_dist, envbom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
spresbom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresbom_dist, ydis = envbom_dist, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2261 
    ##       Significance: 0.205 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.284 0.331 0.391 0.430 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
envbso_dist <- dist(scaled_env[ocean_stratified_sites,c(5:15)], method = "euclidean")
spresbso_mant <- mantel(spresbso_dist, envbso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
spresbso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresbso_dist, ydis = envbso_dist, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5571 
    ##       Significance: 0.393 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.602 0.617 0.628 0.643 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Traits
# All sites
FDpa_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(1)]))
jpresta_mant_s <- mantel(jpres_dist, FDpa_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_s, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5576 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.370 0.408 0.436 0.468 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(2)]))
jpresta_mant_l <- mantel(jpres_dist, FDpa_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_l, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3708 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.261 0.284 0.310 0.327 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_t <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(3)]))
jpresta_mant_t <- mantel(jpres_dist, FDpa_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_t, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4872 
    ##       Significance: 0.028 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.437 0.468 0.489 0.512 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(4)]))
jpresta_mant_dmin <- mantel(jpres_dist, FDpa_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_dmin, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6271 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.558 0.575 0.588 0.599 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_dmax <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(5)]))
jpresta_mant_dmax <- mantel(jpres_dist, FDpa_dist_dmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_dmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_dmax, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2736 
    ##       Significance: 0.176 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.300 0.328 0.344 0.364 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_tmin <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(6)]))
jpresta_mant_tmin <- mantel(jpres_dist, FDpa_dist_tmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_tmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_tmin, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.343 
    ##       Significance: 0.18 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.369 0.390 0.409 0.434 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_tmax <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(7)]))
jpresta_mant_tmax <- mantel(jpres_dist, FDpa_dist_tmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_tmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_tmax, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1955 
    ##       Significance: 0.086 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.189 0.224 0.241 0.267 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_b <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(8:9)]))
jpresta_mant_b <- mantel(jpres_dist, FDpa_dist_b, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_b
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_b, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4724 
    ##       Significance: 0.029 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.408 0.437 0.478 0.503 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_o <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(10)]))
jpresta_mant_o <- mantel(jpres_dist, FDpa_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_o, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4377 
    ##       Significance: 0.092 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.435 0.475 0.484 0.519 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_f <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(11)]))
jpresta_mant_f <- mantel(jpres_dist, FDpa_dist_f, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_f
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_f, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2285 
    ##       Significance: 0.249 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.267 0.276 0.346 0.346 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(12:13)]))
jpresta_mant_ec <- mantel(jpres_dist, FDpa_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_ec, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.401 
    ##       Significance: 0.011 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.310 0.350 0.381 0.392 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(14:15)]))
jpresta_mant_es <- mantel(jpres_dist, FDpa_dist_es, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_es
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_es, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.401 
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.320 0.355 0.381 0.401 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_p <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(16)]))
jpresta_mant_p <- mantel(jpres_dist, FDpa_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_p, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7067 
    ##       Significance: 0.122 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.707 0.707 0.707 0.707 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_w <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(17)]))
jpresta_mant_w <- mantel(jpres_dist, FDpa_dist_w, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
jpresta_mant_w
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_w, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5421 
    ##       Significance: 0.058 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.495 0.542 0.600 0.608 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresta_mant_pv <- rbind(jpresta_mant_s$signif, jpresta_mant_l$signif, jpresta_mant_t$signif, jpresta_mant_dmin$signif, jpresta_mant_dmax$signif, jpresta_mant_tmin$signif, jpresta_mant_tmax$signif, jpresta_mant_b$signif, jpresta_mant_o$signif, jpresta_mant_f$signif, jpresta_mant_ec$signif, jpresta_mant_es$signif, jpresta_mant_p$signif, jpresta_mant_w$signif)
jpresta_mant_pv <- jpresta_mant_pv[,1]
jpresta_mant_pv <- p.adjust(jpresta_mant_pv, method = "bonferroni")
jpresta_mant_pv
```

    ##  [1] 0.014 0.028 0.392 0.014 1.000 1.000 1.000 0.406 1.000 1.000 0.154 0.168
    ## [13] 1.000 0.812

``` r
# Mixed and stratified lakes
FDpms_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(1)]))
jprestms_mant_s <- mantel(jpresms_dist, FDpms_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_s, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.483 
    ##       Significance: 0.004 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.312 0.365 0.417 0.444 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(2)]))
jprestms_mant_l <- mantel(jpresms_dist, FDpms_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_l, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2828 
    ##       Significance: 0.028 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.223 0.257 0.286 0.312 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_t <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(3)]))
jprestms_mant_t <- mantel(jpresms_dist, FDpms_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_t, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3861 
    ##       Significance: 0.037 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.324 0.367 0.408 0.443 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(4)]))
jprestms_mant_dmin <- mantel(jpresms_dist, FDpms_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_dmin, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5503 
    ##       Significance: 0.011 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.476 0.503 0.531 0.550 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_dmax <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(5)]))
jprestms_mant_dmax <- mantel(jpresms_dist, FDpms_dist_dmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_dmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_dmax, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r:  0.28 
    ##       Significance: 0.102 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.281 0.318 0.350 0.372 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_tmin <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(6)]))
jprestms_mant_tmin <- mantel(jpresms_dist, FDpms_dist_tmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_tmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_tmin, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2497 
    ##       Significance: 0.164 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.281 0.321 0.343 0.371 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_tmax <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(7)]))
jprestms_mant_tmax <- mantel(jpresms_dist, FDpms_dist_tmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_tmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_tmax, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1266 
    ##       Significance: 0.152 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.150 0.188 0.221 0.263 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_b <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(8:9)]))
jprestms_mant_b <- mantel(jpresms_dist, FDpms_dist_b, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_b
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_b, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5937 
    ##       Significance: 0.016 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.421 0.481 0.550 0.603 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_o <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(10)]))
jprestms_mant_o <- mantel(jpresms_dist, FDpms_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_o, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4223 
    ##       Significance: 0.086 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.414 0.468 0.483 0.523 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_f <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(11)]))
jprestms_mant_f <- mantel(jpresms_dist, FDpms_dist_f, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_f
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_f, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1818 
    ##       Significance: 0.255 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.224 0.244 0.337 0.337 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(12:13)]))
jprestms_mant_ec <- mantel(jpresms_dist, FDpms_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_ec, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3791 
    ##       Significance: 0.01 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.265 0.316 0.361 0.373 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(14:15)]))
jprestms_mant_es <- mantel(jpresms_dist, FDpms_dist_es, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_es
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_es, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3791 
    ##       Significance: 0.016 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.264 0.316 0.361 0.379 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_p <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(16)]))
jprestms_mant_p <- mantel(jpresms_dist, FDpms_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_p, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7242 
    ##       Significance: 0.131 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.724 0.724 0.724 0.724 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpms_dist_w <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(17)]))
jprestms_mant_w <- mantel(jpresms_dist, FDpms_dist_w, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
jprestms_mant_w
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpms_dist_w, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5217 
    ##       Significance: 0.054 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.454 0.522 0.614 0.641 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jprestms_mant_pv <- rbind(jprestms_mant_s$signif, jprestms_mant_l$signif, jprestms_mant_t$signif, jprestms_mant_dmin$signif, jprestms_mant_dmax$signif, jprestms_mant_tmin$signif, jprestms_mant_tmax$signif, jprestms_mant_b$signif, jprestms_mant_o$signif, jprestms_mant_f$signif, jprestms_mant_ec$signif, jprestms_mant_es$signif, jprestms_mant_p$signif, jprestms_mant_w$signif)
jprestms_mant_pv <- jprestms_mant_pv[,1]
jprestms_mant_pv <- p.adjust(jprestms_mant_pv, method = "bonferroni")
jprestms_mant_pv
```

    ##  [1] 0.056 0.392 0.518 0.154 1.000 1.000 1.000 0.224 1.000 1.000 0.140 0.224
    ## [13] 1.000 0.756

``` r
# Ocean sites and mixed lakes
FDpom_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(1)]))
jprestom_mant_s <- mantel(jpresom_dist, FDpom_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
jprestom_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresom_dist, ydis = FDpom_dist_s, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4413 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.161 0.221 0.262 0.347 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpom_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(2)]))
jprestom_mant_l <- mantel(jpresom_dist, FDpom_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
jprestom_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresom_dist, ydis = FDpom_dist_l, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2336 
    ##       Significance: 0.025 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.143 0.183 0.224 0.283 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpom_dist_t <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(3)]))
jprestom_mant_t <- mantel(jpresom_dist, FDpom_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
jprestom_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresom_dist, ydis = FDpom_dist_t, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2598 
    ##       Significance: 0.083 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.239 0.291 0.351 0.374 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpom_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(4)]))
jprestom_mant_dmin <- mantel(jpresom_dist, FDpom_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
jprestom_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresom_dist, ydis = FDpom_dist_dmin, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.01235 
    ##       Significance: 0.314 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.117 0.212 0.289 0.348 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpom_dist_dmax <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(5)]))
jprestom_mant_dmax <- mantel(jpresom_dist, FDpom_dist_dmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
jprestom_mant_dmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresom_dist, ydis = FDpom_dist_dmax, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.009954 
    ##       Significance: 0.564 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.219 0.277 0.321 0.378 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpom_dist_tmin <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(6)]))
jprestom_mant_tmin <- mantel(jpresom_dist, FDpom_dist_tmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
jprestom_mant_tmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresom_dist, ydis = FDpom_dist_tmin, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1988 
    ##       Significance: 0.04 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.107 0.162 0.227 0.287 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpom_dist_tmax <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(7)]))
jprestom_mant_tmax <- mantel(jpresom_dist, FDpom_dist_tmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
jprestom_mant_tmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresom_dist, ydis = FDpom_dist_tmax, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2742 
    ##       Significance: 0.025 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.157 0.214 0.271 0.346 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# FDpom_dist_b <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(8:9)]))
# jprestom_mant_b <- mantel(jpresom_dist, FDpom_dist_b, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# jprestom_mant_b
# Error in cor(as.vector(xdis), ydis, method = method, use = use) : 
#   no complete element pairs
# FDpom_dist_o <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(10)]))
# jprestom_mant_o <- mantel(jpresom_dist, FDpom_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# jprestom_mant_o
# Error in cor(as.vector(xdis), ydis, method = method, use = use) : 
#   no complete element pairs
# FDpom_dist_f <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(11)]))
# jprestom_mant_f <- mantel(jpresom_dist, FDpom_dist_f, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# jprestom_mant_f
# Error in cor(as.vector(xdis), ydis, method = method, use = use) : 
#   no complete element pairs
# FDpom_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(12:13)]))
# jprestom_mant_ec <- mantel(jpresom_dist, FDpom_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# jprestom_mant_ec
# Error in cor(as.vector(xdis), ydis, method = method, use = use) : 
#   no complete element pairs
# FDpom_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(14:15)]))
# jprestom_mant_es <- mantel(jpresom_dist, FDpom_dist_es, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# jprestom_mant_es
# Error in cor(as.vector(xdis), ydis, method = method, use = use) : 
#   no complete element pairs
# FDpom_dist_p <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(16)]))
# jprestom_mant_p <- mantel(jpresom_dist, FDpom_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# jprestom_mant_p
# Error in cor(as.vector(xdis), ydis, method = method, use = use) : 
#   no complete element pairs
# FDpom_dist_w <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(17)]))
# jprestom_mant_w <- mantel(jpresom_dist, FDpom_dist_w, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# jprestom_mant_w
# Error in cor(as.vector(xdis), ydis, method = method, use = use) : 
#   no complete element pairs
# Adjust p-values
jprestom_mant_pv <- rbind(jprestom_mant_s$signif, jprestom_mant_l$signif, jprestom_mant_t$signif, jprestom_mant_dmin$signif, jprestom_mant_dmax$signif, jprestom_mant_tmin$signif, jprestom_mant_tmax$signif)
jprestom_mant_pv <- jprestom_mant_pv[,1]
jprestom_mant_pv <- p.adjust(jprestom_mant_pv, method = "bonferroni")
jprestom_mant_pv
```

    ## [1] 0.014 0.175 0.581 1.000 1.000 0.280 0.175

``` r
# Stratified lakes and ocean sites
FDpso_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(1)]))
jprestso_mant_s <- mantel(jpresso_dist, FDpso_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_s, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.474 
    ##       Significance: 0.014 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.367 0.410 0.433 0.482 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(2)]))
jprestso_mant_l <- mantel(jpresso_dist, FDpso_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_l, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2556 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.154 0.180 0.199 0.218 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_t <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(3)]))
jprestso_mant_t <- mantel(jpresso_dist, FDpso_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_t, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6994 
    ##       Significance: 0.053 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.679 0.701 0.730 0.750 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(4)]))
jprestso_mant_dmin <- mantel(jpresso_dist, FDpso_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_dmin, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6044 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.512 0.525 0.539 0.557 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_dmax <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(5)]))
jprestso_mant_dmax <- mantel(jpresso_dist, FDpso_dist_dmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_dmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_dmax, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.03229 
    ##       Significance: 0.232 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0647 0.0903 0.1076 0.1279 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_tmin <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(6)]))
jprestso_mant_tmin <- mantel(jpresso_dist, FDpso_dist_tmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_tmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_tmin, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1437 
    ##       Significance: 0.495 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.243 0.282 0.307 0.334 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_tmax <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(7)]))
jprestso_mant_tmax <- mantel(jpresso_dist, FDpso_dist_tmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_tmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_tmax, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.004175 
    ##       Significance: 0.208 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0350 0.0633 0.0857 0.1145 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_b <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(8:9)]))
jprestso_mant_b <- mantel(jpresso_dist, FDpso_dist_b, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_b
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_b, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2052 
    ##       Significance: 0.051 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.163 0.205 0.233 0.261 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_o <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(10)]))
jprestso_mant_o <- mantel(jpresso_dist, FDpso_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_o, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.226 
    ##       Significance: 0.113 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.230 0.267 0.277 0.300 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_f <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(11)]))
jprestso_mant_f <- mantel(jpresso_dist, FDpso_dist_f, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_f
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_f, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.04013 
    ##       Significance: 0.152 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0478 0.0487 0.1161 0.1161 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(12:13)]))
jprestso_mant_ec <- mantel(jpresso_dist, FDpso_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_ec, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.147 
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0743 0.1040 0.1204 0.1470 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(14:15)]))
jprestso_mant_es <- mantel(jpresso_dist, FDpso_dist_es, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_es
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_es, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.147 
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0733 0.1027 0.1194 0.1470 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_p <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(16)]))
jprestso_mant_p <- mantel(jpresso_dist, FDpso_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6938 
    ##       Significance: 0.124 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.694 0.694 0.694 0.694 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpso_dist_w <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(17)]))
jprestso_mant_w <- mantel(jpresso_dist, FDpso_dist_w, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
jprestso_mant_w
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpso_dist_w, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3882 
    ##       Significance: 0.051 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.347 0.362 0.432 0.469 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jprestso_mant_pv <- rbind(jprestso_mant_s$signif, jprestso_mant_l$signif, jprestso_mant_t$signif, jprestso_mant_dmin$signif, jprestso_mant_dmax$signif, jprestso_mant_tmin$signif, jprestso_mant_tmax$signif, jprestso_mant_b$signif, jprestso_mant_o$signif, jprestso_mant_f$signif, jprestso_mant_ec$signif, jprestso_mant_es$signif, jprestso_mant_p$signif, jprestso_mant_w$signif)
jprestso_mant_pv <- jprestso_mant_pv[,1]
jprestso_mant_pv <- p.adjust(jprestso_mant_pv, method = "bonferroni")
jprestso_mant_pv
```

    ##  [1] 0.196 0.014 0.742 0.028 1.000 1.000 1.000 0.714 1.000 1.000 0.168 0.168
    ## [13] 1.000 0.714

``` r
## Dice-Sørensen
# All sites
FDpa_dist <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(1:17)]))
spresta_mant <- mantel(spres_dist, FDpa_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
spresta_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spres_dist, ydis = FDpa_dist, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6454 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.517 0.547 0.564 0.587 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
FDpams_dist <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(1:17)]))
sprestms_mant <- mantel(spresms_dist, FDpams_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
sprestms_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresms_dist, ydis = FDpams_dist, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.664 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.460 0.510 0.532 0.558 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
FDpaom_dist <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(1:7)]))
sprestom_mant <- mantel(spresom_dist, FDpaom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
sprestom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresom_dist, ydis = FDpaom_dist, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3897 
    ##       Significance: 0.006 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.171 0.243 0.275 0.342 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
FDpaso_dist <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(1:17)]))
sprestso_mant <- mantel(spresso_dist, FDpaso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
sprestso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresso_dist, ydis = FDpaso_dist, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4016 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.271 0.299 0.336 0.357 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

## Plot Jaccard incidence NMDS and envfit results

- Includes a plot of correlated environmental variables

``` r
jpres_NMDS_data.scores <- as.data.frame(scores(jpres_NMDS))
jpres_NMDS_data.scores$Stratification <- env[surveyed_sites,19]
jpres_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
jpres_NMDS_data.scores$Stratification <- factor(jpres_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

jprese_ef_coord_cont <- as.data.frame(scores(jprese_ef, "vectors")) * ordiArrowMul(jprese_ef)
jprese_new_row_names <- c("pH")
# Assign the new row names to the data frame
jprese_ef_coord_cont <- data.frame(row.names = jprese_new_row_names, jprese_ef_coord_cont)

# jpresb_ef_coord_cont <- as.data.frame(scores(jpresb_ef, "vectors")) * ordiArrowMul(jpresb_ef)
# jpresb_new_row_names <- c("LA")
# # Assign the new row names to the data frame
# jpresb_ef_coord_cont <- data.frame(row.names = jpresb_new_row_names, jpresb_ef_coord_cont)

jpres_ef_plot <- ggplot(data = jpres_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.95) +
  geom_point(data = jpres_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = 2.5*NMDS1, yend = 2.5*NMDS2),
               data = jprese_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = jprese_ef_coord_cont, aes(x = 2.5*NMDS1, y = 2.5*NMDS2),
            color = env_cont, label = row.names(jprese_ef_coord_cont), size = 5) +
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
  #              data = jpresb_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#698B22") +
  # geom_text(data = jpresb_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
  #           color = "#698B22", label = row.names(jpresb_ef_coord_cont), size = 6) +
    geom_text_repel(data = jpres_NMDS_data.scores, label = jpres_NMDS_data.scores$Lakes, 
                    size = 5, force = 10, point.padding = 5, max.overlaps = 30, nudge_x = -0.01) +
  theme(text = element_text(size = 22), 
        legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank()) +
  annotate("text", x = -1.7, y = 0.92, size = 5,
           label = paste("Stress: ", round(jpres_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ")
jpres_ef_plot <- jpres_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
jpres_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20Jaccard%20incidence%20NMDS%20and%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jpres_ef_plot.png", jpres_ef_plot, width = 6, height = 4, units = "in")

# For 90 percent confidence interval: x = -1.55, y = 0.92

# Determine outliers
ordered_jpres_NMDS1_data.scores <- jpres_NMDS_data.scores[order(jpres_NMDS_data.scores$NMDS1), ]
Q1 <- ordered_jpres_NMDS1_data.scores[6,1]
Q3 <- ordered_jpres_NMDS1_data.scores[18,1]
IQR1 <- IQR(jpres_NMDS_data.scores$NMDS1)
Q1 - 1.5*IQR1
```

    ## [1] -3.153562

``` r
Q3 + 1.5*IQR1
```

    ## [1] 3.061931

``` r
ordered_jpres_NMDS1_data.scores$NMDS1
```

    ##  [1] -1.59317965 -1.53556643 -1.10588844 -0.99472036 -0.95961621 -0.91795845
    ##  [7] -0.53181455 -0.38146764 -0.13652399  0.01126882  0.02219528  0.06217896
    ## [13]  0.46144978  0.56782453  0.57060619  0.65640823  0.67317036  0.82632755
    ## [19]  0.97280768  1.01175808  1.13246109  1.18827918

``` r
boxplot(NMDS1 ~ Stratification, ordered_jpres_NMDS1_data.scores)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20Jaccard%20incidence%20NMDS%20and%20envfit%20results-2.png)<!-- -->

``` r
ordered_jpres_NMDS2_data.scores <- jpres_NMDS_data.scores[order(jpres_NMDS_data.scores$NMDS2), ]
Q1 <- ordered_jpres_NMDS2_data.scores[6,2]
Q3 <- ordered_jpres_NMDS2_data.scores[18,2]
IQR2 <- IQR(jpres_NMDS_data.scores$NMDS2)
Q1 - 1.5*IQR2
```

    ## [1] -0.7755052

``` r
Q3 + 1.5*IQR2
```

    ## [1] 0.7785046

``` r
ordered_jpres_NMDS2_data.scores$NMDS2
```

    ##  [1] -0.89730238 -0.64469542 -0.41056831 -0.26603295 -0.23993447 -0.21680313
    ##  [7] -0.14333603 -0.12925208 -0.03663880 -0.01617379  0.04055497  0.05535307
    ## [13]  0.05755988  0.08552792  0.09947238  0.16900082  0.17570864  0.21980251
    ## [19]  0.37869355  0.43468873  0.53826034  0.74611457

``` r
boxplot(NMDS2 ~ Stratification, ordered_jpres_NMDS2_data.scores)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20Jaccard%20incidence%20NMDS%20and%20envfit%20results-3.png)<!-- -->

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jpres_NMDS_boxplot.tiff", width = 1200, height = 800, res = 150, type = "cairo")

jpres_dissimilarity_matrix <- as.matrix(jpres_dist)

# Set the diagonal elements to NA
diag(jpres_dissimilarity_matrix) <- NA

# Get the labels of the sites
sites <- attr(jpres_dissimilarity_matrix, "Labels")

# Calculate the mean dissimilarity for each group
jpres_mean_dissimilarity <- aggregate(jpres_dissimilarity_matrix, by = list(stratification_group), FUN = mean, na.rm = TRUE)

row.names(jpres_mean_dissimilarity) <- jpres_mean_dissimilarity$Group.1

jpres_mean_dissimilarity <- jpres_mean_dissimilarity[,-1] 

jpres_mean_dissimilarity <- as.data.frame(t(jpres_mean_dissimilarity))

jpres_mean_dissimilarity$Stratification <- env[surveyed_sites,19]

boxplot(Ocean ~ Stratification, jpres_mean_dissimilarity[c(7,9,11,16,18,19),])
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
jpres_mean_dissimilarity$Stratification <- factor(jpres_mean_dissimilarity$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

jpres_ocn_dist_violin_plot <- ggplot(jpres_mean_dissimilarity, mapping = aes(x= Stratification, y= Ocean, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.8,
            width = 0.1) +
  geom_text_repel(data = jpres_mean_dissimilarity, label = row.names(jpres_mean_dissimilarity), size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  ylab("Ocean site distances") +
  xlab("Site type") +
  labs(colour = "Site type:  ", fill = "Site type:  ")
jpres_ocn_dist_violin_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20Jaccard%20incidence%20NMDS%20and%20envfit%20results-4.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jpres_ocn_dist_violin_plot.png", jpres_ocn_dist_violin_plot, width = 8, height = 4, units = "in")

jpres_mix_dist_violin_plot <- ggplot(jpres_mean_dissimilarity, mapping = aes(x= Stratification, y= Mixed, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.8,
            width = 0.1) +
  geom_text_repel(data = jpres_mean_dissimilarity, label = row.names(jpres_mean_dissimilarity), size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  ylab("Mixed lake distances") +
  xlab("Site type") +
  labs(colour = "Site type:  ", fill = "Site type:  ")
jpres_mix_dist_violin_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20Jaccard%20incidence%20NMDS%20and%20envfit%20results-5.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jpres_mix_dist_violin_plot.png", jpres_mix_dist_violin_plot, width = 8, height = 4, units = "in")

jpres_strat_dist_violin_plot <- ggplot(jpres_mean_dissimilarity, mapping = aes(x= Stratification, y= Stratified, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.8,
            width = 0.1) +
  geom_text_repel(data = jpres_mean_dissimilarity, label = row.names(jpres_mean_dissimilarity), size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  ylab("Stratified lake distances") +
  xlab("Site type") +
  labs(colour = "Site type:  ", fill = "Site type:  ")
jpres_strat_dist_violin_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20Jaccard%20incidence%20NMDS%20and%20envfit%20results-6.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jpres_strat_dist_violin_plot.png", jpres_strat_dist_violin_plot, width = 8, height = 4, units = "in")
```

## Boxplot of species intradissimilarity heterogeneity

- The dissimilarity between sites that share the same site type

``` r
jpres_bd_dist <- stratification_group_jpres_bd$distances
jpres_bd_dist <- as.data.frame(jpres_bd_dist)
jpres_bd_dist$X <- row.names(jpres_bd_dist)
jpres_bd_dist_env <- merge(jpres_bd_dist, env[surveyed_sites,], by = "X", sort = F)

jpres_bd_dist_env$Stratification <- factor(jpres_bd_dist_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

jpres_bd_dist_env_plot <- ggplot(jpres_bd_dist_env, aes(x = Stratification, y = jpres_bd_dist, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.8,
            width = 0.1) +
    geom_text_repel(data = jpres_bd_dist_env, label = jpres_bd_dist_env$X, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black")) +
  xlab("Site type") +
  ylab("Distance to Centroid") +
  labs(color = "Stratification", tag = "A")
jpres_bd_dist_env_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Boxplot%20of%20species%20intradissimilarity%20heterogeneity-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jpres_bd_dist_env_plot.png", jpres_bd_dist_env_plot, width = 8.25, height = 4.13, units = "in")
```

## Plot Jaccard incidence NMDS and trait envfit results

Includes a plot of correlated traits

``` r
jprest_ef_coord_cont <- as.data.frame(scores(jprest_ef, "vectors")) * ordiArrowMul(jprest_ef)
jprest_ef_coord_cat = as.data.frame(scores(jprest_ef, "factors")) * ordiArrowMul(jprest_ef)

jprest_ef_plot <- ggplot(data = jpres_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification),type="t", level = 0.95) +
  geom_point(data = jpres_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  geom_text_repel(data = jpres_NMDS_data.scores, label = jpres_NMDS_data.scores$Lakes, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = jprest_ef_coord_cont, linewidth =1, alpha = 0.2, color = trait_cont) +
  geom_text(data = jprest_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
            color = trait_cont, label = row.names(jprest_ef_coord_cont), size = 5) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = jprest_ef_coord_cat, linewidth =1, alpha = 0.2, color = trait_cat) +
  geom_text(data = jprest_ef_coord_cat, aes(x = NMDS1, y = NMDS2), 
            color = trait_cat, label = row.names(jprest_ef_coord_cat), size = 5, check_overlap = T) + 
  theme(text = element_text(size = 22), legend.position = "bottom", legend.text = element_text(size = 16), panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  annotate("text", x = -1.68, y = 0.90, size = 5, 
           label = paste("Stress: ", round(jpres_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "A")
jprest_ef_plot <- jprest_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
jprest_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20Jaccard%20incidence%20NMDS%20and%20trait%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jprest_ef_plot.png", jprest_ef_plot, width = 6, height = 4, units = "in")
```

## Plot Dice-Sørensen incidence NMDS and envfit results

``` r
spres_NMDS_data.scores <- as.data.frame(scores(spres_NMDS))
spres_NMDS_data.scores$Stratification <- env[-c(20),19]
spres_NMDS_data.scores$Stratification <- factor(spres_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

sprese_ef_coord_cont <- as.data.frame(scores(sprese_ef, "vectors")) * ordiArrowMul(sprese_ef)
spresb_ef_coord_cont <- as.data.frame(scores(spresb_ef, "vectors")) * ordiArrowMul(spresb_ef)
sprest_ef_coord_cont <- as.data.frame(scores(sprest_ef, "vectors")) * ordiArrowMul(sprest_ef)
sprest_ef_coord_cat = as.data.frame(scores(sprest_ef, "factors")) * ordiArrowMul(sprest_ef)

spres_ef_plot <- ggplot(data = spres_NMDS_data.scores, aes(x = NMDS1, y = NMDS2)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification),level = 0.50) +
  geom_point(data = spres_NMDS_data.scores, aes(color = Stratification), linewidt = 7, alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sprese_ef_coord_cont, size =1, alpha = 0.5, color = "#8968CD") + 
  geom_text(data = sprese_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
            color = "#8968CD", label = row.names(sprese_ef_coord_cont)) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = spresb_ef_coord_cont, size =1, alpha = 0.5, color = "#CD5555") + 
  geom_text(data = spresb_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
            color = "#CD5555", label = row.names(spresb_ef_coord_cont)) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sprest_ef_coord_cont, size =1, alpha = 0.5, color = "darkolivegreen3") +
  geom_text(data = sprest_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
            color = "darkolivegreen3", label = row.names(sprest_ef_coord_cont)) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sprest_ef_coord_cat, size =1, alpha = 0.5, color = "#EEC900") +
  geom_text(data = sprest_ef_coord_cat, aes(x = NMDS1, y = NMDS2), 
            color = "#EEC900", label = row.names(sprest_ef_coord_cat)) + 
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
        panel.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  annotate("text", x = -1.5, y = 0.82639346, size = 5, 
           label = paste("Stress: ", round(spres_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ")
spres_ef_plot
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/spres_ef_plot.png", spres_ef_plot, width = 8, height = 4, units = "in")
```

## SDis dendrogram

``` r
# cluster communities using average-linkage algorithm
jpres_dist_clust <- hclust(jpres_dist, method = "average")

# Open a PNG device
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jpres_dendro.tiff", width = 1200, height = 800, res = 150, type = "cairo")

# Create your plot using the plot() function
plot(jpres_dist_clust, 
     xlab = "Sites",
     ylab = "Species dissimilarity",
     main = "",
     sub = "",
     pch = 20,
     col= "black")

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jpres_si.tiff", width = 1200, height = 800, res = 150, type = "cairo")

Si <- numeric(nrow(presabs_lake[surveyed_sites,]))
for (k in 2:(nrow(presabs_lake[surveyed_sites,])-1))
{
  sil<-silhouette(cutree(jpres_dist_clust, k=k), jpres_dist)
  Si[k] <- summary(sil)$avg.width
}
k.best<-which.max(Si)
plot(1:nrow(presabs_lake[surveyed_sites,]), Si, type = "h", main = "Silhouette", xlab = "K", ylab = "Width")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col = "red", font = 2, col.axis = "red")

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jpres_pam.tiff", width = 1200, height = 800, res = 150, type = "cairo")
jpres_pam <- pam(jpres_dist,3)
plot(jpres_pam)

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Trait dissimilarity distances

- We use the community weighted means calculated from the FD package. We
  use traits that were found to be significantly correlated with
  stratification and environmental variables.

``` r
### Regular
# DemersPelag has no variation so it will not be included
# RepGuild1 and RepGuild2 are aliases
# RepGuild 1 or 2 and TempPrefMax have too high of vif to be used below
mod <-  lm(FRic_total ~ BodyShapeI + OperculumPresent + MaxLengthTL + Troph + DepthMin + DepthMax + TempPrefMin + FeedingPath + ParentalCare + WaterPref + DorsalSpinesMean, FD_total_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = FRic_total ~ BodyShapeI + OperculumPresent + MaxLengthTL + 
    ##     Troph + DepthMin + DepthMax + TempPrefMin + FeedingPath + 
    ##     ParentalCare + WaterPref + DorsalSpinesMean, data = FD_total_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.29380 -0.05585 -0.00447  0.06114  0.33032 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)          1.3641328  6.9630581   0.196   0.8486  
    ## BodyShapeI3f        -0.0340174  0.1629657  -0.209   0.8388  
    ## BodyShapeI4e        -0.5961886  0.2422615  -2.461   0.0336 *
    ## OperculumPresentyes -0.5577873  0.2485597  -2.244   0.0487 *
    ## MaxLengthTL          0.0159171  0.0104320   1.526   0.1580  
    ## Troph               -0.7876839  0.4560155  -1.727   0.1148  
    ## DepthMin             0.0558398  0.0421961   1.323   0.2152  
    ## DepthMax             0.0002717  0.0090983   0.030   0.9768  
    ## TempPrefMin          0.1047607  0.2619456   0.400   0.6976  
    ## FeedingPathp        -0.5199859  0.4258663  -1.221   0.2501  
    ## ParentalCare4n       0.1371361  0.2987305   0.459   0.6560  
    ## WaterPref3a         -0.5115436  0.2722362  -1.879   0.0897 .
    ## DorsalSpinesMean    -0.0410396  0.0870579  -0.471   0.6475  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.194 on 10 degrees of freedom
    ## Multiple R-squared:  0.9418, Adjusted R-squared:  0.872 
    ## F-statistic: 13.49 on 12 and 10 DF,  p-value: 0.0001313

``` r
rms::vif(mod)
```

    ##        BodyShapeI3f        BodyShapeI4e OperculumPresentyes         MaxLengthTL 
    ##            3.436668            5.153591            6.424385            6.673611 
    ##               Troph            DepthMin            DepthMax         TempPrefMin 
    ##            3.874971            5.268850            8.127159            6.662467 
    ##        FeedingPathp      ParentalCare4n         WaterPref3a    DorsalSpinesMean 
    ##            8.800822           11.547939            7.706577            5.505648

``` r
car::vif(mod)
```

    ##                       GVIF Df GVIF^(1/(2*Df))
    ## BodyShapeI        8.870620  2        1.725792
    ## OperculumPresent  6.424385  1        2.534637
    ## MaxLengthTL       6.673611  1        2.583333
    ## Troph             3.874971  1        1.968495
    ## DepthMin          5.268850  1        2.295398
    ## DepthMax          8.127159  1        2.850817
    ## TempPrefMin       6.662467  1        2.581176
    ## FeedingPath       8.800822  1        2.966618
    ## ParentalCare     11.547939  1        3.398226
    ## WaterPref         7.706577  1        2.776072
    ## DorsalSpinesMean  5.505648  1        2.346412

``` r
# Between all sites
# BodyShapeI 0.30 0.009
# OperculumPresent 0.39 0.021
# Troph 0.69 0.001
# ParentalCare 0.82 0.001
# WaterPref 0.51 0.008
# DorsalSpinesMean 0.46 0.002

# Between stratified and mixed
# Troph 0.44 0.015
# ParentalCare 0.78 0.006
# DorsalSprinesMean 0.38 0.036

# Between stratified and ocean
# Troph 0.75 0.003
# ParentalCare 0.75 0.008
# DorsalSpinesMean 0.48 0.018

# Between mixed and ocean
# Troph 0.53 0.015

# Traits that violate the test of homogeneity of dispersion with individual p-value
# OperculumPresent 2.2e-16, MaxLengthTL 0.004, DepthMax 0.005, WaterPref 0.046
# Traits that individually do not violate the test of homogeneity of dispersion but make the group violate the test with individual p-value 
# TempPrefMin 0.12, FeedingPath 0.16

# Trait combo that does not violate the test of homogeneity
# BodyShapeI 0.08, Troph 0.89, DepthMin 0.50, ParentalCare 0.44, DorsalSpinesMean 0.45
# Dipsersion group p-value is 0.057, M-S p-value is low at 0.046 but others are at 0.05
# PERMANOVA R2 0.65, p-value 0.001, all pairwise are significant with R2 M-S 0.60, S-O 0.59, and M-O 0.29
# mean_traits <- c("BodyShapeI", "Troph", "DepthMin", "ParentalCare", "DorsalSpinesMean")

# When all traits are converted to numerical using factors use this file from dbFD_results.md
# FD_total_n

mean_traits <- c("BodyShapeI", "OperculumPresent", "MaxLengthTL", "Troph", "DepthMin", "DepthMax", "TempPrefMin", "FeedingPath", "ParentalCare", "WaterPref", "DorsalSpinesMean")

# All sites
fd_dist <- gowdis(FD_total_env[surveyed_sites,mean_traits])
# Mixed and stratified lakes
fdms_dist <- gowdis(FD_total_env[mixed_stratified_lakes,mean_traits])
# Ocean sites and mixed lakes
fdom_dist <- gowdis(FD_total_env[ocean_mixed_sites,mean_traits])
# Stratified lakes and ocean sites
fdso_dist <- gowdis(FD_total_env[ocean_stratified_sites,mean_traits])
# Stratified lakes
fds_dist <- gowdis(FD_total_env[stratified_lakes,mean_traits])


### Environmental
# All sites
fde_dist <- gowdis(FD_total_env[surveyed_sites_env,mean_traits])
# Mixed and stratified lakes
fdems_dist <- gowdis(FD_total_env[mixed_stratified_lakes,mean_traits])
# Ocean sites and mixed lakes
fdeom_dist <- gowdis(FD_total_env[ocean_mixed_sites_env,mean_traits])
# Stratified lakes and ocean sites
fdeso_dist <- gowdis(FD_total_env[ocean_stratified_sites_env,mean_traits])


### Geographic
# All sites
fdb_dist <- gowdis(FD_total_env[surveyed_sites_geo,mean_traits])
# Mixed and stratified lakes
fdbms_dist <- gowdis(FD_total_env[mixed_stratified_lakes_geo,mean_traits])
# Ocean sites and mixed lakes
fdbom_dist <- gowdis(FD_total_env[ocean_mixed_sites_geo,mean_traits])
# Stratified lakes and ocean sites
fdbso_dist <- gowdis(FD_total_env[ocean_stratified_sites,mean_traits])
```

## Trait NMDS

- To constrain dissimilarities we perform Nonmetric Multidimensional
  Scaling (NMDS), which tries to find a stable solution using the
  metaMDS package.

``` r
### Regular
# All sites 
fd_NMDS <- metaMDS(fd_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
round(fd_NMDS$stress, digits = 2)
# Mixed and stratified lakes
fdms_NMDS <- metaMDS(fdms_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
fdom_NMDS <- metaMDS(fdom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
fdso_NMDS <- metaMDS(fdso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes
fds_NMDS <- metaMDS(fds_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(fds_dist, try = 1000, parallel = 4, previous.best, trymax =
    ## 500, : stress is (nearly) zero: you may have insufficient data

``` r
### Environmental
# All sites 
fde_NMDS <- metaMDS(fde_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
fdems_NMDS <- metaMDS(fdems_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
fdeom_NMDS <- metaMDS(fdeom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
fdeso_NMDS <- metaMDS(fdeso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)


### Geographic
# All sites 
fdb_NMDS <- metaMDS(fdb_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
fdbms_NMDS <- metaMDS(fdbms_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
fdbom_NMDS <- metaMDS(fdbom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
fdbso_NMDS <- metaMDS(fdbso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
```

## Determine stratification homogeneity using traits

- We use betadisper to determine homogeneity of the trait data based on
  the stratification category. Is the dispersion of traits similar
  between the different stratifications.

``` r
### Stratification
stratification_fd_bd <- betadisper(fd_dist, stratification_group)
```

    ## Warning in betadisper(fd_dist, stratification_group): some squared distances
    ## are negative and changed to zero

``` r
anova(stratification_fd_bd)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq  Mean Sq F value    Pr(>F)    
    ## Groups     2 0.226661 0.113331  23.263 7.802e-06 ***
    ## Residuals 19 0.092562 0.004872                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
permutest(stratification_fd_bd, pairwise = TRUE)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)    
    ## Groups     2 0.226661 0.113331 23.263    999  0.001 ***
    ## Residuals 19 0.092562 0.004872                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##                 Mixed      Ocean Stratified
    ## Mixed                 1.5000e-01      0.001
    ## Ocean      1.4203e-01                 0.001
    ## Stratified 6.6031e-05 1.1054e-03

``` r
(stratification_fd_bd.HSD <- TukeyHSD(stratification_fd_bd))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                        diff         lwr       upr     p adj
    ## Ocean-Mixed      0.02639458 -0.06936777 0.1221569 0.7662467
    ## Stratified-Mixed 0.22120129  0.13254259 0.3098600 0.0000127
    ## Stratified-Ocean 0.19480671  0.09904437 0.2905691 0.0001550

``` r
boxplot(stratification_fd_bd)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Determine%20stratification%20homogeneity%20using%20traits-1.png)<!-- -->

## Determine stratification dissimilarity using traits

- We use adonis to determine if dissimilarities of species traits by
  stratification categories are significant and how much of the
  variation is explained by trait dissimilarities.

``` r
### Stratification
# All sites 
fd_pms <- adonis2(fd_dist ~ env[surveyed_sites,19], permutations = 999)
fd_pms
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = fd_dist ~ env[surveyed_sites, 19], permutations = 999)
    ##                         Df SumOfSqs      R2     F Pr(>F)    
    ## env[surveyed_sites, 19]  2  0.75364 0.50667 9.757  0.001 ***
    ## Residual                19  0.73379 0.49333                 
    ## Total                   21  1.48744 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratifcation
fd_pms_pair <- pairwise.adonis(fd_dist, env[surveyed_sites,19], p.adjust.m = "bonferroni", perm = 999)
fd_pms_pair
```

    ##                 pairs Df  SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 0.54921091 11.097093 0.4421665   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 0.52945100  8.994373 0.4284183   0.003      0.009   *
    ## 3      Mixed vs Ocean  1 0.02717293  4.771947 0.2845196   0.020      0.060

## Envfit environmental influence on traits

- We use envfit to determine significantly correlated environmental
  variables to our trait NMDS. We will use these results, if
  significant, in our figure.

``` r
### Environmental
# All sites for figure
fde_ef <- envfit(fde_NMDS, env[surveyed_sites_env,c(34)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_env,19])
fde_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2   Pr(>r)   
    ## [1,]  0.94663 -0.32231 0.9136 0.005994 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fde_efp <- p.adjust.envfit(fde_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fde_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2   Pr(>r)   
    ## [1,]  0.94663 -0.32231 0.9136 0.005994 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
fdea_ef <- envfit(fde_NMDS, env[surveyed_sites_env,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_env,19])
fdea_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.87420  0.48557 0.0534 0.846154   
    ## salinity_median     0.94663 -0.32231 0.9136 0.003996 **
    ## oxygen_median       0.73955 -0.67310 0.2318 0.998002   
    ## pH_median           0.83888 -0.54431 0.6549 0.107892   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fdea_efp <- p.adjust.envfit(fdea_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fdea_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.87420  0.48557 0.0534 1.00000  
    ## salinity_median     0.94663 -0.32231 0.9136 0.01598 *
    ## oxygen_median       0.73955 -0.67310 0.2318 1.00000  
    ## pH_median           0.83888 -0.54431 0.6549 0.43157  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
fdems_ef <- envfit(fdems_NMDS, env[mixed_stratified_lakes,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
fdems_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.79946  0.60072 0.0773 0.797203   
    ## salinity_median     0.92518 -0.37953 0.9147 0.003996 **
    ## oxygen_median       0.79137 -0.61134 0.1118 1.000000   
    ## pH_median           0.81021 -0.58613 0.5713 0.129870   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fdems_efp <- p.adjust.envfit(fdems_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fdems_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.79946  0.60072 0.0773 1.00000  
    ## salinity_median     0.92518 -0.37953 0.9147 0.01598 *
    ## oxygen_median       0.79137 -0.61134 0.1118 1.00000  
    ## pH_median           0.81021 -0.58613 0.5713 0.51948  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
fdeom_ef <- envfit(fdeom_NMDS, env[ocean_mixed_sites_env,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
fdeom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median  0.24373  0.96984 0.2407 0.30470  
    ## salinity_median    -0.89276  0.45053 0.7838 0.01199 *
    ## oxygen_median      -0.39572  0.91837 0.5589 0.09790 .
    ## pH_median          -0.31837  0.94797 0.4130 0.12388  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fdeom_efp <- p.adjust.envfit(fdeom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fdeom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median  0.24373  0.96984 0.2407 1.00000  
    ## salinity_median    -0.89276  0.45053 0.7838 0.04795 *
    ## oxygen_median      -0.39572  0.91837 0.5589 0.39161  
    ## pH_median          -0.31837  0.94797 0.4130 0.49550  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
fdeso_ef <- envfit(fdeso_NMDS, env[ocean_stratified_sites_env,c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
fdeso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)   
    ## temperature_median  0.74276  0.66956 0.0055 0.97902   
    ## salinity_median    -0.91939 -0.39334 0.8740 0.00999 **
    ## oxygen_median      -0.67775 -0.73529 0.1214 1.00000   
    ## pH_median          -0.79461 -0.60712 0.6924 0.13786   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fdeso_efp <- p.adjust.envfit(fdeso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fdeso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median  0.74276  0.66956 0.0055 1.00000  
    ## salinity_median    -0.91939 -0.39334 0.8740 0.03996 *
    ## oxygen_median      -0.67775 -0.73529 0.1214 1.00000  
    ## pH_median          -0.79461 -0.60712 0.6924 0.55145  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes
fds_ef <- envfit(fds_NMDS, env[stratified_lakes,c(2,6,8,10,22:32)], permutations = 1000, na.rm = TRUE)
fds_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median         -0.28922 -0.95726 0.0255 0.9441
    ## salinity_median             0.45444 -0.89078 0.5940 0.1099
    ## oxygen_median              -0.18095  0.98349 0.5485 0.1548
    ## pH_median                   0.20419 -0.97893 0.3330 0.3856
    ## volume_m3_w_chemocline      0.09370 -0.99560 0.2251 0.4865
    ## volume_m3                  -0.19482 -0.98084 0.3103 0.4386
    ## surface_area_m2            -0.31052 -0.95057 0.1615 0.6843
    ## distance_to_ocean_min_m    -0.99997 -0.00746 0.2102 0.5844
    ## distance_to_ocean_mean_m   -0.59673 -0.80244 0.1947 0.5924
    ## distance_to_ocean_median_m -0.68427 -0.72923 0.1877 0.6114
    ## tidal_lag_time_minutes      0.02960 -0.99956 0.4907 0.1868
    ## tidal_efficiency            0.25521  0.96689 0.0245 0.9860
    ## perimeter_fromSat          -0.97948 -0.20156 0.0411 0.8921
    ## max_depth                  -0.14196 -0.98987 0.2154 0.5634
    ## logArea                    -0.41166 -0.91134 0.0899 0.8212
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fds_efp <- p.adjust.envfit(fds_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fds_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median         -0.28922 -0.95726 0.0255      1
    ## salinity_median             0.45444 -0.89078 0.5940      1
    ## oxygen_median              -0.18095  0.98349 0.5485      1
    ## pH_median                   0.20419 -0.97893 0.3330      1
    ## volume_m3_w_chemocline      0.09370 -0.99560 0.2251      1
    ## volume_m3                  -0.19482 -0.98084 0.3103      1
    ## surface_area_m2            -0.31052 -0.95057 0.1615      1
    ## distance_to_ocean_min_m    -0.99997 -0.00746 0.2102      1
    ## distance_to_ocean_mean_m   -0.59673 -0.80244 0.1947      1
    ## distance_to_ocean_median_m -0.68427 -0.72923 0.1877      1
    ## tidal_lag_time_minutes      0.02960 -0.99956 0.4907      1
    ## tidal_efficiency            0.25521  0.96689 0.0245      1
    ## perimeter_fromSat          -0.97948 -0.20156 0.0411      1
    ## max_depth                  -0.14196 -0.98987 0.2154      1
    ## logArea                    -0.41166 -0.91134 0.0899      1
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Geographic
# All sites for figure
fdb_ef <- envfit(fdb_NMDS, env[surveyed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdb_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.92304 -0.38469 0.0551 0.6204
    ## volume_m3                   0.89185 -0.45234 0.0474 0.6464
    ## surface_area_m2             0.77702 -0.62948 0.0769 0.6903
    ## distance_to_ocean_min_m    -0.90418  0.42715 0.7038 0.1219
    ## distance_to_ocean_mean_m   -0.84975  0.52719 0.5918 0.3217
    ## distance_to_ocean_median_m -0.85295  0.52199 0.5953 0.3097
    ## tidal_lag_time_minutes     -0.69796  0.71613 0.4736 0.8412
    ## tidal_efficiency            0.66946 -0.74285 0.5526 0.2547
    ## perimeter_fromSat           0.76971 -0.63840 0.0615 0.7752
    ## max_depth                  -0.96885  0.24764 0.1087 0.8322
    ## logArea                     0.30262 -0.95311 0.0161 0.7203
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fdb_efp <- p.adjust.envfit(fdb_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fdb_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.92304 -0.38469 0.0551      1
    ## volume_m3                   0.89185 -0.45234 0.0474      1
    ## surface_area_m2             0.77702 -0.62948 0.0769      1
    ## distance_to_ocean_min_m    -0.90418  0.42715 0.7038      1
    ## distance_to_ocean_mean_m   -0.84975  0.52719 0.5918      1
    ## distance_to_ocean_median_m -0.85295  0.52199 0.5953      1
    ## tidal_lag_time_minutes     -0.69796  0.71613 0.4736      1
    ## tidal_efficiency            0.66946 -0.74285 0.5526      1
    ## perimeter_fromSat           0.76971 -0.63840 0.0615      1
    ## max_depth                  -0.96885  0.24764 0.1087      1
    ## logArea                     0.30262 -0.95311 0.0161      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
fdba_ef <- envfit(fdb_NMDS, env[surveyed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.92304 -0.38469 0.0551 0.61039  
    ## volume_m3                   0.89185 -0.45234 0.0474 0.62737  
    ## surface_area_m2             0.77702 -0.62948 0.0769 0.65834  
    ## distance_to_ocean_min_m    -0.90418  0.42715 0.7038 0.08991 .
    ## distance_to_ocean_mean_m   -0.84975  0.52719 0.5918 0.31069  
    ## distance_to_ocean_median_m -0.85295  0.52199 0.5953 0.29670  
    ## tidal_lag_time_minutes     -0.69796  0.71613 0.4736 0.83916  
    ## tidal_efficiency            0.66946 -0.74285 0.5526 0.24775  
    ## perimeter_fromSat           0.76971 -0.63840 0.0615 0.75125  
    ## max_depth                  -0.96885  0.24764 0.1087 0.83616  
    ## logArea                     0.30262 -0.95311 0.0161 0.71628  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fdba_efp <- p.adjust.envfit(fdba_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fdba_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.92304 -0.38469 0.0551  1.000
    ## volume_m3                   0.89185 -0.45234 0.0474  1.000
    ## surface_area_m2             0.77702 -0.62948 0.0769  1.000
    ## distance_to_ocean_min_m    -0.90418  0.42715 0.7038  0.989
    ## distance_to_ocean_mean_m   -0.84975  0.52719 0.5918  1.000
    ## distance_to_ocean_median_m -0.85295  0.52199 0.5953  1.000
    ## tidal_lag_time_minutes     -0.69796  0.71613 0.4736  1.000
    ## tidal_efficiency            0.66946 -0.74285 0.5526  1.000
    ## perimeter_fromSat           0.76971 -0.63840 0.0615  1.000
    ## max_depth                  -0.96885  0.24764 0.1087  1.000
    ## logArea                     0.30262 -0.95311 0.0161  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
fdbms_ef <- envfit(fdbms_NMDS, env[mixed_stratified_lakes_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline     -0.64291  0.76594 0.0050 0.97602  
    ## volume_m3                  -0.92559 -0.37854 0.1644 0.73327  
    ## surface_area_m2            -0.99982  0.01911 0.1561 0.71029  
    ## distance_to_ocean_min_m    -0.92460  0.38093 0.6444 0.09091 .
    ## distance_to_ocean_mean_m   -0.90643  0.42235 0.5078 0.28671  
    ## distance_to_ocean_median_m -0.90864  0.41758 0.5143 0.27672  
    ## tidal_lag_time_minutes     -0.68846  0.72528 0.2897 0.81119  
    ## tidal_efficiency            0.60971 -0.79263 0.3965 0.24675  
    ## perimeter_fromSat          -0.72653  0.68714 0.0981 0.67632  
    ## max_depth                  -0.99980 -0.01979 0.1317 0.77123  
    ## logArea                    -0.94965  0.31330 0.0728 0.78921  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fdbms_efp <- p.adjust.envfit(fdbms_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fdbms_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.64291  0.76594 0.0050      1
    ## volume_m3                  -0.92559 -0.37854 0.1644      1
    ## surface_area_m2            -0.99982  0.01911 0.1561      1
    ## distance_to_ocean_min_m    -0.92460  0.38093 0.6444      1
    ## distance_to_ocean_mean_m   -0.90643  0.42235 0.5078      1
    ## distance_to_ocean_median_m -0.90864  0.41758 0.5143      1
    ## tidal_lag_time_minutes     -0.68846  0.72528 0.2897      1
    ## tidal_efficiency            0.60971 -0.79263 0.3965      1
    ## perimeter_fromSat          -0.72653  0.68714 0.0981      1
    ## max_depth                  -0.99980 -0.01979 0.1317      1
    ## logArea                    -0.94965  0.31330 0.0728      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
fdbom_ef <- envfit(fdbom_NMDS, env[ocean_mixed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.06230  0.99806 0.1381 0.49850  
    ## volume_m3                   0.06102  0.99814 0.1386 0.49750  
    ## surface_area_m2             0.30688  0.95175 0.0856 0.70829  
    ## distance_to_ocean_min_m    -0.87823  0.47823 0.4562 0.13087  
    ## distance_to_ocean_mean_m   -0.83648  0.54799 0.4397 0.14985  
    ## distance_to_ocean_median_m -0.83522  0.54991 0.4273 0.16484  
    ## tidal_lag_time_minutes     -0.80285  0.59618 0.5449 0.04096 *
    ## tidal_efficiency            0.75576 -0.65485 0.4689 0.10689  
    ## perimeter_fromSat           0.19669  0.98047 0.1130 0.67233  
    ## max_depth                   0.16667  0.98601 0.2352 0.20280  
    ## logArea                     0.26008  0.96559 0.1859 0.44356  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fdbom_efp <- p.adjust.envfit(fdbom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fdbom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.06230  0.99806 0.1381 1.0000
    ## volume_m3                   0.06102  0.99814 0.1386 1.0000
    ## surface_area_m2             0.30688  0.95175 0.0856 1.0000
    ## distance_to_ocean_min_m    -0.87823  0.47823 0.4562 1.0000
    ## distance_to_ocean_mean_m   -0.83648  0.54799 0.4397 1.0000
    ## distance_to_ocean_median_m -0.83522  0.54991 0.4273 1.0000
    ## tidal_lag_time_minutes     -0.80285  0.59618 0.5449 0.4505
    ## tidal_efficiency            0.75576 -0.65485 0.4689 1.0000
    ## perimeter_fromSat           0.19669  0.98047 0.1130 1.0000
    ## max_depth                   0.16667  0.98601 0.2352 1.0000
    ## logArea                     0.26008  0.96559 0.1859 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
fdbso_ef <- envfit(fdbso_NMDS, env[ocean_stratified_sites,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline     -0.97186 -0.23554 0.1311 0.63536  
    ## volume_m3                  -0.96137 -0.27525 0.1206 0.68232  
    ## surface_area_m2            -0.91371 -0.40636 0.1947 0.61039  
    ## distance_to_ocean_min_m     0.93664  0.35028 0.7728 0.07892 .
    ## distance_to_ocean_mean_m    0.90284  0.42998 0.6573 0.29171  
    ## distance_to_ocean_median_m  0.90500  0.42541 0.6612 0.27872  
    ## tidal_lag_time_minutes      0.80055  0.59927 0.5558 0.92208  
    ## tidal_efficiency           -0.77188 -0.63577 0.7057 0.34066  
    ## perimeter_fromSat          -0.90588 -0.42354 0.1867 0.73626  
    ## max_depth                   0.88644  0.46284 0.0853 0.82418  
    ## logArea                    -0.77288 -0.63456 0.1003 0.66633  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fdbso_efp <- p.adjust.envfit(fdbso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fdbso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.97186 -0.23554 0.1311 1.0000
    ## volume_m3                  -0.96137 -0.27525 0.1206 1.0000
    ## surface_area_m2            -0.91371 -0.40636 0.1947 1.0000
    ## distance_to_ocean_min_m     0.93664  0.35028 0.7728 0.8681
    ## distance_to_ocean_mean_m    0.90284  0.42998 0.6573 1.0000
    ## distance_to_ocean_median_m  0.90500  0.42541 0.6612 1.0000
    ## tidal_lag_time_minutes      0.80055  0.59927 0.5558 1.0000
    ## tidal_efficiency           -0.77188 -0.63577 0.7057 1.0000
    ## perimeter_fromSat          -0.90588 -0.42354 0.1867 1.0000
    ## max_depth                   0.88644  0.46284 0.0853 1.0000
    ## logArea                    -0.77288 -0.63456 0.1003 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Species
fds_ef <- envfit(fd_NMDS, surveyed_sites_lake, permutations = 1000, na.rm = TRUE)
fds_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                                    NMDS1    NMDS2     r2   Pr(>r)   
    ## Abudefduf_lorenzi               -0.10992 -0.99394 0.0977 0.272727   
    ## Abudefduf_septemfasciatus        0.99507 -0.09914 0.0189 0.782218   
    ## Abudefduf_sexfasciatus          -0.10992 -0.99394 0.0977 0.272727   
    ## Acanthurus_lineatus              0.81508 -0.57934 0.0282 0.654346   
    ## Acanthurus_nigricauda            0.89704 -0.44195 0.0698 0.479520   
    ## Acanthurus_sp                    0.60930 -0.79294 0.0369 0.453546   
    ## Acanthurus_xanthopterus          0.87335  0.48710 0.2458 0.059940 . 
    ## Acentrogobius_janthinopterus    -0.69949  0.71465 0.5665 0.002997 **
    ## Amblyeleotris_gymnocephala       0.62285 -0.78234 0.0230 0.738262   
    ## Amblyglyphidodon_curacao         0.30794 -0.95140 0.1771 0.153846   
    ## Amblygobius_buanensis            0.97524 -0.22117 0.3733 0.006993 **
    ## Amblygobius_decussatus           0.47726 -0.87876 0.0849 0.345654   
    ## Amblygobius_esakiae              0.60930 -0.79294 0.0369 0.453546   
    ## Amblygobius_linki                0.99677  0.08032 0.0184 0.826174   
    ## Amblygobius_phalaena             0.78723 -0.61666 0.0343 0.854146   
    ## Amphiprion_clarkii               0.97551 -0.21996 0.0291 0.558442   
    ## Ancistrogobius_dipus             0.64556  0.76371 0.0286 0.565435   
    ## Apogonichthyoides_melas          0.64556  0.76371 0.0286 0.565435   
    ## Arcygobius_baliurus              0.99496 -0.10025 0.0160 0.916084   
    ## Arothron_reticularis             0.63276  0.77434 0.0173 0.864136   
    ## Asterropteryx_semipunctata       0.85047 -0.52603 0.0633 0.543457   
    ## Asterropteryx_sp                 0.99677  0.08032 0.0184 0.826174   
    ## Atherinomorus_endrachtensis      0.60712 -0.79461 0.1136 0.327672   
    ## Atherinomorus_sp                 0.84878 -0.52874 0.0313 0.757243   
    ## Atule_mate                       0.81508 -0.57934 0.0282 0.654346   
    ## Balistoides_viridescens          0.97754 -0.21073 0.1518 0.185814   
    ## Bolbometopon_muricatum           0.58325 -0.81229 0.0901 0.337662   
    ## Caesio_caerulaurea               0.99620  0.08712 0.0414 0.708292   
    ## Caesio_cuning                    0.56489 -0.82517 0.1710 0.175824   
    ## Caesio_teres                     0.62285 -0.78234 0.0230 0.738262   
    ## Cantherhines_dumerilii           0.97551 -0.21996 0.0291 0.558442   
    ## Canthigaster_bennetti            0.97551 -0.21996 0.0291 0.558442   
    ## Canthigaster_solandri            0.99507 -0.09914 0.0189 0.782218   
    ## Caranx_ignobilis                 0.97015 -0.24250 0.0123 0.943057   
    ## Caranx_melampygus                0.51541 -0.85695 0.0123 0.903097   
    ## Caranx_papuensis                 0.62285 -0.78234 0.0230 0.738262   
    ## Caranx_sexfasciatus              0.66234  0.74921 0.1634 0.166833   
    ## Cephalopholis_argus             -0.44019 -0.89790 0.3364 0.133866   
    ## Cephalopholis_boenak             0.80369 -0.59505 0.0647 0.421578   
    ## Chaetodon_auriga                 0.91053 -0.41344 0.3545 0.013986 * 
    ## Chaetodon_bennetti               0.90905 -0.41668 0.0592 0.469530   
    ## Chaetodon_ephippium              0.69774 -0.71635 0.3357 0.016983 * 
    ## Chaetodon_kleinii                0.82742 -0.56159 0.0851 0.344655   
    ## Chaetodon_lunula                 0.96902 -0.24697 0.1702 0.157842   
    ## Chaetodon_lunulatus              0.58047 -0.81428 0.2239 0.084915 . 
    ## Chaetodon_melannotus             0.69970 -0.71443 0.0258 0.668332   
    ## Chaetodon_ocellicaudus           0.97015 -0.24250 0.0123 0.943057   
    ## Chaetodon_octofasciatus          0.80369 -0.59505 0.0647 0.421578   
    ## Chaetodon_rafflesii              0.84181 -0.53977 0.1148 0.272727   
    ## Chaetodon_semeion                0.86314 -0.50496 0.0554 0.484515   
    ## Chaetodon_trifascialis           0.97551 -0.21996 0.0291 0.558442   
    ## Chaetodon_ulietensis             0.61123 -0.79145 0.1247 0.241758   
    ## Chaetodon_vagabundus             0.85994 -0.51039 0.3542 0.005994 **
    ## Chaetodontoplus_poliourus        0.97551 -0.21996 0.0291 0.558442   
    ## Cheilinus_fasciatus              0.80739 -0.59002 0.0993 0.296703   
    ## Cheilinus_trilobatus             0.81508 -0.57934 0.0282 0.654346   
    ## Cheilinus_undulatus              0.40084 -0.91615 0.0676 0.368631   
    ## Cheilodipterus_artus             0.18546 -0.98265 0.0782 0.499500   
    ## Cheilodipterus_isostigma         0.43580 -0.90004 0.0691 0.370629   
    ## Cheilodipterus_quinquelineatus   0.77642 -0.63021 0.3290 0.015984 * 
    ## Cheilodipterus_singapurensis     0.77208 -0.63553 0.0348 0.864136   
    ## Chlorurus_bleekeri               0.58087 -0.81400 0.1647 0.145854   
    ## Chlorurus_microrhinos            0.97345 -0.22889 0.0416 0.689311   
    ## Chlorurus_spilurus               0.84181 -0.53977 0.1148 0.272727   
    ## Choerodon_anchorago              0.92604 -0.37743 0.4596 0.001998 **
    ## Chromis_viridis                  0.48637 -0.87375 0.0962 0.280719   
    ## Chrysiptera_biocellata           0.99502 -0.09967 0.0366 0.821179   
    ## Chrysiptera_oxycephala           0.30534 -0.95224 0.2330 0.082917 . 
    ## Corythoichthys_ocellatus         0.97551 -0.21996 0.0291 0.558442   
    ## Cryptocentrus_caeruleomaculatus  0.85768  0.51419 0.0763 0.515485   
    ## Cryptocentrus_cyanospilotus      0.07760  0.99698 0.0817 0.337662   
    ## Cryptocentrus_leptocephalus      0.98583 -0.16776 0.0296 0.927073   
    ## Cryptocentrus_strigilliceps      0.59206 -0.80590 0.1206 0.242757   
    ## Ctenochaetus_striatus            0.19139 -0.98151 0.2053 0.086913 . 
    ## Ctenogobiops_pomastictus         0.64556  0.76371 0.0286 0.565435   
    ## Cymbacephalus_beauforti          0.97015 -0.24250 0.0123 0.943057   
    ## Dascyllus_aruanus                0.75557 -0.65507 0.0968 0.354645   
    ## Dascyllus_trimaculatus           0.97551 -0.21996 0.0291 0.558442   
    ## Diplogrammus_goramensis          0.87425  0.48547 0.0363 0.826174   
    ## Diproctacanthus_xanthurus        0.90905 -0.41668 0.0592 0.469530   
    ## Dischistodus_chrysopoecilus      0.80632  0.59148 0.0369 0.802198   
    ## Dischistodus_perspicillatus      0.77943 -0.62649 0.1554 0.183816   
    ## Eleotris_fusca                  -0.29588  0.95523 0.5192 0.014985 * 
    ## Epibulus_brevis                  0.58785 -0.80897 0.2014 0.098901 . 
    ## Epinephelus_coeruleopunctatus    0.99998 -0.00701 0.0360 0.853147   
    ## Epinephelus_merra                0.49995 -0.86605 0.1688 0.124875   
    ## Epinephelus_sp                   0.54640  0.83752 0.0318 0.490509   
    ## Eviota_atriventris               0.99806 -0.06226 0.0696 0.461538   
    ## Eviota_bifasciata                0.71997 -0.69401 0.2562 0.052947 . 
    ## Eviota_fallax                    0.60437 -0.79670 0.2469 0.059940 . 
    ## Eviota_lachdeberei               0.95552 -0.29491 0.0796 0.522478   
    ## Eviota_maculosa                  0.22647 -0.97402 0.0837 0.262737   
    ## Eviota_sigillata                 0.66215 -0.74937 0.0511 0.527473   
    ## Eviota_sp                        0.99677  0.08032 0.0184 0.826174   
    ## Eviota_storthynx                 0.89946 -0.43700 0.0559 0.645355   
    ## Exyrias_belissimus               0.80190 -0.59745 0.1812 0.143856   
    ## Exyrias_puntang                 -0.07025  0.99753 0.3998 0.006993 **
    ## Favonigobius_reichei             0.98583 -0.16776 0.0296 0.927073   
    ## Fibramia_lateralis               0.04521  0.99898 0.0111 0.894106   
    ## Fibramia_thermalis               0.93933 -0.34301 0.0919 0.446553   
    ## Fusigobius_signipinnis           0.87645 -0.48148 0.0656 0.507493   
    ## Gerres_erythrourus               0.97015 -0.24250 0.0123 0.943057   
    ## Gerres_oyena                     0.98642 -0.16424 0.0324 0.890110   
    ## Gladiogobius_ensifer             0.22647 -0.97402 0.0837 0.262737   
    ## Gnathanodon_speciosus            0.25895  0.96589 0.0915 0.303696   
    ## Gobiodon_okinawae                0.81508 -0.57934 0.0282 0.654346   
    ## Gymnothorax_pictus               0.99502 -0.09967 0.0366 0.821179   
    ## Halichoeres_chloropterus         0.42901 -0.90330 0.1679 0.185814   
    ## Halichoeres_leucurus             0.76623 -0.64257 0.1301 0.213786   
    ## Halichoeres_marginatus           0.97015 -0.24250 0.0123 0.943057   
    ## Halichoeres_melanurus            0.16567 -0.98618 0.1794 0.140859   
    ## Halichoeres_papilionaceus        0.93319 -0.35938 0.0594 0.565435   
    ## Halichoeres_trimaculatus         0.98642 -0.16424 0.0324 0.890110   
    ## Hemiglyphidodon_plagiometopon    0.37499 -0.92703 0.2710 0.052947 . 
    ## Hemigymnus_melapterus            0.97345 -0.22889 0.0416 0.689311   
    ## Heniochus_chrysostomus           0.92665 -0.37593 0.0734 0.429570   
    ## Heniochus_monoceros              0.89946 -0.43700 0.0559 0.645355   
    ## Herklotsichthys_quadrimaculatus  0.97015 -0.24250 0.0123 0.943057   
    ## Hipposcarus_longiceps            0.59050 -0.80704 0.1375 0.203796   
    ## Istigobius_decoratus             0.97015 -0.24250 0.0123 0.943057   
    ## Koumansetta_rainfordi            0.97551 -0.21996 0.0291 0.558442   
    ## Kyphosus_cinerascens            -0.12570 -0.99207 0.1021 0.277722   
    ## Lethrinus_erythropterus          0.99831 -0.05809 0.0735 0.518482   
    ## Lethrinus_harak                  0.92136  0.38871 0.1425 0.249750   
    ## Lethrinus_obsoletus              0.64556  0.76371 0.0286 0.565435   
    ## Lethrinus_olivaceus              0.64556  0.76371 0.0286 0.565435   
    ## Lutjanus_argentimaculatus        0.59956  0.80033 0.0950 0.434565   
    ## Lutjanus_biguttatus             -0.13360 -0.99104 0.0695 0.507493   
    ## Lutjanus_bohar                   0.86314 -0.50496 0.0554 0.484515   
    ## Lutjanus_decussatus              0.64556  0.76371 0.0286 0.565435   
    ## Lutjanus_ehrenbergii             0.88272  0.46991 0.1130 0.339660   
    ## Lutjanus_fulviflamma             0.98042 -0.19690 0.0795 0.515485   
    ## Lutjanus_fulvus                  0.74252  0.66982 0.3371 0.017982 * 
    ## Lutjanus_gibbus                  0.99732 -0.07323 0.0916 0.396603   
    ## Lutjanus_monostigma              0.90064  0.43457 0.0269 0.950050   
    ## Lutjanus_russellii               0.99677  0.08032 0.0184 0.826174   
    ## Macrodontogobius_wilburi         0.58047 -0.81428 0.2239 0.084915 . 
    ## Megalops_cyprinoides            -0.21977  0.97555 0.4589 0.048951 * 
    ## Meiacanthus_ditrema              0.60930 -0.79294 0.0369 0.453546   
    ## Meiacanthus_grammistes           0.80369 -0.59505 0.0647 0.421578   
    ## Monodactylus_argenteus           0.50255  0.86455 0.0575 0.606394   
    ## Monotaxis_grandoculis            0.05563 -0.99845 0.0424 0.714286   
    ## Mugilogobius_cavifrons          -0.95857  0.28486 0.1926 0.196803   
    ## Mulloidichthys_flavolineatus     0.97015 -0.24250 0.0123 0.943057   
    ## Myripristis_adusta               0.97925 -0.20266 0.1337 0.265734   
    ## Myripristis_sp                   0.99496 -0.10025 0.0160 0.916084   
    ## Myripristis_violacea             0.97551 -0.21996 0.0291 0.558442   
    ## Naso_brevirostris                0.98962 -0.14373 0.1029 0.323676   
    ## Naso_vlamingii                   0.95166 -0.30714 0.0779 0.402597   
    ## Nectamia_viria                   0.69970 -0.71443 0.0258 0.668332   
    ## Neoglyphidodon_melas            -0.44019 -0.89790 0.3364 0.133866   
    ## Neoniphon_argenteus              0.86892 -0.49495 0.0527 0.656344   
    ## Neoniphon_opercularis           -0.44019 -0.89790 0.3364 0.133866   
    ## Neoniphon_sammara                0.94327 -0.33202 0.1108 0.316683   
    ## Neopomacentrus_nemurus           0.28545 -0.95839 0.1664 0.156843   
    ## Omobranchus_obliquus             0.51251  0.85868 0.0627 0.563437   
    ## Oplopomops_diacanthus            0.97015 -0.24250 0.0123 0.943057   
    ## Oxycheilinus_celebicus           0.17006 -0.98543 0.1401 0.202797   
    ## Parapercis_cylindrica            0.97015 -0.24250 0.0123 0.943057   
    ## Parioglossus_formosus           -0.11512  0.99335 0.1649 0.162837   
    ## Parioglossus_philippinus         0.97303 -0.23068 0.0666 0.489510   
    ## Parioglossus_sp                  0.98551 -0.16963 0.0463 0.593407   
    ## Parupeneus_barberinus            0.93356 -0.35841 0.2063 0.097902 . 
    ## Parupeneus_multifasciatus        0.75567 -0.65496 0.2919 0.025974 * 
    ## Pentapodus_trivittatus           0.79542 -0.60606 0.2723 0.039960 * 
    ## Periophthalmus_argentilineatus   0.99507 -0.09914 0.0189 0.782218   
    ## Pervagor_janthinosoma            0.81508 -0.57934 0.0282 0.654346   
    ## Pervagor_nigrolineatus           0.46280 -0.88646 0.0943 0.312687   
    ## Petroscirtes_breviceps           0.84342 -0.53725 0.1306 0.230769   
    ## Platax_orbicularis               0.86641  0.49933 0.0415 0.691309   
    ## Platax_pinnatus                  0.62285 -0.78234 0.0230 0.738262   
    ## Platax_sp                        0.54640  0.83752 0.0318 0.490509   
    ## Plectorhinchus_albovittatus      0.99677  0.08032 0.0184 0.826174   
    ## Plectorhinchus_chaetodonoides    0.22647 -0.97402 0.0837 0.262737   
    ## Plectorhinchus_gibbosus          0.97015 -0.24250 0.0123 0.943057   
    ## Plectropomus_leopardus           0.97551 -0.21996 0.0291 0.558442   
    ## Pomacanthus_sexstriatus          0.82653 -0.56290 0.0374 0.754246   
    ## Pomacentrus_adelus               0.99507 -0.09914 0.0189 0.782218   
    ## Pomacentrus_amboinensis          0.97551 -0.21996 0.0291 0.558442   
    ## Pomacentrus_burroughi            0.66367 -0.74803 0.2741 0.051948 . 
    ## Pomacentrus_grammorhynchus       0.97015 -0.24250 0.0123 0.943057   
    ## Pomacentrus_nigromanus          -0.09117 -0.99584 0.1122 0.242757   
    ## Pomacentrus_pavo                 0.97551 -0.21996 0.0291 0.558442   
    ## Pomacentrus_simsiang             0.82943 -0.55860 0.3555 0.013986 * 
    ## Pomacentrus_sp                  -0.31964 -0.94754 0.1460 0.186813   
    ## Pomacentrus_taeniometopon        0.99677  0.08032 0.0184 0.826174   
    ## Pseudobalistes_flavimarginatus   0.91882 -0.39467 0.1694 0.146853   
    ## Pseudochromis_fuscus             0.97551 -0.21996 0.0291 0.558442   
    ## Ptereleotris_brachyptera         0.83397 -0.55181 0.0513 0.510490   
    ## Ptereleotris_evides              0.97551 -0.21996 0.0291 0.558442   
    ## Pterocaesio_tile                 0.64556  0.76371 0.0286 0.565435   
    ## Pterocaesio_trilineata           0.69970 -0.71443 0.0258 0.668332   
    ## Pterois_volitans                 0.62285 -0.78234 0.0230 0.738262   
    ## Rastrelliger_kanagurta           0.41320 -0.91064 0.0943 0.287712   
    ## Rhabdamia_gracilis               0.88734  0.46111 0.0849 0.458541   
    ## Rhinecanthus_aculeatus           0.94141 -0.33725 0.0782 0.370629   
    ## Salarias_segmentatus             0.89816 -0.43967 0.2200 0.072927 . 
    ## Sargocentron_diadema             0.75923 -0.65082 0.0564 0.454545   
    ## Sargocentron_spiniferum          0.99496 -0.10025 0.0160 0.916084   
    ## Saurida_gracilis                 0.81508 -0.57934 0.0282 0.654346   
    ## Scarus_dimidiatus                0.84181 -0.53977 0.1148 0.272727   
    ## Scarus_flavipectoralis           0.90905 -0.41668 0.0592 0.469530   
    ## Scarus_ghobban                   0.98551 -0.16963 0.0463 0.593407   
    ## Scarus_hypselopterus             0.81508 -0.57934 0.0282 0.654346   
    ## Scarus_oviceps                   0.60930 -0.79294 0.0369 0.453546   
    ## Scarus_psittacus                 0.85607 -0.51686 0.1350 0.204795   
    ## Scarus_quoyi                     0.97551 -0.21996 0.0291 0.558442   
    ## Scarus_rivulatus                 0.47726 -0.87876 0.0849 0.345654   
    ## Scatophagus_argus                0.63276  0.77434 0.0173 0.864136   
    ## Scolopsis_ciliata                0.78189 -0.62342 0.2106 0.094905 . 
    ## Scolopsis_lineata                0.98583 -0.16776 0.0296 0.927073   
    ## Scolopsis_margaritifera          0.30245 -0.95316 0.2360 0.055944 . 
    ## Scolopsis_trilineata             0.99958  0.02897 0.0611 0.557443   
    ## Siganus_doliatus                 0.58325 -0.81229 0.0901 0.337662   
    ## Siganus_fuscescens               0.99755  0.06989 0.0612 0.553447   
    ## Siganus_lineatus                 0.73021  0.68323 0.2811 0.038961 * 
    ## Siganus_puellus                  0.51858 -0.85503 0.1292 0.217782   
    ## Siganus_punctatissimus           0.97551 -0.21996 0.0291 0.558442   
    ## Siganus_punctatus                0.97551 -0.21996 0.0291 0.558442   
    ## Signigobius_biocellatus          0.81508 -0.57934 0.0282 0.654346   
    ## Silhouettea_sp                   0.07760  0.99698 0.0817 0.337662   
    ## Sphaeramia_nematoptera           0.59124 -0.80650 0.1951 0.119880   
    ## Sphaeramia_orbicularis           0.70741 -0.70680 0.1070 0.383616   
    ## Sphyraena_barracuda              0.88644  0.46285 0.1165 0.306693   
    ## Sphyraena_obtusata               0.88332 -0.46877 0.0097 1.000000   
    ## Sphyraena_qenie                  0.99496 -0.10025 0.0160 0.916084   
    ## Stegastes_nigricans              0.97015 -0.24250 0.0123 0.943057   
    ## Stegastes_punctatus              0.88713 -0.46152 0.0403 0.736264   
    ## Stethojulis_strigiventer         0.98975 -0.14278 0.0518 0.691309   
    ## Strongylura_incisa               0.91608 -0.40099 0.0444 0.627373   
    ## Synodus_variegatus               0.41320 -0.91064 0.0943 0.287712   
    ## Taeniamia_zosterophora           0.46969 -0.88283 0.1098 0.248751   
    ## Taeniurops_meyeni                0.64556  0.76371 0.0286 0.565435   
    ## Toxotes_jaculatrix               0.99567 -0.09299 0.1209 0.285714   
    ## Tylosurus_punctulatus            0.72764 -0.68596 0.0320 0.892108   
    ## Upeneus_taeniopterus             0.99496 -0.10025 0.0160 0.916084   
    ## Valenciennea_muralis             0.99999 -0.00433 0.0565 0.628372   
    ## Valenciennea_parva               0.62285 -0.78234 0.0230 0.738262   
    ## Valenciennea_randalli            0.62285 -0.78234 0.0230 0.738262   
    ## Zanclus_cornutus                 0.87645 -0.48148 0.0656 0.507493   
    ## Zebrasoma_scopas                 0.97551 -0.21996 0.0291 0.558442   
    ## Zebrasoma_velifer                0.64129 -0.76729 0.2697 0.044955 * 
    ## Zenarchopterus_dispar            0.90028 -0.43531 0.3055 0.027972 * 
    ## Zenarchopterus_gilli             0.64801 -0.76163 0.1532 0.225774   
    ## Zoramia_leptacanthus             0.97015 -0.24250 0.0123 0.943057   
    ## Zoramia_viridiventer             0.71256 -0.70161 0.1084 0.278721   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# By stratification
fds_ef <- envfit(fd_NMDS, surveyed_sites_lake, permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites,19])
fds_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                                    NMDS1    NMDS2     r2  Pr(>r)  
    ## Abudefduf_lorenzi               -0.10992 -0.99394 0.0977 0.18781  
    ## Abudefduf_septemfasciatus        0.99507 -0.09914 0.0189 0.63337  
    ## Abudefduf_sexfasciatus          -0.10992 -0.99394 0.0977 0.18781  
    ## Acanthurus_lineatus              0.81508 -0.57934 0.0282 0.82817  
    ## Acanthurus_nigricauda            0.89704 -0.44195 0.0698 0.45355  
    ## Acanthurus_sp                    0.60930 -0.79294 0.0369 0.50150  
    ## Acanthurus_xanthopterus          0.87335  0.48710 0.2458 0.09291 .
    ## Acentrogobius_janthinopterus    -0.69949  0.71465 0.5665 0.11888  
    ## Amblyeleotris_gymnocephala       0.62285 -0.78234 0.0230 0.49351  
    ## Amblyglyphidodon_curacao         0.30794 -0.95140 0.1771 0.23477  
    ## Amblygobius_buanensis            0.97524 -0.22117 0.3733 0.62937  
    ## Amblygobius_decussatus           0.47726 -0.87876 0.0849 0.39860  
    ## Amblygobius_esakiae              0.60930 -0.79294 0.0369 0.50150  
    ## Amblygobius_linki                0.99677  0.08032 0.0184 0.75724  
    ## Amblygobius_phalaena             0.78723 -0.61666 0.0343 0.91708  
    ## Amphiprion_clarkii               0.97551 -0.21996 0.0291 0.67532  
    ## Ancistrogobius_dipus             0.64556  0.76371 0.0286 0.26374  
    ## Apogonichthyoides_melas          0.64556  0.76371 0.0286 0.26374  
    ## Arcygobius_baliurus              0.99496 -0.10025 0.0160 0.86114  
    ## Arothron_reticularis             0.63276  0.77434 0.0173 1.00000  
    ## Asterropteryx_semipunctata       0.85047 -0.52603 0.0633 0.88412  
    ## Asterropteryx_sp                 0.99677  0.08032 0.0184 0.75724  
    ## Atherinomorus_endrachtensis      0.60712 -0.79461 0.1136 0.23177  
    ## Atherinomorus_sp                 0.84878 -0.52874 0.0313 0.69730  
    ## Atule_mate                       0.81508 -0.57934 0.0282 0.82817  
    ## Balistoides_viridescens          0.97754 -0.21073 0.1518 0.79720  
    ## Bolbometopon_muricatum           0.58325 -0.81229 0.0901 0.81219  
    ## Caesio_caerulaurea               0.99620  0.08712 0.0414 0.37962  
    ## Caesio_cuning                    0.56489 -0.82517 0.1710 0.32068  
    ## Caesio_teres                     0.62285 -0.78234 0.0230 0.49351  
    ## Cantherhines_dumerilii           0.97551 -0.21996 0.0291 0.67532  
    ## Canthigaster_bennetti            0.97551 -0.21996 0.0291 0.67532  
    ## Canthigaster_solandri            0.99507 -0.09914 0.0189 0.63337  
    ## Caranx_ignobilis                 0.97015 -0.24250 0.0123 1.00000  
    ## Caranx_melampygus                0.51541 -0.85695 0.0123 0.84316  
    ## Caranx_papuensis                 0.62285 -0.78234 0.0230 0.49351  
    ## Caranx_sexfasciatus              0.66234  0.74921 0.1634 0.04196 *
    ## Cephalopholis_argus             -0.44019 -0.89790 0.3364 0.38861  
    ## Cephalopholis_boenak             0.80369 -0.59505 0.0647 0.65734  
    ## Chaetodon_auriga                 0.91053 -0.41344 0.3545 0.73127  
    ## Chaetodon_bennetti               0.90905 -0.41668 0.0592 0.81219  
    ## Chaetodon_ephippium              0.69774 -0.71635 0.3357 0.33167  
    ## Chaetodon_kleinii                0.82742 -0.56159 0.0851 0.49351  
    ## Chaetodon_lunula                 0.96902 -0.24697 0.1702 0.09690 .
    ## Chaetodon_lunulatus              0.58047 -0.81428 0.2239 0.17882  
    ## Chaetodon_melannotus             0.69970 -0.71443 0.0258 0.37263  
    ## Chaetodon_ocellicaudus           0.97015 -0.24250 0.0123 1.00000  
    ## Chaetodon_octofasciatus          0.80369 -0.59505 0.0647 0.65734  
    ## Chaetodon_rafflesii              0.84181 -0.53977 0.1148 1.00000  
    ## Chaetodon_semeion                0.86314 -0.50496 0.0554 0.34665  
    ## Chaetodon_trifascialis           0.97551 -0.21996 0.0291 0.67532  
    ## Chaetodon_ulietensis             0.61123 -0.79145 0.1247 0.48851  
    ## Chaetodon_vagabundus             0.85994 -0.51039 0.3542 0.36763  
    ## Chaetodontoplus_poliourus        0.97551 -0.21996 0.0291 0.67532  
    ## Cheilinus_fasciatus              0.80739 -0.59002 0.0993 0.64635  
    ## Cheilinus_trilobatus             0.81508 -0.57934 0.0282 0.82817  
    ## Cheilinus_undulatus              0.40084 -0.91615 0.0676 0.55445  
    ## Cheilodipterus_artus             0.18546 -0.98265 0.0782 0.46953  
    ## Cheilodipterus_isostigma         0.43580 -0.90004 0.0691 0.12388  
    ## Cheilodipterus_quinquelineatus   0.77642 -0.63021 0.3290 0.38362  
    ## Cheilodipterus_singapurensis     0.77208 -0.63553 0.0348 0.79520  
    ## Chlorurus_bleekeri               0.58087 -0.81400 0.1647 0.40959  
    ## Chlorurus_microrhinos            0.97345 -0.22889 0.0416 0.93806  
    ## Chlorurus_spilurus               0.84181 -0.53977 0.1148 1.00000  
    ## Choerodon_anchorago              0.92604 -0.37743 0.4596 0.87912  
    ## Chromis_viridis                  0.48637 -0.87375 0.0962 0.31868  
    ## Chrysiptera_biocellata           0.99502 -0.09967 0.0366 0.68032  
    ## Chrysiptera_oxycephala           0.30534 -0.95224 0.2330 0.09990 .
    ## Corythoichthys_ocellatus         0.97551 -0.21996 0.0291 0.67532  
    ## Cryptocentrus_caeruleomaculatus  0.85768  0.51419 0.0763 0.35864  
    ## Cryptocentrus_cyanospilotus      0.07760  0.99698 0.0817 0.74426  
    ## Cryptocentrus_leptocephalus      0.98583 -0.16776 0.0296 0.97502  
    ## Cryptocentrus_strigilliceps      0.59206 -0.80590 0.1206 0.55644  
    ## Ctenochaetus_striatus            0.19139 -0.98151 0.2053 0.21678  
    ## Ctenogobiops_pomastictus         0.64556  0.76371 0.0286 0.26374  
    ## Cymbacephalus_beauforti          0.97015 -0.24250 0.0123 1.00000  
    ## Dascyllus_aruanus                0.75557 -0.65507 0.0968 0.89810  
    ## Dascyllus_trimaculatus           0.97551 -0.21996 0.0291 0.67532  
    ## Diplogrammus_goramensis          0.87425  0.48547 0.0363 0.91009  
    ## Diproctacanthus_xanthurus        0.90905 -0.41668 0.0592 0.81219  
    ## Dischistodus_chrysopoecilus      0.80632  0.59148 0.0369 0.88212  
    ## Dischistodus_perspicillatus      0.77943 -0.62649 0.1554 0.72328  
    ## Eleotris_fusca                  -0.29588  0.95523 0.5192 0.19381  
    ## Epibulus_brevis                  0.58785 -0.80897 0.2014 0.20080  
    ## Epinephelus_coeruleopunctatus    0.99998 -0.00701 0.0360 0.70330  
    ## Epinephelus_merra                0.49995 -0.86605 0.1688 0.10390  
    ## Epinephelus_sp                   0.54640  0.83752 0.0318 0.13886  
    ## Eviota_atriventris               0.99806 -0.06226 0.0696 0.47253  
    ## Eviota_bifasciata                0.71997 -0.69401 0.2562 0.13986  
    ## Eviota_fallax                    0.60437 -0.79670 0.2469 0.05594 .
    ## Eviota_lachdeberei               0.95552 -0.29491 0.0796 0.81019  
    ## Eviota_maculosa                  0.22647 -0.97402 0.0837 0.17083  
    ## Eviota_sigillata                 0.66215 -0.74937 0.0511 0.06893 .
    ## Eviota_sp                        0.99677  0.08032 0.0184 0.75724  
    ## Eviota_storthynx                 0.89946 -0.43700 0.0559 0.86014  
    ## Exyrias_belissimus               0.80190 -0.59745 0.1812 0.82318  
    ## Exyrias_puntang                 -0.07025  0.99753 0.3998 0.04396 *
    ## Favonigobius_reichei             0.98583 -0.16776 0.0296 0.97502  
    ## Fibramia_lateralis               0.04521  0.99898 0.0111 0.84715  
    ## Fibramia_thermalis               0.93933 -0.34301 0.0919 0.32667  
    ## Fusigobius_signipinnis           0.87645 -0.48148 0.0656 0.85315  
    ## Gerres_erythrourus               0.97015 -0.24250 0.0123 1.00000  
    ## Gerres_oyena                     0.98642 -0.16424 0.0324 0.93906  
    ## Gladiogobius_ensifer             0.22647 -0.97402 0.0837 0.17083  
    ## Gnathanodon_speciosus            0.25895  0.96589 0.0915 0.35265  
    ## Gobiodon_okinawae                0.81508 -0.57934 0.0282 0.82817  
    ## Gymnothorax_pictus               0.99502 -0.09967 0.0366 0.68032  
    ## Halichoeres_chloropterus         0.42901 -0.90330 0.1679 0.39660  
    ## Halichoeres_leucurus             0.76623 -0.64257 0.1301 0.41658  
    ## Halichoeres_marginatus           0.97015 -0.24250 0.0123 1.00000  
    ## Halichoeres_melanurus            0.16567 -0.98618 0.1794 0.06893 .
    ## Halichoeres_papilionaceus        0.93319 -0.35938 0.0594 0.94306  
    ## Halichoeres_trimaculatus         0.98642 -0.16424 0.0324 0.93906  
    ## Hemiglyphidodon_plagiometopon    0.37499 -0.92703 0.2710 0.18881  
    ## Hemigymnus_melapterus            0.97345 -0.22889 0.0416 0.93806  
    ## Heniochus_chrysostomus           0.92665 -0.37593 0.0734 1.00000  
    ## Heniochus_monoceros              0.89946 -0.43700 0.0559 0.86014  
    ## Herklotsichthys_quadrimaculatus  0.97015 -0.24250 0.0123 1.00000  
    ## Hipposcarus_longiceps            0.59050 -0.80704 0.1375 0.72827  
    ## Istigobius_decoratus             0.97015 -0.24250 0.0123 1.00000  
    ## Koumansetta_rainfordi            0.97551 -0.21996 0.0291 0.67532  
    ## Kyphosus_cinerascens            -0.12570 -0.99207 0.1021 0.23576  
    ## Lethrinus_erythropterus          0.99831 -0.05809 0.0735 0.93107  
    ## Lethrinus_harak                  0.92136  0.38871 0.1425 0.43556  
    ## Lethrinus_obsoletus              0.64556  0.76371 0.0286 0.26374  
    ## Lethrinus_olivaceus              0.64556  0.76371 0.0286 0.26374  
    ## Lutjanus_argentimaculatus        0.59956  0.80033 0.0950 0.26174  
    ## Lutjanus_biguttatus             -0.13360 -0.99104 0.0695 0.31269  
    ## Lutjanus_bohar                   0.86314 -0.50496 0.0554 0.34665  
    ## Lutjanus_decussatus              0.64556  0.76371 0.0286 0.26374  
    ## Lutjanus_ehrenbergii             0.88272  0.46991 0.1130 0.50450  
    ## Lutjanus_fulviflamma             0.98042 -0.19690 0.0795 0.39960  
    ## Lutjanus_fulvus                  0.74252  0.66982 0.3371 0.04296 *
    ## Lutjanus_gibbus                  0.99732 -0.07323 0.0916 0.71828  
    ## Lutjanus_monostigma              0.90064  0.43457 0.0269 0.83616  
    ## Lutjanus_russellii               0.99677  0.08032 0.0184 0.75724  
    ## Macrodontogobius_wilburi         0.58047 -0.81428 0.2239 0.17882  
    ## Megalops_cyprinoides            -0.21977  0.97555 0.4589 0.11988  
    ## Meiacanthus_ditrema              0.60930 -0.79294 0.0369 0.50150  
    ## Meiacanthus_grammistes           0.80369 -0.59505 0.0647 0.65734  
    ## Monodactylus_argenteus           0.50255  0.86455 0.0575 0.41958  
    ## Monotaxis_grandoculis            0.05563 -0.99845 0.0424 0.59341  
    ## Mugilogobius_cavifrons          -0.95857  0.28486 0.1926 0.52647  
    ## Mulloidichthys_flavolineatus     0.97015 -0.24250 0.0123 1.00000  
    ## Myripristis_adusta               0.97925 -0.20266 0.1337 0.75524  
    ## Myripristis_sp                   0.99496 -0.10025 0.0160 0.86114  
    ## Myripristis_violacea             0.97551 -0.21996 0.0291 0.67532  
    ## Naso_brevirostris                0.98962 -0.14373 0.1029 0.49351  
    ## Naso_vlamingii                   0.95166 -0.30714 0.0779 0.61538  
    ## Nectamia_viria                   0.69970 -0.71443 0.0258 0.37263  
    ## Neoglyphidodon_melas            -0.44019 -0.89790 0.3364 0.38861  
    ## Neoniphon_argenteus              0.86892 -0.49495 0.0527 0.72028  
    ## Neoniphon_opercularis           -0.44019 -0.89790 0.3364 0.38861  
    ## Neoniphon_sammara                0.94327 -0.33202 0.1108 0.82118  
    ## Neopomacentrus_nemurus           0.28545 -0.95839 0.1664 0.11688  
    ## Omobranchus_obliquus             0.51251  0.85868 0.0627 0.37762  
    ## Oplopomops_diacanthus            0.97015 -0.24250 0.0123 1.00000  
    ## Oxycheilinus_celebicus           0.17006 -0.98543 0.1401 0.14585  
    ## Parapercis_cylindrica            0.97015 -0.24250 0.0123 1.00000  
    ## Parioglossus_formosus           -0.11512  0.99335 0.1649 0.12987  
    ## Parioglossus_philippinus         0.97303 -0.23068 0.0666 0.57842  
    ## Parioglossus_sp                  0.98551 -0.16963 0.0463 0.62737  
    ## Parupeneus_barberinus            0.93356 -0.35841 0.2063 0.68731  
    ## Parupeneus_multifasciatus        0.75567 -0.65496 0.2919 0.11788  
    ## Pentapodus_trivittatus           0.79542 -0.60606 0.2723 0.60939  
    ## Periophthalmus_argentilineatus   0.99507 -0.09914 0.0189 0.63337  
    ## Pervagor_janthinosoma            0.81508 -0.57934 0.0282 0.82817  
    ## Pervagor_nigrolineatus           0.46280 -0.88646 0.0943 0.32867  
    ## Petroscirtes_breviceps           0.84342 -0.53725 0.1306 0.70230  
    ## Platax_orbicularis               0.86641  0.49933 0.0415 0.33766  
    ## Platax_pinnatus                  0.62285 -0.78234 0.0230 0.49351  
    ## Platax_sp                        0.54640  0.83752 0.0318 0.13886  
    ## Plectorhinchus_albovittatus      0.99677  0.08032 0.0184 0.75724  
    ## Plectorhinchus_chaetodonoides    0.22647 -0.97402 0.0837 0.17083  
    ## Plectorhinchus_gibbosus          0.97015 -0.24250 0.0123 1.00000  
    ## Plectropomus_leopardus           0.97551 -0.21996 0.0291 0.67532  
    ## Pomacanthus_sexstriatus          0.82653 -0.56290 0.0374 0.85614  
    ## Pomacentrus_adelus               0.99507 -0.09914 0.0189 0.63337  
    ## Pomacentrus_amboinensis          0.97551 -0.21996 0.0291 0.67532  
    ## Pomacentrus_burroughi            0.66367 -0.74803 0.2741 0.05794 .
    ## Pomacentrus_grammorhynchus       0.97015 -0.24250 0.0123 1.00000  
    ## Pomacentrus_nigromanus          -0.09117 -0.99584 0.1122 0.13586  
    ## Pomacentrus_pavo                 0.97551 -0.21996 0.0291 0.67532  
    ## Pomacentrus_simsiang             0.82943 -0.55860 0.3555 0.71928  
    ## Pomacentrus_sp                  -0.31964 -0.94754 0.1460 0.22278  
    ## Pomacentrus_taeniometopon        0.99677  0.08032 0.0184 0.75724  
    ## Pseudobalistes_flavimarginatus   0.91882 -0.39467 0.1694 0.80519  
    ## Pseudochromis_fuscus             0.97551 -0.21996 0.0291 0.67532  
    ## Ptereleotris_brachyptera         0.83397 -0.55181 0.0513 0.45255  
    ## Ptereleotris_evides              0.97551 -0.21996 0.0291 0.67532  
    ## Pterocaesio_tile                 0.64556  0.76371 0.0286 0.26374  
    ## Pterocaesio_trilineata           0.69970 -0.71443 0.0258 0.37263  
    ## Pterois_volitans                 0.62285 -0.78234 0.0230 0.49351  
    ## Rastrelliger_kanagurta           0.41320 -0.91064 0.0943 0.27872  
    ## Rhabdamia_gracilis               0.88734  0.46111 0.0849 0.29471  
    ## Rhinecanthus_aculeatus           0.94141 -0.33725 0.0782 0.61538  
    ## Salarias_segmentatus             0.89816 -0.43967 0.2200 0.20280  
    ## Sargocentron_diadema             0.75923 -0.65082 0.0564 0.32767  
    ## Sargocentron_spiniferum          0.99496 -0.10025 0.0160 0.86114  
    ## Saurida_gracilis                 0.81508 -0.57934 0.0282 0.82817  
    ## Scarus_dimidiatus                0.84181 -0.53977 0.1148 1.00000  
    ## Scarus_flavipectoralis           0.90905 -0.41668 0.0592 0.81219  
    ## Scarus_ghobban                   0.98551 -0.16963 0.0463 0.62737  
    ## Scarus_hypselopterus             0.81508 -0.57934 0.0282 0.82817  
    ## Scarus_oviceps                   0.60930 -0.79294 0.0369 0.50150  
    ## Scarus_psittacus                 0.85607 -0.51686 0.1350 0.64336  
    ## Scarus_quoyi                     0.97551 -0.21996 0.0291 0.67532  
    ## Scarus_rivulatus                 0.47726 -0.87876 0.0849 0.39860  
    ## Scatophagus_argus                0.63276  0.77434 0.0173 1.00000  
    ## Scolopsis_ciliata                0.78189 -0.62342 0.2106 0.62737  
    ## Scolopsis_lineata                0.98583 -0.16776 0.0296 0.97502  
    ## Scolopsis_margaritifera          0.30245 -0.95316 0.2360 0.10689  
    ## Scolopsis_trilineata             0.99958  0.02897 0.0611 0.42058  
    ## Siganus_doliatus                 0.58325 -0.81229 0.0901 0.81219  
    ## Siganus_fuscescens               0.99755  0.06989 0.0612 0.92208  
    ## Siganus_lineatus                 0.73021  0.68323 0.2811 0.01798 *
    ## Siganus_puellus                  0.51858 -0.85503 0.1292 0.33866  
    ## Siganus_punctatissimus           0.97551 -0.21996 0.0291 0.67532  
    ## Siganus_punctatus                0.97551 -0.21996 0.0291 0.67532  
    ## Signigobius_biocellatus          0.81508 -0.57934 0.0282 0.82817  
    ## Silhouettea_sp                   0.07760  0.99698 0.0817 0.74426  
    ## Sphaeramia_nematoptera           0.59124 -0.80650 0.1951 0.46553  
    ## Sphaeramia_orbicularis           0.70741 -0.70680 0.1070 0.44056  
    ## Sphyraena_barracuda              0.88644  0.46285 0.1165 0.40559  
    ## Sphyraena_obtusata               0.88332 -0.46877 0.0097 1.00000  
    ## Sphyraena_qenie                  0.99496 -0.10025 0.0160 0.86114  
    ## Stegastes_nigricans              0.97015 -0.24250 0.0123 1.00000  
    ## Stegastes_punctatus              0.88713 -0.46152 0.0403 1.00000  
    ## Stethojulis_strigiventer         0.98975 -0.14278 0.0518 0.93706  
    ## Strongylura_incisa               0.91608 -0.40099 0.0444 0.69730  
    ## Synodus_variegatus               0.41320 -0.91064 0.0943 0.27872  
    ## Taeniamia_zosterophora           0.46969 -0.88283 0.1098 0.54645  
    ## Taeniurops_meyeni                0.64556  0.76371 0.0286 0.26374  
    ## Toxotes_jaculatrix               0.99567 -0.09299 0.1209 0.97103  
    ## Tylosurus_punctulatus            0.72764 -0.68596 0.0320 0.82318  
    ## Upeneus_taeniopterus             0.99496 -0.10025 0.0160 0.86114  
    ## Valenciennea_muralis             0.99999 -0.00433 0.0565 0.62438  
    ## Valenciennea_parva               0.62285 -0.78234 0.0230 0.49351  
    ## Valenciennea_randalli            0.62285 -0.78234 0.0230 0.49351  
    ## Zanclus_cornutus                 0.87645 -0.48148 0.0656 0.85315  
    ## Zebrasoma_scopas                 0.97551 -0.21996 0.0291 0.67532  
    ## Zebrasoma_velifer                0.64129 -0.76729 0.2697 0.26673  
    ## Zenarchopterus_dispar            0.90028 -0.43531 0.3055 0.85215  
    ## Zenarchopterus_gilli             0.64801 -0.76163 0.1532 0.58641  
    ## Zoramia_leptacanthus             0.97015 -0.24250 0.0123 1.00000  
    ## Zoramia_viridiventer             0.71256 -0.70161 0.1084 0.37063  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
fds_efp <- p.adjust.envfit(fds_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
fds_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                                    NMDS1    NMDS2     r2 Pr(>r)
    ## Abudefduf_lorenzi               -0.10992 -0.99394 0.0977      1
    ## Abudefduf_septemfasciatus        0.99507 -0.09914 0.0189      1
    ## Abudefduf_sexfasciatus          -0.10992 -0.99394 0.0977      1
    ## Acanthurus_lineatus              0.81508 -0.57934 0.0282      1
    ## Acanthurus_nigricauda            0.89704 -0.44195 0.0698      1
    ## Acanthurus_sp                    0.60930 -0.79294 0.0369      1
    ## Acanthurus_xanthopterus          0.87335  0.48710 0.2458      1
    ## Acentrogobius_janthinopterus    -0.69949  0.71465 0.5665      1
    ## Amblyeleotris_gymnocephala       0.62285 -0.78234 0.0230      1
    ## Amblyglyphidodon_curacao         0.30794 -0.95140 0.1771      1
    ## Amblygobius_buanensis            0.97524 -0.22117 0.3733      1
    ## Amblygobius_decussatus           0.47726 -0.87876 0.0849      1
    ## Amblygobius_esakiae              0.60930 -0.79294 0.0369      1
    ## Amblygobius_linki                0.99677  0.08032 0.0184      1
    ## Amblygobius_phalaena             0.78723 -0.61666 0.0343      1
    ## Amphiprion_clarkii               0.97551 -0.21996 0.0291      1
    ## Ancistrogobius_dipus             0.64556  0.76371 0.0286      1
    ## Apogonichthyoides_melas          0.64556  0.76371 0.0286      1
    ## Arcygobius_baliurus              0.99496 -0.10025 0.0160      1
    ## Arothron_reticularis             0.63276  0.77434 0.0173      1
    ## Asterropteryx_semipunctata       0.85047 -0.52603 0.0633      1
    ## Asterropteryx_sp                 0.99677  0.08032 0.0184      1
    ## Atherinomorus_endrachtensis      0.60712 -0.79461 0.1136      1
    ## Atherinomorus_sp                 0.84878 -0.52874 0.0313      1
    ## Atule_mate                       0.81508 -0.57934 0.0282      1
    ## Balistoides_viridescens          0.97754 -0.21073 0.1518      1
    ## Bolbometopon_muricatum           0.58325 -0.81229 0.0901      1
    ## Caesio_caerulaurea               0.99620  0.08712 0.0414      1
    ## Caesio_cuning                    0.56489 -0.82517 0.1710      1
    ## Caesio_teres                     0.62285 -0.78234 0.0230      1
    ## Cantherhines_dumerilii           0.97551 -0.21996 0.0291      1
    ## Canthigaster_bennetti            0.97551 -0.21996 0.0291      1
    ## Canthigaster_solandri            0.99507 -0.09914 0.0189      1
    ## Caranx_ignobilis                 0.97015 -0.24250 0.0123      1
    ## Caranx_melampygus                0.51541 -0.85695 0.0123      1
    ## Caranx_papuensis                 0.62285 -0.78234 0.0230      1
    ## Caranx_sexfasciatus              0.66234  0.74921 0.1634      1
    ## Cephalopholis_argus             -0.44019 -0.89790 0.3364      1
    ## Cephalopholis_boenak             0.80369 -0.59505 0.0647      1
    ## Chaetodon_auriga                 0.91053 -0.41344 0.3545      1
    ## Chaetodon_bennetti               0.90905 -0.41668 0.0592      1
    ## Chaetodon_ephippium              0.69774 -0.71635 0.3357      1
    ## Chaetodon_kleinii                0.82742 -0.56159 0.0851      1
    ## Chaetodon_lunula                 0.96902 -0.24697 0.1702      1
    ## Chaetodon_lunulatus              0.58047 -0.81428 0.2239      1
    ## Chaetodon_melannotus             0.69970 -0.71443 0.0258      1
    ## Chaetodon_ocellicaudus           0.97015 -0.24250 0.0123      1
    ## Chaetodon_octofasciatus          0.80369 -0.59505 0.0647      1
    ## Chaetodon_rafflesii              0.84181 -0.53977 0.1148      1
    ## Chaetodon_semeion                0.86314 -0.50496 0.0554      1
    ## Chaetodon_trifascialis           0.97551 -0.21996 0.0291      1
    ## Chaetodon_ulietensis             0.61123 -0.79145 0.1247      1
    ## Chaetodon_vagabundus             0.85994 -0.51039 0.3542      1
    ## Chaetodontoplus_poliourus        0.97551 -0.21996 0.0291      1
    ## Cheilinus_fasciatus              0.80739 -0.59002 0.0993      1
    ## Cheilinus_trilobatus             0.81508 -0.57934 0.0282      1
    ## Cheilinus_undulatus              0.40084 -0.91615 0.0676      1
    ## Cheilodipterus_artus             0.18546 -0.98265 0.0782      1
    ## Cheilodipterus_isostigma         0.43580 -0.90004 0.0691      1
    ## Cheilodipterus_quinquelineatus   0.77642 -0.63021 0.3290      1
    ## Cheilodipterus_singapurensis     0.77208 -0.63553 0.0348      1
    ## Chlorurus_bleekeri               0.58087 -0.81400 0.1647      1
    ## Chlorurus_microrhinos            0.97345 -0.22889 0.0416      1
    ## Chlorurus_spilurus               0.84181 -0.53977 0.1148      1
    ## Choerodon_anchorago              0.92604 -0.37743 0.4596      1
    ## Chromis_viridis                  0.48637 -0.87375 0.0962      1
    ## Chrysiptera_biocellata           0.99502 -0.09967 0.0366      1
    ## Chrysiptera_oxycephala           0.30534 -0.95224 0.2330      1
    ## Corythoichthys_ocellatus         0.97551 -0.21996 0.0291      1
    ## Cryptocentrus_caeruleomaculatus  0.85768  0.51419 0.0763      1
    ## Cryptocentrus_cyanospilotus      0.07760  0.99698 0.0817      1
    ## Cryptocentrus_leptocephalus      0.98583 -0.16776 0.0296      1
    ## Cryptocentrus_strigilliceps      0.59206 -0.80590 0.1206      1
    ## Ctenochaetus_striatus            0.19139 -0.98151 0.2053      1
    ## Ctenogobiops_pomastictus         0.64556  0.76371 0.0286      1
    ## Cymbacephalus_beauforti          0.97015 -0.24250 0.0123      1
    ## Dascyllus_aruanus                0.75557 -0.65507 0.0968      1
    ## Dascyllus_trimaculatus           0.97551 -0.21996 0.0291      1
    ## Diplogrammus_goramensis          0.87425  0.48547 0.0363      1
    ## Diproctacanthus_xanthurus        0.90905 -0.41668 0.0592      1
    ## Dischistodus_chrysopoecilus      0.80632  0.59148 0.0369      1
    ## Dischistodus_perspicillatus      0.77943 -0.62649 0.1554      1
    ## Eleotris_fusca                  -0.29588  0.95523 0.5192      1
    ## Epibulus_brevis                  0.58785 -0.80897 0.2014      1
    ## Epinephelus_coeruleopunctatus    0.99998 -0.00701 0.0360      1
    ## Epinephelus_merra                0.49995 -0.86605 0.1688      1
    ## Epinephelus_sp                   0.54640  0.83752 0.0318      1
    ## Eviota_atriventris               0.99806 -0.06226 0.0696      1
    ## Eviota_bifasciata                0.71997 -0.69401 0.2562      1
    ## Eviota_fallax                    0.60437 -0.79670 0.2469      1
    ## Eviota_lachdeberei               0.95552 -0.29491 0.0796      1
    ## Eviota_maculosa                  0.22647 -0.97402 0.0837      1
    ## Eviota_sigillata                 0.66215 -0.74937 0.0511      1
    ## Eviota_sp                        0.99677  0.08032 0.0184      1
    ## Eviota_storthynx                 0.89946 -0.43700 0.0559      1
    ## Exyrias_belissimus               0.80190 -0.59745 0.1812      1
    ## Exyrias_puntang                 -0.07025  0.99753 0.3998      1
    ## Favonigobius_reichei             0.98583 -0.16776 0.0296      1
    ## Fibramia_lateralis               0.04521  0.99898 0.0111      1
    ## Fibramia_thermalis               0.93933 -0.34301 0.0919      1
    ## Fusigobius_signipinnis           0.87645 -0.48148 0.0656      1
    ## Gerres_erythrourus               0.97015 -0.24250 0.0123      1
    ## Gerres_oyena                     0.98642 -0.16424 0.0324      1
    ## Gladiogobius_ensifer             0.22647 -0.97402 0.0837      1
    ## Gnathanodon_speciosus            0.25895  0.96589 0.0915      1
    ## Gobiodon_okinawae                0.81508 -0.57934 0.0282      1
    ## Gymnothorax_pictus               0.99502 -0.09967 0.0366      1
    ## Halichoeres_chloropterus         0.42901 -0.90330 0.1679      1
    ## Halichoeres_leucurus             0.76623 -0.64257 0.1301      1
    ## Halichoeres_marginatus           0.97015 -0.24250 0.0123      1
    ## Halichoeres_melanurus            0.16567 -0.98618 0.1794      1
    ## Halichoeres_papilionaceus        0.93319 -0.35938 0.0594      1
    ## Halichoeres_trimaculatus         0.98642 -0.16424 0.0324      1
    ## Hemiglyphidodon_plagiometopon    0.37499 -0.92703 0.2710      1
    ## Hemigymnus_melapterus            0.97345 -0.22889 0.0416      1
    ## Heniochus_chrysostomus           0.92665 -0.37593 0.0734      1
    ## Heniochus_monoceros              0.89946 -0.43700 0.0559      1
    ## Herklotsichthys_quadrimaculatus  0.97015 -0.24250 0.0123      1
    ## Hipposcarus_longiceps            0.59050 -0.80704 0.1375      1
    ## Istigobius_decoratus             0.97015 -0.24250 0.0123      1
    ## Koumansetta_rainfordi            0.97551 -0.21996 0.0291      1
    ## Kyphosus_cinerascens            -0.12570 -0.99207 0.1021      1
    ## Lethrinus_erythropterus          0.99831 -0.05809 0.0735      1
    ## Lethrinus_harak                  0.92136  0.38871 0.1425      1
    ## Lethrinus_obsoletus              0.64556  0.76371 0.0286      1
    ## Lethrinus_olivaceus              0.64556  0.76371 0.0286      1
    ## Lutjanus_argentimaculatus        0.59956  0.80033 0.0950      1
    ## Lutjanus_biguttatus             -0.13360 -0.99104 0.0695      1
    ## Lutjanus_bohar                   0.86314 -0.50496 0.0554      1
    ## Lutjanus_decussatus              0.64556  0.76371 0.0286      1
    ## Lutjanus_ehrenbergii             0.88272  0.46991 0.1130      1
    ## Lutjanus_fulviflamma             0.98042 -0.19690 0.0795      1
    ## Lutjanus_fulvus                  0.74252  0.66982 0.3371      1
    ## Lutjanus_gibbus                  0.99732 -0.07323 0.0916      1
    ## Lutjanus_monostigma              0.90064  0.43457 0.0269      1
    ## Lutjanus_russellii               0.99677  0.08032 0.0184      1
    ## Macrodontogobius_wilburi         0.58047 -0.81428 0.2239      1
    ## Megalops_cyprinoides            -0.21977  0.97555 0.4589      1
    ## Meiacanthus_ditrema              0.60930 -0.79294 0.0369      1
    ## Meiacanthus_grammistes           0.80369 -0.59505 0.0647      1
    ## Monodactylus_argenteus           0.50255  0.86455 0.0575      1
    ## Monotaxis_grandoculis            0.05563 -0.99845 0.0424      1
    ## Mugilogobius_cavifrons          -0.95857  0.28486 0.1926      1
    ## Mulloidichthys_flavolineatus     0.97015 -0.24250 0.0123      1
    ## Myripristis_adusta               0.97925 -0.20266 0.1337      1
    ## Myripristis_sp                   0.99496 -0.10025 0.0160      1
    ## Myripristis_violacea             0.97551 -0.21996 0.0291      1
    ## Naso_brevirostris                0.98962 -0.14373 0.1029      1
    ## Naso_vlamingii                   0.95166 -0.30714 0.0779      1
    ## Nectamia_viria                   0.69970 -0.71443 0.0258      1
    ## Neoglyphidodon_melas            -0.44019 -0.89790 0.3364      1
    ## Neoniphon_argenteus              0.86892 -0.49495 0.0527      1
    ## Neoniphon_opercularis           -0.44019 -0.89790 0.3364      1
    ## Neoniphon_sammara                0.94327 -0.33202 0.1108      1
    ## Neopomacentrus_nemurus           0.28545 -0.95839 0.1664      1
    ## Omobranchus_obliquus             0.51251  0.85868 0.0627      1
    ## Oplopomops_diacanthus            0.97015 -0.24250 0.0123      1
    ## Oxycheilinus_celebicus           0.17006 -0.98543 0.1401      1
    ## Parapercis_cylindrica            0.97015 -0.24250 0.0123      1
    ## Parioglossus_formosus           -0.11512  0.99335 0.1649      1
    ## Parioglossus_philippinus         0.97303 -0.23068 0.0666      1
    ## Parioglossus_sp                  0.98551 -0.16963 0.0463      1
    ## Parupeneus_barberinus            0.93356 -0.35841 0.2063      1
    ## Parupeneus_multifasciatus        0.75567 -0.65496 0.2919      1
    ## Pentapodus_trivittatus           0.79542 -0.60606 0.2723      1
    ## Periophthalmus_argentilineatus   0.99507 -0.09914 0.0189      1
    ## Pervagor_janthinosoma            0.81508 -0.57934 0.0282      1
    ## Pervagor_nigrolineatus           0.46280 -0.88646 0.0943      1
    ## Petroscirtes_breviceps           0.84342 -0.53725 0.1306      1
    ## Platax_orbicularis               0.86641  0.49933 0.0415      1
    ## Platax_pinnatus                  0.62285 -0.78234 0.0230      1
    ## Platax_sp                        0.54640  0.83752 0.0318      1
    ## Plectorhinchus_albovittatus      0.99677  0.08032 0.0184      1
    ## Plectorhinchus_chaetodonoides    0.22647 -0.97402 0.0837      1
    ## Plectorhinchus_gibbosus          0.97015 -0.24250 0.0123      1
    ## Plectropomus_leopardus           0.97551 -0.21996 0.0291      1
    ## Pomacanthus_sexstriatus          0.82653 -0.56290 0.0374      1
    ## Pomacentrus_adelus               0.99507 -0.09914 0.0189      1
    ## Pomacentrus_amboinensis          0.97551 -0.21996 0.0291      1
    ## Pomacentrus_burroughi            0.66367 -0.74803 0.2741      1
    ## Pomacentrus_grammorhynchus       0.97015 -0.24250 0.0123      1
    ## Pomacentrus_nigromanus          -0.09117 -0.99584 0.1122      1
    ## Pomacentrus_pavo                 0.97551 -0.21996 0.0291      1
    ## Pomacentrus_simsiang             0.82943 -0.55860 0.3555      1
    ## Pomacentrus_sp                  -0.31964 -0.94754 0.1460      1
    ## Pomacentrus_taeniometopon        0.99677  0.08032 0.0184      1
    ## Pseudobalistes_flavimarginatus   0.91882 -0.39467 0.1694      1
    ## Pseudochromis_fuscus             0.97551 -0.21996 0.0291      1
    ## Ptereleotris_brachyptera         0.83397 -0.55181 0.0513      1
    ## Ptereleotris_evides              0.97551 -0.21996 0.0291      1
    ## Pterocaesio_tile                 0.64556  0.76371 0.0286      1
    ## Pterocaesio_trilineata           0.69970 -0.71443 0.0258      1
    ## Pterois_volitans                 0.62285 -0.78234 0.0230      1
    ## Rastrelliger_kanagurta           0.41320 -0.91064 0.0943      1
    ## Rhabdamia_gracilis               0.88734  0.46111 0.0849      1
    ## Rhinecanthus_aculeatus           0.94141 -0.33725 0.0782      1
    ## Salarias_segmentatus             0.89816 -0.43967 0.2200      1
    ## Sargocentron_diadema             0.75923 -0.65082 0.0564      1
    ## Sargocentron_spiniferum          0.99496 -0.10025 0.0160      1
    ## Saurida_gracilis                 0.81508 -0.57934 0.0282      1
    ## Scarus_dimidiatus                0.84181 -0.53977 0.1148      1
    ## Scarus_flavipectoralis           0.90905 -0.41668 0.0592      1
    ## Scarus_ghobban                   0.98551 -0.16963 0.0463      1
    ## Scarus_hypselopterus             0.81508 -0.57934 0.0282      1
    ## Scarus_oviceps                   0.60930 -0.79294 0.0369      1
    ## Scarus_psittacus                 0.85607 -0.51686 0.1350      1
    ## Scarus_quoyi                     0.97551 -0.21996 0.0291      1
    ## Scarus_rivulatus                 0.47726 -0.87876 0.0849      1
    ## Scatophagus_argus                0.63276  0.77434 0.0173      1
    ## Scolopsis_ciliata                0.78189 -0.62342 0.2106      1
    ## Scolopsis_lineata                0.98583 -0.16776 0.0296      1
    ## Scolopsis_margaritifera          0.30245 -0.95316 0.2360      1
    ## Scolopsis_trilineata             0.99958  0.02897 0.0611      1
    ## Siganus_doliatus                 0.58325 -0.81229 0.0901      1
    ## Siganus_fuscescens               0.99755  0.06989 0.0612      1
    ## Siganus_lineatus                 0.73021  0.68323 0.2811      1
    ## Siganus_puellus                  0.51858 -0.85503 0.1292      1
    ## Siganus_punctatissimus           0.97551 -0.21996 0.0291      1
    ## Siganus_punctatus                0.97551 -0.21996 0.0291      1
    ## Signigobius_biocellatus          0.81508 -0.57934 0.0282      1
    ## Silhouettea_sp                   0.07760  0.99698 0.0817      1
    ## Sphaeramia_nematoptera           0.59124 -0.80650 0.1951      1
    ## Sphaeramia_orbicularis           0.70741 -0.70680 0.1070      1
    ## Sphyraena_barracuda              0.88644  0.46285 0.1165      1
    ## Sphyraena_obtusata               0.88332 -0.46877 0.0097      1
    ## Sphyraena_qenie                  0.99496 -0.10025 0.0160      1
    ## Stegastes_nigricans              0.97015 -0.24250 0.0123      1
    ## Stegastes_punctatus              0.88713 -0.46152 0.0403      1
    ## Stethojulis_strigiventer         0.98975 -0.14278 0.0518      1
    ## Strongylura_incisa               0.91608 -0.40099 0.0444      1
    ## Synodus_variegatus               0.41320 -0.91064 0.0943      1
    ## Taeniamia_zosterophora           0.46969 -0.88283 0.1098      1
    ## Taeniurops_meyeni                0.64556  0.76371 0.0286      1
    ## Toxotes_jaculatrix               0.99567 -0.09299 0.1209      1
    ## Tylosurus_punctulatus            0.72764 -0.68596 0.0320      1
    ## Upeneus_taeniopterus             0.99496 -0.10025 0.0160      1
    ## Valenciennea_muralis             0.99999 -0.00433 0.0565      1
    ## Valenciennea_parva               0.62285 -0.78234 0.0230      1
    ## Valenciennea_randalli            0.62285 -0.78234 0.0230      1
    ## Zanclus_cornutus                 0.87645 -0.48148 0.0656      1
    ## Zebrasoma_scopas                 0.97551 -0.21996 0.0291      1
    ## Zebrasoma_velifer                0.64129 -0.76729 0.2697      1
    ## Zenarchopterus_dispar            0.90028 -0.43531 0.3055      1
    ## Zenarchopterus_gilli             0.64801 -0.76163 0.1532      1
    ## Zoramia_leptacanthus             0.97015 -0.24250 0.0123      1
    ## Zoramia_viridiventer             0.71256 -0.70161 0.1084      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

## Trait mantel tests

- We used mantel tests to determine significance between the trait and
  environmental distance matrices.

``` r
### Environmental
# All sites
enve_dist_t <- dist(scaled_env[surveyed_sites_env,c(1)], method = "euclidean")
fdea_mant_t <- mantel(fde_dist, enve_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
fdea_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fde_dist, ydis = enve_dist_t, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.347 
    ##       Significance: 0.294 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.443 0.487 0.532 0.572 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_s <- dist(scaled_env[surveyed_sites_env,c(2)], method = "euclidean")
fdea_mant_s <- mantel(fde_dist, enve_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
fdea_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fde_dist, ydis = enve_dist_s, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.8065 
    ##       Significance: 0.017 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.703 0.757 0.794 0.829 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_o <- dist(scaled_env[surveyed_sites_env,c(3)], method = "euclidean")
fdea_mant_o <- mantel(fde_dist, enve_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
fdea_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fde_dist, ydis = enve_dist_o, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1983 
    ##       Significance: 0.961 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.526 0.562 0.601 0.630 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_p <- dist(scaled_env[surveyed_sites_env,c(4)], method = "euclidean")
fdea_mant_p <- mantel(fde_dist, enve_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
fdea_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fde_dist, ydis = enve_dist_p, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4257 
    ##       Significance: 0.158 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.463 0.505 0.534 0.547 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
fdea_mant_pv <- rbind(fdea_mant_t$signif, fdea_mant_s$signif, fdea_mant_o$signif, fdea_mant_p$signif)
fdea_mant_pv <- fdea_mant_pv[,1]
fdea_mant_pv <- p.adjust(fdea_mant_pv, method = "bonferroni")
fdea_mant_pv
```

    ## [1] 1.000 0.068 1.000 0.632

``` r
# Mixed and stratified lakes
envems_dist_t <- dist(scaled_env[mixed_stratified_lakes,c(1)], method = "euclidean")
fdems_mant_t <- mantel(fdems_dist, envems_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
fdems_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdems_dist, ydis = envems_dist_t, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3515 
    ##       Significance: 0.267 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.459 0.517 0.561 0.594 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_s <- dist(scaled_env[mixed_stratified_lakes,c(2)], method = "euclidean")
fdems_mant_s <- mantel(fdems_dist, envems_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
fdems_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdems_dist, ydis = envems_dist_s, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7968 
    ##       Significance: 0.01 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.647 0.706 0.760 0.794 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_o <- dist(scaled_env[mixed_stratified_lakes,c(3)], method = "euclidean")
fdems_mant_o <- mantel(fdems_dist, envems_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
fdems_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdems_dist, ydis = envems_dist_o, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.07334 
    ##       Significance: 0.944 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.482 0.538 0.574 0.623 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_p <- dist(scaled_env[mixed_stratified_lakes,c(4)], method = "euclidean")
fdems_mant_p <- mantel(fdems_dist, envems_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
fdems_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdems_dist, ydis = envems_dist_p, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3022 
    ##       Significance: 0.202 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.363 0.409 0.448 0.476 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
fdems_mant_pv <- rbind(fdems_mant_t$signif, fdems_mant_s$signif, fdems_mant_o$signif, fdems_mant_p$signif)
fdems_mant_pv <- fdems_mant_pv[,1]
fdems_mant_pv <- p.adjust(fdems_mant_pv, method = "bonferroni")
fdems_mant_pv
```

    ## [1] 1.000 0.040 1.000 0.808

``` r
# Ocean sites and mixed lakes
enveom_dist_t <- dist(scaled_env[ocean_mixed_sites_env,c(1)], method = "euclidean")
fdeom_mant_t <- mantel(fdeom_dist, enveom_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
fdeom_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdeom_dist, ydis = enveom_dist_t, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02468 
    ##       Significance: 0.51 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.251 0.316 0.365 0.401 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_dist_s <- dist(scaled_env[ocean_mixed_sites_env,c(2)], method = "euclidean")
fdeom_mant_s <- mantel(fdeom_dist, enveom_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
fdeom_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdeom_dist, ydis = enveom_dist_s, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6041 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.273 0.326 0.373 0.411 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_dist_o <- dist(scaled_env[ocean_mixed_sites_env,c(3)], method = "euclidean")
fdeom_mant_o <- mantel(fdeom_dist, enveom_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
fdeom_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdeom_dist, ydis = enveom_dist_o, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2038 
    ##       Significance: 0.222 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.319 0.402 0.492 0.578 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_dist_p <- dist(scaled_env[ocean_mixed_sites_env,c(4)], method = "euclidean")
fdeom_mant_p <- mantel(fdeom_dist, enveom_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
fdeom_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdeom_dist, ydis = enveom_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.0702 
    ##       Significance: 0.279 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.216 0.287 0.343 0.398 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
fdeom_mant_pv <- rbind(fdeom_mant_t$signif, fdeom_mant_s$signif, fdeom_mant_o$signif, fdeom_mant_p$signif)
fdeom_mant_pv <- fdeom_mant_pv[,1]
fdeom_mant_pv <- p.adjust(fdeom_mant_pv, method = "bonferroni")
fdeom_mant_pv
```

    ## [1] 1.000 0.004 0.888 1.000

``` r
# Stratified lakes and ocean sites
enveso_dist_t <- dist(scaled_env[ocean_stratified_sites_env,c(1)], method = "euclidean")
fdeso_mant_t <- mantel(fdeso_dist, enveso_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
fdeso_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdeso_dist, ydis = enveso_dist_t, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02626 
    ##       Significance: 0.568 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.205 0.255 0.300 0.360 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveso_dist_s <- dist(scaled_env[ocean_stratified_sites_env,c(2)], method = "euclidean")
fdeso_mant_s <- mantel(fdeso_dist, enveso_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
fdeso_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdeso_dist, ydis = enveso_dist_s, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6191 
    ##       Significance: 0.021 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.418 0.495 0.595 0.646 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveso_dist_o <- dist(scaled_env[ocean_stratified_sites_env,c(3)], method = "euclidean")
fdeso_mant_o <- mantel(fdeso_dist, enveso_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
fdeso_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdeso_dist, ydis = enveso_dist_o, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1364 
    ##       Significance: 0.768 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.373 0.437 0.504 0.544 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveso_dist_p <- dist(scaled_env[ocean_stratified_sites_env,c(4)], method = "euclidean")
fdeso_mant_p <- mantel(fdeso_dist, enveso_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
fdeso_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdeso_dist, ydis = enveso_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3058 
    ##       Significance: 0.253 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.413 0.483 0.530 0.590 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
fdeso_mant_pv <- rbind(fdeso_mant_t$signif, fdeso_mant_s$signif, fdeso_mant_o$signif, fdeso_mant_p$signif)
fdeso_mant_pv <- fdeso_mant_pv[,1]
fdeso_mant_pv <- p.adjust(fdeso_mant_pv, method = "bonferroni")
fdeso_mant_pv
```

    ## [1] 1.000 0.084 1.000 1.000

``` r
# Stratified lakes
envs_dist <- dist(scaled_env[stratified_lakes,c(1:15)], method = "euclidean")
fds_mant <- mantel(fds_dist, envs_dist, method = "spearman", permutations = 999, na.rm = TRUE)
fds_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fds_dist, ydis = envs_dist, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.05473 
    ##       Significance: 0.584 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.250 0.333 0.423 0.516 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# All sites
envb_dist_vc <- dist(scaled_env[surveyed_sites_geo,c(5)], method = "euclidean")
fdba_mant_vc <- mantel(fdb_dist, envb_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_vc, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.06623 
    ##       Significance: 0.79 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.104 0.130 0.146 0.171 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_v <- dist(scaled_env[surveyed_sites_geo,c(6)], method = "euclidean")
fdba_mant_v <- mantel(fdb_dist, envb_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_v, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.03113 
    ##       Significance: 0.557 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.170 0.206 0.241 0.276 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_sa <- dist(scaled_env[surveyed_sites_geo,c(7)], method = "euclidean")
fdba_mant_sa <- mantel(fdb_dist, envb_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_sa, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.02536 
    ##       Significance: 0.575 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.154 0.197 0.218 0.263 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_dmin <- dist(scaled_env[surveyed_sites_geo,c(8)], method = "euclidean")
fdba_mant_dmin <- mantel(fdb_dist, envb_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_dmin, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5945 
    ##       Significance: 0.067 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.560 0.611 0.650 0.668 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_dmean <- dist(scaled_env[surveyed_sites_geo,c(9)], method = "euclidean")
fdba_mant_dmean <- mantel(fdb_dist, envb_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_dmean, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5165 
    ##       Significance: 0.161 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.551 0.590 0.614 0.632 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_dmed <- dist(scaled_env[surveyed_sites_geo,c(10)], method = "euclidean")
fdba_mant_dmed <- mantel(fdb_dist, envb_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_dmed, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4911 
    ##       Significance: 0.176 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.527 0.558 0.582 0.630 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_tl <- dist(scaled_env[surveyed_sites_geo,c(11)], method = "euclidean")
fdba_mant_tl <- mantel(fdb_dist, envb_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_tl, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2741 
    ##       Significance: 0.672 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.390 0.414 0.434 0.453 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_te <- dist(scaled_env[surveyed_sites_geo,c(12)], method = "euclidean")
fdba_mant_te <- mantel(fdb_dist, envb_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_te, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3042 
    ##       Significance: 0.236 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.334 0.356 0.376 0.388 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_p <- dist(scaled_env[surveyed_sites_geo,c(13)], method = "euclidean")
fdba_mant_p <- mantel(fdb_dist, envb_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_p, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1227 
    ##       Significance: 0.733 
    ## 
    ## Upper quantiles of permutations (null model):
    ##      90%      95%    97.5%      99% 
    ## -0.00221  0.02048  0.03785  0.05427 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_md <- dist(scaled_env[surveyed_sites_geo,c(14)], method = "euclidean")
fdba_mant_md <- mantel(fdb_dist, envb_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_md, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.0912 
    ##       Significance: 0.581 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.241 0.281 0.308 0.344 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_la <- dist(scaled_env[surveyed_sites_geo,c(15)], method = "euclidean")
fdba_mant_la <- mantel(fdb_dist, envb_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
fdba_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_la, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.09979 
    ##       Significance: 0.73 
    ## 
    ## Upper quantiles of permutations (null model):
    ##      90%      95%    97.5%      99% 
    ## -0.00014  0.02280  0.03871  0.05580 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
fdba__mant_pv <- rbind(fdba_mant_vc$signif, fdba_mant_v$signif, fdba_mant_sa$signif, fdba_mant_dmin$signif, fdba_mant_dmean$signif, fdba_mant_dmed$signif, fdba_mant_tl$signif, fdba_mant_te$signif, fdba_mant_p$signif, fdba_mant_md$signif, fdba_mant_la$signif)
fdba__mant_pv <- fdba__mant_pv[,1]
fdba__mant_pv <- p.adjust(fdba__mant_pv, method = "bonferroni")
fdba__mant_pv
```

    ##  [1] 1.000 1.000 1.000 0.737 1.000 1.000 1.000 1.000 1.000 1.000 1.000

``` r
# Mixed and stratified lakes 
envbms_dist_vc <- dist(scaled_env[mixed_stratified_lakes_geo,c(5)], method = "euclidean")
fdbms_mant_vc <- mantel(fdbms_dist, envbms_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_vc, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02981 
    ##       Significance: 0.818 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.301 0.346 0.392 0.433 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_v <- dist(scaled_env[mixed_stratified_lakes_geo,c(6)], method = "euclidean")
fdbms_mant_v <- mantel(fdbms_dist, envbms_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_v, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1955 
    ##       Significance: 0.514 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.442 0.508 0.564 0.663 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_sa <- dist(scaled_env[mixed_stratified_lakes_geo,c(7)], method = "euclidean")
fdbms_mant_sa <- mantel(fdbms_dist, envbms_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_sa, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1559 
    ##       Significance: 0.55 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.418 0.502 0.567 0.623 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_dmin <- dist(scaled_env[mixed_stratified_lakes_geo,c(8)], method = "euclidean")
fdbms_mant_dmin <- mantel(fdbms_dist, envbms_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_dmin, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4856 
    ##       Significance: 0.087 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.469 0.549 0.615 0.654 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_dmean <- dist(scaled_env[mixed_stratified_lakes_geo,c(9)], method = "euclidean")
fdbms_mant_dmean <- mantel(fdbms_dist, envbms_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_dmean, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4022 
    ##       Significance: 0.217 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.501 0.574 0.614 0.678 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_dmed <- dist(scaled_env[mixed_stratified_lakes_geo,c(10)], method = "euclidean")
fdbms_mant_dmed <- mantel(fdbms_dist, envbms_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_dmed, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3864 
    ##       Significance: 0.22 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.497 0.561 0.606 0.654 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_tl <- dist(scaled_env[mixed_stratified_lakes_geo,c(11)], method = "euclidean")
fdbms_mant_tl <- mantel(fdbms_dist, envbms_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_tl, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.0524 
    ##       Significance: 0.73 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.137 0.188 0.230 0.275 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_te <- dist(scaled_env[mixed_stratified_lakes_geo,c(12)], method = "euclidean")
fdbms_mant_te <- mantel(fdbms_dist, envbms_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_te, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02818 
    ##       Significance: 0.379 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0442 0.0814 0.1057 0.1451 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_p <- dist(scaled_env[mixed_stratified_lakes_geo,c(13)], method = "euclidean")
fdbms_mant_p <- mantel(fdbms_dist, envbms_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_p, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.03561 
    ##       Significance: 0.66 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.165 0.208 0.239 0.287 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_md <- dist(scaled_env[mixed_stratified_lakes_geo,c(14)], method = "euclidean")
fdbms_mant_md <- mantel(fdbms_dist, envbms_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_md, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.04755 
    ##       Significance: 0.493 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.257 0.322 0.382 0.419 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_la <- dist(scaled_env[mixed_stratified_lakes_geo,c(15)], method = "euclidean")
fdbms_mant_la <- mantel(fdbms_dist, envbms_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
fdbms_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist_la, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.07455 
    ##       Significance: 0.839 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0999 0.1286 0.1473 0.1739 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
fdbms__mant_pv <- rbind(fdbms_mant_vc$signif, fdbms_mant_v$signif, fdbms_mant_sa$signif, fdbms_mant_dmin$signif, fdbms_mant_dmean$signif, fdbms_mant_dmed$signif, fdbms_mant_tl$signif, fdbms_mant_te$signif, fdbms_mant_p$signif, fdbms_mant_md$signif, fdbms_mant_la$signif)
fdbms__mant_pv <- fdbms__mant_pv[,1]
fdbms__mant_pv <- p.adjust(fdbms__mant_pv, method = "bonferroni")
fdbms__mant_pv
```

    ##  [1] 1.000 1.000 1.000 0.957 1.000 1.000 1.000 1.000 1.000 1.000 1.000

``` r
# Ocean sites and mixed lakes
envbom_dist_vc <- dist(scaled_env[ocean_mixed_sites_geo,c(5)], method = "euclidean")
fdbom_mant_vc <- mantel(fdbom_dist, envbom_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_vc, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.05755 
    ##       Significance: 0.703 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.273 0.324 0.350 0.388 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_v <- dist(scaled_env[ocean_mixed_sites_geo,c(6)], method = "euclidean")
fdbom_mant_v <- mantel(fdbom_dist, envbom_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_v, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.06941 
    ##       Significance: 0.723 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.261 0.326 0.361 0.387 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_sa <- dist(scaled_env[ocean_mixed_sites_geo,c(7)], method = "euclidean")
fdbom_mant_sa <- mantel(fdbom_dist, envbom_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_sa, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.03054 
    ##       Significance: 0.659 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.276 0.331 0.367 0.402 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_dmin <- dist(scaled_env[ocean_mixed_sites_geo,c(8)], method = "euclidean")
fdbom_mant_dmin <- mantel(fdbom_dist, envbom_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_dmin, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1538 
    ##       Significance: 0.23 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.253 0.329 0.386 0.419 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_dmean <- dist(scaled_env[ocean_mixed_sites_geo,c(9)], method = "euclidean")
fdbom_mant_dmean <- mantel(fdbom_dist, envbom_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_dmean, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1799 
    ##       Significance: 0.241 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.251 0.286 0.321 0.355 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_dmed <- dist(scaled_env[ocean_mixed_sites_geo,c(10)], method = "euclidean")
fdbom_mant_dmed <- mantel(fdbom_dist, envbom_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_dmed, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1365 
    ##       Significance: 0.326 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.255 0.295 0.334 0.373 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_tl <- dist(scaled_env[ocean_mixed_sites_geo,c(11)], method = "euclidean")
fdbom_mant_tl <- mantel(fdbom_dist, envbom_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_tl, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1895 
    ##       Significance: 0.14 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.232 0.292 0.342 0.383 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_te <- dist(scaled_env[ocean_mixed_sites_geo,c(12)], method = "euclidean")
fdbom_mant_te <- mantel(fdbom_dist, envbom_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_te, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1405 
    ##       Significance: 0.26 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.244 0.295 0.329 0.390 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_p <- dist(scaled_env[ocean_mixed_sites_geo,c(13)], method = "euclidean")
fdbom_mant_p <- mantel(fdbom_dist, envbom_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.05181 
    ##       Significance: 0.682 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.282 0.322 0.356 0.382 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_md <- dist(scaled_env[ocean_mixed_sites_geo,c(14)], method = "euclidean")
fdbom_mant_md <- mantel(fdbom_dist, envbom_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_md, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.06813 
    ##       Significance: 0.639 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.178 0.236 0.294 0.373 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_la <- dist(scaled_env[ocean_mixed_sites_geo,c(15)], method = "euclidean")
fdbom_mant_la <- mantel(fdbom_dist, envbom_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
fdbom_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist_la, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.00669 
    ##       Significance: 0.537 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.276 0.363 0.415 0.465 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
fdbom__mant_pv <- rbind(fdbom_mant_vc$signif, fdbom_mant_v$signif, fdbom_mant_sa$signif, fdbom_mant_dmin$signif, fdbom_mant_dmean$signif, fdbom_mant_dmed$signif, fdbom_mant_tl$signif, fdbom_mant_te$signif, fdbom_mant_p$signif, fdbom_mant_md$signif, fdbom_mant_la$signif)
fdbom__mant_pv <- fdbom__mant_pv[,1]
fdbom__mant_pv <- p.adjust(fdbom__mant_pv, method = "bonferroni")
fdbom__mant_pv
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1

``` r
# Stratified lakes and ocean sites
envbso_dist_vc <- dist(scaled_env[ocean_stratified_sites,c(5)], method = "euclidean")
fdbso_mant_vc <- mantel(fdbso_dist, envbso_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_vc, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1871 
    ##       Significance: 0.856 
    ## 
    ## Upper quantiles of permutations (null model):
    ##      90%      95%    97.5%      99% 
    ## -0.05712 -0.03177 -0.00422  0.02222 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_v <- dist(scaled_env[ocean_stratified_sites,c(6)], method = "euclidean")
fdbso_mant_v <- mantel(fdbso_dist, envbso_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_v, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1723 
    ##       Significance: 0.82 
    ## 
    ## Upper quantiles of permutations (null model):
    ##      90%      95%    97.5%      99% 
    ## -0.04855 -0.02637  0.00691  0.03803 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_sa <- dist(scaled_env[ocean_stratified_sites,c(7)], method = "euclidean")
fdbso_mant_sa <- mantel(fdbso_dist, envbso_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_sa, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1542 
    ##       Significance: 0.778 
    ## 
    ## Upper quantiles of permutations (null model):
    ##      90%      95%    97.5%      99% 
    ## -0.03922 -0.00192  0.02420  0.04666 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_dmin <- dist(scaled_env[ocean_stratified_sites,c(8)], method = "euclidean")
fdbso_mant_dmin <- mantel(fdbso_dist, envbso_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_dmin, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6621 
    ##       Significance: 0.052 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.600 0.662 0.714 0.744 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_dmean <- dist(scaled_env[ocean_stratified_sites,c(9)], method = "euclidean")
fdbso_mant_dmean <- mantel(fdbso_dist, envbso_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_dmean, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5043 
    ##       Significance: 0.203 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.573 0.616 0.666 0.702 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_dmed <- dist(scaled_env[ocean_stratified_sites,c(10)], method = "euclidean")
fdbso_mant_dmed <- mantel(fdbso_dist, envbso_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_dmed, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4905 
    ##       Significance: 0.216 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.566 0.620 0.660 0.703 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_tl <- dist(scaled_env[ocean_stratified_sites,c(11)], method = "euclidean")
fdbso_mant_tl <- mantel(fdbso_dist, envbso_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_tl, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3578 
    ##       Significance: 0.849 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.575 0.606 0.640 0.668 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_te <- dist(scaled_env[ocean_stratified_sites,c(12)], method = "euclidean")
fdbso_mant_te <- mantel(fdbso_dist, envbso_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_te, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.517 
    ##       Significance: 0.239 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.577 0.613 0.643 0.672 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_p <- dist(scaled_env[ocean_stratified_sites,c(13)], method = "euclidean")
fdbso_mant_p <- mantel(fdbso_dist, envbso_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1936 
    ##       Significance: 0.948 
    ## 
    ## Upper quantiles of permutations (null model):
    ##     90%     95%   97.5%     99% 
    ## -0.0867 -0.0723 -0.0606 -0.0411 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_md <- dist(scaled_env[ocean_stratified_sites,c(14)], method = "euclidean")
fdbso_mant_md <- mantel(fdbso_dist, envbso_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_md, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.03981 
    ##       Significance: 0.713 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.171 0.235 0.289 0.329 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_la <- dist(scaled_env[ocean_stratified_sites,c(15)], method = "euclidean")
fdbso_mant_la <- mantel(fdbso_dist, envbso_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
fdbso_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist_la, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1062 
    ##       Significance: 0.728 
    ## 
    ## Upper quantiles of permutations (null model):
    ##     90%     95%   97.5%     99% 
    ## -0.0042  0.0224  0.0439  0.0654 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
fdbso__mant_pv <- rbind(fdbso_mant_vc$signif, fdbso_mant_v$signif, fdbso_mant_sa$signif, fdbso_mant_dmin$signif, fdbso_mant_dmean$signif, fdbso_mant_dmed$signif, fdbso_mant_tl$signif, fdbso_mant_te$signif, fdbso_mant_p$signif, fdbso_mant_md$signif, fdbso_mant_la$signif)
fdbso__mant_pv <- fdbso__mant_pv[,1]
fdbso__mant_pv <- p.adjust(fdbso__mant_pv, method = "bonferroni")
fdbso__mant_pv
```

    ##  [1] 1.000 1.000 1.000 0.572 1.000 1.000 1.000 1.000 1.000 1.000 1.000

## Plot trait NMDS and envfit results

- Includes a plot of correlated environmental variables

``` r
fd_NMDS_data.scores <- as.data.frame(scores(fd_NMDS))
fd_NMDS_data.scores$Stratification <- env[surveyed_sites,19]
fd_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
fd_NMDS_data.scores$Stratification <- factor(fd_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

fde_ef_coord_cont <- as.data.frame(scores(fde_ef, "vectors")) * ordiArrowMul(fde_ef)
fde_row_names <- c("S")
# Assign the new row names to the data frame
fde_ef_coord_cont <- data.frame(row.names = fde_row_names, fde_ef_coord_cont)

# fdb_ef_coord_cont <- as.data.frame(scores(fdb_ef, "vectors")) * ordiArrowMul(fdb_ef)
# fdb_row_names <- c("minD")
# # Assign the new row names to the data frame
# fdb_ef_coord_cont <- data.frame(row.names = fdb_row_names, fdb_ef_coord_cont)

fd_ef_plot <- ggplot(data = fd_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type='t', level = 0.95) +
  geom_point(data = fd_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  geom_text_repel(data = fd_NMDS_data.scores, label = fd_NMDS_data.scores$Lakes, size = 5, point.padding = 5, max.overlaps = 30, nudge_x = -0.01) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = fde_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = fde_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
            color = env_cont, label = row.names(fde_ef_coord_cont), size = 5) +
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              data = fdb_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#698B22") +
  # geom_text(data = fdb_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
  #           color = "#698B22", label = row.names(fdb_ef_coord_cont), size = 7) + 
  theme(text = element_text(size = 22), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16), panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  annotate("text", x = 0.55, y = 0.57, size = 5,
            label = paste("Stress: ", round(fd_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ")
fd_ef_plot <- fd_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
fd_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20trait%20NMDS%20and%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_ef_plot.png", fd_ef_plot, width = 6, height = 4, units = "in")

# For 90 percent confidence interval: x = 0.35, y = 0.55
# For combo plot: x = 0.18, y = 0.48
# For num plot: x = 0.29, y = 0.48

# Determine outliers
ordered_fd_NMDS1_data.scores <- fd_NMDS_data.scores[order(fd_NMDS_data.scores$NMDS1), ]
Q1 <- ordered_fd_NMDS1_data.scores[6,1]
Q3 <- ordered_fd_NMDS1_data.scores[18,1]
IQR1 <- IQR(fd_NMDS_data.scores$NMDS1)
Q1 - 1.5*IQR1
```

    ## [1] -0.3259683

``` r
Q3 + 1.5*IQR1
```

    ## [1] 0.4196323

``` r
ordered_fd_NMDS1_data.scores$NMDS1
```

    ##  [1] -0.54471662 -0.46676575 -0.41243841 -0.38993704 -0.27309338 -0.05674459
    ##  [7]  0.04251538  0.10179296  0.11276146  0.11723156  0.11863443  0.11982642
    ## [13]  0.13385200  0.13636556  0.14594387  0.14598171  0.14807661  0.15040854
    ## [19]  0.15218233  0.16742330  0.16829006  0.18240963

``` r
ordered_fd_NMDS2_data.scores <- fd_NMDS_data.scores[order(fd_NMDS_data.scores$NMDS2), ]
Q1 <- ordered_fd_NMDS2_data.scores[6,2]
Q3 <- ordered_fd_NMDS2_data.scores[18,2]
IQR2 <- IQR(fd_NMDS_data.scores$NMDS2)
Q1 - 1.5*IQR2
```

    ## [1] -0.2096355

``` r
Q3 + 1.5*IQR2
```

    ## [1] 0.2270116

``` r
ordered_fd_NMDS2_data.scores$NMDS2
```

    ##  [1] -0.262760710 -0.179342240 -0.160961640 -0.153585891 -0.068051183
    ##  [6] -0.052511068 -0.047965983 -0.037359886 -0.016872227 -0.012846404
    ## [11] -0.009261649 -0.004607925 -0.004291215  0.003673066  0.043323317
    ## [16]  0.044807323  0.056230627  0.069887221  0.076480461  0.166804642
    ## [21]  0.170593494  0.378617872

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_NMDS_boxplot.tiff", width = 1200, height = 800, res = 150, type = "cairo")

fd_dissimilarity_matrix <- as.matrix(fd_dist)

# Set the diagonal elements to NA
diag(fd_dissimilarity_matrix) <- NA

# Get the labels of the sites
sites <- attr(fd_dissimilarity_matrix, "Labels")

# Calculate the mean dissimilarity for each group
fd_mean_dissimilarity <- aggregate(fd_dissimilarity_matrix, by = list(stratification_group), FUN = mean, na.rm = TRUE)

row.names(fd_mean_dissimilarity) <- fd_mean_dissimilarity$Group.1

fd_mean_dissimilarity <- fd_mean_dissimilarity[,-1] 

fd_mean_dissimilarity <- as.data.frame(t(fd_mean_dissimilarity))

fd_mean_dissimilarity$Stratification <- env[surveyed_sites,19]

boxplot(Ocean ~ Stratification, fd_mean_dissimilarity[c(7,9,11,16,18,19),])
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
fd_mean_dissimilarity$Stratification <- factor(fd_mean_dissimilarity$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

fd_ocn_dist_violin_plot <- ggplot(fd_mean_dissimilarity, mapping = aes(x= Stratification, y= Ocean,  color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.8,
            width = 0.1) +
  geom_text_repel(data = fd_mean_dissimilarity, label = row.names(fd_mean_dissimilarity), size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  ylab("Ocean site distances") +
  xlab("Site type") +
  labs(colour = "Site type:  ", fill = "Site type:  ")
fd_ocn_dist_violin_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20trait%20NMDS%20and%20envfit%20results-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_ocn_dist_violin_plot.png", fd_ocn_dist_violin_plot, width = 8, height = 4, units = "in")

fd_mix_dist_violin_plot <- ggplot(fd_mean_dissimilarity, mapping = aes(x= Stratification, y= Mixed,  color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.8,
            width = 0.1) +
  geom_text_repel(data = fd_mean_dissimilarity, label = row.names(fd_mean_dissimilarity), size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  ylab("Mixed lake distances") +
  xlab("Site type") +
  labs(colour = "Site type:  ", fill = "Site type:  ")
fd_mix_dist_violin_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20trait%20NMDS%20and%20envfit%20results-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_mix_dist_violin_plot.png", fd_mix_dist_violin_plot, width = 8, height = 4, units = "in")

fd_strat_dist_violin_plot <- ggplot(fd_mean_dissimilarity, mapping = aes(x= Stratification, y= Stratified, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.8,
            width = 0.1) +
  geom_text_repel(data = fd_mean_dissimilarity, label = row.names(fd_mean_dissimilarity), size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  ylab("Stratified lake distances") +
  xlab("Site type") +
  labs(colour = "Site type:  ", fill = "Site type:  ")
fd_strat_dist_violin_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20trait%20NMDS%20and%20envfit%20results-4.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_strat_dist_violin_plot.png", fd_strat_dist_violin_plot, width = 8, height = 4, units = "in")
```

## Boxplot of trait intradissimilarity heterogeneity

- The dissimilarity between sites that share the same stratification.

``` r
fd_bd_dist <- stratification_fd_bd$distances
fd_bd_dist <- as.data.frame(fd_bd_dist)
fd_bd_dist$X <- row.names(fd_bd_dist)
fd_bd_dist_env <- merge(fd_bd_dist, env[surveyed_sites,], by = "X", sort = F)

fd_bd_dist_env$Stratification <- factor(fd_bd_dist_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

fd_bd_dist_env_plot <- ggplot(fd_bd_dist_env, aes(x = Stratification, y = fd_bd_dist, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.8,
            width = 0.1) +
    geom_text_repel(data = fd_bd_dist_env, label = fd_bd_dist_env$X, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black")) +
  xlab("Site type") +
  ylab("Distance to Centroid") +
  labs(color = "Stratification", tag = "B")
fd_bd_dist_env_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Boxplot%20of%20trait%20intradissimilarity%20heterogeneity-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_bd_dist_env_plot.png", fd_bd_dist_env_plot, width = 8.25, height = 4.13, units = "in")
```

## FDis Dendrogram

``` r
# cluster communities using average-linkage algorithm
fd_dist_clust <- hclust(fd_dist, method = "average")

# Open a PNG device
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_dendro.tiff", width = 1200, height = 800, res = 150, type = "cairo")

# Create your plot using the plot() function
plot(fd_dist_clust, 
     xlab = "Sites",
     ylab = "Functional dissimilarity",
     main = "",
     sub = "",
     pch = 20,
     col= "black")

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_si.tiff", width = 1200, height = 800, res = 150, type = "cairo")

Si <- numeric(nrow(presabs_lake[surveyed_sites,]))
for (k in 2:(nrow(presabs_lake[surveyed_sites,])-1))
{
  sil<-silhouette(cutree(fd_dist_clust, k=k), fd_dist)
  Si[k] <- summary(sil)$avg.width
}
k.best<-which.max(Si)
plot(1:nrow(presabs_lake[surveyed_sites,]), Si, type = "h", main = "Silhouette", xlab = "K", ylab = "Width")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col = "red", font = 2, col.axis = "red")

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_pam.tiff", width = 1200, height = 800, res = 150, type = "cairo")

fd_pam <- pam(fd_dist,2)
plot(fd_pam)

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Phylogeny dissimilarity distances

- We use unifrac to calculate distances in the phylogeny.

``` r
### Regular
pd_dist <- unifrac(presabs_lake[surveyed_sites,], stree)
# Mixed and stratified lakes
pdms_dist <- unifrac(presabs_lake[mixed_stratified_lakes,], stree)
# Ocean sites and mixed lakes
pdom_dist <- unifrac(presabs_lake[ocean_mixed_sites,], stree)
# Stratified lakes and ocean sites
pdso_dist <- unifrac(presabs_lake[ocean_stratified_sites,], stree)
# Stratified lakes
pds_dist <- unifrac(presabs_lake[stratified_lakes,], stree)

### Environmental
pde_dist <- unifrac(presabs_lake[surveyed_sites_env,], etree)
# Mixed and stratified lakes
pdems_dist <- unifrac(presabs_lake[mixed_stratified_lakes,], etree)
# Ocean sites and mixed lakes
pdeom_dist <- unifrac(presabs_lake[ocean_mixed_sites_env,], etree)
# Stratified lakes and ocean sites
pdeso_dist <- unifrac(presabs_lake[ocean_stratified_sites_env,], etree)


### Geographic
pdb_dist <- unifrac(presabs_lake[surveyed_sites_geo,], btree)
# Mixed and stratified lakes
pdbms_dist <- unifrac(presabs_lake[mixed_stratified_lakes_geo,], btree)
# Ocean sites and mixed lakes
pdbom_dist <- unifrac(presabs_lake[ocean_mixed_sites_geo,], btree)
# Stratified lakes and ocean sites
pdbso_dist <- unifrac(presabs_lake[ocean_stratified_sites,], btree)
```

## Phylogeny NMDS

- To constrain dissimilarities we perform Nonmetric Multidimensional
  Scaling (NMDS), which tries to find a stable solution using the
  metaMDS package.

``` r
### Regular
pd_NMDS <- metaMDS(pd_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
round(pd_NMDS$stress, digits = 2)
# Mixed and stratified lakes
pdms_NMDS <- metaMDS(pdms_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
pdom_NMDS <- metaMDS(pdom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
pdso_NMDS <- metaMDS(pdso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes
pds_NMDS <- metaMDS(pds_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)


### Environmental
pde_NMDS <- metaMDS(pde_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
pdems_NMDS <- metaMDS(pdems_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
pdeom_NMDS <- metaMDS(pdeom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
pdeso_NMDS <- metaMDS(pdeso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(pdeso_dist, try = 1000, parallel = 4, previous.best, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
pdb_NMDS <- metaMDS(pdb_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
pdbms_NMDS <- metaMDS(pdbms_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
pdbom_NMDS <- metaMDS(pdbom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
pdbso_NMDS <- metaMDS(pdbso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
```

## Determine stratification homogeneity using phylogenetic distances

- We use betadisper to determine homogeneity of the phylogenetic
  distances based on the stratification category. Is the dispersion of
  phylogenetic distances similar within stratification categories?

``` r
### Stratification
pd_bd <- betadisper(pd_dist, stratification_group)
anova(pd_bd)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value  Pr(>F)  
    ## Groups     2 0.041497 0.0207485  2.8718 0.08133 .
    ## Residuals 19 0.137271 0.0072248                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
permutest(pd_bd, pairwise = TRUE)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
    ## Groups     2 0.041497 0.0207485 2.8718    999  0.071 .
    ## Residuals 19 0.137271 0.0072248                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##               Mixed    Ocean Stratified
    ## Mixed               0.799000      0.048
    ## Ocean      0.786333               0.112
    ## Stratified 0.045108 0.120496

``` r
(pd_bd.HSD <- TukeyHSD(pd_bd))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff        lwr        upr     p adj
    ## Ocean-Mixed       0.01013397 -0.1064845 0.12675242 0.9735283
    ## Stratified-Mixed -0.08555708 -0.1935248 0.02241062 0.1362184
    ## Stratified-Ocean -0.09569106 -0.2123095 0.02092739 0.1198251

``` r
boxplot(pd_bd)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Determine%20stratification%20homogeneity%20using%20phylogenetic%20distances-1.png)<!-- -->

## Determine stratification dissimilarity using phylogenetic distances

- We use adonis to determine if dissimilarities of species phylogenetic
  distances by stratification categories are significant and how much of
  the variation is explained by phylogenetic dissimilarities.

``` r
### Stratification
pd_pms <- adonis2(pd_dist ~ env[surveyed_sites,19], permutations = 999)
pd_pms
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = pd_dist ~ env[surveyed_sites, 19], permutations = 999)
    ##                         Df SumOfSqs     R2     F Pr(>F)    
    ## env[surveyed_sites, 19]  2   1.8501 0.3911 6.102  0.001 ***
    ## Residual                19   2.8804 0.6089                 
    ## Total                   21   4.7306 1.0000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratification
pd_pms_pair <- pairwise.adonis(pd_dist, env[surveyed_sites,19], p.adjust.m = "bonferroni", perm = 999)
pd_pms_pair
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.2940870 9.415922 0.4021162   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.1734121 8.199716 0.4059322   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.2549888 1.443665 0.1073863   0.105      0.315

## Envfit environmental influence on phylogenetic distances

- We use envfit to determine significantly correlated environmental
  variables to our phylogenetic NMDS. We will use these results, if
  significant, in our figure.

``` r
### Environmental
# All sites for figure
pde_ef <- envfit(pde_NMDS, env[surveyed_sites_env,c(34,36)], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
pde_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##       NMDS1    NMDS2     r2 Pr(>r)    
    ## S   0.96421  0.26513 0.7072  0.008 ** 
    ## pH  0.99017 -0.13990 0.8530  0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
pde_efp <- p.adjust.envfit(pde_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pde_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##       NMDS1    NMDS2     r2 Pr(>r)   
    ## S   0.96421  0.26513 0.7072  0.016 * 
    ## pH  0.99017 -0.13990 0.8530  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# All sites by stratification
pdea_ef <- envfit(pde_NMDS, env[surveyed_sites_env,c(2,6,8,10)], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
pdea_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)   
    ## temperature_median -0.73728 -0.67559 0.0772  0.645   
    ## salinity_median     0.96421  0.26513 0.7072  0.011 * 
    ## oxygen_median       0.56896 -0.82236 0.6630  0.168   
    ## pH_median           0.99017 -0.13990 0.8530  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
pdea_efp <- p.adjust.envfit(pdea_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdea_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)   
    ## temperature_median -0.73728 -0.67559 0.0772  1.000   
    ## salinity_median     0.96421  0.26513 0.7072  0.044 * 
    ## oxygen_median       0.56896 -0.82236 0.6630  0.672   
    ## pH_median           0.99017 -0.13990 0.8530  0.008 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
pdems_ef <- envfit(pdems_NMDS, env[mixed_stratified_lakes,c(2,6,8,10)], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdems_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)   
    ## temperature_median -0.89796  0.44007 0.0902  0.725   
    ## salinity_median     0.86015  0.51004 0.6978  0.016 * 
    ## oxygen_median       0.57121 -0.82080 0.4948  0.527   
    ## pH_median           0.99082  0.13521 0.8248  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
pdems_efp <- p.adjust.envfit(pdems_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdems_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)   
    ## temperature_median -0.89796  0.44007 0.0902  1.000   
    ## salinity_median     0.86015  0.51004 0.6978  0.064 . 
    ## oxygen_median       0.57121 -0.82080 0.4948  1.000   
    ## pH_median           0.99082  0.13521 0.8248  0.008 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
pdeom_ef <- envfit(pdeom_NMDS, env[ocean_mixed_sites_env,c(2,6,8,10)], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
pdeom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)   
    ## temperature_median -0.47339 -0.88085 0.0251  0.872   
    ## salinity_median    -0.85759 -0.51433 0.7435  0.008 **
    ## oxygen_median      -0.92981 -0.36805 0.7949  0.011 * 
    ## pH_median          -0.95681  0.29070 0.7286  0.008 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
pdeom_efp <- p.adjust.envfit(pdeom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdeom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.47339 -0.88085 0.0251  1.000  
    ## salinity_median    -0.85759 -0.51433 0.7435  0.032 *
    ## oxygen_median      -0.92981 -0.36805 0.7949  0.044 *
    ## pH_median          -0.95681  0.29070 0.7286  0.032 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
pdeso_ef <- envfit(pdeso_NMDS, env[ocean_stratified_sites_env,c(2,6,8,10)], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
pdeso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median  0.0034001  0.9999900 0.0308  0.918
    ## salinity_median    -0.0028945  1.0000000 0.5080  0.421
    ## oxygen_median      -0.0054402  0.9999900 0.7214  0.553
    ## pH_median          -0.0030174  1.0000000 0.7008  0.283
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
pdeso_efp <- p.adjust.envfit(pdeso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdeso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median  0.0034001  0.9999900 0.0308      1
    ## salinity_median    -0.0028945  1.0000000 0.5080      1
    ## oxygen_median      -0.0054402  0.9999900 0.7214      1
    ## pH_median          -0.0030174  1.0000000 0.7008      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes
pds_ef <- envfit(pds_NMDS, env[stratified_lakes,c(8)], permutations = 999, na.rm = TRUE)
pds_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##       NMDS1  NMDS2    r2 Pr(>r)  
    ## [1,] 0.8374 0.5466 0.745  0.036 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
pds_efp <- p.adjust.envfit(pds_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pds_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##       NMDS1  NMDS2    r2 Pr(>r)  
    ## [1,] 0.8374 0.5466 0.745  0.036 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# All sites for figure
pdb_ef <- envfit(pdb_NMDS, env[surveyed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdb_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.30682 -0.95177 0.2744 0.13786  
    ## volume_m3                   0.29211 -0.95639 0.2676 0.14186  
    ## surface_area_m2             0.24076 -0.97058 0.3748 0.17183  
    ## distance_to_ocean_min_m    -0.92173  0.38784 0.5897 0.05594 .
    ## distance_to_ocean_mean_m   -0.76142  0.64826 0.5393 0.46154  
    ## distance_to_ocean_median_m -0.76966  0.63845 0.5307 0.46953  
    ## tidal_lag_time_minutes     -0.63892  0.76927 0.6134 0.31968  
    ## tidal_efficiency            0.63787 -0.77014 0.6174 0.22178  
    ## perimeter_fromSat           0.31917 -0.94770 0.2856 0.40360  
    ## max_depth                  -0.30928 -0.95097 0.0696 0.95504  
    ## logArea                     0.24301 -0.97002 0.1776 0.40260  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
pdb_efp <- p.adjust.envfit(pdb_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdb_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.30682 -0.95177 0.2744 1.0000
    ## volume_m3                   0.29211 -0.95639 0.2676 1.0000
    ## surface_area_m2             0.24076 -0.97058 0.3748 1.0000
    ## distance_to_ocean_min_m    -0.92173  0.38784 0.5897 0.6154
    ## distance_to_ocean_mean_m   -0.76142  0.64826 0.5393 1.0000
    ## distance_to_ocean_median_m -0.76966  0.63845 0.5307 1.0000
    ## tidal_lag_time_minutes     -0.63892  0.76927 0.6134 1.0000
    ## tidal_efficiency            0.63787 -0.77014 0.6174 1.0000
    ## perimeter_fromSat           0.31917 -0.94770 0.2856 1.0000
    ## max_depth                  -0.30928 -0.95097 0.0696 1.0000
    ## logArea                     0.24301 -0.97002 0.1776 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
pdba_ef <- envfit(pdb_NMDS, env[surveyed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.30682 -0.95177 0.2744 0.14885  
    ## volume_m3                   0.29211 -0.95639 0.2676 0.15185  
    ## surface_area_m2             0.24076 -0.97058 0.3748 0.18482  
    ## distance_to_ocean_min_m    -0.92173  0.38784 0.5897 0.04595 *
    ## distance_to_ocean_mean_m   -0.76142  0.64826 0.5393 0.49850  
    ## distance_to_ocean_median_m -0.76966  0.63845 0.5307 0.51149  
    ## tidal_lag_time_minutes     -0.63892  0.76927 0.6134 0.32268  
    ## tidal_efficiency            0.63787 -0.77014 0.6174 0.24476  
    ## perimeter_fromSat           0.31917 -0.94770 0.2856 0.40260  
    ## max_depth                  -0.30928 -0.95097 0.0696 0.96503  
    ## logArea                     0.24301 -0.97002 0.1776 0.38062  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
pdba_efp <- p.adjust.envfit(pdba_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdba_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.30682 -0.95177 0.2744 1.0000
    ## volume_m3                   0.29211 -0.95639 0.2676 1.0000
    ## surface_area_m2             0.24076 -0.97058 0.3748 1.0000
    ## distance_to_ocean_min_m    -0.92173  0.38784 0.5897 0.5055
    ## distance_to_ocean_mean_m   -0.76142  0.64826 0.5393 1.0000
    ## distance_to_ocean_median_m -0.76966  0.63845 0.5307 1.0000
    ## tidal_lag_time_minutes     -0.63892  0.76927 0.6134 1.0000
    ## tidal_efficiency            0.63787 -0.77014 0.6174 1.0000
    ## perimeter_fromSat           0.31917 -0.94770 0.2856 1.0000
    ## max_depth                  -0.30928 -0.95097 0.0696 1.0000
    ## logArea                     0.24301 -0.97002 0.1776 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
pdbms_ef <- envfit(pdbms_NMDS, env[mixed_stratified_lakes_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline     -0.05752  0.99834 0.0711 0.81419  
    ## volume_m3                  -0.49296  0.87005 0.1314 0.97203  
    ## surface_area_m2            -0.34561  0.93838 0.1047 0.95305  
    ## distance_to_ocean_min_m    -0.83616  0.54849 0.5156 0.08092 .
    ## distance_to_ocean_mean_m   -0.58972  0.80761 0.4501 0.36464  
    ## distance_to_ocean_median_m -0.57832  0.81581 0.4497 0.33367  
    ## tidal_lag_time_minutes     -0.39299  0.91954 0.6050 0.23676  
    ## tidal_efficiency            0.40044 -0.91632 0.6012 0.16783  
    ## perimeter_fromSat          -0.09436  0.99554 0.0784 0.80420  
    ## max_depth                  -0.55340 -0.83292 0.1050 0.95904  
    ## logArea                    -0.50441  0.86346 0.0064 0.99500  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
pdbms_efp <- p.adjust.envfit(pdbms_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdbms_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.05752  0.99834 0.0711 1.0000
    ## volume_m3                  -0.49296  0.87005 0.1314 1.0000
    ## surface_area_m2            -0.34561  0.93838 0.1047 1.0000
    ## distance_to_ocean_min_m    -0.83616  0.54849 0.5156 0.8901
    ## distance_to_ocean_mean_m   -0.58972  0.80761 0.4501 1.0000
    ## distance_to_ocean_median_m -0.57832  0.81581 0.4497 1.0000
    ## tidal_lag_time_minutes     -0.39299  0.91954 0.6050 1.0000
    ## tidal_efficiency            0.40044 -0.91632 0.6012 1.0000
    ## perimeter_fromSat          -0.09436  0.99554 0.0784 1.0000
    ## max_depth                  -0.55340 -0.83292 0.1050 1.0000
    ## logArea                    -0.50441  0.86346 0.0064 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
pdbom_ef <- envfit(pdbom_NMDS, env[ocean_mixed_sites_geo,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2   Pr(>r)   
    ## volume_m3_w_chemocline     -0.95041  0.31100 0.2062 0.492507   
    ## volume_m3                  -0.95024  0.31152 0.2057 0.493506   
    ## surface_area_m2            -0.51540  0.85695 0.2737 0.570430   
    ## distance_to_ocean_min_m     0.56663 -0.82397 0.4592 0.138861   
    ## distance_to_ocean_mean_m    0.29666 -0.95498 0.4458 0.204795   
    ## distance_to_ocean_median_m  0.29140 -0.95660 0.4366 0.226773   
    ## tidal_lag_time_minutes      0.34048 -0.94025 0.5044 0.024975 * 
    ## tidal_efficiency           -0.34008  0.94040 0.4796 0.041958 * 
    ## perimeter_fromSat          -0.74671  0.66515 0.2015 0.623377   
    ## max_depth                  -0.87692 -0.48064 0.6324 0.004995 **
    ## logArea                    -0.72345  0.69038 0.2509 0.475524   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
pdbom_efp <- p.adjust.envfit(pdbom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdbom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline     -0.95041  0.31100 0.2062 1.00000  
    ## volume_m3                  -0.95024  0.31152 0.2057 1.00000  
    ## surface_area_m2            -0.51540  0.85695 0.2737 1.00000  
    ## distance_to_ocean_min_m     0.56663 -0.82397 0.4592 1.00000  
    ## distance_to_ocean_mean_m    0.29666 -0.95498 0.4458 1.00000  
    ## distance_to_ocean_median_m  0.29140 -0.95660 0.4366 1.00000  
    ## tidal_lag_time_minutes      0.34048 -0.94025 0.5044 0.27473  
    ## tidal_efficiency           -0.34008  0.94040 0.4796 0.46154  
    ## perimeter_fromSat          -0.74671  0.66515 0.2015 1.00000  
    ## max_depth                  -0.87692 -0.48064 0.6324 0.05495 .
    ## logArea                    -0.72345  0.69038 0.2509 1.00000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
pdbso_ef <- envfit(pdbso_NMDS, env[ocean_stratified_sites,c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline     -0.32433  0.94594 0.4850 0.03097 *
    ## volume_m3                  -0.31779  0.94816 0.4742 0.04496 *
    ## surface_area_m2            -0.30237  0.95319 0.6099 0.03796 *
    ## distance_to_ocean_min_m     0.99834  0.05760 0.6213 0.25375  
    ## distance_to_ocean_mean_m    0.89864  0.43868 0.6226 0.58042  
    ## distance_to_ocean_median_m  0.87457  0.48490 0.6154 0.57742  
    ## tidal_lag_time_minutes      0.99917  0.04077 0.7125 0.92208  
    ## tidal_efficiency           -0.93913 -0.34357 0.7835 0.58941  
    ## perimeter_fromSat          -0.38917  0.92116 0.5288 0.13387  
    ## max_depth                   0.16212  0.98677 0.1849 0.57642  
    ## logArea                    -0.26902  0.96313 0.4195 0.13387  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
pdbso_efp <- p.adjust.envfit(pdbso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdbso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.32433  0.94594 0.4850 0.3407
    ## volume_m3                  -0.31779  0.94816 0.4742 0.4945
    ## surface_area_m2            -0.30237  0.95319 0.6099 0.4176
    ## distance_to_ocean_min_m     0.99834  0.05760 0.6213 1.0000
    ## distance_to_ocean_mean_m    0.89864  0.43868 0.6226 1.0000
    ## distance_to_ocean_median_m  0.87457  0.48490 0.6154 1.0000
    ## tidal_lag_time_minutes      0.99917  0.04077 0.7125 1.0000
    ## tidal_efficiency           -0.93913 -0.34357 0.7835 1.0000
    ## perimeter_fromSat          -0.38917  0.92116 0.5288 1.0000
    ## max_depth                   0.16212  0.98677 0.1849 1.0000
    ## logArea                    -0.26902  0.96313 0.4195 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Traits
# All sites for figure
pdt_ef <- envfit(pd_NMDS, FD_total_env[surveyed_sites,c(56,62)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites,19])
pdt_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##      NMDS1    NMDS2    r2   Pr(>r)   
    ## T -0.72127  0.69265 0.704 0.004995 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##        NMDS1   NMDS2
    ## EC1b -0.3955 -0.0179
    ## EC2g -0.4873 -0.0630
    ## EC3n  0.0721  0.0076
    ## 
    ## Goodness of fit:
    ##       r2   Pr(>r)   
    ## EC 0.338 0.004995 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
pdt_efp <- p.adjust.envfit(pdt_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdt_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##      NMDS1    NMDS2    r2  Pr(>r)   
    ## T -0.72127  0.69265 0.704 0.00999 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##        NMDS1   NMDS2
    ## EC1b -0.3955 -0.0179
    ## EC2g -0.4873 -0.0630
    ## EC3n  0.0721  0.0076
    ## 
    ## Goodness of fit:
    ##       r2  Pr(>r)   
    ## EC 0.338 0.00999 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
pdta_ef <- envfit(pd_NMDS, FD_total_env[surveyed_sites,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.21079  0.97753 0.2133 0.044955 * 
    ## Troph            -0.72127  0.69265 0.7040 0.003996 **
    ## DepthMin         -0.09925 -0.99506 0.0523 0.122877   
    ## DepthMax          0.72863  0.68490 0.1969 0.160839   
    ## TempPrefMin       0.56071 -0.82801 0.2013 0.398601   
    ## TempPrefMax       0.51263  0.85861 0.0027 0.909091   
    ## DorsalSpinesMean  0.97411 -0.22609 0.5280 0.032967 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.0491 -0.0887
    ## BodyShapeI3f         0.1007  0.0309
    ## BodyShapeI4e        -0.4143 -0.0493
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -0.3899 -0.0284
    ## OperculumPresentyes  0.0866  0.0063
    ## FeedingPathb         0.0382  0.0027
    ## FeedingPathp        -0.3816 -0.0269
    ## RepGuild11b         -0.3955 -0.0179
    ## RepGuild12g         -0.4873 -0.0630
    ## RepGuild13n          0.0721  0.0076
    ## RepGuild22eb        -0.3955 -0.0179
    ## RepGuild23n         -0.4873 -0.0630
    ## RepGuild26s          0.0721  0.0076
    ## ParentalCare3p      -0.3670 -0.0147
    ## ParentalCare4n       0.1713  0.0068
    ## WaterPref1s          0.1177  0.0096
    ## WaterPref3a         -0.4002 -0.0327
    ## 
    ## Goodness of fit:
    ##                      r2   Pr(>r)   
    ## BodyShapeI       0.4083 0.059940 . 
    ## DemersPelag      0.0000 1.000000   
    ## OperculumPresent 0.3415 0.150849   
    ## FeedingPath      0.1471 0.362637   
    ## RepGuild1        0.3380 0.006993 **
    ## RepGuild2        0.3380 0.006993 **
    ## ParentalCare     0.6332 0.238761   
    ## WaterPref        0.4770 0.063936 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
pdta_efp <- p.adjust.envfit(pdta_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdta_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL       0.21079  0.97753 0.2133 0.67433  
    ## Troph            -0.72127  0.69265 0.7040 0.05994 .
    ## DepthMin         -0.09925 -0.99506 0.0523 1.00000  
    ## DepthMax          0.72863  0.68490 0.1969 1.00000  
    ## TempPrefMin       0.56071 -0.82801 0.2013 1.00000  
    ## TempPrefMax       0.51263  0.85861 0.0027 1.00000  
    ## DorsalSpinesMean  0.97411 -0.22609 0.5280 0.49451  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.0491 -0.0887
    ## BodyShapeI3f         0.1007  0.0309
    ## BodyShapeI4e        -0.4143 -0.0493
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -0.3899 -0.0284
    ## OperculumPresentyes  0.0866  0.0063
    ## FeedingPathb         0.0382  0.0027
    ## FeedingPathp        -0.3816 -0.0269
    ## RepGuild11b         -0.3955 -0.0179
    ## RepGuild12g         -0.4873 -0.0630
    ## RepGuild13n          0.0721  0.0076
    ## RepGuild22eb        -0.3955 -0.0179
    ## RepGuild23n         -0.4873 -0.0630
    ## RepGuild26s          0.0721  0.0076
    ## ParentalCare3p      -0.3670 -0.0147
    ## ParentalCare4n       0.1713  0.0068
    ## WaterPref1s          0.1177  0.0096
    ## WaterPref3a         -0.4002 -0.0327
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.4083 0.8991
    ## DemersPelag      0.0000 1.0000
    ## OperculumPresent 0.3415 1.0000
    ## FeedingPath      0.1471 1.0000
    ## RepGuild1        0.3380 0.1049
    ## RepGuild2        0.3380 0.1049
    ## ParentalCare     0.6332 1.0000
    ## WaterPref        0.4770 0.9590
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
pdtms_ef <- envfit(pdms_NMDS, FD_total_env[mixed_stratified_lakes,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL       0.33457  0.94237 0.2150 0.21479  
    ## Troph            -0.81047  0.58578 0.6005 0.05295 .
    ## DepthMin         -0.08091 -0.99672 0.1488 0.20080  
    ## DepthMax          0.58311  0.81240 0.2882 0.12488  
    ## TempPrefMin       0.31951 -0.94758 0.1604 0.44555  
    ## TempPrefMax       0.05947 -0.99823 0.0484 0.49750  
    ## DorsalSpinesMean  0.84071  0.54148 0.5398 0.04995 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.3170 -0.0110
    ## BodyShapeI3f         0.1533  0.0111
    ## BodyShapeI4e        -0.3424 -0.0277
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -0.3184  0.0008
    ## OperculumPresentyes  0.1061 -0.0003
    ## FeedingPathb         0.0442  0.0013
    ## FeedingPathp        -0.3096 -0.0094
    ## RepGuild11b         -0.3170 -0.0110
    ## RepGuild12g         -0.4098 -0.0756
    ## RepGuild13n          0.0874  0.0125
    ## RepGuild22eb        -0.3170 -0.0110
    ## RepGuild23n         -0.4098 -0.0756
    ## RepGuild26s          0.0874  0.0125
    ## ParentalCare3p      -0.2873 -0.0013
    ## ParentalCare4n       0.2234  0.0010
    ## WaterPref1s          0.1472  0.0130
    ## WaterPref3a         -0.3238 -0.0285
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.5321 0.03097 *
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.3456 0.13686  
    ## FeedingPath      0.1402 0.37962  
    ## RepGuild1        0.3511 0.01299 *
    ## RepGuild2        0.3511 0.01299 *
    ## ParentalCare     0.6565 0.29071  
    ## WaterPref        0.4913 0.04995 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
pdtms_efp <- p.adjust.envfit(pdtms_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdtms_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.33457  0.94237 0.2150 1.0000
    ## Troph            -0.81047  0.58578 0.6005 0.7942
    ## DepthMin         -0.08091 -0.99672 0.1488 1.0000
    ## DepthMax          0.58311  0.81240 0.2882 1.0000
    ## TempPrefMin       0.31951 -0.94758 0.1604 1.0000
    ## TempPrefMax       0.05947 -0.99823 0.0484 1.0000
    ## DorsalSpinesMean  0.84071  0.54148 0.5398 0.7493
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.3170 -0.0110
    ## BodyShapeI3f         0.1533  0.0111
    ## BodyShapeI4e        -0.3424 -0.0277
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -0.3184  0.0008
    ## OperculumPresentyes  0.1061 -0.0003
    ## FeedingPathb         0.0442  0.0013
    ## FeedingPathp        -0.3096 -0.0094
    ## RepGuild11b         -0.3170 -0.0110
    ## RepGuild12g         -0.4098 -0.0756
    ## RepGuild13n          0.0874  0.0125
    ## RepGuild22eb        -0.3170 -0.0110
    ## RepGuild23n         -0.4098 -0.0756
    ## RepGuild26s          0.0874  0.0125
    ## ParentalCare3p      -0.2873 -0.0013
    ## ParentalCare4n       0.2234  0.0010
    ## WaterPref1s          0.1472  0.0130
    ## WaterPref3a         -0.3238 -0.0285
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.5321 0.4645
    ## DemersPelag      0.0000 1.0000
    ## OperculumPresent 0.3456 1.0000
    ## FeedingPath      0.1402 1.0000
    ## RepGuild1        0.3511 0.1948
    ## RepGuild2        0.3511 0.1948
    ## ParentalCare     0.6565 1.0000
    ## WaterPref        0.4913 0.7493
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
pdtom_ef <- envfit(pdom_NMDS, FD_total_env[ocean_mixed_sites,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
pdtom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.64050 -0.76796 0.5485 0.008991 **
    ## Troph             0.95598 -0.29345 0.5829 0.039960 * 
    ## DepthMin         -0.89960  0.43672 0.1203 0.401598   
    ## DepthMax         -0.07564 -0.99713 0.0636 0.729271   
    ## TempPrefMin      -0.23517  0.97195 0.4275 0.040959 * 
    ## TempPrefMax       0.38399 -0.92334 0.0023 0.991009   
    ## DorsalSpinesMean -0.51746  0.85571 0.3087 0.113886   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.1684  0.1079
    ## BodyShapeI3f         0.0281 -0.0180
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentyes  0.0000  0.0000
    ## FeedingPathb         0.0000  0.0000
    ## RepGuild13n          0.0000  0.0000
    ## RepGuild26s          0.0000  0.0000
    ## ParentalCare4n       0.0000  0.0000
    ## WaterPref1s          0.0000  0.0000
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.1011 0.7343
    ## DemersPelag      0.0000 1.0000
    ## OperculumPresent 0.0000 1.0000
    ## FeedingPath      0.0000 1.0000
    ## RepGuild1        0.0000 1.0000
    ## RepGuild2        0.0000 1.0000
    ## ParentalCare     0.0000 1.0000
    ## WaterPref        0.0000 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
pdtom_efp <- p.adjust.envfit(pdtom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdtom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.64050 -0.76796 0.5485 0.1349
    ## Troph             0.95598 -0.29345 0.5829 0.5994
    ## DepthMin         -0.89960  0.43672 0.1203 1.0000
    ## DepthMax         -0.07564 -0.99713 0.0636 1.0000
    ## TempPrefMin      -0.23517  0.97195 0.4275 0.6144
    ## TempPrefMax       0.38399 -0.92334 0.0023 1.0000
    ## DorsalSpinesMean -0.51746  0.85571 0.3087 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.1684  0.1079
    ## BodyShapeI3f         0.0281 -0.0180
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentyes  0.0000  0.0000
    ## FeedingPathb         0.0000  0.0000
    ## RepGuild13n          0.0000  0.0000
    ## RepGuild26s          0.0000  0.0000
    ## ParentalCare4n       0.0000  0.0000
    ## WaterPref1s          0.0000  0.0000
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.1011      1
    ## DemersPelag      0.0000      1
    ## OperculumPresent 0.0000      1
    ## FeedingPath      0.0000      1
    ## RepGuild1        0.0000      1
    ## RepGuild2        0.0000      1
    ## ParentalCare     0.0000      1
    ## WaterPref        0.0000      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
pdtso_ef <- envfit(pdso_NMDS, FD_total_env[ocean_stratified_sites,c(6:20)], permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2  Pr(>r)  
    ## MaxLengthTL      -0.48095 -0.87675 0.0372 0.51049  
    ## Troph             0.97238  0.23341 0.8241 0.04795 *
    ## DepthMin          0.12195  0.99254 0.0510 0.26873  
    ## DepthMax         -0.46421  0.88572 0.1642 0.07692 .
    ## TempPrefMin      -0.93579 -0.35256 0.1858 0.66633  
    ## TempPrefMax       0.01813 -0.99984 0.1300 0.06394 .
    ## DorsalSpinesMean -0.98513  0.17179 0.5840 0.06194 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.1725  0.0062
    ## BodyShapeI3f        -0.0947 -0.0187
    ## BodyShapeI4e         0.2950  0.0280
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno   0.2684  0.0032
    ## OperculumPresentyes -0.1073 -0.0013
    ## FeedingPathb        -0.0423  0.0014
    ## FeedingPathp         0.2537 -0.0082
    ## RepGuild11b          0.2680 -0.0175
    ## RepGuild12g          0.3650  0.0291
    ## RepGuild13n         -0.0907 -0.0037
    ## RepGuild22eb         0.2680 -0.0175
    ## RepGuild23n          0.3650  0.0291
    ## RepGuild26s         -0.0907 -0.0037
    ## ParentalCare3p       0.2507  0.0152
    ## ParentalCare4n      -0.2507 -0.0152
    ## WaterPref1s         -0.1581 -0.0050
    ## WaterPref3a          0.2846  0.0089
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.3754 0.07093 .
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.2994 0.16084  
    ## FeedingPath      0.1116 0.34565  
    ## RepGuild1        0.3199 0.01199 *
    ## RepGuild2        0.3199 0.01199 *
    ## ParentalCare     0.6556 0.27373  
    ## WaterPref        0.4681 0.05794 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
pdtso_efp <- p.adjust.envfit(pdtso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
pdtso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL      -0.48095 -0.87675 0.0372 1.0000
    ## Troph             0.97238  0.23341 0.8241 0.7193
    ## DepthMin          0.12195  0.99254 0.0510 1.0000
    ## DepthMax         -0.46421  0.88572 0.1642 1.0000
    ## TempPrefMin      -0.93579 -0.35256 0.1858 1.0000
    ## TempPrefMax       0.01813 -0.99984 0.1300 0.9590
    ## DorsalSpinesMean -0.98513  0.17179 0.5840 0.9291
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.1725  0.0062
    ## BodyShapeI3f        -0.0947 -0.0187
    ## BodyShapeI4e         0.2950  0.0280
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno   0.2684  0.0032
    ## OperculumPresentyes -0.1073 -0.0013
    ## FeedingPathb        -0.0423  0.0014
    ## FeedingPathp         0.2537 -0.0082
    ## RepGuild11b          0.2680 -0.0175
    ## RepGuild12g          0.3650  0.0291
    ## RepGuild13n         -0.0907 -0.0037
    ## RepGuild22eb         0.2680 -0.0175
    ## RepGuild23n          0.3650  0.0291
    ## RepGuild26s         -0.0907 -0.0037
    ## ParentalCare3p       0.2507  0.0152
    ## ParentalCare4n      -0.2507 -0.0152
    ## WaterPref1s         -0.1581 -0.0050
    ## WaterPref3a          0.2846  0.0089
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.3754 1.0000
    ## DemersPelag      0.0000 1.0000
    ## OperculumPresent 0.2994 1.0000
    ## FeedingPath      0.1116 1.0000
    ## RepGuild1        0.3199 0.1798
    ## RepGuild2        0.3199 0.1798
    ## ParentalCare     0.6556 1.0000
    ## WaterPref        0.4681 0.8691
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

## Phylogenetic mantel tests

- We used mantel tests to determine significance between the
  phylogenetic and environmental distance matrices.

``` r
### Environmental
# All sites 
enve_dist_t <- dist(scaled_env[surveyed_sites_env,c(1)], method = "euclidean")
pdea_mant_t <- mantel(pde_dist, enve_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
pdea_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pde_dist, ydis = enve_dist_t, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1742 
    ##       Significance: 0.302 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.230 0.265 0.286 0.299 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_s <- dist(scaled_env[surveyed_sites_env,c(2)], method = "euclidean")
pdea_mant_s <- mantel(pde_dist, enve_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
pdea_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pde_dist, ydis = enve_dist_s, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5397 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.456 0.478 0.500 0.517 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_o <- dist(scaled_env[surveyed_sites_env,c(3)], method = "euclidean")
pdea_mant_o <- mantel(pde_dist, enve_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
pdea_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pde_dist, ydis = enve_dist_o, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4576 
    ##       Significance: 0.288 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.533 0.570 0.605 0.646 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_p <- dist(scaled_env[surveyed_sites_env,c(4)], method = "euclidean")
pdea_mant_p <- mantel(pde_dist, enve_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
pdea_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pde_dist, ydis = enve_dist_p, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7157 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.513 0.542 0.584 0.606 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdea_mant_pv <- rbind(pdea_mant_t$signif, pdea_mant_s$signif, pdea_mant_o$signif, pdea_mant_p$signif)
pdea_mant_pv <- pdea_mant_pv[,1]
pdea_mant_pv <- p.adjust(pdea_mant_pv, method = "bonferroni")
pdea_mant_pv
```

    ## [1] 1.000 0.008 1.000 0.008

``` r
# Mixed and stratified lakes
envems_dist_t <- dist(scaled_env[mixed_stratified_lakes,c(1)], method = "euclidean")
pdems_mant_t <- mantel(pdems_dist, envems_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdems_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdems_dist, ydis = envems_dist_t, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1563 
    ##       Significance: 0.302 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.227 0.252 0.281 0.321 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_s <- dist(scaled_env[mixed_stratified_lakes,c(2)], method = "euclidean")
pdems_mant_s <- mantel(pdems_dist, envems_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdems_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdems_dist, ydis = envems_dist_s, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4842 
    ##       Significance: 0.01 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.396 0.433 0.455 0.483 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_o <- dist(scaled_env[mixed_stratified_lakes,c(3)], method = "euclidean")
pdems_mant_o <- mantel(pdems_dist, envems_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdems_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdems_dist, ydis = envems_dist_o, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3426 
    ##       Significance: 0.289 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.419 0.450 0.480 0.549 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_p <- dist(scaled_env[mixed_stratified_lakes,c(4)], method = "euclidean")
pdems_mant_p <- mantel(pdems_dist, envems_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdems_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdems_dist, ydis = envems_dist_p, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6732 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.444 0.479 0.521 0.555 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdems_mant_pv <- rbind(pdems_mant_t$signif, pdems_mant_s$signif, pdems_mant_o$signif, pdems_mant_p$signif)
pdems_mant_pv <- pdems_mant_pv[,1]
pdems_mant_pv <- p.adjust(pdems_mant_pv, method = "bonferroni")
pdems_mant_pv
```

    ## [1] 1.000 0.040 1.000 0.004

``` r
# Ocean sites and mixed lakes
enveom_dist_t <- dist(scaled_env[ocean_mixed_sites_env,c(1)], method = "euclidean")
pdeom_mant_t <- mantel(pdeom_dist, enveom_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
pdeom_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeom_dist, ydis = enveom_dist_t, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1525 
    ##       Significance: 0.875 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.133 0.167 0.240 0.408 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_dist_s <- dist(scaled_env[ocean_mixed_sites_env,c(2)], method = "euclidean")
pdeom_mant_s <- mantel(pdeom_dist, enveom_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
pdeom_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeom_dist, ydis = enveom_dist_s, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5069 
    ##       Significance: 0.013 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.323 0.394 0.462 0.526 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_dist_o <- dist(scaled_env[ocean_mixed_sites_env,c(3)], method = "euclidean")
pdeom_mant_o <- mantel(pdeom_dist, enveom_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
pdeom_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeom_dist, ydis = enveom_dist_o, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5222 
    ##       Significance: 0.025 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.304 0.401 0.500 0.608 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_dist_p <- dist(scaled_env[ocean_mixed_sites_env,c(4)], method = "euclidean")
pdeom_mant_p <- mantel(pdeom_dist, enveom_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
pdeom_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeom_dist, ydis = enveom_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4543 
    ##       Significance: 0.02 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.269 0.336 0.403 0.549 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdeom_mant_pv <- rbind(pdeom_mant_t$signif, pdeom_mant_s$signif, pdeom_mant_o$signif, pdeom_mant_p$signif)
pdeom_mant_pv <- pdeom_mant_pv[,1]
pdeom_mant_pv <- p.adjust(pdeom_mant_pv, method = "bonferroni")
pdeom_mant_pv
```

    ## [1] 1.000 0.052 0.100 0.080

``` r
# Stratified lakes and ocean sites
enveso_dist_t <- dist(scaled_env[ocean_stratified_sites_env,c(1)], method = "euclidean")
pdeso_mant_t <- mantel(pdeso_dist, enveso_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
pdeso_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeso_dist, ydis = enveso_dist_t, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.05505 
    ##       Significance: 0.372 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0471 0.0719 0.1055 0.1350 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveso_dist_s <- dist(scaled_env[ocean_stratified_sites_env,c(2)], method = "euclidean")
pdeso_mant_s <- mantel(pdeso_dist, enveso_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
pdeso_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeso_dist, ydis = enveso_dist_s, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5084 
    ##       Significance: 0.014 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.412 0.456 0.488 0.527 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveso_dist_o <- dist(scaled_env[ocean_stratified_sites_env,c(3)], method = "euclidean")
pdeso_mant_o <- mantel(pdeso_dist, enveso_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
pdeso_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeso_dist, ydis = enveso_dist_o, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6287 
    ##       Significance: 0.456 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.698 0.718 0.740 0.764 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveso_dist_p <- dist(scaled_env[ocean_stratified_sites_env,c(4)], method = "euclidean")
pdeso_mant_p <- mantel(pdeso_dist, enveso_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
pdeso_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeso_dist, ydis = enveso_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7006 
    ##       Significance: 0.044 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.651 0.692 0.717 0.739 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdeso_mant_pv <- rbind(pdeso_mant_t$signif, pdeso_mant_s$signif, pdeso_mant_o$signif, pdeso_mant_p$signif)
pdeso_mant_pv <- pdeso_mant_pv[,1]
pdeso_mant_pv <- p.adjust(pdeso_mant_pv, method = "bonferroni")
pdeso_mant_pv
```

    ## [1] 1.000 0.056 1.000 0.176

``` r
# Stratified lakes
envs_dist <- dist(scaled_env[stratified_lakes,c(1:15)], method = "euclidean")
pds_mant <- mantel(pds_dist, envs_dist, method = "spearman", permutations = 999, na.rm = TRUE)
pds_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pds_dist, ydis = envs_dist, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1308 
    ##       Significance: 0.313 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.327 0.404 0.463 0.550 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# All sites 
envb_dist_vc <- dist(scaled_env[surveyed_sites_geo,c(5)], method = "euclidean")
pdba_mant_vc <- mantel(pdb_dist, envb_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_vc, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.07599 
    ##       Significance: 0.361 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.125 0.144 0.161 0.178 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_v <- dist(scaled_env[surveyed_sites_geo,c(6)], method = "euclidean")
pdba_mant_v <- mantel(pdb_dist, envb_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_v, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1511 
    ##       Significance: 0.32 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.203 0.225 0.240 0.255 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_sa <- dist(scaled_env[surveyed_sites_geo,c(7)], method = "euclidean")
pdba_mant_sa <- mantel(pdb_dist, envb_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_sa, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1249 
    ##       Significance: 0.439 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.190 0.208 0.223 0.245 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_dmin <- dist(scaled_env[surveyed_sites_geo,c(8)], method = "euclidean")
pdba_mant_dmin <- mantel(pdb_dist, envb_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_dmin, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3981 
    ##       Significance: 0.101 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.398 0.429 0.449 0.473 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_dmean <- dist(scaled_env[surveyed_sites_geo,c(9)], method = "euclidean")
pdba_mant_dmean <- mantel(pdb_dist, envb_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_dmean, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3675 
    ##       Significance: 0.55 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.476 0.505 0.531 0.552 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_dmed <- dist(scaled_env[surveyed_sites_geo,c(10)], method = "euclidean")
pdba_mant_dmed <- mantel(pdb_dist, envb_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_dmed, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3441 
    ##       Significance: 0.592 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.463 0.487 0.509 0.545 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_tl <- dist(scaled_env[surveyed_sites_geo,c(11)], method = "euclidean")
pdba_mant_tl <- mantel(pdb_dist, envb_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_tl, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4216 
    ##       Significance: 0.664 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.557 0.583 0.611 0.631 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_te <- dist(scaled_env[surveyed_sites_geo,c(12)], method = "euclidean")
pdba_mant_te <- mantel(pdb_dist, envb_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_te, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4329 
    ##       Significance: 0.525 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.546 0.576 0.594 0.615 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_p <- dist(scaled_env[surveyed_sites_geo,c(13)], method = "euclidean")
pdba_mant_p <- mantel(pdb_dist, envb_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_p, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.09127 
    ##       Significance: 0.495 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.161 0.179 0.193 0.209 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_md <- dist(scaled_env[surveyed_sites_geo,c(14)], method = "euclidean")
pdba_mant_md <- mantel(pdb_dist, envb_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_md, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.005785 
    ##       Significance: 0.89 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.172 0.209 0.229 0.246 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_dist_la <- dist(scaled_env[surveyed_sites_geo,c(15)], method = "euclidean")
pdba_mant_la <- mantel(pdb_dist, envb_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
pdba_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist_la, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.04335 
    ##       Significance: 0.684 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.164 0.184 0.200 0.213 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdba__mant_pv <- rbind(pdba_mant_vc$signif, pdba_mant_v$signif, pdba_mant_sa$signif, pdba_mant_dmin$signif, pdba_mant_dmean$signif, pdba_mant_dmed$signif, pdba_mant_tl$signif, pdba_mant_te$signif, pdba_mant_p$signif, pdba_mant_md$signif, pdba_mant_la$signif)
pdba__mant_pv <- pdba__mant_pv[,1]
pdba__mant_pv <- p.adjust(pdba__mant_pv, method = "bonferroni")
pdba__mant_pv
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1

``` r
# Mixed and stratified lakes
envbms_dist_vc <- dist(scaled_env[mixed_stratified_lakes_geo,c(5)], method = "euclidean")
pdbms_mant_vc <- mantel(pdbms_dist, envbms_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_vc, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1079 
    ##       Significance: 0.665 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0288 0.0565 0.0756 0.0949 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_v <- dist(scaled_env[mixed_stratified_lakes_geo,c(6)], method = "euclidean")
pdbms_mant_v <- mantel(pdbms_dist, envbms_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_v, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.03613 
    ##       Significance: 0.622 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.178 0.228 0.268 0.314 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_sa <- dist(scaled_env[mixed_stratified_lakes_geo,c(7)], method = "euclidean")
pdbms_mant_sa <- mantel(pdbms_dist, envbms_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_sa, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.01956 
    ##       Significance: 0.771 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.153 0.195 0.231 0.265 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_dmin <- dist(scaled_env[mixed_stratified_lakes_geo,c(8)], method = "euclidean")
pdbms_mant_dmin <- mantel(pdbms_dist, envbms_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_dmin, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1997 
    ##       Significance: 0.106 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.202 0.243 0.274 0.299 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_dmean <- dist(scaled_env[mixed_stratified_lakes_geo,c(9)], method = "euclidean")
pdbms_mant_dmean <- mantel(pdbms_dist, envbms_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_dmean, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1427 
    ##       Significance: 0.608 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.303 0.347 0.388 0.419 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_dmed <- dist(scaled_env[mixed_stratified_lakes_geo,c(10)], method = "euclidean")
pdbms_mant_dmed <- mantel(pdbms_dist, envbms_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_dmed, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1332 
    ##       Significance: 0.531 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.279 0.319 0.373 0.390 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_tl <- dist(scaled_env[mixed_stratified_lakes_geo,c(11)], method = "euclidean")
pdbms_mant_tl <- mantel(pdbms_dist, envbms_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_tl, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2854 
    ##       Significance: 0.544 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.485 0.533 0.586 0.632 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_te <- dist(scaled_env[mixed_stratified_lakes_geo,c(12)], method = "euclidean")
pdbms_mant_te <- mantel(pdbms_dist, envbms_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_te, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2735 
    ##       Significance: 0.632 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.490 0.544 0.589 0.615 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_p <- dist(scaled_env[mixed_stratified_lakes_geo,c(13)], method = "euclidean")
pdbms_mant_p <- mantel(pdbms_dist, envbms_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_p, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.07913 
    ##       Significance: 0.858 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0979 0.1321 0.1519 0.1700 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_md <- dist(scaled_env[mixed_stratified_lakes_geo,c(14)], method = "euclidean")
pdbms_mant_md <- mantel(pdbms_dist, envbms_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_md, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.0662 
    ##       Significance: 0.95 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.174 0.216 0.256 0.287 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_dist_la <- dist(scaled_env[mixed_stratified_lakes_geo,c(15)], method = "euclidean")
pdbms_mant_la <- mantel(pdbms_dist, envbms_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
pdbms_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist_la, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.07261 
    ##       Significance: 0.922 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.176 0.217 0.231 0.257 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdbms__mant_pv <- rbind(pdbms_mant_vc$signif, pdbms_mant_v$signif, pdbms_mant_sa$signif, pdbms_mant_dmin$signif, pdbms_mant_dmean$signif, pdbms_mant_dmed$signif, pdbms_mant_tl$signif, pdbms_mant_te$signif, pdbms_mant_p$signif, pdbms_mant_md$signif, pdbms_mant_la$signif)
pdbms__mant_pv <- pdbms__mant_pv[,1]
pdbms__mant_pv <- p.adjust(pdbms__mant_pv, method = "bonferroni")
pdbms__mant_pv
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1

``` r
# Ocean sites and mixed lakes
envbom_dist_vc <- dist(scaled_env[ocean_mixed_sites_geo,c(5)], method = "euclidean")
pdbom_mant_vc <- mantel(pdbom_dist, envbom_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_vc, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3166 
    ##       Significance: 0.03 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.251 0.286 0.328 0.362 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_v <- dist(scaled_env[ocean_mixed_sites_geo,c(6)], method = "euclidean")
pdbom_mant_v <- mantel(pdbom_dist, envbom_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_v, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2984 
    ##       Significance: 0.035 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.246 0.285 0.316 0.360 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_sa <- dist(scaled_env[ocean_mixed_sites_geo,c(7)], method = "euclidean")
pdbom_mant_sa <- mantel(pdbom_dist, envbom_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_sa, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.213 
    ##       Significance: 0.141 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.231 0.258 0.283 0.310 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_dmin <- dist(scaled_env[ocean_mixed_sites_geo,c(8)], method = "euclidean")
pdbom_mant_dmin <- mantel(pdbom_dist, envbom_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_dmin, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1833 
    ##       Significance: 0.035 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.131 0.166 0.202 0.234 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_dmean <- dist(scaled_env[ocean_mixed_sites_geo,c(9)], method = "euclidean")
pdbom_mant_dmean <- mantel(pdbom_dist, envbom_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_dmean, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.05706 
    ##       Significance: 0.382 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.185 0.228 0.258 0.308 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_dmed <- dist(scaled_env[ocean_mixed_sites_geo,c(10)], method = "euclidean")
pdbom_mant_dmed <- mantel(pdbom_dist, envbom_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_dmed, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.038 
    ##       Significance: 0.43 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.203 0.241 0.291 0.372 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_tl <- dist(scaled_env[ocean_mixed_sites_geo,c(11)], method = "euclidean")
pdbom_mant_tl <- mantel(pdbom_dist, envbom_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_tl, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.06144 
    ##       Significance: 0.118 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0715 0.1342 0.1722 0.2014 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_te <- dist(scaled_env[ocean_mixed_sites_geo,c(12)], method = "euclidean")
pdbom_mant_te <- mantel(pdbom_dist, envbom_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_te, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.06759 
    ##       Significance: 0.102 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0681 0.1250 0.1472 0.1915 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_p <- dist(scaled_env[ocean_mixed_sites_geo,c(13)], method = "euclidean")
pdbom_mant_p <- mantel(pdbom_dist, envbom_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2298 
    ##       Significance: 0.116 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.235 0.264 0.282 0.303 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_md <- dist(scaled_env[ocean_mixed_sites_geo,c(14)], method = "euclidean")
pdbom_mant_md <- mantel(pdbom_dist, envbom_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_md, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2516 
    ##       Significance: 0.033 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.147 0.208 0.264 0.343 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_dist_la <- dist(scaled_env[ocean_mixed_sites_geo,c(15)], method = "euclidean")
pdbom_mant_la <- mantel(pdbom_dist, envbom_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
pdbom_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist_la, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.165 
    ##       Significance: 0.165 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.186 0.215 0.248 0.278 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdbom__mant_pv <- rbind(pdbom_mant_vc$signif, pdbom_mant_v$signif, pdbom_mant_sa$signif, pdbom_mant_dmin$signif, pdbom_mant_dmean$signif, pdbom_mant_dmed$signif, pdbom_mant_tl$signif, pdbom_mant_te$signif, pdbom_mant_p$signif, pdbom_mant_md$signif, pdbom_mant_la$signif)
pdbom__mant_pv <- pdbom__mant_pv[,1]
pdbom__mant_pv <- p.adjust(pdbom__mant_pv, method = "bonferroni")
pdbom__mant_pv
```

    ##  [1] 0.330 0.385 1.000 0.385 1.000 1.000 1.000 1.000 1.000 0.363 1.000

``` r
# Stratified lakes and ocean sites
envbso_dist_vc <- dist(scaled_env[ocean_stratified_sites,c(5)], method = "euclidean")
pdbso_mant_vc <- mantel(pdbso_dist, envbso_dist_vc, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_vc
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_vc, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.289 
    ##       Significance: 0.147 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.313 0.345 0.362 0.390 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_v <- dist(scaled_env[ocean_stratified_sites,c(6)], method = "euclidean")
pdbso_mant_v <- mantel(pdbso_dist, envbso_dist_v, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_v
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_v, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3056 
    ##       Significance: 0.199 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.349 0.372 0.391 0.412 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_sa <- dist(scaled_env[ocean_stratified_sites,c(7)], method = "euclidean")
pdbso_mant_sa <- mantel(pdbso_dist, envbso_dist_sa, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_sa
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_sa, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.343 
    ##       Significance: 0.162 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.373 0.405 0.426 0.452 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_dmin <- dist(scaled_env[ocean_stratified_sites,c(8)], method = "euclidean")
pdbso_mant_dmin <- mantel(pdbso_dist, envbso_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_dmin, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4314 
    ##       Significance: 0.225 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.470 0.500 0.533 0.564 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_dmean <- dist(scaled_env[ocean_stratified_sites,c(9)], method = "euclidean")
pdbso_mant_dmean <- mantel(pdbso_dist, envbso_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_dmean
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_dmean, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5114 
    ##       Significance: 0.53 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.606 0.636 0.654 0.676 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_dmed <- dist(scaled_env[ocean_stratified_sites,c(10)], method = "euclidean")
pdbso_mant_dmed <- mantel(pdbso_dist, envbso_dist_dmed, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_dmed
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_dmed, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.491 
    ##       Significance: 0.57 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.594 0.618 0.650 0.682 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_tl <- dist(scaled_env[ocean_stratified_sites,c(11)], method = "euclidean")
pdbso_mant_tl <- mantel(pdbso_dist, envbso_dist_tl, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_tl
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_tl, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5883 
    ##       Significance: 0.918 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.698 0.716 0.726 0.735 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_te <- dist(scaled_env[ocean_stratified_sites,c(12)], method = "euclidean")
pdbso_mant_te <- mantel(pdbso_dist, envbso_dist_te, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_te
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_te, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6331 
    ##       Significance: 0.523 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.696 0.711 0.722 0.735 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_p <- dist(scaled_env[ocean_stratified_sites,c(13)], method = "euclidean")
pdbso_mant_p <- mantel(pdbso_dist, envbso_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4043 
    ##       Significance: 0.252 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.469 0.498 0.515 0.539 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_md <- dist(scaled_env[ocean_stratified_sites,c(14)], method = "euclidean")
pdbso_mant_md <- mantel(pdbso_dist, envbso_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_md
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_md, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02339 
    ##       Significance: 0.858 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.164 0.198 0.226 0.256 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_dist_la <- dist(scaled_env[ocean_stratified_sites,c(15)], method = "euclidean")
pdbso_mant_la <- mantel(pdbso_dist, envbso_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdbso_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist_la, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3431 
    ##       Significance: 0.173 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.381 0.417 0.430 0.446 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdbso__mant_pv <- rbind(pdbso_mant_vc$signif, pdbso_mant_v$signif, pdbso_mant_sa$signif, pdbso_mant_dmin$signif, pdbso_mant_dmean$signif, pdbso_mant_dmed$signif, pdbso_mant_tl$signif, pdbso_mant_te$signif, pdbso_mant_p$signif, pdbso_mant_md$signif, pdbso_mant_la$signif)
pdbso__mant_pv <- pdbso__mant_pv[,1]
pdbso__mant_pv <- p.adjust(pdbso__mant_pv, method = "bonferroni")
pdbso__mant_pv
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1

``` r
### Traits
# All sites
FDpd_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(1)]))
pdta_mant_s <- mantel(pd_dist, FDpd_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_s, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4625 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.327 0.361 0.377 0.394 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(2)]))
pdta_mant_l <- mantel(pd_dist, FDpd_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_l, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2403 
    ##       Significance: 0.02 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.183 0.208 0.234 0.252 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_t <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(3)]))
pdta_mant_t <- mantel(pd_dist, FDpd_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_t, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4925 
    ##       Significance: 0.011 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.405 0.430 0.454 0.485 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(4)]))
pdta_mant_dmin <- mantel(pd_dist, FDpd_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_dmin, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.597 
    ##       Significance: 0.014 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.529 0.557 0.579 0.606 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_dmax <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(5)]))
pdta_mant_dmax <- mantel(pd_dist, FDpd_dist_dmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_dmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_dmax, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2057 
    ##       Significance: 0.092 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.203 0.232 0.244 0.265 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_tmin <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(6)]))
pdta_mant_tmin <- mantel(pd_dist, FDpd_dist_tmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_tmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_tmin, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2389 
    ##       Significance: 0.236 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.283 0.305 0.328 0.348 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_tmax <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(7)]))
pdta_mant_tmax <- mantel(pd_dist, FDpd_dist_tmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_tmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_tmax, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1368 
    ##       Significance: 0.076 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.124 0.151 0.173 0.196 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_b <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(8:9)]))
pdta_mant_b <- mantel(pd_dist, FDpd_dist_b, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_b
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_b, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3554 
    ##       Significance: 0.046 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.305 0.349 0.370 0.393 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_o <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(10)]))
pdta_mant_o <- mantel(pd_dist, FDpd_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_o, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2912 
    ##       Significance: 0.152 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.304 0.347 0.374 0.409 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_f <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(11)]))
pdta_mant_f <- mantel(pd_dist, FDpd_dist_f, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_f
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_f, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.09556 
    ##       Significance: 0.375 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.176 0.206 0.268 0.268 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(12:13)]))
pdta_mant_ec <- mantel(pd_dist, FDpd_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_ec, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3142 
    ##       Significance: 0.01 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.213 0.236 0.293 0.300 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(14:15)]))
pdta_mant_es <- mantel(pd_dist, FDpd_dist_es, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_es
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_es, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3142 
    ##       Significance: 0.006 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.212 0.236 0.271 0.295 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_p <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(16)]))
pdta_mant_p <- mantel(pd_dist, FDpd_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_p, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6399 
    ##       Significance: 0.245 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.641 0.641 0.641 0.641 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_w <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites,c(17)]))
pdta_mant_w <- mantel(pd_dist, FDpd_dist_w, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19])
pdta_mant_w
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_w, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4477 
    ##       Significance: 0.059 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.397 0.448 0.483 0.484 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdta_mant_pv <- rbind(pdta_mant_s$signif, pdta_mant_l$signif, pdta_mant_t$signif, pdta_mant_dmin$signif, pdta_mant_dmax$signif, pdta_mant_tmin$signif, pdta_mant_tmax$signif, pdta_mant_b$signif, pdta_mant_o$signif, pdta_mant_f$signif, pdta_mant_ec$signif, pdta_mant_es$signif, pdta_mant_p$signif, pdta_mant_w$signif)
pdta_mant_pv <- pdta_mant_pv[,1]
pdta_mant_pv <- p.adjust(pdta_mant_pv, method = "bonferroni")
pdta_mant_pv
```

    ##  [1] 0.028 0.280 0.154 0.196 1.000 1.000 1.000 0.644 1.000 1.000 0.140 0.084
    ## [13] 1.000 0.826

``` r
# Mixed and stratified lakes
FDpdms_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(1)]))
pdtms_mant_s <- mantel(pdms_dist, FDpdms_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_s, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3744 
    ##       Significance: 0.01 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.248 0.296 0.324 0.366 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(2)]))
pdtms_mant_l <- mantel(pdms_dist, FDpdms_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_l, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.09701 
    ##       Significance: 0.095 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0931 0.1283 0.1574 0.1688 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_t <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(3)]))
pdtms_mant_t <- mantel(pdms_dist, FDpdms_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_t, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3895 
    ##       Significance: 0.038 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.327 0.373 0.409 0.462 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(4)]))
pdtms_mant_dmin <- mantel(pdms_dist, FDpdms_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_dmin, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5178 
    ##       Significance: 0.04 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.471 0.506 0.535 0.558 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_dmax <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(5)]))
pdtms_mant_dmax <- mantel(pdms_dist, FDpdms_dist_dmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_dmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_dmax, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1027 
    ##       Significance: 0.123 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.113 0.141 0.168 0.187 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_tmin <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(6)]))
pdtms_mant_tmin <- mantel(pdms_dist, FDpdms_dist_tmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_tmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_tmin, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.07532 
    ##       Significance: 0.32 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.149 0.184 0.212 0.236 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_tmax <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(7)]))
pdtms_mant_tmax <- mantel(pdms_dist, FDpdms_dist_tmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_tmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_tmax, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.007014 
    ##       Significance: 0.165 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0365 0.0711 0.0922 0.1041 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_b <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(8:9)]))
pdtms_mant_b <- mantel(pdms_dist, FDpdms_dist_b, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_b
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_b, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3577 
    ##       Significance: 0.02 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.258 0.320 0.341 0.358 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_o <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(10)]))
pdtms_mant_o <- mantel(pdms_dist, FDpdms_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_o, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1905 
    ##       Significance: 0.121 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.197 0.250 0.278 0.324 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_f <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(11)]))
pdtms_mant_f <- mantel(pdms_dist, FDpdms_dist_f, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_f
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_f, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02332 
    ##       Significance: 0.393 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0745 0.1035 0.1866 0.1866 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(12:13)]))
pdtms_mant_ec <- mantel(pdms_dist, FDpdms_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_ec, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2147 
    ##       Significance: 0.011 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.104 0.157 0.181 0.194 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(14:15)]))
pdtms_mant_es <- mantel(pdms_dist, FDpdms_dist_es, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_es
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_es, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2147 
    ##       Significance: 0.009 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0822 0.1156 0.1635 0.1880 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_p <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(16)]))
pdtms_mant_p <- mantel(pdms_dist, FDpdms_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_p, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6198 
    ##       Significance: 0.271 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ##  0.64  0.64  0.64  0.64 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdms_dist_w <- gowdis(as.data.frame(scaled_FD_total_env[mixed_stratified_lakes,c(17)]))
pdtms_mant_w <- mantel(pdms_dist, FDpdms_dist_w, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
pdtms_mant_w
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_w, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3725 
    ##       Significance: 0.054 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.289 0.372 0.423 0.425 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdtms_mant_pv <- rbind(pdtms_mant_s$signif, pdtms_mant_l$signif, pdtms_mant_t$signif, pdtms_mant_dmin$signif, pdtms_mant_dmax$signif, pdtms_mant_tmin$signif, pdtms_mant_tmax$signif, pdtms_mant_b$signif, pdtms_mant_o$signif, pdtms_mant_f$signif, pdtms_mant_ec$signif, pdtms_mant_es$signif, pdtms_mant_p$signif, pdtms_mant_w$signif)
pdtms_mant_pv <- pdtms_mant_pv[,1]
pdtms_mant_pv <- p.adjust(pdtms_mant_pv, method = "bonferroni")
pdtms_mant_pv
```

    ##  [1] 0.140 1.000 0.532 0.560 1.000 1.000 1.000 0.280 1.000 1.000 0.154 0.126
    ## [13] 1.000 0.756

``` r
# Ocean sites and mixed lakes
FDpdom_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(1)]))
pdtom_mant_s <- mantel(pdom_dist, FDpdom_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
pdtom_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdom_dist, ydis = FDpdom_dist_s, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3424 
    ##       Significance: 0.013 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.155 0.207 0.255 0.351 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdom_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(2)]))
pdtom_mant_l <- mantel(pdom_dist, FDpdom_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
pdtom_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdom_dist, ydis = FDpdom_dist_l, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2404 
    ##       Significance: 0.02 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.105 0.160 0.220 0.284 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdom_dist_t <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(3)]))
pdtom_mant_t <- mantel(pdom_dist, FDpdom_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
pdtom_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdom_dist, ydis = FDpdom_dist_t, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2605 
    ##       Significance: 0.023 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.186 0.228 0.255 0.302 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdom_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(4)]))
pdtom_mant_dmin <- mantel(pdom_dist, FDpdom_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
pdtom_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdom_dist, ydis = FDpdom_dist_dmin, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02977 
    ##       Significance: 0.476 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.112 0.167 0.209 0.305 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdom_dist_dmax <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(5)]))
pdtom_mant_dmax <- mantel(pdom_dist, FDpdom_dist_dmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
pdtom_mant_dmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdom_dist, ydis = FDpdom_dist_dmax, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1447 
    ##       Significance: 0.162 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.170 0.224 0.269 0.319 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdom_dist_tmin <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(6)]))
pdtom_mant_tmin <- mantel(pdom_dist, FDpdom_dist_tmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
pdtom_mant_tmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdom_dist, ydis = FDpdom_dist_tmin, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.177 
    ##       Significance: 0.054 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.113 0.184 0.250 0.297 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdom_dist_tmax <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(7)]))
pdtom_mant_tmax <- mantel(pdom_dist, FDpdom_dist_tmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
pdtom_mant_tmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdom_dist, ydis = FDpdom_dist_tmax, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2323 
    ##       Significance: 0.039 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.141 0.208 0.273 0.332 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# FDpdom_dist_b <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(8:9)]))
# pdtom_mant_b <- mantel(pdom_dist, FDpdom_dist_b, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# pdtom_mant_b
# FDpdom_dist_o <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(10)]))
# pdtom_mant_o <- mantel(pdom_dist, FDpdom_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# pdtom_mant_o
# FDpdom_dist_f <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(11)]))
# pdtom_mant_f <- mantel(pdom_dist, FDpdom_dist_f, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# pdtom_mant_f
# FDpdom_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(12:13)]))
# pdtom_mant_ec <- mantel(pdom_dist, FDpdom_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# pdtom_mant_ec
# FDpdom_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(14:15)]))
# pdtom_mant_es <- mantel(pdom_dist, FDpdom_dist_es, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# pdtom_mant_es
# FDpdom_dist_p <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(16)]))
# pdtom_mant_p <- mantel(pdom_dist, FDpdom_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# pdtom_mant_p
# FDpdom_dist_w <- gowdis(as.data.frame(scaled_FD_total_env[ocean_mixed_sites,c(17)]))
# pdtom_mant_w <- mantel(pdom_dist, FDpdom_dist_w, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19])
# pdtom_mant_w
# Error in cor(as.vector(xdis), ydis, method = method, use = use) : 
#   no complete element pairs
# Adjust p-values
pdtom_mant_pv <- rbind(pdtom_mant_s$signif, pdtom_mant_l$signif, pdtom_mant_t$signif, pdtom_mant_dmin$signif, pdtom_mant_dmax$signif, pdtom_mant_tmin$signif, pdtom_mant_tmax$signif)
pdtom_mant_pv <- pdtom_mant_pv[,1]
pdtom_mant_pv <- p.adjust(pdtom_mant_pv, method = "bonferroni")
pdtom_mant_pv
```

    ## [1] 0.091 0.140 0.161 1.000 1.000 0.378 0.273

``` r
# Stratified lakes and ocean sites
FDpdso_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(1)]))
pdtso_mant_s <- mantel(pdso_dist, FDpdso_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_s, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4399 
    ##       Significance: 0.007 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.336 0.365 0.388 0.433 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(2)]))
pdtso_mant_l <- mantel(pdso_dist, FDpdso_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_l, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1063 
    ##       Significance: 0.027 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0669 0.0877 0.1065 0.1256 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_t <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(3)]))
pdtso_mant_t <- mantel(pdso_dist, FDpdso_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_t, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7366 
    ##       Significance: 0.044 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.696 0.729 0.751 0.786 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(4)]))
pdtso_mant_dmin <- mantel(pdso_dist, FDpdso_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_dmin, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6587 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.515 0.542 0.573 0.589 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_dmax <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(5)]))
pdtso_mant_dmax <- mantel(pdso_dist, FDpdso_dist_dmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_dmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_dmax, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.002262 
    ##       Significance: 0.087 
    ## 
    ## Upper quantiles of permutations (null model):
    ##     90%     95%   97.5%     99% 
    ## -0.0082  0.0195  0.0480  0.0658 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_tmin <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(6)]))
pdtso_mant_tmin <- mantel(pdso_dist, FDpdso_dist_tmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_tmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_tmin, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.0793 
    ##       Significance: 0.425 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.169 0.198 0.225 0.249 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_tmax <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(7)]))
pdtso_mant_tmax <- mantel(pdso_dist, FDpdso_dist_tmax, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_tmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_tmax, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.05174 
    ##       Significance: 0.009 
    ## 
    ## Upper quantiles of permutations (null model):
    ##      90%      95%    97.5%      99% 
    ## -0.02954  0.00469  0.02380  0.04689 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_b <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(8:9)]))
pdtso_mant_b <- mantel(pdso_dist, FDpdso_dist_b, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_b
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_b, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1497 
    ##       Significance: 0.048 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.106 0.144 0.169 0.188 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_o <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(10)]))
pdtso_mant_o <- mantel(pdso_dist, FDpdso_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_o
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_o, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1138 
    ##       Significance: 0.155 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.131 0.197 0.208 0.246 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_f <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(11)]))
pdtso_mant_f <- mantel(pdso_dist, FDpdso_dist_f, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_f
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_f, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.09304 
    ##       Significance: 0.479 
    ## 
    ## Upper quantiles of permutations (null model):
    ##      90%      95%    97.5%      99% 
    ## 0.000949 0.037025 0.115822 0.115822 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(12:13)]))
pdtso_mant_ec <- mantel(pdso_dist, FDpdso_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_ec, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1341 
    ##       Significance: 0.007 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0179 0.0795 0.1052 0.1081 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(14:15)]))
pdtso_mant_es <- mantel(pdso_dist, FDpdso_dist_es, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_es
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_es, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1341 
    ##       Significance: 0.008 
    ## 
    ## Upper quantiles of permutations (null model):
    ##     90%     95%   97.5%     99% 
    ## 0.00809 0.05655 0.08172 0.10523 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_p <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(16)]))
pdtso_mant_p <- mantel(pdso_dist, FDpdso_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_p, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6428 
    ##       Significance: 0.259 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.654 0.654 0.654 0.654 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpdso_dist_w <- gowdis(as.data.frame(scaled_FD_total_env[ocean_stratified_sites,c(17)]))
pdtso_mant_w <- mantel(pdso_dist, FDpdso_dist_w, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
pdtso_mant_w
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_w, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3062 
    ##       Significance: 0.061 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.262 0.306 0.356 0.365 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdtso_mant_pv <- rbind(pdtso_mant_s$signif, pdtso_mant_l$signif, pdtso_mant_t$signif, pdtso_mant_dmin$signif, pdtso_mant_dmax$signif, pdtso_mant_tmin$signif, pdtso_mant_tmax$signif, pdtso_mant_b$signif, pdtso_mant_o$signif, pdtso_mant_f$signif, pdtso_mant_ec$signif, pdtso_mant_es$signif, pdtso_mant_p$signif, pdtso_mant_w$signif)
pdtso_mant_pv <- pdtso_mant_pv[,1]
pdtso_mant_pv <- p.adjust(pdtso_mant_pv, method = "bonferroni")
pdtso_mant_pv
```

    ##  [1] 0.098 0.378 0.616 0.014 1.000 1.000 0.126 0.672 1.000 1.000 0.098 0.112
    ## [13] 1.000 0.854

## Plot phylogenetic NMDS and envfit results

- Includes a plot of correlated environmental variables

``` r
pd_NMDS_data.scores <- as.data.frame(scores(pd_NMDS))
pd_NMDS_data.scores$Stratification <- env[surveyed_sites,19]
pd_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
pd_NMDS_data.scores$Stratification <- factor(pd_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

pde_ef_coord_cont <- as.data.frame(scores(pde_ef, "vectors")) * ordiArrowMul(pde_ef)
# pdb_ef_coord_cont <- as.data.frame(scores(pdb_ef, "vectors")) * ordiArrowMul(pdb_ef)
# pdb_row_names <- c("minD")
# # Assign the new row names to the data frame
# pdb_ef_coord_cont <- data.frame(row.names = pdb_row_names, pdb_ef_coord_cont)

pd_ef_plot <- ggplot(data = pd_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type='t', level=0.95) +
  geom_point(data = pd_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  geom_text_repel(data = pd_NMDS_data.scores, label = pd_NMDS_data.scores$Lakes, size = 5, point.padding = 5, max.overlaps = 30, nudge_x = -0.01) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = pde_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = pde_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
            color = env_cont, label = row.names(pde_ef_coord_cont), size = 5) + 
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
  #              data = pdb_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#698B22") +
  # geom_text(data = pdb_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
  #           color = "#698B22", label = row.names(pdb_ef_coord_cont), size = 7) +
  theme(text = element_text(size = 22), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16), panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  annotate("text", x = -.48, y = .26, size = 5, 
           label = paste("Stress: ", round(pd_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ")
pd_ef_plot <- pd_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
pd_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20phylogenetic%20NMDS%20and%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_ef_plot.png", pd_ef_plot, width = 6, height = 4, units = "in")

# For 90 percent confidence interval: x = -.44, y = .26

# Determine outliers
ordered_pd_NMDS1_data.scores <- pd_NMDS_data.scores[order(pd_NMDS_data.scores$NMDS1), ]
Q1 <- ordered_pd_NMDS1_data.scores[6,1]
Q3 <- ordered_pd_NMDS1_data.scores[18,1]
IQR1 <- IQR(pd_NMDS_data.scores$NMDS1)
Q1 - 1.5*IQR1
```

    ## [1] -1.152046

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.148215

``` r
ordered_pd_NMDS1_data.scores$NMDS1
```

    ##  [1] -0.49293231 -0.48162038 -0.39547839 -0.36765966 -0.31641579 -0.31479471
    ##  [7] -0.20020966 -0.20001606 -0.04866105 -0.01315252  0.01843398  0.02291194
    ## [13]  0.05163823  0.21516810  0.24037173  0.26271219  0.27512126  0.31096335
    ## [19]  0.32748380  0.32776736  0.38786130  0.39050729

``` r
ordered_pd_NMDS2_data.scores <- pd_NMDS_data.scores[order(pd_NMDS_data.scores$NMDS2), ]
Q1 <- ordered_pd_NMDS2_data.scores[6,2]
Q3 <- ordered_pd_NMDS2_data.scores[18,2]
IQR2 <- IQR(pd_NMDS_data.scores$NMDS2)
Q1 - 1.5*IQR2
```

    ## [1] -0.3249817

``` r
Q3 + 1.5*IQR2
```

    ## [1] 0.3395479

``` r
ordered_pd_NMDS2_data.scores$NMDS2
```

    ##  [1] -0.209327936 -0.165818641 -0.156924537 -0.141171643 -0.101569451
    ##  [6] -0.082308191 -0.035953126 -0.035079002 -0.027185827 -0.024496723
    ## [11] -0.017944927  0.008820134  0.015522106  0.032498430  0.041726936
    ## [16]  0.078134930  0.095372222  0.096874460  0.099026995  0.132873992
    ## [21]  0.170573824  0.226355974

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_NMDS_boxplot.tiff", width = 1200, height = 800, res = 150, type = "cairo")

pd_dissimilarity_matrix <- as.matrix(pd_dist)

# Set the diagonal elements to NA
diag(pd_dissimilarity_matrix) <- NA

# Get the labels of the sites
sites <- attr(pd_dissimilarity_matrix, "Labels")

# Calculate the mean dissimilarity for each group
pd_mean_dissimilarity <- aggregate(pd_dissimilarity_matrix, by = list(stratification_group), FUN = mean, na.rm = TRUE)

row.names(pd_mean_dissimilarity) <- pd_mean_dissimilarity$Group.1

pd_mean_dissimilarity <- pd_mean_dissimilarity[,-1] 

pd_mean_dissimilarity <- as.data.frame(t(pd_mean_dissimilarity))

pd_mean_dissimilarity$Stratification <- env[surveyed_sites,19]

boxplot(Ocean ~ Stratification, pd_mean_dissimilarity[c(7,9,11,16,18,19),])
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
pd_mean_dissimilarity$Stratification <- factor(pd_mean_dissimilarity$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

pd_ocn_dist_violin_plot <- ggplot(pd_mean_dissimilarity, mapping = aes(x= Stratification, y= Ocean, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.8,
            width = 0.1) +
  geom_text_repel(data = pd_mean_dissimilarity, label = row.names(pd_mean_dissimilarity), size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  ylab("Ocean site distances") +
  xlab("Site type") +
  labs(colour = "Site type:  ", fill = "Site type:  ")
pd_ocn_dist_violin_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20phylogenetic%20NMDS%20and%20envfit%20results-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_ocn_dist_violin_plot.png", pd_ocn_dist_violin_plot, width = 8, height = 4, units = "in")

pd_mix_dist_violin_plot <- ggplot(pd_mean_dissimilarity, mapping = aes(x= Stratification, y= Mixed, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.8,
            width = 0.1) +
  geom_text_repel(data = pd_mean_dissimilarity, label = row.names(pd_mean_dissimilarity), size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  ylab("Mixed lake distances") +
  xlab("Site type") +
  labs(colour = "Site type:  ", fill = "Site type:  ")
pd_mix_dist_violin_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20phylogenetic%20NMDS%20and%20envfit%20results-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_mix_dist_violin_plot.png", pd_mix_dist_violin_plot, width = 8, height = 4, units = "in")

pd_strat_dist_violin_plot <- ggplot(pd_mean_dissimilarity, mapping = aes(x= Stratification, y= Stratified, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.8,
            width = 0.1) +
  geom_text_repel(data = pd_mean_dissimilarity, label = row.names(pd_mean_dissimilarity), size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  ylab("Stratified lake distances") +
  xlab("Site type") +
  labs(colour = "Site type:  ", fill = "Site type:  ")
pd_strat_dist_violin_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20phylogenetic%20NMDS%20and%20envfit%20results-4.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_strat_dist_violin_plot.png", pd_strat_dist_violin_plot, width = 8, height = 4, units = "in")
```

## Plot phylogenetic NMDS and trait envfit results

- Includes a plot of correlated trait variables

``` r
pdt_ef_coord_cont <- as.data.frame(scores(pdt_ef, "vectors")) * ordiArrowMul(pdt_ef)
pdt_ef_coord_cat = as.data.frame(scores(pdt_ef, "factors")) * ordiArrowMul(pdt_ef)

pdt_ef_plot <- ggplot(data = pd_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification),type="t",level = 0.95) +
  geom_point(data = pd_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  geom_text_repel(data = pd_NMDS_data.scores, label = pd_NMDS_data.scores$Lakes, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = pdt_ef_coord_cont, linewidth =1, alpha = 0.2, color = trait_cont) +
  geom_text(data = pdt_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
            color = trait_cont, label = row.names(pdt_ef_coord_cont), size = 5) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = pdt_ef_coord_cat, linewidth =1, alpha = 0.2, color = trait_cat) +
  geom_text(data = pdt_ef_coord_cat, aes(x = NMDS1, y = NMDS2), 
            color = trait_cat, label = row.names(pdt_ef_coord_cat), size = 5) + 
  theme(text = element_text(size = 22), legend.position = "bottom", legend.text = element_text(size = 16), panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  annotate("text", x = 0.50, y = 0.72, size = 5, 
           label = paste("Stress: ", round(pd_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "B")
pdt_ef_plot <- pdt_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
pdt_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20phylogenetic%20NMDS%20and%20trait%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pdt_ef_plot.png", pdt_ef_plot, width = 6, height = 4, units = "in")
```

## Boxplot of intradissimilarity heterogeneity

- The dissimilarity between sites that share the same stratification.

``` r
pd_bd_dist <- pd_bd$distances
pd_bd_dist <- as.data.frame(pd_bd_dist)
pd_bd_dist$X <- row.names(pd_bd_dist)
pd_bd_dist_env <- merge(pd_bd_dist, env[surveyed_sites,], by = "X", sort = F)

pd_bd_dist_env$Stratification <- factor(pd_bd_dist_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

pd_bd_dist_env_plot <- ggplot(pd_bd_dist_env, aes(x = Stratification, y = pd_bd_dist, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.8,
            width = 0.1) +
    geom_text_repel(data = pd_bd_dist_env, label = pd_bd_dist_env$X, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 16),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black")) +
  xlab("Site type") +
  ylab("Distance to Centroid") +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "C")
pd_bd_dist_env_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Boxplot%20of%20phylogenetic%20intradissimilarity%20heterogeneity-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_bd_dist_env_plot.png", pd_bd_dist_env_plot, width = 8.25, height = 4.13, units = "in")
```

PDis Dendrogram

``` r
# cluster communities using average-linkage algorithm
pd_dist_clust <- hclust(pd_dist, method = "average")

# Open a PNG device
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_dendro.tiff", width = 1200, height = 800, res = 150, type = "cairo")

# Create your plot using the plot() function
plot(pd_dist_clust, 
     xlab = "Sites",
     ylab = "Phylogenetic dissimilarity",
     main = "",
     sub = "",
     pch = 20,
     col= "black")

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_si.tiff", width = 1200, height = 800, res = 150, type = "cairo")

Si <- numeric(nrow(presabs_lake[surveyed_sites,]))
for (k in 2:(nrow(presabs_lake[surveyed_sites,])-1))
{
  sil<-silhouette(cutree(pd_dist_clust, k=k), pd_dist)
  Si[k] <- summary(sil)$avg.width
}
k.best<-which.max(Si)
plot(1:nrow(presabs_lake[surveyed_sites,]), Si, type = "h", main = "Silhouette", xlab = "K", ylab = "Width")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col = "red", font = 2, col.axis = "red")

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_pam.tiff", width = 1200, height = 800, res = 150, type = "cairo")

pd_pam <- pam(pd_dist,3)
plot(pd_pam)

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Environment/biogeography dissimilarity distances

``` r
### Environmental
mod <-  lm(FRic_total ~ temperature_median + salinity_median + oxygen_median + pH_median, FD_total_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = FRic_total ~ temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = FD_total_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.39166 -0.09479 -0.05193  0.15803  0.22816 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -6.56964    2.36826  -2.774 0.014923 *  
    ## temperature_median  0.01128    0.04544   0.248 0.807572    
    ## salinity_median     0.10737    0.02314   4.639 0.000383 ***
    ## oxygen_median       0.06039    0.06267   0.964 0.351589    
    ## pH_median           0.47868    0.40966   1.168 0.262132    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.191 on 14 degrees of freedom
    ##   (4 observations deleted due to missingness)
    ## Multiple R-squared:  0.9032, Adjusted R-squared:  0.8756 
    ## F-statistic: 32.67 on 4 and 14 DF,  p-value: 5.818e-07

``` r
vif(mod)
```

    ## temperature_median    salinity_median      oxygen_median          pH_median 
    ##           1.385456           3.435415           2.209454           4.821892

``` r
#define the variables we want to include in the correlation matrix
data <- env[surveyed_sites_env, c("temperature_median", "salinity_median", "oxygen_median", "pH_median")]
#create correlation matrix
cor(data)
```

    ##                    temperature_median salinity_median oxygen_median  pH_median
    ## temperature_median          1.0000000      -0.3571229    -0.1964940 -0.0995669
    ## salinity_median            -0.3571229       1.0000000     0.4436123  0.7726983
    ## oxygen_median              -0.1964940       0.4436123     1.0000000  0.6930775
    ## pH_median                  -0.0995669       0.7726983     0.6930775  1.0000000

``` r
# All sites
enve_dist <- dist(scaled_env[surveyed_sites_env,c(1:4)], method = "euclidean")
# Mixed and stratified lakes
envems_dist <- dist(scaled_env[mixed_stratified_lakes,c(1:4)], method = "euclidean")
# Ocean sites and mixed lakes
enveom_dist <- dist(scaled_env[ocean_mixed_sites_env,c(1:4)], method = "euclidean")
# Stratified lakes and ocean sites
enveso_dist <- dist(scaled_env[ocean_stratified_sites_env,c(1:4)], method = "euclidean")

### Geographic
mod <-  lm(FRic_total ~ volume_m3 + distance_to_ocean_mean_m + tidal_lag_time_minutes + max_depth + logArea, FD_total_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = FRic_total ~ volume_m3 + distance_to_ocean_mean_m + 
    ##     tidal_lag_time_minutes + max_depth + logArea, data = FD_total_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.73288 -0.18045  0.00905  0.35205  0.51369 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)               1.620e+00  6.946e-01   2.333    0.034 *
    ## volume_m3                -6.106e-09  1.741e-08  -0.351    0.731  
    ## distance_to_ocean_mean_m -2.856e-03  1.845e-03  -1.548    0.142  
    ## tidal_lag_time_minutes   -6.682e-04  2.646e-03  -0.253    0.804  
    ## max_depth                 5.355e-03  1.271e-02   0.421    0.680  
    ## logArea                  -1.990e-02  7.132e-02  -0.279    0.784  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4274 on 15 degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## Multiple R-squared:  0.4775, Adjusted R-squared:  0.3033 
    ## F-statistic: 2.741 on 5 and 15 DF,  p-value: 0.05938

``` r
vif(mod)
```

    ##                volume_m3 distance_to_ocean_mean_m   tidal_lag_time_minutes 
    ##                 2.430634                 5.894264                 4.084926 
    ##                max_depth                  logArea 
    ##                 2.364415                 1.982650

``` r
#define the variables we want to include in the correlation matrix
data <- SR_env[surveyed_sites_geo, c("volume_m3", "distance_to_ocean_mean_m", "tidal_lag_time_minutes", "max_depth", "logArea")]
#create correlation matrix
cor(data)
```

    ##                           volume_m3 distance_to_ocean_mean_m
    ## volume_m3                 1.0000000               -0.3294437
    ## distance_to_ocean_mean_m -0.3294437                1.0000000
    ## tidal_lag_time_minutes   -0.3248531                0.8619003
    ## max_depth                 0.3040825                0.5300135
    ## logArea                   0.6722191               -0.1046599
    ##                          tidal_lag_time_minutes max_depth    logArea
    ## volume_m3                            -0.3248531 0.3040825  0.6722191
    ## distance_to_ocean_mean_m              0.8619003 0.5300135 -0.1046599
    ## tidal_lag_time_minutes                1.0000000 0.3665876 -0.1497090
    ## max_depth                             0.3665876 1.0000000  0.4001014
    ## logArea                              -0.1497090 0.4001014  1.0000000

``` r
# All sites
envb_dist <- dist(scaled_env[surveyed_sites_geo,c(6,8,11,14:15)], method = "euclidean")
# Mixed and stratified lakes
envbms_dist <- dist(scaled_env[mixed_stratified_lakes_geo,c(6,8,11,14:15)], method = "euclidean")
# Ocean sites and mixed lakes
envbom_dist <- dist(scaled_env[ocean_mixed_sites_geo,c(6,8,11,14:15)], method = "euclidean")
# Stratified lakes and ocean sites
envbso_dist <- dist(scaled_env[ocean_stratified_sites,c(6,8,11,14:15)], method = "euclidean")
```

## Environment/biogeography NMDS

- To constrain dissimilarities we perform Nonmetric Multidimensional
  Scaling (NMDS), which tries to find a stable solution using the
  metaMDS package.

``` r
### Environmental
# All sites 
enve_NMDS <- metaMDS(enve_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
envems_NMDS <- metaMDS(envems_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
enveom_NMDS <- metaMDS(enveom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
enveso_NMDS <- metaMDS(enveso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)

### Geographic
# All sites 
envb_NMDS <- metaMDS(envb_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
envbms_NMDS <- metaMDS(envbms_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
envbom_NMDS <- metaMDS(envbom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
envbso_NMDS <- metaMDS(envbso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
```

## Determine stratification homogeneity using environmental/Geographical distances

- We use betadisper to determine homogeneity of the
  environmental/Geographical distances based on the stratification
  category. Is the dispersion of environmental/geographical distances
  similar within stratification categories?

``` r
### Environmental
groupse <- env[surveyed_sites_env,19]
enve_bd <- betadisper(enve_dist, groupse)
anova(enve_bd)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## Groups     2 4.0435 2.02174  4.7926 0.02339 *
    ## Residuals 16 6.7495 0.42184                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
permutest(enve_bd, pairwise = TRUE)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
    ## Groups     2 4.0435 2.02174 4.7926    999  0.021 *
    ## Residuals 16 6.7495 0.42184                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##               Mixed    Ocean Stratified
    ## Mixed               0.295000      0.009
    ## Ocean      0.301869               0.052
    ## Stratified 0.024548 0.066673

``` r
(enve_bd.HSD <- TukeyHSD(enve_bd))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                        diff          lwr       upr     p adj
    ## Ocean-Mixed      -0.3198234 -1.454419946 0.8147731 0.7511538
    ## Stratified-Mixed  0.8209829 -0.016972103 1.6589380 0.0552689
    ## Stratified-Ocean  1.1408064  0.006209893 2.2754029 0.0486582

``` r
boxplot(enve_bd)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Determine%20stratification%20homogeneity%20using%20environmental/Geographical%20distances-1.png)<!-- -->

``` r
### Geographic
groupsb <- env[surveyed_sites_geo,19]
envb_bd <- betadisper(envb_dist, groupsb)
anova(envb_bd)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq Mean Sq F value Pr(>F)
    ## Groups     2  2.6692 1.33459  1.6816  0.214
    ## Residuals 18 14.2853 0.79363

``` r
permutest(envb_bd, pairwise = TRUE)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
    ## Groups     2  2.6692 1.33459 1.6816    999  0.205
    ## Residuals 18 14.2853 0.79363                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              Mixed   Ocean Stratified
    ## Mixed              0.17600      0.524
    ## Ocean      0.18357              0.237
    ## Stratified 0.43965 0.24059

``` r
(envb_bd.HSD <- TukeyHSD(envb_bd))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                        diff        lwr       upr     p adj
    ## Ocean-Mixed       0.8692576 -0.3956640 2.1341791 0.2132322
    ## Stratified-Mixed  0.1866380 -0.9900686 1.3633445 0.9140611
    ## Stratified-Ocean -0.6826196 -1.9105110 0.5452718 0.3525629

``` r
boxplot(envb_bd)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Determine%20stratification%20homogeneity%20using%20environmental/Geographical%20distances-2.png)<!-- -->

## Determine stratification dissimilarity using environmental/geographical distances

- We use adonis to determine if dissimilarities of species
  environmental/geographical distances by stratification categories are
  significant and how much of the variation is explained by
  environmental/Geographical dissimilarities.

``` r
### Environmental
enve_pmc <- adonis2(enve_dist ~ env[surveyed_sites_env,19], permutations = 999)
enve_pmc
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = enve_dist ~ env[surveyed_sites_env, 19], permutations = 999)
    ##                             Df SumOfSqs      R2      F Pr(>F)    
    ## env[surveyed_sites_env, 19]  2   37.503 0.52088 8.6972  0.001 ***
    ## Residual                    16   34.497 0.47912                  
    ## Total                       18   72.000 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratifcation
enve_pmc_pair <- pairwise.adonis(enve_dist, env[surveyed_sites_env,19], p.adjust.m = "bonferroni", perm = 999)
enve_pmc_pair
```

    ##                 pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 27.117826 11.436250 0.4496044   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 23.559911  7.562930 0.4566179   0.009      0.027   .
    ## 3      Mixed vs Ocean  1  1.774351  2.057938 0.1861051   0.144      0.432

``` r
### Geographic
envb_pmc <- adonis2(envb_dist ~ env[surveyed_sites_geo,19], permutations = 999)
envb_pmc
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = envb_dist ~ env[surveyed_sites_geo, 19], permutations = 999)
    ##                             Df SumOfSqs      R2      F Pr(>F)    
    ## env[surveyed_sites_geo, 19]  2   39.138 0.46632 7.8641  0.001 ***
    ## Residual                    18   44.790 0.53368                  
    ## Total                       20   83.928 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratifcation
envb_pmc_pair <- pairwise.adonis(envb_dist, env[surveyed_sites_geo,19], p.adjust.m = "bonferroni", perm = 999)
envb_pmc_pair
```

    ##                 pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1  11.38554  8.044438 0.3822596   0.003      0.009   *
    ## 2 Stratified vs Ocean  1  33.99081 10.935220 0.4767872   0.001      0.003   *
    ## 3      Mixed vs Ocean  1  13.47960  4.376367 0.2846164   0.008      0.024   .

## Envfit trait influence on environmental/Geographical distances

- We use envfit to determine significantly correlated functional
  variables to our environmental/Geographical NMDS. We will use these
  results, if significant, in our figure.

``` r
### Environmental
# All sites for figure
enve_ef <- envfit(enve_NMDS, FD_total_env[surveyed_sites_env,c(56)], permutations = 1000, na.rm = TRUE, strata = env[surveyed_sites_env,19])
enve_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## [1,] -0.99338 -0.11487 0.6303 0.07093 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
enve_efp <- p.adjust.envfit(enve_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
enve_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## [1,] -0.99338 -0.11487 0.6303 0.07093 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
envea_ef <- envfit(enve_NMDS, FD_total_env[surveyed_sites_env,c(6:20)], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
envea_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)  
    ## MaxLengthTL       0.86550 -0.50090 0.1250  0.510  
    ## Troph            -0.99338 -0.11487 0.6303  0.068 .
    ## DepthMin          0.34970 -0.93686 0.0226  0.713  
    ## DepthMax          0.98912 -0.14712 0.2076  0.431  
    ## TempPrefMin       0.67062 -0.74180 0.1311  0.850  
    ## TempPrefMax      -0.98765  0.15670 0.0121  0.885  
    ## DorsalSpinesMean  0.98478 -0.17380 0.5093  0.289  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.2453 -0.9704
    ## BodyShapeI3f         0.5851  0.0547
    ## BodyShapeI4e        -1.7789  0.3075
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.8894 -0.0241
    ## OperculumPresentyes  0.5039  0.0064
    ## FeedingPathb         0.2181  0.1023
    ## FeedingPathp        -1.8541 -0.8695
    ## RepGuild11b         -1.9179 -1.8022
    ## RepGuild12g         -1.7405  1.1871
    ## RepGuild13n          0.3374 -0.0357
    ## RepGuild22eb        -1.9179 -1.8022
    ## RepGuild23n         -1.7405  1.1871
    ## RepGuild26s          0.3374 -0.0357
    ## ParentalCare3p      -1.4585 -0.0847
    ## ParentalCare4n       0.8508  0.0494
    ## WaterPref1s          0.6020  0.1442
    ## WaterPref3a         -1.6855 -0.4036
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.3305  0.214
    ## DemersPelag      0.0000  1.000
    ## OperculumPresent 0.3061  0.207
    ## FeedingPath      0.1586  0.321
    ## RepGuild1        0.2986  0.203
    ## RepGuild2        0.2986  0.203
    ## ParentalCare     0.4003  1.000
    ## WaterPref        0.3449  0.276
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envea_efp <- p.adjust.envfit(envea_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
envea_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.86550 -0.50090 0.1250      1
    ## Troph            -0.99338 -0.11487 0.6303      1
    ## DepthMin          0.34970 -0.93686 0.0226      1
    ## DepthMax          0.98912 -0.14712 0.2076      1
    ## TempPrefMin       0.67062 -0.74180 0.1311      1
    ## TempPrefMax      -0.98765  0.15670 0.0121      1
    ## DorsalSpinesMean  0.98478 -0.17380 0.5093      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.2453 -0.9704
    ## BodyShapeI3f         0.5851  0.0547
    ## BodyShapeI4e        -1.7789  0.3075
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.8894 -0.0241
    ## OperculumPresentyes  0.5039  0.0064
    ## FeedingPathb         0.2181  0.1023
    ## FeedingPathp        -1.8541 -0.8695
    ## RepGuild11b         -1.9179 -1.8022
    ## RepGuild12g         -1.7405  1.1871
    ## RepGuild13n          0.3374 -0.0357
    ## RepGuild22eb        -1.9179 -1.8022
    ## RepGuild23n         -1.7405  1.1871
    ## RepGuild26s          0.3374 -0.0357
    ## ParentalCare3p      -1.4585 -0.0847
    ## ParentalCare4n       0.8508  0.0494
    ## WaterPref1s          0.6020  0.1442
    ## WaterPref3a         -1.6855 -0.4036
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.3305      1
    ## DemersPelag      0.0000      1
    ## OperculumPresent 0.3061      1
    ## FeedingPath      0.1586      1
    ## RepGuild1        0.2986      1
    ## RepGuild2        0.2986      1
    ## ParentalCare     0.4003      1
    ## WaterPref        0.3449      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
envems_ef <- envfit(envems_NMDS, FD_total_env[mixed_stratified_lakes,c(6:20)], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
envems_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)  
    ## MaxLengthTL       0.82417 -0.56634 0.1952  0.443  
    ## Troph            -0.98014 -0.19831 0.6410  0.033 *
    ## DepthMin          0.19540 -0.98072 0.0197  0.816  
    ## DepthMax          0.95523 -0.29585 0.2609  0.392  
    ## TempPrefMin       0.40652 -0.91364 0.1048  0.776  
    ## TempPrefMax      -0.99785  0.06549 0.0140  0.899  
    ## DorsalSpinesMean  0.95399 -0.29984 0.4452  0.266  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -1.9191 -1.6260
    ## BodyShapeI3f         0.7145 -0.0258
    ## BodyShapeI4e        -1.4851  0.4774
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.6458  0.1641
    ## OperculumPresentyes  0.5486 -0.0547
    ## FeedingPathb         0.2459  0.1016
    ## FeedingPathp        -1.7212 -0.7114
    ## RepGuild11b         -1.9191 -1.6260
    ## RepGuild12g         -1.2929  1.3501
    ## RepGuild13n          0.3465 -0.0826
    ## RepGuild22eb        -1.9191 -1.6260
    ## RepGuild23n         -1.2929  1.3501
    ## RepGuild26s          0.3465 -0.0826
    ## ParentalCare3p      -1.2052  0.0324
    ## ParentalCare4n       0.9374 -0.0252
    ## WaterPref1s          0.6760  0.1112
    ## WaterPref3a         -1.4871 -0.2446
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.4299  0.107
    ## DemersPelag      0.0000  1.000
    ## OperculumPresent 0.2892  0.165
    ## FeedingPath      0.1572  0.292
    ## RepGuild1        0.2967  0.182
    ## RepGuild2        0.2967  0.182
    ## ParentalCare     0.3586  0.877
    ## WaterPref        0.3275  0.219
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_efp <- p.adjust.envfit(envems_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
envems_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.82417 -0.56634 0.1952  1.000
    ## Troph            -0.98014 -0.19831 0.6410  0.495
    ## DepthMin          0.19540 -0.98072 0.0197  1.000
    ## DepthMax          0.95523 -0.29585 0.2609  1.000
    ## TempPrefMin       0.40652 -0.91364 0.1048  1.000
    ## TempPrefMax      -0.99785  0.06549 0.0140  1.000
    ## DorsalSpinesMean  0.95399 -0.29984 0.4452  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -1.9191 -1.6260
    ## BodyShapeI3f         0.7145 -0.0258
    ## BodyShapeI4e        -1.4851  0.4774
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.6458  0.1641
    ## OperculumPresentyes  0.5486 -0.0547
    ## FeedingPathb         0.2459  0.1016
    ## FeedingPathp        -1.7212 -0.7114
    ## RepGuild11b         -1.9191 -1.6260
    ## RepGuild12g         -1.2929  1.3501
    ## RepGuild13n          0.3465 -0.0826
    ## RepGuild22eb        -1.9191 -1.6260
    ## RepGuild23n         -1.2929  1.3501
    ## RepGuild26s          0.3465 -0.0826
    ## ParentalCare3p      -1.2052  0.0324
    ## ParentalCare4n       0.9374 -0.0252
    ## WaterPref1s          0.6760  0.1112
    ## WaterPref3a         -1.4871 -0.2446
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.4299      1
    ## DemersPelag      0.0000      1
    ## OperculumPresent 0.2892      1
    ## FeedingPath      0.1572      1
    ## RepGuild1        0.2967      1
    ## RepGuild2        0.2967      1
    ## ParentalCare     0.3586      1
    ## WaterPref        0.3275      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
enveom_ef <- envfit(enveom_NMDS, FD_total_env[ocean_mixed_sites_env,c(6:20)], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
enveom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)  
    ## MaxLengthTL       0.66486  0.74697 0.2866  0.299  
    ## Troph             0.97077  0.24002 0.4898  0.095 .
    ## DepthMin         -0.80338  0.59547 0.4173  0.104  
    ## DepthMax          0.13143  0.99133 0.3472  0.179  
    ## TempPrefMin       0.02081 -0.99978 0.1083  0.642  
    ## TempPrefMax       0.26410 -0.96449 0.3151  0.166  
    ## DorsalSpinesMean -0.91571  0.40185 0.1976  0.436  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.3823 -0.0896
    ## BodyShapeI3f         0.0382  0.0090
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentyes  0.0000  0.0000
    ## FeedingPathb         0.0000  0.0000
    ## RepGuild13n          0.0000  0.0000
    ## RepGuild26s          0.0000  0.0000
    ## ParentalCare4n       0.0000  0.0000
    ## WaterPref1s          0.0000  0.0000
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.0199      1
    ## DemersPelag      0.0000      1
    ## OperculumPresent 0.0000      1
    ## FeedingPath      0.0000      1
    ## RepGuild1        0.0000      1
    ## RepGuild2        0.0000      1
    ## ParentalCare     0.0000      1
    ## WaterPref        0.0000      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_efp <- p.adjust.envfit(enveom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
enveom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.66486  0.74697 0.2866      1
    ## Troph             0.97077  0.24002 0.4898      1
    ## DepthMin         -0.80338  0.59547 0.4173      1
    ## DepthMax          0.13143  0.99133 0.3472      1
    ## TempPrefMin       0.02081 -0.99978 0.1083      1
    ## TempPrefMax       0.26410 -0.96449 0.3151      1
    ## DorsalSpinesMean -0.91571  0.40185 0.1976      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.3823 -0.0896
    ## BodyShapeI3f         0.0382  0.0090
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentyes  0.0000  0.0000
    ## FeedingPathb         0.0000  0.0000
    ## RepGuild13n          0.0000  0.0000
    ## RepGuild26s          0.0000  0.0000
    ## ParentalCare4n       0.0000  0.0000
    ## WaterPref1s          0.0000  0.0000
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.0199      1
    ## DemersPelag      0.0000      1
    ## OperculumPresent 0.0000      1
    ## FeedingPath      0.0000      1
    ## RepGuild1        0.0000      1
    ## RepGuild2        0.0000      1
    ## ParentalCare     0.0000      1
    ## WaterPref        0.0000      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
enveso_ef <- envfit(enveso_NMDS, FD_total_env[ocean_stratified_sites_env,c(6:20)], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
enveso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL      -0.99999 -0.00346 0.1580  0.392
    ## Troph             0.79229  0.61015 0.7642  0.187
    ## DepthMin         -0.45791  0.88900 0.0146  0.834
    ## DepthMax         -0.79449 -0.60727 0.2841  0.208
    ## TempPrefMin      -0.76219  0.64736 0.1318  0.739
    ## TempPrefMax       0.55957  0.82879 0.1005  0.553
    ## DorsalSpinesMean -0.93134 -0.36415 0.5664  0.231
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.6152  0.5766
    ## BodyShapeI3f        -0.7241 -0.3008
    ## BodyShapeI4e         1.2127  0.0877
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno   1.1931  0.3631
    ## OperculumPresentyes -0.6818 -0.2075
    ## FeedingPathb        -0.1433 -0.1861
    ## FeedingPathp         0.6449  0.8373
    ## RepGuild11b          0.5971  1.6435
    ## RepGuild12g          1.8018 -0.3562
    ## RepGuild13n         -0.5251 -0.1164
    ## RepGuild22eb         0.5971  1.6435
    ## RepGuild23n          1.8018 -0.3562
    ## RepGuild26s         -0.5251 -0.1164
    ## ParentalCare3p       0.7000  0.2577
    ## ParentalCare4n      -1.2251 -0.4510
    ## WaterPref1s         -0.7459 -0.5313
    ## WaterPref3a          0.8951  0.6375
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.3098  0.167
    ## DemersPelag      0.0000  1.000
    ## OperculumPresent 0.2910  0.151
    ## FeedingPath      0.0813  0.452
    ## RepGuild1        0.3608  0.167
    ## RepGuild2        0.3608  0.167
    ## ParentalCare     0.3189  0.877
    ## WaterPref        0.3296  0.129
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveso_efp <- p.adjust.envfit(enveso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
enveso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL      -0.99999 -0.00346 0.1580      1
    ## Troph             0.79229  0.61015 0.7642      1
    ## DepthMin         -0.45791  0.88900 0.0146      1
    ## DepthMax         -0.79449 -0.60727 0.2841      1
    ## TempPrefMin      -0.76219  0.64736 0.1318      1
    ## TempPrefMax       0.55957  0.82879 0.1005      1
    ## DorsalSpinesMean -0.93134 -0.36415 0.5664      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.6152  0.5766
    ## BodyShapeI3f        -0.7241 -0.3008
    ## BodyShapeI4e         1.2127  0.0877
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno   1.1931  0.3631
    ## OperculumPresentyes -0.6818 -0.2075
    ## FeedingPathb        -0.1433 -0.1861
    ## FeedingPathp         0.6449  0.8373
    ## RepGuild11b          0.5971  1.6435
    ## RepGuild12g          1.8018 -0.3562
    ## RepGuild13n         -0.5251 -0.1164
    ## RepGuild22eb         0.5971  1.6435
    ## RepGuild23n          1.8018 -0.3562
    ## RepGuild26s         -0.5251 -0.1164
    ## ParentalCare3p       0.7000  0.2577
    ## ParentalCare4n      -1.2251 -0.4510
    ## WaterPref1s         -0.7459 -0.5313
    ## WaterPref3a          0.8951  0.6375
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.3098      1
    ## DemersPelag      0.0000      1
    ## OperculumPresent 0.2910      1
    ## FeedingPath      0.0813      1
    ## RepGuild1        0.3608      1
    ## RepGuild2        0.3608      1
    ## ParentalCare     0.3189      1
    ## WaterPref        0.3296      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# All sites for figure
envb_ef <- envfit(envb_NMDS, FD_total_env[surveyed_sites_geo,c(53,65)], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
envb_ef
```

    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##        NMDS1   NMDS2
    ## Ono  -1.2010 -0.9482
    ## Oyes  0.2826  0.2231
    ## W1s   0.3783  0.2816
    ## W3a  -1.2107 -0.9010
    ## 
    ## Goodness of fit:
    ##       r2 Pr(>r)  
    ## O 0.1993  0.130  
    ## W 0.2575  0.074 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envb_efp <- p.adjust.envfit(envb_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
envb_efp
```

    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##        NMDS1   NMDS2
    ## Ono  -1.2010 -0.9482
    ## Oyes  0.2826  0.2231
    ## W1s   0.3783  0.2816
    ## W3a  -1.2107 -0.9010
    ## 
    ## Goodness of fit:
    ##       r2 Pr(>r)
    ## O 0.1993  0.260
    ## W 0.2575  0.148
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# All sites by stratification
envba_ef <- envfit(envb_NMDS, FD_total_env[surveyed_sites_geo,c(6:20)], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
envba_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.36541  0.93085 0.0296  0.578
    ## Troph            -0.72002 -0.69395 0.6528  0.140
    ## DepthMin          0.50576  0.86268 0.0381  0.366
    ## DepthMax          0.68279  0.73062 0.1312  0.229
    ## TempPrefMin       0.44178  0.89713 0.1583  0.628
    ## TempPrefMax      -0.62465  0.78090 0.0324  0.508
    ## DorsalSpinesMean  0.83957  0.54326 0.2926  0.659
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.0946 -0.0199
    ## BodyShapeI3f         0.2535  0.1731
    ## BodyShapeI4e        -0.9584 -0.5911
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.2010 -0.9482
    ## OperculumPresentyes  0.2826  0.2231
    ## FeedingPathb         0.0951  0.1105
    ## FeedingPathp        -0.9035 -1.0499
    ## RepGuild11b         -1.4866 -1.5154
    ## RepGuild12g         -1.1040 -0.1254
    ## RepGuild13n          0.2052  0.0981
    ## RepGuild22eb        -1.4866 -1.5154
    ## RepGuild23n         -1.1040 -0.1254
    ## RepGuild26s          0.2052  0.0981
    ## ParentalCare3p      -0.9520 -0.8281
    ## ParentalCare4n       0.4760  0.4141
    ## WaterPref1s          0.3783  0.2816
    ## WaterPref3a         -1.2107 -0.9010
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## BodyShapeI       0.1106  0.813  
    ## DemersPelag      0.0000  1.000  
    ## OperculumPresent 0.1993  0.122  
    ## FeedingPath      0.0731  0.349  
    ## RepGuild1        0.1362  0.222  
    ## RepGuild2        0.1362  0.222  
    ## ParentalCare     0.2880  0.479  
    ## WaterPref        0.2575  0.062 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envba_efp <- p.adjust.envfit(envba_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
envba_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.36541  0.93085 0.0296      1
    ## Troph            -0.72002 -0.69395 0.6528      1
    ## DepthMin          0.50576  0.86268 0.0381      1
    ## DepthMax          0.68279  0.73062 0.1312      1
    ## TempPrefMin       0.44178  0.89713 0.1583      1
    ## TempPrefMax      -0.62465  0.78090 0.0324      1
    ## DorsalSpinesMean  0.83957  0.54326 0.2926      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.0946 -0.0199
    ## BodyShapeI3f         0.2535  0.1731
    ## BodyShapeI4e        -0.9584 -0.5911
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.2010 -0.9482
    ## OperculumPresentyes  0.2826  0.2231
    ## FeedingPathb         0.0951  0.1105
    ## FeedingPathp        -0.9035 -1.0499
    ## RepGuild11b         -1.4866 -1.5154
    ## RepGuild12g         -1.1040 -0.1254
    ## RepGuild13n          0.2052  0.0981
    ## RepGuild22eb        -1.4866 -1.5154
    ## RepGuild23n         -1.1040 -0.1254
    ## RepGuild26s          0.2052  0.0981
    ## ParentalCare3p      -0.9520 -0.8281
    ## ParentalCare4n       0.4760  0.4141
    ## WaterPref1s          0.3783  0.2816
    ## WaterPref3a         -1.2107 -0.9010
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.1106   1.00
    ## DemersPelag      0.0000   1.00
    ## OperculumPresent 0.1993   1.00
    ## FeedingPath      0.0731   1.00
    ## RepGuild1        0.1362   1.00
    ## RepGuild2        0.1362   1.00
    ## ParentalCare     0.2880   1.00
    ## WaterPref        0.2575   0.93
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
envbms_ef <- envfit(envbms_NMDS, FD_total_env[mixed_stratified_lakes_geo,c(6:20)], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
envbms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.63646  0.77131 0.0916  0.608
    ## Troph            -0.85662 -0.51595 0.4192  0.200
    ## DepthMin          0.72151  0.69240 0.0646  0.707
    ## DepthMax          0.61982  0.78474 0.3384  0.130
    ## TempPrefMin       0.83214  0.55457 0.0792  0.719
    ## TempPrefMax      -0.10863 -0.99408 0.0869  0.472
    ## DorsalSpinesMean  0.64457  0.76454 0.1359  0.781
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -2.0046 -0.1949
    ## BodyShapeI3f         0.4846  0.0705
    ## BodyShapeI4e        -0.7104 -0.1276
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.2423 -0.2025
    ## OperculumPresentyes  0.4518  0.0736
    ## FeedingPathb         0.1841 -0.0304
    ## FeedingPathp        -1.1967  0.1973
    ## RepGuild11b         -2.0046 -0.1949
    ## RepGuild12g         -0.2485 -0.5350
    ## RepGuild13n          0.2085  0.1054
    ## RepGuild22eb        -2.0046 -0.1949
    ## RepGuild23n         -0.2485 -0.5350
    ## RepGuild26s          0.2085  0.1054
    ## ParentalCare3p      -0.9622  0.0321
    ## ParentalCare4n       0.8420 -0.0281
    ## WaterPref1s          0.5876  0.1076
    ## WaterPref3a         -1.1753 -0.2152
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## BodyShapeI       0.3291  0.245  
    ## DemersPelag      0.0000  1.000  
    ## OperculumPresent 0.3330  0.085 .
    ## FeedingPath      0.1308  0.285  
    ## RepGuild1        0.2084  0.361  
    ## RepGuild2        0.2084  0.361  
    ## ParentalCare     0.4688  0.132  
    ## WaterPref        0.4126  0.070 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbms_efp <- p.adjust.envfit(envbms_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
envbms_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.63646  0.77131 0.0916      1
    ## Troph            -0.85662 -0.51595 0.4192      1
    ## DepthMin          0.72151  0.69240 0.0646      1
    ## DepthMax          0.61982  0.78474 0.3384      1
    ## TempPrefMin       0.83214  0.55457 0.0792      1
    ## TempPrefMax      -0.10863 -0.99408 0.0869      1
    ## DorsalSpinesMean  0.64457  0.76454 0.1359      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -2.0046 -0.1949
    ## BodyShapeI3f         0.4846  0.0705
    ## BodyShapeI4e        -0.7104 -0.1276
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.2423 -0.2025
    ## OperculumPresentyes  0.4518  0.0736
    ## FeedingPathb         0.1841 -0.0304
    ## FeedingPathp        -1.1967  0.1973
    ## RepGuild11b         -2.0046 -0.1949
    ## RepGuild12g         -0.2485 -0.5350
    ## RepGuild13n          0.2085  0.1054
    ## RepGuild22eb        -2.0046 -0.1949
    ## RepGuild23n         -0.2485 -0.5350
    ## RepGuild26s          0.2085  0.1054
    ## ParentalCare3p      -0.9622  0.0321
    ## ParentalCare4n       0.8420 -0.0281
    ## WaterPref1s          0.5876  0.1076
    ## WaterPref3a         -1.1753 -0.2152
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.3291      1
    ## DemersPelag      0.0000      1
    ## OperculumPresent 0.3330      1
    ## FeedingPath      0.1308      1
    ## RepGuild1        0.2084      1
    ## RepGuild2        0.2084      1
    ## ParentalCare     0.4688      1
    ## WaterPref        0.4126      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
envbom_ef <- envfit(envbom_NMDS, FD_total_env[ocean_mixed_sites_geo,c(6:20)], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
envbom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)   
    ## MaxLengthTL       0.14581  0.98931 0.3467  0.151   
    ## Troph             0.38556  0.92268 0.4742  0.291   
    ## DepthMin         -0.40134  0.91593 0.1366  0.332   
    ## DepthMax         -0.04748  0.99887 0.4346  0.098 . 
    ## TempPrefMin       0.03417 -0.99942 0.6417  0.007 **
    ## TempPrefMax       0.12417 -0.99226 0.3988  0.043 * 
    ## DorsalSpinesMean -0.73989  0.67273 0.1931  0.356   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.1541 -0.3693
    ## BodyShapeI3f         0.0280  0.0671
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentyes  0.0000  0.0000
    ## FeedingPathb         0.0000  0.0000
    ## RepGuild13n          0.0000  0.0000
    ## RepGuild26s          0.0000  0.0000
    ## ParentalCare4n       0.0000  0.0000
    ## WaterPref1s          0.0000  0.0000
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.0101   0.94
    ## DemersPelag      0.0000   1.00
    ## OperculumPresent 0.0000   1.00
    ## FeedingPath      0.0000   1.00
    ## RepGuild1        0.0000   1.00
    ## RepGuild2        0.0000   1.00
    ## ParentalCare     0.0000   1.00
    ## WaterPref        0.0000   1.00
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbom_efp <- p.adjust.envfit(envbom_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
envbom_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.14581  0.98931 0.3467  1.000
    ## Troph             0.38556  0.92268 0.4742  1.000
    ## DepthMin         -0.40134  0.91593 0.1366  1.000
    ## DepthMax         -0.04748  0.99887 0.4346  1.000
    ## TempPrefMin       0.03417 -0.99942 0.6417  0.105
    ## TempPrefMax       0.12417 -0.99226 0.3988  0.645
    ## DorsalSpinesMean -0.73989  0.67273 0.1931  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.1541 -0.3693
    ## BodyShapeI3f         0.0280  0.0671
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentyes  0.0000  0.0000
    ## FeedingPathb         0.0000  0.0000
    ## RepGuild13n          0.0000  0.0000
    ## RepGuild26s          0.0000  0.0000
    ## ParentalCare4n       0.0000  0.0000
    ## WaterPref1s          0.0000  0.0000
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.0101      1
    ## DemersPelag      0.0000      1
    ## OperculumPresent 0.0000      1
    ## FeedingPath      0.0000      1
    ## RepGuild1        0.0000      1
    ## RepGuild2        0.0000      1
    ## ParentalCare     0.0000      1
    ## WaterPref        0.0000      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
envbso_ef <- envfit(envbso_NMDS, FD_total_env[ocean_stratified_sites,c(6:20)], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
envbso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)  
    ## MaxLengthTL      -0.65161 -0.75855 0.0682  0.336  
    ## Troph             0.63255  0.77452 0.7239  0.213  
    ## DepthMin         -0.37416 -0.92737 0.0574  0.310  
    ## DepthMax         -0.81515 -0.57925 0.1963  0.084 .
    ## TempPrefMin      -0.67395 -0.73878 0.0963  0.910  
    ## TempPrefMax       0.98595  0.16706 0.0365  0.558  
    ## DorsalSpinesMean -0.91169 -0.41087 0.3870  0.584  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.0038 -0.2135
    ## BodyShapeI3f        -0.6311 -0.0782
    ## BodyShapeI4e         1.1073  0.2970
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno   1.3627  0.6423
    ## OperculumPresentyes -0.5451 -0.2569
    ## FeedingPathb        -0.1750 -0.1201
    ## FeedingPathp         1.0503  0.7208
    ## RepGuild11b          1.6231  1.1809
    ## RepGuild12g          1.2503 -0.1070
    ## RepGuild13n         -0.3749 -0.0879
    ## RepGuild22eb         1.6231  1.1809
    ## RepGuild23n          1.2503 -0.1070
    ## RepGuild26s         -0.3749 -0.0879
    ## ParentalCare3p       1.0986  0.5136
    ## ParentalCare4n      -1.0986 -0.5136
    ## WaterPref1s         -0.7500 -0.3350
    ## WaterPref3a          1.3501  0.6030
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## BodyShapeI       0.1663  0.744  
    ## DemersPelag      0.0000  1.000  
    ## OperculumPresent 0.2570  0.083 .
    ## FeedingPath      0.0766  0.405  
    ## RepGuild1        0.1781  0.195  
    ## RepGuild2        0.1781  0.195  
    ## ParentalCare     0.4163  0.518  
    ## WaterPref        0.3438  0.087 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envbso_efp <- p.adjust.envfit(envbso_ef)
```

    ## Adjustment of significance by bonferroni method

``` r
envbso_efp
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL      -0.65161 -0.75855 0.0682      1
    ## Troph             0.63255  0.77452 0.7239      1
    ## DepthMin         -0.37416 -0.92737 0.0574      1
    ## DepthMax         -0.81515 -0.57925 0.1963      1
    ## TempPrefMin      -0.67395 -0.73878 0.0963      1
    ## TempPrefMax       0.98595  0.16706 0.0365      1
    ## DorsalSpinesMean -0.91169 -0.41087 0.3870      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.0038 -0.2135
    ## BodyShapeI3f        -0.6311 -0.0782
    ## BodyShapeI4e         1.1073  0.2970
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno   1.3627  0.6423
    ## OperculumPresentyes -0.5451 -0.2569
    ## FeedingPathb        -0.1750 -0.1201
    ## FeedingPathp         1.0503  0.7208
    ## RepGuild11b          1.6231  1.1809
    ## RepGuild12g          1.2503 -0.1070
    ## RepGuild13n         -0.3749 -0.0879
    ## RepGuild22eb         1.6231  1.1809
    ## RepGuild23n          1.2503 -0.1070
    ## RepGuild26s         -0.3749 -0.0879
    ## ParentalCare3p       1.0986  0.5136
    ## ParentalCare4n      -1.0986 -0.5136
    ## WaterPref1s         -0.7500 -0.3350
    ## WaterPref3a          1.3501  0.6030
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.1663      1
    ## DemersPelag      0.0000      1
    ## OperculumPresent 0.2570      1
    ## FeedingPath      0.0766      1
    ## RepGuild1        0.1781      1
    ## RepGuild2        0.1781      1
    ## ParentalCare     0.4163      1
    ## WaterPref        0.3438      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

## Environment/biogeography mantel tests

- We used mantel tests to determine significance between the
  environmental/Geographical and trait distance matrices.

``` r
### Environmental
# All sites
FDe_dist <- gowdis(as.data.frame(scaled_FD_total_env[surveyed_sites_env,c(1:17)]))
envea_mant <- mantel(enve_dist, FDe_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19])
envea_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = enve_dist, ydis = FDe_dist, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6119 
    ##       Significance: 0.164 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.630 0.658 0.686 0.701 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
FDems_dist <- gowdis(scaled_FD_total_env[mixed_stratified_lakes,c(1:17)])
envems_mant <- mantel(envems_dist, FDems_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19])
envems_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = envems_dist, ydis = FDems_dist, method = "spearman",      permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5612 
    ##       Significance: 0.103 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.562 0.597 0.622 0.649 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
FDeom_dist <- gowdis(scaled_FD_total_env[ocean_mixed_sites_env,c(1:7)])
enveom_mant <- mantel(enveom_dist, FDeom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19])
enveom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = enveom_dist, ydis = FDeom_dist, method = "spearman",      permutations = 999, strata = env[ocean_mixed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2326 
    ##       Significance: 0.136 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.280 0.392 0.447 0.526 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
FDeso_dist <- gowdis(scaled_FD_total_env[ocean_stratified_sites_env,c(1:17)])
enveso_mant <- mantel(enveso_dist, FDeso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19])
enveso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = enveso_dist, ydis = FDeso_dist, method = "spearman",      permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2606 
    ##       Significance: 0.195 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.316 0.364 0.395 0.434 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjusted p-values
enve_mant_p <- rbind(envems_mant$signif, enveom_mant$signif, enveso_mant$signif)
enve_mant_p <- enve_mant_p[,1]
enve_mant_p <- p.adjust(enve_mant_p, method = "bonferroni")
enve_mant_p
```

    ## [1] 0.309 0.408 0.585

``` r
### Geographic
# All sites for figure
FDb_dist <- gowdis(scaled_FD_total_env[surveyed_sites_geo,c(1:17)])
envba_mant <- mantel(envb_dist, FDb_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19])
envba_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = envb_dist, ydis = FDb_dist, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1251 
    ##       Significance: 0.481 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.227 0.254 0.280 0.318 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
FDbms_dist <- gowdis(scaled_FD_total_env[mixed_stratified_lakes_geo,c(1:17)])
envbms_mant <- mantel(envbms_dist, FDbms_dist, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19])
envbms_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = envbms_dist, ydis = FDbms_dist, method = "spearman",      permutations = 1000, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1603 
    ##       Significance: 0.45554 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.358 0.416 0.450 0.495 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
FDbom_dist <- gowdis(scaled_FD_total_env[ocean_mixed_sites_geo,c(1:7)])
envbom_mant <- mantel(envbom_dist, FDbom_dist, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19])
envbom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = envbom_dist, ydis = FDbom_dist, method = "spearman",      permutations = 1000, strata = env[ocean_mixed_sites_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1571 
    ##       Significance: 0.05994 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.116 0.170 0.202 0.241 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
FDbso_dist <- gowdis(scaled_FD_total_env[ocean_stratified_sites,c(1:17)])
envbso_mant <- mantel(envbso_dist, FDbso_dist, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[ocean_stratified_sites,19])
envbso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = envbso_dist, ydis = FDbso_dist, method = "spearman",      permutations = 1000, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.003759 
    ##       Significance: 0.42557 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.109 0.136 0.171 0.188 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Adjusted p-values
envb_mant_p <- rbind(envbms_mant$signif, envbom_mant$signif, envbso_mant$signif)
envb_mant_p <- envb_mant_p[,1]
envb_mant_p <- p.adjust(envb_mant_p, method = "bonferroni")
envb_mant_p
```

    ## [1] 1.0000000 0.1798202 1.0000000

## Plot environmental/geographical NMDS and envfit results

- Includes a plot of correlated trait variables

``` r
enve_NMDS_data.scores <- as.data.frame(scores(enve_NMDS))
enve_NMDS_data.scores$Stratification <- env[surveyed_sites_env,19]
enve_NMDS_data.scores$Lakes <- env[surveyed_sites_env,1]
enve_NMDS_data.scores$Stratification <- factor(enve_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# enve_ef_coord_cont <- as.data.frame(scores(enve_ef, "vectors")) * ordiArrowMul(enve_ef)
# enve_ef_coord_cat <- as.data.frame(scores(enve_ef, "factors")) * ordiArrowMul(enve_ef)

enve_ef_plot <- ggplot(data = enve_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification),type="t",level = 0.95) +
  geom_point(data = enve_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  geom_text_repel(data = enve_NMDS_data.scores, label = enve_NMDS_data.scores$Lakes, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              data = enve_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#CD3333") +
  # geom_text(data = enve_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
  #           color = "#CD3333", label = row.names(enve_ef_coord_cont), size = 7) + 
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              data = enve_ef_coord_cat, linewidth =1, alpha = 0.2, color = "#8B6508") +
  # geom_text(data = enve_ef_coord_cat, aes(x = NMDS1, y = NMDS2),
  #           color = "#8B6508", label = row.names(enve_ef_coord_cat), size = 7, check_overlap = T) +
  theme(text = element_text(size = 22), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16), panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  annotate("text", x = 1.9, y =3.4, size = 5, 
           label = paste("Stress: ", round(enve_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "A")
enve_ef_plot <- enve_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
enve_ef_plot
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_path()`).

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20environmental/Geographical%20NMDS%20and%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/enve_ef_plot.png", enve_ef_plot, width = 6, height = 4, units = "in")
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_path()`).

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/enve_NMDS_boxplot.tiff", width = 1200, height = 800, res = 150, type = "cairo")
enve_NMDS_data.scores$NMDS <- sqrt((enve_NMDS_data.scores$NMDS1)^2 + (enve_NMDS_data.scores$NMDS2)^2)
ordered_enve_NMDS_data.scores <- enve_NMDS_data.scores[order(enve_NMDS_data.scores$NMDS), ]
Q1 <- ordered_enve_NMDS_data.scores[2,1]
Q3 <- ordered_enve_NMDS_data.scores[5,1]
IQR <- Q3 - Q1
Q1 - 1.5*IQR
```

    ## [1] -0.1354796

``` r
Q3 + 1.5*IQR
```

    ## [1] 1.938206

``` r
ordered_enve_NMDS_data.scores$NMDS
```

    ##  [1] 0.2186056 0.6652153 0.7492830 0.9176676 1.1827237 1.2628071 1.3623662
    ##  [8] 1.4339286 1.5426043 1.5504106 1.7241463 1.7914005 1.8206323 1.9133145
    ## [15] 1.9770652 2.1250499 2.2042851 2.6317930 3.4846341

``` r
boxplot(NMDS ~ Stratification, ordered_enve_NMDS_data.scores)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
envb_NMDS_data.scores <- as.data.frame(scores(envb_NMDS))
envb_NMDS_data.scores$Stratification <- env[surveyed_sites_geo,19]
envb_NMDS_data.scores$Lakes <- env[surveyed_sites_geo,1]
envb_NMDS_data.scores$Stratification <- factor(envb_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# envb_ef_coord_cont <- as.data.frame(scores(envb_ef, "vectors")) * ordiArrowMul(envb_ef)
# envb_ef_coord_cat = as.data.frame(scores(envb_ef, "factors")) * ordiArrowMul(envb_ef)

envb_ef_plot <- ggplot(data = envb_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification),type="t",level = 0.95) +
  geom_point(data = envb_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 3, alpha = 1) + 
  geom_text_repel(data = envb_NMDS_data.scores, label = envb_NMDS_data.scores$Lakes, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              data = envb_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#CD3333") +
  # geom_text(data = envb_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
  #           color = "#CD3333", label = row.names(envb_ef_coord_cont), size = 7) + 
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
  #              data = envb_ef_coord_cat, linewidth =1, alpha = 0.2, color = "#8B6508") +
  # geom_text(data = envb_ef_coord_cat, aes(x = NMDS1, y = NMDS2),
  #           color = "#8B6508", label = row.names(envb_ef_coord_cat), size = 7, check_overlap = T) +
  theme(text = element_text(size = 22), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16), panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  annotate("text", x = 4.5, y = 2.9, size = 5, 
           label = paste("Stress: ", round(envb_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "B")
envb_ef_plot <- envb_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
envb_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/envb_ef_plot.png", envb_ef_plot, width = 6, height = 4, units = "in")

png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/envb_NMDS_boxplot.tiff", width = 1200, height = 800, res = 150, type = "cairo")
envb_NMDS_data.scores$NMDS <- sqrt((envb_NMDS_data.scores$NMDS1)^2 + (envb_NMDS_data.scores$NMDS2)^2)
ordered_envb_NMDS_data.scores <- envb_NMDS_data.scores[order(envb_NMDS_data.scores$NMDS), ]
Q1 <- ordered_envb_NMDS_data.scores[2,1]
Q3 <- ordered_envb_NMDS_data.scores[5,1]
IQR <- Q3 - Q1
Q1 - 1.5*IQR
```

    ## [1] 0.9891062

``` r
Q3 + 1.5*IQR
```

    ## [1] -1.483177

``` r
ordered_envb_NMDS_data.scores$NMDS
```

    ##  [1] 0.2545109 0.2897713 0.5138708 0.5232713 0.5703586 0.6664421 0.7646071
    ##  [8] 0.8314242 1.0801529 1.1064503 1.1103271 1.2476376 1.3676553 1.5232544
    ## [15] 1.6044035 1.6998057 1.8664589 2.0104558 2.1228234 2.2712573 4.8254753

``` r
boxplot(NMDS ~ Stratification, ordered_envb_NMDS_data.scores)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## EDis and BDis Dendrogram

``` r
# cluster communities using average-linkage algorithm
ed_dist_clust <- hclust(enve_dist, method = "average")

# Open a PNG device
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/ed_dendro.tiff", width = 1200, height = 800, res = 150, type = "cairo")

# Create your plot using the plot() function
plot(ed_dist_clust, 
     xlab = "Sites",
     ylab = "Environmental dissimilarity",
     main = "",
     sub = "",
     pch = 20,
     col= "black")

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/ed_si.tiff", width = 1200, height = 800, res = 150, type = "cairo")

Si <- numeric(nrow(presabs_lake[surveyed_sites_env,]))
for (k in 2:(nrow(presabs_lake[surveyed_sites_env,])-1))
{
  sil<-silhouette(cutree(ed_dist_clust, k=k), enve_dist)
  Si[k] <- summary(sil)$avg.width
}
k.best<-which.max(Si)
plot(1:nrow(presabs_lake[surveyed_sites_env,]), Si, type = "h", main = "Silhouette", xlab = "K", ylab = "Width")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col = "red", font = 2, col.axis = "red")

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/ed_pam.tiff", width = 1200, height = 800, res = 150, type = "cairo")

enve_pam <- pam(enve_dist,3)
plot(enve_pam)

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# cluster communities using average-linkage algorithm
bd_dist_clust <- hclust(envb_dist, method = "average")

# Open a PNG device
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/bd_dendro.tiff", width = 1200, height = 800, res = 150, type = "cairo")

# Create your plot using the plot() function
plot(bd_dist_clust, 
     xlab = "Sites",
     ylab = "Geographic dissimilarity",
     main = "",
     sub = "",
     pch = 20,
     col= "black")

 # Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/bd_si.tiff", width = 1200, height = 800, res = 150, type = "cairo")

Si <- numeric(nrow(presabs_lake[surveyed_sites_geo,]))
for (k in 2:(nrow(presabs_lake[surveyed_sites_geo,])-1))
{
  sil<-silhouette(cutree(bd_dist_clust, k=k), envb_dist)
  Si[k] <- summary(sil)$avg.width
}
k.best<-which.max(Si)
plot(1:nrow(presabs_lake[surveyed_sites_geo,]), Si, type = "h", main = "Silhouette", xlab = "K", ylab = "Width")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col = "red", font = 2, col.axis = "red")

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/bd_pam.tiff", width = 1200, height = 800, res = 150, type = "cairo")

envb_pam <- pam(envb_dist,3)
plot(envb_pam)

# Close the PNG device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] rms_6.8-1            Hmisc_5.1-1          car_3.1-2           
    ##  [4] carData_3.0-5        adespatial_0.3-23    picante_1.8.2       
    ##  [7] nlme_3.1-164         pairwiseAdonis_0.4.1 cluster_2.1.6       
    ## [10] FD_1.0-12.3          geometry_0.4.7       ade4_1.7-22         
    ## [13] vegan_2.6-4          lattice_0.22-5       permute_0.9-7       
    ## [16] ggrepel_0.9.4        viridis_0.6.4        viridisLite_0.4.2   
    ## [19] ggplot2_3.5.1        tidyr_1.3.0          phytools_2.0-3      
    ## [22] maps_3.4.1.1         ape_5.7-1            reshape2_1.4.4      
    ## [25] stringr_1.5.1        dplyr_1.1.4          knitr_1.45          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      wk_0.9.1                rstudioapi_0.15.0      
    ##   [4] magrittr_2.0.3          TH.data_1.1-2           farver_2.1.1           
    ##   [7] rmarkdown_2.25          adegraphics_1.0-21      ragg_1.2.6             
    ##  [10] vctrs_0.6.5             spdep_1.3-1             base64enc_0.1-3        
    ##  [13] polspline_1.1.25        htmltools_0.5.7         progress_1.2.3         
    ##  [16] s2_1.1.5                Formula_1.2-5           adegenet_2.1.10        
    ##  [19] spData_2.3.0            KernSmooth_2.23-22      htmlwidgets_1.6.4      
    ##  [22] adephylo_1.1-16         sandwich_3.1-0          plyr_1.8.9             
    ##  [25] zoo_1.8-12              uuid_1.1-1              igraph_1.5.1           
    ##  [28] mime_0.12               lifecycle_1.0.4         iterators_1.0.14       
    ##  [31] pkgconfig_2.0.3         Matrix_1.6-4            R6_2.5.1               
    ##  [34] fastmap_1.1.1           shiny_1.8.0             magic_1.6-1            
    ##  [37] digest_0.6.33           numDeriv_2016.8-1.1     colorspace_2.1-0       
    ##  [40] phylobase_0.8.10        textshaping_0.3.7       labeling_0.4.3         
    ##  [43] clusterGeneration_1.3.8 fansi_1.0.6             httr_1.4.7             
    ##  [46] abind_1.4-5             mgcv_1.9-0              compiler_4.3.1         
    ##  [49] proxy_0.4-27            withr_2.5.2             doParallel_1.0.17      
    ##  [52] backports_1.4.1         htmlTable_2.4.2         optimParallel_1.0-2    
    ##  [55] DBI_1.1.3               highr_0.10              quantreg_5.98          
    ##  [58] MASS_7.3-60             classInt_0.4-10         scatterplot3d_0.3-44   
    ##  [61] tools_4.3.1             units_0.8-5             foreign_0.8-86         
    ##  [64] rncl_0.8.7              httpuv_1.6.13           nnet_7.3-19            
    ##  [67] glue_1.6.2              quadprog_1.5-8          promises_1.2.1         
    ##  [70] grid_4.3.1              sf_1.0-16               checkmate_2.3.1        
    ##  [73] generics_0.1.3          seqinr_4.2-36           gtable_0.3.4           
    ##  [76] class_7.3-22            data.table_1.14.10      hms_1.1.3              
    ##  [79] sp_2.1-2                xml2_1.3.6              utf8_1.2.4             
    ##  [82] foreach_1.5.2           pillar_1.9.0            later_1.3.2            
    ##  [85] splines_4.3.1           survival_3.5-7          deldir_2.0-2           
    ##  [88] SparseM_1.83            tidyselect_1.2.0        gridExtra_2.3          
    ##  [91] xfun_0.41               expm_0.999-8            stringi_1.8.2          
    ##  [94] yaml_2.3.7              boot_1.3-28.1           evaluate_0.23          
    ##  [97] codetools_0.2-19        interp_1.1-5            tibble_3.2.1           
    ## [100] cli_3.6.1               rpart_4.1.23            xtable_1.8-4           
    ## [103] systemfonts_1.0.5       munsell_0.5.0           Rcpp_1.0.11            
    ## [106] coda_0.19-4             png_0.1-8               XML_3.99-0.16          
    ## [109] MatrixModels_0.5-3      ellipsis_0.3.2          RNeXML_2.4.11          
    ## [112] prettyunits_1.2.0       latticeExtra_0.6-30     jpeg_0.1-10            
    ## [115] phangorn_2.11.1         mvtnorm_1.2-4           scales_1.3.0           
    ## [118] e1071_1.7-14            purrr_1.0.2             crayon_1.5.2           
    ## [121] combinat_0.0-8          rlang_1.1.2             multcomp_1.4-25        
    ## [124] fastmatch_1.1-4         mnormt_2.1.1

## Extra unused code

``` r
lakes <- traits_sesmpd_env[,1]
#traits_sesmpd_env$obs.z <- traits_sesmpd_env$mpd.obs.z
t_mpd.obs.z <- traits_sesmpd_env[,7]

#traits_sesmntd_env$obs.z <- traits_sesmntd_env$mntd.obs.z
t_mntd.obs.z <- traits_sesmntd_env[,7]

#traits_sespd_env$obs.z <- traits_sespd_env$pd.obs.z
t_pd.obs.z <- traits_sespd_env[,7]

#phylo_sesmpd_env$obs.z <- phylo_sesmpd_env$mpd.obs.z
p_mpd.obs.z <- phylo_sesmpd_env[,7]

#phylo_sesmntd_env$obs.z <- phylo_sesmntd_env$mntd.obs.z
p_mntd.obs.z <- phylo_sesmntd_env[,7]

#phylo_sespd_env$obs.z <- phylo_sespd_env$pd.obs.z
p_pd.obs.z <- phylo_sespd_env[,7]

p_sespd_env <- phylo_sespd_env[,c(10:40)]

tp_z_env <- cbind(lakes, t_mpd.obs.z, t_mntd.obs.z, t_pd.obs.z, p_mpd.obs.z, p_mntd.obs.z, p_pd.obs.z, p_sespd_env)
row.names(tp_z_env) <- tp_z_env$lakes


lakes_sub <- straits_sesmpd_env[,1]
#straits_sesmpd_env$obs.z <- straits_sesmpd_env$mpd.obs.z
tsub_mpd.obs.z <- straits_sesmpd_env[,7]

#straits_sesmntd_env$obs.z <- straits_sesmntd_env$mntd.obs.z
tsub_mntd.obs.z <- straits_sesmntd_env[,7]

#straits_sespd_env$obs.z <- straits_sespd_env$pd.obs.z
tsub_pd.obs.z <- straits_sespd_env[,7]

#sphylo_sesmpd_env$obs.z <- sphylo_sesmpd_env$mpd.obs.z
psub_mpd.obs.z <- sphylo_sesmpd_env[,7]

#sphylo_sesmntd_env$obs.z <- sphylo_sesmntd_env$mntd.obs.z
psub_mntd.obs.z <- sphylo_sesmntd_env[,7]

#sphylo_sespd_env$obs.z <- sphylo_sespd_env$pd.obs.z
psub_pd.obs.z <- sphylo_sespd_env[,7]

psub_sespd_env <- sphylo_sespd_env[,c(10:40)]

stp_z_env <- cbind(lakes_sub, tsub_mpd.obs.z, tsub_mntd.obs.z, tsub_pd.obs.z, psub_mpd.obs.z, psub_mntd.obs.z, psub_pd.obs.z, psub_sespd_env)
row.names(stp_z_env) <- stp_z_env$lakes_sub

t_z_dist <- vegdist(tp_z_env[,c(2:4)], method = "euclidean")
p_z_dist <- vegdist(tp_z_env[,c(5:7)], method = "euclidean")

st_z_dist <- vegdist(stp_z_env[,c(2:4)], method = "euclidean")
sp_z_dist <- vegdist(stp_z_env[,c(5:7)], method = "euclidean")

groups <- env[,19]
t_z_bd <- betadisper(t_z_dist, groups)
anova(t_z_bd)
permutest(t_z_bd, pairwise = TRUE)
(t_z_bd.HSD <- TukeyHSD(t_z_bd))
boxplot(t_z_bd)

groups <- env[-c(20),19]
st_z_bd <- betadisper(st_z_dist, groups)
anova(st_z_bd)
permutest(st_z_bd, pairwise = TRUE)
(st_z_bd.HSD <- TukeyHSD(st_z_bd))
boxplot(st_z_bd)
```
