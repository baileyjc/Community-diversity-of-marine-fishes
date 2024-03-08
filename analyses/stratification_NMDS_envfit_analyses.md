Stratification NMDS envfit analyses
================

### R Markdown

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
#Graphing
library(ggplot2)
library(viridis)
library(ggrepel)
```

## Load stratification NMDS envfit analyses functions

Need to run this function for downstream analyses in this Markdown file

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
jpres_dist <- vegdist(presabs_lake[-c(20,24:26),], method = "jaccard", binary = TRUE)
# Mixed and stratified lakes
jpresms_dist <- vegdist(presabs_lake[c(1:6,8,10,12:15,17,21:23),], method = "jaccard", binary = TRUE)
# Ocean sites and mixed lakes
jpresom_dist <- vegdist(presabs_lake[c(3,6:11,13:16,18:19,23),], method = "jaccard", binary = TRUE)
# Stratified lakes and ocean sites
jpresso_dist <- vegdist(presabs_lake[c(1:2,4,5,7,9,11:12,16:19,21:22),], method = "jaccard", binary = TRUE)
# Stratified lakes
jpress_dist <- vegdist(presabs_lake[c(1:2,4,5,12,17,21:22),], method = "jaccard", binary = TRUE)

## Dice-Sørensen
# All sites
spres_dist <- vegdist(presabs_lake[-c(20,24:26),], method = "bray", binary = TRUE)
# Mixed and stratified lakes
spresms_dist <- vegdist(presabs_lake[c(1:6,8,10,12:15,17,21:23),], method = "bray", binary = TRUE)
# Ocean sites and mixed lakes
spresom_dist <- vegdist(presabs_lake[c(3,6:11,13:16,18:19,23),], method = "bray", binary = TRUE)
# Stratified lakes and ocean sites
spresso_dist <- vegdist(presabs_lake[c(1:2,4,5,7,9,11:12,16:19,21:22),], method = "bray", binary = TRUE)


### Environmental -c(9,16,18,20)
## Jaccard
# All sites
jprese_dist <- vegdist(presabs_lake[-c(9,16,18,20,24:26),], method = "jaccard", binary = TRUE)
# Mixed and stratified lakes
jpresems_dist <- vegdist(presabs_lake[c(1:6,8,10,12:15,17,21:23),], method = "jaccard", binary = TRUE)
# Ocean sites and mixed lakes
jpreseom_dist <- vegdist(presabs_lake[c(3,6:8,10:11,13:15,19,23),], method = "jaccard", binary = TRUE)
# Stratified lakes and ocean sites
jpreseso_dist <- vegdist(presabs_lake[c(1:2,4,5,7,11:12,17,19,21:22),], method = "jaccard", binary = TRUE)

## Dice-Sørensen
# All sites
sprese_dist <- vegdist(presabs_lake[-c(9,16,18,20,24:26),], method = "bray", binary = TRUE)
# Mixed and stratified lakes
spresems_dist <- vegdist(presabs_lake[c(1:6,8,10,12:15,17,21:23),], method = "bray", binary = TRUE)
# Ocean sites and mixed lakes
spreseom_dist <- vegdist(presabs_lake[c(3,6:8,10:11,13:15,19,23),], method = "bray", binary = TRUE)
# Stratified lakes and ocean sites
spreseso_dist <- vegdist(presabs_lake[c(1:2,4,5,7,11:12,17,19,21:22),], method = "bray", binary = TRUE)


### Biogeographic -c(8,20)
## Jaccard
# All sites
jpresb_dist <- vegdist(presabs_lake[-c(8,20,24:26),], method = "jaccard", binary = TRUE)
# Mixed and stratified lakes
jpresbms_dist <- vegdist(presabs_lake[c(1:6,10,12:15,17,21:23),], method = "jaccard", binary = TRUE)
# Ocean sites and mixed lakes
jpresbom_dist <- vegdist(presabs_lake[c(3,6:7,9:11,13:16,18:19,23),], method = "jaccard", binary = TRUE)
# Stratified lakes and ocean sites
jpresbso_dist <- vegdist(presabs_lake[c(1:2,4,5,7,9,11:12,16:19,21:22),], method = "jaccard", binary = TRUE)

## Dice-Sørensen
# All sites
spresb_dist <- vegdist(presabs_lake[-c(8,20,24:26),], method = "bray", binary = TRUE)
# Mixed and stratified lakes
spresbms_dist <- vegdist(presabs_lake[c(1:6,10,12:15,17,21:23),], method = "bray", binary = TRUE)
# Ocean sites and mixed lakes
spresbom_dist <- vegdist(presabs_lake[c(3,6:7,9:11,13:16,18:19,23),], method = "bray", binary = TRUE)
# Stratified lakes and ocean sites
spresbso_dist <- vegdist(presabs_lake[c(1:2,4,5,7,9,11:12,16:19,21:22),], method = "bray", binary = TRUE)
```

## Incidence NMDS

To constrain dissimilarities we perform Nonmetric Multidimensional
Scaling (NMDS), which tries to find a stable solution using the metaMDS
package.

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
### Biogeographic
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

We use betadisper to determine homogeneity of the presence absence data
based on the stratification category. Is the dispersion of incidences
similar between the different stratifications.

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
    ## Mixed              0.83900      0.433
    ## Ocean      0.83147              0.868
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
    ## Groups     2 0.004075 0.0020373 0.1745    999  0.859
    ## Residuals 19 0.221851 0.0116764                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              Mixed   Ocean Stratified
    ## Mixed              0.74000      0.398
    ## Ocean      0.71901              0.938
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

We use adonis to determine if dissimilarities of species incidences by
stratification categories are significant and how much of the variation
is explained by incidence dissimilarities.

``` r
### Stratification
## Jaccard
# All sites
jpres_pms <- adonis2(jpres_dist ~ env[-20,19], permutations = 999)
jpres_pms
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = jpres_dist ~ env[-20, 19], permutations = 999)
    ##              Df SumOfSqs      R2      F Pr(>F)    
    ## env[-20, 19]  2   2.1042 0.25523 3.2556  0.001 ***
    ## Residual     19   6.1401 0.74477                  
    ## Total        21   8.2443 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratifcation
jpres_pms_pair <- pairwise.adonis(jpres_dist, env[-20,19], p.adjust.m = "bonferroni", perm = 999)
jpres_pms_pair
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.3049256 3.903468 0.2454476   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.4999441 1.560433 0.1150725   0.047      0.141

``` r
## Dice-Sørensen
# All sites
spres_pms <- adonis2(spres_dist ~ env[-20,19], permutations = 999)
spres_pms
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = spres_dist ~ env[-20, 19], permutations = 999)
    ##              Df SumOfSqs     R2      F Pr(>F)    
    ## env[-20, 19]  2   2.5015 0.3572 5.2792  0.001 ***
    ## Residual     19   4.5016 0.6428                  
    ## Total        21   7.0032 1.0000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratifcation
spres_pms_pair <- pairwise.adonis(spres_dist, env[-20,19], p.adjust.m = "bonferroni", perm = 999)
spres_pms_pair
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.6177889 7.163751 0.3384916   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.6459069 6.509883 0.3516977   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.4361964 1.864338 0.1344700   0.077      0.231

## Envfit environmental influence on incidences

We use envfit to determine significantly correlated environmental
variables to our NMDS. We will use these results, if significant, in our
figure.

``` r
### Environmental
## Jaccard
# All sites for figure
jprese_ef <- envfit(jprese_NMDS, env[-c(9,16,18,20),c(36)], permutations = 1000, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
jprese_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##        NMDS1   NMDS2     r2   Pr(>r)   
    ## [1,] 0.94640 0.32299 0.8332 0.001998 **
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
    ## [1,] 0.94640 0.32299 0.8332 0.001998 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
jpresea_ef <- envfit(jprese_NMDS, env[-c(9,16,18,20),c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
jpresea_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.48791 -0.87289 0.0586 0.853147   
    ## salinity_median     0.99783 -0.06579 0.7388 0.030969 * 
    ## oxygen_median       0.79316  0.60902 0.5361 0.825175   
    ## pH_median           0.94640  0.32299 0.8332 0.001998 **
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
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.48791 -0.87289 0.0586 1.000000   
    ## salinity_median     0.99783 -0.06579 0.7388 0.123876   
    ## oxygen_median       0.79316  0.60902 0.5361 1.000000   
    ## pH_median           0.94640  0.32299 0.8332 0.007992 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
jpresems_ef <- envfit(jpresems_NMDS, env[c(1:6,8,10,12:15,17,21:23),c(10)], permutations = 1000, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
jpresems_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##        NMDS1   NMDS2     r2   Pr(>r)   
    ## [1,] 0.96333 0.26831 0.8029 0.001998 **
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
    ##        NMDS1   NMDS2     r2   Pr(>r)   
    ## [1,] 0.96333 0.26831 0.8029 0.001998 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
jpreseom_ef <- envfit(jpreseom_NMDS, env[c(3,6:8,10:11,13:15,19,23),c(6)], permutations = 1000, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
jpreseom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2   Pr(>r)   
    ## [1,]  0.78731 -0.61655 0.7842 0.008991 **
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
    ##         NMDS1    NMDS2     r2   Pr(>r)   
    ## [1,]  0.78731 -0.61655 0.7842 0.008991 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
jpreseso_ef <- envfit(jpreseso_NMDS, env[c(1:2,4,5,7,11:12,17,19,21:22),c(8)], permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,11:12,17,19,21:22),19])
jpreseso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##          NMDS1     NMDS2     r2 Pr(>r)
    ## [1,] 0.0044465 0.9999900 0.7122 0.8871
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
    ##          NMDS1     NMDS2     r2 Pr(>r)
    ## [1,] 0.0044465 0.9999900 0.7122 0.8871
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
## Dice-Sørensen
# All sites 
sprese_ef <- envfit(sprese_NMDS, env[-c(9,16,18,20),c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
sprese_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.49080 -0.87127 0.0577 0.843157   
    ## salinity_median     0.99690 -0.07862 0.7388 0.022977 * 
    ## oxygen_median       0.76133  0.64837 0.5413 0.813187   
    ## pH_median           0.95591  0.29365 0.8277 0.002997 **
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
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.49080 -0.87127 0.0577 1.00000  
    ## salinity_median     0.99690 -0.07862 0.7388 0.09191 .
    ## oxygen_median       0.76133  0.64837 0.5413 1.00000  
    ## pH_median           0.95591  0.29365 0.8277 0.01199 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
spresems_ef <- envfit(spresems_NMDS, env[c(1:6,8,10,12:15,17,21:23),c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
spresems_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2   Pr(>r)   
    ## temperature_median -0.56687 -0.82381 0.0777 0.796204   
    ## salinity_median     0.99178  0.12798 0.7458 0.040959 * 
    ## oxygen_median       0.94372  0.33076 0.3551 0.888112   
    ## pH_median           0.96347  0.26783 0.8029 0.004995 **
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
    ## temperature_median -0.56687 -0.82381 0.0777 1.00000  
    ## salinity_median     0.99178  0.12798 0.7458 0.16384  
    ## oxygen_median       0.94372  0.33076 0.3551 1.00000  
    ## pH_median           0.96347  0.26783 0.8029 0.01998 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
spreseom_ef <- envfit(spreseom_NMDS, env[c(3,6:8,10:11,13:15,19,23),c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
spreseom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median  0.10614  0.99435 0.0704 0.70330  
    ## salinity_median     0.78723 -0.61666 0.7842 0.01698 *
    ## oxygen_median       0.96714  0.25424 0.7188 0.03896 *
    ## pH_median           0.81113  0.58486 0.6421 0.02098 *
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
    ## temperature_median  0.10614  0.99435 0.0704 1.00000  
    ## salinity_median     0.78723 -0.61666 0.7842 0.06793 .
    ## oxygen_median       0.96714  0.25424 0.7188 0.15584  
    ## pH_median           0.81113  0.58486 0.6421 0.08392 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
spreseso_ef <- envfit(spreseso_NMDS, env[c(1:2,4,5,7,11:12,17,19,21:22),c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,11:12,17,19,21:22),19])
spreseso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median -0.0002470 -1.0000000 0.1263 0.3756
    ## salinity_median     0.0078102 -0.9999700 0.4841 0.8731
    ## oxygen_median       0.0031896  0.9999900 0.7249 0.5684
    ## pH_median           0.0239457  0.9997100 0.6686 0.9071
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
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median -0.0002470 -1.0000000 0.1263      1
    ## salinity_median     0.0078102 -0.9999700 0.4841      1
    ## oxygen_median       0.0031896  0.9999900 0.7249      1
    ## pH_median           0.0239457  0.9997100 0.6686      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
jpress_ef <- envfit(jpress_NMDS, env[c(1:2,4,5,12,17,21:22),c(2,6,8,10,22:32)], permutations = 1000, na.rm = TRUE)
jpress_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median         -0.35881 -0.93341 0.1964 0.6334
    ## salinity_median            -0.70036  0.71379 0.5991 0.1279
    ## oxygen_median               0.98543 -0.17006 0.5185 0.1668
    ## pH_median                  -0.98992 -0.14163 0.5058 0.1459
    ## volume_m3_w_chemocline     -0.99942 -0.03420 0.1624 0.6823
    ## volume_m3                  -0.44002 -0.89799 0.1263 0.6434
    ## surface_area_m2            -0.35025 -0.93666 0.2598 0.4386
    ## distance_to_ocean_min_m     0.63373 -0.77355 0.2151 0.5125
    ## distance_to_ocean_mean_m    0.01326 -0.99991 0.1726 0.6014
    ## distance_to_ocean_median_m -0.00456 -0.99999 0.2084 0.5385
    ## tidal_lag_time_minutes     -0.95579 -0.29405 0.3022 0.3906
    ## tidal_efficiency            0.35277  0.93571 0.0640 0.8252
    ## perimeter_fromSat          -0.28273 -0.95920 0.3722 0.3037
    ## max_depth                  -0.46293 -0.88640 0.0126 0.9620
    ## logArea                    -0.26513 -0.96421 0.2144 0.5245
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
    ## volume_m3_w_chemocline     -0.99942 -0.03420 0.1624      1
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
### Biogeographic
## Jaccard
# All sites for figure
jpresb_ef <- envfit(jpresb_NMDS, env[-c(8,20),c(47)], permutations = 1000, na.rm = TRUE, strata = env[-c(8,20),19])
jpresb_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## [1,]  0.20932 -0.97785 0.3086 0.06194 .
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
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## [1,]  0.20932 -0.97785 0.3086 0.06194 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
jpresba_ef <- envfit(jpresb_NMDS, env[-c(8,20),c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[-c(8,20),19])
jpresba_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.68420 -0.72930 0.1623 0.42058  
    ## volume_m3                   0.64631 -0.76308 0.1557 0.38162  
    ## surface_area_m2             0.43194 -0.90190 0.2799 0.18182  
    ## distance_to_ocean_min_m    -0.98169 -0.19050 0.6466 0.11688  
    ## distance_to_ocean_mean_m   -0.95198 -0.30617 0.5726 0.53447  
    ## distance_to_ocean_median_m -0.94460 -0.32822 0.5670 0.51748  
    ## tidal_lag_time_minutes     -0.99902 -0.04415 0.5908 0.70230  
    ## tidal_efficiency            0.99868  0.05130 0.6114 0.41259  
    ## perimeter_fromSat           0.53004 -0.84797 0.2175 0.38162  
    ## max_depth                  -0.22623 -0.97407 0.1435 0.75924  
    ## logArea                     0.20932 -0.97785 0.3086 0.04496 *
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
    ## volume_m3_w_chemocline      0.68420 -0.72930 0.1623 1.0000
    ## volume_m3                   0.64631 -0.76308 0.1557 1.0000
    ## surface_area_m2             0.43194 -0.90190 0.2799 1.0000
    ## distance_to_ocean_min_m    -0.98169 -0.19050 0.6466 1.0000
    ## distance_to_ocean_mean_m   -0.95198 -0.30617 0.5726 1.0000
    ## distance_to_ocean_median_m -0.94460 -0.32822 0.5670 1.0000
    ## tidal_lag_time_minutes     -0.99902 -0.04415 0.5908 1.0000
    ## tidal_efficiency            0.99868  0.05130 0.6114 1.0000
    ## perimeter_fromSat           0.53004 -0.84797 0.2175 1.0000
    ## max_depth                  -0.22623 -0.97407 0.1435 1.0000
    ## logArea                     0.20932 -0.97785 0.3086 0.4945
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
jpresbms_ef <- envfit(jpresbms_NMDS, env[c(1:6,10,12:15,17,21:23),c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[c(1:6,10,12:15,17,21:23),19])
jpresbms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.35181  0.93607 0.0099 0.9830
    ## volume_m3                  -0.68236 -0.73101 0.1272 0.8981
    ## surface_area_m2            -0.93217  0.36201 0.0832 0.9241
    ## distance_to_ocean_min_m    -0.94998  0.31230 0.5323 0.1449
    ## distance_to_ocean_mean_m   -0.85542  0.51793 0.4286 0.4915
    ## distance_to_ocean_median_m -0.81025  0.58609 0.4367 0.4505
    ## tidal_lag_time_minutes     -0.74305  0.66923 0.4175 0.6693
    ## tidal_efficiency            0.58875 -0.80832 0.5328 0.1808
    ## perimeter_fromSat          -0.25657  0.96653 0.0917 0.7353
    ## max_depth                  -0.35004 -0.93673 0.1947 0.7363
    ## logArea                    -0.88135  0.47246 0.0201 0.9740
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
    ## volume_m3_w_chemocline     -0.35181  0.93607 0.0099      1
    ## volume_m3                  -0.68236 -0.73101 0.1272      1
    ## surface_area_m2            -0.93217  0.36201 0.0832      1
    ## distance_to_ocean_min_m    -0.94998  0.31230 0.5323      1
    ## distance_to_ocean_mean_m   -0.85542  0.51793 0.4286      1
    ## distance_to_ocean_median_m -0.81025  0.58609 0.4367      1
    ## tidal_lag_time_minutes     -0.74305  0.66923 0.4175      1
    ## tidal_efficiency            0.58875 -0.80832 0.5328      1
    ## perimeter_fromSat          -0.25657  0.96653 0.0917      1
    ## max_depth                  -0.35004 -0.93673 0.1947      1
    ## logArea                    -0.88135  0.47246 0.0201      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
jpresbom_ef <- envfit(jpresbom_NMDS, env[c(3,6:7,9:11,13:16,18:19,23),c(25)], permutations = 1000, na.rm = TRUE, strata = env[c(3,6:7,9:11,13:16,18:19,23),19])
jpresbom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## [1,] -0.84334  0.53738 0.4301 0.02597 *
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
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## [1,] -0.84334  0.53738 0.4301 0.02597 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
jpresbso_ef <- envfit(jpresbso_NMDS, env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
jpresbso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.94205  0.33547 0.2252 0.7193
    ## volume_m3                   0.95257  0.30431 0.2135 0.7333
    ## surface_area_m2             0.78573 -0.61857 0.3178 0.6613
    ## distance_to_ocean_min_m    -0.89001 -0.45595 0.6701 0.1948
    ## distance_to_ocean_mean_m   -0.89206 -0.45191 0.6471 0.5345
    ## distance_to_ocean_median_m -0.87780 -0.47903 0.6385 0.5315
    ## tidal_lag_time_minutes     -0.98367 -0.18000 0.7191 0.9351
    ## tidal_efficiency            0.97293  0.23112 0.7970 0.5894
    ## perimeter_fromSat           0.99078 -0.13550 0.2877 0.7043
    ## max_depth                  -0.96972 -0.24423 0.0579 0.9700
    ## logArea                     0.54157 -0.84065 0.2156 0.5235
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
    ## volume_m3_w_chemocline      0.94205  0.33547 0.2252      1
    ## volume_m3                   0.95257  0.30431 0.2135      1
    ## surface_area_m2             0.78573 -0.61857 0.3178      1
    ## distance_to_ocean_min_m    -0.89001 -0.45595 0.6701      1
    ## distance_to_ocean_mean_m   -0.89206 -0.45191 0.6471      1
    ## distance_to_ocean_median_m -0.87780 -0.47903 0.6385      1
    ## tidal_lag_time_minutes     -0.98367 -0.18000 0.7191      1
    ## tidal_efficiency            0.97293  0.23112 0.7970      1
    ## perimeter_fromSat           0.99078 -0.13550 0.2877      1
    ## max_depth                  -0.96972 -0.24423 0.0579      1
    ## logArea                     0.54157 -0.84065 0.2156      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
## Dice-Sørensen
# All sites 
spresb_ef <- envfit(spresb_NMDS, env[-c(8,20),c(22,24,26,29,31)], permutations = 1000, na.rm = TRUE, strata = env[-c(8,20),19])
spresb_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline    0.68429 -0.72921 0.1622 0.4426
    ## surface_area_m2           0.43196 -0.90189 0.2799 0.1958
    ## distance_to_ocean_mean_m -0.95195 -0.30626 0.5726 0.5265
    ## tidal_efficiency          0.99868  0.05141 0.6114 0.3886
    ## max_depth                -0.22623 -0.97407 0.1435 0.7612
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
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline    0.68429 -0.72921 0.1622  1.000
    ## surface_area_m2           0.43196 -0.90189 0.2799  0.979
    ## distance_to_ocean_mean_m -0.95195 -0.30626 0.5726  1.000
    ## tidal_efficiency          0.99868  0.05141 0.6114  1.000
    ## max_depth                -0.22623 -0.97407 0.1435  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
spresbms_ef <- envfit(spresbms_NMDS, env[c(1:6,10,12:15,17,21:23),c(22,24,26,29,31)], permutations = 1000, na.rm = TRUE, strata = env[c(1:6,10,12:15,17,21:23),19])
spresbms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline   -0.35176  0.93609 0.0099 0.9830
    ## surface_area_m2          -0.93215  0.36208 0.0832 0.9341
    ## distance_to_ocean_mean_m -0.85542  0.51793 0.4286 0.4875
    ## tidal_efficiency          0.58875 -0.80831 0.5328 0.2078
    ## max_depth                -0.35005 -0.93673 0.1947 0.7453
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
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline   -0.35176  0.93609 0.0099      1
    ## surface_area_m2          -0.93215  0.36208 0.0832      1
    ## distance_to_ocean_mean_m -0.85542  0.51793 0.4286      1
    ## tidal_efficiency          0.58875 -0.80831 0.5328      1
    ## max_depth                -0.35005 -0.93673 0.1947      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
spresbom_ef <- envfit(spresbom_NMDS, env[c(3,6:7,9:11,13:16,18:19,23),c(22,24,26,29,31)], permutations = 1000, na.rm = TRUE, strata = env[c(3,6:7,9:11,13:16,18:19,23),19])
spresbom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                             NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline    0.94026  0.34047 0.1657 0.57143  
    ## surface_area_m2           0.73074 -0.68266 0.2650 0.54246  
    ## distance_to_ocean_mean_m -0.72885  0.68467 0.2299 0.35564  
    ## tidal_efficiency          0.76001 -0.64991 0.2206 0.28472  
    ## max_depth                 0.66523  0.74664 0.5244 0.03497 *
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
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline    0.94026  0.34047 0.1657 1.0000
    ## surface_area_m2           0.73074 -0.68266 0.2650 1.0000
    ## distance_to_ocean_mean_m -0.72885  0.68467 0.2299 1.0000
    ## tidal_efficiency          0.76001 -0.64991 0.2206 1.0000
    ## max_depth                 0.66523  0.74664 0.5244 0.1748
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
spresbso_ef <- envfit(spresbso_NMDS, env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(22,24,26,29,31)], permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
spresbso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline    0.94185  0.33603 0.2252 0.7443
    ## surface_area_m2           0.78585 -0.61842 0.3178 0.6503
    ## distance_to_ocean_mean_m -0.89210 -0.45184 0.6471 0.5534
    ## tidal_efficiency          0.97290  0.23121 0.7970 0.6014
    ## max_depth                -0.96997 -0.24321 0.0579 0.9760
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
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline    0.94185  0.33603 0.2252      1
    ## surface_area_m2           0.78585 -0.61842 0.3178      1
    ## distance_to_ocean_mean_m -0.89210 -0.45184 0.6471      1
    ## tidal_efficiency          0.97290  0.23121 0.7970      1
    ## max_depth                -0.96997 -0.24321 0.0579      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Traits
# All sites for figure 
jprest_ef <- envfit(jpres_NMDS, FD_total_env[-20,c(57:58,64)], permutations = 1000, na.rm = TRUE, strata = env[-20,19])
jprest_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##      NMDS1    NMDS2     r2   Pr(>r)   
    ## L  0.13029  0.99148 0.3305 0.005994 **
    ## T -0.97547 -0.22011 0.7008 0.010989 * 
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
    ## EC1b -1.0503 -0.5277
    ## EC2g -1.5644  0.0244
    ## EC3n  0.2905  0.0559
    ## 
    ## Goodness of fit:
    ##        r2  Pr(>r)  
    ## EC 0.4777 0.01099 *
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
    ##      NMDS1    NMDS2     r2  Pr(>r)  
    ## L  0.13029  0.99148 0.3305 0.01798 *
    ## T -0.97547 -0.22011 0.7008 0.03297 *
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
    ## EC1b -1.0503 -0.5277
    ## EC2g -1.5644  0.0244
    ## EC3n  0.2905  0.0559
    ## 
    ## Goodness of fit:
    ##        r2  Pr(>r)  
    ## EC 0.4777 0.03297 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification 9:10,18:19
jpresta_ef <- envfit(jpres_NMDS, FD_total_env[-20,c(6:22)], permutations = 1000, na.rm = TRUE, strata = env[-20,19])
jpresta_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.13029  0.99148 0.3305 0.001998 **
    ## Troph            -0.97547 -0.22011 0.7008 0.013986 * 
    ## DepthMin         -0.47476  0.88012 0.0104 0.899101   
    ## DepthMax          0.48867  0.87247 0.3001 0.064935 . 
    ## TempPrefMin       0.53929 -0.84212 0.3079 0.143856   
    ## TempPrefMax      -0.06392 -0.99795 0.0801 0.160839   
    ## Weight            0.01274  0.99992 0.1122 0.114885   
    ## CaudalFinLength  -0.13219  0.99123 0.1757 0.014985 * 
    ## DorsalSpinesMean  0.99278 -0.11997 0.6012 0.017982 * 
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
    ## BodyShapeI4e        -1.2708 -0.2090
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.1634 -0.3024
    ## OperculumPresentyes  0.2585  0.0672
    ## FeedingPathb         0.1050  0.0528
    ## FeedingPathp        -1.0503 -0.5277
    ## RepGuild11b         -1.0503 -0.5277
    ## RepGuild12g         -1.5644  0.0244
    ## RepGuild13n          0.2905  0.0559
    ## RepGuild22eb        -1.1059 -0.4106
    ## RepGuild23n         -1.5644  0.0244
    ## RepGuild26s          0.2229  0.0190
    ## ParentalCare3p      -1.1202 -0.1985
    ## ParentalCare4n       0.4201  0.0744
    ## WaterPref1s          0.3596  0.0225
    ## WaterPref3a         -1.2225 -0.0764
    ## 
    ## Goodness of fit:
    ##                      r2   Pr(>r)   
    ## BodyShapeI       0.4385 0.041958 * 
    ## DemersPelag      0.0000 1.000000   
    ## OperculumPresent 0.3654 0.097902 . 
    ## FeedingPath      0.1572 0.340659   
    ## RepGuild1        0.4777 0.009990 **
    ## RepGuild2        0.3744 0.006993 **
    ## ParentalCare     0.5523 0.081918 . 
    ## WaterPref        0.5021 0.055944 . 
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
    ## MaxLengthTL       0.13029  0.99148 0.3305 0.03397 *
    ## Troph            -0.97547 -0.22011 0.7008 0.23776  
    ## DepthMin         -0.47476  0.88012 0.0104 1.00000  
    ## DepthMax          0.48867  0.87247 0.3001 1.00000  
    ## TempPrefMin       0.53929 -0.84212 0.3079 1.00000  
    ## TempPrefMax      -0.06392 -0.99795 0.0801 1.00000  
    ## Weight            0.01274  0.99992 0.1122 1.00000  
    ## CaudalFinLength  -0.13219  0.99123 0.1757 0.25475  
    ## DorsalSpinesMean  0.99278 -0.11997 0.6012 0.30569  
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
    ## BodyShapeI4e        -1.2708 -0.2090
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.1634 -0.3024
    ## OperculumPresentyes  0.2585  0.0672
    ## FeedingPathb         0.1050  0.0528
    ## FeedingPathp        -1.0503 -0.5277
    ## RepGuild11b         -1.0503 -0.5277
    ## RepGuild12g         -1.5644  0.0244
    ## RepGuild13n          0.2905  0.0559
    ## RepGuild22eb        -1.1059 -0.4106
    ## RepGuild23n         -1.5644  0.0244
    ## RepGuild26s          0.2229  0.0190
    ## ParentalCare3p      -1.1202 -0.1985
    ## ParentalCare4n       0.4201  0.0744
    ## WaterPref1s          0.3596  0.0225
    ## WaterPref3a         -1.2225 -0.0764
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.4385 0.7133
    ## DemersPelag      0.0000 1.0000
    ## OperculumPresent 0.3654 1.0000
    ## FeedingPath      0.1572 1.0000
    ## RepGuild1        0.4777 0.1698
    ## RepGuild2        0.3744 0.1189
    ## ParentalCare     0.5523 1.0000
    ## WaterPref        0.5021 0.9510
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes 12:14,18
jprestms_ef <- envfit(jpresms_NMDS, FD_total_env[c(1:6,8,10,12:15,17,21:23),c(12:14,18)], permutations = 1000, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
jprestms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                NMDS1    NMDS2     r2   Pr(>r)    
    ## DepthMax     0.37915  0.92534 0.5650 0.014985 *  
    ## TempPrefMin  0.13433 -0.99094 0.7048 0.000999 ***
    ## TempPrefMax -0.04422 -0.99902 0.4378 0.005994 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##               NMDS1   NMDS2
    ## RepGuild11b -0.8545 -0.3481
    ## RepGuild12g -1.1597  0.1376
    ## RepGuild13n  0.3357  0.0351
    ## 
    ## Goodness of fit:
    ##               r2   Pr(>r)   
    ## RepGuild1 0.5641 0.002997 **
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
    ##                NMDS1    NMDS2     r2   Pr(>r)   
    ## DepthMax     0.37915  0.92534 0.5650 0.059940 . 
    ## TempPrefMin  0.13433 -0.99094 0.7048 0.003996 **
    ## TempPrefMax -0.04422 -0.99902 0.4378 0.023976 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##               NMDS1   NMDS2
    ## RepGuild11b -0.8545 -0.3481
    ## RepGuild12g -1.1597  0.1376
    ## RepGuild13n  0.3357  0.0351
    ## 
    ## Goodness of fit:
    ##               r2  Pr(>r)  
    ## RepGuild1 0.5641 0.01199 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes 9:10
jprestom_ef <- envfit(jpresom_NMDS, FD_total_env[c(3,6:11,13:16,18:19,23),c(9:10)], permutations = 1000, na.rm = TRUE, strata = env[c(3,6:11,13:16,18:19,23),19])
jprestom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL -0.99998 -0.00681 0.4907 0.011988 * 
    ## Troph       -0.66850 -0.74372 0.6816 0.002997 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ##                NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL -0.99998 -0.00681 0.4907 0.023976 * 
    ## Troph       -0.66850 -0.74372 0.6816 0.005994 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites 9:10,18:19
jprestso_ef <- envfit(jpresso_NMDS, FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(9:10,18:19)], permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
jprestso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                NMDS1    NMDS2     r2   Pr(>r)    
    ## MaxLengthTL  0.10154  0.99483 0.3869 0.000999 ***
    ## Troph       -0.79702 -0.60395 0.8666 0.009990 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                NMDS1   NMDS2
    ## RepGuild11b  -0.6890 -0.3566
    ## RepGuild12g  -1.0433  0.0617
    ## RepGuild13n   0.3465  0.0590
    ## RepGuild22eb -0.7282 -0.2733
    ## RepGuild23n  -1.0433  0.0617
    ## RepGuild26s   0.2559  0.0136
    ## 
    ## Goodness of fit:
    ##               r2   Pr(>r)   
    ## RepGuild1 0.4562 0.004995 **
    ## RepGuild2 0.3465 0.010989 * 
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
    ##                NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL  0.10154  0.99483 0.3869 0.003996 **
    ## Troph       -0.79702 -0.60395 0.8666 0.039960 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                NMDS1   NMDS2
    ## RepGuild11b  -0.6890 -0.3566
    ## RepGuild12g  -1.0433  0.0617
    ## RepGuild13n   0.3465  0.0590
    ## RepGuild22eb -0.7282 -0.2733
    ## RepGuild23n  -1.0433  0.0617
    ## RepGuild26s   0.2559  0.0136
    ## 
    ## Goodness of fit:
    ##               r2  Pr(>r)  
    ## RepGuild1 0.4562 0.01998 *
    ## RepGuild2 0.3465 0.04396 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Traits
# All sites 
sprest_ef <- envfit(spres_NMDS, FD_total_env[-20,c(6:22)], permutations = 1000, na.rm = TRUE, strata = env[-20,19])
sprest_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.12761  0.99182 0.3246 0.003996 **
    ## Troph            -0.98932 -0.14574 0.7081 0.005994 **
    ## DepthMin         -0.16081  0.98699 0.0493 0.337662   
    ## DepthMax          0.42707  0.90422 0.3298 0.025974 * 
    ## TempPrefMin       0.37474 -0.92713 0.4123 0.018981 * 
    ## TempPrefMax      -0.04329 -0.99906 0.1597 0.034965 * 
    ## Weight            0.02500  0.99969 0.0790 0.209790   
    ## CaudalFinLength  -0.11136  0.99378 0.2060 0.012987 * 
    ## DorsalSpinesMean  0.91746 -0.39783 0.6344 0.006993 **
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
    ## BodyShapeI2s         0.3322 -0.3011
    ## BodyShapeI3f         0.3370  0.0848
    ## BodyShapeI4e        -1.5131 -0.0921
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.3882 -0.2477
    ## OperculumPresentyes  0.3085  0.0550
    ## FeedingPathb         0.1267  0.0587
    ## FeedingPathp        -1.2673 -0.5874
    ## RepGuild11b         -1.2673 -0.5874
    ## RepGuild12g         -1.8314  0.0235
    ## RepGuild13n          0.3443  0.0627
    ## RepGuild22eb        -1.3282 -0.4187
    ## RepGuild23n         -1.8314  0.0235
    ## RepGuild26s          0.2627  0.0196
    ## ParentalCare3p      -1.3325 -0.1065
    ## ParentalCare4n       0.4997  0.0399
    ## WaterPref1s          0.4232  0.0126
    ## WaterPref3a         -1.4390 -0.0429
    ## 
    ## Goodness of fit:
    ##                      r2   Pr(>r)   
    ## BodyShapeI       0.4394 0.046953 * 
    ## DemersPelag      0.0000 1.000000   
    ## OperculumPresent 0.3680 0.109890   
    ## FeedingPath      0.1625 0.358641   
    ## RepGuild1        0.4851 0.003996 **
    ## RepGuild2        0.3773 0.008991 **
    ## ParentalCare     0.5580 0.081918 . 
    ## WaterPref        0.5077 0.051948 . 
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
    ## MaxLengthTL       0.12761  0.99182 0.3246 0.06793 .
    ## Troph            -0.98932 -0.14574 0.7081 0.10190  
    ## DepthMin         -0.16081  0.98699 0.0493 1.00000  
    ## DepthMax          0.42707  0.90422 0.3298 0.44156  
    ## TempPrefMin       0.37474 -0.92713 0.4123 0.32268  
    ## TempPrefMax      -0.04329 -0.99906 0.1597 0.59441  
    ## Weight            0.02500  0.99969 0.0790 1.00000  
    ## CaudalFinLength  -0.11136  0.99378 0.2060 0.22078  
    ## DorsalSpinesMean  0.91746 -0.39783 0.6344 0.11888  
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
    ## BodyShapeI2s         0.3322 -0.3011
    ## BodyShapeI3f         0.3370  0.0848
    ## BodyShapeI4e        -1.5131 -0.0921
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.3882 -0.2477
    ## OperculumPresentyes  0.3085  0.0550
    ## FeedingPathb         0.1267  0.0587
    ## FeedingPathp        -1.2673 -0.5874
    ## RepGuild11b         -1.2673 -0.5874
    ## RepGuild12g         -1.8314  0.0235
    ## RepGuild13n          0.3443  0.0627
    ## RepGuild22eb        -1.3282 -0.4187
    ## RepGuild23n         -1.8314  0.0235
    ## RepGuild26s          0.2627  0.0196
    ## ParentalCare3p      -1.3325 -0.1065
    ## ParentalCare4n       0.4997  0.0399
    ## WaterPref1s          0.4232  0.0126
    ## WaterPref3a         -1.4390 -0.0429
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.4394 0.79820  
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.3680 1.00000  
    ## FeedingPath      0.1625 1.00000  
    ## RepGuild1        0.4851 0.06793 .
    ## RepGuild2        0.3773 0.15285  
    ## ParentalCare     0.5580 1.00000  
    ## WaterPref        0.5077 0.88312  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
sprestms_ef <- envfit(spresms_NMDS, FD_total_env[c(1:6,8,10,12:15,17,21:23),c(6:22)], permutations = 1000, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
sprestms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)    
    ## MaxLengthTL       0.23258  0.97258 0.4679 0.015984 *  
    ## Troph            -0.99771  0.06768 0.5767 0.057942 .  
    ## DepthMin         -0.11498  0.99337 0.1214 0.274725    
    ## DepthMax          0.37897  0.92541 0.5653 0.007992 ** 
    ## TempPrefMin       0.13434 -0.99094 0.7048 0.000999 ***
    ## TempPrefMax      -0.04422 -0.99902 0.4378 0.011988 *  
    ## Weight            0.13643  0.99065 0.3221 0.036963 *  
    ## CaudalFinLength   0.03300  0.99946 0.3027 0.040959 *  
    ## DorsalSpinesMean  0.98054 -0.19630 0.6172 0.015984 *  
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
    ## BodyShapeI2s        -1.0366 -0.3240
    ## BodyShapeI3f         0.5262  0.0115
    ## BodyShapeI4e        -1.1879  0.0495
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.1084 -0.1168
    ## OperculumPresentyes  0.3695  0.0389
    ## FeedingPathb         0.1496  0.0609
    ## FeedingPathp        -1.0474 -0.4266
    ## RepGuild11b         -1.0474 -0.4266
    ## RepGuild12g         -1.4216  0.1683
    ## RepGuild13n          0.4115  0.0431
    ## RepGuild22eb        -1.0366 -0.3240
    ## RepGuild23n         -1.4216  0.1683
    ## RepGuild26s          0.2984 -0.0010
    ## ParentalCare3p      -1.0069  0.0348
    ## ParentalCare4n       0.6041 -0.0209
    ## WaterPref1s          0.4828 -0.0099
    ## WaterPref3a         -1.0622  0.0219
    ## 
    ## Goodness of fit:
    ##                      r2   Pr(>r)   
    ## BodyShapeI       0.6398 0.016983 * 
    ## DemersPelag      0.0000 1.000000   
    ## OperculumPresent 0.4290 0.078921 . 
    ## FeedingPath      0.1893 0.283716   
    ## RepGuild1        0.5640 0.004995 **
    ## RepGuild2        0.4168 0.019980 * 
    ## ParentalCare     0.6310 0.079920 . 
    ## WaterPref        0.5315 0.061938 . 
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
    ## MaxLengthTL       0.23258  0.97258 0.4679 0.27173  
    ## Troph            -0.99771  0.06768 0.5767 0.98501  
    ## DepthMin         -0.11498  0.99337 0.1214 1.00000  
    ## DepthMax          0.37897  0.92541 0.5653 0.13586  
    ## TempPrefMin       0.13434 -0.99094 0.7048 0.01698 *
    ## TempPrefMax      -0.04422 -0.99902 0.4378 0.20380  
    ## Weight            0.13643  0.99065 0.3221 0.62837  
    ## CaudalFinLength   0.03300  0.99946 0.3027 0.69630  
    ## DorsalSpinesMean  0.98054 -0.19630 0.6172 0.27173  
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
    ## BodyShapeI2s        -1.0366 -0.3240
    ## BodyShapeI3f         0.5262  0.0115
    ## BodyShapeI4e        -1.1879  0.0495
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.1084 -0.1168
    ## OperculumPresentyes  0.3695  0.0389
    ## FeedingPathb         0.1496  0.0609
    ## FeedingPathp        -1.0474 -0.4266
    ## RepGuild11b         -1.0474 -0.4266
    ## RepGuild12g         -1.4216  0.1683
    ## RepGuild13n          0.4115  0.0431
    ## RepGuild22eb        -1.0366 -0.3240
    ## RepGuild23n         -1.4216  0.1683
    ## RepGuild26s          0.2984 -0.0010
    ## ParentalCare3p      -1.0069  0.0348
    ## ParentalCare4n       0.6041 -0.0209
    ## WaterPref1s          0.4828 -0.0099
    ## WaterPref3a         -1.0622  0.0219
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.6398 0.28871  
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.4290 1.00000  
    ## FeedingPath      0.1893 1.00000  
    ## RepGuild1        0.5640 0.08492 .
    ## RepGuild2        0.4168 0.33966  
    ## ParentalCare     0.6310 1.00000  
    ## WaterPref        0.5315 1.00000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
sprestom_ef <- envfit(spresom_NMDS, FD_total_env[c(3,6:11,13:16,18:19,23),c(6:22)], permutations = 1000, na.rm = TRUE, strata = env[c(3,6:11,13:16,18:19,23),19])
sprestom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL      -0.99998 -0.00669 0.4907 0.011988 * 
    ## Troph            -0.66853 -0.74368 0.6816 0.003996 **
    ## DepthMin          0.29446  0.95566 0.1726 0.219780   
    ## DepthMax         -0.10248 -0.99473 0.0077 0.977023   
    ## TempPrefMin       0.63816 -0.76991 0.2869 0.101898   
    ## TempPrefMax       0.34836 -0.93736 0.0734 0.528472   
    ## Weight           -0.99712  0.07580 0.0534 0.661339   
    ## CaudalFinLength  -0.99876  0.04983 0.3784 0.033966 * 
    ## DorsalSpinesMean  0.78594 -0.61830 0.3655 0.065934 . 
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
    ## BodyShapeI2s         0.5434  0.0555
    ## BodyShapeI3f        -0.0906 -0.0092
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
    ## BodyShapeI       0.117 0.6733
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
    ## MaxLengthTL      -0.99998 -0.00669 0.4907 0.20380  
    ## Troph            -0.66853 -0.74368 0.6816 0.06793 .
    ## DepthMin          0.29446  0.95566 0.1726 1.00000  
    ## DepthMax         -0.10248 -0.99473 0.0077 1.00000  
    ## TempPrefMin       0.63816 -0.76991 0.2869 1.00000  
    ## TempPrefMax       0.34836 -0.93736 0.0734 1.00000  
    ## Weight           -0.99712  0.07580 0.0534 1.00000  
    ## CaudalFinLength  -0.99876  0.04983 0.3784 0.57742  
    ## DorsalSpinesMean  0.78594 -0.61830 0.3655 1.00000  
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
    ## BodyShapeI2s         0.5434  0.0555
    ## BodyShapeI3f        -0.0906 -0.0092
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
sprestso_ef <- envfit(spresso_NMDS, FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(6:22)], permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
sprestso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.10156  0.99483 0.3868 0.001998 **
    ## Troph            -0.79698 -0.60400 0.8666 0.009990 **
    ## DepthMin         -0.20026  0.97974 0.0314 0.693307   
    ## DepthMax          0.39587  0.91831 0.2840 0.040959 * 
    ## TempPrefMin       0.39054 -0.92059 0.3671 0.123876   
    ## TempPrefMax      -0.08494 -0.99639 0.1298 0.128871   
    ## Weight           -0.05007  0.99875 0.2028 0.116883   
    ## CaudalFinLength  -0.11160  0.99375 0.2116 0.052947 . 
    ## DorsalSpinesMean  0.99754  0.07007 0.6372 0.058941 . 
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
    ## BodyShapeI2s         0.5873 -0.2448
    ## BodyShapeI3f         0.4117  0.2046
    ## BodyShapeI4e        -1.1609 -0.1744
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.0573 -0.3140
    ## OperculumPresentyes  0.4229  0.1256
    ## FeedingPathb         0.1626  0.0841
    ## FeedingPathp        -0.9753 -0.5047
    ## RepGuild11b         -0.9753 -0.5047
    ## RepGuild12g         -1.4767  0.0876
    ## RepGuild13n          0.4904  0.0834
    ## RepGuild22eb        -1.0308 -0.3868
    ## RepGuild23n         -1.4767  0.0876
    ## RepGuild26s          0.3622  0.0192
    ## ParentalCare3p      -1.0120 -0.1455
    ## ParentalCare4n       0.7590  0.1091
    ## WaterPref1s          0.6224  0.0027
    ## WaterPref3a         -1.1203 -0.0049
    ## 
    ## Goodness of fit:
    ##                      r2   Pr(>r)   
    ## BodyShapeI       0.4042 0.056943 . 
    ## DemersPelag      0.0000 1.000000   
    ## OperculumPresent 0.3356 0.091908 . 
    ## FeedingPath      0.1386 0.369630   
    ## RepGuild1        0.4562 0.003996 **
    ## RepGuild2        0.3465 0.006993 **
    ## ParentalCare     0.5406 0.093906 . 
    ## WaterPref        0.4809 0.032967 * 
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
    ## MaxLengthTL       0.10156  0.99483 0.3868 0.03397 *
    ## Troph            -0.79698 -0.60400 0.8666 0.16983  
    ## DepthMin         -0.20026  0.97974 0.0314 1.00000  
    ## DepthMax          0.39587  0.91831 0.2840 0.69630  
    ## TempPrefMin       0.39054 -0.92059 0.3671 1.00000  
    ## TempPrefMax      -0.08494 -0.99639 0.1298 1.00000  
    ## Weight           -0.05007  0.99875 0.2028 1.00000  
    ## CaudalFinLength  -0.11160  0.99375 0.2116 0.90010  
    ## DorsalSpinesMean  0.99754  0.07007 0.6372 1.00000  
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
    ## BodyShapeI2s         0.5873 -0.2448
    ## BodyShapeI3f         0.4117  0.2046
    ## BodyShapeI4e        -1.1609 -0.1744
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.0573 -0.3140
    ## OperculumPresentyes  0.4229  0.1256
    ## FeedingPathb         0.1626  0.0841
    ## FeedingPathp        -0.9753 -0.5047
    ## RepGuild11b         -0.9753 -0.5047
    ## RepGuild12g         -1.4767  0.0876
    ## RepGuild13n          0.4904  0.0834
    ## RepGuild22eb        -1.0308 -0.3868
    ## RepGuild23n         -1.4767  0.0876
    ## RepGuild26s          0.3622  0.0192
    ## ParentalCare3p      -1.0120 -0.1455
    ## ParentalCare4n       0.7590  0.1091
    ## WaterPref1s          0.6224  0.0027
    ## WaterPref3a         -1.1203 -0.0049
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.4042 0.96803  
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.3356 1.00000  
    ## FeedingPath      0.1386 1.00000  
    ## RepGuild1        0.4562 0.06793 .
    ## RepGuild2        0.3465 0.11888  
    ## ParentalCare     0.5406 1.00000  
    ## WaterPref        0.4809 0.56044  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

## Incidence mantel tests

We used mantel tests to determine significance between the incidence and
environmental distance matrices.

``` r
### Environmental
## Jaccard
# All sites
enve_dist_s <- dist(scaled_env[-c(9,16,18,20),c(2)], method = "euclidean")
jpresea_mant_s <- mantel(jprese_dist, enve_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
jpresea_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jprese_dist, ydis = enve_dist_s, method = "spearman",      permutations = 999, strata = env[-c(9, 16, 18, 20), 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6697 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.560 0.605 0.619 0.642 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_p <- dist(scaled_env[-c(9,16,18,20),c(4)], method = "euclidean")
jpresea_mant_p <- mantel(jprese_dist, enve_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
jpresea_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jprese_dist, ydis = enve_dist_p, method = "spearman",      permutations = 999, strata = env[-c(9, 16, 18, 20), 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7163 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.546 0.579 0.612 0.639 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresea_mant_pv <- rbind(jpresea_mant_s$signif, jpresea_mant_p$signif)
jpresea_mant_pv <- jpresea_mant_pv[,1]
jpresea_mant_pv <- p.adjust(jpresea_mant_pv, method = "bonferroni")
jpresea_mant_pv
```

    ## [1] 0.004 0.002

``` r
# Mixed and stratified lakes
envems_dist_s <- dist(scaled_env[c(1:6,8,10,12:15,17,21:23),c(2)], method = "euclidean")
jpresems_mant_s <- mantel(jpresems_dist, envems_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
jpresems_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresems_dist, ydis = envems_dist_s, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6676 
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.560 0.606 0.631 0.671 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_p <- dist(scaled_env[c(1:6,8,10,12:15,17,21:23),c(4)], method = "euclidean")
jpresems_mant_p <- mantel(jpresems_dist, envems_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
jpresems_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresems_dist, ydis = envems_dist_p, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6357 
    ##       Significance: 0.003 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.420 0.463 0.492 0.542 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresems_mant_pv <- rbind(jpresems_mant_s$signif, jpresems_mant_p$signif)
jpresems_mant_pv <- jpresems_mant_pv[,1]
jpresems_mant_pv <- p.adjust(jpresems_mant_pv, method = "bonferroni")
jpresems_mant_pv
```

    ## [1] 0.024 0.006

``` r
# Ocean sites and mixed lakes
enveom_dist_p <- dist(scaled_env[c(3,6:8,10:11,13:15,19,23),c(4)], method = "euclidean")
jpreseom_mant_p <- mantel(jpreseom_dist, enveom_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
jpreseom_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpreseom_dist, ydis = enveom_dist_p, method = "spearman",      permutations = 999, strata = env[c(3, 6:8, 10:11, 13:15,          19, 23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5793 
    ##       Significance: 0.015 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.298 0.379 0.486 0.617 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
enveso_dist_s <- dist(scaled_env[c(1:2,4,5,7,11:12,17,19,21:22),c(2)], method = "euclidean")
jpreseso_mant_s <- mantel(jpreseso_dist, enveso_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,11:12,17,19,21:22),19])
jpreseso_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpreseso_dist, ydis = enveso_dist_s, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 11:12, 17,          19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5528 
    ##       Significance: 0.016 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.462 0.509 0.539 0.565 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
## Dice-Sørensen
# All sites
enve_dist <- dist(scaled_env[-c(9,16,18,20),c(6,8,10)], method = "euclidean")
sprese_mant <- mantel(sprese_dist, enve_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
sprese_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = sprese_dist, ydis = enve_dist, method = "spearman",      permutations = 999, strata = env[-c(9, 16, 18, 20), 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3836 
    ##       Significance: 0.163 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.417 0.450 0.474 0.503 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
envems_dist <- dist(scaled_env[c(1:6,8,10,12:15,17,21:23),c(6,8,10)], method = "euclidean")
spresems_mant <- mantel(spresems_dist, envems_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
spresems_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresems_dist, ydis = envems_dist, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3463 
    ##       Significance: 0.204 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.409 0.460 0.491 0.528 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
enveom_dist <- dist(scaled_env[c(3,6:8,10:11,13:15,19,23),c(6,8,10)], method = "euclidean")
spreseom_mant <- mantel(spreseom_dist, enveom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
spreseom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spreseom_dist, ydis = enveom_dist, method = "spearman",      permutations = 999, strata = env[c(3, 6:8, 10:11, 13:15,          19, 23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2566 
    ##       Significance: 0.066 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.227 0.274 0.336 0.447 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
enveso_dist <- dist(scaled_env[c(1:2,4,5,7,11:12,17,19,21:22),c(6,8,10)], method = "euclidean")
spreseso_mant <- mantel(spreseso_dist, enveso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,11:12,17,19,21:22),19])
spreseso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spreseso_dist, ydis = enveso_dist, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 11:12, 17,          19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4922 
    ##       Significance: 0.166 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.508 0.535 0.555 0.573 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes
enves_dist <- dist(scaled_env[c(1:2,4,5,12,17,21:22),c(2)], method = "euclidean")
jpreses_mant <- mantel(jpress_dist, enves_dist, method = "spearman", permutations = 999, na.rm = TRUE)
jpreses_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpress_dist, ydis = enves_dist, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4001 
    ##       Significance: 0.021 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.260 0.337 0.371 0.440 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Biogeographic
## Jaccard
# All sites
envb_dist <- dist(scaled_env[-c(8,20),c(5:15)], method = "euclidean")
jpresba_mant <- mantel(jpresb_dist, envb_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(8,20),19])
jpresba_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresb_dist, ydis = envb_dist, method = "spearman",      permutations = 999, strata = env[-c(8, 20), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3881 
    ##       Significance: 0.265 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.440 0.465 0.485 0.509 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
envbms_dist <- dist(scaled_env[c(1:6,10,12:15,17,21:23),c(5:15)], method = "euclidean")
jpresbms_mant <- mantel(jpresbms_dist, envbms_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,10,12:15,17,21:23),19])
jpresbms_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbms_dist, ydis = envbms_dist, method = "spearman",      permutations = 999, strata = env[c(1:6, 10, 12:15, 17, 21:23),          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2892 
    ##       Significance: 0.349 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.412 0.455 0.494 0.535 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
envbom_dist <- dist(scaled_env[c(3,6:7,9:11,13:16,18:19,23),c(5:15)], method = "euclidean")
jpresbom_mant <- mantel(jpresbom_dist, envbom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:7,9:11,13:16,18:19,23),19])
jpresbom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbom_dist, ydis = envbom_dist, method = "spearman",      permutations = 999, strata = env[c(3, 6:7, 9:11, 13:16, 18:19,          23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2261 
    ##       Significance: 0.221 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.305 0.350 0.401 0.457 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
envbso_dist <- dist(scaled_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(5:15)], method = "euclidean")
jpresbso_mant <- mantel(jpresbso_dist, envbso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
jpresbso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresbso_dist, ydis = envbso_dist, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5571 
    ##       Significance: 0.404 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.600 0.618 0.628 0.637 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresb_mant_p <- rbind(jpresbms_mant$signif, jpresbom_mant$signif, jpresbso_mant$signif)
jpresb_mant_p <- jpresb_mant_p[,1]
jpresb_mant_p <- p.adjust(jpresb_mant_p, method = "bonferroni")
jpresb_mant_p
```

    ## [1] 1.000 0.663 1.000

``` r
## Dice-Sørensen
# All sites
spresb_mant <- mantel(spresb_dist, envb_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(8,20),19])
spresb_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresb_dist, ydis = envb_dist, method = "spearman",      permutations = 999, strata = env[-c(8, 20), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3881 
    ##       Significance: 0.285 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.434 0.464 0.483 0.505 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
spresbms_mant <- mantel(spresbms_dist, envbms_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,10,12:15,17,21:23),19])
spresbms_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresbms_dist, ydis = envbms_dist, method = "spearman",      permutations = 999, strata = env[c(1:6, 10, 12:15, 17, 21:23),          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2892 
    ##       Significance: 0.345 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.413 0.467 0.502 0.537 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
spresbom_mant <- mantel(spresbom_dist, envbom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:7,9:11,13:16,18:19,23),19])
spresbom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresbom_dist, ydis = envbom_dist, method = "spearman",      permutations = 999, strata = env[c(3, 6:7, 9:11, 13:16, 18:19,          23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2261 
    ##       Significance: 0.248 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.301 0.339 0.378 0.430 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
spresbso_mant <- mantel(spresbso_dist, envbso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
spresbso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresbso_dist, ydis = envbso_dist, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5571 
    ##       Significance: 0.426 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.600 0.613 0.626 0.645 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Traits
# All sites 1,2,4,12:13
FDpa_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[-20,c(1)]))
jpresta_mant_s <- mantel(jpres_dist, FDpa_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-20,19])
jpresta_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_s, method = "spearman",      permutations = 999, strata = env[-20, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5553 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.370 0.413 0.443 0.472 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[-20,c(2)]))
jpresta_mant_l <- mantel(jpres_dist, FDpa_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-20,19])
jpresta_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_l, method = "spearman",      permutations = 999, strata = env[-20, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3662 
    ##       Significance: 0.004 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.259 0.285 0.303 0.330 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[-20,c(4)]))
jpresta_mant_dmin <- mantel(jpres_dist, FDpa_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-20,19])
jpresta_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_dmin, method = "spearman",      permutations = 999, strata = env[-20, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6215 
    ##       Significance: 0.003 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.555 0.572 0.585 0.604 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpa_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[-20,c(12:13)]))
jpresta_mant_ec <- mantel(jpres_dist, FDpa_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-20,19])
jpresta_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpres_dist, ydis = FDpa_dist_ec, method = "spearman",      permutations = 999, strata = env[-20, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5011 
    ##       Significance: 0.006 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.384 0.419 0.444 0.471 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jpresta_mant_pv <- rbind(jpresta_mant_s$signif, jpresta_mant_l$signif, jpresta_mant_dmin$signif, jpresta_mant_ec$signif)
jpresta_mant_pv <- jpresta_mant_pv[,1]
jpresta_mant_pv <- p.adjust(jpresta_mant_pv, method = "bonferroni")
jpresta_mant_pv
```

    ## [1] 0.008 0.016 0.012 0.024

``` r
# Mixed and stratified lakes 1,4,12:13
FDpams_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[c(1:6,8,10,12:15,17,21:23),c(1)]))
jprestms_mant_s <- mantel(jpresms_dist, FDpams_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
jprestms_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpams_dist_s, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4835 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.312 0.361 0.396 0.428 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpams_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[c(1:6,8,10,12:15,17,21:23),c(4)]))
jprestms_mant_dmin <- mantel(jpresms_dist, FDpams_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
jprestms_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpams_dist_dmin, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5526 
    ##       Significance: 0.009 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.472 0.496 0.526 0.544 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpams_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[c(1:6,8,10,12:15,17,21:23),c(12:13)]))
jprestms_mant_ec <- mantel(jpresms_dist, FDpams_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
jprestms_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresms_dist, ydis = FDpams_dist_ec, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5043 
    ##       Significance: 0.003 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.343 0.380 0.429 0.451 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jprestms_mant_pv <- rbind(jprestms_mant_s$signif, jprestms_mant_dmin$signif, jprestms_mant_ec$signif)
jprestms_mant_pv <- jprestms_mant_pv[,1]
jprestms_mant_pv <- p.adjust(jprestms_mant_pv, method = "bonferroni")
jprestms_mant_pv
```

    ## [1] 0.003 0.027 0.009

``` r
# Ocean sites and mixed lakes 1,2
FDpaom_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[c(3,6:11,13:16,18:19,23),c(1)]))
jprestom_mant_s <- mantel(jpresom_dist, FDpaom_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:11,13:16,18:19,23),19])
jprestom_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresom_dist, ydis = FDpaom_dist_s, method = "spearman",      permutations = 999, strata = env[c(3, 6:11, 13:16, 18:19,          23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4451 
    ##       Significance: 0.003 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.165 0.210 0.263 0.326 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpaom_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[c(3,6:11,13:16,18:19,23),c(2)]))
jprestom_mant_l <- mantel(jpresom_dist, FDpaom_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:11,13:16,18:19,23),19])
jprestom_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresom_dist, ydis = FDpaom_dist_l, method = "spearman",      permutations = 999, strata = env[c(3, 6:11, 13:16, 18:19,          23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2519 
    ##       Significance: 0.022 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.140 0.189 0.225 0.282 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jprestom_mant_pv <- rbind(jprestom_mant_s$signif, jprestom_mant_l$signif)
jprestom_mant_pv <- jprestom_mant_pv[,1]
jprestom_mant_pv <- p.adjust(jprestom_mant_pv, method = "bonferroni")
jprestom_mant_pv
```

    ## [1] 0.006 0.044

``` r
# Stratified lakes and ocean sites 2,4,12:13
FDpaso_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(2)]))
jprestso_mant_l <- mantel(jpresso_dist, FDpaso_dist_l, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
jprestso_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpaso_dist_l, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2324 
    ##       Significance: 0.005 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.142 0.168 0.190 0.208 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpaso_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(4)]))
jprestso_mant_dmin <- mantel(jpresso_dist, FDpaso_dist_dmin, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
jprestso_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpaso_dist_dmin, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5981 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.521 0.537 0.549 0.564 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpaso_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(12:13)]))
jprestso_mant_ec <- mantel(jpresso_dist, FDpaso_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
jprestso_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = jpresso_dist, ydis = FDpaso_dist_ec, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2447 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.149 0.182 0.192 0.217 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
jprestso_mant_pv <- rbind(jprestso_mant_l$signif, jprestso_mant_dmin$signif, jprestso_mant_ec$signif)
jprestso_mant_pv <- jprestso_mant_pv[,1]
jprestso_mant_pv <- p.adjust(jprestso_mant_pv, method = "bonferroni")
jprestso_mant_pv
```

    ## [1] 0.015 0.006 0.006

``` r
## Dice-Sørensen
# All sites
FDpa_dist <- gowdis(as.data.frame(scaled_FD_total_env[-20,c(1:17)]))
spresta_mant <- mantel(spres_dist, FDpa_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-20,19])
spresta_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spres_dist, ydis = FDpa_dist, method = "spearman",      permutations = 999, strata = env[-20, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6331 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.508 0.539 0.562 0.577 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
FDpams_dist <- gowdis(as.data.frame(scaled_FD_total_env[c(1:6,8,10,12:15,17,21:23),c(1:17)]))
sprestms_mant <- mantel(spresms_dist, FDpams_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
sprestms_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresms_dist, ydis = FDpams_dist, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6558 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.460 0.499 0.539 0.570 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
FDpaom_dist <- gowdis(as.data.frame(scaled_FD_total_env[c(3,6:11,13:16,18:19,23),c(1:7)]))
sprestom_mant <- mantel(spresom_dist, FDpaom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:11,13:16,18:19,23),19])
sprestom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresom_dist, ydis = FDpaom_dist, method = "spearman",      permutations = 999, strata = env[c(3, 6:11, 13:16, 18:19,          23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3433 
    ##       Significance: 0.019 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.156 0.225 0.298 0.369 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
FDpaso_dist <- gowdis(as.data.frame(scaled_FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(1:17)]))
sprestso_mant <- mantel(spresso_dist, FDpaso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
sprestso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = spresso_dist, ydis = FDpaso_dist, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3897 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.272 0.310 0.335 0.360 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

## Plot Jaccard incidence NMDS and envfit results

Includes a plot of correlated environmental variables and a plot of
correlated traits

``` r
jpres_NMDS_data.scores <- as.data.frame(scores(jpres_NMDS))
jpres_NMDS_data.scores$Stratification <- env[-20,19]
jpres_NMDS_data.scores$Lakes <- env[-20,1]
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
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), level = 0.50) +
  geom_point(data = jpres_NMDS_data.scores, aes(color = Stratification), size = 4, alpha = 1) + 
  geom_text_repel(data = jpres_NMDS_data.scores, label = jpres_NMDS_data.scores$Lakes, size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = jprese_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#CD661D") +
  geom_text(data = jprese_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
            color = "#CD661D", label = row.names(jprese_ef_coord_cont), size = 7) +
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
  #              data = jpresb_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#698B22") +
  # geom_text(data = jpresb_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
  #           color = "#698B22", label = row.names(jpresb_ef_coord_cont), size = 6) +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
        panel.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  # annotate(geom = "label", x = 1, y = 0.73, size = 6, 
  #          label = paste("Stress: ", round(jpres_NMDS$stress, digits = 2))) +
  labs(color = "Stratification")
jpres_ef_plot <- jpres_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
jpres_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20Jaccard%20incidence%20NMDS%20and%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jpres_ef_plot.png", jpres_ef_plot, width = 8, height = 4, units = "in")


jprest_ef_coord_cont <- as.data.frame(scores(jprest_ef, "vectors")) * ordiArrowMul(jprest_ef)
jprest_ef_coord_cat = as.data.frame(scores(jprest_ef, "factors")) * ordiArrowMul(jprest_ef)

jprest_ef_plot <- ggplot(data = jpres_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification),level = 0.50) +
  geom_point(data = jpres_NMDS_data.scores, aes(color = Stratification), size = 4, alpha = 1) + 
  geom_text_repel(data = jpres_NMDS_data.scores, label = jpres_NMDS_data.scores$Lakes, size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = jprest_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#CD3333") +
  geom_text(data = jprest_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
            color = "#CD3333", label = row.names(jprest_ef_coord_cont), size = 7) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = jprest_ef_coord_cat, linewidth =1, alpha = 0.2, color = "#8B6508") +
  geom_text(data = jprest_ef_coord_cat, aes(x = NMDS1, y = NMDS2), 
            color = "#8B6508", label = row.names(jprest_ef_coord_cat), size = 7, check_overlap = T) + 
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
        panel.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  # annotate(geom = "label", x = 1.05, y = 0.73, size = 6, 
  #          label = paste("Stress: ", round(jpres_NMDS$stress, digits = 2))) +
  labs(color = "Stratification", tag = "B")
jprest_ef_plot <- jprest_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
jprest_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20Jaccard%20incidence%20NMDS%20and%20envfit%20results-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/jprest_ef_plot.png", jprest_ef_plot, width = 8.25, height = 4.13, units = "in")
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
  geom_point(data = spres_NMDS_data.scores, aes(color = Stratification), size = 7, alpha = 1) + 
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
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
        panel.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  annotate(geom = "label", x = -1.5, y = 0.82639346, size = 4, 
           label = paste("Stress: ", round(spres_NMDS$stress, digits = 2))) +
  labs(color = "Stratification")
spres_ef_plot
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/spres_ef_plot.png", spres_ef_plot, width = 8, height = 4, units = "in")
```

## Trait dissimilarity distances

We use the community weighted means calculated from the FD package. We
use traits that were found to be significantly correlated with
stratification and environmental variables.

``` r
### Regular 6,8,10,18,20,21,22
# All sites 6,8,10,18,20,21,22
# Non-siginficant dipersion but significant difference between stratification types: 10,20,22
# BodyShapeI, OperculumPresent, Troph, RepGuild1, ParentalCare, WaterPref, DorsalSpinesMean
fd_dist <- gowdis(as.data.frame(FD_total_env[-20,c(6,8,10,18,20,21,22)]))
# Mixed and stratified lakes
fdms_dist <- gowdis(FD_total_env[c(1:6,8,10,12:15,17,21:23),c(6,8,10,18,20,21,22)])
# Ocean sites and mixed lakes
fdom_dist <- gowdis(FD_total_env[c(3,6:11,13:16,18:19,23),c(6,8,10,18,20,21,22)])
# Stratified lakes and ocean sites
fdso_dist <- gowdis(FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(6,8,10,18,20,21,22)])
# Stratified lakes
fds_dist <- gowdis(FD_total_env[c(1:2,4,5,12,17,21:22),c(6,8,10,18,20,21,22)])


### Environmental -c(9,16,18,20)
# All sites
fde_dist <- gowdis(FD_total_env[-c(9,16,18,20),c(6,8,10,18,20,21,22)])
# Mixed and stratified lakes
fdems_dist <- gowdis(FD_total_env[c(1:6,8,10,12:15,17,21:23),c(6,8,10,18,20,21,22)])
# Ocean sites and mixed lakes
fdeom_dist <- gowdis(FD_total_env[c(3,6:8,10:11,13:15,19,23),c(6,8,10,18,20,21,22)])
# Stratified lakes and ocean sites
fdeso_dist <- gowdis(FD_total_env[c(1:2,4,5,7,11:12,17,19,21:22),c(6,8,10,18,20,21,22)])


### Biogeographic -c(8,20)
# All sites
fdb_dist <- gowdis(FD_total_env[-c(8,20),-c(6,8,10,18,20,21,22)])
# Mixed and stratified lakes
fdbms_dist <- gowdis(FD_total_env[c(1:6,10,12:15,17,21:23),c(6,8,10,18,20,21,22)])
# Ocean sites and mixed lakes
fdbom_dist <- gowdis(FD_total_env[c(3,6:7,9:11,13:16,18:19,23),c(6,8,10,18,20,21,22)])
# Stratified lakes and ocean sites
fdbso_dist <- gowdis(FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(6,8,10,18,20,21,22)])
```

## Trait NMDS

To constrain dissimilarities we perform Nonmetric Multidimensional
Scaling (NMDS), which tries to find a stable solution using the metaMDS
package.

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


### Environmental
# All sites 
fde_NMDS <- metaMDS(fde_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
fdems_NMDS <- metaMDS(fdems_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
fdeom_NMDS <- metaMDS(fdeom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
fdeso_NMDS <- metaMDS(fdeso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)


### Biogeographic
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

We use betadisper to determine homogeneity of the trait data based on
the stratification category. Is the dispersion of traits similar between
the different stratifications.

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
    ##           Df  Sum Sq  Mean Sq F value    Pr(>F)    
    ## Groups     2 0.39642 0.198210   30.89 1.069e-06 ***
    ## Residuals 19 0.12192 0.006417                      
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
    ##           Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)    
    ## Groups     2 0.39642 0.198210 30.89    999  0.001 ***
    ## Residuals 19 0.12192 0.006417                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##                 Mixed      Ocean Stratified
    ## Mixed                 3.7000e-02      0.001
    ## Ocean      4.0346e-02                 0.001
    ## Stratified 6.7195e-06 6.5666e-04

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
    ## Ocean-Mixed      0.05204806 -0.05785458 0.1619507 0.4657824
    ## Stratified-Mixed 0.29806649  0.19631641 0.3998166 0.0000014
    ## Stratified-Ocean 0.24601843  0.13611579 0.3559211 0.0000501

``` r
boxplot(stratification_fd_bd)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Determine%20stratification%20homogeneity%20using%20traits-1.png)<!-- -->

## Determine stratification dissimilarity using traits

We use adonis to determine if dissimilarities of species traits by
stratification categories are significant and how much of the variation
is explained by trait dissimilarities.

``` r
### Stratification
# All sites 
fd_pms <- adonis2(fd_dist ~ env[-20,19], permutations = 999)
fd_pms
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = fd_dist ~ env[-20, 19], permutations = 999)
    ##              Df SumOfSqs      R2     F Pr(>F)   
    ## env[-20, 19]  2   1.2750 0.55876 12.03  0.002 **
    ## Residual     19   1.0068 0.44124                
    ## Total        21   2.2818 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratifcation
fd_pms_pair <- pairwise.adonis(fd_dist, env[-20,19], p.adjust.m = "bonferroni", perm = 999)
fd_pms_pair
```

    ##                 pairs Df  SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 0.91580033 13.423218 0.4894837   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 0.93071542 11.199063 0.4827377   0.005      0.015   .
    ## 3      Mixed vs Ocean  1 0.02619034  5.132228 0.2995657   0.010      0.030   .

## Envfit environmental influence on traits

We use envfit to determine significantly correlated environmental
variables to our trait NMDS. We will use these results, if significant,
in our figure.

``` r
### Environmental
# All sites for figure
fde_ef <- envfit(fde_NMDS, env[-c(9,16,18,20),c(34)], permutations = 1000, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
fde_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##          NMDS1     NMDS2    r2  Pr(>r)  
    ## [1,]  0.998300 -0.058296 0.838 0.02897 *
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
    ##          NMDS1     NMDS2    r2  Pr(>r)  
    ## [1,]  0.998300 -0.058296 0.838 0.02897 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
fdea_ef <- envfit(fde_NMDS, env[-c(9,16,18,20),c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
fdea_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median -0.36833 -0.92969 0.0324 0.89810  
    ## salinity_median     0.99830 -0.05830 0.8380 0.02298 *
    ## oxygen_median       0.85525 -0.51821 0.1609 0.99201  
    ## pH_median           0.63014 -0.77648 0.6622 0.03097 *
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
    ## temperature_median -0.36833 -0.92969 0.0324 1.00000  
    ## salinity_median     0.99830 -0.05830 0.8380 0.09191 .
    ## oxygen_median       0.85525 -0.51821 0.1609 1.00000  
    ## pH_median           0.63014 -0.77648 0.6622 0.12388  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
fdems_ef <- envfit(fdems_NMDS, env[c(1:6,8,10,12:15,17,21:23),c(6)], permutations = 1000, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
fdems_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##        NMDS1   NMDS2     r2  Pr(>r)  
    ## [1,] 0.99444 0.10530 0.8207 0.02198 *
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
    ##        NMDS1   NMDS2     r2  Pr(>r)  
    ## [1,] 0.99444 0.10530 0.8207 0.02198 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
fdeom_ef <- envfit(fdeom_NMDS, env[c(3,6:8,10:11,13:15,19,23),c(2,6,8,10)], permutations = 1000, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
fdeom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2  Pr(>r)  
    ## temperature_median  0.28445 -0.95869 0.0400 0.75924  
    ## salinity_median     0.99869 -0.05117 0.5614 0.09091 .
    ## oxygen_median       0.98315 -0.18282 0.5159 0.21978  
    ## pH_median           0.99360 -0.11295 0.3615 0.09491 .
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
    ##                       NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median  0.28445 -0.95869 0.0400 1.0000
    ## salinity_median     0.99869 -0.05117 0.5614 0.3636
    ## oxygen_median       0.98315 -0.18282 0.5159 0.8791
    ## pH_median           0.99360 -0.11295 0.3615 0.3796
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
fdeso_ef <- envfit(fdeso_NMDS, env[c(1:2,4,5,7,11:12,17,19,21:22),c(6)], permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,11:12,17,19,21:22),19])
fdeso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## [1,]  0.90449 -0.42649 0.8323 0.01798 *
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
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## [1,]  0.90449 -0.42649 0.8323 0.01798 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes
fds_ef <- envfit(fds_NMDS, env[c(1:2,4,5,12,17,21:22),c(28)], permutations = 1000, na.rm = TRUE)
fds_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2   Pr(>r)   
    ## [1,]  0.35484 -0.93493 0.8808 0.002997 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ##         NMDS1    NMDS2     r2   Pr(>r)   
    ## [1,]  0.35484 -0.93493 0.8808 0.002997 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Biogeographic
# All sites for figure
fdb_ef <- envfit(fdb_NMDS, env[-c(8,20),c(39,40)], permutations = 1000, na.rm = TRUE, strata = env[-c(8,20),19])
fdb_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## A     0.79065 -0.61227 0.2147 0.08292 .
    ## minD -0.99118  0.13250 0.7123 0.14286  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ##         NMDS1    NMDS2     r2 Pr(>r)
    ## A     0.79065 -0.61227 0.2147 0.1658
    ## minD -0.99118  0.13250 0.7123 0.2857
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
fdba_ef <- envfit(fdb_NMDS, env[-c(8,20),c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[-c(8,20),19])
fdba_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.94400 -0.32996 0.1475 0.18881  
    ## volume_m3                   0.92109 -0.38934 0.1368 0.20280  
    ## surface_area_m2             0.79065 -0.61227 0.2147 0.08691 .
    ## distance_to_ocean_min_m    -0.99118  0.13250 0.7123 0.16084  
    ## distance_to_ocean_mean_m   -0.91483  0.40385 0.6426 0.42757  
    ## distance_to_ocean_median_m -0.91390  0.40595 0.6420 0.42458  
    ## tidal_lag_time_minutes     -0.81089  0.58520 0.5946 0.78821  
    ## tidal_efficiency            0.74628 -0.66563 0.6686 0.23676  
    ## perimeter_fromSat           0.80691 -0.59068 0.2054 0.10390  
    ## max_depth                  -0.95635  0.29224 0.0859 0.91908  
    ## logArea                     0.58892 -0.80819 0.1041 0.21179  
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
    ## volume_m3_w_chemocline      0.94400 -0.32996 0.1475  1.000
    ## volume_m3                   0.92109 -0.38934 0.1368  1.000
    ## surface_area_m2             0.79065 -0.61227 0.2147  0.956
    ## distance_to_ocean_min_m    -0.99118  0.13250 0.7123  1.000
    ## distance_to_ocean_mean_m   -0.91483  0.40385 0.6426  1.000
    ## distance_to_ocean_median_m -0.91390  0.40595 0.6420  1.000
    ## tidal_lag_time_minutes     -0.81089  0.58520 0.5946  1.000
    ## tidal_efficiency            0.74628 -0.66563 0.6686  1.000
    ## perimeter_fromSat           0.80691 -0.59068 0.2054  1.000
    ## max_depth                  -0.95635  0.29224 0.0859  1.000
    ## logArea                     0.58892 -0.80819 0.1041  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
fdbms_ef <- envfit(fdbms_NMDS, env[c(1:6,10,12:15,17,21:23),c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[c(1:6,10,12:15,17,21:23),19])
fdbms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline     -0.00335 -0.99999 0.1022 0.6284
    ## volume_m3                  -0.23405 -0.97222 0.2640 0.4815
    ## surface_area_m2            -0.30354 -0.95282 0.2595 0.4436
    ## distance_to_ocean_min_m    -0.80427 -0.59426 0.4986 0.1818
    ## distance_to_ocean_mean_m   -0.51381 -0.85790 0.4820 0.2288
    ## distance_to_ocean_median_m -0.52684 -0.84996 0.5015 0.1978
    ## tidal_lag_time_minutes     -0.66990 -0.74246 0.1822 0.8771
    ## tidal_efficiency            0.89055  0.45489 0.2443 0.5524
    ## perimeter_fromSat          -0.47886 -0.87789 0.1380 0.4925
    ## max_depth                  -0.31858 -0.94790 0.1382 0.6773
    ## logArea                    -0.28522 -0.95846 0.1268 0.4955
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
    ## volume_m3_w_chemocline     -0.00335 -0.99999 0.1022      1
    ## volume_m3                  -0.23405 -0.97222 0.2640      1
    ## surface_area_m2            -0.30354 -0.95282 0.2595      1
    ## distance_to_ocean_min_m    -0.80427 -0.59426 0.4986      1
    ## distance_to_ocean_mean_m   -0.51381 -0.85790 0.4820      1
    ## distance_to_ocean_median_m -0.52684 -0.84996 0.5015      1
    ## tidal_lag_time_minutes     -0.66990 -0.74246 0.1822      1
    ## tidal_efficiency            0.89055  0.45489 0.2443      1
    ## perimeter_fromSat          -0.47886 -0.87789 0.1380      1
    ## max_depth                  -0.31858 -0.94790 0.1382      1
    ## logArea                    -0.28522 -0.95846 0.1268      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
fdbom_ef <- envfit(fdbom_NMDS, env[c(3,6:7,9:11,13:16,18:19,23),c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[c(3,6:7,9:11,13:16,18:19,23),19])
fdbom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.13170 -0.99129 0.2733 0.4056
    ## volume_m3                   0.13131 -0.99134 0.2728 0.4066
    ## surface_area_m2             0.15169 -0.98843 0.3040 0.4276
    ## distance_to_ocean_min_m    -0.87057  0.49204 0.3997 0.4096
    ## distance_to_ocean_mean_m   -0.81767  0.57568 0.3986 0.5375
    ## distance_to_ocean_median_m -0.83592  0.54885 0.3857 0.5674
    ## tidal_lag_time_minutes     -0.89858  0.43882 0.2671 0.5425
    ## tidal_efficiency            0.91859 -0.39521 0.2372 0.6314
    ## perimeter_fromSat           0.26271 -0.96487 0.2274 0.5984
    ## max_depth                   0.12579 -0.99206 0.2460 0.1878
    ## logArea                     0.31719 -0.94836 0.2397 0.4815
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
    ## volume_m3_w_chemocline      0.13170 -0.99129 0.2733      1
    ## volume_m3                   0.13131 -0.99134 0.2728      1
    ## surface_area_m2             0.15169 -0.98843 0.3040      1
    ## distance_to_ocean_min_m    -0.87057  0.49204 0.3997      1
    ## distance_to_ocean_mean_m   -0.81767  0.57568 0.3986      1
    ## distance_to_ocean_median_m -0.83592  0.54885 0.3857      1
    ## tidal_lag_time_minutes     -0.89858  0.43882 0.2671      1
    ## tidal_efficiency            0.91859 -0.39521 0.2372      1
    ## perimeter_fromSat           0.26271 -0.96487 0.2274      1
    ## max_depth                   0.12579 -0.99206 0.2460      1
    ## logArea                     0.31719 -0.94836 0.2397      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
fdbso_ef <- envfit(fdbso_NMDS, env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
fdbso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2 Pr(>r)
    ## volume_m3_w_chemocline      0.99010  0.14034 0.1084 0.5275
    ## volume_m3                   0.99740  0.07205 0.1000 0.5924
    ## surface_area_m2             0.99882 -0.04862 0.1525 0.7223
    ## distance_to_ocean_min_m    -0.99978 -0.02093 0.5928 0.1908
    ## distance_to_ocean_mean_m   -0.99024 -0.13936 0.5002 0.4136
    ## distance_to_ocean_median_m -0.98841 -0.15182 0.5156 0.3686
    ## tidal_lag_time_minutes     -0.84934  0.52784 0.3773 0.9770
    ## tidal_efficiency            0.80556 -0.59251 0.5135 0.6374
    ## perimeter_fromSat           0.94537 -0.32599 0.1479 0.7632
    ## max_depth                  -0.82197 -0.56953 0.0417 0.9081
    ## logArea                     0.56914 -0.82224 0.0960 0.5934
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
    ## volume_m3_w_chemocline      0.99010  0.14034 0.1084      1
    ## volume_m3                   0.99740  0.07205 0.1000      1
    ## surface_area_m2             0.99882 -0.04862 0.1525      1
    ## distance_to_ocean_min_m    -0.99978 -0.02093 0.5928      1
    ## distance_to_ocean_mean_m   -0.99024 -0.13936 0.5002      1
    ## distance_to_ocean_median_m -0.98841 -0.15182 0.5156      1
    ## tidal_lag_time_minutes     -0.84934  0.52784 0.3773      1
    ## tidal_efficiency            0.80556 -0.59251 0.5135      1
    ## perimeter_fromSat           0.94537 -0.32599 0.1479      1
    ## max_depth                  -0.82197 -0.56953 0.0417      1
    ## logArea                     0.56914 -0.82224 0.0960      1
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
    ## Abudefduf_lorenzi               -0.23611 -0.97173 0.0149 0.940060   
    ## Abudefduf_septemfasciatus        0.99529  0.09692 0.0169 0.640360   
    ## Abudefduf_sexfasciatus          -0.23611 -0.97173 0.0149 0.940060   
    ## Acanthurus_lineatus              0.63014  0.77648 0.0213 0.414585   
    ## Acanthurus_nigricauda            0.97617  0.21703 0.0556 0.547453   
    ## Acanthurus_sp                    0.90500 -0.42542 0.0199 0.500500   
    ## Acanthurus_xanthopterus          0.75007  0.66135 0.1882 0.145854   
    ## Acentrogobius_janthinopterus    -0.83826  0.54527 0.3818 0.007992 **
    ## Amblyeleotris_gymnocephala       0.95129  0.30831 0.0155 0.724276   
    ## Amblyglyphidodon_curacao         0.60749 -0.79432 0.0382 0.665335   
    ## Amblygobius_buanensis            0.85468  0.51915 0.3321 0.004995 **
    ## Amblygobius_decussatus           0.25879 -0.96593 0.0856 0.219780   
    ## Amblygobius_esakiae              0.90500 -0.42542 0.0199 0.500500   
    ## Amblygobius_linki                0.66600  0.74595 0.0187 0.606394   
    ## Amblygobius_phalaena             0.97146  0.23722 0.0315 0.830170   
    ## Amphiprion_clarkii               0.96122 -0.27578 0.0202 0.476523   
    ## Ancistrogobius_dipus             0.71129  0.70290 0.0134 0.825175   
    ## Apogonichthyoides_melas          0.71129  0.70290 0.0134 0.825175   
    ## Arcygobius_baliurus              0.81538  0.57892 0.0155 0.688312   
    ## Arothron_reticularis             0.74779  0.66393 0.0129 0.842158   
    ## Asterropteryx_semipunctata       0.99446  0.10514 0.0539 0.545455   
    ## Asterropteryx_sp                 0.66600  0.74595 0.0187 0.606394   
    ## Atherinomorus_lacunosus          0.97416  0.22586 0.0160 0.878122   
    ## Atule_mate                       0.63014  0.77648 0.0213 0.414585   
    ## Balistoides_viridescens          0.89129  0.45343 0.1280 0.292707   
    ## Bolbometopon_muricatum           0.36612 -0.93057 0.0885 0.306693   
    ## Caesio_caerulaurea               0.71452  0.69961 0.0341 0.708292   
    ## Caesio_cuning                    0.48565 -0.87415 0.1028 0.375624   
    ## Caesio_teres                     0.95129  0.30831 0.0155 0.724276   
    ## Cantherhines_dumerilii           0.96122 -0.27578 0.0202 0.476523   
    ## Canthigaster_bennetti            0.96122 -0.27578 0.0202 0.476523   
    ## Canthigaster_solandri            0.99529  0.09692 0.0169 0.640360   
    ## Caranx_ignobilis                 0.98712  0.15999 0.0146 0.767233   
    ## Caranx_melampygus                0.38841 -0.92149 0.0053 0.989011   
    ## Caranx_papuensis                 0.95129  0.30831 0.0155 0.724276   
    ## Caranx_sexfasciatus              0.80520  0.59300 0.0939 0.430569   
    ## Cephalopholis_argus             -0.44267 -0.89668 0.2031 0.184815   
    ## Cephalopholis_boenak             0.93543 -0.35352 0.0421 0.561439   
    ## Chaetodon_auriga                 0.96905 -0.24685 0.2867 0.032967 * 
    ## Chaetodon_bennetti               0.90818  0.41857 0.0416 0.578422   
    ## Chaetodon_ephippium              0.89897 -0.43801 0.2449 0.057942 . 
    ## Chaetodon_kleinii                0.92237  0.38631 0.0607 0.495504   
    ## Chaetodon_lunula                 0.87877  0.47725 0.1594 0.166833   
    ## Chaetodon_lunulatus              0.75304 -0.65798 0.1443 0.210789   
    ## Chaetodon_melannotus             0.71723  0.69683 0.0193 0.538462   
    ## Chaetodon_ocellicaudus           0.98712  0.15999 0.0146 0.767233   
    ## Chaetodon_octofasciatus          0.93543 -0.35352 0.0421 0.561439   
    ## Chaetodon_rafflesii              0.98855  0.15092 0.0850 0.421578   
    ## Chaetodon_semeion                0.95129  0.30829 0.0401 0.580420   
    ## Chaetodon_trifascialis           0.96122 -0.27578 0.0202 0.476523   
    ## Chaetodon_ulietensis             0.51780 -0.85550 0.1000 0.267732   
    ## Chaetodon_vagabundus             0.97243 -0.23319 0.2893 0.027972 * 
    ## Chaetodontoplus_poliourus        0.96122 -0.27578 0.0202 0.476523   
    ## Cheilinus_fasciatus              0.98894  0.14830 0.0647 0.466533   
    ## Cheilinus_trilobatus             0.63014  0.77648 0.0213 0.414585   
    ## Cheilinus_undulatus              0.25412 -0.96717 0.0745 0.311688   
    ## Cheilodipterus_artus             0.21337 -0.97697 0.0100 0.873127   
    ## Cheilodipterus_isostigma         0.27332 -0.96192 0.0694 0.363636   
    ## Cheilodipterus_quinquelineatus   0.90937 -0.41599 0.2468 0.049950 * 
    ## Cheilodipterus_singapurensis     0.75449  0.65631 0.0320 0.799201   
    ## Chlorurus_bleekeri               0.51451 -0.85749 0.1093 0.283716   
    ## Chlorurus_microrhinos            0.99678 -0.08021 0.0361 0.686314   
    ## Chlorurus_spilurus               0.98855  0.15092 0.0850 0.421578   
    ## Choerodon_anchorago              0.99634 -0.08543 0.3920 0.002997 **
    ## Chromis_viridis                  0.42127 -0.90693 0.0765 0.372627   
    ## Chrysiptera_biocellata           0.93034  0.36669 0.0338 0.741259   
    ## Chrysiptera_oxycephala           0.35041 -0.93660 0.0784 0.496503   
    ## Corythoichthys_ocellatus         0.96122 -0.27578 0.0202 0.476523   
    ## Cryptocentrus_caeruleomaculatus  0.67547 -0.73739 0.0692 0.554446   
    ## Cryptocentrus_cyanospilotus      0.20811  0.97811 0.0024 1.000000   
    ## Cryptocentrus_leptocephalus      0.91608  0.40100 0.0314 0.819181   
    ## Cryptocentrus_strigilliceps      0.51006 -0.86014 0.0988 0.331668   
    ## Ctenochaetus_striatus            0.16568 -0.98618 0.0700 0.526474   
    ## Ctenogobiops_pomastictus         0.71129  0.70290 0.0134 0.825175   
    ## Cymbacephalus_beauforti          0.98712  0.15999 0.0146 0.767233   
    ## Dascyllus_aruanus                0.59229 -0.80573 0.0848 0.422577   
    ## Dascyllus_trimaculatus           0.96122 -0.27578 0.0202 0.476523   
    ## Diplogrammus_goramensis          0.87289  0.48792 0.0290 0.894106   
    ## Diproctacanthus_xanthurus        0.90818  0.41857 0.0416 0.578422   
    ## Dischistodus_chrysopoecilus      0.88446  0.46662 0.0280 0.922078   
    ## Dischistodus_perspicillatus      0.80269 -0.59639 0.1279 0.297702   
    ## Doboatherina_duodecimalis        0.28051 -0.95985 0.1295 0.283716   
    ## Eleotris_fusca                  -0.30761  0.95151 0.3599 0.030969 * 
    ## Epibulus_brevis                  0.60937 -0.79288 0.1286 0.206793   
    ## Epinephelus_coeruleopunctatus    0.73534  0.67770 0.0358 0.714286   
    ## Epinephelus_merra                0.53707 -0.84354 0.1002 0.307692   
    ## Epinephelus_sp                   0.73002  0.68342 0.0125 0.915085   
    ## Eviota_atriventris               0.96020  0.27931 0.0527 0.561439   
    ## Eviota_bifasciata                0.82172 -0.56989 0.1761 0.133866   
    ## Eviota_fallax                    0.73395 -0.67920 0.1533 0.175824   
    ## Eviota_lachdeberei               0.85288  0.52211 0.0721 0.561439   
    ## Eviota_maculosa                  0.13341 -0.99106 0.0988 0.305694   
    ## Eviota_sigillata                 0.83456  0.55091 0.0361 0.669331   
    ## Eviota_sp                        0.66600  0.74595 0.0187 0.606394   
    ## Eviota_storthynx                 0.84710  0.53144 0.0539 0.540460   
    ## Exyrias_belissimus               0.91234 -0.40944 0.1514 0.197802   
    ## Exyrias_puntang                 -0.04985  0.99876 0.2867 0.027972 * 
    ## Favonigobius_reichei             0.91608  0.40100 0.0314 0.819181   
    ## Fibramia_lateralis               0.04556 -0.99896 0.0232 0.829171   
    ## Fibramia_thermalis               0.66627 -0.74571 0.0970 0.399600   
    ## Fusigobius_signipinnis           0.99878  0.04928 0.0549 0.538462   
    ## Gerres_erythrourus               0.98712  0.15999 0.0146 0.767233   
    ## Gerres_oyena                     0.99184  0.12747 0.0330 0.768232   
    ## Gladiogobius_ensifer             0.13341 -0.99106 0.0988 0.305694   
    ## Gnathanodon_speciosus            0.49547  0.86862 0.0130 0.989011   
    ## Gobiodon_okinawae                0.63014  0.77648 0.0213 0.414585   
    ## Gymnothorax_pictus               0.93034  0.36669 0.0338 0.741259   
    ## Halichoeres_chloropterus         0.77497 -0.63199 0.0536 0.623377   
    ## Halichoeres_leucurus             0.98248  0.18639 0.0862 0.418581   
    ## Halichoeres_marginatus           0.98712  0.15999 0.0146 0.767233   
    ## Halichoeres_melanurus            0.18559 -0.98263 0.0488 0.622378   
    ## Halichoeres_papilionaceus        0.80965  0.58692 0.0555 0.541459   
    ## Halichoeres_trimaculatus         0.99184  0.12747 0.0330 0.768232   
    ## Hemiglyphidodon_plagiometopon    0.45128 -0.89238 0.1069 0.351648   
    ## Hemigymnus_melapterus            0.99678 -0.08021 0.0361 0.686314   
    ## Heniochus_chrysostomus           0.93764  0.34761 0.0595 0.496503   
    ## Heniochus_monoceros              0.84710  0.53144 0.0539 0.540460   
    ## Herklotsichthys_quadrimaculatus  0.98712  0.15999 0.0146 0.767233   
    ## Hipposcarus_longiceps            0.44737 -0.89435 0.1113 0.251748   
    ## Istigobius_decoratus             0.98712  0.15999 0.0146 0.767233   
    ## Koumansetta_rainfordi            0.96122 -0.27578 0.0202 0.476523   
    ## Kyphosus_cinerascens            -0.21320 -0.97701 0.0174 0.925075   
    ## Lethrinus_erythropterus          0.88113  0.47287 0.0682 0.605395   
    ## Lethrinus_harak                  0.83501  0.55023 0.1182 0.365634   
    ## Lethrinus_obsoletus              0.71129  0.70290 0.0134 0.825175   
    ## Lethrinus_olivaceus              0.71129  0.70290 0.0134 0.825175   
    ## Lutjanus_argentimaculatus        0.42544  0.90499 0.1341 0.269730   
    ## Lutjanus_biguttatus             -0.25345 -0.96735 0.0169 0.926074   
    ## Lutjanus_bohar                   0.95129  0.30829 0.0401 0.580420   
    ## Lutjanus_decussatus              0.71129  0.70290 0.0134 0.825175   
    ## Lutjanus_ehrenbergii             0.78768  0.61609 0.0911 0.437562   
    ## Lutjanus_fulviflamma             0.87487  0.48436 0.0859 0.468531   
    ## Lutjanus_fulvus                  0.55440  0.83225 0.3632 0.007992 **
    ## Lutjanus_gibbus                  0.91817  0.39619 0.0762 0.505495   
    ## Lutjanus_monostigma              0.89150  0.45302 0.0285 0.872128   
    ## Lutjanus_russellii               0.66600  0.74595 0.0187 0.606394   
    ## Macrodontogobius_wilburi         0.75304 -0.65798 0.1443 0.210789   
    ## Megalops_cyprinoides            -0.17801  0.98403 0.5178 0.034965 * 
    ## Meiacanthus_ditrema              0.90500 -0.42542 0.0199 0.500500   
    ## Meiacanthus_grammistes           0.93543 -0.35352 0.0421 0.561439   
    ## Monodactylus_argenteus           0.68026  0.73297 0.0316 0.806194   
    ## Monotaxis_grandoculis           -0.02104 -0.99978 0.0121 0.861139   
    ## Mugilogobius_cavifrons          -0.60182  0.79863 0.2146 0.126873   
    ## Mulloidichthys_flavolineatus     0.98712  0.15999 0.0146 0.767233   
    ## Myripristis_adusta               0.81254  0.58291 0.1256 0.317682   
    ## Myripristis_sp                   0.81538  0.57892 0.0155 0.688312   
    ## Myripristis_violacea             0.96122 -0.27578 0.0202 0.476523   
    ## Naso_brevirostris                0.98494  0.17289 0.0776 0.501499   
    ## Naso_vlamingii                   0.87638  0.48163 0.0571 0.509491   
    ## Nectamia_viria                   0.71723  0.69683 0.0193 0.538462   
    ## Neoglyphidodon_melas            -0.44267 -0.89668 0.2031 0.184815   
    ## Neoniphon_argenteus              0.77483  0.63217 0.0508 0.603397   
    ## Neoniphon_opercularis           -0.44267 -0.89668 0.2031 0.184815   
    ## Neoniphon_sammara                0.95550  0.29498 0.0951 0.414585   
    ## Neopomacentrus_nemurus           0.25267 -0.96755 0.0634 0.553447   
    ## Omobranchus_obliquus             0.60665  0.79497 0.0333 0.759241   
    ## Oplopomops_diacanthus            0.98712  0.15999 0.0146 0.767233   
    ## Oxycheilinus_celebicus           0.22029 -0.97544 0.0156 0.816184   
    ## Parapercis_cylindrica            0.98712  0.15999 0.0146 0.767233   
    ## Parioglossus_formosus           -0.10980  0.99395 0.2528 0.083916 . 
    ## Parioglossus_philippinus         0.75917  0.65089 0.0543 0.570430   
    ## Parioglossus_sp                  0.98511  0.17192 0.0366 0.650350   
    ## Parupeneus_barberinus            0.96515  0.26170 0.1647 0.163836   
    ## Parupeneus_multifasciatus        0.88389 -0.46769 0.2084 0.103896   
    ## Pentapodus_trivittatus           0.87597 -0.48237 0.2074 0.106893   
    ## Periophthalmus_argentilineatus   0.99529  0.09692 0.0169 0.640360   
    ## Pervagor_janthinosoma            0.63014  0.77648 0.0213 0.414585   
    ## Pervagor_nigrolineatus           0.37829 -0.92569 0.0786 0.376623   
    ## Petroscirtes_breviceps           0.99811  0.06153 0.1038 0.326673   
    ## Platax_orbicularis               0.76465  0.64444 0.0303 0.868132   
    ## Platax_pinnatus                  0.95129  0.30831 0.0155 0.724276   
    ## Platax_sp                        0.73002  0.68342 0.0125 0.915085   
    ## Plectorhinchus_albovittatus      0.66600  0.74595 0.0187 0.606394   
    ## Plectorhinchus_chaetodonoides    0.13341 -0.99106 0.0988 0.305694   
    ## Plectorhinchus_gibbosus          0.98712  0.15999 0.0146 0.767233   
    ## Plectropomus_leopardus           0.96122 -0.27578 0.0202 0.476523   
    ## Pomacanthus_sexstriatus          0.86214  0.50667 0.0349 0.678322   
    ## Pomacentrus_adelus               0.99529  0.09692 0.0169 0.640360   
    ## Pomacentrus_amboinensis          0.96122 -0.27578 0.0202 0.476523   
    ## Pomacentrus_burroughi            0.84161 -0.54009 0.1823 0.121878   
    ## Pomacentrus_grammorhynchus       0.98712  0.15999 0.0146 0.767233   
    ## Pomacentrus_nigromanus          -0.15416 -0.98805 0.0248 0.876124   
    ## Pomacentrus_pavo                 0.96122 -0.27578 0.0202 0.476523   
    ## Pomacentrus_simsiang             0.92785 -0.37296 0.2872 0.022977 * 
    ## Pomacentrus_sp                  -0.32959 -0.94412 0.0637 0.483516   
    ## Pomacentrus_taeniometopon        0.66600  0.74595 0.0187 0.606394   
    ## Pseudobalistes_flavimarginatus   0.95775  0.28760 0.1323 0.282717   
    ## Pseudochromis_fuscus             0.96122 -0.27578 0.0202 0.476523   
    ## Ptereleotris_brachyptera         1.00000 -0.00184 0.0370 0.638362   
    ## Ptereleotris_evides              0.96122 -0.27578 0.0202 0.476523   
    ## Pterocaesio_tile                 0.71129  0.70290 0.0134 0.825175   
    ## Pterocaesio_trilineata           0.71723  0.69683 0.0193 0.538462   
    ## Pterois_volitans                 0.95129  0.30831 0.0155 0.724276   
    ## Rastrelliger_kanagurta           0.31495 -0.94911 0.0673 0.416583   
    ## Rhabdamia_gracilis               0.83309  0.55314 0.0727 0.536464   
    ## Rhinecanthus_aculeatus           0.88112  0.47290 0.0607 0.483516   
    ## Salarias_segmentatus             0.90427  0.42696 0.1687 0.146853   
    ## Sargocentron_diadema             0.67145  0.74105 0.0425 0.556444   
    ## Sargocentron_spiniferum          0.81538  0.57892 0.0155 0.688312   
    ## Saurida_gracilis                 0.63014  0.77648 0.0213 0.414585   
    ## Scarus_dimidiatus                0.98855  0.15092 0.0850 0.421578   
    ## Scarus_flavipectoralis           0.90818  0.41857 0.0416 0.578422   
    ## Scarus_ghobban                   0.98511  0.17192 0.0366 0.650350   
    ## Scarus_hypselopterus             0.63014  0.77648 0.0213 0.414585   
    ## Scarus_oviceps                   0.90500 -0.42542 0.0199 0.500500   
    ## Scarus_psittacus                 0.98112  0.19339 0.1072 0.318681   
    ## Scarus_quoyi                     0.96122 -0.27578 0.0202 0.476523   
    ## Scarus_rivulatus                 0.25879 -0.96593 0.0856 0.219780   
    ## Scatophagus_argus                0.74779  0.66393 0.0129 0.842158   
    ## Scolopsis_ciliata                0.73156 -0.68177 0.1706 0.139860   
    ## Scolopsis_lineata                0.91608  0.40100 0.0314 0.819181   
    ## Scolopsis_margaritifera          0.35818 -0.93365 0.0807 0.458541   
    ## Scolopsis_trilineata             0.74712  0.66469 0.0530 0.575425   
    ## Siganus_doliatus                 0.36612 -0.93057 0.0885 0.306693   
    ## Siganus_fuscescens               0.77733  0.62909 0.0530 0.568432   
    ## Siganus_lineatus                 0.58173  0.81338 0.2995 0.029970 * 
    ## Siganus_puellus                  0.35551 -0.93467 0.0998 0.298701   
    ## Siganus_punctatissimus           0.96122 -0.27578 0.0202 0.476523   
    ## Siganus_punctatus                0.96122 -0.27578 0.0202 0.476523   
    ## Signigobius_biocellatus          0.63014  0.77648 0.0213 0.414585   
    ## Silhouettea_sp                   0.20811  0.97811 0.0024 1.000000   
    ## Sphaeramia_nematoptera           0.71643 -0.69766 0.1362 0.252747   
    ## Sphaeramia_orbicularis           0.63960 -0.76871 0.1605 0.180819   
    ## Sphyraena_barracuda              0.83899  0.54414 0.0921 0.433566   
    ## Sphyraena_obtusata               0.80335  0.59550 0.0117 0.951049   
    ## Sphyraena_qenie                  0.81538  0.57892 0.0155 0.688312   
    ## Stegastes_nigricans              0.98712  0.15999 0.0146 0.767233   
    ## Stegastes_punctatus              0.80693  0.59064 0.0365 0.652348   
    ## Stethojulis_strigiventer         0.95269  0.30395 0.0516 0.605395   
    ## Strongylura_incisa               0.71190  0.70228 0.0382 0.628372   
    ## Synodus_variegatus               0.31495 -0.94911 0.0673 0.416583   
    ## Taeniamia_zosterophora           0.35694 -0.93413 0.0896 0.315684   
    ## Taeniurops_meyeni                0.71129  0.70290 0.0134 0.825175   
    ## Toxotes_jaculatrix               0.87122  0.49090 0.1170 0.396603   
    ## Tylosurus_punctulatus            0.88861  0.45866 0.0283 0.893107   
    ## Upeneus_taeniopterus             0.81538  0.57892 0.0155 0.688312   
    ## Valenciennea_muralis             0.83346  0.55259 0.0488 0.653347   
    ## Valenciennea_parva               0.95129  0.30831 0.0155 0.724276   
    ## Valenciennea_randalli            0.95129  0.30831 0.0155 0.724276   
    ## Zanclus_cornutus                 0.99878  0.04928 0.0549 0.538462   
    ## Zebrasoma_scopas                 0.96122 -0.27578 0.0202 0.476523   
    ## Zebrasoma_velifer                0.78686 -0.61713 0.1816 0.133866   
    ## Zenarchopterus_dispar            0.98230 -0.18729 0.2631 0.058941 . 
    ## Zenarchopterus_gilli             0.99846 -0.05548 0.0867 0.448551   
    ## Zoramia_leptacanthus             0.98712  0.15999 0.0146 0.767233   
    ## Zoramia_viridiventer             0.51116 -0.85949 0.0923 0.352647   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# By stratification
fds_ef <- envfit(fd_NMDS, surveyed_sites_lake, permutations = 1000, na.rm = TRUE, strata = env[-20,19])
fds_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                                    NMDS1    NMDS2     r2  Pr(>r)  
    ## Abudefduf_lorenzi               -0.23611 -0.97173 0.0149 0.89411  
    ## Abudefduf_septemfasciatus        0.99529  0.09692 0.0169 0.34965  
    ## Abudefduf_sexfasciatus          -0.23611 -0.97173 0.0149 0.89411  
    ## Acanthurus_lineatus              0.63014  0.77648 0.0213 0.51648  
    ## Acanthurus_nigricauda            0.97617  0.21703 0.0556 0.53746  
    ## Acanthurus_sp                    0.90500 -0.42542 0.0199 0.83716  
    ## Acanthurus_xanthopterus          0.75007  0.66135 0.1882 0.24376  
    ## Acentrogobius_janthinopterus    -0.83826  0.54527 0.3818 0.35864  
    ## Amblyeleotris_gymnocephala       0.95129  0.30831 0.0155 0.62438  
    ## Amblyglyphidodon_curacao         0.60749 -0.79432 0.0382 0.84316  
    ## Amblygobius_buanensis            0.85468  0.51915 0.3321 0.38262  
    ## Amblygobius_decussatus           0.25879 -0.96593 0.0856 0.34965  
    ## Amblygobius_esakiae              0.90500 -0.42542 0.0199 0.83716  
    ## Amblygobius_linki                0.66600  0.74595 0.0187 0.27672  
    ## Amblygobius_phalaena             0.97146  0.23722 0.0315 0.92108  
    ## Amphiprion_clarkii               0.96122 -0.27578 0.0202 0.68232  
    ## Ancistrogobius_dipus             0.71129  0.70290 0.0134 0.77423  
    ## Apogonichthyoides_melas          0.71129  0.70290 0.0134 0.77423  
    ## Arcygobius_baliurus              0.81538  0.57892 0.0155 0.49151  
    ## Arothron_reticularis             0.74779  0.66393 0.0129 0.86214  
    ## Asterropteryx_semipunctata       0.99446  0.10514 0.0539 0.91209  
    ## Asterropteryx_sp                 0.66600  0.74595 0.0187 0.27672  
    ## Atherinomorus_lacunosus          0.97416  0.22586 0.0160 0.86813  
    ## Atule_mate                       0.63014  0.77648 0.0213 0.51648  
    ## Balistoides_viridescens          0.89129  0.45343 0.1280 0.86314  
    ## Bolbometopon_muricatum           0.36612 -0.93057 0.0885 0.51249  
    ## Caesio_caerulaurea               0.71452  0.69961 0.0341 0.27772  
    ## Caesio_cuning                    0.48565 -0.87415 0.1028 0.56144  
    ## Caesio_teres                     0.95129  0.30831 0.0155 0.62438  
    ## Cantherhines_dumerilii           0.96122 -0.27578 0.0202 0.68232  
    ## Canthigaster_bennetti            0.96122 -0.27578 0.0202 0.68232  
    ## Canthigaster_solandri            0.99529  0.09692 0.0169 0.34965  
    ## Caranx_ignobilis                 0.98712  0.15999 0.0146 1.00000  
    ## Caranx_melampygus                0.38841 -0.92149 0.0053 0.98501  
    ## Caranx_papuensis                 0.95129  0.30831 0.0155 0.62438  
    ## Caranx_sexfasciatus              0.80520  0.59300 0.0939 0.48751  
    ## Cephalopholis_argus             -0.44267 -0.89668 0.2031 0.50649  
    ## Cephalopholis_boenak             0.93543 -0.35352 0.0421 0.69730  
    ## Chaetodon_auriga                 0.96905 -0.24685 0.2867 0.84316  
    ## Chaetodon_bennetti               0.90818  0.41857 0.0416 0.73726  
    ## Chaetodon_ephippium              0.89897 -0.43801 0.2449 0.83117  
    ## Chaetodon_kleinii                0.92237  0.38631 0.0607 0.69930  
    ## Chaetodon_lunula                 0.87877  0.47725 0.1594 0.08891 .
    ## Chaetodon_lunulatus              0.75304 -0.65798 0.1443 0.64436  
    ## Chaetodon_melannotus             0.71723  0.69683 0.0193 0.12488  
    ## Chaetodon_ocellicaudus           0.98712  0.15999 0.0146 1.00000  
    ## Chaetodon_octofasciatus          0.93543 -0.35352 0.0421 0.69730  
    ## Chaetodon_rafflesii              0.98855  0.15092 0.0850 1.00000  
    ## Chaetodon_semeion                0.95129  0.30829 0.0401 0.40559  
    ## Chaetodon_trifascialis           0.96122 -0.27578 0.0202 0.68232  
    ## Chaetodon_ulietensis             0.51780 -0.85550 0.1000 0.46054  
    ## Chaetodon_vagabundus             0.97243 -0.23319 0.2893 0.43856  
    ## Chaetodontoplus_poliourus        0.96122 -0.27578 0.0202 0.68232  
    ## Cheilinus_fasciatus              0.98894  0.14830 0.0647 0.86014  
    ## Cheilinus_trilobatus             0.63014  0.77648 0.0213 0.51648  
    ## Cheilinus_undulatus              0.25412 -0.96717 0.0745 0.46753  
    ## Cheilodipterus_artus             0.21337 -0.97697 0.0100 0.86713  
    ## Cheilodipterus_isostigma         0.27332 -0.96192 0.0694 0.17283  
    ## Cheilodipterus_quinquelineatus   0.90937 -0.41599 0.2468 0.78821  
    ## Cheilodipterus_singapurensis     0.75449  0.65631 0.0320 0.55644  
    ## Chlorurus_bleekeri               0.51451 -0.85749 0.1093 0.67433  
    ## Chlorurus_microrhinos            0.99678 -0.08021 0.0361 0.93207  
    ## Chlorurus_spilurus               0.98855  0.15092 0.0850 1.00000  
    ## Choerodon_anchorago              0.99634 -0.08543 0.3920 0.78621  
    ## Chromis_viridis                  0.42127 -0.90693 0.0765 0.43656  
    ## Chrysiptera_biocellata           0.93034  0.36669 0.0338 0.31768  
    ## Chrysiptera_oxycephala           0.35041 -0.93660 0.0784 0.71828  
    ## Corythoichthys_ocellatus         0.96122 -0.27578 0.0202 0.68232  
    ## Cryptocentrus_caeruleomaculatus  0.67547 -0.73739 0.0692 0.43357  
    ## Cryptocentrus_cyanospilotus      0.20811  0.97811 0.0024 1.00000  
    ## Cryptocentrus_leptocephalus      0.91608  0.40100 0.0314 0.92807  
    ## Cryptocentrus_strigilliceps      0.51006 -0.86014 0.0988 0.51748  
    ## Ctenochaetus_striatus            0.16568 -0.98618 0.0700 0.70829  
    ## Ctenogobiops_pomastictus         0.71129  0.70290 0.0134 0.77423  
    ## Cymbacephalus_beauforti          0.98712  0.15999 0.0146 1.00000  
    ## Dascyllus_aruanus                0.59229 -0.80573 0.0848 0.83716  
    ## Dascyllus_trimaculatus           0.96122 -0.27578 0.0202 0.68232  
    ## Diplogrammus_goramensis          0.87289  0.48792 0.0290 0.96204  
    ## Diproctacanthus_xanthurus        0.90818  0.41857 0.0416 0.73726  
    ## Dischistodus_chrysopoecilus      0.88446  0.46662 0.0280 0.97802  
    ## Dischistodus_perspicillatus      0.80269 -0.59639 0.1279 0.89411  
    ## Doboatherina_duodecimalis        0.28051 -0.95985 0.1295 0.18681  
    ## Eleotris_fusca                  -0.30761  0.95151 0.3599 0.31069  
    ## Epibulus_brevis                  0.60937 -0.79288 0.1286 0.49451  
    ## Epinephelus_coeruleopunctatus    0.73534  0.67770 0.0358 0.21878  
    ## Epinephelus_merra                0.53707 -0.84354 0.1002 0.45055  
    ## Epinephelus_sp                   0.73002  0.68342 0.0125 0.87313  
    ## Eviota_atriventris               0.96020  0.27931 0.0527 0.73227  
    ## Eviota_bifasciata                0.82172 -0.56989 0.1761 0.47752  
    ## Eviota_fallax                    0.73395 -0.67920 0.1533 0.44555  
    ## Eviota_lachdeberei               0.85288  0.52211 0.0721 0.79221  
    ## Eviota_maculosa                  0.13341 -0.99106 0.0988 0.32967  
    ## Eviota_sigillata                 0.83456  0.55091 0.0361 0.18681  
    ## Eviota_sp                        0.66600  0.74595 0.0187 0.27672  
    ## Eviota_storthynx                 0.84710  0.53144 0.0539 0.66234  
    ## Exyrias_belissimus               0.91234 -0.40944 0.1514 0.95604  
    ## Exyrias_puntang                 -0.04985  0.99876 0.2867 0.19980  
    ## Favonigobius_reichei             0.91608  0.40100 0.0314 0.92807  
    ## Fibramia_lateralis               0.04556 -0.99896 0.0232 0.78921  
    ## Fibramia_thermalis               0.66627 -0.74571 0.0970 0.27872  
    ## Fusigobius_signipinnis           0.99878  0.04928 0.0549 0.90609  
    ## Gerres_erythrourus               0.98712  0.15999 0.0146 1.00000  
    ## Gerres_oyena                     0.99184  0.12747 0.0330 0.81019  
    ## Gladiogobius_ensifer             0.13341 -0.99106 0.0988 0.32967  
    ## Gnathanodon_speciosus            0.49547  0.86862 0.0130 0.97902  
    ## Gobiodon_okinawae                0.63014  0.77648 0.0213 0.51648  
    ## Gymnothorax_pictus               0.93034  0.36669 0.0338 0.31768  
    ## Halichoeres_chloropterus         0.77497 -0.63199 0.0536 0.85614  
    ## Halichoeres_leucurus             0.98248  0.18639 0.0862 0.79920  
    ## Halichoeres_marginatus           0.98712  0.15999 0.0146 1.00000  
    ## Halichoeres_melanurus            0.18559 -0.98263 0.0488 0.68931  
    ## Halichoeres_papilionaceus        0.80965  0.58692 0.0555 0.86014  
    ## Halichoeres_trimaculatus         0.99184  0.12747 0.0330 0.81019  
    ## Hemiglyphidodon_plagiometopon    0.45128 -0.89238 0.1069 0.70729  
    ## Hemigymnus_melapterus            0.99678 -0.08021 0.0361 0.93207  
    ## Heniochus_chrysostomus           0.93764  0.34761 0.0595 0.95005  
    ## Heniochus_monoceros              0.84710  0.53144 0.0539 0.66234  
    ## Herklotsichthys_quadrimaculatus  0.98712  0.15999 0.0146 1.00000  
    ## Hipposcarus_longiceps            0.44737 -0.89435 0.1113 0.55145  
    ## Istigobius_decoratus             0.98712  0.15999 0.0146 1.00000  
    ## Koumansetta_rainfordi            0.96122 -0.27578 0.0202 0.68232  
    ## Kyphosus_cinerascens            -0.21320 -0.97701 0.0174 0.89910  
    ## Lethrinus_erythropterus          0.88113  0.47287 0.0682 0.91608  
    ## Lethrinus_harak                  0.83501  0.55023 0.1182 0.80619  
    ## Lethrinus_obsoletus              0.71129  0.70290 0.0134 0.77423  
    ## Lethrinus_olivaceus              0.71129  0.70290 0.0134 0.77423  
    ## Lutjanus_argentimaculatus        0.42544  0.90499 0.1341 0.12388  
    ## Lutjanus_biguttatus             -0.25345 -0.96735 0.0169 0.86314  
    ## Lutjanus_bohar                   0.95129  0.30829 0.0401 0.40559  
    ## Lutjanus_decussatus              0.71129  0.70290 0.0134 0.77423  
    ## Lutjanus_ehrenbergii             0.78768  0.61609 0.0911 0.90310  
    ## Lutjanus_fulviflamma             0.87487  0.48436 0.0859 0.31868  
    ## Lutjanus_fulvus                  0.55440  0.83225 0.3632 0.02498 *
    ## Lutjanus_gibbus                  0.91817  0.39619 0.0762 0.85514  
    ## Lutjanus_monostigma              0.89150  0.45302 0.0285 0.89311  
    ## Lutjanus_russellii               0.66600  0.74595 0.0187 0.27672  
    ## Macrodontogobius_wilburi         0.75304 -0.65798 0.1443 0.64436  
    ## Megalops_cyprinoides            -0.17801  0.98403 0.5178 0.11189  
    ## Meiacanthus_ditrema              0.90500 -0.42542 0.0199 0.83716  
    ## Meiacanthus_grammistes           0.93543 -0.35352 0.0421 0.69730  
    ## Monodactylus_argenteus           0.68026  0.73297 0.0316 0.64535  
    ## Monotaxis_grandoculis           -0.02104 -0.99978 0.0121 0.81818  
    ## Mugilogobius_cavifrons          -0.60182  0.79863 0.2146 0.36364  
    ## Mulloidichthys_flavolineatus     0.98712  0.15999 0.0146 1.00000  
    ## Myripristis_adusta               0.81254  0.58291 0.1256 0.35165  
    ## Myripristis_sp                   0.81538  0.57892 0.0155 0.49151  
    ## Myripristis_violacea             0.96122 -0.27578 0.0202 0.68232  
    ## Naso_brevirostris                0.98494  0.17289 0.0776 0.82617  
    ## Naso_vlamingii                   0.87638  0.48163 0.0571 0.77522  
    ## Nectamia_viria                   0.71723  0.69683 0.0193 0.12488  
    ## Neoglyphidodon_melas            -0.44267 -0.89668 0.2031 0.50649  
    ## Neoniphon_argenteus              0.77483  0.63217 0.0508 0.48452  
    ## Neoniphon_opercularis           -0.44267 -0.89668 0.2031 0.50649  
    ## Neoniphon_sammara                0.95550  0.29498 0.0951 0.93307  
    ## Neopomacentrus_nemurus           0.25267 -0.96755 0.0634 0.60040  
    ## Omobranchus_obliquus             0.60665  0.79497 0.0333 0.54545  
    ## Oplopomops_diacanthus            0.98712  0.15999 0.0146 1.00000  
    ## Oxycheilinus_celebicus           0.22029 -0.97544 0.0156 0.82218  
    ## Parapercis_cylindrica            0.98712  0.15999 0.0146 1.00000  
    ## Parioglossus_formosus           -0.10980  0.99395 0.2528 0.05395 .
    ## Parioglossus_philippinus         0.75917  0.65089 0.0543 0.63836  
    ## Parioglossus_sp                  0.98511  0.17192 0.0366 0.59441  
    ## Parupeneus_barberinus            0.96515  0.26170 0.1647 0.81419  
    ## Parupeneus_multifasciatus        0.88389 -0.46769 0.2084 0.47353  
    ## Pentapodus_trivittatus           0.87597 -0.48237 0.2074 0.82517  
    ## Periophthalmus_argentilineatus   0.99529  0.09692 0.0169 0.34965  
    ## Pervagor_janthinosoma            0.63014  0.77648 0.0213 0.51648  
    ## Pervagor_nigrolineatus           0.37829 -0.92569 0.0786 0.37463  
    ## Petroscirtes_breviceps           0.99811  0.06153 0.1038 0.88811  
    ## Platax_orbicularis               0.76465  0.64444 0.0303 0.65534  
    ## Platax_pinnatus                  0.95129  0.30831 0.0155 0.62438  
    ## Platax_sp                        0.73002  0.68342 0.0125 0.87313  
    ## Plectorhinchus_albovittatus      0.66600  0.74595 0.0187 0.27672  
    ## Plectorhinchus_chaetodonoides    0.13341 -0.99106 0.0988 0.32967  
    ## Plectorhinchus_gibbosus          0.98712  0.15999 0.0146 1.00000  
    ## Plectropomus_leopardus           0.96122 -0.27578 0.0202 0.68232  
    ## Pomacanthus_sexstriatus          0.86214  0.50667 0.0349 0.71029  
    ## Pomacentrus_adelus               0.99529  0.09692 0.0169 0.34965  
    ## Pomacentrus_amboinensis          0.96122 -0.27578 0.0202 0.68232  
    ## Pomacentrus_burroughi            0.84161 -0.54009 0.1823 0.42258  
    ## Pomacentrus_grammorhynchus       0.98712  0.15999 0.0146 1.00000  
    ## Pomacentrus_nigromanus          -0.15416 -0.98805 0.0248 0.77822  
    ## Pomacentrus_pavo                 0.96122 -0.27578 0.0202 0.68232  
    ## Pomacentrus_simsiang             0.92785 -0.37296 0.2872 0.82318  
    ## Pomacentrus_sp                  -0.32959 -0.94412 0.0637 0.57343  
    ## Pomacentrus_taeniometopon        0.66600  0.74595 0.0187 0.27672  
    ## Pseudobalistes_flavimarginatus   0.95775  0.28760 0.1323 0.94605  
    ## Pseudochromis_fuscus             0.96122 -0.27578 0.0202 0.68232  
    ## Ptereleotris_brachyptera         1.00000 -0.00184 0.0370 0.61039  
    ## Ptereleotris_evides              0.96122 -0.27578 0.0202 0.68232  
    ## Pterocaesio_tile                 0.71129  0.70290 0.0134 0.77423  
    ## Pterocaesio_trilineata           0.71723  0.69683 0.0193 0.12488  
    ## Pterois_volitans                 0.95129  0.30831 0.0155 0.62438  
    ## Rastrelliger_kanagurta           0.31495 -0.94911 0.0673 0.60739  
    ## Rhabdamia_gracilis               0.83309  0.55314 0.0727 0.40559  
    ## Rhinecanthus_aculeatus           0.88112  0.47290 0.0607 0.67932  
    ## Salarias_segmentatus             0.90427  0.42696 0.1687 0.34765  
    ## Sargocentron_diadema             0.67145  0.74105 0.0425 0.35764  
    ## Sargocentron_spiniferum          0.81538  0.57892 0.0155 0.49151  
    ## Saurida_gracilis                 0.63014  0.77648 0.0213 0.51648  
    ## Scarus_dimidiatus                0.98855  0.15092 0.0850 1.00000  
    ## Scarus_flavipectoralis           0.90818  0.41857 0.0416 0.73726  
    ## Scarus_ghobban                   0.98511  0.17192 0.0366 0.59441  
    ## Scarus_hypselopterus             0.63014  0.77648 0.0213 0.51648  
    ## Scarus_oviceps                   0.90500 -0.42542 0.0199 0.83716  
    ## Scarus_psittacus                 0.98112  0.19339 0.1072 0.78022  
    ## Scarus_quoyi                     0.96122 -0.27578 0.0202 0.68232  
    ## Scarus_rivulatus                 0.25879 -0.96593 0.0856 0.34965  
    ## Scatophagus_argus                0.74779  0.66393 0.0129 0.86214  
    ## Scolopsis_ciliata                0.73156 -0.68177 0.1706 0.62937  
    ## Scolopsis_lineata                0.91608  0.40100 0.0314 0.92807  
    ## Scolopsis_margaritifera          0.35818 -0.93365 0.0807 0.69830  
    ## Scolopsis_trilineata             0.74712  0.66469 0.0530 0.27173  
    ## Siganus_doliatus                 0.36612 -0.93057 0.0885 0.51249  
    ## Siganus_fuscescens               0.77733  0.62909 0.0530 0.93007  
    ## Siganus_lineatus                 0.58173  0.81338 0.2995 0.02298 *
    ## Siganus_puellus                  0.35551 -0.93467 0.0998 0.30270  
    ## Siganus_punctatissimus           0.96122 -0.27578 0.0202 0.68232  
    ## Siganus_punctatus                0.96122 -0.27578 0.0202 0.68232  
    ## Signigobius_biocellatus          0.63014  0.77648 0.0213 0.51648  
    ## Silhouettea_sp                   0.20811  0.97811 0.0024 1.00000  
    ## Sphaeramia_nematoptera           0.71643 -0.69766 0.1362 0.88412  
    ## Sphaeramia_orbicularis           0.63960 -0.76871 0.1605 0.20779  
    ## Sphyraena_barracuda              0.83899  0.54414 0.0921 0.86713  
    ## Sphyraena_obtusata               0.80335  0.59550 0.0117 1.00000  
    ## Sphyraena_qenie                  0.81538  0.57892 0.0155 0.49151  
    ## Stegastes_nigricans              0.98712  0.15999 0.0146 1.00000  
    ## Stegastes_punctatus              0.80693  0.59064 0.0365 0.87213  
    ## Stethojulis_strigiventer         0.95269  0.30395 0.0516 0.76723  
    ## Strongylura_incisa               0.71190  0.70228 0.0382 0.53746  
    ## Synodus_variegatus               0.31495 -0.94911 0.0673 0.60739  
    ## Taeniamia_zosterophora           0.35694 -0.93413 0.0896 0.47053  
    ## Taeniurops_meyeni                0.71129  0.70290 0.0134 0.77423  
    ## Toxotes_jaculatrix               0.87122  0.49090 0.1170 0.85714  
    ## Tylosurus_punctulatus            0.88861  0.45866 0.0283 0.89111  
    ## Upeneus_taeniopterus             0.81538  0.57892 0.0155 0.49151  
    ## Valenciennea_muralis             0.83346  0.55259 0.0488 0.66833  
    ## Valenciennea_parva               0.95129  0.30831 0.0155 0.62438  
    ## Valenciennea_randalli            0.95129  0.30831 0.0155 0.62438  
    ## Zanclus_cornutus                 0.99878  0.04928 0.0549 0.90609  
    ## Zebrasoma_scopas                 0.96122 -0.27578 0.0202 0.68232  
    ## Zebrasoma_velifer                0.78686 -0.61713 0.1816 0.73127  
    ## Zenarchopterus_dispar            0.98230 -0.18729 0.2631 0.97602  
    ## Zenarchopterus_gilli             0.99846 -0.05548 0.0867 0.93107  
    ## Zoramia_leptacanthus             0.98712  0.15999 0.0146 1.00000  
    ## Zoramia_viridiventer             0.51116 -0.85949 0.0923 0.31668  
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
    ## Abudefduf_lorenzi               -0.23611 -0.97173 0.0149      1
    ## Abudefduf_septemfasciatus        0.99529  0.09692 0.0169      1
    ## Abudefduf_sexfasciatus          -0.23611 -0.97173 0.0149      1
    ## Acanthurus_lineatus              0.63014  0.77648 0.0213      1
    ## Acanthurus_nigricauda            0.97617  0.21703 0.0556      1
    ## Acanthurus_sp                    0.90500 -0.42542 0.0199      1
    ## Acanthurus_xanthopterus          0.75007  0.66135 0.1882      1
    ## Acentrogobius_janthinopterus    -0.83826  0.54527 0.3818      1
    ## Amblyeleotris_gymnocephala       0.95129  0.30831 0.0155      1
    ## Amblyglyphidodon_curacao         0.60749 -0.79432 0.0382      1
    ## Amblygobius_buanensis            0.85468  0.51915 0.3321      1
    ## Amblygobius_decussatus           0.25879 -0.96593 0.0856      1
    ## Amblygobius_esakiae              0.90500 -0.42542 0.0199      1
    ## Amblygobius_linki                0.66600  0.74595 0.0187      1
    ## Amblygobius_phalaena             0.97146  0.23722 0.0315      1
    ## Amphiprion_clarkii               0.96122 -0.27578 0.0202      1
    ## Ancistrogobius_dipus             0.71129  0.70290 0.0134      1
    ## Apogonichthyoides_melas          0.71129  0.70290 0.0134      1
    ## Arcygobius_baliurus              0.81538  0.57892 0.0155      1
    ## Arothron_reticularis             0.74779  0.66393 0.0129      1
    ## Asterropteryx_semipunctata       0.99446  0.10514 0.0539      1
    ## Asterropteryx_sp                 0.66600  0.74595 0.0187      1
    ## Atherinomorus_lacunosus          0.97416  0.22586 0.0160      1
    ## Atule_mate                       0.63014  0.77648 0.0213      1
    ## Balistoides_viridescens          0.89129  0.45343 0.1280      1
    ## Bolbometopon_muricatum           0.36612 -0.93057 0.0885      1
    ## Caesio_caerulaurea               0.71452  0.69961 0.0341      1
    ## Caesio_cuning                    0.48565 -0.87415 0.1028      1
    ## Caesio_teres                     0.95129  0.30831 0.0155      1
    ## Cantherhines_dumerilii           0.96122 -0.27578 0.0202      1
    ## Canthigaster_bennetti            0.96122 -0.27578 0.0202      1
    ## Canthigaster_solandri            0.99529  0.09692 0.0169      1
    ## Caranx_ignobilis                 0.98712  0.15999 0.0146      1
    ## Caranx_melampygus                0.38841 -0.92149 0.0053      1
    ## Caranx_papuensis                 0.95129  0.30831 0.0155      1
    ## Caranx_sexfasciatus              0.80520  0.59300 0.0939      1
    ## Cephalopholis_argus             -0.44267 -0.89668 0.2031      1
    ## Cephalopholis_boenak             0.93543 -0.35352 0.0421      1
    ## Chaetodon_auriga                 0.96905 -0.24685 0.2867      1
    ## Chaetodon_bennetti               0.90818  0.41857 0.0416      1
    ## Chaetodon_ephippium              0.89897 -0.43801 0.2449      1
    ## Chaetodon_kleinii                0.92237  0.38631 0.0607      1
    ## Chaetodon_lunula                 0.87877  0.47725 0.1594      1
    ## Chaetodon_lunulatus              0.75304 -0.65798 0.1443      1
    ## Chaetodon_melannotus             0.71723  0.69683 0.0193      1
    ## Chaetodon_ocellicaudus           0.98712  0.15999 0.0146      1
    ## Chaetodon_octofasciatus          0.93543 -0.35352 0.0421      1
    ## Chaetodon_rafflesii              0.98855  0.15092 0.0850      1
    ## Chaetodon_semeion                0.95129  0.30829 0.0401      1
    ## Chaetodon_trifascialis           0.96122 -0.27578 0.0202      1
    ## Chaetodon_ulietensis             0.51780 -0.85550 0.1000      1
    ## Chaetodon_vagabundus             0.97243 -0.23319 0.2893      1
    ## Chaetodontoplus_poliourus        0.96122 -0.27578 0.0202      1
    ## Cheilinus_fasciatus              0.98894  0.14830 0.0647      1
    ## Cheilinus_trilobatus             0.63014  0.77648 0.0213      1
    ## Cheilinus_undulatus              0.25412 -0.96717 0.0745      1
    ## Cheilodipterus_artus             0.21337 -0.97697 0.0100      1
    ## Cheilodipterus_isostigma         0.27332 -0.96192 0.0694      1
    ## Cheilodipterus_quinquelineatus   0.90937 -0.41599 0.2468      1
    ## Cheilodipterus_singapurensis     0.75449  0.65631 0.0320      1
    ## Chlorurus_bleekeri               0.51451 -0.85749 0.1093      1
    ## Chlorurus_microrhinos            0.99678 -0.08021 0.0361      1
    ## Chlorurus_spilurus               0.98855  0.15092 0.0850      1
    ## Choerodon_anchorago              0.99634 -0.08543 0.3920      1
    ## Chromis_viridis                  0.42127 -0.90693 0.0765      1
    ## Chrysiptera_biocellata           0.93034  0.36669 0.0338      1
    ## Chrysiptera_oxycephala           0.35041 -0.93660 0.0784      1
    ## Corythoichthys_ocellatus         0.96122 -0.27578 0.0202      1
    ## Cryptocentrus_caeruleomaculatus  0.67547 -0.73739 0.0692      1
    ## Cryptocentrus_cyanospilotus      0.20811  0.97811 0.0024      1
    ## Cryptocentrus_leptocephalus      0.91608  0.40100 0.0314      1
    ## Cryptocentrus_strigilliceps      0.51006 -0.86014 0.0988      1
    ## Ctenochaetus_striatus            0.16568 -0.98618 0.0700      1
    ## Ctenogobiops_pomastictus         0.71129  0.70290 0.0134      1
    ## Cymbacephalus_beauforti          0.98712  0.15999 0.0146      1
    ## Dascyllus_aruanus                0.59229 -0.80573 0.0848      1
    ## Dascyllus_trimaculatus           0.96122 -0.27578 0.0202      1
    ## Diplogrammus_goramensis          0.87289  0.48792 0.0290      1
    ## Diproctacanthus_xanthurus        0.90818  0.41857 0.0416      1
    ## Dischistodus_chrysopoecilus      0.88446  0.46662 0.0280      1
    ## Dischistodus_perspicillatus      0.80269 -0.59639 0.1279      1
    ## Doboatherina_duodecimalis        0.28051 -0.95985 0.1295      1
    ## Eleotris_fusca                  -0.30761  0.95151 0.3599      1
    ## Epibulus_brevis                  0.60937 -0.79288 0.1286      1
    ## Epinephelus_coeruleopunctatus    0.73534  0.67770 0.0358      1
    ## Epinephelus_merra                0.53707 -0.84354 0.1002      1
    ## Epinephelus_sp                   0.73002  0.68342 0.0125      1
    ## Eviota_atriventris               0.96020  0.27931 0.0527      1
    ## Eviota_bifasciata                0.82172 -0.56989 0.1761      1
    ## Eviota_fallax                    0.73395 -0.67920 0.1533      1
    ## Eviota_lachdeberei               0.85288  0.52211 0.0721      1
    ## Eviota_maculosa                  0.13341 -0.99106 0.0988      1
    ## Eviota_sigillata                 0.83456  0.55091 0.0361      1
    ## Eviota_sp                        0.66600  0.74595 0.0187      1
    ## Eviota_storthynx                 0.84710  0.53144 0.0539      1
    ## Exyrias_belissimus               0.91234 -0.40944 0.1514      1
    ## Exyrias_puntang                 -0.04985  0.99876 0.2867      1
    ## Favonigobius_reichei             0.91608  0.40100 0.0314      1
    ## Fibramia_lateralis               0.04556 -0.99896 0.0232      1
    ## Fibramia_thermalis               0.66627 -0.74571 0.0970      1
    ## Fusigobius_signipinnis           0.99878  0.04928 0.0549      1
    ## Gerres_erythrourus               0.98712  0.15999 0.0146      1
    ## Gerres_oyena                     0.99184  0.12747 0.0330      1
    ## Gladiogobius_ensifer             0.13341 -0.99106 0.0988      1
    ## Gnathanodon_speciosus            0.49547  0.86862 0.0130      1
    ## Gobiodon_okinawae                0.63014  0.77648 0.0213      1
    ## Gymnothorax_pictus               0.93034  0.36669 0.0338      1
    ## Halichoeres_chloropterus         0.77497 -0.63199 0.0536      1
    ## Halichoeres_leucurus             0.98248  0.18639 0.0862      1
    ## Halichoeres_marginatus           0.98712  0.15999 0.0146      1
    ## Halichoeres_melanurus            0.18559 -0.98263 0.0488      1
    ## Halichoeres_papilionaceus        0.80965  0.58692 0.0555      1
    ## Halichoeres_trimaculatus         0.99184  0.12747 0.0330      1
    ## Hemiglyphidodon_plagiometopon    0.45128 -0.89238 0.1069      1
    ## Hemigymnus_melapterus            0.99678 -0.08021 0.0361      1
    ## Heniochus_chrysostomus           0.93764  0.34761 0.0595      1
    ## Heniochus_monoceros              0.84710  0.53144 0.0539      1
    ## Herklotsichthys_quadrimaculatus  0.98712  0.15999 0.0146      1
    ## Hipposcarus_longiceps            0.44737 -0.89435 0.1113      1
    ## Istigobius_decoratus             0.98712  0.15999 0.0146      1
    ## Koumansetta_rainfordi            0.96122 -0.27578 0.0202      1
    ## Kyphosus_cinerascens            -0.21320 -0.97701 0.0174      1
    ## Lethrinus_erythropterus          0.88113  0.47287 0.0682      1
    ## Lethrinus_harak                  0.83501  0.55023 0.1182      1
    ## Lethrinus_obsoletus              0.71129  0.70290 0.0134      1
    ## Lethrinus_olivaceus              0.71129  0.70290 0.0134      1
    ## Lutjanus_argentimaculatus        0.42544  0.90499 0.1341      1
    ## Lutjanus_biguttatus             -0.25345 -0.96735 0.0169      1
    ## Lutjanus_bohar                   0.95129  0.30829 0.0401      1
    ## Lutjanus_decussatus              0.71129  0.70290 0.0134      1
    ## Lutjanus_ehrenbergii             0.78768  0.61609 0.0911      1
    ## Lutjanus_fulviflamma             0.87487  0.48436 0.0859      1
    ## Lutjanus_fulvus                  0.55440  0.83225 0.3632      1
    ## Lutjanus_gibbus                  0.91817  0.39619 0.0762      1
    ## Lutjanus_monostigma              0.89150  0.45302 0.0285      1
    ## Lutjanus_russellii               0.66600  0.74595 0.0187      1
    ## Macrodontogobius_wilburi         0.75304 -0.65798 0.1443      1
    ## Megalops_cyprinoides            -0.17801  0.98403 0.5178      1
    ## Meiacanthus_ditrema              0.90500 -0.42542 0.0199      1
    ## Meiacanthus_grammistes           0.93543 -0.35352 0.0421      1
    ## Monodactylus_argenteus           0.68026  0.73297 0.0316      1
    ## Monotaxis_grandoculis           -0.02104 -0.99978 0.0121      1
    ## Mugilogobius_cavifrons          -0.60182  0.79863 0.2146      1
    ## Mulloidichthys_flavolineatus     0.98712  0.15999 0.0146      1
    ## Myripristis_adusta               0.81254  0.58291 0.1256      1
    ## Myripristis_sp                   0.81538  0.57892 0.0155      1
    ## Myripristis_violacea             0.96122 -0.27578 0.0202      1
    ## Naso_brevirostris                0.98494  0.17289 0.0776      1
    ## Naso_vlamingii                   0.87638  0.48163 0.0571      1
    ## Nectamia_viria                   0.71723  0.69683 0.0193      1
    ## Neoglyphidodon_melas            -0.44267 -0.89668 0.2031      1
    ## Neoniphon_argenteus              0.77483  0.63217 0.0508      1
    ## Neoniphon_opercularis           -0.44267 -0.89668 0.2031      1
    ## Neoniphon_sammara                0.95550  0.29498 0.0951      1
    ## Neopomacentrus_nemurus           0.25267 -0.96755 0.0634      1
    ## Omobranchus_obliquus             0.60665  0.79497 0.0333      1
    ## Oplopomops_diacanthus            0.98712  0.15999 0.0146      1
    ## Oxycheilinus_celebicus           0.22029 -0.97544 0.0156      1
    ## Parapercis_cylindrica            0.98712  0.15999 0.0146      1
    ## Parioglossus_formosus           -0.10980  0.99395 0.2528      1
    ## Parioglossus_philippinus         0.75917  0.65089 0.0543      1
    ## Parioglossus_sp                  0.98511  0.17192 0.0366      1
    ## Parupeneus_barberinus            0.96515  0.26170 0.1647      1
    ## Parupeneus_multifasciatus        0.88389 -0.46769 0.2084      1
    ## Pentapodus_trivittatus           0.87597 -0.48237 0.2074      1
    ## Periophthalmus_argentilineatus   0.99529  0.09692 0.0169      1
    ## Pervagor_janthinosoma            0.63014  0.77648 0.0213      1
    ## Pervagor_nigrolineatus           0.37829 -0.92569 0.0786      1
    ## Petroscirtes_breviceps           0.99811  0.06153 0.1038      1
    ## Platax_orbicularis               0.76465  0.64444 0.0303      1
    ## Platax_pinnatus                  0.95129  0.30831 0.0155      1
    ## Platax_sp                        0.73002  0.68342 0.0125      1
    ## Plectorhinchus_albovittatus      0.66600  0.74595 0.0187      1
    ## Plectorhinchus_chaetodonoides    0.13341 -0.99106 0.0988      1
    ## Plectorhinchus_gibbosus          0.98712  0.15999 0.0146      1
    ## Plectropomus_leopardus           0.96122 -0.27578 0.0202      1
    ## Pomacanthus_sexstriatus          0.86214  0.50667 0.0349      1
    ## Pomacentrus_adelus               0.99529  0.09692 0.0169      1
    ## Pomacentrus_amboinensis          0.96122 -0.27578 0.0202      1
    ## Pomacentrus_burroughi            0.84161 -0.54009 0.1823      1
    ## Pomacentrus_grammorhynchus       0.98712  0.15999 0.0146      1
    ## Pomacentrus_nigromanus          -0.15416 -0.98805 0.0248      1
    ## Pomacentrus_pavo                 0.96122 -0.27578 0.0202      1
    ## Pomacentrus_simsiang             0.92785 -0.37296 0.2872      1
    ## Pomacentrus_sp                  -0.32959 -0.94412 0.0637      1
    ## Pomacentrus_taeniometopon        0.66600  0.74595 0.0187      1
    ## Pseudobalistes_flavimarginatus   0.95775  0.28760 0.1323      1
    ## Pseudochromis_fuscus             0.96122 -0.27578 0.0202      1
    ## Ptereleotris_brachyptera         1.00000 -0.00184 0.0370      1
    ## Ptereleotris_evides              0.96122 -0.27578 0.0202      1
    ## Pterocaesio_tile                 0.71129  0.70290 0.0134      1
    ## Pterocaesio_trilineata           0.71723  0.69683 0.0193      1
    ## Pterois_volitans                 0.95129  0.30831 0.0155      1
    ## Rastrelliger_kanagurta           0.31495 -0.94911 0.0673      1
    ## Rhabdamia_gracilis               0.83309  0.55314 0.0727      1
    ## Rhinecanthus_aculeatus           0.88112  0.47290 0.0607      1
    ## Salarias_segmentatus             0.90427  0.42696 0.1687      1
    ## Sargocentron_diadema             0.67145  0.74105 0.0425      1
    ## Sargocentron_spiniferum          0.81538  0.57892 0.0155      1
    ## Saurida_gracilis                 0.63014  0.77648 0.0213      1
    ## Scarus_dimidiatus                0.98855  0.15092 0.0850      1
    ## Scarus_flavipectoralis           0.90818  0.41857 0.0416      1
    ## Scarus_ghobban                   0.98511  0.17192 0.0366      1
    ## Scarus_hypselopterus             0.63014  0.77648 0.0213      1
    ## Scarus_oviceps                   0.90500 -0.42542 0.0199      1
    ## Scarus_psittacus                 0.98112  0.19339 0.1072      1
    ## Scarus_quoyi                     0.96122 -0.27578 0.0202      1
    ## Scarus_rivulatus                 0.25879 -0.96593 0.0856      1
    ## Scatophagus_argus                0.74779  0.66393 0.0129      1
    ## Scolopsis_ciliata                0.73156 -0.68177 0.1706      1
    ## Scolopsis_lineata                0.91608  0.40100 0.0314      1
    ## Scolopsis_margaritifera          0.35818 -0.93365 0.0807      1
    ## Scolopsis_trilineata             0.74712  0.66469 0.0530      1
    ## Siganus_doliatus                 0.36612 -0.93057 0.0885      1
    ## Siganus_fuscescens               0.77733  0.62909 0.0530      1
    ## Siganus_lineatus                 0.58173  0.81338 0.2995      1
    ## Siganus_puellus                  0.35551 -0.93467 0.0998      1
    ## Siganus_punctatissimus           0.96122 -0.27578 0.0202      1
    ## Siganus_punctatus                0.96122 -0.27578 0.0202      1
    ## Signigobius_biocellatus          0.63014  0.77648 0.0213      1
    ## Silhouettea_sp                   0.20811  0.97811 0.0024      1
    ## Sphaeramia_nematoptera           0.71643 -0.69766 0.1362      1
    ## Sphaeramia_orbicularis           0.63960 -0.76871 0.1605      1
    ## Sphyraena_barracuda              0.83899  0.54414 0.0921      1
    ## Sphyraena_obtusata               0.80335  0.59550 0.0117      1
    ## Sphyraena_qenie                  0.81538  0.57892 0.0155      1
    ## Stegastes_nigricans              0.98712  0.15999 0.0146      1
    ## Stegastes_punctatus              0.80693  0.59064 0.0365      1
    ## Stethojulis_strigiventer         0.95269  0.30395 0.0516      1
    ## Strongylura_incisa               0.71190  0.70228 0.0382      1
    ## Synodus_variegatus               0.31495 -0.94911 0.0673      1
    ## Taeniamia_zosterophora           0.35694 -0.93413 0.0896      1
    ## Taeniurops_meyeni                0.71129  0.70290 0.0134      1
    ## Toxotes_jaculatrix               0.87122  0.49090 0.1170      1
    ## Tylosurus_punctulatus            0.88861  0.45866 0.0283      1
    ## Upeneus_taeniopterus             0.81538  0.57892 0.0155      1
    ## Valenciennea_muralis             0.83346  0.55259 0.0488      1
    ## Valenciennea_parva               0.95129  0.30831 0.0155      1
    ## Valenciennea_randalli            0.95129  0.30831 0.0155      1
    ## Zanclus_cornutus                 0.99878  0.04928 0.0549      1
    ## Zebrasoma_scopas                 0.96122 -0.27578 0.0202      1
    ## Zebrasoma_velifer                0.78686 -0.61713 0.1816      1
    ## Zenarchopterus_dispar            0.98230 -0.18729 0.2631      1
    ## Zenarchopterus_gilli             0.99846 -0.05548 0.0867      1
    ## Zoramia_leptacanthus             0.98712  0.15999 0.0146      1
    ## Zoramia_viridiventer             0.51116 -0.85949 0.0923      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

## Trait mantel tests

We used mantel tests to determine significance between the trait and
environmental distance matrices.

``` r
### Environmental
# All sites
enve_dist <- dist(scaled_env[-c(9,16,18,20),c(2)], method = "euclidean")
fdea_mant <- mantel(fde_dist, enve_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
fdea_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fde_dist, ydis = enve_dist, method = "spearman",      permutations = 999, strata = env[-c(9, 16, 18, 20), 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.8057 
    ##       Significance: 0.015 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.702 0.743 0.774 0.816 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
envems_dist <- dist(scaled_env[c(1:6,8,10,12:15,17,21:23),c(2)], method = "euclidean")
fdems_mant <- mantel(fdems_dist, envems_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
fdems_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdems_dist, ydis = envems_dist, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.8214 
    ##       Significance: 0.011 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.684 0.745 0.793 0.817 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
enveom_dist <- dist(scaled_env[c(3,6:8,10:11,13:15,19,23),c(2:4)], method = "euclidean")
fdeom_mant <- mantel(fdeom_dist, enveom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
fdeom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdeom_dist, ydis = enveom_dist, method = "spearman",      permutations = 999, strata = env[c(3, 6:8, 10:11, 13:15,          19, 23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1356 
    ##       Significance: 0.142 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.170 0.227 0.275 0.369 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
enveso_dist <- dist(scaled_env[c(1:2,4,5,7,11:12,17,19,21:22),c(2)], method = "euclidean")
fdeso_mant <- mantel(fdeso_dist, enveso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,11:12,17,19,21:22),19])
fdeso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdeso_dist, ydis = enveso_dist, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 11:12, 17,          19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6976 
    ##       Significance: 0.01 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.462 0.527 0.593 0.670 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes
envs_dist <- dist(scaled_env[c(1:2,4,5,12,17,21:22),c(1:15)], method = "euclidean")
fds_mant <- mantel(fds_dist, envs_dist, method = "spearman", permutations = 999, na.rm = TRUE)
fds_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fds_dist, ydis = envs_dist, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.03996 
    ##       Significance: 0.409 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.258 0.348 0.415 0.474 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Biogeographic
# All sites
envb_dist_la <- dist(scaled_env[-c(8,20),c(15)], method = "euclidean")
fdba_mant_la <- mantel(fdb_dist, envb_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(8,20),19])
fdba_mant_la
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdb_dist, ydis = envb_dist_la, method = "spearman",      permutations = 999, strata = env[-c(8, 20), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1117 
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0426 0.0669 0.0821 0.1143 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes 
envbms_dist <- dist(scaled_env[c(1:6,10,12:15,17,21:23),c(5:15)], method = "euclidean")
fdbms_mant <- mantel(fdbms_dist, envbms_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,10,12:15,17,21:23),19])
fdbms_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbms_dist, ydis = envbms_dist, method = "spearman",      permutations = 999, strata = env[c(1:6, 10, 12:15, 17, 21:23),          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2809 
    ##       Significance: 0.283 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.417 0.473 0.526 0.578 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
envbom_dist <- dist(scaled_env[c(3,6:7,9:11,13:16,18:19,23),c(5:15)], method = "euclidean")
fdbom_mant <- mantel(fdbom_dist, envbom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:7,9:11,13:16,18:19,23),19])
fdbom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbom_dist, ydis = envbom_dist, method = "spearman",      permutations = 999, strata = env[c(3, 6:7, 9:11, 13:16, 18:19,          23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1645 
    ##       Significance: 0.603 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.487 0.582 0.617 0.669 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
envbso_dist <- dist(scaled_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(5:15)], method = "euclidean")
fdbso_mant <- mantel(fdbso_dist, envbso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
fdbso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = fdbso_dist, ydis = envbso_dist, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1656 
    ##       Significance: 0.317 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.242 0.270 0.295 0.318 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

## Plot trait NMDS and envfit results

Includes a plot of correlated environmental variables

``` r
fd_NMDS_data.scores <- as.data.frame(scores(fd_NMDS))
fd_NMDS_data.scores$Stratification <- env[-20,19]
fd_NMDS_data.scores$Lakes <- env[-20,1]
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
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), level = 0.50) +
  geom_point(data = fd_NMDS_data.scores, aes(color = Stratification), size = 4, alpha = 1) + 
  geom_text_repel(data = fd_NMDS_data.scores, label = fd_NMDS_data.scores$Lakes, size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = fde_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#CD661D") +
  geom_text(data = fde_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
            color = "#CD661D", label = row.names(fde_ef_coord_cont), size = 7) + 
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              data = fdb_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#698B22") +
  # geom_text(data = fdb_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
  #           color = "#698B22", label = row.names(fdb_ef_coord_cont), size = 7) + 
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
        panel.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  # annotate(geom = "label", x = 0.19, y = -0.3, size = 6,
  #          label = paste("Stress: ", round(fd_NMDS$stress, digits = 2))) +
  labs(color = "Stratification")
fd_ef_plot <- fd_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
fd_ef_plot
```

    ## Warning: ggrepel: 14 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20trait%20NMDS%20and%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_ef_plot.png", fd_ef_plot, width = 8, height = 4, units = "in")
```

    ## Warning: ggrepel: 14 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

## Boxplot of trait intradissimilarity heterogeneity

The dissimilarity between sites that share the same stratification.

``` r
fd_bd_dist <- stratification_fd_bd$distances
fd_bd_dist <- as.data.frame(fd_bd_dist)
fd_bd_dist$X <- row.names(fd_bd_dist)
fd_bd_dist_env <- merge(fd_bd_dist, env[-20,], by = "X", sort = F)

fd_bd_dist_env_plot <- ggplot(fd_bd_dist_env, aes(x = Stratification, y = fd_bd_dist, 
                                                  color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.75,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black")) +
  ylab("Distance to Centroid") +
  labs(color = "Stratification", tag = "A")
fd_bd_dist_env_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Boxplot%20of%20trait%20intradissimilarity%20heterogeneity-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/fd_bd_dist_env_plot.png", fd_bd_dist_env_plot, width = 8.25, height = 4.13, units = "in")
```

## Phylogeny dissimilarity distances

We use unifrac to calculate distances in the phylogeny.

``` r
### Regular
pd_dist <- unifrac(presabs_lake[-20,], stree)
# Mixed and stratified lakes
pdms_dist <- unifrac(presabs_lake[c(1:6,8,10,12:15,17,21:23),], stree)
# Ocean sites and mixed lakes
pdom_dist <- unifrac(presabs_lake[c(3,6:11,13:16,18:19,23),], stree)
# Stratified lakes and ocean sites
pdso_dist <- unifrac(presabs_lake[c(1:2,4,5,7,9,11:12,16:19,21:22),], stree)
# Stratified lakes
pds_dist <- unifrac(presabs_lake[c(1:2,4,5,12,17,21:22),], stree)

### Environmental
pde_dist <- unifrac(presabs_lake[-c(9,16,18,20),], etree)
# Mixed and stratified lakes
pdems_dist <- unifrac(presabs_lake[c(1:6,8,10,12:15,17,21:23),], etree)
# Ocean sites and mixed lakes
pdeom_dist <- unifrac(presabs_lake[c(3,6:8,10:11,13:15,19,23),], etree)
# Stratified lakes and ocean sites
pdeso_dist <- unifrac(presabs_lake[c(1:2,4,5,7,11:12,17,19,21:22),], etree)


### Biogeographic
pdb_dist <- unifrac(presabs_lake[-c(8,20),], btree)
# Mixed and stratified lakes
pdbms_dist <- unifrac(presabs_lake[c(1:6,10,12:15,17,21:23),], btree)
# Ocean sites and mixed lakes
pdbom_dist <- unifrac(presabs_lake[c(3,6:7,9:11,13:16,18:19,23),], btree)
# Stratified lakes and ocean sites
pdbso_dist <- unifrac(presabs_lake[c(1:2,4,5,7,9,11:12,16:19,21:22),], btree)
```

## Phylogeny NMDS

To constrain dissimilarities we perform Nonmetric Multidimensional
Scaling (NMDS), which tries to find a stable solution using the metaMDS
package.

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
### Biogeographic
pdb_NMDS <- metaMDS(pdb_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
pdbms_NMDS <- metaMDS(pdbms_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
pdbom_NMDS <- metaMDS(pdbom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
pdbso_NMDS <- metaMDS(pdbso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
```

## Determine stratification homogeneity using phylogenetic distances

We use betadisper to determine homogeneity of the phylogenetic distances
based on the stratification category. Is the dispersion of phylogenetic
distances similar within stratification categories?

``` r
### Stratification
pd_bd <- betadisper(pd_dist, stratification_group)
anova(pd_bd)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value  Pr(>F)  
    ## Groups     2 0.041384 0.0206918  2.8636 0.08185 .
    ## Residuals 19 0.137291 0.0072258                  
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
    ## Groups     2 0.041384 0.0206918 2.8636    999  0.084 .
    ## Residuals 19 0.137291 0.0072258                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##               Mixed    Ocean Stratified
    ## Mixed               0.792000      0.041
    ## Ocean      0.785767               0.134
    ## Stratified 0.045423 0.120915

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
    ## Ocean-Mixed       0.01016167 -0.1064651 0.12678849 0.9733893
    ## Stratified-Mixed -0.08541898 -0.1933944 0.02255647 0.1370417
    ## Stratified-Ocean -0.09558066 -0.2122075 0.02104616 0.1203813

``` r
boxplot(pd_bd)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Determine%20stratification%20homogeneity%20using%20phylogenetic%20distances-1.png)<!-- -->

## Determine stratification dissimilarity using phylogenetic distances

We use adonis to determine if dissimilarities of species phylogenetic
distances by stratification categories are significant and how much of
the variation is explained by phylogenetic dissimilarities.

``` r
### Stratification
pd_pms <- adonis2(pd_dist ~ env[-20,19], permutations = 999)
pd_pms
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = pd_dist ~ env[-20, 19], permutations = 999)
    ##              Df SumOfSqs      R2      F Pr(>F)    
    ## env[-20, 19]  2   1.8501 0.39099 6.0992  0.001 ***
    ## Residual     19   2.8817 0.60901                  
    ## Total        21   4.7317 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratification
pd_pms_pair <- pairwise.adonis(pd_dist, env[-20,19], p.adjust.m = "bonferroni", perm = 999)
pd_pms_pair
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.2941262 9.411578 0.4020053   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.1732986 8.193920 0.4057617   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.2549694 1.443238 0.1073579   0.106      0.318

## Envfit environmental influence on phylogenetic distances

We use envfit to determine significantly correlated environmental
variables to our phylogenetic NMDS. We will use these results, if
significant, in our figure.

``` r
### Environmental
# All sites for figure
pde_ef <- envfit(pde_NMDS, env[-c(9,16,18,20),c(34,36)], permutations = 999, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
pde_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##       NMDS1    NMDS2     r2 Pr(>r)    
    ## S   0.96421  0.26513 0.7072  0.017 *  
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
    ## S   0.96421  0.26513 0.7072  0.034 * 
    ## pH  0.99017 -0.13990 0.8530  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# All sites by stratification
pdea_ef <- envfit(pde_NMDS, env[-c(9,16,18,20),c(2,6,8,10)], permutations = 999, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
pdea_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)    
    ## temperature_median -0.73729 -0.67558 0.0772  0.622    
    ## salinity_median     0.96421  0.26513 0.7072  0.011 *  
    ## oxygen_median       0.56896 -0.82236 0.6630  0.191    
    ## pH_median           0.99017 -0.13990 0.8530  0.001 ***
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
    ## temperature_median -0.73729 -0.67558 0.0772  1.000   
    ## salinity_median     0.96421  0.26513 0.7072  0.044 * 
    ## oxygen_median       0.56896 -0.82236 0.6630  0.764   
    ## pH_median           0.99017 -0.13990 0.8530  0.004 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
pdems_ef <- envfit(pdems_NMDS, env[c(1:6,8,10,12:15,17,21:23),c(6,10)], permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
pdems_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                   NMDS1   NMDS2     r2 Pr(>r)   
    ## salinity_median 0.86013 0.51008 0.6978  0.013 * 
    ## pH_median       0.99082 0.13519 0.8248  0.002 **
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
    ##                   NMDS1   NMDS2     r2 Pr(>r)   
    ## salinity_median 0.86013 0.51008 0.6978  0.026 * 
    ## pH_median       0.99082 0.13519 0.8248  0.004 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
pdeom_ef <- envfit(pdeom_NMDS, env[c(3,6:8,10:11,13:15,19,23),c(8,10)], permutations = 999, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
pdeom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                  NMDS1    NMDS2     r2 Pr(>r)  
    ## oxygen_median  0.93268  0.36070 0.7964  0.014 *
    ## pH_median      0.95280 -0.30359 0.7303  0.011 *
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
    ##                  NMDS1    NMDS2     r2 Pr(>r)  
    ## oxygen_median  0.93268  0.36070 0.7964  0.028 *
    ## pH_median      0.95280 -0.30359 0.7303  0.022 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
pdeso_ef <- envfit(pdeso_NMDS, env[c(1:2,4,5,7,11:12,17,19,21:22),c(2,6,8,10)], permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,11:12,17,19,21:22),19])
pdeso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)  
    ## temperature_median -0.00030110  1.00000000 0.0570  0.637  
    ## salinity_median     0.00041739  1.00000000 0.7126  0.053 .
    ## oxygen_median       0.00077290 -1.00000000 0.8105  0.077 .
    ## pH_median           0.00055988  1.00000000 0.8454  0.018 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ##                          NMDS1       NMDS2     r2 Pr(>r)  
    ## temperature_median -0.00030110  1.00000000 0.0570  1.000  
    ## salinity_median     0.00041739  1.00000000 0.7126  0.212  
    ## oxygen_median       0.00077290 -1.00000000 0.8105  0.308  
    ## pH_median           0.00055988  1.00000000 0.8454  0.072 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes
pds_ef <- envfit(pds_NMDS, env[c(1:2,4,5,12,17,21:22),c(8)], permutations = 999, na.rm = TRUE)
pds_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##        NMDS1   NMDS2    r2 Pr(>r)  
    ## [1,] 0.83740 0.54659 0.745  0.035 *
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
    ##        NMDS1   NMDS2    r2 Pr(>r)  
    ## [1,] 0.83740 0.54659 0.745  0.035 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Biogeographic
# All sites for figure
pdb_ef <- envfit(pdb_NMDS, env[-c(8,20),c(40)], permutations = 1000, na.rm = TRUE, strata = env[-c(8,20),19])
pdb_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## [1,] -0.92159  0.38816 0.5898 0.06094 .
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
    ##         NMDS1    NMDS2     r2  Pr(>r)  
    ## [1,] -0.92159  0.38816 0.5898 0.06094 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
pdba_ef <- envfit(pdb_NMDS, env[-c(8,20),c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[-c(8,20),19])
pdba_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.30704 -0.95170 0.2742 0.14186  
    ## volume_m3                   0.29231 -0.95632 0.2674 0.14186  
    ## surface_area_m2             0.24090 -0.97055 0.3746 0.17383  
    ## distance_to_ocean_min_m    -0.92159  0.38816 0.5898 0.07393 .
    ## distance_to_ocean_mean_m   -0.76126  0.64845 0.5394 0.48951  
    ## distance_to_ocean_median_m -0.76948  0.63867 0.5308 0.49251  
    ## tidal_lag_time_minutes     -0.63880  0.76938 0.6134 0.32468  
    ## tidal_efficiency            0.63770 -0.77028 0.6175 0.22577  
    ## perimeter_fromSat           0.31940 -0.94762 0.2854 0.40460  
    ## max_depth                  -0.30949 -0.95090 0.0695 0.96703  
    ## logArea                     0.24320 -0.96998 0.1774 0.40060  
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
    ## volume_m3_w_chemocline      0.30704 -0.95170 0.2742 1.0000
    ## volume_m3                   0.29231 -0.95632 0.2674 1.0000
    ## surface_area_m2             0.24090 -0.97055 0.3746 1.0000
    ## distance_to_ocean_min_m    -0.92159  0.38816 0.5898 0.8132
    ## distance_to_ocean_mean_m   -0.76126  0.64845 0.5394 1.0000
    ## distance_to_ocean_median_m -0.76948  0.63867 0.5308 1.0000
    ## tidal_lag_time_minutes     -0.63880  0.76938 0.6134 1.0000
    ## tidal_efficiency            0.63770 -0.77028 0.6175 1.0000
    ## perimeter_fromSat           0.31940 -0.94762 0.2854 1.0000
    ## max_depth                  -0.30949 -0.95090 0.0695 1.0000
    ## logArea                     0.24320 -0.96998 0.1774 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
pdbms_ef <- envfit(pdbms_NMDS, env[c(1:6,10,12:15,17,21:23),c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[c(1:6,10,12:15,17,21:23),19])
pdbms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline     -0.06129  0.99812 0.0721 0.82318  
    ## volume_m3                  -0.49381  0.86957 0.1345 0.96503  
    ## surface_area_m2            -0.34954  0.93692 0.1072 0.93706  
    ## distance_to_ocean_min_m    -0.83987  0.54279 0.5147 0.07692 .
    ## distance_to_ocean_mean_m   -0.59314  0.80510 0.4517 0.37363  
    ## distance_to_ocean_median_m -0.58207  0.81314 0.4510 0.34565  
    ## tidal_lag_time_minutes     -0.39660  0.91799 0.6069 0.24775  
    ## tidal_efficiency            0.40488 -0.91437 0.6010 0.19281  
    ## perimeter_fromSat          -0.09999  0.99499 0.0799 0.79520  
    ## max_depth                  -0.56729 -0.82352 0.1064 0.96104  
    ## logArea                    -0.52584  0.85059 0.0070 0.99700  
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
    ## volume_m3_w_chemocline     -0.06129  0.99812 0.0721 1.0000
    ## volume_m3                  -0.49381  0.86957 0.1345 1.0000
    ## surface_area_m2            -0.34954  0.93692 0.1072 1.0000
    ## distance_to_ocean_min_m    -0.83987  0.54279 0.5147 0.8462
    ## distance_to_ocean_mean_m   -0.59314  0.80510 0.4517 1.0000
    ## distance_to_ocean_median_m -0.58207  0.81314 0.4510 1.0000
    ## tidal_lag_time_minutes     -0.39660  0.91799 0.6069 1.0000
    ## tidal_efficiency            0.40488 -0.91437 0.6010 1.0000
    ## perimeter_fromSat          -0.09999  0.99499 0.0799 1.0000
    ## max_depth                  -0.56729 -0.82352 0.1064 1.0000
    ## logArea                    -0.52584  0.85059 0.0070 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
pdbom_ef <- envfit(pdbom_NMDS, env[c(3,6:7,9:11,13:16,18:19,23),c(28,31)], permutations = 1000, na.rm = TRUE, strata = env[c(3,6:7,9:11,13:16,18:19,23),19])
pdbom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                           NMDS1    NMDS2     r2   Pr(>r)   
    ## tidal_lag_time_minutes -0.34047 -0.94025 0.5044 0.011988 * 
    ## max_depth               0.87692 -0.48064 0.6324 0.004995 **
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
    ##                           NMDS1    NMDS2     r2  Pr(>r)   
    ## tidal_lag_time_minutes -0.34047 -0.94025 0.5044 0.02398 * 
    ## max_depth               0.87692 -0.48064 0.6324 0.00999 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
pdbso_ef <- envfit(pdbso_NMDS, env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(22:32)], permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
pdbso_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1    NMDS2     r2  Pr(>r)  
    ## volume_m3_w_chemocline      0.32433  0.94594 0.4850 0.03197 *
    ## volume_m3                   0.31779  0.94816 0.4742 0.04096 *
    ## surface_area_m2             0.30237  0.95319 0.6099 0.03596 *
    ## distance_to_ocean_min_m    -0.99834  0.05759 0.6213 0.26474  
    ## distance_to_ocean_mean_m   -0.89864  0.43869 0.6226 0.56543  
    ## distance_to_ocean_median_m -0.87456  0.48491 0.6154 0.55644  
    ## tidal_lag_time_minutes     -0.99917  0.04079 0.7125 0.93606  
    ## tidal_efficiency            0.93912 -0.34358 0.7835 0.56843  
    ## perimeter_fromSat           0.38918  0.92116 0.5288 0.11988  
    ## max_depth                  -0.16212  0.98677 0.1849 0.55744  
    ## logArea                     0.26902  0.96314 0.4195 0.11289  
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
    ## volume_m3_w_chemocline      0.32433  0.94594 0.4850 0.3516
    ## volume_m3                   0.31779  0.94816 0.4742 0.4505
    ## surface_area_m2             0.30237  0.95319 0.6099 0.3956
    ## distance_to_ocean_min_m    -0.99834  0.05759 0.6213 1.0000
    ## distance_to_ocean_mean_m   -0.89864  0.43869 0.6226 1.0000
    ## distance_to_ocean_median_m -0.87456  0.48491 0.6154 1.0000
    ## tidal_lag_time_minutes     -0.99917  0.04079 0.7125 1.0000
    ## tidal_efficiency            0.93912 -0.34358 0.7835 1.0000
    ## perimeter_fromSat           0.38918  0.92116 0.5288 1.0000
    ## max_depth                  -0.16212  0.98677 0.1849 1.0000
    ## logArea                     0.26902  0.96314 0.4195 1.0000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Traits
# All sites for figure
pdt_ef <- envfit(pd_NMDS, FD_total_env[-20,c(58,64)], permutations = 1000, na.rm = TRUE, strata = env[-20,19])
pdt_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##      NMDS1    NMDS2     r2   Pr(>r)   
    ## T -0.71490  0.69923 0.7028 0.004995 **
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
    ## EC1b -0.3816 -0.0270
    ## EC2g -0.4873 -0.0630
    ## EC3n  0.0965  0.0100
    ## 
    ## Goodness of fit:
    ##       r2   Pr(>r)   
    ## EC 0.432 0.006993 **
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
    ##      NMDS1    NMDS2     r2  Pr(>r)   
    ## T -0.71490  0.69923 0.7028 0.00999 **
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
    ## EC1b -0.3816 -0.0270
    ## EC2g -0.4873 -0.0630
    ## EC3n  0.0965  0.0100
    ## 
    ## Goodness of fit:
    ##       r2  Pr(>r)  
    ## EC 0.432 0.01399 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
pdta_ef <- envfit(pd_NMDS, FD_total_env[-20,c(6:22)], permutations = 1000, na.rm = TRUE, strata = env[-20,19])
pdta_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2   Pr(>r)   
    ## MaxLengthTL       0.19595  0.98061 0.2103 0.050949 . 
    ## Troph            -0.71490  0.69923 0.7028 0.007992 **
    ## DepthMin         -0.10558 -0.99441 0.0497 0.177822   
    ## DepthMax          0.74848  0.66316 0.2306 0.161838   
    ## TempPrefMin       0.56814 -0.82293 0.1991 0.395604   
    ## TempPrefMax      -0.39503 -0.91867 0.0052 0.820180   
    ## Weight            0.08261  0.99658 0.0847 0.207792   
    ## CaudalFinLength  -0.22744  0.97379 0.0329 0.433566   
    ## DorsalSpinesMean  0.97324 -0.22979 0.5299 0.042957 * 
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
    ## FeedingPathp        -0.3816 -0.0270
    ## RepGuild11b         -0.3816 -0.0270
    ## RepGuild12g         -0.4873 -0.0630
    ## RepGuild13n          0.0965  0.0100
    ## RepGuild22eb        -0.3955 -0.0180
    ## RepGuild23n         -0.4873 -0.0630
    ## RepGuild26s          0.0721  0.0076
    ## ParentalCare3p      -0.3755 -0.0197
    ## ParentalCare4n       0.1408  0.0074
    ## WaterPref1s          0.1177  0.0096
    ## WaterPref3a         -0.4003 -0.0327
    ## 
    ## Goodness of fit:
    ##                      r2   Pr(>r)   
    ## BodyShapeI       0.4083 0.041958 * 
    ## DemersPelag      0.0000 1.000000   
    ## OperculumPresent 0.3415 0.155844   
    ## FeedingPath      0.1471 0.350649   
    ## RepGuild1        0.4320 0.001998 **
    ## RepGuild2        0.3380 0.005994 **
    ## ParentalCare     0.5331 0.135864   
    ## WaterPref        0.4770 0.054945 . 
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
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.19595  0.98061 0.2103 0.8661
    ## Troph            -0.71490  0.69923 0.7028 0.1359
    ## DepthMin         -0.10558 -0.99441 0.0497 1.0000
    ## DepthMax          0.74848  0.66316 0.2306 1.0000
    ## TempPrefMin       0.56814 -0.82293 0.1991 1.0000
    ## TempPrefMax      -0.39503 -0.91867 0.0052 1.0000
    ## Weight            0.08261  0.99658 0.0847 1.0000
    ## CaudalFinLength  -0.22744  0.97379 0.0329 1.0000
    ## DorsalSpinesMean  0.97324 -0.22979 0.5299 0.7303
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
    ## FeedingPathp        -0.3816 -0.0270
    ## RepGuild11b         -0.3816 -0.0270
    ## RepGuild12g         -0.4873 -0.0630
    ## RepGuild13n          0.0965  0.0100
    ## RepGuild22eb        -0.3955 -0.0180
    ## RepGuild23n         -0.4873 -0.0630
    ## RepGuild26s          0.0721  0.0076
    ## ParentalCare3p      -0.3755 -0.0197
    ## ParentalCare4n       0.1408  0.0074
    ## WaterPref1s          0.1177  0.0096
    ## WaterPref3a         -0.4003 -0.0327
    ## 
    ## Goodness of fit:
    ##                      r2  Pr(>r)  
    ## BodyShapeI       0.4083 0.71329  
    ## DemersPelag      0.0000 1.00000  
    ## OperculumPresent 0.3415 1.00000  
    ## FeedingPath      0.1471 1.00000  
    ## RepGuild1        0.4320 0.03397 *
    ## RepGuild2        0.3380 0.10190  
    ## ParentalCare     0.5331 1.00000  
    ## WaterPref        0.4770 0.93407  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Mixed and stratified lakes
pdtms_ef <- envfit(pdms_NMDS, FD_total_env[c(1:6,8,10,12:15,17,21:23),c(18,19)], permutations = 1000, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
pdtms_ef
```

    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                NMDS1   NMDS2
    ## RepGuild11b  -0.3099 -0.0094
    ## RepGuild12g  -0.4096 -0.0753
    ## RepGuild13n   0.1199  0.0141
    ## RepGuild22eb -0.3167 -0.0114
    ## RepGuild23n  -0.4096 -0.0753
    ## RepGuild26s   0.0874  0.0125
    ## 
    ## Goodness of fit:
    ##               r2   Pr(>r)   
    ## RepGuild1 0.4570 0.007992 **
    ## RepGuild2 0.3511 0.012987 * 
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
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                NMDS1   NMDS2
    ## RepGuild11b  -0.3099 -0.0094
    ## RepGuild12g  -0.4096 -0.0753
    ## RepGuild13n   0.1199  0.0141
    ## RepGuild22eb -0.3167 -0.0114
    ## RepGuild23n  -0.4096 -0.0753
    ## RepGuild26s   0.0874  0.0125
    ## 
    ## Goodness of fit:
    ##               r2  Pr(>r)  
    ## RepGuild1 0.4570 0.01598 *
    ## RepGuild2 0.3511 0.02597 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
pdtom_ef <- envfit(pdom_NMDS, FD_total_env[c(3,6:11,13:16,18:19,23),c(9)], permutations = 1000, na.rm = TRUE, strata = env[c(3,6:11,13:16,18:19,23),19])
pdtom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##         NMDS1    NMDS2   r2   Pr(>r)   
    ## [1,] -0.67324  0.73943 0.56 0.005994 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ##         NMDS1    NMDS2   r2   Pr(>r)   
    ## [1,] -0.67324  0.73943 0.56 0.005994 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
pdtso_ef <- envfit(pdso_NMDS, FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(18,19)], permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
pdtso_ef
```

    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                NMDS1   NMDS2
    ## RepGuild11b  -0.2537 -0.0082
    ## RepGuild12g  -0.3650  0.0291
    ## RepGuild13n   0.1237 -0.0042
    ## RepGuild22eb -0.2680 -0.0175
    ## RepGuild23n  -0.3650  0.0291
    ## RepGuild26s   0.0907 -0.0037
    ## 
    ## Goodness of fit:
    ##               r2   Pr(>r)   
    ## RepGuild1 0.4085 0.006993 **
    ## RepGuild2 0.3199 0.012987 * 
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
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                NMDS1   NMDS2
    ## RepGuild11b  -0.2537 -0.0082
    ## RepGuild12g  -0.3650  0.0291
    ## RepGuild13n   0.1237 -0.0042
    ## RepGuild22eb -0.2680 -0.0175
    ## RepGuild23n  -0.3650  0.0291
    ## RepGuild26s   0.0907 -0.0037
    ## 
    ## Goodness of fit:
    ##               r2  Pr(>r)  
    ## RepGuild1 0.4085 0.01399 *
    ## RepGuild2 0.3199 0.02597 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

## Phylogenetic mantel tests

We used mantel tests to determine significance between the phylogenetic
and environmental distance matrices.

``` r
### Environmental
# All sites 
enve_dist_s <- dist(scaled_env[-c(9,16,18,20),c(2)], method = "euclidean")
pdea_mant_s <- mantel(pde_dist, enve_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
pdea_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pde_dist, ydis = enve_dist_s, method = "spearman",      permutations = 999, strata = env[-c(9, 16, 18, 20), 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5397 
    ##       Significance: 0.004 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.452 0.480 0.496 0.515 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enve_dist_p <- dist(scaled_env[-c(9,16,18,20),c(4)], method = "euclidean")
pdea_mant_p <- mantel(pde_dist, enve_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
pdea_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pde_dist, ydis = enve_dist_p, method = "spearman",      permutations = 999, strata = env[-c(9, 16, 18, 20), 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7157 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.509 0.544 0.582 0.612 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdea_mant_pv <- rbind(pdea_mant_s$signif, pdea_mant_p$signif)
pdea_mant_pv <- pdea_mant_pv[,1]
pdea_mant_pv <- p.adjust(pdea_mant_pv, method = "bonferroni")
pdea_mant_pv
```

    ## [1] 0.008 0.002

``` r
# Mixed and stratified lakes
envems_dist_s <- dist(scaled_env[c(1:6,8,10,12:15,17,21:23),c(2)], method = "euclidean")
pdems_mant_s <- mantel(pdems_dist, envems_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
pdems_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdems_dist, ydis = envems_dist_s, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4842 
    ##       Significance: 0.011 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.394 0.419 0.444 0.481 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
envems_dist_p <- dist(scaled_env[c(1:6,8,10,12:15,17,21:23),c(4)], method = "euclidean")
pdems_mant_p <- mantel(pdems_dist, envems_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
pdems_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdems_dist, ydis = envems_dist_p, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6732 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.435 0.472 0.519 0.576 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdems_mant_pv <- rbind(pdems_mant_s$signif, pdems_mant_p$signif)
pdems_mant_pv <- pdems_mant_pv[,1]
pdems_mant_pv <- p.adjust(pdems_mant_pv, method = "bonferroni")
pdems_mant_pv
```

    ## [1] 0.022 0.002

``` r
# Ocean sites and mixed lakes
enveom_dist_s <- dist(scaled_env[c(3,6:8,10:11,13:15,19,23),c(2)], method = "euclidean")
pdeom_mant_s <- mantel(pdeom_dist, enveom_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
pdeom_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeom_dist, ydis = enveom_dist_s, method = "spearman",      permutations = 999, strata = env[c(3, 6:8, 10:11, 13:15,          19, 23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5084 
    ##       Significance: 0.01 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.313 0.383 0.443 0.507 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
enveom_dist_p <- dist(scaled_env[c(3,6:8,10:11,13:15,19,23),c(4)], method = "euclidean")
pdeom_mant_p <- mantel(pdeom_dist, enveom_dist_p, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
pdeom_mant_p
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeom_dist, ydis = enveom_dist_p, method = "spearman",      permutations = 999, strata = env[c(3, 6:8, 10:11, 13:15,          19, 23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4535 
    ##       Significance: 0.018 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.261 0.319 0.388 0.516 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdeom_mant_pv <- rbind(pdeom_mant_s$signif, pdeom_mant_p$signif)
pdeom_mant_pv <- pdeom_mant_pv[,1]
pdeom_mant_pv <- p.adjust(pdeom_mant_pv, method = "bonferroni")
pdeom_mant_pv
```

    ## [1] 0.020 0.036

``` r
# Stratified lakes and ocean sites
enveso_dist_s <- dist(scaled_env[c(1:2,4,5,7,11:12,17,19,21:22),c(2)], method = "euclidean")
pdeso_mant_s <- mantel(pdeso_dist, enveso_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,11:12,17,19,21:22),19])
pdeso_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdeso_dist, ydis = enveso_dist_s, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 11:12, 17,          19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5084 
    ##       Significance: 0.024 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.418 0.468 0.503 0.543 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes
envs_dist <- dist(scaled_env[c(1:2,4,5,12,17,21:22),c(1:15)], method = "euclidean")
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
    ##       Significance: 0.291 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.337 0.406 0.464 0.586 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Biogeographic
# All sites 
envb_dist <- dist(scaled_env[-c(8,20),c(5:15)], method = "euclidean")
pdb_mant <- mantel(pdb_dist, envb_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(8,20),19])
pdb_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdb_dist, ydis = envb_dist, method = "spearman",      permutations = 999, strata = env[-c(8, 20), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3595 
    ##       Significance: 0.336 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.413 0.438 0.459 0.483 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
envbms_dist <- dist(scaled_env[c(1:6,10,12:15,17,21:23),c(5:15)], method = "euclidean")
pdbms_mant <- mantel(pdbms_dist, envbms_dist, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:6,10,12:15,17,21:23),19])
pdbms_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbms_dist, ydis = envbms_dist, method = "spearman",      permutations = 1000, strata = env[c(1:6, 10, 12:15, 17, 21:23),          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2029 
    ##       Significance: 0.58042 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.384 0.432 0.464 0.496 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
envbom_dist <- dist(scaled_env[c(3,6:7,9:11,13:16,18:19,23),c(6)], method = "euclidean")
pdbom_mant <- mantel(pdbom_dist, envbom_dist, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(3,6:7,9:11,13:16,18:19,23),19])
pdbom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbom_dist, ydis = envbom_dist, method = "spearman",      permutations = 1000, strata = env[c(3, 6:7, 9:11, 13:16,          18:19, 23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2984 
    ##       Significance: 0.032967 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.233 0.268 0.305 0.343 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
envbso_dist <- dist(scaled_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(5:15)], method = "euclidean")
pdbso_mant <- mantel(pdbso_dist, envbso_dist, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
pdbso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdbso_dist, ydis = envbso_dist, method = "spearman",      permutations = 1000, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6341 
    ##       Significance: 0.25075 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.673 0.689 0.709 0.735 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
### Traits
# All sites 1,3,12:13,14:15
FDpd_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[-c(20),c(1)]))
pdt_mant_s <- mantel(pd_dist, FDpd_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(20),19])
pdt_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_s, method = "spearman",      permutations = 999, strata = env[-c(20), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4615 
    ##       Significance: 0.004 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.316 0.347 0.370 0.416 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_t <- gowdis(as.data.frame(scaled_FD_total_env[-c(20),c(3)]))
pdt_mant_t <- mantel(pd_dist, FDpd_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(20),19])
pdt_mant_t
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_t, method = "spearman",      permutations = 999, strata = env[-c(20), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4874 
    ##       Significance: 0.005 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.396 0.427 0.442 0.463 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[-c(20),c(12:13)]))
pdt_mant_ec <- mantel(pd_dist, FDpd_dist_ec, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(20),19])
pdt_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_ec, method = "spearman",      permutations = 999, strata = env[-c(20), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3755 
    ##       Significance: 0.005 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.264 0.302 0.329 0.351 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
FDpd_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[-c(20),c(14:15)]))
pdt_mant_es <- mantel(pd_dist, FDpd_dist_es, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(20),19])
pdt_mant_es
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pd_dist, ydis = FDpd_dist_es, method = "spearman",      permutations = 999, strata = env[-c(20), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3137 
    ##       Significance: 0.009 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.211 0.235 0.271 0.299 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
pdt_mant_pv <- rbind(pdt_mant_s$signif, pdt_mant_t$signif, pdt_mant_ec$signif, pdt_mant_es$signif)
pdt_mant_pv <- pdt_mant_pv[,1]
pdt_mant_pv <- p.adjust(pdt_mant_pv, method = "bonferroni")
pdt_mant_pv
```

    ## [1] 0.016 0.020 0.020 0.036

``` r
# Mixed and stratified lakes 1,12:13,14:15
FDpdms_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[c(1:6,8,10,12:15,17,21:23),c(1)]))
pdtms_mant_s <- mantel(pdms_dist, FDpdms_dist_s, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
pdtms_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_s, method = "spearman",      permutations = 1000, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3773 
    ##       Significance: 0.00999 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.252 0.297 0.330 0.368 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
FDpdms_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[c(1:6,8,10,12:15,17,21:23),c(12:13)]))
pdtms_mant_ec <- mantel(pdms_dist, FDpdms_dist_ec, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
pdtms_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_ec, method = "spearman",      permutations = 1000, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2809 
    ##       Significance: 0.000999 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.138 0.185 0.216 0.247 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
FDpdms_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[c(1:6,8,10,12:15,17,21:23),c(14:15)]))
pdtms_mant_es <- mantel(pdms_dist, FDpdms_dist_es, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
pdtms_mant_es
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdms_dist, ydis = FDpdms_dist_es, method = "spearman",      permutations = 1000, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2143 
    ##       Significance: 0.003996 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0845 0.1184 0.1635 0.1874 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Adjust p-values
pdtms_mant_pv <- rbind(pdtms_mant_s$signif, pdtms_mant_ec$signif, pdtms_mant_es$signif)
pdtms_mant_pv <- pdtms_mant_pv[,1]
pdtms_mant_pv <- p.adjust(pdtms_mant_pv, method = "bonferroni")
pdtms_mant_pv
```

    ## [1] 0.029970030 0.002997003 0.011988012

``` r
# Ocean sites and mixed lakes 1,2,6
FDpdom_dist_s <- gowdis(as.data.frame(scaled_FD_total_env[c(3,6:11,13:16,18:19,23),c(1)]))
pdtom_mant_s <- mantel(pdom_dist, FDpdom_dist_s, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(3,6:11,13:16,18:19,23),19])
pdtom_mant_s
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdom_dist, ydis = FDpdom_dist_s, method = "spearman",      permutations = 1000, strata = env[c(3, 6:11, 13:16, 18:19,          23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3449 
    ##       Significance: 0.006993 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.152 0.199 0.250 0.318 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
FDpdom_dist_l <- gowdis(as.data.frame(scaled_FD_total_env[c(3,6:11,13:16,18:19,23),c(2)]))
pdtom_mant_l <- mantel(pdom_dist, FDpdom_dist_l, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(3,6:11,13:16,18:19,23),19])
pdtom_mant_l
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdom_dist, ydis = FDpdom_dist_l, method = "spearman",      permutations = 1000, strata = env[c(3, 6:11, 13:16, 18:19,          23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2589 
    ##       Significance: 0.00999 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0966 0.1464 0.2071 0.2507 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Adjust p-values
pdtom_mant_pv <- rbind(pdtom_mant_s$signif, pdtom_mant_l$signif)
pdtom_mant_pv <- pdtom_mant_pv[,1]
pdtom_mant_pv <- p.adjust(pdtom_mant_pv, method = "bonferroni")
pdtom_mant_pv
```

    ## [1] 0.01398601 0.01998002

``` r
# Stratified lakes and ocean sites 1,3,4,7,12:13,14:15
FDpdso_dist_dmin <- gowdis(as.data.frame(scaled_FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(4)]))
pdtso_mant_dmin <- mantel(pdso_dist, FDpdso_dist_dmin, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
pdtso_mant_dmin
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_dmin, method = "spearman",      permutations = 1000, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6478 
    ##       Significance: 0.000999 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.508 0.534 0.550 0.568 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
FDpdso_dist_tmax <- gowdis(as.data.frame(scaled_FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(7)]))
pdtso_mant_tmax <- mantel(pdso_dist, FDpdso_dist_tmax, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
pdtso_mant_tmax
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_tmax, method = "spearman",      permutations = 1000, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.117 
    ##       Significance: 0.001998 
    ## 
    ## Upper quantiles of permutations (null model):
    ##     90%     95%   97.5%     99% 
    ## 0.00632 0.02959 0.04652 0.06397 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
FDpdso_dist_ec <- gowdis(as.data.frame(scaled_FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(12:13)]))
pdtso_mant_ec <- mantel(pdso_dist, FDpdso_dist_ec, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
pdtso_mant_ec
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_ec, method = "spearman",      permutations = 1000, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1961 
    ##       Significance: 0.001998 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0548 0.0883 0.1159 0.1431 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
FDpdso_dist_es <- gowdis(as.data.frame(scaled_FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(14:15)]))
pdtso_mant_es <- mantel(pdso_dist, FDpdso_dist_es, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
pdtso_mant_es
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = pdso_dist, ydis = FDpdso_dist_es, method = "spearman",      permutations = 1000, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1341 
    ##       Significance: 0.008991 
    ## 
    ## Upper quantiles of permutations (null model):
    ##     90%     95%   97.5%     99% 
    ## 0.00968 0.06506 0.09480 0.10807 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Adjust p-values
pdtso_mant_pv <- rbind(pdtso_mant_dmin$signif, pdtso_mant_tmax$signif, pdtso_mant_ec$signif, pdtso_mant_es$signif)
pdtso_mant_pv <- pdtso_mant_pv[,1]
pdtso_mant_pv <- p.adjust(pdtso_mant_pv, method = "bonferroni")
pdtso_mant_pv
```

    ## [1] 0.003996004 0.007992008 0.007992008 0.035964036

## Plot phylogenetic NMDS and envfit results

Includes a plot of correlated environmental variables

``` r
pd_NMDS_data.scores <- as.data.frame(scores(pd_NMDS))
pd_NMDS_data.scores$Stratification <- env[-c(20),19]
pd_NMDS_data.scores$Lakes <- env[-c(20),1]
pd_NMDS_data.scores$Stratification <- factor(pd_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

pde_ef_coord_cont <- as.data.frame(scores(pde_ef, "vectors")) * ordiArrowMul(pde_ef)
# pdb_ef_coord_cont <- as.data.frame(scores(pdb_ef, "vectors")) * ordiArrowMul(pdb_ef)
# pdb_row_names <- c("minD")
# # Assign the new row names to the data frame
# pdb_ef_coord_cont <- data.frame(row.names = pdb_row_names, pdb_ef_coord_cont)

pd_ef_plot <- ggplot(data = pd_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification),level = 0.50) +
  geom_point(data = pd_NMDS_data.scores, aes(color = Stratification), size = 4, alpha = 1) + 
  geom_text_repel(data = pd_NMDS_data.scores, label = pd_NMDS_data.scores$Lakes, size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = pde_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#CD661D") +
  geom_text(data = pde_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
            color = "#CD661D", label = row.names(pde_ef_coord_cont), size = 7) + 
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
  #              data = pdb_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#698B22") +
  # geom_text(data = pdb_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
  #           color = "#698B22", label = row.names(pdb_ef_coord_cont), size = 7) +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
        panel.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  # annotate(geom = "label", x = .29, y = .22, size = 6, 
  #          label = paste("Stress: ", round(pd_NMDS$stress, digits = 2))) +
  labs(color = "Stratification")
pd_ef_plot <- pd_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
pd_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20phylogenetic%20NMDS%20and%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_ef_plot.png", pd_ef_plot, width = 8, height = 4, units = "in")
```

## Boxplot of intradissimilarity heterogeneity

The dissimilarity between sites that share the same stratification.

``` r
pd_bd_dist <- pd_bd$distances
pd_bd_dist <- as.data.frame(pd_bd_dist)
pd_bd_dist$X <- row.names(pd_bd_dist)
pd_bd_dist_env <- merge(pd_bd_dist, env[-20,], by = "X", sort = F)

pd_bd_dist_env_plot <- ggplot(pd_bd_dist_env, aes(x = Stratification, y = pd_bd_dist, 
                                                  color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.75,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black")) +
  ylab("Distance to Centroid") +
  labs(color = "Stratification", tag = "B")
pd_bd_dist_env_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Boxplot%20of%20phylogenetic%20intradissimilarity%20heterogeneity-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pd_bd_dist_env_plot.png", pd_bd_dist_env_plot, width = 8.25, height = 4.13, units = "in")
```

## Plot phylogenetic NMDS and trait envfit results

Includes a plot of correlated trait variables

``` r
pdt_ef_coord_cont <- as.data.frame(scores(pdt_ef, "vectors")) * ordiArrowMul(pdt_ef)
pdt_ef_coord_cat = as.data.frame(scores(pdt_ef, "factors")) * ordiArrowMul(pdt_ef)

pdt_ef_plot <- ggplot(data = pd_NMDS_data.scores, aes(x = NMDS1, y = NMDS2)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification),level = 0.50) +
  geom_point(data = pd_NMDS_data.scores, aes(color = Stratification), size = 4, alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = pdt_ef_coord_cont, size =1, alpha = 0.5, color = "#CD3333") +
  geom_text(data = pdt_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
            color = "#CD3333", label = row.names(pdt_ef_coord_cont)) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = pdt_ef_coord_cat, size =1, alpha = 0.5, color = "#8B6508") +
  geom_text(data = pdt_ef_coord_cat, aes(x = NMDS1, y = NMDS2), 
            color = "#8B6508", label = row.names(pdt_ef_coord_cat)) + 
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
        panel.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  # annotate(geom = "label", x = .29, y = .22, size = 6, 
  #          label = paste("Stress: ", round(pd_NMDS$stress, digits = 2))) +
  labs(color = "Stratification")
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
pdt_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20phylogenetic%20NMDS%20and%20trait%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/pdt_ef_plot.png", pdt_ef_plot, width = 8, height = 4, units = "in")
```

## Environment/biogeography dissimilarity distances

``` r
### Environmental
# All sites
enve_dist <- dist(scaled_env[-c(9,16,18,20),c(2:4)], method = "euclidean")
# Mixed and stratified lakes
envems_dist <- dist(scaled_env[c(1:6,8,10,12:15,17,21:23),c(2:4)], method = "euclidean")
# Ocean sites and mixed lakes
enveom_dist <- dist(scaled_env[c(3,6:8,10:11,13:15,19,23),c(2:4)], method = "euclidean")
# Stratified lakes and ocean sites
enveso_dist <- dist(scaled_env[c(1:2,4,5,7,11:12,17,19,21:22),c(2:4)], method = "euclidean")

### Biogeographic
# All sites
envb_dist <- dist(scaled_env[-c(8,20),c(5:6,8,12:15)], method = "euclidean")
# Mixed and stratified lakes
envbms_dist <- dist(scaled_env[c(1:6,10,12:15,17,21:23),c(5:6,8,12:15)], method = "euclidean")
# Ocean sites and mixed lakes
envbom_dist <- dist(scaled_env[c(3,6:7,9:11,13:16,18:19,23),c(5:6,8,12:15)], method = "euclidean")
# Stratified lakes and ocean sites
envbso_dist <- dist(scaled_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(5:6,8,12:15)], method = "euclidean")
```

## Environment/biogeography NMDS

To constrain dissimilarities we perform Nonmetric Multidimensional
Scaling (NMDS), which tries to find a stable solution using the metaMDS
package.

``` r
### Environmental
# All sites 
enve_NMDS <- metaMDS(enve_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
envems_NMDS <- metaMDS(envems_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
enveom_NMDS <- metaMDS(enveom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(enveom_dist, try = 1000, parallel = 4, previous.best, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Stratified lakes and ocean sites
enveso_NMDS <- metaMDS(enveso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)

### Biogeographic
# All sites 
envb_NMDS <- metaMDS(envb_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Mixed and stratified lakes
envbms_NMDS <- metaMDS(envbms_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
envbom_NMDS <- metaMDS(envbom_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
envbso_NMDS <- metaMDS(envbso_dist, try = 1000, parallel = 4, previous.best, trymax = 500, maxit = 500)
```

## Determine stratification homogeneity using environmental/biogeographical distances

We use betadisper to determine homogeneity of the
environmental/biogeographical distances based on the stratification
category. Is the dispersion of environmental/biogeographical distances
similar within stratification categories?

``` r
### Environmental
groupse <- env[-c(9,16,18,20),19]
enve_bd <- betadisper(enve_dist, groupse)
anova(enve_bd)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value Pr(>F)  
    ## Groups     2 2.6595  1.3298  6.0582  0.011 *
    ## Residuals 16 3.5120  0.2195                 
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
    ## Groups     2 2.6595  1.3298 6.0582    999  0.009 **
    ## Residuals 16 3.5120  0.2195                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##               Mixed    Ocean Stratified
    ## Mixed               0.139000      0.029
    ## Ocean      0.146462               0.016
    ## Stratified 0.034379 0.013298

``` r
(enve_bd.HSD <- TukeyHSD(enve_bd))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                        diff         lwr       upr     p adj
    ## Ocean-Mixed      -0.4426265 -1.26105575 0.3758028 0.3665113
    ## Stratified-Mixed  0.5734890 -0.03096109 1.1779390 0.0643192
    ## Stratified-Ocean  1.0161155  0.19768622 1.8345447 0.0144650

``` r
boxplot(enve_bd)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Determine%20stratification%20homogeneity%20using%20environmental/biogeographical%20distances-1.png)<!-- -->

``` r
### Biogeographic
groupsb <- env[-c(8,20),19]
envb_bd <- betadisper(envb_dist, groupsb)
anova(envb_bd)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## Groups     2  9.293  4.6467  2.6437 0.09848 .
    ## Residuals 18 31.638  1.7577                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
permutest(envb_bd, pairwise = TRUE)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
    ## Groups     2  9.293  4.6467 2.6437    999  0.062 .
    ## Residuals 18 31.638  1.7577                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              Mixed   Ocean Stratified
    ## Mixed              0.06700      0.695
    ## Ocean      0.12207              0.081
    ## Stratified 0.53863 0.12721

``` r
(envb_bd.HSD <- TukeyHSD(envb_bd))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                      diff        lwr       upr     p adj
    ## Ocean-Mixed       1.55171 -0.3307344 3.4341546 0.1171621
    ## Stratified-Mixed  0.16319 -1.5879737 1.9143538 0.9693514
    ## Stratified-Ocean -1.38852 -3.2158566 0.4388165 0.1565488

``` r
boxplot(envb_bd)
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Determine%20stratification%20homogeneity%20using%20environmental/biogeographical%20distances-2.png)<!-- -->

## Determine stratification dissimilarity using environmental/biogeographical distances

We use adonis to determine if dissimilarities of species
environmental/biogeographical distances by stratification categories are
significant and how much of the variation is explained by
environmental/biogeographical dissimilarities.

``` r
### Environmental
enve_pmc <- adonis2(enve_dist ~ env[-c(9,16,18,20),19], permutations = 999)
enve_pmc
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = enve_dist ~ env[-c(9, 16, 18, 20), 19], permutations = 999)
    ##                            Df SumOfSqs      R2      F Pr(>F)    
    ## env[-c(9, 16, 18, 20), 19]  2   35.341 0.65446 15.152  0.001 ***
    ## Residual                   16   18.659 0.34554                  
    ## Total                      18   54.000 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratifcation
enve_pmc_pair <- pairwise.adonis(enve_dist, env[-c(9,16,18,20),19], p.adjust.m = "bonferroni", perm = 999)
enve_pmc_pair
```

    ##                 pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 24.969657 18.950602 0.5751216   0.003      0.009   *
    ## 2 Stratified vs Ocean  1 23.134197 14.727147 0.6206877   0.003      0.009   *
    ## 3      Mixed vs Ocean  1  1.589454  3.021701 0.2513539   0.086      0.258

``` r
### Biogeographic
envb_pmc <- adonis2(envb_dist ~ env[-c(8,20),19], permutations = 999)
envb_pmc
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = envb_dist ~ env[-c(8, 20), 19], permutations = 999)
    ##                    Df SumOfSqs      R2      F Pr(>F)    
    ## env[-c(8, 20), 19]  2   48.891 0.39926 5.9816  0.001 ***
    ## Residual           18   73.562 0.60074                  
    ## Total              20  122.454 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# By stratifcation
envb_pmc_pair <- pairwise.adonis(envb_dist, env[-c(8,20),19], p.adjust.m = "bonferroni", perm = 999)
envb_pmc_pair
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1  10.20361 7.086195 0.3527893   0.006      0.018   .
    ## 2 Stratified vs Ocean  1  40.90496 7.473867 0.3837896   0.001      0.003   *
    ## 3      Mixed vs Ocean  1  23.15310 4.060068 0.2695916   0.009      0.027   .

## Envfit trait influence on environmental/biogeographical distances

We use envfit to determine significantly correlated functional variables
to our environmental/biogeographical NMDS. We will use these results, if
significant, in our figure.

``` r
### Environmental
# All sites for figure
enve_ef <- envfit(enve_NMDS, FD_total_env[-c(9,16,18,20),c(54:55)], permutations = 1000, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
enve_ef
```

    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##        NMDS1   NMDS2
    ## B2s  -0.0712 -0.8103
    ## B3f   0.6100  0.2397
    ## B4e  -1.9470 -0.3737
    ## Ono  -1.9370 -0.6037
    ## Oyes  0.5165  0.1610
    ## 
    ## Goodness of fit:
    ##       r2   Pr(>r)   
    ## B 0.4513 0.006993 **
    ## O 0.4159 0.030969 * 
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
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##        NMDS1   NMDS2
    ## B2s  -0.0712 -0.8103
    ## B3f   0.6100  0.2397
    ## B4e  -1.9470 -0.3737
    ## Ono  -1.9370 -0.6037
    ## Oyes  0.5165  0.1610
    ## 
    ## Goodness of fit:
    ##       r2  Pr(>r)  
    ## B 0.4513 0.01399 *
    ## O 0.4159 0.06194 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# All sites by stratification
envea_ef <- envfit(enve_NMDS, FD_total_env[-c(9,16,18,20),c(6:22)], permutations = 999, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
envea_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)  
    ## MaxLengthTL       0.38726  0.92197 0.2545  0.160  
    ## Troph            -0.98542 -0.17016 0.6098  0.138  
    ## DepthMin          0.11233 -0.99367 0.0343  0.770  
    ## DepthMax          0.55073  0.83469 0.4131  0.101  
    ## TempPrefMin       0.72291 -0.69095 0.1232  0.879  
    ## TempPrefMax      -0.37699 -0.92622 0.1176  0.386  
    ## Weight            0.33233  0.94316 0.1191  0.268  
    ## CaudalFinLength   0.08773  0.99614 0.0656  0.529  
    ## DorsalSpinesMean  0.85784  0.51391 0.6120  0.076 .
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
    ## BodyShapeI2s        -0.0712 -0.8103
    ## BodyShapeI3f         0.6100  0.2397
    ## BodyShapeI4e        -1.9470 -0.3737
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.9370 -0.6037
    ## OperculumPresentyes  0.5165  0.1610
    ## FeedingPathb         0.2063  0.0528
    ## FeedingPathp        -1.7534 -0.4488
    ## RepGuild11b         -1.7534 -0.4488
    ## RepGuild12g         -2.1183 -0.4960
    ## RepGuild13n          0.5162  0.1260
    ## RepGuild22eb        -1.6415 -1.4438
    ## RepGuild23n         -2.1183 -0.4960
    ## RepGuild26s          0.3674  0.1522
    ## ParentalCare3p      -1.5812 -0.2305
    ## ParentalCare4n       0.7298  0.1064
    ## WaterPref1s          0.6173  0.2224
    ## WaterPref3a         -1.7286 -0.6226
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)   
    ## BodyShapeI       0.4513  0.008 **
    ## DemersPelag      0.0000  1.000   
    ## OperculumPresent 0.4159  0.014 * 
    ## FeedingPath      0.1460  0.379   
    ## RepGuild1        0.4039  0.098 . 
    ## RepGuild2        0.3345  0.078 . 
    ## ParentalCare     0.4465  0.323   
    ## WaterPref        0.4568  0.115   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ## MaxLengthTL       0.38726  0.92197 0.2545      1
    ## Troph            -0.98542 -0.17016 0.6098      1
    ## DepthMin          0.11233 -0.99367 0.0343      1
    ## DepthMax          0.55073  0.83469 0.4131      1
    ## TempPrefMin       0.72291 -0.69095 0.1232      1
    ## TempPrefMax      -0.37699 -0.92622 0.1176      1
    ## Weight            0.33233  0.94316 0.1191      1
    ## CaudalFinLength   0.08773  0.99614 0.0656      1
    ## DorsalSpinesMean  0.85784  0.51391 0.6120      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -0.0712 -0.8103
    ## BodyShapeI3f         0.6100  0.2397
    ## BodyShapeI4e        -1.9470 -0.3737
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.9370 -0.6037
    ## OperculumPresentyes  0.5165  0.1610
    ## FeedingPathb         0.2063  0.0528
    ## FeedingPathp        -1.7534 -0.4488
    ## RepGuild11b         -1.7534 -0.4488
    ## RepGuild12g         -2.1183 -0.4960
    ## RepGuild13n          0.5162  0.1260
    ## RepGuild22eb        -1.6415 -1.4438
    ## RepGuild23n         -2.1183 -0.4960
    ## RepGuild26s          0.3674  0.1522
    ## ParentalCare3p      -1.5812 -0.2305
    ## ParentalCare4n       0.7298  0.1064
    ## WaterPref1s          0.6173  0.2224
    ## WaterPref3a         -1.7286 -0.6226
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.4513  0.136
    ## DemersPelag      0.0000  1.000
    ## OperculumPresent 0.4159  0.238
    ## FeedingPath      0.1460  1.000
    ## RepGuild1        0.4039  1.000
    ## RepGuild2        0.3345  1.000
    ## ParentalCare     0.4465  1.000
    ## WaterPref        0.4568  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
envems_ef <- envfit(envems_NMDS, FD_total_env[c(1:6,8,10,12:15,17,21:23),c(6,8,10,19)], permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
envems_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##          NMDS1    NMDS2    r2 Pr(>r)  
    ## Troph -0.86893 -0.49494 0.608  0.039 *
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
    ## BodyShapeI2s        -1.4705 -1.3732
    ## BodyShapeI3f         0.7485  0.2253
    ## BodyShapeI4e        -1.6908 -0.2764
    ## OperculumPresentno  -1.6936 -0.5096
    ## OperculumPresentyes  0.5645  0.1699
    ## RepGuild22eb        -1.4705 -1.3732
    ## RepGuild23n         -1.8914 -0.3776
    ## RepGuild26s          0.4041  0.1637
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## BodyShapeI       0.5560  0.011 *
    ## OperculumPresent 0.4121  0.029 *
    ## RepGuild2        0.3448  0.045 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ##          NMDS1    NMDS2    r2 Pr(>r)
    ## Troph -0.86893 -0.49494 0.608  0.156
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -1.4705 -1.3732
    ## BodyShapeI3f         0.7485  0.2253
    ## BodyShapeI4e        -1.6908 -0.2764
    ## OperculumPresentno  -1.6936 -0.5096
    ## OperculumPresentyes  0.5645  0.1699
    ## RepGuild22eb        -1.4705 -1.3732
    ## RepGuild23n         -1.8914 -0.3776
    ## RepGuild26s          0.4041  0.1637
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## BodyShapeI       0.5560  0.044 *
    ## OperculumPresent 0.4121  0.116  
    ## RepGuild2        0.3448  0.180  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
enveom_ef <- envfit(enveom_NMDS, FD_total_env[c(3,6:8,10:11,13:15,19,23),c(6:22)], permutations = 999, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
enveom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL      -0.88843  0.45902 0.2417  0.429
    ## Troph            -0.75114  0.66014 0.5404  0.190
    ## DepthMin          0.80672 -0.59093 0.2325  0.309
    ## DepthMax         -0.99176  0.12808 0.0417  0.853
    ## TempPrefMin       0.17837 -0.98396 0.0508  0.792
    ## TempPrefMax       0.00397  0.99999 0.0211  0.920
    ## Weight           -0.18882  0.98201 0.1201  0.534
    ## CaudalFinLength  -0.97758 -0.21058 0.3166  0.276
    ## DorsalSpinesMean  0.86406 -0.50339 0.1736  0.549
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.4375 -0.0563
    ## BodyShapeI3f        -0.0438  0.0056
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
    ## BodyShapeI       0.0351      1
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
    ## MaxLengthTL      -0.88843  0.45902 0.2417      1
    ## Troph            -0.75114  0.66014 0.5404      1
    ## DepthMin          0.80672 -0.59093 0.2325      1
    ## DepthMax         -0.99176  0.12808 0.0417      1
    ## TempPrefMin       0.17837 -0.98396 0.0508      1
    ## TempPrefMax       0.00397  0.99999 0.0211      1
    ## Weight           -0.18882  0.98201 0.1201      1
    ## CaudalFinLength  -0.97758 -0.21058 0.3166      1
    ## DorsalSpinesMean  0.86406 -0.50339 0.1736      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.4375 -0.0563
    ## BodyShapeI3f        -0.0438  0.0056
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
    ## BodyShapeI       0.0351      1
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
enveso_ef <- envfit(enveso_NMDS, FD_total_env[c(1:2,4,5,7,11:12,17,19,21:22),c(6,8,19)], permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,11:12,17,19,21:22),19])
enveso_ef
```

    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.6527 -0.6941
    ## BodyShapeI3f         0.7703  0.5322
    ## BodyShapeI4e        -1.2893 -0.3182
    ## OperculumPresentno  -1.2190 -0.5445
    ## OperculumPresentyes  0.6966  0.3111
    ## RepGuild22eb        -0.8320 -1.3628
    ## RepGuild23n         -1.5060 -0.4656
    ## RepGuild26s          0.4805  0.2868
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## BodyShapeI       0.3964  0.027 *
    ## OperculumPresent 0.3351  0.020 *
    ## RepGuild2        0.2998  0.038 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.6527 -0.6941
    ## BodyShapeI3f         0.7703  0.5322
    ## BodyShapeI4e        -1.2893 -0.3182
    ## OperculumPresentno  -1.2190 -0.5445
    ## OperculumPresentyes  0.6966  0.3111
    ## RepGuild22eb        -0.8320 -1.3628
    ## RepGuild23n         -1.5060 -0.4656
    ## RepGuild26s          0.4805  0.2868
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## BodyShapeI       0.3964  0.081 .
    ## OperculumPresent 0.3351  0.060 .
    ## RepGuild2        0.2998  0.114  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Biogeographic
# All sites for figure
envb_ef <- envfit(envb_NMDS, FD_total_env[-c(8,20),c(58,67)], permutations = 999, na.rm = TRUE, strata = env[-c(8,20),19])
envb_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##      NMDS1    NMDS2     r2 Pr(>r)
    ## T -0.44636 -0.89485 0.6512  0.158
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##       NMDS1   NMDS2
    ## W1s  0.3441  0.3821
    ## W3a -1.1012 -1.2227
    ## 
    ## Goodness of fit:
    ##       r2 Pr(>r)  
    ## W 0.1855  0.024 *
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
    ## ***VECTORS
    ## 
    ##      NMDS1    NMDS2     r2 Pr(>r)
    ## T -0.44636 -0.89485 0.6512  0.316
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##       NMDS1   NMDS2
    ## W1s  0.3441  0.3821
    ## W3a -1.1012 -1.2227
    ## 
    ## Goodness of fit:
    ##       r2 Pr(>r)  
    ## W 0.1855  0.048 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# All sites by stratification
envba_ef <- envfit(envb_NMDS, FD_total_env[-c(8,20),c(6:22)], permutations = 999, na.rm = TRUE, strata = env[-c(8,20),19])
envba_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.18650  0.98246 0.0236  0.640
    ## Troph            -0.44636 -0.89485 0.6512  0.143
    ## DepthMin          0.26612  0.96394 0.0213  0.600
    ## DepthMax          0.41126  0.91152 0.1295  0.316
    ## TempPrefMin       0.17523  0.98453 0.1889  0.447
    ## TempPrefMax      -0.97425 -0.22545 0.0225  0.517
    ## Weight           -0.96923  0.24618 0.0132  0.823
    ## CaudalFinLength  -0.62812 -0.77812 0.0030  0.969
    ## DorsalSpinesMean  0.53852  0.84261 0.3114  0.605
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.0089  0.1071
    ## BodyShapeI3f         0.2695  0.2241
    ## BodyShapeI4e        -0.9500 -0.8646
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.0607 -1.2425
    ## OperculumPresentyes  0.2496  0.2924
    ## FeedingPathb         0.0779  0.1295
    ## FeedingPathp        -0.7403 -1.2303
    ## RepGuild11b         -0.7403 -1.2303
    ## RepGuild12g         -1.2944 -0.4556
    ## RepGuild13n          0.2394  0.1983
    ## RepGuild22eb        -1.1859 -1.8665
    ## RepGuild23n         -1.2944 -0.4556
    ## RepGuild26s          0.2097  0.1543
    ## ParentalCare3p      -0.8603 -1.0218
    ## ParentalCare4n       0.3441  0.4087
    ## WaterPref1s          0.3441  0.3821
    ## WaterPref3a         -1.1012 -1.2227
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## BodyShapeI       0.0872  0.781  
    ## DemersPelag      0.0000  1.000  
    ## OperculumPresent 0.1377  0.100 .
    ## FeedingPath      0.0476  0.355  
    ## RepGuild1        0.0995  0.425  
    ## RepGuild2        0.1031  0.141  
    ## ParentalCare     0.1565  0.293  
    ## WaterPref        0.1855  0.016 *
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
    ## MaxLengthTL       0.18650  0.98246 0.0236      1
    ## Troph            -0.44636 -0.89485 0.6512      1
    ## DepthMin          0.26612  0.96394 0.0213      1
    ## DepthMax          0.41126  0.91152 0.1295      1
    ## TempPrefMin       0.17523  0.98453 0.1889      1
    ## TempPrefMax      -0.97425 -0.22545 0.0225      1
    ## Weight           -0.96923  0.24618 0.0132      1
    ## CaudalFinLength  -0.62812 -0.77812 0.0030      1
    ## DorsalSpinesMean  0.53852  0.84261 0.3114      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s         0.0089  0.1071
    ## BodyShapeI3f         0.2695  0.2241
    ## BodyShapeI4e        -0.9500 -0.8646
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.0607 -1.2425
    ## OperculumPresentyes  0.2496  0.2924
    ## FeedingPathb         0.0779  0.1295
    ## FeedingPathp        -0.7403 -1.2303
    ## RepGuild11b         -0.7403 -1.2303
    ## RepGuild12g         -1.2944 -0.4556
    ## RepGuild13n          0.2394  0.1983
    ## RepGuild22eb        -1.1859 -1.8665
    ## RepGuild23n         -1.2944 -0.4556
    ## RepGuild26s          0.2097  0.1543
    ## ParentalCare3p      -0.8603 -1.0218
    ## ParentalCare4n       0.3441  0.4087
    ## WaterPref1s          0.3441  0.3821
    ## WaterPref3a         -1.1012 -1.2227
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.0872  1.000
    ## DemersPelag      0.0000  1.000
    ## OperculumPresent 0.1377  1.000
    ## FeedingPath      0.0476  1.000
    ## RepGuild1        0.0995  1.000
    ## RepGuild2        0.1031  1.000
    ## ParentalCare     0.1565  1.000
    ## WaterPref        0.1855  0.272
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
envbms_ef <- envfit(envbms_NMDS, FD_total_env[c(1:6,10,12:15,17,21:23),c(6:22)], permutations = 999, na.rm = TRUE, strata = env[c(1:6,10,12:15,17,21:23),19])
envbms_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                     NMDS1    NMDS2     r2 Pr(>r)
    ## MaxLengthTL       0.57614  0.81735 0.0924  0.596
    ## Troph            -0.84701 -0.53158 0.4169  0.173
    ## DepthMin          0.73775  0.67507 0.0336  0.819
    ## DepthMax          0.62925  0.77720 0.2900  0.216
    ## TempPrefMin       0.94201  0.33559 0.0941  0.653
    ## TempPrefMax      -0.27317 -0.96197 0.0783  0.561
    ## Weight            0.10132  0.99485 0.0437  0.710
    ## CaudalFinLength   0.53417  0.84537 0.0180  0.874
    ## DorsalSpinesMean  0.66895  0.74331 0.1653  0.659
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -1.9564 -0.2898
    ## BodyShapeI3f         0.4935  0.0919
    ## BodyShapeI4e        -0.7447 -0.1573
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.2058 -0.2432
    ## OperculumPresentyes  0.4385  0.0884
    ## FeedingPathb         0.1780 -0.0238
    ## FeedingPathp        -1.1571  0.1545
    ## RepGuild11b         -1.1571  0.1545
    ## RepGuild12g         -0.2955 -0.5918
    ## RepGuild13n          0.2641  0.0795
    ## RepGuild22eb        -1.9564 -0.2898
    ## RepGuild23n         -0.2955 -0.5918
    ## RepGuild26s          0.2123  0.1228
    ## ParentalCare3p      -0.9083 -0.0211
    ## ParentalCare4n       0.6056  0.0141
    ## WaterPref1s          0.5920  0.1443
    ## WaterPref3a         -1.1839 -0.2887
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## BodyShapeI       0.3495  0.240  
    ## DemersPelag      0.0000  1.000  
    ## OperculumPresent 0.3297  0.108  
    ## FeedingPath      0.1256  0.301  
    ## RepGuild1        0.1773  0.529  
    ## RepGuild2        0.2200  0.292  
    ## ParentalCare     0.3298  0.322  
    ## WaterPref        0.4449  0.073 .
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
    ## MaxLengthTL       0.57614  0.81735 0.0924      1
    ## Troph            -0.84701 -0.53158 0.4169      1
    ## DepthMin          0.73775  0.67507 0.0336      1
    ## DepthMax          0.62925  0.77720 0.2900      1
    ## TempPrefMin       0.94201  0.33559 0.0941      1
    ## TempPrefMax      -0.27317 -0.96197 0.0783      1
    ## Weight            0.10132  0.99485 0.0437      1
    ## CaudalFinLength   0.53417  0.84537 0.0180      1
    ## DorsalSpinesMean  0.66895  0.74331 0.1653      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## BodyShapeI2s        -1.9564 -0.2898
    ## BodyShapeI3f         0.4935  0.0919
    ## BodyShapeI4e        -0.7447 -0.1573
    ## DemersPelag1r        0.0000  0.0000
    ## OperculumPresentno  -1.2058 -0.2432
    ## OperculumPresentyes  0.4385  0.0884
    ## FeedingPathb         0.1780 -0.0238
    ## FeedingPathp        -1.1571  0.1545
    ## RepGuild11b         -1.1571  0.1545
    ## RepGuild12g         -0.2955 -0.5918
    ## RepGuild13n          0.2641  0.0795
    ## RepGuild22eb        -1.9564 -0.2898
    ## RepGuild23n         -0.2955 -0.5918
    ## RepGuild26s          0.2123  0.1228
    ## ParentalCare3p      -0.9083 -0.0211
    ## ParentalCare4n       0.6056  0.0141
    ## WaterPref1s          0.5920  0.1443
    ## WaterPref3a         -1.1839 -0.2887
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)
    ## BodyShapeI       0.3495      1
    ## DemersPelag      0.0000      1
    ## OperculumPresent 0.3297      1
    ## FeedingPath      0.1256      1
    ## RepGuild1        0.1773      1
    ## RepGuild2        0.2200      1
    ## ParentalCare     0.3298      1
    ## WaterPref        0.4449      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
envbom_ef <- envfit(envbom_NMDS, FD_total_env[c(3,6:7,9:11,13:16,18:19,23),c(13:14)], permutations = 999, na.rm = TRUE, strata = env[c(3,6:7,9:11,13:16,18:19,23),19])
envbom_ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                 NMDS1     NMDS2     r2 Pr(>r)   
    ## TempPrefMin -0.050506 -0.998720 0.6062  0.006 **
    ## TempPrefMax -0.092281 -0.995730 0.4263  0.013 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ##                 NMDS1     NMDS2     r2 Pr(>r)  
    ## TempPrefMin -0.050506 -0.998720 0.6062  0.012 *
    ## TempPrefMax -0.092281 -0.995730 0.4263  0.026 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
envbso_ef <- envfit(envbso_NMDS, FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(8,21)], permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
envbso_ef
```

    ## 
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## OperculumPresentno  -1.4615 -0.9116
    ## OperculumPresentyes  0.5846  0.3647
    ## WaterPref1s          0.8538  0.4932
    ## WaterPref3a         -1.5368 -0.8878
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## OperculumPresent 0.1931  0.083 .
    ## WaterPref        0.2848  0.020 *
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
    ## ***FACTORS:
    ## 
    ## Centroids:
    ##                       NMDS1   NMDS2
    ## OperculumPresentno  -1.4615 -0.9116
    ## OperculumPresentyes  0.5846  0.3647
    ## WaterPref1s          0.8538  0.4932
    ## WaterPref3a         -1.5368 -0.8878
    ## 
    ## Goodness of fit:
    ##                      r2 Pr(>r)  
    ## OperculumPresent 0.1931  0.166  
    ## WaterPref        0.2848  0.040 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

## Environment/biogeography mantel tests

We used mantel tests to determine significance between the
environmental/biogeographical and trait distance matrices.

``` r
### Environmental
# All sites
FDe_dist <- gowdis(scaled_FD_total_env[-c(9,16,18,20),c(1:17)])
envea_mant <- mantel(enve_dist, FDe_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(9,16,18,20),19])
envea_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = enve_dist, ydis = FDe_dist, method = "spearman",      permutations = 999, strata = env[-c(9, 16, 18, 20), 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5905 
    ##       Significance: 0.093 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.585 0.615 0.629 0.648 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
FDems_dist <- gowdis(scaled_FD_total_env[c(1:6,8,10,12:15,17,21:23),c(1:17)])
envems_mant <- mantel(envems_dist, FDems_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:6,8,10,12:15,17,21:23),19])
envems_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = envems_dist, ydis = FDems_dist, method = "spearman",      permutations = 999, strata = env[c(1:6, 8, 10, 12:15, 17,          21:23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5362 
    ##       Significance: 0.078 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.521 0.558 0.577 0.599 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
FDeom_dist <- gowdis(scaled_FD_total_env[c(3,6:8,10:11,13:15,19,23),c(1:7)])
enveom_mant <- mantel(enveom_dist, FDeom_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(3,6:8,10:11,13:15,19,23),19])
enveom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = enveom_dist, ydis = FDeom_dist, method = "spearman",      permutations = 999, strata = env[c(3, 6:8, 10:11, 13:15,          19, 23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3447 
    ##       Significance: 0.094 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.328 0.445 0.525 0.566 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
FDeso_dist <- gowdis(scaled_FD_total_env[c(1:2,4,5,7,11:12,17,19,21:22),c(1:17)])
enveso_mant <- mantel(enveso_dist, FDeso_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[c(1:2,4,5,7,11:12,17,19,21:22),19])
enveso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = enveso_dist, ydis = FDeso_dist, method = "spearman",      permutations = 999, strata = env[c(1:2, 4, 5, 7, 11:12, 17,          19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2766 
    ##       Significance: 0.067 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.249 0.297 0.327 0.364 
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

    ## [1] 0.234 0.282 0.201

``` r
### Biogeographic
# All sites for figure
FDb_dist <- gowdis(scaled_FD_total_env[-c(8,20),c(1:17)])
envba_mant <- mantel(envb_dist, FDb_dist, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[-c(8,20),19])
envba_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = envb_dist, ydis = FDb_dist, method = "spearman",      permutations = 999, strata = env[-c(8, 20), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.08245 
    ##       Significance: 0.337 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.159 0.193 0.226 0.251 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
FDbms_dist <- gowdis(scaled_FD_total_env[c(1:6,10,12:15,17,21:23),c(1:17)])
envbms_mant <- mantel(envbms_dist, FDbms_dist, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:6,10,12:15,17,21:23),19])
envbms_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = envbms_dist, ydis = FDbms_dist, method = "spearman",      permutations = 1000, strata = env[c(1:6, 10, 12:15, 17, 21:23),          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1734 
    ##       Significance: 0.34366 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.340 0.413 0.441 0.483 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Ocean sites and mixed lakes
FDbom_dist <- gowdis(scaled_FD_total_env[c(3,6:7,9:11,13:16,18:19,23),c(1:7)])
envbom_mant <- mantel(envbom_dist, FDbom_dist, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(3,6:7,9:11,13:16,18:19,23),19])
envbom_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = envbom_dist, ydis = FDbom_dist, method = "spearman",      permutations = 1000, strata = env[c(3, 6:7, 9:11, 13:16,          18:19, 23), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.06714 
    ##       Significance: 0.12887 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0834 0.1254 0.1621 0.1890 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1000

``` r
# Stratified lakes and ocean sites
FDbso_dist <- gowdis(scaled_FD_total_env[c(1:2,4,5,7,9,11:12,16:19,21:22),c(1:17)])
envbso_mant <- mantel(envbso_dist, FDbso_dist, method = "spearman", permutations = 1000, na.rm = TRUE, strata = env[c(1:2,4,5,7,9,11:12,16:19,21:22),19])
envbso_mant
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = envbso_dist, ydis = FDbso_dist, method = "spearman",      permutations = 1000, strata = env[c(1:2, 4, 5, 7, 9, 11:12,          16:19, 21:22), 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.05808 
    ##       Significance: 0.30669 
    ## 
    ## Upper quantiles of permutations (null model):
    ##     90%     95%   97.5%     99% 
    ## 0.00961 0.03714 0.06229 0.08660 
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

    ## [1] 1.0000000 0.3866134 0.9200799

## Plot environmental/biogeographical NMDS and envfit results

Includes a plot of correlated trait variables

``` r
enve_NMDS_data.scores <- as.data.frame(scores(enve_NMDS))
enve_NMDS_data.scores$Stratification <- env[-c(9,16,18,20),19]
enve_NMDS_data.scores$Lakes <- env[-c(9,16,18,20),1]
enve_NMDS_data.scores$Stratification <- factor(enve_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# enve_ef_coord_cont <- as.data.frame(scores(enve_ef, "vectors")) * ordiArrowMul(enve_ef)
# enve_ef_coord_cat = as.data.frame(scores(enve_ef, "factors")) * ordiArrowMul(enve_ef)

enve_ef_plot <- ggplot(data = enve_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification),level = 0.50) +
  geom_point(data = enve_NMDS_data.scores, aes(color = Stratification), size = 4, alpha = 1) + 
  geom_text_repel(data = enve_NMDS_data.scores, label = enve_NMDS_data.scores$Lakes, size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              data = enve_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#CD3333") +
  # geom_text(data = enve_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
  #           color = "#CD3333", label = row.names(enve_ef_coord_cont), size = 7) + 
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              data = enve_ef_coord_cat, linewidth =1, alpha = 0.2, color = "#8B6508") +
  # geom_text(data = enve_ef_coord_cat, aes(x = NMDS1, y = NMDS2), 
  #           color = "#8B6508", label = row.names(enve_ef_coord_cat), size = 7, check_overlap = T) + 
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
        panel.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  # annotate(geom = "label", x = 1.35, y = 1.55, size = 6, 
  #          label = paste("Stress: ", round(enve_NMDS$stress, digits = 2))) +
  labs(color = "Stratification", tag = "A")
enve_ef_plot <- enve_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
enve_ef_plot
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values (`geom_path()`).

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20environmental/biogeographical%20NMDS%20and%20envfit%20results-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/enve_ef_plot.png", enve_ef_plot, width = 8, height = 4, units = "in")
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values (`geom_path()`).

    ## Warning: ggrepel: 1 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
envb_NMDS_data.scores <- as.data.frame(scores(envb_NMDS))
envb_NMDS_data.scores$Stratification <- env[-c(8,20),19]
envb_NMDS_data.scores$Lakes <- env[-c(8,20),1]
envb_NMDS_data.scores$Stratification <- factor(envb_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# envb_ef_coord_cont <- as.data.frame(scores(envb_ef, "vectors")) * ordiArrowMul(envb_ef)
# envb_ef_coord_cat = as.data.frame(scores(envb_ef, "factors")) * ordiArrowMul(envb_ef)

envb_ef_plot <- ggplot(data = envb_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification),level = 0.50) +
  geom_point(data = envb_NMDS_data.scores, aes(color = Stratification), size = 3, alpha = 1) + 
  geom_text_repel(data = envb_NMDS_data.scores, label = envb_NMDS_data.scores$Lakes, size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              data = envb_ef_coord_cont, linewidth =1, alpha = 0.2, color = "#CD3333") +
  # geom_text(data = envb_ef_coord_cont, aes(x = NMDS1, y = NMDS2), 
  #           color = "#CD3333", label = row.names(envb_ef_coord_cont), size = 7) + 
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  #              data = envb_ef_coord_cat, linewidth =1, alpha = 0.2, color = "#8B6508") +
  # geom_text(data = envb_ef_coord_cat, aes(x = NMDS1, y = NMDS2), 
  #           color = "#8B6508", label = row.names(envb_ef_coord_cat), size = 7, check_overlap = T) + 
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
        panel.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank()) +
  # annotate(geom = "label", x = 5.9, y = 1.59, size = 6, 
  #          label = paste("Stress: ", round(envb_NMDS$stress, digits = 2))) +
  labs(color = "Stratification", tag = "B")
envb_ef_plot <- envb_ef_plot + guides(color = guide_legend(override.aes = list(label = "")))
envb_ef_plot
```

![](stratification_NMDS_envfit_analyses_files/figure-gfm/Plot%20environmental/biogeographical%20NMDS%20and%20envfit%20results-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_NMDS_envfit_plots/envb_ef_plot.png", envb_ef_plot, width = 8, height = 4, units = "in")
```

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.4.1
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
    ##  [1] adespatial_0.3-23    picante_1.8.2        nlme_3.1-164        
    ##  [4] pairwiseAdonis_0.4.1 cluster_2.1.6        FD_1.0-12.3         
    ##  [7] geometry_0.4.7       ade4_1.7-22          vegan_2.6-4         
    ## [10] lattice_0.22-5       permute_0.9-7        ggrepel_0.9.4       
    ## [13] viridis_0.6.4        viridisLite_0.4.2    ggplot2_3.4.4       
    ## [16] tidyr_1.3.0          phytools_2.0-3       maps_3.4.1.1        
    ## [19] ape_5.7-1            reshape2_1.4.4       stringr_1.5.1       
    ## [22] dplyr_1.1.4          knitr_1.45          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      wk_0.9.1                rstudioapi_0.15.0      
    ##   [4] magrittr_2.0.3          farver_2.1.1            rmarkdown_2.25         
    ##   [7] adegraphics_1.0-21      ragg_1.2.6              vctrs_0.6.5            
    ##  [10] spdep_1.3-1             htmltools_0.5.7         progress_1.2.3         
    ##  [13] s2_1.1.5                adegenet_2.1.10         spData_2.3.0           
    ##  [16] KernSmooth_2.23-22      adephylo_1.1-16         plyr_1.8.9             
    ##  [19] uuid_1.1-1              igraph_1.5.1            mime_0.12              
    ##  [22] lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3        
    ##  [25] Matrix_1.6-4            R6_2.5.1                fastmap_1.1.1          
    ##  [28] shiny_1.8.0             magic_1.6-1             digest_0.6.33          
    ##  [31] numDeriv_2016.8-1.1     colorspace_2.1-0        phylobase_0.8.10       
    ##  [34] textshaping_0.3.7       labeling_0.4.3          clusterGeneration_1.3.8
    ##  [37] fansi_1.0.6             httr_1.4.7              abind_1.4-5            
    ##  [40] mgcv_1.9-0              compiler_4.3.1          proxy_0.4-27           
    ##  [43] withr_2.5.2             doParallel_1.0.17       optimParallel_1.0-2    
    ##  [46] DBI_1.1.3               highr_0.10              MASS_7.3-60            
    ##  [49] classInt_0.4-10         scatterplot3d_0.3-44    tools_4.3.1            
    ##  [52] units_0.8-5             rncl_0.8.7              httpuv_1.6.13          
    ##  [55] glue_1.6.2              quadprog_1.5-8          promises_1.2.1         
    ##  [58] grid_4.3.1              sf_1.0-14               generics_0.1.3         
    ##  [61] seqinr_4.2-36           gtable_0.3.4            class_7.3-22           
    ##  [64] hms_1.1.3               sp_2.1-2                xml2_1.3.6             
    ##  [67] utf8_1.2.4              foreach_1.5.2           pillar_1.9.0           
    ##  [70] later_1.3.2             splines_4.3.1           deldir_2.0-2           
    ##  [73] tidyselect_1.2.0        gridExtra_2.3           xfun_0.41              
    ##  [76] expm_0.999-8            stringi_1.8.2           yaml_2.3.7             
    ##  [79] boot_1.3-28.1           evaluate_0.23           codetools_0.2-19       
    ##  [82] interp_1.1-5            tibble_3.2.1            cli_3.6.1              
    ##  [85] xtable_1.8-4            systemfonts_1.0.5       munsell_0.5.0          
    ##  [88] Rcpp_1.0.11             coda_0.19-4             png_0.1-8              
    ##  [91] XML_3.99-0.16           ellipsis_0.3.2          RNeXML_2.4.11          
    ##  [94] prettyunits_1.2.0       latticeExtra_0.6-30     jpeg_0.1-10            
    ##  [97] phangorn_2.11.1         scales_1.3.0            e1071_1.7-14           
    ## [100] purrr_1.0.2             crayon_1.5.2            combinat_0.0-8         
    ## [103] rlang_1.1.2             fastmatch_1.1-4         mnormt_2.1.1

# Extra unused code

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
