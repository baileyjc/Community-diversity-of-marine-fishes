Functional diversity analyses
================

### Load packages and files

``` r
# Load the knitr package if not already loaded
library(knitr)

# Source the R Markdown file
knit("/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.Rmd", output = "/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md")
```

    ## 
    ## 
    ## processing file: /Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.Rmd

    ##   |                  |          |   0%  |                  |          |   3%                                                                          |                  |.         |   7% [Bringing everything together load modifying files packages]             |                  |.         |  10%                                                                          |                  |.         |  13% [Bringing everything together load in modifying files]                   |                  |..        |  17%                                                                          |                  |..        |  20% [Check species names across files]                                       |                  |..        |  23%                                                                          |                  |...       |  27% [Modify environment data]                                                |                  |...       |  30%                                                                          |                  |...       |  33% [Modify incidence matrices]                                              |                  |....      |  37%                                                                          |                  |....      |  40% [Modify phylogeny]                                                       |                  |....      |  43%                                                                          |                  |.....     |  47% [Modify trait data]                                                      |                  |.....     |  50%                                                                          |                  |.....     |  53% [Modify Site_type data frames]                                           |                  |......    |  57%                                                                          |                  |......    |  60% [Modify Site_type trait data]                                            |                  |......    |  63%                                                                          |                  |.......   |  67% [Site_type trait data tests]                                             |                  |.......   |  70%                                                                          |                  |.......   |  73% [Modify site trait data frames]                                          |                  |........  |  77%                                                                          |                  |........  |  80% [Site trait data tests]                                                  |                  |........  |  83%                                                                          |                  |......... |  87% [unnamed-chunk-7]                                                        |                  |......... |  90%                                                                          |                  |......... |  93% [unnamed-chunk-8]                                                        |                  |..........|  97%                                                                          |                  |..........| 100% [Bringing everything together load out modified files and session info]

    ## output file: /Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md

    ## [1] "/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md"

``` r
# Packages
library(picante)
```

    ## Loading required package: vegan

    ## Loading required package: permute

    ## Loading required package: lattice

    ## 
    ## Attaching package: 'vegan'

    ## The following object is masked from 'package:phytools':
    ## 
    ##     scores

    ## Loading required package: nlme

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

``` r
library(vegan)
library(phytools)
library(ggplot2)
library(ggrepel)
library(caret)
```

    ## 
    ## Attaching package: 'caret'

    ## The following object is masked from 'package:vegan':
    ## 
    ##     tolerance

``` r
library(BAT)
```

    ## 
    ## Attaching package: 'BAT'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     alpha

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     fill

    ## The following objects are masked from 'package:base':
    ## 
    ##     beta, gamma

``` r
library(pairwiseAdonis)
```

    ## Loading required package: cluster

    ## 
    ## Attaching package: 'cluster'

    ## The following object is masked from 'package:maps':
    ## 
    ##     votes.repub

``` r
library(emmeans)
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

``` r
# Site vectors
surveyed_sites <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LCN", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOM", "OOO", "OTM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_wo_LCN <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOM", "OOO", "OTM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_wo_TLN_HLM <- c("BCM", "CLM", "FLK", "GLK", "HLO", "IBK", "LCN", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOM", "OOO", "OTM", "RCA", "SLN", "ULN")

surveyed_sites_env <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OTM", "RCA", "SLN", "TLN", "ULN")

ocean_mixed_sites <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA", "FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")

ocean_mixed_sites_env <- c("IBK", "NCN", "RCA", "FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")

ocean_stratified_sites <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA", "BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")

ocean_stratified_sites_env <- c("IBK", "NCN", "RCA", "BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")

mixed_stratified_lakes <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "LLN", "MLN", "NLK", "NLN", "NLU", "OLO", "OTM", "SLN", "TLN", "ULN")

ocean_sites <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA")

ocean_sites_env <- c("IBK", "NCN", "RCA")

mixed_lakes <- c("FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")

stratified_lakes <- c("BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")

# Define your custom colors
custom_colors <- c("Reference" = "black", "Ocean" = "#EE6363", "Mixed" = "#87CEFA", "Stratified" = "#6E8B3D")

env_cont <- "#FFB90F"
```

### Edit input files

``` r
env <- env[order(env$X),]
presabs_lake <- presabs_lake[order(row.names(presabs_lake)),]
surveyed_sites_lake <- surveyed_sites_lake[order(row.names(surveyed_sites_lake)),]

Site_type_group_ref <- env[,19]
Site_type_group <- env[surveyed_sites,19]

straits_na <- which(rowSums(is.na(straits)) > 0)
straits_clean <- straits[-straits_na, ]
```

# Functional Diversity

### Identify traits with variation

``` r
## Reference
# Determine which traits have variance
nzv <- nearZeroVar(traits)
nzv
```

    ## integer(0)

``` r
# All traits have variation

traits$presence_sum <- rowSums(presabs_species)  # Sum occurrences across all sites
# Ensure traits_numeric is a data frame
mod <-  lm(presence_sum ~ BodyShapeI + DemersPelag + OperculumPresent + MaxLengthTL + Troph + DepthMin + DepthMax + TempPrefMin + TempPrefMax + FeedingPath + RepGuild2 + ParentalCare + WaterPref + DorsalSpinesMean, traits)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = presence_sum ~ BodyShapeI + DemersPelag + OperculumPresent + 
    ##     MaxLengthTL + Troph + DepthMin + DepthMax + TempPrefMin + 
    ##     TempPrefMax + FeedingPath + RepGuild2 + ParentalCare + WaterPref + 
    ##     DorsalSpinesMean, data = traits)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6137 -0.6768 -0.3531  0.0526 12.2027 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          2.6908970  1.6929838   1.589 0.112151    
    ## BodyShapeI2s         0.7134603  0.4509856   1.582 0.113839    
    ## BodyShapeI3f         0.7184052  0.4508056   1.594 0.111214    
    ## BodyShapeI4e         0.6939366  0.4520646   1.535 0.124964    
    ## BodyShapeI5l         0.5260132  0.4721222   1.114 0.265377    
    ## DemersPelag2pn      -0.2705386  0.2737268  -0.988 0.323123    
    ## DemersPelag3p       -0.7505488  1.2534547  -0.599 0.549398    
    ## DemersPelag4po      -0.1211062  0.3705341  -0.327 0.743828    
    ## DemersPelag5d       -0.4626085  0.1364634  -3.390 0.000715 ***
    ## DemersPelag6bp      -0.3574284  0.2538203  -1.408 0.159260    
    ## DemersPelag7bd       0.2451175  1.2234277   0.200 0.841229    
    ## OperculumPresentyes  0.3711498  0.0876401   4.235 2.41e-05 ***
    ## MaxLengthTL         -0.0000217  0.0007159  -0.030 0.975821    
    ## Troph               -0.1086047  0.0817110  -1.329 0.183986    
    ## DepthMin            -0.0036201  0.0034454  -1.051 0.293553    
    ## DepthMax            -0.0000868  0.0003638  -0.239 0.811471    
    ## TempPrefMin          0.0519679  0.0222752   2.333 0.019767 *  
    ## TempPrefMax         -0.0614448  0.0599027  -1.026 0.305161    
    ## FeedingPathp        -0.1307345  0.0925640  -1.412 0.158028    
    ## RepGuild22eb        -0.4328537  0.4713257  -0.918 0.358556    
    ## RepGuild23n         -0.6764647  0.4471733  -1.513 0.130531    
    ## RepGuild24t         -0.6146815  0.4810022  -1.278 0.201456    
    ## RepGuild25h         -0.9890426  1.1956561  -0.827 0.408244    
    ## RepGuild26s         -0.6392609  0.4987003  -1.282 0.200072    
    ## ParentalCare2m      -0.5531592  0.5453568  -1.014 0.310584    
    ## ParentalCare3p      -0.7564641  0.4608894  -1.641 0.100921    
    ## ParentalCare4n      -0.6356969  0.5166541  -1.230 0.218717    
    ## WaterPref2bs         0.2840775  0.1229165   2.311 0.020946 *  
    ## WaterPref3a          0.7526872  0.1862863   4.040 5.58e-05 ***
    ## WaterPref4b          0.0318851  1.2526535   0.025 0.979696    
    ## WaterPref5fb         0.1306558  0.5272456   0.248 0.804313    
    ## WaterPref6f          0.2861852  0.6483562   0.441 0.658980    
    ## DorsalSpinesMean     0.0221649  0.0125085   1.772 0.076579 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.528 on 1667 degrees of freedom
    ##   (24 observations deleted due to missingness)
    ## Multiple R-squared:  0.068,  Adjusted R-squared:  0.0501 
    ## F-statistic: 3.801 on 32 and 1667 DF,  p-value: 7.522e-12

``` r
rms::vif(mod)
```

    ##        BodyShapeI2s        BodyShapeI3f        BodyShapeI4e        BodyShapeI5l 
    ##           26.649284           35.090143           30.268831           12.663803 
    ##      DemersPelag2pn       DemersPelag3p      DemersPelag4po       DemersPelag5d 
    ##            1.069413            1.344441            1.162392            1.352421 
    ##      DemersPelag6bp      DemersPelag7bd OperculumPresentyes         MaxLengthTL 
    ##            1.156647            1.280800            1.363029            1.895574 
    ##               Troph            DepthMin            DepthMax         TempPrefMin 
    ##            1.522357            2.620143            1.943623            2.688081 
    ##         TempPrefMax        FeedingPathp        RepGuild22eb         RepGuild23n 
    ##            3.121962            1.369398           12.225085           31.308560 
    ##         RepGuild24t         RepGuild25h         RepGuild26s      ParentalCare2m 
    ##            7.105026            1.223312           44.986852            6.541658 
    ##      ParentalCare3p      ParentalCare4n        WaterPref2bs         WaterPref3a 
    ##           37.780554           48.319826            1.097236            1.267030 
    ##         WaterPref4b        WaterPref5fb         WaterPref6f    DorsalSpinesMean 
    ##            1.342723            1.066029            1.076586            2.088104

``` r
car::vif(mod)
```

    ##                       GVIF Df GVIF^(1/(2*Df))
    ## BodyShapeI        3.411923  4        1.165802
    ## DemersPelag       3.119262  6        1.099439
    ## OperculumPresent  1.363029  1        1.167489
    ## MaxLengthTL       1.895574  1        1.376798
    ## Troph             1.522357  1        1.233838
    ## DepthMin          2.620143  1        1.618686
    ## DepthMax          1.943623  1        1.394139
    ## TempPrefMin       2.688081  1        1.639537
    ## TempPrefMax       3.121962  1        1.766908
    ## FeedingPath       1.369398  1        1.170213
    ## RepGuild2        52.640048  5        1.486386
    ## ParentalCare     37.128063  3        1.826489
    ## WaterPref         2.034620  5        1.073614
    ## DorsalSpinesMean  2.088104  1        1.445027

``` r
# RepGuild1 and RepGuild2 are aliases

# Trait distance calculation
traits_keep <- c("BodyShapeI", "DemersPelag", "OperculumPresent", "MaxLengthTL", "Troph", "DepthMin", "DepthMax", "TempPrefMin", "TempPrefMax", "FeedingPath", "RepGuild2", "ParentalCare", "WaterPref", "DorsalSpinesMean")
# traits_vif <- traits[,traits_keep]
# traits_dist <- gower(traits_vif)
# traits_clust <- hclust(traits_dist,"average")
# traits_tree <- as.phylo(traits_clust)


## Surveyed sites
# Determine which traits have variance
nzv <- nearZeroVar(straits)
nzv
```

    ## [1] 2

``` r
# DemersPelag has no variation

straits$presence_sum <- rowSums(surveyed_sites_species)  # Sum occurrences across all sites
# Ensure straits_numeric is a data frame
mod <-  lm(presence_sum ~ BodyShapeI + OperculumPresent + MaxLengthTL + Troph + DepthMin + DepthMax + TempPrefMin + TempPrefMax + FeedingPath + RepGuild2 + ParentalCare + WaterPref + DorsalSpinesMean, straits)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = presence_sum ~ BodyShapeI + OperculumPresent + MaxLengthTL + 
    ##     Troph + DepthMin + DepthMax + TempPrefMin + TempPrefMax + 
    ##     FeedingPath + RepGuild2 + ParentalCare + WaterPref + DorsalSpinesMean, 
    ##     data = straits)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.1223 -1.7162 -0.5149  0.9708  9.2324 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)          2.747921  26.383107   0.104  0.91714   
    ## BodyShapeI2s         8.865859   4.501084   1.970  0.05011 . 
    ## BodyShapeI3f         9.678157   4.469116   2.166  0.03141 * 
    ## BodyShapeI4e         9.371918   4.451885   2.105  0.03640 * 
    ## BodyShapeI5l         8.718795   4.694057   1.857  0.06457 . 
    ## OperculumPresentyes  0.886427   0.459851   1.928  0.05517 . 
    ## MaxLengthTL         -0.005259   0.006449  -0.815  0.41569   
    ## Troph               -0.488821   0.334417  -1.462  0.14523   
    ## DepthMin            -0.121188   0.046855  -2.586  0.01033 * 
    ## DepthMax             0.013385   0.006337   2.112  0.03578 * 
    ## TempPrefMin          0.458240   0.175164   2.616  0.00950 **
    ## TempPrefMax         -0.401754   0.913469  -0.440  0.66050   
    ## FeedingPathp        -0.388418   0.479110  -0.811  0.41840   
    ## RepGuild22eb        -2.713529   2.782956  -0.975  0.33059   
    ## RepGuild23n         -4.866405   2.613964  -1.862  0.06396 . 
    ## RepGuild24t         -4.396258   2.708204  -1.623  0.10594   
    ## RepGuild26s         -4.769057   2.860485  -1.667  0.09687 . 
    ## ParentalCare2m      -3.516499   2.439742  -1.441  0.15089   
    ## ParentalCare3p      -4.130168   2.006203  -2.059  0.04069 * 
    ## ParentalCare4n      -3.924847   2.239042  -1.753  0.08099 . 
    ## WaterPref2bs         0.390157   0.479530   0.814  0.41673   
    ## WaterPref3a          2.156145   0.782921   2.754  0.00637 **
    ## DorsalSpinesMean     0.079604   0.073033   1.090  0.27690   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.623 on 223 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.2132, Adjusted R-squared:  0.1355 
    ## F-statistic: 2.746 on 22 and 223 DF,  p-value: 9.16e-05

``` r
rms::vif(mod)
```

    ##        BodyShapeI2s        BodyShapeI3f        BodyShapeI4e        BodyShapeI5l 
    ##          159.945958          172.783698          130.653845            9.488429 
    ## OperculumPresentyes         MaxLengthTL               Troph            DepthMin 
    ##            1.777374            2.525748            1.714540            1.583170 
    ##            DepthMax         TempPrefMin         TempPrefMax        FeedingPathp 
    ##            3.586064            2.123964            2.447730            1.460756 
    ##        RepGuild22eb         RepGuild23n         RepGuild24t         RepGuild26s 
    ##           15.852072           52.522651           10.224568           71.557916 
    ##      ParentalCare2m      ParentalCare3p      ParentalCare4n        WaterPref2bs 
    ##            5.063143           34.361958           43.843286            1.141666 
    ##         WaterPref3a    DorsalSpinesMean 
    ##            1.332457            2.175458

``` r
car::vif(mod)
```

    ##                       GVIF Df GVIF^(1/(2*Df))
    ## BodyShapeI        8.454690  4        1.305832
    ## OperculumPresent  1.777374  1        1.333182
    ## MaxLengthTL       2.525748  1        1.589260
    ## Troph             1.714540  1        1.309405
    ## DepthMin          1.583170  1        1.258241
    ## DepthMax          3.586064  1        1.893690
    ## TempPrefMin       2.123964  1        1.457383
    ## TempPrefMax       2.447730  1        1.564522
    ## FeedingPath       1.460756  1        1.208617
    ## RepGuild2        46.456285  4        1.615774
    ## ParentalCare     29.516419  3        1.757967
    ## WaterPref         1.465432  2        1.100250
    ## DorsalSpinesMean  2.175458  1        1.474943

``` r
# RepGuild1 and RepGuild2 are aliases

# Trait distance calculation
straits_keep <- c("BodyShapeI", "OperculumPresent", "MaxLengthTL", "Troph", "DepthMin", "DepthMax", "TempPrefMin", "TempPrefMax", "FeedingPath", "RepGuild2", "ParentalCare", "WaterPref", "DorsalSpinesMean")
straits_vif <- straits[,straits_keep]
straits_dist <- gower(straits_vif)
straits_clust <- hclust(straits_dist,"average")
straits_tree <- as.phylo(straits_clust)
```

## FD alpha Diversity

### FD alpha mntd, mpd, and pd calculations of traits and files to be saved

``` r
## Reference
# Dispersion
traits_disp <- dispersion(comm = presabs_lake, distance = traits_dist, abund = F)
write.csv(traits_disp, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/traits_disp.csv")

#Mean pairwise distances
traits_sesmpd <- ses.mpd(presabs_lake, traits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)
write.csv(traits_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/traits_sesmpd.csv")

#Mean nearest taxon distances
traits_sesmntd <- ses.mntd(presabs_lake, traits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)
write.csv(traits_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/traits_sesmntd.csv")

#Faith's PD
traits_sespd <- ses.pd(presabs_lake, traits_tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)
write.csv(traits_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/traits_sespd.csv")


## Surveyed sites
# Dispersion
straits_disp <- dispersion(comm = surveyed_sites_lake, distance = straits_dist, abund = F)
write.csv(straits_disp, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/straits_disp.csv")

#Mean pairwise distances
straits_sesmpd <- ses.mpd(surveyed_sites_lake, straits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)
write.csv(straits_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/straits_sesmpd.csv")

#Mean nearest taxon distances
straits_sesmntd <- ses.mntd(surveyed_sites_lake, straits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)
write.csv(straits_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/straits_sesmntd.csv")

#Faith's PD
straits_sespd <- ses.pd(surveyed_sites_lake, straits_tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)
write.csv(straits_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/straits_sespd.csv")
```

### FD alpha files mntd, mpd, and pd to read into R

- Read in Functional Diversity files and combine with env data

``` r
## Reference
# Dispersion 
traits_disp <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/traits_disp.csv")
traits_disp_env <- merge(traits_disp, env, by = "X", sort = F)
traits_disp_env$measure <- "traits_disp"
traits_disp_env$Site_type <- factor(traits_disp_env$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
row.names(traits_disp_env) <- traits_disp_env$X

# Mean pairwise distances
traits_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/traits_sesmpd.csv")
traits_sesmpd_env <- merge(traits_sesmpd, env, by = "X", sort = F)
traits_sesmpd_env$measure <- "traits_sesmpd"
traits_sesmpd_env$Site_type <- factor(traits_sesmpd_env$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
row.names(traits_sesmpd_env) <- traits_sesmpd_env$X

# Mean nearest taxon distances
traits_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/traits_sesmntd.csv")
traits_sesmntd_env <- merge(traits_sesmntd, env, by = "X", sort = F)
traits_sesmntd_env$measure <- "traits_sesmntd"
traits_sesmntd_env$Site_type <- factor(traits_sesmntd_env$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
row.names(traits_sesmntd_env) <- traits_sesmntd_env$X

#Faith's PD
traits_sespd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/traits_sespd.csv")
traits_sespd_env <- merge(traits_sespd, env, by = "X", sort = F)
traits_sespd_env$measure <- "traits_sespd"
traits_sespd_env$Site_type <- factor(traits_sespd_env$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
straits_sespd_env <- merge(traits_sespd_env, traits_disp, by = "X", sort = F)
row.names(traits_sespd_env) <- traits_sespd_env$X


## Surveyed sites
# Dispersion 
straits_disp <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/straits_disp.csv")
straits_disp_env <- merge(straits_disp, env, by = "X", sort = F)
straits_disp_env$measure <- "straits_disp"
straits_disp_env$Site_type <- factor(straits_disp_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))
row.names(straits_disp_env) <- straits_disp_env$X

# Mean pairwise distances
straits_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/straits_sesmpd.csv")
straits_sesmpd_env <- merge(straits_sesmpd, env, by = "X", sort = F)
straits_sesmpd_env$measure <- "straits_sesmpd"
straits_sesmpd_env$Site_type <- factor(straits_sesmpd_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))
row.names(straits_sesmpd_env) <- straits_sesmpd_env$X

# Mean nearest taxon distances
straits_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/straits_sesmntd.csv")
straits_sesmntd_env <- merge(straits_sesmntd, env, by = "X", sort = F)
straits_sesmntd_env$measure <- "straits_sesmntd"
straits_sesmntd_env$Site_type <- factor(straits_sesmntd_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))
row.names(straits_sesmntd_env) <- straits_sesmntd_env$X

#Faith's PD
straits_sespd <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FD/straits_sespd.csv")
straits_sespd_env <- merge(straits_sespd, env, by = "X", sort = F)
straits_sespd_env$measure <- "straits_sespd"
straits_sespd_env$Site_type <- factor(straits_sespd_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))
straits_sespd_env <- merge(straits_sespd_env, straits_disp, by = "X", sort = F)
row.names(straits_sespd_env) <- straits_sespd_env$X
```

### FD alpha outliers for each site type

``` r
outlier_FD_alpha <- straits_sespd_env %>%
  group_by(Site_type) %>%
  mutate(
    Q1 = quantile(pd.obs.z, 0.25),
    Q3 = quantile(pd.obs.z, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = pd.obs.z < lower_bound | pd.obs.z > upper_bound
  )

outlier_FD_alpha <- as.data.frame(outlier_FD_alpha)
row.names(outlier_FD_alpha) <- outlier_FD_alpha$X

# Create the plot
(outlier_FD_alpha_plot <- ggplot(outlier_FD_alpha, aes(x = Site_type, y = pd.obs.z, fill = Site_type)) +
  geom_violin(trim = FALSE) +
  geom_jitter(aes(color = is_outlier), width = 0.2, size = 4, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "purple")) +
  scale_fill_manual(values = custom_colors) +
  theme(text = element_text(size = 16),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(color = "black", size = 16)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "FD Alpha z-scores", x = "Site type", color = "Outlier:", tag = "c")) +
  guides(fill = "none")
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/outlier_FD_alpha.jpg", outlier_FD_alpha_plot, width = 6.26, height = 6, units = "in")
```

### Identify VIF of environmental and geographic variables

``` r
library(MASS)
```

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
# Determine which env variables have variance
environment <- c("temperature_median", "salinity_median", "oxygen_median", "pH_median")
nzv <- nearZeroVar(env[,environment])
nzv
```

    ## integer(0)

``` r
# All env variables have variance

mod <-  lm(pd.obs ~ temperature_median + salinity_median + oxygen_median + pH_median, straits_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = straits_sespd_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.83746 -0.61357 -0.05812  0.58753  1.75741 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        -43.09612   13.56172  -3.178  0.00671 **
    ## temperature_median  -0.41881    0.26022  -1.609  0.12983   
    ## salinity_median      0.08583    0.13253   0.648  0.52774   
    ## oxygen_median        0.23667    0.35888   0.659  0.52030   
    ## pH_median            7.27367    2.34591   3.101  0.00782 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.094 on 14 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.8514, Adjusted R-squared:  0.8089 
    ## F-statistic: 20.05 on 4 and 14 DF,  p-value: 1.115e-05

``` r
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                     Sum Sq Df F value   Pr(>F)   
    ## (Intercept)        12.0798  1 10.0983 0.006711 **
    ## temperature_median  3.0986  1  2.5903 0.129831   
    ## salinity_median     0.5016  1  0.4194 0.527738   
    ## oxygen_median       0.5202  1  0.4349 0.520302   
    ## pH_median          11.5000  1  9.6136 0.007823 **
    ## Residuals          16.7472 14                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
rms::vif(mod)
```

    ## temperature_median    salinity_median      oxygen_median          pH_median 
    ##           1.385456           3.435415           2.209454           4.821892

``` r
car::vif(mod)
```

    ## temperature_median    salinity_median      oxygen_median          pH_median 
    ##           1.385456           3.435415           2.209454           4.821892

``` r
mod <-  lm(pd.obs ~ temperature_median + salinity_median + oxygen_median, straits_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = straits_sespd_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9810 -0.7654 -0.4725  0.5645  3.1780 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        -10.17836   10.58758  -0.961   0.3516   
    ## temperature_median  -0.08517    0.29727  -0.287   0.7784   
    ## salinity_median      0.40409    0.10519   3.841   0.0016 **
    ## oxygen_median        0.97086    0.33836   2.869   0.0117 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.372 on 15 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.7493, Adjusted R-squared:  0.6992 
    ## F-statistic: 14.94 on 3 and 15 DF,  p-value: 8.913e-05

``` r
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                     Sum Sq Df F value   Pr(>F)   
    ## (Intercept)         1.7404  1  0.9242 0.351619   
    ## temperature_median  0.1546  1  0.0821 0.778402   
    ## salinity_median    27.7885  1 14.7564 0.001602 **
    ## oxygen_median      15.5038  1  8.2329 0.011700 * 
    ## Residuals          28.2472 15                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
rms::vif(mod)
```

    ## temperature_median    salinity_median      oxygen_median 
    ##           1.148556           1.374749           1.247588

``` r
car::vif(mod)
```

    ## temperature_median    salinity_median      oxygen_median 
    ##           1.148556           1.374749           1.247588

``` r
environment <- c("temperature_median", "salinity_median", "oxygen_median")

mod <- aov(salinity_median ~ Site_type, straits_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2 147.34   73.67   13.61 0.000353 ***
    ## Residuals   16  86.62    5.41                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 3 observations deleted due to missingness

``` r
mod <- aov(pd.obs ~ salinity_median + Site_type, straits_sespd_env)
car::vif(mod)
```

    ##                     GVIF Df GVIF^(1/(2*Df))
    ## salinity_median 2.701125  1        1.643510
    ## Site_type       2.701125  2        1.281995

``` r
mod <- aov(oxygen_median ~ Site_type, straits_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2 14.440    7.22      19 5.95e-05 ***
    ## Residuals   16  6.081    0.38                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 3 observations deleted due to missingness

``` r
mod <- aov(pd.obs ~ oxygen_median + Site_type, straits_sespd_env)
car::vif(mod)
```

    ##                   GVIF Df GVIF^(1/(2*Df))
    ## oxygen_median 3.374424  1        1.836961
    ## Site_type     3.374424  2        1.355345

``` r
mod <- aov(temperature_median ~ Site_type, straits_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2   2.94   1.470   1.092  0.359
    ## Residuals   16  21.54   1.346               
    ## 3 observations deleted due to missingness

``` r
mod <- aov(pd.obs ~ temperature_median + Site_type, straits_sespd_env)
car::vif(mod)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## temperature_median 1.136543  1        1.066088
    ## Site_type          1.136543  2        1.032515

``` r
SR_ss_env <- straits_sespd_env[surveyed_sites_env,]
SR_ss_env$Site_type <- factor(SR_ss_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

# Run linear discriminant analysis of environmental variables
LDA <- lda(SR_ss_env[,environment], SR_ss_env$Site_type)
spe.class <- predict(LDA)$class
spe.post <- predict(LDA)$posterior
(spe.table <- table(SR_ss_env$Site_type, spe.class))
```

    ##             spe.class
    ##              Ocean Mixed Stratified
    ##   Ocean          1     2          0
    ##   Mixed          0     8          0
    ##   Stratified     0     0          8

``` r
diag(prop.table(spe.table, 1))
```

    ##      Ocean      Mixed Stratified 
    ##  0.3333333  1.0000000  1.0000000

``` r
# Determine which geo variables have variance
geography <- c("volume_m3_w_chemocline", "volume_m3", "surface_area_m2", "distance_to_ocean_min_m", "distance_to_ocean_mean_m", "distance_to_ocean_median_m", "tidal_lag_time_minutes", "tidal_efficiency", "perimeter_fromSat", "max_depth", "logArea")
nzv <- nearZeroVar(env[,geography])
nzv
```

    ## integer(0)

``` r
# All geo variables have variance

mod <-  lm(pd.obs ~ volume_m3_w_chemocline + volume_m3 + surface_area_m2 + distance_to_ocean_min_m + distance_to_ocean_mean_m + distance_to_ocean_median_m + tidal_lag_time_minutes + tidal_efficiency + perimeter_fromSat + max_depth + logArea, straits_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ volume_m3_w_chemocline + volume_m3 + surface_area_m2 + 
    ##     distance_to_ocean_min_m + distance_to_ocean_mean_m + distance_to_ocean_median_m + 
    ##     tidal_lag_time_minutes + tidal_efficiency + perimeter_fromSat + 
    ##     max_depth + logArea, data = straits_sespd_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.31770 -0.79585  0.08146  0.80992  2.25461 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                -7.002e+00  1.093e+01  -0.641    0.538
    ## volume_m3_w_chemocline      3.602e-06  3.994e-06   0.902    0.391
    ## volume_m3                  -3.063e-06  4.137e-06  -0.740    0.478
    ## surface_area_m2            -7.011e-06  3.981e-06  -1.761    0.112
    ## distance_to_ocean_min_m    -7.433e-03  2.183e-02  -0.340    0.741
    ## distance_to_ocean_mean_m    3.557e-02  9.273e-02   0.384    0.710
    ## distance_to_ocean_median_m -3.767e-02  8.523e-02  -0.442    0.669
    ## tidal_lag_time_minutes     -1.020e-02  2.559e-02  -0.399    0.699
    ## tidal_efficiency           -1.327e-01  6.817e+00  -0.019    0.985
    ## perimeter_fromSat          -1.336e-03  1.685e-03  -0.793    0.448
    ## max_depth                  -5.807e-02  1.123e-01  -0.517    0.618
    ## logArea                     1.527e+00  1.091e+00   1.400    0.195
    ## 
    ## Residual standard error: 1.815 on 9 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.6945, Adjusted R-squared:  0.3212 
    ## F-statistic:  1.86 on 11 and 9 DF,  p-value: 0.1802

``` r
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                             Sum Sq Df F value Pr(>F)
    ## (Intercept)                 1.3526  1  0.4106 0.5377
    ## volume_m3_w_chemocline      2.6794  1  0.8133 0.3907
    ## volume_m3                   1.8054  1  0.5480 0.4780
    ## surface_area_m2            10.2193  1  3.1018 0.1120
    ## distance_to_ocean_min_m     0.3819  1  0.1159 0.7413
    ## distance_to_ocean_mean_m    0.4847  1  0.1471 0.7102
    ## distance_to_ocean_median_m  0.6436  1  0.1954 0.6689
    ## tidal_lag_time_minutes      0.5234  1  0.1589 0.6995
    ## tidal_efficiency            0.0012  1  0.0004 0.9849
    ## perimeter_fromSat           2.0710  1  0.6286 0.4483
    ## max_depth                   0.8810  1  0.2674 0.6176
    ## logArea                     6.4561  1  1.9596 0.1951
    ## Residuals                  29.6515  9

``` r
rms::vif(mod)
```

    ##     volume_m3_w_chemocline                  volume_m3 
    ##                 6433.74411                 6855.30224 
    ##            surface_area_m2    distance_to_ocean_min_m 
    ##                   22.87840                   18.42491 
    ##   distance_to_ocean_mean_m distance_to_ocean_median_m 
    ##                  825.68615                  687.45478 
    ##     tidal_lag_time_minutes           tidal_efficiency 
    ##                   21.17746                   27.90908 
    ##          perimeter_fromSat                  max_depth 
    ##                   94.39757                   10.74366 
    ##                    logArea 
    ##                   25.70304

``` r
car::vif(mod)
```

    ##     volume_m3_w_chemocline                  volume_m3 
    ##                 6433.74411                 6855.30224 
    ##            surface_area_m2    distance_to_ocean_min_m 
    ##                   22.87840                   18.42491 
    ##   distance_to_ocean_mean_m distance_to_ocean_median_m 
    ##                  825.68615                  687.45478 
    ##     tidal_lag_time_minutes           tidal_efficiency 
    ##                   21.17746                   27.90908 
    ##          perimeter_fromSat                  max_depth 
    ##                   94.39757                   10.74366 
    ##                    logArea 
    ##                   25.70304

``` r
mod <-  lm(pd.obs ~ distance_to_ocean_min_m + max_depth + logArea, straits_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ distance_to_ocean_min_m + max_depth + logArea, 
    ##     data = straits_sespd_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3790 -0.8852 -0.0909  0.8312  5.0222 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)              3.613318   2.826549   1.278   0.2174   
    ## distance_to_ocean_min_m -0.018558   0.006373  -2.912   0.0093 **
    ## max_depth                0.000272   0.046853   0.006   0.9954   
    ## logArea                  0.154562   0.280858   0.550   0.5889   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.029 on 18 degrees of freedom
    ## Multiple R-squared:  0.3957, Adjusted R-squared:  0.295 
    ## F-statistic: 3.929 on 3 and 18 DF,  p-value: 0.02553

``` r
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                         Sum Sq Df F value   Pr(>F)   
    ## (Intercept)              6.730  1  1.6342 0.217365   
    ## distance_to_ocean_min_m 34.918  1  8.4791 0.009303 **
    ## max_depth                0.000  1  0.0000 0.995432   
    ## logArea                  1.247  1  0.3029 0.588865   
    ## Residuals               74.127 18                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
rms::vif(mod)
```

    ## distance_to_ocean_min_m               max_depth                 logArea 
    ##                1.258244                1.497112                1.392300

``` r
car::vif(mod)
```

    ## distance_to_ocean_min_m               max_depth                 logArea 
    ##                1.258244                1.497112                1.392300

``` r
geography <- c("distance_to_ocean_min_m", "max_depth", "logArea")

mod <- aov(distance_to_ocean_min_m ~ Site_type, straits_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2  80935   40468   16.49 7.05e-05 ***
    ## Residuals   19  46637    2455                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mod <- aov(pd.obs ~ distance_to_ocean_min_m + Site_type, straits_sespd_env)
car::vif(mod)
```

    ##                             GVIF Df GVIF^(1/(2*Df))
    ## distance_to_ocean_min_m 2.735421  1        1.653911
    ## Site_type               2.735421  2        1.286045

``` r
mod <- aov(max_depth ~ Site_type, straits_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2    535   267.5   2.235  0.134
    ## Residuals   19   2274   119.7

``` r
mod <- aov(pd.obs ~ max_depth + Site_type, straits_sespd_env)
car::vif(mod)
```

    ##               GVIF Df GVIF^(1/(2*Df))
    ## max_depth 1.235315  1        1.111447
    ## Site_type 1.235315  2        1.054252

``` r
mod <- aov(logArea ~ Site_type, straits_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2  14.37   7.184   2.341  0.123
    ## Residuals   19  58.32   3.069

``` r
mod <- aov(pd.obs ~ logArea + Site_type, straits_sespd_env)
car::vif(mod)
```

    ##               GVIF Df GVIF^(1/(2*Df))
    ## logArea   1.246371  1        1.116410
    ## Site_type 1.246371  2        1.056603

``` r
SR_ss_geo <- straits_sespd_env[surveyed_sites,]
SR_ss_geo$Site_type <- factor(SR_ss_geo$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

# Run linear discriminant analysis of geographical variables
LDA <- lda(SR_ss_geo[,geography], SR_ss_geo$Site_type)
spe.class <- predict(LDA)$class
spe.post <- predict(LDA)$posterior
(spe.table <- table(SR_ss_geo$Site_type, spe.class))
```

    ##             spe.class
    ##              Ocean Mixed Stratified
    ##   Ocean          4     2          0
    ##   Mixed          0     8          0
    ##   Stratified     0     2          6

``` r
diag(prop.table(spe.table, 1))
```

    ##      Ocean      Mixed Stratified 
    ##  0.6666667  1.0000000  0.7500000

``` r
# Combine environment and geography
mod <-  lm(pd.obs ~ salinity_median + oxygen_median + temperature_median + distance_to_ocean_min_m + max_depth + logArea, straits_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ salinity_median + oxygen_median + temperature_median + 
    ##     distance_to_ocean_min_m + max_depth + logArea, data = straits_sespd_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.41579 -0.57810 -0.05263  0.66856  1.91595 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -1.192e+01  1.216e+01  -0.980   0.3463  
    ## salinity_median          4.180e-01  1.596e-01   2.619   0.0224 *
    ## oxygen_median            6.673e-01  3.416e-01   1.954   0.0745 .
    ## temperature_median      -1.993e-01  3.277e-01  -0.608   0.5543  
    ## distance_to_ocean_min_m  6.698e-04  7.531e-03   0.089   0.9306  
    ## max_depth               -6.971e-02  4.434e-02  -1.572   0.1419  
    ## logArea                  7.276e-01  3.036e-01   2.396   0.0338 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.259 on 12 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.8311, Adjusted R-squared:  0.7466 
    ## F-statistic: 9.841 on 6 and 12 DF,  p-value: 0.0004759

``` r
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                          Sum Sq Df F value  Pr(>F)  
    ## (Intercept)              1.5243  1  0.9611 0.34627  
    ## salinity_median         10.8775  1  6.8584 0.02243 *
    ## oxygen_median            6.0534  1  3.8167 0.07445 .
    ## temperature_median       0.5870  1  0.3701 0.55428  
    ## distance_to_ocean_min_m  0.0125  1  0.0079 0.93059  
    ## max_depth                3.9193  1  2.4711 0.14193  
    ## logArea                  9.1064  1  5.7417 0.03375 *
    ## Residuals               19.0321 12                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
rms::vif(mod)
```

    ##         salinity_median           oxygen_median      temperature_median 
    ##                3.758774                1.509635                1.656753 
    ## distance_to_ocean_min_m               max_depth                 logArea 
    ##                3.720344                2.979442                2.810937

``` r
car::vif(mod)
```

    ##         salinity_median           oxygen_median      temperature_median 
    ##                3.758774                1.509635                1.656753 
    ## distance_to_ocean_min_m               max_depth                 logArea 
    ##                3.720344                2.979442                2.810937

``` r
mod <-  lm(pd.obs ~ oxygen_median + temperature_median + distance_to_ocean_min_m + logArea, straits_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ oxygen_median + temperature_median + distance_to_ocean_min_m + 
    ##     logArea, data = straits_sespd_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1507 -0.9751 -0.2223  1.0147  2.9702 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             12.84425   10.78255   1.191   0.2534  
    ## oxygen_median            0.85968    0.40231   2.137   0.0507 .
    ## temperature_median      -0.51124    0.37257  -1.372   0.1916  
    ## distance_to_ocean_min_m -0.01499    0.00556  -2.696   0.0174 *
    ## logArea                  0.43729    0.25945   1.685   0.1141  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.536 on 14 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.707,  Adjusted R-squared:  0.6233 
    ## F-statistic: 8.447 on 4 and 14 DF,  p-value: 0.001102

``` r
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                         Sum Sq Df F value  Pr(>F)  
    ## (Intercept)              3.346  1  1.4190 0.25338  
    ## oxygen_median           10.767  1  4.5663 0.05074 .
    ## temperature_median       4.440  1  1.8829 0.19159  
    ## distance_to_ocean_min_m 17.138  1  7.2685 0.01739 *
    ## logArea                  6.698  1  2.8406 0.11406  
    ## Residuals               33.009 14                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
rms::vif(mod)
```

    ##           oxygen_median      temperature_median distance_to_ocean_min_m 
    ##                1.408629                1.440908                1.364232 
    ##                 logArea 
    ##                1.380536

``` r
car::vif(mod)
```

    ##           oxygen_median      temperature_median distance_to_ocean_min_m 
    ##                1.408629                1.440908                1.364232 
    ##                 logArea 
    ##                1.380536

``` r
envgeo <- c("oxygen_median", "temperature_median", "distance_to_ocean_min_m", "logArea")


# Run linear discriminant analysis of geographical variables
LDA <- lda(SR_ss_env[,envgeo], SR_ss_env$Site_type)
spe.class <- predict(LDA)$class
spe.post <- predict(LDA)$posterior
(spe.table <- table(SR_ss_env$Site_type, spe.class))
```

    ##             spe.class
    ##              Ocean Mixed Stratified
    ##   Ocean          3     0          0
    ##   Mixed          1     7          0
    ##   Stratified     0     1          7

``` r
diag(prop.table(spe.table, 1))
```

    ##      Ocean      Mixed Stratified 
    ##      1.000      0.875      0.875

### FD alpha & site type

``` r
### Site_type
# Run ANOVA on the original data
anova_FRic_result <- aov(pd.obs ~ Site_type, data = straits_sespd_env[surveyed_sites,])
summary(anova_FRic_result)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2  74.23   37.12   14.56 0.000147 ***
    ## Residuals   19  48.44    2.55                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Perform Tukey's HSD test
tukey_FRic_result <- TukeyHSD(anova_FRic_result)
print(tukey_FRic_result)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = pd.obs ~ Site_type, data = straits_sespd_env[surveyed_sites, ])
    ## 
    ## $Site_type
    ##                        diff       lwr       upr     p adj
    ## Mixed-Ocean       0.4304859 -1.760127  2.621099 0.8725386
    ## Stratified-Ocean -3.5561361 -5.746749 -1.365523 0.0015963
    ## Stratified-Mixed -3.9866221 -6.014735 -1.958509 0.0002278

``` r
# Sites
N <- length(anova_FRic_result$residuals)
# Categories
k <- length(anova_FRic_result$coefficients)
# Sum of squares
SS <- sum(resid(anova_FRic_result)^2) # from anova, the residual SS
# Mean squares
MS <- SS/(N-k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

# Get differences
Site_type <- tukey_FRic_result$Site_type

# q-values
(OM_q <- (abs(Site_type[1,1]))/USE)
```

    ## [1] 0.7060233

``` r
(MS_q <- (abs(Site_type[3,1]))/BSE)
```

    ## [1] 7.062175

``` r
(SO_q <- (abs(Site_type[2,1]))/USE)
```

    ## [1] 5.832281

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx

# Run ANOVA on the original data
anova_FRicz_result <- aov(pd.obs.z ~ Site_type, data = straits_sespd_env[surveyed_sites,])
summary(anova_FRicz_result)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type    2  13.83   6.915   5.008 0.0179 *
    ## Residuals   19  26.24   1.381                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Perform Tukey's HSD test
tukey_FRicz_result <- TukeyHSD(anova_FRicz_result)
print(tukey_FRicz_result)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = pd.obs.z ~ Site_type, data = straits_sespd_env[surveyed_sites, ])
    ## 
    ## $Site_type
    ##                       diff        lwr      upr     p adj
    ## Mixed-Ocean      1.5069087 -0.1053032 3.119121 0.0694252
    ## Stratified-Ocean 1.9506386  0.3384267 3.562851 0.0164204
    ## Stratified-Mixed 0.4437299 -1.0488884 1.936348 0.7341934

``` r
# Sites
N <- length(anova_FRicz_result$residuals)
# Categories
k <- length(anova_FRicz_result$coefficients)
# Sum of squares
SS <- sum(resid(anova_FRicz_result)^2) # from anova, the residual SS
# Mean squares
MS <- SS/(N-k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

# Get differences
Site_type <- tukey_FRicz_result$Site_type

# q-values
(OM_q <- (abs(Site_type[1,1]))/USE)
```

    ## [1] 3.358076

``` r
(MS_q <- (abs(Site_type[3,1]))/BSE)
```

    ## [1] 1.06806

``` r
(SO_q <- (abs(Site_type[2,1]))/USE)
```

    ## [1] 4.346907

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx

# Run ANOVA on the original data
anova_FDisp_result <- aov(Dispersion ~ Site_type, data = straits_sespd_env[surveyed_sites,])
summary(anova_FDisp_result)
```

    ##             Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Site_type    2 0.00826 0.004132   2.121  0.147
    ## Residuals   19 0.03701 0.001948

``` r
# Perform Tukey's HSD test
tukey_FDisp_result <- TukeyHSD(anova_FDisp_result)
print(tukey_FDisp_result)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = Dispersion ~ Site_type, data = straits_sespd_env[surveyed_sites, ])
    ## 
    ## $Site_type
    ##                         diff         lwr        upr     p adj
    ## Mixed-Ocean      0.041543913 -0.01901261 0.10210044 0.2156444
    ## Stratified-Ocean 0.045209326 -0.01534720 0.10576585 0.1669220
    ## Stratified-Mixed 0.003665412 -0.05239903 0.05972986 0.9849189

``` r
# Sites
N <- length(anova_FDisp_result$residuals)
# Categories
k <- length(anova_FDisp_result$coefficients)
# Sum of squares
SS <- sum(resid(anova_FDisp_result)^2) # from anova, the residual SS
# Mean squares
MS <- SS/(N-k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

# Get differences
Site_type <- tukey_FDisp_result$Site_type

# q-values
(OM_q <- (abs(Site_type[1,1]))/USE)
```

    ## [1] 2.464746

``` r
(MS_q <- (abs(Site_type[3,1]))/BSE)
```

    ## [1] 0.2348881

``` r
(SO_q <- (abs(Site_type[2,1]))/USE)
```

    ## [1] 2.68221

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx
```

### FD alpha with env and geo linear models & ANOVAs

``` r
par(mfrow=c(2,2)) 
##### FRicz ANOVAs
#### Environmental
### Surveyed sites
## Temperature
# Interaction
FD_FRicz_am_temp <- aov(pd.obs.z ~ Site_type * temperature_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_temp)
    ## W = 0.88867, p-value = 0.03049

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_temp) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.3137 0.7351
    ##       16

``` r
# Check for outliers
outlierTest(FD_FRicz_am_temp)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -3.991594          0.0017886     0.033983

``` r
# Plot residuals
plot(FD_FRicz_am_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-1.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FRicz_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                     2  9.003   4.502   3.088 0.0799 .
    ## temperature_median            1  0.017   0.017   0.012 0.9155  
    ## Site_type:temperature_median  2  3.730   1.865   1.280 0.3110  
    ## Residuals                    13 18.949   1.458                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(FD_FRicz_am_temp)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.2397285 1.0000000 0.9329753

``` r
# ANCOVA
FD_FRicz_am_temp <- aov(pd.obs.z ~ Site_type + temperature_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_temp)
    ## W = 0.93906, p-value = 0.2536

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_temp) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.6434 0.5386
    ##       16

``` r
# Check for outliers
outlierTest(FD_FRicz_am_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLO -3.27426          0.0055388      0.10524

``` r
# Plot residuals
plot(FD_FRicz_am_temp)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-2.png)<!-- -->

``` r
## Salinity
# Interaction
FD_FRicz_am_sal <- aov(pd.obs.z ~ Site_type * salinity_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_sal)
    ## W = 0.90983, p-value = 0.07351

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_sal) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.6958 0.5131
    ##       16

``` r
# Check for outliers
outlierTest(FD_FRicz_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -3.174753          0.0079987      0.15198

``` r
# Plot residuals
plot(FD_FRicz_am_sal)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-3.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FRicz_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                  2  9.003   4.502   3.081 0.0803 .
    ## salinity_median            1  0.866   0.866   0.593 0.4551  
    ## Site_type:salinity_median  2  2.838   1.419   0.971 0.4044  
    ## Residuals                 13 18.991   1.461                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
sal_p_values <- summary(FD_FRicz_am_sal)[[1]][, "Pr(>F)"]
p_values <- c(sal_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.2408546 1.0000000 1.0000000

``` r
# ANCOVA
FD_FRicz_am_sal <- aov(pd.obs.z ~ Site_type + salinity_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_sal)
    ## W = 0.93044, p-value = 0.1764

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_sal) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  1.1097 0.3537
    ##       16

``` r
# Check for outliers
outlierTest(FD_FRicz_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -3.423837          0.0041132      0.07815

``` r
# Plot residuals
plot(FD_FRicz_am_sal)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-4.png)<!-- -->

``` r
## Oxygen
# Interaction
FD_FRicz_am_oxy <- aov(pd.obs.z ~ Site_type * oxygen_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_oxy)
    ## W = 0.95447, p-value = 0.4691

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_oxy) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  1.1446 0.3431
    ##       16

``` r
# Check for outliers
outlierTest(FD_FRicz_am_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -3.677548           0.003163     0.060097

``` r
# Plot residuals
plot(FD_FRicz_am_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-5.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FRicz_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value  Pr(>F)   
    ## Site_type                2  9.003   4.502   6.378 0.01174 * 
    ## oxygen_median            1  8.732   8.732  12.373 0.00378 **
    ## Site_type:oxygen_median  2  4.789   2.394   3.392 0.06523 . 
    ## Residuals               13  9.175   0.706                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
oxy_p_values <- summary(FD_FRicz_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.03523476 0.01135418 0.19569759

``` r
# ANCOVA
FD_FRicz_am_oxy <- aov(pd.obs.z ~ Site_type + oxygen_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_oxy)
    ## W = 0.94666, p-value = 0.3461

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_oxy) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.5998 0.5608
    ##       16

``` r
# Check for outliers
outlierTest(FD_FRicz_am_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLO -3.34533          0.0048084      0.09136

``` r
# Plot residuals
plot(FD_FRicz_am_oxy)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-6.png)<!-- -->

``` r
# Summarize ANOVA results and calculate pairwise comparisons
summary(FD_FRicz_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type           2  9.003   4.502   2.977 0.0815 .
    ## temperature_median  1  0.017   0.017   0.011 0.9168  
    ## Residuals          15 22.679   1.512                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_FRicz_amp_temp <- emmeans(FD_FRicz_am_temp, pairwise ~ Site_type, adjust = "bonferroni")
FD_FRicz_amp_temp$contrasts
```

    ##  contrast           estimate    SE df t.ratio p.value
    ##  Ocean - Mixed        -1.572 0.837 15  -1.877  0.2402
    ##  Ocean - Stratified   -2.040 0.844 15  -2.418  0.0864
    ##  Mixed - Stratified   -0.468 0.655 15  -0.714  1.0000
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(FD_FRicz_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type        2  9.003   4.502   3.093  0.075 .
    ## salinity_median  1  0.866   0.866   0.595  0.452  
    ## Residuals       15 21.830   1.455                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_FRicz_amp_sal <- emmeans(FD_FRicz_am_sal, pairwise ~ Site_type, adjust = "bonferroni")
FD_FRicz_amp_sal$contrasts
```

    ##  contrast           estimate    SE df t.ratio p.value
    ##  Ocean - Mixed        -1.638 0.820 15  -1.997  0.1928
    ##  Ocean - Stratified   -2.629 1.130 15  -2.324  0.1037
    ##  Mixed - Stratified   -0.991 0.931 15  -1.064  0.9122
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(FD_FRicz_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value  Pr(>F)   
    ## Site_type      2  9.003   4.502   4.836 0.02395 * 
    ## oxygen_median  1  8.732   8.732   9.381 0.00789 **
    ## Residuals     15 13.964   0.931                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_FRicz_amp_oxy <- emmeans(FD_FRicz_am_oxy, pairwise ~ Site_type, adjust = "bonferroni")
FD_FRicz_amp_oxy$contrasts
```

    ##  contrast           estimate    SE df t.ratio p.value
    ##  Ocean - Mixed        -0.632 0.723 15  -0.874  1.0000
    ##  Ocean - Stratified    0.678 1.100 15   0.618  1.0000
    ##  Mixed - Stratified    1.310 0.749 15   1.750  0.3017
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
# p-values
temp_p_values <- summary(FD_FRicz_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(FD_FRicz_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(FD_FRicz_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.48888764 1.00000000         NA 0.45019983 1.00000000         NA 0.14367016
    ## [8] 0.04736936         NA

``` r
### Mixed lakes
FD_FRicz_env_M_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, straits_sespd_env[mixed_lakes,])
summary(FD_FRicz_env_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = straits_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##      FLK      HLO      LLN      MLN      NLN      NLU      OLO      ULN 
    ## -0.61048 -1.77239 -0.22376  1.46783  0.22579  0.03092  1.03129 -0.14921 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         29.7425    61.7605   0.482    0.655
    ## salinity_median     -0.5172     1.4358  -0.360    0.737
    ## oxygen_median       -1.8203     1.1614  -1.567    0.192
    ## temperature_median  -0.1500     0.8912  -0.168    0.874
    ## 
    ## Residual standard error: 1.309 on 4 degrees of freedom
    ## Multiple R-squared:  0.4868, Adjusted R-squared:  0.1019 
    ## F-statistic: 1.265 on 3 and 4 DF,  p-value: 0.3989

``` r
Anova(FD_FRicz_env_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                    Sum Sq Df F value Pr(>F)
    ## (Intercept)        0.3975  1  0.2319 0.6553
    ## salinity_median    0.2225  1  0.1298 0.7369
    ## oxygen_median      4.2107  1  2.4565 0.1921
    ## temperature_median 0.0486  1  0.0283 0.8745
    ## Residuals          6.8564  4

``` r
p_values <- summary(FD_FRicz_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          1.0000000          0.7684344          1.0000000

``` r
### Stratified lakes
FD_FRicz_env_S_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, straits_sespd_env[stratified_lakes,])
summary(FD_FRicz_env_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = straits_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##      BCM      CLM      GLK      HLM      NLK      OTM      SLN      TLN 
    ## -0.09011 -0.19557  0.68363 -0.12195 -0.57115  0.60031 -0.33905  0.03389 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.930946   6.378917  -0.146    0.891
    ## salinity_median     0.006006   0.099066   0.061    0.955
    ## oxygen_median      -0.708241   0.436772  -1.622    0.180
    ## temperature_median  0.100618   0.141284   0.712    0.516
    ## 
    ## Residual standard error: 0.5769 on 4 degrees of freedom
    ## Multiple R-squared:  0.6384, Adjusted R-squared:  0.3673 
    ## F-statistic: 2.354 on 3 and 4 DF,  p-value: 0.2132

``` r
Anova(FD_FRicz_env_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                     Sum Sq Df F value Pr(>F)
    ## (Intercept)        0.00709  1  0.0213 0.8910
    ## salinity_median    0.00122  1  0.0037 0.9546
    ## oxygen_median      0.87511  1  2.6294 0.1802
    ## temperature_median 0.16880  1  0.5072 0.5157
    ## Residuals          1.33128  4

``` r
p_values <- summary(FD_FRicz_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          1.0000000          0.7208861          1.0000000

``` r
#### Geographical
### Surveyed sites
## Distance
# Interaction
FD_FRicz_am_dist <- aov(pd.obs.z ~ Site_type * distance_to_ocean_min_m, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_dist)
    ## W = 0.95485, p-value = 0.3927

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_dist) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  2.8964 0.07981 .
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(FD_FRicz_am_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -2.614257           0.019536      0.41026

``` r
# Plot residuals
plot(FD_FRicz_am_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   19

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-7.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FRicz_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                          2 13.829   6.915   5.495 0.0152 *
    ## distance_to_ocean_min_m            1  0.502   0.502   0.399 0.5365  
    ## Site_type:distance_to_ocean_min_m  2  5.601   2.801   2.226 0.1403  
    ## Residuals                         16 20.132   1.258                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
dist_p_values <- summary(FD_FRicz_am_dist)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.04574943 1.00000000 0.42099277

``` r
# ANCOVA
FD_FRicz_am_dist <- aov(pd.obs.z ~ Site_type + distance_to_ocean_min_m, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_dist)
    ## W = 0.92361, p-value = 0.09034

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_dist) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  1.5084 0.2466
    ##       19

``` r
# Check for outliers
outlierTest(FD_FRicz_am_dist)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -3.629915            0.00207     0.045541

``` r
# Plot residuals
plot(FD_FRicz_am_dist)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-8.png)<!-- -->

``` r
## Max Depth
# Interaction
FD_FRicz_am_mxd <- aov(pd.obs.z ~ Site_type * max_depth, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_mxd)
    ## W = 0.97462, p-value = 0.8148

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_mxd) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.2106  0.812
    ##       19

``` r
# Check for outliers
outlierTest(FD_FRicz_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -2.591101           0.020458      0.45008

``` r
# Plot residuals
plot(FD_FRicz_am_mxd)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-9.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FRicz_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value  Pr(>F)   
    ## Site_type            2 13.829   6.915   6.775 0.00739 **
    ## max_depth            1  5.070   5.070   4.968 0.04051 * 
    ## Site_type:max_depth  2  4.836   2.418   2.369 0.12552   
    ## Residuals           16 16.329   1.021                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
mxd_p_values <- summary(FD_FRicz_am_mxd)[[1]][, "Pr(>F)"]
p_values <- c(mxd_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.02216163 0.12152254 0.37657012

``` r
# ANCOVA
FD_FRicz_am_mxd <- aov(pd.obs.z ~ Site_type + max_depth, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_mxd)
    ## W = 0.91886, p-value = 0.072

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_mxd) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.0353 0.9654
    ##       19

``` r
# Check for outliers
outlierTest(FD_FRicz_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -2.919671          0.0095555      0.21022

``` r
# Plot residuals
plot(FD_FRicz_am_mxd)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-10.png)<!-- -->

``` r
## Log Area
# Interaction
FD_FRicz_am_lga <- aov(pd.obs.z ~ Site_type * logArea, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_lga)
    ## W = 0.91956, p-value = 0.07443

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_lga) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.0837   0.92
    ##       19

``` r
# Check for outliers
outlierTest(FD_FRicz_am_lga)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -3.909146           0.001395     0.030689

``` r
# Plot residuals
plot(FD_FRicz_am_lga)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-11.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FRicz_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type          2 13.829   6.915   6.156 0.0104 *
    ## logArea            1  6.496   6.496   5.784 0.0286 *
    ## Site_type:logArea  2  1.769   0.885   0.788 0.4717  
    ## Residuals         16 17.970   1.123                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
lga_p_values <- summary(FD_FRicz_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.03120403 0.08590842 1.00000000

``` r
# ANCOVA
FD_FRicz_am_lga <- aov(pd.obs.z ~ Site_type + logArea, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FRicz_am_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FRicz_am_lga)
    ## W = 0.91604, p-value = 0.06296

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FRicz_am_lga) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.1025 0.9031
    ##       19

``` r
# Check for outliers
outlierTest(FD_FRicz_am_lga)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -3.859162          0.0012584     0.027685

``` r
# Plot residuals
plot(FD_FRicz_am_lga)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-12.png)<!-- -->

``` r
# Summarize ANOVA results and calculate pairwise comparisons
summary(FD_FRicz_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                2 13.829   6.915   4.837 0.0208 *
    ## distance_to_ocean_min_m  1  0.502   0.502   0.351 0.5608  
    ## Residuals               18 25.733   1.430                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_FRicz_amp_dist <- emmeans(FD_FRicz_am_dist, pairwise ~ Site_type, adjust = "bonferroni")
FD_FRicz_amp_dist$contrasts
```

    ##  contrast           estimate    SE df t.ratio p.value
    ##  Ocean - Mixed        -1.718 0.738 18  -2.329  0.0951
    ##  Ocean - Stratified   -2.447 1.060 18  -2.314  0.0981
    ##  Mixed - Stratified   -0.729 0.767 18  -0.950  1.0000
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(FD_FRicz_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type    2  13.83   6.915   5.880 0.0108 *
    ## max_depth    1   5.07   5.070   4.312 0.0524 .
    ## Residuals   18  21.16   1.176                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_FRicz_amp_mxd <- emmeans(FD_FRicz_am_mxd, pairwise ~ Site_type, adjust = "bonferroni")
FD_FRicz_amp_mxd$contrasts
```

    ##  contrast           estimate    SE df t.ratio p.value
    ##  Ocean - Mixed        -1.446 0.586 18  -2.466  0.0718
    ##  Ocean - Stratified   -2.397 0.624 18  -3.843  0.0036
    ##  Mixed - Stratified   -0.951 0.595 18  -1.600  0.3813
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(FD_FRicz_am_lga)
```

    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## Site_type    2 13.829   6.915   6.305 0.00841 **
    ## logArea      1  6.496   6.496   5.923 0.02558 * 
    ## Residuals   18 19.740   1.097                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_FRicz_amp_lga <- emmeans(FD_FRicz_am_lga, pairwise ~ Site_type, adjust = "bonferroni")
FD_FRicz_amp_lga$contrasts
```

    ##  contrast           estimate    SE df t.ratio p.value
    ##  Ocean - Mixed        -0.832 0.630 18  -1.321  0.6090
    ##  Ocean - Stratified   -1.477 0.598 18  -2.470  0.0712
    ##  Mixed - Stratified   -0.645 0.530 18  -1.217  0.7179
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
# p-values
dist_p_values <- summary(FD_FRicz_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(FD_FRicz_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(FD_FRicz_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.12504911 1.00000000         NA 0.06498293 0.31466408         NA 0.05044252
    ## [8] 0.15350312         NA

``` r
### Mixed lakes
FD_FRicz_geo_M_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, straits_sespd_env[mixed_lakes,])
summary(FD_FRicz_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = straits_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##     FLK     HLO     LLN     MLN     NLN     NLU     OLO     ULN 
    ## -0.6273 -1.7255 -0.4627  1.2569  1.1317  0.3203  0.2481 -0.1415 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              1.70380    3.58957   0.475    0.660
    ## distance_to_ocean_min_m  0.01580    0.02123   0.744    0.498
    ## max_depth               -0.04308    0.09275  -0.464    0.666
    ## logArea                 -0.26433    0.47610  -0.555    0.608
    ## 
    ## Residual standard error: 1.287 on 4 degrees of freedom
    ## Multiple R-squared:  0.5038, Adjusted R-squared:  0.1317 
    ## F-statistic: 1.354 on 3 and 4 DF,  p-value: 0.3762

``` r
Anova(FD_FRicz_geo_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                         Sum Sq Df F value Pr(>F)
    ## (Intercept)             0.3734  1  0.2253 0.6598
    ## distance_to_ocean_min_m 0.9186  1  0.5542 0.4979
    ## max_depth               0.3575  1  0.2157 0.6665
    ## logArea                 0.5109  1  0.3082 0.6083
    ## Residuals               6.6294  4

``` r
p_values <- summary(FD_FRicz_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                       1                       1                       1 
    ##                 logArea 
    ##                       1

``` r
### Stratified lakes
FD_FRicz_geo_S_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, straits_sespd_env[stratified_lakes,])
summary(FD_FRicz_geo_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = straits_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##      BCM      CLM      GLK      HLM      NLK      OTM      SLN      TLN 
    ##  0.11092  0.13801 -0.41897 -0.46751 -0.01302  0.40110 -0.19225  0.44172 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -1.739347   2.785779  -0.624   0.5662  
    ## distance_to_ocean_min_m -0.008809   0.002367  -3.722   0.0204 *
    ## max_depth               -0.018654   0.033622  -0.555   0.6086  
    ## logArea                  0.357353   0.344625   1.037   0.3583  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4524 on 4 degrees of freedom
    ## Multiple R-squared:  0.7777, Adjusted R-squared:  0.611 
    ## F-statistic: 4.664 on 3 and 4 DF,  p-value: 0.08549

``` r
Anova(FD_FRicz_geo_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                          Sum Sq Df F value  Pr(>F)  
    ## (Intercept)             0.07978  1  0.3898 0.56623  
    ## distance_to_ocean_min_m 2.83428  1 13.8498 0.02045 *
    ## max_depth               0.06300  1  0.3078 0.60858  
    ## logArea                 0.22004  1  1.0752 0.35834  
    ## Residuals               0.81858  4                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(FD_FRicz_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              1.00000000              0.08178872              1.00000000 
    ##                 logArea 
    ##              1.00000000

``` r
##### PDisp ANOVAs
#### Environmental
### Surveyed sites
## Temperature
# Interaction
FD_FDisp_am_temp <- aov(Dispersion ~ Site_type * temperature_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_temp)
    ## W = 0.93005, p-value = 0.1736

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_temp) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2   1.938 0.1763
    ##       16

``` r
# Check for outliers
outlierTest(FD_FDisp_am_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## CLM -3.241058          0.0070727      0.13438

``` r
# Plot residuals
plot(FD_FDisp_am_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-13.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FDisp_am_temp)
```

    ##                              Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Site_type                     2 0.00550 0.002748   1.113  0.358
    ## temperature_median            1 0.00234 0.002344   0.950  0.348
    ## Site_type:temperature_median  2 0.00042 0.000208   0.084  0.920
    ## Residuals                    13 0.03209 0.002469

``` r
# ANCOVA
FD_FDisp_am_temp <- aov(Dispersion ~ Site_type + temperature_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_temp)
    ## W = 0.94458, p-value = 0.3182

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_temp) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  1.6097 0.2307
    ##       16

``` r
# Check for outliers
outlierTest(FD_FDisp_am_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## CLM -3.313187           0.005126     0.097394

``` r
# Plot residuals
plot(FD_FDisp_am_temp)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-14.png)<!-- -->

``` r
## Salinity
# Interaction
FD_FDisp_am_sal <- aov(Dispersion ~ Site_type * salinity_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_sal)
    ## W = 0.9724, p-value = 0.8232

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_sal) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  3.8228 0.04395 *
    ##       16                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(FD_FDisp_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## OTM 2.816173            0.01557      0.29583

``` r
# Plot residuals
plot(FD_FDisp_am_sal)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-15.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FDisp_am_sal)
```

    ##                           Df   Sum Sq  Mean Sq F value Pr(>F)  
    ## Site_type                  2 0.005495 0.002748   1.453 0.2695  
    ## salinity_median            1 0.007347 0.007347   3.885 0.0704 .
    ## Site_type:salinity_median  2 0.002921 0.001460   0.772 0.4821  
    ## Residuals                 13 0.024586 0.001891                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# ANCOVA
FD_FDisp_am_sal <- aov(Dispersion ~ Site_type + salinity_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_sal)
    ## W = 0.9827, p-value = 0.969

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_sal) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  2.6953 0.09799 .
    ##       16                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(FD_FDisp_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## OTM 2.807424           0.013975      0.26552

``` r
# Plot residuals
plot(FD_FDisp_am_sal)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-16.png)<!-- -->

``` r
## Oxygen
# Interaction
FD_FDisp_am_oxy <- aov(Dispersion ~ Site_type * oxygen_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_oxy)
    ## W = 0.96558, p-value = 0.6859

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_oxy) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  4.5341 0.02754 *
    ##       16                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(FD_FDisp_am_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## OTM 2.351261           0.036629      0.69594

``` r
# Plot residuals
plot(FD_FDisp_am_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-17.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FDisp_am_oxy)
```

    ##                         Df   Sum Sq  Mean Sq F value Pr(>F)  
    ## Site_type                2 0.005495 0.002748   1.350 0.2933  
    ## oxygen_median            1 0.008330 0.008330   4.092 0.0641 .
    ## Site_type:oxygen_median  2 0.000062 0.000031   0.015 0.9850  
    ## Residuals               13 0.026462 0.002036                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# ANCOVA
FD_FDisp_am_oxy <- aov(Dispersion ~ Site_type + oxygen_median, data = straits_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_oxy)
    ## W = 0.96607, p-value = 0.6961

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_oxy) ~ straits_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)  
    ## group  2   4.486 0.0284 *
    ##       16                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(FD_FDisp_am_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## OTM 2.490991           0.025914      0.49236

``` r
# Plot residuals
plot(FD_FDisp_am_oxy)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-18.png)<!-- -->

``` r
# Summarize ANOVA results and calculate pairwise comparisons
summary(FD_FDisp_am_temp)
```

    ##                    Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Site_type           2 0.00550 0.002748   1.268  0.310
    ## temperature_median  1 0.00234 0.002344   1.082  0.315
    ## Residuals          15 0.03251 0.002167

``` r
FD_FDisp_amp_temp <- emmeans(FD_FDisp_am_temp, pairwise ~ Site_type, adjust = "bonferroni")
FD_FDisp_amp_temp$contrasts
```

    ##  contrast           estimate     SE df t.ratio p.value
    ##  Ocean - Mixed      -0.04812 0.0317 15  -1.518  0.4495
    ##  Ocean - Stratified -0.04287 0.0319 15  -1.342  0.5984
    ##  Mixed - Stratified  0.00525 0.0248 15   0.212  1.0000
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(FD_FDisp_am_sal)
```

    ##                 Df   Sum Sq  Mean Sq F value Pr(>F)  
    ## Site_type        2 0.005495 0.002748   1.498 0.2551  
    ## salinity_median  1 0.007347 0.007347   4.007 0.0637 .
    ## Residuals       15 0.027507 0.001834                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_FDisp_amp_sal <- emmeans(FD_FDisp_am_sal, pairwise ~ Site_type, adjust = "bonferroni")
FD_FDisp_amp_sal$contrasts
```

    ##  contrast           estimate     SE df t.ratio p.value
    ##  Ocean - Mixed       -0.0498 0.0291 15  -1.710  0.3238
    ##  Ocean - Stratified  -0.1038 0.0401 15  -2.586  0.0620
    ##  Mixed - Stratified  -0.0541 0.0331 15  -1.636  0.3680
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(FD_FDisp_am_oxy)
```

    ##               Df   Sum Sq  Mean Sq F value Pr(>F)  
    ## Site_type      2 0.005495 0.002748   1.554 0.2436  
    ## oxygen_median  1 0.008330 0.008330   4.711 0.0464 *
    ## Residuals     15 0.026524 0.001768                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_FDisp_amp_oxy <- emmeans(FD_FDisp_am_oxy, pairwise ~ Site_type, adjust = "bonferroni")
FD_FDisp_amp_oxy$contrasts
```

    ##  contrast           estimate     SE df t.ratio p.value
    ##  Ocean - Mixed       -0.0153 0.0315 15  -0.484  1.0000
    ##  Ocean - Stratified   0.0353 0.0479 15   0.737  1.0000
    ##  Mixed - Stratified   0.0505 0.0326 15   1.548  0.4276
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
# p-values
temp_p_values <- summary(FD_FDisp_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(FD_FDisp_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(FD_FDisp_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000        NA 1.0000000 0.3824842        NA 1.0000000
    ## [8] 0.2786215        NA

``` r
### Mixed lakes
FD_FDisp_env_M_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, straits_sespd_env[mixed_lakes,])
summary(FD_FDisp_env_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = straits_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##       FLK       HLO       LLN       MLN       NLN       NLU       OLO       ULN 
    ## -0.015916 -0.013987 -0.018774  0.022246 -0.008384  0.016583  0.020369 -0.002138 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         1.173514   1.070801   1.096    0.335
    ## salinity_median    -0.019235   0.024894  -0.773    0.483
    ## oxygen_median      -0.031717   0.020136  -1.575    0.190
    ## temperature_median  0.005769   0.015451   0.373    0.728
    ## 
    ## Residual standard error: 0.0227 on 4 degrees of freedom
    ## Multiple R-squared:  0.6021, Adjusted R-squared:  0.3037 
    ## F-statistic: 2.018 on 3 and 4 DF,  p-value: 0.2539

``` r
Anova(FD_FDisp_env_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                        Sum Sq Df F value Pr(>F)
    ## (Intercept)        0.00061886  1  1.2010 0.3347
    ## salinity_median    0.00030763  1  0.5970 0.4828
    ## oxygen_median      0.00127837  1  2.4810 0.1904
    ## temperature_median 0.00007183  1  0.1394 0.7278
    ## Residuals          0.00206108  4

``` r
p_values <- summary(FD_FDisp_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          1.0000000          0.7614088          1.0000000

``` r
### Stratified lakes
FD_FDisp_env_S_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, straits_sespd_env[stratified_lakes,])
summary(FD_FDisp_env_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = straits_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##      BCM      CLM      GLK      HLM      NLK      OTM      SLN      TLN 
    ##  0.03742 -0.06283  0.03624 -0.01857 -0.03218  0.08796 -0.02234 -0.02569 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.206216   0.719516  -0.287    0.789
    ## salinity_median     0.010809   0.011174   0.967    0.388
    ## oxygen_median      -0.001923   0.049266  -0.039    0.971
    ## temperature_median  0.015583   0.015936   0.978    0.384
    ## 
    ## Residual standard error: 0.06507 on 4 degrees of freedom
    ## Multiple R-squared:  0.4238, Adjusted R-squared:  -0.008303 
    ## F-statistic: 0.9808 on 3 and 4 DF,  p-value: 0.4856

``` r
Anova(FD_FDisp_env_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                       Sum Sq Df F value Pr(>F)
    ## (Intercept)        0.0003478  1  0.0821 0.7886
    ## salinity_median    0.0039618  1  0.9356 0.3882
    ## oxygen_median      0.0000065  1  0.0015 0.9707
    ## temperature_median 0.0040490  1  0.9562 0.3835
    ## Residuals          0.0169377  4

``` r
p_values <- summary(FD_FDisp_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
#### Geographical
### Surveyed sites
## Distance
# Interaction
FD_FDisp_am_dist <- aov(Dispersion ~ Site_type * distance_to_ocean_min_m, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_dist)
    ## W = 0.98183, p-value = 0.9422

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_dist) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  3.3094 0.05846 .
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(FD_FDisp_am_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## OTM 2.315624           0.035146      0.73806

``` r
# Plot residuals
plot(FD_FDisp_am_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   19

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-19.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FDisp_am_dist)
```

    ##                                   Df   Sum Sq  Mean Sq F value Pr(>F)  
    ## Site_type                          2 0.008264 0.004132   2.471 0.1161  
    ## distance_to_ocean_min_m            1 0.006815 0.006815   4.075 0.0606 .
    ## Site_type:distance_to_ocean_min_m  2 0.003445 0.001723   1.030 0.3795  
    ## Residuals                         16 0.026754 0.001672                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# ANCOVA
FD_FDisp_am_dist <- aov(Dispersion ~ Site_type + distance_to_ocean_min_m, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_dist)
    ## W = 0.99307, p-value = 0.9998

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_dist) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2   2.139 0.1453
    ##       19

``` r
# Check for outliers
outlierTest(FD_FDisp_am_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## OTM 2.488257           0.023506      0.51712

``` r
# Plot residuals
plot(FD_FDisp_am_dist)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-20.png)<!-- -->

``` r
## Max Depth
# Interaction
FD_FDisp_am_mxd <- aov(Dispersion ~ Site_type * max_depth, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_mxd)
    ## W = 0.92965, p-value = 0.1207

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_mxd) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  3.4962 0.05095 .
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(FD_FDisp_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## CLM -3.596628          0.0026439     0.058165

``` r
# Plot residuals
plot(FD_FDisp_am_mxd)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-21.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FDisp_am_mxd)
```

    ##                     Df   Sum Sq  Mean Sq F value Pr(>F)
    ## Site_type            2 0.008264 0.004132   2.291  0.133
    ## max_depth            1 0.000061 0.000061   0.034  0.856
    ## Site_type:max_depth  2 0.008090 0.004045   2.242  0.139
    ## Residuals           16 0.028863 0.001804

``` r
# ANCOVA
FD_FDisp_am_mxd <- aov(Dispersion ~ Site_type + max_depth, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_mxd)
    ## W = 0.91994, p-value = 0.07581

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_mxd) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  1.0582 0.3667
    ##       19

``` r
# Check for outliers
outlierTest(FD_FDisp_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -2.825612           0.011659       0.2565

``` r
# Plot residuals
plot(FD_FDisp_am_mxd)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-22.png)<!-- -->

``` r
## Log Area
# Interaction
FD_FDisp_am_lga <- aov(Dispersion ~ Site_type * logArea, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_lga)
    ## W = 0.91287, p-value = 0.05419

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_lga) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  3.0074 0.07333 .
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(FD_FDisp_am_lga)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## CLM -3.955096          0.0012703     0.027946

``` r
# Plot residuals
plot(FD_FDisp_am_lga)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-23.png)<!-- -->

``` r
# Summarize the ANOVA results
summary(FD_FDisp_am_lga)
```

    ##                   Df   Sum Sq  Mean Sq F value Pr(>F)  
    ## Site_type          2 0.008264 0.004132   2.434 0.1195  
    ## logArea            1 0.000176 0.000176   0.104 0.7514  
    ## Site_type:logArea  2 0.009671 0.004835   2.848 0.0875 .
    ## Residuals         16 0.027167 0.001698                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# ANCOVA
FD_FDisp_am_lga <- aov(Dispersion ~ Site_type + logArea, data = straits_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(FD_FDisp_am_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(FD_FDisp_am_lga)
    ## W = 0.90519, p-value = 0.03781

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(FD_FDisp_am_lga) ~ straits_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  1.4745 0.2539
    ##       19

``` r
# Check for outliers
outlierTest(FD_FDisp_am_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## SLN -2.92446          0.0094589       0.2081

``` r
# Plot residuals
plot(FD_FDisp_am_lga)
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20&%20ANOVA-24.png)<!-- -->

``` r
# Summarize ANOVA results and calculate pairwise comparisons
summary(FD_FDisp_am_dist)
```

    ##                         Df   Sum Sq  Mean Sq F value Pr(>F)  
    ## Site_type                2 0.008264 0.004132   2.463  0.113  
    ## distance_to_ocean_min_m  1 0.006815 0.006815   4.062  0.059 .
    ## Residuals               18 0.030199 0.001678                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_FDisp_amp_dist <- emmeans(FD_FDisp_am_dist, pairwise ~ Site_type, adjust = "bonferroni")
FD_FDisp_amp_dist$contrasts
```

    ##  contrast           estimate     SE df t.ratio p.value
    ##  Ocean - Mixed       -0.0662 0.0253 18  -2.619  0.0522
    ##  Ocean - Stratified  -0.1030 0.0362 18  -2.844  0.0323
    ##  Mixed - Stratified  -0.0369 0.0263 18  -1.402  0.5335
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(FD_FDisp_am_mxd)
```

    ##             Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Site_type    2 0.00826 0.004132   2.013  0.163
    ## max_depth    1 0.00006 0.000061   0.030  0.865
    ## Residuals   18 0.03695 0.002053

``` r
FD_FDisp_amp_mxd <- emmeans(FD_FDisp_am_mxd, pairwise ~ Site_type, adjust = "bonferroni")
FD_FDisp_amp_mxd$contrasts
```

    ##  contrast           estimate     SE df t.ratio p.value
    ##  Ocean - Mixed       -0.0418 0.0245 18  -1.704  0.3166
    ##  Ocean - Stratified  -0.0437 0.0261 18  -1.675  0.3338
    ##  Mixed - Stratified  -0.0019 0.0249 18  -0.077  1.0000
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(FD_FDisp_am_lga)
```

    ##             Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Site_type    2 0.00826 0.004132   2.019  0.162
    ## logArea      1 0.00018 0.000176   0.086  0.772
    ## Residuals   18 0.03684 0.002047

``` r
FD_FDisp_amp_lga <- emmeans(FD_FDisp_am_lga, pairwise ~ Site_type, adjust = "bonferroni")
FD_FDisp_amp_lga$contrasts
```

    ##  contrast           estimate     SE df t.ratio p.value
    ##  Ocean - Mixed      -0.03803 0.0272 18  -1.398  0.5377
    ##  Ocean - Stratified -0.04274 0.0258 18  -1.654  0.3461
    ##  Mixed - Stratified -0.00472 0.0229 18  -0.206  1.0000
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
# p-values
dist_p_values <- summary(FD_FDisp_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(FD_FDisp_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(FD_FDisp_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.6803005 0.3542655        NA 0.9756095 1.0000000        NA 0.9706062
    ## [8] 1.0000000        NA

``` r
### Mixed lakes
FD_FDisp_geo_M_lm <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, straits_sespd_env[mixed_lakes,])
summary(FD_FDisp_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = straits_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##       FLK       HLO       LLN       MLN       NLN       NLU       OLO       ULN 
    ## -0.012313 -0.014805 -0.016547  0.014215  0.019358 -0.002100  0.020937 -0.008746 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              0.6670211  0.0581647  11.468  0.00033 ***
    ## distance_to_ocean_min_m  0.0001073  0.0003440   0.312  0.77074    
    ## max_depth               -0.0009766  0.0015029  -0.650  0.55123    
    ## logArea                 -0.0097234  0.0077146  -1.260  0.27604    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02086 on 4 degrees of freedom
    ## Multiple R-squared:  0.664,  Adjusted R-squared:  0.412 
    ## F-statistic: 2.635 on 3 and 4 DF,  p-value: 0.1862

``` r
Anova(FD_FDisp_geo_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                           Sum Sq Df  F value  Pr(>F)    
    ## (Intercept)             0.057228  1 131.5103 0.00033 ***
    ## distance_to_ocean_min_m 0.000042  1   0.0972 0.77074    
    ## max_depth               0.000184  1   0.4223 0.55123    
    ## logArea                 0.000691  1   1.5886 0.27604    
    ## Residuals               0.001741  4                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(FD_FDisp_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             0.001320052             1.000000000             1.000000000 
    ##                 logArea 
    ##             1.000000000

``` r
### Stratified lakes
FD_FDisp_geo_S_lm <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, straits_sespd_env[stratified_lakes,])
summary(FD_FDisp_geo_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = straits_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##       BCM       CLM       GLK       HLM       NLK       OTM       SLN       TLN 
    ##  0.028173 -0.051918 -0.005307 -0.040354  0.026877  0.042165  0.007181 -0.006817 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              0.2113637  0.2709585   0.780   0.4789  
    ## distance_to_ocean_min_m -0.0006456  0.0002302  -2.804   0.0486 *
    ## max_depth               -0.0013738  0.0032702  -0.420   0.6960  
    ## logArea                  0.0480282  0.0335199   1.433   0.2252  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.044 on 4 degrees of freedom
    ## Multiple R-squared:  0.7366, Adjusted R-squared:  0.539 
    ## F-statistic: 3.728 on 3 and 4 DF,  p-value: 0.1181

``` r
Anova(FD_FDisp_geo_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                            Sum Sq Df F value  Pr(>F)  
    ## (Intercept)             0.0011781  1  0.6085 0.47894  
    ## distance_to_ocean_min_m 0.0152233  1  7.8632 0.04861 *
    ## max_depth               0.0003417  1  0.1765 0.69600  
    ## logArea                 0.0039746  1  2.0530 0.22519  
    ## Residuals               0.0077441  4                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(FD_FDisp_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               0.1944256               1.0000000 
    ##                 logArea 
    ##               0.9007568

``` r
par(mfrow=c(1,1))
```

### FD alpha with env and geo correlated variables

``` r
# Log Area
FD_alpha_LA_plot <- ggplot(data = straits_sespd_env, mapping = aes(y = pd.obs.z, x = logArea, color = Site_type, fill = Site_type)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = straits_sespd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 16),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Log Area (m"^"2"~")", y="FD alpha z-score", colour = "Site type:", fill = "Site type:", tag = "c")
(FD_alpha_LA_plot <- FD_alpha_LA_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_alpha_LA_plot.jpg", plot = FD_alpha_LA_plot, width = 6.26, height = 6, units = "in")

# Distance from the ocean mean
FD_alpha_D_plot <- ggplot(data = straits_sespd_env, mapping = aes(y = pd.obs.z, x = distance_to_ocean_min_m, color = Site_type, fill = Site_type)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = straits_sespd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 16),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Isolation (m)", y="FD alpha z-score", colour = "Site type:", fill = "Site type:", tag = "b")
(FD_alpha_D_plot <- FD_alpha_D_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_alpha_D_plot.jpg", plot = FD_alpha_D_plot, width = 6.26, height = 6, units = "in")

# Max depth
FD_alpha_MD_plot <- ggplot(data = straits_sespd_env, mapping = aes(y = pd.obs.z, x = max_depth, color = Site_type, fill = Site_type)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = straits_sespd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 16),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Age (m)", y="FD alpha z-score", colour = "Site type:", fill = "Site type:", tag = "a")
(FD_alpha_MD_plot <- FD_alpha_MD_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](FD_analyses_files/figure-gfm/FD%20alpha%20lm%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_alpha_MD_plot.jpg", plot = FD_alpha_MD_plot, width = 6.26, height = 6, units = "in")
```

### FD alpha z-scores & dispersion

``` r
# Functional richness and dispersion
FD_alpha_plot <- ggplot(data = straits_sespd_env, mapping = aes(y = pd.obs.z, x = Dispersion, color = Site_type, fill = Site_type)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = straits_sespd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 14),
    legend.position = "bottom", 
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 14), 
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="FD alpha Dispersion", y="FD alpha z-score", colour = "Site type:", fill = "Site type:", tag = "b")
(FD_alpha_plot <- FD_alpha_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](FD_analyses_files/figure-gfm/z-scores%20&%20dispersion-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_alpha_plot.jpg", plot = FD_alpha_plot, width = 4.88, height = 6, units = "in")
```

## FD beta diversity

### FD beta diversity of individual traits, dispersion, and PERMANOVA

``` r
# Run individual trait Functional Diversity analysis
FD_beta_BS <- BAT::beta(surveyed_sites_lake, straits[, "BodyShapeI", drop = FALSE], abund = FALSE)
(FD_beta_BS_BD <- betadisper(FD_beta_BS$Btotal, Site_type_group))
```

    ## Warning in betadisper(FD_beta_BS$Btotal, Site_type_group): some squared
    ## distances are negative and changed to zero

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_BS$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 3
    ## No. of Negative Eigenvalues: 1
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##    0.08333    0.13057    0.17625 
    ## 
    ## Eigenvalues for PCoA axes:
    ##    PCoA1    PCoA2    PCoA3    PCoA4 
    ##  0.44444  0.35835  0.10304 -0.06896

``` r
(FD_beta_BS_AOV <- anova(FD_beta_BS_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Groups     2 0.02977 0.014887   0.604 0.5568
    ## Residuals 19 0.46828 0.024646

``` r
(FD_beta_BS_THSD <- TukeyHSD(FD_beta_BS_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                        diff        lwr       upr     p adj
    ## Mixed-Ocean      0.04723377 -0.1681593 0.2626269 0.8441142
    ## Stratified-Ocean 0.09291369 -0.1224794 0.3083068 0.5281069
    ## Stratified-Mixed 0.04567992 -0.1537353 0.2450952 0.8313064

``` r
(FD_beta_BS_PM <- adonis2(FD_beta_BS$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_BS$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## Model     2  0.12409 0.14828 1.6539  0.156
    ## Residual 19  0.71278 0.85172              
    ## Total    21  0.83687 1.00000

``` r
FD_beta_OP <- BAT::beta(surveyed_sites_lake, straits[, "OperculumPresent", drop = FALSE], abund = FALSE)
(FD_beta_OP_BD <- betadisper(FD_beta_OP$Btotal, Site_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_OP$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 1
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##  8.095e-18  1.041e-17  6.250e-02 
    ## 
    ## Eigenvalues for PCoA axes:
    ##  PCoA1 
    ## 0.2386

``` r
(FD_beta_OP_AOV <- anova(FD_beta_OP_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.019886 0.0099432  0.8636 0.4375
    ## Residuals 19 0.218750 0.0115132

``` r
(FD_beta_OP_THSD <- TukeyHSD(FD_beta_OP_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                          diff         lwr       upr     p adj
    ## Mixed-Ocean      6.158268e-17 -0.14721474 0.1472147 1.0000000
    ## Stratified-Ocean 6.250000e-02 -0.08471474 0.2097147 0.5384199
    ## Stratified-Mixed 6.250000e-02 -0.07379436 0.1987944 0.4876159

``` r
(FD_beta_OP_PM <- adonis2(FD_beta_OP$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_OP$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## Model     2 0.019886 0.08333 0.8636      1
    ## Residual 19 0.218750 0.91667              
    ## Total    21 0.238636 1.00000

``` r
FD_beta_ML <- BAT::beta(surveyed_sites_lake, straits[, "MaxLengthTL", drop = FALSE], abund = FALSE)
(FD_beta_ML_BD <- betadisper(FD_beta_ML$Btotal, Site_type_group))
```

    ## Warning in betadisper(FD_beta_ML$Btotal, Site_type_group): some squared
    ## distances are negative and changed to zero

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_ML$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 11
    ## No. of Negative Eigenvalues: 4
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##    0.10794    0.16442    0.01531 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 15 eigenvalues)
    ##     PCoA1     PCoA2     PCoA3     PCoA4     PCoA5     PCoA6     PCoA7     PCoA8 
    ## 5.234e-01 4.977e-02 1.766e-02 6.025e-03 1.101e-03 7.605e-04 4.890e-05 1.092e-05

``` r
(FD_beta_ML_AOV <- anova(FD_beta_ML_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Groups     2 0.09036 0.045180  2.5959 0.1008
    ## Residuals 19 0.33068 0.017404

``` r
(FD_beta_ML_THSD <- TukeyHSD(FD_beta_ML_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff        lwr        upr     p adj
    ## Mixed-Ocean       0.05648372 -0.1245180 0.23748540 0.7118635
    ## Stratified-Ocean -0.09262672 -0.2736284 0.08837496 0.4123271
    ## Stratified-Mixed -0.14911044 -0.3166854 0.01846455 0.0864317

``` r
(FD_beta_ML_PM <- adonis2(FD_beta_ML$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_ML$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## Model     2  0.11024 0.18889 2.2124  0.101
    ## Residual 19  0.47339 0.81111              
    ## Total    21  0.58364 1.00000

``` r
FD_beta_T <- BAT::beta(surveyed_sites_lake, straits[, "Troph", drop = FALSE], abund = FALSE)
(FD_beta_T_BD <- betadisper(FD_beta_T$Btotal, Site_type_group))
```

    ## Warning in betadisper(FD_beta_T$Btotal, Site_type_group): some squared
    ## distances are negative and changed to zero

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_T$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 11
    ## No. of Negative Eigenvalues: 1
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##   0.056047   0.008992   0.419319 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 12 eigenvalues)
    ##     PCoA1     PCoA2     PCoA3     PCoA4     PCoA5     PCoA6     PCoA7     PCoA8 
    ## 1.885e+00 3.486e-01 1.322e-01 7.932e-02 2.705e-02 1.592e-02 5.227e-04 5.605e-05

``` r
(FD_beta_T_AOV <- anova(FD_beta_T_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## Groups     2 0.78256 0.39128  92.951 1.543e-10 ***
    ## Residuals 19 0.07998 0.00421                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(FD_beta_T_THSD <- TukeyHSD(FD_beta_T_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff        lwr        upr     p adj
    ## Mixed-Ocean      -0.04705512 -0.1360714 0.04196116 0.3897448
    ## Stratified-Ocean  0.36327195  0.2742557 0.45228823 0.0000000
    ## Stratified-Mixed  0.41032707  0.3279140 0.49274014 0.0000000

``` r
(FD_beta_T_PM <- adonis2(FD_beta_T$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_T$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs     R2      F Pr(>F)   
    ## Model     2  0.98446 0.4003 6.3411  0.002 **
    ## Residual 19  1.47488 0.5997                 
    ## Total    21  2.45934 1.0000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_beta_DMin <- BAT::beta(surveyed_sites_lake, straits[, "DepthMin", drop = FALSE], abund = FALSE)
(FD_beta_DMin_BD <- betadisper(FD_beta_DMin$Btotal, Site_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_DMin$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 2
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##  1.017e-15  8.333e-03  6.250e-02 
    ## 
    ## Eigenvalues for PCoA axes:
    ##    PCoA1    PCoA2 
    ## 0.238315 0.002797

``` r
(FD_beta_DMin_AOV <- anova(FD_beta_DMin_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.01721 0.0086048  0.7343  0.493
    ## Residuals 19 0.22264 0.0117178

``` r
(FD_beta_DMin_THSD <- TukeyHSD(FD_beta_DMin_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff         lwr       upr     p adj
    ## Mixed-Ocean      0.008333333 -0.14018421 0.1568509 0.9888671
    ## Stratified-Ocean 0.062500000 -0.08601755 0.2110175 0.5440766
    ## Stratified-Mixed 0.054166667 -0.08333386 0.1916672 0.5853185

``` r
(FD_beta_DMin_PM <- adonis2(FD_beta_DMin$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_DMin$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## Model     2 0.018472 0.07661 0.7882      1
    ## Residual 19 0.222639 0.92339              
    ## Total    21 0.241111 1.00000

``` r
FD_beta_DMax <- BAT::beta(surveyed_sites_lake, straits_clean[, "DepthMax", drop = FALSE], abund = FALSE)
(FD_beta_DMax_BD <- betadisper(FD_beta_DMax$Btotal, Site_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_DMax$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 5
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##   0.010820   0.062500   0.004024 
    ## 
    ## Eigenvalues for PCoA axes:
    ##     PCoA1     PCoA2     PCoA3     PCoA4     PCoA5 
    ## 2.411e-01 1.029e-03 6.243e-06 2.087e-07 7.113e-08

``` r
(FD_beta_DMax_AOV <- anova(FD_beta_DMax_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.015875 0.0079377  0.6868 0.5153
    ## Residuals 19 0.219595 0.0115576

``` r
(FD_beta_DMax_THSD <- TukeyHSD(FD_beta_DMax_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                          diff         lwr        upr     p adj
    ## Mixed-Ocean       0.051679517 -0.09581938 0.19917841 0.6528984
    ## Stratified-Ocean -0.006796339 -0.15429524 0.14070256 0.9924772
    ## Stratified-Mixed -0.058475855 -0.19503330 0.07808159 0.5328699

``` r
(FD_beta_DMax_PM <- adonis2(FD_beta_DMax$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_DMax$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2    F Pr(>F)
    ## Model     2 0.022436 0.09265 0.97  0.268
    ## Residual 19 0.219727 0.90735            
    ## Total    21 0.242163 1.00000

``` r
FD_beta_TMin <- BAT::beta(surveyed_sites_lake, straits[, "TempPrefMin", drop = FALSE], abund = FALSE)
(FD_beta_TMin_BD <- betadisper(FD_beta_TMin$Btotal, Site_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_TMin$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 10
    ## No. of Negative Eigenvalues: 6
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.1444     0.1580     0.3001 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 16 eigenvalues)
    ##    PCoA1    PCoA2    PCoA3    PCoA4    PCoA5    PCoA6    PCoA7    PCoA8 
    ## 1.497465 0.742113 0.221929 0.080618 0.029435 0.015949 0.003644 0.001314

``` r
(FD_beta_TMin_AOV <- anova(FD_beta_TMin_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value  Pr(>F)  
    ## Groups     2 0.11194 0.055969  4.0674 0.03386 *
    ## Residuals 19 0.26145 0.013760                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(FD_beta_TMin_THSD <- TukeyHSD(FD_beta_TMin_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                        diff          lwr       upr     p adj
    ## Mixed-Ocean      0.01360189 -0.147339848 0.1745436 0.9749408
    ## Stratified-Ocean 0.15563432 -0.005307425 0.3165761 0.0590889
    ## Stratified-Mixed 0.14203242 -0.006970676 0.2910355 0.0633188

``` r
(FD_beta_TMin_PM <- adonis2(FD_beta_TMin$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_TMin$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   1.1730 0.48719 9.0255  0.001 ***
    ## Residual 19   1.2346 0.51281                  
    ## Total    21   2.4076 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_beta_TMax <- BAT::beta(surveyed_sites_lake, straits[, "TempPrefMax", drop = FALSE], abund = FALSE)
(FD_beta_TMax_BD <- betadisper(FD_beta_TMax$Btotal, Site_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_TMax$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 7
    ## No. of Negative Eigenvalues: 4
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.1617     0.2430     0.1507 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 11 eigenvalues)
    ##      PCoA1      PCoA2      PCoA3      PCoA4      PCoA5      PCoA6      PCoA7 
    ##  1.175e+00  7.863e-01  1.565e-01  8.697e-02  2.044e-02  4.424e-03  4.134e-05 
    ##      PCoA8 
    ## -1.723e-02

``` r
(FD_beta_TMax_AOV <- anova(FD_beta_TMax_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Groups     2 0.03950 0.019752  0.3815  0.688
    ## Residuals 19 0.98377 0.051777

``` r
(FD_beta_TMax_THSD <- TukeyHSD(FD_beta_TMax_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff        lwr       upr     p adj
    ## Mixed-Ocean       0.08136072 -0.2308322 0.3935537 0.7878800
    ## Stratified-Ocean -0.01096956 -0.3231625 0.3012234 0.9956177
    ## Stratified-Mixed -0.09233029 -0.3813648 0.1967042 0.7006019

``` r
(FD_beta_TMax_PM <- adonis2(FD_beta_TMax$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_TMax$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)   
    ## Model     2  0.63483 0.30603 4.1894  0.007 **
    ## Residual 19  1.43955 0.69397                 
    ## Total    21  2.07438 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_beta_FP <- BAT::beta(surveyed_sites_lake, straits[, "FeedingPath", drop = FALSE], abund = FALSE)
(FD_beta_FP_BD <- betadisper(FD_beta_FP$Btotal, Site_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_FP$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 1
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##  6.939e-18  1.214e-17  1.250e-01 
    ## 
    ## Eigenvalues for PCoA axes:
    ##  PCoA1 
    ## 0.4545

``` r
(FD_beta_FP_AOV <- anova(FD_beta_FP_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Groups     2 0.07955 0.039773  2.0152 0.1608
    ## Residuals 19 0.37500 0.019737

``` r
(FD_beta_FP_THSD <- TukeyHSD(FD_beta_FP_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                           diff         lwr       upr     p adj
    ## Mixed-Ocean      -1.561251e-17 -0.19274933 0.1927493 1.0000000
    ## Stratified-Ocean  1.250000e-01 -0.06774933 0.3177493 0.2508645
    ## Stratified-Mixed  1.250000e-01 -0.05345121 0.3034512 0.2031343

``` r
(FD_beta_FP_PM <- adonis2(FD_beta_FP$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_FP$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs    R2      F Pr(>F)
    ## Model     2  0.07955 0.175 2.0152  0.313
    ## Residual 19  0.37500 0.825              
    ## Total    21  0.45455 1.000

``` r
FD_beta_RG <- BAT::beta(surveyed_sites_lake, straits[, "RepGuild2", drop = FALSE], abund = FALSE)
(FD_beta_RG_BD <- betadisper(FD_beta_RG$Btotal, Site_type_group))
```

    ## Warning in betadisper(FD_beta_RG$Btotal, Site_type_group): some squared
    ## distances are negative and changed to zero

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_RG$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 4
    ## No. of Negative Eigenvalues: 2
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.1000     0.0750     0.2412 
    ## 
    ## Eigenvalues for PCoA axes:
    ##    PCoA1    PCoA2    PCoA3    PCoA4    PCoA5    PCoA6 
    ##  0.72496  0.51693  0.19241  0.03141 -0.07329 -0.12917

``` r
(FD_beta_RG_AOV <- anova(FD_beta_RG_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value  Pr(>F)  
    ## Groups     2 0.12514 0.062570  3.2564 0.06081 .
    ## Residuals 19 0.36508 0.019215                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(FD_beta_RG_THSD <- TukeyHSD(FD_beta_RG_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff          lwr       upr     p adj
    ## Mixed-Ocean      -0.02500186 -0.215184397 0.1651807 0.9405594
    ## Stratified-Ocean  0.14114897 -0.049033563 0.3313315 0.1701310
    ## Stratified-Mixed  0.16615083 -0.009923979 0.3422256 0.0664123

``` r
(FD_beta_RG_PM <- adonis2(FD_beta_RG$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_RG$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)   
    ## Model     2  0.36130 0.28601 3.8055  0.004 **
    ## Residual 19  0.90194 0.71399                 
    ## Total    21  1.26324 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_beta_PC <- BAT::beta(surveyed_sites_lake, straits[, "ParentalCare", drop = FALSE], abund = FALSE)
(FD_beta_PC_BD <- betadisper(FD_beta_PC$Btotal, Site_type_group))
```

    ## Warning in betadisper(FD_beta_PC$Btotal, Site_type_group): some squared
    ## distances are negative and changed to zero

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_PC$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 4
    ## No. of Negative Eigenvalues: 2
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.1662     0.1562     0.1464 
    ## 
    ## Eigenvalues for PCoA axes:
    ##    PCoA1    PCoA2    PCoA3    PCoA4    PCoA5    PCoA6 
    ##  1.15653  0.38132  0.29734  0.07416 -0.04217 -0.12291

``` r
(FD_beta_PC_AOV <- anova(FD_beta_PC_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Groups     2 0.00134 0.000671  0.0191 0.9811
    ## Residuals 19 0.66900 0.035211

``` r
(FD_beta_PC_THSD <- TukeyHSD(FD_beta_PC_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                          diff        lwr       upr     p adj
    ## Mixed-Ocean      -0.009914027 -0.2673625 0.2475345 0.9947389
    ## Stratified-Ocean -0.019717416 -0.2771659 0.2377311 0.9793691
    ## Stratified-Mixed -0.009803389 -0.2481544 0.2285476 0.9940008

``` r
(FD_beta_PC_PM <- adonis2(FD_beta_PC$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_PC$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs     R2      F Pr(>F)    
    ## Model     2  0.83237 0.4772 8.6714  0.001 ***
    ## Residual 19  0.91191 0.5228                  
    ## Total    21  1.74428 1.0000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
straits_augmented <- rbind(straits, straits[rownames(straits) == "Toxotes_jaculatrix", ])
rownames(straits_augmented)[nrow(straits_augmented)] <- "Dummy_Species"
FD_beta_WP <- BAT::beta(surveyed_sites_lake, straits_augmented[, "WaterPref", drop = FALSE], abund = FALSE)
(FD_beta_WP_BD <- betadisper(FD_beta_WP$Btotal, Site_type_group))
```

    ## Warning in betadisper(FD_beta_WP$Btotal, Site_type_group): some squared
    ## distances are negative and changed to zero

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_WP$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 3
    ## No. of Negative Eigenvalues: 1
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##    0.09722    0.09375    0.18454 
    ## 
    ## Eigenvalues for PCoA axes:
    ##   PCoA1   PCoA2   PCoA3   PCoA4 
    ##  0.4444  0.4017  0.2640 -0.1423

``` r
(FD_beta_WP_AOV <- anova(FD_beta_WP_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Groups     2 0.04064 0.020319  0.7622 0.4804
    ## Residuals 19 0.50652 0.026659

``` r
(FD_beta_WP_THSD <- TukeyHSD(FD_beta_WP_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                          diff        lwr       upr     p adj
    ## Mixed-Ocean      -0.003472222 -0.2274862 0.2205418 0.9991455
    ## Stratified-Ocean  0.087315534 -0.1366985 0.3113295 0.5917495
    ## Stratified-Mixed  0.090787756 -0.1166089 0.2981844 0.5185089

``` r
(FD_beta_WP_PM <- adonis2(FD_beta_WP$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_WP$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)   
    ## Model     2  0.24876 0.25704 3.2867  0.008 **
    ## Residual 19  0.71904 0.74296                 
    ## Total    21  0.96780 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_beta_DS <- BAT::beta(surveyed_sites_lake, straits[, "DorsalSpinesMean", drop = FALSE], abund = FALSE)
(FD_beta_DS_BD <- betadisper(FD_beta_DS$Btotal, Site_type_group))
```

    ## Warning in betadisper(FD_beta_DS$Btotal, Site_type_group): some squared
    ## distances are negative and changed to zero

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_DS$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 10
    ## No. of Negative Eigenvalues: 3
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##    0.06886    0.05665    0.47184 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 13 eigenvalues)
    ##    PCoA1    PCoA2    PCoA3    PCoA4    PCoA5    PCoA6    PCoA7    PCoA8 
    ## 2.174590 0.714139 0.225799 0.095939 0.049400 0.029249 0.012557 0.004629

``` r
(FD_beta_DS_AOV <- anova(FD_beta_DS_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## Groups     2 0.85614 0.42807  38.097 2.246e-07 ***
    ## Residuals 19 0.21349 0.01124                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(FD_beta_DS_THSD <- TukeyHSD(FD_beta_DS_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff        lwr       upr     p adj
    ## Mixed-Ocean      -0.01220565 -0.1576405 0.1332292 0.9752841
    ## Stratified-Ocean  0.40298954  0.2575547 0.5484243 0.0000031
    ## Stratified-Mixed  0.41519519  0.2805487 0.5498417 0.0000007

``` r
(FD_beta_DS_PM <- adonis2(FD_beta_DS$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_DS$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)   
    ## Model     2   1.1927 0.37257 5.6411  0.002 **
    ## Residual 19   2.0086 0.62743                 
    ## Total    21   3.2014 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- c(FD_beta_BS_PM$`Pr(>F)`[1], FD_beta_OP_PM$`Pr(>F)`[1], FD_beta_ML_PM$`Pr(>F)`[1], FD_beta_T_PM$`Pr(>F)`[1], FD_beta_DMin_PM$`Pr(>F)`[1], FD_beta_DMax_PM$`Pr(>F)`[1], FD_beta_TMin_PM$`Pr(>F)`[1], FD_beta_TMax_PM$`Pr(>F)`[1], FD_beta_FP_PM$`Pr(>F)`[1], FD_beta_RG_PM$`Pr(>F)`[1], FD_beta_PC_PM$`Pr(>F)`[1], FD_beta_WP_PM$`Pr(>F)`[1], FD_beta_DS_PM$`Pr(>F)`[1])
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##  [1] 1.000 1.000 1.000 0.026 1.000 1.000 0.013 0.091 1.000 0.052 0.013 0.104
    ## [13] 0.026

``` r
straits_sub_keep <- c("Troph", "TempPrefMin", "ParentalCare", "DorsalSpinesMean")
```

### FD beta diversity distance calculations plus dispersion and PERMANOVA

``` r
### Regular
# Reference
FD_beta_ref_dist <- BAT::beta(presabs_lake, traits[,traits_keep], abund = F)
(FD_beta_ref_BD <- betadisper(FD_beta_ref_dist$Btotal, Site_type_group_ref))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_ref_dist$Btotal, group =
    ## Site_type_group_ref)
    ## 
    ## No. of Positive Eigenvalues: 22
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##  Reference      Ocean      Mixed Stratified 
    ##     0.0000     0.4115     0.3806     0.4737 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 22 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.9723 0.8426 0.6132 0.5455 0.3872 0.3347 0.2472 0.2355

``` r
(FD_beta_ref_AOV <- anova(FD_beta_ref_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value    Pr(>F)    
    ## Groups     3 0.20686 0.068952  9.5639 0.0004609 ***
    ## Residuals 19 0.13698 0.007210                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(FD_beta_ref_THSD <- TukeyHSD(FD_beta_ref_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                             diff         lwr        upr     p adj
    ## Ocean-Reference       0.41148304  0.15360159 0.66936450 0.0013252
    ## Mixed-Reference       0.38055273  0.12731817 0.63378729 0.0023684
    ## Stratified-Reference  0.47374107  0.22050651 0.72697563 0.0002408
    ## Mixed-Ocean          -0.03093031 -0.15987104 0.09801042 0.9054389
    ## Stratified-Ocean      0.06225803 -0.06668270 0.19119876 0.5395373
    ## Stratified-Mixed      0.09318834 -0.02618758 0.21256426 0.1604943

``` r
(FD_beta_ref_PM <- adonis2(FD_beta_ref_dist$Btotal ~ env[,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_ref_dist$Btotal ~ env[, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     3   2.5608 0.38604 3.9823  0.001 ***
    ## Residual 19   4.0726 0.61396                  
    ## Total    22   6.6333 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(FD_beta_ref_PM_pair <- pairwise.adonis(FD_beta_ref_dist$Btotal, env[,19], p.adjust.m = "bonferroni", perm = 999))
```

    ## Set of permutations < 'minperm'. Generating entire set.

    ##                     pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted
    ## 1     Stratified vs Mixed  1 1.3638470 6.386376 0.3132669   0.001      0.006
    ## 2     Stratified vs Ocean  1 1.2705416 5.257780 0.3046614   0.001      0.006
    ## 3 Stratified vs Reference  1 0.6491367 2.500782 0.2632186   0.103      0.618
    ## 4          Mixed vs Ocean  1 0.2773129 1.475364 0.1094860   0.109      0.654
    ## 5      Mixed vs Reference  1 0.6155581 3.674142 0.3442096   0.120      0.720
    ## 6      Ocean vs Reference  1 0.5747831 2.654192 0.3467632   0.141      0.846
    ##   sig
    ## 1   *
    ## 2   *
    ## 3    
    ## 4    
    ## 5    
    ## 6

``` r
# Surveyed sites
FD_beta_dist <- BAT::beta(surveyed_sites_lake, straits[,straits_keep], abund = F)
(FD_beta_BD <- betadisper(FD_beta_dist$Btotal, Site_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_dist$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 21
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.4299     0.4045     0.4490 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 21 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.9071 0.8322 0.5418 0.3840 0.3575 0.2883 0.2370 0.1843

``` r
(FD_beta_AOV <- anova(FD_beta_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.007963 0.0039814   0.465 0.6351
    ## Residuals 19 0.162669 0.0085615

``` r
(FD_beta_THSD <- TukeyHSD(FD_beta_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff        lwr       upr     p adj
    ## Mixed-Ocean      -0.02541832 -0.1523674 0.1015308 0.8680513
    ## Stratified-Ocean  0.01907573 -0.1078734 0.1460248 0.9231386
    ## Stratified-Mixed  0.04449404 -0.0730380 0.1620261 0.6091140

``` r
(FD_beta_PM <- adonis2(FD_beta_dist$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_dist$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2     F Pr(>F)    
    ## Model     2   1.9240 0.31665 4.402  0.001 ***
    ## Residual 19   4.1522 0.68335                 
    ## Total    21   6.0761 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(FD_beta_PM_pair <- pairwise.adonis(FD_beta_dist$Btotal, env[surveyed_sites,19], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.2933447 6.065453 0.3022834   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.2283632 5.202464 0.3024255   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.3169502 1.530099 0.1130885   0.088      0.264

``` r
# A subset of traits that are significant
FD_beta_sub_dist <- BAT::beta(surveyed_sites_lake, straits[,straits_sub_keep], abund = F)
(FD_beta_sub_BD <- betadisper(FD_beta_sub_dist$Btotal, Site_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_sub_dist$Btotal, group = Site_type_group)
    ## 
    ## No. of Positive Eigenvalues: 21
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.3481     0.3321     0.3741 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 21 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 2.2098 0.6058 0.2889 0.2655 0.2315 0.2114 0.1761 0.1560

``` r
(FD_beta_sub_AOV <- anova(FD_beta_sub_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.007171 0.0035854  0.5039  0.612
    ## Residuals 19 0.135205 0.0071161

``` r
(FD_beta_sub_THSD <- TukeyHSD(FD_beta_sub_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff         lwr        upr     p adj
    ## Mixed-Ocean      -0.01603858 -0.13177610 0.09969893 0.9341982
    ## Stratified-Ocean  0.02598228 -0.08975523 0.14171980 0.8373426
    ## Stratified-Mixed  0.04202087 -0.06513125 0.14917298 0.5880639

``` r
(FD_beta_sub_PM <- adonis2(FD_beta_sub_dist$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_sub_dist$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   2.1756 0.43516 7.3189  0.001 ***
    ## Residual 19   2.8239 0.56484                  
    ## Total    21   4.9995 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(FD_beta_sub_PM_pair <- pairwise.adonis(FD_beta_sub_dist$Btotal, env[surveyed_sites,19], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.5753973 10.669131 0.43248912   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.4358446  8.923915 0.42649356   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.1824607  1.327119 0.09958037   0.148      0.444

``` r
# Without LCN
FD_beta_wo_LCN_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_LCN,], straits[,straits_keep], abund = F)
(FD_beta_wo_LCN_BD <- betadisper(FD_beta_wo_LCN_dist$Btotal, Site_type_group[-c(8)]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_wo_LCN_dist$Btotal, group =
    ## Site_type_group[-c(8)])
    ## 
    ## No. of Positive Eigenvalues: 20
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.3890     0.4045     0.4490 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 20 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.9062 0.7709 0.5320 0.3840 0.3128 0.2719 0.2065 0.1842

``` r
(FD_beta_wo_LCN_AOV <- anova(FD_beta_wo_LCN_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq  Mean Sq F value Pr(>F)
    ## Groups     2 0.013338 0.006669  0.8285 0.4527
    ## Residuals 18 0.144899 0.008050

``` r
(FD_beta_wo_LCN_THSD <- TukeyHSD(FD_beta_wo_LCN_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                        diff         lwr       upr     p adj
    ## Mixed-Ocean      0.01546925 -0.11507171 0.1460102 0.9509663
    ## Stratified-Ocean 0.05996329 -0.07057767 0.1905043 0.4841014
    ## Stratified-Mixed 0.04449404 -0.06999796 0.1589860 0.5912362

``` r
(FD_beta_wo_LCN_PM <- adonis2(FD_beta_wo_LCN_dist$Btotal ~ env[surveyed_sites_wo_LCN,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_wo_LCN_dist$Btotal ~ env[surveyed_sites_wo_LCN, 19], permutations = 999)
    ##          Df SumOfSqs    R2      F Pr(>F)    
    ## Model     2   2.0313 0.349 4.8249  0.001 ***
    ## Residual 18   3.7892 0.651                  
    ## Total    20   5.8205 1.000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(FD_beta_wo_LCN_PM_pair <- pairwise.adonis(FD_beta_wo_LCN_dist$Btotal, env[surveyed_sites_wo_LCN,19], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.2933447 6.065453 0.3022834   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.2747220 5.676102 0.3403734   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.4148743 2.149883 0.1634907   0.021      0.063

``` r
# Without TLN and HLM
FD_beta_wo_TLN_HLM_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_TLN_HLM,], straits[,straits_keep], abund = F)
(FD_beta_wo_TLN_HLM_BD <- betadisper(FD_beta_wo_TLN_HLM_dist$Btotal, Site_type_group[-c(5,21)]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = FD_beta_wo_TLN_HLM_dist$Btotal, group =
    ## Site_type_group[-c(5, 21)])
    ## 
    ## No. of Positive Eigenvalues: 19
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.4299     0.4045     0.4050 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 19 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.8781 0.7491 0.5344 0.3839 0.2847 0.2321 0.1873 0.1757

``` r
(FD_beta_wo_TLN_HLM_AOV <- anova(FD_beta_wo_TLN_HLM_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.002671 0.0013353  0.1431 0.8677
    ## Residuals 17 0.158633 0.0093313

``` r
(FD_beta_wo_TLN_HLM_THSD <- TukeyHSD(FD_beta_wo_TLN_HLM_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                           diff        lwr       upr     p adj
    ## Mixed-Ocean      -0.0254183161 -0.1592513 0.1084147 0.8782607
    ## Stratified-Ocean -0.0249373574 -0.1680108 0.1181361 0.8963108
    ## Stratified-Mixed  0.0004809588 -0.1333520 0.1343140 0.9999531

``` r
(FD_beta_wo_TLN_HLM_PM <- adonis2(FD_beta_wo_TLN_HLM_dist$Btotal ~ env[surveyed_sites_wo_TLN_HLM,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = FD_beta_wo_TLN_HLM_dist$Btotal ~ env[surveyed_sites_wo_TLN_HLM, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   2.0491 0.36875 4.9653  0.001 ***
    ## Residual 17   3.5078 0.63125                  
    ## Total    19   5.5569 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(FD_beta_wo_TLN_HLM_PM_pair <- pairwise.adonis(FD_beta_wo_TLN_HLM_dist$Btotal, env[surveyed_sites_wo_TLN_HLM,19], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.4720254 7.546050 0.3860652   0.002      0.006   *
    ## 2 Stratified vs Ocean  1 1.3280216 6.066877 0.3776015   0.004      0.012   .
    ## 3      Mixed vs Ocean  1 0.3169502 1.530099 0.1130885   0.084      0.252

``` r
# Mixed and stratified lakes
FD_beta_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], straits[,straits_keep], abund = F)
# Ocean sites and mixed lakes
FD_beta_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites,], straits[,straits_keep], abund = F)
# Stratified lakes and ocean sites
FD_beta_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites,], straits[,straits_keep], abund = F)
# Stratified lakes
FD_beta_S_dist<- BAT::beta(surveyed_sites_lake[stratified_lakes,], straits[,straits_keep], abund = F)


### Environmental
# Surveyed sites
FD_beta_env_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_env,], straits[,straits_keep], abund = F)
# A subset of traits that are significant
FD_beta_sub_env_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_env,], straits[,straits_sub_keep], abund = F)
# Mixed and stratified lakes
FD_beta_env_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], straits[,straits_keep], abund = F)
# Ocean sites and mixed lakes
FD_beta_env_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites_env,], straits[,straits_keep], abund = F)
# Stratified lakes and ocean sites
FD_beta_env_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites_env,], straits[,straits_keep], abund = F)
# Mixed lakes
FD_beta_env_M_dist<- BAT::beta(surveyed_sites_lake[mixed_lakes,], straits[,straits_keep], abund = F)


### Geographic
# Surveyed sites
FD_beta_geo_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites,], straits[,straits_keep], abund = F)
# A subset of traits that are significant
FD_beta_sub_geo_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites,], straits[,straits_sub_keep], abund = F)
# Mixed and stratified lakes
FD_beta_geo_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], straits[,straits_keep], abund = F)
# Ocean sites and mixed lakes
FD_beta_geo_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites,], straits[,straits_keep], abund = F)
# Stratified lakes and ocean sites
FD_beta_geo_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites,], straits[,straits_keep], abund = F)
# Mixed lakes
FD_beta_geo_M_dist<- BAT::beta(surveyed_sites_lake[mixed_lakes,], straits[,straits_keep], abund = F)
```

### FD beta q values of dispersion statistics

``` r
# Sites
N <- 22
# Categories
k <- 3
# Sum of squares
SS <- FD_beta_AOV$`Sum Sq`[2]
# Mean of squares
MS <- FD_beta_AOV$`Mean Sq`[2]
# q value
q.value <- qtukey(p = 0.95, nmeans = k, df = N - k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

group <- FD_beta_THSD$group

(OM_q <- (abs(group[1,1]))/USE)
```

    ## [1] 0.7193541

``` r
(MS_q <- (abs(group[2,1]))/BSE)
```

    ## [1] 0.58311

``` r
(SO_q <- (abs(group[3,1]))/USE)
```

    ## [1] 1.259209

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx
```

### FD beta trait dispersions

- The dissimilarity between sites of the same site type

``` r
# Dispersion within-groups
FD_beta_BD_dist <- FD_beta_BD$distances
FD_beta_BD_dist <- as.data.frame(FD_beta_BD_dist)
FD_beta_BD_dist$X <- row.names(FD_beta_BD_dist)
FD_beta_BD_dist_env <- merge(FD_beta_BD_dist, env[surveyed_sites,], by = "X", sort = F)

FD_beta_BD_dist_env$Site_type <- factor(FD_beta_BD_dist_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

(FD_beta_BD_dist_env_plot <- ggplot(FD_beta_BD_dist_env, aes(x = Site_type, y = FD_beta_BD_dist, color = Site_type, fill = Site_type)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Site_type, fill = Site_type)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.8,
            width = 0.1) +
    geom_text_repel(data = FD_beta_BD_dist_env, label = FD_beta_BD_dist_env$X, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 16), legend.text = element_text(size = 16),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black")) +
  xlab("Site type") +
  ylab("Distance to Centroid") +
  labs(color = "Site_type", tag = "c"))
```

![](FD_analyses_files/figure-gfm/FD%20beta%20dispersion%20plot-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_beta_BD_dist_plot.jpg", FD_beta_BD_dist_env_plot, width = 4.88, height = 6, units = "in")
```

### FD beta NMDS

``` r
### Regular
# Reference
FD_beta_ref_NMDS <- metaMDS(FD_beta_ref_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Surveyed sites 
FD_beta_NMDS <- metaMDS(FD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
round(FD_beta_NMDS$stress, digits = 2)
FD_beta_rep_NMDS <- metaMDS(FD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
FD_beta_ric_NMDS <- metaMDS(FD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# A subset of traits that are significant
FD_beta_sub_NMDS <- metaMDS(FD_beta_sub_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
FD_beta_MS_NMDS <- metaMDS(FD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
FD_beta_OM_NMDS <- metaMDS(FD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
FD_beta_SO_NMDS <- metaMDS(FD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes
FD_beta_S_NMDS <- metaMDS(FD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(FD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax =
    ## 500, : stress is (nearly) zero: you may have insufficient data

``` r
### Environmental
# Surveyed sites 
FD_beta_env_NMDS <- metaMDS(FD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# A subset of traits that are significant
FD_beta_sub_env_NMDS <- metaMDS(FD_beta_sub_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
FD_beta_env_MS_NMDS <- metaMDS(FD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
FD_beta_env_OM_NMDS <- metaMDS(FD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
FD_beta_env_SO_NMDS <- metaMDS(FD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed lakes
FD_beta_env_M_NMDS <- metaMDS(FD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(FD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
FD_beta_geo_NMDS <- metaMDS(FD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# A subset of traits that are significant
FD_beta_sub_geo_NMDS <- metaMDS(FD_beta_sub_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
FD_beta_geo_MS_NMDS <- metaMDS(FD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
FD_beta_geo_OM_NMDS <- metaMDS(FD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
FD_beta_geo_SO_NMDS <- metaMDS(FD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed lakes
FD_beta_geo_M_NMDS <- metaMDS(FD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(FD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

### Load an envfit analysis function for adjusting p-values

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

### FD beta env and geo correlated variables using envfit

``` r
### Environmental
# Surveyed sites 
# For figure
(FD_beta_env_ef <- envfit(FD_beta_env_NMDS, env[surveyed_sites_env,c(34)], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                                  NMDS1   NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)] 0.99992 0.01303 0.7819  0.031 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_env_efp <- p.adjust.envfit(FD_beta_env_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                                  NMDS1   NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)] 0.99992 0.01303 0.7819  0.031 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# A subset of traits that are significant
# For figure
(FD_beta_sub_env_ef <- envfit(FD_beta_sub_env_NMDS, env[surveyed_sites_env,c(34)], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                                   NMDS1    NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)]  0.91531 -0.40275 0.8141  0.012 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_sub_env_efp <- p.adjust.envfit(FD_beta_sub_env_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                                   NMDS1    NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)]  0.91531 -0.40275 0.8141  0.012 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by Site_type
(FD_beta_env_A_ef <- envfit(FD_beta_env_NMDS, env[surveyed_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.27074 -0.96265 0.1290  0.533  
    ## salinity_median     0.99992  0.01303 0.7819  0.033 *
    ## oxygen_median       0.56289  0.82653 0.5918  0.544  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_env_A_efp <- p.adjust.envfit(FD_beta_env_A_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.27074 -0.96265 0.1290  1.000  
    ## salinity_median     0.99992  0.01303 0.7819  0.099 .
    ## oxygen_median       0.56289  0.82653 0.5918  1.000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
(FD_beta_env_MS_ef <- envfit(FD_beta_env_MS_NMDS, env[mixed_stratified_lakes,environment], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.24071 -0.97060 0.2281  0.343  
    ## salinity_median     0.98126  0.19270 0.7759  0.020 *
    ## oxygen_median       0.47911  0.87776 0.4513  0.624  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_env_MS_efp <- p.adjust.envfit(FD_beta_env_MS_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.24071 -0.97060 0.2281   1.00  
    ## salinity_median     0.98126  0.19270 0.7759   0.06 .
    ## oxygen_median       0.47911  0.87776 0.4513   1.00  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
(FD_beta_env_OM_ef <- envfit(FD_beta_env_OM_NMDS, env[ocean_mixed_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.13800  0.99043 0.1165  0.553  
    ## salinity_median    -0.78183 -0.62349 0.7761  0.017 *
    ## oxygen_median      -0.98751 -0.15753 0.7007  0.058 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_env_OM_efp <- p.adjust.envfit(FD_beta_env_OM_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.13800  0.99043 0.1165  1.000  
    ## salinity_median    -0.78183 -0.62349 0.7761  0.051 .
    ## oxygen_median      -0.98751 -0.15753 0.7007  0.174  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
(FD_beta_env_SO_ef <- envfit(FD_beta_env_SO_NMDS, env[ocean_stratified_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median  0.11111 -0.99381 0.2427  0.227  
    ## salinity_median    -0.99549  0.09487 0.7716  0.027 *
    ## oxygen_median      -0.54284  0.83984 0.5013  0.832  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_env_SO_efp <- p.adjust.envfit(FD_beta_env_SO_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median  0.11111 -0.99381 0.2427  0.681  
    ## salinity_median    -0.99549  0.09487 0.7716  0.081 .
    ## oxygen_median      -0.54284  0.83984 0.5013  1.000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed lakes
(FD_beta_env_M_ef <- envfit(FD_beta_env_M_NMDS, env[mixed_lakes,environment], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)  
    ## temperature_median -0.00016273  1.00000000 0.1062  0.760  
    ## salinity_median     0.00055983 -1.00000000 0.5164  0.185  
    ## oxygen_median       0.00136234 -1.00000000 0.6166  0.077 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_env_M_efp <- p.adjust.envfit(FD_beta_env_M_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median -0.00016273  1.00000000 0.1062  1.000
    ## salinity_median     0.00055983 -1.00000000 0.5164  0.555
    ## oxygen_median       0.00136234 -1.00000000 0.6166  0.231
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes
(FD_beta_S_ef <- envfit(FD_beta_S_NMDS, env[stratified_lakes,c(2,6,8,26,31:32)], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                             NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median       -0.92000  0.39191 0.5767  0.075 .
    ## salinity_median           0.45813  0.88888 0.6026  0.087 .
    ## oxygen_median            -0.24369 -0.96985 0.6771  0.049 *
    ## distance_to_ocean_mean_m -0.98380  0.17928 0.3436  0.413  
    ## max_depth                -0.82907  0.55914 0.5586  0.136  
    ## logArea                  -0.85185  0.52378 0.5527  0.121  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_S_efp <- p.adjust.envfit(FD_beta_S_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median       -0.92000  0.39191 0.5767  0.450
    ## salinity_median           0.45813  0.88888 0.6026  0.522
    ## oxygen_median            -0.24369 -0.96985 0.6771  0.294
    ## distance_to_ocean_mean_m -0.98380  0.17928 0.3436  1.000
    ## max_depth                -0.82907  0.55914 0.5586  0.816
    ## logArea                  -0.85185  0.52378 0.5527  0.726
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
# For figure
(FD_beta_geo_ef <- envfit(FD_beta_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.81970 -0.57280 0.6240  0.084 .
    ## max_depth               -0.87040  0.49234 0.0216  0.995  
    ## logArea                  0.28188  0.95945 0.1375  0.245  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_geo_efp <- p.adjust.envfit(FD_beta_geo_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.81970 -0.57280 0.6240  0.252
    ## max_depth               -0.87040  0.49234 0.0216  1.000
    ## logArea                  0.28188  0.95945 0.1375  0.735
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# A subset of traits that are significant
# For figure
(FD_beta_sub_geo_ef <- envfit(FD_beta_sub_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.89447 -0.44713 0.6680  0.030 *
    ## max_depth               -0.39699  0.91782 0.1065  0.920  
    ## logArea                  0.13969  0.99020 0.1972  0.137  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_sub_geo_efp <- p.adjust.envfit(FD_beta_sub_geo_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.89447 -0.44713 0.6680  0.090 .
    ## max_depth               -0.39699  0.91782 0.1065  1.000  
    ## logArea                  0.13969  0.99020 0.1972  0.411  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by Site_type
(FD_beta_geo_A_ef <- envfit(FD_beta_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.81970 -0.57280 0.6240  0.094 .
    ## max_depth               -0.87040  0.49234 0.0216  0.996  
    ## logArea                  0.28188  0.95945 0.1375  0.225  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_geo_A_efp <- p.adjust.envfit(FD_beta_geo_A_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.81970 -0.57280 0.6240  0.282
    ## max_depth               -0.87040  0.49234 0.0216  1.000
    ## logArea                  0.28188  0.95945 0.1375  0.675
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
(FD_beta_geo_MS_ef <- envfit(FD_beta_geo_MS_NMDS, env[mixed_stratified_lakes,geography], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.83844 -0.54499 0.5385  0.141
    ## max_depth               -0.99994 -0.01098 0.0552  0.970
    ## logArea                 -0.06131  0.99812 0.0122  0.946
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_geo_MS_efp <- p.adjust.envfit(FD_beta_geo_MS_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.83844 -0.54499 0.5385  0.423
    ## max_depth               -0.99994 -0.01098 0.0552  1.000
    ## logArea                 -0.06131  0.99812 0.0122  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
(FD_beta_geo_OM_ef <- envfit(FD_beta_geo_OM_NMDS, env[ocean_mixed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m  0.49144 -0.87091 0.4051  0.046 *
    ## max_depth               -0.97142 -0.23737 0.5902  0.018 *
    ## logArea                 -0.69406  0.71992 0.3102  0.249  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_geo_OM_efp <- p.adjust.envfit(FD_beta_geo_OM_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m  0.49144 -0.87091 0.4051  0.138  
    ## max_depth               -0.97142 -0.23737 0.5902  0.054 .
    ## logArea                 -0.69406  0.71992 0.3102  0.747  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
(FD_beta_geo_SO_ef <- envfit(FD_beta_geo_SO_NMDS, env[ocean_stratified_sites,geography], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m  0.83952 -0.54333 0.7036  0.136
    ## max_depth                0.22386 -0.97462 0.0872  0.852
    ## logArea                 -0.95841  0.28541 0.1510  0.623
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_geo_SO_efp <- p.adjust.envfit(FD_beta_geo_SO_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m  0.83952 -0.54333 0.7036  0.408
    ## max_depth                0.22386 -0.97462 0.0872  1.000
    ## logArea                 -0.95841  0.28541 0.1510  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed lakes
(FD_beta_geo_M_ef <- envfit(FD_beta_geo_M_NMDS, env[mixed_lakes,geography], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                              NMDS1      NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.0036097 -0.9999900 0.3193  0.379  
    ## max_depth                0.0149725  0.9998900 0.7042  0.066 .
    ## logArea                  0.0030400 -1.0000000 0.6749  0.063 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
(FD_beta_geo_M_efp <- p.adjust.envfit(FD_beta_geo_M_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                              NMDS1      NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.0036097 -0.9999900 0.3193  1.000
    ## max_depth                0.0149725  0.9998900 0.7042  0.198
    ## logArea                  0.0030400 -1.0000000 0.6749  0.189
    ## Permutation: free
    ## Number of permutations: 999

### FD beta Mantel correlation tests

``` r
### Environmental
# Surveyed sites
env_dist_t <- dist(scaled_env[surveyed_sites_env,c(1)], method = "euclidean")
(FD_beta_env_mant_t <- mantel(FD_beta_env_dist$Btotal, env_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_dist$Btotal, ydis = env_dist_t, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2958 
    ##       Significance: 0.213 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.352 0.381 0.404 0.433 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_dist_s <- dist(scaled_env[surveyed_sites_env,c(2)], method = "euclidean")
(FD_beta_env_mant_s <- mantel(FD_beta_env_dist$Btotal, env_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_dist$Btotal, ydis = env_dist_s, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.7174 
    ##       Significance: 0.007 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.614 0.658 0.684 0.706 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_dist_o <- dist(scaled_env[surveyed_sites_env,c(3)], method = "euclidean")
(FD_beta_env_mant_o <- mantel(FD_beta_env_dist$Btotal, env_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_dist$Btotal, ydis = env_dist_o, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4174 
    ##       Significance: 0.578 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.572 0.618 0.656 0.702 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
FD_beta_env_mant_pv <- rbind(FD_beta_env_mant_t$signif, FD_beta_env_mant_s$signif, FD_beta_env_mant_o$signif)
FD_beta_env_mant_pv <- FD_beta_env_mant_pv[,1]
(FD_beta_env_mant_pv <- p.adjust(FD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.639 0.021 1.000

``` r
# Mixed and stratified lakes
env_MS_dist_t <- dist(scaled_env[mixed_stratified_lakes,c(1)], method = "euclidean")
(FD_beta_env_MS_mant_t <- mantel(FD_beta_env_MS_dist$Btotal, env_MS_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_t,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2965 
    ##       Significance: 0.231 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.362 0.393 0.440 0.470 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_MS_dist_s <- dist(scaled_env[mixed_stratified_lakes,c(2)], method = "euclidean")
(FD_beta_env_MS_mant_s <- mantel(FD_beta_env_MS_dist$Btotal, env_MS_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_s,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6985 
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.599 0.648 0.680 0.699 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_MS_dist_o <- dist(scaled_env[mixed_stratified_lakes,c(3)], method = "euclidean")
(FD_beta_env_MS_mant_o <- mantel(FD_beta_env_MS_dist$Btotal, env_MS_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_o,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2429 
    ##       Significance: 0.73 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.478 0.534 0.571 0.612 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
FD_beta_env_MS_mant_pv <- rbind(FD_beta_env_MS_mant_t$signif, FD_beta_env_MS_mant_s$signif, FD_beta_env_MS_mant_o$signif)
FD_beta_env_MS_mant_pv <- FD_beta_env_MS_mant_pv[,1]
(FD_beta_env_MS_mant_pv <- p.adjust(FD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.693 0.036 1.000

``` r
# Ocean sites and mixed lakes
env_OM_dist_t <- dist(scaled_env[ocean_mixed_sites_env,c(1)], method = "euclidean")
(FD_beta_env_OM_mant_t <- mantel(FD_beta_env_OM_dist$Btotal, env_OM_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_t,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.07165 
    ##       Significance: 0.635 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.170 0.240 0.316 0.402 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_OM_dist_s <- dist(scaled_env[ocean_mixed_sites_env,c(2)], method = "euclidean")
(FD_beta_env_OM_mant_s <- mantel(FD_beta_env_OM_dist$Btotal, env_OM_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_s,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4251 
    ##       Significance: 0.028 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.314 0.366 0.436 0.521 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_OM_dist_o <- dist(scaled_env[ocean_mixed_sites_env,c(3)], method = "euclidean")
(FD_beta_env_OM_mant_o <- mantel(FD_beta_env_OM_dist$Btotal, env_OM_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_o,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5302 
    ##       Significance: 0.03 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.336 0.420 0.542 0.633 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
FD_beta_env_OM_mant_pv <- rbind(FD_beta_env_OM_mant_t$signif, FD_beta_env_OM_mant_s$signif, FD_beta_env_OM_mant_o$signif)
FD_beta_env_OM_mant_pv <- FD_beta_env_OM_mant_pv[,1]
(FD_beta_env_OM_mant_pv <- p.adjust(FD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.084 0.090

``` r
# Stratified lakes and ocean sites
env_SO_dist_t <- dist(scaled_env[ocean_stratified_sites_env,c(1)], method = "euclidean")
(FD_beta_env_SO_mant_t <- mantel(FD_beta_env_SO_dist$Btotal, env_SO_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_t,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.005267 
    ##       Significance: 0.319 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0997 0.1393 0.1648 0.1938 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_SO_dist_s <- dist(scaled_env[ocean_stratified_sites_env,c(2)], method = "euclidean")
(FD_beta_env_SO_mant_s <- mantel(FD_beta_env_SO_dist$Btotal, env_SO_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_s,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5652 
    ##       Significance: 0.025 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.466 0.526 0.564 0.613 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_SO_dist_o <- dist(scaled_env[ocean_stratified_sites_env,c(3)], method = "euclidean")
(FD_beta_env_SO_mant_o <- mantel(FD_beta_env_SO_dist$Btotal, env_SO_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_o,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6221 
    ##       Significance: 0.269 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.688 0.719 0.758 0.799 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
FD_beta_env_SO_mant_pv <- rbind(FD_beta_env_SO_mant_t$signif, FD_beta_env_SO_mant_s$signif, FD_beta_env_SO_mant_o$signif)
FD_beta_env_SO_mant_pv <- FD_beta_env_SO_mant_pv[,1]
(FD_beta_env_SO_mant_pv <- p.adjust(FD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 0.957 0.075 0.807

``` r
# Mixed lakes
env_M_dist_t <- dist(scaled_env[mixed_lakes,c(1)], method = "euclidean")
(FD_beta_env_M_mant_t <- mantel(FD_beta_env_M_dist$Btotal, env_M_dist_t, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_M_dist$Btotal, ydis = env_M_dist_t,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1609 
    ##       Significance: 0.785 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.222 0.321 0.418 0.670 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_M_dist_s <- dist(scaled_env[mixed_lakes,c(2)], method = "euclidean")
(FD_beta_env_M_mant_s <- mantel(FD_beta_env_M_dist$Btotal, env_M_dist_s, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_M_dist$Btotal, ydis = env_M_dist_s,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1713 
    ##       Significance: 0.193 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.253 0.340 0.429 0.486 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_M_dist_o <- dist(scaled_env[mixed_lakes,c(3)], method = "euclidean")
(FD_beta_env_M_mant_o <- mantel(FD_beta_env_M_dist$Btotal, env_M_dist_o, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_env_M_dist$Btotal, ydis = env_M_dist_o,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1719 
    ##       Significance: 0.136 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.213 0.399 0.506 0.564 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
FD_beta_env_M_mant_pv <- rbind(FD_beta_env_M_mant_t$signif, FD_beta_env_M_mant_s$signif, FD_beta_env_M_mant_o$signif)
FD_beta_env_M_mant_pv <- FD_beta_env_M_mant_pv[,1]
(FD_beta_env_M_mant_pv <- p.adjust(FD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.579 0.408

``` r
# Stratified lakes
env_geo_S_dist <- dist(scaled_env[stratified_lakes,c(1:3,9,14:15)], method = "euclidean")
(FD_beta_S_mant <- mantel(FD_beta_S_dist$Btotal, env_geo_S_dist, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_S_dist$Btotal, ydis = env_geo_S_dist, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2501 
    ##       Significance: 0.122 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.269 0.371 0.424 0.502 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
geo_dist_dmean <- dist(scaled_env[surveyed_sites,c(9)], method = "euclidean")
(FD_beta_geo_mant_dmean <- mantel(FD_beta_geo_dist$Btotal, geo_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_dist$Btotal, ydis = geo_dist_dmean,      method = "spearman", permutations = 999, strata = env[surveyed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4088 
    ##       Significance: 0.485 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.541 0.581 0.614 0.643 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_md <- dist(scaled_env[surveyed_sites,c(14)], method = "euclidean")
(FD_beta_geo_mant_md <- mantel(FD_beta_geo_dist$Btotal, geo_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_dist$Btotal, ydis = geo_dist_md, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.08001 
    ##       Significance: 0.603 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.223 0.255 0.278 0.298 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_la <- dist(scaled_env[surveyed_sites,c(15)], method = "euclidean")
(FD_beta_geo_mant_la <- mantel(FD_beta_geo_dist$Btotal, geo_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_dist$Btotal, ydis = geo_dist_la, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.04732 
    ##       Significance: 0.124 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0523 0.0702 0.0831 0.1019 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
FD_beta_geo_mant_pv <- rbind(FD_beta_geo_mant_dmean$signif, FD_beta_geo_mant_md$signif, FD_beta_geo_mant_la$signif)
FD_beta_geo_mant_pv <- FD_beta_geo_mant_pv[,1]
(FD_beta_geo_mant_pv <- p.adjust(FD_beta_geo_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.372

``` r
# Mixed and stratified lakes 
geo_MS_dist_dmean <- dist(scaled_env[mixed_stratified_lakes,c(9)], method = "euclidean")
(FD_beta_geo_MS_mant_dmean <- mantel(FD_beta_geo_MS_dist$Btotal, geo_MS_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_dmean,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2878 
    ##       Significance: 0.425 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.453 0.508 0.552 0.606 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_md <- dist(scaled_env[mixed_stratified_lakes,c(14)], method = "euclidean")
(FD_beta_geo_MS_mant_md <- mantel(FD_beta_geo_MS_dist$Btotal, geo_MS_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_md,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.0198 
    ##       Significance: 0.701 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.231 0.288 0.328 0.368 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_la <- dist(scaled_env[mixed_stratified_lakes,c(15)], method = "euclidean")
(FD_beta_geo_MS_mant_la <- mantel(FD_beta_geo_MS_dist$Btotal, geo_MS_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_la,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.001 
    ##       Significance: 0.365 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0716 0.1073 0.1421 0.1663 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
FD_beta_geo_MS_mant_pv <- rbind(FD_beta_geo_MS_mant_dmean$signif, FD_beta_geo_MS_mant_md$signif, FD_beta_geo_MS_mant_la$signif)
FD_beta_geo_MS_mant_pv <- FD_beta_geo_MS_mant_pv[,1]
(FD_beta_geo_MS_mant_pv <- p.adjust(FD_beta_geo_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 1 1 1

``` r
# Ocean sites and mixed lakes
geo_OM_dist_dmean <- dist(scaled_env[ocean_mixed_sites,c(9)], method = "euclidean")
(FD_beta_geo_OM_mant_dmean <- mantel(FD_beta_geo_OM_dist$Btotal, geo_OM_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_dmean,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.04424 
    ##       Significance: 0.576 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.175 0.207 0.234 0.258 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_md <- dist(scaled_env[ocean_mixed_sites,c(14)], method = "euclidean")
(FD_beta_geo_OM_mant_md <- mantel(FD_beta_geo_OM_dist$Btotal, geo_OM_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_md,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2662 
    ##       Significance: 0.05 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.181 0.255 0.322 0.381 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_la <- dist(scaled_env[ocean_mixed_sites,c(15)], method = "euclidean")
(FD_beta_geo_OM_mant_la <- mantel(FD_beta_geo_OM_dist$Btotal, geo_OM_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_la,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2378 
    ##       Significance: 0.094 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.235 0.276 0.309 0.354 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
FD_beta_geo_OM_mant_pv <- rbind( FD_beta_geo_OM_mant_dmean$signif, FD_beta_geo_OM_mant_md$signif, FD_beta_geo_OM_mant_la$signif)
FD_beta_geo_OM_mant_pv <- FD_beta_geo_OM_mant_pv[,1]
(FD_beta_geo_OM_mant_pv <- p.adjust(FD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.150 0.282

``` r
# Stratified lakes and ocean sites
geo_SO_dist_dmean <- dist(scaled_env[ocean_stratified_sites,c(9)], method = "euclidean")
(FD_beta_geo_SO_mant_dmean <- mantel(FD_beta_geo_SO_dist$Btotal, geo_SO_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_dmean,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5579 
    ##       Significance: 0.411 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.671 0.704 0.731 0.752 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_SO_dist_md <- dist(scaled_env[ocean_stratified_sites,c(14)], method = "euclidean")
(FD_beta_geo_SO_mant_md <- mantel(FD_beta_geo_SO_dist$Btotal, geo_SO_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_md,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02172 
    ##       Significance: 0.683 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.141 0.194 0.224 0.254 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_SO_dist_la <- dist(scaled_env[ocean_stratified_sites,c(15)], method = "euclidean")
(FD_beta_geo_SO_mant_la <- mantel(FD_beta_geo_SO_dist$Btotal, geo_SO_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_la,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2541 
    ##       Significance: 0.134 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.266 0.296 0.319 0.357 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
FD_beta_geo_SO_mant_pv <- rbind( FD_beta_geo_SO_mant_dmean$signif, FD_beta_geo_SO_mant_md$signif, FD_beta_geo_SO_mant_la$signif)
FD_beta_geo_SO_mant_pv <- FD_beta_geo_SO_mant_pv[,1]
(FD_beta_geo_SO_mant_pv <- p.adjust(FD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.402

``` r
# Stratified lakes and ocean sites
geo_M_dist_dmean <- dist(scaled_env[mixed_lakes,c(9)], method = "euclidean")
(FD_beta_geo_M_mant_dmean <- mantel(FD_beta_geo_M_dist$Btotal, geo_M_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_dmean,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1262 
    ##       Significance: 0.197 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.230 0.399 0.582 0.730 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_M_dist_md <- dist(scaled_env[mixed_lakes,c(14)], method = "euclidean")
(FD_beta_geo_M_mant_md <- mantel(FD_beta_geo_M_dist$Btotal, geo_M_dist_md, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_md,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6274 
    ##       Significance: 0.006 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.278 0.397 0.483 0.559 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_M_dist_la <- dist(scaled_env[mixed_lakes,c(15)], method = "euclidean")
(FD_beta_geo_M_mant_la <- mantel(FD_beta_geo_M_dist$Btotal, geo_M_dist_la, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = FD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_la,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4045 
    ##       Significance: 0.042 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.240 0.373 0.461 0.570 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
FD_beta_geo_M_mant_pv <- rbind( FD_beta_geo_M_mant_dmean$signif, FD_beta_geo_M_mant_md$signif, FD_beta_geo_M_mant_la$signif)
FD_beta_geo_M_mant_pv <- FD_beta_geo_M_mant_pv[,1]
(FD_beta_geo_M_mant_pv <- p.adjust(FD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 0.591 0.018 0.126

### FD beta NMDS ordination plots

``` r
# FD beta total NMDS scores
FD_beta_NMDS_data.scores <- as.data.frame(scores(FD_beta_NMDS))
FD_beta_NMDS_data.scores$Site_type <- env[surveyed_sites,19]
FD_beta_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
FD_beta_NMDS_data.scores$Site_type <- factor(FD_beta_NMDS_data.scores$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

# Get significantly correlated environmental variables
FD_beta_env_ef_coord_cont <- as.data.frame(scores(FD_beta_env_ef, "vectors")) * ordiArrowMul(FD_beta_env_ef)
FD_beta_env_row_names <- c("S")
# Assign the new row names to the data frame
FD_beta_env_ef_coord_cont <- data.frame(row.names = FD_beta_env_row_names, FD_beta_env_ef_coord_cont)

# FD_beta_geo_ef_coord_cont <- as.data.frame(scores(FD_beta_geo_ef, "vectors")) * ordiArrowMul(FD_beta_geo_ef)
# FD_beta_geo_row_names <- c("minD")
# # Assign the new row names to the data frame
# FD_beta_geo_ef_coord_cont <- data.frame(row.names = FD_beta_geo_row_names, FD_beta_geo_ef_coord_cont)

# Plot NMDS ordination of FD beta total with CI = 0.95
FD_beta_ef_plot <- ggplot(data = FD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = FD_beta_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = FD_beta_env_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = FD_beta_env_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
            color = env_cont, label = row.names(FD_beta_env_ef_coord_cont), size = 5) +
  geom_text_repel(data = FD_beta_NMDS_data.scores, label = FD_beta_NMDS_data.scores$Lakes, 
                  size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 16), 
        legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_rect(fill = NA), 
        axis.text = element_text(color = "black", size = 16), 
        legend.key = element_blank()) +
  annotate("text", x = -0.6, y = 0.5, size = 5,
           label = paste("Stress: ", round(FD_beta_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81))
(FD_beta_ef_plot <- FD_beta_ef_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](FD_analyses_files/figure-gfm/FD%20beta%20NMDS%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_beta_NMDS.jpg", FD_beta_ef_plot, width = 6, height = 6, units = "in")


# Plot NMDS ordination of FD beta total with CI = 0.90
FD_beta_ef_plot_CI90 <- ggplot(data = FD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.90) +
  geom_point(data = FD_beta_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = FD_beta_env_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = FD_beta_env_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
            color = env_cont, label = row.names(FD_beta_env_ef_coord_cont), size = 5) +
  geom_text_repel(data = FD_beta_NMDS_data.scores, label = FD_beta_NMDS_data.scores$Lakes, 
                  size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 16), 
        legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_rect(fill = NA), 
        axis.text = element_text(color = "black", size = 16), 
        legend.key = element_blank()) +
  annotate("text", x = -0.6, y = 0.5, size = 5,
           label = paste("Stress: ", round(FD_beta_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81))
(FD_beta_ef_plot_CI90 <- FD_beta_ef_plot_CI90 + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](FD_analyses_files/figure-gfm/FD%20beta%20NMDS%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_beta_NMDS_CI90.jpg", FD_beta_ef_plot_CI90, width = 6, height = 6, units = "in")


# FD beta total NMDS scores for the subset of traits
FD_beta_sub_NMDS_data.scores <- as.data.frame(scores(FD_beta_sub_NMDS))
FD_beta_sub_NMDS_data.scores$Site_type <- env[surveyed_sites,19]
FD_beta_sub_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
FD_beta_sub_NMDS_data.scores$Site_type <- factor(FD_beta_sub_NMDS_data.scores$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

# Get significantly correlated environmental variables
FD_beta_sub_env_ef_coord_cont <- as.data.frame(scores(FD_beta_sub_env_ef, "vectors")) * ordiArrowMul(FD_beta_sub_env_ef)
FD_beta_sub_env_row_names <- c("S")
# Assign the new row names to the data frame
FD_beta_sub_env_ef_coord_cont <- data.frame(row.names = FD_beta_sub_env_row_names, FD_beta_sub_env_ef_coord_cont)

# Plot NMDS ordination of FD beta total with the subset of traits and with CI = 0.95
FD_beta_sub_ef_plot <- ggplot(data = FD_beta_sub_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = FD_beta_sub_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = FD_beta_sub_env_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = FD_beta_sub_env_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
            color = env_cont, label = row.names(FD_beta_sub_env_ef_coord_cont), size = 5) +
  geom_text_repel(data = FD_beta_sub_NMDS_data.scores, label = FD_beta_sub_NMDS_data.scores$Lakes, 
                  size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 16), 
        legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_rect(fill = NA), 
        axis.text = element_text(color = "black", size = 16), 
        legend.key = element_blank()) +
  annotate("text", x = -0.6, y = 0.5, size = 5,
           label = paste("Stress: ", round(FD_beta_sub_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1,1)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1,1))
(FD_beta_sub_ef_plot <- FD_beta_sub_ef_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](FD_analyses_files/figure-gfm/FD%20beta%20NMDS%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_beta_sub_NMDS.jpg", FD_beta_sub_ef_plot, width = 6, height = 6, units = "in")


# FD beta replacement NMDS scores
FD_beta_rep_NMDS_data.scores <- as.data.frame(scores(FD_beta_rep_NMDS))
FD_beta_rep_NMDS_data.scores$Site_type <- env[surveyed_sites,19]
FD_beta_rep_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
FD_beta_rep_NMDS_data.scores$Site_type <- factor(FD_beta_rep_NMDS_data.scores$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of FD beta replacement with CI = 0.95
FD_beta_rep_plot <- ggplot(data = FD_beta_rep_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = FD_beta_rep_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = FD_beta_rep_NMDS_data.scores, label = FD_beta_rep_NMDS_data.scores$Lakes, 
                  size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 14), 
        legend.position = "bottom", 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_rect(fill = NA), 
        axis.text = element_text(color = "black", size = 14), 
        legend.key = element_blank()) +
  annotate("text", x = -0.8, y = 0.5, size = 5,
           label = paste("Stress: ", round(FD_beta_rep_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "a") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(FD_beta_rep_plot <- FD_beta_rep_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](FD_analyses_files/figure-gfm/FD%20beta%20NMDS%20plots-4.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_beta_rep_NMDS.jpg", FD_beta_rep_plot, width = 4.88, height = 6, units = "in")


# FD beta richness NMDS scores
FD_beta_ric_NMDS_data.scores <- as.data.frame(scores(FD_beta_ric_NMDS))
FD_beta_ric_NMDS_data.scores$Site_type <- env[surveyed_sites,19]
FD_beta_ric_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
FD_beta_ric_NMDS_data.scores$Site_type <- factor(FD_beta_ric_NMDS_data.scores$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of FD beta richness with CI = 0.95
FD_beta_ric_plot <- ggplot(data = FD_beta_ric_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = FD_beta_ric_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = FD_beta_ric_NMDS_data.scores, label = FD_beta_ric_NMDS_data.scores$Lakes, 
                  size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 14), 
        legend.position = "bottom", 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_rect(fill = NA), 
        axis.text = element_text(color = "black", size = 14), 
        legend.key = element_blank()) +
  annotate("text", x = -0.8, y = 0.5, size = 5,
           label = paste("Stress: ", round(FD_beta_ric_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "b") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(FD_beta_ric_plot <- FD_beta_ric_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](FD_analyses_files/figure-gfm/FD%20beta%20NMDS%20plots-5.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_beta_ric_NMDS.jpg", FD_beta_ric_plot, width = 4.88, height = 6, units = "in")


# FD beta ref NMDS scores
FD_beta_ref_NMDS_data.scores <- as.data.frame(scores(FD_beta_ref_NMDS))
FD_beta_ref_NMDS_data.scores$Site_type <- env[,19]
FD_beta_ref_NMDS_data.scores$Lakes <- env[,1]
FD_beta_ref_NMDS_data.scores$Site_type <- factor(FD_beta_ref_NMDS_data.scores$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of FD beta ref with CI = 0.95
FD_beta_ref_plot <- ggplot(data = FD_beta_ref_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = FD_beta_ref_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = FD_beta_ref_NMDS_data.scores, label = FD_beta_ref_NMDS_data.scores$Lakes, 
                  size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 14), 
        legend.position = "right", 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_rect(fill = NA), 
        axis.text = element_text(color = "black", size = 14), 
        legend.key = element_blank()) +
  annotate("text", x = -0.8, y = 0.5, size = 5,
           label = paste("Stress: ", round(FD_beta_ref_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "c") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(FD_beta_ref_plot <- FD_beta_ref_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_path()`).

![](FD_analyses_files/figure-gfm/FD%20beta%20NMDS%20plots-6.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_beta_ref_NMDS.jpg", FD_beta_ref_plot, width = 6.26, height = 6, units = "in")
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_path()`).

### FD beta dendrogram

``` r
# Cluster communities using average-linkage algorithm
FD_beta_dist_clust <- hclust(FD_beta_dist$Btotal, method = "average")

# The following will create a dendrogram of FD_beta putting things in one dimensional space
# Open a jpg device
png("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_beta_dendrogram.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

# Create your plot using the plot() function
plot(FD_beta_dist_clust, 
     xlab = "Sites",
     ylab = "FD beta",
     main = "",
     sub = "",
     pch = 20,
     col= "black")

# Close the jpg device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# The following will graph and identify the best supported number of clusters of the dendrogram
# Open a jpg device
png("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_beta_dg_silhouette.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

Si <- numeric(nrow(presabs_lake[surveyed_sites,]))
for (k in 2:(nrow(presabs_lake[surveyed_sites,])-1))
{
  sil<-silhouette(cutree(FD_beta_dist_clust, k=k), FD_beta_dist$Btotal)
  Si[k] <- summary(sil)$avg.width
}
k.best<-which.max(Si)
plot(1:nrow(presabs_lake[surveyed_sites,]), Si, type = "h", main = "Silhouette", xlab = "K", ylab = "Width")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col = "red", font = 2, col.axis = "red")

# Close the jpg device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# Partitioning around mediods uses distances matrix to found how points according to silhouette are grouped around a centroid. The greater the number, the greater the distance from the centroid.
# Open a jpg device
png("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/FD_beta_dg_partitioning.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

FD_beta_pam <- pam(FD_beta_dist$Btotal,2)
plot(FD_beta_pam)

# Close the jpg device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### FD beta outliers

``` r
# Determine NMDS1 outliers
ordered_FD_beta_NMDS1_data.scores <- FD_beta_NMDS_data.scores[order(FD_beta_NMDS_data.scores$NMDS1), ]

outlier_FD_beta_NMDS1 <- ordered_FD_beta_NMDS1_data.scores %>%
  group_by(Site_type) %>%
  mutate(
    Q1 = quantile(NMDS1, 0.25),
    Q3 = quantile(NMDS1, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = NMDS1 < lower_bound | NMDS1 > upper_bound
  )

outlier_FD_beta_NMDS1 <- as.data.frame(outlier_FD_beta_NMDS1)

row.names(outlier_FD_beta_NMDS1) <- outlier_FD_beta_NMDS1$X

# Create the plot
(outlier_FD_beta_NMDS1_plot <- ggplot(outlier_FD_beta_NMDS1, aes(x = Site_type, y = NMDS1, fill = Site_type)) +
  geom_violin(trim = FALSE) +
  geom_jitter(aes(color = is_outlier), width = 0.2, size = 4, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "purple")) +
  scale_fill_manual(values = custom_colors) +
  theme(text = element_text(size = 16),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(color = "black", size = 16)) +
  guides(fill = "none") + 
  labs(y = "NMDS1 Distances", x = "Site type", color = "Outlier:", tag = "a"))
```

![](FD_analyses_files/figure-gfm/FD%20beta%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/outlier_FD_beta_NMDS1.jpg", outlier_FD_beta_NMDS1_plot, width = 6.26, height = 6, units = "in")


# Determine NMDS2 outliers
ordered_FD_beta_NMDS2_data.scores <- FD_beta_NMDS_data.scores[order(FD_beta_NMDS_data.scores$NMDS2), ]

outlier_FD_beta_NMDS2 <- ordered_FD_beta_NMDS2_data.scores %>%
  group_by(Site_type) %>%
  mutate(
    Q1 = quantile(NMDS2, 0.25),
    Q3 = quantile(NMDS2, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = NMDS2 < lower_bound | NMDS2 > upper_bound
  )

outlier_FD_beta_NMDS2 <- as.data.frame(outlier_FD_beta_NMDS2)

row.names(outlier_FD_beta_NMDS2) <- outlier_FD_beta_NMDS2$X

# Create the plot
(outlier_FD_beta_NMDS2_plot <- ggplot(outlier_FD_beta_NMDS2, aes(x = Site_type, y = NMDS2, fill = Site_type)) +
  geom_violin(trim = FALSE) +
  geom_jitter(aes(color = is_outlier), width = 0.2, size = 4, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "purple")) +
  scale_fill_manual(values = custom_colors) +
  theme(text = element_text(size = 16),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(color = "black", size = 16)) +
  guides(fill = "none") + 
  labs(y = "NMDS2 Distances", x = "Site type", color = "Outlier:", tag = "b"))
```

![](FD_analyses_files/figure-gfm/FD%20beta%20outliers-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/outlier_FD_beta_NMDS2.jpg", outlier_FD_beta_NMDS2_plot, width = 6.26, height = 6, units = "in")


# Make FD_beta distances into a matrix
FD_beta_dist_matrix <- as.matrix(FD_beta_dist$Btotal)

# Set the diagonal elements to NA
diag(FD_beta_dist_matrix) <- NA

# Get the labels of the sites
sites <- attr(FD_beta_dist_matrix, "Labels")

# Calculate the mean dist for each group
FD_beta_mean_dist <- aggregate(FD_beta_dist_matrix, by = list(Site_type_group), FUN = mean, na.rm = TRUE)

# Rename rows, delete column, and transpose matrix and make it a data frame
row.names(FD_beta_mean_dist) <- FD_beta_mean_dist$Group.1
FD_beta_mean_dist <- FD_beta_mean_dist[,-1] 
FD_beta_mean_dist <- as.data.frame(t(FD_beta_mean_dist))

# Add Site_type column
FD_beta_mean_dist$Site_type <- env[surveyed_sites,19]
FD_beta_mean_dist$Site_type <- factor(FD_beta_mean_dist$Site_type, levels = c("Ocean", "Mixed", "Stratified"))


# Determine ocean sites outliers based on distance from other ocean sites
outlier_FD_beta_mean_dist_ocean <- FD_beta_mean_dist[ocean_sites,] %>%
  group_by(Site_type) %>%
  mutate(
    Q1 = quantile(Ocean, 0.25),
    Q3 = quantile(Ocean, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = Ocean < lower_bound | Ocean > upper_bound
  )

outlier_FD_beta_mean_dist_ocean <- as.data.frame(outlier_FD_beta_mean_dist_ocean)
row.names(outlier_FD_beta_mean_dist_ocean) <- outlier_FD_beta_mean_dist_ocean$X
outlier_FD_beta_mean_dist_ocean$Distances <- outlier_FD_beta_mean_dist_ocean$Ocean
outlier_FD_beta_mean_dist_ocean <- outlier_FD_beta_mean_dist_ocean[,-c(1:3)]


# Determine mixed lake outliers based on distance from other mixed lakes
outlier_FD_beta_mean_dist_mixed <- FD_beta_mean_dist[mixed_lakes,] %>%
  group_by(Site_type) %>%
  mutate(
    Q1 = quantile(Mixed, 0.25),
    Q3 = quantile(Mixed, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = Mixed < lower_bound | Mixed > upper_bound
  )

outlier_FD_beta_mean_dist_mixed <- as.data.frame(outlier_FD_beta_mean_dist_mixed)
row.names(outlier_FD_beta_mean_dist_mixed) <- outlier_FD_beta_mean_dist_mixed$X
outlier_FD_beta_mean_dist_mixed$Distances <- outlier_FD_beta_mean_dist_mixed$Mixed
outlier_FD_beta_mean_dist_mixed <- outlier_FD_beta_mean_dist_mixed[,-c(1:3)]


# Determine stratified lake outliers based on distance from other stratified lakes
outlier_FD_beta_mean_dist_stratified <- FD_beta_mean_dist[stratified_lakes,] %>%
  group_by(Site_type) %>%
  mutate(
    Q1 = quantile(Stratified, 0.25),
    Q3 = quantile(Stratified, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = Stratified < lower_bound | Stratified > upper_bound
  )

outlier_FD_beta_mean_dist_stratified <- as.data.frame(outlier_FD_beta_mean_dist_stratified)
row.names(outlier_FD_beta_mean_dist_stratified) <- outlier_FD_beta_mean_dist_stratified$X
outlier_FD_beta_mean_dist_stratified$Distances <- outlier_FD_beta_mean_dist_stratified$Stratified
outlier_FD_beta_mean_dist_stratified <- outlier_FD_beta_mean_dist_stratified[,-c(1:3)]

outlier_FD_beta_mean_dist <- rbind(outlier_FD_beta_mean_dist_ocean, outlier_FD_beta_mean_dist_mixed, outlier_FD_beta_mean_dist_stratified)

(outlier_FD_beta_mean_dist_plot <- ggplot(outlier_FD_beta_mean_dist, aes(x = Site_type, y = Distances, fill = Site_type)) +
  geom_violin(trim = FALSE) +
  geom_jitter(aes(color = is_outlier), width = 0.2, size = 4, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "purple")) +
  scale_fill_manual(values = custom_colors) +
  theme(text = element_text(size = 16),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(color = "black", size = 16)) +
  guides(fill = "none") + 
  labs(y = "Average Distance from other sites by site type", x = "Site type", color = "Outlier:", tag = "c"))
```

![](FD_analyses_files/figure-gfm/FD%20beta%20outliers-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/outlier_FD_beta_mean_dist.jpg", outlier_FD_beta_mean_dist_plot, width = 6.26, height = 6, units = "in")
```

### Package and version info

``` r
sessionInfo()
```

    ## R version 4.4.3 (2025-02-28)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Ventura 13.6.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
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
    ##  [1] MASS_7.3-65          emmeans_1.10.7       pairwiseAdonis_0.4.1
    ##  [4] cluster_2.1.8        BAT_2.9.6            caret_7.0-1         
    ##  [7] ggrepel_0.9.6        ggplot2_3.5.1        picante_1.8.2       
    ## [10] nlme_3.1-167         vegan_2.6-10         lattice_0.22-6      
    ## [13] permute_0.9-7        car_3.1-3            carData_3.0-5       
    ## [16] tidyr_1.3.1          phytools_2.4-4       maps_3.4.2.1        
    ## [19] ape_5.8-1            reshape2_1.4.4       stringr_1.5.1       
    ## [22] dplyr_1.1.4          knitr_1.49          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rstudioapi_0.17.1       magrittr_2.0.3          TH.data_1.1-3          
    ##   [4] estimability_1.5.1      farver_2.1.2            rmarkdown_2.29         
    ##   [7] ragg_1.3.3              vctrs_0.6.5             base64enc_0.1-3        
    ##  [10] terra_1.8-29            polspline_1.1.25        htmltools_0.5.8.1      
    ##  [13] progress_1.2.3          DEoptim_2.2-8           Formula_1.2-5          
    ##  [16] pROC_1.18.5             parallelly_1.42.0       pracma_2.4.4           
    ##  [19] KernSmooth_2.23-26      htmlwidgets_1.6.4       plyr_1.8.9             
    ##  [22] sandwich_3.1-1          palmerpenguins_0.1.1    zoo_1.8-13             
    ##  [25] lubridate_1.9.4         igraph_2.1.4            lifecycle_1.0.4        
    ##  [28] iterators_1.0.14        pkgconfig_2.0.3         Matrix_1.7-2           
    ##  [31] R6_2.6.1                fastmap_1.2.0           future_1.34.0          
    ##  [34] magic_1.6-1             digest_0.6.37           numDeriv_2016.8-1.1    
    ##  [37] colorspace_2.1-1        textshaping_1.0.0       Hmisc_5.2-2            
    ##  [40] pdist_1.2.1             labeling_0.4.3          clusterGeneration_1.3.8
    ##  [43] timechange_0.3.0        abind_1.4-8             mgcv_1.9-1             
    ##  [46] compiler_4.4.3          proxy_0.4-27            withr_3.0.2            
    ##  [49] doParallel_1.0.17       backports_1.5.0         htmlTable_2.4.3        
    ##  [52] optimParallel_1.0-2     quantreg_6.00           lava_1.8.1             
    ##  [55] scatterplot3d_0.3-44    ModelMetrics_1.2.2.2    tools_4.4.3            
    ##  [58] foreign_0.8-88          future.apply_1.11.3     nnet_7.3-20            
    ##  [61] glue_1.8.0              quadprog_1.5-8          grid_4.4.3             
    ##  [64] checkmate_2.3.2         generics_0.1.3          recipes_1.1.1          
    ##  [67] gtable_0.3.6            class_7.3-23            data.table_1.17.0      
    ##  [70] hms_1.1.3               foreach_1.5.2           pillar_1.10.1          
    ##  [73] splines_4.4.3           survival_3.8-3          SparseM_1.84-2         
    ##  [76] ks_1.14.3               tidyselect_1.2.1        rms_7.0-0              
    ##  [79] gridExtra_2.3           stats4_4.4.3            xfun_0.51              
    ##  [82] expm_1.0-0              hardhat_1.4.1           timeDate_4041.110      
    ##  [85] proto_1.0.0             stringi_1.8.4           yaml_2.3.10            
    ##  [88] evaluate_1.0.3          codetools_0.2-20        tibble_3.2.1           
    ##  [91] cli_3.6.4               rpart_4.1.24            nls2_0.3-4             
    ##  [94] systemfonts_1.2.1       xtable_1.8-4            geometry_0.5.2         
    ##  [97] munsell_0.5.1           Rcpp_1.0.14             globals_0.16.3         
    ## [100] coda_0.19-4.1           fastcluster_1.2.6       MatrixModels_0.5-3     
    ## [103] gower_1.0.2             prettyunits_1.2.0       mclust_6.1.1           
    ## [106] listenv_0.9.1           phangorn_2.12.1         mvtnorm_1.3-3          
    ## [109] ipred_0.9-15            scales_1.3.0            prodlim_2024.06.25     
    ## [112] e1071_1.7-16            purrr_1.0.4             crayon_1.5.3           
    ## [115] combinat_0.0-8          rlang_1.1.5             fastmatch_1.1-6        
    ## [118] multcomp_1.4-28         mnormt_2.1.1            hypervolume_3.1.5

#### straits z score vs p value plots

``` r
straits_sesmpd_plot <- ggplot(data = straits_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Site_type, fill = Site_type)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = straits_sesmpd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 16),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "p-value", x= "z-score", colour = "Site type:", fill = "Site type:", tag = "D") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmpd_plot <- straits_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmpd_plot
```

![](FD_analyses_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/straits_sesmpd_plot_z.jpg", straits_sesmpd_plot, width = 5, height = 4, units = "in")
```

    ## Warning: ggrepel: 3 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
straits_sesmntd_plot <- ggplot(data = straits_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Site_type, fill = Site_type)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = straits_sesmntd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 16),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "p-value", x= "z-score", colour = "Site type:", fill = "Site type:", tag = "E") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmntd_plot <- straits_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmntd_plot
```

![](FD_analyses_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/straits_sesmntd_plot_z.jpg", straits_sesmntd_plot, width = 5, height = 4, units = "in")

straits_sespd_plot <- ggplot(data = straits_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Site_type, fill = Site_type)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = straits_sespd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 16),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "p-value", x= "z-score", colour = "Site type:", fill = "Site type:", tag = "F") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sespd_plot <- straits_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sespd_plot
```

![](FD_analyses_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/straits_sespd_plot_z.jpg", straits_sespd_plot, width = 5, height = 4, units = "in")
```

#### straits z score vs obs value plots

``` r
straits_sesmpd_plot <- ggplot(data = straits_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Site_type, fill = Site_type)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 16),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed MPD", x= "z-score", colour = "Site type:", fill = "Site type:") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmpd_plot <- straits_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmpd_plot
```

![](FD_analyses_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/straits_sesmpd_plot_obs+z.jpg", straits_sesmpd_plot, width = 5, height = 4, units = "in")

straits_sesmntd_plot <- ggplot(data = straits_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Site_type, fill = Site_type)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 16),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed MNTD", x= "z-score", colour = "Site type:", fill = "Site type:") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmntd_plot <- straits_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmntd_plot
```

![](FD_analyses_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/straits_sesmntd_plot_obs+z.jpg", straits_sesmntd_plot, width = 5, height = 4, units = "in")

straits_sespd_plot <- ggplot(data = straits_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Site_type, fill = Site_type)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 16),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed FD", x= "z-score", colour = "Site type:", fill = "Site type:") +
  ylim(c(0,6)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=6, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sespd_plot <- straits_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sespd_plot
```

    ## Warning: Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](FD_analyses_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD/straits_sespd_plot_obs+z.jpg", straits_sespd_plot, width = 5, height = 4, units = "in")
```

    ## Warning: Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

#### Phylogenetic Signal

``` r
phylo_tree <- tree

check_diff <- setdiff(phylo_tree$tip.label, row.names(traits))
check_diff

check_diff <- setdiff(row.names(traits), phylo_tree$tip.label)
check_diff

# Get the tip labels of the phylogenetic tree
tip_labels <- phylo_tree$tip.label

# Order the rows of the trait matrix to match the tip labels
ordered_trait_matrix <- traits[tip_labels, ]

#Determine phylogenetic signal for each max length
ML <- ordered_trait_matrix$MaxLengthTL
ML_psig <- phylosig(phylo_tree, ML, method = "lambda",test=TRUE)
ML_psig

#Determine phylogenetic signal for each trophic level
TL <- ordered_trait_matrix$Troph
TL_psig <- phylosig(phylo_tree, TL, method = "lambda",test=TRUE)
TL_psig

#Determine phylogenetic signal for each depth min
DMin <- ordered_trait_matrix$DepthMin
DMin_psig <- phylosig(phylo_tree, DMin , method = "lambda",test=TRUE)
DMin_psig

#Determine phylogenetic signal for each depth max
DMax <- ordered_trait_matrix$DepthMax
DMax_psig <- phylosig(phylo_tree, DMax , method = "lambda",test=TRUE)
DMax_psig

#Determine phylogenetic signal for each temperature min
TMin <- ordered_trait_matrix$TempPrefMin
TMin_psig <- phylosig(phylo_tree, TMin , method = "lambda",test=TRUE)
TMin_psig

#Determine phylogenetic signal for each temperature max
TMax <- ordered_trait_matrix$TempPrefMax
TMax_psig <- phylosig(phylo_tree, TMax , method = "lambda",test=TRUE)
TMax_psig

#Determine phylogenetic signal for each dorsal spines
DS <- ordered_trait_matrix$DorsalSpinesMean
DS_psig <- phylosig(phylo_tree, DS , method = "lambda",test=TRUE)
DS_psig
```

#### Subset of traits Functional Diversity MNTD, MPD, and PD

``` r
# Trait distance calculation
traits_subset <- traits[,straits_sub_keep]
traits_sub_dist <- gowdis(traits_subset)
traits_sub_clust <- hclust(traits_sub_dist,"average")
traits_sub_tree <- as.phylo(traits_sub_clust)

#Mean pairwise differences
traits_sub_sesmpd <- ses.mpd(presabs_lake, traits_sub_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

#Mean nearest taxon distance
traits_sub_sesmntd <- ses.mntd(presabs_lake, traits_sub_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

#Faith's PD
traits_sub_sespd <- ses.pd(presabs_lake, traits_sub_tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)


# Trait distance calculation
straits_subset <- straits[,straits_sub_keep]
straits_sub_dist <- gowdis(straits_subset)
straits_sub_clust <- hclust(straits_sub_dist,"average")
straits_sub_tree <- as.phylo(straits_sub_clust)

#Mean pairwise differences
straits_sub_sesmpd <- ses.mpd(surveyed_sites_lake, straits_sub_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

#Mean nearest taxon distance
straits_sub_sesmntd <- ses.mntd(surveyed_sites_lake, straits_sub_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

#Faith's PD
straits_sub_sespd <- ses.pd(surveyed_sites_lake, straits_sub_tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)
```

#### Fourth corner approach

``` r
library(mvabund)
library(lattice)
rows_with_na <- which(rowSums(is.na(straits_vif)) > 0)
rows_with_na

ft1=traitglm(surveyed_sites_lake[surveyed_sites_env,], env[surveyed_sites_env,environment], straits_vif[,c(1:13)])
ft1$fourth

a = max(abs(ft1$fourth.corner))
colort = colorRampPalette(c("blue","purple","red")) 
plot.4th = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list(x = list(rot = 45)))
print(plot.4th)

ft <- fourthcorner(env[surveyed_sites_env,environment], surveyed_sites_lake[surveyed_sites_env,], straits_vif[,c(1:5,7:10,12:13)], tr01 = T,modeltype=6)

straits_comp <- straits_vif[complete.cases(straits_vif),]
ssl_comp <- surveyed_sites_lake[,complete.cases(straits_vif)]

ft <- fourthcorner(env[surveyed_sites_env,environment], ssl_comp[surveyed_sites_env,], straits_comp, tr01 = T,modeltype=6)

ft <- fourthcorner(env[surveyed_sites_env,environment], surveyed_sites_lake[surveyed_sites_env,], straits_vif[,c(1:5,7:10,12:13)], tr01 = T,modeltype=6)

straits_comp <- straits_vif[complete.cases(straits_vif),]
ssl_comp <- surveyed_sites_lake[,complete.cases(straits_vif)]

ft <- fourthcorner(env[surveyed_sites_env,environment], ssl_comp[surveyed_sites_env,], straits_comp, tr01 = T,modeltype=6)
```
