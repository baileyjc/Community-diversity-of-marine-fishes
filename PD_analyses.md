Phylogenetic diversity analyses
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

    ##   |                  |          |   0%  |                  |          |   3%                                                                          |                  |.         |   7% [Bringing everything together load modifying files packages]             |                  |.         |  10%                                                                          |                  |.         |  13% [Bringing everything together load in modifying files]                   |                  |..        |  17%                                                                          |                  |..        |  20% [Check species names across files]                                       |                  |..        |  23%                                                                          |                  |...       |  27% [Modify environment data]                                                |                  |...       |  30%                                                                          |                  |...       |  33% [Modify incidence matrices]                                              |                  |....      |  37%                                                                          |                  |....      |  40% [Modify phylogeny]                                                       |                  |....      |  43%                                                                          |                  |.....     |  47% [Modify trait data]                                                      |                  |.....     |  50%                                                                          |                  |.....     |  53% [Modify stratification data frames]                                      |                  |......    |  57%                                                                          |                  |......    |  60% [Modify stratification trait data]                                       |                  |......    |  63%                                                                          |                  |.......   |  67% [Stratification trait data tests]                                        |                  |.......   |  70%                                                                          |                  |.......   |  73% [Modify site trait data frames]                                          |                  |........  |  77%                                                                          |                  |........  |  80% [Site trait data tests]                                                  |                  |........  |  83%                                                                          |                  |......... |  87% [unnamed-chunk-5]                                                        |                  |......... |  90%                                                                          |                  |......... |  93% [unnamed-chunk-6]                                                        |                  |..........|  97%                                                                          |                  |..........| 100% [Bringing everything together load out modified files and session info]

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
library(viridis)
```

    ## Loading required package: viridisLite

    ## 
    ## Attaching package: 'viridis'

    ## The following object is masked from 'package:maps':
    ## 
    ##     unemp

``` r
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
# Site vectors
surveyed_sites <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LCN", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOM", "OOO", "OTM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_wo_LCN <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOM", "OOO", "OTM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_wo_TLN_HLM <- c("BCM", "CLM", "FLK", "GLK", "HLO", "IBK", "LCN", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOM", "OOO", "OTM", "RCA", "SLN", "ULN")

surveyed_sites_env <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OTM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_geo <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LCN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOM", "OOO", "OTM", "RCA", "SLN", "TLN", "ULN")

ocean_mixed_sites <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA", "FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")

ocean_mixed_sites_env <- c("IBK", "NCN", "RCA", "FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")

ocean_mixed_sites_geo <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA", "FLK", "HLO", "MLN", "NLN", "NLU", "OLO", "ULN")

ocean_stratified_sites <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA", "BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")

ocean_stratified_sites_env <- c("IBK", "NCN", "RCA", "BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")

mixed_stratified_lakes <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "LLN", "MLN", "NLK", "NLN", "NLU", "OLO", "OTM", "SLN", "TLN", "ULN")

mixed_stratified_lakes_geo <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "MLN", "NLK", "NLN", "NLU", "OLO", "OTM", "SLN", "TLN", "ULN")

ocean_sites <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA")

ocean_sites_env <- c("IBK", "NCN", "RCA")

mixed_lakes <- c("FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")

mixed_lakes_geo <- c("FLK", "HLO", "MLN", "NLN", "NLU", "OLO", "ULN")

stratified_lakes <- c("BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")

# Define your custom colors
custom_colors <- c("Reference" = "black", "Ocean" = "#EE6363", "Mixed" = "#87CEFA", "Stratified" = "#6E8B3D")

env_cont <- "orange"

trait_cont <- "#68228B"

trait_cat <- "#DEB887"
```

### Edit input files

``` r
env <- env[order(env$X),]
presabs_lake <- presabs_lake[order(row.names(presabs_lake)),]
surveyed_sites_lake <- surveyed_sites_lake[order(row.names(surveyed_sites_lake)),]

stratification_group_ref <- env[,19]
stratification_group <- env[surveyed_sites,19]
```

# Phylogenetic Diversity

## PD alpha Diversity

### PD alpha mntd, mpd, and pd calculations of traits and files to be saved

``` r
## Reference
# Dispersion
tree_disp <- dispersion(comm = presabs_lake, tree = tree, abund = F)
write.csv(tree_disp, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_disp.csv")

#Mean pairwise distances
tree_sesmpd <- ses.mpd(presabs_lake, cophenetic(tree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)
write.csv(tree_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_sesmpd.csv")

#Mean nearest taxon distances
tree_sesmntd <- ses.mntd(presabs_lake, cophenetic(tree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)
write.csv(tree_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_sesmntd.csv")

#Faith's PD
tree_sespd <- ses.pd(presabs_lake, tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)
write.csv(tree_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_sespd.csv")


## Surveyed sites
# Dispersion
stree_disp <- dispersion(comm = surveyed_sites_lake, tree = stree, abund = F)
write.csv(stree_disp, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_disp.csv")

#Mean pairwise distances
stree_sesmpd <- ses.mpd(surveyed_sites_lake, cophenetic(stree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)
write.csv(stree_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_sesmpd.csv")

#Mean nearest taxon distances
stree_sesmntd <- ses.mntd(surveyed_sites_lake, cophenetic(stree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)
write.csv(stree_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_sesmntd.csv")

#Faith's PD
stree_sespd <- ses.pd(surveyed_sites_lake, stree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)
write.csv(stree_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_sespd.csv")
```

### PD alpha files mntd, mpd, and pd to read into R

- Read in Phylogenetic Diversity files and combine with env data

``` r
## Reference
# Dispersion 
tree_disp <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_disp.csv")
tree_disp_env <- merge(tree_disp, env, by = "X", sort = F)
tree_disp_env$measure <- "tree_disp"
tree_disp_env$Stratification <- factor(tree_disp_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
row.names(tree_disp_env) <- tree_disp_env$X

# Mean pairwise distances
tree_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_sesmpd.csv")
tree_sesmpd_env <- merge(tree_sesmpd, env, by = "X", sort = F)
tree_sesmpd_env$measure <- "tree_sesmpd"
tree_sesmpd_env$Stratification <- factor(tree_sesmpd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
row.names(tree_sesmpd_env) <- tree_sesmpd_env$X

# Mean nearest taxon distances
tree_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_sesmntd.csv")
tree_sesmntd_env <- merge(tree_sesmntd, env, by = "X", sort = F)
tree_sesmntd_env$measure <- "tree_sesmntd"
tree_sesmntd_env$Stratification <- factor(tree_sesmntd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
row.names(tree_sesmntd_env) <- tree_sesmntd_env$X

#Faith's PD
tree_sespd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_sespd.csv")
tree_sespd_env <- merge(tree_sespd, env, by = "X", sort = F)
tree_sespd_env$measure <- "tree_sespd"
tree_sespd_env$Stratification <- factor(tree_sespd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
stree_sespd_env <- merge(tree_sespd_env, tree_disp, by = "X", sort = F)
row.names(tree_sespd_env) <- tree_sespd_env$X


## Surveyed sites
# Dispersion 
stree_disp <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_disp.csv")
stree_disp_env <- merge(stree_disp, env, by = "X", sort = F)
stree_disp_env$measure <- "stree_disp"
stree_disp_env$Stratification <- factor(stree_disp_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))
row.names(stree_disp_env) <- stree_disp_env$X

# Mean pairwise distances
stree_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_sesmpd.csv")
stree_sesmpd_env <- merge(stree_sesmpd, env, by = "X", sort = F)
stree_sesmpd_env$measure <- "stree_sesmpd"
stree_sesmpd_env$Stratification <- factor(stree_sesmpd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))
row.names(stree_sesmpd_env) <- stree_sesmpd_env$X

# Mean nearest taxon distances
stree_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_sesmntd.csv")
stree_sesmntd_env <- merge(stree_sesmntd, env, by = "X", sort = F)
stree_sesmntd_env$measure <- "stree_sesmntd"
stree_sesmntd_env$Stratification <- factor(stree_sesmntd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))
row.names(stree_sesmntd_env) <- stree_sesmntd_env$X

#Faith's PD
stree_sespd <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_sespd.csv")
stree_sespd_env <- merge(stree_sespd, env, by = "X", sort = F)
stree_sespd_env$measure <- "stree_sespd"
stree_sespd_env$Stratification <- factor(stree_sespd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))
stree_sespd_env <- merge(stree_sespd_env, stree_disp, by = "X", sort = F)
row.names(stree_sespd_env) <- stree_sespd_env$X
```

### PD alpha outliers for each stratification type based on z-scores

``` r
outlier_PD_alpha <- stree_sespd_env %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(pd.obs.z, 0.25),
    Q3 = quantile(pd.obs.z, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = pd.obs.z < lower_bound | pd.obs.z > upper_bound
  )

outlier_PD_alpha <- as.data.frame(outlier_PD_alpha)
row.names(outlier_PD_alpha) <- outlier_PD_alpha$X

# Create the plot
(outlier_PD_alpha_plot <- ggplot(outlier_PD_alpha, aes(x = Stratification, y = pd.obs.z, fill = Stratification)) +
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
  labs(y = "PD alpha z-scores", x = "Site type", color = "Outlier", fill = "Site type", tag = "b")) +
  guides(fill = "none")
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/outlier_PD_alpha.jpg", outlier_PD_alpha_plot, width = 6.26, height = 6, units = "in")
```

### Identify VIF of env and geo variables

``` r
# Determine which env variables have variance
environment <- c("temperature_median", "salinity_median", "oxygen_median", "pH_median")
nzv <- nearZeroVar(env[,environment])
nzv
```

    ## integer(0)

``` r
# All env variables have variance

mod <-  lm(pd.obs ~ temperature_median + salinity_median + oxygen_median + pH_median, stree_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = stree_sespd_env)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -960.6 -344.4    5.4  266.2 1158.6 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -23046.44    7912.77  -2.913   0.0114 *
    ## temperature_median   -226.09     151.83  -1.489   0.1586  
    ## salinity_median        29.70      77.33   0.384   0.7067  
    ## oxygen_median         187.18     209.39   0.894   0.3865  
    ## pH_median            3973.56    1368.75   2.903   0.0116 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 638.1 on 14 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.8326, Adjusted R-squared:  0.7847 
    ## F-statistic:  17.4 on 4 and 14 DF,  p-value: 2.52e-05

``` r
anova(mod)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##                    Df   Sum Sq  Mean Sq F value    Pr(>F)    
    ## temperature_median  1  3431415  3431415  8.4262  0.011578 *  
    ## salinity_median     1 15790866 15790866 38.7761 2.211e-05 ***
    ## oxygen_median       1  5692040  5692040 13.9774  0.002202 ** 
    ## pH_median           1  3432036  3432036  8.4277  0.011572 *  
    ## Residuals          14  5701242   407232                      
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
mod <-  lm(pd.obs ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1060.5  -388.2  -180.5   295.4  1777.7 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        -5063.68    6020.36  -0.841  0.41350   
    ## temperature_median   -43.83     169.04  -0.259  0.79895   
    ## salinity_median      203.57      59.81   3.403  0.00393 **
    ## oxygen_median        588.26     192.40   3.057  0.00798 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 780.3 on 15 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.7317, Adjusted R-squared:  0.6781 
    ## F-statistic: 13.64 on 3 and 15 DF,  p-value: 0.0001466

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

# Determine which geo variables have variance
geography <- c("volume_m3_w_chemocline", "volume_m3", "surface_area_m2", "distance_to_ocean_min_m", "distance_to_ocean_mean_m", "distance_to_ocean_median_m", "tidal_lag_time_minutes", "tidal_efficiency", "perimeter_fromSat", "max_depth", "logArea")
nzv <- nearZeroVar(env[,geography])
nzv
```

    ## integer(0)

``` r
# All geo variables have variance

mod <-  lm(pd.obs ~ volume_m3_w_chemocline + volume_m3 + surface_area_m2 + distance_to_ocean_min_m + distance_to_ocean_mean_m + distance_to_ocean_median_m + tidal_lag_time_minutes + tidal_efficiency + perimeter_fromSat + max_depth + logArea, stree_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ volume_m3_w_chemocline + volume_m3 + surface_area_m2 + 
    ##     distance_to_ocean_min_m + distance_to_ocean_mean_m + distance_to_ocean_median_m + 
    ##     tidal_lag_time_minutes + tidal_efficiency + perimeter_fromSat + 
    ##     max_depth + logArea, data = stree_sespd_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1236.90  -373.04    26.01   421.08  1540.38 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                -3.481e+03  5.821e+03  -0.598    0.565  
    ## volume_m3_w_chemocline      1.813e-03  2.127e-03   0.852    0.416  
    ## volume_m3                  -1.536e-03  2.204e-03  -0.697    0.503  
    ## surface_area_m2            -3.984e-03  2.120e-03  -1.879    0.093 .
    ## distance_to_ocean_min_m    -3.792e+00  1.163e+01  -0.326    0.752  
    ## distance_to_ocean_mean_m    2.534e+01  4.939e+01   0.513    0.620  
    ## distance_to_ocean_median_m -2.641e+01  4.540e+01  -0.582    0.575  
    ## tidal_lag_time_minutes     -8.253e+00  1.363e+01  -0.606    0.560  
    ## tidal_efficiency           -5.195e+02  3.631e+03  -0.143    0.889  
    ## perimeter_fromSat          -6.617e-01  8.975e-01  -0.737    0.480  
    ## max_depth                  -3.843e+01  5.981e+01  -0.643    0.536  
    ## logArea                     8.836e+02  5.809e+02   1.521    0.163  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 966.8 on 9 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.7095, Adjusted R-squared:  0.3544 
    ## F-statistic: 1.998 on 11 and 9 DF,  p-value: 0.1541

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
mod <-  lm(pd.obs ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ distance_to_ocean_min_m + max_depth + logArea, 
    ##     data = stree_sespd_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1739.8  -528.6    69.2   274.2  2683.2 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             1894.996   1521.690   1.245    0.229  
    ## distance_to_ocean_min_m   -9.726      3.431  -2.835    0.011 *
    ## max_depth                 -6.612     25.223  -0.262    0.796  
    ## logArea                  131.270    151.202   0.868    0.397  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1093 on 18 degrees of freedom
    ## Multiple R-squared:  0.4137, Adjusted R-squared:  0.3159 
    ## F-statistic: 4.233 on 3 and 18 DF,  p-value: 0.01981

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
```

### PD alpha PRic and z-score & stratification

``` r
### Stratification
# Run ANOVA on the original data
anova_PRic_result <- aov(pd.obs ~ Stratification, data = stree_sespd_env[surveyed_sites,])
summary(anova_PRic_result)
```

    ##                Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Stratification  2 23597133 11798567   17.19 5.48e-05 ***
    ## Residuals      19 13044225   686538                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Perform Tukey's HSD test
tukey_PRic_result <- TukeyHSD(anova_PRic_result)
print(tukey_PRic_result)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = pd.obs ~ Stratification, data = stree_sespd_env[surveyed_sites, ])
    ## 
    ## $Stratification
    ##                        diff        lwr        upr     p adj
    ## Mixed-Ocean        218.8183  -917.9877  1355.6244 0.8773441
    ## Stratified-Ocean -2020.3990 -3157.2051  -883.5929 0.0006627
    ## Stratified-Mixed -2239.2173 -3291.6953 -1186.7394 0.0000923

``` r
# Sites
N <- length(anova_PRic_result$residuals)
# Categories
k <- length(anova_PRic_result$coefficients)
# Sum of squares
SS <- sum(resid(anova_PRic_result)^2) # from anova, the residual SS
# Mean squares
MS <- SS/(N-k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

# Get differences
Stratification <- tukey_PRic_result$Stratification

# q-values
(OM_q <- (abs(Stratification[1,1]))/USE)
```

    ## [1] 0.6915491

``` r
(MS_q <- (abs(Stratification[3,1]))/BSE)
```

    ## [1] 7.643793

``` r
(SO_q <- (abs(Stratification[2,1]))/USE)
```

    ## [1] 6.385228

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx

# Run ANOVA on the original data
anova_zscore_result <- aov(pd.obs.z ~ Stratification, data = stree_sespd_env[surveyed_sites,])
summary(anova_zscore_result)
```

    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## Stratification  2  1.838   0.919   0.888  0.428
    ## Residuals      19 19.668   1.035

``` r
# Perform Tukey's HSD test
tukey_zscore_result <- TukeyHSD(anova_zscore_result)
print(tukey_zscore_result)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = pd.obs.z ~ Stratification, data = stree_sespd_env[surveyed_sites, ])
    ## 
    ## $Stratification
    ##                          diff        lwr      upr     p adj
    ## Mixed-Ocean      0.6488785535 -0.7470411 2.044798 0.4784457
    ## Stratified-Ocean 0.6491655128 -0.7467542 2.045085 0.4781464
    ## Stratified-Mixed 0.0002869593 -1.2920836 1.292657 0.9999998

``` r
# Sites
N <- length(anova_zscore_result$residuals)
# Categories
k <- length(anova_zscore_result$coefficients)
# Sum of squares
SS <- sum(resid(anova_zscore_result)^2) # from anova, the residual SS
# Mean squares
MS <- SS/(N-k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

# Get differences
Stratification <- tukey_zscore_result$Stratification

# q-values
(OM_q <- (abs(Stratification[1,1]))/USE)
```

    ## [1] 1.670047

``` r
(MS_q <- (abs(Stratification[3,1]))/BSE)
```

    ## [1] 0.0007977355

``` r
(SO_q <- (abs(Stratification[2,1]))/USE)
```

    ## [1] 1.670785

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx

# Run ANOVA on the original data
anova_PDisp_result <- aov(Dispersion ~ Stratification, data = stree_sespd_env[surveyed_sites,])
summary(anova_PDisp_result)
```

    ##                Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Stratification  2 0.000707 0.0003537   0.703  0.507
    ## Residuals      19 0.009556 0.0005030

``` r
# Perform Tukey's HSD test
tukey_PDisp_result <- TukeyHSD(anova_PDisp_result)
print(tukey_PDisp_result)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = Dispersion ~ Stratification, data = stree_sespd_env[surveyed_sites, ])
    ## 
    ## $Stratification
    ##                          diff         lwr        upr     p adj
    ## Mixed-Ocean       0.007071665 -0.02369775 0.03784108 0.8302990
    ## Stratified-Ocean -0.006219734 -0.03698915 0.02454968 0.8657120
    ## Stratified-Mixed -0.013291399 -0.04177834 0.01519555 0.4759158

``` r
# Sites
N <- length(anova_PDisp_result$residuals)
# Categories
k <- length(anova_PDisp_result$coefficients)
# Sum of squares
SS <- sum(resid(anova_PDisp_result)^2) # from anova, the residual SS
# Mean squares
MS <- SS/(N-k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

# Get differences
Stratification <- tukey_PDisp_result$Stratification

# q-values
(OM_q <- (abs(Stratification[1,1]))/USE)
```

    ## [1] 0.825711

``` r
(MS_q <- (abs(Stratification[3,1]))/BSE)
```

    ## [1] 1.676295

``` r
(SO_q <- (abs(Stratification[2,1]))/USE)
```

    ## [1] 0.7262367

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx
```

### PD alpha PRic and z-score & env and geo linear models

``` r
### Environmental
# Surveyed sites
PD_PRic_env_model <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[surveyed_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.34125 -0.17839 -0.09832  0.17872  0.51241 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        2.210239   2.191942   1.008  0.32928    
    ## temperature_median 0.009972   0.061544   0.162  0.87345    
    ## salinity_median    0.130178   0.021778   5.978 2.53e-05 ***
    ## oxygen_median      0.244896   0.070051   3.496  0.00325 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2841 on 15 degrees of freedom
    ## Multiple R-squared:  0.8537, Adjusted R-squared:  0.8245 
    ## F-statistic: 29.18 on 3 and 15 DF,  p-value: 1.658e-06

``` r
anova(PD_PRic_env_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                    Df Sum Sq Mean Sq F value    Pr(>F)    
    ## temperature_median  1 0.7740  0.7740  9.5888  0.007369 ** 
    ## salinity_median     1 5.3056  5.3056 65.7336 7.301e-07 ***
    ## oxygen_median       1 0.9865  0.9865 12.2219  0.003250 ** 
    ## Residuals          15 1.2107  0.0807                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##       1.0000000000       1.0000000000       0.0001013439       0.0130007650

``` r
PD_z_env_model <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[surveyed_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.0491 -0.6604  0.2329  0.6865  1.6970 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.05264    7.98903  -0.007    0.995
    ## temperature_median  0.05069    0.22431   0.226    0.824
    ## salinity_median    -0.05571    0.07937  -0.702    0.493
    ## oxygen_median      -0.20248    0.25532  -0.793    0.440
    ## 
    ## Residual standard error: 1.035 on 15 degrees of freedom
    ## Multiple R-squared:  0.1376, Adjusted R-squared:  -0.03483 
    ## F-statistic: 0.798 on 3 and 15 DF,  p-value: 0.514

``` r
anova(PD_z_env_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                    Df  Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1  0.5408 0.54076  0.5043 0.4885
    ## salinity_median     1  1.3519 1.35190  1.2609 0.2791
    ## oxygen_median       1  0.6743 0.67435  0.6289 0.4401
    ## Residuals          15 16.0831 1.07221

``` r
p_values <- summary(PD_z_env_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_model <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[surveyed_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.032704 -0.010212 -0.000927  0.012572  0.048834 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -0.2235120  0.1681181  -1.329   0.2036  
    ## temperature_median  0.0097192  0.0047203   2.059   0.0573 .
    ## salinity_median     0.0030130  0.0016703   1.804   0.0914 .
    ## oxygen_median      -0.0002191  0.0053728  -0.041   0.9680  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02179 on 15 degrees of freedom
    ## Multiple R-squared:  0.2845, Adjusted R-squared:  0.1414 
    ## F-statistic: 1.988 on 3 and 15 DF,  p-value: 0.1591

``` r
anova(PD_PDisp_env_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                    Df    Sum Sq    Mean Sq F value  Pr(>F)  
    ## temperature_median  1 0.0010125 0.00101250  2.1324 0.16484  
    ## salinity_median     1 0.0018190 0.00181899  3.8310 0.06919 .
    ## oxygen_median       1 0.0000008 0.00000079  0.0017 0.96801  
    ## Residuals          15 0.0071221 0.00047481                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_env_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          0.8142048          0.2291871          0.3655132          1.0000000

``` r
# Mixed and stratified lakes
PD_PRic_env_MS_model <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_env_MS_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[mixed_stratified_lakes, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.40218 -0.20251 -0.06598  0.23118  0.45102 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         1.72817    2.52835   0.684 0.507264    
    ## temperature_median  0.01795    0.06991   0.257 0.801749    
    ## salinity_median     0.13399    0.02390   5.606 0.000115 ***
    ## oxygen_median       0.28179    0.08638   3.262 0.006801 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3062 on 12 degrees of freedom
    ## Multiple R-squared:  0.8431, Adjusted R-squared:  0.8039 
    ## F-statistic: 21.49 on 3 and 12 DF,  p-value: 4.073e-05

``` r
anova(PD_PRic_env_MS_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                    Df Sum Sq Mean Sq F value    Pr(>F)    
    ## temperature_median  1 0.8524  0.8524  9.0901  0.010763 *  
    ## salinity_median     1 4.1956  4.1956 44.7418 2.231e-05 ***
    ## oxygen_median       1 0.9979  0.9979 10.6420  0.006801 ** 
    ## Residuals          12 1.1253  0.0938                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_MS_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##       1.0000000000       1.0000000000       0.0004600793       0.0272028786

``` r
PD_z_env_MS_model <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_z_env_MS_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[mixed_stratified_lakes, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.40930 -0.48650  0.04376  0.57984  1.20044 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -7.01848    6.96162  -1.008    0.333
    ## temperature_median  0.20879    0.19249   1.085    0.299
    ## salinity_median    -0.02443    0.06581  -0.371    0.717
    ## oxygen_median       0.13066    0.23784   0.549    0.593
    ## 
    ## Residual standard error: 0.8432 on 12 degrees of freedom
    ## Multiple R-squared:  0.1298, Adjusted R-squared:  -0.08774 
    ## F-statistic: 0.5967 on 3 and 12 DF,  p-value: 0.6292

``` r
anova(PD_z_env_MS_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                    Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1 1.0185 1.01846  1.4326 0.2545
    ## salinity_median     1 0.0396 0.03961  0.0557 0.8174
    ## oxygen_median       1 0.2146 0.21456  0.3018 0.5928
    ## Residuals          12 8.5312 0.71093

``` r
p_values <- summary(PD_z_env_MS_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_MS_model <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_env_MS_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[mixed_stratified_lakes, 
    ##     ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.026782 -0.011687 -0.000974  0.010088  0.046961 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -0.325673   0.182766  -1.782   0.1001  
    ## temperature_median  0.012123   0.005054   2.399   0.0336 *
    ## salinity_median     0.003451   0.001728   1.997   0.0690 .
    ## oxygen_median       0.004063   0.006244   0.651   0.5275  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02214 on 12 degrees of freedom
    ## Multiple R-squared:  0.3997, Adjusted R-squared:  0.2497 
    ## F-statistic: 2.664 on 3 and 12 DF,  p-value: 0.09534

``` r
anova(PD_PDisp_env_MS_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                    Df    Sum Sq    Mean Sq F value  Pr(>F)  
    ## temperature_median  1 0.0012390 0.00123905  2.5286 0.13778  
    ## salinity_median     1 0.0024691 0.00246905  5.0388 0.04442 *
    ## oxygen_median       1 0.0002075 0.00020749  0.4234 0.52749  
    ## Residuals          12 0.0058801 0.00049001                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_env_MS_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          0.4002584          0.1343561          0.2759753          1.0000000

``` r
# Ocean sites and mixed lakes
PD_PRic_env_OM_model <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PRic_env_OM_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[ocean_mixed_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.34221 -0.18181  0.00482  0.10561  0.49457 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         3.67640   12.32519   0.298    0.774
    ## temperature_median -0.07279    0.17016  -0.428    0.682
    ## salinity_median     0.16910    0.31857   0.531    0.612
    ## oxygen_median       0.20225    0.25697   0.787    0.457
    ## 
    ## Residual standard error: 0.3262 on 7 degrees of freedom
    ## Multiple R-squared:  0.2849, Adjusted R-squared:  -0.02157 
    ## F-statistic: 0.9296 on 3 and 7 DF,  p-value: 0.4751

``` r
anova(PD_PRic_env_OM_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                    Df  Sum Sq  Mean Sq F value Pr(>F)
    ## temperature_median  1 0.00704 0.007044  0.0662 0.8044
    ## salinity_median     1 0.22387 0.223866  2.1033 0.1903
    ## oxygen_median       1 0.06593 0.065929  0.6194 0.4571
    ## Residuals           7 0.74507 0.106438

``` r
p_values <- summary(PD_PRic_env_OM_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_z_env_OM_model <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_z_env_OM_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[ocean_mixed_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.8318 -0.4833 -0.2593  0.5477  1.0164 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        31.53378   30.82399   1.023    0.340
    ## temperature_median -0.05226    0.42556  -0.123    0.906
    ## salinity_median    -0.77195    0.79671  -0.969    0.365
    ## oxygen_median      -1.15354    0.64267  -1.795    0.116
    ## 
    ## Residual standard error: 0.8159 on 7 degrees of freedom
    ## Multiple R-squared:  0.641,  Adjusted R-squared:  0.4872 
    ## F-statistic: 4.167 on 3 and 7 DF,  p-value: 0.05472

``` r
anova(PD_z_env_OM_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                    Df Sum Sq Mean Sq F value  Pr(>F)  
    ## temperature_median  1 0.2716  0.2716  0.4080 0.54331  
    ## salinity_median     1 5.9057  5.9057  8.8712 0.02056 *
    ## oxygen_median       1 2.1448  2.1448  3.2218 0.11574  
    ## Residuals           7 4.6600  0.6657                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_env_OM_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          1.0000000          1.0000000          0.4629584

``` r
PD_PDisp_env_OM_model <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PDisp_env_OM_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[ocean_mixed_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.009814 -0.004982 -0.001190  0.001979  0.019442 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         2.457e-01  3.643e-01   0.674    0.522
    ## temperature_median -5.727e-05  5.030e-03  -0.011    0.991
    ## salinity_median    -8.417e-04  9.417e-03  -0.089    0.931
    ## oxygen_median      -9.043e-03  7.596e-03  -1.190    0.273
    ## 
    ## Residual standard error: 0.009644 on 7 degrees of freedom
    ## Multiple R-squared:  0.3056, Adjusted R-squared:  0.007994 
    ## F-statistic: 1.027 on 3 and 7 DF,  p-value: 0.4369

``` r
anova(PD_PDisp_env_OM_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                    Df     Sum Sq    Mean Sq F value Pr(>F)
    ## temperature_median  1 0.00001581 1.5809e-05  0.1700 0.6925
    ## salinity_median     1 0.00013890 1.3890e-04  1.4934 0.2612
    ## oxygen_median       1 0.00013182 1.3181e-04  1.4172 0.2727
    ## Residuals           7 0.00065107 9.3010e-05

``` r
p_values <- summary(PD_PDisp_env_OM_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
# Stratified lakes and ocean sites
PD_PRic_env_SO_model <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PRic_env_SO_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[ocean_stratified_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.33792 -0.07625  0.00662  0.06006  0.33720 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         1.68269    1.85239   0.908 0.393871    
    ## temperature_median  0.03913    0.05299   0.738 0.484238    
    ## salinity_median     0.11951    0.01902   6.283 0.000411 ***
    ## oxygen_median       0.21178    0.06018   3.519 0.009743 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.227 on 7 degrees of freedom
    ## Multiple R-squared:  0.9111, Adjusted R-squared:  0.873 
    ## F-statistic: 23.92 on 3 and 7 DF,  p-value: 0.0004701

``` r
anova(PD_PRic_env_SO_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                    Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## temperature_median  1 0.06831 0.06831  1.3253 0.2874402    
    ## salinity_median     1 2.99225 2.99225 58.0564 0.0001243 ***
    ## oxygen_median       1 0.63818 0.63818 12.3822 0.0097429 ** 
    ## Residuals           7 0.36078 0.05154                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_SO_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##        1.000000000        1.000000000        0.001642753        0.038971567

``` r
PD_z_env_SO_model <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_z_env_SO_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[ocean_stratified_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.2602 -0.6625 -0.1193  0.7336  1.1685 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         1.47357    8.55414   0.172    0.868
    ## temperature_median  0.06716    0.24472   0.274    0.792
    ## salinity_median    -0.13008    0.08783  -1.481    0.182
    ## oxygen_median      -0.26149    0.27793  -0.941    0.378
    ## 
    ## Residual standard error: 1.048 on 7 degrees of freedom
    ## Multiple R-squared:  0.4059, Adjusted R-squared:  0.1514 
    ## F-statistic: 1.594 on 3 and 7 DF,  p-value: 0.2747

``` r
anova(PD_z_env_SO_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                    Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1 0.6060  0.6060  0.5514 0.4819
    ## salinity_median     1 3.6786  3.6786  3.3469 0.1100
    ## oxygen_median       1 0.9729  0.9729  0.8852 0.3781
    ## Residuals           7 7.6937  1.0991

``` r
p_values <- summary(PD_z_env_SO_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          1.0000000          0.7285996          1.0000000

``` r
PD_PDisp_env_SO_model <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PDisp_env_SO_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[ocean_stratified_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.027570 -0.022235  0.006409  0.012965  0.050260 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.249416   0.239944  -1.039    0.333
    ## temperature_median  0.011174   0.006864   1.628    0.148
    ## salinity_median     0.002342   0.002464   0.951    0.373
    ## oxygen_median      -0.001076   0.007796  -0.138    0.894
    ## 
    ## Residual standard error: 0.02941 on 7 degrees of freedom
    ## Multiple R-squared:  0.3004, Adjusted R-squared:  0.0005571 
    ## F-statistic: 1.002 on 3 and 7 DF,  p-value: 0.4464

``` r
anova(PD_PDisp_env_SO_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                    Df    Sum Sq    Mean Sq F value Pr(>F)
    ## temperature_median  1 0.0017999 0.00179990  2.0814 0.1923
    ## salinity_median     1 0.0007828 0.00078275  0.9052 0.3731
    ## oxygen_median       1 0.0000165 0.00001647  0.0190 0.8941
    ## Residuals           7 0.0060534 0.00086477

``` r
p_values <- summary(PD_PDisp_env_SO_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          0.5902869          1.0000000          1.0000000

``` r
# Mixed lakes
PD_PRic_env_M_model <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_lakes,])
summary(PD_PRic_env_M_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##      FLK      HLO      LLN      MLN      NLN      NLU      OLO      ULN 
    ##  0.04358 -0.23393  0.33491 -0.15413  0.48318 -0.31577 -0.28003  0.12220 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.72550   18.55796  -0.039    0.971
    ## temperature_median -0.06134    0.26779  -0.229    0.830
    ## salinity_median     0.28019    0.43143   0.649    0.551
    ## oxygen_median       0.29681    0.34898   0.850    0.443
    ## 
    ## Residual standard error: 0.3934 on 4 degrees of freedom
    ## Multiple R-squared:  0.3807, Adjusted R-squared:  -0.08376 
    ## F-statistic: 0.8197 on 3 and 4 DF,  p-value: 0.5469

``` r
anova(PD_PRic_env_M_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                    Df  Sum Sq  Mean Sq F value Pr(>F)
    ## temperature_median  1 0.07422 0.074216  0.4795 0.5267
    ## salinity_median     1 0.19440 0.194402  1.2561 0.3251
    ## oxygen_median       1 0.11195 0.111950  0.7233 0.4430
    ## Residuals           4 0.61907 0.154766

``` r
p_values <- summary(PD_PRic_env_M_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_z_env_M_model <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_lakes,])
summary(PD_z_env_M_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##     FLK     HLO     LLN     MLN     NLN     NLU     OLO     ULN 
    ## -0.6311 -0.9440 -0.1008  0.4796  0.7484  0.1079  0.1433  0.1967 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         18.8693    34.6682   0.544    0.615
    ## temperature_median   0.3083     0.5003   0.616    0.571
    ## salinity_median     -0.7909     0.8060  -0.981    0.382
    ## oxygen_median       -0.6458     0.6519  -0.991    0.378
    ## 
    ## Residual standard error: 0.7349 on 4 degrees of freedom
    ## Multiple R-squared:  0.5607, Adjusted R-squared:  0.2313 
    ## F-statistic: 1.702 on 3 and 4 DF,  p-value: 0.3034

``` r
anova(PD_z_env_M_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                    Df  Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1 0.92351 0.92351  1.7099 0.2611
    ## salinity_median     1 1.30448 1.30448  2.4152 0.1951
    ## oxygen_median       1 0.52992 0.52992  0.9811 0.3780
    ## Residuals           4 2.16043 0.54011

``` r
p_values <- summary(PD_z_env_M_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_M_model <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_lakes,])
summary(PD_PDisp_env_M_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##        FLK        HLO        LLN        MLN        NLN        NLU        OLO 
    ## -8.607e-03 -6.537e-03  8.692e-05 -3.098e-03  1.640e-02  7.300e-04 -6.838e-03 
    ##        ULN 
    ##  7.864e-03 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.016539   0.529784  -0.031    0.977
    ## temperature_median  0.004574   0.007645   0.598    0.582
    ## salinity_median     0.002243   0.012316   0.182    0.864
    ## oxygen_median      -0.004627   0.009963  -0.464    0.666
    ## 
    ## Residual standard error: 0.01123 on 4 degrees of freedom
    ## Multiple R-squared:  0.1322, Adjusted R-squared:  -0.5187 
    ## F-statistic: 0.2031 on 3 and 4 DF,  p-value: 0.8894

``` r
anova(PD_PDisp_env_M_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                    Df     Sum Sq    Mean Sq F value Pr(>F)
    ## temperature_median  1 0.00004959 4.9595e-05  0.3932 0.5646
    ## salinity_median     1 0.00000004 4.2000e-08  0.0003 0.9864
    ## oxygen_median       1 0.00002720 2.7204e-05  0.2157 0.6665
    ## Residuals           4 0.00050452 1.2613e-04

``` r
p_values <- summary(PD_PDisp_env_M_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
# Stratified lakes
PD_PRic_env_S_model <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[stratified_lakes,])
summary(PD_PRic_env_S_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##       BCM       CLM       GLK       HLM       NLK       OTM       SLN       TLN 
    ##  0.222256  0.002541  0.093548  0.144056 -0.223573 -0.278494 -0.109216  0.148881 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         4.17718    2.71151   1.541    0.198
    ## temperature_median  0.02913    0.06006   0.485    0.653
    ## salinity_median     0.06782    0.04211   1.610    0.183
    ## oxygen_median      -0.03224    0.18566  -0.174    0.871
    ## 
    ## Residual standard error: 0.2452 on 4 degrees of freedom
    ## Multiple R-squared:  0.6503, Adjusted R-squared:  0.388 
    ## F-statistic: 2.479 on 3 and 4 DF,  p-value: 0.2006

``` r
anova(PD_PRic_env_S_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                    Df  Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median  1 0.00002 0.00002  0.0003 0.9862  
    ## salinity_median     1 0.44541 0.44541  7.4067 0.0529 .
    ## oxygen_median       1 0.00181 0.00181  0.0302 0.8706  
    ## Residuals           4 0.24054 0.06014                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_S_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          0.7931036          1.0000000          0.7303422          1.0000000

``` r
PD_z_env_S_model <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[stratified_lakes,])
summary(PD_z_env_S_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##     BCM     CLM     GLK     HLM     NLK     OTM     SLN     TLN 
    ##  0.6329 -0.5528  0.2931  0.0925 -0.7839  1.2305 -0.4533 -0.4590 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -9.36668   10.11196  -0.926    0.407
    ## temperature_median  0.18240    0.22397   0.814    0.461
    ## salinity_median     0.03683    0.15704   0.235    0.826
    ## oxygen_median       0.57919    0.69238   0.837    0.450
    ## 
    ## Residual standard error: 0.9145 on 4 degrees of freedom
    ## Multiple R-squared:  0.3152, Adjusted R-squared:  -0.1983 
    ## F-statistic: 0.6138 on 3 and 4 DF,  p-value: 0.6412

``` r
anova(PD_z_env_S_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                    Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1 0.6162 0.61619  0.7368 0.4391
    ## salinity_median     1 0.3387 0.33869  0.4050 0.5591
    ## oxygen_median       1 0.5853 0.58526  0.6998 0.4499
    ## Residuals           4 3.3454 0.83635

``` r
p_values <- summary(PD_z_env_S_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_S_model <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[stratified_lakes,])
summary(PD_PDisp_env_S_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##       BCM       CLM       GLK       HLM       NLK       OTM       SLN       TLN 
    ##  0.018547 -0.026500  0.037634 -0.004772 -0.034031  0.026157 -0.020451  0.003415 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.491926   0.381330  -1.290    0.267
    ## temperature_median  0.013502   0.008446   1.599    0.185
    ## salinity_median     0.006309   0.005922   1.065    0.347
    ## oxygen_median       0.018246   0.026110   0.699    0.523
    ## 
    ## Residual standard error: 0.03449 on 4 degrees of freedom
    ## Multiple R-squared:  0.4408, Adjusted R-squared:  0.0214 
    ## F-statistic: 1.051 on 3 and 4 DF,  p-value: 0.4619

``` r
anova(PD_PDisp_env_S_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                    Df    Sum Sq    Mean Sq F value Pr(>F)
    ## temperature_median  1 0.0023592 0.00235923  1.9836 0.2318
    ## salinity_median     1 0.0008101 0.00081014  0.6812 0.4556
    ## oxygen_median       1 0.0005808 0.00058081  0.4883 0.5232
    ## Residuals           4 0.0047575 0.00118937

``` r
p_values <- summary(PD_PDisp_env_S_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          0.7405587          1.0000000          1.0000000

``` r
### Geographical
# Surveyed sites
PD_PRic_geo_model <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[surveyed_sites_geo,])
summary(PD_PRic_geo_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[surveyed_sites_geo, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.75378 -0.25663  0.06493  0.19185  0.60705 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              7.9434974  0.5821945  13.644 1.38e-10 ***
    ## distance_to_ocean_min_m -0.0062107  0.0013054  -4.758 0.000182 ***
    ## max_depth               -0.0001038  0.0096171  -0.011 0.991513    
    ## logArea                  0.0098393  0.0582911   0.169 0.867949    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4128 on 17 degrees of freedom
    ## Multiple R-squared:  0.6328, Adjusted R-squared:  0.568 
    ## F-statistic: 9.765 on 3 and 17 DF,  p-value: 0.0005634

``` r
anova(PD_PRic_geo_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                         Df Sum Sq Mean Sq F value  Pr(>F)    
    ## distance_to_ocean_min_m  1 4.9847  4.9847 29.2578 4.7e-05 ***
    ## max_depth                1 0.0014  0.0014  0.0083  0.9286    
    ## logArea                  1 0.0049  0.0049  0.0285  0.8679    
    ## Residuals               17 2.8963  0.1704                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##            5.514245e-10            7.296449e-04            1.000000e+00 
    ##                 logArea 
    ##            1.000000e+00

``` r
PD_z_geo_model <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[surveyed_sites_geo,])
summary(PD_z_geo_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[surveyed_sites_geo, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.68687 -0.38327 -0.03268  0.42764  1.56108 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)             -0.899618   1.406672  -0.640    0.531
    ## distance_to_ocean_min_m  0.004960   0.003154   1.573    0.134
    ## max_depth               -0.035436   0.023236  -1.525    0.146
    ## logArea                  0.013862   0.140840   0.098    0.923
    ## 
    ## Residual standard error: 0.9973 on 17 degrees of freedom
    ## Multiple R-squared:  0.2044, Adjusted R-squared:  0.06396 
    ## F-statistic: 1.456 on 3 and 17 DF,  p-value: 0.2619

``` r
anova(PD_z_geo_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1  1.3601 1.36006  1.3674 0.2584
    ## max_depth                1  2.9733 2.97327  2.9894 0.1019
    ## logArea                  1  0.0096 0.00964  0.0097 0.9227
    ## Residuals               17 16.9082 0.99460

``` r
p_values <- summary(PD_z_geo_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               0.5369031               0.5825424 
    ##                 logArea 
    ##               1.0000000

``` r
PD_PDisp_geo_model <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[surveyed_sites_geo,])
summary(PD_PDisp_geo_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[surveyed_sites_geo, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.060645 -0.009683 -0.001792  0.010261  0.041420 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              1.812e-01  3.226e-02   5.618 3.07e-05 ***
    ## distance_to_ocean_min_m -1.148e-04  7.232e-05  -1.588    0.131    
    ## max_depth                3.460e-04  5.328e-04   0.649    0.525    
    ## logArea                 -9.409e-04  3.230e-03  -0.291    0.774    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02287 on 17 degrees of freedom
    ## Multiple R-squared:  0.1336, Adjusted R-squared:  -0.01926 
    ## F-statistic: 0.874 on 3 and 17 DF,  p-value: 0.4739

``` r
anova(PD_PDisp_geo_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                         Df    Sum Sq    Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.0011490 0.00114904  2.1972 0.1566
    ## max_depth                1 0.0001778 0.00017777  0.3399 0.5675
    ## logArea                  1 0.0000444 0.00004439  0.0849 0.7743
    ## Residuals               17 0.0088904 0.00052296

``` r
p_values <- summary(PD_PDisp_geo_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##            0.0001229092            0.5230886216            1.0000000000 
    ##                 logArea 
    ##            1.0000000000

``` r
# Mixed and stratified lakes
PD_PRic_geo_MS_model <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_stratified_lakes_geo,])
summary(PD_PRic_geo_MS_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_stratified_lakes_geo, 
    ##     ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7367 -0.1808  0.1454  0.1928  0.5426 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              6.198472   1.260178   4.919 0.000458 ***
    ## distance_to_ocean_min_m -0.006718   0.001643  -4.090 0.001789 ** 
    ## max_depth               -0.023129   0.016392  -1.411 0.185904    
    ## logArea                  0.238012   0.151228   1.574 0.143823    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4308 on 11 degrees of freedom
    ## Multiple R-squared:  0.6533, Adjusted R-squared:  0.5587 
    ## F-statistic: 6.908 on 3 and 11 DF,  p-value: 0.006996

``` r
anova(PD_PRic_geo_MS_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                         Df Sum Sq Mean Sq F value   Pr(>F)   
    ## distance_to_ocean_min_m  1 3.3741  3.3741  18.177 0.001335 **
    ## max_depth                1 0.0132  0.0132   0.071 0.794867   
    ## logArea                  1 0.4598  0.4598   2.477 0.143823   
    ## Residuals               11 2.0418  0.1856                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_MS_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             0.001831138             0.007155424             0.743617082 
    ##                 logArea 
    ##             0.575290766

``` r
PD_z_geo_MS_model <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_stratified_lakes_geo,])
summary(PD_z_geo_MS_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_stratified_lakes_geo, 
    ##     ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.0495 -0.5087 -0.1472  0.4868  1.1265 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -4.558494   2.406723  -1.894   0.0848 .
    ## distance_to_ocean_min_m  0.001378   0.003137   0.439   0.6690  
    ## max_depth               -0.044080   0.031307  -1.408   0.1868  
    ## logArea                  0.452107   0.288820   1.565   0.1458  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8228 on 11 degrees of freedom
    ## Multiple R-squared:  0.2013, Adjusted R-squared:  -0.01655 
    ## F-statistic: 0.924 on 3 and 11 DF,  p-value: 0.4613

``` r
anova(PD_z_geo_MS_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                         Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.1675 0.16747  0.2474 0.6287
    ## max_depth                1 0.0503 0.05031  0.0743 0.7902
    ## logArea                  1 1.6590 1.65899  2.4503 0.1458
    ## Residuals               11 7.4475 0.67704

``` r
p_values <- summary(PD_z_geo_MS_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.3391891               1.0000000               0.7470347 
    ##                 logArea 
    ##               0.5831770

``` r
PD_PDisp_geo_MS_model <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_stratified_lakes_geo,])
summary(PD_PDisp_geo_MS_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_stratified_lakes_geo, 
    ##     ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.043795 -0.011847 -0.004155  0.011385  0.035405 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              7.284e-02  6.703e-02   1.087   0.3004  
    ## distance_to_ocean_min_m -2.014e-04  8.737e-05  -2.305   0.0417 *
    ## max_depth               -5.016e-04  8.720e-04  -0.575   0.5767  
    ## logArea                  1.292e-02  8.044e-03   1.606   0.1366  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02292 on 11 degrees of freedom
    ## Multiple R-squared:  0.4101, Adjusted R-squared:  0.2492 
    ## F-statistic: 2.549 on 3 and 11 DF,  p-value: 0.1093

``` r
anova(PD_PDisp_geo_MS_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                         Df    Sum Sq    Mean Sq F value  Pr(>F)  
    ## distance_to_ocean_min_m  1 0.0019691 0.00196906  3.7491 0.07894 .
    ## max_depth                1 0.0006931 0.00069308  1.3196 0.27503  
    ## logArea                  1 0.0013545 0.00135452  2.5790 0.13659  
    ## Residuals               11 0.0057773 0.00052521                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_geo_MS_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               0.1666132               1.0000000 
    ##                 logArea 
    ##               0.5463676

``` r
# Ocean sites and mixed lakes
PD_PRic_geo_OM_model <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_mixed_sites_geo,])
summary(PD_PRic_geo_OM_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_mixed_sites_geo, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.40077 -0.14310 -0.06238  0.12398  0.37652 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              7.731682   0.512281  15.093 1.07e-07 ***
    ## distance_to_ocean_min_m -0.000930   0.002511  -0.370    0.720    
    ## max_depth                0.014052   0.008747   1.607    0.143    
    ## logArea                  0.007430   0.045343   0.164    0.873    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.282 on 9 degrees of freedom
    ## Multiple R-squared:  0.3188, Adjusted R-squared:  0.09179 
    ## F-statistic: 1.404 on 3 and 9 DF,  p-value: 0.3039

``` r
anova(PD_PRic_geo_OM_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                         Df  Sum Sq  Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.09150 0.091502  1.1504 0.3114
    ## max_depth                1 0.24144 0.241442  3.0355 0.1154
    ## logArea                  1 0.00214 0.002136  0.0268 0.8735
    ## Residuals                9 0.71585 0.079538

``` r
p_values <- summary(PD_PRic_geo_OM_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##            4.277349e-07            1.000000e+00            5.704839e-01 
    ##                 logArea 
    ##            1.000000e+00

``` r
PD_z_geo_OM_model <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_mixed_sites_geo,])
summary(PD_z_geo_OM_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_mixed_sites_geo, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.74075 -0.53419 -0.21181  0.06012  1.70884 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              0.124303   1.552195   0.080   0.9379  
    ## distance_to_ocean_min_m  0.003932   0.007610   0.517   0.6178  
    ## max_depth               -0.078263   0.026502  -2.953   0.0161 *
    ## logArea                 -0.027456   0.137387  -0.200   0.8460  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8545 on 9 degrees of freedom
    ## Multiple R-squared:  0.5908, Adjusted R-squared:  0.4544 
    ## F-statistic: 4.331 on 3 and 9 DF,  p-value: 0.0378

``` r
anova(PD_z_geo_OM_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                         Df Sum Sq Mean Sq F value  Pr(>F)  
    ## distance_to_ocean_min_m  1 2.1228  2.1228  2.9071 0.12238  
    ## max_depth                1 7.3366  7.3366 10.0471 0.01137 *
    ## logArea                  1 0.0292  0.0292  0.0399 0.84604  
    ## Residuals                9 6.5720  0.7302                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_geo_OM_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               1.0000000               0.0645594 
    ##                 logArea 
    ##               1.0000000

``` r
PD_PDisp_geo_OM_model <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_mixed_sites_geo,])
summary(PD_PDisp_geo_OM_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_mixed_sites_geo, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.008414 -0.005258 -0.004180  0.002201  0.023940 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              1.911e-01  1.858e-02  10.283 2.84e-06 ***
    ## distance_to_ocean_min_m -9.184e-06  9.111e-05  -0.101    0.922    
    ## max_depth               -2.449e-04  3.173e-04  -0.772    0.460    
    ## logArea                 -1.487e-03  1.645e-03  -0.904    0.390    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01023 on 9 degrees of freedom
    ## Multiple R-squared:  0.224,  Adjusted R-squared:  -0.0347 
    ## F-statistic: 0.8658 on 3 and 9 DF,  p-value: 0.4935

``` r
anova(PD_PDisp_geo_OM_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                         Df     Sum Sq    Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.00005863 5.8627e-05  0.5601 0.4733
    ## max_depth                1 0.00012780 1.2780e-04  1.2208 0.2979
    ## logArea                  1 0.00008549 8.5485e-05  0.8166 0.3897
    ## Residuals                9 0.00094212 1.0468e-04

``` r
p_values <- summary(PD_PDisp_geo_OM_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##            1.134274e-05            1.000000e+00            1.000000e+00 
    ##                 logArea 
    ##            1.000000e+00

``` r
# Stratified lakes and ocean sites
PD_PRic_geo_SO_model <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_stratified_sites,])
summary(PD_PRic_geo_SO_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_stratified_sites, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.63353 -0.16913  0.02623  0.29476  0.39100 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              7.316228   0.637554  11.475 4.44e-07 ***
    ## distance_to_ocean_min_m -0.005662   0.001275  -4.440  0.00125 ** 
    ## max_depth                0.001600   0.009906   0.162  0.87491    
    ## logArea                  0.045681   0.059535   0.767  0.46064    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3613 on 10 degrees of freedom
    ## Multiple R-squared:  0.7599, Adjusted R-squared:  0.6878 
    ## F-statistic: 10.55 on 3 and 10 DF,  p-value: 0.001933

``` r
anova(PD_PRic_geo_SO_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                         Df Sum Sq Mean Sq F value    Pr(>F)    
    ## distance_to_ocean_min_m  1 4.0120  4.0120 30.7268 0.0002464 ***
    ## max_depth                1 0.0428  0.0428  0.3280 0.5795115    
    ## logArea                  1 0.0769  0.0769  0.5887 0.4606420    
    ## Residuals               10 1.3057  0.1306                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_SO_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##            1.777047e-06            5.014832e-03            1.000000e+00 
    ##                 logArea 
    ##            1.000000e+00

``` r
PD_z_geo_SO_model <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_stratified_sites,])
summary(PD_z_geo_SO_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_stratified_sites, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.59838 -0.72653  0.01201  0.55790  1.68035 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)             -1.125010   2.047266  -0.550    0.595
    ## distance_to_ocean_min_m  0.004991   0.004095   1.219    0.251
    ## max_depth               -0.030543   0.031811  -0.960    0.360
    ## logArea                  0.017241   0.191176   0.090    0.930
    ## 
    ## Residual standard error: 1.16 on 10 degrees of freedom
    ## Multiple R-squared:  0.1687, Adjusted R-squared:  -0.08075 
    ## F-statistic: 0.6762 on 3 and 10 DF,  p-value: 0.5861

``` r
anova(PD_z_geo_SO_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1  1.2869 1.28693  0.9559 0.3513
    ## max_depth                1  1.4334 1.43341  1.0647 0.3265
    ## logArea                  1  0.0110 0.01095  0.0081 0.9299
    ## Residuals               10 13.4635 1.34635

``` r
p_values <- summary(PD_z_geo_SO_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                       1                       1                       1 
    ##                 logArea 
    ##                       1

``` r
PD_PDisp_geo_SO_model <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_stratified_sites,])
summary(PD_PDisp_geo_SO_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_stratified_sites, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.050938 -0.008060 -0.003566  0.013064  0.046248 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)              1.671e-01  4.861e-02   3.438  0.00635 **
    ## distance_to_ocean_min_m -1.223e-04  9.723e-05  -1.257  0.23719   
    ## max_depth                6.412e-04  7.554e-04   0.849  0.41585   
    ## logArea                 -4.052e-04  4.540e-03  -0.089  0.93065   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02755 on 10 degrees of freedom
    ## Multiple R-squared:  0.1664, Adjusted R-squared:  -0.08362 
    ## F-statistic: 0.6656 on 3 and 10 DF,  p-value: 0.592

``` r
anova(PD_PDisp_geo_SO_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                         Df    Sum Sq    Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.0008849 0.00088486  1.1656 0.3057
    ## max_depth                1 0.0006250 0.00062497  0.8232 0.3856
    ## logArea                  1 0.0000060 0.00000605  0.0080 0.9306
    ## Residuals               10 0.0075916 0.00075916

``` r
p_values <- summary(PD_PDisp_geo_SO_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.0254187               0.9487510               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
# Mixed lakes
PD_PRic_geo_M_model <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_lakes_geo,])
summary(PD_PRic_geo_M_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_lakes_geo, ])
    ## 
    ## Residuals:
    ##       FLK       HLO       MLN       NLN       NLU       OLO       ULN 
    ##  0.173264 -0.200974  0.079429  0.134824  0.004932 -0.301321  0.109845 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)              6.5196514  0.8802290   7.407  0.00509 **
    ## distance_to_ocean_min_m -0.0036623  0.0045323  -0.808  0.47820   
    ## max_depth               -0.0006091  0.0185026  -0.033  0.97581   
    ## logArea                  0.1779275  0.1054215   1.688  0.19004   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2568 on 3 degrees of freedom
    ## Multiple R-squared:  0.7179, Adjusted R-squared:  0.4358 
    ## F-statistic: 2.545 on 3 and 3 DF,  p-value: 0.2316

``` r
anova(PD_PRic_geo_M_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                         Df  Sum Sq  Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.19091 0.190908  2.8957 0.1874
    ## max_depth                1 0.12468 0.124676  1.8911 0.2628
    ## logArea                  1 0.18780 0.187799  2.8486 0.1900
    ## Residuals                3 0.19778 0.065927

``` r
p_values <- summary(PD_PRic_geo_M_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              0.02036363              1.00000000              1.00000000 
    ##                 logArea 
    ##              0.76015675

``` r
PD_z_geo_M_model <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_lakes_geo,])
summary(PD_z_geo_M_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_lakes_geo, ])
    ## 
    ## Residuals:
    ##      FLK      HLO      MLN      NLN      NLU      OLO      ULN 
    ## -0.39164 -0.73397  0.23629  1.16555  0.34108  0.07677 -0.69408 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)             -2.378054   3.258512  -0.730    0.518
    ## distance_to_ocean_min_m  0.001892   0.016778   0.113    0.917
    ## max_depth               -0.065705   0.068495  -0.959    0.408
    ## logArea                  0.251911   0.390259   0.645    0.565
    ## 
    ## Residual standard error: 0.9505 on 3 degrees of freedom
    ## Multiple R-squared:  0.3847, Adjusted R-squared:  -0.2307 
    ## F-statistic: 0.6251 on 3 and 3 DF,  p-value: 0.6456

``` r
anova(PD_z_geo_M_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.86165 0.86165  0.9537 0.4008
    ## max_depth                1 0.45618 0.45618  0.5049 0.5286
    ## logArea                  1 0.37645 0.37645  0.4167 0.5646
    ## Residuals                3 2.71041 0.90347

``` r
p_values <- summary(PD_z_geo_M_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                       1                       1                       1 
    ##                 logArea 
    ##                       1

``` r
PD_PDisp_geo_M_model <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_lakes_geo,])
summary(PD_PDisp_geo_M_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_lakes_geo, ])
    ## 
    ## Residuals:
    ##       FLK       HLO       MLN       NLN       NLU       OLO       ULN 
    ##  0.002631 -0.008280  0.001062  0.014305  0.002095 -0.007256 -0.004558 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              0.1546430  0.0375045   4.123   0.0259 *
    ## distance_to_ocean_min_m -0.0001868  0.0001931  -0.967   0.4048  
    ## max_depth               -0.0007955  0.0007884  -1.009   0.3872  
    ## logArea                  0.0045731  0.0044918   1.018   0.3836  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01094 on 3 degrees of freedom
    ## Multiple R-squared:  0.3451, Adjusted R-squared:  -0.3097 
    ## F-statistic: 0.5271 on 3 and 3 DF,  p-value: 0.694

``` r
anova(PD_PDisp_geo_M_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                         Df     Sum Sq    Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.00003957 3.9573e-05  0.3306 0.6056
    ## max_depth                1 0.00002561 2.5612e-05  0.2140 0.6751
    ## logArea                  1 0.00012406 1.2406e-04  1.0365 0.3836
    ## Residuals                3 0.00035906 1.1969e-04

``` r
p_values <- summary(PD_PDisp_geo_M_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.1034443               1.0000000               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
# Stratified lakes
PD_PRic_geo_S_model <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[stratified_lakes,])
summary(PD_PRic_geo_S_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##      BCM      CLM      GLK      HLM      NLK      OTM      SLN      TLN 
    ##  0.11905 -0.05492 -0.18213  0.29714  0.07181 -0.40328 -0.05091  0.20324 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              5.821097   1.822485   3.194   0.0331 *
    ## distance_to_ocean_min_m -0.002957   0.001549  -1.909   0.1289  
    ## max_depth               -0.008760   0.021996  -0.398   0.7108  
    ## logArea                  0.164865   0.225457   0.731   0.5052  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2959 on 4 degrees of freedom
    ## Multiple R-squared:  0.4906, Adjusted R-squared:  0.1086 
    ## F-statistic: 1.284 on 3 and 4 DF,  p-value: 0.3938

``` r
anova(PD_PRic_geo_S_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                         Df  Sum Sq  Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.25797 0.257971  2.9454 0.1613
    ## max_depth                1 0.03264 0.032640  0.3727 0.5745
    ## logArea                  1 0.04683 0.046834  0.5347 0.5052
    ## Residuals                4 0.35034 0.087586

``` r
p_values <- summary(PD_PRic_geo_S_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.1323530               0.5154496               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
PD_z_geo_S_model <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[stratified_lakes,])
summary(PD_z_geo_S_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##      BCM      CLM      GLK      HLM      NLK      OTM      SLN      TLN 
    ## -0.12449 -0.24694  0.59937 -0.58482  0.07139  0.72921  0.18406 -0.62777 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -10.973624   4.064548  -2.700   0.0541 .
    ## distance_to_ocean_min_m  -0.000481   0.003454  -0.139   0.8960  
    ## max_depth                -0.079262   0.049056  -1.616   0.1815  
    ## logArea                   1.180681   0.502819   2.348   0.0787 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.66 on 4 degrees of freedom
    ## Multiple R-squared:  0.6433, Adjusted R-squared:  0.3758 
    ## F-statistic: 2.405 on 3 and 4 DF,  p-value: 0.2079

``` r
anova(PD_z_geo_S_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##                         Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## distance_to_ocean_min_m  1 0.13873 0.13873  0.3185 0.60268  
    ## max_depth                1 0.60222 0.60222  1.3824 0.30490  
    ## logArea                  1 2.40199 2.40199  5.5137 0.07868 .
    ## Residuals                4 1.74257 0.43564                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_geo_S_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.2164135               1.0000000               0.7258090 
    ##                 logArea 
    ##               0.3147069

``` r
PD_PDisp_geo_S_model <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[stratified_lakes,])
summary(PD_PDisp_geo_S_model)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##       BCM       CLM       GLK       HLM       NLK       OTM       SLN       TLN 
    ##  0.002563 -0.008611  0.028940 -0.025547  0.020337 -0.002319 -0.011485 -0.003878 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -0.1881649  0.1423307  -1.322   0.2567  
    ## distance_to_ocean_min_m -0.0002779  0.0001209  -2.298   0.0831 .
    ## max_depth               -0.0025469  0.0017178  -1.483   0.2123  
    ## logArea                  0.0441519  0.0176075   2.508   0.0662 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02311 on 4 degrees of freedom
    ## Multiple R-squared:  0.7488, Adjusted R-squared:  0.5605 
    ## F-statistic: 3.975 on 3 and 4 DF,  p-value: 0.1079

``` r
anova(PD_PDisp_geo_S_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dispersion
    ##                         Df    Sum Sq   Mean Sq F value  Pr(>F)  
    ## distance_to_ocean_min_m  1 0.0012404 0.0012404  2.3220 0.20223  
    ## max_depth                1 0.0017715 0.0017715  3.3161 0.14271  
    ## logArea                  1 0.0033590 0.0033590  6.2879 0.06623 .
    ## Residuals                4 0.0021368 0.0005342                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_geo_S_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               0.3324746               0.8492304 
    ##                 logArea 
    ##               0.2649182

### PD alpha z-scores with env and geo correlated variables

``` r
# Log Area
PD_alpha_LA_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.z, x = logArea, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = stree_sespd_env[,1], size = 5, point.padding = 3) +
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
  labs(x="Log Area (m"^"2"~")", y="PD alpha z-score", colour = "Site type", fill = "Site type", tag = "c")
(PD_alpha_LA_plot <- PD_alpha_LA_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_alpha_LA_plot.jpg", plot = PD_alpha_LA_plot, width = 6.26, height = 6, units = "in")

# Distance from the ocean mean
PD_alpha_D_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.z, x = distance_to_ocean_min_m, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = stree_sespd_env[,1], size = 5, point.padding = 3) +
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
  labs(x="Isolation (m)", y="PD alpha z-score", colour = "Site type", fill = "Site type", tag = "b")
(PD_alpha_D_plot <- PD_alpha_D_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_alpha_D_plot.jpg", plot = PD_alpha_D_plot, width = 6.26, height = 6, units = "in")

# Max depth
PD_alpha_MD_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.z, x = max_depth, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = stree_sespd_env[,1], size = 5, point.padding = 3) +
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
  labs(x="Age (m)", y="PD alpha z-score", colour = "Site type", fill = "Site type", tag = "a")
(PD_alpha_MD_plot <- PD_alpha_MD_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_alpha_MD_plot.jpg", plot = PD_alpha_MD_plot, width = 6.26, height = 6, units = "in")

# Phylogenetic richness and dispersion
PD_alpha_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.z, x = Dispersion, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = stree_sespd_env[,1], size = 5, point.padding = 3) +
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
  labs(x="PD Dispersion", y="PD alpha z-score", colour = "Site type", fill = "Site type", tag = "a")
(PD_alpha_plot <- PD_alpha_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20plots-4.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_alpha_plot.jpg", plot = PD_alpha_plot, width = 4.88, height = 6, units = "in")
```

    ## Warning: ggrepel: 1 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

## PD beta diversity

### PD beta diversity distance calculations plus dispersion and PERMANOVA

``` r
### Regular
# Reference
PD_beta_ref_dist<- BAT::beta(presabs_lake, tree, abund = F)
(PD_beta_ref_BD <- betadisper(PD_beta_ref_dist$Btotal, stratification_group_ref))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = PD_beta_ref_dist$Btotal, group =
    ## stratification_group_ref)
    ## 
    ## No. of Positive Eigenvalues: 21
    ## No. of Negative Eigenvalues: 1
    ## 
    ## Average distance to median:
    ##  Reference      Ocean      Mixed Stratified 
    ##     0.0000     0.3888     0.3810     0.2954 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 22 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.8527 0.7216 0.5556 0.2855 0.2604 0.2334 0.1974 0.1937

``` r
(PD_beta_ref_AOV <- anova(PD_beta_ref_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value   Pr(>F)   
    ## Groups     3 0.15894 0.052981  7.5847 0.001563 **
    ## Residuals 19 0.13272 0.006985                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(PD_beta_ref_THSD <- TukeyHSD(PD_beta_ref_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                              diff         lwr        upr     p adj
    ## Ocean-Reference       0.388759022  0.13492092 0.64259712 0.0019786
    ## Mixed-Reference       0.380972871  0.13170881 0.63023693 0.0020177
    ## Stratified-Reference  0.295414817  0.04615075 0.54467888 0.0168613
    ## Mixed-Ocean          -0.007786151 -0.13470520 0.11913290 0.9981093
    ## Stratified-Ocean     -0.093344205 -0.22026325 0.03357484 0.1992503
    ## Stratified-Mixed     -0.085558054 -0.20306226 0.03194615 0.2061814

``` r
(PD_beta_ref_PM <- adonis2(PD_beta_ref_dist$Btotal ~ env[,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = PD_beta_ref_dist$Btotal ~ env[, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     3   2.5315 0.46901 5.5941  0.001 ***
    ## Residual 19   2.8660 0.53099                  
    ## Total    22   5.3974 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(PD_beta_ref_PM_pair <- pairwise.adonis(PD_beta_ref_dist$Btotal, env[,19], p.adjust.m = "bonferroni", perm = 999))
```

    ## Set of permutations < 'minperm'. Generating entire set.

    ##                     pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted
    ## 1     Stratified vs Mixed  1 1.2940870 9.415922 0.4021162   0.002      0.012
    ## 2     Stratified vs Ocean  1 1.1858559 8.357071 0.4105242   0.001      0.006
    ## 3 Stratified vs Reference  1 0.7729479 7.110647 0.5039207   0.112      0.672
    ## 4          Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.104      0.624
    ## 5      Mixed vs Reference  1 0.6592270 3.967203 0.3617334   0.133      0.798
    ## 6      Ocean vs Reference  1 0.6319662 3.354877 0.4015472   0.145      0.870
    ##   sig
    ## 1   .
    ## 2   *
    ## 3    
    ## 4    
    ## 5    
    ## 6

``` r
# Surveyed sites
PD_beta_dist<- BAT::beta(surveyed_sites_lake, stree, abund = F)
(PD_beta_BD <- betadisper(PD_beta_dist$Btotal, stratification_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = PD_beta_dist$Btotal, group = stratification_group)
    ## 
    ## No. of Positive Eigenvalues: 20
    ## No. of Negative Eigenvalues: 1
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.3888     0.3810     0.2954 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 21 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.8451 0.5956 0.2940 0.2611 0.2343 0.2009 0.1951 0.1567

``` r
(PD_beta_AOV <- anova(PD_beta_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value  Pr(>F)  
    ## Groups     2 0.040437 0.0202187  2.8944 0.07993 .
    ## Residuals 19 0.132722 0.0069854                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(PD_beta_THSD <- TukeyHSD(PD_beta_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                          diff        lwr        upr     p adj
    ## Mixed-Ocean      -0.007786117 -0.1224560 0.10688378 0.9837439
    ## Stratified-Ocean -0.093343777 -0.2080137 0.02132612 0.1234745
    ## Stratified-Mixed -0.085557660 -0.1917214 0.02060604 0.1281442

``` r
(PD_beta_PM <- adonis2(PD_beta_dist$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = PD_beta_dist$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2     F Pr(>F)    
    ## Model     2   1.8596 0.39351 6.164  0.001 ***
    ## Residual 19   2.8660 0.60649                 
    ## Total    21   4.7255 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(PD_beta_PM_pair <- pairwise.adonis(PD_beta_dist$Btotal, env[surveyed_sites,19], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.2940870 9.415922 0.4021162   0.002      0.006   *
    ## 2 Stratified vs Ocean  1 1.1858559 8.357071 0.4105242   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.106      0.318

``` r
# Without LCN
PD_beta_wo_LCN_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_LCN,], stree, abund = F)
(PD_beta_wo_LCN_BD <- betadisper(PD_beta_wo_LCN_dist$Btotal, stratification_group[-8]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = PD_beta_wo_LCN_dist$Btotal, group =
    ## stratification_group[-8])
    ## 
    ## No. of Positive Eigenvalues: 19
    ## No. of Negative Eigenvalues: 1
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.3528     0.3810     0.2954 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 20 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.8416 0.5421 0.2862 0.2610 0.2193 0.2008 0.1669 0.1561

``` r
(PD_beta_wo_LCN_AOV <- anova(PD_beta_wo_LCN_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.030094 0.0150469  2.3046 0.1285
    ## Residuals 18 0.117526 0.0065292

``` r
(PD_beta_wo_LCN_THSD <- TukeyHSD(PD_beta_wo_LCN_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff         lwr        upr     p adj
    ## Mixed-Ocean       0.02816227 -0.08940337 0.14572791 0.8158143
    ## Stratified-Ocean -0.05739477 -0.17496041 0.06017087 0.4426199
    ## Stratified-Mixed -0.08555704 -0.18866893 0.01755485 0.1142581

``` r
(PD_beta_wo_LCN_PM <- adonis2(PD_beta_wo_LCN_dist$Btotal ~ env[surveyed_sites_wo_LCN,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = PD_beta_wo_LCN_dist$Btotal ~ env[surveyed_sites_wo_LCN, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   1.9412 0.43018 6.7946  0.001 ***
    ## Residual 18   2.5713 0.56982                  
    ## Total    20   4.5125 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(PD_beta_wo_LCN_PM_pair <- pairwise.adonis(PD_beta_wo_LCN_dist$Btotal, env[surveyed_sites_wo_LCN,19], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.2940870 9.415922 0.4021162   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.2083230 9.439158 0.4618174   0.003      0.009   *
    ## 3      Mixed vs Ocean  1 0.3347638 2.034034 0.1560556   0.022      0.066

``` r
# Without TLN and HLM
PD_beta_wo_TLN_HLM_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_TLN_HLM,], stree, abund = F)
(PD_beta_wo_TLN_HLM_BD <- betadisper(PD_beta_wo_TLN_HLM_dist$Btotal, stratification_group[-c(5,21)]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = PD_beta_wo_TLN_HLM_dist$Btotal, group =
    ## stratification_group[-c(5, 21)])
    ## 
    ## No. of Positive Eigenvalues: 19
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.3888     0.3810     0.2405 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 19 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.7597 0.5513 0.2751 0.2560 0.2076 0.1856 0.1688 0.1485

``` r
(PD_beta_wo_TLN_HLM_AOV <- anova(PD_beta_wo_TLN_HLM_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq  Mean Sq F value   Pr(>F)   
    ## Groups     2 0.087067 0.043533  7.6495 0.004273 **
    ## Residuals 17 0.096748 0.005691                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(PD_beta_wo_TLN_HLM_THSD <- TukeyHSD(PD_beta_wo_TLN_HLM_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                          diff        lwr         upr     p adj
    ## Mixed-Ocean      -0.007785914 -0.1123030  0.09673115 0.9800914
    ## Stratified-Ocean -0.148257109 -0.2599906 -0.03652367 0.0089809
    ## Stratified-Mixed -0.140471196 -0.2449883 -0.03595413 0.0081854

``` r
(PD_beta_wo_TLN_HLM_PM <- adonis2(PD_beta_wo_TLN_HLM_dist$Btotal ~ env[surveyed_sites_wo_TLN_HLM,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = PD_beta_wo_TLN_HLM_dist$Btotal ~ env[surveyed_sites_wo_TLN_HLM, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   1.8539 0.42732 6.3425  0.001 ***
    ## Residual 17   2.4845 0.57268                  
    ## Total    19   4.3383 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(PD_beta_wo_TLN_HLM_PM_pair <- pairwise.adonis(PD_beta_wo_TLN_HLM_dist$Btotal, env[surveyed_sites_wo_TLN_HLM,19], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.3499033 10.500850 0.4666868   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.2146330  9.192717 0.4789690   0.003      0.009   *
    ## 3      Mixed vs Ocean  1 0.2573592  1.467099 0.1089395   0.105      0.315

``` r
# Mixed and stratified lakes
PD_beta_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], stree, abund = F)
# Ocean sites and mixed lakes
PD_beta_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites,], stree, abund = F)
# Stratified lakes and ocean sites
PD_beta_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites,], stree, abund = F)
# Stratified lakes
PD_beta_S_dist<- BAT::beta(surveyed_sites_lake[stratified_lakes,], stree, abund = F)


### Environmental
# Surveyed sites
PD_beta_env_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_env,], stree, abund = F)
# Mixed and stratified lakes
PD_beta_env_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], stree, abund = F)
# Ocean sites and mixed lakes
PD_beta_env_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites_env,], stree, abund = F)
# Stratified lakes and ocean sites
PD_beta_env_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites_env,], stree, abund = F)
# Mixed lakes
PD_beta_env_M_dist<- BAT::beta(surveyed_sites_lake[mixed_lakes,], stree, abund = F)


### Geographic
# Surveyed sites
PD_beta_geo_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_geo,], stree, abund = F)
# Mixed and stratified lakes
PD_beta_geo_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes_geo,], stree, abund = F)
# Ocean sites and mixed lakes
PD_beta_geo_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites_geo,], stree, abund = F)
# Stratified lakes and ocean sites
PD_beta_geo_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites,], stree, abund = F)
# Mixed lakes
PD_beta_geo_M_dist<- BAT::beta(surveyed_sites_lake[mixed_lakes_geo,], stree, abund = F)
```

### PD beta q values of dispersion statistics

``` r
# Sites
N <- 22
# Categories
k <- 3
# Sum of squares
SS <- PD_beta_AOV$`Sum Sq`[2]
# Mean of squares
MS <- PD_beta_AOV$`Mean Sq`[2]
# q value
q.value <- qtukey(p = 0.95, nmeans = k, df = N - k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

group <- PD_beta_THSD$group

(OM_q <- (abs(group[1,1]))/USE)
```

    ## [1] 0.2439479

``` r
(MS_q <- (abs(group[2,1]))/BSE)
```

    ## [1] 3.158893

``` r
(SO_q <- (abs(group[3,1]))/USE)
```

    ## [1] 2.680619

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx
```

### PD beta trait dispersions

- The dissimilarity between sites of the same site type

``` r
# Dispersion within-groups
PD_beta_BD_dist <- PD_beta_BD$distances
PD_beta_BD_dist <- as.data.frame(PD_beta_BD_dist)
PD_beta_BD_dist$X <- row.names(PD_beta_BD_dist)
PD_beta_BD_dist_env <- merge(PD_beta_BD_dist, env[surveyed_sites,], by = "X", sort = F)

PD_beta_BD_dist_env$Stratification <- factor(PD_beta_BD_dist_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

(PD_beta_BD_dist_env_plot <- ggplot(PD_beta_BD_dist_env, aes(x = Stratification, y = PD_beta_BD_dist, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.8,
            width = 0.1) +
    geom_text_repel(data = PD_beta_BD_dist_env, label = PD_beta_BD_dist_env$X, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
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
  labs(color = "Stratification", tag = "b"))
```

![](PD_analyses_files/figure-gfm/PD%20beta%20dispersion%20plot-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_beta_BD_dist_plot.jpg", PD_beta_BD_dist_env_plot, width = 4.88, height = 6, units = "in")
```

### PD beta NMDS

``` r
### Regular
PD_beta_ref_NMDS <- metaMDS(PD_beta_ref_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.905255e-05 
    ## Run 1 stress 9.914421e-05 
    ## ... Procrustes: rmse 0.0001715228  max resid 0.0004002377 
    ## ... Similar to previous best
    ## Run 2 stress 9.361446e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001438732  max resid 0.0003320175 
    ## ... Similar to previous best
    ## Run 3 stress 5.84523e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002013244  max resid 0.0003360322 
    ## ... Similar to previous best
    ## Run 4 stress 9.766058e-05 
    ## ... Procrustes: rmse 0.0001936725  max resid 0.0003377915 
    ## ... Similar to previous best
    ## Run 5 stress 9.988971e-05 
    ## ... Procrustes: rmse 0.0002496749  max resid 0.0004137558 
    ## ... Similar to previous best
    ## Run 6 stress 9.64562e-05 
    ## ... Procrustes: rmse 0.0002421935  max resid 0.0003983401 
    ## ... Similar to previous best
    ## Run 7 stress 9.677921e-05 
    ## ... Procrustes: rmse 0.0002065084  max resid 0.0003532497 
    ## ... Similar to previous best
    ## Run 8 stress 8.923641e-05 
    ## ... Procrustes: rmse 0.0001479651  max resid 0.0002621281 
    ## ... Similar to previous best
    ## Run 9 stress 9.707986e-05 
    ## ... Procrustes: rmse 0.0002104232  max resid 0.0003647147 
    ## ... Similar to previous best
    ## Run 10 stress 9.963734e-05 
    ## ... Procrustes: rmse 0.0002488593  max resid 0.0004197714 
    ## ... Similar to previous best
    ## Run 11 stress 9.643945e-05 
    ## ... Procrustes: rmse 0.0002136752  max resid 0.0003560257 
    ## ... Similar to previous best
    ## Run 12 stress 9.91623e-05 
    ## ... Procrustes: rmse 0.0002476373  max resid 0.0004063934 
    ## ... Similar to previous best
    ## Run 13 stress 9.915241e-05 
    ## ... Procrustes: rmse 0.0002053806  max resid 0.0003511837 
    ## ... Similar to previous best
    ## Run 14 stress 9.938781e-05 
    ## ... Procrustes: rmse 0.0002434591  max resid 0.0004008459 
    ## ... Similar to previous best
    ## Run 15 stress 9.918533e-05 
    ## ... Procrustes: rmse 0.0002359281  max resid 0.000401734 
    ## ... Similar to previous best
    ## Run 16 stress 9.546519e-05 
    ## ... Procrustes: rmse 0.0001544039  max resid 0.000261405 
    ## ... Similar to previous best
    ## Run 17 stress 9.730872e-05 
    ## ... Procrustes: rmse 0.0001965195  max resid 0.0003403949 
    ## ... Similar to previous best
    ## Run 18 stress 9.836727e-05 
    ## ... Procrustes: rmse 0.0001667585  max resid 0.0002835428 
    ## ... Similar to previous best
    ## Run 19 stress 8.355124e-05 
    ## ... Procrustes: rmse 4.868127e-05  max resid 8.492851e-05 
    ## ... Similar to previous best
    ## Run 20 stress 9.821285e-05 
    ## ... Procrustes: rmse 0.0002450422  max resid 0.0004124304 
    ## ... Similar to previous best
    ## Run 21 stress 9.692356e-05 
    ## ... Procrustes: rmse 6.66187e-05  max resid 0.0001243904 
    ## ... Similar to previous best
    ## Run 22 stress 9.994645e-05 
    ## ... Procrustes: rmse 6.170422e-05  max resid 0.0001119867 
    ## ... Similar to previous best
    ## Run 23 stress 9.384668e-05 
    ## ... Procrustes: rmse 0.0001734701  max resid 0.0003008543 
    ## ... Similar to previous best
    ## Run 24 stress 9.958735e-05 
    ## ... Procrustes: rmse 0.0002273447  max resid 0.0003813083 
    ## ... Similar to previous best
    ## Run 25 stress 9.958829e-05 
    ## ... Procrustes: rmse 0.0002213123  max resid 0.0003663487 
    ## ... Similar to previous best
    ## Run 26 stress 8.156469e-05 
    ## ... Procrustes: rmse 6.513962e-05  max resid 0.0001372599 
    ## ... Similar to previous best
    ## Run 27 stress 9.904479e-05 
    ## ... Procrustes: rmse 0.0002154225  max resid 0.0003626759 
    ## ... Similar to previous best
    ## Run 28 stress 9.782073e-05 
    ## ... Procrustes: rmse 0.0002147196  max resid 0.0003548778 
    ## ... Similar to previous best
    ## Run 29 stress 9.358353e-05 
    ## ... Procrustes: rmse 0.0002064258  max resid 0.0003436962 
    ## ... Similar to previous best
    ## Run 30 stress 9.981673e-05 
    ## ... Procrustes: rmse 0.000192443  max resid 0.0003134164 
    ## ... Similar to previous best
    ## Run 31 stress 9.863863e-05 
    ## ... Procrustes: rmse 0.0002469356  max resid 0.0004046323 
    ## ... Similar to previous best
    ## Run 32 stress 9.490975e-05 
    ## ... Procrustes: rmse 0.0002340336  max resid 0.0003962394 
    ## ... Similar to previous best
    ## Run 33 stress 9.871587e-05 
    ## ... Procrustes: rmse 0.0002451729  max resid 0.0004110005 
    ## ... Similar to previous best
    ## Run 34 stress 9.846503e-05 
    ## ... Procrustes: rmse 0.0002171865  max resid 0.0003705307 
    ## ... Similar to previous best
    ## Run 35 stress 9.988114e-05 
    ## ... Procrustes: rmse 0.000213733  max resid 0.0003569253 
    ## ... Similar to previous best
    ## Run 36 stress 9.547489e-05 
    ## ... Procrustes: rmse 5.219156e-05  max resid 0.0001174633 
    ## ... Similar to previous best
    ## Run 37 stress 6.518133e-05 
    ## ... Procrustes: rmse 4.038745e-05  max resid 7.325786e-05 
    ## ... Similar to previous best
    ## Run 38 stress 9.082894e-05 
    ## ... Procrustes: rmse 0.000187093  max resid 0.000308894 
    ## ... Similar to previous best
    ## Run 39 stress 9.588606e-05 
    ## ... Procrustes: rmse 0.0002121672  max resid 0.0003492806 
    ## ... Similar to previous best
    ## Run 40 stress 9.767007e-05 
    ## ... Procrustes: rmse 0.0002006211  max resid 0.0003426165 
    ## ... Similar to previous best
    ## Run 41 stress 9.926881e-05 
    ## ... Procrustes: rmse 0.0002494778  max resid 0.0004195758 
    ## ... Similar to previous best
    ## Run 42 stress 9.867849e-05 
    ## ... Procrustes: rmse 0.0002473001  max resid 0.0004104708 
    ## ... Similar to previous best
    ## Run 43 stress 9.682353e-05 
    ## ... Procrustes: rmse 0.0002421676  max resid 0.0003964214 
    ## ... Similar to previous best
    ## Run 44 stress 9.930458e-05 
    ## ... Procrustes: rmse 0.0002127006  max resid 0.0003587805 
    ## ... Similar to previous best
    ## Run 45 stress 9.776227e-05 
    ## ... Procrustes: rmse 0.0001334414  max resid 0.000234034 
    ## ... Similar to previous best
    ## Run 46 stress 4.76979e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 4.165038e-05  max resid 8.085781e-05 
    ## ... Similar to previous best
    ## Run 47 stress 6.907203e-05 
    ## ... Procrustes: rmse 4.198503e-05  max resid 7.168064e-05 
    ## ... Similar to previous best
    ## Run 48 stress 9.845017e-05 
    ## ... Procrustes: rmse 0.0002402233  max resid 0.0003918805 
    ## ... Similar to previous best
    ## Run 49 stress 6.236974e-05 
    ## ... Procrustes: rmse 8.137366e-05  max resid 0.0001166118 
    ## ... Similar to previous best
    ## Run 50 stress 9.948697e-05 
    ## ... Procrustes: rmse 0.0002457807  max resid 0.0004070562 
    ## ... Similar to previous best
    ## Run 51 stress 6.843759e-05 
    ## ... Procrustes: rmse 4.375428e-05  max resid 8.793526e-05 
    ## ... Similar to previous best
    ## Run 52 stress 9.83296e-05 
    ## ... Procrustes: rmse 0.0002064209  max resid 0.0003367308 
    ## ... Similar to previous best
    ## Run 53 stress 0.3793789 
    ## Run 54 stress 9.736428e-05 
    ## ... Procrustes: rmse 0.0002009939  max resid 0.0003215811 
    ## ... Similar to previous best
    ## Run 55 stress 9.931047e-05 
    ## ... Procrustes: rmse 0.0002203507  max resid 0.0003619711 
    ## ... Similar to previous best
    ## Run 56 stress 9.369296e-05 
    ## ... Procrustes: rmse 0.0001667117  max resid 0.0002717258 
    ## ... Similar to previous best
    ## Run 57 stress 4.203045e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 3.818258e-05  max resid 7.923687e-05 
    ## ... Similar to previous best
    ## Run 58 stress 9.932502e-05 
    ## ... Procrustes: rmse 0.0002234944  max resid 0.00036371 
    ## ... Similar to previous best
    ## Run 59 stress 9.838986e-05 
    ## ... Procrustes: rmse 0.0002114705  max resid 0.0003425311 
    ## ... Similar to previous best
    ## Run 60 stress 9.802339e-05 
    ## ... Procrustes: rmse 0.0002462227  max resid 0.0003912306 
    ## ... Similar to previous best
    ## Run 61 stress 9.784309e-05 
    ## ... Procrustes: rmse 0.0001958931  max resid 0.0003252639 
    ## ... Similar to previous best
    ## Run 62 stress 7.87867e-05 
    ## ... Procrustes: rmse 4.42767e-05  max resid 9.190704e-05 
    ## ... Similar to previous best
    ## Run 63 stress 9.737927e-05 
    ## ... Procrustes: rmse 0.0002152636  max resid 0.0003473428 
    ## ... Similar to previous best
    ## Run 64 stress 9.890984e-05 
    ## ... Procrustes: rmse 0.0002036596  max resid 0.0003474409 
    ## ... Similar to previous best
    ## Run 65 stress 9.865846e-05 
    ## ... Procrustes: rmse 0.0002508691  max resid 0.0003965625 
    ## ... Similar to previous best
    ## Run 66 stress 9.446465e-05 
    ## ... Procrustes: rmse 0.0001819542  max resid 0.0002893788 
    ## ... Similar to previous best
    ## Run 67 stress 9.882968e-05 
    ## ... Procrustes: rmse 0.000221074  max resid 0.0003580746 
    ## ... Similar to previous best
    ## Run 68 stress 9.476672e-05 
    ## ... Procrustes: rmse 0.0002378659  max resid 0.000383047 
    ## ... Similar to previous best
    ## Run 69 stress 7.933899e-05 
    ## ... Procrustes: rmse 4.258488e-05  max resid 6.887356e-05 
    ## ... Similar to previous best
    ## Run 70 stress 7.365152e-05 
    ## ... Procrustes: rmse 4.530727e-05  max resid 8.323035e-05 
    ## ... Similar to previous best
    ## Run 71 stress 9.935208e-05 
    ## ... Procrustes: rmse 0.000222433  max resid 0.0003590896 
    ## ... Similar to previous best
    ## Run 72 stress 8.981047e-05 
    ## ... Procrustes: rmse 4.439954e-05  max resid 7.510655e-05 
    ## ... Similar to previous best
    ## Run 73 stress 9.923248e-05 
    ## ... Procrustes: rmse 0.0002529554  max resid 0.0004027427 
    ## ... Similar to previous best
    ## Run 74 stress 8.92066e-05 
    ## ... Procrustes: rmse 0.0002252656  max resid 0.0003600202 
    ## ... Similar to previous best
    ## Run 75 stress 7.417478e-05 
    ## ... Procrustes: rmse 4.278034e-05  max resid 7.797964e-05 
    ## ... Similar to previous best
    ## Run 76 stress 5.951654e-05 
    ## ... Procrustes: rmse 2.924422e-05  max resid 5.445281e-05 
    ## ... Similar to previous best
    ## Run 77 stress 9.832933e-05 
    ## ... Procrustes: rmse 0.0002121705  max resid 0.0003423426 
    ## ... Similar to previous best
    ## Run 78 stress 9.117933e-05 
    ## ... Procrustes: rmse 0.000201951  max resid 0.0003240602 
    ## ... Similar to previous best
    ## Run 79 stress 9.651254e-05 
    ## ... Procrustes: rmse 0.0001617363  max resid 0.0002686254 
    ## ... Similar to previous best
    ## Run 80 stress 9.853172e-05 
    ## ... Procrustes: rmse 0.0002506433  max resid 0.0003992689 
    ## ... Similar to previous best
    ## Run 81 stress 6.630217e-05 
    ## ... Procrustes: rmse 6.447715e-05  max resid 0.0001168618 
    ## ... Similar to previous best
    ## Run 82 stress 9.129911e-05 
    ## ... Procrustes: rmse 0.0002018754  max resid 0.000330834 
    ## ... Similar to previous best
    ## Run 83 stress 9.938491e-05 
    ## ... Procrustes: rmse 0.0002167564  max resid 0.0003536507 
    ## ... Similar to previous best
    ## Run 84 stress 9.723269e-05 
    ## ... Procrustes: rmse 0.0002444183  max resid 0.0003932312 
    ## ... Similar to previous best
    ## Run 85 stress 9.843559e-05 
    ## ... Procrustes: rmse 0.0002177598  max resid 0.0003510621 
    ## ... Similar to previous best
    ## Run 86 stress 7.905685e-05 
    ## ... Procrustes: rmse 4.45648e-05  max resid 7.331565e-05 
    ## ... Similar to previous best
    ## Run 87 stress 9.847326e-05 
    ## ... Procrustes: rmse 0.0001982465  max resid 0.0003279245 
    ## ... Similar to previous best
    ## Run 88 stress 8.535775e-05 
    ## ... Procrustes: rmse 0.0001303292  max resid 0.0001990116 
    ## ... Similar to previous best
    ## Run 89 stress 9.594435e-05 
    ## ... Procrustes: rmse 8.148882e-05  max resid 0.0001389792 
    ## ... Similar to previous best
    ## Run 90 stress 6.093271e-05 
    ## ... Procrustes: rmse 3.357867e-05  max resid 6.713035e-05 
    ## ... Similar to previous best
    ## Run 91 stress 9.981803e-05 
    ## ... Procrustes: rmse 0.0002545098  max resid 0.0004062419 
    ## ... Similar to previous best
    ## Run 92 stress 9.374361e-05 
    ## ... Procrustes: rmse 5.172149e-05  max resid 8.598557e-05 
    ## ... Similar to previous best
    ## Run 93 stress 9.859541e-05 
    ## ... Procrustes: rmse 0.0002096988  max resid 0.0003435302 
    ## ... Similar to previous best
    ## Run 94 stress 9.977382e-05 
    ## ... Procrustes: rmse 0.0002317458  max resid 0.0003752088 
    ## ... Similar to previous best
    ## Run 95 stress 6.321976e-05 
    ## ... Procrustes: rmse 7.565796e-05  max resid 0.0001269929 
    ## ... Similar to previous best
    ## Run 96 stress 7.571675e-05 
    ## ... Procrustes: rmse 7.986176e-05  max resid 0.0001217647 
    ## ... Similar to previous best
    ## Run 97 stress 9.865346e-05 
    ## ... Procrustes: rmse 0.0002381252  max resid 0.0003742874 
    ## ... Similar to previous best
    ## Run 98 stress 9.790745e-05 
    ## ... Procrustes: rmse 0.0002468994  max resid 0.0003946442 
    ## ... Similar to previous best
    ## Run 99 stress 9.962767e-05 
    ## ... Procrustes: rmse 0.0002171889  max resid 0.0003507573 
    ## ... Similar to previous best
    ## Run 100 stress 8.250156e-05 
    ## ... Procrustes: rmse 4.959875e-05  max resid 0.0001008255 
    ## ... Similar to previous best
    ## Run 101 stress 9.649699e-05 
    ## ... Procrustes: rmse 0.0002131345  max resid 0.0003467906 
    ## ... Similar to previous best
    ## Run 102 stress 6.897299e-05 
    ## ... Procrustes: rmse 3.437696e-05  max resid 6.061815e-05 
    ## ... Similar to previous best
    ## Run 103 stress 9.883655e-05 
    ## ... Procrustes: rmse 0.0002174715  max resid 0.0003522897 
    ## ... Similar to previous best
    ## Run 104 stress 9.563675e-05 
    ## ... Procrustes: rmse 0.0002268601  max resid 0.0003716602 
    ## ... Similar to previous best
    ## Run 105 stress 9.933568e-05 
    ## ... Procrustes: rmse 0.0002517841  max resid 0.0004004324 
    ## ... Similar to previous best
    ## Run 106 stress 9.48978e-05 
    ## ... Procrustes: rmse 0.0002418719  max resid 0.0003846268 
    ## ... Similar to previous best
    ## Run 107 stress 9.73504e-05 
    ## ... Procrustes: rmse 0.0002092933  max resid 0.0003483515 
    ## ... Similar to previous best
    ## Run 108 stress 8.488543e-05 
    ## ... Procrustes: rmse 4.025917e-05  max resid 7.600229e-05 
    ## ... Similar to previous best
    ## Run 109 stress 8.472539e-05 
    ## ... Procrustes: rmse 3.935441e-05  max resid 6.846837e-05 
    ## ... Similar to previous best
    ## Run 110 stress 9.783359e-05 
    ## ... Procrustes: rmse 0.0002479707  max resid 0.0004007319 
    ## ... Similar to previous best
    ## Run 111 stress 9.143692e-05 
    ## ... Procrustes: rmse 0.0001796036  max resid 0.0002996866 
    ## ... Similar to previous best
    ## Run 112 stress 9.890157e-05 
    ## ... Procrustes: rmse 0.0002431739  max resid 0.0003788798 
    ## ... Similar to previous best
    ## Run 113 stress 6.41207e-05 
    ## ... Procrustes: rmse 3.254855e-05  max resid 5.705943e-05 
    ## ... Similar to previous best
    ## Run 114 stress 9.855147e-05 
    ## ... Procrustes: rmse 0.0002469398  max resid 0.0003957316 
    ## ... Similar to previous best
    ## Run 115 stress 9.764312e-05 
    ## ... Procrustes: rmse 0.0002187379  max resid 0.0003563914 
    ## ... Similar to previous best
    ## Run 116 stress 9.613234e-05 
    ## ... Procrustes: rmse 0.0001886444  max resid 0.0003106416 
    ## ... Similar to previous best
    ## Run 117 stress 9.827652e-05 
    ## ... Procrustes: rmse 0.0002019119  max resid 0.00032998 
    ## ... Similar to previous best
    ## Run 118 stress 9.910481e-05 
    ## ... Procrustes: rmse 0.0002504791  max resid 0.0004006533 
    ## ... Similar to previous best
    ## Run 119 stress 6.402113e-05 
    ## ... Procrustes: rmse 3.516926e-05  max resid 5.612888e-05 
    ## ... Similar to previous best
    ## Run 120 stress 9.872777e-05 
    ## ... Procrustes: rmse 0.0002493912  max resid 0.0003978233 
    ## ... Similar to previous best
    ## Run 121 stress 9.727931e-05 
    ## ... Procrustes: rmse 0.0002169086  max resid 0.0003503667 
    ## ... Similar to previous best
    ## Run 122 stress 9.845412e-05 
    ## ... Procrustes: rmse 0.0002449584  max resid 0.0003906089 
    ## ... Similar to previous best
    ## Run 123 stress 9.92206e-05 
    ## ... Procrustes: rmse 0.0002526424  max resid 0.0004020821 
    ## ... Similar to previous best
    ## Run 124 stress 9.715375e-05 
    ## ... Procrustes: rmse 0.0002448854  max resid 0.0003894695 
    ## ... Similar to previous best
    ## Run 125 stress 9.779554e-05 
    ## ... Procrustes: rmse 0.0002476619  max resid 0.0003964502 
    ## ... Similar to previous best
    ## Run 126 stress 9.900956e-05 
    ## ... Procrustes: rmse 0.0002486052  max resid 0.0003981261 
    ## ... Similar to previous best
    ## Run 127 stress 9.941214e-05 
    ## ... Procrustes: rmse 0.0002201365  max resid 0.0003604376 
    ## ... Similar to previous best
    ## Run 128 stress 9.883131e-05 
    ## ... Procrustes: rmse 0.0002226302  max resid 0.000361287 
    ## ... Similar to previous best
    ## Run 129 stress 9.568571e-05 
    ## ... Procrustes: rmse 0.0001794334  max resid 0.0002959836 
    ## ... Similar to previous best
    ## Run 130 stress 9.965068e-05 
    ## ... Procrustes: rmse 0.000217021  max resid 0.0003503269 
    ## ... Similar to previous best
    ## Run 131 stress 9.616887e-05 
    ## ... Procrustes: rmse 0.0002069009  max resid 0.000337211 
    ## ... Similar to previous best
    ## Run 132 stress 9.862349e-05 
    ## ... Procrustes: rmse 0.0002503788  max resid 0.000398649 
    ## ... Similar to previous best
    ## Run 133 stress 9.846093e-05 
    ## ... Procrustes: rmse 0.0001749007  max resid 0.0002838006 
    ## ... Similar to previous best
    ## Run 134 stress 9.995517e-05 
    ## ... Procrustes: rmse 4.870107e-05  max resid 9.241647e-05 
    ## ... Similar to previous best
    ## Run 135 stress 9.93033e-05 
    ## ... Procrustes: rmse 0.0002511912  max resid 0.0003997459 
    ## ... Similar to previous best
    ## Run 136 stress 9.963701e-05 
    ## ... Procrustes: rmse 0.0002444268  max resid 0.0003908108 
    ## ... Similar to previous best
    ## Run 137 stress 9.774546e-05 
    ## ... Procrustes: rmse 0.0002490303  max resid 0.0003993713 
    ## ... Similar to previous best
    ## Run 138 stress 9.995089e-05 
    ## ... Procrustes: rmse 0.000180485  max resid 0.0002952008 
    ## ... Similar to previous best
    ## Run 139 stress 9.505209e-05 
    ## ... Procrustes: rmse 0.0001709433  max resid 0.000261058 
    ## ... Similar to previous best
    ## Run 140 stress 9.900963e-05 
    ## ... Procrustes: rmse 0.0002199169  max resid 0.0003594887 
    ## ... Similar to previous best
    ## Run 141 stress 8.753121e-05 
    ## ... Procrustes: rmse 9.905779e-05  max resid 0.0001655418 
    ## ... Similar to previous best
    ## Run 142 stress 6.397792e-05 
    ## ... Procrustes: rmse 5.39425e-05  max resid 9.495857e-05 
    ## ... Similar to previous best
    ## Run 143 stress 9.939812e-05 
    ## ... Procrustes: rmse 0.0001211577  max resid 0.0002135037 
    ## ... Similar to previous best
    ## Run 144 stress 8.160155e-05 
    ## ... Procrustes: rmse 6.6004e-05  max resid 0.0001175193 
    ## ... Similar to previous best
    ## Run 145 stress 5.579229e-05 
    ## ... Procrustes: rmse 3.282227e-05  max resid 5.853397e-05 
    ## ... Similar to previous best
    ## Run 146 stress 9.63318e-05 
    ## ... Procrustes: rmse 0.0001946502  max resid 0.0003177099 
    ## ... Similar to previous best
    ## Run 147 stress 9.497027e-05 
    ## ... Procrustes: rmse 0.0002104006  max resid 0.0003390364 
    ## ... Similar to previous best
    ## Run 148 stress 8.780072e-05 
    ## ... Procrustes: rmse 0.0001943142  max resid 0.0003174835 
    ## ... Similar to previous best
    ## Run 149 stress 9.291791e-05 
    ## ... Procrustes: rmse 7.973474e-05  max resid 0.0001323763 
    ## ... Similar to previous best
    ## Run 150 stress 9.201299e-05 
    ## ... Procrustes: rmse 0.0001691299  max resid 0.0002733414 
    ## ... Similar to previous best
    ## Run 151 stress 9.230145e-05 
    ## ... Procrustes: rmse 5.259749e-05  max resid 9.389933e-05 
    ## ... Similar to previous best
    ## Run 152 stress 8.685239e-05 
    ## ... Procrustes: rmse 3.973477e-05  max resid 8.409432e-05 
    ## ... Similar to previous best
    ## Run 153 stress 9.619264e-05 
    ## ... Procrustes: rmse 0.0002432447  max resid 0.0003838453 
    ## ... Similar to previous best
    ## Run 154 stress 9.954041e-05 
    ## ... Procrustes: rmse 0.0002096262  max resid 0.0003384183 
    ## ... Similar to previous best
    ## Run 155 stress 5.656338e-05 
    ## ... Procrustes: rmse 2.92669e-05  max resid 5.279402e-05 
    ## ... Similar to previous best
    ## Run 156 stress 8.503983e-05 
    ## ... Procrustes: rmse 4.537451e-05  max resid 0.000137927 
    ## ... Similar to previous best
    ## Run 157 stress 9.801304e-05 
    ## ... Procrustes: rmse 0.0002129604  max resid 0.0003441443 
    ## ... Similar to previous best
    ## Run 158 stress 9.996332e-05 
    ## ... Procrustes: rmse 0.0002084527  max resid 0.0003448673 
    ## ... Similar to previous best
    ## Run 159 stress 9.999876e-05 
    ## ... Procrustes: rmse 0.0002127581  max resid 0.0003498057 
    ## ... Similar to previous best
    ## Run 160 stress 9.93062e-05 
    ## ... Procrustes: rmse 0.0002089971  max resid 0.0003390316 
    ## ... Similar to previous best
    ## Run 161 stress 9.287873e-05 
    ## ... Procrustes: rmse 0.0002269235  max resid 0.0003644704 
    ## ... Similar to previous best
    ## Run 162 stress 9.818652e-05 
    ## ... Procrustes: rmse 0.0002400037  max resid 0.0003810404 
    ## ... Similar to previous best
    ## Run 163 stress 7.73054e-05 
    ## ... Procrustes: rmse 4.601318e-05  max resid 0.000100102 
    ## ... Similar to previous best
    ## Run 164 stress 9.858006e-05 
    ## ... Procrustes: rmse 0.000247823  max resid 0.0003980351 
    ## ... Similar to previous best
    ## Run 165 stress 9.780254e-05 
    ## ... Procrustes: rmse 0.0002483989  max resid 0.0003981354 
    ## ... Similar to previous best
    ## Run 166 stress 9.641995e-05 
    ## ... Procrustes: rmse 0.0001744481  max resid 0.0002952752 
    ## ... Similar to previous best
    ## Run 167 stress 9.690062e-05 
    ## ... Procrustes: rmse 0.0002157226  max resid 0.0003478455 
    ## ... Similar to previous best
    ## Run 168 stress 9.766252e-05 
    ## ... Procrustes: rmse 0.0002360707  max resid 0.0003820481 
    ## ... Similar to previous best
    ## Run 169 stress 8.522932e-05 
    ## ... Procrustes: rmse 9.667774e-05  max resid 0.0001569532 
    ## ... Similar to previous best
    ## Run 170 stress 9.960654e-05 
    ## ... Procrustes: rmse 0.0001731554  max resid 0.0002989666 
    ## ... Similar to previous best
    ## Run 171 stress 9.831267e-05 
    ## ... Procrustes: rmse 0.0002142739  max resid 0.0003448251 
    ## ... Similar to previous best
    ## Run 172 stress 9.897279e-05 
    ## ... Procrustes: rmse 0.0002060367  max resid 0.000337303 
    ## ... Similar to previous best
    ## Run 173 stress 9.785362e-05 
    ## ... Procrustes: rmse 0.0002471301  max resid 0.0003963338 
    ## ... Similar to previous best
    ## Run 174 stress 9.907229e-05 
    ## ... Procrustes: rmse 0.0002155835  max resid 0.0003521104 
    ## ... Similar to previous best
    ## Run 175 stress 9.922479e-05 
    ## ... Procrustes: rmse 0.0002025786  max resid 0.0003289831 
    ## ... Similar to previous best
    ## Run 176 stress 9.869283e-05 
    ## ... Procrustes: rmse 0.0002232534  max resid 0.0003565729 
    ## ... Similar to previous best
    ## Run 177 stress 8.955053e-05 
    ## ... Procrustes: rmse 4.961706e-05  max resid 7.716003e-05 
    ## ... Similar to previous best
    ## Run 178 stress 8.759881e-05 
    ## ... Procrustes: rmse 9.149592e-05  max resid 0.0001334681 
    ## ... Similar to previous best
    ## Run 179 stress 8.370821e-05 
    ## ... Procrustes: rmse 4.210189e-05  max resid 8.893368e-05 
    ## ... Similar to previous best
    ## Run 180 stress 9.4353e-05 
    ## ... Procrustes: rmse 0.0001616369  max resid 0.0002687321 
    ## ... Similar to previous best
    ## Run 181 stress 9.871262e-05 
    ## ... Procrustes: rmse 0.0002107631  max resid 0.0003399742 
    ## ... Similar to previous best
    ## Run 182 stress 9.082328e-05 
    ## ... Procrustes: rmse 4.883723e-05  max resid 0.0001120046 
    ## ... Similar to previous best
    ## Run 183 stress 9.985234e-05 
    ## ... Procrustes: rmse 0.0002133953  max resid 0.0003473846 
    ## ... Similar to previous best
    ## Run 184 stress 9.823613e-05 
    ## ... Procrustes: rmse 0.0002504368  max resid 0.0003993194 
    ## ... Similar to previous best
    ## Run 185 stress 9.901223e-05 
    ## ... Procrustes: rmse 0.0002034435  max resid 0.0003422297 
    ## ... Similar to previous best
    ## Run 186 stress 4.789079e-05 
    ## ... Procrustes: rmse 2.637703e-05  max resid 4.464734e-05 
    ## ... Similar to previous best
    ## Run 187 stress 9.895466e-05 
    ## ... Procrustes: rmse 0.0002444056  max resid 0.0003901511 
    ## ... Similar to previous best
    ## Run 188 stress 7.880924e-05 
    ## ... Procrustes: rmse 6.890229e-05  max resid 0.0001102591 
    ## ... Similar to previous best
    ## Run 189 stress 9.780377e-05 
    ## ... Procrustes: rmse 0.0002204954  max resid 0.0003582753 
    ## ... Similar to previous best
    ## Run 190 stress 9.793427e-05 
    ## ... Procrustes: rmse 0.0002448359  max resid 0.0003869279 
    ## ... Similar to previous best
    ## Run 191 stress 9.765193e-05 
    ## ... Procrustes: rmse 0.0002469954  max resid 0.0003922557 
    ## ... Similar to previous best
    ## Run 192 stress 9.774626e-05 
    ## ... Procrustes: rmse 0.0002436395  max resid 0.0003836721 
    ## ... Similar to previous best
    ## Run 193 stress 0.3824882 
    ## Run 194 stress 4.180056e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 2.527196e-05  max resid 4.968737e-05 
    ## ... Similar to previous best
    ## Run 195 stress 5.464814e-05 
    ## ... Procrustes: rmse 2.755296e-05  max resid 6.370575e-05 
    ## ... Similar to previous best
    ## Run 196 stress 9.567199e-05 
    ## ... Procrustes: rmse 0.0002324441  max resid 0.0003944946 
    ## ... Similar to previous best
    ## Run 197 stress 9.677824e-05 
    ## ... Procrustes: rmse 0.0002106632  max resid 0.0003544276 
    ## ... Similar to previous best
    ## Run 198 stress 9.955788e-05 
    ## ... Procrustes: rmse 0.0002491484  max resid 0.0004054125 
    ## ... Similar to previous best
    ## Run 199 stress 9.665009e-05 
    ## ... Procrustes: rmse 0.000213041  max resid 0.0003626917 
    ## ... Similar to previous best
    ## Run 200 stress 9.866584e-05 
    ## ... Procrustes: rmse 5.606115e-05  max resid 8.965671e-05 
    ## ... Similar to previous best
    ## Run 201 stress 9.854966e-05 
    ## ... Procrustes: rmse 0.0002198564  max resid 0.0003760838 
    ## ... Similar to previous best
    ## Run 202 stress 9.839956e-05 
    ## ... Procrustes: rmse 0.0002205017  max resid 0.0003636513 
    ## ... Similar to previous best
    ## Run 203 stress 9.610592e-05 
    ## ... Procrustes: rmse 0.0002381935  max resid 0.0003906672 
    ## ... Similar to previous best
    ## Run 204 stress 9.686642e-05 
    ## ... Procrustes: rmse 0.0001996155  max resid 0.0003482057 
    ## ... Similar to previous best
    ## Run 205 stress 9.861216e-05 
    ## ... Procrustes: rmse 0.0002165511  max resid 0.0003606316 
    ## ... Similar to previous best
    ## Run 206 stress 9.9451e-05 
    ## ... Procrustes: rmse 0.0002186267  max resid 0.00036318 
    ## ... Similar to previous best
    ## Run 207 stress 9.663912e-05 
    ## ... Procrustes: rmse 0.0002219017  max resid 0.0003831429 
    ## ... Similar to previous best
    ## Run 208 stress 9.915558e-05 
    ## ... Procrustes: rmse 0.0002481286  max resid 0.0004099103 
    ## ... Similar to previous best
    ## Run 209 stress 9.116956e-05 
    ## ... Procrustes: rmse 0.0001245921  max resid 0.0002133284 
    ## ... Similar to previous best
    ## Run 210 stress 9.867265e-05 
    ## ... Procrustes: rmse 0.0002410097  max resid 0.0003934483 
    ## ... Similar to previous best
    ## Run 211 stress 9.979484e-05 
    ## ... Procrustes: rmse 0.0002510835  max resid 0.0004089848 
    ## ... Similar to previous best
    ## Run 212 stress 9.866299e-05 
    ## ... Procrustes: rmse 0.0002190375  max resid 0.0003582289 
    ## ... Similar to previous best
    ## Run 213 stress 9.918422e-05 
    ## ... Procrustes: rmse 0.0002053194  max resid 0.0003548881 
    ## ... Similar to previous best
    ## Run 214 stress 8.535108e-05 
    ## ... Procrustes: rmse 4.714579e-05  max resid 7.499395e-05 
    ## ... Similar to previous best
    ## Run 215 stress 6.897979e-05 
    ## ... Procrustes: rmse 4.290014e-05  max resid 7.23647e-05 
    ## ... Similar to previous best
    ## Run 216 stress 9.900806e-05 
    ## ... Procrustes: rmse 0.0002135436  max resid 0.0003654682 
    ## ... Similar to previous best
    ## Run 217 stress 5.139696e-05 
    ## ... Procrustes: rmse 2.865684e-05  max resid 4.721031e-05 
    ## ... Similar to previous best
    ## Run 218 stress 0.3793776 
    ## Run 219 stress 5.842507e-05 
    ## ... Procrustes: rmse 3.388664e-05  max resid 5.632476e-05 
    ## ... Similar to previous best
    ## Run 220 stress 9.928913e-05 
    ## ... Procrustes: rmse 0.0002220403  max resid 0.000366391 
    ## ... Similar to previous best
    ## Run 221 stress 9.899302e-05 
    ## ... Procrustes: rmse 0.0002502236  max resid 0.0004068848 
    ## ... Similar to previous best
    ## Run 222 stress 9.702508e-05 
    ## ... Procrustes: rmse 0.0001537773  max resid 0.0002607847 
    ## ... Similar to previous best
    ## Run 223 stress 9.909686e-05 
    ## ... Procrustes: rmse 0.000226261  max resid 0.000379019 
    ## ... Similar to previous best
    ## Run 224 stress 9.819382e-05 
    ## ... Procrustes: rmse 0.0002194129  max resid 0.0003645547 
    ## ... Similar to previous best
    ## Run 225 stress 6.020998e-05 
    ## ... Procrustes: rmse 3.231795e-05  max resid 6.29543e-05 
    ## ... Similar to previous best
    ## Run 226 stress 7.582441e-05 
    ## ... Procrustes: rmse 0.0001202034  max resid 0.0002100665 
    ## ... Similar to previous best
    ## Run 227 stress 9.537694e-05 
    ## ... Procrustes: rmse 0.000147322  max resid 0.0002612225 
    ## ... Similar to previous best
    ## Run 228 stress 8.984686e-05 
    ## ... Procrustes: rmse 0.0001103516  max resid 0.0002008346 
    ## ... Similar to previous best
    ## Run 229 stress 9.999353e-05 
    ## ... Procrustes: rmse 0.0002228633  max resid 0.0003687927 
    ## ... Similar to previous best
    ## Run 230 stress 9.971577e-05 
    ## ... Procrustes: rmse 0.0002470562  max resid 0.0004023144 
    ## ... Similar to previous best
    ## Run 231 stress 9.382157e-05 
    ## ... Procrustes: rmse 3.915302e-05  max resid 7.977584e-05 
    ## ... Similar to previous best
    ## Run 232 stress 9.804591e-05 
    ## ... Procrustes: rmse 0.0002154025  max resid 0.0003578533 
    ## ... Similar to previous best
    ## Run 233 stress 9.531641e-05 
    ## ... Procrustes: rmse 5.313861e-05  max resid 0.0001171431 
    ## ... Similar to previous best
    ## Run 234 stress 9.934776e-05 
    ## ... Procrustes: rmse 0.0002453688  max resid 0.0003964184 
    ## ... Similar to previous best
    ## Run 235 stress 7.331269e-05 
    ## ... Procrustes: rmse 3.763816e-05  max resid 6.879726e-05 
    ## ... Similar to previous best
    ## Run 236 stress 9.784195e-05 
    ## ... Procrustes: rmse 0.0002093752  max resid 0.0003527529 
    ## ... Similar to previous best
    ## Run 237 stress 9.306122e-05 
    ## ... Procrustes: rmse 4.125051e-05  max resid 8.587265e-05 
    ## ... Similar to previous best
    ## Run 238 stress 8.564407e-05 
    ## ... Procrustes: rmse 4.966475e-05  max resid 8.121047e-05 
    ## ... Similar to previous best
    ## Run 239 stress 9.777189e-05 
    ## ... Procrustes: rmse 0.0002477618  max resid 0.0004041906 
    ## ... Similar to previous best
    ## Run 240 stress 9.492658e-05 
    ## ... Procrustes: rmse 6.262646e-05  max resid 0.0001235255 
    ## ... Similar to previous best
    ## Run 241 stress 9.763212e-05 
    ## ... Procrustes: rmse 0.0002044789  max resid 0.0003405615 
    ## ... Similar to previous best
    ## Run 242 stress 9.967706e-05 
    ## ... Procrustes: rmse 0.0002430733  max resid 0.0004091239 
    ## ... Similar to previous best
    ## Run 243 stress 4.108328e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 2.608123e-05  max resid 4.810216e-05 
    ## ... Similar to previous best
    ## Run 244 stress 9.333089e-05 
    ## ... Procrustes: rmse 4.479936e-05  max resid 9.174324e-05 
    ## ... Similar to previous best
    ## Run 245 stress 9.924361e-05 
    ## ... Procrustes: rmse 0.0002233074  max resid 0.0003686095 
    ## ... Similar to previous best
    ## Run 246 stress 9.832258e-05 
    ## ... Procrustes: rmse 0.0002370615  max resid 0.0003862502 
    ## ... Similar to previous best
    ## Run 247 stress 7.920357e-05 
    ## ... Procrustes: rmse 4.398014e-05  max resid 6.829234e-05 
    ## ... Similar to previous best
    ## Run 248 stress 9.848104e-05 
    ## ... Procrustes: rmse 0.0002546054  max resid 0.0004184163 
    ## ... Similar to previous best
    ## Run 249 stress 8.100175e-05 
    ## ... Procrustes: rmse 5.698666e-05  max resid 0.000111604 
    ## ... Similar to previous best
    ## Run 250 stress 7.967305e-05 
    ## ... Procrustes: rmse 4.696953e-05  max resid 7.668427e-05 
    ## ... Similar to previous best
    ## Run 251 stress 9.53254e-05 
    ## ... Procrustes: rmse 0.0002244058  max resid 0.0003761417 
    ## ... Similar to previous best
    ## Run 252 stress 9.788712e-05 
    ## ... Procrustes: rmse 0.0002525827  max resid 0.000414429 
    ## ... Similar to previous best
    ## Run 253 stress 9.910845e-05 
    ## ... Procrustes: rmse 0.00025941  max resid 0.000422838 
    ## ... Similar to previous best
    ## Run 254 stress 9.495009e-05 
    ## ... Procrustes: rmse 0.000204102  max resid 0.0003418032 
    ## ... Similar to previous best
    ## Run 255 stress 9.824176e-05 
    ## ... Procrustes: rmse 0.0002257727  max resid 0.0003779212 
    ## ... Similar to previous best
    ## Run 256 stress 9.451148e-05 
    ## ... Procrustes: rmse 0.0001739354  max resid 0.0003046213 
    ## ... Similar to previous best
    ## Run 257 stress 9.937382e-05 
    ## ... Procrustes: rmse 0.0002070449  max resid 0.0003494929 
    ## ... Similar to previous best
    ## Run 258 stress 9.727726e-05 
    ## ... Procrustes: rmse 0.0001883935  max resid 0.0003153418 
    ## ... Similar to previous best
    ## Run 259 stress 6.109908e-05 
    ## ... Procrustes: rmse 8.963384e-05  max resid 0.0001473584 
    ## ... Similar to previous best
    ## Run 260 stress 9.704889e-05 
    ## ... Procrustes: rmse 0.0002554752  max resid 0.0004179659 
    ## ... Similar to previous best
    ## Run 261 stress 9.578216e-05 
    ## ... Procrustes: rmse 0.000218184  max resid 0.0003683417 
    ## ... Similar to previous best
    ## Run 262 stress 9.757346e-05 
    ## ... Procrustes: rmse 0.0002572798  max resid 0.0004191695 
    ## ... Similar to previous best
    ## Run 263 stress 9.995287e-05 
    ## ... Procrustes: rmse 0.0002237675  max resid 0.0003783791 
    ## ... Similar to previous best
    ## Run 264 stress 9.635177e-05 
    ## ... Procrustes: rmse 0.0002555879  max resid 0.0004180856 
    ## ... Similar to previous best
    ## Run 265 stress 9.567011e-05 
    ## ... Procrustes: rmse 0.0001884476  max resid 0.0003127639 
    ## ... Similar to previous best
    ## Run 266 stress 9.504568e-05 
    ## ... Procrustes: rmse 0.0002040773  max resid 0.0003325924 
    ## ... Similar to previous best
    ## Run 267 stress 9.954307e-05 
    ## ... Procrustes: rmse 0.0002285016  max resid 0.0003824144 
    ## ... Similar to previous best
    ## Run 268 stress 9.741223e-05 
    ## ... Procrustes: rmse 0.0002495121  max resid 0.0004122483 
    ## ... Similar to previous best
    ## Run 269 stress 9.550616e-05 
    ## ... Procrustes: rmse 0.0002229965  max resid 0.0003736325 
    ## ... Similar to previous best
    ## Run 270 stress 9.850189e-05 
    ## ... Procrustes: rmse 0.0002476888  max resid 0.0004059188 
    ## ... Similar to previous best
    ## Run 271 stress 8.933948e-05 
    ## ... Procrustes: rmse 5.083331e-05  max resid 8.012154e-05 
    ## ... Similar to previous best
    ## Run 272 stress 9.964555e-05 
    ## ... Procrustes: rmse 0.0002512076  max resid 0.0004122138 
    ## ... Similar to previous best
    ## Run 273 stress 9.853997e-05 
    ## ... Procrustes: rmse 0.0002433171  max resid 0.0004003431 
    ## ... Similar to previous best
    ## Run 274 stress 9.896556e-05 
    ## ... Procrustes: rmse 0.0002338845  max resid 0.0003925054 
    ## ... Similar to previous best
    ## Run 275 stress 9.640157e-05 
    ## ... Procrustes: rmse 0.000243333  max resid 0.0003993063 
    ## ... Similar to previous best
    ## Run 276 stress 7.812533e-05 
    ## ... Procrustes: rmse 0.0001356737  max resid 0.0002345501 
    ## ... Similar to previous best
    ## Run 277 stress 9.657081e-05 
    ## ... Procrustes: rmse 0.0002190326  max resid 0.0003648157 
    ## ... Similar to previous best
    ## Run 278 stress 9.887502e-05 
    ## ... Procrustes: rmse 0.0001943682  max resid 0.0003235815 
    ## ... Similar to previous best
    ## Run 279 stress 9.9292e-05 
    ## ... Procrustes: rmse 0.0002084328  max resid 0.0003403425 
    ## ... Similar to previous best
    ## Run 280 stress 9.834114e-05 
    ## ... Procrustes: rmse 0.000227759  max resid 0.0003762387 
    ## ... Similar to previous best
    ## Run 281 stress 6.102505e-05 
    ## ... Procrustes: rmse 2.758634e-05  max resid 5.538996e-05 
    ## ... Similar to previous best
    ## Run 282 stress 9.897064e-05 
    ## ... Procrustes: rmse 0.0002294076  max resid 0.0003852501 
    ## ... Similar to previous best
    ## Run 283 stress 9.802085e-05 
    ## ... Procrustes: rmse 0.0002538219  max resid 0.0004143746 
    ## ... Similar to previous best
    ## Run 284 stress 2.355037e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 2.19266e-05  max resid 3.552613e-05 
    ## ... Similar to previous best
    ## Run 285 stress 9.928003e-05 
    ## ... Procrustes: rmse 0.0002151648  max resid 0.0003541612 
    ## ... Similar to previous best
    ## Run 286 stress 9.861845e-05 
    ## ... Procrustes: rmse 0.0002174754  max resid 0.0003543028 
    ## ... Similar to previous best
    ## Run 287 stress 9.75559e-05 
    ## ... Procrustes: rmse 0.0002535622  max resid 0.0004073804 
    ## ... Similar to previous best
    ## Run 288 stress 9.746229e-05 
    ## ... Procrustes: rmse 0.000250109  max resid 0.0004074853 
    ## ... Similar to previous best
    ## Run 289 stress 7.40406e-05 
    ## ... Procrustes: rmse 3.552644e-05  max resid 7.163877e-05 
    ## ... Similar to previous best
    ## Run 290 stress 9.82789e-05 
    ## ... Procrustes: rmse 0.000252843  max resid 0.0004079571 
    ## ... Similar to previous best
    ## Run 291 stress 9.636496e-05 
    ## ... Procrustes: rmse 0.0002511745  max resid 0.0004017988 
    ## ... Similar to previous best
    ## Run 292 stress 8.641409e-05 
    ## ... Procrustes: rmse 4.170074e-05  max resid 6.68812e-05 
    ## ... Similar to previous best
    ## Run 293 stress 9.776819e-05 
    ## ... Procrustes: rmse 0.0002517293  max resid 0.0004077498 
    ## ... Similar to previous best
    ## Run 294 stress 9.833033e-05 
    ## ... Procrustes: rmse 5.125489e-05  max resid 8.968777e-05 
    ## ... Similar to previous best
    ## Run 295 stress 7.888579e-05 
    ## ... Procrustes: rmse 3.904431e-05  max resid 6.711898e-05 
    ## ... Similar to previous best
    ## Run 296 stress 9.948373e-05 
    ## ... Procrustes: rmse 0.0002298868  max resid 0.0003744668 
    ## ... Similar to previous best
    ## Run 297 stress 6.367144e-05 
    ## ... Procrustes: rmse 3.124402e-05  max resid 4.952828e-05 
    ## ... Similar to previous best
    ## Run 298 stress 0.382594 
    ## Run 299 stress 9.923849e-05 
    ## ... Procrustes: rmse 0.000223174  max resid 0.0003673635 
    ## ... Similar to previous best
    ## Run 300 stress 9.866482e-05 
    ## ... Procrustes: rmse 0.0002266229  max resid 0.0003704571 
    ## ... Similar to previous best
    ## Run 301 stress 9.517965e-05 
    ## ... Procrustes: rmse 0.0001806627  max resid 0.0003062179 
    ## ... Similar to previous best
    ## Run 302 stress 9.910555e-05 
    ## ... Procrustes: rmse 4.259207e-05  max resid 9.394094e-05 
    ## ... Similar to previous best
    ## Run 303 stress 8.45854e-05 
    ## ... Procrustes: rmse 4.919124e-05  max resid 9.849449e-05 
    ## ... Similar to previous best
    ## Run 304 stress 4.672695e-05 
    ## ... Procrustes: rmse 2.438745e-05  max resid 5.067153e-05 
    ## ... Similar to previous best
    ## Run 305 stress 9.928043e-05 
    ## ... Procrustes: rmse 0.0002300688  max resid 0.0003737378 
    ## ... Similar to previous best
    ## Run 306 stress 9.464725e-05 
    ## ... Procrustes: rmse 0.0002459892  max resid 0.0003943591 
    ## ... Similar to previous best
    ## Run 307 stress 9.783947e-05 
    ## ... Procrustes: rmse 0.0002515792  max resid 0.0004065229 
    ## ... Similar to previous best
    ## Run 308 stress 8.097173e-05 
    ## ... Procrustes: rmse 0.0001213305  max resid 0.0001950263 
    ## ... Similar to previous best
    ## Run 309 stress 6.679408e-05 
    ## ... Procrustes: rmse 7.195493e-05  max resid 0.0001398609 
    ## ... Similar to previous best
    ## Run 310 stress 9.31135e-05 
    ## ... Procrustes: rmse 8.862884e-05  max resid 0.0001553209 
    ## ... Similar to previous best
    ## Run 311 stress 9.568804e-05 
    ## ... Procrustes: rmse 0.0002174857  max resid 0.0003537983 
    ## ... Similar to previous best
    ## Run 312 stress 8.152539e-05 
    ## ... Procrustes: rmse 3.606919e-05  max resid 7.105049e-05 
    ## ... Similar to previous best
    ## Run 313 stress 9.976599e-05 
    ## ... Procrustes: rmse 0.000252954  max resid 0.0004069709 
    ## ... Similar to previous best
    ## Run 314 stress 9.587349e-05 
    ## ... Procrustes: rmse 5.303545e-05  max resid 8.947578e-05 
    ## ... Similar to previous best
    ## Run 315 stress 9.640986e-05 
    ## ... Procrustes: rmse 0.0002508203  max resid 0.0003999256 
    ## ... Similar to previous best
    ## Run 316 stress 9.491392e-05 
    ## ... Procrustes: rmse 0.000244759  max resid 0.0003945709 
    ## ... Similar to previous best
    ## Run 317 stress 9.976479e-05 
    ## ... Procrustes: rmse 0.0002086088  max resid 0.0003520538 
    ## ... Similar to previous best
    ## Run 318 stress 9.796422e-05 
    ## ... Procrustes: rmse 0.0002144145  max resid 0.0003469274 
    ## ... Similar to previous best
    ## Run 319 stress 9.982136e-05 
    ## ... Procrustes: rmse 0.0001844423  max resid 0.0003082568 
    ## ... Similar to previous best
    ## Run 320 stress 9.771163e-05 
    ## ... Procrustes: rmse 0.0002464322  max resid 0.0004001142 
    ## ... Similar to previous best
    ## Run 321 stress 9.754539e-05 
    ## ... Procrustes: rmse 0.0002203992  max resid 0.0003743732 
    ## ... Similar to previous best
    ## Run 322 stress 9.908612e-05 
    ## ... Procrustes: rmse 0.0001849955  max resid 0.0003088254 
    ## ... Similar to previous best
    ## Run 323 stress 8.614361e-05 
    ## ... Procrustes: rmse 0.0001559861  max resid 0.0002642564 
    ## ... Similar to previous best
    ## Run 324 stress 9.901683e-05 
    ## ... Procrustes: rmse 0.0002407031  max resid 0.0004039078 
    ## ... Similar to previous best
    ## Run 325 stress 4.341487e-05 
    ## ... Procrustes: rmse 2.121076e-05  max resid 4.715602e-05 
    ## ... Similar to previous best
    ## Run 326 stress 9.969003e-05 
    ## ... Procrustes: rmse 0.0001954047  max resid 0.0003231307 
    ## ... Similar to previous best
    ## Run 327 stress 7.312507e-05 
    ## ... Procrustes: rmse 0.0001338661  max resid 0.000217497 
    ## ... Similar to previous best
    ## Run 328 stress 9.410805e-05 
    ## ... Procrustes: rmse 0.0001730863  max resid 0.0003010036 
    ## ... Similar to previous best
    ## Run 329 stress 9.703676e-05 
    ## ... Procrustes: rmse 0.0002139457  max resid 0.0003573511 
    ## ... Similar to previous best
    ## Run 330 stress 9.875404e-05 
    ## ... Procrustes: rmse 0.000199117  max resid 0.0003252561 
    ## ... Similar to previous best
    ## Run 331 stress 9.243784e-05 
    ## ... Procrustes: rmse 0.0001862927  max resid 0.0002987573 
    ## ... Similar to previous best
    ## Run 332 stress 9.571029e-05 
    ## ... Procrustes: rmse 0.0002294754  max resid 0.0003820262 
    ## ... Similar to previous best
    ## Run 333 stress 8.838961e-05 
    ## ... Procrustes: rmse 0.0001202822  max resid 0.0002119947 
    ## ... Similar to previous best
    ## Run 334 stress 9.801548e-05 
    ## ... Procrustes: rmse 0.0002114945  max resid 0.0003511615 
    ## ... Similar to previous best
    ## Run 335 stress 9.76262e-05 
    ## ... Procrustes: rmse 0.0002528924  max resid 0.0004098619 
    ## ... Similar to previous best
    ## Run 336 stress 9.821139e-05 
    ## ... Procrustes: rmse 0.0002543039  max resid 0.0004107944 
    ## ... Similar to previous best
    ## Run 337 stress 9.887938e-05 
    ## ... Procrustes: rmse 0.0002051729  max resid 0.0003465602 
    ## ... Similar to previous best
    ## Run 338 stress 9.690361e-05 
    ## ... Procrustes: rmse 0.0002031694  max resid 0.0003435816 
    ## ... Similar to previous best
    ## Run 339 stress 9.90111e-05 
    ## ... Procrustes: rmse 0.0002352724  max resid 0.0003904622 
    ## ... Similar to previous best
    ## Run 340 stress 6.646784e-05 
    ## ... Procrustes: rmse 2.980764e-05  max resid 5.938085e-05 
    ## ... Similar to previous best
    ## Run 341 stress 9.733825e-05 
    ## ... Procrustes: rmse 0.0002002559  max resid 0.000327123 
    ## ... Similar to previous best
    ## Run 342 stress 9.69384e-05 
    ## ... Procrustes: rmse 0.0001575228  max resid 0.000257357 
    ## ... Similar to previous best
    ## Run 343 stress 9.732546e-05 
    ## ... Procrustes: rmse 0.0002213364  max resid 0.0003657237 
    ## ... Similar to previous best
    ## Run 344 stress 9.76069e-05 
    ## ... Procrustes: rmse 0.0001940444  max resid 0.000316571 
    ## ... Similar to previous best
    ## Run 345 stress 9.827129e-05 
    ## ... Procrustes: rmse 0.0001811038  max resid 0.000300993 
    ## ... Similar to previous best
    ## Run 346 stress 9.648667e-05 
    ## ... Procrustes: rmse 0.0002341622  max resid 0.0003809356 
    ## ... Similar to previous best
    ## Run 347 stress 9.972835e-05 
    ## ... Procrustes: rmse 0.0002560831  max resid 0.000414796 
    ## ... Similar to previous best
    ## Run 348 stress 9.763375e-05 
    ## ... Procrustes: rmse 0.000218757  max resid 0.0003601337 
    ## ... Similar to previous best
    ## Run 349 stress 9.847393e-05 
    ## ... Procrustes: rmse 0.0002567891  max resid 0.0004129108 
    ## ... Similar to previous best
    ## Run 350 stress 9.317403e-05 
    ## ... Procrustes: rmse 5.289894e-05  max resid 8.766934e-05 
    ## ... Similar to previous best
    ## Run 351 stress 8.838582e-05 
    ## ... Procrustes: rmse 0.0001151964  max resid 0.0001949066 
    ## ... Similar to previous best
    ## Run 352 stress 7.24999e-05 
    ## ... Procrustes: rmse 4.823006e-05  max resid 7.709388e-05 
    ## ... Similar to previous best
    ## Run 353 stress 9.869382e-05 
    ## ... Procrustes: rmse 0.0002133337  max resid 0.000353896 
    ## ... Similar to previous best
    ## Run 354 stress 5.711172e-05 
    ## ... Procrustes: rmse 3.392865e-05  max resid 5.900194e-05 
    ## ... Similar to previous best
    ## Run 355 stress 9.933259e-05 
    ## ... Procrustes: rmse 0.0002190732  max resid 0.0003624039 
    ## ... Similar to previous best
    ## Run 356 stress 9.699839e-05 
    ## ... Procrustes: rmse 0.0002491677  max resid 0.0004061354 
    ## ... Similar to previous best
    ## Run 357 stress 9.714444e-05 
    ## ... Procrustes: rmse 0.0002179474  max resid 0.0003473043 
    ## ... Similar to previous best
    ## Run 358 stress 9.829585e-05 
    ## ... Procrustes: rmse 0.0002215906  max resid 0.0003623222 
    ## ... Similar to previous best
    ## Run 359 stress 9.275901e-05 
    ## ... Procrustes: rmse 0.0001772714  max resid 0.0003028384 
    ## ... Similar to previous best
    ## Run 360 stress 9.340962e-05 
    ## ... Procrustes: rmse 4.60654e-05  max resid 8.006728e-05 
    ## ... Similar to previous best
    ## Run 361 stress 9.749877e-05 
    ## ... Procrustes: rmse 0.0002543662  max resid 0.0004072383 
    ## ... Similar to previous best
    ## Run 362 stress 9.675204e-05 
    ## ... Procrustes: rmse 0.0002122548  max resid 0.0003564901 
    ## ... Similar to previous best
    ## Run 363 stress 7.911216e-05 
    ## ... Procrustes: rmse 4.229449e-05  max resid 7.225377e-05 
    ## ... Similar to previous best
    ## Run 364 stress 9.622809e-05 
    ## ... Procrustes: rmse 0.0002094319  max resid 0.0003492324 
    ## ... Similar to previous best
    ## Run 365 stress 9.752916e-05 
    ## ... Procrustes: rmse 0.000250712  max resid 0.0004083964 
    ## ... Similar to previous best
    ## Run 366 stress 9.925673e-05 
    ## ... Procrustes: rmse 0.0002486818  max resid 0.0004043969 
    ## ... Similar to previous best
    ## Run 367 stress 6.86072e-05 
    ## ... Procrustes: rmse 0.0001166071  max resid 0.0001900552 
    ## ... Similar to previous best
    ## Run 368 stress 9.674711e-05 
    ## ... Procrustes: rmse 0.0002197237  max resid 0.0003649718 
    ## ... Similar to previous best
    ## Run 369 stress 9.792582e-05 
    ## ... Procrustes: rmse 0.0002246411  max resid 0.0003681408 
    ## ... Similar to previous best
    ## Run 370 stress 9.90065e-05 
    ## ... Procrustes: rmse 0.0002506462  max resid 0.0004092409 
    ## ... Similar to previous best
    ## Run 371 stress 9.982113e-05 
    ## ... Procrustes: rmse 0.0002582998  max resid 0.0004148271 
    ## ... Similar to previous best
    ## Run 372 stress 9.800524e-05 
    ## ... Procrustes: rmse 0.0002322445  max resid 0.0003866085 
    ## ... Similar to previous best
    ## Run 373 stress 9.978611e-05 
    ## ... Procrustes: rmse 0.0002244596  max resid 0.0003671887 
    ## ... Similar to previous best
    ## Run 374 stress 9.850362e-05 
    ## ... Procrustes: rmse 0.0002527802  max resid 0.0004096225 
    ## ... Similar to previous best
    ## Run 375 stress 9.958174e-05 
    ## ... Procrustes: rmse 0.0002502706  max resid 0.0004066347 
    ## ... Similar to previous best
    ## Run 376 stress 4.710863e-05 
    ## ... Procrustes: rmse 2.585825e-05  max resid 5.854558e-05 
    ## ... Similar to previous best
    ## Run 377 stress 7.622116e-05 
    ## ... Procrustes: rmse 3.480526e-05  max resid 7.372974e-05 
    ## ... Similar to previous best
    ## Run 378 stress 9.978805e-05 
    ## ... Procrustes: rmse 0.0002597105  max resid 0.0004183423 
    ## ... Similar to previous best
    ## Run 379 stress 9.75305e-05 
    ## ... Procrustes: rmse 0.0002136253  max resid 0.0003493515 
    ## ... Similar to previous best
    ## Run 380 stress 9.584348e-05 
    ## ... Procrustes: rmse 0.0002044788  max resid 0.0003358194 
    ## ... Similar to previous best
    ## Run 381 stress 9.296307e-05 
    ## ... Procrustes: rmse 0.0002019892  max resid 0.0003333506 
    ## ... Similar to previous best
    ## Run 382 stress 9.980367e-05 
    ## ... Procrustes: rmse 4.976268e-05  max resid 8.790814e-05 
    ## ... Similar to previous best
    ## Run 383 stress 9.500864e-05 
    ## ... Procrustes: rmse 0.0002262342  max resid 0.0003733881 
    ## ... Similar to previous best
    ## Run 384 stress 9.743971e-05 
    ## ... Procrustes: rmse 0.0002225832  max resid 0.0003616429 
    ## ... Similar to previous best
    ## Run 385 stress 9.736813e-05 
    ## ... Procrustes: rmse 0.0002203801  max resid 0.0003658757 
    ## ... Similar to previous best
    ## Run 386 stress 8.146562e-05 
    ## ... Procrustes: rmse 0.0001536208  max resid 0.0002520457 
    ## ... Similar to previous best
    ## Run 387 stress 9.132036e-05 
    ## ... Procrustes: rmse 0.0001593354  max resid 0.0002792687 
    ## ... Similar to previous best
    ## Run 388 stress 9.815372e-05 
    ## ... Procrustes: rmse 0.0002222492  max resid 0.0003656534 
    ## ... Similar to previous best
    ## Run 389 stress 9.492632e-05 
    ## ... Procrustes: rmse 0.0002423358  max resid 0.0003936911 
    ## ... Similar to previous best
    ## Run 390 stress 9.356185e-05 
    ## ... Procrustes: rmse 4.351291e-05  max resid 8.565733e-05 
    ## ... Similar to previous best
    ## Run 391 stress 9.634663e-05 
    ## ... Procrustes: rmse 0.0002473525  max resid 0.0003976169 
    ## ... Similar to previous best
    ## Run 392 stress 9.912448e-05 
    ## ... Procrustes: rmse 0.0002225275  max resid 0.0003672544 
    ## ... Similar to previous best
    ## Run 393 stress 7.892141e-05 
    ## ... Procrustes: rmse 4.834943e-05  max resid 9.206432e-05 
    ## ... Similar to previous best
    ## Run 394 stress 7.51158e-05 
    ## ... Procrustes: rmse 5.778304e-05  max resid 0.000119311 
    ## ... Similar to previous best
    ## Run 395 stress 9.918035e-05 
    ## ... Procrustes: rmse 0.0002512753  max resid 0.0004077879 
    ## ... Similar to previous best
    ## Run 396 stress 9.937709e-05 
    ## ... Procrustes: rmse 0.0001239457  max resid 0.0002112772 
    ## ... Similar to previous best
    ## Run 397 stress 9.882719e-05 
    ## ... Procrustes: rmse 0.0002560945  max resid 0.0004143417 
    ## ... Similar to previous best
    ## Run 398 stress 9.697137e-05 
    ## ... Procrustes: rmse 0.0002132885  max resid 0.0003476037 
    ## ... Similar to previous best
    ## Run 399 stress 9.672653e-05 
    ## ... Procrustes: rmse 0.0002515564  max resid 0.0004009609 
    ## ... Similar to previous best
    ## Run 400 stress 9.939498e-05 
    ## ... Procrustes: rmse 0.000257309  max resid 0.0004121385 
    ## ... Similar to previous best
    ## Run 401 stress 9.628241e-05 
    ## ... Procrustes: rmse 0.0001964388  max resid 0.0003165123 
    ## ... Similar to previous best
    ## Run 402 stress 9.009172e-05 
    ## ... Procrustes: rmse 3.924138e-05  max resid 8.360379e-05 
    ## ... Similar to previous best
    ## Run 403 stress 9.65667e-05 
    ## ... Procrustes: rmse 0.0002178358  max resid 0.0003595101 
    ## ... Similar to previous best
    ## Run 404 stress 9.994866e-05 
    ## ... Procrustes: rmse 0.0002427924  max resid 0.0003906403 
    ## ... Similar to previous best
    ## Run 405 stress 9.407615e-05 
    ## ... Procrustes: rmse 0.0001978538  max resid 0.0003295183 
    ## ... Similar to previous best
    ## Run 406 stress 9.500405e-05 
    ## ... Procrustes: rmse 0.0002314279  max resid 0.0003801882 
    ## ... Similar to previous best
    ## Run 407 stress 9.963909e-05 
    ## ... Procrustes: rmse 0.0002574037  max resid 0.0004147368 
    ## ... Similar to previous best
    ## Run 408 stress 9.722469e-05 
    ## ... Procrustes: rmse 0.0001761344  max resid 0.000301577 
    ## ... Similar to previous best
    ## Run 409 stress 4.742527e-05 
    ## ... Procrustes: rmse 2.632451e-05  max resid 7.346564e-05 
    ## ... Similar to previous best
    ## Run 410 stress 9.963373e-05 
    ## ... Procrustes: rmse 0.0002583067  max resid 0.0004094536 
    ## ... Similar to previous best
    ## Run 411 stress 8.658475e-05 
    ## ... Procrustes: rmse 9.600794e-05  max resid 0.0001534341 
    ## ... Similar to previous best
    ## Run 412 stress 9.8155e-05 
    ## ... Procrustes: rmse 0.0002126464  max resid 0.000348424 
    ## ... Similar to previous best
    ## Run 413 stress 7.262402e-05 
    ## ... Procrustes: rmse 0.0001207511  max resid 0.0002082417 
    ## ... Similar to previous best
    ## Run 414 stress 0.3824891 
    ## Run 415 stress 6.841032e-05 
    ## ... Procrustes: rmse 3.353255e-05  max resid 5.008731e-05 
    ## ... Similar to previous best
    ## Run 416 stress 8.102252e-05 
    ## ... Procrustes: rmse 3.915108e-05  max resid 5.915503e-05 
    ## ... Similar to previous best
    ## Run 417 stress 9.795951e-05 
    ## ... Procrustes: rmse 0.0002450728  max resid 0.0003950571 
    ## ... Similar to previous best
    ## Run 418 stress 9.947152e-05 
    ## ... Procrustes: rmse 0.0002089013  max resid 0.0003409253 
    ## ... Similar to previous best
    ## Run 419 stress 8.491889e-05 
    ## ... Procrustes: rmse 4.524254e-05  max resid 8.642643e-05 
    ## ... Similar to previous best
    ## Run 420 stress 9.730541e-05 
    ## ... Procrustes: rmse 0.0002027458  max resid 0.0003371324 
    ## ... Similar to previous best
    ## Run 421 stress 9.73293e-05 
    ## ... Procrustes: rmse 0.0001950714  max resid 0.0003259026 
    ## ... Similar to previous best
    ## Run 422 stress 9.336299e-05 
    ## ... Procrustes: rmse 0.0001937487  max resid 0.000328416 
    ## ... Similar to previous best
    ## Run 423 stress 9.94104e-05 
    ## ... Procrustes: rmse 0.0002165283  max resid 0.0003546285 
    ## ... Similar to previous best
    ## Run 424 stress 9.825454e-05 
    ## ... Procrustes: rmse 0.0002541821  max resid 0.0004117377 
    ## ... Similar to previous best
    ## Run 425 stress 9.699466e-05 
    ## ... Procrustes: rmse 0.0002029954  max resid 0.0003351585 
    ## ... Similar to previous best
    ## Run 426 stress 9.864303e-05 
    ## ... Procrustes: rmse 0.0002367441  max resid 0.0003937654 
    ## ... Similar to previous best
    ## Run 427 stress 9.515313e-05 
    ## ... Procrustes: rmse 0.0001606118  max resid 0.0002668104 
    ## ... Similar to previous best
    ## Run 428 stress 9.862335e-05 
    ## ... Procrustes: rmse 0.0002274824  max resid 0.0003681312 
    ## ... Similar to previous best
    ## Run 429 stress 8.76934e-05 
    ## ... Procrustes: rmse 4.277893e-05  max resid 7.996635e-05 
    ## ... Similar to previous best
    ## Run 430 stress 9.752315e-05 
    ## ... Procrustes: rmse 5.110308e-05  max resid 8.117888e-05 
    ## ... Similar to previous best
    ## Run 431 stress 9.827203e-05 
    ## ... Procrustes: rmse 0.0002243749  max resid 0.0003651475 
    ## ... Similar to previous best
    ## Run 432 stress 9.805033e-05 
    ## ... Procrustes: rmse 0.0002232338  max resid 0.0003626176 
    ## ... Similar to previous best
    ## Run 433 stress 9.658913e-05 
    ## ... Procrustes: rmse 0.0002519756  max resid 0.0004041235 
    ## ... Similar to previous best
    ## Run 434 stress 9.883694e-05 
    ## ... Procrustes: rmse 0.0002241249  max resid 0.0003707687 
    ## ... Similar to previous best
    ## Run 435 stress 9.840651e-05 
    ## ... Procrustes: rmse 0.000256294  max resid 0.0004133056 
    ## ... Similar to previous best
    ## Run 436 stress 9.687145e-05 
    ## ... Procrustes: rmse 0.0002472474  max resid 0.00040251 
    ## ... Similar to previous best
    ## Run 437 stress 8.821423e-05 
    ## ... Procrustes: rmse 0.0001464763  max resid 0.0002339673 
    ## ... Similar to previous best
    ## Run 438 stress 9.130117e-05 
    ## ... Procrustes: rmse 0.0001179213  max resid 0.0002002315 
    ## ... Similar to previous best
    ## Run 439 stress 5.899897e-05 
    ## ... Procrustes: rmse 6.470642e-05  max resid 0.0001143712 
    ## ... Similar to previous best
    ## Run 440 stress 8.687935e-05 
    ## ... Procrustes: rmse 8.067796e-05  max resid 0.0001584455 
    ## ... Similar to previous best
    ## Run 441 stress 9.337621e-05 
    ## ... Procrustes: rmse 0.0001066214  max resid 0.0001809562 
    ## ... Similar to previous best
    ## Run 442 stress 9.541174e-05 
    ## ... Procrustes: rmse 0.0002345079  max resid 0.0003695206 
    ## ... Similar to previous best
    ## Run 443 stress 9.909893e-05 
    ## ... Procrustes: rmse 0.0002564807  max resid 0.0004145111 
    ## ... Similar to previous best
    ## Run 444 stress 9.920543e-05 
    ## ... Procrustes: rmse 0.0002581491  max resid 0.0004166184 
    ## ... Similar to previous best
    ## Run 445 stress 9.716965e-05 
    ## ... Procrustes: rmse 0.000206807  max resid 0.0003355293 
    ## ... Similar to previous best
    ## Run 446 stress 9.600521e-05 
    ## ... Procrustes: rmse 0.0002463185  max resid 0.0003994324 
    ## ... Similar to previous best
    ## Run 447 stress 0.3793354 
    ## Run 448 stress 9.655641e-05 
    ## ... Procrustes: rmse 0.0002467522  max resid 0.000399437 
    ## ... Similar to previous best
    ## Run 449 stress 9.783352e-05 
    ## ... Procrustes: rmse 0.0002253807  max resid 0.0003804571 
    ## ... Similar to previous best
    ## Run 450 stress 7.07802e-05 
    ## ... Procrustes: rmse 0.000108349  max resid 0.0001892398 
    ## ... Similar to previous best
    ## Run 451 stress 9.555021e-05 
    ## ... Procrustes: rmse 0.0002140755  max resid 0.0003527303 
    ## ... Similar to previous best
    ## Run 452 stress 9.564291e-05 
    ## ... Procrustes: rmse 0.0001522317  max resid 0.0002572377 
    ## ... Similar to previous best
    ## Run 453 stress 9.827286e-05 
    ## ... Procrustes: rmse 0.0002479684  max resid 0.0004041039 
    ## ... Similar to previous best
    ## Run 454 stress 5.780819e-05 
    ## ... Procrustes: rmse 3.041633e-05  max resid 5.058904e-05 
    ## ... Similar to previous best
    ## Run 455 stress 9.765445e-05 
    ## ... Procrustes: rmse 0.0002260548  max resid 0.0003808101 
    ## ... Similar to previous best
    ## Run 456 stress 9.807647e-05 
    ## ... Procrustes: rmse 0.0002143801  max resid 0.0003478103 
    ## ... Similar to previous best
    ## Run 457 stress 9.874549e-05 
    ## ... Procrustes: rmse 0.0002192732  max resid 0.000360792 
    ## ... Similar to previous best
    ## Run 458 stress 9.962746e-05 
    ## ... Procrustes: rmse 0.0002241735  max resid 0.0003644613 
    ## ... Similar to previous best
    ## Run 459 stress 7.443358e-05 
    ## ... Procrustes: rmse 0.0001070857  max resid 0.0001818357 
    ## ... Similar to previous best
    ## Run 460 stress 9.90428e-05 
    ## ... Procrustes: rmse 0.0002504576  max resid 0.0004077756 
    ## ... Similar to previous best
    ## Run 461 stress 8.291957e-05 
    ## ... Procrustes: rmse 5.132569e-05  max resid 9.155757e-05 
    ## ... Similar to previous best
    ## Run 462 stress 6.620611e-05 
    ## ... Procrustes: rmse 3.497648e-05  max resid 6.630692e-05 
    ## ... Similar to previous best
    ## Run 463 stress 9.812221e-05 
    ## ... Procrustes: rmse 0.0002205535  max resid 0.0003594751 
    ## ... Similar to previous best
    ## Run 464 stress 9.522599e-05 
    ## ... Procrustes: rmse 0.0002315695  max resid 0.0003927391 
    ## ... Similar to previous best
    ## Run 465 stress 8.954969e-05 
    ## ... Procrustes: rmse 0.000189909  max resid 0.0003200294 
    ## ... Similar to previous best
    ## Run 466 stress 9.705652e-05 
    ## ... Procrustes: rmse 0.0002496519  max resid 0.0004058727 
    ## ... Similar to previous best
    ## Run 467 stress 9.980044e-05 
    ## ... Procrustes: rmse 0.0002526988  max resid 0.0004090804 
    ## ... Similar to previous best
    ## Run 468 stress 9.802817e-05 
    ## ... Procrustes: rmse 0.0002059108  max resid 0.0003493187 
    ## ... Similar to previous best
    ## Run 469 stress 8.402564e-05 
    ## ... Procrustes: rmse 4.501633e-05  max resid 7.761513e-05 
    ## ... Similar to previous best
    ## Run 470 stress 9.323408e-05 
    ## ... Procrustes: rmse 8.005635e-05  max resid 0.0001366249 
    ## ... Similar to previous best
    ## Run 471 stress 9.840884e-05 
    ## ... Procrustes: rmse 0.000250929  max resid 0.0004076891 
    ## ... Similar to previous best
    ## Run 472 stress 9.799315e-05 
    ## ... Procrustes: rmse 0.0002516576  max resid 0.0004044604 
    ## ... Similar to previous best
    ## Run 473 stress 9.953532e-05 
    ## ... Procrustes: rmse 0.000193759  max resid 0.0003249045 
    ## ... Similar to previous best
    ## Run 474 stress 9.789878e-05 
    ## ... Procrustes: rmse 0.0002496738  max resid 0.0004084815 
    ## ... Similar to previous best
    ## Run 475 stress 9.448082e-05 
    ## ... Procrustes: rmse 0.0002185709  max resid 0.0003548666 
    ## ... Similar to previous best
    ## Run 476 stress 9.763815e-05 
    ## ... Procrustes: rmse 0.0001787965  max resid 0.0002943593 
    ## ... Similar to previous best
    ## Run 477 stress 7.342411e-05 
    ## ... Procrustes: rmse 8.044212e-05  max resid 0.0001367167 
    ## ... Similar to previous best
    ## Run 478 stress 9.865867e-05 
    ## ... Procrustes: rmse 0.0001972812  max resid 0.0003280082 
    ## ... Similar to previous best
    ## Run 479 stress 9.785046e-05 
    ## ... Procrustes: rmse 0.0002331088  max resid 0.0003895808 
    ## ... Similar to previous best
    ## Run 480 stress 9.895174e-05 
    ## ... Procrustes: rmse 0.0002492654  max resid 0.0004019441 
    ## ... Similar to previous best
    ## Run 481 stress 9.03581e-05 
    ## ... Procrustes: rmse 4.829648e-05  max resid 8.196962e-05 
    ## ... Similar to previous best
    ## Run 482 stress 9.819323e-05 
    ## ... Procrustes: rmse 0.0001993931  max resid 0.0003303882 
    ## ... Similar to previous best
    ## Run 483 stress 9.734135e-05 
    ## ... Procrustes: rmse 0.0002282638  max resid 0.0003839941 
    ## ... Similar to previous best
    ## Run 484 stress 9.677864e-05 
    ## ... Procrustes: rmse 0.0001877151  max resid 0.0003044897 
    ## ... Similar to previous best
    ## Run 485 stress 9.989192e-05 
    ## ... Procrustes: rmse 0.0002127495  max resid 0.0003545693 
    ## ... Similar to previous best
    ## Run 486 stress 9.913994e-05 
    ## ... Procrustes: rmse 0.0002539027  max resid 0.0004137051 
    ## ... Similar to previous best
    ## Run 487 stress 9.347196e-05 
    ## ... Procrustes: rmse 0.000238965  max resid 0.000384147 
    ## ... Similar to previous best
    ## Run 488 stress 9.980424e-05 
    ## ... Procrustes: rmse 0.0002580467  max resid 0.0004150564 
    ## ... Similar to previous best
    ## Run 489 stress 9.974521e-05 
    ## ... Procrustes: rmse 0.0002577797  max resid 0.0004164755 
    ## ... Similar to previous best
    ## Run 490 stress 9.730971e-05 
    ## ... Procrustes: rmse 0.0002182915  max resid 0.0003588431 
    ## ... Similar to previous best
    ## Run 491 stress 9.718783e-05 
    ## ... Procrustes: rmse 0.0002518349  max resid 0.0004088262 
    ## ... Similar to previous best
    ## Run 492 stress 9.511715e-05 
    ## ... Procrustes: rmse 0.0001759555  max resid 0.0002932729 
    ## ... Similar to previous best
    ## Run 493 stress 9.973856e-05 
    ## ... Procrustes: rmse 0.0002588958  max resid 0.0004157649 
    ## ... Similar to previous best
    ## Run 494 stress 7.032393e-05 
    ## ... Procrustes: rmse 3.689417e-05  max resid 5.626976e-05 
    ## ... Similar to previous best
    ## Run 495 stress 9.928091e-05 
    ## ... Procrustes: rmse 0.0002278864  max resid 0.0003699359 
    ## ... Similar to previous best
    ## Run 496 stress 8.800003e-05 
    ## ... Procrustes: rmse 4.306021e-05  max resid 8.999927e-05 
    ## ... Similar to previous best
    ## Run 497 stress 9.644667e-05 
    ## ... Procrustes: rmse 7.13641e-05  max resid 0.0001409266 
    ## ... Similar to previous best
    ## Run 498 stress 9.79136e-05 
    ## ... Procrustes: rmse 0.000202124  max resid 0.0003262237 
    ## ... Similar to previous best
    ## Run 499 stress 9.663308e-05 
    ## ... Procrustes: rmse 0.0002509083  max resid 0.0004032565 
    ## ... Similar to previous best
    ## Run 500 stress 9.895529e-05 
    ## ... Procrustes: rmse 0.000211028  max resid 0.0003478031 
    ## ... Similar to previous best
    ## *** Best solution repeated 214 times

    ## Warning in metaMDS(PD_beta_ref_dist$Btotal, try = 1000, parallel = 4, trymax =
    ## 500, : stress is (nearly) zero: you may have insufficient data

``` r
# Surveyed sites 
PD_beta_NMDS <- metaMDS(PD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07927535 
    ## Run 1 stress 0.09290716 
    ## Run 2 stress 0.0796926 
    ## ... Procrustes: rmse 0.01235944  max resid 0.04963438 
    ## Run 3 stress 0.09106027 
    ## Run 4 stress 0.0795144 
    ## ... Procrustes: rmse 0.01733858  max resid 0.05637792 
    ## Run 5 stress 0.09247102 
    ## Run 6 stress 0.07936951 
    ## ... Procrustes: rmse 0.009031072  max resid 0.03238542 
    ## Run 7 stress 0.09072897 
    ## Run 8 stress 0.07936951 
    ## ... Procrustes: rmse 0.009033295  max resid 0.03239898 
    ## Run 9 stress 0.09290693 
    ## Run 10 stress 0.08985606 
    ## Run 11 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733028  max resid 0.05632434 
    ## Run 12 stress 0.08985603 
    ## Run 13 stress 0.07969258 
    ## ... Procrustes: rmse 0.01236611  max resid 0.0496941 
    ## Run 14 stress 0.08985607 
    ## Run 15 stress 0.09252226 
    ## Run 16 stress 0.08985605 
    ## Run 17 stress 0.0938206 
    ## Run 18 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 2.169723e-05  max resid 5.307107e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.07927536 
    ## ... Procrustes: rmse 3.789644e-05  max resid 8.596156e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.0938836 
    ## Run 21 stress 0.07927535 
    ## ... Procrustes: rmse 3.266754e-05  max resid 8.788966e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734402  max resid 0.05637668 
    ## Run 23 stress 0.07927536 
    ## ... Procrustes: rmse 2.679626e-05  max resid 6.034761e-05 
    ## ... Similar to previous best
    ## Run 24 stress 0.09544397 
    ## Run 25 stress 0.07936954 
    ## ... Procrustes: rmse 0.009026598  max resid 0.03235175 
    ## Run 26 stress 0.09106032 
    ## Run 27 stress 0.08985614 
    ## Run 28 stress 0.08985608 
    ## Run 29 stress 0.09544393 
    ## Run 30 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237238  max resid 0.04971676 
    ## Run 31 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240219  max resid 0.04987646 
    ## Run 32 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237826  max resid 0.04973477 
    ## Run 33 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735672  max resid 0.05649109 
    ## Run 34 stress 0.0904376 
    ## Run 35 stress 0.07936953 
    ## ... Procrustes: rmse 0.00902888  max resid 0.03236515 
    ## Run 36 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734773  max resid 0.05639929 
    ## Run 37 stress 0.07927535 
    ## ... Procrustes: rmse 3.246901e-05  max resid 7.713916e-05 
    ## ... Similar to previous best
    ## Run 38 stress 0.07936954 
    ## ... Procrustes: rmse 0.009029941  max resid 0.03237211 
    ## Run 39 stress 0.07927534 
    ## ... Procrustes: rmse 1.392356e-05  max resid 3.18197e-05 
    ## ... Similar to previous best
    ## Run 40 stress 0.09043771 
    ## Run 41 stress 0.07936956 
    ## ... Procrustes: rmse 0.009022461  max resid 0.03232548 
    ## Run 42 stress 0.07927551 
    ## ... Procrustes: rmse 6.435041e-05  max resid 0.0001647472 
    ## ... Similar to previous best
    ## Run 43 stress 0.09296252 
    ## Run 44 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734589  max resid 0.05635991 
    ## Run 45 stress 0.09043777 
    ## Run 46 stress 0.09106031 
    ## Run 47 stress 0.0902278 
    ## Run 48 stress 0.07936949 
    ## ... Procrustes: rmse 0.009044727  max resid 0.03246719 
    ## Run 49 stress 0.1033494 
    ## Run 50 stress 0.09344088 
    ## Run 51 stress 0.07936949 
    ## ... Procrustes: rmse 0.009044809  max resid 0.03247392 
    ## Run 52 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041243  max resid 0.03244403 
    ## Run 53 stress 0.09043749 
    ## Run 54 stress 0.07969268 
    ## ... Procrustes: rmse 0.01235974  max resid 0.04952268 
    ## Run 55 stress 0.08985604 
    ## Run 56 stress 0.08985606 
    ## Run 57 stress 0.07936951 
    ## ... Procrustes: rmse 0.00904444  max resid 0.03246568 
    ## Run 58 stress 0.08985607 
    ## Run 59 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735265  max resid 0.05643234 
    ## Run 60 stress 0.09106039 
    ## Run 61 stress 0.09565367 
    ## Run 62 stress 0.09043776 
    ## Run 63 stress 0.07927534 
    ## ... Procrustes: rmse 1.477421e-05  max resid 3.901e-05 
    ## ... Similar to previous best
    ## Run 64 stress 0.09252214 
    ## Run 65 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733225  max resid 0.05628602 
    ## Run 66 stress 0.08985603 
    ## Run 67 stress 0.09106042 
    ## Run 68 stress 0.07927537 
    ## ... Procrustes: rmse 7.416181e-05  max resid 0.0002081285 
    ## ... Similar to previous best
    ## Run 69 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045964  max resid 0.03247328 
    ## Run 70 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732104  max resid 0.05625997 
    ## Run 71 stress 0.09072898 
    ## Run 72 stress 0.09382094 
    ## Run 73 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049262  max resid 0.03249119 
    ## Run 74 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239715  max resid 0.04985295 
    ## Run 75 stress 0.09485701 
    ## Run 76 stress 0.1033494 
    ## Run 77 stress 0.09022778 
    ## Run 78 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 1.111528e-05  max resid 2.753847e-05 
    ## ... Similar to previous best
    ## Run 79 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 4.702055e-06  max resid 1.308928e-05 
    ## ... Similar to previous best
    ## Run 80 stress 0.07927535 
    ## ... Procrustes: rmse 3.378938e-05  max resid 7.74486e-05 
    ## ... Similar to previous best
    ## Run 81 stress 0.0795144 
    ## ... Procrustes: rmse 0.01736038  max resid 0.05645521 
    ## Run 82 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238706  max resid 0.04976859 
    ## Run 83 stress 0.08985616 
    ## Run 84 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734949  max resid 0.05639443 
    ## Run 85 stress 0.07927535 
    ## ... Procrustes: rmse 1.961823e-05  max resid 4.723458e-05 
    ## ... Similar to previous best
    ## Run 86 stress 0.09485683 
    ## Run 87 stress 0.09296189 
    ## Run 88 stress 0.07927534 
    ## ... Procrustes: rmse 1.917181e-05  max resid 5.485996e-05 
    ## ... Similar to previous best
    ## Run 89 stress 0.09252214 
    ## Run 90 stress 0.1052402 
    ## Run 91 stress 0.09072895 
    ## Run 92 stress 0.07927534 
    ## ... Procrustes: rmse 6.08071e-06  max resid 1.403294e-05 
    ## ... Similar to previous best
    ## Run 93 stress 0.09296214 
    ## Run 94 stress 0.1067421 
    ## Run 95 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736359  max resid 0.05649263 
    ## Run 96 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734584  max resid 0.05637814 
    ## Run 97 stress 0.108875 
    ## Run 98 stress 0.09382091 
    ## Run 99 stress 0.09043774 
    ## Run 100 stress 0.0947407 
    ## Run 101 stress 0.09344088 
    ## Run 102 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237775  max resid 0.04971199 
    ## Run 103 stress 0.07936951 
    ## ... Procrustes: rmse 0.00904509  max resid 0.03246517 
    ## Run 104 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239578  max resid 0.04979723 
    ## Run 105 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735136  max resid 0.05641154 
    ## Run 106 stress 0.09344088 
    ## Run 107 stress 0.08985612 
    ## Run 108 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237305  max resid 0.04966297 
    ## Run 109 stress 0.0796927 
    ## ... Procrustes: rmse 0.01236475  max resid 0.04951728 
    ## Run 110 stress 0.07927537 
    ## ... Procrustes: rmse 4.185483e-05  max resid 0.0001101101 
    ## ... Similar to previous best
    ## Run 111 stress 0.09043748 
    ## Run 112 stress 0.07927535 
    ## ... Procrustes: rmse 3.488426e-05  max resid 9.152656e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039855  max resid 0.03243303 
    ## Run 114 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237724  max resid 0.049709 
    ## Run 115 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735841  max resid 0.05649574 
    ## Run 116 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238981  max resid 0.04979305 
    ## Run 117 stress 0.09172781 
    ## Run 118 stress 0.08985607 
    ## Run 119 stress 0.07927534 
    ## ... Procrustes: rmse 1.20755e-05  max resid 2.95436e-05 
    ## ... Similar to previous best
    ## Run 120 stress 0.07927536 
    ## ... Procrustes: rmse 5.832059e-05  max resid 0.0001613527 
    ## ... Similar to previous best
    ## Run 121 stress 0.09043779 
    ## Run 122 stress 0.07927535 
    ## ... Procrustes: rmse 2.7679e-05  max resid 7.684004e-05 
    ## ... Similar to previous best
    ## Run 123 stress 0.09022783 
    ## Run 124 stress 0.07951445 
    ## ... Procrustes: rmse 0.01737325  max resid 0.0565578 
    ## Run 125 stress 0.1039004 
    ## Run 126 stress 0.07927534 
    ## ... Procrustes: rmse 1.305901e-05  max resid 3.132983e-05 
    ## ... Similar to previous best
    ## Run 127 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237477  max resid 0.04968192 
    ## Run 128 stress 0.07927534 
    ## ... Procrustes: rmse 7.979604e-06  max resid 1.808918e-05 
    ## ... Similar to previous best
    ## Run 129 stress 0.09022796 
    ## Run 130 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041701  max resid 0.0324357 
    ## Run 131 stress 0.09382055 
    ## Run 132 stress 0.07951446 
    ## ... Procrustes: rmse 0.01737494  max resid 0.05659714 
    ## Run 133 stress 0.09106047 
    ## Run 134 stress 0.07951448 
    ## ... Procrustes: rmse 0.01740096  max resid 0.05671865 
    ## Run 135 stress 0.09290721 
    ## Run 136 stress 0.08985619 
    ## Run 137 stress 0.07936952 
    ## ... Procrustes: rmse 0.009031963  max resid 0.03238253 
    ## Run 138 stress 0.0796926 
    ## ... Procrustes: rmse 0.01239642  max resid 0.04982782 
    ## Run 139 stress 0.09106035 
    ## Run 140 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240477  max resid 0.04988147 
    ## Run 141 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045208  max resid 0.03246976 
    ## Run 142 stress 0.07951442 
    ## ... Procrustes: rmse 0.01734861  max resid 0.05638152 
    ## Run 143 stress 0.09344086 
    ## Run 144 stress 0.09072909 
    ## Run 145 stress 0.0796926 
    ## ... Procrustes: rmse 0.01236977  max resid 0.04961838 
    ## Run 146 stress 0.09072898 
    ## Run 147 stress 0.0907291 
    ## Run 148 stress 0.0796927 
    ## ... Procrustes: rmse 0.01235756  max resid 0.04948575 
    ## Run 149 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237401  max resid 0.04970055 
    ## Run 150 stress 0.09296203 
    ## Run 151 stress 0.08985615 
    ## Run 152 stress 0.08985603 
    ## Run 153 stress 0.09106027 
    ## Run 154 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039976  max resid 0.03243366 
    ## Run 155 stress 0.1069268 
    ## Run 156 stress 0.08985603 
    ## Run 157 stress 0.09388352 
    ## Run 158 stress 0.07927534 
    ## ... Procrustes: rmse 1.274315e-05  max resid 2.707189e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.08985604 
    ## Run 160 stress 0.07936953 
    ## ... Procrustes: rmse 0.009055364  max resid 0.03251372 
    ## Run 161 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735339  max resid 0.05643587 
    ## Run 162 stress 0.07927535 
    ## ... Procrustes: rmse 1.131164e-05  max resid 2.560888e-05 
    ## ... Similar to previous best
    ## Run 163 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237221  max resid 0.04965253 
    ## Run 164 stress 0.1039004 
    ## Run 165 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047412  max resid 0.03248289 
    ## Run 166 stress 0.07927535 
    ## ... Procrustes: rmse 3.724622e-05  max resid 9.039398e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047722  max resid 0.03248305 
    ## Run 168 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237416  max resid 0.0496872 
    ## Run 169 stress 0.09565379 
    ## Run 170 stress 0.09247113 
    ## Run 171 stress 0.07936952 
    ## ... Procrustes: rmse 0.009040544  max resid 0.03245058 
    ## Run 172 stress 0.07969268 
    ## ... Procrustes: rmse 0.01242602  max resid 0.05001986 
    ## Run 173 stress 0.07969262 
    ## ... Procrustes: rmse 0.01237122  max resid 0.04960689 
    ## Run 174 stress 0.09072892 
    ## Run 175 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043736  max resid 0.03246172 
    ## Run 176 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732084  max resid 0.05624022 
    ## Run 177 stress 0.09388352 
    ## Run 178 stress 0.09072893 
    ## Run 179 stress 0.07969262 
    ## ... Procrustes: rmse 0.01240755  max resid 0.04989794 
    ## Run 180 stress 0.07936955 
    ## ... Procrustes: rmse 0.009023483  max resid 0.0323307 
    ## Run 181 stress 0.07927534 
    ## ... Procrustes: rmse 3.304958e-06  max resid 9.733208e-06 
    ## ... Similar to previous best
    ## Run 182 stress 0.0796926 
    ## ... Procrustes: rmse 0.01239132  max resid 0.04979784 
    ## Run 183 stress 0.0910603 
    ## Run 184 stress 0.08985604 
    ## Run 185 stress 0.08985611 
    ## Run 186 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904943  max resid 0.03248547 
    ## Run 187 stress 0.07927534 
    ## ... Procrustes: rmse 1.616319e-05  max resid 4.442432e-05 
    ## ... Similar to previous best
    ## Run 188 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735572  max resid 0.05647663 
    ## Run 189 stress 0.09296193 
    ## Run 190 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173668  max resid 0.05650312 
    ## Run 191 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047471  max resid 0.03247581 
    ## Run 192 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237227  max resid 0.04967216 
    ## Run 193 stress 0.1077077 
    ## Run 194 stress 0.09106036 
    ## Run 195 stress 0.07936954 
    ## ... Procrustes: rmse 0.009027416  max resid 0.03235427 
    ## Run 196 stress 0.07951442 
    ## ... Procrustes: rmse 0.01736673  max resid 0.05651194 
    ## Run 197 stress 0.08985606 
    ## Run 198 stress 0.08985603 
    ## Run 199 stress 0.09043758 
    ## Run 200 stress 0.07951443 
    ## ... Procrustes: rmse 0.01735277  max resid 0.05651451 
    ## Run 201 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047542  max resid 0.03248717 
    ## Run 202 stress 0.107708 
    ## Run 203 stress 0.1105179 
    ## Run 204 stress 0.09106028 
    ## Run 205 stress 0.1110283 
    ## Run 206 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735498  max resid 0.05645847 
    ## Run 207 stress 0.07969265 
    ## ... Procrustes: rmse 0.01241895  max resid 0.04996194 
    ## Run 208 stress 0.07936955 
    ## ... Procrustes: rmse 0.009025517  max resid 0.03234342 
    ## Run 209 stress 0.07927535 
    ## ... Procrustes: rmse 4.486407e-05  max resid 0.0001067225 
    ## ... Similar to previous best
    ## Run 210 stress 0.09565393 
    ## Run 211 stress 0.07969265 
    ## ... Procrustes: rmse 0.01240949  max resid 0.04990934 
    ## Run 212 stress 0.09388338 
    ## Run 213 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237533  max resid 0.04970241 
    ## Run 214 stress 0.09474061 
    ## Run 215 stress 0.09544387 
    ## Run 216 stress 0.09043758 
    ## Run 217 stress 0.07927534 
    ## ... Procrustes: rmse 6.467138e-06  max resid 1.466652e-05 
    ## ... Similar to previous best
    ## Run 218 stress 0.07927534 
    ## ... Procrustes: rmse 2.163027e-05  max resid 5.622564e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.1033494 
    ## Run 220 stress 0.09043736 
    ## Run 221 stress 0.09022789 
    ## Run 222 stress 0.09043756 
    ## Run 223 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237677  max resid 0.04973344 
    ## Run 224 stress 0.09072891 
    ## Run 225 stress 0.08985612 
    ## Run 226 stress 0.08985608 
    ## Run 227 stress 0.07936952 
    ## ... Procrustes: rmse 0.009033974  max resid 0.03239612 
    ## Run 228 stress 0.09106028 
    ## Run 229 stress 0.104507 
    ## Run 230 stress 0.07951442 
    ## ... Procrustes: rmse 0.01734368  max resid 0.05641019 
    ## Run 231 stress 0.08985606 
    ## Run 232 stress 0.0938834 
    ## Run 233 stress 0.09043766 
    ## Run 234 stress 0.104507 
    ## Run 235 stress 0.09072892 
    ## Run 236 stress 0.09485683 
    ## Run 237 stress 0.08985612 
    ## Run 238 stress 0.07927535 
    ## ... Procrustes: rmse 4.411312e-05  max resid 0.0001160829 
    ## ... Similar to previous best
    ## Run 239 stress 0.08985619 
    ## Run 240 stress 0.07927534 
    ## ... Procrustes: rmse 1.845499e-05  max resid 4.754803e-05 
    ## ... Similar to previous best
    ## Run 241 stress 0.07936955 
    ## ... Procrustes: rmse 0.009024524  max resid 0.03233672 
    ## Run 242 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735165  max resid 0.05641202 
    ## Run 243 stress 0.07969261 
    ## ... Procrustes: rmse 0.01238309  max resid 0.04968019 
    ## Run 244 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238487  max resid 0.04977727 
    ## Run 245 stress 0.08985603 
    ## Run 246 stress 0.09290699 
    ## Run 247 stress 0.07969263 
    ## ... Procrustes: rmse 0.01238489  max resid 0.04975466 
    ## Run 248 stress 0.07927537 
    ## ... Procrustes: rmse 6.438024e-05  max resid 0.0001528284 
    ## ... Similar to previous best
    ## Run 249 stress 0.09344088 
    ## Run 250 stress 0.08985607 
    ## Run 251 stress 0.07927534 
    ## ... Procrustes: rmse 6.859366e-06  max resid 1.885776e-05 
    ## ... Similar to previous best
    ## Run 252 stress 0.09106047 
    ## Run 253 stress 0.09043755 
    ## Run 254 stress 0.09172783 
    ## Run 255 stress 0.07927536 
    ## ... Procrustes: rmse 4.01921e-05  max resid 0.0001049652 
    ## ... Similar to previous best
    ## Run 256 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735611  max resid 0.05645735 
    ## Run 257 stress 0.07969266 
    ## ... Procrustes: rmse 0.012369  max resid 0.04957828 
    ## Run 258 stress 0.09022785 
    ## Run 259 stress 0.09290698 
    ## Run 260 stress 0.09106039 
    ## Run 261 stress 0.09252227 
    ## Run 262 stress 0.1067418 
    ## Run 263 stress 0.09043749 
    ## Run 264 stress 0.07969258 
    ## ... Procrustes: rmse 0.0123728  max resid 0.04968612 
    ## Run 265 stress 0.07927534 
    ## ... Procrustes: rmse 8.717416e-06  max resid 2.021801e-05 
    ## ... Similar to previous best
    ## Run 266 stress 0.08985611 
    ## Run 267 stress 0.08985606 
    ## Run 268 stress 0.07951443 
    ## ... Procrustes: rmse 0.01735393  max resid 0.05651521 
    ## Run 269 stress 0.09043747 
    ## Run 270 stress 0.09252216 
    ## Run 271 stress 0.07951442 
    ## ... Procrustes: rmse 0.017331  max resid 0.05627011 
    ## Run 272 stress 0.09043761 
    ## Run 273 stress 0.1056806 
    ## Run 274 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735465  max resid 0.05650578 
    ## Run 275 stress 0.07927535 
    ## ... Procrustes: rmse 4.275908e-05  max resid 9.304567e-05 
    ## ... Similar to previous best
    ## Run 276 stress 0.09172797 
    ## Run 277 stress 0.09296251 
    ## Run 278 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733099  max resid 0.05627232 
    ## Run 279 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734876  max resid 0.05639308 
    ## Run 280 stress 0.09474049 
    ## Run 281 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735341  max resid 0.05643192 
    ## Run 282 stress 0.07951444 
    ## ... Procrustes: rmse 0.0173606  max resid 0.05651716 
    ## Run 283 stress 0.07927536 
    ## ... Procrustes: rmse 2.683636e-05  max resid 6.299209e-05 
    ## ... Similar to previous best
    ## Run 284 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904191  max resid 0.03244603 
    ## Run 285 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732431  max resid 0.05626696 
    ## Run 286 stress 0.08985604 
    ## Run 287 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237033  max resid 0.04961143 
    ## Run 288 stress 0.09544357 
    ## Run 289 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733999  max resid 0.05632646 
    ## Run 290 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736199  max resid 0.05647793 
    ## Run 291 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735653  max resid 0.05647904 
    ## Run 292 stress 0.08985612 
    ## Run 293 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734902  max resid 0.05637029 
    ## Run 294 stress 0.07936954 
    ## ... Procrustes: rmse 0.009029726  max resid 0.03236988 
    ## Run 295 stress 0.08985603 
    ## Run 296 stress 0.09022785 
    ## Run 297 stress 0.07951442 
    ## ... Procrustes: rmse 0.01732456  max resid 0.05623462 
    ## Run 298 stress 0.08985603 
    ## Run 299 stress 0.09382101 
    ## Run 300 stress 0.07927535 
    ## ... Procrustes: rmse 3.117924e-05  max resid 7.3681e-05 
    ## ... Similar to previous best
    ## Run 301 stress 0.08985606 
    ## Run 302 stress 0.07927538 
    ## ... Procrustes: rmse 7.332549e-05  max resid 0.0002031056 
    ## ... Similar to previous best
    ## Run 303 stress 0.09252227 
    ## Run 304 stress 0.07936951 
    ## ... Procrustes: rmse 0.009053462  max resid 0.03250021 
    ## Run 305 stress 0.0793695 
    ## ... Procrustes: rmse 0.009042793  max resid 0.03245354 
    ## Run 306 stress 0.07927536 
    ## ... Procrustes: rmse 4.062748e-05  max resid 0.0001087245 
    ## ... Similar to previous best
    ## Run 307 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735465  max resid 0.0564323 
    ## Run 308 stress 0.0793695 
    ## ... Procrustes: rmse 0.009042584  max resid 0.03245075 
    ## Run 309 stress 0.07936952 
    ## ... Procrustes: rmse 0.009050714  max resid 0.03249424 
    ## Run 310 stress 0.07951448 
    ## ... Procrustes: rmse 0.0173336  max resid 0.05628646 
    ## Run 311 stress 0.07936949 
    ## ... Procrustes: rmse 0.009044358  max resid 0.03246244 
    ## Run 312 stress 0.09172797 
    ## Run 313 stress 0.0934409 
    ## Run 314 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041961  max resid 0.03244491 
    ## Run 315 stress 0.1061856 
    ## Run 316 stress 0.0793695 
    ## ... Procrustes: rmse 0.009042374  max resid 0.03246599 
    ## Run 317 stress 0.07936957 
    ## ... Procrustes: rmse 0.009021652  max resid 0.03231948 
    ## Run 318 stress 0.09356313 
    ## Run 319 stress 0.090729 
    ## Run 320 stress 0.07927535 
    ## ... Procrustes: rmse 3.210922e-05  max resid 8.827651e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.07936949 
    ## ... Procrustes: rmse 0.00904778  max resid 0.03248373 
    ## Run 322 stress 0.09172795 
    ## Run 323 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733355  max resid 0.05629235 
    ## Run 324 stress 0.09382085 
    ## Run 325 stress 0.0792754 
    ## ... Procrustes: rmse 9.416834e-05  max resid 0.0002659567 
    ## ... Similar to previous best
    ## Run 326 stress 0.07927535 
    ## ... Procrustes: rmse 2.85925e-05  max resid 7.517633e-05 
    ## ... Similar to previous best
    ## Run 327 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734797  max resid 0.056388 
    ## Run 328 stress 0.0910603 
    ## Run 329 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733544  max resid 0.05630398 
    ## Run 330 stress 0.07927534 
    ## ... Procrustes: rmse 5.630181e-06  max resid 1.375724e-05 
    ## ... Similar to previous best
    ## Run 331 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237814  max resid 0.04972297 
    ## Run 332 stress 0.0898561 
    ## Run 333 stress 0.07951443 
    ## ... Procrustes: rmse 0.01737794  max resid 0.05656931 
    ## Run 334 stress 0.07951443 
    ## ... Procrustes: rmse 0.01734106  max resid 0.0563316 
    ## Run 335 stress 0.07927535 
    ## ... Procrustes: rmse 3.242409e-05  max resid 8.121836e-05 
    ## ... Similar to previous best
    ## Run 336 stress 0.09382074 
    ## Run 337 stress 0.07969258 
    ## ... Procrustes: rmse 0.0123755  max resid 0.04969822 
    ## Run 338 stress 0.0792754 
    ## ... Procrustes: rmse 9.040469e-05  max resid 0.0002542763 
    ## ... Similar to previous best
    ## Run 339 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039761  max resid 0.03242795 
    ## Run 340 stress 0.1108444 
    ## Run 341 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734942  max resid 0.05639975 
    ## Run 342 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173441  max resid 0.05634802 
    ## Run 343 stress 0.07936949 
    ## ... Procrustes: rmse 0.009043689  max resid 0.03246392 
    ## Run 344 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734618  max resid 0.05637865 
    ## Run 345 stress 0.09290718 
    ## Run 346 stress 0.09356315 
    ## Run 347 stress 0.1058826 
    ## Run 348 stress 0.07951446 
    ## ... Procrustes: rmse 0.01739268  max resid 0.05665706 
    ## Run 349 stress 0.0795144 
    ## ... Procrustes: rmse 0.01733594  max resid 0.05631511 
    ## Run 350 stress 0.09072903 
    ## Run 351 stress 0.09388327 
    ## Run 352 stress 0.07927534 
    ## ... Procrustes: rmse 1.456676e-05  max resid 3.065484e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.07951442 
    ## ... Procrustes: rmse 0.017356  max resid 0.05650134 
    ## Run 354 stress 0.1108127 
    ## Run 355 stress 0.07927534 
    ## ... Procrustes: rmse 5.243543e-06  max resid 1.214465e-05 
    ## ... Similar to previous best
    ## Run 356 stress 0.07936951 
    ## ... Procrustes: rmse 0.009033483  max resid 0.03239244 
    ## Run 357 stress 0.09290691 
    ## Run 358 stress 0.08985605 
    ## Run 359 stress 0.09474081 
    ## Run 360 stress 0.07927534 
    ## ... Procrustes: rmse 1.235234e-05  max resid 3.831185e-05 
    ## ... Similar to previous best
    ## Run 361 stress 0.1108468 
    ## Run 362 stress 0.09043738 
    ## Run 363 stress 0.1058826 
    ## Run 364 stress 0.104507 
    ## Run 365 stress 0.09072909 
    ## Run 366 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904486  max resid 0.03246596 
    ## Run 367 stress 0.1089998 
    ## Run 368 stress 0.07936954 
    ## ... Procrustes: rmse 0.009025992  max resid 0.03234588 
    ## Run 369 stress 0.09043765 
    ## Run 370 stress 0.09022782 
    ## Run 371 stress 0.08985606 
    ## Run 372 stress 0.1108463 
    ## Run 373 stress 0.1064572 
    ## Run 374 stress 0.07927534 
    ## ... Procrustes: rmse 1.034262e-05  max resid 2.204064e-05 
    ## ... Similar to previous best
    ## Run 375 stress 0.07927535 
    ## ... Procrustes: rmse 2.921397e-05  max resid 7.74725e-05 
    ## ... Similar to previous best
    ## Run 376 stress 0.07969261 
    ## ... Procrustes: rmse 0.01238926  max resid 0.04971588 
    ## Run 377 stress 0.07927535 
    ## ... Procrustes: rmse 3.74946e-05  max resid 9.731698e-05 
    ## ... Similar to previous best
    ## Run 378 stress 0.09485696 
    ## Run 379 stress 0.07951442 
    ## ... Procrustes: rmse 0.01732916  max resid 0.05626042 
    ## Run 380 stress 0.09106041 
    ## Run 381 stress 0.08985606 
    ## Run 382 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734058  max resid 0.05633576 
    ## Run 383 stress 0.09344092 
    ## Run 384 stress 0.09072892 
    ## Run 385 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734493  max resid 0.05636782 
    ## Run 386 stress 0.07927534 
    ## ... Procrustes: rmse 1.596553e-05  max resid 4.073998e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.09388354 
    ## Run 388 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237859  max resid 0.04972365 
    ## Run 389 stress 0.09106047 
    ## Run 390 stress 0.08985605 
    ## Run 391 stress 0.0917279 
    ## Run 392 stress 0.09544383 
    ## Run 393 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733622  max resid 0.05631033 
    ## Run 394 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733675  max resid 0.056318 
    ## Run 395 stress 0.09172805 
    ## Run 396 stress 0.09388356 
    ## Run 397 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735221  max resid 0.05641701 
    ## Run 398 stress 0.07927534 
    ## ... Procrustes: rmse 1.347557e-05  max resid 4.762202e-05 
    ## ... Similar to previous best
    ## Run 399 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173483  max resid 0.05641575 
    ## Run 400 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173403  max resid 0.05633883 
    ## Run 401 stress 0.07969258 
    ## ... Procrustes: rmse 0.0123749  max resid 0.04970104 
    ## Run 402 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735277  max resid 0.0564086 
    ## Run 403 stress 0.07927534 
    ## ... Procrustes: rmse 8.053676e-06  max resid 2.427859e-05 
    ## ... Similar to previous best
    ## Run 404 stress 0.09106027 
    ## Run 405 stress 0.09290698 
    ## Run 406 stress 0.07936951 
    ## ... Procrustes: rmse 0.00904294  max resid 0.0324515 
    ## Run 407 stress 0.08985603 
    ## Run 408 stress 0.09485703 
    ## Run 409 stress 0.09022785 
    ## Run 410 stress 0.09172797 
    ## Run 411 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734557  max resid 0.05637487 
    ## Run 412 stress 0.08985612 
    ## Run 413 stress 0.07927534 
    ## ... Procrustes: rmse 1.651409e-05  max resid 4.459801e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.08985607 
    ## Run 415 stress 0.07969262 
    ## ... Procrustes: rmse 0.01240682  max resid 0.04989789 
    ## Run 416 stress 0.09290717 
    ## Run 417 stress 0.09043746 
    ## Run 418 stress 0.07927534 
    ## ... Procrustes: rmse 1.371037e-05  max resid 3.281989e-05 
    ## ... Similar to previous best
    ## Run 419 stress 0.09296252 
    ## Run 420 stress 0.106491 
    ## Run 421 stress 0.09072891 
    ## Run 422 stress 0.09172789 
    ## Run 423 stress 0.09022777 
    ## Run 424 stress 0.07927534 
    ## ... Procrustes: rmse 1.015972e-05  max resid 2.938668e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.08985607 
    ## Run 426 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736175  max resid 0.05647638 
    ## Run 427 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041348  max resid 0.03244308 
    ## Run 428 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735034  max resid 0.05640475 
    ## Run 429 stress 0.07927534 
    ## ... Procrustes: rmse 1.088498e-05  max resid 3.081835e-05 
    ## ... Similar to previous best
    ## Run 430 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049734  max resid 0.03248921 
    ## Run 431 stress 0.07927538 
    ## ... Procrustes: rmse 2.082081e-05  max resid 6.808381e-05 
    ## ... Similar to previous best
    ## Run 432 stress 0.09072915 
    ## Run 433 stress 0.07927534 
    ## ... Procrustes: rmse 1.446351e-05  max resid 3.791441e-05 
    ## ... Similar to previous best
    ## Run 434 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736712  max resid 0.05650463 
    ## Run 435 stress 0.1067418 
    ## Run 436 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735284  max resid 0.05643158 
    ## Run 437 stress 0.1044027 
    ## Run 438 stress 0.0904377 
    ## Run 439 stress 0.09247105 
    ## Run 440 stress 0.09296237 
    ## Run 441 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733105  max resid 0.05627998 
    ## Run 442 stress 0.09247103 
    ## Run 443 stress 0.07927534 
    ## ... Procrustes: rmse 5.380757e-06  max resid 1.567569e-05 
    ## ... Similar to previous best
    ## Run 444 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045234  max resid 0.03246829 
    ## Run 445 stress 0.07927534 
    ## ... Procrustes: rmse 1.814248e-05  max resid 4.82553e-05 
    ## ... Similar to previous best
    ## Run 446 stress 0.09344096 
    ## Run 447 stress 0.09344099 
    ## Run 448 stress 0.07927534 
    ## ... Procrustes: rmse 1.087635e-05  max resid 3.19102e-05 
    ## ... Similar to previous best
    ## Run 449 stress 0.1033494 
    ## Run 450 stress 0.09043741 
    ## Run 451 stress 0.07927535 
    ## ... Procrustes: rmse 2.105475e-05  max resid 6.020717e-05 
    ## ... Similar to previous best
    ## Run 452 stress 0.07927534 
    ## ... Procrustes: rmse 1.31457e-05  max resid 3.783705e-05 
    ## ... Similar to previous best
    ## Run 453 stress 0.09022778 
    ## Run 454 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047406  max resid 0.03248231 
    ## Run 455 stress 0.07936955 
    ## ... Procrustes: rmse 0.009024576  max resid 0.03233743 
    ## Run 456 stress 0.07936949 
    ## ... Procrustes: rmse 0.009043721  max resid 0.03245788 
    ## Run 457 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237225  max resid 0.0496718 
    ## Run 458 stress 0.07951442 
    ## ... Procrustes: rmse 0.0173583  max resid 0.05645226 
    ## Run 459 stress 0.09072896 
    ## Run 460 stress 0.106141 
    ## Run 461 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735041  max resid 0.05640542 
    ## Run 462 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238506  max resid 0.04973238 
    ## Run 463 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240072  max resid 0.04985699 
    ## Run 464 stress 0.07969263 
    ## ... Procrustes: rmse 0.01239421  max resid 0.04985624 
    ## Run 465 stress 0.07927536 
    ## ... Procrustes: rmse 5.925857e-05  max resid 0.0001412242 
    ## ... Similar to previous best
    ## Run 466 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735191  max resid 0.05641399 
    ## Run 467 stress 0.09485786 
    ## Run 468 stress 0.07927536 
    ## ... Procrustes: rmse 2.17755e-05  max resid 5.057724e-05 
    ## ... Similar to previous best
    ## Run 469 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734905  max resid 0.05639623 
    ## Run 470 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734155  max resid 0.05634747 
    ## Run 471 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237223  max resid 0.0496475 
    ## Run 472 stress 0.09072898 
    ## Run 473 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733662  max resid 0.05631759 
    ## Run 474 stress 0.090729 
    ## Run 475 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041045  max resid 0.03244029 
    ## Run 476 stress 0.07927536 
    ## ... Procrustes: rmse 5.033717e-05  max resid 0.0001210979 
    ## ... Similar to previous best
    ## Run 477 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736253  max resid 0.05648765 
    ## Run 478 stress 0.07951442 
    ## ... Procrustes: rmse 0.01732959  max resid 0.05626643 
    ## Run 479 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173345  max resid 0.0563024 
    ## Run 480 stress 0.09072913 
    ## Run 481 stress 0.09290693 
    ## Run 482 stress 0.07936951 
    ## ... Procrustes: rmse 0.009033194  max resid 0.0323909 
    ## Run 483 stress 0.1033495 
    ## Run 484 stress 0.09043737 
    ## Run 485 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735561  max resid 0.0564164 
    ## Run 486 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238152  max resid 0.04973841 
    ## Run 487 stress 0.09382079 
    ## Run 488 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735425  max resid 0.05643532 
    ## Run 489 stress 0.07936954 
    ## ... Procrustes: rmse 0.009029971  max resid 0.03237071 
    ## Run 490 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735518  max resid 0.05646313 
    ## Run 491 stress 0.09022781 
    ## Run 492 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238082  max resid 0.04974567 
    ## Run 493 stress 0.09290701 
    ## Run 494 stress 0.09290705 
    ## Run 495 stress 0.0956536 
    ## Run 496 stress 0.1117676 
    ## Run 497 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733679  max resid 0.05634394 
    ## Run 498 stress 0.07927534 
    ## ... Procrustes: rmse 5.066645e-06  max resid 1.209182e-05 
    ## ... Similar to previous best
    ## Run 499 stress 0.0904376 
    ## Run 500 stress 0.07927534 
    ## ... Procrustes: rmse 1.567446e-05  max resid 3.744727e-05 
    ## ... Similar to previous best
    ## *** Best solution repeated 62 times

``` r
round(PD_beta_NMDS$stress, digits = 2)
```

    ## [1] 0.08

``` r
PD_beta_rep_NMDS <- metaMDS(PD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.3173923 
    ## Run 1 stress 0.3322988 
    ## Run 2 stress 0.3427269 
    ## Run 3 stress 0.3306792 
    ## Run 4 stress 0.3180892 
    ## Run 5 stress 0.319705 
    ## Run 6 stress 0.3239742 
    ## Run 7 stress 0.3155733 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1620467  max resid 0.2901415 
    ## Run 8 stress 0.3178759 
    ## Run 9 stress 0.3353028 
    ## Run 10 stress 0.3258488 
    ## Run 11 stress 0.3403249 
    ## Run 12 stress 0.3145029 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1404084  max resid 0.3668316 
    ## Run 13 stress 0.337815 
    ## Run 14 stress 0.337242 
    ## Run 15 stress 0.3339303 
    ## Run 16 stress 0.3412496 
    ## Run 17 stress 0.3279736 
    ## Run 18 stress 0.3287229 
    ## Run 19 stress 0.3391044 
    ## Run 20 stress 0.3336656 
    ## Run 21 stress 0.318106 
    ## Run 22 stress 0.3326232 
    ## Run 23 stress 0.328546 
    ## Run 24 stress 0.3249052 
    ## Run 25 stress 0.322298 
    ## Run 26 stress 0.3198603 
    ## Run 27 stress 0.31708 
    ## Run 28 stress 0.320144 
    ## Run 29 stress 0.3375357 
    ## Run 30 stress 0.3132855 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1901384  max resid 0.304453 
    ## Run 31 stress 0.3276564 
    ## Run 32 stress 0.3396109 
    ## Run 33 stress 0.3221426 
    ## Run 34 stress 0.3225219 
    ## Run 35 stress 0.3310708 
    ## Run 36 stress 0.3306469 
    ## Run 37 stress 0.3183461 
    ## Run 38 stress 0.3218412 
    ## Run 39 stress 0.3380031 
    ## Run 40 stress 0.330313 
    ## Run 41 stress 0.3396594 
    ## Run 42 stress 0.3315276 
    ## Run 43 stress 0.3154265 
    ## Run 44 stress 0.3193572 
    ## Run 45 stress 0.3348483 
    ## Run 46 stress 0.3176677 
    ## Run 47 stress 0.3170649 
    ## Run 48 stress 0.3180889 
    ## Run 49 stress 0.3189063 
    ## Run 50 stress 0.331389 
    ## Run 51 stress 0.3142239 
    ## Run 52 stress 0.3203537 
    ## Run 53 stress 0.3321956 
    ## Run 54 stress 0.3204027 
    ## Run 55 stress 0.3219694 
    ## Run 56 stress 0.3203325 
    ## Run 57 stress 0.329386 
    ## Run 58 stress 0.3381236 
    ## Run 59 stress 0.3324636 
    ## Run 60 stress 0.3173344 
    ## Run 61 stress 0.3173063 
    ## Run 62 stress 0.3322698 
    ## Run 63 stress 0.3258988 
    ## Run 64 stress 0.3228891 
    ## Run 65 stress 0.3286575 
    ## Run 66 stress 0.3161783 
    ## Run 67 stress 0.3366513 
    ## Run 68 stress 0.316077 
    ## Run 69 stress 0.3317899 
    ## Run 70 stress 0.3236698 
    ## Run 71 stress 0.3334303 
    ## Run 72 stress 0.3391191 
    ## Run 73 stress 0.3346212 
    ## Run 74 stress 0.3352311 
    ## Run 75 stress 0.3251544 
    ## Run 76 stress 0.3218803 
    ## Run 77 stress 0.3209812 
    ## Run 78 stress 0.3174694 
    ## Run 79 stress 0.3262357 
    ## Run 80 stress 0.350655 
    ## Run 81 stress 0.324755 
    ## Run 82 stress 0.3226601 
    ## Run 83 stress 0.3287336 
    ## Run 84 stress 0.3149565 
    ## Run 85 stress 0.3168405 
    ## Run 86 stress 0.3407252 
    ## Run 87 stress 0.3318456 
    ## Run 88 stress 0.3191628 
    ## Run 89 stress 0.3340988 
    ## Run 90 stress 0.3227052 
    ## Run 91 stress 0.3207856 
    ## Run 92 stress 0.3260505 
    ## Run 93 stress 0.3212453 
    ## Run 94 stress 0.3328038 
    ## Run 95 stress 0.3285335 
    ## Run 96 stress 0.330948 
    ## Run 97 stress 0.3317852 
    ## Run 98 stress 0.3287225 
    ## Run 99 stress 0.3299457 
    ## Run 100 stress 0.3313981 
    ## Run 101 stress 0.327495 
    ## Run 102 stress 0.3201849 
    ## Run 103 stress 0.3205968 
    ## Run 104 stress 0.3265402 
    ## Run 105 stress 0.3219321 
    ## Run 106 stress 0.3210406 
    ## Run 107 stress 0.3162139 
    ## Run 108 stress 0.3320626 
    ## Run 109 stress 0.3165161 
    ## Run 110 stress 0.3272218 
    ## Run 111 stress 0.3264906 
    ## Run 112 stress 0.3255794 
    ## Run 113 stress 0.3278001 
    ## Run 114 stress 0.3162543 
    ## Run 115 stress 0.3327844 
    ## Run 116 stress 0.3202837 
    ## Run 117 stress 0.3143221 
    ## Run 118 stress 0.3213747 
    ## Run 119 stress 0.3164712 
    ## Run 120 stress 0.3304627 
    ## Run 121 stress 0.3235437 
    ## Run 122 stress 0.33216 
    ## Run 123 stress 0.3288352 
    ## Run 124 stress 0.3378333 
    ## Run 125 stress 0.3253854 
    ## Run 126 stress 0.3236073 
    ## Run 127 stress 0.319947 
    ## Run 128 stress 0.3371901 
    ## Run 129 stress 0.3326663 
    ## Run 130 stress 0.3219044 
    ## Run 131 stress 0.3316697 
    ## Run 132 stress 0.3179464 
    ## Run 133 stress 0.3250941 
    ## Run 134 stress 0.337319 
    ## Run 135 stress 0.3123683 
    ## ... New best solution
    ## ... Procrustes: rmse 0.145526  max resid 0.2822 
    ## Run 136 stress 0.331995 
    ## Run 137 stress 0.3341931 
    ## Run 138 stress 0.3258278 
    ## Run 139 stress 0.326534 
    ## Run 140 stress 0.3336229 
    ## Run 141 stress 0.3248705 
    ## Run 142 stress 0.3177189 
    ## Run 143 stress 0.3156598 
    ## Run 144 stress 0.3367971 
    ## Run 145 stress 0.3308131 
    ## Run 146 stress 0.3223904 
    ## Run 147 stress 0.3202451 
    ## Run 148 stress 0.3196444 
    ## Run 149 stress 0.3338749 
    ## Run 150 stress 0.3273081 
    ## Run 151 stress 0.3266167 
    ## Run 152 stress 0.3206419 
    ## Run 153 stress 0.3296067 
    ## Run 154 stress 0.3229911 
    ## Run 155 stress 0.3359087 
    ## Run 156 stress 0.3211936 
    ## Run 157 stress 0.3169759 
    ## Run 158 stress 0.3333803 
    ## Run 159 stress 0.3183342 
    ## Run 160 stress 0.3370704 
    ## Run 161 stress 0.3113196 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1404822  max resid 0.3391763 
    ## Run 162 stress 0.3251533 
    ## Run 163 stress 0.3310952 
    ## Run 164 stress 0.3245905 
    ## Run 165 stress 0.3367571 
    ## Run 166 stress 0.3243648 
    ## Run 167 stress 0.3306907 
    ## Run 168 stress 0.3308694 
    ## Run 169 stress 0.3253829 
    ## Run 170 stress 0.3194207 
    ## Run 171 stress 0.3210934 
    ## Run 172 stress 0.3253607 
    ## Run 173 stress 0.3298636 
    ## Run 174 stress 0.3220035 
    ## Run 175 stress 0.3183424 
    ## Run 176 stress 0.3245744 
    ## Run 177 stress 0.3188635 
    ## Run 178 stress 0.3358183 
    ## Run 179 stress 0.315741 
    ## Run 180 stress 0.3335869 
    ## Run 181 stress 0.3370335 
    ## Run 182 stress 0.3216145 
    ## Run 183 stress 0.3218857 
    ## Run 184 stress 0.3132676 
    ## Run 185 stress 0.3235683 
    ## Run 186 stress 0.3320918 
    ## Run 187 stress 0.3286269 
    ## Run 188 stress 0.3255551 
    ## Run 189 stress 0.3422136 
    ## Run 190 stress 0.3310946 
    ## Run 191 stress 0.322497 
    ## Run 192 stress 0.3149857 
    ## Run 193 stress 0.3167892 
    ## Run 194 stress 0.3291789 
    ## Run 195 stress 0.324835 
    ## Run 196 stress 0.3170441 
    ## Run 197 stress 0.3165024 
    ## Run 198 stress 0.3335059 
    ## Run 199 stress 0.3328261 
    ## Run 200 stress 0.3232435 
    ## Run 201 stress 0.3300464 
    ## Run 202 stress 0.3306342 
    ## Run 203 stress 0.3202488 
    ## Run 204 stress 0.3289099 
    ## Run 205 stress 0.3095954 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0910219  max resid 0.239616 
    ## Run 206 stress 0.315386 
    ## Run 207 stress 0.332076 
    ## Run 208 stress 0.3246737 
    ## Run 209 stress 0.3329134 
    ## Run 210 stress 0.323129 
    ## Run 211 stress 0.33727 
    ## Run 212 stress 0.3327072 
    ## Run 213 stress 0.3182866 
    ## Run 214 stress 0.3300169 
    ## Run 215 stress 0.3347279 
    ## Run 216 stress 0.3262241 
    ## Run 217 stress 0.3198784 
    ## Run 218 stress 0.3265885 
    ## Run 219 stress 0.3214354 
    ## Run 220 stress 0.324018 
    ## Run 221 stress 0.3154769 
    ## Run 222 stress 0.3352714 
    ## Run 223 stress 0.332781 
    ## Run 224 stress 0.3252269 
    ## Run 225 stress 0.3268245 
    ## Run 226 stress 0.3248908 
    ## Run 227 stress 0.3194973 
    ## Run 228 stress 0.3199729 
    ## Run 229 stress 0.3227505 
    ## Run 230 stress 0.3188086 
    ## Run 231 stress 0.3179882 
    ## Run 232 stress 0.3374711 
    ## Run 233 stress 0.3235127 
    ## Run 234 stress 0.3357145 
    ## Run 235 stress 0.3279562 
    ## Run 236 stress 0.3276201 
    ## Run 237 stress 0.334536 
    ## Run 238 stress 0.3246522 
    ## Run 239 stress 0.3313416 
    ## Run 240 stress 0.3146755 
    ## Run 241 stress 0.3808263 
    ## Run 242 stress 0.321902 
    ## Run 243 stress 0.3332715 
    ## Run 244 stress 0.327872 
    ## Run 245 stress 0.3208861 
    ## Run 246 stress 0.3386641 
    ## Run 247 stress 0.3144303 
    ## Run 248 stress 0.3202046 
    ## Run 249 stress 0.3193472 
    ## Run 250 stress 0.3323246 
    ## Run 251 stress 0.3230002 
    ## Run 252 stress 0.3303845 
    ## Run 253 stress 0.3211793 
    ## Run 254 stress 0.3188125 
    ## Run 255 stress 0.3237842 
    ## Run 256 stress 0.3170776 
    ## Run 257 stress 0.3341022 
    ## Run 258 stress 0.3349642 
    ## Run 259 stress 0.3200974 
    ## Run 260 stress 0.3224823 
    ## Run 261 stress 0.3208487 
    ## Run 262 stress 0.3362889 
    ## Run 263 stress 0.3241447 
    ## Run 264 stress 0.3232978 
    ## Run 265 stress 0.3175027 
    ## Run 266 stress 0.3222255 
    ## Run 267 stress 0.3301758 
    ## Run 268 stress 0.3172747 
    ## Run 269 stress 0.3281922 
    ## Run 270 stress 0.3401919 
    ## Run 271 stress 0.3291765 
    ## Run 272 stress 0.3255182 
    ## Run 273 stress 0.318717 
    ## Run 274 stress 0.3327463 
    ## Run 275 stress 0.3277947 
    ## Run 276 stress 0.3337313 
    ## Run 277 stress 0.3308405 
    ## Run 278 stress 0.3272854 
    ## Run 279 stress 0.3203856 
    ## Run 280 stress 0.3318865 
    ## Run 281 stress 0.3219576 
    ## Run 282 stress 0.3368783 
    ## Run 283 stress 0.321402 
    ## Run 284 stress 0.3269698 
    ## Run 285 stress 0.3350553 
    ## Run 286 stress 0.3330057 
    ## Run 287 stress 0.3283205 
    ## Run 288 stress 0.3169238 
    ## Run 289 stress 0.3356722 
    ## Run 290 stress 0.336331 
    ## Run 291 stress 0.3290185 
    ## Run 292 stress 0.3261331 
    ## Run 293 stress 0.3307822 
    ## Run 294 stress 0.3228925 
    ## Run 295 stress 0.3310021 
    ## Run 296 stress 0.3342068 
    ## Run 297 stress 0.3264534 
    ## Run 298 stress 0.3175519 
    ## Run 299 stress 0.3278358 
    ## Run 300 stress 0.3266885 
    ## Run 301 stress 0.3356301 
    ## Run 302 stress 0.3317406 
    ## Run 303 stress 0.3357753 
    ## Run 304 stress 0.3395668 
    ## Run 305 stress 0.31385 
    ## Run 306 stress 0.3330738 
    ## Run 307 stress 0.3182474 
    ## Run 308 stress 0.3296252 
    ## Run 309 stress 0.3208508 
    ## Run 310 stress 0.3366562 
    ## Run 311 stress 0.3194903 
    ## Run 312 stress 0.3278788 
    ## Run 313 stress 0.3258078 
    ## Run 314 stress 0.3139029 
    ## Run 315 stress 0.3206895 
    ## Run 316 stress 0.3154882 
    ## Run 317 stress 0.3208358 
    ## Run 318 stress 0.3231396 
    ## Run 319 stress 0.3347094 
    ## Run 320 stress 0.3286226 
    ## Run 321 stress 0.3253077 
    ## Run 322 stress 0.314503 
    ## Run 323 stress 0.3219924 
    ## Run 324 stress 0.3175029 
    ## Run 325 stress 0.3227793 
    ## Run 326 stress 0.329326 
    ## Run 327 stress 0.3246513 
    ## Run 328 stress 0.3343032 
    ## Run 329 stress 0.3207798 
    ## Run 330 stress 0.3217676 
    ## Run 331 stress 0.3261122 
    ## Run 332 stress 0.319325 
    ## Run 333 stress 0.3304712 
    ## Run 334 stress 0.322432 
    ## Run 335 stress 0.3386216 
    ## Run 336 stress 0.3254053 
    ## Run 337 stress 0.3217976 
    ## Run 338 stress 0.3289614 
    ## Run 339 stress 0.3366606 
    ## Run 340 stress 0.3292818 
    ## Run 341 stress 0.3134163 
    ## Run 342 stress 0.3303645 
    ## Run 343 stress 0.3249074 
    ## Run 344 stress 0.3223335 
    ## Run 345 stress 0.3273307 
    ## Run 346 stress 0.3214538 
    ## Run 347 stress 0.3573876 
    ## Run 348 stress 0.3378695 
    ## Run 349 stress 0.3221583 
    ## Run 350 stress 0.3134762 
    ## Run 351 stress 0.3399396 
    ## Run 352 stress 0.3233118 
    ## Run 353 stress 0.3167788 
    ## Run 354 stress 0.3175295 
    ## Run 355 stress 0.3162794 
    ## Run 356 stress 0.3345682 
    ## Run 357 stress 0.3329109 
    ## Run 358 stress 0.3186089 
    ## Run 359 stress 0.3293305 
    ## Run 360 stress 0.3250883 
    ## Run 361 stress 0.3258126 
    ## Run 362 stress 0.3297467 
    ## Run 363 stress 0.3359582 
    ## Run 364 stress 0.3203515 
    ## Run 365 stress 0.3372872 
    ## Run 366 stress 0.3288611 
    ## Run 367 stress 0.3350303 
    ## Run 368 stress 0.3198222 
    ## Run 369 stress 0.3166933 
    ## Run 370 stress 0.3225584 
    ## Run 371 stress 0.330472 
    ## Run 372 stress 0.332479 
    ## Run 373 stress 0.314256 
    ## Run 374 stress 0.3219792 
    ## Run 375 stress 0.3257399 
    ## Run 376 stress 0.3226422 
    ## Run 377 stress 0.3266382 
    ## Run 378 stress 0.3252573 
    ## Run 379 stress 0.3168409 
    ## Run 380 stress 0.3288346 
    ## Run 381 stress 0.3346904 
    ## Run 382 stress 0.324026 
    ## Run 383 stress 0.3186148 
    ## Run 384 stress 0.3773351 
    ## Run 385 stress 0.3179977 
    ## Run 386 stress 0.3121692 
    ## Run 387 stress 0.3316393 
    ## Run 388 stress 0.3215163 
    ## Run 389 stress 0.3266617 
    ## Run 390 stress 0.3231316 
    ## Run 391 stress 0.3288164 
    ## Run 392 stress 0.3193297 
    ## Run 393 stress 0.3386561 
    ## Run 394 stress 0.33396 
    ## Run 395 stress 0.3179045 
    ## Run 396 stress 0.3252367 
    ## Run 397 stress 0.3259782 
    ## Run 398 stress 0.3140382 
    ## Run 399 stress 0.3296453 
    ## Run 400 stress 0.3262998 
    ## Run 401 stress 0.3245308 
    ## Run 402 stress 0.3170205 
    ## Run 403 stress 0.3209848 
    ## Run 404 stress 0.3158899 
    ## Run 405 stress 0.3348369 
    ## Run 406 stress 0.3229763 
    ## Run 407 stress 0.3298628 
    ## Run 408 stress 0.3242306 
    ## Run 409 stress 0.3269688 
    ## Run 410 stress 0.3305753 
    ## Run 411 stress 0.3299727 
    ## Run 412 stress 0.3317328 
    ## Run 413 stress 0.3226989 
    ## Run 414 stress 0.3315044 
    ## Run 415 stress 0.3348065 
    ## Run 416 stress 0.3233819 
    ## Run 417 stress 0.3299358 
    ## Run 418 stress 0.3367748 
    ## Run 419 stress 0.3246876 
    ## Run 420 stress 0.3262325 
    ## Run 421 stress 0.3238811 
    ## Run 422 stress 0.3261184 
    ## Run 423 stress 0.3355517 
    ## Run 424 stress 0.3304241 
    ## Run 425 stress 0.3124763 
    ## Run 426 stress 0.326165 
    ## Run 427 stress 0.3294632 
    ## Run 428 stress 0.3144183 
    ## Run 429 stress 0.3293881 
    ## Run 430 stress 0.3311556 
    ## Run 431 stress 0.3253636 
    ## Run 432 stress 0.3329236 
    ## Run 433 stress 0.3397078 
    ## Run 434 stress 0.3246789 
    ## Run 435 stress 0.3145709 
    ## Run 436 stress 0.3252987 
    ## Run 437 stress 0.3159956 
    ## Run 438 stress 0.3210813 
    ## Run 439 stress 0.3361765 
    ## Run 440 stress 0.3300249 
    ## Run 441 stress 0.318326 
    ## Run 442 stress 0.3218791 
    ## Run 443 stress 0.3263737 
    ## Run 444 stress 0.3509927 
    ## Run 445 stress 0.3151603 
    ## Run 446 stress 0.3192619 
    ## Run 447 stress 0.3393202 
    ## Run 448 stress 0.3355247 
    ## Run 449 stress 0.3192395 
    ## Run 450 stress 0.3187857 
    ## Run 451 stress 0.3260075 
    ## Run 452 stress 0.3197286 
    ## Run 453 stress 0.3155374 
    ## Run 454 stress 0.3122497 
    ## Run 455 stress 0.3206702 
    ## Run 456 stress 0.3234951 
    ## Run 457 stress 0.3226032 
    ## Run 458 stress 0.3197161 
    ## Run 459 stress 0.3221439 
    ## Run 460 stress 0.3373759 
    ## Run 461 stress 0.3175135 
    ## Run 462 stress 0.3263423 
    ## Run 463 stress 0.335817 
    ## Run 464 stress 0.3329205 
    ## Run 465 stress 0.3269546 
    ## Run 466 stress 0.3337935 
    ## Run 467 stress 0.326733 
    ## Run 468 stress 0.3247692 
    ## Run 469 stress 0.3241143 
    ## Run 470 stress 0.3284956 
    ## Run 471 stress 0.3312619 
    ## Run 472 stress 0.3249441 
    ## Run 473 stress 0.3220263 
    ## Run 474 stress 0.3328042 
    ## Run 475 stress 0.321236 
    ## Run 476 stress 0.3250967 
    ## Run 477 stress 0.3177798 
    ## Run 478 stress 0.3284424 
    ## Run 479 stress 0.336047 
    ## Run 480 stress 0.3190025 
    ## Run 481 stress 0.3301114 
    ## Run 482 stress 0.3253761 
    ## Run 483 stress 0.3320912 
    ## Run 484 stress 0.3164018 
    ## Run 485 stress 0.3360118 
    ## Run 486 stress 0.3330361 
    ## Run 487 stress 0.3216888 
    ## Run 488 stress 0.3255734 
    ## Run 489 stress 0.3389702 
    ## Run 490 stress 0.3243904 
    ## Run 491 stress 0.3217745 
    ## Run 492 stress 0.3263399 
    ## Run 493 stress 0.3277397 
    ## Run 494 stress 0.3340399 
    ## Run 495 stress 0.3290904 
    ## Run 496 stress 0.3352693 
    ## Run 497 stress 0.3215114 
    ## Run 498 stress 0.3140585 
    ## Run 499 stress 0.3153243 
    ## Run 500 stress 0.3177117 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##    500: stress ratio > sratmax

``` r
PD_beta_ric_NMDS <- metaMDS(PD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.01115177 
    ## Run 1 stress 0.01130817 
    ## ... Procrustes: rmse 0.01032205  max resid 0.02320562 
    ## Run 2 stress 0.01517348 
    ## Run 3 stress 0.01536148 
    ## Run 4 stress 0.01528189 
    ## Run 5 stress 0.01542801 
    ## Run 6 stress 0.01115177 
    ## ... Procrustes: rmse 3.22524e-06  max resid 7.044025e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.01659949 
    ## Run 8 stress 0.01558333 
    ## Run 9 stress 0.01515719 
    ## Run 10 stress 0.01648718 
    ## Run 11 stress 0.01115239 
    ## ... Procrustes: rmse 0.001611141  max resid 0.003632785 
    ## ... Similar to previous best
    ## Run 12 stress 0.01115188 
    ## ... Procrustes: rmse 3.810472e-05  max resid 8.599449e-05 
    ## ... Similar to previous best
    ## Run 13 stress 0.01515723 
    ## Run 14 stress 0.01115183 
    ## ... Procrustes: rmse 2.528346e-05  max resid 5.723187e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.01499443 
    ## Run 16 stress 0.01557929 
    ## Run 17 stress 0.01528211 
    ## Run 18 stress 0.01547425 
    ## Run 19 stress 0.01116612 
    ## ... Procrustes: rmse 0.003687103  max resid 0.008310261 
    ## ... Similar to previous best
    ## Run 20 stress 0.01547434 
    ## Run 21 stress 0.01533788 
    ## Run 22 stress 0.01493159 
    ## Run 23 stress 0.01641872 
    ## Run 24 stress 0.01557665 
    ## Run 25 stress 0.01147956 
    ## ... Procrustes: rmse 0.0131094  max resid 0.02863804 
    ## Run 26 stress 0.01554831 
    ## Run 27 stress 0.01561832 
    ## Run 28 stress 0.01581909 
    ## Run 29 stress 0.01572108 
    ## Run 30 stress 0.01499095 
    ## Run 31 stress 0.01529089 
    ## Run 32 stress 0.01115174 
    ## ... New best solution
    ## ... Procrustes: rmse 1.183984e-05  max resid 2.678197e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.01571541 
    ## Run 34 stress 0.01498599 
    ## Run 35 stress 0.01115261 
    ## ... Procrustes: rmse 0.001663788  max resid 0.003751525 
    ## ... Similar to previous best
    ## Run 36 stress 0.01558377 
    ## Run 37 stress 0.01648729 
    ## Run 38 stress 0.01515712 
    ## Run 39 stress 0.01556111 
    ## Run 40 stress 0.01556015 
    ## Run 41 stress 0.01515674 
    ## Run 42 stress 0.01641897 
    ## Run 43 stress 0.01562233 
    ## Run 44 stress 0.01525787 
    ## Run 45 stress 0.01493179 
    ## Run 46 stress 0.01533763 
    ## Run 47 stress 0.01115227 
    ## ... Procrustes: rmse 0.001559852  max resid 0.003516957 
    ## ... Similar to previous best
    ## Run 48 stress 0.01115164 
    ## ... New best solution
    ## ... Procrustes: rmse 4.555051e-05  max resid 0.0001026186 
    ## ... Similar to previous best
    ## Run 49 stress 0.01115163 
    ## ... New best solution
    ## ... Procrustes: rmse 3.272592e-06  max resid 7.377869e-06 
    ## ... Similar to previous best
    ## Run 50 stress 0.01556729 
    ## Run 51 stress 0.01498589 
    ## Run 52 stress 0.0154994 
    ## Run 53 stress 0.01562086 
    ## Run 54 stress 0.01499272 
    ## Run 55 stress 0.01528214 
    ## Run 56 stress 0.01527924 
    ## Run 57 stress 0.01549927 
    ## Run 58 stress 0.01556982 
    ## Run 59 stress 0.01658168 
    ## Run 60 stress 0.01576556 
    ## Run 61 stress 0.01527903 
    ## Run 62 stress 0.01382445 
    ## Run 63 stress 0.01498592 
    ## Run 64 stress 0.01557881 
    ## Run 65 stress 0.01528202 
    ## Run 66 stress 0.01115175 
    ## ... Procrustes: rmse 5.236919e-05  max resid 0.0001184473 
    ## ... Similar to previous best
    ## Run 67 stress 0.0157153 
    ## Run 68 stress 0.01528183 
    ## Run 69 stress 0.01498603 
    ## Run 70 stress 0.01513836 
    ## Run 71 stress 0.2963757 
    ## Run 72 stress 0.01498582 
    ## Run 73 stress 0.01498583 
    ## Run 74 stress 0.01576336 
    ## Run 75 stress 0.01498587 
    ## Run 76 stress 0.01528189 
    ## Run 77 stress 0.01499131 
    ## Run 78 stress 0.01648735 
    ## Run 79 stress 0.01527931 
    ## Run 80 stress 0.01515661 
    ## Run 81 stress 0.01529486 
    ## Run 82 stress 0.01718695 
    ## Run 83 stress 0.0111519 
    ## ... Procrustes: rmse 0.001372299  max resid 0.003094482 
    ## ... Similar to previous best
    ## Run 84 stress 0.01115183 
    ## ... Procrustes: rmse 8.574901e-05  max resid 0.0001938273 
    ## ... Similar to previous best
    ## Run 85 stress 0.01549943 
    ## Run 86 stress 0.0155609 
    ## Run 87 stress 0.01557925 
    ## Run 88 stress 0.01687631 
    ## Run 89 stress 0.01567356 
    ## Run 90 stress 0.01498584 
    ## Run 91 stress 0.01549942 
    ## Run 92 stress 0.01152831 
    ## ... Procrustes: rmse 0.01427101  max resid 0.03094784 
    ## Run 93 stress 0.0154518 
    ## Run 94 stress 0.015157 
    ## Run 95 stress 0.01594376 
    ## Run 96 stress 0.01556101 
    ## Run 97 stress 0.01545043 
    ## Run 98 stress 0.01564541 
    ## Run 99 stress 0.01715732 
    ## Run 100 stress 0.3807932 
    ## Run 101 stress 0.01529067 
    ## Run 102 stress 0.0154547 
    ## Run 103 stress 0.01116897 
    ## ... Procrustes: rmse 0.003889223  max resid 0.008763748 
    ## ... Similar to previous best
    ## Run 104 stress 0.0151566 
    ## Run 105 stress 0.01557 
    ## Run 106 stress 0.01552954 
    ## Run 107 stress 0.01577177 
    ## Run 108 stress 0.01117119 
    ## ... Procrustes: rmse 0.004090091  max resid 0.009217888 
    ## ... Similar to previous best
    ## Run 109 stress 0.01514972 
    ## Run 110 stress 0.01551346 
    ## Run 111 stress 0.01498599 
    ## Run 112 stress 0.01572372 
    ## Run 113 stress 0.0156456 
    ## Run 114 stress 0.01115204 
    ## ... Procrustes: rmse 0.0001610557  max resid 0.0003646003 
    ## ... Similar to previous best
    ## Run 115 stress 0.015439 
    ## Run 116 stress 0.0154819 
    ## Run 117 stress 0.01529077 
    ## Run 118 stress 0.01128808 
    ## ... Procrustes: rmse 0.009658659  max resid 0.02172263 
    ## Run 119 stress 0.01515887 
    ## Run 120 stress 0.01115184 
    ## ... Procrustes: rmse 0.001339772  max resid 0.003021314 
    ## ... Similar to previous best
    ## Run 121 stress 0.01545345 
    ## Run 122 stress 0.01181268 
    ## Run 123 stress 0.01515078 
    ## Run 124 stress 0.01686732 
    ## Run 125 stress 0.01648626 
    ## Run 126 stress 0.01562148 
    ## Run 127 stress 0.01528184 
    ## Run 128 stress 0.01543175 
    ## Run 129 stress 0.01514942 
    ## Run 130 stress 0.01493159 
    ## Run 131 stress 0.01556099 
    ## Run 132 stress 0.01528189 
    ## Run 133 stress 0.01574912 
    ## Run 134 stress 0.01637134 
    ## Run 135 stress 0.01545773 
    ## Run 136 stress 0.01498601 
    ## Run 137 stress 0.01115181 
    ## ... Procrustes: rmse 0.001338711  max resid 0.003018317 
    ## ... Similar to previous best
    ## Run 138 stress 0.01547454 
    ## Run 139 stress 0.01582993 
    ## Run 140 stress 0.01493179 
    ## Run 141 stress 0.01571563 
    ## Run 142 stress 0.01694539 
    ## Run 143 stress 0.01127533 
    ## ... Procrustes: rmse 0.00923767  max resid 0.02077225 
    ## Run 144 stress 0.0111518 
    ## ... Procrustes: rmse 7.470584e-05  max resid 0.0001687936 
    ## ... Similar to previous best
    ## Run 145 stress 0.01117871 
    ## ... Procrustes: rmse 0.004684189  max resid 0.01055517 
    ## Run 146 stress 0.01556202 
    ## Run 147 stress 0.01115198 
    ## ... Procrustes: rmse 0.0014073  max resid 0.003173059 
    ## ... Similar to previous best
    ## Run 148 stress 0.01115149 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001183934  max resid 0.002669284 
    ## ... Similar to previous best
    ## Run 149 stress 0.01571573 
    ## Run 150 stress 0.01499279 
    ## Run 151 stress 0.0111518 
    ## ... Procrustes: rmse 0.0001478896  max resid 0.000334162 
    ## ... Similar to previous best
    ## Run 152 stress 0.01675859 
    ## Run 153 stress 0.01648721 
    ## Run 154 stress 0.01550291 
    ## Run 155 stress 0.01663487 
    ## Run 156 stress 0.01528181 
    ## Run 157 stress 0.01115203 
    ## ... Procrustes: rmse 0.001339541  max resid 0.003024874 
    ## ... Similar to previous best
    ## Run 158 stress 0.01533726 
    ## Run 159 stress 0.0155395 
    ## Run 160 stress 0.01493154 
    ## Run 161 stress 0.01560034 
    ## Run 162 stress 0.01525856 
    ## Run 163 stress 0.01115665 
    ## ... Procrustes: rmse 0.001282952  max resid 0.00289516 
    ## ... Similar to previous best
    ## Run 164 stress 0.01493162 
    ## Run 165 stress 0.0151573 
    ## Run 166 stress 0.01558404 
    ## Run 167 stress 0.01547398 
    ## Run 168 stress 0.01120004 
    ## ... Procrustes: rmse 0.004836785  max resid 0.01089805 
    ## Run 169 stress 0.01116353 
    ## ... Procrustes: rmse 0.0021772  max resid 0.00491156 
    ## ... Similar to previous best
    ## Run 170 stress 0.01124223 
    ## ... Procrustes: rmse 0.005619886  max resid 0.01266728 
    ## Run 171 stress 0.01549952 
    ## Run 172 stress 0.01694528 
    ## Run 173 stress 0.01115181 
    ## ... Procrustes: rmse 0.001257373  max resid 0.002837975 
    ## ... Similar to previous best
    ## Run 174 stress 0.01499296 
    ## Run 175 stress 0.01121265 
    ## ... Procrustes: rmse 0.005492372  max resid 0.01236936 
    ## Run 176 stress 0.01142391 
    ## ... Procrustes: rmse 0.01071352  max resid 0.02361565 
    ## Run 177 stress 0.01641869 
    ## Run 178 stress 0.01545148 
    ## Run 179 stress 0.01515676 
    ## Run 180 stress 0.01303576 
    ## Run 181 stress 0.01498585 
    ## Run 182 stress 0.01542977 
    ## Run 183 stress 0.01583101 
    ## Run 184 stress 0.0164862 
    ## Run 185 stress 0.01663516 
    ## Run 186 stress 0.01528164 
    ## Run 187 stress 0.0164872 
    ## Run 188 stress 0.01115154 
    ## ... Procrustes: rmse 0.001137889  max resid 0.002567974 
    ## ... Similar to previous best
    ## Run 189 stress 0.01556105 
    ## Run 190 stress 0.0154913 
    ## Run 191 stress 0.01498596 
    ## Run 192 stress 0.01556956 
    ## Run 193 stress 0.01529117 
    ## Run 194 stress 0.01115165 
    ## ... Procrustes: rmse 0.001182428  max resid 0.002668907 
    ## ... Similar to previous best
    ## Run 195 stress 0.0153613 
    ## Run 196 stress 0.01663479 
    ## Run 197 stress 0.01561048 
    ## Run 198 stress 0.0158195 
    ## Run 199 stress 0.01581983 
    ## Run 200 stress 0.01549931 
    ## Run 201 stress 0.01499269 
    ## Run 202 stress 0.01658353 
    ## Run 203 stress 0.01164339 
    ## ... Procrustes: rmse 0.01661525  max resid 0.0357829 
    ## Run 204 stress 0.01547432 
    ## Run 205 stress 0.01658143 
    ## Run 206 stress 0.01528194 
    ## Run 207 stress 0.01594434 
    ## Run 208 stress 0.0154991 
    ## Run 209 stress 0.01581925 
    ## Run 210 stress 0.01562179 
    ## Run 211 stress 0.01120161 
    ## ... Procrustes: rmse 0.004922425  max resid 0.01109158 
    ## Run 212 stress 0.01115196 
    ## ... Procrustes: rmse 0.0002083188  max resid 0.0004706613 
    ## ... Similar to previous best
    ## Run 213 stress 0.01549945 
    ## Run 214 stress 0.01115248 
    ## ... Procrustes: rmse 0.0003930052  max resid 0.0008870135 
    ## ... Similar to previous best
    ## Run 215 stress 0.01115194 
    ## ... Procrustes: rmse 0.001309986  max resid 0.002957863 
    ## ... Similar to previous best
    ## Run 216 stress 0.0155599 
    ## Run 217 stress 0.01498586 
    ## Run 218 stress 0.01583003 
    ## Run 219 stress 0.01561962 
    ## Run 220 stress 0.01574951 
    ## Run 221 stress 0.01148923 
    ## ... Procrustes: rmse 0.01212585  max resid 0.02631301 
    ## Run 222 stress 0.01527908 
    ## Run 223 stress 0.01557922 
    ## Run 224 stress 0.01528215 
    ## Run 225 stress 0.01525818 
    ## Run 226 stress 0.01564594 
    ## Run 227 stress 0.01528203 
    ## Run 228 stress 0.01169108 
    ## Run 229 stress 0.0125828 
    ## Run 230 stress 0.01567205 
    ## Run 231 stress 0.0170028 
    ## Run 232 stress 0.01545344 
    ## Run 233 stress 0.01115256 
    ## ... Procrustes: rmse 0.0004197408  max resid 0.0009473866 
    ## ... Similar to previous best
    ## Run 234 stress 0.01584389 
    ## Run 235 stress 0.01536139 
    ## Run 236 stress 0.0156225 
    ## Run 237 stress 0.01529934 
    ## Run 238 stress 0.01121011 
    ## ... Procrustes: rmse 0.005367208  max resid 0.01208892 
    ## Run 239 stress 0.01116725 
    ## ... Procrustes: rmse 0.002552176  max resid 0.005757744 
    ## ... Similar to previous best
    ## Run 240 stress 0.01545399 
    ## Run 241 stress 0.0164984 
    ## Run 242 stress 0.0112854 
    ## ... Procrustes: rmse 0.008396098  max resid 0.01889346 
    ## Run 243 stress 0.01122389 
    ## ... Procrustes: rmse 0.005973248  max resid 0.0134543 
    ## Run 244 stress 0.01516063 
    ## Run 245 stress 0.01498595 
    ## Run 246 stress 0.01152877 
    ## ... Procrustes: rmse 0.01314945  max resid 0.02837703 
    ## Run 247 stress 0.01515712 
    ## Run 248 stress 0.01493234 
    ## Run 249 stress 0.01115156 
    ## ... Procrustes: rmse 0.001147563  max resid 0.002589822 
    ## ... Similar to previous best
    ## Run 250 stress 0.01556909 
    ## Run 251 stress 0.01547207 
    ## Run 252 stress 0.01649821 
    ## Run 253 stress 0.01659301 
    ## Run 254 stress 0.01556989 
    ## Run 255 stress 0.01562119 
    ## Run 256 stress 0.01638825 
    ## Run 257 stress 0.01515725 
    ## Run 258 stress 0.01515705 
    ## Run 259 stress 0.01676273 
    ## Run 260 stress 0.01558253 
    ## Run 261 stress 0.01140888 
    ## ... Procrustes: rmse 0.01050952  max resid 0.02327096 
    ## Run 262 stress 0.0152816 
    ## Run 263 stress 0.01637113 
    ## Run 264 stress 0.01178886 
    ## Run 265 stress 0.01698899 
    ## Run 266 stress 0.01115179 
    ## ... Procrustes: rmse 0.0001443611  max resid 0.0003259993 
    ## ... Similar to previous best
    ## Run 267 stress 0.01529479 
    ## Run 268 stress 0.01517258 
    ## Run 269 stress 0.01118347 
    ## ... Procrustes: rmse 0.003833808  max resid 0.008644205 
    ## ... Similar to previous best
    ## Run 270 stress 0.01544627 
    ## Run 271 stress 0.01545402 
    ## Run 272 stress 0.01547265 
    ## Run 273 stress 0.01533771 
    ## Run 274 stress 0.01120098 
    ## ... Procrustes: rmse 0.00489223  max resid 0.01102297 
    ## Run 275 stress 0.01529116 
    ## Run 276 stress 0.01655614 
    ## Run 277 stress 0.01641859 
    ## Run 278 stress 0.01115306 
    ## ... Procrustes: rmse 0.0005623244  max resid 0.001269146 
    ## ... Similar to previous best
    ## Run 279 stress 0.01115188 
    ## ... Procrustes: rmse 0.000182308  max resid 0.0004116186 
    ## ... Similar to previous best
    ## Run 280 stress 0.01554263 
    ## Run 281 stress 0.01529066 
    ## Run 282 stress 0.01515768 
    ## Run 283 stress 0.0151497 
    ## Run 284 stress 0.01545378 
    ## Run 285 stress 0.01545355 
    ## Run 286 stress 0.01557885 
    ## Run 287 stress 0.01567103 
    ## Run 288 stress 0.01115164 
    ## ... Procrustes: rmse 0.00118522  max resid 0.002674933 
    ## ... Similar to previous best
    ## Run 289 stress 0.01115194 
    ## ... Procrustes: rmse 0.001307138  max resid 0.002951255 
    ## ... Similar to previous best
    ## Run 290 stress 0.0166349 
    ## Run 291 stress 0.0149858 
    ## Run 292 stress 0.01115841 
    ## ... Procrustes: rmse 0.001544472  max resid 0.003485451 
    ## ... Similar to previous best
    ## Run 293 stress 0.0111519 
    ## ... Procrustes: rmse 0.001290921  max resid 0.002914641 
    ## ... Similar to previous best
    ## Run 294 stress 0.01567947 
    ## Run 295 stress 0.01145381 
    ## ... Procrustes: rmse 0.01126389  max resid 0.02468372 
    ## Run 296 stress 0.01562167 
    ## Run 297 stress 0.01116495 
    ## ... Procrustes: rmse 0.002326285  max resid 0.005247481 
    ## ... Similar to previous best
    ## Run 298 stress 0.01562351 
    ## Run 299 stress 0.01670564 
    ## Run 300 stress 0.01556117 
    ## Run 301 stress 0.01120224 
    ## ... Procrustes: rmse 0.004910919  max resid 0.01106304 
    ## Run 302 stress 0.01384224 
    ## Run 303 stress 0.0124869 
    ## Run 304 stress 0.01527577 
    ## Run 305 stress 0.01498595 
    ## Run 306 stress 0.01557143 
    ## Run 307 stress 0.01549948 
    ## Run 308 stress 0.01528201 
    ## Run 309 stress 0.01529101 
    ## Run 310 stress 0.01515731 
    ## Run 311 stress 0.0151571 
    ## Run 312 stress 0.01547236 
    ## Run 313 stress 0.01498584 
    ## Run 314 stress 0.0155608 
    ## Run 315 stress 0.01129467 
    ## ... Procrustes: rmse 0.008632992  max resid 0.0193709 
    ## Run 316 stress 0.01557932 
    ## Run 317 stress 0.01543188 
    ## Run 318 stress 0.01528209 
    ## Run 319 stress 0.01648647 
    ## Run 320 stress 0.01533764 
    ## Run 321 stress 0.01122932 
    ## ... Procrustes: rmse 0.005583648  max resid 0.01258668 
    ## Run 322 stress 0.01515734 
    ## Run 323 stress 0.01637135 
    ## Run 324 stress 0.01515896 
    ## Run 325 stress 0.0157382 
    ## Run 326 stress 0.01562203 
    ## Run 327 stress 0.2711612 
    ## Run 328 stress 0.01537282 
    ## Run 329 stress 0.01544863 
    ## Run 330 stress 0.01558384 
    ## Run 331 stress 0.01556982 
    ## Run 332 stress 0.01124409 
    ## ... Procrustes: rmse 0.006879032  max resid 0.01548933 
    ## Run 333 stress 0.01528188 
    ## Run 334 stress 0.0156723 
    ## Run 335 stress 0.01557886 
    ## Run 336 stress 0.01115158 
    ## ... Procrustes: rmse 0.001155575  max resid 0.002608115 
    ## ... Similar to previous best
    ## Run 337 stress 0.0111519 
    ## ... Procrustes: rmse 0.001297304  max resid 0.002928883 
    ## ... Similar to previous best
    ## Run 338 stress 0.01649825 
    ## Run 339 stress 0.01115214 
    ## ... Procrustes: rmse 0.0002800811  max resid 0.0006323861 
    ## ... Similar to previous best
    ## Run 340 stress 0.01542913 
    ## Run 341 stress 0.01564574 
    ## Run 342 stress 0.01119453 
    ## ... Procrustes: rmse 0.004468722  max resid 0.01006997 
    ## Run 343 stress 0.01117367 
    ## ... Procrustes: rmse 0.003071657  max resid 0.006925216 
    ## ... Similar to previous best
    ## Run 344 stress 0.01576328 
    ## Run 345 stress 0.01499108 
    ## Run 346 stress 0.01499204 
    ## Run 347 stress 0.01498585 
    ## Run 348 stress 0.01115211 
    ## ... Procrustes: rmse 0.001374147  max resid 0.003102593 
    ## ... Similar to previous best
    ## Run 349 stress 0.0153379 
    ## Run 350 stress 0.0149859 
    ## Run 351 stress 0.01498603 
    ## Run 352 stress 0.01637122 
    ## Run 353 stress 0.01547441 
    ## Run 354 stress 0.01641868 
    ## Run 355 stress 0.01527569 
    ## Run 356 stress 0.01498565 
    ## Run 357 stress 0.0152579 
    ## Run 358 stress 0.01558237 
    ## Run 359 stress 0.01169131 
    ## Run 360 stress 0.01564563 
    ## Run 361 stress 0.01577125 
    ## Run 362 stress 0.01504126 
    ## Run 363 stress 0.01116063 
    ## ... Procrustes: rmse 0.001834036  max resid 0.004137158 
    ## ... Similar to previous best
    ## Run 364 stress 0.01582974 
    ## Run 365 stress 0.01115188 
    ## ... Procrustes: rmse 0.001292388  max resid 0.002917714 
    ## ... Similar to previous best
    ## Run 366 stress 0.01528222 
    ## Run 367 stress 0.0158304 
    ## Run 368 stress 0.01614177 
    ## Run 369 stress 0.01121714 
    ## ... Procrustes: rmse 0.005648735  max resid 0.01271782 
    ## Run 370 stress 0.01594321 
    ## Run 371 stress 0.01496318 
    ## Run 372 stress 0.01115188 
    ## ... Procrustes: rmse 0.00128394  max resid 0.002898404 
    ## ... Similar to previous best
    ## Run 373 stress 0.01115323 
    ## ... Procrustes: rmse 0.0006066391  max resid 0.00136937 
    ## ... Similar to previous best
    ## Run 374 stress 0.01115196 
    ## ... Procrustes: rmse 0.001323869  max resid 0.002989042 
    ## ... Similar to previous best
    ## Run 375 stress 0.01556079 
    ## Run 376 stress 0.01115173 
    ## ... Procrustes: rmse 0.001230232  max resid 0.002776888 
    ## ... Similar to previous best
    ## Run 377 stress 0.01115175 
    ## ... Procrustes: rmse 0.00123776  max resid 0.002793843 
    ## ... Similar to previous best
    ## Run 378 stress 0.01115178 
    ## ... Procrustes: rmse 0.001248972  max resid 0.002819349 
    ## ... Similar to previous best
    ## Run 379 stress 0.01556968 
    ## Run 380 stress 0.01576983 
    ## Run 381 stress 0.01115297 
    ## ... Procrustes: rmse 0.0005383888  max resid 0.001215255 
    ## ... Similar to previous best
    ## Run 382 stress 0.01649813 
    ## Run 383 stress 0.01660054 
    ## Run 384 stress 0.01649846 
    ## Run 385 stress 0.01121612 
    ## ... Procrustes: rmse 0.005661027  max resid 0.01274805 
    ## Run 386 stress 0.01557919 
    ## Run 387 stress 0.01545669 
    ## Run 388 stress 0.01536181 
    ## Run 389 stress 0.01528201 
    ## Run 390 stress 0.01115197 
    ## ... Procrustes: rmse 0.001316908  max resid 0.002972798 
    ## ... Similar to previous best
    ## Run 391 stress 0.01649844 
    ## Run 392 stress 0.01116865 
    ## ... Procrustes: rmse 0.002504554  max resid 0.005653307 
    ## ... Similar to previous best
    ## Run 393 stress 0.01498573 
    ## Run 394 stress 0.0151572 
    ## Run 395 stress 0.01557292 
    ## Run 396 stress 0.01154181 
    ## ... Procrustes: rmse 0.01373807  max resid 0.02959762 
    ## Run 397 stress 0.01567225 
    ## Run 398 stress 0.01515687 
    ## Run 399 stress 0.01547341 
    ## Run 400 stress 0.01679832 
    ## Run 401 stress 0.01536064 
    ## Run 402 stress 0.01546044 
    ## Run 403 stress 0.01533786 
    ## Run 404 stress 0.01648713 
    ## Run 405 stress 0.01515699 
    ## Run 406 stress 0.01515729 
    ## Run 407 stress 0.01582979 
    ## Run 408 stress 0.01547436 
    ## Run 409 stress 0.01528174 
    ## Run 410 stress 0.01499221 
    ## Run 411 stress 0.01583007 
    ## Run 412 stress 0.01555972 
    ## Run 413 stress 0.01556962 
    ## Run 414 stress 0.01677926 
    ## Run 415 stress 0.01665941 
    ## Run 416 stress 0.0155616 
    ## Run 417 stress 0.01515702 
    ## Run 418 stress 0.01556992 
    ## Run 419 stress 0.01545359 
    ## Run 420 stress 0.01498564 
    ## Run 421 stress 0.01515196 
    ## Run 422 stress 0.01115196 
    ## ... Procrustes: rmse 0.001322093  max resid 0.002985006 
    ## ... Similar to previous best
    ## Run 423 stress 0.01515706 
    ## Run 424 stress 0.3809694 
    ## Run 425 stress 0.01582269 
    ## Run 426 stress 0.01571544 
    ## Run 427 stress 0.01116344 
    ## ... Procrustes: rmse 0.002167471  max resid 0.004889623 
    ## ... Similar to previous best
    ## Run 428 stress 0.01576881 
    ## Run 429 stress 0.01574198 
    ## Run 430 stress 0.01546181 
    ## Run 431 stress 0.01498587 
    ## Run 432 stress 0.01648703 
    ## Run 433 stress 0.01493169 
    ## Run 434 stress 0.01549963 
    ## Run 435 stress 0.01515739 
    ## Run 436 stress 0.01546471 
    ## Run 437 stress 0.01115188 
    ## ... Procrustes: rmse 0.001290004  max resid 0.002912346 
    ## ... Similar to previous best
    ## Run 438 stress 0.01547237 
    ## Run 439 stress 0.01547879 
    ## Run 440 stress 0.01282443 
    ## Run 441 stress 0.01115192 
    ## ... Procrustes: rmse 0.0001918788  max resid 0.000433527 
    ## ... Similar to previous best
    ## Run 442 stress 0.01116311 
    ## ... Procrustes: rmse 0.002115846  max resid 0.004772155 
    ## ... Similar to previous best
    ## Run 443 stress 0.01493169 
    ## Run 444 stress 0.01529491 
    ## Run 445 stress 0.01564575 
    ## Run 446 stress 0.01115193 
    ## ... Procrustes: rmse 0.001310725  max resid 0.002959169 
    ## ... Similar to previous best
    ## Run 447 stress 0.01498573 
    ## Run 448 stress 0.01499245 
    ## Run 449 stress 0.01711397 
    ## Run 450 stress 0.01515701 
    ## Run 451 stress 0.0149922 
    ## Run 452 stress 0.01545072 
    ## Run 453 stress 0.01115304 
    ## ... Procrustes: rmse 0.0005572139  max resid 0.001257801 
    ## ... Similar to previous best
    ## Run 454 stress 0.01669216 
    ## Run 455 stress 0.01545288 
    ## Run 456 stress 0.01666971 
    ## Run 457 stress 0.01536184 
    ## Run 458 stress 0.0154321 
    ## Run 459 stress 0.0155693 
    ## Run 460 stress 0.01135942 
    ## ... Procrustes: rmse 0.01036564  max resid 0.02314017 
    ## Run 461 stress 0.01554869 
    ## Run 462 stress 0.01515859 
    ## Run 463 stress 0.01498605 
    ## Run 464 stress 0.01663496 
    ## Run 465 stress 0.01545338 
    ## Run 466 stress 0.01548125 
    ## Run 467 stress 0.01564595 
    ## Run 468 stress 0.0168681 
    ## Run 469 stress 0.01546488 
    ## Run 470 stress 0.01556146 
    ## Run 471 stress 0.01119351 
    ## ... Procrustes: rmse 0.004466501  max resid 0.01006988 
    ## Run 472 stress 0.01185521 
    ## Run 473 stress 0.01196246 
    ## Run 474 stress 0.01115329 
    ## ... Procrustes: rmse 0.0006101789  max resid 0.001377779 
    ## ... Similar to previous best
    ## Run 475 stress 0.01115195 
    ## ... Procrustes: rmse 0.001318218  max resid 0.002976244 
    ## ... Similar to previous best
    ## Run 476 stress 0.01549943 
    ## Run 477 stress 0.01555698 
    ## Run 478 stress 0.01493166 
    ## Run 479 stress 0.01545346 
    ## Run 480 stress 0.01564827 
    ## Run 481 stress 0.01694536 
    ## Run 482 stress 0.01137115 
    ## ... Procrustes: rmse 0.009536715  max resid 0.02135713 
    ## Run 483 stress 0.01498594 
    ## Run 484 stress 0.01527904 
    ## Run 485 stress 0.01576307 
    ## Run 486 stress 0.01172655 
    ## Run 487 stress 0.01115201 
    ## ... Procrustes: rmse 0.0002287565  max resid 0.0005167867 
    ## ... Similar to previous best
    ## Run 488 stress 0.01525775 
    ## Run 489 stress 0.01549939 
    ## Run 490 stress 0.0155367 
    ## Run 491 stress 0.01149918 
    ## ... Procrustes: rmse 0.01245464  max resid 0.02695796 
    ## Run 492 stress 0.01685388 
    ## Run 493 stress 0.01550897 
    ## Run 494 stress 0.01528198 
    ## Run 495 stress 0.01545332 
    ## Run 496 stress 0.01562925 
    ## Run 497 stress 0.01564557 
    ## Run 498 stress 0.01115188 
    ## ... Procrustes: rmse 0.001290593  max resid 0.002913642 
    ## ... Similar to previous best
    ## Run 499 stress 0.01151167 
    ## ... Procrustes: rmse 0.01256156  max resid 0.02716891 
    ## Run 500 stress 0.01533798 
    ## *** Best solution repeated 50 times

``` r
# Mixed and stratified lakes
PD_beta_MS_NMDS <- metaMDS(PD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.05171059 
    ## Run 2 stress 0.05466703 
    ## Run 3 stress 0.05403606 
    ## Run 4 stress 0.05968643 
    ## Run 5 stress 0.05900149 
    ## Run 6 stress 0.05036794 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002175539  max resid 0.0005532499 
    ## ... Similar to previous best
    ## Run 7 stress 0.06362778 
    ## Run 8 stress 0.05337578 
    ## Run 9 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001601898  max resid 0.0004082782 
    ## ... Similar to previous best
    ## Run 10 stress 0.0517106 
    ## Run 11 stress 0.05403617 
    ## Run 12 stress 0.05337577 
    ## Run 13 stress 0.05829931 
    ## Run 14 stress 0.05036792 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001189267  max resid 0.0002856475 
    ## ... Similar to previous best
    ## Run 15 stress 0.05466701 
    ## Run 16 stress 0.06655724 
    ## Run 17 stress 0.05171059 
    ## Run 18 stress 0.05051332 
    ## ... Procrustes: rmse 0.03574779  max resid 0.1221322 
    ## Run 19 stress 0.05051305 
    ## ... Procrustes: rmse 0.03572644  max resid 0.1219358 
    ## Run 20 stress 0.05337579 
    ## Run 21 stress 0.05602049 
    ## Run 22 stress 0.06826693 
    ## Run 23 stress 0.05171059 
    ## Run 24 stress 0.05036796 
    ## ... Procrustes: rmse 2.198215e-05  max resid 5.297272e-05 
    ## ... Similar to previous best
    ## Run 25 stress 0.06476238 
    ## Run 26 stress 0.05036799 
    ## ... Procrustes: rmse 7.839575e-05  max resid 0.0002107381 
    ## ... Similar to previous best
    ## Run 27 stress 0.05979236 
    ## Run 28 stress 0.05051304 
    ## ... Procrustes: rmse 0.03572679  max resid 0.1220586 
    ## Run 29 stress 0.05036792 
    ## ... New best solution
    ## ... Procrustes: rmse 8.033181e-06  max resid 1.536379e-05 
    ## ... Similar to previous best
    ## Run 30 stress 0.05051349 
    ## ... Procrustes: rmse 0.03574  max resid 0.1219442 
    ## Run 31 stress 0.05171064 
    ## Run 32 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001937897  max resid 0.0004532081 
    ## ... Similar to previous best
    ## Run 33 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001525087  max resid 0.0003728733 
    ## ... Similar to previous best
    ## Run 34 stress 0.0533758 
    ## Run 35 stress 0.06743836 
    ## Run 36 stress 0.05171058 
    ## Run 37 stress 0.05602054 
    ## Run 38 stress 0.05036798 
    ## ... Procrustes: rmse 7.660083e-05  max resid 0.0002206096 
    ## ... Similar to previous best
    ## Run 39 stress 0.05928219 
    ## Run 40 stress 0.05171058 
    ## Run 41 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001736833  max resid 0.0003577561 
    ## ... Similar to previous best
    ## Run 42 stress 0.050368 
    ## ... Procrustes: rmse 0.0001711096  max resid 0.0004217409 
    ## ... Similar to previous best
    ## Run 43 stress 0.06845227 
    ## Run 44 stress 0.05051353 
    ## ... Procrustes: rmse 0.03576413  max resid 0.1221868 
    ## Run 45 stress 0.05337589 
    ## Run 46 stress 0.0505133 
    ## ... Procrustes: rmse 0.03571131  max resid 0.1218632 
    ## Run 47 stress 0.06609073 
    ## Run 48 stress 0.365735 
    ## Run 49 stress 0.05036794 
    ## ... Procrustes: rmse 0.0001425305  max resid 0.0003368169 
    ## ... Similar to previous best
    ## Run 50 stress 0.05979228 
    ## Run 51 stress 0.05171061 
    ## Run 52 stress 0.05171069 
    ## Run 53 stress 0.05171059 
    ## Run 54 stress 0.05171059 
    ## Run 55 stress 0.05337583 
    ## Run 56 stress 0.05171065 
    ## Run 57 stress 0.05051319 
    ## ... Procrustes: rmse 0.03571434  max resid 0.1220277 
    ## Run 58 stress 0.050368 
    ## ... Procrustes: rmse 0.0001921611  max resid 0.0004699539 
    ## ... Similar to previous best
    ## Run 59 stress 0.05337592 
    ## Run 60 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001318113  max resid 0.0003182515 
    ## ... Similar to previous best
    ## Run 61 stress 0.05051341 
    ## ... Procrustes: rmse 0.03575924  max resid 0.1221705 
    ## Run 62 stress 0.06136679 
    ## Run 63 stress 0.0560207 
    ## Run 64 stress 0.05051299 
    ## ... Procrustes: rmse 0.03571832  max resid 0.1220283 
    ## Run 65 stress 0.05829932 
    ## Run 66 stress 0.062313 
    ## Run 67 stress 0.05051307 
    ## ... Procrustes: rmse 0.03572964  max resid 0.1220682 
    ## Run 68 stress 0.05466701 
    ## Run 69 stress 0.05171062 
    ## Run 70 stress 0.05337579 
    ## Run 71 stress 0.05466701 
    ## Run 72 stress 0.05466714 
    ## Run 73 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001751891  max resid 0.0004014476 
    ## ... Similar to previous best
    ## Run 74 stress 0.05036793 
    ## ... Procrustes: rmse 2.117606e-05  max resid 5.78606e-05 
    ## ... Similar to previous best
    ## Run 75 stress 0.06309836 
    ## Run 76 stress 0.05403593 
    ## Run 77 stress 0.05051337 
    ## ... Procrustes: rmse 0.035755  max resid 0.1221563 
    ## Run 78 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 7.875698e-05  max resid 0.0001933801 
    ## ... Similar to previous best
    ## Run 79 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001306374  max resid 0.0003828958 
    ## ... Similar to previous best
    ## Run 80 stress 0.05171064 
    ## Run 81 stress 0.05051338 
    ## ... Procrustes: rmse 0.03568711  max resid 0.1219401 
    ## Run 82 stress 0.05829932 
    ## Run 83 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001337376  max resid 0.0003779968 
    ## ... Similar to previous best
    ## Run 84 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001342319  max resid 0.0003378474 
    ## ... Similar to previous best
    ## Run 85 stress 0.06409281 
    ## Run 86 stress 0.05171059 
    ## Run 87 stress 0.05466704 
    ## Run 88 stress 0.05928242 
    ## Run 89 stress 0.06362785 
    ## Run 90 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001411524  max resid 0.000391851 
    ## ... Similar to previous best
    ## Run 91 stress 0.05337577 
    ## Run 92 stress 0.05403581 
    ## Run 93 stress 0.06111766 
    ## Run 94 stress 0.06950981 
    ## Run 95 stress 0.05602043 
    ## Run 96 stress 0.05171061 
    ## Run 97 stress 0.05602074 
    ## Run 98 stress 0.06730251 
    ## Run 99 stress 0.06136677 
    ## Run 100 stress 0.0540359 
    ## Run 101 stress 0.05171066 
    ## Run 102 stress 0.05403569 
    ## Run 103 stress 0.05829932 
    ## Run 104 stress 0.06309854 
    ## Run 105 stress 0.05979225 
    ## Run 106 stress 0.05051347 
    ## ... Procrustes: rmse 0.03573081  max resid 0.1220738 
    ## Run 107 stress 0.05602056 
    ## Run 108 stress 0.06499236 
    ## Run 109 stress 0.05968634 
    ## Run 110 stress 0.0533759 
    ## Run 111 stress 0.0517106 
    ## Run 112 stress 0.06231302 
    ## Run 113 stress 0.05602088 
    ## Run 114 stress 0.06136685 
    ## Run 115 stress 0.06826704 
    ## Run 116 stress 0.05337591 
    ## Run 117 stress 0.05979229 
    ## Run 118 stress 0.05968634 
    ## Run 119 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001110445  max resid 0.0002707657 
    ## ... Similar to previous best
    ## Run 120 stress 0.05171064 
    ## Run 121 stress 0.05829934 
    ## Run 122 stress 0.05968627 
    ## Run 123 stress 0.06231287 
    ## Run 124 stress 0.05337579 
    ## Run 125 stress 0.05051359 
    ## ... Procrustes: rmse 0.03570857  max resid 0.1218141 
    ## Run 126 stress 0.05051325 
    ## ... Procrustes: rmse 0.03565664  max resid 0.1216849 
    ## Run 127 stress 0.05051342 
    ## ... Procrustes: rmse 0.03568772  max resid 0.1219433 
    ## Run 128 stress 0.05466703 
    ## Run 129 stress 0.05829942 
    ## Run 130 stress 0.05171062 
    ## Run 131 stress 0.05979228 
    ## Run 132 stress 0.06362794 
    ## Run 133 stress 0.05036796 
    ## ... Procrustes: rmse 8.655473e-05  max resid 0.0002141631 
    ## ... Similar to previous best
    ## Run 134 stress 0.0505133 
    ## ... Procrustes: rmse 0.03565718  max resid 0.1216893 
    ## Run 135 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001519101  max resid 0.0003954179 
    ## ... Similar to previous best
    ## Run 136 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 2.54073e-05  max resid 7.089861e-05 
    ## ... Similar to previous best
    ## Run 137 stress 0.05051331 
    ## ... Procrustes: rmse 0.03565742  max resid 0.121688 
    ## Run 138 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001530957  max resid 0.0003838424 
    ## ... Similar to previous best
    ## Run 139 stress 0.05337579 
    ## Run 140 stress 0.05051325 
    ## ... Procrustes: rmse 0.03569054  max resid 0.1217907 
    ## Run 141 stress 0.05051323 
    ## ... Procrustes: rmse 0.03569114  max resid 0.1219509 
    ## Run 142 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001225672  max resid 0.0003006884 
    ## ... Similar to previous best
    ## Run 143 stress 0.05403573 
    ## Run 144 stress 0.0582993 
    ## Run 145 stress 0.06231302 
    ## Run 146 stress 0.05036798 
    ## ... Procrustes: rmse 0.000122851  max resid 0.0003180058 
    ## ... Similar to previous best
    ## Run 147 stress 0.05602051 
    ## Run 148 stress 0.05051301 
    ## ... Procrustes: rmse 0.03581053  max resid 0.1222839 
    ## Run 149 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001552087  max resid 0.0003588982 
    ## ... Similar to previous best
    ## Run 150 stress 0.05602047 
    ## Run 151 stress 0.06309844 
    ## Run 152 stress 0.05337578 
    ## Run 153 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001415964  max resid 0.0003861105 
    ## ... Similar to previous best
    ## Run 154 stress 0.06231288 
    ## Run 155 stress 0.06309848 
    ## Run 156 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.624107e-05  max resid 3.230942e-05 
    ## ... Similar to previous best
    ## Run 157 stress 0.05466704 
    ## Run 158 stress 0.05051347 
    ## ... Procrustes: rmse 0.03571916  max resid 0.1220419 
    ## Run 159 stress 0.05036794 
    ## ... Procrustes: rmse 8.772499e-05  max resid 0.0002062369 
    ## ... Similar to previous best
    ## Run 160 stress 0.06716796 
    ## Run 161 stress 0.06362799 
    ## Run 162 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001410151  max resid 0.0003305578 
    ## ... Similar to previous best
    ## Run 163 stress 0.06136693 
    ## Run 164 stress 0.06309857 
    ## Run 165 stress 0.05036792 
    ## ... Procrustes: rmse 1.283495e-05  max resid 2.846495e-05 
    ## ... Similar to previous best
    ## Run 166 stress 0.06111764 
    ## Run 167 stress 0.05051343 
    ## ... Procrustes: rmse 0.03573809  max resid 0.1219233 
    ## Run 168 stress 0.05051337 
    ## ... Procrustes: rmse 0.03570356  max resid 0.1219934 
    ## Run 169 stress 0.06038942 
    ## Run 170 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001440552  max resid 0.0003374171 
    ## ... Similar to previous best
    ## Run 171 stress 0.0505129 
    ## ... Procrustes: rmse 0.03569302  max resid 0.1218414 
    ## Run 172 stress 0.05051305 
    ## ... Procrustes: rmse 0.03570851  max resid 0.1219943 
    ## Run 173 stress 0.05171063 
    ## Run 174 stress 0.05036792 
    ## ... Procrustes: rmse 2.931127e-05  max resid 7.836762e-05 
    ## ... Similar to previous best
    ## Run 175 stress 0.06478515 
    ## Run 176 stress 0.05928264 
    ## Run 177 stress 0.05171061 
    ## Run 178 stress 0.05829932 
    ## Run 179 stress 0.05968626 
    ## Run 180 stress 0.05171062 
    ## Run 181 stress 0.05036806 
    ## ... Procrustes: rmse 0.0001800461  max resid 0.0005262529 
    ## ... Similar to previous best
    ## Run 182 stress 0.06136677 
    ## Run 183 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001322607  max resid 0.0003195766 
    ## ... Similar to previous best
    ## Run 184 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001394806  max resid 0.0003365585 
    ## ... Similar to previous best
    ## Run 185 stress 0.05036793 
    ## ... Procrustes: rmse 8.089141e-05  max resid 0.0001868517 
    ## ... Similar to previous best
    ## Run 186 stress 0.05337583 
    ## Run 187 stress 0.06499275 
    ## Run 188 stress 0.05036795 
    ## ... Procrustes: rmse 9.978868e-05  max resid 0.0002161271 
    ## ... Similar to previous best
    ## Run 189 stress 0.05036795 
    ## ... Procrustes: rmse 9.688607e-05  max resid 0.0002274583 
    ## ... Similar to previous best
    ## Run 190 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001598714  max resid 0.0004500634 
    ## ... Similar to previous best
    ## Run 191 stress 0.05051335 
    ## ... Procrustes: rmse 0.03573264  max resid 0.1220765 
    ## Run 192 stress 0.05051307 
    ## ... Procrustes: rmse 0.03570684  max resid 0.1219893 
    ## Run 193 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001617744  max resid 0.0003876228 
    ## ... Similar to previous best
    ## Run 194 stress 0.06309852 
    ## Run 195 stress 0.05928219 
    ## Run 196 stress 0.05036805 
    ## ... Procrustes: rmse 0.000175683  max resid 0.0004366716 
    ## ... Similar to previous best
    ## Run 197 stress 0.05171059 
    ## Run 198 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001685338  max resid 0.0004508449 
    ## ... Similar to previous best
    ## Run 199 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001449233  max resid 0.0003849734 
    ## ... Similar to previous best
    ## Run 200 stress 0.0505133 
    ## ... Procrustes: rmse 0.03574686  max resid 0.1221187 
    ## Run 201 stress 0.05171058 
    ## Run 202 stress 0.05036791 
    ## ... Procrustes: rmse 4.06925e-05  max resid 8.084264e-05 
    ## ... Similar to previous best
    ## Run 203 stress 0.05928218 
    ## Run 204 stress 0.06640348 
    ## Run 205 stress 0.05968632 
    ## Run 206 stress 0.0582993 
    ## Run 207 stress 0.05403598 
    ## Run 208 stress 0.05036794 
    ## ... Procrustes: rmse 9.042097e-05  max resid 0.0002138995 
    ## ... Similar to previous best
    ## Run 209 stress 0.05900151 
    ## Run 210 stress 0.05968628 
    ## Run 211 stress 0.05051345 
    ## ... Procrustes: rmse 0.03574911  max resid 0.1221323 
    ## Run 212 stress 0.05602063 
    ## Run 213 stress 0.05968629 
    ## Run 214 stress 0.06231297 
    ## Run 215 stress 0.05928216 
    ## Run 216 stress 0.05602057 
    ## Run 217 stress 0.05968628 
    ## Run 218 stress 0.05602057 
    ## Run 219 stress 0.05403584 
    ## Run 220 stress 0.0505132 
    ## ... Procrustes: rmse 0.03570666  max resid 0.1218489 
    ## Run 221 stress 0.05051289 
    ## ... Procrustes: rmse 0.03569765  max resid 0.1219494 
    ## Run 222 stress 0.05829935 
    ## Run 223 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001466505  max resid 0.0004055619 
    ## ... Similar to previous best
    ## Run 224 stress 0.05036794 
    ## ... Procrustes: rmse 7.581715e-05  max resid 0.0001744249 
    ## ... Similar to previous best
    ## Run 225 stress 0.3312177 
    ## Run 226 stress 0.05928221 
    ## Run 227 stress 0.05403599 
    ## Run 228 stress 0.0688511 
    ## Run 229 stress 0.05403571 
    ## Run 230 stress 0.06231295 
    ## Run 231 stress 0.06478525 
    ## Run 232 stress 0.06309856 
    ## Run 233 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001535267  max resid 0.0003699005 
    ## ... Similar to previous best
    ## Run 234 stress 0.05171059 
    ## Run 235 stress 0.06111766 
    ## Run 236 stress 0.05171058 
    ## Run 237 stress 0.05968633 
    ## Run 238 stress 0.05036794 
    ## ... Procrustes: rmse 7.247607e-05  max resid 0.0001722907 
    ## ... Similar to previous best
    ## Run 239 stress 0.05051298 
    ## ... Procrustes: rmse 0.03567751  max resid 0.1217833 
    ## Run 240 stress 0.06478538 
    ## Run 241 stress 0.06943819 
    ## Run 242 stress 0.05171063 
    ## Run 243 stress 0.05466702 
    ## Run 244 stress 0.05171058 
    ## Run 245 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001579992  max resid 0.0003821617 
    ## ... Similar to previous best
    ## Run 246 stress 0.05337579 
    ## Run 247 stress 0.05337586 
    ## Run 248 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001509992  max resid 0.0003905053 
    ## ... Similar to previous best
    ## Run 249 stress 0.05928216 
    ## Run 250 stress 0.05900148 
    ## Run 251 stress 0.0540357 
    ## Run 252 stress 0.05036793 
    ## ... Procrustes: rmse 7.902599e-05  max resid 0.0001791365 
    ## ... Similar to previous best
    ## Run 253 stress 0.05051349 
    ## ... Procrustes: rmse 0.03574979  max resid 0.1221351 
    ## Run 254 stress 0.06309837 
    ## Run 255 stress 0.050368 
    ## ... Procrustes: rmse 0.0001438931  max resid 0.0004012551 
    ## ... Similar to previous best
    ## Run 256 stress 0.05171066 
    ## Run 257 stress 0.06136675 
    ## Run 258 stress 0.05171059 
    ## Run 259 stress 0.05337577 
    ## Run 260 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001040546  max resid 0.0002316836 
    ## ... Similar to previous best
    ## Run 261 stress 0.05036796 
    ## ... Procrustes: rmse 8.713128e-05  max resid 0.0002409782 
    ## ... Similar to previous best
    ## Run 262 stress 0.06476217 
    ## Run 263 stress 0.0620303 
    ## Run 264 stress 0.06136681 
    ## Run 265 stress 0.05602051 
    ## Run 266 stress 0.05403583 
    ## Run 267 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001425582  max resid 0.000325168 
    ## ... Similar to previous best
    ## Run 268 stress 0.06609082 
    ## Run 269 stress 0.05979225 
    ## Run 270 stress 0.05036809 
    ## ... Procrustes: rmse 0.0001831451  max resid 0.0004045496 
    ## ... Similar to previous best
    ## Run 271 stress 0.06476233 
    ## Run 272 stress 0.06630524 
    ## Run 273 stress 0.05829929 
    ## Run 274 stress 0.05171058 
    ## Run 275 stress 0.05979229 
    ## Run 276 stress 0.05036793 
    ## ... Procrustes: rmse 7.095901e-05  max resid 0.0002072042 
    ## ... Similar to previous best
    ## Run 277 stress 0.05403603 
    ## Run 278 stress 0.05036791 
    ## ... Procrustes: rmse 8.485603e-06  max resid 2.456227e-05 
    ## ... Similar to previous best
    ## Run 279 stress 0.05051318 
    ## ... Procrustes: rmse 0.03569704  max resid 0.1219663 
    ## Run 280 stress 0.05036795 
    ## ... Procrustes: rmse 9.454039e-05  max resid 0.0002549992 
    ## ... Similar to previous best
    ## Run 281 stress 0.05829936 
    ## Run 282 stress 0.06027433 
    ## Run 283 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001057213  max resid 0.0002456476 
    ## ... Similar to previous best
    ## Run 284 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001273185  max resid 0.0003070854 
    ## ... Similar to previous best
    ## Run 285 stress 0.050513 
    ## ... Procrustes: rmse 0.03576295  max resid 0.1220429 
    ## Run 286 stress 0.0517106 
    ## Run 287 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001052024  max resid 0.0003067416 
    ## ... Similar to previous best
    ## Run 288 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001404841  max resid 0.0003398198 
    ## ... Similar to previous best
    ## Run 289 stress 0.05928222 
    ## Run 290 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001772927  max resid 0.0004328264 
    ## ... Similar to previous best
    ## Run 291 stress 0.06231315 
    ## Run 292 stress 0.0505134 
    ## ... Procrustes: rmse 0.03574441  max resid 0.1221159 
    ## Run 293 stress 0.05602072 
    ## Run 294 stress 0.05968626 
    ## Run 295 stress 0.05051333 
    ## ... Procrustes: rmse 0.03573251  max resid 0.1220787 
    ## Run 296 stress 0.05036795 
    ## ... Procrustes: rmse 9.063927e-05  max resid 0.0002500568 
    ## ... Similar to previous best
    ## Run 297 stress 0.05036793 
    ## ... Procrustes: rmse 3.784485e-05  max resid 9.408714e-05 
    ## ... Similar to previous best
    ## Run 298 stress 0.05036797 
    ## ... Procrustes: rmse 9.94947e-05  max resid 0.0001999011 
    ## ... Similar to previous best
    ## Run 299 stress 0.05466703 
    ## Run 300 stress 0.06409302 
    ## Run 301 stress 0.05036793 
    ## ... Procrustes: rmse 7.253338e-05  max resid 0.0001639618 
    ## ... Similar to previous best
    ## Run 302 stress 0.0582994 
    ## Run 303 stress 0.05051314 
    ## ... Procrustes: rmse 0.03571518  max resid 0.1220189 
    ## Run 304 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001908756  max resid 0.0004664193 
    ## ... Similar to previous best
    ## Run 305 stress 0.05337577 
    ## Run 306 stress 0.3658208 
    ## Run 307 stress 0.05036795 
    ## ... Procrustes: rmse 7.756097e-05  max resid 0.0002011095 
    ## ... Similar to previous best
    ## Run 308 stress 0.05829936 
    ## Run 309 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001186461  max resid 0.0002622659 
    ## ... Similar to previous best
    ## Run 310 stress 0.05979225 
    ## Run 311 stress 0.06027445 
    ## Run 312 stress 0.05171067 
    ## Run 313 stress 0.05466702 
    ## Run 314 stress 0.05051325 
    ## ... Procrustes: rmse 0.0357339  max resid 0.1220799 
    ## Run 315 stress 0.05403572 
    ## Run 316 stress 0.05051325 
    ## ... Procrustes: rmse 0.03564566  max resid 0.1216603 
    ## Run 317 stress 0.05829935 
    ## Run 318 stress 0.05051345 
    ## ... Procrustes: rmse 0.03574893  max resid 0.1221318 
    ## Run 319 stress 0.05466703 
    ## Run 320 stress 0.06621981 
    ## Run 321 stress 0.06231287 
    ## Run 322 stress 0.05829931 
    ## Run 323 stress 0.05171058 
    ## Run 324 stress 0.06309851 
    ## Run 325 stress 0.06476244 
    ## Run 326 stress 0.05051334 
    ## ... Procrustes: rmse 0.03565064  max resid 0.1216656 
    ## Run 327 stress 0.0503681 
    ## ... Procrustes: rmse 0.0002161648  max resid 0.0004700076 
    ## ... Similar to previous best
    ## Run 328 stress 0.05036793 
    ## ... Procrustes: rmse 6.902899e-05  max resid 0.0001595039 
    ## ... Similar to previous best
    ## Run 329 stress 0.0611178 
    ## Run 330 stress 0.06478524 
    ## Run 331 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001178851  max resid 0.0003396006 
    ## ... Similar to previous best
    ## Run 332 stress 0.05051315 
    ## ... Procrustes: rmse 0.03569605  max resid 0.1219621 
    ## Run 333 stress 0.05337584 
    ## Run 334 stress 0.06136674 
    ## Run 335 stress 0.06537738 
    ## Run 336 stress 0.05051338 
    ## ... Procrustes: rmse 0.0357152  max resid 0.122029 
    ## Run 337 stress 0.0623132 
    ## Run 338 stress 0.05602056 
    ## Run 339 stress 0.0505131 
    ## ... Procrustes: rmse 0.03574503  max resid 0.1219739 
    ## Run 340 stress 0.05466707 
    ## Run 341 stress 0.05979225 
    ## Run 342 stress 0.05036791 
    ## ... Procrustes: rmse 3.054675e-05  max resid 6.319014e-05 
    ## ... Similar to previous best
    ## Run 343 stress 0.07221107 
    ## Run 344 stress 0.05036792 
    ## ... Procrustes: rmse 4.561775e-05  max resid 0.0001256751 
    ## ... Similar to previous best
    ## Run 345 stress 0.05051308 
    ## ... Procrustes: rmse 0.03568225  max resid 0.1217855 
    ## Run 346 stress 0.05403572 
    ## Run 347 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001580171  max resid 0.0003861859 
    ## ... Similar to previous best
    ## Run 348 stress 0.05036792 
    ## ... Procrustes: rmse 4.958777e-05  max resid 0.0001303951 
    ## ... Similar to previous best
    ## Run 349 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001060361  max resid 0.000249366 
    ## ... Similar to previous best
    ## Run 350 stress 0.05171059 
    ## Run 351 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001477376  max resid 0.0004077855 
    ## ... Similar to previous best
    ## Run 352 stress 0.06493695 
    ## Run 353 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001391538  max resid 0.0003886864 
    ## ... Similar to previous best
    ## Run 354 stress 0.05171063 
    ## Run 355 stress 0.05051327 
    ## ... Procrustes: rmse 0.03566374  max resid 0.1217153 
    ## Run 356 stress 0.05337589 
    ## Run 357 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001230872  max resid 0.000339854 
    ## ... Similar to previous best
    ## Run 358 stress 0.05829941 
    ## Run 359 stress 0.06493701 
    ## Run 360 stress 0.05968632 
    ## Run 361 stress 0.05466703 
    ## Run 362 stress 0.05829935 
    ## Run 363 stress 0.05968632 
    ## Run 364 stress 0.05466703 
    ## Run 365 stress 0.05337589 
    ## Run 366 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001753333  max resid 0.000481503 
    ## ... Similar to previous best
    ## Run 367 stress 0.05171063 
    ## Run 368 stress 0.06231305 
    ## Run 369 stress 0.05036795 
    ## ... Procrustes: rmse 9.723954e-05  max resid 0.000228641 
    ## ... Similar to previous best
    ## Run 370 stress 0.05171058 
    ## Run 371 stress 0.05036793 
    ## ... Procrustes: rmse 6.35547e-05  max resid 0.0001622965 
    ## ... Similar to previous best
    ## Run 372 stress 0.05928234 
    ## Run 373 stress 0.05036795 
    ## ... Procrustes: rmse 8.342655e-05  max resid 0.0002365218 
    ## ... Similar to previous best
    ## Run 374 stress 0.06111767 
    ## Run 375 stress 0.0596864 
    ## Run 376 stress 0.06231307 
    ## Run 377 stress 0.05979225 
    ## Run 378 stress 0.06476214 
    ## Run 379 stress 0.05968634 
    ## Run 380 stress 0.05968627 
    ## Run 381 stress 0.05466702 
    ## Run 382 stress 0.05171059 
    ## Run 383 stress 0.05466707 
    ## Run 384 stress 0.05602052 
    ## Run 385 stress 0.05337576 
    ## Run 386 stress 0.05337591 
    ## Run 387 stress 0.05602055 
    ## Run 388 stress 0.06136685 
    ## Run 389 stress 0.05036792 
    ## ... Procrustes: rmse 4.434059e-05  max resid 9.364877e-05 
    ## ... Similar to previous best
    ## Run 390 stress 0.05171064 
    ## Run 391 stress 0.0582993 
    ## Run 392 stress 0.05337577 
    ## Run 393 stress 0.05337583 
    ## Run 394 stress 0.05968632 
    ## Run 395 stress 0.05968637 
    ## Run 396 stress 0.06609052 
    ## Run 397 stress 0.06362778 
    ## Run 398 stress 0.3594412 
    ## Run 399 stress 0.05036795 
    ## ... Procrustes: rmse 9.900239e-05  max resid 0.0002327377 
    ## ... Similar to previous best
    ## Run 400 stress 0.05051273 
    ## ... Procrustes: rmse 0.0356688  max resid 0.1218378 
    ## Run 401 stress 0.05171058 
    ## Run 402 stress 0.05051284 
    ## ... Procrustes: rmse 0.03570708  max resid 0.1218921 
    ## Run 403 stress 0.06027433 
    ## Run 404 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001561379  max resid 0.0004494389 
    ## ... Similar to previous best
    ## Run 405 stress 0.06493696 
    ## Run 406 stress 0.0546671 
    ## Run 407 stress 0.05337584 
    ## Run 408 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001161084  max resid 0.0002702475 
    ## ... Similar to previous best
    ## Run 409 stress 0.05036798 
    ## ... Procrustes: rmse 8.201004e-05  max resid 0.0001804617 
    ## ... Similar to previous best
    ## Run 410 stress 0.05466705 
    ## Run 411 stress 0.06027439 
    ## Run 412 stress 0.05051318 
    ## ... Procrustes: rmse 0.03571876  max resid 0.1220313 
    ## Run 413 stress 0.06362795 
    ## Run 414 stress 0.05051325 
    ## ... Procrustes: rmse 0.03565657  max resid 0.121692 
    ## Run 415 stress 0.06478519 
    ## Run 416 stress 0.05979225 
    ## Run 417 stress 0.05602048 
    ## Run 418 stress 0.05602047 
    ## Run 419 stress 0.05036791 
    ## ... Procrustes: rmse 3.076305e-05  max resid 8.70901e-05 
    ## ... Similar to previous best
    ## Run 420 stress 0.06038956 
    ## Run 421 stress 0.05036793 
    ## ... Procrustes: rmse 7.533952e-05  max resid 0.0001697883 
    ## ... Similar to previous best
    ## Run 422 stress 0.0592822 
    ## Run 423 stress 0.05829936 
    ## Run 424 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001793944  max resid 0.0005059102 
    ## ... Similar to previous best
    ## Run 425 stress 0.06844316 
    ## Run 426 stress 0.06743843 
    ## Run 427 stress 0.05403619 
    ## Run 428 stress 0.05829949 
    ## Run 429 stress 0.05051327 
    ## ... Procrustes: rmse 0.03572193  max resid 0.1220446 
    ## Run 430 stress 0.05337594 
    ## Run 431 stress 0.06950966 
    ## Run 432 stress 0.05036792 
    ## ... Procrustes: rmse 5.556963e-05  max resid 0.0001235002 
    ## ... Similar to previous best
    ## Run 433 stress 0.05829932 
    ## Run 434 stress 0.05337578 
    ## Run 435 stress 0.06309849 
    ## Run 436 stress 0.05036794 
    ## ... Procrustes: rmse 8.714363e-05  max resid 0.0002028488 
    ## ... Similar to previous best
    ## Run 437 stress 0.05051335 
    ## ... Procrustes: rmse 0.03572467  max resid 0.122056 
    ## Run 438 stress 0.06309863 
    ## Run 439 stress 0.05928225 
    ## Run 440 stress 0.06038948 
    ## Run 441 stress 0.05051304 
    ## ... Procrustes: rmse 0.03567682  max resid 0.1217737 
    ## Run 442 stress 0.05403561 
    ## Run 443 stress 0.05829937 
    ## Run 444 stress 0.05968634 
    ## Run 445 stress 0.05968627 
    ## Run 446 stress 0.05968633 
    ## Run 447 stress 0.0582994 
    ## Run 448 stress 0.05928246 
    ## Run 449 stress 0.05829944 
    ## Run 450 stress 0.05036794 
    ## ... Procrustes: rmse 9.146709e-05  max resid 0.0002033554 
    ## ... Similar to previous best
    ## Run 451 stress 0.05928217 
    ## Run 452 stress 0.06136677 
    ## Run 453 stress 0.05036792 
    ## ... Procrustes: rmse 6.361173e-05  max resid 0.0001374572 
    ## ... Similar to previous best
    ## Run 454 stress 0.05036791 
    ## ... Procrustes: rmse 2.760171e-05  max resid 7.832958e-05 
    ## ... Similar to previous best
    ## Run 455 stress 0.05036793 
    ## ... Procrustes: rmse 2.852625e-05  max resid 7.205806e-05 
    ## ... Similar to previous best
    ## Run 456 stress 0.05466705 
    ## Run 457 stress 0.05036806 
    ## ... Procrustes: rmse 0.0001816176  max resid 0.0004906767 
    ## ... Similar to previous best
    ## Run 458 stress 0.0533758 
    ## Run 459 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001128778  max resid 0.0002200527 
    ## ... Similar to previous best
    ## Run 460 stress 0.05968627 
    ## Run 461 stress 0.06309878 
    ## Run 462 stress 0.062313 
    ## Run 463 stress 0.06478527 
    ## Run 464 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001576543  max resid 0.0004290249 
    ## ... Similar to previous best
    ## Run 465 stress 0.05602045 
    ## Run 466 stress 0.05968629 
    ## Run 467 stress 0.0505132 
    ## ... Procrustes: rmse 0.03579943  max resid 0.1221412 
    ## Run 468 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001857019  max resid 0.0004134382 
    ## ... Similar to previous best
    ## Run 469 stress 0.06231324 
    ## Run 470 stress 0.05051308 
    ## ... Procrustes: rmse 0.03570841  max resid 0.1219952 
    ## Run 471 stress 0.05403581 
    ## Run 472 stress 0.06136678 
    ## Run 473 stress 0.06231302 
    ## Run 474 stress 0.05051321 
    ## ... Procrustes: rmse 0.03572453  max resid 0.1220499 
    ## Run 475 stress 0.3657893 
    ## Run 476 stress 0.06231305 
    ## Run 477 stress 0.05036796 
    ## ... Procrustes: rmse 9.197023e-05  max resid 0.0002103698 
    ## ... Similar to previous best
    ## Run 478 stress 0.05051316 
    ## ... Procrustes: rmse 0.03569469  max resid 0.1219585 
    ## Run 479 stress 0.05036808 
    ## ... Procrustes: rmse 0.0002038571  max resid 0.000488914 
    ## ... Similar to previous best
    ## Run 480 stress 0.06309845 
    ## Run 481 stress 0.06362783 
    ## Run 482 stress 0.05051333 
    ## ... Procrustes: rmse 0.03570802  max resid 0.1220041 
    ## Run 483 stress 0.05051327 
    ## ... Procrustes: rmse 0.03565781  max resid 0.1216939 
    ## Run 484 stress 0.05928228 
    ## Run 485 stress 0.0597923 
    ## Run 486 stress 0.06231309 
    ## Run 487 stress 0.05051322 
    ## ... Procrustes: rmse 0.03565437  max resid 0.1216877 
    ## Run 488 stress 0.05036807 
    ## ... Procrustes: rmse 0.0001769504  max resid 0.0003913299 
    ## ... Similar to previous best
    ## Run 489 stress 0.05051339 
    ## ... Procrustes: rmse 0.03572333  max resid 0.1218816 
    ## Run 490 stress 0.05171058 
    ## Run 491 stress 0.05051344 
    ## ... Procrustes: rmse 0.03574077  max resid 0.1221062 
    ## Run 492 stress 0.05036795 
    ## ... Procrustes: rmse 9.748704e-05  max resid 0.0002184061 
    ## ... Similar to previous best
    ## Run 493 stress 0.0517106 
    ## Run 494 stress 0.05337593 
    ## Run 495 stress 0.05466701 
    ## Run 496 stress 0.06504224 
    ## Run 497 stress 0.0505131 
    ## ... Procrustes: rmse 0.03567795  max resid 0.1217702 
    ## Run 498 stress 0.06478535 
    ## Run 499 stress 0.05036791 
    ## ... Procrustes: rmse 1.503809e-05  max resid 2.890716e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.0613668 
    ## *** Best solution repeated 84 times

``` r
# Ocean sites and mixed lakes
PD_beta_OM_NMDS <- metaMDS(PD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1016168 
    ## Run 1 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0005628893  max resid 0.001573775 
    ## ... Similar to previous best
    ## Run 2 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005288972  max resid 0.001475712 
    ## ... Similar to previous best
    ## Run 3 stress 0.135615 
    ## Run 4 stress 0.1356145 
    ## Run 5 stress 0.101617 
    ## ... Procrustes: rmse 0.0006325507  max resid 0.001770142 
    ## ... Similar to previous best
    ## Run 6 stress 0.1016164 
    ## ... Procrustes: rmse 3.505939e-05  max resid 9.8139e-05 
    ## ... Similar to previous best
    ## Run 7 stress 0.1016164 
    ## ... Procrustes: rmse 1.414625e-05  max resid 3.911194e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.1341868 
    ## Run 9 stress 0.1016164 
    ## ... Procrustes: rmse 1.533116e-05  max resid 4.281264e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.1356147 
    ## Run 11 stress 0.1016164 
    ## ... Procrustes: rmse 4.172077e-05  max resid 0.0001172944 
    ## ... Similar to previous best
    ## Run 12 stress 0.1341868 
    ## Run 13 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002876666  max resid 0.0008064082 
    ## ... Similar to previous best
    ## Run 14 stress 0.1016164 
    ## ... Procrustes: rmse 2.602708e-05  max resid 7.011612e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.1537554 
    ## Run 16 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 7.590634e-06  max resid 1.739196e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.1016164 
    ## ... Procrustes: rmse 1.754824e-05  max resid 4.951363e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005362768  max resid 0.001499699 
    ## ... Similar to previous best
    ## Run 19 stress 0.1519033 
    ## Run 20 stress 0.1341868 
    ## Run 21 stress 0.1016164 
    ## ... Procrustes: rmse 2.782919e-05  max resid 7.835611e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.1016164 
    ## ... Procrustes: rmse 1.252422e-05  max resid 3.531558e-05 
    ## ... Similar to previous best
    ## Run 23 stress 0.1016164 
    ## ... Procrustes: rmse 6.620379e-05  max resid 0.0001865664 
    ## ... Similar to previous best
    ## Run 24 stress 0.1016168 
    ## ... Procrustes: rmse 0.00050004  max resid 0.001398572 
    ## ... Similar to previous best
    ## Run 25 stress 0.1016164 
    ## ... Procrustes: rmse 3.859142e-05  max resid 0.0001088935 
    ## ... Similar to previous best
    ## Run 26 stress 0.1356148 
    ## Run 27 stress 0.1356145 
    ## Run 28 stress 0.3106949 
    ## Run 29 stress 0.135615 
    ## Run 30 stress 0.1356152 
    ## Run 31 stress 0.1016164 
    ## ... Procrustes: rmse 2.112259e-05  max resid 5.67157e-05 
    ## ... Similar to previous best
    ## Run 32 stress 0.1016164 
    ## ... Procrustes: rmse 3.279674e-05  max resid 9.21301e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.1341868 
    ## Run 34 stress 0.1016164 
    ## ... Procrustes: rmse 5.872982e-05  max resid 0.0001652995 
    ## ... Similar to previous best
    ## Run 35 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006107349  max resid 0.00170559 
    ## ... Similar to previous best
    ## Run 36 stress 0.1016164 
    ## ... Procrustes: rmse 5.721133e-05  max resid 0.0001592287 
    ## ... Similar to previous best
    ## Run 37 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005349381  max resid 0.001495779 
    ## ... Similar to previous best
    ## Run 38 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002351507  max resid 0.0006587198 
    ## ... Similar to previous best
    ## Run 39 stress 0.135615 
    ## Run 40 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003283413  max resid 0.0009172759 
    ## ... Similar to previous best
    ## Run 41 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004375244  max resid 0.001222922 
    ## ... Similar to previous best
    ## Run 42 stress 0.1016164 
    ## ... Procrustes: rmse 8.087908e-05  max resid 0.0002258669 
    ## ... Similar to previous best
    ## Run 43 stress 0.1016169 
    ## ... Procrustes: rmse 0.000611851  max resid 0.001710132 
    ## ... Similar to previous best
    ## Run 44 stress 0.1016164 
    ## ... Procrustes: rmse 7.930975e-06  max resid 2.113378e-05 
    ## ... Similar to previous best
    ## Run 45 stress 0.1016164 
    ## ... Procrustes: rmse 2.119516e-05  max resid 5.865481e-05 
    ## ... Similar to previous best
    ## Run 46 stress 0.1016164 
    ## ... Procrustes: rmse 7.746403e-06  max resid 2.166146e-05 
    ## ... Similar to previous best
    ## Run 47 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003635321  max resid 0.001017299 
    ## ... Similar to previous best
    ## Run 48 stress 0.1341868 
    ## Run 49 stress 0.1341868 
    ## Run 50 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001333309  max resid 0.0003751624 
    ## ... Similar to previous best
    ## Run 51 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006228068  max resid 0.001738596 
    ## ... Similar to previous best
    ## Run 52 stress 0.1016164 
    ## ... Procrustes: rmse 1.667626e-05  max resid 4.510059e-05 
    ## ... Similar to previous best
    ## Run 53 stress 0.1016164 
    ## ... Procrustes: rmse 3.715783e-05  max resid 0.0001042618 
    ## ... Similar to previous best
    ## Run 54 stress 0.1356148 
    ## Run 55 stress 0.1341868 
    ## Run 56 stress 0.1341868 
    ## Run 57 stress 0.1341868 
    ## Run 58 stress 0.1016164 
    ## ... Procrustes: rmse 9.290306e-05  max resid 0.0002615425 
    ## ... Similar to previous best
    ## Run 59 stress 0.1016164 
    ## ... Procrustes: rmse 1.846068e-05  max resid 5.201386e-05 
    ## ... Similar to previous best
    ## Run 60 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003240596  max resid 0.0009083249 
    ## ... Similar to previous best
    ## Run 61 stress 0.1341868 
    ## Run 62 stress 0.1356145 
    ## Run 63 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003601109  max resid 0.001009375 
    ## ... Similar to previous best
    ## Run 64 stress 0.1519033 
    ## Run 65 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003974644  max resid 0.001112923 
    ## ... Similar to previous best
    ## Run 66 stress 0.135615 
    ## Run 67 stress 0.1341868 
    ## Run 68 stress 0.1341868 
    ## Run 69 stress 0.1356146 
    ## Run 70 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003432189  max resid 0.0009576384 
    ## ... Similar to previous best
    ## Run 71 stress 0.135615 
    ## Run 72 stress 0.1341868 
    ## Run 73 stress 0.1016164 
    ## ... Procrustes: rmse 2.349449e-05  max resid 3.904915e-05 
    ## ... Similar to previous best
    ## Run 74 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006182842  max resid 0.001727175 
    ## ... Similar to previous best
    ## Run 75 stress 0.3495731 
    ## Run 76 stress 0.1016164 
    ## ... Procrustes: rmse 8.343282e-06  max resid 1.69304e-05 
    ## ... Similar to previous best
    ## Run 77 stress 0.1341868 
    ## Run 78 stress 0.1016164 
    ## ... Procrustes: rmse 1.31846e-05  max resid 3.73076e-05 
    ## ... Similar to previous best
    ## Run 79 stress 0.1016165 
    ## ... Procrustes: rmse 0.000284461  max resid 0.000797831 
    ## ... Similar to previous best
    ## Run 80 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004982233  max resid 0.001393445 
    ## ... Similar to previous best
    ## Run 81 stress 0.1530424 
    ## Run 82 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002085023  max resid 0.000585085 
    ## ... Similar to previous best
    ## Run 83 stress 0.1016164 
    ## ... Procrustes: rmse 1.733628e-05  max resid 4.879143e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.1016164 
    ## ... Procrustes: rmse 6.673976e-05  max resid 0.000188112 
    ## ... Similar to previous best
    ## Run 85 stress 0.1341868 
    ## Run 86 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002208157  max resid 0.0006163135 
    ## ... Similar to previous best
    ## Run 87 stress 0.1356148 
    ## Run 88 stress 0.3065089 
    ## Run 89 stress 0.1016164 
    ## ... Procrustes: rmse 2.060898e-06  max resid 5.301574e-06 
    ## ... Similar to previous best
    ## Run 90 stress 0.1341868 
    ## Run 91 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004515706  max resid 0.001263809 
    ## ... Similar to previous best
    ## Run 92 stress 0.1016164 
    ## ... Procrustes: rmse 6.591612e-05  max resid 0.0001835384 
    ## ... Similar to previous best
    ## Run 93 stress 0.3186006 
    ## Run 94 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005484  max resid 0.001529588 
    ## ... Similar to previous best
    ## Run 95 stress 0.1341868 
    ## Run 96 stress 0.1016164 
    ## ... Procrustes: rmse 1.712688e-05  max resid 4.824701e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.1016164 
    ## ... Procrustes: rmse 7.366856e-05  max resid 0.000207759 
    ## ... Similar to previous best
    ## Run 98 stress 0.1016164 
    ## ... Procrustes: rmse 4.229422e-05  max resid 0.0001193254 
    ## ... Similar to previous best
    ## Run 99 stress 0.1016164 
    ## ... Procrustes: rmse 2.829126e-05  max resid 7.938486e-05 
    ## ... Similar to previous best
    ## Run 100 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004076773  max resid 0.001141583 
    ## ... Similar to previous best
    ## Run 101 stress 0.1016164 
    ## ... Procrustes: rmse 2.551769e-05  max resid 7.150252e-05 
    ## ... Similar to previous best
    ## Run 102 stress 0.1016164 
    ## ... Procrustes: rmse 4.988834e-05  max resid 0.0001407324 
    ## ... Similar to previous best
    ## Run 103 stress 0.1016165 
    ## ... Procrustes: rmse 0.000265475  max resid 0.0007446074 
    ## ... Similar to previous best
    ## Run 104 stress 0.1016164 
    ## ... Procrustes: rmse 1.859224e-05  max resid 5.257566e-05 
    ## ... Similar to previous best
    ## Run 105 stress 0.1530426 
    ## Run 106 stress 0.1341868 
    ## Run 107 stress 0.1016169 
    ## ... Procrustes: rmse 0.000608294  max resid 0.001696209 
    ## ... Similar to previous best
    ## Run 108 stress 0.1356151 
    ## Run 109 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002656662  max resid 0.0007454947 
    ## ... Similar to previous best
    ## Run 110 stress 0.1341868 
    ## Run 111 stress 0.1356147 
    ## Run 112 stress 0.1016164 
    ## ... Procrustes: rmse 1.088774e-05  max resid 2.403483e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.1016164 
    ## ... Procrustes: rmse 8.084662e-05  max resid 0.0002266658 
    ## ... Similar to previous best
    ## Run 114 stress 0.1016164 
    ## ... Procrustes: rmse 4.231847e-05  max resid 0.0001192601 
    ## ... Similar to previous best
    ## Run 115 stress 0.1016164 
    ## ... Procrustes: rmse 2.100876e-05  max resid 4.298123e-05 
    ## ... Similar to previous best
    ## Run 116 stress 0.1356149 
    ## Run 117 stress 0.1016164 
    ## ... Procrustes: rmse 6.03706e-05  max resid 0.0001698307 
    ## ... Similar to previous best
    ## Run 118 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004833032  max resid 0.001352017 
    ## ... Similar to previous best
    ## Run 119 stress 0.1016164 
    ## ... Procrustes: rmse 1.149562e-05  max resid 3.242211e-05 
    ## ... Similar to previous best
    ## Run 120 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006038553  max resid 0.001688391 
    ## ... Similar to previous best
    ## Run 121 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003601657  max resid 0.001009524 
    ## ... Similar to previous best
    ## Run 122 stress 0.1356145 
    ## Run 123 stress 0.1519033 
    ## Run 124 stress 0.1016165 
    ## ... Procrustes: rmse 0.000182575  max resid 0.0005122715 
    ## ... Similar to previous best
    ## Run 125 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003162399  max resid 0.0008857949 
    ## ... Similar to previous best
    ## Run 126 stress 0.1341868 
    ## Run 127 stress 0.1341868 
    ## Run 128 stress 0.1016164 
    ## ... Procrustes: rmse 6.526737e-06  max resid 1.80456e-05 
    ## ... Similar to previous best
    ## Run 129 stress 0.1341868 
    ## Run 130 stress 0.101617 
    ## ... Procrustes: rmse 0.0007056463  max resid 0.001972503 
    ## ... Similar to previous best
    ## Run 131 stress 0.1016164 
    ## ... Procrustes: rmse 3.009837e-05  max resid 8.491209e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.1356147 
    ## Run 133 stress 0.1530425 
    ## Run 134 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003005588  max resid 0.0008411236 
    ## ... Similar to previous best
    ## Run 135 stress 0.1016164 
    ## ... Procrustes: rmse 1.702052e-05  max resid 4.788686e-05 
    ## ... Similar to previous best
    ## Run 136 stress 0.1356145 
    ## Run 137 stress 0.1530423 
    ## Run 138 stress 0.1016164 
    ## ... Procrustes: rmse 3.196794e-06  max resid 5.764484e-06 
    ## ... Similar to previous best
    ## Run 139 stress 0.1016164 
    ## ... Procrustes: rmse 4.40426e-05  max resid 0.0001241263 
    ## ... Similar to previous best
    ## Run 140 stress 0.1356149 
    ## Run 141 stress 0.1341868 
    ## Run 142 stress 0.1341868 
    ## Run 143 stress 0.1341868 
    ## Run 144 stress 0.1016164 
    ## ... Procrustes: rmse 5.200295e-06  max resid 1.462761e-05 
    ## ... Similar to previous best
    ## Run 145 stress 0.1016164 
    ## ... Procrustes: rmse 8.762347e-05  max resid 0.0002469341 
    ## ... Similar to previous best
    ## Run 146 stress 0.1356151 
    ## Run 147 stress 0.1016164 
    ## ... Procrustes: rmse 2.475073e-05  max resid 6.968665e-05 
    ## ... Similar to previous best
    ## Run 148 stress 0.1016164 
    ## ... Procrustes: rmse 2.31914e-05  max resid 6.496788e-05 
    ## ... Similar to previous best
    ## Run 149 stress 0.1341868 
    ## Run 150 stress 0.1016164 
    ## ... Procrustes: rmse 2.32655e-05  max resid 6.557495e-05 
    ## ... Similar to previous best
    ## Run 151 stress 0.1016164 
    ## ... Procrustes: rmse 1.626503e-05  max resid 4.502573e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.1016164 
    ## ... Procrustes: rmse 4.415024e-05  max resid 0.0001243004 
    ## ... Similar to previous best
    ## Run 153 stress 0.3316359 
    ## Run 154 stress 0.1016164 
    ## ... Procrustes: rmse 6.317125e-05  max resid 0.0001784703 
    ## ... Similar to previous best
    ## Run 155 stress 0.1341868 
    ## Run 156 stress 0.1016164 
    ## ... Procrustes: rmse 2.516961e-05  max resid 7.059353e-05 
    ## ... Similar to previous best
    ## Run 157 stress 0.1341868 
    ## Run 158 stress 0.1016164 
    ## ... Procrustes: rmse 6.143475e-05  max resid 0.0001732573 
    ## ... Similar to previous best
    ## Run 159 stress 0.1341868 
    ## Run 160 stress 0.1016164 
    ## ... Procrustes: rmse 6.337342e-05  max resid 0.0001788333 
    ## ... Similar to previous best
    ## Run 161 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005487059  max resid 0.001530821 
    ## ... Similar to previous best
    ## Run 162 stress 0.3578489 
    ## Run 163 stress 0.1016164 
    ## ... Procrustes: rmse 1.455042e-06  max resid 3.606423e-06 
    ## ... Similar to previous best
    ## Run 164 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006199603  max resid 0.001732873 
    ## ... Similar to previous best
    ## Run 165 stress 0.1016167 
    ## ... Procrustes: rmse 0.000425051  max resid 0.001191233 
    ## ... Similar to previous best
    ## Run 166 stress 0.1016164 
    ## ... Procrustes: rmse 2.506104e-05  max resid 7.043818e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.1016164 
    ## ... Procrustes: rmse 4.068992e-05  max resid 0.0001146641 
    ## ... Similar to previous best
    ## Run 168 stress 0.1016164 
    ## ... Procrustes: rmse 2.85271e-05  max resid 7.007238e-05 
    ## ... Similar to previous best
    ## Run 169 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005649089  max resid 0.00157915 
    ## ... Similar to previous best
    ## Run 170 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001295472  max resid 0.0003645995 
    ## ... Similar to previous best
    ## Run 171 stress 0.1341868 
    ## Run 172 stress 0.1016164 
    ## ... Procrustes: rmse 5.978067e-06  max resid 1.561292e-05 
    ## ... Similar to previous best
    ## Run 173 stress 0.1356147 
    ## Run 174 stress 0.1016164 
    ## ... Procrustes: rmse 1.841205e-05  max resid 5.192349e-05 
    ## ... Similar to previous best
    ## Run 175 stress 0.1016166 
    ## ... Procrustes: rmse 0.000308866  max resid 0.0008651329 
    ## ... Similar to previous best
    ## Run 176 stress 0.1016164 
    ## ... Procrustes: rmse 4.475782e-05  max resid 0.0001218534 
    ## ... Similar to previous best
    ## Run 177 stress 0.1016164 
    ## ... Procrustes: rmse 5.055414e-05  max resid 0.0001421257 
    ## ... Similar to previous best
    ## Run 178 stress 0.1016164 
    ## ... Procrustes: rmse 5.598304e-05  max resid 0.0001573301 
    ## ... Similar to previous best
    ## Run 179 stress 0.1341868 
    ## Run 180 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 9.315804e-06  max resid 2.581977e-05 
    ## ... Similar to previous best
    ## Run 181 stress 0.1016166 
    ## ... Procrustes: rmse 0.00037747  max resid 0.001056142 
    ## ... Similar to previous best
    ## Run 182 stress 0.1341868 
    ## Run 183 stress 0.1016164 
    ## ... Procrustes: rmse 1.795308e-05  max resid 5.156238e-05 
    ## ... Similar to previous best
    ## Run 184 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004228236  max resid 0.00118179 
    ## ... Similar to previous best
    ## Run 185 stress 0.1356145 
    ## Run 186 stress 0.1016164 
    ## ... Procrustes: rmse 2.64604e-05  max resid 7.463258e-05 
    ## ... Similar to previous best
    ## Run 187 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005763601  max resid 0.001610571 
    ## ... Similar to previous best
    ## Run 188 stress 0.1016164 
    ## ... Procrustes: rmse 7.435569e-06  max resid 1.966278e-05 
    ## ... Similar to previous best
    ## Run 189 stress 0.1016164 
    ## ... Procrustes: rmse 2.266407e-05  max resid 6.349735e-05 
    ## ... Similar to previous best
    ## Run 190 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003591325  max resid 0.001005126 
    ## ... Similar to previous best
    ## Run 191 stress 0.3411883 
    ## Run 192 stress 0.1356148 
    ## Run 193 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004423426  max resid 0.001237898 
    ## ... Similar to previous best
    ## Run 194 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003201043  max resid 0.00089672 
    ## ... Similar to previous best
    ## Run 195 stress 0.1016164 
    ## ... Procrustes: rmse 1.089148e-05  max resid 2.67688e-05 
    ## ... Similar to previous best
    ## Run 196 stress 0.1356143 
    ## Run 197 stress 0.1341868 
    ## Run 198 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003965077  max resid 0.001109926 
    ## ... Similar to previous best
    ## Run 199 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005700734  max resid 0.001592911 
    ## ... Similar to previous best
    ## Run 200 stress 0.1016164 
    ## ... Procrustes: rmse 9.713505e-05  max resid 0.0002727943 
    ## ... Similar to previous best
    ## Run 201 stress 0.1356147 
    ## Run 202 stress 0.1016164 
    ## ... Procrustes: rmse 2.604529e-05  max resid 7.316791e-05 
    ## ... Similar to previous best
    ## Run 203 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001059002  max resid 0.0002979419 
    ## ... Similar to previous best
    ## Run 204 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001030652  max resid 0.0002901933 
    ## ... Similar to previous best
    ## Run 205 stress 0.1341868 
    ## Run 206 stress 0.1341868 
    ## Run 207 stress 0.1016164 
    ## ... Procrustes: rmse 4.831464e-06  max resid 1.286647e-05 
    ## ... Similar to previous best
    ## Run 208 stress 0.1341868 
    ## Run 209 stress 0.1530423 
    ## Run 210 stress 0.1356146 
    ## Run 211 stress 0.1016164 
    ## ... Procrustes: rmse 5.893043e-06  max resid 1.679762e-05 
    ## ... Similar to previous best
    ## Run 212 stress 0.1016164 
    ## ... Procrustes: rmse 8.674266e-05  max resid 0.0002439477 
    ## ... Similar to previous best
    ## Run 213 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005110253  max resid 0.001429723 
    ## ... Similar to previous best
    ## Run 214 stress 0.1016164 
    ## ... Procrustes: rmse 6.982118e-05  max resid 0.0001962761 
    ## ... Similar to previous best
    ## Run 215 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003757968  max resid 0.001051342 
    ## ... Similar to previous best
    ## Run 216 stress 0.101617 
    ## ... Procrustes: rmse 0.0006905514  max resid 0.001930162 
    ## ... Similar to previous best
    ## Run 217 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005853294  max resid 0.001635911 
    ## ... Similar to previous best
    ## Run 218 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002822136  max resid 0.0007910597 
    ## ... Similar to previous best
    ## Run 219 stress 0.1016164 
    ## ... Procrustes: rmse 1.499511e-05  max resid 4.144534e-05 
    ## ... Similar to previous best
    ## Run 220 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004332631  max resid 0.001212094 
    ## ... Similar to previous best
    ## Run 221 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 1.6475e-06  max resid 3.513878e-06 
    ## ... Similar to previous best
    ## Run 222 stress 0.1016164 
    ## ... Procrustes: rmse 2.246258e-05  max resid 5.798058e-05 
    ## ... Similar to previous best
    ## Run 223 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001052723  max resid 0.0002971279 
    ## ... Similar to previous best
    ## Run 224 stress 0.1016166 
    ## ... Procrustes: rmse 0.000362334  max resid 0.001013712 
    ## ... Similar to previous best
    ## Run 225 stress 0.1341868 
    ## Run 226 stress 0.1341868 
    ## Run 227 stress 0.1016164 
    ## ... Procrustes: rmse 8.412382e-05  max resid 0.000237686 
    ## ... Similar to previous best
    ## Run 228 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001210331  max resid 0.000339868 
    ## ... Similar to previous best
    ## Run 229 stress 0.1537554 
    ## Run 230 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002669509  max resid 0.0007480545 
    ## ... Similar to previous best
    ## Run 231 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006141033  max resid 0.001713787 
    ## ... Similar to previous best
    ## Run 232 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004610942  max resid 0.00128798 
    ## ... Similar to previous best
    ## Run 233 stress 0.1356143 
    ## Run 234 stress 0.1016164 
    ## ... Procrustes: rmse 1.902573e-05  max resid 5.366542e-05 
    ## ... Similar to previous best
    ## Run 235 stress 0.1356147 
    ## Run 236 stress 0.1341868 
    ## Run 237 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001027765  max resid 0.0002894725 
    ## ... Similar to previous best
    ## Run 238 stress 0.1016164 
    ## ... Procrustes: rmse 1.806529e-06  max resid 4.293336e-06 
    ## ... Similar to previous best
    ## Run 239 stress 0.1530429 
    ## Run 240 stress 0.1341868 
    ## Run 241 stress 0.1519033 
    ## Run 242 stress 0.1356145 
    ## Run 243 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003393145  max resid 0.000939369 
    ## ... Similar to previous best
    ## Run 244 stress 0.1016164 
    ## ... Procrustes: rmse 5.781227e-06  max resid 1.104674e-05 
    ## ... Similar to previous best
    ## Run 245 stress 0.1016164 
    ## ... Procrustes: rmse 0.000118822  max resid 0.0003334423 
    ## ... Similar to previous best
    ## Run 246 stress 0.1016164 
    ## ... Procrustes: rmse 7.36946e-07  max resid 1.39605e-06 
    ## ... Similar to previous best
    ## Run 247 stress 0.1341868 
    ## Run 248 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003279188  max resid 0.0009166518 
    ## ... Similar to previous best
    ## Run 249 stress 0.1016164 
    ## ... Procrustes: rmse 8.740343e-06  max resid 2.27606e-05 
    ## ... Similar to previous best
    ## Run 250 stress 0.1537554 
    ## Run 251 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002597862  max resid 0.0007277757 
    ## ... Similar to previous best
    ## Run 252 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003362443  max resid 0.0009415773 
    ## ... Similar to previous best
    ## Run 253 stress 0.1341868 
    ## Run 254 stress 0.1016164 
    ## ... Procrustes: rmse 1.24747e-05  max resid 3.514405e-05 
    ## ... Similar to previous best
    ## Run 255 stress 0.1016164 
    ## ... Procrustes: rmse 6.063748e-05  max resid 0.0001709681 
    ## ... Similar to previous best
    ## Run 256 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004878842  max resid 0.001363327 
    ## ... Similar to previous best
    ## Run 257 stress 0.1016164 
    ## ... Procrustes: rmse 5.512125e-05  max resid 0.0001544669 
    ## ... Similar to previous best
    ## Run 258 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002281972  max resid 0.0006398518 
    ## ... Similar to previous best
    ## Run 259 stress 0.1341868 
    ## Run 260 stress 0.1016164 
    ## ... Procrustes: rmse 1.355164e-05  max resid 3.820531e-05 
    ## ... Similar to previous best
    ## Run 261 stress 0.1341868 
    ## Run 262 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005520039  max resid 0.001542424 
    ## ... Similar to previous best
    ## Run 263 stress 0.1016164 
    ## ... Procrustes: rmse 3.460585e-05  max resid 9.788013e-05 
    ## ... Similar to previous best
    ## Run 264 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001293094  max resid 0.00036365 
    ## ... Similar to previous best
    ## Run 265 stress 0.1016164 
    ## ... Procrustes: rmse 6.033107e-05  max resid 0.0001702066 
    ## ... Similar to previous best
    ## Run 266 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002182119  max resid 0.0006094294 
    ## ... Similar to previous best
    ## Run 267 stress 0.1016164 
    ## ... Procrustes: rmse 1.717082e-06  max resid 3.60591e-06 
    ## ... Similar to previous best
    ## Run 268 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005281984  max resid 0.001475907 
    ## ... Similar to previous best
    ## Run 269 stress 0.1016164 
    ## ... Procrustes: rmse 3.575173e-05  max resid 0.0001007478 
    ## ... Similar to previous best
    ## Run 270 stress 0.1519033 
    ## Run 271 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005067286  max resid 0.001413521 
    ## ... Similar to previous best
    ## Run 272 stress 0.1341868 
    ## Run 273 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004708924  max resid 0.001316958 
    ## ... Similar to previous best
    ## Run 274 stress 0.1356152 
    ## Run 275 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004586862  max resid 0.001282429 
    ## ... Similar to previous best
    ## Run 276 stress 0.1016164 
    ## ... Procrustes: rmse 3.34915e-06  max resid 8.848663e-06 
    ## ... Similar to previous best
    ## Run 277 stress 0.1016164 
    ## ... Procrustes: rmse 2.820257e-05  max resid 7.946084e-05 
    ## ... Similar to previous best
    ## Run 278 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002109138  max resid 0.0005906447 
    ## ... Similar to previous best
    ## Run 279 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001665684  max resid 0.0004659566 
    ## ... Similar to previous best
    ## Run 280 stress 0.1530423 
    ## Run 281 stress 0.1341868 
    ## Run 282 stress 0.1356145 
    ## Run 283 stress 0.1356144 
    ## Run 284 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004772093  max resid 0.001333454 
    ## ... Similar to previous best
    ## Run 285 stress 0.1356149 
    ## Run 286 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004705866  max resid 0.001315018 
    ## ... Similar to previous best
    ## Run 287 stress 0.1016164 
    ## ... Procrustes: rmse 8.982081e-06  max resid 2.166512e-05 
    ## ... Similar to previous best
    ## Run 288 stress 0.3578485 
    ## Run 289 stress 0.1016164 
    ## ... Procrustes: rmse 6.320633e-05  max resid 0.000177892 
    ## ... Similar to previous best
    ## Run 290 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007304067  max resid 0.002040605 
    ## ... Similar to previous best
    ## Run 291 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003984424  max resid 0.001114552 
    ## ... Similar to previous best
    ## Run 292 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001104599  max resid 0.000312119 
    ## ... Similar to previous best
    ## Run 293 stress 0.1530428 
    ## Run 294 stress 0.1519033 
    ## Run 295 stress 0.1016164 
    ## ... Procrustes: rmse 7.991682e-05  max resid 0.0002188314 
    ## ... Similar to previous best
    ## Run 296 stress 0.1016164 
    ## ... Procrustes: rmse 6.320126e-05  max resid 0.0001783208 
    ## ... Similar to previous best
    ## Run 297 stress 0.1016164 
    ## ... Procrustes: rmse 1.193207e-05  max resid 3.182577e-05 
    ## ... Similar to previous best
    ## Run 298 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005355751  max resid 0.001496954 
    ## ... Similar to previous best
    ## Run 299 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005447324  max resid 0.001522993 
    ## ... Similar to previous best
    ## Run 300 stress 0.1016164 
    ## ... Procrustes: rmse 4.515588e-05  max resid 0.0001274807 
    ## ... Similar to previous best
    ## Run 301 stress 0.1519033 
    ## Run 302 stress 0.1016164 
    ## ... Procrustes: rmse 9.028946e-06  max resid 1.714915e-05 
    ## ... Similar to previous best
    ## Run 303 stress 0.1016164 
    ## ... Procrustes: rmse 8.009319e-05  max resid 0.0002257126 
    ## ... Similar to previous best
    ## Run 304 stress 0.1356148 
    ## Run 305 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005542212  max resid 0.001547286 
    ## ... Similar to previous best
    ## Run 306 stress 0.1016164 
    ## ... Procrustes: rmse 9.155151e-06  max resid 2.582963e-05 
    ## ... Similar to previous best
    ## Run 307 stress 0.101617 
    ## ... Procrustes: rmse 0.0006626732  max resid 0.001851986 
    ## ... Similar to previous best
    ## Run 308 stress 0.1341868 
    ## Run 309 stress 0.1341868 
    ## Run 310 stress 0.1341868 
    ## Run 311 stress 0.1341868 
    ## Run 312 stress 0.1016164 
    ## ... Procrustes: rmse 8.847994e-05  max resid 0.0002498577 
    ## ... Similar to previous best
    ## Run 313 stress 0.1356148 
    ## Run 314 stress 0.212628 
    ## Run 315 stress 0.1016164 
    ## ... Procrustes: rmse 7.026732e-05  max resid 0.0001981596 
    ## ... Similar to previous best
    ## Run 316 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005499911  max resid 0.001535257 
    ## ... Similar to previous best
    ## Run 317 stress 0.1356149 
    ## Run 318 stress 0.1016164 
    ## ... Procrustes: rmse 5.767373e-05  max resid 0.0001617394 
    ## ... Similar to previous best
    ## Run 319 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006350071  max resid 0.001773205 
    ## ... Similar to previous best
    ## Run 320 stress 0.1016164 
    ## ... Procrustes: rmse 1.274037e-05  max resid 3.601607e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.1016164 
    ## ... Procrustes: rmse 1.012959e-05  max resid 2.855851e-05 
    ## ... Similar to previous best
    ## Run 322 stress 0.1341868 
    ## Run 323 stress 0.3411844 
    ## Run 324 stress 0.1341868 
    ## Run 325 stress 0.1016164 
    ## ... Procrustes: rmse 1.473366e-05  max resid 4.161135e-05 
    ## ... Similar to previous best
    ## Run 326 stress 0.1016164 
    ## ... Procrustes: rmse 4.136574e-05  max resid 0.0001159702 
    ## ... Similar to previous best
    ## Run 327 stress 0.1016164 
    ## ... Procrustes: rmse 1.195305e-05  max resid 3.378846e-05 
    ## ... Similar to previous best
    ## Run 328 stress 0.1016164 
    ## ... Procrustes: rmse 2.116254e-05  max resid 3.948604e-05 
    ## ... Similar to previous best
    ## Run 329 stress 0.1016164 
    ## ... Procrustes: rmse 6.715929e-06  max resid 1.893981e-05 
    ## ... Similar to previous best
    ## Run 330 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001128171  max resid 0.0003155752 
    ## ... Similar to previous best
    ## Run 331 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003466757  max resid 0.0009703784 
    ## ... Similar to previous best
    ## Run 332 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001773508  max resid 0.000497639 
    ## ... Similar to previous best
    ## Run 333 stress 0.1016164 
    ## ... Procrustes: rmse 9.183736e-05  max resid 0.0002589316 
    ## ... Similar to previous best
    ## Run 334 stress 0.1016164 
    ## ... Procrustes: rmse 7.836182e-05  max resid 0.0002205005 
    ## ... Similar to previous best
    ## Run 335 stress 0.1356143 
    ## Run 336 stress 0.1530427 
    ## Run 337 stress 0.1016171 
    ## ... Procrustes: rmse 0.000718259  max resid 0.002006815 
    ## ... Similar to previous best
    ## Run 338 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005479081  max resid 0.001531742 
    ## ... Similar to previous best
    ## Run 339 stress 0.1341868 
    ## Run 340 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003032092  max resid 0.0008491104 
    ## ... Similar to previous best
    ## Run 341 stress 0.1016164 
    ## ... Procrustes: rmse 2.986018e-05  max resid 8.334855e-05 
    ## ... Similar to previous best
    ## Run 342 stress 0.1016167 
    ## ... Procrustes: rmse 0.000438317  max resid 0.001223866 
    ## ... Similar to previous best
    ## Run 343 stress 0.1016164 
    ## ... Procrustes: rmse 5.074058e-05  max resid 0.0001421448 
    ## ... Similar to previous best
    ## Run 344 stress 0.1016164 
    ## ... Procrustes: rmse 5.906331e-05  max resid 0.0001664425 
    ## ... Similar to previous best
    ## Run 345 stress 0.1016164 
    ## ... Procrustes: rmse 3.172647e-05  max resid 8.911268e-05 
    ## ... Similar to previous best
    ## Run 346 stress 0.1356147 
    ## Run 347 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002446579  max resid 0.0006856521 
    ## ... Similar to previous best
    ## Run 348 stress 0.1341868 
    ## Run 349 stress 0.1341868 
    ## Run 350 stress 0.101617 
    ## ... Procrustes: rmse 0.0006691971  max resid 0.001870242 
    ## ... Similar to previous best
    ## Run 351 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002289828  max resid 0.0006429241 
    ## ... Similar to previous best
    ## Run 352 stress 0.1016164 
    ## ... Procrustes: rmse 1.322142e-05  max resid 3.726206e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.1519033 
    ## Run 354 stress 0.1016164 
    ## ... Procrustes: rmse 8.807728e-06  max resid 2.489373e-05 
    ## ... Similar to previous best
    ## Run 355 stress 0.1016164 
    ## ... Procrustes: rmse 1.20609e-05  max resid 3.401346e-05 
    ## ... Similar to previous best
    ## Run 356 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004387225  max resid 0.001227445 
    ## ... Similar to previous best
    ## Run 357 stress 0.1519033 
    ## Run 358 stress 0.1356146 
    ## Run 359 stress 0.1341868 
    ## Run 360 stress 0.1016164 
    ## ... Procrustes: rmse 5.050826e-05  max resid 0.0001332622 
    ## ... Similar to previous best
    ## Run 361 stress 0.1016164 
    ## ... Procrustes: rmse 3.228513e-05  max resid 9.086234e-05 
    ## ... Similar to previous best
    ## Run 362 stress 0.1016164 
    ## ... Procrustes: rmse 7.627752e-06  max resid 1.931603e-05 
    ## ... Similar to previous best
    ## Run 363 stress 0.1016164 
    ## ... Procrustes: rmse 9.306512e-05  max resid 0.0002627471 
    ## ... Similar to previous best
    ## Run 364 stress 0.1341868 
    ## Run 365 stress 0.1016164 
    ## ... Procrustes: rmse 7.511991e-05  max resid 0.0002113683 
    ## ... Similar to previous best
    ## Run 366 stress 0.1341868 
    ## Run 367 stress 0.1016164 
    ## ... Procrustes: rmse 1.200827e-05  max resid 3.37789e-05 
    ## ... Similar to previous best
    ## Run 368 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003979915  max resid 0.001114686 
    ## ... Similar to previous best
    ## Run 369 stress 0.1016164 
    ## ... Procrustes: rmse 1.759564e-05  max resid 2.937189e-05 
    ## ... Similar to previous best
    ## Run 370 stress 0.1537556 
    ## Run 371 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004246667  max resid 0.001187802 
    ## ... Similar to previous best
    ## Run 372 stress 0.1016164 
    ## ... Procrustes: rmse 3.354978e-05  max resid 9.480102e-05 
    ## ... Similar to previous best
    ## Run 373 stress 0.1016164 
    ## ... Procrustes: rmse 1.420726e-05  max resid 3.786058e-05 
    ## ... Similar to previous best
    ## Run 374 stress 0.1016164 
    ## ... Procrustes: rmse 4.966904e-05  max resid 0.0001398958 
    ## ... Similar to previous best
    ## Run 375 stress 0.1016164 
    ## ... Procrustes: rmse 6.226136e-05  max resid 0.0001752948 
    ## ... Similar to previous best
    ## Run 376 stress 0.1016164 
    ## ... Procrustes: rmse 4.293826e-05  max resid 0.0001211718 
    ## ... Similar to previous best
    ## Run 377 stress 0.1356145 
    ## Run 378 stress 0.1016164 
    ## ... Procrustes: rmse 4.094417e-05  max resid 0.000115073 
    ## ... Similar to previous best
    ## Run 379 stress 0.1016164 
    ## ... Procrustes: rmse 5.80025e-05  max resid 0.0001637591 
    ## ... Similar to previous best
    ## Run 380 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005902786  max resid 0.001647058 
    ## ... Similar to previous best
    ## Run 381 stress 0.1016164 
    ## ... Procrustes: rmse 5.508723e-05  max resid 0.0001550536 
    ## ... Similar to previous best
    ## Run 382 stress 0.1356151 
    ## Run 383 stress 0.1016164 
    ## ... Procrustes: rmse 2.066996e-05  max resid 5.810088e-05 
    ## ... Similar to previous best
    ## Run 384 stress 0.1341868 
    ## Run 385 stress 0.1016164 
    ## ... Procrustes: rmse 2.391787e-05  max resid 5.847048e-05 
    ## ... Similar to previous best
    ## Run 386 stress 0.1016164 
    ## ... Procrustes: rmse 5.87852e-06  max resid 1.538555e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.1341868 
    ## Run 388 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003333073  max resid 0.0009310357 
    ## ... Similar to previous best
    ## Run 389 stress 0.1016164 
    ## ... Procrustes: rmse 4.797292e-05  max resid 0.0001349037 
    ## ... Similar to previous best
    ## Run 390 stress 0.1341868 
    ## Run 391 stress 0.1519033 
    ## Run 392 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002651002  max resid 0.0007435329 
    ## ... Similar to previous best
    ## Run 393 stress 0.1530423 
    ## Run 394 stress 0.1016164 
    ## ... Procrustes: rmse 2.211646e-05  max resid 4.85809e-05 
    ## ... Similar to previous best
    ## Run 395 stress 0.1016164 
    ## ... Procrustes: rmse 1.906248e-05  max resid 3.724417e-05 
    ## ... Similar to previous best
    ## Run 396 stress 0.1341868 
    ## Run 397 stress 0.1016164 
    ## ... Procrustes: rmse 3.749995e-05  max resid 0.0001031336 
    ## ... Similar to previous best
    ## Run 398 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001047262  max resid 0.0002951189 
    ## ... Similar to previous best
    ## Run 399 stress 0.1356149 
    ## Run 400 stress 0.1530425 
    ## Run 401 stress 0.1016171 
    ## ... Procrustes: rmse 0.0006378995  max resid 0.001784258 
    ## ... Similar to previous best
    ## Run 402 stress 0.1016164 
    ## ... Procrustes: rmse 5.779348e-05  max resid 0.0001629946 
    ## ... Similar to previous best
    ## Run 403 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007216735  max resid 0.002017079 
    ## ... Similar to previous best
    ## Run 404 stress 0.1356144 
    ## Run 405 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004573087  max resid 0.001279962 
    ## ... Similar to previous best
    ## Run 406 stress 0.1356146 
    ## Run 407 stress 0.1016164 
    ## ... Procrustes: rmse 9.701412e-05  max resid 0.0002701626 
    ## ... Similar to previous best
    ## Run 408 stress 0.1016164 
    ## ... Procrustes: rmse 2.438443e-06  max resid 6.715463e-06 
    ## ... Similar to previous best
    ## Run 409 stress 0.135615 
    ## Run 410 stress 0.1341868 
    ## Run 411 stress 0.101617 
    ## ... Procrustes: rmse 0.0006587407  max resid 0.001840594 
    ## ... Similar to previous best
    ## Run 412 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002485155  max resid 0.0006967641 
    ## ... Similar to previous best
    ## Run 413 stress 0.1016164 
    ## ... Procrustes: rmse 1.089709e-05  max resid 2.044761e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.1016164 
    ## ... Procrustes: rmse 2.63484e-05  max resid 7.415302e-05 
    ## ... Similar to previous best
    ## Run 415 stress 0.1356148 
    ## Run 416 stress 0.1016164 
    ## ... Procrustes: rmse 3.568607e-06  max resid 1.03818e-05 
    ## ... Similar to previous best
    ## Run 417 stress 0.1016164 
    ## ... Procrustes: rmse 7.549084e-05  max resid 0.0002129278 
    ## ... Similar to previous best
    ## Run 418 stress 0.1356147 
    ## Run 419 stress 0.1341868 
    ## Run 420 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005788895  max resid 0.00161687 
    ## ... Similar to previous best
    ## Run 421 stress 0.1356146 
    ## Run 422 stress 0.1016164 
    ## ... Procrustes: rmse 3.497757e-05  max resid 9.867933e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.1341868 
    ## Run 424 stress 0.1016164 
    ## ... Procrustes: rmse 1.247474e-05  max resid 3.449038e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.1341868 
    ## Run 426 stress 0.1016164 
    ## ... Procrustes: rmse 3.022532e-05  max resid 8.51334e-05 
    ## ... Similar to previous best
    ## Run 427 stress 0.101617 
    ## ... Procrustes: rmse 0.0006487275  max resid 0.00181199 
    ## ... Similar to previous best
    ## Run 428 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003781751  max resid 0.001058066 
    ## ... Similar to previous best
    ## Run 429 stress 0.1016164 
    ## ... Procrustes: rmse 6.203842e-05  max resid 0.0001717991 
    ## ... Similar to previous best
    ## Run 430 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003763842  max resid 0.001051378 
    ## ... Similar to previous best
    ## Run 431 stress 0.1016164 
    ## ... Procrustes: rmse 8.879438e-06  max resid 2.240764e-05 
    ## ... Similar to previous best
    ## Run 432 stress 0.1356146 
    ## Run 433 stress 0.1016164 
    ## ... Procrustes: rmse 7.037102e-05  max resid 0.0001978692 
    ## ... Similar to previous best
    ## Run 434 stress 0.1016164 
    ## ... Procrustes: rmse 0.000112767  max resid 0.0003183681 
    ## ... Similar to previous best
    ## Run 435 stress 0.1341868 
    ## Run 436 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004001417  max resid 0.001116209 
    ## ... Similar to previous best
    ## Run 437 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001078695  max resid 0.0003045759 
    ## ... Similar to previous best
    ## Run 438 stress 0.1356149 
    ## Run 439 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005379989  max resid 0.001504084 
    ## ... Similar to previous best
    ## Run 440 stress 0.1016164 
    ## ... Procrustes: rmse 5.545807e-05  max resid 0.0001557077 
    ## ... Similar to previous best
    ## Run 441 stress 0.1341868 
    ## Run 442 stress 0.1016164 
    ## ... Procrustes: rmse 4.197741e-05  max resid 0.0001185573 
    ## ... Similar to previous best
    ## Run 443 stress 0.1341868 
    ## Run 444 stress 0.1016164 
    ## ... Procrustes: rmse 8.075351e-05  max resid 0.0002277369 
    ## ... Similar to previous best
    ## Run 445 stress 0.1016164 
    ## ... Procrustes: rmse 7.868642e-06  max resid 2.200515e-05 
    ## ... Similar to previous best
    ## Run 446 stress 0.1016164 
    ## ... Procrustes: rmse 4.886213e-05  max resid 0.0001372574 
    ## ... Similar to previous best
    ## Run 447 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003255962  max resid 0.0009116246 
    ## ... Similar to previous best
    ## Run 448 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001159826  max resid 0.0003269759 
    ## ... Similar to previous best
    ## Run 449 stress 0.1356145 
    ## Run 450 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001262524  max resid 0.0003546222 
    ## ... Similar to previous best
    ## Run 451 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001545176  max resid 0.0004333133 
    ## ... Similar to previous best
    ## Run 452 stress 0.1016164 
    ## ... Procrustes: rmse 3.51613e-05  max resid 9.90422e-05 
    ## ... Similar to previous best
    ## Run 453 stress 0.1356145 
    ## Run 454 stress 0.1341868 
    ## Run 455 stress 0.1016164 
    ## ... Procrustes: rmse 6.588934e-05  max resid 0.0001859926 
    ## ... Similar to previous best
    ## Run 456 stress 0.1341868 
    ## Run 457 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003158897  max resid 0.0008843749 
    ## ... Similar to previous best
    ## Run 458 stress 0.1016164 
    ## ... Procrustes: rmse 1.132169e-05  max resid 3.159821e-05 
    ## ... Similar to previous best
    ## Run 459 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003207319  max resid 0.000897971 
    ## ... Similar to previous best
    ## Run 460 stress 0.1016164 
    ## ... Procrustes: rmse 6.744917e-05  max resid 0.000190459 
    ## ... Similar to previous best
    ## Run 461 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004815287  max resid 0.001346022 
    ## ... Similar to previous best
    ## Run 462 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004896916  max resid 0.001369359 
    ## ... Similar to previous best
    ## Run 463 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006324812  max resid 0.001768285 
    ## ... Similar to previous best
    ## Run 464 stress 0.1341868 
    ## Run 465 stress 0.1016164 
    ## ... Procrustes: rmse 4.751398e-05  max resid 0.0001339974 
    ## ... Similar to previous best
    ## Run 466 stress 0.135615 
    ## Run 467 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001518507  max resid 0.0004266353 
    ## ... Similar to previous best
    ## Run 468 stress 0.1016164 
    ## ... Procrustes: rmse 4.881653e-05  max resid 0.000137329 
    ## ... Similar to previous best
    ## Run 469 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004883451  max resid 0.001366673 
    ## ... Similar to previous best
    ## Run 470 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001158417  max resid 0.0003269897 
    ## ... Similar to previous best
    ## Run 471 stress 0.1016164 
    ## ... Procrustes: rmse 0.000104463  max resid 0.0002930887 
    ## ... Similar to previous best
    ## Run 472 stress 0.3578479 
    ## Run 473 stress 0.101617 
    ## ... Procrustes: rmse 0.0007038603  max resid 0.001967133 
    ## ... Similar to previous best
    ## Run 474 stress 0.1356151 
    ## Run 475 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004608879  max resid 0.001289299 
    ## ... Similar to previous best
    ## Run 476 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004771394  max resid 0.001334833 
    ## ... Similar to previous best
    ## Run 477 stress 0.1341868 
    ## Run 478 stress 0.1341868 
    ## Run 479 stress 0.1016169 
    ## ... Procrustes: rmse 0.0004230962  max resid 0.001173253 
    ## ... Similar to previous best
    ## Run 480 stress 0.1341868 
    ## Run 481 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005185412  max resid 0.001449442 
    ## ... Similar to previous best
    ## Run 482 stress 0.3060031 
    ## Run 483 stress 0.1519033 
    ## Run 484 stress 0.1341868 
    ## Run 485 stress 0.3237824 
    ## Run 486 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003718373  max resid 0.001041203 
    ## ... Similar to previous best
    ## Run 487 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003775418  max resid 0.001056775 
    ## ... Similar to previous best
    ## Run 488 stress 0.1356148 
    ## Run 489 stress 0.1016164 
    ## ... Procrustes: rmse 1.834616e-05  max resid 5.166298e-05 
    ## ... Similar to previous best
    ## Run 490 stress 0.1356143 
    ## Run 491 stress 0.1016164 
    ## ... Procrustes: rmse 8.89714e-05  max resid 0.0002489133 
    ## ... Similar to previous best
    ## Run 492 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002915778  max resid 0.0008148573 
    ## ... Similar to previous best
    ## Run 493 stress 0.1341868 
    ## Run 494 stress 0.1016164 
    ## ... Procrustes: rmse 1.447359e-05  max resid 3.680558e-05 
    ## ... Similar to previous best
    ## Run 495 stress 0.1016167 
    ## ... Procrustes: rmse 0.000481461  max resid 0.001346093 
    ## ... Similar to previous best
    ## Run 496 stress 0.1356149 
    ## Run 497 stress 0.1530427 
    ## Run 498 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004906466  max resid 0.001372084 
    ## ... Similar to previous best
    ## Run 499 stress 0.1016164 
    ## ... Procrustes: rmse 3.175039e-05  max resid 8.925365e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.101617 
    ## ... Procrustes: rmse 0.0006647374  max resid 0.00185696 
    ## ... Similar to previous best
    ## *** Best solution repeated 184 times

``` r
# Stratified lakes and ocean sites
PD_beta_SO_NMDS <- metaMDS(PD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.04851417 
    ## Run 1 stress 0.0527641 
    ## Run 2 stress 0.04851417 
    ## ... New best solution
    ## ... Procrustes: rmse 2.667887e-05  max resid 4.479348e-05 
    ## ... Similar to previous best
    ## Run 3 stress 0.04938615 
    ## Run 4 stress 0.06400747 
    ## Run 5 stress 0.05605153 
    ## Run 6 stress 0.05474817 
    ## Run 7 stress 0.05216138 
    ## Run 8 stress 0.04851417 
    ## ... New best solution
    ## ... Procrustes: rmse 2.40589e-05  max resid 4.652575e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.05605148 
    ## Run 10 stress 0.0544779 
    ## Run 11 stress 0.05728063 
    ## Run 12 stress 0.04935943 
    ## Run 13 stress 0.04612768 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04986046  max resid 0.1396998 
    ## Run 14 stress 0.05607344 
    ## Run 15 stress 0.05607344 
    ## Run 16 stress 0.04938616 
    ## Run 17 stress 0.04938618 
    ## Run 18 stress 0.05643262 
    ## Run 19 stress 0.0561224 
    ## Run 20 stress 0.04938616 
    ## Run 21 stress 0.05593458 
    ## Run 22 stress 0.04612763 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002767382  max resid 0.0005530267 
    ## ... Similar to previous best
    ## Run 23 stress 0.05572688 
    ## Run 24 stress 0.04851419 
    ## Run 25 stress 0.05438942 
    ## Run 26 stress 0.04851418 
    ## Run 27 stress 0.05607347 
    ## Run 28 stress 0.04935944 
    ## Run 29 stress 0.05607345 
    ## Run 30 stress 0.05607344 
    ## Run 31 stress 0.04938617 
    ## Run 32 stress 0.05728055 
    ## Run 33 stress 0.05260056 
    ## Run 34 stress 0.05451493 
    ## Run 35 stress 0.04935943 
    ## Run 36 stress 0.04860423 
    ## Run 37 stress 0.06088973 
    ## Run 38 stress 0.05439285 
    ## Run 39 stress 0.04851417 
    ## Run 40 stress 0.05474821 
    ## Run 41 stress 0.05216138 
    ## Run 42 stress 0.05695856 
    ## Run 43 stress 0.05474815 
    ## Run 44 stress 0.0557269 
    ## Run 45 stress 0.05474817 
    ## Run 46 stress 0.05216141 
    ## Run 47 stress 0.0543894 
    ## Run 48 stress 0.04938616 
    ## Run 49 stress 0.05276412 
    ## Run 50 stress 0.04860422 
    ## Run 51 stress 0.05607344 
    ## Run 52 stress 0.04935951 
    ## Run 53 stress 0.05607344 
    ## Run 54 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001173806  max resid 0.0002310564 
    ## ... Similar to previous best
    ## Run 55 stress 0.05695855 
    ## Run 56 stress 0.05593458 
    ## Run 57 stress 0.06933824 
    ## Run 58 stress 0.04612757 
    ## ... Procrustes: rmse 2.168563e-05  max resid 3.737477e-05 
    ## ... Similar to previous best
    ## Run 59 stress 0.04612757 
    ## ... Procrustes: rmse 3.387089e-05  max resid 6.422863e-05 
    ## ... Similar to previous best
    ## Run 60 stress 0.04851417 
    ## Run 61 stress 0.05605154 
    ## Run 62 stress 0.05728081 
    ## Run 63 stress 0.06291429 
    ## Run 64 stress 0.04938618 
    ## Run 65 stress 0.05216128 
    ## Run 66 stress 0.05643261 
    ## Run 67 stress 0.04851421 
    ## Run 68 stress 0.05612247 
    ## Run 69 stress 0.05474819 
    ## Run 70 stress 0.05474819 
    ## Run 71 stress 0.06291425 
    ## Run 72 stress 0.05137964 
    ## Run 73 stress 0.04860423 
    ## Run 74 stress 0.05959898 
    ## Run 75 stress 0.04938616 
    ## Run 76 stress 0.05607344 
    ## Run 77 stress 0.0547482 
    ## Run 78 stress 0.04612758 
    ## ... Procrustes: rmse 4.812829e-05  max resid 9.495737e-05 
    ## ... Similar to previous best
    ## Run 79 stress 0.05855035 
    ## Run 80 stress 0.04612759 
    ## ... Procrustes: rmse 7.095144e-05  max resid 0.0001409655 
    ## ... Similar to previous best
    ## Run 81 stress 0.0618505 
    ## Run 82 stress 0.04851417 
    ## Run 83 stress 0.04612757 
    ## ... Procrustes: rmse 2.102574e-05  max resid 3.724526e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.04612759 
    ## ... Procrustes: rmse 7.58998e-05  max resid 0.0001526463 
    ## ... Similar to previous best
    ## Run 85 stress 0.05474817 
    ## Run 86 stress 0.06205754 
    ## Run 87 stress 0.0557269 
    ## Run 88 stress 0.04938616 
    ## Run 89 stress 0.05842785 
    ## Run 90 stress 0.0543928 
    ## Run 91 stress 0.04860423 
    ## Run 92 stress 0.04938615 
    ## Run 93 stress 0.05433178 
    ## Run 94 stress 0.05605147 
    ## Run 95 stress 0.04612759 
    ## ... Procrustes: rmse 6.111241e-05  max resid 0.0001171685 
    ## ... Similar to previous best
    ## Run 96 stress 0.05474824 
    ## Run 97 stress 0.06933848 
    ## Run 98 stress 0.05474814 
    ## Run 99 stress 0.06205755 
    ## Run 100 stress 0.04935944 
    ## Run 101 stress 0.04938615 
    ## Run 102 stress 0.05474817 
    ## Run 103 stress 0.04935947 
    ## Run 104 stress 0.0543928 
    ## Run 105 stress 0.05695871 
    ## Run 106 stress 0.05433172 
    ## Run 107 stress 0.04935946 
    ## Run 108 stress 0.05593458 
    ## Run 109 stress 0.05474818 
    ## Run 110 stress 0.04851417 
    ## Run 111 stress 0.05728068 
    ## Run 112 stress 0.06185052 
    ## Run 113 stress 0.04938617 
    ## Run 114 stress 0.05593469 
    ## Run 115 stress 0.04612761 
    ## ... Procrustes: rmse 0.0001057155  max resid 0.0002110325 
    ## ... Similar to previous best
    ## Run 116 stress 0.05474813 
    ## Run 117 stress 0.04860422 
    ## Run 118 stress 0.04851417 
    ## Run 119 stress 0.05695867 
    ## Run 120 stress 0.05695864 
    ## Run 121 stress 0.0543928 
    ## Run 122 stress 0.04938616 
    ## Run 123 stress 0.05216161 
    ## Run 124 stress 0.05728053 
    ## Run 125 stress 0.05216133 
    ## Run 126 stress 0.05605148 
    ## Run 127 stress 0.04938616 
    ## Run 128 stress 0.04851417 
    ## Run 129 stress 0.05593454 
    ## Run 130 stress 0.06400749 
    ## Run 131 stress 0.05276411 
    ## Run 132 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001227209  max resid 0.0002498136 
    ## ... Similar to previous best
    ## Run 133 stress 0.05474821 
    ## Run 134 stress 0.04851417 
    ## Run 135 stress 0.06185048 
    ## Run 136 stress 0.05605145 
    ## Run 137 stress 0.05605149 
    ## Run 138 stress 0.04938616 
    ## Run 139 stress 0.05695858 
    ## Run 140 stress 0.05643261 
    ## Run 141 stress 0.05488104 
    ## Run 142 stress 0.05695861 
    ## Run 143 stress 0.05605146 
    ## Run 144 stress 0.05572688 
    ## Run 145 stress 0.05474814 
    ## Run 146 stress 0.05572688 
    ## Run 147 stress 0.05959896 
    ## Run 148 stress 0.04935946 
    ## Run 149 stress 0.05643264 
    ## Run 150 stress 0.05474817 
    ## Run 151 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 1.357543e-05  max resid 3.02077e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.05216163 
    ## Run 153 stress 0.04938616 
    ## Run 154 stress 0.05137965 
    ## Run 155 stress 0.05474815 
    ## Run 156 stress 0.0560515 
    ## Run 157 stress 0.05605152 
    ## Run 158 stress 0.04938617 
    ## Run 159 stress 0.04938615 
    ## Run 160 stress 0.05216162 
    ## Run 161 stress 0.05643262 
    ## Run 162 stress 0.0526005 
    ## Run 163 stress 0.04938616 
    ## Run 164 stress 0.05959896 
    ## Run 165 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001186723  max resid 0.00023568 
    ## ... Similar to previous best
    ## Run 166 stress 0.05438938 
    ## Run 167 stress 0.06185051 
    ## Run 168 stress 0.05728055 
    ## Run 169 stress 0.04935949 
    ## Run 170 stress 0.06185042 
    ## Run 171 stress 0.05137966 
    ## Run 172 stress 0.04860423 
    ## Run 173 stress 0.05605153 
    ## Run 174 stress 0.04851417 
    ## Run 175 stress 0.05605147 
    ## Run 176 stress 0.04851418 
    ## Run 177 stress 0.05842782 
    ## Run 178 stress 0.05855038 
    ## Run 179 stress 0.05643264 
    ## Run 180 stress 0.04612759 
    ## ... Procrustes: rmse 7.273964e-05  max resid 0.0001453915 
    ## ... Similar to previous best
    ## Run 181 stress 0.06291426 
    ## Run 182 stress 0.04938616 
    ## Run 183 stress 0.04851417 
    ## Run 184 stress 0.05959897 
    ## Run 185 stress 0.05728068 
    ## Run 186 stress 0.05607344 
    ## Run 187 stress 0.04612758 
    ## ... Procrustes: rmse 5.727032e-05  max resid 0.0001142861 
    ## ... Similar to previous best
    ## Run 188 stress 0.04851417 
    ## Run 189 stress 0.06205756 
    ## Run 190 stress 0.05855034 
    ## Run 191 stress 0.05695862 
    ## Run 192 stress 0.06291425 
    ## Run 193 stress 0.06933818 
    ## Run 194 stress 0.05855034 
    ## Run 195 stress 0.0527641 
    ## Run 196 stress 0.04851417 
    ## Run 197 stress 0.05451497 
    ## Run 198 stress 0.04935943 
    ## Run 199 stress 0.3578482 
    ## Run 200 stress 0.05607345 
    ## Run 201 stress 0.04851418 
    ## Run 202 stress 0.05438944 
    ## Run 203 stress 0.05959895 
    ## Run 204 stress 0.04938616 
    ## Run 205 stress 0.05731516 
    ## Run 206 stress 0.05451494 
    ## Run 207 stress 0.05842784 
    ## Run 208 stress 0.05438947 
    ## Run 209 stress 0.0543928 
    ## Run 210 stress 0.05605147 
    ## Run 211 stress 0.06400729 
    ## Run 212 stress 0.05728067 
    ## Run 213 stress 0.05605145 
    ## Run 214 stress 0.0461276 
    ## ... Procrustes: rmse 7.52865e-05  max resid 0.0001469153 
    ## ... Similar to previous best
    ## Run 215 stress 0.04612765 
    ## ... Procrustes: rmse 0.0001317445  max resid 0.0002676666 
    ## ... Similar to previous best
    ## Run 216 stress 0.05451493 
    ## Run 217 stress 0.04851417 
    ## Run 218 stress 0.0527641 
    ## Run 219 stress 0.04860422 
    ## Run 220 stress 0.05137965 
    ## Run 221 stress 0.05438941 
    ## Run 222 stress 0.05572688 
    ## Run 223 stress 0.0585503 
    ## Run 224 stress 0.04938616 
    ## Run 225 stress 0.05605157 
    ## Run 226 stress 0.05607346 
    ## Run 227 stress 0.04851417 
    ## Run 228 stress 0.05451493 
    ## Run 229 stress 0.05474814 
    ## Run 230 stress 0.05605146 
    ## Run 231 stress 0.04612762 
    ## ... Procrustes: rmse 0.0001103859  max resid 0.0002224206 
    ## ... Similar to previous best
    ## Run 232 stress 0.04938616 
    ## Run 233 stress 0.05728091 
    ## Run 234 stress 0.0526005 
    ## Run 235 stress 0.04938616 
    ## Run 236 stress 0.05438945 
    ## Run 237 stress 0.0485142 
    ## Run 238 stress 0.04860422 
    ## Run 239 stress 0.05488106 
    ## Run 240 stress 0.04851417 
    ## Run 241 stress 0.05216125 
    ## Run 242 stress 0.05474814 
    ## Run 243 stress 0.04938617 
    ## Run 244 stress 0.05643261 
    ## Run 245 stress 0.05605145 
    ## Run 246 stress 0.04938618 
    ## Run 247 stress 0.3125127 
    ## Run 248 stress 0.05612244 
    ## Run 249 stress 0.04851417 
    ## Run 250 stress 0.05959898 
    ## Run 251 stress 0.06088974 
    ## Run 252 stress 0.05612248 
    ## Run 253 stress 0.05433175 
    ## Run 254 stress 0.05728053 
    ## Run 255 stress 0.05842784 
    ## Run 256 stress 0.04938617 
    ## Run 257 stress 0.05433178 
    ## Run 258 stress 0.05605146 
    ## Run 259 stress 0.05260061 
    ## Run 260 stress 0.06088978 
    ## Run 261 stress 0.04851418 
    ## Run 262 stress 0.05607348 
    ## Run 263 stress 0.05276411 
    ## Run 264 stress 0.05643261 
    ## Run 265 stress 0.04935948 
    ## Run 266 stress 0.04860422 
    ## Run 267 stress 0.04851417 
    ## Run 268 stress 0.04938616 
    ## Run 269 stress 0.05216128 
    ## Run 270 stress 0.05605157 
    ## Run 271 stress 0.05842794 
    ## Run 272 stress 0.05855033 
    ## Run 273 stress 0.05451494 
    ## Run 274 stress 0.04938615 
    ## Run 275 stress 0.05728058 
    ## Run 276 stress 0.05438952 
    ## Run 277 stress 0.05607344 
    ## Run 278 stress 0.04612765 
    ## ... Procrustes: rmse 0.0001440068  max resid 0.0002910602 
    ## ... Similar to previous best
    ## Run 279 stress 0.05474823 
    ## Run 280 stress 0.05643261 
    ## Run 281 stress 0.04938615 
    ## Run 282 stress 0.0461276 
    ## ... Procrustes: rmse 9.186597e-05  max resid 0.0001875852 
    ## ... Similar to previous best
    ## Run 283 stress 0.04938615 
    ## Run 284 stress 0.04938617 
    ## Run 285 stress 0.04851418 
    ## Run 286 stress 0.05474813 
    ## Run 287 stress 0.04851417 
    ## Run 288 stress 0.05607344 
    ## Run 289 stress 0.05842785 
    ## Run 290 stress 0.05260047 
    ## Run 291 stress 0.05643261 
    ## Run 292 stress 0.05855029 
    ## Run 293 stress 0.04935948 
    ## Run 294 stress 0.04860422 
    ## Run 295 stress 0.05612243 
    ## Run 296 stress 0.05959895 
    ## Run 297 stress 0.05260047 
    ## Run 298 stress 0.05433177 
    ## Run 299 stress 0.06400728 
    ## Run 300 stress 0.05605147 
    ## Run 301 stress 0.04860423 
    ## Run 302 stress 0.05605145 
    ## Run 303 stress 0.05643261 
    ## Run 304 stress 0.05572688 
    ## Run 305 stress 0.05438946 
    ## Run 306 stress 0.05474822 
    ## Run 307 stress 0.04851417 
    ## Run 308 stress 0.05605149 
    ## Run 309 stress 0.04938616 
    ## Run 310 stress 0.0547482 
    ## Run 311 stress 0.0585503 
    ## Run 312 stress 0.06185047 
    ## Run 313 stress 0.05488104 
    ## Run 314 stress 0.04851417 
    ## Run 315 stress 0.05260053 
    ## Run 316 stress 0.05474816 
    ## Run 317 stress 0.0527641 
    ## Run 318 stress 0.04851418 
    ## Run 319 stress 0.04938615 
    ## Run 320 stress 0.05728056 
    ## Run 321 stress 0.05572692 
    ## Run 322 stress 0.05438942 
    ## Run 323 stress 0.04851417 
    ## Run 324 stress 0.05572688 
    ## Run 325 stress 0.04851417 
    ## Run 326 stress 0.05488104 
    ## Run 327 stress 0.05695875 
    ## Run 328 stress 0.04851417 
    ## Run 329 stress 0.05605145 
    ## Run 330 stress 0.06291426 
    ## Run 331 stress 0.0493862 
    ## Run 332 stress 0.04938616 
    ## Run 333 stress 0.05842781 
    ## Run 334 stress 0.06088973 
    ## Run 335 stress 0.05643263 
    ## Run 336 stress 0.05451493 
    ## Run 337 stress 0.05439279 
    ## Run 338 stress 0.06185045 
    ## Run 339 stress 0.04938617 
    ## Run 340 stress 0.04938617 
    ## Run 341 stress 0.06291426 
    ## Run 342 stress 0.06291425 
    ## Run 343 stress 0.05433183 
    ## Run 344 stress 0.05488103 
    ## Run 345 stress 0.05276411 
    ## Run 346 stress 0.04938616 
    ## Run 347 stress 0.05439282 
    ## Run 348 stress 0.05451493 
    ## Run 349 stress 0.04851417 
    ## Run 350 stress 0.0561224 
    ## Run 351 stress 0.0561224 
    ## Run 352 stress 0.05137965 
    ## Run 353 stress 0.05643264 
    ## Run 354 stress 0.04612759 
    ## ... Procrustes: rmse 6.253543e-05  max resid 0.0001244445 
    ## ... Similar to previous best
    ## Run 355 stress 0.05474822 
    ## Run 356 stress 0.05605152 
    ## Run 357 stress 0.05695882 
    ## Run 358 stress 0.05474822 
    ## Run 359 stress 0.05643261 
    ## Run 360 stress 0.05728055 
    ## Run 361 stress 0.04938616 
    ## Run 362 stress 0.05728066 
    ## Run 363 stress 0.05572688 
    ## Run 364 stress 0.05438938 
    ## Run 365 stress 0.05488107 
    ## Run 366 stress 0.0543895 
    ## Run 367 stress 0.05728078 
    ## Run 368 stress 0.0461276 
    ## ... Procrustes: rmse 8.683153e-05  max resid 0.0001731122 
    ## ... Similar to previous best
    ## Run 369 stress 0.05607344 
    ## Run 370 stress 0.05855035 
    ## Run 371 stress 0.05607346 
    ## Run 372 stress 0.0561225 
    ## Run 373 stress 0.05433176 
    ## Run 374 stress 0.05959898 
    ## Run 375 stress 0.05474823 
    ## Run 376 stress 0.05216142 
    ## Run 377 stress 0.04860423 
    ## Run 378 stress 0.04612761 
    ## ... Procrustes: rmse 9.927795e-05  max resid 0.0002007141 
    ## ... Similar to previous best
    ## Run 379 stress 0.05474819 
    ## Run 380 stress 0.05612244 
    ## Run 381 stress 0.05643261 
    ## Run 382 stress 0.05137964 
    ## Run 383 stress 0.04851418 
    ## Run 384 stress 0.04938621 
    ## Run 385 stress 0.04851419 
    ## Run 386 stress 0.04935946 
    ## Run 387 stress 0.04851417 
    ## Run 388 stress 0.0585503 
    ## Run 389 stress 0.04935944 
    ## Run 390 stress 0.04612757 
    ## ... Procrustes: rmse 1.479319e-05  max resid 2.187482e-05 
    ## ... Similar to previous best
    ## Run 391 stress 0.05731513 
    ## Run 392 stress 0.05605152 
    ## Run 393 stress 0.04612757 
    ## ... Procrustes: rmse 1.14367e-05  max resid 1.650716e-05 
    ## ... Similar to previous best
    ## Run 394 stress 0.05260049 
    ## Run 395 stress 0.05438943 
    ## Run 396 stress 0.05137971 
    ## Run 397 stress 0.0608897 
    ## Run 398 stress 0.05728063 
    ## Run 399 stress 0.05607344 
    ## Run 400 stress 0.04851417 
    ## Run 401 stress 0.0608897 
    ## Run 402 stress 0.05451494 
    ## Run 403 stress 0.05605145 
    ## Run 404 stress 0.04851418 
    ## Run 405 stress 0.04938619 
    ## Run 406 stress 0.06400727 
    ## Run 407 stress 0.05643261 
    ## Run 408 stress 0.05276412 
    ## Run 409 stress 0.0461277 
    ## ... Procrustes: rmse 0.000187846  max resid 0.0003822884 
    ## ... Similar to previous best
    ## Run 410 stress 0.04938619 
    ## Run 411 stress 0.04938616 
    ## Run 412 stress 0.05612242 
    ## Run 413 stress 0.05607345 
    ## Run 414 stress 0.05216128 
    ## Run 415 stress 0.05728053 
    ## Run 416 stress 0.05959897 
    ## Run 417 stress 0.05605153 
    ## Run 418 stress 0.05276412 
    ## Run 419 stress 0.04938617 
    ## Run 420 stress 0.04612757 
    ## ... Procrustes: rmse 1.604625e-05  max resid 2.453299e-05 
    ## ... Similar to previous best
    ## Run 421 stress 0.05643262 
    ## Run 422 stress 0.04851417 
    ## Run 423 stress 0.06205753 
    ## Run 424 stress 0.05216142 
    ## Run 425 stress 0.05451493 
    ## Run 426 stress 0.04851417 
    ## Run 427 stress 0.05438948 
    ## Run 428 stress 0.0461276 
    ## ... Procrustes: rmse 6.943861e-05  max resid 0.000134563 
    ## ... Similar to previous best
    ## Run 429 stress 0.05607347 
    ## Run 430 stress 0.0527641 
    ## Run 431 stress 0.04612758 
    ## ... Procrustes: rmse 6.363922e-05  max resid 0.0001294575 
    ## ... Similar to previous best
    ## Run 432 stress 0.05438942 
    ## Run 433 stress 0.05488103 
    ## Run 434 stress 0.05216141 
    ## Run 435 stress 0.04938615 
    ## Run 436 stress 0.05593456 
    ## Run 437 stress 0.05728054 
    ## Run 438 stress 0.0585503 
    ## Run 439 stress 0.05607345 
    ## Run 440 stress 0.05728058 
    ## Run 441 stress 0.04612757 
    ## ... Procrustes: rmse 3.53468e-05  max resid 6.855623e-05 
    ## ... Similar to previous best
    ## Run 442 stress 0.05474814 
    ## Run 443 stress 0.04935942 
    ## Run 444 stress 0.05607344 
    ## Run 445 stress 0.04612759 
    ## ... Procrustes: rmse 8.442655e-05  max resid 0.0001731988 
    ## ... Similar to previous best
    ## Run 446 stress 0.04938618 
    ## Run 447 stress 0.04938615 
    ## Run 448 stress 0.05474814 
    ## Run 449 stress 0.0561224 
    ## Run 450 stress 0.05605146 
    ## Run 451 stress 0.04860422 
    ## Run 452 stress 0.04851418 
    ## Run 453 stress 0.06291429 
    ## Run 454 stress 0.0569586 
    ## Run 455 stress 0.05605149 
    ## Run 456 stress 0.05605145 
    ## Run 457 stress 0.05216128 
    ## Run 458 stress 0.05728058 
    ## Run 459 stress 0.05605151 
    ## Run 460 stress 0.04938616 
    ## Run 461 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001097302  max resid 0.0002179311 
    ## ... Similar to previous best
    ## Run 462 stress 0.04938615 
    ## Run 463 stress 0.04612757 
    ## ... Procrustes: rmse 2.334291e-05  max resid 4.860158e-05 
    ## ... Similar to previous best
    ## Run 464 stress 0.0547482 
    ## Run 465 stress 0.0572806 
    ## Run 466 stress 0.04851417 
    ## Run 467 stress 0.05695875 
    ## Run 468 stress 0.04612771 
    ## ... Procrustes: rmse 0.0001857178  max resid 0.0003772882 
    ## ... Similar to previous best
    ## Run 469 stress 0.05855034 
    ## Run 470 stress 0.06400741 
    ## Run 471 stress 0.05607347 
    ## Run 472 stress 0.05607345 
    ## Run 473 stress 0.04860423 
    ## Run 474 stress 0.04851418 
    ## Run 475 stress 0.05695858 
    ## Run 476 stress 0.05593461 
    ## Run 477 stress 0.04938618 
    ## Run 478 stress 0.05643264 
    ## Run 479 stress 0.05728052 
    ## Run 480 stress 0.04612757 
    ## ... Procrustes: rmse 1.367335e-05  max resid 2.679597e-05 
    ## ... Similar to previous best
    ## Run 481 stress 0.04935942 
    ## Run 482 stress 0.04860423 
    ## Run 483 stress 0.05433182 
    ## Run 484 stress 0.05474824 
    ## Run 485 stress 0.05451493 
    ## Run 486 stress 0.05474817 
    ## Run 487 stress 0.05488105 
    ## Run 488 stress 0.06291431 
    ## Run 489 stress 0.05593451 
    ## Run 490 stress 0.05137964 
    ## Run 491 stress 0.05439282 
    ## Run 492 stress 0.06205752 
    ## Run 493 stress 0.05451493 
    ## Run 494 stress 0.05842787 
    ## Run 495 stress 0.05216152 
    ## Run 496 stress 0.04851418 
    ## Run 497 stress 0.05728065 
    ## Run 498 stress 0.05137965 
    ## Run 499 stress 0.04851417 
    ## Run 500 stress 0.05137965 
    ## *** Best solution repeated 24 times

``` r
# Stratified lakes
PD_beta_S_NMDS <- metaMDS(PD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.02693829 
    ## Run 1 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001402344  max resid 0.0002380092 
    ## ... Similar to previous best
    ## Run 2 stress 0.07547299 
    ## Run 3 stress 0.02693838 
    ## ... Procrustes: rmse 0.000196303  max resid 0.000335739 
    ## ... Similar to previous best
    ## Run 4 stress 0.07374808 
    ## Run 5 stress 0.08382828 
    ## Run 6 stress 0.1862386 
    ## Run 7 stress 0.02693832 
    ## ... Procrustes: rmse 7.491169e-05  max resid 0.0001279208 
    ## ... Similar to previous best
    ## Run 8 stress 0.1756177 
    ## Run 9 stress 0.0269383 
    ## ... Procrustes: rmse 2.840091e-05  max resid 4.862505e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001883563  max resid 0.0003220647 
    ## ... Similar to previous best
    ## Run 11 stress 0.073748 
    ## Run 12 stress 0.07374804 
    ## Run 13 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001598703  max resid 0.0002741502 
    ## ... Similar to previous best
    ## Run 14 stress 0.2714151 
    ## Run 15 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001000529  max resid 0.0001704148 
    ## ... Similar to previous best
    ## Run 16 stress 0.2852158 
    ## Run 17 stress 0.02693828 
    ## ... New best solution
    ## ... Procrustes: rmse 1.653729e-05  max resid 2.796864e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.02693841 
    ## ... Procrustes: rmse 0.000201537  max resid 0.000344398 
    ## ... Similar to previous best
    ## Run 19 stress 0.07374801 
    ## Run 20 stress 0.07374802 
    ## Run 21 stress 0.07547299 
    ## Run 22 stress 0.179905 
    ## Run 23 stress 0.0737481 
    ## Run 24 stress 0.179905 
    ## Run 25 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001677653  max resid 0.0002878945 
    ## ... Similar to previous best
    ## Run 26 stress 0.02693829 
    ## ... Procrustes: rmse 4.24885e-05  max resid 7.253951e-05 
    ## ... Similar to previous best
    ## Run 27 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001709668  max resid 0.0002923885 
    ## ... Similar to previous best
    ## Run 28 stress 0.075473 
    ## Run 29 stress 0.1889141 
    ## Run 30 stress 0.08382832 
    ## Run 31 stress 0.08382872 
    ## Run 32 stress 0.1756177 
    ## Run 33 stress 0.07547299 
    ## Run 34 stress 0.02693832 
    ## ... Procrustes: rmse 9.170727e-05  max resid 0.000156475 
    ## ... Similar to previous best
    ## Run 35 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001489474  max resid 0.0002545607 
    ## ... Similar to previous best
    ## Run 36 stress 0.02693838 
    ## ... Procrustes: rmse 0.000153046  max resid 0.0002624568 
    ## ... Similar to previous best
    ## Run 37 stress 0.1862392 
    ## Run 38 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001602821  max resid 0.000275026 
    ## ... Similar to previous best
    ## Run 39 stress 0.07374802 
    ## Run 40 stress 0.02693833 
    ## ... Procrustes: rmse 9.703249e-05  max resid 0.0001663526 
    ## ... Similar to previous best
    ## Run 41 stress 0.07374805 
    ## Run 42 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001308574  max resid 0.0002236381 
    ## ... Similar to previous best
    ## Run 43 stress 0.07374808 
    ## Run 44 stress 0.07547302 
    ## Run 45 stress 0.02693832 
    ## ... Procrustes: rmse 8.846149e-05  max resid 0.0001494531 
    ## ... Similar to previous best
    ## Run 46 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001531256  max resid 0.0002625755 
    ## ... Similar to previous best
    ## Run 47 stress 0.07374801 
    ## Run 48 stress 0.1862387 
    ## Run 49 stress 0.0838285 
    ## Run 50 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001597871  max resid 0.0002719564 
    ## ... Similar to previous best
    ## Run 51 stress 0.075473 
    ## Run 52 stress 0.07374804 
    ## Run 53 stress 0.02693839 
    ## ... Procrustes: rmse 0.000159635  max resid 0.0002694649 
    ## ... Similar to previous best
    ## Run 54 stress 0.08382837 
    ## Run 55 stress 0.1756177 
    ## Run 56 stress 0.1862387 
    ## Run 57 stress 0.1900871 
    ## Run 58 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001339388  max resid 0.0002283637 
    ## ... Similar to previous best
    ## Run 59 stress 0.2714043 
    ## Run 60 stress 0.02693833 
    ## ... Procrustes: rmse 8.654782e-05  max resid 0.0001464399 
    ## ... Similar to previous best
    ## Run 61 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001893376  max resid 0.0003252669 
    ## ... Similar to previous best
    ## Run 62 stress 0.02693831 
    ## ... Procrustes: rmse 6.749363e-05  max resid 0.0001144321 
    ## ... Similar to previous best
    ## Run 63 stress 0.07547299 
    ## Run 64 stress 0.07374803 
    ## Run 65 stress 0.1756177 
    ## Run 66 stress 0.08382835 
    ## Run 67 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001443239  max resid 0.0002470061 
    ## ... Similar to previous best
    ## Run 68 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001251816  max resid 0.0002146034 
    ## ... Similar to previous best
    ## Run 69 stress 0.179905 
    ## Run 70 stress 0.02693828 
    ## ... New best solution
    ## ... Procrustes: rmse 1.315764e-05  max resid 2.262523e-05 
    ## ... Similar to previous best
    ## Run 71 stress 0.09391896 
    ## Run 72 stress 0.0269383 
    ## ... Procrustes: rmse 6.851737e-05  max resid 0.0001169724 
    ## ... Similar to previous best
    ## Run 73 stress 0.02693837 
    ## ... Procrustes: rmse 0.000146221  max resid 0.0002478909 
    ## ... Similar to previous best
    ## Run 74 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001669853  max resid 0.0002848471 
    ## ... Similar to previous best
    ## Run 75 stress 0.07374805 
    ## Run 76 stress 0.179905 
    ## Run 77 stress 0.02693831 
    ## ... Procrustes: rmse 8.994272e-05  max resid 0.000153578 
    ## ... Similar to previous best
    ## Run 78 stress 0.07374804 
    ## Run 79 stress 0.0737481 
    ## Run 80 stress 0.02693831 
    ## ... Procrustes: rmse 9.177315e-05  max resid 0.0001574311 
    ## ... Similar to previous best
    ## Run 81 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001061083  max resid 0.0001819079 
    ## ... Similar to previous best
    ## Run 82 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001475501  max resid 0.0002509002 
    ## ... Similar to previous best
    ## Run 83 stress 0.07547299 
    ## Run 84 stress 0.0269383 
    ## ... Procrustes: rmse 7.012084e-05  max resid 0.0001206137 
    ## ... Similar to previous best
    ## Run 85 stress 0.07547299 
    ## Run 86 stress 0.07547299 
    ## Run 87 stress 0.2273964 
    ## Run 88 stress 0.07374807 
    ## Run 89 stress 0.1756177 
    ## Run 90 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001233253  max resid 0.0002107733 
    ## ... Similar to previous best
    ## Run 91 stress 0.075473 
    ## Run 92 stress 0.07547299 
    ## Run 93 stress 0.07374807 
    ## Run 94 stress 0.07374803 
    ## Run 95 stress 0.08382871 
    ## Run 96 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001266372  max resid 0.0002193814 
    ## ... Similar to previous best
    ## Run 97 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001192189  max resid 0.0002028029 
    ## ... Similar to previous best
    ## Run 98 stress 0.179905 
    ## Run 99 stress 0.07547301 
    ## Run 100 stress 0.07374802 
    ## Run 101 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001106587  max resid 0.0001876862 
    ## ... Similar to previous best
    ## Run 102 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001505418  max resid 0.0002583842 
    ## ... Similar to previous best
    ## Run 103 stress 0.07374809 
    ## Run 104 stress 0.09391925 
    ## Run 105 stress 0.075473 
    ## Run 106 stress 0.1756177 
    ## Run 107 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001216548  max resid 0.000206849 
    ## ... Similar to previous best
    ## Run 108 stress 0.07374802 
    ## Run 109 stress 0.02693832 
    ## ... Procrustes: rmse 9.80675e-05  max resid 0.0001687998 
    ## ... Similar to previous best
    ## Run 110 stress 0.179905 
    ## Run 111 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001378716  max resid 0.0002360899 
    ## ... Similar to previous best
    ## Run 112 stress 0.08382851 
    ## Run 113 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001039193  max resid 0.0001944682 
    ## ... Similar to previous best
    ## Run 114 stress 0.07547302 
    ## Run 115 stress 0.02693831 
    ## ... Procrustes: rmse 6.391661e-05  max resid 0.0001078197 
    ## ... Similar to previous best
    ## Run 116 stress 0.1889141 
    ## Run 117 stress 0.1889141 
    ## Run 118 stress 0.07547299 
    ## Run 119 stress 0.0269383 
    ## ... Procrustes: rmse 2.710557e-05  max resid 4.675034e-05 
    ## ... Similar to previous best
    ## Run 120 stress 0.08382865 
    ## Run 121 stress 0.02693829 
    ## ... Procrustes: rmse 1.65156e-05  max resid 2.704048e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.02693832 
    ## ... Procrustes: rmse 9.703513e-05  max resid 0.0001665591 
    ## ... Similar to previous best
    ## Run 123 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001407395  max resid 0.0002406604 
    ## ... Similar to previous best
    ## Run 124 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001171074  max resid 0.0002008386 
    ## ... Similar to previous best
    ## Run 125 stress 0.1889141 
    ## Run 126 stress 0.07374804 
    ## Run 127 stress 0.1889143 
    ## Run 128 stress 0.08382841 
    ## Run 129 stress 0.0269383 
    ## ... Procrustes: rmse 4.495536e-05  max resid 7.584867e-05 
    ## ... Similar to previous best
    ## Run 130 stress 0.07374803 
    ## Run 131 stress 0.07374802 
    ## Run 132 stress 0.07547299 
    ## Run 133 stress 0.07374801 
    ## Run 134 stress 0.07374803 
    ## Run 135 stress 0.08382849 
    ## Run 136 stress 0.1889141 
    ## Run 137 stress 0.179905 
    ## Run 138 stress 0.07547299 
    ## Run 139 stress 0.07374804 
    ## Run 140 stress 0.02693831 
    ## ... Procrustes: rmse 8.590079e-05  max resid 0.0001465002 
    ## ... Similar to previous best
    ## Run 141 stress 0.07374811 
    ## Run 142 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001767069  max resid 0.0003019169 
    ## ... Similar to previous best
    ## Run 143 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001657731  max resid 0.0002820403 
    ## ... Similar to previous best
    ## Run 144 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001207379  max resid 0.0002070993 
    ## ... Similar to previous best
    ## Run 145 stress 0.073748 
    ## Run 146 stress 0.08382827 
    ## Run 147 stress 0.07374801 
    ## Run 148 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001153534  max resid 0.0001979659 
    ## ... Similar to previous best
    ## Run 149 stress 0.07374802 
    ## Run 150 stress 0.09391918 
    ## Run 151 stress 0.02693829 
    ## ... Procrustes: rmse 6.866998e-06  max resid 1.144482e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.07547299 
    ## Run 153 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001228183  max resid 0.0002107356 
    ## ... Similar to previous best
    ## Run 154 stress 0.07547304 
    ## Run 155 stress 0.07547299 
    ## Run 156 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001629242  max resid 0.0002804016 
    ## ... Similar to previous best
    ## Run 157 stress 0.07547301 
    ## Run 158 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001492647  max resid 0.0002557843 
    ## ... Similar to previous best
    ## Run 159 stress 0.2852146 
    ## Run 160 stress 0.179905 
    ## Run 161 stress 0.08382834 
    ## Run 162 stress 0.1756177 
    ## Run 163 stress 0.0269383 
    ## ... Procrustes: rmse 6.24257e-05  max resid 0.0001069344 
    ## ... Similar to previous best
    ## Run 164 stress 0.1756177 
    ## Run 165 stress 0.075473 
    ## Run 166 stress 0.179905 
    ## Run 167 stress 0.1862389 
    ## Run 168 stress 0.07374806 
    ## Run 169 stress 0.1862385 
    ## Run 170 stress 0.02693831 
    ## ... Procrustes: rmse 8.082565e-05  max resid 0.0001373358 
    ## ... Similar to previous best
    ## Run 171 stress 0.02693841 
    ## ... Procrustes: rmse 0.0001894122  max resid 0.0003256741 
    ## ... Similar to previous best
    ## Run 172 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001267987  max resid 0.0002178808 
    ## ... Similar to previous best
    ## Run 173 stress 0.07547299 
    ## Run 174 stress 0.1799052 
    ## Run 175 stress 0.179905 
    ## Run 176 stress 0.179905 
    ## Run 177 stress 0.02693832 
    ## ... Procrustes: rmse 7.864215e-05  max resid 0.0001356791 
    ## ... Similar to previous best
    ## Run 178 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001826568  max resid 0.0003128338 
    ## ... Similar to previous best
    ## Run 179 stress 0.07547301 
    ## Run 180 stress 0.179905 
    ## Run 181 stress 0.1862388 
    ## Run 182 stress 0.1889142 
    ## Run 183 stress 0.1756177 
    ## Run 184 stress 0.1756177 
    ## Run 185 stress 0.0269383 
    ## ... Procrustes: rmse 5.355711e-05  max resid 9.081556e-05 
    ## ... Similar to previous best
    ## Run 186 stress 0.07374803 
    ## Run 187 stress 0.07374804 
    ## Run 188 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001467421  max resid 0.0002517239 
    ## ... Similar to previous best
    ## Run 189 stress 0.07374806 
    ## Run 190 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001072481  max resid 0.0001843478 
    ## ... Similar to previous best
    ## Run 191 stress 0.1889141 
    ## Run 192 stress 0.0737481 
    ## Run 193 stress 0.07374804 
    ## Run 194 stress 0.08382856 
    ## Run 195 stress 0.1799051 
    ## Run 196 stress 0.02693829 
    ## ... Procrustes: rmse 3.952415e-05  max resid 6.804804e-05 
    ## ... Similar to previous best
    ## Run 197 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001633391  max resid 0.0002790365 
    ## ... Similar to previous best
    ## Run 198 stress 0.179905 
    ## Run 199 stress 0.1799051 
    ## Run 200 stress 0.0838283 
    ## Run 201 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001319513  max resid 0.0002249094 
    ## ... Similar to previous best
    ## Run 202 stress 0.075473 
    ## Run 203 stress 0.07374801 
    ## Run 204 stress 0.07374806 
    ## Run 205 stress 0.08382852 
    ## Run 206 stress 0.07374808 
    ## Run 207 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001000721  max resid 0.0001708521 
    ## ... Similar to previous best
    ## Run 208 stress 0.07547301 
    ## Run 209 stress 0.0269383 
    ## ... Procrustes: rmse 7.61909e-05  max resid 0.0001306262 
    ## ... Similar to previous best
    ## Run 210 stress 0.07374803 
    ## Run 211 stress 0.07374805 
    ## Run 212 stress 0.0269383 
    ## ... Procrustes: rmse 7.702756e-05  max resid 0.0001317932 
    ## ... Similar to previous best
    ## Run 213 stress 0.02693829 
    ## ... Procrustes: rmse 5.226856e-05  max resid 7.4058e-05 
    ## ... Similar to previous best
    ## Run 214 stress 0.179905 
    ## Run 215 stress 0.07547299 
    ## Run 216 stress 0.07374804 
    ## Run 217 stress 0.1862383 
    ## Run 218 stress 0.07374801 
    ## Run 219 stress 0.179905 
    ## Run 220 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001078707  max resid 0.0001839145 
    ## ... Similar to previous best
    ## Run 221 stress 0.179905 
    ## Run 222 stress 0.1862384 
    ## Run 223 stress 0.07374809 
    ## Run 224 stress 0.07547299 
    ## Run 225 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001389665  max resid 0.0002378594 
    ## ... Similar to previous best
    ## Run 226 stress 0.07374809 
    ## Run 227 stress 0.07547299 
    ## Run 228 stress 0.1756177 
    ## Run 229 stress 0.1889143 
    ## Run 230 stress 0.08382843 
    ## Run 231 stress 0.08382828 
    ## Run 232 stress 0.02693829 
    ## ... Procrustes: rmse 2.51512e-05  max resid 4.31226e-05 
    ## ... Similar to previous best
    ## Run 233 stress 0.02693834 
    ## ... Procrustes: rmse 9.802613e-05  max resid 0.0001691525 
    ## ... Similar to previous best
    ## Run 234 stress 0.179905 
    ## Run 235 stress 0.07374801 
    ## Run 236 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001716102  max resid 0.0002893821 
    ## ... Similar to previous best
    ## Run 237 stress 0.1756177 
    ## Run 238 stress 0.07374807 
    ## Run 239 stress 0.075473 
    ## Run 240 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001477529  max resid 0.000252233 
    ## ... Similar to previous best
    ## Run 241 stress 0.07547301 
    ## Run 242 stress 0.075473 
    ## Run 243 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001598581  max resid 0.0002619409 
    ## ... Similar to previous best
    ## Run 244 stress 0.1889141 
    ## Run 245 stress 0.179905 
    ## Run 246 stress 0.02693829 
    ## ... Procrustes: rmse 4.655722e-05  max resid 7.908946e-05 
    ## ... Similar to previous best
    ## Run 247 stress 0.07374801 
    ## Run 248 stress 0.07374802 
    ## Run 249 stress 0.07374803 
    ## Run 250 stress 0.02693829 
    ## ... Procrustes: rmse 1.844198e-05  max resid 3.168965e-05 
    ## ... Similar to previous best
    ## Run 251 stress 0.07374803 
    ## Run 252 stress 0.0269383 
    ## ... Procrustes: rmse 6.761757e-05  max resid 0.0001151001 
    ## ... Similar to previous best
    ## Run 253 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001494223  max resid 0.0002549051 
    ## ... Similar to previous best
    ## Run 254 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001679839  max resid 0.0002872302 
    ## ... Similar to previous best
    ## Run 255 stress 0.07374806 
    ## Run 256 stress 0.073748 
    ## Run 257 stress 0.07374804 
    ## Run 258 stress 0.08382843 
    ## Run 259 stress 0.07547299 
    ## Run 260 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001465555  max resid 0.0002511065 
    ## ... Similar to previous best
    ## Run 261 stress 0.07374805 
    ## Run 262 stress 0.07374804 
    ## Run 263 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001576725  max resid 0.0002706841 
    ## ... Similar to previous best
    ## Run 264 stress 0.07374803 
    ## Run 265 stress 0.1862387 
    ## Run 266 stress 0.08382841 
    ## Run 267 stress 0.07374801 
    ## Run 268 stress 0.0269383 
    ## ... Procrustes: rmse 7.265661e-05  max resid 0.0001356751 
    ## ... Similar to previous best
    ## Run 269 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001176583  max resid 0.0001990025 
    ## ... Similar to previous best
    ## Run 270 stress 0.07547302 
    ## Run 271 stress 0.08382836 
    ## Run 272 stress 0.08382854 
    ## Run 273 stress 0.1889143 
    ## Run 274 stress 0.07547301 
    ## Run 275 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001302553  max resid 0.000223211 
    ## ... Similar to previous best
    ## Run 276 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001266568  max resid 0.0002170191 
    ## ... Similar to previous best
    ## Run 277 stress 0.0269383 
    ## ... Procrustes: rmse 7.073319e-05  max resid 0.0001213806 
    ## ... Similar to previous best
    ## Run 278 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001470149  max resid 0.0002508597 
    ## ... Similar to previous best
    ## Run 279 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001688613  max resid 0.0002896267 
    ## ... Similar to previous best
    ## Run 280 stress 0.07374805 
    ## Run 281 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001788443  max resid 0.0003055762 
    ## ... Similar to previous best
    ## Run 282 stress 0.30831 
    ## Run 283 stress 0.07374807 
    ## Run 284 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001079687  max resid 0.0001855318 
    ## ... Similar to previous best
    ## Run 285 stress 0.0269383 
    ## ... Procrustes: rmse 5.985911e-05  max resid 0.0001021027 
    ## ... Similar to previous best
    ## Run 286 stress 0.07547299 
    ## Run 287 stress 0.075473 
    ## Run 288 stress 0.179905 
    ## Run 289 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001380766  max resid 0.0002376399 
    ## ... Similar to previous best
    ## Run 290 stress 0.1889142 
    ## Run 291 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001035934  max resid 0.0001775284 
    ## ... Similar to previous best
    ## Run 292 stress 0.02693829 
    ## ... Procrustes: rmse 3.5253e-05  max resid 6.057884e-05 
    ## ... Similar to previous best
    ## Run 293 stress 0.0269383 
    ## ... Procrustes: rmse 5.753343e-05  max resid 0.00011069 
    ## ... Similar to previous best
    ## Run 294 stress 0.179905 
    ## Run 295 stress 0.08382832 
    ## Run 296 stress 0.075473 
    ## Run 297 stress 0.1862387 
    ## Run 298 stress 0.1756177 
    ## Run 299 stress 0.1862382 
    ## Run 300 stress 0.02693828 
    ## ... Procrustes: rmse 1.266299e-05  max resid 2.170645e-05 
    ## ... Similar to previous best
    ## Run 301 stress 0.1889141 
    ## Run 302 stress 0.1889141 
    ## Run 303 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001561574  max resid 0.000266911 
    ## ... Similar to previous best
    ## Run 304 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001338664  max resid 0.000229234 
    ## ... Similar to previous best
    ## Run 305 stress 0.1799051 
    ## Run 306 stress 0.07374807 
    ## Run 307 stress 0.1756177 
    ## Run 308 stress 0.179905 
    ## Run 309 stress 0.02693839 
    ## ... Procrustes: rmse 0.000177879  max resid 0.0003047699 
    ## ... Similar to previous best
    ## Run 310 stress 0.07374804 
    ## Run 311 stress 0.0269383 
    ## ... Procrustes: rmse 6.714957e-05  max resid 0.0001143963 
    ## ... Similar to previous best
    ## Run 312 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001243957  max resid 0.0002129042 
    ## ... Similar to previous best
    ## Run 313 stress 0.0269383 
    ## ... Procrustes: rmse 6.287877e-05  max resid 0.000108039 
    ## ... Similar to previous best
    ## Run 314 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001178015  max resid 0.0002028423 
    ## ... Similar to previous best
    ## Run 315 stress 0.07374803 
    ## Run 316 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001116331  max resid 0.0001913657 
    ## ... Similar to previous best
    ## Run 317 stress 0.1756177 
    ## Run 318 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001347031  max resid 0.0002300361 
    ## ... Similar to previous best
    ## Run 319 stress 0.08382842 
    ## Run 320 stress 0.07547299 
    ## Run 321 stress 0.02693829 
    ## ... Procrustes: rmse 5.538551e-05  max resid 9.52485e-05 
    ## ... Similar to previous best
    ## Run 322 stress 0.1756177 
    ## Run 323 stress 0.08382829 
    ## Run 324 stress 0.1756177 
    ## Run 325 stress 0.02693833 
    ## ... Procrustes: rmse 0.000108847  max resid 0.0001849259 
    ## ... Similar to previous best
    ## Run 326 stress 0.02693842 
    ## ... Procrustes: rmse 0.0001803264  max resid 0.0003105402 
    ## ... Similar to previous best
    ## Run 327 stress 0.07374802 
    ## Run 328 stress 0.07374803 
    ## Run 329 stress 0.07374801 
    ## Run 330 stress 0.1799051 
    ## Run 331 stress 0.08382848 
    ## Run 332 stress 0.07374801 
    ## Run 333 stress 0.1889141 
    ## Run 334 stress 0.1862386 
    ## Run 335 stress 0.07374815 
    ## Run 336 stress 0.02693829 
    ## ... Procrustes: rmse 4.279213e-05  max resid 7.346202e-05 
    ## ... Similar to previous best
    ## Run 337 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001146821  max resid 0.0001950255 
    ## ... Similar to previous best
    ## Run 338 stress 0.07374804 
    ## Run 339 stress 0.2273965 
    ## Run 340 stress 0.3083098 
    ## Run 341 stress 0.2273964 
    ## Run 342 stress 0.1756177 
    ## Run 343 stress 0.3083098 
    ## Run 344 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001320299  max resid 0.0002268418 
    ## ... Similar to previous best
    ## Run 345 stress 0.07374807 
    ## Run 346 stress 0.0269383 
    ## ... Procrustes: rmse 5.564269e-05  max resid 9.480302e-05 
    ## ... Similar to previous best
    ## Run 347 stress 0.1756177 
    ## Run 348 stress 0.1756177 
    ## Run 349 stress 0.0269383 
    ## ... Procrustes: rmse 5.73515e-05  max resid 9.729746e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.0269383 
    ## ... Procrustes: rmse 5.682535e-05  max resid 9.713447e-05 
    ## ... Similar to previous best
    ## Run 351 stress 0.07374807 
    ## Run 352 stress 0.02693832 
    ## ... Procrustes: rmse 9.441278e-05  max resid 0.0001539592 
    ## ... Similar to previous best
    ## Run 353 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001189993  max resid 0.0002041955 
    ## ... Similar to previous best
    ## Run 354 stress 0.09391929 
    ## Run 355 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001714155  max resid 0.0002927871 
    ## ... Similar to previous best
    ## Run 356 stress 0.02693831 
    ## ... Procrustes: rmse 9.443092e-05  max resid 0.0001641112 
    ## ... Similar to previous best
    ## Run 357 stress 0.07374801 
    ## Run 358 stress 0.0269383 
    ## ... Procrustes: rmse 7.48198e-05  max resid 0.000128067 
    ## ... Similar to previous best
    ## Run 359 stress 0.07374804 
    ## Run 360 stress 0.07374801 
    ## Run 361 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001144407  max resid 0.0001961285 
    ## ... Similar to previous best
    ## Run 362 stress 0.0737481 
    ## Run 363 stress 0.07374805 
    ## Run 364 stress 0.07374803 
    ## Run 365 stress 0.07547299 
    ## Run 366 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001773289  max resid 0.0003038063 
    ## ... Similar to previous best
    ## Run 367 stress 0.07547301 
    ## Run 368 stress 0.3083098 
    ## Run 369 stress 0.07547299 
    ## Run 370 stress 0.075473 
    ## Run 371 stress 0.179905 
    ## Run 372 stress 0.0737481 
    ## Run 373 stress 0.07547299 
    ## Run 374 stress 0.0269383 
    ## ... Procrustes: rmse 5.245773e-05  max resid 9.008936e-05 
    ## ... Similar to previous best
    ## Run 375 stress 0.07374806 
    ## Run 376 stress 0.1889142 
    ## Run 377 stress 0.08382828 
    ## Run 378 stress 0.02693829 
    ## ... Procrustes: rmse 5.160772e-05  max resid 9.129337e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.07374805 
    ## Run 380 stress 0.02693833 
    ## ... Procrustes: rmse 7.399776e-05  max resid 0.0001234131 
    ## ... Similar to previous best
    ## Run 381 stress 0.1889141 
    ## Run 382 stress 0.07547299 
    ## Run 383 stress 0.0737481 
    ## Run 384 stress 0.0269383 
    ## ... Procrustes: rmse 6.649274e-05  max resid 0.0001135731 
    ## ... Similar to previous best
    ## Run 385 stress 0.1799052 
    ## Run 386 stress 0.08382822 
    ## Run 387 stress 0.07374813 
    ## Run 388 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001628494  max resid 0.0002789145 
    ## ... Similar to previous best
    ## Run 389 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001253713  max resid 0.0002149856 
    ## ... Similar to previous best
    ## Run 390 stress 0.09391919 
    ## Run 391 stress 0.1799051 
    ## Run 392 stress 0.07374804 
    ## Run 393 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001458575  max resid 0.0002494746 
    ## ... Similar to previous best
    ## Run 394 stress 0.075473 
    ## Run 395 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001625087  max resid 0.0002798314 
    ## ... Similar to previous best
    ## Run 396 stress 0.1889141 
    ## Run 397 stress 0.075473 
    ## Run 398 stress 0.07374818 
    ## Run 399 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001799884  max resid 0.0003084549 
    ## ... Similar to previous best
    ## Run 400 stress 0.08382825 
    ## Run 401 stress 0.073748 
    ## Run 402 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001378959  max resid 0.0002357399 
    ## ... Similar to previous best
    ## Run 403 stress 0.0269383 
    ## ... Procrustes: rmse 6.23679e-05  max resid 0.0001067407 
    ## ... Similar to previous best
    ## Run 404 stress 0.02693828 
    ## ... Procrustes: rmse 1.546541e-05  max resid 2.656471e-05 
    ## ... Similar to previous best
    ## Run 405 stress 0.02693832 
    ## ... Procrustes: rmse 9.756031e-05  max resid 0.0001674875 
    ## ... Similar to previous best
    ## Run 406 stress 0.07374807 
    ## Run 407 stress 0.07547299 
    ## Run 408 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001460347  max resid 0.0002503063 
    ## ... Similar to previous best
    ## Run 409 stress 0.1756177 
    ## Run 410 stress 0.1756177 
    ## Run 411 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001396617  max resid 0.000239407 
    ## ... Similar to previous best
    ## Run 412 stress 0.0737481 
    ## Run 413 stress 0.179905 
    ## Run 414 stress 0.07374808 
    ## Run 415 stress 0.02693829 
    ## ... Procrustes: rmse 2.0606e-05  max resid 3.288778e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.08382821 
    ## Run 417 stress 0.02693829 
    ## ... Procrustes: rmse 4.024126e-05  max resid 6.895924e-05 
    ## ... Similar to previous best
    ## Run 418 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001332278  max resid 0.0002305408 
    ## ... Similar to previous best
    ## Run 419 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001082406  max resid 0.0001852039 
    ## ... Similar to previous best
    ## Run 420 stress 0.02693832 
    ## ... Procrustes: rmse 9.707342e-05  max resid 0.0001663351 
    ## ... Similar to previous best
    ## Run 421 stress 0.02693829 
    ## ... Procrustes: rmse 3.503406e-05  max resid 5.996808e-05 
    ## ... Similar to previous best
    ## Run 422 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001015631  max resid 0.0001746886 
    ## ... Similar to previous best
    ## Run 423 stress 0.02693839 
    ## ... Procrustes: rmse 0.000177783  max resid 0.0003038563 
    ## ... Similar to previous best
    ## Run 424 stress 0.02693829 
    ## ... Procrustes: rmse 3.502462e-05  max resid 6.045735e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.02693829 
    ## ... Procrustes: rmse 3.130645e-05  max resid 5.279542e-05 
    ## ... Similar to previous best
    ## Run 426 stress 0.07374805 
    ## Run 427 stress 0.07547301 
    ## Run 428 stress 0.1862387 
    ## Run 429 stress 0.07374802 
    ## Run 430 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001113368  max resid 0.0001907727 
    ## ... Similar to previous best
    ## Run 431 stress 0.08382853 
    ## Run 432 stress 0.07547302 
    ## Run 433 stress 0.07374803 
    ## Run 434 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001732206  max resid 0.0002964457 
    ## ... Similar to previous best
    ## Run 435 stress 0.179905 
    ## Run 436 stress 0.07374804 
    ## Run 437 stress 0.07547299 
    ## Run 438 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001312585  max resid 0.0002238264 
    ## ... Similar to previous best
    ## Run 439 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001630691  max resid 0.0002803884 
    ## ... Similar to previous best
    ## Run 440 stress 0.09391917 
    ## Run 441 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001850824  max resid 0.0003161616 
    ## ... Similar to previous best
    ## Run 442 stress 0.07374802 
    ## Run 443 stress 0.07547299 
    ## Run 444 stress 0.09391913 
    ## Run 445 stress 0.07547299 
    ## Run 446 stress 0.07374809 
    ## Run 447 stress 0.07374809 
    ## Run 448 stress 0.02693829 
    ## ... Procrustes: rmse 3.197413e-05  max resid 5.473549e-05 
    ## ... Similar to previous best
    ## Run 449 stress 0.08382821 
    ## Run 450 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001654284  max resid 0.0002825278 
    ## ... Similar to previous best
    ## Run 451 stress 0.07374804 
    ## Run 452 stress 0.02693832 
    ## ... Procrustes: rmse 9.532261e-05  max resid 0.0001639998 
    ## ... Similar to previous best
    ## Run 453 stress 0.07374801 
    ## Run 454 stress 0.0269383 
    ## ... Procrustes: rmse 7.532528e-05  max resid 0.0001289945 
    ## ... Similar to previous best
    ## Run 455 stress 0.0269383 
    ## ... Procrustes: rmse 5.989128e-05  max resid 9.720992e-05 
    ## ... Similar to previous best
    ## Run 456 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001763924  max resid 0.0003023986 
    ## ... Similar to previous best
    ## Run 457 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001345224  max resid 0.0002291267 
    ## ... Similar to previous best
    ## Run 458 stress 0.02693829 
    ## ... Procrustes: rmse 2.091977e-05  max resid 3.500599e-05 
    ## ... Similar to previous best
    ## Run 459 stress 0.07374807 
    ## Run 460 stress 0.07547299 
    ## Run 461 stress 0.07374804 
    ## Run 462 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001559307  max resid 0.000264609 
    ## ... Similar to previous best
    ## Run 463 stress 0.02693837 
    ## ... Procrustes: rmse 0.000155923  max resid 0.0002671176 
    ## ... Similar to previous best
    ## Run 464 stress 0.02693837 
    ## ... Procrustes: rmse 9.557279e-05  max resid 0.000164971 
    ## ... Similar to previous best
    ## Run 465 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001428796  max resid 0.0002439214 
    ## ... Similar to previous best
    ## Run 466 stress 0.07374804 
    ## Run 467 stress 0.07374808 
    ## Run 468 stress 0.07374804 
    ## Run 469 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001278798  max resid 0.0002191723 
    ## ... Similar to previous best
    ## Run 470 stress 0.07374801 
    ## Run 471 stress 0.07374801 
    ## Run 472 stress 0.02693831 
    ## ... Procrustes: rmse 7.914678e-05  max resid 0.0001358047 
    ## ... Similar to previous best
    ## Run 473 stress 0.07547302 
    ## Run 474 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001828221  max resid 0.0003128159 
    ## ... Similar to previous best
    ## Run 475 stress 0.02693832 
    ## ... Procrustes: rmse 9.589619e-05  max resid 0.0001605422 
    ## ... Similar to previous best
    ## Run 476 stress 0.02693831 
    ## ... Procrustes: rmse 8.09151e-05  max resid 0.0001400342 
    ## ... Similar to previous best
    ## Run 477 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001600383  max resid 0.0002727386 
    ## ... Similar to previous best
    ## Run 478 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001370777  max resid 0.0002343283 
    ## ... Similar to previous best
    ## Run 479 stress 0.07547299 
    ## Run 480 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001276238  max resid 0.000218904 
    ## ... Similar to previous best
    ## Run 481 stress 0.075473 
    ## Run 482 stress 0.02693829 
    ## ... Procrustes: rmse 2.398027e-05  max resid 4.102838e-05 
    ## ... Similar to previous best
    ## Run 483 stress 0.1889141 
    ## Run 484 stress 0.07374801 
    ## Run 485 stress 0.0269383 
    ## ... Procrustes: rmse 7.647408e-05  max resid 0.0001311408 
    ## ... Similar to previous best
    ## Run 486 stress 0.07374803 
    ## Run 487 stress 0.07374802 
    ## Run 488 stress 0.09391897 
    ## Run 489 stress 0.0269383 
    ## ... Procrustes: rmse 7.048052e-05  max resid 0.0001217771 
    ## ... Similar to previous best
    ## Run 490 stress 0.07374804 
    ## Run 491 stress 0.07374808 
    ## Run 492 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001033763  max resid 0.0001761471 
    ## ... Similar to previous best
    ## Run 493 stress 0.075473 
    ## Run 494 stress 0.07374803 
    ## Run 495 stress 0.07374805 
    ## Run 496 stress 0.179905 
    ## Run 497 stress 0.09391927 
    ## Run 498 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001076465  max resid 0.0001855077 
    ## ... Similar to previous best
    ## Run 499 stress 0.179905 
    ## Run 500 stress 0.07374809 
    ## *** Best solution repeated 159 times

``` r
### Environmental
# Surveyed sites 
PD_beta_env_NMDS <- metaMDS(PD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06330868 
    ## Run 1 stress 0.08315618 
    ## Run 2 stress 0.06463315 
    ## Run 3 stress 0.06330869 
    ## ... Procrustes: rmse 2.719912e-05  max resid 5.979736e-05 
    ## ... Similar to previous best
    ## Run 4 stress 0.08089037 
    ## Run 5 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 8.437874e-06  max resid 1.514469e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 1.674317e-05  max resid 3.625804e-05 
    ## ... Similar to previous best
    ## Run 7 stress 0.06463313 
    ## Run 8 stress 0.06463313 
    ## Run 9 stress 0.06463315 
    ## Run 10 stress 0.0640262 
    ## Run 11 stress 0.06463313 
    ## Run 12 stress 0.06463321 
    ## Run 13 stress 0.09232133 
    ## Run 14 stress 0.08089036 
    ## Run 15 stress 0.06463316 
    ## Run 16 stress 0.09563687 
    ## Run 17 stress 0.0640263 
    ## Run 18 stress 0.09106964 
    ## Run 19 stress 0.06330868 
    ## ... Procrustes: rmse 1.518146e-05  max resid 4.086173e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.06402633 
    ## Run 21 stress 0.08246336 
    ## Run 22 stress 0.09336664 
    ## Run 23 stress 0.09571429 
    ## Run 24 stress 0.0646332 
    ## Run 25 stress 0.0933668 
    ## Run 26 stress 0.06402624 
    ## Run 27 stress 0.08315609 
    ## Run 28 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 9.554149e-06  max resid 1.720351e-05 
    ## ... Similar to previous best
    ## Run 29 stress 0.06402639 
    ## Run 30 stress 0.08089066 
    ## Run 31 stress 0.09233441 
    ## Run 32 stress 0.08344983 
    ## Run 33 stress 0.06402621 
    ## Run 34 stress 0.06330869 
    ## ... Procrustes: rmse 3.147539e-05  max resid 7.042305e-05 
    ## ... Similar to previous best
    ## Run 35 stress 0.06463315 
    ## Run 36 stress 0.06463313 
    ## Run 37 stress 0.06402616 
    ## Run 38 stress 0.06463314 
    ## Run 39 stress 0.06330876 
    ## ... Procrustes: rmse 0.000107527  max resid 0.0002420748 
    ## ... Similar to previous best
    ## Run 40 stress 0.06402631 
    ## Run 41 stress 0.06330871 
    ## ... Procrustes: rmse 5.430866e-05  max resid 0.0001182796 
    ## ... Similar to previous best
    ## Run 42 stress 0.06402627 
    ## Run 43 stress 0.06463314 
    ## Run 44 stress 0.06402618 
    ## Run 45 stress 0.06402634 
    ## Run 46 stress 0.08315611 
    ## Run 47 stress 0.0633087 
    ## ... Procrustes: rmse 5.559104e-05  max resid 0.000122917 
    ## ... Similar to previous best
    ## Run 48 stress 0.06463315 
    ## Run 49 stress 0.06330871 
    ## ... Procrustes: rmse 1.634726e-05  max resid 5.619403e-05 
    ## ... Similar to previous best
    ## Run 50 stress 0.08315615 
    ## Run 51 stress 0.08357051 
    ## Run 52 stress 0.06463316 
    ## Run 53 stress 0.06463315 
    ## Run 54 stress 0.06330868 
    ## ... Procrustes: rmse 1.462622e-05  max resid 5.450502e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.06402617 
    ## Run 56 stress 0.06402616 
    ## Run 57 stress 0.06463314 
    ## Run 58 stress 0.06330868 
    ## ... Procrustes: rmse 1.008579e-05  max resid 2.97635e-05 
    ## ... Similar to previous best
    ## Run 59 stress 0.06402621 
    ## Run 60 stress 0.06463313 
    ## Run 61 stress 0.08089042 
    ## Run 62 stress 0.06330869 
    ## ... Procrustes: rmse 4.719666e-05  max resid 0.0001051772 
    ## ... Similar to previous best
    ## Run 63 stress 0.06330868 
    ## ... Procrustes: rmse 1.476508e-05  max resid 4.659415e-05 
    ## ... Similar to previous best
    ## Run 64 stress 0.06463313 
    ## Run 65 stress 0.06463313 
    ## Run 66 stress 0.09788382 
    ## Run 67 stress 0.06463313 
    ## Run 68 stress 0.09429113 
    ## Run 69 stress 0.06330869 
    ## ... Procrustes: rmse 1.482158e-05  max resid 4.165537e-05 
    ## ... Similar to previous best
    ## Run 70 stress 0.06330868 
    ## ... Procrustes: rmse 4.401236e-06  max resid 9.83783e-06 
    ## ... Similar to previous best
    ## Run 71 stress 0.0633087 
    ## ... Procrustes: rmse 4.415462e-05  max resid 0.0001551888 
    ## ... Similar to previous best
    ## Run 72 stress 0.3666548 
    ## Run 73 stress 0.06463313 
    ## Run 74 stress 0.06330869 
    ## ... Procrustes: rmse 2.81368e-05  max resid 8.590752e-05 
    ## ... Similar to previous best
    ## Run 75 stress 0.06463315 
    ## Run 76 stress 0.09270094 
    ## Run 77 stress 0.06463314 
    ## Run 78 stress 0.06330871 
    ## ... Procrustes: rmse 6.261452e-05  max resid 0.0001088055 
    ## ... Similar to previous best
    ## Run 79 stress 0.08389992 
    ## Run 80 stress 0.06402617 
    ## Run 81 stress 0.06330875 
    ## ... Procrustes: rmse 0.0001068611  max resid 0.0002386947 
    ## ... Similar to previous best
    ## Run 82 stress 0.06463314 
    ## Run 83 stress 0.06330868 
    ## ... Procrustes: rmse 1.657153e-05  max resid 3.576797e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.06463313 
    ## Run 85 stress 0.06463313 
    ## Run 86 stress 0.06463314 
    ## Run 87 stress 0.09551 
    ## Run 88 stress 0.06330868 
    ## ... Procrustes: rmse 6.121547e-06  max resid 1.865153e-05 
    ## ... Similar to previous best
    ## Run 89 stress 0.06402619 
    ## Run 90 stress 0.0633087 
    ## ... Procrustes: rmse 5.494714e-05  max resid 9.274546e-05 
    ## ... Similar to previous best
    ## Run 91 stress 0.06330873 
    ## ... Procrustes: rmse 9.821422e-05  max resid 0.0002200728 
    ## ... Similar to previous best
    ## Run 92 stress 0.06402633 
    ## Run 93 stress 0.06330869 
    ## ... Procrustes: rmse 2.829714e-05  max resid 0.0001025547 
    ## ... Similar to previous best
    ## Run 94 stress 0.06330872 
    ## ... Procrustes: rmse 8.589422e-05  max resid 0.000191533 
    ## ... Similar to previous best
    ## Run 95 stress 0.08233737 
    ## Run 96 stress 0.06463322 
    ## Run 97 stress 0.06463313 
    ## Run 98 stress 0.06402627 
    ## Run 99 stress 0.08246323 
    ## Run 100 stress 0.06463315 
    ## Run 101 stress 0.06330874 
    ## ... Procrustes: rmse 0.000104172  max resid 0.0002337451 
    ## ... Similar to previous best
    ## Run 102 stress 0.06330878 
    ## ... Procrustes: rmse 0.0001322907  max resid 0.0002951178 
    ## ... Similar to previous best
    ## Run 103 stress 0.0960229 
    ## Run 104 stress 0.08233735 
    ## Run 105 stress 0.09581748 
    ## Run 106 stress 0.06402632 
    ## Run 107 stress 0.0633087 
    ## ... Procrustes: rmse 5.722117e-05  max resid 0.0001270044 
    ## ... Similar to previous best
    ## Run 108 stress 0.08089033 
    ## Run 109 stress 0.09504617 
    ## Run 110 stress 0.06330868 
    ## ... Procrustes: rmse 1.669246e-05  max resid 3.247245e-05 
    ## ... Similar to previous best
    ## Run 111 stress 0.06463318 
    ## Run 112 stress 0.06463314 
    ## Run 113 stress 0.3576493 
    ## Run 114 stress 0.0976747 
    ## Run 115 stress 0.09283188 
    ## Run 116 stress 0.06463313 
    ## Run 117 stress 0.06402627 
    ## Run 118 stress 0.06402628 
    ## Run 119 stress 0.09245977 
    ## Run 120 stress 0.06463313 
    ## Run 121 stress 0.08126163 
    ## Run 122 stress 0.06402628 
    ## Run 123 stress 0.06463313 
    ## Run 124 stress 0.06402616 
    ## Run 125 stress 0.06402626 
    ## Run 126 stress 0.06330872 
    ## ... Procrustes: rmse 6.442017e-05  max resid 0.0001154832 
    ## ... Similar to previous best
    ## Run 127 stress 0.08089066 
    ## Run 128 stress 0.06402622 
    ## Run 129 stress 0.06463313 
    ## Run 130 stress 0.06402616 
    ## Run 131 stress 0.06402626 
    ## Run 132 stress 0.06402617 
    ## Run 133 stress 0.09283209 
    ## Run 134 stress 0.06402616 
    ## Run 135 stress 0.06463313 
    ## Run 136 stress 0.06463315 
    ## Run 137 stress 0.09206658 
    ## Run 138 stress 0.06402621 
    ## Run 139 stress 0.0633087 
    ## ... Procrustes: rmse 6.2111e-05  max resid 0.0001400531 
    ## ... Similar to previous best
    ## Run 140 stress 0.06330869 
    ## ... Procrustes: rmse 4.417229e-05  max resid 9.633572e-05 
    ## ... Similar to previous best
    ## Run 141 stress 0.06402616 
    ## Run 142 stress 0.06463313 
    ## Run 143 stress 0.06463315 
    ## Run 144 stress 0.08089044 
    ## Run 145 stress 0.08089023 
    ## Run 146 stress 0.06463313 
    ## Run 147 stress 0.06463316 
    ## Run 148 stress 0.06463313 
    ## Run 149 stress 0.06402623 
    ## Run 150 stress 0.06402619 
    ## Run 151 stress 0.09270065 
    ## Run 152 stress 0.09290985 
    ## Run 153 stress 0.06330868 
    ## ... Procrustes: rmse 9.472455e-06  max resid 2.381924e-05 
    ## ... Similar to previous best
    ## Run 154 stress 0.08344971 
    ## Run 155 stress 0.09885892 
    ## Run 156 stress 0.06463314 
    ## Run 157 stress 0.09550928 
    ## Run 158 stress 0.06330871 
    ## ... Procrustes: rmse 6.321148e-05  max resid 0.0001414296 
    ## ... Similar to previous best
    ## Run 159 stress 0.08315605 
    ## Run 160 stress 0.06463315 
    ## Run 161 stress 0.06402625 
    ## Run 162 stress 0.09615034 
    ## Run 163 stress 0.08233739 
    ## Run 164 stress 0.06463313 
    ## Run 165 stress 0.08246357 
    ## Run 166 stress 0.06330868 
    ## ... Procrustes: rmse 4.439679e-06  max resid 8.170834e-06 
    ## ... Similar to previous best
    ## Run 167 stress 0.09473003 
    ## Run 168 stress 0.06463318 
    ## Run 169 stress 0.06402618 
    ## Run 170 stress 0.06402631 
    ## Run 171 stress 0.06463313 
    ## Run 172 stress 0.06402624 
    ## Run 173 stress 0.06402619 
    ## Run 174 stress 0.06402616 
    ## Run 175 stress 0.06330873 
    ## ... Procrustes: rmse 7.284611e-05  max resid 0.0002652463 
    ## ... Similar to previous best
    ## Run 176 stress 0.06463316 
    ## Run 177 stress 0.08089007 
    ## Run 178 stress 0.0640262 
    ## Run 179 stress 0.0640263 
    ## Run 180 stress 0.06402636 
    ## Run 181 stress 0.06463313 
    ## Run 182 stress 0.06402619 
    ## Run 183 stress 0.06402641 
    ## Run 184 stress 0.06330868 
    ## ... Procrustes: rmse 2.421588e-05  max resid 8.268896e-05 
    ## ... Similar to previous best
    ## Run 185 stress 0.06463313 
    ## Run 186 stress 0.06463317 
    ## Run 187 stress 0.06330872 
    ## ... Procrustes: rmse 7.236131e-05  max resid 0.0001598396 
    ## ... Similar to previous best
    ## Run 188 stress 0.06330868 
    ## ... Procrustes: rmse 1.362002e-05  max resid 2.847662e-05 
    ## ... Similar to previous best
    ## Run 189 stress 0.06463313 
    ## Run 190 stress 0.06330868 
    ## ... Procrustes: rmse 1.448658e-05  max resid 3.757819e-05 
    ## ... Similar to previous best
    ## Run 191 stress 0.06330868 
    ## ... Procrustes: rmse 9.782853e-06  max resid 3.554255e-05 
    ## ... Similar to previous best
    ## Run 192 stress 0.06330871 
    ## ... Procrustes: rmse 7.713054e-05  max resid 0.0001732971 
    ## ... Similar to previous best
    ## Run 193 stress 0.0835705 
    ## Run 194 stress 0.0646332 
    ## Run 195 stress 0.08089059 
    ## Run 196 stress 0.06463313 
    ## Run 197 stress 0.06330876 
    ## ... Procrustes: rmse 0.0001184262  max resid 0.0002640794 
    ## ... Similar to previous best
    ## Run 198 stress 0.06463317 
    ## Run 199 stress 0.06402619 
    ## Run 200 stress 0.0831561 
    ## Run 201 stress 0.06463316 
    ## Run 202 stress 0.06463314 
    ## Run 203 stress 0.06402623 
    ## Run 204 stress 0.06402627 
    ## Run 205 stress 0.06402634 
    ## Run 206 stress 0.097675 
    ## Run 207 stress 0.06402626 
    ## Run 208 stress 0.0633087 
    ## ... Procrustes: rmse 4.745356e-05  max resid 0.0001095908 
    ## ... Similar to previous best
    ## Run 209 stress 0.06330872 
    ## ... Procrustes: rmse 7.32735e-05  max resid 0.0001620313 
    ## ... Similar to previous best
    ## Run 210 stress 0.06463316 
    ## Run 211 stress 0.06402626 
    ## Run 212 stress 0.06330868 
    ## ... Procrustes: rmse 9.192362e-06  max resid 2.02842e-05 
    ## ... Similar to previous best
    ## Run 213 stress 0.06330868 
    ## ... Procrustes: rmse 1.408317e-05  max resid 3.25563e-05 
    ## ... Similar to previous best
    ## Run 214 stress 0.06330874 
    ## ... Procrustes: rmse 9.753119e-05  max resid 0.0002149972 
    ## ... Similar to previous best
    ## Run 215 stress 0.0633087 
    ## ... Procrustes: rmse 4.583838e-05  max resid 9.778714e-05 
    ## ... Similar to previous best
    ## Run 216 stress 0.06330871 
    ## ... Procrustes: rmse 6.633474e-05  max resid 0.000145717 
    ## ... Similar to previous best
    ## Run 217 stress 0.0823373 
    ## Run 218 stress 0.06463313 
    ## Run 219 stress 0.06402622 
    ## Run 220 stress 0.06330872 
    ## ... Procrustes: rmse 8.863075e-05  max resid 0.0001995493 
    ## ... Similar to previous best
    ## Run 221 stress 0.06402616 
    ## Run 222 stress 0.06330868 
    ## ... Procrustes: rmse 9.699223e-06  max resid 1.981756e-05 
    ## ... Similar to previous best
    ## Run 223 stress 0.06463313 
    ## Run 224 stress 0.06463315 
    ## Run 225 stress 0.08327913 
    ## Run 226 stress 0.06330868 
    ## ... Procrustes: rmse 5.989977e-06  max resid 1.616994e-05 
    ## ... Similar to previous best
    ## Run 227 stress 0.06330868 
    ## ... Procrustes: rmse 4.67988e-06  max resid 8.620281e-06 
    ## ... Similar to previous best
    ## Run 228 stress 0.06330869 
    ## ... Procrustes: rmse 2.767126e-05  max resid 9.931676e-05 
    ## ... Similar to previous best
    ## Run 229 stress 0.06402638 
    ## Run 230 stress 0.06330868 
    ## ... Procrustes: rmse 5.236585e-06  max resid 1.207108e-05 
    ## ... Similar to previous best
    ## Run 231 stress 0.3549387 
    ## Run 232 stress 0.06402619 
    ## Run 233 stress 0.06402618 
    ## Run 234 stress 0.06463313 
    ## Run 235 stress 0.06330869 
    ## ... Procrustes: rmse 2.739847e-05  max resid 9.996805e-05 
    ## ... Similar to previous best
    ## Run 236 stress 0.06463313 
    ## Run 237 stress 0.06463314 
    ## Run 238 stress 0.06330869 
    ## ... Procrustes: rmse 4.113786e-05  max resid 9.339791e-05 
    ## ... Similar to previous best
    ## Run 239 stress 0.06402622 
    ## Run 240 stress 0.06330868 
    ## ... Procrustes: rmse 6.24532e-06  max resid 2.049366e-05 
    ## ... Similar to previous best
    ## Run 241 stress 0.06330869 
    ## ... Procrustes: rmse 2.262499e-05  max resid 7.219669e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.08344968 
    ## Run 243 stress 0.06330869 
    ## ... Procrustes: rmse 3.813063e-05  max resid 8.543209e-05 
    ## ... Similar to previous best
    ## Run 244 stress 0.06463314 
    ## Run 245 stress 0.06463313 
    ## Run 246 stress 0.08233728 
    ## Run 247 stress 0.06330869 
    ## ... Procrustes: rmse 3.272027e-05  max resid 0.0001184789 
    ## ... Similar to previous best
    ## Run 248 stress 0.06463313 
    ## Run 249 stress 0.06463316 
    ## Run 250 stress 0.06463316 
    ## Run 251 stress 0.0633087 
    ## ... Procrustes: rmse 4.665915e-05  max resid 0.0001015793 
    ## ... Similar to previous best
    ## Run 252 stress 0.09463828 
    ## Run 253 stress 0.06330871 
    ## ... Procrustes: rmse 5.85807e-05  max resid 0.0002135335 
    ## ... Similar to previous best
    ## Run 254 stress 0.06463316 
    ## Run 255 stress 0.06330878 
    ## ... Procrustes: rmse 0.0001273812  max resid 0.0002786801 
    ## ... Similar to previous best
    ## Run 256 stress 0.06330871 
    ## ... Procrustes: rmse 6.216222e-05  max resid 0.0001372596 
    ## ... Similar to previous best
    ## Run 257 stress 0.09485781 
    ## Run 258 stress 0.06330869 
    ## ... Procrustes: rmse 1.391221e-05  max resid 3.082689e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.06330871 
    ## ... Procrustes: rmse 6.736657e-05  max resid 0.0001483907 
    ## ... Similar to previous best
    ## Run 260 stress 0.0640262 
    ## Run 261 stress 0.09291007 
    ## Run 262 stress 0.06402625 
    ## Run 263 stress 0.0640263 
    ## Run 264 stress 0.06330868 
    ## ... Procrustes: rmse 9.865756e-06  max resid 3.172041e-05 
    ## ... Similar to previous best
    ## Run 265 stress 0.06330869 
    ## ... Procrustes: rmse 3.509443e-05  max resid 6.228186e-05 
    ## ... Similar to previous best
    ## Run 266 stress 0.06330872 
    ## ... Procrustes: rmse 8.391301e-05  max resid 0.0001881307 
    ## ... Similar to previous best
    ## Run 267 stress 0.06402633 
    ## Run 268 stress 0.09283193 
    ## Run 269 stress 0.06402621 
    ## Run 270 stress 0.06330868 
    ## ... Procrustes: rmse 6.289339e-06  max resid 1.323719e-05 
    ## ... Similar to previous best
    ## Run 271 stress 0.06402618 
    ## Run 272 stress 0.09245973 
    ## Run 273 stress 0.06402623 
    ## Run 274 stress 0.06402617 
    ## Run 275 stress 0.06330871 
    ## ... Procrustes: rmse 6.61705e-05  max resid 0.0001400349 
    ## ... Similar to previous best
    ## Run 276 stress 0.08126184 
    ## Run 277 stress 0.06402618 
    ## Run 278 stress 0.0633087 
    ## ... Procrustes: rmse 5.032481e-05  max resid 0.0001109605 
    ## ... Similar to previous best
    ## Run 279 stress 0.06463315 
    ## Run 280 stress 0.08327917 
    ## Run 281 stress 0.0835706 
    ## Run 282 stress 0.06330868 
    ## ... Procrustes: rmse 4.995998e-06  max resid 1.242924e-05 
    ## ... Similar to previous best
    ## Run 283 stress 0.06330868 
    ## ... Procrustes: rmse 9.074955e-06  max resid 1.99442e-05 
    ## ... Similar to previous best
    ## Run 284 stress 0.06330869 
    ## ... Procrustes: rmse 2.394703e-05  max resid 5.269816e-05 
    ## ... Similar to previous best
    ## Run 285 stress 0.06330877 
    ## ... Procrustes: rmse 0.0001133751  max resid 0.0002457361 
    ## ... Similar to previous best
    ## Run 286 stress 0.06402616 
    ## Run 287 stress 0.06402628 
    ## Run 288 stress 0.08089049 
    ## Run 289 stress 0.06402624 
    ## Run 290 stress 0.06330869 
    ## ... Procrustes: rmse 3.5995e-05  max resid 0.0001251056 
    ## ... Similar to previous best
    ## Run 291 stress 0.06463315 
    ## Run 292 stress 0.06402623 
    ## Run 293 stress 0.06402624 
    ## Run 294 stress 0.06463314 
    ## Run 295 stress 0.06402636 
    ## Run 296 stress 0.06330869 
    ## ... Procrustes: rmse 3.075084e-05  max resid 0.0001068133 
    ## ... Similar to previous best
    ## Run 297 stress 0.08389989 
    ## Run 298 stress 0.08605892 
    ## Run 299 stress 0.06402624 
    ## Run 300 stress 0.06463314 
    ## Run 301 stress 0.09378696 
    ## Run 302 stress 0.06330868 
    ## ... Procrustes: rmse 3.411965e-06  max resid 7.750705e-06 
    ## ... Similar to previous best
    ## Run 303 stress 0.06330869 
    ## ... Procrustes: rmse 2.517214e-05  max resid 6.870507e-05 
    ## ... Similar to previous best
    ## Run 304 stress 0.06463314 
    ## Run 305 stress 0.09157928 
    ## Run 306 stress 0.06330868 
    ## ... Procrustes: rmse 4.681067e-06  max resid 1.177572e-05 
    ## ... Similar to previous best
    ## Run 307 stress 0.0640263 
    ## Run 308 stress 0.08089042 
    ## Run 309 stress 0.06402618 
    ## Run 310 stress 0.08327879 
    ## Run 311 stress 0.06330872 
    ## ... Procrustes: rmse 3.609175e-05  max resid 7.279436e-05 
    ## ... Similar to previous best
    ## Run 312 stress 0.06402641 
    ## Run 313 stress 0.0910696 
    ## Run 314 stress 0.0633087 
    ## ... Procrustes: rmse 5.170338e-05  max resid 0.000114442 
    ## ... Similar to previous best
    ## Run 315 stress 0.06402617 
    ## Run 316 stress 0.06463313 
    ## Run 317 stress 0.0640263 
    ## Run 318 stress 0.09106924 
    ## Run 319 stress 0.06330868 
    ## ... Procrustes: rmse 1.602911e-05  max resid 5.926294e-05 
    ## ... Similar to previous best
    ## Run 320 stress 0.06330868 
    ## ... Procrustes: rmse 5.500925e-06  max resid 1.267975e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.06330869 
    ## ... Procrustes: rmse 2.543024e-05  max resid 8.832267e-05 
    ## ... Similar to previous best
    ## Run 322 stress 0.06463314 
    ## Run 323 stress 0.06402627 
    ## Run 324 stress 0.06463316 
    ## Run 325 stress 0.06463316 
    ## Run 326 stress 0.06330869 
    ## ... Procrustes: rmse 3.949673e-05  max resid 8.878845e-05 
    ## ... Similar to previous best
    ## Run 327 stress 0.06402629 
    ## Run 328 stress 0.06463314 
    ## Run 329 stress 0.06463313 
    ## Run 330 stress 0.06330868 
    ## ... Procrustes: rmse 2.929168e-06  max resid 9.336211e-06 
    ## ... Similar to previous best
    ## Run 331 stress 0.06463314 
    ## Run 332 stress 0.0633087 
    ## ... Procrustes: rmse 6.16916e-05  max resid 0.0001370651 
    ## ... Similar to previous best
    ## Run 333 stress 0.06402627 
    ## Run 334 stress 0.06330868 
    ## ... Procrustes: rmse 4.356672e-06  max resid 9.690527e-06 
    ## ... Similar to previous best
    ## Run 335 stress 0.06402619 
    ## Run 336 stress 0.06402626 
    ## Run 337 stress 0.09283236 
    ## Run 338 stress 0.06330868 
    ## ... Procrustes: rmse 2.479649e-05  max resid 4.394024e-05 
    ## ... Similar to previous best
    ## Run 339 stress 0.08233748 
    ## Run 340 stress 0.06402617 
    ## Run 341 stress 0.09106937 
    ## Run 342 stress 0.06402622 
    ## Run 343 stress 0.08389993 
    ## Run 344 stress 0.06330868 
    ## ... Procrustes: rmse 2.087837e-05  max resid 4.769681e-05 
    ## ... Similar to previous best
    ## Run 345 stress 0.06330869 
    ## ... Procrustes: rmse 2.645538e-05  max resid 9.436897e-05 
    ## ... Similar to previous best
    ## Run 346 stress 0.06463313 
    ## Run 347 stress 0.06463314 
    ## Run 348 stress 0.06463313 
    ## Run 349 stress 0.0633087 
    ## ... Procrustes: rmse 3.761406e-05  max resid 8.078795e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.06402619 
    ## Run 351 stress 0.06402616 
    ## Run 352 stress 0.06463315 
    ## Run 353 stress 0.06402619 
    ## Run 354 stress 0.06463313 
    ## Run 355 stress 0.06402624 
    ## Run 356 stress 0.08301931 
    ## Run 357 stress 0.08357062 
    ## Run 358 stress 0.0633087 
    ## ... Procrustes: rmse 5.70562e-05  max resid 0.0001242579 
    ## ... Similar to previous best
    ## Run 359 stress 0.06402631 
    ## Run 360 stress 0.06463313 
    ## Run 361 stress 0.06330871 
    ## ... Procrustes: rmse 7.563048e-05  max resid 0.000169974 
    ## ... Similar to previous best
    ## Run 362 stress 0.06463313 
    ## Run 363 stress 0.06463313 
    ## Run 364 stress 0.06463314 
    ## Run 365 stress 0.09290993 
    ## Run 366 stress 0.08089038 
    ## Run 367 stress 0.06330868 
    ## ... Procrustes: rmse 7.290906e-06  max resid 1.992894e-05 
    ## ... Similar to previous best
    ## Run 368 stress 0.06330869 
    ## ... Procrustes: rmse 3.442383e-05  max resid 0.0001064237 
    ## ... Similar to previous best
    ## Run 369 stress 0.06463313 
    ## Run 370 stress 0.06463313 
    ## Run 371 stress 0.08089038 
    ## Run 372 stress 0.0812618 
    ## Run 373 stress 0.06402615 
    ## Run 374 stress 0.06330868 
    ## ... Procrustes: rmse 5.066709e-06  max resid 1.681548e-05 
    ## ... Similar to previous best
    ## Run 375 stress 0.06330874 
    ## ... Procrustes: rmse 0.0001029173  max resid 0.0002327484 
    ## ... Similar to previous best
    ## Run 376 stress 0.06402616 
    ## Run 377 stress 0.06463313 
    ## Run 378 stress 0.08315607 
    ## Run 379 stress 0.08315617 
    ## Run 380 stress 0.06463315 
    ## Run 381 stress 0.09472756 
    ## Run 382 stress 0.06330869 
    ## ... Procrustes: rmse 4.215835e-05  max resid 0.0001525668 
    ## ... Similar to previous best
    ## Run 383 stress 0.09581749 
    ## Run 384 stress 0.06402621 
    ## Run 385 stress 0.06330868 
    ## ... Procrustes: rmse 8.618019e-06  max resid 3.172355e-05 
    ## ... Similar to previous best
    ## Run 386 stress 0.09563703 
    ## Run 387 stress 0.08315615 
    ## Run 388 stress 0.06402616 
    ## Run 389 stress 0.06402638 
    ## Run 390 stress 0.06402634 
    ## Run 391 stress 0.09217179 
    ## Run 392 stress 0.0640262 
    ## Run 393 stress 0.06463313 
    ## Run 394 stress 0.08487995 
    ## Run 395 stress 0.0812618 
    ## Run 396 stress 0.06463314 
    ## Run 397 stress 0.06330876 
    ## ... Procrustes: rmse 0.0001178082  max resid 0.0002623979 
    ## ... Similar to previous best
    ## Run 398 stress 0.08246304 
    ## Run 399 stress 0.06330869 
    ## ... Procrustes: rmse 4.373349e-05  max resid 9.723774e-05 
    ## ... Similar to previous best
    ## Run 400 stress 0.06402622 
    ## Run 401 stress 0.06463316 
    ## Run 402 stress 0.06330877 
    ## ... Procrustes: rmse 3.267974e-05  max resid 8.937387e-05 
    ## ... Similar to previous best
    ## Run 403 stress 0.06330874 
    ## ... Procrustes: rmse 0.0001016065  max resid 0.0002336106 
    ## ... Similar to previous best
    ## Run 404 stress 0.06330868 
    ## ... Procrustes: rmse 1.102712e-05  max resid 2.283602e-05 
    ## ... Similar to previous best
    ## Run 405 stress 0.08126161 
    ## Run 406 stress 0.06463313 
    ## Run 407 stress 0.06463313 
    ## Run 408 stress 0.06463315 
    ## Run 409 stress 0.06463313 
    ## Run 410 stress 0.09547946 
    ## Run 411 stress 0.06330869 
    ## ... Procrustes: rmse 3.05964e-05  max resid 6.930462e-05 
    ## ... Similar to previous best
    ## Run 412 stress 0.0640262 
    ## Run 413 stress 0.08327878 
    ## Run 414 stress 0.06330872 
    ## ... Procrustes: rmse 8.023804e-05  max resid 0.000183097 
    ## ... Similar to previous best
    ## Run 415 stress 0.06402618 
    ## Run 416 stress 0.06402626 
    ## Run 417 stress 0.08126182 
    ## Run 418 stress 0.09106987 
    ## Run 419 stress 0.06330869 
    ## ... Procrustes: rmse 3.744594e-05  max resid 0.0001245505 
    ## ... Similar to previous best
    ## Run 420 stress 0.09106924 
    ## Run 421 stress 0.06463314 
    ## Run 422 stress 0.06330869 
    ## ... Procrustes: rmse 3.020849e-05  max resid 6.737799e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.09106937 
    ## Run 424 stress 0.06463314 
    ## Run 425 stress 0.06402623 
    ## Run 426 stress 0.06402619 
    ## Run 427 stress 0.06463313 
    ## Run 428 stress 0.08089053 
    ## Run 429 stress 0.08089057 
    ## Run 430 stress 0.09846445 
    ## Run 431 stress 0.06463314 
    ## Run 432 stress 0.06463321 
    ## Run 433 stress 0.06330868 
    ## ... Procrustes: rmse 2.246778e-05  max resid 5.029752e-05 
    ## ... Similar to previous best
    ## Run 434 stress 0.06463314 
    ## Run 435 stress 0.06330871 
    ## ... Procrustes: rmse 6.475822e-05  max resid 0.0002369418 
    ## ... Similar to previous best
    ## Run 436 stress 0.06463313 
    ## Run 437 stress 0.09612804 
    ## Run 438 stress 0.09232166 
    ## Run 439 stress 0.08344975 
    ## Run 440 stress 0.06330868 
    ## ... Procrustes: rmse 1.445065e-05  max resid 3.064965e-05 
    ## ... Similar to previous best
    ## Run 441 stress 0.06463314 
    ## Run 442 stress 0.06463313 
    ## Run 443 stress 0.06463313 
    ## Run 444 stress 0.09283241 
    ## Run 445 stress 0.06402624 
    ## Run 446 stress 0.06463315 
    ## Run 447 stress 0.06402625 
    ## Run 448 stress 0.06330871 
    ## ... Procrustes: rmse 6.225002e-05  max resid 0.0001077707 
    ## ... Similar to previous best
    ## Run 449 stress 0.0633087 
    ## ... Procrustes: rmse 2.249081e-05  max resid 4.349023e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.06330869 
    ## ... Procrustes: rmse 3.103537e-05  max resid 6.734973e-05 
    ## ... Similar to previous best
    ## Run 451 stress 0.06463313 
    ## Run 452 stress 0.09245931 
    ## Run 453 stress 0.0640262 
    ## Run 454 stress 0.06402621 
    ## Run 455 stress 0.06330869 
    ## ... Procrustes: rmse 3.969627e-05  max resid 9.030952e-05 
    ## ... Similar to previous best
    ## Run 456 stress 0.09428664 
    ## Run 457 stress 0.06402623 
    ## Run 458 stress 0.06402638 
    ## Run 459 stress 0.06463313 
    ## Run 460 stress 0.0640262 
    ## Run 461 stress 0.06402618 
    ## Run 462 stress 0.06463315 
    ## Run 463 stress 0.06463313 
    ## Run 464 stress 0.06463313 
    ## Run 465 stress 0.08126184 
    ## Run 466 stress 0.06402626 
    ## Run 467 stress 0.06463316 
    ## Run 468 stress 0.0808904 
    ## Run 469 stress 0.06402616 
    ## Run 470 stress 0.06402628 
    ## Run 471 stress 0.06463313 
    ## Run 472 stress 0.08301952 
    ## Run 473 stress 0.06402621 
    ## Run 474 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 3.389637e-06  max resid 8.114043e-06 
    ## ... Similar to previous best
    ## Run 475 stress 0.06402623 
    ## Run 476 stress 0.09245956 
    ## Run 477 stress 0.06463314 
    ## Run 478 stress 0.09106937 
    ## Run 479 stress 0.06402622 
    ## Run 480 stress 0.06330868 
    ## ... Procrustes: rmse 6.662707e-06  max resid 2.387659e-05 
    ## ... Similar to previous best
    ## Run 481 stress 0.06330868 
    ## ... Procrustes: rmse 2.692464e-05  max resid 9.007199e-05 
    ## ... Similar to previous best
    ## Run 482 stress 0.06463313 
    ## Run 483 stress 0.08089016 
    ## Run 484 stress 0.08246373 
    ## Run 485 stress 0.06463316 
    ## Run 486 stress 0.06402618 
    ## Run 487 stress 0.06463315 
    ## Run 488 stress 0.06402639 
    ## Run 489 stress 0.06402622 
    ## Run 490 stress 0.06402621 
    ## Run 491 stress 0.0633087 
    ## ... Procrustes: rmse 5.268771e-05  max resid 0.0001193103 
    ## ... Similar to previous best
    ## Run 492 stress 0.06402622 
    ## Run 493 stress 0.06463315 
    ## Run 494 stress 0.08315604 
    ## Run 495 stress 0.06463315 
    ## Run 496 stress 0.06463319 
    ## Run 497 stress 0.06463316 
    ## Run 498 stress 0.09245964 
    ## Run 499 stress 0.06330876 
    ## ... Procrustes: rmse 0.0001068923  max resid 0.0002325115 
    ## ... Similar to previous best
    ## Run 500 stress 0.06463313 
    ## *** Best solution repeated 5 times

``` r
# Mixed and stratified lakes
PD_beta_env_MS_NMDS <- metaMDS(PD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.05979235 
    ## Run 2 stress 0.05171062 
    ## Run 3 stress 0.06136675 
    ## Run 4 stress 0.05829933 
    ## Run 5 stress 0.05051299 
    ## ... Procrustes: rmse 0.03575913  max resid 0.1220654 
    ## Run 6 stress 0.05171062 
    ## Run 7 stress 0.05036801 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002929221  max resid 0.0007457058 
    ## ... Similar to previous best
    ## Run 8 stress 0.06704952 
    ## Run 9 stress 0.05036795 
    ## ... New best solution
    ## ... Procrustes: rmse 4.000882e-05  max resid 9.156621e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.05051336 
    ## ... Procrustes: rmse 0.03565739  max resid 0.1218382 
    ## Run 11 stress 0.06832183 
    ## Run 12 stress 0.05171063 
    ## Run 13 stress 0.05036797 
    ## ... Procrustes: rmse 2.146403e-05  max resid 5.329471e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.05036796 
    ## ... Procrustes: rmse 0.000211728  max resid 0.0005384954 
    ## ... Similar to previous best
    ## Run 15 stress 0.05171058 
    ## Run 16 stress 0.05051349 
    ## ... Procrustes: rmse 0.03562859  max resid 0.1215697 
    ## Run 17 stress 0.05036792 
    ## ... New best solution
    ## ... Procrustes: rmse 9.456244e-05  max resid 0.0002381264 
    ## ... Similar to previous best
    ## Run 18 stress 0.06647197 
    ## Run 19 stress 0.05337577 
    ## Run 20 stress 0.05403581 
    ## Run 21 stress 0.06640365 
    ## Run 22 stress 0.06630527 
    ## Run 23 stress 0.06027449 
    ## Run 24 stress 0.06362792 
    ## Run 25 stress 0.05829932 
    ## Run 26 stress 0.05051295 
    ## ... Procrustes: rmse 0.03567486  max resid 0.1217796 
    ## Run 27 stress 0.05051339 
    ## ... Procrustes: rmse 0.03574338  max resid 0.1219526 
    ## Run 28 stress 0.05928237 
    ## Run 29 stress 0.06478549 
    ## Run 30 stress 0.06027455 
    ## Run 31 stress 0.05900152 
    ## Run 32 stress 0.05968627 
    ## Run 33 stress 0.06027442 
    ## Run 34 stress 0.050368 
    ## ... Procrustes: rmse 0.0001473546  max resid 0.0003745776 
    ## ... Similar to previous best
    ## Run 35 stress 0.06027426 
    ## Run 36 stress 0.05968634 
    ## Run 37 stress 0.05036807 
    ## ... Procrustes: rmse 0.0001898388  max resid 0.0004795833 
    ## ... Similar to previous best
    ## Run 38 stress 0.05051348 
    ## ... Procrustes: rmse 0.03574139  max resid 0.1221104 
    ## Run 39 stress 0.0596864 
    ## Run 40 stress 0.06362779 
    ## Run 41 stress 0.06478517 
    ## Run 42 stress 0.0582993 
    ## Run 43 stress 0.05337592 
    ## Run 44 stress 0.05051346 
    ## ... Procrustes: rmse 0.03563818  max resid 0.1216219 
    ## Run 45 stress 0.05968629 
    ## Run 46 stress 0.05602047 
    ## Run 47 stress 0.07026165 
    ## Run 48 stress 0.05466706 
    ## Run 49 stress 0.05466705 
    ## Run 50 stress 0.06478527 
    ## Run 51 stress 0.06111778 
    ## Run 52 stress 0.06027428 
    ## Run 53 stress 0.05979232 
    ## Run 54 stress 0.06231303 
    ## Run 55 stress 0.05928244 
    ## Run 56 stress 0.05829929 
    ## Run 57 stress 0.05337592 
    ## Run 58 stress 0.06309836 
    ## Run 59 stress 0.05900148 
    ## Run 60 stress 0.05051294 
    ## ... Procrustes: rmse 0.03568217  max resid 0.1218023 
    ## Run 61 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.654558e-05  max resid 3.351872e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.05968629 
    ## Run 63 stress 0.05829928 
    ## Run 64 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001209151  max resid 0.000327834 
    ## ... Similar to previous best
    ## Run 65 stress 0.05337583 
    ## Run 66 stress 0.05051326 
    ## ... Procrustes: rmse 0.03570643  max resid 0.1218419 
    ## Run 67 stress 0.05968632 
    ## Run 68 stress 0.05051347 
    ## ... Procrustes: rmse 0.03565971  max resid 0.1216848 
    ## Run 69 stress 0.05051344 
    ## ... Procrustes: rmse 0.03574424  max resid 0.1221186 
    ## Run 70 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001305126  max resid 0.0003181678 
    ## ... Similar to previous best
    ## Run 71 stress 0.05051279 
    ## ... Procrustes: rmse 0.03568758  max resid 0.1219095 
    ## Run 72 stress 0.05036796 
    ## ... Procrustes: rmse 9.314329e-05  max resid 0.0002357138 
    ## ... Similar to previous best
    ## Run 73 stress 0.05466708 
    ## Run 74 stress 0.05829929 
    ## Run 75 stress 0.05171059 
    ## Run 76 stress 0.06630527 
    ## Run 77 stress 0.05051317 
    ## ... Procrustes: rmse 0.03569667  max resid 0.1219663 
    ## Run 78 stress 0.05602044 
    ## Run 79 stress 0.06730244 
    ## Run 80 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001050334  max resid 0.0002631168 
    ## ... Similar to previous best
    ## Run 81 stress 0.05466704 
    ## Run 82 stress 0.05051324 
    ## ... Procrustes: rmse 0.03569873  max resid 0.1219761 
    ## Run 83 stress 0.05968638 
    ## Run 84 stress 0.050368 
    ## ... Procrustes: rmse 0.0001545708  max resid 0.0003875208 
    ## ... Similar to previous best
    ## Run 85 stress 0.05036797 
    ## ... Procrustes: rmse 0.000129461  max resid 0.0003273179 
    ## ... Similar to previous best
    ## Run 86 stress 0.05051292 
    ## ... Procrustes: rmse 0.03575835  max resid 0.1221322 
    ## Run 87 stress 0.3657341 
    ## Run 88 stress 0.05829931 
    ## Run 89 stress 0.06904468 
    ## Run 90 stress 0.05928222 
    ## Run 91 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001658707  max resid 0.0004048407 
    ## ... Similar to previous best
    ## Run 92 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001539432  max resid 0.0004120418 
    ## ... Similar to previous best
    ## Run 93 stress 0.06231318 
    ## Run 94 stress 0.05171061 
    ## Run 95 stress 0.05051309 
    ## ... Procrustes: rmse 0.03572195  max resid 0.1220381 
    ## Run 96 stress 0.06111767 
    ## Run 97 stress 0.05979236 
    ## Run 98 stress 0.05051322 
    ## ... Procrustes: rmse 0.03570359  max resid 0.121847 
    ## Run 99 stress 0.05051332 
    ## ... Procrustes: rmse 0.03572014  max resid 0.1220431 
    ## Run 100 stress 0.05051307 
    ## ... Procrustes: rmse 0.035701  max resid 0.1219745 
    ## Run 101 stress 0.0533758 
    ## Run 102 stress 0.06231303 
    ## Run 103 stress 0.05337586 
    ## Run 104 stress 0.05036794 
    ## ... Procrustes: rmse 8.99101e-05  max resid 0.0002190923 
    ## ... Similar to previous best
    ## Run 105 stress 0.05036795 
    ## ... Procrustes: rmse 8.581519e-05  max resid 0.0002202508 
    ## ... Similar to previous best
    ## Run 106 stress 0.05051335 
    ## ... Procrustes: rmse 0.03573665  max resid 0.122093 
    ## Run 107 stress 0.06837525 
    ## Run 108 stress 0.07221116 
    ## Run 109 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001060378  max resid 0.0002716564 
    ## ... Similar to previous best
    ## Run 110 stress 0.05171065 
    ## Run 111 stress 0.06630534 
    ## Run 112 stress 0.05171062 
    ## Run 113 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001345936  max resid 0.0003347515 
    ## ... Similar to previous best
    ## Run 114 stress 0.05171059 
    ## Run 115 stress 0.05051355 
    ## ... Procrustes: rmse 0.03572322  max resid 0.1218667 
    ## Run 116 stress 0.06309864 
    ## Run 117 stress 0.05036793 
    ## ... Procrustes: rmse 4.825516e-05  max resid 9.774475e-05 
    ## ... Similar to previous best
    ## Run 118 stress 0.05403584 
    ## Run 119 stress 0.05171059 
    ## Run 120 stress 0.06309862 
    ## Run 121 stress 0.05403554 
    ## Run 122 stress 0.05337586 
    ## Run 123 stress 0.05051325 
    ## ... Procrustes: rmse 0.03573899  max resid 0.1220966 
    ## Run 124 stress 0.05337578 
    ## Run 125 stress 0.05403654 
    ## Run 126 stress 0.05602053 
    ## Run 127 stress 0.06027429 
    ## Run 128 stress 0.05051314 
    ## ... Procrustes: rmse 0.03568404  max resid 0.1219287 
    ## Run 129 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001243735  max resid 0.0002924475 
    ## ... Similar to previous best
    ## Run 130 stress 0.0613668 
    ## Run 131 stress 0.05051316 
    ## ... Procrustes: rmse 0.03572609  max resid 0.1220541 
    ## Run 132 stress 0.0505133 
    ## ... Procrustes: rmse 0.03572635  max resid 0.1220602 
    ## Run 133 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001125662  max resid 0.0002745063 
    ## ... Similar to previous best
    ## Run 134 stress 0.05051284 
    ## ... Procrustes: rmse 0.03569075  max resid 0.1219257 
    ## Run 135 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001802341  max resid 0.0004513132 
    ## ... Similar to previous best
    ## Run 136 stress 0.05051345 
    ## ... Procrustes: rmse 0.03568036  max resid 0.1217529 
    ## Run 137 stress 0.05051352 
    ## ... Procrustes: rmse 0.0357551  max resid 0.1221531 
    ## Run 138 stress 0.06231312 
    ## Run 139 stress 0.05829935 
    ## Run 140 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001163015  max resid 0.0002847274 
    ## ... Similar to previous best
    ## Run 141 stress 0.05928241 
    ## Run 142 stress 0.05051338 
    ## ... Procrustes: rmse 0.03570751  max resid 0.1220073 
    ## Run 143 stress 0.06478532 
    ## Run 144 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001117053  max resid 0.0003136873 
    ## ... Similar to previous best
    ## Run 145 stress 0.06231299 
    ## Run 146 stress 0.054667 
    ## Run 147 stress 0.05602052 
    ## Run 148 stress 0.05968641 
    ## Run 149 stress 0.05171063 
    ## Run 150 stress 0.05051356 
    ## ... Procrustes: rmse 0.03570377  max resid 0.1219966 
    ## Run 151 stress 0.06950997 
    ## Run 152 stress 0.05051295 
    ## ... Procrustes: rmse 0.03569  max resid 0.1218282 
    ## Run 153 stress 0.05900154 
    ## Run 154 stress 0.05036793 
    ## ... Procrustes: rmse 7.640115e-05  max resid 0.0001838537 
    ## ... Similar to previous best
    ## Run 155 stress 0.05403573 
    ## Run 156 stress 0.05036792 
    ## ... Procrustes: rmse 3.733069e-05  max resid 7.546448e-05 
    ## ... Similar to previous best
    ## Run 157 stress 0.0517106 
    ## Run 158 stress 0.05171065 
    ## Run 159 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001710824  max resid 0.0004139837 
    ## ... Similar to previous best
    ## Run 160 stress 0.05171061 
    ## Run 161 stress 0.05979231 
    ## Run 162 stress 0.06478521 
    ## Run 163 stress 0.05337588 
    ## Run 164 stress 0.05466711 
    ## Run 165 stress 0.05051296 
    ## ... Procrustes: rmse 0.03567347  max resid 0.1218847 
    ## Run 166 stress 0.05171062 
    ## Run 167 stress 0.05968631 
    ## Run 168 stress 0.05979226 
    ## Run 169 stress 0.06826705 
    ## Run 170 stress 0.06038946 
    ## Run 171 stress 0.06231302 
    ## Run 172 stress 0.05036795 
    ## ... Procrustes: rmse 7.252801e-05  max resid 0.0002060508 
    ## ... Similar to previous best
    ## Run 173 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001354269  max resid 0.0003403574 
    ## ... Similar to previous best
    ## Run 174 stress 0.06111765 
    ## Run 175 stress 0.06730223 
    ## Run 176 stress 0.06537627 
    ## Run 177 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001678104  max resid 0.0004474869 
    ## ... Similar to previous best
    ## Run 178 stress 0.05171058 
    ## Run 179 stress 0.06478532 
    ## Run 180 stress 0.05979241 
    ## Run 181 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001650307  max resid 0.0003952936 
    ## ... Similar to previous best
    ## Run 182 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001801024  max resid 0.0004432944 
    ## ... Similar to previous best
    ## Run 183 stress 0.05051327 
    ## ... Procrustes: rmse 0.03566401  max resid 0.1218738 
    ## Run 184 stress 0.06038958 
    ## Run 185 stress 0.05466701 
    ## Run 186 stress 0.0517106 
    ## Run 187 stress 0.050368 
    ## ... Procrustes: rmse 0.0001361482  max resid 0.0003376985 
    ## ... Similar to previous best
    ## Run 188 stress 0.05602056 
    ## Run 189 stress 0.05036793 
    ## ... Procrustes: rmse 6.5314e-05  max resid 0.0001613677 
    ## ... Similar to previous best
    ## Run 190 stress 0.06309864 
    ## Run 191 stress 0.05051296 
    ## ... Procrustes: rmse 0.0355533  max resid 0.1214481 
    ## Run 192 stress 0.06231296 
    ## Run 193 stress 0.05051282 
    ## ... Procrustes: rmse 0.03559323  max resid 0.1216179 
    ## Run 194 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001487621  max resid 0.0003696588 
    ## ... Similar to previous best
    ## Run 195 stress 0.05051305 
    ## ... Procrustes: rmse 0.03568334  max resid 0.1217943 
    ## Run 196 stress 0.06038953 
    ## Run 197 stress 0.0517106 
    ## Run 198 stress 0.05051297 
    ## ... Procrustes: rmse 0.03569151  max resid 0.1218315 
    ## Run 199 stress 0.05171061 
    ## Run 200 stress 0.05829936 
    ## Run 201 stress 0.050368 
    ## ... Procrustes: rmse 0.0001565037  max resid 0.0003915615 
    ## ... Similar to previous best
    ## Run 202 stress 0.05979228 
    ## Run 203 stress 0.06038968 
    ## Run 204 stress 0.0590015 
    ## Run 205 stress 0.06478531 
    ## Run 206 stress 0.05051294 
    ## ... Procrustes: rmse 0.03569227  max resid 0.1219376 
    ## Run 207 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001445334  max resid 0.000350767 
    ## ... Similar to previous best
    ## Run 208 stress 0.05036811 
    ## ... Procrustes: rmse 0.0001914615  max resid 0.0005287151 
    ## ... Similar to previous best
    ## Run 209 stress 0.05337581 
    ## Run 210 stress 0.06362788 
    ## Run 211 stress 0.05036801 
    ## ... Procrustes: rmse 8.063531e-05  max resid 0.0002063135 
    ## ... Similar to previous best
    ## Run 212 stress 0.05036794 
    ## ... Procrustes: rmse 5.981509e-05  max resid 0.0001633281 
    ## ... Similar to previous best
    ## Run 213 stress 0.05829932 
    ## Run 214 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001710156  max resid 0.0004200432 
    ## ... Similar to previous best
    ## Run 215 stress 0.06730232 
    ## Run 216 stress 0.05051284 
    ## ... Procrustes: rmse 0.03578951  max resid 0.1221622 
    ## Run 217 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001344335  max resid 0.0003238502 
    ## ... Similar to previous best
    ## Run 218 stress 0.05171061 
    ## Run 219 stress 0.05036794 
    ## ... Procrustes: rmse 9.352236e-05  max resid 0.0002293614 
    ## ... Similar to previous best
    ## Run 220 stress 0.05602062 
    ## Run 221 stress 0.06630525 
    ## Run 222 stress 0.05036796 
    ## ... Procrustes: rmse 9.518081e-05  max resid 0.0002318939 
    ## ... Similar to previous best
    ## Run 223 stress 0.05602056 
    ## Run 224 stress 0.06362791 
    ## Run 225 stress 0.06111779 
    ## Run 226 stress 0.05968628 
    ## Run 227 stress 0.05403581 
    ## Run 228 stress 0.05171065 
    ## Run 229 stress 0.05968629 
    ## Run 230 stress 0.06309861 
    ## Run 231 stress 0.05403595 
    ## Run 232 stress 0.0582993 
    ## Run 233 stress 0.05051347 
    ## ... Procrustes: rmse 0.03571397  max resid 0.1220297 
    ## Run 234 stress 0.05928225 
    ## Run 235 stress 0.05968638 
    ## Run 236 stress 0.05036798 
    ## ... Procrustes: rmse 7.069183e-05  max resid 0.0001802096 
    ## ... Similar to previous best
    ## Run 237 stress 0.05051307 
    ## ... Procrustes: rmse 0.03569427  max resid 0.1219538 
    ## Run 238 stress 0.05829931 
    ## Run 239 stress 0.05051294 
    ## ... Procrustes: rmse 0.03568507  max resid 0.1218139 
    ## Run 240 stress 0.05051359 
    ## ... Procrustes: rmse 0.03570058  max resid 0.1217974 
    ## Run 241 stress 0.05051311 
    ## ... Procrustes: rmse 0.03573085  max resid 0.1219327 
    ## Run 242 stress 0.05466702 
    ## Run 243 stress 0.05051298 
    ## ... Procrustes: rmse 0.03569374  max resid 0.1218346 
    ## Run 244 stress 0.05602051 
    ## Run 245 stress 0.05036796 
    ## ... Procrustes: rmse 9.725902e-05  max resid 0.0002601944 
    ## ... Similar to previous best
    ## Run 246 stress 0.05171059 
    ## Run 247 stress 0.05036792 
    ## ... Procrustes: rmse 3.931982e-05  max resid 0.0001058892 
    ## ... Similar to previous best
    ## Run 248 stress 0.05036795 
    ## ... Procrustes: rmse 8.954481e-05  max resid 0.000223468 
    ## ... Similar to previous best
    ## Run 249 stress 0.06203063 
    ## Run 250 stress 0.05171065 
    ## Run 251 stress 0.0582993 
    ## Run 252 stress 0.0546671 
    ## Run 253 stress 0.0517107 
    ## Run 254 stress 0.05602051 
    ## Run 255 stress 0.050368 
    ## ... Procrustes: rmse 0.0001245592  max resid 0.0003264647 
    ## ... Similar to previous best
    ## Run 256 stress 0.06136679 
    ## Run 257 stress 0.0592825 
    ## Run 258 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001624549  max resid 0.000408578 
    ## ... Similar to previous best
    ## Run 259 stress 0.05036794 
    ## ... Procrustes: rmse 9.535474e-05  max resid 0.000223809 
    ## ... Similar to previous best
    ## Run 260 stress 0.06111765 
    ## Run 261 stress 0.06027433 
    ## Run 262 stress 0.050513 
    ## ... Procrustes: rmse 0.03573292  max resid 0.1219511 
    ## Run 263 stress 0.05403555 
    ## Run 264 stress 0.0613668 
    ## Run 265 stress 0.06938764 
    ## Run 266 stress 0.05829934 
    ## Run 267 stress 0.05466704 
    ## Run 268 stress 0.06640331 
    ## Run 269 stress 0.05051296 
    ## ... Procrustes: rmse 0.0356837  max resid 0.1218068 
    ## Run 270 stress 0.06938769 
    ## Run 271 stress 0.0517106 
    ## Run 272 stress 0.06478516 
    ## Run 273 stress 0.06136676 
    ## Run 274 stress 0.0505133 
    ## ... Procrustes: rmse 0.03566502  max resid 0.1217196 
    ## Run 275 stress 0.05036794 
    ## ... Procrustes: rmse 5.667322e-05  max resid 0.0001400208 
    ## ... Similar to previous best
    ## Run 276 stress 0.05051297 
    ## ... Procrustes: rmse 0.03577027  max resid 0.1221715 
    ## Run 277 stress 0.06038955 
    ## Run 278 stress 0.06730212 
    ## Run 279 stress 0.05466702 
    ## Run 280 stress 0.05051328 
    ## ... Procrustes: rmse 0.03570127  max resid 0.1219847 
    ## Run 281 stress 0.05051295 
    ## ... Procrustes: rmse 0.03571359  max resid 0.1219014 
    ## Run 282 stress 0.06647169 
    ## Run 283 stress 0.07204087 
    ## Run 284 stress 0.050368 
    ## ... Procrustes: rmse 0.0001598441  max resid 0.0004015886 
    ## ... Similar to previous best
    ## Run 285 stress 0.05051322 
    ## ... Procrustes: rmse 0.03573142  max resid 0.121921 
    ## Run 286 stress 0.05036792 
    ## ... Procrustes: rmse 5.774692e-05  max resid 0.000136831 
    ## ... Similar to previous best
    ## Run 287 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001038459  max resid 0.0002554638 
    ## ... Similar to previous best
    ## Run 288 stress 0.0664034 
    ## Run 289 stress 0.05036792 
    ## ... Procrustes: rmse 6.27422e-05  max resid 0.0001513671 
    ## ... Similar to previous best
    ## Run 290 stress 0.05051322 
    ## ... Procrustes: rmse 0.03569682  max resid 0.1219695 
    ## Run 291 stress 0.05051332 
    ## ... Procrustes: rmse 0.03569974  max resid 0.1219802 
    ## Run 292 stress 0.05171058 
    ## Run 293 stress 0.05466705 
    ## Run 294 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001186856  max resid 0.0003146213 
    ## ... Similar to previous best
    ## Run 295 stress 0.05036796 
    ## ... Procrustes: rmse 9.415626e-05  max resid 0.0002451699 
    ## ... Similar to previous best
    ## Run 296 stress 0.06038947 
    ## Run 297 stress 0.05051339 
    ## ... Procrustes: rmse 0.0357383  max resid 0.1220994 
    ## Run 298 stress 0.05466701 
    ## Run 299 stress 0.05051278 
    ## ... Procrustes: rmse 0.03564313  max resid 0.1217715 
    ## Run 300 stress 0.05466709 
    ## Run 301 stress 0.06478552 
    ## Run 302 stress 0.06136683 
    ## Run 303 stress 0.05036811 
    ## ... Procrustes: rmse 0.0002099662  max resid 0.0005698243 
    ## ... Similar to previous best
    ## Run 304 stress 0.06478531 
    ## Run 305 stress 0.05466706 
    ## Run 306 stress 0.0560205 
    ## Run 307 stress 0.05051329 
    ## ... Procrustes: rmse 0.03573437  max resid 0.1220846 
    ## Run 308 stress 0.05466702 
    ## Run 309 stress 0.0533759 
    ## Run 310 stress 0.3657347 
    ## Run 311 stress 0.0517106 
    ## Run 312 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001434437  max resid 0.000352216 
    ## ... Similar to previous best
    ## Run 313 stress 0.0540362 
    ## Run 314 stress 0.05036793 
    ## ... Procrustes: rmse 4.722083e-05  max resid 0.0001160787 
    ## ... Similar to previous best
    ## Run 315 stress 0.0590015 
    ## Run 316 stress 0.06136677 
    ## Run 317 stress 0.06136689 
    ## Run 318 stress 0.05466708 
    ## Run 319 stress 0.06027435 
    ## Run 320 stress 0.0505132 
    ## ... Procrustes: rmse 0.03566673  max resid 0.1217287 
    ## Run 321 stress 0.0630986 
    ## Run 322 stress 0.05979227 
    ## Run 323 stress 0.06476247 
    ## Run 324 stress 0.05602043 
    ## Run 325 stress 0.05051337 
    ## ... Procrustes: rmse 0.03570098  max resid 0.1219876 
    ## Run 326 stress 0.05051335 
    ## ... Procrustes: rmse 0.03573694  max resid 0.1219254 
    ## Run 327 stress 0.05466706 
    ## Run 328 stress 0.05403568 
    ## Run 329 stress 0.06743847 
    ## Run 330 stress 0.05171058 
    ## Run 331 stress 0.05171058 
    ## Run 332 stress 0.05829931 
    ## Run 333 stress 0.06027445 
    ## Run 334 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001467673  max resid 0.0003805519 
    ## ... Similar to previous best
    ## Run 335 stress 0.050368 
    ## ... Procrustes: rmse 0.0001378592  max resid 0.0003148498 
    ## ... Similar to previous best
    ## Run 336 stress 0.05968627 
    ## Run 337 stress 0.05403596 
    ## Run 338 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001269039  max resid 0.0002958293 
    ## ... Similar to previous best
    ## Run 339 stress 0.06231325 
    ## Run 340 stress 0.3594388 
    ## Run 341 stress 0.05403574 
    ## Run 342 stress 0.06136678 
    ## Run 343 stress 0.05171058 
    ## Run 344 stress 0.05036797 
    ## ... Procrustes: rmse 0.000126379  max resid 0.0003091482 
    ## ... Similar to previous best
    ## Run 345 stress 0.05337579 
    ## Run 346 stress 0.06674556 
    ## Run 347 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001722057  max resid 0.0004413843 
    ## ... Similar to previous best
    ## Run 348 stress 0.0517106 
    ## Run 349 stress 0.06231291 
    ## Run 350 stress 0.05968634 
    ## Run 351 stress 0.05036792 
    ## ... Procrustes: rmse 3.053462e-05  max resid 8.418151e-05 
    ## ... Similar to previous best
    ## Run 352 stress 0.0602743 
    ## Run 353 stress 0.07272462 
    ## Run 354 stress 0.05602054 
    ## Run 355 stress 0.05466704 
    ## Run 356 stress 0.05051341 
    ## ... Procrustes: rmse 0.03573265  max resid 0.1220831 
    ## Run 357 stress 0.06609066 
    ## Run 358 stress 0.05466704 
    ## Run 359 stress 0.050368 
    ## ... Procrustes: rmse 0.000131532  max resid 0.0003620548 
    ## ... Similar to previous best
    ## Run 360 stress 0.05036793 
    ## ... Procrustes: rmse 7.113326e-05  max resid 0.0001683454 
    ## ... Similar to previous best
    ## Run 361 stress 0.06231294 
    ## Run 362 stress 0.0603895 
    ## Run 363 stress 0.05602053 
    ## Run 364 stress 0.05466704 
    ## Run 365 stress 0.05051311 
    ## ... Procrustes: rmse 0.03570576  max resid 0.1219909 
    ## Run 366 stress 0.06609071 
    ## Run 367 stress 0.05036792 
    ## ... Procrustes: rmse 6.942745e-05  max resid 0.0001668077 
    ## ... Similar to previous best
    ## Run 368 stress 0.05900149 
    ## Run 369 stress 0.05036791 
    ## ... Procrustes: rmse 3.074539e-05  max resid 7.436383e-05 
    ## ... Similar to previous best
    ## Run 370 stress 0.06136684 
    ## Run 371 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001306002  max resid 0.0003380566 
    ## ... Similar to previous best
    ## Run 372 stress 0.05036804 
    ## ... Procrustes: rmse 0.000168169  max resid 0.0003972651 
    ## ... Similar to previous best
    ## Run 373 stress 0.0590015 
    ## Run 374 stress 0.05036793 
    ## ... Procrustes: rmse 5.573065e-05  max resid 0.0001214514 
    ## ... Similar to previous best
    ## Run 375 stress 0.05171062 
    ## Run 376 stress 0.05051341 
    ## ... Procrustes: rmse 0.03571209  max resid 0.121852 
    ## Run 377 stress 0.05466711 
    ## Run 378 stress 0.05051331 
    ## ... Procrustes: rmse 0.03571064  max resid 0.121853 
    ## Run 379 stress 0.05602055 
    ## Run 380 stress 0.05337584 
    ## Run 381 stress 0.05171059 
    ## Run 382 stress 0.05602044 
    ## Run 383 stress 0.06136686 
    ## Run 384 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001165171  max resid 0.0003144951 
    ## ... Similar to previous best
    ## Run 385 stress 0.06476225 
    ## Run 386 stress 0.06478538 
    ## Run 387 stress 0.05051311 
    ## ... Procrustes: rmse 0.03567903  max resid 0.1217749 
    ## Run 388 stress 0.0517106 
    ## Run 389 stress 0.0582993 
    ## Run 390 stress 0.05602043 
    ## Run 391 stress 0.05171064 
    ## Run 392 stress 0.050368 
    ## ... Procrustes: rmse 0.0001586101  max resid 0.0003933667 
    ## ... Similar to previous best
    ## Run 393 stress 0.05036797 
    ## ... Procrustes: rmse 8.84285e-05  max resid 0.0001779825 
    ## ... Similar to previous best
    ## Run 394 stress 0.05051299 
    ## ... Procrustes: rmse 0.03569405  max resid 0.1219487 
    ## Run 395 stress 0.05051313 
    ## ... Procrustes: rmse 0.03571174  max resid 0.1220099 
    ## Run 396 stress 0.07026201 
    ## Run 397 stress 0.0613668 
    ## Run 398 stress 0.05036791 
    ## ... Procrustes: rmse 4.264558e-05  max resid 0.0001031603 
    ## ... Similar to previous best
    ## Run 399 stress 0.0533758 
    ## Run 400 stress 0.050368 
    ## ... Procrustes: rmse 0.0001526457  max resid 0.0003651225 
    ## ... Similar to previous best
    ## Run 401 stress 0.05602062 
    ## Run 402 stress 0.05968646 
    ## Run 403 stress 0.05051346 
    ## ... Procrustes: rmse 0.03568126  max resid 0.1217625 
    ## Run 404 stress 0.05036794 
    ## ... Procrustes: rmse 9.85813e-05  max resid 0.0002420145 
    ## ... Similar to previous best
    ## Run 405 stress 0.05171065 
    ## Run 406 stress 0.05036792 
    ## ... Procrustes: rmse 4.861232e-05  max resid 0.000101141 
    ## ... Similar to previous best
    ## Run 407 stress 0.05829947 
    ## Run 408 stress 0.05829941 
    ## Run 409 stress 0.05051339 
    ## ... Procrustes: rmse 0.03570412  max resid 0.1219977 
    ## Run 410 stress 0.05051332 
    ## ... Procrustes: rmse 0.03570562  max resid 0.1219998 
    ## Run 411 stress 0.05051318 
    ## ... Procrustes: rmse 0.03576644  max resid 0.1220355 
    ## Run 412 stress 0.05036794 
    ## ... Procrustes: rmse 8.963445e-05  max resid 0.0002168462 
    ## ... Similar to previous best
    ## Run 413 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001371705  max resid 0.000336254 
    ## ... Similar to previous best
    ## Run 414 stress 0.05171061 
    ## Run 415 stress 0.07461542 
    ## Run 416 stress 0.06309853 
    ## Run 417 stress 0.05403594 
    ## Run 418 stress 0.05051322 
    ## ... Procrustes: rmse 0.03577953  max resid 0.1220743 
    ## Run 419 stress 0.05036794 
    ## ... Procrustes: rmse 7.298622e-05  max resid 0.000142848 
    ## ... Similar to previous best
    ## Run 420 stress 0.05036796 
    ## ... Procrustes: rmse 8.006281e-05  max resid 0.0001500516 
    ## ... Similar to previous best
    ## Run 421 stress 0.05051344 
    ## ... Procrustes: rmse 0.03577052  max resid 0.122197 
    ## Run 422 stress 0.05051306 
    ## ... Procrustes: rmse 0.03569459  max resid 0.1218269 
    ## Run 423 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.585248e-05  max resid 4.357567e-05 
    ## ... Similar to previous best
    ## Run 424 stress 0.05602047 
    ## Run 425 stress 0.05051303 
    ## ... Procrustes: rmse 0.03569442  max resid 0.1218266 
    ## Run 426 stress 0.05403604 
    ## Run 427 stress 0.05171061 
    ## Run 428 stress 0.050368 
    ## ... Procrustes: rmse 0.0001455284  max resid 0.0003561594 
    ## ... Similar to previous best
    ## Run 429 stress 0.05051326 
    ## ... Procrustes: rmse 0.03569851  max resid 0.1218165 
    ## Run 430 stress 0.05602058 
    ## Run 431 stress 0.05051323 
    ## ... Procrustes: rmse 0.03569892  max resid 0.1218184 
    ## Run 432 stress 0.05051364 
    ## ... Procrustes: rmse 0.03570495  max resid 0.1220052 
    ## Run 433 stress 0.05602057 
    ## Run 434 stress 0.0560205 
    ## Run 435 stress 0.05051339 
    ## ... Procrustes: rmse 0.03572448  max resid 0.1220534 
    ## Run 436 stress 0.05036796 
    ## ... Procrustes: rmse 9.958046e-05  max resid 0.0002416598 
    ## ... Similar to previous best
    ## Run 437 stress 0.05602048 
    ## Run 438 stress 0.05337584 
    ## Run 439 stress 0.05171065 
    ## Run 440 stress 0.05337578 
    ## Run 441 stress 0.06027427 
    ## Run 442 stress 0.05051339 
    ## ... Procrustes: rmse 0.03567163  max resid 0.1217241 
    ## Run 443 stress 0.05968638 
    ## Run 444 stress 0.05171059 
    ## Run 445 stress 0.05051325 
    ## ... Procrustes: rmse 0.03571264  max resid 0.1220156 
    ## Run 446 stress 0.0690449 
    ## Run 447 stress 0.06027441 
    ## Run 448 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001470312  max resid 0.0003608859 
    ## ... Similar to previous best
    ## Run 449 stress 0.05036795 
    ## ... Procrustes: rmse 8.825837e-05  max resid 0.0002081221 
    ## ... Similar to previous best
    ## Run 450 stress 0.05403617 
    ## Run 451 stress 0.05171058 
    ## Run 452 stress 0.06951031 
    ## Run 453 stress 0.05602046 
    ## Run 454 stress 0.05051315 
    ## ... Procrustes: rmse 0.03565695  max resid 0.1218438 
    ## Run 455 stress 0.05036793 
    ## ... Procrustes: rmse 7.722861e-05  max resid 0.000199168 
    ## ... Similar to previous best
    ## Run 456 stress 0.05900147 
    ## Run 457 stress 0.05171058 
    ## Run 458 stress 0.05036795 
    ## ... Procrustes: rmse 9.889529e-05  max resid 0.0002791519 
    ## ... Similar to previous best
    ## Run 459 stress 0.05466707 
    ## Run 460 stress 0.06640336 
    ## Run 461 stress 0.0505132 
    ## ... Procrustes: rmse 0.03561555  max resid 0.1217179 
    ## Run 462 stress 0.05171061 
    ## Run 463 stress 0.05051327 
    ## ... Procrustes: rmse 0.03573836  max resid 0.1220913 
    ## Run 464 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001325791  max resid 0.0003493032 
    ## ... Similar to previous best
    ## Run 465 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001173838  max resid 0.0003111161 
    ## ... Similar to previous best
    ## Run 466 stress 0.06027428 
    ## Run 467 stress 0.05051347 
    ## ... Procrustes: rmse 0.03573361  max resid 0.1220852 
    ## Run 468 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001418313  max resid 0.0003873524 
    ## ... Similar to previous best
    ## Run 469 stress 0.0503681 
    ## ... Procrustes: rmse 0.0001818867  max resid 0.0004583025 
    ## ... Similar to previous best
    ## Run 470 stress 0.05171059 
    ## Run 471 stress 0.05968641 
    ## Run 472 stress 0.05171063 
    ## Run 473 stress 0.05051322 
    ## ... Procrustes: rmse 0.03571279  max resid 0.1220132 
    ## Run 474 stress 0.05036802 
    ## ... Procrustes: rmse 0.000166989  max resid 0.0004355347 
    ## ... Similar to previous best
    ## Run 475 stress 0.06027436 
    ## Run 476 stress 0.05051338 
    ## ... Procrustes: rmse 0.03568849  max resid 0.1217719 
    ## Run 477 stress 0.05036793 
    ## ... Procrustes: rmse 5.649257e-05  max resid 0.000126559 
    ## ... Similar to previous best
    ## Run 478 stress 0.0505135 
    ## ... Procrustes: rmse 0.03572056  max resid 0.1218612 
    ## Run 479 stress 0.05602049 
    ## Run 480 stress 0.05928229 
    ## Run 481 stress 0.05829937 
    ## Run 482 stress 0.06478537 
    ## Run 483 stress 0.05051304 
    ## ... Procrustes: rmse 0.03569011  max resid 0.1219377 
    ## Run 484 stress 0.05051308 
    ## ... Procrustes: rmse 0.03569543  max resid 0.1219561 
    ## Run 485 stress 0.3555094 
    ## Run 486 stress 0.06476217 
    ## Run 487 stress 0.05036793 
    ## ... Procrustes: rmse 7.169582e-05  max resid 0.0001680439 
    ## ... Similar to previous best
    ## Run 488 stress 0.0582993 
    ## Run 489 stress 0.05051329 
    ## ... Procrustes: rmse 0.03572345  max resid 0.1220479 
    ## Run 490 stress 0.07026245 
    ## Run 491 stress 0.05928219 
    ## Run 492 stress 0.05171061 
    ## Run 493 stress 0.05337595 
    ## Run 494 stress 0.06499232 
    ## Run 495 stress 0.06674557 
    ## Run 496 stress 0.05171058 
    ## Run 497 stress 0.05466701 
    ## Run 498 stress 0.06478512 
    ## Run 499 stress 0.05051355 
    ## ... Procrustes: rmse 0.0357528  max resid 0.1219635 
    ## Run 500 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001287126  max resid 0.0003126197 
    ## ... Similar to previous best
    ## *** Best solution repeated 15 times

``` r
# Ocean sites and mixed lakes
PD_beta_env_OM_NMDS <- metaMDS(PD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08352634 
    ## Run 1 stress 0.08352634 
    ## ... Procrustes: rmse 7.733823e-06  max resid 1.727751e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.1134014 
    ## Run 3 stress 0.1134014 
    ## Run 4 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 2.305599e-06  max resid 5.051656e-06 
    ## ... Similar to previous best
    ## Run 5 stress 0.1134014 
    ## Run 6 stress 0.08352634 
    ## ... Procrustes: rmse 1.119888e-06  max resid 2.199118e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.08352634 
    ## ... Procrustes: rmse 1.707024e-06  max resid 3.293156e-06 
    ## ... Similar to previous best
    ## Run 8 stress 0.08352634 
    ## ... Procrustes: rmse 1.865666e-06  max resid 3.351341e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.08352634 
    ## ... Procrustes: rmse 1.234264e-06  max resid 2.543835e-06 
    ## ... Similar to previous best
    ## Run 10 stress 0.1134014 
    ## Run 11 stress 0.127931 
    ## Run 12 stress 0.1134015 
    ## Run 13 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.671679e-06  max resid 3.587133e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.08352634 
    ## ... Procrustes: rmse 3.740105e-06  max resid 8.218878e-06 
    ## ... Similar to previous best
    ## Run 15 stress 0.08352634 
    ## ... Procrustes: rmse 5.713204e-06  max resid 1.263301e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.08352634 
    ## ... Procrustes: rmse 2.459279e-06  max resid 5.30656e-06 
    ## ... Similar to previous best
    ## Run 17 stress 0.08352634 
    ## ... Procrustes: rmse 1.417698e-06  max resid 3.056717e-06 
    ## ... Similar to previous best
    ## Run 18 stress 0.1134014 
    ## Run 19 stress 0.08352634 
    ## ... Procrustes: rmse 1.197213e-06  max resid 2.476641e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.08352634 
    ## ... Procrustes: rmse 1.989859e-06  max resid 4.222867e-06 
    ## ... Similar to previous best
    ## Run 21 stress 0.08352634 
    ## ... Procrustes: rmse 1.260112e-06  max resid 2.635235e-06 
    ## ... Similar to previous best
    ## Run 22 stress 0.1134015 
    ## Run 23 stress 0.08352634 
    ## ... Procrustes: rmse 3.472703e-06  max resid 7.662238e-06 
    ## ... Similar to previous best
    ## Run 24 stress 0.1134015 
    ## Run 25 stress 0.08352634 
    ## ... Procrustes: rmse 8.263949e-07  max resid 1.274947e-06 
    ## ... Similar to previous best
    ## Run 26 stress 0.1279309 
    ## Run 27 stress 0.08352634 
    ## ... Procrustes: rmse 1.764966e-06  max resid 3.917643e-06 
    ## ... Similar to previous best
    ## Run 28 stress 0.08352634 
    ## ... Procrustes: rmse 3.476317e-06  max resid 7.63949e-06 
    ## ... Similar to previous best
    ## Run 29 stress 0.1279309 
    ## Run 30 stress 0.08352634 
    ## ... Procrustes: rmse 1.057372e-06  max resid 2.245276e-06 
    ## ... Similar to previous best
    ## Run 31 stress 0.3262791 
    ## Run 32 stress 0.1279309 
    ## Run 33 stress 0.08352634 
    ## ... Procrustes: rmse 6.582015e-07  max resid 1.292397e-06 
    ## ... Similar to previous best
    ## Run 34 stress 0.08352634 
    ## ... Procrustes: rmse 1.783646e-06  max resid 3.990808e-06 
    ## ... Similar to previous best
    ## Run 35 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 6.576339e-07  max resid 1.192561e-06 
    ## ... Similar to previous best
    ## Run 36 stress 0.1134015 
    ## Run 37 stress 0.08352634 
    ## ... Procrustes: rmse 2.282347e-06  max resid 4.92345e-06 
    ## ... Similar to previous best
    ## Run 38 stress 0.08352634 
    ## ... Procrustes: rmse 3.367351e-06  max resid 7.185343e-06 
    ## ... Similar to previous best
    ## Run 39 stress 0.1134014 
    ## Run 40 stress 0.1134014 
    ## Run 41 stress 0.08352634 
    ## ... Procrustes: rmse 3.109448e-06  max resid 6.956743e-06 
    ## ... Similar to previous best
    ## Run 42 stress 0.08352634 
    ## ... Procrustes: rmse 6.226739e-06  max resid 1.383852e-05 
    ## ... Similar to previous best
    ## Run 43 stress 0.08352634 
    ## ... Procrustes: rmse 3.448285e-06  max resid 7.601849e-06 
    ## ... Similar to previous best
    ## Run 44 stress 0.08352634 
    ## ... Procrustes: rmse 8.914921e-07  max resid 1.924407e-06 
    ## ... Similar to previous best
    ## Run 45 stress 0.08352634 
    ## ... Procrustes: rmse 5.112996e-07  max resid 9.523741e-07 
    ## ... Similar to previous best
    ## Run 46 stress 0.08352634 
    ## ... Procrustes: rmse 1.300031e-05  max resid 2.874308e-05 
    ## ... Similar to previous best
    ## Run 47 stress 0.1279309 
    ## Run 48 stress 0.08352634 
    ## ... Procrustes: rmse 2.513555e-06  max resid 5.523931e-06 
    ## ... Similar to previous best
    ## Run 49 stress 0.08352634 
    ## ... Procrustes: rmse 8.299899e-07  max resid 1.775955e-06 
    ## ... Similar to previous best
    ## Run 50 stress 0.08352634 
    ## ... Procrustes: rmse 3.755687e-06  max resid 7.497211e-06 
    ## ... Similar to previous best
    ## Run 51 stress 0.08352634 
    ## ... Procrustes: rmse 2.974446e-06  max resid 5.803355e-06 
    ## ... Similar to previous best
    ## Run 52 stress 0.1134015 
    ## Run 53 stress 0.127931 
    ## Run 54 stress 0.1279309 
    ## Run 55 stress 0.1134014 
    ## Run 56 stress 0.08352634 
    ## ... Procrustes: rmse 1.270338e-06  max resid 2.675727e-06 
    ## ... Similar to previous best
    ## Run 57 stress 0.08352634 
    ## ... Procrustes: rmse 3.000931e-06  max resid 6.208494e-06 
    ## ... Similar to previous best
    ## Run 58 stress 0.1134014 
    ## Run 59 stress 0.127931 
    ## Run 60 stress 0.08352634 
    ## ... Procrustes: rmse 2.385798e-06  max resid 5.305082e-06 
    ## ... Similar to previous best
    ## Run 61 stress 0.08352634 
    ## ... Procrustes: rmse 2.636124e-06  max resid 5.509212e-06 
    ## ... Similar to previous best
    ## Run 62 stress 0.08352634 
    ## ... Procrustes: rmse 1.814646e-06  max resid 3.917434e-06 
    ## ... Similar to previous best
    ## Run 63 stress 0.1989488 
    ## Run 64 stress 0.08352634 
    ## ... Procrustes: rmse 2.215299e-06  max resid 4.30578e-06 
    ## ... Similar to previous best
    ## Run 65 stress 0.1134014 
    ## Run 66 stress 0.08352634 
    ## ... Procrustes: rmse 1.663642e-06  max resid 3.777261e-06 
    ## ... Similar to previous best
    ## Run 67 stress 0.1783142 
    ## Run 68 stress 0.08352634 
    ## ... Procrustes: rmse 4.616081e-07  max resid 7.491856e-07 
    ## ... Similar to previous best
    ## Run 69 stress 0.08352634 
    ## ... Procrustes: rmse 3.059274e-06  max resid 6.912392e-06 
    ## ... Similar to previous best
    ## Run 70 stress 0.08352634 
    ## ... Procrustes: rmse 3.09606e-06  max resid 6.95058e-06 
    ## ... Similar to previous best
    ## Run 71 stress 0.08352634 
    ## ... Procrustes: rmse 3.39117e-06  max resid 7.558041e-06 
    ## ... Similar to previous best
    ## Run 72 stress 0.08352634 
    ## ... Procrustes: rmse 6.050095e-06  max resid 1.340005e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.1134014 
    ## Run 74 stress 0.1279309 
    ## Run 75 stress 0.08352634 
    ## ... Procrustes: rmse 2.144403e-06  max resid 4.758277e-06 
    ## ... Similar to previous best
    ## Run 76 stress 0.1279309 
    ## Run 77 stress 0.08352634 
    ## ... Procrustes: rmse 1.527827e-06  max resid 3.361948e-06 
    ## ... Similar to previous best
    ## Run 78 stress 0.127931 
    ## Run 79 stress 0.1134015 
    ## Run 80 stress 0.1134015 
    ## Run 81 stress 0.08352634 
    ## ... Procrustes: rmse 3.609075e-06  max resid 8.100256e-06 
    ## ... Similar to previous best
    ## Run 82 stress 0.08352634 
    ## ... Procrustes: rmse 2.330232e-06  max resid 5.131211e-06 
    ## ... Similar to previous best
    ## Run 83 stress 0.1134014 
    ## Run 84 stress 0.08352634 
    ## ... Procrustes: rmse 3.050211e-06  max resid 6.730502e-06 
    ## ... Similar to previous best
    ## Run 85 stress 0.08352634 
    ## ... Procrustes: rmse 7.075948e-07  max resid 1.222995e-06 
    ## ... Similar to previous best
    ## Run 86 stress 0.08352634 
    ## ... Procrustes: rmse 3.230242e-06  max resid 7.107085e-06 
    ## ... Similar to previous best
    ## Run 87 stress 0.1279309 
    ## Run 88 stress 0.08352634 
    ## ... Procrustes: rmse 8.382422e-07  max resid 1.871555e-06 
    ## ... Similar to previous best
    ## Run 89 stress 0.127931 
    ## Run 90 stress 0.08352634 
    ## ... Procrustes: rmse 3.136184e-06  max resid 7.114628e-06 
    ## ... Similar to previous best
    ## Run 91 stress 0.1279309 
    ## Run 92 stress 0.08352634 
    ## ... Procrustes: rmse 2.351437e-06  max resid 5.309414e-06 
    ## ... Similar to previous best
    ## Run 93 stress 0.08352634 
    ## ... Procrustes: rmse 1.736568e-06  max resid 3.844307e-06 
    ## ... Similar to previous best
    ## Run 94 stress 0.1134015 
    ## Run 95 stress 0.08352634 
    ## ... Procrustes: rmse 1.493925e-06  max resid 2.360087e-06 
    ## ... Similar to previous best
    ## Run 96 stress 0.1134014 
    ## Run 97 stress 0.08352634 
    ## ... Procrustes: rmse 1.780293e-06  max resid 3.944098e-06 
    ## ... Similar to previous best
    ## Run 98 stress 0.08352634 
    ## ... Procrustes: rmse 2.195965e-06  max resid 4.680906e-06 
    ## ... Similar to previous best
    ## Run 99 stress 0.08352634 
    ## ... Procrustes: rmse 7.737002e-07  max resid 1.274026e-06 
    ## ... Similar to previous best
    ## Run 100 stress 0.08352634 
    ## ... Procrustes: rmse 2.625785e-06  max resid 5.853663e-06 
    ## ... Similar to previous best
    ## Run 101 stress 0.08352634 
    ## ... Procrustes: rmse 3.502724e-06  max resid 6.91011e-06 
    ## ... Similar to previous best
    ## Run 102 stress 0.1134014 
    ## Run 103 stress 0.1134015 
    ## Run 104 stress 0.127931 
    ## Run 105 stress 0.1279309 
    ## Run 106 stress 0.1279309 
    ## Run 107 stress 0.08352634 
    ## ... Procrustes: rmse 2.653438e-06  max resid 5.746006e-06 
    ## ... Similar to previous best
    ## Run 108 stress 0.1134014 
    ## Run 109 stress 0.08352634 
    ## ... Procrustes: rmse 3.998651e-06  max resid 8.711515e-06 
    ## ... Similar to previous best
    ## Run 110 stress 0.08352634 
    ## ... Procrustes: rmse 3.544225e-06  max resid 6.595682e-06 
    ## ... Similar to previous best
    ## Run 111 stress 0.1279309 
    ## Run 112 stress 0.1279309 
    ## Run 113 stress 0.08352634 
    ## ... Procrustes: rmse 2.616767e-06  max resid 5.670971e-06 
    ## ... Similar to previous best
    ## Run 114 stress 0.1134015 
    ## Run 115 stress 0.08352634 
    ## ... Procrustes: rmse 3.715912e-06  max resid 8.294748e-06 
    ## ... Similar to previous best
    ## Run 116 stress 0.1134014 
    ## Run 117 stress 0.08352634 
    ## ... Procrustes: rmse 2.983546e-06  max resid 5.587903e-06 
    ## ... Similar to previous best
    ## Run 118 stress 0.1134014 
    ## Run 119 stress 0.08352634 
    ## ... Procrustes: rmse 8.41259e-07  max resid 1.878531e-06 
    ## ... Similar to previous best
    ## Run 120 stress 0.08352634 
    ## ... Procrustes: rmse 1.60893e-06  max resid 3.232274e-06 
    ## ... Similar to previous best
    ## Run 121 stress 0.08352634 
    ## ... Procrustes: rmse 2.733899e-06  max resid 6.154329e-06 
    ## ... Similar to previous best
    ## Run 122 stress 0.08352634 
    ## ... Procrustes: rmse 3.775882e-07  max resid 7.549474e-07 
    ## ... Similar to previous best
    ## Run 123 stress 0.08352634 
    ## ... Procrustes: rmse 3.459045e-06  max resid 7.868324e-06 
    ## ... Similar to previous best
    ## Run 124 stress 0.08352634 
    ## ... Procrustes: rmse 8.332668e-07  max resid 1.773205e-06 
    ## ... Similar to previous best
    ## Run 125 stress 0.08352634 
    ## ... Procrustes: rmse 2.157274e-06  max resid 4.730164e-06 
    ## ... Similar to previous best
    ## Run 126 stress 0.08352634 
    ## ... Procrustes: rmse 1.273592e-06  max resid 2.320904e-06 
    ## ... Similar to previous best
    ## Run 127 stress 0.08352634 
    ## ... Procrustes: rmse 3.186784e-06  max resid 7.061616e-06 
    ## ... Similar to previous best
    ## Run 128 stress 0.08352634 
    ## ... Procrustes: rmse 9.057673e-07  max resid 1.822529e-06 
    ## ... Similar to previous best
    ## Run 129 stress 0.08352634 
    ## ... Procrustes: rmse 2.67382e-06  max resid 5.876967e-06 
    ## ... Similar to previous best
    ## Run 130 stress 0.1279309 
    ## Run 131 stress 0.08352634 
    ## ... Procrustes: rmse 2.498076e-06  max resid 5.34394e-06 
    ## ... Similar to previous best
    ## Run 132 stress 0.1134015 
    ## Run 133 stress 0.08352634 
    ## ... Procrustes: rmse 3.497164e-06  max resid 7.817363e-06 
    ## ... Similar to previous best
    ## Run 134 stress 0.1279309 
    ## Run 135 stress 0.1134014 
    ## Run 136 stress 0.1134014 
    ## Run 137 stress 0.08352634 
    ## ... Procrustes: rmse 1.267475e-06  max resid 2.692612e-06 
    ## ... Similar to previous best
    ## Run 138 stress 0.08352634 
    ## ... Procrustes: rmse 6.573422e-07  max resid 1.03407e-06 
    ## ... Similar to previous best
    ## Run 139 stress 0.08352634 
    ## ... Procrustes: rmse 2.360326e-06  max resid 5.340321e-06 
    ## ... Similar to previous best
    ## Run 140 stress 0.08352634 
    ## ... Procrustes: rmse 3.065839e-06  max resid 6.357144e-06 
    ## ... Similar to previous best
    ## Run 141 stress 0.1279309 
    ## Run 142 stress 0.08352634 
    ## ... Procrustes: rmse 6.532916e-07  max resid 1.299322e-06 
    ## ... Similar to previous best
    ## Run 143 stress 0.08352634 
    ## ... Procrustes: rmse 1.820472e-06  max resid 4.027328e-06 
    ## ... Similar to previous best
    ## Run 144 stress 0.1134014 
    ## Run 145 stress 0.08352634 
    ## ... Procrustes: rmse 1.914816e-06  max resid 4.171283e-06 
    ## ... Similar to previous best
    ## Run 146 stress 0.08352634 
    ## ... Procrustes: rmse 5.41123e-06  max resid 1.232227e-05 
    ## ... Similar to previous best
    ## Run 147 stress 0.1279309 
    ## Run 148 stress 0.08352634 
    ## ... Procrustes: rmse 2.526587e-06  max resid 5.66211e-06 
    ## ... Similar to previous best
    ## Run 149 stress 0.08352634 
    ## ... Procrustes: rmse 1.035087e-06  max resid 2.29133e-06 
    ## ... Similar to previous best
    ## Run 150 stress 0.1783104 
    ## Run 151 stress 0.08352634 
    ## ... Procrustes: rmse 2.526986e-06  max resid 5.473644e-06 
    ## ... Similar to previous best
    ## Run 152 stress 0.1134014 
    ## Run 153 stress 0.127931 
    ## Run 154 stress 0.1134014 
    ## Run 155 stress 0.08352634 
    ## ... Procrustes: rmse 7.010664e-06  max resid 1.502604e-05 
    ## ... Similar to previous best
    ## Run 156 stress 0.08352634 
    ## ... Procrustes: rmse 2.830226e-06  max resid 6.243757e-06 
    ## ... Similar to previous best
    ## Run 157 stress 0.1134015 
    ## Run 158 stress 0.1989154 
    ## Run 159 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 2.69235e-07  max resid 4.575343e-07 
    ## ... Similar to previous best
    ## Run 160 stress 0.08352634 
    ## ... Procrustes: rmse 3.536273e-06  max resid 7.848421e-06 
    ## ... Similar to previous best
    ## Run 161 stress 0.1279309 
    ## Run 162 stress 0.1279309 
    ## Run 163 stress 0.08352634 
    ## ... Procrustes: rmse 1.432183e-06  max resid 2.755174e-06 
    ## ... Similar to previous best
    ## Run 164 stress 0.1134014 
    ## Run 165 stress 0.08352634 
    ## ... Procrustes: rmse 2.230034e-06  max resid 4.818549e-06 
    ## ... Similar to previous best
    ## Run 166 stress 0.08352634 
    ## ... Procrustes: rmse 2.525256e-06  max resid 5.392486e-06 
    ## ... Similar to previous best
    ## Run 167 stress 0.1134014 
    ## Run 168 stress 0.1134014 
    ## Run 169 stress 0.1279309 
    ## Run 170 stress 0.08352634 
    ## ... Procrustes: rmse 1.775695e-06  max resid 3.671883e-06 
    ## ... Similar to previous best
    ## Run 171 stress 0.08352634 
    ## ... Procrustes: rmse 1.487601e-06  max resid 3.109978e-06 
    ## ... Similar to previous best
    ## Run 172 stress 0.08352634 
    ## ... Procrustes: rmse 1.359692e-06  max resid 2.257096e-06 
    ## ... Similar to previous best
    ## Run 173 stress 0.1134015 
    ## Run 174 stress 0.08352634 
    ## ... Procrustes: rmse 2.979582e-06  max resid 6.613172e-06 
    ## ... Similar to previous best
    ## Run 175 stress 0.1279309 
    ## Run 176 stress 0.08352634 
    ## ... Procrustes: rmse 1.953458e-06  max resid 4.346483e-06 
    ## ... Similar to previous best
    ## Run 177 stress 0.1279309 
    ## Run 178 stress 0.08352634 
    ## ... Procrustes: rmse 9.200295e-07  max resid 2.020632e-06 
    ## ... Similar to previous best
    ## Run 179 stress 0.1989159 
    ## Run 180 stress 0.08352634 
    ## ... Procrustes: rmse 2.629511e-06  max resid 5.706445e-06 
    ## ... Similar to previous best
    ## Run 181 stress 0.1134014 
    ## Run 182 stress 0.1279309 
    ## Run 183 stress 0.1134014 
    ## Run 184 stress 0.1279309 
    ## Run 185 stress 0.1279309 
    ## Run 186 stress 0.08352634 
    ## ... Procrustes: rmse 4.619852e-06  max resid 1.017887e-05 
    ## ... Similar to previous best
    ## Run 187 stress 0.08352634 
    ## ... Procrustes: rmse 1.68497e-06  max resid 2.865858e-06 
    ## ... Similar to previous best
    ## Run 188 stress 0.08352634 
    ## ... Procrustes: rmse 2.764656e-06  max resid 5.418399e-06 
    ## ... Similar to previous best
    ## Run 189 stress 0.1279309 
    ## Run 190 stress 0.1279309 
    ## Run 191 stress 0.127931 
    ## Run 192 stress 0.08352634 
    ## ... Procrustes: rmse 1.376175e-06  max resid 2.502039e-06 
    ## ... Similar to previous best
    ## Run 193 stress 0.1134015 
    ## Run 194 stress 0.08352634 
    ## ... Procrustes: rmse 1.82204e-06  max resid 4.011032e-06 
    ## ... Similar to previous best
    ## Run 195 stress 0.08352634 
    ## ... Procrustes: rmse 1.710456e-06  max resid 3.112253e-06 
    ## ... Similar to previous best
    ## Run 196 stress 0.1279309 
    ## Run 197 stress 0.08352634 
    ## ... Procrustes: rmse 1.371591e-06  max resid 2.440551e-06 
    ## ... Similar to previous best
    ## Run 198 stress 0.08352634 
    ## ... Procrustes: rmse 1.361947e-06  max resid 2.845503e-06 
    ## ... Similar to previous best
    ## Run 199 stress 0.1134015 
    ## Run 200 stress 0.08352634 
    ## ... Procrustes: rmse 3.971429e-06  max resid 7.383899e-06 
    ## ... Similar to previous best
    ## Run 201 stress 0.08352634 
    ## ... Procrustes: rmse 6.912305e-07  max resid 1.369633e-06 
    ## ... Similar to previous best
    ## Run 202 stress 0.08352634 
    ## ... Procrustes: rmse 9.323356e-06  max resid 2.004287e-05 
    ## ... Similar to previous best
    ## Run 203 stress 0.08352634 
    ## ... Procrustes: rmse 2.925689e-06  max resid 5.684547e-06 
    ## ... Similar to previous best
    ## Run 204 stress 0.1134015 
    ## Run 205 stress 0.1279309 
    ## Run 206 stress 0.08352634 
    ## ... Procrustes: rmse 5.545373e-07  max resid 1.10493e-06 
    ## ... Similar to previous best
    ## Run 207 stress 0.08352634 
    ## ... Procrustes: rmse 2.699687e-06  max resid 5.988832e-06 
    ## ... Similar to previous best
    ## Run 208 stress 0.08352634 
    ## ... Procrustes: rmse 2.932429e-06  max resid 6.500173e-06 
    ## ... Similar to previous best
    ## Run 209 stress 0.127931 
    ## Run 210 stress 0.08352634 
    ## ... Procrustes: rmse 6.047709e-07  max resid 1.318248e-06 
    ## ... Similar to previous best
    ## Run 211 stress 0.1134015 
    ## Run 212 stress 0.08352634 
    ## ... Procrustes: rmse 7.116342e-07  max resid 1.45399e-06 
    ## ... Similar to previous best
    ## Run 213 stress 0.08352634 
    ## ... Procrustes: rmse 1.286891e-06  max resid 2.57032e-06 
    ## ... Similar to previous best
    ## Run 214 stress 0.08352634 
    ## ... Procrustes: rmse 2.322283e-06  max resid 5.077077e-06 
    ## ... Similar to previous best
    ## Run 215 stress 0.08352634 
    ## ... Procrustes: rmse 1.962448e-06  max resid 4.34428e-06 
    ## ... Similar to previous best
    ## Run 216 stress 0.08352634 
    ## ... Procrustes: rmse 5.740036e-06  max resid 1.293069e-05 
    ## ... Similar to previous best
    ## Run 217 stress 0.08352634 
    ## ... Procrustes: rmse 3.275006e-06  max resid 7.247768e-06 
    ## ... Similar to previous best
    ## Run 218 stress 0.08352634 
    ## ... Procrustes: rmse 2.710541e-06  max resid 6.551448e-06 
    ## ... Similar to previous best
    ## Run 219 stress 0.1134014 
    ## Run 220 stress 0.08352634 
    ## ... Procrustes: rmse 3.21787e-06  max resid 6.983886e-06 
    ## ... Similar to previous best
    ## Run 221 stress 0.1134015 
    ## Run 222 stress 0.1783138 
    ## Run 223 stress 0.08352634 
    ## ... Procrustes: rmse 3.308189e-06  max resid 7.348768e-06 
    ## ... Similar to previous best
    ## Run 224 stress 0.1279309 
    ## Run 225 stress 0.1989489 
    ## Run 226 stress 0.08352634 
    ## ... Procrustes: rmse 5.479275e-07  max resid 1.158406e-06 
    ## ... Similar to previous best
    ## Run 227 stress 0.08352634 
    ## ... Procrustes: rmse 1.354486e-06  max resid 3.002861e-06 
    ## ... Similar to previous best
    ## Run 228 stress 0.1134015 
    ## Run 229 stress 0.08352634 
    ## ... Procrustes: rmse 4.926867e-06  max resid 1.034479e-05 
    ## ... Similar to previous best
    ## Run 230 stress 0.1134015 
    ## Run 231 stress 0.08352634 
    ## ... Procrustes: rmse 3.134416e-06  max resid 7.226461e-06 
    ## ... Similar to previous best
    ## Run 232 stress 0.1989159 
    ## Run 233 stress 0.1134014 
    ## Run 234 stress 0.08352634 
    ## ... Procrustes: rmse 2.216424e-06  max resid 4.862779e-06 
    ## ... Similar to previous best
    ## Run 235 stress 0.1279309 
    ## Run 236 stress 0.08352634 
    ## ... Procrustes: rmse 5.903979e-07  max resid 1.104736e-06 
    ## ... Similar to previous best
    ## Run 237 stress 0.1134014 
    ## Run 238 stress 0.08352634 
    ## ... Procrustes: rmse 1.5245e-06  max resid 3.42322e-06 
    ## ... Similar to previous best
    ## Run 239 stress 0.08352634 
    ## ... Procrustes: rmse 9.65277e-07  max resid 1.864829e-06 
    ## ... Similar to previous best
    ## Run 240 stress 0.127931 
    ## Run 241 stress 0.08352634 
    ## ... Procrustes: rmse 2.944705e-06  max resid 6.521032e-06 
    ## ... Similar to previous best
    ## Run 242 stress 0.08352634 
    ## ... Procrustes: rmse 4.336089e-06  max resid 9.680459e-06 
    ## ... Similar to previous best
    ## Run 243 stress 0.08352634 
    ## ... Procrustes: rmse 6.560974e-07  max resid 1.464195e-06 
    ## ... Similar to previous best
    ## Run 244 stress 0.08352634 
    ## ... Procrustes: rmse 3.059502e-06  max resid 6.725915e-06 
    ## ... Similar to previous best
    ## Run 245 stress 0.1989156 
    ## Run 246 stress 0.08352634 
    ## ... Procrustes: rmse 4.034982e-06  max resid 8.683852e-06 
    ## ... Similar to previous best
    ## Run 247 stress 0.1134014 
    ## Run 248 stress 0.08352634 
    ## ... Procrustes: rmse 2.921046e-06  max resid 6.479869e-06 
    ## ... Similar to previous best
    ## Run 249 stress 0.08352634 
    ## ... Procrustes: rmse 3.450491e-06  max resid 6.395109e-06 
    ## ... Similar to previous best
    ## Run 250 stress 0.1134015 
    ## Run 251 stress 0.08352634 
    ## ... Procrustes: rmse 8.536451e-07  max resid 1.485762e-06 
    ## ... Similar to previous best
    ## Run 252 stress 0.08352634 
    ## ... Procrustes: rmse 1.688739e-06  max resid 3.630897e-06 
    ## ... Similar to previous best
    ## Run 253 stress 0.08352634 
    ## ... Procrustes: rmse 1.07739e-06  max resid 2.162092e-06 
    ## ... Similar to previous best
    ## Run 254 stress 0.1134014 
    ## Run 255 stress 0.1279309 
    ## Run 256 stress 0.08352634 
    ## ... Procrustes: rmse 1.324947e-06  max resid 2.102636e-06 
    ## ... Similar to previous best
    ## Run 257 stress 0.08352634 
    ## ... Procrustes: rmse 6.254777e-07  max resid 1.199128e-06 
    ## ... Similar to previous best
    ## Run 258 stress 0.1134014 
    ## Run 259 stress 0.08352634 
    ## ... Procrustes: rmse 1.368763e-06  max resid 2.93679e-06 
    ## ... Similar to previous best
    ## Run 260 stress 0.08352634 
    ## ... Procrustes: rmse 3.027917e-06  max resid 6.708378e-06 
    ## ... Similar to previous best
    ## Run 261 stress 0.08352634 
    ## ... Procrustes: rmse 4.641162e-06  max resid 9.1923e-06 
    ## ... Similar to previous best
    ## Run 262 stress 0.08352634 
    ## ... Procrustes: rmse 1.677718e-06  max resid 3.635735e-06 
    ## ... Similar to previous best
    ## Run 263 stress 0.1134015 
    ## Run 264 stress 0.08352634 
    ## ... Procrustes: rmse 1.181095e-06  max resid 2.458168e-06 
    ## ... Similar to previous best
    ## Run 265 stress 0.08352634 
    ## ... Procrustes: rmse 2.591751e-06  max resid 5.744462e-06 
    ## ... Similar to previous best
    ## Run 266 stress 0.08352634 
    ## ... Procrustes: rmse 1.828377e-06  max resid 3.998514e-06 
    ## ... Similar to previous best
    ## Run 267 stress 0.08352634 
    ## ... Procrustes: rmse 2.210382e-06  max resid 4.84201e-06 
    ## ... Similar to previous best
    ## Run 268 stress 0.1134014 
    ## Run 269 stress 0.08352634 
    ## ... Procrustes: rmse 1.884452e-06  max resid 4.185148e-06 
    ## ... Similar to previous best
    ## Run 270 stress 0.1279309 
    ## Run 271 stress 0.178312 
    ## Run 272 stress 0.08352634 
    ## ... Procrustes: rmse 3.893786e-06  max resid 8.545403e-06 
    ## ... Similar to previous best
    ## Run 273 stress 0.08352634 
    ## ... Procrustes: rmse 2.056841e-06  max resid 4.498255e-06 
    ## ... Similar to previous best
    ## Run 274 stress 0.08352634 
    ## ... Procrustes: rmse 1.329328e-06  max resid 2.932492e-06 
    ## ... Similar to previous best
    ## Run 275 stress 0.1134014 
    ## Run 276 stress 0.3121569 
    ## Run 277 stress 0.1279309 
    ## Run 278 stress 0.08352634 
    ## ... Procrustes: rmse 3.419233e-06  max resid 7.575767e-06 
    ## ... Similar to previous best
    ## Run 279 stress 0.08352634 
    ## ... Procrustes: rmse 1.400313e-06  max resid 2.953613e-06 
    ## ... Similar to previous best
    ## Run 280 stress 0.1134015 
    ## Run 281 stress 0.1134015 
    ## Run 282 stress 0.127931 
    ## Run 283 stress 0.08352634 
    ## ... Procrustes: rmse 1.119798e-06  max resid 2.335967e-06 
    ## ... Similar to previous best
    ## Run 284 stress 0.08352634 
    ## ... Procrustes: rmse 3.883978e-06  max resid 8.900927e-06 
    ## ... Similar to previous best
    ## Run 285 stress 0.08352634 
    ## ... Procrustes: rmse 2.937392e-06  max resid 6.505487e-06 
    ## ... Similar to previous best
    ## Run 286 stress 0.08352634 
    ## ... Procrustes: rmse 2.106389e-06  max resid 4.677098e-06 
    ## ... Similar to previous best
    ## Run 287 stress 0.1134014 
    ## Run 288 stress 0.08352634 
    ## ... Procrustes: rmse 8.900872e-07  max resid 2.111449e-06 
    ## ... Similar to previous best
    ## Run 289 stress 0.1989485 
    ## Run 290 stress 0.08352634 
    ## ... Procrustes: rmse 4.124767e-06  max resid 9.170373e-06 
    ## ... Similar to previous best
    ## Run 291 stress 0.1989154 
    ## Run 292 stress 0.1134014 
    ## Run 293 stress 0.08352634 
    ## ... Procrustes: rmse 2.292671e-06  max resid 5.085286e-06 
    ## ... Similar to previous best
    ## Run 294 stress 0.1989158 
    ## Run 295 stress 0.127931 
    ## Run 296 stress 0.08352634 
    ## ... Procrustes: rmse 5.070975e-07  max resid 8.140519e-07 
    ## ... Similar to previous best
    ## Run 297 stress 0.1134014 
    ## Run 298 stress 0.1134014 
    ## Run 299 stress 0.08352634 
    ## ... Procrustes: rmse 4.623946e-07  max resid 7.906876e-07 
    ## ... Similar to previous best
    ## Run 300 stress 0.1279309 
    ## Run 301 stress 0.08352634 
    ## ... Procrustes: rmse 1.890799e-06  max resid 4.170452e-06 
    ## ... Similar to previous best
    ## Run 302 stress 0.08352634 
    ## ... Procrustes: rmse 3.744911e-06  max resid 8.259911e-06 
    ## ... Similar to previous best
    ## Run 303 stress 0.08352634 
    ## ... Procrustes: rmse 2.531284e-06  max resid 5.56958e-06 
    ## ... Similar to previous best
    ## Run 304 stress 0.08352634 
    ## ... Procrustes: rmse 2.308876e-06  max resid 4.050486e-06 
    ## ... Similar to previous best
    ## Run 305 stress 0.1989489 
    ## Run 306 stress 0.1279309 
    ## Run 307 stress 0.08352634 
    ## ... Procrustes: rmse 1.807156e-06  max resid 4.025547e-06 
    ## ... Similar to previous best
    ## Run 308 stress 0.3270104 
    ## Run 309 stress 0.08352634 
    ## ... Procrustes: rmse 1.609294e-06  max resid 3.548406e-06 
    ## ... Similar to previous best
    ## Run 310 stress 0.08352634 
    ## ... Procrustes: rmse 1.338747e-06  max resid 2.929776e-06 
    ## ... Similar to previous best
    ## Run 311 stress 0.08352634 
    ## ... Procrustes: rmse 2.04226e-06  max resid 4.61052e-06 
    ## ... Similar to previous best
    ## Run 312 stress 0.08352634 
    ## ... Procrustes: rmse 2.425983e-06  max resid 5.28259e-06 
    ## ... Similar to previous best
    ## Run 313 stress 0.1989156 
    ## Run 314 stress 0.1134015 
    ## Run 315 stress 0.08352634 
    ## ... Procrustes: rmse 4.541915e-06  max resid 9.862297e-06 
    ## ... Similar to previous best
    ## Run 316 stress 0.08352634 
    ## ... Procrustes: rmse 5.938411e-07  max resid 9.727615e-07 
    ## ... Similar to previous best
    ## Run 317 stress 0.1134014 
    ## Run 318 stress 0.1134014 
    ## Run 319 stress 0.3266727 
    ## Run 320 stress 0.08352634 
    ## ... Procrustes: rmse 9.870722e-07  max resid 1.511905e-06 
    ## ... Similar to previous best
    ## Run 321 stress 0.1134014 
    ## Run 322 stress 0.08352634 
    ## ... Procrustes: rmse 2.410025e-06  max resid 5.403344e-06 
    ## ... Similar to previous best
    ## Run 323 stress 0.1279309 
    ## Run 324 stress 0.1279309 
    ## Run 325 stress 0.1279309 
    ## Run 326 stress 0.1134014 
    ## Run 327 stress 0.1279309 
    ## Run 328 stress 0.08352634 
    ## ... Procrustes: rmse 2.292896e-06  max resid 4.97357e-06 
    ## ... Similar to previous best
    ## Run 329 stress 0.08352634 
    ## ... Procrustes: rmse 3.732266e-06  max resid 8.524759e-06 
    ## ... Similar to previous best
    ## Run 330 stress 0.1989156 
    ## Run 331 stress 0.08352634 
    ## ... Procrustes: rmse 3.132194e-06  max resid 7.009598e-06 
    ## ... Similar to previous best
    ## Run 332 stress 0.08352634 
    ## ... Procrustes: rmse 2.849595e-06  max resid 6.250204e-06 
    ## ... Similar to previous best
    ## Run 333 stress 0.08352634 
    ## ... Procrustes: rmse 2.853446e-06  max resid 6.288228e-06 
    ## ... Similar to previous best
    ## Run 334 stress 0.1279309 
    ## Run 335 stress 0.2031894 
    ## Run 336 stress 0.1134014 
    ## Run 337 stress 0.08352634 
    ## ... Procrustes: rmse 1.2442e-06  max resid 2.499711e-06 
    ## ... Similar to previous best
    ## Run 338 stress 0.1134015 
    ## Run 339 stress 0.1989154 
    ## Run 340 stress 0.08352634 
    ## ... Procrustes: rmse 6.29301e-06  max resid 1.398251e-05 
    ## ... Similar to previous best
    ## Run 341 stress 0.1279309 
    ## Run 342 stress 0.08352634 
    ## ... Procrustes: rmse 8.634378e-07  max resid 1.761373e-06 
    ## ... Similar to previous best
    ## Run 343 stress 0.1279309 
    ## Run 344 stress 0.1279309 
    ## Run 345 stress 0.1134014 
    ## Run 346 stress 0.08352634 
    ## ... Procrustes: rmse 2.524193e-06  max resid 5.570217e-06 
    ## ... Similar to previous best
    ## Run 347 stress 0.1989154 
    ## Run 348 stress 0.1134014 
    ## Run 349 stress 0.08352634 
    ## ... Procrustes: rmse 2.396498e-06  max resid 5.381756e-06 
    ## ... Similar to previous best
    ## Run 350 stress 0.1134014 
    ## Run 351 stress 0.2034292 
    ## Run 352 stress 0.08352634 
    ## ... Procrustes: rmse 4.352058e-07  max resid 9.291812e-07 
    ## ... Similar to previous best
    ## Run 353 stress 0.1134015 
    ## Run 354 stress 0.1279309 
    ## Run 355 stress 0.08352634 
    ## ... Procrustes: rmse 1.924978e-06  max resid 4.180667e-06 
    ## ... Similar to previous best
    ## Run 356 stress 0.08352634 
    ## ... Procrustes: rmse 1.478953e-06  max resid 3.27805e-06 
    ## ... Similar to previous best
    ## Run 357 stress 0.08352634 
    ## ... Procrustes: rmse 2.121591e-06  max resid 4.705266e-06 
    ## ... Similar to previous best
    ## Run 358 stress 0.1134014 
    ## Run 359 stress 0.08352634 
    ## ... Procrustes: rmse 4.087269e-06  max resid 8.997017e-06 
    ## ... Similar to previous best
    ## Run 360 stress 0.08352634 
    ## ... Procrustes: rmse 2.951076e-06  max resid 6.501449e-06 
    ## ... Similar to previous best
    ## Run 361 stress 0.1134014 
    ## Run 362 stress 0.1989155 
    ## Run 363 stress 0.08352634 
    ## ... Procrustes: rmse 1.82661e-06  max resid 4.086688e-06 
    ## ... Similar to previous best
    ## Run 364 stress 0.1279309 
    ## Run 365 stress 0.127931 
    ## Run 366 stress 0.1134015 
    ## Run 367 stress 0.08352634 
    ## ... Procrustes: rmse 2.092428e-06  max resid 4.544947e-06 
    ## ... Similar to previous best
    ## Run 368 stress 0.1134014 
    ## Run 369 stress 0.08352634 
    ## ... Procrustes: rmse 2.360989e-06  max resid 5.391288e-06 
    ## ... Similar to previous best
    ## Run 370 stress 0.08352634 
    ## ... Procrustes: rmse 1.413084e-06  max resid 3.094587e-06 
    ## ... Similar to previous best
    ## Run 371 stress 0.08352634 
    ## ... Procrustes: rmse 6.812144e-07  max resid 1.208594e-06 
    ## ... Similar to previous best
    ## Run 372 stress 0.08352634 
    ## ... Procrustes: rmse 1.748292e-06  max resid 3.772672e-06 
    ## ... Similar to previous best
    ## Run 373 stress 0.1279309 
    ## Run 374 stress 0.08352634 
    ## ... Procrustes: rmse 3.150603e-06  max resid 6.638524e-06 
    ## ... Similar to previous best
    ## Run 375 stress 0.1279309 
    ## Run 376 stress 0.1134015 
    ## Run 377 stress 0.1134014 
    ## Run 378 stress 0.1279309 
    ## Run 379 stress 0.08352634 
    ## ... Procrustes: rmse 1.999973e-06  max resid 4.411938e-06 
    ## ... Similar to previous best
    ## Run 380 stress 0.1279309 
    ## Run 381 stress 0.08352634 
    ## ... Procrustes: rmse 1.788867e-06  max resid 3.993106e-06 
    ## ... Similar to previous best
    ## Run 382 stress 0.08352634 
    ## ... Procrustes: rmse 3.869613e-06  max resid 8.642525e-06 
    ## ... Similar to previous best
    ## Run 383 stress 0.08352634 
    ## ... Procrustes: rmse 3.914341e-06  max resid 8.841098e-06 
    ## ... Similar to previous best
    ## Run 384 stress 0.08352634 
    ## ... Procrustes: rmse 4.223839e-06  max resid 9.451957e-06 
    ## ... Similar to previous best
    ## Run 385 stress 0.08352634 
    ## ... Procrustes: rmse 1.520987e-06  max resid 2.605692e-06 
    ## ... Similar to previous best
    ## Run 386 stress 0.08352634 
    ## ... Procrustes: rmse 2.977725e-06  max resid 6.388961e-06 
    ## ... Similar to previous best
    ## Run 387 stress 0.1134014 
    ## Run 388 stress 0.08352634 
    ## ... Procrustes: rmse 1.062789e-06  max resid 2.267716e-06 
    ## ... Similar to previous best
    ## Run 389 stress 0.1989158 
    ## Run 390 stress 0.08352634 
    ## ... Procrustes: rmse 9.915305e-07  max resid 2.182455e-06 
    ## ... Similar to previous best
    ## Run 391 stress 0.1134014 
    ## Run 392 stress 0.08352634 
    ## ... Procrustes: rmse 2.321879e-06  max resid 5.209907e-06 
    ## ... Similar to previous best
    ## Run 393 stress 0.08352634 
    ## ... Procrustes: rmse 2.525316e-06  max resid 5.650197e-06 
    ## ... Similar to previous best
    ## Run 394 stress 0.1134014 
    ## Run 395 stress 0.1134014 
    ## Run 396 stress 0.1279309 
    ## Run 397 stress 0.1279309 
    ## Run 398 stress 0.1134014 
    ## Run 399 stress 0.1134015 
    ## Run 400 stress 0.08352634 
    ## ... Procrustes: rmse 5.589606e-07  max resid 9.440222e-07 
    ## ... Similar to previous best
    ## Run 401 stress 0.1279309 
    ## Run 402 stress 0.1134015 
    ## Run 403 stress 0.08352634 
    ## ... Procrustes: rmse 3.126453e-06  max resid 6.006705e-06 
    ## ... Similar to previous best
    ## Run 404 stress 0.08352634 
    ## ... Procrustes: rmse 1.120852e-06  max resid 2.24027e-06 
    ## ... Similar to previous best
    ## Run 405 stress 0.08352634 
    ## ... Procrustes: rmse 1.557406e-06  max resid 3.434461e-06 
    ## ... Similar to previous best
    ## Run 406 stress 0.08352634 
    ## ... Procrustes: rmse 1.768896e-06  max resid 3.78379e-06 
    ## ... Similar to previous best
    ## Run 407 stress 0.08352634 
    ## ... Procrustes: rmse 2.857903e-06  max resid 6.379324e-06 
    ## ... Similar to previous best
    ## Run 408 stress 0.1279309 
    ## Run 409 stress 0.08352634 
    ## ... Procrustes: rmse 1.792396e-06  max resid 3.655565e-06 
    ## ... Similar to previous best
    ## Run 410 stress 0.08352634 
    ## ... Procrustes: rmse 1.777218e-06  max resid 3.617376e-06 
    ## ... Similar to previous best
    ## Run 411 stress 0.1134014 
    ## Run 412 stress 0.08352634 
    ## ... Procrustes: rmse 3.33175e-06  max resid 6.979932e-06 
    ## ... Similar to previous best
    ## Run 413 stress 0.08352634 
    ## ... Procrustes: rmse 1.248403e-06  max resid 2.767551e-06 
    ## ... Similar to previous best
    ## Run 414 stress 0.1134014 
    ## Run 415 stress 0.340734 
    ## Run 416 stress 0.08352634 
    ## ... Procrustes: rmse 1.109102e-06  max resid 2.380238e-06 
    ## ... Similar to previous best
    ## Run 417 stress 0.08352634 
    ## ... Procrustes: rmse 3.116474e-06  max resid 6.910129e-06 
    ## ... Similar to previous best
    ## Run 418 stress 0.08352634 
    ## ... Procrustes: rmse 1.072837e-06  max resid 2.408731e-06 
    ## ... Similar to previous best
    ## Run 419 stress 0.08352634 
    ## ... Procrustes: rmse 7.571518e-07  max resid 1.495127e-06 
    ## ... Similar to previous best
    ## Run 420 stress 0.08352634 
    ## ... Procrustes: rmse 5.989372e-06  max resid 1.207619e-05 
    ## ... Similar to previous best
    ## Run 421 stress 0.08352634 
    ## ... Procrustes: rmse 1.27333e-06  max resid 2.705213e-06 
    ## ... Similar to previous best
    ## Run 422 stress 0.1279309 
    ## Run 423 stress 0.08352634 
    ## ... Procrustes: rmse 2.266883e-06  max resid 5.088331e-06 
    ## ... Similar to previous best
    ## Run 424 stress 0.08352634 
    ## ... Procrustes: rmse 1.629729e-05  max resid 3.662769e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.08352634 
    ## ... Procrustes: rmse 6.671268e-07  max resid 1.336231e-06 
    ## ... Similar to previous best
    ## Run 426 stress 0.1134014 
    ## Run 427 stress 0.08352634 
    ## ... Procrustes: rmse 2.511188e-07  max resid 4.546327e-07 
    ## ... Similar to previous best
    ## Run 428 stress 0.1134014 
    ## Run 429 stress 0.1134014 
    ## Run 430 stress 0.08352634 
    ## ... Procrustes: rmse 1.136563e-06  max resid 2.520965e-06 
    ## ... Similar to previous best
    ## Run 431 stress 0.08352634 
    ## ... Procrustes: rmse 1.340495e-06  max resid 2.271349e-06 
    ## ... Similar to previous best
    ## Run 432 stress 0.1279309 
    ## Run 433 stress 0.08352634 
    ## ... Procrustes: rmse 3.215755e-06  max resid 7.120695e-06 
    ## ... Similar to previous best
    ## Run 434 stress 0.1134015 
    ## Run 435 stress 0.08352634 
    ## ... Procrustes: rmse 1.537014e-06  max resid 3.281825e-06 
    ## ... Similar to previous best
    ## Run 436 stress 0.08352634 
    ## ... Procrustes: rmse 1.501082e-06  max resid 2.71843e-06 
    ## ... Similar to previous best
    ## Run 437 stress 0.08352634 
    ## ... Procrustes: rmse 1.277753e-06  max resid 2.806408e-06 
    ## ... Similar to previous best
    ## Run 438 stress 0.08352634 
    ## ... Procrustes: rmse 2.977986e-06  max resid 6.042407e-06 
    ## ... Similar to previous best
    ## Run 439 stress 0.08352634 
    ## ... Procrustes: rmse 4.716079e-06  max resid 1.115054e-05 
    ## ... Similar to previous best
    ## Run 440 stress 0.1134015 
    ## Run 441 stress 0.1279309 
    ## Run 442 stress 0.08352634 
    ## ... Procrustes: rmse 3.858245e-06  max resid 8.653211e-06 
    ## ... Similar to previous best
    ## Run 443 stress 0.1279309 
    ## Run 444 stress 0.1279309 
    ## Run 445 stress 0.08352634 
    ## ... Procrustes: rmse 3.286069e-06  max resid 7.220603e-06 
    ## ... Similar to previous best
    ## Run 446 stress 0.08352634 
    ## ... Procrustes: rmse 2.423702e-06  max resid 4.637519e-06 
    ## ... Similar to previous best
    ## Run 447 stress 0.08352634 
    ## ... Procrustes: rmse 2.756792e-06  max resid 6.08645e-06 
    ## ... Similar to previous best
    ## Run 448 stress 0.08352634 
    ## ... Procrustes: rmse 3.543282e-06  max resid 7.946456e-06 
    ## ... Similar to previous best
    ## Run 449 stress 0.08352634 
    ## ... Procrustes: rmse 3.97494e-06  max resid 8.806247e-06 
    ## ... Similar to previous best
    ## Run 450 stress 0.08352634 
    ## ... Procrustes: rmse 2.95031e-06  max resid 6.556259e-06 
    ## ... Similar to previous best
    ## Run 451 stress 0.3267553 
    ## Run 452 stress 0.08352634 
    ## ... Procrustes: rmse 1.762557e-06  max resid 3.378765e-06 
    ## ... Similar to previous best
    ## Run 453 stress 0.08352634 
    ## ... Procrustes: rmse 1.809416e-06  max resid 3.996955e-06 
    ## ... Similar to previous best
    ## Run 454 stress 0.1279309 
    ## Run 455 stress 0.1134015 
    ## Run 456 stress 0.08352634 
    ## ... Procrustes: rmse 2.362856e-06  max resid 5.109447e-06 
    ## ... Similar to previous best
    ## Run 457 stress 0.08352634 
    ## ... Procrustes: rmse 3.084459e-06  max resid 6.394798e-06 
    ## ... Similar to previous best
    ## Run 458 stress 0.08352634 
    ## ... Procrustes: rmse 2.614224e-06  max resid 5.789777e-06 
    ## ... Similar to previous best
    ## Run 459 stress 0.1279309 
    ## Run 460 stress 0.08352634 
    ## ... Procrustes: rmse 1.237596e-06  max resid 2.426016e-06 
    ## ... Similar to previous best
    ## Run 461 stress 0.1134014 
    ## Run 462 stress 0.08352634 
    ## ... Procrustes: rmse 2.264834e-06  max resid 5.213417e-06 
    ## ... Similar to previous best
    ## Run 463 stress 0.1279309 
    ## Run 464 stress 0.127931 
    ## Run 465 stress 0.08352634 
    ## ... Procrustes: rmse 8.28773e-07  max resid 1.556475e-06 
    ## ... Similar to previous best
    ## Run 466 stress 0.1134014 
    ## Run 467 stress 0.08352634 
    ## ... Procrustes: rmse 1.612539e-06  max resid 3.584007e-06 
    ## ... Similar to previous best
    ## Run 468 stress 0.3403702 
    ## Run 469 stress 0.1279309 
    ## Run 470 stress 0.08352634 
    ## ... Procrustes: rmse 1.92764e-06  max resid 4.414003e-06 
    ## ... Similar to previous best
    ## Run 471 stress 0.08352634 
    ## ... Procrustes: rmse 8.517395e-07  max resid 1.401044e-06 
    ## ... Similar to previous best
    ## Run 472 stress 0.08352634 
    ## ... Procrustes: rmse 4.021936e-06  max resid 8.944966e-06 
    ## ... Similar to previous best
    ## Run 473 stress 0.08352634 
    ## ... Procrustes: rmse 2.615996e-06  max resid 5.838873e-06 
    ## ... Similar to previous best
    ## Run 474 stress 0.08352634 
    ## ... Procrustes: rmse 2.248083e-06  max resid 5.01343e-06 
    ## ... Similar to previous best
    ## Run 475 stress 0.2031893 
    ## Run 476 stress 0.08352634 
    ## ... Procrustes: rmse 1.790818e-06  max resid 3.942635e-06 
    ## ... Similar to previous best
    ## Run 477 stress 0.08352634 
    ## ... Procrustes: rmse 2.233435e-06  max resid 4.828595e-06 
    ## ... Similar to previous best
    ## Run 478 stress 0.1134014 
    ## Run 479 stress 0.1134015 
    ## Run 480 stress 0.08352634 
    ## ... Procrustes: rmse 6.605028e-06  max resid 1.439183e-05 
    ## ... Similar to previous best
    ## Run 481 stress 0.08352634 
    ## ... Procrustes: rmse 1.709098e-06  max resid 3.715689e-06 
    ## ... Similar to previous best
    ## Run 482 stress 0.08352634 
    ## ... Procrustes: rmse 1.721678e-06  max resid 3.353407e-06 
    ## ... Similar to previous best
    ## Run 483 stress 0.08352634 
    ## ... Procrustes: rmse 6.028988e-07  max resid 1.036083e-06 
    ## ... Similar to previous best
    ## Run 484 stress 0.1989159 
    ## Run 485 stress 0.08352634 
    ## ... Procrustes: rmse 5.07731e-06  max resid 8.114959e-06 
    ## ... Similar to previous best
    ## Run 486 stress 0.1134015 
    ## Run 487 stress 0.08352634 
    ## ... Procrustes: rmse 2.491429e-06  max resid 5.47867e-06 
    ## ... Similar to previous best
    ## Run 488 stress 0.08352634 
    ## ... Procrustes: rmse 1.496314e-06  max resid 3.379051e-06 
    ## ... Similar to previous best
    ## Run 489 stress 0.08352634 
    ## ... Procrustes: rmse 2.877191e-07  max resid 4.057559e-07 
    ## ... Similar to previous best
    ## Run 490 stress 0.08352634 
    ## ... Procrustes: rmse 7.161301e-07  max resid 1.614635e-06 
    ## ... Similar to previous best
    ## Run 491 stress 0.1279309 
    ## Run 492 stress 0.08352634 
    ## ... Procrustes: rmse 2.612096e-06  max resid 5.704297e-06 
    ## ... Similar to previous best
    ## Run 493 stress 0.08352634 
    ## ... Procrustes: rmse 2.582864e-06  max resid 5.728095e-06 
    ## ... Similar to previous best
    ## Run 494 stress 0.127931 
    ## Run 495 stress 0.1134015 
    ## Run 496 stress 0.08352634 
    ## ... Procrustes: rmse 3.420759e-06  max resid 6.969997e-06 
    ## ... Similar to previous best
    ## Run 497 stress 0.08352634 
    ## ... Procrustes: rmse 2.657835e-06  max resid 6.161213e-06 
    ## ... Similar to previous best
    ## Run 498 stress 0.127931 
    ## Run 499 stress 0.08352634 
    ## ... Procrustes: rmse 1.870264e-06  max resid 3.491958e-06 
    ## ... Similar to previous best
    ## Run 500 stress 0.1134014 
    ## *** Best solution repeated 192 times

``` r
# Stratified lakes and ocean sites
PD_beta_env_SO_NMDS <- metaMDS(PD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.92006e-05 
    ## Run 1 stress 0.0005331308 
    ## ... Procrustes: rmse 0.004390472  max resid 0.007054716 
    ## ... Similar to previous best
    ## Run 2 stress 0.0001069238 
    ## ... Procrustes: rmse 0.0006744225  max resid 0.00135416 
    ## ... Similar to previous best
    ## Run 3 stress 0.001069164 
    ## Run 4 stress 9.774129e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002089218  max resid 0.0005997602 
    ## ... Similar to previous best
    ## Run 5 stress 0.0001585666 
    ## ... Procrustes: rmse 0.001140299  max resid 0.001920752 
    ## ... Similar to previous best
    ## Run 6 stress 0.0008103545 
    ## Run 7 stress 9.995858e-05 
    ## ... Procrustes: rmse 0.000209677  max resid 0.0006011573 
    ## ... Similar to previous best
    ## Run 8 stress 0.000220052 
    ## ... Procrustes: rmse 0.001668506  max resid 0.002730733 
    ## ... Similar to previous best
    ## Run 9 stress 9.919671e-05 
    ## ... Procrustes: rmse 0.0003133718  max resid 0.0007230929 
    ## ... Similar to previous best
    ## Run 10 stress 9.980578e-05 
    ## ... Procrustes: rmse 0.0003331654  max resid 0.000693143 
    ## ... Similar to previous best
    ## Run 11 stress 0.0007357797 
    ## Run 12 stress 9.833004e-05 
    ## ... Procrustes: rmse 0.0002517117  max resid 0.0006083248 
    ## ... Similar to previous best
    ## Run 13 stress 0.0007579538 
    ## Run 14 stress 9.885879e-05 
    ## ... Procrustes: rmse 1.054103e-05  max resid 1.837798e-05 
    ## ... Similar to previous best
    ## Run 15 stress 9.97708e-05 
    ## ... Procrustes: rmse 0.0002514165  max resid 0.0007120216 
    ## ... Similar to previous best
    ## Run 16 stress 0.001198275 
    ## Run 17 stress 9.381731e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004294431  max resid 0.0006243536 
    ## ... Similar to previous best
    ## Run 18 stress 0.000795812 
    ## Run 19 stress 0.001208218 
    ## Run 20 stress 9.99646e-05 
    ## ... Procrustes: rmse 0.00049015  max resid 0.001037702 
    ## ... Similar to previous best
    ## Run 21 stress 9.979478e-05 
    ## ... Procrustes: rmse 0.0004913463  max resid 0.00103631 
    ## ... Similar to previous best
    ## Run 22 stress 0.0008543977 
    ## Run 23 stress 9.735302e-05 
    ## ... Procrustes: rmse 0.0004657245  max resid 0.001020487 
    ## ... Similar to previous best
    ## Run 24 stress 9.838973e-05 
    ## ... Procrustes: rmse 0.0005316975  max resid 0.001121886 
    ## ... Similar to previous best
    ## Run 25 stress 9.958302e-05 
    ## ... Procrustes: rmse 0.0005388002  max resid 0.001006565 
    ## ... Similar to previous best
    ## Run 26 stress 9.734914e-05 
    ## ... Procrustes: rmse 0.0004541867  max resid 0.001021951 
    ## ... Similar to previous best
    ## Run 27 stress 9.662641e-05 
    ## ... Procrustes: rmse 0.0004651086  max resid 0.00101962 
    ## ... Similar to previous best
    ## Run 28 stress 0.000711189 
    ## Run 29 stress 0.001050591 
    ## Run 30 stress 9.624404e-05 
    ## ... Procrustes: rmse 0.0004311912  max resid 0.0006219045 
    ## ... Similar to previous best
    ## Run 31 stress 0.0007502064 
    ## Run 32 stress 9.981893e-05 
    ## ... Procrustes: rmse 0.0005386793  max resid 0.001009241 
    ## ... Similar to previous best
    ## Run 33 stress 0.3407336 
    ## Run 34 stress 0.001022864 
    ## Run 35 stress 9.883506e-05 
    ## ... Procrustes: rmse 0.0005649293  max resid 0.001029944 
    ## ... Similar to previous best
    ## Run 36 stress 0.001073211 
    ## Run 37 stress 0.00080461 
    ## Run 38 stress 9.950102e-05 
    ## ... Procrustes: rmse 0.0004246571  max resid 0.0006143481 
    ## ... Similar to previous best
    ## Run 39 stress 0.0006209792 
    ## Run 40 stress 9.977705e-05 
    ## ... Procrustes: rmse 0.0004900833  max resid 0.001037336 
    ## ... Similar to previous best
    ## Run 41 stress 9.93928e-05 
    ## ... Procrustes: rmse 0.0004239007  max resid 0.0006052418 
    ## ... Similar to previous best
    ## Run 42 stress 9.825166e-05 
    ## ... Procrustes: rmse 0.0004514945  max resid 0.001021859 
    ## ... Similar to previous best
    ## Run 43 stress 9.634195e-05 
    ## ... Procrustes: rmse 0.0004557513  max resid 0.001021667 
    ## ... Similar to previous best
    ## Run 44 stress 9.90966e-05 
    ## ... Procrustes: rmse 0.0005594327  max resid 0.001138143 
    ## ... Similar to previous best
    ## Run 45 stress 9.799276e-05 
    ## ... Procrustes: rmse 0.0005632802  max resid 0.001130909 
    ## ... Similar to previous best
    ## Run 46 stress 0.0005459613 
    ## ... Procrustes: rmse 0.004277453  max resid 0.006282052 
    ## ... Similar to previous best
    ## Run 47 stress 9.995349e-05 
    ## ... Procrustes: rmse 0.0005385176  max resid 0.001009368 
    ## ... Similar to previous best
    ## Run 48 stress 0.0006600545 
    ## Run 49 stress 0.3192024 
    ## Run 50 stress 9.673643e-05 
    ## ... Procrustes: rmse 0.0004342346  max resid 0.0006340402 
    ## ... Similar to previous best
    ## Run 51 stress 9.857702e-05 
    ## ... Procrustes: rmse 0.0005655116  max resid 0.001032156 
    ## ... Similar to previous best
    ## Run 52 stress 9.947065e-05 
    ## ... Procrustes: rmse 0.0005301891  max resid 0.001120315 
    ## ... Similar to previous best
    ## Run 53 stress 9.978883e-05 
    ## ... Procrustes: rmse 0.0005384753  max resid 0.001008963 
    ## ... Similar to previous best
    ## Run 54 stress 0.0006129965 
    ## Run 55 stress 0.2348704 
    ## Run 56 stress 9.585426e-05 
    ## ... Procrustes: rmse 0.0004356088  max resid 0.0006336455 
    ## ... Similar to previous best
    ## Run 57 stress 9.96002e-05 
    ## ... Procrustes: rmse 0.0004225862  max resid 0.0006083991 
    ## ... Similar to previous best
    ## Run 58 stress 0.0007287779 
    ## Run 59 stress 0.0007732636 
    ## Run 60 stress 0.0009188519 
    ## Run 61 stress 0.0006417615 
    ## Run 62 stress 0.0008388654 
    ## Run 63 stress 0.0009373567 
    ## Run 64 stress 9.844854e-05 
    ## ... Procrustes: rmse 0.0004948057  max resid 0.001034493 
    ## ... Similar to previous best
    ## Run 65 stress 9.77126e-05 
    ## ... Procrustes: rmse 0.0005232734  max resid 0.0009967391 
    ## ... Similar to previous best
    ## Run 66 stress 9.878169e-05 
    ## ... Procrustes: rmse 0.0005434363  max resid 0.001185941 
    ## ... Similar to previous best
    ## Run 67 stress 0.0006987909 
    ## Run 68 stress 9.887675e-05 
    ## ... Procrustes: rmse 0.000566543  max resid 0.001030311 
    ## ... Similar to previous best
    ## Run 69 stress 9.591716e-05 
    ## ... Procrustes: rmse 0.000433862  max resid 0.0006302581 
    ## ... Similar to previous best
    ## Run 70 stress 9.997701e-05 
    ## ... Procrustes: rmse 0.0004892899  max resid 0.001036208 
    ## ... Similar to previous best
    ## Run 71 stress 9.822047e-05 
    ## ... Procrustes: rmse 0.0005679261  max resid 0.001029624 
    ## ... Similar to previous best
    ## Run 72 stress 9.981748e-05 
    ## ... Procrustes: rmse 0.000489915  max resid 0.001037197 
    ## ... Similar to previous best
    ## Run 73 stress 0.0007048726 
    ## Run 74 stress 0.0009129374 
    ## Run 75 stress 9.880338e-05 
    ## ... Procrustes: rmse 0.0005311445  max resid 0.001118933 
    ## ... Similar to previous best
    ## Run 76 stress 0.0006199485 
    ## Run 77 stress 9.746795e-05 
    ## ... Procrustes: rmse 0.0004632917  max resid 0.001020062 
    ## ... Similar to previous best
    ## Run 78 stress 9.89276e-05 
    ## ... Procrustes: rmse 0.0005647997  max resid 0.001029712 
    ## ... Similar to previous best
    ## Run 79 stress 0.00106286 
    ## Run 80 stress 9.840321e-05 
    ## ... Procrustes: rmse 0.0004283997  max resid 0.0006132447 
    ## ... Similar to previous best
    ## Run 81 stress 0.0007049573 
    ## Run 82 stress 0.001082442 
    ## Run 83 stress 0.0006886495 
    ## Run 84 stress 9.723614e-05 
    ## ... Procrustes: rmse 0.0004948528  max resid 0.001032383 
    ## ... Similar to previous best
    ## Run 85 stress 0.0007442949 
    ## Run 86 stress 9.850007e-05 
    ## ... Procrustes: rmse 0.0004517908  max resid 0.001021737 
    ## ... Similar to previous best
    ## Run 87 stress 0.0005580945 
    ## ... Procrustes: rmse 0.00440754  max resid 0.006437643 
    ## ... Similar to previous best
    ## Run 88 stress 9.315522e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0005381197  max resid 0.001115771 
    ## ... Similar to previous best
    ## Run 89 stress 9.813537e-05 
    ## ... Procrustes: rmse 0.0003120156  max resid 0.0007103595 
    ## ... Similar to previous best
    ## Run 90 stress 0.0008049259 
    ## Run 91 stress 9.953407e-05 
    ## ... Procrustes: rmse 0.000202981  max resid 0.000585577 
    ## ... Similar to previous best
    ## Run 92 stress 9.717469e-05 
    ## ... Procrustes: rmse 0.0004401373  max resid 0.0007564001 
    ## ... Similar to previous best
    ## Run 93 stress 0.00106632 
    ## Run 94 stress 0.0004801882 
    ## ... Procrustes: rmse 0.003865056  max resid 0.005892554 
    ## ... Similar to previous best
    ## Run 95 stress 9.965222e-05 
    ## ... Procrustes: rmse 0.0002036085  max resid 0.0005873035 
    ## ... Similar to previous best
    ## Run 96 stress 9.946255e-05 
    ## ... Procrustes: rmse 0.0002032975  max resid 0.0005867918 
    ## ... Similar to previous best
    ## Run 97 stress 9.985789e-05 
    ## ... Procrustes: rmse 0.0003344817  max resid 0.0007240182 
    ## ... Similar to previous best
    ## Run 98 stress 0.0008877452 
    ## Run 99 stress 9.972669e-05 
    ## ... Procrustes: rmse 0.0003171516  max resid 0.0007190073 
    ## ... Similar to previous best
    ## Run 100 stress 9.892388e-05 
    ## ... Procrustes: rmse 0.0003230304  max resid 0.0008716092 
    ## ... Similar to previous best
    ## Run 101 stress 9.953711e-05 
    ## ... Procrustes: rmse 0.0002701971  max resid 0.0006929907 
    ## ... Similar to previous best
    ## Run 102 stress 9.891277e-05 
    ## ... Procrustes: rmse 0.0003230692  max resid 0.0008717582 
    ## ... Similar to previous best
    ## Run 103 stress 9.97384e-05 
    ## ... Procrustes: rmse 0.0002690853  max resid 0.0006914944 
    ## ... Similar to previous best
    ## Run 104 stress 9.865066e-05 
    ## ... Procrustes: rmse 0.0003229951  max resid 0.0008716265 
    ## ... Similar to previous best
    ## Run 105 stress 0.3275214 
    ## Run 106 stress 9.940714e-05 
    ## ... Procrustes: rmse 0.0002256787  max resid 0.0004527486 
    ## ... Similar to previous best
    ## Run 107 stress 9.950314e-05 
    ## ... Procrustes: rmse 0.000347974  max resid 0.0007081631 
    ## ... Similar to previous best
    ## Run 108 stress 9.989077e-05 
    ## ... Procrustes: rmse 0.0003253052  max resid 0.0008767412 
    ## ... Similar to previous best
    ## Run 109 stress 0.0007878386 
    ## Run 110 stress 9.957476e-05 
    ## ... Procrustes: rmse 0.0003338055  max resid 0.0007229055 
    ## ... Similar to previous best
    ## Run 111 stress 9.777302e-05 
    ## ... Procrustes: rmse 0.0002011155  max resid 0.0005807848 
    ## ... Similar to previous best
    ## Run 112 stress 0.000742696 
    ## Run 113 stress 9.948883e-05 
    ## ... Procrustes: rmse 0.0002020898  max resid 0.0004591959 
    ## ... Similar to previous best
    ## Run 114 stress 9.975001e-05 
    ## ... Procrustes: rmse 0.0003164708  max resid 0.000718408 
    ## ... Similar to previous best
    ## Run 115 stress 9.750816e-05 
    ## ... Procrustes: rmse 0.0003283864  max resid 0.0007125019 
    ## ... Similar to previous best
    ## Run 116 stress 0.001013917 
    ## Run 117 stress 0.3131388 
    ## Run 118 stress 9.939293e-05 
    ## ... Procrustes: rmse 0.0003477791  max resid 0.0007073037 
    ## ... Similar to previous best
    ## Run 119 stress 9.865372e-05 
    ## ... Procrustes: rmse 0.0002680741  max resid 0.0006872513 
    ## ... Similar to previous best
    ## Run 120 stress 9.834692e-05 
    ## ... Procrustes: rmse 0.0003132616  max resid 0.0007125299 
    ## ... Similar to previous best
    ## Run 121 stress 0.0009990513 
    ## Run 122 stress 0.0006651401 
    ## Run 123 stress 0.001160058 
    ## Run 124 stress 0.001261722 
    ## Run 125 stress 9.601473e-05 
    ## ... Procrustes: rmse 0.0003167717  max resid 0.0008563761 
    ## ... Similar to previous best
    ## Run 126 stress 9.878982e-05 
    ## ... Procrustes: rmse 0.0001455419  max resid 0.0003878154 
    ## ... Similar to previous best
    ## Run 127 stress 9.846864e-05 
    ## ... Procrustes: rmse 0.0002021674  max resid 0.0005837181 
    ## ... Similar to previous best
    ## Run 128 stress 0.0009543522 
    ## Run 129 stress 0.0009719955 
    ## Run 130 stress 9.929273e-05 
    ## ... Procrustes: rmse 0.0002031203  max resid 0.0005860914 
    ## ... Similar to previous best
    ## Run 131 stress 9.917496e-05 
    ## ... Procrustes: rmse 0.0002695933  max resid 0.0006914961 
    ## ... Similar to previous best
    ## Run 132 stress 0.000752784 
    ## Run 133 stress 9.707867e-05 
    ## ... Procrustes: rmse 0.0002002361  max resid 0.0005783151 
    ## ... Similar to previous best
    ## Run 134 stress 0.0007436068 
    ## Run 135 stress 0.0005636746 
    ## ... Procrustes: rmse 0.004585881  max resid 0.006998813 
    ## ... Similar to previous best
    ## Run 136 stress 0.000747055 
    ## Run 137 stress 0.0005858985 
    ## ... Procrustes: rmse 0.00478703  max resid 0.007319655 
    ## ... Similar to previous best
    ## Run 138 stress 9.896877e-05 
    ## ... Procrustes: rmse 0.0001478847  max resid 0.0003894684 
    ## ... Similar to previous best
    ## Run 139 stress 9.903132e-05 
    ## ... Procrustes: rmse 0.0002016379  max resid 0.0004583691 
    ## ... Similar to previous best
    ## Run 140 stress 0.0003350821 
    ## ... Procrustes: rmse 0.002643401  max resid 0.004348669 
    ## ... Similar to previous best
    ## Run 141 stress 0.001053539 
    ## Run 142 stress 0.0008982309 
    ## Run 143 stress 9.995012e-05 
    ## ... Procrustes: rmse 0.0003492166  max resid 0.0007110135 
    ## ... Similar to previous best
    ## Run 144 stress 9.98353e-05 
    ## ... Procrustes: rmse 0.0002702626  max resid 0.0006935852 
    ## ... Similar to previous best
    ## Run 145 stress 0.0006821288 
    ## Run 146 stress 9.904675e-05 
    ## ... Procrustes: rmse 0.0001466493  max resid 0.0003895064 
    ## ... Similar to previous best
    ## Run 147 stress 9.901596e-05 
    ## ... Procrustes: rmse 2.789633e-05  max resid 4.223016e-05 
    ## ... Similar to previous best
    ## Run 148 stress 9.93488e-05 
    ## ... Procrustes: rmse 3.002186e-05  max resid 4.59329e-05 
    ## ... Similar to previous best
    ## Run 149 stress 9.355857e-05 
    ## ... Procrustes: rmse 0.0001959306  max resid 0.0005648343 
    ## ... Similar to previous best
    ## Run 150 stress 0.3407334 
    ## Run 151 stress 0.0007383965 
    ## Run 152 stress 9.524822e-05 
    ## ... Procrustes: rmse 0.0001975607  max resid 0.000569949 
    ## ... Similar to previous best
    ## Run 153 stress 9.337246e-05 
    ## ... Procrustes: rmse 0.0004119087  max resid 0.0006742968 
    ## ... Similar to previous best
    ## Run 154 stress 0.0008403896 
    ## Run 155 stress 9.921869e-05 
    ## ... Procrustes: rmse 0.0002694533  max resid 0.0006908225 
    ## ... Similar to previous best
    ## Run 156 stress 0.3266727 
    ## Run 157 stress 0.0004881736 
    ## ... Procrustes: rmse 0.003934948  max resid 0.006003567 
    ## ... Similar to previous best
    ## Run 158 stress 9.998754e-05 
    ## ... Procrustes: rmse 3.292355e-05  max resid 4.878329e-05 
    ## ... Similar to previous best
    ## Run 159 stress 9.575692e-05 
    ## ... Procrustes: rmse 0.0001454428  max resid 0.0003882849 
    ## ... Similar to previous best
    ## Run 160 stress 0.0007202344 
    ## Run 161 stress 9.796645e-05 
    ## ... Procrustes: rmse 0.000146314  max resid 0.0003871636 
    ## ... Similar to previous best
    ## Run 162 stress 9.48251e-05 
    ## ... Procrustes: rmse 0.0004117328  max resid 0.000673054 
    ## ... Similar to previous best
    ## Run 163 stress 0.0006616438 
    ## Run 164 stress 0.0005795202 
    ## ... Procrustes: rmse 0.004696657  max resid 0.007158121 
    ## ... Similar to previous best
    ## Run 165 stress 0.0006904957 
    ## Run 166 stress 9.081451e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004389742  max resid 0.0007451079 
    ## ... Similar to previous best
    ## Run 167 stress 0.001237247 
    ## Run 168 stress 9.916876e-05 
    ## ... Procrustes: rmse 0.0004634719  max resid 0.0007838797 
    ## ... Similar to previous best
    ## Run 169 stress 9.953843e-05 
    ## ... Procrustes: rmse 0.00048184  max resid 0.0007095913 
    ## ... Similar to previous best
    ## Run 170 stress 0.001090757 
    ## Run 171 stress 0.0007251494 
    ## Run 172 stress 9.915016e-05 
    ## ... Procrustes: rmse 0.000425333  max resid 0.0006657268 
    ## ... Similar to previous best
    ## Run 173 stress 9.739181e-05 
    ## ... Procrustes: rmse 0.0004255054  max resid 0.0006367394 
    ## ... Similar to previous best
    ## Run 174 stress 0.001346021 
    ## Run 175 stress 9.354918e-05 
    ## ... Procrustes: rmse 0.0003817187  max resid 0.0006261519 
    ## ... Similar to previous best
    ## Run 176 stress 9.887479e-05 
    ## ... Procrustes: rmse 0.0004403977  max resid 0.0008157174 
    ## ... Similar to previous best
    ## Run 177 stress 9.432175e-05 
    ## ... Procrustes: rmse 0.0004173697  max resid 0.0006284869 
    ## ... Similar to previous best
    ## Run 178 stress 9.894222e-05 
    ## ... Procrustes: rmse 0.00048935  max resid 0.0007213498 
    ## ... Similar to previous best
    ## Run 179 stress 9.769589e-05 
    ## ... Procrustes: rmse 0.0003951481  max resid 0.0006462127 
    ## ... Similar to previous best
    ## Run 180 stress 0.0008824865 
    ## Run 181 stress 9.944842e-05 
    ## ... Procrustes: rmse 0.0004032675  max resid 0.0006556412 
    ## ... Similar to previous best
    ## Run 182 stress 9.808345e-05 
    ## ... Procrustes: rmse 0.0004378281  max resid 0.0008125756 
    ## ... Similar to previous best
    ## Run 183 stress 0.0005920101 
    ## Run 184 stress 0.0007610143 
    ## Run 185 stress 0.0005297098 
    ## ... Procrustes: rmse 0.0046465  max resid 0.007359948 
    ## ... Similar to previous best
    ## Run 186 stress 0.001195939 
    ## Run 187 stress 0.001042978 
    ## Run 188 stress 9.886098e-05 
    ## ... Procrustes: rmse 0.0003991483  max resid 0.0006443685 
    ## ... Similar to previous best
    ## Run 189 stress 0.001335439 
    ## Run 190 stress 9.910296e-05 
    ## ... Procrustes: rmse 0.0004251646  max resid 0.0006656263 
    ## ... Similar to previous best
    ## Run 191 stress 9.837371e-05 
    ## ... Procrustes: rmse 0.0004394591  max resid 0.0008130508 
    ## ... Similar to previous best
    ## Run 192 stress 0.0009347218 
    ## Run 193 stress 9.897898e-05 
    ## ... Procrustes: rmse 0.0004020441  max resid 0.000652857 
    ## ... Similar to previous best
    ## Run 194 stress 0.001024642 
    ## Run 195 stress 9.986145e-05 
    ## ... Procrustes: rmse 0.0005314042  max resid 0.0008320135 
    ## ... Similar to previous best
    ## Run 196 stress 0.001083596 
    ## Run 197 stress 9.92995e-05 
    ## ... Procrustes: rmse 0.0004356119  max resid 0.000653071 
    ## ... Similar to previous best
    ## Run 198 stress 9.804187e-05 
    ## ... Procrustes: rmse 0.0004379699  max resid 0.0008108726 
    ## ... Similar to previous best
    ## Run 199 stress 9.511812e-05 
    ## ... Procrustes: rmse 0.000180599  max resid 0.0003637186 
    ## ... Similar to previous best
    ## Run 200 stress 9.676461e-05 
    ## ... Procrustes: rmse 1.854695e-05  max resid 4.248808e-05 
    ## ... Similar to previous best
    ## Run 201 stress 0.001008098 
    ## Run 202 stress 9.928267e-05 
    ## ... Procrustes: rmse 0.0004855153  max resid 0.0007135397 
    ## ... Similar to previous best
    ## Run 203 stress 9.856584e-05 
    ## ... Procrustes: rmse 0.0004391767  max resid 0.0008126882 
    ## ... Similar to previous best
    ## Run 204 stress 0.0005987336 
    ## Run 205 stress 9.89357e-05 
    ## ... Procrustes: rmse 0.0004342175  max resid 0.0006528055 
    ## ... Similar to previous best
    ## Run 206 stress 9.573147e-05 
    ## ... Procrustes: rmse 0.0004548188  max resid 0.0008218068 
    ## ... Similar to previous best
    ## Run 207 stress 9.96343e-05 
    ## ... Procrustes: rmse 0.0004688959  max resid 0.0008411398 
    ## ... Similar to previous best
    ## Run 208 stress 9.742102e-05 
    ## ... Procrustes: rmse 0.0004650033  max resid 0.0007006617 
    ## ... Similar to previous best
    ## Run 209 stress 0.001095071 
    ## Run 210 stress 9.802512e-05 
    ## ... Procrustes: rmse 0.0004597211  max resid 0.0007760759 
    ## ... Similar to previous best
    ## Run 211 stress 0.0007122844 
    ## Run 212 stress 9.348835e-05 
    ## ... Procrustes: rmse 0.0004209615  max resid 0.0007846651 
    ## ... Similar to previous best
    ## Run 213 stress 9.718345e-05 
    ## ... Procrustes: rmse 0.0004812833  max resid 0.000712736 
    ## ... Similar to previous best
    ## Run 214 stress 9.985771e-05 
    ## ... Procrustes: rmse 0.0004670303  max resid 0.0007874801 
    ## ... Similar to previous best
    ## Run 215 stress 8.205782e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003743915  max resid 0.0005597264 
    ## ... Similar to previous best
    ## Run 216 stress 9.851799e-05 
    ## ... Procrustes: rmse 0.0001058462  max resid 0.0001605145 
    ## ... Similar to previous best
    ## Run 217 stress 0.001085005 
    ## Run 218 stress 0.001236768 
    ## Run 219 stress 9.985487e-05 
    ## ... Procrustes: rmse 0.0001914159  max resid 0.0004190382 
    ## ... Similar to previous best
    ## Run 220 stress 9.981166e-05 
    ## ... Procrustes: rmse 0.0002881368  max resid 0.0007014328 
    ## ... Similar to previous best
    ## Run 221 stress 9.857559e-05 
    ## ... Procrustes: rmse 0.0002843077  max resid 0.0006956124 
    ## ... Similar to previous best
    ## Run 222 stress 9.627243e-05 
    ## ... Procrustes: rmse 0.0001820251  max resid 0.0004082376 
    ## ... Similar to previous best
    ## Run 223 stress 9.853207e-05 
    ## ... Procrustes: rmse 0.0002193391  max resid 0.0004088841 
    ## ... Similar to previous best
    ## Run 224 stress 0.0008806301 
    ## Run 225 stress 9.740219e-05 
    ## ... Procrustes: rmse 0.0002569459  max resid 0.00047434 
    ## ... Similar to previous best
    ## Run 226 stress 0.001296737 
    ## Run 227 stress 0.0007781359 
    ## Run 228 stress 0.000531859 
    ## ... Procrustes: rmse 0.004417317  max resid 0.006639244 
    ## ... Similar to previous best
    ## Run 229 stress 9.85469e-05 
    ## ... Procrustes: rmse 0.0002610249  max resid 0.0004788753 
    ## ... Similar to previous best
    ## Run 230 stress 9.954347e-05 
    ## ... Procrustes: rmse 0.000291382  max resid 0.0006625977 
    ## ... Similar to previous best
    ## Run 231 stress 9.845314e-05 
    ## ... Procrustes: rmse 0.0001877078  max resid 0.0004148214 
    ## ... Similar to previous best
    ## Run 232 stress 9.976195e-05 
    ## ... Procrustes: rmse 0.0002498371  max resid 0.0006287809 
    ## ... Similar to previous best
    ## Run 233 stress 9.921561e-05 
    ## ... Procrustes: rmse 0.0001870454  max resid 0.0004151635 
    ## ... Similar to previous best
    ## Run 234 stress 9.942533e-05 
    ## ... Procrustes: rmse 0.0002492359  max resid 0.0006287401 
    ## ... Similar to previous best
    ## Run 235 stress 9.809575e-05 
    ## ... Procrustes: rmse 0.000175166  max resid 0.0003879147 
    ## ... Similar to previous best
    ## Run 236 stress 9.934208e-05 
    ## ... Procrustes: rmse 0.0002479772  max resid 0.0006236645 
    ## ... Similar to previous best
    ## Run 237 stress 0.001015281 
    ## Run 238 stress 9.926057e-05 
    ## ... Procrustes: rmse 0.0002483507  max resid 0.00062686 
    ## ... Similar to previous best
    ## Run 239 stress 9.818295e-05 
    ## ... Procrustes: rmse 0.0002837275  max resid 0.0006949808 
    ## ... Similar to previous best
    ## Run 240 stress 9.980474e-05 
    ## ... Procrustes: rmse 0.0002924999  max resid 0.00066545 
    ## ... Similar to previous best
    ## Run 241 stress 9.690397e-05 
    ## ... Procrustes: rmse 0.0001840377  max resid 0.0004114264 
    ## ... Similar to previous best
    ## Run 242 stress 0.0006392072 
    ## Run 243 stress 0.001246043 
    ## Run 244 stress 0.0008384203 
    ## Run 245 stress 0.0007779108 
    ## Run 246 stress 9.997495e-05 
    ## ... Procrustes: rmse 0.0001829252  max resid 0.0003973963 
    ## ... Similar to previous best
    ## Run 247 stress 0.000894039 
    ## Run 248 stress 9.77898e-05 
    ## ... Procrustes: rmse 0.0002162317  max resid 0.0004081218 
    ## ... Similar to previous best
    ## Run 249 stress 9.784545e-05 
    ## ... Procrustes: rmse 0.0001861866  max resid 0.0004126602 
    ## ... Similar to previous best
    ## Run 250 stress 9.818482e-05 
    ## ... Procrustes: rmse 0.000288195  max resid 0.0006573151 
    ## ... Similar to previous best
    ## Run 251 stress 0.001155029 
    ## Run 252 stress 9.782232e-05 
    ## ... Procrustes: rmse 0.0003388908  max resid 0.0008669076 
    ## ... Similar to previous best
    ## Run 253 stress 9.984603e-05 
    ## ... Procrustes: rmse 0.0002641373  max resid 0.0004837468 
    ## ... Similar to previous best
    ## Run 254 stress 9.774803e-05 
    ## ... Procrustes: rmse 0.0003078222  max resid 0.0006403033 
    ## ... Similar to previous best
    ## Run 255 stress 9.915154e-05 
    ## ... Procrustes: rmse 0.0001763636  max resid 0.0003848855 
    ## ... Similar to previous best
    ## Run 256 stress 0.001241703 
    ## Run 257 stress 9.957396e-05 
    ## ... Procrustes: rmse 0.0001902683  max resid 0.0004171268 
    ## ... Similar to previous best
    ## Run 258 stress 0.0009019767 
    ## Run 259 stress 9.876048e-05 
    ## ... Procrustes: rmse 0.0002615326  max resid 0.0004804555 
    ## ... Similar to previous best
    ## Run 260 stress 9.791747e-05 
    ## ... Procrustes: rmse 9.73658e-05  max resid 0.000153422 
    ## ... Similar to previous best
    ## Run 261 stress 0.0009870676 
    ## Run 262 stress 9.867745e-05 
    ## ... Procrustes: rmse 0.0002866265  max resid 0.0006533363 
    ## ... Similar to previous best
    ## Run 263 stress 9.746224e-05 
    ## ... Procrustes: rmse 0.000183934  max resid 0.0004132487 
    ## ... Similar to previous best
    ## Run 264 stress 8.756539e-05 
    ## ... Procrustes: rmse 0.0003526701  max resid 0.000631546 
    ## ... Similar to previous best
    ## Run 265 stress 9.942086e-05 
    ## ... Procrustes: rmse 0.0002618143  max resid 0.0004820196 
    ## ... Similar to previous best
    ## Run 266 stress 0.0005473417 
    ## ... Procrustes: rmse 0.004557667  max resid 0.007210523 
    ## ... Similar to previous best
    ## Run 267 stress 0.0006256317 
    ## Run 268 stress 9.98208e-05 
    ## ... Procrustes: rmse 0.0003453024  max resid 0.0008777813 
    ## ... Similar to previous best
    ## Run 269 stress 0.0007669247 
    ## Run 270 stress 9.918919e-05 
    ## ... Procrustes: rmse 0.0002890008  max resid 0.0006597687 
    ## ... Similar to previous best
    ## Run 271 stress 9.851049e-05 
    ## ... Procrustes: rmse 0.0003089302  max resid 0.0006383434 
    ## ... Similar to previous best
    ## Run 272 stress 9.930577e-05 
    ## ... Procrustes: rmse 0.0002869821  max resid 0.0007012137 
    ## ... Similar to previous best
    ## Run 273 stress 9.966553e-05 
    ## ... Procrustes: rmse 0.0002883461  max resid 0.0007022099 
    ## ... Similar to previous best
    ## Run 274 stress 9.803444e-05 
    ## ... Procrustes: rmse 0.0002449933  max resid 0.000620358 
    ## ... Similar to previous best
    ## Run 275 stress 9.755624e-05 
    ## ... Procrustes: rmse 0.0002860216  max resid 0.0006519703 
    ## ... Similar to previous best
    ## Run 276 stress 0.001019925 
    ## Run 277 stress 0.0006252627 
    ## Run 278 stress 9.982367e-05 
    ## ... Procrustes: rmse 0.0002609124  max resid 0.0004776761 
    ## ... Similar to previous best
    ## Run 279 stress 9.927493e-05 
    ## ... Procrustes: rmse 0.0002464201  max resid 0.0006220616 
    ## ... Similar to previous best
    ## Run 280 stress 9.908883e-05 
    ## ... Procrustes: rmse 0.0002592342  max resid 0.0004785011 
    ## ... Similar to previous best
    ## Run 281 stress 0.001100474 
    ## Run 282 stress 9.752233e-05 
    ## ... Procrustes: rmse 0.0001010712  max resid 0.0001601794 
    ## ... Similar to previous best
    ## Run 283 stress 0.0008255717 
    ## Run 284 stress 9.565814e-05 
    ## ... Procrustes: rmse 0.0001657912  max resid 0.0003711442 
    ## ... Similar to previous best
    ## Run 285 stress 0.0001760484 
    ## ... Procrustes: rmse 0.001323908  max resid 0.00191464 
    ## ... Similar to previous best
    ## Run 286 stress 9.887374e-05 
    ## ... Procrustes: rmse 0.0002900458  max resid 0.0006596464 
    ## ... Similar to previous best
    ## Run 287 stress 0.0006245091 
    ## Run 288 stress 9.900596e-05 
    ## ... Procrustes: rmse 0.0001807276  max resid 0.0003924442 
    ## ... Similar to previous best
    ## Run 289 stress 9.937593e-05 
    ## ... Procrustes: rmse 0.0002630119  max resid 0.0004821967 
    ## ... Similar to previous best
    ## Run 290 stress 9.9746e-05 
    ## ... Procrustes: rmse 0.000190661  max resid 0.0004170844 
    ## ... Similar to previous best
    ## Run 291 stress 0.00113702 
    ## Run 292 stress 0.0005943907 
    ## Run 293 stress 0.0004358331 
    ## ... Procrustes: rmse 0.003566253  max resid 0.005667887 
    ## ... Similar to previous best
    ## Run 294 stress 9.952825e-05 
    ## ... Procrustes: rmse 0.0002633762  max resid 0.0004828578 
    ## ... Similar to previous best
    ## Run 295 stress 9.892541e-05 
    ## ... Procrustes: rmse 0.0001120504  max resid 0.0001738569 
    ## ... Similar to previous best
    ## Run 296 stress 9.588557e-05 
    ## ... Procrustes: rmse 0.0002119758  max resid 0.000404827 
    ## ... Similar to previous best
    ## Run 297 stress 9.959411e-05 
    ## ... Procrustes: rmse 0.0002880569  max resid 0.0007008666 
    ## ... Similar to previous best
    ## Run 298 stress 0.001280845 
    ## Run 299 stress 9.900797e-05 
    ## ... Procrustes: rmse 0.0001806428  max resid 0.0003919847 
    ## ... Similar to previous best
    ## Run 300 stress 9.932326e-05 
    ## ... Procrustes: rmse 0.0003121993  max resid 0.0006463746 
    ## ... Similar to previous best
    ## Run 301 stress 0.0009150503 
    ## Run 302 stress 0.0008899947 
    ## Run 303 stress 9.809313e-05 
    ## ... Procrustes: rmse 0.0002598855  max resid 0.0004780964 
    ## ... Similar to previous best
    ## Run 304 stress 9.570206e-05 
    ## ... Procrustes: rmse 0.0003502135  max resid 0.000634015 
    ## ... Similar to previous best
    ## Run 305 stress 8.686825e-05 
    ## ... Procrustes: rmse 0.0003468557  max resid 0.0006189968 
    ## ... Similar to previous best
    ## Run 306 stress 9.977487e-05 
    ## ... Procrustes: rmse 0.0002923573  max resid 0.0006643407 
    ## ... Similar to previous best
    ## Run 307 stress 9.944872e-05 
    ## ... Procrustes: rmse 0.000311952  max resid 0.000647157 
    ## ... Similar to previous best
    ## Run 308 stress 0.0003189806 
    ## ... Procrustes: rmse 0.002527314  max resid 0.003742018 
    ## ... Similar to previous best
    ## Run 309 stress 0.0005961396 
    ## Run 310 stress 0.0005971576 
    ## Run 311 stress 9.884388e-05 
    ## ... Procrustes: rmse 9.503171e-05  max resid 0.0001930804 
    ## ... Similar to previous best
    ## Run 312 stress 9.974983e-05 
    ## ... Procrustes: rmse 0.0002190088  max resid 0.0004108514 
    ## ... Similar to previous best
    ## Run 313 stress 0.001222809 
    ## Run 314 stress 9.892463e-05 
    ## ... Procrustes: rmse 0.000289887  max resid 0.0006619713 
    ## ... Similar to previous best
    ## Run 315 stress 9.65973e-05 
    ## ... Procrustes: rmse 0.0003757546  max resid 0.0006419073 
    ## ... Similar to previous best
    ## Run 316 stress 9.9948e-05 
    ## ... Procrustes: rmse 0.000261271  max resid 0.0004808075 
    ## ... Similar to previous best
    ## Run 317 stress 9.972757e-05 
    ## ... Procrustes: rmse 0.0003463863  max resid 0.000880318 
    ## ... Similar to previous best
    ## Run 318 stress 9.868382e-05 
    ## ... Procrustes: rmse 0.0002852924  max resid 0.0006970504 
    ## ... Similar to previous best
    ## Run 319 stress 9.730475e-05 
    ## ... Procrustes: rmse 0.0001848289  max resid 0.0004111135 
    ## ... Similar to previous best
    ## Run 320 stress 9.814242e-05 
    ## ... Procrustes: rmse 0.0003088186  max resid 0.0006407202 
    ## ... Similar to previous best
    ## Run 321 stress 9.686464e-05 
    ## ... Procrustes: rmse 0.000233781  max resid 0.0006042187 
    ## ... Similar to previous best
    ## Run 322 stress 9.912429e-05 
    ## ... Procrustes: rmse 0.0002892977  max resid 0.0006600996 
    ## ... Similar to previous best
    ## Run 323 stress 9.957608e-05 
    ## ... Procrustes: rmse 0.0002438719  max resid 0.0006166698 
    ## ... Similar to previous best
    ## Run 324 stress 0.001114134 
    ## Run 325 stress 0.001069758 
    ## Run 326 stress 8.528404e-05 
    ## ... Procrustes: rmse 0.0003573349  max resid 0.0006390179 
    ## ... Similar to previous best
    ## Run 327 stress 0.001090721 
    ## Run 328 stress 9.805899e-05 
    ## ... Procrustes: rmse 0.0003496039  max resid 0.0006322143 
    ## ... Similar to previous best
    ## Run 329 stress 9.931951e-05 
    ## ... Procrustes: rmse 0.0002471882  max resid 0.0006235919 
    ## ... Similar to previous best
    ## Run 330 stress 0.001158917 
    ## Run 331 stress 9.444308e-05 
    ## ... Procrustes: rmse 0.000327604  max resid 0.0008473658 
    ## ... Similar to previous best
    ## Run 332 stress 9.84526e-05 
    ## ... Procrustes: rmse 0.0002458788  max resid 0.0006211496 
    ## ... Similar to previous best
    ## Run 333 stress 9.937904e-05 
    ## ... Procrustes: rmse 0.0002471832  max resid 0.0006218434 
    ## ... Similar to previous best
    ## Run 334 stress 0.0006734091 
    ## Run 335 stress 9.885009e-05 
    ## ... Procrustes: rmse 0.0001113718  max resid 0.000172435 
    ## ... Similar to previous best
    ## Run 336 stress 0.0007584769 
    ## Run 337 stress 0.0009742124 
    ## Run 338 stress 9.737661e-05 
    ## ... Procrustes: rmse 0.0002868548  max resid 0.000655018 
    ## ... Similar to previous best
    ## Run 339 stress 0.0006933664 
    ## Run 340 stress 9.993099e-05 
    ## ... Procrustes: rmse 0.0002494779  max resid 0.0006275251 
    ## ... Similar to previous best
    ## Run 341 stress 0.001244754 
    ## Run 342 stress 9.84578e-05 
    ## ... Procrustes: rmse 0.0002593554  max resid 0.0004758342 
    ## ... Similar to previous best
    ## Run 343 stress 0.0006558869 
    ## Run 344 stress 0.001189808 
    ## Run 345 stress 9.836321e-05 
    ## ... Procrustes: rmse 0.0001864951  max resid 0.0004132177 
    ## ... Similar to previous best
    ## Run 346 stress 9.855261e-05 
    ## ... Procrustes: rmse 0.0002197556  max resid 0.0004091273 
    ## ... Similar to previous best
    ## Run 347 stress 0.0005804785 
    ## ... Procrustes: rmse 0.004850131  max resid 0.007666094 
    ## ... Similar to previous best
    ## Run 348 stress 0.0006248273 
    ## Run 349 stress 0.001093493 
    ## Run 350 stress 0.0005588149 
    ## ... Procrustes: rmse 0.004660225  max resid 0.007024644 
    ## ... Similar to previous best
    ## Run 351 stress 9.84908e-05 
    ## ... Procrustes: rmse 0.0002192741  max resid 0.0004074202 
    ## ... Similar to previous best
    ## Run 352 stress 0.001074761 
    ## Run 353 stress 9.851989e-05 
    ## ... Procrustes: rmse 0.0001013534  max resid 0.0001564913 
    ## ... Similar to previous best
    ## Run 354 stress 9.26283e-05 
    ## ... Procrustes: rmse 0.0002008134  max resid 0.0003893344 
    ## ... Similar to previous best
    ## Run 355 stress 9.381588e-05 
    ## ... Procrustes: rmse 0.000176443  max resid 0.0004050022 
    ## ... Similar to previous best
    ## Run 356 stress 0.0007289725 
    ## Run 357 stress 9.987459e-05 
    ## ... Procrustes: rmse 0.0003449556  max resid 0.0008781075 
    ## ... Similar to previous best
    ## Run 358 stress 9.868197e-05 
    ## ... Procrustes: rmse 0.0002897871  max resid 0.0006598177 
    ## ... Similar to previous best
    ## Run 359 stress 9.618544e-05 
    ## ... Procrustes: rmse 0.0002095495  max resid 0.0003977324 
    ## ... Similar to previous best
    ## Run 360 stress 9.911502e-05 
    ## ... Procrustes: rmse 0.0002215163  max resid 0.000409638 
    ## ... Similar to previous best
    ## Run 361 stress 9.455909e-05 
    ## ... Procrustes: rmse 0.0002518392  max resid 0.0004684088 
    ## ... Similar to previous best
    ## Run 362 stress 9.977414e-05 
    ## ... Procrustes: rmse 0.0001907805  max resid 0.0004179208 
    ## ... Similar to previous best
    ## Run 363 stress 9.873592e-05 
    ## ... Procrustes: rmse 0.0002897877  max resid 0.0006602925 
    ## ... Similar to previous best
    ## Run 364 stress 9.786628e-05 
    ## ... Procrustes: rmse 0.0002591633  max resid 0.0004775168 
    ## ... Similar to previous best
    ## Run 365 stress 9.880588e-05 
    ## ... Procrustes: rmse 0.0001799145  max resid 0.0003914805 
    ## ... Similar to previous best
    ## Run 366 stress 9.866397e-05 
    ## ... Procrustes: rmse 0.000186123  max resid 0.0004140103 
    ## ... Similar to previous best
    ## Run 367 stress 9.892636e-05 
    ## ... Procrustes: rmse 0.0002595876  max resid 0.000475679 
    ## ... Similar to previous best
    ## Run 368 stress 9.676562e-05 
    ## ... Procrustes: rmse 0.0002916627  max resid 0.000621824 
    ## ... Similar to previous best
    ## Run 369 stress 0.0009159234 
    ## Run 370 stress 0.0006910815 
    ## Run 371 stress 9.830553e-05 
    ## ... Procrustes: rmse 0.0003083453  max resid 0.0006400123 
    ## ... Similar to previous best
    ## Run 372 stress 0.001228336 
    ## Run 373 stress 9.859073e-05 
    ## ... Procrustes: rmse 0.0002889237  max resid 0.0006586318 
    ## ... Similar to previous best
    ## Run 374 stress 9.854279e-05 
    ## ... Procrustes: rmse 0.0003755171  max resid 0.0006430503 
    ## ... Similar to previous best
    ## Run 375 stress 0.0008892877 
    ## Run 376 stress 9.930073e-05 
    ## ... Procrustes: rmse 0.0002216761  max resid 0.0004092636 
    ## ... Similar to previous best
    ## Run 377 stress 9.854454e-05 
    ## ... Procrustes: rmse 0.0003808226  max resid 0.0007972321 
    ## ... Similar to previous best
    ## Run 378 stress 9.980961e-05 
    ## ... Procrustes: rmse 0.0001841165  max resid 0.0003975352 
    ## ... Similar to previous best
    ## Run 379 stress 9.947139e-05 
    ## ... Procrustes: rmse 0.0002625843  max resid 0.0004812518 
    ## ... Similar to previous best
    ## Run 380 stress 0.2629798 
    ## Run 381 stress 0.3407343 
    ## Run 382 stress 0.0005615176 
    ## ... Procrustes: rmse 0.004678391  max resid 0.007063615 
    ## ... Similar to previous best
    ## Run 383 stress 0.001277547 
    ## Run 384 stress 9.91798e-05 
    ## ... Procrustes: rmse 0.0002889755  max resid 0.0006562477 
    ## ... Similar to previous best
    ## Run 385 stress 0.0006937257 
    ## Run 386 stress 9.841615e-05 
    ## ... Procrustes: rmse 0.0002846062  max resid 0.0006966508 
    ## ... Similar to previous best
    ## Run 387 stress 0.0007914374 
    ## Run 388 stress 0.0007288724 
    ## Run 389 stress 9.862411e-05 
    ## ... Procrustes: rmse 0.0002898458  max resid 0.0006603601 
    ## ... Similar to previous best
    ## Run 390 stress 9.960318e-05 
    ## ... Procrustes: rmse 0.0001143984  max resid 0.000178486 
    ## ... Similar to previous best
    ## Run 391 stress 9.877686e-05 
    ## ... Procrustes: rmse 0.0002613704  max resid 0.0004796711 
    ## ... Similar to previous best
    ## Run 392 stress 9.85391e-05 
    ## ... Procrustes: rmse 0.0003092662  max resid 0.0006428054 
    ## ... Similar to previous best
    ## Run 393 stress 9.981224e-05 
    ## ... Procrustes: rmse 0.0002920294  max resid 0.0006652693 
    ## ... Similar to previous best
    ## Run 394 stress 0.0006412299 
    ## Run 395 stress 0.0005561002 
    ## ... Procrustes: rmse 0.004619429  max resid 0.006943825 
    ## ... Similar to previous best
    ## Run 396 stress 0.0007763468 
    ## Run 397 stress 0.0007400251 
    ## Run 398 stress 0.0002769215 
    ## ... Procrustes: rmse 0.002161742  max resid 0.003518509 
    ## ... Similar to previous best
    ## Run 399 stress 0.0007383143 
    ## Run 400 stress 0.0006681157 
    ## Run 401 stress 0.0006289436 
    ## Run 402 stress 0.3266727 
    ## Run 403 stress 0.0008060724 
    ## Run 404 stress 9.780155e-05 
    ## ... Procrustes: rmse 0.0002452245  max resid 0.0006217491 
    ## ... Similar to previous best
    ## Run 405 stress 9.795373e-05 
    ## ... Procrustes: rmse 0.000259467  max resid 0.0004772736 
    ## ... Similar to previous best
    ## Run 406 stress 9.856457e-05 
    ## ... Procrustes: rmse 0.0002468623  max resid 0.0006245861 
    ## ... Similar to previous best
    ## Run 407 stress 0.0008640787 
    ## Run 408 stress 9.88192e-05 
    ## ... Procrustes: rmse 0.0002609516  max resid 0.0004784315 
    ## ... Similar to previous best
    ## Run 409 stress 0.0006326913 
    ## Run 410 stress 9.804395e-05 
    ## ... Procrustes: rmse 0.0002581067  max resid 0.0004743818 
    ## ... Similar to previous best
    ## Run 411 stress 0.001046881 
    ## Run 412 stress 9.795394e-05 
    ## ... Procrustes: rmse 0.0002588838  max resid 0.0004759032 
    ## ... Similar to previous best
    ## Run 413 stress 0.001021632 
    ## Run 414 stress 9.945137e-05 
    ## ... Procrustes: rmse 0.0004577667  max resid 0.0008816962 
    ## ... Similar to previous best
    ## Run 415 stress 0.0004265106 
    ## ... Procrustes: rmse 0.003478781  max resid 0.005529814 
    ## ... Similar to previous best
    ## Run 416 stress 9.549611e-05 
    ## ... Procrustes: rmse 0.0002822298  max resid 0.0006460416 
    ## ... Similar to previous best
    ## Run 417 stress 9.940283e-05 
    ## ... Procrustes: rmse 0.0002210826  max resid 0.0004087962 
    ## ... Similar to previous best
    ## Run 418 stress 0.0002072933 
    ## ... Procrustes: rmse 0.001544189  max resid 0.002574138 
    ## ... Similar to previous best
    ## Run 419 stress 0.000997487 
    ## Run 420 stress 9.912583e-05 
    ## ... Procrustes: rmse 0.0003114545  max resid 0.0006457957 
    ## ... Similar to previous best
    ## Run 421 stress 9.866515e-05 
    ## ... Procrustes: rmse 0.0002470712  max resid 0.0006239494 
    ## ... Similar to previous best
    ## Run 422 stress 9.93442e-05 
    ## ... Procrustes: rmse 0.0002157302  max resid 0.0004086964 
    ## ... Similar to previous best
    ## Run 423 stress 9.881455e-05 
    ## ... Procrustes: rmse 0.0002471115  max resid 0.0006255659 
    ## ... Similar to previous best
    ## Run 424 stress 9.9973e-05 
    ## ... Procrustes: rmse 0.0003136852  max resid 0.000649215 
    ## ... Similar to previous best
    ## Run 425 stress 0.0007219705 
    ## Run 426 stress 0.0007394167 
    ## Run 427 stress 9.745355e-05 
    ## ... Procrustes: rmse 0.0001833705  max resid 0.0004105137 
    ## ... Similar to previous best
    ## Run 428 stress 0.0007034475 
    ## Run 429 stress 0.0008386753 
    ## Run 430 stress 0.0003512599 
    ## ... Procrustes: rmse 0.002846969  max resid 0.004213257 
    ## ... Similar to previous best
    ## Run 431 stress 9.535469e-05 
    ## ... Procrustes: rmse 0.000364231  max resid 0.0006541907 
    ## ... Similar to previous best
    ## Run 432 stress 9.946379e-05 
    ## ... Procrustes: rmse 0.000258994  max resid 0.0004788587 
    ## ... Similar to previous best
    ## Run 433 stress 0.001200544 
    ## Run 434 stress 0.0007311285 
    ## Run 435 stress 9.885521e-05 
    ## ... Procrustes: rmse 0.0001881619  max resid 0.0004139636 
    ## ... Similar to previous best
    ## Run 436 stress 9.929267e-05 
    ## ... Procrustes: rmse 0.0002799885  max resid 0.0006889332 
    ## ... Similar to previous best
    ## Run 437 stress 0.0006591849 
    ## Run 438 stress 9.088895e-05 
    ## ... Procrustes: rmse 0.000363628  max resid 0.0005866986 
    ## ... Similar to previous best
    ## Run 439 stress 9.918707e-05 
    ## ... Procrustes: rmse 0.0001131557  max resid 0.0001759174 
    ## ... Similar to previous best
    ## Run 440 stress 9.191889e-05 
    ## ... Procrustes: rmse 0.0002451134  max resid 0.0004588715 
    ## ... Similar to previous best
    ## Run 441 stress 0.001124464 
    ## Run 442 stress 9.719128e-05 
    ## ... Procrustes: rmse 0.000356269  max resid 0.0005393957 
    ## ... Similar to previous best
    ## Run 443 stress 9.973634e-05 
    ## ... Procrustes: rmse 0.0002494335  max resid 0.000629083 
    ## ... Similar to previous best
    ## Run 444 stress 9.916581e-05 
    ## ... Procrustes: rmse 0.0003436041  max resid 0.0008751229 
    ## ... Similar to previous best
    ## Run 445 stress 9.991084e-05 
    ## ... Procrustes: rmse 0.0002892183  max resid 0.0007042922 
    ## ... Similar to previous best
    ## Run 446 stress 0.0009413395 
    ## Run 447 stress 0.001091028 
    ## Run 448 stress 0.0008757913 
    ## Run 449 stress 0.0008332521 
    ## Run 450 stress 9.886901e-05 
    ## ... Procrustes: rmse 0.000247389  max resid 0.0006240773 
    ## ... Similar to previous best
    ## Run 451 stress 9.8831e-05 
    ## ... Procrustes: rmse 0.0003418184  max resid 0.0008730246 
    ## ... Similar to previous best
    ## Run 452 stress 9.219052e-05 
    ## ... Procrustes: rmse 0.0003596846  max resid 0.0005880553 
    ## ... Similar to previous best
    ## Run 453 stress 9.992369e-05 
    ## ... Procrustes: rmse 0.0001908247  max resid 0.0004187259 
    ## ... Similar to previous best
    ## Run 454 stress 0.001228478 
    ## Run 455 stress 9.895519e-05 
    ## ... Procrustes: rmse 0.0003088809  max resid 0.0006437194 
    ## ... Similar to previous best
    ## Run 456 stress 0.0009317585 
    ## Run 457 stress 0.001012496 
    ## Run 458 stress 9.815157e-05 
    ## ... Procrustes: rmse 0.0002175118  max resid 0.000405971 
    ## ... Similar to previous best
    ## Run 459 stress 0.00108183 
    ## Run 460 stress 0.0007381887 
    ## Run 461 stress 0.0006555245 
    ## Run 462 stress 9.604431e-05 
    ## ... Procrustes: rmse 0.0002825482  max resid 0.0006492087 
    ## ... Similar to previous best
    ## Run 463 stress 9.983382e-05 
    ## ... Procrustes: rmse 0.0003466862  max resid 0.0008809492 
    ## ... Similar to previous best
    ## Run 464 stress 0.0006140382 
    ## Run 465 stress 9.97705e-05 
    ## ... Procrustes: rmse 0.0002231486  max resid 0.0004118596 
    ## ... Similar to previous best
    ## Run 466 stress 9.872621e-05 
    ## ... Procrustes: rmse 0.000559841  max resid 0.001014492 
    ## ... Similar to previous best
    ## Run 467 stress 0.0005404645 
    ## ... Procrustes: rmse 0.004455613  max resid 0.006683897 
    ## ... Similar to previous best
    ## Run 468 stress 9.938292e-05 
    ## ... Procrustes: rmse 0.0003108869  max resid 0.0006445866 
    ## ... Similar to previous best
    ## Run 469 stress 9.680202e-05 
    ## ... Procrustes: rmse 0.0001798181  max resid 0.0004101861 
    ## ... Similar to previous best
    ## Run 470 stress 9.884736e-05 
    ## ... Procrustes: rmse 0.0003108791  max resid 0.0006448598 
    ## ... Similar to previous best
    ## Run 471 stress 9.907117e-05 
    ## ... Procrustes: rmse 0.0002614387  max resid 0.0004818183 
    ## ... Similar to previous best
    ## Run 472 stress 9.981043e-05 
    ## ... Procrustes: rmse 0.0002400714  max resid 0.0006093029 
    ## ... Similar to previous best
    ## Run 473 stress 9.974293e-05 
    ## ... Procrustes: rmse 0.0001157441  max resid 0.0001788619 
    ## ... Similar to previous best
    ## Run 474 stress 0.000354874 
    ## ... Procrustes: rmse 0.00288148  max resid 0.004271906 
    ## ... Similar to previous best
    ## Run 475 stress 9.526848e-05 
    ## ... Procrustes: rmse 0.0002744497  max resid 0.0006298337 
    ## ... Similar to previous best
    ## Run 476 stress 9.823353e-05 
    ## ... Procrustes: rmse 0.0002167244  max resid 0.0004049241 
    ## ... Similar to previous best
    ## Run 477 stress 9.949608e-05 
    ## ... Procrustes: rmse 0.0001901466  max resid 0.0004157408 
    ## ... Similar to previous best
    ## Run 478 stress 9.929197e-05 
    ## ... Procrustes: rmse 0.000344319  max resid 0.0008764261 
    ## ... Similar to previous best
    ## Run 479 stress 0.0008742512 
    ## Run 480 stress 9.892829e-05 
    ## ... Procrustes: rmse 0.0002613435  max resid 0.0004782543 
    ## ... Similar to previous best
    ## Run 481 stress 9.990107e-05 
    ## ... Procrustes: rmse 0.000114294  max resid 0.0001782997 
    ## ... Similar to previous best
    ## Run 482 stress 9.958234e-05 
    ## ... Procrustes: rmse 0.0002228229  max resid 0.0004110559 
    ## ... Similar to previous best
    ## Run 483 stress 9.916002e-05 
    ## ... Procrustes: rmse 0.0001895127  max resid 0.0004172686 
    ## ... Similar to previous best
    ## Run 484 stress 9.827887e-05 
    ## ... Procrustes: rmse 0.0003411673  max resid 0.0008715717 
    ## ... Similar to previous best
    ## Run 485 stress 0.0008193841 
    ## Run 486 stress 9.930468e-05 
    ## ... Procrustes: rmse 0.000246009  max resid 0.000462467 
    ## ... Similar to previous best
    ## Run 487 stress 9.9335e-05 
    ## ... Procrustes: rmse 0.0002215822  max resid 0.000410741 
    ## ... Similar to previous best
    ## Run 488 stress 9.963787e-05 
    ## ... Procrustes: rmse 0.0006208921  max resid 0.001174403 
    ## ... Similar to previous best
    ## Run 489 stress 9.88886e-05 
    ## ... Procrustes: rmse 0.0002207918  max resid 0.0004099273 
    ## ... Similar to previous best
    ## Run 490 stress 0.001130324 
    ## Run 491 stress 9.737943e-05 
    ## ... Procrustes: rmse 0.0003068566  max resid 0.0006375554 
    ## ... Similar to previous best
    ## Run 492 stress 9.966963e-05 
    ## ... Procrustes: rmse 0.0002223306  max resid 0.0004106806 
    ## ... Similar to previous best
    ## Run 493 stress 9.978559e-05 
    ## ... Procrustes: rmse 0.0002388654  max resid 0.00061263 
    ## ... Similar to previous best
    ## Run 494 stress 0.0007543041 
    ## Run 495 stress 9.849605e-05 
    ## ... Procrustes: rmse 0.0003088093  max resid 0.0006411619 
    ## ... Similar to previous best
    ## Run 496 stress 0.0008096528 
    ## Run 497 stress 9.795001e-05 
    ## ... Procrustes: rmse 0.0002448294  max resid 0.0006179247 
    ## ... Similar to previous best
    ## Run 498 stress 9.752699e-05 
    ## ... Procrustes: rmse 0.0002585557  max resid 0.0004762635 
    ## ... Similar to previous best
    ## Run 499 stress 9.885986e-05 
    ## ... Procrustes: rmse 0.0001819542  max resid 0.0004045354 
    ## ... Similar to previous best
    ## Run 500 stress 9.792867e-05 
    ## ... Procrustes: rmse 0.0002173481  max resid 0.000406608 
    ## ... Similar to previous best
    ## *** Best solution repeated 194 times

    ## Warning in metaMDS(PD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
PD_beta_env_M_NMDS <- metaMDS(PD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.019724e-05 
    ## Run 1 stress 9.621882e-05 
    ## ... Procrustes: rmse 8.902152e-05  max resid 0.0001357669 
    ## ... Similar to previous best
    ## Run 2 stress 9.437719e-05 
    ## ... Procrustes: rmse 0.0001223353  max resid 0.0001979837 
    ## ... Similar to previous best
    ## Run 3 stress 9.400875e-05 
    ## ... Procrustes: rmse 0.0002459147  max resid 0.0003480021 
    ## ... Similar to previous best
    ## Run 4 stress 9.171644e-05 
    ## ... Procrustes: rmse 0.0002473902  max resid 0.000351768 
    ## ... Similar to previous best
    ## Run 5 stress 8.966151e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001153806  max resid 0.0002390181 
    ## ... Similar to previous best
    ## Run 6 stress 9.281535e-05 
    ## ... Procrustes: rmse 6.878908e-05  max resid 0.0001512173 
    ## ... Similar to previous best
    ## Run 7 stress 9.573857e-05 
    ## ... Procrustes: rmse 0.0002006568  max resid 0.0003598541 
    ## ... Similar to previous best
    ## Run 8 stress 9.327151e-05 
    ## ... Procrustes: rmse 0.000138378  max resid 0.0003105261 
    ## ... Similar to previous best
    ## Run 9 stress 8.727443e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001919229  max resid 0.000318249 
    ## ... Similar to previous best
    ## Run 10 stress 8.780459e-05 
    ## ... Procrustes: rmse 0.0002201873  max resid 0.0005198957 
    ## ... Similar to previous best
    ## Run 11 stress 9.242285e-05 
    ## ... Procrustes: rmse 0.000213715  max resid 0.000384886 
    ## ... Similar to previous best
    ## Run 12 stress 9.527807e-05 
    ## ... Procrustes: rmse 0.0001947399  max resid 0.00031448 
    ## ... Similar to previous best
    ## Run 13 stress 9.276917e-05 
    ## ... Procrustes: rmse 0.0001783907  max resid 0.0002978636 
    ## ... Similar to previous best
    ## Run 14 stress 9.502126e-05 
    ## ... Procrustes: rmse 0.0002401353  max resid 0.0005608943 
    ## ... Similar to previous best
    ## Run 15 stress 9.214801e-05 
    ## ... Procrustes: rmse 0.0001698631  max resid 0.0002946771 
    ## ... Similar to previous best
    ## Run 16 stress 8.401621e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001604519  max resid 0.0002457722 
    ## ... Similar to previous best
    ## Run 17 stress 8.873756e-05 
    ## ... Procrustes: rmse 0.0002000763  max resid 0.0004474317 
    ## ... Similar to previous best
    ## Run 18 stress 9.85958e-05 
    ## ... Procrustes: rmse 0.0001717618  max resid 0.0003273587 
    ## ... Similar to previous best
    ## Run 19 stress 0.3078749 
    ## Run 20 stress 9.171182e-05 
    ## ... Procrustes: rmse 0.0001324692  max resid 0.0002686797 
    ## ... Similar to previous best
    ## Run 21 stress 9.520232e-05 
    ## ... Procrustes: rmse 0.0002097732  max resid 0.0004642653 
    ## ... Similar to previous best
    ## Run 22 stress 9.146062e-05 
    ## ... Procrustes: rmse 7.449418e-05  max resid 0.0001552929 
    ## ... Similar to previous best
    ## Run 23 stress 9.675322e-05 
    ## ... Procrustes: rmse 0.0002136805  max resid 0.0004710834 
    ## ... Similar to previous best
    ## Run 24 stress 9.38626e-05 
    ## ... Procrustes: rmse 0.0001633945  max resid 0.0002441831 
    ## ... Similar to previous best
    ## Run 25 stress 8.674815e-05 
    ## ... Procrustes: rmse 0.0001497043  max resid 0.0003014395 
    ## ... Similar to previous best
    ## Run 26 stress 9.175285e-05 
    ## ... Procrustes: rmse 0.0002062488  max resid 0.0003154147 
    ## ... Similar to previous best
    ## Run 27 stress 9.315091e-05 
    ## ... Procrustes: rmse 0.0002046323  max resid 0.000304939 
    ## ... Similar to previous best
    ## Run 28 stress 9.211078e-05 
    ## ... Procrustes: rmse 0.0001586539  max resid 0.0003044364 
    ## ... Similar to previous best
    ## Run 29 stress 9.969045e-05 
    ## ... Procrustes: rmse 0.0002210955  max resid 0.000352569 
    ## ... Similar to previous best
    ## Run 30 stress 8.921816e-05 
    ## ... Procrustes: rmse 0.0001491257  max resid 0.0002824631 
    ## ... Similar to previous best
    ## Run 31 stress 9.174181e-05 
    ## ... Procrustes: rmse 0.0001854526  max resid 0.0003146319 
    ## ... Similar to previous best
    ## Run 32 stress 9.79395e-05 
    ## ... Procrustes: rmse 0.0002151139  max resid 0.000473972 
    ## ... Similar to previous best
    ## Run 33 stress 9.434801e-05 
    ## ... Procrustes: rmse 0.0001659143  max resid 0.0003200763 
    ## ... Similar to previous best
    ## Run 34 stress 9.852367e-05 
    ## ... Procrustes: rmse 0.0001203236  max resid 0.0002243688 
    ## ... Similar to previous best
    ## Run 35 stress 9.194528e-05 
    ## ... Procrustes: rmse 0.0002111203  max resid 0.000376231 
    ## ... Similar to previous best
    ## Run 36 stress 9.726318e-05 
    ## ... Procrustes: rmse 0.0001738798  max resid 0.0003369352 
    ## ... Similar to previous best
    ## Run 37 stress 9.485902e-05 
    ## ... Procrustes: rmse 0.0001667292  max resid 0.0003333217 
    ## ... Similar to previous best
    ## Run 38 stress 9.522106e-05 
    ## ... Procrustes: rmse 0.0001672596  max resid 0.0003227558 
    ## ... Similar to previous best
    ## Run 39 stress 8.644418e-05 
    ## ... Procrustes: rmse 0.0001547367  max resid 0.0002388539 
    ## ... Similar to previous best
    ## Run 40 stress 9.647045e-05 
    ## ... Procrustes: rmse 0.0001714984  max resid 0.0003323428 
    ## ... Similar to previous best
    ## Run 41 stress 9.377638e-05 
    ## ... Procrustes: rmse 0.0002083855  max resid 0.0003255235 
    ## ... Similar to previous best
    ## Run 42 stress 9.372567e-05 
    ## ... Procrustes: rmse 0.0001565805  max resid 0.0002735593 
    ## ... Similar to previous best
    ## Run 43 stress 9.052992e-05 
    ## ... Procrustes: rmse 0.0001472919  max resid 0.0002520495 
    ## ... Similar to previous best
    ## Run 44 stress 9.210934e-05 
    ## ... Procrustes: rmse 0.0002072122  max resid 0.0003760878 
    ## ... Similar to previous best
    ## Run 45 stress 9.839935e-05 
    ## ... Procrustes: rmse 0.0002092195  max resid 0.0003322726 
    ## ... Similar to previous best
    ## Run 46 stress 9.693247e-05 
    ## ... Procrustes: rmse 0.0002108566  max resid 0.0003641109 
    ## ... Similar to previous best
    ## Run 47 stress 9.850725e-05 
    ## ... Procrustes: rmse 0.0001541236  max resid 0.0002377673 
    ## ... Similar to previous best
    ## Run 48 stress 9.734421e-05 
    ## ... Procrustes: rmse 0.0001713173  max resid 0.000341248 
    ## ... Similar to previous best
    ## Run 49 stress 9.896767e-05 
    ## ... Procrustes: rmse 0.0002149983  max resid 0.0004737702 
    ## ... Similar to previous best
    ## Run 50 stress 8.846424e-05 
    ## ... Procrustes: rmse 7.155741e-05  max resid 0.0001549283 
    ## ... Similar to previous best
    ## Run 51 stress 9.450222e-05 
    ## ... Procrustes: rmse 0.0001124364  max resid 0.0002626242 
    ## ... Similar to previous best
    ## Run 52 stress 8.554214e-05 
    ## ... Procrustes: rmse 0.0001933405  max resid 0.0004409896 
    ## ... Similar to previous best
    ## Run 53 stress 9.499212e-05 
    ## ... Procrustes: rmse 0.0002135473  max resid 0.0003719871 
    ## ... Similar to previous best
    ## Run 54 stress 9.916184e-05 
    ## ... Procrustes: rmse 0.0002198393  max resid 0.0003823375 
    ## ... Similar to previous best
    ## Run 55 stress 9.990364e-05 
    ## ... Procrustes: rmse 0.000220467  max resid 0.0003483639 
    ## ... Similar to previous best
    ## Run 56 stress 9.063679e-05 
    ## ... Procrustes: rmse 0.0001917144  max resid 0.0004381208 
    ## ... Similar to previous best
    ## Run 57 stress 9.069813e-05 
    ## ... Procrustes: rmse 0.0001508093  max resid 0.0002889722 
    ## ... Similar to previous best
    ## Run 58 stress 9.503571e-05 
    ## ... Procrustes: rmse 0.0001690162  max resid 0.0003268712 
    ## ... Similar to previous best
    ## Run 59 stress 9.397832e-05 
    ## ... Procrustes: rmse 0.0001637721  max resid 0.0003124948 
    ## ... Similar to previous best
    ## Run 60 stress 8.533383e-05 
    ## ... Procrustes: rmse 0.0001478045  max resid 0.000295653 
    ## ... Similar to previous best
    ## Run 61 stress 9.967077e-05 
    ## ... Procrustes: rmse 0.0002163616  max resid 0.0004760214 
    ## ... Similar to previous best
    ## Run 62 stress 9.411916e-05 
    ## ... Procrustes: rmse 0.0002049695  max resid 0.0004545488 
    ## ... Similar to previous best
    ## Run 63 stress 8.855721e-05 
    ## ... Procrustes: rmse 0.0001816806  max resid 0.0003232615 
    ## ... Similar to previous best
    ## Run 64 stress 0.2319215 
    ## Run 65 stress 9.916221e-05 
    ## ... Procrustes: rmse 0.0002209361  max resid 0.0003405106 
    ## ... Similar to previous best
    ## Run 66 stress 9.472346e-05 
    ## ... Procrustes: rmse 0.0002124811  max resid 0.000325873 
    ## ... Similar to previous best
    ## Run 67 stress 9.35446e-05 
    ## ... Procrustes: rmse 0.000207366  max resid 0.000459909 
    ## ... Similar to previous best
    ## Run 68 stress 9.85964e-05 
    ## ... Procrustes: rmse 0.0002207616  max resid 0.0003414751 
    ## ... Similar to previous best
    ## Run 69 stress 9.313932e-05 
    ## ... Procrustes: rmse 0.0002069876  max resid 0.0003237965 
    ## ... Similar to previous best
    ## Run 70 stress 9.639729e-05 
    ## ... Procrustes: rmse 0.0001728626  max resid 0.0003390456 
    ## ... Similar to previous best
    ## Run 71 stress 9.822447e-05 
    ## ... Procrustes: rmse 0.0002161  max resid 0.0003297597 
    ## ... Similar to previous best
    ## Run 72 stress 9.084688e-05 
    ## ... Procrustes: rmse 0.000205203  max resid 0.0004558801 
    ## ... Similar to previous best
    ## Run 73 stress 9.689928e-05 
    ## ... Procrustes: rmse 0.0002205187  max resid 0.0003125241 
    ## ... Similar to previous best
    ## Run 74 stress 9.674663e-05 
    ## ... Procrustes: rmse 0.0001719706  max resid 0.0003361413 
    ## ... Similar to previous best
    ## Run 75 stress 9.560838e-05 
    ## ... Procrustes: rmse 0.0002120582  max resid 0.0003203761 
    ## ... Similar to previous best
    ## Run 76 stress 9.171876e-05 
    ## ... Procrustes: rmse 0.000197644  max resid 0.0004453216 
    ## ... Similar to previous best
    ## Run 77 stress 9.966432e-05 
    ## ... Procrustes: rmse 0.000200798  max resid 0.0002911593 
    ## ... Similar to previous best
    ## Run 78 stress 9.909838e-05 
    ## ... Procrustes: rmse 0.0002186803  max resid 0.0003338411 
    ## ... Similar to previous best
    ## Run 79 stress 9.61451e-05 
    ## ... Procrustes: rmse 0.0001634895  max resid 0.0003216715 
    ## ... Similar to previous best
    ## Run 80 stress 9.461403e-05 
    ## ... Procrustes: rmse 0.0002040802  max resid 0.0003192168 
    ## ... Similar to previous best
    ## Run 81 stress 9.869202e-05 
    ## ... Procrustes: rmse 0.0001677335  max resid 0.0002516223 
    ## ... Similar to previous best
    ## Run 82 stress 9.192401e-05 
    ## ... Procrustes: rmse 0.000113645  max resid 0.0001739061 
    ## ... Similar to previous best
    ## Run 83 stress 9.778971e-05 
    ## ... Procrustes: rmse 0.0001753332  max resid 0.0003406585 
    ## ... Similar to previous best
    ## Run 84 stress 9.006905e-05 
    ## ... Procrustes: rmse 0.0001478732  max resid 0.0002268371 
    ## ... Similar to previous best
    ## Run 85 stress 8.600845e-05 
    ## ... Procrustes: rmse 0.0001506912  max resid 0.000300427 
    ## ... Similar to previous best
    ## Run 86 stress 9.421851e-05 
    ## ... Procrustes: rmse 0.0001605433  max resid 0.0003061946 
    ## ... Similar to previous best
    ## Run 87 stress 8.60113e-05 
    ## ... Procrustes: rmse 0.0001499992  max resid 0.0002948578 
    ## ... Similar to previous best
    ## Run 88 stress 9.525522e-05 
    ## ... Procrustes: rmse 0.0001499892  max resid 0.0002646135 
    ## ... Similar to previous best
    ## Run 89 stress 9.766926e-05 
    ## ... Procrustes: rmse 0.0001739936  max resid 0.0003370611 
    ## ... Similar to previous best
    ## Run 90 stress 9.22833e-05 
    ## ... Procrustes: rmse 0.0001943029  max resid 0.0003215449 
    ## ... Similar to previous best
    ## Run 91 stress 9.541423e-05 
    ## ... Procrustes: rmse 0.0001644602  max resid 0.000314474 
    ## ... Similar to previous best
    ## Run 92 stress 9.837692e-05 
    ## ... Procrustes: rmse 0.000174956  max resid 0.0003382995 
    ## ... Similar to previous best
    ## Run 93 stress 9.314879e-05 
    ## ... Procrustes: rmse 0.0001586408  max resid 0.0003153591 
    ## ... Similar to previous best
    ## Run 94 stress 9.809996e-05 
    ## ... Procrustes: rmse 0.0001738849  max resid 0.0003450481 
    ## ... Similar to previous best
    ## Run 95 stress 9.695529e-05 
    ## ... Procrustes: rmse 0.0002121243  max resid 0.0004684741 
    ## ... Similar to previous best
    ## Run 96 stress 9.408184e-05 
    ## ... Procrustes: rmse 0.0001672891  max resid 0.0003273298 
    ## ... Similar to previous best
    ## Run 97 stress 9.880425e-05 
    ## ... Procrustes: rmse 0.0001568263  max resid 0.0002901183 
    ## ... Similar to previous best
    ## Run 98 stress 9.203588e-05 
    ## ... Procrustes: rmse 0.0001892936  max resid 0.0002777764 
    ## ... Similar to previous best
    ## Run 99 stress 9.299402e-05 
    ## ... Procrustes: rmse 0.0002053893  max resid 0.0004556325 
    ## ... Similar to previous best
    ## Run 100 stress 9.777262e-05 
    ## ... Procrustes: rmse 0.000214032  max resid 0.0003380849 
    ## ... Similar to previous best
    ## Run 101 stress 9.756564e-05 
    ## ... Procrustes: rmse 0.0001913037  max resid 0.0003164297 
    ## ... Similar to previous best
    ## Run 102 stress 9.353448e-05 
    ## ... Procrustes: rmse 0.0001115003  max resid 0.0001705952 
    ## ... Similar to previous best
    ## Run 103 stress 9.929958e-05 
    ## ... Procrustes: rmse 0.0001777205  max resid 0.0003479524 
    ## ... Similar to previous best
    ## Run 104 stress 9.323961e-05 
    ## ... Procrustes: rmse 0.0001576296  max resid 0.000301657 
    ## ... Similar to previous best
    ## Run 105 stress 9.530293e-05 
    ## ... Procrustes: rmse 0.0001072225  max resid 0.0001664918 
    ## ... Similar to previous best
    ## Run 106 stress 8.931821e-05 
    ## ... Procrustes: rmse 0.0001761687  max resid 0.0002471125 
    ## ... Similar to previous best
    ## Run 107 stress 9.398632e-05 
    ## ... Procrustes: rmse 0.0001020808  max resid 0.0001536957 
    ## ... Similar to previous best
    ## Run 108 stress 8.975632e-05 
    ## ... Procrustes: rmse 0.0001009492  max resid 0.0001506656 
    ## ... Similar to previous best
    ## Run 109 stress 9.396827e-05 
    ## ... Procrustes: rmse 0.000115277  max resid 0.0001752227 
    ## ... Similar to previous best
    ## Run 110 stress 9.281698e-05 
    ## ... Procrustes: rmse 0.0001641233  max resid 0.0003195619 
    ## ... Similar to previous best
    ## Run 111 stress 8.565928e-05 
    ## ... Procrustes: rmse 0.0001881312  max resid 0.0002740771 
    ## ... Similar to previous best
    ## Run 112 stress 8.676967e-05 
    ## ... Procrustes: rmse 0.0001487995  max resid 0.0002910067 
    ## ... Similar to previous best
    ## Run 113 stress 9.461636e-05 
    ## ... Procrustes: rmse 0.0001684435  max resid 0.0003274978 
    ## ... Similar to previous best
    ## Run 114 stress 8.775589e-05 
    ## ... Procrustes: rmse 0.0001525728  max resid 0.0002378277 
    ## ... Similar to previous best
    ## Run 115 stress 9.110016e-05 
    ## ... Procrustes: rmse 0.0001578425  max resid 0.000311935 
    ## ... Similar to previous best
    ## Run 116 stress 9.729175e-05 
    ## ... Procrustes: rmse 0.0001152506  max resid 0.0002175208 
    ## ... Similar to previous best
    ## Run 117 stress 0.3083098 
    ## Run 118 stress 9.864493e-05 
    ## ... Procrustes: rmse 0.0002003067  max resid 0.0004435313 
    ## ... Similar to previous best
    ## Run 119 stress 9.470903e-05 
    ## ... Procrustes: rmse 0.0002093698  max resid 0.0004690306 
    ## ... Similar to previous best
    ## Run 120 stress 9.751367e-05 
    ## ... Procrustes: rmse 0.0002115976  max resid 0.0004667107 
    ## ... Similar to previous best
    ## Run 121 stress 9.923105e-05 
    ## ... Procrustes: rmse 0.0001765397  max resid 0.0003460555 
    ## ... Similar to previous best
    ## Run 122 stress 8.858569e-05 
    ## ... Procrustes: rmse 0.0002008704  max resid 0.0004524039 
    ## ... Similar to previous best
    ## Run 123 stress 8.988244e-05 
    ## ... Procrustes: rmse 0.000154662  max resid 0.0002995643 
    ## ... Similar to previous best
    ## Run 124 stress 9.157104e-05 
    ## ... Procrustes: rmse 0.0001480406  max resid 0.0002123598 
    ## ... Similar to previous best
    ## Run 125 stress 9.510135e-05 
    ## ... Procrustes: rmse 0.000111507  max resid 0.0001709775 
    ## ... Similar to previous best
    ## Run 126 stress 8.861204e-05 
    ## ... Procrustes: rmse 0.0001485159  max resid 0.0002863518 
    ## ... Similar to previous best
    ## Run 127 stress 9.684418e-05 
    ## ... Procrustes: rmse 0.0002152555  max resid 0.0003386571 
    ## ... Similar to previous best
    ## Run 128 stress 9.425402e-05 
    ## ... Procrustes: rmse 0.0002096356  max resid 0.0003271181 
    ## ... Similar to previous best
    ## Run 129 stress 6.319104e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001123746  max resid 0.0001779557 
    ## ... Similar to previous best
    ## Run 130 stress 7.494647e-05 
    ## ... Procrustes: rmse 0.0001514591  max resid 0.0003039078 
    ## ... Similar to previous best
    ## Run 131 stress 9.786403e-05 
    ## ... Procrustes: rmse 0.0001795095  max resid 0.0003169236 
    ## ... Similar to previous best
    ## Run 132 stress 9.85563e-05 
    ## ... Procrustes: rmse 0.0001804388  max resid 0.0003189048 
    ## ... Similar to previous best
    ## Run 133 stress 8.786877e-05 
    ## ... Procrustes: rmse 7.300143e-05  max resid 0.0001405668 
    ## ... Similar to previous best
    ## Run 134 stress 9.830654e-05 
    ## ... Procrustes: rmse 0.0001727938  max resid 0.0003165252 
    ## ... Similar to previous best
    ## Run 135 stress 9.440056e-05 
    ## ... Procrustes: rmse 0.000238848  max resid 0.0003888419 
    ## ... Similar to previous best
    ## Run 136 stress 9.818524e-05 
    ## ... Procrustes: rmse 0.0001801584  max resid 0.0003181998 
    ## ... Similar to previous best
    ## Run 137 stress 0.2842783 
    ## Run 138 stress 9.404356e-05 
    ## ... Procrustes: rmse 0.0002353381  max resid 0.0003865432 
    ## ... Similar to previous best
    ## Run 139 stress 9.611181e-05 
    ## ... Procrustes: rmse 0.0002105945  max resid 0.0003939395 
    ## ... Similar to previous best
    ## Run 140 stress 9.342841e-05 
    ## ... Procrustes: rmse 0.0001717406  max resid 0.0003011283 
    ## ... Similar to previous best
    ## Run 141 stress 9.167638e-05 
    ## ... Procrustes: rmse 0.0001719574  max resid 0.0003021415 
    ## ... Similar to previous best
    ## Run 142 stress 9.210485e-05 
    ## ... Procrustes: rmse 0.0002298734  max resid 0.0003824551 
    ## ... Similar to previous best
    ## Run 143 stress 0.2842783 
    ## Run 144 stress 9.982486e-05 
    ## ... Procrustes: rmse 0.0002474567  max resid 0.0004085839 
    ## ... Similar to previous best
    ## Run 145 stress 9.372667e-05 
    ## ... Procrustes: rmse 0.0002049082  max resid 0.0003887182 
    ## ... Similar to previous best
    ## Run 146 stress 9.267107e-05 
    ## ... Procrustes: rmse 0.0001024262  max resid 0.0002163748 
    ## ... Similar to previous best
    ## Run 147 stress 8.465169e-05 
    ## ... Procrustes: rmse 0.0001716282  max resid 0.00026964 
    ## ... Similar to previous best
    ## Run 148 stress 9.625998e-05 
    ## ... Procrustes: rmse 0.000108579  max resid 0.0001628817 
    ## ... Similar to previous best
    ## Run 149 stress 9.965176e-05 
    ## ... Procrustes: rmse 0.0002461969  max resid 0.0004050756 
    ## ... Similar to previous best
    ## Run 150 stress 0.3080295 
    ## Run 151 stress 9.80123e-05 
    ## ... Procrustes: rmse 0.0002066724  max resid 0.0003780259 
    ## ... Similar to previous best
    ## Run 152 stress 9.660091e-05 
    ## ... Procrustes: rmse 8.491447e-05  max resid 0.0001687675 
    ## ... Similar to previous best
    ## Run 153 stress 9.701223e-05 
    ## ... Procrustes: rmse 0.0002115374  max resid 0.0003920629 
    ## ... Similar to previous best
    ## Run 154 stress 9.902915e-05 
    ## ... Procrustes: rmse 0.0001124348  max resid 0.0001531596 
    ## ... Similar to previous best
    ## Run 155 stress 9.87131e-05 
    ## ... Procrustes: rmse 0.000177739  max resid 0.0003200913 
    ## ... Similar to previous best
    ## Run 156 stress 8.819289e-05 
    ## ... Procrustes: rmse 7.247734e-05  max resid 0.0001395003 
    ## ... Similar to previous best
    ## Run 157 stress 9.698917e-05 
    ## ... Procrustes: rmse 0.0002090949  max resid 0.0003911305 
    ## ... Similar to previous best
    ## Run 158 stress 9.31234e-05 
    ## ... Procrustes: rmse 0.0001970953  max resid 0.0003711191 
    ## ... Similar to previous best
    ## Run 159 stress 9.545686e-05 
    ## ... Procrustes: rmse 0.0002354089  max resid 0.0003694878 
    ## ... Similar to previous best
    ## Run 160 stress 9.29712e-05 
    ## ... Procrustes: rmse 0.0001670363  max resid 0.0003059249 
    ## ... Similar to previous best
    ## Run 161 stress 0.2848382 
    ## Run 162 stress 9.205504e-05 
    ## ... Procrustes: rmse 0.0001669555  max resid 0.0003259259 
    ## ... Similar to previous best
    ## Run 163 stress 9.347614e-05 
    ## ... Procrustes: rmse 8.808117e-05  max resid 0.0002003889 
    ## ... Similar to previous best
    ## Run 164 stress 9.891666e-05 
    ## ... Procrustes: rmse 0.0001390021  max resid 0.0002272262 
    ## ... Similar to previous best
    ## Run 165 stress 9.804432e-05 
    ## ... Procrustes: rmse 0.0002064588  max resid 0.0003863667 
    ## ... Similar to previous best
    ## Run 166 stress 9.738228e-05 
    ## ... Procrustes: rmse 0.0002394243  max resid 0.0003757731 
    ## ... Similar to previous best
    ## Run 167 stress 9.755806e-05 
    ## ... Procrustes: rmse 0.0002466402  max resid 0.0003932934 
    ## ... Similar to previous best
    ## Run 168 stress 9.53702e-05 
    ## ... Procrustes: rmse 0.0002079857  max resid 0.0003934759 
    ## ... Similar to previous best
    ## Run 169 stress 9.954211e-05 
    ## ... Procrustes: rmse 0.0002496061  max resid 0.0003925772 
    ## ... Similar to previous best
    ## Run 170 stress 8.775891e-05 
    ## ... Procrustes: rmse 0.0001567497  max resid 0.0003141242 
    ## ... Similar to previous best
    ## Run 171 stress 8.60028e-05 
    ## ... Procrustes: rmse 6.78997e-05  max resid 0.0001287143 
    ## ... Similar to previous best
    ## Run 172 stress 8.069887e-05 
    ## ... Procrustes: rmse 7.124715e-05  max resid 0.0001466104 
    ## ... Similar to previous best
    ## Run 173 stress 9.108598e-05 
    ## ... Procrustes: rmse 0.0001609994  max resid 0.0003203836 
    ## ... Similar to previous best
    ## Run 174 stress 9.533452e-05 
    ## ... Procrustes: rmse 0.0001053009  max resid 0.0002220042 
    ## ... Similar to previous best
    ## Run 175 stress 8.646106e-05 
    ## ... Procrustes: rmse 9.211547e-05  max resid 0.0001938588 
    ## ... Similar to previous best
    ## Run 176 stress 8.63088e-05 
    ## ... Procrustes: rmse 0.0002207088  max resid 0.000359337 
    ## ... Similar to previous best
    ## Run 177 stress 8.997135e-05 
    ## ... Procrustes: rmse 0.0002273644  max resid 0.0003606556 
    ## ... Similar to previous best
    ## Run 178 stress 6.811067e-05 
    ## ... Procrustes: rmse 4.408283e-05  max resid 6.028766e-05 
    ## ... Similar to previous best
    ## Run 179 stress 8.685332e-05 
    ## ... Procrustes: rmse 0.0001619931  max resid 0.0002849177 
    ## ... Similar to previous best
    ## Run 180 stress 9.042313e-05 
    ## ... Procrustes: rmse 0.0001823636  max resid 0.0003559238 
    ## ... Similar to previous best
    ## Run 181 stress 9.385478e-05 
    ## ... Procrustes: rmse 0.0001096592  max resid 0.0001952782 
    ## ... Similar to previous best
    ## Run 182 stress 9.804826e-05 
    ## ... Procrustes: rmse 0.0002466162  max resid 0.0003920867 
    ## ... Similar to previous best
    ## Run 183 stress 9.782902e-05 
    ## ... Procrustes: rmse 0.0002420885  max resid 0.0003830208 
    ## ... Similar to previous best
    ## Run 184 stress 8.734918e-05 
    ## ... Procrustes: rmse 7.698905e-05  max resid 0.0001701196 
    ## ... Similar to previous best
    ## Run 185 stress 9.500588e-05 
    ## ... Procrustes: rmse 8.002845e-05  max resid 0.0001571257 
    ## ... Similar to previous best
    ## Run 186 stress 7.784249e-05 
    ## ... Procrustes: rmse 0.0001935928  max resid 0.000316064 
    ## ... Similar to previous best
    ## Run 187 stress 9.631507e-05 
    ## ... Procrustes: rmse 0.0002439026  max resid 0.0004000337 
    ## ... Similar to previous best
    ## Run 188 stress 9.910752e-05 
    ## ... Procrustes: rmse 0.0002363369  max resid 0.0003940867 
    ## ... Similar to previous best
    ## Run 189 stress 9.076674e-05 
    ## ... Procrustes: rmse 9.480381e-05  max resid 0.0001964309 
    ## ... Similar to previous best
    ## Run 190 stress 9.748139e-05 
    ## ... Procrustes: rmse 0.0001789528  max resid 0.0003157541 
    ## ... Similar to previous best
    ## Run 191 stress 0.3120129 
    ## Run 192 stress 9.093954e-05 
    ## ... Procrustes: rmse 0.0001596823  max resid 0.0002982031 
    ## ... Similar to previous best
    ## Run 193 stress 9.552356e-05 
    ## ... Procrustes: rmse 0.0001793343  max resid 0.0003087554 
    ## ... Similar to previous best
    ## Run 194 stress 8.895749e-05 
    ## ... Procrustes: rmse 0.0001609527  max resid 0.0003152061 
    ## ... Similar to previous best
    ## Run 195 stress 9.311929e-05 
    ## ... Procrustes: rmse 0.0002332953  max resid 0.0003727792 
    ## ... Similar to previous best
    ## Run 196 stress 9.628907e-05 
    ## ... Procrustes: rmse 0.0001696016  max resid 0.0003372653 
    ## ... Similar to previous best
    ## Run 197 stress 8.699072e-05 
    ## ... Procrustes: rmse 0.0002165986  max resid 0.000355767 
    ## ... Similar to previous best
    ## Run 198 stress 9.919887e-05 
    ## ... Procrustes: rmse 0.0002520341  max resid 0.0004068626 
    ## ... Similar to previous best
    ## Run 199 stress 9.448843e-05 
    ## ... Procrustes: rmse 0.00023534  max resid 0.000374914 
    ## ... Similar to previous best
    ## Run 200 stress 9.305727e-05 
    ## ... Procrustes: rmse 0.0002010677  max resid 0.0003217098 
    ## ... Similar to previous best
    ## Run 201 stress 9.836717e-05 
    ## ... Procrustes: rmse 0.0002485995  max resid 0.000398569 
    ## ... Similar to previous best
    ## Run 202 stress 9.521171e-05 
    ## ... Procrustes: rmse 0.000241047  max resid 0.0003931999 
    ## ... Similar to previous best
    ## Run 203 stress 9.624691e-05 
    ## ... Procrustes: rmse 0.0002246837  max resid 0.0003882228 
    ## ... Similar to previous best
    ## Run 204 stress 7.17525e-05 
    ## ... Procrustes: rmse 8.426705e-05  max resid 0.0001482035 
    ## ... Similar to previous best
    ## Run 205 stress 9.169341e-05 
    ## ... Procrustes: rmse 0.0001563854  max resid 0.0003005237 
    ## ... Similar to previous best
    ## Run 206 stress 9.816972e-05 
    ## ... Procrustes: rmse 0.0002299495  max resid 0.0003536528 
    ## ... Similar to previous best
    ## Run 207 stress 9.366977e-05 
    ## ... Procrustes: rmse 0.0001073203  max resid 0.0001625433 
    ## ... Similar to previous best
    ## Run 208 stress 8.403076e-05 
    ## ... Procrustes: rmse 0.0001336169  max resid 0.0002021267 
    ## ... Similar to previous best
    ## Run 209 stress 9.761477e-05 
    ## ... Procrustes: rmse 0.0001082937  max resid 0.0001622597 
    ## ... Similar to previous best
    ## Run 210 stress 9.96895e-05 
    ## ... Procrustes: rmse 0.0002200709  max resid 0.0003374695 
    ## ... Similar to previous best
    ## Run 211 stress 9.497138e-05 
    ## ... Procrustes: rmse 0.0001772715  max resid 0.0003077211 
    ## ... Similar to previous best
    ## Run 212 stress 0.2852155 
    ## Run 213 stress 9.215224e-05 
    ## ... Procrustes: rmse 0.0002010625  max resid 0.0003787268 
    ## ... Similar to previous best
    ## Run 214 stress 9.581168e-05 
    ## ... Procrustes: rmse 0.0001773683  max resid 0.0003093316 
    ## ... Similar to previous best
    ## Run 215 stress 7.773226e-05 
    ## ... Procrustes: rmse 0.0001461604  max resid 0.0002591159 
    ## ... Similar to previous best
    ## Run 216 stress 9.09725e-05 
    ## ... Procrustes: rmse 0.0001981785  max resid 0.000373913 
    ## ... Similar to previous best
    ## Run 217 stress 9.246776e-05 
    ## ... Procrustes: rmse 0.0001517826  max resid 0.0002996354 
    ## ... Similar to previous best
    ## Run 218 stress 9.325086e-05 
    ## ... Procrustes: rmse 0.0001887165  max resid 0.0004175529 
    ## ... Similar to previous best
    ## Run 219 stress 9.14084e-05 
    ## ... Procrustes: rmse 0.0002200482  max resid 0.0003744047 
    ## ... Similar to previous best
    ## Run 220 stress 9.535407e-05 
    ## ... Procrustes: rmse 7.617345e-05  max resid 0.0001440107 
    ## ... Similar to previous best
    ## Run 221 stress 9.337542e-05 
    ## ... Procrustes: rmse 0.0001486606  max resid 0.0002504577 
    ## ... Similar to previous best
    ## Run 222 stress 9.926504e-05 
    ## ... Procrustes: rmse 0.0002508368  max resid 0.0004074855 
    ## ... Similar to previous best
    ## Run 223 stress 8.977173e-05 
    ## ... Procrustes: rmse 0.0001396523  max resid 0.0002104491 
    ## ... Similar to previous best
    ## Run 224 stress 9.646551e-05 
    ## ... Procrustes: rmse 0.0002440965  max resid 0.0003898663 
    ## ... Similar to previous best
    ## Run 225 stress 9.236299e-05 
    ## ... Procrustes: rmse 0.0001545118  max resid 0.0002279178 
    ## ... Similar to previous best
    ## Run 226 stress 9.369059e-05 
    ## ... Procrustes: rmse 0.0001058767  max resid 0.0001604371 
    ## ... Similar to previous best
    ## Run 227 stress 9.055353e-05 
    ## ... Procrustes: rmse 0.0001738939  max resid 0.0002902296 
    ## ... Similar to previous best
    ## Run 228 stress 7.240638e-05 
    ## ... Procrustes: rmse 0.0001094984  max resid 0.0001688187 
    ## ... Similar to previous best
    ## Run 229 stress 9.747671e-05 
    ## ... Procrustes: rmse 0.0002118257  max resid 0.0003923629 
    ## ... Similar to previous best
    ## Run 230 stress 0.2854316 
    ## Run 231 stress 9.488098e-05 
    ## ... Procrustes: rmse 0.0002404066  max resid 0.0003840478 
    ## ... Similar to previous best
    ## Run 232 stress 9.786645e-05 
    ## ... Procrustes: rmse 0.000245128  max resid 0.0004016824 
    ## ... Similar to previous best
    ## Run 233 stress 8.872848e-05 
    ## ... Procrustes: rmse 9.674089e-05  max resid 0.0002066192 
    ## ... Similar to previous best
    ## Run 234 stress 9.46154e-05 
    ## ... Procrustes: rmse 0.0002007047  max resid 0.0003701221 
    ## ... Similar to previous best
    ## Run 235 stress 9.837893e-05 
    ## ... Procrustes: rmse 0.0001831197  max resid 0.0003139364 
    ## ... Similar to previous best
    ## Run 236 stress 9.528847e-05 
    ## ... Procrustes: rmse 0.0002064962  max resid 0.0003912332 
    ## ... Similar to previous best
    ## Run 237 stress 8.872324e-05 
    ## ... Procrustes: rmse 0.0001003614  max resid 0.0001487976 
    ## ... Similar to previous best
    ## Run 238 stress 9.634218e-05 
    ## ... Procrustes: rmse 0.0002432984  max resid 0.0003962597 
    ## ... Similar to previous best
    ## Run 239 stress 8.491754e-05 
    ## ... Procrustes: rmse 7.859249e-05  max resid 0.000165947 
    ## ... Similar to previous best
    ## Run 240 stress 8.7317e-05 
    ## ... Procrustes: rmse 6.569544e-05  max resid 0.0001216359 
    ## ... Similar to previous best
    ## Run 241 stress 9.894669e-05 
    ## ... Procrustes: rmse 0.000248076  max resid 0.0003931213 
    ## ... Similar to previous best
    ## Run 242 stress 9.344726e-05 
    ## ... Procrustes: rmse 0.0001754825  max resid 0.000305902 
    ## ... Similar to previous best
    ## Run 243 stress 9.981583e-05 
    ## ... Procrustes: rmse 0.0002484523  max resid 0.0004083389 
    ## ... Similar to previous best
    ## Run 244 stress 9.639873e-05 
    ## ... Procrustes: rmse 0.0002344025  max resid 0.0003868551 
    ## ... Similar to previous best
    ## Run 245 stress 9.905087e-05 
    ## ... Procrustes: rmse 0.0002468917  max resid 0.0003897492 
    ## ... Similar to previous best
    ## Run 246 stress 9.124359e-05 
    ## ... Procrustes: rmse 0.0001710704  max resid 0.0003012435 
    ## ... Similar to previous best
    ## Run 247 stress 9.682525e-05 
    ## ... Procrustes: rmse 0.0002353081  max resid 0.000371102 
    ## ... Similar to previous best
    ## Run 248 stress 9.383114e-05 
    ## ... Procrustes: rmse 0.0001754104  max resid 0.0003058322 
    ## ... Similar to previous best
    ## Run 249 stress 9.740138e-05 
    ## ... Procrustes: rmse 0.0002457604  max resid 0.0003990472 
    ## ... Similar to previous best
    ## Run 250 stress 9.442755e-05 
    ## ... Procrustes: rmse 0.0001496028  max resid 0.0002260681 
    ## ... Similar to previous best
    ## Run 251 stress 9.257969e-05 
    ## ... Procrustes: rmse 0.0002359842  max resid 0.0003787023 
    ## ... Similar to previous best
    ## Run 252 stress 9.876832e-05 
    ## ... Procrustes: rmse 0.0001570859  max resid 0.0003193605 
    ## ... Similar to previous best
    ## Run 253 stress 9.713602e-05 
    ## ... Procrustes: rmse 0.0002422115  max resid 0.0003953222 
    ## ... Similar to previous best
    ## Run 254 stress 8.756134e-05 
    ## ... Procrustes: rmse 9.289858e-05  max resid 0.000198544 
    ## ... Similar to previous best
    ## Run 255 stress 8.283849e-05 
    ## ... Procrustes: rmse 0.0001553889  max resid 0.0002843673 
    ## ... Similar to previous best
    ## Run 256 stress 8.519662e-05 
    ## ... Procrustes: rmse 0.0001402483  max resid 0.0002394784 
    ## ... Similar to previous best
    ## Run 257 stress 9.310269e-05 
    ## ... Procrustes: rmse 9.948851e-05  max resid 0.0001563175 
    ## ... Similar to previous best
    ## Run 258 stress 9.362166e-05 
    ## ... Procrustes: rmse 0.0002043477  max resid 0.0003829617 
    ## ... Similar to previous best
    ## Run 259 stress 8.973217e-05 
    ## ... Procrustes: rmse 0.0001004215  max resid 0.000214619 
    ## ... Similar to previous best
    ## Run 260 stress 9.703798e-05 
    ## ... Procrustes: rmse 0.0001998952  max resid 0.0004287535 
    ## ... Similar to previous best
    ## Run 261 stress 9.994626e-05 
    ## ... Procrustes: rmse 0.0002190391  max resid 0.0004040189 
    ## ... Similar to previous best
    ## Run 262 stress 0.2848512 
    ## Run 263 stress 0.2479794 
    ## Run 264 stress 9.663922e-05 
    ## ... Procrustes: rmse 0.000238251  max resid 0.0003745572 
    ## ... Similar to previous best
    ## Run 265 stress 9.586416e-05 
    ## ... Procrustes: rmse 0.0002427384  max resid 0.0003948091 
    ## ... Similar to previous best
    ## Run 266 stress 9.547246e-05 
    ## ... Procrustes: rmse 0.0001514167  max resid 0.0002921788 
    ## ... Similar to previous best
    ## Run 267 stress 9.016387e-05 
    ## ... Procrustes: rmse 0.0001960724  max resid 0.0003706453 
    ## ... Similar to previous best
    ## Run 268 stress 8.918501e-05 
    ## ... Procrustes: rmse 0.0001461284  max resid 0.0002846271 
    ## ... Similar to previous best
    ## Run 269 stress 0.2762713 
    ## Run 270 stress 9.91967e-05 
    ## ... Procrustes: rmse 0.0002009651  max resid 0.0003187205 
    ## ... Similar to previous best
    ## Run 271 stress 9.599846e-05 
    ## ... Procrustes: rmse 0.0002091626  max resid 0.000395279 
    ## ... Similar to previous best
    ## Run 272 stress 9.749907e-05 
    ## ... Procrustes: rmse 0.0002471446  max resid 0.0004012428 
    ## ... Similar to previous best
    ## Run 273 stress 9.845899e-05 
    ## ... Procrustes: rmse 0.0002423564  max resid 0.000381476 
    ## ... Similar to previous best
    ## Run 274 stress 9.385362e-05 
    ## ... Procrustes: rmse 0.000239094  max resid 0.0003843533 
    ## ... Similar to previous best
    ## Run 275 stress 9.687786e-05 
    ## ... Procrustes: rmse 0.0002458864  max resid 0.0004010613 
    ## ... Similar to previous best
    ## Run 276 stress 8.667917e-05 
    ## ... Procrustes: rmse 6.684867e-05  max resid 0.0001253361 
    ## ... Similar to previous best
    ## Run 277 stress 9.44607e-05 
    ## ... Procrustes: rmse 0.0001980898  max resid 0.0003644566 
    ## ... Similar to previous best
    ## Run 278 stress 9.41611e-05 
    ## ... Procrustes: rmse 0.0002385145  max resid 0.0003805459 
    ## ... Similar to previous best
    ## Run 279 stress 9.909426e-05 
    ## ... Procrustes: rmse 0.0002157873  max resid 0.0004026873 
    ## ... Similar to previous best
    ## Run 280 stress 9.289569e-05 
    ## ... Procrustes: rmse 0.0002032279  max resid 0.000384742 
    ## ... Similar to previous best
    ## Run 281 stress 9.566883e-05 
    ## ... Procrustes: rmse 0.0002341959  max resid 0.0003791912 
    ## ... Similar to previous best
    ## Run 282 stress 9.829949e-05 
    ## ... Procrustes: rmse 0.000182983  max resid 0.0003137963 
    ## ... Similar to previous best
    ## Run 283 stress 9.330763e-05 
    ## ... Procrustes: rmse 0.0002347743  max resid 0.0003702904 
    ## ... Similar to previous best
    ## Run 284 stress 9.825772e-05 
    ## ... Procrustes: rmse 0.0002159765  max resid 0.0004039396 
    ## ... Similar to previous best
    ## Run 285 stress 0.3082791 
    ## Run 286 stress 9.781628e-05 
    ## ... Procrustes: rmse 0.0002134042  max resid 0.0004016535 
    ## ... Similar to previous best
    ## Run 287 stress 8.610516e-05 
    ## ... Procrustes: rmse 9.752309e-05  max resid 0.000145583 
    ## ... Similar to previous best
    ## Run 288 stress 9.98388e-05 
    ## ... Procrustes: rmse 0.0002440618  max resid 0.0004094002 
    ## ... Similar to previous best
    ## Run 289 stress 8.974194e-05 
    ## ... Procrustes: rmse 0.0001697257  max resid 0.0003262924 
    ## ... Similar to previous best
    ## Run 290 stress 9.796117e-05 
    ## ... Procrustes: rmse 0.0002097985  max resid 0.0003853424 
    ## ... Similar to previous best
    ## Run 291 stress 9.91471e-05 
    ## ... Procrustes: rmse 0.0002088021  max resid 0.0003846232 
    ## ... Similar to previous best
    ## Run 292 stress 9.152698e-05 
    ## ... Procrustes: rmse 0.000197488  max resid 0.0003239425 
    ## ... Similar to previous best
    ## Run 293 stress 9.541948e-05 
    ## ... Procrustes: rmse 0.0002040508  max resid 0.0003774023 
    ## ... Similar to previous best
    ## Run 294 stress 9.110953e-05 
    ## ... Procrustes: rmse 0.0002321614  max resid 0.0003771435 
    ## ... Similar to previous best
    ## Run 295 stress 9.649208e-05 
    ## ... Procrustes: rmse 0.0002342319  max resid 0.0003869437 
    ## ... Similar to previous best
    ## Run 296 stress 9.472622e-05 
    ## ... Procrustes: rmse 0.0002396793  max resid 0.0003858475 
    ## ... Similar to previous best
    ## Run 297 stress 9.631726e-05 
    ## ... Procrustes: rmse 0.00024026  max resid 0.0003862825 
    ## ... Similar to previous best
    ## Run 298 stress 9.619847e-05 
    ## ... Procrustes: rmse 0.0001738888  max resid 0.0003120215 
    ## ... Similar to previous best
    ## Run 299 stress 9.464948e-05 
    ## ... Procrustes: rmse 0.0001741369  max resid 0.0003059168 
    ## ... Similar to previous best
    ## Run 300 stress 8.568588e-05 
    ## ... Procrustes: rmse 8.379466e-05  max resid 0.000170804 
    ## ... Similar to previous best
    ## Run 301 stress 9.789986e-05 
    ## ... Procrustes: rmse 0.0002118679  max resid 0.000400298 
    ## ... Similar to previous best
    ## Run 302 stress 6.323222e-05 
    ## ... Procrustes: rmse 4.60644e-05  max resid 6.910659e-05 
    ## ... Similar to previous best
    ## Run 303 stress 9.794056e-05 
    ## ... Procrustes: rmse 0.0001080747  max resid 0.0001630622 
    ## ... Similar to previous best
    ## Run 304 stress 9.802816e-05 
    ## ... Procrustes: rmse 0.0002454491  max resid 0.000389174 
    ## ... Similar to previous best
    ## Run 305 stress 9.818593e-05 
    ## ... Procrustes: rmse 0.0002485308  max resid 0.0003982521 
    ## ... Similar to previous best
    ## Run 306 stress 9.995853e-05 
    ## ... Procrustes: rmse 0.0002513219  max resid 0.0004089586 
    ## ... Similar to previous best
    ## Run 307 stress 9.516022e-05 
    ## ... Procrustes: rmse 0.0001071728  max resid 0.0001558809 
    ## ... Similar to previous best
    ## Run 308 stress 9.275484e-05 
    ## ... Procrustes: rmse 0.0001995964  max resid 0.0004005645 
    ## ... Similar to previous best
    ## Run 309 stress 9.626864e-05 
    ## ... Procrustes: rmse 0.0001175015  max resid 0.0001695979 
    ## ... Similar to previous best
    ## Run 310 stress 9.594052e-05 
    ## ... Procrustes: rmse 0.0002094692  max resid 0.0003901797 
    ## ... Similar to previous best
    ## Run 311 stress 9.890966e-05 
    ## ... Procrustes: rmse 0.0002079608  max resid 0.0003810934 
    ## ... Similar to previous best
    ## Run 312 stress 8.720547e-05 
    ## ... Procrustes: rmse 5.766917e-05  max resid 9.088687e-05 
    ## ... Similar to previous best
    ## Run 313 stress 8.481773e-05 
    ## ... Procrustes: rmse 6.949271e-05  max resid 0.0001454579 
    ## ... Similar to previous best
    ## Run 314 stress 9.224122e-05 
    ## ... Procrustes: rmse 0.000161101  max resid 0.000316077 
    ## ... Similar to previous best
    ## Run 315 stress 9.527591e-05 
    ## ... Procrustes: rmse 0.0002367894  max resid 0.0003745622 
    ## ... Similar to previous best
    ## Run 316 stress 9.250971e-05 
    ## ... Procrustes: rmse 0.0002003067  max resid 0.0003254064 
    ## ... Similar to previous best
    ## Run 317 stress 9.212352e-05 
    ## ... Procrustes: rmse 0.0002321372  max resid 0.0003723427 
    ## ... Similar to previous best
    ## Run 318 stress 8.662001e-05 
    ## ... Procrustes: rmse 9.398827e-05  max resid 0.0001998812 
    ## ... Similar to previous best
    ## Run 319 stress 9.047136e-05 
    ## ... Procrustes: rmse 0.0002304333  max resid 0.0003761472 
    ## ... Similar to previous best
    ## Run 320 stress 6.399233e-05 
    ## ... Procrustes: rmse 7.782733e-05  max resid 0.000146325 
    ## ... Similar to previous best
    ## Run 321 stress 9.497105e-05 
    ## ... Procrustes: rmse 0.0002411074  max resid 0.0003878508 
    ## ... Similar to previous best
    ## Run 322 stress 0.2848516 
    ## Run 323 stress 8.988978e-05 
    ## ... Procrustes: rmse 7.33572e-05  max resid 0.0001182059 
    ## ... Similar to previous best
    ## Run 324 stress 8.651666e-05 
    ## ... Procrustes: rmse 0.000186909  max resid 0.0003573358 
    ## ... Similar to previous best
    ## Run 325 stress 8.280876e-05 
    ## ... Procrustes: rmse 9.963843e-05  max resid 0.0001690807 
    ## ... Similar to previous best
    ## Run 326 stress 9.692954e-05 
    ## ... Procrustes: rmse 0.0002077076  max resid 0.000394999 
    ## ... Similar to previous best
    ## Run 327 stress 9.480663e-05 
    ## ... Procrustes: rmse 0.0002372326  max resid 0.0003886854 
    ## ... Similar to previous best
    ## Run 328 stress 7.90169e-05 
    ## ... Procrustes: rmse 8.777156e-05  max resid 0.0001891705 
    ## ... Similar to previous best
    ## Run 329 stress 9.629481e-05 
    ## ... Procrustes: rmse 0.0001804372  max resid 0.0003169083 
    ## ... Similar to previous best
    ## Run 330 stress 9.301009e-05 
    ## ... Procrustes: rmse 0.0001733819  max resid 0.000303023 
    ## ... Similar to previous best
    ## Run 331 stress 9.406121e-05 
    ## ... Procrustes: rmse 0.0002365735  max resid 0.0003918807 
    ## ... Similar to previous best
    ## Run 332 stress 9.499743e-05 
    ## ... Procrustes: rmse 0.0002037297  max resid 0.0003247977 
    ## ... Similar to previous best
    ## Run 333 stress 8.976978e-05 
    ## ... Procrustes: rmse 7.713775e-05  max resid 0.0001498359 
    ## ... Similar to previous best
    ## Run 334 stress 9.530578e-05 
    ## ... Procrustes: rmse 0.0002077649  max resid 0.0003918982 
    ## ... Similar to previous best
    ## Run 335 stress 0.2846007 
    ## Run 336 stress 9.930039e-05 
    ## ... Procrustes: rmse 0.0001844606  max resid 0.000315829 
    ## ... Similar to previous best
    ## Run 337 stress 9.376408e-05 
    ## ... Procrustes: rmse 0.0002340737  max resid 0.0003761151 
    ## ... Similar to previous best
    ## Run 338 stress 9.396604e-05 
    ## ... Procrustes: rmse 0.0001728483  max resid 0.0003033235 
    ## ... Similar to previous best
    ## Run 339 stress 9.541934e-05 
    ## ... Procrustes: rmse 8.122758e-05  max resid 0.0001574839 
    ## ... Similar to previous best
    ## Run 340 stress 8.858977e-05 
    ## ... Procrustes: rmse 0.0002265637  max resid 0.0003673818 
    ## ... Similar to previous best
    ## Run 341 stress 9.35717e-05 
    ## ... Procrustes: rmse 0.0002365006  max resid 0.00038312 
    ## ... Similar to previous best
    ## Run 342 stress 9.867327e-05 
    ## ... Procrustes: rmse 0.0002165641  max resid 0.0004057935 
    ## ... Similar to previous best
    ## Run 343 stress 8.979622e-05 
    ## ... Procrustes: rmse 9.04512e-05  max resid 0.0001555873 
    ## ... Similar to previous best
    ## Run 344 stress 8.697959e-05 
    ## ... Procrustes: rmse 9.196616e-05  max resid 0.0001984965 
    ## ... Similar to previous best
    ## Run 345 stress 9.524124e-05 
    ## ... Procrustes: rmse 0.0001052648  max resid 0.0001514172 
    ## ... Similar to previous best
    ## Run 346 stress 9.719505e-05 
    ## ... Procrustes: rmse 0.0002119325  max resid 0.0003934925 
    ## ... Similar to previous best
    ## Run 347 stress 9.358168e-05 
    ## ... Procrustes: rmse 0.0001656011  max resid 0.0002904917 
    ## ... Similar to previous best
    ## Run 348 stress 9.971941e-05 
    ## ... Procrustes: rmse 0.0001518733  max resid 0.0002540473 
    ## ... Similar to previous best
    ## Run 349 stress 9.597338e-05 
    ## ... Procrustes: rmse 0.0001533078  max resid 0.0002598896 
    ## ... Similar to previous best
    ## Run 350 stress 9.172531e-05 
    ## ... Procrustes: rmse 0.0001998406  max resid 0.0003744788 
    ## ... Similar to previous best
    ## Run 351 stress 9.96187e-05 
    ## ... Procrustes: rmse 0.0001791795  max resid 0.0003206484 
    ## ... Similar to previous best
    ## Run 352 stress 9.443845e-05 
    ## ... Procrustes: rmse 0.0002355401  max resid 0.0003720354 
    ## ... Similar to previous best
    ## Run 353 stress 9.750775e-05 
    ## ... Procrustes: rmse 0.0001822264  max resid 0.0003128652 
    ## ... Similar to previous best
    ## Run 354 stress 8.904551e-05 
    ## ... Procrustes: rmse 9.099192e-05  max resid 0.0001888355 
    ## ... Similar to previous best
    ## Run 355 stress 9.732031e-05 
    ## ... Procrustes: rmse 0.0001783457  max resid 0.0003146273 
    ## ... Similar to previous best
    ## Run 356 stress 9.153624e-05 
    ## ... Procrustes: rmse 0.0001935021  max resid 0.0003899706 
    ## ... Similar to previous best
    ## Run 357 stress 9.836846e-05 
    ## ... Procrustes: rmse 0.0001827224  max resid 0.0003435121 
    ## ... Similar to previous best
    ## Run 358 stress 9.683023e-05 
    ## ... Procrustes: rmse 0.0002450777  max resid 0.0003984915 
    ## ... Similar to previous best
    ## Run 359 stress 9.954662e-05 
    ## ... Procrustes: rmse 0.0001797326  max resid 0.0003516538 
    ## ... Similar to previous best
    ## Run 360 stress 9.995066e-05 
    ## ... Procrustes: rmse 0.0002134505  max resid 0.0003530507 
    ## ... Similar to previous best
    ## Run 361 stress 9.622638e-05 
    ## ... Procrustes: rmse 0.0001794581  max resid 0.000310104 
    ## ... Similar to previous best
    ## Run 362 stress 9.361937e-05 
    ## ... Procrustes: rmse 0.0002361547  max resid 0.0003859635 
    ## ... Similar to previous best
    ## Run 363 stress 9.823937e-05 
    ## ... Procrustes: rmse 0.0002111846  max resid 0.0003975131 
    ## ... Similar to previous best
    ## Run 364 stress 8.242523e-05 
    ## ... Procrustes: rmse 9.139549e-05  max resid 0.0001951607 
    ## ... Similar to previous best
    ## Run 365 stress 0.3083098 
    ## Run 366 stress 9.364998e-05 
    ## ... Procrustes: rmse 0.0001730674  max resid 0.0003036338 
    ## ... Similar to previous best
    ## Run 367 stress 9.451504e-05 
    ## ... Procrustes: rmse 0.0002356257  max resid 0.000385821 
    ## ... Similar to previous best
    ## Run 368 stress 9.819842e-05 
    ## ... Procrustes: rmse 8.744298e-05  max resid 0.0001706351 
    ## ... Similar to previous best
    ## Run 369 stress 8.897634e-05 
    ## ... Procrustes: rmse 0.0002256013  max resid 0.0003598824 
    ## ... Similar to previous best
    ## Run 370 stress 0.2319225 
    ## Run 371 stress 9.518043e-05 
    ## ... Procrustes: rmse 0.0002370131  max resid 0.0003780228 
    ## ... Similar to previous best
    ## Run 372 stress 9.774998e-05 
    ## ... Procrustes: rmse 0.0002483286  max resid 0.0004013734 
    ## ... Similar to previous best
    ## Run 373 stress 9.620142e-05 
    ## ... Procrustes: rmse 0.0002018096  max resid 0.0004214057 
    ## ... Similar to previous best
    ## Run 374 stress 9.108885e-05 
    ## ... Procrustes: rmse 0.0001892038  max resid 0.0003722687 
    ## ... Similar to previous best
    ## Run 375 stress 8.828345e-05 
    ## ... Procrustes: rmse 9.298985e-05  max resid 0.0001937412 
    ## ... Similar to previous best
    ## Run 376 stress 9.204831e-05 
    ## ... Procrustes: rmse 0.0002342688  max resid 0.0003771238 
    ## ... Similar to previous best
    ## Run 377 stress 8.986305e-05 
    ## ... Procrustes: rmse 0.0001798684  max resid 0.0003384858 
    ## ... Similar to previous best
    ## Run 378 stress 9.843324e-05 
    ## ... Procrustes: rmse 0.0002373624  max resid 0.0003936087 
    ## ... Similar to previous best
    ## Run 379 stress 9.587377e-05 
    ## ... Procrustes: rmse 0.0002103833  max resid 0.000394388 
    ## ... Similar to previous best
    ## Run 380 stress 8.683325e-05 
    ## ... Procrustes: rmse 9.478532e-05  max resid 0.0001337732 
    ## ... Similar to previous best
    ## Run 381 stress 8.971776e-05 
    ## ... Procrustes: rmse 6.362318e-05  max resid 0.0001061828 
    ## ... Similar to previous best
    ## Run 382 stress 9.644785e-05 
    ## ... Procrustes: rmse 0.0002449663  max resid 0.0003980627 
    ## ... Similar to previous best
    ## Run 383 stress 9.637322e-05 
    ## ... Procrustes: rmse 0.0002102256  max resid 0.0003956758 
    ## ... Similar to previous best
    ## Run 384 stress 9.851436e-05 
    ## ... Procrustes: rmse 0.0002499349  max resid 0.000404453 
    ## ... Similar to previous best
    ## Run 385 stress 9.897242e-05 
    ## ... Procrustes: rmse 0.0001782029  max resid 0.0003132352 
    ## ... Similar to previous best
    ## Run 386 stress 9.831174e-05 
    ## ... Procrustes: rmse 8.227649e-05  max resid 0.0001459135 
    ## ... Similar to previous best
    ## Run 387 stress 9.775621e-05 
    ## ... Procrustes: rmse 0.0001966651  max resid 0.0003641331 
    ## ... Similar to previous best
    ## Run 388 stress 7.723371e-05 
    ## ... Procrustes: rmse 8.932496e-05  max resid 0.0001534024 
    ## ... Similar to previous best
    ## Run 389 stress 9.378251e-05 
    ## ... Procrustes: rmse 0.0002372447  max resid 0.0003860868 
    ## ... Similar to previous best
    ## Run 390 stress 0.2319166 
    ## Run 391 stress 0.2842805 
    ## Run 392 stress 9.713942e-05 
    ## ... Procrustes: rmse 0.0002464737  max resid 0.0003999429 
    ## ... Similar to previous best
    ## Run 393 stress 9.992389e-05 
    ## ... Procrustes: rmse 0.0002434817  max resid 0.0003812598 
    ## ... Similar to previous best
    ## Run 394 stress 9.752478e-05 
    ## ... Procrustes: rmse 0.000227284  max resid 0.0003679019 
    ## ... Similar to previous best
    ## Run 395 stress 9.533061e-05 
    ## ... Procrustes: rmse 0.0001749884  max resid 0.0003076755 
    ## ... Similar to previous best
    ## Run 396 stress 9.995806e-05 
    ## ... Procrustes: rmse 0.0002189004  max resid 0.0004131459 
    ## ... Similar to previous best
    ## Run 397 stress 9.53376e-05 
    ## ... Procrustes: rmse 0.0002300993  max resid 0.0003643889 
    ## ... Similar to previous best
    ## Run 398 stress 9.0833e-05 
    ## ... Procrustes: rmse 0.0001894995  max resid 0.0003679723 
    ## ... Similar to previous best
    ## Run 399 stress 9.999451e-05 
    ## ... Procrustes: rmse 0.0002167585  max resid 0.0004094517 
    ## ... Similar to previous best
    ## Run 400 stress 9.022288e-05 
    ## ... Procrustes: rmse 0.0001672033  max resid 0.0002951522 
    ## ... Similar to previous best
    ## Run 401 stress 9.221386e-05 
    ## ... Procrustes: rmse 0.0001038289  max resid 0.0001505538 
    ## ... Similar to previous best
    ## Run 402 stress 9.750704e-05 
    ## ... Procrustes: rmse 0.0002052215  max resid 0.0004296237 
    ## ... Similar to previous best
    ## Run 403 stress 9.517657e-05 
    ## ... Procrustes: rmse 0.0001748197  max resid 0.00030738 
    ## ... Similar to previous best
    ## Run 404 stress 9.747265e-05 
    ## ... Procrustes: rmse 0.000212053  max resid 0.0003937648 
    ## ... Similar to previous best
    ## Run 405 stress 9.113848e-05 
    ## ... Procrustes: rmse 0.0001895474  max resid 0.0003706721 
    ## ... Similar to previous best
    ## Run 406 stress 9.651863e-05 
    ## ... Procrustes: rmse 0.0002020168  max resid 0.0003708891 
    ## ... Similar to previous best
    ## Run 407 stress 9.823068e-05 
    ## ... Procrustes: rmse 0.0002149861  max resid 0.0004038527 
    ## ... Similar to previous best
    ## Run 408 stress 9.599793e-05 
    ## ... Procrustes: rmse 0.0002099  max resid 0.0003991688 
    ## ... Similar to previous best
    ## Run 409 stress 9.896143e-05 
    ## ... Procrustes: rmse 0.000247764  max resid 0.0003875231 
    ## ... Similar to previous best
    ## Run 410 stress 9.675517e-05 
    ## ... Procrustes: rmse 0.0002070661  max resid 0.0003811912 
    ## ... Similar to previous best
    ## Run 411 stress 9.578445e-05 
    ## ... Procrustes: rmse 0.0002436996  max resid 0.0003930189 
    ## ... Similar to previous best
    ## Run 412 stress 9.533229e-05 
    ## ... Procrustes: rmse 0.000152955  max resid 0.0002595777 
    ## ... Similar to previous best
    ## Run 413 stress 9.786945e-05 
    ## ... Procrustes: rmse 0.0002344439  max resid 0.0003869226 
    ## ... Similar to previous best
    ## Run 414 stress 9.805437e-05 
    ## ... Procrustes: rmse 0.0001826246  max resid 0.0003136647 
    ## ... Similar to previous best
    ## Run 415 stress 8.332526e-05 
    ## ... Procrustes: rmse 0.0001765326  max resid 0.0003396498 
    ## ... Similar to previous best
    ## Run 416 stress 9.224865e-05 
    ## ... Procrustes: rmse 0.0002324696  max resid 0.0003755664 
    ## ... Similar to previous best
    ## Run 417 stress 9.536886e-05 
    ## ... Procrustes: rmse 0.0002059014  max resid 0.0003838103 
    ## ... Similar to previous best
    ## Run 418 stress 8.938007e-05 
    ## ... Procrustes: rmse 6.857102e-05  max resid 0.0001047992 
    ## ... Similar to previous best
    ## Run 419 stress 9.947359e-05 
    ## ... Procrustes: rmse 0.0001742192  max resid 0.0003453742 
    ## ... Similar to previous best
    ## Run 420 stress 9.223817e-05 
    ## ... Procrustes: rmse 0.0002338064  max resid 0.0003709827 
    ## ... Similar to previous best
    ## Run 421 stress 9.419583e-05 
    ## ... Procrustes: rmse 0.0002047022  max resid 0.0003805067 
    ## ... Similar to previous best
    ## Run 422 stress 9.490838e-05 
    ## ... Procrustes: rmse 0.0002369843  max resid 0.0003894211 
    ## ... Similar to previous best
    ## Run 423 stress 9.811497e-05 
    ## ... Procrustes: rmse 0.0001749483  max resid 0.000314998 
    ## ... Similar to previous best
    ## Run 424 stress 9.603006e-05 
    ## ... Procrustes: rmse 0.0002415901  max resid 0.0003842528 
    ## ... Similar to previous best
    ## Run 425 stress 9.387571e-05 
    ## ... Procrustes: rmse 0.000233962  max resid 0.0003841897 
    ## ... Similar to previous best
    ## Run 426 stress 9.645042e-05 
    ## ... Procrustes: rmse 0.0002427317  max resid 0.0003928095 
    ## ... Similar to previous best
    ## Run 427 stress 9.33546e-05 
    ## ... Procrustes: rmse 0.0001971796  max resid 0.0003959887 
    ## ... Similar to previous best
    ## Run 428 stress 8.921284e-05 
    ## ... Procrustes: rmse 0.0002278952  max resid 0.0003700501 
    ## ... Similar to previous best
    ## Run 429 stress 9.773374e-05 
    ## ... Procrustes: rmse 0.0002036389  max resid 0.0003779983 
    ## ... Similar to previous best
    ## Run 430 stress 8.92935e-05 
    ## ... Procrustes: rmse 7.242382e-05  max resid 0.0001386323 
    ## ... Similar to previous best
    ## Run 431 stress 8.839908e-05 
    ## ... Procrustes: rmse 0.0001247785  max resid 0.0001933617 
    ## ... Similar to previous best
    ## Run 432 stress 9.669593e-05 
    ## ... Procrustes: rmse 0.0001701133  max resid 0.0003319853 
    ## ... Similar to previous best
    ## Run 433 stress 9.479934e-05 
    ## ... Procrustes: rmse 0.0001742286  max resid 0.0003061596 
    ## ... Similar to previous best
    ## Run 434 stress 9.602532e-05 
    ## ... Procrustes: rmse 0.0002007949  max resid 0.0003690336 
    ## ... Similar to previous best
    ## Run 435 stress 9.440626e-05 
    ## ... Procrustes: rmse 0.0001639641  max resid 0.0002975474 
    ## ... Similar to previous best
    ## Run 436 stress 9.651738e-05 
    ## ... Procrustes: rmse 0.0002110535  max resid 0.000397018 
    ## ... Similar to previous best
    ## Run 437 stress 9.608946e-05 
    ## ... Procrustes: rmse 0.0002430506  max resid 0.000395273 
    ## ... Similar to previous best
    ## Run 438 stress 9.5646e-05 
    ## ... Procrustes: rmse 0.0002412881  max resid 0.0003872275 
    ## ... Similar to previous best
    ## Run 439 stress 9.94451e-05 
    ## ... Procrustes: rmse 0.0001762413  max resid 0.0003490832 
    ## ... Similar to previous best
    ## Run 440 stress 9.381479e-05 
    ## ... Procrustes: rmse 0.0002034115  max resid 0.000313255 
    ## ... Similar to previous best
    ## Run 441 stress 0.2319172 
    ## Run 442 stress 9.885708e-05 
    ## ... Procrustes: rmse 0.0002479098  max resid 0.0004054058 
    ## ... Similar to previous best
    ## Run 443 stress 9.336824e-05 
    ## ... Procrustes: rmse 0.0001756423  max resid 0.0003061266 
    ## ... Similar to previous best
    ## Run 444 stress 9.756932e-05 
    ## ... Procrustes: rmse 0.0002091411  max resid 0.0003348615 
    ## ... Similar to previous best
    ## Run 445 stress 9.717445e-05 
    ## ... Procrustes: rmse 0.0001823463  max resid 0.0003138342 
    ## ... Similar to previous best
    ## Run 446 stress 9.490929e-05 
    ## ... Procrustes: rmse 0.0001767321  max resid 0.0003077298 
    ## ... Similar to previous best
    ## Run 447 stress 9.397786e-05 
    ## ... Procrustes: rmse 0.0002043564  max resid 0.000387238 
    ## ... Similar to previous best
    ## Run 448 stress 9.845534e-05 
    ## ... Procrustes: rmse 0.0001074736  max resid 0.0001520252 
    ## ... Similar to previous best
    ## Run 449 stress 9.25946e-05 
    ## ... Procrustes: rmse 0.0001733323  max resid 0.000303645 
    ## ... Similar to previous best
    ## Run 450 stress 9.533256e-05 
    ## ... Procrustes: rmse 0.0001048116  max resid 0.0001495573 
    ## ... Similar to previous best
    ## Run 451 stress 9.414882e-05 
    ## ... Procrustes: rmse 0.0002354028  max resid 0.0003870891 
    ## ... Similar to previous best
    ## Run 452 stress 9.290713e-05 
    ## ... Procrustes: rmse 0.0001741049  max resid 0.0003044061 
    ## ... Similar to previous best
    ## Run 453 stress 8.78625e-05 
    ## ... Procrustes: rmse 0.0002172225  max resid 0.0003451122 
    ## ... Similar to previous best
    ## Run 454 stress 8.422302e-05 
    ## ... Procrustes: rmse 0.000146022  max resid 0.000294038 
    ## ... Similar to previous best
    ## Run 455 stress 9.683055e-05 
    ## ... Procrustes: rmse 0.0002383923  max resid 0.0003746558 
    ## ... Similar to previous best
    ## Run 456 stress 9.286814e-05 
    ## ... Procrustes: rmse 0.0001964003  max resid 0.000315652 
    ## ... Similar to previous best
    ## Run 457 stress 6.541472e-05 
    ## ... Procrustes: rmse 7.243756e-05  max resid 0.0001240794 
    ## ... Similar to previous best
    ## Run 458 stress 8.780023e-05 
    ## ... Procrustes: rmse 0.000223507  max resid 0.0003606625 
    ## ... Similar to previous best
    ## Run 459 stress 0.3083098 
    ## Run 460 stress 9.138105e-05 
    ## ... Procrustes: rmse 0.000165931  max resid 0.0003242371 
    ## ... Similar to previous best
    ## Run 461 stress 9.675813e-05 
    ## ... Procrustes: rmse 5.698518e-05  max resid 8.87146e-05 
    ## ... Similar to previous best
    ## Run 462 stress 9.896114e-05 
    ## ... Procrustes: rmse 0.0001537843  max resid 0.0002325912 
    ## ... Similar to previous best
    ## Run 463 stress 9.126005e-05 
    ## ... Procrustes: rmse 0.0002274032  max resid 0.0003709659 
    ## ... Similar to previous best
    ## Run 464 stress 9.588888e-05 
    ## ... Procrustes: rmse 0.0002431459  max resid 0.000395865 
    ## ... Similar to previous best
    ## Run 465 stress 9.434101e-05 
    ## ... Procrustes: rmse 0.0002353404  max resid 0.0003717133 
    ## ... Similar to previous best
    ## Run 466 stress 8.744046e-05 
    ## ... Procrustes: rmse 0.0001275781  max resid 0.0002033113 
    ## ... Similar to previous best
    ## Run 467 stress 9.503198e-05 
    ## ... Procrustes: rmse 0.0002404318  max resid 0.000381649 
    ## ... Similar to previous best
    ## Run 468 stress 9.036352e-05 
    ## ... Procrustes: rmse 0.0001591223  max resid 0.0003039605 
    ## ... Similar to previous best
    ## Run 469 stress 9.865452e-05 
    ## ... Procrustes: rmse 0.0001977071  max resid 0.000383338 
    ## ... Similar to previous best
    ## Run 470 stress 8.395004e-05 
    ## ... Procrustes: rmse 9.184568e-05  max resid 0.0001986811 
    ## ... Similar to previous best
    ## Run 471 stress 9.882219e-05 
    ## ... Procrustes: rmse 0.0001954322  max resid 0.0003810218 
    ## ... Similar to previous best
    ## Run 472 stress 0.3083098 
    ## Run 473 stress 9.731456e-05 
    ## ... Procrustes: rmse 0.000245428  max resid 0.0004005534 
    ## ... Similar to previous best
    ## Run 474 stress 9.160219e-05 
    ## ... Procrustes: rmse 0.0001604429  max resid 0.0002942453 
    ## ... Similar to previous best
    ## Run 475 stress 8.818047e-05 
    ## ... Procrustes: rmse 0.0001422918  max resid 0.0002261182 
    ## ... Similar to previous best
    ## Run 476 stress 9.539252e-05 
    ## ... Procrustes: rmse 0.0002067642  max resid 0.0003921304 
    ## ... Similar to previous best
    ## Run 477 stress 9.333163e-05 
    ## ... Procrustes: rmse 0.0001944318  max resid 0.0003089408 
    ## ... Similar to previous best
    ## Run 478 stress 9.053644e-05 
    ## ... Procrustes: rmse 0.0001718197  max resid 0.0003017284 
    ## ... Similar to previous best
    ## Run 479 stress 8.611418e-05 
    ## ... Procrustes: rmse 0.0001808801  max resid 0.0003544422 
    ## ... Similar to previous best
    ## Run 480 stress 8.524785e-05 
    ## ... Procrustes: rmse 6.07813e-05  max resid 0.0001214684 
    ## ... Similar to previous best
    ## Run 481 stress 9.689219e-05 
    ## ... Procrustes: rmse 0.0001806974  max resid 0.0003113777 
    ## ... Similar to previous best
    ## Run 482 stress 9.931631e-05 
    ## ... Procrustes: rmse 0.0002447789  max resid 0.0004060269 
    ## ... Similar to previous best
    ## Run 483 stress 9.354815e-05 
    ## ... Procrustes: rmse 0.0002365857  max resid 0.0003886423 
    ## ... Similar to previous best
    ## Run 484 stress 8.546875e-05 
    ## ... Procrustes: rmse 6.836052e-05  max resid 0.0001299655 
    ## ... Similar to previous best
    ## Run 485 stress 9.671087e-05 
    ## ... Procrustes: rmse 0.0002453884  max resid 0.0003981385 
    ## ... Similar to previous best
    ## Run 486 stress 9.310917e-05 
    ## ... Procrustes: rmse 0.0002281393  max resid 0.0003544902 
    ## ... Similar to previous best
    ## Run 487 stress 0.283906 
    ## Run 488 stress 9.993718e-05 
    ## ... Procrustes: rmse 0.0001281817  max resid 0.000182693 
    ## ... Similar to previous best
    ## Run 489 stress 8.4591e-05 
    ## ... Procrustes: rmse 9.116582e-05  max resid 0.0001913179 
    ## ... Similar to previous best
    ## Run 490 stress 9.62339e-05 
    ## ... Procrustes: rmse 0.0002433371  max resid 0.0003873812 
    ## ... Similar to previous best
    ## Run 491 stress 9.840575e-05 
    ## ... Procrustes: rmse 0.0001827936  max resid 0.0003222682 
    ## ... Similar to previous best
    ## Run 492 stress 9.465332e-05 
    ## ... Procrustes: rmse 0.0001942479  max resid 0.000382741 
    ## ... Similar to previous best
    ## Run 493 stress 9.941205e-05 
    ## ... Procrustes: rmse 9.065609e-05  max resid 0.0001806845 
    ## ... Similar to previous best
    ## Run 494 stress 9.918479e-05 
    ## ... Procrustes: rmse 0.0001530266  max resid 0.0002285505 
    ## ... Similar to previous best
    ## Run 495 stress 9.870255e-05 
    ## ... Procrustes: rmse 0.0002506874  max resid 0.0004055328 
    ## ... Similar to previous best
    ## Run 496 stress 9.276212e-05 
    ## ... Procrustes: rmse 0.0002024787  max resid 0.0003808835 
    ## ... Similar to previous best
    ## Run 497 stress 9.011369e-05 
    ## ... Procrustes: rmse 0.0002298648  max resid 0.0003726641 
    ## ... Similar to previous best
    ## Run 498 stress 9.67157e-05 
    ## ... Procrustes: rmse 0.0001789733  max resid 0.0003149138 
    ## ... Similar to previous best
    ## Run 499 stress 9.138397e-05 
    ## ... Procrustes: rmse 0.0002297997  max resid 0.0003796049 
    ## ... Similar to previous best
    ## Run 500 stress 9.475445e-05 
    ## ... Procrustes: rmse 7.466487e-05  max resid 0.0001469053 
    ## ... Similar to previous best
    ## *** Best solution repeated 351 times

    ## Warning in metaMDS(PD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
PD_beta_geo_NMDS <- metaMDS(PD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08080643 
    ## Run 1 stress 0.08121396 
    ## ... Procrustes: rmse 0.02027908  max resid 0.06580645 
    ## Run 2 stress 0.081318 
    ## Run 3 stress 0.09563364 
    ## Run 4 stress 0.08121402 
    ## ... Procrustes: rmse 0.02026566  max resid 0.06577583 
    ## Run 5 stress 0.0929675 
    ## Run 6 stress 0.09124926 
    ## Run 7 stress 0.08152294 
    ## Run 8 stress 0.09270878 
    ## Run 9 stress 0.09563388 
    ## Run 10 stress 0.08131795 
    ## Run 11 stress 0.08080641 
    ## ... New best solution
    ## ... Procrustes: rmse 3.129413e-05  max resid 0.0001020405 
    ## ... Similar to previous best
    ## Run 12 stress 0.08085076 
    ## ... Procrustes: rmse 0.009594926  max resid 0.03305794 
    ## Run 13 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027546  max resid 0.06584248 
    ## Run 14 stress 0.09180031 
    ## Run 15 stress 0.08085077 
    ## ... Procrustes: rmse 0.009611593  max resid 0.03313903 
    ## Run 16 stress 0.09296726 
    ## Run 17 stress 0.08085079 
    ## ... Procrustes: rmse 0.009574452  max resid 0.03293812 
    ## Run 18 stress 0.09563365 
    ## Run 19 stress 0.08152292 
    ## Run 20 stress 0.1029169 
    ## Run 21 stress 0.09559571 
    ## Run 22 stress 0.08152299 
    ## Run 23 stress 0.08080643 
    ## ... Procrustes: rmse 4.178836e-05  max resid 0.0001561572 
    ## ... Similar to previous best
    ## Run 24 stress 0.08152291 
    ## Run 25 stress 0.08080642 
    ## ... Procrustes: rmse 1.225352e-05  max resid 2.398662e-05 
    ## ... Similar to previous best
    ## Run 26 stress 0.08152296 
    ## Run 27 stress 0.0813179 
    ## Run 28 stress 0.09129242 
    ## Run 29 stress 0.08085074 
    ## ... Procrustes: rmse 0.009596391  max resid 0.033071 
    ## Run 30 stress 0.08080643 
    ## ... Procrustes: rmse 3.159221e-05  max resid 0.0001121091 
    ## ... Similar to previous best
    ## Run 31 stress 0.08121397 
    ## ... Procrustes: rmse 0.0202715  max resid 0.06577171 
    ## Run 32 stress 0.09457609 
    ## Run 33 stress 0.08080655 
    ## ... Procrustes: rmse 0.0002181002  max resid 0.0007848341 
    ## ... Similar to previous best
    ## Run 34 stress 0.09120964 
    ## Run 35 stress 0.08085077 
    ## ... Procrustes: rmse 0.009604889  max resid 0.03313459 
    ## Run 36 stress 0.1043987 
    ## Run 37 stress 0.09129243 
    ## Run 38 stress 0.1021228 
    ## Run 39 stress 0.1051522 
    ## Run 40 stress 0.08080643 
    ## ... Procrustes: rmse 0.0001109547  max resid 0.0004028407 
    ## ... Similar to previous best
    ## Run 41 stress 0.09129239 
    ## Run 42 stress 0.08121398 
    ## ... Procrustes: rmse 0.02024674  max resid 0.06572407 
    ## Run 43 stress 0.08080652 
    ## ... Procrustes: rmse 0.0001346887  max resid 0.0004848097 
    ## ... Similar to previous best
    ## Run 44 stress 0.08152291 
    ## Run 45 stress 0.08085074 
    ## ... Procrustes: rmse 0.009602617  max resid 0.0331029 
    ## Run 46 stress 0.09180031 
    ## Run 47 stress 0.08080648 
    ## ... Procrustes: rmse 9.428892e-05  max resid 0.000340725 
    ## ... Similar to previous best
    ## Run 48 stress 0.105976 
    ## Run 49 stress 0.08080646 
    ## ... Procrustes: rmse 7.727643e-05  max resid 0.0002782969 
    ## ... Similar to previous best
    ## Run 50 stress 0.09563366 
    ## Run 51 stress 0.09178837 
    ## Run 52 stress 0.1067885 
    ## Run 53 stress 0.09124912 
    ## Run 54 stress 0.09283153 
    ## Run 55 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027278  max resid 0.06579776 
    ## Run 56 stress 0.09236463 
    ## Run 57 stress 0.09180034 
    ## Run 58 stress 0.09145261 
    ## Run 59 stress 0.08121396 
    ## ... Procrustes: rmse 0.02024996  max resid 0.06563746 
    ## Run 60 stress 0.09563378 
    ## Run 61 stress 0.09120965 
    ## Run 62 stress 0.08121402 
    ## ... Procrustes: rmse 0.02031665  max resid 0.06608109 
    ## Run 63 stress 0.08085075 
    ## ... Procrustes: rmse 0.009595477  max resid 0.03307365 
    ## Run 64 stress 0.09236468 
    ## Run 65 stress 0.08121401 
    ## ... Procrustes: rmse 0.02029315  max resid 0.0660417 
    ## Run 66 stress 0.1090857 
    ## Run 67 stress 0.08085078 
    ## ... Procrustes: rmse 0.009609069  max resid 0.03312174 
    ## Run 68 stress 0.09236468 
    ## Run 69 stress 0.09261215 
    ## Run 70 stress 0.08131789 
    ## Run 71 stress 0.08152291 
    ## Run 72 stress 0.08131806 
    ## Run 73 stress 0.08121402 
    ## ... Procrustes: rmse 0.0202404  max resid 0.06551825 
    ## Run 74 stress 0.1061419 
    ## Run 75 stress 0.09129242 
    ## Run 76 stress 0.08121395 
    ## ... Procrustes: rmse 0.02026505  max resid 0.06575368 
    ## Run 77 stress 0.08085074 
    ## ... Procrustes: rmse 0.009592936  max resid 0.03304966 
    ## Run 78 stress 0.1096608 
    ## Run 79 stress 0.09387655 
    ## Run 80 stress 0.09145255 
    ## Run 81 stress 0.0928316 
    ## Run 82 stress 0.08121399 
    ## ... Procrustes: rmse 0.02030451  max resid 0.06594998 
    ## Run 83 stress 0.09236466 
    ## Run 84 stress 0.09120968 
    ## Run 85 stress 0.08085076 
    ## ... Procrustes: rmse 0.009585992  max resid 0.03300662 
    ## Run 86 stress 0.09178824 
    ## Run 87 stress 0.08080641 
    ## ... Procrustes: rmse 1.905617e-05  max resid 5.084513e-05 
    ## ... Similar to previous best
    ## Run 88 stress 0.08121397 
    ## ... Procrustes: rmse 0.02024717  max resid 0.06561033 
    ## Run 89 stress 0.08080642 
    ## ... Procrustes: rmse 1.116033e-05  max resid 3.338014e-05 
    ## ... Similar to previous best
    ## Run 90 stress 0.09314321 
    ## Run 91 stress 0.08152291 
    ## Run 92 stress 0.09124914 
    ## Run 93 stress 0.08131794 
    ## Run 94 stress 0.08080646 
    ## ... Procrustes: rmse 0.0001416685  max resid 0.0005109715 
    ## ... Similar to previous best
    ## Run 95 stress 0.08080643 
    ## ... Procrustes: rmse 0.0001014211  max resid 0.0003584823 
    ## ... Similar to previous best
    ## Run 96 stress 0.09314329 
    ## Run 97 stress 0.08121395 
    ## ... Procrustes: rmse 0.02026907  max resid 0.06577886 
    ## Run 98 stress 0.08085074 
    ## ... Procrustes: rmse 0.00959287  max resid 0.03304861 
    ## Run 99 stress 0.09120961 
    ## Run 100 stress 0.09124913 
    ## Run 101 stress 0.09457621 
    ## Run 102 stress 0.09124916 
    ## Run 103 stress 0.1047493 
    ## Run 104 stress 0.09124919 
    ## Run 105 stress 0.09283159 
    ## Run 106 stress 0.08121401 
    ## ... Procrustes: rmse 0.02028835  max resid 0.06603463 
    ## Run 107 stress 0.09178831 
    ## Run 108 stress 0.09129245 
    ## Run 109 stress 0.1014979 
    ## Run 110 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027694  max resid 0.06582567 
    ## Run 111 stress 0.081214 
    ## ... Procrustes: rmse 0.02031555  max resid 0.06608514 
    ## Run 112 stress 0.09178819 
    ## Run 113 stress 0.09178818 
    ## Run 114 stress 0.106788 
    ## Run 115 stress 0.08121396 
    ## ... Procrustes: rmse 0.02028248  max resid 0.06589005 
    ## Run 116 stress 0.09124933 
    ## Run 117 stress 0.1062402 
    ## Run 118 stress 0.09457589 
    ## Run 119 stress 0.09145239 
    ## Run 120 stress 0.09120968 
    ## Run 121 stress 0.08085074 
    ## ... Procrustes: rmse 0.009594508  max resid 0.0330581 
    ## Run 122 stress 0.08152292 
    ## Run 123 stress 0.0931432 
    ## Run 124 stress 0.09559565 
    ## Run 125 stress 0.0912097 
    ## Run 126 stress 0.09236464 
    ## Run 127 stress 0.0945763 
    ## Run 128 stress 0.09563364 
    ## Run 129 stress 0.08085073 
    ## ... Procrustes: rmse 0.009599131  max resid 0.03308482 
    ## Run 130 stress 0.08152297 
    ## Run 131 stress 0.08152298 
    ## Run 132 stress 0.08085074 
    ## ... Procrustes: rmse 0.009597603  max resid 0.03307086 
    ## Run 133 stress 0.08121395 
    ## ... Procrustes: rmse 0.02026284  max resid 0.06574451 
    ## Run 134 stress 0.0912096 
    ## Run 135 stress 0.09124911 
    ## Run 136 stress 0.08080643 
    ## ... Procrustes: rmse 4.847233e-05  max resid 0.0001744813 
    ## ... Similar to previous best
    ## Run 137 stress 0.09120962 
    ## Run 138 stress 0.08131825 
    ## Run 139 stress 0.08085077 
    ## ... Procrustes: rmse 0.009580996  max resid 0.03297681 
    ## Run 140 stress 0.08152292 
    ## Run 141 stress 0.09178832 
    ## Run 142 stress 0.08131812 
    ## Run 143 stress 0.08085077 
    ## ... Procrustes: rmse 0.009577579  max resid 0.03295722 
    ## Run 144 stress 0.08080642 
    ## ... Procrustes: rmse 2.268397e-05  max resid 8.071924e-05 
    ## ... Similar to previous best
    ## Run 145 stress 0.08121397 
    ## ... Procrustes: rmse 0.02026263  max resid 0.06575943 
    ## Run 146 stress 0.09236467 
    ## Run 147 stress 0.09563361 
    ## Run 148 stress 0.1047357 
    ## Run 149 stress 0.3624763 
    ## Run 150 stress 0.09120963 
    ## Run 151 stress 0.08121399 
    ## ... Procrustes: rmse 0.02029666  max resid 0.06602439 
    ## Run 152 stress 0.09296729 
    ## Run 153 stress 0.0912097 
    ## Run 154 stress 0.1060198 
    ## Run 155 stress 0.08080641 
    ## ... New best solution
    ## ... Procrustes: rmse 6.183876e-05  max resid 0.0002252786 
    ## ... Similar to previous best
    ## Run 156 stress 0.08085074 
    ## ... Procrustes: rmse 0.009599389  max resid 0.0330941 
    ## Run 157 stress 0.09124911 
    ## Run 158 stress 0.08085077 
    ## ... Procrustes: rmse 0.00958215  max resid 0.03297564 
    ## Run 159 stress 0.1010535 
    ## Run 160 stress 0.09559562 
    ## Run 161 stress 0.08152292 
    ## Run 162 stress 0.1029169 
    ## Run 163 stress 0.0912096 
    ## Run 164 stress 0.09178815 
    ## Run 165 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027277  max resid 0.0657391 
    ## Run 166 stress 0.09129238 
    ## Run 167 stress 0.08080648 
    ## ... Procrustes: rmse 0.0001054119  max resid 0.0003780029 
    ## ... Similar to previous best
    ## Run 168 stress 0.1028979 
    ## Run 169 stress 0.09270889 
    ## Run 170 stress 0.08080641 
    ## ... New best solution
    ## ... Procrustes: rmse 4.885378e-05  max resid 0.0001816786 
    ## ... Similar to previous best
    ## Run 171 stress 0.08085074 
    ## ... Procrustes: rmse 0.00959802  max resid 0.03308574 
    ## Run 172 stress 0.08085074 
    ## ... Procrustes: rmse 0.009600111  max resid 0.03309476 
    ## Run 173 stress 0.08085082 
    ## ... Procrustes: rmse 0.009613926  max resid 0.03314106 
    ## Run 174 stress 0.09178822 
    ## Run 175 stress 0.08080645 
    ## ... Procrustes: rmse 7.409517e-05  max resid 0.0002413117 
    ## ... Similar to previous best
    ## Run 176 stress 0.1062566 
    ## Run 177 stress 0.09270877 
    ## Run 178 stress 0.08121396 
    ## ... Procrustes: rmse 0.02028784  max resid 0.06590524 
    ## Run 179 stress 0.08131805 
    ## Run 180 stress 0.0917882 
    ## Run 181 stress 0.0928316 
    ## Run 182 stress 0.09180031 
    ## Run 183 stress 0.1051522 
    ## Run 184 stress 0.09387652 
    ## Run 185 stress 0.09236472 
    ## Run 186 stress 0.08085083 
    ## ... Procrustes: rmse 0.009614981  max resid 0.03313896 
    ## Run 187 stress 0.09180031 
    ## Run 188 stress 0.09563386 
    ## Run 189 stress 0.09129239 
    ## Run 190 stress 0.105895 
    ## Run 191 stress 0.08121395 
    ## ... Procrustes: rmse 0.02026378  max resid 0.06572029 
    ## Run 192 stress 0.081214 
    ## ... Procrustes: rmse 0.02027446  max resid 0.06587543 
    ## Run 193 stress 0.08152294 
    ## Run 194 stress 0.1045884 
    ## Run 195 stress 0.0808508 
    ## ... Procrustes: rmse 0.00957123  max resid 0.03291861 
    ## Run 196 stress 0.08085074 
    ## ... Procrustes: rmse 0.009595704  max resid 0.03306057 
    ## Run 197 stress 0.08121398 
    ## ... Procrustes: rmse 0.02027855  max resid 0.06583766 
    ## Run 198 stress 0.08131813 
    ## Run 199 stress 0.09145246 
    ## Run 200 stress 0.08121397 
    ## ... Procrustes: rmse 0.02024351  max resid 0.06562432 
    ## Run 201 stress 0.0927088 
    ## Run 202 stress 0.08080658 
    ## ... Procrustes: rmse 0.0001642918  max resid 0.0005796359 
    ## ... Similar to previous best
    ## Run 203 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027102  max resid 0.06578185 
    ## Run 204 stress 0.08131822 
    ## Run 205 stress 0.08152291 
    ## Run 206 stress 0.1068324 
    ## Run 207 stress 0.09296739 
    ## Run 208 stress 0.09283153 
    ## Run 209 stress 0.08121396 
    ## ... Procrustes: rmse 0.02026051  max resid 0.06572909 
    ## Run 210 stress 0.09124912 
    ## Run 211 stress 0.09129238 
    ## Run 212 stress 0.08080642 
    ## ... Procrustes: rmse 6.942679e-05  max resid 0.0002553938 
    ## ... Similar to previous best
    ## Run 213 stress 0.08131798 
    ## Run 214 stress 0.09236467 
    ## Run 215 stress 0.09124911 
    ## Run 216 stress 0.08131814 
    ## Run 217 stress 0.1055546 
    ## Run 218 stress 0.08121397 
    ## ... Procrustes: rmse 0.02029957  max resid 0.06598134 
    ## Run 219 stress 0.08085074 
    ## ... Procrustes: rmse 0.009604553  max resid 0.03310853 
    ## Run 220 stress 0.08121396 
    ## ... Procrustes: rmse 0.02028679  max resid 0.06591179 
    ## Run 221 stress 0.08131798 
    ## Run 222 stress 0.09124913 
    ## Run 223 stress 0.08121396 
    ## ... Procrustes: rmse 0.02027562  max resid 0.06581259 
    ## Run 224 stress 0.09270896 
    ## Run 225 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027974  max resid 0.06583933 
    ## Run 226 stress 0.09236469 
    ## Run 227 stress 0.08131806 
    ## Run 228 stress 0.0912493 
    ## Run 229 stress 0.08131804 
    ## Run 230 stress 0.08085078 
    ## ... Procrustes: rmse 0.00958774  max resid 0.03301542 
    ## Run 231 stress 0.09129239 
    ## Run 232 stress 0.08121396 
    ## ... Procrustes: rmse 0.0202646  max resid 0.06578336 
    ## Run 233 stress 0.1068383 
    ## Run 234 stress 0.09180031 
    ## Run 235 stress 0.08085083 
    ## ... Procrustes: rmse 0.009568039  max resid 0.03289595 
    ## Run 236 stress 0.08080645 
    ## ... Procrustes: rmse 7.704894e-05  max resid 0.000273863 
    ## ... Similar to previous best
    ## Run 237 stress 0.08121398 
    ## ... Procrustes: rmse 0.02025191  max resid 0.06565754 
    ## Run 238 stress 0.08080651 
    ## ... Procrustes: rmse 0.0001402432  max resid 0.0004978439 
    ## ... Similar to previous best
    ## Run 239 stress 0.1085157 
    ## Run 240 stress 0.08152293 
    ## Run 241 stress 0.09178819 
    ## Run 242 stress 0.08152292 
    ## Run 243 stress 0.08121395 
    ## ... Procrustes: rmse 0.0202763  max resid 0.06584078 
    ## Run 244 stress 0.09236464 
    ## Run 245 stress 0.0912097 
    ## Run 246 stress 0.08121397 
    ## ... Procrustes: rmse 0.02026161  max resid 0.06575182 
    ## Run 247 stress 0.08080647 
    ## ... Procrustes: rmse 0.0001043332  max resid 0.0003726909 
    ## ... Similar to previous best
    ## Run 248 stress 0.08152291 
    ## Run 249 stress 0.08152291 
    ## Run 250 stress 0.1048181 
    ## Run 251 stress 0.102716 
    ## Run 252 stress 0.09457613 
    ## Run 253 stress 0.09296718 
    ## Run 254 stress 0.09283155 
    ## Run 255 stress 0.08085089 
    ## ... Procrustes: rmse 0.00958799  max resid 0.03304028 
    ## Run 256 stress 0.08131808 
    ## Run 257 stress 0.0812142 
    ## ... Procrustes: rmse 0.02027584  max resid 0.06579424 
    ## Run 258 stress 0.08152292 
    ## Run 259 stress 0.09178823 
    ## Run 260 stress 0.08085083 
    ## ... Procrustes: rmse 0.009589615  max resid 0.03303141 
    ## Run 261 stress 0.08121398 
    ## ... Procrustes: rmse 0.02025053  max resid 0.06574302 
    ## Run 262 stress 0.09236465 
    ## Run 263 stress 0.08080642 
    ## ... Procrustes: rmse 2.751387e-05  max resid 9.259843e-05 
    ## ... Similar to previous best
    ## Run 264 stress 0.09270885 
    ## Run 265 stress 0.08121396 
    ## ... Procrustes: rmse 0.02029164  max resid 0.06588119 
    ## Run 266 stress 0.09178818 
    ## Run 267 stress 0.08152296 
    ## Run 268 stress 0.08085074 
    ## ... Procrustes: rmse 0.009603525  max resid 0.03310417 
    ## Run 269 stress 0.09124912 
    ## Run 270 stress 0.08085073 
    ## ... Procrustes: rmse 0.009601653  max resid 0.03309907 
    ## Run 271 stress 0.08121397 
    ## ... Procrustes: rmse 0.02027789  max resid 0.06581041 
    ## Run 272 stress 0.08121395 
    ## ... Procrustes: rmse 0.0202679  max resid 0.06577624 
    ## Run 273 stress 0.1053118 
    ## Run 274 stress 0.08085077 
    ## ... Procrustes: rmse 0.009579478  max resid 0.03296831 
    ## Run 275 stress 0.08085073 
    ## ... Procrustes: rmse 0.009602807  max resid 0.0331068 
    ## Run 276 stress 0.09314346 
    ## Run 277 stress 0.08121395 
    ## ... Procrustes: rmse 0.0202745  max resid 0.06581093 
    ## Run 278 stress 0.08080648 
    ## ... Procrustes: rmse 0.0001157801  max resid 0.0004030443 
    ## ... Similar to previous best
    ## Run 279 stress 0.08121396 
    ## ... Procrustes: rmse 0.02025844  max resid 0.06567173 
    ## Run 280 stress 0.08152301 
    ## Run 281 stress 0.08131801 
    ## Run 282 stress 0.08121396 
    ## ... Procrustes: rmse 0.02027602  max resid 0.06583764 
    ## Run 283 stress 0.09120962 
    ## Run 284 stress 0.0926121 
    ## Run 285 stress 0.0955957 
    ## Run 286 stress 0.09383882 
    ## Run 287 stress 0.0813179 
    ## Run 288 stress 0.09120967 
    ## Run 289 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027187  max resid 0.06575534 
    ## Run 290 stress 0.09270876 
    ## Run 291 stress 0.08085073 
    ## ... Procrustes: rmse 0.009601577  max resid 0.03309791 
    ## Run 292 stress 0.1035563 
    ## Run 293 stress 0.08080642 
    ## ... Procrustes: rmse 6.148374e-05  max resid 0.0002248574 
    ## ... Similar to previous best
    ## Run 294 stress 0.08121399 
    ## ... Procrustes: rmse 0.020295  max resid 0.06600615 
    ## Run 295 stress 0.081214 
    ## ... Procrustes: rmse 0.02024401  max resid 0.06555796 
    ## Run 296 stress 0.09383869 
    ## Run 297 stress 0.08152296 
    ## Run 298 stress 0.09180031 
    ## Run 299 stress 0.08152291 
    ## Run 300 stress 0.1007649 
    ## Run 301 stress 0.0808508 
    ## ... Procrustes: rmse 0.009572891  max resid 0.03292868 
    ## Run 302 stress 0.08080644 
    ## ... Procrustes: rmse 6.230945e-05  max resid 0.0002201863 
    ## ... Similar to previous best
    ## Run 303 stress 0.09124912 
    ## Run 304 stress 0.08152291 
    ## Run 305 stress 0.09314319 
    ## Run 306 stress 0.081523 
    ## Run 307 stress 0.1051646 
    ## Run 308 stress 0.09383872 
    ## Run 309 stress 0.08085082 
    ## ... Procrustes: rmse 0.0095916  max resid 0.0330579 
    ## Run 310 stress 0.08131812 
    ## Run 311 stress 0.08085074 
    ## ... Procrustes: rmse 0.009592596  max resid 0.03304582 
    ## Run 312 stress 0.09457656 
    ## Run 313 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027274  max resid 0.06576199 
    ## Run 314 stress 0.08152292 
    ## Run 315 stress 0.08131816 
    ## Run 316 stress 0.08131802 
    ## Run 317 stress 0.09129249 
    ## Run 318 stress 0.08080641 
    ## ... Procrustes: rmse 6.445564e-05  max resid 0.0002356737 
    ## ... Similar to previous best
    ## Run 319 stress 0.08131791 
    ## Run 320 stress 0.08131801 
    ## Run 321 stress 0.09296741 
    ## Run 322 stress 0.08152299 
    ## Run 323 stress 0.08121398 
    ## ... Procrustes: rmse 0.02029465  max resid 0.0659837 
    ## Run 324 stress 0.08121399 
    ## ... Procrustes: rmse 0.02024425  max resid 0.06554754 
    ## Run 325 stress 0.08080648 
    ## ... Procrustes: rmse 9.760058e-05  max resid 0.0003445629 
    ## ... Similar to previous best
    ## Run 326 stress 0.08080641 
    ## ... Procrustes: rmse 1.280223e-05  max resid 3.483031e-05 
    ## ... Similar to previous best
    ## Run 327 stress 0.08152291 
    ## Run 328 stress 0.08121396 
    ## ... Procrustes: rmse 0.02025879  max resid 0.06575846 
    ## Run 329 stress 0.08121398 
    ## ... Procrustes: rmse 0.02029432  max resid 0.06597519 
    ## Run 330 stress 0.09180037 
    ## Run 331 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027685  max resid 0.06581146 
    ## Run 332 stress 0.1035431 
    ## Run 333 stress 0.1060473 
    ## Run 334 stress 0.08152294 
    ## Run 335 stress 0.08080641 
    ## ... New best solution
    ## ... Procrustes: rmse 1.550404e-05  max resid 3.892401e-05 
    ## ... Similar to previous best
    ## Run 336 stress 0.08131809 
    ## Run 337 stress 0.08085076 
    ## ... Procrustes: rmse 0.009583144  max resid 0.03298871 
    ## Run 338 stress 0.1042223 
    ## Run 339 stress 0.08121399 
    ## ... Procrustes: rmse 0.0202935  max resid 0.06597995 
    ## Run 340 stress 0.09261232 
    ## Run 341 stress 0.08121399 
    ## ... Procrustes: rmse 0.02031638  max resid 0.0660456 
    ## Run 342 stress 0.09145237 
    ## Run 343 stress 0.09120962 
    ## Run 344 stress 0.08121397 
    ## ... Procrustes: rmse 0.02028455  max resid 0.06591475 
    ## Run 345 stress 0.08152292 
    ## Run 346 stress 0.08080651 
    ## ... Procrustes: rmse 0.0001745823  max resid 0.0006306853 
    ## ... Similar to previous best
    ## Run 347 stress 0.08152294 
    ## Run 348 stress 0.09261238 
    ## Run 349 stress 0.08121396 
    ## ... Procrustes: rmse 0.02029557  max resid 0.06586167 
    ## Run 350 stress 0.09124911 
    ## Run 351 stress 0.0938387 
    ## Run 352 stress 0.09457625 
    ## Run 353 stress 0.09145276 
    ## Run 354 stress 0.08085077 
    ## ... Procrustes: rmse 0.009593451  max resid 0.03303355 
    ## Run 355 stress 0.08152295 
    ## Run 356 stress 0.08085074 
    ## ... Procrustes: rmse 0.009594375  max resid 0.0330554 
    ## Run 357 stress 0.08131799 
    ## Run 358 stress 0.1030879 
    ## Run 359 stress 0.08080649 
    ## ... Procrustes: rmse 0.0001599282  max resid 0.0005787823 
    ## ... Similar to previous best
    ## Run 360 stress 0.0808508 
    ## ... Procrustes: rmse 0.009576227  max resid 0.03294775 
    ## Run 361 stress 0.08085077 
    ## ... Procrustes: rmse 0.009580374  max resid 0.03297309 
    ## Run 362 stress 0.09403915 
    ## Run 363 stress 0.0912492 
    ## Run 364 stress 0.08131811 
    ## Run 365 stress 0.08121395 
    ## ... Procrustes: rmse 0.02028437  max resid 0.06582767 
    ## Run 366 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027638  max resid 0.06579447 
    ## Run 367 stress 0.08085085 
    ## ... Procrustes: rmse 0.009563344  max resid 0.03287196 
    ## Run 368 stress 0.09563402 
    ## Run 369 stress 0.08152292 
    ## Run 370 stress 0.09314317 
    ## Run 371 stress 0.09178816 
    ## Run 372 stress 0.08085082 
    ## ... Procrustes: rmse 0.009601703  max resid 0.03306903 
    ## Run 373 stress 0.08121396 
    ## ... Procrustes: rmse 0.02029858  max resid 0.06592587 
    ## Run 374 stress 0.102123 
    ## Run 375 stress 0.0912924 
    ## Run 376 stress 0.09403877 
    ## Run 377 stress 0.08080649 
    ## ... Procrustes: rmse 0.000133392  max resid 0.0004812603 
    ## ... Similar to previous best
    ## Run 378 stress 0.08085073 
    ## ... Procrustes: rmse 0.009601156  max resid 0.03309784 
    ## Run 379 stress 0.1007647 
    ## Run 380 stress 0.09178815 
    ## Run 381 stress 0.09124911 
    ## Run 382 stress 0.102716 
    ## Run 383 stress 0.08152291 
    ## Run 384 stress 0.09457615 
    ## Run 385 stress 0.08121397 
    ## ... Procrustes: rmse 0.02029034  max resid 0.06592176 
    ## Run 386 stress 0.09283157 
    ## Run 387 stress 0.09180031 
    ## Run 388 stress 0.08121396 
    ## ... Procrustes: rmse 0.02025907  max resid 0.0656419 
    ## Run 389 stress 0.08085076 
    ## ... Procrustes: rmse 0.009582755  max resid 0.03298449 
    ## Run 390 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027621  max resid 0.06577265 
    ## Run 391 stress 0.08080643 
    ## ... Procrustes: rmse 6.15719e-05  max resid 0.0002157831 
    ## ... Similar to previous best
    ## Run 392 stress 0.08152296 
    ## Run 393 stress 0.08121396 
    ## ... Procrustes: rmse 0.02026131  max resid 0.06567035 
    ## Run 394 stress 0.09178819 
    ## Run 395 stress 0.08085079 
    ## ... Procrustes: rmse 0.00959713  max resid 0.0330866 
    ## Run 396 stress 0.08152292 
    ## Run 397 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027177  max resid 0.06574853 
    ## Run 398 stress 0.08121395 
    ## ... Procrustes: rmse 0.02026511  max resid 0.06570765 
    ## Run 399 stress 0.08080652 
    ## ... Procrustes: rmse 0.0001801389  max resid 0.0006503052 
    ## ... Similar to previous best
    ## Run 400 stress 0.0912924 
    ## Run 401 stress 0.08121404 
    ## ... Procrustes: rmse 0.02025062  max resid 0.06578178 
    ## Run 402 stress 0.09180031 
    ## Run 403 stress 0.09129247 
    ## Run 404 stress 0.09383883 
    ## Run 405 stress 0.08080641 
    ## ... Procrustes: rmse 2.168682e-05  max resid 4.680259e-05 
    ## ... Similar to previous best
    ## Run 406 stress 0.09559583 
    ## Run 407 stress 0.08080641 
    ## ... Procrustes: rmse 1.237937e-05  max resid 2.900641e-05 
    ## ... Similar to previous best
    ## Run 408 stress 0.09261213 
    ## Run 409 stress 0.08152291 
    ## Run 410 stress 0.09261222 
    ## Run 411 stress 0.09261242 
    ## Run 412 stress 0.09129238 
    ## Run 413 stress 0.08152294 
    ## Run 414 stress 0.09314339 
    ## Run 415 stress 0.08121404 
    ## ... Procrustes: rmse 0.02034164  max resid 0.06618281 
    ## Run 416 stress 0.09124925 
    ## Run 417 stress 0.09314324 
    ## Run 418 stress 0.1007647 
    ## Run 419 stress 0.09387661 
    ## Run 420 stress 0.08131803 
    ## Run 421 stress 0.09178817 
    ## Run 422 stress 0.08131807 
    ## Run 423 stress 0.08121395 
    ## ... Procrustes: rmse 0.02026941  max resid 0.06573652 
    ## Run 424 stress 0.08085076 
    ## ... Procrustes: rmse 0.009606095  max resid 0.03311476 
    ## Run 425 stress 0.09236468 
    ## Run 426 stress 0.08121396 
    ## ... Procrustes: rmse 0.02027477  max resid 0.06572398 
    ## Run 427 stress 0.08152291 
    ## Run 428 stress 0.08131803 
    ## Run 429 stress 0.0938388 
    ## Run 430 stress 0.08121397 
    ## ... Procrustes: rmse 0.02030453  max resid 0.0659683 
    ## Run 431 stress 0.09236466 
    ## Run 432 stress 0.09124913 
    ## Run 433 stress 0.08121395 
    ## ... Procrustes: rmse 0.02026463  max resid 0.06569174 
    ## Run 434 stress 0.08121396 
    ## ... Procrustes: rmse 0.02025624  max resid 0.06563289 
    ## Run 435 stress 0.09283172 
    ## Run 436 stress 0.09403873 
    ## Run 437 stress 0.08152298 
    ## Run 438 stress 0.09129238 
    ## Run 439 stress 0.08085074 
    ## ... Procrustes: rmse 0.009602548  max resid 0.03309971 
    ## Run 440 stress 0.08080643 
    ## ... Procrustes: rmse 6.175442e-05  max resid 0.000218634 
    ## ... Similar to previous best
    ## Run 441 stress 0.08121395 
    ## ... Procrustes: rmse 0.02026146  max resid 0.0656765 
    ## Run 442 stress 0.08121396 
    ## ... Procrustes: rmse 0.02028024  max resid 0.06585296 
    ## Run 443 stress 0.09120961 
    ## Run 444 stress 0.08121395 
    ## ... Procrustes: rmse 0.02026894  max resid 0.06572382 
    ## Run 445 stress 0.081523 
    ## Run 446 stress 0.09477979 
    ## Run 447 stress 0.08085073 
    ## ... Procrustes: rmse 0.009599787  max resid 0.03308764 
    ## Run 448 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027172  max resid 0.06576315 
    ## Run 449 stress 0.1068323 
    ## Run 450 stress 0.08121395 
    ## ... Procrustes: rmse 0.0202649  max resid 0.06574273 
    ## Run 451 stress 0.09283168 
    ## Run 452 stress 0.08085074 
    ## ... Procrustes: rmse 0.009603423  max resid 0.03310199 
    ## Run 453 stress 0.09120968 
    ## Run 454 stress 0.09383872 
    ## Run 455 stress 0.09178819 
    ## Run 456 stress 0.09124912 
    ## Run 457 stress 0.08085075 
    ## ... Procrustes: rmse 0.009597532  max resid 0.03308632 
    ## Run 458 stress 0.09178829 
    ## Run 459 stress 0.08121396 
    ## ... Procrustes: rmse 0.02026335  max resid 0.06568068 
    ## Run 460 stress 0.09283161 
    ## Run 461 stress 0.1045585 
    ## Run 462 stress 0.08121396 
    ## ... Procrustes: rmse 0.02028455  max resid 0.06586908 
    ## Run 463 stress 0.08121395 
    ## ... Procrustes: rmse 0.02026863  max resid 0.06572442 
    ## Run 464 stress 0.08085076 
    ## ... Procrustes: rmse 0.009609016  max resid 0.03312621 
    ## Run 465 stress 0.09178815 
    ## Run 466 stress 0.09180037 
    ## Run 467 stress 0.08085077 
    ## ... Procrustes: rmse 0.009608746  max resid 0.03311875 
    ## Run 468 stress 0.08080651 
    ## ... Procrustes: rmse 0.0001461854  max resid 0.0005191291 
    ## ... Similar to previous best
    ## Run 469 stress 0.09178818 
    ## Run 470 stress 0.09120964 
    ## Run 471 stress 0.08121399 
    ## ... Procrustes: rmse 0.02030137  max resid 0.06600216 
    ## Run 472 stress 0.0808508 
    ## ... Procrustes: rmse 0.009589791  max resid 0.03304585 
    ## Run 473 stress 0.09457643 
    ## Run 474 stress 0.0912096 
    ## Run 475 stress 0.105976 
    ## Run 476 stress 0.08085081 
    ## ... Procrustes: rmse 0.009597731  max resid 0.03308931 
    ## Run 477 stress 0.09178827 
    ## Run 478 stress 0.08152291 
    ## Run 479 stress 0.09178819 
    ## Run 480 stress 0.09124911 
    ## Run 481 stress 0.09261211 
    ## Run 482 stress 0.08152291 
    ## Run 483 stress 0.08085075 
    ## ... Procrustes: rmse 0.00960466  max resid 0.0331106 
    ## Run 484 stress 0.08121395 
    ## ... Procrustes: rmse 0.0202784  max resid 0.06580277 
    ## Run 485 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027127  max resid 0.06572964 
    ## Run 486 stress 0.08121397 
    ## ... Procrustes: rmse 0.0202923  max resid 0.06591875 
    ## Run 487 stress 0.08080643 
    ## ... Procrustes: rmse 9.630668e-05  max resid 0.0003484623 
    ## ... Similar to previous best
    ## Run 488 stress 0.09178821 
    ## Run 489 stress 0.1051521 
    ## Run 490 stress 0.09124915 
    ## Run 491 stress 0.08121396 
    ## ... Procrustes: rmse 0.02026039  max resid 0.065737 
    ## Run 492 stress 0.08121395 
    ## ... Procrustes: rmse 0.02027915  max resid 0.0657926 
    ## Run 493 stress 0.09236465 
    ## Run 494 stress 0.09563372 
    ## Run 495 stress 0.09270876 
    ## Run 496 stress 0.08085076 
    ## ... Procrustes: rmse 0.009595646  max resid 0.03307326 
    ## Run 497 stress 0.09477998 
    ## Run 498 stress 0.08080642 
    ## ... Procrustes: rmse 3.801439e-05  max resid 0.0001266613 
    ## ... Similar to previous best
    ## Run 499 stress 0.08085075 
    ## ... Procrustes: rmse 0.009591346  max resid 0.03303651 
    ## Run 500 stress 0.09387651 
    ## *** Best solution repeated 12 times

``` r
# Mixed and stratified lakes
PD_beta_geo_MS_NMDS <- metaMDS(PD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.04845889 
    ## Run 1 stress 0.0620353 
    ## Run 2 stress 0.05333362 
    ## Run 3 stress 0.059333 
    ## Run 4 stress 0.06317301 
    ## Run 5 stress 0.05992872 
    ## Run 6 stress 0.04845894 
    ## ... Procrustes: rmse 0.0001076822  max resid 0.0001867216 
    ## ... Similar to previous best
    ## Run 7 stress 0.06036745 
    ## Run 8 stress 0.05401009 
    ## Run 9 stress 0.04875109 
    ## ... Procrustes: rmse 0.04119801  max resid 0.138538 
    ## Run 10 stress 0.05333361 
    ## Run 11 stress 0.04875108 
    ## ... Procrustes: rmse 0.04108895  max resid 0.1382016 
    ## Run 12 stress 0.05401012 
    ## Run 13 stress 0.06614146 
    ## Run 14 stress 0.06850977 
    ## Run 15 stress 0.04845892 
    ## ... Procrustes: rmse 9.761443e-05  max resid 0.0001454544 
    ## ... Similar to previous best
    ## Run 16 stress 0.06617394 
    ## Run 17 stress 0.05333361 
    ## Run 18 stress 0.06096275 
    ## Run 19 stress 0.04875115 
    ## ... Procrustes: rmse 0.04105599  max resid 0.1380999 
    ## Run 20 stress 0.04845891 
    ## ... Procrustes: rmse 4.389065e-05  max resid 7.254226e-05 
    ## ... Similar to previous best
    ## Run 21 stress 0.05118024 
    ## Run 22 stress 0.07023414 
    ## Run 23 stress 0.05849524 
    ## Run 24 stress 0.04875113 
    ## ... Procrustes: rmse 0.04121195  max resid 0.1385809 
    ## Run 25 stress 0.05456969 
    ## Run 26 stress 0.06096274 
    ## Run 27 stress 0.05118021 
    ## Run 28 stress 0.05333364 
    ## Run 29 stress 0.06036745 
    ## Run 30 stress 0.04845891 
    ## ... Procrustes: rmse 4.783896e-05  max resid 7.217234e-05 
    ## ... Similar to previous best
    ## Run 31 stress 0.07117142 
    ## Run 32 stress 0.04845892 
    ## ... Procrustes: rmse 6.010066e-05  max resid 9.27492e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.04845894 
    ## ... Procrustes: rmse 0.0001215333  max resid 0.000203761 
    ## ... Similar to previous best
    ## Run 34 stress 0.04875117 
    ## ... Procrustes: rmse 0.04105028  max resid 0.1380817 
    ## Run 35 stress 0.04845889 
    ## ... New best solution
    ## ... Procrustes: rmse 1.158696e-05  max resid 1.686677e-05 
    ## ... Similar to previous best
    ## Run 36 stress 0.05348165 
    ## Run 37 stress 0.05118027 
    ## Run 38 stress 0.06995447 
    ## Run 39 stress 0.04875107 
    ## ... Procrustes: rmse 0.04116267  max resid 0.1384289 
    ## Run 40 stress 0.06086811 
    ## Run 41 stress 0.05992866 
    ## Run 42 stress 0.04845896 
    ## ... Procrustes: rmse 0.000119371  max resid 0.0001815917 
    ## ... Similar to previous best
    ## Run 43 stress 0.04845892 
    ## ... Procrustes: rmse 7.805865e-05  max resid 0.0001258381 
    ## ... Similar to previous best
    ## Run 44 stress 0.06814641 
    ## Run 45 stress 0.05991816 
    ## Run 46 stress 0.3473218 
    ## Run 47 stress 0.06036746 
    ## Run 48 stress 0.05865982 
    ## Run 49 stress 0.06086822 
    ## Run 50 stress 0.04845895 
    ## ... Procrustes: rmse 0.0001231049  max resid 0.0002151627 
    ## ... Similar to previous best
    ## Run 51 stress 0.0484589 
    ## ... Procrustes: rmse 2.130774e-05  max resid 3.760733e-05 
    ## ... Similar to previous best
    ## Run 52 stress 0.04845895 
    ## ... Procrustes: rmse 0.0001150977  max resid 0.0002070238 
    ## ... Similar to previous best
    ## Run 53 stress 0.05933306 
    ## Run 54 stress 0.05118022 
    ## Run 55 stress 0.05667836 
    ## Run 56 stress 0.05348158 
    ## Run 57 stress 0.05748881 
    ## Run 58 stress 0.0534816 
    ## Run 59 stress 0.05456978 
    ## Run 60 stress 0.0593314 
    ## Run 61 stress 0.0534817 
    ## Run 62 stress 0.04845889 
    ## ... Procrustes: rmse 2.407104e-05  max resid 4.781742e-05 
    ## ... Similar to previous best
    ## Run 63 stress 0.3425938 
    ## Run 64 stress 0.0545697 
    ## Run 65 stress 0.06041479 
    ## Run 66 stress 0.05333374 
    ## Run 67 stress 0.06236433 
    ## Run 68 stress 0.04875107 
    ## ... Procrustes: rmse 0.04110101  max resid 0.138239 
    ## Run 69 stress 0.06203521 
    ## Run 70 stress 0.05849522 
    ## Run 71 stress 0.06649432 
    ## Run 72 stress 0.05748878 
    ## Run 73 stress 0.05865992 
    ## Run 74 stress 0.04845891 
    ## ... Procrustes: rmse 6.181313e-05  max resid 9.229162e-05 
    ## ... Similar to previous best
    ## Run 75 stress 0.05118021 
    ## Run 76 stress 0.05748877 
    ## Run 77 stress 0.06622334 
    ## Run 78 stress 0.05333367 
    ## Run 79 stress 0.07117148 
    ## Run 80 stress 0.06096276 
    ## Run 81 stress 0.04875119 
    ## ... Procrustes: rmse 0.04105039  max resid 0.138081 
    ## Run 82 stress 0.06086812 
    ## Run 83 stress 0.04845892 
    ## ... Procrustes: rmse 6.725488e-05  max resid 9.791472e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.06036747 
    ## Run 85 stress 0.05401 
    ## Run 86 stress 0.06086822 
    ## Run 87 stress 0.05348175 
    ## Run 88 stress 0.04845891 
    ## ... Procrustes: rmse 5.784424e-05  max resid 9.038743e-05 
    ## ... Similar to previous best
    ## Run 89 stress 0.05991817 
    ## Run 90 stress 0.0662235 
    ## Run 91 stress 0.05400987 
    ## Run 92 stress 0.0623643 
    ## Run 93 stress 0.06036759 
    ## Run 94 stress 0.05849521 
    ## Run 95 stress 0.05992868 
    ## Run 96 stress 0.04875117 
    ## ... Procrustes: rmse 0.04122444  max resid 0.138618 
    ## Run 97 stress 0.0545698 
    ## Run 98 stress 0.05118022 
    ## Run 99 stress 0.04875104 
    ## ... Procrustes: rmse 0.04112453  max resid 0.1383112 
    ## Run 100 stress 0.05118023 
    ## Run 101 stress 0.04845893 
    ## ... Procrustes: rmse 9.387276e-05  max resid 0.0001381467 
    ## ... Similar to previous best
    ## Run 102 stress 0.05938472 
    ## Run 103 stress 0.0673981 
    ## Run 104 stress 0.06614167 
    ## Run 105 stress 0.0593847 
    ## Run 106 stress 0.06622327 
    ## Run 107 stress 0.05748877 
    ## Run 108 stress 0.07023432 
    ## Run 109 stress 0.06774299 
    ## Run 110 stress 0.05118022 
    ## Run 111 stress 0.04845897 
    ## ... Procrustes: rmse 0.0001252121  max resid 0.0001835761 
    ## ... Similar to previous best
    ## Run 112 stress 0.0484589 
    ## ... Procrustes: rmse 5.346944e-05  max resid 8.02219e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.04875111 
    ## ... Procrustes: rmse 0.04107511  max resid 0.1381576 
    ## Run 114 stress 0.05933147 
    ## Run 115 stress 0.06236428 
    ## Run 116 stress 0.04845891 
    ## ... Procrustes: rmse 5.96128e-05  max resid 9.436548e-05 
    ## ... Similar to previous best
    ## Run 117 stress 0.05400997 
    ## Run 118 stress 0.04875109 
    ## ... Procrustes: rmse 0.04108977  max resid 0.1382044 
    ## Run 119 stress 0.06086807 
    ## Run 120 stress 0.05118022 
    ## Run 121 stress 0.05667835 
    ## Run 122 stress 0.05348169 
    ## Run 123 stress 0.05991824 
    ## Run 124 stress 0.05865985 
    ## Run 125 stress 0.05865985 
    ## Run 126 stress 0.3620352 
    ## Run 127 stress 0.0484589 
    ## ... Procrustes: rmse 1.209586e-05  max resid 2.445785e-05 
    ## ... Similar to previous best
    ## Run 128 stress 0.05456968 
    ## Run 129 stress 0.0484589 
    ## ... Procrustes: rmse 5.105892e-05  max resid 7.657928e-05 
    ## ... Similar to previous best
    ## Run 130 stress 0.04875112 
    ## ... Procrustes: rmse 0.04120451  max resid 0.1385576 
    ## Run 131 stress 0.05938478 
    ## Run 132 stress 0.04845897 
    ## ... Procrustes: rmse 0.0001155431  max resid 0.000167996 
    ## ... Similar to previous best
    ## Run 133 stress 0.05333369 
    ## Run 134 stress 0.05992877 
    ## Run 135 stress 0.0487511 
    ## ... Procrustes: rmse 0.04118614  max resid 0.1385007 
    ## Run 136 stress 0.07504564 
    ## Run 137 stress 0.05849515 
    ## Run 138 stress 0.06160425 
    ## Run 139 stress 0.04845896 
    ## ... Procrustes: rmse 0.0001081527  max resid 0.0001956765 
    ## ... Similar to previous best
    ## Run 140 stress 0.05933181 
    ## Run 141 stress 0.05938475 
    ## Run 142 stress 0.0661348 
    ## Run 143 stress 0.05333362 
    ## Run 144 stress 0.05456974 
    ## Run 145 stress 0.04875107 
    ## ... Procrustes: rmse 0.0410943  max resid 0.1382165 
    ## Run 146 stress 0.05992867 
    ## Run 147 stress 0.05933295 
    ## Run 148 stress 0.04875104 
    ## ... Procrustes: rmse 0.04111877  max resid 0.1382931 
    ## Run 149 stress 0.05865983 
    ## Run 150 stress 0.05118021 
    ## Run 151 stress 0.05401009 
    ## Run 152 stress 0.04845894 
    ## ... Procrustes: rmse 9.848241e-05  max resid 0.0001585608 
    ## ... Similar to previous best
    ## Run 153 stress 0.06036746 
    ## Run 154 stress 0.06236425 
    ## Run 155 stress 0.06036744 
    ## Run 156 stress 0.06236427 
    ## Run 157 stress 0.04845889 
    ## ... Procrustes: rmse 1.969029e-05  max resid 3.103345e-05 
    ## ... Similar to previous best
    ## Run 158 stress 0.04845891 
    ## ... Procrustes: rmse 2.482943e-05  max resid 3.366926e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.05667831 
    ## Run 160 stress 0.05667844 
    ## Run 161 stress 0.06036747 
    ## Run 162 stress 0.06851041 
    ## Run 163 stress 0.05933309 
    ## Run 164 stress 0.06036745 
    ## Run 165 stress 0.04845891 
    ## ... Procrustes: rmse 6.54975e-05  max resid 9.770814e-05 
    ## ... Similar to previous best
    ## Run 166 stress 0.04845889 
    ## ... Procrustes: rmse 9.752573e-06  max resid 1.42411e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.0533336 
    ## Run 168 stress 0.07023429 
    ## Run 169 stress 0.04875107 
    ## ... Procrustes: rmse 0.04117602  max resid 0.1384698 
    ## Run 170 stress 0.04875104 
    ## ... Procrustes: rmse 0.04114769  max resid 0.1383824 
    ## Run 171 stress 0.05118024 
    ## Run 172 stress 0.06617393 
    ## Run 173 stress 0.0484589 
    ## ... Procrustes: rmse 2.459597e-05  max resid 3.51089e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.04845891 
    ## ... Procrustes: rmse 5.489588e-05  max resid 9.159509e-05 
    ## ... Similar to previous best
    ## Run 175 stress 0.06851018 
    ## Run 176 stress 0.06624126 
    ## Run 177 stress 0.06096274 
    ## Run 178 stress 0.05400998 
    ## Run 179 stress 0.04845899 
    ## ... Procrustes: rmse 0.0001444234  max resid 0.0002623567 
    ## ... Similar to previous best
    ## Run 180 stress 0.06203526 
    ## Run 181 stress 0.05667813 
    ## Run 182 stress 0.06041486 
    ## Run 183 stress 0.362036 
    ## Run 184 stress 0.04845892 
    ## ... Procrustes: rmse 8.007665e-05  max resid 0.0001180342 
    ## ... Similar to previous best
    ## Run 185 stress 0.05938475 
    ## Run 186 stress 0.04875109 
    ## ... Procrustes: rmse 0.04118968  max resid 0.138512 
    ## Run 187 stress 0.04875117 
    ## ... Procrustes: rmse 0.04105239  max resid 0.1380868 
    ## Run 188 stress 0.05933312 
    ## Run 189 stress 0.04845898 
    ## ... Procrustes: rmse 0.0001413221  max resid 0.000256715 
    ## ... Similar to previous best
    ## Run 190 stress 0.04875111 
    ## ... Procrustes: rmse 0.0410717  max resid 0.1381476 
    ## Run 191 stress 0.05933139 
    ## Run 192 stress 0.05667823 
    ## Run 193 stress 0.05333369 
    ## Run 194 stress 0.04845889 
    ## ... Procrustes: rmse 2.458631e-05  max resid 3.710605e-05 
    ## ... Similar to previous best
    ## Run 195 stress 0.05933162 
    ## Run 196 stress 0.0599288 
    ## Run 197 stress 0.068311 
    ## Run 198 stress 0.05401002 
    ## Run 199 stress 0.06036746 
    ## Run 200 stress 0.06041472 
    ## Run 201 stress 0.05933123 
    ## Run 202 stress 0.04875105 
    ## ... Procrustes: rmse 0.04110628  max resid 0.1382546 
    ## Run 203 stress 0.06614132 
    ## Run 204 stress 0.06616889 
    ## Run 205 stress 0.05748878 
    ## Run 206 stress 0.06878466 
    ## Run 207 stress 0.05748884 
    ## Run 208 stress 0.04875114 
    ## ... Procrustes: rmse 0.04120202  max resid 0.1385504 
    ## Run 209 stress 0.05933321 
    ## Run 210 stress 0.05333362 
    ## Run 211 stress 0.04845894 
    ## ... Procrustes: rmse 9.813438e-05  max resid 0.00014679 
    ## ... Similar to previous best
    ## Run 212 stress 0.06624126 
    ## Run 213 stress 0.06845765 
    ## Run 214 stress 0.04875118 
    ## ... Procrustes: rmse 0.04104423  max resid 0.138062 
    ## Run 215 stress 0.04845899 
    ## ... Procrustes: rmse 0.0001317967  max resid 0.0002122613 
    ## ... Similar to previous best
    ## Run 216 stress 0.05333368 
    ## Run 217 stress 0.05118035 
    ## Run 218 stress 0.05865991 
    ## Run 219 stress 0.04845894 
    ## ... Procrustes: rmse 0.0001097435  max resid 0.000186368 
    ## ... Similar to previous best
    ## Run 220 stress 0.0593847 
    ## Run 221 stress 0.0484589 
    ## ... Procrustes: rmse 4.518594e-05  max resid 8.544735e-05 
    ## ... Similar to previous best
    ## Run 222 stress 0.05865988 
    ## Run 223 stress 0.04845896 
    ## ... Procrustes: rmse 0.0001285901  max resid 0.0002335218 
    ## ... Similar to previous best
    ## Run 224 stress 0.04845892 
    ## ... Procrustes: rmse 8.827313e-05  max resid 0.0001581337 
    ## ... Similar to previous best
    ## Run 225 stress 0.05456969 
    ## Run 226 stress 0.06617482 
    ## Run 227 stress 0.05748877 
    ## Run 228 stress 0.06831103 
    ## Run 229 stress 0.362035 
    ## Run 230 stress 0.05865987 
    ## Run 231 stress 0.0487511 
    ## ... Procrustes: rmse 0.04110115  max resid 0.1382401 
    ## Run 232 stress 0.05748879 
    ## Run 233 stress 0.05348166 
    ## Run 234 stress 0.05456976 
    ## Run 235 stress 0.05933164 
    ## Run 236 stress 0.04845894 
    ## ... Procrustes: rmse 0.000106332  max resid 0.0001910178 
    ## ... Similar to previous best
    ## Run 237 stress 0.06614159 
    ## Run 238 stress 0.05748879 
    ## Run 239 stress 0.06649437 
    ## Run 240 stress 0.05748878 
    ## Run 241 stress 0.05933295 
    ## Run 242 stress 0.06096274 
    ## Run 243 stress 0.0750454 
    ## Run 244 stress 0.05401 
    ## Run 245 stress 0.06036748 
    ## Run 246 stress 0.06203523 
    ## Run 247 stress 0.05400985 
    ## Run 248 stress 0.05933154 
    ## Run 249 stress 0.04875107 
    ## ... Procrustes: rmse 0.04110486  max resid 0.1382489 
    ## Run 250 stress 0.0593331 
    ## Run 251 stress 0.070195 
    ## Run 252 stress 0.05748887 
    ## Run 253 stress 0.04875115 
    ## ... Procrustes: rmse 0.04105685  max resid 0.138101 
    ## Run 254 stress 0.05748886 
    ## Run 255 stress 0.05400994 
    ## Run 256 stress 0.0659327 
    ## Run 257 stress 0.06785974 
    ## Run 258 stress 0.05401005 
    ## Run 259 stress 0.06845778 
    ## Run 260 stress 0.0484589 
    ## ... Procrustes: rmse 3.390474e-05  max resid 4.962104e-05 
    ## ... Similar to previous best
    ## Run 261 stress 0.05992866 
    ## Run 262 stress 0.05400998 
    ## Run 263 stress 0.05992878 
    ## Run 264 stress 0.05400985 
    ## Run 265 stress 0.04875107 
    ## ... Procrustes: rmse 0.04117608  max resid 0.1384699 
    ## Run 266 stress 0.04845893 
    ## ... Procrustes: rmse 8.674049e-05  max resid 0.0001397872 
    ## ... Similar to previous best
    ## Run 267 stress 0.04845893 
    ## ... Procrustes: rmse 9.040278e-05  max resid 0.0001446916 
    ## ... Similar to previous best
    ## Run 268 stress 0.04845901 
    ## ... Procrustes: rmse 3.457049e-05  max resid 8.184184e-05 
    ## ... Similar to previous best
    ## Run 269 stress 0.05933315 
    ## Run 270 stress 0.04875111 
    ## ... Procrustes: rmse 0.0410692  max resid 0.13814 
    ## Run 271 stress 0.05118029 
    ## Run 272 stress 0.06160427 
    ## Run 273 stress 0.06036755 
    ## Run 274 stress 0.05933157 
    ## Run 275 stress 0.0620352 
    ## Run 276 stress 0.05118034 
    ## Run 277 stress 0.0593847 
    ## Run 278 stress 0.05938473 
    ## Run 279 stress 0.05849527 
    ## Run 280 stress 0.06160425 
    ## Run 281 stress 0.05401009 
    ## Run 282 stress 0.05849517 
    ## Run 283 stress 0.05933305 
    ## Run 284 stress 0.06930436 
    ## Run 285 stress 0.07600337 
    ## Run 286 stress 0.04875116 
    ## ... Procrustes: rmse 0.04120931  max resid 0.1385727 
    ## Run 287 stress 0.0484589 
    ## ... Procrustes: rmse 5.7215e-05  max resid 0.0001021978 
    ## ... Similar to previous best
    ## Run 288 stress 0.04845891 
    ## ... Procrustes: rmse 6.529271e-05  max resid 0.0001146072 
    ## ... Similar to previous best
    ## Run 289 stress 0.05748887 
    ## Run 290 stress 0.04875105 
    ## ... Procrustes: rmse 0.04111907  max resid 0.1382938 
    ## Run 291 stress 0.05333363 
    ## Run 292 stress 0.05849515 
    ## Run 293 stress 0.06160425 
    ## Run 294 stress 0.04875107 
    ## ... Procrustes: rmse 0.04109266  max resid 0.1382121 
    ## Run 295 stress 0.05748878 
    ## Run 296 stress 0.05938482 
    ## Run 297 stress 0.0484589 
    ## ... Procrustes: rmse 4.666245e-05  max resid 7.101439e-05 
    ## ... Similar to previous best
    ## Run 298 stress 0.04845889 
    ## ... Procrustes: rmse 2.03373e-05  max resid 3.060196e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.04875106 
    ## ... Procrustes: rmse 0.04110917  max resid 0.1382642 
    ## Run 300 stress 0.05748882 
    ## Run 301 stress 0.06995442 
    ## Run 302 stress 0.05933131 
    ## Run 303 stress 0.04845892 
    ## ... Procrustes: rmse 7.883161e-05  max resid 0.0001384004 
    ## ... Similar to previous best
    ## Run 304 stress 0.04845891 
    ## ... Procrustes: rmse 4.781164e-05  max resid 7.175876e-05 
    ## ... Similar to previous best
    ## Run 305 stress 0.07131955 
    ## Run 306 stress 0.05849525 
    ## Run 307 stress 0.04845902 
    ## ... Procrustes: rmse 0.000163343  max resid 0.0002641646 
    ## ... Similar to previous best
    ## Run 308 stress 0.06096276 
    ## Run 309 stress 0.05933134 
    ## Run 310 stress 0.05748883 
    ## Run 311 stress 0.06160433 
    ## Run 312 stress 0.3620352 
    ## Run 313 stress 0.05992864 
    ## Run 314 stress 0.06897918 
    ## Run 315 stress 0.04845891 
    ## ... Procrustes: rmse 5.882514e-05  max resid 9.046557e-05 
    ## ... Similar to previous best
    ## Run 316 stress 0.06086825 
    ## Run 317 stress 0.04845895 
    ## ... Procrustes: rmse 0.0001074588  max resid 0.0001712334 
    ## ... Similar to previous best
    ## Run 318 stress 0.06203524 
    ## Run 319 stress 0.04845892 
    ## ... Procrustes: rmse 8.843441e-05  max resid 0.0001592638 
    ## ... Similar to previous best
    ## Run 320 stress 0.05849516 
    ## Run 321 stress 0.04875108 
    ## ... Procrustes: rmse 0.04117972  max resid 0.1384812 
    ## Run 322 stress 0.04875104 
    ## ... Procrustes: rmse 0.04115387  max resid 0.1384015 
    ## Run 323 stress 0.05933295 
    ## Run 324 stress 0.04845895 
    ## ... Procrustes: rmse 0.0001172122  max resid 0.0002126555 
    ## ... Similar to previous best
    ## Run 325 stress 0.04845889 
    ## ... Procrustes: rmse 3.551965e-05  max resid 6.269968e-05 
    ## ... Similar to previous best
    ## Run 326 stress 0.05933303 
    ## Run 327 stress 0.05401 
    ## Run 328 stress 0.05748878 
    ## Run 329 stress 0.05748881 
    ## Run 330 stress 0.05849521 
    ## Run 331 stress 0.05849518 
    ## Run 332 stress 0.07117135 
    ## Run 333 stress 0.04875108 
    ## ... Procrustes: rmse 0.04108583  max resid 0.1381912 
    ## Run 334 stress 0.05991818 
    ## Run 335 stress 0.06203521 
    ## Run 336 stress 0.0484589 
    ## ... Procrustes: rmse 3.417273e-05  max resid 5.232204e-05 
    ## ... Similar to previous best
    ## Run 337 stress 0.05849523 
    ## Run 338 stress 0.04875107 
    ## ... Procrustes: rmse 0.04115745  max resid 0.1384133 
    ## Run 339 stress 0.06041486 
    ## Run 340 stress 0.05118022 
    ## Run 341 stress 0.05933299 
    ## Run 342 stress 0.04875106 
    ## ... Procrustes: rmse 0.04109653  max resid 0.1382243 
    ## Run 343 stress 0.05849516 
    ## Run 344 stress 0.04875106 
    ## ... Procrustes: rmse 0.04110143  max resid 0.1382398 
    ## Run 345 stress 0.05118021 
    ## Run 346 stress 0.06614145 
    ## Run 347 stress 0.04875114 
    ## ... Procrustes: rmse 0.04105646  max resid 0.1381005 
    ## Run 348 stress 0.05849524 
    ## Run 349 stress 0.06041474 
    ## Run 350 stress 0.06203521 
    ## Run 351 stress 0.04845893 
    ## ... Procrustes: rmse 9.359307e-05  max resid 0.0001390415 
    ## ... Similar to previous best
    ## Run 352 stress 0.04845896 
    ## ... Procrustes: rmse 0.0001236133  max resid 0.0002133007 
    ## ... Similar to previous best
    ## Run 353 stress 0.0750464 
    ## Run 354 stress 0.04845891 
    ## ... Procrustes: rmse 7.522959e-05  max resid 0.000128162 
    ## ... Similar to previous best
    ## Run 355 stress 0.05333371 
    ## Run 356 stress 0.05456973 
    ## Run 357 stress 0.06969978 
    ## Run 358 stress 0.05400994 
    ## Run 359 stress 0.05333364 
    ## Run 360 stress 0.05849528 
    ## Run 361 stress 0.05991813 
    ## Run 362 stress 0.05348169 
    ## Run 363 stress 0.05933297 
    ## Run 364 stress 0.06939826 
    ## Run 365 stress 0.05667815 
    ## Run 366 stress 0.04845893 
    ## ... Procrustes: rmse 9.379135e-05  max resid 0.0001668368 
    ## ... Similar to previous best
    ## Run 367 stress 0.06036746 
    ## Run 368 stress 0.06036744 
    ## Run 369 stress 0.05938469 
    ## Run 370 stress 0.048459 
    ## ... Procrustes: rmse 3.573124e-05  max resid 7.622858e-05 
    ## ... Similar to previous best
    ## Run 371 stress 0.3620358 
    ## Run 372 stress 0.0593847 
    ## Run 373 stress 0.05333361 
    ## Run 374 stress 0.06174022 
    ## Run 375 stress 0.05992875 
    ## Run 376 stress 0.05991814 
    ## Run 377 stress 0.06236432 
    ## Run 378 stress 0.05933174 
    ## Run 379 stress 0.05991819 
    ## Run 380 stress 0.05938476 
    ## Run 381 stress 0.04875115 
    ## ... Procrustes: rmse 0.0410599  max resid 0.1381102 
    ## Run 382 stress 0.04845892 
    ## ... Procrustes: rmse 9.178333e-05  max resid 0.0001600898 
    ## ... Similar to previous best
    ## Run 383 stress 0.04875116 
    ## ... Procrustes: rmse 0.04105051  max resid 0.1380816 
    ## Run 384 stress 0.06036754 
    ## Run 385 stress 0.06174019 
    ## Run 386 stress 0.07271913 
    ## Run 387 stress 0.06096274 
    ## Run 388 stress 0.06617474 
    ## Run 389 stress 0.05333362 
    ## Run 390 stress 0.06614124 
    ## Run 391 stress 0.0584952 
    ## Run 392 stress 0.04875113 
    ## ... Procrustes: rmse 0.04120865  max resid 0.1385699 
    ## Run 393 stress 0.06897877 
    ## Run 394 stress 0.06174032 
    ## Run 395 stress 0.04845892 
    ## ... Procrustes: rmse 7.988457e-05  max resid 0.0001180883 
    ## ... Similar to previous best
    ## Run 396 stress 0.0484589 
    ## ... Procrustes: rmse 3.489454e-05  max resid 5.246108e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.04845891 
    ## ... Procrustes: rmse 7.032433e-05  max resid 0.0001081153 
    ## ... Similar to previous best
    ## Run 398 stress 0.05849517 
    ## Run 399 stress 0.04845895 
    ## ... Procrustes: rmse 0.0001228469  max resid 0.0002239167 
    ## ... Similar to previous best
    ## Run 400 stress 0.05992871 
    ## Run 401 stress 0.06174017 
    ## Run 402 stress 0.048459 
    ## ... Procrustes: rmse 0.0001536651  max resid 0.0002591455 
    ## ... Similar to previous best
    ## Run 403 stress 0.0484589 
    ## ... Procrustes: rmse 4.908369e-05  max resid 8.56329e-05 
    ## ... Similar to previous best
    ## Run 404 stress 0.0620353 
    ## Run 405 stress 0.05748889 
    ## Run 406 stress 0.06850994 
    ## Run 407 stress 0.05748879 
    ## Run 408 stress 0.06160425 
    ## Run 409 stress 0.05748879 
    ## Run 410 stress 0.06850963 
    ## Run 411 stress 0.06236427 
    ## Run 412 stress 0.05849517 
    ## Run 413 stress 0.04875107 
    ## ... Procrustes: rmse 0.04109327  max resid 0.1382137 
    ## Run 414 stress 0.04875111 
    ## ... Procrustes: rmse 0.04107381  max resid 0.1381545 
    ## Run 415 stress 0.04875114 
    ## ... Procrustes: rmse 0.04105605  max resid 0.1380989 
    ## Run 416 stress 0.0484589 
    ## ... Procrustes: rmse 5.021805e-05  max resid 9.081496e-05 
    ## ... Similar to previous best
    ## Run 417 stress 0.04875105 
    ## ... Procrustes: rmse 0.0411149  max resid 0.1382814 
    ## Run 418 stress 0.04845893 
    ## ... Procrustes: rmse 8.78343e-05  max resid 0.0001295442 
    ## ... Similar to previous best
    ## Run 419 stress 0.04845891 
    ## ... Procrustes: rmse 7.231169e-05  max resid 0.0001200169 
    ## ... Similar to previous best
    ## Run 420 stress 0.04875114 
    ## ... Procrustes: rmse 0.04105778  max resid 0.1381042 
    ## Run 421 stress 0.05991816 
    ## Run 422 stress 0.05118021 
    ## Run 423 stress 0.05933311 
    ## Run 424 stress 0.07023468 
    ## Run 425 stress 0.07600306 
    ## Run 426 stress 0.0586599 
    ## Run 427 stress 0.0683387 
    ## Run 428 stress 0.06739811 
    ## Run 429 stress 0.05933297 
    ## Run 430 stress 0.04845893 
    ## ... Procrustes: rmse 9.583956e-05  max resid 0.0001681379 
    ## ... Similar to previous best
    ## Run 431 stress 0.04875113 
    ## ... Procrustes: rmse 0.04106042  max resid 0.1381124 
    ## Run 432 stress 0.05992865 
    ## Run 433 stress 0.07061484 
    ## Run 434 stress 0.04845891 
    ## ... Procrustes: rmse 4.843681e-05  max resid 7.209197e-05 
    ## ... Similar to previous best
    ## Run 435 stress 0.07580057 
    ## Run 436 stress 0.0484589 
    ## ... Procrustes: rmse 4.644179e-05  max resid 8.239485e-05 
    ## ... Similar to previous best
    ## Run 437 stress 0.05865984 
    ## Run 438 stress 0.06614128 
    ## Run 439 stress 0.06878405 
    ## Run 440 stress 0.05992877 
    ## Run 441 stress 0.05333362 
    ## Run 442 stress 0.0484589 
    ## ... Procrustes: rmse 4.714168e-05  max resid 8.284557e-05 
    ## ... Similar to previous best
    ## Run 443 stress 0.05667814 
    ## Run 444 stress 0.06036744 
    ## Run 445 stress 0.05333364 
    ## Run 446 stress 0.05865982 
    ## Run 447 stress 0.04845891 
    ## ... Procrustes: rmse 6.115193e-05  max resid 9.12737e-05 
    ## ... Similar to previous best
    ## Run 448 stress 0.06096276 
    ## Run 449 stress 0.04845891 
    ## ... Procrustes: rmse 6.404043e-05  max resid 9.635011e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.04845889 
    ## ... Procrustes: rmse 3.050372e-05  max resid 5.604316e-05 
    ## ... Similar to previous best
    ## Run 451 stress 0.05333363 
    ## Run 452 stress 0.05992869 
    ## Run 453 stress 0.06041492 
    ## Run 454 stress 0.04875108 
    ## ... Procrustes: rmse 0.0410999  max resid 0.138234 
    ## Run 455 stress 0.06160419 
    ## Run 456 stress 0.05933297 
    ## Run 457 stress 0.06617254 
    ## Run 458 stress 0.06086814 
    ## Run 459 stress 0.04875106 
    ## ... Procrustes: rmse 0.04109722  max resid 0.1382264 
    ## Run 460 stress 0.0484589 
    ## ... Procrustes: rmse 3.752267e-05  max resid 6.168832e-05 
    ## ... Similar to previous best
    ## Run 461 stress 0.05333369 
    ## Run 462 stress 0.05865996 
    ## Run 463 stress 0.04875108 
    ## ... Procrustes: rmse 0.04118334  max resid 0.1384924 
    ## Run 464 stress 0.04875104 
    ## ... Procrustes: rmse 0.04111735  max resid 0.1382887 
    ## Run 465 stress 0.0484589 
    ## ... Procrustes: rmse 4.27969e-05  max resid 6.630996e-05 
    ## ... Similar to previous best
    ## Run 466 stress 0.06041479 
    ## Run 467 stress 0.04875111 
    ## ... Procrustes: rmse 0.04107034  max resid 0.1381432 
    ## Run 468 stress 0.05992866 
    ## Run 469 stress 0.06599252 
    ## Run 470 stress 0.06614151 
    ## Run 471 stress 0.05933312 
    ## Run 472 stress 0.04845895 
    ## ... Procrustes: rmse 0.0001147484  max resid 0.0001846497 
    ## ... Similar to previous best
    ## Run 473 stress 0.06317371 
    ## Run 474 stress 0.05748877 
    ## Run 475 stress 0.05748886 
    ## Run 476 stress 0.05400987 
    ## Run 477 stress 0.06911399 
    ## Run 478 stress 0.06614139 
    ## Run 479 stress 0.04845891 
    ## ... Procrustes: rmse 5.930849e-05  max resid 8.83691e-05 
    ## ... Similar to previous best
    ## Run 480 stress 0.06160424 
    ## Run 481 stress 0.05456973 
    ## Run 482 stress 0.04875105 
    ## ... Procrustes: rmse 0.04112757  max resid 0.1383209 
    ## Run 483 stress 0.0484589 
    ## ... Procrustes: rmse 4.469481e-05  max resid 7.804938e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.05748885 
    ## Run 485 stress 0.06236436 
    ## Run 486 stress 0.05118021 
    ## Run 487 stress 0.05348162 
    ## Run 488 stress 0.06036752 
    ## Run 489 stress 0.05865981 
    ## Run 490 stress 0.05933143 
    ## Run 491 stress 0.04845895 
    ## ... Procrustes: rmse 0.000108026  max resid 0.0001583117 
    ## ... Similar to previous best
    ## Run 492 stress 0.05748879 
    ## Run 493 stress 0.06036743 
    ## Run 494 stress 0.06971999 
    ## Run 495 stress 0.05991814 
    ## Run 496 stress 0.0484589 
    ## ... Procrustes: rmse 3.803685e-05  max resid 6.292225e-05 
    ## ... Similar to previous best
    ## Run 497 stress 0.05938473 
    ## Run 498 stress 0.04875114 
    ## ... Procrustes: rmse 0.04121229  max resid 0.1385809 
    ## Run 499 stress 0.05849516 
    ## Run 500 stress 0.06174021 
    ## *** Best solution repeated 82 times

``` r
# Ocean sites and mixed lakes
PD_beta_geo_OM_NMDS <- metaMDS(PD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09495903 
    ## Run 1 stress 0.138422 
    ## Run 2 stress 0.09495903 
    ## ... New best solution
    ## ... Procrustes: rmse 6.315405e-07  max resid 1.536231e-06 
    ## ... Similar to previous best
    ## Run 3 stress 0.09495903 
    ## ... Procrustes: rmse 1.724954e-06  max resid 3.438241e-06 
    ## ... Similar to previous best
    ## Run 4 stress 0.1264617 
    ## Run 5 stress 0.09495903 
    ## ... Procrustes: rmse 5.839111e-06  max resid 1.277622e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.09495903 
    ## ... Procrustes: rmse 7.582537e-06  max resid 1.780697e-05 
    ## ... Similar to previous best
    ## Run 7 stress 0.09495903 
    ## ... New best solution
    ## ... Procrustes: rmse 4.1516e-07  max resid 1.001376e-06 
    ## ... Similar to previous best
    ## Run 8 stress 0.09495903 
    ## ... Procrustes: rmse 1.688405e-06  max resid 4.654073e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.1264617 
    ## Run 10 stress 0.09495903 
    ## ... Procrustes: rmse 1.02855e-05  max resid 2.164176e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.1264617 
    ## Run 12 stress 0.1501836 
    ## Run 13 stress 0.09495903 
    ## ... Procrustes: rmse 7.442489e-06  max resid 1.777907e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.09495903 
    ## ... Procrustes: rmse 1.439747e-06  max resid 3.627402e-06 
    ## ... Similar to previous best
    ## Run 15 stress 0.09495903 
    ## ... Procrustes: rmse 8.141505e-06  max resid 1.92058e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.1264617 
    ## Run 17 stress 0.09495903 
    ## ... Procrustes: rmse 1.262242e-06  max resid 3.238334e-06 
    ## ... Similar to previous best
    ## Run 18 stress 0.09495903 
    ## ... Procrustes: rmse 1.412143e-06  max resid 3.36067e-06 
    ## ... Similar to previous best
    ## Run 19 stress 0.09495903 
    ## ... Procrustes: rmse 5.938486e-06  max resid 1.155898e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.09495903 
    ## ... Procrustes: rmse 5.013111e-06  max resid 1.073291e-05 
    ## ... Similar to previous best
    ## Run 21 stress 0.09495903 
    ## ... Procrustes: rmse 9.908091e-06  max resid 2.26959e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.09495903 
    ## ... Procrustes: rmse 9.189353e-06  max resid 2.259234e-05 
    ## ... Similar to previous best
    ## Run 23 stress 0.1384223 
    ## Run 24 stress 0.09495903 
    ## ... Procrustes: rmse 5.247891e-07  max resid 1.126014e-06 
    ## ... Similar to previous best
    ## Run 25 stress 0.1264617 
    ## Run 26 stress 0.09495903 
    ## ... New best solution
    ## ... Procrustes: rmse 3.486313e-07  max resid 7.23668e-07 
    ## ... Similar to previous best
    ## Run 27 stress 0.09495903 
    ## ... New best solution
    ## ... Procrustes: rmse 3.435729e-07  max resid 9.441217e-07 
    ## ... Similar to previous best
    ## Run 28 stress 0.09495903 
    ## ... Procrustes: rmse 8.522931e-06  max resid 2.001917e-05 
    ## ... Similar to previous best
    ## Run 29 stress 0.09495903 
    ## ... Procrustes: rmse 3.714714e-06  max resid 6.613256e-06 
    ## ... Similar to previous best
    ## Run 30 stress 0.09495903 
    ## ... Procrustes: rmse 2.10426e-06  max resid 5.170189e-06 
    ## ... Similar to previous best
    ## Run 31 stress 0.09495903 
    ## ... Procrustes: rmse 9.73358e-06  max resid 2.272351e-05 
    ## ... Similar to previous best
    ## Run 32 stress 0.09495903 
    ## ... Procrustes: rmse 5.600873e-07  max resid 1.341333e-06 
    ## ... Similar to previous best
    ## Run 33 stress 0.1264617 
    ## Run 34 stress 0.1521592 
    ## Run 35 stress 0.325977 
    ## Run 36 stress 0.1264617 
    ## Run 37 stress 0.09495903 
    ## ... Procrustes: rmse 8.39831e-07  max resid 1.944382e-06 
    ## ... Similar to previous best
    ## Run 38 stress 0.09495903 
    ## ... New best solution
    ## ... Procrustes: rmse 1.470396e-06  max resid 3.54182e-06 
    ## ... Similar to previous best
    ## Run 39 stress 0.09495903 
    ## ... Procrustes: rmse 1.590005e-06  max resid 3.654402e-06 
    ## ... Similar to previous best
    ## Run 40 stress 0.1521592 
    ## Run 41 stress 0.1521592 
    ## Run 42 stress 0.09495903 
    ## ... Procrustes: rmse 3.838644e-06  max resid 9.449694e-06 
    ## ... Similar to previous best
    ## Run 43 stress 0.09495903 
    ## ... Procrustes: rmse 1.598153e-06  max resid 3.23014e-06 
    ## ... Similar to previous best
    ## Run 44 stress 0.09495903 
    ## ... Procrustes: rmse 2.830421e-06  max resid 6.577634e-06 
    ## ... Similar to previous best
    ## Run 45 stress 0.09495903 
    ## ... Procrustes: rmse 3.5754e-06  max resid 8.915352e-06 
    ## ... Similar to previous best
    ## Run 46 stress 0.09495903 
    ## ... Procrustes: rmse 3.295149e-05  max resid 7.758065e-05 
    ## ... Similar to previous best
    ## Run 47 stress 0.09495903 
    ## ... Procrustes: rmse 2.209925e-06  max resid 4.372977e-06 
    ## ... Similar to previous best
    ## Run 48 stress 0.09495903 
    ## ... Procrustes: rmse 1.193052e-06  max resid 2.733676e-06 
    ## ... Similar to previous best
    ## Run 49 stress 0.1264617 
    ## Run 50 stress 0.09495903 
    ## ... Procrustes: rmse 3.113824e-06  max resid 7.990341e-06 
    ## ... Similar to previous best
    ## Run 51 stress 0.1264617 
    ## Run 52 stress 0.1264617 
    ## Run 53 stress 0.09495903 
    ## ... Procrustes: rmse 1.948812e-06  max resid 4.412306e-06 
    ## ... Similar to previous best
    ## Run 54 stress 0.09495903 
    ## ... Procrustes: rmse 2.856153e-06  max resid 7.009842e-06 
    ## ... Similar to previous best
    ## Run 55 stress 0.09495903 
    ## ... Procrustes: rmse 1.620162e-06  max resid 2.972022e-06 
    ## ... Similar to previous best
    ## Run 56 stress 0.1264617 
    ## Run 57 stress 0.09495903 
    ## ... Procrustes: rmse 1.797918e-06  max resid 3.055663e-06 
    ## ... Similar to previous best
    ## Run 58 stress 0.09495903 
    ## ... Procrustes: rmse 2.195701e-06  max resid 5.767633e-06 
    ## ... Similar to previous best
    ## Run 59 stress 0.09495903 
    ## ... Procrustes: rmse 6.188714e-06  max resid 1.566049e-05 
    ## ... Similar to previous best
    ## Run 60 stress 0.1264617 
    ## Run 61 stress 0.09495903 
    ## ... Procrustes: rmse 2.961308e-06  max resid 7.029599e-06 
    ## ... Similar to previous best
    ## Run 62 stress 0.09495903 
    ## ... Procrustes: rmse 4.720342e-06  max resid 1.154289e-05 
    ## ... Similar to previous best
    ## Run 63 stress 0.1264617 
    ## Run 64 stress 0.1264617 
    ## Run 65 stress 0.09495903 
    ## ... New best solution
    ## ... Procrustes: rmse 5.790253e-07  max resid 1.229114e-06 
    ## ... Similar to previous best
    ## Run 66 stress 0.1264617 
    ## Run 67 stress 0.09495903 
    ## ... Procrustes: rmse 9.262077e-07  max resid 2.519904e-06 
    ## ... Similar to previous best
    ## Run 68 stress 0.09495903 
    ## ... Procrustes: rmse 7.082727e-06  max resid 1.702672e-05 
    ## ... Similar to previous best
    ## Run 69 stress 0.1501835 
    ## Run 70 stress 0.09495903 
    ## ... Procrustes: rmse 1.440047e-05  max resid 3.42486e-05 
    ## ... Similar to previous best
    ## Run 71 stress 0.1264617 
    ## Run 72 stress 0.09495903 
    ## ... Procrustes: rmse 5.755574e-06  max resid 1.468102e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.09495903 
    ## ... Procrustes: rmse 1.046841e-06  max resid 2.489485e-06 
    ## ... Similar to previous best
    ## Run 74 stress 0.09495903 
    ## ... Procrustes: rmse 7.422445e-07  max resid 2.018767e-06 
    ## ... Similar to previous best
    ## Run 75 stress 0.1384222 
    ## Run 76 stress 0.09495903 
    ## ... Procrustes: rmse 2.032638e-06  max resid 4.857109e-06 
    ## ... Similar to previous best
    ## Run 77 stress 0.09495903 
    ## ... Procrustes: rmse 1.777071e-05  max resid 4.227873e-05 
    ## ... Similar to previous best
    ## Run 78 stress 0.09495903 
    ## ... Procrustes: rmse 1.032402e-05  max resid 2.464393e-05 
    ## ... Similar to previous best
    ## Run 79 stress 0.09495903 
    ## ... Procrustes: rmse 2.643932e-06  max resid 6.417605e-06 
    ## ... Similar to previous best
    ## Run 80 stress 0.09495903 
    ## ... Procrustes: rmse 2.124195e-06  max resid 3.819674e-06 
    ## ... Similar to previous best
    ## Run 81 stress 0.09495903 
    ## ... New best solution
    ## ... Procrustes: rmse 6.356918e-07  max resid 1.44068e-06 
    ## ... Similar to previous best
    ## Run 82 stress 0.09495903 
    ## ... Procrustes: rmse 3.535918e-06  max resid 8.253246e-06 
    ## ... Similar to previous best
    ## Run 83 stress 0.09495903 
    ## ... Procrustes: rmse 2.533642e-06  max resid 5.374645e-06 
    ## ... Similar to previous best
    ## Run 84 stress 0.09495903 
    ## ... Procrustes: rmse 4.221298e-06  max resid 1.020206e-05 
    ## ... Similar to previous best
    ## Run 85 stress 0.1384219 
    ## Run 86 stress 0.09495903 
    ## ... Procrustes: rmse 2.169012e-06  max resid 3.535863e-06 
    ## ... Similar to previous best
    ## Run 87 stress 0.09495903 
    ## ... Procrustes: rmse 1.64515e-06  max resid 3.741461e-06 
    ## ... Similar to previous best
    ## Run 88 stress 0.09495903 
    ## ... Procrustes: rmse 2.914229e-06  max resid 7.943698e-06 
    ## ... Similar to previous best
    ## Run 89 stress 0.09495903 
    ## ... Procrustes: rmse 2.978977e-06  max resid 7.030517e-06 
    ## ... Similar to previous best
    ## Run 90 stress 0.1264617 
    ## Run 91 stress 0.09495903 
    ## ... Procrustes: rmse 1.454612e-06  max resid 4.247188e-06 
    ## ... Similar to previous best
    ## Run 92 stress 0.09495903 
    ## ... Procrustes: rmse 5.050268e-06  max resid 1.342478e-05 
    ## ... Similar to previous best
    ## Run 93 stress 0.1264617 
    ## Run 94 stress 0.09495903 
    ## ... Procrustes: rmse 1.47573e-06  max resid 4.010274e-06 
    ## ... Similar to previous best
    ## Run 95 stress 0.138422 
    ## Run 96 stress 0.1264617 
    ## Run 97 stress 0.09495903 
    ## ... Procrustes: rmse 6.692407e-06  max resid 1.557671e-05 
    ## ... Similar to previous best
    ## Run 98 stress 0.09495903 
    ## ... Procrustes: rmse 2.60038e-06  max resid 6.130034e-06 
    ## ... Similar to previous best
    ## Run 99 stress 0.1264617 
    ## Run 100 stress 0.09495903 
    ## ... Procrustes: rmse 1.204513e-05  max resid 2.894419e-05 
    ## ... Similar to previous best
    ## Run 101 stress 0.1521592 
    ## Run 102 stress 0.09495903 
    ## ... New best solution
    ## ... Procrustes: rmse 6.006424e-07  max resid 1.075531e-06 
    ## ... Similar to previous best
    ## Run 103 stress 0.09495903 
    ## ... Procrustes: rmse 7.056391e-07  max resid 1.809303e-06 
    ## ... Similar to previous best
    ## Run 104 stress 0.09495903 
    ## ... Procrustes: rmse 2.537132e-06  max resid 6.109726e-06 
    ## ... Similar to previous best
    ## Run 105 stress 0.1264617 
    ## Run 106 stress 0.09495903 
    ## ... Procrustes: rmse 1.911522e-06  max resid 4.490803e-06 
    ## ... Similar to previous best
    ## Run 107 stress 0.09495903 
    ## ... Procrustes: rmse 1.709042e-06  max resid 3.787926e-06 
    ## ... Similar to previous best
    ## Run 108 stress 0.09495903 
    ## ... Procrustes: rmse 3.110125e-06  max resid 7.603055e-06 
    ## ... Similar to previous best
    ## Run 109 stress 0.09495903 
    ## ... Procrustes: rmse 6.235247e-07  max resid 9.36071e-07 
    ## ... Similar to previous best
    ## Run 110 stress 0.09495903 
    ## ... Procrustes: rmse 2.573315e-06  max resid 6.18897e-06 
    ## ... Similar to previous best
    ## Run 111 stress 0.09495903 
    ## ... Procrustes: rmse 2.072036e-06  max resid 5.259195e-06 
    ## ... Similar to previous best
    ## Run 112 stress 0.1264617 
    ## Run 113 stress 0.1264617 
    ## Run 114 stress 0.09495903 
    ## ... Procrustes: rmse 2.0822e-06  max resid 5.426832e-06 
    ## ... Similar to previous best
    ## Run 115 stress 0.1264617 
    ## Run 116 stress 0.1264617 
    ## Run 117 stress 0.09495903 
    ## ... Procrustes: rmse 5.957065e-06  max resid 1.056778e-05 
    ## ... Similar to previous best
    ## Run 118 stress 0.09495903 
    ## ... Procrustes: rmse 1.842695e-06  max resid 4.326048e-06 
    ## ... Similar to previous best
    ## Run 119 stress 0.09495903 
    ## ... Procrustes: rmse 1.692315e-06  max resid 3.976569e-06 
    ## ... Similar to previous best
    ## Run 120 stress 0.1264617 
    ## Run 121 stress 0.09495903 
    ## ... Procrustes: rmse 1.420083e-05  max resid 3.144644e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.09495903 
    ## ... Procrustes: rmse 3.056152e-06  max resid 7.311e-06 
    ## ... Similar to previous best
    ## Run 123 stress 0.09495903 
    ## ... Procrustes: rmse 3.017493e-06  max resid 7.369071e-06 
    ## ... Similar to previous best
    ## Run 124 stress 0.304369 
    ## Run 125 stress 0.09495903 
    ## ... Procrustes: rmse 2.163595e-06  max resid 4.183149e-06 
    ## ... Similar to previous best
    ## Run 126 stress 0.1264617 
    ## Run 127 stress 0.09495903 
    ## ... Procrustes: rmse 1.744014e-06  max resid 4.950964e-06 
    ## ... Similar to previous best
    ## Run 128 stress 0.09495905 
    ## ... Procrustes: rmse 3.132191e-05  max resid 7.077899e-05 
    ## ... Similar to previous best
    ## Run 129 stress 0.09495903 
    ## ... Procrustes: rmse 2.884706e-06  max resid 7.014147e-06 
    ## ... Similar to previous best
    ## Run 130 stress 0.09495903 
    ## ... Procrustes: rmse 1.218115e-06  max resid 2.359994e-06 
    ## ... Similar to previous best
    ## Run 131 stress 0.09495903 
    ## ... Procrustes: rmse 1.11453e-05  max resid 2.632861e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.09495903 
    ## ... Procrustes: rmse 2.363846e-06  max resid 5.567868e-06 
    ## ... Similar to previous best
    ## Run 133 stress 0.1501835 
    ## Run 134 stress 0.09495903 
    ## ... Procrustes: rmse 5.029379e-06  max resid 1.140039e-05 
    ## ... Similar to previous best
    ## Run 135 stress 0.09495903 
    ## ... Procrustes: rmse 1.614089e-06  max resid 3.908726e-06 
    ## ... Similar to previous best
    ## Run 136 stress 0.1264617 
    ## Run 137 stress 0.09495903 
    ## ... Procrustes: rmse 9.628454e-06  max resid 2.280838e-05 
    ## ... Similar to previous best
    ## Run 138 stress 0.138422 
    ## Run 139 stress 0.09495903 
    ## ... Procrustes: rmse 3.925946e-06  max resid 9.100615e-06 
    ## ... Similar to previous best
    ## Run 140 stress 0.09495903 
    ## ... Procrustes: rmse 2.776366e-06  max resid 6.793882e-06 
    ## ... Similar to previous best
    ## Run 141 stress 0.1264617 
    ## Run 142 stress 0.138422 
    ## Run 143 stress 0.09495903 
    ## ... Procrustes: rmse 1.550872e-05  max resid 3.674771e-05 
    ## ... Similar to previous best
    ## Run 144 stress 0.2972401 
    ## Run 145 stress 0.09495903 
    ## ... Procrustes: rmse 1.66246e-06  max resid 3.928758e-06 
    ## ... Similar to previous best
    ## Run 146 stress 0.1264617 
    ## Run 147 stress 0.09495903 
    ## ... Procrustes: rmse 1.766947e-06  max resid 4.246417e-06 
    ## ... Similar to previous best
    ## Run 148 stress 0.09495903 
    ## ... New best solution
    ## ... Procrustes: rmse 3.368961e-07  max resid 7.491669e-07 
    ## ... Similar to previous best
    ## Run 149 stress 0.1264617 
    ## Run 150 stress 0.1264617 
    ## Run 151 stress 0.09495903 
    ## ... Procrustes: rmse 1.161519e-06  max resid 2.856094e-06 
    ## ... Similar to previous best
    ## Run 152 stress 0.1264617 
    ## Run 153 stress 0.1264617 
    ## Run 154 stress 0.09495903 
    ## ... Procrustes: rmse 2.093959e-06  max resid 5.48392e-06 
    ## ... Similar to previous best
    ## Run 155 stress 0.09495903 
    ## ... Procrustes: rmse 3.876636e-06  max resid 9.552002e-06 
    ## ... Similar to previous best
    ## Run 156 stress 0.1264617 
    ## Run 157 stress 0.1384221 
    ## Run 158 stress 0.09495903 
    ## ... Procrustes: rmse 7.909263e-06  max resid 1.792809e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.1264617 
    ## Run 160 stress 0.09495903 
    ## ... Procrustes: rmse 3.965841e-06  max resid 9.653309e-06 
    ## ... Similar to previous best
    ## Run 161 stress 0.09495903 
    ## ... Procrustes: rmse 8.656102e-07  max resid 1.932954e-06 
    ## ... Similar to previous best
    ## Run 162 stress 0.09495903 
    ## ... Procrustes: rmse 1.464545e-05  max resid 3.526161e-05 
    ## ... Similar to previous best
    ## Run 163 stress 0.09495903 
    ## ... Procrustes: rmse 2.567983e-06  max resid 6.110256e-06 
    ## ... Similar to previous best
    ## Run 164 stress 0.1264617 
    ## Run 165 stress 0.09495903 
    ## ... Procrustes: rmse 8.395971e-07  max resid 2.106742e-06 
    ## ... Similar to previous best
    ## Run 166 stress 0.09495903 
    ## ... Procrustes: rmse 2.192076e-06  max resid 5.320115e-06 
    ## ... Similar to previous best
    ## Run 167 stress 0.09495903 
    ## ... Procrustes: rmse 1.029357e-05  max resid 2.360744e-05 
    ## ... Similar to previous best
    ## Run 168 stress 0.1264617 
    ## Run 169 stress 0.1521592 
    ## Run 170 stress 0.1264617 
    ## Run 171 stress 0.1264617 
    ## Run 172 stress 0.09495903 
    ## ... Procrustes: rmse 1.385487e-06  max resid 3.489383e-06 
    ## ... Similar to previous best
    ## Run 173 stress 0.09495903 
    ## ... Procrustes: rmse 2.821155e-06  max resid 6.839591e-06 
    ## ... Similar to previous best
    ## Run 174 stress 0.09495903 
    ## ... Procrustes: rmse 2.027366e-05  max resid 4.802419e-05 
    ## ... Similar to previous best
    ## Run 175 stress 0.3158094 
    ## Run 176 stress 0.1264617 
    ## Run 177 stress 0.138422 
    ## Run 178 stress 0.09495903 
    ## ... Procrustes: rmse 2.374764e-06  max resid 5.799341e-06 
    ## ... Similar to previous best
    ## Run 179 stress 0.3105733 
    ## Run 180 stress 0.09495903 
    ## ... Procrustes: rmse 1.806986e-06  max resid 2.676482e-06 
    ## ... Similar to previous best
    ## Run 181 stress 0.09495903 
    ## ... Procrustes: rmse 4.609736e-07  max resid 9.61038e-07 
    ## ... Similar to previous best
    ## Run 182 stress 0.09495903 
    ## ... Procrustes: rmse 1.559915e-06  max resid 3.241213e-06 
    ## ... Similar to previous best
    ## Run 183 stress 0.1264617 
    ## Run 184 stress 0.09495903 
    ## ... Procrustes: rmse 3.266773e-07  max resid 8.052779e-07 
    ## ... Similar to previous best
    ## Run 185 stress 0.09495903 
    ## ... Procrustes: rmse 2.906655e-06  max resid 7.082315e-06 
    ## ... Similar to previous best
    ## Run 186 stress 0.09495903 
    ## ... Procrustes: rmse 4.284649e-06  max resid 1.026395e-05 
    ## ... Similar to previous best
    ## Run 187 stress 0.1264617 
    ## Run 188 stress 0.09495903 
    ## ... Procrustes: rmse 3.368961e-06  max resid 8.122005e-06 
    ## ... Similar to previous best
    ## Run 189 stress 0.09495903 
    ## ... Procrustes: rmse 3.230027e-06  max resid 7.81436e-06 
    ## ... Similar to previous best
    ## Run 190 stress 0.3332219 
    ## Run 191 stress 0.09495903 
    ## ... Procrustes: rmse 1.35873e-06  max resid 3.196207e-06 
    ## ... Similar to previous best
    ## Run 192 stress 0.3240618 
    ## Run 193 stress 0.1264617 
    ## Run 194 stress 0.09495903 
    ## ... Procrustes: rmse 2.656952e-06  max resid 6.282263e-06 
    ## ... Similar to previous best
    ## Run 195 stress 0.1384219 
    ## Run 196 stress 0.09495903 
    ## ... Procrustes: rmse 2.518439e-06  max resid 5.983021e-06 
    ## ... Similar to previous best
    ## Run 197 stress 0.09495903 
    ## ... Procrustes: rmse 2.84844e-06  max resid 7.007131e-06 
    ## ... Similar to previous best
    ## Run 198 stress 0.09495903 
    ## ... Procrustes: rmse 8.686166e-06  max resid 1.644385e-05 
    ## ... Similar to previous best
    ## Run 199 stress 0.09495903 
    ## ... Procrustes: rmse 5.814614e-06  max resid 1.342154e-05 
    ## ... Similar to previous best
    ## Run 200 stress 0.1384219 
    ## Run 201 stress 0.09495903 
    ## ... Procrustes: rmse 6.739949e-07  max resid 1.300162e-06 
    ## ... Similar to previous best
    ## Run 202 stress 0.09495903 
    ## ... Procrustes: rmse 3.640935e-06  max resid 8.919514e-06 
    ## ... Similar to previous best
    ## Run 203 stress 0.1264617 
    ## Run 204 stress 0.1501836 
    ## Run 205 stress 0.1501838 
    ## Run 206 stress 0.1384221 
    ## Run 207 stress 0.09495903 
    ## ... Procrustes: rmse 8.76532e-06  max resid 2.115943e-05 
    ## ... Similar to previous best
    ## Run 208 stress 0.09495903 
    ## ... Procrustes: rmse 1.794395e-06  max resid 3.602986e-06 
    ## ... Similar to previous best
    ## Run 209 stress 0.1264617 
    ## Run 210 stress 0.1264617 
    ## Run 211 stress 0.09495903 
    ## ... Procrustes: rmse 1.318194e-06  max resid 2.963004e-06 
    ## ... Similar to previous best
    ## Run 212 stress 0.1264617 
    ## Run 213 stress 0.1521592 
    ## Run 214 stress 0.09495903 
    ## ... Procrustes: rmse 5.252121e-07  max resid 8.112329e-07 
    ## ... Similar to previous best
    ## Run 215 stress 0.1384222 
    ## Run 216 stress 0.09495903 
    ## ... Procrustes: rmse 1.091867e-06  max resid 2.619586e-06 
    ## ... Similar to previous best
    ## Run 217 stress 0.1264617 
    ## Run 218 stress 0.09495903 
    ## ... Procrustes: rmse 6.903572e-07  max resid 1.500693e-06 
    ## ... Similar to previous best
    ## Run 219 stress 0.09495903 
    ## ... Procrustes: rmse 4.200621e-06  max resid 9.934764e-06 
    ## ... Similar to previous best
    ## Run 220 stress 0.1264617 
    ## Run 221 stress 0.09495903 
    ## ... Procrustes: rmse 1.967769e-07  max resid 4.79916e-07 
    ## ... Similar to previous best
    ## Run 222 stress 0.09495903 
    ## ... Procrustes: rmse 1.369342e-06  max resid 3.164544e-06 
    ## ... Similar to previous best
    ## Run 223 stress 0.1264617 
    ## Run 224 stress 0.09495903 
    ## ... Procrustes: rmse 2.222414e-06  max resid 5.425268e-06 
    ## ... Similar to previous best
    ## Run 225 stress 0.09495903 
    ## ... Procrustes: rmse 5.143516e-07  max resid 8.761965e-07 
    ## ... Similar to previous best
    ## Run 226 stress 0.09495903 
    ## ... Procrustes: rmse 8.264747e-07  max resid 1.536446e-06 
    ## ... Similar to previous best
    ## Run 227 stress 0.09495903 
    ## ... Procrustes: rmse 1.549054e-06  max resid 3.877506e-06 
    ## ... Similar to previous best
    ## Run 228 stress 0.09495903 
    ## ... Procrustes: rmse 3.488603e-06  max resid 7.732422e-06 
    ## ... Similar to previous best
    ## Run 229 stress 0.1384223 
    ## Run 230 stress 0.09495903 
    ## ... Procrustes: rmse 9.598265e-06  max resid 2.209333e-05 
    ## ... Similar to previous best
    ## Run 231 stress 0.1264617 
    ## Run 232 stress 0.1521592 
    ## Run 233 stress 0.1264617 
    ## Run 234 stress 0.09495903 
    ## ... Procrustes: rmse 2.521855e-06  max resid 5.962787e-06 
    ## ... Similar to previous best
    ## Run 235 stress 0.1384223 
    ## Run 236 stress 0.09495903 
    ## ... Procrustes: rmse 3.777503e-06  max resid 9.055023e-06 
    ## ... Similar to previous best
    ## Run 237 stress 0.09495903 
    ## ... Procrustes: rmse 2.117147e-06  max resid 5.908055e-06 
    ## ... Similar to previous best
    ## Run 238 stress 0.1264617 
    ## Run 239 stress 0.1384224 
    ## Run 240 stress 0.1264617 
    ## Run 241 stress 0.09495903 
    ## ... Procrustes: rmse 5.095996e-06  max resid 1.298432e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.1264617 
    ## Run 243 stress 0.09495903 
    ## ... Procrustes: rmse 2.09726e-06  max resid 5.083324e-06 
    ## ... Similar to previous best
    ## Run 244 stress 0.09495903 
    ## ... Procrustes: rmse 4.633881e-06  max resid 1.164545e-05 
    ## ... Similar to previous best
    ## Run 245 stress 0.09495903 
    ## ... Procrustes: rmse 3.67703e-06  max resid 8.788024e-06 
    ## ... Similar to previous best
    ## Run 246 stress 0.1264617 
    ## Run 247 stress 0.09495903 
    ## ... Procrustes: rmse 4.859315e-06  max resid 1.13424e-05 
    ## ... Similar to previous best
    ## Run 248 stress 0.09495903 
    ## ... Procrustes: rmse 1.027795e-05  max resid 2.475146e-05 
    ## ... Similar to previous best
    ## Run 249 stress 0.09495903 
    ## ... Procrustes: rmse 1.427458e-06  max resid 2.843213e-06 
    ## ... Similar to previous best
    ## Run 250 stress 0.1264617 
    ## Run 251 stress 0.1264617 
    ## Run 252 stress 0.1501835 
    ## Run 253 stress 0.1384221 
    ## Run 254 stress 0.1264617 
    ## Run 255 stress 0.1264617 
    ## Run 256 stress 0.09495903 
    ## ... Procrustes: rmse 3.241616e-06  max resid 7.944499e-06 
    ## ... Similar to previous best
    ## Run 257 stress 0.09495903 
    ## ... Procrustes: rmse 2.07056e-06  max resid 3.731617e-06 
    ## ... Similar to previous best
    ## Run 258 stress 0.09495903 
    ## ... Procrustes: rmse 5.3211e-07  max resid 1.411683e-06 
    ## ... Similar to previous best
    ## Run 259 stress 0.09495903 
    ## ... Procrustes: rmse 1.433958e-06  max resid 3.431174e-06 
    ## ... Similar to previous best
    ## Run 260 stress 0.09495903 
    ## ... Procrustes: rmse 1.266808e-06  max resid 3.022179e-06 
    ## ... Similar to previous best
    ## Run 261 stress 0.1264617 
    ## Run 262 stress 0.09495903 
    ## ... Procrustes: rmse 1.039364e-05  max resid 2.954199e-05 
    ## ... Similar to previous best
    ## Run 263 stress 0.09495903 
    ## ... Procrustes: rmse 2.456864e-06  max resid 6.1224e-06 
    ## ... Similar to previous best
    ## Run 264 stress 0.09495903 
    ## ... Procrustes: rmse 5.827302e-07  max resid 9.799013e-07 
    ## ... Similar to previous best
    ## Run 265 stress 0.09495903 
    ## ... Procrustes: rmse 4.523712e-07  max resid 1.064438e-06 
    ## ... Similar to previous best
    ## Run 266 stress 0.09495903 
    ## ... Procrustes: rmse 1.728288e-06  max resid 3.979131e-06 
    ## ... Similar to previous best
    ## Run 267 stress 0.09495903 
    ## ... Procrustes: rmse 1.002814e-05  max resid 2.362209e-05 
    ## ... Similar to previous best
    ## Run 268 stress 0.1264617 
    ## Run 269 stress 0.09495903 
    ## ... Procrustes: rmse 5.058718e-06  max resid 1.185553e-05 
    ## ... Similar to previous best
    ## Run 270 stress 0.09495903 
    ## ... Procrustes: rmse 1.056092e-06  max resid 2.443785e-06 
    ## ... Similar to previous best
    ## Run 271 stress 0.09495903 
    ## ... Procrustes: rmse 4.256047e-06  max resid 1.038778e-05 
    ## ... Similar to previous best
    ## Run 272 stress 0.09495903 
    ## ... Procrustes: rmse 2.848272e-06  max resid 7.711602e-06 
    ## ... Similar to previous best
    ## Run 273 stress 0.1264617 
    ## Run 274 stress 0.1264617 
    ## Run 275 stress 0.09495903 
    ## ... Procrustes: rmse 1.434363e-06  max resid 3.50544e-06 
    ## ... Similar to previous best
    ## Run 276 stress 0.09495903 
    ## ... Procrustes: rmse 1.696851e-05  max resid 4.031826e-05 
    ## ... Similar to previous best
    ## Run 277 stress 0.09495903 
    ## ... Procrustes: rmse 9.441795e-07  max resid 2.063272e-06 
    ## ... Similar to previous best
    ## Run 278 stress 0.09495903 
    ## ... Procrustes: rmse 3.062721e-06  max resid 7.945568e-06 
    ## ... Similar to previous best
    ## Run 279 stress 0.09495903 
    ## ... Procrustes: rmse 4.105357e-06  max resid 1.011799e-05 
    ## ... Similar to previous best
    ## Run 280 stress 0.1501838 
    ## Run 281 stress 0.09495903 
    ## ... Procrustes: rmse 1.816851e-06  max resid 3.435851e-06 
    ## ... Similar to previous best
    ## Run 282 stress 0.1501836 
    ## Run 283 stress 0.09495903 
    ## ... Procrustes: rmse 3.301101e-06  max resid 7.649185e-06 
    ## ... Similar to previous best
    ## Run 284 stress 0.09495903 
    ## ... Procrustes: rmse 5.972449e-06  max resid 1.581467e-05 
    ## ... Similar to previous best
    ## Run 285 stress 0.09495903 
    ## ... Procrustes: rmse 2.228252e-06  max resid 5.43154e-06 
    ## ... Similar to previous best
    ## Run 286 stress 0.1264617 
    ## Run 287 stress 0.09495903 
    ## ... Procrustes: rmse 2.492434e-06  max resid 6.198928e-06 
    ## ... Similar to previous best
    ## Run 288 stress 0.09495903 
    ## ... Procrustes: rmse 1.070325e-06  max resid 2.584615e-06 
    ## ... Similar to previous best
    ## Run 289 stress 0.1384223 
    ## Run 290 stress 0.09495903 
    ## ... Procrustes: rmse 1.540745e-06  max resid 3.746263e-06 
    ## ... Similar to previous best
    ## Run 291 stress 0.1264617 
    ## Run 292 stress 0.1264617 
    ## Run 293 stress 0.09495903 
    ## ... Procrustes: rmse 1.021655e-05  max resid 1.708182e-05 
    ## ... Similar to previous best
    ## Run 294 stress 0.09495903 
    ## ... Procrustes: rmse 3.499797e-06  max resid 8.711954e-06 
    ## ... Similar to previous best
    ## Run 295 stress 0.1264617 
    ## Run 296 stress 0.09495903 
    ## ... Procrustes: rmse 6.029494e-06  max resid 1.373548e-05 
    ## ... Similar to previous best
    ## Run 297 stress 0.3531928 
    ## Run 298 stress 0.1264617 
    ## Run 299 stress 0.09495903 
    ## ... Procrustes: rmse 2.687168e-06  max resid 5.220911e-06 
    ## ... Similar to previous best
    ## Run 300 stress 0.09495903 
    ## ... Procrustes: rmse 3.025023e-06  max resid 7.579022e-06 
    ## ... Similar to previous best
    ## Run 301 stress 0.09495903 
    ## ... Procrustes: rmse 1.592489e-06  max resid 3.533876e-06 
    ## ... Similar to previous best
    ## Run 302 stress 0.1264617 
    ## Run 303 stress 0.09495903 
    ## ... Procrustes: rmse 1.968798e-06  max resid 4.650267e-06 
    ## ... Similar to previous best
    ## Run 304 stress 0.1264617 
    ## Run 305 stress 0.09495903 
    ## ... Procrustes: rmse 5.274998e-07  max resid 1.23934e-06 
    ## ... Similar to previous best
    ## Run 306 stress 0.09495903 
    ## ... Procrustes: rmse 1.467262e-06  max resid 3.44475e-06 
    ## ... Similar to previous best
    ## Run 307 stress 0.09495903 
    ## ... Procrustes: rmse 9.883621e-07  max resid 2.25588e-06 
    ## ... Similar to previous best
    ## Run 308 stress 0.09495903 
    ## ... Procrustes: rmse 1.295406e-06  max resid 3.700754e-06 
    ## ... Similar to previous best
    ## Run 309 stress 0.09495903 
    ## ... Procrustes: rmse 4.656638e-06  max resid 1.112239e-05 
    ## ... Similar to previous best
    ## Run 310 stress 0.09495903 
    ## ... Procrustes: rmse 1.397046e-06  max resid 3.18161e-06 
    ## ... Similar to previous best
    ## Run 311 stress 0.09495903 
    ## ... Procrustes: rmse 2.063698e-06  max resid 4.22687e-06 
    ## ... Similar to previous best
    ## Run 312 stress 0.1264617 
    ## Run 313 stress 0.09495903 
    ## ... Procrustes: rmse 1.567352e-05  max resid 3.404936e-05 
    ## ... Similar to previous best
    ## Run 314 stress 0.09495903 
    ## ... Procrustes: rmse 7.157599e-06  max resid 1.729186e-05 
    ## ... Similar to previous best
    ## Run 315 stress 0.09495903 
    ## ... Procrustes: rmse 3.624028e-06  max resid 8.662936e-06 
    ## ... Similar to previous best
    ## Run 316 stress 0.09495903 
    ## ... Procrustes: rmse 3.370499e-06  max resid 8.353004e-06 
    ## ... Similar to previous best
    ## Run 317 stress 0.09495903 
    ## ... Procrustes: rmse 3.404403e-06  max resid 8.739953e-06 
    ## ... Similar to previous best
    ## Run 318 stress 0.1501837 
    ## Run 319 stress 0.09495903 
    ## ... Procrustes: rmse 5.016289e-06  max resid 1.160901e-05 
    ## ... Similar to previous best
    ## Run 320 stress 0.09495903 
    ## ... Procrustes: rmse 2.368737e-06  max resid 5.840967e-06 
    ## ... Similar to previous best
    ## Run 321 stress 0.09495903 
    ## ... Procrustes: rmse 1.9906e-06  max resid 4.127425e-06 
    ## ... Similar to previous best
    ## Run 322 stress 0.09495903 
    ## ... Procrustes: rmse 7.842599e-07  max resid 1.700619e-06 
    ## ... Similar to previous best
    ## Run 323 stress 0.1521592 
    ## Run 324 stress 0.09495903 
    ## ... Procrustes: rmse 5.040199e-07  max resid 1.156227e-06 
    ## ... Similar to previous best
    ## Run 325 stress 0.1264617 
    ## Run 326 stress 0.09495903 
    ## ... Procrustes: rmse 1.721366e-06  max resid 4.154307e-06 
    ## ... Similar to previous best
    ## Run 327 stress 0.09495903 
    ## ... Procrustes: rmse 8.665393e-06  max resid 1.570877e-05 
    ## ... Similar to previous best
    ## Run 328 stress 0.1264617 
    ## Run 329 stress 0.1384221 
    ## Run 330 stress 0.09495903 
    ## ... Procrustes: rmse 8.572807e-06  max resid 1.982686e-05 
    ## ... Similar to previous best
    ## Run 331 stress 0.09495903 
    ## ... Procrustes: rmse 2.895413e-06  max resid 7.083367e-06 
    ## ... Similar to previous best
    ## Run 332 stress 0.09495903 
    ## ... Procrustes: rmse 1.15405e-05  max resid 2.07655e-05 
    ## ... Similar to previous best
    ## Run 333 stress 0.09495903 
    ## ... Procrustes: rmse 2.697623e-06  max resid 6.436411e-06 
    ## ... Similar to previous best
    ## Run 334 stress 0.1264617 
    ## Run 335 stress 0.09495903 
    ## ... Procrustes: rmse 1.572641e-06  max resid 3.627684e-06 
    ## ... Similar to previous best
    ## Run 336 stress 0.09495903 
    ## ... Procrustes: rmse 1.835506e-06  max resid 4.179779e-06 
    ## ... Similar to previous best
    ## Run 337 stress 0.09495903 
    ## ... Procrustes: rmse 3.370337e-06  max resid 8.10732e-06 
    ## ... Similar to previous best
    ## Run 338 stress 0.09495903 
    ## ... Procrustes: rmse 2.01127e-06  max resid 4.87509e-06 
    ## ... Similar to previous best
    ## Run 339 stress 0.1264617 
    ## Run 340 stress 0.1264617 
    ## Run 341 stress 0.09495903 
    ## ... Procrustes: rmse 1.444431e-06  max resid 3.448432e-06 
    ## ... Similar to previous best
    ## Run 342 stress 0.09495903 
    ## ... Procrustes: rmse 5.660782e-07  max resid 9.406857e-07 
    ## ... Similar to previous best
    ## Run 343 stress 0.09495903 
    ## ... Procrustes: rmse 9.069967e-06  max resid 2.063283e-05 
    ## ... Similar to previous best
    ## Run 344 stress 0.1264617 
    ## Run 345 stress 0.09495903 
    ## ... Procrustes: rmse 4.528358e-06  max resid 7.460908e-06 
    ## ... Similar to previous best
    ## Run 346 stress 0.09495903 
    ## ... Procrustes: rmse 2.027316e-06  max resid 4.6392e-06 
    ## ... Similar to previous best
    ## Run 347 stress 0.1264617 
    ## Run 348 stress 0.1384219 
    ## Run 349 stress 0.09495903 
    ## ... Procrustes: rmse 7.617672e-06  max resid 1.836341e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.09495903 
    ## ... Procrustes: rmse 5.59967e-06  max resid 1.37598e-05 
    ## ... Similar to previous best
    ## Run 351 stress 0.09495903 
    ## ... Procrustes: rmse 3.483059e-06  max resid 8.503956e-06 
    ## ... Similar to previous best
    ## Run 352 stress 0.09495903 
    ## ... Procrustes: rmse 3.132239e-06  max resid 7.616705e-06 
    ## ... Similar to previous best
    ## Run 353 stress 0.09495903 
    ## ... Procrustes: rmse 7.71639e-06  max resid 1.833181e-05 
    ## ... Similar to previous best
    ## Run 354 stress 0.09495903 
    ## ... Procrustes: rmse 3.536614e-06  max resid 8.707593e-06 
    ## ... Similar to previous best
    ## Run 355 stress 0.09495903 
    ## ... Procrustes: rmse 5.113035e-06  max resid 1.196799e-05 
    ## ... Similar to previous best
    ## Run 356 stress 0.1384219 
    ## Run 357 stress 0.1264617 
    ## Run 358 stress 0.09495903 
    ## ... Procrustes: rmse 3.125604e-06  max resid 7.525835e-06 
    ## ... Similar to previous best
    ## Run 359 stress 0.1264617 
    ## Run 360 stress 0.09495903 
    ## ... Procrustes: rmse 1.492422e-05  max resid 3.583445e-05 
    ## ... Similar to previous best
    ## Run 361 stress 0.1264617 
    ## Run 362 stress 0.1521592 
    ## Run 363 stress 0.09495903 
    ## ... Procrustes: rmse 2.84406e-06  max resid 3.852486e-06 
    ## ... Similar to previous best
    ## Run 364 stress 0.09495903 
    ## ... Procrustes: rmse 1.324876e-06  max resid 2.418925e-06 
    ## ... Similar to previous best
    ## Run 365 stress 0.1264617 
    ## Run 366 stress 0.09495903 
    ## ... Procrustes: rmse 9.78497e-07  max resid 2.388968e-06 
    ## ... Similar to previous best
    ## Run 367 stress 0.09495903 
    ## ... Procrustes: rmse 6.016787e-06  max resid 1.507227e-05 
    ## ... Similar to previous best
    ## Run 368 stress 0.09495903 
    ## ... Procrustes: rmse 2.906802e-06  max resid 7.139767e-06 
    ## ... Similar to previous best
    ## Run 369 stress 0.1264617 
    ## Run 370 stress 0.09495903 
    ## ... Procrustes: rmse 1.76554e-06  max resid 3.755019e-06 
    ## ... Similar to previous best
    ## Run 371 stress 0.09495903 
    ## ... Procrustes: rmse 4.562875e-06  max resid 1.082012e-05 
    ## ... Similar to previous best
    ## Run 372 stress 0.1264617 
    ## Run 373 stress 0.09495903 
    ## ... Procrustes: rmse 1.882389e-06  max resid 4.190557e-06 
    ## ... Similar to previous best
    ## Run 374 stress 0.1384221 
    ## Run 375 stress 0.09495903 
    ## ... Procrustes: rmse 3.852482e-06  max resid 9.075898e-06 
    ## ... Similar to previous best
    ## Run 376 stress 0.1264617 
    ## Run 377 stress 0.1384222 
    ## Run 378 stress 0.1264617 
    ## Run 379 stress 0.09495903 
    ## ... Procrustes: rmse 7.42728e-06  max resid 1.771985e-05 
    ## ... Similar to previous best
    ## Run 380 stress 0.1264617 
    ## Run 381 stress 0.09495903 
    ## ... Procrustes: rmse 3.027216e-06  max resid 8.235125e-06 
    ## ... Similar to previous best
    ## Run 382 stress 0.09495903 
    ## ... Procrustes: rmse 3.232044e-06  max resid 7.8818e-06 
    ## ... Similar to previous best
    ## Run 383 stress 0.09495903 
    ## ... Procrustes: rmse 2.373854e-06  max resid 6.293915e-06 
    ## ... Similar to previous best
    ## Run 384 stress 0.09495903 
    ## ... Procrustes: rmse 1.902228e-06  max resid 3.760477e-06 
    ## ... Similar to previous best
    ## Run 385 stress 0.09495903 
    ## ... Procrustes: rmse 4.472782e-06  max resid 1.050588e-05 
    ## ... Similar to previous best
    ## Run 386 stress 0.09495903 
    ## ... Procrustes: rmse 9.43012e-07  max resid 1.934308e-06 
    ## ... Similar to previous best
    ## Run 387 stress 0.09495903 
    ## ... Procrustes: rmse 2.534909e-06  max resid 6.241309e-06 
    ## ... Similar to previous best
    ## Run 388 stress 0.09495903 
    ## ... Procrustes: rmse 6.026364e-06  max resid 1.445697e-05 
    ## ... Similar to previous best
    ## Run 389 stress 0.09495903 
    ## ... Procrustes: rmse 1.26943e-05  max resid 2.157675e-05 
    ## ... Similar to previous best
    ## Run 390 stress 0.1264617 
    ## Run 391 stress 0.09495903 
    ## ... Procrustes: rmse 5.826782e-06  max resid 1.356523e-05 
    ## ... Similar to previous best
    ## Run 392 stress 0.1264617 
    ## Run 393 stress 0.1521592 
    ## Run 394 stress 0.09495903 
    ## ... Procrustes: rmse 2.068427e-06  max resid 4.384536e-06 
    ## ... Similar to previous best
    ## Run 395 stress 0.09495903 
    ## ... Procrustes: rmse 1.905279e-06  max resid 4.487726e-06 
    ## ... Similar to previous best
    ## Run 396 stress 0.09495903 
    ## ... Procrustes: rmse 3.450992e-06  max resid 8.318174e-06 
    ## ... Similar to previous best
    ## Run 397 stress 0.09495903 
    ## ... Procrustes: rmse 1.522877e-06  max resid 3.781815e-06 
    ## ... Similar to previous best
    ## Run 398 stress 0.09495903 
    ## ... Procrustes: rmse 3.82714e-06  max resid 9.157794e-06 
    ## ... Similar to previous best
    ## Run 399 stress 0.09495903 
    ## ... Procrustes: rmse 9.162523e-06  max resid 2.179972e-05 
    ## ... Similar to previous best
    ## Run 400 stress 0.1264617 
    ## Run 401 stress 0.09495903 
    ## ... Procrustes: rmse 3.635872e-06  max resid 8.659876e-06 
    ## ... Similar to previous best
    ## Run 402 stress 0.1264617 
    ## Run 403 stress 0.09495903 
    ## ... Procrustes: rmse 3.13421e-06  max resid 7.723861e-06 
    ## ... Similar to previous best
    ## Run 404 stress 0.09495903 
    ## ... Procrustes: rmse 1.11884e-06  max resid 2.684814e-06 
    ## ... Similar to previous best
    ## Run 405 stress 0.1264617 
    ## Run 406 stress 0.09495903 
    ## ... Procrustes: rmse 5.826561e-06  max resid 1.391699e-05 
    ## ... Similar to previous best
    ## Run 407 stress 0.09495903 
    ## ... Procrustes: rmse 6.202425e-07  max resid 1.377654e-06 
    ## ... Similar to previous best
    ## Run 408 stress 0.1264617 
    ## Run 409 stress 0.09495903 
    ## ... Procrustes: rmse 3.505746e-06  max resid 9.953219e-06 
    ## ... Similar to previous best
    ## Run 410 stress 0.09495903 
    ## ... Procrustes: rmse 4.328647e-06  max resid 1.038365e-05 
    ## ... Similar to previous best
    ## Run 411 stress 0.09495903 
    ## ... Procrustes: rmse 7.358606e-06  max resid 1.757149e-05 
    ## ... Similar to previous best
    ## Run 412 stress 0.09495903 
    ## ... Procrustes: rmse 1.77385e-06  max resid 4.156493e-06 
    ## ... Similar to previous best
    ## Run 413 stress 0.1384219 
    ## Run 414 stress 0.09495903 
    ## ... Procrustes: rmse 1.24452e-06  max resid 2.732375e-06 
    ## ... Similar to previous best
    ## Run 415 stress 0.09495903 
    ## ... Procrustes: rmse 1.067469e-05  max resid 2.541997e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.09495903 
    ## ... Procrustes: rmse 6.142933e-07  max resid 1.362336e-06 
    ## ... Similar to previous best
    ## Run 417 stress 0.09495903 
    ## ... Procrustes: rmse 3.512913e-07  max resid 5.423629e-07 
    ## ... Similar to previous best
    ## Run 418 stress 0.09495903 
    ## ... Procrustes: rmse 1.176651e-05  max resid 2.142788e-05 
    ## ... Similar to previous best
    ## Run 419 stress 0.09495903 
    ## ... Procrustes: rmse 1.437319e-06  max resid 3.402576e-06 
    ## ... Similar to previous best
    ## Run 420 stress 0.09495903 
    ## ... Procrustes: rmse 6.67129e-06  max resid 1.574516e-05 
    ## ... Similar to previous best
    ## Run 421 stress 0.09495903 
    ## ... Procrustes: rmse 2.171011e-06  max resid 3.90758e-06 
    ## ... Similar to previous best
    ## Run 422 stress 0.09495903 
    ## ... Procrustes: rmse 6.308864e-06  max resid 1.527421e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.09495903 
    ## ... Procrustes: rmse 4.55797e-06  max resid 1.0791e-05 
    ## ... Similar to previous best
    ## Run 424 stress 0.138422 
    ## Run 425 stress 0.1264617 
    ## Run 426 stress 0.09495903 
    ## ... Procrustes: rmse 4.009648e-06  max resid 9.579959e-06 
    ## ... Similar to previous best
    ## Run 427 stress 0.09495903 
    ## ... Procrustes: rmse 3.686925e-06  max resid 7.479658e-06 
    ## ... Similar to previous best
    ## Run 428 stress 0.1264617 
    ## Run 429 stress 0.09495903 
    ## ... Procrustes: rmse 2.83787e-06  max resid 7.281914e-06 
    ## ... Similar to previous best
    ## Run 430 stress 0.1264617 
    ## Run 431 stress 0.3106309 
    ## Run 432 stress 0.09495903 
    ## ... Procrustes: rmse 1.188402e-05  max resid 3.040316e-05 
    ## ... Similar to previous best
    ## Run 433 stress 0.3168452 
    ## Run 434 stress 0.09495903 
    ## ... Procrustes: rmse 4.502295e-06  max resid 8.887435e-06 
    ## ... Similar to previous best
    ## Run 435 stress 0.09495903 
    ## ... Procrustes: rmse 3.117884e-06  max resid 7.629502e-06 
    ## ... Similar to previous best
    ## Run 436 stress 0.1384219 
    ## Run 437 stress 0.1264617 
    ## Run 438 stress 0.09495903 
    ## ... Procrustes: rmse 5.987947e-06  max resid 1.443872e-05 
    ## ... Similar to previous best
    ## Run 439 stress 0.1264617 
    ## Run 440 stress 0.09495903 
    ## ... Procrustes: rmse 5.306957e-07  max resid 1.144502e-06 
    ## ... Similar to previous best
    ## Run 441 stress 0.09495903 
    ## ... Procrustes: rmse 6.093063e-06  max resid 1.455072e-05 
    ## ... Similar to previous best
    ## Run 442 stress 0.09495903 
    ## ... Procrustes: rmse 1.953246e-06  max resid 4.872718e-06 
    ## ... Similar to previous best
    ## Run 443 stress 0.1384221 
    ## Run 444 stress 0.09495903 
    ## ... Procrustes: rmse 1.311725e-06  max resid 3.1777e-06 
    ## ... Similar to previous best
    ## Run 445 stress 0.09495903 
    ## ... Procrustes: rmse 1.845612e-06  max resid 4.298777e-06 
    ## ... Similar to previous best
    ## Run 446 stress 0.09495903 
    ## ... Procrustes: rmse 2.058232e-06  max resid 4.874928e-06 
    ## ... Similar to previous best
    ## Run 447 stress 0.09495903 
    ## ... Procrustes: rmse 3.540714e-06  max resid 9.062521e-06 
    ## ... Similar to previous best
    ## Run 448 stress 0.1264617 
    ## Run 449 stress 0.09495903 
    ## ... Procrustes: rmse 3.219964e-06  max resid 7.710898e-06 
    ## ... Similar to previous best
    ## Run 450 stress 0.09495903 
    ## ... Procrustes: rmse 9.449842e-07  max resid 1.821338e-06 
    ## ... Similar to previous best
    ## Run 451 stress 0.09495903 
    ## ... Procrustes: rmse 4.233528e-06  max resid 9.994881e-06 
    ## ... Similar to previous best
    ## Run 452 stress 0.09495903 
    ## ... Procrustes: rmse 8.809056e-07  max resid 1.587188e-06 
    ## ... Similar to previous best
    ## Run 453 stress 0.1384221 
    ## Run 454 stress 0.09495903 
    ## ... Procrustes: rmse 2.314663e-06  max resid 5.664527e-06 
    ## ... Similar to previous best
    ## Run 455 stress 0.1264617 
    ## Run 456 stress 0.1264617 
    ## Run 457 stress 0.09495903 
    ## ... Procrustes: rmse 2.622171e-06  max resid 5.969013e-06 
    ## ... Similar to previous best
    ## Run 458 stress 0.1264617 
    ## Run 459 stress 0.09495903 
    ## ... Procrustes: rmse 8.1109e-07  max resid 1.961411e-06 
    ## ... Similar to previous best
    ## Run 460 stress 0.09495903 
    ## ... Procrustes: rmse 5.034878e-06  max resid 1.200137e-05 
    ## ... Similar to previous best
    ## Run 461 stress 0.09495903 
    ## ... Procrustes: rmse 1.490993e-06  max resid 3.225116e-06 
    ## ... Similar to previous best
    ## Run 462 stress 0.09495903 
    ## ... Procrustes: rmse 4.281984e-06  max resid 1.051544e-05 
    ## ... Similar to previous best
    ## Run 463 stress 0.09495903 
    ## ... Procrustes: rmse 5.90598e-06  max resid 1.437792e-05 
    ## ... Similar to previous best
    ## Run 464 stress 0.09495903 
    ## ... Procrustes: rmse 5.91256e-07  max resid 1.196938e-06 
    ## ... Similar to previous best
    ## Run 465 stress 0.09495903 
    ## ... Procrustes: rmse 9.765137e-07  max resid 2.429266e-06 
    ## ... Similar to previous best
    ## Run 466 stress 0.09495903 
    ## ... Procrustes: rmse 8.966304e-06  max resid 2.139038e-05 
    ## ... Similar to previous best
    ## Run 467 stress 0.09495903 
    ## ... Procrustes: rmse 3.681993e-06  max resid 8.945778e-06 
    ## ... Similar to previous best
    ## Run 468 stress 0.138422 
    ## Run 469 stress 0.1264617 
    ## Run 470 stress 0.1264617 
    ## Run 471 stress 0.09495903 
    ## ... Procrustes: rmse 5.388925e-06  max resid 1.175494e-05 
    ## ... Similar to previous best
    ## Run 472 stress 0.09495903 
    ## ... Procrustes: rmse 4.706024e-06  max resid 9.358765e-06 
    ## ... Similar to previous best
    ## Run 473 stress 0.09495903 
    ## ... Procrustes: rmse 2.700098e-06  max resid 6.296089e-06 
    ## ... Similar to previous best
    ## Run 474 stress 0.09495903 
    ## ... Procrustes: rmse 1.210498e-06  max resid 2.958733e-06 
    ## ... Similar to previous best
    ## Run 475 stress 0.1264617 
    ## Run 476 stress 0.09495903 
    ## ... Procrustes: rmse 6.739197e-06  max resid 1.256282e-05 
    ## ... Similar to previous best
    ## Run 477 stress 0.1264617 
    ## Run 478 stress 0.1264617 
    ## Run 479 stress 0.1264617 
    ## Run 480 stress 0.09495903 
    ## ... Procrustes: rmse 1.456361e-06  max resid 3.763674e-06 
    ## ... Similar to previous best
    ## Run 481 stress 0.09495903 
    ## ... Procrustes: rmse 6.345354e-07  max resid 1.465017e-06 
    ## ... Similar to previous best
    ## Run 482 stress 0.1264617 
    ## Run 483 stress 0.09495903 
    ## ... Procrustes: rmse 9.949804e-06  max resid 2.352026e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.09495903 
    ## ... Procrustes: rmse 2.883442e-06  max resid 8.199665e-06 
    ## ... Similar to previous best
    ## Run 485 stress 0.09495903 
    ## ... Procrustes: rmse 2.283246e-06  max resid 5.607689e-06 
    ## ... Similar to previous best
    ## Run 486 stress 0.1264617 
    ## Run 487 stress 0.1264617 
    ## Run 488 stress 0.09495903 
    ## ... Procrustes: rmse 2.282655e-06  max resid 5.611618e-06 
    ## ... Similar to previous best
    ## Run 489 stress 0.09495903 
    ## ... Procrustes: rmse 2.410703e-06  max resid 5.73537e-06 
    ## ... Similar to previous best
    ## Run 490 stress 0.09495903 
    ## ... Procrustes: rmse 6.336408e-06  max resid 1.528345e-05 
    ## ... Similar to previous best
    ## Run 491 stress 0.09495903 
    ## ... Procrustes: rmse 2.440713e-06  max resid 6.000989e-06 
    ## ... Similar to previous best
    ## Run 492 stress 0.09495903 
    ## ... Procrustes: rmse 1.966172e-06  max resid 4.692825e-06 
    ## ... Similar to previous best
    ## Run 493 stress 0.1264617 
    ## Run 494 stress 0.09495903 
    ## ... Procrustes: rmse 6.149342e-06  max resid 1.434691e-05 
    ## ... Similar to previous best
    ## Run 495 stress 0.09495903 
    ## ... Procrustes: rmse 8.223437e-06  max resid 1.914088e-05 
    ## ... Similar to previous best
    ## Run 496 stress 0.09495903 
    ## ... Procrustes: rmse 3.64619e-06  max resid 9.447244e-06 
    ## ... Similar to previous best
    ## Run 497 stress 0.09495903 
    ## ... Procrustes: rmse 1.710671e-06  max resid 4.171021e-06 
    ## ... Similar to previous best
    ## Run 498 stress 0.09495903 
    ## ... Procrustes: rmse 2.099014e-06  max resid 4.397152e-06 
    ## ... Similar to previous best
    ## Run 499 stress 0.1264617 
    ## Run 500 stress 0.09495903 
    ## ... Procrustes: rmse 2.993582e-06  max resid 7.628662e-06 
    ## ... Similar to previous best
    ## *** Best solution repeated 227 times

``` r
# Stratified lakes and ocean sites
PD_beta_geo_SO_NMDS <- metaMDS(PD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.04851417 
    ## Run 1 stress 0.05959896 
    ## Run 2 stress 0.05731513 
    ## Run 3 stress 0.05728055 
    ## Run 4 stress 0.05572688 
    ## Run 5 stress 0.05276414 
    ## Run 6 stress 0.04612761 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04982063  max resid 0.1395299 
    ## Run 7 stress 0.04935951 
    ## Run 8 stress 0.04612763 
    ## ... Procrustes: rmse 0.0002015856  max resid 0.0003999025 
    ## ... Similar to previous best
    ## Run 9 stress 0.05607345 
    ## Run 10 stress 0.04860423 
    ## Run 11 stress 0.06088971 
    ## Run 12 stress 0.05643263 
    ## Run 13 stress 0.0560516 
    ## Run 14 stress 0.0527641 
    ## Run 15 stress 0.05260047 
    ## Run 16 stress 0.05593459 
    ## Run 17 stress 0.04938618 
    ## Run 18 stress 0.05612241 
    ## Run 19 stress 0.05474814 
    ## Run 20 stress 0.04860423 
    ## Run 21 stress 0.04935943 
    ## Run 22 stress 0.0559346 
    ## Run 23 stress 0.05605147 
    ## Run 24 stress 0.0557269 
    ## Run 25 stress 0.04851417 
    ## Run 26 stress 0.04612766 
    ## ... Procrustes: rmse 6.031795e-05  max resid 9.964202e-05 
    ## ... Similar to previous best
    ## Run 27 stress 0.05607344 
    ## Run 28 stress 0.05728057 
    ## Run 29 stress 0.05276412 
    ## Run 30 stress 0.06291426 
    ## Run 31 stress 0.04851419 
    ## Run 32 stress 0.05276416 
    ## Run 33 stress 0.05612256 
    ## Run 34 stress 0.04935946 
    ## Run 35 stress 0.04851419 
    ## Run 36 stress 0.05695857 
    ## Run 37 stress 0.04851417 
    ## Run 38 stress 0.05728054 
    ## Run 39 stress 0.05643262 
    ## Run 40 stress 0.04938616 
    ## Run 41 stress 0.05643262 
    ## Run 42 stress 0.05605149 
    ## Run 43 stress 0.05216165 
    ## Run 44 stress 0.04851419 
    ## Run 45 stress 0.05260054 
    ## Run 46 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 6.659111e-05  max resid 0.0001133602 
    ## ... Similar to previous best
    ## Run 47 stress 0.04938616 
    ## Run 48 stress 0.05439282 
    ## Run 49 stress 0.0629143 
    ## Run 50 stress 0.05728065 
    ## Run 51 stress 0.04938618 
    ## Run 52 stress 0.04851417 
    ## Run 53 stress 0.04612762 
    ## ... Procrustes: rmse 0.0001411833  max resid 0.0002871082 
    ## ... Similar to previous best
    ## Run 54 stress 0.06400728 
    ## Run 55 stress 0.05216151 
    ## Run 56 stress 0.0527641 
    ## Run 57 stress 0.05643261 
    ## Run 58 stress 0.04851417 
    ## Run 59 stress 0.04851419 
    ## Run 60 stress 0.04860423 
    ## Run 61 stress 0.05728071 
    ## Run 62 stress 0.05137973 
    ## Run 63 stress 0.05855033 
    ## Run 64 stress 0.06291428 
    ## Run 65 stress 0.05474815 
    ## Run 66 stress 0.05855029 
    ## Run 67 stress 0.05438954 
    ## Run 68 stress 0.05643263 
    ## Run 69 stress 0.04938616 
    ## Run 70 stress 0.04938616 
    ## Run 71 stress 0.05612242 
    ## Run 72 stress 0.04938619 
    ## Run 73 stress 0.05607344 
    ## Run 74 stress 0.05216161 
    ## Run 75 stress 0.05728063 
    ## Run 76 stress 0.05695889 
    ## Run 77 stress 0.05216137 
    ## Run 78 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001600011  max resid 0.0003235312 
    ## ... Similar to previous best
    ## Run 79 stress 0.05728062 
    ## Run 80 stress 0.05572688 
    ## Run 81 stress 0.05728057 
    ## Run 82 stress 0.05842809 
    ## Run 83 stress 0.04851417 
    ## Run 84 stress 0.05260049 
    ## Run 85 stress 0.04860423 
    ## Run 86 stress 0.04938615 
    ## Run 87 stress 0.05439283 
    ## Run 88 stress 0.05451493 
    ## Run 89 stress 0.04938615 
    ## Run 90 stress 0.05216179 
    ## Run 91 stress 0.05216126 
    ## Run 92 stress 0.06400746 
    ## Run 93 stress 0.05643263 
    ## Run 94 stress 0.05216152 
    ## Run 95 stress 0.06185053 
    ## Run 96 stress 0.0485142 
    ## Run 97 stress 0.05695883 
    ## Run 98 stress 0.05260046 
    ## Run 99 stress 0.04938615 
    ## Run 100 stress 0.04860422 
    ## Run 101 stress 0.04851417 
    ## Run 102 stress 0.05216157 
    ## Run 103 stress 0.05451494 
    ## Run 104 stress 0.05137963 
    ## Run 105 stress 0.0640074 
    ## Run 106 stress 0.04938615 
    ## Run 107 stress 0.05605147 
    ## Run 108 stress 0.05439281 
    ## Run 109 stress 0.05439282 
    ## Run 110 stress 0.05438936 
    ## Run 111 stress 0.05643261 
    ## Run 112 stress 0.0513797 
    ## Run 113 stress 0.04935947 
    ## Run 114 stress 0.05607344 
    ## Run 115 stress 0.05438955 
    ## Run 116 stress 0.04851418 
    ## Run 117 stress 0.05488104 
    ## Run 118 stress 0.05605155 
    ## Run 119 stress 0.0543894 
    ## Run 120 stress 0.04938616 
    ## Run 121 stress 0.06185053 
    ## Run 122 stress 0.05728055 
    ## Run 123 stress 0.04938616 
    ## Run 124 stress 0.04935944 
    ## Run 125 stress 0.05728077 
    ## Run 126 stress 0.0527641 
    ## Run 127 stress 0.05855031 
    ## Run 128 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 5.299066e-05  max resid 0.0001065521 
    ## ... Similar to previous best
    ## Run 129 stress 0.0561225 
    ## Run 130 stress 0.05612249 
    ## Run 131 stress 0.05728071 
    ## Run 132 stress 0.04860423 
    ## Run 133 stress 0.06205752 
    ## Run 134 stress 0.05572688 
    ## Run 135 stress 0.04612759 
    ## ... Procrustes: rmse 5.558245e-05  max resid 0.0001112033 
    ## ... Similar to previous best
    ## Run 136 stress 0.05607344 
    ## Run 137 stress 0.05260049 
    ## Run 138 stress 0.04938616 
    ## Run 139 stress 0.05474814 
    ## Run 140 stress 0.04851418 
    ## Run 141 stress 0.05855031 
    ## Run 142 stress 0.05451493 
    ## Run 143 stress 0.05260048 
    ## Run 144 stress 0.04612758 
    ## ... Procrustes: rmse 8.707657e-05  max resid 0.0001776244 
    ## ... Similar to previous best
    ## Run 145 stress 0.05728053 
    ## Run 146 stress 0.04938618 
    ## Run 147 stress 0.05607348 
    ## Run 148 stress 0.04851418 
    ## Run 149 stress 0.04851418 
    ## Run 150 stress 0.04851418 
    ## Run 151 stress 0.05695857 
    ## Run 152 stress 0.04938617 
    ## Run 153 stress 0.05731512 
    ## Run 154 stress 0.06185049 
    ## Run 155 stress 0.04938616 
    ## Run 156 stress 0.05855029 
    ## Run 157 stress 0.05855029 
    ## Run 158 stress 0.0526005 
    ## Run 159 stress 0.05643261 
    ## Run 160 stress 0.3257554 
    ## Run 161 stress 0.0527641 
    ## Run 162 stress 0.05607344 
    ## Run 163 stress 0.0493595 
    ## Run 164 stress 0.05728053 
    ## Run 165 stress 0.05607347 
    ## Run 166 stress 0.06185052 
    ## Run 167 stress 0.04938619 
    ## Run 168 stress 0.05728056 
    ## Run 169 stress 0.05137965 
    ## Run 170 stress 0.04851417 
    ## Run 171 stress 0.05593465 
    ## Run 172 stress 0.04851417 
    ## Run 173 stress 0.05605154 
    ## Run 174 stress 0.04935953 
    ## Run 175 stress 0.04938616 
    ## Run 176 stress 0.04935951 
    ## Run 177 stress 0.05474815 
    ## Run 178 stress 0.04851417 
    ## Run 179 stress 0.05842817 
    ## Run 180 stress 0.0521613 
    ## Run 181 stress 0.05451497 
    ## Run 182 stress 0.06205754 
    ## Run 183 stress 0.05260054 
    ## Run 184 stress 0.0543895 
    ## Run 185 stress 0.05438952 
    ## Run 186 stress 0.0543928 
    ## Run 187 stress 0.05438948 
    ## Run 188 stress 0.0493862 
    ## Run 189 stress 0.04851418 
    ## Run 190 stress 0.04860424 
    ## Run 191 stress 0.059599 
    ## Run 192 stress 0.05438941 
    ## Run 193 stress 0.04935949 
    ## Run 194 stress 0.06400754 
    ## Run 195 stress 0.05438941 
    ## Run 196 stress 0.04935955 
    ## Run 197 stress 0.04851417 
    ## Run 198 stress 0.04851417 
    ## Run 199 stress 0.0561224 
    ## Run 200 stress 0.04938615 
    ## Run 201 stress 0.04851418 
    ## Run 202 stress 0.04938616 
    ## Run 203 stress 0.05451494 
    ## Run 204 stress 0.0543928 
    ## Run 205 stress 0.05695877 
    ## Run 206 stress 0.04938616 
    ## Run 207 stress 0.04851417 
    ## Run 208 stress 0.05433173 
    ## Run 209 stress 0.05474814 
    ## Run 210 stress 0.04938615 
    ## Run 211 stress 0.05593459 
    ## Run 212 stress 0.04612759 
    ## ... Procrustes: rmse 4.474595e-05  max resid 8.974817e-05 
    ## ... Similar to previous best
    ## Run 213 stress 0.04612758 
    ## ... Procrustes: rmse 2.669201e-05  max resid 4.557727e-05 
    ## ... Similar to previous best
    ## Run 214 stress 0.05842829 
    ## Run 215 stress 0.05855033 
    ## Run 216 stress 0.05643271 
    ## Run 217 stress 0.0521613 
    ## Run 218 stress 0.05728053 
    ## Run 219 stress 0.0572806 
    ## Run 220 stress 0.05607344 
    ## Run 221 stress 0.04612765 
    ## ... Procrustes: rmse 0.0001601455  max resid 0.0003232421 
    ## ... Similar to previous best
    ## Run 222 stress 0.06400733 
    ## Run 223 stress 0.06088971 
    ## Run 224 stress 0.05607346 
    ## Run 225 stress 0.05593468 
    ## Run 226 stress 0.05728054 
    ## Run 227 stress 0.05612241 
    ## Run 228 stress 0.06185054 
    ## Run 229 stress 0.06400759 
    ## Run 230 stress 0.06205752 
    ## Run 231 stress 0.04851417 
    ## Run 232 stress 0.05607345 
    ## Run 233 stress 0.05607345 
    ## Run 234 stress 0.05488106 
    ## Run 235 stress 0.05695853 
    ## Run 236 stress 0.0521615 
    ## Run 237 stress 0.05607344 
    ## Run 238 stress 0.04612758 
    ## ... Procrustes: rmse 5.929353e-05  max resid 0.0001120689 
    ## ... Similar to previous best
    ## Run 239 stress 0.05572688 
    ## Run 240 stress 0.05605148 
    ## Run 241 stress 0.05959896 
    ## Run 242 stress 0.05643262 
    ## Run 243 stress 0.05959895 
    ## Run 244 stress 0.0569586 
    ## Run 245 stress 0.05607347 
    ## Run 246 stress 0.059599 
    ## Run 247 stress 0.04612758 
    ## ... Procrustes: rmse 3.367494e-05  max resid 6.735865e-05 
    ## ... Similar to previous best
    ## Run 248 stress 0.05695883 
    ## Run 249 stress 0.04938616 
    ## Run 250 stress 0.05137966 
    ## Run 251 stress 0.05842827 
    ## Run 252 stress 0.0547482 
    ## Run 253 stress 0.05695871 
    ## Run 254 stress 0.05260049 
    ## Run 255 stress 0.05607347 
    ## Run 256 stress 0.05842802 
    ## Run 257 stress 0.0543928 
    ## Run 258 stress 0.06205755 
    ## Run 259 stress 0.05605155 
    ## Run 260 stress 0.04935945 
    ## Run 261 stress 0.04612759 
    ## ... Procrustes: rmse 9.985581e-05  max resid 0.0002031251 
    ## ... Similar to previous best
    ## Run 262 stress 0.06088972 
    ## Run 263 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001553854  max resid 0.0003127105 
    ## ... Similar to previous best
    ## Run 264 stress 0.05260049 
    ## Run 265 stress 0.05438945 
    ## Run 266 stress 0.04851417 
    ## Run 267 stress 0.06185053 
    ## Run 268 stress 0.05728057 
    ## Run 269 stress 0.05593458 
    ## Run 270 stress 0.06205753 
    ## Run 271 stress 0.05276412 
    ## Run 272 stress 0.04851418 
    ## Run 273 stress 0.06185045 
    ## Run 274 stress 0.05216165 
    ## Run 275 stress 0.04935944 
    ## Run 276 stress 0.05474819 
    ## Run 277 stress 0.04938616 
    ## Run 278 stress 0.04851417 
    ## Run 279 stress 0.05216134 
    ## Run 280 stress 0.04851417 
    ## Run 281 stress 0.05474815 
    ## Run 282 stress 0.05474817 
    ## Run 283 stress 0.04935952 
    ## Run 284 stress 0.05260048 
    ## Run 285 stress 0.04860423 
    ## Run 286 stress 0.04938616 
    ## Run 287 stress 0.04851417 
    ## Run 288 stress 0.06185067 
    ## Run 289 stress 0.05728063 
    ## Run 290 stress 0.04938615 
    ## Run 291 stress 0.05137964 
    ## Run 292 stress 0.06205753 
    ## Run 293 stress 0.05607345 
    ## Run 294 stress 0.06291435 
    ## Run 295 stress 0.07046361 
    ## Run 296 stress 0.0561224 
    ## Run 297 stress 0.05728055 
    ## Run 298 stress 0.05488111 
    ## Run 299 stress 0.06291431 
    ## Run 300 stress 0.04938621 
    ## Run 301 stress 0.05612247 
    ## Run 302 stress 0.0461277 
    ## ... Procrustes: rmse 0.0001536691  max resid 0.0003134713 
    ## ... Similar to previous best
    ## Run 303 stress 0.04851417 
    ## Run 304 stress 0.04851418 
    ## Run 305 stress 0.0527641 
    ## Run 306 stress 0.04938616 
    ## Run 307 stress 0.05728067 
    ## Run 308 stress 0.05438943 
    ## Run 309 stress 0.05451497 
    ## Run 310 stress 0.0527641 
    ## Run 311 stress 0.04938616 
    ## Run 312 stress 0.05488103 
    ## Run 313 stress 0.05607345 
    ## Run 314 stress 0.04612758 
    ## ... Procrustes: rmse 7.242677e-05  max resid 0.0001493829 
    ## ... Similar to previous best
    ## Run 315 stress 0.3578429 
    ## Run 316 stress 0.05474815 
    ## Run 317 stress 0.04938615 
    ## Run 318 stress 0.0543928 
    ## Run 319 stress 0.05439287 
    ## Run 320 stress 0.06205756 
    ## Run 321 stress 0.05260047 
    ## Run 322 stress 0.05451493 
    ## Run 323 stress 0.05855029 
    ## Run 324 stress 0.05643265 
    ## Run 325 stress 0.05605156 
    ## Run 326 stress 0.06291425 
    ## Run 327 stress 0.04938617 
    ## Run 328 stress 0.05607345 
    ## Run 329 stress 0.04935949 
    ## Run 330 stress 0.0560515 
    ## Run 331 stress 0.04935943 
    ## Run 332 stress 0.06291424 
    ## Run 333 stress 0.05260049 
    ## Run 334 stress 0.05607344 
    ## Run 335 stress 0.06088973 
    ## Run 336 stress 0.05216155 
    ## Run 337 stress 0.05643262 
    ## Run 338 stress 0.05695872 
    ## Run 339 stress 0.04851417 
    ## Run 340 stress 0.0543928 
    ## Run 341 stress 0.05439283 
    ## Run 342 stress 0.05612255 
    ## Run 343 stress 0.04938616 
    ## Run 344 stress 0.04851417 
    ## Run 345 stress 0.04851419 
    ## Run 346 stress 0.04938615 
    ## Run 347 stress 0.05260052 
    ## Run 348 stress 0.04612758 
    ## ... Procrustes: rmse 7.416153e-05  max resid 0.0001512788 
    ## ... Similar to previous best
    ## Run 349 stress 0.05607344 
    ## Run 350 stress 0.0559345 
    ## Run 351 stress 0.05605155 
    ## Run 352 stress 0.05451495 
    ## Run 353 stress 0.05593449 
    ## Run 354 stress 0.05607346 
    ## Run 355 stress 0.04851417 
    ## Run 356 stress 0.04860422 
    ## Run 357 stress 0.05842817 
    ## Run 358 stress 0.04612759 
    ## ... Procrustes: rmse 4.40118e-05  max resid 8.625212e-05 
    ## ... Similar to previous best
    ## Run 359 stress 0.0461276 
    ## ... Procrustes: rmse 0.0001094317  max resid 0.0002236564 
    ## ... Similar to previous best
    ## Run 360 stress 0.05438938 
    ## Run 361 stress 0.05607346 
    ## Run 362 stress 0.05137968 
    ## Run 363 stress 0.05438945 
    ## Run 364 stress 0.05438939 
    ## Run 365 stress 0.04851418 
    ## Run 366 stress 0.05216144 
    ## Run 367 stress 0.04860422 
    ## Run 368 stress 0.05137969 
    ## Run 369 stress 0.04935943 
    ## Run 370 stress 0.05593455 
    ## Run 371 stress 0.0608897 
    ## Run 372 stress 0.3348997 
    ## Run 373 stress 0.04938618 
    ## Run 374 stress 0.04938616 
    ## Run 375 stress 0.05959898 
    ## Run 376 stress 0.05643261 
    ## Run 377 stress 0.05438944 
    ## Run 378 stress 0.05605152 
    ## Run 379 stress 0.05438951 
    ## Run 380 stress 0.05728066 
    ## Run 381 stress 0.05695884 
    ## Run 382 stress 0.04851418 
    ## Run 383 stress 0.05728057 
    ## Run 384 stress 0.05260048 
    ## Run 385 stress 0.05607345 
    ## Run 386 stress 0.05451493 
    ## Run 387 stress 0.0585503 
    ## Run 388 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 4.968043e-06  max resid 9.309495e-06 
    ## ... Similar to previous best
    ## Run 389 stress 0.04612759 
    ## ... Procrustes: rmse 8.87973e-05  max resid 0.0001767372 
    ## ... Similar to previous best
    ## Run 390 stress 0.04938616 
    ## Run 391 stress 0.04851421 
    ## Run 392 stress 0.05451494 
    ## Run 393 stress 0.05605158 
    ## Run 394 stress 0.06400742 
    ## Run 395 stress 0.04938615 
    ## Run 396 stress 0.05842787 
    ## Run 397 stress 0.05260051 
    ## Run 398 stress 0.04851417 
    ## Run 399 stress 0.0543929 
    ## Run 400 stress 0.05643261 
    ## Run 401 stress 0.0572806 
    ## Run 402 stress 0.06088973 
    ## Run 403 stress 0.05695874 
    ## Run 404 stress 0.05474816 
    ## Run 405 stress 0.05855031 
    ## Run 406 stress 0.04938618 
    ## Run 407 stress 0.05572693 
    ## Run 408 stress 0.05728064 
    ## Run 409 stress 0.05216159 
    ## Run 410 stress 0.05438944 
    ## Run 411 stress 0.04860422 
    ## Run 412 stress 0.05451494 
    ## Run 413 stress 0.05593466 
    ## Run 414 stress 0.05612244 
    ## Run 415 stress 0.04938616 
    ## Run 416 stress 0.04851418 
    ## Run 417 stress 0.04851418 
    ## Run 418 stress 0.04851417 
    ## Run 419 stress 0.05605149 
    ## Run 420 stress 0.04938615 
    ## Run 421 stress 0.06185047 
    ## Run 422 stress 0.07046373 
    ## Run 423 stress 0.0629143 
    ## Run 424 stress 0.06205756 
    ## Run 425 stress 0.04612757 
    ## ... Procrustes: rmse 5.245538e-05  max resid 0.0001063533 
    ## ... Similar to previous best
    ## Run 426 stress 0.05855033 
    ## Run 427 stress 0.05451493 
    ## Run 428 stress 0.05728053 
    ## Run 429 stress 0.05607346 
    ## Run 430 stress 0.04938616 
    ## Run 431 stress 0.05572694 
    ## Run 432 stress 0.05605151 
    ## Run 433 stress 0.05728053 
    ## Run 434 stress 0.05474816 
    ## Run 435 stress 0.05695882 
    ## Run 436 stress 0.0559346 
    ## Run 437 stress 0.04851417 
    ## Run 438 stress 0.05695894 
    ## Run 439 stress 0.06205753 
    ## Run 440 stress 0.04851418 
    ## Run 441 stress 0.04851417 
    ## Run 442 stress 0.05276412 
    ## Run 443 stress 0.05593456 
    ## Run 444 stress 0.05439284 
    ## Run 445 stress 0.05605151 
    ## Run 446 stress 0.04851417 
    ## Run 447 stress 0.05260048 
    ## Run 448 stress 0.04938616 
    ## Run 449 stress 0.05643261 
    ## Run 450 stress 0.05451493 
    ## Run 451 stress 0.0461276 
    ## ... Procrustes: rmse 0.0001138627  max resid 0.000230799 
    ## ... Similar to previous best
    ## Run 452 stress 0.04938619 
    ## Run 453 stress 0.04938619 
    ## Run 454 stress 0.05643269 
    ## Run 455 stress 0.05728054 
    ## Run 456 stress 0.0543318 
    ## Run 457 stress 0.04851417 
    ## Run 458 stress 0.05959896 
    ## Run 459 stress 0.04851417 
    ## Run 460 stress 0.04851417 
    ## Run 461 stress 0.04938616 
    ## Run 462 stress 0.04938617 
    ## Run 463 stress 0.04938616 
    ## Run 464 stress 0.05842786 
    ## Run 465 stress 0.05474815 
    ## Run 466 stress 0.05728056 
    ## Run 467 stress 0.06205754 
    ## Run 468 stress 0.05451493 
    ## Run 469 stress 0.05447786 
    ## Run 470 stress 0.04851417 
    ## Run 471 stress 0.05488103 
    ## Run 472 stress 0.04938615 
    ## Run 473 stress 0.0461276 
    ## ... Procrustes: rmse 0.0001041956  max resid 0.0002121098 
    ## ... Similar to previous best
    ## Run 474 stress 0.04938616 
    ## Run 475 stress 0.04612763 
    ## ... Procrustes: rmse 0.000127061  max resid 0.0002555236 
    ## ... Similar to previous best
    ## Run 476 stress 0.05607347 
    ## Run 477 stress 0.05474816 
    ## Run 478 stress 0.05451494 
    ## Run 479 stress 0.04935947 
    ## Run 480 stress 0.05260054 
    ## Run 481 stress 0.05959901 
    ## Run 482 stress 0.05137965 
    ## Run 483 stress 0.06205752 
    ## Run 484 stress 0.05607344 
    ## Run 485 stress 0.05216135 
    ## Run 486 stress 0.04935951 
    ## Run 487 stress 0.04860422 
    ## Run 488 stress 0.05607347 
    ## Run 489 stress 0.04851417 
    ## Run 490 stress 0.05451493 
    ## Run 491 stress 0.04935944 
    ## Run 492 stress 0.05855029 
    ## Run 493 stress 0.04935942 
    ## Run 494 stress 0.05855029 
    ## Run 495 stress 0.05439286 
    ## Run 496 stress 0.05474813 
    ## Run 497 stress 0.05607348 
    ## Run 498 stress 0.05728058 
    ## Run 499 stress 0.05959898 
    ## Run 500 stress 0.05643263 
    ## *** Best solution repeated 6 times

``` r
# Mixed lakes
PD_beta_geo_M_NMDS <- metaMDS(PD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.953481e-05 
    ## Run 1 stress 9.869257e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.202134  max resid 0.318843 
    ## Run 2 stress 8.996246e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1978375  max resid 0.3255758 
    ## Run 3 stress 9.435038e-05 
    ## ... Procrustes: rmse 0.1978758  max resid 0.3547605 
    ## Run 4 stress 9.029539e-05 
    ## ... Procrustes: rmse 0.0205847  max resid 0.04605052 
    ## Run 5 stress 9.818079e-05 
    ## ... Procrustes: rmse 0.06116588  max resid 0.1308831 
    ## Run 6 stress 9.836492e-05 
    ## ... Procrustes: rmse 0.03383692  max resid 0.05677174 
    ## Run 7 stress 9.674111e-05 
    ## ... Procrustes: rmse 0.1051452  max resid 0.168685 
    ## Run 8 stress 0.2570537 
    ## Run 9 stress 8.745797e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.06922479  max resid 0.1465717 
    ## Run 10 stress 9.63136e-05 
    ## ... Procrustes: rmse 0.2345418  max resid 0.3315968 
    ## Run 11 stress 9.967046e-05 
    ## ... Procrustes: rmse 0.09528764  max resid 0.1767152 
    ## Run 12 stress 9.917994e-05 
    ## ... Procrustes: rmse 0.05543524  max resid 0.07978323 
    ## Run 13 stress 9.451289e-05 
    ## ... Procrustes: rmse 0.06206284  max resid 0.09919253 
    ## Run 14 stress 9.358828e-05 
    ## ... Procrustes: rmse 0.1680007  max resid 0.2520839 
    ## Run 15 stress 0.2570537 
    ## Run 16 stress 9.634975e-05 
    ## ... Procrustes: rmse 0.04564588  max resid 0.07772253 
    ## Run 17 stress 9.651849e-05 
    ## ... Procrustes: rmse 0.04679628  max resid 0.07094662 
    ## Run 18 stress 8.98925e-05 
    ## ... Procrustes: rmse 0.04276556  max resid 0.07528706 
    ## Run 19 stress 9.809901e-05 
    ## ... Procrustes: rmse 0.07190767  max resid 0.1112987 
    ## Run 20 stress 9.6388e-05 
    ## ... Procrustes: rmse 0.04274003  max resid 0.08861251 
    ## Run 21 stress 9.895304e-05 
    ## ... Procrustes: rmse 0.110294  max resid 0.1465314 
    ## Run 22 stress 9.347527e-05 
    ## ... Procrustes: rmse 0.1104402  max resid 0.1893554 
    ## Run 23 stress 8.8264e-05 
    ## ... Procrustes: rmse 0.2351022  max resid 0.3325901 
    ## Run 24 stress 9.583716e-05 
    ## ... Procrustes: rmse 0.0529213  max resid 0.08723266 
    ## Run 25 stress 8.923641e-05 
    ## ... Procrustes: rmse 0.2351063  max resid 0.3327724 
    ## Run 26 stress 9.736827e-05 
    ## ... Procrustes: rmse 0.2351062  max resid 0.3327673 
    ## Run 27 stress 0.2290535 
    ## Run 28 stress 9.357511e-05 
    ## ... Procrustes: rmse 0.04725686  max resid 0.07982909 
    ## Run 29 stress 9.185206e-05 
    ## ... Procrustes: rmse 0.2350986  max resid 0.3325788 
    ## Run 30 stress 9.859513e-05 
    ## ... Procrustes: rmse 0.05070575  max resid 0.08933033 
    ## Run 31 stress 9.460612e-05 
    ## ... Procrustes: rmse 0.05715584  max resid 0.1226592 
    ## Run 32 stress 9.873372e-05 
    ## ... Procrustes: rmse 0.2351015  max resid 0.3325674 
    ## Run 33 stress 0.2570537 
    ## Run 34 stress 9.698968e-05 
    ## ... Procrustes: rmse 0.03302769  max resid 0.06738822 
    ## Run 35 stress 8.841807e-05 
    ## ... Procrustes: rmse 0.2350979  max resid 0.3325817 
    ## Run 36 stress 0.2570537 
    ## Run 37 stress 8.636205e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2351064  max resid 0.3327771 
    ## Run 38 stress 8.375162e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001992302  max resid 0.0003622008 
    ## ... Similar to previous best
    ## Run 39 stress 9.857514e-05 
    ## ... Procrustes: rmse 0.0760133  max resid 0.1465856 
    ## Run 40 stress 0.2570537 
    ## Run 41 stress 8.410022e-05 
    ## ... Procrustes: rmse 0.1849761  max resid 0.2275304 
    ## Run 42 stress 9.355322e-05 
    ## ... Procrustes: rmse 0.03306617  max resid 0.06124728 
    ## Run 43 stress 9.565801e-05 
    ## ... Procrustes: rmse 0.131667  max resid 0.179074 
    ## Run 44 stress 0.2290618 
    ## Run 45 stress 9.597655e-05 
    ## ... Procrustes: rmse 0.1410786  max resid 0.26175 
    ## Run 46 stress 9.64911e-05 
    ## ... Procrustes: rmse 0.1940298  max resid 0.2621917 
    ## Run 47 stress 9.563016e-05 
    ## ... Procrustes: rmse 0.1687604  max resid 0.3029135 
    ## Run 48 stress 9.923418e-05 
    ## ... Procrustes: rmse 0.03324785  max resid 0.04484771 
    ## Run 49 stress 9.520536e-05 
    ## ... Procrustes: rmse 0.002060276  max resid 0.004069095 
    ## ... Similar to previous best
    ## Run 50 stress 9.616319e-05 
    ## ... Procrustes: rmse 0.1751797  max resid 0.2172519 
    ## Run 51 stress 9.038016e-05 
    ## ... Procrustes: rmse 0.0001157209  max resid 0.0001816986 
    ## ... Similar to previous best
    ## Run 52 stress 9.408191e-05 
    ## ... Procrustes: rmse 0.0001162518  max resid 0.000190998 
    ## ... Similar to previous best
    ## Run 53 stress 0.293306 
    ## Run 54 stress 9.652548e-05 
    ## ... Procrustes: rmse 0.05667169  max resid 0.07888786 
    ## Run 55 stress 8.727126e-05 
    ## ... Procrustes: rmse 0.0001849472  max resid 0.0003172566 
    ## ... Similar to previous best
    ## Run 56 stress 8.76389e-05 
    ## ... Procrustes: rmse 0.1439178  max resid 0.2689654 
    ## Run 57 stress 9.464878e-05 
    ## ... Procrustes: rmse 0.2105332  max resid 0.3135928 
    ## Run 58 stress 9.164216e-05 
    ## ... Procrustes: rmse 0.2131919  max resid 0.3599159 
    ## Run 59 stress 8.661813e-05 
    ## ... Procrustes: rmse 0.1409328  max resid 0.2600314 
    ## Run 60 stress 9.859264e-05 
    ## ... Procrustes: rmse 0.0001993487  max resid 0.0002661956 
    ## ... Similar to previous best
    ## Run 61 stress 9.947273e-05 
    ## ... Procrustes: rmse 0.2159726  max resid 0.2578771 
    ## Run 62 stress 9.301759e-05 
    ## ... Procrustes: rmse 0.0001819054  max resid 0.0002351079 
    ## ... Similar to previous best
    ## Run 63 stress 0.2932328 
    ## Run 64 stress 9.655143e-05 
    ## ... Procrustes: rmse 0.0001975299  max resid 0.0003336778 
    ## ... Similar to previous best
    ## Run 65 stress 9.634606e-05 
    ## ... Procrustes: rmse 0.0001947051  max resid 0.000275223 
    ## ... Similar to previous best
    ## Run 66 stress 9.76166e-05 
    ## ... Procrustes: rmse 0.1536249  max resid 0.2768738 
    ## Run 67 stress 8.538774e-05 
    ## ... Procrustes: rmse 0.0002352938  max resid 0.0005045361 
    ## ... Similar to previous best
    ## Run 68 stress 9.018418e-05 
    ## ... Procrustes: rmse 0.1054294  max resid 0.1967931 
    ## Run 69 stress 9.876766e-05 
    ## ... Procrustes: rmse 0.0001222369  max resid 0.0002033895 
    ## ... Similar to previous best
    ## Run 70 stress 0.2601497 
    ## Run 71 stress 9.16079e-05 
    ## ... Procrustes: rmse 0.0001924687  max resid 0.0002605096 
    ## ... Similar to previous best
    ## Run 72 stress 0.2601482 
    ## Run 73 stress 9.939008e-05 
    ## ... Procrustes: rmse 0.0002499492  max resid 0.0005261677 
    ## ... Similar to previous best
    ## Run 74 stress 9.549476e-05 
    ## ... Procrustes: rmse 0.07586948  max resid 0.1371388 
    ## Run 75 stress 9.164697e-05 
    ## ... Procrustes: rmse 0.1589907  max resid 0.2596004 
    ## Run 76 stress 9.208367e-05 
    ## ... Procrustes: rmse 0.1593044  max resid 0.2829034 
    ## Run 77 stress 9.75241e-05 
    ## ... Procrustes: rmse 0.2132416  max resid 0.2551824 
    ## Run 78 stress 8.731949e-05 
    ## ... Procrustes: rmse 0.1968652  max resid 0.3493331 
    ## Run 79 stress 9.237559e-05 
    ## ... Procrustes: rmse 0.01160031  max resid 0.0225347 
    ## Run 80 stress 0.2601502 
    ## Run 81 stress 0.2635436 
    ## Run 82 stress 0.2570537 
    ## Run 83 stress 9.11391e-05 
    ## ... Procrustes: rmse 0.0001165859  max resid 0.0001805786 
    ## ... Similar to previous best
    ## Run 84 stress 9.627494e-05 
    ## ... Procrustes: rmse 0.1540206  max resid 0.2110361 
    ## Run 85 stress 0.2933059 
    ## Run 86 stress 9.920122e-05 
    ## ... Procrustes: rmse 0.2180708  max resid 0.259799 
    ## Run 87 stress 9.824385e-05 
    ## ... Procrustes: rmse 0.1921646  max resid 0.2348897 
    ## Run 88 stress 9.972522e-05 
    ## ... Procrustes: rmse 0.1856122  max resid 0.3001809 
    ## Run 89 stress 8.915807e-05 
    ## ... Procrustes: rmse 0.0001254104  max resid 0.0002163274 
    ## ... Similar to previous best
    ## Run 90 stress 9.88242e-05 
    ## ... Procrustes: rmse 0.1875391  max resid 0.2931933 
    ## Run 91 stress 0.2570537 
    ## Run 92 stress 9.391885e-05 
    ## ... Procrustes: rmse 0.1568995  max resid 0.2642106 
    ## Run 93 stress 9.305524e-05 
    ## ... Procrustes: rmse 0.000117775  max resid 0.0001850107 
    ## ... Similar to previous best
    ## Run 94 stress 9.283508e-05 
    ## ... Procrustes: rmse 0.01778034  max resid 0.03504135 
    ## Run 95 stress 9.242998e-05 
    ## ... Procrustes: rmse 0.0001974897  max resid 0.000256637 
    ## ... Similar to previous best
    ## Run 96 stress 9.449863e-05 
    ## ... Procrustes: rmse 0.2021093  max resid 0.2447493 
    ## Run 97 stress 9.864373e-05 
    ## ... Procrustes: rmse 0.2139311  max resid 0.3048712 
    ## Run 98 stress 9.227364e-05 
    ## ... Procrustes: rmse 0.1863738  max resid 0.2899347 
    ## Run 99 stress 9.950797e-05 
    ## ... Procrustes: rmse 0.027928  max resid 0.05400187 
    ## Run 100 stress 8.882875e-05 
    ## ... Procrustes: rmse 0.0001261478  max resid 0.0002091718 
    ## ... Similar to previous best
    ## Run 101 stress 9.475745e-05 
    ## ... Procrustes: rmse 0.1879099  max resid 0.2336143 
    ## Run 102 stress 0.2570537 
    ## Run 103 stress 0.2570537 
    ## Run 104 stress 9.871006e-05 
    ## ... Procrustes: rmse 0.0001249053  max resid 0.0002009131 
    ## ... Similar to previous best
    ## Run 105 stress 9.704307e-05 
    ## ... Procrustes: rmse 0.1875634  max resid 0.3302944 
    ## Run 106 stress 9.29121e-05 
    ## ... Procrustes: rmse 0.0001157443  max resid 0.0001882314 
    ## ... Similar to previous best
    ## Run 107 stress 9.149728e-05 
    ## ... Procrustes: rmse 0.0001135445  max resid 0.0001859132 
    ## ... Similar to previous best
    ## Run 108 stress 9.483413e-05 
    ## ... Procrustes: rmse 0.02420714  max resid 0.04799821 
    ## Run 109 stress 9.246955e-05 
    ## ... Procrustes: rmse 0.02488058  max resid 0.04922908 
    ## Run 110 stress 9.900829e-05 
    ## ... Procrustes: rmse 0.0001122317  max resid 0.000197837 
    ## ... Similar to previous best
    ## Run 111 stress 9.2297e-05 
    ## ... Procrustes: rmse 0.0001935018  max resid 0.0002707979 
    ## ... Similar to previous best
    ## Run 112 stress 9.762602e-05 
    ## ... Procrustes: rmse 0.0978587  max resid 0.177174 
    ## Run 113 stress 8.809511e-05 
    ## ... Procrustes: rmse 0.1252543  max resid 0.2337274 
    ## Run 114 stress 0.2290538 
    ## Run 115 stress 9.537204e-05 
    ## ... Procrustes: rmse 0.2080695  max resid 0.3056954 
    ## Run 116 stress 9.16845e-05 
    ## ... Procrustes: rmse 0.0001910807  max resid 0.000324853 
    ## ... Similar to previous best
    ## Run 117 stress 9.147962e-05 
    ## ... Procrustes: rmse 0.0001890996  max resid 0.0003243729 
    ## ... Similar to previous best
    ## Run 118 stress 9.551595e-05 
    ## ... Procrustes: rmse 0.0001974941  max resid 0.0002767388 
    ## ... Similar to previous best
    ## Run 119 stress 9.438295e-05 
    ## ... Procrustes: rmse 0.0001205244  max resid 0.0001917265 
    ## ... Similar to previous best
    ## Run 120 stress 8.964526e-05 
    ## ... Procrustes: rmse 0.0001974969  max resid 0.000360438 
    ## ... Similar to previous best
    ## Run 121 stress 9.51022e-05 
    ## ... Procrustes: rmse 0.0001174672  max resid 0.0001936543 
    ## ... Similar to previous best
    ## Run 122 stress 9.768673e-05 
    ## ... Procrustes: rmse 0.0002016717  max resid 0.0002820654 
    ## ... Similar to previous best
    ## Run 123 stress 9.747881e-05 
    ## ... Procrustes: rmse 0.08321459  max resid 0.1599816 
    ## Run 124 stress 9.260973e-05 
    ## ... Procrustes: rmse 0.0001156997  max resid 0.0001880616 
    ## ... Similar to previous best
    ## Run 125 stress 9.338373e-05 
    ## ... Procrustes: rmse 0.1358675  max resid 0.2356968 
    ## Run 126 stress 9.470393e-05 
    ## ... Procrustes: rmse 0.0001959071  max resid 0.0002640981 
    ## ... Similar to previous best
    ## Run 127 stress 9.982693e-05 
    ## ... Procrustes: rmse 0.1943038  max resid 0.3125055 
    ## Run 128 stress 9.559683e-05 
    ## ... Procrustes: rmse 0.1478026  max resid 0.2592272 
    ## Run 129 stress 9.570946e-05 
    ## ... Procrustes: rmse 0.1033905  max resid 0.1704092 
    ## Run 130 stress 9.568052e-05 
    ## ... Procrustes: rmse 0.2120458  max resid 0.3288521 
    ## Run 131 stress 9.296678e-05 
    ## ... Procrustes: rmse 0.1712565  max resid 0.2582863 
    ## Run 132 stress 9.787583e-05 
    ## ... Procrustes: rmse 0.1844147  max resid 0.3208488 
    ## Run 133 stress 0.2927707 
    ## Run 134 stress 9.928414e-05 
    ## ... Procrustes: rmse 0.0001238134  max resid 0.0002058136 
    ## ... Similar to previous best
    ## Run 135 stress 0.293122 
    ## Run 136 stress 9.445083e-05 
    ## ... Procrustes: rmse 0.0001877888  max resid 0.0003149239 
    ## ... Similar to previous best
    ## Run 137 stress 9.484912e-05 
    ## ... Procrustes: rmse 0.0001220441  max resid 0.0002007779 
    ## ... Similar to previous best
    ## Run 138 stress 9.089503e-05 
    ## ... Procrustes: rmse 0.0001904957  max resid 0.0003226242 
    ## ... Similar to previous best
    ## Run 139 stress 9.798774e-05 
    ## ... Procrustes: rmse 0.0001939244  max resid 0.0003205346 
    ## ... Similar to previous best
    ## Run 140 stress 7.594718e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000135415  max resid 0.000193878 
    ## ... Similar to previous best
    ## Run 141 stress 9.319093e-05 
    ## ... Procrustes: rmse 0.1498239  max resid 0.2720539 
    ## Run 142 stress 9.925921e-05 
    ## ... Procrustes: rmse 0.04778716  max resid 0.08546012 
    ## Run 143 stress 9.613953e-05 
    ## ... Procrustes: rmse 0.0001690138  max resid 0.0002200556 
    ## ... Similar to previous best
    ## Run 144 stress 9.431765e-05 
    ## ... Procrustes: rmse 0.002599501  max resid 0.00477445 
    ## ... Similar to previous best
    ## Run 145 stress 6.855747e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.07877467  max resid 0.1368944 
    ## Run 146 stress 9.300422e-05 
    ## ... Procrustes: rmse 0.07879467  max resid 0.1386575 
    ## Run 147 stress 9.56444e-05 
    ## ... Procrustes: rmse 0.05914412  max resid 0.1228678 
    ## Run 148 stress 9.441825e-05 
    ## ... Procrustes: rmse 0.04002065  max resid 0.06260509 
    ## Run 149 stress 8.684325e-05 
    ## ... Procrustes: rmse 0.07879462  max resid 0.1386589 
    ## Run 150 stress 9.043873e-05 
    ## ... Procrustes: rmse 0.0787964  max resid 0.1386777 
    ## Run 151 stress 9.581612e-05 
    ## ... Procrustes: rmse 0.1209733  max resid 0.169409 
    ## Run 152 stress 9.626516e-05 
    ## ... Procrustes: rmse 0.07879436  max resid 0.1386612 
    ## Run 153 stress 9.636113e-05 
    ## ... Procrustes: rmse 0.154598  max resid 0.1853948 
    ## Run 154 stress 0.2933059 
    ## Run 155 stress 8.565382e-05 
    ## ... Procrustes: rmse 0.07879674  max resid 0.1386795 
    ## Run 156 stress 9.234178e-05 
    ## ... Procrustes: rmse 0.1414514  max resid 0.2020286 
    ## Run 157 stress 9.112674e-05 
    ## ... Procrustes: rmse 0.1028308  max resid 0.2002008 
    ## Run 158 stress 8.848667e-05 
    ## ... Procrustes: rmse 0.07876654  max resid 0.1384692 
    ## Run 159 stress 9.491665e-05 
    ## ... Procrustes: rmse 0.0787613  max resid 0.1384541 
    ## Run 160 stress 0.2601508 
    ## Run 161 stress 7.78236e-05 
    ## ... Procrustes: rmse 0.05988773  max resid 0.08581098 
    ## Run 162 stress 9.622437e-05 
    ## ... Procrustes: rmse 0.07876283  max resid 0.1384485 
    ## Run 163 stress 9.429535e-05 
    ## ... Procrustes: rmse 0.07879501  max resid 0.138659 
    ## Run 164 stress 9.132849e-05 
    ## ... Procrustes: rmse 0.07876549  max resid 0.1384628 
    ## Run 165 stress 8.483217e-05 
    ## ... Procrustes: rmse 0.1449139  max resid 0.2007595 
    ## Run 166 stress 8.525046e-05 
    ## ... Procrustes: rmse 0.09923835  max resid 0.162796 
    ## Run 167 stress 9.563762e-05 
    ## ... Procrustes: rmse 0.1640088  max resid 0.2281528 
    ## Run 168 stress 0.2932707 
    ## Run 169 stress 0.293306 
    ## Run 170 stress 9.513782e-05 
    ## ... Procrustes: rmse 0.07879535  max resid 0.1386606 
    ## Run 171 stress 9.165324e-05 
    ## ... Procrustes: rmse 0.0733409  max resid 0.128244 
    ## Run 172 stress 9.474189e-05 
    ## ... Procrustes: rmse 0.07876385  max resid 0.138452 
    ## Run 173 stress 9.326828e-05 
    ## ... Procrustes: rmse 0.150018  max resid 0.2125325 
    ## Run 174 stress 9.695257e-05 
    ## ... Procrustes: rmse 0.149087  max resid 0.1821596 
    ## Run 175 stress 9.870699e-05 
    ## ... Procrustes: rmse 0.1225942  max resid 0.1505806 
    ## Run 176 stress 8.757476e-05 
    ## ... Procrustes: rmse 0.05123284  max resid 0.101871 
    ## Run 177 stress 9.252091e-05 
    ## ... Procrustes: rmse 0.01040805  max resid 0.02038425 
    ## Run 178 stress 9.961978e-05 
    ## ... Procrustes: rmse 0.0632833  max resid 0.08547212 
    ## Run 179 stress 9.343114e-05 
    ## ... Procrustes: rmse 0.1126396  max resid 0.1586151 
    ## Run 180 stress 9.445874e-05 
    ## ... Procrustes: rmse 0.05778012  max resid 0.0946854 
    ## Run 181 stress 0.2931216 
    ## Run 182 stress 9.927942e-05 
    ## ... Procrustes: rmse 0.04086509  max resid 0.07438748 
    ## Run 183 stress 8.940084e-05 
    ## ... Procrustes: rmse 0.07876214  max resid 0.1384649 
    ## Run 184 stress 9.302212e-05 
    ## ... Procrustes: rmse 0.07879532  max resid 0.1386603 
    ## Run 185 stress 9.087292e-05 
    ## ... Procrustes: rmse 0.01859125  max resid 0.04141153 
    ## Run 186 stress 9.441477e-05 
    ## ... Procrustes: rmse 0.1462437  max resid 0.2142954 
    ## Run 187 stress 9.693269e-05 
    ## ... Procrustes: rmse 0.05900076  max resid 0.1011606 
    ## Run 188 stress 9.696957e-05 
    ## ... Procrustes: rmse 0.01680061  max resid 0.03310962 
    ## Run 189 stress 9.229104e-05 
    ## ... Procrustes: rmse 0.07876458  max resid 0.138458 
    ## Run 190 stress 8.64467e-05 
    ## ... Procrustes: rmse 0.09422985  max resid 0.1485728 
    ## Run 191 stress 8.797009e-05 
    ## ... Procrustes: rmse 0.07876607  max resid 0.1384692 
    ## Run 192 stress 9.446791e-05 
    ## ... Procrustes: rmse 0.07741907  max resid 0.1473106 
    ## Run 193 stress 5.563959e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1317186  max resid 0.1729132 
    ## Run 194 stress 9.836984e-05 
    ## ... Procrustes: rmse 0.05138407  max resid 0.1089168 
    ## Run 195 stress 0.2933059 
    ## Run 196 stress 9.354385e-05 
    ## ... Procrustes: rmse 0.01763816  max resid 0.03367582 
    ## Run 197 stress 9.517952e-05 
    ## ... Procrustes: rmse 0.1339077  max resid 0.1993746 
    ## Run 198 stress 8.886306e-05 
    ## ... Procrustes: rmse 0.1213618  max resid 0.1652837 
    ## Run 199 stress 9.22039e-05 
    ## ... Procrustes: rmse 0.02098021  max resid 0.04311191 
    ## Run 200 stress 9.601441e-05 
    ## ... Procrustes: rmse 0.2018394  max resid 0.3150198 
    ## Run 201 stress 9.204395e-05 
    ## ... Procrustes: rmse 0.1434419  max resid 0.1948616 
    ## Run 202 stress 0.2570537 
    ## Run 203 stress 9.785588e-05 
    ## ... Procrustes: rmse 0.1970134  max resid 0.3080924 
    ## Run 204 stress 9.27389e-05 
    ## ... Procrustes: rmse 0.07080243  max resid 0.1435815 
    ## Run 205 stress 8.574965e-05 
    ## ... Procrustes: rmse 0.2005423  max resid 0.3131947 
    ## Run 206 stress 0.2601495 
    ## Run 207 stress 9.269636e-05 
    ## ... Procrustes: rmse 0.05256885  max resid 0.06776522 
    ## Run 208 stress 8.902555e-05 
    ## ... Procrustes: rmse 0.1331838  max resid 0.2206572 
    ## Run 209 stress 9.932691e-05 
    ## ... Procrustes: rmse 0.06949973  max resid 0.1325094 
    ## Run 210 stress 9.196322e-05 
    ## ... Procrustes: rmse 0.03179719  max resid 0.06995227 
    ## Run 211 stress 9.943006e-05 
    ## ... Procrustes: rmse 0.02901341  max resid 0.05672161 
    ## Run 212 stress 9.341193e-05 
    ## ... Procrustes: rmse 0.03341861  max resid 0.07466064 
    ## Run 213 stress 9.705727e-05 
    ## ... Procrustes: rmse 0.04621895  max resid 0.08768673 
    ## Run 214 stress 9.502394e-05 
    ## ... Procrustes: rmse 0.2018393  max resid 0.3150573 
    ## Run 215 stress 9.254645e-05 
    ## ... Procrustes: rmse 0.200752  max resid 0.312833 
    ## Run 216 stress 9.076409e-05 
    ## ... Procrustes: rmse 0.03768944  max resid 0.08393808 
    ## Run 217 stress 9.773218e-05 
    ## ... Procrustes: rmse 0.1150762  max resid 0.1819284 
    ## Run 218 stress 8.775194e-05 
    ## ... Procrustes: rmse 0.03774219  max resid 0.08290541 
    ## Run 219 stress 0.2598985 
    ## Run 220 stress 9.571464e-05 
    ## ... Procrustes: rmse 0.2018372  max resid 0.3150177 
    ## Run 221 stress 0.2927708 
    ## Run 222 stress 9.574475e-05 
    ## ... Procrustes: rmse 0.2006677  max resid 0.3126485 
    ## Run 223 stress 7.977609e-05 
    ## ... Procrustes: rmse 0.2009374  max resid 0.3138106 
    ## Run 224 stress 9.501919e-05 
    ## ... Procrustes: rmse 0.1620186  max resid 0.2346657 
    ## Run 225 stress 9.168693e-05 
    ## ... Procrustes: rmse 0.2018393  max resid 0.3150209 
    ## Run 226 stress 9.154924e-05 
    ## ... Procrustes: rmse 0.05298439  max resid 0.1096142 
    ## Run 227 stress 9.693227e-05 
    ## ... Procrustes: rmse 0.2018236  max resid 0.3148407 
    ## Run 228 stress 0.2294239 
    ## Run 229 stress 9.46469e-05 
    ## ... Procrustes: rmse 0.20182  max resid 0.3148446 
    ## Run 230 stress 9.865176e-05 
    ## ... Procrustes: rmse 0.2018377  max resid 0.3150168 
    ## Run 231 stress 9.799335e-05 
    ## ... Procrustes: rmse 0.2018227  max resid 0.3148383 
    ## Run 232 stress 9.324448e-05 
    ## ... Procrustes: rmse 0.2018235  max resid 0.3148489 
    ## Run 233 stress 9.852338e-05 
    ## ... Procrustes: rmse 0.2018397  max resid 0.3150313 
    ## Run 234 stress 9.021323e-05 
    ## ... Procrustes: rmse 0.07719405  max resid 0.1434362 
    ## Run 235 stress 8.968108e-05 
    ## ... Procrustes: rmse 0.06337435  max resid 0.1033353 
    ## Run 236 stress 8.9137e-05 
    ## ... Procrustes: rmse 0.08881718  max resid 0.1472428 
    ## Run 237 stress 9.829895e-05 
    ## ... Procrustes: rmse 0.008034389  max resid 0.01528731 
    ## Run 238 stress 8.759109e-05 
    ## ... Procrustes: rmse 0.2018243  max resid 0.3148602 
    ## Run 239 stress 0.2601509 
    ## Run 240 stress 8.252565e-05 
    ## ... Procrustes: rmse 0.04177099  max resid 0.09145816 
    ## Run 241 stress 9.042232e-05 
    ## ... Procrustes: rmse 0.2018398  max resid 0.3150232 
    ## Run 242 stress 9.217559e-05 
    ## ... Procrustes: rmse 0.2018188  max resid 0.3148474 
    ## Run 243 stress 9.214302e-05 
    ## ... Procrustes: rmse 0.06135217  max resid 0.1119444 
    ## Run 244 stress 9.319008e-05 
    ## ... Procrustes: rmse 0.2018239  max resid 0.3148497 
    ## Run 245 stress 0.2570537 
    ## Run 246 stress 8.972671e-05 
    ## ... Procrustes: rmse 0.2018392  max resid 0.3150214 
    ## Run 247 stress 8.871112e-05 
    ## ... Procrustes: rmse 0.2018393  max resid 0.3150314 
    ## Run 248 stress 9.585201e-05 
    ## ... Procrustes: rmse 0.03782862  max resid 0.08420884 
    ## Run 249 stress 9.943809e-05 
    ## ... Procrustes: rmse 0.2018227  max resid 0.3148354 
    ## Run 250 stress 8.204944e-05 
    ## ... Procrustes: rmse 0.01962945  max resid 0.0297036 
    ## Run 251 stress 0.2933061 
    ## Run 252 stress 9.819991e-05 
    ## ... Procrustes: rmse 0.03919451  max resid 0.08513646 
    ## Run 253 stress 8.873671e-05 
    ## ... Procrustes: rmse 0.0431255  max resid 0.05671025 
    ## Run 254 stress 9.897157e-05 
    ## ... Procrustes: rmse 0.03004666  max resid 0.05825016 
    ## Run 255 stress 9.68618e-05 
    ## ... Procrustes: rmse 0.1140315  max resid 0.1790949 
    ## Run 256 stress 8.913916e-05 
    ## ... Procrustes: rmse 0.2018395  max resid 0.3150224 
    ## Run 257 stress 9.2335e-05 
    ## ... Procrustes: rmse 0.05678909  max resid 0.08387784 
    ## Run 258 stress 8.98048e-05 
    ## ... Procrustes: rmse 0.06722845  max resid 0.1411927 
    ## Run 259 stress 8.9049e-05 
    ## ... Procrustes: rmse 0.2018247  max resid 0.3148607 
    ## Run 260 stress 9.754117e-05 
    ## ... Procrustes: rmse 0.05190332  max resid 0.07180483 
    ## Run 261 stress 9.652994e-05 
    ## ... Procrustes: rmse 0.0286201  max resid 0.06275234 
    ## Run 262 stress 9.439345e-05 
    ## ... Procrustes: rmse 0.2018391  max resid 0.3150267 
    ## Run 263 stress 8.586239e-05 
    ## ... Procrustes: rmse 0.2018284  max resid 0.31486 
    ## Run 264 stress 9.060724e-05 
    ## ... Procrustes: rmse 0.2018238  max resid 0.3148543 
    ## Run 265 stress 9.49741e-05 
    ## ... Procrustes: rmse 0.1696366  max resid 0.2478776 
    ## Run 266 stress 9.891822e-05 
    ## ... Procrustes: rmse 0.03859482  max resid 0.08568006 
    ## Run 267 stress 9.379896e-05 
    ## ... Procrustes: rmse 0.09724663  max resid 0.151834 
    ## Run 268 stress 8.494632e-05 
    ## ... Procrustes: rmse 0.03727397  max resid 0.08377059 
    ## Run 269 stress 9.566314e-05 
    ## ... Procrustes: rmse 0.04083124  max resid 0.07372386 
    ## Run 270 stress 9.803518e-05 
    ## ... Procrustes: rmse 0.03030186  max resid 0.06763347 
    ## Run 271 stress 9.698082e-05 
    ## ... Procrustes: rmse 0.2018398  max resid 0.3150203 
    ## Run 272 stress 7.526523e-05 
    ## ... Procrustes: rmse 0.07937412  max resid 0.1323139 
    ## Run 273 stress 0.2601532 
    ## Run 274 stress 0.2601669 
    ## Run 275 stress 9.176538e-05 
    ## ... Procrustes: rmse 0.01988298  max resid 0.0282101 
    ## Run 276 stress 0.2865417 
    ## Run 277 stress 9.167176e-05 
    ## ... Procrustes: rmse 0.1875132  max resid 0.2853463 
    ## Run 278 stress 9.787721e-05 
    ## ... Procrustes: rmse 0.02973804  max resid 0.06240995 
    ## Run 279 stress 9.532366e-05 
    ## ... Procrustes: rmse 0.00872657  max resid 0.01746695 
    ## Run 280 stress 0.2831454 
    ## Run 281 stress 9.996594e-05 
    ## ... Procrustes: rmse 0.03422382  max resid 0.07595635 
    ## Run 282 stress 9.496075e-05 
    ## ... Procrustes: rmse 0.2018232  max resid 0.3148458 
    ## Run 283 stress 0.2933059 
    ## Run 284 stress 9.647152e-05 
    ## ... Procrustes: rmse 0.2018391  max resid 0.315026 
    ## Run 285 stress 0.2294232 
    ## Run 286 stress 0.2927399 
    ## Run 287 stress 8.839324e-05 
    ## ... Procrustes: rmse 0.2018246  max resid 0.3148611 
    ## Run 288 stress 9.825484e-05 
    ## ... Procrustes: rmse 0.1902232  max resid 0.291678 
    ## Run 289 stress 9.025452e-05 
    ## ... Procrustes: rmse 0.2018208  max resid 0.314852 
    ## Run 290 stress 9.751169e-05 
    ## ... Procrustes: rmse 0.121354  max resid 0.1753702 
    ## Run 291 stress 8.681584e-05 
    ## ... Procrustes: rmse 0.05812228  max resid 0.1124108 
    ## Run 292 stress 9.569632e-05 
    ## ... Procrustes: rmse 0.2018243  max resid 0.3148464 
    ## Run 293 stress 0.2601495 
    ## Run 294 stress 9.781842e-05 
    ## ... Procrustes: rmse 0.04807769  max resid 0.0933465 
    ## Run 295 stress 9.676522e-05 
    ## ... Procrustes: rmse 0.06622912  max resid 0.1378508 
    ## Run 296 stress 9.904482e-05 
    ## ... Procrustes: rmse 0.2018196  max resid 0.3148364 
    ## Run 297 stress 9.93497e-05 
    ## ... Procrustes: rmse 0.1860865  max resid 0.2835898 
    ## Run 298 stress 9.466788e-05 
    ## ... Procrustes: rmse 0.03992181  max resid 0.05266497 
    ## Run 299 stress 9.875585e-05 
    ## ... Procrustes: rmse 0.03093345  max resid 0.06728026 
    ## Run 300 stress 8.880099e-05 
    ## ... Procrustes: rmse 0.1836914  max resid 0.2783589 
    ## Run 301 stress 8.678893e-05 
    ## ... Procrustes: rmse 0.1472754  max resid 0.2415477 
    ## Run 302 stress 8.979305e-05 
    ## ... Procrustes: rmse 0.08966216  max resid 0.126591 
    ## Run 303 stress 9.881688e-05 
    ## ... Procrustes: rmse 0.03074643  max resid 0.06683811 
    ## Run 304 stress 5.413186e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03406447  max resid 0.07622744 
    ## Run 305 stress 0.2598893 
    ## Run 306 stress 9.428242e-05 
    ## ... Procrustes: rmse 0.0692149  max resid 0.1562561 
    ## Run 307 stress 9.782024e-05 
    ## ... Procrustes: rmse 0.1848214  max resid 0.2447301 
    ## Run 308 stress 9.219139e-05 
    ## ... Procrustes: rmse 0.05835775  max resid 0.1258154 
    ## Run 309 stress 9.32991e-05 
    ## ... Procrustes: rmse 0.1794632  max resid 0.24197 
    ## Run 310 stress 8.341454e-05 
    ## ... Procrustes: rmse 0.03314613  max resid 0.07175862 
    ## Run 311 stress 9.887068e-05 
    ## ... Procrustes: rmse 0.1506572  max resid 0.2257341 
    ## Run 312 stress 9.233139e-05 
    ## ... Procrustes: rmse 0.003479641  max resid 0.004441 
    ## ... Similar to previous best
    ## Run 313 stress 9.257096e-05 
    ## ... Procrustes: rmse 0.1848239  max resid 0.2447014 
    ## Run 314 stress 8.748904e-05 
    ## ... Procrustes: rmse 0.06238753  max resid 0.1419543 
    ## Run 315 stress 0.2601497 
    ## Run 316 stress 9.488877e-05 
    ## ... Procrustes: rmse 0.04397519  max resid 0.09744305 
    ## Run 317 stress 0.2933059 
    ## Run 318 stress 0.2294233 
    ## Run 319 stress 9.729562e-05 
    ## ... Procrustes: rmse 0.1457368  max resid 0.2242916 
    ## Run 320 stress 9.037545e-05 
    ## ... Procrustes: rmse 0.1848247  max resid 0.2447125 
    ## Run 321 stress 8.716302e-05 
    ## ... Procrustes: rmse 0.1848239  max resid 0.2447259 
    ## Run 322 stress 0.293122 
    ## Run 323 stress 9.561164e-05 
    ## ... Procrustes: rmse 0.1424392  max resid 0.2223132 
    ## Run 324 stress 0.293122 
    ## Run 325 stress 9.998341e-05 
    ## ... Procrustes: rmse 0.1558781  max resid 0.2059191 
    ## Run 326 stress 8.875895e-05 
    ## ... Procrustes: rmse 0.1770571  max resid 0.2399533 
    ## Run 327 stress 8.722415e-05 
    ## ... Procrustes: rmse 0.08509535  max resid 0.1796277 
    ## Run 328 stress 8.423532e-05 
    ## ... Procrustes: rmse 0.04022051  max resid 0.04801834 
    ## Run 329 stress 9.599167e-05 
    ## ... Procrustes: rmse 0.007720732  max resid 0.01719309 
    ## Run 330 stress 9.690439e-05 
    ## ... Procrustes: rmse 0.04241485  max resid 0.08466037 
    ## Run 331 stress 8.906247e-05 
    ## ... Procrustes: rmse 0.1844615  max resid 0.2443359 
    ## Run 332 stress 9.548588e-05 
    ## ... Procrustes: rmse 0.01695047  max resid 0.02074301 
    ## Run 333 stress 9.301358e-05 
    ## ... Procrustes: rmse 0.05916419  max resid 0.1255875 
    ## Run 334 stress 9.25081e-05 
    ## ... Procrustes: rmse 0.03954515  max resid 0.08908011 
    ## Run 335 stress 9.040159e-05 
    ## ... Procrustes: rmse 0.1848248  max resid 0.2447228 
    ## Run 336 stress 9.382543e-05 
    ## ... Procrustes: rmse 0.06388377  max resid 0.1447593 
    ## Run 337 stress 8.791809e-05 
    ## ... Procrustes: rmse 0.05480448  max resid 0.06451491 
    ## Run 338 stress 9.174289e-05 
    ## ... Procrustes: rmse 0.1834466  max resid 0.24361 
    ## Run 339 stress 9.56608e-05 
    ## ... Procrustes: rmse 0.04898519  max resid 0.05800688 
    ## Run 340 stress 9.416356e-05 
    ## ... Procrustes: rmse 0.04585804  max resid 0.1026337 
    ## Run 341 stress 9.216196e-05 
    ## ... Procrustes: rmse 0.05518935  max resid 0.06494938 
    ## Run 342 stress 9.152709e-05 
    ## ... Procrustes: rmse 0.04534864  max resid 0.09735244 
    ## Run 343 stress 9.100729e-05 
    ## ... Procrustes: rmse 0.1848244  max resid 0.2447345 
    ## Run 344 stress 0.2294234 
    ## Run 345 stress 9.653073e-05 
    ## ... Procrustes: rmse 0.1619524  max resid 0.2322259 
    ## Run 346 stress 0.293122 
    ## Run 347 stress 9.919607e-05 
    ## ... Procrustes: rmse 0.1068187  max resid 0.1395711 
    ## Run 348 stress 9.720785e-05 
    ## ... Procrustes: rmse 0.03782867  max resid 0.04528505 
    ## Run 349 stress 9.606758e-05 
    ## ... Procrustes: rmse 0.08380786  max resid 0.1766787 
    ## Run 350 stress 9.151608e-05 
    ## ... Procrustes: rmse 0.01374289  max resid 0.0168978 
    ## Run 351 stress 9.569948e-05 
    ## ... Procrustes: rmse 0.1848231  max resid 0.2447033 
    ## Run 352 stress 9.475678e-05 
    ## ... Procrustes: rmse 0.1848225  max resid 0.2447277 
    ## Run 353 stress 9.177389e-05 
    ## ... Procrustes: rmse 0.0417702  max resid 0.08005078 
    ## Run 354 stress 9.676489e-05 
    ## ... Procrustes: rmse 0.02023976  max resid 0.03847653 
    ## Run 355 stress 9.726931e-05 
    ## ... Procrustes: rmse 0.05295236  max resid 0.1153679 
    ## Run 356 stress 9.875865e-05 
    ## ... Procrustes: rmse 0.1067157  max resid 0.154351 
    ## Run 357 stress 9.3249e-05 
    ## ... Procrustes: rmse 0.1848222  max resid 0.2447214 
    ## Run 358 stress 9.744829e-05 
    ## ... Procrustes: rmse 0.04211348  max resid 0.08952181 
    ## Run 359 stress 0.2570537 
    ## Run 360 stress 8.204311e-05 
    ## ... Procrustes: rmse 0.0353642  max resid 0.08014597 
    ## Run 361 stress 9.568418e-05 
    ## ... Procrustes: rmse 0.03284922  max resid 0.071887 
    ## Run 362 stress 9.140964e-05 
    ## ... Procrustes: rmse 0.1804082  max resid 0.2392819 
    ## Run 363 stress 9.61248e-05 
    ## ... Procrustes: rmse 0.0532337  max resid 0.1044997 
    ## Run 364 stress 9.716724e-05 
    ## ... Procrustes: rmse 0.1794244  max resid 0.242052 
    ## Run 365 stress 9.89038e-05 
    ## ... Procrustes: rmse 0.02854913  max resid 0.06344484 
    ## Run 366 stress 9.583013e-05 
    ## ... Procrustes: rmse 0.1848237  max resid 0.2447448 
    ## Run 367 stress 9.071231e-05 
    ## ... Procrustes: rmse 0.09072754  max resid 0.1955195 
    ## Run 368 stress 9.231949e-05 
    ## ... Procrustes: rmse 0.008282369  max resid 0.01029546 
    ## Run 369 stress 9.815891e-05 
    ## ... Procrustes: rmse 0.03851981  max resid 0.07839088 
    ## Run 370 stress 8.062333e-05 
    ## ... Procrustes: rmse 0.04506811  max resid 0.102214 
    ## Run 371 stress 9.253705e-05 
    ## ... Procrustes: rmse 0.08651251  max resid 0.1783622 
    ## Run 372 stress 9.412655e-05 
    ## ... Procrustes: rmse 0.184824  max resid 0.2447413 
    ## Run 373 stress 9.539387e-05 
    ## ... Procrustes: rmse 0.1848227  max resid 0.2447274 
    ## Run 374 stress 8.082794e-05 
    ## ... Procrustes: rmse 0.1846135  max resid 0.2444023 
    ## Run 375 stress 9.336522e-05 
    ## ... Procrustes: rmse 0.01841434  max resid 0.03444402 
    ## Run 376 stress 9.085125e-05 
    ## ... Procrustes: rmse 0.05394665  max resid 0.1212553 
    ## Run 377 stress 9.898327e-05 
    ## ... Procrustes: rmse 0.04903289  max resid 0.05794536 
    ## Run 378 stress 9.27144e-05 
    ## ... Procrustes: rmse 0.1848235  max resid 0.2447318 
    ## Run 379 stress 9.006062e-05 
    ## ... Procrustes: rmse 0.1009587  max resid 0.1976236 
    ## Run 380 stress 8.570596e-05 
    ## ... Procrustes: rmse 0.04622496  max resid 0.09939867 
    ## Run 381 stress 9.316505e-05 
    ## ... Procrustes: rmse 0.03346513  max resid 0.04225138 
    ## Run 382 stress 9.981e-05 
    ## ... Procrustes: rmse 0.1848227  max resid 0.2447407 
    ## Run 383 stress 8.789792e-05 
    ## ... Procrustes: rmse 0.04175143  max resid 0.09022483 
    ## Run 384 stress 9.057341e-05 
    ## ... Procrustes: rmse 0.1827158  max resid 0.2435024 
    ## Run 385 stress 9.578422e-05 
    ## ... Procrustes: rmse 0.08239133  max resid 0.1740318 
    ## Run 386 stress 9.114504e-05 
    ## ... Procrustes: rmse 0.1848233  max resid 0.244722 
    ## Run 387 stress 8.200233e-05 
    ## ... Procrustes: rmse 0.03459463  max resid 0.07613134 
    ## Run 388 stress 9.523741e-05 
    ## ... Procrustes: rmse 0.1847327  max resid 0.2446004 
    ## Run 389 stress 9.299435e-05 
    ## ... Procrustes: rmse 0.1036182  max resid 0.1872708 
    ## Run 390 stress 9.411772e-05 
    ## ... Procrustes: rmse 0.05582473  max resid 0.110046 
    ## Run 391 stress 9.467164e-05 
    ## ... Procrustes: rmse 0.1848238  max resid 0.2447053 
    ## Run 392 stress 9.443313e-05 
    ## ... Procrustes: rmse 0.06334987  max resid 0.07383606 
    ## Run 393 stress 8.88667e-05 
    ## ... Procrustes: rmse 0.1768363  max resid 0.2396149 
    ## Run 394 stress 9.839197e-05 
    ## ... Procrustes: rmse 0.109968  max resid 0.2002874 
    ## Run 395 stress 0.2601669 
    ## Run 396 stress 9.6043e-05 
    ## ... Procrustes: rmse 0.05007714  max resid 0.06017059 
    ## Run 397 stress 9.812052e-05 
    ## ... Procrustes: rmse 0.1848237  max resid 0.2447475 
    ## Run 398 stress 0.293306 
    ## Run 399 stress 9.327188e-05 
    ## ... Procrustes: rmse 0.02104132  max resid 0.02956669 
    ## Run 400 stress 0.2601507 
    ## Run 401 stress 9.550793e-05 
    ## ... Procrustes: rmse 0.1628835  max resid 0.2179171 
    ## Run 402 stress 9.134868e-05 
    ## ... Procrustes: rmse 0.07391039  max resid 0.1594814 
    ## Run 403 stress 9.191109e-05 
    ## ... Procrustes: rmse 0.05769368  max resid 0.1290605 
    ## Run 404 stress 9.329584e-05 
    ## ... Procrustes: rmse 0.1810853  max resid 0.242717 
    ## Run 405 stress 8.801744e-05 
    ## ... Procrustes: rmse 0.0151665  max resid 0.01860117 
    ## Run 406 stress 9.569596e-05 
    ## ... Procrustes: rmse 0.09702584  max resid 0.1343956 
    ## Run 407 stress 0.2570537 
    ## Run 408 stress 9.856967e-05 
    ## ... Procrustes: rmse 0.1535126  max resid 0.2264221 
    ## Run 409 stress 0.2812919 
    ## Run 410 stress 9.276311e-05 
    ## ... Procrustes: rmse 0.1848235  max resid 0.2447357 
    ## Run 411 stress 9.40407e-05 
    ## ... Procrustes: rmse 0.1489854  max resid 0.2258519 
    ## Run 412 stress 9.227044e-05 
    ## ... Procrustes: rmse 0.01331664  max resid 0.0269354 
    ## Run 413 stress 0.293306 
    ## Run 414 stress 9.562027e-05 
    ## ... Procrustes: rmse 0.1848218  max resid 0.2447254 
    ## Run 415 stress 9.905431e-05 
    ## ... Procrustes: rmse 0.134565  max resid 0.2173536 
    ## Run 416 stress 9.43909e-05 
    ## ... Procrustes: rmse 0.1377155  max resid 0.2066642 
    ## Run 417 stress 9.4576e-05 
    ## ... Procrustes: rmse 0.04834348  max resid 0.1083701 
    ## Run 418 stress 9.78667e-05 
    ## ... Procrustes: rmse 0.06736689  max resid 0.1415012 
    ## Run 419 stress 9.342297e-05 
    ## ... Procrustes: rmse 0.1847006  max resid 0.2445045 
    ## Run 420 stress 8.806008e-05 
    ## ... Procrustes: rmse 0.1831053  max resid 0.2434226 
    ## Run 421 stress 9.87256e-05 
    ## ... Procrustes: rmse 0.1848227  max resid 0.2447373 
    ## Run 422 stress 0.293122 
    ## Run 423 stress 9.189721e-05 
    ## ... Procrustes: rmse 0.1464338  max resid 0.2167946 
    ## Run 424 stress 9.323621e-05 
    ## ... Procrustes: rmse 0.02837859  max resid 0.03427429 
    ## Run 425 stress 9.489819e-05 
    ## ... Procrustes: rmse 0.1848234  max resid 0.244736 
    ## Run 426 stress 9.71539e-05 
    ## ... Procrustes: rmse 0.0396428  max resid 0.04734616 
    ## Run 427 stress 9.768403e-05 
    ## ... Procrustes: rmse 0.05475239  max resid 0.1093545 
    ## Run 428 stress 9.019908e-05 
    ## ... Procrustes: rmse 0.05940521  max resid 0.1245026 
    ## Run 429 stress 9.740394e-05 
    ## ... Procrustes: rmse 0.06329613  max resid 0.1351989 
    ## Run 430 stress 9.589478e-05 
    ## ... Procrustes: rmse 0.02941599  max resid 0.06095244 
    ## Run 431 stress 8.457989e-05 
    ## ... Procrustes: rmse 0.02673119  max resid 0.0323488 
    ## Run 432 stress 9.071495e-05 
    ## ... Procrustes: rmse 0.07549927  max resid 0.1725207 
    ## Run 433 stress 9.395767e-05 
    ## ... Procrustes: rmse 0.1848239  max resid 0.2447046 
    ## Run 434 stress 8.630725e-05 
    ## ... Procrustes: rmse 0.1848212  max resid 0.2447067 
    ## Run 435 stress 9.247338e-05 
    ## ... Procrustes: rmse 0.06307516  max resid 0.1417 
    ## Run 436 stress 7.595403e-05 
    ## ... Procrustes: rmse 0.1848246  max resid 0.244704 
    ## Run 437 stress 0.2601501 
    ## Run 438 stress 9.146755e-05 
    ## ... Procrustes: rmse 0.03931497  max resid 0.08836444 
    ## Run 439 stress 9.918929e-05 
    ## ... Procrustes: rmse 0.1785546  max resid 0.2409051 
    ## Run 440 stress 9.364329e-05 
    ## ... Procrustes: rmse 0.05929072  max resid 0.1335107 
    ## Run 441 stress 8.679277e-05 
    ## ... Procrustes: rmse 0.1321403  max resid 0.1949112 
    ## Run 442 stress 8.910656e-05 
    ## ... Procrustes: rmse 0.1829954  max resid 0.2436312 
    ## Run 443 stress 9.827126e-05 
    ## ... Procrustes: rmse 0.05802928  max resid 0.111146 
    ## Run 444 stress 9.603636e-05 
    ## ... Procrustes: rmse 0.1424232  max resid 0.2180606 
    ## Run 445 stress 9.108259e-05 
    ## ... Procrustes: rmse 0.0243759  max resid 0.04657388 
    ## Run 446 stress 9.697855e-05 
    ## ... Procrustes: rmse 0.006383614  max resid 0.007992179 
    ## Run 447 stress 9.486215e-05 
    ## ... Procrustes: rmse 0.1848238  max resid 0.2447067 
    ## Run 448 stress 9.394404e-05 
    ## ... Procrustes: rmse 0.1848238  max resid 0.244704 
    ## Run 449 stress 9.494745e-05 
    ## ... Procrustes: rmse 0.1847294  max resid 0.2445798 
    ## Run 450 stress 0.2570537 
    ## Run 451 stress 9.587672e-05 
    ## ... Procrustes: rmse 0.1848214  max resid 0.2447274 
    ## Run 452 stress 8.972849e-05 
    ## ... Procrustes: rmse 0.1848225  max resid 0.2447167 
    ## Run 453 stress 9.551388e-05 
    ## ... Procrustes: rmse 0.07348898  max resid 0.154007 
    ## Run 454 stress 8.946718e-05 
    ## ... Procrustes: rmse 0.1847293  max resid 0.2445737 
    ## Run 455 stress 9.027794e-05 
    ## ... Procrustes: rmse 0.06754208  max resid 0.1430204 
    ## Run 456 stress 0.2570537 
    ## Run 457 stress 9.547399e-05 
    ## ... Procrustes: rmse 0.1012499  max resid 0.1576614 
    ## Run 458 stress 8.585509e-05 
    ## ... Procrustes: rmse 0.08106053  max resid 0.1787451 
    ## Run 459 stress 9.056416e-05 
    ## ... Procrustes: rmse 0.03458891  max resid 0.073187 
    ## Run 460 stress 0.2930582 
    ## Run 461 stress 9.95954e-05 
    ## ... Procrustes: rmse 0.06270308  max resid 0.09859369 
    ## Run 462 stress 9.561619e-05 
    ## ... Procrustes: rmse 0.08415063  max resid 0.1456497 
    ## Run 463 stress 9.777322e-05 
    ## ... Procrustes: rmse 0.1848234  max resid 0.2447458 
    ## Run 464 stress 9.976885e-05 
    ## ... Procrustes: rmse 0.1231334  max resid 0.215762 
    ## Run 465 stress 9.076169e-05 
    ## ... Procrustes: rmse 0.1819838  max resid 0.2411982 
    ## Run 466 stress 9.610317e-05 
    ## ... Procrustes: rmse 0.1459748  max resid 0.2250022 
    ## Run 467 stress 9.549557e-05 
    ## ... Procrustes: rmse 0.07732238  max resid 0.1225534 
    ## Run 468 stress 0.2570537 
    ## Run 469 stress 9.44219e-05 
    ## ... Procrustes: rmse 0.1049742  max resid 0.2033131 
    ## Run 470 stress 9.661267e-05 
    ## ... Procrustes: rmse 0.1848228  max resid 0.2447048 
    ## Run 471 stress 9.971701e-05 
    ## ... Procrustes: rmse 0.1656589  max resid 0.2356954 
    ## Run 472 stress 8.728978e-05 
    ## ... Procrustes: rmse 0.1848238  max resid 0.2447285 
    ## Run 473 stress 9.095284e-05 
    ## ... Procrustes: rmse 0.1037343  max resid 0.2006917 
    ## Run 474 stress 9.817655e-05 
    ## ... Procrustes: rmse 0.1846906  max resid 0.244521 
    ## Run 475 stress 9.28238e-05 
    ## ... Procrustes: rmse 0.02774687  max resid 0.03633525 
    ## Run 476 stress 9.960605e-05 
    ## ... Procrustes: rmse 0.1116857  max resid 0.1872683 
    ## Run 477 stress 8.728772e-05 
    ## ... Procrustes: rmse 0.1848237  max resid 0.2446878 
    ## Run 478 stress 0.2933059 
    ## Run 479 stress 9.773058e-05 
    ## ... Procrustes: rmse 0.1848231  max resid 0.2447077 
    ## Run 480 stress 9.084169e-05 
    ## ... Procrustes: rmse 0.1256733  max resid 0.166247 
    ## Run 481 stress 9.285362e-05 
    ## ... Procrustes: rmse 0.05937791  max resid 0.08495967 
    ## Run 482 stress 0.260149 
    ## Run 483 stress 9.077611e-05 
    ## ... Procrustes: rmse 0.02862461  max resid 0.04004718 
    ## Run 484 stress 9.496711e-05 
    ## ... Procrustes: rmse 0.04342879  max resid 0.06939688 
    ## Run 485 stress 0.2290537 
    ## Run 486 stress 9.498826e-05 
    ## ... Procrustes: rmse 0.1467311  max resid 0.2230551 
    ## Run 487 stress 8.708581e-05 
    ## ... Procrustes: rmse 0.1848239  max resid 0.2446885 
    ## Run 488 stress 9.054992e-05 
    ## ... Procrustes: rmse 0.1848221  max resid 0.2447162 
    ## Run 489 stress 9.252885e-05 
    ## ... Procrustes: rmse 0.05925095  max resid 0.1342146 
    ## Run 490 stress 9.490515e-05 
    ## ... Procrustes: rmse 0.1848235  max resid 0.2447026 
    ## Run 491 stress 0.2601503 
    ## Run 492 stress 9.559127e-05 
    ## ... Procrustes: rmse 0.1848225  max resid 0.2447272 
    ## Run 493 stress 9.911865e-05 
    ## ... Procrustes: rmse 0.06734886  max resid 0.1428602 
    ## Run 494 stress 9.595664e-05 
    ## ... Procrustes: rmse 0.1848228  max resid 0.2447317 
    ## Run 495 stress 9.82986e-05 
    ## ... Procrustes: rmse 0.1457584  max resid 0.2255464 
    ## Run 496 stress 9.340235e-05 
    ## ... Procrustes: rmse 0.1848238  max resid 0.2447276 
    ## Run 497 stress 9.804689e-05 
    ## ... Procrustes: rmse 0.1848233  max resid 0.2447422 
    ## Run 498 stress 0.293122 
    ## Run 499 stress 8.455115e-05 
    ## ... Procrustes: rmse 0.184823  max resid 0.2447181 
    ## Run 500 stress 9.130749e-05 
    ## ... Procrustes: rmse 0.1848213  max resid 0.2447164 
    ## *** Best solution repeated 1 times

    ## Warning in metaMDS(PD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, :
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

### PD beta env and geo correlated variables using envfit

``` r
### Environmental
# Surveyed sites 
# For figure
(PD_beta_env_ef <- envfit(PD_beta_env_NMDS, env[surveyed_sites_env,c(34)], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                                  NMDS1   NMDS2    r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31834 0.699  0.016 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_env_efp <- p.adjust.envfit(PD_beta_env_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                                  NMDS1   NMDS2    r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31834 0.699  0.016 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by stratification
(PD_beta_env_A_ef <- envfit(PD_beta_env_NMDS, env[surveyed_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.75816 -0.65207 0.0724  0.680  
    ## salinity_median     0.94798  0.31834 0.6990  0.012 *
    ## oxygen_median       0.57184 -0.82036 0.6836  0.110  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_env_A_efp <- p.adjust.envfit(PD_beta_env_A_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.75816 -0.65207 0.0724  1.000  
    ## salinity_median     0.94798  0.31834 0.6990  0.036 *
    ## oxygen_median       0.57184 -0.82036 0.6836  0.330  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
(PD_beta_env_MS_ef <- envfit(PD_beta_env_MS_NMDS, env[mixed_stratified_lakes,environment], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.92677  0.37564 0.0954  0.751  
    ## salinity_median     0.90474  0.42596 0.6840  0.029 *
    ## oxygen_median       0.44478 -0.89564 0.5837  0.209  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_env_MS_efp <- p.adjust.envfit(PD_beta_env_MS_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.92677  0.37564 0.0954  1.000  
    ## salinity_median     0.90474  0.42596 0.6840  0.087 .
    ## oxygen_median       0.44478 -0.89564 0.5837  0.627  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
(PD_beta_env_OM_ef <- envfit(PD_beta_env_OM_NMDS, env[ocean_mixed_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.52519 -0.85098 0.0223  0.905  
    ## salinity_median    -0.79785 -0.60286 0.7269  0.016 *
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.019 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_env_OM_efp <- p.adjust.envfit(PD_beta_env_OM_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.52519 -0.85098 0.0223  1.000  
    ## salinity_median    -0.79785 -0.60286 0.7269  0.048 *
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.057 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
(PD_beta_env_SO_ef <- envfit(PD_beta_env_SO_NMDS, env[ocean_stratified_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median  0.0007130  1.0000000 0.0512  0.656
    ## salinity_median    -0.0021384 -1.0000000 0.5220  0.378
    ## oxygen_median      -0.0056266  0.9999800 0.7192  0.642
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_env_SO_efp <- p.adjust.envfit(PD_beta_env_SO_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median  0.0007130  1.0000000 0.0512      1
    ## salinity_median    -0.0021384 -1.0000000 0.5220      1
    ## oxygen_median      -0.0056266  0.9999800 0.7192      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed lakes
(PD_beta_env_M_ef <- envfit(PD_beta_env_M_NMDS, env[mixed_lakes,environment], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)  
    ## temperature_median -0.00007896 -1.00000000 0.2066  0.570  
    ## salinity_median     0.00053353  1.00000000 0.4558  0.224  
    ## oxygen_median       0.00110373  1.00000000 0.6086  0.097 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_env_M_efp <- p.adjust.envfit(PD_beta_env_M_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median -0.00007896 -1.00000000 0.2066  1.000
    ## salinity_median     0.00053353  1.00000000 0.4558  0.672
    ## oxygen_median       0.00110373  1.00000000 0.6086  0.291
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes
(PD_beta_S_ef <- envfit(PD_beta_S_NMDS, env[stratified_lakes,c(2,6,8,26,31:32)], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                             NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median       -0.39651 -0.91803 0.0381  0.882  
    ## salinity_median          -0.96164  0.27432 0.5830  0.129  
    ## oxygen_median             0.83740  0.54660 0.7450  0.042 *
    ## distance_to_ocean_mean_m -0.09986  0.99500 0.0319  0.925  
    ## max_depth                -0.87449  0.48505 0.1022  0.741  
    ## logArea                  -0.40926  0.91242 0.2148  0.538  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_S_efp <- p.adjust.envfit(PD_beta_S_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median       -0.39651 -0.91803 0.0381  1.000
    ## salinity_median          -0.96164  0.27432 0.5830  0.774
    ## oxygen_median             0.83740  0.54660 0.7450  0.252
    ## distance_to_ocean_mean_m -0.09986  0.99500 0.0319  1.000
    ## max_depth                -0.87449  0.48505 0.1022  1.000
    ## logArea                  -0.40926  0.91242 0.2148  1.000
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
# For figure
(PD_beta_geo_ef <- envfit(PD_beta_geo_NMDS, env[surveyed_sites_geo,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.89767  0.44067 0.6055  0.044 *
    ## max_depth               -0.21091 -0.97751 0.0713  0.929  
    ## logArea                  0.23451 -0.97211 0.1882  0.376  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_geo_efp <- p.adjust.envfit(PD_beta_geo_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.89767  0.44067 0.6055  0.132
    ## max_depth               -0.21091 -0.97751 0.0713  1.000
    ## logArea                  0.23451 -0.97211 0.1882  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by stratification
(PD_beta_geo_A_ef <- envfit(PD_beta_geo_NMDS, env[surveyed_sites_geo,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.89767  0.44067 0.6055  0.041 *
    ## max_depth               -0.21091 -0.97751 0.0713  0.923  
    ## logArea                  0.23451 -0.97211 0.1882  0.365  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_geo_A_efp <- p.adjust.envfit(PD_beta_geo_A_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.89767  0.44067 0.6055  0.123
    ## max_depth               -0.21091 -0.97751 0.0713  1.000
    ## logArea                  0.23451 -0.97211 0.1882  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
(PD_beta_geo_MS_ef <- envfit(PD_beta_geo_MS_NMDS, env[mixed_stratified_lakes_geo,geography], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.83618  0.54846 0.5156  0.074 .
    ## max_depth               -0.40273 -0.91532 0.0797  0.960  
    ## logArea                 -0.50467  0.86331 0.0064  0.994  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_geo_MS_efp <- p.adjust.envfit(PD_beta_geo_MS_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.83618  0.54846 0.5156  0.222
    ## max_depth               -0.40273 -0.91532 0.0797  1.000
    ## logArea                 -0.50467  0.86331 0.0064  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
(PD_beta_geo_OM_ef <- envfit(PD_beta_geo_OM_NMDS, env[ocean_mixed_sites_geo,geography], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m  0.59015 -0.80729 0.4463  0.136   
    ## max_depth               -0.89003 -0.45591 0.6830  0.004 **
    ## logArea                 -0.73947  0.67319 0.2464  0.488   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_geo_OM_efp <- p.adjust.envfit(PD_beta_geo_OM_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m  0.59015 -0.80729 0.4463  0.408  
    ## max_depth               -0.89003 -0.45591 0.6830  0.012 *
    ## logArea                 -0.73947  0.67319 0.2464  1.000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
(PD_beta_geo_SO_ef <- envfit(PD_beta_geo_SO_NMDS, env[ocean_stratified_sites,geography], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m  0.99474  0.10248 0.6338  0.212
    ## max_depth                0.73003  0.68342 0.0374  0.975
    ## logArea                 -0.39131  0.92026 0.3079  0.319
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_geo_SO_efp <- p.adjust.envfit(PD_beta_geo_SO_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m  0.99474  0.10248 0.6338  0.636
    ## max_depth                0.73003  0.68342 0.0374  1.000
    ## logArea                 -0.39131  0.92026 0.3079  0.957
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed lakes
(PD_beta_geo_M_ef <- envfit(PD_beta_geo_M_NMDS, env[mixed_lakes_geo,geography], permutations = 999, na.rm = TRUE))
```

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m -0.96293  0.26974 0.5712  0.154   
    ## max_depth                0.68418 -0.72931 0.9636  0.003 **
    ## logArea                  0.96987  0.24360 0.6605  0.115   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 5039

``` r
(PD_beta_geo_M_efp <- p.adjust.envfit(PD_beta_geo_M_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m -0.96293  0.26974 0.5712  0.462   
    ## max_depth                0.68418 -0.72931 0.9636  0.009 **
    ## logArea                  0.96987  0.24360 0.6605  0.345   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 5039

### PD beta Mantel correlation tests

``` r
### Environmental
# Surveyed sites
env_dist_t <- dist(scaled_env[surveyed_sites_env,c(1)], method = "euclidean")
(PD_beta_env_mant_t <- mantel(PD_beta_env_dist$Btotal, env_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_dist$Btotal, ydis = env_dist_t, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1792 
    ##       Significance: 0.212 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.215 0.239 0.269 0.306 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_dist_s <- dist(scaled_env[surveyed_sites_env,c(2)], method = "euclidean")
(PD_beta_env_mant_s <- mantel(PD_beta_env_dist$Btotal, env_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_dist$Btotal, ydis = env_dist_s, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5142 
    ##       Significance: 0.005 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.431 0.457 0.483 0.501 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_dist_o <- dist(scaled_env[surveyed_sites_env,c(3)], method = "euclidean")
(PD_beta_env_mant_o <- mantel(PD_beta_env_dist$Btotal, env_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_dist$Btotal, ydis = env_dist_o, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.462 
    ##       Significance: 0.256 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.529 0.568 0.596 0.638 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_mant_pv <- rbind(PD_beta_env_mant_t$signif, PD_beta_env_mant_s$signif, PD_beta_env_mant_o$signif)
PD_beta_env_mant_pv <- PD_beta_env_mant_pv[,1]
(PD_beta_env_mant_pv <- p.adjust(PD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.636 0.015 0.768

``` r
# Mixed and stratified lakes
env_MS_dist_t <- dist(scaled_env[mixed_stratified_lakes,c(1)], method = "euclidean")
(PD_beta_env_MS_mant_t <- mantel(PD_beta_env_MS_dist$Btotal, env_MS_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_t,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1644 
    ##       Significance: 0.247 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.214 0.253 0.286 0.332 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_MS_dist_s <- dist(scaled_env[mixed_stratified_lakes,c(2)], method = "euclidean")
(PD_beta_env_MS_mant_s <- mantel(PD_beta_env_MS_dist$Btotal, env_MS_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_s,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4638 
    ##       Significance: 0.007 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.375 0.420 0.438 0.455 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_MS_dist_o <- dist(scaled_env[mixed_stratified_lakes,c(3)], method = "euclidean")
(PD_beta_env_MS_mant_o <- mantel(PD_beta_env_MS_dist$Btotal, env_MS_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_o,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3305 
    ##       Significance: 0.289 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.406 0.449 0.497 0.552 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_MS_mant_pv <- rbind(PD_beta_env_MS_mant_t$signif, PD_beta_env_MS_mant_s$signif, PD_beta_env_MS_mant_o$signif)
PD_beta_env_MS_mant_pv <- PD_beta_env_MS_mant_pv[,1]
(PD_beta_env_MS_mant_pv <- p.adjust(PD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.741 0.021 0.867

``` r
# Ocean sites and mixed lakes
env_OM_dist_t <- dist(scaled_env[ocean_mixed_sites_env,c(1)], method = "euclidean")
(PD_beta_env_OM_mant_t <- mantel(PD_beta_env_OM_dist$Btotal, env_OM_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_t,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.03824 
    ##       Significance: 0.528 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.144 0.194 0.251 0.392 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_OM_dist_s <- dist(scaled_env[ocean_mixed_sites_env,c(2)], method = "euclidean")
(PD_beta_env_OM_mant_s <- mantel(PD_beta_env_OM_dist$Btotal, env_OM_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_s,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4565 
    ##       Significance: 0.016 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.288 0.344 0.405 0.486 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_OM_dist_o <- dist(scaled_env[ocean_mixed_sites_env,c(3)], method = "euclidean")
(PD_beta_env_OM_mant_o <- mantel(PD_beta_env_OM_dist$Btotal, env_OM_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_o,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5527 
    ##       Significance: 0.016 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.306 0.384 0.484 0.572 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_OM_mant_pv <- rbind(PD_beta_env_OM_mant_t$signif, PD_beta_env_OM_mant_s$signif, PD_beta_env_OM_mant_o$signif)
PD_beta_env_OM_mant_pv <- PD_beta_env_OM_mant_pv[,1]
(PD_beta_env_OM_mant_pv <- p.adjust(PD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.048 0.048

``` r
# Stratified lakes and ocean sites
env_SO_dist_t <- dist(scaled_env[ocean_stratified_sites_env,c(1)], method = "euclidean")
(PD_beta_env_SO_mant_t <- mantel(PD_beta_env_SO_dist$Btotal, env_SO_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_t,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.07569 
    ##       Significance: 0.466 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0414 0.0742 0.1041 0.1407 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_SO_dist_s <- dist(scaled_env[ocean_stratified_sites_env,c(2)], method = "euclidean")
(PD_beta_env_SO_mant_s <- mantel(PD_beta_env_SO_dist$Btotal, env_SO_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_s,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4755 
    ##       Significance: 0.018 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.394 0.432 0.467 0.492 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_SO_dist_o <- dist(scaled_env[ocean_stratified_sites_env,c(3)], method = "euclidean")
(PD_beta_env_SO_mant_o <- mantel(PD_beta_env_SO_dist$Btotal, env_SO_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_o,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6587 
    ##       Significance: 0.207 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.689 0.708 0.727 0.764 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_SO_mant_pv <- rbind(PD_beta_env_SO_mant_t$signif, PD_beta_env_SO_mant_s$signif, PD_beta_env_SO_mant_o$signif)
PD_beta_env_SO_mant_pv <- PD_beta_env_SO_mant_pv[,1]
(PD_beta_env_SO_mant_pv <- p.adjust(PD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.054 0.621

``` r
# Mixed lakes
env_M_dist_t <- dist(scaled_env[mixed_lakes,c(1)], method = "euclidean")
(PD_beta_env_M_mant_t <- mantel(PD_beta_env_M_dist$Btotal, env_M_dist_t, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_M_dist$Btotal, ydis = env_M_dist_t,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1368 
    ##       Significance: 0.773 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.224 0.316 0.430 0.541 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_M_dist_s <- dist(scaled_env[mixed_lakes,c(2)], method = "euclidean")
(PD_beta_env_M_mant_s <- mantel(PD_beta_env_M_dist$Btotal, env_M_dist_s, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_M_dist$Btotal, ydis = env_M_dist_s,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1095 
    ##       Significance: 0.263 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.226 0.289 0.326 0.371 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_M_dist_o <- dist(scaled_env[mixed_lakes,c(3)], method = "euclidean")
(PD_beta_env_M_mant_o <- mantel(PD_beta_env_M_dist$Btotal, env_M_dist_o, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_env_M_dist$Btotal, ydis = env_M_dist_o,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2748 
    ##       Significance: 0.083 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.237 0.354 0.455 0.566 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_M_mant_pv <- rbind(PD_beta_env_M_mant_t$signif, PD_beta_env_M_mant_s$signif, PD_beta_env_M_mant_o$signif)
PD_beta_env_M_mant_pv <- PD_beta_env_M_mant_pv[,1]
(PD_beta_env_M_mant_pv <- p.adjust(PD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.789 0.249

``` r
# Stratified lakes
env_geo_S_dist <- dist(scaled_env[stratified_lakes,c(1:3,9,14:15)], method = "euclidean")
(PD_beta_S_mant <- mantel(PD_beta_S_dist$Btotal, env_geo_S_dist, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_S_dist$Btotal, ydis = env_geo_S_dist, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1511 
    ##       Significance: 0.28 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.370 0.460 0.538 0.611 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
geo_dist_dmean <- dist(scaled_env[surveyed_sites_geo,c(9)], method = "euclidean")
(PD_beta_geo_mant_dmean <- mantel(PD_beta_geo_dist$Btotal, geo_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_dist$Btotal, ydis = geo_dist_dmean,      method = "spearman", permutations = 999, strata = env[surveyed_sites_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3791 
    ##       Significance: 0.551 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.490 0.516 0.539 0.569 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_md <- dist(scaled_env[surveyed_sites_geo,c(14)], method = "euclidean")
(PD_beta_geo_mant_md <- mantel(PD_beta_geo_dist$Btotal, geo_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_dist$Btotal, ydis = geo_dist_md, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.01139 
    ##       Significance: 0.81 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.159 0.184 0.211 0.248 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_la <- dist(scaled_env[surveyed_sites_geo,c(15)], method = "euclidean")
(PD_beta_geo_mant_la <- mantel(PD_beta_geo_dist$Btotal, geo_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_dist$Btotal, ydis = geo_dist_la, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.03412 
    ##       Significance: 0.763 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.159 0.184 0.202 0.216 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_mant_pv <- rbind(PD_beta_geo_mant_dmean$signif, PD_beta_geo_mant_md$signif, PD_beta_geo_mant_la$signif)
PD_beta_geo_mant_pv <- PD_beta_geo_mant_pv[,1]
(PD_beta_geo_mant_pv <- p.adjust(PD_beta_geo_mant_pv, method = "bonferroni"))
```

    ## [1] 1 1 1

``` r
# Mixed and stratified lakes 
geo_MS_dist_dmean <- dist(scaled_env[mixed_stratified_lakes_geo,c(9)], method = "euclidean")
(PD_beta_geo_MS_mant_dmean <- mantel(PD_beta_geo_MS_dist$Btotal, geo_MS_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_dmean,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1427 
    ##       Significance: 0.588 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.295 0.345 0.376 0.422 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_md <- dist(scaled_env[mixed_stratified_lakes_geo,c(14)], method = "euclidean")
(PD_beta_geo_MS_mant_md <- mantel(PD_beta_geo_MS_dist$Btotal, geo_MS_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_md,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.06657 
    ##       Significance: 0.938 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.156 0.198 0.224 0.256 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_la <- dist(scaled_env[mixed_stratified_lakes_geo,c(15)], method = "euclidean")
(PD_beta_geo_MS_mant_la <- mantel(PD_beta_geo_MS_dist$Btotal, geo_MS_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_la,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.07261 
    ##       Significance: 0.92 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.179 0.213 0.243 0.272 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_MS_mant_pv <- rbind(PD_beta_geo_MS_mant_dmean$signif, PD_beta_geo_MS_mant_md$signif, PD_beta_geo_MS_mant_la$signif)
PD_beta_geo_MS_mant_pv <- PD_beta_geo_MS_mant_pv[,1]
(PD_beta_geo_MS_mant_pv <- p.adjust(PD_beta_geo_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 1 1 1

``` r
# Ocean sites and mixed lakes
geo_OM_dist_dmean <- dist(scaled_env[ocean_mixed_sites_geo,c(9)], method = "euclidean")
(PD_beta_geo_OM_mant_dmean <- mantel(PD_beta_geo_OM_dist$Btotal, geo_OM_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_dmean,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.06768 
    ##       Significance: 0.406 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.198 0.243 0.282 0.333 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_md <- dist(scaled_env[ocean_mixed_sites_geo,c(14)], method = "euclidean")
(PD_beta_geo_OM_mant_md <- mantel(PD_beta_geo_OM_dist$Btotal, geo_OM_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_md,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3254 
    ##       Significance: 0.014 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.163 0.221 0.285 0.356 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_la <- dist(scaled_env[ocean_mixed_sites_geo,c(15)], method = "euclidean")
(PD_beta_geo_OM_mant_la <- mantel(PD_beta_geo_OM_dist$Btotal, geo_OM_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_la,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1664 
    ##       Significance: 0.192 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.199 0.224 0.274 0.316 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_OM_mant_pv <- rbind( PD_beta_geo_OM_mant_dmean$signif, PD_beta_geo_OM_mant_md$signif, PD_beta_geo_OM_mant_la$signif)
PD_beta_geo_OM_mant_pv <- PD_beta_geo_OM_mant_pv[,1]
(PD_beta_geo_OM_mant_pv <- p.adjust(PD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.042 0.576

``` r
# Stratified lakes and ocean sites
geo_SO_dist_dmean <- dist(scaled_env[ocean_stratified_sites,c(9)], method = "euclidean")
(PD_beta_geo_SO_mant_dmean <- mantel(PD_beta_geo_SO_dist$Btotal, geo_SO_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_dmean,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5262 
    ##       Significance: 0.534 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.622 0.655 0.674 0.690 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_SO_dist_md <- dist(scaled_env[ocean_stratified_sites,c(14)], method = "euclidean")
(PD_beta_geo_SO_mant_md <- mantel(PD_beta_geo_SO_dist$Btotal, geo_SO_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_md,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.04194 
    ##       Significance: 0.78 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.106 0.135 0.164 0.193 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_SO_dist_la <- dist(scaled_env[ocean_stratified_sites,c(15)], method = "euclidean")
(PD_beta_geo_SO_mant_la <- mantel(PD_beta_geo_SO_dist$Btotal, geo_SO_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_la,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3193 
    ##       Significance: 0.199 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.367 0.398 0.425 0.443 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_SO_mant_pv <- rbind( PD_beta_geo_SO_mant_dmean$signif, PD_beta_geo_SO_mant_md$signif, PD_beta_geo_SO_mant_la$signif)
PD_beta_geo_SO_mant_pv <- PD_beta_geo_SO_mant_pv[,1]
(PD_beta_geo_SO_mant_pv <- p.adjust(PD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.597

``` r
# Mixed lakes
geo_M_dist_dmean <- dist(scaled_env[mixed_lakes_geo,c(9)], method = "euclidean")
(PD_beta_geo_M_mant_dmean <- mantel(PD_beta_geo_M_dist$Btotal, geo_M_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_dmean,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1429 
    ##       Significance: 0.7 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.346 0.455 0.606 0.770 
    ## Permutation: free
    ## Number of permutations: 5039

``` r
geo_M_dist_md <- dist(scaled_env[mixed_lakes_geo,c(14)], method = "euclidean")
(PD_beta_geo_M_mant_md <- mantel(PD_beta_geo_M_dist$Btotal, geo_M_dist_md, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_md,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6367 
    ##       Significance: 0.019 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.288 0.483 0.607 0.746 
    ## Permutation: free
    ## Number of permutations: 5039

``` r
geo_M_dist_la <- dist(scaled_env[mixed_lakes_geo,c(15)], method = "euclidean")
(PD_beta_geo_M_mant_la <- mantel(PD_beta_geo_M_dist$Btotal, geo_M_dist_la, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_la,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2532 
    ##       Significance: 0.083 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.233 0.316 0.369 0.429 
    ## Permutation: free
    ## Number of permutations: 5039

``` r
# Adjust p-values
PD_beta_geo_M_mant_pv <- rbind( PD_beta_geo_M_mant_dmean$signif, PD_beta_geo_M_mant_md$signif, PD_beta_geo_M_mant_la$signif)
PD_beta_geo_M_mant_pv <- PD_beta_geo_M_mant_pv[,1]
(PD_beta_geo_M_mant_pv <- p.adjust(PD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.057 0.249

### PD beta NMDS ordination plots

``` r
# PD beta total NMDS scores
PD_beta_NMDS_data.scores <- as.data.frame(scores(PD_beta_NMDS))
PD_beta_NMDS_data.scores$Stratification <- env[surveyed_sites,19]
PD_beta_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
PD_beta_NMDS_data.scores$Stratification <- factor(PD_beta_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# Get significantly correlated environmental variables
PD_beta_env_ef_coord_cont <- as.data.frame(scores(PD_beta_env_ef, "vectors")) * ordiArrowMul(PD_beta_env_ef)
PD_beta_env_row_names <- c("S")
# Assign the new row names to the data frame
PD_beta_env_ef_coord_cont <- data.frame(row.names = PD_beta_env_row_names, PD_beta_env_ef_coord_cont)

# PD_beta_geo_ef_coord_cont <- as.data.frame(scores(PD_beta_geo_ef, "vectors")) * ordiArrowMul(PD_beta_geo_ef)
# PD_beta_geo_row_names <- c("minD")
# # Assign the new row names to the data frame
# PD_beta_geo_ef_coord_cont <- data.frame(row.names = PD_beta_geo_row_names, PD_beta_geo_ef_coord_cont)

# Plot NMDS ordination of PD beta total with CI = 0.95
PD_beta_ef_plot <- ggplot(data = PD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.95) +
  geom_point(data = PD_beta_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = PD_beta_env_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = PD_beta_env_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
            color = env_cont, label = row.names(PD_beta_env_ef_coord_cont), size = 5) +
  geom_text_repel(data = PD_beta_NMDS_data.scores, label = PD_beta_NMDS_data.scores$Lakes, 
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
           label = paste("Stress: ", round(PD_beta_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81))
(PD_beta_ef_plot <- PD_beta_ef_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20beta%20NMDS%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_beta_NMDS.jpg", PD_beta_ef_plot, width = 6, height = 6, units = "in")


# Plot NMDS ordination of PD beta total with CI = 0.90
PD_beta_ef_plot_CI90 <- ggplot(data = PD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.90) +
  geom_point(data = PD_beta_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = PD_beta_env_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = PD_beta_env_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
            color = env_cont, label = row.names(PD_beta_env_ef_coord_cont), size = 5) +
  geom_text_repel(data = PD_beta_NMDS_data.scores, label = PD_beta_NMDS_data.scores$Lakes, 
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
           label = paste("Stress: ", round(PD_beta_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81))
(PD_beta_ef_plot_CI90 <- PD_beta_ef_plot_CI90 + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20beta%20NMDS%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_beta_NMDS_CI90.jpg", PD_beta_ef_plot_CI90, width = 6, height = 6, units = "in")


# PD beta replacement NMDS scores
PD_beta_rep_NMDS_data.scores <- as.data.frame(scores(PD_beta_rep_NMDS))
PD_beta_rep_NMDS_data.scores$Stratification <- env[surveyed_sites,19]
PD_beta_rep_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
PD_beta_rep_NMDS_data.scores$Stratification <- factor(PD_beta_rep_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of PD beta replacement with CI = 0.95
PD_beta_rep_plot <- ggplot(data = PD_beta_rep_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.95) +
  geom_point(data = PD_beta_rep_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = PD_beta_rep_NMDS_data.scores, label = PD_beta_rep_NMDS_data.scores$Lakes, 
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
           label = paste("Stress: ", round(PD_beta_rep_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "a") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(PD_beta_rep_plot <- PD_beta_rep_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20beta%20NMDS%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_beta_rep_NMDS.jpg", PD_beta_rep_plot, width = 4.88, height = 6, units = "in")


# PD beta richness NMDS scores
PD_beta_ric_NMDS_data.scores <- as.data.frame(scores(PD_beta_ric_NMDS))
PD_beta_ric_NMDS_data.scores$Stratification <- env[surveyed_sites,19]
PD_beta_ric_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
PD_beta_ric_NMDS_data.scores$Stratification <- factor(PD_beta_ric_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of PD beta richness with CI = 0.95
PD_beta_ric_plot <- ggplot(data = PD_beta_ric_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.95) +
  geom_point(data = PD_beta_ric_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = PD_beta_ric_NMDS_data.scores, label = PD_beta_ric_NMDS_data.scores$Lakes, 
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
           label = paste("Stress: ", round(PD_beta_ric_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "b") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(PD_beta_ric_plot <- PD_beta_ric_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20beta%20NMDS%20plots-4.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_beta_ric_NMDS.jpg", PD_beta_ric_plot, width = 4.88, height = 6, units = "in")


# PD beta ref NMDS scores
PD_beta_ref_NMDS_data.scores <- as.data.frame(scores(PD_beta_ref_NMDS))
PD_beta_ref_NMDS_data.scores$Stratification <- env[,19]
PD_beta_ref_NMDS_data.scores$Lakes <- env[,1]
PD_beta_ref_NMDS_data.scores$Stratification <- factor(PD_beta_ref_NMDS_data.scores$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of PD beta ref with CI = 0.95
PD_beta_ref_plot <- ggplot(data = PD_beta_ref_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.95) +
  geom_point(data = PD_beta_ref_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = PD_beta_ref_NMDS_data.scores, label = PD_beta_ref_NMDS_data.scores$Lakes, 
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
           label = paste("Stress: ", round(PD_beta_ref_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "b") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(PD_beta_ref_plot <- PD_beta_ref_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_path()`).

    ## Warning: ggrepel: 13 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](PD_analyses_files/figure-gfm/PD%20beta%20NMDS%20plots-5.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_beta_ref_NMDS.jpg", PD_beta_ref_plot, width = 4.88, height = 6, units = "in")
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_path()`).
    ## ggrepel: 13 unlabeled data points (too many overlaps). Consider increasing max.overlaps

### PD beta dendrogram

``` r
# Cluster communities using average-linkage algorithm
PD_beta_dist_clust <- hclust(PD_beta_dist$Btotal, method = "average")

# The following will create a dendrogram of PD_beta putting things in one dimensional space
# Open a jpg device
png("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_beta_dendrogram.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

# Create your plot using the plot() function
plot(PD_beta_dist_clust, 
     xlab = "Sites",
     ylab = "PD beta",
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
png("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_beta_dg_silhouette.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

Si <- numeric(nrow(presabs_lake[surveyed_sites,]))
for (k in 2:(nrow(presabs_lake[surveyed_sites,])-1))
{
  sil<-silhouette(cutree(PD_beta_dist_clust, k=k), PD_beta_dist$Btotal)
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
png("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_beta_dg_partitioning.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

PD_beta_pam <- pam(PD_beta_dist$Btotal,2)
plot(PD_beta_pam)

# Close the jpg device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### PD beta outliers

``` r
# Determine NMDS1 outliers
ordered_PD_beta_NMDS1_data.scores <- PD_beta_NMDS_data.scores[order(PD_beta_NMDS_data.scores$NMDS1), ]

outlier_PD_beta_NMDS1 <- ordered_PD_beta_NMDS1_data.scores %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(NMDS1, 0.25),
    Q3 = quantile(NMDS1, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = NMDS1 < lower_bound | NMDS1 > upper_bound
  )

outlier_PD_beta_NMDS1 <- as.data.frame(outlier_PD_beta_NMDS1)

row.names(outlier_PD_beta_NMDS1) <- outlier_PD_beta_NMDS1$X

# Create the plot
(outlier_PD_beta_NMDS1_plot <- ggplot(outlier_PD_beta_NMDS1, aes(x = Stratification, y = NMDS1, fill = Stratification)) +
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
  labs(y = "NMDS1 Distances", x = "Site type", color = "Outlier", fill = "Site type", tag = "a"))
```

![](PD_analyses_files/figure-gfm/PD%20beta%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/outlier_PD_beta_NMDS1.jpg", outlier_PD_beta_NMDS1_plot, width = 6.26, height = 6, units = "in")


# Determine NMDS2 outliers
ordered_PD_beta_NMDS2_data.scores <- PD_beta_NMDS_data.scores[order(PD_beta_NMDS_data.scores$NMDS2), ]

outlier_PD_beta_NMDS2 <- ordered_PD_beta_NMDS2_data.scores %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(NMDS2, 0.25),
    Q3 = quantile(NMDS2, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = NMDS2 < lower_bound | NMDS2 > upper_bound
  )

outlier_PD_beta_NMDS2 <- as.data.frame(outlier_PD_beta_NMDS2)

row.names(outlier_PD_beta_NMDS2) <- outlier_PD_beta_NMDS2$X

# Create the plot
(outlier_PD_beta_NMDS2_plot <- ggplot(outlier_PD_beta_NMDS2, aes(x = Stratification, y = NMDS2, fill = Stratification)) +
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
  labs(y = "NMDS2 Distances", x = "Site type", color = "Outlier", fill = "Site type", tag = "b"))
```

![](PD_analyses_files/figure-gfm/PD%20beta%20outliers-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/outlier_PD_beta_NMDS2.jpg", outlier_PD_beta_NMDS2_plot, width = 6.26, height = 6, units = "in")


# Make PD_beta distances into a matrix
PD_beta_dist_matrix <- as.matrix(PD_beta_dist$Btotal)

# Set the diagonal elements to NA
diag(PD_beta_dist_matrix) <- NA

# Get the labels of the sites
sites <- attr(PD_beta_dist_matrix, "Labels")

# Calculate the mean dist for each group
PD_beta_mean_dist <- aggregate(PD_beta_dist_matrix, by = list(stratification_group), FUN = mean, na.rm = TRUE)

# Rename rows, delete column, and transpose matrix and make it a data frame
row.names(PD_beta_mean_dist) <- PD_beta_mean_dist$Group.1
PD_beta_mean_dist <- PD_beta_mean_dist[,-1] 
PD_beta_mean_dist <- as.data.frame(t(PD_beta_mean_dist))

# Add Stratification column
PD_beta_mean_dist$Stratification <- env[surveyed_sites,19]
PD_beta_mean_dist$Stratification <- factor(PD_beta_mean_dist$Stratification, levels = c("Ocean", "Mixed", "Stratified"))


# Determine ocean sites outliers based on distance from other ocean sites
outlier_PD_beta_mean_dist_ocean <- PD_beta_mean_dist[ocean_sites,] %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(Ocean, 0.25),
    Q3 = quantile(Ocean, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = Ocean < lower_bound | Ocean > upper_bound
  )

outlier_PD_beta_mean_dist_ocean <- as.data.frame(outlier_PD_beta_mean_dist_ocean)
row.names(outlier_PD_beta_mean_dist_ocean) <- outlier_PD_beta_mean_dist_ocean$X
outlier_PD_beta_mean_dist_ocean$Distances <- outlier_PD_beta_mean_dist_ocean$Ocean
outlier_PD_beta_mean_dist_ocean <- outlier_PD_beta_mean_dist_ocean[,-c(1:3)]


# Determine mixed lake outliers based on distance from other mixed lakes
outlier_PD_beta_mean_dist_mixed <- PD_beta_mean_dist[mixed_lakes,] %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(Mixed, 0.25),
    Q3 = quantile(Mixed, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = Mixed < lower_bound | Mixed > upper_bound
  )

outlier_PD_beta_mean_dist_mixed <- as.data.frame(outlier_PD_beta_mean_dist_mixed)
row.names(outlier_PD_beta_mean_dist_mixed) <- outlier_PD_beta_mean_dist_mixed$X
outlier_PD_beta_mean_dist_mixed$Distances <- outlier_PD_beta_mean_dist_mixed$Mixed
outlier_PD_beta_mean_dist_mixed <- outlier_PD_beta_mean_dist_mixed[,-c(1:3)]


# Determine stratified lake outliers based on distance from other stratified lakes
outlier_PD_beta_mean_dist_stratified <- PD_beta_mean_dist[stratified_lakes,] %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(Stratified, 0.25),
    Q3 = quantile(Stratified, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = Stratified < lower_bound | Stratified > upper_bound
  )

outlier_PD_beta_mean_dist_stratified <- as.data.frame(outlier_PD_beta_mean_dist_stratified)
row.names(outlier_PD_beta_mean_dist_stratified) <- outlier_PD_beta_mean_dist_stratified$X
outlier_PD_beta_mean_dist_stratified$Distances <- outlier_PD_beta_mean_dist_stratified$Stratified
outlier_PD_beta_mean_dist_stratified <- outlier_PD_beta_mean_dist_stratified[,-c(1:3)]

outlier_PD_beta_mean_dist <- rbind(outlier_PD_beta_mean_dist_ocean, outlier_PD_beta_mean_dist_mixed, outlier_PD_beta_mean_dist_stratified)

(outlier_PD_beta_mean_dist_plot <- ggplot(outlier_PD_beta_mean_dist, aes(x = Stratification, y = Distances, fill = Stratification)) +
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
  labs(y = "Average Distance from other sites by site type", x = "Site type", color = "Outlier", fill = "Site type", tag = "c"))
```

![](PD_analyses_files/figure-gfm/PD%20beta%20outliers-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/outlier_PD_beta_mean_dist.jpg", outlier_PD_beta_mean_dist_plot, width = 6.26, height = 6, units = "in")
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
    ##  [1] pairwiseAdonis_0.4.1 cluster_2.1.8        BAT_2.9.6           
    ##  [4] caret_7.0-1          ggrepel_0.9.6        viridis_0.6.5       
    ##  [7] viridisLite_0.4.2    ggplot2_3.5.1        picante_1.8.2       
    ## [10] nlme_3.1-167         vegan_2.6-10         lattice_0.22-6      
    ## [13] permute_0.9-7        car_3.1-3            carData_3.0-5       
    ## [16] tidyr_1.3.1          phytools_2.4-4       maps_3.4.2.1        
    ## [19] ape_5.8-1            reshape2_1.4.4       stringr_1.5.1       
    ## [22] dplyr_1.1.4          knitr_1.49          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rstudioapi_0.17.1       magrittr_2.0.3          TH.data_1.1-3          
    ##   [4] farver_2.1.2            rmarkdown_2.29          ragg_1.3.3             
    ##   [7] vctrs_0.6.5             base64enc_0.1-3         terra_1.8-29           
    ##  [10] polspline_1.1.25        htmltools_0.5.8.1       progress_1.2.3         
    ##  [13] DEoptim_2.2-8           Formula_1.2-5           pROC_1.18.5            
    ##  [16] parallelly_1.42.0       pracma_2.4.4            KernSmooth_2.23-26     
    ##  [19] htmlwidgets_1.6.4       sandwich_3.1-1          plyr_1.8.9             
    ##  [22] zoo_1.8-13              palmerpenguins_0.1.1    lubridate_1.9.4        
    ##  [25] igraph_2.1.4            lifecycle_1.0.4         iterators_1.0.14       
    ##  [28] pkgconfig_2.0.3         Matrix_1.7-2            R6_2.6.1               
    ##  [31] fastmap_1.2.0           future_1.34.0           magic_1.6-1            
    ##  [34] digest_0.6.37           numDeriv_2016.8-1.1     colorspace_2.1-1       
    ##  [37] textshaping_1.0.0       Hmisc_5.2-2             pdist_1.2.1            
    ##  [40] labeling_0.4.3          clusterGeneration_1.3.8 timechange_0.3.0       
    ##  [43] abind_1.4-8             mgcv_1.9-1              compiler_4.4.3         
    ##  [46] proxy_0.4-27            withr_3.0.2             doParallel_1.0.17      
    ##  [49] backports_1.5.0         htmlTable_2.4.3         optimParallel_1.0-2    
    ##  [52] quantreg_6.00           MASS_7.3-64             lava_1.8.1             
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
    ##  [94] geometry_0.5.2          systemfonts_1.2.1       munsell_0.5.1          
    ##  [97] Rcpp_1.0.14             globals_0.16.3          coda_0.19-4.1          
    ## [100] fastcluster_1.2.6       MatrixModels_0.5-3      gower_1.0.2            
    ## [103] prettyunits_1.2.0       mclust_6.1.1            listenv_0.9.1          
    ## [106] phangorn_2.12.1         mvtnorm_1.3-3           ipred_0.9-15           
    ## [109] scales_1.3.0            prodlim_2024.06.25      e1071_1.7-16           
    ## [112] purrr_1.0.4             crayon_1.5.3            combinat_0.0-8         
    ## [115] rlang_1.1.5             multcomp_1.4-28         fastmatch_1.1-6        
    ## [118] mnormt_2.1.1            hypervolume_3.1.5

#### stree Z score vs p value plots

``` r
stree_sesmpd_plot <- ggplot(data = stree_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = stree_sesmpd_env[,1], size = 5, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "D") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
stree_sesmpd_plot <- stree_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
stree_sesmpd_plot
```

    ## Warning: ggrepel: 3 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](PD_analyses_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/stree_sesmpd_plot_z.jpg", stree_sesmpd_plot, width = 5, height = 4, units = "in")
```

    ## Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
stree_sesmntd_plot <- ggplot(data = stree_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = stree_sesmntd_env[,1], size = 5, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "E") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
stree_sesmntd_plot <- stree_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
stree_sesmntd_plot
```

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](PD_analyses_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/stree_sesmntd_plot_z.jpg", stree_sesmntd_plot, width = 5, height = 4, units = "in")
```

    ## Warning: ggrepel: 11 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
stree_sespd_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = stree_sespd_env[,1], size = 5, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "F") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
stree_sespd_plot <- stree_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
stree_sespd_plot
```

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](PD_analyses_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/stree_sespd_plot_z.jpg", stree_sespd_plot, width = 5, height = 4, units = "in")
```

    ## Warning: ggrepel: 13 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

#### stree Z score vs obs value plots

``` r
stree_sesmpd_plot <- ggplot(data = stree_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
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
  labs(y= "Observed MPD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
stree_sesmpd_plot <- stree_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
stree_sesmpd_plot
```

    ## Warning: Removed 22 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](PD_analyses_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/stree_sesmpd_plot_obs+z.jpg", stree_sesmpd_plot, width = 5, height = 4, units = "in")
```

    ## Warning: Removed 22 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

``` r
stree_sesmntd_plot <- ggplot(data = stree_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
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
  labs(y= "Observed MNTD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
stree_sesmntd_plot <- stree_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
stree_sesmntd_plot
```

    ## Warning: Removed 22 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](PD_analyses_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/stree_sesmntd_plot_obs+z.jpg", stree_sesmntd_plot, width = 5, height = 4, units = "in")
```

    ## Warning: Removed 22 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

``` r
stree_sespd_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
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
  labs(y= "Observed PD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,6)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=6, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
stree_sespd_plot <- stree_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
stree_sespd_plot
```

    ## Warning: Removed 22 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](PD_analyses_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/stree_sespd_plot_obs+z.jpg", stree_sespd_plot, width = 5, height = 4, units = "in")
```

    ## Warning: Removed 22 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

### Fourth corner approach

``` r
library(mvabund)
library(lattice)
rows_with_na <- which(rowSums(is.na(stree_vif)) > 0)
rows_with_na

ft1=traitglm(surveyed_sites_lake[surveyed_sites_env,], env[surveyed_sites_env,environment], stree_vif[,c(1:13)])
ft1$fourth

a = max(abs(ft1$fourth.corner))
colort = colorRampPalette(c("blue","purple","red")) 
plot.4th = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list(x = list(rot = 45)))
print(plot.4th)

ft <- fourthcorner(env[surveyed_sites_env,environment], surveyed_sites_lake[surveyed_sites_env,], stree_vif[,c(1:5,7:10,12:13)], tr01 = T,modeltype=6)

stree_comp <- stree_vif[complete.cases(stree_vif),]
ssl_comp <- surveyed_sites_lake[,complete.cases(stree_vif)]

ft <- fourthcorner(env[surveyed_sites_env,environment], ssl_comp[surveyed_sites_env,], stree_comp, tr01 = T,modeltype=6)

ft <- fourthcorner(env[surveyed_sites_env,environment], surveyed_sites_lake[surveyed_sites_env,], stree_vif[,c(1:5,7:10,12:13)], tr01 = T,modeltype=6)

stree_comp <- stree_vif[complete.cases(stree_vif),]
ssl_comp <- surveyed_sites_lake[,complete.cases(stree_vif)]

ft <- fourthcorner(env[surveyed_sites_env,environment], ssl_comp[surveyed_sites_env,], stree_comp, tr01 = T,modeltype=6)
```
