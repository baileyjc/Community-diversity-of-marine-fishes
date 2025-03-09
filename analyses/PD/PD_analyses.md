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
PD_PRic_env_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_am)
```

    ##                                        Df Sum Sq Mean Sq F value  Pr(>F)    
    ## distance_to_ocean_min_m                 1  4.688   4.688  63.091 4.6e-05 ***
    ## max_depth                               1  0.013   0.013   0.175 0.68693    
    ## logArea                                 1  0.906   0.906  12.198 0.00817 ** 
    ## Stratification                          2  1.775   0.888  11.944 0.00396 ** 
    ## distance_to_ocean_min_m:Stratification  2  0.025   0.012   0.168 0.84841    
    ## max_depth:Stratification                2  0.271   0.135   1.822 0.22287    
    ## logArea:Stratification                  1  0.004   0.004   0.050 0.82917    
    ## Residuals                               8  0.594   0.074                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_env_lm <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_lm)
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
p_values <- summary(PD_PRic_env_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##       1.0000000000       1.0000000000       0.0001013439       0.0130007650

``` r
PD_z_env_am <- aov(pd.obs.z ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_am)
```

    ##                                        Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1  2.589  2.5886   4.043 0.0792 .
    ## max_depth                               1  1.596  1.5957   2.492 0.1531  
    ## logArea                                 1  0.163  0.1634   0.255 0.6270  
    ## Stratification                          2  2.381  1.1906   1.860 0.2172  
    ## distance_to_ocean_min_m:Stratification  2  0.159  0.0794   0.124 0.8850  
    ## max_depth:Stratification                2  5.030  2.5151   3.928 0.0648 .
    ## logArea:Stratification                  1  1.610  1.6104   2.515 0.1514  
    ## Residuals                               8  5.122  0.6402                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_env_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_env_lm <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_lm)
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
p_values <- summary(PD_z_env_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_am <- aov(Dispersion ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_am)
```

    ##                                        Df    Sum Sq   Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1 0.0012153 0.0012153   3.852 0.0853 .
    ## max_depth                               1 0.0001601 0.0001601   0.507 0.4965  
    ## logArea                                 1 0.0000480 0.0000480   0.152 0.7068  
    ## Stratification                          2 0.0010606 0.0005303   1.681 0.2458  
    ## distance_to_ocean_min_m:Stratification  2 0.0002363 0.0001181   0.374 0.6991  
    ## max_depth:Stratification                2 0.0024191 0.0012096   3.833 0.0680 .
    ## logArea:Stratification                  1 0.0022907 0.0022907   7.260 0.0273 *
    ## Residuals                               8 0.0025243 0.0003155                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_env_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_env_lm <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_lm)
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
p_values <- summary(PD_PDisp_env_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          0.8142048          0.2291871          0.3655132          1.0000000

``` r
# Mixed and stratified lakes
PD_PRic_env_MS_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_MS_am)
```

    ##                                        Df Sum Sq Mean Sq F value  Pr(>F)    
    ## distance_to_ocean_min_m                 1  4.688   4.688  63.091 4.6e-05 ***
    ## max_depth                               1  0.013   0.013   0.175 0.68693    
    ## logArea                                 1  0.906   0.906  12.198 0.00817 ** 
    ## Stratification                          2  1.775   0.888  11.944 0.00396 ** 
    ## distance_to_ocean_min_m:Stratification  2  0.025   0.012   0.168 0.84841    
    ## max_depth:Stratification                2  0.271   0.135   1.822 0.22287    
    ## logArea:Stratification                  1  0.004   0.004   0.050 0.82917    
    ## Residuals                               8  0.594   0.074                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_env_MS_lm <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_env_MS_lm)
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
p_values <- summary(PD_PRic_env_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##       1.0000000000       1.0000000000       0.0004600793       0.0272028786

``` r
PD_z_env_MS_am <- aov(pd.obs.z ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_MS_am)
```

    ##                                        Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1  2.589  2.5886   4.043 0.0792 .
    ## max_depth                               1  1.596  1.5957   2.492 0.1531  
    ## logArea                                 1  0.163  0.1634   0.255 0.6270  
    ## Stratification                          2  2.381  1.1906   1.860 0.2172  
    ## distance_to_ocean_min_m:Stratification  2  0.159  0.0794   0.124 0.8850  
    ## max_depth:Stratification                2  5.030  2.5151   3.928 0.0648 .
    ## logArea:Stratification                  1  1.610  1.6104   2.515 0.1514  
    ## Residuals                               8  5.122  0.6402                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_env_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_env_MS_lm <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_z_env_MS_lm)
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
p_values <- summary(PD_z_env_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_MS_am <- aov(Dispersion ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_MS_am)
```

    ##                                        Df    Sum Sq   Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1 0.0012153 0.0012153   3.852 0.0853 .
    ## max_depth                               1 0.0001601 0.0001601   0.507 0.4965  
    ## logArea                                 1 0.0000480 0.0000480   0.152 0.7068  
    ## Stratification                          2 0.0010606 0.0005303   1.681 0.2458  
    ## distance_to_ocean_min_m:Stratification  2 0.0002363 0.0001181   0.374 0.6991  
    ## max_depth:Stratification                2 0.0024191 0.0012096   3.833 0.0680 .
    ## logArea:Stratification                  1 0.0022907 0.0022907   7.260 0.0273 *
    ## Residuals                               8 0.0025243 0.0003155                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_env_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_env_MS_lm <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_env_MS_lm)
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
p_values <- summary(PD_PDisp_env_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          0.4002584          0.1343561          0.2759753          1.0000000

``` r
# Ocean sites and mixed lakes
PD_PRic_env_OM_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_OM_am)
```

    ##                                        Df Sum Sq Mean Sq F value  Pr(>F)    
    ## distance_to_ocean_min_m                 1  4.688   4.688  63.091 4.6e-05 ***
    ## max_depth                               1  0.013   0.013   0.175 0.68693    
    ## logArea                                 1  0.906   0.906  12.198 0.00817 ** 
    ## Stratification                          2  1.775   0.888  11.944 0.00396 ** 
    ## distance_to_ocean_min_m:Stratification  2  0.025   0.012   0.168 0.84841    
    ## max_depth:Stratification                2  0.271   0.135   1.822 0.22287    
    ## logArea:Stratification                  1  0.004   0.004   0.050 0.82917    
    ## Residuals                               8  0.594   0.074                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_env_OM_lm <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PRic_env_OM_lm)
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
p_values <- summary(PD_PRic_env_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_z_env_OM_am <- aov(pd.obs.z ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_OM_am)
```

    ##                                        Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1  2.589  2.5886   4.043 0.0792 .
    ## max_depth                               1  1.596  1.5957   2.492 0.1531  
    ## logArea                                 1  0.163  0.1634   0.255 0.6270  
    ## Stratification                          2  2.381  1.1906   1.860 0.2172  
    ## distance_to_ocean_min_m:Stratification  2  0.159  0.0794   0.124 0.8850  
    ## max_depth:Stratification                2  5.030  2.5151   3.928 0.0648 .
    ## logArea:Stratification                  1  1.610  1.6104   2.515 0.1514  
    ## Residuals                               8  5.122  0.6402                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_env_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_env_OM_lm <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_z_env_OM_lm)
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
p_values <- summary(PD_z_env_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          1.0000000          1.0000000          0.4629584

``` r
PD_PDisp_env_OM_am <- aov(Dispersion ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_OM_am)
```

    ##                                        Df    Sum Sq   Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1 0.0012153 0.0012153   3.852 0.0853 .
    ## max_depth                               1 0.0001601 0.0001601   0.507 0.4965  
    ## logArea                                 1 0.0000480 0.0000480   0.152 0.7068  
    ## Stratification                          2 0.0010606 0.0005303   1.681 0.2458  
    ## distance_to_ocean_min_m:Stratification  2 0.0002363 0.0001181   0.374 0.6991  
    ## max_depth:Stratification                2 0.0024191 0.0012096   3.833 0.0680 .
    ## logArea:Stratification                  1 0.0022907 0.0022907   7.260 0.0273 *
    ## Residuals                               8 0.0025243 0.0003155                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_env_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_env_OM_lm <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PDisp_env_OM_lm)
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
p_values <- summary(PD_PDisp_env_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
# Stratified lakes and ocean sites
PD_PRic_env_SO_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_SO_am)
```

    ##                                        Df Sum Sq Mean Sq F value  Pr(>F)    
    ## distance_to_ocean_min_m                 1  4.688   4.688  63.091 4.6e-05 ***
    ## max_depth                               1  0.013   0.013   0.175 0.68693    
    ## logArea                                 1  0.906   0.906  12.198 0.00817 ** 
    ## Stratification                          2  1.775   0.888  11.944 0.00396 ** 
    ## distance_to_ocean_min_m:Stratification  2  0.025   0.012   0.168 0.84841    
    ## max_depth:Stratification                2  0.271   0.135   1.822 0.22287    
    ## logArea:Stratification                  1  0.004   0.004   0.050 0.82917    
    ## Residuals                               8  0.594   0.074                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_env_SO_lm <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PRic_env_SO_lm)
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
p_values <- summary(PD_PRic_env_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##        1.000000000        1.000000000        0.001642753        0.038971567

``` r
PD_z_env_SO_am <- aov(pd.obs.z ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_SO_am)
```

    ##                                        Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1  2.589  2.5886   4.043 0.0792 .
    ## max_depth                               1  1.596  1.5957   2.492 0.1531  
    ## logArea                                 1  0.163  0.1634   0.255 0.6270  
    ## Stratification                          2  2.381  1.1906   1.860 0.2172  
    ## distance_to_ocean_min_m:Stratification  2  0.159  0.0794   0.124 0.8850  
    ## max_depth:Stratification                2  5.030  2.5151   3.928 0.0648 .
    ## logArea:Stratification                  1  1.610  1.6104   2.515 0.1514  
    ## Residuals                               8  5.122  0.6402                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_env_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_env_SO_lm <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_z_env_SO_lm)
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
p_values <- summary(PD_z_env_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          1.0000000          0.7285996          1.0000000

``` r
PD_PDisp_env_SO_am <- aov(Dispersion ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_SO_am)
```

    ##                                        Df    Sum Sq   Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1 0.0012153 0.0012153   3.852 0.0853 .
    ## max_depth                               1 0.0001601 0.0001601   0.507 0.4965  
    ## logArea                                 1 0.0000480 0.0000480   0.152 0.7068  
    ## Stratification                          2 0.0010606 0.0005303   1.681 0.2458  
    ## distance_to_ocean_min_m:Stratification  2 0.0002363 0.0001181   0.374 0.6991  
    ## max_depth:Stratification                2 0.0024191 0.0012096   3.833 0.0680 .
    ## logArea:Stratification                  1 0.0022907 0.0022907   7.260 0.0273 *
    ## Residuals                               8 0.0025243 0.0003155                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_env_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_env_SO_lm <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PDisp_env_SO_lm)
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
p_values <- summary(PD_PDisp_env_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          0.5902869          1.0000000          1.0000000

``` r
# Ocean sites
PD_PRic_env_O_lm <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_sites,])
summary(PD_PRic_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         -0.8962        NaN     NaN      NaN
    ## temperature_median   0.1391        NaN     NaN      NaN
    ## salinity_median      0.1405        NaN     NaN      NaN
    ## oxygen_median            NA         NA      NA       NA
    ## 
    ## Residual standard error: NaN on 0 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
p_values <- summary(PD_PRic_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median 
    ##                NaN                NaN                NaN

``` r
PD_z_env_O_lm <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_sites,])
summary(PD_z_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -197.484        NaN     NaN      NaN
    ## temperature_median   -1.463        NaN     NaN      NaN
    ## salinity_median       7.170        NaN     NaN      NaN
    ## oxygen_median            NA         NA      NA       NA
    ## 
    ## Residual standard error: NaN on 0 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
p_values <- summary(PD_z_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median 
    ##                NaN                NaN                NaN

``` r
PD_PDisp_env_O_lm <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[ocean_sites,])
summary(PD_PDisp_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.518321        NaN     NaN      NaN
    ## temperature_median -0.009953        NaN     NaN      NaN
    ## salinity_median     0.029523        NaN     NaN      NaN
    ## oxygen_median             NA         NA      NA       NA
    ## 
    ## Residual standard error: NaN on 0 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
p_values <- summary(PD_PDisp_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median 
    ##                NaN                NaN                NaN

``` r
# Mixed lakes
PD_PRic_env_M_lm <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_lakes,])
summary(PD_PRic_env_M_lm)
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
p_values <- summary(PD_PRic_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_z_env_M_lm <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_lakes,])
summary(PD_z_env_M_lm)
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
p_values <- summary(PD_z_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_M_lm <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[mixed_lakes,])
summary(PD_PDisp_env_M_lm)
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
p_values <- summary(PD_PDisp_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
# Stratified lakes
PD_PRic_env_S_lm <- lm(log(pd.obs) ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[stratified_lakes,])
summary(PD_PRic_env_S_lm)
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
p_values <- summary(PD_PRic_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          0.7931036          1.0000000          0.7303422          1.0000000

``` r
PD_z_env_S_lm <- lm(pd.obs.z ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[stratified_lakes,])
summary(PD_z_env_S_lm)
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
p_values <- summary(PD_z_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_S_lm <- lm(Dispersion ~ temperature_median + salinity_median + oxygen_median, stree_sespd_env[stratified_lakes,])
summary(PD_PDisp_env_S_lm)
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
p_values <- summary(PD_PDisp_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          0.7405587          1.0000000          1.0000000

``` r
### Geographical
# Surveyed sites
PD_PRic_geo_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites,])
summary(PD_PRic_geo_am)
```

    ##                                        Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m                 1  4.802   4.802  56.041  2.1e-05 ***
    ## max_depth                               1  0.000   0.000   0.000 0.997311    
    ## logArea                                 1  0.078   0.078   0.909 0.362755    
    ## Stratification                          2  2.571   1.286  15.005 0.000975 ***
    ## distance_to_ocean_min_m:Stratification  2  0.034   0.017   0.197 0.824545    
    ## max_depth:Stratification                2  0.138   0.069   0.808 0.472782    
    ## logArea:Stratification                  2  0.346   0.173   2.022 0.183079    
    ## Residuals                              10  0.857   0.086                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_geo_lm <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[surveyed_sites,])
summary(PD_PRic_geo_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[surveyed_sites, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.78051 -0.26900  0.04979  0.20082  0.97999 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              7.715598   0.652234  11.829 6.35e-10 ***
    ## distance_to_ocean_min_m -0.005829   0.001471  -3.964  0.00091 ***
    ## max_depth               -0.003314   0.010811  -0.307  0.76270    
    ## logArea                  0.038634   0.064809   0.596  0.55852    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4683 on 18 degrees of freedom
    ## Multiple R-squared:  0.5528, Adjusted R-squared:  0.4783 
    ## F-statistic: 7.418 on 3 and 18 DF,  p-value: 0.001946

``` r
p_values <- summary(PD_PRic_geo_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##            2.540718e-09            3.640987e-03            1.000000e+00 
    ##                 logArea 
    ##            1.000000e+00

``` r
PD_z_geo_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites,])
summary(PD_z_geo_am)
```

    ##                                        Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m                 1  4.802   4.802  56.041  2.1e-05 ***
    ## max_depth                               1  0.000   0.000   0.000 0.997311    
    ## logArea                                 1  0.078   0.078   0.909 0.362755    
    ## Stratification                          2  2.571   1.286  15.005 0.000975 ***
    ## distance_to_ocean_min_m:Stratification  2  0.034   0.017   0.197 0.824545    
    ## max_depth:Stratification                2  0.138   0.069   0.808 0.472782    
    ## logArea:Stratification                  2  0.346   0.173   2.022 0.183079    
    ## Residuals                              10  0.857   0.086                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_geo_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_geo_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[surveyed_sites,])
summary(PD_z_geo_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[surveyed_sites, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.64778 -0.52435  0.00465  0.45465  1.58573 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)             -0.761374   1.365296  -0.558    0.584
    ## distance_to_ocean_min_m  0.004729   0.003078   1.536    0.142
    ## max_depth               -0.033489   0.022631  -1.480    0.156
    ## logArea                 -0.003605   0.135662  -0.027    0.979
    ## 
    ## Residual standard error: 0.9802 on 18 degrees of freedom
    ## Multiple R-squared:  0.1958, Adjusted R-squared:  0.0618 
    ## F-statistic: 1.461 on 3 and 18 DF,  p-value: 0.2585

``` r
p_values <- summary(PD_z_geo_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               0.5676289               0.6248880 
    ##                 logArea 
    ##               1.0000000

``` r
PD_PDisp_geo_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[surveyed_sites,])
summary(PD_PDisp_geo_am)
```

    ##                                        Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m                 1  4.802   4.802  56.041  2.1e-05 ***
    ## max_depth                               1  0.000   0.000   0.000 0.997311    
    ## logArea                                 1  0.078   0.078   0.909 0.362755    
    ## Stratification                          2  2.571   1.286  15.005 0.000975 ***
    ## distance_to_ocean_min_m:Stratification  2  0.034   0.017   0.197 0.824545    
    ## max_depth:Stratification                2  0.138   0.069   0.808 0.472782    
    ## logArea:Stratification                  2  0.346   0.173   2.022 0.183079    
    ## Residuals                              10  0.857   0.086                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_geo_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_geo_lm <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[surveyed_sites,])
summary(PD_PDisp_geo_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[surveyed_sites, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.060930 -0.009537 -0.001303  0.009854  0.041231 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              1.802e-01  3.099e-02   5.815 1.65e-05 ***
    ## distance_to_ocean_min_m -1.131e-04  6.988e-05  -1.619    0.123    
    ## max_depth                3.318e-04  5.137e-04   0.646    0.527    
    ## logArea                 -8.135e-04  3.079e-03  -0.264    0.795    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02225 on 18 degrees of freedom
    ## Multiple R-squared:  0.1318, Adjusted R-squared:  -0.01291 
    ## F-statistic: 0.9108 on 3 and 18 DF,  p-value: 0.4554

``` r
p_values <- summary(PD_PDisp_geo_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##            6.601688e-05            4.912372e-01            1.000000e+00 
    ##                 logArea 
    ##            1.000000e+00

``` r
# Mixed and stratified lakes
PD_PRic_geo_MS_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_geo_MS_am)
```

    ##                                        Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m                 1  3.604   3.604  48.493 0.000117 ***
    ## max_depth                               1  0.023   0.023   0.308 0.593953    
    ## logArea                                 1  1.285   1.285  17.295 0.003171 ** 
    ## Stratification                          1  1.610   1.610  21.671 0.001634 ** 
    ## distance_to_ocean_min_m:Stratification  1  0.003   0.003   0.039 0.848475    
    ## max_depth:Stratification                1  0.048   0.048   0.648 0.444062    
    ## logArea:Stratification                  1  0.004   0.004   0.050 0.829172    
    ## Residuals                               8  0.594   0.074                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_geo_MS_lm <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_geo_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_stratified_lakes, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7599 -0.2009  0.1229  0.2318  0.5418 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              5.486777   1.083074   5.066 0.000277 ***
    ## distance_to_ocean_min_m -0.006856   0.001649  -4.157 0.001331 ** 
    ## max_depth               -0.031337   0.014641  -2.140 0.053553 .  
    ## logArea                  0.329719   0.126209   2.612 0.022697 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4339 on 12 degrees of freedom
    ## Multiple R-squared:  0.6849, Adjusted R-squared:  0.6061 
    ## F-statistic: 8.695 on 3 and 12 DF,  p-value: 0.002451

``` r
p_values <- summary(PD_PRic_geo_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             0.001108257             0.005324119             0.214210399 
    ##                 logArea 
    ##             0.090789723

``` r
PD_z_geo_MS_am <- aov(pd.obs.z ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_z_geo_MS_am)
```

    ##                                        Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1  0.201  0.2007   0.313 0.5909  
    ## max_depth                               1  0.040  0.0401   0.063 0.8086  
    ## logArea                                 1  0.475  0.4747   0.741 0.4143  
    ## Stratification                          1  0.014  0.0144   0.022 0.8847  
    ## distance_to_ocean_min_m:Stratification  1  0.054  0.0545   0.085 0.7780  
    ## max_depth:Stratification                1  2.287  2.2872   3.572 0.0954 .
    ## logArea:Stratification                  1  1.610  1.6104   2.515 0.1514  
    ## Residuals                               8  5.122  0.6402                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_geo_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_geo_MS_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_z_geo_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_stratified_lakes, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.0506 -0.5132 -0.1355  0.6584  1.2369 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)             -2.604956   2.172122  -1.199    0.254
    ## distance_to_ocean_min_m  0.001755   0.003308   0.531    0.605
    ## max_depth               -0.021551   0.029362  -0.734    0.477
    ## logArea                  0.200379   0.253114   0.792    0.444
    ## 
    ## Residual standard error: 0.8703 on 12 degrees of freedom
    ## Multiple R-squared:  0.07298,    Adjusted R-squared:  -0.1588 
    ## F-statistic: 0.3149 on 3 and 12 DF,  p-value: 0.8144

``` r
p_values <- summary(PD_z_geo_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                       1                       1                       1 
    ##                 logArea 
    ##                       1

``` r
PD_PDisp_geo_MS_am <- aov(Dispersion ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_geo_MS_am)
```

    ##                                        Df    Sum Sq   Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1 0.0019693 0.0019693   6.241 0.0370 *
    ## max_depth                               1 0.0006945 0.0006945   2.201 0.1762  
    ## logArea                                 1 0.0009206 0.0009206   2.917 0.1260  
    ## Stratification                          1 0.0000004 0.0000004   0.001 0.9731  
    ## distance_to_ocean_min_m:Stratification  1 0.0000779 0.0000779   0.247 0.6325  
    ## max_depth:Stratification                1 0.0013180 0.0013180   4.177 0.0752 .
    ## logArea:Stratification                  1 0.0022907 0.0022907   7.260 0.0273 *
    ## Residuals                               8 0.0025243 0.0003155                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_geo_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_geo_MS_lm <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_geo_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_stratified_lakes, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.045321 -0.011717 -0.001919  0.010822  0.037509 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              1.046e-01  5.678e-02   1.842   0.0903 .
    ## distance_to_ocean_min_m -1.953e-04  8.648e-05  -2.258   0.0434 *
    ## max_depth               -1.352e-04  7.676e-04  -0.176   0.8632  
    ## logArea                  8.825e-03  6.617e-03   1.334   0.2071  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02275 on 12 degrees of freedom
    ## Multiple R-squared:  0.3659, Adjusted R-squared:  0.2074 
    ## F-statistic: 2.308 on 3 and 12 DF,  p-value: 0.1284

``` r
p_values <- summary(PD_PDisp_geo_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.3610680               0.1734679               1.0000000 
    ##                 logArea 
    ##               0.8283874

``` r
# Ocean sites and mixed lakes
PD_PRic_geo_OM_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PRic_geo_OM_am)
```

    ##                                        Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1 0.0038  0.0038   0.045 0.8385  
    ## max_depth                               1 0.3603  0.3603   4.268 0.0844 .
    ## logArea                                 1 0.0498  0.0498   0.590 0.4715  
    ## Stratification                          1 0.0243  0.0243   0.288 0.6107  
    ## distance_to_ocean_min_m:Stratification  1 0.0007  0.0007   0.008 0.9300  
    ## max_depth:Stratification                1 0.0842  0.0842   0.997 0.3566  
    ## logArea:Stratification                  1 0.3205  0.3205   3.796 0.0993 .
    ## Residuals                               6 0.5065  0.0844                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_geo_OM_lm <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_mixed_sites,])
summary(PD_PRic_geo_OM_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_mixed_sites, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.44690 -0.17366 -0.03444  0.26718  0.38264 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             7.396661   0.511085  14.472 4.93e-08 ***
    ## distance_to_ocean_min_m 0.001238   0.002330   0.532    0.607    
    ## max_depth               0.014949   0.009471   1.578    0.146    
    ## logArea                 0.033648   0.046129   0.729    0.482    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.306 on 10 degrees of freedom
    ## Multiple R-squared:  0.3066, Adjusted R-squared:  0.09858 
    ## F-statistic: 1.474 on 3 and 10 DF,  p-value: 0.2802

``` r
p_values <- summary(PD_PRic_geo_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##            1.971830e-07            1.000000e+00            5.823038e-01 
    ##                 logArea 
    ##            1.000000e+00

``` r
PD_z_geo_OM_am <- aov(pd.obs.z ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_z_geo_OM_am)
```

    ##                                        Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1  1.387   1.387   1.871 0.2204  
    ## max_depth                               1  8.062   8.062  10.874 0.0165 *
    ## logArea                                 1  0.094   0.094   0.126 0.7346  
    ## Stratification                          1  1.077   1.077   1.453 0.2734  
    ## distance_to_ocean_min_m:Stratification  1  0.077   0.077   0.104 0.7583  
    ## max_depth:Stratification                1  0.858   0.858   1.157 0.3234  
    ## logArea:Stratification                  1  0.223   0.223   0.301 0.6032  
    ## Residuals                               6  4.448   0.741                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_geo_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_geo_OM_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_mixed_sites,])
summary(PD_z_geo_OM_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_mixed_sites, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7581 -0.4528 -0.2380  0.1176  1.7563 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              0.362625   1.365557   0.266   0.7960  
    ## distance_to_ocean_min_m  0.002389   0.006224   0.384   0.7091  
    ## max_depth               -0.078901   0.025306  -3.118   0.0109 *
    ## logArea                 -0.046107   0.123252  -0.374   0.7161  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8175 on 10 degrees of freedom
    ## Multiple R-squared:  0.5881, Adjusted R-squared:  0.4645 
    ## F-statistic: 4.759 on 3 and 10 DF,  p-value: 0.02598

``` r
p_values <- summary(PD_z_geo_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              1.00000000              1.00000000              0.04365647 
    ##                 logArea 
    ##              1.00000000

``` r
PD_PDisp_geo_OM_am <- aov(Dispersion ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PDisp_geo_OM_am)
```

    ##                                        Df    Sum Sq   Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m                 1 0.0000369 3.687e-05   0.434  0.534
    ## max_depth                               1 0.0001463 1.463e-04   1.723  0.237
    ## logArea                                 1 0.0000945 9.450e-05   1.113  0.332
    ## Stratification                          1 0.0002191 2.191e-04   2.579  0.159
    ## distance_to_ocean_min_m:Stratification  1 0.0000172 1.719e-05   0.202  0.669
    ## max_depth:Stratification                1 0.0000351 3.512e-05   0.413  0.544
    ## logArea:Stratification                  1 0.0001612 1.612e-04   1.898  0.217
    ## Residuals                               6 0.0005096 8.494e-05

``` r
p_values <- summary(PD_PDisp_geo_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_geo_OM_lm <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_mixed_sites,])
summary(PD_PDisp_geo_OM_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_mixed_sites, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.008451 -0.005301 -0.003319  0.001950  0.023887 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              1.908e-01  1.621e-02  11.770  3.5e-07 ***
    ## distance_to_ocean_min_m -7.455e-06  7.390e-05  -0.101    0.922    
    ## max_depth               -2.442e-04  3.005e-04  -0.813    0.435    
    ## logArea                 -1.466e-03  1.463e-03  -1.001    0.340    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.009707 on 10 degrees of freedom
    ## Multiple R-squared:  0.2276, Adjusted R-squared:  -0.004077 
    ## F-statistic: 0.9824 on 3 and 10 DF,  p-value: 0.4395

``` r
p_values <- summary(PD_PDisp_geo_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##            1.401993e-06            1.000000e+00            1.000000e+00 
    ##                 logArea 
    ##            1.000000e+00

``` r
# Stratified lakes and ocean sites
PD_PRic_geo_SO_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PRic_geo_SO_am)
```

    ##                                        Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m                 1  4.012   4.012  39.287 0.000766 ***
    ## max_depth                               1  0.043   0.043   0.419 0.541244    
    ## logArea                                 1  0.077   0.077   0.753 0.418950    
    ## Stratification                          1  0.624   0.624   6.109 0.048356 *  
    ## distance_to_ocean_min_m:Stratification  1  0.001   0.001   0.011 0.919385    
    ## max_depth:Stratification                1  0.008   0.008   0.081 0.785905    
    ## logArea:Stratification                  1  0.060   0.060   0.586 0.473152    
    ## Residuals                               6  0.613   0.102                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_geo_SO_lm <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_stratified_sites,])
summary(PD_PRic_geo_SO_lm)
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
p_values <- summary(PD_PRic_geo_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##            1.777047e-06            5.014832e-03            1.000000e+00 
    ##                 logArea 
    ##            1.000000e+00

``` r
PD_z_geo_SO_am <- aov(pd.obs.z ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_z_geo_SO_am)
```

    ##                                        Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m                 1  1.287   1.287   2.746 0.14855   
    ## max_depth                               1  1.433   1.433   3.059 0.13088   
    ## logArea                                 1  0.011   0.011   0.023 0.88351   
    ## Stratification                          1  0.520   0.520   1.111 0.33254   
    ## distance_to_ocean_min_m:Stratification  1  0.582   0.582   1.242 0.30765   
    ## max_depth:Stratification                1  6.959   6.959  14.851 0.00842 **
    ## logArea:Stratification                  1  2.590   2.590   5.527 0.05698 . 
    ## Residuals                               6  2.812   0.469                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_geo_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_geo_SO_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_stratified_sites,])
summary(PD_z_geo_SO_lm)
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
p_values <- summary(PD_z_geo_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                       1                       1                       1 
    ##                 logArea 
    ##                       1

``` r
PD_PDisp_geo_SO_am <- aov(Dispersion ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PDisp_geo_SO_am)
```

    ##                                        Df   Sum Sq  Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1 0.000885 0.000885   2.350 0.1761  
    ## max_depth                               1 0.000625 0.000625   1.660 0.2451  
    ## logArea                                 1 0.000006 0.000006   0.016 0.9033  
    ## Stratification                          1 0.000385 0.000385   1.022 0.3511  
    ## distance_to_ocean_min_m:Stratification  1 0.000011 0.000011   0.031 0.8671  
    ## max_depth:Stratification                1 0.001583 0.001583   4.205 0.0862 .
    ## logArea:Stratification                  1 0.003353 0.003353   8.907 0.0245 *
    ## Residuals                               6 0.002259 0.000376                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_geo_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_geo_SO_lm <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_stratified_sites,])
summary(PD_PDisp_geo_SO_lm)
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
p_values <- summary(PD_PDisp_geo_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.0254187               0.9487510               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
# Ocean sites
PD_PRic_geo_O_lm <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_sites,])
summary(PD_PRic_geo_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ##        IBK        LCN        NCN        OOO        OOM        RCA 
    ##  1.469e-02 -1.162e-01  1.374e-02 -3.058e-01  3.936e-01  2.398e-17 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              8.170431   0.927594   8.808   0.0126 *
    ## distance_to_ocean_min_m -0.008315   0.022033  -0.377   0.7422  
    ## max_depth                0.013225   0.016645   0.795   0.5102  
    ## logArea                 -0.029584   0.082546  -0.358   0.7543  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3622 on 2 degrees of freedom
    ## Multiple R-squared:  0.241,  Adjusted R-squared:  -0.8975 
    ## F-statistic: 0.2117 on 3 and 2 DF,  p-value: 0.8817

``` r
p_values <- summary(PD_PRic_geo_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              0.05058098              1.00000000              1.00000000 
    ##                 logArea 
    ##              1.00000000

``` r
PD_z_geo_O_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_sites,])
summary(PD_z_geo_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ##        IBK        LCN        NCN        OOO        OOM        RCA 
    ##  1.413e-01  7.442e-01 -4.173e-01  9.104e-02 -5.592e-01 -1.244e-16 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              1.21947    1.87239   0.651   0.5817  
    ## distance_to_ocean_min_m -0.02868    0.04447  -0.645   0.5851  
    ## max_depth               -0.09940    0.03360  -2.959   0.0978 .
    ## logArea                 -0.09907    0.16662  -0.595   0.6124  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7311 on 2 degrees of freedom
    ## Multiple R-squared:  0.8916, Adjusted R-squared:  0.7291 
    ## F-statistic: 5.485 on 3 and 2 DF,  p-value: 0.1581

``` r
p_values <- summary(PD_z_geo_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               1.0000000               0.3910896 
    ##                 logArea 
    ##               1.0000000

``` r
PD_PDisp_geo_O_lm <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_sites,])
summary(PD_PDisp_geo_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ##        IBK        LCN        NCN        OOO        OOM        RCA 
    ##  2.187e-03  6.768e-03 -5.057e-03 -6.328e-03  2.430e-03 -1.204e-18 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              0.1969825  0.0200118   9.843   0.0102 *
    ## distance_to_ocean_min_m -0.0001342  0.0004753  -0.282   0.8042  
    ## max_depth               -0.0004631  0.0003591  -1.290   0.3262  
    ## logArea                 -0.0018941  0.0017808  -1.064   0.3989  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.007814 on 2 degrees of freedom
    ## Multiple R-squared:  0.7386, Adjusted R-squared:  0.3465 
    ## F-statistic: 1.884 on 3 and 2 DF,  p-value: 0.3653

``` r
p_values <- summary(PD_PDisp_geo_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              0.04065527              1.00000000              1.00000000 
    ##                 logArea 
    ##              1.00000000

``` r
# Mixed lakes
PD_PRic_geo_M_lm <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_lakes,])
summary(PD_PRic_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##      FLK      HLO      LLN      MLN      NLN      NLU      OLO      ULN 
    ##  0.12605 -0.19765  0.15209  0.03516  0.08471  0.10265 -0.37645  0.07345 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              6.0902122  0.6888556   8.841 0.000904 ***
    ## distance_to_ocean_min_m -0.0023065  0.0040741  -0.566 0.601547    
    ## max_depth               -0.0002914  0.0177993  -0.016 0.987722    
    ## logArea                  0.2163289  0.0913654   2.368 0.077011 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2471 on 4 degrees of freedom
    ## Multiple R-squared:  0.7558, Adjusted R-squared:  0.5726 
    ## F-statistic: 4.126 on 3 and 4 DF,  p-value: 0.1023

``` r
p_values <- summary(PD_PRic_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             0.003614367             1.000000000             1.000000000 
    ##                 logArea 
    ##             0.308043921

``` r
PD_z_geo_M_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_lakes,])
summary(PD_z_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##     FLK     HLO     LLN     MLN     NLN     NLU     OLO     ULN 
    ## -0.2123 -0.7466 -0.5777  0.4045  1.3559 -0.0301  0.3621 -0.5558 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)             -0.746784   2.562853  -0.291    0.785
    ## distance_to_ocean_min_m -0.003258   0.015157  -0.215    0.840
    ## max_depth               -0.066912   0.066221  -1.010    0.369
    ## logArea                  0.106039   0.339921   0.312    0.771
    ## 
    ## Residual standard error: 0.9192 on 4 degrees of freedom
    ## Multiple R-squared:  0.3129, Adjusted R-squared:  -0.2024 
    ## F-statistic: 0.6072 on 3 and 4 DF,  p-value: 0.6446

``` r
p_values <- summary(PD_z_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                       1                       1                       1 
    ##                 logArea 
    ##                       1

``` r
PD_PDisp_geo_M_lm <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_lakes,])
summary(PD_PDisp_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##        FLK        HLO        LLN        MLN        NLN        NLU        OLO 
    ##  0.0038012 -0.0083620 -0.0037686  0.0021595  0.0155470 -0.0003265 -0.0053946 
    ##        ULN 
    ## -0.0036561 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)              0.1652839  0.0274444   6.023  0.00383 **
    ## distance_to_ocean_min_m -0.0002204  0.0001623  -1.358  0.24611   
    ## max_depth               -0.0008034  0.0007091  -1.133  0.32054   
    ## logArea                  0.0036215  0.0036400   0.995  0.37609   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.009843 on 4 degrees of freedom
    ## Multiple R-squared:  0.3334, Adjusted R-squared:  -0.1665 
    ## F-statistic: 0.6669 on 3 and 4 DF,  p-value: 0.615

``` r
p_values <- summary(PD_PDisp_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              0.01531847              0.98443978              1.00000000 
    ##                 logArea 
    ##              1.00000000

``` r
# Stratified lakes
PD_PRic_geo_S_lm <- lm(log(pd.obs) ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[stratified_lakes,])
summary(PD_PRic_geo_S_lm)
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
p_values <- summary(PD_PRic_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.1323530               0.5154496               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
PD_z_geo_S_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[stratified_lakes,])
summary(PD_z_geo_S_lm)
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
p_values <- summary(PD_z_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.2164135               1.0000000               0.7258090 
    ##                 logArea 
    ##               0.3147069

``` r
PD_PDisp_geo_S_lm <- lm(Dispersion ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[stratified_lakes,])
summary(PD_PDisp_geo_S_lm)
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
p_values <- summary(PD_PDisp_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
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
PD_beta_ref_dist <- BAT::beta(presabs_lake, tree, abund = F)
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
    ## 3 Stratified vs Reference  1 0.7729479 7.110647 0.5039207   0.121      0.726
    ## 4          Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.097      0.582
    ## 5      Mixed vs Reference  1 0.6592270 3.967203 0.3617334   0.117      0.702
    ## 6      Ocean vs Reference  1 0.6319662 3.354877 0.4015472   0.137      0.822
    ##   sig
    ## 1   .
    ## 2   *
    ## 3    
    ## 4    
    ## 5    
    ## 6

``` r
# Surveyed sites
PD_beta_dist <- BAT::beta(surveyed_sites_lake, stree, abund = F)
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
    ## 1 Stratified vs Mixed  1 1.2940870 9.415922 0.4021162   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.1858559 8.357071 0.4105242   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.108      0.324

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
    ## 1 Stratified vs Mixed  1 1.2940870 9.415922 0.4021162   0.002      0.006   *
    ## 2 Stratified vs Ocean  1 1.2083230 9.439158 0.4618174   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.3347638 2.034034 0.1560556   0.019      0.057

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
    ## 3      Mixed vs Ocean  1 0.2573592  1.467099 0.1089395   0.094      0.282

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
PD_beta_geo_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites,], stree, abund = F)
# Mixed and stratified lakes
PD_beta_geo_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], stree, abund = F)
# Ocean sites and mixed lakes
PD_beta_geo_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites,], stree, abund = F)
# Stratified lakes and ocean sites
PD_beta_geo_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites,], stree, abund = F)
# Mixed lakes
PD_beta_geo_M_dist<- BAT::beta(surveyed_sites_lake[mixed_lakes,], stree, abund = F)
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
    ## Run 1 stress 9.360011e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001535056  max resid 0.0003589049 
    ## ... Similar to previous best
    ## Run 2 stress 9.961151e-05 
    ## ... Procrustes: rmse 0.0001605938  max resid 0.0003085181 
    ## ... Similar to previous best
    ## Run 3 stress 9.669407e-05 
    ## ... Procrustes: rmse 6.370208e-05  max resid 0.000137317 
    ## ... Similar to previous best
    ## Run 4 stress 9.907165e-05 
    ## ... Procrustes: rmse 0.0001791006  max resid 0.000317618 
    ## ... Similar to previous best
    ## Run 5 stress 9.845563e-05 
    ## ... Procrustes: rmse 0.0001566071  max resid 0.0002906076 
    ## ... Similar to previous best
    ## Run 6 stress 9.919638e-05 
    ## ... Procrustes: rmse 0.0001845691  max resid 0.0004141959 
    ## ... Similar to previous best
    ## Run 7 stress 8.494268e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001434376  max resid 0.0002499643 
    ## ... Similar to previous best
    ## Run 8 stress 9.867027e-05 
    ## ... Procrustes: rmse 8.201786e-05  max resid 0.0001330559 
    ## ... Similar to previous best
    ## Run 9 stress 7.764609e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 6.521635e-05  max resid 0.0001101225 
    ## ... Similar to previous best
    ## Run 10 stress 3.694074e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 4.088427e-05  max resid 7.792976e-05 
    ## ... Similar to previous best
    ## Run 11 stress 9.79299e-05 
    ## ... Procrustes: rmse 0.0002018331  max resid 0.0003461666 
    ## ... Similar to previous best
    ## Run 12 stress 9.854933e-05 
    ## ... Procrustes: rmse 0.0002503834  max resid 0.0004271661 
    ## ... Similar to previous best
    ## Run 13 stress 9.998871e-05 
    ## ... Procrustes: rmse 0.0001918585  max resid 0.0003224949 
    ## ... Similar to previous best
    ## Run 14 stress 9.866826e-05 
    ## ... Procrustes: rmse 0.0002563008  max resid 0.000438526 
    ## ... Similar to previous best
    ## Run 15 stress 9.580757e-05 
    ## ... Procrustes: rmse 0.0001984122  max resid 0.0003399463 
    ## ... Similar to previous best
    ## Run 16 stress 9.578798e-05 
    ## ... Procrustes: rmse 0.0001739413  max resid 0.0002961196 
    ## ... Similar to previous best
    ## Run 17 stress 9.710668e-05 
    ## ... Procrustes: rmse 0.0002094315  max resid 0.000349388 
    ## ... Similar to previous best
    ## Run 18 stress 9.927488e-05 
    ## ... Procrustes: rmse 0.0002552335  max resid 0.0004373452 
    ## ... Similar to previous best
    ## Run 19 stress 9.854536e-05 
    ## ... Procrustes: rmse 0.0002316042  max resid 0.0003882935 
    ## ... Similar to previous best
    ## Run 20 stress 7.489959e-05 
    ## ... Procrustes: rmse 9.911946e-05  max resid 0.0001963151 
    ## ... Similar to previous best
    ## Run 21 stress 9.87025e-05 
    ## ... Procrustes: rmse 0.0002569169  max resid 0.0004355935 
    ## ... Similar to previous best
    ## Run 22 stress 9.503837e-05 
    ## ... Procrustes: rmse 0.0002177298  max resid 0.0003702095 
    ## ... Similar to previous best
    ## Run 23 stress 9.811979e-05 
    ## ... Procrustes: rmse 0.0002497933  max resid 0.0004353706 
    ## ... Similar to previous best
    ## Run 24 stress 7.532461e-05 
    ## ... Procrustes: rmse 3.588584e-05  max resid 7.6888e-05 
    ## ... Similar to previous best
    ## Run 25 stress 9.624196e-05 
    ## ... Procrustes: rmse 0.0002329979  max resid 0.000397415 
    ## ... Similar to previous best
    ## Run 26 stress 9.704299e-05 
    ## ... Procrustes: rmse 0.0002013806  max resid 0.0003356635 
    ## ... Similar to previous best
    ## Run 27 stress 9.802555e-05 
    ## ... Procrustes: rmse 0.000208071  max resid 0.0003617311 
    ## ... Similar to previous best
    ## Run 28 stress 5.757503e-05 
    ## ... Procrustes: rmse 3.10491e-05  max resid 5.599121e-05 
    ## ... Similar to previous best
    ## Run 29 stress 9.640476e-05 
    ## ... Procrustes: rmse 0.0002070022  max resid 0.0003528157 
    ## ... Similar to previous best
    ## Run 30 stress 6.068776e-05 
    ## ... Procrustes: rmse 8.371764e-05  max resid 0.0001397508 
    ## ... Similar to previous best
    ## Run 31 stress 9.908971e-05 
    ## ... Procrustes: rmse 0.0002583135  max resid 0.0004388237 
    ## ... Similar to previous best
    ## Run 32 stress 9.368417e-05 
    ## ... Procrustes: rmse 0.0002406349  max resid 0.0004111136 
    ## ... Similar to previous best
    ## Run 33 stress 7.805511e-05 
    ## ... Procrustes: rmse 4.872592e-05  max resid 6.775139e-05 
    ## ... Similar to previous best
    ## Run 34 stress 9.857582e-05 
    ## ... Procrustes: rmse 0.0002550636  max resid 0.0004282435 
    ## ... Similar to previous best
    ## Run 35 stress 9.75207e-05 
    ## ... Procrustes: rmse 4.990306e-05  max resid 0.0001017851 
    ## ... Similar to previous best
    ## Run 36 stress 9.893627e-05 
    ## ... Procrustes: rmse 0.0002577646  max resid 0.000437316 
    ## ... Similar to previous best
    ## Run 37 stress 9.931417e-05 
    ## ... Procrustes: rmse 0.0002228621  max resid 0.0003818998 
    ## ... Similar to previous best
    ## Run 38 stress 9.805747e-05 
    ## ... Procrustes: rmse 0.000182414  max resid 0.0003179777 
    ## ... Similar to previous best
    ## Run 39 stress 9.434553e-05 
    ## ... Procrustes: rmse 0.0001815026  max resid 0.00031029 
    ## ... Similar to previous best
    ## Run 40 stress 9.880687e-05 
    ## ... Procrustes: rmse 0.0002475055  max resid 0.0004368858 
    ## ... Similar to previous best
    ## Run 41 stress 9.711129e-05 
    ## ... Procrustes: rmse 0.0002498048  max resid 0.0004256832 
    ## ... Similar to previous best
    ## Run 42 stress 9.661648e-05 
    ## ... Procrustes: rmse 0.0002273203  max resid 0.0003868393 
    ## ... Similar to previous best
    ## Run 43 stress 9.928648e-05 
    ## ... Procrustes: rmse 0.0002204858  max resid 0.0003776856 
    ## ... Similar to previous best
    ## Run 44 stress 9.242537e-05 
    ## ... Procrustes: rmse 4.507052e-05  max resid 8.482818e-05 
    ## ... Similar to previous best
    ## Run 45 stress 9.92419e-05 
    ## ... Procrustes: rmse 0.0002407892  max resid 0.0004051126 
    ## ... Similar to previous best
    ## Run 46 stress 9.888563e-05 
    ## ... Procrustes: rmse 0.0002254232  max resid 0.0003887568 
    ## ... Similar to previous best
    ## Run 47 stress 9.277652e-05 
    ## ... Procrustes: rmse 4.298332e-05  max resid 7.091003e-05 
    ## ... Similar to previous best
    ## Run 48 stress 9.782366e-05 
    ## ... Procrustes: rmse 0.0002515322  max resid 0.000426862 
    ## ... Similar to previous best
    ## Run 49 stress 9.437131e-05 
    ## ... Procrustes: rmse 0.0001735259  max resid 0.0002981528 
    ## ... Similar to previous best
    ## Run 50 stress 9.890657e-05 
    ## ... Procrustes: rmse 0.0002112609  max resid 0.0003646327 
    ## ... Similar to previous best
    ## Run 51 stress 9.941657e-05 
    ## ... Procrustes: rmse 0.0002507878  max resid 0.0004225268 
    ## ... Similar to previous best
    ## Run 52 stress 4.79734e-05 
    ## ... Procrustes: rmse 2.885124e-05  max resid 5.609834e-05 
    ## ... Similar to previous best
    ## Run 53 stress 7.202386e-05 
    ## ... Procrustes: rmse 5.665758e-05  max resid 0.0001001583 
    ## ... Similar to previous best
    ## Run 54 stress 9.219945e-05 
    ## ... Procrustes: rmse 9.48731e-05  max resid 0.0001438523 
    ## ... Similar to previous best
    ## Run 55 stress 9.909513e-05 
    ## ... Procrustes: rmse 0.0002236095  max resid 0.0003823384 
    ## ... Similar to previous best
    ## Run 56 stress 9.7724e-05 
    ## ... Procrustes: rmse 0.0002516127  max resid 0.0004299803 
    ## ... Similar to previous best
    ## Run 57 stress 9.800587e-05 
    ## ... Procrustes: rmse 0.0002348827  max resid 0.0004045597 
    ## ... Similar to previous best
    ## Run 58 stress 6.927116e-05 
    ## ... Procrustes: rmse 8.820241e-05  max resid 0.0001668315 
    ## ... Similar to previous best
    ## Run 59 stress 7.653566e-05 
    ## ... Procrustes: rmse 5.752564e-05  max resid 0.0001040362 
    ## ... Similar to previous best
    ## Run 60 stress 9.930912e-05 
    ## ... Procrustes: rmse 0.0002543032  max resid 0.0004332378 
    ## ... Similar to previous best
    ## Run 61 stress 9.788304e-05 
    ## ... Procrustes: rmse 0.0002427356  max resid 0.0004101248 
    ## ... Similar to previous best
    ## Run 62 stress 9.933059e-05 
    ## ... Procrustes: rmse 0.0002157197  max resid 0.0003617311 
    ## ... Similar to previous best
    ## Run 63 stress 9.558761e-05 
    ## ... Procrustes: rmse 0.0002004659  max resid 0.0003347722 
    ## ... Similar to previous best
    ## Run 64 stress 9.915814e-05 
    ## ... Procrustes: rmse 0.0002504648  max resid 0.0004215827 
    ## ... Similar to previous best
    ## Run 65 stress 9.679469e-05 
    ## ... Procrustes: rmse 0.0001314188  max resid 0.0002113032 
    ## ... Similar to previous best
    ## Run 66 stress 9.96127e-05 
    ## ... Procrustes: rmse 0.0002594583  max resid 0.0004382513 
    ## ... Similar to previous best
    ## Run 67 stress 9.839493e-05 
    ## ... Procrustes: rmse 0.000216532  max resid 0.0003719319 
    ## ... Similar to previous best
    ## Run 68 stress 9.870292e-05 
    ## ... Procrustes: rmse 0.000243265  max resid 0.0004118332 
    ## ... Similar to previous best
    ## Run 69 stress 7.831848e-05 
    ## ... Procrustes: rmse 4.400806e-05  max resid 8.851835e-05 
    ## ... Similar to previous best
    ## Run 70 stress 9.830509e-05 
    ## ... Procrustes: rmse 0.0001889089  max resid 0.0003275916 
    ## ... Similar to previous best
    ## Run 71 stress 6.735998e-05 
    ## ... Procrustes: rmse 6.019015e-05  max resid 9.948216e-05 
    ## ... Similar to previous best
    ## Run 72 stress 9.696191e-05 
    ## ... Procrustes: rmse 0.0002528241  max resid 0.0004307472 
    ## ... Similar to previous best
    ## Run 73 stress 9.604667e-05 
    ## ... Procrustes: rmse 0.0002448678  max resid 0.0004174913 
    ## ... Similar to previous best
    ## Run 74 stress 9.853763e-05 
    ## ... Procrustes: rmse 0.0002105043  max resid 0.0003724175 
    ## ... Similar to previous best
    ## Run 75 stress 9.991237e-05 
    ## ... Procrustes: rmse 0.0002241485  max resid 0.0003816639 
    ## ... Similar to previous best
    ## Run 76 stress 9.931075e-05 
    ## ... Procrustes: rmse 0.0002508806  max resid 0.0004307072 
    ## ... Similar to previous best
    ## Run 77 stress 4.615499e-05 
    ## ... Procrustes: rmse 3.108399e-05  max resid 6.410344e-05 
    ## ... Similar to previous best
    ## Run 78 stress 9.867509e-05 
    ## ... Procrustes: rmse 0.0002431369  max resid 0.0004190717 
    ## ... Similar to previous best
    ## Run 79 stress 9.987071e-05 
    ## ... Procrustes: rmse 0.0002146997  max resid 0.0003627796 
    ## ... Similar to previous best
    ## Run 80 stress 9.932065e-05 
    ## ... Procrustes: rmse 0.0002233637  max resid 0.0003792842 
    ## ... Similar to previous best
    ## Run 81 stress 8.679784e-05 
    ## ... Procrustes: rmse 4.229714e-05  max resid 9.746594e-05 
    ## ... Similar to previous best
    ## Run 82 stress 9.888929e-05 
    ## ... Procrustes: rmse 0.0002220693  max resid 0.0003797887 
    ## ... Similar to previous best
    ## Run 83 stress 4.201462e-05 
    ## ... Procrustes: rmse 2.759324e-05  max resid 5.70354e-05 
    ## ... Similar to previous best
    ## Run 84 stress 9.608236e-05 
    ## ... Procrustes: rmse 0.0002378822  max resid 0.0004044281 
    ## ... Similar to previous best
    ## Run 85 stress 7.666896e-05 
    ## ... Procrustes: rmse 6.53727e-05  max resid 0.0001269312 
    ## ... Similar to previous best
    ## Run 86 stress 9.876483e-05 
    ## ... Procrustes: rmse 0.0002523723  max resid 0.0004305345 
    ## ... Similar to previous best
    ## Run 87 stress 9.815518e-05 
    ## ... Procrustes: rmse 9.571456e-05  max resid 0.0001621476 
    ## ... Similar to previous best
    ## Run 88 stress 9.93521e-05 
    ## ... Procrustes: rmse 0.0002573734  max resid 0.0004408369 
    ## ... Similar to previous best
    ## Run 89 stress 8.241511e-05 
    ## ... Procrustes: rmse 4.088425e-05  max resid 7.962027e-05 
    ## ... Similar to previous best
    ## Run 90 stress 9.067765e-05 
    ## ... Procrustes: rmse 0.0001965234  max resid 0.0003431024 
    ## ... Similar to previous best
    ## Run 91 stress 9.442952e-05 
    ## ... Procrustes: rmse 0.0001682501  max resid 0.0002894496 
    ## ... Similar to previous best
    ## Run 92 stress 9.954817e-05 
    ## ... Procrustes: rmse 0.0002066092  max resid 0.0003632674 
    ## ... Similar to previous best
    ## Run 93 stress 9.787946e-05 
    ## ... Procrustes: rmse 0.0002035767  max resid 0.0003401648 
    ## ... Similar to previous best
    ## Run 94 stress 9.327302e-05 
    ## ... Procrustes: rmse 0.0001355532  max resid 0.0002159519 
    ## ... Similar to previous best
    ## Run 95 stress 9.925606e-05 
    ## ... Procrustes: rmse 0.0002088095  max resid 0.0003482846 
    ## ... Similar to previous best
    ## Run 96 stress 9.420215e-05 
    ## ... Procrustes: rmse 0.0001987082  max resid 0.0003385859 
    ## ... Similar to previous best
    ## Run 97 stress 9.69911e-05 
    ## ... Procrustes: rmse 0.0001206472  max resid 0.0002435071 
    ## ... Similar to previous best
    ## Run 98 stress 9.950524e-05 
    ## ... Procrustes: rmse 0.0002262692  max resid 0.0003971705 
    ## ... Similar to previous best
    ## Run 99 stress 9.88838e-05 
    ## ... Procrustes: rmse 0.0002561029  max resid 0.0004354249 
    ## ... Similar to previous best
    ## Run 100 stress 9.930244e-05 
    ## ... Procrustes: rmse 0.0002118774  max resid 0.0003559989 
    ## ... Similar to previous best
    ## Run 101 stress 9.667147e-05 
    ## ... Procrustes: rmse 0.0002186585  max resid 0.0003707913 
    ## ... Similar to previous best
    ## Run 102 stress 9.056034e-05 
    ## ... Procrustes: rmse 8.559067e-05  max resid 0.0001370308 
    ## ... Similar to previous best
    ## Run 103 stress 9.588604e-05 
    ## ... Procrustes: rmse 0.0002466477  max resid 0.0004133483 
    ## ... Similar to previous best
    ## Run 104 stress 9.979683e-05 
    ## ... Procrustes: rmse 0.0002500979  max resid 0.0004275982 
    ## ... Similar to previous best
    ## Run 105 stress 8.165455e-05 
    ## ... Procrustes: rmse 4.27838e-05  max resid 6.645586e-05 
    ## ... Similar to previous best
    ## Run 106 stress 6.096582e-05 
    ## ... Procrustes: rmse 3.692733e-05  max resid 9.550348e-05 
    ## ... Similar to previous best
    ## Run 107 stress 9.599346e-05 
    ## ... Procrustes: rmse 0.0002051149  max resid 0.0003473815 
    ## ... Similar to previous best
    ## Run 108 stress 9.951288e-05 
    ## ... Procrustes: rmse 0.0002581012  max resid 0.0004397622 
    ## ... Similar to previous best
    ## Run 109 stress 9.91633e-05 
    ## ... Procrustes: rmse 0.0001913668  max resid 0.0003282596 
    ## ... Similar to previous best
    ## Run 110 stress 9.855874e-05 
    ## ... Procrustes: rmse 0.0002550208  max resid 0.0004310075 
    ## ... Similar to previous best
    ## Run 111 stress 6.318526e-05 
    ## ... Procrustes: rmse 7.091678e-05  max resid 0.0001238568 
    ## ... Similar to previous best
    ## Run 112 stress 9.707325e-05 
    ## ... Procrustes: rmse 0.0002514925  max resid 0.0004282724 
    ## ... Similar to previous best
    ## Run 113 stress 8.362147e-05 
    ## ... Procrustes: rmse 0.0001245083  max resid 0.000206727 
    ## ... Similar to previous best
    ## Run 114 stress 9.569921e-05 
    ## ... Procrustes: rmse 0.0002147382  max resid 0.0003651557 
    ## ... Similar to previous best
    ## Run 115 stress 9.702849e-05 
    ## ... Procrustes: rmse 0.0002137653  max resid 0.0003738484 
    ## ... Similar to previous best
    ## Run 116 stress 9.605158e-05 
    ## ... Procrustes: rmse 0.0001974178  max resid 0.0003285701 
    ## ... Similar to previous best
    ## Run 117 stress 6.457605e-05 
    ## ... Procrustes: rmse 3.117717e-05  max resid 5.262663e-05 
    ## ... Similar to previous best
    ## Run 118 stress 8.441914e-05 
    ## ... Procrustes: rmse 3.904146e-05  max resid 7.431923e-05 
    ## ... Similar to previous best
    ## Run 119 stress 9.863361e-05 
    ## ... Procrustes: rmse 0.0002212784  max resid 0.0003785022 
    ## ... Similar to previous best
    ## Run 120 stress 9.50561e-05 
    ## ... Procrustes: rmse 9.825023e-05  max resid 0.0001614475 
    ## ... Similar to previous best
    ## Run 121 stress 9.927541e-05 
    ## ... Procrustes: rmse 0.000227861  max resid 0.0003933437 
    ## ... Similar to previous best
    ## Run 122 stress 9.059118e-05 
    ## ... Procrustes: rmse 0.0001826712  max resid 0.000310435 
    ## ... Similar to previous best
    ## Run 123 stress 9.84848e-05 
    ## ... Procrustes: rmse 0.0002543654  max resid 0.0004380348 
    ## ... Similar to previous best
    ## Run 124 stress 9.930197e-05 
    ## ... Procrustes: rmse 0.0002544023  max resid 0.0004374046 
    ## ... Similar to previous best
    ## Run 125 stress 9.569066e-05 
    ## ... Procrustes: rmse 0.0001695202  max resid 0.0002919415 
    ## ... Similar to previous best
    ## Run 126 stress 7.922389e-05 
    ## ... Procrustes: rmse 0.000103412  max resid 0.0001671691 
    ## ... Similar to previous best
    ## Run 127 stress 0.3824858 
    ## Run 128 stress 4.592965e-05 
    ## ... Procrustes: rmse 2.770609e-05  max resid 5.4668e-05 
    ## ... Similar to previous best
    ## Run 129 stress 9.815867e-05 
    ## ... Procrustes: rmse 0.0002495618  max resid 0.0004340181 
    ## ... Similar to previous best
    ## Run 130 stress 9.774978e-05 
    ## ... Procrustes: rmse 0.0002503652  max resid 0.0004300552 
    ## ... Similar to previous best
    ## Run 131 stress 9.593315e-05 
    ## ... Procrustes: rmse 0.0002461878  max resid 0.0004204544 
    ## ... Similar to previous best
    ## Run 132 stress 9.26616e-05 
    ## ... Procrustes: rmse 0.0002332613  max resid 0.0003924423 
    ## ... Similar to previous best
    ## Run 133 stress 7.789668e-05 
    ## ... Procrustes: rmse 4.2937e-05  max resid 7.768734e-05 
    ## ... Similar to previous best
    ## Run 134 stress 7.200821e-05 
    ## ... Procrustes: rmse 3.388747e-05  max resid 7.134744e-05 
    ## ... Similar to previous best
    ## Run 135 stress 9.936517e-05 
    ## ... Procrustes: rmse 0.0002561614  max resid 0.0004376265 
    ## ... Similar to previous best
    ## Run 136 stress 6.437709e-05 
    ## ... Procrustes: rmse 3.549036e-05  max resid 5.52625e-05 
    ## ... Similar to previous best
    ## Run 137 stress 6.458653e-05 
    ## ... Procrustes: rmse 3.76855e-05  max resid 7.157878e-05 
    ## ... Similar to previous best
    ## Run 138 stress 9.985286e-05 
    ## ... Procrustes: rmse 0.0002567106  max resid 0.0004378449 
    ## ... Similar to previous best
    ## Run 139 stress 9.685149e-05 
    ## ... Procrustes: rmse 9.511141e-05  max resid 0.0001476527 
    ## ... Similar to previous best
    ## Run 140 stress 7.475896e-05 
    ## ... Procrustes: rmse 3.992238e-05  max resid 7.980844e-05 
    ## ... Similar to previous best
    ## Run 141 stress 9.914071e-05 
    ## ... Procrustes: rmse 0.0002399046  max resid 0.0004031199 
    ## ... Similar to previous best
    ## Run 142 stress 9.760775e-05 
    ## ... Procrustes: rmse 0.0002170258  max resid 0.0003673 
    ## ... Similar to previous best
    ## Run 143 stress 8.428316e-05 
    ## ... Procrustes: rmse 5.063353e-05  max resid 9.128942e-05 
    ## ... Similar to previous best
    ## Run 144 stress 9.330992e-05 
    ## ... Procrustes: rmse 0.0001213676  max resid 0.0001983273 
    ## ... Similar to previous best
    ## Run 145 stress 8.95018e-05 
    ## ... Procrustes: rmse 4.682343e-05  max resid 6.785441e-05 
    ## ... Similar to previous best
    ## Run 146 stress 9.900088e-05 
    ## ... Procrustes: rmse 0.000218925  max resid 0.0003725511 
    ## ... Similar to previous best
    ## Run 147 stress 9.782947e-05 
    ## ... Procrustes: rmse 0.0002100959  max resid 0.0003627898 
    ## ... Similar to previous best
    ## Run 148 stress 9.80279e-05 
    ## ... Procrustes: rmse 0.0002147439  max resid 0.0003747943 
    ## ... Similar to previous best
    ## Run 149 stress 9.997104e-05 
    ## ... Procrustes: rmse 0.000203466  max resid 0.000347 
    ## ... Similar to previous best
    ## Run 150 stress 9.926102e-05 
    ## ... Procrustes: rmse 0.0002498192  max resid 0.0004252888 
    ## ... Similar to previous best
    ## Run 151 stress 9.490423e-05 
    ## ... Procrustes: rmse 4.678535e-05  max resid 8.237461e-05 
    ## ... Similar to previous best
    ## Run 152 stress 7.827685e-05 
    ## ... Procrustes: rmse 4.389653e-05  max resid 7.273757e-05 
    ## ... Similar to previous best
    ## Run 153 stress 9.636699e-05 
    ## ... Procrustes: rmse 0.0002171013  max resid 0.0003708538 
    ## ... Similar to previous best
    ## Run 154 stress 9.992651e-05 
    ## ... Procrustes: rmse 0.000248793  max resid 0.0004120952 
    ## ... Similar to previous best
    ## Run 155 stress 0.3825725 
    ## Run 156 stress 9.668725e-05 
    ## ... Procrustes: rmse 0.0002149017  max resid 0.0003658662 
    ## ... Similar to previous best
    ## Run 157 stress 9.905314e-05 
    ## ... Procrustes: rmse 0.0002177074  max resid 0.0003741212 
    ## ... Similar to previous best
    ## Run 158 stress 9.936714e-05 
    ## ... Procrustes: rmse 0.0002221386  max resid 0.0003770663 
    ## ... Similar to previous best
    ## Run 159 stress 9.495844e-05 
    ## ... Procrustes: rmse 5.143632e-05  max resid 7.683242e-05 
    ## ... Similar to previous best
    ## Run 160 stress 9.73585e-05 
    ## ... Procrustes: rmse 0.0002323826  max resid 0.0003921593 
    ## ... Similar to previous best
    ## Run 161 stress 0.379374 
    ## Run 162 stress 9.699927e-05 
    ## ... Procrustes: rmse 0.0002450638  max resid 0.0004174334 
    ## ... Similar to previous best
    ## Run 163 stress 9.822515e-05 
    ## ... Procrustes: rmse 0.0002501694  max resid 0.0004230112 
    ## ... Similar to previous best
    ## Run 164 stress 9.823293e-05 
    ## ... Procrustes: rmse 0.0002442658  max resid 0.0004105502 
    ## ... Similar to previous best
    ## Run 165 stress 5.03157e-05 
    ## ... Procrustes: rmse 2.861582e-05  max resid 4.749896e-05 
    ## ... Similar to previous best
    ## Run 166 stress 9.917754e-05 
    ## ... Procrustes: rmse 0.0002226329  max resid 0.0003842569 
    ## ... Similar to previous best
    ## Run 167 stress 9.814159e-05 
    ## ... Procrustes: rmse 0.000219861  max resid 0.0003775888 
    ## ... Similar to previous best
    ## Run 168 stress 9.823194e-05 
    ## ... Procrustes: rmse 0.0002409219  max resid 0.0004101569 
    ## ... Similar to previous best
    ## Run 169 stress 9.755901e-05 
    ## ... Procrustes: rmse 0.000249127  max resid 0.00043133 
    ## ... Similar to previous best
    ## Run 170 stress 8.637616e-05 
    ## ... Procrustes: rmse 4.994969e-05  max resid 9.442219e-05 
    ## ... Similar to previous best
    ## Run 171 stress 9.952872e-05 
    ## ... Procrustes: rmse 0.0002584529  max resid 0.000434101 
    ## ... Similar to previous best
    ## Run 172 stress 9.650287e-05 
    ## ... Procrustes: rmse 0.0002500637  max resid 0.0004222336 
    ## ... Similar to previous best
    ## Run 173 stress 7.451423e-05 
    ## ... Procrustes: rmse 8.408999e-05  max resid 0.000175114 
    ## ... Similar to previous best
    ## Run 174 stress 9.774508e-05 
    ## ... Procrustes: rmse 0.0002536882  max resid 0.0004348809 
    ## ... Similar to previous best
    ## Run 175 stress 9.038597e-05 
    ## ... Procrustes: rmse 8.474026e-05  max resid 0.0001500634 
    ## ... Similar to previous best
    ## Run 176 stress 6.300928e-05 
    ## ... Procrustes: rmse 8.107653e-05  max resid 0.0001321155 
    ## ... Similar to previous best
    ## Run 177 stress 9.047546e-05 
    ## ... Procrustes: rmse 0.0002156623  max resid 0.0003628019 
    ## ... Similar to previous best
    ## Run 178 stress 5.200824e-05 
    ## ... Procrustes: rmse 4.221106e-05  max resid 9.451375e-05 
    ## ... Similar to previous best
    ## Run 179 stress 9.999057e-05 
    ## ... Procrustes: rmse 5.214272e-05  max resid 8.423346e-05 
    ## ... Similar to previous best
    ## Run 180 stress 8.215097e-05 
    ## ... Procrustes: rmse 4.304171e-05  max resid 7.318375e-05 
    ## ... Similar to previous best
    ## Run 181 stress 9.342082e-05 
    ## ... Procrustes: rmse 0.0002010784  max resid 0.0003462126 
    ## ... Similar to previous best
    ## Run 182 stress 8.446534e-05 
    ## ... Procrustes: rmse 0.0001269386  max resid 0.0002199803 
    ## ... Similar to previous best
    ## Run 183 stress 9.438611e-05 
    ## ... Procrustes: rmse 0.0001043314  max resid 0.0002000143 
    ## ... Similar to previous best
    ## Run 184 stress 9.410516e-05 
    ## ... Procrustes: rmse 0.0002108017  max resid 0.0003560832 
    ## ... Similar to previous best
    ## Run 185 stress 9.94579e-05 
    ## ... Procrustes: rmse 0.0001874141  max resid 0.0003279439 
    ## ... Similar to previous best
    ## Run 186 stress 3.86908e-05 
    ## ... Procrustes: rmse 2.343915e-05  max resid 3.944036e-05 
    ## ... Similar to previous best
    ## Run 187 stress 9.648538e-05 
    ## ... Procrustes: rmse 0.0002118671  max resid 0.0003630053 
    ## ... Similar to previous best
    ## Run 188 stress 9.770942e-05 
    ## ... Procrustes: rmse 0.0002022187  max resid 0.000339834 
    ## ... Similar to previous best
    ## Run 189 stress 9.97722e-05 
    ## ... Procrustes: rmse 0.0002213132  max resid 0.0003799301 
    ## ... Similar to previous best
    ## Run 190 stress 9.318852e-05 
    ## ... Procrustes: rmse 0.0002141782  max resid 0.0003643913 
    ## ... Similar to previous best
    ## Run 191 stress 9.673536e-05 
    ## ... Procrustes: rmse 0.0002518847  max resid 0.0004277522 
    ## ... Similar to previous best
    ## Run 192 stress 9.943632e-05 
    ## ... Procrustes: rmse 0.0002511517  max resid 0.0004289561 
    ## ... Similar to previous best
    ## Run 193 stress 9.832591e-05 
    ## ... Procrustes: rmse 0.0002207943  max resid 0.0003849899 
    ## ... Similar to previous best
    ## Run 194 stress 9.83244e-05 
    ## ... Procrustes: rmse 0.0002203045  max resid 0.0003786136 
    ## ... Similar to previous best
    ## Run 195 stress 7.513916e-05 
    ## ... Procrustes: rmse 3.554693e-05  max resid 7.422714e-05 
    ## ... Similar to previous best
    ## Run 196 stress 9.909737e-05 
    ## ... Procrustes: rmse 0.0002548609  max resid 0.0004262807 
    ## ... Similar to previous best
    ## Run 197 stress 9.991418e-05 
    ## ... Procrustes: rmse 0.0002564433  max resid 0.0004369105 
    ## ... Similar to previous best
    ## Run 198 stress 9.456518e-05 
    ## ... Procrustes: rmse 0.0002228167  max resid 0.0003856774 
    ## ... Similar to previous best
    ## Run 199 stress 9.683775e-05 
    ## ... Procrustes: rmse 0.00023449  max resid 0.0004045101 
    ## ... Similar to previous best
    ## Run 200 stress 8.738009e-05 
    ## ... Procrustes: rmse 5.207757e-05  max resid 9.061565e-05 
    ## ... Similar to previous best
    ## Run 201 stress 9.670283e-05 
    ## ... Procrustes: rmse 0.0002106986  max resid 0.0003616303 
    ## ... Similar to previous best
    ## Run 202 stress 9.982167e-05 
    ## ... Procrustes: rmse 0.0002119068  max resid 0.0003581607 
    ## ... Similar to previous best
    ## Run 203 stress 9.640857e-05 
    ## ... Procrustes: rmse 5.233961e-05  max resid 7.540841e-05 
    ## ... Similar to previous best
    ## Run 204 stress 9.880212e-05 
    ## ... Procrustes: rmse 0.0002517329  max resid 0.0004305944 
    ## ... Similar to previous best
    ## Run 205 stress 5.806031e-05 
    ## ... Procrustes: rmse 3.067432e-05  max resid 5.425634e-05 
    ## ... Similar to previous best
    ## Run 206 stress 9.465845e-05 
    ## ... Procrustes: rmse 0.000119129  max resid 0.0001785382 
    ## ... Similar to previous best
    ## Run 207 stress 9.545209e-05 
    ## ... Procrustes: rmse 0.0002475027  max resid 0.0004245251 
    ## ... Similar to previous best
    ## Run 208 stress 8.940663e-05 
    ## ... Procrustes: rmse 0.0002110646  max resid 0.0003654296 
    ## ... Similar to previous best
    ## Run 209 stress 7.30817e-05 
    ## ... Procrustes: rmse 3.956612e-05  max resid 7.173354e-05 
    ## ... Similar to previous best
    ## Run 210 stress 9.112806e-05 
    ## ... Procrustes: rmse 0.0001808027  max resid 0.0003024456 
    ## ... Similar to previous best
    ## Run 211 stress 9.826456e-05 
    ## ... Procrustes: rmse 5.800713e-05  max resid 0.0001049579 
    ## ... Similar to previous best
    ## Run 212 stress 0.3824866 
    ## Run 213 stress 8.981198e-05 
    ## ... Procrustes: rmse 0.0001385316  max resid 0.0002408135 
    ## ... Similar to previous best
    ## Run 214 stress 9.723795e-05 
    ## ... Procrustes: rmse 0.0002516233  max resid 0.0004317634 
    ## ... Similar to previous best
    ## Run 215 stress 9.628354e-05 
    ## ... Procrustes: rmse 0.0002507634  max resid 0.0004256898 
    ## ... Similar to previous best
    ## Run 216 stress 9.987893e-05 
    ## ... Procrustes: rmse 0.0002461493  max resid 0.000408419 
    ## ... Similar to previous best
    ## Run 217 stress 9.614376e-05 
    ## ... Procrustes: rmse 0.0002002042  max resid 0.0003359154 
    ## ... Similar to previous best
    ## Run 218 stress 9.78662e-05 
    ## ... Procrustes: rmse 0.0002540821  max resid 0.0004358914 
    ## ... Similar to previous best
    ## Run 219 stress 4.189585e-05 
    ## ... Procrustes: rmse 2.411342e-05  max resid 5.038697e-05 
    ## ... Similar to previous best
    ## Run 220 stress 9.625142e-05 
    ## ... Procrustes: rmse 0.0002175472  max resid 0.0003673861 
    ## ... Similar to previous best
    ## Run 221 stress 9.994536e-05 
    ## ... Procrustes: rmse 0.0002599241  max resid 0.0004441364 
    ## ... Similar to previous best
    ## Run 222 stress 9.957392e-05 
    ## ... Procrustes: rmse 0.0002497567  max resid 0.0004189665 
    ## ... Similar to previous best
    ## Run 223 stress 9.78777e-05 
    ## ... Procrustes: rmse 0.0002181603  max resid 0.0003750342 
    ## ... Similar to previous best
    ## Run 224 stress 6.259636e-05 
    ## ... Procrustes: rmse 3.69236e-05  max resid 5.837608e-05 
    ## ... Similar to previous best
    ## Run 225 stress 5.883572e-05 
    ## ... Procrustes: rmse 3.743581e-05  max resid 7.430803e-05 
    ## ... Similar to previous best
    ## Run 226 stress 9.90178e-05 
    ## ... Procrustes: rmse 0.0002113516  max resid 0.000361221 
    ## ... Similar to previous best
    ## Run 227 stress 9.855055e-05 
    ## ... Procrustes: rmse 0.0002533798  max resid 0.0004321459 
    ## ... Similar to previous best
    ## Run 228 stress 9.895819e-05 
    ## ... Procrustes: rmse 0.0002548487  max resid 0.0004303301 
    ## ... Similar to previous best
    ## Run 229 stress 9.859202e-05 
    ## ... Procrustes: rmse 0.0002251641  max resid 0.000386524 
    ## ... Similar to previous best
    ## Run 230 stress 5.637411e-05 
    ## ... Procrustes: rmse 4.202261e-05  max resid 6.953963e-05 
    ## ... Similar to previous best
    ## Run 231 stress 9.717046e-05 
    ## ... Procrustes: rmse 0.0002188593  max resid 0.0003766916 
    ## ... Similar to previous best
    ## Run 232 stress 9.914957e-05 
    ## ... Procrustes: rmse 0.0002029969  max resid 0.0003393384 
    ## ... Similar to previous best
    ## Run 233 stress 9.958994e-05 
    ## ... Procrustes: rmse 0.0001971251  max resid 0.000321023 
    ## ... Similar to previous best
    ## Run 234 stress 9.938297e-05 
    ## ... Procrustes: rmse 0.0002169622  max resid 0.0003762691 
    ## ... Similar to previous best
    ## Run 235 stress 9.933764e-05 
    ## ... Procrustes: rmse 0.0002273625  max resid 0.0003895094 
    ## ... Similar to previous best
    ## Run 236 stress 9.767599e-05 
    ## ... Procrustes: rmse 0.0002544468  max resid 0.0004325455 
    ## ... Similar to previous best
    ## Run 237 stress 6.410364e-05 
    ## ... Procrustes: rmse 3.941202e-05  max resid 6.525675e-05 
    ## ... Similar to previous best
    ## Run 238 stress 8.767643e-05 
    ## ... Procrustes: rmse 0.0001981836  max resid 0.0003323794 
    ## ... Similar to previous best
    ## Run 239 stress 9.935296e-05 
    ## ... Procrustes: rmse 0.0002586218  max resid 0.0004397087 
    ## ... Similar to previous best
    ## Run 240 stress 5.273878e-05 
    ## ... Procrustes: rmse 2.800744e-05  max resid 5.756668e-05 
    ## ... Similar to previous best
    ## Run 241 stress 7.109115e-05 
    ## ... Procrustes: rmse 9.266727e-05  max resid 0.0001549748 
    ## ... Similar to previous best
    ## Run 242 stress 9.916325e-05 
    ## ... Procrustes: rmse 0.0002584978  max resid 0.0004407956 
    ## ... Similar to previous best
    ## Run 243 stress 5.761356e-05 
    ## ... Procrustes: rmse 3.528665e-05  max resid 6.43623e-05 
    ## ... Similar to previous best
    ## Run 244 stress 6.292414e-05 
    ## ... Procrustes: rmse 3.557971e-05  max resid 6.757418e-05 
    ## ... Similar to previous best
    ## Run 245 stress 9.878863e-05 
    ## ... Procrustes: rmse 0.000221827  max resid 0.0003813605 
    ## ... Similar to previous best
    ## Run 246 stress 5.668739e-05 
    ## ... Procrustes: rmse 5.97278e-05  max resid 9.105964e-05 
    ## ... Similar to previous best
    ## Run 247 stress 9.926078e-05 
    ## ... Procrustes: rmse 0.0002218766  max resid 0.0003636038 
    ## ... Similar to previous best
    ## Run 248 stress 8.13296e-05 
    ## ... Procrustes: rmse 3.606959e-05  max resid 8.854282e-05 
    ## ... Similar to previous best
    ## Run 249 stress 7.725756e-05 
    ## ... Procrustes: rmse 3.810898e-05  max resid 8.47285e-05 
    ## ... Similar to previous best
    ## Run 250 stress 9.676255e-05 
    ## ... Procrustes: rmse 0.000228989  max resid 0.0003788012 
    ## ... Similar to previous best
    ## Run 251 stress 9.898778e-05 
    ## ... Procrustes: rmse 0.0001802209  max resid 0.0003236758 
    ## ... Similar to previous best
    ## Run 252 stress 5.369138e-05 
    ## ... Procrustes: rmse 3.427886e-05  max resid 6.937907e-05 
    ## ... Similar to previous best
    ## Run 253 stress 9.991141e-05 
    ## ... Procrustes: rmse 0.0002570413  max resid 0.000427705 
    ## ... Similar to previous best
    ## Run 254 stress 8.689193e-05 
    ## ... Procrustes: rmse 9.189055e-05  max resid 0.0001498835 
    ## ... Similar to previous best
    ## Run 255 stress 7.779894e-05 
    ## ... Procrustes: rmse 5.795479e-05  max resid 9.058224e-05 
    ## ... Similar to previous best
    ## Run 256 stress 9.484913e-05 
    ## ... Procrustes: rmse 0.0002188183  max resid 0.0003712765 
    ## ... Similar to previous best
    ## Run 257 stress 9.968288e-05 
    ## ... Procrustes: rmse 0.0002188925  max resid 0.0003683929 
    ## ... Similar to previous best
    ## Run 258 stress 7.533363e-05 
    ## ... Procrustes: rmse 4.628613e-05  max resid 7.796296e-05 
    ## ... Similar to previous best
    ## Run 259 stress 6.482759e-05 
    ## ... Procrustes: rmse 4.337704e-05  max resid 7.496384e-05 
    ## ... Similar to previous best
    ## Run 260 stress 9.822737e-05 
    ## ... Procrustes: rmse 0.0002219933  max resid 0.0003781972 
    ## ... Similar to previous best
    ## Run 261 stress 9.801519e-05 
    ## ... Procrustes: rmse 0.0002216091  max resid 0.0003827232 
    ## ... Similar to previous best
    ## Run 262 stress 9.964626e-05 
    ## ... Procrustes: rmse 0.0002587006  max resid 0.0004408955 
    ## ... Similar to previous best
    ## Run 263 stress 8.533412e-05 
    ## ... Procrustes: rmse 6.782554e-05  max resid 0.0001137622 
    ## ... Similar to previous best
    ## Run 264 stress 9.629801e-05 
    ## ... Procrustes: rmse 0.0002201943  max resid 0.0003785789 
    ## ... Similar to previous best
    ## Run 265 stress 9.261166e-05 
    ## ... Procrustes: rmse 0.0001958273  max resid 0.0003367299 
    ## ... Similar to previous best
    ## Run 266 stress 7.48272e-05 
    ## ... Procrustes: rmse 0.0001024525  max resid 0.0001770529 
    ## ... Similar to previous best
    ## Run 267 stress 9.877581e-05 
    ## ... Procrustes: rmse 0.0001827519  max resid 0.0003226032 
    ## ... Similar to previous best
    ## Run 268 stress 9.125912e-05 
    ## ... Procrustes: rmse 0.0002089743  max resid 0.0003579297 
    ## ... Similar to previous best
    ## Run 269 stress 9.881783e-05 
    ## ... Procrustes: rmse 0.0002506671  max resid 0.0004290276 
    ## ... Similar to previous best
    ## Run 270 stress 9.884678e-05 
    ## ... Procrustes: rmse 0.0001738164  max resid 0.0002932788 
    ## ... Similar to previous best
    ## Run 271 stress 8.989931e-05 
    ## ... Procrustes: rmse 0.0002344119  max resid 0.0004012846 
    ## ... Similar to previous best
    ## Run 272 stress 9.364456e-05 
    ## ... Procrustes: rmse 0.0001585114  max resid 0.0002611252 
    ## ... Similar to previous best
    ## Run 273 stress 8.472735e-05 
    ## ... Procrustes: rmse 5.24252e-05  max resid 8.45672e-05 
    ## ... Similar to previous best
    ## Run 274 stress 9.719984e-05 
    ## ... Procrustes: rmse 0.0002280048  max resid 0.0003891825 
    ## ... Similar to previous best
    ## Run 275 stress 9.978615e-05 
    ## ... Procrustes: rmse 0.0002204156  max resid 0.0003782455 
    ## ... Similar to previous best
    ## Run 276 stress 6.900095e-05 
    ## ... Procrustes: rmse 7.867822e-05  max resid 0.0001426504 
    ## ... Similar to previous best
    ## Run 277 stress 6.36934e-05 
    ## ... Procrustes: rmse 5.692482e-05  max resid 9.922093e-05 
    ## ... Similar to previous best
    ## Run 278 stress 5.65238e-05 
    ## ... Procrustes: rmse 3.220989e-05  max resid 5.952241e-05 
    ## ... Similar to previous best
    ## Run 279 stress 9.641257e-05 
    ## ... Procrustes: rmse 0.0002127562  max resid 0.0003665713 
    ## ... Similar to previous best
    ## Run 280 stress 9.886315e-05 
    ## ... Procrustes: rmse 4.996671e-05  max resid 0.0001055529 
    ## ... Similar to previous best
    ## Run 281 stress 9.852877e-05 
    ## ... Procrustes: rmse 0.000254569  max resid 0.0004306184 
    ## ... Similar to previous best
    ## Run 282 stress 7.118891e-05 
    ## ... Procrustes: rmse 3.507742e-05  max resid 7.222012e-05 
    ## ... Similar to previous best
    ## Run 283 stress 8.418657e-05 
    ## ... Procrustes: rmse 4.819052e-05  max resid 7.734934e-05 
    ## ... Similar to previous best
    ## Run 284 stress 6.952893e-05 
    ## ... Procrustes: rmse 4.30693e-05  max resid 7.256192e-05 
    ## ... Similar to previous best
    ## Run 285 stress 9.397569e-05 
    ## ... Procrustes: rmse 0.0001694856  max resid 0.0002913566 
    ## ... Similar to previous best
    ## Run 286 stress 7.154575e-05 
    ## ... Procrustes: rmse 4.358818e-05  max resid 6.435621e-05 
    ## ... Similar to previous best
    ## Run 287 stress 9.316355e-05 
    ## ... Procrustes: rmse 0.0001840842  max resid 0.0003092712 
    ## ... Similar to previous best
    ## Run 288 stress 0.3770684 
    ## Run 289 stress 9.957461e-05 
    ## ... Procrustes: rmse 0.0002107347  max resid 0.0003420153 
    ## ... Similar to previous best
    ## Run 290 stress 9.445996e-05 
    ## ... Procrustes: rmse 0.000214004  max resid 0.0003693695 
    ## ... Similar to previous best
    ## Run 291 stress 9.876463e-05 
    ## ... Procrustes: rmse 0.0002087299  max resid 0.0003535588 
    ## ... Similar to previous best
    ## Run 292 stress 9.768535e-05 
    ## ... Procrustes: rmse 0.0001884198  max resid 0.0003444252 
    ## ... Similar to previous best
    ## Run 293 stress 8.657666e-05 
    ## ... Procrustes: rmse 4.120743e-05  max resid 6.995689e-05 
    ## ... Similar to previous best
    ## Run 294 stress 9.853854e-05 
    ## ... Procrustes: rmse 0.0002529901  max resid 0.0004308363 
    ## ... Similar to previous best
    ## Run 295 stress 6.715428e-05 
    ## ... Procrustes: rmse 3.05729e-05  max resid 6.739124e-05 
    ## ... Similar to previous best
    ## Run 296 stress 6.431834e-05 
    ## ... Procrustes: rmse 3.605753e-05  max resid 7.423996e-05 
    ## ... Similar to previous best
    ## Run 297 stress 8.523819e-05 
    ## ... Procrustes: rmse 4.683167e-05  max resid 7.021155e-05 
    ## ... Similar to previous best
    ## Run 298 stress 9.5352e-05 
    ## ... Procrustes: rmse 0.000199346  max resid 0.0003338916 
    ## ... Similar to previous best
    ## Run 299 stress 9.636711e-05 
    ## ... Procrustes: rmse 0.0002298062  max resid 0.0003947585 
    ## ... Similar to previous best
    ## Run 300 stress 7.458342e-05 
    ## ... Procrustes: rmse 4.95797e-05  max resid 8.569102e-05 
    ## ... Similar to previous best
    ## Run 301 stress 9.852231e-05 
    ## ... Procrustes: rmse 0.0001499552  max resid 0.0002401067 
    ## ... Similar to previous best
    ## Run 302 stress 9.567064e-05 
    ## ... Procrustes: rmse 0.0001978575  max resid 0.0003346474 
    ## ... Similar to previous best
    ## Run 303 stress 7.830459e-05 
    ## ... Procrustes: rmse 4.773872e-05  max resid 8.505674e-05 
    ## ... Similar to previous best
    ## Run 304 stress 9.844468e-05 
    ## ... Procrustes: rmse 0.0002523539  max resid 0.0004292843 
    ## ... Similar to previous best
    ## Run 305 stress 9.810467e-05 
    ## ... Procrustes: rmse 0.0002211472  max resid 0.0003828182 
    ## ... Similar to previous best
    ## Run 306 stress 9.796896e-05 
    ## ... Procrustes: rmse 0.0002024992  max resid 0.0003349894 
    ## ... Similar to previous best
    ## Run 307 stress 9.801301e-05 
    ## ... Procrustes: rmse 0.0002185516  max resid 0.000368411 
    ## ... Similar to previous best
    ## Run 308 stress 9.895074e-05 
    ## ... Procrustes: rmse 0.0002520823  max resid 0.0004312059 
    ## ... Similar to previous best
    ## Run 309 stress 9.904726e-05 
    ## ... Procrustes: rmse 0.0002519248  max resid 0.0004286317 
    ## ... Similar to previous best
    ## Run 310 stress 4.863784e-05 
    ## ... Procrustes: rmse 2.907213e-05  max resid 5.96238e-05 
    ## ... Similar to previous best
    ## Run 311 stress 7.369557e-05 
    ## ... Procrustes: rmse 4.090496e-05  max resid 8.395109e-05 
    ## ... Similar to previous best
    ## Run 312 stress 9.833069e-05 
    ## ... Procrustes: rmse 0.0002527178  max resid 0.0004331883 
    ## ... Similar to previous best
    ## Run 313 stress 9.922392e-05 
    ## ... Procrustes: rmse 0.0002271554  max resid 0.0003906848 
    ## ... Similar to previous best
    ## Run 314 stress 9.721934e-05 
    ## ... Procrustes: rmse 0.0002513371  max resid 0.0004288997 
    ## ... Similar to previous best
    ## Run 315 stress 5.801237e-05 
    ## ... Procrustes: rmse 3.296488e-05  max resid 5.454913e-05 
    ## ... Similar to previous best
    ## Run 316 stress 9.93995e-05 
    ## ... Procrustes: rmse 0.0002545491  max resid 0.0004262248 
    ## ... Similar to previous best
    ## Run 317 stress 9.883015e-05 
    ## ... Procrustes: rmse 0.0002168193  max resid 0.000358989 
    ## ... Similar to previous best
    ## Run 318 stress 9.892951e-05 
    ## ... Procrustes: rmse 0.0002540742  max resid 0.0004357756 
    ## ... Similar to previous best
    ## Run 319 stress 8.201789e-05 
    ## ... Procrustes: rmse 4.075144e-05  max resid 8.100668e-05 
    ## ... Similar to previous best
    ## Run 320 stress 9.866524e-05 
    ## ... Procrustes: rmse 0.0002552328  max resid 0.0004361946 
    ## ... Similar to previous best
    ## Run 321 stress 9.807578e-05 
    ## ... Procrustes: rmse 0.0002540553  max resid 0.0004313886 
    ## ... Similar to previous best
    ## Run 322 stress 9.641064e-05 
    ## ... Procrustes: rmse 0.0001646065  max resid 0.0002772786 
    ## ... Similar to previous best
    ## Run 323 stress 9.851034e-05 
    ## ... Procrustes: rmse 0.0002524912  max resid 0.0004294292 
    ## ... Similar to previous best
    ## Run 324 stress 7.432181e-05 
    ## ... Procrustes: rmse 4.234787e-05  max resid 0.0001049762 
    ## ... Similar to previous best
    ## Run 325 stress 9.691308e-05 
    ## ... Procrustes: rmse 0.0002192182  max resid 0.0003773547 
    ## ... Similar to previous best
    ## Run 326 stress 8.951881e-05 
    ## ... Procrustes: rmse 4.378822e-05  max resid 7.262215e-05 
    ## ... Similar to previous best
    ## Run 327 stress 9.846739e-05 
    ## ... Procrustes: rmse 0.000251351  max resid 0.0004328315 
    ## ... Similar to previous best
    ## Run 328 stress 9.851e-05 
    ## ... Procrustes: rmse 0.0002480785  max resid 0.000421187 
    ## ... Similar to previous best
    ## Run 329 stress 9.980983e-05 
    ## ... Procrustes: rmse 0.0002540928  max resid 0.0004370537 
    ## ... Similar to previous best
    ## Run 330 stress 9.879475e-05 
    ## ... Procrustes: rmse 0.0002487855  max resid 0.0004060209 
    ## ... Similar to previous best
    ## Run 331 stress 9.986376e-05 
    ## ... Procrustes: rmse 0.0002188914  max resid 0.0003748899 
    ## ... Similar to previous best
    ## Run 332 stress 9.784626e-05 
    ## ... Procrustes: rmse 0.0001322419  max resid 0.0002028506 
    ## ... Similar to previous best
    ## Run 333 stress 9.463472e-05 
    ## ... Procrustes: rmse 4.631034e-05  max resid 9.324523e-05 
    ## ... Similar to previous best
    ## Run 334 stress 9.831709e-05 
    ## ... Procrustes: rmse 0.0002229695  max resid 0.0003807288 
    ## ... Similar to previous best
    ## Run 335 stress 9.420128e-05 
    ## ... Procrustes: rmse 0.0001926523  max resid 0.0003389484 
    ## ... Similar to previous best
    ## Run 336 stress 9.471953e-05 
    ## ... Procrustes: rmse 0.0002157017  max resid 0.0003723231 
    ## ... Similar to previous best
    ## Run 337 stress 9.958509e-05 
    ## ... Procrustes: rmse 0.0002571607  max resid 0.0004425419 
    ## ... Similar to previous best
    ## Run 338 stress 9.80645e-05 
    ## ... Procrustes: rmse 0.0002539109  max resid 0.0004326001 
    ## ... Similar to previous best
    ## Run 339 stress 9.907023e-05 
    ## ... Procrustes: rmse 0.0002509415  max resid 0.0004384434 
    ## ... Similar to previous best
    ## Run 340 stress 9.9407e-05 
    ## ... Procrustes: rmse 0.0001882503  max resid 0.0003253861 
    ## ... Similar to previous best
    ## Run 341 stress 9.725262e-05 
    ## ... Procrustes: rmse 0.000221146  max resid 0.0003752063 
    ## ... Similar to previous best
    ## Run 342 stress 9.766493e-05 
    ## ... Procrustes: rmse 0.0001946346  max resid 0.0003388325 
    ## ... Similar to previous best
    ## Run 343 stress 9.659274e-05 
    ## ... Procrustes: rmse 0.0002072689  max resid 0.0003414379 
    ## ... Similar to previous best
    ## Run 344 stress 9.964515e-05 
    ## ... Procrustes: rmse 0.0002581479  max resid 0.000438986 
    ## ... Similar to previous best
    ## Run 345 stress 9.864055e-05 
    ## ... Procrustes: rmse 0.0002569414  max resid 0.0004383489 
    ## ... Similar to previous best
    ## Run 346 stress 9.943426e-05 
    ## ... Procrustes: rmse 0.0002203993  max resid 0.0003718681 
    ## ... Similar to previous best
    ## Run 347 stress 9.76489e-05 
    ## ... Procrustes: rmse 0.0002425041  max resid 0.0004153196 
    ## ... Similar to previous best
    ## Run 348 stress 9.546366e-05 
    ## ... Procrustes: rmse 0.0002078129  max resid 0.0003496951 
    ## ... Similar to previous best
    ## Run 349 stress 9.536226e-05 
    ## ... Procrustes: rmse 0.0002078983  max resid 0.000346932 
    ## ... Similar to previous best
    ## Run 350 stress 8.09146e-05 
    ## ... Procrustes: rmse 3.869079e-05  max resid 6.9427e-05 
    ## ... Similar to previous best
    ## Run 351 stress 9.867446e-05 
    ## ... Procrustes: rmse 0.0002540868  max resid 0.0004334936 
    ## ... Similar to previous best
    ## Run 352 stress 9.497756e-05 
    ## ... Procrustes: rmse 0.000223134  max resid 0.00037776 
    ## ... Similar to previous best
    ## Run 353 stress 9.887677e-05 
    ## ... Procrustes: rmse 0.0002564819  max resid 0.0004310008 
    ## ... Similar to previous best
    ## Run 354 stress 9.143881e-05 
    ## ... Procrustes: rmse 5.068217e-05  max resid 9.149547e-05 
    ## ... Similar to previous best
    ## Run 355 stress 9.821959e-05 
    ## ... Procrustes: rmse 0.0002167641  max resid 0.0003710189 
    ## ... Similar to previous best
    ## Run 356 stress 9.044546e-05 
    ## ... Procrustes: rmse 4.650629e-05  max resid 7.338069e-05 
    ## ... Similar to previous best
    ## Run 357 stress 9.806164e-05 
    ## ... Procrustes: rmse 0.0002089378  max resid 0.0003537203 
    ## ... Similar to previous best
    ## Run 358 stress 9.735474e-05 
    ## ... Procrustes: rmse 0.0002205379  max resid 0.0003722337 
    ## ... Similar to previous best
    ## Run 359 stress 8.941136e-05 
    ## ... Procrustes: rmse 4.623575e-05  max resid 7.59901e-05 
    ## ... Similar to previous best
    ## Run 360 stress 9.451684e-05 
    ## ... Procrustes: rmse 0.0001677031  max resid 0.0002958336 
    ## ... Similar to previous best
    ## Run 361 stress 9.845896e-05 
    ## ... Procrustes: rmse 0.000189857  max resid 0.0003144141 
    ## ... Similar to previous best
    ## Run 362 stress 9.909154e-05 
    ## ... Procrustes: rmse 0.0002538919  max resid 0.0004342976 
    ## ... Similar to previous best
    ## Run 363 stress 8.473021e-05 
    ## ... Procrustes: rmse 4.334405e-05  max resid 9.565993e-05 
    ## ... Similar to previous best
    ## Run 364 stress 9.940211e-05 
    ## ... Procrustes: rmse 0.0002518656  max resid 0.0004240805 
    ## ... Similar to previous best
    ## Run 365 stress 9.828729e-05 
    ## ... Procrustes: rmse 0.000253123  max resid 0.0004294773 
    ## ... Similar to previous best
    ## Run 366 stress 9.909125e-05 
    ## ... Procrustes: rmse 0.0002435275  max resid 0.000396814 
    ## ... Similar to previous best
    ## Run 367 stress 9.90953e-05 
    ## ... Procrustes: rmse 0.0002559783  max resid 0.0004356966 
    ## ... Similar to previous best
    ## Run 368 stress 9.966744e-05 
    ## ... Procrustes: rmse 0.000210469  max resid 0.0003556311 
    ## ... Similar to previous best
    ## Run 369 stress 9.740207e-05 
    ## ... Procrustes: rmse 6.567567e-05  max resid 0.0001115064 
    ## ... Similar to previous best
    ## Run 370 stress 9.959797e-05 
    ## ... Procrustes: rmse 0.0002213477  max resid 0.0003751419 
    ## ... Similar to previous best
    ## Run 371 stress 8.951768e-05 
    ## ... Procrustes: rmse 4.802794e-05  max resid 0.0001117698 
    ## ... Similar to previous best
    ## Run 372 stress 9.419981e-05 
    ## ... Procrustes: rmse 0.0002055487  max resid 0.0003542855 
    ## ... Similar to previous best
    ## Run 373 stress 9.906182e-05 
    ## ... Procrustes: rmse 0.000218334  max resid 0.0003689966 
    ## ... Similar to previous best
    ## Run 374 stress 9.899464e-05 
    ## ... Procrustes: rmse 4.94812e-05  max resid 0.0001000156 
    ## ... Similar to previous best
    ## Run 375 stress 5.053115e-05 
    ## ... Procrustes: rmse 3.056713e-05  max resid 5.661414e-05 
    ## ... Similar to previous best
    ## Run 376 stress 9.868411e-05 
    ## ... Procrustes: rmse 0.0002571684  max resid 0.0004353447 
    ## ... Similar to previous best
    ## Run 377 stress 9.835684e-05 
    ## ... Procrustes: rmse 0.0002541636  max resid 0.0004294129 
    ## ... Similar to previous best
    ## Run 378 stress 9.933363e-05 
    ## ... Procrustes: rmse 0.0002289211  max resid 0.0003941127 
    ## ... Similar to previous best
    ## Run 379 stress 9.884343e-05 
    ## ... Procrustes: rmse 0.0002571984  max resid 0.0004342567 
    ## ... Similar to previous best
    ## Run 380 stress 3.30059e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 2.210728e-05  max resid 4.235463e-05 
    ## ... Similar to previous best
    ## Run 381 stress 7.933303e-05 
    ## ... Procrustes: rmse 4.388772e-05  max resid 6.815215e-05 
    ## ... Similar to previous best
    ## Run 382 stress 6.976967e-05 
    ## ... Procrustes: rmse 6.909375e-05  max resid 0.0001619776 
    ## ... Similar to previous best
    ## Run 383 stress 9.950216e-05 
    ## ... Procrustes: rmse 0.0002484358  max resid 0.0004187828 
    ## ... Similar to previous best
    ## Run 384 stress 9.528713e-05 
    ## ... Procrustes: rmse 0.000199673  max resid 0.0003522005 
    ## ... Similar to previous best
    ## Run 385 stress 9.701975e-05 
    ## ... Procrustes: rmse 0.0002106314  max resid 0.0003733293 
    ## ... Similar to previous best
    ## Run 386 stress 4.380033e-05 
    ## ... Procrustes: rmse 2.293988e-05  max resid 4.045573e-05 
    ## ... Similar to previous best
    ## Run 387 stress 9.384356e-05 
    ## ... Procrustes: rmse 3.785829e-05  max resid 7.621964e-05 
    ## ... Similar to previous best
    ## Run 388 stress 9.719377e-05 
    ## ... Procrustes: rmse 0.0002482846  max resid 0.0004178277 
    ## ... Similar to previous best
    ## Run 389 stress 9.876599e-05 
    ## ... Procrustes: rmse 0.0002202701  max resid 0.000374823 
    ## ... Similar to previous best
    ## Run 390 stress 9.699022e-05 
    ## ... Procrustes: rmse 0.0001708392  max resid 0.0003092142 
    ## ... Similar to previous best
    ## Run 391 stress 9.915586e-05 
    ## ... Procrustes: rmse 0.0002520325  max resid 0.000421893 
    ## ... Similar to previous best
    ## Run 392 stress 9.910132e-05 
    ## ... Procrustes: rmse 0.0002250451  max resid 0.0003795307 
    ## ... Similar to previous best
    ## Run 393 stress 4.731195e-05 
    ## ... Procrustes: rmse 2.649783e-05  max resid 4.933496e-05 
    ## ... Similar to previous best
    ## Run 394 stress 6.673589e-05 
    ## ... Procrustes: rmse 8.263028e-05  max resid 0.0001528093 
    ## ... Similar to previous best
    ## Run 395 stress 9.524257e-05 
    ## ... Procrustes: rmse 0.0002131915  max resid 0.0003633513 
    ## ... Similar to previous best
    ## Run 396 stress 9.765884e-05 
    ## ... Procrustes: rmse 0.0001989548  max resid 0.0003557656 
    ## ... Similar to previous best
    ## Run 397 stress 9.879186e-05 
    ## ... Procrustes: rmse 0.0002109133  max resid 0.0003733548 
    ## ... Similar to previous best
    ## Run 398 stress 9.119597e-05 
    ## ... Procrustes: rmse 3.951922e-05  max resid 7.005671e-05 
    ## ... Similar to previous best
    ## Run 399 stress 8.511259e-05 
    ## ... Procrustes: rmse 4.509258e-05  max resid 6.846099e-05 
    ## ... Similar to previous best
    ## Run 400 stress 9.353152e-05 
    ## ... Procrustes: rmse 0.0002059049  max resid 0.0003589351 
    ## ... Similar to previous best
    ## Run 401 stress 9.930348e-05 
    ## ... Procrustes: rmse 0.0001766267  max resid 0.0003222954 
    ## ... Similar to previous best
    ## Run 402 stress 6.770794e-05 
    ## ... Procrustes: rmse 3.548146e-05  max resid 6.745979e-05 
    ## ... Similar to previous best
    ## Run 403 stress 8.845644e-05 
    ## ... Procrustes: rmse 0.000134464  max resid 0.000211862 
    ## ... Similar to previous best
    ## Run 404 stress 9.682153e-05 
    ## ... Procrustes: rmse 0.0001674527  max resid 0.000300492 
    ## ... Similar to previous best
    ## Run 405 stress 9.930458e-05 
    ## ... Procrustes: rmse 0.0002202868  max resid 0.000386011 
    ## ... Similar to previous best
    ## Run 406 stress 8.69955e-05 
    ## ... Procrustes: rmse 4.945562e-05  max resid 7.092879e-05 
    ## ... Similar to previous best
    ## Run 407 stress 9.881348e-05 
    ## ... Procrustes: rmse 0.0002139142  max resid 0.0003743156 
    ## ... Similar to previous best
    ## Run 408 stress 9.775001e-05 
    ## ... Procrustes: rmse 0.00022098  max resid 0.0003854713 
    ## ... Similar to previous best
    ## Run 409 stress 7.961084e-05 
    ## ... Procrustes: rmse 0.000104745  max resid 0.0001765404 
    ## ... Similar to previous best
    ## Run 410 stress 6.731274e-05 
    ## ... Procrustes: rmse 3.864864e-05  max resid 6.422554e-05 
    ## ... Similar to previous best
    ## Run 411 stress 9.182113e-05 
    ## ... Procrustes: rmse 4.445241e-05  max resid 0.0001084313 
    ## ... Similar to previous best
    ## Run 412 stress 9.9622e-05 
    ## ... Procrustes: rmse 0.0002186243  max resid 0.0003812622 
    ## ... Similar to previous best
    ## Run 413 stress 9.79385e-05 
    ## ... Procrustes: rmse 0.0002425741  max resid 0.0004154877 
    ## ... Similar to previous best
    ## Run 414 stress 9.702396e-05 
    ## ... Procrustes: rmse 0.0001854776  max resid 0.0003336215 
    ## ... Similar to previous best
    ## Run 415 stress 9.848143e-05 
    ## ... Procrustes: rmse 0.0002451051  max resid 0.0004298686 
    ## ... Similar to previous best
    ## Run 416 stress 6.854695e-05 
    ## ... Procrustes: rmse 3.541611e-05  max resid 5.423219e-05 
    ## ... Similar to previous best
    ## Run 417 stress 9.785601e-05 
    ## ... Procrustes: rmse 0.0001901605  max resid 0.0003355823 
    ## ... Similar to previous best
    ## Run 418 stress 8.137991e-05 
    ## ... Procrustes: rmse 6.222102e-05  max resid 9.676244e-05 
    ## ... Similar to previous best
    ## Run 419 stress 9.904574e-05 
    ## ... Procrustes: rmse 0.0002519398  max resid 0.0004299315 
    ## ... Similar to previous best
    ## Run 420 stress 6.761438e-05 
    ## ... Procrustes: rmse 3.028435e-05  max resid 6.202659e-05 
    ## ... Similar to previous best
    ## Run 421 stress 9.90946e-05 
    ## ... Procrustes: rmse 0.0002106268  max resid 0.0003749484 
    ## ... Similar to previous best
    ## Run 422 stress 9.86499e-05 
    ## ... Procrustes: rmse 0.0002496184  max resid 0.0004363106 
    ## ... Similar to previous best
    ## Run 423 stress 9.803781e-05 
    ## ... Procrustes: rmse 0.0002024074  max resid 0.0003576914 
    ## ... Similar to previous best
    ## Run 424 stress 7.295741e-05 
    ## ... Procrustes: rmse 4.091233e-05  max resid 9.145651e-05 
    ## ... Similar to previous best
    ## Run 425 stress 9.824888e-05 
    ## ... Procrustes: rmse 0.0002495828  max resid 0.0004154531 
    ## ... Similar to previous best
    ## Run 426 stress 9.933334e-05 
    ## ... Procrustes: rmse 0.0002518891  max resid 0.0004395472 
    ## ... Similar to previous best
    ## Run 427 stress 9.815725e-05 
    ## ... Procrustes: rmse 0.0002478836  max resid 0.0004145856 
    ## ... Similar to previous best
    ## Run 428 stress 6.100118e-05 
    ## ... Procrustes: rmse 3.847826e-05  max resid 7.183669e-05 
    ## ... Similar to previous best
    ## Run 429 stress 7.721863e-05 
    ## ... Procrustes: rmse 0.0001531908  max resid 0.0002841943 
    ## ... Similar to previous best
    ## Run 430 stress 6.543785e-05 
    ## ... Procrustes: rmse 4.952063e-05  max resid 9.013044e-05 
    ## ... Similar to previous best
    ## Run 431 stress 6.491756e-05 
    ## ... Procrustes: rmse 7.677543e-05  max resid 0.0001446214 
    ## ... Similar to previous best
    ## Run 432 stress 9.847948e-05 
    ## ... Procrustes: rmse 0.0002413614  max resid 0.0004254618 
    ## ... Similar to previous best
    ## Run 433 stress 7.692742e-05 
    ## ... Procrustes: rmse 0.0001264944  max resid 0.0002541868 
    ## ... Similar to previous best
    ## Run 434 stress 9.912616e-05 
    ## ... Procrustes: rmse 4.34326e-05  max resid 8.729153e-05 
    ## ... Similar to previous best
    ## Run 435 stress 9.377618e-05 
    ## ... Procrustes: rmse 0.0002317346  max resid 0.0004079981 
    ## ... Similar to previous best
    ## Run 436 stress 9.587588e-05 
    ## ... Procrustes: rmse 5.207305e-05  max resid 9.959508e-05 
    ## ... Similar to previous best
    ## Run 437 stress 9.978598e-05 
    ## ... Procrustes: rmse 0.0001960183  max resid 0.0003385594 
    ## ... Similar to previous best
    ## Run 438 stress 8.613356e-05 
    ## ... Procrustes: rmse 3.973666e-05  max resid 6.745622e-05 
    ## ... Similar to previous best
    ## Run 439 stress 9.790459e-05 
    ## ... Procrustes: rmse 0.0001846626  max resid 0.0003373171 
    ## ... Similar to previous best
    ## Run 440 stress 6.801879e-05 
    ## ... Procrustes: rmse 3.917158e-05  max resid 7.315979e-05 
    ## ... Similar to previous best
    ## Run 441 stress 9.877395e-05 
    ## ... Procrustes: rmse 0.0002079166  max resid 0.0003637791 
    ## ... Similar to previous best
    ## Run 442 stress 6.339864e-05 
    ## ... Procrustes: rmse 3.373781e-05  max resid 5.750557e-05 
    ## ... Similar to previous best
    ## Run 443 stress 9.855401e-05 
    ## ... Procrustes: rmse 0.0002461481  max resid 0.0004314124 
    ## ... Similar to previous best
    ## Run 444 stress 9.893473e-05 
    ## ... Procrustes: rmse 0.0002524169  max resid 0.0004214784 
    ## ... Similar to previous best
    ## Run 445 stress 9.692337e-05 
    ## ... Procrustes: rmse 0.0002148867  max resid 0.0003796365 
    ## ... Similar to previous best
    ## Run 446 stress 9.903315e-05 
    ## ... Procrustes: rmse 0.0002090207  max resid 0.0003659389 
    ## ... Similar to previous best
    ## Run 447 stress 9.843963e-05 
    ## ... Procrustes: rmse 0.0002086794  max resid 0.000366731 
    ## ... Similar to previous best
    ## Run 448 stress 9.972732e-05 
    ## ... Procrustes: rmse 0.0002119488  max resid 0.0003747018 
    ## ... Similar to previous best
    ## Run 449 stress 8.733797e-05 
    ## ... Procrustes: rmse 0.0001628818  max resid 0.0002962644 
    ## ... Similar to previous best
    ## Run 450 stress 9.840114e-05 
    ## ... Procrustes: rmse 0.0002183667  max resid 0.0003759689 
    ## ... Similar to previous best
    ## Run 451 stress 9.284316e-05 
    ## ... Procrustes: rmse 0.0001360837  max resid 0.0002210308 
    ## ... Similar to previous best
    ## Run 452 stress 9.776913e-05 
    ## ... Procrustes: rmse 0.0002058627  max resid 0.0003635684 
    ## ... Similar to previous best
    ## Run 453 stress 9.583414e-05 
    ## ... Procrustes: rmse 0.0002081167  max resid 0.0003659159 
    ## ... Similar to previous best
    ## Run 454 stress 9.823932e-05 
    ## ... Procrustes: rmse 0.0002508815  max resid 0.0004264654 
    ## ... Similar to previous best
    ## Run 455 stress 6.932564e-05 
    ## ... Procrustes: rmse 3.734825e-05  max resid 5.813175e-05 
    ## ... Similar to previous best
    ## Run 456 stress 9.935014e-05 
    ## ... Procrustes: rmse 0.0002544941  max resid 0.0004301366 
    ## ... Similar to previous best
    ## Run 457 stress 9.927116e-05 
    ## ... Procrustes: rmse 0.0002460764  max resid 0.000418634 
    ## ... Similar to previous best
    ## Run 458 stress 5.359726e-05 
    ## ... Procrustes: rmse 6.863343e-05  max resid 0.0001545371 
    ## ... Similar to previous best
    ## Run 459 stress 9.568842e-05 
    ## ... Procrustes: rmse 0.0002163068  max resid 0.0003737911 
    ## ... Similar to previous best
    ## Run 460 stress 9.827042e-05 
    ## ... Procrustes: rmse 0.0002171173  max resid 0.0003798083 
    ## ... Similar to previous best
    ## Run 461 stress 9.709013e-05 
    ## ... Procrustes: rmse 0.0002485929  max resid 0.0004186312 
    ## ... Similar to previous best
    ## Run 462 stress 6.635358e-05 
    ## ... Procrustes: rmse 3.064255e-05  max resid 6.073809e-05 
    ## ... Similar to previous best
    ## Run 463 stress 9.582343e-05 
    ## ... Procrustes: rmse 0.000211662  max resid 0.0003665758 
    ## ... Similar to previous best
    ## Run 464 stress 6.620287e-05 
    ## ... Procrustes: rmse 3.873723e-05  max resid 6.258727e-05 
    ## ... Similar to previous best
    ## Run 465 stress 9.550409e-05 
    ## ... Procrustes: rmse 0.0002414942  max resid 0.0004221288 
    ## ... Similar to previous best
    ## Run 466 stress 9.632556e-05 
    ## ... Procrustes: rmse 0.0001742729  max resid 0.0003131834 
    ## ... Similar to previous best
    ## Run 467 stress 9.787548e-05 
    ## ... Procrustes: rmse 0.0002086232  max resid 0.0003639146 
    ## ... Similar to previous best
    ## Run 468 stress 8.831348e-05 
    ## ... Procrustes: rmse 5.385639e-05  max resid 8.164704e-05 
    ## ... Similar to previous best
    ## Run 469 stress 9.097934e-05 
    ## ... Procrustes: rmse 0.0002223119  max resid 0.0003763421 
    ## ... Similar to previous best
    ## Run 470 stress 9.429047e-05 
    ## ... Procrustes: rmse 5.256672e-05  max resid 0.0001154439 
    ## ... Similar to previous best
    ## Run 471 stress 0.377012 
    ## Run 472 stress 8.002671e-05 
    ## ... Procrustes: rmse 3.841481e-05  max resid 6.908859e-05 
    ## ... Similar to previous best
    ## Run 473 stress 9.859665e-05 
    ## ... Procrustes: rmse 0.0002475924  max resid 0.0004331949 
    ## ... Similar to previous best
    ## Run 474 stress 9.83012e-05 
    ## ... Procrustes: rmse 0.0002194695  max resid 0.0003751706 
    ## ... Similar to previous best
    ## Run 475 stress 9.847605e-05 
    ## ... Procrustes: rmse 0.0002455535  max resid 0.0004262978 
    ## ... Similar to previous best
    ## Run 476 stress 9.830806e-05 
    ## ... Procrustes: rmse 0.0002152714  max resid 0.0003716743 
    ## ... Similar to previous best
    ## Run 477 stress 9.843979e-05 
    ## ... Procrustes: rmse 0.0002457604  max resid 0.0004272323 
    ## ... Similar to previous best
    ## Run 478 stress 9.893672e-05 
    ## ... Procrustes: rmse 0.0001018879  max resid 0.0001812853 
    ## ... Similar to previous best
    ## Run 479 stress 9.36746e-05 
    ## ... Procrustes: rmse 0.0001855557  max resid 0.0003304927 
    ## ... Similar to previous best
    ## Run 480 stress 9.824787e-05 
    ## ... Procrustes: rmse 0.0002504726  max resid 0.0004224416 
    ## ... Similar to previous best
    ## Run 481 stress 9.717324e-05 
    ## ... Procrustes: rmse 0.000248556  max resid 0.0004230236 
    ## ... Similar to previous best
    ## Run 482 stress 8.303726e-05 
    ## ... Procrustes: rmse 4.67116e-05  max resid 7.810636e-05 
    ## ... Similar to previous best
    ## Run 483 stress 9.430396e-05 
    ## ... Procrustes: rmse 6.540057e-05  max resid 0.0001369664 
    ## ... Similar to previous best
    ## Run 484 stress 9.762199e-05 
    ## ... Procrustes: rmse 0.0001463212  max resid 0.0002570212 
    ## ... Similar to previous best
    ## Run 485 stress 9.953559e-05 
    ## ... Procrustes: rmse 0.0002067263  max resid 0.0003673995 
    ## ... Similar to previous best
    ## Run 486 stress 9.730441e-05 
    ## ... Procrustes: rmse 0.0002440855  max resid 0.0004284646 
    ## ... Similar to previous best
    ## Run 487 stress 9.891969e-05 
    ## ... Procrustes: rmse 0.0002132849  max resid 0.0003716001 
    ## ... Similar to previous best
    ## Run 488 stress 9.837281e-05 
    ## ... Procrustes: rmse 0.0002495656  max resid 0.0004288718 
    ## ... Similar to previous best
    ## Run 489 stress 9.869393e-05 
    ## ... Procrustes: rmse 0.0001880018  max resid 0.000337852 
    ## ... Similar to previous best
    ## Run 490 stress 7.166478e-05 
    ## ... Procrustes: rmse 4.052657e-05  max resid 6.124181e-05 
    ## ... Similar to previous best
    ## Run 491 stress 9.933572e-05 
    ## ... Procrustes: rmse 0.0002462439  max resid 0.0004327672 
    ## ... Similar to previous best
    ## Run 492 stress 9.856933e-05 
    ## ... Procrustes: rmse 0.0002481335  max resid 0.0004227204 
    ## ... Similar to previous best
    ## Run 493 stress 8.587433e-05 
    ## ... Procrustes: rmse 8.664368e-05  max resid 0.0001806218 
    ## ... Similar to previous best
    ## Run 494 stress 9.21079e-05 
    ## ... Procrustes: rmse 4.670612e-05  max resid 9.384548e-05 
    ## ... Similar to previous best
    ## Run 495 stress 9.805764e-05 
    ## ... Procrustes: rmse 0.0002176599  max resid 0.000382884 
    ## ... Similar to previous best
    ## Run 496 stress 9.737166e-05 
    ## ... Procrustes: rmse 0.0002208575  max resid 0.0003769522 
    ## ... Similar to previous best
    ## Run 497 stress 9.652697e-05 
    ## ... Procrustes: rmse 0.0002077115  max resid 0.0003677362 
    ## ... Similar to previous best
    ## Run 498 stress 6.097251e-05 
    ## ... Procrustes: rmse 3.329012e-05  max resid 4.637419e-05 
    ## ... Similar to previous best
    ## Run 499 stress 9.832023e-05 
    ## ... Procrustes: rmse 5.309117e-05  max resid 0.0001033956 
    ## ... Similar to previous best
    ## Run 500 stress 9.697961e-05 
    ## ... Procrustes: rmse 0.0001808179  max resid 0.000327178 
    ## ... Similar to previous best
    ## *** Best solution repeated 120 times

    ## Warning in metaMDS(PD_beta_ref_dist$Btotal, try = 1000, parallel = 4, trymax =
    ## 500, : stress is (nearly) zero: you may have insufficient data

``` r
# Surveyed sites 
PD_beta_NMDS <- metaMDS(PD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07927535 
    ## Run 1 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734411  max resid 0.0564073 
    ## Run 2 stress 0.07927535 
    ## ... Procrustes: rmse 6.886207e-05  max resid 0.0001864682 
    ## ... Similar to previous best
    ## Run 3 stress 0.0793695 
    ## ... Procrustes: rmse 0.009038379  max resid 0.03242695 
    ## Run 4 stress 0.09252216 
    ## Run 5 stress 0.09252214 
    ## Run 6 stress 0.09388354 
    ## Run 7 stress 0.09072901 
    ## Run 8 stress 0.104507 
    ## Run 9 stress 0.0796926 
    ## ... Procrustes: rmse 0.01236603  max resid 0.04966287 
    ## Run 10 stress 0.08985605 
    ## Run 11 stress 0.09106041 
    ## Run 12 stress 0.07969266 
    ## ... Procrustes: rmse 0.01236182  max resid 0.04957102 
    ## Run 13 stress 0.0907291 
    ## Run 14 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237206  max resid 0.04975039 
    ## Run 15 stress 0.09485686 
    ## Run 16 stress 0.09247106 
    ## Run 17 stress 0.09544387 
    ## Run 18 stress 0.08985604 
    ## Run 19 stress 0.09043765 
    ## Run 20 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735202  max resid 0.05646367 
    ## Run 21 stress 0.09247112 
    ## Run 22 stress 0.1077585 
    ## Run 23 stress 0.07951445 
    ## ... Procrustes: rmse 0.0173177  max resid 0.05622108 
    ## Run 24 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735002  max resid 0.05644444 
    ## Run 25 stress 0.08985608 
    ## Run 26 stress 0.0898561 
    ## Run 27 stress 0.09106049 
    ## Run 28 stress 0.07969259 
    ## ... Procrustes: rmse 0.01236044  max resid 0.0496572 
    ## Run 29 stress 0.07951445 
    ## ... Procrustes: rmse 0.01734397  max resid 0.05640078 
    ## Run 30 stress 0.09344086 
    ## Run 31 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 4.203221e-05  max resid 0.0001108871 
    ## ... Similar to previous best
    ## Run 32 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 1.056191e-05  max resid 2.601851e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.09565378 
    ## Run 34 stress 0.08985604 
    ## Run 35 stress 0.08985612 
    ## Run 36 stress 0.07936951 
    ## ... Procrustes: rmse 0.009035393  max resid 0.03240219 
    ## Run 37 stress 0.07936951 
    ## ... Procrustes: rmse 0.009037495  max resid 0.03241634 
    ## Run 38 stress 0.07969271 
    ## ... Procrustes: rmse 0.0124286  max resid 0.05003604 
    ## Run 39 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733776  max resid 0.05631232 
    ## Run 40 stress 0.09022776 
    ## Run 41 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044983  max resid 0.03245997 
    ## Run 42 stress 0.09344085 
    ## Run 43 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238691  max resid 0.04975732 
    ## Run 44 stress 0.07969263 
    ## ... Procrustes: rmse 0.01241467  max resid 0.04993384 
    ## Run 45 stress 0.07927535 
    ## ... Procrustes: rmse 2.361806e-05  max resid 6.209458e-05 
    ## ... Similar to previous best
    ## Run 46 stress 0.08985607 
    ## Run 47 stress 0.07936951 
    ## ... Procrustes: rmse 0.009039544  max resid 0.03242826 
    ## Run 48 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735012  max resid 0.05639519 
    ## Run 49 stress 0.09072901 
    ## Run 50 stress 0.0910603 
    ## Run 51 stress 0.07927534 
    ## ... Procrustes: rmse 2.017172e-05  max resid 5.586869e-05 
    ## ... Similar to previous best
    ## Run 52 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043671  max resid 0.03245605 
    ## Run 53 stress 0.07969259 
    ## ... Procrustes: rmse 0.0123716  max resid 0.04969412 
    ## Run 54 stress 0.07927534 
    ## ... Procrustes: rmse 7.612227e-06  max resid 1.955614e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.09043757 
    ## Run 56 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735176  max resid 0.05641288 
    ## Run 57 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046872  max resid 0.03247551 
    ## Run 58 stress 0.07927536 
    ## ... Procrustes: rmse 4.617404e-05  max resid 0.0001271442 
    ## ... Similar to previous best
    ## Run 59 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734701  max resid 0.05637653 
    ## Run 60 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044301  max resid 0.03245722 
    ## Run 61 stress 0.07936951 
    ## ... Procrustes: rmse 0.00903515  max resid 0.0324011 
    ## Run 62 stress 0.08985605 
    ## Run 63 stress 0.07927535 
    ## ... Procrustes: rmse 2.255041e-05  max resid 5.699109e-05 
    ## ... Similar to previous best
    ## Run 64 stress 0.09072901 
    ## Run 65 stress 0.08985606 
    ## Run 66 stress 0.09043753 
    ## Run 67 stress 0.09043738 
    ## Run 68 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173641  max resid 0.05648552 
    ## Run 69 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045427  max resid 0.03247554 
    ## Run 70 stress 0.07927537 
    ## ... Procrustes: rmse 5.643572e-05  max resid 0.00015268 
    ## ... Similar to previous best
    ## Run 71 stress 0.08985604 
    ## Run 72 stress 0.07927534 
    ## ... Procrustes: rmse 1.915475e-05  max resid 5.250204e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.09474053 
    ## Run 74 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734428  max resid 0.05635857 
    ## Run 75 stress 0.07927537 
    ## ... Procrustes: rmse 6.511259e-05  max resid 0.0001545326 
    ## ... Similar to previous best
    ## Run 76 stress 0.07936956 
    ## ... Procrustes: rmse 0.009024575  max resid 0.03233593 
    ## Run 77 stress 0.07927536 
    ## ... Procrustes: rmse 1.98667e-05  max resid 3.446746e-05 
    ## ... Similar to previous best
    ## Run 78 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735448  max resid 0.05644095 
    ## Run 79 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 3.558866e-06  max resid 1.094296e-05 
    ## ... Similar to previous best
    ## Run 80 stress 0.0935632 
    ## Run 81 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733548  max resid 0.05630624 
    ## Run 82 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041074  max resid 0.03243935 
    ## Run 83 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237635  max resid 0.04974017 
    ## Run 84 stress 0.07927536 
    ## ... Procrustes: rmse 4.961562e-05  max resid 0.0001375319 
    ## ... Similar to previous best
    ## Run 85 stress 0.09247101 
    ## Run 86 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238513  max resid 0.04975937 
    ## Run 87 stress 0.07936952 
    ## ... Procrustes: rmse 0.009039277  max resid 0.03244105 
    ## Run 88 stress 0.08985606 
    ## Run 89 stress 0.09290699 
    ## Run 90 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734894  max resid 0.0563933 
    ## Run 91 stress 0.09356322 
    ## Run 92 stress 0.08985615 
    ## Run 93 stress 0.07927534 
    ## ... Procrustes: rmse 2.805028e-06  max resid 6.793097e-06 
    ## ... Similar to previous best
    ## Run 94 stress 0.09106034 
    ## Run 95 stress 0.08985625 
    ## Run 96 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239016  max resid 0.04979105 
    ## Run 97 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046676  max resid 0.03247812 
    ## Run 98 stress 0.09485688 
    ## Run 99 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735681  max resid 0.05648701 
    ## Run 100 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735894  max resid 0.05648087 
    ## Run 101 stress 0.09022775 
    ## Run 102 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735317  max resid 0.05643978 
    ## Run 103 stress 0.09388326 
    ## Run 104 stress 0.09290691 
    ## Run 105 stress 0.1088749 
    ## Run 106 stress 0.07927534 
    ## ... Procrustes: rmse 1.399979e-05  max resid 3.222634e-05 
    ## ... Similar to previous best
    ## Run 107 stress 0.09022794 
    ## Run 108 stress 0.09296243 
    ## Run 109 stress 0.09043735 
    ## Run 110 stress 0.09106045 
    ## Run 111 stress 0.09072893 
    ## Run 112 stress 0.09072897 
    ## Run 113 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904447  max resid 0.03246297 
    ## Run 114 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735018  max resid 0.05645888 
    ## Run 115 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904145  max resid 0.03244197 
    ## Run 116 stress 0.08985609 
    ## Run 117 stress 0.09290709 
    ## Run 118 stress 0.07951444 
    ## ... Procrustes: rmse 0.01737515  max resid 0.05655844 
    ## Run 119 stress 0.09296252 
    ## Run 120 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734959  max resid 0.05638734 
    ## Run 121 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735953  max resid 0.05645862 
    ## Run 122 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734908  max resid 0.05638549 
    ## Run 123 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238815  max resid 0.04975004 
    ## Run 124 stress 0.09072892 
    ## Run 125 stress 0.07927534 
    ## ... Procrustes: rmse 1.066773e-05  max resid 2.777598e-05 
    ## ... Similar to previous best
    ## Run 126 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733707  max resid 0.05632161 
    ## Run 127 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735572  max resid 0.05646816 
    ## Run 128 stress 0.09106029 
    ## Run 129 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238214  max resid 0.04974664 
    ## Run 130 stress 0.09022775 
    ## Run 131 stress 0.07951446 
    ## ... Procrustes: rmse 0.01738721  max resid 0.05663935 
    ## Run 132 stress 0.07927536 
    ## ... Procrustes: rmse 5.217195e-05  max resid 0.0001235044 
    ## ... Similar to previous best
    ## Run 133 stress 0.08985608 
    ## Run 134 stress 0.0929069 
    ## Run 135 stress 0.08985605 
    ## Run 136 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736109  max resid 0.05647393 
    ## Run 137 stress 0.09544374 
    ## Run 138 stress 0.09290691 
    ## Run 139 stress 0.07927534 
    ## ... Procrustes: rmse 5.905926e-06  max resid 2.223606e-05 
    ## ... Similar to previous best
    ## Run 140 stress 0.1061856 
    ## Run 141 stress 0.07927534 
    ## ... Procrustes: rmse 2.145543e-05  max resid 4.820631e-05 
    ## ... Similar to previous best
    ## Run 142 stress 0.07927534 
    ## ... Procrustes: rmse 6.557472e-06  max resid 1.750738e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238107  max resid 0.04975052 
    ## Run 144 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239773  max resid 0.049831 
    ## Run 145 stress 0.09043746 
    ## Run 146 stress 0.09072907 
    ## Run 147 stress 0.08985604 
    ## Run 148 stress 0.07927534 
    ## ... Procrustes: rmse 2.001448e-05  max resid 5.333865e-05 
    ## ... Similar to previous best
    ## Run 149 stress 0.0948573 
    ## Run 150 stress 0.07927534 
    ## ... Procrustes: rmse 1.323335e-06  max resid 4.121582e-06 
    ## ... Similar to previous best
    ## Run 151 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238965  max resid 0.04978287 
    ## Run 152 stress 0.07951444 
    ## ... Procrustes: rmse 0.01734892  max resid 0.05639089 
    ## Run 153 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734994  max resid 0.05640598 
    ## Run 154 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735302  max resid 0.05640665 
    ## Run 155 stress 0.07951445 
    ## ... Procrustes: rmse 0.01734884  max resid 0.05640025 
    ## Run 156 stress 0.09565386 
    ## Run 157 stress 0.07951447 
    ## ... Procrustes: rmse 0.01739288  max resid 0.0566844 
    ## Run 158 stress 0.09072892 
    ## Run 159 stress 0.08985606 
    ## Run 160 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735347  max resid 0.05642201 
    ## Run 161 stress 0.07927535 
    ## ... Procrustes: rmse 1.377129e-05  max resid 3.506731e-05 
    ## ... Similar to previous best
    ## Run 162 stress 0.07927541 
    ## ... Procrustes: rmse 8.475911e-05  max resid 0.0002285606 
    ## ... Similar to previous best
    ## Run 163 stress 0.09290715 
    ## Run 164 stress 0.09296276 
    ## Run 165 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735575  max resid 0.05643607 
    ## Run 166 stress 0.07951444 
    ## ... Procrustes: rmse 0.01735895  max resid 0.05653866 
    ## Run 167 stress 0.07951445 
    ## ... Procrustes: rmse 0.01736314  max resid 0.05656569 
    ## Run 168 stress 0.07927534 
    ## ... Procrustes: rmse 9.489327e-06  max resid 2.422153e-05 
    ## ... Similar to previous best
    ## Run 169 stress 0.09252213 
    ## Run 170 stress 0.07927537 
    ## ... Procrustes: rmse 2.564236e-05  max resid 6.754901e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.07969259 
    ## ... Procrustes: rmse 0.0123871  max resid 0.04977188 
    ## Run 172 stress 0.07951443 
    ## ... Procrustes: rmse 0.01737447  max resid 0.05655513 
    ## Run 173 stress 0.09043763 
    ## Run 174 stress 0.09247102 
    ## Run 175 stress 0.08985603 
    ## Run 176 stress 0.0793695 
    ## ... Procrustes: rmse 0.009038915  max resid 0.03242649 
    ## Run 177 stress 0.09485708 
    ## Run 178 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237089  max resid 0.04961614 
    ## Run 179 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734849  max resid 0.05637449 
    ## Run 180 stress 0.0925222 
    ## Run 181 stress 0.07927534 
    ## ... Procrustes: rmse 2.909675e-06  max resid 6.552352e-06 
    ## ... Similar to previous best
    ## Run 182 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735711  max resid 0.05644325 
    ## Run 183 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734755  max resid 0.05638294 
    ## Run 184 stress 0.08985616 
    ## Run 185 stress 0.07936951 
    ## ... Procrustes: rmse 0.009037367  max resid 0.03241363 
    ## Run 186 stress 0.09565379 
    ## Run 187 stress 0.07969268 
    ## ... Procrustes: rmse 0.01236569  max resid 0.04953587 
    ## Run 188 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733351  max resid 0.05628933 
    ## Run 189 stress 0.09106046 
    ## Run 190 stress 0.0793695 
    ## ... Procrustes: rmse 0.009037004  max resid 0.0324132 
    ## Run 191 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735444  max resid 0.05649328 
    ## Run 192 stress 0.09388363 
    ## Run 193 stress 0.09043737 
    ## Run 194 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237609  max resid 0.04971904 
    ## Run 195 stress 0.07969258 
    ## ... Procrustes: rmse 0.0123752  max resid 0.04969343 
    ## Run 196 stress 0.07951445 
    ## ... Procrustes: rmse 0.0173285  max resid 0.0562334 
    ## Run 197 stress 0.09356313 
    ## Run 198 stress 0.07927536 
    ## ... Procrustes: rmse 3.56129e-05  max resid 8.897254e-05 
    ## ... Similar to previous best
    ## Run 199 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734785  max resid 0.05638671 
    ## Run 200 stress 0.0929626 
    ## Run 201 stress 0.07969265 
    ## ... Procrustes: rmse 0.01240299  max resid 0.04985675 
    ## Run 202 stress 0.07927534 
    ## ... Procrustes: rmse 6.440895e-06  max resid 1.750314e-05 
    ## ... Similar to previous best
    ## Run 203 stress 0.1033494 
    ## Run 204 stress 0.08985606 
    ## Run 205 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735557  max resid 0.05646945 
    ## Run 206 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238211  max resid 0.04976746 
    ## Run 207 stress 0.09344089 
    ## Run 208 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735411  max resid 0.0564691 
    ## Run 209 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237738  max resid 0.04971392 
    ## Run 210 stress 0.07951442 
    ## ... Procrustes: rmse 0.0173698  max resid 0.0565304 
    ## Run 211 stress 0.07927537 
    ## ... Procrustes: rmse 5.29089e-05  max resid 0.0001469987 
    ## ... Similar to previous best
    ## Run 212 stress 0.0907291 
    ## Run 213 stress 0.0795145 
    ## ... Procrustes: rmse 0.01731635  max resid 0.05615824 
    ## Run 214 stress 0.07936956 
    ## ... Procrustes: rmse 0.009022726  max resid 0.03232492 
    ## Run 215 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735568  max resid 0.05642996 
    ## Run 216 stress 0.09388328 
    ## Run 217 stress 0.07936949 
    ## ... Procrustes: rmse 0.009044922  max resid 0.03246436 
    ## Run 218 stress 0.0938833 
    ## Run 219 stress 0.08985615 
    ## Run 220 stress 0.07927536 
    ## ... Procrustes: rmse 5.440064e-05  max resid 0.0001517828 
    ## ... Similar to previous best
    ## Run 221 stress 0.07927535 
    ## ... Procrustes: rmse 3.487702e-05  max resid 9.369741e-05 
    ## ... Similar to previous best
    ## Run 222 stress 0.1066421 
    ## Run 223 stress 0.09106029 
    ## Run 224 stress 0.09544363 
    ## Run 225 stress 0.09172795 
    ## Run 226 stress 0.09252236 
    ## Run 227 stress 0.09106037 
    ## Run 228 stress 0.07927535 
    ## ... Procrustes: rmse 2.248445e-05  max resid 5.803823e-05 
    ## ... Similar to previous best
    ## Run 229 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039729  max resid 0.03243107 
    ## Run 230 stress 0.07936956 
    ## ... Procrustes: rmse 0.009022975  max resid 0.03232629 
    ## Run 231 stress 0.07969262 
    ## ... Procrustes: rmse 0.01239912  max resid 0.04985622 
    ## Run 232 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237377  max resid 0.0496857 
    ## Run 233 stress 0.09106029 
    ## Run 234 stress 0.07951443 
    ## ... Procrustes: rmse 0.01737142  max resid 0.05654529 
    ## Run 235 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039962  max resid 0.03242986 
    ## Run 236 stress 0.07969259 
    ## ... Procrustes: rmse 0.0123746  max resid 0.04968659 
    ## Run 237 stress 0.09344086 
    ## Run 238 stress 0.07951443 
    ## ... Procrustes: rmse 0.01735131  max resid 0.05641847 
    ## Run 239 stress 0.07927537 
    ## ... Procrustes: rmse 5.705742e-05  max resid 0.0001545842 
    ## ... Similar to previous best
    ## Run 240 stress 0.09106029 
    ## Run 241 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238773  max resid 0.04977315 
    ## Run 242 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239785  max resid 0.04985365 
    ## Run 243 stress 0.07969263 
    ## ... Procrustes: rmse 0.01237754  max resid 0.04974808 
    ## Run 244 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240833  max resid 0.04987025 
    ## Run 245 stress 0.07951445 
    ## ... Procrustes: rmse 0.01736079  max resid 0.05656842 
    ## Run 246 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735074  max resid 0.05639931 
    ## Run 247 stress 0.07951448 
    ## ... Procrustes: rmse 0.01732235  max resid 0.05619413 
    ## Run 248 stress 0.09022776 
    ## Run 249 stress 0.09252227 
    ## Run 250 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239377  max resid 0.04977523 
    ## Run 251 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733642  max resid 0.05630922 
    ## Run 252 stress 0.09247104 
    ## Run 253 stress 0.07936953 
    ## ... Procrustes: rmse 0.009029103  max resid 0.03236459 
    ## Run 254 stress 0.07936951 
    ## ... Procrustes: rmse 0.009046188  max resid 0.03247546 
    ## Run 255 stress 0.09296259 
    ## Run 256 stress 0.07927534 
    ## ... Procrustes: rmse 3.487293e-06  max resid 7.892382e-06 
    ## ... Similar to previous best
    ## Run 257 stress 0.09106036 
    ## Run 258 stress 0.07927534 
    ## ... Procrustes: rmse 1.975423e-05  max resid 5.382668e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.08985613 
    ## Run 260 stress 0.07936952 
    ## ... Procrustes: rmse 0.009031998  max resid 0.03238133 
    ## Run 261 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 7.851152e-07  max resid 1.725764e-06 
    ## ... Similar to previous best
    ## Run 262 stress 0.07969268 
    ## ... Procrustes: rmse 0.01242799  max resid 0.0500217 
    ## Run 263 stress 0.09296227 
    ## Run 264 stress 0.08985603 
    ## Run 265 stress 0.0793695 
    ## ... Procrustes: rmse 0.009048235  max resid 0.03248025 
    ## Run 266 stress 0.09172787 
    ## Run 267 stress 0.08985605 
    ## Run 268 stress 0.09106036 
    ## Run 269 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735738  max resid 0.0564448 
    ## Run 270 stress 0.09544359 
    ## Run 271 stress 0.07969262 
    ## ... Procrustes: rmse 0.01237782  max resid 0.04964553 
    ## Run 272 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237698  max resid 0.04972358 
    ## Run 273 stress 0.08985607 
    ## Run 274 stress 0.09043744 
    ## Run 275 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735077  max resid 0.05639718 
    ## Run 276 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173572  max resid 0.05644582 
    ## Run 277 stress 0.07927534 
    ## ... Procrustes: rmse 7.739149e-06  max resid 2.747355e-05 
    ## ... Similar to previous best
    ## Run 278 stress 0.07927534 
    ## ... Procrustes: rmse 1.279727e-05  max resid 3.066364e-05 
    ## ... Similar to previous best
    ## Run 279 stress 0.09022799 
    ## Run 280 stress 0.106141 
    ## Run 281 stress 0.07927536 
    ## ... Procrustes: rmse 4.907047e-05  max resid 0.0001167026 
    ## ... Similar to previous best
    ## Run 282 stress 0.07927534 
    ## ... Procrustes: rmse 7.752168e-06  max resid 1.38055e-05 
    ## ... Similar to previous best
    ## Run 283 stress 0.09565391 
    ## Run 284 stress 0.08985611 
    ## Run 285 stress 0.07927538 
    ## ... Procrustes: rmse 7.598958e-05  max resid 0.0002133309 
    ## ... Similar to previous best
    ## Run 286 stress 0.07927535 
    ## ... Procrustes: rmse 1.216857e-05  max resid 2.274712e-05 
    ## ... Similar to previous best
    ## Run 287 stress 0.07927534 
    ## ... Procrustes: rmse 5.730894e-06  max resid 1.371904e-05 
    ## ... Similar to previous best
    ## Run 288 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238342  max resid 0.04974685 
    ## Run 289 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734822  max resid 0.05639403 
    ## Run 290 stress 0.08985603 
    ## Run 291 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237209  max resid 0.04964781 
    ## Run 292 stress 0.09485734 
    ## Run 293 stress 0.09474085 
    ## Run 294 stress 0.07927536 
    ## ... Procrustes: rmse 5.175374e-05  max resid 0.0001346564 
    ## ... Similar to previous best
    ## Run 295 stress 0.09106031 
    ## Run 296 stress 0.0925222 
    ## Run 297 stress 0.09388348 
    ## Run 298 stress 0.09106041 
    ## Run 299 stress 0.07951443 
    ## ... Procrustes: rmse 0.01735358  max resid 0.05643487 
    ## Run 300 stress 0.09474065 
    ## Run 301 stress 0.09072904 
    ## Run 302 stress 0.09296206 
    ## Run 303 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238157  max resid 0.04975152 
    ## Run 304 stress 0.09252217 
    ## Run 305 stress 0.08985604 
    ## Run 306 stress 0.0910605 
    ## Run 307 stress 0.08985605 
    ## Run 308 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237268  max resid 0.04966096 
    ## Run 309 stress 0.0792754 
    ## ... Procrustes: rmse 7.816513e-05  max resid 0.0002191593 
    ## ... Similar to previous best
    ## Run 310 stress 0.07927537 
    ## ... Procrustes: rmse 6.810544e-05  max resid 0.000191473 
    ## ... Similar to previous best
    ## Run 311 stress 0.09388358 
    ## Run 312 stress 0.09356315 
    ## Run 313 stress 0.09565375 
    ## Run 314 stress 0.09296235 
    ## Run 315 stress 0.09544372 
    ## Run 316 stress 0.07927534 
    ## ... Procrustes: rmse 7.33902e-06  max resid 1.463161e-05 
    ## ... Similar to previous best
    ## Run 317 stress 0.07927534 
    ## ... Procrustes: rmse 1.376143e-05  max resid 3.636666e-05 
    ## ... Similar to previous best
    ## Run 318 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735684  max resid 0.05645603 
    ## Run 319 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045547  max resid 0.03246968 
    ## Run 320 stress 0.09565362 
    ## Run 321 stress 0.07927535 
    ## ... Procrustes: rmse 3.572863e-05  max resid 9.929664e-05 
    ## ... Similar to previous best
    ## Run 322 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734838  max resid 0.05637951 
    ## Run 323 stress 0.09106054 
    ## Run 324 stress 0.08985609 
    ## Run 325 stress 0.07927537 
    ## ... Procrustes: rmse 6.297645e-05  max resid 0.0001762928 
    ## ... Similar to previous best
    ## Run 326 stress 0.09296216 
    ## Run 327 stress 0.1108125 
    ## Run 328 stress 0.09043757 
    ## Run 329 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238141  max resid 0.04976502 
    ## Run 330 stress 0.1110285 
    ## Run 331 stress 0.07927534 
    ## ... Procrustes: rmse 9.132466e-06  max resid 2.737838e-05 
    ## ... Similar to previous best
    ## Run 332 stress 0.09043777 
    ## Run 333 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734866  max resid 0.05639187 
    ## Run 334 stress 0.07951444 
    ## ... Procrustes: rmse 0.0173255  max resid 0.05626008 
    ## Run 335 stress 0.07927539 
    ## ... Procrustes: rmse 8.359776e-05  max resid 0.0002360337 
    ## ... Similar to previous best
    ## Run 336 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734988  max resid 0.05641485 
    ## Run 337 stress 0.09544386 
    ## Run 338 stress 0.0795144 
    ## ... Procrustes: rmse 0.01733751  max resid 0.05635309 
    ## Run 339 stress 0.07927534 
    ## ... Procrustes: rmse 1.182446e-05  max resid 3.194257e-05 
    ## ... Similar to previous best
    ## Run 340 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049868  max resid 0.03248958 
    ## Run 341 stress 0.08985604 
    ## Run 342 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733967  max resid 0.05632602 
    ## Run 343 stress 0.0796926 
    ## ... Procrustes: rmse 0.01236287  max resid 0.04958277 
    ## Run 344 stress 0.09043766 
    ## Run 345 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173541  max resid 0.05640965 
    ## Run 346 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734413  max resid 0.05635859 
    ## Run 347 stress 0.09106054 
    ## Run 348 stress 0.09072898 
    ## Run 349 stress 0.09485676 
    ## Run 350 stress 0.09106033 
    ## Run 351 stress 0.07927534 
    ## ... Procrustes: rmse 2.092837e-05  max resid 4.752733e-05 
    ## ... Similar to previous best
    ## Run 352 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045213  max resid 0.03246978 
    ## Run 353 stress 0.08985607 
    ## Run 354 stress 0.092907 
    ## Run 355 stress 0.07936956 
    ## ... Procrustes: rmse 0.00902368  max resid 0.03233044 
    ## Run 356 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173514  max resid 0.05640872 
    ## Run 357 stress 0.09106041 
    ## Run 358 stress 0.07927535 
    ## ... Procrustes: rmse 3.001584e-05  max resid 7.482853e-05 
    ## ... Similar to previous best
    ## Run 359 stress 0.07927534 
    ## ... Procrustes: rmse 8.560098e-06  max resid 2.011591e-05 
    ## ... Similar to previous best
    ## Run 360 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044073  max resid 0.03246501 
    ## Run 361 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734723  max resid 0.05638226 
    ## Run 362 stress 0.07936952 
    ## ... Procrustes: rmse 0.009033858  max resid 0.0323942 
    ## Run 363 stress 0.09106043 
    ## Run 364 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734565  max resid 0.05634144 
    ## Run 365 stress 0.09106029 
    ## Run 366 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237441  max resid 0.04969461 
    ## Run 367 stress 0.09072893 
    ## Run 368 stress 0.09296251 
    ## Run 369 stress 0.1047632 
    ## Run 370 stress 0.09296249 
    ## Run 371 stress 0.1074957 
    ## Run 372 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735265  max resid 0.05643677 
    ## Run 373 stress 0.09247105 
    ## Run 374 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735209  max resid 0.05640201 
    ## Run 375 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173606  max resid 0.05649504 
    ## Run 376 stress 0.09544375 
    ## Run 377 stress 0.09106036 
    ## Run 378 stress 0.07936951 
    ## ... Procrustes: rmse 0.009033786  max resid 0.03239414 
    ## Run 379 stress 0.08985604 
    ## Run 380 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239715  max resid 0.04981374 
    ## Run 381 stress 0.09247109 
    ## Run 382 stress 0.09043764 
    ## Run 383 stress 0.1038073 
    ## Run 384 stress 0.07969263 
    ## ... Procrustes: rmse 0.01241029  max resid 0.04991814 
    ## Run 385 stress 0.09072902 
    ## Run 386 stress 0.09072892 
    ## Run 387 stress 0.08985605 
    ## Run 388 stress 0.07927537 
    ## ... Procrustes: rmse 6.508052e-05  max resid 0.0001826911 
    ## ... Similar to previous best
    ## Run 389 stress 0.1033494 
    ## Run 390 stress 0.09252216 
    ## Run 391 stress 0.08985611 
    ## Run 392 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734198  max resid 0.05633737 
    ## Run 393 stress 0.1047992 
    ## Run 394 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734529  max resid 0.05637279 
    ## Run 395 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735325  max resid 0.05642345 
    ## Run 396 stress 0.08985603 
    ## Run 397 stress 0.07969263 
    ## ... Procrustes: rmse 0.01238074  max resid 0.04968794 
    ## Run 398 stress 0.07927536 
    ## ... Procrustes: rmse 4.971286e-05  max resid 0.0001382076 
    ## ... Similar to previous best
    ## Run 399 stress 0.07927534 
    ## ... Procrustes: rmse 1.729848e-06  max resid 3.719871e-06 
    ## ... Similar to previous best
    ## Run 400 stress 0.07969258 
    ## ... Procrustes: rmse 0.0123786  max resid 0.04971973 
    ## Run 401 stress 0.08985605 
    ## Run 402 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734006  max resid 0.05633526 
    ## Run 403 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734397  max resid 0.0563619 
    ## Run 404 stress 0.09043739 
    ## Run 405 stress 0.07969265 
    ## ... Procrustes: rmse 0.01236519  max resid 0.04955505 
    ## Run 406 stress 0.07927536 
    ## ... Procrustes: rmse 3.883795e-05  max resid 0.0001046021 
    ## ... Similar to previous best
    ## Run 407 stress 0.09296218 
    ## Run 408 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238326  max resid 0.04975313 
    ## Run 409 stress 0.09296218 
    ## Run 410 stress 0.104507 
    ## Run 411 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734603  max resid 0.05637402 
    ## Run 412 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047544  max resid 0.03248182 
    ## Run 413 stress 0.07927534 
    ## ... Procrustes: rmse 2.304586e-05  max resid 6.408458e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.1089998 
    ## Run 415 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047324  max resid 0.03248227 
    ## Run 416 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735438  max resid 0.05644915 
    ## Run 417 stress 0.09022805 
    ## Run 418 stress 0.07927535 
    ## ... Procrustes: rmse 3.294455e-05  max resid 9.138596e-05 
    ## ... Similar to previous best
    ## Run 419 stress 0.0934409 
    ## Run 420 stress 0.09172795 
    ## Run 421 stress 0.07969264 
    ## ... Procrustes: rmse 0.01236811  max resid 0.04957852 
    ## Run 422 stress 0.09172792 
    ## Run 423 stress 0.07927534 
    ## ... Procrustes: rmse 5.159238e-06  max resid 1.231108e-05 
    ## ... Similar to previous best
    ## Run 424 stress 0.09106034 
    ## Run 425 stress 0.08985606 
    ## Run 426 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046299  max resid 0.03247479 
    ## Run 427 stress 0.07927534 
    ## ... Procrustes: rmse 7.348634e-06  max resid 1.781849e-05 
    ## ... Similar to previous best
    ## Run 428 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733613  max resid 0.05635503 
    ## Run 429 stress 0.07927534 
    ## ... Procrustes: rmse 2.216604e-05  max resid 5.572833e-05 
    ## ... Similar to previous best
    ## Run 430 stress 0.07927534 
    ## ... Procrustes: rmse 1.740994e-05  max resid 4.929135e-05 
    ## ... Similar to previous best
    ## Run 431 stress 0.07969261 
    ## ... Procrustes: rmse 0.0123976  max resid 0.04984066 
    ## Run 432 stress 0.0793695 
    ## ... Procrustes: rmse 0.009050055  max resid 0.03249562 
    ## Run 433 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734564  max resid 0.05637221 
    ## Run 434 stress 0.07936955 
    ## ... Procrustes: rmse 0.009025308  max resid 0.03234066 
    ## Run 435 stress 0.08985603 
    ## Run 436 stress 0.07927535 
    ## ... Procrustes: rmse 1.187636e-05  max resid 3.927522e-05 
    ## ... Similar to previous best
    ## Run 437 stress 0.09290692 
    ## Run 438 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046664  max resid 0.03247818 
    ## Run 439 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239236  max resid 0.04980183 
    ## Run 440 stress 0.09106029 
    ## Run 441 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237866  max resid 0.04975087 
    ## Run 442 stress 0.1142079 
    ## Run 443 stress 0.09106032 
    ## Run 444 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237626  max resid 0.04970667 
    ## Run 445 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237465  max resid 0.04969643 
    ## Run 446 stress 0.07927534 
    ## ... Procrustes: rmse 1.118148e-05  max resid 3.066293e-05 
    ## ... Similar to previous best
    ## Run 447 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039408  max resid 0.03242878 
    ## Run 448 stress 0.07969262 
    ## ... Procrustes: rmse 0.01238228  max resid 0.04977884 
    ## Run 449 stress 0.07927535 
    ## ... Procrustes: rmse 1.746527e-05  max resid 4.25535e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.07936951 
    ## ... Procrustes: rmse 0.009052355  max resid 0.03249911 
    ## Run 451 stress 0.08985604 
    ## Run 452 stress 0.07927536 
    ## ... Procrustes: rmse 4.977858e-05  max resid 0.0001285044 
    ## ... Similar to previous best
    ## Run 453 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043328  max resid 0.0324525 
    ## Run 454 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239236  max resid 0.04976155 
    ## Run 455 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045014  max resid 0.03246923 
    ## Run 456 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735017  max resid 0.05640187 
    ## Run 457 stress 0.09474067 
    ## Run 458 stress 0.07927534 
    ## ... Procrustes: rmse 3.24205e-06  max resid 9.29765e-06 
    ## ... Similar to previous best
    ## Run 459 stress 0.09388344 
    ## Run 460 stress 0.09022778 
    ## Run 461 stress 0.09022783 
    ## Run 462 stress 0.09247104 
    ## Run 463 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047761  max resid 0.03248596 
    ## Run 464 stress 0.09106035 
    ## Run 465 stress 0.1097212 
    ## Run 466 stress 0.07927535 
    ## ... Procrustes: rmse 3.212774e-05  max resid 8.778509e-05 
    ## ... Similar to previous best
    ## Run 467 stress 0.0796926 
    ## ... Procrustes: rmse 0.0123988  max resid 0.04980644 
    ## Run 468 stress 0.09382087 
    ## Run 469 stress 0.09296253 
    ## Run 470 stress 0.09106028 
    ## Run 471 stress 0.09290692 
    ## Run 472 stress 0.08985611 
    ## Run 473 stress 0.07936952 
    ## ... Procrustes: rmse 0.009033395  max resid 0.0323909 
    ## Run 474 stress 0.08985612 
    ## Run 475 stress 0.1089997 
    ## Run 476 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237309  max resid 0.04968274 
    ## Run 477 stress 0.08985609 
    ## Run 478 stress 0.07927538 
    ## ... Procrustes: rmse 1.98872e-05  max resid 6.128814e-05 
    ## ... Similar to previous best
    ## Run 479 stress 0.09344097 
    ## Run 480 stress 0.09565367 
    ## Run 481 stress 0.08985614 
    ## Run 482 stress 0.07969263 
    ## ... Procrustes: rmse 0.01236763  max resid 0.0495776 
    ## Run 483 stress 0.07927535 
    ## ... Procrustes: rmse 1.261413e-05  max resid 2.849784e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.09043771 
    ## Run 485 stress 0.0929622 
    ## Run 486 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237646  max resid 0.04968559 
    ## Run 487 stress 0.09072911 
    ## Run 488 stress 0.09247103 
    ## Run 489 stress 0.09106052 
    ## Run 490 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734816  max resid 0.05637792 
    ## Run 491 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049751  max resid 0.03250196 
    ## Run 492 stress 0.08985609 
    ## Run 493 stress 0.07927534 
    ## ... Procrustes: rmse 2.527134e-06  max resid 5.427614e-06 
    ## ... Similar to previous best
    ## Run 494 stress 0.0795144 
    ## ... Procrustes: rmse 0.017347  max resid 0.056377 
    ## Run 495 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044084  max resid 0.03246556 
    ## Run 496 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735023  max resid 0.05640097 
    ## Run 497 stress 0.08985609 
    ## Run 498 stress 0.08985605 
    ## Run 499 stress 0.1125337 
    ## Run 500 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734974  max resid 0.05639827 
    ## *** Best solution repeated 40 times

``` r
round(PD_beta_NMDS$stress, digits = 2)
```

    ## [1] 0.08

``` r
PD_beta_rep_NMDS <- metaMDS(PD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.3173923 
    ## Run 1 stress 0.317187 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1823507  max resid 0.4269161 
    ## Run 2 stress 0.3247522 
    ## Run 3 stress 0.3364812 
    ## Run 4 stress 0.3366198 
    ## Run 5 stress 0.3332025 
    ## Run 6 stress 0.3409185 
    ## Run 7 stress 0.3198697 
    ## Run 8 stress 0.320794 
    ## Run 9 stress 0.3199274 
    ## Run 10 stress 0.3313277 
    ## Run 11 stress 0.3232835 
    ## Run 12 stress 0.3302472 
    ## Run 13 stress 0.3223389 
    ## Run 14 stress 0.3336998 
    ## Run 15 stress 0.3389413 
    ## Run 16 stress 0.3403138 
    ## Run 17 stress 0.3174219 
    ## ... Procrustes: rmse 0.1937749  max resid 0.2958065 
    ## Run 18 stress 0.3198177 
    ## Run 19 stress 0.3237462 
    ## Run 20 stress 0.3300829 
    ## Run 21 stress 0.3410338 
    ## Run 22 stress 0.3266951 
    ## Run 23 stress 0.3336668 
    ## Run 24 stress 0.3205935 
    ## Run 25 stress 0.3310028 
    ## Run 26 stress 0.3239098 
    ## Run 27 stress 0.3412461 
    ## Run 28 stress 0.3159145 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1572167  max resid 0.3652745 
    ## Run 29 stress 0.3374269 
    ## Run 30 stress 0.3142529 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1228915  max resid 0.3979812 
    ## Run 31 stress 0.3236333 
    ## Run 32 stress 0.3314693 
    ## Run 33 stress 0.3196885 
    ## Run 34 stress 0.3170154 
    ## Run 35 stress 0.3327346 
    ## Run 36 stress 0.3281987 
    ## Run 37 stress 0.3297494 
    ## Run 38 stress 0.3370925 
    ## Run 39 stress 0.3207298 
    ## Run 40 stress 0.3198202 
    ## Run 41 stress 0.3353005 
    ## Run 42 stress 0.3185171 
    ## Run 43 stress 0.3318145 
    ## Run 44 stress 0.3265244 
    ## Run 45 stress 0.3299046 
    ## Run 46 stress 0.3164207 
    ## Run 47 stress 0.3376705 
    ## Run 48 stress 0.3304891 
    ## Run 49 stress 0.3284022 
    ## Run 50 stress 0.3193402 
    ## Run 51 stress 0.3242652 
    ## Run 52 stress 0.3254061 
    ## Run 53 stress 0.33366 
    ## Run 54 stress 0.3253655 
    ## Run 55 stress 0.3163635 
    ## Run 56 stress 0.327258 
    ## Run 57 stress 0.3215366 
    ## Run 58 stress 0.3177451 
    ## Run 59 stress 0.3220484 
    ## Run 60 stress 0.3271592 
    ## Run 61 stress 0.3414369 
    ## Run 62 stress 0.3257825 
    ## Run 63 stress 0.3294163 
    ## Run 64 stress 0.3224749 
    ## Run 65 stress 0.3272938 
    ## Run 66 stress 0.3301684 
    ## Run 67 stress 0.3218125 
    ## Run 68 stress 0.3206041 
    ## Run 69 stress 0.3203071 
    ## Run 70 stress 0.3376841 
    ## Run 71 stress 0.3338398 
    ## Run 72 stress 0.3372768 
    ## Run 73 stress 0.328279 
    ## Run 74 stress 0.3137695 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1149402  max resid 0.27336 
    ## Run 75 stress 0.3773896 
    ## Run 76 stress 0.3338181 
    ## Run 77 stress 0.3306892 
    ## Run 78 stress 0.3362366 
    ## Run 79 stress 0.3318176 
    ## Run 80 stress 0.3300291 
    ## Run 81 stress 0.3370387 
    ## Run 82 stress 0.3154256 
    ## Run 83 stress 0.3390357 
    ## Run 84 stress 0.3220651 
    ## Run 85 stress 0.3308418 
    ## Run 86 stress 0.3310107 
    ## Run 87 stress 0.3483742 
    ## Run 88 stress 0.3170689 
    ## Run 89 stress 0.3377526 
    ## Run 90 stress 0.3329868 
    ## Run 91 stress 0.3773515 
    ## Run 92 stress 0.3116016 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1810788  max resid 0.3215802 
    ## Run 93 stress 0.3289757 
    ## Run 94 stress 0.3242085 
    ## Run 95 stress 0.3303415 
    ## Run 96 stress 0.3172899 
    ## Run 97 stress 0.3305898 
    ## Run 98 stress 0.3301445 
    ## Run 99 stress 0.3316563 
    ## Run 100 stress 0.3259528 
    ## Run 101 stress 0.3321738 
    ## Run 102 stress 0.3284103 
    ## Run 103 stress 0.332509 
    ## Run 104 stress 0.3212351 
    ## Run 105 stress 0.3159459 
    ## Run 106 stress 0.3288759 
    ## Run 107 stress 0.3197076 
    ## Run 108 stress 0.3222866 
    ## Run 109 stress 0.3186841 
    ## Run 110 stress 0.3228457 
    ## Run 111 stress 0.3274408 
    ## Run 112 stress 0.3206715 
    ## Run 113 stress 0.321321 
    ## Run 114 stress 0.3266407 
    ## Run 115 stress 0.3310526 
    ## Run 116 stress 0.3263564 
    ## Run 117 stress 0.3351822 
    ## Run 118 stress 0.3341148 
    ## Run 119 stress 0.3183997 
    ## Run 120 stress 0.3297444 
    ## Run 121 stress 0.3200156 
    ## Run 122 stress 0.3184208 
    ## Run 123 stress 0.3229455 
    ## Run 124 stress 0.3288851 
    ## Run 125 stress 0.338916 
    ## Run 126 stress 0.3267079 
    ## Run 127 stress 0.3282768 
    ## Run 128 stress 0.3170971 
    ## Run 129 stress 0.3376632 
    ## Run 130 stress 0.3235783 
    ## Run 131 stress 0.3363545 
    ## Run 132 stress 0.3291057 
    ## Run 133 stress 0.3381073 
    ## Run 134 stress 0.3264197 
    ## Run 135 stress 0.3235454 
    ## Run 136 stress 0.3338358 
    ## Run 137 stress 0.3239422 
    ## Run 138 stress 0.333357 
    ## Run 139 stress 0.3290793 
    ## Run 140 stress 0.3321603 
    ## Run 141 stress 0.3318935 
    ## Run 142 stress 0.3292916 
    ## Run 143 stress 0.327982 
    ## Run 144 stress 0.3195956 
    ## Run 145 stress 0.3295349 
    ## Run 146 stress 0.3205499 
    ## Run 147 stress 0.3328325 
    ## Run 148 stress 0.3322671 
    ## Run 149 stress 0.3208405 
    ## Run 150 stress 0.3286485 
    ## Run 151 stress 0.3353119 
    ## Run 152 stress 0.3236531 
    ## Run 153 stress 0.3155503 
    ## Run 154 stress 0.3273803 
    ## Run 155 stress 0.320252 
    ## Run 156 stress 0.3316833 
    ## Run 157 stress 0.3299862 
    ## Run 158 stress 0.3239008 
    ## Run 159 stress 0.333908 
    ## Run 160 stress 0.3144 
    ## Run 161 stress 0.3282286 
    ## Run 162 stress 0.3190966 
    ## Run 163 stress 0.3333971 
    ## Run 164 stress 0.3217292 
    ## Run 165 stress 0.3215227 
    ## Run 166 stress 0.3310271 
    ## Run 167 stress 0.3242492 
    ## Run 168 stress 0.3210324 
    ## Run 169 stress 0.3349316 
    ## Run 170 stress 0.3369603 
    ## Run 171 stress 0.330683 
    ## Run 172 stress 0.3257975 
    ## Run 173 stress 0.3271683 
    ## Run 174 stress 0.3328448 
    ## Run 175 stress 0.3183245 
    ## Run 176 stress 0.3252822 
    ## Run 177 stress 0.339323 
    ## Run 178 stress 0.3157797 
    ## Run 179 stress 0.3279358 
    ## Run 180 stress 0.3179303 
    ## Run 181 stress 0.3133642 
    ## Run 182 stress 0.3274188 
    ## Run 183 stress 0.326203 
    ## Run 184 stress 0.3303948 
    ## Run 185 stress 0.3351748 
    ## Run 186 stress 0.3388019 
    ## Run 187 stress 0.341734 
    ## Run 188 stress 0.3332243 
    ## Run 189 stress 0.3236501 
    ## Run 190 stress 0.3157422 
    ## Run 191 stress 0.3188614 
    ## Run 192 stress 0.3350443 
    ## Run 193 stress 0.3121251 
    ## Run 194 stress 0.3198298 
    ## Run 195 stress 0.3367609 
    ## Run 196 stress 0.3148087 
    ## Run 197 stress 0.3309687 
    ## Run 198 stress 0.3321328 
    ## Run 199 stress 0.3345358 
    ## Run 200 stress 0.3353846 
    ## Run 201 stress 0.3387555 
    ## Run 202 stress 0.3244758 
    ## Run 203 stress 0.3162139 
    ## Run 204 stress 0.3182126 
    ## Run 205 stress 0.3223991 
    ## Run 206 stress 0.3323211 
    ## Run 207 stress 0.3332159 
    ## Run 208 stress 0.333574 
    ## Run 209 stress 0.3299975 
    ## Run 210 stress 0.3278384 
    ## Run 211 stress 0.3289356 
    ## Run 212 stress 0.313718 
    ## Run 213 stress 0.3280116 
    ## Run 214 stress 0.3363302 
    ## Run 215 stress 0.3200825 
    ## Run 216 stress 0.3208468 
    ## Run 217 stress 0.3311361 
    ## Run 218 stress 0.3168054 
    ## Run 219 stress 0.3263547 
    ## Run 220 stress 0.3253134 
    ## Run 221 stress 0.3234555 
    ## Run 222 stress 0.3227186 
    ## Run 223 stress 0.3300168 
    ## Run 224 stress 0.3283633 
    ## Run 225 stress 0.3276281 
    ## Run 226 stress 0.3362194 
    ## Run 227 stress 0.3215282 
    ## Run 228 stress 0.3191265 
    ## Run 229 stress 0.3164224 
    ## Run 230 stress 0.3159937 
    ## Run 231 stress 0.3134403 
    ## Run 232 stress 0.3205277 
    ## Run 233 stress 0.3151557 
    ## Run 234 stress 0.3279143 
    ## Run 235 stress 0.3240165 
    ## Run 236 stress 0.3239895 
    ## Run 237 stress 0.3373546 
    ## Run 238 stress 0.33375 
    ## Run 239 stress 0.3260881 
    ## Run 240 stress 0.3154305 
    ## Run 241 stress 0.3275926 
    ## Run 242 stress 0.317904 
    ## Run 243 stress 0.3337546 
    ## Run 244 stress 0.3301427 
    ## Run 245 stress 0.3365824 
    ## Run 246 stress 0.3384613 
    ## Run 247 stress 0.3334389 
    ## Run 248 stress 0.326968 
    ## Run 249 stress 0.3299367 
    ## Run 250 stress 0.3234106 
    ## Run 251 stress 0.3141809 
    ## Run 252 stress 0.3211966 
    ## Run 253 stress 0.330865 
    ## Run 254 stress 0.3265077 
    ## Run 255 stress 0.3192049 
    ## Run 256 stress 0.3340613 
    ## Run 257 stress 0.3263121 
    ## Run 258 stress 0.3198728 
    ## Run 259 stress 0.3250737 
    ## Run 260 stress 0.3373659 
    ## Run 261 stress 0.3342677 
    ## Run 262 stress 0.3216801 
    ## Run 263 stress 0.3146589 
    ## Run 264 stress 0.3357334 
    ## Run 265 stress 0.3155085 
    ## Run 266 stress 0.3279867 
    ## Run 267 stress 0.3144332 
    ## Run 268 stress 0.3209029 
    ## Run 269 stress 0.3221943 
    ## Run 270 stress 0.3306958 
    ## Run 271 stress 0.3316676 
    ## Run 272 stress 0.3159941 
    ## Run 273 stress 0.316438 
    ## Run 274 stress 0.3316318 
    ## Run 275 stress 0.3271567 
    ## Run 276 stress 0.3355963 
    ## Run 277 stress 0.3396976 
    ## Run 278 stress 0.3311063 
    ## Run 279 stress 0.3258284 
    ## Run 280 stress 0.3135339 
    ## Run 281 stress 0.3240667 
    ## Run 282 stress 0.3215131 
    ## Run 283 stress 0.3150447 
    ## Run 284 stress 0.3362261 
    ## Run 285 stress 0.3240946 
    ## Run 286 stress 0.3310611 
    ## Run 287 stress 0.3182775 
    ## Run 288 stress 0.3313377 
    ## Run 289 stress 0.3370848 
    ## Run 290 stress 0.3191024 
    ## Run 291 stress 0.3240533 
    ## Run 292 stress 0.3381603 
    ## Run 293 stress 0.3397747 
    ## Run 294 stress 0.3249535 
    ## Run 295 stress 0.3358005 
    ## Run 296 stress 0.334093 
    ## Run 297 stress 0.3260729 
    ## Run 298 stress 0.3202738 
    ## Run 299 stress 0.3216525 
    ## Run 300 stress 0.3216415 
    ## Run 301 stress 0.3195049 
    ## Run 302 stress 0.3294689 
    ## Run 303 stress 0.3164968 
    ## Run 304 stress 0.328135 
    ## Run 305 stress 0.3411275 
    ## Run 306 stress 0.3271465 
    ## Run 307 stress 0.315874 
    ## Run 308 stress 0.3187453 
    ## Run 309 stress 0.3212949 
    ## Run 310 stress 0.329125 
    ## Run 311 stress 0.3237038 
    ## Run 312 stress 0.3166962 
    ## Run 313 stress 0.3255723 
    ## Run 314 stress 0.3337464 
    ## Run 315 stress 0.3226692 
    ## Run 316 stress 0.3237349 
    ## Run 317 stress 0.3295486 
    ## Run 318 stress 0.329226 
    ## Run 319 stress 0.3283694 
    ## Run 320 stress 0.324744 
    ## Run 321 stress 0.3157421 
    ## Run 322 stress 0.3224816 
    ## Run 323 stress 0.322053 
    ## Run 324 stress 0.3339416 
    ## Run 325 stress 0.3256991 
    ## Run 326 stress 0.3205676 
    ## Run 327 stress 0.3225403 
    ## Run 328 stress 0.3185361 
    ## Run 329 stress 0.3357257 
    ## Run 330 stress 0.3206023 
    ## Run 331 stress 0.3310359 
    ## Run 332 stress 0.3375943 
    ## Run 333 stress 0.3317682 
    ## Run 334 stress 0.3332276 
    ## Run 335 stress 0.3180939 
    ## Run 336 stress 0.3352403 
    ## Run 337 stress 0.3205059 
    ## Run 338 stress 0.3282836 
    ## Run 339 stress 0.329139 
    ## Run 340 stress 0.3424487 
    ## Run 341 stress 0.3162962 
    ## Run 342 stress 0.3206543 
    ## Run 343 stress 0.3379235 
    ## Run 344 stress 0.3301145 
    ## Run 345 stress 0.3192395 
    ## Run 346 stress 0.3385774 
    ## Run 347 stress 0.3264283 
    ## Run 348 stress 0.3369078 
    ## Run 349 stress 0.3303038 
    ## Run 350 stress 0.3206224 
    ## Run 351 stress 0.3151443 
    ## Run 352 stress 0.3225397 
    ## Run 353 stress 0.3210237 
    ## Run 354 stress 0.3419319 
    ## Run 355 stress 0.3454631 
    ## Run 356 stress 0.3203194 
    ## Run 357 stress 0.3268482 
    ## Run 358 stress 0.3354762 
    ## Run 359 stress 0.3359545 
    ## Run 360 stress 0.3372843 
    ## Run 361 stress 0.333812 
    ## Run 362 stress 0.335371 
    ## Run 363 stress 0.3174221 
    ## Run 364 stress 0.3321186 
    ## Run 365 stress 0.3243905 
    ## Run 366 stress 0.3244764 
    ## Run 367 stress 0.3366272 
    ## Run 368 stress 0.3209841 
    ## Run 369 stress 0.337669 
    ## Run 370 stress 0.3288402 
    ## Run 371 stress 0.3251209 
    ## Run 372 stress 0.3393954 
    ## Run 373 stress 0.3255395 
    ## Run 374 stress 0.320657 
    ## Run 375 stress 0.3294823 
    ## Run 376 stress 0.3310366 
    ## Run 377 stress 0.3270823 
    ## Run 378 stress 0.331413 
    ## Run 379 stress 0.3193724 
    ## Run 380 stress 0.3162443 
    ## Run 381 stress 0.3152909 
    ## Run 382 stress 0.3185155 
    ## Run 383 stress 0.3158907 
    ## Run 384 stress 0.3353318 
    ## Run 385 stress 0.3182773 
    ## Run 386 stress 0.3361393 
    ## Run 387 stress 0.3295081 
    ## Run 388 stress 0.3371007 
    ## Run 389 stress 0.3255891 
    ## Run 390 stress 0.3382175 
    ## Run 391 stress 0.3333989 
    ## Run 392 stress 0.3279199 
    ## Run 393 stress 0.3280922 
    ## Run 394 stress 0.3239941 
    ## Run 395 stress 0.3307714 
    ## Run 396 stress 0.3422017 
    ## Run 397 stress 0.3289868 
    ## Run 398 stress 0.3192354 
    ## Run 399 stress 0.3420399 
    ## Run 400 stress 0.3283356 
    ## Run 401 stress 0.3289906 
    ## Run 402 stress 0.3167883 
    ## Run 403 stress 0.3207606 
    ## Run 404 stress 0.3180175 
    ## Run 405 stress 0.3388847 
    ## Run 406 stress 0.3284138 
    ## Run 407 stress 0.3303182 
    ## Run 408 stress 0.3380142 
    ## Run 409 stress 0.3329487 
    ## Run 410 stress 0.3194513 
    ## Run 411 stress 0.3201348 
    ## Run 412 stress 0.3266024 
    ## Run 413 stress 0.3384617 
    ## Run 414 stress 0.3308703 
    ## Run 415 stress 0.3425648 
    ## Run 416 stress 0.3330154 
    ## Run 417 stress 0.3272838 
    ## Run 418 stress 0.3229949 
    ## Run 419 stress 0.3325584 
    ## Run 420 stress 0.3299934 
    ## Run 421 stress 0.3213845 
    ## Run 422 stress 0.3298176 
    ## Run 423 stress 0.3330309 
    ## Run 424 stress 0.3344867 
    ## Run 425 stress 0.327845 
    ## Run 426 stress 0.3234448 
    ## Run 427 stress 0.3341723 
    ## Run 428 stress 0.3298378 
    ## Run 429 stress 0.3191714 
    ## Run 430 stress 0.3199913 
    ## Run 431 stress 0.3381511 
    ## Run 432 stress 0.3222459 
    ## Run 433 stress 0.3229461 
    ## Run 434 stress 0.3216257 
    ## Run 435 stress 0.3229924 
    ## Run 436 stress 0.3302756 
    ## Run 437 stress 0.337193 
    ## Run 438 stress 0.3287392 
    ## Run 439 stress 0.3329966 
    ## Run 440 stress 0.3189527 
    ## Run 441 stress 0.3294051 
    ## Run 442 stress 0.334908 
    ## Run 443 stress 0.3241666 
    ## Run 444 stress 0.3279217 
    ## Run 445 stress 0.3249616 
    ## Run 446 stress 0.3326274 
    ## Run 447 stress 0.3329308 
    ## Run 448 stress 0.3315654 
    ## Run 449 stress 0.3261203 
    ## Run 450 stress 0.3334908 
    ## Run 451 stress 0.3222604 
    ## Run 452 stress 0.3391906 
    ## Run 453 stress 0.338644 
    ## Run 454 stress 0.3267392 
    ## Run 455 stress 0.3196053 
    ## Run 456 stress 0.3324216 
    ## Run 457 stress 0.316899 
    ## Run 458 stress 0.3199159 
    ## Run 459 stress 0.3184071 
    ## Run 460 stress 0.326899 
    ## Run 461 stress 0.3371942 
    ## Run 462 stress 0.3216368 
    ## Run 463 stress 0.3371303 
    ## Run 464 stress 0.3349934 
    ## Run 465 stress 0.3313947 
    ## Run 466 stress 0.336037 
    ## Run 467 stress 0.312013 
    ## ... Procrustes: rmse 0.1062699  max resid 0.2820896 
    ## Run 468 stress 0.322492 
    ## Run 469 stress 0.3207067 
    ## Run 470 stress 0.3236695 
    ## Run 471 stress 0.3203623 
    ## Run 472 stress 0.3177049 
    ## Run 473 stress 0.3264661 
    ## Run 474 stress 0.3284562 
    ## Run 475 stress 0.3260948 
    ## Run 476 stress 0.3298568 
    ## Run 477 stress 0.3170301 
    ## Run 478 stress 0.3195224 
    ## Run 479 stress 0.3231544 
    ## Run 480 stress 0.3368399 
    ## Run 481 stress 0.319266 
    ## Run 482 stress 0.3165763 
    ## Run 483 stress 0.3136102 
    ## Run 484 stress 0.3280835 
    ## Run 485 stress 0.3216294 
    ## Run 486 stress 0.3194799 
    ## Run 487 stress 0.3302471 
    ## Run 488 stress 0.3351692 
    ## Run 489 stress 0.3195283 
    ## Run 490 stress 0.3191434 
    ## Run 491 stress 0.3350188 
    ## Run 492 stress 0.3201873 
    ## Run 493 stress 0.3259292 
    ## Run 494 stress 0.3372827 
    ## Run 495 stress 0.3351318 
    ## Run 496 stress 0.3227499 
    ## Run 497 stress 0.3325364 
    ## Run 498 stress 0.339675 
    ## Run 499 stress 0.324018 
    ## Run 500 stress 0.31579 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##      1: no. of iterations >= maxit
    ##    499: stress ratio > sratmax

``` r
PD_beta_ric_NMDS <- metaMDS(PD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.01115177 
    ## Run 1 stress 0.01117184 
    ## ... Procrustes: rmse 0.004210336  max resid 0.00948797 
    ## ... Similar to previous best
    ## Run 2 stress 0.01545425 
    ## Run 3 stress 0.01115452 
    ## ... Procrustes: rmse 0.002143607  max resid 0.004832637 
    ## ... Similar to previous best
    ## Run 4 stress 0.01528187 
    ## Run 5 stress 0.01118475 
    ## ... Procrustes: rmse 0.005163155  max resid 0.01163142 
    ## Run 6 stress 0.01547442 
    ## Run 7 stress 0.01547548 
    ## Run 8 stress 0.01117701 
    ## ... Procrustes: rmse 0.004621212  max resid 0.01041261 
    ## Run 9 stress 0.01561833 
    ## Run 10 stress 0.01529722 
    ## Run 11 stress 0.01546104 
    ## Run 12 stress 0.01527908 
    ## Run 13 stress 0.01528208 
    ## Run 14 stress 0.01718687 
    ## Run 15 stress 0.01553814 
    ## Run 16 stress 0.01686832 
    ## Run 17 stress 0.01557914 
    ## Run 18 stress 0.01498576 
    ## Run 19 stress 0.0111519 
    ## ... Procrustes: rmse 5.078247e-05  max resid 0.0001148014 
    ## ... Similar to previous best
    ## Run 20 stress 0.01498578 
    ## Run 21 stress 0.01542829 
    ## Run 22 stress 0.0152719 
    ## Run 23 stress 0.01115144 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001220334  max resid 0.00275158 
    ## ... Similar to previous best
    ## Run 24 stress 0.01567585 
    ## Run 25 stress 0.01515719 
    ## Run 26 stress 0.01147228 
    ## ... Procrustes: rmse 0.0134113  max resid 0.02973911 
    ## Run 27 stress 0.01528218 
    ## Run 28 stress 0.01542788 
    ## Run 29 stress 0.01549913 
    ## Run 30 stress 0.01557887 
    ## Run 31 stress 0.01548291 
    ## Run 32 stress 0.01558336 
    ## Run 33 stress 0.01571533 
    ## Run 34 stress 0.01498584 
    ## Run 35 stress 0.01545378 
    ## Run 36 stress 0.01499504 
    ## Run 37 stress 0.01499237 
    ## Run 38 stress 0.01556368 
    ## Run 39 stress 0.01577185 
    ## Run 40 stress 0.01159861 
    ## ... Procrustes: rmse 0.01541845  max resid 0.03319352 
    ## Run 41 stress 0.01115194 
    ## ... Procrustes: rmse 0.00128476  max resid 0.002900922 
    ## ... Similar to previous best
    ## Run 42 stress 0.01542807 
    ## Run 43 stress 0.0149922 
    ## Run 44 stress 0.01556106 
    ## Run 45 stress 0.01528201 
    ## Run 46 stress 0.01670173 
    ## Run 47 stress 0.01515699 
    ## Run 48 stress 0.01115175 
    ## ... Procrustes: rmse 0.001210063  max resid 0.002731124 
    ## ... Similar to previous best
    ## Run 49 stress 0.0158938 
    ## Run 50 stress 0.01564566 
    ## Run 51 stress 0.01115189 
    ## ... Procrustes: rmse 0.001266889  max resid 0.002860305 
    ## ... Similar to previous best
    ## Run 52 stress 0.0167439 
    ## Run 53 stress 0.0153616 
    ## Run 54 stress 0.01550467 
    ## Run 55 stress 0.3807941 
    ## Run 56 stress 0.01525794 
    ## Run 57 stress 0.01641888 
    ## Run 58 stress 0.01665895 
    ## Run 59 stress 0.01252527 
    ## Run 60 stress 0.01115199 
    ## ... Procrustes: rmse 0.001308859  max resid 0.00295511 
    ## ... Similar to previous best
    ## Run 61 stress 0.01115619 
    ## ... Procrustes: rmse 0.001230846  max resid 0.002777257 
    ## ... Similar to previous best
    ## Run 62 stress 0.01115324 
    ## ... Procrustes: rmse 0.0006331283  max resid 0.001428861 
    ## ... Similar to previous best
    ## Run 63 stress 0.0156153 
    ## Run 64 stress 0.01529492 
    ## Run 65 stress 0.01115193 
    ## ... Procrustes: rmse 0.001287526  max resid 0.002906939 
    ## ... Similar to previous best
    ## Run 66 stress 0.0154793 
    ## Run 67 stress 0.01115626 
    ## ... Procrustes: rmse 0.001244009  max resid 0.002807241 
    ## ... Similar to previous best
    ## Run 68 stress 0.01132092 
    ## ... Procrustes: rmse 0.009168337  max resid 0.02062988 
    ## Run 69 stress 0.01556947 
    ## Run 70 stress 0.01124814 
    ## ... Procrustes: rmse 0.007043197  max resid 0.01586224 
    ## Run 71 stress 0.01673564 
    ## Run 72 stress 0.01493542 
    ## Run 73 stress 0.0152906 
    ## Run 74 stress 0.01515715 
    ## Run 75 stress 0.01546852 
    ## Run 76 stress 0.01147662 
    ## ... Procrustes: rmse 0.01183465  max resid 0.02575568 
    ## Run 77 stress 0.01561938 
    ## Run 78 stress 0.01558344 
    ## Run 79 stress 0.01552957 
    ## Run 80 stress 0.01515707 
    ## Run 81 stress 0.0112203 
    ## ... Procrustes: rmse 0.005883306  max resid 0.0132493 
    ## Run 82 stress 0.01554175 
    ## Run 83 stress 0.01663507 
    ## Run 84 stress 0.01544843 
    ## Run 85 stress 0.01115176 
    ## ... Procrustes: rmse 0.001214865  max resid 0.002742342 
    ## ... Similar to previous best
    ## Run 86 stress 0.01515732 
    ## Run 87 stress 0.01493163 
    ## Run 88 stress 0.01122796 
    ## ... Procrustes: rmse 0.006116515  max resid 0.01376639 
    ## Run 89 stress 0.01576947 
    ## Run 90 stress 0.01569323 
    ## Run 91 stress 0.01556656 
    ## Run 92 stress 0.01536254 
    ## Run 93 stress 0.01549966 
    ## Run 94 stress 0.01515745 
    ## Run 95 stress 0.0151574 
    ## Run 96 stress 0.2867431 
    ## Run 97 stress 0.01582 
    ## Run 98 stress 0.01529522 
    ## Run 99 stress 0.01547439 
    ## Run 100 stress 0.01548355 
    ## Run 101 stress 0.015338 
    ## Run 102 stress 0.01695065 
    ## Run 103 stress 0.01166451 
    ## Run 104 stress 0.01536151 
    ## Run 105 stress 0.01543096 
    ## Run 106 stress 0.01557899 
    ## Run 107 stress 0.01221332 
    ## Run 108 stress 0.01515142 
    ## Run 109 stress 0.01557919 
    ## Run 110 stress 0.0156722 
    ## Run 111 stress 0.01549944 
    ## Run 112 stress 0.01658387 
    ## Run 113 stress 0.01564592 
    ## Run 114 stress 0.0155696 
    ## Run 115 stress 0.01564567 
    ## Run 116 stress 0.01685172 
    ## Run 117 stress 0.01533781 
    ## Run 118 stress 0.01498579 
    ## Run 119 stress 0.01549962 
    ## Run 120 stress 0.01552203 
    ## Run 121 stress 0.01134879 
    ## ... Procrustes: rmse 0.009718625  max resid 0.0217394 
    ## Run 122 stress 0.0111703 
    ## ... Procrustes: rmse 0.002839521  max resid 0.006406154 
    ## ... Similar to previous best
    ## Run 123 stress 0.01577187 
    ## Run 124 stress 0.2816576 
    ## Run 125 stress 0.01182438 
    ## Run 126 stress 0.01641888 
    ## Run 127 stress 0.01663494 
    ## Run 128 stress 0.01549926 
    ## Run 129 stress 0.01115215 
    ## ... Procrustes: rmse 0.0003083922  max resid 0.0006961793 
    ## ... Similar to previous best
    ## Run 130 stress 0.01548168 
    ## Run 131 stress 0.0113098 
    ## ... Procrustes: rmse 0.009131141  max resid 0.02054191 
    ## Run 132 stress 0.380794 
    ## Run 133 stress 0.01493188 
    ## Run 134 stress 0.01116049 
    ## ... Procrustes: rmse 0.001848331  max resid 0.004170352 
    ## ... Similar to previous best
    ## Run 135 stress 0.01528186 
    ## Run 136 stress 0.01663503 
    ## Run 137 stress 0.01568979 
    ## Run 138 stress 0.01558482 
    ## Run 139 stress 0.01527551 
    ## Run 140 stress 0.01641897 
    ## Run 141 stress 0.01555453 
    ## Run 142 stress 0.0111517 
    ## ... Procrustes: rmse 0.0001314958  max resid 0.0002970515 
    ## ... Similar to previous best
    ## Run 143 stress 0.01498583 
    ## Run 144 stress 0.01525763 
    ## Run 145 stress 0.01574939 
    ## Run 146 stress 0.01515736 
    ## Run 147 stress 0.01549934 
    ## Run 148 stress 0.01574952 
    ## Run 149 stress 0.01641869 
    ## Run 150 stress 0.01562163 
    ## Run 151 stress 0.01559147 
    ## Run 152 stress 0.01549938 
    ## Run 153 stress 0.01528214 
    ## Run 154 stress 0.0111518 
    ## ... Procrustes: rmse 0.00123332  max resid 0.002784027 
    ## ... Similar to previous best
    ## Run 155 stress 0.01115192 
    ## ... Procrustes: rmse 0.001276568  max resid 0.002881666 
    ## ... Similar to previous best
    ## Run 156 stress 0.0164984 
    ## Run 157 stress 0.01536127 
    ## Run 158 stress 0.01648743 
    ## Run 159 stress 0.01556977 
    ## Run 160 stress 0.01131818 
    ## ... Procrustes: rmse 0.009449461  max resid 0.02114819 
    ## Run 161 stress 0.2816595 
    ## Run 162 stress 0.2819805 
    ## Run 163 stress 0.01115179 
    ## ... Procrustes: rmse 0.001229997  max resid 0.002776429 
    ## ... Similar to previous best
    ## Run 164 stress 0.01545379 
    ## Run 165 stress 0.01556119 
    ## Run 166 stress 0.01554501 
    ## Run 167 stress 0.01656446 
    ## Run 168 stress 0.01115569 
    ## ... Procrustes: rmse 0.001145662  max resid 0.002585471 
    ## ... Similar to previous best
    ## Run 169 stress 0.01115199 
    ## ... Procrustes: rmse 0.0002437669  max resid 0.0005499203 
    ## ... Similar to previous best
    ## Run 170 stress 0.01116863 
    ## ... Procrustes: rmse 0.002661036  max resid 0.006005342 
    ## ... Similar to previous best
    ## Run 171 stress 0.01115193 
    ## ... Procrustes: rmse 0.001276261  max resid 0.002880911 
    ## ... Similar to previous best
    ## Run 172 stress 0.01542942 
    ## Run 173 stress 0.0155297 
    ## Run 174 stress 0.01116259 
    ## ... Procrustes: rmse 0.00209725  max resid 0.004731628 
    ## ... Similar to previous best
    ## Run 175 stress 0.01550068 
    ## Run 176 stress 0.01515358 
    ## Run 177 stress 0.0111519 
    ## ... Procrustes: rmse 0.001275318  max resid 0.002879217 
    ## ... Similar to previous best
    ## Run 178 stress 0.01582946 
    ## Run 179 stress 0.01574952 
    ## Run 180 stress 0.01556961 
    ## Run 181 stress 0.01546466 
    ## Run 182 stress 0.01493169 
    ## Run 183 stress 0.01528165 
    ## Run 184 stress 0.01499402 
    ## Run 185 stress 0.01556099 
    ## Run 186 stress 0.01641866 
    ## Run 187 stress 0.01685501 
    ## Run 188 stress 0.0155463 
    ## Run 189 stress 0.01515675 
    ## Run 190 stress 0.01518944 
    ## Run 191 stress 0.0153615 
    ## Run 192 stress 0.01115165 
    ## ... Procrustes: rmse 0.001167615  max resid 0.00263523 
    ## ... Similar to previous best
    ## Run 193 stress 0.01522993 
    ## Run 194 stress 0.0111516 
    ## ... Procrustes: rmse 0.001145962  max resid 0.002586252 
    ## ... Similar to previous best
    ## Run 195 stress 0.01566026 
    ## Run 196 stress 0.01121256 
    ## ... Procrustes: rmse 0.005513838  max resid 0.01241841 
    ## Run 197 stress 0.01135495 
    ## ... Procrustes: rmse 0.009822508  max resid 0.02196291 
    ## Run 198 stress 0.01115493 
    ## ... Procrustes: rmse 0.001002154  max resid 0.002261888 
    ## ... Similar to previous best
    ## Run 199 stress 0.0156835 
    ## Run 200 stress 0.01141805 
    ## ... Procrustes: rmse 0.01214332  max resid 0.02699744 
    ## Run 201 stress 0.01577018 
    ## Run 202 stress 0.01536178 
    ## Run 203 stress 0.01514952 
    ## Run 204 stress 0.01543899 
    ## Run 205 stress 0.01115179 
    ## ... Procrustes: rmse 0.001227288  max resid 0.002770511 
    ## ... Similar to previous best
    ## Run 206 stress 0.01546884 
    ## Run 207 stress 0.01533795 
    ## Run 208 stress 0.01545365 
    ## Run 209 stress 0.01525789 
    ## Run 210 stress 0.01574953 
    ## Run 211 stress 0.01552916 
    ## Run 212 stress 0.01549097 
    ## Run 213 stress 0.01493154 
    ## Run 214 stress 0.014992 
    ## Run 215 stress 0.01527577 
    ## Run 216 stress 0.01515713 
    ## Run 217 stress 0.01115395 
    ## ... Procrustes: rmse 0.0007999795  max resid 0.001805272 
    ## ... Similar to previous best
    ## Run 218 stress 0.01529099 
    ## Run 219 stress 0.01555535 
    ## Run 220 stress 0.01567382 
    ## Run 221 stress 0.01515705 
    ## Run 222 stress 0.01498606 
    ## Run 223 stress 0.01546048 
    ## Run 224 stress 0.01557901 
    ## Run 225 stress 0.01555419 
    ## Run 226 stress 0.01556921 
    ## Run 227 stress 0.01116125 
    ## ... Procrustes: rmse 0.001868201  max resid 0.004217089 
    ## ... Similar to previous best
    ## Run 228 stress 0.0149314 
    ## Run 229 stress 0.01545275 
    ## Run 230 stress 0.01574922 
    ## Run 231 stress 0.01594408 
    ## Run 232 stress 0.0111529 
    ## ... Procrustes: rmse 0.0005425054  max resid 0.001224435 
    ## ... Similar to previous best
    ## Run 233 stress 0.01115184 
    ## ... Procrustes: rmse 0.0001820989  max resid 0.0004106356 
    ## ... Similar to previous best
    ## Run 234 stress 0.01518657 
    ## Run 235 stress 0.01556188 
    ## Run 236 stress 0.01498606 
    ## Run 237 stress 0.01536174 
    ## Run 238 stress 0.01525751 
    ## Run 239 stress 0.01562246 
    ## Run 240 stress 0.01564577 
    ## Run 241 stress 0.01329501 
    ## Run 242 stress 0.01545296 
    ## Run 243 stress 0.01658205 
    ## Run 244 stress 0.01161268 
    ## ... Procrustes: rmse 0.01576842  max resid 0.03394904 
    ## Run 245 stress 0.01150291 
    ## ... Procrustes: rmse 0.01247122  max resid 0.02700284 
    ## Run 246 stress 0.01670621 
    ## Run 247 stress 0.0158199 
    ## Run 248 stress 0.01118737 
    ## ... Procrustes: rmse 0.004108136  max resid 0.009263457 
    ## ... Similar to previous best
    ## Run 249 stress 0.01700456 
    ## Run 250 stress 0.0155793 
    ## Run 251 stress 0.3742621 
    ## Run 252 stress 0.01528177 
    ## Run 253 stress 0.01542821 
    ## Run 254 stress 0.01527645 
    ## Run 255 stress 0.016582 
    ## Run 256 stress 0.01116259 
    ## ... Procrustes: rmse 0.002091961  max resid 0.004720205 
    ## ... Similar to previous best
    ## Run 257 stress 0.01570974 
    ## Run 258 stress 0.01515705 
    ## Run 259 stress 0.01581945 
    ## Run 260 stress 0.01545314 
    ## Run 261 stress 0.01649858 
    ## Run 262 stress 0.01499238 
    ## Run 263 stress 0.01703363 
    ## Run 264 stress 0.01115197 
    ## ... Procrustes: rmse 0.001303926  max resid 0.002944151 
    ## ... Similar to previous best
    ## Run 265 stress 0.01542778 
    ## Run 266 stress 0.01549938 
    ## Run 267 stress 0.01115303 
    ## ... Procrustes: rmse 0.0005738432  max resid 0.001295541 
    ## ... Similar to previous best
    ## Run 268 stress 0.01554267 
    ## Run 269 stress 0.01118181 
    ## ... Procrustes: rmse 0.003744898  max resid 0.008444261 
    ## ... Similar to previous best
    ## Run 270 stress 0.01147522 
    ## ... Procrustes: rmse 0.01347924  max resid 0.02988666 
    ## Run 271 stress 0.01498594 
    ## Run 272 stress 0.01641883 
    ## Run 273 stress 0.01515714 
    ## Run 274 stress 0.01583022 
    ## Run 275 stress 0.01119083 
    ## ... Procrustes: rmse 0.004330348  max resid 0.009763785 
    ## ... Similar to previous best
    ## Run 276 stress 0.01549935 
    ## Run 277 stress 0.01115204 
    ## ... Procrustes: rmse 0.0002696073  max resid 0.000608651 
    ## ... Similar to previous best
    ## Run 278 stress 0.01571572 
    ## Run 279 stress 0.01575658 
    ## Run 280 stress 0.01515705 
    ## Run 281 stress 0.01115175 
    ## ... Procrustes: rmse 0.0001542726  max resid 0.0003484419 
    ## ... Similar to previous best
    ## Run 282 stress 0.0111516 
    ## ... Procrustes: rmse 8.251388e-05  max resid 0.0001861937 
    ## ... Similar to previous best
    ## Run 283 stress 0.01551044 
    ## Run 284 stress 0.01115161 
    ## ... Procrustes: rmse 9.135223e-05  max resid 0.0002063682 
    ## ... Similar to previous best
    ## Run 285 stress 0.01520049 
    ## Run 286 stress 0.01493163 
    ## Run 287 stress 0.01533805 
    ## Run 288 stress 0.01523604 
    ## Run 289 stress 0.0149921 
    ## Run 290 stress 0.01564582 
    ## Run 291 stress 0.01161704 
    ## ... Procrustes: rmse 0.01579835  max resid 0.034012 
    ## Run 292 stress 0.01658223 
    ## Run 293 stress 0.01556934 
    ## Run 294 stress 0.01116788 
    ## ... Procrustes: rmse 0.002572977  max resid 0.005801176 
    ## ... Similar to previous best
    ## Run 295 stress 0.01498599 
    ## Run 296 stress 0.01694506 
    ## Run 297 stress 0.01552589 
    ## Run 298 stress 0.01533776 
    ## Run 299 stress 0.01543898 
    ## Run 300 stress 0.01548182 
    ## Run 301 stress 0.01533794 
    ## Run 302 stress 0.01648595 
    ## Run 303 stress 0.01556958 
    ## Run 304 stress 0.01515676 
    ## Run 305 stress 0.01654969 
    ## Run 306 stress 0.01680936 
    ## Run 307 stress 0.01641866 
    ## Run 308 stress 0.016419 
    ## Run 309 stress 0.01544859 
    ## Run 310 stress 0.01115579 
    ## ... Procrustes: rmse 0.0003489156  max resid 0.0007766696 
    ## ... Similar to previous best
    ## Run 311 stress 0.01117516 
    ## ... Procrustes: rmse 0.003258142  max resid 0.007347874 
    ## ... Similar to previous best
    ## Run 312 stress 0.0154678 
    ## Run 313 stress 0.0151496 
    ## Run 314 stress 0.01498579 
    ## Run 315 stress 0.01138507 
    ## ... Procrustes: rmse 0.01019744  max resid 0.02269126 
    ## Run 316 stress 0.01552985 
    ## Run 317 stress 0.0155296 
    ## Run 318 stress 0.01556095 
    ## Run 319 stress 0.01700774 
    ## Run 320 stress 0.01154024 
    ## ... Procrustes: rmse 0.01350019  max resid 0.02992498 
    ## Run 321 stress 0.01115165 
    ## ... Procrustes: rmse 0.0001069937  max resid 0.0002416568 
    ## ... Similar to previous best
    ## Run 322 stress 0.01641883 
    ## Run 323 stress 0.01515659 
    ## Run 324 stress 0.01115835 
    ## ... Procrustes: rmse 0.001558833  max resid 0.003517839 
    ## ... Similar to previous best
    ## Run 325 stress 0.01115154 
    ## ... Procrustes: rmse 0.001116233  max resid 0.002519134 
    ## ... Similar to previous best
    ## Run 326 stress 0.015766 
    ## Run 327 stress 0.01493175 
    ## Run 328 stress 0.01515738 
    ## Run 329 stress 0.01547103 
    ## Run 330 stress 0.01528217 
    ## Run 331 stress 0.01116969 
    ## ... Procrustes: rmse 0.002799497  max resid 0.006315239 
    ## ... Similar to previous best
    ## Run 332 stress 0.01499228 
    ## Run 333 stress 0.01542731 
    ## Run 334 stress 0.01552923 
    ## Run 335 stress 0.0113349 
    ## ... Procrustes: rmse 0.009948535  max resid 0.0222355 
    ## Run 336 stress 0.01115221 
    ## ... Procrustes: rmse 0.001390796  max resid 0.003140918 
    ## ... Similar to previous best
    ## Run 337 stress 0.01549933 
    ## Run 338 stress 0.01558321 
    ## Run 339 stress 0.01514968 
    ## Run 340 stress 0.3380803 
    ## Run 341 stress 0.01648823 
    ## Run 342 stress 0.01498604 
    ## Run 343 stress 0.01115915 
    ## ... Procrustes: rmse 0.0016736  max resid 0.003776313 
    ## ... Similar to previous best
    ## Run 344 stress 0.01117368 
    ## ... Procrustes: rmse 0.002835242  max resid 0.006400294 
    ## ... Similar to previous best
    ## Run 345 stress 0.01536244 
    ## Run 346 stress 0.01552963 
    ## Run 347 stress 0.01536057 
    ## Run 348 stress 0.0115196 
    ## ... Procrustes: rmse 0.01308045  max resid 0.02823833 
    ## Run 349 stress 0.015561 
    ## Run 350 stress 0.01515728 
    ## Run 351 stress 0.01528173 
    ## Run 352 stress 0.01121141 
    ## ... Procrustes: rmse 0.005299722  max resid 0.01194375 
    ## Run 353 stress 0.01583074 
    ## Run 354 stress 0.01548174 
    ## Run 355 stress 0.01546982 
    ## Run 356 stress 0.01115197 
    ## ... Procrustes: rmse 0.001295884  max resid 0.002926134 
    ## ... Similar to previous best
    ## Run 357 stress 0.01562094 
    ## Run 358 stress 0.01149206 
    ## ... Procrustes: rmse 0.01227052  max resid 0.02660498 
    ## Run 359 stress 0.01685492 
    ## Run 360 stress 0.01542809 
    ## Run 361 stress 0.01115187 
    ## ... Procrustes: rmse 0.0002018005  max resid 0.0004555502 
    ## ... Similar to previous best
    ## Run 362 stress 0.01564579 
    ## Run 363 stress 0.01115193 
    ## ... Procrustes: rmse 0.0002229548  max resid 0.0005036416 
    ## ... Similar to previous best
    ## Run 364 stress 0.01549935 
    ## Run 365 stress 0.01688325 
    ## Run 366 stress 0.01515722 
    ## Run 367 stress 0.01536147 
    ## Run 368 stress 0.01546138 
    ## Run 369 stress 0.01546447 
    ## Run 370 stress 0.01498603 
    ## Run 371 stress 0.01115222 
    ## ... Procrustes: rmse 0.0003227852  max resid 0.0007281229 
    ## ... Similar to previous best
    ## Run 372 stress 0.01525805 
    ## Run 373 stress 0.01543127 
    ## Run 374 stress 0.01370521 
    ## Run 375 stress 0.01633912 
    ## Run 376 stress 0.0111568 
    ## ... Procrustes: rmse 0.001322837  max resid 0.002984167 
    ## ... Similar to previous best
    ## Run 377 stress 0.01150317 
    ## ... Procrustes: rmse 0.01260022  max resid 0.02726298 
    ## Run 378 stress 0.01122227 
    ## ... Procrustes: rmse 0.004969907  max resid 0.01117751 
    ## Run 379 stress 0.01115175 
    ## ... Procrustes: rmse 0.0001448543  max resid 0.0003266197 
    ## ... Similar to previous best
    ## Run 380 stress 0.0149924 
    ## Run 381 stress 0.0155553 
    ## Run 382 stress 0.01498592 
    ## Run 383 stress 0.01574946 
    ## Run 384 stress 0.0152946 
    ## Run 385 stress 0.01719163 
    ## Run 386 stress 0.01116874 
    ## ... Procrustes: rmse 0.002676198  max resid 0.006038677 
    ## ... Similar to previous best
    ## Run 387 stress 0.01115924 
    ## ... Procrustes: rmse 0.001686983  max resid 0.00380612 
    ## ... Similar to previous best
    ## Run 388 stress 0.01572961 
    ## Run 389 stress 0.01557901 
    ## Run 390 stress 0.01576299 
    ## Run 391 stress 0.01528161 
    ## Run 392 stress 0.01719159 
    ## Run 393 stress 0.0156224 
    ## Run 394 stress 0.01649802 
    ## Run 395 stress 0.0123981 
    ## Run 396 stress 0.01557895 
    ## Run 397 stress 0.01547403 
    ## Run 398 stress 0.01499248 
    ## Run 399 stress 0.01558323 
    ## Run 400 stress 0.01499301 
    ## Run 401 stress 0.01120678 
    ## ... Procrustes: rmse 0.005199221  max resid 0.01171463 
    ## Run 402 stress 0.01115168 
    ## ... Procrustes: rmse 0.001182382  max resid 0.002668519 
    ## ... Similar to previous best
    ## Run 403 stress 0.0152818 
    ## Run 404 stress 0.01515031 
    ## Run 405 stress 0.01121494 
    ## ... Procrustes: rmse 0.005582393  max resid 0.01256731 
    ## Run 406 stress 0.01123469 
    ## ... Procrustes: rmse 0.006519026  max resid 0.01467776 
    ## Run 407 stress 0.01532595 
    ## Run 408 stress 0.01203768 
    ## Run 409 stress 0.01536131 
    ## Run 410 stress 0.015157 
    ## Run 411 stress 0.01557943 
    ## Run 412 stress 0.01556311 
    ## Run 413 stress 0.01523271 
    ## Run 414 stress 0.01556866 
    ## Run 415 stress 0.01557933 
    ## Run 416 stress 0.01583077 
    ## Run 417 stress 0.01542796 
    ## Run 418 stress 0.01685457 
    ## Run 419 stress 0.0151493 
    ## Run 420 stress 0.01686769 
    ## Run 421 stress 0.01533789 
    ## Run 422 stress 0.01115182 
    ## ... Procrustes: rmse 0.001240656  max resid 0.002800583 
    ## ... Similar to previous best
    ## Run 423 stress 0.01558331 
    ## Run 424 stress 0.01137252 
    ## ... Procrustes: rmse 0.01000178  max resid 0.02229166 
    ## Run 425 stress 0.01115181 
    ## ... Procrustes: rmse 0.0001771813  max resid 0.0004002592 
    ## ... Similar to previous best
    ## Run 426 stress 0.2940868 
    ## Run 427 stress 0.0153377 
    ## Run 428 stress 0.01542822 
    ## Run 429 stress 0.01515734 
    ## Run 430 stress 0.0155836 
    ## Run 431 stress 0.01116545 
    ## ... Procrustes: rmse 0.002402846  max resid 0.005420521 
    ## ... Similar to previous best
    ## Run 432 stress 0.01115334 
    ## ... Procrustes: rmse 0.0006450164  max resid 0.001456423 
    ## ... Similar to previous best
    ## Run 433 stress 0.01558325 
    ## Run 434 stress 0.01115168 
    ## ... Procrustes: rmse 0.001179291  max resid 0.002661386 
    ## ... Similar to previous best
    ## Run 435 stress 0.014986 
    ## Run 436 stress 0.01115529 
    ## ... Procrustes: rmse 0.001059996  max resid 0.002391193 
    ## ... Similar to previous best
    ## Run 437 stress 0.01147801 
    ## ... Procrustes: rmse 0.01144272  max resid 0.02500969 
    ## Run 438 stress 0.01493154 
    ## Run 439 stress 0.01547428 
    ## Run 440 stress 0.0151571 
    ## Run 441 stress 0.01542785 
    ## Run 442 stress 0.01543121 
    ## Run 443 stress 0.01499231 
    ## Run 444 stress 0.01517658 
    ## Run 445 stress 0.01498598 
    ## Run 446 stress 0.01115165 
    ## ... Procrustes: rmse 0.000108212  max resid 0.0002444346 
    ## ... Similar to previous best
    ## Run 447 stress 0.01600065 
    ## Run 448 stress 0.01571537 
    ## Run 449 stress 0.01123908 
    ## ... Procrustes: rmse 0.00659245  max resid 0.01483708 
    ## Run 450 stress 0.01567624 
    ## Run 451 stress 0.01115417 
    ## ... Procrustes: rmse 0.0008492006  max resid 0.001916425 
    ## ... Similar to previous best
    ## Run 452 stress 0.01515711 
    ## Run 453 stress 0.01557926 
    ## Run 454 stress 0.0166379 
    ## Run 455 stress 0.01498603 
    ## Run 456 stress 0.01493154 
    ## Run 457 stress 0.01116114 
    ## ... Procrustes: rmse 0.001908151  max resid 0.004304197 
    ## ... Similar to previous best
    ## Run 458 stress 0.01514914 
    ## Run 459 stress 0.01117846 
    ## ... Procrustes: rmse 0.003425702  max resid 0.007720647 
    ## ... Similar to previous best
    ## Run 460 stress 0.01515713 
    ## Run 461 stress 0.01549919 
    ## Run 462 stress 0.01533781 
    ## Run 463 stress 0.01582012 
    ## Run 464 stress 0.01581988 
    ## Run 465 stress 0.01115183 
    ## ... Procrustes: rmse 0.001245773  max resid 0.002812196 
    ## ... Similar to previous best
    ## Run 466 stress 0.01549943 
    ## Run 467 stress 0.01655638 
    ## Run 468 stress 0.01115201 
    ## ... Procrustes: rmse 0.0002530614  max resid 0.000571613 
    ## ... Similar to previous best
    ## Run 469 stress 0.0111945 
    ## ... Procrustes: rmse 0.004548276  max resid 0.01025265 
    ## Run 470 stress 0.01115177 
    ## ... Procrustes: rmse 0.001221294  max resid 0.002756659 
    ## ... Similar to previous best
    ## Run 471 stress 0.01514998 
    ## Run 472 stress 0.01170254 
    ## Run 473 stress 0.01525797 
    ## Run 474 stress 0.01545299 
    ## Run 475 stress 0.01536098 
    ## Run 476 stress 0.3040633 
    ## Run 477 stress 0.01498601 
    ## Run 478 stress 0.01116881 
    ## ... Procrustes: rmse 0.002718184  max resid 0.006132132 
    ## ... Similar to previous best
    ## Run 479 stress 0.01584383 
    ## Run 480 stress 0.0157683 
    ## Run 481 stress 0.01115242 
    ## ... Procrustes: rmse 0.0003997912  max resid 0.0009023774 
    ## ... Similar to previous best
    ## Run 482 stress 0.01499251 
    ## Run 483 stress 0.01529108 
    ## Run 484 stress 0.01498581 
    ## Run 485 stress 0.01543922 
    ## Run 486 stress 0.0112922 
    ## ... Procrustes: rmse 0.00864128  max resid 0.01945204 
    ## Run 487 stress 0.01560269 
    ## Run 488 stress 0.01581957 
    ## Run 489 stress 0.01525743 
    ## Run 490 stress 0.01536183 
    ## Run 491 stress 0.01499498 
    ## Run 492 stress 0.01545298 
    ## Run 493 stress 0.01390287 
    ## Run 494 stress 0.01545297 
    ## Run 495 stress 0.01567246 
    ## Run 496 stress 0.01557165 
    ## Run 497 stress 0.01115328 
    ## ... Procrustes: rmse 0.0006297803  max resid 0.00142065 
    ## ... Similar to previous best
    ## Run 498 stress 0.01568109 
    ## Run 499 stress 0.01562182 
    ## Run 500 stress 0.01439298 
    ## *** Best solution repeated 76 times

``` r
# Mixed and stratified lakes
PD_beta_MS_NMDS <- metaMDS(PD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.05036792 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002043154  max resid 0.0005214254 
    ## ... Similar to previous best
    ## Run 2 stress 0.05051317 
    ## ... Procrustes: rmse 0.03567173  max resid 0.1218823 
    ## Run 3 stress 0.05928218 
    ## Run 4 stress 0.05036807 
    ## ... Procrustes: rmse 0.0001472172  max resid 0.0003610559 
    ## ... Similar to previous best
    ## Run 5 stress 0.05036799 
    ## ... Procrustes: rmse 9.799864e-05  max resid 0.0002548983 
    ## ... Similar to previous best
    ## Run 6 stress 0.05051293 
    ## ... Procrustes: rmse 0.03565913  max resid 0.1217249 
    ## Run 7 stress 0.05036808 
    ## ... Procrustes: rmse 0.0001546594  max resid 0.0003778773 
    ## ... Similar to previous best
    ## Run 8 stress 0.05051333 
    ## ... Procrustes: rmse 0.03570209  max resid 0.1218172 
    ## Run 9 stress 0.05036802 
    ## ... Procrustes: rmse 0.000114781  max resid 0.0002735839 
    ## ... Similar to previous best
    ## Run 10 stress 0.05051306 
    ## ... Procrustes: rmse 0.03565124  max resid 0.1216864 
    ## Run 11 stress 0.05036795 
    ## ... Procrustes: rmse 5.394576e-05  max resid 0.0001339348 
    ## ... Similar to previous best
    ## Run 12 stress 0.05403599 
    ## Run 13 stress 0.0517106 
    ## Run 14 stress 0.05602044 
    ## Run 15 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001149541  max resid 0.0002858371 
    ## ... Similar to previous best
    ## Run 16 stress 0.05036801 
    ## ... Procrustes: rmse 0.00010736  max resid 0.0002374046 
    ## ... Similar to previous best
    ## Run 17 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 5.308515e-05  max resid 0.0001152688 
    ## ... Similar to previous best
    ## Run 18 stress 0.05171058 
    ## Run 19 stress 0.0582995 
    ## Run 20 stress 0.05051305 
    ## ... Procrustes: rmse 0.03565311  max resid 0.1216981 
    ## Run 21 stress 0.05171059 
    ## Run 22 stress 0.0517106 
    ## Run 23 stress 0.05968625 
    ## Run 24 stress 0.05466702 
    ## Run 25 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001072912  max resid 0.000265875 
    ## ... Similar to previous best
    ## Run 26 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001004583  max resid 0.0002157251 
    ## ... Similar to previous best
    ## Run 27 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001398471  max resid 0.0003114657 
    ## ... Similar to previous best
    ## Run 28 stress 0.05171059 
    ## Run 29 stress 0.06136687 
    ## Run 30 stress 0.05968628 
    ## Run 31 stress 0.05036793 
    ## ... Procrustes: rmse 6.670787e-05  max resid 0.0001532234 
    ## ... Similar to previous best
    ## Run 32 stress 0.05602053 
    ## Run 33 stress 0.05602048 
    ## Run 34 stress 0.05036793 
    ## ... Procrustes: rmse 6.289261e-05  max resid 0.0001397414 
    ## ... Similar to previous best
    ## Run 35 stress 0.050368 
    ## ... Procrustes: rmse 0.0001317299  max resid 0.0002909661 
    ## ... Similar to previous best
    ## Run 36 stress 0.06499236 
    ## Run 37 stress 0.05036792 
    ## ... Procrustes: rmse 6.748001e-05  max resid 0.0001454213 
    ## ... Similar to previous best
    ## Run 38 stress 0.05403619 
    ## Run 39 stress 0.05036795 
    ## ... Procrustes: rmse 9.630984e-05  max resid 0.0002607282 
    ## ... Similar to previous best
    ## Run 40 stress 0.054667 
    ## Run 41 stress 0.06904478 
    ## Run 42 stress 0.05337593 
    ## Run 43 stress 0.05036794 
    ## ... Procrustes: rmse 6.163163e-05  max resid 0.0001611584 
    ## ... Similar to previous best
    ## Run 44 stress 0.06476214 
    ## Run 45 stress 0.06111767 
    ## Run 46 stress 0.06027447 
    ## Run 47 stress 0.06826717 
    ## Run 48 stress 0.06231298 
    ## Run 49 stress 0.05036794 
    ## ... Procrustes: rmse 7.825338e-05  max resid 0.0002086529 
    ## ... Similar to previous best
    ## Run 50 stress 0.0517106 
    ## Run 51 stress 0.05171061 
    ## Run 52 stress 0.0517106 
    ## Run 53 stress 0.05036794 
    ## ... Procrustes: rmse 6.95383e-05  max resid 0.0001854423 
    ## ... Similar to previous best
    ## Run 54 stress 0.05171063 
    ## Run 55 stress 0.0560207 
    ## Run 56 stress 0.05036792 
    ## ... Procrustes: rmse 6.098128e-05  max resid 0.0001367821 
    ## ... Similar to previous best
    ## Run 57 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001564618  max resid 0.0004580144 
    ## ... Similar to previous best
    ## Run 58 stress 0.05928228 
    ## Run 59 stress 0.050368 
    ## ... Procrustes: rmse 0.0001376496  max resid 0.0003715162 
    ## ... Similar to previous best
    ## Run 60 stress 0.06478539 
    ## Run 61 stress 0.05602062 
    ## Run 62 stress 0.06609077 
    ## Run 63 stress 0.05466705 
    ## Run 64 stress 0.05171063 
    ## Run 65 stress 0.0611177 
    ## Run 66 stress 0.0592826 
    ## Run 67 stress 0.05171059 
    ## Run 68 stress 0.05403622 
    ## Run 69 stress 0.05051286 
    ## ... Procrustes: rmse 0.03569513  max resid 0.1219406 
    ## Run 70 stress 0.05337577 
    ## Run 71 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001347907  max resid 0.0003191314 
    ## ... Similar to previous best
    ## Run 72 stress 0.05036808 
    ## ... Procrustes: rmse 0.000196557  max resid 0.0005004856 
    ## ... Similar to previous best
    ## Run 73 stress 0.06231297 
    ## Run 74 stress 0.05036794 
    ## ... Procrustes: rmse 9.621511e-05  max resid 0.0002219995 
    ## ... Similar to previous best
    ## Run 75 stress 0.05979234 
    ## Run 76 stress 0.05900152 
    ## Run 77 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001055496  max resid 0.0002434483 
    ## ... Similar to previous best
    ## Run 78 stress 0.05602062 
    ## Run 79 stress 0.05036802 
    ## ... Procrustes: rmse 0.000155723  max resid 0.0004305444 
    ## ... Similar to previous best
    ## Run 80 stress 0.05051321 
    ## ... Procrustes: rmse 0.03566563  max resid 0.121724 
    ## Run 81 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001483865  max resid 0.0004228579 
    ## ... Similar to previous best
    ## Run 82 stress 0.05171062 
    ## Run 83 stress 0.05036794 
    ## ... Procrustes: rmse 8.948974e-05  max resid 0.0002009884 
    ## ... Similar to previous best
    ## Run 84 stress 0.05051327 
    ## ... Procrustes: rmse 0.03573141  max resid 0.1220739 
    ## Run 85 stress 0.05036792 
    ## ... Procrustes: rmse 4.836397e-05  max resid 0.000141161 
    ## ... Similar to previous best
    ## Run 86 stress 0.0560205 
    ## Run 87 stress 0.05602043 
    ## Run 88 stress 0.05051297 
    ## ... Procrustes: rmse 0.03569435  max resid 0.1219456 
    ## Run 89 stress 0.05051329 
    ## ... Procrustes: rmse 0.03570613  max resid 0.1218357 
    ## Run 90 stress 0.05829946 
    ## Run 91 stress 0.05051325 
    ## ... Procrustes: rmse 0.035772  max resid 0.1220436 
    ## Run 92 stress 0.05171058 
    ## Run 93 stress 0.06027437 
    ## Run 94 stress 0.05466701 
    ## Run 95 stress 0.05171059 
    ## Run 96 stress 0.05051324 
    ## ... Procrustes: rmse 0.03575439  max resid 0.121989 
    ## Run 97 stress 0.05036792 
    ## ... Procrustes: rmse 2.334486e-05  max resid 6.069079e-05 
    ## ... Similar to previous best
    ## Run 98 stress 0.05051324 
    ## ... Procrustes: rmse 0.03572347  max resid 0.1220486 
    ## Run 99 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001385039  max resid 0.0003345034 
    ## ... Similar to previous best
    ## Run 100 stress 0.07221112 
    ## Run 101 stress 0.06309848 
    ## Run 102 stress 0.05337579 
    ## Run 103 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001020806  max resid 0.0002777923 
    ## ... Similar to previous best
    ## Run 104 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001702293  max resid 0.0004172753 
    ## ... Similar to previous best
    ## Run 105 stress 0.0611178 
    ## Run 106 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001136338  max resid 0.0002372422 
    ## ... Similar to previous best
    ## Run 107 stress 0.05051336 
    ## ... Procrustes: rmse 0.03570243  max resid 0.1219904 
    ## Run 108 stress 0.06231293 
    ## Run 109 stress 0.06851912 
    ## Run 110 stress 0.05171062 
    ## Run 111 stress 0.05928225 
    ## Run 112 stress 0.0630986 
    ## Run 113 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001052736  max resid 0.0003091444 
    ## ... Similar to previous best
    ## Run 114 stress 0.05036792 
    ## ... Procrustes: rmse 4.914698e-05  max resid 0.0001076386 
    ## ... Similar to previous best
    ## Run 115 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001793201  max resid 0.0005069275 
    ## ... Similar to previous best
    ## Run 116 stress 0.05968626 
    ## Run 117 stress 0.06231297 
    ## Run 118 stress 0.0505131 
    ## ... Procrustes: rmse 0.03570032  max resid 0.1218384 
    ## Run 119 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001601862  max resid 0.0003931594 
    ## ... Similar to previous best
    ## Run 120 stress 0.07104253 
    ## Run 121 stress 0.06362799 
    ## Run 122 stress 0.05337578 
    ## Run 123 stress 0.05051317 
    ## ... Procrustes: rmse 0.03568423  max resid 0.1217822 
    ## Run 124 stress 0.05602046 
    ## Run 125 stress 0.05968631 
    ## Run 126 stress 0.05466701 
    ## Run 127 stress 0.05051311 
    ## ... Procrustes: rmse 0.03572981  max resid 0.1220619 
    ## Run 128 stress 0.06309847 
    ## Run 129 stress 0.05036794 
    ## ... Procrustes: rmse 8.895772e-05  max resid 0.0002031379 
    ## ... Similar to previous best
    ## Run 130 stress 0.05337579 
    ## Run 131 stress 0.05051303 
    ## ... Procrustes: rmse 0.03570694  max resid 0.1219886 
    ## Run 132 stress 0.05829928 
    ## Run 133 stress 0.05602053 
    ## Run 134 stress 0.05602054 
    ## Run 135 stress 0.0505134 
    ## ... Procrustes: rmse 0.03568583  max resid 0.1217651 
    ## Run 136 stress 0.0505132 
    ## ... Procrustes: rmse 0.03571928  max resid 0.1220353 
    ## Run 137 stress 0.05337577 
    ## Run 138 stress 0.05051325 
    ## ... Procrustes: rmse 0.03573219  max resid 0.1220756 
    ## Run 139 stress 0.05337584 
    ## Run 140 stress 0.06362798 
    ## Run 141 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001253961  max resid 0.0003034459 
    ## ... Similar to previous best
    ## Run 142 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001041964  max resid 0.0002439821 
    ## ... Similar to previous best
    ## Run 143 stress 0.06885079 
    ## Run 144 stress 0.05171064 
    ## Run 145 stress 0.05036809 
    ## ... Procrustes: rmse 0.0002167198  max resid 0.000521603 
    ## ... Similar to previous best
    ## Run 146 stress 0.06027428 
    ## Run 147 stress 0.05036794 
    ## ... Procrustes: rmse 6.026667e-05  max resid 0.0001619287 
    ## ... Similar to previous best
    ## Run 148 stress 0.0505133 
    ## ... Procrustes: rmse 0.03571109  max resid 0.1218511 
    ## Run 149 stress 0.06027436 
    ## Run 150 stress 0.05829935 
    ## Run 151 stress 0.050368 
    ## ... Procrustes: rmse 0.0001222323  max resid 0.0003167843 
    ## ... Similar to previous best
    ## Run 152 stress 0.05337578 
    ## Run 153 stress 0.05466704 
    ## Run 154 stress 0.05036794 
    ## ... Procrustes: rmse 9.423028e-05  max resid 0.0002213158 
    ## ... Similar to previous best
    ## Run 155 stress 0.0673021 
    ## Run 156 stress 0.05051341 
    ## ... Procrustes: rmse 0.03570534  max resid 0.1220015 
    ## Run 157 stress 0.06476213 
    ## Run 158 stress 0.06309855 
    ## Run 159 stress 0.05337583 
    ## Run 160 stress 0.05968627 
    ## Run 161 stress 0.05602046 
    ## Run 162 stress 0.05968626 
    ## Run 163 stress 0.05403558 
    ## Run 164 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001725103  max resid 0.0004051365 
    ## ... Similar to previous best
    ## Run 165 stress 0.05051348 
    ## ... Procrustes: rmse 0.03574834  max resid 0.1221322 
    ## Run 166 stress 0.05968638 
    ## Run 167 stress 0.05602062 
    ## Run 168 stress 0.05171059 
    ## Run 169 stress 0.05337584 
    ## Run 170 stress 0.05036792 
    ## ... Procrustes: rmse 4.364734e-05  max resid 9.518488e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.05051298 
    ## ... Procrustes: rmse 0.03569542  max resid 0.1218386 
    ## Run 172 stress 0.05051297 
    ## ... Procrustes: rmse 0.03557545  max resid 0.1215707 
    ## Run 173 stress 0.0517106 
    ## Run 174 stress 0.05928228 
    ## Run 175 stress 0.05036791 
    ## ... Procrustes: rmse 2.17711e-05  max resid 4.470759e-05 
    ## ... Similar to previous best
    ## Run 176 stress 0.05171058 
    ## Run 177 stress 0.05171059 
    ## Run 178 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001476385  max resid 0.0003532077 
    ## ... Similar to previous best
    ## Run 179 stress 0.05036793 
    ## ... Procrustes: rmse 7.768589e-05  max resid 0.000180776 
    ## ... Similar to previous best
    ## Run 180 stress 0.05171059 
    ## Run 181 stress 0.05036794 
    ## ... Procrustes: rmse 8.700962e-05  max resid 0.0001962227 
    ## ... Similar to previous best
    ## Run 182 stress 0.05602052 
    ## Run 183 stress 0.06896879 
    ## Run 184 stress 0.05051274 
    ## ... Procrustes: rmse 0.03569034  max resid 0.1219088 
    ## Run 185 stress 0.05036808 
    ## ... Procrustes: rmse 0.0002056679  max resid 0.0005095372 
    ## ... Similar to previous best
    ## Run 186 stress 0.06476235 
    ## Run 187 stress 0.05928231 
    ## Run 188 stress 0.05900148 
    ## Run 189 stress 0.05036795 
    ## ... Procrustes: rmse 8.10755e-05  max resid 0.0002209852 
    ## ... Similar to previous best
    ## Run 190 stress 0.05036793 
    ## ... Procrustes: rmse 6.871364e-05  max resid 0.0001946615 
    ## ... Similar to previous best
    ## Run 191 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001436661  max resid 0.0004201516 
    ## ... Similar to previous best
    ## Run 192 stress 0.0647623 
    ## Run 193 stress 0.0582993 
    ## Run 194 stress 0.05171058 
    ## Run 195 stress 0.07015953 
    ## Run 196 stress 0.06231308 
    ## Run 197 stress 0.05829929 
    ## Run 198 stress 0.05051338 
    ## ... Procrustes: rmse 0.03564407  max resid 0.1216449 
    ## Run 199 stress 0.050513 
    ## ... Procrustes: rmse 0.03570228  max resid 0.121973 
    ## Run 200 stress 0.06231298 
    ## Run 201 stress 0.06231291 
    ## Run 202 stress 0.05051278 
    ## ... Procrustes: rmse 0.03562349  max resid 0.1217076 
    ## Run 203 stress 0.05051334 
    ## ... Procrustes: rmse 0.03570171  max resid 0.1218178 
    ## Run 204 stress 0.06640347 
    ## Run 205 stress 0.05171059 
    ## Run 206 stress 0.05403607 
    ## Run 207 stress 0.05171059 
    ## Run 208 stress 0.06478541 
    ## Run 209 stress 0.05968632 
    ## Run 210 stress 0.05051337 
    ## ... Procrustes: rmse 0.03570232  max resid 0.1219906 
    ## Run 211 stress 0.05829938 
    ## Run 212 stress 0.05602054 
    ## Run 213 stress 0.05602056 
    ## Run 214 stress 0.050368 
    ## ... Procrustes: rmse 0.0001290807  max resid 0.0003480665 
    ## ... Similar to previous best
    ## Run 215 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001294829  max resid 0.0003200314 
    ## ... Similar to previous best
    ## Run 216 stress 0.05928219 
    ## Run 217 stress 0.06478524 
    ## Run 218 stress 0.06309856 
    ## Run 219 stress 0.05036793 
    ## ... Procrustes: rmse 5.787848e-05  max resid 0.0001692259 
    ## ... Similar to previous best
    ## Run 220 stress 0.05968625 
    ## Run 221 stress 0.06896923 
    ## Run 222 stress 0.05051366 
    ## ... Procrustes: rmse 0.03570845  max resid 0.1220147 
    ## Run 223 stress 0.05602047 
    ## Run 224 stress 0.05051295 
    ## ... Procrustes: rmse 0.03569562  max resid 0.1218422 
    ## Run 225 stress 0.06499243 
    ## Run 226 stress 0.05051327 
    ## ... Procrustes: rmse 0.0358098  max resid 0.1221624 
    ## Run 227 stress 0.06655674 
    ## Run 228 stress 0.05829929 
    ## Run 229 stress 0.05051304 
    ## ... Procrustes: rmse 0.03568917  max resid 0.1218131 
    ## Run 230 stress 0.0596864 
    ## Run 231 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001777249  max resid 0.000440521 
    ## ... Similar to previous best
    ## Run 232 stress 0.05051334 
    ## ... Procrustes: rmse 0.03570322  max resid 0.1219926 
    ## Run 233 stress 0.0517106 
    ## Run 234 stress 0.06027433 
    ## Run 235 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001599232  max resid 0.0003848054 
    ## ... Similar to previous best
    ## Run 236 stress 0.05036791 
    ## ... Procrustes: rmse 9.579261e-06  max resid 2.59926e-05 
    ## ... Similar to previous best
    ## Run 237 stress 0.05829939 
    ## Run 238 stress 0.050368 
    ## ... Procrustes: rmse 0.0001090891  max resid 0.0001933766 
    ## ... Similar to previous best
    ## Run 239 stress 0.05602053 
    ## Run 240 stress 0.05036794 
    ## ... Procrustes: rmse 1.82828e-05  max resid 3.582363e-05 
    ## ... Similar to previous best
    ## Run 241 stress 0.06231287 
    ## Run 242 stress 0.06640335 
    ## Run 243 stress 0.05051365 
    ## ... Procrustes: rmse 0.03574806  max resid 0.122133 
    ## Run 244 stress 0.05337589 
    ## Run 245 stress 0.05403598 
    ## Run 246 stress 0.05829929 
    ## Run 247 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001772239  max resid 0.0004300818 
    ## ... Similar to previous best
    ## Run 248 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001639423  max resid 0.0003801604 
    ## ... Similar to previous best
    ## Run 249 stress 0.05466705 
    ## Run 250 stress 0.06136677 
    ## Run 251 stress 0.06478536 
    ## Run 252 stress 0.05337577 
    ## Run 253 stress 0.05829928 
    ## Run 254 stress 0.05337599 
    ## Run 255 stress 0.05829931 
    ## Run 256 stress 0.05829931 
    ## Run 257 stress 0.0630985 
    ## Run 258 stress 0.05602055 
    ## Run 259 stress 0.0596863 
    ## Run 260 stress 0.05466704 
    ## Run 261 stress 0.05337578 
    ## Run 262 stress 0.06826716 
    ## Run 263 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001097638  max resid 0.0002185635 
    ## ... Similar to previous best
    ## Run 264 stress 0.05466704 
    ## Run 265 stress 0.05602066 
    ## Run 266 stress 0.06476228 
    ## Run 267 stress 0.05036793 
    ## ... Procrustes: rmse 4.919692e-05  max resid 0.000130767 
    ## ... Similar to previous best
    ## Run 268 stress 0.05051315 
    ## ... Procrustes: rmse 0.03569176  max resid 0.1218074 
    ## Run 269 stress 0.05979232 
    ## Run 270 stress 0.05051321 
    ## ... Procrustes: rmse 0.03570846  max resid 0.1218513 
    ## Run 271 stress 0.05829932 
    ## Run 272 stress 0.05051326 
    ## ... Procrustes: rmse 0.03571767  max resid 0.1220324 
    ## Run 273 stress 0.05051306 
    ## ... Procrustes: rmse 0.03572131  max resid 0.122034 
    ## Run 274 stress 0.05466709 
    ## Run 275 stress 0.0602743 
    ## Run 276 stress 0.0505129 
    ## ... Procrustes: rmse 0.03569741  max resid 0.1219511 
    ## Run 277 stress 0.0613668 
    ## Run 278 stress 0.05036794 
    ## ... Procrustes: rmse 7.452825e-05  max resid 0.0001934898 
    ## ... Similar to previous best
    ## Run 279 stress 0.06609059 
    ## Run 280 stress 0.05171059 
    ## Run 281 stress 0.0533759 
    ## Run 282 stress 0.05171059 
    ## Run 283 stress 0.05928245 
    ## Run 284 stress 0.0505131 
    ## ... Procrustes: rmse 0.03570709  max resid 0.1219938 
    ## Run 285 stress 0.05051349 
    ## ... Procrustes: rmse 0.03573732  max resid 0.1220988 
    ## Run 286 stress 0.05051339 
    ## ... Procrustes: rmse 0.03574255  max resid 0.1221115 
    ## Run 287 stress 0.3592584 
    ## Run 288 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001262729  max resid 0.000338933 
    ## ... Similar to previous best
    ## Run 289 stress 0.05036795 
    ## ... Procrustes: rmse 9.740816e-05  max resid 0.0002654284 
    ## ... Similar to previous best
    ## Run 290 stress 0.05036801 
    ## ... Procrustes: rmse 7.786652e-05  max resid 0.0001961242 
    ## ... Similar to previous best
    ## Run 291 stress 0.05403599 
    ## Run 292 stress 0.06478512 
    ## Run 293 stress 0.06478514 
    ## Run 294 stress 0.05900147 
    ## Run 295 stress 0.05051296 
    ## ... Procrustes: rmse 0.03569405  max resid 0.1219458 
    ## Run 296 stress 0.05036802 
    ## ... Procrustes: rmse 0.000147357  max resid 0.0004030999 
    ## ... Similar to previous best
    ## Run 297 stress 0.0505131 
    ## ... Procrustes: rmse 0.03571813  max resid 0.121894 
    ## Run 298 stress 0.05602048 
    ## Run 299 stress 0.05928217 
    ## Run 300 stress 0.05171059 
    ## Run 301 stress 0.05466701 
    ## Run 302 stress 0.05900148 
    ## Run 303 stress 0.05829938 
    ## Run 304 stress 0.05171059 
    ## Run 305 stress 0.06504236 
    ## Run 306 stress 0.06647183 
    ## Run 307 stress 0.05171059 
    ## Run 308 stress 0.0533759 
    ## Run 309 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001262377  max resid 0.0003710398 
    ## ... Similar to previous best
    ## Run 310 stress 0.05928212 
    ## Run 311 stress 0.06038952 
    ## Run 312 stress 0.05602064 
    ## Run 313 stress 0.05051284 
    ## ... Procrustes: rmse 0.03569359  max resid 0.1218526 
    ## Run 314 stress 0.0517106 
    ## Run 315 stress 0.05466702 
    ## Run 316 stress 0.06478537 
    ## Run 317 stress 0.05466707 
    ## Run 318 stress 0.05466707 
    ## Run 319 stress 0.05829932 
    ## Run 320 stress 0.05602055 
    ## Run 321 stress 0.05051345 
    ## ... Procrustes: rmse 0.03577805  max resid 0.1222184 
    ## Run 322 stress 0.06309849 
    ## Run 323 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001739421  max resid 0.0004203977 
    ## ... Similar to previous best
    ## Run 324 stress 0.05051342 
    ## ... Procrustes: rmse 0.03572738  max resid 0.122067 
    ## Run 325 stress 0.05403615 
    ## Run 326 stress 0.06647186 
    ## Run 327 stress 0.0592822 
    ## Run 328 stress 0.05337594 
    ## Run 329 stress 0.05036792 
    ## ... Procrustes: rmse 4.533258e-05  max resid 9.739142e-05 
    ## ... Similar to previous best
    ## Run 330 stress 0.05968635 
    ## Run 331 stress 0.06309843 
    ## Run 332 stress 0.05968627 
    ## Run 333 stress 0.06499234 
    ## Run 334 stress 0.05466702 
    ## Run 335 stress 0.05051311 
    ## ... Procrustes: rmse 0.03567445  max resid 0.1217632 
    ## Run 336 stress 0.05036791 
    ## ... Procrustes: rmse 2.870608e-05  max resid 6.247199e-05 
    ## ... Similar to previous best
    ## Run 337 stress 0.05602078 
    ## Run 338 stress 0.06499232 
    ## Run 339 stress 0.05051336 
    ## ... Procrustes: rmse 0.03570532  max resid 0.1218279 
    ## Run 340 stress 0.05036795 
    ## ... Procrustes: rmse 8.838431e-05  max resid 0.0002594258 
    ## ... Similar to previous best
    ## Run 341 stress 0.06476246 
    ## Run 342 stress 0.05602065 
    ## Run 343 stress 0.06136683 
    ## Run 344 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001637072  max resid 0.0004319861 
    ## ... Similar to previous best
    ## Run 345 stress 0.05900153 
    ## Run 346 stress 0.06362779 
    ## Run 347 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001057145  max resid 0.0002463827 
    ## ... Similar to previous best
    ## Run 348 stress 0.05171066 
    ## Run 349 stress 0.0649924 
    ## Run 350 stress 0.06231288 
    ## Run 351 stress 0.05036798 
    ## ... Procrustes: rmse 0.000122804  max resid 0.0003289265 
    ## ... Similar to previous best
    ## Run 352 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001603254  max resid 0.0003778668 
    ## ... Similar to previous best
    ## Run 353 stress 0.05036793 
    ## ... Procrustes: rmse 6.937258e-05  max resid 0.0001574767 
    ## ... Similar to previous best
    ## Run 354 stress 0.05403599 
    ## Run 355 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001026756  max resid 0.0002410186 
    ## ... Similar to previous best
    ## Run 356 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001148759  max resid 0.0002725322 
    ## ... Similar to previous best
    ## Run 357 stress 0.05036792 
    ## ... Procrustes: rmse 3.88301e-05  max resid 0.0001062358 
    ## ... Similar to previous best
    ## Run 358 stress 0.06504219 
    ## Run 359 stress 0.05337588 
    ## Run 360 stress 0.06896914 
    ## Run 361 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001600035  max resid 0.0003889921 
    ## ... Similar to previous best
    ## Run 362 stress 0.05602045 
    ## Run 363 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001695404  max resid 0.0003910283 
    ## ... Similar to previous best
    ## Run 364 stress 0.05051287 
    ## ... Procrustes: rmse 0.03558611  max resid 0.1215427 
    ## Run 365 stress 0.05051274 
    ## ... Procrustes: rmse 0.03575363  max resid 0.1220703 
    ## Run 366 stress 0.05036801 
    ## ... Procrustes: rmse 0.000164839  max resid 0.000396395 
    ## ... Similar to previous best
    ## Run 367 stress 0.05036794 
    ## ... Procrustes: rmse 9.289215e-05  max resid 0.0002160349 
    ## ... Similar to previous best
    ## Run 368 stress 0.0505132 
    ## ... Procrustes: rmse 0.03567198  max resid 0.1217442 
    ## Run 369 stress 0.05466701 
    ## Run 370 stress 0.05051294 
    ## ... Procrustes: rmse 0.03568812  max resid 0.1218213 
    ## Run 371 stress 0.05337577 
    ## Run 372 stress 0.06640332 
    ## Run 373 stress 0.05036793 
    ## ... Procrustes: rmse 6.285283e-05  max resid 0.000183427 
    ## ... Similar to previous best
    ## Run 374 stress 0.05051297 
    ## ... Procrustes: rmse 0.03567907  max resid 0.12179 
    ## Run 375 stress 0.05051295 
    ## ... Procrustes: rmse 0.03570429  max resid 0.1219752 
    ## Run 376 stress 0.05036794 
    ## ... Procrustes: rmse 7.083713e-05  max resid 0.0001555404 
    ## ... Similar to previous best
    ## Run 377 stress 0.05171061 
    ## Run 378 stress 0.06038957 
    ## Run 379 stress 0.05928217 
    ## Run 380 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001129159  max resid 0.0002974039 
    ## ... Similar to previous best
    ## Run 381 stress 0.05403608 
    ## Run 382 stress 0.3444382 
    ## Run 383 stress 0.06476216 
    ## Run 384 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001296377  max resid 0.0003492177 
    ## ... Similar to previous best
    ## Run 385 stress 0.05036793 
    ## ... Procrustes: rmse 5.607675e-05  max resid 0.0001399957 
    ## ... Similar to previous best
    ## Run 386 stress 0.06609072 
    ## Run 387 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001377448  max resid 0.0003293374 
    ## ... Similar to previous best
    ## Run 388 stress 0.07462514 
    ## Run 389 stress 0.05171063 
    ## Run 390 stress 0.05051301 
    ## ... Procrustes: rmse 0.03570281  max resid 0.1219757 
    ## Run 391 stress 0.06231295 
    ## Run 392 stress 0.05928224 
    ## Run 393 stress 0.05051297 
    ## ... Procrustes: rmse 0.03572259  max resid 0.1219219 
    ## Run 394 stress 0.05968636 
    ## Run 395 stress 0.05051315 
    ## ... Procrustes: rmse 0.03570402  max resid 0.1219813 
    ## Run 396 stress 0.05051318 
    ## ... Procrustes: rmse 0.03574086  max resid 0.1220987 
    ## Run 397 stress 0.05466701 
    ## Run 398 stress 0.05466701 
    ## Run 399 stress 0.06478515 
    ## Run 400 stress 0.05602068 
    ## Run 401 stress 0.0694471 
    ## Run 402 stress 0.0596863 
    ## Run 403 stress 0.05403604 
    ## Run 404 stress 0.05829931 
    ## Run 405 stress 0.05171063 
    ## Run 406 stress 0.05036796 
    ## ... Procrustes: rmse 9.265828e-05  max resid 0.0002690518 
    ## ... Similar to previous best
    ## Run 407 stress 0.05171059 
    ## Run 408 stress 0.06478513 
    ## Run 409 stress 0.05171058 
    ## Run 410 stress 0.05171066 
    ## Run 411 stress 0.05051303 
    ## ... Procrustes: rmse 0.03577424  max resid 0.122187 
    ## Run 412 stress 0.05036791 
    ## ... Procrustes: rmse 3.045236e-05  max resid 6.117827e-05 
    ## ... Similar to previous best
    ## Run 413 stress 0.0517106 
    ## Run 414 stress 0.05051296 
    ## ... Procrustes: rmse 0.03568515  max resid 0.1218094 
    ## Run 415 stress 0.05602071 
    ## Run 416 stress 0.05928245 
    ## Run 417 stress 0.05036795 
    ## ... Procrustes: rmse 7.875723e-05  max resid 0.0002289233 
    ## ... Similar to previous best
    ## Run 418 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001257011  max resid 0.0003677572 
    ## ... Similar to previous best
    ## Run 419 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001000891  max resid 0.0002364334 
    ## ... Similar to previous best
    ## Run 420 stress 0.0596864 
    ## Run 421 stress 0.06409297 
    ## Run 422 stress 0.06499238 
    ## Run 423 stress 0.05051327 
    ## ... Procrustes: rmse 0.03575678  max resid 0.1221494 
    ## Run 424 stress 0.05051343 
    ## ... Procrustes: rmse 0.03573404  max resid 0.122088 
    ## Run 425 stress 0.05466703 
    ## Run 426 stress 0.05036799 
    ## ... Procrustes: rmse 0.000140758  max resid 0.0003349821 
    ## ... Similar to previous best
    ## Run 427 stress 0.06136676 
    ## Run 428 stress 0.05051352 
    ## ... Procrustes: rmse 0.03570536  max resid 0.1220045 
    ## Run 429 stress 0.05051351 
    ## ... Procrustes: rmse 0.0357034  max resid 0.1219958 
    ## Run 430 stress 0.05171059 
    ## Run 431 stress 0.05829938 
    ## Run 432 stress 0.0505128 
    ## ... Procrustes: rmse 0.03571478  max resid 0.1219257 
    ## Run 433 stress 0.05403598 
    ## Run 434 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001185557  max resid 0.0002836773 
    ## ... Similar to previous best
    ## Run 435 stress 0.050368 
    ## ... Procrustes: rmse 0.0001528945  max resid 0.000364027 
    ## ... Similar to previous best
    ## Run 436 stress 0.05602063 
    ## Run 437 stress 0.05928237 
    ## Run 438 stress 0.06655663 
    ## Run 439 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001590867  max resid 0.0003899316 
    ## ... Similar to previous best
    ## Run 440 stress 0.06136684 
    ## Run 441 stress 0.05171059 
    ## Run 442 stress 0.06231313 
    ## Run 443 stress 0.05171058 
    ## Run 444 stress 0.05036794 
    ## ... Procrustes: rmse 3.45162e-05  max resid 7.005055e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.05036794 
    ## ... Procrustes: rmse 8.538059e-05  max resid 0.0002286492 
    ## ... Similar to previous best
    ## Run 446 stress 0.05036799 
    ## ... Procrustes: rmse 0.000133042  max resid 0.0003174336 
    ## ... Similar to previous best
    ## Run 447 stress 0.05928225 
    ## Run 448 stress 0.0602743 
    ## Run 449 stress 0.06027434 
    ## Run 450 stress 0.0560206 
    ## Run 451 stress 0.06478512 
    ## Run 452 stress 0.06027436 
    ## Run 453 stress 0.05036807 
    ## ... Procrustes: rmse 0.0001926575  max resid 0.0005121027 
    ## ... Similar to previous best
    ## Run 454 stress 0.05051281 
    ## ... Procrustes: rmse 0.03560168  max resid 0.1216331 
    ## Run 455 stress 0.05829931 
    ## Run 456 stress 0.05968629 
    ## Run 457 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001178319  max resid 0.00031626 
    ## ... Similar to previous best
    ## Run 458 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001719197  max resid 0.0004285447 
    ## ... Similar to previous best
    ## Run 459 stress 0.05928236 
    ## Run 460 stress 0.06730241 
    ## Run 461 stress 0.05051357 
    ## ... Procrustes: rmse 0.03574642  max resid 0.1221271 
    ## Run 462 stress 0.05968626 
    ## Run 463 stress 0.05051288 
    ## ... Procrustes: rmse 0.03569599  max resid 0.1219448 
    ## Run 464 stress 0.05403593 
    ## Run 465 stress 0.05968636 
    ## Run 466 stress 0.05928223 
    ## Run 467 stress 0.06027427 
    ## Run 468 stress 0.06478536 
    ## Run 469 stress 0.05337577 
    ## Run 470 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001293791  max resid 0.0003692649 
    ## ... Similar to previous best
    ## Run 471 stress 0.06478524 
    ## Run 472 stress 0.05979229 
    ## Run 473 stress 0.05968629 
    ## Run 474 stress 0.05051354 
    ## ... Procrustes: rmse 0.03570522  max resid 0.1220047 
    ## Run 475 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001096341  max resid 0.000244294 
    ## ... Similar to previous best
    ## Run 476 stress 0.3662723 
    ## Run 477 stress 0.06309844 
    ## Run 478 stress 0.05337582 
    ## Run 479 stress 0.05466707 
    ## Run 480 stress 0.05171059 
    ## Run 481 stress 0.05829929 
    ## Run 482 stress 0.06730207 
    ## Run 483 stress 0.0533758 
    ## Run 484 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001887015  max resid 0.0004617346 
    ## ... Similar to previous best
    ## Run 485 stress 0.05928247 
    ## Run 486 stress 0.0505129 
    ## ... Procrustes: rmse 0.03562657  max resid 0.1217345 
    ## Run 487 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001090915  max resid 0.0002580214 
    ## ... Similar to previous best
    ## Run 488 stress 0.05968627 
    ## Run 489 stress 0.05051287 
    ## ... Procrustes: rmse 0.03562079  max resid 0.1216467 
    ## Run 490 stress 0.05051355 
    ## ... Procrustes: rmse 0.03565496  max resid 0.1216649 
    ## Run 491 stress 0.05403552 
    ## Run 492 stress 0.0623129 
    ## Run 493 stress 0.05829931 
    ## Run 494 stress 0.050368 
    ## ... Procrustes: rmse 0.0001423655  max resid 0.000371748 
    ## ... Similar to previous best
    ## Run 495 stress 0.06630526 
    ## Run 496 stress 0.06730224 
    ## Run 497 stress 0.05602052 
    ## Run 498 stress 0.05466704 
    ## Run 499 stress 0.05829929 
    ## Run 500 stress 0.05466708 
    ## *** Best solution repeated 109 times

``` r
# Ocean sites and mixed lakes
PD_beta_OM_NMDS <- metaMDS(PD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1016168 
    ## Run 1 stress 0.1016167 
    ## ... New best solution
    ## ... Procrustes: rmse 8.002264e-05  max resid 0.0002212783 
    ## ... Similar to previous best
    ## Run 2 stress 0.1356147 
    ## Run 3 stress 0.1016166 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001373519  max resid 0.0003829742 
    ## ... Similar to previous best
    ## Run 4 stress 0.1016169 
    ## ... Procrustes: rmse 0.0002328239  max resid 0.000648225 
    ## ... Similar to previous best
    ## Run 5 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004067302  max resid 0.00113933 
    ## ... Similar to previous best
    ## Run 6 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007825569  max resid 0.00219037 
    ## ... Similar to previous best
    ## Run 7 stress 0.1016165 
    ## ... Procrustes: rmse 0.0003528457  max resid 0.0009914708 
    ## ... Similar to previous best
    ## Run 8 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 6.224864e-05  max resid 0.0001714221 
    ## ... Similar to previous best
    ## Run 9 stress 0.1016164 
    ## ... Procrustes: rmse 1.227977e-05  max resid 3.148635e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.1356146 
    ## Run 11 stress 0.1016164 
    ## ... Procrustes: rmse 5.343278e-05  max resid 0.0001506705 
    ## ... Similar to previous best
    ## Run 12 stress 0.1016164 
    ## ... Procrustes: rmse 3.398209e-05  max resid 8.663372e-05 
    ## ... Similar to previous best
    ## Run 13 stress 0.1016164 
    ## ... Procrustes: rmse 5.266953e-06  max resid 1.500052e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001233934  max resid 0.0003485535 
    ## ... Similar to previous best
    ## Run 15 stress 0.1016164 
    ## ... Procrustes: rmse 2.426077e-05  max resid 6.82715e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003228124  max resid 0.0009039225 
    ## ... Similar to previous best
    ## Run 17 stress 0.1341868 
    ## Run 18 stress 0.1530428 
    ## Run 19 stress 0.1016164 
    ## ... Procrustes: rmse 8.413137e-05  max resid 0.000236913 
    ## ... Similar to previous best
    ## Run 20 stress 0.1016164 
    ## ... Procrustes: rmse 2.036625e-05  max resid 5.745925e-05 
    ## ... Similar to previous best
    ## Run 21 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005804518  max resid 0.001623072 
    ## ... Similar to previous best
    ## Run 22 stress 0.1530424 
    ## Run 23 stress 0.1016164 
    ## ... Procrustes: rmse 1.476139e-05  max resid 3.739545e-05 
    ## ... Similar to previous best
    ## Run 24 stress 0.1356146 
    ## Run 25 stress 0.1016164 
    ## ... Procrustes: rmse 6.804151e-06  max resid 1.394305e-05 
    ## ... Similar to previous best
    ## Run 26 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 3.824909e-06  max resid 9.73723e-06 
    ## ... Similar to previous best
    ## Run 27 stress 0.3411821 
    ## Run 28 stress 0.1016164 
    ## ... Procrustes: rmse 7.50083e-05  max resid 0.000211587 
    ## ... Similar to previous best
    ## Run 29 stress 0.1016164 
    ## ... Procrustes: rmse 1.195873e-05  max resid 3.321119e-05 
    ## ... Similar to previous best
    ## Run 30 stress 0.1519033 
    ## Run 31 stress 0.1016164 
    ## ... Procrustes: rmse 3.60649e-05  max resid 0.0001014574 
    ## ... Similar to previous best
    ## Run 32 stress 0.1530425 
    ## Run 33 stress 0.1341868 
    ## Run 34 stress 0.1016164 
    ## ... Procrustes: rmse 2.356524e-05  max resid 6.637366e-05 
    ## ... Similar to previous best
    ## Run 35 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005283115  max resid 0.001477199 
    ## ... Similar to previous best
    ## Run 36 stress 0.1341868 
    ## Run 37 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005562788  max resid 0.001555367 
    ## ... Similar to previous best
    ## Run 38 stress 0.1341868 
    ## Run 39 stress 0.1016164 
    ## ... Procrustes: rmse 4.447316e-05  max resid 0.0001050121 
    ## ... Similar to previous best
    ## Run 40 stress 0.1016164 
    ## ... Procrustes: rmse 7.603119e-05  max resid 0.000214105 
    ## ... Similar to previous best
    ## Run 41 stress 0.1341868 
    ## Run 42 stress 0.1356149 
    ## Run 43 stress 0.1016164 
    ## ... Procrustes: rmse 7.023035e-05  max resid 0.0001974725 
    ## ... Similar to previous best
    ## Run 44 stress 0.1356148 
    ## Run 45 stress 0.1341868 
    ## Run 46 stress 0.1016164 
    ## ... Procrustes: rmse 6.534511e-05  max resid 0.000183776 
    ## ... Similar to previous best
    ## Run 47 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004152508  max resid 0.001160685 
    ## ... Similar to previous best
    ## Run 48 stress 0.1356148 
    ## Run 49 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 5.012187e-06  max resid 1.291096e-05 
    ## ... Similar to previous best
    ## Run 50 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004504676  max resid 0.00125983 
    ## ... Similar to previous best
    ## Run 51 stress 0.1341868 
    ## Run 52 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003908148  max resid 0.001093569 
    ## ... Similar to previous best
    ## Run 53 stress 0.1016164 
    ## ... Procrustes: rmse 1.111666e-05  max resid 3.161061e-05 
    ## ... Similar to previous best
    ## Run 54 stress 0.1016164 
    ## ... Procrustes: rmse 9.554131e-06  max resid 2.467243e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.1016164 
    ## ... Procrustes: rmse 3.356903e-06  max resid 9.073096e-06 
    ## ... Similar to previous best
    ## Run 56 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001473958  max resid 0.000414445 
    ## ... Similar to previous best
    ## Run 57 stress 0.1016164 
    ## ... Procrustes: rmse 1.686918e-05  max resid 4.642001e-05 
    ## ... Similar to previous best
    ## Run 58 stress 0.101617 
    ## ... Procrustes: rmse 0.0006527442  max resid 0.001824535 
    ## ... Similar to previous best
    ## Run 59 stress 0.1341868 
    ## Run 60 stress 0.1016164 
    ## ... Procrustes: rmse 8.693242e-05  max resid 0.0002435807 
    ## ... Similar to previous best
    ## Run 61 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004578822  max resid 0.001282517 
    ## ... Similar to previous best
    ## Run 62 stress 0.1341868 
    ## Run 63 stress 0.1016164 
    ## ... Procrustes: rmse 9.994142e-05  max resid 0.000281123 
    ## ... Similar to previous best
    ## Run 64 stress 0.1356145 
    ## Run 65 stress 0.1016164 
    ## ... Procrustes: rmse 7.575872e-05  max resid 0.000213625 
    ## ... Similar to previous best
    ## Run 66 stress 0.1016164 
    ## ... Procrustes: rmse 5.195346e-06  max resid 1.431293e-05 
    ## ... Similar to previous best
    ## Run 67 stress 0.1016168 
    ## ... Procrustes: rmse 0.000534156  max resid 0.001493832 
    ## ... Similar to previous best
    ## Run 68 stress 0.1016164 
    ## ... Procrustes: rmse 6.29751e-05  max resid 0.0001774002 
    ## ... Similar to previous best
    ## Run 69 stress 0.1341868 
    ## Run 70 stress 0.3195522 
    ## Run 71 stress 0.1341868 
    ## Run 72 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006202729  max resid 0.001733522 
    ## ... Similar to previous best
    ## Run 73 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003261946  max resid 0.0009122613 
    ## ... Similar to previous best
    ## Run 74 stress 0.135615 
    ## Run 75 stress 0.1016164 
    ## ... Procrustes: rmse 1.827511e-05  max resid 5.084571e-05 
    ## ... Similar to previous best
    ## Run 76 stress 0.1016164 
    ## ... Procrustes: rmse 5.332013e-05  max resid 0.0001488942 
    ## ... Similar to previous best
    ## Run 77 stress 0.1016164 
    ## ... Procrustes: rmse 1.892394e-05  max resid 4.887541e-05 
    ## ... Similar to previous best
    ## Run 78 stress 0.1016164 
    ## ... Procrustes: rmse 8.125795e-05  max resid 0.0002289793 
    ## ... Similar to previous best
    ## Run 79 stress 0.1356147 
    ## Run 80 stress 0.1016164 
    ## ... Procrustes: rmse 2.19512e-05  max resid 6.1443e-05 
    ## ... Similar to previous best
    ## Run 81 stress 0.1016164 
    ## ... Procrustes: rmse 5.782336e-05  max resid 0.0001626743 
    ## ... Similar to previous best
    ## Run 82 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002003781  max resid 0.0005619839 
    ## ... Similar to previous best
    ## Run 83 stress 0.1016164 
    ## ... Procrustes: rmse 2.939265e-05  max resid 8.246882e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.101617 
    ## ... Procrustes: rmse 0.0006532907  max resid 0.001824412 
    ## ... Similar to previous best
    ## Run 85 stress 0.1016164 
    ## ... Procrustes: rmse 2.63651e-05  max resid 7.397258e-05 
    ## ... Similar to previous best
    ## Run 86 stress 0.1016169 
    ## ... Procrustes: rmse 0.0004915233  max resid 0.001366964 
    ## ... Similar to previous best
    ## Run 87 stress 0.1016164 
    ## ... Procrustes: rmse 9.234498e-06  max resid 2.618233e-05 
    ## ... Similar to previous best
    ## Run 88 stress 0.1341868 
    ## Run 89 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003465344  max resid 0.0009708407 
    ## ... Similar to previous best
    ## Run 90 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001049819  max resid 0.0002926527 
    ## ... Similar to previous best
    ## Run 91 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005197075  max resid 0.001451988 
    ## ... Similar to previous best
    ## Run 92 stress 0.1016164 
    ## ... Procrustes: rmse 5.961284e-06  max resid 1.218219e-05 
    ## ... Similar to previous best
    ## Run 93 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004159509  max resid 0.001163302 
    ## ... Similar to previous best
    ## Run 94 stress 0.101617 
    ## ... Procrustes: rmse 0.0007051406  max resid 0.001970726 
    ## ... Similar to previous best
    ## Run 95 stress 0.1356146 
    ## Run 96 stress 0.1341868 
    ## Run 97 stress 0.1356148 
    ## Run 98 stress 0.1016169 
    ## ... Procrustes: rmse 0.000581547  max resid 0.00162575 
    ## ... Similar to previous best
    ## Run 99 stress 0.1519033 
    ## Run 100 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 4.797906e-06  max resid 1.37651e-05 
    ## ... Similar to previous best
    ## Run 101 stress 0.1356148 
    ## Run 102 stress 0.1016164 
    ## ... Procrustes: rmse 3.701394e-05  max resid 0.0001041105 
    ## ... Similar to previous best
    ## Run 103 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002892595  max resid 0.0008114652 
    ## ... Similar to previous best
    ## Run 104 stress 0.1016165 
    ## ... Procrustes: rmse 0.000239962  max resid 0.0006724163 
    ## ... Similar to previous best
    ## Run 105 stress 0.1356145 
    ## Run 106 stress 0.1530423 
    ## Run 107 stress 0.1356146 
    ## Run 108 stress 0.1016164 
    ## ... Procrustes: rmse 6.461671e-05  max resid 0.0001782536 
    ## ... Similar to previous best
    ## Run 109 stress 0.1356144 
    ## Run 110 stress 0.1016164 
    ## ... Procrustes: rmse 1.790671e-05  max resid 4.052314e-05 
    ## ... Similar to previous best
    ## Run 111 stress 0.101617 
    ## ... Procrustes: rmse 0.0006653792  max resid 0.001858296 
    ## ... Similar to previous best
    ## Run 112 stress 0.3495546 
    ## Run 113 stress 0.1016164 
    ## ... Procrustes: rmse 8.69237e-05  max resid 0.0002449668 
    ## ... Similar to previous best
    ## Run 114 stress 0.1016164 
    ## ... Procrustes: rmse 4.348424e-05  max resid 0.000118494 
    ## ... Similar to previous best
    ## Run 115 stress 0.1016164 
    ## ... Procrustes: rmse 1.092683e-05  max resid 2.890674e-05 
    ## ... Similar to previous best
    ## Run 116 stress 0.1341868 
    ## Run 117 stress 0.1356144 
    ## Run 118 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004090189  max resid 0.001141197 
    ## ... Similar to previous best
    ## Run 119 stress 0.1016164 
    ## ... Procrustes: rmse 2.993894e-05  max resid 8.438994e-05 
    ## ... Similar to previous best
    ## Run 120 stress 0.1016164 
    ## ... Procrustes: rmse 2.045869e-05  max resid 5.743271e-05 
    ## ... Similar to previous best
    ## Run 121 stress 0.1016164 
    ## ... Procrustes: rmse 3.119047e-05  max resid 8.777916e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.1356147 
    ## Run 123 stress 0.1016166 
    ## ... Procrustes: rmse 0.000329487  max resid 0.0009218625 
    ## ... Similar to previous best
    ## Run 124 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004283104  max resid 0.001198122 
    ## ... Similar to previous best
    ## Run 125 stress 0.1341868 
    ## Run 126 stress 0.1016164 
    ## ... Procrustes: rmse 2.057115e-05  max resid 5.804068e-05 
    ## ... Similar to previous best
    ## Run 127 stress 0.1356148 
    ## Run 128 stress 0.1016164 
    ## ... Procrustes: rmse 6.340104e-05  max resid 0.0001775468 
    ## ... Similar to previous best
    ## Run 129 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003513468  max resid 0.0009833093 
    ## ... Similar to previous best
    ## Run 130 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003940641  max resid 0.001103795 
    ## ... Similar to previous best
    ## Run 131 stress 0.135615 
    ## Run 132 stress 0.1016164 
    ## ... Procrustes: rmse 1.970471e-05  max resid 5.49627e-05 
    ## ... Similar to previous best
    ## Run 133 stress 0.135615 
    ## Run 134 stress 0.1356146 
    ## Run 135 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001875339  max resid 0.0005263432 
    ## ... Similar to previous best
    ## Run 136 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005033994  max resid 0.001407402 
    ## ... Similar to previous best
    ## Run 137 stress 0.1530425 
    ## Run 138 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001056949  max resid 0.0002977555 
    ## ... Similar to previous best
    ## Run 139 stress 0.1016164 
    ## ... Procrustes: rmse 1.053919e-05  max resid 2.676864e-05 
    ## ... Similar to previous best
    ## Run 140 stress 0.1016164 
    ## ... Procrustes: rmse 1.393621e-05  max resid 3.743802e-05 
    ## ... Similar to previous best
    ## Run 141 stress 0.1016164 
    ## ... Procrustes: rmse 7.734251e-06  max resid 2.149047e-05 
    ## ... Similar to previous best
    ## Run 142 stress 0.1016164 
    ## ... Procrustes: rmse 1.453311e-05  max resid 4.084377e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002911003  max resid 0.000815292 
    ## ... Similar to previous best
    ## Run 144 stress 0.1341868 
    ## Run 145 stress 0.1519033 
    ## Run 146 stress 0.1356148 
    ## Run 147 stress 0.1530427 
    ## Run 148 stress 0.3495709 
    ## Run 149 stress 0.1016166 
    ## ... Procrustes: rmse 0.000373435  max resid 0.001045285 
    ## ... Similar to previous best
    ## Run 150 stress 0.1356147 
    ## Run 151 stress 0.1016164 
    ## ... Procrustes: rmse 1.151392e-05  max resid 3.242346e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.1016164 
    ## ... Procrustes: rmse 2.493581e-05  max resid 6.007183e-05 
    ## ... Similar to previous best
    ## Run 153 stress 0.1016164 
    ## ... Procrustes: rmse 1.899533e-05  max resid 5.350426e-05 
    ## ... Similar to previous best
    ## Run 154 stress 0.1341868 
    ## Run 155 stress 0.1016164 
    ## ... Procrustes: rmse 6.092089e-05  max resid 0.0001713924 
    ## ... Similar to previous best
    ## Run 156 stress 0.317577 
    ## Run 157 stress 0.1356147 
    ## Run 158 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002637379  max resid 0.0007390069 
    ## ... Similar to previous best
    ## Run 159 stress 0.1016164 
    ## ... Procrustes: rmse 8.839203e-05  max resid 0.0002491724 
    ## ... Similar to previous best
    ## Run 160 stress 0.1341868 
    ## Run 161 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006027673  max resid 0.001683856 
    ## ... Similar to previous best
    ## Run 162 stress 0.1016164 
    ## ... Procrustes: rmse 8.390865e-06  max resid 2.014879e-05 
    ## ... Similar to previous best
    ## Run 163 stress 0.1530425 
    ## Run 164 stress 0.1016164 
    ## ... Procrustes: rmse 9.971142e-05  max resid 0.0002796946 
    ## ... Similar to previous best
    ## Run 165 stress 0.1341868 
    ## Run 166 stress 0.1016164 
    ## ... Procrustes: rmse 2.670064e-05  max resid 7.538047e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003124831  max resid 0.000874754 
    ## ... Similar to previous best
    ## Run 168 stress 0.1016164 
    ## ... Procrustes: rmse 1.552083e-05  max resid 4.356146e-05 
    ## ... Similar to previous best
    ## Run 169 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006340837  max resid 0.001772486 
    ## ... Similar to previous best
    ## Run 170 stress 0.1016164 
    ## ... Procrustes: rmse 3.194779e-05  max resid 9.011741e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.1016164 
    ## ... Procrustes: rmse 8.309869e-06  max resid 2.344057e-05 
    ## ... Similar to previous best
    ## Run 172 stress 0.1341868 
    ## Run 173 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001990034  max resid 0.0005547735 
    ## ... Similar to previous best
    ## Run 174 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002809773  max resid 0.0007802181 
    ## ... Similar to previous best
    ## Run 175 stress 0.1016164 
    ## ... Procrustes: rmse 4.113879e-05  max resid 0.0001138736 
    ## ... Similar to previous best
    ## Run 176 stress 0.1016164 
    ## ... Procrustes: rmse 1.475489e-06  max resid 3.1305e-06 
    ## ... Similar to previous best
    ## Run 177 stress 0.1016164 
    ## ... Procrustes: rmse 2.764277e-05  max resid 6.267155e-05 
    ## ... Similar to previous best
    ## Run 178 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003877297  max resid 0.001085109 
    ## ... Similar to previous best
    ## Run 179 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005858714  max resid 0.001635234 
    ## ... Similar to previous best
    ## Run 180 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001017248  max resid 0.0002868544 
    ## ... Similar to previous best
    ## Run 181 stress 0.1016164 
    ## ... Procrustes: rmse 7.625892e-06  max resid 2.145661e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.1016164 
    ## ... Procrustes: rmse 2.741145e-06  max resid 6.492408e-06 
    ## ... Similar to previous best
    ## Run 183 stress 0.1356143 
    ## Run 184 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005753136  max resid 0.001608929 
    ## ... Similar to previous best
    ## Run 185 stress 0.1356146 
    ## Run 186 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004405283  max resid 0.001230695 
    ## ... Similar to previous best
    ## Run 187 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001323789  max resid 0.0003672115 
    ## ... Similar to previous best
    ## Run 188 stress 0.1016164 
    ## ... Procrustes: rmse 5.037396e-05  max resid 0.0001419743 
    ## ... Similar to previous best
    ## Run 189 stress 0.1016164 
    ## ... Procrustes: rmse 2.473359e-05  max resid 6.877857e-05 
    ## ... Similar to previous best
    ## Run 190 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001742498  max resid 0.0004892997 
    ## ... Similar to previous best
    ## Run 191 stress 0.1016164 
    ## ... Procrustes: rmse 2.209665e-05  max resid 3.948438e-05 
    ## ... Similar to previous best
    ## Run 192 stress 0.1016164 
    ## ... Procrustes: rmse 3.740256e-05  max resid 0.0001022299 
    ## ... Similar to previous best
    ## Run 193 stress 0.1016164 
    ## ... Procrustes: rmse 1.076066e-05  max resid 3.04896e-05 
    ## ... Similar to previous best
    ## Run 194 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002135988  max resid 0.000598383 
    ## ... Similar to previous best
    ## Run 195 stress 0.1016164 
    ## ... Procrustes: rmse 2.376661e-05  max resid 6.369123e-05 
    ## ... Similar to previous best
    ## Run 196 stress 0.1356145 
    ## Run 197 stress 0.1016164 
    ## ... Procrustes: rmse 1.1157e-05  max resid 3.118586e-05 
    ## ... Similar to previous best
    ## Run 198 stress 0.1341868 
    ## Run 199 stress 0.1016164 
    ## ... Procrustes: rmse 1.070583e-05  max resid 3.006761e-05 
    ## ... Similar to previous best
    ## Run 200 stress 0.3398622 
    ## Run 201 stress 0.1341868 
    ## Run 202 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005436616  max resid 0.001519987 
    ## ... Similar to previous best
    ## Run 203 stress 0.1016164 
    ## ... Procrustes: rmse 9.363824e-05  max resid 0.0002638825 
    ## ... Similar to previous best
    ## Run 204 stress 0.1016171 
    ## ... Procrustes: rmse 0.0006945367  max resid 0.001935073 
    ## ... Similar to previous best
    ## Run 205 stress 0.1341868 
    ## Run 206 stress 0.1356145 
    ## Run 207 stress 0.101617 
    ## ... Procrustes: rmse 0.0006660267  max resid 0.001861871 
    ## ... Similar to previous best
    ## Run 208 stress 0.1016164 
    ## ... Procrustes: rmse 1.380425e-05  max resid 3.895372e-05 
    ## ... Similar to previous best
    ## Run 209 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002081447  max resid 0.0005831899 
    ## ... Similar to previous best
    ## Run 210 stress 0.1356147 
    ## Run 211 stress 0.1016164 
    ## ... Procrustes: rmse 5.221312e-05  max resid 0.0001474142 
    ## ... Similar to previous best
    ## Run 212 stress 0.1016164 
    ## ... Procrustes: rmse 1.745774e-05  max resid 4.875472e-05 
    ## ... Similar to previous best
    ## Run 213 stress 0.1016164 
    ## ... Procrustes: rmse 1.631803e-05  max resid 4.318026e-05 
    ## ... Similar to previous best
    ## Run 214 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003732217  max resid 0.001042164 
    ## ... Similar to previous best
    ## Run 215 stress 0.1356148 
    ## Run 216 stress 0.2074952 
    ## Run 217 stress 0.1341868 
    ## Run 218 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005670969  max resid 0.001583379 
    ## ... Similar to previous best
    ## Run 219 stress 0.1530423 
    ## Run 220 stress 0.1016164 
    ## ... Procrustes: rmse 2.14404e-05  max resid 6.011951e-05 
    ## ... Similar to previous best
    ## Run 221 stress 0.1016164 
    ## ... Procrustes: rmse 4.106357e-05  max resid 0.0001155766 
    ## ... Similar to previous best
    ## Run 222 stress 0.3495732 
    ## Run 223 stress 0.1016164 
    ## ... Procrustes: rmse 5.347058e-05  max resid 0.0001510628 
    ## ... Similar to previous best
    ## Run 224 stress 0.1341868 
    ## Run 225 stress 0.1356147 
    ## Run 226 stress 0.1356148 
    ## Run 227 stress 0.1016164 
    ## ... Procrustes: rmse 7.165454e-05  max resid 0.0002010582 
    ## ... Similar to previous best
    ## Run 228 stress 0.1341868 
    ## Run 229 stress 0.1016164 
    ## ... Procrustes: rmse 1.001953e-05  max resid 2.280611e-05 
    ## ... Similar to previous best
    ## Run 230 stress 0.1016164 
    ## ... Procrustes: rmse 7.26519e-05  max resid 0.0002045076 
    ## ... Similar to previous best
    ## Run 231 stress 0.1016166 
    ## ... Procrustes: rmse 0.000330914  max resid 0.0009259001 
    ## ... Similar to previous best
    ## Run 232 stress 0.1341868 
    ## Run 233 stress 0.1356149 
    ## Run 234 stress 0.1356145 
    ## Run 235 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004652596  max resid 0.001302537 
    ## ... Similar to previous best
    ## Run 236 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004870954  max resid 0.001358366 
    ## ... Similar to previous best
    ## Run 237 stress 0.1341868 
    ## Run 238 stress 0.1530422 
    ## Run 239 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005243291  max resid 0.001466466 
    ## ... Similar to previous best
    ## Run 240 stress 0.1016164 
    ## ... Procrustes: rmse 2.596665e-05  max resid 5.972196e-05 
    ## ... Similar to previous best
    ## Run 241 stress 0.1341868 
    ## Run 242 stress 0.1016164 
    ## ... Procrustes: rmse 1.073963e-05  max resid 3.028714e-05 
    ## ... Similar to previous best
    ## Run 243 stress 0.1341868 
    ## Run 244 stress 0.1341868 
    ## Run 245 stress 0.1016164 
    ## ... Procrustes: rmse 9.358389e-06  max resid 2.634387e-05 
    ## ... Similar to previous best
    ## Run 246 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004512367  max resid 0.001260996 
    ## ... Similar to previous best
    ## Run 247 stress 0.1341868 
    ## Run 248 stress 0.1016164 
    ## ... Procrustes: rmse 9.094564e-06  max resid 2.480842e-05 
    ## ... Similar to previous best
    ## Run 249 stress 0.1016164 
    ## ... Procrustes: rmse 7.496533e-05  max resid 0.0002114743 
    ## ... Similar to previous best
    ## Run 250 stress 0.1016164 
    ## ... Procrustes: rmse 1.438086e-05  max resid 4.041435e-05 
    ## ... Similar to previous best
    ## Run 251 stress 0.1341868 
    ## Run 252 stress 0.1016164 
    ## ... Procrustes: rmse 2.199311e-05  max resid 6.158769e-05 
    ## ... Similar to previous best
    ## Run 253 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005137492  max resid 0.001436759 
    ## ... Similar to previous best
    ## Run 254 stress 0.1356145 
    ## Run 255 stress 0.1356151 
    ## Run 256 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002894864  max resid 0.0008116536 
    ## ... Similar to previous best
    ## Run 257 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005022493  max resid 0.001405587 
    ## ... Similar to previous best
    ## Run 258 stress 0.1341868 
    ## Run 259 stress 0.101617 
    ## ... Procrustes: rmse 0.0006436571  max resid 0.001797383 
    ## ... Similar to previous best
    ## Run 260 stress 0.1341868 
    ## Run 261 stress 0.1530427 
    ## Run 262 stress 0.1016164 
    ## ... Procrustes: rmse 9.04259e-06  max resid 2.544266e-05 
    ## ... Similar to previous best
    ## Run 263 stress 0.1016164 
    ## ... Procrustes: rmse 8.013473e-05  max resid 0.0002253826 
    ## ... Similar to previous best
    ## Run 264 stress 0.1016164 
    ## ... Procrustes: rmse 8.621692e-06  max resid 2.412067e-05 
    ## ... Similar to previous best
    ## Run 265 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003661439  max resid 0.001024787 
    ## ... Similar to previous best
    ## Run 266 stress 0.1016164 
    ## ... Procrustes: rmse 2.337083e-05  max resid 6.591595e-05 
    ## ... Similar to previous best
    ## Run 267 stress 0.1016164 
    ## ... Procrustes: rmse 9.323603e-05  max resid 0.0002630689 
    ## ... Similar to previous best
    ## Run 268 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005687431  max resid 0.001585984 
    ## ... Similar to previous best
    ## Run 269 stress 0.1356144 
    ## Run 270 stress 0.1016164 
    ## ... Procrustes: rmse 9.493968e-06  max resid 2.654637e-05 
    ## ... Similar to previous best
    ## Run 271 stress 0.1016164 
    ## ... Procrustes: rmse 5.000361e-05  max resid 0.0001395804 
    ## ... Similar to previous best
    ## Run 272 stress 0.1016166 
    ## ... Procrustes: rmse 0.000354868  max resid 0.000994813 
    ## ... Similar to previous best
    ## Run 273 stress 0.1016164 
    ## ... Procrustes: rmse 4.323628e-05  max resid 0.0001155279 
    ## ... Similar to previous best
    ## Run 274 stress 0.1341868 
    ## Run 275 stress 0.1016164 
    ## ... Procrustes: rmse 4.358492e-05  max resid 0.0001230168 
    ## ... Similar to previous best
    ## Run 276 stress 0.1016164 
    ## ... Procrustes: rmse 6.232451e-05  max resid 0.0001755326 
    ## ... Similar to previous best
    ## Run 277 stress 0.1016164 
    ## ... Procrustes: rmse 1.705061e-05  max resid 4.806467e-05 
    ## ... Similar to previous best
    ## Run 278 stress 0.1519033 
    ## Run 279 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005099199  max resid 0.001423598 
    ## ... Similar to previous best
    ## Run 280 stress 0.1016164 
    ## ... Procrustes: rmse 9.43499e-05  max resid 0.0002658413 
    ## ... Similar to previous best
    ## Run 281 stress 0.1341868 
    ## Run 282 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002224831  max resid 0.0006234043 
    ## ... Similar to previous best
    ## Run 283 stress 0.1016165 
    ## ... Procrustes: rmse 0.000215152  max resid 0.0006044614 
    ## ... Similar to previous best
    ## Run 284 stress 0.1341868 
    ## Run 285 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001144895  max resid 0.0003216171 
    ## ... Similar to previous best
    ## Run 286 stress 0.1016164 
    ## ... Procrustes: rmse 7.080768e-05  max resid 0.0001977401 
    ## ... Similar to previous best
    ## Run 287 stress 0.1530423 
    ## Run 288 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002449615  max resid 0.0006860996 
    ## ... Similar to previous best
    ## Run 289 stress 0.1356147 
    ## Run 290 stress 0.1016164 
    ## ... Procrustes: rmse 7.16039e-05  max resid 0.0002019062 
    ## ... Similar to previous best
    ## Run 291 stress 0.1016164 
    ## ... Procrustes: rmse 9.450512e-05  max resid 0.0002659124 
    ## ... Similar to previous best
    ## Run 292 stress 0.1016164 
    ## ... Procrustes: rmse 6.410188e-05  max resid 0.000177598 
    ## ... Similar to previous best
    ## Run 293 stress 0.1016164 
    ## ... Procrustes: rmse 3.038302e-05  max resid 8.546038e-05 
    ## ... Similar to previous best
    ## Run 294 stress 0.1016164 
    ## ... Procrustes: rmse 7.870786e-06  max resid 2.215954e-05 
    ## ... Similar to previous best
    ## Run 295 stress 0.1016164 
    ## ... Procrustes: rmse 4.43109e-05  max resid 0.0001250463 
    ## ... Similar to previous best
    ## Run 296 stress 0.1530422 
    ## Run 297 stress 0.3027329 
    ## Run 298 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001028246  max resid 0.0002897522 
    ## ... Similar to previous best
    ## Run 299 stress 0.1016164 
    ## ... Procrustes: rmse 5.970512e-06  max resid 1.671603e-05 
    ## ... Similar to previous best
    ## Run 300 stress 0.1016164 
    ## ... Procrustes: rmse 7.475421e-05  max resid 0.0002086576 
    ## ... Similar to previous best
    ## Run 301 stress 0.1530428 
    ## Run 302 stress 0.1016164 
    ## ... Procrustes: rmse 3.536799e-06  max resid 9.819002e-06 
    ## ... Similar to previous best
    ## Run 303 stress 0.1356149 
    ## Run 304 stress 0.101617 
    ## ... Procrustes: rmse 0.0006931264  max resid 0.001936353 
    ## ... Similar to previous best
    ## Run 305 stress 0.1016164 
    ## ... Procrustes: rmse 2.783142e-05  max resid 7.843586e-05 
    ## ... Similar to previous best
    ## Run 306 stress 0.1016164 
    ## ... Procrustes: rmse 1.582791e-05  max resid 4.389414e-05 
    ## ... Similar to previous best
    ## Run 307 stress 0.1016164 
    ## ... Procrustes: rmse 3.330645e-05  max resid 9.388125e-05 
    ## ... Similar to previous best
    ## Run 308 stress 0.1016164 
    ## ... Procrustes: rmse 7.051983e-05  max resid 0.0001980827 
    ## ... Similar to previous best
    ## Run 309 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004373434  max resid 0.001220493 
    ## ... Similar to previous best
    ## Run 310 stress 0.1016164 
    ## ... Procrustes: rmse 6.254739e-05  max resid 0.0001762111 
    ## ... Similar to previous best
    ## Run 311 stress 0.299814 
    ## Run 312 stress 0.1341868 
    ## Run 313 stress 0.1016164 
    ## ... Procrustes: rmse 1.254264e-05  max resid 3.504642e-05 
    ## ... Similar to previous best
    ## Run 314 stress 0.101617 
    ## ... Procrustes: rmse 0.0006424419  max resid 0.001794468 
    ## ... Similar to previous best
    ## Run 315 stress 0.1341868 
    ## Run 316 stress 0.1356147 
    ## Run 317 stress 0.1356148 
    ## Run 318 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004293172  max resid 0.001200659 
    ## ... Similar to previous best
    ## Run 319 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004795128  max resid 0.0013422 
    ## ... Similar to previous best
    ## Run 320 stress 0.1341868 
    ## Run 321 stress 0.1356148 
    ## Run 322 stress 0.1341868 
    ## Run 323 stress 0.1016164 
    ## ... Procrustes: rmse 4.627255e-05  max resid 0.00012029 
    ## ... Similar to previous best
    ## Run 324 stress 0.1356147 
    ## Run 325 stress 0.1016171 
    ## ... Procrustes: rmse 0.0005778758  max resid 0.001604293 
    ## ... Similar to previous best
    ## Run 326 stress 0.1341868 
    ## Run 327 stress 0.1016164 
    ## ... Procrustes: rmse 6.506525e-05  max resid 0.0001829059 
    ## ... Similar to previous best
    ## Run 328 stress 0.1016164 
    ## ... Procrustes: rmse 8.300508e-05  max resid 0.0002317926 
    ## ... Similar to previous best
    ## Run 329 stress 0.1356147 
    ## Run 330 stress 0.1341868 
    ## Run 331 stress 0.1016164 
    ## ... Procrustes: rmse 1.771372e-05  max resid 4.944113e-05 
    ## ... Similar to previous best
    ## Run 332 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004399087  max resid 0.001229302 
    ## ... Similar to previous best
    ## Run 333 stress 0.1356148 
    ## Run 334 stress 0.1016164 
    ## ... Procrustes: rmse 1.53366e-05  max resid 2.309466e-05 
    ## ... Similar to previous best
    ## Run 335 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002993957  max resid 0.0008381891 
    ## ... Similar to previous best
    ## Run 336 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003731379  max resid 0.001044383 
    ## ... Similar to previous best
    ## Run 337 stress 0.1341868 
    ## Run 338 stress 0.1016164 
    ## ... Procrustes: rmse 8.236756e-07  max resid 1.759746e-06 
    ## ... Similar to previous best
    ## Run 339 stress 0.1016169 
    ## ... Procrustes: rmse 0.000568547  max resid 0.00158825 
    ## ... Similar to previous best
    ## Run 340 stress 0.1356147 
    ## Run 341 stress 0.1341868 
    ## Run 342 stress 0.1356151 
    ## Run 343 stress 0.1356149 
    ## Run 344 stress 0.1530427 
    ## Run 345 stress 0.1356147 
    ## Run 346 stress 0.1356148 
    ## Run 347 stress 0.1341868 
    ## Run 348 stress 0.1356147 
    ## Run 349 stress 0.1016164 
    ## ... Procrustes: rmse 2.151721e-05  max resid 6.066301e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003894783  max resid 0.001089757 
    ## ... Similar to previous best
    ## Run 351 stress 0.1016164 
    ## ... Procrustes: rmse 3.560989e-05  max resid 0.0001005968 
    ## ... Similar to previous best
    ## Run 352 stress 0.1016164 
    ## ... Procrustes: rmse 1.783748e-05  max resid 5.030167e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.1016164 
    ## ... Procrustes: rmse 8.637818e-05  max resid 0.0002422555 
    ## ... Similar to previous best
    ## Run 354 stress 0.1016164 
    ## ... Procrustes: rmse 4.047154e-05  max resid 0.0001098251 
    ## ... Similar to previous best
    ## Run 355 stress 0.1016164 
    ## ... Procrustes: rmse 5.438545e-05  max resid 0.0001535338 
    ## ... Similar to previous best
    ## Run 356 stress 0.1016164 
    ## ... Procrustes: rmse 3.068915e-05  max resid 8.585371e-05 
    ## ... Similar to previous best
    ## Run 357 stress 0.1530424 
    ## Run 358 stress 0.1341868 
    ## Run 359 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003473267  max resid 0.0009716524 
    ## ... Similar to previous best
    ## Run 360 stress 0.1341868 
    ## Run 361 stress 0.1016164 
    ## ... Procrustes: rmse 6.882331e-06  max resid 1.791875e-05 
    ## ... Similar to previous best
    ## Run 362 stress 0.1016164 
    ## ... Procrustes: rmse 0.000157198  max resid 0.0004415295 
    ## ... Similar to previous best
    ## Run 363 stress 0.1356144 
    ## Run 364 stress 0.1016164 
    ## ... Procrustes: rmse 3.424979e-05  max resid 9.638145e-05 
    ## ... Similar to previous best
    ## Run 365 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005430989  max resid 0.001519802 
    ## ... Similar to previous best
    ## Run 366 stress 0.1016164 
    ## ... Procrustes: rmse 1.553469e-05  max resid 4.236428e-05 
    ## ... Similar to previous best
    ## Run 367 stress 0.1016164 
    ## ... Procrustes: rmse 2.788363e-05  max resid 7.828027e-05 
    ## ... Similar to previous best
    ## Run 368 stress 0.1356148 
    ## Run 369 stress 0.1356148 
    ## Run 370 stress 0.1016164 
    ## ... Procrustes: rmse 3.451125e-05  max resid 8.544753e-05 
    ## ... Similar to previous best
    ## Run 371 stress 0.1341868 
    ## Run 372 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001993888  max resid 0.0005581806 
    ## ... Similar to previous best
    ## Run 373 stress 0.1341868 
    ## Run 374 stress 0.1356147 
    ## Run 375 stress 0.1016164 
    ## ... Procrustes: rmse 3.804559e-05  max resid 0.0001043943 
    ## ... Similar to previous best
    ## Run 376 stress 0.1341868 
    ## Run 377 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003451785  max resid 0.0009662248 
    ## ... Similar to previous best
    ## Run 378 stress 0.1016164 
    ## ... Procrustes: rmse 3.637441e-05  max resid 0.0001020575 
    ## ... Similar to previous best
    ## Run 379 stress 0.1016164 
    ## ... Procrustes: rmse 3.85868e-05  max resid 0.0001086521 
    ## ... Similar to previous best
    ## Run 380 stress 0.1356151 
    ## Run 381 stress 0.1016164 
    ## ... Procrustes: rmse 6.586596e-06  max resid 1.729013e-05 
    ## ... Similar to previous best
    ## Run 382 stress 0.1356144 
    ## Run 383 stress 0.135615 
    ## Run 384 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004080696  max resid 0.001141313 
    ## ... Similar to previous best
    ## Run 385 stress 0.1016164 
    ## ... Procrustes: rmse 1.695917e-05  max resid 4.697919e-05 
    ## ... Similar to previous best
    ## Run 386 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001719139  max resid 0.000482385 
    ## ... Similar to previous best
    ## Run 387 stress 0.1016164 
    ## ... Procrustes: rmse 1.898928e-05  max resid 5.348315e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.1016164 
    ## ... Procrustes: rmse 5.658302e-05  max resid 0.0001598553 
    ## ... Similar to previous best
    ## Run 389 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005396572  max resid 0.00150711 
    ## ... Similar to previous best
    ## Run 390 stress 0.1341868 
    ## Run 391 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005282137  max resid 0.00147805 
    ## ... Similar to previous best
    ## Run 392 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003833371  max resid 0.001071113 
    ## ... Similar to previous best
    ## Run 393 stress 0.1341868 
    ## Run 394 stress 0.1016164 
    ## ... Procrustes: rmse 1.485009e-06  max resid 4.169041e-06 
    ## ... Similar to previous best
    ## Run 395 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006154503  max resid 0.001718922 
    ## ... Similar to previous best
    ## Run 396 stress 0.1016164 
    ## ... Procrustes: rmse 2.090442e-06  max resid 5.68955e-06 
    ## ... Similar to previous best
    ## Run 397 stress 0.1016164 
    ## ... Procrustes: rmse 2.345422e-05  max resid 6.620333e-05 
    ## ... Similar to previous best
    ## Run 398 stress 0.1016164 
    ## ... Procrustes: rmse 3.841028e-05  max resid 0.00010842 
    ## ... Similar to previous best
    ## Run 399 stress 0.1016171 
    ## ... Procrustes: rmse 0.0006831821  max resid 0.001907069 
    ## ... Similar to previous best
    ## Run 400 stress 0.1016164 
    ## ... Procrustes: rmse 5.113427e-05  max resid 0.0001415603 
    ## ... Similar to previous best
    ## Run 401 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005018839  max resid 0.001403175 
    ## ... Similar to previous best
    ## Run 402 stress 0.1341868 
    ## Run 403 stress 0.1341868 
    ## Run 404 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003285022  max resid 0.0009158778 
    ## ... Similar to previous best
    ## Run 405 stress 0.1356147 
    ## Run 406 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002965853  max resid 0.0008310009 
    ## ... Similar to previous best
    ## Run 407 stress 0.1341868 
    ## Run 408 stress 0.1341868 
    ## Run 409 stress 0.1356145 
    ## Run 410 stress 0.101617 
    ## ... Procrustes: rmse 0.0007012612  max resid 0.001958769 
    ## ... Similar to previous best
    ## Run 411 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001185431  max resid 0.000334054 
    ## ... Similar to previous best
    ## Run 412 stress 0.1016164 
    ## ... Procrustes: rmse 4.873776e-05  max resid 0.0001374034 
    ## ... Similar to previous best
    ## Run 413 stress 0.3202131 
    ## Run 414 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004946111  max resid 0.001383058 
    ## ... Similar to previous best
    ## Run 415 stress 0.1016164 
    ## ... Procrustes: rmse 3.146747e-05  max resid 8.806581e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006219563  max resid 0.001737848 
    ## ... Similar to previous best
    ## Run 417 stress 0.1016164 
    ## ... Procrustes: rmse 7.4636e-06  max resid 2.097898e-05 
    ## ... Similar to previous best
    ## Run 418 stress 0.1016164 
    ## ... Procrustes: rmse 6.422508e-05  max resid 0.0001814746 
    ## ... Similar to previous best
    ## Run 419 stress 0.1016164 
    ## ... Procrustes: rmse 3.994505e-05  max resid 0.0001120584 
    ## ... Similar to previous best
    ## Run 420 stress 0.1016164 
    ## ... Procrustes: rmse 9.036176e-05  max resid 0.0002551186 
    ## ... Similar to previous best
    ## Run 421 stress 0.1016164 
    ## ... Procrustes: rmse 3.999319e-05  max resid 0.0001121689 
    ## ... Similar to previous best
    ## Run 422 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004329576  max resid 0.001211355 
    ## ... Similar to previous best
    ## Run 423 stress 0.1016164 
    ## ... Procrustes: rmse 3.146441e-05  max resid 8.801338e-05 
    ## ... Similar to previous best
    ## Run 424 stress 0.1356151 
    ## Run 425 stress 0.1016164 
    ## ... Procrustes: rmse 3.045843e-05  max resid 8.402927e-05 
    ## ... Similar to previous best
    ## Run 426 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004180186  max resid 0.001165615 
    ## ... Similar to previous best
    ## Run 427 stress 0.1016164 
    ## ... Procrustes: rmse 3.341694e-05  max resid 9.420057e-05 
    ## ... Similar to previous best
    ## Run 428 stress 0.1341868 
    ## Run 429 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006303683  max resid 0.00176143 
    ## ... Similar to previous best
    ## Run 430 stress 0.1016164 
    ## ... Procrustes: rmse 3.850714e-06  max resid 1.066808e-05 
    ## ... Similar to previous best
    ## Run 431 stress 0.1016164 
    ## ... Procrustes: rmse 5.559404e-05  max resid 0.0001564339 
    ## ... Similar to previous best
    ## Run 432 stress 0.1016164 
    ## ... Procrustes: rmse 6.539366e-05  max resid 0.0001821207 
    ## ... Similar to previous best
    ## Run 433 stress 0.1016164 
    ## ... Procrustes: rmse 8.831888e-05  max resid 0.0002487795 
    ## ... Similar to previous best
    ## Run 434 stress 0.1016164 
    ## ... Procrustes: rmse 1.350403e-06  max resid 2.755201e-06 
    ## ... Similar to previous best
    ## Run 435 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004207683  max resid 0.001176524 
    ## ... Similar to previous best
    ## Run 436 stress 0.1016164 
    ## ... Procrustes: rmse 9.179728e-05  max resid 0.0002571004 
    ## ... Similar to previous best
    ## Run 437 stress 0.1016164 
    ## ... Procrustes: rmse 9.996341e-05  max resid 0.0002810149 
    ## ... Similar to previous best
    ## Run 438 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004901438  max resid 0.001369122 
    ## ... Similar to previous best
    ## Run 439 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006322935  max resid 0.001766435 
    ## ... Similar to previous best
    ## Run 440 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003107692  max resid 0.0008665749 
    ## ... Similar to previous best
    ## Run 441 stress 0.1016164 
    ## ... Procrustes: rmse 4.651159e-05  max resid 0.0001293808 
    ## ... Similar to previous best
    ## Run 442 stress 0.1016164 
    ## ... Procrustes: rmse 8.91527e-05  max resid 0.0002512122 
    ## ... Similar to previous best
    ## Run 443 stress 0.101617 
    ## ... Procrustes: rmse 0.0006662576  max resid 0.001862161 
    ## ... Similar to previous best
    ## Run 444 stress 0.1341868 
    ## Run 445 stress 0.1016164 
    ## ... Procrustes: rmse 4.696656e-05  max resid 0.0001321817 
    ## ... Similar to previous best
    ## Run 446 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003067653  max resid 0.0008591611 
    ## ... Similar to previous best
    ## Run 447 stress 0.1016164 
    ## ... Procrustes: rmse 8.199518e-05  max resid 0.0002297419 
    ## ... Similar to previous best
    ## Run 448 stress 0.1356147 
    ## Run 449 stress 0.1341868 
    ## Run 450 stress 0.1356148 
    ## Run 451 stress 0.1016164 
    ## ... Procrustes: rmse 3.027511e-05  max resid 8.479767e-05 
    ## ... Similar to previous best
    ## Run 452 stress 0.1341868 
    ## Run 453 stress 0.1341868 
    ## Run 454 stress 0.1016164 
    ## ... Procrustes: rmse 8.696998e-06  max resid 2.131773e-05 
    ## ... Similar to previous best
    ## Run 455 stress 0.1016164 
    ## ... Procrustes: rmse 1.202665e-05  max resid 3.383627e-05 
    ## ... Similar to previous best
    ## Run 456 stress 0.1016164 
    ## ... Procrustes: rmse 3.51578e-05  max resid 9.887078e-05 
    ## ... Similar to previous best
    ## Run 457 stress 0.1356149 
    ## Run 458 stress 0.1341868 
    ## Run 459 stress 0.1016164 
    ## ... Procrustes: rmse 4.386422e-05  max resid 0.0001221489 
    ## ... Similar to previous best
    ## Run 460 stress 0.1016164 
    ## ... Procrustes: rmse 7.220102e-05  max resid 0.000203004 
    ## ... Similar to previous best
    ## Run 461 stress 0.1016164 
    ## ... Procrustes: rmse 6.103914e-05  max resid 0.0001714774 
    ## ... Similar to previous best
    ## Run 462 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002401659  max resid 0.0006744558 
    ## ... Similar to previous best
    ## Run 463 stress 0.1016164 
    ## ... Procrustes: rmse 1.035636e-05  max resid 2.925911e-05 
    ## ... Similar to previous best
    ## Run 464 stress 0.1016164 
    ## ... Procrustes: rmse 5.327063e-05  max resid 0.0001497252 
    ## ... Similar to previous best
    ## Run 465 stress 0.1016164 
    ## ... Procrustes: rmse 3.634043e-05  max resid 0.0001023774 
    ## ... Similar to previous best
    ## Run 466 stress 0.1356143 
    ## Run 467 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005716694  max resid 0.001593399 
    ## ... Similar to previous best
    ## Run 468 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001019204  max resid 0.0002866113 
    ## ... Similar to previous best
    ## Run 469 stress 0.1016164 
    ## ... Procrustes: rmse 2.36822e-05  max resid 6.678056e-05 
    ## ... Similar to previous best
    ## Run 470 stress 0.1016164 
    ## ... Procrustes: rmse 1.578685e-05  max resid 3.422783e-05 
    ## ... Similar to previous best
    ## Run 471 stress 0.1356146 
    ## Run 472 stress 0.135615 
    ## Run 473 stress 0.1016164 
    ## ... Procrustes: rmse 9.440517e-05  max resid 0.0002646114 
    ## ... Similar to previous best
    ## Run 474 stress 0.1016164 
    ## ... Procrustes: rmse 4.101775e-06  max resid 1.160849e-05 
    ## ... Similar to previous best
    ## Run 475 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003513868  max resid 0.0009838365 
    ## ... Similar to previous best
    ## Run 476 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004453181  max resid 0.001245327 
    ## ... Similar to previous best
    ## Run 477 stress 0.1016164 
    ## ... Procrustes: rmse 9.236566e-05  max resid 0.0002608731 
    ## ... Similar to previous best
    ## Run 478 stress 0.1016164 
    ## ... Procrustes: rmse 1.994843e-05  max resid 4.822888e-05 
    ## ... Similar to previous best
    ## Run 479 stress 0.1016164 
    ## ... Procrustes: rmse 7.65934e-06  max resid 2.152637e-05 
    ## ... Similar to previous best
    ## Run 480 stress 0.1530422 
    ## Run 481 stress 0.101617 
    ## ... Procrustes: rmse 0.0006624571  max resid 0.001851513 
    ## ... Similar to previous best
    ## Run 482 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002659289  max resid 0.0007453671 
    ## ... Similar to previous best
    ## Run 483 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003886169  max resid 0.001084672 
    ## ... Similar to previous best
    ## Run 484 stress 0.1341868 
    ## Run 485 stress 0.1341868 
    ## Run 486 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004997138  max resid 0.001397332 
    ## ... Similar to previous best
    ## Run 487 stress 0.1016164 
    ## ... Procrustes: rmse 2.377311e-05  max resid 6.694734e-05 
    ## ... Similar to previous best
    ## Run 488 stress 0.1016169 
    ## ... Procrustes: rmse 0.000552247  max resid 0.001541297 
    ## ... Similar to previous best
    ## Run 489 stress 0.1016164 
    ## ... Procrustes: rmse 1.618499e-05  max resid 4.550789e-05 
    ## ... Similar to previous best
    ## Run 490 stress 0.1341868 
    ## Run 491 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001338341  max resid 0.000375736 
    ## ... Similar to previous best
    ## Run 492 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006024157  max resid 0.0016839 
    ## ... Similar to previous best
    ## Run 493 stress 0.3495547 
    ## Run 494 stress 0.1016164 
    ## ... Procrustes: rmse 1.890191e-05  max resid 5.324978e-05 
    ## ... Similar to previous best
    ## Run 495 stress 0.1356151 
    ## Run 496 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001105414  max resid 0.0003119913 
    ## ... Similar to previous best
    ## Run 497 stress 0.1016164 
    ## ... Procrustes: rmse 3.458504e-05  max resid 9.745845e-05 
    ## ... Similar to previous best
    ## Run 498 stress 0.1016164 
    ## ... Procrustes: rmse 7.168532e-05  max resid 0.0002024648 
    ## ... Similar to previous best
    ## Run 499 stress 0.1016164 
    ## ... Procrustes: rmse 2.801607e-05  max resid 7.767857e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.1341868 
    ## *** Best solution repeated 264 times

``` r
# Stratified lakes and ocean sites
PD_beta_SO_NMDS <- metaMDS(PD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.04851417 
    ## Run 1 stress 0.05643261 
    ## Run 2 stress 0.05607344 
    ## Run 3 stress 0.05643261 
    ## Run 4 stress 0.04851417 
    ## ... New best solution
    ## ... Procrustes: rmse 2.093446e-05  max resid 5.789737e-05 
    ## ... Similar to previous best
    ## Run 5 stress 0.04851417 
    ## ... Procrustes: rmse 2.384185e-05  max resid 6.292011e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.05959895 
    ## Run 7 stress 0.0561224 
    ## Run 8 stress 0.06400741 
    ## Run 9 stress 0.05488103 
    ## Run 10 stress 0.05605147 
    ## Run 11 stress 0.05728074 
    ## Run 12 stress 0.05612241 
    ## Run 13 stress 0.05605149 
    ## Run 14 stress 0.06205753 
    ## Run 15 stress 0.06088971 
    ## Run 16 stress 0.04935948 
    ## Run 17 stress 0.04851417 
    ## ... Procrustes: rmse 2.659224e-05  max resid 7.353441e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.06291424 
    ## Run 19 stress 0.05216149 
    ## Run 20 stress 0.05216132 
    ## Run 21 stress 0.04938616 
    ## Run 22 stress 0.05607345 
    ## Run 23 stress 0.05276411 
    ## Run 24 stress 0.05607347 
    ## Run 25 stress 0.05605155 
    ## Run 26 stress 0.05728067 
    ## Run 27 stress 0.05607344 
    ## Run 28 stress 0.06400732 
    ## Run 29 stress 0.05695861 
    ## Run 30 stress 0.05643262 
    ## Run 31 stress 0.06205753 
    ## Run 32 stress 0.04851419 
    ## ... Procrustes: rmse 5.263368e-05  max resid 0.0001435711 
    ## ... Similar to previous best
    ## Run 33 stress 0.04851417 
    ## ... New best solution
    ## ... Procrustes: rmse 3.369357e-06  max resid 6.297015e-06 
    ## ... Similar to previous best
    ## Run 34 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04982653  max resid 0.139517 
    ## Run 35 stress 0.04935942 
    ## Run 36 stress 0.05572695 
    ## Run 37 stress 0.0485142 
    ## Run 38 stress 0.04938617 
    ## Run 39 stress 0.04935942 
    ## Run 40 stress 0.05643261 
    ## Run 41 stress 0.05572688 
    ## Run 42 stress 0.0543895 
    ## Run 43 stress 0.05643262 
    ## Run 44 stress 0.05855031 
    ## Run 45 stress 0.0543928 
    ## Run 46 stress 0.05433178 
    ## Run 47 stress 0.06400727 
    ## Run 48 stress 0.05216147 
    ## Run 49 stress 0.05855039 
    ## Run 50 stress 0.04851418 
    ## Run 51 stress 0.05959897 
    ## Run 52 stress 0.05572689 
    ## Run 53 stress 0.04851418 
    ## Run 54 stress 0.05959897 
    ## Run 55 stress 0.05593462 
    ## Run 56 stress 0.04938615 
    ## Run 57 stress 0.06933831 
    ## Run 58 stress 0.04938615 
    ## Run 59 stress 0.05474817 
    ## Run 60 stress 0.04935945 
    ## Run 61 stress 0.04851417 
    ## Run 62 stress 0.05607348 
    ## Run 63 stress 0.06185044 
    ## Run 64 stress 0.0543894 
    ## Run 65 stress 0.05643263 
    ## Run 66 stress 0.06291427 
    ## Run 67 stress 0.05438939 
    ## Run 68 stress 0.04851417 
    ## Run 69 stress 0.05728063 
    ## Run 70 stress 0.06291424 
    ## Run 71 stress 0.05643271 
    ## Run 72 stress 0.05605149 
    ## Run 73 stress 0.0521615 
    ## Run 74 stress 0.05607344 
    ## Run 75 stress 0.05607345 
    ## Run 76 stress 0.0569587 
    ## Run 77 stress 0.04860422 
    ## Run 78 stress 0.05216152 
    ## Run 79 stress 0.05447775 
    ## Run 80 stress 0.04938619 
    ## Run 81 stress 0.05643261 
    ## Run 82 stress 0.06088972 
    ## Run 83 stress 0.05137967 
    ## Run 84 stress 0.04938616 
    ## Run 85 stress 0.05842815 
    ## Run 86 stress 0.04938616 
    ## Run 87 stress 0.05605146 
    ## Run 88 stress 0.05728052 
    ## Run 89 stress 0.05451493 
    ## Run 90 stress 0.05643263 
    ## Run 91 stress 0.04851417 
    ## Run 92 stress 0.05842786 
    ## Run 93 stress 0.05695869 
    ## Run 94 stress 0.0493862 
    ## Run 95 stress 0.06933824 
    ## Run 96 stress 0.05216143 
    ## Run 97 stress 0.05439279 
    ## Run 98 stress 0.05216132 
    ## Run 99 stress 0.04935953 
    ## Run 100 stress 0.059599 
    ## Run 101 stress 0.349556 
    ## Run 102 stress 0.05593453 
    ## Run 103 stress 0.05276411 
    ## Run 104 stress 0.05607347 
    ## Run 105 stress 0.06185047 
    ## Run 106 stress 0.05474816 
    ## Run 107 stress 0.05728061 
    ## Run 108 stress 0.05216143 
    ## Run 109 stress 0.06205753 
    ## Run 110 stress 0.05438941 
    ## Run 111 stress 0.04851417 
    ## Run 112 stress 0.04851417 
    ## Run 113 stress 0.05216136 
    ## Run 114 stress 0.0493596 
    ## Run 115 stress 0.06185057 
    ## Run 116 stress 0.05855038 
    ## Run 117 stress 0.04851417 
    ## Run 118 stress 0.06205753 
    ## Run 119 stress 0.05438944 
    ## Run 120 stress 0.05728053 
    ## Run 121 stress 0.05593451 
    ## Run 122 stress 0.04851417 
    ## Run 123 stress 0.0559346 
    ## Run 124 stress 0.05607346 
    ## Run 125 stress 0.05855031 
    ## Run 126 stress 0.04935955 
    ## Run 127 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 4.225813e-05  max resid 7.142555e-05 
    ## ... Similar to previous best
    ## Run 128 stress 0.05959897 
    ## Run 129 stress 0.05695861 
    ## Run 130 stress 0.05607344 
    ## Run 131 stress 0.05137965 
    ## Run 132 stress 0.05605147 
    ## Run 133 stress 0.04938615 
    ## Run 134 stress 0.05643263 
    ## Run 135 stress 0.05474817 
    ## Run 136 stress 0.06088971 
    ## Run 137 stress 0.05216152 
    ## Run 138 stress 0.06400744 
    ## Run 139 stress 0.05643265 
    ## Run 140 stress 0.06291427 
    ## Run 141 stress 0.0527641 
    ## Run 142 stress 0.05451494 
    ## Run 143 stress 0.3398639 
    ## Run 144 stress 0.06088971 
    ## Run 145 stress 0.05572689 
    ## Run 146 stress 0.05438943 
    ## Run 147 stress 0.0543895 
    ## Run 148 stress 0.05607347 
    ## Run 149 stress 0.04938616 
    ## Run 150 stress 0.05488104 
    ## Run 151 stress 0.04938616 
    ## Run 152 stress 0.06088974 
    ## Run 153 stress 0.04935943 
    ## Run 154 stress 0.0527641 
    ## Run 155 stress 0.04938616 
    ## Run 156 stress 0.06205754 
    ## Run 157 stress 0.05276413 
    ## Run 158 stress 0.05612246 
    ## Run 159 stress 0.05612246 
    ## Run 160 stress 0.04938616 
    ## Run 161 stress 0.04935946 
    ## Run 162 stress 0.04938617 
    ## Run 163 stress 0.05572689 
    ## Run 164 stress 0.06185047 
    ## Run 165 stress 0.06291425 
    ## Run 166 stress 0.04612758 
    ## ... Procrustes: rmse 6.713554e-05  max resid 0.0001275938 
    ## ... Similar to previous best
    ## Run 167 stress 0.05216127 
    ## Run 168 stress 0.05643261 
    ## Run 169 stress 0.0559345 
    ## Run 170 stress 0.05855033 
    ## Run 171 stress 0.04612765 
    ## ... Procrustes: rmse 0.0001417915  max resid 0.0002814563 
    ## ... Similar to previous best
    ## Run 172 stress 0.06088979 
    ## Run 173 stress 0.04935946 
    ## Run 174 stress 0.04851417 
    ## Run 175 stress 0.04860423 
    ## Run 176 stress 0.05728054 
    ## Run 177 stress 0.05607344 
    ## Run 178 stress 0.0493862 
    ## Run 179 stress 0.05643265 
    ## Run 180 stress 0.3582285 
    ## Run 181 stress 0.04938617 
    ## Run 182 stress 0.06205752 
    ## Run 183 stress 0.05260049 
    ## Run 184 stress 0.05728056 
    ## Run 185 stress 0.05216129 
    ## Run 186 stress 0.04851417 
    ## Run 187 stress 0.05137966 
    ## Run 188 stress 0.05728057 
    ## Run 189 stress 0.05728096 
    ## Run 190 stress 0.05137963 
    ## Run 191 stress 0.04851417 
    ## Run 192 stress 0.3495702 
    ## Run 193 stress 0.05695871 
    ## Run 194 stress 0.05451493 
    ## Run 195 stress 0.04851418 
    ## Run 196 stress 0.04938617 
    ## Run 197 stress 0.05451498 
    ## Run 198 stress 0.04612765 
    ## ... Procrustes: rmse 0.0001424714  max resid 0.0002845256 
    ## ... Similar to previous best
    ## Run 199 stress 0.05488104 
    ## Run 200 stress 0.05433177 
    ## Run 201 stress 0.04860422 
    ## Run 202 stress 0.04612759 
    ## ... Procrustes: rmse 7.527016e-05  max resid 0.0001429649 
    ## ... Similar to previous best
    ## Run 203 stress 0.0557269 
    ## Run 204 stress 0.04851417 
    ## Run 205 stress 0.04938617 
    ## Run 206 stress 0.04851418 
    ## Run 207 stress 0.06291426 
    ## Run 208 stress 0.04938616 
    ## Run 209 stress 0.04851421 
    ## Run 210 stress 0.0585503 
    ## Run 211 stress 0.0527641 
    ## Run 212 stress 0.05607351 
    ## Run 213 stress 0.04938615 
    ## Run 214 stress 0.05607345 
    ## Run 215 stress 0.05572689 
    ## Run 216 stress 0.04938617 
    ## Run 217 stress 0.05260057 
    ## Run 218 stress 0.04938615 
    ## Run 219 stress 0.05695875 
    ## Run 220 stress 0.05216141 
    ## Run 221 stress 0.05643262 
    ## Run 222 stress 0.04851418 
    ## Run 223 stress 0.04935947 
    ## Run 224 stress 0.05728055 
    ## Run 225 stress 0.3305548 
    ## Run 226 stress 0.05488108 
    ## Run 227 stress 0.0527641 
    ## Run 228 stress 0.04851418 
    ## Run 229 stress 0.06088972 
    ## Run 230 stress 0.05855031 
    ## Run 231 stress 0.04938618 
    ## Run 232 stress 0.04938615 
    ## Run 233 stress 0.04938617 
    ## Run 234 stress 0.04860423 
    ## Run 235 stress 0.04851417 
    ## Run 236 stress 0.05260048 
    ## Run 237 stress 0.05728052 
    ## Run 238 stress 0.04851417 
    ## Run 239 stress 0.05474815 
    ## Run 240 stress 0.05216132 
    ## Run 241 stress 0.06291431 
    ## Run 242 stress 0.04938616 
    ## Run 243 stress 0.05216143 
    ## Run 244 stress 0.0693381 
    ## Run 245 stress 0.0527641 
    ## Run 246 stress 0.05855039 
    ## Run 247 stress 0.0526005 
    ## Run 248 stress 0.04851417 
    ## Run 249 stress 0.04851417 
    ## Run 250 stress 0.06205754 
    ## Run 251 stress 0.05137965 
    ## Run 252 stress 0.05474814 
    ## Run 253 stress 0.05474817 
    ## Run 254 stress 0.05607349 
    ## Run 255 stress 0.0527641 
    ## Run 256 stress 0.04938619 
    ## Run 257 stress 0.05695875 
    ## Run 258 stress 0.04860423 
    ## Run 259 stress 0.05605152 
    ## Run 260 stress 0.05216131 
    ## Run 261 stress 0.04860422 
    ## Run 262 stress 0.04851418 
    ## Run 263 stress 0.05728057 
    ## Run 264 stress 0.05612241 
    ## Run 265 stress 0.05728069 
    ## Run 266 stress 0.04851418 
    ## Run 267 stress 0.05216162 
    ## Run 268 stress 0.05855034 
    ## Run 269 stress 0.04935951 
    ## Run 270 stress 0.05260054 
    ## Run 271 stress 0.05572688 
    ## Run 272 stress 0.05643261 
    ## Run 273 stress 0.04938615 
    ## Run 274 stress 0.0560515 
    ## Run 275 stress 0.0543928 
    ## Run 276 stress 0.04851417 
    ## Run 277 stress 0.05959901 
    ## Run 278 stress 0.04612758 
    ## ... Procrustes: rmse 5.608822e-05  max resid 0.0001021161 
    ## ... Similar to previous best
    ## Run 279 stress 0.04938616 
    ## Run 280 stress 0.05474828 
    ## Run 281 stress 0.06185066 
    ## Run 282 stress 0.04851417 
    ## Run 283 stress 0.05643263 
    ## Run 284 stress 0.05216157 
    ## Run 285 stress 0.3314337 
    ## Run 286 stress 0.05593456 
    ## Run 287 stress 0.04938616 
    ## Run 288 stress 0.05260049 
    ## Run 289 stress 0.05607345 
    ## Run 290 stress 0.05607345 
    ## Run 291 stress 0.06185067 
    ## Run 292 stress 0.05137965 
    ## Run 293 stress 0.05488104 
    ## Run 294 stress 0.04851417 
    ## Run 295 stress 0.05433178 
    ## Run 296 stress 0.04851418 
    ## Run 297 stress 0.05607345 
    ## Run 298 stress 0.05695876 
    ## Run 299 stress 0.06291427 
    ## Run 300 stress 0.05695855 
    ## Run 301 stress 0.06088972 
    ## Run 302 stress 0.04938615 
    ## Run 303 stress 0.04851417 
    ## Run 304 stress 0.05607344 
    ## Run 305 stress 0.05451493 
    ## Run 306 stress 0.04860422 
    ## Run 307 stress 0.05605148 
    ## Run 308 stress 0.05593461 
    ## Run 309 stress 0.0527641 
    ## Run 310 stress 0.05276411 
    ## Run 311 stress 0.05695869 
    ## Run 312 stress 0.05216136 
    ## Run 313 stress 0.04938616 
    ## Run 314 stress 0.04938616 
    ## Run 315 stress 0.06185046 
    ## Run 316 stress 0.06185048 
    ## Run 317 stress 0.04935946 
    ## Run 318 stress 0.07046337 
    ## Run 319 stress 0.0585503 
    ## Run 320 stress 0.05842787 
    ## Run 321 stress 0.06400736 
    ## Run 322 stress 0.0485142 
    ## Run 323 stress 0.06205753 
    ## Run 324 stress 0.05137965 
    ## Run 325 stress 0.05607347 
    ## Run 326 stress 0.05607346 
    ## Run 327 stress 0.05593452 
    ## Run 328 stress 0.05612243 
    ## Run 329 stress 0.05605148 
    ## Run 330 stress 0.3398585 
    ## Run 331 stress 0.05137965 
    ## Run 332 stress 0.04938616 
    ## Run 333 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 3.408226e-05  max resid 6.132798e-05 
    ## ... Similar to previous best
    ## Run 334 stress 0.05695854 
    ## Run 335 stress 0.04938615 
    ## Run 336 stress 0.05728053 
    ## Run 337 stress 0.05695853 
    ## Run 338 stress 0.05643261 
    ## Run 339 stress 0.05137965 
    ## Run 340 stress 0.05728058 
    ## Run 341 stress 0.0493862 
    ## Run 342 stress 0.04860422 
    ## Run 343 stress 0.05137974 
    ## Run 344 stress 0.04851417 
    ## Run 345 stress 0.05605152 
    ## Run 346 stress 0.04851417 
    ## Run 347 stress 0.05572688 
    ## Run 348 stress 0.04938615 
    ## Run 349 stress 0.05474818 
    ## Run 350 stress 0.04860422 
    ## Run 351 stress 0.06088974 
    ## Run 352 stress 0.05842802 
    ## Run 353 stress 0.05438944 
    ## Run 354 stress 0.05216133 
    ## Run 355 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001044764  max resid 0.0002100922 
    ## ... Similar to previous best
    ## Run 356 stress 0.05728058 
    ## Run 357 stress 0.05216127 
    ## Run 358 stress 0.05728062 
    ## Run 359 stress 0.0543928 
    ## Run 360 stress 0.05855034 
    ## Run 361 stress 0.04935945 
    ## Run 362 stress 0.0543894 
    ## Run 363 stress 0.05612243 
    ## Run 364 stress 0.05643261 
    ## Run 365 stress 0.04938618 
    ## Run 366 stress 0.0585503 
    ## Run 367 stress 0.04612762 
    ## ... Procrustes: rmse 4.822e-05  max resid 8.431061e-05 
    ## ... Similar to previous best
    ## Run 368 stress 0.04851418 
    ## Run 369 stress 0.05474813 
    ## Run 370 stress 0.05959895 
    ## Run 371 stress 0.05842782 
    ## Run 372 stress 0.04851417 
    ## Run 373 stress 0.04860423 
    ## Run 374 stress 0.06933817 
    ## Run 375 stress 0.05855036 
    ## Run 376 stress 0.06400751 
    ## Run 377 stress 0.04938616 
    ## Run 378 stress 0.04612762 
    ## ... Procrustes: rmse 9.640881e-05  max resid 0.0001833685 
    ## ... Similar to previous best
    ## Run 379 stress 0.05643264 
    ## Run 380 stress 0.05728082 
    ## Run 381 stress 0.04612768 
    ## ... Procrustes: rmse 0.0001294061  max resid 0.0002540531 
    ## ... Similar to previous best
    ## Run 382 stress 0.05959895 
    ## Run 383 stress 0.06205753 
    ## Run 384 stress 0.04938616 
    ## Run 385 stress 0.06205755 
    ## Run 386 stress 0.05605145 
    ## Run 387 stress 0.04938615 
    ## Run 388 stress 0.04938615 
    ## Run 389 stress 0.04935946 
    ## Run 390 stress 0.04860423 
    ## Run 391 stress 0.05276412 
    ## Run 392 stress 0.05855032 
    ## Run 393 stress 0.06205753 
    ## Run 394 stress 0.05605145 
    ## Run 395 stress 0.04851417 
    ## Run 396 stress 0.05643261 
    ## Run 397 stress 0.05260048 
    ## Run 398 stress 0.04851417 
    ## Run 399 stress 0.05260052 
    ## Run 400 stress 0.05474827 
    ## Run 401 stress 0.05959898 
    ## Run 402 stress 0.05216126 
    ## Run 403 stress 0.05474815 
    ## Run 404 stress 0.05728057 
    ## Run 405 stress 0.04938617 
    ## Run 406 stress 0.05276414 
    ## Run 407 stress 0.0521614 
    ## Run 408 stress 0.0461276 
    ## ... Procrustes: rmse 7.650063e-05  max resid 0.0001530286 
    ## ... Similar to previous best
    ## Run 409 stress 0.05260047 
    ## Run 410 stress 0.05695857 
    ## Run 411 stress 0.04938617 
    ## Run 412 stress 0.05643261 
    ## Run 413 stress 0.05438942 
    ## Run 414 stress 0.04938617 
    ## Run 415 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 1.455344e-05  max resid 3.099466e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.0521615 
    ## Run 417 stress 0.05451494 
    ## Run 418 stress 0.04935949 
    ## Run 419 stress 0.05474815 
    ## Run 420 stress 0.04851417 
    ## Run 421 stress 0.04860423 
    ## Run 422 stress 0.05728085 
    ## Run 423 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001259464  max resid 0.0002589016 
    ## ... Similar to previous best
    ## Run 424 stress 0.0543928 
    ## Run 425 stress 0.04935946 
    ## Run 426 stress 0.06205752 
    ## Run 427 stress 0.05959898 
    ## Run 428 stress 0.05260047 
    ## Run 429 stress 0.04612757 
    ## ... Procrustes: rmse 3.92532e-05  max resid 8.013225e-05 
    ## ... Similar to previous best
    ## Run 430 stress 0.05474814 
    ## Run 431 stress 0.04860423 
    ## Run 432 stress 0.05695873 
    ## Run 433 stress 0.06291425 
    ## Run 434 stress 0.0608897 
    ## Run 435 stress 0.05728063 
    ## Run 436 stress 0.04851417 
    ## Run 437 stress 0.05605154 
    ## Run 438 stress 0.05728061 
    ## Run 439 stress 0.0608897 
    ## Run 440 stress 0.06185044 
    ## Run 441 stress 0.04851417 
    ## Run 442 stress 0.05474815 
    ## Run 443 stress 0.0572808 
    ## Run 444 stress 0.05695868 
    ## Run 445 stress 0.05643261 
    ## Run 446 stress 0.04938617 
    ## Run 447 stress 0.06400739 
    ## Run 448 stress 0.05137963 
    ## Run 449 stress 0.04612758 
    ## ... Procrustes: rmse 3.611902e-05  max resid 7.062196e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.05451494 
    ## Run 451 stress 0.05605154 
    ## Run 452 stress 0.05607345 
    ## Run 453 stress 0.05728058 
    ## Run 454 stress 0.04938616 
    ## Run 455 stress 0.04860422 
    ## Run 456 stress 0.04860422 
    ## Run 457 stress 0.05216133 
    ## Run 458 stress 0.04851417 
    ## Run 459 stress 0.04612757 
    ## ... Procrustes: rmse 4.208651e-05  max resid 8.654233e-05 
    ## ... Similar to previous best
    ## Run 460 stress 0.04938619 
    ## Run 461 stress 0.04851417 
    ## Run 462 stress 0.06400744 
    ## Run 463 stress 0.04938616 
    ## Run 464 stress 0.04938615 
    ## Run 465 stress 0.05216156 
    ## Run 466 stress 0.05612247 
    ## Run 467 stress 0.05643261 
    ## Run 468 stress 0.06205755 
    ## Run 469 stress 0.0560515 
    ## Run 470 stress 0.06291437 
    ## Run 471 stress 0.0526006 
    ## Run 472 stress 0.0461276 
    ## ... Procrustes: rmse 9.317295e-05  max resid 0.0001878662 
    ## ... Similar to previous best
    ## Run 473 stress 0.05728071 
    ## Run 474 stress 0.05276411 
    ## Run 475 stress 0.05605151 
    ## Run 476 stress 0.05612241 
    ## Run 477 stress 0.05451493 
    ## Run 478 stress 0.05607344 
    ## Run 479 stress 0.04851417 
    ## Run 480 stress 0.05572689 
    ## Run 481 stress 0.04851418 
    ## Run 482 stress 0.06205753 
    ## Run 483 stress 0.05959899 
    ## Run 484 stress 0.04938618 
    ## Run 485 stress 0.0618505 
    ## Run 486 stress 0.0561225 
    ## Run 487 stress 0.06088976 
    ## Run 488 stress 0.06205756 
    ## Run 489 stress 0.06400728 
    ## Run 490 stress 0.05572688 
    ## Run 491 stress 0.04938615 
    ## Run 492 stress 0.05842793 
    ## Run 493 stress 0.05855032 
    ## Run 494 stress 0.05438937 
    ## Run 495 stress 0.06088979 
    ## Run 496 stress 0.05855043 
    ## Run 497 stress 0.05216158 
    ## Run 498 stress 0.05607344 
    ## Run 499 stress 0.05451497 
    ## Run 500 stress 0.05607345 
    ## *** Best solution repeated 6 times

``` r
# Stratified lakes
PD_beta_S_NMDS <- metaMDS(PD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.02693829 
    ## Run 1 stress 0.02693831 
    ## ... Procrustes: rmse 5.024783e-05  max resid 8.56228e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.075473 
    ## Run 3 stress 0.02693839 
    ## ... Procrustes: rmse 0.0002055306  max resid 0.000352324 
    ## ... Similar to previous best
    ## Run 4 stress 0.1756177 
    ## Run 5 stress 0.02693833 
    ## ... Procrustes: rmse 8.77519e-05  max resid 0.0001495885 
    ## ... Similar to previous best
    ## Run 6 stress 0.02693829 
    ## ... Procrustes: rmse 6.744687e-05  max resid 0.0001156995 
    ## ... Similar to previous best
    ## Run 7 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001818531  max resid 0.0003161972 
    ## ... Similar to previous best
    ## Run 8 stress 0.07374802 
    ## Run 9 stress 0.08382861 
    ## Run 10 stress 0.073748 
    ## Run 11 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001440904  max resid 0.0002461028 
    ## ... Similar to previous best
    ## Run 12 stress 0.07374811 
    ## Run 13 stress 0.07374803 
    ## Run 14 stress 0.07547299 
    ## Run 15 stress 0.1799051 
    ## Run 16 stress 0.02693832 
    ## ... Procrustes: rmse 0.00012713  max resid 0.000218895 
    ## ... Similar to previous best
    ## Run 17 stress 0.073748 
    ## Run 18 stress 0.07374806 
    ## Run 19 stress 0.07547299 
    ## Run 20 stress 0.02693831 
    ## ... Procrustes: rmse 4.282698e-05  max resid 7.386016e-05 
    ## ... Similar to previous best
    ## Run 21 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001424447  max resid 0.0002441652 
    ## ... Similar to previous best
    ## Run 22 stress 0.08382861 
    ## Run 23 stress 0.07374809 
    ## Run 24 stress 0.07374804 
    ## Run 25 stress 0.08382825 
    ## Run 26 stress 0.1756177 
    ## Run 27 stress 0.07374802 
    ## Run 28 stress 0.02693833 
    ## ... Procrustes: rmse 8.108494e-05  max resid 0.0001389356 
    ## ... Similar to previous best
    ## Run 29 stress 0.07374806 
    ## Run 30 stress 0.02693831 
    ## ... Procrustes: rmse 0.0001148837  max resid 0.0001960359 
    ## ... Similar to previous best
    ## Run 31 stress 0.02693828 
    ## ... New best solution
    ## ... Procrustes: rmse 4.376548e-05  max resid 7.464732e-05 
    ## ... Similar to previous best
    ## Run 32 stress 0.08382865 
    ## Run 33 stress 0.073748 
    ## Run 34 stress 0.07374801 
    ## Run 35 stress 0.1862384 
    ## Run 36 stress 0.1756177 
    ## Run 37 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001901078  max resid 0.0003253787 
    ## ... Similar to previous best
    ## Run 38 stress 0.2852159 
    ## Run 39 stress 0.2852156 
    ## Run 40 stress 0.07374805 
    ## Run 41 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001235006  max resid 0.0002091632 
    ## ... Similar to previous best
    ## Run 42 stress 0.07374807 
    ## Run 43 stress 0.07374803 
    ## Run 44 stress 0.2273964 
    ## Run 45 stress 0.07547299 
    ## Run 46 stress 0.07547299 
    ## Run 47 stress 0.07374803 
    ## Run 48 stress 0.08382863 
    ## Run 49 stress 0.07374804 
    ## Run 50 stress 0.02693842 
    ## ... Procrustes: rmse 0.0001643439  max resid 0.0002790616 
    ## ... Similar to previous best
    ## Run 51 stress 0.1756177 
    ## Run 52 stress 0.07374808 
    ## Run 53 stress 0.0269383 
    ## ... Procrustes: rmse 5.003099e-05  max resid 8.842263e-05 
    ## ... Similar to previous best
    ## Run 54 stress 0.07374802 
    ## Run 55 stress 0.08382826 
    ## Run 56 stress 0.0269383 
    ## ... Procrustes: rmse 4.133206e-05  max resid 7.066521e-05 
    ## ... Similar to previous best
    ## Run 57 stress 0.07374808 
    ## Run 58 stress 0.02693837 
    ## ... Procrustes: rmse 0.000175735  max resid 0.0002999879 
    ## ... Similar to previous best
    ## Run 59 stress 0.3083098 
    ## Run 60 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001645596  max resid 0.0002812909 
    ## ... Similar to previous best
    ## Run 61 stress 0.07374809 
    ## Run 62 stress 0.179905 
    ## Run 63 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001778774  max resid 0.0003060218 
    ## ... Similar to previous best
    ## Run 64 stress 0.02693842 
    ## ... Procrustes: rmse 0.000197633  max resid 0.0003396727 
    ## ... Similar to previous best
    ## Run 65 stress 0.3083098 
    ## Run 66 stress 0.2273964 
    ## Run 67 stress 0.02693833 
    ## ... Procrustes: rmse 9.173539e-05  max resid 0.0001558127 
    ## ... Similar to previous best
    ## Run 68 stress 0.07547299 
    ## Run 69 stress 0.07374801 
    ## Run 70 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001307164  max resid 0.0002231215 
    ## ... Similar to previous best
    ## Run 71 stress 0.179905 
    ## Run 72 stress 0.07374801 
    ## Run 73 stress 0.07374808 
    ## Run 74 stress 0.07374801 
    ## Run 75 stress 0.02693829 
    ## ... Procrustes: rmse 4.680348e-05  max resid 7.961884e-05 
    ## ... Similar to previous best
    ## Run 76 stress 0.07374804 
    ## Run 77 stress 0.07547299 
    ## Run 78 stress 0.0269383 
    ## ... Procrustes: rmse 8.426673e-05  max resid 0.0001442731 
    ## ... Similar to previous best
    ## Run 79 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001041034  max resid 0.0001789416 
    ## ... Similar to previous best
    ## Run 80 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001374653  max resid 0.0002357705 
    ## ... Similar to previous best
    ## Run 81 stress 0.1756177 
    ## Run 82 stress 0.02693841 
    ## ... Procrustes: rmse 0.0001785794  max resid 0.0003067475 
    ## ... Similar to previous best
    ## Run 83 stress 0.179905 
    ## Run 84 stress 0.07374808 
    ## Run 85 stress 0.08382854 
    ## Run 86 stress 0.02693829 
    ## ... Procrustes: rmse 3.046265e-05  max resid 5.192076e-05 
    ## ... Similar to previous best
    ## Run 87 stress 0.179905 
    ## Run 88 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001555319  max resid 0.0002654703 
    ## ... Similar to previous best
    ## Run 89 stress 0.1862384 
    ## Run 90 stress 0.07374803 
    ## Run 91 stress 0.075473 
    ## Run 92 stress 0.1889141 
    ## Run 93 stress 0.073748 
    ## Run 94 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001347555  max resid 0.0002289026 
    ## ... Similar to previous best
    ## Run 95 stress 0.07547299 
    ## Run 96 stress 0.1756177 
    ## Run 97 stress 0.2852156 
    ## Run 98 stress 0.075473 
    ## Run 99 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001651127  max resid 0.0002833842 
    ## ... Similar to previous best
    ## Run 100 stress 0.07547299 
    ## Run 101 stress 0.07374804 
    ## Run 102 stress 0.0269383 
    ## ... Procrustes: rmse 8.270644e-05  max resid 0.0001420568 
    ## ... Similar to previous best
    ## Run 103 stress 0.07374816 
    ## Run 104 stress 0.1889141 
    ## Run 105 stress 0.02693837 
    ## ... Procrustes: rmse 0.000172516  max resid 0.0002948771 
    ## ... Similar to previous best
    ## Run 106 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001464581  max resid 0.0002506707 
    ## ... Similar to previous best
    ## Run 107 stress 0.1756177 
    ## Run 108 stress 0.075473 
    ## Run 109 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001469669  max resid 0.0002515428 
    ## ... Similar to previous best
    ## Run 110 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001045473  max resid 0.0001808893 
    ## ... Similar to previous best
    ## Run 111 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001046521  max resid 0.0001789996 
    ## ... Similar to previous best
    ## Run 112 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001637476  max resid 0.0002807049 
    ## ... Similar to previous best
    ## Run 113 stress 0.02693829 
    ## ... Procrustes: rmse 1.819368e-05  max resid 3.148615e-05 
    ## ... Similar to previous best
    ## Run 114 stress 0.07547299 
    ## Run 115 stress 0.07547302 
    ## Run 116 stress 0.2626666 
    ## Run 117 stress 0.09391927 
    ## Run 118 stress 0.0737481 
    ## Run 119 stress 0.02693842 
    ## ... Procrustes: rmse 0.0002011541  max resid 0.0003423331 
    ## ... Similar to previous best
    ## Run 120 stress 0.1756177 
    ## Run 121 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001646977  max resid 0.0002825267 
    ## ... Similar to previous best
    ## Run 122 stress 0.075473 
    ## Run 123 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001160022  max resid 0.0002008231 
    ## ... Similar to previous best
    ## Run 124 stress 0.075473 
    ## Run 125 stress 0.075473 
    ## Run 126 stress 0.07547299 
    ## Run 127 stress 0.07374801 
    ## Run 128 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001375387  max resid 0.0002352433 
    ## ... Similar to previous best
    ## Run 129 stress 0.075473 
    ## Run 130 stress 0.02693832 
    ## ... Procrustes: rmse 2.226575e-05  max resid 3.216206e-05 
    ## ... Similar to previous best
    ## Run 131 stress 0.02693831 
    ## ... Procrustes: rmse 5.607722e-05  max resid 9.472211e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.07547299 
    ## Run 133 stress 0.075473 
    ## Run 134 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001008238  max resid 0.0001733445 
    ## ... Similar to previous best
    ## Run 135 stress 0.02693841 
    ## ... Procrustes: rmse 0.0001789806  max resid 0.0003063762 
    ## ... Similar to previous best
    ## Run 136 stress 0.08382829 
    ## Run 137 stress 0.07374801 
    ## Run 138 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001317083  max resid 0.0002252149 
    ## ... Similar to previous best
    ## Run 139 stress 0.02693831 
    ## ... Procrustes: rmse 8.893628e-05  max resid 0.0001344853 
    ## ... Similar to previous best
    ## Run 140 stress 0.09391896 
    ## Run 141 stress 0.02693831 
    ## ... Procrustes: rmse 7.177465e-05  max resid 0.0001232521 
    ## ... Similar to previous best
    ## Run 142 stress 0.08382844 
    ## Run 143 stress 0.08382853 
    ## Run 144 stress 0.073748 
    ## Run 145 stress 0.07374808 
    ## Run 146 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001702639  max resid 0.0002893285 
    ## ... Similar to previous best
    ## Run 147 stress 0.07374803 
    ## Run 148 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001110803  max resid 0.0001899101 
    ## ... Similar to previous best
    ## Run 149 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001042092  max resid 0.0001785168 
    ## ... Similar to previous best
    ## Run 150 stress 0.1756177 
    ## Run 151 stress 0.02693829 
    ## ... Procrustes: rmse 5.839492e-05  max resid 9.968953e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.07374802 
    ## Run 153 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001521034  max resid 0.0002609435 
    ## ... Similar to previous best
    ## Run 154 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001208412  max resid 0.0002071619 
    ## ... Similar to previous best
    ## Run 155 stress 0.1862385 
    ## Run 156 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001719035  max resid 0.000293889 
    ## ... Similar to previous best
    ## Run 157 stress 0.1889141 
    ## Run 158 stress 0.08382829 
    ## Run 159 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001137746  max resid 0.0001940097 
    ## ... Similar to previous best
    ## Run 160 stress 0.07374808 
    ## Run 161 stress 0.07374809 
    ## Run 162 stress 0.07374801 
    ## Run 163 stress 0.179905 
    ## Run 164 stress 0.07547299 
    ## Run 165 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001540293  max resid 0.000263542 
    ## ... Similar to previous best
    ## Run 166 stress 0.07374802 
    ## Run 167 stress 0.02693833 
    ## ... Procrustes: rmse 9.108835e-05  max resid 0.0001551808 
    ## ... Similar to previous best
    ## Run 168 stress 0.02693829 
    ## ... Procrustes: rmse 6.581567e-05  max resid 0.0001121097 
    ## ... Similar to previous best
    ## Run 169 stress 0.1889141 
    ## Run 170 stress 0.02693829 
    ## ... Procrustes: rmse 6.85589e-05  max resid 0.0001178687 
    ## ... Similar to previous best
    ## Run 171 stress 0.02693834 
    ## ... Procrustes: rmse 0.00011225  max resid 0.0001919781 
    ## ... Similar to previous best
    ## Run 172 stress 0.1889141 
    ## Run 173 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001544185  max resid 0.0002685406 
    ## ... Similar to previous best
    ## Run 174 stress 0.1756177 
    ## Run 175 stress 0.07374802 
    ## Run 176 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001000709  max resid 0.0001719459 
    ## ... Similar to previous best
    ## Run 177 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001290864  max resid 0.0002218548 
    ## ... Similar to previous best
    ## Run 178 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001829461  max resid 0.0003120012 
    ## ... Similar to previous best
    ## Run 179 stress 0.08382831 
    ## Run 180 stress 0.08382822 
    ## Run 181 stress 0.179905 
    ## Run 182 stress 0.1756177 
    ## Run 183 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001188195  max resid 0.0002029834 
    ## ... Similar to previous best
    ## Run 184 stress 0.1756177 
    ## Run 185 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001640096  max resid 0.0002811472 
    ## ... Similar to previous best
    ## Run 186 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001434805  max resid 0.0002464572 
    ## ... Similar to previous best
    ## Run 187 stress 0.08382841 
    ## Run 188 stress 0.08382824 
    ## Run 189 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001136689  max resid 0.0001943116 
    ## ... Similar to previous best
    ## Run 190 stress 0.02693841 
    ## ... Procrustes: rmse 0.0001912169  max resid 0.0003249266 
    ## ... Similar to previous best
    ## Run 191 stress 0.1756177 
    ## Run 192 stress 0.073748 
    ## Run 193 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001233954  max resid 0.0002119799 
    ## ... Similar to previous best
    ## Run 194 stress 0.02693829 
    ## ... Procrustes: rmse 3.89399e-05  max resid 6.668608e-05 
    ## ... Similar to previous best
    ## Run 195 stress 0.07374803 
    ## Run 196 stress 0.02693829 
    ## ... Procrustes: rmse 5.196996e-05  max resid 8.797431e-05 
    ## ... Similar to previous best
    ## Run 197 stress 0.08382843 
    ## Run 198 stress 0.07547299 
    ## Run 199 stress 0.02693831 
    ## ... Procrustes: rmse 8.417317e-05  max resid 0.0001429143 
    ## ... Similar to previous best
    ## Run 200 stress 0.09391919 
    ## Run 201 stress 0.07374802 
    ## Run 202 stress 0.07374803 
    ## Run 203 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001162587  max resid 0.0001984996 
    ## ... Similar to previous best
    ## Run 204 stress 0.073748 
    ## Run 205 stress 0.07374809 
    ## Run 206 stress 0.07374801 
    ## Run 207 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001599057  max resid 0.0002725915 
    ## ... Similar to previous best
    ## Run 208 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001708514  max resid 0.0002922894 
    ## ... Similar to previous best
    ## Run 209 stress 0.179905 
    ## Run 210 stress 0.07374812 
    ## Run 211 stress 0.075473 
    ## Run 212 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001152616  max resid 0.0001968675 
    ## ... Similar to previous best
    ## Run 213 stress 0.07374808 
    ## Run 214 stress 0.02693829 
    ## ... Procrustes: rmse 4.091685e-05  max resid 7.063996e-05 
    ## ... Similar to previous best
    ## Run 215 stress 0.07374807 
    ## Run 216 stress 0.07374808 
    ## Run 217 stress 0.07374804 
    ## Run 218 stress 0.07547299 
    ## Run 219 stress 0.07374806 
    ## Run 220 stress 0.02693832 
    ## ... Procrustes: rmse 8.642641e-05  max resid 0.0001484711 
    ## ... Similar to previous best
    ## Run 221 stress 0.179905 
    ## Run 222 stress 0.1862385 
    ## Run 223 stress 0.07374802 
    ## Run 224 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001072948  max resid 0.0001812642 
    ## ... Similar to previous best
    ## Run 225 stress 0.07547299 
    ## Run 226 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001506506  max resid 0.0002585502 
    ## ... Similar to previous best
    ## Run 227 stress 0.07547299 
    ## Run 228 stress 0.02693831 
    ## ... Procrustes: rmse 7.253758e-05  max resid 0.0001217185 
    ## ... Similar to previous best
    ## Run 229 stress 0.02693831 
    ## ... Procrustes: rmse 9.915532e-05  max resid 0.0001687124 
    ## ... Similar to previous best
    ## Run 230 stress 0.0269383 
    ## ... Procrustes: rmse 8.617433e-05  max resid 0.0001473493 
    ## ... Similar to previous best
    ## Run 231 stress 0.07374805 
    ## Run 232 stress 0.07374801 
    ## Run 233 stress 0.0269383 
    ## ... Procrustes: rmse 6.127906e-05  max resid 0.0001049245 
    ## ... Similar to previous best
    ## Run 234 stress 0.1889142 
    ## Run 235 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001059132  max resid 0.0001822963 
    ## ... Similar to previous best
    ## Run 236 stress 0.02693832 
    ## ... Procrustes: rmse 0.000107655  max resid 0.0001834429 
    ## ... Similar to previous best
    ## Run 237 stress 0.07374811 
    ## Run 238 stress 0.075473 
    ## Run 239 stress 0.1889141 
    ## Run 240 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001575552  max resid 0.0002694128 
    ## ... Similar to previous best
    ## Run 241 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001188139  max resid 0.0002034959 
    ## ... Similar to previous best
    ## Run 242 stress 0.179905 
    ## Run 243 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001046424  max resid 0.0001793106 
    ## ... Similar to previous best
    ## Run 244 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001787527  max resid 0.0003058768 
    ## ... Similar to previous best
    ## Run 245 stress 0.02693843 
    ## ... Procrustes: rmse 0.000212978  max resid 0.0003631606 
    ## ... Similar to previous best
    ## Run 246 stress 0.075473 
    ## Run 247 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001087021  max resid 0.0001769071 
    ## ... Similar to previous best
    ## Run 248 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001133531  max resid 0.000193679 
    ## ... Similar to previous best
    ## Run 249 stress 0.02693838 
    ## ... Procrustes: rmse 9.938548e-05  max resid 0.0001660671 
    ## ... Similar to previous best
    ## Run 250 stress 0.1889142 
    ## Run 251 stress 0.07374807 
    ## Run 252 stress 0.07374802 
    ## Run 253 stress 0.0269383 
    ## ... Procrustes: rmse 7.618975e-05  max resid 0.0001293489 
    ## ... Similar to previous best
    ## Run 254 stress 0.02693828 
    ## ... New best solution
    ## ... Procrustes: rmse 1.002996e-06  max resid 1.435918e-06 
    ## ... Similar to previous best
    ## Run 255 stress 0.07374803 
    ## Run 256 stress 0.02693829 
    ## ... Procrustes: rmse 3.096278e-05  max resid 5.293232e-05 
    ## ... Similar to previous best
    ## Run 257 stress 0.1862389 
    ## Run 258 stress 0.179905 
    ## Run 259 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001460929  max resid 0.0002509217 
    ## ... Similar to previous best
    ## Run 260 stress 0.07547301 
    ## Run 261 stress 0.0737481 
    ## Run 262 stress 0.08382846 
    ## Run 263 stress 0.0939189 
    ## Run 264 stress 0.07547299 
    ## Run 265 stress 0.08382878 
    ## Run 266 stress 0.07374801 
    ## Run 267 stress 0.02693833 
    ## ... Procrustes: rmse 0.000126239  max resid 0.0002161717 
    ## ... Similar to previous best
    ## Run 268 stress 0.02693832 
    ## ... Procrustes: rmse 6.96441e-05  max resid 0.0001175871 
    ## ... Similar to previous best
    ## Run 269 stress 0.1756177 
    ## Run 270 stress 0.02693831 
    ## ... Procrustes: rmse 9.329816e-05  max resid 0.0001594971 
    ## ... Similar to previous best
    ## Run 271 stress 0.09391898 
    ## Run 272 stress 0.0269383 
    ## ... Procrustes: rmse 7.811862e-05  max resid 0.0001335008 
    ## ... Similar to previous best
    ## Run 273 stress 0.02693831 
    ## ... Procrustes: rmse 8.215763e-05  max resid 0.0001406574 
    ## ... Similar to previous best
    ## Run 274 stress 0.07374804 
    ## Run 275 stress 0.07374802 
    ## Run 276 stress 0.07374808 
    ## Run 277 stress 0.08382824 
    ## Run 278 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001047941  max resid 0.0001797905 
    ## ... Similar to previous best
    ## Run 279 stress 0.073748 
    ## Run 280 stress 0.1889143 
    ## Run 281 stress 0.02693831 
    ## ... Procrustes: rmse 7.615806e-05  max resid 0.0001306092 
    ## ... Similar to previous best
    ## Run 282 stress 0.07374801 
    ## Run 283 stress 0.07374802 
    ## Run 284 stress 0.07374801 
    ## Run 285 stress 0.07374806 
    ## Run 286 stress 0.07374813 
    ## Run 287 stress 0.08382843 
    ## Run 288 stress 0.0269383 
    ## ... Procrustes: rmse 5.635521e-05  max resid 9.525767e-05 
    ## ... Similar to previous best
    ## Run 289 stress 0.02693829 
    ## ... Procrustes: rmse 1.884112e-05  max resid 3.235417e-05 
    ## ... Similar to previous best
    ## Run 290 stress 0.02693831 
    ## ... Procrustes: rmse 0.0001057458  max resid 0.0001804992 
    ## ... Similar to previous best
    ## Run 291 stress 0.07547302 
    ## Run 292 stress 0.07547299 
    ## Run 293 stress 0.07374801 
    ## Run 294 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001656996  max resid 0.0002840017 
    ## ... Similar to previous best
    ## Run 295 stress 0.07547299 
    ## Run 296 stress 0.07374801 
    ## Run 297 stress 0.07374805 
    ## Run 298 stress 0.1862386 
    ## Run 299 stress 0.1889141 
    ## Run 300 stress 0.08382824 
    ## Run 301 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001333681  max resid 0.0002298763 
    ## ... Similar to previous best
    ## Run 302 stress 0.2273964 
    ## Run 303 stress 0.1862386 
    ## Run 304 stress 0.07547301 
    ## Run 305 stress 0.3083098 
    ## Run 306 stress 0.07374803 
    ## Run 307 stress 0.08382884 
    ## Run 308 stress 0.07547299 
    ## Run 309 stress 0.07374802 
    ## Run 310 stress 0.07374808 
    ## Run 311 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001357031  max resid 0.0002315102 
    ## ... Similar to previous best
    ## Run 312 stress 0.02693832 
    ## ... Procrustes: rmse 9.406788e-05  max resid 0.0001611946 
    ## ... Similar to previous best
    ## Run 313 stress 0.1862385 
    ## Run 314 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001373787  max resid 0.0002352356 
    ## ... Similar to previous best
    ## Run 315 stress 0.02693831 
    ## ... Procrustes: rmse 7.268595e-05  max resid 0.0001247787 
    ## ... Similar to previous best
    ## Run 316 stress 0.02693829 
    ## ... Procrustes: rmse 2.527997e-05  max resid 4.065809e-05 
    ## ... Similar to previous best
    ## Run 317 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001467687  max resid 0.0002510024 
    ## ... Similar to previous best
    ## Run 318 stress 0.179905 
    ## Run 319 stress 0.07547299 
    ## Run 320 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001021967  max resid 0.000175076 
    ## ... Similar to previous best
    ## Run 321 stress 0.07374802 
    ## Run 322 stress 0.07374809 
    ## Run 323 stress 0.07374807 
    ## Run 324 stress 0.07374803 
    ## Run 325 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001153254  max resid 0.0002000983 
    ## ... Similar to previous best
    ## Run 326 stress 0.08382843 
    ## Run 327 stress 0.07374803 
    ## Run 328 stress 0.179905 
    ## Run 329 stress 0.08382855 
    ## Run 330 stress 0.1799051 
    ## Run 331 stress 0.1889141 
    ## Run 332 stress 0.07374801 
    ## Run 333 stress 0.02693829 
    ## ... Procrustes: rmse 3.782072e-05  max resid 6.622266e-05 
    ## ... Similar to previous best
    ## Run 334 stress 0.1862386 
    ## Run 335 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001439352  max resid 0.0002474401 
    ## ... Similar to previous best
    ## Run 336 stress 0.08382847 
    ## Run 337 stress 0.02693831 
    ## ... Procrustes: rmse 9.280274e-05  max resid 0.0001592682 
    ## ... Similar to previous best
    ## Run 338 stress 0.07374801 
    ## Run 339 stress 0.07374803 
    ## Run 340 stress 0.1862387 
    ## Run 341 stress 0.07547299 
    ## Run 342 stress 0.02693829 
    ## ... Procrustes: rmse 3.68405e-06  max resid 6.365592e-06 
    ## ... Similar to previous best
    ## Run 343 stress 0.02693829 
    ## ... Procrustes: rmse 4.081489e-05  max resid 6.839339e-05 
    ## ... Similar to previous best
    ## Run 344 stress 0.09391908 
    ## Run 345 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001753524  max resid 0.000307414 
    ## ... Similar to previous best
    ## Run 346 stress 0.02693833 
    ## ... Procrustes: rmse 9.681256e-05  max resid 0.000165562 
    ## ... Similar to previous best
    ## Run 347 stress 0.07374802 
    ## Run 348 stress 0.08382831 
    ## Run 349 stress 0.07374801 
    ## Run 350 stress 0.07374808 
    ## Run 351 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001532442  max resid 0.0002618025 
    ## ... Similar to previous best
    ## Run 352 stress 0.07374805 
    ## Run 353 stress 0.1900871 
    ## Run 354 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001903938  max resid 0.0003257588 
    ## ... Similar to previous best
    ## Run 355 stress 0.07374811 
    ## Run 356 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001505863  max resid 0.0002625346 
    ## ... Similar to previous best
    ## Run 357 stress 0.1889141 
    ## Run 358 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001493852  max resid 0.0002559012 
    ## ... Similar to previous best
    ## Run 359 stress 0.0269383 
    ## ... Procrustes: rmse 6.152894e-05  max resid 0.0001049124 
    ## ... Similar to previous best
    ## Run 360 stress 0.07374808 
    ## Run 361 stress 0.1756177 
    ## Run 362 stress 0.07374812 
    ## Run 363 stress 0.07374813 
    ## Run 364 stress 0.07374802 
    ## Run 365 stress 0.0269383 
    ## ... Procrustes: rmse 5.293327e-05  max resid 8.912237e-05 
    ## ... Similar to previous best
    ## Run 366 stress 0.1756177 
    ## Run 367 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001683152  max resid 0.0002880446 
    ## ... Similar to previous best
    ## Run 368 stress 0.02693829 
    ## ... Procrustes: rmse 4.307978e-05  max resid 7.419197e-05 
    ## ... Similar to previous best
    ## Run 369 stress 0.075473 
    ## Run 370 stress 0.07374803 
    ## Run 371 stress 0.02693831 
    ## ... Procrustes: rmse 9.079025e-05  max resid 0.000153794 
    ## ... Similar to previous best
    ## Run 372 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001382854  max resid 0.0002365346 
    ## ... Similar to previous best
    ## Run 373 stress 0.3083098 
    ## Run 374 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001618608  max resid 0.0002770803 
    ## ... Similar to previous best
    ## Run 375 stress 0.08382833 
    ## Run 376 stress 0.179905 
    ## Run 377 stress 0.07547299 
    ## Run 378 stress 0.07547299 
    ## Run 379 stress 0.07374806 
    ## Run 380 stress 0.07374801 
    ## Run 381 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001147238  max resid 0.0001961487 
    ## ... Similar to previous best
    ## Run 382 stress 0.07547303 
    ## Run 383 stress 0.1756177 
    ## Run 384 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001162574  max resid 0.0001986836 
    ## ... Similar to previous best
    ## Run 385 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001345515  max resid 0.0002294498 
    ## ... Similar to previous best
    ## Run 386 stress 0.0939192 
    ## Run 387 stress 0.0737481 
    ## Run 388 stress 0.07374802 
    ## Run 389 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001359727  max resid 0.0002321063 
    ## ... Similar to previous best
    ## Run 390 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001425464  max resid 0.0002450584 
    ## ... Similar to previous best
    ## Run 391 stress 0.07547299 
    ## Run 392 stress 0.07374803 
    ## Run 393 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001244271  max resid 0.0002131071 
    ## ... Similar to previous best
    ## Run 394 stress 0.1756177 
    ## Run 395 stress 0.1756177 
    ## Run 396 stress 0.179905 
    ## Run 397 stress 0.07374806 
    ## Run 398 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001137213  max resid 0.0001945694 
    ## ... Similar to previous best
    ## Run 399 stress 0.07374806 
    ## Run 400 stress 0.0737481 
    ## Run 401 stress 0.02693835 
    ## ... Procrustes: rmse 0.000143376  max resid 0.0002444453 
    ## ... Similar to previous best
    ## Run 402 stress 0.07374802 
    ## Run 403 stress 0.0838285 
    ## Run 404 stress 0.073748 
    ## Run 405 stress 0.07547302 
    ## Run 406 stress 0.07374802 
    ## Run 407 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001482589  max resid 0.0002520443 
    ## ... Similar to previous best
    ## Run 408 stress 0.0269383 
    ## ... Procrustes: rmse 6.379485e-05  max resid 0.0001094303 
    ## ... Similar to previous best
    ## Run 409 stress 0.073748 
    ## Run 410 stress 0.2852148 
    ## Run 411 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001670809  max resid 0.0002859369 
    ## ... Similar to previous best
    ## Run 412 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001520449  max resid 0.0002589274 
    ## ... Similar to previous best
    ## Run 413 stress 0.02693829 
    ## ... Procrustes: rmse 3.98797e-05  max resid 6.850016e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.07547299 
    ## Run 415 stress 0.07374804 
    ## Run 416 stress 0.08382854 
    ## Run 417 stress 0.02693831 
    ## ... Procrustes: rmse 9.802176e-05  max resid 0.0001672263 
    ## ... Similar to previous best
    ## Run 418 stress 0.02693831 
    ## ... Procrustes: rmse 0.0001040034  max resid 0.0001784103 
    ## ... Similar to previous best
    ## Run 419 stress 0.07374803 
    ## Run 420 stress 0.1862392 
    ## Run 421 stress 0.02693833 
    ## ... Procrustes: rmse 9.956616e-05  max resid 0.0001702863 
    ## ... Similar to previous best
    ## Run 422 stress 0.07374802 
    ## Run 423 stress 0.02693829 
    ## ... Procrustes: rmse 3.299243e-05  max resid 5.636835e-05 
    ## ... Similar to previous best
    ## Run 424 stress 0.07547299 
    ## Run 425 stress 0.07547299 
    ## Run 426 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001636141  max resid 0.0002803382 
    ## ... Similar to previous best
    ## Run 427 stress 0.07374814 
    ## Run 428 stress 0.1862389 
    ## Run 429 stress 0.02693831 
    ## ... Procrustes: rmse 4.27191e-05  max resid 7.005612e-05 
    ## ... Similar to previous best
    ## Run 430 stress 0.02693829 
    ## ... Procrustes: rmse 1.56872e-05  max resid 2.525018e-05 
    ## ... Similar to previous best
    ## Run 431 stress 0.1756177 
    ## Run 432 stress 0.075473 
    ## Run 433 stress 0.07374803 
    ## Run 434 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001990202  max resid 0.0003416077 
    ## ... Similar to previous best
    ## Run 435 stress 0.02693829 
    ## ... Procrustes: rmse 5.69836e-05  max resid 9.73847e-05 
    ## ... Similar to previous best
    ## Run 436 stress 0.09391926 
    ## Run 437 stress 0.07374804 
    ## Run 438 stress 0.02693841 
    ## ... Procrustes: rmse 0.0001764235  max resid 0.0003025946 
    ## ... Similar to previous best
    ## Run 439 stress 0.07374812 
    ## Run 440 stress 0.1799051 
    ## Run 441 stress 0.07374805 
    ## Run 442 stress 0.07374802 
    ## Run 443 stress 0.1889141 
    ## Run 444 stress 0.1889141 
    ## Run 445 stress 0.02693829 
    ## ... Procrustes: rmse 3.773878e-05  max resid 6.370892e-05 
    ## ... Similar to previous best
    ## Run 446 stress 0.08382848 
    ## Run 447 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001422751  max resid 0.0002422961 
    ## ... Similar to previous best
    ## Run 448 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001142251  max resid 0.000194188 
    ## ... Similar to previous best
    ## Run 449 stress 0.07374805 
    ## Run 450 stress 0.07374807 
    ## Run 451 stress 0.02693831 
    ## ... Procrustes: rmse 0.0001057495  max resid 0.0001811147 
    ## ... Similar to previous best
    ## Run 452 stress 0.0269383 
    ## ... Procrustes: rmse 7.192549e-05  max resid 0.0001232015 
    ## ... Similar to previous best
    ## Run 453 stress 0.179905 
    ## Run 454 stress 0.179905 
    ## Run 455 stress 0.1862383 
    ## Run 456 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001280904  max resid 0.0002199794 
    ## ... Similar to previous best
    ## Run 457 stress 0.1889143 
    ## Run 458 stress 0.073748 
    ## Run 459 stress 0.07374802 
    ## Run 460 stress 0.1862386 
    ## Run 461 stress 0.02693834 
    ## ... Procrustes: rmse 0.000136998  max resid 0.0002346169 
    ## ... Similar to previous best
    ## Run 462 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001425072  max resid 0.0002446026 
    ## ... Similar to previous best
    ## Run 463 stress 0.02693831 
    ## ... Procrustes: rmse 6.19826e-05  max resid 0.0001071936 
    ## ... Similar to previous best
    ## Run 464 stress 0.285217 
    ## Run 465 stress 0.08382838 
    ## Run 466 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001776989  max resid 0.0003033759 
    ## ... Similar to previous best
    ## Run 467 stress 0.02693831 
    ## ... Procrustes: rmse 6.947739e-05  max resid 0.0001186611 
    ## ... Similar to previous best
    ## Run 468 stress 0.179905 
    ## Run 469 stress 0.02693829 
    ## ... Procrustes: rmse 4.152335e-05  max resid 7.073323e-05 
    ## ... Similar to previous best
    ## Run 470 stress 0.08382842 
    ## Run 471 stress 0.07374805 
    ## Run 472 stress 0.07547299 
    ## Run 473 stress 0.08382836 
    ## Run 474 stress 0.09391909 
    ## Run 475 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001076378  max resid 0.0001851219 
    ## ... Similar to previous best
    ## Run 476 stress 0.07547301 
    ## Run 477 stress 0.07374804 
    ## Run 478 stress 0.07374805 
    ## Run 479 stress 0.02693832 
    ## ... Procrustes: rmse 8.608638e-05  max resid 0.0001475994 
    ## ... Similar to previous best
    ## Run 480 stress 0.07374808 
    ## Run 481 stress 0.07547301 
    ## Run 482 stress 0.08382859 
    ## Run 483 stress 0.07374802 
    ## Run 484 stress 0.08382853 
    ## Run 485 stress 0.07547299 
    ## Run 486 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001001238  max resid 0.0001719632 
    ## ... Similar to previous best
    ## Run 487 stress 0.1889141 
    ## Run 488 stress 0.07374812 
    ## Run 489 stress 0.1889143 
    ## Run 490 stress 0.179905 
    ## Run 491 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001967775  max resid 0.0003361801 
    ## ... Similar to previous best
    ## Run 492 stress 0.07547299 
    ## Run 493 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001700122  max resid 0.0002896224 
    ## ... Similar to previous best
    ## Run 494 stress 0.1889141 
    ## Run 495 stress 0.0269383 
    ## ... Procrustes: rmse 8.370014e-05  max resid 0.000143033 
    ## ... Similar to previous best
    ## Run 496 stress 0.07374808 
    ## Run 497 stress 0.07374804 
    ## Run 498 stress 0.1756177 
    ## Run 499 stress 0.1889141 
    ## Run 500 stress 0.08382821 
    ## *** Best solution repeated 82 times

``` r
### Environmental
# Surveyed sites 
PD_beta_env_NMDS <- metaMDS(PD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06330868 
    ## Run 1 stress 0.09217199 
    ## Run 2 stress 0.06330871 
    ## ... Procrustes: rmse 7.667043e-05  max resid 0.0001504282 
    ## ... Similar to previous best
    ## Run 3 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 1.70334e-05  max resid 3.526278e-05 
    ## ... Similar to previous best
    ## Run 4 stress 0.06402616 
    ## Run 5 stress 0.06402639 
    ## Run 6 stress 0.06402623 
    ## Run 7 stress 0.06463313 
    ## Run 8 stress 0.06463314 
    ## Run 9 stress 0.06463313 
    ## Run 10 stress 0.0640263 
    ## Run 11 stress 0.0640263 
    ## Run 12 stress 0.06330868 
    ## ... Procrustes: rmse 1.203864e-05  max resid 2.194453e-05 
    ## ... Similar to previous best
    ## Run 13 stress 0.3699075 
    ## Run 14 stress 0.06330871 
    ## ... Procrustes: rmse 5.588322e-05  max resid 9.903335e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 9.30119e-06  max resid 3.404455e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.1004383 
    ## Run 17 stress 0.06402623 
    ## Run 18 stress 0.06330871 
    ## ... Procrustes: rmse 7.650807e-05  max resid 0.0002078442 
    ## ... Similar to previous best
    ## Run 19 stress 0.06463313 
    ## Run 20 stress 0.06402631 
    ## Run 21 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 8.639698e-07  max resid 1.94083e-06 
    ## ... Similar to previous best
    ## Run 22 stress 0.08357057 
    ## Run 23 stress 0.06330869 
    ## ... Procrustes: rmse 2.518162e-05  max resid 5.995167e-05 
    ## ... Similar to previous best
    ## Run 24 stress 0.06402623 
    ## Run 25 stress 0.0923208 
    ## Run 26 stress 0.0823374 
    ## Run 27 stress 0.0640263 
    ## Run 28 stress 0.06330872 
    ## ... Procrustes: rmse 8.906948e-05  max resid 0.000195983 
    ## ... Similar to previous best
    ## Run 29 stress 0.06402616 
    ## Run 30 stress 0.08431279 
    ## Run 31 stress 0.08301942 
    ## Run 32 stress 0.06463313 
    ## Run 33 stress 0.06330869 
    ## ... Procrustes: rmse 3.174493e-05  max resid 0.0001004939 
    ## ... Similar to previous best
    ## Run 34 stress 0.08126181 
    ## Run 35 stress 0.06463313 
    ## Run 36 stress 0.06402638 
    ## Run 37 stress 0.06463313 
    ## Run 38 stress 0.06402615 
    ## Run 39 stress 0.06463313 
    ## Run 40 stress 0.06463316 
    ## Run 41 stress 0.06463314 
    ## Run 42 stress 0.06402631 
    ## Run 43 stress 0.06330869 
    ## ... Procrustes: rmse 2.363122e-05  max resid 5.38845e-05 
    ## ... Similar to previous best
    ## Run 44 stress 0.06402631 
    ## Run 45 stress 0.06330875 
    ## ... Procrustes: rmse 0.0001050894  max resid 0.000236276 
    ## ... Similar to previous best
    ## Run 46 stress 0.06330869 
    ## ... Procrustes: rmse 3.005664e-05  max resid 9.889193e-05 
    ## ... Similar to previous best
    ## Run 47 stress 0.06330868 
    ## ... Procrustes: rmse 2.707571e-06  max resid 4.653775e-06 
    ## ... Similar to previous best
    ## Run 48 stress 0.0633087 
    ## ... Procrustes: rmse 5.432163e-05  max resid 0.0001207492 
    ## ... Similar to previous best
    ## Run 49 stress 0.06402616 
    ## Run 50 stress 0.06402624 
    ## Run 51 stress 0.06463313 
    ## Run 52 stress 0.08126217 
    ## Run 53 stress 0.06402621 
    ## Run 54 stress 0.06463313 
    ## Run 55 stress 0.081262 
    ## Run 56 stress 0.08233728 
    ## Run 57 stress 0.08233747 
    ## Run 58 stress 0.3744265 
    ## Run 59 stress 0.06330868 
    ## ... Procrustes: rmse 1.859959e-05  max resid 4.255147e-05 
    ## ... Similar to previous best
    ## Run 60 stress 0.06330868 
    ## ... Procrustes: rmse 6.154168e-06  max resid 1.677465e-05 
    ## ... Similar to previous best
    ## Run 61 stress 0.06330871 
    ## ... Procrustes: rmse 6.914129e-05  max resid 0.0001538815 
    ## ... Similar to previous best
    ## Run 62 stress 0.06330869 
    ## ... Procrustes: rmse 1.822307e-05  max resid 4.472381e-05 
    ## ... Similar to previous best
    ## Run 63 stress 0.06330869 
    ## ... Procrustes: rmse 4.7626e-05  max resid 0.000106149 
    ## ... Similar to previous best
    ## Run 64 stress 0.06463313 
    ## Run 65 stress 0.0640263 
    ## Run 66 stress 0.06330875 
    ## ... Procrustes: rmse 8.068132e-05  max resid 0.0001759089 
    ## ... Similar to previous best
    ## Run 67 stress 0.06463317 
    ## Run 68 stress 0.06330869 
    ## ... Procrustes: rmse 3.005576e-05  max resid 0.0001025648 
    ## ... Similar to previous best
    ## Run 69 stress 0.06330869 
    ## ... Procrustes: rmse 3.424247e-05  max resid 7.596076e-05 
    ## ... Similar to previous best
    ## Run 70 stress 0.08126202 
    ## Run 71 stress 0.06402625 
    ## Run 72 stress 0.0633087 
    ## ... Procrustes: rmse 5.599232e-05  max resid 9.726777e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.0808905 
    ## Run 74 stress 0.06463316 
    ## Run 75 stress 0.09270079 
    ## Run 76 stress 0.06463314 
    ## Run 77 stress 0.08403369 
    ## Run 78 stress 0.06402622 
    ## Run 79 stress 0.06330872 
    ## ... Procrustes: rmse 8.037047e-05  max resid 0.00017194 
    ## ... Similar to previous best
    ## Run 80 stress 0.06330868 
    ## ... Procrustes: rmse 1.199413e-05  max resid 2.66753e-05 
    ## ... Similar to previous best
    ## Run 81 stress 0.09106921 
    ## Run 82 stress 0.08089066 
    ## Run 83 stress 0.06463314 
    ## Run 84 stress 0.06402629 
    ## Run 85 stress 0.06402618 
    ## Run 86 stress 0.06463315 
    ## Run 87 stress 0.06463313 
    ## Run 88 stress 0.06330869 
    ## ... Procrustes: rmse 2.485814e-05  max resid 5.375808e-05 
    ## ... Similar to previous best
    ## Run 89 stress 0.06402617 
    ## Run 90 stress 0.0640262 
    ## Run 91 stress 0.06330869 
    ## ... Procrustes: rmse 2.470516e-05  max resid 5.998001e-05 
    ## ... Similar to previous best
    ## Run 92 stress 0.08089032 
    ## Run 93 stress 0.06402622 
    ## Run 94 stress 0.06463314 
    ## Run 95 stress 0.08089068 
    ## Run 96 stress 0.06463314 
    ## Run 97 stress 0.06402628 
    ## Run 98 stress 0.0937868 
    ## Run 99 stress 0.06330869 
    ## ... Procrustes: rmse 3.524153e-05  max resid 6.272712e-05 
    ## ... Similar to previous best
    ## Run 100 stress 0.06463313 
    ## Run 101 stress 0.06330868 
    ## ... Procrustes: rmse 7.640124e-06  max resid 1.545447e-05 
    ## ... Similar to previous best
    ## Run 102 stress 0.06330871 
    ## ... Procrustes: rmse 7.35889e-05  max resid 0.0001641381 
    ## ... Similar to previous best
    ## Run 103 stress 0.06402631 
    ## Run 104 stress 0.06402632 
    ## Run 105 stress 0.0633087 
    ## ... Procrustes: rmse 4.269382e-05  max resid 9.512067e-05 
    ## ... Similar to previous best
    ## Run 106 stress 0.06330868 
    ## ... Procrustes: rmse 2.444419e-05  max resid 9.047838e-05 
    ## ... Similar to previous best
    ## Run 107 stress 0.09788388 
    ## Run 108 stress 0.06330869 
    ## ... Procrustes: rmse 3.619417e-05  max resid 8.154198e-05 
    ## ... Similar to previous best
    ## Run 109 stress 0.06402634 
    ## Run 110 stress 0.06463313 
    ## Run 111 stress 0.06463316 
    ## Run 112 stress 0.06463315 
    ## Run 113 stress 0.06402623 
    ## Run 114 stress 0.06330872 
    ## ... Procrustes: rmse 6.20005e-05  max resid 0.0001346344 
    ## ... Similar to previous best
    ## Run 115 stress 0.06463313 
    ## Run 116 stress 0.06402625 
    ## Run 117 stress 0.06330868 
    ## ... Procrustes: rmse 9.371793e-06  max resid 1.976532e-05 
    ## ... Similar to previous best
    ## Run 118 stress 0.06463313 
    ## Run 119 stress 0.06463318 
    ## Run 120 stress 0.06330873 
    ## ... Procrustes: rmse 8.041482e-05  max resid 0.0001761185 
    ## ... Similar to previous best
    ## Run 121 stress 0.06330869 
    ## ... Procrustes: rmse 4.97031e-05  max resid 0.0001073163 
    ## ... Similar to previous best
    ## Run 122 stress 0.06463313 
    ## Run 123 stress 0.06402616 
    ## Run 124 stress 0.06330868 
    ## ... Procrustes: rmse 1.942644e-05  max resid 5.54115e-05 
    ## ... Similar to previous best
    ## Run 125 stress 0.06402615 
    ## Run 126 stress 0.06330872 
    ## ... Procrustes: rmse 7.578793e-05  max resid 0.0001695872 
    ## ... Similar to previous best
    ## Run 127 stress 0.08357059 
    ## Run 128 stress 0.06402622 
    ## Run 129 stress 0.06402636 
    ## Run 130 stress 0.09245944 
    ## Run 131 stress 0.06463315 
    ## Run 132 stress 0.06463315 
    ## Run 133 stress 0.06402621 
    ## Run 134 stress 0.06402622 
    ## Run 135 stress 0.08089064 
    ## Run 136 stress 0.06330868 
    ## ... Procrustes: rmse 3.906631e-06  max resid 1.415495e-05 
    ## ... Similar to previous best
    ## Run 137 stress 0.0633087 
    ## ... Procrustes: rmse 4.879697e-05  max resid 9.246487e-05 
    ## ... Similar to previous best
    ## Run 138 stress 0.06330868 
    ## ... Procrustes: rmse 1.291319e-05  max resid 2.85656e-05 
    ## ... Similar to previous best
    ## Run 139 stress 0.06402626 
    ## Run 140 stress 0.06330878 
    ## ... Procrustes: rmse 0.0001162275  max resid 0.0002345067 
    ## ... Similar to previous best
    ## Run 141 stress 0.06402616 
    ## Run 142 stress 0.0843124 
    ## Run 143 stress 0.08233734 
    ## Run 144 stress 0.06330869 
    ## ... Procrustes: rmse 1.728366e-05  max resid 5.639341e-05 
    ## ... Similar to previous best
    ## Run 145 stress 0.06463313 
    ## Run 146 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 2.727992e-06  max resid 5.513651e-06 
    ## ... Similar to previous best
    ## Run 147 stress 0.06330872 
    ## ... Procrustes: rmse 7.339936e-05  max resid 0.0001645279 
    ## ... Similar to previous best
    ## Run 148 stress 0.06330869 
    ## ... Procrustes: rmse 3.749539e-05  max resid 7.77591e-05 
    ## ... Similar to previous best
    ## Run 149 stress 0.09429108 
    ## Run 150 stress 0.08315616 
    ## Run 151 stress 0.06330868 
    ## ... Procrustes: rmse 5.813923e-06  max resid 1.840429e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.06402616 
    ## Run 153 stress 0.06330869 
    ## ... Procrustes: rmse 2.711225e-05  max resid 6.13075e-05 
    ## ... Similar to previous best
    ## Run 154 stress 0.06463313 
    ## Run 155 stress 0.06402618 
    ## Run 156 stress 0.06402617 
    ## Run 157 stress 0.09106935 
    ## Run 158 stress 0.06330869 
    ## ... Procrustes: rmse 4.573256e-05  max resid 0.0001022536 
    ## ... Similar to previous best
    ## Run 159 stress 0.08126215 
    ## Run 160 stress 0.06330869 
    ## ... Procrustes: rmse 2.553544e-05  max resid 5.653883e-05 
    ## ... Similar to previous best
    ## Run 161 stress 0.06463313 
    ## Run 162 stress 0.06402616 
    ## Run 163 stress 0.06463316 
    ## Run 164 stress 0.06463313 
    ## Run 165 stress 0.06402619 
    ## Run 166 stress 0.09533682 
    ## Run 167 stress 0.06402616 
    ## Run 168 stress 0.06463313 
    ## Run 169 stress 0.06402632 
    ## Run 170 stress 0.06402631 
    ## Run 171 stress 0.06330872 
    ## ... Procrustes: rmse 2.034741e-05  max resid 4.395917e-05 
    ## ... Similar to previous best
    ## Run 172 stress 0.06463317 
    ## Run 173 stress 0.06330868 
    ## ... Procrustes: rmse 7.729751e-06  max resid 2.766777e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.08126208 
    ## Run 175 stress 0.06463313 
    ## Run 176 stress 0.06463315 
    ## Run 177 stress 0.06402622 
    ## Run 178 stress 0.06463314 
    ## Run 179 stress 0.06463313 
    ## Run 180 stress 0.06402617 
    ## Run 181 stress 0.06330868 
    ## ... Procrustes: rmse 5.545951e-06  max resid 1.729824e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.06463314 
    ## Run 183 stress 0.06402626 
    ## Run 184 stress 0.06402629 
    ## Run 185 stress 0.08126168 
    ## Run 186 stress 0.09336679 
    ## Run 187 stress 0.06463314 
    ## Run 188 stress 0.06463315 
    ## Run 189 stress 0.09270092 
    ## Run 190 stress 0.06330874 
    ## ... Procrustes: rmse 9.038243e-05  max resid 0.0002003439 
    ## ... Similar to previous best
    ## Run 191 stress 0.06402617 
    ## Run 192 stress 0.08089068 
    ## Run 193 stress 0.06402623 
    ## Run 194 stress 0.06463314 
    ## Run 195 stress 0.06463315 
    ## Run 196 stress 0.06330869 
    ## ... Procrustes: rmse 4.610117e-05  max resid 0.0001665344 
    ## ... Similar to previous best
    ## Run 197 stress 0.06402617 
    ## Run 198 stress 0.06463315 
    ## Run 199 stress 0.06330868 
    ## ... Procrustes: rmse 1.092736e-05  max resid 2.681268e-05 
    ## ... Similar to previous best
    ## Run 200 stress 0.06463314 
    ## Run 201 stress 0.0993032 
    ## Run 202 stress 0.06402619 
    ## Run 203 stress 0.06463313 
    ## Run 204 stress 0.06330868 
    ## ... Procrustes: rmse 6.762834e-06  max resid 1.432398e-05 
    ## ... Similar to previous best
    ## Run 205 stress 0.06463315 
    ## Run 206 stress 0.06402622 
    ## Run 207 stress 0.06463315 
    ## Run 208 stress 0.09504611 
    ## Run 209 stress 0.06402621 
    ## Run 210 stress 0.06330869 
    ## ... Procrustes: rmse 3.762046e-05  max resid 8.118113e-05 
    ## ... Similar to previous best
    ## Run 211 stress 0.0640262 
    ## Run 212 stress 0.06402618 
    ## Run 213 stress 0.08315607 
    ## Run 214 stress 0.06330871 
    ## ... Procrustes: rmse 6.450611e-05  max resid 0.0001131798 
    ## ... Similar to previous best
    ## Run 215 stress 0.06463317 
    ## Run 216 stress 0.06463313 
    ## Run 217 stress 0.08315606 
    ## Run 218 stress 0.06330869 
    ## ... Procrustes: rmse 3.848926e-05  max resid 8.67698e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.06402627 
    ## Run 220 stress 0.06402617 
    ## Run 221 stress 0.06463313 
    ## Run 222 stress 0.08315621 
    ## Run 223 stress 0.06402617 
    ## Run 224 stress 0.06463314 
    ## Run 225 stress 0.06402616 
    ## Run 226 stress 0.06463313 
    ## Run 227 stress 0.06330869 
    ## ... Procrustes: rmse 3.663545e-05  max resid 8.167798e-05 
    ## ... Similar to previous best
    ## Run 228 stress 0.06402619 
    ## Run 229 stress 0.09246065 
    ## Run 230 stress 0.0633087 
    ## ... Procrustes: rmse 6.047285e-05  max resid 0.0001347205 
    ## ... Similar to previous best
    ## Run 231 stress 0.06330868 
    ## ... Procrustes: rmse 1.394065e-05  max resid 5.086839e-05 
    ## ... Similar to previous best
    ## Run 232 stress 0.06463313 
    ## Run 233 stress 0.06330869 
    ## ... Procrustes: rmse 4.191505e-05  max resid 0.0001549685 
    ## ... Similar to previous best
    ## Run 234 stress 0.06330871 
    ## ... Procrustes: rmse 6.368251e-05  max resid 0.0001422544 
    ## ... Similar to previous best
    ## Run 235 stress 0.06330868 
    ## ... Procrustes: rmse 2.054229e-05  max resid 4.639501e-05 
    ## ... Similar to previous best
    ## Run 236 stress 0.08246342 
    ## Run 237 stress 0.06402616 
    ## Run 238 stress 0.06463314 
    ## Run 239 stress 0.06330877 
    ## ... Procrustes: rmse 0.0001214743  max resid 0.0002644284 
    ## ... Similar to previous best
    ## Run 240 stress 0.06463314 
    ## Run 241 stress 0.06330868 
    ## ... Procrustes: rmse 1.410449e-05  max resid 4.5488e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.06330869 
    ## ... Procrustes: rmse 4.560842e-05  max resid 0.0001006141 
    ## ... Similar to previous best
    ## Run 243 stress 0.06463316 
    ## Run 244 stress 0.06402618 
    ## Run 245 stress 0.06463315 
    ## Run 246 stress 0.06402616 
    ## Run 247 stress 0.0633087 
    ## ... Procrustes: rmse 2.042028e-05  max resid 4.635322e-05 
    ## ... Similar to previous best
    ## Run 248 stress 0.09336651 
    ## Run 249 stress 0.08089061 
    ## Run 250 stress 0.06463315 
    ## Run 251 stress 0.094638 
    ## Run 252 stress 0.06463314 
    ## Run 253 stress 0.06463313 
    ## Run 254 stress 0.06402624 
    ## Run 255 stress 0.0633087 
    ## ... Procrustes: rmse 4.699316e-05  max resid 0.000107239 
    ## ... Similar to previous best
    ## Run 256 stress 0.09999193 
    ## Run 257 stress 0.06402621 
    ## Run 258 stress 0.06463313 
    ## Run 259 stress 0.06330869 
    ## ... Procrustes: rmse 1.7959e-05  max resid 3.395826e-05 
    ## ... Similar to previous best
    ## Run 260 stress 0.06402624 
    ## Run 261 stress 0.09682895 
    ## Run 262 stress 0.06402623 
    ## Run 263 stress 0.09462051 
    ## Run 264 stress 0.06402619 
    ## Run 265 stress 0.06402616 
    ## Run 266 stress 0.06463313 
    ## Run 267 stress 0.06463315 
    ## Run 268 stress 0.06463314 
    ## Run 269 stress 0.08089067 
    ## Run 270 stress 0.06402634 
    ## Run 271 stress 0.06402616 
    ## Run 272 stress 0.06402628 
    ## Run 273 stress 0.06330869 
    ## ... Procrustes: rmse 2.148993e-05  max resid 5.004125e-05 
    ## ... Similar to previous best
    ## Run 274 stress 0.09336685 
    ## Run 275 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 2.811134e-06  max resid 6.849376e-06 
    ## ... Similar to previous best
    ## Run 276 stress 0.06402619 
    ## Run 277 stress 0.06330874 
    ## ... Procrustes: rmse 9.566225e-05  max resid 0.0002156144 
    ## ... Similar to previous best
    ## Run 278 stress 0.08089007 
    ## Run 279 stress 0.08089035 
    ## Run 280 stress 0.06463314 
    ## Run 281 stress 0.06330873 
    ## ... Procrustes: rmse 8.399346e-05  max resid 0.0001427902 
    ## ... Similar to previous best
    ## Run 282 stress 0.06330868 
    ## ... Procrustes: rmse 1.263817e-05  max resid 2.770276e-05 
    ## ... Similar to previous best
    ## Run 283 stress 0.06402639 
    ## Run 284 stress 0.09157919 
    ## Run 285 stress 0.06463313 
    ## Run 286 stress 0.06330868 
    ## ... Procrustes: rmse 2.478427e-05  max resid 5.518103e-05 
    ## ... Similar to previous best
    ## Run 287 stress 0.06463313 
    ## Run 288 stress 0.06402617 
    ## Run 289 stress 0.0640264 
    ## Run 290 stress 0.08327884 
    ## Run 291 stress 0.06402641 
    ## Run 292 stress 0.06402621 
    ## Run 293 stress 0.06330869 
    ## ... Procrustes: rmse 1.782169e-05  max resid 6.043875e-05 
    ## ... Similar to previous best
    ## Run 294 stress 0.06463313 
    ## Run 295 stress 0.0923346 
    ## Run 296 stress 0.06330869 
    ## ... Procrustes: rmse 3.39325e-05  max resid 6.64163e-05 
    ## ... Similar to previous best
    ## Run 297 stress 0.06330869 
    ## ... Procrustes: rmse 2.933286e-05  max resid 6.420608e-05 
    ## ... Similar to previous best
    ## Run 298 stress 0.06402622 
    ## Run 299 stress 0.09217185 
    ## Run 300 stress 0.06402631 
    ## Run 301 stress 0.06402616 
    ## Run 302 stress 0.06463317 
    ## Run 303 stress 0.06463315 
    ## Run 304 stress 0.06463313 
    ## Run 305 stress 0.06402633 
    ## Run 306 stress 0.06463313 
    ## Run 307 stress 0.06463314 
    ## Run 308 stress 0.06402633 
    ## Run 309 stress 0.06402625 
    ## Run 310 stress 0.08390019 
    ## Run 311 stress 0.06402617 
    ## Run 312 stress 0.06463313 
    ## Run 313 stress 0.06463313 
    ## Run 314 stress 0.06463313 
    ## Run 315 stress 0.06463313 
    ## Run 316 stress 0.06330874 
    ## ... Procrustes: rmse 9.829636e-05  max resid 0.0002202271 
    ## ... Similar to previous best
    ## Run 317 stress 0.06463313 
    ## Run 318 stress 0.06330868 
    ## ... Procrustes: rmse 8.215183e-06  max resid 1.650965e-05 
    ## ... Similar to previous best
    ## Run 319 stress 0.09106994 
    ## Run 320 stress 0.09524012 
    ## Run 321 stress 0.06463318 
    ## Run 322 stress 0.06463314 
    ## Run 323 stress 0.06402621 
    ## Run 324 stress 0.09157947 
    ## Run 325 stress 0.06463315 
    ## Run 326 stress 0.06463313 
    ## Run 327 stress 0.06330869 
    ## ... Procrustes: rmse 3.765542e-05  max resid 8.340255e-05 
    ## ... Similar to previous best
    ## Run 328 stress 0.0640262 
    ## Run 329 stress 0.06463313 
    ## Run 330 stress 0.09551327 
    ## Run 331 stress 0.06330871 
    ## ... Procrustes: rmse 6.052986e-05  max resid 0.0001023622 
    ## ... Similar to previous best
    ## Run 332 stress 0.06402643 
    ## Run 333 stress 0.06402616 
    ## Run 334 stress 0.06402616 
    ## Run 335 stress 0.06463313 
    ## Run 336 stress 0.06402626 
    ## Run 337 stress 0.06330876 
    ## ... Procrustes: rmse 0.000103935  max resid 0.0002333009 
    ## ... Similar to previous best
    ## Run 338 stress 0.06330876 
    ## ... Procrustes: rmse 9.970737e-05  max resid 0.0002151349 
    ## ... Similar to previous best
    ## Run 339 stress 0.06402626 
    ## Run 340 stress 0.06402627 
    ## Run 341 stress 0.09369689 
    ## Run 342 stress 0.08315609 
    ## Run 343 stress 0.06402631 
    ## Run 344 stress 0.06402626 
    ## Run 345 stress 0.06402627 
    ## Run 346 stress 0.08403577 
    ## Run 347 stress 0.06463313 
    ## Run 348 stress 0.09523987 
    ## Run 349 stress 0.06330869 
    ## ... Procrustes: rmse 3.002588e-05  max resid 5.916996e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.08357055 
    ## Run 351 stress 0.0640262 
    ## Run 352 stress 0.06463313 
    ## Run 353 stress 0.09106914 
    ## Run 354 stress 0.06463313 
    ## Run 355 stress 0.0832797 
    ## Run 356 stress 0.06402635 
    ## Run 357 stress 0.06463313 
    ## Run 358 stress 0.06330874 
    ## ... Procrustes: rmse 9.816945e-05  max resid 0.0002211673 
    ## ... Similar to previous best
    ## Run 359 stress 0.06402623 
    ## Run 360 stress 0.06402616 
    ## Run 361 stress 0.06330871 
    ## ... Procrustes: rmse 6.664608e-05  max resid 0.0001463184 
    ## ... Similar to previous best
    ## Run 362 stress 0.06463317 
    ## Run 363 stress 0.0633087 
    ## ... Procrustes: rmse 3.868769e-05  max resid 8.215541e-05 
    ## ... Similar to previous best
    ## Run 364 stress 0.09824418 
    ## Run 365 stress 0.06402615 
    ## Run 366 stress 0.06402638 
    ## Run 367 stress 0.08089012 
    ## Run 368 stress 0.06402619 
    ## Run 369 stress 0.08605812 
    ## Run 370 stress 0.06463313 
    ## Run 371 stress 0.06463314 
    ## Run 372 stress 0.0640263 
    ## Run 373 stress 0.06463313 
    ## Run 374 stress 0.06463313 
    ## Run 375 stress 0.08351229 
    ## Run 376 stress 0.06330879 
    ## ... Procrustes: rmse 0.0001347514  max resid 0.0003005335 
    ## ... Similar to previous best
    ## Run 377 stress 0.09504616 
    ## Run 378 stress 0.06330871 
    ## ... Procrustes: rmse 5.88028e-05  max resid 0.0001332536 
    ## ... Similar to previous best
    ## Run 379 stress 0.06463315 
    ## Run 380 stress 0.06463313 
    ## Run 381 stress 0.06463314 
    ## Run 382 stress 0.06402615 
    ## Run 383 stress 0.08357052 
    ## Run 384 stress 0.06330869 
    ## ... Procrustes: rmse 3.238516e-05  max resid 9.939006e-05 
    ## ... Similar to previous best
    ## Run 385 stress 0.06402626 
    ## Run 386 stress 0.06463313 
    ## Run 387 stress 0.06402616 
    ## Run 388 stress 0.06402622 
    ## Run 389 stress 0.06402638 
    ## Run 390 stress 0.06463313 
    ## Run 391 stress 0.06463313 
    ## Run 392 stress 0.06330873 
    ## ... Procrustes: rmse 8.645268e-05  max resid 0.0001902117 
    ## ... Similar to previous best
    ## Run 393 stress 0.06402618 
    ## Run 394 stress 0.06330869 
    ## ... Procrustes: rmse 2.384624e-05  max resid 7.768232e-05 
    ## ... Similar to previous best
    ## Run 395 stress 0.3698488 
    ## Run 396 stress 0.06330868 
    ## ... Procrustes: rmse 6.801639e-06  max resid 1.654589e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.0640262 
    ## Run 398 stress 0.09283183 
    ## Run 399 stress 0.06402616 
    ## Run 400 stress 0.06402634 
    ## Run 401 stress 0.09416694 
    ## Run 402 stress 0.09876246 
    ## Run 403 stress 0.08126176 
    ## Run 404 stress 0.06402619 
    ## Run 405 stress 0.09504635 
    ## Run 406 stress 0.06463313 
    ## Run 407 stress 0.08089032 
    ## Run 408 stress 0.06463313 
    ## Run 409 stress 0.06402633 
    ## Run 410 stress 0.08126184 
    ## Run 411 stress 0.06330868 
    ## ... Procrustes: rmse 1.646424e-05  max resid 3.785977e-05 
    ## ... Similar to previous best
    ## Run 412 stress 0.08403581 
    ## Run 413 stress 0.06330869 
    ## ... Procrustes: rmse 4.649447e-05  max resid 0.000103295 
    ## ... Similar to previous best
    ## Run 414 stress 0.06463313 
    ## Run 415 stress 0.06463313 
    ## Run 416 stress 0.06330869 
    ## ... Procrustes: rmse 4.141742e-05  max resid 9.217697e-05 
    ## ... Similar to previous best
    ## Run 417 stress 0.06330869 
    ## ... Procrustes: rmse 2.900634e-05  max resid 4.861881e-05 
    ## ... Similar to previous best
    ## Run 418 stress 0.06463316 
    ## Run 419 stress 0.06463315 
    ## Run 420 stress 0.06330872 
    ## ... Procrustes: rmse 7.183303e-05  max resid 0.000153474 
    ## ... Similar to previous best
    ## Run 421 stress 0.06402642 
    ## Run 422 stress 0.06463313 
    ## Run 423 stress 0.06402617 
    ## Run 424 stress 0.06402634 
    ## Run 425 stress 0.06330869 
    ## ... Procrustes: rmse 4.644643e-05  max resid 0.0001041542 
    ## ... Similar to previous best
    ## Run 426 stress 0.06330868 
    ## ... Procrustes: rmse 2.310052e-05  max resid 4.039813e-05 
    ## ... Similar to previous best
    ## Run 427 stress 0.06330868 
    ## ... Procrustes: rmse 6.315059e-06  max resid 1.188068e-05 
    ## ... Similar to previous best
    ## Run 428 stress 0.08089011 
    ## Run 429 stress 0.06402617 
    ## Run 430 stress 0.06463314 
    ## Run 431 stress 0.08126186 
    ## Run 432 stress 0.06402617 
    ## Run 433 stress 0.06463314 
    ## Run 434 stress 0.06463317 
    ## Run 435 stress 0.06463313 
    ## Run 436 stress 0.06330868 
    ## ... Procrustes: rmse 3.507855e-06  max resid 1.08792e-05 
    ## ... Similar to previous best
    ## Run 437 stress 0.06463317 
    ## Run 438 stress 0.3744184 
    ## Run 439 stress 0.0640262 
    ## Run 440 stress 0.06463314 
    ## Run 441 stress 0.081262 
    ## Run 442 stress 0.08126182 
    ## Run 443 stress 0.08089049 
    ## Run 444 stress 0.06330871 
    ## ... Procrustes: rmse 7.183638e-05  max resid 0.0001604853 
    ## ... Similar to previous best
    ## Run 445 stress 0.06402627 
    ## Run 446 stress 0.06330871 
    ## ... Procrustes: rmse 4.700904e-05  max resid 9.949623e-05 
    ## ... Similar to previous best
    ## Run 447 stress 0.06463314 
    ## Run 448 stress 0.09503217 
    ## Run 449 stress 0.08126174 
    ## Run 450 stress 0.06463313 
    ## Run 451 stress 0.06463314 
    ## Run 452 stress 0.0633087 
    ## ... Procrustes: rmse 5.247811e-05  max resid 0.0001220988 
    ## ... Similar to previous best
    ## Run 453 stress 0.06330869 
    ## ... Procrustes: rmse 3.310996e-05  max resid 6.698972e-05 
    ## ... Similar to previous best
    ## Run 454 stress 0.06402623 
    ## Run 455 stress 0.06463314 
    ## Run 456 stress 0.06330869 
    ## ... Procrustes: rmse 2.537705e-05  max resid 5.63088e-05 
    ## ... Similar to previous best
    ## Run 457 stress 0.06402635 
    ## Run 458 stress 0.06402617 
    ## Run 459 stress 0.06330868 
    ## ... Procrustes: rmse 1.288371e-05  max resid 2.852554e-05 
    ## ... Similar to previous best
    ## Run 460 stress 0.06402616 
    ## Run 461 stress 0.06402622 
    ## Run 462 stress 0.08315604 
    ## Run 463 stress 0.06330868 
    ## ... Procrustes: rmse 2.110105e-05  max resid 4.097496e-05 
    ## ... Similar to previous best
    ## Run 464 stress 0.06402621 
    ## Run 465 stress 0.06330868 
    ## ... Procrustes: rmse 8.97101e-06  max resid 2.10751e-05 
    ## ... Similar to previous best
    ## Run 466 stress 0.06330868 
    ## ... Procrustes: rmse 4.292625e-06  max resid 8.740485e-06 
    ## ... Similar to previous best
    ## Run 467 stress 0.06463315 
    ## Run 468 stress 0.09533279 
    ## Run 469 stress 0.06463313 
    ## Run 470 stress 0.06463313 
    ## Run 471 stress 0.06402625 
    ## Run 472 stress 0.06463314 
    ## Run 473 stress 0.09416678 
    ## Run 474 stress 0.0633087 
    ## ... Procrustes: rmse 4.571167e-05  max resid 0.0001243146 
    ## ... Similar to previous best
    ## Run 475 stress 0.06330869 
    ## ... Procrustes: rmse 4.354794e-05  max resid 9.56596e-05 
    ## ... Similar to previous best
    ## Run 476 stress 0.06463313 
    ## Run 477 stress 0.06330869 
    ## ... Procrustes: rmse 3.556477e-05  max resid 0.0001080177 
    ## ... Similar to previous best
    ## Run 478 stress 0.06330868 
    ## ... Procrustes: rmse 3.544063e-06  max resid 1.164146e-05 
    ## ... Similar to previous best
    ## Run 479 stress 0.08327879 
    ## Run 480 stress 0.08379064 
    ## Run 481 stress 0.06463314 
    ## Run 482 stress 0.06402624 
    ## Run 483 stress 0.06463313 
    ## Run 484 stress 0.09206634 
    ## Run 485 stress 0.06463313 
    ## Run 486 stress 0.06330875 
    ## ... Procrustes: rmse 0.000106249  max resid 0.0002379039 
    ## ... Similar to previous best
    ## Run 487 stress 0.06402619 
    ## Run 488 stress 0.06330868 
    ## ... Procrustes: rmse 4.126754e-06  max resid 9.188246e-06 
    ## ... Similar to previous best
    ## Run 489 stress 0.06330868 
    ## ... Procrustes: rmse 8.651199e-06  max resid 1.717811e-05 
    ## ... Similar to previous best
    ## Run 490 stress 0.0633087 
    ## ... Procrustes: rmse 3.100715e-05  max resid 6.626147e-05 
    ## ... Similar to previous best
    ## Run 491 stress 0.06463313 
    ## Run 492 stress 0.06330873 
    ## ... Procrustes: rmse 9.454529e-05  max resid 0.000211181 
    ## ... Similar to previous best
    ## Run 493 stress 0.06402617 
    ## Run 494 stress 0.06402623 
    ## Run 495 stress 0.06463313 
    ## Run 496 stress 0.06402623 
    ## Run 497 stress 0.06402616 
    ## Run 498 stress 0.0923344 
    ## Run 499 stress 0.06330868 
    ## ... Procrustes: rmse 1.276939e-05  max resid 2.843927e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.06463313 
    ## *** Best solution repeated 52 times

``` r
# Mixed and stratified lakes
PD_beta_env_MS_NMDS <- metaMDS(PD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.06409272 
    ## Run 2 stress 0.05036799 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002906372  max resid 0.0007417951 
    ## ... Similar to previous best
    ## Run 3 stress 0.05171059 
    ## Run 4 stress 0.06136676 
    ## Run 5 stress 0.06231308 
    ## Run 6 stress 0.05036793 
    ## ... New best solution
    ## ... Procrustes: rmse 7.823109e-05  max resid 0.0001944265 
    ## ... Similar to previous best
    ## Run 7 stress 0.05602044 
    ## Run 8 stress 0.05036792 
    ## ... New best solution
    ## ... Procrustes: rmse 8.16458e-05  max resid 0.0001937319 
    ## ... Similar to previous best
    ## Run 9 stress 0.05036796 
    ## ... Procrustes: rmse 8.041976e-05  max resid 0.0002140609 
    ## ... Similar to previous best
    ## Run 10 stress 0.0582994 
    ## Run 11 stress 0.06730219 
    ## Run 12 stress 0.0670504 
    ## Run 13 stress 0.05171063 
    ## Run 14 stress 0.05171059 
    ## Run 15 stress 0.05337578 
    ## Run 16 stress 0.05337578 
    ## Run 17 stress 0.3592577 
    ## Run 18 stress 0.06478524 
    ## Run 19 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001315086  max resid 0.0003586773 
    ## ... Similar to previous best
    ## Run 20 stress 0.05051335 
    ## ... Procrustes: rmse 0.03573655  max resid 0.1220955 
    ## Run 21 stress 0.05829937 
    ## Run 22 stress 0.05036792 
    ## ... New best solution
    ## ... Procrustes: rmse 1.674291e-05  max resid 2.714565e-05 
    ## ... Similar to previous best
    ## Run 23 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001432754  max resid 0.0003491434 
    ## ... Similar to previous best
    ## Run 24 stress 0.05051281 
    ## ... Procrustes: rmse 0.03577893  max resid 0.1221338 
    ## Run 25 stress 0.05337597 
    ## Run 26 stress 0.06640333 
    ## Run 27 stress 0.05171067 
    ## Run 28 stress 0.05403572 
    ## Run 29 stress 0.05337586 
    ## Run 30 stress 0.05051332 
    ## ... Procrustes: rmse 0.03566606  max resid 0.1217217 
    ## Run 31 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001890065  max resid 0.0004612658 
    ## ... Similar to previous best
    ## Run 32 stress 0.053376 
    ## Run 33 stress 0.06640355 
    ## Run 34 stress 0.06111772 
    ## Run 35 stress 0.05602048 
    ## Run 36 stress 0.06231311 
    ## Run 37 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 4.666737e-05  max resid 0.0001092646 
    ## ... Similar to previous best
    ## Run 38 stress 0.05602069 
    ## Run 39 stress 0.050368 
    ## ... Procrustes: rmse 0.0001226724  max resid 0.0002994271 
    ## ... Similar to previous best
    ## Run 40 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001015677  max resid 0.0002459259 
    ## ... Similar to previous best
    ## Run 41 stress 0.06885076 
    ## Run 42 stress 0.06478537 
    ## Run 43 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001098312  max resid 0.0002741295 
    ## ... Similar to previous best
    ## Run 44 stress 0.06362793 
    ## Run 45 stress 0.06478526 
    ## Run 46 stress 0.0517106 
    ## Run 47 stress 0.05051332 
    ## ... Procrustes: rmse 0.03572026  max resid 0.1220391 
    ## Run 48 stress 0.05968632 
    ## Run 49 stress 0.05051321 
    ## ... Procrustes: rmse 0.03566204  max resid 0.1217091 
    ## Run 50 stress 0.05928248 
    ## Run 51 stress 0.05466703 
    ## Run 52 stress 0.05466701 
    ## Run 53 stress 0.06904514 
    ## Run 54 stress 0.05051296 
    ## ... Procrustes: rmse 0.0356844  max resid 0.1219139 
    ## Run 55 stress 0.05051322 
    ## ... Procrustes: rmse 0.03566843  max resid 0.1217268 
    ## Run 56 stress 0.05968629 
    ## Run 57 stress 0.06027432 
    ## Run 58 stress 0.0582993 
    ## Run 59 stress 0.05403622 
    ## Run 60 stress 0.06231286 
    ## Run 61 stress 0.05036791 
    ## ... Procrustes: rmse 2.958542e-05  max resid 5.824988e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.05968626 
    ## Run 63 stress 0.05036793 
    ## ... Procrustes: rmse 5.55805e-05  max resid 0.0001392987 
    ## ... Similar to previous best
    ## Run 64 stress 0.05466701 
    ## Run 65 stress 0.05036798 
    ## ... Procrustes: rmse 0.00013155  max resid 0.0003687896 
    ## ... Similar to previous best
    ## Run 66 stress 0.06705124 
    ## Run 67 stress 0.06231311 
    ## Run 68 stress 0.05051304 
    ## ... Procrustes: rmse 0.03569482  max resid 0.1219502 
    ## Run 69 stress 0.05979226 
    ## Run 70 stress 0.06730248 
    ## Run 71 stress 0.05051317 
    ## ... Procrustes: rmse 0.0356512  max resid 0.121681 
    ## Run 72 stress 0.06478534 
    ## Run 73 stress 0.0533758 
    ## Run 74 stress 0.05171059 
    ## Run 75 stress 0.05036794 
    ## ... Procrustes: rmse 6.963087e-05  max resid 0.000175223 
    ## ... Similar to previous best
    ## Run 76 stress 0.05051328 
    ## ... Procrustes: rmse 0.03565676  max resid 0.1216861 
    ## Run 77 stress 0.05466701 
    ## Run 78 stress 0.05829939 
    ## Run 79 stress 0.06309846 
    ## Run 80 stress 0.05036793 
    ## ... Procrustes: rmse 5.717811e-05  max resid 0.0001406489 
    ## ... Similar to previous best
    ## Run 81 stress 0.05466707 
    ## Run 82 stress 0.06027434 
    ## Run 83 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001882069  max resid 0.0005077471 
    ## ... Similar to previous best
    ## Run 84 stress 0.06111766 
    ## Run 85 stress 0.05403591 
    ## Run 86 stress 0.05829929 
    ## Run 87 stress 0.05051299 
    ## ... Procrustes: rmse 0.03567946  max resid 0.1217852 
    ## Run 88 stress 0.05979225 
    ## Run 89 stress 0.06493686 
    ## Run 90 stress 0.0596863 
    ## Run 91 stress 0.06111766 
    ## Run 92 stress 0.06362797 
    ## Run 93 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001085738  max resid 0.0002624565 
    ## ... Similar to previous best
    ## Run 94 stress 0.05171059 
    ## Run 95 stress 0.0505135 
    ## ... Procrustes: rmse 0.03567374  max resid 0.121716 
    ## Run 96 stress 0.050368 
    ## ... Procrustes: rmse 0.0001375678  max resid 0.0003085471 
    ## ... Similar to previous best
    ## Run 97 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001185993  max resid 0.0002852831 
    ## ... Similar to previous best
    ## Run 98 stress 0.05171061 
    ## Run 99 stress 0.06647182 
    ## Run 100 stress 0.06362798 
    ## Run 101 stress 0.05036795 
    ## ... Procrustes: rmse 8.602845e-05  max resid 0.0002101049 
    ## ... Similar to previous best
    ## Run 102 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001442363  max resid 0.00040761 
    ## ... Similar to previous best
    ## Run 103 stress 0.05403588 
    ## Run 104 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001513011  max resid 0.0003762148 
    ## ... Similar to previous best
    ## Run 105 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001484091  max resid 0.0003875464 
    ## ... Similar to previous best
    ## Run 106 stress 0.05051294 
    ## ... Procrustes: rmse 0.0356819  max resid 0.1217993 
    ## Run 107 stress 0.05036794 
    ## ... Procrustes: rmse 7.994843e-05  max resid 0.000196284 
    ## ... Similar to previous best
    ## Run 108 stress 0.06231322 
    ## Run 109 stress 0.05051301 
    ## ... Procrustes: rmse 0.03570469  max resid 0.1219773 
    ## Run 110 stress 0.05928222 
    ## Run 111 stress 0.05928219 
    ## Run 112 stress 0.06038952 
    ## Run 113 stress 0.05051338 
    ## ... Procrustes: rmse 0.03570088  max resid 0.1218087 
    ## Run 114 stress 0.05337578 
    ## Run 115 stress 0.05968631 
    ## Run 116 stress 0.05602042 
    ## Run 117 stress 0.05403576 
    ## Run 118 stress 0.05051315 
    ## ... Procrustes: rmse 0.03568785  max resid 0.1219352 
    ## Run 119 stress 0.06231311 
    ## Run 120 stress 0.05829928 
    ## Run 121 stress 0.05171059 
    ## Run 122 stress 0.05928255 
    ## Run 123 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001324202  max resid 0.000336647 
    ## ... Similar to previous best
    ## Run 124 stress 0.0592822 
    ## Run 125 stress 0.05337579 
    ## Run 126 stress 0.05403587 
    ## Run 127 stress 0.050513 
    ## ... Procrustes: rmse 0.03568667  max resid 0.1218051 
    ## Run 128 stress 0.06904459 
    ## Run 129 stress 0.05036793 
    ## ... Procrustes: rmse 7.456191e-05  max resid 0.0001862682 
    ## ... Similar to previous best
    ## Run 130 stress 0.05171062 
    ## Run 131 stress 0.06309851 
    ## Run 132 stress 0.05036791 
    ## ... Procrustes: rmse 2.138875e-05  max resid 5.398154e-05 
    ## ... Similar to previous best
    ## Run 133 stress 0.05051307 
    ## ... Procrustes: rmse 0.03568449  max resid 0.1217905 
    ## Run 134 stress 0.07104212 
    ## Run 135 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001920678  max resid 0.0005022668 
    ## ... Similar to previous best
    ## Run 136 stress 0.05171061 
    ## Run 137 stress 0.05403615 
    ## Run 138 stress 0.05051345 
    ## ... Procrustes: rmse 0.03572632  max resid 0.1220607 
    ## Run 139 stress 0.05171058 
    ## Run 140 stress 0.05337585 
    ## Run 141 stress 0.06478535 
    ## Run 142 stress 0.05036798 
    ## ... Procrustes: rmse 0.000117473  max resid 0.0002878891 
    ## ... Similar to previous best
    ## Run 143 stress 0.05036797 
    ## ... Procrustes: rmse 9.848652e-05  max resid 0.0002442591 
    ## ... Similar to previous best
    ## Run 144 stress 0.06231315 
    ## Run 145 stress 0.05051328 
    ## ... Procrustes: rmse 0.03572263  max resid 0.1220452 
    ## Run 146 stress 0.05829931 
    ## Run 147 stress 0.06609056 
    ## Run 148 stress 0.05968641 
    ## Run 149 stress 0.05466701 
    ## Run 150 stress 0.06743851 
    ## Run 151 stress 0.07026125 
    ## Run 152 stress 0.05051301 
    ## ... Procrustes: rmse 0.03569574  max resid 0.12195 
    ## Run 153 stress 0.06309875 
    ## Run 154 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001157596  max resid 0.0002856891 
    ## ... Similar to previous best
    ## Run 155 stress 0.05403584 
    ## Run 156 stress 0.06231309 
    ## Run 157 stress 0.05403567 
    ## Run 158 stress 0.06309846 
    ## Run 159 stress 0.05036811 
    ## ... Procrustes: rmse 0.0001875708  max resid 0.0004672329 
    ## ... Similar to previous best
    ## Run 160 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001067672  max resid 0.000275823 
    ## ... Similar to previous best
    ## Run 161 stress 0.05171059 
    ## Run 162 stress 0.05602055 
    ## Run 163 stress 0.05337577 
    ## Run 164 stress 0.06640339 
    ## Run 165 stress 0.05337579 
    ## Run 166 stress 0.0603895 
    ## Run 167 stress 0.06111764 
    ## Run 168 stress 0.06730244 
    ## Run 169 stress 0.05466702 
    ## Run 170 stress 0.05051284 
    ## ... Procrustes: rmse 0.03559433  max resid 0.121622 
    ## Run 171 stress 0.05928231 
    ## Run 172 stress 0.06231293 
    ## Run 173 stress 0.05171058 
    ## Run 174 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001576808  max resid 0.0003930532 
    ## ... Similar to previous best
    ## Run 175 stress 0.05403614 
    ## Run 176 stress 0.05036795 
    ## ... Procrustes: rmse 5.24969e-05  max resid 0.0001405676 
    ## ... Similar to previous best
    ## Run 177 stress 0.05968646 
    ## Run 178 stress 0.05979232 
    ## Run 179 stress 0.05051344 
    ## ... Procrustes: rmse 0.0356677  max resid 0.1217131 
    ## Run 180 stress 0.06630524 
    ## Run 181 stress 0.05036792 
    ## ... Procrustes: rmse 3.640323e-05  max resid 8.583763e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.0517106 
    ## Run 183 stress 0.05900148 
    ## Run 184 stress 0.06309836 
    ## Run 185 stress 0.05036807 
    ## ... Procrustes: rmse 0.0001850268  max resid 0.0004597732 
    ## ... Similar to previous best
    ## Run 186 stress 0.05466707 
    ## Run 187 stress 0.05036795 
    ## ... Procrustes: rmse 9.240239e-05  max resid 0.0002275675 
    ## ... Similar to previous best
    ## Run 188 stress 0.05171058 
    ## Run 189 stress 0.05036793 
    ## ... Procrustes: rmse 6.813516e-05  max resid 0.0001972192 
    ## ... Similar to previous best
    ## Run 190 stress 0.05051337 
    ## ... Procrustes: rmse 0.03569028  max resid 0.1219511 
    ## Run 191 stress 0.05051336 
    ## ... Procrustes: rmse 0.03568921  max resid 0.1219487 
    ## Run 192 stress 0.05036794 
    ## ... Procrustes: rmse 8.154069e-05  max resid 0.0001631119 
    ## ... Similar to previous best
    ## Run 193 stress 0.05036794 
    ## ... Procrustes: rmse 7.675839e-05  max resid 0.0001887619 
    ## ... Similar to previous best
    ## Run 194 stress 0.05602048 
    ## Run 195 stress 0.05051348 
    ## ... Procrustes: rmse 0.03569548  max resid 0.1219715 
    ## Run 196 stress 0.05403589 
    ## Run 197 stress 0.05036808 
    ## ... Procrustes: rmse 0.0001899184  max resid 0.0003606643 
    ## ... Similar to previous best
    ## Run 198 stress 0.05036795 
    ## ... Procrustes: rmse 7.115779e-05  max resid 0.0001975839 
    ## ... Similar to previous best
    ## Run 199 stress 0.06309844 
    ## Run 200 stress 0.05036791 
    ## ... Procrustes: rmse 1.522204e-05  max resid 3.341225e-05 
    ## ... Similar to previous best
    ## Run 201 stress 0.05051301 
    ## ... Procrustes: rmse 0.03569348  max resid 0.1219443 
    ## Run 202 stress 0.05171059 
    ## Run 203 stress 0.0592825 
    ## Run 204 stress 0.06309837 
    ## Run 205 stress 0.05968639 
    ## Run 206 stress 0.05403554 
    ## Run 207 stress 0.05979225 
    ## Run 208 stress 0.05036798 
    ## ... Procrustes: rmse 9.026386e-05  max resid 0.0002252146 
    ## ... Similar to previous best
    ## Run 209 stress 0.05968627 
    ## Run 210 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001277095  max resid 0.0003525861 
    ## ... Similar to previous best
    ## Run 211 stress 0.0597923 
    ## Run 212 stress 0.05036794 
    ## ... Procrustes: rmse 7.774929e-05  max resid 0.000190406 
    ## ... Similar to previous best
    ## Run 213 stress 0.05036794 
    ## ... Procrustes: rmse 6.597406e-05  max resid 0.0001657174 
    ## ... Similar to previous best
    ## Run 214 stress 0.05829947 
    ## Run 215 stress 0.05036792 
    ## ... Procrustes: rmse 4.923582e-05  max resid 0.0001017091 
    ## ... Similar to previous best
    ## Run 216 stress 0.05036794 
    ## ... Procrustes: rmse 5.626399e-05  max resid 0.0001245378 
    ## ... Similar to previous best
    ## Run 217 stress 0.05466703 
    ## Run 218 stress 0.06478512 
    ## Run 219 stress 0.05036797 
    ## ... Procrustes: rmse 9.133664e-05  max resid 0.0002223058 
    ## ... Similar to previous best
    ## Run 220 stress 0.06111771 
    ## Run 221 stress 0.06038943 
    ## Run 222 stress 0.06476233 
    ## Run 223 stress 0.06027435 
    ## Run 224 stress 0.06655657 
    ## Run 225 stress 0.05036791 
    ## ... Procrustes: rmse 4.407895e-05  max resid 8.53306e-05 
    ## ... Similar to previous best
    ## Run 226 stress 0.05171061 
    ## Run 227 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001607761  max resid 0.0004215881 
    ## ... Similar to previous best
    ## Run 228 stress 0.05968628 
    ## Run 229 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001388513  max resid 0.0003440053 
    ## ... Similar to previous best
    ## Run 230 stress 0.05829929 
    ## Run 231 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001019715  max resid 0.0002689729 
    ## ... Similar to previous best
    ## Run 232 stress 0.05829941 
    ## Run 233 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001680979  max resid 0.0004369435 
    ## ... Similar to previous best
    ## Run 234 stress 0.05051302 
    ## ... Procrustes: rmse 0.03572714  max resid 0.1219286 
    ## Run 235 stress 0.05403584 
    ## Run 236 stress 0.05900154 
    ## Run 237 stress 0.05403582 
    ## Run 238 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001384391  max resid 0.0003659226 
    ## ... Similar to previous best
    ## Run 239 stress 0.05171058 
    ## Run 240 stress 0.06716809 
    ## Run 241 stress 0.0667457 
    ## Run 242 stress 0.05602062 
    ## Run 243 stress 0.05036795 
    ## ... Procrustes: rmse 7.45485e-05  max resid 0.0001707388 
    ## ... Similar to previous best
    ## Run 244 stress 0.05602052 
    ## Run 245 stress 0.06309861 
    ## Run 246 stress 0.05051325 
    ## ... Procrustes: rmse 0.03567647  max resid 0.1217474 
    ## Run 247 stress 0.06231286 
    ## Run 248 stress 0.050368 
    ## ... Procrustes: rmse 0.0001551009  max resid 0.0004442945 
    ## ... Similar to previous best
    ## Run 249 stress 0.05968637 
    ## Run 250 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001015599  max resid 0.0002506414 
    ## ... Similar to previous best
    ## Run 251 stress 0.05171059 
    ## Run 252 stress 0.05036792 
    ## ... Procrustes: rmse 2.673402e-05  max resid 6.64463e-05 
    ## ... Similar to previous best
    ## Run 253 stress 0.05036796 
    ## ... Procrustes: rmse 9.851404e-05  max resid 0.0002403817 
    ## ... Similar to previous best
    ## Run 254 stress 0.05036792 
    ## ... Procrustes: rmse 4.630665e-05  max resid 0.0001105304 
    ## ... Similar to previous best
    ## Run 255 stress 0.05466701 
    ## Run 256 stress 0.05051299 
    ## ... Procrustes: rmse 0.03582318  max resid 0.1223039 
    ## Run 257 stress 0.0613668 
    ## Run 258 stress 0.05036794 
    ## ... Procrustes: rmse 9.369725e-05  max resid 0.0002482382 
    ## ... Similar to previous best
    ## Run 259 stress 0.05051308 
    ## ... Procrustes: rmse 0.03569808  max resid 0.1219623 
    ## Run 260 stress 0.05466701 
    ## Run 261 stress 0.05036791 
    ## ... Procrustes: rmse 3.114092e-05  max resid 8.832023e-05 
    ## ... Similar to previous best
    ## Run 262 stress 0.06027448 
    ## Run 263 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001193182  max resid 0.0003419566 
    ## ... Similar to previous best
    ## Run 264 stress 0.05051355 
    ## ... Procrustes: rmse 0.03569656  max resid 0.1219761 
    ## Run 265 stress 0.06231306 
    ## Run 266 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001692864  max resid 0.000466472 
    ## ... Similar to previous best
    ## Run 267 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001510042  max resid 0.0003661761 
    ## ... Similar to previous best
    ## Run 268 stress 0.06537641 
    ## Run 269 stress 0.05403564 
    ## Run 270 stress 0.05036796 
    ## ... Procrustes: rmse 9.8704e-05  max resid 0.0002891763 
    ## ... Similar to previous best
    ## Run 271 stress 0.06136681 
    ## Run 272 stress 0.06027436 
    ## Run 273 stress 0.05928214 
    ## Run 274 stress 0.05051335 
    ## ... Procrustes: rmse 0.03569821  max resid 0.121805 
    ## Run 275 stress 0.05337584 
    ## Run 276 stress 0.05036791 
    ## ... Procrustes: rmse 1.207272e-05  max resid 2.355267e-05 
    ## ... Similar to previous best
    ## Run 277 stress 0.05968632 
    ## Run 278 stress 0.05036807 
    ## ... Procrustes: rmse 0.0001845848  max resid 0.0004513185 
    ## ... Similar to previous best
    ## Run 279 stress 0.06231287 
    ## Run 280 stress 0.05337585 
    ## Run 281 stress 0.05051296 
    ## ... Procrustes: rmse 0.03568625  max resid 0.1219152 
    ## Run 282 stress 0.05466704 
    ## Run 283 stress 0.05171059 
    ## Run 284 stress 0.05171059 
    ## Run 285 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001680009  max resid 0.0004210545 
    ## ... Similar to previous best
    ## Run 286 stress 0.06231309 
    ## Run 287 stress 0.0533758 
    ## Run 288 stress 0.05051298 
    ## ... Procrustes: rmse 0.03567277  max resid 0.1217662 
    ## Run 289 stress 0.05602049 
    ## Run 290 stress 0.0582994 
    ## Run 291 stress 0.06203052 
    ## Run 292 stress 0.05036792 
    ## ... Procrustes: rmse 3.526335e-05  max resid 8.355534e-05 
    ## ... Similar to previous best
    ## Run 293 stress 0.05466701 
    ## Run 294 stress 0.06647181 
    ## Run 295 stress 0.05403597 
    ## Run 296 stress 0.05337589 
    ## Run 297 stress 0.06203032 
    ## Run 298 stress 0.05036818 
    ## ... Procrustes: rmse 0.0002129457  max resid 0.0005403725 
    ## ... Similar to previous best
    ## Run 299 stress 0.05171062 
    ## Run 300 stress 0.0505135 
    ## ... Procrustes: rmse 0.03573978  max resid 0.1221043 
    ## Run 301 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001742682  max resid 0.0004519178 
    ## ... Similar to previous best
    ## Run 302 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001785365  max resid 0.0004613965 
    ## ... Similar to previous best
    ## Run 303 stress 0.05337577 
    ## Run 304 stress 0.05171059 
    ## Run 305 stress 0.05036796 
    ## ... Procrustes: rmse 9.151299e-05  max resid 0.0002265414 
    ## ... Similar to previous best
    ## Run 306 stress 0.06504237 
    ## Run 307 stress 0.06478527 
    ## Run 308 stress 0.05051323 
    ## ... Procrustes: rmse 0.03569875  max resid 0.1219717 
    ## Run 309 stress 0.05051346 
    ## ... Procrustes: rmse 0.035705  max resid 0.1218163 
    ## Run 310 stress 0.05928219 
    ## Run 311 stress 0.05928221 
    ## Run 312 stress 0.05036793 
    ## ... Procrustes: rmse 4.901686e-05  max resid 0.0001165849 
    ## ... Similar to previous best
    ## Run 313 stress 0.05051343 
    ## ... Procrustes: rmse 0.03571799  max resid 0.1220357 
    ## Run 314 stress 0.05403606 
    ## Run 315 stress 0.05900148 
    ## Run 316 stress 0.06837497 
    ## Run 317 stress 0.0533758 
    ## Run 318 stress 0.05036793 
    ## ... Procrustes: rmse 5.327069e-05  max resid 0.0001272605 
    ## ... Similar to previous best
    ## Run 319 stress 0.0505132 
    ## ... Procrustes: rmse 0.03566765  max resid 0.1217279 
    ## Run 320 stress 0.05466704 
    ## Run 321 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001425784  max resid 0.0003637946 
    ## ... Similar to previous best
    ## Run 322 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001204611  max resid 0.0003224935 
    ## ... Similar to previous best
    ## Run 323 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001394072  max resid 0.0003610936 
    ## ... Similar to previous best
    ## Run 324 stress 0.07033192 
    ## Run 325 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001200605  max resid 0.0002474779 
    ## ... Similar to previous best
    ## Run 326 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001725387  max resid 0.0004899304 
    ## ... Similar to previous best
    ## Run 327 stress 0.05036795 
    ## ... Procrustes: rmse 7.937264e-05  max resid 0.0001879843 
    ## ... Similar to previous best
    ## Run 328 stress 0.0596864 
    ## Run 329 stress 0.05171064 
    ## Run 330 stress 0.06231293 
    ## Run 331 stress 0.0560205 
    ## Run 332 stress 0.05466702 
    ## Run 333 stress 0.05036798 
    ## ... Procrustes: rmse 0.000117307  max resid 0.0002911501 
    ## ... Similar to previous best
    ## Run 334 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001617279  max resid 0.0004340535 
    ## ... Similar to previous best
    ## Run 335 stress 0.06136678 
    ## Run 336 stress 0.06231294 
    ## Run 337 stress 0.05036791 
    ## ... Procrustes: rmse 2.640276e-05  max resid 6.212629e-05 
    ## ... Similar to previous best
    ## Run 338 stress 0.05466702 
    ## Run 339 stress 0.05403578 
    ## Run 340 stress 0.05036793 
    ## ... Procrustes: rmse 6.314952e-05  max resid 0.0001535411 
    ## ... Similar to previous best
    ## Run 341 stress 0.05051322 
    ## ... Procrustes: rmse 0.03570508  max resid 0.1219896 
    ## Run 342 stress 0.05928241 
    ## Run 343 stress 0.05602072 
    ## Run 344 stress 0.05036794 
    ## ... Procrustes: rmse 9.770515e-05  max resid 0.0002565432 
    ## ... Similar to previous best
    ## Run 345 stress 0.05337586 
    ## Run 346 stress 0.05337581 
    ## Run 347 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001024272  max resid 0.0002147211 
    ## ... Similar to previous best
    ## Run 348 stress 0.05602055 
    ## Run 349 stress 0.06231287 
    ## Run 350 stress 0.05036808 
    ## ... Procrustes: rmse 0.0001901709  max resid 0.0004754237 
    ## ... Similar to previous best
    ## Run 351 stress 0.05051335 
    ## ... Procrustes: rmse 0.03574486  max resid 0.1221141 
    ## Run 352 stress 0.05051282 
    ## ... Procrustes: rmse 0.03568272  max resid 0.1218947 
    ## Run 353 stress 0.06038943 
    ## Run 354 stress 0.05466703 
    ## Run 355 stress 0.06136681 
    ## Run 356 stress 0.05051314 
    ## ... Procrustes: rmse 0.03567994  max resid 0.1217693 
    ## Run 357 stress 0.06904485 
    ## Run 358 stress 0.05051365 
    ## ... Procrustes: rmse 0.03569714  max resid 0.1219792 
    ## Run 359 stress 0.05602046 
    ## Run 360 stress 0.05829943 
    ## Run 361 stress 0.06027436 
    ## Run 362 stress 0.06111769 
    ## Run 363 stress 0.05466706 
    ## Run 364 stress 0.0620305 
    ## Run 365 stress 0.05337585 
    ## Run 366 stress 0.05337584 
    ## Run 367 stress 0.05051288 
    ## ... Procrustes: rmse 0.0356748  max resid 0.1217881 
    ## Run 368 stress 0.05051311 
    ## ... Procrustes: rmse 0.03573217  max resid 0.1219316 
    ## Run 369 stress 0.0582994 
    ## Run 370 stress 0.05036794 
    ## ... Procrustes: rmse 2.547866e-05  max resid 5.503832e-05 
    ## ... Similar to previous best
    ## Run 371 stress 0.05051338 
    ## ... Procrustes: rmse 0.03564192  max resid 0.1216343 
    ## Run 372 stress 0.05403563 
    ## Run 373 stress 0.06478529 
    ## Run 374 stress 0.05968631 
    ## Run 375 stress 0.05968626 
    ## Run 376 stress 0.06027429 
    ## Run 377 stress 0.05036794 
    ## ... Procrustes: rmse 7.146038e-05  max resid 0.000153689 
    ## ... Similar to previous best
    ## Run 378 stress 0.05051337 
    ## ... Procrustes: rmse 0.03569719  max resid 0.1219724 
    ## Run 379 stress 0.05051306 
    ## ... Procrustes: rmse 0.0357282  max resid 0.121927 
    ## Run 380 stress 0.05403611 
    ## Run 381 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001225762  max resid 0.0003062965 
    ## ... Similar to previous best
    ## Run 382 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001323876  max resid 0.0003362269 
    ## ... Similar to previous best
    ## Run 383 stress 0.062313 
    ## Run 384 stress 0.06362787 
    ## Run 385 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001175195  max resid 0.0002922253 
    ## ... Similar to previous best
    ## Run 386 stress 0.05036793 
    ## ... Procrustes: rmse 8.754708e-05  max resid 0.0002189992 
    ## ... Similar to previous best
    ## Run 387 stress 0.0613668 
    ## Run 388 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001669418  max resid 0.000448395 
    ## ... Similar to previous best
    ## Run 389 stress 0.05968628 
    ## Run 390 stress 0.05036794 
    ## ... Procrustes: rmse 8.230412e-05  max resid 0.0002293023 
    ## ... Similar to previous best
    ## Run 391 stress 0.06027428 
    ## Run 392 stress 0.05337585 
    ## Run 393 stress 0.05051309 
    ## ... Procrustes: rmse 0.0356663  max resid 0.121736 
    ## Run 394 stress 0.05968633 
    ## Run 395 stress 0.0560205 
    ## Run 396 stress 0.05928224 
    ## Run 397 stress 0.05051333 
    ## ... Procrustes: rmse 0.03573197  max resid 0.1220752 
    ## Run 398 stress 0.05171058 
    ## Run 399 stress 0.06231305 
    ## Run 400 stress 0.05829939 
    ## Run 401 stress 0.05403594 
    ## Run 402 stress 0.05171059 
    ## Run 403 stress 0.06609057 
    ## Run 404 stress 0.05602049 
    ## Run 405 stress 0.05466701 
    ## Run 406 stress 0.06478516 
    ## Run 407 stress 0.0602744 
    ## Run 408 stress 0.05337586 
    ## Run 409 stress 0.05602071 
    ## Run 410 stress 0.06478516 
    ## Run 411 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001004748  max resid 0.0002435995 
    ## ... Similar to previous best
    ## Run 412 stress 0.06478551 
    ## Run 413 stress 0.05036794 
    ## ... Procrustes: rmse 0.000100859  max resid 0.000266251 
    ## ... Similar to previous best
    ## Run 414 stress 0.05036792 
    ## ... Procrustes: rmse 1.704067e-05  max resid 3.237861e-05 
    ## ... Similar to previous best
    ## Run 415 stress 0.05602052 
    ## Run 416 stress 0.06674549 
    ## Run 417 stress 0.05928245 
    ## Run 418 stress 0.05968633 
    ## Run 419 stress 0.05403602 
    ## Run 420 stress 0.06309845 
    ## Run 421 stress 0.06478538 
    ## Run 422 stress 0.05036806 
    ## ... Procrustes: rmse 0.0001766951  max resid 0.0004473792 
    ## ... Similar to previous best
    ## Run 423 stress 0.05466704 
    ## Run 424 stress 0.0517106 
    ## Run 425 stress 0.05051314 
    ## ... Procrustes: rmse 0.03568995  max resid 0.1219399 
    ## Run 426 stress 0.05036791 
    ## ... Procrustes: rmse 4.040956e-05  max resid 9.151818e-05 
    ## ... Similar to previous best
    ## Run 427 stress 0.05036808 
    ## ... Procrustes: rmse 0.0001903155  max resid 0.0004798905 
    ## ... Similar to previous best
    ## Run 428 stress 0.05051293 
    ## ... Procrustes: rmse 0.03569488  max resid 0.1218395 
    ## Run 429 stress 0.05036796 
    ## ... Procrustes: rmse 8.993702e-05  max resid 0.0002143729 
    ## ... Similar to previous best
    ## Run 430 stress 0.05051316 
    ## ... Procrustes: rmse 0.0356843  max resid 0.1219274 
    ## Run 431 stress 0.06674546 
    ## Run 432 stress 0.05337581 
    ## Run 433 stress 0.05403602 
    ## Run 434 stress 0.06478524 
    ## Run 435 stress 0.0582993 
    ## Run 436 stress 0.06309849 
    ## Run 437 stress 0.05036808 
    ## ... Procrustes: rmse 0.0001903248  max resid 0.0004918223 
    ## ... Similar to previous best
    ## Run 438 stress 0.06111766 
    ## Run 439 stress 0.05337584 
    ## Run 440 stress 0.05051299 
    ## ... Procrustes: rmse 0.03568415  max resid 0.121915 
    ## Run 441 stress 0.05337578 
    ## Run 442 stress 0.06476217 
    ## Run 443 stress 0.06111768 
    ## Run 444 stress 0.05051326 
    ## ... Procrustes: rmse 0.03571656  max resid 0.1220254 
    ## Run 445 stress 0.05602049 
    ## Run 446 stress 0.05051308 
    ## ... Procrustes: rmse 0.03569976  max resid 0.1219677 
    ## Run 447 stress 0.06826696 
    ## Run 448 stress 0.05036795 
    ## ... Procrustes: rmse 8.842001e-05  max resid 0.0002176274 
    ## ... Similar to previous best
    ## Run 449 stress 0.05051345 
    ## ... Procrustes: rmse 0.03570338  max resid 0.1218114 
    ## Run 450 stress 0.0596863 
    ## Run 451 stress 0.05466703 
    ## Run 452 stress 0.06136678 
    ## Run 453 stress 0.05602068 
    ## Run 454 stress 0.0623129 
    ## Run 455 stress 0.050368 
    ## ... Procrustes: rmse 0.0001471072  max resid 0.0003819231 
    ## ... Similar to previous best
    ## Run 456 stress 0.06409318 
    ## Run 457 stress 0.06309853 
    ## Run 458 stress 0.05337591 
    ## Run 459 stress 0.05171059 
    ## Run 460 stress 0.06476214 
    ## Run 461 stress 0.06743841 
    ## Run 462 stress 0.05051294 
    ## ... Procrustes: rmse 0.0356703  max resid 0.1217645 
    ## Run 463 stress 0.05602046 
    ## Run 464 stress 0.0560206 
    ## Run 465 stress 0.05968636 
    ## Run 466 stress 0.05051327 
    ## ... Procrustes: rmse 0.03569256  max resid 0.1219551 
    ## Run 467 stress 0.05036794 
    ## ... Procrustes: rmse 8.527582e-05  max resid 0.0002348669 
    ## ... Similar to previous best
    ## Run 468 stress 0.05829931 
    ## Run 469 stress 0.05036792 
    ## ... Procrustes: rmse 3.003965e-05  max resid 7.239185e-05 
    ## ... Similar to previous best
    ## Run 470 stress 0.06478536 
    ## Run 471 stress 0.05602063 
    ## Run 472 stress 0.05051286 
    ## ... Procrustes: rmse 0.03568266  max resid 0.1218137 
    ## Run 473 stress 0.05051302 
    ## ... Procrustes: rmse 0.03568776  max resid 0.121927 
    ## Run 474 stress 0.05036792 
    ## ... Procrustes: rmse 5.837522e-05  max resid 0.0001519278 
    ## ... Similar to previous best
    ## Run 475 stress 0.05051338 
    ## ... Procrustes: rmse 0.03570327  max resid 0.1218172 
    ## Run 476 stress 0.05337582 
    ## Run 477 stress 0.06309853 
    ## Run 478 stress 0.05051335 
    ## ... Procrustes: rmse 0.03566624  max resid 0.1217079 
    ## Run 479 stress 0.05403552 
    ## Run 480 stress 0.06136678 
    ## Run 481 stress 0.05051296 
    ## ... Procrustes: rmse 0.03557674  max resid 0.1215794 
    ## Run 482 stress 0.05829929 
    ## Run 483 stress 0.05171058 
    ## Run 484 stress 0.05036794 
    ## ... Procrustes: rmse 7.45221e-05  max resid 0.0001836765 
    ## ... Similar to previous best
    ## Run 485 stress 0.05403557 
    ## Run 486 stress 0.06478544 
    ## Run 487 stress 0.05403573 
    ## Run 488 stress 0.05968636 
    ## Run 489 stress 0.05051286 
    ## ... Procrustes: rmse 0.03555126  max resid 0.1214656 
    ## Run 490 stress 0.05051296 
    ## ... Procrustes: rmse 0.03570206  max resid 0.1219667 
    ## Run 491 stress 0.06111766 
    ## Run 492 stress 0.06647178 
    ## Run 493 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001633174  max resid 0.0004014101 
    ## ... Similar to previous best
    ## Run 494 stress 0.05403587 
    ## Run 495 stress 0.05979227 
    ## Run 496 stress 0.0640926 
    ## Run 497 stress 0.05829932 
    ## Run 498 stress 0.06630535 
    ## Run 499 stress 0.06478536 
    ## Run 500 stress 0.06938756 
    ## *** Best solution repeated 109 times

``` r
# Ocean sites and mixed lakes
PD_beta_env_OM_NMDS <- metaMDS(PD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08352634 
    ## Run 1 stress 0.1134015 
    ## Run 2 stress 0.08352634 
    ## ... Procrustes: rmse 1.018336e-06  max resid 2.155041e-06 
    ## ... Similar to previous best
    ## Run 3 stress 0.1134014 
    ## Run 4 stress 0.1134014 
    ## Run 5 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 2.392178e-06  max resid 5.305114e-06 
    ## ... Similar to previous best
    ## Run 6 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.098146e-06  max resid 1.786005e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.08352634 
    ## ... Procrustes: rmse 1.554448e-06  max resid 3.337342e-06 
    ## ... Similar to previous best
    ## Run 8 stress 0.1134014 
    ## Run 9 stress 0.1134015 
    ## Run 10 stress 0.2031893 
    ## Run 11 stress 0.08352634 
    ## ... Procrustes: rmse 2.32057e-06  max resid 5.219433e-06 
    ## ... Similar to previous best
    ## Run 12 stress 0.1279309 
    ## Run 13 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.359069e-06  max resid 2.698155e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.08352634 
    ## ... Procrustes: rmse 1.941341e-06  max resid 3.835297e-06 
    ## ... Similar to previous best
    ## Run 15 stress 0.08352634 
    ## ... Procrustes: rmse 1.285832e-06  max resid 2.579018e-06 
    ## ... Similar to previous best
    ## Run 16 stress 0.1134014 
    ## Run 17 stress 0.1989154 
    ## Run 18 stress 0.08352634 
    ## ... Procrustes: rmse 4.793037e-06  max resid 1.04933e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.08352634 
    ## ... Procrustes: rmse 2.209615e-06  max resid 4.639465e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.1134014 
    ## Run 21 stress 0.1279309 
    ## Run 22 stress 0.1279309 
    ## Run 23 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 6.048376e-07  max resid 9.509763e-07 
    ## ... Similar to previous best
    ## Run 24 stress 0.08352634 
    ## ... Procrustes: rmse 1.725156e-06  max resid 3.386472e-06 
    ## ... Similar to previous best
    ## Run 25 stress 0.08352634 
    ## ... Procrustes: rmse 2.009413e-06  max resid 4.511762e-06 
    ## ... Similar to previous best
    ## Run 26 stress 0.08352634 
    ## ... Procrustes: rmse 1.56823e-06  max resid 3.463953e-06 
    ## ... Similar to previous best
    ## Run 27 stress 0.08352634 
    ## ... Procrustes: rmse 2.060126e-06  max resid 4.55789e-06 
    ## ... Similar to previous best
    ## Run 28 stress 0.08352634 
    ## ... Procrustes: rmse 1.179889e-06  max resid 2.20678e-06 
    ## ... Similar to previous best
    ## Run 29 stress 0.1134014 
    ## Run 30 stress 0.1134014 
    ## Run 31 stress 0.08352634 
    ## ... Procrustes: rmse 2.406123e-06  max resid 5.344467e-06 
    ## ... Similar to previous best
    ## Run 32 stress 0.1134015 
    ## Run 33 stress 0.08352634 
    ## ... Procrustes: rmse 2.943683e-06  max resid 6.476171e-06 
    ## ... Similar to previous best
    ## Run 34 stress 0.08352634 
    ## ... Procrustes: rmse 1.986021e-06  max resid 4.422774e-06 
    ## ... Similar to previous best
    ## Run 35 stress 0.1134014 
    ## Run 36 stress 0.08352634 
    ## ... Procrustes: rmse 2.81133e-06  max resid 6.241856e-06 
    ## ... Similar to previous best
    ## Run 37 stress 0.08352634 
    ## ... Procrustes: rmse 3.468833e-07  max resid 5.185558e-07 
    ## ... Similar to previous best
    ## Run 38 stress 0.1134015 
    ## Run 39 stress 0.08352634 
    ## ... Procrustes: rmse 1.586511e-06  max resid 3.618893e-06 
    ## ... Similar to previous best
    ## Run 40 stress 0.127931 
    ## Run 41 stress 0.127931 
    ## Run 42 stress 0.08352634 
    ## ... Procrustes: rmse 1.439262e-06  max resid 3.054904e-06 
    ## ... Similar to previous best
    ## Run 43 stress 0.08352634 
    ## ... Procrustes: rmse 1.955684e-06  max resid 4.111039e-06 
    ## ... Similar to previous best
    ## Run 44 stress 0.08352634 
    ## ... Procrustes: rmse 1.709647e-06  max resid 2.771171e-06 
    ## ... Similar to previous best
    ## Run 45 stress 0.1134014 
    ## Run 46 stress 0.08352634 
    ## ... Procrustes: rmse 1.641441e-06  max resid 3.614073e-06 
    ## ... Similar to previous best
    ## Run 47 stress 0.08352634 
    ## ... Procrustes: rmse 1.378476e-06  max resid 2.666118e-06 
    ## ... Similar to previous best
    ## Run 48 stress 0.08352634 
    ## ... Procrustes: rmse 2.278744e-06  max resid 5.034715e-06 
    ## ... Similar to previous best
    ## Run 49 stress 0.1279309 
    ## Run 50 stress 0.1134014 
    ## Run 51 stress 0.08352634 
    ## ... Procrustes: rmse 1.483232e-06  max resid 3.111e-06 
    ## ... Similar to previous best
    ## Run 52 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 9.147155e-07  max resid 2.079815e-06 
    ## ... Similar to previous best
    ## Run 53 stress 0.08352634 
    ## ... Procrustes: rmse 3.056832e-06  max resid 6.885762e-06 
    ## ... Similar to previous best
    ## Run 54 stress 0.1134015 
    ## Run 55 stress 0.1279309 
    ## Run 56 stress 0.08352634 
    ## ... Procrustes: rmse 9.976901e-07  max resid 2.196893e-06 
    ## ... Similar to previous best
    ## Run 57 stress 0.3271573 
    ## Run 58 stress 0.08352634 
    ## ... Procrustes: rmse 4.519969e-07  max resid 8.15918e-07 
    ## ... Similar to previous best
    ## Run 59 stress 0.1134014 
    ## Run 60 stress 0.1279309 
    ## Run 61 stress 0.08352634 
    ## ... Procrustes: rmse 3.173251e-06  max resid 7.08663e-06 
    ## ... Similar to previous best
    ## Run 62 stress 0.08352634 
    ## ... Procrustes: rmse 3.389232e-06  max resid 7.552647e-06 
    ## ... Similar to previous best
    ## Run 63 stress 0.08352634 
    ## ... Procrustes: rmse 1.756712e-06  max resid 3.328475e-06 
    ## ... Similar to previous best
    ## Run 64 stress 0.08352634 
    ## ... Procrustes: rmse 9.530179e-07  max resid 2.16248e-06 
    ## ... Similar to previous best
    ## Run 65 stress 0.08352634 
    ## ... Procrustes: rmse 1.678933e-06  max resid 3.785436e-06 
    ## ... Similar to previous best
    ## Run 66 stress 0.3266727 
    ## Run 67 stress 0.1279311 
    ## Run 68 stress 0.08352634 
    ## ... Procrustes: rmse 8.87805e-07  max resid 1.975944e-06 
    ## ... Similar to previous best
    ## Run 69 stress 0.1134015 
    ## Run 70 stress 0.08352634 
    ## ... Procrustes: rmse 2.714128e-06  max resid 6.006111e-06 
    ## ... Similar to previous best
    ## Run 71 stress 0.1134015 
    ## Run 72 stress 0.2031893 
    ## Run 73 stress 0.1783127 
    ## Run 74 stress 0.08352634 
    ## ... Procrustes: rmse 3.491659e-06  max resid 6.938822e-06 
    ## ... Similar to previous best
    ## Run 75 stress 0.08352634 
    ## ... Procrustes: rmse 2.397951e-06  max resid 3.647123e-06 
    ## ... Similar to previous best
    ## Run 76 stress 0.08352634 
    ## ... Procrustes: rmse 3.078733e-06  max resid 6.892404e-06 
    ## ... Similar to previous best
    ## Run 77 stress 0.1134014 
    ## Run 78 stress 0.08352634 
    ## ... Procrustes: rmse 1.417698e-06  max resid 2.990822e-06 
    ## ... Similar to previous best
    ## Run 79 stress 0.1134015 
    ## Run 80 stress 0.08352634 
    ## ... Procrustes: rmse 1.013e-06  max resid 2.379356e-06 
    ## ... Similar to previous best
    ## Run 81 stress 0.1134015 
    ## Run 82 stress 0.08352634 
    ## ... Procrustes: rmse 6.131738e-07  max resid 1.142707e-06 
    ## ... Similar to previous best
    ## Run 83 stress 0.1134014 
    ## Run 84 stress 0.08352634 
    ## ... Procrustes: rmse 1.569555e-06  max resid 2.864101e-06 
    ## ... Similar to previous best
    ## Run 85 stress 0.08352634 
    ## ... Procrustes: rmse 9.96921e-07  max resid 2.22569e-06 
    ## ... Similar to previous best
    ## Run 86 stress 0.08352634 
    ## ... Procrustes: rmse 5.374568e-07  max resid 1.001811e-06 
    ## ... Similar to previous best
    ## Run 87 stress 0.127931 
    ## Run 88 stress 0.1989485 
    ## Run 89 stress 0.1134015 
    ## Run 90 stress 0.08352634 
    ## ... Procrustes: rmse 3.169753e-06  max resid 7.349423e-06 
    ## ... Similar to previous best
    ## Run 91 stress 0.08352634 
    ## ... Procrustes: rmse 1.078485e-06  max resid 2.443778e-06 
    ## ... Similar to previous best
    ## Run 92 stress 0.1279309 
    ## Run 93 stress 0.1134014 
    ## Run 94 stress 0.1989486 
    ## Run 95 stress 0.08352634 
    ## ... Procrustes: rmse 3.020829e-06  max resid 6.770578e-06 
    ## ... Similar to previous best
    ## Run 96 stress 0.1279309 
    ## Run 97 stress 0.08352634 
    ## ... Procrustes: rmse 2.516227e-06  max resid 5.627665e-06 
    ## ... Similar to previous best
    ## Run 98 stress 0.08352634 
    ## ... Procrustes: rmse 1.951499e-06  max resid 4.387199e-06 
    ## ... Similar to previous best
    ## Run 99 stress 0.08352634 
    ## ... Procrustes: rmse 3.466053e-06  max resid 7.75286e-06 
    ## ... Similar to previous best
    ## Run 100 stress 0.3407345 
    ## Run 101 stress 0.1279309 
    ## Run 102 stress 0.08352634 
    ## ... Procrustes: rmse 3.765711e-07  max resid 7.477513e-07 
    ## ... Similar to previous best
    ## Run 103 stress 0.1279309 
    ## Run 104 stress 0.1134014 
    ## Run 105 stress 0.1279309 
    ## Run 106 stress 0.1134014 
    ## Run 107 stress 0.08352634 
    ## ... Procrustes: rmse 1.094611e-06  max resid 2.418983e-06 
    ## ... Similar to previous best
    ## Run 108 stress 0.3126576 
    ## Run 109 stress 0.08352634 
    ## ... Procrustes: rmse 2.48213e-06  max resid 5.578089e-06 
    ## ... Similar to previous best
    ## Run 110 stress 0.1279309 
    ## Run 111 stress 0.08352634 
    ## ... Procrustes: rmse 2.861016e-06  max resid 6.372948e-06 
    ## ... Similar to previous best
    ## Run 112 stress 0.08352634 
    ## ... Procrustes: rmse 3.422511e-06  max resid 7.615106e-06 
    ## ... Similar to previous best
    ## Run 113 stress 0.08352634 
    ## ... Procrustes: rmse 8.034325e-07  max resid 1.61291e-06 
    ## ... Similar to previous best
    ## Run 114 stress 0.08352634 
    ## ... Procrustes: rmse 3.03258e-06  max resid 6.836503e-06 
    ## ... Similar to previous best
    ## Run 115 stress 0.3266727 
    ## Run 116 stress 0.08352634 
    ## ... Procrustes: rmse 1.217833e-06  max resid 2.590116e-06 
    ## ... Similar to previous best
    ## Run 117 stress 0.08352634 
    ## ... Procrustes: rmse 2.712536e-06  max resid 4.301306e-06 
    ## ... Similar to previous best
    ## Run 118 stress 0.08352634 
    ## ... Procrustes: rmse 7.699079e-07  max resid 1.56379e-06 
    ## ... Similar to previous best
    ## Run 119 stress 0.1134014 
    ## Run 120 stress 0.1134015 
    ## Run 121 stress 0.3154655 
    ## Run 122 stress 0.1134014 
    ## Run 123 stress 0.1279309 
    ## Run 124 stress 0.08352634 
    ## ... Procrustes: rmse 2.520275e-06  max resid 5.172436e-06 
    ## ... Similar to previous best
    ## Run 125 stress 0.1279309 
    ## Run 126 stress 0.08352634 
    ## ... Procrustes: rmse 3.088245e-06  max resid 6.972291e-06 
    ## ... Similar to previous best
    ## Run 127 stress 0.08352634 
    ## ... Procrustes: rmse 1.349515e-06  max resid 2.716566e-06 
    ## ... Similar to previous best
    ## Run 128 stress 0.08352634 
    ## ... Procrustes: rmse 2.356946e-06  max resid 5.095682e-06 
    ## ... Similar to previous best
    ## Run 129 stress 0.08352634 
    ## ... Procrustes: rmse 1.608742e-06  max resid 3.542057e-06 
    ## ... Similar to previous best
    ## Run 130 stress 0.08352634 
    ## ... Procrustes: rmse 4.634995e-06  max resid 1.002024e-05 
    ## ... Similar to previous best
    ## Run 131 stress 0.127931 
    ## Run 132 stress 0.08352634 
    ## ... Procrustes: rmse 1.74527e-06  max resid 3.630041e-06 
    ## ... Similar to previous best
    ## Run 133 stress 0.08352634 
    ## ... Procrustes: rmse 1.058461e-05  max resid 2.350132e-05 
    ## ... Similar to previous best
    ## Run 134 stress 0.08352634 
    ## ... Procrustes: rmse 2.554096e-06  max resid 5.75418e-06 
    ## ... Similar to previous best
    ## Run 135 stress 0.08352634 
    ## ... Procrustes: rmse 1.332944e-06  max resid 2.907676e-06 
    ## ... Similar to previous best
    ## Run 136 stress 0.08352634 
    ## ... Procrustes: rmse 1.227063e-06  max resid 2.509447e-06 
    ## ... Similar to previous best
    ## Run 137 stress 0.08352634 
    ## ... Procrustes: rmse 3.497808e-07  max resid 7.430201e-07 
    ## ... Similar to previous best
    ## Run 138 stress 0.08352634 
    ## ... Procrustes: rmse 1.348954e-06  max resid 2.855188e-06 
    ## ... Similar to previous best
    ## Run 139 stress 0.127931 
    ## Run 140 stress 0.1279309 
    ## Run 141 stress 0.1279309 
    ## Run 142 stress 0.08352634 
    ## ... Procrustes: rmse 6.499384e-06  max resid 1.478405e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.08352634 
    ## ... Procrustes: rmse 1.849076e-06  max resid 4.157453e-06 
    ## ... Similar to previous best
    ## Run 144 stress 0.08352634 
    ## ... Procrustes: rmse 7.387867e-07  max resid 8.981924e-07 
    ## ... Similar to previous best
    ## Run 145 stress 0.08352634 
    ## ... Procrustes: rmse 4.941748e-07  max resid 9.234822e-07 
    ## ... Similar to previous best
    ## Run 146 stress 0.1134015 
    ## Run 147 stress 0.08352634 
    ## ... Procrustes: rmse 2.088194e-05  max resid 4.61297e-05 
    ## ... Similar to previous best
    ## Run 148 stress 0.08352634 
    ## ... Procrustes: rmse 1.031767e-06  max resid 1.898693e-06 
    ## ... Similar to previous best
    ## Run 149 stress 0.1134014 
    ## Run 150 stress 0.1134014 
    ## Run 151 stress 0.1279309 
    ## Run 152 stress 0.1134014 
    ## Run 153 stress 0.08352634 
    ## ... Procrustes: rmse 3.203227e-06  max resid 7.040039e-06 
    ## ... Similar to previous best
    ## Run 154 stress 0.08352634 
    ## ... Procrustes: rmse 6.781852e-07  max resid 1.753641e-06 
    ## ... Similar to previous best
    ## Run 155 stress 0.1279309 
    ## Run 156 stress 0.1134015 
    ## Run 157 stress 0.1134015 
    ## Run 158 stress 0.1279309 
    ## Run 159 stress 0.1279311 
    ## Run 160 stress 0.08352634 
    ## ... Procrustes: rmse 3.466024e-06  max resid 7.762677e-06 
    ## ... Similar to previous best
    ## Run 161 stress 0.08352634 
    ## ... Procrustes: rmse 1.237914e-06  max resid 2.803993e-06 
    ## ... Similar to previous best
    ## Run 162 stress 0.08352634 
    ## ... Procrustes: rmse 2.60037e-06  max resid 5.758392e-06 
    ## ... Similar to previous best
    ## Run 163 stress 0.08352634 
    ## ... Procrustes: rmse 1.992313e-06  max resid 5.072059e-06 
    ## ... Similar to previous best
    ## Run 164 stress 0.08352634 
    ## ... Procrustes: rmse 3.647258e-07  max resid 6.204447e-07 
    ## ... Similar to previous best
    ## Run 165 stress 0.3266727 
    ## Run 166 stress 0.1134015 
    ## Run 167 stress 0.08352634 
    ## ... Procrustes: rmse 4.131557e-06  max resid 9.232747e-06 
    ## ... Similar to previous best
    ## Run 168 stress 0.08352634 
    ## ... Procrustes: rmse 1.939218e-06  max resid 4.154894e-06 
    ## ... Similar to previous best
    ## Run 169 stress 0.08352634 
    ## ... Procrustes: rmse 2.995929e-06  max resid 6.630558e-06 
    ## ... Similar to previous best
    ## Run 170 stress 0.08352634 
    ## ... Procrustes: rmse 3.503862e-07  max resid 6.742547e-07 
    ## ... Similar to previous best
    ## Run 171 stress 0.1134014 
    ## Run 172 stress 0.127931 
    ## Run 173 stress 0.1279309 
    ## Run 174 stress 0.08352634 
    ## ... Procrustes: rmse 2.667955e-06  max resid 5.592081e-06 
    ## ... Similar to previous best
    ## Run 175 stress 0.08352634 
    ## ... Procrustes: rmse 2.5893e-06  max resid 6.125548e-06 
    ## ... Similar to previous best
    ## Run 176 stress 0.1279309 
    ## Run 177 stress 0.1134015 
    ## Run 178 stress 0.08352634 
    ## ... Procrustes: rmse 8.733518e-07  max resid 1.967663e-06 
    ## ... Similar to previous best
    ## Run 179 stress 0.1134014 
    ## Run 180 stress 0.08352634 
    ## ... Procrustes: rmse 1.7449e-06  max resid 3.133467e-06 
    ## ... Similar to previous best
    ## Run 181 stress 0.08352634 
    ## ... Procrustes: rmse 2.880006e-06  max resid 6.530925e-06 
    ## ... Similar to previous best
    ## Run 182 stress 0.1279309 
    ## Run 183 stress 0.1134014 
    ## Run 184 stress 0.08352634 
    ## ... Procrustes: rmse 9.789198e-07  max resid 2.207155e-06 
    ## ... Similar to previous best
    ## Run 185 stress 0.1134015 
    ## Run 186 stress 0.08352634 
    ## ... Procrustes: rmse 6.343225e-07  max resid 1.2105e-06 
    ## ... Similar to previous best
    ## Run 187 stress 0.08352634 
    ## ... Procrustes: rmse 2.339034e-06  max resid 5.196285e-06 
    ## ... Similar to previous best
    ## Run 188 stress 0.08352634 
    ## ... Procrustes: rmse 1.135132e-06  max resid 2.534173e-06 
    ## ... Similar to previous best
    ## Run 189 stress 0.08352634 
    ## ... Procrustes: rmse 6.858509e-07  max resid 1.390949e-06 
    ## ... Similar to previous best
    ## Run 190 stress 0.1134015 
    ## Run 191 stress 0.1134015 
    ## Run 192 stress 0.1279309 
    ## Run 193 stress 0.08352634 
    ## ... Procrustes: rmse 2.118177e-06  max resid 4.912709e-06 
    ## ... Similar to previous best
    ## Run 194 stress 0.08352634 
    ## ... Procrustes: rmse 3.759541e-07  max resid 7.650466e-07 
    ## ... Similar to previous best
    ## Run 195 stress 0.08352634 
    ## ... Procrustes: rmse 3.937097e-07  max resid 6.529707e-07 
    ## ... Similar to previous best
    ## Run 196 stress 0.1134015 
    ## Run 197 stress 0.08352634 
    ## ... Procrustes: rmse 2.172371e-06  max resid 4.828742e-06 
    ## ... Similar to previous best
    ## Run 198 stress 0.08352634 
    ## ... Procrustes: rmse 1.900802e-06  max resid 3.857507e-06 
    ## ... Similar to previous best
    ## Run 199 stress 0.1783127 
    ## Run 200 stress 0.08352634 
    ## ... Procrustes: rmse 2.414628e-06  max resid 4.795423e-06 
    ## ... Similar to previous best
    ## Run 201 stress 0.08352634 
    ## ... Procrustes: rmse 2.055261e-06  max resid 4.735947e-06 
    ## ... Similar to previous best
    ## Run 202 stress 0.1783098 
    ## Run 203 stress 0.08352634 
    ## ... Procrustes: rmse 3.051416e-06  max resid 6.447249e-06 
    ## ... Similar to previous best
    ## Run 204 stress 0.08352634 
    ## ... Procrustes: rmse 1.442715e-06  max resid 2.112195e-06 
    ## ... Similar to previous best
    ## Run 205 stress 0.1989484 
    ## Run 206 stress 0.08352634 
    ## ... Procrustes: rmse 9.56137e-07  max resid 1.629202e-06 
    ## ... Similar to previous best
    ## Run 207 stress 0.08352634 
    ## ... Procrustes: rmse 6.330005e-07  max resid 1.152417e-06 
    ## ... Similar to previous best
    ## Run 208 stress 0.08352634 
    ## ... Procrustes: rmse 2.443212e-06  max resid 5.46174e-06 
    ## ... Similar to previous best
    ## Run 209 stress 0.127931 
    ## Run 210 stress 0.1134015 
    ## Run 211 stress 0.08352634 
    ## ... Procrustes: rmse 3.047924e-06  max resid 6.691219e-06 
    ## ... Similar to previous best
    ## Run 212 stress 0.08352634 
    ## ... Procrustes: rmse 2.038936e-06  max resid 4.336432e-06 
    ## ... Similar to previous best
    ## Run 213 stress 0.08352634 
    ## ... Procrustes: rmse 2.224093e-06  max resid 4.023327e-06 
    ## ... Similar to previous best
    ## Run 214 stress 0.08352634 
    ## ... Procrustes: rmse 1.291331e-06  max resid 2.802977e-06 
    ## ... Similar to previous best
    ## Run 215 stress 0.1134014 
    ## Run 216 stress 0.08352634 
    ## ... Procrustes: rmse 1.99287e-06  max resid 4.53742e-06 
    ## ... Similar to previous best
    ## Run 217 stress 0.1134014 
    ## Run 218 stress 0.08352634 
    ## ... Procrustes: rmse 1.050352e-06  max resid 2.004837e-06 
    ## ... Similar to previous best
    ## Run 219 stress 0.08352634 
    ## ... Procrustes: rmse 1.811273e-06  max resid 3.834156e-06 
    ## ... Similar to previous best
    ## Run 220 stress 0.08352634 
    ## ... Procrustes: rmse 1.779097e-06  max resid 3.854791e-06 
    ## ... Similar to previous best
    ## Run 221 stress 0.08352634 
    ## ... Procrustes: rmse 8.903706e-07  max resid 1.983234e-06 
    ## ... Similar to previous best
    ## Run 222 stress 0.08352634 
    ## ... Procrustes: rmse 4.070418e-06  max resid 9.087008e-06 
    ## ... Similar to previous best
    ## Run 223 stress 0.1134014 
    ## Run 224 stress 0.08352634 
    ## ... Procrustes: rmse 2.167757e-06  max resid 4.842036e-06 
    ## ... Similar to previous best
    ## Run 225 stress 0.08352634 
    ## ... Procrustes: rmse 1.051687e-06  max resid 1.963603e-06 
    ## ... Similar to previous best
    ## Run 226 stress 0.1134014 
    ## Run 227 stress 0.1134014 
    ## Run 228 stress 0.08352634 
    ## ... Procrustes: rmse 3.548112e-06  max resid 7.974132e-06 
    ## ... Similar to previous best
    ## Run 229 stress 0.3407337 
    ## Run 230 stress 0.127931 
    ## Run 231 stress 0.127931 
    ## Run 232 stress 0.08352634 
    ## ... Procrustes: rmse 3.792131e-06  max resid 8.475189e-06 
    ## ... Similar to previous best
    ## Run 233 stress 0.2031893 
    ## Run 234 stress 0.1134014 
    ## Run 235 stress 0.08352634 
    ## ... Procrustes: rmse 3.400906e-06  max resid 7.342648e-06 
    ## ... Similar to previous best
    ## Run 236 stress 0.08352634 
    ## ... Procrustes: rmse 2.135342e-06  max resid 4.735322e-06 
    ## ... Similar to previous best
    ## Run 237 stress 0.08352634 
    ## ... Procrustes: rmse 1.724876e-06  max resid 3.909499e-06 
    ## ... Similar to previous best
    ## Run 238 stress 0.1134014 
    ## Run 239 stress 0.08352634 
    ## ... Procrustes: rmse 3.061121e-06  max resid 6.886631e-06 
    ## ... Similar to previous best
    ## Run 240 stress 0.1989154 
    ## Run 241 stress 0.3266729 
    ## Run 242 stress 0.2050058 
    ## Run 243 stress 0.08352634 
    ## ... Procrustes: rmse 1.310586e-06  max resid 2.84775e-06 
    ## ... Similar to previous best
    ## Run 244 stress 0.08352634 
    ## ... Procrustes: rmse 1.589295e-06  max resid 3.20514e-06 
    ## ... Similar to previous best
    ## Run 245 stress 0.1134015 
    ## Run 246 stress 0.08352634 
    ## ... Procrustes: rmse 1.773989e-06  max resid 3.458393e-06 
    ## ... Similar to previous best
    ## Run 247 stress 0.1134015 
    ## Run 248 stress 0.1279309 
    ## Run 249 stress 0.08352634 
    ## ... Procrustes: rmse 2.308775e-06  max resid 5.204643e-06 
    ## ... Similar to previous best
    ## Run 250 stress 0.127931 
    ## Run 251 stress 0.08352634 
    ## ... Procrustes: rmse 3.349786e-06  max resid 4.77303e-06 
    ## ... Similar to previous best
    ## Run 252 stress 0.08352634 
    ## ... Procrustes: rmse 4.659162e-06  max resid 9.421718e-06 
    ## ... Similar to previous best
    ## Run 253 stress 0.08352634 
    ## ... Procrustes: rmse 2.960419e-06  max resid 6.661513e-06 
    ## ... Similar to previous best
    ## Run 254 stress 0.1134015 
    ## Run 255 stress 0.08352634 
    ## ... Procrustes: rmse 2.996697e-06  max resid 6.851933e-06 
    ## ... Similar to previous best
    ## Run 256 stress 0.08352634 
    ## ... Procrustes: rmse 7.363929e-06  max resid 1.634795e-05 
    ## ... Similar to previous best
    ## Run 257 stress 0.1279309 
    ## Run 258 stress 0.1279309 
    ## Run 259 stress 0.1134015 
    ## Run 260 stress 0.1134014 
    ## Run 261 stress 0.08352634 
    ## ... Procrustes: rmse 7.205559e-07  max resid 1.452047e-06 
    ## ... Similar to previous best
    ## Run 262 stress 0.1279309 
    ## Run 263 stress 0.08352634 
    ## ... Procrustes: rmse 5.780282e-07  max resid 9.440525e-07 
    ## ... Similar to previous best
    ## Run 264 stress 0.1134015 
    ## Run 265 stress 0.3403702 
    ## Run 266 stress 0.08352634 
    ## ... Procrustes: rmse 1.421913e-06  max resid 2.922878e-06 
    ## ... Similar to previous best
    ## Run 267 stress 0.08352634 
    ## ... Procrustes: rmse 2.481915e-06  max resid 5.588445e-06 
    ## ... Similar to previous best
    ## Run 268 stress 0.08352634 
    ## ... Procrustes: rmse 3.025106e-06  max resid 6.552124e-06 
    ## ... Similar to previous best
    ## Run 269 stress 0.1134014 
    ## Run 270 stress 0.08352634 
    ## ... Procrustes: rmse 1.848296e-06  max resid 4.117393e-06 
    ## ... Similar to previous best
    ## Run 271 stress 0.08352634 
    ## ... Procrustes: rmse 5.502435e-07  max resid 9.991934e-07 
    ## ... Similar to previous best
    ## Run 272 stress 0.08352634 
    ## ... Procrustes: rmse 1.290009e-06  max resid 2.414142e-06 
    ## ... Similar to previous best
    ## Run 273 stress 0.1134014 
    ## Run 274 stress 0.1279309 
    ## Run 275 stress 0.08352634 
    ## ... Procrustes: rmse 2.927652e-06  max resid 6.587904e-06 
    ## ... Similar to previous best
    ## Run 276 stress 0.08352634 
    ## ... Procrustes: rmse 1.873289e-06  max resid 3.932719e-06 
    ## ... Similar to previous best
    ## Run 277 stress 0.08352634 
    ## ... Procrustes: rmse 5.193298e-07  max resid 1.00357e-06 
    ## ... Similar to previous best
    ## Run 278 stress 0.1989486 
    ## Run 279 stress 0.08352634 
    ## ... Procrustes: rmse 1.069643e-06  max resid 2.331115e-06 
    ## ... Similar to previous best
    ## Run 280 stress 0.08352634 
    ## ... Procrustes: rmse 9.135672e-07  max resid 2.01465e-06 
    ## ... Similar to previous best
    ## Run 281 stress 0.08352634 
    ## ... Procrustes: rmse 1.417997e-06  max resid 3.124491e-06 
    ## ... Similar to previous best
    ## Run 282 stress 0.127931 
    ## Run 283 stress 0.08352634 
    ## ... Procrustes: rmse 3.303066e-06  max resid 7.360386e-06 
    ## ... Similar to previous best
    ## Run 284 stress 0.3154698 
    ## Run 285 stress 0.3070017 
    ## Run 286 stress 0.1279309 
    ## Run 287 stress 0.08352634 
    ## ... Procrustes: rmse 2.262225e-06  max resid 4.825022e-06 
    ## ... Similar to previous best
    ## Run 288 stress 0.1134015 
    ## Run 289 stress 0.1279309 
    ## Run 290 stress 0.1134015 
    ## Run 291 stress 0.1134015 
    ## Run 292 stress 0.08352634 
    ## ... Procrustes: rmse 6.61078e-07  max resid 9.21934e-07 
    ## ... Similar to previous best
    ## Run 293 stress 0.08352634 
    ## ... Procrustes: rmse 2.788425e-06  max resid 6.055583e-06 
    ## ... Similar to previous best
    ## Run 294 stress 0.08352634 
    ## ... Procrustes: rmse 3.211465e-06  max resid 6.97714e-06 
    ## ... Similar to previous best
    ## Run 295 stress 0.08352634 
    ## ... Procrustes: rmse 1.661803e-06  max resid 3.688417e-06 
    ## ... Similar to previous best
    ## Run 296 stress 0.1134015 
    ## Run 297 stress 0.08352634 
    ## ... Procrustes: rmse 1.050578e-05  max resid 2.33078e-05 
    ## ... Similar to previous best
    ## Run 298 stress 0.08352634 
    ## ... Procrustes: rmse 3.111888e-06  max resid 6.849901e-06 
    ## ... Similar to previous best
    ## Run 299 stress 0.1134015 
    ## Run 300 stress 0.08352634 
    ## ... Procrustes: rmse 1.545384e-06  max resid 3.352661e-06 
    ## ... Similar to previous best
    ## Run 301 stress 0.3270102 
    ## Run 302 stress 0.08352634 
    ## ... Procrustes: rmse 1.043845e-06  max resid 2.274251e-06 
    ## ... Similar to previous best
    ## Run 303 stress 0.08352634 
    ## ... Procrustes: rmse 2.240005e-06  max resid 4.963183e-06 
    ## ... Similar to previous best
    ## Run 304 stress 0.1134014 
    ## Run 305 stress 0.08352634 
    ## ... Procrustes: rmse 3.742695e-06  max resid 6.203526e-06 
    ## ... Similar to previous best
    ## Run 306 stress 0.1279309 
    ## Run 307 stress 0.08352634 
    ## ... Procrustes: rmse 3.249575e-06  max resid 7.087386e-06 
    ## ... Similar to previous best
    ## Run 308 stress 0.08352634 
    ## ... Procrustes: rmse 1.786462e-06  max resid 3.705906e-06 
    ## ... Similar to previous best
    ## Run 309 stress 0.1134014 
    ## Run 310 stress 0.08352634 
    ## ... Procrustes: rmse 6.762477e-07  max resid 8.800622e-07 
    ## ... Similar to previous best
    ## Run 311 stress 0.08352634 
    ## ... Procrustes: rmse 5.036613e-06  max resid 1.138696e-05 
    ## ... Similar to previous best
    ## Run 312 stress 0.1134014 
    ## Run 313 stress 0.08352634 
    ## ... Procrustes: rmse 3.850405e-06  max resid 8.682192e-06 
    ## ... Similar to previous best
    ## Run 314 stress 0.1279309 
    ## Run 315 stress 0.08352634 
    ## ... Procrustes: rmse 2.415029e-06  max resid 5.312528e-06 
    ## ... Similar to previous best
    ## Run 316 stress 0.08352634 
    ## ... Procrustes: rmse 8.738948e-07  max resid 1.921377e-06 
    ## ... Similar to previous best
    ## Run 317 stress 0.1134014 
    ## Run 318 stress 0.1279309 
    ## Run 319 stress 0.127931 
    ## Run 320 stress 0.1989487 
    ## Run 321 stress 0.08352634 
    ## ... Procrustes: rmse 2.396532e-06  max resid 5.380149e-06 
    ## ... Similar to previous best
    ## Run 322 stress 0.1134015 
    ## Run 323 stress 0.1134015 
    ## Run 324 stress 0.1134015 
    ## Run 325 stress 0.08352634 
    ## ... Procrustes: rmse 1.023171e-06  max resid 1.83124e-06 
    ## ... Similar to previous best
    ## Run 326 stress 0.1134014 
    ## Run 327 stress 0.08352634 
    ## ... Procrustes: rmse 1.225384e-06  max resid 2.759903e-06 
    ## ... Similar to previous best
    ## Run 328 stress 0.08352634 
    ## ... Procrustes: rmse 1.0905e-06  max resid 2.233377e-06 
    ## ... Similar to previous best
    ## Run 329 stress 0.1279309 
    ## Run 330 stress 0.08352634 
    ## ... Procrustes: rmse 4.96335e-07  max resid 7.548503e-07 
    ## ... Similar to previous best
    ## Run 331 stress 0.08352634 
    ## ... Procrustes: rmse 3.377876e-06  max resid 7.565202e-06 
    ## ... Similar to previous best
    ## Run 332 stress 0.08352634 
    ## ... Procrustes: rmse 4.788016e-07  max resid 1.026692e-06 
    ## ... Similar to previous best
    ## Run 333 stress 0.08352634 
    ## ... Procrustes: rmse 1.361955e-06  max resid 3.023214e-06 
    ## ... Similar to previous best
    ## Run 334 stress 0.1279309 
    ## Run 335 stress 0.08352634 
    ## ... Procrustes: rmse 3.574026e-06  max resid 8.00592e-06 
    ## ... Similar to previous best
    ## Run 336 stress 0.08352634 
    ## ... Procrustes: rmse 1.623555e-06  max resid 3.601718e-06 
    ## ... Similar to previous best
    ## Run 337 stress 0.1134015 
    ## Run 338 stress 0.08352634 
    ## ... Procrustes: rmse 2.347282e-06  max resid 4.928481e-06 
    ## ... Similar to previous best
    ## Run 339 stress 0.08352634 
    ## ... Procrustes: rmse 3.475814e-06  max resid 7.701208e-06 
    ## ... Similar to previous best
    ## Run 340 stress 0.08352634 
    ## ... Procrustes: rmse 1.230792e-06  max resid 2.341524e-06 
    ## ... Similar to previous best
    ## Run 341 stress 0.08352634 
    ## ... Procrustes: rmse 2.129629e-06  max resid 4.381443e-06 
    ## ... Similar to previous best
    ## Run 342 stress 0.1134014 
    ## Run 343 stress 0.08352634 
    ## ... Procrustes: rmse 2.584025e-06  max resid 5.94363e-06 
    ## ... Similar to previous best
    ## Run 344 stress 0.08352634 
    ## ... Procrustes: rmse 1.393992e-06  max resid 2.957616e-06 
    ## ... Similar to previous best
    ## Run 345 stress 0.1279309 
    ## Run 346 stress 0.1989487 
    ## Run 347 stress 0.08352634 
    ## ... Procrustes: rmse 2.73152e-06  max resid 4.583959e-06 
    ## ... Similar to previous best
    ## Run 348 stress 0.08352634 
    ## ... Procrustes: rmse 2.19534e-06  max resid 4.8024e-06 
    ## ... Similar to previous best
    ## Run 349 stress 0.08352634 
    ## ... Procrustes: rmse 2.447488e-06  max resid 5.926996e-06 
    ## ... Similar to previous best
    ## Run 350 stress 0.1989155 
    ## Run 351 stress 0.1134014 
    ## Run 352 stress 0.08352634 
    ## ... Procrustes: rmse 5.464389e-06  max resid 1.195899e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.08352634 
    ## ... Procrustes: rmse 1.171967e-06  max resid 2.542255e-06 
    ## ... Similar to previous best
    ## Run 354 stress 0.1134014 
    ## Run 355 stress 0.08352634 
    ## ... Procrustes: rmse 1.007505e-06  max resid 1.887153e-06 
    ## ... Similar to previous best
    ## Run 356 stress 0.1134015 
    ## Run 357 stress 0.1279309 
    ## Run 358 stress 0.08352634 
    ## ... Procrustes: rmse 3.014307e-06  max resid 6.728691e-06 
    ## ... Similar to previous best
    ## Run 359 stress 0.08352634 
    ## ... Procrustes: rmse 3.539833e-06  max resid 7.946672e-06 
    ## ... Similar to previous best
    ## Run 360 stress 0.08352634 
    ## ... Procrustes: rmse 1.217294e-06  max resid 2.480944e-06 
    ## ... Similar to previous best
    ## Run 361 stress 0.1134014 
    ## Run 362 stress 0.08352634 
    ## ... Procrustes: rmse 2.988224e-06  max resid 6.377543e-06 
    ## ... Similar to previous best
    ## Run 363 stress 0.1134015 
    ## Run 364 stress 0.08352634 
    ## ... Procrustes: rmse 8.882256e-07  max resid 1.817725e-06 
    ## ... Similar to previous best
    ## Run 365 stress 0.1134015 
    ## Run 366 stress 0.1134015 
    ## Run 367 stress 0.08352634 
    ## ... Procrustes: rmse 1.275682e-06  max resid 2.833052e-06 
    ## ... Similar to previous best
    ## Run 368 stress 0.08352634 
    ## ... Procrustes: rmse 1.150034e-06  max resid 2.402611e-06 
    ## ... Similar to previous best
    ## Run 369 stress 0.08352634 
    ## ... Procrustes: rmse 2.129302e-06  max resid 4.961155e-06 
    ## ... Similar to previous best
    ## Run 370 stress 0.08352634 
    ## ... Procrustes: rmse 5.876048e-07  max resid 9.627146e-07 
    ## ... Similar to previous best
    ## Run 371 stress 0.08352634 
    ## ... Procrustes: rmse 1.592689e-06  max resid 3.103101e-06 
    ## ... Similar to previous best
    ## Run 372 stress 0.1134015 
    ## Run 373 stress 0.08352634 
    ## ... Procrustes: rmse 1.947227e-06  max resid 4.378908e-06 
    ## ... Similar to previous best
    ## Run 374 stress 0.08352634 
    ## ... Procrustes: rmse 1.884864e-06  max resid 4.275624e-06 
    ## ... Similar to previous best
    ## Run 375 stress 0.08352634 
    ## ... Procrustes: rmse 2.068097e-06  max resid 4.553453e-06 
    ## ... Similar to previous best
    ## Run 376 stress 0.1279309 
    ## Run 377 stress 0.1134014 
    ## Run 378 stress 0.1279309 
    ## Run 379 stress 0.08352634 
    ## ... Procrustes: rmse 9.314336e-07  max resid 2.019819e-06 
    ## ... Similar to previous best
    ## Run 380 stress 0.08352634 
    ## ... Procrustes: rmse 2.416825e-06  max resid 5.468526e-06 
    ## ... Similar to previous best
    ## Run 381 stress 0.08352634 
    ## ... Procrustes: rmse 3.026674e-06  max resid 6.72937e-06 
    ## ... Similar to previous best
    ## Run 382 stress 0.1134015 
    ## Run 383 stress 0.08352634 
    ## ... Procrustes: rmse 1.170399e-06  max resid 2.528322e-06 
    ## ... Similar to previous best
    ## Run 384 stress 0.08352634 
    ## ... Procrustes: rmse 3.904663e-07  max resid 6.458451e-07 
    ## ... Similar to previous best
    ## Run 385 stress 0.1989156 
    ## Run 386 stress 0.1279309 
    ## Run 387 stress 0.08352634 
    ## ... Procrustes: rmse 4.435354e-06  max resid 1.002094e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.2050061 
    ## Run 389 stress 0.1279309 
    ## Run 390 stress 0.08352634 
    ## ... Procrustes: rmse 1.041464e-06  max resid 2.338867e-06 
    ## ... Similar to previous best
    ## Run 391 stress 0.08352635 
    ## ... Procrustes: rmse 8.793232e-06  max resid 1.401621e-05 
    ## ... Similar to previous best
    ## Run 392 stress 0.08352634 
    ## ... Procrustes: rmse 3.022883e-06  max resid 6.611023e-06 
    ## ... Similar to previous best
    ## Run 393 stress 0.08352634 
    ## ... Procrustes: rmse 3.707826e-06  max resid 8.280247e-06 
    ## ... Similar to previous best
    ## Run 394 stress 0.08352634 
    ## ... Procrustes: rmse 2.143442e-06  max resid 4.842588e-06 
    ## ... Similar to previous best
    ## Run 395 stress 0.08352634 
    ## ... Procrustes: rmse 8.363135e-07  max resid 1.547336e-06 
    ## ... Similar to previous best
    ## Run 396 stress 0.1134015 
    ## Run 397 stress 0.08352634 
    ## ... Procrustes: rmse 5.092441e-06  max resid 1.137855e-05 
    ## ... Similar to previous best
    ## Run 398 stress 0.08352634 
    ## ... Procrustes: rmse 1.427319e-06  max resid 2.861625e-06 
    ## ... Similar to previous best
    ## Run 399 stress 0.08352634 
    ## ... Procrustes: rmse 3.066893e-06  max resid 6.757584e-06 
    ## ... Similar to previous best
    ## Run 400 stress 0.08352634 
    ## ... Procrustes: rmse 1.777939e-06  max resid 3.920196e-06 
    ## ... Similar to previous best
    ## Run 401 stress 0.08352634 
    ## ... Procrustes: rmse 5.728246e-06  max resid 1.268346e-05 
    ## ... Similar to previous best
    ## Run 402 stress 0.1134015 
    ## Run 403 stress 0.1279309 
    ## Run 404 stress 0.1134014 
    ## Run 405 stress 0.08352634 
    ## ... Procrustes: rmse 6.701909e-07  max resid 1.274429e-06 
    ## ... Similar to previous best
    ## Run 406 stress 0.1134015 
    ## Run 407 stress 0.08352634 
    ## ... Procrustes: rmse 1.551016e-06  max resid 2.859351e-06 
    ## ... Similar to previous best
    ## Run 408 stress 0.1279309 
    ## Run 409 stress 0.08352634 
    ## ... Procrustes: rmse 2.443914e-06  max resid 5.321984e-06 
    ## ... Similar to previous best
    ## Run 410 stress 0.1134015 
    ## Run 411 stress 0.1134014 
    ## Run 412 stress 0.08352634 
    ## ... Procrustes: rmse 4.237531e-06  max resid 9.620057e-06 
    ## ... Similar to previous best
    ## Run 413 stress 0.1134014 
    ## Run 414 stress 0.1279309 
    ## Run 415 stress 0.1279309 
    ## Run 416 stress 0.08352634 
    ## ... Procrustes: rmse 6.081161e-07  max resid 1.187981e-06 
    ## ... Similar to previous best
    ## Run 417 stress 0.08352634 
    ## ... Procrustes: rmse 3.627565e-06  max resid 8.043936e-06 
    ## ... Similar to previous best
    ## Run 418 stress 0.08352634 
    ## ... Procrustes: rmse 3.793333e-06  max resid 8.347462e-06 
    ## ... Similar to previous best
    ## Run 419 stress 0.1134014 
    ## Run 420 stress 0.08352634 
    ## ... Procrustes: rmse 4.497483e-06  max resid 1.010265e-05 
    ## ... Similar to previous best
    ## Run 421 stress 0.08352634 
    ## ... Procrustes: rmse 2.538471e-06  max resid 5.695518e-06 
    ## ... Similar to previous best
    ## Run 422 stress 0.08352634 
    ## ... Procrustes: rmse 6.61307e-07  max resid 1.268944e-06 
    ## ... Similar to previous best
    ## Run 423 stress 0.1279309 
    ## Run 424 stress 0.1279309 
    ## Run 425 stress 0.08352634 
    ## ... Procrustes: rmse 8.856424e-07  max resid 1.924832e-06 
    ## ... Similar to previous best
    ## Run 426 stress 0.08352634 
    ## ... Procrustes: rmse 9.292639e-07  max resid 1.772276e-06 
    ## ... Similar to previous best
    ## Run 427 stress 0.08352634 
    ## ... Procrustes: rmse 2.384689e-06  max resid 4.131107e-06 
    ## ... Similar to previous best
    ## Run 428 stress 0.08352634 
    ## ... Procrustes: rmse 2.404399e-06  max resid 5.317406e-06 
    ## ... Similar to previous best
    ## Run 429 stress 0.1279309 
    ## Run 430 stress 0.1989485 
    ## Run 431 stress 0.08352634 
    ## ... Procrustes: rmse 1.437451e-05  max resid 3.079865e-05 
    ## ... Similar to previous best
    ## Run 432 stress 0.08352634 
    ## ... Procrustes: rmse 3.729539e-06  max resid 8.036529e-06 
    ## ... Similar to previous best
    ## Run 433 stress 0.08352634 
    ## ... Procrustes: rmse 1.97627e-06  max resid 4.397179e-06 
    ## ... Similar to previous best
    ## Run 434 stress 0.08352634 
    ## ... Procrustes: rmse 1.07899e-06  max resid 2.659474e-06 
    ## ... Similar to previous best
    ## Run 435 stress 0.08352634 
    ## ... Procrustes: rmse 2.729813e-06  max resid 5.979175e-06 
    ## ... Similar to previous best
    ## Run 436 stress 0.1134014 
    ## Run 437 stress 0.1279309 
    ## Run 438 stress 0.1279309 
    ## Run 439 stress 0.2050052 
    ## Run 440 stress 0.1279309 
    ## Run 441 stress 0.08352634 
    ## ... Procrustes: rmse 2.179426e-06  max resid 4.639061e-06 
    ## ... Similar to previous best
    ## Run 442 stress 0.1134014 
    ## Run 443 stress 0.1989484 
    ## Run 444 stress 0.08352634 
    ## ... Procrustes: rmse 1.680093e-06  max resid 3.783942e-06 
    ## ... Similar to previous best
    ## Run 445 stress 0.08352634 
    ## ... Procrustes: rmse 2.41206e-06  max resid 5.414025e-06 
    ## ... Similar to previous best
    ## Run 446 stress 0.1134014 
    ## Run 447 stress 0.08352634 
    ## ... Procrustes: rmse 1.574627e-06  max resid 3.48636e-06 
    ## ... Similar to previous best
    ## Run 448 stress 0.1279309 
    ## Run 449 stress 0.1134015 
    ## Run 450 stress 0.08352634 
    ## ... Procrustes: rmse 2.589783e-06  max resid 5.549033e-06 
    ## ... Similar to previous best
    ## Run 451 stress 0.08352634 
    ## ... Procrustes: rmse 3.140191e-07  max resid 4.492853e-07 
    ## ... Similar to previous best
    ## Run 452 stress 0.08352634 
    ## ... Procrustes: rmse 8.828802e-07  max resid 1.486511e-06 
    ## ... Similar to previous best
    ## Run 453 stress 0.08352634 
    ## ... Procrustes: rmse 1.265666e-06  max resid 2.43412e-06 
    ## ... Similar to previous best
    ## Run 454 stress 0.08352634 
    ## ... Procrustes: rmse 4.763429e-06  max resid 1.05959e-05 
    ## ... Similar to previous best
    ## Run 455 stress 0.08352634 
    ## ... Procrustes: rmse 9.622924e-07  max resid 2.035808e-06 
    ## ... Similar to previous best
    ## Run 456 stress 0.08352634 
    ## ... Procrustes: rmse 6.01491e-07  max resid 1.32535e-06 
    ## ... Similar to previous best
    ## Run 457 stress 0.1134014 
    ## Run 458 stress 0.1279309 
    ## Run 459 stress 0.1279309 
    ## Run 460 stress 0.08352634 
    ## ... Procrustes: rmse 1.321599e-06  max resid 2.424942e-06 
    ## ... Similar to previous best
    ## Run 461 stress 0.08352634 
    ## ... Procrustes: rmse 6.52797e-07  max resid 1.340525e-06 
    ## ... Similar to previous best
    ## Run 462 stress 0.08352634 
    ## ... Procrustes: rmse 2.519582e-06  max resid 5.723206e-06 
    ## ... Similar to previous best
    ## Run 463 stress 0.1279309 
    ## Run 464 stress 0.1134014 
    ## Run 465 stress 0.08352634 
    ## ... Procrustes: rmse 1.141923e-06  max resid 2.288348e-06 
    ## ... Similar to previous best
    ## Run 466 stress 0.1989156 
    ## Run 467 stress 0.08352634 
    ## ... Procrustes: rmse 1.004927e-06  max resid 2.198289e-06 
    ## ... Similar to previous best
    ## Run 468 stress 0.1134015 
    ## Run 469 stress 0.08352634 
    ## ... Procrustes: rmse 1.989134e-06  max resid 4.44555e-06 
    ## ... Similar to previous best
    ## Run 470 stress 0.1279309 
    ## Run 471 stress 0.1134014 
    ## Run 472 stress 0.08352634 
    ## ... Procrustes: rmse 1.19463e-06  max resid 2.715001e-06 
    ## ... Similar to previous best
    ## Run 473 stress 0.1134014 
    ## Run 474 stress 0.08352634 
    ## ... Procrustes: rmse 4.064117e-06  max resid 9.213573e-06 
    ## ... Similar to previous best
    ## Run 475 stress 0.08352634 
    ## ... Procrustes: rmse 1.733124e-06  max resid 3.734493e-06 
    ## ... Similar to previous best
    ## Run 476 stress 0.127931 
    ## Run 477 stress 0.1134015 
    ## Run 478 stress 0.1134014 
    ## Run 479 stress 0.1279309 
    ## Run 480 stress 0.1279309 
    ## Run 481 stress 0.08352634 
    ## ... Procrustes: rmse 3.092011e-06  max resid 6.912763e-06 
    ## ... Similar to previous best
    ## Run 482 stress 0.1134015 
    ## Run 483 stress 0.08352634 
    ## ... Procrustes: rmse 3.139577e-06  max resid 7.070095e-06 
    ## ... Similar to previous best
    ## Run 484 stress 0.1134014 
    ## Run 485 stress 0.08352634 
    ## ... Procrustes: rmse 1.182973e-06  max resid 2.493976e-06 
    ## ... Similar to previous best
    ## Run 486 stress 0.1279309 
    ## Run 487 stress 0.1989485 
    ## Run 488 stress 0.3084871 
    ## Run 489 stress 0.08352634 
    ## ... Procrustes: rmse 1.177055e-06  max resid 2.628765e-06 
    ## ... Similar to previous best
    ## Run 490 stress 0.08352634 
    ## ... Procrustes: rmse 4.229928e-06  max resid 9.425608e-06 
    ## ... Similar to previous best
    ## Run 491 stress 0.08352634 
    ## ... Procrustes: rmse 1.752467e-06  max resid 3.922403e-06 
    ## ... Similar to previous best
    ## Run 492 stress 0.1989154 
    ## Run 493 stress 0.08352634 
    ## ... Procrustes: rmse 1.591909e-06  max resid 3.07551e-06 
    ## ... Similar to previous best
    ## Run 494 stress 0.08352634 
    ## ... Procrustes: rmse 1.976112e-06  max resid 4.278494e-06 
    ## ... Similar to previous best
    ## Run 495 stress 0.08352634 
    ## ... Procrustes: rmse 2.96479e-06  max resid 6.650306e-06 
    ## ... Similar to previous best
    ## Run 496 stress 0.08352634 
    ## ... Procrustes: rmse 3.402081e-06  max resid 7.516925e-06 
    ## ... Similar to previous best
    ## Run 497 stress 0.1134015 
    ## Run 498 stress 0.1134014 
    ## Run 499 stress 0.08352634 
    ## ... Procrustes: rmse 1.348767e-06  max resid 3.08856e-06 
    ## ... Similar to previous best
    ## Run 500 stress 0.1134015 
    ## *** Best solution repeated 250 times

``` r
# Stratified lakes and ocean sites
PD_beta_env_SO_NMDS <- metaMDS(PD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.92006e-05 
    ## Run 1 stress 9.956341e-05 
    ## ... Procrustes: rmse 0.0003181812  max resid 0.0007309883 
    ## ... Similar to previous best
    ## Run 2 stress 9.901888e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000276791  max resid 0.0007085606 
    ## ... Similar to previous best
    ## Run 3 stress 9.238321e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003194426  max resid 0.0006628984 
    ## ... Similar to previous best
    ## Run 4 stress 0.000728062 
    ## Run 5 stress 9.978078e-05 
    ## ... Procrustes: rmse 0.0003060536  max resid 0.0007063616 
    ## ... Similar to previous best
    ## Run 6 stress 9.880904e-05 
    ## ... Procrustes: rmse 0.0003041332  max resid 0.0007018839 
    ## ... Similar to previous best
    ## Run 7 stress 9.760439e-05 
    ## ... Procrustes: rmse 0.0003023038  max resid 0.0006980452 
    ## ... Similar to previous best
    ## Run 8 stress 9.615273e-05 
    ## ... Procrustes: rmse 0.0002435384  max resid 0.0005917881 
    ## ... Similar to previous best
    ## Run 9 stress 9.776918e-05 
    ## ... Procrustes: rmse 0.0003978603  max resid 0.0005913997 
    ## ... Similar to previous best
    ## Run 10 stress 0.0005012829 
    ## ... Procrustes: rmse 0.004131279  max resid 0.006516341 
    ## ... Similar to previous best
    ## Run 11 stress 9.972084e-05 
    ## ... Procrustes: rmse 3.3972e-05  max resid 5.250801e-05 
    ## ... Similar to previous best
    ## Run 12 stress 9.658441e-05 
    ## ... Procrustes: rmse 0.0002994403  max resid 0.0006955395 
    ## ... Similar to previous best
    ## Run 13 stress 9.829055e-05 
    ## ... Procrustes: rmse 2.410912e-05  max resid 3.893499e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.0006207721 
    ## Run 15 stress 9.910108e-05 
    ## ... Procrustes: rmse 0.0002808727  max resid 0.0005317108 
    ## ... Similar to previous best
    ## Run 16 stress 9.926339e-05 
    ## ... Procrustes: rmse 0.0002488048  max resid 0.0006042433 
    ## ... Similar to previous best
    ## Run 17 stress 0.0006008358 
    ## Run 18 stress 9.918813e-05 
    ## ... Procrustes: rmse 0.0002673445  max resid 0.000699041 
    ## ... Similar to previous best
    ## Run 19 stress 9.80503e-05 
    ## ... Procrustes: rmse 0.0002671015  max resid 0.000700039 
    ## ... Similar to previous best
    ## Run 20 stress 9.846588e-05 
    ## ... Procrustes: rmse 0.0002066209  max resid 0.000585387 
    ## ... Similar to previous best
    ## Run 21 stress 8.906376e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004348295  max resid 0.0007736746 
    ## ... Similar to previous best
    ## Run 22 stress 9.996035e-05 
    ## ... Procrustes: rmse 0.0004679627  max resid 0.0008270971 
    ## ... Similar to previous best
    ## Run 23 stress 9.816065e-05 
    ## ... Procrustes: rmse 0.0004602411  max resid 0.000815426 
    ## ... Similar to previous best
    ## Run 24 stress 0.001030871 
    ## Run 25 stress 0.001161657 
    ## Run 26 stress 9.865286e-05 
    ## ... Procrustes: rmse 0.0004225266  max resid 0.0006658678 
    ## ... Similar to previous best
    ## Run 27 stress 0.001181308 
    ## Run 28 stress 0.00136176 
    ## Run 29 stress 9.808434e-05 
    ## ... Procrustes: rmse 0.0004396094  max resid 0.0008147996 
    ## ... Similar to previous best
    ## Run 30 stress 0.001117907 
    ## Run 31 stress 0.0003261331 
    ## ... Procrustes: rmse 0.00285939  max resid 0.004241667 
    ## ... Similar to previous best
    ## Run 32 stress 9.92372e-05 
    ## ... Procrustes: rmse 0.0004432159  max resid 0.0008192631 
    ## ... Similar to previous best
    ## Run 33 stress 9.940378e-05 
    ## ... Procrustes: rmse 0.0004635316  max resid 0.0007788513 
    ## ... Similar to previous best
    ## Run 34 stress 9.799995e-05 
    ## ... Procrustes: rmse 0.0004578211  max resid 0.000812567 
    ## ... Similar to previous best
    ## Run 35 stress 0.001174344 
    ## Run 36 stress 0.000465663 
    ## ... Procrustes: rmse 0.004090051  max resid 0.006496424 
    ## ... Similar to previous best
    ## Run 37 stress 9.868001e-05 
    ## ... Procrustes: rmse 0.0004627503  max resid 0.0008197661 
    ## ... Similar to previous best
    ## Run 38 stress 9.863633e-05 
    ## ... Procrustes: rmse 0.0004584698  max resid 0.0008108693 
    ## ... Similar to previous best
    ## Run 39 stress 0.0007200112 
    ## Run 40 stress 9.670913e-05 
    ## ... Procrustes: rmse 0.0004550125  max resid 0.0008083265 
    ## ... Similar to previous best
    ## Run 41 stress 9.870264e-05 
    ## ... Procrustes: rmse 0.000399994  max resid 0.000654674 
    ## ... Similar to previous best
    ## Run 42 stress 9.841809e-05 
    ## ... Procrustes: rmse 0.0005193724  max resid 0.0008156275 
    ## ... Similar to previous best
    ## Run 43 stress 9.86766e-05 
    ## ... Procrustes: rmse 0.00039817  max resid 0.0006547638 
    ## ... Similar to previous best
    ## Run 44 stress 0.0007808466 
    ## Run 45 stress 9.943896e-05 
    ## ... Procrustes: rmse 0.0003974663  max resid 0.0006590058 
    ## ... Similar to previous best
    ## Run 46 stress 0.001300362 
    ## Run 47 stress 9.916687e-05 
    ## ... Procrustes: rmse 0.0004253238  max resid 0.0006684993 
    ## ... Similar to previous best
    ## Run 48 stress 9.922574e-05 
    ## ... Procrustes: rmse 0.0003924877  max resid 0.0006346183 
    ## ... Similar to previous best
    ## Run 49 stress 9.929373e-05 
    ## ... Procrustes: rmse 0.000441888  max resid 0.0008218127 
    ## ... Similar to previous best
    ## Run 50 stress 0.0006940575 
    ## Run 51 stress 0.0004952876 
    ## ... Procrustes: rmse 0.004350722  max resid 0.006903915 
    ## ... Similar to previous best
    ## Run 52 stress 0.0003466657 
    ## ... Procrustes: rmse 0.003052126  max resid 0.004894698 
    ## ... Similar to previous best
    ## Run 53 stress 9.937822e-05 
    ## ... Procrustes: rmse 0.0004022272  max resid 0.0006570627 
    ## ... Similar to previous best
    ## Run 54 stress 9.916593e-05 
    ## ... Procrustes: rmse 0.0004615353  max resid 0.0008132 
    ## ... Similar to previous best
    ## Run 55 stress 0.0005602291 
    ## ... Procrustes: rmse 0.004916196  max resid 0.00778421 
    ## ... Similar to previous best
    ## Run 56 stress 0.0007998259 
    ## Run 57 stress 9.825855e-05 
    ## ... Procrustes: rmse 0.0004216349  max resid 0.0006639907 
    ## ... Similar to previous best
    ## Run 58 stress 9.996249e-05 
    ## ... Procrustes: rmse 0.0008760706  max resid 0.001380217 
    ## ... Similar to previous best
    ## Run 59 stress 0.0006649613 
    ## Run 60 stress 9.712273e-05 
    ## ... Procrustes: rmse 0.0004770085  max resid 0.0007127774 
    ## ... Similar to previous best
    ## Run 61 stress 0.001181871 
    ## Run 62 stress 9.896227e-05 
    ## ... Procrustes: rmse 0.000470666  max resid 0.0007123889 
    ## ... Similar to previous best
    ## Run 63 stress 9.833252e-05 
    ## ... Procrustes: rmse 0.0003977661  max resid 0.0006472954 
    ## ... Similar to previous best
    ## Run 64 stress 0.0005084739 
    ## ... Procrustes: rmse 0.004462975  max resid 0.00706734 
    ## ... Similar to previous best
    ## Run 65 stress 9.78117e-05 
    ## ... Procrustes: rmse 0.0004785669  max resid 0.0007021689 
    ## ... Similar to previous best
    ## Run 66 stress 9.972692e-05 
    ## ... Procrustes: rmse 0.0004259237  max resid 0.0006682403 
    ## ... Similar to previous best
    ## Run 67 stress 9.931104e-05 
    ## ... Procrustes: rmse 0.0004614528  max resid 0.0008143519 
    ## ... Similar to previous best
    ## Run 68 stress 9.865632e-05 
    ## ... Procrustes: rmse 0.0004234287  max resid 0.0006672791 
    ## ... Similar to previous best
    ## Run 69 stress 9.264469e-05 
    ## ... Procrustes: rmse 0.0004127619  max resid 0.0007726822 
    ## ... Similar to previous best
    ## Run 70 stress 9.931853e-05 
    ## ... Procrustes: rmse 0.0004637635  max resid 0.0008175679 
    ## ... Similar to previous best
    ## Run 71 stress 0.0009878358 
    ## Run 72 stress 9.917355e-05 
    ## ... Procrustes: rmse 0.0004905835  max resid 0.0007260842 
    ## ... Similar to previous best
    ## Run 73 stress 9.984284e-05 
    ## ... Procrustes: rmse 0.0004153685  max resid 0.0006618661 
    ## ... Similar to previous best
    ## Run 74 stress 9.926118e-05 
    ## ... Procrustes: rmse 0.0001462093  max resid 0.0003984022 
    ## ... Similar to previous best
    ## Run 75 stress 0.0003992457 
    ## ... Procrustes: rmse 0.003503979  max resid 0.005240949 
    ## ... Similar to previous best
    ## Run 76 stress 9.974899e-05 
    ## ... Procrustes: rmse 0.000445868  max resid 0.0008246191 
    ## ... Similar to previous best
    ## Run 77 stress 0.0006607034 
    ## Run 78 stress 9.661978e-05 
    ## ... Procrustes: rmse 0.000423419  max resid 0.0006368922 
    ## ... Similar to previous best
    ## Run 79 stress 0.3275204 
    ## Run 80 stress 0.0005463365 
    ## ... Procrustes: rmse 0.004804057  max resid 0.007263532 
    ## ... Similar to previous best
    ## Run 81 stress 9.947596e-05 
    ## ... Procrustes: rmse 0.0004624871  max resid 0.0007743094 
    ## ... Similar to previous best
    ## Run 82 stress 0.000588128 
    ## ... Procrustes: rmse 0.005177538  max resid 0.0078275 
    ## Run 83 stress 9.744509e-05 
    ## ... Procrustes: rmse 0.0004561639  max resid 0.0008108455 
    ## ... Similar to previous best
    ## Run 84 stress 9.895488e-05 
    ## ... Procrustes: rmse 0.0004336298  max resid 0.0006552932 
    ## ... Similar to previous best
    ## Run 85 stress 0.0009744089 
    ## Run 86 stress 9.82321e-05 
    ## ... Procrustes: rmse 0.0003980894  max resid 0.0006522867 
    ## ... Similar to previous best
    ## Run 87 stress 9.94899e-05 
    ## ... Procrustes: rmse 0.0004834759  max resid 0.0007169223 
    ## ... Similar to previous best
    ## Run 88 stress 9.978497e-05 
    ## ... Procrustes: rmse 0.000529532  max resid 0.0008246984 
    ## ... Similar to previous best
    ## Run 89 stress 9.966124e-05 
    ## ... Procrustes: rmse 0.0008734698  max resid 0.001376868 
    ## ... Similar to previous best
    ## Run 90 stress 0.001019094 
    ## Run 91 stress 0.0007834982 
    ## Run 92 stress 9.561108e-05 
    ## ... Procrustes: rmse 0.0003756289  max resid 0.0006254085 
    ## ... Similar to previous best
    ## Run 93 stress 0.000915342 
    ## Run 94 stress 0.001099316 
    ## Run 95 stress 9.923359e-05 
    ## ... Procrustes: rmse 0.000462031  max resid 0.0008166255 
    ## ... Similar to previous best
    ## Run 96 stress 0.0009342591 
    ## Run 97 stress 9.939899e-05 
    ## ... Procrustes: rmse 0.0004625824  max resid 0.0007751982 
    ## ... Similar to previous best
    ## Run 98 stress 9.971159e-05 
    ## ... Procrustes: rmse 0.0004671207  max resid 0.0008249166 
    ## ... Similar to previous best
    ## Run 99 stress 9.809137e-05 
    ## ... Procrustes: rmse 0.0004671968  max resid 0.000708743 
    ## ... Similar to previous best
    ## Run 100 stress 0.0008197625 
    ## Run 101 stress 0.0008028585 
    ## Run 102 stress 9.916579e-05 
    ## ... Procrustes: rmse 0.0004865526  max resid 0.0007117679 
    ## ... Similar to previous best
    ## Run 103 stress 9.977359e-05 
    ## ... Procrustes: rmse 0.000391963  max resid 0.0006522168 
    ## ... Similar to previous best
    ## Run 104 stress 0.0007760765 
    ## Run 105 stress 0.0007700975 
    ## Run 106 stress 0.001142638 
    ## Run 107 stress 9.814835e-05 
    ## ... Procrustes: rmse 0.0004215495  max resid 0.000663359 
    ## ... Similar to previous best
    ## Run 108 stress 9.930378e-05 
    ## ... Procrustes: rmse 0.0004708552  max resid 0.000708754 
    ## ... Similar to previous best
    ## Run 109 stress 9.948216e-05 
    ## ... Procrustes: rmse 0.0004897411  max resid 0.0007180753 
    ## ... Similar to previous best
    ## Run 110 stress 9.865919e-05 
    ## ... Procrustes: rmse 0.0004233614  max resid 0.0006664469 
    ## ... Similar to previous best
    ## Run 111 stress 9.976385e-05 
    ## ... Procrustes: rmse 0.0004685898  max resid 0.000706189 
    ## ... Similar to previous best
    ## Run 112 stress 9.556285e-05 
    ## ... Procrustes: rmse 0.0003883966  max resid 0.0006380361 
    ## ... Similar to previous best
    ## Run 113 stress 0.001130889 
    ## Run 114 stress 0.0005472944 
    ## ... Procrustes: rmse 0.004804294  max resid 0.007607148 
    ## ... Similar to previous best
    ## Run 115 stress 9.762916e-05 
    ## ... Procrustes: rmse 0.0004739418  max resid 0.0006985692 
    ## ... Similar to previous best
    ## Run 116 stress 0.0005025723 
    ## ... Procrustes: rmse 0.004408771  max resid 0.006985426 
    ## ... Similar to previous best
    ## Run 117 stress 9.752122e-05 
    ## ... Procrustes: rmse 0.0004280184  max resid 0.0006462474 
    ## ... Similar to previous best
    ## Run 118 stress 0.001175415 
    ## Run 119 stress 9.378429e-05 
    ## ... Procrustes: rmse 0.0004196073  max resid 0.0007856925 
    ## ... Similar to previous best
    ## Run 120 stress 9.8716e-05 
    ## ... Procrustes: rmse 0.0004234717  max resid 0.0006659623 
    ## ... Similar to previous best
    ## Run 121 stress 9.88461e-05 
    ## ... Procrustes: rmse 0.0004005836  max resid 0.000654445 
    ## ... Similar to previous best
    ## Run 122 stress 9.953809e-05 
    ## ... Procrustes: rmse 0.000402781  max resid 0.000655842 
    ## ... Similar to previous best
    ## Run 123 stress 0.001343783 
    ## Run 124 stress 9.895901e-05 
    ## ... Procrustes: rmse 0.0004501348  max resid 0.0008025343 
    ## ... Similar to previous best
    ## Run 125 stress 9.957145e-05 
    ## ... Procrustes: rmse 0.0003949085  max resid 0.0006541018 
    ## ... Similar to previous best
    ## Run 126 stress 9.962571e-05 
    ## ... Procrustes: rmse 0.0004415782  max resid 0.0008220769 
    ## ... Similar to previous best
    ## Run 127 stress 9.726954e-05 
    ## ... Procrustes: rmse 0.0004573675  max resid 0.0008119446 
    ## ... Similar to previous best
    ## Run 128 stress 9.570988e-05 
    ## ... Procrustes: rmse 0.000366151  max resid 0.0005979 
    ## ... Similar to previous best
    ## Run 129 stress 9.998759e-05 
    ## ... Procrustes: rmse 0.0004655147  max resid 0.0007798813 
    ## ... Similar to previous best
    ## Run 130 stress 9.844291e-05 
    ## ... Procrustes: rmse 0.0004613001  max resid 0.0008193799 
    ## ... Similar to previous best
    ## Run 131 stress 0.3266727 
    ## Run 132 stress 9.985654e-05 
    ## ... Procrustes: rmse 0.0004739599  max resid 0.0007200791 
    ## ... Similar to previous best
    ## Run 133 stress 9.883451e-05 
    ## ... Procrustes: rmse 0.0004227942  max resid 0.0006641193 
    ## ... Similar to previous best
    ## Run 134 stress 9.9433e-05 
    ## ... Procrustes: rmse 0.0004250294  max resid 0.0006687402 
    ## ... Similar to previous best
    ## Run 135 stress 9.986561e-05 
    ## ... Procrustes: rmse 0.000428146  max resid 0.0006718059 
    ## ... Similar to previous best
    ## Run 136 stress 9.997469e-05 
    ## ... Procrustes: rmse 0.0004287988  max resid 0.0006734605 
    ## ... Similar to previous best
    ## Run 137 stress 0.00065746 
    ## Run 138 stress 8.613602e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001890321  max resid 0.0004714183 
    ## ... Similar to previous best
    ## Run 139 stress 9.974803e-05 
    ## ... Procrustes: rmse 0.0004916669  max resid 0.0008067042 
    ## ... Similar to previous best
    ## Run 140 stress 9.907926e-05 
    ## ... Procrustes: rmse 0.0004402394  max resid 0.0006777255 
    ## ... Similar to previous best
    ## Run 141 stress 9.993485e-05 
    ## ... Procrustes: rmse 0.0004265855  max resid 0.0007656669 
    ## ... Similar to previous best
    ## Run 142 stress 9.704797e-05 
    ## ... Procrustes: rmse 0.0004050938  max resid 0.0007301916 
    ## ... Similar to previous best
    ## Run 143 stress 9.998855e-05 
    ## ... Procrustes: rmse 0.000465935  max resid 0.0007937684 
    ## ... Similar to previous best
    ## Run 144 stress 9.992825e-05 
    ## ... Procrustes: rmse 0.0004031612  max resid 0.0006324746 
    ## ... Similar to previous best
    ## Run 145 stress 9.785639e-05 
    ## ... Procrustes: rmse 0.0004417626  max resid 0.0007402611 
    ## ... Similar to previous best
    ## Run 146 stress 9.94804e-05 
    ## ... Procrustes: rmse 0.0001340501  max resid 0.0002916371 
    ## ... Similar to previous best
    ## Run 147 stress 0.00112867 
    ## Run 148 stress 9.665425e-05 
    ## ... Procrustes: rmse 0.0003862572  max resid 0.0006016244 
    ## ... Similar to previous best
    ## Run 149 stress 9.411926e-05 
    ## ... Procrustes: rmse 0.0004211507  max resid 0.0007193051 
    ## ... Similar to previous best
    ## Run 150 stress 0.001250927 
    ## Run 151 stress 9.728824e-05 
    ## ... Procrustes: rmse 0.0004539291  max resid 0.0007781948 
    ## ... Similar to previous best
    ## Run 152 stress 9.891907e-05 
    ## ... Procrustes: rmse 0.0004446417  max resid 0.0007544591 
    ## ... Similar to previous best
    ## Run 153 stress 9.73004e-05 
    ## ... Procrustes: rmse 0.000446178  max resid 0.0007652187 
    ## ... Similar to previous best
    ## Run 154 stress 9.863853e-05 
    ## ... Procrustes: rmse 0.0004377645  max resid 0.0006750021 
    ## ... Similar to previous best
    ## Run 155 stress 0.0004565613 
    ## ... Procrustes: rmse 0.004004966  max resid 0.006199141 
    ## ... Similar to previous best
    ## Run 156 stress 0.0006927904 
    ## Run 157 stress 9.562724e-05 
    ## ... Procrustes: rmse 0.0004095401  max resid 0.0007399432 
    ## ... Similar to previous best
    ## Run 158 stress 9.889641e-05 
    ## ... Procrustes: rmse 0.000399299  max resid 0.000747247 
    ## ... Similar to previous best
    ## Run 159 stress 9.876764e-05 
    ## ... Procrustes: rmse 0.000457648  max resid 0.0007960072 
    ## ... Similar to previous best
    ## Run 160 stress 9.763751e-05 
    ## ... Procrustes: rmse 0.0004298289  max resid 0.000739424 
    ## ... Similar to previous best
    ## Run 161 stress 0.0005955541 
    ## Run 162 stress 9.94947e-05 
    ## ... Procrustes: rmse 0.0004594011  max resid 0.0007967227 
    ## ... Similar to previous best
    ## Run 163 stress 0.2570691 
    ## Run 164 stress 9.852726e-05 
    ## ... Procrustes: rmse 0.0004215087  max resid 0.0007560906 
    ## ... Similar to previous best
    ## Run 165 stress 9.920689e-05 
    ## ... Procrustes: rmse 0.0003990409  max resid 0.000625136 
    ## ... Similar to previous best
    ## Run 166 stress 9.844285e-05 
    ## ... Procrustes: rmse 0.0004194189  max resid 0.0007533727 
    ## ... Similar to previous best
    ## Run 167 stress 0.001242467 
    ## Run 168 stress 9.983436e-05 
    ## ... Procrustes: rmse 0.0004216261  max resid 0.0007572648 
    ## ... Similar to previous best
    ## Run 169 stress 9.89362e-05 
    ## ... Procrustes: rmse 0.0003958373  max resid 0.0006146932 
    ## ... Similar to previous best
    ## Run 170 stress 9.742266e-05 
    ## ... Procrustes: rmse 0.0004100428  max resid 0.000739394 
    ## ... Similar to previous best
    ## Run 171 stress 9.97579e-05 
    ## ... Procrustes: rmse 0.0004420942  max resid 0.0006758706 
    ## ... Similar to previous best
    ## Run 172 stress 0.0007107939 
    ## Run 173 stress 0.0004968999 
    ## ... Procrustes: rmse 0.004363779  max resid 0.006695902 
    ## ... Similar to previous best
    ## Run 174 stress 0.0005343961 
    ## ... Procrustes: rmse 0.004690071  max resid 0.007202169 
    ## ... Similar to previous best
    ## Run 175 stress 0.001045214 
    ## Run 176 stress 0.0007653798 
    ## Run 177 stress 9.526999e-05 
    ## ... Procrustes: rmse 0.0004276181  max resid 0.0007201028 
    ## ... Similar to previous best
    ## Run 178 stress 9.920561e-05 
    ## ... Procrustes: rmse 0.0004021515  max resid 0.0006293505 
    ## ... Similar to previous best
    ## Run 179 stress 0.0008189545 
    ## Run 180 stress 9.457852e-05 
    ## ... Procrustes: rmse 0.0004219187  max resid 0.000649919 
    ## ... Similar to previous best
    ## Run 181 stress 9.960521e-05 
    ## ... Procrustes: rmse 0.0004220527  max resid 0.0007645749 
    ## ... Similar to previous best
    ## Run 182 stress 0.0007829732 
    ## Run 183 stress 9.972829e-05 
    ## ... Procrustes: rmse 0.0004045736  max resid 0.0006326868 
    ## ... Similar to previous best
    ## Run 184 stress 0.0007884079 
    ## Run 185 stress 9.988603e-05 
    ## ... Procrustes: rmse 0.0004931392  max resid 0.0008085606 
    ## ... Similar to previous best
    ## Run 186 stress 0.0004556081 
    ## ... Procrustes: rmse 0.004006996  max resid 0.006132958 
    ## ... Similar to previous best
    ## Run 187 stress 9.823272e-05 
    ## ... Procrustes: rmse 0.0004788724  max resid 0.0007887659 
    ## ... Similar to previous best
    ## Run 188 stress 0.0003302244 
    ## ... Procrustes: rmse 0.002888465  max resid 0.004467522 
    ## ... Similar to previous best
    ## Run 189 stress 0.0004890764 
    ## ... Procrustes: rmse 0.004299904  max resid 0.006590803 
    ## ... Similar to previous best
    ## Run 190 stress 9.973673e-05 
    ## ... Procrustes: rmse 0.0004047638  max resid 0.0006323662 
    ## ... Similar to previous best
    ## Run 191 stress 9.549405e-05 
    ## ... Procrustes: rmse 0.0004319026  max resid 0.0007236291 
    ## ... Similar to previous best
    ## Run 192 stress 0.0008148307 
    ## Run 193 stress 9.953309e-05 
    ## ... Procrustes: rmse 0.0004673826  max resid 0.0007971546 
    ## ... Similar to previous best
    ## Run 194 stress 9.86707e-05 
    ## ... Procrustes: rmse 0.0004431905  max resid 0.0007521339 
    ## ... Similar to previous best
    ## Run 195 stress 0.0006840721 
    ## Run 196 stress 9.411251e-05 
    ## ... Procrustes: rmse 0.0004139494  max resid 0.000641313 
    ## ... Similar to previous best
    ## Run 197 stress 0.0005680102 
    ## ... Procrustes: rmse 0.004988637  max resid 0.00765739 
    ## ... Similar to previous best
    ## Run 198 stress 9.81326e-05 
    ## ... Procrustes: rmse 0.000418242  max resid 0.0007474822 
    ## ... Similar to previous best
    ## Run 199 stress 9.854972e-05 
    ## ... Procrustes: rmse 0.0004603424  max resid 0.0007850089 
    ## ... Similar to previous best
    ## Run 200 stress 0.0004182202 
    ## ... Procrustes: rmse 0.003680543  max resid 0.005629003 
    ## ... Similar to previous best
    ## Run 201 stress 0.0006952094 
    ## Run 202 stress 9.849094e-05 
    ## ... Procrustes: rmse 0.00043772  max resid 0.0006715371 
    ## ... Similar to previous best
    ## Run 203 stress 0.0006048988 
    ## Run 204 stress 9.697213e-05 
    ## ... Procrustes: rmse 0.0004793961  max resid 0.0007890694 
    ## ... Similar to previous best
    ## Run 205 stress 9.827046e-05 
    ## ... Procrustes: rmse 0.0004413287  max resid 0.0007504881 
    ## ... Similar to previous best
    ## Run 206 stress 0.0007263879 
    ## Run 207 stress 9.951525e-05 
    ## ... Procrustes: rmse 0.0004258874  max resid 0.0007634923 
    ## ... Similar to previous best
    ## Run 208 stress 9.830396e-05 
    ## ... Procrustes: rmse 0.0004517617  max resid 0.0007862487 
    ## ... Similar to previous best
    ## Run 209 stress 0.0007698892 
    ## Run 210 stress 9.912872e-05 
    ## ... Procrustes: rmse 0.0004446098  max resid 0.000754334 
    ## ... Similar to previous best
    ## Run 211 stress 9.8937e-05 
    ## ... Procrustes: rmse 0.0004930898  max resid 0.000796053 
    ## ... Similar to previous best
    ## Run 212 stress 9.959979e-05 
    ## ... Procrustes: rmse 0.0004068019  max resid 0.0007531881 
    ## ... Similar to previous best
    ## Run 213 stress 9.885795e-05 
    ## ... Procrustes: rmse 0.0004176729  max resid 0.0007545102 
    ## ... Similar to previous best
    ## Run 214 stress 9.320579e-05 
    ## ... Procrustes: rmse 0.0004197413  max resid 0.0007148172 
    ## ... Similar to previous best
    ## Run 215 stress 9.836205e-05 
    ## ... Procrustes: rmse 0.0004861731  max resid 0.0007990349 
    ## ... Similar to previous best
    ## Run 216 stress 9.981259e-05 
    ## ... Procrustes: rmse 0.0004919216  max resid 0.0008082621 
    ## ... Similar to previous best
    ## Run 217 stress 9.865903e-05 
    ## ... Procrustes: rmse 0.000400159  max resid 0.0006246664 
    ## ... Similar to previous best
    ## Run 218 stress 9.956172e-05 
    ## ... Procrustes: rmse 0.0004625273  max resid 0.0008030798 
    ## ... Similar to previous best
    ## Run 219 stress 9.684078e-05 
    ## ... Procrustes: rmse 0.0006417955  max resid 0.001002735 
    ## ... Similar to previous best
    ## Run 220 stress 0.001241101 
    ## Run 221 stress 0.0009487459 
    ## Run 222 stress 9.343097e-05 
    ## ... Procrustes: rmse 0.0003924619  max resid 0.0007214419 
    ## ... Similar to previous best
    ## Run 223 stress 9.85219e-05 
    ## ... Procrustes: rmse 0.0004578724  max resid 0.0007956319 
    ## ... Similar to previous best
    ## Run 224 stress 9.890317e-05 
    ## ... Procrustes: rmse 0.0004011824  max resid 0.0006261989 
    ## ... Similar to previous best
    ## Run 225 stress 9.669121e-05 
    ## ... Procrustes: rmse 0.0004147578  max resid 0.0007471315 
    ## ... Similar to previous best
    ## Run 226 stress 0.0005198296 
    ## ... Procrustes: rmse 0.004566899  max resid 0.006997955 
    ## ... Similar to previous best
    ## Run 227 stress 0.0005959245 
    ## Run 228 stress 9.802026e-05 
    ## ... Procrustes: rmse 0.000397614  max resid 0.0006216218 
    ## ... Similar to previous best
    ## Run 229 stress 9.945786e-05 
    ## ... Procrustes: rmse 0.0002028797  max resid 0.0004836829 
    ## ... Similar to previous best
    ## Run 230 stress 0.3407338 
    ## Run 231 stress 9.221907e-05 
    ## ... Procrustes: rmse 0.0001412036  max resid 0.000279239 
    ## ... Similar to previous best
    ## Run 232 stress 9.850996e-05 
    ## ... Procrustes: rmse 0.0004220555  max resid 0.0007583006 
    ## ... Similar to previous best
    ## Run 233 stress 9.883171e-05 
    ## ... Procrustes: rmse 0.0004602194  max resid 0.0007847419 
    ## ... Similar to previous best
    ## Run 234 stress 0.0008492356 
    ## Run 235 stress 0.001164807 
    ## Run 236 stress 9.567699e-05 
    ## ... Procrustes: rmse 0.0001176972  max resid 0.0003169472 
    ## ... Similar to previous best
    ## Run 237 stress 9.645915e-05 
    ## ... Procrustes: rmse 0.000391133  max resid 0.0006981199 
    ## ... Similar to previous best
    ## Run 238 stress 0.0004790511 
    ## ... Procrustes: rmse 0.004203617  max resid 0.006429414 
    ## ... Similar to previous best
    ## Run 239 stress 9.98255e-05 
    ## ... Procrustes: rmse 0.0004931625  max resid 0.0008090142 
    ## ... Similar to previous best
    ## Run 240 stress 0.000738917 
    ## Run 241 stress 9.953109e-05 
    ## ... Procrustes: rmse 0.0004238696  max resid 0.0007596958 
    ## ... Similar to previous best
    ## Run 242 stress 9.933409e-05 
    ## ... Procrustes: rmse 0.0004557953  max resid 0.0007931651 
    ## ... Similar to previous best
    ## Run 243 stress 9.486489e-05 
    ## ... Procrustes: rmse 0.0003997406  max resid 0.0007230516 
    ## ... Similar to previous best
    ## Run 244 stress 9.979485e-05 
    ## ... Procrustes: rmse 0.00047677  max resid 0.0007815479 
    ## ... Similar to previous best
    ## Run 245 stress 9.347039e-05 
    ## ... Procrustes: rmse 0.0003774018  max resid 0.0005935969 
    ## ... Similar to previous best
    ## Run 246 stress 9.715467e-05 
    ## ... Procrustes: rmse 0.0004365185  max resid 0.0007370754 
    ## ... Similar to previous best
    ## Run 247 stress 0.3407335 
    ## Run 248 stress 0.0007413625 
    ## Run 249 stress 0.001109941 
    ## Run 250 stress 9.737917e-05 
    ## ... Procrustes: rmse 0.0004399052  max resid 0.0007350489 
    ## ... Similar to previous best
    ## Run 251 stress 9.858542e-05 
    ## ... Procrustes: rmse 0.0004176087  max resid 0.0007591955 
    ## ... Similar to previous best
    ## Run 252 stress 0.001005761 
    ## Run 253 stress 9.837256e-05 
    ## ... Procrustes: rmse 0.0004826018  max resid 0.0007923293 
    ## ... Similar to previous best
    ## Run 254 stress 9.862396e-05 
    ## ... Procrustes: rmse 0.0004390038  max resid 0.0006752367 
    ## ... Similar to previous best
    ## Run 255 stress 8.689253e-05 
    ## ... Procrustes: rmse 6.353691e-05  max resid 0.0001497417 
    ## ... Similar to previous best
    ## Run 256 stress 0.0001264407 
    ## ... Procrustes: rmse 0.001089677  max resid 0.0016915 
    ## ... Similar to previous best
    ## Run 257 stress 9.893757e-05 
    ## ... Procrustes: rmse 0.00045612  max resid 0.0007909464 
    ## ... Similar to previous best
    ## Run 258 stress 0.0001658369 
    ## ... Procrustes: rmse 0.001433736  max resid 0.002226546 
    ## ... Similar to previous best
    ## Run 259 stress 0.0007726666 
    ## Run 260 stress 9.858395e-05 
    ## ... Procrustes: rmse 0.0004174197  max resid 0.0007565205 
    ## ... Similar to previous best
    ## Run 261 stress 9.863693e-05 
    ## ... Procrustes: rmse 0.0004180684  max resid 0.0007602214 
    ## ... Similar to previous best
    ## Run 262 stress 0.0007687166 
    ## Run 263 stress 9.733119e-05 
    ## ... Procrustes: rmse 0.0004167546  max resid 0.0007513466 
    ## ... Similar to previous best
    ## Run 264 stress 9.863699e-05 
    ## ... Procrustes: rmse 0.0003995546  max resid 0.0006237725 
    ## ... Similar to previous best
    ## Run 265 stress 9.856516e-05 
    ## ... Procrustes: rmse 0.0004044823  max resid 0.000626352 
    ## ... Similar to previous best
    ## Run 266 stress 9.958722e-05 
    ## ... Procrustes: rmse 0.0004216667  max resid 0.0007670122 
    ## ... Similar to previous best
    ## Run 267 stress 9.771342e-05 
    ## ... Procrustes: rmse 0.0005964486  max resid 0.0008694918 
    ## ... Similar to previous best
    ## Run 268 stress 0.0009192079 
    ## Run 269 stress 0.0006220727 
    ## Run 270 stress 0.001327914 
    ## Run 271 stress 9.959482e-05 
    ## ... Procrustes: rmse 0.0004262888  max resid 0.0007654323 
    ## ... Similar to previous best
    ## Run 272 stress 0.0008141963 
    ## Run 273 stress 0.0008290386 
    ## Run 274 stress 9.881265e-05 
    ## ... Procrustes: rmse 0.0004625079  max resid 0.0007900428 
    ## ... Similar to previous best
    ## Run 275 stress 0.0003149961 
    ## ... Procrustes: rmse 0.002694115  max resid 0.004143304 
    ## ... Similar to previous best
    ## Run 276 stress 0.0006510614 
    ## Run 277 stress 0.0002994286 
    ## ... Procrustes: rmse 0.002640936  max resid 0.004027321 
    ## ... Similar to previous best
    ## Run 278 stress 0.000196444 
    ## ... Procrustes: rmse 0.001706828  max resid 0.00264398 
    ## ... Similar to previous best
    ## Run 279 stress 9.823142e-05 
    ## ... Procrustes: rmse 0.0004619748  max resid 0.0007887083 
    ## ... Similar to previous best
    ## Run 280 stress 9.932849e-05 
    ## ... Procrustes: rmse 0.0004163653  max resid 0.0007477461 
    ## ... Similar to previous best
    ## Run 281 stress 0.0008791149 
    ## Run 282 stress 9.925826e-05 
    ## ... Procrustes: rmse 0.0004317685  max resid 0.0007351496 
    ## ... Similar to previous best
    ## Run 283 stress 0.0005903601 
    ## Run 284 stress 0.0006168519 
    ## Run 285 stress 9.991046e-05 
    ## ... Procrustes: rmse 0.0004440354  max resid 0.0006800397 
    ## ... Similar to previous best
    ## Run 286 stress 9.760088e-05 
    ## ... Procrustes: rmse 0.0004576465  max resid 0.0007813323 
    ## ... Similar to previous best
    ## Run 287 stress 0.0006521923 
    ## Run 288 stress 9.962531e-05 
    ## ... Procrustes: rmse 0.0004252818  max resid 0.0007660141 
    ## ... Similar to previous best
    ## Run 289 stress 9.966911e-05 
    ## ... Procrustes: rmse 0.0004913493  max resid 0.0008071147 
    ## ... Similar to previous best
    ## Run 290 stress 9.921907e-05 
    ## ... Procrustes: rmse 0.0004453554  max resid 0.0007423188 
    ## ... Similar to previous best
    ## Run 291 stress 9.816297e-05 
    ## ... Procrustes: rmse 0.0004161711  max resid 0.0007535261 
    ## ... Similar to previous best
    ## Run 292 stress 0.0008487832 
    ## Run 293 stress 0.0006830488 
    ## Run 294 stress 9.849473e-05 
    ## ... Procrustes: rmse 0.0004170085  max resid 0.0007548187 
    ## ... Similar to previous best
    ## Run 295 stress 9.849763e-05 
    ## ... Procrustes: rmse 0.0004379904  max resid 0.0006713685 
    ## ... Similar to previous best
    ## Run 296 stress 9.82966e-05 
    ## ... Procrustes: rmse 0.0004619299  max resid 0.0007891152 
    ## ... Similar to previous best
    ## Run 297 stress 0.0009777189 
    ## Run 298 stress 9.561811e-05 
    ## ... Procrustes: rmse 0.0004305787  max resid 0.0007422987 
    ## ... Similar to previous best
    ## Run 299 stress 9.850877e-05 
    ## ... Procrustes: rmse 0.0004863868  max resid 0.0007985146 
    ## ... Similar to previous best
    ## Run 300 stress 9.964758e-05 
    ## ... Procrustes: rmse 0.0004626663  max resid 0.0008030401 
    ## ... Similar to previous best
    ## Run 301 stress 0.0007977539 
    ## Run 302 stress 9.960094e-05 
    ## ... Procrustes: rmse 0.0004680566  max resid 0.0007988167 
    ## ... Similar to previous best
    ## Run 303 stress 9.882151e-05 
    ## ... Procrustes: rmse 0.0004457213  max resid 0.0007440188 
    ## ... Similar to previous best
    ## Run 304 stress 9.985187e-05 
    ## ... Procrustes: rmse 0.0004276064  max resid 0.0007664603 
    ## ... Similar to previous best
    ## Run 305 stress 0.0005790562 
    ## ... Procrustes: rmse 0.005083969  max resid 0.007798686 
    ## Run 306 stress 9.674505e-05 
    ## ... Procrustes: rmse 0.0004123905  max resid 0.0007481318 
    ## ... Similar to previous best
    ## Run 307 stress 9.9575e-05 
    ## ... Procrustes: rmse 0.0003812299  max resid 0.0006257616 
    ## ... Similar to previous best
    ## Run 308 stress 0.001262408 
    ## Run 309 stress 9.875675e-05 
    ## ... Procrustes: rmse 0.0004186074  max resid 0.0007600881 
    ## ... Similar to previous best
    ## Run 310 stress 9.954205e-05 
    ## ... Procrustes: rmse 0.0004440643  max resid 0.0007518695 
    ## ... Similar to previous best
    ## Run 311 stress 9.874263e-05 
    ## ... Procrustes: rmse 0.000416155  max resid 0.0007603144 
    ## ... Similar to previous best
    ## Run 312 stress 0.0007272256 
    ## Run 313 stress 9.811989e-05 
    ## ... Procrustes: rmse 0.0004155602  max resid 0.0007582504 
    ## ... Similar to previous best
    ## Run 314 stress 9.544119e-05 
    ## ... Procrustes: rmse 0.0003869158  max resid 0.0007202832 
    ## ... Similar to previous best
    ## Run 315 stress 9.826946e-05 
    ## ... Procrustes: rmse 0.0003975786  max resid 0.0006202177 
    ## ... Similar to previous best
    ## Run 316 stress 9.981581e-05 
    ## ... Procrustes: rmse 0.0008563284  max resid 0.001333057 
    ## ... Similar to previous best
    ## Run 317 stress 9.858783e-05 
    ## ... Procrustes: rmse 0.0004409979  max resid 0.0007477289 
    ## ... Similar to previous best
    ## Run 318 stress 9.820824e-05 
    ## ... Procrustes: rmse 0.0004377702  max resid 0.000672518 
    ## ... Similar to previous best
    ## Run 319 stress 9.974641e-05 
    ## ... Procrustes: rmse 0.0004622669  max resid 0.0008021584 
    ## ... Similar to previous best
    ## Run 320 stress 0.001170576 
    ## Run 321 stress 9.914825e-05 
    ## ... Procrustes: rmse 0.0004955608  max resid 0.0007974699 
    ## ... Similar to previous best
    ## Run 322 stress 9.895307e-05 
    ## ... Procrustes: rmse 0.0004651607  max resid 0.0007950087 
    ## ... Similar to previous best
    ## Run 323 stress 0.001073024 
    ## Run 324 stress 0.0008206202 
    ## Run 325 stress 0.0009146263 
    ## Run 326 stress 9.853912e-05 
    ## ... Procrustes: rmse 0.0004632523  max resid 0.0007906699 
    ## ... Similar to previous best
    ## Run 327 stress 9.95859e-05 
    ## ... Procrustes: rmse 0.0004577542  max resid 0.0007965884 
    ## ... Similar to previous best
    ## Run 328 stress 0.0006466919 
    ## Run 329 stress 0.0002959441 
    ## ... Procrustes: rmse 0.002564557  max resid 0.003979195 
    ## ... Similar to previous best
    ## Run 330 stress 0.0006189984 
    ## Run 331 stress 9.770277e-05 
    ## ... Procrustes: rmse 0.0004176691  max resid 0.0007535578 
    ## ... Similar to previous best
    ## Run 332 stress 9.967935e-05 
    ## ... Procrustes: rmse 0.0004128683  max resid 0.0007407078 
    ## ... Similar to previous best
    ## Run 333 stress 9.709822e-05 
    ## ... Procrustes: rmse 0.0004156117  max resid 0.0007488419 
    ## ... Similar to previous best
    ## Run 334 stress 0.0007393404 
    ## Run 335 stress 9.732421e-05 
    ## ... Procrustes: rmse 0.0004090974  max resid 0.0007425891 
    ## ... Similar to previous best
    ## Run 336 stress 0.0008772884 
    ## Run 337 stress 0.0007084096 
    ## Run 338 stress 0.3407342 
    ## Run 339 stress 0.000625451 
    ## Run 340 stress 9.877997e-05 
    ## ... Procrustes: rmse 0.0004000669  max resid 0.0006275113 
    ## ... Similar to previous best
    ## Run 341 stress 9.903366e-05 
    ## ... Procrustes: rmse 0.000486811  max resid 0.0008022413 
    ## ... Similar to previous best
    ## Run 342 stress 0.001139502 
    ## Run 343 stress 9.822437e-05 
    ## ... Procrustes: rmse 0.0004840708  max resid 0.0007951646 
    ## ... Similar to previous best
    ## Run 344 stress 9.834508e-05 
    ## ... Procrustes: rmse 0.0004859984  max resid 0.0007986186 
    ## ... Similar to previous best
    ## Run 345 stress 0.0006942164 
    ## Run 346 stress 0.0004817191 
    ## ... Procrustes: rmse 0.003999209  max resid 0.006102148 
    ## ... Similar to previous best
    ## Run 347 stress 0.000908809 
    ## Run 348 stress 0.0007168706 
    ## Run 349 stress 9.964012e-05 
    ## ... Procrustes: rmse 0.0004627558  max resid 0.0008031768 
    ## ... Similar to previous best
    ## Run 350 stress 0.000568911 
    ## ... Procrustes: rmse 0.004996444  max resid 0.007669445 
    ## ... Similar to previous best
    ## Run 351 stress 9.98513e-05 
    ## ... Procrustes: rmse 0.0003969072  max resid 0.0006145371 
    ## ... Similar to previous best
    ## Run 352 stress 0.0008662753 
    ## Run 353 stress 0.0008893499 
    ## Run 354 stress 0.0005258392 
    ## ... Procrustes: rmse 0.004621006  max resid 0.007086279 
    ## ... Similar to previous best
    ## Run 355 stress 9.973397e-05 
    ## ... Procrustes: rmse 0.0004226202  max resid 0.0007648437 
    ## ... Similar to previous best
    ## Run 356 stress 0.0004617178 
    ## ... Procrustes: rmse 0.004046774  max resid 0.006267389 
    ## ... Similar to previous best
    ## Run 357 stress 9.903773e-05 
    ## ... Procrustes: rmse 0.0004473273  max resid 0.0007456812 
    ## ... Similar to previous best
    ## Run 358 stress 9.861639e-05 
    ## ... Procrustes: rmse 0.0003996155  max resid 0.0006231398 
    ## ... Similar to previous best
    ## Run 359 stress 9.757236e-05 
    ## ... Procrustes: rmse 0.0004281786  max resid 0.0006556044 
    ## ... Similar to previous best
    ## Run 360 stress 0.0008405217 
    ## Run 361 stress 9.93108e-05 
    ## ... Procrustes: rmse 0.0004170416  max resid 0.0007478539 
    ## ... Similar to previous best
    ## Run 362 stress 9.849683e-05 
    ## ... Procrustes: rmse 0.0003995314  max resid 0.000623822 
    ## ... Similar to previous best
    ## Run 363 stress 9.933296e-05 
    ## ... Procrustes: rmse 0.0004963024  max resid 0.0008023441 
    ## ... Similar to previous best
    ## Run 364 stress 0.0003637227 
    ## ... Procrustes: rmse 0.003183122  max resid 0.004917014 
    ## ... Similar to previous best
    ## Run 365 stress 9.910604e-05 
    ## ... Procrustes: rmse 0.0004194399  max resid 0.0007571739 
    ## ... Similar to previous best
    ## Run 366 stress 0.000634632 
    ## Run 367 stress 0.0006490635 
    ## Run 368 stress 0.001118545 
    ## Run 369 stress 9.864069e-05 
    ## ... Procrustes: rmse 0.0004178897  max resid 0.0007565521 
    ## ... Similar to previous best
    ## Run 370 stress 9.887381e-05 
    ## ... Procrustes: rmse 0.0004461639  max resid 0.000745334 
    ## ... Similar to previous best
    ## Run 371 stress 9.922181e-05 
    ## ... Procrustes: rmse 0.0004441586  max resid 0.0007456081 
    ## ... Similar to previous best
    ## Run 372 stress 9.794526e-05 
    ## ... Procrustes: rmse 0.0004409795  max resid 0.0007327087 
    ## ... Similar to previous best
    ## Run 373 stress 9.648265e-05 
    ## ... Procrustes: rmse 0.0004062277  max resid 0.000745693 
    ## ... Similar to previous best
    ## Run 374 stress 0.001036947 
    ## Run 375 stress 9.971304e-05 
    ## ... Procrustes: rmse 0.0004360037  max resid 0.000666005 
    ## ... Similar to previous best
    ## Run 376 stress 9.815116e-05 
    ## ... Procrustes: rmse 0.0004361439  max resid 0.0006723681 
    ## ... Similar to previous best
    ## Run 377 stress 0.0006727751 
    ## Run 378 stress 9.931451e-05 
    ## ... Procrustes: rmse 0.0004459508  max resid 0.0007567243 
    ## ... Similar to previous best
    ## Run 379 stress 9.877982e-05 
    ## ... Procrustes: rmse 0.0003998913  max resid 0.0006243871 
    ## ... Similar to previous best
    ## Run 380 stress 0.001137083 
    ## Run 381 stress 9.996448e-05 
    ## ... Procrustes: rmse 0.0004221758  max resid 0.0007542765 
    ## ... Similar to previous best
    ## Run 382 stress 0.0002887739 
    ## ... Procrustes: rmse 0.002522473  max resid 0.003901006 
    ## ... Similar to previous best
    ## Run 383 stress 0.0006582383 
    ## Run 384 stress 0.0005700019 
    ## ... Procrustes: rmse 0.004998337  max resid 0.007681509 
    ## ... Similar to previous best
    ## Run 385 stress 9.871541e-05 
    ## ... Procrustes: rmse 0.000418049  max resid 0.0007613934 
    ## ... Similar to previous best
    ## Run 386 stress 0.246707 
    ## Run 387 stress 9.887193e-05 
    ## ... Procrustes: rmse 0.000400685  max resid 0.0006247519 
    ## ... Similar to previous best
    ## Run 388 stress 9.858585e-05 
    ## ... Procrustes: rmse 0.0003988523  max resid 0.0006211412 
    ## ... Similar to previous best
    ## Run 389 stress 0.001252274 
    ## Run 390 stress 9.687448e-05 
    ## ... Procrustes: rmse 0.0004122226  max resid 0.0007459562 
    ## ... Similar to previous best
    ## Run 391 stress 0.0007773385 
    ## Run 392 stress 0.001076961 
    ## Run 393 stress 0.001180913 
    ## Run 394 stress 0.00101735 
    ## Run 395 stress 0.0002806052 
    ## ... Procrustes: rmse 0.002438358  max resid 0.003698451 
    ## ... Similar to previous best
    ## Run 396 stress 9.92206e-05 
    ## ... Procrustes: rmse 0.0004023528  max resid 0.0006304203 
    ## ... Similar to previous best
    ## Run 397 stress 0.001069564 
    ## Run 398 stress 0.00106756 
    ## Run 399 stress 9.767008e-05 
    ## ... Procrustes: rmse 0.0004124812  max resid 0.0007580489 
    ## ... Similar to previous best
    ## Run 400 stress 9.710859e-05 
    ## ... Procrustes: rmse 0.0004064596  max resid 0.0007311797 
    ## ... Similar to previous best
    ## Run 401 stress 9.720258e-05 
    ## ... Procrustes: rmse 0.0003940832  max resid 0.0006155582 
    ## ... Similar to previous best
    ## Run 402 stress 0.001091444 
    ## Run 403 stress 9.764388e-05 
    ## ... Procrustes: rmse 0.0004309947  max resid 0.0006686683 
    ## ... Similar to previous best
    ## Run 404 stress 9.934586e-05 
    ## ... Procrustes: rmse 0.000461684  max resid 0.0008018043 
    ## ... Similar to previous best
    ## Run 405 stress 0.0007226087 
    ## Run 406 stress 9.822331e-05 
    ## ... Procrustes: rmse 0.0004202747  max resid 0.0007535125 
    ## ... Similar to previous best
    ## Run 407 stress 0.001270163 
    ## Run 408 stress 0.001374763 
    ## Run 409 stress 9.828666e-05 
    ## ... Procrustes: rmse 0.000392107  max resid 0.0006682634 
    ## ... Similar to previous best
    ## Run 410 stress 9.641653e-05 
    ## ... Procrustes: rmse 0.0004345064  max resid 0.0007290488 
    ## ... Similar to previous best
    ## Run 411 stress 9.928941e-05 
    ## ... Procrustes: rmse 0.000424016  max resid 0.0007607711 
    ## ... Similar to previous best
    ## Run 412 stress 9.893133e-05 
    ## ... Procrustes: rmse 0.0004596771  max resid 0.0007983125 
    ## ... Similar to previous best
    ## Run 413 stress 9.849911e-05 
    ## ... Procrustes: rmse 0.0003994289  max resid 0.0006247652 
    ## ... Similar to previous best
    ## Run 414 stress 0.0008599878 
    ## Run 415 stress 9.811245e-05 
    ## ... Procrustes: rmse 0.000414315  max resid 0.0007484683 
    ## ... Similar to previous best
    ## Run 416 stress 9.964846e-05 
    ## ... Procrustes: rmse 0.0004437141  max resid 0.0006797203 
    ## ... Similar to previous best
    ## Run 417 stress 0.0007507855 
    ## Run 418 stress 9.895181e-05 
    ## ... Procrustes: rmse 0.0004392635  max resid 0.0006758382 
    ## ... Similar to previous best
    ## Run 419 stress 0.0007519838 
    ## Run 420 stress 9.947095e-05 
    ## ... Procrustes: rmse 0.000462188  max resid 0.0008027102 
    ## ... Similar to previous best
    ## Run 421 stress 9.533508e-05 
    ## ... Procrustes: rmse 0.0004205366  max resid 0.0007279055 
    ## ... Similar to previous best
    ## Run 422 stress 9.82439e-05 
    ## ... Procrustes: rmse 0.0003974774  max resid 0.0006229838 
    ## ... Similar to previous best
    ## Run 423 stress 0.0005079635 
    ## ... Procrustes: rmse 0.004464346  max resid 0.006840551 
    ## ... Similar to previous best
    ## Run 424 stress 0.00102744 
    ## Run 425 stress 9.818524e-05 
    ## ... Procrustes: rmse 0.0004196072  max resid 0.0007554106 
    ## ... Similar to previous best
    ## Run 426 stress 9.787705e-05 
    ## ... Procrustes: rmse 0.0004605824  max resid 0.0007871617 
    ## ... Similar to previous best
    ## Run 427 stress 9.75316e-05 
    ## ... Procrustes: rmse 0.0004342919  max resid 0.0006681171 
    ## ... Similar to previous best
    ## Run 428 stress 0.0004061535 
    ## ... Procrustes: rmse 0.003574553  max resid 0.005464461 
    ## ... Similar to previous best
    ## Run 429 stress 9.963296e-05 
    ## ... Procrustes: rmse 0.000496838  max resid 0.0007955645 
    ## ... Similar to previous best
    ## Run 430 stress 9.93209e-05 
    ## ... Procrustes: rmse 0.0004419498  max resid 0.0006790289 
    ## ... Similar to previous best
    ## Run 431 stress 9.948562e-05 
    ## ... Procrustes: rmse 0.0004194564  max resid 0.0007559804 
    ## ... Similar to previous best
    ## Run 432 stress 9.99086e-05 
    ## ... Procrustes: rmse 0.0004005137  max resid 0.000621574 
    ## ... Similar to previous best
    ## Run 433 stress 9.816582e-05 
    ## ... Procrustes: rmse 0.0003361819  max resid 0.000693022 
    ## ... Similar to previous best
    ## Run 434 stress 9.983021e-05 
    ## ... Procrustes: rmse 0.0004258487  max resid 0.0007628994 
    ## ... Similar to previous best
    ## Run 435 stress 0.0008172185 
    ## Run 436 stress 9.833018e-05 
    ## ... Procrustes: rmse 0.0004416844  max resid 0.0007502022 
    ## ... Similar to previous best
    ## Run 437 stress 0.3275209 
    ## Run 438 stress 9.6336e-05 
    ## ... Procrustes: rmse 0.0004294734  max resid 0.0006604023 
    ## ... Similar to previous best
    ## Run 439 stress 9.937268e-05 
    ## ... Procrustes: rmse 0.0004461973  max resid 0.0007567795 
    ## ... Similar to previous best
    ## Run 440 stress 0.0008292393 
    ## Run 441 stress 9.562902e-05 
    ## ... Procrustes: rmse 0.0004268449  max resid 0.0006582192 
    ## ... Similar to previous best
    ## Run 442 stress 9.761381e-05 
    ## ... Procrustes: rmse 0.0003958562  max resid 0.0006191485 
    ## ... Similar to previous best
    ## Run 443 stress 9.913817e-05 
    ## ... Procrustes: rmse 0.000402319  max resid 0.000628066 
    ## ... Similar to previous best
    ## Run 444 stress 9.993682e-05 
    ## ... Procrustes: rmse 0.0005002804  max resid 0.0008042684 
    ## ... Similar to previous best
    ## Run 445 stress 0.000703423 
    ## Run 446 stress 9.818858e-05 
    ## ... Procrustes: rmse 0.0004150566  max resid 0.000757711 
    ## ... Similar to previous best
    ## Run 447 stress 0.0006097926 
    ## Run 448 stress 0.0006082116 
    ## Run 449 stress 0.0003986191 
    ## ... Procrustes: rmse 0.003509  max resid 0.005368573 
    ## ... Similar to previous best
    ## Run 450 stress 9.920637e-05 
    ## ... Procrustes: rmse 0.0004017285  max resid 0.0006294488 
    ## ... Similar to previous best
    ## Run 451 stress 9.863294e-05 
    ## ... Procrustes: rmse 0.000418159  max resid 0.0007580298 
    ## ... Similar to previous best
    ## Run 452 stress 0.001129468 
    ## Run 453 stress 9.983767e-05 
    ## ... Procrustes: rmse 0.0004693766  max resid 0.0008010886 
    ## ... Similar to previous best
    ## Run 454 stress 9.748595e-05 
    ## ... Procrustes: rmse 0.0004362113  max resid 0.0007383236 
    ## ... Similar to previous best
    ## Run 455 stress 9.487411e-05 
    ## ... Procrustes: rmse 0.0004057703  max resid 0.0007330144 
    ## ... Similar to previous best
    ## Run 456 stress 9.945287e-05 
    ## ... Procrustes: rmse 0.000402385  max resid 0.0006307758 
    ## ... Similar to previous best
    ## Run 457 stress 9.943923e-05 
    ## ... Procrustes: rmse 0.0004521116  max resid 0.0007884321 
    ## ... Similar to previous best
    ## Run 458 stress 9.913019e-05 
    ## ... Procrustes: rmse 0.00042466  max resid 0.0007618965 
    ## ... Similar to previous best
    ## Run 459 stress 9.955463e-05 
    ## ... Procrustes: rmse 0.0004448248  max resid 0.0007430646 
    ## ... Similar to previous best
    ## Run 460 stress 9.639425e-05 
    ## ... Procrustes: rmse 0.0004334869  max resid 0.0007372069 
    ## ... Similar to previous best
    ## Run 461 stress 9.992679e-05 
    ## ... Procrustes: rmse 0.0008096682  max resid 0.001288407 
    ## ... Similar to previous best
    ## Run 462 stress 9.787431e-05 
    ## ... Procrustes: rmse 0.000396997  max resid 0.0006212654 
    ## ... Similar to previous best
    ## Run 463 stress 9.886494e-05 
    ## ... Procrustes: rmse 0.000396845  max resid 0.0006243279 
    ## ... Similar to previous best
    ## Run 464 stress 0.000646285 
    ## Run 465 stress 9.982588e-05 
    ## ... Procrustes: rmse 0.0004484096  max resid 0.0007603184 
    ## ... Similar to previous best
    ## Run 466 stress 9.821928e-05 
    ## ... Procrustes: rmse 0.0004161465  max resid 0.0007587307 
    ## ... Similar to previous best
    ## Run 467 stress 9.976575e-05 
    ## ... Procrustes: rmse 0.0004471533  max resid 0.0007581673 
    ## ... Similar to previous best
    ## Run 468 stress 9.93282e-05 
    ## ... Procrustes: rmse 0.0003104216  max resid 0.0006659303 
    ## ... Similar to previous best
    ## Run 469 stress 0.00116003 
    ## Run 470 stress 0.001093503 
    ## Run 471 stress 9.95909e-05 
    ## ... Procrustes: rmse 0.0004494265  max resid 0.0007495745 
    ## ... Similar to previous best
    ## Run 472 stress 9.547716e-05 
    ## ... Procrustes: rmse 0.0001302906  max resid 0.0002928738 
    ## ... Similar to previous best
    ## Run 473 stress 9.843661e-05 
    ## ... Procrustes: rmse 0.0004185968  max resid 0.0007536346 
    ## ... Similar to previous best
    ## Run 474 stress 0.0009876051 
    ## Run 475 stress 9.953356e-05 
    ## ... Procrustes: rmse 0.0004029609  max resid 0.0006313621 
    ## ... Similar to previous best
    ## Run 476 stress 9.967999e-05 
    ## ... Procrustes: rmse 0.0001982529  max resid 0.0004942244 
    ## ... Similar to previous best
    ## Run 477 stress 9.897164e-05 
    ## ... Procrustes: rmse 0.0004890439  max resid 0.0008031802 
    ## ... Similar to previous best
    ## Run 478 stress 9.770759e-05 
    ## ... Procrustes: rmse 0.000435685  max resid 0.0006693601 
    ## ... Similar to previous best
    ## Run 479 stress 0.0008704278 
    ## Run 480 stress 0.00139757 
    ## Run 481 stress 9.95608e-05 
    ## ... Procrustes: rmse 0.0004211871  max resid 0.000766975 
    ## ... Similar to previous best
    ## Run 482 stress 9.79618e-05 
    ## ... Procrustes: rmse 0.0004363653  max resid 0.0006697056 
    ## ... Similar to previous best
    ## Run 483 stress 9.808604e-05 
    ## ... Procrustes: rmse 0.0004432911  max resid 0.0007405282 
    ## ... Similar to previous best
    ## Run 484 stress 9.855461e-05 
    ## ... Procrustes: rmse 0.0004380367  max resid 0.0007460774 
    ## ... Similar to previous best
    ## Run 485 stress 0.3266727 
    ## Run 486 stress 0.000216061 
    ## ... Procrustes: rmse 0.001879859  max resid 0.002911512 
    ## ... Similar to previous best
    ## Run 487 stress 9.867339e-05 
    ## ... Procrustes: rmse 0.0004180754  max resid 0.000757395 
    ## ... Similar to previous best
    ## Run 488 stress 0.001191152 
    ## Run 489 stress 0.0005882491 
    ## Run 490 stress 0.0005792878 
    ## ... Procrustes: rmse 0.00508699  max resid 0.00780607 
    ## Run 491 stress 9.920314e-05 
    ## ... Procrustes: rmse 0.0004841674  max resid 0.000797481 
    ## ... Similar to previous best
    ## Run 492 stress 9.863806e-05 
    ## ... Procrustes: rmse 0.0004578121  max resid 0.0007966174 
    ## ... Similar to previous best
    ## Run 493 stress 9.986527e-05 
    ## ... Procrustes: rmse 0.0004987648  max resid 0.0008003233 
    ## ... Similar to previous best
    ## Run 494 stress 0.0007480492 
    ## Run 495 stress 0.0003275279 
    ## ... Procrustes: rmse 0.002886786  max resid 0.004407512 
    ## ... Similar to previous best
    ## Run 496 stress 9.967448e-05 
    ## ... Procrustes: rmse 0.000447445  max resid 0.0007487576 
    ## ... Similar to previous best
    ## Run 497 stress 9.97225e-05 
    ## ... Procrustes: rmse 0.0004268809  max resid 0.0007661382 
    ## ... Similar to previous best
    ## Run 498 stress 9.991181e-05 
    ## ... Procrustes: rmse 0.0004421628  max resid 0.0006760296 
    ## ... Similar to previous best
    ## Run 499 stress 9.945585e-05 
    ## ... Procrustes: rmse 0.0004666956  max resid 0.00079574 
    ## ... Similar to previous best
    ## Run 500 stress 9.986954e-05 
    ## ... Procrustes: rmse 0.0004232  max resid 0.0007677959 
    ## ... Similar to previous best
    ## *** Best solution repeated 256 times

    ## Warning in metaMDS(PD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
PD_beta_env_M_NMDS <- metaMDS(PD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.019724e-05 
    ## Run 1 stress 9.627152e-05 
    ## ... Procrustes: rmse 0.0002334003  max resid 0.0003798012 
    ## ... Similar to previous best
    ## Run 2 stress 9.240751e-05 
    ## ... Procrustes: rmse 0.0001723021  max resid 0.0003493094 
    ## ... Similar to previous best
    ## Run 3 stress 9.5742e-05 
    ## ... Procrustes: rmse 1.450529e-05  max resid 1.902703e-05 
    ## ... Similar to previous best
    ## Run 4 stress 8.297251e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001203791  max resid 0.0002400659 
    ## ... Similar to previous best
    ## Run 5 stress 8.816652e-05 
    ## ... Procrustes: rmse 0.0001180955  max resid 0.0002194197 
    ## ... Similar to previous best
    ## Run 6 stress 9.340238e-05 
    ## ... Procrustes: rmse 0.0001347601  max resid 0.0003067964 
    ## ... Similar to previous best
    ## Run 7 stress 9.263928e-05 
    ## ... Procrustes: rmse 0.0001981148  max resid 0.0003452707 
    ## ... Similar to previous best
    ## Run 8 stress 9.478553e-05 
    ## ... Procrustes: rmse 0.000128665  max resid 0.0002471699 
    ## ... Similar to previous best
    ## Run 9 stress 9.105154e-05 
    ## ... Procrustes: rmse 0.0001345059  max resid 0.0002961751 
    ## ... Similar to previous best
    ## Run 10 stress 9.415231e-05 
    ## ... Procrustes: rmse 0.0001282059  max resid 0.0002475273 
    ## ... Similar to previous best
    ## Run 11 stress 9.767961e-05 
    ## ... Procrustes: rmse 0.0001580653  max resid 0.0003569856 
    ## ... Similar to previous best
    ## Run 12 stress 0.2852157 
    ## Run 13 stress 8.658612e-05 
    ## ... Procrustes: rmse 0.000180112  max resid 0.0003125779 
    ## ... Similar to previous best
    ## Run 14 stress 8.515755e-05 
    ## ... Procrustes: rmse 0.0001728349  max resid 0.0003166229 
    ## ... Similar to previous best
    ## Run 15 stress 9.402398e-05 
    ## ... Procrustes: rmse 0.0001398794  max resid 0.0003312163 
    ## ... Similar to previous best
    ## Run 16 stress 9.770648e-05 
    ## ... Procrustes: rmse 0.0001876476  max resid 0.0003124592 
    ## ... Similar to previous best
    ## Run 17 stress 9.703414e-05 
    ## ... Procrustes: rmse 0.0001938323  max resid 0.0003152714 
    ## ... Similar to previous best
    ## Run 18 stress 9.482891e-05 
    ## ... Procrustes: rmse 0.0001976081  max resid 0.0003441053 
    ## ... Similar to previous best
    ## Run 19 stress 9.394369e-05 
    ## ... Procrustes: rmse 0.0001444035  max resid 0.0002945466 
    ## ... Similar to previous best
    ## Run 20 stress 9.994486e-05 
    ## ... Procrustes: rmse 0.000120703  max resid 0.0002436694 
    ## ... Similar to previous best
    ## Run 21 stress 9.81733e-05 
    ## ... Procrustes: rmse 0.0001858245  max resid 0.0003157557 
    ## ... Similar to previous best
    ## Run 22 stress 0.2680462 
    ## Run 23 stress 9.533112e-05 
    ## ... Procrustes: rmse 0.0001309482  max resid 0.0002594649 
    ## ... Similar to previous best
    ## Run 24 stress 9.951907e-05 
    ## ... Procrustes: rmse 0.0001375622  max resid 0.000267302 
    ## ... Similar to previous best
    ## Run 25 stress 9.221579e-05 
    ## ... Procrustes: rmse 0.0002009615  max resid 0.0003476149 
    ## ... Similar to previous best
    ## Run 26 stress 9.699069e-05 
    ## ... Procrustes: rmse 0.0001275156  max resid 0.0002445073 
    ## ... Similar to previous best
    ## Run 27 stress 8.98743e-05 
    ## ... Procrustes: rmse 0.0001208362  max resid 0.0002424432 
    ## ... Similar to previous best
    ## Run 28 stress 9.240512e-05 
    ## ... Procrustes: rmse 0.0001978434  max resid 0.0003417045 
    ## ... Similar to previous best
    ## Run 29 stress 9.762635e-05 
    ## ... Procrustes: rmse 0.0001982654  max resid 0.0003452579 
    ## ... Similar to previous best
    ## Run 30 stress 9.828549e-05 
    ## ... Procrustes: rmse 0.0001338155  max resid 0.0002619485 
    ## ... Similar to previous best
    ## Run 31 stress 6.953052e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 5.258198e-05  max resid 6.837806e-05 
    ## ... Similar to previous best
    ## Run 32 stress 9.013523e-05 
    ## ... Procrustes: rmse 0.0001257581  max resid 0.0002803654 
    ## ... Similar to previous best
    ## Run 33 stress 9.968868e-05 
    ## ... Procrustes: rmse 0.0001412674  max resid 0.0002266782 
    ## ... Similar to previous best
    ## Run 34 stress 9.678884e-05 
    ## ... Procrustes: rmse 0.000199065  max resid 0.0002926953 
    ## ... Similar to previous best
    ## Run 35 stress 9.757067e-05 
    ## ... Procrustes: rmse 0.0001448391  max resid 0.000221124 
    ## ... Similar to previous best
    ## Run 36 stress 9.031908e-05 
    ## ... Procrustes: rmse 0.0001171975  max resid 0.0001827801 
    ## ... Similar to previous best
    ## Run 37 stress 0.2971485 
    ## Run 38 stress 9.300844e-05 
    ## ... Procrustes: rmse 0.0001347787  max resid 0.0002116854 
    ## ... Similar to previous best
    ## Run 39 stress 8.074924e-05 
    ## ... Procrustes: rmse 9.066996e-05  max resid 0.0001845212 
    ## ... Similar to previous best
    ## Run 40 stress 9.96425e-05 
    ## ... Procrustes: rmse 0.0002082771  max resid 0.000302176 
    ## ... Similar to previous best
    ## Run 41 stress 9.666659e-05 
    ## ... Procrustes: rmse 0.0001430461  max resid 0.000227136 
    ## ... Similar to previous best
    ## Run 42 stress 9.695368e-05 
    ## ... Procrustes: rmse 0.0001381889  max resid 0.0002202745 
    ## ... Similar to previous best
    ## Run 43 stress 9.84867e-05 
    ## ... Procrustes: rmse 0.000140629  max resid 0.000224915 
    ## ... Similar to previous best
    ## Run 44 stress 9.745938e-05 
    ## ... Procrustes: rmse 0.0001435754  max resid 0.0002240677 
    ## ... Similar to previous best
    ## Run 45 stress 9.291763e-05 
    ## ... Procrustes: rmse 0.0001196876  max resid 0.0002766687 
    ## ... Similar to previous best
    ## Run 46 stress 9.163597e-05 
    ## ... Procrustes: rmse 0.0001211311  max resid 0.0002611003 
    ## ... Similar to previous best
    ## Run 47 stress 8.740848e-05 
    ## ... Procrustes: rmse 0.0001163375  max resid 0.0002007884 
    ## ... Similar to previous best
    ## Run 48 stress 9.904555e-05 
    ## ... Procrustes: rmse 0.0001319907  max resid 0.000292883 
    ## ... Similar to previous best
    ## Run 49 stress 9.330333e-05 
    ## ... Procrustes: rmse 0.0001353675  max resid 0.0002115682 
    ## ... Similar to previous best
    ## Run 50 stress 8.529196e-05 
    ## ... Procrustes: rmse 0.0001637019  max resid 0.000289339 
    ## ... Similar to previous best
    ## Run 51 stress 9.761873e-05 
    ## ... Procrustes: rmse 0.0001449308  max resid 0.0002305867 
    ## ... Similar to previous best
    ## Run 52 stress 9.189125e-05 
    ## ... Procrustes: rmse 0.0001306644  max resid 0.0002068915 
    ## ... Similar to previous best
    ## Run 53 stress 9.320568e-05 
    ## ... Procrustes: rmse 0.0001344304  max resid 0.0002135044 
    ## ... Similar to previous best
    ## Run 54 stress 9.572746e-05 
    ## ... Procrustes: rmse 0.0001686651  max resid 0.0002597268 
    ## ... Similar to previous best
    ## Run 55 stress 9.097413e-05 
    ## ... Procrustes: rmse 0.0001267086  max resid 0.0002760074 
    ## ... Similar to previous best
    ## Run 56 stress 8.621943e-05 
    ## ... Procrustes: rmse 0.0001113722  max resid 0.00018403 
    ## ... Similar to previous best
    ## Run 57 stress 9.930941e-05 
    ## ... Procrustes: rmse 0.0001489432  max resid 0.0002353052 
    ## ... Similar to previous best
    ## Run 58 stress 9.623449e-05 
    ## ... Procrustes: rmse 0.0001379578  max resid 0.0002207238 
    ## ... Similar to previous best
    ## Run 59 stress 9.843922e-05 
    ## ... Procrustes: rmse 0.0001992451  max resid 0.0002927917 
    ## ... Similar to previous best
    ## Run 60 stress 9.355916e-05 
    ## ... Procrustes: rmse 0.0001359029  max resid 0.000211964 
    ## ... Similar to previous best
    ## Run 61 stress 0.3082927 
    ## Run 62 stress 8.592328e-05 
    ## ... Procrustes: rmse 0.0001150957  max resid 0.0001986281 
    ## ... Similar to previous best
    ## Run 63 stress 9.687003e-05 
    ## ... Procrustes: rmse 0.0001459943  max resid 0.0002703294 
    ## ... Similar to previous best
    ## Run 64 stress 9.102464e-05 
    ## ... Procrustes: rmse 0.0001441578  max resid 0.0003003393 
    ## ... Similar to previous best
    ## Run 65 stress 9.376274e-05 
    ## ... Procrustes: rmse 0.0001943834  max resid 0.0002897809 
    ## ... Similar to previous best
    ## Run 66 stress 9.627721e-05 
    ## ... Procrustes: rmse 0.0001642603  max resid 0.0002492556 
    ## ... Similar to previous best
    ## Run 67 stress 8.050709e-05 
    ## ... Procrustes: rmse 0.0001507741  max resid 0.0002644904 
    ## ... Similar to previous best
    ## Run 68 stress 9.765478e-05 
    ## ... Procrustes: rmse 0.0002037554  max resid 0.0002983005 
    ## ... Similar to previous best
    ## Run 69 stress 8.935528e-05 
    ## ... Procrustes: rmse 0.000113263  max resid 0.0002702529 
    ## ... Similar to previous best
    ## Run 70 stress 9.271021e-05 
    ## ... Procrustes: rmse 0.000135207  max resid 0.0002572518 
    ## ... Similar to previous best
    ## Run 71 stress 9.423272e-05 
    ## ... Procrustes: rmse 0.0001352034  max resid 0.0002164176 
    ## ... Similar to previous best
    ## Run 72 stress 9.771411e-05 
    ## ... Procrustes: rmse 0.0002003871  max resid 0.0002939645 
    ## ... Similar to previous best
    ## Run 73 stress 9.411969e-05 
    ## ... Procrustes: rmse 0.0001260993  max resid 0.0002504292 
    ## ... Similar to previous best
    ## Run 74 stress 9.913941e-05 
    ## ... Procrustes: rmse 0.0002073023  max resid 0.0002995979 
    ## ... Similar to previous best
    ## Run 75 stress 9.291352e-05 
    ## ... Procrustes: rmse 0.0001338824  max resid 0.0002110973 
    ## ... Similar to previous best
    ## Run 76 stress 9.981329e-05 
    ## ... Procrustes: rmse 0.0001483009  max resid 0.0002371422 
    ## ... Similar to previous best
    ## Run 77 stress 9.596023e-05 
    ## ... Procrustes: rmse 0.0001339769  max resid 0.0002961629 
    ## ... Similar to previous best
    ## Run 78 stress 9.094094e-05 
    ## ... Procrustes: rmse 0.0001426337  max resid 0.0002265683 
    ## ... Similar to previous best
    ## Run 79 stress 9.546181e-05 
    ## ... Procrustes: rmse 0.0001255153  max resid 0.0002724329 
    ## ... Similar to previous best
    ## Run 80 stress 9.962686e-05 
    ## ... Procrustes: rmse 0.0002019561  max resid 0.0002962834 
    ## ... Similar to previous best
    ## Run 81 stress 9.339591e-05 
    ## ... Procrustes: rmse 0.000138901  max resid 0.0002938198 
    ## ... Similar to previous best
    ## Run 82 stress 9.820574e-05 
    ## ... Procrustes: rmse 0.0001512221  max resid 0.0002770746 
    ## ... Similar to previous best
    ## Run 83 stress 9.293383e-05 
    ## ... Procrustes: rmse 0.0001345303  max resid 0.0002132654 
    ## ... Similar to previous best
    ## Run 84 stress 8.819263e-05 
    ## ... Procrustes: rmse 9.829487e-05  max resid 0.0002129509 
    ## ... Similar to previous best
    ## Run 85 stress 9.881514e-05 
    ## ... Procrustes: rmse 0.0001885877  max resid 0.0002966898 
    ## ... Similar to previous best
    ## Run 86 stress 9.548478e-05 
    ## ... Procrustes: rmse 0.0001383747  max resid 0.0002155131 
    ## ... Similar to previous best
    ## Run 87 stress 9.177513e-05 
    ## ... Procrustes: rmse 0.000171419  max resid 0.0003088039 
    ## ... Similar to previous best
    ## Run 88 stress 9.519693e-05 
    ## ... Procrustes: rmse 0.0001395488  max resid 0.0002199278 
    ## ... Similar to previous best
    ## Run 89 stress 9.110514e-05 
    ## ... Procrustes: rmse 0.0001863674  max resid 0.0002800167 
    ## ... Similar to previous best
    ## Run 90 stress 9.89596e-05 
    ## ... Procrustes: rmse 0.0001479543  max resid 0.000236123 
    ## ... Similar to previous best
    ## Run 91 stress 9.698588e-05 
    ## ... Procrustes: rmse 0.0001384779  max resid 0.0002216384 
    ## ... Similar to previous best
    ## Run 92 stress 8.620471e-05 
    ## ... Procrustes: rmse 0.0001263133  max resid 0.0002596332 
    ## ... Similar to previous best
    ## Run 93 stress 9.348839e-05 
    ## ... Procrustes: rmse 0.000193756  max resid 0.0002904383 
    ## ... Similar to previous best
    ## Run 94 stress 9.404576e-05 
    ## ... Procrustes: rmse 7.189591e-05  max resid 0.0001086343 
    ## ... Similar to previous best
    ## Run 95 stress 9.533792e-05 
    ## ... Procrustes: rmse 0.0001741188  max resid 0.0003108921 
    ## ... Similar to previous best
    ## Run 96 stress 7.13438e-05 
    ## ... Procrustes: rmse 0.000114906  max resid 0.0002271025 
    ## ... Similar to previous best
    ## Run 97 stress 9.796216e-05 
    ## ... Procrustes: rmse 0.0001546453  max resid 0.0003135465 
    ## ... Similar to previous best
    ## Run 98 stress 9.807309e-05 
    ## ... Procrustes: rmse 0.0002042681  max resid 0.000298752 
    ## ... Similar to previous best
    ## Run 99 stress 9.968691e-05 
    ## ... Procrustes: rmse 0.0001413775  max resid 0.0002194158 
    ## ... Similar to previous best
    ## Run 100 stress 9.695968e-05 
    ## ... Procrustes: rmse 0.000195362  max resid 0.0002914926 
    ## ... Similar to previous best
    ## Run 101 stress 9.706374e-05 
    ## ... Procrustes: rmse 0.0001203428  max resid 0.0002054305 
    ## ... Similar to previous best
    ## Run 102 stress 9.191196e-05 
    ## ... Procrustes: rmse 0.0001714458  max resid 0.0003084318 
    ## ... Similar to previous best
    ## Run 103 stress 9.486969e-05 
    ## ... Procrustes: rmse 0.000128628  max resid 0.000219077 
    ## ... Similar to previous best
    ## Run 104 stress 8.732799e-05 
    ## ... Procrustes: rmse 0.000119825  max resid 0.0001966531 
    ## ... Similar to previous best
    ## Run 105 stress 8.419942e-05 
    ## ... Procrustes: rmse 0.0001630607  max resid 0.0002785411 
    ## ... Similar to previous best
    ## Run 106 stress 9.010639e-05 
    ## ... Procrustes: rmse 0.0001517747  max resid 0.0002393916 
    ## ... Similar to previous best
    ## Run 107 stress 9.592856e-05 
    ## ... Procrustes: rmse 0.0001382786  max resid 0.0002159486 
    ## ... Similar to previous best
    ## Run 108 stress 9.18481e-05 
    ## ... Procrustes: rmse 0.0001861149  max resid 0.0002819685 
    ## ... Similar to previous best
    ## Run 109 stress 9.999357e-05 
    ## ... Procrustes: rmse 0.0001503705  max resid 0.0002428912 
    ## ... Similar to previous best
    ## Run 110 stress 8.520056e-05 
    ## ... Procrustes: rmse 0.000173923  max resid 0.0002709421 
    ## ... Similar to previous best
    ## Run 111 stress 9.717439e-05 
    ## ... Procrustes: rmse 0.0001427046  max resid 0.0002294808 
    ## ... Similar to previous best
    ## Run 112 stress 8.504777e-05 
    ## ... Procrustes: rmse 0.0001308144  max resid 0.0002530288 
    ## ... Similar to previous best
    ## Run 113 stress 9.516819e-05 
    ## ... Procrustes: rmse 0.0001081893  max resid 0.0002330312 
    ## ... Similar to previous best
    ## Run 114 stress 9.717744e-05 
    ## ... Procrustes: rmse 0.0001447207  max resid 0.0002303095 
    ## ... Similar to previous best
    ## Run 115 stress 9.293618e-05 
    ## ... Procrustes: rmse 0.0001354547  max resid 0.0002152515 
    ## ... Similar to previous best
    ## Run 116 stress 9.067676e-05 
    ## ... Procrustes: rmse 0.0001150264  max resid 0.000176733 
    ## ... Similar to previous best
    ## Run 117 stress 9.175029e-05 
    ## ... Procrustes: rmse 0.0001285588  max resid 0.0002817261 
    ## ... Similar to previous best
    ## Run 118 stress 8.539532e-05 
    ## ... Procrustes: rmse 5.171703e-05  max resid 7.312573e-05 
    ## ... Similar to previous best
    ## Run 119 stress 9.884418e-05 
    ## ... Procrustes: rmse 0.0001874342  max resid 0.0003039171 
    ## ... Similar to previous best
    ## Run 120 stress 8.843565e-05 
    ## ... Procrustes: rmse 0.0001581829  max resid 0.0002412293 
    ## ... Similar to previous best
    ## Run 121 stress 9.450083e-05 
    ## ... Procrustes: rmse 0.0001280015  max resid 0.000286225 
    ## ... Similar to previous best
    ## Run 122 stress 9.788185e-05 
    ## ... Procrustes: rmse 0.0001426937  max resid 0.0002218611 
    ## ... Similar to previous best
    ## Run 123 stress 9.841469e-05 
    ## ... Procrustes: rmse 0.0001265232  max resid 0.0001921934 
    ## ... Similar to previous best
    ## Run 124 stress 9.952243e-05 
    ## ... Procrustes: rmse 0.0001833935  max resid 0.0002853472 
    ## ... Similar to previous best
    ## Run 125 stress 9.801062e-05 
    ## ... Procrustes: rmse 0.0001322697  max resid 0.0002839514 
    ## ... Similar to previous best
    ## Run 126 stress 9.573677e-05 
    ## ... Procrustes: rmse 0.0001991682  max resid 0.0002946846 
    ## ... Similar to previous best
    ## Run 127 stress 9.532838e-05 
    ## ... Procrustes: rmse 0.000198075  max resid 0.0002954094 
    ## ... Similar to previous best
    ## Run 128 stress 9.56092e-05 
    ## ... Procrustes: rmse 0.0001978974  max resid 0.0002929973 
    ## ... Similar to previous best
    ## Run 129 stress 8.737623e-05 
    ## ... Procrustes: rmse 0.0001467821  max resid 0.0002318563 
    ## ... Similar to previous best
    ## Run 130 stress 0.2971485 
    ## Run 131 stress 8.588022e-05 
    ## ... Procrustes: rmse 0.0001380516  max resid 0.0002933856 
    ## ... Similar to previous best
    ## Run 132 stress 9.369168e-05 
    ## ... Procrustes: rmse 0.0001948132  max resid 0.000290951 
    ## ... Similar to previous best
    ## Run 133 stress 9.103766e-05 
    ## ... Procrustes: rmse 0.0001596454  max resid 0.0002486314 
    ## ... Similar to previous best
    ## Run 134 stress 9.774384e-05 
    ## ... Procrustes: rmse 0.0001424482  max resid 0.0002260405 
    ## ... Similar to previous best
    ## Run 135 stress 9.301652e-05 
    ## ... Procrustes: rmse 0.0001279338  max resid 0.0002844745 
    ## ... Similar to previous best
    ## Run 136 stress 9.833438e-05 
    ## ... Procrustes: rmse 0.0001513647  max resid 0.0003112856 
    ## ... Similar to previous best
    ## Run 137 stress 0.3078749 
    ## Run 138 stress 9.536794e-05 
    ## ... Procrustes: rmse 0.0001736058  max resid 0.0002873661 
    ## ... Similar to previous best
    ## Run 139 stress 9.509733e-05 
    ## ... Procrustes: rmse 0.000140021  max resid 0.0002175922 
    ## ... Similar to previous best
    ## Run 140 stress 9.938297e-05 
    ## ... Procrustes: rmse 0.0002020202  max resid 0.0002982649 
    ## ... Similar to previous best
    ## Run 141 stress 9.407801e-05 
    ## ... Procrustes: rmse 0.0001364637  max resid 0.0002127456 
    ## ... Similar to previous best
    ## Run 142 stress 9.74892e-05 
    ## ... Procrustes: rmse 0.0002021773  max resid 0.0002946476 
    ## ... Similar to previous best
    ## Run 143 stress 9.240812e-05 
    ## ... Procrustes: rmse 0.0001165531  max resid 0.0002653164 
    ## ... Similar to previous best
    ## Run 144 stress 8.856251e-05 
    ## ... Procrustes: rmse 0.0001561901  max resid 0.0002457067 
    ## ... Similar to previous best
    ## Run 145 stress 9.059296e-05 
    ## ... Procrustes: rmse 0.0001620674  max resid 0.0002836286 
    ## ... Similar to previous best
    ## Run 146 stress 9.12912e-05 
    ## ... Procrustes: rmse 0.0001312081  max resid 0.0002510672 
    ## ... Similar to previous best
    ## Run 147 stress 9.389857e-05 
    ## ... Procrustes: rmse 0.0001308789  max resid 0.0002131288 
    ## ... Similar to previous best
    ## Run 148 stress 9.646756e-05 
    ## ... Procrustes: rmse 0.0001397792  max resid 0.0002200079 
    ## ... Similar to previous best
    ## Run 149 stress 9.639275e-05 
    ## ... Procrustes: rmse 6.535904e-05  max resid 0.0001059141 
    ## ... Similar to previous best
    ## Run 150 stress 9.432912e-05 
    ## ... Procrustes: rmse 0.0001185991  max resid 0.0002371672 
    ## ... Similar to previous best
    ## Run 151 stress 9.316138e-05 
    ## ... Procrustes: rmse 0.00012845  max resid 0.000285175 
    ## ... Similar to previous best
    ## Run 152 stress 9.668591e-05 
    ## ... Procrustes: rmse 0.0001171706  max resid 0.0002785549 
    ## ... Similar to previous best
    ## Run 153 stress 9.89898e-05 
    ## ... Procrustes: rmse 0.0001480092  max resid 0.0002362232 
    ## ... Similar to previous best
    ## Run 154 stress 8.754003e-05 
    ## ... Procrustes: rmse 0.0001389164  max resid 0.0002880155 
    ## ... Similar to previous best
    ## Run 155 stress 9.840746e-05 
    ## ... Procrustes: rmse 0.000200962  max resid 0.0002975476 
    ## ... Similar to previous best
    ## Run 156 stress 9.951989e-05 
    ## ... Procrustes: rmse 0.0001417624  max resid 0.0002183209 
    ## ... Similar to previous best
    ## Run 157 stress 8.308866e-05 
    ## ... Procrustes: rmse 0.0001486872  max resid 0.0002325788 
    ## ... Similar to previous best
    ## Run 158 stress 7.887375e-05 
    ## ... Procrustes: rmse 0.0001111688  max resid 0.0002502405 
    ## ... Similar to previous best
    ## Run 159 stress 9.974426e-05 
    ## ... Procrustes: rmse 0.0001376666  max resid 0.0003021976 
    ## ... Similar to previous best
    ## Run 160 stress 8.476099e-05 
    ## ... Procrustes: rmse 0.0001110964  max resid 0.0002589091 
    ## ... Similar to previous best
    ## Run 161 stress 9.70417e-05 
    ## ... Procrustes: rmse 0.00020277  max resid 0.0002957922 
    ## ... Similar to previous best
    ## Run 162 stress 9.053401e-05 
    ## ... Procrustes: rmse 0.0001260913  max resid 0.0002019217 
    ## ... Similar to previous best
    ## Run 163 stress 8.937525e-05 
    ## ... Procrustes: rmse 0.0001264917  max resid 0.0002057712 
    ## ... Similar to previous best
    ## Run 164 stress 0.3082926 
    ## Run 165 stress 9.478205e-05 
    ## ... Procrustes: rmse 0.0001660673  max resid 0.0002558996 
    ## ... Similar to previous best
    ## Run 166 stress 9.223195e-05 
    ## ... Procrustes: rmse 0.0001917989  max resid 0.0002877963 
    ## ... Similar to previous best
    ## Run 167 stress 9.667356e-05 
    ## ... Procrustes: rmse 0.0001343128  max resid 0.0002902287 
    ## ... Similar to previous best
    ## Run 168 stress 0.2985159 
    ## Run 169 stress 8.002184e-05 
    ## ... Procrustes: rmse 0.0001515234  max resid 0.0002404551 
    ## ... Similar to previous best
    ## Run 170 stress 0.2854215 
    ## Run 171 stress 9.457249e-05 
    ## ... Procrustes: rmse 0.0001299565  max resid 0.0002883157 
    ## ... Similar to previous best
    ## Run 172 stress 0.284852 
    ## Run 173 stress 8.46966e-05 
    ## ... Procrustes: rmse 0.0001110767  max resid 0.0002590726 
    ## ... Similar to previous best
    ## Run 174 stress 8.318279e-05 
    ## ... Procrustes: rmse 0.0001404812  max resid 0.0002215718 
    ## ... Similar to previous best
    ## Run 175 stress 9.620193e-05 
    ## ... Procrustes: rmse 0.0001435813  max resid 0.0002944553 
    ## ... Similar to previous best
    ## Run 176 stress 9.570936e-05 
    ## ... Procrustes: rmse 0.0001406277  max resid 0.0002221721 
    ## ... Similar to previous best
    ## Run 177 stress 9.674939e-05 
    ## ... Procrustes: rmse 0.0001420513  max resid 0.0002255952 
    ## ... Similar to previous best
    ## Run 178 stress 9.617687e-05 
    ## ... Procrustes: rmse 0.0001404679  max resid 0.0002221542 
    ## ... Similar to previous best
    ## Run 179 stress 9.994303e-05 
    ## ... Procrustes: rmse 0.0001465644  max resid 0.0002331569 
    ## ... Similar to previous best
    ## Run 180 stress 9.494545e-05 
    ## ... Procrustes: rmse 0.0001975987  max resid 0.0002909497 
    ## ... Similar to previous best
    ## Run 181 stress 8.868927e-05 
    ## ... Procrustes: rmse 4.476315e-05  max resid 8.123355e-05 
    ## ... Similar to previous best
    ## Run 182 stress 9.716769e-05 
    ## ... Procrustes: rmse 6.64329e-05  max resid 0.0001052894 
    ## ... Similar to previous best
    ## Run 183 stress 8.049498e-05 
    ## ... Procrustes: rmse 0.0001629687  max resid 0.0002664887 
    ## ... Similar to previous best
    ## Run 184 stress 9.792155e-05 
    ## ... Procrustes: rmse 0.0001402732  max resid 0.000227875 
    ## ... Similar to previous best
    ## Run 185 stress 9.821254e-05 
    ## ... Procrustes: rmse 0.0002017617  max resid 0.0002978233 
    ## ... Similar to previous best
    ## Run 186 stress 9.807366e-05 
    ## ... Procrustes: rmse 7.50377e-05  max resid 0.000115716 
    ## ... Similar to previous best
    ## Run 187 stress 9.419697e-05 
    ## ... Procrustes: rmse 0.000129454  max resid 0.000287376 
    ## ... Similar to previous best
    ## Run 188 stress 9.472181e-05 
    ## ... Procrustes: rmse 0.000144844  max resid 0.0002636685 
    ## ... Similar to previous best
    ## Run 189 stress 9.504376e-05 
    ## ... Procrustes: rmse 0.0001301604  max resid 0.0002900722 
    ## ... Similar to previous best
    ## Run 190 stress 9.082188e-05 
    ## ... Procrustes: rmse 0.0001253082  max resid 0.0002792692 
    ## ... Similar to previous best
    ## Run 191 stress 9.058536e-05 
    ## ... Procrustes: rmse 0.0001215589  max resid 0.0002087425 
    ## ... Similar to previous best
    ## Run 192 stress 9.711559e-05 
    ## ... Procrustes: rmse 0.0001774447  max resid 0.0002400604 
    ## ... Similar to previous best
    ## Run 193 stress 8.836427e-05 
    ## ... Procrustes: rmse 0.0001601834  max resid 0.0002513836 
    ## ... Similar to previous best
    ## Run 194 stress 9.92702e-05 
    ## ... Procrustes: rmse 0.0001358079  max resid 0.0002912379 
    ## ... Similar to previous best
    ## Run 195 stress 7.918963e-05 
    ## ... Procrustes: rmse 0.0001096258  max resid 0.0002489997 
    ## ... Similar to previous best
    ## Run 196 stress 9.657925e-05 
    ## ... Procrustes: rmse 0.0001282311  max resid 0.0002013522 
    ## ... Similar to previous best
    ## Run 197 stress 0.2319196 
    ## Run 198 stress 9.847497e-05 
    ## ... Procrustes: rmse 0.0001247389  max resid 0.0001874489 
    ## ... Similar to previous best
    ## Run 199 stress 9.359532e-05 
    ## ... Procrustes: rmse 0.0001637513  max resid 0.0002524185 
    ## ... Similar to previous best
    ## Run 200 stress 9.57381e-05 
    ## ... Procrustes: rmse 0.0001628034  max resid 0.0002410911 
    ## ... Similar to previous best
    ## Run 201 stress 8.957512e-05 
    ## ... Procrustes: rmse 0.0001259008  max resid 0.0002025746 
    ## ... Similar to previous best
    ## Run 202 stress 8.852768e-05 
    ## ... Procrustes: rmse 0.0001139198  max resid 0.0002627158 
    ## ... Similar to previous best
    ## Run 203 stress 9.159108e-05 
    ## ... Procrustes: rmse 0.000129671  max resid 0.0002860899 
    ## ... Similar to previous best
    ## Run 204 stress 9.536348e-05 
    ## ... Procrustes: rmse 0.0001299071  max resid 0.0002000368 
    ## ... Similar to previous best
    ## Run 205 stress 9.175023e-05 
    ## ... Procrustes: rmse 0.0001749882  max resid 0.0002854886 
    ## ... Similar to previous best
    ## Run 206 stress 8.832131e-05 
    ## ... Procrustes: rmse 9.421571e-05  max resid 0.0002068486 
    ## ... Similar to previous best
    ## Run 207 stress 9.738089e-05 
    ## ... Procrustes: rmse 0.0002028698  max resid 0.0002977275 
    ## ... Similar to previous best
    ## Run 208 stress 9.27042e-05 
    ## ... Procrustes: rmse 0.0001594316  max resid 0.0002372728 
    ## ... Similar to previous best
    ## Run 209 stress 9.154806e-05 
    ## ... Procrustes: rmse 0.0001548247  max resid 0.0002211306 
    ## ... Similar to previous best
    ## Run 210 stress 9.975928e-05 
    ## ... Procrustes: rmse 0.0001386095  max resid 0.0002979583 
    ## ... Similar to previous best
    ## Run 211 stress 9.674718e-05 
    ## ... Procrustes: rmse 0.0001436104  max resid 0.0002285169 
    ## ... Similar to previous best
    ## Run 212 stress 9.259734e-05 
    ## ... Procrustes: rmse 0.0001547579  max resid 0.0002341858 
    ## ... Similar to previous best
    ## Run 213 stress 8.066243e-05 
    ## ... Procrustes: rmse 0.0001641093  max resid 0.0002622387 
    ## ... Similar to previous best
    ## Run 214 stress 5.612822e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 6.226539e-05  max resid 0.0001031501 
    ## ... Similar to previous best
    ## Run 215 stress 9.778513e-05 
    ## ... Procrustes: rmse 0.000182136  max resid 0.0002892696 
    ## ... Similar to previous best
    ## Run 216 stress 8.965911e-05 
    ## ... Procrustes: rmse 0.000109371  max resid 0.0001848057 
    ## ... Similar to previous best
    ## Run 217 stress 9.244012e-05 
    ## ... Procrustes: rmse 0.0001699023  max resid 0.0002655275 
    ## ... Similar to previous best
    ## Run 218 stress 9.526277e-05 
    ## ... Procrustes: rmse 0.0001743406  max resid 0.0002605598 
    ## ... Similar to previous best
    ## Run 219 stress 9.362187e-05 
    ## ... Procrustes: rmse 0.000183953  max resid 0.0002645896 
    ## ... Similar to previous best
    ## Run 220 stress 9.331396e-05 
    ## ... Procrustes: rmse 0.000146868  max resid 0.0002314478 
    ## ... Similar to previous best
    ## Run 221 stress 9.261258e-05 
    ## ... Procrustes: rmse 0.0001672009  max resid 0.0002684672 
    ## ... Similar to previous best
    ## Run 222 stress 8.550716e-05 
    ## ... Procrustes: rmse 0.0001144472  max resid 0.000173845 
    ## ... Similar to previous best
    ## Run 223 stress 9.676869e-05 
    ## ... Procrustes: rmse 0.0001971341  max resid 0.0002877797 
    ## ... Similar to previous best
    ## Run 224 stress 9.505562e-05 
    ## ... Procrustes: rmse 0.0001455545  max resid 0.0003147787 
    ## ... Similar to previous best
    ## Run 225 stress 8.826027e-05 
    ## ... Procrustes: rmse 0.0001545618  max resid 0.0002365235 
    ## ... Similar to previous best
    ## Run 226 stress 8.160041e-05 
    ## ... Procrustes: rmse 0.0001191793  max resid 0.0002459164 
    ## ... Similar to previous best
    ## Run 227 stress 9.362983e-05 
    ## ... Procrustes: rmse 0.0001419823  max resid 0.0002430873 
    ## ... Similar to previous best
    ## Run 228 stress 9.819664e-05 
    ## ... Procrustes: rmse 0.0001559695  max resid 0.000250783 
    ## ... Similar to previous best
    ## Run 229 stress 9.660354e-05 
    ## ... Procrustes: rmse 0.0001781366  max resid 0.0002760013 
    ## ... Similar to previous best
    ## Run 230 stress 9.172239e-05 
    ## ... Procrustes: rmse 0.0001152794  max resid 0.0001956788 
    ## ... Similar to previous best
    ## Run 231 stress 9.144404e-05 
    ## ... Procrustes: rmse 0.0001839247  max resid 0.0002638168 
    ## ... Similar to previous best
    ## Run 232 stress 9.97756e-05 
    ## ... Procrustes: rmse 0.0001671175  max resid 0.0002662824 
    ## ... Similar to previous best
    ## Run 233 stress 9.824653e-05 
    ## ... Procrustes: rmse 0.0001630677  max resid 0.0002587415 
    ## ... Similar to previous best
    ## Run 234 stress 9.645588e-05 
    ## ... Procrustes: rmse 0.0001798002  max resid 0.0002837119 
    ## ... Similar to previous best
    ## Run 235 stress 9.108563e-05 
    ## ... Procrustes: rmse 0.000125623  max resid 0.0001878543 
    ## ... Similar to previous best
    ## Run 236 stress 9.785296e-05 
    ## ... Procrustes: rmse 0.0001827867  max resid 0.0002944 
    ## ... Similar to previous best
    ## Run 237 stress 9.204793e-05 
    ## ... Procrustes: rmse 0.0001079465  max resid 0.0001579749 
    ## ... Similar to previous best
    ## Run 238 stress 9.182196e-05 
    ## ... Procrustes: rmse 0.0001127669  max resid 0.0002067217 
    ## ... Similar to previous best
    ## Run 239 stress 9.66528e-05 
    ## ... Procrustes: rmse 0.0001481512  max resid 0.0002749811 
    ## ... Similar to previous best
    ## Run 240 stress 7.17653e-05 
    ## ... Procrustes: rmse 6.670491e-05  max resid 0.0001048398 
    ## ... Similar to previous best
    ## Run 241 stress 8.918577e-05 
    ## ... Procrustes: rmse 0.0001109549  max resid 0.0002006456 
    ## ... Similar to previous best
    ## Run 242 stress 8.106177e-05 
    ## ... Procrustes: rmse 0.0001247725  max resid 0.0002579227 
    ## ... Similar to previous best
    ## Run 243 stress 9.937964e-05 
    ## ... Procrustes: rmse 0.0001940111  max resid 0.000282028 
    ## ... Similar to previous best
    ## Run 244 stress 9.431417e-05 
    ## ... Procrustes: rmse 0.0001746116  max resid 0.0002782752 
    ## ... Similar to previous best
    ## Run 245 stress 9.195799e-05 
    ## ... Procrustes: rmse 0.0001847536  max resid 0.0002586728 
    ## ... Similar to previous best
    ## Run 246 stress 9.529188e-05 
    ## ... Procrustes: rmse 0.0001755812  max resid 0.0002839782 
    ## ... Similar to previous best
    ## Run 247 stress 9.262714e-05 
    ## ... Procrustes: rmse 0.0001511365  max resid 0.0002434531 
    ## ... Similar to previous best
    ## Run 248 stress 0.3083099 
    ## Run 249 stress 8.710776e-05 
    ## ... Procrustes: rmse 0.0001255391  max resid 0.0002897565 
    ## ... Similar to previous best
    ## Run 250 stress 6.765252e-05 
    ## ... Procrustes: rmse 8.983548e-05  max resid 0.0001335376 
    ## ... Similar to previous best
    ## Run 251 stress 0.3083098 
    ## Run 252 stress 9.408991e-05 
    ## ... Procrustes: rmse 0.000165882  max resid 0.0002649575 
    ## ... Similar to previous best
    ## Run 253 stress 9.705498e-05 
    ## ... Procrustes: rmse 0.0001364502  max resid 0.0003195206 
    ## ... Similar to previous best
    ## Run 254 stress 9.948143e-05 
    ## ... Procrustes: rmse 0.000185103  max resid 0.0002943834 
    ## ... Similar to previous best
    ## Run 255 stress 9.576572e-05 
    ## ... Procrustes: rmse 0.0001764973  max resid 0.0002836051 
    ## ... Similar to previous best
    ## Run 256 stress 9.388102e-05 
    ## ... Procrustes: rmse 0.0001269396  max resid 0.0002341926 
    ## ... Similar to previous best
    ## Run 257 stress 9.459445e-05 
    ## ... Procrustes: rmse 0.0001702583  max resid 0.0002611089 
    ## ... Similar to previous best
    ## Run 258 stress 9.769691e-05 
    ## ... Procrustes: rmse 0.0001750465  max resid 0.0002824694 
    ## ... Similar to previous best
    ## Run 259 stress 9.467411e-05 
    ## ... Procrustes: rmse 0.000149764  max resid 0.0003146759 
    ## ... Similar to previous best
    ## Run 260 stress 9.539199e-05 
    ## ... Procrustes: rmse 0.0001497455  max resid 0.0003219973 
    ## ... Similar to previous best
    ## Run 261 stress 9.534314e-05 
    ## ... Procrustes: rmse 0.0001310304  max resid 0.0001889484 
    ## ... Similar to previous best
    ## Run 262 stress 9.092342e-05 
    ## ... Procrustes: rmse 0.0001358133  max resid 0.0001955611 
    ## ... Similar to previous best
    ## Run 263 stress 8.955212e-05 
    ## ... Procrustes: rmse 0.0001608637  max resid 0.0002535234 
    ## ... Similar to previous best
    ## Run 264 stress 0.3120129 
    ## Run 265 stress 9.373314e-05 
    ## ... Procrustes: rmse 0.0001456018  max resid 0.0002728047 
    ## ... Similar to previous best
    ## Run 266 stress 9.203381e-05 
    ## ... Procrustes: rmse 0.000102793  max resid 0.000171617 
    ## ... Similar to previous best
    ## Run 267 stress 9.563344e-05 
    ## ... Procrustes: rmse 0.0001524058  max resid 0.0002189761 
    ## ... Similar to previous best
    ## Run 268 stress 9.412157e-05 
    ## ... Procrustes: rmse 0.0001659037  max resid 0.0002667205 
    ## ... Similar to previous best
    ## Run 269 stress 0.2680462 
    ## Run 270 stress 9.696169e-05 
    ## ... Procrustes: rmse 0.0001763873  max resid 0.0002767708 
    ## ... Similar to previous best
    ## Run 271 stress 9.439824e-05 
    ## ... Procrustes: rmse 0.0001733313  max resid 0.0002709547 
    ## ... Similar to previous best
    ## Run 272 stress 9.07424e-05 
    ## ... Procrustes: rmse 0.0001359454  max resid 0.0002593187 
    ## ... Similar to previous best
    ## Run 273 stress 9.937389e-05 
    ## ... Procrustes: rmse 0.0001847674  max resid 0.0002906087 
    ## ... Similar to previous best
    ## Run 274 stress 9.560063e-05 
    ## ... Procrustes: rmse 0.0001892423  max resid 0.0002739908 
    ## ... Similar to previous best
    ## Run 275 stress 9.674472e-05 
    ## ... Procrustes: rmse 0.0001482066  max resid 0.0002758391 
    ## ... Similar to previous best
    ## Run 276 stress 9.853778e-05 
    ## ... Procrustes: rmse 0.0001822871  max resid 0.0002822545 
    ## ... Similar to previous best
    ## Run 277 stress 9.627901e-05 
    ## ... Procrustes: rmse 0.0001932158  max resid 0.000281611 
    ## ... Similar to previous best
    ## Run 278 stress 9.651561e-05 
    ## ... Procrustes: rmse 0.0001789251  max resid 0.0002857826 
    ## ... Similar to previous best
    ## Run 279 stress 9.636373e-05 
    ## ... Procrustes: rmse 0.0001785126  max resid 0.0002800332 
    ## ... Similar to previous best
    ## Run 280 stress 9.463184e-05 
    ## ... Procrustes: rmse 0.0001864628  max resid 0.0002597855 
    ## ... Similar to previous best
    ## Run 281 stress 9.967561e-05 
    ## ... Procrustes: rmse 0.0002058014  max resid 0.0003037044 
    ## ... Similar to previous best
    ## Run 282 stress 9.256477e-05 
    ## ... Procrustes: rmse 0.0001656537  max resid 0.0002615494 
    ## ... Similar to previous best
    ## Run 283 stress 9.199202e-05 
    ## ... Procrustes: rmse 0.0001377406  max resid 0.0002551207 
    ## ... Similar to previous best
    ## Run 284 stress 9.781563e-05 
    ## ... Procrustes: rmse 0.0001810979  max resid 0.0002879153 
    ## ... Similar to previous best
    ## Run 285 stress 9.37402e-05 
    ## ... Procrustes: rmse 0.000185109  max resid 0.0002692089 
    ## ... Similar to previous best
    ## Run 286 stress 9.562613e-05 
    ## ... Procrustes: rmse 0.0001771445  max resid 0.0002837946 
    ## ... Similar to previous best
    ## Run 287 stress 9.752068e-05 
    ## ... Procrustes: rmse 0.0001417481  max resid 0.0003284868 
    ## ... Similar to previous best
    ## Run 288 stress 9.396978e-05 
    ## ... Procrustes: rmse 0.0001430287  max resid 0.0002093019 
    ## ... Similar to previous best
    ## Run 289 stress 9.989781e-05 
    ## ... Procrustes: rmse 0.0001876412  max resid 0.000304523 
    ## ... Similar to previous best
    ## Run 290 stress 9.822678e-05 
    ## ... Procrustes: rmse 0.000151632  max resid 0.000277895 
    ## ... Similar to previous best
    ## Run 291 stress 9.822267e-05 
    ## ... Procrustes: rmse 0.0001489623  max resid 0.0002798601 
    ## ... Similar to previous best
    ## Run 292 stress 9.86198e-05 
    ## ... Procrustes: rmse 0.0002025289  max resid 0.0002965661 
    ## ... Similar to previous best
    ## Run 293 stress 9.428385e-05 
    ## ... Procrustes: rmse 0.0001342331  max resid 0.0002496215 
    ## ... Similar to previous best
    ## Run 294 stress 9.854045e-05 
    ## ... Procrustes: rmse 0.0001325101  max resid 0.0002362654 
    ## ... Similar to previous best
    ## Run 295 stress 9.785322e-05 
    ## ... Procrustes: rmse 0.0001392195  max resid 0.0003258424 
    ## ... Similar to previous best
    ## Run 296 stress 9.451628e-05 
    ## ... Procrustes: rmse 0.0001746309  max resid 0.000276516 
    ## ... Similar to previous best
    ## Run 297 stress 8.937784e-05 
    ## ... Procrustes: rmse 0.0001367501  max resid 0.0002954231 
    ## ... Similar to previous best
    ## Run 298 stress 9.341382e-05 
    ## ... Procrustes: rmse 9.160362e-05  max resid 0.0001427051 
    ## ... Similar to previous best
    ## Run 299 stress 0.3083098 
    ## Run 300 stress 9.648039e-05 
    ## ... Procrustes: rmse 0.0001955171  max resid 0.0002859552 
    ## ... Similar to previous best
    ## Run 301 stress 9.990202e-05 
    ## ... Procrustes: rmse 0.0001275134  max resid 0.0001807765 
    ## ... Similar to previous best
    ## Run 302 stress 9.057239e-05 
    ## ... Procrustes: rmse 0.0001394539  max resid 0.0002820302 
    ## ... Similar to previous best
    ## Run 303 stress 9.429587e-05 
    ## ... Procrustes: rmse 0.0001740728  max resid 0.0002731679 
    ## ... Similar to previous best
    ## Run 304 stress 9.951227e-05 
    ## ... Procrustes: rmse 0.0001829112  max resid 0.0002905313 
    ## ... Similar to previous best
    ## Run 305 stress 9.234343e-05 
    ## ... Procrustes: rmse 0.0001428043  max resid 0.0002350474 
    ## ... Similar to previous best
    ## Run 306 stress 9.964394e-05 
    ## ... Procrustes: rmse 0.0001323995  max resid 0.0002120508 
    ## ... Similar to previous best
    ## Run 307 stress 9.751465e-05 
    ## ... Procrustes: rmse 0.0001332738  max resid 0.0001976239 
    ## ... Similar to previous best
    ## Run 308 stress 9.513286e-05 
    ## ... Procrustes: rmse 0.000175839  max resid 0.0002825884 
    ## ... Similar to previous best
    ## Run 309 stress 9.109872e-05 
    ## ... Procrustes: rmse 0.0001364663  max resid 0.0002949248 
    ## ... Similar to previous best
    ## Run 310 stress 8.733709e-05 
    ## ... Procrustes: rmse 0.0001324142  max resid 0.0002532538 
    ## ... Similar to previous best
    ## Run 311 stress 9.124503e-05 
    ## ... Procrustes: rmse 0.0001252622  max resid 0.0002251691 
    ## ... Similar to previous best
    ## Run 312 stress 0.3120129 
    ## Run 313 stress 9.775595e-05 
    ## ... Procrustes: rmse 0.0001493219  max resid 0.0002767563 
    ## ... Similar to previous best
    ## Run 314 stress 9.177172e-05 
    ## ... Procrustes: rmse 0.0001648196  max resid 0.0002579877 
    ## ... Similar to previous best
    ## Run 315 stress 9.088009e-05 
    ## ... Procrustes: rmse 0.0001236724  max resid 0.0002193109 
    ## ... Similar to previous best
    ## Run 316 stress 0.3081116 
    ## Run 317 stress 9.960483e-05 
    ## ... Procrustes: rmse 0.0001071713  max resid 0.0001859384 
    ## ... Similar to previous best
    ## Run 318 stress 0.2842084 
    ## Run 319 stress 9.794559e-05 
    ## ... Procrustes: rmse 0.0001737371  max resid 0.0002810819 
    ## ... Similar to previous best
    ## Run 320 stress 9.763847e-05 
    ## ... Procrustes: rmse 0.000127746  max resid 0.0001834443 
    ## ... Similar to previous best
    ## Run 321 stress 7.919588e-05 
    ## ... Procrustes: rmse 0.0001123339  max resid 0.00024824 
    ## ... Similar to previous best
    ## Run 322 stress 9.428306e-05 
    ## ... Procrustes: rmse 0.0001096177  max resid 0.000185134 
    ## ... Similar to previous best
    ## Run 323 stress 9.246971e-05 
    ## ... Procrustes: rmse 0.0001790604  max resid 0.0002577585 
    ## ... Similar to previous best
    ## Run 324 stress 9.719888e-05 
    ## ... Procrustes: rmse 0.0001229333  max resid 0.000195771 
    ## ... Similar to previous best
    ## Run 325 stress 9.737776e-05 
    ## ... Procrustes: rmse 0.0001457294  max resid 0.0002686089 
    ## ... Similar to previous best
    ## Run 326 stress 9.845736e-05 
    ## ... Procrustes: rmse 0.0002032135  max resid 0.0002932636 
    ## ... Similar to previous best
    ## Run 327 stress 8.591229e-05 
    ## ... Procrustes: rmse 0.0001306944  max resid 0.0002823728 
    ## ... Similar to previous best
    ## Run 328 stress 9.412414e-05 
    ## ... Procrustes: rmse 0.0001924417  max resid 0.0002775625 
    ## ... Similar to previous best
    ## Run 329 stress 9.397642e-05 
    ## ... Procrustes: rmse 0.0001696668  max resid 0.0002714946 
    ## ... Similar to previous best
    ## Run 330 stress 8.665114e-05 
    ## ... Procrustes: rmse 0.0001158065  max resid 0.0001578739 
    ## ... Similar to previous best
    ## Run 331 stress 9.732238e-05 
    ## ... Procrustes: rmse 0.0001296099  max resid 0.0002373504 
    ## ... Similar to previous best
    ## Run 332 stress 8.256161e-05 
    ## ... Procrustes: rmse 0.0001124161  max resid 0.0001635276 
    ## ... Similar to previous best
    ## Run 333 stress 0.3078066 
    ## Run 334 stress 9.768612e-05 
    ## ... Procrustes: rmse 0.0001497247  max resid 0.0002773531 
    ## ... Similar to previous best
    ## Run 335 stress 9.177228e-05 
    ## ... Procrustes: rmse 0.0001379536  max resid 0.0002580797 
    ## ... Similar to previous best
    ## Run 336 stress 9.262258e-05 
    ## ... Procrustes: rmse 0.0001629394  max resid 0.0002178179 
    ## ... Similar to previous best
    ## Run 337 stress 9.697785e-05 
    ## ... Procrustes: rmse 0.0001804955  max resid 0.0002901766 
    ## ... Similar to previous best
    ## Run 338 stress 9.571453e-05 
    ## ... Procrustes: rmse 0.0001195397  max resid 0.0001805502 
    ## ... Similar to previous best
    ## Run 339 stress 8.959764e-05 
    ## ... Procrustes: rmse 0.0001774513  max resid 0.0002404648 
    ## ... Similar to previous best
    ## Run 340 stress 8.047842e-05 
    ## ... Procrustes: rmse 0.0001131522  max resid 0.0001718509 
    ## ... Similar to previous best
    ## Run 341 stress 9.317144e-05 
    ## ... Procrustes: rmse 0.0001217809  max resid 0.0002029624 
    ## ... Similar to previous best
    ## Run 342 stress 9.961293e-05 
    ## ... Procrustes: rmse 0.0002029044  max resid 0.0002935397 
    ## ... Similar to previous best
    ## Run 343 stress 9.454337e-05 
    ## ... Procrustes: rmse 0.0001741663  max resid 0.0002760326 
    ## ... Similar to previous best
    ## Run 344 stress 9.735425e-05 
    ## ... Procrustes: rmse 0.0001771259  max resid 0.0002801928 
    ## ... Similar to previous best
    ## Run 345 stress 9.559259e-05 
    ## ... Procrustes: rmse 0.0001446911  max resid 0.0002680595 
    ## ... Similar to previous best
    ## Run 346 stress 9.145395e-05 
    ## ... Procrustes: rmse 0.0001401008  max resid 0.0002633667 
    ## ... Similar to previous best
    ## Run 347 stress 9.298501e-05 
    ## ... Procrustes: rmse 0.0001813934  max resid 0.0002488017 
    ## ... Similar to previous best
    ## Run 348 stress 8.594035e-05 
    ## ... Procrustes: rmse 0.0001507492  max resid 0.0002398099 
    ## ... Similar to previous best
    ## Run 349 stress 9.799937e-05 
    ## ... Procrustes: rmse 0.0002000417  max resid 0.0002860314 
    ## ... Similar to previous best
    ## Run 350 stress 0.2848522 
    ## Run 351 stress 9.182393e-05 
    ## ... Procrustes: rmse 0.0001621432  max resid 0.0002589704 
    ## ... Similar to previous best
    ## Run 352 stress 9.909919e-05 
    ## ... Procrustes: rmse 0.0001511459  max resid 0.0002172976 
    ## ... Similar to previous best
    ## Run 353 stress 9.961418e-05 
    ## ... Procrustes: rmse 0.0001800894  max resid 0.000287427 
    ## ... Similar to previous best
    ## Run 354 stress 9.98468e-05 
    ## ... Procrustes: rmse 0.0001850518  max resid 0.0002932199 
    ## ... Similar to previous best
    ## Run 355 stress 8.588128e-05 
    ## ... Procrustes: rmse 0.0001187529  max resid 0.0002136137 
    ## ... Similar to previous best
    ## Run 356 stress 9.568287e-05 
    ## ... Procrustes: rmse 0.0001347205  max resid 0.0003066539 
    ## ... Similar to previous best
    ## Run 357 stress 9.07855e-05 
    ## ... Procrustes: rmse 0.0001416605  max resid 0.0002396086 
    ## ... Similar to previous best
    ## Run 358 stress 8.486665e-05 
    ## ... Procrustes: rmse 0.0001385027  max resid 0.0002164801 
    ## ... Similar to previous best
    ## Run 359 stress 9.168511e-05 
    ## ... Procrustes: rmse 0.0001227553  max resid 0.0002299996 
    ## ... Similar to previous best
    ## Run 360 stress 8.823746e-05 
    ## ... Procrustes: rmse 0.0001187466  max resid 0.0002049945 
    ## ... Similar to previous best
    ## Run 361 stress 9.882684e-05 
    ## ... Procrustes: rmse 0.0001856487  max resid 0.0002989877 
    ## ... Similar to previous best
    ## Run 362 stress 9.503446e-05 
    ## ... Procrustes: rmse 0.000145096  max resid 0.0002704973 
    ## ... Similar to previous best
    ## Run 363 stress 9.96857e-05 
    ## ... Procrustes: rmse 0.0001875284  max resid 0.0002975741 
    ## ... Similar to previous best
    ## Run 364 stress 9.477976e-05 
    ## ... Procrustes: rmse 0.0001874007  max resid 0.0002654028 
    ## ... Similar to previous best
    ## Run 365 stress 8.743192e-05 
    ## ... Procrustes: rmse 0.0001656638  max resid 0.0002324293 
    ## ... Similar to previous best
    ## Run 366 stress 8.976814e-05 
    ## ... Procrustes: rmse 0.0001394872  max resid 0.0002645498 
    ## ... Similar to previous best
    ## Run 367 stress 9.112298e-05 
    ## ... Procrustes: rmse 0.0001177329  max resid 0.0001989587 
    ## ... Similar to previous best
    ## Run 368 stress 7.302138e-05 
    ## ... Procrustes: rmse 8.24262e-05  max resid 0.0001251029 
    ## ... Similar to previous best
    ## Run 369 stress 9.916584e-05 
    ## ... Procrustes: rmse 0.0001865092  max resid 0.0003005746 
    ## ... Similar to previous best
    ## Run 370 stress 9.5525e-05 
    ## ... Procrustes: rmse 0.0001759975  max resid 0.0002751596 
    ## ... Similar to previous best
    ## Run 371 stress 9.868606e-05 
    ## ... Procrustes: rmse 0.0001813619  max resid 0.0002908236 
    ## ... Similar to previous best
    ## Run 372 stress 9.025764e-05 
    ## ... Procrustes: rmse 0.000136785  max resid 0.0002563559 
    ## ... Similar to previous best
    ## Run 373 stress 9.641823e-05 
    ## ... Procrustes: rmse 0.0001463824  max resid 0.0002705636 
    ## ... Similar to previous best
    ## Run 374 stress 8.446143e-05 
    ## ... Procrustes: rmse 0.0001176009  max resid 0.0002513517 
    ## ... Similar to previous best
    ## Run 375 stress 9.721906e-05 
    ## ... Procrustes: rmse 0.00013077  max resid 0.0002248505 
    ## ... Similar to previous best
    ## Run 376 stress 8.92045e-05 
    ## ... Procrustes: rmse 0.0001526071  max resid 0.0002405044 
    ## ... Similar to previous best
    ## Run 377 stress 8.680334e-05 
    ## ... Procrustes: rmse 0.0001190116  max resid 0.0001780309 
    ## ... Similar to previous best
    ## Run 378 stress 9.583242e-05 
    ## ... Procrustes: rmse 0.0001648525  max resid 0.0002684963 
    ## ... Similar to previous best
    ## Run 379 stress 0.3061064 
    ## Run 380 stress 0.2854402 
    ## Run 381 stress 9.561284e-05 
    ## ... Procrustes: rmse 0.0001476098  max resid 0.0002346951 
    ## ... Similar to previous best
    ## Run 382 stress 9.772333e-05 
    ## ... Procrustes: rmse 0.0001312392  max resid 0.0002252155 
    ## ... Similar to previous best
    ## Run 383 stress 0.285215 
    ## Run 384 stress 9.700953e-05 
    ## ... Procrustes: rmse 0.0001317819  max resid 0.0002107203 
    ## ... Similar to previous best
    ## Run 385 stress 8.511519e-05 
    ## ... Procrustes: rmse 0.0001139806  max resid 0.0002011401 
    ## ... Similar to previous best
    ## Run 386 stress 9.966938e-05 
    ## ... Procrustes: rmse 0.0001794634  max resid 0.000291164 
    ## ... Similar to previous best
    ## Run 387 stress 9.628884e-05 
    ## ... Procrustes: rmse 0.0001917795  max resid 0.0002757093 
    ## ... Similar to previous best
    ## Run 388 stress 0.3082926 
    ## Run 389 stress 0.2848518 
    ## Run 390 stress 9.830649e-05 
    ## ... Procrustes: rmse 0.000191616  max resid 0.0002813812 
    ## ... Similar to previous best
    ## Run 391 stress 9.21602e-05 
    ## ... Procrustes: rmse 0.0001588299  max resid 0.0002553194 
    ## ... Similar to previous best
    ## Run 392 stress 7.195194e-05 
    ## ... Procrustes: rmse 0.0001033982  max resid 0.0001590657 
    ## ... Similar to previous best
    ## Run 393 stress 9.367809e-05 
    ## ... Procrustes: rmse 0.0001903069  max resid 0.0002695797 
    ## ... Similar to previous best
    ## Run 394 stress 8.89732e-05 
    ## ... Procrustes: rmse 0.0001077448  max resid 0.0002376319 
    ## ... Similar to previous best
    ## Run 395 stress 9.538892e-05 
    ## ... Procrustes: rmse 0.0001312905  max resid 0.0002103518 
    ## ... Similar to previous best
    ## Run 396 stress 9.243508e-05 
    ## ... Procrustes: rmse 0.0001585822  max resid 0.0002548831 
    ## ... Similar to previous best
    ## Run 397 stress 9.947172e-05 
    ## ... Procrustes: rmse 0.0002056434  max resid 0.0003022544 
    ## ... Similar to previous best
    ## Run 398 stress 7.972531e-05 
    ## ... Procrustes: rmse 0.0001136883  max resid 0.0002002937 
    ## ... Similar to previous best
    ## Run 399 stress 9.622807e-05 
    ## ... Procrustes: rmse 0.0001973916  max resid 0.0002870235 
    ## ... Similar to previous best
    ## Run 400 stress 9.397265e-05 
    ## ... Procrustes: rmse 0.0001912711  max resid 0.0002699255 
    ## ... Similar to previous best
    ## Run 401 stress 9.331223e-05 
    ## ... Procrustes: rmse 0.0001294059  max resid 0.0001701009 
    ## ... Similar to previous best
    ## Run 402 stress 9.174909e-05 
    ## ... Procrustes: rmse 0.0001864222  max resid 0.0002672501 
    ## ... Similar to previous best
    ## Run 403 stress 9.389469e-05 
    ## ... Procrustes: rmse 0.0001850555  max resid 0.000266493 
    ## ... Similar to previous best
    ## Run 404 stress 0.3083098 
    ## Run 405 stress 9.2042e-05 
    ## ... Procrustes: rmse 9.983324e-05  max resid 0.0001246047 
    ## ... Similar to previous best
    ## Run 406 stress 9.386114e-05 
    ## ... Procrustes: rmse 0.0001418288  max resid 0.0002638513 
    ## ... Similar to previous best
    ## Run 407 stress 7.134987e-05 
    ## ... Procrustes: rmse 9.493229e-05  max resid 0.000206398 
    ## ... Similar to previous best
    ## Run 408 stress 9.851194e-05 
    ## ... Procrustes: rmse 0.0001613631  max resid 0.0002563565 
    ## ... Similar to previous best
    ## Run 409 stress 9.906042e-05 
    ## ... Procrustes: rmse 0.0001870179  max resid 0.0002970589 
    ## ... Similar to previous best
    ## Run 410 stress 9.647027e-05 
    ## ... Procrustes: rmse 0.000144943  max resid 0.0003119578 
    ## ... Similar to previous best
    ## Run 411 stress 9.293384e-05 
    ## ... Procrustes: rmse 0.0001475217  max resid 0.0002443869 
    ## ... Similar to previous best
    ## Run 412 stress 9.262021e-05 
    ## ... Procrustes: rmse 0.0001422484  max resid 0.000264992 
    ## ... Similar to previous best
    ## Run 413 stress 8.990592e-05 
    ## ... Procrustes: rmse 0.0001624769  max resid 0.0002592088 
    ## ... Similar to previous best
    ## Run 414 stress 9.45406e-05 
    ## ... Procrustes: rmse 0.0001430531  max resid 0.0002657365 
    ## ... Similar to previous best
    ## Run 415 stress 9.612872e-05 
    ## ... Procrustes: rmse 0.000195377  max resid 0.0002752444 
    ## ... Similar to previous best
    ## Run 416 stress 9.046832e-05 
    ## ... Procrustes: rmse 9.800233e-05  max resid 0.0001602042 
    ## ... Similar to previous best
    ## Run 417 stress 9.731689e-05 
    ## ... Procrustes: rmse 0.0001910029  max resid 0.0002779162 
    ## ... Similar to previous best
    ## Run 418 stress 9.466781e-05 
    ## ... Procrustes: rmse 0.0001749884  max resid 0.0002749546 
    ## ... Similar to previous best
    ## Run 419 stress 9.533508e-05 
    ## ... Procrustes: rmse 0.0001311433  max resid 0.0002526056 
    ## ... Similar to previous best
    ## Run 420 stress 9.842499e-05 
    ## ... Procrustes: rmse 0.000150844  max resid 0.0003261788 
    ## ... Similar to previous best
    ## Run 421 stress 9.924648e-05 
    ## ... Procrustes: rmse 0.0002034875  max resid 0.0002929606 
    ## ... Similar to previous best
    ## Run 422 stress 0.3083098 
    ## Run 423 stress 9.620214e-05 
    ## ... Procrustes: rmse 0.0001970191  max resid 0.0002863469 
    ## ... Similar to previous best
    ## Run 424 stress 8.338473e-05 
    ## ... Procrustes: rmse 0.0001293594  max resid 0.0002015116 
    ## ... Similar to previous best
    ## Run 425 stress 9.132959e-05 
    ## ... Procrustes: rmse 0.0001742685  max resid 0.0002508612 
    ## ... Similar to previous best
    ## Run 426 stress 9.605192e-05 
    ## ... Procrustes: rmse 0.0001921348  max resid 0.0002817925 
    ## ... Similar to previous best
    ## Run 427 stress 9.589038e-05 
    ## ... Procrustes: rmse 0.0001769105  max resid 0.0002777063 
    ## ... Similar to previous best
    ## Run 428 stress 9.457186e-05 
    ## ... Procrustes: rmse 0.0001286581  max resid 0.0003024977 
    ## ... Similar to previous best
    ## Run 429 stress 8.764078e-05 
    ## ... Procrustes: rmse 0.0001573061  max resid 0.0002451385 
    ## ... Similar to previous best
    ## Run 430 stress 9.460414e-05 
    ## ... Procrustes: rmse 0.0001443375  max resid 0.0002693694 
    ## ... Similar to previous best
    ## Run 431 stress 9.529468e-05 
    ## ... Procrustes: rmse 0.0001437087  max resid 0.0002716194 
    ## ... Similar to previous best
    ## Run 432 stress 8.874935e-05 
    ## ... Procrustes: rmse 0.0001245558  max resid 0.0002929099 
    ## ... Similar to previous best
    ## Run 433 stress 9.906562e-05 
    ## ... Procrustes: rmse 0.0001507053  max resid 0.0002465109 
    ## ... Similar to previous best
    ## Run 434 stress 9.70215e-05 
    ## ... Procrustes: rmse 0.0001726325  max resid 0.0002782791 
    ## ... Similar to previous best
    ## Run 435 stress 7.865971e-05 
    ## ... Procrustes: rmse 9.703723e-05  max resid 0.0001560663 
    ## ... Similar to previous best
    ## Run 436 stress 9.847983e-05 
    ## ... Procrustes: rmse 0.0001784259  max resid 0.000288138 
    ## ... Similar to previous best
    ## Run 437 stress 9.296335e-05 
    ## ... Procrustes: rmse 0.0001414719  max resid 0.0002651582 
    ## ... Similar to previous best
    ## Run 438 stress 8.955666e-05 
    ## ... Procrustes: rmse 0.0001582553  max resid 0.0002449978 
    ## ... Similar to previous best
    ## Run 439 stress 9.211603e-05 
    ## ... Procrustes: rmse 0.0001237152  max resid 0.00021944 
    ## ... Similar to previous best
    ## Run 440 stress 9.332341e-05 
    ## ... Procrustes: rmse 0.0001718446  max resid 0.000274317 
    ## ... Similar to previous best
    ## Run 441 stress 8.908139e-05 
    ## ... Procrustes: rmse 0.0001794513  max resid 0.0002543597 
    ## ... Similar to previous best
    ## Run 442 stress 9.940104e-05 
    ## ... Procrustes: rmse 0.0001487672  max resid 0.0003214771 
    ## ... Similar to previous best
    ## Run 443 stress 9.666153e-05 
    ## ... Procrustes: rmse 0.0001804647  max resid 0.0002894272 
    ## ... Similar to previous best
    ## Run 444 stress 9.599716e-05 
    ## ... Procrustes: rmse 0.0001674917  max resid 0.0002705313 
    ## ... Similar to previous best
    ## Run 445 stress 9.072493e-05 
    ## ... Procrustes: rmse 0.0001828307  max resid 0.0002593369 
    ## ... Similar to previous best
    ## Run 446 stress 9.891615e-05 
    ## ... Procrustes: rmse 0.0001986718  max resid 0.0002934985 
    ## ... Similar to previous best
    ## Run 447 stress 9.022701e-05 
    ## ... Procrustes: rmse 0.0001373619  max resid 0.0002581282 
    ## ... Similar to previous best
    ## Run 448 stress 9.427522e-05 
    ## ... Procrustes: rmse 0.0001310282  max resid 0.0003011879 
    ## ... Similar to previous best
    ## Run 449 stress 9.608732e-05 
    ## ... Procrustes: rmse 0.0001480816  max resid 0.0002751185 
    ## ... Similar to previous best
    ## Run 450 stress 9.712662e-05 
    ## ... Procrustes: rmse 0.0001497815  max resid 0.0002776907 
    ## ... Similar to previous best
    ## Run 451 stress 9.943144e-05 
    ## ... Procrustes: rmse 0.0001518113  max resid 0.0003252429 
    ## ... Similar to previous best
    ## Run 452 stress 9.364046e-05 
    ## ... Procrustes: rmse 0.0001732892  max resid 0.000272793 
    ## ... Similar to previous best
    ## Run 453 stress 9.690539e-05 
    ## ... Procrustes: rmse 0.0001636026  max resid 0.0002723559 
    ## ... Similar to previous best
    ## Run 454 stress 9.940248e-05 
    ## ... Procrustes: rmse 0.0001863164  max resid 0.0002998304 
    ## ... Similar to previous best
    ## Run 455 stress 9.436972e-05 
    ## ... Procrustes: rmse 0.0001911924  max resid 0.000277381 
    ## ... Similar to previous best
    ## Run 456 stress 9.241059e-05 
    ## ... Procrustes: rmse 0.0001390656  max resid 0.0002597066 
    ## ... Similar to previous best
    ## Run 457 stress 9.231561e-05 
    ## ... Procrustes: rmse 0.0001352022  max resid 0.0002557859 
    ## ... Similar to previous best
    ## Run 458 stress 8.623686e-05 
    ## ... Procrustes: rmse 0.0001427577  max resid 0.0002195429 
    ## ... Similar to previous best
    ## Run 459 stress 9.966028e-05 
    ## ... Procrustes: rmse 0.0001868868  max resid 0.0002980794 
    ## ... Similar to previous best
    ## Run 460 stress 9.692034e-05 
    ## ... Procrustes: rmse 0.000195884  max resid 0.0002839247 
    ## ... Similar to previous best
    ## Run 461 stress 9.375111e-05 
    ## ... Procrustes: rmse 0.0001420554  max resid 0.0002017029 
    ## ... Similar to previous best
    ## Run 462 stress 8.964356e-05 
    ## ... Procrustes: rmse 0.0001166337  max resid 0.000186786 
    ## ... Similar to previous best
    ## Run 463 stress 9.251561e-05 
    ## ... Procrustes: rmse 0.000163641  max resid 0.0002618283 
    ## ... Similar to previous best
    ## Run 464 stress 4.528323e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 7.876004e-05  max resid 0.0001121431 
    ## ... Similar to previous best
    ## Run 465 stress 0.2848519 
    ## Run 466 stress 8.834392e-05 
    ## ... Procrustes: rmse 9.283146e-05  max resid 0.0001526498 
    ## ... Similar to previous best
    ## Run 467 stress 9.911243e-05 
    ## ... Procrustes: rmse 0.0002506635  max resid 0.0003607773 
    ## ... Similar to previous best
    ## Run 468 stress 9.83036e-05 
    ## ... Procrustes: rmse 0.000187738  max resid 0.0003773464 
    ## ... Similar to previous best
    ## Run 469 stress 9.353364e-05 
    ## ... Procrustes: rmse 0.0001800264  max resid 0.0003562759 
    ## ... Similar to previous best
    ## Run 470 stress 9.40367e-05 
    ## ... Procrustes: rmse 0.0001002141  max resid 0.0001829222 
    ## ... Similar to previous best
    ## Run 471 stress 9.86385e-05 
    ## ... Procrustes: rmse 0.0002272008  max resid 0.0003588095 
    ## ... Similar to previous best
    ## Run 472 stress 9.725972e-05 
    ## ... Procrustes: rmse 0.0001423342  max resid 0.0002281655 
    ## ... Similar to previous best
    ## Run 473 stress 8.867001e-05 
    ## ... Procrustes: rmse 0.000112881  max resid 0.0002111496 
    ## ... Similar to previous best
    ## Run 474 stress 9.223543e-05 
    ## ... Procrustes: rmse 0.000172819  max resid 0.0002614436 
    ## ... Similar to previous best
    ## Run 475 stress 9.254344e-05 
    ## ... Procrustes: rmse 0.0001697856  max resid 0.0003446984 
    ## ... Similar to previous best
    ## Run 476 stress 8.730643e-05 
    ## ... Procrustes: rmse 0.0002190534  max resid 0.0003074028 
    ## ... Similar to previous best
    ## Run 477 stress 8.752558e-05 
    ## ... Procrustes: rmse 0.0001228396  max resid 0.0001996785 
    ## ... Similar to previous best
    ## Run 478 stress 0.3083098 
    ## Run 479 stress 9.332564e-05 
    ## ... Procrustes: rmse 0.0001805353  max resid 0.0002846052 
    ## ... Similar to previous best
    ## Run 480 stress 9.45511e-05 
    ## ... Procrustes: rmse 0.0002382595  max resid 0.0003382588 
    ## ... Similar to previous best
    ## Run 481 stress 8.844684e-05 
    ## ... Procrustes: rmse 0.0001461017  max resid 0.0002586971 
    ## ... Similar to previous best
    ## Run 482 stress 9.717075e-05 
    ## ... Procrustes: rmse 0.0001845138  max resid 0.0003710646 
    ## ... Similar to previous best
    ## Run 483 stress 9.968051e-05 
    ## ... Procrustes: rmse 0.0002300564  max resid 0.0003655399 
    ## ... Similar to previous best
    ## Run 484 stress 9.101172e-05 
    ## ... Procrustes: rmse 0.00013667  max resid 0.0002229848 
    ## ... Similar to previous best
    ## Run 485 stress 9.102436e-05 
    ## ... Procrustes: rmse 0.0001754472  max resid 0.000347775 
    ## ... Similar to previous best
    ## Run 486 stress 9.686198e-05 
    ## ... Procrustes: rmse 0.0002107205  max resid 0.0003385522 
    ## ... Similar to previous best
    ## Run 487 stress 9.531036e-05 
    ## ... Procrustes: rmse 0.0002110618  max resid 0.0003345798 
    ## ... Similar to previous best
    ## Run 488 stress 9.604668e-05 
    ## ... Procrustes: rmse 0.0002274828  max resid 0.0003338294 
    ## ... Similar to previous best
    ## Run 489 stress 9.306523e-05 
    ## ... Procrustes: rmse 0.0002117871  max resid 0.0003310269 
    ## ... Similar to previous best
    ## Run 490 stress 8.917082e-05 
    ## ... Procrustes: rmse 0.0001983037  max resid 0.0003097887 
    ## ... Similar to previous best
    ## Run 491 stress 9.869921e-05 
    ## ... Procrustes: rmse 0.0002483277  max resid 0.0003607257 
    ## ... Similar to previous best
    ## Run 492 stress 9.010538e-05 
    ## ... Procrustes: rmse 0.0001721227  max resid 0.0003482029 
    ## ... Similar to previous best
    ## Run 493 stress 9.457158e-05 
    ## ... Procrustes: rmse 0.0002367797  max resid 0.0003444953 
    ## ... Similar to previous best
    ## Run 494 stress 0.3083098 
    ## Run 495 stress 9.975412e-05 
    ## ... Procrustes: rmse 0.0001889062  max resid 0.0003718283 
    ## ... Similar to previous best
    ## Run 496 stress 9.841179e-05 
    ## ... Procrustes: rmse 0.0002200946  max resid 0.0003498608 
    ## ... Similar to previous best
    ## Run 497 stress 9.875366e-05 
    ## ... Procrustes: rmse 0.0002218856  max resid 0.0003589288 
    ## ... Similar to previous best
    ## Run 498 stress 9.208946e-05 
    ## ... Procrustes: rmse 9.568308e-05  max resid 0.0001418817 
    ## ... Similar to previous best
    ## Run 499 stress 9.88181e-05 
    ## ... Procrustes: rmse 0.000248778  max resid 0.0003576355 
    ## ... Similar to previous best
    ## Run 500 stress 9.461868e-05 
    ## ... Procrustes: rmse 0.0002119165  max resid 0.0003370552 
    ## ... Similar to previous best
    ## *** Best solution repeated 34 times

    ## Warning in metaMDS(PD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
PD_beta_geo_NMDS <- metaMDS(PD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07927535 
    ## Run 1 stress 0.1047632 
    ## Run 2 stress 0.08985616 
    ## Run 3 stress 0.08985605 
    ## Run 4 stress 0.09344087 
    ## Run 5 stress 0.07951441 
    ## ... Procrustes: rmse 0.01732471  max resid 0.05628178 
    ## Run 6 stress 0.09247107 
    ## Run 7 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734583  max resid 0.05645003 
    ## Run 8 stress 0.09106031 
    ## Run 9 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 4.238408e-05  max resid 0.0001149878 
    ## ... Similar to previous best
    ## Run 10 stress 0.0904376 
    ## Run 11 stress 0.07927535 
    ## ... Procrustes: rmse 3.383557e-05  max resid 9.358122e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.07936955 
    ## ... Procrustes: rmse 0.009024217  max resid 0.0323331 
    ## Run 13 stress 0.0956538 
    ## Run 14 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735267  max resid 0.05639622 
    ## Run 15 stress 0.08985606 
    ## Run 16 stress 0.09296198 
    ## Run 17 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 9.434154e-06  max resid 2.59402e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.09382093 
    ## Run 19 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237369  max resid 0.04965063 
    ## Run 20 stress 0.07969265 
    ## ... Procrustes: rmse 0.01237002  max resid 0.04956505 
    ## Run 21 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239108  max resid 0.04973282 
    ## Run 22 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734898  max resid 0.05638418 
    ## Run 23 stress 0.07927537 
    ## ... Procrustes: rmse 2.311849e-05  max resid 5.196821e-05 
    ## ... Similar to previous best
    ## Run 24 stress 0.08985617 
    ## Run 25 stress 0.08985606 
    ## Run 26 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046238  max resid 0.0324615 
    ## Run 27 stress 0.07927534 
    ## ... Procrustes: rmse 4.413846e-06  max resid 1.437259e-05 
    ## ... Similar to previous best
    ## Run 28 stress 0.07936949 
    ## ... Procrustes: rmse 0.009048498  max resid 0.03248423 
    ## Run 29 stress 0.09043774 
    ## Run 30 stress 0.07969263 
    ## ... Procrustes: rmse 0.0124034  max resid 0.04987215 
    ## Run 31 stress 0.07969267 
    ## ... Procrustes: rmse 0.01238854  max resid 0.04983614 
    ## Run 32 stress 0.09172782 
    ## Run 33 stress 0.08985606 
    ## Run 34 stress 0.0925223 
    ## Run 35 stress 0.0929071 
    ## Run 36 stress 0.09043764 
    ## Run 37 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 9.488034e-06  max resid 2.807564e-05 
    ## ... Similar to previous best
    ## Run 38 stress 0.0795144 
    ## ... Procrustes: rmse 0.01733852  max resid 0.05633519 
    ## Run 39 stress 0.09290703 
    ## Run 40 stress 0.0793695 
    ## ... Procrustes: rmse 0.009038142  max resid 0.03242142 
    ## Run 41 stress 0.0925222 
    ## Run 42 stress 0.09296257 
    ## Run 43 stress 0.09106048 
    ## Run 44 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734458  max resid 0.05636995 
    ## Run 45 stress 0.1076849 
    ## Run 46 stress 0.09043744 
    ## Run 47 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735451  max resid 0.05649965 
    ## Run 48 stress 0.07927537 
    ## ... Procrustes: rmse 5.891149e-05  max resid 0.0001620783 
    ## ... Similar to previous best
    ## Run 49 stress 0.09072902 
    ## Run 50 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734841  max resid 0.05639551 
    ## Run 51 stress 0.07969259 
    ## ... Procrustes: rmse 0.0123816  max resid 0.04971041 
    ## Run 52 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046209  max resid 0.03247517 
    ## Run 53 stress 0.0938834 
    ## Run 54 stress 0.07927534 
    ## ... Procrustes: rmse 1.612365e-05  max resid 4.667937e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.09565365 
    ## Run 56 stress 0.07936949 
    ## ... Procrustes: rmse 0.009044399  max resid 0.03246108 
    ## Run 57 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732828  max resid 0.0562435 
    ## Run 58 stress 0.07927534 
    ## ... Procrustes: rmse 4.638757e-06  max resid 1.38537e-05 
    ## ... Similar to previous best
    ## Run 59 stress 0.09106044 
    ## Run 60 stress 0.09106044 
    ## Run 61 stress 0.1097215 
    ## Run 62 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735451  max resid 0.0564748 
    ## Run 63 stress 0.09072893 
    ## Run 64 stress 0.08985614 
    ## Run 65 stress 0.07936955 
    ## ... Procrustes: rmse 0.009024084  max resid 0.03233385 
    ## Run 66 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735428  max resid 0.05645222 
    ## Run 67 stress 0.09022792 
    ## Run 68 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733086  max resid 0.05627475 
    ## Run 69 stress 0.09296236 
    ## Run 70 stress 0.07927536 
    ## ... Procrustes: rmse 5.146518e-05  max resid 0.0001448896 
    ## ... Similar to previous best
    ## Run 71 stress 0.07951443 
    ## ... Procrustes: rmse 0.01735874  max resid 0.05652212 
    ## Run 72 stress 0.07927534 
    ## ... Procrustes: rmse 1.160726e-05  max resid 4.606138e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.090729 
    ## Run 74 stress 0.07927535 
    ## ... Procrustes: rmse 3.496263e-05  max resid 8.05398e-05 
    ## ... Similar to previous best
    ## Run 75 stress 0.09252224 
    ## Run 76 stress 0.07969259 
    ## ... Procrustes: rmse 0.0123788  max resid 0.04972386 
    ## Run 77 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237398  max resid 0.04967044 
    ## Run 78 stress 0.09172784 
    ## Run 79 stress 0.07936953 
    ## ... Procrustes: rmse 0.009026834  max resid 0.03236309 
    ## Run 80 stress 0.09382079 
    ## Run 81 stress 0.07969265 
    ## ... Procrustes: rmse 0.01237194  max resid 0.04976338 
    ## Run 82 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173595  max resid 0.05647659 
    ## Run 83 stress 0.09106036 
    ## Run 84 stress 0.09565381 
    ## Run 85 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173399  max resid 0.05634067 
    ## Run 86 stress 0.09106059 
    ## Run 87 stress 0.07936953 
    ## ... Procrustes: rmse 0.009029257  max resid 0.03236468 
    ## Run 88 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045845  max resid 0.0324696 
    ## Run 89 stress 0.07951441 
    ## ... Procrustes: rmse 0.017359  max resid 0.05642683 
    ## Run 90 stress 0.08985604 
    ## Run 91 stress 0.07936951 
    ## ... Procrustes: rmse 0.009047697  max resid 0.0324806 
    ## Run 92 stress 0.07936951 
    ## ... Procrustes: rmse 0.009038771  max resid 0.03242825 
    ## Run 93 stress 0.0929627 
    ## Run 94 stress 0.07969262 
    ## ... Procrustes: rmse 0.01238873  max resid 0.04970065 
    ## Run 95 stress 0.07936953 
    ## ... Procrustes: rmse 0.009033163  max resid 0.03238965 
    ## Run 96 stress 0.0910603 
    ## Run 97 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173503  max resid 0.05638431 
    ## Run 98 stress 0.07936952 
    ## ... Procrustes: rmse 0.009051942  max resid 0.03249269 
    ## Run 99 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735245  max resid 0.05644457 
    ## Run 100 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734496  max resid 0.05637455 
    ## Run 101 stress 0.07927534 
    ## ... Procrustes: rmse 4.861973e-06  max resid 1.427254e-05 
    ## ... Similar to previous best
    ## Run 102 stress 0.08985603 
    ## Run 103 stress 0.07927535 
    ## ... Procrustes: rmse 3.944829e-05  max resid 0.0001090623 
    ## ... Similar to previous best
    ## Run 104 stress 0.07927536 
    ## ... Procrustes: rmse 5.72989e-05  max resid 0.00016217 
    ## ... Similar to previous best
    ## Run 105 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240373  max resid 0.04987417 
    ## Run 106 stress 0.0792754 
    ## ... Procrustes: rmse 9.83318e-05  max resid 0.0002702916 
    ## ... Similar to previous best
    ## Run 107 stress 0.0796926 
    ## ... Procrustes: rmse 0.01239026  max resid 0.04980977 
    ## Run 108 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044547  max resid 0.03246317 
    ## Run 109 stress 0.09043736 
    ## Run 110 stress 0.08985603 
    ## Run 111 stress 0.09072894 
    ## Run 112 stress 0.09247103 
    ## Run 113 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173636  max resid 0.05647629 
    ## Run 114 stress 0.09565368 
    ## Run 115 stress 0.09290693 
    ## Run 116 stress 0.09296215 
    ## Run 117 stress 0.09544395 
    ## Run 118 stress 0.07936953 
    ## ... Procrustes: rmse 0.009035786  max resid 0.03240808 
    ## Run 119 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734265  max resid 0.05631178 
    ## Run 120 stress 0.07969268 
    ## ... Procrustes: rmse 0.0123671  max resid 0.04956373 
    ## Run 121 stress 0.09252221 
    ## Run 122 stress 0.08985604 
    ## Run 123 stress 0.0796926 
    ## ... Procrustes: rmse 0.01236692  max resid 0.0496127 
    ## Run 124 stress 0.07927534 
    ## ... Procrustes: rmse 1.257429e-05  max resid 4.019855e-05 
    ## ... Similar to previous best
    ## Run 125 stress 0.09388354 
    ## Run 126 stress 0.09290717 
    ## Run 127 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237311  max resid 0.04967378 
    ## Run 128 stress 0.09290695 
    ## Run 129 stress 0.07927536 
    ## ... Procrustes: rmse 5.860345e-05  max resid 0.0001656856 
    ## ... Similar to previous best
    ## Run 130 stress 0.07969258 
    ## ... Procrustes: rmse 0.01238082  max resid 0.0497322 
    ## Run 131 stress 0.07927534 
    ## ... Procrustes: rmse 7.498143e-06  max resid 1.321172e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.0792754 
    ## ... Procrustes: rmse 9.723882e-05  max resid 0.0002782296 
    ## ... Similar to previous best
    ## Run 133 stress 0.07951444 
    ## ... Procrustes: rmse 0.01738285  max resid 0.0566158 
    ## Run 134 stress 0.07927535 
    ## ... Procrustes: rmse 2.83928e-05  max resid 7.985252e-05 
    ## ... Similar to previous best
    ## Run 135 stress 0.07936954 
    ## ... Procrustes: rmse 0.009030872  max resid 0.03237425 
    ## Run 136 stress 0.07927535 
    ## ... Procrustes: rmse 4.034287e-05  max resid 0.0001128553 
    ## ... Similar to previous best
    ## Run 137 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046129  max resid 0.03247443 
    ## Run 138 stress 0.07936952 
    ## ... Procrustes: rmse 0.009053447  max resid 0.03249982 
    ## Run 139 stress 0.09172783 
    ## Run 140 stress 0.09106048 
    ## Run 141 stress 0.07927555 
    ## ... Procrustes: rmse 4.531758e-05  max resid 0.0001069964 
    ## ... Similar to previous best
    ## Run 142 stress 0.08985605 
    ## Run 143 stress 0.07969258 
    ## ... Procrustes: rmse 0.01236768  max resid 0.04966084 
    ## Run 144 stress 0.07969269 
    ## ... Procrustes: rmse 0.01236357  max resid 0.04953146 
    ## Run 145 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735196  max resid 0.05642935 
    ## Run 146 stress 0.08985607 
    ## Run 147 stress 0.08985611 
    ## Run 148 stress 0.08985604 
    ## Run 149 stress 0.09022808 
    ## Run 150 stress 0.08985606 
    ## Run 151 stress 0.09290695 
    ## Run 152 stress 0.09043758 
    ## Run 153 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047635  max resid 0.03248313 
    ## Run 154 stress 0.08985608 
    ## Run 155 stress 0.09485767 
    ## Run 156 stress 0.07927538 
    ## ... Procrustes: rmse 6.342653e-05  max resid 0.0001773359 
    ## ... Similar to previous best
    ## Run 157 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045861  max resid 0.03246089 
    ## Run 158 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238783  max resid 0.0497828 
    ## Run 159 stress 0.08985603 
    ## Run 160 stress 0.09072905 
    ## Run 161 stress 0.09043781 
    ## Run 162 stress 0.0898562 
    ## Run 163 stress 0.09356327 
    ## Run 164 stress 0.0910603 
    ## Run 165 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237543  max resid 0.04970364 
    ## Run 166 stress 0.07927534 
    ## ... Procrustes: rmse 1.148858e-05  max resid 2.903113e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.07927538 
    ## ... Procrustes: rmse 8.174085e-05  max resid 0.00022535 
    ## ... Similar to previous best
    ## Run 168 stress 0.07969263 
    ## ... Procrustes: rmse 0.01241011  max resid 0.0499174 
    ## Run 169 stress 0.09022784 
    ## Run 170 stress 0.07969263 
    ## ... Procrustes: rmse 0.0123695  max resid 0.04959588 
    ## Run 171 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735023  max resid 0.05640366 
    ## Run 172 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735538  max resid 0.05644035 
    ## Run 173 stress 0.106141 
    ## Run 174 stress 0.07936955 
    ## ... Procrustes: rmse 0.009043592  max resid 0.03245983 
    ## Run 175 stress 0.07927539 
    ## ... Procrustes: rmse 7.372559e-05  max resid 0.000207658 
    ## ... Similar to previous best
    ## Run 176 stress 0.09290716 
    ## Run 177 stress 0.09072894 
    ## Run 178 stress 0.08985606 
    ## Run 179 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735071  max resid 0.05643908 
    ## Run 180 stress 0.07936949 
    ## ... Procrustes: rmse 0.00904481  max resid 0.03246571 
    ## Run 181 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049126  max resid 0.0324893 
    ## Run 182 stress 0.07951443 
    ## ... Procrustes: rmse 0.01734699  max resid 0.05638255 
    ## Run 183 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237885  max resid 0.04970445 
    ## Run 184 stress 0.09565387 
    ## Run 185 stress 0.08985605 
    ## Run 186 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046153  max resid 0.03247229 
    ## Run 187 stress 0.09247107 
    ## Run 188 stress 0.07936949 
    ## ... Procrustes: rmse 0.00904544  max resid 0.0324693 
    ## Run 189 stress 0.07927535 
    ## ... Procrustes: rmse 2.790941e-05  max resid 7.967973e-05 
    ## ... Similar to previous best
    ## Run 190 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237577  max resid 0.04970671 
    ## Run 191 stress 0.09072891 
    ## Run 192 stress 0.09290716 
    ## Run 193 stress 0.07927537 
    ## ... Procrustes: rmse 4.41849e-05  max resid 9.827205e-05 
    ## ... Similar to previous best
    ## Run 194 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237461  max resid 0.04969489 
    ## Run 195 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734419  max resid 0.05636926 
    ## Run 196 stress 0.09290706 
    ## Run 197 stress 0.09247104 
    ## Run 198 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735286  max resid 0.0564138 
    ## Run 199 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735394  max resid 0.05644132 
    ## Run 200 stress 0.1033494 
    ## Run 201 stress 0.09072899 
    ## Run 202 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237569  max resid 0.04971107 
    ## Run 203 stress 0.09382056 
    ## Run 204 stress 0.07927536 
    ## ... Procrustes: rmse 4.176487e-05  max resid 0.0001003603 
    ## ... Similar to previous best
    ## Run 205 stress 0.104507 
    ## Run 206 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237001  max resid 0.04963812 
    ## Run 207 stress 0.07969262 
    ## ... Procrustes: rmse 0.01240228  max resid 0.04987562 
    ## Run 208 stress 0.09043759 
    ## Run 209 stress 0.09565352 
    ## Run 210 stress 0.07927534 
    ## ... Procrustes: rmse 1.138509e-05  max resid 2.560091e-05 
    ## ... Similar to previous best
    ## Run 211 stress 0.07927535 
    ## ... Procrustes: rmse 3.566926e-05  max resid 9.972743e-05 
    ## ... Similar to previous best
    ## Run 212 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045682  max resid 0.03247275 
    ## Run 213 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735267  max resid 0.05642205 
    ## Run 214 stress 0.07927534 
    ## ... Procrustes: rmse 1.589745e-05  max resid 3.805309e-05 
    ## ... Similar to previous best
    ## Run 215 stress 0.07927536 
    ## ... Procrustes: rmse 5.15599e-05  max resid 0.0001372557 
    ## ... Similar to previous best
    ## Run 216 stress 0.09022784 
    ## Run 217 stress 0.08985603 
    ## Run 218 stress 0.09072892 
    ## Run 219 stress 0.08985612 
    ## Run 220 stress 0.08985606 
    ## Run 221 stress 0.07969263 
    ## ... Procrustes: rmse 0.01241511  max resid 0.04993308 
    ## Run 222 stress 0.09072891 
    ## Run 223 stress 0.08985605 
    ## Run 224 stress 0.07969264 
    ## ... Procrustes: rmse 0.01241157  max resid 0.04993869 
    ## Run 225 stress 0.07927538 
    ## ... Procrustes: rmse 8.267096e-05  max resid 0.0002358222 
    ## ... Similar to previous best
    ## Run 226 stress 0.07936953 
    ## ... Procrustes: rmse 0.009029789  max resid 0.03236776 
    ## Run 227 stress 0.09022779 
    ## Run 228 stress 0.07927534 
    ## ... Procrustes: rmse 1.768376e-05  max resid 4.644002e-05 
    ## ... Similar to previous best
    ## Run 229 stress 0.0796926 
    ## ... Procrustes: rmse 0.01239155  max resid 0.04980733 
    ## Run 230 stress 0.08985614 
    ## Run 231 stress 0.07927535 
    ## ... Procrustes: rmse 3.819375e-05  max resid 0.0001036413 
    ## ... Similar to previous best
    ## Run 232 stress 0.09388344 
    ## Run 233 stress 0.09072891 
    ## Run 234 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045433  max resid 0.03247061 
    ## Run 235 stress 0.1108444 
    ## Run 236 stress 0.1064573 
    ## Run 237 stress 0.07936951 
    ## ... Procrustes: rmse 0.009036096  max resid 0.03240807 
    ## Run 238 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 3.667033e-06  max resid 1.033809e-05 
    ## ... Similar to previous best
    ## Run 239 stress 0.1105014 
    ## Run 240 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733556  max resid 0.05630338 
    ## Run 241 stress 0.07927535 
    ## ... Procrustes: rmse 2.852842e-05  max resid 7.339877e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.07927534 
    ## ... Procrustes: rmse 2.174448e-05  max resid 5.893118e-05 
    ## ... Similar to previous best
    ## Run 243 stress 0.07969261 
    ## ... Procrustes: rmse 0.01241164  max resid 0.04986082 
    ## Run 244 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237972  max resid 0.0497399 
    ## Run 245 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237285  max resid 0.04966838 
    ## Run 246 stress 0.09043738 
    ## Run 247 stress 0.0956536 
    ## Run 248 stress 0.09072913 
    ## Run 249 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237493  max resid 0.04969484 
    ## Run 250 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237472  max resid 0.04973582 
    ## Run 251 stress 0.07927538 
    ## ... Procrustes: rmse 7.568687e-05  max resid 0.0001787447 
    ## ... Similar to previous best
    ## Run 252 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735036  max resid 0.05638383 
    ## Run 253 stress 0.08985614 
    ## Run 254 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734954  max resid 0.05640981 
    ## Run 255 stress 0.07936949 
    ## ... Procrustes: rmse 0.009044519  max resid 0.03246331 
    ## Run 256 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735333  max resid 0.05645098 
    ## Run 257 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240322  max resid 0.04987152 
    ## Run 258 stress 0.09485696 
    ## Run 259 stress 0.07927538 
    ## ... Procrustes: rmse 6.829388e-05  max resid 0.0001473893 
    ## ... Similar to previous best
    ## Run 260 stress 0.07936957 
    ## ... Procrustes: rmse 0.00902151  max resid 0.03231607 
    ## Run 261 stress 0.07951447 
    ## ... Procrustes: rmse 0.01731855  max resid 0.05617116 
    ## Run 262 stress 0.09474097 
    ## Run 263 stress 0.09022777 
    ## Run 264 stress 0.09106036 
    ## Run 265 stress 0.08985605 
    ## Run 266 stress 0.09388328 
    ## Run 267 stress 0.09072897 
    ## Run 268 stress 0.07951444 
    ## ... Procrustes: rmse 0.0173757  max resid 0.05657254 
    ## Run 269 stress 0.0796926 
    ## ... Procrustes: rmse 0.01239804  max resid 0.04982426 
    ## Run 270 stress 0.106141 
    ## Run 271 stress 0.09022791 
    ## Run 272 stress 0.09247104 
    ## Run 273 stress 0.08985606 
    ## Run 274 stress 0.09106032 
    ## Run 275 stress 0.07951442 
    ## ... Procrustes: rmse 0.01736119  max resid 0.05646977 
    ## Run 276 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734081  max resid 0.05634089 
    ## Run 277 stress 0.09106047 
    ## Run 278 stress 0.09388355 
    ## Run 279 stress 0.08985604 
    ## Run 280 stress 0.1083561 
    ## Run 281 stress 0.0796926 
    ## ... Procrustes: rmse 0.01239354  max resid 0.04981041 
    ## Run 282 stress 0.1044029 
    ## Run 283 stress 0.07927535 
    ## ... Procrustes: rmse 1.451336e-05  max resid 3.572265e-05 
    ## ... Similar to previous best
    ## Run 284 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047954  max resid 0.03248464 
    ## Run 285 stress 0.09043771 
    ## Run 286 stress 0.07969262 
    ## ... Procrustes: rmse 0.01237656  max resid 0.04975798 
    ## Run 287 stress 0.09388335 
    ## Run 288 stress 0.07927538 
    ## ... Procrustes: rmse 6.950126e-05  max resid 0.0001852332 
    ## ... Similar to previous best
    ## Run 289 stress 0.09106039 
    ## Run 290 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039401  max resid 0.03243377 
    ## Run 291 stress 0.09544394 
    ## Run 292 stress 0.09382065 
    ## Run 293 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045051  max resid 0.03246394 
    ## Run 294 stress 0.09290699 
    ## Run 295 stress 0.07936952 
    ## ... Procrustes: rmse 0.009031276  max resid 0.03237344 
    ## Run 296 stress 0.08985605 
    ## Run 297 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237726  max resid 0.04972782 
    ## Run 298 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735656  max resid 0.05644241 
    ## Run 299 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041944  max resid 0.03244927 
    ## Run 300 stress 0.07951445 
    ## ... Procrustes: rmse 0.01738025  max resid 0.05662857 
    ## Run 301 stress 0.07927538 
    ## ... Procrustes: rmse 4.924018e-05  max resid 0.0001058627 
    ## ... Similar to previous best
    ## Run 302 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173535  max resid 0.05642152 
    ## Run 303 stress 0.07927534 
    ## ... Procrustes: rmse 1.367557e-05  max resid 3.471951e-05 
    ## ... Similar to previous best
    ## Run 304 stress 0.09022776 
    ## Run 305 stress 0.1041926 
    ## Run 306 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044128  max resid 0.03245459 
    ## Run 307 stress 0.09290703 
    ## Run 308 stress 0.07936951 
    ## ... Procrustes: rmse 0.009040285  max resid 0.03243645 
    ## Run 309 stress 0.1046274 
    ## Run 310 stress 0.07951446 
    ## ... Procrustes: rmse 0.01732577  max resid 0.05622529 
    ## Run 311 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237541  max resid 0.04970488 
    ## Run 312 stress 0.0902279 
    ## Run 313 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039413  max resid 0.03242157 
    ## Run 314 stress 0.07927535 
    ## ... Procrustes: rmse 3.251168e-05  max resid 8.301345e-05 
    ## ... Similar to previous best
    ## Run 315 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044395  max resid 0.032457 
    ## Run 316 stress 0.09565351 
    ## Run 317 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041484  max resid 0.03244008 
    ## Run 318 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734468  max resid 0.05642481 
    ## Run 319 stress 0.08985603 
    ## Run 320 stress 0.09072894 
    ## Run 321 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735449  max resid 0.05643579 
    ## Run 322 stress 0.08985607 
    ## Run 323 stress 0.09106038 
    ## Run 324 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239724  max resid 0.04983684 
    ## Run 325 stress 0.0929072 
    ## Run 326 stress 0.09043739 
    ## Run 327 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238918  max resid 0.04979048 
    ## Run 328 stress 0.09565386 
    ## Run 329 stress 0.07927534 
    ## ... Procrustes: rmse 6.699674e-06  max resid 1.903099e-05 
    ## ... Similar to previous best
    ## Run 330 stress 0.09043737 
    ## Run 331 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046771  max resid 0.03247458 
    ## Run 332 stress 0.09022787 
    ## Run 333 stress 0.09072895 
    ## Run 334 stress 0.09252216 
    ## Run 335 stress 0.07969267 
    ## ... Procrustes: rmse 0.01235889  max resid 0.04950578 
    ## Run 336 stress 0.07951443 
    ## ... Procrustes: rmse 0.01733006  max resid 0.05627955 
    ## Run 337 stress 0.0793695 
    ## ... Procrustes: rmse 0.009054015  max resid 0.03251466 
    ## Run 338 stress 0.07936953 
    ## ... Procrustes: rmse 0.009029852  max resid 0.03236744 
    ## Run 339 stress 0.09172784 
    ## Run 340 stress 0.07927535 
    ## ... Procrustes: rmse 2.366785e-05  max resid 6.363552e-05 
    ## ... Similar to previous best
    ## Run 341 stress 0.1047992 
    ## Run 342 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734922  max resid 0.05639541 
    ## Run 343 stress 0.07951442 
    ## ... Procrustes: rmse 0.01732817  max resid 0.05626835 
    ## Run 344 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237053  max resid 0.04962347 
    ## Run 345 stress 0.09388361 
    ## Run 346 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043476  max resid 0.03244818 
    ## Run 347 stress 0.07936949 
    ## ... Procrustes: rmse 0.009048087  max resid 0.03248479 
    ## Run 348 stress 0.08985609 
    ## Run 349 stress 0.09474024 
    ## Run 350 stress 0.09474008 
    ## Run 351 stress 0.1097213 
    ## Run 352 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238803  max resid 0.04978057 
    ## Run 353 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733249  max resid 0.0563535 
    ## Run 354 stress 0.07969261 
    ## ... Procrustes: rmse 0.0123894  max resid 0.04977569 
    ## Run 355 stress 0.07936951 
    ## ... Procrustes: rmse 0.009041771  max resid 0.03244199 
    ## Run 356 stress 0.09072896 
    ## Run 357 stress 0.08985605 
    ## Run 358 stress 0.09344108 
    ## Run 359 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237864  max resid 0.04971916 
    ## Run 360 stress 0.08985613 
    ## Run 361 stress 0.07936952 
    ## ... Procrustes: rmse 0.009045583  max resid 0.03246854 
    ## Run 362 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 3.058991e-06  max resid 9.309266e-06 
    ## ... Similar to previous best
    ## Run 363 stress 0.09106035 
    ## Run 364 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735952  max resid 0.05645527 
    ## Run 365 stress 0.0898561 
    ## Run 366 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732371  max resid 0.05621961 
    ## Run 367 stress 0.08985611 
    ## Run 368 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734647  max resid 0.05636661 
    ## Run 369 stress 0.07936954 
    ## ... Procrustes: rmse 0.009030985  max resid 0.03237565 
    ## Run 370 stress 0.09252226 
    ## Run 371 stress 0.08985606 
    ## Run 372 stress 0.09106036 
    ## Run 373 stress 0.07936951 
    ## ... Procrustes: rmse 0.009041832  max resid 0.03244317 
    ## Run 374 stress 0.09043735 
    ## Run 375 stress 0.09072894 
    ## Run 376 stress 0.0793695 
    ## ... Procrustes: rmse 0.009051643  max resid 0.03250284 
    ## Run 377 stress 0.09290692 
    ## Run 378 stress 0.07927535 
    ## ... Procrustes: rmse 2.016125e-05  max resid 5.454281e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.08985606 
    ## Run 380 stress 0.07936952 
    ## ... Procrustes: rmse 0.009034111  max resid 0.03239566 
    ## Run 381 stress 0.09296245 
    ## Run 382 stress 0.07951444 
    ## ... Procrustes: rmse 0.0173615  max resid 0.05653656 
    ## Run 383 stress 0.09247104 
    ## Run 384 stress 0.07951443 
    ## ... Procrustes: rmse 0.01736147  max resid 0.05653905 
    ## Run 385 stress 0.09544375 
    ## Run 386 stress 0.09388336 
    ## Run 387 stress 0.09247104 
    ## Run 388 stress 0.3637375 
    ## Run 389 stress 0.08985604 
    ## Run 390 stress 0.07927534 
    ## ... Procrustes: rmse 2.225974e-05  max resid 6.119973e-05 
    ## ... Similar to previous best
    ## Run 391 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735282  max resid 0.05641431 
    ## Run 392 stress 0.09022783 
    ## Run 393 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733372  max resid 0.05629237 
    ## Run 394 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736791  max resid 0.05649613 
    ## Run 395 stress 0.07927534 
    ## ... Procrustes: rmse 1.007901e-05  max resid 2.026714e-05 
    ## ... Similar to previous best
    ## Run 396 stress 0.07969266 
    ## ... Procrustes: rmse 0.01236788  max resid 0.04957359 
    ## Run 397 stress 0.07969268 
    ## ... Procrustes: rmse 0.01236873  max resid 0.04954963 
    ## Run 398 stress 0.090228 
    ## Run 399 stress 0.09172784 
    ## Run 400 stress 0.09252217 
    ## Run 401 stress 0.07927535 
    ## ... Procrustes: rmse 2.62548e-05  max resid 5.835344e-05 
    ## ... Similar to previous best
    ## Run 402 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045122  max resid 0.032465 
    ## Run 403 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239735  max resid 0.04975734 
    ## Run 404 stress 0.07927535 
    ## ... Procrustes: rmse 1.147292e-05  max resid 2.641186e-05 
    ## ... Similar to previous best
    ## Run 405 stress 0.07936949 
    ## ... Procrustes: rmse 0.009044032  max resid 0.03245863 
    ## Run 406 stress 0.07936949 
    ## ... Procrustes: rmse 0.00904538  max resid 0.03246895 
    ## Run 407 stress 0.08985614 
    ## Run 408 stress 0.07927534 
    ## ... Procrustes: rmse 2.191844e-05  max resid 5.980683e-05 
    ## ... Similar to previous best
    ## Run 409 stress 0.09290722 
    ## Run 410 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237286  max resid 0.04964946 
    ## Run 411 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173597  max resid 0.05645767 
    ## Run 412 stress 0.07936952 
    ## ... Procrustes: rmse 0.009034702  max resid 0.03239957 
    ## Run 413 stress 0.09072904 
    ## Run 414 stress 0.07927534 
    ## ... Procrustes: rmse 1.241842e-05  max resid 3.491445e-05 
    ## ... Similar to previous best
    ## Run 415 stress 0.0793695 
    ## ... Procrustes: rmse 0.009048103  max resid 0.03248023 
    ## Run 416 stress 0.07927534 
    ## ... Procrustes: rmse 1.962705e-06  max resid 4.476242e-06 
    ## ... Similar to previous best
    ## Run 417 stress 0.07927534 
    ## ... Procrustes: rmse 1.630113e-05  max resid 3.837597e-05 
    ## ... Similar to previous best
    ## Run 418 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734603  max resid 0.05635547 
    ## Run 419 stress 0.09106027 
    ## Run 420 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237872  max resid 0.04971471 
    ## Run 421 stress 0.07969262 
    ## ... Procrustes: rmse 0.01237266  max resid 0.04972582 
    ## Run 422 stress 0.09252219 
    ## Run 423 stress 0.09252215 
    ## Run 424 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237193  max resid 0.04963865 
    ## Run 425 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736025  max resid 0.05647985 
    ## Run 426 stress 0.0934409 
    ## Run 427 stress 0.07927534 
    ## ... Procrustes: rmse 5.562173e-06  max resid 1.476708e-05 
    ## ... Similar to previous best
    ## Run 428 stress 0.09290716 
    ## Run 429 stress 0.0954437 
    ## Run 430 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045905  max resid 0.0324709 
    ## Run 431 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734826  max resid 0.05637337 
    ## Run 432 stress 0.1117676 
    ## Run 433 stress 0.09043735 
    ## Run 434 stress 0.106141 
    ## Run 435 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044609  max resid 0.032462 
    ## Run 436 stress 0.09356314 
    ## Run 437 stress 0.07927535 
    ## ... Procrustes: rmse 3.438456e-05  max resid 8.596772e-05 
    ## ... Similar to previous best
    ## Run 438 stress 0.07927534 
    ## ... Procrustes: rmse 7.957542e-06  max resid 1.419165e-05 
    ## ... Similar to previous best
    ## Run 439 stress 0.07936951 
    ## ... Procrustes: rmse 0.009051248  max resid 0.03248968 
    ## Run 440 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734188  max resid 0.05634407 
    ## Run 441 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736081  max resid 0.056466 
    ## Run 442 stress 0.09296222 
    ## Run 443 stress 0.08985608 
    ## Run 444 stress 0.09474053 
    ## Run 445 stress 0.09106029 
    ## Run 446 stress 0.09072896 
    ## Run 447 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239781  max resid 0.04984217 
    ## Run 448 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237496  max resid 0.04963734 
    ## Run 449 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173413  max resid 0.05631536 
    ## Run 450 stress 0.09106034 
    ## Run 451 stress 0.07951444 
    ## ... Procrustes: rmse 0.0173408  max resid 0.05633574 
    ## Run 452 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045514  max resid 0.03246867 
    ## Run 453 stress 0.0902279 
    ## Run 454 stress 0.08985606 
    ## Run 455 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735935  max resid 0.0564597 
    ## Run 456 stress 0.0898561 
    ## Run 457 stress 0.09290688 
    ## Run 458 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047353  max resid 0.03248031 
    ## Run 459 stress 0.07927535 
    ## ... Procrustes: rmse 3.200905e-05  max resid 8.894627e-05 
    ## ... Similar to previous best
    ## Run 460 stress 0.07927534 
    ## ... Procrustes: rmse 2.546742e-06  max resid 7.032294e-06 
    ## ... Similar to previous best
    ## Run 461 stress 0.07936951 
    ## ... Procrustes: rmse 0.009037871  max resid 0.03242289 
    ## Run 462 stress 0.07969267 
    ## ... Procrustes: rmse 0.01236884  max resid 0.04955625 
    ## Run 463 stress 0.08985607 
    ## Run 464 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 3.172667e-06  max resid 6.923583e-06 
    ## ... Similar to previous best
    ## Run 465 stress 0.07927536 
    ## ... Procrustes: rmse 4.529914e-05  max resid 0.0001251592 
    ## ... Similar to previous best
    ## Run 466 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734917  max resid 0.05639385 
    ## Run 467 stress 0.08985608 
    ## Run 468 stress 0.108875 
    ## Run 469 stress 0.09022778 
    ## Run 470 stress 0.07951445 
    ## ... Procrustes: rmse 0.01734846  max resid 0.05649105 
    ## Run 471 stress 0.07927534 
    ## ... Procrustes: rmse 2.860536e-06  max resid 6.49152e-06 
    ## ... Similar to previous best
    ## Run 472 stress 0.07927534 
    ## ... Procrustes: rmse 4.846387e-06  max resid 1.656729e-05 
    ## ... Similar to previous best
    ## Run 473 stress 0.07927534 
    ## ... Procrustes: rmse 3.545047e-06  max resid 1.109731e-05 
    ## ... Similar to previous best
    ## Run 474 stress 0.07936963 
    ## ... Procrustes: rmse 0.009047613  max resid 0.03247621 
    ## Run 475 stress 0.09106035 
    ## Run 476 stress 0.08985605 
    ## Run 477 stress 0.08985608 
    ## Run 478 stress 0.09172788 
    ## Run 479 stress 0.09022788 
    ## Run 480 stress 0.07927537 
    ## ... Procrustes: rmse 6.061376e-05  max resid 0.0001616788 
    ## ... Similar to previous best
    ## Run 481 stress 0.07927535 
    ## ... Procrustes: rmse 1.433178e-05  max resid 3.471261e-05 
    ## ... Similar to previous best
    ## Run 482 stress 0.07927535 
    ## ... Procrustes: rmse 2.868823e-05  max resid 7.362787e-05 
    ## ... Similar to previous best
    ## Run 483 stress 0.09106029 
    ## Run 484 stress 0.1108459 
    ## Run 485 stress 0.07951443 
    ## ... Procrustes: rmse 0.01735193  max resid 0.05644566 
    ## Run 486 stress 0.07927535 
    ## ... Procrustes: rmse 3.527429e-05  max resid 9.752775e-05 
    ## ... Similar to previous best
    ## Run 487 stress 0.09382055 
    ## Run 488 stress 0.07927535 
    ## ... Procrustes: rmse 1.06016e-05  max resid 2.479138e-05 
    ## ... Similar to previous best
    ## Run 489 stress 0.0793695 
    ## ... Procrustes: rmse 0.009048736  max resid 0.03248639 
    ## Run 490 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043058  max resid 0.03245061 
    ## Run 491 stress 0.07927536 
    ## ... Procrustes: rmse 4.780161e-05  max resid 0.0001187153 
    ## ... Similar to previous best
    ## Run 492 stress 0.09388347 
    ## Run 493 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734042  max resid 0.05637612 
    ## Run 494 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238439  max resid 0.0497305 
    ## Run 495 stress 0.07927534 
    ## ... Procrustes: rmse 1.102196e-05  max resid 3.078141e-05 
    ## ... Similar to previous best
    ## Run 496 stress 0.07936951 
    ## ... Procrustes: rmse 0.009050007  max resid 0.03249349 
    ## Run 497 stress 0.09043767 
    ## Run 498 stress 0.07927535 
    ## ... Procrustes: rmse 4.497847e-05  max resid 0.0001203157 
    ## ... Similar to previous best
    ## Run 499 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733782  max resid 0.05632038 
    ## Run 500 stress 0.0793695 
    ## ... Procrustes: rmse 0.009037993  max resid 0.0324212 
    ## *** Best solution repeated 13 times

``` r
# Mixed and stratified lakes
PD_beta_geo_MS_NMDS <- metaMDS(PD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.05968634 
    ## Run 2 stress 0.054667 
    ## Run 3 stress 0.06027439 
    ## Run 4 stress 0.06136676 
    ## Run 5 stress 0.05051336 
    ## ... Procrustes: rmse 0.03578581  max resid 0.1221097 
    ## Run 6 stress 0.05171058 
    ## Run 7 stress 0.05036798 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002882235  max resid 0.0007343884 
    ## ... Similar to previous best
    ## Run 8 stress 0.05051319 
    ## ... Procrustes: rmse 0.03566987  max resid 0.1218623 
    ## Run 9 stress 0.05036795 
    ## ... New best solution
    ## ... Procrustes: rmse 4.231699e-05  max resid 0.0001050754 
    ## ... Similar to previous best
    ## Run 10 stress 0.05829939 
    ## Run 11 stress 0.05036797 
    ## ... Procrustes: rmse 3.707951e-05  max resid 8.660985e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.05968627 
    ## Run 13 stress 0.0517106 
    ## Run 14 stress 0.05403567 
    ## Run 15 stress 0.05036795 
    ## ... Procrustes: rmse 3.642396e-05  max resid 6.030471e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.05403566 
    ## Run 17 stress 0.05171062 
    ## Run 18 stress 0.06136692 
    ## Run 19 stress 0.05466707 
    ## Run 20 stress 0.05829937 
    ## Run 21 stress 0.06136677 
    ## Run 22 stress 0.050368 
    ## ... Procrustes: rmse 5.582244e-05  max resid 0.0001377112 
    ## ... Similar to previous best
    ## Run 23 stress 0.05036803 
    ## ... Procrustes: rmse 8.179903e-05  max resid 0.0002296078 
    ## ... Similar to previous best
    ## Run 24 stress 0.05036799 
    ## ... Procrustes: rmse 0.0002251823  max resid 0.000587124 
    ## ... Similar to previous best
    ## Run 25 stress 0.05036799 
    ## ... Procrustes: rmse 0.0002150197  max resid 0.0005949204 
    ## ... Similar to previous best
    ## Run 26 stress 0.06027428 
    ## Run 27 stress 0.06231307 
    ## Run 28 stress 0.05900152 
    ## Run 29 stress 0.06647173 
    ## Run 30 stress 0.05051334 
    ## ... Procrustes: rmse 0.03566649  max resid 0.1216893 
    ## Run 31 stress 0.06640352 
    ## Run 32 stress 0.05036793 
    ## ... New best solution
    ## ... Procrustes: rmse 2.592598e-05  max resid 6.450103e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.3657332 
    ## Run 34 stress 0.050368 
    ## ... Procrustes: rmse 0.0002076969  max resid 0.0005541372 
    ## ... Similar to previous best
    ## Run 35 stress 0.05171061 
    ## Run 36 stress 0.05051313 
    ## ... Procrustes: rmse 0.03566831  max resid 0.1217219 
    ## Run 37 stress 0.05051345 
    ## ... Procrustes: rmse 0.03564787  max resid 0.1216355 
    ## Run 38 stress 0.06655703 
    ## Run 39 stress 0.0590015 
    ## Run 40 stress 0.05036797 
    ## ... Procrustes: rmse 5.186325e-05  max resid 0.000126888 
    ## ... Similar to previous best
    ## Run 41 stress 0.06309863 
    ## Run 42 stress 0.05171065 
    ## Run 43 stress 0.05051274 
    ## ... Procrustes: rmse 0.03561428  max resid 0.1216566 
    ## Run 44 stress 0.05602044 
    ## Run 45 stress 0.06362779 
    ## Run 46 stress 0.05466701 
    ## Run 47 stress 0.05036801 
    ## ... Procrustes: rmse 0.0002088144  max resid 0.0005486076 
    ## ... Similar to previous best
    ## Run 48 stress 0.05403599 
    ## Run 49 stress 0.06609058 
    ## Run 50 stress 0.05466703 
    ## Run 51 stress 0.06640338 
    ## Run 52 stress 0.05337581 
    ## Run 53 stress 0.05602053 
    ## Run 54 stress 0.05928222 
    ## Run 55 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001591812  max resid 0.0004210044 
    ## ... Similar to previous best
    ## Run 56 stress 0.0517106 
    ## Run 57 stress 0.05036794 
    ## ... Procrustes: rmse 2.027869e-05  max resid 4.372387e-05 
    ## ... Similar to previous best
    ## Run 58 stress 0.05051295 
    ## ... Procrustes: rmse 0.03565871  max resid 0.1217141 
    ## Run 59 stress 0.06478536 
    ## Run 60 stress 0.05337581 
    ## Run 61 stress 0.05036795 
    ## ... Procrustes: rmse 3.451537e-05  max resid 8.401518e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.0505132 
    ## ... Procrustes: rmse 0.03568488  max resid 0.1219181 
    ## Run 63 stress 0.05051351 
    ## ... Procrustes: rmse 0.03570763  max resid 0.1219977 
    ## Run 64 stress 0.06704955 
    ## Run 65 stress 0.06027438 
    ## Run 66 stress 0.05466707 
    ## Run 67 stress 0.05051296 
    ## ... Procrustes: rmse 0.03566648  max resid 0.1218462 
    ## Run 68 stress 0.05036792 
    ## ... New best solution
    ## ... Procrustes: rmse 8.437097e-06  max resid 1.664182e-05 
    ## ... Similar to previous best
    ## Run 69 stress 0.05036794 
    ## ... Procrustes: rmse 2.834126e-05  max resid 5.530705e-05 
    ## ... Similar to previous best
    ## Run 70 stress 0.0505133 
    ## ... Procrustes: rmse 0.03570329  max resid 0.1219795 
    ## Run 71 stress 0.05036795 
    ## ... Procrustes: rmse 4.158894e-05  max resid 8.415548e-05 
    ## ... Similar to previous best
    ## Run 72 stress 0.0671676 
    ## Run 73 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001091815  max resid 0.0002731927 
    ## ... Similar to previous best
    ## Run 74 stress 0.3380554 
    ## Run 75 stress 0.06743837 
    ## Run 76 stress 0.05403584 
    ## Run 77 stress 0.05036796 
    ## ... Procrustes: rmse 5.098639e-05  max resid 0.000131557 
    ## ... Similar to previous best
    ## Run 78 stress 0.0660906 
    ## Run 79 stress 0.05337579 
    ## Run 80 stress 0.050368 
    ## ... Procrustes: rmse 0.0002021766  max resid 0.0005268793 
    ## ... Similar to previous best
    ## Run 81 stress 0.05968626 
    ## Run 82 stress 0.05036798 
    ## ... Procrustes: rmse 7.036755e-05  max resid 0.0001682552 
    ## ... Similar to previous best
    ## Run 83 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001680748  max resid 0.0004317372 
    ## ... Similar to previous best
    ## Run 84 stress 0.05171058 
    ## Run 85 stress 0.05051324 
    ## ... Procrustes: rmse 0.03570919  max resid 0.1218361 
    ## Run 86 stress 0.05051355 
    ## ... Procrustes: rmse 0.03567702  max resid 0.1219102 
    ## Run 87 stress 0.06640341 
    ## Run 88 stress 0.05829934 
    ## Run 89 stress 0.05928252 
    ## Run 90 stress 0.05403601 
    ## Run 91 stress 0.05829948 
    ## Run 92 stress 0.05602047 
    ## Run 93 stress 0.06851853 
    ## Run 94 stress 0.06231287 
    ## Run 95 stress 0.05036803 
    ## ... Procrustes: rmse 0.0002278681  max resid 0.0005934171 
    ## ... Similar to previous best
    ## Run 96 stress 0.05829937 
    ## Run 97 stress 0.05036806 
    ## ... Procrustes: rmse 0.0002471658  max resid 0.0006383972 
    ## ... Similar to previous best
    ## Run 98 stress 0.06136681 
    ## Run 99 stress 0.05036798 
    ## ... Procrustes: rmse 7.484037e-05  max resid 0.0001847654 
    ## ... Similar to previous best
    ## Run 100 stress 0.0517106 
    ## Run 101 stress 0.3462972 
    ## Run 102 stress 0.06136677 
    ## Run 103 stress 0.050368 
    ## ... Procrustes: rmse 7.465486e-05  max resid 0.000183847 
    ## ... Similar to previous best
    ## Run 104 stress 0.05171058 
    ## Run 105 stress 0.05051336 
    ## ... Procrustes: rmse 0.03567082  max resid 0.1218851 
    ## Run 106 stress 0.05036795 
    ## ... Procrustes: rmse 3.995978e-05  max resid 9.881008e-05 
    ## ... Similar to previous best
    ## Run 107 stress 0.05036797 
    ## ... Procrustes: rmse 6.144769e-05  max resid 0.0001529019 
    ## ... Similar to previous best
    ## Run 108 stress 0.05466704 
    ## Run 109 stress 0.05036811 
    ## ... Procrustes: rmse 0.0001477782  max resid 0.0003774332 
    ## ... Similar to previous best
    ## Run 110 stress 0.05051322 
    ## ... Procrustes: rmse 0.03568042  max resid 0.121751 
    ## Run 111 stress 0.06111765 
    ## Run 112 stress 0.05171063 
    ## Run 113 stress 0.05466702 
    ## Run 114 stress 0.05036793 
    ## ... Procrustes: rmse 0.0001116158  max resid 0.0002751223 
    ## ... Similar to previous best
    ## Run 115 stress 0.05171058 
    ## Run 116 stress 0.05171059 
    ## Run 117 stress 0.05171059 
    ## Run 118 stress 0.05051344 
    ## ... Procrustes: rmse 0.03564285  max resid 0.1216194 
    ## Run 119 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001061483  max resid 0.0002719092 
    ## ... Similar to previous best
    ## Run 120 stress 0.05171058 
    ## Run 121 stress 0.05171058 
    ## Run 122 stress 0.06231288 
    ## Run 123 stress 0.05051324 
    ## ... Procrustes: rmse 0.03562559  max resid 0.121589 
    ## Run 124 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001174987  max resid 0.0002966773 
    ## ... Similar to previous best
    ## Run 125 stress 0.05968643 
    ## Run 126 stress 0.05036796 
    ## ... Procrustes: rmse 4.568729e-05  max resid 0.0001099701 
    ## ... Similar to previous best
    ## Run 127 stress 0.05036795 
    ## ... Procrustes: rmse 4.225104e-05  max resid 0.0001036046 
    ## ... Similar to previous best
    ## Run 128 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001537098  max resid 0.0003903745 
    ## ... Similar to previous best
    ## Run 129 stress 0.05036802 
    ## ... Procrustes: rmse 8.925336e-05  max resid 0.000220066 
    ## ... Similar to previous best
    ## Run 130 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001949207  max resid 0.0005123188 
    ## ... Similar to previous best
    ## Run 131 stress 0.05466705 
    ## Run 132 stress 0.05036797 
    ## ... Procrustes: rmse 6.776478e-05  max resid 0.0001795249 
    ## ... Similar to previous best
    ## Run 133 stress 0.0533758 
    ## Run 134 stress 0.05171058 
    ## Run 135 stress 0.05036814 
    ## ... Procrustes: rmse 0.0001496593  max resid 0.0003860946 
    ## ... Similar to previous best
    ## Run 136 stress 0.05051357 
    ## ... Procrustes: rmse 0.03572861  max resid 0.1220646 
    ## Run 137 stress 0.05337596 
    ## Run 138 stress 0.05968631 
    ## Run 139 stress 0.05051305 
    ## ... Procrustes: rmse 0.03566746  max resid 0.1218584 
    ## Run 140 stress 0.05036805 
    ## ... Procrustes: rmse 0.0002363182  max resid 0.0006281779 
    ## ... Similar to previous best
    ## Run 141 stress 0.050368 
    ## ... Procrustes: rmse 0.0002010205  max resid 0.0005098096 
    ## ... Similar to previous best
    ## Run 142 stress 0.06136682 
    ## Run 143 stress 0.05051348 
    ## ... Procrustes: rmse 0.03568236  max resid 0.1217364 
    ## Run 144 stress 0.06038958 
    ## Run 145 stress 0.05036795 
    ## ... Procrustes: rmse 0.000138437  max resid 0.0003488599 
    ## ... Similar to previous best
    ## Run 146 stress 0.05337585 
    ## Run 147 stress 0.05051317 
    ## ... Procrustes: rmse 0.03563635  max resid 0.1216245 
    ## Run 148 stress 0.05171059 
    ## Run 149 stress 0.05051281 
    ## ... Procrustes: rmse 0.03566673  max resid 0.1217657 
    ## Run 150 stress 0.05036792 
    ## ... New best solution
    ## ... Procrustes: rmse 1.933153e-05  max resid 4.555293e-05 
    ## ... Similar to previous best
    ## Run 151 stress 0.05171062 
    ## Run 152 stress 0.05036796 
    ## ... Procrustes: rmse 7.184235e-05  max resid 0.0001579512 
    ## ... Similar to previous best
    ## Run 153 stress 0.05928232 
    ## Run 154 stress 0.05051355 
    ## ... Procrustes: rmse 0.03565155  max resid 0.1218346 
    ## Run 155 stress 0.05051323 
    ## ... Procrustes: rmse 0.03567538  max resid 0.1218944 
    ## Run 156 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001163671  max resid 0.0002905615 
    ## ... Similar to previous best
    ## Run 157 stress 0.0560207 
    ## Run 158 stress 0.0560205 
    ## Run 159 stress 0.05466702 
    ## Run 160 stress 0.06038949 
    ## Run 161 stress 0.05979227 
    ## Run 162 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 6.344592e-05  max resid 0.00014598 
    ## ... Similar to previous best
    ## Run 163 stress 0.05051316 
    ## ... Procrustes: rmse 0.03571584  max resid 0.1220255 
    ## Run 164 stress 0.0505128 
    ## ... Procrustes: rmse 0.03569332  max resid 0.1218947 
    ## Run 165 stress 0.05900147 
    ## Run 166 stress 0.0505135 
    ## ... Procrustes: rmse 0.03575042  max resid 0.1221429 
    ## Run 167 stress 0.05036796 
    ## ... Procrustes: rmse 8.267894e-05  max resid 0.0002318398 
    ## ... Similar to previous best
    ## Run 168 stress 0.0582993 
    ## Run 169 stress 0.06647191 
    ## Run 170 stress 0.05036792 
    ## ... Procrustes: rmse 3.714936e-05  max resid 7.896959e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.0517106 
    ## Run 172 stress 0.07272426 
    ## Run 173 stress 0.05900149 
    ## Run 174 stress 0.05051334 
    ## ... Procrustes: rmse 0.03578035  max resid 0.1222255 
    ## Run 175 stress 0.06231315 
    ## Run 176 stress 0.06844313 
    ## Run 177 stress 0.0505129 
    ## ... Procrustes: rmse 0.03570079  max resid 0.1219662 
    ## Run 178 stress 0.06674573 
    ## Run 179 stress 0.05968629 
    ## Run 180 stress 0.05171063 
    ## Run 181 stress 0.05466707 
    ## Run 182 stress 0.05829931 
    ## Run 183 stress 0.06231302 
    ## Run 184 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001172341  max resid 0.0002804014 
    ## ... Similar to previous best
    ## Run 185 stress 0.06478528 
    ## Run 186 stress 0.05403556 
    ## Run 187 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001245509  max resid 0.0003500379 
    ## ... Similar to previous best
    ## Run 188 stress 0.05051316 
    ## ... Procrustes: rmse 0.03571874  max resid 0.1220351 
    ## Run 189 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001850213  max resid 0.0004208365 
    ## ... Similar to previous best
    ## Run 190 stress 0.05403585 
    ## Run 191 stress 0.06136676 
    ## Run 192 stress 0.05602054 
    ## Run 193 stress 0.05171058 
    ## Run 194 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001564838  max resid 0.0003837463 
    ## ... Similar to previous best
    ## Run 195 stress 0.06894572 
    ## Run 196 stress 0.05036793 
    ## ... Procrustes: rmse 3.907302e-05  max resid 8.199077e-05 
    ## ... Similar to previous best
    ## Run 197 stress 0.06309852 
    ## Run 198 stress 0.06027436 
    ## Run 199 stress 0.05051339 
    ## ... Procrustes: rmse 0.03574563  max resid 0.1221262 
    ## Run 200 stress 0.0533758 
    ## Run 201 stress 0.05171058 
    ## Run 202 stress 0.05466702 
    ## Run 203 stress 0.05602054 
    ## Run 204 stress 0.050513 
    ## ... Procrustes: rmse 0.03571087  max resid 0.1220044 
    ## Run 205 stress 0.0582993 
    ## Run 206 stress 0.05979228 
    ## Run 207 stress 0.06478538 
    ## Run 208 stress 0.05900159 
    ## Run 209 stress 0.06027437 
    ## Run 210 stress 0.06826724 
    ## Run 211 stress 0.05051291 
    ## ... Procrustes: rmse 0.03571251  max resid 0.1220026 
    ## Run 212 stress 0.05968629 
    ## Run 213 stress 0.0517106 
    ## Run 214 stress 0.05403579 
    ## Run 215 stress 0.06111771 
    ## Run 216 stress 0.05051331 
    ## ... Procrustes: rmse 0.03569356  max resid 0.1219674 
    ## Run 217 stress 0.05051314 
    ## ... Procrustes: rmse 0.03572236  max resid 0.1220463 
    ## Run 218 stress 0.05051343 
    ## ... Procrustes: rmse 0.03567498  max resid 0.1217365 
    ## Run 219 stress 0.0505133 
    ## ... Procrustes: rmse 0.03574523  max resid 0.1221216 
    ## Run 220 stress 0.05036795 
    ## ... Procrustes: rmse 0.000112259  max resid 0.0003147257 
    ## ... Similar to previous best
    ## Run 221 stress 0.05036799 
    ## ... Procrustes: rmse 6.583447e-05  max resid 0.0001711398 
    ## ... Similar to previous best
    ## Run 222 stress 0.05051321 
    ## ... Procrustes: rmse 0.03570143  max resid 0.1218374 
    ## Run 223 stress 0.05051311 
    ## ... Procrustes: rmse 0.03574089  max resid 0.1219674 
    ## Run 224 stress 0.05171059 
    ## Run 225 stress 0.05602048 
    ## Run 226 stress 0.0636279 
    ## Run 227 stress 0.05051341 
    ## ... Procrustes: rmse 0.03574386  max resid 0.1221213 
    ## Run 228 stress 0.0517106 
    ## Run 229 stress 0.05171061 
    ## Run 230 stress 0.05602072 
    ## Run 231 stress 0.05171064 
    ## Run 232 stress 0.05051294 
    ## ... Procrustes: rmse 0.0357078  max resid 0.1219909 
    ## Run 233 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 7.097338e-06  max resid 1.380113e-05 
    ## ... Similar to previous best
    ## Run 234 stress 0.05051313 
    ## ... Procrustes: rmse 0.03568899  max resid 0.1218073 
    ## Run 235 stress 0.06499242 
    ## Run 236 stress 0.06309838 
    ## Run 237 stress 0.06478534 
    ## Run 238 stress 0.06851848 
    ## Run 239 stress 0.05171058 
    ## Run 240 stress 0.0647852 
    ## Run 241 stress 0.05602048 
    ## Run 242 stress 0.06111768 
    ## Run 243 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001520182  max resid 0.0003719102 
    ## ... Similar to previous best
    ## Run 244 stress 0.05466701 
    ## Run 245 stress 0.05051325 
    ## ... Procrustes: rmse 0.03572277  max resid 0.1220508 
    ## Run 246 stress 0.06027436 
    ## Run 247 stress 0.05036798 
    ## ... Procrustes: rmse 6.952773e-05  max resid 0.00017943 
    ## ... Similar to previous best
    ## Run 248 stress 0.05051295 
    ## ... Procrustes: rmse 0.03569083  max resid 0.121833 
    ## Run 249 stress 0.06231288 
    ## Run 250 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001191317  max resid 0.0002920164 
    ## ... Similar to previous best
    ## Run 251 stress 0.054667 
    ## Run 252 stress 0.05337586 
    ## Run 253 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001161634  max resid 0.0002646763 
    ## ... Similar to previous best
    ## Run 254 stress 0.05051324 
    ## ... Procrustes: rmse 0.03570711  max resid 0.121849 
    ## Run 255 stress 0.05051325 
    ## ... Procrustes: rmse 0.03573093  max resid 0.1220752 
    ## Run 256 stress 0.0560205 
    ## Run 257 stress 0.06111768 
    ## Run 258 stress 0.06027439 
    ## Run 259 stress 0.05829929 
    ## Run 260 stress 0.06309836 
    ## Run 261 stress 0.0560205 
    ## Run 262 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001251918  max resid 0.0003511345 
    ## ... Similar to previous best
    ## Run 263 stress 0.05968639 
    ## Run 264 stress 0.05403582 
    ## Run 265 stress 0.05968631 
    ## Run 266 stress 0.05051335 
    ## ... Procrustes: rmse 0.03573875  max resid 0.1221024 
    ## Run 267 stress 0.05928214 
    ## Run 268 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001299448  max resid 0.0003042333 
    ## ... Similar to previous best
    ## Run 269 stress 0.05036796 
    ## ... Procrustes: rmse 7.909128e-05  max resid 0.000205675 
    ## ... Similar to previous best
    ## Run 270 stress 0.05051305 
    ## ... Procrustes: rmse 0.03568355  max resid 0.121799 
    ## Run 271 stress 0.05051315 
    ## ... Procrustes: rmse 0.03576373  max resid 0.1220316 
    ## Run 272 stress 0.05337576 
    ## Run 273 stress 0.05036793 
    ## ... Procrustes: rmse 9.965332e-05  max resid 0.0002418689 
    ## ... Similar to previous best
    ## Run 274 stress 0.05171058 
    ## Run 275 stress 0.05337584 
    ## Run 276 stress 0.05051305 
    ## ... Procrustes: rmse 0.03576511  max resid 0.1220554 
    ## Run 277 stress 0.05051331 
    ## ... Procrustes: rmse 0.03568584  max resid 0.1217776 
    ## Run 278 stress 0.05051314 
    ## ... Procrustes: rmse 0.03567373  max resid 0.1217599 
    ## Run 279 stress 0.0682671 
    ## Run 280 stress 0.05829939 
    ## Run 281 stress 0.05036795 
    ## ... Procrustes: rmse 7.134299e-05  max resid 0.0001988797 
    ## ... Similar to previous best
    ## Run 282 stress 0.05602068 
    ## Run 283 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001459832  max resid 0.0003541914 
    ## ... Similar to previous best
    ## Run 284 stress 0.06309836 
    ## Run 285 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001542594  max resid 0.0003745954 
    ## ... Similar to previous best
    ## Run 286 stress 0.06136676 
    ## Run 287 stress 0.05979227 
    ## Run 288 stress 0.05036793 
    ## ... Procrustes: rmse 4.436614e-05  max resid 0.0001263073 
    ## ... Similar to previous best
    ## Run 289 stress 0.05051333 
    ## ... Procrustes: rmse 0.0356897  max resid 0.1217876 
    ## Run 290 stress 0.06136678 
    ## Run 291 stress 0.05036808 
    ## ... Procrustes: rmse 0.0001995518  max resid 0.0004873518 
    ## ... Similar to previous best
    ## Run 292 stress 0.05051338 
    ## ... Procrustes: rmse 0.03569886  max resid 0.12181 
    ## Run 293 stress 0.05036804 
    ## ... Procrustes: rmse 0.0002011557  max resid 0.0004635478 
    ## ... Similar to previous best
    ## Run 294 stress 0.0694384 
    ## Run 295 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001307126  max resid 0.0003616618 
    ## ... Similar to previous best
    ## Run 296 stress 0.05602063 
    ## Run 297 stress 0.05968629 
    ## Run 298 stress 0.05337578 
    ## Run 299 stress 0.05403603 
    ## Run 300 stress 0.05171059 
    ## Run 301 stress 0.05928219 
    ## Run 302 stress 0.06038955 
    ## Run 303 stress 0.05171058 
    ## Run 304 stress 0.0664036 
    ## Run 305 stress 0.05051329 
    ## ... Procrustes: rmse 0.03573116  max resid 0.1220771 
    ## Run 306 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001054606  max resid 0.0002980572 
    ## ... Similar to previous best
    ## Run 307 stress 0.05829934 
    ## Run 308 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001613657  max resid 0.0003975487 
    ## ... Similar to previous best
    ## Run 309 stress 0.05968641 
    ## Run 310 stress 0.05968636 
    ## Run 311 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001768144  max resid 0.0004370325 
    ## ... Similar to previous best
    ## Run 312 stress 0.05171059 
    ## Run 313 stress 0.05968629 
    ## Run 314 stress 0.0505133 
    ## ... Procrustes: rmse 0.03579358  max resid 0.1221104 
    ## Run 315 stress 0.05051345 
    ## ... Procrustes: rmse 0.03573737  max resid 0.1221013 
    ## Run 316 stress 0.0647854 
    ## Run 317 stress 0.05051327 
    ## ... Procrustes: rmse 0.03572911  max resid 0.12207 
    ## Run 318 stress 0.07217903 
    ## Run 319 stress 0.05051296 
    ## ... Procrustes: rmse 0.03569713  max resid 0.1218512 
    ## Run 320 stress 0.05928242 
    ## Run 321 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001096905  max resid 0.0003175051 
    ## ... Similar to previous best
    ## Run 322 stress 0.05602061 
    ## Run 323 stress 0.07221091 
    ## Run 324 stress 0.05051334 
    ## ... Procrustes: rmse 0.03573827  max resid 0.1221004 
    ## Run 325 stress 0.05036793 
    ## ... Procrustes: rmse 5.482338e-05  max resid 0.0001284844 
    ## ... Similar to previous best
    ## Run 326 stress 0.05171058 
    ## Run 327 stress 0.06499233 
    ## Run 328 stress 0.05337588 
    ## Run 329 stress 0.05466702 
    ## Run 330 stress 0.05036792 
    ## ... Procrustes: rmse 4.864711e-05  max resid 0.0001148282 
    ## ... Similar to previous best
    ## Run 331 stress 0.0636279 
    ## Run 332 stress 0.05051343 
    ## ... Procrustes: rmse 0.03570877  max resid 0.1220145 
    ## Run 333 stress 0.06609055 
    ## Run 334 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001174703  max resid 0.0002748487 
    ## ... Similar to previous best
    ## Run 335 stress 0.05171063 
    ## Run 336 stress 0.05968633 
    ## Run 337 stress 0.06789814 
    ## Run 338 stress 0.05968629 
    ## Run 339 stress 0.05928246 
    ## Run 340 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001413574  max resid 0.0003342442 
    ## ... Similar to previous best
    ## Run 341 stress 0.05051285 
    ## ... Procrustes: rmse 0.03569702  max resid 0.1219487 
    ## Run 342 stress 0.05337585 
    ## Run 343 stress 0.05968626 
    ## Run 344 stress 0.06504227 
    ## Run 345 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001034543  max resid 0.0002546842 
    ## ... Similar to previous best
    ## Run 346 stress 0.05051345 
    ## ... Procrustes: rmse 0.03572604  max resid 0.1218892 
    ## Run 347 stress 0.05051316 
    ## ... Procrustes: rmse 0.03572648  max resid 0.1220556 
    ## Run 348 stress 0.05051305 
    ## ... Procrustes: rmse 0.03569459  max resid 0.1218354 
    ## Run 349 stress 0.06478525 
    ## Run 350 stress 0.06478515 
    ## Run 351 stress 0.06027426 
    ## Run 352 stress 0.05051313 
    ## ... Procrustes: rmse 0.0357056  max resid 0.1218559 
    ## Run 353 stress 0.05051306 
    ## ... Procrustes: rmse 0.03575408  max resid 0.1220127 
    ## Run 354 stress 0.05337578 
    ## Run 355 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001409915  max resid 0.0003454169 
    ## ... Similar to previous best
    ## Run 356 stress 0.05928228 
    ## Run 357 stress 0.05602048 
    ## Run 358 stress 0.05900153 
    ## Run 359 stress 0.0623129 
    ## Run 360 stress 0.0560205 
    ## Run 361 stress 0.06231294 
    ## Run 362 stress 0.06231313 
    ## Run 363 stress 0.06027442 
    ## Run 364 stress 0.05036792 
    ## ... Procrustes: rmse 8.021452e-05  max resid 0.0001935285 
    ## ... Similar to previous best
    ## Run 365 stress 0.05171059 
    ## Run 366 stress 0.3658211 
    ## Run 367 stress 0.3563931 
    ## Run 368 stress 0.06851844 
    ## Run 369 stress 0.0546671 
    ## Run 370 stress 0.05171065 
    ## Run 371 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001866992  max resid 0.0004597906 
    ## ... Similar to previous best
    ## Run 372 stress 0.06231303 
    ## Run 373 stress 0.05036795 
    ## ... Procrustes: rmse 9.81308e-05  max resid 0.0002302082 
    ## ... Similar to previous best
    ## Run 374 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001365627  max resid 0.0003158364 
    ## ... Similar to previous best
    ## Run 375 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001087068  max resid 0.0003171307 
    ## ... Similar to previous best
    ## Run 376 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001127664  max resid 0.0003055275 
    ## ... Similar to previous best
    ## Run 377 stress 0.06111766 
    ## Run 378 stress 0.05171058 
    ## Run 379 stress 0.05829931 
    ## Run 380 stress 0.05051314 
    ## ... Procrustes: rmse 0.03572412  max resid 0.12205 
    ## Run 381 stress 0.05036794 
    ## ... Procrustes: rmse 6.438246e-05  max resid 0.0001663576 
    ## ... Similar to previous best
    ## Run 382 stress 0.06111778 
    ## Run 383 stress 0.05171062 
    ## Run 384 stress 0.05968628 
    ## Run 385 stress 0.06027429 
    ## Run 386 stress 0.054667 
    ## Run 387 stress 0.0613669 
    ## Run 388 stress 0.05602051 
    ## Run 389 stress 0.05036793 
    ## ... Procrustes: rmse 7.600973e-05  max resid 0.000183386 
    ## ... Similar to previous best
    ## Run 390 stress 0.05051346 
    ## ... Procrustes: rmse 0.03573025  max resid 0.1219001 
    ## Run 391 stress 0.06362778 
    ## Run 392 stress 0.06674551 
    ## Run 393 stress 0.05900151 
    ## Run 394 stress 0.06111771 
    ## Run 395 stress 0.05171059 
    ## Run 396 stress 0.05928211 
    ## Run 397 stress 0.05036794 
    ## ... Procrustes: rmse 6.42182e-05  max resid 0.0001679523 
    ## ... Similar to previous best
    ## Run 398 stress 0.05051292 
    ## ... Procrustes: rmse 0.03577809  max resid 0.1221922 
    ## Run 399 stress 0.05403607 
    ## Run 400 stress 0.05051327 
    ## ... Procrustes: rmse 0.03570909  max resid 0.1218524 
    ## Run 401 stress 0.05829935 
    ## Run 402 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001470471  max resid 0.0003568423 
    ## ... Similar to previous best
    ## Run 403 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001372328  max resid 0.0003370683 
    ## ... Similar to previous best
    ## Run 404 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001254707  max resid 0.0003119197 
    ## ... Similar to previous best
    ## Run 405 stress 0.05171068 
    ## Run 406 stress 0.05051309 
    ## ... Procrustes: rmse 0.03570083  max resid 0.1219773 
    ## Run 407 stress 0.05900149 
    ## Run 408 stress 0.05900148 
    ## Run 409 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001334606  max resid 0.0003613751 
    ## ... Similar to previous best
    ## Run 410 stress 0.0582993 
    ## Run 411 stress 0.05051313 
    ## ... Procrustes: rmse 0.03570717  max resid 0.1219992 
    ## Run 412 stress 0.05051333 
    ## ... Procrustes: rmse 0.03562656  max resid 0.121765 
    ## Run 413 stress 0.05602062 
    ## Run 414 stress 0.05928218 
    ## Run 415 stress 0.05171058 
    ## Run 416 stress 0.05602055 
    ## Run 417 stress 0.05928234 
    ## Run 418 stress 0.05171058 
    ## Run 419 stress 0.06896927 
    ## Run 420 stress 0.06478545 
    ## Run 421 stress 0.05051302 
    ## ... Procrustes: rmse 0.03571125  max resid 0.1220046 
    ## Run 422 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001386154  max resid 0.0003440828 
    ## ... Similar to previous best
    ## Run 423 stress 0.05051339 
    ## ... Procrustes: rmse 0.03574066  max resid 0.1221096 
    ## Run 424 stress 0.06743841 
    ## Run 425 stress 0.05968632 
    ## Run 426 stress 0.06111767 
    ## Run 427 stress 0.05051331 
    ## ... Procrustes: rmse 0.0356426  max resid 0.1218148 
    ## Run 428 stress 0.05406543 
    ## Run 429 stress 0.06743852 
    ## Run 430 stress 0.05036797 
    ## ... Procrustes: rmse 6.577844e-05  max resid 0.0001627921 
    ## ... Similar to previous best
    ## Run 431 stress 0.05466706 
    ## Run 432 stress 0.05829943 
    ## Run 433 stress 0.05036792 
    ## ... Procrustes: rmse 5.578696e-05  max resid 0.0001349729 
    ## ... Similar to previous best
    ## Run 434 stress 0.06478538 
    ## Run 435 stress 0.05337587 
    ## Run 436 stress 0.05171061 
    ## Run 437 stress 0.05036792 
    ## ... Procrustes: rmse 6.704884e-05  max resid 0.000161176 
    ## ... Similar to previous best
    ## Run 438 stress 0.06844299 
    ## Run 439 stress 0.05928233 
    ## Run 440 stress 0.06362791 
    ## Run 441 stress 0.05036794 
    ## ... Procrustes: rmse 6.019912e-05  max resid 0.0001518187 
    ## ... Similar to previous best
    ## Run 442 stress 0.05337587 
    ## Run 443 stress 0.0694379 
    ## Run 444 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.563283e-05  max resid 3.207981e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.05602061 
    ## Run 446 stress 0.06951026 
    ## Run 447 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001409601  max resid 0.0003477117 
    ## ... Similar to previous best
    ## Run 448 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001356653  max resid 0.0003415058 
    ## ... Similar to previous best
    ## Run 449 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001442355  max resid 0.000362296 
    ## ... Similar to previous best
    ## Run 450 stress 0.05036797 
    ## ... Procrustes: rmse 9.913695e-05  max resid 0.000289346 
    ## ... Similar to previous best
    ## Run 451 stress 0.06674556 
    ## Run 452 stress 0.06896892 
    ## Run 453 stress 0.05403591 
    ## Run 454 stress 0.06640336 
    ## Run 455 stress 0.0533758 
    ## Run 456 stress 0.05051346 
    ## ... Procrustes: rmse 0.03573841  max resid 0.1221012 
    ## Run 457 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001497405  max resid 0.000394166 
    ## ... Similar to previous best
    ## Run 458 stress 0.06640334 
    ## Run 459 stress 0.050513 
    ## ... Procrustes: rmse 0.03570476  max resid 0.1219811 
    ## Run 460 stress 0.05171058 
    ## Run 461 stress 0.06309845 
    ## Run 462 stress 0.05036799 
    ## ... Procrustes: rmse 9.120126e-05  max resid 0.0002290153 
    ## ... Similar to previous best
    ## Run 463 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001060687  max resid 0.0002601608 
    ## ... Similar to previous best
    ## Run 464 stress 0.05036809 
    ## ... Procrustes: rmse 0.0002176904  max resid 0.0005492137 
    ## ... Similar to previous best
    ## Run 465 stress 0.050368 
    ## ... Procrustes: rmse 0.0001371413  max resid 0.0003735201 
    ## ... Similar to previous best
    ## Run 466 stress 0.05979236 
    ## Run 467 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001072595  max resid 0.0002712611 
    ## ... Similar to previous best
    ## Run 468 stress 0.06730253 
    ## Run 469 stress 0.06362797 
    ## Run 470 stress 0.05466701 
    ## Run 471 stress 0.3549658 
    ## Run 472 stress 0.05829932 
    ## Run 473 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001221901  max resid 0.0003067639 
    ## ... Similar to previous best
    ## Run 474 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001465304  max resid 0.0003568721 
    ## ... Similar to previous best
    ## Run 475 stress 0.06904524 
    ## Run 476 stress 0.05036792 
    ## ... Procrustes: rmse 4.769744e-05  max resid 0.0001218014 
    ## ... Similar to previous best
    ## Run 477 stress 0.06493696 
    ## Run 478 stress 0.06478521 
    ## Run 479 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001073204  max resid 0.0002674886 
    ## ... Similar to previous best
    ## Run 480 stress 0.05968641 
    ## Run 481 stress 0.05900154 
    ## Run 482 stress 0.05036792 
    ## ... Procrustes: rmse 5.204143e-05  max resid 0.0001300146 
    ## ... Similar to previous best
    ## Run 483 stress 0.07211138 
    ## Run 484 stress 0.05051325 
    ## ... Procrustes: rmse 0.03567211  max resid 0.1217412 
    ## Run 485 stress 0.05968634 
    ## Run 486 stress 0.05171063 
    ## Run 487 stress 0.05466706 
    ## Run 488 stress 0.05171064 
    ## Run 489 stress 0.0689688 
    ## Run 490 stress 0.0560206 
    ## Run 491 stress 0.06231299 
    ## Run 492 stress 0.06309854 
    ## Run 493 stress 0.06309843 
    ## Run 494 stress 0.06478527 
    ## Run 495 stress 0.05051334 
    ## ... Procrustes: rmse 0.03572441  max resid 0.1220566 
    ## Run 496 stress 0.05337578 
    ## Run 497 stress 0.050513 
    ## ... Procrustes: rmse 0.03571118  max resid 0.1218842 
    ## Run 498 stress 0.05051321 
    ## ... Procrustes: rmse 0.03566683  max resid 0.1217347 
    ## Run 499 stress 0.05171062 
    ## Run 500 stress 0.05602063 
    ## *** Best solution repeated 16 times

``` r
# Ocean sites and mixed lakes
PD_beta_geo_OM_NMDS <- metaMDS(PD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1016168 
    ## Run 1 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0006420421  max resid 0.001797906 
    ## ... Similar to previous best
    ## Run 2 stress 0.1016167 
    ## ... Procrustes: rmse 0.000576117  max resid 0.001615366 
    ## ... Similar to previous best
    ## Run 3 stress 0.1341868 
    ## Run 4 stress 0.1341868 
    ## Run 5 stress 0.1016166 
    ## ... Procrustes: rmse 0.0004895504  max resid 0.001372236 
    ## ... Similar to previous best
    ## Run 6 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 8.870729e-05  max resid 0.0002477039 
    ## ... Similar to previous best
    ## Run 7 stress 0.1341868 
    ## Run 8 stress 0.1016164 
    ## ... Procrustes: rmse 1.033161e-05  max resid 2.898301e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.1356149 
    ## Run 10 stress 0.1530429 
    ## Run 11 stress 0.1016164 
    ## ... Procrustes: rmse 3.850179e-05  max resid 0.0001073892 
    ## ... Similar to previous best
    ## Run 12 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006178623  max resid 0.00172774 
    ## ... Similar to previous best
    ## Run 13 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002141332  max resid 0.000599918 
    ## ... Similar to previous best
    ## Run 14 stress 0.1016164 
    ## ... Procrustes: rmse 1.081889e-05  max resid 2.72046e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.1356144 
    ## Run 16 stress 0.1016164 
    ## ... Procrustes: rmse 7.858791e-05  max resid 0.0002219646 
    ## ... Similar to previous best
    ## Run 17 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 3.741368e-06  max resid 9.990511e-06 
    ## ... Similar to previous best
    ## Run 18 stress 0.1530424 
    ## Run 19 stress 0.1016164 
    ## ... Procrustes: rmse 1.15896e-05  max resid 3.117971e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.1341868 
    ## Run 21 stress 0.1016164 
    ## ... Procrustes: rmse 6.714375e-05  max resid 0.0001882715 
    ## ... Similar to previous best
    ## Run 22 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005948472  max resid 0.001660509 
    ## ... Similar to previous best
    ## Run 23 stress 0.1341868 
    ## Run 24 stress 0.1341868 
    ## Run 25 stress 0.1356149 
    ## Run 26 stress 0.1341868 
    ## Run 27 stress 0.1016164 
    ## ... Procrustes: rmse 5.145323e-05  max resid 0.0001432574 
    ## ... Similar to previous best
    ## Run 28 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004566838  max resid 0.001275106 
    ## ... Similar to previous best
    ## Run 29 stress 0.1016164 
    ## ... Procrustes: rmse 4.452667e-05  max resid 0.0001246263 
    ## ... Similar to previous best
    ## Run 30 stress 0.1016164 
    ## ... Procrustes: rmse 2.055692e-05  max resid 5.556863e-05 
    ## ... Similar to previous best
    ## Run 31 stress 0.135615 
    ## Run 32 stress 0.1016164 
    ## ... Procrustes: rmse 2.328626e-05  max resid 6.588717e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.1341868 
    ## Run 34 stress 0.1016164 
    ## ... Procrustes: rmse 2.61096e-05  max resid 7.414884e-05 
    ## ... Similar to previous best
    ## Run 35 stress 0.1341868 
    ## Run 36 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001197634  max resid 0.0003327584 
    ## ... Similar to previous best
    ## Run 37 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007354815  max resid 0.002051363 
    ## ... Similar to previous best
    ## Run 38 stress 0.1356145 
    ## Run 39 stress 0.1016164 
    ## ... Procrustes: rmse 2.801742e-05  max resid 7.998359e-05 
    ## ... Similar to previous best
    ## Run 40 stress 0.3193072 
    ## Run 41 stress 0.1016164 
    ## ... Procrustes: rmse 2.750558e-05  max resid 7.738528e-05 
    ## ... Similar to previous best
    ## Run 42 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004282199  max resid 0.001194273 
    ## ... Similar to previous best
    ## Run 43 stress 0.1016164 
    ## ... Procrustes: rmse 3.693122e-05  max resid 0.0001049813 
    ## ... Similar to previous best
    ## Run 44 stress 0.1341868 
    ## Run 45 stress 0.1016164 
    ## ... Procrustes: rmse 8.472302e-05  max resid 0.0002395374 
    ## ... Similar to previous best
    ## Run 46 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006212067  max resid 0.001735807 
    ## ... Similar to previous best
    ## Run 47 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005285712  max resid 0.001476371 
    ## ... Similar to previous best
    ## Run 48 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004683596  max resid 0.001309149 
    ## ... Similar to previous best
    ## Run 49 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002946095  max resid 0.0008241835 
    ## ... Similar to previous best
    ## Run 50 stress 0.1016164 
    ## ... Procrustes: rmse 6.494846e-05  max resid 0.0001836187 
    ## ... Similar to previous best
    ## Run 51 stress 0.1341868 
    ## Run 52 stress 0.1016164 
    ## ... Procrustes: rmse 3.555227e-05  max resid 9.878831e-05 
    ## ... Similar to previous best
    ## Run 53 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002353238  max resid 0.0006586364 
    ## ... Similar to previous best
    ## Run 54 stress 0.1016164 
    ## ... Procrustes: rmse 1.827562e-05  max resid 5.244462e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.1016164 
    ## ... Procrustes: rmse 2.402582e-06  max resid 6.510302e-06 
    ## ... Similar to previous best
    ## Run 56 stress 0.1016164 
    ## ... Procrustes: rmse 8.147723e-06  max resid 2.114213e-05 
    ## ... Similar to previous best
    ## Run 57 stress 0.1341868 
    ## Run 58 stress 0.1341868 
    ## Run 59 stress 0.1341868 
    ## Run 60 stress 0.1016164 
    ## ... Procrustes: rmse 3.712538e-05  max resid 0.0001051878 
    ## ... Similar to previous best
    ## Run 61 stress 0.1016164 
    ## ... Procrustes: rmse 8.904592e-06  max resid 2.335439e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001935716  max resid 0.0005401352 
    ## ... Similar to previous best
    ## Run 63 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002185208  max resid 0.0006097428 
    ## ... Similar to previous best
    ## Run 64 stress 0.1519033 
    ## Run 65 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005751132  max resid 0.001605917 
    ## ... Similar to previous best
    ## Run 66 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004294719  max resid 0.001196385 
    ## ... Similar to previous best
    ## Run 67 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003801961  max resid 0.001063672 
    ## ... Similar to previous best
    ## Run 68 stress 0.1016164 
    ## ... Procrustes: rmse 2.311739e-05  max resid 6.604629e-05 
    ## ... Similar to previous best
    ## Run 69 stress 0.1016164 
    ## ... Procrustes: rmse 1.554633e-05  max resid 4.051329e-05 
    ## ... Similar to previous best
    ## Run 70 stress 0.1016164 
    ## ... Procrustes: rmse 1.150473e-05  max resid 3.112632e-05 
    ## ... Similar to previous best
    ## Run 71 stress 0.101617 
    ## ... Procrustes: rmse 0.000692217  max resid 0.001933699 
    ## ... Similar to previous best
    ## Run 72 stress 0.1016164 
    ## ... Procrustes: rmse 6.496145e-05  max resid 0.0001842179 
    ## ... Similar to previous best
    ## Run 73 stress 0.1016164 
    ## ... Procrustes: rmse 7.05541e-05  max resid 0.0001997983 
    ## ... Similar to previous best
    ## Run 74 stress 0.1016164 
    ## ... Procrustes: rmse 4.853243e-06  max resid 1.412891e-05 
    ## ... Similar to previous best
    ## Run 75 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004227176  max resid 0.001180959 
    ## ... Similar to previous best
    ## Run 76 stress 0.1016164 
    ## ... Procrustes: rmse 7.261717e-05  max resid 0.0002053619 
    ## ... Similar to previous best
    ## Run 77 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003818101  max resid 0.001067524 
    ## ... Similar to previous best
    ## Run 78 stress 0.1016164 
    ## ... Procrustes: rmse 2.130588e-05  max resid 5.889625e-05 
    ## ... Similar to previous best
    ## Run 79 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003603775  max resid 0.001007271 
    ## ... Similar to previous best
    ## Run 80 stress 0.1016164 
    ## ... Procrustes: rmse 6.067713e-06  max resid 1.56407e-05 
    ## ... Similar to previous best
    ## Run 81 stress 0.1016164 
    ## ... Procrustes: rmse 2.101694e-06  max resid 5.474451e-06 
    ## ... Similar to previous best
    ## Run 82 stress 0.1356152 
    ## Run 83 stress 0.1016164 
    ## ... Procrustes: rmse 8.545877e-05  max resid 0.0002407234 
    ## ... Similar to previous best
    ## Run 84 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002879843  max resid 0.000806682 
    ## ... Similar to previous best
    ## Run 85 stress 0.1016164 
    ## ... Procrustes: rmse 4.912784e-05  max resid 0.0001396133 
    ## ... Similar to previous best
    ## Run 86 stress 0.1016164 
    ## ... Procrustes: rmse 2.496063e-05  max resid 7.134753e-05 
    ## ... Similar to previous best
    ## Run 87 stress 0.1016164 
    ## ... Procrustes: rmse 8.596816e-06  max resid 2.483677e-05 
    ## ... Similar to previous best
    ## Run 88 stress 0.1341868 
    ## Run 89 stress 0.302167 
    ## Run 90 stress 0.1016164 
    ## ... Procrustes: rmse 1.569545e-05  max resid 4.519716e-05 
    ## ... Similar to previous best
    ## Run 91 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003390457  max resid 0.0009465778 
    ## ... Similar to previous best
    ## Run 92 stress 0.1016164 
    ## ... Procrustes: rmse 9.118385e-05  max resid 0.000258177 
    ## ... Similar to previous best
    ## Run 93 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002812442  max resid 0.0007867872 
    ## ... Similar to previous best
    ## Run 94 stress 0.1341868 
    ## Run 95 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001900289  max resid 0.0005330287 
    ## ... Similar to previous best
    ## Run 96 stress 0.1341868 
    ## Run 97 stress 0.1341868 
    ## Run 98 stress 0.135615 
    ## Run 99 stress 0.212628 
    ## Run 100 stress 0.1341868 
    ## Run 101 stress 0.1341868 
    ## Run 102 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004886238  max resid 0.00136547 
    ## ... Similar to previous best
    ## Run 103 stress 0.1016164 
    ## ... Procrustes: rmse 7.111884e-06  max resid 1.907301e-05 
    ## ... Similar to previous best
    ## Run 104 stress 0.1341868 
    ## Run 105 stress 0.1016164 
    ## ... Procrustes: rmse 2.983596e-05  max resid 8.285562e-05 
    ## ... Similar to previous best
    ## Run 106 stress 0.1537555 
    ## Run 107 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004470452  max resid 0.001231505 
    ## ... Similar to previous best
    ## Run 108 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003860175  max resid 0.0010788 
    ## ... Similar to previous best
    ## Run 109 stress 0.1016164 
    ## ... Procrustes: rmse 0.000138945  max resid 0.0003927344 
    ## ... Similar to previous best
    ## Run 110 stress 0.1530422 
    ## Run 111 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002622454  max resid 0.0007338566 
    ## ... Similar to previous best
    ## Run 112 stress 0.1016164 
    ## ... Procrustes: rmse 4.048847e-05  max resid 0.00011505 
    ## ... Similar to previous best
    ## Run 113 stress 0.1016164 
    ## ... Procrustes: rmse 2.001793e-05  max resid 5.568127e-05 
    ## ... Similar to previous best
    ## Run 114 stress 0.101617 
    ## ... Procrustes: rmse 0.0006339483  max resid 0.001767184 
    ## ... Similar to previous best
    ## Run 115 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001495914  max resid 0.0004214925 
    ## ... Similar to previous best
    ## Run 116 stress 0.1016164 
    ## ... Procrustes: rmse 1.166464e-05  max resid 3.369475e-05 
    ## ... Similar to previous best
    ## Run 117 stress 0.1016164 
    ## ... Procrustes: rmse 6.850745e-05  max resid 0.0001938611 
    ## ... Similar to previous best
    ## Run 118 stress 0.2972302 
    ## Run 119 stress 0.1356144 
    ## Run 120 stress 0.1016164 
    ## ... Procrustes: rmse 5.99514e-05  max resid 0.0001690944 
    ## ... Similar to previous best
    ## Run 121 stress 0.1016164 
    ## ... Procrustes: rmse 1.257641e-05  max resid 3.522666e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.1016164 
    ## ... Procrustes: rmse 3.807529e-06  max resid 1.09384e-05 
    ## ... Similar to previous best
    ## Run 123 stress 0.3495561 
    ## Run 124 stress 0.1356152 
    ## Run 125 stress 0.1016164 
    ## ... Procrustes: rmse 1.875498e-05  max resid 5.380248e-05 
    ## ... Similar to previous best
    ## Run 126 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004167915  max resid 0.001163484 
    ## ... Similar to previous best
    ## Run 127 stress 0.1016164 
    ## ... Procrustes: rmse 1.852618e-05  max resid 5.246079e-05 
    ## ... Similar to previous best
    ## Run 128 stress 0.1356148 
    ## Run 129 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001047416  max resid 0.0002965516 
    ## ... Similar to previous best
    ## Run 130 stress 0.1016164 
    ## ... Procrustes: rmse 9.160934e-05  max resid 0.0002598502 
    ## ... Similar to previous best
    ## Run 131 stress 0.1356148 
    ## Run 132 stress 0.2984513 
    ## Run 133 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002721186  max resid 0.0007617935 
    ## ... Similar to previous best
    ## Run 134 stress 0.1016164 
    ## ... Procrustes: rmse 2.426494e-05  max resid 6.942099e-05 
    ## ... Similar to previous best
    ## Run 135 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005960637  max resid 0.001664533 
    ## ... Similar to previous best
    ## Run 136 stress 0.1016164 
    ## ... Procrustes: rmse 6.47806e-05  max resid 0.0001837778 
    ## ... Similar to previous best
    ## Run 137 stress 0.1341868 
    ## Run 138 stress 0.1016164 
    ## ... Procrustes: rmse 3.690727e-05  max resid 0.0001052183 
    ## ... Similar to previous best
    ## Run 139 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005082649  max resid 0.001420286 
    ## ... Similar to previous best
    ## Run 140 stress 0.1530423 
    ## Run 141 stress 0.1341868 
    ## Run 142 stress 0.1016169 
    ## ... Procrustes: rmse 0.000576332  max resid 0.001609587 
    ## ... Similar to previous best
    ## Run 143 stress 0.1341868 
    ## Run 144 stress 0.1341868 
    ## Run 145 stress 0.1356152 
    ## Run 146 stress 0.1016164 
    ## ... Procrustes: rmse 9.22655e-05  max resid 0.0002610981 
    ## ... Similar to previous best
    ## Run 147 stress 0.1341868 
    ## Run 148 stress 0.1356144 
    ## Run 149 stress 0.1356148 
    ## Run 150 stress 0.1356148 
    ## Run 151 stress 0.1016164 
    ## ... Procrustes: rmse 7.428389e-05  max resid 0.0002104877 
    ## ... Similar to previous best
    ## Run 152 stress 0.1341868 
    ## Run 153 stress 0.1016164 
    ## ... Procrustes: rmse 1.713571e-05  max resid 4.832496e-05 
    ## ... Similar to previous best
    ## Run 154 stress 0.1537556 
    ## Run 155 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005068953  max resid 0.001416134 
    ## ... Similar to previous best
    ## Run 156 stress 0.1016164 
    ## ... Procrustes: rmse 4.757128e-05  max resid 0.0001351368 
    ## ... Similar to previous best
    ## Run 157 stress 0.1356148 
    ## Run 158 stress 0.1016164 
    ## ... Procrustes: rmse 7.246519e-06  max resid 2.04706e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.1530428 
    ## Run 160 stress 0.1016164 
    ## ... Procrustes: rmse 9.06906e-06  max resid 2.632185e-05 
    ## ... Similar to previous best
    ## Run 161 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005145855  max resid 0.001437062 
    ## ... Similar to previous best
    ## Run 162 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004582781  max resid 0.001280021 
    ## ... Similar to previous best
    ## Run 163 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003718254  max resid 0.001037862 
    ## ... Similar to previous best
    ## Run 164 stress 0.1016164 
    ## ... Procrustes: rmse 4.949011e-06  max resid 1.432563e-05 
    ## ... Similar to previous best
    ## Run 165 stress 0.101617 
    ## ... Procrustes: rmse 0.0006459345  max resid 0.001804126 
    ## ... Similar to previous best
    ## Run 166 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005020391  max resid 0.001402323 
    ## ... Similar to previous best
    ## Run 167 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005124432  max resid 0.001431002 
    ## ... Similar to previous best
    ## Run 168 stress 0.1016164 
    ## ... Procrustes: rmse 2.779586e-05  max resid 7.959469e-05 
    ## ... Similar to previous best
    ## Run 169 stress 0.349553 
    ## Run 170 stress 0.1016164 
    ## ... Procrustes: rmse 2.322861e-05  max resid 6.421915e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.1016164 
    ## ... Procrustes: rmse 6.109051e-05  max resid 0.0001734214 
    ## ... Similar to previous best
    ## Run 172 stress 0.1016164 
    ## ... Procrustes: rmse 2.048309e-05  max resid 5.876232e-05 
    ## ... Similar to previous best
    ## Run 173 stress 0.1016164 
    ## ... Procrustes: rmse 7.914135e-06  max resid 2.299952e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.135615 
    ## Run 175 stress 0.1016164 
    ## ... Procrustes: rmse 4.824077e-05  max resid 0.0001323079 
    ## ... Similar to previous best
    ## Run 176 stress 0.1016164 
    ## ... Procrustes: rmse 8.030933e-05  max resid 0.0002238966 
    ## ... Similar to previous best
    ## Run 177 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002767409  max resid 0.0007746391 
    ## ... Similar to previous best
    ## Run 178 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005721746  max resid 0.001596237 
    ## ... Similar to previous best
    ## Run 179 stress 0.1016164 
    ## ... Procrustes: rmse 4.621202e-05  max resid 0.0001313675 
    ## ... Similar to previous best
    ## Run 180 stress 0.1016164 
    ## ... Procrustes: rmse 2.241222e-05  max resid 6.323082e-05 
    ## ... Similar to previous best
    ## Run 181 stress 0.1356149 
    ## Run 182 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001116948  max resid 0.0003162057 
    ## ... Similar to previous best
    ## Run 183 stress 0.1341868 
    ## Run 184 stress 0.1356146 
    ## Run 185 stress 0.1016164 
    ## ... Procrustes: rmse 2.339318e-05  max resid 6.695089e-05 
    ## ... Similar to previous best
    ## Run 186 stress 0.1356145 
    ## Run 187 stress 0.1016164 
    ## ... Procrustes: rmse 6.460706e-05  max resid 0.0001831676 
    ## ... Similar to previous best
    ## Run 188 stress 0.1016164 
    ## ... Procrustes: rmse 3.310864e-05  max resid 9.435836e-05 
    ## ... Similar to previous best
    ## Run 189 stress 0.2109396 
    ## Run 190 stress 0.1016164 
    ## ... Procrustes: rmse 2.472915e-05  max resid 6.968345e-05 
    ## ... Similar to previous best
    ## Run 191 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003353872  max resid 0.0009367398 
    ## ... Similar to previous best
    ## Run 192 stress 0.1016164 
    ## ... Procrustes: rmse 2.397902e-05  max resid 6.854306e-05 
    ## ... Similar to previous best
    ## Run 193 stress 0.1016164 
    ## ... Procrustes: rmse 3.66564e-05  max resid 0.0001035708 
    ## ... Similar to previous best
    ## Run 194 stress 0.1016164 
    ## ... Procrustes: rmse 5.343279e-05  max resid 0.0001479729 
    ## ... Similar to previous best
    ## Run 195 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004596592  max resid 0.001284707 
    ## ... Similar to previous best
    ## Run 196 stress 0.1016164 
    ## ... Procrustes: rmse 4.529725e-05  max resid 0.0001286131 
    ## ... Similar to previous best
    ## Run 197 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005889077  max resid 0.001645574 
    ## ... Similar to previous best
    ## Run 198 stress 0.1356148 
    ## Run 199 stress 0.1341868 
    ## Run 200 stress 0.1016164 
    ## ... Procrustes: rmse 1.877514e-05  max resid 4.321491e-05 
    ## ... Similar to previous best
    ## Run 201 stress 0.1016164 
    ## ... Procrustes: rmse 5.471277e-05  max resid 0.0001520819 
    ## ... Similar to previous best
    ## Run 202 stress 0.1016164 
    ## ... Procrustes: rmse 2.904202e-05  max resid 8.283765e-05 
    ## ... Similar to previous best
    ## Run 203 stress 0.1356149 
    ## Run 204 stress 0.1016164 
    ## ... Procrustes: rmse 3.390073e-05  max resid 9.662443e-05 
    ## ... Similar to previous best
    ## Run 205 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003512791  max resid 0.0009715596 
    ## ... Similar to previous best
    ## Run 206 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004203259  max resid 0.001174447 
    ## ... Similar to previous best
    ## Run 207 stress 0.1016164 
    ## ... Procrustes: rmse 6.290229e-05  max resid 0.0001784887 
    ## ... Similar to previous best
    ## Run 208 stress 0.1341868 
    ## Run 209 stress 0.1016164 
    ## ... Procrustes: rmse 7.578536e-05  max resid 0.0002139211 
    ## ... Similar to previous best
    ## Run 210 stress 0.1016164 
    ## ... Procrustes: rmse 4.38827e-05  max resid 0.000124265 
    ## ... Similar to previous best
    ## Run 211 stress 0.1016164 
    ## ... Procrustes: rmse 8.844395e-05  max resid 0.0002495753 
    ## ... Similar to previous best
    ## Run 212 stress 0.1356145 
    ## Run 213 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006100409  max resid 0.001702981 
    ## ... Similar to previous best
    ## Run 214 stress 0.1341868 
    ## Run 215 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004681909  max resid 0.001308756 
    ## ... Similar to previous best
    ## Run 216 stress 0.1016164 
    ## ... Procrustes: rmse 5.272563e-05  max resid 0.0001448658 
    ## ... Similar to previous best
    ## Run 217 stress 0.1519033 
    ## Run 218 stress 0.1016164 
    ## ... Procrustes: rmse 6.109515e-05  max resid 0.000173048 
    ## ... Similar to previous best
    ## Run 219 stress 0.1016164 
    ## ... Procrustes: rmse 5.939273e-05  max resid 0.0001686226 
    ## ... Similar to previous best
    ## Run 220 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001921842  max resid 0.0005376358 
    ## ... Similar to previous best
    ## Run 221 stress 0.1016164 
    ## ... Procrustes: rmse 9.53714e-06  max resid 2.759978e-05 
    ## ... Similar to previous best
    ## Run 222 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004734667  max resid 0.001324227 
    ## ... Similar to previous best
    ## Run 223 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001440236  max resid 0.0004065455 
    ## ... Similar to previous best
    ## Run 224 stress 0.1016164 
    ## ... Procrustes: rmse 2.98719e-05  max resid 8.522824e-05 
    ## ... Similar to previous best
    ## Run 225 stress 0.1016164 
    ## ... Procrustes: rmse 5.029703e-05  max resid 0.0001430895 
    ## ... Similar to previous best
    ## Run 226 stress 0.101617 
    ## ... Procrustes: rmse 0.0006536745  max resid 0.001825313 
    ## ... Similar to previous best
    ## Run 227 stress 0.1016164 
    ## ... Procrustes: rmse 3.386931e-05  max resid 9.655755e-05 
    ## ... Similar to previous best
    ## Run 228 stress 0.1016164 
    ## ... Procrustes: rmse 6.595367e-06  max resid 1.365787e-05 
    ## ... Similar to previous best
    ## Run 229 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005572009  max resid 0.001555204 
    ## ... Similar to previous best
    ## Run 230 stress 0.1341868 
    ## Run 231 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004358952  max resid 0.001212914 
    ## ... Similar to previous best
    ## Run 232 stress 0.1356147 
    ## Run 233 stress 0.1341868 
    ## Run 234 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002891178  max resid 0.0008090948 
    ## ... Similar to previous best
    ## Run 235 stress 0.1016164 
    ## ... Procrustes: rmse 9.284035e-06  max resid 2.65463e-05 
    ## ... Similar to previous best
    ## Run 236 stress 0.1016164 
    ## ... Procrustes: rmse 3.144576e-05  max resid 8.871821e-05 
    ## ... Similar to previous best
    ## Run 237 stress 0.1341868 
    ## Run 238 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001016092  max resid 0.0002876282 
    ## ... Similar to previous best
    ## Run 239 stress 0.135615 
    ## Run 240 stress 0.1016164 
    ## ... Procrustes: rmse 3.601686e-05  max resid 0.0001024727 
    ## ... Similar to previous best
    ## Run 241 stress 0.1341868 
    ## Run 242 stress 0.1016164 
    ## ... Procrustes: rmse 2.090119e-05  max resid 5.986533e-05 
    ## ... Similar to previous best
    ## Run 243 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001199154  max resid 0.0003390578 
    ## ... Similar to previous best
    ## Run 244 stress 0.1356148 
    ## Run 245 stress 0.1356148 
    ## Run 246 stress 0.1016164 
    ## ... Procrustes: rmse 2.241556e-05  max resid 6.417256e-05 
    ## ... Similar to previous best
    ## Run 247 stress 0.1016166 
    ## ... Procrustes: rmse 0.000351285  max resid 0.0009834009 
    ## ... Similar to previous best
    ## Run 248 stress 0.1016164 
    ## ... Procrustes: rmse 6.590875e-05  max resid 0.0001869751 
    ## ... Similar to previous best
    ## Run 249 stress 0.1016164 
    ## ... Procrustes: rmse 9.718489e-06  max resid 2.819265e-05 
    ## ... Similar to previous best
    ## Run 250 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001347552  max resid 0.0003770922 
    ## ... Similar to previous best
    ## Run 251 stress 0.1016164 
    ## ... Procrustes: rmse 1.681033e-05  max resid 4.809394e-05 
    ## ... Similar to previous best
    ## Run 252 stress 0.1016164 
    ## ... Procrustes: rmse 4.50502e-05  max resid 0.0001279749 
    ## ... Similar to previous best
    ## Run 253 stress 0.1356147 
    ## Run 254 stress 0.1016164 
    ## ... Procrustes: rmse 2.690372e-05  max resid 7.678541e-05 
    ## ... Similar to previous best
    ## Run 255 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004768032  max resid 0.001331842 
    ## ... Similar to previous best
    ## Run 256 stress 0.1016164 
    ## ... Procrustes: rmse 4.416317e-05  max resid 0.0001255813 
    ## ... Similar to previous best
    ## Run 257 stress 0.1016164 
    ## ... Procrustes: rmse 2.422639e-05  max resid 6.490019e-05 
    ## ... Similar to previous best
    ## Run 258 stress 0.1016164 
    ## ... Procrustes: rmse 1.61769e-05  max resid 4.628948e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003809551  max resid 0.001064198 
    ## ... Similar to previous best
    ## Run 260 stress 0.1016164 
    ## ... Procrustes: rmse 1.565703e-05  max resid 3.938664e-05 
    ## ... Similar to previous best
    ## Run 261 stress 0.1519033 
    ## Run 262 stress 0.1016164 
    ## ... Procrustes: rmse 3.84679e-05  max resid 0.0001052919 
    ## ... Similar to previous best
    ## Run 263 stress 0.1341868 
    ## Run 264 stress 0.1341868 
    ## Run 265 stress 0.1519033 
    ## Run 266 stress 0.1341868 
    ## Run 267 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001133722  max resid 0.0003176836 
    ## ... Similar to previous best
    ## Run 268 stress 0.2109396 
    ## Run 269 stress 0.101617 
    ## ... Procrustes: rmse 0.0006168312  max resid 0.001721701 
    ## ... Similar to previous best
    ## Run 270 stress 0.1537555 
    ## Run 271 stress 0.1530428 
    ## Run 272 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003352127  max resid 0.0009351551 
    ## ... Similar to previous best
    ## Run 273 stress 0.1016164 
    ## ... Procrustes: rmse 3.238958e-06  max resid 9.607538e-06 
    ## ... Similar to previous best
    ## Run 274 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005771161  max resid 0.001611661 
    ## ... Similar to previous best
    ## Run 275 stress 0.1356144 
    ## Run 276 stress 0.1016165 
    ## ... Procrustes: rmse 0.000171947  max resid 0.0004817772 
    ## ... Similar to previous best
    ## Run 277 stress 0.1016164 
    ## ... Procrustes: rmse 2.325979e-05  max resid 6.645008e-05 
    ## ... Similar to previous best
    ## Run 278 stress 0.1341868 
    ## Run 279 stress 0.1341868 
    ## Run 280 stress 0.135615 
    ## Run 281 stress 0.1341868 
    ## Run 282 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004880959  max resid 0.001363107 
    ## ... Similar to previous best
    ## Run 283 stress 0.1356146 
    ## Run 284 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001069425  max resid 0.0003029358 
    ## ... Similar to previous best
    ## Run 285 stress 0.1016164 
    ## ... Procrustes: rmse 2.370194e-05  max resid 6.786423e-05 
    ## ... Similar to previous best
    ## Run 286 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003604346  max resid 0.001008107 
    ## ... Similar to previous best
    ## Run 287 stress 0.1016164 
    ## ... Procrustes: rmse 1.898272e-05  max resid 5.157205e-05 
    ## ... Similar to previous best
    ## Run 288 stress 0.1341868 
    ## Run 289 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003483516  max resid 0.0009734656 
    ## ... Similar to previous best
    ## Run 290 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005501238  max resid 0.001536267 
    ## ... Similar to previous best
    ## Run 291 stress 0.1519033 
    ## Run 292 stress 0.1016164 
    ## ... Procrustes: rmse 4.336158e-05  max resid 0.0001215697 
    ## ... Similar to previous best
    ## Run 293 stress 0.1530428 
    ## Run 294 stress 0.1341868 
    ## Run 295 stress 0.1016164 
    ## ... Procrustes: rmse 7.595978e-05  max resid 0.0002145177 
    ## ... Similar to previous best
    ## Run 296 stress 0.1016164 
    ## ... Procrustes: rmse 6.683437e-05  max resid 0.0001895991 
    ## ... Similar to previous best
    ## Run 297 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003645121  max resid 0.001018796 
    ## ... Similar to previous best
    ## Run 298 stress 0.1530425 
    ## Run 299 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004214653  max resid 0.001176801 
    ## ... Similar to previous best
    ## Run 300 stress 0.1016164 
    ## ... Procrustes: rmse 2.907318e-05  max resid 8.291299e-05 
    ## ... Similar to previous best
    ## Run 301 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001032147  max resid 0.0002889387 
    ## ... Similar to previous best
    ## Run 302 stress 0.1530428 
    ## Run 303 stress 0.101617 
    ## ... Procrustes: rmse 0.0006452058  max resid 0.001802081 
    ## ... Similar to previous best
    ## Run 304 stress 0.1016164 
    ## ... Procrustes: rmse 6.108415e-05  max resid 0.0001718773 
    ## ... Similar to previous best
    ## Run 305 stress 0.3394536 
    ## Run 306 stress 0.1356151 
    ## Run 307 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005663326  max resid 0.00157894 
    ## ... Similar to previous best
    ## Run 308 stress 0.1016164 
    ## ... Procrustes: rmse 3.183351e-05  max resid 8.075927e-05 
    ## ... Similar to previous best
    ## Run 309 stress 0.1016164 
    ## ... Procrustes: rmse 1.531221e-05  max resid 4.184042e-05 
    ## ... Similar to previous best
    ## Run 310 stress 0.135615 
    ## Run 311 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001352444  max resid 0.0003786063 
    ## ... Similar to previous best
    ## Run 312 stress 0.1016164 
    ## ... Procrustes: rmse 1.260347e-05  max resid 2.240691e-05 
    ## ... Similar to previous best
    ## Run 313 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002764801  max resid 0.0007739284 
    ## ... Similar to previous best
    ## Run 314 stress 0.1356148 
    ## Run 315 stress 0.1356143 
    ## Run 316 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004741169  max resid 0.001323852 
    ## ... Similar to previous best
    ## Run 317 stress 0.135615 
    ## Run 318 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004821759  max resid 0.001347252 
    ## ... Similar to previous best
    ## Run 319 stress 0.153043 
    ## Run 320 stress 0.1356148 
    ## Run 321 stress 0.1341868 
    ## Run 322 stress 0.1016164 
    ## ... Procrustes: rmse 9.995623e-06  max resid 2.891816e-05 
    ## ... Similar to previous best
    ## Run 323 stress 0.1016164 
    ## ... Procrustes: rmse 5.711974e-06  max resid 1.251256e-05 
    ## ... Similar to previous best
    ## Run 324 stress 0.1530427 
    ## Run 325 stress 0.1016164 
    ## ... Procrustes: rmse 9.896112e-05  max resid 0.000280194 
    ## ... Similar to previous best
    ## Run 326 stress 0.1341868 
    ## Run 327 stress 0.1341868 
    ## Run 328 stress 0.1016164 
    ## ... Procrustes: rmse 7.374458e-05  max resid 0.0002091674 
    ## ... Similar to previous best
    ## Run 329 stress 0.1016164 
    ## ... Procrustes: rmse 4.196742e-06  max resid 1.200351e-05 
    ## ... Similar to previous best
    ## Run 330 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003246123  max resid 0.0009078864 
    ## ... Similar to previous best
    ## Run 331 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001016443  max resid 0.0002822329 
    ## ... Similar to previous best
    ## Run 332 stress 0.1016164 
    ## ... Procrustes: rmse 6.045384e-05  max resid 0.0001712502 
    ## ... Similar to previous best
    ## Run 333 stress 0.1016164 
    ## ... Procrustes: rmse 3.84914e-06  max resid 1.111612e-05 
    ## ... Similar to previous best
    ## Run 334 stress 0.1016165 
    ## ... Procrustes: rmse 4.74122e-05  max resid 0.000123353 
    ## ... Similar to previous best
    ## Run 335 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003262639  max resid 0.0009124178 
    ## ... Similar to previous best
    ## Run 336 stress 0.1016164 
    ## ... Procrustes: rmse 2.659355e-05  max resid 7.445329e-05 
    ## ... Similar to previous best
    ## Run 337 stress 0.1016164 
    ## ... Procrustes: rmse 4.68562e-05  max resid 0.0001294013 
    ## ... Similar to previous best
    ## Run 338 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004067377  max resid 0.001136584 
    ## ... Similar to previous best
    ## Run 339 stress 0.1016164 
    ## ... Procrustes: rmse 1.996745e-05  max resid 5.724814e-05 
    ## ... Similar to previous best
    ## Run 340 stress 0.1016164 
    ## ... Procrustes: rmse 8.705209e-06  max resid 2.522401e-05 
    ## ... Similar to previous best
    ## Run 341 stress 0.135615 
    ## Run 342 stress 0.1016164 
    ## ... Procrustes: rmse 4.672565e-05  max resid 0.0001309274 
    ## ... Similar to previous best
    ## Run 343 stress 0.1537554 
    ## Run 344 stress 0.1016164 
    ## ... Procrustes: rmse 2.353928e-05  max resid 6.720683e-05 
    ## ... Similar to previous best
    ## Run 345 stress 0.1519033 
    ## Run 346 stress 0.1016164 
    ## ... Procrustes: rmse 9.861045e-05  max resid 0.0002664087 
    ## ... Similar to previous best
    ## Run 347 stress 0.2074952 
    ## Run 348 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004578148  max resid 0.001280105 
    ## ... Similar to previous best
    ## Run 349 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001053098  max resid 0.0002975235 
    ## ... Similar to previous best
    ## Run 350 stress 0.1016164 
    ## ... Procrustes: rmse 5.708063e-05  max resid 0.0001598998 
    ## ... Similar to previous best
    ## Run 351 stress 0.1530425 
    ## Run 352 stress 0.1016164 
    ## ... Procrustes: rmse 2.056139e-05  max resid 4.942833e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.1530425 
    ## Run 354 stress 0.1341868 
    ## Run 355 stress 0.1016164 
    ## ... Procrustes: rmse 9.794813e-05  max resid 0.0002740222 
    ## ... Similar to previous best
    ## Run 356 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002682673  max resid 0.0007513313 
    ## ... Similar to previous best
    ## Run 357 stress 0.1016164 
    ## ... Procrustes: rmse 1.022013e-05  max resid 2.857502e-05 
    ## ... Similar to previous best
    ## Run 358 stress 0.1016169 
    ## ... Procrustes: rmse 0.000576732  max resid 0.001610429 
    ## ... Similar to previous best
    ## Run 359 stress 0.1016164 
    ## ... Procrustes: rmse 1.772569e-05  max resid 5.071081e-05 
    ## ... Similar to previous best
    ## Run 360 stress 0.1341868 
    ## Run 361 stress 0.101617 
    ## ... Procrustes: rmse 0.0006518141  max resid 0.001820855 
    ## ... Similar to previous best
    ## Run 362 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004733349  max resid 0.001322614 
    ## ... Similar to previous best
    ## Run 363 stress 0.1356148 
    ## Run 364 stress 0.1341868 
    ## Run 365 stress 0.1341868 
    ## Run 366 stress 0.1016164 
    ## ... Procrustes: rmse 3.241399e-05  max resid 9.248819e-05 
    ## ... Similar to previous best
    ## Run 367 stress 0.3198734 
    ## Run 368 stress 0.1016164 
    ## ... Procrustes: rmse 2.717784e-05  max resid 7.76494e-05 
    ## ... Similar to previous best
    ## Run 369 stress 0.1016164 
    ## ... Procrustes: rmse 2.947606e-05  max resid 7.823172e-05 
    ## ... Similar to previous best
    ## Run 370 stress 0.1016164 
    ## ... Procrustes: rmse 5.023117e-05  max resid 0.000142282 
    ## ... Similar to previous best
    ## Run 371 stress 0.1016164 
    ## ... Procrustes: rmse 1.346168e-05  max resid 3.887492e-05 
    ## ... Similar to previous best
    ## Run 372 stress 0.1341868 
    ## Run 373 stress 0.1341868 
    ## Run 374 stress 0.1341868 
    ## Run 375 stress 0.1016164 
    ## ... Procrustes: rmse 3.124459e-05  max resid 8.847156e-05 
    ## ... Similar to previous best
    ## Run 376 stress 0.1341868 
    ## Run 377 stress 0.1341868 
    ## Run 378 stress 0.1341868 
    ## Run 379 stress 0.1341868 
    ## Run 380 stress 0.1341868 
    ## Run 381 stress 0.1530423 
    ## Run 382 stress 0.135615 
    ## Run 383 stress 0.1356151 
    ## Run 384 stress 0.1356148 
    ## Run 385 stress 0.1016164 
    ## ... Procrustes: rmse 4.99592e-05  max resid 0.0001419826 
    ## ... Similar to previous best
    ## Run 386 stress 0.101617 
    ## ... Procrustes: rmse 0.0006845037  max resid 0.001909954 
    ## ... Similar to previous best
    ## Run 387 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003329849  max resid 0.0009323353 
    ## ... Similar to previous best
    ## Run 388 stress 0.1341868 
    ## Run 389 stress 0.1341868 
    ## Run 390 stress 0.1016164 
    ## ... Procrustes: rmse 4.076172e-05  max resid 0.0001153914 
    ## ... Similar to previous best
    ## Run 391 stress 0.1341868 
    ## Run 392 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002835299  max resid 0.0007929468 
    ## ... Similar to previous best
    ## Run 393 stress 0.1341868 
    ## Run 394 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002981349  max resid 0.0008330312 
    ## ... Similar to previous best
    ## Run 395 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002737746  max resid 0.0007657987 
    ## ... Similar to previous best
    ## Run 396 stress 0.1016167 
    ## ... Procrustes: rmse 0.000490014  max resid 0.001365585 
    ## ... Similar to previous best
    ## Run 397 stress 0.1341868 
    ## Run 398 stress 0.1341868 
    ## Run 399 stress 0.1519033 
    ## Run 400 stress 0.101617 
    ## ... Procrustes: rmse 0.0006447399  max resid 0.001795574 
    ## ... Similar to previous best
    ## Run 401 stress 0.1341868 
    ## Run 402 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001197832  max resid 0.0003388361 
    ## ... Similar to previous best
    ## Run 403 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003160406  max resid 0.0008849228 
    ## ... Similar to previous best
    ## Run 404 stress 0.1016164 
    ## ... Procrustes: rmse 4.665614e-06  max resid 1.360977e-05 
    ## ... Similar to previous best
    ## Run 405 stress 0.1016164 
    ## ... Procrustes: rmse 7.48043e-05  max resid 0.0002117712 
    ## ... Similar to previous best
    ## Run 406 stress 0.1530426 
    ## Run 407 stress 0.1016164 
    ## ... Procrustes: rmse 1.933151e-06  max resid 5.484412e-06 
    ## ... Similar to previous best
    ## Run 408 stress 0.1016164 
    ## ... Procrustes: rmse 1.297327e-05  max resid 3.742849e-05 
    ## ... Similar to previous best
    ## Run 409 stress 0.1356149 
    ## Run 410 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003771788  max resid 0.001053072 
    ## ... Similar to previous best
    ## Run 411 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 2.084002e-06  max resid 5.254721e-06 
    ## ... Similar to previous best
    ## Run 412 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006331259  max resid 0.001768769 
    ## ... Similar to previous best
    ## Run 413 stress 0.1016164 
    ## ... Procrustes: rmse 8.410831e-06  max resid 2.064519e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.1016164 
    ## ... Procrustes: rmse 0.000100139  max resid 0.0002787508 
    ## ... Similar to previous best
    ## Run 415 stress 0.1341868 
    ## Run 416 stress 0.1016164 
    ## ... Procrustes: rmse 4.749616e-05  max resid 0.0001242739 
    ## ... Similar to previous best
    ## Run 417 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003168331  max resid 0.0008867555 
    ## ... Similar to previous best
    ## Run 418 stress 0.1016164 
    ## ... Procrustes: rmse 5.882631e-06  max resid 1.560583e-05 
    ## ... Similar to previous best
    ## Run 419 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005994496  max resid 0.001674447 
    ## ... Similar to previous best
    ## Run 420 stress 0.1016164 
    ## ... Procrustes: rmse 6.687596e-05  max resid 0.0001880349 
    ## ... Similar to previous best
    ## Run 421 stress 0.101617 
    ## ... Procrustes: rmse 0.0006129874  max resid 0.001707865 
    ## ... Similar to previous best
    ## Run 422 stress 0.1016164 
    ## ... Procrustes: rmse 1.965082e-05  max resid 5.532036e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.1016164 
    ## ... Procrustes: rmse 1.349883e-05  max resid 2.924201e-05 
    ## ... Similar to previous best
    ## Run 424 stress 0.1016164 
    ## ... Procrustes: rmse 2.786798e-05  max resid 7.865621e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006024321  max resid 0.00168179 
    ## ... Similar to previous best
    ## Run 426 stress 0.1530428 
    ## Run 427 stress 0.1016164 
    ## ... Procrustes: rmse 3.681022e-05  max resid 0.0001036321 
    ## ... Similar to previous best
    ## Run 428 stress 0.1341868 
    ## Run 429 stress 0.1530423 
    ## Run 430 stress 0.1016164 
    ## ... Procrustes: rmse 4.028481e-05  max resid 0.000105349 
    ## ... Similar to previous best
    ## Run 431 stress 0.1016164 
    ## ... Procrustes: rmse 9.445251e-05  max resid 0.0002651989 
    ## ... Similar to previous best
    ## Run 432 stress 0.1356147 
    ## Run 433 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004483268  max resid 0.001253972 
    ## ... Similar to previous best
    ## Run 434 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002664849  max resid 0.0007469363 
    ## ... Similar to previous best
    ## Run 435 stress 0.1016164 
    ## ... Procrustes: rmse 1.728439e-05  max resid 4.818909e-05 
    ## ... Similar to previous best
    ## Run 436 stress 0.1341868 
    ## Run 437 stress 0.1341868 
    ## Run 438 stress 0.1016164 
    ## ... Procrustes: rmse 6.786767e-06  max resid 1.914888e-05 
    ## ... Similar to previous best
    ## Run 439 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005337268  max resid 0.001492852 
    ## ... Similar to previous best
    ## Run 440 stress 0.1530422 
    ## Run 441 stress 0.1016168 
    ## ... Procrustes: rmse 0.00051275  max resid 0.001433854 
    ## ... Similar to previous best
    ## Run 442 stress 0.1016164 
    ## ... Procrustes: rmse 2.963199e-05  max resid 8.358133e-05 
    ## ... Similar to previous best
    ## Run 443 stress 0.1341868 
    ## Run 444 stress 0.1537555 
    ## Run 445 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003957777  max resid 0.001107239 
    ## ... Similar to previous best
    ## Run 446 stress 0.1341868 
    ## Run 447 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005666294  max resid 0.001582974 
    ## ... Similar to previous best
    ## Run 448 stress 0.1341868 
    ## Run 449 stress 0.1016164 
    ## ... Procrustes: rmse 3.895457e-05  max resid 0.0001096481 
    ## ... Similar to previous best
    ## Run 450 stress 0.1356149 
    ## Run 451 stress 0.1016165 
    ## ... Procrustes: rmse 0.000257465  max resid 0.0007209956 
    ## ... Similar to previous best
    ## Run 452 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004304105  max resid 0.001202184 
    ## ... Similar to previous best
    ## Run 453 stress 0.1341868 
    ## Run 454 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005703895  max resid 0.001593309 
    ## ... Similar to previous best
    ## Run 455 stress 0.1356148 
    ## Run 456 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006154708  max resid 0.001720199 
    ## ... Similar to previous best
    ## Run 457 stress 0.1016164 
    ## ... Procrustes: rmse 9.946827e-05  max resid 0.0002798443 
    ## ... Similar to previous best
    ## Run 458 stress 0.1341868 
    ## Run 459 stress 0.1016164 
    ## ... Procrustes: rmse 4.080657e-07  max resid 6.324753e-07 
    ## ... Similar to previous best
    ## Run 460 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007314565  max resid 0.002042626 
    ## ... Similar to previous best
    ## Run 461 stress 0.1341868 
    ## Run 462 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005226384  max resid 0.001460842 
    ## ... Similar to previous best
    ## Run 463 stress 0.1341868 
    ## Run 464 stress 0.135615 
    ## Run 465 stress 0.1341868 
    ## Run 466 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005735409  max resid 0.001601955 
    ## ... Similar to previous best
    ## Run 467 stress 0.1016164 
    ## ... Procrustes: rmse 1.968884e-06  max resid 5.198959e-06 
    ## ... Similar to previous best
    ## Run 468 stress 0.1016164 
    ## ... Procrustes: rmse 0.000152879  max resid 0.0004287161 
    ## ... Similar to previous best
    ## Run 469 stress 0.1016164 
    ## ... Procrustes: rmse 3.762148e-05  max resid 0.0001052164 
    ## ... Similar to previous best
    ## Run 470 stress 0.1016164 
    ## ... Procrustes: rmse 2.322219e-06  max resid 5.615235e-06 
    ## ... Similar to previous best
    ## Run 471 stress 0.1016169 
    ## ... Procrustes: rmse 0.0004757027  max resid 0.001324679 
    ## ... Similar to previous best
    ## Run 472 stress 0.1016165 
    ## ... Procrustes: rmse 0.000251085  max resid 0.0007044355 
    ## ... Similar to previous best
    ## Run 473 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003575399  max resid 0.0009988432 
    ## ... Similar to previous best
    ## Run 474 stress 0.1016166 
    ## ... Procrustes: rmse 0.0001519884  max resid 0.0004095987 
    ## ... Similar to previous best
    ## Run 475 stress 0.101617 
    ## ... Procrustes: rmse 0.0006643545  max resid 0.001856745 
    ## ... Similar to previous best
    ## Run 476 stress 0.1016164 
    ## ... Procrustes: rmse 2.199973e-05  max resid 3.966364e-05 
    ## ... Similar to previous best
    ## Run 477 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003796761  max resid 0.001062146 
    ## ... Similar to previous best
    ## Run 478 stress 0.1016164 
    ## ... Procrustes: rmse 1.205528e-05  max resid 3.389234e-05 
    ## ... Similar to previous best
    ## Run 479 stress 0.1016164 
    ## ... Procrustes: rmse 5.494106e-06  max resid 1.374874e-05 
    ## ... Similar to previous best
    ## Run 480 stress 0.1016164 
    ## ... Procrustes: rmse 5.274816e-05  max resid 0.0001486077 
    ## ... Similar to previous best
    ## Run 481 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003502837  max resid 0.0009799854 
    ## ... Similar to previous best
    ## Run 482 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004657922  max resid 0.001302301 
    ## ... Similar to previous best
    ## Run 483 stress 0.1016164 
    ## ... Procrustes: rmse 2.248129e-05  max resid 6.32889e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.1016164 
    ## ... Procrustes: rmse 5.380785e-05  max resid 0.0001515209 
    ## ... Similar to previous best
    ## Run 485 stress 0.1341868 
    ## Run 486 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003934305  max resid 0.001100217 
    ## ... Similar to previous best
    ## Run 487 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003703545  max resid 0.001035756 
    ## ... Similar to previous best
    ## Run 488 stress 0.1016164 
    ## ... Procrustes: rmse 5.431103e-05  max resid 0.0001528431 
    ## ... Similar to previous best
    ## Run 489 stress 0.1356147 
    ## Run 490 stress 0.1016164 
    ## ... Procrustes: rmse 5.121215e-05  max resid 0.000141232 
    ## ... Similar to previous best
    ## Run 491 stress 0.1016164 
    ## ... Procrustes: rmse 3.021386e-05  max resid 8.500253e-05 
    ## ... Similar to previous best
    ## Run 492 stress 0.1016164 
    ## ... Procrustes: rmse 1.754798e-05  max resid 4.664707e-05 
    ## ... Similar to previous best
    ## Run 493 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001578466  max resid 0.0004427598 
    ## ... Similar to previous best
    ## Run 494 stress 0.1016164 
    ## ... Procrustes: rmse 5.46699e-06  max resid 1.020689e-05 
    ## ... Similar to previous best
    ## Run 495 stress 0.1016164 
    ## ... Procrustes: rmse 3.612631e-05  max resid 0.0001016923 
    ## ... Similar to previous best
    ## Run 496 stress 0.1356148 
    ## Run 497 stress 0.1341868 
    ## Run 498 stress 0.1016164 
    ## ... Procrustes: rmse 1.59005e-05  max resid 4.461235e-05 
    ## ... Similar to previous best
    ## Run 499 stress 0.1341868 
    ## Run 500 stress 0.1016164 
    ## ... Procrustes: rmse 1.38621e-05  max resid 3.192627e-05 
    ## ... Similar to previous best
    ## *** Best solution repeated 65 times

``` r
# Stratified lakes and ocean sites
PD_beta_geo_SO_NMDS <- metaMDS(PD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.04851417 
    ## Run 1 stress 0.04851417 
    ## ... New best solution
    ## ... Procrustes: rmse 5.022817e-05  max resid 0.0001409192 
    ## ... Similar to previous best
    ## Run 2 stress 0.04935947 
    ## Run 3 stress 0.05728059 
    ## Run 4 stress 0.04935947 
    ## Run 5 stress 0.05728094 
    ## Run 6 stress 0.04938615 
    ## Run 7 stress 0.05855031 
    ## Run 8 stress 0.05643265 
    ## Run 9 stress 0.04851417 
    ## ... New best solution
    ## ... Procrustes: rmse 2.018848e-05  max resid 5.207196e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.05276412 
    ## Run 11 stress 0.05607344 
    ## Run 12 stress 0.05607344 
    ## Run 13 stress 0.05439286 
    ## Run 14 stress 0.04938617 
    ## Run 15 stress 0.06185071 
    ## Run 16 stress 0.05612239 
    ## Run 17 stress 0.05728066 
    ## Run 18 stress 0.04935953 
    ## Run 19 stress 0.05433173 
    ## Run 20 stress 0.05572688 
    ## Run 21 stress 0.04938615 
    ## Run 22 stress 0.04612762 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04985707  max resid 0.1396717 
    ## Run 23 stress 0.04935944 
    ## Run 24 stress 0.05728066 
    ## Run 25 stress 0.05605146 
    ## Run 26 stress 0.04612758 
    ## ... New best solution
    ## ... Procrustes: rmse 5.061711e-05  max resid 0.0001035205 
    ## ... Similar to previous best
    ## Run 27 stress 0.05572688 
    ## Run 28 stress 0.06205759 
    ## Run 29 stress 0.04851417 
    ## Run 30 stress 0.05216138 
    ## Run 31 stress 0.05607346 
    ## Run 32 stress 0.05728073 
    ## Run 33 stress 0.0561224 
    ## Run 34 stress 0.05605145 
    ## Run 35 stress 0.05695856 
    ## Run 36 stress 0.06205754 
    ## Run 37 stress 0.05260046 
    ## Run 38 stress 0.05276411 
    ## Run 39 stress 0.06088979 
    ## Run 40 stress 0.05643261 
    ## Run 41 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 5.985126e-05  max resid 0.000119249 
    ## ... Similar to previous best
    ## Run 42 stress 0.04612766 
    ## ... Procrustes: rmse 0.0001391456  max resid 0.0002828072 
    ## ... Similar to previous best
    ## Run 43 stress 0.05643261 
    ## Run 44 stress 0.04860422 
    ## Run 45 stress 0.06205752 
    ## Run 46 stress 0.054515 
    ## Run 47 stress 0.0526005 
    ## Run 48 stress 0.04938621 
    ## Run 49 stress 0.04938619 
    ## Run 50 stress 0.05137969 
    ## Run 51 stress 0.0526005 
    ## Run 52 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001232608  max resid 0.0002476004 
    ## ... Similar to previous best
    ## Run 53 stress 0.04938616 
    ## Run 54 stress 0.06291433 
    ## Run 55 stress 0.06088979 
    ## Run 56 stress 0.3374347 
    ## Run 57 stress 0.05643261 
    ## Run 58 stress 0.05605148 
    ## Run 59 stress 0.04851417 
    ## Run 60 stress 0.04612759 
    ## ... Procrustes: rmse 7.680468e-05  max resid 0.0001561289 
    ## ... Similar to previous best
    ## Run 61 stress 0.04935947 
    ## Run 62 stress 0.04851417 
    ## Run 63 stress 0.06205761 
    ## Run 64 stress 0.04938616 
    ## Run 65 stress 0.04851417 
    ## Run 66 stress 0.0569586 
    ## Run 67 stress 0.07046354 
    ## Run 68 stress 0.05728057 
    ## Run 69 stress 0.04860423 
    ## Run 70 stress 0.04851417 
    ## Run 71 stress 0.05612244 
    ## Run 72 stress 0.04851418 
    ## Run 73 stress 0.05451493 
    ## Run 74 stress 0.0585503 
    ## Run 75 stress 0.04612758 
    ## ... Procrustes: rmse 6.131917e-05  max resid 0.0001219805 
    ## ... Similar to previous best
    ## Run 76 stress 0.05855029 
    ## Run 77 stress 0.05855031 
    ## Run 78 stress 0.05607346 
    ## Run 79 stress 0.05260056 
    ## Run 80 stress 0.05643261 
    ## Run 81 stress 0.06291425 
    ## Run 82 stress 0.04851417 
    ## Run 83 stress 0.05855033 
    ## Run 84 stress 0.04938616 
    ## Run 85 stress 0.05842818 
    ## Run 86 stress 0.07046351 
    ## Run 87 stress 0.05643261 
    ## Run 88 stress 0.04851417 
    ## Run 89 stress 0.06088978 
    ## Run 90 stress 0.04860422 
    ## Run 91 stress 0.04612757 
    ## ... Procrustes: rmse 2.360823e-05  max resid 4.446015e-05 
    ## ... Similar to previous best
    ## Run 92 stress 0.05439286 
    ## Run 93 stress 0.05572694 
    ## Run 94 stress 0.05438939 
    ## Run 95 stress 0.05643261 
    ## Run 96 stress 0.0608897 
    ## Run 97 stress 0.05474817 
    ## Run 98 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001361123  max resid 0.0002775828 
    ## ... Similar to previous best
    ## Run 99 stress 0.04851419 
    ## Run 100 stress 0.0620576 
    ## Run 101 stress 0.05643261 
    ## Run 102 stress 0.06291425 
    ## Run 103 stress 0.04851417 
    ## Run 104 stress 0.06400747 
    ## Run 105 stress 0.06205755 
    ## Run 106 stress 0.05643261 
    ## Run 107 stress 0.05605148 
    ## Run 108 stress 0.05728063 
    ## Run 109 stress 0.05605159 
    ## Run 110 stress 0.04851417 
    ## Run 111 stress 0.05728083 
    ## Run 112 stress 0.04612757 
    ## ... Procrustes: rmse 4.203634e-05  max resid 8.329616e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.04851419 
    ## Run 114 stress 0.04851417 
    ## Run 115 stress 0.05643265 
    ## Run 116 stress 0.04860422 
    ## Run 117 stress 0.0608897 
    ## Run 118 stress 0.05593468 
    ## Run 119 stress 0.05276411 
    ## Run 120 stress 0.04851418 
    ## Run 121 stress 0.06205756 
    ## Run 122 stress 0.04938616 
    ## Run 123 stress 0.05438943 
    ## Run 124 stress 0.05695879 
    ## Run 125 stress 0.05643261 
    ## Run 126 stress 0.04938617 
    ## Run 127 stress 0.05451493 
    ## Run 128 stress 0.04938615 
    ## Run 129 stress 0.05605152 
    ## Run 130 stress 0.05260051 
    ## Run 131 stress 0.04938615 
    ## Run 132 stress 0.04860422 
    ## Run 133 stress 0.04851418 
    ## Run 134 stress 0.04612758 
    ## ... Procrustes: rmse 1.736205e-05  max resid 3.164367e-05 
    ## ... Similar to previous best
    ## Run 135 stress 0.05572689 
    ## Run 136 stress 0.3411951 
    ## Run 137 stress 0.05728056 
    ## Run 138 stress 0.05605146 
    ## Run 139 stress 0.05728068 
    ## Run 140 stress 0.05572688 
    ## Run 141 stress 0.05607345 
    ## Run 142 stress 0.04612767 
    ## ... Procrustes: rmse 0.0001545151  max resid 0.0003116211 
    ## ... Similar to previous best
    ## Run 143 stress 0.05488107 
    ## Run 144 stress 0.05695872 
    ## Run 145 stress 0.05593452 
    ## Run 146 stress 0.05842781 
    ## Run 147 stress 0.04860423 
    ## Run 148 stress 0.06088973 
    ## Run 149 stress 0.05607348 
    ## Run 150 stress 0.05605147 
    ## Run 151 stress 0.05695857 
    ## Run 152 stress 0.04860423 
    ## Run 153 stress 0.04851417 
    ## Run 154 stress 0.05728071 
    ## Run 155 stress 0.05605148 
    ## Run 156 stress 0.04851417 
    ## Run 157 stress 0.05439281 
    ## Run 158 stress 0.05216134 
    ## Run 159 stress 0.0521615 
    ## Run 160 stress 0.05605158 
    ## Run 161 stress 0.04851417 
    ## Run 162 stress 0.06185055 
    ## Run 163 stress 0.0527641 
    ## Run 164 stress 0.04938617 
    ## Run 165 stress 0.05605152 
    ## Run 166 stress 0.05216135 
    ## Run 167 stress 0.0560735 
    ## Run 168 stress 0.06088977 
    ## Run 169 stress 0.05643261 
    ## Run 170 stress 0.05612239 
    ## Run 171 stress 0.04938616 
    ## Run 172 stress 0.05474814 
    ## Run 173 stress 0.0640074 
    ## Run 174 stress 0.06291425 
    ## Run 175 stress 0.0585503 
    ## Run 176 stress 0.05439285 
    ## Run 177 stress 0.05607344 
    ## Run 178 stress 0.05137964 
    ## Run 179 stress 0.04851417 
    ## Run 180 stress 0.0543895 
    ## Run 181 stress 0.04935946 
    ## Run 182 stress 0.0620576 
    ## Run 183 stress 0.05607344 
    ## Run 184 stress 0.05643268 
    ## Run 185 stress 0.05695859 
    ## Run 186 stress 0.04935944 
    ## Run 187 stress 0.05276412 
    ## Run 188 stress 0.04612758 
    ## ... Procrustes: rmse 5.764233e-05  max resid 0.0001138022 
    ## ... Similar to previous best
    ## Run 189 stress 0.05643266 
    ## Run 190 stress 0.05433185 
    ## Run 191 stress 0.04851417 
    ## Run 192 stress 0.05643261 
    ## Run 193 stress 0.05695865 
    ## Run 194 stress 0.06291425 
    ## Run 195 stress 0.05572689 
    ## Run 196 stress 0.05605147 
    ## Run 197 stress 0.06088971 
    ## Run 198 stress 0.04851418 
    ## Run 199 stress 0.04612759 
    ## ... Procrustes: rmse 6.517971e-05  max resid 0.0001323319 
    ## ... Similar to previous best
    ## Run 200 stress 0.04851418 
    ## Run 201 stress 0.05433186 
    ## Run 202 stress 0.05728067 
    ## Run 203 stress 0.05488108 
    ## Run 204 stress 0.05728082 
    ## Run 205 stress 0.06185055 
    ## Run 206 stress 0.05488105 
    ## Run 207 stress 0.04851417 
    ## Run 208 stress 0.06205752 
    ## Run 209 stress 0.04851417 
    ## Run 210 stress 0.04935947 
    ## Run 211 stress 0.05216138 
    ## Run 212 stress 0.05216146 
    ## Run 213 stress 0.04938616 
    ## Run 214 stress 0.05605148 
    ## Run 215 stress 0.04612761 
    ## ... Procrustes: rmse 0.000102917  max resid 0.0002075051 
    ## ... Similar to previous best
    ## Run 216 stress 0.05607345 
    ## Run 217 stress 0.05474823 
    ## Run 218 stress 0.05612243 
    ## Run 219 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001396611  max resid 0.0002858442 
    ## ... Similar to previous best
    ## Run 220 stress 0.05474817 
    ## Run 221 stress 0.04938616 
    ## Run 222 stress 0.04860422 
    ## Run 223 stress 0.04935946 
    ## Run 224 stress 0.05959902 
    ## Run 225 stress 0.04851417 
    ## Run 226 stress 0.06205753 
    ## Run 227 stress 0.04612765 
    ## ... Procrustes: rmse 0.0001338283  max resid 0.0002702356 
    ## ... Similar to previous best
    ## Run 228 stress 0.0585503 
    ## Run 229 stress 0.05260057 
    ## Run 230 stress 0.05728075 
    ## Run 231 stress 0.05605151 
    ## Run 232 stress 0.04860422 
    ## Run 233 stress 0.04938615 
    ## Run 234 stress 0.05728057 
    ## Run 235 stress 0.04612758 
    ## ... Procrustes: rmse 4.452939e-05  max resid 8.850995e-05 
    ## ... Similar to previous best
    ## Run 236 stress 0.05451493 
    ## Run 237 stress 0.05451496 
    ## Run 238 stress 0.04938615 
    ## Run 239 stress 0.05612241 
    ## Run 240 stress 0.04612759 
    ## ... Procrustes: rmse 7.738576e-05  max resid 0.0001586628 
    ## ... Similar to previous best
    ## Run 241 stress 0.3495706 
    ## Run 242 stress 0.06185058 
    ## Run 243 stress 0.04938617 
    ## Run 244 stress 0.05474814 
    ## Run 245 stress 0.05643262 
    ## Run 246 stress 0.04851417 
    ## Run 247 stress 0.04851417 
    ## Run 248 stress 0.04851417 
    ## Run 249 stress 0.06205754 
    ## Run 250 stress 0.05728068 
    ## Run 251 stress 0.06185052 
    ## Run 252 stress 0.05276415 
    ## Run 253 stress 0.04851417 
    ## Run 254 stress 0.05728057 
    ## Run 255 stress 0.05728064 
    ## Run 256 stress 0.04851417 
    ## Run 257 stress 0.05474818 
    ## Run 258 stress 0.0561224 
    ## Run 259 stress 0.05572688 
    ## Run 260 stress 0.04851419 
    ## Run 261 stress 0.05216133 
    ## Run 262 stress 0.05605152 
    ## Run 263 stress 0.05438951 
    ## Run 264 stress 0.06185054 
    ## Run 265 stress 0.05607344 
    ## Run 266 stress 0.05728067 
    ## Run 267 stress 0.0461276 
    ## ... Procrustes: rmse 8.62046e-05  max resid 0.00017309 
    ## ... Similar to previous best
    ## Run 268 stress 0.04612761 
    ## ... Procrustes: rmse 9.43784e-05  max resid 0.0001891488 
    ## ... Similar to previous best
    ## Run 269 stress 0.04612758 
    ## ... Procrustes: rmse 4.809633e-05  max resid 9.526423e-05 
    ## ... Similar to previous best
    ## Run 270 stress 0.05643261 
    ## Run 271 stress 0.04938616 
    ## Run 272 stress 0.04851417 
    ## Run 273 stress 0.05137963 
    ## Run 274 stress 0.05612243 
    ## Run 275 stress 0.05607346 
    ## Run 276 stress 0.05855038 
    ## Run 277 stress 0.05216129 
    ## Run 278 stress 0.06400735 
    ## Run 279 stress 0.05728062 
    ## Run 280 stress 0.05728078 
    ## Run 281 stress 0.05137972 
    ## Run 282 stress 0.05605146 
    ## Run 283 stress 0.05643261 
    ## Run 284 stress 0.05137965 
    ## Run 285 stress 0.05855029 
    ## Run 286 stress 0.05607344 
    ## Run 287 stress 0.06088974 
    ## Run 288 stress 0.05842794 
    ## Run 289 stress 0.05612241 
    ## Run 290 stress 0.05728058 
    ## Run 291 stress 0.05855029 
    ## Run 292 stress 0.06291425 
    ## Run 293 stress 0.05842803 
    ## Run 294 stress 0.05451493 
    ## Run 295 stress 0.05488106 
    ## Run 296 stress 0.05959902 
    ## Run 297 stress 0.05643262 
    ## Run 298 stress 0.05728056 
    ## Run 299 stress 0.05695854 
    ## Run 300 stress 0.0527641 
    ## Run 301 stress 0.05438947 
    ## Run 302 stress 0.04851419 
    ## Run 303 stress 0.05607344 
    ## Run 304 stress 0.05728055 
    ## Run 305 stress 0.05728063 
    ## Run 306 stress 0.04935965 
    ## Run 307 stress 0.06185062 
    ## Run 308 stress 0.04938616 
    ## Run 309 stress 0.05607344 
    ## Run 310 stress 0.05605151 
    ## Run 311 stress 0.05605148 
    ## Run 312 stress 0.0493862 
    ## Run 313 stress 0.06185042 
    ## Run 314 stress 0.0585503 
    ## Run 315 stress 0.05607348 
    ## Run 316 stress 0.04612757 
    ## ... Procrustes: rmse 2.65422e-05  max resid 5.12567e-05 
    ## ... Similar to previous best
    ## Run 317 stress 0.05959897 
    ## Run 318 stress 0.05605154 
    ## Run 319 stress 0.04938615 
    ## Run 320 stress 0.05695874 
    ## Run 321 stress 0.05276413 
    ## Run 322 stress 0.0526005 
    ## Run 323 stress 0.05728058 
    ## Run 324 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 1.07579e-05  max resid 1.863384e-05 
    ## ... Similar to previous best
    ## Run 325 stress 0.0572806 
    ## Run 326 stress 0.06291425 
    ## Run 327 stress 0.0543928 
    ## Run 328 stress 0.05605147 
    ## Run 329 stress 0.04860422 
    ## Run 330 stress 0.05643261 
    ## Run 331 stress 0.05612244 
    ## Run 332 stress 0.05842818 
    ## Run 333 stress 0.04851418 
    ## Run 334 stress 0.0527641 
    ## Run 335 stress 0.05451494 
    ## Run 336 stress 0.05216134 
    ## Run 337 stress 0.05607345 
    ## Run 338 stress 0.0527641 
    ## Run 339 stress 0.04938615 
    ## Run 340 stress 0.05607344 
    ## Run 341 stress 0.05612248 
    ## Run 342 stress 0.06185045 
    ## Run 343 stress 0.04851419 
    ## Run 344 stress 0.05728066 
    ## Run 345 stress 0.04860423 
    ## Run 346 stress 0.05728061 
    ## Run 347 stress 0.04935957 
    ## Run 348 stress 0.05842793 
    ## Run 349 stress 0.05474816 
    ## Run 350 stress 0.06088977 
    ## Run 351 stress 0.05607345 
    ## Run 352 stress 0.06185045 
    ## Run 353 stress 0.05260053 
    ## Run 354 stress 0.05474817 
    ## Run 355 stress 0.05605148 
    ## Run 356 stress 0.05695875 
    ## Run 357 stress 0.04851417 
    ## Run 358 stress 0.05612241 
    ## Run 359 stress 0.05855033 
    ## Run 360 stress 0.0485142 
    ## Run 361 stress 0.05451493 
    ## Run 362 stress 0.05728063 
    ## Run 363 stress 0.04938618 
    ## Run 364 stress 0.04612757 
    ## ... Procrustes: rmse 1.985145e-05  max resid 4.071526e-05 
    ## ... Similar to previous best
    ## Run 365 stress 0.04938615 
    ## Run 366 stress 0.04851418 
    ## Run 367 stress 0.04938616 
    ## Run 368 stress 0.05612247 
    ## Run 369 stress 0.04851418 
    ## Run 370 stress 0.05572688 
    ## Run 371 stress 0.04935957 
    ## Run 372 stress 0.05731512 
    ## Run 373 stress 0.04938616 
    ## Run 374 stress 0.05438936 
    ## Run 375 stress 0.05474816 
    ## Run 376 stress 0.04938615 
    ## Run 377 stress 0.04851417 
    ## Run 378 stress 0.05607344 
    ## Run 379 stress 0.05260053 
    ## Run 380 stress 0.05216139 
    ## Run 381 stress 0.05607345 
    ## Run 382 stress 0.04938618 
    ## Run 383 stress 0.05607345 
    ## Run 384 stress 0.06291426 
    ## Run 385 stress 0.05433172 
    ## Run 386 stress 0.05728078 
    ## Run 387 stress 0.05260049 
    ## Run 388 stress 0.05612239 
    ## Run 389 stress 0.04938621 
    ## Run 390 stress 0.05474819 
    ## Run 391 stress 0.04612758 
    ## ... Procrustes: rmse 3.522189e-05  max resid 7.544554e-05 
    ## ... Similar to previous best
    ## Run 392 stress 0.0560516 
    ## Run 393 stress 0.04851418 
    ## Run 394 stress 0.05605149 
    ## Run 395 stress 0.06291434 
    ## Run 396 stress 0.05216134 
    ## Run 397 stress 0.05451495 
    ## Run 398 stress 0.05474814 
    ## Run 399 stress 0.33777 
    ## Run 400 stress 0.04938615 
    ## Run 401 stress 0.05451493 
    ## Run 402 stress 0.04860422 
    ## Run 403 stress 0.05695877 
    ## Run 404 stress 0.04860422 
    ## Run 405 stress 0.04938615 
    ## Run 406 stress 0.04851417 
    ## Run 407 stress 0.05643262 
    ## Run 408 stress 0.05216138 
    ## Run 409 stress 0.05643262 
    ## Run 410 stress 0.04851419 
    ## Run 411 stress 0.0572806 
    ## Run 412 stress 0.05276412 
    ## Run 413 stress 0.06400737 
    ## Run 414 stress 0.05474814 
    ## Run 415 stress 0.04860422 
    ## Run 416 stress 0.04851417 
    ## Run 417 stress 0.06205755 
    ## Run 418 stress 0.3578484 
    ## Run 419 stress 0.05572688 
    ## Run 420 stress 0.05728061 
    ## Run 421 stress 0.05572692 
    ## Run 422 stress 0.05855029 
    ## Run 423 stress 0.05433175 
    ## Run 424 stress 0.05607348 
    ## Run 425 stress 0.0493595 
    ## Run 426 stress 0.05474821 
    ## Run 427 stress 0.04860423 
    ## Run 428 stress 0.04938615 
    ## Run 429 stress 0.04851418 
    ## Run 430 stress 0.04938616 
    ## Run 431 stress 0.05593456 
    ## Run 432 stress 0.0560515 
    ## Run 433 stress 0.06291435 
    ## Run 434 stress 0.04938615 
    ## Run 435 stress 0.05612255 
    ## Run 436 stress 0.04851417 
    ## Run 437 stress 0.0585503 
    ## Run 438 stress 0.05260056 
    ## Run 439 stress 0.0572807 
    ## Run 440 stress 0.05855029 
    ## Run 441 stress 0.04851417 
    ## Run 442 stress 0.05276412 
    ## Run 443 stress 0.04612757 
    ## ... Procrustes: rmse 1.689025e-05  max resid 3.362976e-05 
    ## ... Similar to previous best
    ## Run 444 stress 0.05728053 
    ## Run 445 stress 0.06185044 
    ## Run 446 stress 0.04851418 
    ## Run 447 stress 0.05643261 
    ## Run 448 stress 0.05695867 
    ## Run 449 stress 0.04935944 
    ## Run 450 stress 0.04851418 
    ## Run 451 stress 0.06185061 
    ## Run 452 stress 0.05695867 
    ## Run 453 stress 0.05855035 
    ## Run 454 stress 0.05593452 
    ## Run 455 stress 0.05137968 
    ## Run 456 stress 0.05607345 
    ## Run 457 stress 0.05474818 
    ## Run 458 stress 0.05612247 
    ## Run 459 stress 0.04851417 
    ## Run 460 stress 0.05695869 
    ## Run 461 stress 0.05842813 
    ## Run 462 stress 0.05474814 
    ## Run 463 stress 0.05607344 
    ## Run 464 stress 0.04612757 
    ## ... Procrustes: rmse 4.943024e-05  max resid 9.924952e-05 
    ## ... Similar to previous best
    ## Run 465 stress 0.05260048 
    ## Run 466 stress 0.04935943 
    ## Run 467 stress 0.05451493 
    ## Run 468 stress 0.05572689 
    ## Run 469 stress 0.05474815 
    ## Run 470 stress 0.05728077 
    ## Run 471 stress 0.05439281 
    ## Run 472 stress 0.05731509 
    ## Run 473 stress 0.0493595 
    ## Run 474 stress 0.05612244 
    ## Run 475 stress 0.06400744 
    ## Run 476 stress 0.05842782 
    ## Run 477 stress 0.05216136 
    ## Run 478 stress 0.05842789 
    ## Run 479 stress 0.04938618 
    ## Run 480 stress 0.04851419 
    ## Run 481 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001250505  max resid 0.000251717 
    ## ... Similar to previous best
    ## Run 482 stress 0.05643261 
    ## Run 483 stress 0.05643261 
    ## Run 484 stress 0.04851417 
    ## Run 485 stress 0.05607344 
    ## Run 486 stress 0.05605151 
    ## Run 487 stress 0.05643264 
    ## Run 488 stress 0.0527641 
    ## Run 489 stress 0.05607345 
    ## Run 490 stress 0.04612757 
    ## ... Procrustes: rmse 3.313365e-05  max resid 6.680457e-05 
    ## ... Similar to previous best
    ## Run 491 stress 0.05474813 
    ## Run 492 stress 0.04938616 
    ## Run 493 stress 0.04938617 
    ## Run 494 stress 0.04851417 
    ## Run 495 stress 0.04612762 
    ## ... Procrustes: rmse 8.923498e-05  max resid 0.0001839294 
    ## ... Similar to previous best
    ## Run 496 stress 0.04612758 
    ## ... Procrustes: rmse 6.263514e-05  max resid 0.0001274459 
    ## ... Similar to previous best
    ## Run 497 stress 0.06205752 
    ## Run 498 stress 0.04935946 
    ## Run 499 stress 0.05612243 
    ## Run 500 stress 0.05695859 
    ## *** Best solution repeated 9 times

``` r
# Mixed lakes
PD_beta_geo_M_NMDS <- metaMDS(PD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.019724e-05 
    ## Run 1 stress 0.3080294 
    ## Run 2 stress 9.877091e-05 
    ## ... Procrustes: rmse 0.0001716695  max resid 0.0002985787 
    ## ... Similar to previous best
    ## Run 3 stress 8.984402e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002233403  max resid 0.0003515279 
    ## ... Similar to previous best
    ## Run 4 stress 9.202079e-05 
    ## ... Procrustes: rmse 8.331384e-05  max resid 0.000154281 
    ## ... Similar to previous best
    ## Run 5 stress 9.580871e-05 
    ## ... Procrustes: rmse 0.0002346592  max resid 0.0003542922 
    ## ... Similar to previous best
    ## Run 6 stress 9.920981e-05 
    ## ... Procrustes: rmse 0.0002076836  max resid 0.0003651905 
    ## ... Similar to previous best
    ## Run 7 stress 9.192201e-05 
    ## ... Procrustes: rmse 0.000167853  max resid 0.0002953308 
    ## ... Similar to previous best
    ## Run 8 stress 9.305833e-05 
    ## ... Procrustes: rmse 0.0001963786  max resid 0.0003472111 
    ## ... Similar to previous best
    ## Run 9 stress 9.460889e-05 
    ## ... Procrustes: rmse 0.0002300511  max resid 0.0003526777 
    ## ... Similar to previous best
    ## Run 10 stress 9.354591e-05 
    ## ... Procrustes: rmse 0.0002330372  max resid 0.0004411834 
    ## ... Similar to previous best
    ## Run 11 stress 0.2842783 
    ## Run 12 stress 9.125287e-05 
    ## ... Procrustes: rmse 0.0002009888  max resid 0.0004018609 
    ## ... Similar to previous best
    ## Run 13 stress 0.2848516 
    ## Run 14 stress 9.769001e-05 
    ## ... Procrustes: rmse 0.0001859105  max resid 0.0002619384 
    ## ... Similar to previous best
    ## Run 15 stress 9.289134e-05 
    ## ... Procrustes: rmse 0.0001821187  max resid 0.0003264508 
    ## ... Similar to previous best
    ## Run 16 stress 9.838609e-05 
    ## ... Procrustes: rmse 0.0002305501  max resid 0.0003458246 
    ## ... Similar to previous best
    ## Run 17 stress 9.92604e-05 
    ## ... Procrustes: rmse 0.0001907014  max resid 0.0003953594 
    ## ... Similar to previous best
    ## Run 18 stress 9.576265e-05 
    ## ... Procrustes: rmse 0.0002283832  max resid 0.0003591867 
    ## ... Similar to previous best
    ## Run 19 stress 9.572692e-05 
    ## ... Procrustes: rmse 0.0002048105  max resid 0.000296498 
    ## ... Similar to previous best
    ## Run 20 stress 9.696687e-05 
    ## ... Procrustes: rmse 0.0002376584  max resid 0.0003575486 
    ## ... Similar to previous best
    ## Run 21 stress 7.740591e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001206421  max resid 0.0002247499 
    ## ... Similar to previous best
    ## Run 22 stress 9.298519e-05 
    ## ... Procrustes: rmse 0.0001311102  max resid 0.0002973925 
    ## ... Similar to previous best
    ## Run 23 stress 9.613874e-05 
    ## ... Procrustes: rmse 0.0002524621  max resid 0.0005828623 
    ## ... Similar to previous best
    ## Run 24 stress 9.781328e-05 
    ## ... Procrustes: rmse 0.0002608856  max resid 0.0006030666 
    ## ... Similar to previous best
    ## Run 25 stress 9.228678e-05 
    ## ... Procrustes: rmse 0.0001853045  max resid 0.0003152468 
    ## ... Similar to previous best
    ## Run 26 stress 9.928965e-05 
    ## ... Procrustes: rmse 0.0001242574  max resid 0.0002306661 
    ## ... Similar to previous best
    ## Run 27 stress 9.898937e-05 
    ## ... Procrustes: rmse 0.0002134905  max resid 0.0003874508 
    ## ... Similar to previous best
    ## Run 28 stress 9.832545e-05 
    ## ... Procrustes: rmse 0.0002578395  max resid 0.0005995746 
    ## ... Similar to previous best
    ## Run 29 stress 8.999119e-05 
    ## ... Procrustes: rmse 0.0002368198  max resid 0.0003402902 
    ## ... Similar to previous best
    ## Run 30 stress 9.850533e-05 
    ## ... Procrustes: rmse 0.0001300307  max resid 0.000253713 
    ## ... Similar to previous best
    ## Run 31 stress 9.842074e-05 
    ## ... Procrustes: rmse 0.0002618017  max resid 0.0006043202 
    ## ... Similar to previous best
    ## Run 32 stress 9.836204e-05 
    ## ... Procrustes: rmse 0.0001987102  max resid 0.0003304265 
    ## ... Similar to previous best
    ## Run 33 stress 9.035296e-05 
    ## ... Procrustes: rmse 0.0002439975  max resid 0.0003571982 
    ## ... Similar to previous best
    ## Run 34 stress 9.231707e-05 
    ## ... Procrustes: rmse 0.0001294565  max resid 0.0002425462 
    ## ... Similar to previous best
    ## Run 35 stress 9.870816e-05 
    ## ... Procrustes: rmse 0.0002612364  max resid 0.0006024916 
    ## ... Similar to previous best
    ## Run 36 stress 0.3079616 
    ## Run 37 stress 9.58018e-05 
    ## ... Procrustes: rmse 0.0001919941  max resid 0.0003282917 
    ## ... Similar to previous best
    ## Run 38 stress 8.290847e-05 
    ## ... Procrustes: rmse 5.629958e-05  max resid 9.589249e-05 
    ## ... Similar to previous best
    ## Run 39 stress 9.942909e-05 
    ## ... Procrustes: rmse 5.071249e-05  max resid 7.316962e-05 
    ## ... Similar to previous best
    ## Run 40 stress 8.693031e-05 
    ## ... Procrustes: rmse 0.0001542204  max resid 0.0002286075 
    ## ... Similar to previous best
    ## Run 41 stress 9.819571e-05 
    ## ... Procrustes: rmse 0.000257197  max resid 0.0005918928 
    ## ... Similar to previous best
    ## Run 42 stress 7.350332e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001120679  max resid 0.0001645557 
    ## ... Similar to previous best
    ## Run 43 stress 9.727791e-05 
    ## ... Procrustes: rmse 0.0002468075  max resid 0.0004510982 
    ## ... Similar to previous best
    ## Run 44 stress 9.218717e-05 
    ## ... Procrustes: rmse 0.0001523795  max resid 0.0002239603 
    ## ... Similar to previous best
    ## Run 45 stress 9.596755e-05 
    ## ... Procrustes: rmse 0.0002074567  max resid 0.0004163368 
    ## ... Similar to previous best
    ## Run 46 stress 8.787392e-05 
    ## ... Procrustes: rmse 0.0001987785  max resid 0.0004047317 
    ## ... Similar to previous best
    ## Run 47 stress 9.736586e-05 
    ## ... Procrustes: rmse 0.000160045  max resid 0.0003137046 
    ## ... Similar to previous best
    ## Run 48 stress 9.610341e-05 
    ## ... Procrustes: rmse 0.0002400354  max resid 0.0004315027 
    ## ... Similar to previous best
    ## Run 49 stress 8.351963e-05 
    ## ... Procrustes: rmse 0.0001830049  max resid 0.0003204238 
    ## ... Similar to previous best
    ## Run 50 stress 8.570051e-05 
    ## ... Procrustes: rmse 0.0001472501  max resid 0.0002993989 
    ## ... Similar to previous best
    ## Run 51 stress 9.912863e-05 
    ## ... Procrustes: rmse 0.000172738  max resid 0.0003042099 
    ## ... Similar to previous best
    ## Run 52 stress 8.986104e-05 
    ## ... Procrustes: rmse 0.0001501988  max resid 0.0002215365 
    ## ... Similar to previous best
    ## Run 53 stress 9.602986e-05 
    ## ... Procrustes: rmse 0.000215809  max resid 0.0004349994 
    ## ... Similar to previous best
    ## Run 54 stress 9.692654e-05 
    ## ... Procrustes: rmse 0.0002210296  max resid 0.0004438312 
    ## ... Similar to previous best
    ## Run 55 stress 9.580384e-05 
    ## ... Procrustes: rmse 0.0001444442  max resid 0.00020875 
    ## ... Similar to previous best
    ## Run 56 stress 9.080204e-05 
    ## ... Procrustes: rmse 0.0002077216  max resid 0.0004224368 
    ## ... Similar to previous best
    ## Run 57 stress 9.664891e-05 
    ## ... Procrustes: rmse 0.0002408943  max resid 0.0004835216 
    ## ... Similar to previous best
    ## Run 58 stress 9.740698e-05 
    ## ... Procrustes: rmse 0.0001701869  max resid 0.0002689069 
    ## ... Similar to previous best
    ## Run 59 stress 8.865087e-05 
    ## ... Procrustes: rmse 0.0002253277  max resid 0.0004621288 
    ## ... Similar to previous best
    ## Run 60 stress 9.526678e-05 
    ## ... Procrustes: rmse 0.0002354456  max resid 0.0004273217 
    ## ... Similar to previous best
    ## Run 61 stress 9.493063e-05 
    ## ... Procrustes: rmse 0.0002356319  max resid 0.0004741474 
    ## ... Similar to previous best
    ## Run 62 stress 9.146936e-05 
    ## ... Procrustes: rmse 5.698593e-05  max resid 9.814282e-05 
    ## ... Similar to previous best
    ## Run 63 stress 9.850013e-05 
    ## ... Procrustes: rmse 0.0002246116  max resid 0.0004545008 
    ## ... Similar to previous best
    ## Run 64 stress 9.782408e-05 
    ## ... Procrustes: rmse 0.0001420393  max resid 0.0002410083 
    ## ... Similar to previous best
    ## Run 65 stress 9.911411e-05 
    ## ... Procrustes: rmse 0.0002273955  max resid 0.0004591043 
    ## ... Similar to previous best
    ## Run 66 stress 9.704676e-05 
    ## ... Procrustes: rmse 0.0002470854  max resid 0.0004984062 
    ## ... Similar to previous best
    ## Run 67 stress 8.345554e-05 
    ## ... Procrustes: rmse 0.000101748  max resid 0.0002027338 
    ## ... Similar to previous best
    ## Run 68 stress 8.584278e-05 
    ## ... Procrustes: rmse 0.0001566671  max resid 0.0002380628 
    ## ... Similar to previous best
    ## Run 69 stress 9.3341e-05 
    ## ... Procrustes: rmse 0.0002312247  max resid 0.0004182441 
    ## ... Similar to previous best
    ## Run 70 stress 9.978221e-05 
    ## ... Procrustes: rmse 0.0001216967  max resid 0.0002368854 
    ## ... Similar to previous best
    ## Run 71 stress 9.856765e-05 
    ## ... Procrustes: rmse 0.0001879381  max resid 0.0003440204 
    ## ... Similar to previous best
    ## Run 72 stress 9.819833e-05 
    ## ... Procrustes: rmse 0.0002240887  max resid 0.0004531689 
    ## ... Similar to previous best
    ## Run 73 stress 9.532068e-05 
    ## ... Procrustes: rmse 0.0002098173  max resid 0.0004271155 
    ## ... Similar to previous best
    ## Run 74 stress 9.34798e-05 
    ## ... Procrustes: rmse 0.0001548818  max resid 0.000227516 
    ## ... Similar to previous best
    ## Run 75 stress 9.732152e-05 
    ## ... Procrustes: rmse 0.000130532  max resid 0.0001919734 
    ## ... Similar to previous best
    ## Run 76 stress 9.186958e-05 
    ## ... Procrustes: rmse 0.0001400452  max resid 0.0003159765 
    ## ... Similar to previous best
    ## Run 77 stress 9.598455e-05 
    ## ... Procrustes: rmse 0.0002171613  max resid 0.0004411894 
    ## ... Similar to previous best
    ## Run 78 stress 9.540013e-05 
    ## ... Procrustes: rmse 0.0002377617  max resid 0.0004243799 
    ## ... Similar to previous best
    ## Run 79 stress 3.284437e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 8.071968e-05  max resid 0.0001429614 
    ## ... Similar to previous best
    ## Run 80 stress 9.748527e-05 
    ## ... Procrustes: rmse 0.0002034818  max resid 0.0003700431 
    ## ... Similar to previous best
    ## Run 81 stress 9.49473e-05 
    ## ... Procrustes: rmse 0.000143775  max resid 0.0002918032 
    ## ... Similar to previous best
    ## Run 82 stress 9.095393e-05 
    ## ... Procrustes: rmse 0.0001345364  max resid 0.0002222356 
    ## ... Similar to previous best
    ## Run 83 stress 9.606997e-05 
    ## ... Procrustes: rmse 0.0002035895  max resid 0.0003740582 
    ## ... Similar to previous best
    ## Run 84 stress 0.3083098 
    ## Run 85 stress 9.929352e-05 
    ## ... Procrustes: rmse 0.0001319604  max resid 0.0002421228 
    ## ... Similar to previous best
    ## Run 86 stress 9.479228e-05 
    ## ... Procrustes: rmse 0.0001173316  max resid 0.0002038634 
    ## ... Similar to previous best
    ## Run 87 stress 9.266122e-05 
    ## ... Procrustes: rmse 0.0001969162  max resid 0.000363481 
    ## ... Similar to previous best
    ## Run 88 stress 9.605067e-05 
    ## ... Procrustes: rmse 0.0002059789  max resid 0.0003788004 
    ## ... Similar to previous best
    ## Run 89 stress 9.360991e-05 
    ## ... Procrustes: rmse 0.000146966  max resid 0.0002365201 
    ## ... Similar to previous best
    ## Run 90 stress 9.364655e-05 
    ## ... Procrustes: rmse 0.0001187688  max resid 0.0002079328 
    ## ... Similar to previous best
    ## Run 91 stress 9.862315e-05 
    ## ... Procrustes: rmse 0.0001946837  max resid 0.0003533159 
    ## ... Similar to previous best
    ## Run 92 stress 9.81897e-05 
    ## ... Procrustes: rmse 0.0001218538  max resid 0.0001907061 
    ## ... Similar to previous best
    ## Run 93 stress 9.439795e-05 
    ## ... Procrustes: rmse 0.0001211268  max resid 0.0002159851 
    ## ... Similar to previous best
    ## Run 94 stress 9.122066e-05 
    ## ... Procrustes: rmse 0.0002076022  max resid 0.0003513723 
    ## ... Similar to previous best
    ## Run 95 stress 9.825562e-05 
    ## ... Procrustes: rmse 0.0002138745  max resid 0.000360219 
    ## ... Similar to previous best
    ## Run 96 stress 9.326165e-05 
    ## ... Procrustes: rmse 0.0001518962  max resid 0.0003151184 
    ## ... Similar to previous best
    ## Run 97 stress 9.819648e-05 
    ## ... Procrustes: rmse 0.0001207301  max resid 0.0001887953 
    ## ... Similar to previous best
    ## Run 98 stress 9.566776e-05 
    ## ... Procrustes: rmse 0.000182713  max resid 0.0002690947 
    ## ... Similar to previous best
    ## Run 99 stress 8.740808e-05 
    ## ... Procrustes: rmse 0.0001023631  max resid 0.0001557093 
    ## ... Similar to previous best
    ## Run 100 stress 9.438084e-05 
    ## ... Procrustes: rmse 0.0001981081  max resid 0.0003696766 
    ## ... Similar to previous best
    ## Run 101 stress 9.756799e-05 
    ## ... Procrustes: rmse 0.0001948658  max resid 0.0003415238 
    ## ... Similar to previous best
    ## Run 102 stress 9.786422e-05 
    ## ... Procrustes: rmse 0.0002047881  max resid 0.000378573 
    ## ... Similar to previous best
    ## Run 103 stress 9.540971e-05 
    ## ... Procrustes: rmse 0.0001974868  max resid 0.0003582823 
    ## ... Similar to previous best
    ## Run 104 stress 9.921485e-05 
    ## ... Procrustes: rmse 0.0002136648  max resid 0.0003899754 
    ## ... Similar to previous best
    ## Run 105 stress 7.11147e-05 
    ## ... Procrustes: rmse 0.000126718  max resid 0.0001930641 
    ## ... Similar to previous best
    ## Run 106 stress 9.233468e-05 
    ## ... Procrustes: rmse 0.0002075539  max resid 0.0003496814 
    ## ... Similar to previous best
    ## Run 107 stress 9.216063e-05 
    ## ... Procrustes: rmse 0.0001471453  max resid 0.0003008625 
    ## ... Similar to previous best
    ## Run 108 stress 9.809935e-05 
    ## ... Procrustes: rmse 0.0002037897  max resid 0.0003768797 
    ## ... Similar to previous best
    ## Run 109 stress 9.585247e-05 
    ## ... Procrustes: rmse 0.0001217835  max resid 0.0002082406 
    ## ... Similar to previous best
    ## Run 110 stress 0.2479796 
    ## Run 111 stress 9.25235e-05 
    ## ... Procrustes: rmse 0.000116237  max resid 0.0002074548 
    ## ... Similar to previous best
    ## Run 112 stress 8.532531e-05 
    ## ... Procrustes: rmse 0.0001269147  max resid 0.000229014 
    ## ... Similar to previous best
    ## Run 113 stress 9.181921e-05 
    ## ... Procrustes: rmse 0.0001366577  max resid 0.0002177181 
    ## ... Similar to previous best
    ## Run 114 stress 9.927568e-05 
    ## ... Procrustes: rmse 0.0002305677  max resid 0.0003923929 
    ## ... Similar to previous best
    ## Run 115 stress 9.370927e-05 
    ## ... Procrustes: rmse 0.0001933439  max resid 0.0003582852 
    ## ... Similar to previous best
    ## Run 116 stress 9.169623e-05 
    ## ... Procrustes: rmse 0.0001650566  max resid 0.0002508243 
    ## ... Similar to previous best
    ## Run 117 stress 9.664197e-05 
    ## ... Procrustes: rmse 0.0001515431  max resid 0.0003063252 
    ## ... Similar to previous best
    ## Run 118 stress 9.04044e-05 
    ## ... Procrustes: rmse 0.0002055639  max resid 0.0003467319 
    ## ... Similar to previous best
    ## Run 119 stress 0.2842076 
    ## Run 120 stress 9.890137e-05 
    ## ... Procrustes: rmse 0.0001413215  max resid 0.0002356724 
    ## ... Similar to previous best
    ## Run 121 stress 8.819577e-05 
    ## ... Procrustes: rmse 0.0002006183  max resid 0.0003401521 
    ## ... Similar to previous best
    ## Run 122 stress 9.533332e-05 
    ## ... Procrustes: rmse 0.0002026951  max resid 0.0003699404 
    ## ... Similar to previous best
    ## Run 123 stress 9.435142e-05 
    ## ... Procrustes: rmse 0.0001918114  max resid 0.0003408742 
    ## ... Similar to previous best
    ## Run 124 stress 9.013734e-05 
    ## ... Procrustes: rmse 0.0002013157  max resid 0.0003386197 
    ## ... Similar to previous best
    ## Run 125 stress 9.713915e-05 
    ## ... Procrustes: rmse 0.0002014679  max resid 0.0003675219 
    ## ... Similar to previous best
    ## Run 126 stress 9.285629e-05 
    ## ... Procrustes: rmse 0.0001939059  max resid 0.0003596313 
    ## ... Similar to previous best
    ## Run 127 stress 9.250687e-05 
    ## ... Procrustes: rmse 0.0001964727  max resid 0.0003631688 
    ## ... Similar to previous best
    ## Run 128 stress 9.307363e-05 
    ## ... Procrustes: rmse 0.000114657  max resid 0.0001995674 
    ## ... Similar to previous best
    ## Run 129 stress 9.979884e-05 
    ## ... Procrustes: rmse 0.000130531  max resid 0.0002332203 
    ## ... Similar to previous best
    ## Run 130 stress 9.907804e-05 
    ## ... Procrustes: rmse 0.0001291431  max resid 0.0002307192 
    ## ... Similar to previous best
    ## Run 131 stress 7.412977e-05 
    ## ... Procrustes: rmse 0.000106889  max resid 0.0002007777 
    ## ... Similar to previous best
    ## Run 132 stress 9.015352e-05 
    ## ... Procrustes: rmse 0.0001897378  max resid 0.0003475491 
    ## ... Similar to previous best
    ## Run 133 stress 9.8982e-05 
    ## ... Procrustes: rmse 0.0002031133  max resid 0.0003769952 
    ## ... Similar to previous best
    ## Run 134 stress 9.253517e-05 
    ## ... Procrustes: rmse 0.0001946813  max resid 0.0003576219 
    ## ... Similar to previous best
    ## Run 135 stress 9.672747e-05 
    ## ... Procrustes: rmse 0.0001386254  max resid 0.0002319594 
    ## ... Similar to previous best
    ## Run 136 stress 9.265011e-05 
    ## ... Procrustes: rmse 0.0001160808  max resid 0.0002028101 
    ## ... Similar to previous best
    ## Run 137 stress 9.342095e-05 
    ## ... Procrustes: rmse 0.0001445097  max resid 0.000225726 
    ## ... Similar to previous best
    ## Run 138 stress 9.361372e-05 
    ## ... Procrustes: rmse 0.0001607429  max resid 0.0002318932 
    ## ... Similar to previous best
    ## Run 139 stress 8.736837e-05 
    ## ... Procrustes: rmse 0.0001313575  max resid 0.000242157 
    ## ... Similar to previous best
    ## Run 140 stress 9.752218e-05 
    ## ... Procrustes: rmse 0.0002076214  max resid 0.0003815039 
    ## ... Similar to previous best
    ## Run 141 stress 8.995675e-05 
    ## ... Procrustes: rmse 0.0001855816  max resid 0.0003444562 
    ## ... Similar to previous best
    ## Run 142 stress 9.087117e-05 
    ## ... Procrustes: rmse 0.0001418108  max resid 0.0003025012 
    ## ... Similar to previous best
    ## Run 143 stress 9.606945e-05 
    ## ... Procrustes: rmse 0.0001588206  max resid 0.0003272317 
    ## ... Similar to previous best
    ## Run 144 stress 9.383468e-05 
    ## ... Procrustes: rmse 0.0001903864  max resid 0.0003451595 
    ## ... Similar to previous best
    ## Run 145 stress 8.887701e-05 
    ## ... Procrustes: rmse 0.0001102493  max resid 0.0001762256 
    ## ... Similar to previous best
    ## Run 146 stress 8.869318e-05 
    ## ... Procrustes: rmse 0.0001651871  max resid 0.0002277749 
    ## ... Similar to previous best
    ## Run 147 stress 9.306133e-05 
    ## ... Procrustes: rmse 0.0001575727  max resid 0.0002454381 
    ## ... Similar to previous best
    ## Run 148 stress 9.718647e-05 
    ## ... Procrustes: rmse 0.0002054687  max resid 0.0003738245 
    ## ... Similar to previous best
    ## Run 149 stress 9.675978e-05 
    ## ... Procrustes: rmse 0.0001453114  max resid 0.0002938088 
    ## ... Similar to previous best
    ## Run 150 stress 0.3076939 
    ## Run 151 stress 9.774365e-05 
    ## ... Procrustes: rmse 0.0001258313  max resid 0.0002248737 
    ## ... Similar to previous best
    ## Run 152 stress 9.197854e-05 
    ## ... Procrustes: rmse 0.0002058735  max resid 0.0003502588 
    ## ... Similar to previous best
    ## Run 153 stress 9.537279e-05 
    ## ... Procrustes: rmse 0.0002026767  max resid 0.0003699916 
    ## ... Similar to previous best
    ## Run 154 stress 8.859996e-05 
    ## ... Procrustes: rmse 0.0001322307  max resid 0.000268868 
    ## ... Similar to previous best
    ## Run 155 stress 9.974816e-05 
    ## ... Procrustes: rmse 0.0001699506  max resid 0.0002753071 
    ## ... Similar to previous best
    ## Run 156 stress 0.231921 
    ## Run 157 stress 9.398687e-05 
    ## ... Procrustes: rmse 0.0001631637  max resid 0.0003261075 
    ## ... Similar to previous best
    ## Run 158 stress 9.354317e-05 
    ## ... Procrustes: rmse 0.0001989839  max resid 0.0003677691 
    ## ... Similar to previous best
    ## Run 159 stress 9.693027e-05 
    ## ... Procrustes: rmse 0.0002057179  max resid 0.0003795376 
    ## ... Similar to previous best
    ## Run 160 stress 9.492021e-05 
    ## ... Procrustes: rmse 0.0001204984  max resid 0.0002105334 
    ## ... Similar to previous best
    ## Run 161 stress 9.825602e-05 
    ## ... Procrustes: rmse 0.0002026341  max resid 0.0003084902 
    ## ... Similar to previous best
    ## Run 162 stress 9.576047e-05 
    ## ... Procrustes: rmse 0.0002029525  max resid 0.0003745799 
    ## ... Similar to previous best
    ## Run 163 stress 9.396779e-05 
    ## ... Procrustes: rmse 0.0001192839  max resid 0.0002129345 
    ## ... Similar to previous best
    ## Run 164 stress 9.939551e-05 
    ## ... Procrustes: rmse 0.0001284336  max resid 0.0002278679 
    ## ... Similar to previous best
    ## Run 165 stress 9.526255e-05 
    ## ... Procrustes: rmse 0.0001218627  max resid 0.0002134786 
    ## ... Similar to previous best
    ## Run 166 stress 9.109041e-05 
    ## ... Procrustes: rmse 0.0001134242  max resid 0.0002023829 
    ## ... Similar to previous best
    ## Run 167 stress 0.3120129 
    ## Run 168 stress 9.995307e-05 
    ## ... Procrustes: rmse 0.0001307158  max resid 0.0002335683 
    ## ... Similar to previous best
    ## Run 169 stress 9.053485e-05 
    ## ... Procrustes: rmse 0.0001576306  max resid 0.0003144223 
    ## ... Similar to previous best
    ## Run 170 stress 9.074218e-05 
    ## ... Procrustes: rmse 0.0001132627  max resid 0.000202105 
    ## ... Similar to previous best
    ## Run 171 stress 9.647928e-05 
    ## ... Procrustes: rmse 0.0001241284  max resid 0.0002216745 
    ## ... Similar to previous best
    ## Run 172 stress 8.956623e-05 
    ## ... Procrustes: rmse 0.0001168741  max resid 0.0001657611 
    ## ... Similar to previous best
    ## Run 173 stress 9.818913e-05 
    ## ... Procrustes: rmse 0.0001665644  max resid 0.0002622787 
    ## ... Similar to previous best
    ## Run 174 stress 6.916564e-05 
    ## ... Procrustes: rmse 0.0001174704  max resid 0.0002503863 
    ## ... Similar to previous best
    ## Run 175 stress 8.633128e-05 
    ## ... Procrustes: rmse 0.0001061421  max resid 0.0001709689 
    ## ... Similar to previous best
    ## Run 176 stress 9.794048e-05 
    ## ... Procrustes: rmse 0.00020521  max resid 0.0003795125 
    ## ... Similar to previous best
    ## Run 177 stress 8.745173e-05 
    ## ... Procrustes: rmse 0.0001233058  max resid 0.0002342495 
    ## ... Similar to previous best
    ## Run 178 stress 9.164992e-05 
    ## ... Procrustes: rmse 0.0001492989  max resid 0.0003094628 
    ## ... Similar to previous best
    ## Run 179 stress 9.858929e-05 
    ## ... Procrustes: rmse 0.0002278532  max resid 0.0003876957 
    ## ... Similar to previous best
    ## Run 180 stress 0.2848515 
    ## Run 181 stress 9.401179e-05 
    ## ... Procrustes: rmse 0.0001484133  max resid 0.0003012325 
    ## ... Similar to previous best
    ## Run 182 stress 9.818375e-05 
    ## ... Procrustes: rmse 0.0002107472  max resid 0.0003887339 
    ## ... Similar to previous best
    ## Run 183 stress 8.61062e-05 
    ## ... Procrustes: rmse 0.0001410573  max resid 0.0003036341 
    ## ... Similar to previous best
    ## Run 184 stress 9.667229e-05 
    ## ... Procrustes: rmse 0.0001327755  max resid 0.0002382516 
    ## ... Similar to previous best
    ## Run 185 stress 9.141597e-05 
    ## ... Procrustes: rmse 0.0001349037  max resid 0.000250302 
    ## ... Similar to previous best
    ## Run 186 stress 9.851245e-05 
    ## ... Procrustes: rmse 0.0002096494  max resid 0.0003868236 
    ## ... Similar to previous best
    ## Run 187 stress 9.841436e-05 
    ## ... Procrustes: rmse 0.0002042936  max resid 0.000371562 
    ## ... Similar to previous best
    ## Run 188 stress 9.056799e-05 
    ## ... Procrustes: rmse 0.0001502571  max resid 0.0002504111 
    ## ... Similar to previous best
    ## Run 189 stress 9.651844e-05 
    ## ... Procrustes: rmse 0.0002010439  max resid 0.0003655426 
    ## ... Similar to previous best
    ## Run 190 stress 9.435942e-05 
    ## ... Procrustes: rmse 0.0002170916  max resid 0.0003684109 
    ## ... Similar to previous best
    ## Run 191 stress 9.049859e-05 
    ## ... Procrustes: rmse 0.0001252286  max resid 0.0002433246 
    ## ... Similar to previous best
    ## Run 192 stress 9.772265e-05 
    ## ... Procrustes: rmse 0.0001403158  max resid 0.0002171261 
    ## ... Similar to previous best
    ## Run 193 stress 9.168091e-05 
    ## ... Procrustes: rmse 0.0001152781  max resid 0.0001912076 
    ## ... Similar to previous best
    ## Run 194 stress 8.474646e-05 
    ## ... Procrustes: rmse 0.0001045264  max resid 0.000168579 
    ## ... Similar to previous best
    ## Run 195 stress 8.817562e-05 
    ## ... Procrustes: rmse 0.0001342579  max resid 0.0002430252 
    ## ... Similar to previous best
    ## Run 196 stress 9.067025e-05 
    ## ... Procrustes: rmse 0.0001793463  max resid 0.0003244215 
    ## ... Similar to previous best
    ## Run 197 stress 9.568061e-05 
    ## ... Procrustes: rmse 0.0002028121  max resid 0.0003696858 
    ## ... Similar to previous best
    ## Run 198 stress 9.808298e-05 
    ## ... Procrustes: rmse 0.0001503621  max resid 0.0003195034 
    ## ... Similar to previous best
    ## Run 199 stress 9.018877e-05 
    ## ... Procrustes: rmse 0.0001718269  max resid 0.000250486 
    ## ... Similar to previous best
    ## Run 200 stress 8.706282e-05 
    ## ... Procrustes: rmse 0.000142314  max resid 0.0003039435 
    ## ... Similar to previous best
    ## Run 201 stress 9.883291e-05 
    ## ... Procrustes: rmse 0.000206546  max resid 0.0003812189 
    ## ... Similar to previous best
    ## Run 202 stress 9.950035e-05 
    ## ... Procrustes: rmse 0.0001859949  max resid 0.0003322379 
    ## ... Similar to previous best
    ## Run 203 stress 9.7828e-05 
    ## ... Procrustes: rmse 0.0001204639  max resid 0.0002183464 
    ## ... Similar to previous best
    ## Run 204 stress 8.460113e-05 
    ## ... Procrustes: rmse 0.0001009942  max resid 0.0001780967 
    ## ... Similar to previous best
    ## Run 205 stress 9.640655e-05 
    ## ... Procrustes: rmse 0.0002108402  max resid 0.0003506361 
    ## ... Similar to previous best
    ## Run 206 stress 9.678562e-05 
    ## ... Procrustes: rmse 0.0002030558  max resid 0.00036968 
    ## ... Similar to previous best
    ## Run 207 stress 8.817864e-05 
    ## ... Procrustes: rmse 0.0001372168  max resid 0.0002800315 
    ## ... Similar to previous best
    ## Run 208 stress 9.42456e-05 
    ## ... Procrustes: rmse 0.0001483984  max resid 0.0003090627 
    ## ... Similar to previous best
    ## Run 209 stress 0.3083098 
    ## Run 210 stress 9.167707e-05 
    ## ... Procrustes: rmse 0.0001664266  max resid 0.000249284 
    ## ... Similar to previous best
    ## Run 211 stress 9.972163e-05 
    ## ... Procrustes: rmse 0.0001262344  max resid 0.0002191456 
    ## ... Similar to previous best
    ## Run 212 stress 9.470667e-05 
    ## ... Procrustes: rmse 0.0002187276  max resid 0.0003689821 
    ## ... Similar to previous best
    ## Run 213 stress 9.253223e-05 
    ## ... Procrustes: rmse 0.0001503873  max resid 0.0003218525 
    ## ... Similar to previous best
    ## Run 214 stress 8.890032e-05 
    ## ... Procrustes: rmse 0.0001993642  max resid 0.0003341329 
    ## ... Similar to previous best
    ## Run 215 stress 0.231921 
    ## Run 216 stress 0.3081919 
    ## Run 217 stress 7.737875e-05 
    ## ... Procrustes: rmse 0.0001243512  max resid 0.0002662319 
    ## ... Similar to previous best
    ## Run 218 stress 8.588199e-05 
    ## ... Procrustes: rmse 0.000132746  max resid 0.0002184652 
    ## ... Similar to previous best
    ## Run 219 stress 9.7281e-05 
    ## ... Procrustes: rmse 0.0001818026  max resid 0.0002501082 
    ## ... Similar to previous best
    ## Run 220 stress 9.200691e-05 
    ## ... Procrustes: rmse 0.0001317965  max resid 0.0002217353 
    ## ... Similar to previous best
    ## Run 221 stress 9.490961e-05 
    ## ... Procrustes: rmse 0.00021885  max resid 0.0003664039 
    ## ... Similar to previous best
    ## Run 222 stress 8.670288e-05 
    ## ... Procrustes: rmse 0.0001077374  max resid 0.0001893389 
    ## ... Similar to previous best
    ## Run 223 stress 9.537783e-05 
    ## ... Procrustes: rmse 0.0001882159  max resid 0.0003394052 
    ## ... Similar to previous best
    ## Run 224 stress 9.716263e-05 
    ## ... Procrustes: rmse 0.0002044026  max resid 0.0003721094 
    ## ... Similar to previous best
    ## Run 225 stress 9.255131e-05 
    ## ... Procrustes: rmse 0.0001729074  max resid 0.0002479015 
    ## ... Similar to previous best
    ## Run 226 stress 9.589552e-05 
    ## ... Procrustes: rmse 0.000200375  max resid 0.0003711539 
    ## ... Similar to previous best
    ## Run 227 stress 8.389937e-05 
    ## ... Procrustes: rmse 0.0001246017  max resid 0.0001967835 
    ## ... Similar to previous best
    ## Run 228 stress 9.765914e-05 
    ## ... Procrustes: rmse 0.0002017113  max resid 0.0003653646 
    ## ... Similar to previous best
    ## Run 229 stress 9.488489e-05 
    ## ... Procrustes: rmse 0.0001208816  max resid 0.0002115742 
    ## ... Similar to previous best
    ## Run 230 stress 9.802717e-05 
    ## ... Procrustes: rmse 0.0001985133  max resid 0.0003371151 
    ## ... Similar to previous best
    ## Run 231 stress 9.244929e-05 
    ## ... Procrustes: rmse 0.0001361947  max resid 0.0002737849 
    ## ... Similar to previous best
    ## Run 232 stress 9.841722e-05 
    ## ... Procrustes: rmse 0.0002204991  max resid 0.0003686833 
    ## ... Similar to previous best
    ## Run 233 stress 9.127215e-05 
    ## ... Procrustes: rmse 0.0001820914  max resid 0.0003356927 
    ## ... Similar to previous best
    ## Run 234 stress 9.96188e-05 
    ## ... Procrustes: rmse 0.0001531679  max resid 0.0002590914 
    ## ... Similar to previous best
    ## Run 235 stress 9.945729e-05 
    ## ... Procrustes: rmse 0.0001950757  max resid 0.0003554936 
    ## ... Similar to previous best
    ## Run 236 stress 9.721315e-05 
    ## ... Procrustes: rmse 0.000124842  max resid 0.0002183983 
    ## ... Similar to previous best
    ## Run 237 stress 9.759095e-05 
    ## ... Procrustes: rmse 0.0002012604  max resid 0.0003757231 
    ## ... Similar to previous best
    ## Run 238 stress 9.322299e-05 
    ## ... Procrustes: rmse 0.0001205445  max resid 0.0002106997 
    ## ... Similar to previous best
    ## Run 239 stress 9.649733e-05 
    ## ... Procrustes: rmse 0.0002038482  max resid 0.0003766754 
    ## ... Similar to previous best
    ## Run 240 stress 8.994297e-05 
    ## ... Procrustes: rmse 0.0001892251  max resid 0.0003494287 
    ## ... Similar to previous best
    ## Run 241 stress 7.485135e-05 
    ## ... Procrustes: rmse 0.0001138215  max resid 0.0002426177 
    ## ... Similar to previous best
    ## Run 242 stress 9.933084e-05 
    ## ... Procrustes: rmse 0.0002257226  max resid 0.000379694 
    ## ... Similar to previous best
    ## Run 243 stress 0.3081116 
    ## Run 244 stress 9.025717e-05 
    ## ... Procrustes: rmse 0.0001127853  max resid 0.000179768 
    ## ... Similar to previous best
    ## Run 245 stress 9.87131e-05 
    ## ... Procrustes: rmse 0.0001997789  max resid 0.0003739438 
    ## ... Similar to previous best
    ## Run 246 stress 9.484583e-05 
    ## ... Procrustes: rmse 0.000215175  max resid 0.0003644975 
    ## ... Similar to previous best
    ## Run 247 stress 8.964958e-05 
    ## ... Procrustes: rmse 0.0001446255  max resid 0.0002979546 
    ## ... Similar to previous best
    ## Run 248 stress 8.877207e-05 
    ## ... Procrustes: rmse 0.0001051122  max resid 0.0001792325 
    ## ... Similar to previous best
    ## Run 249 stress 9.841412e-05 
    ## ... Procrustes: rmse 0.0002023983  max resid 0.0003565679 
    ## ... Similar to previous best
    ## Run 250 stress 9.884496e-05 
    ## ... Procrustes: rmse 0.000230098  max resid 0.0003917931 
    ## ... Similar to previous best
    ## Run 251 stress 9.219751e-05 
    ## ... Procrustes: rmse 0.0001958712  max resid 0.0003621286 
    ## ... Similar to previous best
    ## Run 252 stress 9.776369e-05 
    ## ... Procrustes: rmse 0.0001194317  max resid 0.0002074698 
    ## ... Similar to previous best
    ## Run 253 stress 9.348823e-05 
    ## ... Procrustes: rmse 0.0001155124  max resid 0.0002053518 
    ## ... Similar to previous best
    ## Run 254 stress 9.588075e-05 
    ## ... Procrustes: rmse 0.0001215474  max resid 0.0002043085 
    ## ... Similar to previous best
    ## Run 255 stress 0.231917 
    ## Run 256 stress 9.617142e-05 
    ## ... Procrustes: rmse 0.0002057125  max resid 0.0003796839 
    ## ... Similar to previous best
    ## Run 257 stress 7.223799e-05 
    ## ... Procrustes: rmse 7.844351e-05  max resid 0.0001342378 
    ## ... Similar to previous best
    ## Run 258 stress 9.70533e-05 
    ## ... Procrustes: rmse 0.0001284371  max resid 0.0002090236 
    ## ... Similar to previous best
    ## Run 259 stress 9.980922e-05 
    ## ... Procrustes: rmse 0.0001285386  max resid 0.0002274309 
    ## ... Similar to previous best
    ## Run 260 stress 9.828396e-05 
    ## ... Procrustes: rmse 0.0001644658  max resid 0.0002692127 
    ## ... Similar to previous best
    ## Run 261 stress 9.69315e-05 
    ## ... Procrustes: rmse 0.0001589806  max resid 0.0002399895 
    ## ... Similar to previous best
    ## Run 262 stress 8.058692e-05 
    ## ... Procrustes: rmse 0.0001123898  max resid 0.0001950943 
    ## ... Similar to previous best
    ## Run 263 stress 9.656827e-05 
    ## ... Procrustes: rmse 0.0002171155  max resid 0.0003697232 
    ## ... Similar to previous best
    ## Run 264 stress 8.54145e-05 
    ## ... Procrustes: rmse 0.0001389146  max resid 0.0002882212 
    ## ... Similar to previous best
    ## Run 265 stress 9.162699e-05 
    ## ... Procrustes: rmse 0.000118329  max resid 0.0001947244 
    ## ... Similar to previous best
    ## Run 266 stress 9.890602e-05 
    ## ... Procrustes: rmse 0.0001110709  max resid 0.0001871358 
    ## ... Similar to previous best
    ## Run 267 stress 9.727454e-05 
    ## ... Procrustes: rmse 0.0002179906  max resid 0.0003741278 
    ## ... Similar to previous best
    ## Run 268 stress 9.773857e-05 
    ## ... Procrustes: rmse 0.0001060108  max resid 0.0001644923 
    ## ... Similar to previous best
    ## Run 269 stress 9.646546e-05 
    ## ... Procrustes: rmse 0.000145063  max resid 0.0002864465 
    ## ... Similar to previous best
    ## Run 270 stress 9.935788e-05 
    ## ... Procrustes: rmse 0.0002141876  max resid 0.0003922274 
    ## ... Similar to previous best
    ## Run 271 stress 9.742403e-05 
    ## ... Procrustes: rmse 0.0002033778  max resid 0.0003739108 
    ## ... Similar to previous best
    ## Run 272 stress 9.476396e-05 
    ## ... Procrustes: rmse 0.000200299  max resid 0.0003640989 
    ## ... Similar to previous best
    ## Run 273 stress 9.685894e-05 
    ## ... Procrustes: rmse 0.0002062861  max resid 0.0003760697 
    ## ... Similar to previous best
    ## Run 274 stress 9.669767e-05 
    ## ... Procrustes: rmse 0.0002034541  max resid 0.0003774268 
    ## ... Similar to previous best
    ## Run 275 stress 9.53241e-05 
    ## ... Procrustes: rmse 0.0001986191  max resid 0.0003257533 
    ## ... Similar to previous best
    ## Run 276 stress 8.99273e-05 
    ## ... Procrustes: rmse 0.0001843729  max resid 0.0003438896 
    ## ... Similar to previous best
    ## Run 277 stress 5.922406e-05 
    ## ... Procrustes: rmse 8.940275e-05  max resid 0.0001991788 
    ## ... Similar to previous best
    ## Run 278 stress 9.425401e-05 
    ## ... Procrustes: rmse 0.0001922555  max resid 0.0003547731 
    ## ... Similar to previous best
    ## Run 279 stress 9.857906e-05 
    ## ... Procrustes: rmse 0.0001456671  max resid 0.0002448486 
    ## ... Similar to previous best
    ## Run 280 stress 8.857475e-05 
    ## ... Procrustes: rmse 0.0001419199  max resid 0.0003011057 
    ## ... Similar to previous best
    ## Run 281 stress 9.299155e-05 
    ## ... Procrustes: rmse 0.0001955257  max resid 0.0003531893 
    ## ... Similar to previous best
    ## Run 282 stress 9.479936e-05 
    ## ... Procrustes: rmse 0.000153922  max resid 0.0003286105 
    ## ... Similar to previous best
    ## Run 283 stress 9.432159e-05 
    ## ... Procrustes: rmse 0.0002003661  max resid 0.0003688255 
    ## ... Similar to previous best
    ## Run 284 stress 9.142472e-05 
    ## ... Procrustes: rmse 0.0001669556  max resid 0.0002278383 
    ## ... Similar to previous best
    ## Run 285 stress 9.366424e-05 
    ## ... Procrustes: rmse 0.0001993682  max resid 0.0003677927 
    ## ... Similar to previous best
    ## Run 286 stress 9.987551e-05 
    ## ... Procrustes: rmse 0.0002153128  max resid 0.0003937632 
    ## ... Similar to previous best
    ## Run 287 stress 9.04447e-05 
    ## ... Procrustes: rmse 0.0001913177  max resid 0.0003536263 
    ## ... Similar to previous best
    ## Run 288 stress 9.876887e-05 
    ## ... Procrustes: rmse 0.0002112285  max resid 0.00038491 
    ## ... Similar to previous best
    ## Run 289 stress 9.788661e-05 
    ## ... Procrustes: rmse 0.0001228815  max resid 0.0001971943 
    ## ... Similar to previous best
    ## Run 290 stress 9.95968e-05 
    ## ... Procrustes: rmse 0.0001581107  max resid 0.0003209108 
    ## ... Similar to previous best
    ## Run 291 stress 8.753335e-05 
    ## ... Procrustes: rmse 0.0001404191  max resid 0.0003027026 
    ## ... Similar to previous best
    ## Run 292 stress 0.3120129 
    ## Run 293 stress 9.492241e-05 
    ## ... Procrustes: rmse 0.0001332902  max resid 0.0002226146 
    ## ... Similar to previous best
    ## Run 294 stress 8.030551e-05 
    ## ... Procrustes: rmse 0.0001082367  max resid 0.0001989496 
    ## ... Similar to previous best
    ## Run 295 stress 8.968479e-05 
    ## ... Procrustes: rmse 0.0001266261  max resid 0.0002350608 
    ## ... Similar to previous best
    ## Run 296 stress 9.291213e-05 
    ## ... Procrustes: rmse 0.0001957271  max resid 0.0003570045 
    ## ... Similar to previous best
    ## Run 297 stress 7.95295e-05 
    ## ... Procrustes: rmse 0.0001331752  max resid 0.0002835173 
    ## ... Similar to previous best
    ## Run 298 stress 9.369221e-05 
    ## ... Procrustes: rmse 0.000155889  max resid 0.0002575857 
    ## ... Similar to previous best
    ## Run 299 stress 9.550893e-05 
    ## ... Procrustes: rmse 0.0001145279  max resid 0.0002020364 
    ## ... Similar to previous best
    ## Run 300 stress 9.789367e-05 
    ## ... Procrustes: rmse 0.0001289259  max resid 0.0002299905 
    ## ... Similar to previous best
    ## Run 301 stress 9.581348e-05 
    ## ... Procrustes: rmse 0.0001790571  max resid 0.0002562973 
    ## ... Similar to previous best
    ## Run 302 stress 9.536543e-05 
    ## ... Procrustes: rmse 0.0001217137  max resid 0.0002173358 
    ## ... Similar to previous best
    ## Run 303 stress 9.564395e-05 
    ## ... Procrustes: rmse 0.0002191607  max resid 0.0003682274 
    ## ... Similar to previous best
    ## Run 304 stress 0.2479793 
    ## Run 305 stress 8.628068e-05 
    ## ... Procrustes: rmse 0.0001412571  max resid 0.000304119 
    ## ... Similar to previous best
    ## Run 306 stress 9.501862e-05 
    ## ... Procrustes: rmse 0.0002004776  max resid 0.0003618953 
    ## ... Similar to previous best
    ## Run 307 stress 9.349278e-05 
    ## ... Procrustes: rmse 0.0001470125  max resid 0.0002468207 
    ## ... Similar to previous best
    ## Run 308 stress 9.754102e-05 
    ## ... Procrustes: rmse 0.0002037481  max resid 0.0003765622 
    ## ... Similar to previous best
    ## Run 309 stress 9.627554e-05 
    ## ... Procrustes: rmse 0.000196094  max resid 0.0003548702 
    ## ... Similar to previous best
    ## Run 310 stress 8.931345e-05 
    ## ... Procrustes: rmse 0.0001249464  max resid 0.0001750038 
    ## ... Similar to previous best
    ## Run 311 stress 9.763909e-05 
    ## ... Procrustes: rmse 0.0002021822  max resid 0.0003742705 
    ## ... Similar to previous best
    ## Run 312 stress 9.669119e-05 
    ## ... Procrustes: rmse 0.000147778  max resid 0.0002997782 
    ## ... Similar to previous best
    ## Run 313 stress 9.431466e-05 
    ## ... Procrustes: rmse 0.0001226058  max resid 0.0002009472 
    ## ... Similar to previous best
    ## Run 314 stress 9.135869e-05 
    ## ... Procrustes: rmse 0.0001500314  max resid 0.0003188597 
    ## ... Similar to previous best
    ## Run 315 stress 9.972366e-05 
    ## ... Procrustes: rmse 0.000162825  max resid 0.0003427273 
    ## ... Similar to previous best
    ## Run 316 stress 9.015557e-05 
    ## ... Procrustes: rmse 0.0002011189  max resid 0.0003417712 
    ## ... Similar to previous best
    ## Run 317 stress 8.626633e-05 
    ## ... Procrustes: rmse 0.0001314766  max resid 0.0002308079 
    ## ... Similar to previous best
    ## Run 318 stress 9.066457e-05 
    ## ... Procrustes: rmse 0.000126456  max resid 0.0002097046 
    ## ... Similar to previous best
    ## Run 319 stress 9.189147e-05 
    ## ... Procrustes: rmse 0.0001135829  max resid 0.0001803708 
    ## ... Similar to previous best
    ## Run 320 stress 8.629769e-05 
    ## ... Procrustes: rmse 0.0001409423  max resid 0.0003007349 
    ## ... Similar to previous best
    ## Run 321 stress 9.423969e-05 
    ## ... Procrustes: rmse 0.0001998909  max resid 0.0003642569 
    ## ... Similar to previous best
    ## Run 322 stress 9.095966e-05 
    ## ... Procrustes: rmse 0.0001335482  max resid 0.0002208495 
    ## ... Similar to previous best
    ## Run 323 stress 9.634341e-05 
    ## ... Procrustes: rmse 0.0001239719  max resid 0.0002171987 
    ## ... Similar to previous best
    ## Run 324 stress 9.820002e-05 
    ## ... Procrustes: rmse 0.0002041927  max resid 0.0003772917 
    ## ... Similar to previous best
    ## Run 325 stress 9.985176e-05 
    ## ... Procrustes: rmse 0.0001297457  max resid 0.0002269494 
    ## ... Similar to previous best
    ## Run 326 stress 9.465373e-05 
    ## ... Procrustes: rmse 0.0001564811  max resid 0.0002415213 
    ## ... Similar to previous best
    ## Run 327 stress 9.406521e-05 
    ## ... Procrustes: rmse 0.0001225277  max resid 0.000218615 
    ## ... Similar to previous best
    ## Run 328 stress 9.60016e-05 
    ## ... Procrustes: rmse 0.0002210947  max resid 0.0003711995 
    ## ... Similar to previous best
    ## Run 329 stress 9.269267e-05 
    ## ... Procrustes: rmse 0.0001163612  max resid 0.0002033984 
    ## ... Similar to previous best
    ## Run 330 stress 9.889179e-05 
    ## ... Procrustes: rmse 0.0002109986  max resid 0.0003840858 
    ## ... Similar to previous best
    ## Run 331 stress 9.638921e-05 
    ## ... Procrustes: rmse 0.0002064698  max resid 0.00038069 
    ## ... Similar to previous best
    ## Run 332 stress 9.355137e-05 
    ## ... Procrustes: rmse 0.0001555842  max resid 0.0002361559 
    ## ... Similar to previous best
    ## Run 333 stress 0.2848516 
    ## Run 334 stress 9.752328e-05 
    ## ... Procrustes: rmse 0.0001195104  max resid 0.0002116235 
    ## ... Similar to previous best
    ## Run 335 stress 9.694616e-05 
    ## ... Procrustes: rmse 0.0002166018  max resid 0.0003610114 
    ## ... Similar to previous best
    ## Run 336 stress 9.052145e-05 
    ## ... Procrustes: rmse 0.0001352155  max resid 0.0002730171 
    ## ... Similar to previous best
    ## Run 337 stress 9.397518e-05 
    ## ... Procrustes: rmse 0.000122808  max resid 0.0002183294 
    ## ... Similar to previous best
    ## Run 338 stress 8.278741e-05 
    ## ... Procrustes: rmse 0.0001017297  max resid 0.000223042 
    ## ... Similar to previous best
    ## Run 339 stress 9.57679e-05 
    ## ... Procrustes: rmse 0.0001218522  max resid 0.000212558 
    ## ... Similar to previous best
    ## Run 340 stress 9.775722e-05 
    ## ... Procrustes: rmse 0.0001261295  max resid 0.0002207944 
    ## ... Similar to previous best
    ## Run 341 stress 9.664222e-05 
    ## ... Procrustes: rmse 0.0001248607  max resid 0.000223053 
    ## ... Similar to previous best
    ## Run 342 stress 9.643815e-05 
    ## ... Procrustes: rmse 0.000221007  max resid 0.0003711193 
    ## ... Similar to previous best
    ## Run 343 stress 9.387938e-05 
    ## ... Procrustes: rmse 0.0001405874  max resid 0.0002109238 
    ## ... Similar to previous best
    ## Run 344 stress 9.969497e-05 
    ## ... Procrustes: rmse 0.0002021746  max resid 0.0003659355 
    ## ... Similar to previous best
    ## Run 345 stress 9.905665e-05 
    ## ... Procrustes: rmse 0.0002124515  max resid 0.000387455 
    ## ... Similar to previous best
    ## Run 346 stress 9.371931e-05 
    ## ... Procrustes: rmse 0.0001215831  max resid 0.000214101 
    ## ... Similar to previous best
    ## Run 347 stress 9.9023e-05 
    ## ... Procrustes: rmse 0.0002118906  max resid 0.0003889587 
    ## ... Similar to previous best
    ## Run 348 stress 9.861281e-05 
    ## ... Procrustes: rmse 0.0001529819  max resid 0.0002592291 
    ## ... Similar to previous best
    ## Run 349 stress 9.562988e-05 
    ## ... Procrustes: rmse 0.0001541293  max resid 0.0003162336 
    ## ... Similar to previous best
    ## Run 350 stress 9.389019e-05 
    ## ... Procrustes: rmse 0.0001097799  max resid 0.0001870942 
    ## ... Similar to previous best
    ## Run 351 stress 9.916735e-05 
    ## ... Procrustes: rmse 0.000128613  max resid 0.0002251138 
    ## ... Similar to previous best
    ## Run 352 stress 9.819317e-05 
    ## ... Procrustes: rmse 0.0002079946  max resid 0.000378228 
    ## ... Similar to previous best
    ## Run 353 stress 4.949832e-05 
    ## ... Procrustes: rmse 8.626472e-05  max resid 0.0001888523 
    ## ... Similar to previous best
    ## Run 354 stress 9.923577e-05 
    ## ... Procrustes: rmse 0.0002117637  max resid 0.0003915779 
    ## ... Similar to previous best
    ## Run 355 stress 9.706861e-05 
    ## ... Procrustes: rmse 0.0002073653  max resid 0.0003835649 
    ## ... Similar to previous best
    ## Run 356 stress 9.803468e-05 
    ## ... Procrustes: rmse 0.0001559964  max resid 0.0002637385 
    ## ... Similar to previous best
    ## Run 357 stress 0.2479794 
    ## Run 358 stress 9.787677e-05 
    ## ... Procrustes: rmse 0.0001278334  max resid 0.0001811958 
    ## ... Similar to previous best
    ## Run 359 stress 9.465728e-05 
    ## ... Procrustes: rmse 0.0001191012  max resid 0.00020781 
    ## ... Similar to previous best
    ## Run 360 stress 9.941228e-05 
    ## ... Procrustes: rmse 0.0002117931  max resid 0.0003903921 
    ## ... Similar to previous best
    ## Run 361 stress 9.920745e-05 
    ## ... Procrustes: rmse 0.0002113569  max resid 0.0003916672 
    ## ... Similar to previous best
    ## Run 362 stress 8.872886e-05 
    ## ... Procrustes: rmse 0.0002029321  max resid 0.000339113 
    ## ... Similar to previous best
    ## Run 363 stress 9.404659e-05 
    ## ... Procrustes: rmse 0.0001198185  max resid 0.0002128031 
    ## ... Similar to previous best
    ## Run 364 stress 9.121951e-05 
    ## ... Procrustes: rmse 0.0001379559  max resid 0.0002506189 
    ## ... Similar to previous best
    ## Run 365 stress 9.612578e-05 
    ## ... Procrustes: rmse 0.0001225352  max resid 0.00027029 
    ## ... Similar to previous best
    ## Run 366 stress 9.255427e-05 
    ## ... Procrustes: rmse 0.000143083  max resid 0.0002893778 
    ## ... Similar to previous best
    ## Run 367 stress 9.389069e-05 
    ## ... Procrustes: rmse 0.0001972369  max resid 0.0003591058 
    ## ... Similar to previous best
    ## Run 368 stress 9.674095e-05 
    ## ... Procrustes: rmse 0.0002069247  max resid 0.0003806078 
    ## ... Similar to previous best
    ## Run 369 stress 9.983901e-05 
    ## ... Procrustes: rmse 0.000188945  max resid 0.0002618779 
    ## ... Similar to previous best
    ## Run 370 stress 9.688162e-05 
    ## ... Procrustes: rmse 0.0002006107  max resid 0.0003673343 
    ## ... Similar to previous best
    ## Run 371 stress 8.766252e-05 
    ## ... Procrustes: rmse 0.0001202476  max resid 0.0001924031 
    ## ... Similar to previous best
    ## Run 372 stress 8.80607e-05 
    ## ... Procrustes: rmse 0.0001191403  max resid 0.0001726873 
    ## ... Similar to previous best
    ## Run 373 stress 9.452298e-05 
    ## ... Procrustes: rmse 0.0001653282  max resid 0.0002487291 
    ## ... Similar to previous best
    ## Run 374 stress 9.435937e-05 
    ## ... Procrustes: rmse 0.0001219014  max resid 0.0002144292 
    ## ... Similar to previous best
    ## Run 375 stress 8.664335e-05 
    ## ... Procrustes: rmse 0.0001311012  max resid 0.0002399178 
    ## ... Similar to previous best
    ## Run 376 stress 9.680978e-05 
    ## ... Procrustes: rmse 0.0001226374  max resid 0.0002213079 
    ## ... Similar to previous best
    ## Run 377 stress 9.42461e-05 
    ## ... Procrustes: rmse 0.0001983915  max resid 0.0003642787 
    ## ... Similar to previous best
    ## Run 378 stress 9.843762e-05 
    ## ... Procrustes: rmse 0.0002290824  max resid 0.0003899873 
    ## ... Similar to previous best
    ## Run 379 stress 9.354567e-05 
    ## ... Procrustes: rmse 0.0001907785  max resid 0.0003542538 
    ## ... Similar to previous best
    ## Run 380 stress 9.56062e-05 
    ## ... Procrustes: rmse 0.0001615143  max resid 0.0003003933 
    ## ... Similar to previous best
    ## Run 381 stress 9.80047e-05 
    ## ... Procrustes: rmse 0.0002071231  max resid 0.000372243 
    ## ... Similar to previous best
    ## Run 382 stress 9.14165e-05 
    ## ... Procrustes: rmse 0.0001926896  max resid 0.000355462 
    ## ... Similar to previous best
    ## Run 383 stress 9.537722e-05 
    ## ... Procrustes: rmse 0.0002038005  max resid 0.0003725302 
    ## ... Similar to previous best
    ## Run 384 stress 9.653858e-05 
    ## ... Procrustes: rmse 0.0002174915  max resid 0.0003627396 
    ## ... Similar to previous best
    ## Run 385 stress 0.3051802 
    ## Run 386 stress 9.054154e-05 
    ## ... Procrustes: rmse 0.0002083813  max resid 0.0003534977 
    ## ... Similar to previous best
    ## Run 387 stress 9.878583e-05 
    ## ... Procrustes: rmse 0.0001262027  max resid 0.0002189272 
    ## ... Similar to previous best
    ## Run 388 stress 9.167208e-05 
    ## ... Procrustes: rmse 0.0001345634  max resid 0.0001966887 
    ## ... Similar to previous best
    ## Run 389 stress 9.597332e-05 
    ## ... Procrustes: rmse 0.0001634331  max resid 0.0003201415 
    ## ... Similar to previous best
    ## Run 390 stress 9.566732e-05 
    ## ... Procrustes: rmse 0.0002152981  max resid 0.0003679667 
    ## ... Similar to previous best
    ## Run 391 stress 9.62609e-05 
    ## ... Procrustes: rmse 0.0001445582  max resid 0.0002371229 
    ## ... Similar to previous best
    ## Run 392 stress 9.931891e-05 
    ## ... Procrustes: rmse 0.0002060099  max resid 0.000381115 
    ## ... Similar to previous best
    ## Run 393 stress 9.949541e-05 
    ## ... Procrustes: rmse 0.0002148287  max resid 0.0003919772 
    ## ... Similar to previous best
    ## Run 394 stress 9.102044e-05 
    ## ... Procrustes: rmse 0.0001447699  max resid 0.0002958245 
    ## ... Similar to previous best
    ## Run 395 stress 9.966301e-05 
    ## ... Procrustes: rmse 0.0002111831  max resid 0.0003898733 
    ## ... Similar to previous best
    ## Run 396 stress 8.450545e-05 
    ## ... Procrustes: rmse 0.0001727288  max resid 0.0003240968 
    ## ... Similar to previous best
    ## Run 397 stress 9.497791e-05 
    ## ... Procrustes: rmse 0.0002159804  max resid 0.0003618468 
    ## ... Similar to previous best
    ## Run 398 stress 8.416118e-05 
    ## ... Procrustes: rmse 0.0001157806  max resid 0.0001988819 
    ## ... Similar to previous best
    ## Run 399 stress 9.847842e-05 
    ## ... Procrustes: rmse 0.0002077896  max resid 0.0003772019 
    ## ... Similar to previous best
    ## Run 400 stress 9.63255e-05 
    ## ... Procrustes: rmse 0.000222189  max resid 0.0003757665 
    ## ... Similar to previous best
    ## Run 401 stress 9.405246e-05 
    ## ... Procrustes: rmse 0.0001512019  max resid 0.0003099071 
    ## ... Similar to previous best
    ## Run 402 stress 9.869135e-05 
    ## ... Procrustes: rmse 0.000211143  max resid 0.0003882865 
    ## ... Similar to previous best
    ## Run 403 stress 9.524987e-05 
    ## ... Procrustes: rmse 0.0002213117  max resid 0.0003730784 
    ## ... Similar to previous best
    ## Run 404 stress 9.894442e-05 
    ## ... Procrustes: rmse 0.0001518401  max resid 0.0003055943 
    ## ... Similar to previous best
    ## Run 405 stress 9.981513e-05 
    ## ... Procrustes: rmse 0.0001306019  max resid 0.0002289729 
    ## ... Similar to previous best
    ## Run 406 stress 9.278973e-05 
    ## ... Procrustes: rmse 0.0001896337  max resid 0.0003519214 
    ## ... Similar to previous best
    ## Run 407 stress 9.395886e-05 
    ## ... Procrustes: rmse 0.0001971455  max resid 0.0003652035 
    ## ... Similar to previous best
    ## Run 408 stress 9.78412e-05 
    ## ... Procrustes: rmse 0.000187438  max resid 0.0002718815 
    ## ... Similar to previous best
    ## Run 409 stress 9.108217e-05 
    ## ... Procrustes: rmse 0.0001584907  max resid 0.0002349102 
    ## ... Similar to previous best
    ## Run 410 stress 9.431904e-05 
    ## ... Procrustes: rmse 0.0001194849  max resid 0.0002133704 
    ## ... Similar to previous best
    ## Run 411 stress 9.835304e-05 
    ## ... Procrustes: rmse 9.245601e-05  max resid 0.0001539168 
    ## ... Similar to previous best
    ## Run 412 stress 9.859609e-05 
    ## ... Procrustes: rmse 0.0002113088  max resid 0.000389006 
    ## ... Similar to previous best
    ## Run 413 stress 9.734038e-05 
    ## ... Procrustes: rmse 0.0001973286  max resid 0.0003614579 
    ## ... Similar to previous best
    ## Run 414 stress 9.727136e-05 
    ## ... Procrustes: rmse 0.0001574949  max resid 0.0003358725 
    ## ... Similar to previous best
    ## Run 415 stress 9.875408e-05 
    ## ... Procrustes: rmse 0.0001947102  max resid 0.000363142 
    ## ... Similar to previous best
    ## Run 416 stress 9.208063e-05 
    ## ... Procrustes: rmse 0.0001157189  max resid 0.0002023872 
    ## ... Similar to previous best
    ## Run 417 stress 9.988293e-05 
    ## ... Procrustes: rmse 0.000213902  max resid 0.0003896087 
    ## ... Similar to previous best
    ## Run 418 stress 9.906742e-05 
    ## ... Procrustes: rmse 0.0001252541  max resid 0.0001955544 
    ## ... Similar to previous best
    ## Run 419 stress 9.609297e-05 
    ## ... Procrustes: rmse 0.0002053755  max resid 0.0003784795 
    ## ... Similar to previous best
    ## Run 420 stress 9.958452e-05 
    ## ... Procrustes: rmse 0.0001204093  max resid 0.0002089189 
    ## ... Similar to previous best
    ## Run 421 stress 9.969021e-05 
    ## ... Procrustes: rmse 0.0001486995  max resid 0.0002598806 
    ## ... Similar to previous best
    ## Run 422 stress 0.3083098 
    ## Run 423 stress 8.727603e-05 
    ## ... Procrustes: rmse 0.0001096655  max resid 0.0001889066 
    ## ... Similar to previous best
    ## Run 424 stress 9.698973e-05 
    ## ... Procrustes: rmse 0.0002014239  max resid 0.0003338465 
    ## ... Similar to previous best
    ## Run 425 stress 9.934804e-05 
    ## ... Procrustes: rmse 0.0001843012  max resid 0.0002560216 
    ## ... Similar to previous best
    ## Run 426 stress 9.943592e-05 
    ## ... Procrustes: rmse 0.0001320789  max resid 0.0002557658 
    ## ... Similar to previous best
    ## Run 427 stress 9.158482e-05 
    ## ... Procrustes: rmse 0.000114396  max resid 0.0001999243 
    ## ... Similar to previous best
    ## Run 428 stress 9.667613e-05 
    ## ... Procrustes: rmse 0.0002046551  max resid 0.0003732897 
    ## ... Similar to previous best
    ## Run 429 stress 9.359514e-05 
    ## ... Procrustes: rmse 0.0001708295  max resid 0.0002540231 
    ## ... Similar to previous best
    ## Run 430 stress 9.486055e-05 
    ## ... Procrustes: rmse 0.0001467505  max resid 0.0003105376 
    ## ... Similar to previous best
    ## Run 431 stress 9.868903e-05 
    ## ... Procrustes: rmse 0.0001283927  max resid 0.0002293656 
    ## ... Similar to previous best
    ## Run 432 stress 9.480474e-05 
    ## ... Procrustes: rmse 0.0001983323  max resid 0.0003604361 
    ## ... Similar to previous best
    ## Run 433 stress 9.943983e-05 
    ## ... Procrustes: rmse 0.0001557673  max resid 0.0002525801 
    ## ... Similar to previous best
    ## Run 434 stress 9.126744e-05 
    ## ... Procrustes: rmse 0.000134025  max resid 0.0002892169 
    ## ... Similar to previous best
    ## Run 435 stress 8.95589e-05 
    ## ... Procrustes: rmse 0.00018584  max resid 0.0003473383 
    ## ... Similar to previous best
    ## Run 436 stress 9.510438e-05 
    ## ... Procrustes: rmse 0.0001211731  max resid 0.0002163624 
    ## ... Similar to previous best
    ## Run 437 stress 9.934594e-05 
    ## ... Procrustes: rmse 0.0001795043  max resid 0.000328202 
    ## ... Similar to previous best
    ## Run 438 stress 9.253871e-05 
    ## ... Procrustes: rmse 0.0001525795  max resid 0.0002381902 
    ## ... Similar to previous best
    ## Run 439 stress 8.973489e-05 
    ## ... Procrustes: rmse 0.0001928741  max resid 0.000331633 
    ## ... Similar to previous best
    ## Run 440 stress 9.209995e-05 
    ## ... Procrustes: rmse 0.0001144387  max resid 0.0002025402 
    ## ... Similar to previous best
    ## Run 441 stress 8.966447e-05 
    ## ... Procrustes: rmse 0.0001301961  max resid 0.0002222686 
    ## ... Similar to previous best
    ## Run 442 stress 9.992974e-05 
    ## ... Procrustes: rmse 0.0002034783  max resid 0.0003335234 
    ## ... Similar to previous best
    ## Run 443 stress 9.348028e-05 
    ## ... Procrustes: rmse 0.0002132939  max resid 0.0003603277 
    ## ... Similar to previous best
    ## Run 444 stress 9.952587e-05 
    ## ... Procrustes: rmse 0.0001415553  max resid 0.0002354704 
    ## ... Similar to previous best
    ## Run 445 stress 9.323348e-05 
    ## ... Procrustes: rmse 0.0001232576  max resid 0.0002073796 
    ## ... Similar to previous best
    ## Run 446 stress 9.927801e-05 
    ## ... Procrustes: rmse 0.0001763769  max resid 0.0003096734 
    ## ... Similar to previous best
    ## Run 447 stress 9.337191e-05 
    ## ... Procrustes: rmse 0.0001661683  max resid 0.0003243218 
    ## ... Similar to previous best
    ## Run 448 stress 0.3083098 
    ## Run 449 stress 9.256331e-05 
    ## ... Procrustes: rmse 0.0001894389  max resid 0.0003537353 
    ## ... Similar to previous best
    ## Run 450 stress 9.794962e-05 
    ## ... Procrustes: rmse 0.0001451857  max resid 0.0002181469 
    ## ... Similar to previous best
    ## Run 451 stress 9.094373e-05 
    ## ... Procrustes: rmse 0.0001131481  max resid 0.0002018806 
    ## ... Similar to previous best
    ## Run 452 stress 9.158301e-05 
    ## ... Procrustes: rmse 0.000147387  max resid 0.0003107413 
    ## ... Similar to previous best
    ## Run 453 stress 9.626669e-05 
    ## ... Procrustes: rmse 0.000161276  max resid 0.0002632012 
    ## ... Similar to previous best
    ## Run 454 stress 9.739905e-05 
    ## ... Procrustes: rmse 0.0002246679  max resid 0.0003778521 
    ## ... Similar to previous best
    ## Run 455 stress 9.732552e-05 
    ## ... Procrustes: rmse 0.0001489086  max resid 0.0002534399 
    ## ... Similar to previous best
    ## Run 456 stress 7.389172e-05 
    ## ... Procrustes: rmse 0.0001063904  max resid 0.0001480271 
    ## ... Similar to previous best
    ## Run 457 stress 8.580265e-05 
    ## ... Procrustes: rmse 0.0001213882  max resid 0.0002081523 
    ## ... Similar to previous best
    ## Run 458 stress 9.64024e-05 
    ## ... Procrustes: rmse 0.0001623234  max resid 0.0002687295 
    ## ... Similar to previous best
    ## Run 459 stress 9.568781e-05 
    ## ... Procrustes: rmse 0.0001571973  max resid 0.00023509 
    ## ... Similar to previous best
    ## Run 460 stress 9.636651e-05 
    ## ... Procrustes: rmse 0.0002215843  max resid 0.0003729055 
    ## ... Similar to previous best
    ## Run 461 stress 9.432143e-05 
    ## ... Procrustes: rmse 0.0001891267  max resid 0.0003471503 
    ## ... Similar to previous best
    ## Run 462 stress 7.846388e-05 
    ## ... Procrustes: rmse 8.853721e-05  max resid 0.0001578583 
    ## ... Similar to previous best
    ## Run 463 stress 9.330841e-05 
    ## ... Procrustes: rmse 0.0001385953  max resid 0.0002287969 
    ## ... Similar to previous best
    ## Run 464 stress 0.307194 
    ## Run 465 stress 8.960752e-05 
    ## ... Procrustes: rmse 0.0001565431  max resid 0.0003104411 
    ## ... Similar to previous best
    ## Run 466 stress 9.628099e-05 
    ## ... Procrustes: rmse 0.0001245366  max resid 0.0002225737 
    ## ... Similar to previous best
    ## Run 467 stress 9.347689e-05 
    ## ... Procrustes: rmse 0.0001941454  max resid 0.000359899 
    ## ... Similar to previous best
    ## Run 468 stress 9.566321e-05 
    ## ... Procrustes: rmse 0.0001976239  max resid 0.0003597709 
    ## ... Similar to previous best
    ## Run 469 stress 9.622481e-05 
    ## ... Procrustes: rmse 0.000222427  max resid 0.0003794738 
    ## ... Similar to previous best
    ## Run 470 stress 9.065273e-05 
    ## ... Procrustes: rmse 0.0001759307  max resid 0.0003103115 
    ## ... Similar to previous best
    ## Run 471 stress 9.955378e-05 
    ## ... Procrustes: rmse 0.0002316611  max resid 0.0003948054 
    ## ... Similar to previous best
    ## Run 472 stress 9.071179e-05 
    ## ... Procrustes: rmse 0.0001481255  max resid 0.0003142143 
    ## ... Similar to previous best
    ## Run 473 stress 9.17135e-05 
    ## ... Procrustes: rmse 0.0001149341  max resid 0.0002050886 
    ## ... Similar to previous best
    ## Run 474 stress 9.781315e-05 
    ## ... Procrustes: rmse 0.000222237  max resid 0.0003783674 
    ## ... Similar to previous best
    ## Run 475 stress 8.877825e-05 
    ## ... Procrustes: rmse 0.0001323225  max resid 0.0002282367 
    ## ... Similar to previous best
    ## Run 476 stress 9.966463e-05 
    ## ... Procrustes: rmse 0.000214029  max resid 0.0003947158 
    ## ... Similar to previous best
    ## Run 477 stress 9.395248e-05 
    ## ... Procrustes: rmse 0.0002003036  max resid 0.0003661735 
    ## ... Similar to previous best
    ## Run 478 stress 9.996923e-05 
    ## ... Procrustes: rmse 0.0002027004  max resid 0.0003694592 
    ## ... Similar to previous best
    ## Run 479 stress 8.05881e-05 
    ## ... Procrustes: rmse 8.59819e-05  max resid 0.0001356954 
    ## ... Similar to previous best
    ## Run 480 stress 9.325777e-05 
    ## ... Procrustes: rmse 0.0001474726  max resid 0.0002179935 
    ## ... Similar to previous best
    ## Run 481 stress 9.259097e-05 
    ## ... Procrustes: rmse 0.0001652422  max resid 0.0002298463 
    ## ... Similar to previous best
    ## Run 482 stress 9.644621e-05 
    ## ... Procrustes: rmse 0.0001769376  max resid 0.0002372785 
    ## ... Similar to previous best
    ## Run 483 stress 5.55999e-05 
    ## ... Procrustes: rmse 8.898691e-05  max resid 0.0001582767 
    ## ... Similar to previous best
    ## Run 484 stress 9.252866e-05 
    ## ... Procrustes: rmse 0.0002100954  max resid 0.0003577539 
    ## ... Similar to previous best
    ## Run 485 stress 9.062302e-05 
    ## ... Procrustes: rmse 0.0001404613  max resid 0.0002189865 
    ## ... Similar to previous best
    ## Run 486 stress 8.926554e-05 
    ## ... Procrustes: rmse 0.0001080082  max resid 0.0001918372 
    ## ... Similar to previous best
    ## Run 487 stress 9.064706e-05 
    ## ... Procrustes: rmse 0.0001078747  max resid 0.0001855943 
    ## ... Similar to previous best
    ## Run 488 stress 9.673782e-05 
    ## ... Procrustes: rmse 0.000164374  max resid 0.0002641776 
    ## ... Similar to previous best
    ## Run 489 stress 9.682716e-05 
    ## ... Procrustes: rmse 0.0002060312  max resid 0.0003806672 
    ## ... Similar to previous best
    ## Run 490 stress 9.982684e-05 
    ## ... Procrustes: rmse 0.0002152218  max resid 0.0003942349 
    ## ... Similar to previous best
    ## Run 491 stress 9.094275e-05 
    ## ... Procrustes: rmse 0.0001302933  max resid 0.0002204631 
    ## ... Similar to previous best
    ## Run 492 stress 8.726208e-05 
    ## ... Procrustes: rmse 0.0001225014  max resid 0.0002047898 
    ## ... Similar to previous best
    ## Run 493 stress 0.3040432 
    ## Run 494 stress 9.214432e-05 
    ## ... Procrustes: rmse 0.0001776986  max resid 0.0002504941 
    ## ... Similar to previous best
    ## Run 495 stress 9.69016e-05 
    ## ... Procrustes: rmse 0.0002067787  max resid 0.0003770324 
    ## ... Similar to previous best
    ## Run 496 stress 9.364455e-05 
    ## ... Procrustes: rmse 0.0002122452  max resid 0.0003607509 
    ## ... Similar to previous best
    ## Run 497 stress 9.185512e-05 
    ## ... Procrustes: rmse 0.0001905143  max resid 0.0003486785 
    ## ... Similar to previous best
    ## Run 498 stress 9.962767e-05 
    ## ... Procrustes: rmse 0.0001319606  max resid 0.0002349784 
    ## ... Similar to previous best
    ## Run 499 stress 9.397977e-05 
    ## ... Procrustes: rmse 0.0001185566  max resid 0.0002122887 
    ## ... Similar to previous best
    ## Run 500 stress 9.408372e-05 
    ## ... Procrustes: rmse 0.0001991288  max resid 0.0003672403 
    ## ... Similar to previous best
    ## *** Best solution repeated 401 times

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
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31832 0.699  0.022 *
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
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31832 0.699  0.022 *
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
    ## temperature_median -0.75818 -0.65204 0.0724  0.699  
    ## salinity_median     0.94798  0.31832 0.6990  0.023 *
    ## oxygen_median       0.57185 -0.82036 0.6836  0.110  
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
    ## temperature_median -0.75818 -0.65204 0.0724  1.000  
    ## salinity_median     0.94798  0.31832 0.6990  0.069 .
    ## oxygen_median       0.57185 -0.82036 0.6836  0.330  
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
    ## temperature_median -0.92689  0.37534 0.0954  0.748  
    ## salinity_median     0.90474  0.42596 0.6840  0.038 *
    ## oxygen_median       0.44478 -0.89564 0.5836  0.218  
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
    ## temperature_median -0.92689  0.37534 0.0954  1.000
    ## salinity_median     0.90474  0.42596 0.6840  0.114
    ## oxygen_median       0.44478 -0.89564 0.5836  0.654
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
    ## temperature_median -0.52519 -0.85098 0.0223  0.902  
    ## salinity_median    -0.79785 -0.60286 0.7269  0.016 *
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.015 *
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
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.045 *
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
    ## temperature_median  0.0001722 -1.0000000 0.1579  0.264
    ## salinity_median    -0.0063343  0.9999800 0.4836  0.873
    ## oxygen_median      -0.0014073  1.0000000 0.7574  0.228
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
    ## temperature_median  0.0001722 -1.0000000 0.1579  0.792
    ## salinity_median    -0.0063343  0.9999800 0.4836  1.000
    ## oxygen_median      -0.0014073  1.0000000 0.7574  0.684
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
    ## temperature_median -0.00007504  1.00000000 0.0973  0.765  
    ## salinity_median     0.00161315  1.00000000 0.3817  0.271  
    ## oxygen_median       0.00048263 -1.00000000 0.6352  0.077 .
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
    ## temperature_median -0.00007504  1.00000000 0.0973  1.000
    ## salinity_median     0.00161315  1.00000000 0.3817  0.813
    ## oxygen_median       0.00048263 -1.00000000 0.6352  0.231
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
    ## temperature_median       -0.39655 -0.91801 0.0381  0.915  
    ## salinity_median          -0.96163  0.27434 0.5830  0.109  
    ## oxygen_median             0.83739  0.54660 0.7450  0.026 *
    ## distance_to_ocean_mean_m -0.09989  0.99500 0.0319  0.936  
    ## max_depth                -0.87439  0.48522 0.1022  0.773  
    ## logArea                  -0.40921  0.91244 0.2149  0.576  
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
    ## temperature_median       -0.39655 -0.91801 0.0381  1.000
    ## salinity_median          -0.96163  0.27434 0.5830  0.654
    ## oxygen_median             0.83739  0.54660 0.7450  0.156
    ## distance_to_ocean_mean_m -0.09989  0.99500 0.0319  1.000
    ## max_depth                -0.87439  0.48522 0.1022  1.000
    ## logArea                  -0.40921  0.91244 0.2149  1.000
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
# For figure
(PD_beta_geo_ef <- envfit(PD_beta_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.82791  0.56086 0.5558  0.084 .
    ## max_depth               -0.20869 -0.97798 0.0859  0.902  
    ## logArea                  0.25301 -0.96746 0.2012  0.306  
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
    ## distance_to_ocean_min_m -0.82791  0.56086 0.5558  0.252
    ## max_depth               -0.20869 -0.97798 0.0859  1.000
    ## logArea                  0.25301 -0.96746 0.2012  0.918
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by stratification
(PD_beta_geo_A_ef <- envfit(PD_beta_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.82791  0.56086 0.5558  0.068 .
    ## max_depth               -0.20869 -0.97798 0.0859  0.900  
    ## logArea                  0.25301 -0.96746 0.2012  0.280  
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
    ## distance_to_ocean_min_m -0.82791  0.56086 0.5558  0.204
    ## max_depth               -0.20869 -0.97798 0.0859  1.000
    ## logArea                  0.25301 -0.96746 0.2012  0.840
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
(PD_beta_geo_MS_ef <- envfit(PD_beta_geo_MS_NMDS, env[mixed_stratified_lakes,geography], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.76307  0.64632 0.4671  0.108
    ## max_depth               -0.34716 -0.93780 0.0777  0.966
    ## logArea                  0.55445  0.83222 0.0074  0.974
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
    ## distance_to_ocean_min_m -0.76307  0.64632 0.4671  0.324
    ## max_depth               -0.34716 -0.93780 0.0777  1.000
    ## logArea                  0.55445  0.83222 0.0074  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
(PD_beta_geo_OM_ef <- envfit(PD_beta_geo_OM_NMDS, env[ocean_mixed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m  0.45615 -0.88990 0.4273  0.162   
    ## max_depth               -0.94118 -0.33789 0.6682  0.005 **
    ## logArea                 -0.86397  0.50354 0.2339  0.453   
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
    ## distance_to_ocean_min_m  0.45615 -0.88990 0.4273  0.486  
    ## max_depth               -0.94118 -0.33789 0.6682  0.015 *
    ## logArea                 -0.86397  0.50354 0.2339  1.000  
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
    ## distance_to_ocean_min_m  0.99474  0.10244 0.6338  0.233
    ## max_depth                0.73022  0.68321 0.0374  0.973
    ## logArea                 -0.39139  0.92023 0.3078  0.322
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
    ## distance_to_ocean_min_m  0.99474  0.10244 0.6338  0.699
    ## max_depth                0.73022  0.68321 0.0374  1.000
    ## logArea                 -0.39139  0.92023 0.3078  0.966
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed lakes
(PD_beta_geo_M_ef <- envfit(PD_beta_geo_M_NMDS, env[mixed_lakes,geography], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1       NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m -0.00013720  1.00000000 0.7961  0.030 * 
    ## max_depth                0.00054677 -1.00000000 0.7721  0.043 * 
    ## logArea                  0.00031450  1.00000000 0.8439  0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
(PD_beta_geo_M_efp <- p.adjust.envfit(PD_beta_geo_M_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1       NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.00013720  1.00000000 0.7961  0.090 .
    ## max_depth                0.00054677 -1.00000000 0.7721  0.129  
    ## logArea                  0.00031450  1.00000000 0.8439  0.018 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

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
    ##       Significance: 0.216 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.219 0.252 0.279 0.301 
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
    ## 0.438 0.463 0.483 0.504 
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
    ##       Significance: 0.244 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.522 0.559 0.581 0.611 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_mant_pv <- rbind(PD_beta_env_mant_t$signif, PD_beta_env_mant_s$signif, PD_beta_env_mant_o$signif)
PD_beta_env_mant_pv <- PD_beta_env_mant_pv[,1]
(PD_beta_env_mant_pv <- p.adjust(PD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.648 0.015 0.732

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
    ##       Significance: 0.236 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.221 0.258 0.287 0.318 
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
    ##       Significance: 0.013 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.381 0.420 0.446 0.468 
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
    ##       Significance: 0.294 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.405 0.445 0.476 0.520 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_MS_mant_pv <- rbind(PD_beta_env_MS_mant_t$signif, PD_beta_env_MS_mant_s$signif, PD_beta_env_MS_mant_o$signif)
PD_beta_env_MS_mant_pv <- PD_beta_env_MS_mant_pv[,1]
(PD_beta_env_MS_mant_pv <- p.adjust(PD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.708 0.039 0.882

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
    ##       Significance: 0.532 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.137 0.190 0.277 0.376 
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
    ## 0.305 0.357 0.424 0.483 
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
    ##       Significance: 0.013 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.284 0.379 0.464 0.578 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_OM_mant_pv <- rbind(PD_beta_env_OM_mant_t$signif, PD_beta_env_OM_mant_s$signif, PD_beta_env_OM_mant_o$signif)
PD_beta_env_OM_mant_pv <- PD_beta_env_OM_mant_pv[,1]
(PD_beta_env_OM_mant_pv <- p.adjust(PD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.048 0.039

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
    ##       Significance: 0.46 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0314 0.0651 0.0948 0.1290 
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
    ##       Significance: 0.021 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.387 0.428 0.462 0.493 
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
    ##       Significance: 0.227 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.690 0.718 0.737 0.751 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_SO_mant_pv <- rbind(PD_beta_env_SO_mant_t$signif, PD_beta_env_SO_mant_s$signif, PD_beta_env_SO_mant_o$signif)
PD_beta_env_SO_mant_pv <- PD_beta_env_SO_mant_pv[,1]
(PD_beta_env_SO_mant_pv <- p.adjust(PD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.063 0.681

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
    ##       Significance: 0.781 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.252 0.330 0.394 0.600 
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
    ##       Significance: 0.234 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.201 0.259 0.321 0.390 
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
    ##       Significance: 0.074 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.228 0.335 0.417 0.556 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_M_mant_pv <- rbind(PD_beta_env_M_mant_t$signif, PD_beta_env_M_mant_s$signif, PD_beta_env_M_mant_o$signif)
PD_beta_env_M_mant_pv <- PD_beta_env_M_mant_pv[,1]
(PD_beta_env_M_mant_pv <- p.adjust(PD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.702 0.222

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
    ##       Significance: 0.301 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.348 0.440 0.521 0.607 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
geo_dist_dmean <- dist(scaled_env[surveyed_sites,c(9)], method = "euclidean")
(PD_beta_geo_mant_dmean <- mantel(PD_beta_geo_dist$Btotal, geo_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_dist$Btotal, ydis = geo_dist_dmean,      method = "spearman", permutations = 999, strata = env[surveyed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3483 
    ##       Significance: 0.596 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.459 0.486 0.507 0.526 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_md <- dist(scaled_env[surveyed_sites,c(14)], method = "euclidean")
(PD_beta_geo_mant_md <- mantel(PD_beta_geo_dist$Btotal, geo_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_dist$Btotal, ydis = geo_dist_md, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.0242 
    ##       Significance: 0.782 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.167 0.194 0.217 0.242 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_la <- dist(scaled_env[surveyed_sites,c(15)], method = "euclidean")
(PD_beta_geo_mant_la <- mantel(PD_beta_geo_dist$Btotal, geo_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_dist$Btotal, ydis = geo_dist_la, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.01132 
    ##       Significance: 0.711 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.114 0.128 0.148 0.165 
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
geo_MS_dist_dmean <- dist(scaled_env[mixed_stratified_lakes,c(9)], method = "euclidean")
(PD_beta_geo_MS_mant_dmean <- mantel(PD_beta_geo_MS_dist$Btotal, geo_MS_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_dmean,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1391 
    ##       Significance: 0.639 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.292 0.335 0.365 0.409 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_md <- dist(scaled_env[mixed_stratified_lakes,c(14)], method = "euclidean")
(PD_beta_geo_MS_mant_md <- mantel(PD_beta_geo_MS_dist$Btotal, geo_MS_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_md,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.05502 
    ##       Significance: 0.901 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.149 0.190 0.225 0.258 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_la <- dist(scaled_env[mixed_stratified_lakes,c(15)], method = "euclidean")
(PD_beta_geo_MS_mant_la <- mantel(PD_beta_geo_MS_dist$Btotal, geo_MS_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_la,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.06883 
    ##       Significance: 0.897 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.129 0.163 0.182 0.227 
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
geo_OM_dist_dmean <- dist(scaled_env[ocean_mixed_sites,c(9)], method = "euclidean")
(PD_beta_geo_OM_mant_dmean <- mantel(PD_beta_geo_OM_dist$Btotal, geo_OM_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_dmean,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.06927 
    ##       Significance: 0.38 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.169 0.197 0.219 0.253 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_md <- dist(scaled_env[ocean_mixed_sites,c(14)], method = "euclidean")
(PD_beta_geo_OM_mant_md <- mantel(PD_beta_geo_OM_dist$Btotal, geo_OM_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_md,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3135 
    ##       Significance: 0.016 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.153 0.218 0.264 0.333 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_la <- dist(scaled_env[ocean_mixed_sites,c(15)], method = "euclidean")
(PD_beta_geo_OM_mant_la <- mantel(PD_beta_geo_OM_dist$Btotal, geo_OM_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_la,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1943 
    ##       Significance: 0.072 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.173 0.205 0.234 0.274 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_OM_mant_pv <- rbind( PD_beta_geo_OM_mant_dmean$signif, PD_beta_geo_OM_mant_md$signif, PD_beta_geo_OM_mant_la$signif)
PD_beta_geo_OM_mant_pv <- PD_beta_geo_OM_mant_pv[,1]
(PD_beta_geo_OM_mant_pv <- p.adjust(PD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.048 0.216

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
    ##       Significance: 0.519 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.619 0.649 0.678 0.694 
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
    ##       Significance: 0.796 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.113 0.134 0.169 0.198 
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
    ##       Significance: 0.22 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.376 0.404 0.424 0.440 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_SO_mant_pv <- rbind( PD_beta_geo_SO_mant_dmean$signif, PD_beta_geo_SO_mant_md$signif, PD_beta_geo_SO_mant_la$signif)
PD_beta_geo_SO_mant_pv <- PD_beta_geo_SO_mant_pv[,1]
(PD_beta_geo_SO_mant_pv <- p.adjust(PD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.00 1.00 0.66

``` r
# Mixed lakes
geo_M_dist_dmean <- dist(scaled_env[mixed_lakes,c(9)], method = "euclidean")
(PD_beta_geo_M_mant_dmean <- mantel(PD_beta_geo_M_dist$Btotal, geo_M_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_dmean,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1161 
    ##       Significance: 0.71 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.222 0.310 0.390 0.692 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_M_dist_md <- dist(scaled_env[mixed_lakes,c(14)], method = "euclidean")
(PD_beta_geo_M_mant_md <- mantel(PD_beta_geo_M_dist$Btotal, geo_M_dist_md, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_md,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5985 
    ##       Significance: 0.011 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.293 0.407 0.506 0.583 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_M_dist_la <- dist(scaled_env[mixed_lakes,c(15)], method = "euclidean")
(PD_beta_geo_M_mant_la <- mantel(PD_beta_geo_M_dist$Btotal, geo_M_dist_la, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = PD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_la,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4007 
    ##       Significance: 0.017 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.204 0.290 0.353 0.427 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_M_mant_pv <- rbind( PD_beta_geo_M_mant_dmean$signif, PD_beta_geo_M_mant_md$signif, PD_beta_geo_M_mant_la$signif)
PD_beta_geo_M_mant_pv <- PD_beta_geo_M_mant_pv[,1]
(PD_beta_geo_M_mant_pv <- p.adjust(PD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.033 0.051

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
