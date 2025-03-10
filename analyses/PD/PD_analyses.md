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

### PD alpha mntd, mpd, and pd calculations of phylogeny and files to be saved

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

- Read in Phylogenetic Diversity files and combine with env data frame

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

### PD alpha outliers for each site type

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

### Identify VIF of environmental and geographic variables

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

### PD alpha & site type

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

### PD alpha with env and geo linear models & ANOVAs

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

### PD alpha with env and geo correlated variables

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
```

### PD alpha z-scores & dispersion

``` r
PD_alpha_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.z, x = Dispersion, color = Stratification, fill = Stratification)) +
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = stree_sespd_env[,1], size = 5, point.padding = 3) +
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
  labs(x="PD alpha Dispersion", y="PD alpha z-score", colour = "Site type", fill = "Site type", tag = "a")
(PD_alpha_plot <- PD_alpha_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/z-scores%20&%20dispersion-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_alpha_plot.jpg", plot = PD_alpha_plot, width = 4.88, height = 6, units = "in")
```

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
    ## 3 Stratified vs Reference  1 0.7729479 7.110647 0.5039207   0.102      0.612
    ## 4          Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.093      0.558
    ## 5      Mixed vs Reference  1 0.6592270 3.967203 0.3617334   0.102      0.612
    ## 6      Ocean vs Reference  1 0.6319662 3.354877 0.4015472   0.147      0.882
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
    ## 3      Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.092      0.276

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
    ## 2 Stratified vs Ocean  1 1.2083230 9.439158 0.4618174   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.3347638 2.034034 0.1560556   0.023      0.069

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
    ## 1 Stratified vs Mixed  1 1.3499033 10.500850 0.4666868   0.002      0.006   *
    ## 2 Stratified vs Ocean  1 1.2146330  9.192717 0.4789690   0.004      0.012   .
    ## 3      Mixed vs Ocean  1 0.2573592  1.467099 0.1089395   0.093      0.279

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

### PD beta phylogenetic dispersions

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
    ## Run 1 stress 9.906576e-05 
    ## ... Procrustes: rmse 0.0001133848  max resid 0.0003491934 
    ## ... Similar to previous best
    ## Run 2 stress 4.738901e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002321375  max resid 0.000383672 
    ## ... Similar to previous best
    ## Run 3 stress 9.548172e-05 
    ## ... Procrustes: rmse 0.0002019027  max resid 0.0003318656 
    ## ... Similar to previous best
    ## Run 4 stress 9.666749e-05 
    ## ... Procrustes: rmse 0.0002066786  max resid 0.000340836 
    ## ... Similar to previous best
    ## Run 5 stress 9.222035e-05 
    ## ... Procrustes: rmse 0.0002294947  max resid 0.0003845618 
    ## ... Similar to previous best
    ## Run 6 stress 9.91004e-05 
    ## ... Procrustes: rmse 0.0002451891  max resid 0.0004019308 
    ## ... Similar to previous best
    ## Run 7 stress 9.942373e-05 
    ## ... Procrustes: rmse 0.000195071  max resid 0.0003245132 
    ## ... Similar to previous best
    ## Run 8 stress 8.116102e-05 
    ## ... Procrustes: rmse 5.441835e-05  max resid 0.0001139688 
    ## ... Similar to previous best
    ## Run 9 stress 9.82123e-05 
    ## ... Procrustes: rmse 0.0002154673  max resid 0.0003539293 
    ## ... Similar to previous best
    ## Run 10 stress 5.068774e-05 
    ## ... Procrustes: rmse 3.192693e-05  max resid 6.396154e-05 
    ## ... Similar to previous best
    ## Run 11 stress 9.277893e-05 
    ## ... Procrustes: rmse 0.0001069345  max resid 0.0001913484 
    ## ... Similar to previous best
    ## Run 12 stress 9.512345e-05 
    ## ... Procrustes: rmse 0.0002387063  max resid 0.0003893015 
    ## ... Similar to previous best
    ## Run 13 stress 9.947612e-05 
    ## ... Procrustes: rmse 0.0002494419  max resid 0.0004082945 
    ## ... Similar to previous best
    ## Run 14 stress 9.9342e-05 
    ## ... Procrustes: rmse 0.0002492274  max resid 0.0004038786 
    ## ... Similar to previous best
    ## Run 15 stress 8.95296e-05 
    ## ... Procrustes: rmse 0.0001651472  max resid 0.0002743321 
    ## ... Similar to previous best
    ## Run 16 stress 9.456465e-05 
    ## ... Procrustes: rmse 0.0001593204  max resid 0.0002529274 
    ## ... Similar to previous best
    ## Run 17 stress 9.93254e-05 
    ## ... Procrustes: rmse 0.0002106116  max resid 0.0003462624 
    ## ... Similar to previous best
    ## Run 18 stress 8.314897e-05 
    ## ... Procrustes: rmse 0.0001596862  max resid 0.0002601827 
    ## ... Similar to previous best
    ## Run 19 stress 6.158655e-05 
    ## ... Procrustes: rmse 3.157101e-05  max resid 6.429141e-05 
    ## ... Similar to previous best
    ## Run 20 stress 8.610362e-05 
    ## ... Procrustes: rmse 0.0001137254  max resid 0.0001988245 
    ## ... Similar to previous best
    ## Run 21 stress 9.653441e-05 
    ## ... Procrustes: rmse 0.0001381285  max resid 0.0002158293 
    ## ... Similar to previous best
    ## Run 22 stress 9.804031e-05 
    ## ... Procrustes: rmse 0.0002040825  max resid 0.0003339943 
    ## ... Similar to previous best
    ## Run 23 stress 9.80191e-05 
    ## ... Procrustes: rmse 0.0002170617  max resid 0.0003564733 
    ## ... Similar to previous best
    ## Run 24 stress 9.79745e-05 
    ## ... Procrustes: rmse 0.0002415843  max resid 0.0003992988 
    ## ... Similar to previous best
    ## Run 25 stress 9.288398e-05 
    ## ... Procrustes: rmse 0.0002037209  max resid 0.0003358522 
    ## ... Similar to previous best
    ## Run 26 stress 9.289194e-05 
    ## ... Procrustes: rmse 5.284828e-05  max resid 7.963111e-05 
    ## ... Similar to previous best
    ## Run 27 stress 9.825688e-05 
    ## ... Procrustes: rmse 0.0002337668  max resid 0.0003845612 
    ## ... Similar to previous best
    ## Run 28 stress 9.644873e-05 
    ## ... Procrustes: rmse 0.0001280076  max resid 0.0002147438 
    ## ... Similar to previous best
    ## Run 29 stress 9.779982e-05 
    ## ... Procrustes: rmse 0.0002323003  max resid 0.0003866525 
    ## ... Similar to previous best
    ## Run 30 stress 9.019179e-05 
    ## ... Procrustes: rmse 6.262557e-05  max resid 0.0001022857 
    ## ... Similar to previous best
    ## Run 31 stress 0.382571 
    ## Run 32 stress 8.283099e-05 
    ## ... Procrustes: rmse 4.847527e-05  max resid 8.631527e-05 
    ## ... Similar to previous best
    ## Run 33 stress 9.71858e-05 
    ## ... Procrustes: rmse 0.000236938  max resid 0.0003916341 
    ## ... Similar to previous best
    ## Run 34 stress 9.578373e-05 
    ## ... Procrustes: rmse 0.0002384526  max resid 0.0003953894 
    ## ... Similar to previous best
    ## Run 35 stress 8.996675e-05 
    ## ... Procrustes: rmse 0.0001242929  max resid 0.0002003942 
    ## ... Similar to previous best
    ## Run 36 stress 9.770962e-05 
    ## ... Procrustes: rmse 0.0002084342  max resid 0.0003442304 
    ## ... Similar to previous best
    ## Run 37 stress 9.799722e-05 
    ## ... Procrustes: rmse 0.0002443383  max resid 0.0003970018 
    ## ... Similar to previous best
    ## Run 38 stress 9.854132e-05 
    ## ... Procrustes: rmse 0.0002411529  max resid 0.0003977975 
    ## ... Similar to previous best
    ## Run 39 stress 9.411281e-05 
    ## ... Procrustes: rmse 0.000204972  max resid 0.000336529 
    ## ... Similar to previous best
    ## Run 40 stress 9.334009e-05 
    ## ... Procrustes: rmse 0.0002182388  max resid 0.0003604349 
    ## ... Similar to previous best
    ## Run 41 stress 9.303652e-05 
    ## ... Procrustes: rmse 0.0001699987  max resid 0.0002867005 
    ## ... Similar to previous best
    ## Run 42 stress 9.868557e-05 
    ## ... Procrustes: rmse 0.0002476274  max resid 0.0004089223 
    ## ... Similar to previous best
    ## Run 43 stress 9.888411e-05 
    ## ... Procrustes: rmse 0.0002464604  max resid 0.0004113609 
    ## ... Similar to previous best
    ## Run 44 stress 9.636929e-05 
    ## ... Procrustes: rmse 0.0002274353  max resid 0.000383505 
    ## ... Similar to previous best
    ## Run 45 stress 9.850976e-05 
    ## ... Procrustes: rmse 0.0002104055  max resid 0.0003464751 
    ## ... Similar to previous best
    ## Run 46 stress 6.694871e-05 
    ## ... Procrustes: rmse 3.1313e-05  max resid 6.364058e-05 
    ## ... Similar to previous best
    ## Run 47 stress 4.915884e-05 
    ## ... Procrustes: rmse 3.08672e-05  max resid 6.688709e-05 
    ## ... Similar to previous best
    ## Run 48 stress 9.829225e-05 
    ## ... Procrustes: rmse 0.0002024027  max resid 0.0003364112 
    ## ... Similar to previous best
    ## Run 49 stress 9.969631e-05 
    ## ... Procrustes: rmse 0.0002483357  max resid 0.0004058804 
    ## ... Similar to previous best
    ## Run 50 stress 9.880902e-05 
    ## ... Procrustes: rmse 0.0002369094  max resid 0.0003917151 
    ## ... Similar to previous best
    ## Run 51 stress 8.152615e-05 
    ## ... Procrustes: rmse 6.466118e-05  max resid 9.801308e-05 
    ## ... Similar to previous best
    ## Run 52 stress 9.590737e-05 
    ## ... Procrustes: rmse 0.000238876  max resid 0.0004027543 
    ## ... Similar to previous best
    ## Run 53 stress 9.982493e-05 
    ## ... Procrustes: rmse 0.0002214121  max resid 0.0003641636 
    ## ... Similar to previous best
    ## Run 54 stress 9.921925e-05 
    ## ... Procrustes: rmse 0.0002426776  max resid 0.0003987094 
    ## ... Similar to previous best
    ## Run 55 stress 9.208381e-05 
    ## ... Procrustes: rmse 5.110111e-05  max resid 9.422281e-05 
    ## ... Similar to previous best
    ## Run 56 stress 9.72805e-05 
    ## ... Procrustes: rmse 0.0002422203  max resid 0.0004049872 
    ## ... Similar to previous best
    ## Run 57 stress 4.928227e-05 
    ## ... Procrustes: rmse 3.032951e-05  max resid 5.287216e-05 
    ## ... Similar to previous best
    ## Run 58 stress 9.350142e-05 
    ## ... Procrustes: rmse 0.0001811402  max resid 0.000306386 
    ## ... Similar to previous best
    ## Run 59 stress 9.756577e-05 
    ## ... Procrustes: rmse 0.0002458884  max resid 0.0004044289 
    ## ... Similar to previous best
    ## Run 60 stress 9.610247e-05 
    ## ... Procrustes: rmse 0.0002130002  max resid 0.0003498974 
    ## ... Similar to previous best
    ## Run 61 stress 9.911358e-05 
    ## ... Procrustes: rmse 0.0002418744  max resid 0.0004005053 
    ## ... Similar to previous best
    ## Run 62 stress 9.522548e-05 
    ## ... Procrustes: rmse 0.0002392477  max resid 0.0003911721 
    ## ... Similar to previous best
    ## Run 63 stress 9.865158e-05 
    ## ... Procrustes: rmse 0.0002174773  max resid 0.0003603947 
    ## ... Similar to previous best
    ## Run 64 stress 4.670718e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 3.21358e-05  max resid 6.165586e-05 
    ## ... Similar to previous best
    ## Run 65 stress 6.403231e-05 
    ## ... Procrustes: rmse 6.268175e-05  max resid 0.0001008038 
    ## ... Similar to previous best
    ## Run 66 stress 9.836075e-05 
    ## ... Procrustes: rmse 0.0002566714  max resid 0.0004334765 
    ## ... Similar to previous best
    ## Run 67 stress 7.716629e-05 
    ## ... Procrustes: rmse 4.384435e-05  max resid 8.853041e-05 
    ## ... Similar to previous best
    ## Run 68 stress 9.78516e-05 
    ## ... Procrustes: rmse 0.0002413771  max resid 0.0004257394 
    ## ... Similar to previous best
    ## Run 69 stress 7.740055e-05 
    ## ... Procrustes: rmse 4.399837e-05  max resid 0.0001062123 
    ## ... Similar to previous best
    ## Run 70 stress 9.882456e-05 
    ## ... Procrustes: rmse 0.0002249247  max resid 0.0003975749 
    ## ... Similar to previous best
    ## Run 71 stress 9.980833e-05 
    ## ... Procrustes: rmse 0.0002608725  max resid 0.0004437512 
    ## ... Similar to previous best
    ## Run 72 stress 9.909477e-05 
    ## ... Procrustes: rmse 0.0002594594  max resid 0.0004355752 
    ## ... Similar to previous best
    ## Run 73 stress 9.180564e-05 
    ## ... Procrustes: rmse 0.0001956956  max resid 0.0003445909 
    ## ... Similar to previous best
    ## Run 74 stress 9.465256e-05 
    ## ... Procrustes: rmse 0.0002285572  max resid 0.0003992803 
    ## ... Similar to previous best
    ## Run 75 stress 9.685446e-05 
    ## ... Procrustes: rmse 0.0001995371  max resid 0.0003517028 
    ## ... Similar to previous best
    ## Run 76 stress 9.807635e-05 
    ## ... Procrustes: rmse 0.0002466607  max resid 0.0004320201 
    ## ... Similar to previous best
    ## Run 77 stress 9.87442e-05 
    ## ... Procrustes: rmse 0.0002570125  max resid 0.0004363043 
    ## ... Similar to previous best
    ## Run 78 stress 9.878691e-05 
    ## ... Procrustes: rmse 0.0002416867  max resid 0.0004206346 
    ## ... Similar to previous best
    ## Run 79 stress 9.859136e-05 
    ## ... Procrustes: rmse 0.000249534  max resid 0.0004353744 
    ## ... Similar to previous best
    ## Run 80 stress 9.963339e-05 
    ## ... Procrustes: rmse 0.0002581438  max resid 0.000453291 
    ## ... Similar to previous best
    ## Run 81 stress 9.985759e-05 
    ## ... Procrustes: rmse 0.0002587768  max resid 0.0004439172 
    ## ... Similar to previous best
    ## Run 82 stress 9.946917e-05 
    ## ... Procrustes: rmse 0.0002294503  max resid 0.0003967486 
    ## ... Similar to previous best
    ## Run 83 stress 9.429905e-05 
    ## ... Procrustes: rmse 4.930101e-05  max resid 0.0001048142 
    ## ... Similar to previous best
    ## Run 84 stress 9.737895e-05 
    ## ... Procrustes: rmse 0.0001845068  max resid 0.0003287742 
    ## ... Similar to previous best
    ## Run 85 stress 9.601633e-05 
    ## ... Procrustes: rmse 0.000214939  max resid 0.0003846337 
    ## ... Similar to previous best
    ## Run 86 stress 9.726004e-05 
    ## ... Procrustes: rmse 0.0002483104  max resid 0.0004360258 
    ## ... Similar to previous best
    ## Run 87 stress 9.565289e-05 
    ## ... Procrustes: rmse 0.0002124615  max resid 0.000378355 
    ## ... Similar to previous best
    ## Run 88 stress 5.880901e-05 
    ## ... Procrustes: rmse 3.930128e-05  max resid 7.584394e-05 
    ## ... Similar to previous best
    ## Run 89 stress 9.760037e-05 
    ## ... Procrustes: rmse 0.0001995873  max resid 0.0003536429 
    ## ... Similar to previous best
    ## Run 90 stress 9.811044e-05 
    ## ... Procrustes: rmse 0.0002571653  max resid 0.0004355441 
    ## ... Similar to previous best
    ## Run 91 stress 9.834241e-05 
    ## ... Procrustes: rmse 0.0002010077  max resid 0.000359677 
    ## ... Similar to previous best
    ## Run 92 stress 9.888637e-05 
    ## ... Procrustes: rmse 0.00025073  max resid 0.0004271219 
    ## ... Similar to previous best
    ## Run 93 stress 9.823125e-05 
    ## ... Procrustes: rmse 0.0002572466  max resid 0.0004379506 
    ## ... Similar to previous best
    ## Run 94 stress 7.511997e-05 
    ## ... Procrustes: rmse 4.130082e-05  max resid 8.012228e-05 
    ## ... Similar to previous best
    ## Run 95 stress 9.741706e-05 
    ## ... Procrustes: rmse 0.0002509082  max resid 0.0004243394 
    ## ... Similar to previous best
    ## Run 96 stress 5.522683e-05 
    ## ... Procrustes: rmse 3.988077e-05  max resid 7.099384e-05 
    ## ... Similar to previous best
    ## Run 97 stress 8.551389e-05 
    ## ... Procrustes: rmse 5.026503e-05  max resid 0.0001073858 
    ## ... Similar to previous best
    ## Run 98 stress 9.65177e-05 
    ## ... Procrustes: rmse 0.0002488717  max resid 0.0004190922 
    ## ... Similar to previous best
    ## Run 99 stress 9.557674e-05 
    ## ... Procrustes: rmse 0.0002452357  max resid 0.0004330193 
    ## ... Similar to previous best
    ## Run 100 stress 9.51602e-05 
    ## ... Procrustes: rmse 0.0002133642  max resid 0.0003555534 
    ## ... Similar to previous best
    ## Run 101 stress 9.814086e-05 
    ## ... Procrustes: rmse 0.0002108121  max resid 0.0003606075 
    ## ... Similar to previous best
    ## Run 102 stress 0.3770104 
    ## Run 103 stress 9.550809e-05 
    ## ... Procrustes: rmse 0.0001936029  max resid 0.0003476997 
    ## ... Similar to previous best
    ## Run 104 stress 9.898804e-05 
    ## ... Procrustes: rmse 0.0002524873  max resid 0.000444628 
    ## ... Similar to previous best
    ## Run 105 stress 8.388021e-05 
    ## ... Procrustes: rmse 4.83629e-05  max resid 8.456333e-05 
    ## ... Similar to previous best
    ## Run 106 stress 7.687174e-05 
    ## ... Procrustes: rmse 5.1766e-05  max resid 0.0001030321 
    ## ... Similar to previous best
    ## Run 107 stress 9.859119e-05 
    ## ... Procrustes: rmse 0.0002585932  max resid 0.0004368222 
    ## ... Similar to previous best
    ## Run 108 stress 9.846418e-05 
    ## ... Procrustes: rmse 0.0002245975  max resid 0.0003996099 
    ## ... Similar to previous best
    ## Run 109 stress 9.861284e-05 
    ## ... Procrustes: rmse 0.0002193364  max resid 0.0003761155 
    ## ... Similar to previous best
    ## Run 110 stress 9.707218e-05 
    ## ... Procrustes: rmse 0.000171178  max resid 0.000306487 
    ## ... Similar to previous best
    ## Run 111 stress 9.807134e-05 
    ## ... Procrustes: rmse 0.0002568764  max resid 0.0004328448 
    ## ... Similar to previous best
    ## Run 112 stress 6.962187e-05 
    ## ... Procrustes: rmse 4.04939e-05  max resid 7.215878e-05 
    ## ... Similar to previous best
    ## Run 113 stress 8.056759e-05 
    ## ... Procrustes: rmse 6.376575e-05  max resid 0.0001140706 
    ## ... Similar to previous best
    ## Run 114 stress 9.844968e-05 
    ## ... Procrustes: rmse 0.0002221481  max resid 0.0003771669 
    ## ... Similar to previous best
    ## Run 115 stress 9.973603e-05 
    ## ... Procrustes: rmse 0.0002564857  max resid 0.0004469164 
    ## ... Similar to previous best
    ## Run 116 stress 9.696873e-05 
    ## ... Procrustes: rmse 0.0002287849  max resid 0.0004021855 
    ## ... Similar to previous best
    ## Run 117 stress 8.855888e-05 
    ## ... Procrustes: rmse 0.0001003026  max resid 0.0001664709 
    ## ... Similar to previous best
    ## Run 118 stress 5.903356e-05 
    ## ... Procrustes: rmse 6.976667e-05  max resid 0.0001267489 
    ## ... Similar to previous best
    ## Run 119 stress 9.731573e-05 
    ## ... Procrustes: rmse 0.000195281  max resid 0.000340048 
    ## ... Similar to previous best
    ## Run 120 stress 9.857359e-05 
    ## ... Procrustes: rmse 0.000254553  max resid 0.0004453929 
    ## ... Similar to previous best
    ## Run 121 stress 9.943013e-05 
    ## ... Procrustes: rmse 0.0002203638  max resid 0.0003943105 
    ## ... Similar to previous best
    ## Run 122 stress 9.603309e-05 
    ## ... Procrustes: rmse 0.0001877253  max resid 0.0003359938 
    ## ... Similar to previous best
    ## Run 123 stress 9.932358e-05 
    ## ... Procrustes: rmse 0.0002563509  max resid 0.0004306863 
    ## ... Similar to previous best
    ## Run 124 stress 9.995381e-05 
    ## ... Procrustes: rmse 0.0002112764  max resid 0.0003781663 
    ## ... Similar to previous best
    ## Run 125 stress 9.761174e-05 
    ## ... Procrustes: rmse 0.0002219965  max resid 0.0003947124 
    ## ... Similar to previous best
    ## Run 126 stress 9.686751e-05 
    ## ... Procrustes: rmse 0.0001575471  max resid 0.0002821265 
    ## ... Similar to previous best
    ## Run 127 stress 9.55538e-05 
    ## ... Procrustes: rmse 0.0001641682  max resid 0.0003050679 
    ## ... Similar to previous best
    ## Run 128 stress 9.939685e-05 
    ## ... Procrustes: rmse 0.0002144649  max resid 0.0003798575 
    ## ... Similar to previous best
    ## Run 129 stress 5.909996e-05 
    ## ... Procrustes: rmse 3.911955e-05  max resid 6.744582e-05 
    ## ... Similar to previous best
    ## Run 130 stress 8.38755e-05 
    ## ... Procrustes: rmse 4.536289e-05  max resid 9.013188e-05 
    ## ... Similar to previous best
    ## Run 131 stress 9.718928e-05 
    ## ... Procrustes: rmse 0.0002546365  max resid 0.0004300298 
    ## ... Similar to previous best
    ## Run 132 stress 7.703313e-05 
    ## ... Procrustes: rmse 4.52292e-05  max resid 7.877665e-05 
    ## ... Similar to previous best
    ## Run 133 stress 9.813246e-05 
    ## ... Procrustes: rmse 0.0002214112  max resid 0.0003952587 
    ## ... Similar to previous best
    ## Run 134 stress 6.113877e-05 
    ## ... Procrustes: rmse 5.491079e-05  max resid 8.746675e-05 
    ## ... Similar to previous best
    ## Run 135 stress 9.943104e-05 
    ## ... Procrustes: rmse 0.0002601179  max resid 0.0004436412 
    ## ... Similar to previous best
    ## Run 136 stress 9.931355e-05 
    ## ... Procrustes: rmse 0.0002242245  max resid 0.0003975575 
    ## ... Similar to previous best
    ## Run 137 stress 9.708663e-05 
    ## ... Procrustes: rmse 0.0002444595  max resid 0.0004335056 
    ## ... Similar to previous best
    ## Run 138 stress 9.80084e-05 
    ## ... Procrustes: rmse 0.0002549798  max resid 0.0004473317 
    ## ... Similar to previous best
    ## Run 139 stress 9.944275e-05 
    ## ... Procrustes: rmse 0.0002213853  max resid 0.0003651079 
    ## ... Similar to previous best
    ## Run 140 stress 5.990919e-05 
    ## ... Procrustes: rmse 3.702166e-05  max resid 6.704365e-05 
    ## ... Similar to previous best
    ## Run 141 stress 7.340613e-05 
    ## ... Procrustes: rmse 4.890033e-05  max resid 9.065803e-05 
    ## ... Similar to previous best
    ## Run 142 stress 9.878734e-05 
    ## ... Procrustes: rmse 0.0002588941  max resid 0.0004369493 
    ## ... Similar to previous best
    ## Run 143 stress 9.301556e-05 
    ## ... Procrustes: rmse 0.000187768  max resid 0.0003373304 
    ## ... Similar to previous best
    ## Run 144 stress 9.934202e-05 
    ## ... Procrustes: rmse 0.0002563175  max resid 0.0004290375 
    ## ... Similar to previous best
    ## Run 145 stress 9.592646e-05 
    ## ... Procrustes: rmse 0.0001749797  max resid 0.0003168563 
    ## ... Similar to previous best
    ## Run 146 stress 9.602719e-05 
    ## ... Procrustes: rmse 0.0002448978  max resid 0.000423059 
    ## ... Similar to previous best
    ## Run 147 stress 9.934218e-05 
    ## ... Procrustes: rmse 0.0002244827  max resid 0.0003785994 
    ## ... Similar to previous best
    ## Run 148 stress 9.996457e-05 
    ## ... Procrustes: rmse 0.0002523099  max resid 0.0004376918 
    ## ... Similar to previous best
    ## Run 149 stress 7.215193e-05 
    ## ... Procrustes: rmse 4.357929e-05  max resid 8.176583e-05 
    ## ... Similar to previous best
    ## Run 150 stress 7.550848e-05 
    ## ... Procrustes: rmse 8.197857e-05  max resid 0.0001468989 
    ## ... Similar to previous best
    ## Run 151 stress 9.854622e-05 
    ## ... Procrustes: rmse 0.0002577044  max resid 0.0004432866 
    ## ... Similar to previous best
    ## Run 152 stress 9.995382e-05 
    ## ... Procrustes: rmse 0.0002267965  max resid 0.0003870847 
    ## ... Similar to previous best
    ## Run 153 stress 9.900392e-05 
    ## ... Procrustes: rmse 0.0002118373  max resid 0.0003635937 
    ## ... Similar to previous best
    ## Run 154 stress 9.943277e-05 
    ## ... Procrustes: rmse 0.0002212734  max resid 0.0003949168 
    ## ... Similar to previous best
    ## Run 155 stress 9.807671e-05 
    ## ... Procrustes: rmse 0.0002553247  max resid 0.0004325249 
    ## ... Similar to previous best
    ## Run 156 stress 5.539466e-05 
    ## ... Procrustes: rmse 3.533081e-05  max resid 6.551665e-05 
    ## ... Similar to previous best
    ## Run 157 stress 9.878899e-05 
    ## ... Procrustes: rmse 0.0002064168  max resid 0.000373055 
    ## ... Similar to previous best
    ## Run 158 stress 9.966338e-05 
    ## ... Procrustes: rmse 0.0002077876  max resid 0.0003569257 
    ## ... Similar to previous best
    ## Run 159 stress 9.933098e-05 
    ## ... Procrustes: rmse 0.0002575357  max resid 0.00044822 
    ## ... Similar to previous best
    ## Run 160 stress 9.906722e-05 
    ## ... Procrustes: rmse 0.0002588823  max resid 0.0004412815 
    ## ... Similar to previous best
    ## Run 161 stress 9.686723e-05 
    ## ... Procrustes: rmse 0.0002055353  max resid 0.0003762854 
    ## ... Similar to previous best
    ## Run 162 stress 6.43772e-05 
    ## ... Procrustes: rmse 8.666334e-05  max resid 0.0001717725 
    ## ... Similar to previous best
    ## Run 163 stress 9.823649e-05 
    ## ... Procrustes: rmse 0.0002152717  max resid 0.0003786362 
    ## ... Similar to previous best
    ## Run 164 stress 9.682216e-05 
    ## ... Procrustes: rmse 0.0002485547  max resid 0.0004346759 
    ## ... Similar to previous best
    ## Run 165 stress 8.387031e-05 
    ## ... Procrustes: rmse 8.742151e-05  max resid 0.0001235212 
    ## ... Similar to previous best
    ## Run 166 stress 9.267736e-05 
    ## ... Procrustes: rmse 0.0001680677  max resid 0.0003011885 
    ## ... Similar to previous best
    ## Run 167 stress 9.888087e-05 
    ## ... Procrustes: rmse 0.0002584442  max resid 0.0004484839 
    ## ... Similar to previous best
    ## Run 168 stress 9.923047e-05 
    ## ... Procrustes: rmse 0.0002137919  max resid 0.0003599254 
    ## ... Similar to previous best
    ## Run 169 stress 9.557526e-05 
    ## ... Procrustes: rmse 0.0002068001  max resid 0.000346566 
    ## ... Similar to previous best
    ## Run 170 stress 9.056797e-05 
    ## ... Procrustes: rmse 0.0001577897  max resid 0.0002908924 
    ## ... Similar to previous best
    ## Run 171 stress 9.765512e-05 
    ## ... Procrustes: rmse 0.0002456618  max resid 0.0004236243 
    ## ... Similar to previous best
    ## Run 172 stress 9.419824e-05 
    ## ... Procrustes: rmse 0.0002093796  max resid 0.0003634082 
    ## ... Similar to previous best
    ## Run 173 stress 8.669694e-05 
    ## ... Procrustes: rmse 0.0001367334  max resid 0.0002346411 
    ## ... Similar to previous best
    ## Run 174 stress 9.725989e-05 
    ## ... Procrustes: rmse 0.0002533014  max resid 0.0004353046 
    ## ... Similar to previous best
    ## Run 175 stress 9.383793e-05 
    ## ... Procrustes: rmse 5.943068e-05  max resid 9.910399e-05 
    ## ... Similar to previous best
    ## Run 176 stress 9.913526e-05 
    ## ... Procrustes: rmse 0.0002128208  max resid 0.0003674422 
    ## ... Similar to previous best
    ## Run 177 stress 9.957397e-05 
    ## ... Procrustes: rmse 0.0002472855  max resid 0.0004364281 
    ## ... Similar to previous best
    ## Run 178 stress 9.831987e-05 
    ## ... Procrustes: rmse 0.000226365  max resid 0.0003969464 
    ## ... Similar to previous best
    ## Run 179 stress 9.962846e-05 
    ## ... Procrustes: rmse 0.0002274853  max resid 0.0003962058 
    ## ... Similar to previous best
    ## Run 180 stress 9.903127e-05 
    ## ... Procrustes: rmse 0.0002580684  max resid 0.0004480666 
    ## ... Similar to previous best
    ## Run 181 stress 9.781175e-05 
    ## ... Procrustes: rmse 0.0002536692  max resid 0.0004226842 
    ## ... Similar to previous best
    ## Run 182 stress 8.037894e-05 
    ## ... Procrustes: rmse 5.029399e-05  max resid 8.152678e-05 
    ## ... Similar to previous best
    ## Run 183 stress 9.967943e-05 
    ## ... Procrustes: rmse 0.0002145058  max resid 0.0003573007 
    ## ... Similar to previous best
    ## Run 184 stress 9.733638e-05 
    ## ... Procrustes: rmse 0.0001978235  max resid 0.0003538913 
    ## ... Similar to previous best
    ## Run 185 stress 9.881314e-05 
    ## ... Procrustes: rmse 0.0002578043  max resid 0.0004415613 
    ## ... Similar to previous best
    ## Run 186 stress 7.592188e-05 
    ## ... Procrustes: rmse 8.478123e-05  max resid 0.0001323907 
    ## ... Similar to previous best
    ## Run 187 stress 9.950595e-05 
    ## ... Procrustes: rmse 0.0002519155  max resid 0.0004332871 
    ## ... Similar to previous best
    ## Run 188 stress 9.570871e-05 
    ## ... Procrustes: rmse 0.0002158603  max resid 0.0003637656 
    ## ... Similar to previous best
    ## Run 189 stress 9.817619e-05 
    ## ... Procrustes: rmse 0.0002478506  max resid 0.0004390225 
    ## ... Similar to previous best
    ## Run 190 stress 9.72468e-05 
    ## ... Procrustes: rmse 0.0002551362  max resid 0.0004309721 
    ## ... Similar to previous best
    ## Run 191 stress 9.819019e-05 
    ## ... Procrustes: rmse 5.06255e-05  max resid 8.184417e-05 
    ## ... Similar to previous best
    ## Run 192 stress 9.690747e-05 
    ## ... Procrustes: rmse 0.0002236518  max resid 0.0003941266 
    ## ... Similar to previous best
    ## Run 193 stress 8.503247e-05 
    ## ... Procrustes: rmse 8.36177e-05  max resid 0.0001413921 
    ## ... Similar to previous best
    ## Run 194 stress 9.709991e-05 
    ## ... Procrustes: rmse 0.0002534619  max resid 0.0004256063 
    ## ... Similar to previous best
    ## Run 195 stress 9.712996e-05 
    ## ... Procrustes: rmse 0.0001971569  max resid 0.0003580301 
    ## ... Similar to previous best
    ## Run 196 stress 9.969798e-05 
    ## ... Procrustes: rmse 0.0002602834  max resid 0.0004385056 
    ## ... Similar to previous best
    ## Run 197 stress 9.836061e-05 
    ## ... Procrustes: rmse 0.00025541  max resid 0.0004439155 
    ## ... Similar to previous best
    ## Run 198 stress 9.968554e-05 
    ## ... Procrustes: rmse 0.0002231005  max resid 0.000397298 
    ## ... Similar to previous best
    ## Run 199 stress 5.67463e-05 
    ## ... Procrustes: rmse 6.954236e-05  max resid 0.0001162474 
    ## ... Similar to previous best
    ## Run 200 stress 9.703843e-05 
    ## ... Procrustes: rmse 0.0002529011  max resid 0.0004261642 
    ## ... Similar to previous best
    ## Run 201 stress 9.795382e-05 
    ## ... Procrustes: rmse 0.0002201863  max resid 0.0003726898 
    ## ... Similar to previous best
    ## Run 202 stress 9.67353e-05 
    ## ... Procrustes: rmse 0.0001992141  max resid 0.0003343544 
    ## ... Similar to previous best
    ## Run 203 stress 9.994795e-05 
    ## ... Procrustes: rmse 0.0002603779  max resid 0.0004492169 
    ## ... Similar to previous best
    ## Run 204 stress 8.837481e-05 
    ## ... Procrustes: rmse 0.0001503652  max resid 0.0002587942 
    ## ... Similar to previous best
    ## Run 205 stress 9.023054e-05 
    ## ... Procrustes: rmse 0.0001666975  max resid 0.000311573 
    ## ... Similar to previous best
    ## Run 206 stress 9.023029e-05 
    ## ... Procrustes: rmse 0.0001944557  max resid 0.0003532113 
    ## ... Similar to previous best
    ## Run 207 stress 9.605039e-05 
    ## ... Procrustes: rmse 0.0001937106  max resid 0.000333902 
    ## ... Similar to previous best
    ## Run 208 stress 9.994215e-05 
    ## ... Procrustes: rmse 0.0002615626  max resid 0.0004399697 
    ## ... Similar to previous best
    ## Run 209 stress 9.632618e-05 
    ## ... Procrustes: rmse 0.000219078  max resid 0.0003850966 
    ## ... Similar to previous best
    ## Run 210 stress 9.850129e-05 
    ## ... Procrustes: rmse 0.0001695204  max resid 0.0003081427 
    ## ... Similar to previous best
    ## Run 211 stress 9.163974e-05 
    ## ... Procrustes: rmse 0.0002093869  max resid 0.0003801113 
    ## ... Similar to previous best
    ## Run 212 stress 9.673721e-05 
    ## ... Procrustes: rmse 0.000254011  max resid 0.0004317099 
    ## ... Similar to previous best
    ## Run 213 stress 9.938891e-05 
    ## ... Procrustes: rmse 0.0001136319  max resid 0.0002087844 
    ## ... Similar to previous best
    ## Run 214 stress 9.510335e-05 
    ## ... Procrustes: rmse 5.441832e-05  max resid 8.231511e-05 
    ## ... Similar to previous best
    ## Run 215 stress 9.940937e-05 
    ## ... Procrustes: rmse 0.0002222808  max resid 0.0003950191 
    ## ... Similar to previous best
    ## Run 216 stress 6.202669e-05 
    ## ... Procrustes: rmse 4.150588e-05  max resid 6.991864e-05 
    ## ... Similar to previous best
    ## Run 217 stress 6.879424e-05 
    ## ... Procrustes: rmse 4.723914e-05  max resid 8.632826e-05 
    ## ... Similar to previous best
    ## Run 218 stress 6.430925e-05 
    ## ... Procrustes: rmse 6.73709e-05  max resid 0.0001330892 
    ## ... Similar to previous best
    ## Run 219 stress 9.880544e-05 
    ## ... Procrustes: rmse 0.000213772  max resid 0.0003579188 
    ## ... Similar to previous best
    ## Run 220 stress 9.569596e-05 
    ## ... Procrustes: rmse 0.0002192402  max resid 0.0003836126 
    ## ... Similar to previous best
    ## Run 221 stress 9.677146e-05 
    ## ... Procrustes: rmse 0.0002065543  max resid 0.0003656449 
    ## ... Similar to previous best
    ## Run 222 stress 9.815793e-05 
    ## ... Procrustes: rmse 0.00025704  max resid 0.0004346981 
    ## ... Similar to previous best
    ## Run 223 stress 9.908008e-05 
    ## ... Procrustes: rmse 0.0002572045  max resid 0.0004513575 
    ## ... Similar to previous best
    ## Run 224 stress 9.876191e-05 
    ## ... Procrustes: rmse 5.104463e-05  max resid 8.985512e-05 
    ## ... Similar to previous best
    ## Run 225 stress 9.994792e-05 
    ## ... Procrustes: rmse 0.0002166367  max resid 0.0003908278 
    ## ... Similar to previous best
    ## Run 226 stress 6.310355e-05 
    ## ... Procrustes: rmse 3.900483e-05  max resid 6.082809e-05 
    ## ... Similar to previous best
    ## Run 227 stress 9.975593e-05 
    ## ... Procrustes: rmse 0.0002232689  max resid 0.0003894702 
    ## ... Similar to previous best
    ## Run 228 stress 7.59481e-05 
    ## ... Procrustes: rmse 3.843388e-05  max resid 6.352807e-05 
    ## ... Similar to previous best
    ## Run 229 stress 9.725718e-05 
    ## ... Procrustes: rmse 0.0002215938  max resid 0.000394819 
    ## ... Similar to previous best
    ## Run 230 stress 7.582671e-05 
    ## ... Procrustes: rmse 4.193062e-05  max resid 0.0001181248 
    ## ... Similar to previous best
    ## Run 231 stress 9.602167e-05 
    ## ... Procrustes: rmse 0.0002194378  max resid 0.0003853399 
    ## ... Similar to previous best
    ## Run 232 stress 9.969134e-05 
    ## ... Procrustes: rmse 0.0002481444  max resid 0.0004339032 
    ## ... Similar to previous best
    ## Run 233 stress 9.928153e-05 
    ## ... Procrustes: rmse 0.0002144119  max resid 0.00036805 
    ## ... Similar to previous best
    ## Run 234 stress 9.949504e-05 
    ## ... Procrustes: rmse 0.0001891856  max resid 0.0003346019 
    ## ... Similar to previous best
    ## Run 235 stress 9.911011e-05 
    ## ... Procrustes: rmse 0.0002413888  max resid 0.0004196213 
    ## ... Similar to previous best
    ## Run 236 stress 7.993128e-05 
    ## ... Procrustes: rmse 4.30849e-05  max resid 7.270719e-05 
    ## ... Similar to previous best
    ## Run 237 stress 9.764139e-05 
    ## ... Procrustes: rmse 0.0002198197  max resid 0.0003645655 
    ## ... Similar to previous best
    ## Run 238 stress 9.866197e-05 
    ## ... Procrustes: rmse 0.0002280913  max resid 0.0003874429 
    ## ... Similar to previous best
    ## Run 239 stress 9.863927e-05 
    ## ... Procrustes: rmse 0.0002365819  max resid 0.0004161302 
    ## ... Similar to previous best
    ## Run 240 stress 9.609283e-05 
    ## ... Procrustes: rmse 0.0002512502  max resid 0.0004260269 
    ## ... Similar to previous best
    ## Run 241 stress 7.432569e-05 
    ## ... Procrustes: rmse 4.800247e-05  max resid 6.841433e-05 
    ## ... Similar to previous best
    ## Run 242 stress 7.577824e-05 
    ## ... Procrustes: rmse 6.18176e-05  max resid 0.000127209 
    ## ... Similar to previous best
    ## Run 243 stress 8.641611e-05 
    ## ... Procrustes: rmse 4.901176e-05  max resid 8.239058e-05 
    ## ... Similar to previous best
    ## Run 244 stress 9.968226e-05 
    ## ... Procrustes: rmse 0.0001561026  max resid 0.0002700721 
    ## ... Similar to previous best
    ## Run 245 stress 9.03956e-05 
    ## ... Procrustes: rmse 0.0001851183  max resid 0.0003340402 
    ## ... Similar to previous best
    ## Run 246 stress 8.87167e-05 
    ## ... Procrustes: rmse 0.0002307239  max resid 0.000397451 
    ## ... Similar to previous best
    ## Run 247 stress 9.170042e-05 
    ## ... Procrustes: rmse 0.0001589122  max resid 0.0003034347 
    ## ... Similar to previous best
    ## Run 248 stress 9.950704e-05 
    ## ... Procrustes: rmse 0.0002574152  max resid 0.0004515978 
    ## ... Similar to previous best
    ## Run 249 stress 9.623267e-05 
    ## ... Procrustes: rmse 0.0002070018  max resid 0.0003600023 
    ## ... Similar to previous best
    ## Run 250 stress 9.749621e-05 
    ## ... Procrustes: rmse 0.0002396473  max resid 0.0004185898 
    ## ... Similar to previous best
    ## Run 251 stress 9.659263e-05 
    ## ... Procrustes: rmse 0.0002159978  max resid 0.0003625362 
    ## ... Similar to previous best
    ## Run 252 stress 9.788103e-05 
    ## ... Procrustes: rmse 0.0002191044  max resid 0.0003666039 
    ## ... Similar to previous best
    ## Run 253 stress 9.78223e-05 
    ## ... Procrustes: rmse 0.0002521012  max resid 0.0004350382 
    ## ... Similar to previous best
    ## Run 254 stress 9.990106e-05 
    ## ... Procrustes: rmse 0.0002600217  max resid 0.0004555094 
    ## ... Similar to previous best
    ## Run 255 stress 9.789138e-05 
    ## ... Procrustes: rmse 0.0002222106  max resid 0.0003864211 
    ## ... Similar to previous best
    ## Run 256 stress 9.942958e-05 
    ## ... Procrustes: rmse 0.0002511304  max resid 0.0004429514 
    ## ... Similar to previous best
    ## Run 257 stress 9.837621e-05 
    ## ... Procrustes: rmse 0.0002043359  max resid 0.0003677696 
    ## ... Similar to previous best
    ## Run 258 stress 9.886161e-05 
    ## ... Procrustes: rmse 0.0002586769  max resid 0.0004356922 
    ## ... Similar to previous best
    ## Run 259 stress 6.565811e-05 
    ## ... Procrustes: rmse 4.153231e-05  max resid 7.543764e-05 
    ## ... Similar to previous best
    ## Run 260 stress 8.401524e-05 
    ## ... Procrustes: rmse 5.681668e-05  max resid 0.0001340183 
    ## ... Similar to previous best
    ## Run 261 stress 9.317465e-05 
    ## ... Procrustes: rmse 0.0001870175  max resid 0.0003348516 
    ## ... Similar to previous best
    ## Run 262 stress 7.761954e-05 
    ## ... Procrustes: rmse 6.005836e-05  max resid 0.0001023901 
    ## ... Similar to previous best
    ## Run 263 stress 8.657161e-05 
    ## ... Procrustes: rmse 0.00014003  max resid 0.0002437297 
    ## ... Similar to previous best
    ## Run 264 stress 8.055827e-05 
    ## ... Procrustes: rmse 0.0001306807  max resid 0.0002549207 
    ## ... Similar to previous best
    ## Run 265 stress 9.780904e-05 
    ## ... Procrustes: rmse 0.0002225134  max resid 0.0003923807 
    ## ... Similar to previous best
    ## Run 266 stress 9.861057e-05 
    ## ... Procrustes: rmse 0.0002192273  max resid 0.0003901795 
    ## ... Similar to previous best
    ## Run 267 stress 8.808451e-05 
    ## ... Procrustes: rmse 7.182749e-05  max resid 0.0001278048 
    ## ... Similar to previous best
    ## Run 268 stress 7.339694e-05 
    ## ... Procrustes: rmse 8.029184e-05  max resid 0.0001187952 
    ## ... Similar to previous best
    ## Run 269 stress 8.685897e-05 
    ## ... Procrustes: rmse 4.380537e-05  max resid 8.979871e-05 
    ## ... Similar to previous best
    ## Run 270 stress 5.937778e-05 
    ## ... Procrustes: rmse 3.179013e-05  max resid 5.437949e-05 
    ## ... Similar to previous best
    ## Run 271 stress 9.871594e-05 
    ## ... Procrustes: rmse 0.0002575422  max resid 0.0004406223 
    ## ... Similar to previous best
    ## Run 272 stress 9.949051e-05 
    ## ... Procrustes: rmse 0.0002145198  max resid 0.0003773889 
    ## ... Similar to previous best
    ## Run 273 stress 8.893269e-05 
    ## ... Procrustes: rmse 0.0001482577  max resid 0.0002603998 
    ## ... Similar to previous best
    ## Run 274 stress 9.700546e-05 
    ## ... Procrustes: rmse 0.0002190734  max resid 0.0003899908 
    ## ... Similar to previous best
    ## Run 275 stress 9.266889e-05 
    ## ... Procrustes: rmse 0.0001388742  max resid 0.000246147 
    ## ... Similar to previous best
    ## Run 276 stress 9.264591e-05 
    ## ... Procrustes: rmse 0.0001641388  max resid 0.0003186073 
    ## ... Similar to previous best
    ## Run 277 stress 9.783546e-05 
    ## ... Procrustes: rmse 0.0002552499  max resid 0.0004411681 
    ## ... Similar to previous best
    ## Run 278 stress 9.997996e-05 
    ## ... Procrustes: rmse 0.000221478  max resid 0.0003889012 
    ## ... Similar to previous best
    ## Run 279 stress 9.606504e-05 
    ## ... Procrustes: rmse 0.0002494709  max resid 0.0004137777 
    ## ... Similar to previous best
    ## Run 280 stress 5.448803e-05 
    ## ... Procrustes: rmse 4.752158e-05  max resid 9.638087e-05 
    ## ... Similar to previous best
    ## Run 281 stress 9.09735e-05 
    ## ... Procrustes: rmse 4.455767e-05  max resid 0.0001043929 
    ## ... Similar to previous best
    ## Run 282 stress 0.3825773 
    ## Run 283 stress 9.673864e-05 
    ## ... Procrustes: rmse 0.0002165007  max resid 0.0003753119 
    ## ... Similar to previous best
    ## Run 284 stress 9.823071e-05 
    ## ... Procrustes: rmse 0.0002126633  max resid 0.0003765148 
    ## ... Similar to previous best
    ## Run 285 stress 6.776027e-05 
    ## ... Procrustes: rmse 4.56848e-05  max resid 7.387633e-05 
    ## ... Similar to previous best
    ## Run 286 stress 8.458451e-05 
    ## ... Procrustes: rmse 4.764067e-05  max resid 0.0001040114 
    ## ... Similar to previous best
    ## Run 287 stress 9.756042e-05 
    ## ... Procrustes: rmse 0.0002515972  max resid 0.0004377503 
    ## ... Similar to previous best
    ## Run 288 stress 9.28727e-05 
    ## ... Procrustes: rmse 0.0002293804  max resid 0.0004093301 
    ## ... Similar to previous best
    ## Run 289 stress 5.84685e-05 
    ## ... Procrustes: rmse 4.038108e-05  max resid 6.937308e-05 
    ## ... Similar to previous best
    ## Run 290 stress 7.075935e-05 
    ## ... Procrustes: rmse 9.570252e-05  max resid 0.0001501293 
    ## ... Similar to previous best
    ## Run 291 stress 9.92388e-05 
    ## ... Procrustes: rmse 0.0002469667  max resid 0.000399243 
    ## ... Similar to previous best
    ## Run 292 stress 9.707912e-05 
    ## ... Procrustes: rmse 0.0001729418  max resid 0.0003117222 
    ## ... Similar to previous best
    ## Run 293 stress 9.840706e-05 
    ## ... Procrustes: rmse 0.0002240676  max resid 0.0003988231 
    ## ... Similar to previous best
    ## Run 294 stress 7.758099e-05 
    ## ... Procrustes: rmse 4.658846e-05  max resid 9.792215e-05 
    ## ... Similar to previous best
    ## Run 295 stress 9.423729e-05 
    ## ... Procrustes: rmse 0.0002149954  max resid 0.0003760429 
    ## ... Similar to previous best
    ## Run 296 stress 8.190521e-05 
    ## ... Procrustes: rmse 4.36036e-05  max resid 7.884917e-05 
    ## ... Similar to previous best
    ## Run 297 stress 9.983798e-05 
    ## ... Procrustes: rmse 8.160454e-05  max resid 0.000154202 
    ## ... Similar to previous best
    ## Run 298 stress 5.317563e-05 
    ## ... Procrustes: rmse 5.278575e-05  max resid 8.332836e-05 
    ## ... Similar to previous best
    ## Run 299 stress 4.976132e-05 
    ## ... Procrustes: rmse 3.716701e-05  max resid 7.158333e-05 
    ## ... Similar to previous best
    ## Run 300 stress 9.996133e-05 
    ## ... Procrustes: rmse 0.0002598584  max resid 0.000448656 
    ## ... Similar to previous best
    ## Run 301 stress 9.622121e-05 
    ## ... Procrustes: rmse 0.0001755884  max resid 0.0003127776 
    ## ... Similar to previous best
    ## Run 302 stress 9.458457e-05 
    ## ... Procrustes: rmse 0.0002209409  max resid 0.0003923645 
    ## ... Similar to previous best
    ## Run 303 stress 9.500912e-05 
    ## ... Procrustes: rmse 0.0002491325  max resid 0.0004206436 
    ## ... Similar to previous best
    ## Run 304 stress 9.878587e-05 
    ## ... Procrustes: rmse 6.628966e-05  max resid 0.0001099147 
    ## ... Similar to previous best
    ## Run 305 stress 7.168417e-05 
    ## ... Procrustes: rmse 4.658063e-05  max resid 0.0001006162 
    ## ... Similar to previous best
    ## Run 306 stress 9.491175e-05 
    ## ... Procrustes: rmse 0.0001970893  max resid 0.0003498529 
    ## ... Similar to previous best
    ## Run 307 stress 7.929651e-05 
    ## ... Procrustes: rmse 4.719964e-05  max resid 7.666656e-05 
    ## ... Similar to previous best
    ## Run 308 stress 9.750022e-05 
    ## ... Procrustes: rmse 0.0001920212  max resid 0.0003502477 
    ## ... Similar to previous best
    ## Run 309 stress 9.560657e-05 
    ## ... Procrustes: rmse 0.0002160975  max resid 0.0003825818 
    ## ... Similar to previous best
    ## Run 310 stress 9.475691e-05 
    ## ... Procrustes: rmse 0.0002114335  max resid 0.0003772801 
    ## ... Similar to previous best
    ## Run 311 stress 4.493296e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 3.366757e-05  max resid 6.308311e-05 
    ## ... Similar to previous best
    ## Run 312 stress 6.602878e-05 
    ## ... Procrustes: rmse 3.614024e-05  max resid 6.710984e-05 
    ## ... Similar to previous best
    ## Run 313 stress 9.53551e-05 
    ## ... Procrustes: rmse 0.0001717169  max resid 0.000283106 
    ## ... Similar to previous best
    ## Run 314 stress 9.62421e-05 
    ## ... Procrustes: rmse 0.0002098465  max resid 0.0003189539 
    ## ... Similar to previous best
    ## Run 315 stress 9.80344e-05 
    ## ... Procrustes: rmse 0.0002432483  max resid 0.000415968 
    ## ... Similar to previous best
    ## Run 316 stress 9.71415e-05 
    ## ... Procrustes: rmse 0.0001947687  max resid 0.000303031 
    ## ... Similar to previous best
    ## Run 317 stress 9.984555e-05 
    ## ... Procrustes: rmse 0.0002344836  max resid 0.0003798685 
    ## ... Similar to previous best
    ## Run 318 stress 9.561711e-05 
    ## ... Procrustes: rmse 0.0002216927  max resid 0.0003542265 
    ## ... Similar to previous best
    ## Run 319 stress 9.693373e-05 
    ## ... Procrustes: rmse 0.0001968956  max resid 0.0002983834 
    ## ... Similar to previous best
    ## Run 320 stress 9.237142e-05 
    ## ... Procrustes: rmse 0.0002278734  max resid 0.0003912665 
    ## ... Similar to previous best
    ## Run 321 stress 8.994123e-05 
    ## ... Procrustes: rmse 7.482784e-05  max resid 0.0001300969 
    ## ... Similar to previous best
    ## Run 322 stress 9.682603e-05 
    ## ... Procrustes: rmse 0.0002426153  max resid 0.0004157472 
    ## ... Similar to previous best
    ## Run 323 stress 9.94771e-05 
    ## ... Procrustes: rmse 0.0002098441  max resid 0.0003218272 
    ## ... Similar to previous best
    ## Run 324 stress 7.24956e-05 
    ## ... Procrustes: rmse 3.190925e-05  max resid 5.534683e-05 
    ## ... Similar to previous best
    ## Run 325 stress 8.496871e-05 
    ## ... Procrustes: rmse 4.470537e-05  max resid 9.07105e-05 
    ## ... Similar to previous best
    ## Run 326 stress 9.937639e-05 
    ## ... Procrustes: rmse 0.0001920799  max resid 0.0003383625 
    ## ... Similar to previous best
    ## Run 327 stress 9.993679e-05 
    ## ... Procrustes: rmse 0.0002150828  max resid 0.0003290582 
    ## ... Similar to previous best
    ## Run 328 stress 9.597191e-05 
    ## ... Procrustes: rmse 0.0002388435  max resid 0.0004081271 
    ## ... Similar to previous best
    ## Run 329 stress 9.886007e-05 
    ## ... Procrustes: rmse 0.0002140716  max resid 0.0003432774 
    ## ... Similar to previous best
    ## Run 330 stress 9.882823e-05 
    ## ... Procrustes: rmse 0.0002048015  max resid 0.0003467235 
    ## ... Similar to previous best
    ## Run 331 stress 9.968714e-05 
    ## ... Procrustes: rmse 0.0001904214  max resid 0.0002957168 
    ## ... Similar to previous best
    ## Run 332 stress 9.541225e-05 
    ## ... Procrustes: rmse 0.0001899374  max resid 0.000290341 
    ## ... Similar to previous best
    ## Run 333 stress 9.613571e-05 
    ## ... Procrustes: rmse 0.0002331652  max resid 0.0003820766 
    ## ... Similar to previous best
    ## Run 334 stress 8.381097e-05 
    ## ... Procrustes: rmse 4.405655e-05  max resid 7.307088e-05 
    ## ... Similar to previous best
    ## Run 335 stress 9.945099e-05 
    ## ... Procrustes: rmse 0.0002113473  max resid 0.0003362247 
    ## ... Similar to previous best
    ## Run 336 stress 9.472787e-05 
    ## ... Procrustes: rmse 0.0002373352  max resid 0.0004114039 
    ## ... Similar to previous best
    ## Run 337 stress 9.441641e-05 
    ## ... Procrustes: rmse 5.367934e-05  max resid 0.0001043239 
    ## ... Similar to previous best
    ## Run 338 stress 9.954967e-05 
    ## ... Procrustes: rmse 0.0002401578  max resid 0.0003958449 
    ## ... Similar to previous best
    ## Run 339 stress 9.834887e-05 
    ## ... Procrustes: rmse 0.0002116115  max resid 0.0003403619 
    ## ... Similar to previous best
    ## Run 340 stress 9.75485e-05 
    ## ... Procrustes: rmse 0.0002038766  max resid 0.0003288625 
    ## ... Similar to previous best
    ## Run 341 stress 8.972914e-05 
    ## ... Procrustes: rmse 0.0001854792  max resid 0.0003147364 
    ## ... Similar to previous best
    ## Run 342 stress 9.819098e-05 
    ## ... Procrustes: rmse 0.0001753111  max resid 0.0002708611 
    ## ... Similar to previous best
    ## Run 343 stress 8.498618e-05 
    ## ... Procrustes: rmse 6.763051e-05  max resid 0.0001098752 
    ## ... Similar to previous best
    ## Run 344 stress 9.814475e-05 
    ## ... Procrustes: rmse 0.0001965775  max resid 0.000295382 
    ## ... Similar to previous best
    ## Run 345 stress 9.696372e-05 
    ## ... Procrustes: rmse 0.0001745201  max resid 0.0002814249 
    ## ... Similar to previous best
    ## Run 346 stress 9.604313e-05 
    ## ... Procrustes: rmse 0.0001834634  max resid 0.000289909 
    ## ... Similar to previous best
    ## Run 347 stress 9.590694e-05 
    ## ... Procrustes: rmse 0.0002005891  max resid 0.0003338158 
    ## ... Similar to previous best
    ## Run 348 stress 9.950752e-05 
    ## ... Procrustes: rmse 0.0002132912  max resid 0.0003349357 
    ## ... Similar to previous best
    ## Run 349 stress 9.905422e-05 
    ## ... Procrustes: rmse 0.0001791961  max resid 0.0003094817 
    ## ... Similar to previous best
    ## Run 350 stress 9.825911e-05 
    ## ... Procrustes: rmse 0.0002436267  max resid 0.0004129775 
    ## ... Similar to previous best
    ## Run 351 stress 9.855778e-05 
    ## ... Procrustes: rmse 0.0002445154  max resid 0.0004212184 
    ## ... Similar to previous best
    ## Run 352 stress 9.485096e-05 
    ## ... Procrustes: rmse 0.0002304124  max resid 0.0003820654 
    ## ... Similar to previous best
    ## Run 353 stress 9.895608e-05 
    ## ... Procrustes: rmse 0.0002324004  max resid 0.0003626339 
    ## ... Similar to previous best
    ## Run 354 stress 9.811375e-05 
    ## ... Procrustes: rmse 0.0002107872  max resid 0.0003323907 
    ## ... Similar to previous best
    ## Run 355 stress 9.868847e-05 
    ## ... Procrustes: rmse 0.0002147314  max resid 0.0003566877 
    ## ... Similar to previous best
    ## Run 356 stress 9.748221e-05 
    ## ... Procrustes: rmse 0.0001881423  max resid 0.0003016474 
    ## ... Similar to previous best
    ## Run 357 stress 9.928614e-05 
    ## ... Procrustes: rmse 0.0002425397  max resid 0.0004217328 
    ## ... Similar to previous best
    ## Run 358 stress 9.995836e-05 
    ## ... Procrustes: rmse 0.0002337165  max resid 0.0003708073 
    ## ... Similar to previous best
    ## Run 359 stress 7.93858e-05 
    ## ... Procrustes: rmse 4.292754e-05  max resid 8.780383e-05 
    ## ... Similar to previous best
    ## Run 360 stress 7.15484e-05 
    ## ... Procrustes: rmse 4.135661e-05  max resid 8.097454e-05 
    ## ... Similar to previous best
    ## Run 361 stress 9.963074e-05 
    ## ... Procrustes: rmse 0.0002352808  max resid 0.0003714114 
    ## ... Similar to previous best
    ## Run 362 stress 9.867749e-05 
    ## ... Procrustes: rmse 0.0002355998  max resid 0.000382054 
    ## ... Similar to previous best
    ## Run 363 stress 5.733583e-05 
    ## ... Procrustes: rmse 3.546633e-05  max resid 6.78607e-05 
    ## ... Similar to previous best
    ## Run 364 stress 8.795395e-05 
    ## ... Procrustes: rmse 0.0001156435  max resid 0.0001977916 
    ## ... Similar to previous best
    ## Run 365 stress 7.988189e-05 
    ## ... Procrustes: rmse 7.74246e-05  max resid 0.0001339206 
    ## ... Similar to previous best
    ## Run 366 stress 9.739656e-05 
    ## ... Procrustes: rmse 0.0002437733  max resid 0.0004146774 
    ## ... Similar to previous best
    ## Run 367 stress 9.684742e-05 
    ## ... Procrustes: rmse 0.000203636  max resid 0.0003253402 
    ## ... Similar to previous best
    ## Run 368 stress 9.603535e-05 
    ## ... Procrustes: rmse 0.0001812372  max resid 0.0003006163 
    ## ... Similar to previous best
    ## Run 369 stress 9.509869e-05 
    ## ... Procrustes: rmse 0.0002364048  max resid 0.0004042096 
    ## ... Similar to previous best
    ## Run 370 stress 9.781255e-05 
    ## ... Procrustes: rmse 0.0002043246  max resid 0.0003296123 
    ## ... Similar to previous best
    ## Run 371 stress 5.621972e-05 
    ## ... Procrustes: rmse 2.758835e-05  max resid 5.181271e-05 
    ## ... Similar to previous best
    ## Run 372 stress 9.967829e-05 
    ## ... Procrustes: rmse 0.0002011514  max resid 0.000311621 
    ## ... Similar to previous best
    ## Run 373 stress 9.627455e-05 
    ## ... Procrustes: rmse 0.0002353878  max resid 0.0003868594 
    ## ... Similar to previous best
    ## Run 374 stress 9.893116e-05 
    ## ... Procrustes: rmse 0.0002427754  max resid 0.0004004203 
    ## ... Similar to previous best
    ## Run 375 stress 5.800372e-05 
    ## ... Procrustes: rmse 3.443358e-05  max resid 8.02639e-05 
    ## ... Similar to previous best
    ## Run 376 stress 9.954843e-05 
    ## ... Procrustes: rmse 0.0001984784  max resid 0.0003026701 
    ## ... Similar to previous best
    ## Run 377 stress 9.141957e-05 
    ## ... Procrustes: rmse 0.0001150768  max resid 0.0001980786 
    ## ... Similar to previous best
    ## Run 378 stress 9.827421e-05 
    ## ... Procrustes: rmse 0.0001887056  max resid 0.0002862776 
    ## ... Similar to previous best
    ## Run 379 stress 7.627635e-05 
    ## ... Procrustes: rmse 4.441537e-05  max resid 7.242685e-05 
    ## ... Similar to previous best
    ## Run 380 stress 9.35896e-05 
    ## ... Procrustes: rmse 0.0002240642  max resid 0.0003530638 
    ## ... Similar to previous best
    ## Run 381 stress 8.177569e-05 
    ## ... Procrustes: rmse 0.0001139491  max resid 0.0001685093 
    ## ... Similar to previous best
    ## Run 382 stress 6.937018e-05 
    ## ... Procrustes: rmse 3.699008e-05  max resid 7.633524e-05 
    ## ... Similar to previous best
    ## Run 383 stress 9.926353e-05 
    ## ... Procrustes: rmse 0.0002104845  max resid 0.0003417124 
    ## ... Similar to previous best
    ## Run 384 stress 9.882351e-05 
    ## ... Procrustes: rmse 0.0002131867  max resid 0.0003604134 
    ## ... Similar to previous best
    ## Run 385 stress 6.270527e-05 
    ## ... Procrustes: rmse 3.431562e-05  max resid 7.096775e-05 
    ## ... Similar to previous best
    ## Run 386 stress 9.566066e-05 
    ## ... Procrustes: rmse 0.0001915561  max resid 0.0002914305 
    ## ... Similar to previous best
    ## Run 387 stress 9.798838e-05 
    ## ... Procrustes: rmse 0.0001682618  max resid 0.0002592075 
    ## ... Similar to previous best
    ## Run 388 stress 9.745037e-05 
    ## ... Procrustes: rmse 0.0002423448  max resid 0.0004174226 
    ## ... Similar to previous best
    ## Run 389 stress 9.998659e-05 
    ## ... Procrustes: rmse 0.0001815143  max resid 0.0002826328 
    ## ... Similar to previous best
    ## Run 390 stress 9.989947e-05 
    ## ... Procrustes: rmse 0.0002465298  max resid 0.0004030764 
    ## ... Similar to previous best
    ## Run 391 stress 7.158105e-05 
    ## ... Procrustes: rmse 9.304818e-05  max resid 0.0001589221 
    ## ... Similar to previous best
    ## Run 392 stress 9.501889e-05 
    ## ... Procrustes: rmse 4.919827e-05  max resid 7.79114e-05 
    ## ... Similar to previous best
    ## Run 393 stress 9.878269e-05 
    ## ... Procrustes: rmse 0.0002448836  max resid 0.000424489 
    ## ... Similar to previous best
    ## Run 394 stress 9.972546e-05 
    ## ... Procrustes: rmse 0.0002168648  max resid 0.0003566229 
    ## ... Similar to previous best
    ## Run 395 stress 9.9957e-05 
    ## ... Procrustes: rmse 0.0002138293  max resid 0.0003560842 
    ## ... Similar to previous best
    ## Run 396 stress 9.197934e-05 
    ## ... Procrustes: rmse 0.0001948207  max resid 0.0003185729 
    ## ... Similar to previous best
    ## Run 397 stress 9.955002e-05 
    ## ... Procrustes: rmse 0.0001847391  max resid 0.0002996074 
    ## ... Similar to previous best
    ## Run 398 stress 8.51019e-05 
    ## ... Procrustes: rmse 4.769941e-05  max resid 9.434931e-05 
    ## ... Similar to previous best
    ## Run 399 stress 9.691427e-05 
    ## ... Procrustes: rmse 0.0002022948  max resid 0.0003211839 
    ## ... Similar to previous best
    ## Run 400 stress 9.810066e-05 
    ## ... Procrustes: rmse 0.0001987083  max resid 0.0002960323 
    ## ... Similar to previous best
    ## Run 401 stress 9.967342e-05 
    ## ... Procrustes: rmse 0.0002015401  max resid 0.0003049501 
    ## ... Similar to previous best
    ## Run 402 stress 9.897474e-05 
    ## ... Procrustes: rmse 0.0002478533  max resid 0.0004229499 
    ## ... Similar to previous best
    ## Run 403 stress 9.901089e-05 
    ## ... Procrustes: rmse 0.0002413044  max resid 0.0003966845 
    ## ... Similar to previous best
    ## Run 404 stress 9.77192e-05 
    ## ... Procrustes: rmse 0.0002083375  max resid 0.0003412583 
    ## ... Similar to previous best
    ## Run 405 stress 9.636217e-05 
    ## ... Procrustes: rmse 0.0002404977  max resid 0.0004095764 
    ## ... Similar to previous best
    ## Run 406 stress 9.676585e-05 
    ## ... Procrustes: rmse 0.0002084223  max resid 0.0003373099 
    ## ... Similar to previous best
    ## Run 407 stress 9.502952e-05 
    ## ... Procrustes: rmse 0.0002330278  max resid 0.0004132035 
    ## ... Similar to previous best
    ## Run 408 stress 7.25261e-05 
    ## ... Procrustes: rmse 3.865008e-05  max resid 7.297507e-05 
    ## ... Similar to previous best
    ## Run 409 stress 9.716469e-05 
    ## ... Procrustes: rmse 0.0002400276  max resid 0.000381273 
    ## ... Similar to previous best
    ## Run 410 stress 9.756323e-05 
    ## ... Procrustes: rmse 0.0002114719  max resid 0.0003586631 
    ## ... Similar to previous best
    ## Run 411 stress 5.869124e-05 
    ## ... Procrustes: rmse 7.504489e-05  max resid 0.000118758 
    ## ... Similar to previous best
    ## Run 412 stress 9.780517e-05 
    ## ... Procrustes: rmse 0.0002136888  max resid 0.0003573747 
    ## ... Similar to previous best
    ## Run 413 stress 5.754542e-05 
    ## ... Procrustes: rmse 3.465246e-05  max resid 6.791698e-05 
    ## ... Similar to previous best
    ## Run 414 stress 9.438027e-05 
    ## ... Procrustes: rmse 0.0001725311  max resid 0.0002591216 
    ## ... Similar to previous best
    ## Run 415 stress 8.781582e-05 
    ## ... Procrustes: rmse 0.0001409563  max resid 0.0002347633 
    ## ... Similar to previous best
    ## Run 416 stress 5.481724e-05 
    ## ... Procrustes: rmse 3.49165e-05  max resid 7.915027e-05 
    ## ... Similar to previous best
    ## Run 417 stress 9.823672e-05 
    ## ... Procrustes: rmse 0.000237646  max resid 0.0003811789 
    ## ... Similar to previous best
    ## Run 418 stress 8.591527e-05 
    ## ... Procrustes: rmse 5.262614e-05  max resid 8.535883e-05 
    ## ... Similar to previous best
    ## Run 419 stress 9.502304e-05 
    ## ... Procrustes: rmse 0.0002276979  max resid 0.0003696829 
    ## ... Similar to previous best
    ## Run 420 stress 9.52364e-05 
    ## ... Procrustes: rmse 0.0002054631  max resid 0.0003405836 
    ## ... Similar to previous best
    ## Run 421 stress 9.730689e-05 
    ## ... Procrustes: rmse 0.0002125199  max resid 0.0003251822 
    ## ... Similar to previous best
    ## Run 422 stress 9.700643e-05 
    ## ... Procrustes: rmse 0.0002406555  max resid 0.0004159618 
    ## ... Similar to previous best
    ## Run 423 stress 9.6567e-05 
    ## ... Procrustes: rmse 0.0002047065  max resid 0.0003138751 
    ## ... Similar to previous best
    ## Run 424 stress 9.70147e-05 
    ## ... Procrustes: rmse 0.0001847925  max resid 0.0003136805 
    ## ... Similar to previous best
    ## Run 425 stress 9.765234e-05 
    ## ... Procrustes: rmse 0.0002384238  max resid 0.000384997 
    ## ... Similar to previous best
    ## Run 426 stress 9.750615e-05 
    ## ... Procrustes: rmse 0.0002317213  max resid 0.0003717408 
    ## ... Similar to previous best
    ## Run 427 stress 9.841111e-05 
    ## ... Procrustes: rmse 0.0002123717  max resid 0.0003369007 
    ## ... Similar to previous best
    ## Run 428 stress 8.218722e-05 
    ## ... Procrustes: rmse 0.0001123225  max resid 0.0001889446 
    ## ... Similar to previous best
    ## Run 429 stress 5.920849e-05 
    ## ... Procrustes: rmse 3.154128e-05  max resid 7.909086e-05 
    ## ... Similar to previous best
    ## Run 430 stress 9.948573e-05 
    ## ... Procrustes: rmse 6.567557e-05  max resid 0.0001413555 
    ## ... Similar to previous best
    ## Run 431 stress 9.817558e-05 
    ## ... Procrustes: rmse 0.0001805144  max resid 0.0002836029 
    ## ... Similar to previous best
    ## Run 432 stress 6.110957e-05 
    ## ... Procrustes: rmse 3.656817e-05  max resid 8.316653e-05 
    ## ... Similar to previous best
    ## Run 433 stress 9.630569e-05 
    ## ... Procrustes: rmse 0.0002064166  max resid 0.000345827 
    ## ... Similar to previous best
    ## Run 434 stress 9.859882e-05 
    ## ... Procrustes: rmse 0.0002188304  max resid 0.000345953 
    ## ... Similar to previous best
    ## Run 435 stress 9.488549e-05 
    ## ... Procrustes: rmse 0.0001409609  max resid 0.000222411 
    ## ... Similar to previous best
    ## Run 436 stress 9.98647e-05 
    ## ... Procrustes: rmse 0.0002161619  max resid 0.0003438814 
    ## ... Similar to previous best
    ## Run 437 stress 9.939822e-05 
    ## ... Procrustes: rmse 0.0002455538  max resid 0.0003900938 
    ## ... Similar to previous best
    ## Run 438 stress 9.97065e-05 
    ## ... Procrustes: rmse 0.0002496377  max resid 0.0004275293 
    ## ... Similar to previous best
    ## Run 439 stress 6.210671e-05 
    ## ... Procrustes: rmse 3.45641e-05  max resid 6.247753e-05 
    ## ... Similar to previous best
    ## Run 440 stress 9.931379e-05 
    ## ... Procrustes: rmse 0.0002085115  max resid 0.0003534015 
    ## ... Similar to previous best
    ## Run 441 stress 9.71744e-05 
    ## ... Procrustes: rmse 0.0001850475  max resid 0.000278394 
    ## ... Similar to previous best
    ## Run 442 stress 9.572023e-05 
    ## ... Procrustes: rmse 0.0002098754  max resid 0.0003461678 
    ## ... Similar to previous best
    ## Run 443 stress 9.490103e-05 
    ## ... Procrustes: rmse 0.0002376831  max resid 0.0004068626 
    ## ... Similar to previous best
    ## Run 444 stress 9.865126e-05 
    ## ... Procrustes: rmse 0.0002432126  max resid 0.000408484 
    ## ... Similar to previous best
    ## Run 445 stress 9.867731e-05 
    ## ... Procrustes: rmse 0.0001919013  max resid 0.0003130179 
    ## ... Similar to previous best
    ## Run 446 stress 9.871663e-05 
    ## ... Procrustes: rmse 0.0002430259  max resid 0.0004026385 
    ## ... Similar to previous best
    ## Run 447 stress 8.859491e-05 
    ## ... Procrustes: rmse 0.0001337649  max resid 0.0002171907 
    ## ... Similar to previous best
    ## Run 448 stress 9.909584e-05 
    ## ... Procrustes: rmse 0.0001905354  max resid 0.0002864557 
    ## ... Similar to previous best
    ## Run 449 stress 9.899908e-05 
    ## ... Procrustes: rmse 0.0001963131  max resid 0.0003005989 
    ## ... Similar to previous best
    ## Run 450 stress 8.256219e-05 
    ## ... Procrustes: rmse 4.930962e-05  max resid 9.525708e-05 
    ## ... Similar to previous best
    ## Run 451 stress 8.297121e-05 
    ## ... Procrustes: rmse 4.476962e-05  max resid 8.246858e-05 
    ## ... Similar to previous best
    ## Run 452 stress 9.897296e-05 
    ## ... Procrustes: rmse 0.000173462  max resid 0.0003054788 
    ## ... Similar to previous best
    ## Run 453 stress 8.179952e-05 
    ## ... Procrustes: rmse 7.20618e-05  max resid 0.0001071941 
    ## ... Similar to previous best
    ## Run 454 stress 9.845332e-05 
    ## ... Procrustes: rmse 0.0002130033  max resid 0.0003434054 
    ## ... Similar to previous best
    ## Run 455 stress 9.483464e-05 
    ## ... Procrustes: rmse 0.0001987023  max resid 0.0003034692 
    ## ... Similar to previous best
    ## Run 456 stress 9.791858e-05 
    ## ... Procrustes: rmse 0.0002092904  max resid 0.0003185766 
    ## ... Similar to previous best
    ## Run 457 stress 9.859351e-05 
    ## ... Procrustes: rmse 0.0002450627  max resid 0.0003894708 
    ## ... Similar to previous best
    ## Run 458 stress 8.619512e-05 
    ## ... Procrustes: rmse 4.640907e-05  max resid 0.0001052784 
    ## ... Similar to previous best
    ## Run 459 stress 9.828034e-05 
    ## ... Procrustes: rmse 0.0002010194  max resid 0.0003019418 
    ## ... Similar to previous best
    ## Run 460 stress 9.955751e-05 
    ## ... Procrustes: rmse 0.0001914599  max resid 0.0002951227 
    ## ... Similar to previous best
    ## Run 461 stress 9.647971e-05 
    ## ... Procrustes: rmse 0.0002258937  max resid 0.0003621622 
    ## ... Similar to previous best
    ## Run 462 stress 9.968548e-05 
    ## ... Procrustes: rmse 0.0002036835  max resid 0.000329648 
    ## ... Similar to previous best
    ## Run 463 stress 9.256331e-05 
    ## ... Procrustes: rmse 0.0001438689  max resid 0.0002415817 
    ## ... Similar to previous best
    ## Run 464 stress 9.880847e-05 
    ## ... Procrustes: rmse 0.0002154786  max resid 0.0003690397 
    ## ... Similar to previous best
    ## Run 465 stress 9.815386e-05 
    ## ... Procrustes: rmse 0.0002075852  max resid 0.0003207931 
    ## ... Similar to previous best
    ## Run 466 stress 9.849672e-05 
    ## ... Procrustes: rmse 0.0002439778  max resid 0.000415852 
    ## ... Similar to previous best
    ## Run 467 stress 9.958611e-05 
    ## ... Procrustes: rmse 0.0002093821  max resid 0.0003270895 
    ## ... Similar to previous best
    ## Run 468 stress 9.801008e-05 
    ## ... Procrustes: rmse 0.0001976583  max resid 0.0003012884 
    ## ... Similar to previous best
    ## Run 469 stress 9.750918e-05 
    ## ... Procrustes: rmse 0.0002391324  max resid 0.0003860915 
    ## ... Similar to previous best
    ## Run 470 stress 9.852297e-05 
    ## ... Procrustes: rmse 0.0002272874  max resid 0.0004020429 
    ## ... Similar to previous best
    ## Run 471 stress 9.644906e-05 
    ## ... Procrustes: rmse 0.0001922755  max resid 0.0002977491 
    ## ... Similar to previous best
    ## Run 472 stress 8.375287e-05 
    ## ... Procrustes: rmse 4.411572e-05  max resid 0.0001004989 
    ## ... Similar to previous best
    ## Run 473 stress 5.463379e-05 
    ## ... Procrustes: rmse 3.731262e-05  max resid 6.309656e-05 
    ## ... Similar to previous best
    ## Run 474 stress 9.83013e-05 
    ## ... Procrustes: rmse 0.0002147339  max resid 0.0003539184 
    ## ... Similar to previous best
    ## Run 475 stress 9.284276e-05 
    ## ... Procrustes: rmse 0.0002026758  max resid 0.0003304469 
    ## ... Similar to previous best
    ## Run 476 stress 9.809585e-05 
    ## ... Procrustes: rmse 0.0002380049  max resid 0.0003772306 
    ## ... Similar to previous best
    ## Run 477 stress 6.370692e-05 
    ## ... Procrustes: rmse 3.131668e-05  max resid 5.317645e-05 
    ## ... Similar to previous best
    ## Run 478 stress 9.387145e-05 
    ## ... Procrustes: rmse 0.0001975065  max resid 0.0003103508 
    ## ... Similar to previous best
    ## Run 479 stress 9.956531e-05 
    ## ... Procrustes: rmse 0.0002054602  max resid 0.0003363649 
    ## ... Similar to previous best
    ## Run 480 stress 7.457765e-05 
    ## ... Procrustes: rmse 4.086191e-05  max resid 0.000104065 
    ## ... Similar to previous best
    ## Run 481 stress 9.82634e-05 
    ## ... Procrustes: rmse 0.0002453341  max resid 0.0004186977 
    ## ... Similar to previous best
    ## Run 482 stress 9.622735e-05 
    ## ... Procrustes: rmse 0.0002086572  max resid 0.0003232388 
    ## ... Similar to previous best
    ## Run 483 stress 8.065755e-05 
    ## ... Procrustes: rmse 4.62792e-05  max resid 8.023276e-05 
    ## ... Similar to previous best
    ## Run 484 stress 6.012957e-05 
    ## ... Procrustes: rmse 4.043442e-05  max resid 8.157681e-05 
    ## ... Similar to previous best
    ## Run 485 stress 9.964946e-05 
    ## ... Procrustes: rmse 0.0002145596  max resid 0.0003602993 
    ## ... Similar to previous best
    ## Run 486 stress 9.592726e-05 
    ## ... Procrustes: rmse 0.0002087697  max resid 0.00035224 
    ## ... Similar to previous best
    ## Run 487 stress 9.985124e-05 
    ## ... Procrustes: rmse 0.0002493366  max resid 0.0004314097 
    ## ... Similar to previous best
    ## Run 488 stress 9.609566e-05 
    ## ... Procrustes: rmse 0.0002351182  max resid 0.0003833244 
    ## ... Similar to previous best
    ## Run 489 stress 9.980256e-05 
    ## ... Procrustes: rmse 0.0002307939  max resid 0.0003695068 
    ## ... Similar to previous best
    ## Run 490 stress 9.805914e-05 
    ## ... Procrustes: rmse 0.0001952707  max resid 0.0003043097 
    ## ... Similar to previous best
    ## Run 491 stress 7.824882e-05 
    ## ... Procrustes: rmse 4.813709e-05  max resid 9.649706e-05 
    ## ... Similar to previous best
    ## Run 492 stress 9.983882e-05 
    ## ... Procrustes: rmse 0.0002024813  max resid 0.0003157682 
    ## ... Similar to previous best
    ## Run 493 stress 9.985274e-05 
    ## ... Procrustes: rmse 0.000247645  max resid 0.0004146826 
    ## ... Similar to previous best
    ## Run 494 stress 7.398366e-05 
    ## ... Procrustes: rmse 8.883437e-05  max resid 0.0001446943 
    ## ... Similar to previous best
    ## Run 495 stress 9.713341e-05 
    ## ... Procrustes: rmse 0.0002154593  max resid 0.0003522129 
    ## ... Similar to previous best
    ## Run 496 stress 9.542838e-05 
    ## ... Procrustes: rmse 0.0002107258  max resid 0.0003429614 
    ## ... Similar to previous best
    ## Run 497 stress 9.890094e-05 
    ## ... Procrustes: rmse 0.0002044791  max resid 0.0003290006 
    ## ... Similar to previous best
    ## Run 498 stress 7.363897e-05 
    ## ... Procrustes: rmse 3.999769e-05  max resid 6.310565e-05 
    ## ... Similar to previous best
    ## Run 499 stress 9.771682e-05 
    ## ... Procrustes: rmse 0.0002300213  max resid 0.0003685198 
    ## ... Similar to previous best
    ## Run 500 stress 9.241315e-05 
    ## ... Procrustes: rmse 0.0001838945  max resid 0.0002751226 
    ## ... Similar to previous best
    ## *** Best solution repeated 190 times

    ## Warning in metaMDS(PD_beta_ref_dist$Btotal, try = 1000, parallel = 4, trymax =
    ## 500, : stress is (nearly) zero: you may have insufficient data

``` r
# Surveyed sites 
PD_beta_NMDS <- metaMDS(PD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07927535 
    ## Run 1 stress 0.08985621 
    ## Run 2 stress 0.07936951 
    ## ... Procrustes: rmse 0.009039891  max resid 0.0324453 
    ## Run 3 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735682  max resid 0.05649241 
    ## Run 4 stress 0.09072905 
    ## Run 5 stress 0.07951443 
    ## ... Procrustes: rmse 0.01733519  max resid 0.05635657 
    ## Run 6 stress 0.09022777 
    ## Run 7 stress 0.07969259 
    ## ... Procrustes: rmse 0.01236448  max resid 0.04964809 
    ## Run 8 stress 0.07969258 
    ## ... Procrustes: rmse 0.01236878  max resid 0.04970381 
    ## Run 9 stress 0.07969259 
    ## ... Procrustes: rmse 0.01236553  max resid 0.04966544 
    ## Run 10 stress 0.07969261 
    ## ... Procrustes: rmse 0.01233031  max resid 0.04949824 
    ## Run 11 stress 0.0796926 
    ## ... Procrustes: rmse 0.01236314  max resid 0.04963283 
    ## Run 12 stress 0.07927535 
    ## ... New best solution
    ## ... Procrustes: rmse 4.886762e-05  max resid 0.0001072178 
    ## ... Similar to previous best
    ## Run 13 stress 0.07969271 
    ## ... Procrustes: rmse 0.01243677  max resid 0.05005605 
    ## Run 14 stress 0.09388349 
    ## Run 15 stress 0.07927537 
    ## ... Procrustes: rmse 5.544001e-05  max resid 0.0001456211 
    ## ... Similar to previous best
    ## Run 16 stress 0.1135254 
    ## Run 17 stress 0.07927537 
    ## ... Procrustes: rmse 8.660146e-05  max resid 0.0001957356 
    ## ... Similar to previous best
    ## Run 18 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734784  max resid 0.05637406 
    ## Run 19 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 1.94548e-05  max resid 4.976392e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.07927534 
    ## ... Procrustes: rmse 2.435213e-05  max resid 6.717145e-05 
    ## ... Similar to previous best
    ## Run 21 stress 0.07969259 
    ## ... Procrustes: rmse 0.0123663  max resid 0.04963315 
    ## Run 22 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735864  max resid 0.05650531 
    ## Run 23 stress 0.09022782 
    ## Run 24 stress 0.09290714 
    ## Run 25 stress 0.07927537 
    ## ... Procrustes: rmse 5.913343e-05  max resid 0.0001660407 
    ## ... Similar to previous best
    ## Run 26 stress 0.07927535 
    ## ... Procrustes: rmse 2.402483e-05  max resid 6.689949e-05 
    ## ... Similar to previous best
    ## Run 27 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237888  max resid 0.04971312 
    ## Run 28 stress 0.08985609 
    ## Run 29 stress 0.09072895 
    ## Run 30 stress 0.1138787 
    ## Run 31 stress 0.09022775 
    ## Run 32 stress 0.09072897 
    ## Run 33 stress 0.09022785 
    ## Run 34 stress 0.08985605 
    ## Run 35 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173411  max resid 0.05633206 
    ## Run 36 stress 0.09247103 
    ## Run 37 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237338  max resid 0.04964145 
    ## Run 38 stress 0.07969262 
    ## ... Procrustes: rmse 0.0123983  max resid 0.04975889 
    ## Run 39 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044684  max resid 0.03246239 
    ## Run 40 stress 0.1108443 
    ## Run 41 stress 0.07927537 
    ## ... Procrustes: rmse 6.829886e-05  max resid 0.0001701144 
    ## ... Similar to previous best
    ## Run 42 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173541  max resid 0.05646486 
    ## Run 43 stress 0.0793695 
    ## ... Procrustes: rmse 0.00905016  max resid 0.03249381 
    ## Run 44 stress 0.08985603 
    ## Run 45 stress 0.08985611 
    ## Run 46 stress 0.08985609 
    ## Run 47 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049128  max resid 0.03249118 
    ## Run 48 stress 0.09172785 
    ## Run 49 stress 0.07927536 
    ## ... Procrustes: rmse 4.254047e-05  max resid 0.0001183063 
    ## ... Similar to previous best
    ## Run 50 stress 0.0796926 
    ## ... Procrustes: rmse 0.01234497  max resid 0.04951406 
    ## Run 51 stress 0.08985606 
    ## Run 52 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046365  max resid 0.0324751 
    ## Run 53 stress 0.07969259 
    ## ... Procrustes: rmse 0.0123767  max resid 0.04968588 
    ## Run 54 stress 0.1058826 
    ## Run 55 stress 0.09247116 
    ## Run 56 stress 0.09022802 
    ## Run 57 stress 0.09247114 
    ## Run 58 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904714  max resid 0.03248168 
    ## Run 59 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045078  max resid 0.03246443 
    ## Run 60 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173539  max resid 0.0563942 
    ## Run 61 stress 0.08985609 
    ## Run 62 stress 0.1061856 
    ## Run 63 stress 0.09388332 
    ## Run 64 stress 0.1045071 
    ## Run 65 stress 0.09072913 
    ## Run 66 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046518  max resid 0.03247587 
    ## Run 67 stress 0.1077585 
    ## Run 68 stress 0.09247105 
    ## Run 69 stress 0.07927534 
    ## ... Procrustes: rmse 6.56477e-06  max resid 1.748145e-05 
    ## ... Similar to previous best
    ## Run 70 stress 0.07936951 
    ## ... Procrustes: rmse 0.009041523  max resid 0.03243872 
    ## Run 71 stress 0.1138779 
    ## Run 72 stress 0.1047992 
    ## Run 73 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735437  max resid 0.05647208 
    ## Run 74 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237791  max resid 0.04969627 
    ## Run 75 stress 0.08985605 
    ## Run 76 stress 0.1047633 
    ## Run 77 stress 0.07927538 
    ## ... Procrustes: rmse 7.732133e-05  max resid 0.0001703 
    ## ... Similar to previous best
    ## Run 78 stress 0.09296222 
    ## Run 79 stress 0.09043745 
    ## Run 80 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734824  max resid 0.05639185 
    ## Run 81 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237991  max resid 0.04971182 
    ## Run 82 stress 0.09106028 
    ## Run 83 stress 0.09072913 
    ## Run 84 stress 0.09172789 
    ## Run 85 stress 0.07936951 
    ## ... Procrustes: rmse 0.009035676  max resid 0.03240301 
    ## Run 86 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734859  max resid 0.05637624 
    ## Run 87 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240317  max resid 0.04986032 
    ## Run 88 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 1.232286e-05  max resid 2.672204e-05 
    ## ... Similar to previous best
    ## Run 89 stress 0.0902278 
    ## Run 90 stress 0.0793695 
    ## ... Procrustes: rmse 0.009048951  max resid 0.03248814 
    ## Run 91 stress 0.0793695 
    ## ... Procrustes: rmse 0.00903977  max resid 0.03243113 
    ## Run 92 stress 0.07927534 
    ## ... Procrustes: rmse 8.942328e-06  max resid 1.847766e-05 
    ## ... Similar to previous best
    ## Run 93 stress 0.07951445 
    ## ... Procrustes: rmse 0.01736639  max resid 0.0565928 
    ## Run 94 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735314  max resid 0.05644327 
    ## Run 95 stress 0.07951444 
    ## ... Procrustes: rmse 0.0173714  max resid 0.05654163 
    ## Run 96 stress 0.09043736 
    ## Run 97 stress 0.08985608 
    ## Run 98 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046127  max resid 0.03247214 
    ## Run 99 stress 0.07927535 
    ## ... Procrustes: rmse 1.869507e-05  max resid 4.387114e-05 
    ## ... Similar to previous best
    ## Run 100 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733668  max resid 0.05636482 
    ## Run 101 stress 0.09247104 
    ## Run 102 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734918  max resid 0.05637632 
    ## Run 103 stress 0.09106047 
    ## Run 104 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733199  max resid 0.05628046 
    ## Run 105 stress 0.07927536 
    ## ... Procrustes: rmse 5.876343e-05  max resid 0.0001632277 
    ## ... Similar to previous best
    ## Run 106 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041486  max resid 0.03244301 
    ## Run 107 stress 0.08985612 
    ## Run 108 stress 0.09290691 
    ## Run 109 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173504  max resid 0.05641293 
    ## Run 110 stress 0.0910603 
    ## Run 111 stress 0.09296226 
    ## Run 112 stress 0.07927534 
    ## ... Procrustes: rmse 1.204569e-05  max resid 3.570724e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.07936952 
    ## ... Procrustes: rmse 0.00903469  max resid 0.03240645 
    ## Run 114 stress 0.07927535 
    ## ... Procrustes: rmse 2.630856e-05  max resid 7.122692e-05 
    ## ... Similar to previous best
    ## Run 115 stress 0.07927534 
    ## ... Procrustes: rmse 5.604495e-06  max resid 1.015081e-05 
    ## ... Similar to previous best
    ## Run 116 stress 0.07969269 
    ## ... Procrustes: rmse 0.01236724  max resid 0.04956003 
    ## Run 117 stress 0.09565357 
    ## Run 118 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173536  max resid 0.05642526 
    ## Run 119 stress 0.07936952 
    ## ... Procrustes: rmse 0.009031806  max resid 0.03237982 
    ## Run 120 stress 0.09356315 
    ## Run 121 stress 0.07951445 
    ## ... Procrustes: rmse 0.01733034  max resid 0.05627074 
    ## Run 122 stress 0.09544363 
    ## Run 123 stress 0.07927534 
    ## ... Procrustes: rmse 8.461907e-06  max resid 2.466864e-05 
    ## ... Similar to previous best
    ## Run 124 stress 0.07951444 
    ## ... Procrustes: rmse 0.01735981  max resid 0.05646745 
    ## Run 125 stress 0.09290693 
    ## Run 126 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735283  max resid 0.05641664 
    ## Run 127 stress 0.09043762 
    ## Run 128 stress 0.09565389 
    ## Run 129 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734514  max resid 0.05637429 
    ## Run 130 stress 0.07927534 
    ## ... Procrustes: rmse 1.23454e-05  max resid 3.148728e-05 
    ## ... Similar to previous best
    ## Run 131 stress 0.07927534 
    ## ... Procrustes: rmse 2.439529e-05  max resid 6.490656e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.07936954 
    ## ... Procrustes: rmse 0.009026131  max resid 0.03234597 
    ## Run 133 stress 0.109712 
    ## Run 134 stress 0.0796926 
    ## ... Procrustes: rmse 0.01236981  max resid 0.0496339 
    ## Run 135 stress 0.07927537 
    ## ... Procrustes: rmse 2.876552e-05  max resid 0.0001002921 
    ## ... Similar to previous best
    ## Run 136 stress 0.09474103 
    ## Run 137 stress 0.109491 
    ## Run 138 stress 0.09072892 
    ## Run 139 stress 0.07936953 
    ## ... Procrustes: rmse 0.009031282  max resid 0.03237847 
    ## Run 140 stress 0.07927535 
    ## ... Procrustes: rmse 3.012824e-05  max resid 7.916523e-05 
    ## ... Similar to previous best
    ## Run 141 stress 0.09388347 
    ## Run 142 stress 0.106141 
    ## Run 143 stress 0.07927534 
    ## ... Procrustes: rmse 2.75768e-05  max resid 6.004677e-05 
    ## ... Similar to previous best
    ## Run 144 stress 0.08985613 
    ## Run 145 stress 0.09388329 
    ## Run 146 stress 0.09290705 
    ## Run 147 stress 0.1088087 
    ## Run 148 stress 0.09247102 
    ## Run 149 stress 0.104507 
    ## Run 150 stress 0.08985611 
    ## Run 151 stress 0.09022781 
    ## Run 152 stress 0.09072891 
    ## Run 153 stress 0.1047108 
    ## Run 154 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240118  max resid 0.04986481 
    ## Run 155 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237506  max resid 0.04969175 
    ## Run 156 stress 0.1088751 
    ## Run 157 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237521  max resid 0.04966572 
    ## Run 158 stress 0.106642 
    ## Run 159 stress 0.09106031 
    ## Run 160 stress 0.1061411 
    ## Run 161 stress 0.07936952 
    ## ... Procrustes: rmse 0.00903087  max resid 0.03237552 
    ## Run 162 stress 0.07927534 
    ## ... Procrustes: rmse 1.958031e-05  max resid 4.676134e-05 
    ## ... Similar to previous best
    ## Run 163 stress 0.08985603 
    ## Run 164 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045464  max resid 0.0324599 
    ## Run 165 stress 0.0910603 
    ## Run 166 stress 0.07951445 
    ## ... Procrustes: rmse 0.01737828  max resid 0.05659114 
    ## Run 167 stress 0.1047992 
    ## Run 168 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734904  max resid 0.05640128 
    ## Run 169 stress 0.09296263 
    ## Run 170 stress 0.07927534 
    ## ... Procrustes: rmse 1.246128e-05  max resid 3.077114e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.07927534 
    ## ... Procrustes: rmse 6.434902e-06  max resid 1.603613e-05 
    ## ... Similar to previous best
    ## Run 172 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237352  max resid 0.04969399 
    ## Run 173 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734055  max resid 0.05634904 
    ## Run 174 stress 0.0796926 
    ## ... Procrustes: rmse 0.01240052  max resid 0.04985321 
    ## Run 175 stress 0.07927534 
    ## ... Procrustes: rmse 8.37782e-06  max resid 2.55563e-05 
    ## ... Similar to previous best
    ## Run 176 stress 0.08985614 
    ## Run 177 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237467  max resid 0.04970376 
    ## Run 178 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239444  max resid 0.0497635 
    ## Run 179 stress 0.09043737 
    ## Run 180 stress 0.09247102 
    ## Run 181 stress 0.09474084 
    ## Run 182 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735156  max resid 0.0564216 
    ## Run 183 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237455  max resid 0.04969833 
    ## Run 184 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735101  max resid 0.05641191 
    ## Run 185 stress 0.08985604 
    ## Run 186 stress 0.07969262 
    ## ... Procrustes: rmse 0.01236839  max resid 0.04959984 
    ## Run 187 stress 0.09043755 
    ## Run 188 stress 0.08985615 
    ## Run 189 stress 0.08985604 
    ## Run 190 stress 0.07951447 
    ## ... Procrustes: rmse 0.01739241  max resid 0.05667407 
    ## Run 191 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735459  max resid 0.05650509 
    ## Run 192 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735052  max resid 0.05640833 
    ## Run 193 stress 0.08985604 
    ## Run 194 stress 0.09382088 
    ## Run 195 stress 0.09296245 
    ## Run 196 stress 0.07927534 
    ## ... Procrustes: rmse 2.547355e-05  max resid 6.847333e-05 
    ## ... Similar to previous best
    ## Run 197 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735887  max resid 0.05648843 
    ## Run 198 stress 0.09356319 
    ## Run 199 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 8.0654e-06  max resid 1.564959e-05 
    ## ... Similar to previous best
    ## Run 200 stress 0.0796927 
    ## ... Procrustes: rmse 0.01235425  max resid 0.04948595 
    ## Run 201 stress 0.07927536 
    ## ... Procrustes: rmse 5.856364e-05  max resid 0.0001372201 
    ## ... Similar to previous best
    ## Run 202 stress 0.07969262 
    ## ... Procrustes: rmse 0.01239432  max resid 0.04984466 
    ## Run 203 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735787  max resid 0.05645077 
    ## Run 204 stress 0.07969261 
    ## ... Procrustes: rmse 0.01238679  max resid 0.04973305 
    ## Run 205 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238038  max resid 0.04973877 
    ## Run 206 stress 0.09247111 
    ## Run 207 stress 0.08985616 
    ## Run 208 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733002  max resid 0.05626443 
    ## Run 209 stress 0.08985605 
    ## Run 210 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735816  max resid 0.05644554 
    ## Run 211 stress 0.07951444 
    ## ... Procrustes: rmse 0.01735907  max resid 0.05653481 
    ## Run 212 stress 0.07951441 
    ## ... Procrustes: rmse 0.017357  max resid 0.05644595 
    ## Run 213 stress 0.07936951 
    ## ... Procrustes: rmse 0.009037504  max resid 0.03241132 
    ## Run 214 stress 0.07927534 
    ## ... Procrustes: rmse 4.435971e-06  max resid 1.168174e-05 
    ## ... Similar to previous best
    ## Run 215 stress 0.07927536 
    ## ... Procrustes: rmse 3.967323e-05  max resid 0.0001022102 
    ## ... Similar to previous best
    ## Run 216 stress 0.07927534 
    ## ... Procrustes: rmse 8.630688e-06  max resid 2.389301e-05 
    ## ... Similar to previous best
    ## Run 217 stress 0.07969268 
    ## ... Procrustes: rmse 0.01236899  max resid 0.04954908 
    ## Run 218 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734099  max resid 0.05633859 
    ## Run 219 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732381  max resid 0.05621542 
    ## Run 220 stress 0.09043768 
    ## Run 221 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238433  max resid 0.04974722 
    ## Run 222 stress 0.09072893 
    ## Run 223 stress 0.09290711 
    ## Run 224 stress 0.1058826 
    ## Run 225 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735377  max resid 0.05643776 
    ## Run 226 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735381  max resid 0.05640172 
    ## Run 227 stress 0.09106028 
    ## Run 228 stress 0.07936949 
    ## ... Procrustes: rmse 0.00904897  max resid 0.0324876 
    ## Run 229 stress 0.07927535 
    ## ... Procrustes: rmse 3.332809e-05  max resid 9.166367e-05 
    ## ... Similar to previous best
    ## Run 230 stress 0.07951442 
    ## ... Procrustes: rmse 0.01732343  max resid 0.05627551 
    ## Run 231 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734006  max resid 0.05632957 
    ## Run 232 stress 0.09544397 
    ## Run 233 stress 0.09072895 
    ## Run 234 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237816  max resid 0.04974124 
    ## Run 235 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732574  max resid 0.05621838 
    ## Run 236 stress 0.09344106 
    ## Run 237 stress 0.09344087 
    ## Run 238 stress 0.09290694 
    ## Run 239 stress 0.09106031 
    ## Run 240 stress 0.07927535 
    ## ... Procrustes: rmse 1.211787e-05  max resid 2.672127e-05 
    ## ... Similar to previous best
    ## Run 241 stress 0.1047993 
    ## Run 242 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047585  max resid 0.03247976 
    ## Run 243 stress 0.09072897 
    ## Run 244 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173539  max resid 0.05643441 
    ## Run 245 stress 0.3773318 
    ## Run 246 stress 0.09072899 
    ## Run 247 stress 0.07927535 
    ## ... Procrustes: rmse 2.18482e-05  max resid 5.404202e-05 
    ## ... Similar to previous best
    ## Run 248 stress 0.07969262 
    ## ... Procrustes: rmse 0.01240579  max resid 0.04988062 
    ## Run 249 stress 0.08985607 
    ## Run 250 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237602  max resid 0.04969265 
    ## Run 251 stress 0.07969265 
    ## ... Procrustes: rmse 0.0123974  max resid 0.04987312 
    ## Run 252 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238913  max resid 0.04977703 
    ## Run 253 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735289  max resid 0.05641026 
    ## Run 254 stress 0.09382091 
    ## Run 255 stress 0.09106028 
    ## Run 256 stress 0.07936951 
    ## ... Procrustes: rmse 0.009037397  max resid 0.03242032 
    ## Run 257 stress 0.08985606 
    ## Run 258 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735759  max resid 0.0564426 
    ## Run 259 stress 0.07969263 
    ## ... Procrustes: rmse 0.01237163  max resid 0.04959622 
    ## Run 260 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237201  max resid 0.04963326 
    ## Run 261 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047707  max resid 0.03248019 
    ## Run 262 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046809  max resid 0.03248584 
    ## Run 263 stress 0.07969263 
    ## ... Procrustes: rmse 0.01236372  max resid 0.0495579 
    ## Run 264 stress 0.07969258 
    ## ... Procrustes: rmse 0.0123762  max resid 0.04970104 
    ## Run 265 stress 0.07927535 
    ## ... Procrustes: rmse 3.64506e-05  max resid 0.0001020189 
    ## ... Similar to previous best
    ## Run 266 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904878  max resid 0.03248421 
    ## Run 267 stress 0.07927534 
    ## ... Procrustes: rmse 8.947531e-06  max resid 2.184679e-05 
    ## ... Similar to previous best
    ## Run 268 stress 0.08985613 
    ## Run 269 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047307  max resid 0.03247745 
    ## Run 270 stress 0.08985604 
    ## Run 271 stress 0.07927534 
    ## ... Procrustes: rmse 2.015864e-05  max resid 5.110987e-05 
    ## ... Similar to previous best
    ## Run 272 stress 0.0793695 
    ## ... Procrustes: rmse 0.009048466  max resid 0.03248636 
    ## Run 273 stress 0.09344096 
    ## Run 274 stress 0.1047992 
    ## Run 275 stress 0.07969263 
    ## ... Procrustes: rmse 0.0123756  max resid 0.04962062 
    ## Run 276 stress 0.07927538 
    ## ... Procrustes: rmse 7.049054e-05  max resid 0.0001994669 
    ## ... Similar to previous best
    ## Run 277 stress 0.07969262 
    ## ... Procrustes: rmse 0.01240706  max resid 0.04988999 
    ## Run 278 stress 0.09252224 
    ## Run 279 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237741  max resid 0.04969853 
    ## Run 280 stress 0.07927534 
    ## ... Procrustes: rmse 2.170418e-05  max resid 5.359947e-05 
    ## ... Similar to previous best
    ## Run 281 stress 0.08985617 
    ## Run 282 stress 0.1041898 
    ## Run 283 stress 0.09344089 
    ## Run 284 stress 0.08985605 
    ## Run 285 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735218  max resid 0.05640722 
    ## Run 286 stress 0.07927535 
    ## ... Procrustes: rmse 2.817414e-05  max resid 6.081872e-05 
    ## ... Similar to previous best
    ## Run 287 stress 0.09106033 
    ## Run 288 stress 0.07951445 
    ## ... Procrustes: rmse 0.01732504  max resid 0.05621935 
    ## Run 289 stress 0.07927534 
    ## ... Procrustes: rmse 1.871898e-05  max resid 4.9898e-05 
    ## ... Similar to previous best
    ## Run 290 stress 0.07927536 
    ## ... Procrustes: rmse 5.543628e-05  max resid 0.0001540511 
    ## ... Similar to previous best
    ## Run 291 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237854  max resid 0.04971958 
    ## Run 292 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237234  max resid 0.04965133 
    ## Run 293 stress 0.07951445 
    ## ... Procrustes: rmse 0.01732668  max resid 0.05622596 
    ## Run 294 stress 0.07927534 
    ## ... Procrustes: rmse 6.32798e-06  max resid 1.906116e-05 
    ## ... Similar to previous best
    ## Run 295 stress 0.07927535 
    ## ... Procrustes: rmse 3.119476e-05  max resid 6.591604e-05 
    ## ... Similar to previous best
    ## Run 296 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734541  max resid 0.05636252 
    ## Run 297 stress 0.09043764 
    ## Run 298 stress 0.09106033 
    ## Run 299 stress 0.09043736 
    ## Run 300 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239017  max resid 0.04978597 
    ## Run 301 stress 0.07927536 
    ## ... Procrustes: rmse 3.88113e-05  max resid 0.0001057343 
    ## ... Similar to previous best
    ## Run 302 stress 0.07951444 
    ## ... Procrustes: rmse 0.0173261  max resid 0.05622389 
    ## Run 303 stress 0.09072897 
    ## Run 304 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904457  max resid 0.03245908 
    ## Run 305 stress 0.0793695 
    ## ... Procrustes: rmse 0.009051368  max resid 0.03249171 
    ## Run 306 stress 0.09290692 
    ## Run 307 stress 0.09296229 
    ## Run 308 stress 0.07927536 
    ## ... Procrustes: rmse 4.668723e-05  max resid 0.0001223097 
    ## ... Similar to previous best
    ## Run 309 stress 0.09022817 
    ## Run 310 stress 0.09296231 
    ## Run 311 stress 0.1088749 
    ## Run 312 stress 0.09290716 
    ## Run 313 stress 0.09290691 
    ## Run 314 stress 0.07927535 
    ## ... Procrustes: rmse 2.354049e-05  max resid 6.62335e-05 
    ## ... Similar to previous best
    ## Run 315 stress 0.07927534 
    ## ... Procrustes: rmse 6.26904e-06  max resid 1.856961e-05 
    ## ... Similar to previous best
    ## Run 316 stress 0.09290696 
    ## Run 317 stress 0.09022786 
    ## Run 318 stress 0.07927535 
    ## ... Procrustes: rmse 1.433527e-05  max resid 3.676923e-05 
    ## ... Similar to previous best
    ## Run 319 stress 0.07936955 
    ## ... Procrustes: rmse 0.009026728  max resid 0.03234812 
    ## Run 320 stress 0.07927543 
    ## ... Procrustes: rmse 3.222423e-05  max resid 0.0001130678 
    ## ... Similar to previous best
    ## Run 321 stress 0.07927534 
    ## ... Procrustes: rmse 1.100628e-05  max resid 2.497308e-05 
    ## ... Similar to previous best
    ## Run 322 stress 0.09106038 
    ## Run 323 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043739  max resid 0.0324546 
    ## Run 324 stress 0.08985605 
    ## Run 325 stress 0.07927535 
    ## ... Procrustes: rmse 4.479827e-05  max resid 0.000106184 
    ## ... Similar to previous best
    ## Run 326 stress 0.07927535 
    ## ... Procrustes: rmse 3.69945e-05  max resid 9.091525e-05 
    ## ... Similar to previous best
    ## Run 327 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173417  max resid 0.05640093 
    ## Run 328 stress 0.07936952 
    ## ... Procrustes: rmse 0.009039165  max resid 0.03243832 
    ## Run 329 stress 0.07927534 
    ## ... Procrustes: rmse 1.27682e-05  max resid 3.255361e-05 
    ## ... Similar to previous best
    ## Run 330 stress 0.09388345 
    ## Run 331 stress 0.07927535 
    ## ... Procrustes: rmse 1.340511e-05  max resid 3.494116e-05 
    ## ... Similar to previous best
    ## Run 332 stress 0.07927535 
    ## ... Procrustes: rmse 2.028811e-05  max resid 4.366378e-05 
    ## ... Similar to previous best
    ## Run 333 stress 0.07969268 
    ## ... Procrustes: rmse 0.01237186  max resid 0.04956821 
    ## Run 334 stress 0.1047994 
    ## Run 335 stress 0.1039003 
    ## Run 336 stress 0.09172786 
    ## Run 337 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238407  max resid 0.04971383 
    ## Run 338 stress 0.1117824 
    ## Run 339 stress 0.106141 
    ## Run 340 stress 0.07951442 
    ## ... Procrustes: rmse 0.01736975  max resid 0.05651546 
    ## Run 341 stress 0.09072906 
    ## Run 342 stress 0.09072895 
    ## Run 343 stress 0.09022777 
    ## Run 344 stress 0.09072895 
    ## Run 345 stress 0.07951443 
    ## ... Procrustes: rmse 0.01737668  max resid 0.05655647 
    ## Run 346 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733807  max resid 0.05631363 
    ## Run 347 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734462  max resid 0.05635569 
    ## Run 348 stress 0.09290693 
    ## Run 349 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735729  max resid 0.05644366 
    ## Run 350 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734694  max resid 0.05633295 
    ## Run 351 stress 0.09388342 
    ## Run 352 stress 0.09544381 
    ## Run 353 stress 0.07927537 
    ## ... Procrustes: rmse 6.126192e-05  max resid 0.0001590223 
    ## ... Similar to previous best
    ## Run 354 stress 0.08985609 
    ## Run 355 stress 0.08985613 
    ## Run 356 stress 0.07969268 
    ## ... Procrustes: rmse 0.01243186  max resid 0.0500299 
    ## Run 357 stress 0.07927534 
    ## ... Procrustes: rmse 6.677791e-06  max resid 1.717685e-05 
    ## ... Similar to previous best
    ## Run 358 stress 0.09388332 
    ## Run 359 stress 0.08985603 
    ## Run 360 stress 0.0793695 
    ## ... Procrustes: rmse 0.009038913  max resid 0.03242318 
    ## Run 361 stress 0.09252213 
    ## Run 362 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733811  max resid 0.05631377 
    ## Run 363 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734588  max resid 0.05636908 
    ## Run 364 stress 0.1056806 
    ## Run 365 stress 0.07969264 
    ## ... Procrustes: rmse 0.01238684  max resid 0.04967585 
    ## Run 366 stress 0.09072908 
    ## Run 367 stress 0.09252217 
    ## Run 368 stress 0.09022794 
    ## Run 369 stress 0.09172794 
    ## Run 370 stress 0.08985607 
    ## Run 371 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735194  max resid 0.05641421 
    ## Run 372 stress 0.08985604 
    ## Run 373 stress 0.09565358 
    ## Run 374 stress 0.09474046 
    ## Run 375 stress 0.09290689 
    ## Run 376 stress 0.09290713 
    ## Run 377 stress 0.07927539 
    ## ... Procrustes: rmse 8.673908e-05  max resid 0.0001943056 
    ## ... Similar to previous best
    ## Run 378 stress 0.07927534 
    ## ... Procrustes: rmse 1.235509e-05  max resid 3.486357e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.08985607 
    ## Run 380 stress 0.0793695 
    ## ... Procrustes: rmse 0.009042213  max resid 0.03244404 
    ## Run 381 stress 0.090729 
    ## Run 382 stress 0.09106063 
    ## Run 383 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735625  max resid 0.05645638 
    ## Run 384 stress 0.09252226 
    ## Run 385 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045254  max resid 0.03246493 
    ## Run 386 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039896  max resid 0.03242895 
    ## Run 387 stress 0.08985611 
    ## Run 388 stress 0.09388326 
    ## Run 389 stress 0.0793695 
    ## ... Procrustes: rmse 0.009052842  max resid 0.0325106 
    ## Run 390 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237301  max resid 0.04963967 
    ## Run 391 stress 0.07936951 
    ## ... Procrustes: rmse 0.009036757  max resid 0.03240819 
    ## Run 392 stress 0.08985609 
    ## Run 393 stress 0.07951442 
    ## ... Procrustes: rmse 0.01736299  max resid 0.05652783 
    ## Run 394 stress 0.07927534 
    ## ... Procrustes: rmse 3.245842e-06  max resid 8.933744e-06 
    ## ... Similar to previous best
    ## Run 395 stress 0.07927538 
    ## ... Procrustes: rmse 7.514318e-05  max resid 0.0001828494 
    ## ... Similar to previous best
    ## Run 396 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173558  max resid 0.0564311 
    ## Run 397 stress 0.09252216 
    ## Run 398 stress 0.07951442 
    ## ... Procrustes: rmse 0.01734882  max resid 0.05635717 
    ## Run 399 stress 0.07951444 
    ## ... Procrustes: rmse 0.0173239  max resid 0.05621105 
    ## Run 400 stress 0.09106033 
    ## Run 401 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045967  max resid 0.0324739 
    ## Run 402 stress 0.08985605 
    ## Run 403 stress 0.104507 
    ## Run 404 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238215  max resid 0.04973571 
    ## Run 405 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734446  max resid 0.05635908 
    ## Run 406 stress 0.07951444 
    ## ... Procrustes: rmse 0.01736427  max resid 0.0564888 
    ## Run 407 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 2.646277e-06  max resid 8.262232e-06 
    ## ... Similar to previous best
    ## Run 408 stress 0.0938837 
    ## Run 409 stress 0.09296224 
    ## Run 410 stress 0.0902279 
    ## Run 411 stress 0.07927534 
    ## ... Procrustes: rmse 2.052487e-05  max resid 5.304871e-05 
    ## ... Similar to previous best
    ## Run 412 stress 0.07927534 
    ## ... Procrustes: rmse 4.670195e-06  max resid 1.297737e-05 
    ## ... Similar to previous best
    ## Run 413 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240966  max resid 0.0498875 
    ## Run 414 stress 0.07927534 
    ## ... Procrustes: rmse 4.253273e-06  max resid 1.197006e-05 
    ## ... Similar to previous best
    ## Run 415 stress 0.08985605 
    ## Run 416 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735123  max resid 0.05640641 
    ## Run 417 stress 0.09106042 
    ## Run 418 stress 0.09072895 
    ## Run 419 stress 0.09544371 
    ## Run 420 stress 0.07927534 
    ## ... Procrustes: rmse 9.588855e-06  max resid 2.325527e-05 
    ## ... Similar to previous best
    ## Run 421 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045713  max resid 0.03247275 
    ## Run 422 stress 0.07927534 
    ## ... Procrustes: rmse 3.165683e-06  max resid 6.645846e-06 
    ## ... Similar to previous best
    ## Run 423 stress 0.0795144 
    ## ... Procrustes: rmse 0.01733905  max resid 0.05632623 
    ## Run 424 stress 0.09485765 
    ## Run 425 stress 0.07927535 
    ## ... Procrustes: rmse 4.577216e-05  max resid 0.0001155184 
    ## ... Similar to previous best
    ## Run 426 stress 0.0947403 
    ## Run 427 stress 0.09106027 
    ## Run 428 stress 0.1097216 
    ## Run 429 stress 0.07951442 
    ## ... Procrustes: rmse 0.01736919  max resid 0.05651907 
    ## Run 430 stress 0.07936949 
    ## ... Procrustes: rmse 0.009042817  max resid 0.03245712 
    ## Run 431 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734778  max resid 0.05640546 
    ## Run 432 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734419  max resid 0.05636248 
    ## Run 433 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237677  max resid 0.0496995 
    ## Run 434 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239288  max resid 0.0498043 
    ## Run 435 stress 0.09072893 
    ## Run 436 stress 0.09544356 
    ## Run 437 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046492  max resid 0.03247625 
    ## Run 438 stress 0.09296227 
    ## Run 439 stress 0.07927534 
    ## ... Procrustes: rmse 7.921834e-06  max resid 1.961861e-05 
    ## ... Similar to previous best
    ## Run 440 stress 0.1058826 
    ## Run 441 stress 0.09474077 
    ## Run 442 stress 0.09290708 
    ## Run 443 stress 0.07927535 
    ## ... Procrustes: rmse 3.657254e-05  max resid 0.0001007457 
    ## ... Similar to previous best
    ## Run 444 stress 0.08985603 
    ## Run 445 stress 0.07936957 
    ## ... Procrustes: rmse 0.009020965  max resid 0.03231388 
    ## Run 446 stress 0.07927534 
    ## ... Procrustes: rmse 1.278407e-05  max resid 3.442496e-05 
    ## ... Similar to previous best
    ## Run 447 stress 0.08985605 
    ## Run 448 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237414  max resid 0.04967647 
    ## Run 449 stress 0.07969269 
    ## ... Procrustes: rmse 0.01236986  max resid 0.04956835 
    ## Run 450 stress 0.0902279 
    ## Run 451 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734022  max resid 0.05633525 
    ## Run 452 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237842  max resid 0.04969621 
    ## Run 453 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732947  max resid 0.05625899 
    ## Run 454 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735526  max resid 0.05645894 
    ## Run 455 stress 0.08985612 
    ## Run 456 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733098  max resid 0.05627291 
    ## Run 457 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239031  max resid 0.049772 
    ## Run 458 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237586  max resid 0.04969259 
    ## Run 459 stress 0.1047631 
    ## Run 460 stress 0.1039003 
    ## Run 461 stress 0.106642 
    ## Run 462 stress 0.07936951 
    ## ... Procrustes: rmse 0.009035569  max resid 0.03240495 
    ## Run 463 stress 0.104507 
    ## Run 464 stress 0.0796926 
    ## ... Procrustes: rmse 0.01240259  max resid 0.04983925 
    ## Run 465 stress 0.07927544 
    ## ... Procrustes: rmse 4.411821e-05  max resid 0.0001074208 
    ## ... Similar to previous best
    ## Run 466 stress 0.09072898 
    ## Run 467 stress 0.08985606 
    ## Run 468 stress 0.07927534 
    ## ... Procrustes: rmse 3.86099e-06  max resid 1.277571e-05 
    ## ... Similar to previous best
    ## Run 469 stress 0.0796927 
    ## ... Procrustes: rmse 0.01243868  max resid 0.05007915 
    ## Run 470 stress 0.08985613 
    ## Run 471 stress 0.0793695 
    ## ... Procrustes: rmse 0.009040203  max resid 0.0324331 
    ## Run 472 stress 0.08985604 
    ## Run 473 stress 0.09072891 
    ## Run 474 stress 0.07927537 
    ## ... Procrustes: rmse 5.313869e-05  max resid 0.0001463618 
    ## ... Similar to previous best
    ## Run 475 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734914  max resid 0.05639683 
    ## Run 476 stress 0.09043735 
    ## Run 477 stress 0.09106029 
    ## Run 478 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735701  max resid 0.05644302 
    ## Run 479 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733794  max resid 0.05631852 
    ## Run 480 stress 0.09252215 
    ## Run 481 stress 0.09290725 
    ## Run 482 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237077  max resid 0.0496869 
    ## Run 483 stress 0.09252236 
    ## Run 484 stress 0.07927534 
    ## ... Procrustes: rmse 1.153049e-06  max resid 2.436316e-06 
    ## ... Similar to previous best
    ## Run 485 stress 0.09106032 
    ## Run 486 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173482  max resid 0.0563947 
    ## Run 487 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237616  max resid 0.04970859 
    ## Run 488 stress 0.07927534 
    ## ... Procrustes: rmse 2.115425e-05  max resid 5.776624e-05 
    ## ... Similar to previous best
    ## Run 489 stress 0.0793695 
    ## ... Procrustes: rmse 0.009042151  max resid 0.03244619 
    ## Run 490 stress 0.09043756 
    ## Run 491 stress 0.3598401 
    ## Run 492 stress 0.07951445 
    ## ... Procrustes: rmse 0.01732452  max resid 0.05621229 
    ## Run 493 stress 0.07927534 
    ## ... Procrustes: rmse 2.055221e-05  max resid 5.809949e-05 
    ## ... Similar to previous best
    ## Run 494 stress 0.07936952 
    ## ... Procrustes: rmse 0.009036819  max resid 0.03241221 
    ## Run 495 stress 0.0938207 
    ## Run 496 stress 0.07927534 
    ## ... Procrustes: rmse 4.696215e-06  max resid 1.085843e-05 
    ## ... Similar to previous best
    ## Run 497 stress 0.07936954 
    ## ... Procrustes: rmse 0.00902695  max resid 0.03235101 
    ## Run 498 stress 0.08985616 
    ## Run 499 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734376  max resid 0.05636674 
    ## Run 500 stress 0.07936955 
    ## ... Procrustes: rmse 0.009026788  max resid 0.03234999 
    ## *** Best solution repeated 17 times

``` r
round(PD_beta_NMDS$stress, digits = 2)
```

    ## [1] 0.08

``` r
PD_beta_rep_NMDS <- metaMDS(PD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.3173923 
    ## Run 1 stress 0.3158987 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1883019  max resid 0.3444579 
    ## Run 2 stress 0.3265745 
    ## Run 3 stress 0.3216895 
    ## Run 4 stress 0.3283094 
    ## Run 5 stress 0.3113012 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1564313  max resid 0.3982677 
    ## Run 6 stress 0.3364185 
    ## Run 7 stress 0.3337168 
    ## Run 8 stress 0.3262254 
    ## Run 9 stress 0.3329835 
    ## Run 10 stress 0.3379119 
    ## Run 11 stress 0.3295595 
    ## Run 12 stress 0.3267556 
    ## Run 13 stress 0.3250023 
    ## Run 14 stress 0.3217957 
    ## Run 15 stress 0.3232131 
    ## Run 16 stress 0.3357178 
    ## Run 17 stress 0.315906 
    ## Run 18 stress 0.313378 
    ## Run 19 stress 0.3383201 
    ## Run 20 stress 0.327212 
    ## Run 21 stress 0.3259204 
    ## Run 22 stress 0.3325703 
    ## Run 23 stress 0.3176454 
    ## Run 24 stress 0.3331113 
    ## Run 25 stress 0.341323 
    ## Run 26 stress 0.3344258 
    ## Run 27 stress 0.3272814 
    ## Run 28 stress 0.3247217 
    ## Run 29 stress 0.3259927 
    ## Run 30 stress 0.3211777 
    ## Run 31 stress 0.3273559 
    ## Run 32 stress 0.3338879 
    ## Run 33 stress 0.3167804 
    ## Run 34 stress 0.320453 
    ## Run 35 stress 0.3191107 
    ## Run 36 stress 0.3304977 
    ## Run 37 stress 0.3365809 
    ## Run 38 stress 0.3196782 
    ## Run 39 stress 0.325884 
    ## Run 40 stress 0.3146088 
    ## Run 41 stress 0.3203884 
    ## Run 42 stress 0.3273754 
    ## Run 43 stress 0.3278487 
    ## Run 44 stress 0.3317621 
    ## Run 45 stress 0.3374606 
    ## Run 46 stress 0.316828 
    ## Run 47 stress 0.3345902 
    ## Run 48 stress 0.324214 
    ## Run 49 stress 0.3160488 
    ## Run 50 stress 0.3184898 
    ## Run 51 stress 0.3166637 
    ## Run 52 stress 0.3171067 
    ## Run 53 stress 0.3219804 
    ## Run 54 stress 0.3373114 
    ## Run 55 stress 0.3162096 
    ## Run 56 stress 0.3343177 
    ## Run 57 stress 0.3235261 
    ## Run 58 stress 0.316889 
    ## Run 59 stress 0.3284648 
    ## Run 60 stress 0.3351137 
    ## Run 61 stress 0.3216358 
    ## Run 62 stress 0.31874 
    ## Run 63 stress 0.3186555 
    ## Run 64 stress 0.3164568 
    ## Run 65 stress 0.3208775 
    ## Run 66 stress 0.3220609 
    ## Run 67 stress 0.3245318 
    ## Run 68 stress 0.3177215 
    ## Run 69 stress 0.3280998 
    ## Run 70 stress 0.3086972 
    ## ... New best solution
    ## ... Procrustes: rmse 0.09514246  max resid 0.2087958 
    ## Run 71 stress 0.3198136 
    ## Run 72 stress 0.3255264 
    ## Run 73 stress 0.3170433 
    ## Run 74 stress 0.3395217 
    ## Run 75 stress 0.3262521 
    ## Run 76 stress 0.3254104 
    ## Run 77 stress 0.3289707 
    ## Run 78 stress 0.3316622 
    ## Run 79 stress 0.3329286 
    ## Run 80 stress 0.3275361 
    ## Run 81 stress 0.3266756 
    ## Run 82 stress 0.3274209 
    ## Run 83 stress 0.3224418 
    ## Run 84 stress 0.3250993 
    ## Run 85 stress 0.3330244 
    ## Run 86 stress 0.3207402 
    ## Run 87 stress 0.3200781 
    ## Run 88 stress 0.3259014 
    ## Run 89 stress 0.3230483 
    ## Run 90 stress 0.3150063 
    ## Run 91 stress 0.3226762 
    ## Run 92 stress 0.331681 
    ## Run 93 stress 0.3357596 
    ## Run 94 stress 0.3207139 
    ## Run 95 stress 0.3385349 
    ## Run 96 stress 0.3340213 
    ## Run 97 stress 0.3292769 
    ## Run 98 stress 0.3141469 
    ## Run 99 stress 0.3347985 
    ## Run 100 stress 0.3384437 
    ## Run 101 stress 0.3193025 
    ## Run 102 stress 0.3155191 
    ## Run 103 stress 0.3232646 
    ## Run 104 stress 0.3317897 
    ## Run 105 stress 0.3345538 
    ## Run 106 stress 0.3189205 
    ## Run 107 stress 0.3191868 
    ## Run 108 stress 0.3274565 
    ## Run 109 stress 0.3169275 
    ## Run 110 stress 0.3342016 
    ## Run 111 stress 0.3308063 
    ## Run 112 stress 0.340713 
    ## Run 113 stress 0.3189766 
    ## Run 114 stress 0.3242426 
    ## Run 115 stress 0.3152905 
    ## Run 116 stress 0.3275817 
    ## Run 117 stress 0.3348656 
    ## Run 118 stress 0.3299152 
    ## Run 119 stress 0.3402493 
    ## Run 120 stress 0.3213274 
    ## Run 121 stress 0.3223166 
    ## Run 122 stress 0.3382306 
    ## Run 123 stress 0.3313362 
    ## Run 124 stress 0.333901 
    ## Run 125 stress 0.3296486 
    ## Run 126 stress 0.3376336 
    ## Run 127 stress 0.3384674 
    ## Run 128 stress 0.340329 
    ## Run 129 stress 0.3249895 
    ## Run 130 stress 0.3147571 
    ## Run 131 stress 0.3187741 
    ## Run 132 stress 0.334726 
    ## Run 133 stress 0.329028 
    ## Run 134 stress 0.3319885 
    ## Run 135 stress 0.3240638 
    ## Run 136 stress 0.3395988 
    ## Run 137 stress 0.3383321 
    ## Run 138 stress 0.3350026 
    ## Run 139 stress 0.3193118 
    ## Run 140 stress 0.3265196 
    ## Run 141 stress 0.3349474 
    ## Run 142 stress 0.3199465 
    ## Run 143 stress 0.324986 
    ## Run 144 stress 0.3259002 
    ## Run 145 stress 0.360255 
    ## Run 146 stress 0.3174142 
    ## Run 147 stress 0.3239373 
    ## Run 148 stress 0.3310988 
    ## Run 149 stress 0.3180162 
    ## Run 150 stress 0.3196204 
    ## Run 151 stress 0.3253572 
    ## Run 152 stress 0.3273864 
    ## Run 153 stress 0.3327176 
    ## Run 154 stress 0.3146727 
    ## Run 155 stress 0.3412181 
    ## Run 156 stress 0.3291083 
    ## Run 157 stress 0.319596 
    ## Run 158 stress 0.3230555 
    ## Run 159 stress 0.3166625 
    ## Run 160 stress 0.3223316 
    ## Run 161 stress 0.3262183 
    ## Run 162 stress 0.3192773 
    ## Run 163 stress 0.3398143 
    ## Run 164 stress 0.3219207 
    ## Run 165 stress 0.3145217 
    ## Run 166 stress 0.3269364 
    ## Run 167 stress 0.3116333 
    ## Run 168 stress 0.3184116 
    ## Run 169 stress 0.3264739 
    ## Run 170 stress 0.3201819 
    ## Run 171 stress 0.3419135 
    ## Run 172 stress 0.3368865 
    ## Run 173 stress 0.3417212 
    ## Run 174 stress 0.3188734 
    ## Run 175 stress 0.3270266 
    ## Run 176 stress 0.3331808 
    ## Run 177 stress 0.3280972 
    ## Run 178 stress 0.3133472 
    ## Run 179 stress 0.3265969 
    ## Run 180 stress 0.3292327 
    ## Run 181 stress 0.3205453 
    ## Run 182 stress 0.3307545 
    ## Run 183 stress 0.32244 
    ## Run 184 stress 0.3147647 
    ## Run 185 stress 0.3163024 
    ## Run 186 stress 0.3272584 
    ## Run 187 stress 0.3132577 
    ## Run 188 stress 0.3377117 
    ## Run 189 stress 0.334733 
    ## Run 190 stress 0.3215583 
    ## Run 191 stress 0.3296492 
    ## Run 192 stress 0.3373617 
    ## Run 193 stress 0.3380696 
    ## Run 194 stress 0.3349537 
    ## Run 195 stress 0.3178324 
    ## Run 196 stress 0.3416113 
    ## Run 197 stress 0.3213046 
    ## Run 198 stress 0.3360955 
    ## Run 199 stress 0.3224806 
    ## Run 200 stress 0.3215668 
    ## Run 201 stress 0.3453638 
    ## Run 202 stress 0.3227069 
    ## Run 203 stress 0.3219283 
    ## Run 204 stress 0.3267165 
    ## Run 205 stress 0.3153713 
    ## Run 206 stress 0.3272888 
    ## Run 207 stress 0.3390877 
    ## Run 208 stress 0.3299735 
    ## Run 209 stress 0.3412105 
    ## Run 210 stress 0.3116713 
    ## Run 211 stress 0.3179289 
    ## Run 212 stress 0.3193651 
    ## Run 213 stress 0.3161385 
    ## Run 214 stress 0.3248796 
    ## Run 215 stress 0.3350961 
    ## Run 216 stress 0.3244254 
    ## Run 217 stress 0.3249795 
    ## Run 218 stress 0.333392 
    ## Run 219 stress 0.3328436 
    ## Run 220 stress 0.3285358 
    ## Run 221 stress 0.3243848 
    ## Run 222 stress 0.3333468 
    ## Run 223 stress 0.3230279 
    ## Run 224 stress 0.3334916 
    ## Run 225 stress 0.3354913 
    ## Run 226 stress 0.3275555 
    ## Run 227 stress 0.3289093 
    ## Run 228 stress 0.3368321 
    ## Run 229 stress 0.338413 
    ## Run 230 stress 0.3315999 
    ## Run 231 stress 0.3268997 
    ## Run 232 stress 0.3321666 
    ## Run 233 stress 0.3150233 
    ## Run 234 stress 0.3161742 
    ## Run 235 stress 0.3207786 
    ## Run 236 stress 0.3375565 
    ## Run 237 stress 0.3274759 
    ## Run 238 stress 0.3248293 
    ## Run 239 stress 0.3209557 
    ## Run 240 stress 0.3298302 
    ## Run 241 stress 0.3204071 
    ## Run 242 stress 0.3186728 
    ## Run 243 stress 0.3218245 
    ## Run 244 stress 0.3285592 
    ## Run 245 stress 0.3140772 
    ## Run 246 stress 0.3202725 
    ## Run 247 stress 0.3201786 
    ## Run 248 stress 0.324928 
    ## Run 249 stress 0.3315806 
    ## Run 250 stress 0.3247651 
    ## Run 251 stress 0.3258724 
    ## Run 252 stress 0.3222495 
    ## Run 253 stress 0.3404032 
    ## Run 254 stress 0.3197552 
    ## Run 255 stress 0.3176554 
    ## Run 256 stress 0.3252618 
    ## Run 257 stress 0.3371874 
    ## Run 258 stress 0.324307 
    ## Run 259 stress 0.3296994 
    ## Run 260 stress 0.323143 
    ## Run 261 stress 0.3201479 
    ## Run 262 stress 0.3334114 
    ## Run 263 stress 0.3361361 
    ## Run 264 stress 0.3219637 
    ## Run 265 stress 0.3287277 
    ## Run 266 stress 0.3345241 
    ## Run 267 stress 0.330013 
    ## Run 268 stress 0.3149699 
    ## Run 269 stress 0.3254593 
    ## Run 270 stress 0.327623 
    ## Run 271 stress 0.3330582 
    ## Run 272 stress 0.3281831 
    ## Run 273 stress 0.3255956 
    ## Run 274 stress 0.3279709 
    ## Run 275 stress 0.3363469 
    ## Run 276 stress 0.3244611 
    ## Run 277 stress 0.3367335 
    ## Run 278 stress 0.3323109 
    ## Run 279 stress 0.3208506 
    ## Run 280 stress 0.3315932 
    ## Run 281 stress 0.3230515 
    ## Run 282 stress 0.3315824 
    ## Run 283 stress 0.3243865 
    ## Run 284 stress 0.3315064 
    ## Run 285 stress 0.3220512 
    ## Run 286 stress 0.335994 
    ## Run 287 stress 0.3214731 
    ## Run 288 stress 0.327202 
    ## Run 289 stress 0.3229245 
    ## Run 290 stress 0.3174726 
    ## Run 291 stress 0.3206697 
    ## Run 292 stress 0.3307317 
    ## Run 293 stress 0.3199576 
    ## Run 294 stress 0.3306099 
    ## Run 295 stress 0.3313877 
    ## Run 296 stress 0.339575 
    ## Run 297 stress 0.3320093 
    ## Run 298 stress 0.3190126 
    ## Run 299 stress 0.3280581 
    ## Run 300 stress 0.3147193 
    ## Run 301 stress 0.3179673 
    ## Run 302 stress 0.3331967 
    ## Run 303 stress 0.3357426 
    ## Run 304 stress 0.3241431 
    ## Run 305 stress 0.3396317 
    ## Run 306 stress 0.3304021 
    ## Run 307 stress 0.3213273 
    ## Run 308 stress 0.3137231 
    ## Run 309 stress 0.3328551 
    ## Run 310 stress 0.3238467 
    ## Run 311 stress 0.3169904 
    ## Run 312 stress 0.3259351 
    ## Run 313 stress 0.3177494 
    ## Run 314 stress 0.3423487 
    ## Run 315 stress 0.3276742 
    ## Run 316 stress 0.3194362 
    ## Run 317 stress 0.3221081 
    ## Run 318 stress 0.3213864 
    ## Run 319 stress 0.3385883 
    ## Run 320 stress 0.3170923 
    ## Run 321 stress 0.3372429 
    ## Run 322 stress 0.3341651 
    ## Run 323 stress 0.321625 
    ## Run 324 stress 0.3346115 
    ## Run 325 stress 0.3320163 
    ## Run 326 stress 0.3227841 
    ## Run 327 stress 0.3404236 
    ## Run 328 stress 0.3210551 
    ## Run 329 stress 0.3335274 
    ## Run 330 stress 0.3356787 
    ## Run 331 stress 0.3241145 
    ## Run 332 stress 0.3280364 
    ## Run 333 stress 0.3193322 
    ## Run 334 stress 0.3220411 
    ## Run 335 stress 0.3191782 
    ## Run 336 stress 0.3333034 
    ## Run 337 stress 0.3358201 
    ## Run 338 stress 0.3182555 
    ## Run 339 stress 0.3254414 
    ## Run 340 stress 0.3389164 
    ## Run 341 stress 0.317078 
    ## Run 342 stress 0.3192461 
    ## Run 343 stress 0.3236041 
    ## Run 344 stress 0.3301511 
    ## Run 345 stress 0.333528 
    ## Run 346 stress 0.320651 
    ## Run 347 stress 0.3229587 
    ## Run 348 stress 0.3159198 
    ## Run 349 stress 0.3411604 
    ## Run 350 stress 0.3340029 
    ## Run 351 stress 0.3187014 
    ## Run 352 stress 0.3166267 
    ## Run 353 stress 0.3239176 
    ## Run 354 stress 0.3258504 
    ## Run 355 stress 0.3235161 
    ## Run 356 stress 0.3299985 
    ## Run 357 stress 0.3182677 
    ## Run 358 stress 0.3327363 
    ## Run 359 stress 0.3368784 
    ## Run 360 stress 0.3225143 
    ## Run 361 stress 0.3178979 
    ## Run 362 stress 0.3332817 
    ## Run 363 stress 0.3300723 
    ## Run 364 stress 0.3363201 
    ## Run 365 stress 0.3254097 
    ## Run 366 stress 0.3252047 
    ## Run 367 stress 0.3581146 
    ## Run 368 stress 0.3335392 
    ## Run 369 stress 0.3196105 
    ## Run 370 stress 0.3408997 
    ## Run 371 stress 0.3183101 
    ## Run 372 stress 0.3262008 
    ## Run 373 stress 0.3210571 
    ## Run 374 stress 0.3313027 
    ## Run 375 stress 0.3343703 
    ## Run 376 stress 0.3199801 
    ## Run 377 stress 0.3215768 
    ## Run 378 stress 0.3255836 
    ## Run 379 stress 0.3180647 
    ## Run 380 stress 0.3134472 
    ## Run 381 stress 0.3134408 
    ## Run 382 stress 0.3408656 
    ## Run 383 stress 0.3274016 
    ## Run 384 stress 0.3380532 
    ## Run 385 stress 0.3335178 
    ## Run 386 stress 0.317934 
    ## Run 387 stress 0.3182519 
    ## Run 388 stress 0.3329303 
    ## Run 389 stress 0.3158156 
    ## Run 390 stress 0.3176109 
    ## Run 391 stress 0.3363087 
    ## Run 392 stress 0.3198979 
    ## Run 393 stress 0.3240492 
    ## Run 394 stress 0.3258161 
    ## Run 395 stress 0.3384592 
    ## Run 396 stress 0.3340877 
    ## Run 397 stress 0.3310713 
    ## Run 398 stress 0.3277724 
    ## Run 399 stress 0.3201892 
    ## Run 400 stress 0.332395 
    ## Run 401 stress 0.3284189 
    ## Run 402 stress 0.3377229 
    ## Run 403 stress 0.318271 
    ## Run 404 stress 0.318469 
    ## Run 405 stress 0.334677 
    ## Run 406 stress 0.3438005 
    ## Run 407 stress 0.3307154 
    ## Run 408 stress 0.3248726 
    ## Run 409 stress 0.3214951 
    ## Run 410 stress 0.3184037 
    ## Run 411 stress 0.3265265 
    ## Run 412 stress 0.3217718 
    ## Run 413 stress 0.3146164 
    ## Run 414 stress 0.3282306 
    ## Run 415 stress 0.3289489 
    ## Run 416 stress 0.3283686 
    ## Run 417 stress 0.3323642 
    ## Run 418 stress 0.3222844 
    ## Run 419 stress 0.332098 
    ## Run 420 stress 0.3185804 
    ## Run 421 stress 0.3283341 
    ## Run 422 stress 0.3232123 
    ## Run 423 stress 0.3368358 
    ## Run 424 stress 0.32459 
    ## Run 425 stress 0.3182324 
    ## Run 426 stress 0.3196871 
    ## Run 427 stress 0.3106117 
    ## Run 428 stress 0.3293726 
    ## Run 429 stress 0.3317838 
    ## Run 430 stress 0.3184043 
    ## Run 431 stress 0.3371667 
    ## Run 432 stress 0.316608 
    ## Run 433 stress 0.3153599 
    ## Run 434 stress 0.3253113 
    ## Run 435 stress 0.3191195 
    ## Run 436 stress 0.313997 
    ## Run 437 stress 0.3162017 
    ## Run 438 stress 0.3273813 
    ## Run 439 stress 0.327878 
    ## Run 440 stress 0.3277555 
    ## Run 441 stress 0.3292299 
    ## Run 442 stress 0.3186437 
    ## Run 443 stress 0.3368732 
    ## Run 444 stress 0.3348824 
    ## Run 445 stress 0.3318176 
    ## Run 446 stress 0.3320801 
    ## Run 447 stress 0.3295035 
    ## Run 448 stress 0.3253759 
    ## Run 449 stress 0.3223845 
    ## Run 450 stress 0.3198649 
    ## Run 451 stress 0.3348174 
    ## Run 452 stress 0.3235746 
    ## Run 453 stress 0.3218259 
    ## Run 454 stress 0.324535 
    ## Run 455 stress 0.3181066 
    ## Run 456 stress 0.3246625 
    ## Run 457 stress 0.312873 
    ## Run 458 stress 0.3163912 
    ## Run 459 stress 0.3227454 
    ## Run 460 stress 0.3151851 
    ## Run 461 stress 0.3247664 
    ## Run 462 stress 0.3381527 
    ## Run 463 stress 0.3197476 
    ## Run 464 stress 0.3179173 
    ## Run 465 stress 0.3186461 
    ## Run 466 stress 0.3353221 
    ## Run 467 stress 0.3171385 
    ## Run 468 stress 0.3308096 
    ## Run 469 stress 0.3364448 
    ## Run 470 stress 0.3404941 
    ## Run 471 stress 0.3341338 
    ## Run 472 stress 0.3279266 
    ## Run 473 stress 0.3365979 
    ## Run 474 stress 0.3367801 
    ## Run 475 stress 0.3348062 
    ## Run 476 stress 0.33497 
    ## Run 477 stress 0.3308027 
    ## Run 478 stress 0.3215827 
    ## Run 479 stress 0.3320776 
    ## Run 480 stress 0.324743 
    ## Run 481 stress 0.3191483 
    ## Run 482 stress 0.3213965 
    ## Run 483 stress 0.3403449 
    ## Run 484 stress 0.3293499 
    ## Run 485 stress 0.3163439 
    ## Run 486 stress 0.3154625 
    ## Run 487 stress 0.3240808 
    ## Run 488 stress 0.3299453 
    ## Run 489 stress 0.3159606 
    ## Run 490 stress 0.321134 
    ## Run 491 stress 0.3305223 
    ## Run 492 stress 0.3345737 
    ## Run 493 stress 0.3398123 
    ## Run 494 stress 0.338685 
    ## Run 495 stress 0.3176346 
    ## Run 496 stress 0.3394136 
    ## Run 497 stress 0.3343832 
    ## Run 498 stress 0.3305841 
    ## Run 499 stress 0.3232544 
    ## Run 500 stress 0.3250726 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##    499: stress ratio > sratmax
    ##      1: scale factor of the gradient < sfgrmin

``` r
PD_beta_ric_NMDS <- metaMDS(PD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.01115177 
    ## Run 1 stress 0.01115407 
    ## ... Procrustes: rmse 0.0006329619  max resid 0.001432375 
    ## ... Similar to previous best
    ## Run 2 stress 0.01115195 
    ## ... Procrustes: rmse 6.824454e-05  max resid 0.0001544799 
    ## ... Similar to previous best
    ## Run 3 stress 0.01556108 
    ## Run 4 stress 0.01686719 
    ## Run 5 stress 0.01583078 
    ## Run 6 stress 0.01117361 
    ## ... Procrustes: rmse 0.004344701  max resid 0.009791387 
    ## ... Similar to previous best
    ## Run 7 stress 0.01115173 
    ## ... New best solution
    ## ... Procrustes: rmse 2.274843e-05  max resid 5.043236e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.01542698 
    ## Run 9 stress 0.01557017 
    ## Run 10 stress 0.01528213 
    ## Run 11 stress 0.01191575 
    ## Run 12 stress 0.01549932 
    ## Run 13 stress 0.01173042 
    ## Run 14 stress 0.0152941 
    ## Run 15 stress 0.01119966 
    ## ... Procrustes: rmse 0.006022541  max resid 0.01356229 
    ## Run 16 stress 0.3808264 
    ## Run 17 stress 0.01247437 
    ## Run 18 stress 0.01658079 
    ## Run 19 stress 0.01498594 
    ## Run 20 stress 0.01529508 
    ## Run 21 stress 0.01515711 
    ## Run 22 stress 0.0111794 
    ## ... Procrustes: rmse 0.00477406  max resid 0.01075627 
    ## Run 23 stress 0.01555305 
    ## Run 24 stress 0.01574923 
    ## Run 25 stress 0.01528181 
    ## Run 26 stress 0.01553529 
    ## Run 27 stress 0.01658134 
    ## Run 28 stress 0.01533785 
    ## Run 29 stress 0.01685484 
    ## Run 30 stress 0.01493145 
    ## Run 31 stress 0.01515735 
    ## Run 32 stress 0.01543249 
    ## Run 33 stress 0.01545126 
    ## Run 34 stress 0.01527557 
    ## Run 35 stress 0.01525887 
    ## Run 36 stress 0.0111517 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001330045  max resid 0.002999163 
    ## ... Similar to previous best
    ## Run 37 stress 0.01567186 
    ## Run 38 stress 0.01543156 
    ## Run 39 stress 0.01528179 
    ## Run 40 stress 0.01536151 
    ## Run 41 stress 0.0112455 
    ## ... Procrustes: rmse 0.006828552  max resid 0.01537694 
    ## Run 42 stress 0.0154743 
    ## Run 43 stress 0.0149936 
    ## Run 44 stress 0.01637126 
    ## Run 45 stress 0.01115141 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001175259  max resid 0.002652665 
    ## ... Similar to previous best
    ## Run 46 stress 0.01493137 
    ## Run 47 stress 0.01525798 
    ## Run 48 stress 0.01641872 
    ## Run 49 stress 0.01581989 
    ## Run 50 stress 0.01498799 
    ## Run 51 stress 0.0151574 
    ## Run 52 stress 0.0155796 
    ## Run 53 stress 0.01567251 
    ## Run 54 stress 0.01115155 
    ## ... Procrustes: rmse 0.001101189  max resid 0.002482673 
    ## ... Similar to previous best
    ## Run 55 stress 0.01556153 
    ## Run 56 stress 0.01558539 
    ## Run 57 stress 0.01648714 
    ## Run 58 stress 0.01557923 
    ## Run 59 stress 0.01115156 
    ## ... Procrustes: rmse 8.202324e-05  max resid 0.0001851767 
    ## ... Similar to previous best
    ## Run 60 stress 0.01115176 
    ## ... Procrustes: rmse 0.0001691725  max resid 0.0003817367 
    ## ... Similar to previous best
    ## Run 61 stress 0.01115167 
    ## ... Procrustes: rmse 0.0001346989  max resid 0.0003038391 
    ## ... Similar to previous best
    ## Run 62 stress 0.01528178 
    ## Run 63 stress 0.01117722 
    ## ... Procrustes: rmse 0.004459917  max resid 0.01005014 
    ## Run 64 stress 0.01515705 
    ## Run 65 stress 0.01115203 
    ## ... Procrustes: rmse 0.001308633  max resid 0.002950747 
    ## ... Similar to previous best
    ## Run 66 stress 0.01128009 
    ## ... Procrustes: rmse 0.009151395  max resid 0.02057794 
    ## Run 67 stress 0.01115188 
    ## ... Procrustes: rmse 0.0002197141  max resid 0.0004966262 
    ## ... Similar to previous best
    ## Run 68 stress 0.01567204 
    ## Run 69 stress 0.01528205 
    ## Run 70 stress 0.01493152 
    ## Run 71 stress 0.01182316 
    ## Run 72 stress 0.01516452 
    ## Run 73 stress 0.01545619 
    ## Run 74 stress 0.01527542 
    ## Run 75 stress 0.0153744 
    ## Run 76 stress 0.01516099 
    ## Run 77 stress 0.01115666 
    ## ... Procrustes: rmse 0.002350854  max resid 0.005300222 
    ## ... Similar to previous best
    ## Run 78 stress 0.01115173 
    ## ... Procrustes: rmse 0.001183188  max resid 0.002667621 
    ## ... Similar to previous best
    ## Run 79 stress 0.01115159 
    ## ... Procrustes: rmse 9.668771e-05  max resid 0.0002179091 
    ## ... Similar to previous best
    ## Run 80 stress 0.01571528 
    ## Run 81 stress 0.01498603 
    ## Run 82 stress 0.01115161 
    ## ... Procrustes: rmse 0.0001062435  max resid 0.000239584 
    ## ... Similar to previous best
    ## Run 83 stress 0.01115197 
    ## ... Procrustes: rmse 0.001277887  max resid 0.002880868 
    ## ... Similar to previous best
    ## Run 84 stress 0.01509227 
    ## Run 85 stress 0.01685404 
    ## Run 86 stress 0.0155699 
    ## Run 87 stress 0.01528205 
    ## Run 88 stress 0.0155793 
    ## Run 89 stress 0.01558563 
    ## Run 90 stress 0.0155143 
    ## Run 91 stress 0.01553486 
    ## Run 92 stress 0.01515756 
    ## Run 93 stress 0.01546464 
    ## Run 94 stress 0.01546465 
    ## Run 95 stress 0.0114667 
    ## ... Procrustes: rmse 0.01259589  max resid 0.02759781 
    ## Run 96 stress 0.01544934 
    ## Run 97 stress 0.01499115 
    ## Run 98 stress 0.01515712 
    ## Run 99 stress 0.01581959 
    ## Run 100 stress 0.01563872 
    ## Run 101 stress 0.01527561 
    ## Run 102 stress 0.01533809 
    ## Run 103 stress 0.01510441 
    ## Run 104 stress 0.01498603 
    ## Run 105 stress 0.01673803 
    ## Run 106 stress 0.01498571 
    ## Run 107 stress 0.01533801 
    ## Run 108 stress 0.01499311 
    ## Run 109 stress 0.0156716 
    ## Run 110 stress 0.01116918 
    ## ... Procrustes: rmse 0.003801336  max resid 0.008567544 
    ## ... Similar to previous best
    ## Run 111 stress 0.01515937 
    ## Run 112 stress 0.01514985 
    ## Run 113 stress 0.01536154 
    ## Run 114 stress 0.01571582 
    ## Run 115 stress 0.01562181 
    ## Run 116 stress 0.01567762 
    ## Run 117 stress 0.0149859 
    ## Run 118 stress 0.01536317 
    ## Run 119 stress 0.01116504 
    ## ... Procrustes: rmse 0.003402752  max resid 0.007669258 
    ## ... Similar to previous best
    ## Run 120 stress 0.01568496 
    ## Run 121 stress 0.01514921 
    ## Run 122 stress 0.0153617 
    ## Run 123 stress 0.01542817 
    ## Run 124 stress 0.01663523 
    ## Run 125 stress 0.01116813 
    ## ... Procrustes: rmse 0.003698695  max resid 0.008335325 
    ## ... Similar to previous best
    ## Run 126 stress 0.01550027 
    ## Run 127 stress 0.01123934 
    ## ... Procrustes: rmse 0.007754941  max resid 0.01744854 
    ## Run 128 stress 0.01641893 
    ## Run 129 stress 0.0116852 
    ## Run 130 stress 0.01527926 
    ## Run 131 stress 0.01562244 
    ## Run 132 stress 0.01515707 
    ## Run 133 stress 0.01576232 
    ## Run 134 stress 0.01499427 
    ## Run 135 stress 0.0164871 
    ## Run 136 stress 0.01562215 
    ## Run 137 stress 0.01513758 
    ## Run 138 stress 0.01582004 
    ## Run 139 stress 0.01557892 
    ## Run 140 stress 0.01115655 
    ## ... Procrustes: rmse 0.00233417  max resid 0.005262366 
    ## ... Similar to previous best
    ## Run 141 stress 0.0159432 
    ## Run 142 stress 0.01515728 
    ## Run 143 stress 0.01115212 
    ## ... Procrustes: rmse 0.001340295  max resid 0.003021916 
    ## ... Similar to previous best
    ## Run 144 stress 0.01546902 
    ## Run 145 stress 0.01564598 
    ## Run 146 stress 0.01527537 
    ## Run 147 stress 0.01115183 
    ## ... Procrustes: rmse 0.000201025  max resid 0.0004542615 
    ## ... Similar to previous best
    ## Run 148 stress 0.01567816 
    ## Run 149 stress 0.01558364 
    ## Run 150 stress 0.01115195 
    ## ... Procrustes: rmse 0.0002352538  max resid 0.0005309109 
    ## ... Similar to previous best
    ## Run 151 stress 0.01146774 
    ## ... Procrustes: rmse 0.01262208  max resid 0.02765471 
    ## Run 152 stress 0.01546011 
    ## Run 153 stress 0.01558332 
    ## Run 154 stress 0.01546482 
    ## Run 155 stress 0.01557041 
    ## Run 156 stress 0.0151571 
    ## Run 157 stress 0.01498597 
    ## Run 158 stress 0.0154321 
    ## Run 159 stress 0.01150555 
    ## ... Procrustes: rmse 0.01367857  max resid 0.02971 
    ## Run 160 stress 0.01558539 
    ## Run 161 stress 0.01145673 
    ## ... Procrustes: rmse 0.01236607  max resid 0.02716877 
    ## Run 162 stress 0.01115204 
    ## ... Procrustes: rmse 0.0002722409  max resid 0.00061563 
    ## ... Similar to previous best
    ## Run 163 stress 0.01528217 
    ## Run 164 stress 0.01515703 
    ## Run 165 stress 0.015157 
    ## Run 166 stress 0.01115354 
    ## ... Procrustes: rmse 0.00174777  max resid 0.003941015 
    ## ... Similar to previous best
    ## Run 167 stress 0.01498585 
    ## Run 168 stress 0.01115172 
    ## ... Procrustes: rmse 0.001183035  max resid 0.002667415 
    ## ... Similar to previous best
    ## Run 169 stress 0.01515666 
    ## Run 170 stress 0.01515668 
    ## Run 171 stress 0.01557903 
    ## Run 172 stress 0.01556984 
    ## Run 173 stress 0.3773861 
    ## Run 174 stress 0.01493173 
    ## Run 175 stress 0.01587247 
    ## Run 176 stress 0.01544451 
    ## Run 177 stress 0.3707802 
    ## Run 178 stress 0.01528179 
    ## Run 179 stress 0.01515742 
    ## Run 180 stress 0.01515721 
    ## Run 181 stress 0.01519793 
    ## Run 182 stress 0.01515712 
    ## Run 183 stress 0.01670614 
    ## Run 184 stress 0.0156195 
    ## Run 185 stress 0.01115184 
    ## ... Procrustes: rmse 0.0002097053  max resid 0.0004737526 
    ## ... Similar to previous best
    ## Run 186 stress 0.01498574 
    ## Run 187 stress 0.01566102 
    ## Run 188 stress 0.01138495 
    ## ... Procrustes: rmse 0.01121753  max resid 0.0249823 
    ## Run 189 stress 0.01115182 
    ## ... Procrustes: rmse 0.0002005829  max resid 0.0004530994 
    ## ... Similar to previous best
    ## Run 190 stress 0.01515712 
    ## Run 191 stress 0.01549942 
    ## Run 192 stress 0.01115151 
    ## ... Procrustes: rmse 0.001073158  max resid 0.002420065 
    ## ... Similar to previous best
    ## Run 193 stress 0.01543206 
    ## Run 194 stress 0.01515723 
    ## Run 195 stress 0.01543918 
    ## Run 196 stress 0.01493159 
    ## Run 197 stress 0.01123031 
    ## ... Procrustes: rmse 0.006546364  max resid 0.01472195 
    ## Run 198 stress 0.01529489 
    ## Run 199 stress 0.0156192 
    ## Run 200 stress 0.01498572 
    ## Run 201 stress 0.01115208 
    ## ... Procrustes: rmse 0.0002946276  max resid 0.0006656433 
    ## ... Similar to previous best
    ## Run 202 stress 0.01498609 
    ## Run 203 stress 0.01115563 
    ## ... Procrustes: rmse 0.002173942  max resid 0.004900728 
    ## ... Similar to previous best
    ## Run 204 stress 0.0156708 
    ## Run 205 stress 0.01498614 
    ## Run 206 stress 0.01119627 
    ## ... Procrustes: rmse 0.005695755  max resid 0.01282828 
    ## Run 207 stress 0.01493169 
    ## Run 208 stress 0.01116259 
    ## ... Procrustes: rmse 0.003139791  max resid 0.007077062 
    ## ... Similar to previous best
    ## Run 209 stress 0.01115811 
    ## ... Procrustes: rmse 0.002573918  max resid 0.005801471 
    ## ... Similar to previous best
    ## Run 210 stress 0.01562242 
    ## Run 211 stress 0.01510794 
    ## Run 212 stress 0.01663518 
    ## Run 213 stress 0.01515702 
    ## Run 214 stress 0.01556071 
    ## Run 215 stress 0.01548148 
    ## Run 216 stress 0.01584409 
    ## Run 217 stress 0.01536177 
    ## Run 218 stress 0.01700773 
    ## Run 219 stress 0.01701596 
    ## Run 220 stress 0.01690495 
    ## Run 221 stress 0.01560346 
    ## Run 222 stress 0.01659124 
    ## Run 223 stress 0.01115195 
    ## ... Procrustes: rmse 0.001271926  max resid 0.00286825 
    ## ... Similar to previous best
    ## Run 224 stress 0.01582953 
    ## Run 225 stress 0.0111517 
    ## ... Procrustes: rmse 0.001174702  max resid 0.002648706 
    ## ... Similar to previous best
    ## Run 226 stress 0.01672241 
    ## Run 227 stress 0.01549931 
    ## Run 228 stress 0.01665908 
    ## Run 229 stress 0.0111768 
    ## ... Procrustes: rmse 0.004375138  max resid 0.009856314 
    ## ... Similar to previous best
    ## Run 230 stress 0.01543176 
    ## Run 231 stress 0.01525794 
    ## Run 232 stress 0.01529088 
    ## Run 233 stress 0.01712915 
    ## Run 234 stress 0.01528197 
    ## Run 235 stress 0.01582156 
    ## Run 236 stress 0.01124237 
    ## ... Procrustes: rmse 0.007838406  max resid 0.01763728 
    ## Run 237 stress 0.01493154 
    ## Run 238 stress 0.0111521 
    ## ... Procrustes: rmse 0.0003068312  max resid 0.0006936178 
    ## ... Similar to previous best
    ## Run 239 stress 0.01543213 
    ## Run 240 stress 0.01117534 
    ## ... Procrustes: rmse 0.004312759  max resid 0.009719436 
    ## ... Similar to previous best
    ## Run 241 stress 0.01498609 
    ## Run 242 stress 0.01544625 
    ## Run 243 stress 0.01138952 
    ## ... Procrustes: rmse 0.01129651  max resid 0.02515091 
    ## Run 244 stress 0.01499201 
    ## Run 245 stress 0.01498582 
    ## Run 246 stress 0.01546936 
    ## Run 247 stress 0.0154533 
    ## Run 248 stress 0.01549941 
    ## Run 249 stress 0.01117125 
    ## ... Procrustes: rmse 0.003572456  max resid 0.008043524 
    ## ... Similar to previous best
    ## Run 250 stress 0.01515702 
    ## Run 251 stress 0.01115195 
    ## ... Procrustes: rmse 0.001273951  max resid 0.002872047 
    ## ... Similar to previous best
    ## Run 252 stress 0.01554233 
    ## Run 253 stress 0.01549238 
    ## Run 254 stress 0.0152821 
    ## Run 255 stress 0.01542816 
    ## Run 256 stress 0.01668216 
    ## Run 257 stress 0.01551165 
    ## Run 258 stress 0.01545385 
    ## Run 259 stress 0.0154995 
    ## Run 260 stress 0.01157652 
    ## ... Procrustes: rmse 0.01659272  max resid 0.03668768 
    ## Run 261 stress 0.01526975 
    ## Run 262 stress 0.01556133 
    ## Run 263 stress 0.01649959 
    ## Run 264 stress 0.01115195 
    ## ... Procrustes: rmse 0.0002503897  max resid 0.0005659356 
    ## ... Similar to previous best
    ## Run 265 stress 0.01556956 
    ## Run 266 stress 0.01545347 
    ## Run 267 stress 0.01569013 
    ## Run 268 stress 0.01121006 
    ## ... Procrustes: rmse 0.006426558  max resid 0.01446462 
    ## Run 269 stress 0.01669203 
    ## Run 270 stress 0.01528203 
    ## Run 271 stress 0.01544813 
    ## Run 272 stress 0.01500393 
    ## Run 273 stress 0.01557902 
    ## Run 274 stress 0.0111614 
    ## ... Procrustes: rmse 0.00278151  max resid 0.006273858 
    ## ... Similar to previous best
    ## Run 275 stress 0.01558396 
    ## Run 276 stress 0.01548193 
    ## Run 277 stress 0.01115192 
    ## ... Procrustes: rmse 0.0002331422  max resid 0.0005270569 
    ## ... Similar to previous best
    ## Run 278 stress 0.01684182 
    ## Run 279 stress 0.01163983 
    ## ... Procrustes: rmse 0.0175609  max resid 0.03794067 
    ## Run 280 stress 0.01556123 
    ## Run 281 stress 0.01552937 
    ## Run 282 stress 0.01564572 
    ## Run 283 stress 0.011212 
    ## ... Procrustes: rmse 0.006527089  max resid 0.01469037 
    ## Run 284 stress 0.01499372 
    ## Run 285 stress 0.01115174 
    ## ... Procrustes: rmse 0.001185524  max resid 0.002672629 
    ## ... Similar to previous best
    ## Run 286 stress 0.01553249 
    ## Run 287 stress 0.01556102 
    ## Run 288 stress 0.01649862 
    ## Run 289 stress 0.01559322 
    ## Run 290 stress 0.01641909 
    ## Run 291 stress 0.0151566 
    ## Run 292 stress 0.01545086 
    ## Run 293 stress 0.01663625 
    ## Run 294 stress 0.01557877 
    ## Run 295 stress 0.01116646 
    ## ... Procrustes: rmse 0.003546476  max resid 0.007993337 
    ## ... Similar to previous best
    ## Run 296 stress 0.01554615 
    ## Run 297 stress 0.01543218 
    ## Run 298 stress 0.01556145 
    ## Run 299 stress 0.01571532 
    ## Run 300 stress 0.01604362 
    ## Run 301 stress 0.01554169 
    ## Run 302 stress 0.01558329 
    ## Run 303 stress 0.01658126 
    ## Run 304 stress 0.01257964 
    ## Run 305 stress 0.01115175 
    ## ... Procrustes: rmse 0.0001706936  max resid 0.0003854847 
    ## ... Similar to previous best
    ## Run 306 stress 0.01543934 
    ## Run 307 stress 0.01553014 
    ## Run 308 stress 0.0151566 
    ## Run 309 stress 0.0152909 
    ## Run 310 stress 0.01543933 
    ## Run 311 stress 0.01577388 
    ## Run 312 stress 0.01515709 
    ## Run 313 stress 0.01542694 
    ## Run 314 stress 0.0149858 
    ## Run 315 stress 0.01515717 
    ## Run 316 stress 0.01115877 
    ## ... Procrustes: rmse 0.002661191  max resid 0.005999705 
    ## ... Similar to previous best
    ## Run 317 stress 0.0111517 
    ## ... Procrustes: rmse 0.0001462004  max resid 0.000330031 
    ## ... Similar to previous best
    ## Run 318 stress 0.01498587 
    ## Run 319 stress 0.01515742 
    ## Run 320 stress 0.01493154 
    ## Run 321 stress 0.01553789 
    ## Run 322 stress 0.0154499 
    ## Run 323 stress 0.01529082 
    ## Run 324 stress 0.01498609 
    ## Run 325 stress 0.01576731 
    ## Run 326 stress 0.01115189 
    ## ... Procrustes: rmse 0.0002234563  max resid 0.0005045436 
    ## ... Similar to previous best
    ## Run 327 stress 0.0168547 
    ## Run 328 stress 0.01150183 
    ## ... Procrustes: rmse 0.01509538  max resid 0.03344144 
    ## Run 329 stress 0.01498588 
    ## Run 330 stress 0.01558339 
    ## Run 331 stress 0.01515697 
    ## Run 332 stress 0.01700753 
    ## Run 333 stress 0.0152821 
    ## Run 334 stress 0.01557878 
    ## Run 335 stress 0.0111526 
    ## ... Procrustes: rmse 0.001494315  max resid 0.003369607 
    ## ... Similar to previous best
    ## Run 336 stress 0.01533778 
    ## Run 337 stress 0.01525846 
    ## Run 338 stress 0.01689043 
    ## Run 339 stress 0.01558302 
    ## Run 340 stress 0.01582043 
    ## Run 341 stress 0.01545348 
    ## Run 342 stress 0.01555701 
    ## Run 343 stress 0.01527902 
    ## Run 344 stress 0.0158252 
    ## Run 345 stress 0.01499146 
    ## Run 346 stress 0.01115226 
    ## ... Procrustes: rmse 0.001390026  max resid 0.003134086 
    ## ... Similar to previous best
    ## Run 347 stress 0.01116803 
    ## ... Procrustes: rmse 0.003693056  max resid 0.008322829 
    ## ... Similar to previous best
    ## Run 348 stress 0.01515696 
    ## Run 349 stress 0.01574944 
    ## Run 350 stress 0.01499129 
    ## Run 351 stress 0.0112173 
    ## ... Procrustes: rmse 0.006787046  max resid 0.01527473 
    ## Run 352 stress 0.01123872 
    ## ... Procrustes: rmse 0.007730597  max resid 0.01739378 
    ## Run 353 stress 0.01137374 
    ## ... Procrustes: rmse 0.01088899  max resid 0.02434572 
    ## Run 354 stress 0.015561 
    ## Run 355 stress 0.0111518 
    ## ... Procrustes: rmse 0.000190401  max resid 0.0004300561 
    ## ... Similar to previous best
    ## Run 356 stress 0.01115238 
    ## ... Procrustes: rmse 0.00143036  max resid 0.003225037 
    ## ... Similar to previous best
    ## Run 357 stress 0.01529412 
    ## Run 358 stress 0.01547234 
    ## Run 359 stress 0.01451961 
    ## Run 360 stress 0.01552985 
    ## Run 361 stress 0.01115195 
    ## ... Procrustes: rmse 0.001277253  max resid 0.002879831 
    ## ... Similar to previous best
    ## Run 362 stress 0.01499217 
    ## Run 363 stress 0.01553804 
    ## Run 364 stress 0.01576821 
    ## Run 365 stress 0.01561993 
    ## Run 366 stress 0.01493156 
    ## Run 367 stress 0.01516107 
    ## Run 368 stress 0.01533795 
    ## Run 369 stress 0.01546701 
    ## Run 370 stress 0.01564576 
    ## Run 371 stress 0.0120077 
    ## Run 372 stress 0.01556118 
    ## Run 373 stress 0.01553403 
    ## Run 374 stress 0.01116923 
    ## ... Procrustes: rmse 0.00380576  max resid 0.008577174 
    ## ... Similar to previous best
    ## Run 375 stress 0.01545683 
    ## Run 376 stress 0.01515731 
    ## Run 377 stress 0.01498597 
    ## Run 378 stress 0.01536191 
    ## Run 379 stress 0.01536108 
    ## Run 380 stress 0.01134843 
    ## ... Procrustes: rmse 0.01078571  max resid 0.02414513 
    ## Run 381 stress 0.01119101 
    ## ... Procrustes: rmse 0.005386097  max resid 0.01213338 
    ## Run 382 stress 0.01115197 
    ## ... Procrustes: rmse 0.001283807  max resid 0.002892733 
    ## ... Similar to previous best
    ## Run 383 stress 0.01528198 
    ## Run 384 stress 0.01562201 
    ## Run 385 stress 0.01662281 
    ## Run 386 stress 0.01528194 
    ## Run 387 stress 0.01242762 
    ## Run 388 stress 0.01552927 
    ## Run 389 stress 0.01515705 
    ## Run 390 stress 0.01668192 
    ## Run 391 stress 0.0157494 
    ## Run 392 stress 0.01115555 
    ## ... Procrustes: rmse 0.002164172  max resid 0.004879095 
    ## ... Similar to previous best
    ## Run 393 stress 0.01545236 
    ## Run 394 stress 0.01127988 
    ## ... Procrustes: rmse 0.009279982  max resid 0.02083226 
    ## Run 395 stress 0.0155699 
    ## Run 396 stress 0.01536162 
    ## Run 397 stress 0.2782669 
    ## Run 398 stress 0.011226 
    ## ... Procrustes: rmse 0.007149444  max resid 0.0160875 
    ## Run 399 stress 0.01701592 
    ## Run 400 stress 0.0149911 
    ## Run 401 stress 0.01119659 
    ## ... Procrustes: rmse 0.005704457  max resid 0.01284547 
    ## Run 402 stress 0.01558431 
    ## Run 403 stress 0.01118556 
    ## ... Procrustes: rmse 0.005043374  max resid 0.01136298 
    ## Run 404 stress 0.01547203 
    ## Run 405 stress 0.0111515 
    ## ... Procrustes: rmse 0.001078877  max resid 0.002432848 
    ## ... Similar to previous best
    ## Run 406 stress 0.01575609 
    ## Run 407 stress 0.01550027 
    ## Run 408 stress 0.01648704 
    ## Run 409 stress 0.01115165 
    ## ... Procrustes: rmse 0.000116153  max resid 0.0002614374 
    ## ... Similar to previous best
    ## Run 410 stress 0.01571544 
    ## Run 411 stress 0.0111519 
    ## ... Procrustes: rmse 0.0002237769  max resid 0.000505851 
    ## ... Similar to previous best
    ## Run 412 stress 0.01515786 
    ## Run 413 stress 0.01515719 
    ## Run 414 stress 0.0167356 
    ## Run 415 stress 0.0156211 
    ## Run 416 stress 0.373834 
    ## Run 417 stress 0.01529499 
    ## Run 418 stress 0.01547397 
    ## Run 419 stress 0.01547378 
    ## Run 420 stress 0.01115192 
    ## ... Procrustes: rmse 0.0002386277  max resid 0.0005392755 
    ## ... Similar to previous best
    ## Run 421 stress 0.01562004 
    ## Run 422 stress 0.01499247 
    ## Run 423 stress 0.01536114 
    ## Run 424 stress 0.01115185 
    ## ... Procrustes: rmse 0.0001944853  max resid 0.0004385525 
    ## ... Similar to previous best
    ## Run 425 stress 0.01583118 
    ## Run 426 stress 0.01556975 
    ## Run 427 stress 0.01560082 
    ## Run 428 stress 0.01544606 
    ## Run 429 stress 0.01116472 
    ## ... Procrustes: rmse 0.003222863  max resid 0.007259844 
    ## ... Similar to previous best
    ## Run 430 stress 0.01428124 
    ## Run 431 stress 0.3093895 
    ## Run 432 stress 0.01641884 
    ## Run 433 stress 0.01115308 
    ## ... Procrustes: rmse 0.001633463  max resid 0.003683029 
    ## ... Similar to previous best
    ## Run 434 stress 0.01156814 
    ## ... Procrustes: rmse 0.0155507  max resid 0.0336152 
    ## Run 435 stress 0.2750948 
    ## Run 436 stress 0.01558344 
    ## Run 437 stress 0.0168162 
    ## Run 438 stress 0.01115812 
    ## ... Procrustes: rmse 0.002573654  max resid 0.005801131 
    ## ... Similar to previous best
    ## Run 439 stress 0.01117591 
    ## ... Procrustes: rmse 0.004329705  max resid 0.009758788 
    ## ... Similar to previous best
    ## Run 440 stress 0.0111859 
    ## ... Procrustes: rmse 0.005036882  max resid 0.01134609 
    ## Run 441 stress 0.01143747 
    ## ... Procrustes: rmse 0.01197186  max resid 0.0264006 
    ## Run 442 stress 0.01670562 
    ## Run 443 stress 0.01498597 
    ## Run 444 stress 0.01119267 
    ## ... Procrustes: rmse 0.005486061  max resid 0.01235845 
    ## Run 445 stress 0.01536192 
    ## Run 446 stress 0.0112923 
    ## ... Procrustes: rmse 0.00957337  max resid 0.02146748 
    ## Run 447 stress 0.01493144 
    ## Run 448 stress 0.01117801 
    ## ... Procrustes: rmse 0.004519378  max resid 0.01018393 
    ## Run 449 stress 0.3773686 
    ## Run 450 stress 0.01543168 
    ## Run 451 stress 0.01561946 
    ## Run 452 stress 0.01574923 
    ## Run 453 stress 0.01115184 
    ## ... Procrustes: rmse 0.0002060057  max resid 0.0004651024 
    ## ... Similar to previous best
    ## Run 454 stress 0.01360985 
    ## Run 455 stress 0.0154991 
    ## Run 456 stress 0.01515653 
    ## Run 457 stress 0.01499281 
    ## Run 458 stress 0.01493158 
    ## Run 459 stress 0.01528204 
    ## Run 460 stress 0.01562222 
    ## Run 461 stress 0.0111679 
    ## ... Procrustes: rmse 0.003656815  max resid 0.008240134 
    ## ... Similar to previous best
    ## Run 462 stress 0.01515733 
    ## Run 463 stress 0.01516735 
    ## Run 464 stress 0.0154537 
    ## Run 465 stress 0.01528212 
    ## Run 466 stress 0.01533756 
    ## Run 467 stress 0.01543221 
    ## Run 468 stress 0.01564581 
    ## Run 469 stress 0.01545523 
    ## Run 470 stress 0.01536189 
    ## Run 471 stress 0.01498593 
    ## Run 472 stress 0.01515747 
    ## Run 473 stress 0.01580045 
    ## Run 474 stress 0.01498578 
    ## Run 475 stress 0.01528148 
    ## Run 476 stress 0.01558336 
    ## Run 477 stress 0.01557901 
    ## Run 478 stress 0.01558466 
    ## Run 479 stress 0.01655648 
    ## Run 480 stress 0.01286647 
    ## Run 481 stress 0.01115203 
    ## ... Procrustes: rmse 0.0002751501  max resid 0.0006221864 
    ## ... Similar to previous best
    ## Run 482 stress 0.01130072 
    ## ... Procrustes: rmse 0.009939691  max resid 0.02235807 
    ## Run 483 stress 0.01167491 
    ## Run 484 stress 0.01115207 
    ## ... Procrustes: rmse 0.001313493  max resid 0.002961086 
    ## ... Similar to previous best
    ## Run 485 stress 0.01571574 
    ## Run 486 stress 0.01545799 
    ## Run 487 stress 0.01115176 
    ## ... Procrustes: rmse 0.0001753982  max resid 0.0003961148 
    ## ... Similar to previous best
    ## Run 488 stress 0.01583092 
    ## Run 489 stress 0.297888 
    ## Run 490 stress 0.01545267 
    ## Run 491 stress 0.01641863 
    ## Run 492 stress 0.01514966 
    ## Run 493 stress 0.01498587 
    ## Run 494 stress 0.01558353 
    ## Run 495 stress 0.0111519 
    ## ... Procrustes: rmse 0.001257352  max resid 0.002834981 
    ## ... Similar to previous best
    ## Run 496 stress 0.01515045 
    ## Run 497 stress 0.01219237 
    ## Run 498 stress 0.01122343 
    ## ... Procrustes: rmse 0.007039447  max resid 0.01584214 
    ## Run 499 stress 0.01555684 
    ## Run 500 stress 0.01498584 
    ## *** Best solution repeated 69 times

``` r
# Mixed and stratified lakes
PD_beta_MS_NMDS <- metaMDS(PD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.05036795 
    ## ... New best solution
    ## ... Procrustes: rmse 6.537218e-05  max resid 0.0001626756 
    ## ... Similar to previous best
    ## Run 2 stress 0.06231323 
    ## Run 3 stress 0.06609053 
    ## Run 4 stress 0.0540362 
    ## Run 5 stress 0.05036799 
    ## ... Procrustes: rmse 2.871045e-05  max resid 7.59588e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.05968629 
    ## Run 7 stress 0.05036797 
    ## ... Procrustes: rmse 0.0002121013  max resid 0.0005367999 
    ## ... Similar to previous best
    ## Run 8 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001270769  max resid 0.0003374637 
    ## ... Similar to previous best
    ## Run 9 stress 0.05171067 
    ## Run 10 stress 0.06938774 
    ## Run 11 stress 0.05036793 
    ## ... Procrustes: rmse 3.895708e-05  max resid 8.143654e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001244041  max resid 0.0003033615 
    ## ... Similar to previous best
    ## Run 13 stress 0.0505131 
    ## ... Procrustes: rmse 0.03569984  max resid 0.121965 
    ## Run 14 stress 0.05337578 
    ## Run 15 stress 0.0505131 
    ## ... Procrustes: rmse 0.03568025  max resid 0.1219059 
    ## Run 16 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001657627  max resid 0.0004505572 
    ## ... Similar to previous best
    ## Run 17 stress 0.05036793 
    ## ... Procrustes: rmse 4.727065e-05  max resid 8.280889e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.05036798 
    ## ... Procrustes: rmse 9.718158e-05  max resid 0.0002354673 
    ## ... Similar to previous best
    ## Run 19 stress 0.05036793 
    ## ... Procrustes: rmse 0.0001065673  max resid 0.0002722437 
    ## ... Similar to previous best
    ## Run 20 stress 0.06027432 
    ## Run 21 stress 0.3545079 
    ## Run 22 stress 0.06655711 
    ## Run 23 stress 0.05036793 
    ## ... Procrustes: rmse 0.0001000026  max resid 0.0002294896 
    ## ... Similar to previous best
    ## Run 24 stress 0.0682669 
    ## Run 25 stress 0.05051343 
    ## ... Procrustes: rmse 0.03564664  max resid 0.1216365 
    ## Run 26 stress 0.050368 
    ## ... Procrustes: rmse 0.0001068614  max resid 0.0002539223 
    ## ... Similar to previous best
    ## Run 27 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001705178  max resid 0.0004530708 
    ## ... Similar to previous best
    ## Run 28 stress 0.05171061 
    ## Run 29 stress 0.05602054 
    ## Run 30 stress 0.05928213 
    ## Run 31 stress 0.05928236 
    ## Run 32 stress 0.05968637 
    ## Run 33 stress 0.05602055 
    ## Run 34 stress 0.06478538 
    ## Run 35 stress 0.05403613 
    ## Run 36 stress 0.0517106 
    ## Run 37 stress 0.05829931 
    ## Run 38 stress 0.06231302 
    ## Run 39 stress 0.05466706 
    ## Run 40 stress 0.06844342 
    ## Run 41 stress 0.05337583 
    ## Run 42 stress 0.07204143 
    ## Run 43 stress 0.05968625 
    ## Run 44 stress 0.06362788 
    ## Run 45 stress 0.05051299 
    ## ... Procrustes: rmse 0.03566472  max resid 0.1217367 
    ## Run 46 stress 0.05466707 
    ## Run 47 stress 0.05466705 
    ## Run 48 stress 0.06640337 
    ## Run 49 stress 0.06904487 
    ## Run 50 stress 0.05036792 
    ## ... Procrustes: rmse 4.621972e-05  max resid 9.476566e-05 
    ## ... Similar to previous best
    ## Run 51 stress 0.06493688 
    ## Run 52 stress 0.05051343 
    ## ... Procrustes: rmse 0.035641  max resid 0.1216193 
    ## Run 53 stress 0.05900155 
    ## Run 54 stress 0.05829929 
    ## Run 55 stress 0.05171059 
    ## Run 56 stress 0.05403586 
    ## Run 57 stress 0.05829928 
    ## Run 58 stress 0.05403606 
    ## Run 59 stress 0.05051323 
    ## ... Procrustes: rmse 0.03571624  max resid 0.122018 
    ## Run 60 stress 0.05968626 
    ## Run 61 stress 0.05466713 
    ## Run 62 stress 0.05051346 
    ## ... Procrustes: rmse 0.03566984  max resid 0.1217059 
    ## Run 63 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001757572  max resid 0.000467711 
    ## ... Similar to previous best
    ## Run 64 stress 0.06231298 
    ## Run 65 stress 0.05337591 
    ## Run 66 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001252345  max resid 0.0003040868 
    ## ... Similar to previous best
    ## Run 67 stress 0.05036811 
    ## ... Procrustes: rmse 0.000164674  max resid 0.0003933072 
    ## ... Similar to previous best
    ## Run 68 stress 0.05051317 
    ## ... Procrustes: rmse 0.03565121  max resid 0.121677 
    ## Run 69 stress 0.06027449 
    ## Run 70 stress 0.05036794 
    ## ... Procrustes: rmse 9.86468e-05  max resid 0.000287641 
    ## ... Similar to previous best
    ## Run 71 stress 0.05466703 
    ## Run 72 stress 0.05602058 
    ## Run 73 stress 0.0560205 
    ## Run 74 stress 0.05968635 
    ## Run 75 stress 0.05036795 
    ## ... Procrustes: rmse 6.552209e-05  max resid 0.0001522364 
    ## ... Similar to previous best
    ## Run 76 stress 0.05602053 
    ## Run 77 stress 0.05900148 
    ## Run 78 stress 0.05403569 
    ## Run 79 stress 0.05829936 
    ## Run 80 stress 0.06136677 
    ## Run 81 stress 0.05829941 
    ## Run 82 stress 0.05829933 
    ## Run 83 stress 0.06950987 
    ## Run 84 stress 0.06362792 
    ## Run 85 stress 0.0688508 
    ## Run 86 stress 0.05968626 
    ## Run 87 stress 0.06896923 
    ## Run 88 stress 0.05036793 
    ## ... Procrustes: rmse 3.5391e-05  max resid 7.692839e-05 
    ## ... Similar to previous best
    ## Run 89 stress 0.05171059 
    ## Run 90 stress 0.05036792 
    ## ... Procrustes: rmse 8.786682e-05  max resid 0.0002302718 
    ## ... Similar to previous best
    ## Run 91 stress 0.05051335 
    ## ... Procrustes: rmse 0.03571568  max resid 0.1220233 
    ## Run 92 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001590593  max resid 0.0003943781 
    ## ... Similar to previous best
    ## Run 93 stress 0.05602077 
    ## Run 94 stress 0.05036797 
    ## ... Procrustes: rmse 8.768526e-05  max resid 0.0001756362 
    ## ... Similar to previous best
    ## Run 95 stress 0.05051307 
    ## ... Procrustes: rmse 0.03566517  max resid 0.1217347 
    ## Run 96 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001899334  max resid 0.000521161 
    ## ... Similar to previous best
    ## Run 97 stress 0.05171059 
    ## Run 98 stress 0.06309864 
    ## Run 99 stress 0.050368 
    ## ... Procrustes: rmse 0.0001130234  max resid 0.0002753852 
    ## ... Similar to previous best
    ## Run 100 stress 0.05051322 
    ## ... Procrustes: rmse 0.03568443  max resid 0.1219249 
    ## Run 101 stress 0.05051319 
    ## ... Procrustes: rmse 0.03568934  max resid 0.1217891 
    ## Run 102 stress 0.05036797 
    ## ... Procrustes: rmse 9.108906e-05  max resid 0.0001902122 
    ## ... Similar to previous best
    ## Run 103 stress 0.06730254 
    ## Run 104 stress 0.05051348 
    ## ... Procrustes: rmse 0.03568702  max resid 0.1219418 
    ## Run 105 stress 0.05036796 
    ## ... Procrustes: rmse 8.030916e-05  max resid 0.0001875754 
    ## ... Similar to previous best
    ## Run 106 stress 0.05337589 
    ## Run 107 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001353145  max resid 0.000332798 
    ## ... Similar to previous best
    ## Run 108 stress 0.05337583 
    ## Run 109 stress 0.05051304 
    ## ... Procrustes: rmse 0.03570775  max resid 0.12186 
    ## Run 110 stress 0.05051323 
    ## ... Procrustes: rmse 0.035642  max resid 0.1216411 
    ## Run 111 stress 0.0517106 
    ## Run 112 stress 0.05403596 
    ## Run 113 stress 0.06231289 
    ## Run 114 stress 0.050368 
    ## ... Procrustes: rmse 0.0001806066  max resid 0.0004770568 
    ## ... Similar to previous best
    ## Run 115 stress 0.05968626 
    ## Run 116 stress 0.06038941 
    ## Run 117 stress 0.05171059 
    ## Run 118 stress 0.05036796 
    ## ... Procrustes: rmse 6.068624e-05  max resid 0.0001437675 
    ## ... Similar to previous best
    ## Run 119 stress 0.05036803 
    ## ... Procrustes: rmse 0.000189138  max resid 0.0005231763 
    ## ... Similar to previous best
    ## Run 120 stress 0.0649924 
    ## Run 121 stress 0.05466705 
    ## Run 122 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001298388  max resid 0.0003657767 
    ## ... Similar to previous best
    ## Run 123 stress 0.05036792 
    ## ... Procrustes: rmse 1.546708e-05  max resid 3.193063e-05 
    ## ... Similar to previous best
    ## Run 124 stress 0.05051304 
    ## ... Procrustes: rmse 0.03573046  max resid 0.12193 
    ## Run 125 stress 0.054036 
    ## Run 126 stress 0.05829931 
    ## Run 127 stress 0.05036812 
    ## ... Procrustes: rmse 0.0001594569  max resid 0.0003939454 
    ## ... Similar to previous best
    ## Run 128 stress 0.05466713 
    ## Run 129 stress 0.05051325 
    ## ... Procrustes: rmse 0.03570407  max resid 0.121985 
    ## Run 130 stress 0.06478542 
    ## Run 131 stress 0.05171059 
    ## Run 132 stress 0.05171058 
    ## Run 133 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001321146  max resid 0.000345311 
    ## ... Similar to previous best
    ## Run 134 stress 0.05051302 
    ## ... Procrustes: rmse 0.03568938  max resid 0.1219288 
    ## Run 135 stress 0.05171063 
    ## Run 136 stress 0.05171059 
    ## Run 137 stress 0.05051347 
    ## ... Procrustes: rmse 0.03572375  max resid 0.1220513 
    ## Run 138 stress 0.06309837 
    ## Run 139 stress 0.0647855 
    ## Run 140 stress 0.0640928 
    ## Run 141 stress 0.05036794 
    ## ... Procrustes: rmse 3.129177e-05  max resid 6.323028e-05 
    ## ... Similar to previous best
    ## Run 142 stress 0.05036798 
    ## ... Procrustes: rmse 9.639399e-05  max resid 0.0002294779 
    ## ... Similar to previous best
    ## Run 143 stress 0.06476249 
    ## Run 144 stress 0.05051281 
    ## ... Procrustes: rmse 0.03567251  max resid 0.1218588 
    ## Run 145 stress 0.05171058 
    ## Run 146 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001543868  max resid 0.0004118244 
    ## ... Similar to previous best
    ## Run 147 stress 0.0533758 
    ## Run 148 stress 0.05171058 
    ## Run 149 stress 0.05171059 
    ## Run 150 stress 0.06136678 
    ## Run 151 stress 0.05171058 
    ## Run 152 stress 0.05968639 
    ## Run 153 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001010192  max resid 0.0002517122 
    ## ... Similar to previous best
    ## Run 154 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001268123  max resid 0.0003557772 
    ## ... Similar to previous best
    ## Run 155 stress 0.05337581 
    ## Run 156 stress 0.05466708 
    ## Run 157 stress 0.05036799 
    ## ... Procrustes: rmse 9.056129e-05  max resid 0.00021838 
    ## ... Similar to previous best
    ## Run 158 stress 0.05928235 
    ## Run 159 stress 0.06231306 
    ## Run 160 stress 0.05036806 
    ## ... Procrustes: rmse 0.0001478466  max resid 0.0003628513 
    ## ... Similar to previous best
    ## Run 161 stress 0.05051342 
    ## ... Procrustes: rmse 0.03569023  max resid 0.1219504 
    ## Run 162 stress 0.05036797 
    ## ... Procrustes: rmse 7.478269e-05  max resid 0.0001627913 
    ## ... Similar to previous best
    ## Run 163 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001439424  max resid 0.0004008944 
    ## ... Similar to previous best
    ## Run 164 stress 0.05403605 
    ## Run 165 stress 0.05602061 
    ## Run 166 stress 0.05036794 
    ## ... Procrustes: rmse 4.451615e-05  max resid 9.814412e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.05602057 
    ## Run 168 stress 0.05829946 
    ## Run 169 stress 0.05466703 
    ## Run 170 stress 0.0693876 
    ## Run 171 stress 0.06904479 
    ## Run 172 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 4.576145e-05  max resid 0.0001320794 
    ## ... Similar to previous best
    ## Run 173 stress 0.05829932 
    ## Run 174 stress 0.06478525 
    ## Run 175 stress 0.0517106 
    ## Run 176 stress 0.05036801 
    ## ... Procrustes: rmse 0.000155626  max resid 0.000407451 
    ## ... Similar to previous best
    ## Run 177 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001421334  max resid 0.0003694097 
    ## ... Similar to previous best
    ## Run 178 stress 0.05171058 
    ## Run 179 stress 0.0597923 
    ## Run 180 stress 0.05171059 
    ## Run 181 stress 0.05171058 
    ## Run 182 stress 0.05051305 
    ## ... Procrustes: rmse 0.03570761  max resid 0.1218704 
    ## Run 183 stress 0.06478548 
    ## Run 184 stress 0.05171059 
    ## Run 185 stress 0.05051311 
    ## ... Procrustes: rmse 0.03574061  max resid 0.1219641 
    ## Run 186 stress 0.05051312 
    ## ... Procrustes: rmse 0.0356832  max resid 0.1217883 
    ## Run 187 stress 0.05466708 
    ## Run 188 stress 0.05171061 
    ## Run 189 stress 0.05337583 
    ## Run 190 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001004352  max resid 0.0002519684 
    ## ... Similar to previous best
    ## Run 191 stress 0.05829933 
    ## Run 192 stress 0.05051328 
    ## ... Procrustes: rmse 0.03571538  max resid 0.1220279 
    ## Run 193 stress 0.0540362 
    ## Run 194 stress 0.05829934 
    ## Run 195 stress 0.05051313 
    ## ... Procrustes: rmse 0.03571727  max resid 0.1220269 
    ## Run 196 stress 0.05051319 
    ## ... Procrustes: rmse 0.03569174  max resid 0.1219536 
    ## Run 197 stress 0.07024607 
    ## Run 198 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001066999  max resid 0.0002636719 
    ## ... Similar to previous best
    ## Run 199 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.85813e-05  max resid 3.746372e-05 
    ## ... Similar to previous best
    ## Run 200 stress 0.05979227 
    ## Run 201 stress 0.05466703 
    ## Run 202 stress 0.0533759 
    ## Run 203 stress 0.05466701 
    ## Run 204 stress 0.06136691 
    ## Run 205 stress 0.05602052 
    ## Run 206 stress 0.05051369 
    ## ... Procrustes: rmse 0.03570958  max resid 0.1220198 
    ## Run 207 stress 0.05171058 
    ## Run 208 stress 0.05036791 
    ## ... Procrustes: rmse 2.152431e-05  max resid 4.569861e-05 
    ## ... Similar to previous best
    ## Run 209 stress 0.05171058 
    ## Run 210 stress 0.05036792 
    ## ... Procrustes: rmse 5.65062e-05  max resid 0.0001331803 
    ## ... Similar to previous best
    ## Run 211 stress 0.050368 
    ## ... Procrustes: rmse 0.0001367601  max resid 0.0003877416 
    ## ... Similar to previous best
    ## Run 212 stress 0.05171062 
    ## Run 213 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001766644  max resid 0.0004383282 
    ## ... Similar to previous best
    ## Run 214 stress 0.06872637 
    ## Run 215 stress 0.05602051 
    ## Run 216 stress 0.05051325 
    ## ... Procrustes: rmse 0.03573592  max resid 0.1220852 
    ## Run 217 stress 0.05337585 
    ## Run 218 stress 0.06478541 
    ## Run 219 stress 0.05829935 
    ## Run 220 stress 0.05403599 
    ## Run 221 stress 0.05051342 
    ## ... Procrustes: rmse 0.03577248  max resid 0.1222009 
    ## Run 222 stress 0.05036794 
    ## ... Procrustes: rmse 7.551667e-05  max resid 0.0001392707 
    ## ... Similar to previous best
    ## Run 223 stress 0.05051303 
    ## ... Procrustes: rmse 0.03557218  max resid 0.1215728 
    ## Run 224 stress 0.05171061 
    ## Run 225 stress 0.05051314 
    ## ... Procrustes: rmse 0.03569602  max resid 0.1219627 
    ## Run 226 stress 0.05928213 
    ## Run 227 stress 0.05900152 
    ## Run 228 stress 0.05051339 
    ## ... Procrustes: rmse 0.03572203  max resid 0.1218747 
    ## Run 229 stress 0.06309862 
    ## Run 230 stress 0.06027442 
    ## Run 231 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001593924  max resid 0.0003901844 
    ## ... Similar to previous best
    ## Run 232 stress 0.06231299 
    ## Run 233 stress 0.05036792 
    ## ... Procrustes: rmse 1.107389e-05  max resid 1.974043e-05 
    ## ... Similar to previous best
    ## Run 234 stress 0.05337583 
    ## Run 235 stress 0.05036793 
    ## ... Procrustes: rmse 7.830559e-05  max resid 0.0001611788 
    ## ... Similar to previous best
    ## Run 236 stress 0.05979233 
    ## Run 237 stress 0.06309841 
    ## Run 238 stress 0.05036804 
    ## ... Procrustes: rmse 0.000171938  max resid 0.0004125661 
    ## ... Similar to previous best
    ## Run 239 stress 0.0505132 
    ## ... Procrustes: rmse 0.03567612  max resid 0.1217548 
    ## Run 240 stress 0.05036794 
    ## ... Procrustes: rmse 4.811063e-05  max resid 0.0001095171 
    ## ... Similar to previous best
    ## Run 241 stress 0.05403604 
    ## Run 242 stress 0.05466706 
    ## Run 243 stress 0.05829931 
    ## Run 244 stress 0.05928234 
    ## Run 245 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001292955  max resid 0.0003346593 
    ## ... Similar to previous best
    ## Run 246 stress 0.05036812 
    ## ... Procrustes: rmse 0.0002155999  max resid 0.0005333559 
    ## ... Similar to previous best
    ## Run 247 stress 0.05051344 
    ## ... Procrustes: rmse 0.03571241  max resid 0.1220222 
    ## Run 248 stress 0.05036791 
    ## ... Procrustes: rmse 1.150133e-05  max resid 2.384033e-05 
    ## ... Similar to previous best
    ## Run 249 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001898346  max resid 0.0004729866 
    ## ... Similar to previous best
    ## Run 250 stress 0.06730247 
    ## Run 251 stress 0.05051293 
    ## ... Procrustes: rmse 0.03569021  max resid 0.1219317 
    ## Run 252 stress 0.06478513 
    ## Run 253 stress 0.06309841 
    ## Run 254 stress 0.06537664 
    ## Run 255 stress 0.05051332 
    ## ... Procrustes: rmse 0.03573281  max resid 0.1220803 
    ## Run 256 stress 0.05829937 
    ## Run 257 stress 0.0582993 
    ## Run 258 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001344675  max resid 0.0003855466 
    ## ... Similar to previous best
    ## Run 259 stress 0.05337579 
    ## Run 260 stress 0.06136674 
    ## Run 261 stress 0.05036792 
    ## ... Procrustes: rmse 5.879715e-05  max resid 0.0001376108 
    ## ... Similar to previous best
    ## Run 262 stress 0.05051329 
    ## ... Procrustes: rmse 0.03572948  max resid 0.1220682 
    ## Run 263 stress 0.05900148 
    ## Run 264 stress 0.05036808 
    ## ... Procrustes: rmse 0.0002030332  max resid 0.0005018838 
    ## ... Similar to previous best
    ## Run 265 stress 0.05036793 
    ## ... Procrustes: rmse 7.75689e-05  max resid 0.0001846387 
    ## ... Similar to previous best
    ## Run 266 stress 0.05171067 
    ## Run 267 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001094741  max resid 0.0002609291 
    ## ... Similar to previous best
    ## Run 268 stress 0.05036794 
    ## ... Procrustes: rmse 5.837662e-05  max resid 0.0001666563 
    ## ... Similar to previous best
    ## Run 269 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001754012  max resid 0.0004346821 
    ## ... Similar to previous best
    ## Run 270 stress 0.0582993 
    ## Run 271 stress 0.06851879 
    ## Run 272 stress 0.06111764 
    ## Run 273 stress 0.0517106 
    ## Run 274 stress 0.05337577 
    ## Run 275 stress 0.06231291 
    ## Run 276 stress 0.06362783 
    ## Run 277 stress 0.05051309 
    ## ... Procrustes: rmse 0.0357243  max resid 0.1219168 
    ## Run 278 stress 0.06231303 
    ## Run 279 stress 0.06674555 
    ## Run 280 stress 0.05968635 
    ## Run 281 stress 0.05036797 
    ## ... Procrustes: rmse 9.02058e-05  max resid 0.0002412955 
    ## ... Similar to previous best
    ## Run 282 stress 0.06504229 
    ## Run 283 stress 0.05602057 
    ## Run 284 stress 0.05051334 
    ## ... Procrustes: rmse 0.03573105  max resid 0.1220759 
    ## Run 285 stress 0.05829933 
    ## Run 286 stress 0.05829933 
    ## Run 287 stress 0.05979232 
    ## Run 288 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001159703  max resid 0.0003021921 
    ## ... Similar to previous best
    ## Run 289 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001323044  max resid 0.000319005 
    ## ... Similar to previous best
    ## Run 290 stress 0.05337586 
    ## Run 291 stress 0.06136687 
    ## Run 292 stress 0.05051338 
    ## ... Procrustes: rmse 0.03570529  max resid 0.1220003 
    ## Run 293 stress 0.05051321 
    ## ... Procrustes: rmse 0.03574792  max resid 0.1221202 
    ## Run 294 stress 0.05051307 
    ## ... Procrustes: rmse 0.03569441  max resid 0.1219542 
    ## Run 295 stress 0.06309852 
    ## Run 296 stress 0.0590016 
    ## Run 297 stress 0.05051291 
    ## ... Procrustes: rmse 0.03569121  max resid 0.1219322 
    ## Run 298 stress 0.05829936 
    ## Run 299 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001602311  max resid 0.0004404034 
    ## ... Similar to previous best
    ## Run 300 stress 0.06826707 
    ## Run 301 stress 0.05051315 
    ## ... Procrustes: rmse 0.03568015  max resid 0.1217733 
    ## Run 302 stress 0.0546671 
    ## Run 303 stress 0.05051298 
    ## ... Procrustes: rmse 0.03570077  max resid 0.1219674 
    ## Run 304 stress 0.05602061 
    ## Run 305 stress 0.05036794 
    ## ... Procrustes: rmse 8.45794e-05  max resid 0.0001849573 
    ## ... Similar to previous best
    ## Run 306 stress 0.06309852 
    ## Run 307 stress 0.05829932 
    ## Run 308 stress 0.05036803 
    ## ... Procrustes: rmse 0.000168881  max resid 0.0004135521 
    ## ... Similar to previous best
    ## Run 309 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001185553  max resid 0.0002662052 
    ## ... Similar to previous best
    ## Run 310 stress 0.0517106 
    ## Run 311 stress 0.05051325 
    ## ... Procrustes: rmse 0.03575669  max resid 0.1220031 
    ## Run 312 stress 0.05036792 
    ## ... Procrustes: rmse 4.993168e-05  max resid 0.0001286272 
    ## ... Similar to previous best
    ## Run 313 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001622753  max resid 0.0004312793 
    ## ... Similar to previous best
    ## Run 314 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001402409  max resid 0.0003847813 
    ## ... Similar to previous best
    ## Run 315 stress 0.050368 
    ## ... Procrustes: rmse 0.0001558266  max resid 0.0003828167 
    ## ... Similar to previous best
    ## Run 316 stress 0.05337584 
    ## Run 317 stress 0.05602048 
    ## Run 318 stress 0.0592823 
    ## Run 319 stress 0.06231316 
    ## Run 320 stress 0.05829936 
    ## Run 321 stress 0.05051332 
    ## ... Procrustes: rmse 0.03567127  max resid 0.1217297 
    ## Run 322 stress 0.0693878 
    ## Run 323 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001477216  max resid 0.0003748191 
    ## ... Similar to previous best
    ## Run 324 stress 0.05171058 
    ## Run 325 stress 0.05171059 
    ## Run 326 stress 0.06136677 
    ## Run 327 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001554533  max resid 0.0003715214 
    ## ... Similar to previous best
    ## Run 328 stress 0.06851893 
    ## Run 329 stress 0.05466703 
    ## Run 330 stress 0.050368 
    ## ... Procrustes: rmse 0.0001540406  max resid 0.0003776382 
    ## ... Similar to previous best
    ## Run 331 stress 0.05466702 
    ## Run 332 stress 0.05900148 
    ## Run 333 stress 0.05051282 
    ## ... Procrustes: rmse 0.03569205  max resid 0.1219267 
    ## Run 334 stress 0.0582993 
    ## Run 335 stress 0.05051298 
    ## ... Procrustes: rmse 0.03567654  max resid 0.1217814 
    ## Run 336 stress 0.06896903 
    ## Run 337 stress 0.05051305 
    ## ... Procrustes: rmse 0.03569932  max resid 0.1218417 
    ## Run 338 stress 0.06655695 
    ## Run 339 stress 0.05829931 
    ## Run 340 stress 0.06027431 
    ## Run 341 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001470999  max resid 0.0004302702 
    ## ... Similar to previous best
    ## Run 342 stress 0.0546671 
    ## Run 343 stress 0.06038952 
    ## Run 344 stress 0.0505129 
    ## ... Procrustes: rmse 0.03568258  max resid 0.1218118 
    ## Run 345 stress 0.05337584 
    ## Run 346 stress 0.06499236 
    ## Run 347 stress 0.06231301 
    ## Run 348 stress 0.06111768 
    ## Run 349 stress 0.0505135 
    ## ... Procrustes: rmse 0.03575322  max resid 0.1221468 
    ## Run 350 stress 0.05979237 
    ## Run 351 stress 0.05051301 
    ## ... Procrustes: rmse 0.03569716  max resid 0.1219583 
    ## Run 352 stress 0.05466705 
    ## Run 353 stress 0.06309866 
    ## Run 354 stress 0.05036802 
    ## ... Procrustes: rmse 0.000169289  max resid 0.0004191062 
    ## ... Similar to previous best
    ## Run 355 stress 0.05968635 
    ## Run 356 stress 0.05036795 
    ## ... Procrustes: rmse 9.191376e-05  max resid 0.0002219869 
    ## ... Similar to previous best
    ## Run 357 stress 0.05968636 
    ## Run 358 stress 0.3222158 
    ## Run 359 stress 0.05900151 
    ## Run 360 stress 0.05036792 
    ## ... Procrustes: rmse 5.667618e-05  max resid 0.0001311325 
    ## ... Similar to previous best
    ## Run 361 stress 0.05171058 
    ## Run 362 stress 0.05928214 
    ## Run 363 stress 0.05051307 
    ## ... Procrustes: rmse 0.0357361  max resid 0.1219527 
    ## Run 364 stress 0.06409296 
    ## Run 365 stress 0.05036794 
    ## ... Procrustes: rmse 7.768026e-05  max resid 0.0001889715 
    ## ... Similar to previous best
    ## Run 366 stress 0.05051352 
    ## ... Procrustes: rmse 0.03570292  max resid 0.1219951 
    ## Run 367 stress 0.05403574 
    ## Run 368 stress 0.05051336 
    ## ... Procrustes: rmse 0.03566499  max resid 0.1217062 
    ## Run 369 stress 0.05051338 
    ## ... Procrustes: rmse 0.03572249  max resid 0.1220492 
    ## Run 370 stress 0.05051302 
    ## ... Procrustes: rmse 0.03570772  max resid 0.1219842 
    ## Run 371 stress 0.06362785 
    ## Run 372 stress 0.05171063 
    ## Run 373 stress 0.05829934 
    ## Run 374 stress 0.05602058 
    ## Run 375 stress 0.06476213 
    ## Run 376 stress 0.05403585 
    ## Run 377 stress 0.06478524 
    ## Run 378 stress 0.05051328 
    ## ... Procrustes: rmse 0.03572529  max resid 0.1220562 
    ## Run 379 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001650167  max resid 0.0003949596 
    ## ... Similar to previous best
    ## Run 380 stress 0.05968634 
    ## Run 381 stress 0.05968629 
    ## Run 382 stress 0.05051346 
    ## ... Procrustes: rmse 0.03566574  max resid 0.1217058 
    ## Run 383 stress 0.05051326 
    ## ... Procrustes: rmse 0.03572235  max resid 0.1220464 
    ## Run 384 stress 0.05051296 
    ## ... Procrustes: rmse 0.03567822  max resid 0.1217903 
    ## Run 385 stress 0.05466705 
    ## Run 386 stress 0.0647854 
    ## Run 387 stress 0.05337589 
    ## Run 388 stress 0.05900147 
    ## Run 389 stress 0.05051275 
    ## ... Procrustes: rmse 0.03569141  max resid 0.1219127 
    ## Run 390 stress 0.06409288 
    ## Run 391 stress 0.05829939 
    ## Run 392 stress 0.06231313 
    ## Run 393 stress 0.05051341 
    ## ... Procrustes: rmse 0.03573723  max resid 0.1220967 
    ## Run 394 stress 0.05337581 
    ## Run 395 stress 0.05403566 
    ## Run 396 stress 0.06309848 
    ## Run 397 stress 0.0630986 
    ## Run 398 stress 0.05051285 
    ## ... Procrustes: rmse 0.03567604  max resid 0.1217993 
    ## Run 399 stress 0.06609067 
    ## Run 400 stress 0.05051297 
    ## ... Procrustes: rmse 0.03569666  max resid 0.1218427 
    ## Run 401 stress 0.05036795 
    ## ... Procrustes: rmse 8.716745e-05  max resid 0.0002303445 
    ## ... Similar to previous best
    ## Run 402 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001038922  max resid 0.0003022285 
    ## ... Similar to previous best
    ## Run 403 stress 0.05403601 
    ## Run 404 stress 0.05051293 
    ## ... Procrustes: rmse 0.03569939  max resid 0.1219582 
    ## Run 405 stress 0.05036793 
    ## ... Procrustes: rmse 7.035236e-05  max resid 0.0002056205 
    ## ... Similar to previous best
    ## Run 406 stress 0.06231292 
    ## Run 407 stress 0.06136686 
    ## Run 408 stress 0.06136685 
    ## Run 409 stress 0.06640332 
    ## Run 410 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001373828  max resid 0.0003312627 
    ## ... Similar to previous best
    ## Run 411 stress 0.05466701 
    ## Run 412 stress 0.05171059 
    ## Run 413 stress 0.05602057 
    ## Run 414 stress 0.05928239 
    ## Run 415 stress 0.05928228 
    ## Run 416 stress 0.06499234 
    ## Run 417 stress 0.05928225 
    ## Run 418 stress 0.05337592 
    ## Run 419 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001034815  max resid 0.000266373 
    ## ... Similar to previous best
    ## Run 420 stress 0.05602072 
    ## Run 421 stress 0.06027435 
    ## Run 422 stress 0.05829932 
    ## Run 423 stress 0.05036799 
    ## ... Procrustes: rmse 9.136928e-05  max resid 0.0002307043 
    ## ... Similar to previous best
    ## Run 424 stress 0.05602043 
    ## Run 425 stress 0.05036793 
    ## ... Procrustes: rmse 7.529159e-05  max resid 0.000178514 
    ## ... Similar to previous best
    ## Run 426 stress 0.06872658 
    ## Run 427 stress 0.05466701 
    ## Run 428 stress 0.0517106 
    ## Run 429 stress 0.05171059 
    ## Run 430 stress 0.06537673 
    ## Run 431 stress 0.05036794 
    ## ... Procrustes: rmse 6.704298e-05  max resid 0.0001490325 
    ## ... Similar to previous best
    ## Run 432 stress 0.050513 
    ## ... Procrustes: rmse 0.03569468  max resid 0.1219503 
    ## Run 433 stress 0.0533758 
    ## Run 434 stress 0.05602049 
    ## Run 435 stress 0.05036795 
    ## ... Procrustes: rmse 8.848769e-05  max resid 0.0001943002 
    ## ... Similar to previous best
    ## Run 436 stress 0.05036793 
    ## ... Procrustes: rmse 7.04959e-05  max resid 0.0002017255 
    ## ... Similar to previous best
    ## Run 437 stress 0.05036792 
    ## ... Procrustes: rmse 4.294216e-05  max resid 9.875944e-05 
    ## ... Similar to previous best
    ## Run 438 stress 0.05337586 
    ## Run 439 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001237476  max resid 0.0002943012 
    ## ... Similar to previous best
    ## Run 440 stress 0.06231287 
    ## Run 441 stress 0.05403605 
    ## Run 442 stress 0.05337584 
    ## Run 443 stress 0.050368 
    ## ... Procrustes: rmse 0.0001363562  max resid 0.0003725933 
    ## ... Similar to previous best
    ## Run 444 stress 0.05036794 
    ## ... Procrustes: rmse 8.675838e-05  max resid 0.0002309915 
    ## ... Similar to previous best
    ## Run 445 stress 0.05928218 
    ## Run 446 stress 0.05051297 
    ## ... Procrustes: rmse 0.03569445  max resid 0.1218366 
    ## Run 447 stress 0.05602053 
    ## Run 448 stress 0.0582994 
    ## Run 449 stress 0.06111782 
    ## Run 450 stress 0.05466711 
    ## Run 451 stress 0.05051308 
    ## ... Procrustes: rmse 0.03570336  max resid 0.1219819 
    ## Run 452 stress 0.05337577 
    ## Run 453 stress 0.06231287 
    ## Run 454 stress 0.06362786 
    ## Run 455 stress 0.05466702 
    ## Run 456 stress 0.0582993 
    ## Run 457 stress 0.06630526 
    ## Run 458 stress 0.05602051 
    ## Run 459 stress 0.05171059 
    ## Run 460 stress 0.0597923 
    ## Run 461 stress 0.06231287 
    ## Run 462 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001320201  max resid 0.0003503063 
    ## ... Similar to previous best
    ## Run 463 stress 0.05036795 
    ## ... Procrustes: rmse 9.941446e-05  max resid 0.0002373497 
    ## ... Similar to previous best
    ## Run 464 stress 0.06309843 
    ## Run 465 stress 0.05171058 
    ## Run 466 stress 0.05171059 
    ## Run 467 stress 0.050368 
    ## ... Procrustes: rmse 0.0001505161  max resid 0.0003690943 
    ## ... Similar to previous best
    ## Run 468 stress 0.05171059 
    ## Run 469 stress 0.05829928 
    ## Run 470 stress 0.05928218 
    ## Run 471 stress 0.05051346 
    ## ... Procrustes: rmse 0.03576439  max resid 0.1221783 
    ## Run 472 stress 0.05337579 
    ## Run 473 stress 0.05466708 
    ## Run 474 stress 0.050368 
    ## ... Procrustes: rmse 0.0001477441  max resid 0.0003566589 
    ## ... Similar to previous best
    ## Run 475 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001302965  max resid 0.0003821109 
    ## ... Similar to previous best
    ## Run 476 stress 0.05171059 
    ## Run 477 stress 0.06038952 
    ## Run 478 stress 0.0517106 
    ## Run 479 stress 0.05051301 
    ## ... Procrustes: rmse 0.03569241  max resid 0.1219445 
    ## Run 480 stress 0.05051296 
    ## ... Procrustes: rmse 0.0356864  max resid 0.1218142 
    ## Run 481 stress 0.05337588 
    ## Run 482 stress 0.0647852 
    ## Run 483 stress 0.05036792 
    ## ... Procrustes: rmse 3.681079e-05  max resid 9.357165e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.05928244 
    ## Run 485 stress 0.06499244 
    ## Run 486 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001637183  max resid 0.0004498518 
    ## ... Similar to previous best
    ## Run 487 stress 0.05466712 
    ## Run 488 stress 0.05171061 
    ## Run 489 stress 0.05968633 
    ## Run 490 stress 0.06038947 
    ## Run 491 stress 0.05051327 
    ## ... Procrustes: rmse 0.03572474  max resid 0.122054 
    ## Run 492 stress 0.05403601 
    ## Run 493 stress 0.06478514 
    ## Run 494 stress 0.05051304 
    ## ... Procrustes: rmse 0.03567676  max resid 0.1217752 
    ## Run 495 stress 0.05602049 
    ## Run 496 stress 0.06674562 
    ## Run 497 stress 0.05403584 
    ## Run 498 stress 0.06027449 
    ## Run 499 stress 0.06409281 
    ## Run 500 stress 0.05036793 
    ## ... Procrustes: rmse 7.848217e-05  max resid 0.0001860766 
    ## ... Similar to previous best
    ## *** Best solution repeated 64 times

``` r
# Ocean sites and mixed lakes
PD_beta_OM_NMDS <- metaMDS(PD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1016168 
    ## Run 1 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0006446856  max resid 0.001807324 
    ## ... Similar to previous best
    ## Run 2 stress 0.1530425 
    ## Run 3 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004996837  max resid 0.001402316 
    ## ... Similar to previous best
    ## Run 4 stress 0.1341868 
    ## Run 5 stress 0.1519033 
    ## Run 6 stress 0.1016167 
    ## ... Procrustes: rmse 0.0005479064  max resid 0.001539412 
    ## ... Similar to previous best
    ## Run 7 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 5.230149e-05  max resid 0.0001475586 
    ## ... Similar to previous best
    ## Run 8 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 2.219296e-05  max resid 6.254232e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.1341868 
    ## Run 10 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004398265  max resid 0.001229981 
    ## ... Similar to previous best
    ## Run 11 stress 0.1356152 
    ## Run 12 stress 0.1016168 
    ## ... Procrustes: rmse 0.000541518  max resid 0.001515007 
    ## ... Similar to previous best
    ## Run 13 stress 0.1356144 
    ## Run 14 stress 0.1016164 
    ## ... Procrustes: rmse 3.411846e-05  max resid 9.579165e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.1016164 
    ## ... Procrustes: rmse 6.731548e-05  max resid 0.0001898105 
    ## ... Similar to previous best
    ## Run 16 stress 0.1341868 
    ## Run 17 stress 0.101617 
    ## ... Procrustes: rmse 0.0006642187  max resid 0.001856857 
    ## ... Similar to previous best
    ## Run 18 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 8.160824e-06  max resid 2.28068e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.1356151 
    ## Run 20 stress 0.1356147 
    ## Run 21 stress 0.1356151 
    ## Run 22 stress 0.1016164 
    ## ... Procrustes: rmse 5.004625e-05  max resid 0.0001411772 
    ## ... Similar to previous best
    ## Run 23 stress 0.1341868 
    ## Run 24 stress 0.1016164 
    ## ... Procrustes: rmse 1.124209e-05  max resid 3.116807e-05 
    ## ... Similar to previous best
    ## Run 25 stress 0.1356149 
    ## Run 26 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006103959  max resid 0.00170556 
    ## ... Similar to previous best
    ## Run 27 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005974254  max resid 0.001668632 
    ## ... Similar to previous best
    ## Run 28 stress 0.1530422 
    ## Run 29 stress 0.1016164 
    ## ... Procrustes: rmse 5.63542e-06  max resid 1.410872e-05 
    ## ... Similar to previous best
    ## Run 30 stress 0.1016164 
    ## ... Procrustes: rmse 8.557969e-05  max resid 0.0002416346 
    ## ... Similar to previous best
    ## Run 31 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002808174  max resid 0.0007878125 
    ## ... Similar to previous best
    ## Run 32 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002125343  max resid 0.0005969564 
    ## ... Similar to previous best
    ## Run 33 stress 0.1016164 
    ## ... Procrustes: rmse 1.946886e-06  max resid 4.613133e-06 
    ## ... Similar to previous best
    ## Run 34 stress 0.1016165 
    ## ... Procrustes: rmse 0.000173769  max resid 0.0004887451 
    ## ... Similar to previous best
    ## Run 35 stress 0.1519033 
    ## Run 36 stress 0.1016164 
    ## ... Procrustes: rmse 4.931168e-06  max resid 1.395564e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 4.012225e-06  max resid 1.090772e-05 
    ## ... Similar to previous best
    ## Run 38 stress 0.1341868 
    ## Run 39 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005848497  max resid 0.001633748 
    ## ... Similar to previous best
    ## Run 40 stress 0.1016164 
    ## ... Procrustes: rmse 5.001219e-05  max resid 0.0001410958 
    ## ... Similar to previous best
    ## Run 41 stress 0.1356144 
    ## Run 42 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006336312  max resid 0.00177088 
    ## ... Similar to previous best
    ## Run 43 stress 0.1016164 
    ## ... Procrustes: rmse 2.438167e-06  max resid 6.778e-06 
    ## ... Similar to previous best
    ## Run 44 stress 0.1016164 
    ## ... Procrustes: rmse 1.598495e-05  max resid 4.495766e-05 
    ## ... Similar to previous best
    ## Run 45 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001843504  max resid 0.0005141305 
    ## ... Similar to previous best
    ## Run 46 stress 0.1530426 
    ## Run 47 stress 0.1341868 
    ## Run 48 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 1.557087e-06  max resid 3.944489e-06 
    ## ... Similar to previous best
    ## Run 49 stress 0.1341868 
    ## Run 50 stress 0.1016164 
    ## ... Procrustes: rmse 5.400439e-05  max resid 0.0001507116 
    ## ... Similar to previous best
    ## Run 51 stress 0.1356146 
    ## Run 52 stress 0.135615 
    ## Run 53 stress 0.1016164 
    ## ... Procrustes: rmse 6.849507e-05  max resid 0.0001887959 
    ## ... Similar to previous best
    ## Run 54 stress 0.1537555 
    ## Run 55 stress 0.1356146 
    ## Run 56 stress 0.1356151 
    ## Run 57 stress 0.1341868 
    ## Run 58 stress 0.1530429 
    ## Run 59 stress 0.1537555 
    ## Run 60 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003538047  max resid 0.0009911712 
    ## ... Similar to previous best
    ## Run 61 stress 0.1341868 
    ## Run 62 stress 0.1356146 
    ## Run 63 stress 0.1016164 
    ## ... Procrustes: rmse 8.713835e-06  max resid 2.366755e-05 
    ## ... Similar to previous best
    ## Run 64 stress 0.1356152 
    ## Run 65 stress 0.1356149 
    ## Run 66 stress 0.1016164 
    ## ... Procrustes: rmse 5.750627e-06  max resid 1.586028e-05 
    ## ... Similar to previous best
    ## Run 67 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006441183  max resid 0.00179963 
    ## ... Similar to previous best
    ## Run 68 stress 0.3411858 
    ## Run 69 stress 0.1016164 
    ## ... Procrustes: rmse 2.721803e-05  max resid 7.522325e-05 
    ## ... Similar to previous best
    ## Run 70 stress 0.1341868 
    ## Run 71 stress 0.1016164 
    ## ... Procrustes: rmse 3.101296e-06  max resid 7.894011e-06 
    ## ... Similar to previous best
    ## Run 72 stress 0.3495557 
    ## Run 73 stress 0.1530424 
    ## Run 74 stress 0.1356147 
    ## Run 75 stress 0.1356145 
    ## Run 76 stress 0.1016164 
    ## ... Procrustes: rmse 3.191389e-05  max resid 8.888348e-05 
    ## ... Similar to previous best
    ## Run 77 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001502315  max resid 0.0004199322 
    ## ... Similar to previous best
    ## Run 78 stress 0.1016164 
    ## ... Procrustes: rmse 2.420812e-05  max resid 6.72076e-05 
    ## ... Similar to previous best
    ## Run 79 stress 0.1016164 
    ## ... Procrustes: rmse 8.783179e-05  max resid 0.0002466826 
    ## ... Similar to previous best
    ## Run 80 stress 0.1341868 
    ## Run 81 stress 0.1016164 
    ## ... Procrustes: rmse 7.591302e-06  max resid 1.643437e-05 
    ## ... Similar to previous best
    ## Run 82 stress 0.1016164 
    ## ... Procrustes: rmse 0.000130796  max resid 0.0003689801 
    ## ... Similar to previous best
    ## Run 83 stress 0.1016164 
    ## ... Procrustes: rmse 2.39945e-05  max resid 6.729461e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.1016164 
    ## ... Procrustes: rmse 4.810024e-06  max resid 1.185605e-05 
    ## ... Similar to previous best
    ## Run 85 stress 0.1016164 
    ## ... Procrustes: rmse 2.224164e-05  max resid 6.249143e-05 
    ## ... Similar to previous best
    ## Run 86 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003811877  max resid 0.001067251 
    ## ... Similar to previous best
    ## Run 87 stress 0.1356148 
    ## Run 88 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004612111  max resid 0.00129105 
    ## ... Similar to previous best
    ## Run 89 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003495073  max resid 0.0009788534 
    ## ... Similar to previous best
    ## Run 90 stress 0.1016164 
    ## ... Procrustes: rmse 8.065023e-06  max resid 2.271427e-05 
    ## ... Similar to previous best
    ## Run 91 stress 0.1016164 
    ## ... Procrustes: rmse 0.000119301  max resid 0.0003370961 
    ## ... Similar to previous best
    ## Run 92 stress 0.1341868 
    ## Run 93 stress 0.1016168 
    ## ... Procrustes: rmse 0.000511925  max resid 0.001429043 
    ## ... Similar to previous best
    ## Run 94 stress 0.1016164 
    ## ... Procrustes: rmse 1.590128e-05  max resid 4.36358e-05 
    ## ... Similar to previous best
    ## Run 95 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005953354  max resid 0.001663724 
    ## ... Similar to previous best
    ## Run 96 stress 0.1016164 
    ## ... Procrustes: rmse 5.81705e-05  max resid 0.0001639639 
    ## ... Similar to previous best
    ## Run 97 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003709043  max resid 0.00103845 
    ## ... Similar to previous best
    ## Run 98 stress 0.1016164 
    ## ... Procrustes: rmse 2.541945e-05  max resid 6.012955e-05 
    ## ... Similar to previous best
    ## Run 99 stress 0.1341868 
    ## Run 100 stress 0.2074952 
    ## Run 101 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004056031  max resid 0.001134607 
    ## ... Similar to previous best
    ## Run 102 stress 0.101617 
    ## ... Procrustes: rmse 0.0006121101  max resid 0.001704513 
    ## ... Similar to previous best
    ## Run 103 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003260663  max resid 0.0009129781 
    ## ... Similar to previous best
    ## Run 104 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002768329  max resid 0.0007782238 
    ## ... Similar to previous best
    ## Run 105 stress 0.1016164 
    ## ... Procrustes: rmse 2.231795e-05  max resid 6.297142e-05 
    ## ... Similar to previous best
    ## Run 106 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004717708  max resid 0.00132 
    ## ... Similar to previous best
    ## Run 107 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 2.422468e-06  max resid 6.576871e-06 
    ## ... Similar to previous best
    ## Run 108 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 7.952828e-06  max resid 2.240458e-05 
    ## ... Similar to previous best
    ## Run 109 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005802469  max resid 0.001621805 
    ## ... Similar to previous best
    ## Run 110 stress 0.1016164 
    ## ... Procrustes: rmse 3.373427e-06  max resid 9.36514e-06 
    ## ... Similar to previous best
    ## Run 111 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003313798  max resid 0.0009268273 
    ## ... Similar to previous best
    ## Run 112 stress 0.1016164 
    ## ... Procrustes: rmse 3.038407e-05  max resid 8.567329e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003190413  max resid 0.00089445 
    ## ... Similar to previous best
    ## Run 114 stress 0.1341868 
    ## Run 115 stress 0.1016164 
    ## ... Procrustes: rmse 2.179877e-05  max resid 6.116632e-05 
    ## ... Similar to previous best
    ## Run 116 stress 0.1016164 
    ## ... Procrustes: rmse 8.198035e-05  max resid 0.0002309811 
    ## ... Similar to previous best
    ## Run 117 stress 0.1016164 
    ## ... Procrustes: rmse 9.508728e-06  max resid 2.679045e-05 
    ## ... Similar to previous best
    ## Run 118 stress 0.1530423 
    ## Run 119 stress 0.1356148 
    ## Run 120 stress 0.1341868 
    ## Run 121 stress 0.1016164 
    ## ... Procrustes: rmse 2.801147e-05  max resid 7.852259e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003263474  max resid 0.0009127601 
    ## ... Similar to previous best
    ## Run 123 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004991294  max resid 0.001395092 
    ## ... Similar to previous best
    ## Run 124 stress 0.1016164 
    ## ... Procrustes: rmse 4.353947e-06  max resid 1.030507e-05 
    ## ... Similar to previous best
    ## Run 125 stress 0.1016164 
    ## ... Procrustes: rmse 3.125615e-05  max resid 8.653676e-05 
    ## ... Similar to previous best
    ## Run 126 stress 0.1016171 
    ## ... Procrustes: rmse 0.0005842087  max resid 0.001623383 
    ## ... Similar to previous best
    ## Run 127 stress 0.1341868 
    ## Run 128 stress 0.1016164 
    ## ... Procrustes: rmse 4.842074e-05  max resid 0.0001362235 
    ## ... Similar to previous best
    ## Run 129 stress 0.1341868 
    ## Run 130 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001361875  max resid 0.0003838397 
    ## ... Similar to previous best
    ## Run 131 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005491312  max resid 0.001534822 
    ## ... Similar to previous best
    ## Run 132 stress 0.1356151 
    ## Run 133 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003885001  max resid 0.00108756 
    ## ... Similar to previous best
    ## Run 134 stress 0.1519033 
    ## Run 135 stress 0.1341868 
    ## Run 136 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005846569  max resid 0.001632628 
    ## ... Similar to previous best
    ## Run 137 stress 0.1341868 
    ## Run 138 stress 0.1016171 
    ## ... Procrustes: rmse 0.0006617719  max resid 0.001844537 
    ## ... Similar to previous best
    ## Run 139 stress 0.1356147 
    ## Run 140 stress 0.1341868 
    ## Run 141 stress 0.1341868 
    ## Run 142 stress 0.1016164 
    ## ... Procrustes: rmse 2.156623e-05  max resid 5.993183e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002085116  max resid 0.0005809417 
    ## ... Similar to previous best
    ## Run 144 stress 0.1341868 
    ## Run 145 stress 0.1341868 
    ## Run 146 stress 0.1016169 
    ## ... Procrustes: rmse 0.0004655774  max resid 0.001301923 
    ## ... Similar to previous best
    ## Run 147 stress 0.1016164 
    ## ... Procrustes: rmse 3.143169e-05  max resid 8.878629e-05 
    ## ... Similar to previous best
    ## Run 148 stress 0.1341868 
    ## Run 149 stress 0.1016164 
    ## ... Procrustes: rmse 5.61264e-05  max resid 0.0001581828 
    ## ... Similar to previous best
    ## Run 150 stress 0.1016164 
    ## ... Procrustes: rmse 3.160652e-05  max resid 8.856426e-05 
    ## ... Similar to previous best
    ## Run 151 stress 0.1016164 
    ## ... Procrustes: rmse 2.117913e-06  max resid 5.469493e-06 
    ## ... Similar to previous best
    ## Run 152 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005218472  max resid 0.001459148 
    ## ... Similar to previous best
    ## Run 153 stress 0.1016164 
    ## ... Procrustes: rmse 1.964705e-05  max resid 5.092478e-05 
    ## ... Similar to previous best
    ## Run 154 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005311216  max resid 0.001485403 
    ## ... Similar to previous best
    ## Run 155 stress 0.1356148 
    ## Run 156 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001833401  max resid 0.0005129921 
    ## ... Similar to previous best
    ## Run 157 stress 0.1341868 
    ## Run 158 stress 0.1016164 
    ## ... Procrustes: rmse 1.045128e-05  max resid 2.935561e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005339563  max resid 0.001492448 
    ## ... Similar to previous best
    ## Run 160 stress 0.1519033 
    ## Run 161 stress 0.1016164 
    ## ... Procrustes: rmse 9.636776e-05  max resid 0.0002716733 
    ## ... Similar to previous best
    ## Run 162 stress 0.1341868 
    ## Run 163 stress 0.1016164 
    ## ... Procrustes: rmse 8.971925e-05  max resid 0.0002535025 
    ## ... Similar to previous best
    ## Run 164 stress 0.1016164 
    ## ... Procrustes: rmse 4.6361e-05  max resid 0.0001308391 
    ## ... Similar to previous best
    ## Run 165 stress 0.1016164 
    ## ... Procrustes: rmse 4.100583e-05  max resid 0.0001155988 
    ## ... Similar to previous best
    ## Run 166 stress 0.1016164 
    ## ... Procrustes: rmse 9.207847e-05  max resid 0.0002591579 
    ## ... Similar to previous best
    ## Run 167 stress 0.1016164 
    ## ... Procrustes: rmse 6.380892e-05  max resid 0.0001801218 
    ## ... Similar to previous best
    ## Run 168 stress 0.1341868 
    ## Run 169 stress 0.1341868 
    ## Run 170 stress 0.1016164 
    ## ... Procrustes: rmse 7.820511e-05  max resid 0.000220442 
    ## ... Similar to previous best
    ## Run 171 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002317033  max resid 0.0006479343 
    ## ... Similar to previous best
    ## Run 172 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002042318  max resid 0.000573319 
    ## ... Similar to previous best
    ## Run 173 stress 0.1016164 
    ## ... Procrustes: rmse 1.918237e-05  max resid 3.754131e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004842293  max resid 0.001355193 
    ## ... Similar to previous best
    ## Run 175 stress 0.1016164 
    ## ... Procrustes: rmse 3.513841e-05  max resid 7.974134e-05 
    ## ... Similar to previous best
    ## Run 176 stress 0.1356148 
    ## Run 177 stress 0.1341868 
    ## Run 178 stress 0.1016164 
    ## ... Procrustes: rmse 8.911848e-05  max resid 0.0002510352 
    ## ... Similar to previous best
    ## Run 179 stress 0.135615 
    ## Run 180 stress 0.1537555 
    ## Run 181 stress 0.3495737 
    ## Run 182 stress 0.1341868 
    ## Run 183 stress 0.1016164 
    ## ... Procrustes: rmse 2.934291e-05  max resid 8.289316e-05 
    ## ... Similar to previous best
    ## Run 184 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004435258  max resid 0.001240978 
    ## ... Similar to previous best
    ## Run 185 stress 0.1356143 
    ## Run 186 stress 0.1016164 
    ## ... Procrustes: rmse 1.191368e-05  max resid 3.018883e-05 
    ## ... Similar to previous best
    ## Run 187 stress 0.1016164 
    ## ... Procrustes: rmse 0.000123759  max resid 0.0003466064 
    ## ... Similar to previous best
    ## Run 188 stress 0.1016164 
    ## ... Procrustes: rmse 2.65414e-05  max resid 4.109201e-05 
    ## ... Similar to previous best
    ## Run 189 stress 0.1341868 
    ## Run 190 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005716869  max resid 0.001592498 
    ## ... Similar to previous best
    ## Run 191 stress 0.1016164 
    ## ... Procrustes: rmse 5.493878e-05  max resid 0.0001549419 
    ## ... Similar to previous best
    ## Run 192 stress 0.1530423 
    ## Run 193 stress 0.1356146 
    ## Run 194 stress 0.1016164 
    ## ... Procrustes: rmse 1.231322e-05  max resid 3.452405e-05 
    ## ... Similar to previous best
    ## Run 195 stress 0.1356145 
    ## Run 196 stress 0.1016164 
    ## ... Procrustes: rmse 4.162542e-05  max resid 0.0001172166 
    ## ... Similar to previous best
    ## Run 197 stress 0.1016164 
    ## ... Procrustes: rmse 8.501048e-05  max resid 0.0002387429 
    ## ... Similar to previous best
    ## Run 198 stress 0.1016164 
    ## ... Procrustes: rmse 1.124803e-05  max resid 3.173405e-05 
    ## ... Similar to previous best
    ## Run 199 stress 0.1341868 
    ## Run 200 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003334578  max resid 0.0009327233 
    ## ... Similar to previous best
    ## Run 201 stress 0.2074952 
    ## Run 202 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005065605  max resid 0.001415925 
    ## ... Similar to previous best
    ## Run 203 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004810471  max resid 0.001345071 
    ## ... Similar to previous best
    ## Run 204 stress 0.1356144 
    ## Run 205 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005364538  max resid 0.001499221 
    ## ... Similar to previous best
    ## Run 206 stress 0.1341868 
    ## Run 207 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003920548  max resid 0.001096857 
    ## ... Similar to previous best
    ## Run 208 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002238706  max resid 0.0006278519 
    ## ... Similar to previous best
    ## Run 209 stress 0.1356151 
    ## Run 210 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 7.782609e-07  max resid 1.844671e-06 
    ## ... Similar to previous best
    ## Run 211 stress 0.1016164 
    ## ... Procrustes: rmse 2.231449e-05  max resid 6.189603e-05 
    ## ... Similar to previous best
    ## Run 212 stress 0.1016164 
    ## ... Procrustes: rmse 4.69599e-05  max resid 0.0001319955 
    ## ... Similar to previous best
    ## Run 213 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005472542  max resid 0.001530406 
    ## ... Similar to previous best
    ## Run 214 stress 0.1016164 
    ## ... Procrustes: rmse 4.360099e-05  max resid 0.0001228771 
    ## ... Similar to previous best
    ## Run 215 stress 0.3021671 
    ## Run 216 stress 0.1016164 
    ## ... Procrustes: rmse 6.980495e-06  max resid 1.059493e-05 
    ## ... Similar to previous best
    ## Run 217 stress 0.1341868 
    ## Run 218 stress 0.1341868 
    ## Run 219 stress 0.3056547 
    ## Run 220 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004951408  max resid 0.001378234 
    ## ... Similar to previous best
    ## Run 221 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003269026  max resid 0.0009131694 
    ## ... Similar to previous best
    ## Run 222 stress 0.1016164 
    ## ... Procrustes: rmse 7.134466e-05  max resid 0.0002014609 
    ## ... Similar to previous best
    ## Run 223 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004110401  max resid 0.001149892 
    ## ... Similar to previous best
    ## Run 224 stress 0.1341868 
    ## Run 225 stress 0.1016164 
    ## ... Procrustes: rmse 6.319644e-05  max resid 0.0001779333 
    ## ... Similar to previous best
    ## Run 226 stress 0.1016166 
    ## ... Procrustes: rmse 0.000376988  max resid 0.001053358 
    ## ... Similar to previous best
    ## Run 227 stress 0.1356148 
    ## Run 228 stress 0.1016164 
    ## ... Procrustes: rmse 6.659678e-06  max resid 1.884496e-05 
    ## ... Similar to previous best
    ## Run 229 stress 0.1341868 
    ## Run 230 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002162942  max resid 0.0006041472 
    ## ... Similar to previous best
    ## Run 231 stress 0.1016164 
    ## ... Procrustes: rmse 5.782962e-06  max resid 1.628299e-05 
    ## ... Similar to previous best
    ## Run 232 stress 0.1530422 
    ## Run 233 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005682816  max resid 0.001588354 
    ## ... Similar to previous best
    ## Run 234 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005610482  max resid 0.001567153 
    ## ... Similar to previous best
    ## Run 235 stress 0.1537555 
    ## Run 236 stress 0.1016164 
    ## ... Procrustes: rmse 8.37469e-05  max resid 0.0002352142 
    ## ... Similar to previous best
    ## Run 237 stress 0.1016164 
    ## ... Procrustes: rmse 9.26501e-06  max resid 2.503238e-05 
    ## ... Similar to previous best
    ## Run 238 stress 0.1016164 
    ## ... Procrustes: rmse 2.512658e-05  max resid 7.061768e-05 
    ## ... Similar to previous best
    ## Run 239 stress 0.3073147 
    ## Run 240 stress 0.1356145 
    ## Run 241 stress 0.1016164 
    ## ... Procrustes: rmse 1.954584e-05  max resid 5.215957e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.1016164 
    ## ... Procrustes: rmse 8.37196e-06  max resid 1.672749e-05 
    ## ... Similar to previous best
    ## Run 243 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002499505  max resid 0.0007004539 
    ## ... Similar to previous best
    ## Run 244 stress 0.1016164 
    ## ... Procrustes: rmse 2.898542e-05  max resid 8.182377e-05 
    ## ... Similar to previous best
    ## Run 245 stress 0.1530424 
    ## Run 246 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005308943  max resid 0.001483604 
    ## ... Similar to previous best
    ## Run 247 stress 0.1016167 
    ## ... Procrustes: rmse 0.000442874  max resid 0.001238194 
    ## ... Similar to previous best
    ## Run 248 stress 0.1356152 
    ## Run 249 stress 0.1356149 
    ## Run 250 stress 0.1341868 
    ## Run 251 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002421213  max resid 0.0006780193 
    ## ... Similar to previous best
    ## Run 252 stress 0.1016164 
    ## ... Procrustes: rmse 3.612924e-05  max resid 0.0001017422 
    ## ... Similar to previous best
    ## Run 253 stress 0.1016164 
    ## ... Procrustes: rmse 4.975913e-05  max resid 0.000139716 
    ## ... Similar to previous best
    ## Run 254 stress 0.1016164 
    ## ... Procrustes: rmse 1.704774e-05  max resid 4.77429e-05 
    ## ... Similar to previous best
    ## Run 255 stress 0.1530425 
    ## Run 256 stress 0.1016164 
    ## ... Procrustes: rmse 5.640317e-05  max resid 0.000148238 
    ## ... Similar to previous best
    ## Run 257 stress 0.1530428 
    ## Run 258 stress 0.1016164 
    ## ... Procrustes: rmse 1.979471e-05  max resid 5.561694e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.1519033 
    ## Run 260 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002182371  max resid 0.0006116182 
    ## ... Similar to previous best
    ## Run 261 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006086222  max resid 0.001700667 
    ## ... Similar to previous best
    ## Run 262 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001227083  max resid 0.0003450037 
    ## ... Similar to previous best
    ## Run 263 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003232781  max resid 0.0009031151 
    ## ... Similar to previous best
    ## Run 264 stress 0.1016167 
    ## ... Procrustes: rmse 0.000458891  max resid 0.001279393 
    ## ... Similar to previous best
    ## Run 265 stress 0.1341868 
    ## Run 266 stress 0.1016164 
    ## ... Procrustes: rmse 2.276468e-05  max resid 6.398932e-05 
    ## ... Similar to previous best
    ## Run 267 stress 0.1530422 
    ## Run 268 stress 0.1016164 
    ## ... Procrustes: rmse 1.387098e-05  max resid 3.731948e-05 
    ## ... Similar to previous best
    ## Run 269 stress 0.1530424 
    ## Run 270 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005561771  max resid 0.001554371 
    ## ... Similar to previous best
    ## Run 271 stress 0.3495732 
    ## Run 272 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003350848  max resid 0.0009372913 
    ## ... Similar to previous best
    ## Run 273 stress 0.1341868 
    ## Run 274 stress 0.1016164 
    ## ... Procrustes: rmse 6.388512e-05  max resid 0.0001801062 
    ## ... Similar to previous best
    ## Run 275 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003432875  max resid 0.0009602433 
    ## ... Similar to previous best
    ## Run 276 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003789099  max resid 0.001059142 
    ## ... Similar to previous best
    ## Run 277 stress 0.1016164 
    ## ... Procrustes: rmse 2.585519e-05  max resid 6.606572e-05 
    ## ... Similar to previous best
    ## Run 278 stress 0.1341868 
    ## Run 279 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005138996  max resid 0.001438001 
    ## ... Similar to previous best
    ## Run 280 stress 0.1356146 
    ## Run 281 stress 0.1016164 
    ## ... Procrustes: rmse 3.923126e-05  max resid 0.0001099545 
    ## ... Similar to previous best
    ## Run 282 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004816673  max resid 0.00134608 
    ## ... Similar to previous best
    ## Run 283 stress 0.1016164 
    ## ... Procrustes: rmse 6.151731e-05  max resid 0.0001727369 
    ## ... Similar to previous best
    ## Run 284 stress 0.1016164 
    ## ... Procrustes: rmse 9.708512e-06  max resid 2.739827e-05 
    ## ... Similar to previous best
    ## Run 285 stress 0.1356149 
    ## Run 286 stress 0.1016164 
    ## ... Procrustes: rmse 4.451808e-05  max resid 0.0001254183 
    ## ... Similar to previous best
    ## Run 287 stress 0.1341868 
    ## Run 288 stress 0.1016164 
    ## ... Procrustes: rmse 3.909613e-05  max resid 9.789005e-05 
    ## ... Similar to previous best
    ## Run 289 stress 0.1016164 
    ## ... Procrustes: rmse 4.269107e-05  max resid 0.000118709 
    ## ... Similar to previous best
    ## Run 290 stress 0.1341868 
    ## Run 291 stress 0.1016164 
    ## ... Procrustes: rmse 7.632639e-05  max resid 0.0002130986 
    ## ... Similar to previous best
    ## Run 292 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004145219  max resid 0.001158855 
    ## ... Similar to previous best
    ## Run 293 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003848329  max resid 0.001078907 
    ## ... Similar to previous best
    ## Run 294 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001497006  max resid 0.0004200933 
    ## ... Similar to previous best
    ## Run 295 stress 0.1341868 
    ## Run 296 stress 0.1016164 
    ## ... Procrustes: rmse 4.486535e-05  max resid 0.0001251811 
    ## ... Similar to previous best
    ## Run 297 stress 0.1356147 
    ## Run 298 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003302765  max resid 0.0009226917 
    ## ... Similar to previous best
    ## Run 299 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004802002  max resid 0.001342593 
    ## ... Similar to previous best
    ## Run 300 stress 0.1016164 
    ## ... Procrustes: rmse 2.497847e-05  max resid 6.983499e-05 
    ## ... Similar to previous best
    ## Run 301 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006004926  max resid 0.001676279 
    ## ... Similar to previous best
    ## Run 302 stress 0.1016164 
    ## ... Procrustes: rmse 1.368378e-05  max resid 3.81348e-05 
    ## ... Similar to previous best
    ## Run 303 stress 0.1016164 
    ## ... Procrustes: rmse 1.185105e-05  max resid 3.269146e-05 
    ## ... Similar to previous best
    ## Run 304 stress 0.1016164 
    ## ... Procrustes: rmse 1.74195e-05  max resid 4.84122e-05 
    ## ... Similar to previous best
    ## Run 305 stress 0.1016164 
    ## ... Procrustes: rmse 2.620372e-06  max resid 7.322833e-06 
    ## ... Similar to previous best
    ## Run 306 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003501599  max resid 0.0009718523 
    ## ... Similar to previous best
    ## Run 307 stress 0.1341868 
    ## Run 308 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005254125  max resid 0.001469122 
    ## ... Similar to previous best
    ## Run 309 stress 0.1341868 
    ## Run 310 stress 0.1016164 
    ## ... Procrustes: rmse 4.418044e-05  max resid 0.0001222945 
    ## ... Similar to previous best
    ## Run 311 stress 0.1016164 
    ## ... Procrustes: rmse 4.534951e-05  max resid 0.000127195 
    ## ... Similar to previous best
    ## Run 312 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 9.657549e-07  max resid 2.693804e-06 
    ## ... Similar to previous best
    ## Run 313 stress 0.1016164 
    ## ... Procrustes: rmse 9.381054e-05  max resid 0.0002642854 
    ## ... Similar to previous best
    ## Run 314 stress 0.101617 
    ## ... Procrustes: rmse 0.0005378244  max resid 0.001491964 
    ## ... Similar to previous best
    ## Run 315 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005689727  max resid 0.001591013 
    ## ... Similar to previous best
    ## Run 316 stress 0.1530424 
    ## Run 317 stress 0.1016164 
    ## ... Procrustes: rmse 2.422375e-05  max resid 6.806478e-05 
    ## ... Similar to previous best
    ## Run 318 stress 0.1016164 
    ## ... Procrustes: rmse 2.985705e-05  max resid 8.328194e-05 
    ## ... Similar to previous best
    ## Run 319 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004977783  max resid 0.001389579 
    ## ... Similar to previous best
    ## Run 320 stress 0.1016164 
    ## ... Procrustes: rmse 2.828928e-05  max resid 7.946481e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.1016164 
    ## ... Procrustes: rmse 1.158351e-05  max resid 3.290092e-05 
    ## ... Similar to previous best
    ## Run 322 stress 0.1016164 
    ## ... Procrustes: rmse 3.381288e-06  max resid 7.627949e-06 
    ## ... Similar to previous best
    ## Run 323 stress 0.1016164 
    ## ... Procrustes: rmse 2.109449e-05  max resid 5.316159e-05 
    ## ... Similar to previous best
    ## Run 324 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005005739  max resid 0.001393599 
    ## ... Similar to previous best
    ## Run 325 stress 0.1016164 
    ## ... Procrustes: rmse 1.109704e-06  max resid 2.567086e-06 
    ## ... Similar to previous best
    ## Run 326 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003783585  max resid 0.001053779 
    ## ... Similar to previous best
    ## Run 327 stress 0.1356145 
    ## Run 328 stress 0.1016166 
    ## ... Procrustes: rmse 0.000304269  max resid 0.0008521015 
    ## ... Similar to previous best
    ## Run 329 stress 0.1016164 
    ## ... Procrustes: rmse 7.316921e-05  max resid 0.0002066061 
    ## ... Similar to previous best
    ## Run 330 stress 0.1016164 
    ## ... Procrustes: rmse 1.523107e-05  max resid 4.256243e-05 
    ## ... Similar to previous best
    ## Run 331 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005302741  max resid 0.00147925 
    ## ... Similar to previous best
    ## Run 332 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006081059  max resid 0.001698955 
    ## ... Similar to previous best
    ## Run 333 stress 0.1016164 
    ## ... Procrustes: rmse 2.221821e-05  max resid 6.164259e-05 
    ## ... Similar to previous best
    ## Run 334 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003524433  max resid 0.0009868147 
    ## ... Similar to previous best
    ## Run 335 stress 0.1016164 
    ## ... Procrustes: rmse 4.136266e-06  max resid 8.465063e-06 
    ## ... Similar to previous best
    ## Run 336 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003598617  max resid 0.001008018 
    ## ... Similar to previous best
    ## Run 337 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001366839  max resid 0.0003819364 
    ## ... Similar to previous best
    ## Run 338 stress 0.1341868 
    ## Run 339 stress 0.1016164 
    ## ... Procrustes: rmse 3.216509e-05  max resid 7.435809e-05 
    ## ... Similar to previous best
    ## Run 340 stress 0.1016164 
    ## ... Procrustes: rmse 3.701528e-05  max resid 0.0001014349 
    ## ... Similar to previous best
    ## Run 341 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001019199  max resid 0.0002875112 
    ## ... Similar to previous best
    ## Run 342 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004082616  max resid 0.001143203 
    ## ... Similar to previous best
    ## Run 343 stress 0.1356151 
    ## Run 344 stress 0.1530425 
    ## Run 345 stress 0.1341868 
    ## Run 346 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004211498  max resid 0.001178683 
    ## ... Similar to previous best
    ## Run 347 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004312069  max resid 0.001205588 
    ## ... Similar to previous best
    ## Run 348 stress 0.1356147 
    ## Run 349 stress 0.1016164 
    ## ... Procrustes: rmse 3.88066e-05  max resid 0.0001079929 
    ## ... Similar to previous best
    ## Run 350 stress 0.3578487 
    ## Run 351 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005211566  max resid 0.001456314 
    ## ... Similar to previous best
    ## Run 352 stress 0.1016164 
    ## ... Procrustes: rmse 2.006261e-06  max resid 5.217642e-06 
    ## ... Similar to previous best
    ## Run 353 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005030065  max resid 0.001406693 
    ## ... Similar to previous best
    ## Run 354 stress 0.1530424 
    ## Run 355 stress 0.1356153 
    ## Run 356 stress 0.1356148 
    ## Run 357 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002017074  max resid 0.0005643663 
    ## ... Similar to previous best
    ## Run 358 stress 0.1016164 
    ## ... Procrustes: rmse 7.028932e-06  max resid 1.953499e-05 
    ## ... Similar to previous best
    ## Run 359 stress 0.1016164 
    ## ... Procrustes: rmse 2.813765e-05  max resid 7.868362e-05 
    ## ... Similar to previous best
    ## Run 360 stress 0.1341868 
    ## Run 361 stress 0.1016164 
    ## ... Procrustes: rmse 7.986706e-05  max resid 0.0002248577 
    ## ... Similar to previous best
    ## Run 362 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005794318  max resid 0.00161766 
    ## ... Similar to previous best
    ## Run 363 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003497398  max resid 0.0009781037 
    ## ... Similar to previous best
    ## Run 364 stress 0.1016164 
    ## ... Procrustes: rmse 8.866862e-05  max resid 0.0002496965 
    ## ... Similar to previous best
    ## Run 365 stress 0.1016164 
    ## ... Procrustes: rmse 1.19263e-05  max resid 3.361785e-05 
    ## ... Similar to previous best
    ## Run 366 stress 0.1341868 
    ## Run 367 stress 0.1016164 
    ## ... Procrustes: rmse 7.88194e-06  max resid 1.80052e-05 
    ## ... Similar to previous best
    ## Run 368 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004955092  max resid 0.001384569 
    ## ... Similar to previous best
    ## Run 369 stress 0.1016164 
    ## ... Procrustes: rmse 1.100231e-05  max resid 2.986461e-05 
    ## ... Similar to previous best
    ## Run 370 stress 0.1341868 
    ## Run 371 stress 0.1016167 
    ## ... Procrustes: rmse 0.000407004  max resid 0.00113644 
    ## ... Similar to previous best
    ## Run 372 stress 0.1341868 
    ## Run 373 stress 0.1016164 
    ## ... Procrustes: rmse 8.020594e-06  max resid 2.269939e-05 
    ## ... Similar to previous best
    ## Run 374 stress 0.1016164 
    ## ... Procrustes: rmse 4.919301e-05  max resid 0.0001385239 
    ## ... Similar to previous best
    ## Run 375 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001029923  max resid 0.000288617 
    ## ... Similar to previous best
    ## Run 376 stress 0.1016164 
    ## ... Procrustes: rmse 6.892489e-05  max resid 0.0001937573 
    ## ... Similar to previous best
    ## Run 377 stress 0.1016164 
    ## ... Procrustes: rmse 2.50358e-05  max resid 7.007016e-05 
    ## ... Similar to previous best
    ## Run 378 stress 0.1016164 
    ## ... Procrustes: rmse 1.55152e-05  max resid 4.371542e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004349791  max resid 0.00121537 
    ## ... Similar to previous best
    ## Run 380 stress 0.1016164 
    ## ... Procrustes: rmse 9.243988e-05  max resid 0.0002591642 
    ## ... Similar to previous best
    ## Run 381 stress 0.1341868 
    ## Run 382 stress 0.1016164 
    ## ... Procrustes: rmse 6.432372e-06  max resid 1.163132e-05 
    ## ... Similar to previous best
    ## Run 383 stress 0.1016164 
    ## ... Procrustes: rmse 4.285005e-05  max resid 0.0001207041 
    ## ... Similar to previous best
    ## Run 384 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004496369  max resid 0.001257426 
    ## ... Similar to previous best
    ## Run 385 stress 0.1341868 
    ## Run 386 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004576593  max resid 0.00127877 
    ## ... Similar to previous best
    ## Run 387 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005923471  max resid 0.001655757 
    ## ... Similar to previous best
    ## Run 388 stress 0.1016164 
    ## ... Procrustes: rmse 4.121555e-05  max resid 0.0001163747 
    ## ... Similar to previous best
    ## Run 389 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004401628  max resid 0.001231685 
    ## ... Similar to previous best
    ## Run 390 stress 0.1341868 
    ## Run 391 stress 0.1016164 
    ## ... Procrustes: rmse 8.080862e-05  max resid 0.0002278583 
    ## ... Similar to previous best
    ## Run 392 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002123106  max resid 0.0005900662 
    ## ... Similar to previous best
    ## Run 393 stress 0.1356144 
    ## Run 394 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004572698  max resid 0.001278161 
    ## ... Similar to previous best
    ## Run 395 stress 0.1016164 
    ## ... Procrustes: rmse 6.454883e-06  max resid 1.809815e-05 
    ## ... Similar to previous best
    ## Run 396 stress 0.1016164 
    ## ... Procrustes: rmse 2.43264e-05  max resid 6.865506e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005825728  max resid 0.001628705 
    ## ... Similar to previous best
    ## Run 398 stress 0.1016164 
    ## ... Procrustes: rmse 2.052393e-05  max resid 4.726276e-05 
    ## ... Similar to previous best
    ## Run 399 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006343954  max resid 0.001772671 
    ## ... Similar to previous best
    ## Run 400 stress 0.1016164 
    ## ... Procrustes: rmse 7.429757e-05  max resid 0.0002091672 
    ## ... Similar to previous best
    ## Run 401 stress 0.1530426 
    ## Run 402 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001442935  max resid 0.000402673 
    ## ... Similar to previous best
    ## Run 403 stress 0.1341868 
    ## Run 404 stress 0.1341868 
    ## Run 405 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004432131  max resid 0.001237376 
    ## ... Similar to previous best
    ## Run 406 stress 0.1016164 
    ## ... Procrustes: rmse 1.639245e-05  max resid 4.18897e-05 
    ## ... Similar to previous best
    ## Run 407 stress 0.1016168 
    ## ... Procrustes: rmse 0.000539451  max resid 0.001508721 
    ## ... Similar to previous best
    ## Run 408 stress 0.1356152 
    ## Run 409 stress 0.1016164 
    ## ... Procrustes: rmse 6.230023e-05  max resid 0.0001753052 
    ## ... Similar to previous best
    ## Run 410 stress 0.1016164 
    ## ... Procrustes: rmse 5.67416e-06  max resid 1.564791e-05 
    ## ... Similar to previous best
    ## Run 411 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004448582  max resid 0.001244017 
    ## ... Similar to previous best
    ## Run 412 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005031962  max resid 0.001405098 
    ## ... Similar to previous best
    ## Run 413 stress 0.1016164 
    ## ... Procrustes: rmse 2.578129e-05  max resid 7.221866e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.212628 
    ## Run 415 stress 0.1016164 
    ## ... Procrustes: rmse 5.925176e-06  max resid 1.65466e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.101617 
    ## ... Procrustes: rmse 0.0006013779  max resid 0.001675821 
    ## ... Similar to previous best
    ## Run 417 stress 0.1341868 
    ## Run 418 stress 0.135615 
    ## Run 419 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001487842  max resid 0.0004188 
    ## ... Similar to previous best
    ## Run 420 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004175971  max resid 0.001167896 
    ## ... Similar to previous best
    ## Run 421 stress 0.1016164 
    ## ... Procrustes: rmse 4.261468e-06  max resid 1.19725e-05 
    ## ... Similar to previous best
    ## Run 422 stress 0.1341868 
    ## Run 423 stress 0.1016164 
    ## ... Procrustes: rmse 7.242456e-05  max resid 0.0002042249 
    ## ... Similar to previous best
    ## Run 424 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005094994  max resid 0.001425053 
    ## ... Similar to previous best
    ## Run 425 stress 0.1016164 
    ## ... Procrustes: rmse 3.404349e-05  max resid 9.58305e-05 
    ## ... Similar to previous best
    ## Run 426 stress 0.1341868 
    ## Run 427 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001083472  max resid 0.0003049676 
    ## ... Similar to previous best
    ## Run 428 stress 0.1016164 
    ## ... Procrustes: rmse 5.141895e-05  max resid 0.0001388442 
    ## ... Similar to previous best
    ## Run 429 stress 0.1341868 
    ## Run 430 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005164551  max resid 0.00144401 
    ## ... Similar to previous best
    ## Run 431 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003318441  max resid 0.0009285174 
    ## ... Similar to previous best
    ## Run 432 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003398221  max resid 0.0009513821 
    ## ... Similar to previous best
    ## Run 433 stress 0.1016164 
    ## ... Procrustes: rmse 6.880169e-05  max resid 0.0001935029 
    ## ... Similar to previous best
    ## Run 434 stress 0.1016164 
    ## ... Procrustes: rmse 3.455796e-05  max resid 9.093476e-05 
    ## ... Similar to previous best
    ## Run 435 stress 0.1016164 
    ## ... Procrustes: rmse 9.044387e-06  max resid 2.516373e-05 
    ## ... Similar to previous best
    ## Run 436 stress 0.1341868 
    ## Run 437 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004872865  max resid 0.001362165 
    ## ... Similar to previous best
    ## Run 438 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004189107  max resid 0.001172398 
    ## ... Similar to previous best
    ## Run 439 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003285314  max resid 0.0009201297 
    ## ... Similar to previous best
    ## Run 440 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003332134  max resid 0.0009330865 
    ## ... Similar to previous best
    ## Run 441 stress 0.1016164 
    ## ... Procrustes: rmse 4.409749e-05  max resid 0.000117168 
    ## ... Similar to previous best
    ## Run 442 stress 0.1016164 
    ## ... Procrustes: rmse 1.729462e-05  max resid 4.862356e-05 
    ## ... Similar to previous best
    ## Run 443 stress 0.1016164 
    ## ... Procrustes: rmse 3.791843e-06  max resid 1.075908e-05 
    ## ... Similar to previous best
    ## Run 444 stress 0.1356147 
    ## Run 445 stress 0.1016164 
    ## ... Procrustes: rmse 1.874971e-05  max resid 5.288506e-05 
    ## ... Similar to previous best
    ## Run 446 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004956095  max resid 0.001382746 
    ## ... Similar to previous best
    ## Run 447 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002388137  max resid 0.0006702123 
    ## ... Similar to previous best
    ## Run 448 stress 0.1356147 
    ## Run 449 stress 0.1356145 
    ## Run 450 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002087405  max resid 0.000585818 
    ## ... Similar to previous best
    ## Run 451 stress 0.1356146 
    ## Run 452 stress 0.1016164 
    ## ... Procrustes: rmse 1.377305e-05  max resid 3.933126e-05 
    ## ... Similar to previous best
    ## Run 453 stress 0.1016164 
    ## ... Procrustes: rmse 4.482492e-05  max resid 0.0001258977 
    ## ... Similar to previous best
    ## Run 454 stress 0.1016164 
    ## ... Procrustes: rmse 5.043654e-05  max resid 0.0001417515 
    ## ... Similar to previous best
    ## Run 455 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004332398  max resid 0.001211605 
    ## ... Similar to previous best
    ## Run 456 stress 0.1016164 
    ## ... Procrustes: rmse 1.299296e-05  max resid 3.630887e-05 
    ## ... Similar to previous best
    ## Run 457 stress 0.1016164 
    ## ... Procrustes: rmse 1.877372e-05  max resid 5.295208e-05 
    ## ... Similar to previous best
    ## Run 458 stress 0.1016164 
    ## ... Procrustes: rmse 3.93082e-05  max resid 0.0001107561 
    ## ... Similar to previous best
    ## Run 459 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004520143  max resid 0.001264326 
    ## ... Similar to previous best
    ## Run 460 stress 0.1356146 
    ## Run 461 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004434911  max resid 0.001241811 
    ## ... Similar to previous best
    ## Run 462 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002725569  max resid 0.0007639679 
    ## ... Similar to previous best
    ## Run 463 stress 0.1341868 
    ## Run 464 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003217536  max resid 0.0009004632 
    ## ... Similar to previous best
    ## Run 465 stress 0.1356147 
    ## Run 466 stress 0.1341868 
    ## Run 467 stress 0.1016164 
    ## ... Procrustes: rmse 4.763365e-05  max resid 0.0001345201 
    ## ... Similar to previous best
    ## Run 468 stress 0.1016164 
    ## ... Procrustes: rmse 6.090845e-06  max resid 1.715169e-05 
    ## ... Similar to previous best
    ## Run 469 stress 0.1016164 
    ## ... Procrustes: rmse 3.147003e-05  max resid 8.733605e-05 
    ## ... Similar to previous best
    ## Run 470 stress 0.101617 
    ## ... Procrustes: rmse 0.0006716774  max resid 0.001876506 
    ## ... Similar to previous best
    ## Run 471 stress 0.1356144 
    ## Run 472 stress 0.1530424 
    ## Run 473 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004932172  max resid 0.001379406 
    ## ... Similar to previous best
    ## Run 474 stress 0.1016164 
    ## ... Procrustes: rmse 1.750414e-05  max resid 4.915623e-05 
    ## ... Similar to previous best
    ## Run 475 stress 0.1341868 
    ## Run 476 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005659254  max resid 0.001581728 
    ## ... Similar to previous best
    ## Run 477 stress 0.1341868 
    ## Run 478 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003342932  max resid 0.0009341875 
    ## ... Similar to previous best
    ## Run 479 stress 0.1341868 
    ## Run 480 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003875326  max resid 0.001084295 
    ## ... Similar to previous best
    ## Run 481 stress 0.1356149 
    ## Run 482 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001454675  max resid 0.0004072711 
    ## ... Similar to previous best
    ## Run 483 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005024854  max resid 0.001406597 
    ## ... Similar to previous best
    ## Run 484 stress 0.1016169 
    ## ... Procrustes: rmse 0.0004612486  max resid 0.001288084 
    ## ... Similar to previous best
    ## Run 485 stress 0.1356145 
    ## Run 486 stress 0.1016164 
    ## ... Procrustes: rmse 2.084855e-05  max resid 5.881817e-05 
    ## ... Similar to previous best
    ## Run 487 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006237335  max resid 0.001744153 
    ## ... Similar to previous best
    ## Run 488 stress 0.1341868 
    ## Run 489 stress 0.1016164 
    ## ... Procrustes: rmse 1.914405e-05  max resid 5.387751e-05 
    ## ... Similar to previous best
    ## Run 490 stress 0.1356144 
    ## Run 491 stress 0.1341868 
    ## Run 492 stress 0.1341868 
    ## Run 493 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004244926  max resid 0.001188261 
    ## ... Similar to previous best
    ## Run 494 stress 0.1356144 
    ## Run 495 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004598436  max resid 0.001286395 
    ## ... Similar to previous best
    ## Run 496 stress 0.1341868 
    ## Run 497 stress 0.1016164 
    ## ... Procrustes: rmse 0.000100305  max resid 0.0002816958 
    ## ... Similar to previous best
    ## Run 498 stress 0.1016167 
    ## ... Procrustes: rmse 0.000473559  max resid 0.001324646 
    ## ... Similar to previous best
    ## Run 499 stress 0.1341868 
    ## Run 500 stress 0.1016164 
    ## ... Procrustes: rmse 7.612745e-05  max resid 0.0002094876 
    ## ... Similar to previous best
    ## *** Best solution repeated 137 times

``` r
# Stratified lakes and ocean sites
PD_beta_SO_NMDS <- metaMDS(PD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.04851417 
    ## Run 1 stress 0.05607346 
    ## Run 2 stress 0.04938616 
    ## Run 3 stress 0.04860423 
    ## ... Procrustes: rmse 0.05619241  max resid 0.1295428 
    ## Run 4 stress 0.04612762 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04983454  max resid 0.1395702 
    ## Run 5 stress 0.04851417 
    ## Run 6 stress 0.05451494 
    ## Run 7 stress 0.04938616 
    ## Run 8 stress 0.0608897 
    ## Run 9 stress 0.04612759 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001855547  max resid 0.0003792552 
    ## ... Similar to previous best
    ## Run 10 stress 0.04935945 
    ## Run 11 stress 0.05605146 
    ## Run 12 stress 0.05593455 
    ## Run 13 stress 0.05260054 
    ## Run 14 stress 0.05607345 
    ## Run 15 stress 0.05695875 
    ## Run 16 stress 0.05643261 
    ## Run 17 stress 0.05728066 
    ## Run 18 stress 0.0521613 
    ## Run 19 stress 0.05593456 
    ## Run 20 stress 0.058428 
    ## Run 21 stress 0.05728082 
    ## Run 22 stress 0.0557269 
    ## Run 23 stress 0.05842783 
    ## Run 24 stress 0.05137965 
    ## Run 25 stress 0.05260053 
    ## Run 26 stress 0.04612762 
    ## ... Procrustes: rmse 0.0001800292  max resid 0.0003572982 
    ## ... Similar to previous best
    ## Run 27 stress 0.05447778 
    ## Run 28 stress 0.05643261 
    ## Run 29 stress 0.04851417 
    ## Run 30 stress 0.05695864 
    ## Run 31 stress 0.05695868 
    ## Run 32 stress 0.06933845 
    ## Run 33 stress 0.04612762 
    ## ... Procrustes: rmse 0.0001870188  max resid 0.0003798978 
    ## ... Similar to previous best
    ## Run 34 stress 0.05728055 
    ## Run 35 stress 0.07046342 
    ## Run 36 stress 0.05474821 
    ## Run 37 stress 0.04938615 
    ## Run 38 stress 0.05572688 
    ## Run 39 stress 0.05572688 
    ## Run 40 stress 0.05612254 
    ## Run 41 stress 0.0527641 
    ## Run 42 stress 0.05439281 
    ## Run 43 stress 0.0585503 
    ## Run 44 stress 0.04938618 
    ## Run 45 stress 0.0527641 
    ## Run 46 stress 0.05643261 
    ## Run 47 stress 0.04935946 
    ## Run 48 stress 0.04851417 
    ## Run 49 stress 0.05612252 
    ## Run 50 stress 0.05605153 
    ## Run 51 stress 0.04851417 
    ## Run 52 stress 0.05612239 
    ## Run 53 stress 0.05605152 
    ## Run 54 stress 0.04851417 
    ## Run 55 stress 0.0543894 
    ## Run 56 stress 0.04851419 
    ## Run 57 stress 0.0543894 
    ## Run 58 stress 0.0640074 
    ## Run 59 stress 0.06088975 
    ## Run 60 stress 0.05728059 
    ## Run 61 stress 0.05855038 
    ## Run 62 stress 0.05474819 
    ## Run 63 stress 0.05216126 
    ## Run 64 stress 0.05216156 
    ## Run 65 stress 0.05607344 
    ## Run 66 stress 0.05612243 
    ## Run 67 stress 0.04860422 
    ## Run 68 stress 0.05959896 
    ## Run 69 stress 0.04938616 
    ## Run 70 stress 0.357849 
    ## Run 71 stress 0.05607345 
    ## Run 72 stress 0.05607349 
    ## Run 73 stress 0.05612246 
    ## Run 74 stress 0.05474815 
    ## Run 75 stress 0.04938616 
    ## Run 76 stress 0.05451496 
    ## Run 77 stress 0.05959896 
    ## Run 78 stress 0.0547482 
    ## Run 79 stress 0.05593449 
    ## Run 80 stress 0.05605148 
    ## Run 81 stress 0.05643264 
    ## Run 82 stress 0.05643261 
    ## Run 83 stress 0.05643261 
    ## Run 84 stress 0.05260057 
    ## Run 85 stress 0.04935945 
    ## Run 86 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 8.621064e-05  max resid 0.000175482 
    ## ... Similar to previous best
    ## Run 87 stress 0.05607344 
    ## Run 88 stress 0.05643263 
    ## Run 89 stress 0.05643261 
    ## Run 90 stress 0.04851418 
    ## Run 91 stress 0.04860423 
    ## Run 92 stress 0.0572806 
    ## Run 93 stress 0.05593449 
    ## Run 94 stress 0.04860422 
    ## Run 95 stress 0.06933817 
    ## Run 96 stress 0.05605147 
    ## Run 97 stress 0.04851417 
    ## Run 98 stress 0.04612766 
    ## ... Procrustes: rmse 0.0001390434  max resid 0.0002803391 
    ## ... Similar to previous best
    ## Run 99 stress 0.04851417 
    ## Run 100 stress 0.05842796 
    ## Run 101 stress 0.05433173 
    ## Run 102 stress 0.05572688 
    ## Run 103 stress 0.05137974 
    ## Run 104 stress 0.05451493 
    ## Run 105 stress 0.0608897 
    ## Run 106 stress 0.05605151 
    ## Run 107 stress 0.04612759 
    ## ... Procrustes: rmse 8.235906e-05  max resid 0.0001662118 
    ## ... Similar to previous best
    ## Run 108 stress 0.04851417 
    ## Run 109 stress 0.05439281 
    ## Run 110 stress 0.06205755 
    ## Run 111 stress 0.05216171 
    ## Run 112 stress 0.05855031 
    ## Run 113 stress 0.0559346 
    ## Run 114 stress 0.05276411 
    ## Run 115 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001297581  max resid 0.0002610396 
    ## ... Similar to previous best
    ## Run 116 stress 0.04851417 
    ## Run 117 stress 0.04938616 
    ## Run 118 stress 0.04851418 
    ## Run 119 stress 0.0543928 
    ## Run 120 stress 0.04938615 
    ## Run 121 stress 0.04860423 
    ## Run 122 stress 0.05216128 
    ## Run 123 stress 0.05572688 
    ## Run 124 stress 0.04938615 
    ## Run 125 stress 0.04851417 
    ## Run 126 stress 0.06205753 
    ## Run 127 stress 0.04938616 
    ## Run 128 stress 0.04860422 
    ## Run 129 stress 0.05260053 
    ## Run 130 stress 0.06205752 
    ## Run 131 stress 0.0521614 
    ## Run 132 stress 0.05439287 
    ## Run 133 stress 0.05612242 
    ## Run 134 stress 0.04860423 
    ## Run 135 stress 0.05593456 
    ## Run 136 stress 0.04935947 
    ## Run 137 stress 0.05728053 
    ## Run 138 stress 0.05607344 
    ## Run 139 stress 0.05605145 
    ## Run 140 stress 0.0527641 
    ## Run 141 stress 0.0526005 
    ## Run 142 stress 0.05643263 
    ## Run 143 stress 0.05728055 
    ## Run 144 stress 0.05605152 
    ## Run 145 stress 0.05488105 
    ## Run 146 stress 0.05451493 
    ## Run 147 stress 0.0526005 
    ## Run 148 stress 0.05137966 
    ## Run 149 stress 0.04851419 
    ## Run 150 stress 0.05728053 
    ## Run 151 stress 0.05612242 
    ## Run 152 stress 0.05607344 
    ## Run 153 stress 0.04938616 
    ## Run 154 stress 0.05216138 
    ## Run 155 stress 0.06400732 
    ## Run 156 stress 0.04938616 
    ## Run 157 stress 0.3291585 
    ## Run 158 stress 0.04938617 
    ## Run 159 stress 0.04851417 
    ## Run 160 stress 0.04851418 
    ## Run 161 stress 0.04851418 
    ## Run 162 stress 0.05488105 
    ## Run 163 stress 0.06185054 
    ## Run 164 stress 0.05439283 
    ## Run 165 stress 0.05451498 
    ## Run 166 stress 0.05728052 
    ## Run 167 stress 0.04938617 
    ## Run 168 stress 0.05607349 
    ## Run 169 stress 0.05593454 
    ## Run 170 stress 0.05216134 
    ## Run 171 stress 0.05842786 
    ## Run 172 stress 0.3254264 
    ## Run 173 stress 0.04938615 
    ## Run 174 stress 0.05474814 
    ## Run 175 stress 0.05216161 
    ## Run 176 stress 0.05605157 
    ## Run 177 stress 0.05451493 
    ## Run 178 stress 0.04612766 
    ## ... Procrustes: rmse 0.0001358056  max resid 0.0002770948 
    ## ... Similar to previous best
    ## Run 179 stress 0.04938616 
    ## Run 180 stress 0.05728058 
    ## Run 181 stress 0.04851417 
    ## Run 182 stress 0.05572688 
    ## Run 183 stress 0.0640073 
    ## Run 184 stress 0.05488111 
    ## Run 185 stress 0.05488104 
    ## Run 186 stress 0.04612761 
    ## ... Procrustes: rmse 9.152544e-05  max resid 0.0001863133 
    ## ... Similar to previous best
    ## Run 187 stress 0.05605148 
    ## Run 188 stress 0.04851417 
    ## Run 189 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001185665  max resid 0.000238113 
    ## ... Similar to previous best
    ## Run 190 stress 0.05474817 
    ## Run 191 stress 0.3397485 
    ## Run 192 stress 0.0543928 
    ## Run 193 stress 0.04860422 
    ## Run 194 stress 0.04851417 
    ## Run 195 stress 0.04851417 
    ## Run 196 stress 0.05728053 
    ## Run 197 stress 0.05433173 
    ## Run 198 stress 0.04938617 
    ## Run 199 stress 0.06933831 
    ## Run 200 stress 0.04938616 
    ## Run 201 stress 0.05260055 
    ## Run 202 stress 0.05728056 
    ## Run 203 stress 0.05216135 
    ## Run 204 stress 0.05474816 
    ## Run 205 stress 0.04938616 
    ## Run 206 stress 0.06400748 
    ## Run 207 stress 0.05451493 
    ## Run 208 stress 0.04851417 
    ## Run 209 stress 0.05842805 
    ## Run 210 stress 0.05137964 
    ## Run 211 stress 0.04938616 
    ## Run 212 stress 0.04612759 
    ## ... Procrustes: rmse 8.195877e-05  max resid 0.000165196 
    ## ... Similar to previous best
    ## Run 213 stress 0.05728053 
    ## Run 214 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001242088  max resid 0.000250984 
    ## ... Similar to previous best
    ## Run 215 stress 0.05593466 
    ## Run 216 stress 0.05451496 
    ## Run 217 stress 0.05216143 
    ## Run 218 stress 0.05474814 
    ## Run 219 stress 0.04935953 
    ## Run 220 stress 0.04612757 
    ## ... Procrustes: rmse 7.749857e-06  max resid 1.46354e-05 
    ## ... Similar to previous best
    ## Run 221 stress 0.04938616 
    ## Run 222 stress 0.04860422 
    ## Run 223 stress 0.04938616 
    ## Run 224 stress 0.04851417 
    ## Run 225 stress 0.04612757 
    ## ... Procrustes: rmse 5.005341e-05  max resid 0.0001013373 
    ## ... Similar to previous best
    ## Run 226 stress 0.05216154 
    ## Run 227 stress 0.04938617 
    ## Run 228 stress 0.06291424 
    ## Run 229 stress 0.06205754 
    ## Run 230 stress 0.05605145 
    ## Run 231 stress 0.0527641 
    ## Run 232 stress 0.05439281 
    ## Run 233 stress 0.05695881 
    ## Run 234 stress 0.05612241 
    ## Run 235 stress 0.3305477 
    ## Run 236 stress 0.05607345 
    ## Run 237 stress 0.05643261 
    ## Run 238 stress 0.05607344 
    ## Run 239 stress 0.05474813 
    ## Run 240 stress 0.04938616 
    ## Run 241 stress 0.06088972 
    ## Run 242 stress 0.05607344 
    ## Run 243 stress 0.05572688 
    ## Run 244 stress 0.04935947 
    ## Run 245 stress 0.05572691 
    ## Run 246 stress 0.05593462 
    ## Run 247 stress 0.05607346 
    ## Run 248 stress 0.04860422 
    ## Run 249 stress 0.05451493 
    ## Run 250 stress 0.05643262 
    ## Run 251 stress 0.0585503 
    ## Run 252 stress 0.05447819 
    ## Run 253 stress 0.05474813 
    ## Run 254 stress 0.04612764 
    ## ... Procrustes: rmse 0.000128663  max resid 0.0002589984 
    ## ... Similar to previous best
    ## Run 255 stress 0.05643263 
    ## Run 256 stress 0.05607345 
    ## Run 257 stress 0.05607344 
    ## Run 258 stress 0.06291425 
    ## Run 259 stress 0.0585503 
    ## Run 260 stress 0.04612759 
    ## ... Procrustes: rmse 8.745867e-05  max resid 0.000176844 
    ## ... Similar to previous best
    ## Run 261 stress 0.0493862 
    ## Run 262 stress 0.05728057 
    ## Run 263 stress 0.0561224 
    ## Run 264 stress 0.0585503 
    ## Run 265 stress 0.04851417 
    ## Run 266 stress 0.05260049 
    ## Run 267 stress 0.05855032 
    ## Run 268 stress 0.05728068 
    ## Run 269 stress 0.04938616 
    ## Run 270 stress 0.05607345 
    ## Run 271 stress 0.0521615 
    ## Run 272 stress 0.05607345 
    ## Run 273 stress 0.04938616 
    ## Run 274 stress 0.04612767 
    ## ... Procrustes: rmse 0.0001611675  max resid 0.0003292013 
    ## ... Similar to previous best
    ## Run 275 stress 0.05451497 
    ## Run 276 stress 0.04938625 
    ## Run 277 stress 0.04938615 
    ## Run 278 stress 0.05438938 
    ## Run 279 stress 0.05260057 
    ## Run 280 stress 0.05276411 
    ## Run 281 stress 0.06185068 
    ## Run 282 stress 0.05474814 
    ## Run 283 stress 0.04935949 
    ## Run 284 stress 0.05451493 
    ## Run 285 stress 0.0543928 
    ## Run 286 stress 0.04938616 
    ## Run 287 stress 0.05605145 
    ## Run 288 stress 0.04860423 
    ## Run 289 stress 0.05593458 
    ## Run 290 stress 0.05695863 
    ## Run 291 stress 0.0527641 
    ## Run 292 stress 0.05474815 
    ## Run 293 stress 0.06291424 
    ## Run 294 stress 0.04935949 
    ## Run 295 stress 0.0493862 
    ## Run 296 stress 0.04938615 
    ## Run 297 stress 0.04935943 
    ## Run 298 stress 0.04612761 
    ## ... Procrustes: rmse 0.0001006727  max resid 0.0001982518 
    ## ... Similar to previous best
    ## Run 299 stress 0.05643261 
    ## Run 300 stress 0.05451493 
    ## Run 301 stress 0.04938615 
    ## Run 302 stress 0.05433173 
    ## Run 303 stress 0.04851417 
    ## Run 304 stress 0.0560515 
    ## Run 305 stress 0.05451493 
    ## Run 306 stress 0.04851417 
    ## Run 307 stress 0.05728083 
    ## Run 308 stress 0.05728052 
    ## Run 309 stress 0.04851417 
    ## Run 310 stress 0.04938616 
    ## Run 311 stress 0.05216142 
    ## Run 312 stress 0.05572688 
    ## Run 313 stress 0.05572688 
    ## Run 314 stress 0.05451493 
    ## Run 315 stress 0.05695859 
    ## Run 316 stress 0.05137973 
    ## Run 317 stress 0.05728064 
    ## Run 318 stress 0.05959896 
    ## Run 319 stress 0.05728052 
    ## Run 320 stress 0.05572688 
    ## Run 321 stress 0.05451493 
    ## Run 322 stress 0.04851418 
    ## Run 323 stress 0.04851417 
    ## Run 324 stress 0.05137968 
    ## Run 325 stress 0.05643261 
    ## Run 326 stress 0.0461276 
    ## ... Procrustes: rmse 6.91154e-05  max resid 0.0001414513 
    ## ... Similar to previous best
    ## Run 327 stress 0.05643267 
    ## Run 328 stress 0.04935949 
    ## Run 329 stress 0.05607345 
    ## Run 330 stress 0.04612759 
    ## ... Procrustes: rmse 6.708757e-05  max resid 0.0001350821 
    ## ... Similar to previous best
    ## Run 331 stress 0.05607347 
    ## Run 332 stress 0.0557269 
    ## Run 333 stress 0.04938625 
    ## Run 334 stress 0.04612761 
    ## ... Procrustes: rmse 9.473696e-05  max resid 0.0001920338 
    ## ... Similar to previous best
    ## Run 335 stress 0.05607347 
    ## Run 336 stress 0.05842793 
    ## Run 337 stress 0.04938617 
    ## Run 338 stress 0.05643262 
    ## Run 339 stress 0.04612757 
    ## ... Procrustes: rmse 3.788296e-05  max resid 7.69245e-05 
    ## ... Similar to previous best
    ## Run 340 stress 0.05842789 
    ## Run 341 stress 0.05572688 
    ## Run 342 stress 0.06933805 
    ## Run 343 stress 0.05593453 
    ## Run 344 stress 0.05451494 
    ## Run 345 stress 0.05695864 
    ## Run 346 stress 0.05605146 
    ## Run 347 stress 0.05855041 
    ## Run 348 stress 0.05488103 
    ## Run 349 stress 0.05728058 
    ## Run 350 stress 0.06400726 
    ## Run 351 stress 0.06205756 
    ## Run 352 stress 0.06205752 
    ## Run 353 stress 0.04851419 
    ## Run 354 stress 0.05260053 
    ## Run 355 stress 0.05572688 
    ## Run 356 stress 0.04935946 
    ## Run 357 stress 0.05607344 
    ## Run 358 stress 0.05137965 
    ## Run 359 stress 0.05605157 
    ## Run 360 stress 0.0543928 
    ## Run 361 stress 0.04935944 
    ## Run 362 stress 0.04938617 
    ## Run 363 stress 0.04938615 
    ## Run 364 stress 0.04612762 
    ## ... Procrustes: rmse 0.0001058187  max resid 0.0002139966 
    ## ... Similar to previous best
    ## Run 365 stress 0.04851419 
    ## Run 366 stress 0.05276411 
    ## Run 367 stress 0.0585503 
    ## Run 368 stress 0.05728061 
    ## Run 369 stress 0.05488104 
    ## Run 370 stress 0.05216158 
    ## Run 371 stress 0.06291426 
    ## Run 372 stress 0.06185068 
    ## Run 373 stress 0.05260048 
    ## Run 374 stress 0.05607344 
    ## Run 375 stress 0.06205754 
    ## Run 376 stress 0.06291425 
    ## Run 377 stress 0.04860422 
    ## Run 378 stress 0.05474814 
    ## Run 379 stress 0.04860422 
    ## Run 380 stress 0.04612765 
    ## ... Procrustes: rmse 0.0001427917  max resid 0.0002931882 
    ## ... Similar to previous best
    ## Run 381 stress 0.05605157 
    ## Run 382 stress 0.05474813 
    ## Run 383 stress 0.3493323 
    ## Run 384 stress 0.05260046 
    ## Run 385 stress 0.05607345 
    ## Run 386 stress 0.05451495 
    ## Run 387 stress 0.05728057 
    ## Run 388 stress 0.05855039 
    ## Run 389 stress 0.05842798 
    ## Run 390 stress 0.06205753 
    ## Run 391 stress 0.05572689 
    ## Run 392 stress 0.06088971 
    ## Run 393 stress 0.05959896 
    ## Run 394 stress 0.3582284 
    ## Run 395 stress 0.3493365 
    ## Run 396 stress 0.04612771 
    ## ... Procrustes: rmse 0.000179378  max resid 0.0003597283 
    ## ... Similar to previous best
    ## Run 397 stress 0.05137969 
    ## Run 398 stress 0.05842781 
    ## Run 399 stress 0.05488106 
    ## Run 400 stress 0.04935958 
    ## Run 401 stress 0.04851417 
    ## Run 402 stress 0.04938616 
    ## Run 403 stress 0.05433173 
    ## Run 404 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001324146  max resid 0.0002714978 
    ## ... Similar to previous best
    ## Run 405 stress 0.04851417 
    ## Run 406 stress 0.0526005 
    ## Run 407 stress 0.05728051 
    ## Run 408 stress 0.3388069 
    ## Run 409 stress 0.05728052 
    ## Run 410 stress 0.05855032 
    ## Run 411 stress 0.0557269 
    ## Run 412 stress 0.05643261 
    ## Run 413 stress 0.05438953 
    ## Run 414 stress 0.05728066 
    ## Run 415 stress 0.04938619 
    ## Run 416 stress 0.05605153 
    ## Run 417 stress 0.06088973 
    ## Run 418 stress 0.06185044 
    ## Run 419 stress 0.05728062 
    ## Run 420 stress 0.04938615 
    ## Run 421 stress 0.05605149 
    ## Run 422 stress 0.04938616 
    ## Run 423 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001114592  max resid 0.0002228166 
    ## ... Similar to previous best
    ## Run 424 stress 0.06291425 
    ## Run 425 stress 0.05607345 
    ## Run 426 stress 0.05572688 
    ## Run 427 stress 0.05216153 
    ## Run 428 stress 0.05959897 
    ## Run 429 stress 0.04612766 
    ## ... Procrustes: rmse 0.0001454608  max resid 0.0002907802 
    ## ... Similar to previous best
    ## Run 430 stress 0.05607345 
    ## Run 431 stress 0.05959902 
    ## Run 432 stress 0.05607344 
    ## Run 433 stress 0.04935944 
    ## Run 434 stress 0.3418612 
    ## Run 435 stress 0.05605147 
    ## Run 436 stress 0.04935947 
    ## Run 437 stress 0.05474818 
    ## Run 438 stress 0.06205754 
    ## Run 439 stress 0.0521614 
    ## Run 440 stress 0.0572808 
    ## Run 441 stress 0.04860423 
    ## Run 442 stress 0.05643261 
    ## Run 443 stress 0.05607344 
    ## Run 444 stress 0.06088971 
    ## Run 445 stress 0.05605147 
    ## Run 446 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001406135  max resid 0.0002858281 
    ## ... Similar to previous best
    ## Run 447 stress 0.05643261 
    ## Run 448 stress 0.05607346 
    ## Run 449 stress 0.05607344 
    ## Run 450 stress 0.04938616 
    ## Run 451 stress 0.05643262 
    ## Run 452 stress 0.05607346 
    ## Run 453 stress 0.05451496 
    ## Run 454 stress 0.04935944 
    ## Run 455 stress 0.05643266 
    ## Run 456 stress 0.05137972 
    ## Run 457 stress 0.05605147 
    ## Run 458 stress 0.04851417 
    ## Run 459 stress 0.0572808 
    ## Run 460 stress 0.05216134 
    ## Run 461 stress 0.04860422 
    ## Run 462 stress 0.05605149 
    ## Run 463 stress 0.05439282 
    ## Run 464 stress 0.05695868 
    ## Run 465 stress 0.05216147 
    ## Run 466 stress 0.04851419 
    ## Run 467 stress 0.05855031 
    ## Run 468 stress 0.05728057 
    ## Run 469 stress 0.04860422 
    ## Run 470 stress 0.0585503 
    ## Run 471 stress 0.04938617 
    ## Run 472 stress 0.06291434 
    ## Run 473 stress 0.05216148 
    ## Run 474 stress 0.04938617 
    ## Run 475 stress 0.04935944 
    ## Run 476 stress 0.05959901 
    ## Run 477 stress 0.04935943 
    ## Run 478 stress 0.05728053 
    ## Run 479 stress 0.04935949 
    ## Run 480 stress 0.04851418 
    ## Run 481 stress 0.05474814 
    ## Run 482 stress 0.05959913 
    ## Run 483 stress 0.05607345 
    ## Run 484 stress 0.0561224 
    ## Run 485 stress 0.05276412 
    ## Run 486 stress 0.05607346 
    ## Run 487 stress 0.06933839 
    ## Run 488 stress 0.04938616 
    ## Run 489 stress 0.05607348 
    ## Run 490 stress 0.05728053 
    ## Run 491 stress 0.04851417 
    ## Run 492 stress 0.04938615 
    ## Run 493 stress 0.05728062 
    ## Run 494 stress 0.05612239 
    ## Run 495 stress 0.05474827 
    ## Run 496 stress 0.04938625 
    ## Run 497 stress 0.05137967 
    ## Run 498 stress 0.0560515 
    ## Run 499 stress 0.04938617 
    ## Run 500 stress 0.04938615 
    ## *** Best solution repeated 26 times

``` r
# Stratified lakes
PD_beta_S_NMDS <- metaMDS(PD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.02693829 
    ## Run 1 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001841296  max resid 0.0003146492 
    ## ... Similar to previous best
    ## Run 2 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001235581  max resid 0.0002106226 
    ## ... Similar to previous best
    ## Run 3 stress 0.07374801 
    ## Run 4 stress 0.08382827 
    ## Run 5 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001345633  max resid 0.0002306588 
    ## ... Similar to previous best
    ## Run 6 stress 0.07374805 
    ## Run 7 stress 0.1862388 
    ## Run 8 stress 0.07547299 
    ## Run 9 stress 0.02693832 
    ## ... Procrustes: rmse 7.975641e-05  max resid 0.0001364646 
    ## ... Similar to previous best
    ## Run 10 stress 0.075473 
    ## Run 11 stress 0.1756177 
    ## Run 12 stress 0.07374812 
    ## Run 13 stress 0.07374801 
    ## Run 14 stress 0.07547299 
    ## Run 15 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001409598  max resid 0.0002415975 
    ## ... Similar to previous best
    ## Run 16 stress 0.07547299 
    ## Run 17 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001495789  max resid 0.0002560481 
    ## ... Similar to previous best
    ## Run 18 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001421381  max resid 0.0002435694 
    ## ... Similar to previous best
    ## Run 19 stress 0.0269383 
    ## ... Procrustes: rmse 0.0001057893  max resid 0.0001810431 
    ## ... Similar to previous best
    ## Run 20 stress 0.08382844 
    ## Run 21 stress 0.08382856 
    ## Run 22 stress 0.0737481 
    ## Run 23 stress 0.2714151 
    ## Run 24 stress 0.075473 
    ## Run 25 stress 0.07547299 
    ## Run 26 stress 0.07374802 
    ## Run 27 stress 0.07374805 
    ## Run 28 stress 0.07374803 
    ## Run 29 stress 0.1756177 
    ## Run 30 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001311859  max resid 0.0002248311 
    ## ... Similar to previous best
    ## Run 31 stress 0.179905 
    ## Run 32 stress 0.179905 
    ## Run 33 stress 0.073748 
    ## Run 34 stress 0.0269384 
    ## ... Procrustes: rmse 0.0002129821  max resid 0.00036538 
    ## ... Similar to previous best
    ## Run 35 stress 0.07374814 
    ## Run 36 stress 0.09391898 
    ## Run 37 stress 0.2273964 
    ## Run 38 stress 0.07374804 
    ## Run 39 stress 0.07374804 
    ## Run 40 stress 0.07374808 
    ## Run 41 stress 0.07374803 
    ## Run 42 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001644212  max resid 0.0002822155 
    ## ... Similar to previous best
    ## Run 43 stress 0.07547301 
    ## Run 44 stress 0.0269383 
    ## ... Procrustes: rmse 4.028595e-05  max resid 6.757925e-05 
    ## ... Similar to previous best
    ## Run 45 stress 0.07374805 
    ## Run 46 stress 0.075473 
    ## Run 47 stress 0.075473 
    ## Run 48 stress 0.07547299 
    ## Run 49 stress 0.07547299 
    ## Run 50 stress 0.1756177 
    ## Run 51 stress 0.02693831 
    ## ... Procrustes: rmse 5.787083e-05  max resid 9.889591e-05 
    ## ... Similar to previous best
    ## Run 52 stress 0.07374801 
    ## Run 53 stress 0.1756177 
    ## Run 54 stress 0.08382851 
    ## Run 55 stress 0.02693833 
    ## ... Procrustes: rmse 9.111271e-05  max resid 0.0001557872 
    ## ... Similar to previous best
    ## Run 56 stress 0.08382848 
    ## Run 57 stress 0.08382854 
    ## Run 58 stress 0.07547299 
    ## Run 59 stress 0.07374806 
    ## Run 60 stress 0.07547302 
    ## Run 61 stress 0.09391919 
    ## Run 62 stress 0.07374808 
    ## Run 63 stress 0.0269383 
    ## ... Procrustes: rmse 4.748346e-05  max resid 8.122811e-05 
    ## ... Similar to previous best
    ## Run 64 stress 0.07374804 
    ## Run 65 stress 0.08382833 
    ## Run 66 stress 0.07374804 
    ## Run 67 stress 0.07374804 
    ## Run 68 stress 0.07374806 
    ## Run 69 stress 0.08382867 
    ## Run 70 stress 0.07547301 
    ## Run 71 stress 0.07547299 
    ## Run 72 stress 0.07374802 
    ## Run 73 stress 0.02693835 
    ## ... Procrustes: rmse 0.000106674  max resid 0.0001824672 
    ## ... Similar to previous best
    ## Run 74 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001616832  max resid 0.0002780173 
    ## ... Similar to previous best
    ## Run 75 stress 0.0269383 
    ## ... Procrustes: rmse 4.43743e-05  max resid 7.486454e-05 
    ## ... Similar to previous best
    ## Run 76 stress 0.02693834 
    ## ... Procrustes: rmse 9.301725e-05  max resid 0.0001595562 
    ## ... Similar to previous best
    ## Run 77 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001258629  max resid 0.0002165273 
    ## ... Similar to previous best
    ## Run 78 stress 0.02693841 
    ## ... Procrustes: rmse 0.0001504249  max resid 0.0002553024 
    ## ... Similar to previous best
    ## Run 79 stress 0.02693839 
    ## ... Procrustes: rmse 0.000142116  max resid 0.0002419525 
    ## ... Similar to previous best
    ## Run 80 stress 0.1889142 
    ## Run 81 stress 0.08382863 
    ## Run 82 stress 0.02693833 
    ## ... Procrustes: rmse 7.622683e-05  max resid 0.0001288435 
    ## ... Similar to previous best
    ## Run 83 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001173587  max resid 0.0002001259 
    ## ... Similar to previous best
    ## Run 84 stress 0.02693835 
    ## ... Procrustes: rmse 6.38796e-05  max resid 0.0001046061 
    ## ... Similar to previous best
    ## Run 85 stress 0.0269383 
    ## ... Procrustes: rmse 3.755693e-05  max resid 6.459605e-05 
    ## ... Similar to previous best
    ## Run 86 stress 0.07547299 
    ## Run 87 stress 0.2650138 
    ## Run 88 stress 0.0838286 
    ## Run 89 stress 0.07374803 
    ## Run 90 stress 0.08382861 
    ## Run 91 stress 0.07374807 
    ## Run 92 stress 0.08382865 
    ## Run 93 stress 0.1862385 
    ## Run 94 stress 0.07374811 
    ## Run 95 stress 0.07374805 
    ## Run 96 stress 0.08382845 
    ## Run 97 stress 0.075473 
    ## Run 98 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001561446  max resid 0.000267001 
    ## ... Similar to previous best
    ## Run 99 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001402033  max resid 0.0002399065 
    ## ... Similar to previous best
    ## Run 100 stress 0.1756177 
    ## Run 101 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001487256  max resid 0.0002546104 
    ## ... Similar to previous best
    ## Run 102 stress 0.075473 
    ## Run 103 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001370965  max resid 0.0002343494 
    ## ... Similar to previous best
    ## Run 104 stress 0.1756177 
    ## Run 105 stress 0.07547301 
    ## Run 106 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001231464  max resid 0.0002119153 
    ## ... Similar to previous best
    ## Run 107 stress 0.0269383 
    ## ... Procrustes: rmse 9.506146e-05  max resid 0.0001631866 
    ## ... Similar to previous best
    ## Run 108 stress 0.02693834 
    ## ... Procrustes: rmse 6.675529e-05  max resid 0.0001108398 
    ## ... Similar to previous best
    ## Run 109 stress 0.08382852 
    ## Run 110 stress 0.1756177 
    ## Run 111 stress 0.07374809 
    ## Run 112 stress 0.02693828 
    ## ... New best solution
    ## ... Procrustes: rmse 2.437814e-05  max resid 4.160545e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001129327  max resid 0.0001927918 
    ## ... Similar to previous best
    ## Run 114 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001137119  max resid 0.00019384 
    ## ... Similar to previous best
    ## Run 115 stress 0.02693829 
    ## ... Procrustes: rmse 2.144502e-05  max resid 3.651137e-05 
    ## ... Similar to previous best
    ## Run 116 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001443438  max resid 0.000263196 
    ## ... Similar to previous best
    ## Run 117 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001568135  max resid 0.000269594 
    ## ... Similar to previous best
    ## Run 118 stress 0.07547301 
    ## Run 119 stress 0.08382876 
    ## Run 120 stress 0.02693831 
    ## ... Procrustes: rmse 8.762235e-05  max resid 0.0001499 
    ## ... Similar to previous best
    ## Run 121 stress 0.2273964 
    ## Run 122 stress 0.1862385 
    ## Run 123 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001515632  max resid 0.0002569387 
    ## ... Similar to previous best
    ## Run 124 stress 0.02693829 
    ## ... Procrustes: rmse 2.985033e-05  max resid 5.113728e-05 
    ## ... Similar to previous best
    ## Run 125 stress 0.02693829 
    ## ... Procrustes: rmse 2.522218e-05  max resid 4.337169e-05 
    ## ... Similar to previous best
    ## Run 126 stress 0.0269383 
    ## ... Procrustes: rmse 7.081535e-05  max resid 0.0001234907 
    ## ... Similar to previous best
    ## Run 127 stress 0.02693833 
    ## ... Procrustes: rmse 9.765582e-05  max resid 0.0001656604 
    ## ... Similar to previous best
    ## Run 128 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001642785  max resid 0.0002811498 
    ## ... Similar to previous best
    ## Run 129 stress 0.1756177 
    ## Run 130 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001429046  max resid 0.0002451635 
    ## ... Similar to previous best
    ## Run 131 stress 0.08382823 
    ## Run 132 stress 0.02693841 
    ## ... Procrustes: rmse 0.0001749175  max resid 0.0002978724 
    ## ... Similar to previous best
    ## Run 133 stress 0.1756177 
    ## Run 134 stress 0.075473 
    ## Run 135 stress 0.02693828 
    ## ... New best solution
    ## ... Procrustes: rmse 7.246926e-06  max resid 1.246408e-05 
    ## ... Similar to previous best
    ## Run 136 stress 0.07374807 
    ## Run 137 stress 0.07547302 
    ## Run 138 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001543602  max resid 0.000263633 
    ## ... Similar to previous best
    ## Run 139 stress 0.02693831 
    ## ... Procrustes: rmse 8.674195e-05  max resid 0.0001259132 
    ## ... Similar to previous best
    ## Run 140 stress 0.07374803 
    ## Run 141 stress 0.02693828 
    ## ... Procrustes: rmse 5.585583e-06  max resid 8.739269e-06 
    ## ... Similar to previous best
    ## Run 142 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001633479  max resid 0.0002808392 
    ## ... Similar to previous best
    ## Run 143 stress 0.02693832 
    ## ... Procrustes: rmse 9.646037e-05  max resid 0.0001654507 
    ## ... Similar to previous best
    ## Run 144 stress 0.179905 
    ## Run 145 stress 0.07374807 
    ## Run 146 stress 0.07547299 
    ## Run 147 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001874869  max resid 0.0003205452 
    ## ... Similar to previous best
    ## Run 148 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001804635  max resid 0.0003081553 
    ## ... Similar to previous best
    ## Run 149 stress 0.02693831 
    ## ... Procrustes: rmse 9.257868e-05  max resid 0.0001584969 
    ## ... Similar to previous best
    ## Run 150 stress 0.1756177 
    ## Run 151 stress 0.07547299 
    ## Run 152 stress 0.07374803 
    ## Run 153 stress 0.07547299 
    ## Run 154 stress 0.07374809 
    ## Run 155 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001206714  max resid 0.0002062125 
    ## ... Similar to previous best
    ## Run 156 stress 0.02693842 
    ## ... Procrustes: rmse 0.0002029711  max resid 0.0003474502 
    ## ... Similar to previous best
    ## Run 157 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001371758  max resid 0.0002345929 
    ## ... Similar to previous best
    ## Run 158 stress 0.1756177 
    ## Run 159 stress 0.08382851 
    ## Run 160 stress 0.07374801 
    ## Run 161 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001403814  max resid 0.0002226731 
    ## ... Similar to previous best
    ## Run 162 stress 0.08382852 
    ## Run 163 stress 0.07374804 
    ## Run 164 stress 0.07374802 
    ## Run 165 stress 0.08382863 
    ## Run 166 stress 0.07374808 
    ## Run 167 stress 0.1889141 
    ## Run 168 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001392091  max resid 0.000238708 
    ## ... Similar to previous best
    ## Run 169 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001275779  max resid 0.0002184117 
    ## ... Similar to previous best
    ## Run 170 stress 0.02693831 
    ## ... Procrustes: rmse 8.53095e-05  max resid 0.0001461148 
    ## ... Similar to previous best
    ## Run 171 stress 0.07547299 
    ## Run 172 stress 0.08382868 
    ## Run 173 stress 0.07374802 
    ## Run 174 stress 0.07374812 
    ## Run 175 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001113311  max resid 0.0001910815 
    ## ... Similar to previous best
    ## Run 176 stress 0.02693832 
    ## ... Procrustes: rmse 9.460544e-05  max resid 0.0001609134 
    ## ... Similar to previous best
    ## Run 177 stress 0.07374809 
    ## Run 178 stress 0.02693828 
    ## ... Procrustes: rmse 8.278113e-06  max resid 1.1434e-05 
    ## ... Similar to previous best
    ## Run 179 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001218797  max resid 0.0002087861 
    ## ... Similar to previous best
    ## Run 180 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001309488  max resid 0.0002243216 
    ## ... Similar to previous best
    ## Run 181 stress 0.0838284 
    ## Run 182 stress 0.07547302 
    ## Run 183 stress 0.02693829 
    ## ... Procrustes: rmse 2.767368e-05  max resid 4.781319e-05 
    ## ... Similar to previous best
    ## Run 184 stress 0.3083098 
    ## Run 185 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001615059  max resid 0.0002754608 
    ## ... Similar to previous best
    ## Run 186 stress 0.073748 
    ## Run 187 stress 0.1889141 
    ## Run 188 stress 0.02693829 
    ## ... Procrustes: rmse 4.905673e-05  max resid 8.393995e-05 
    ## ... Similar to previous best
    ## Run 189 stress 0.02693833 
    ## ... Procrustes: rmse 8.470788e-05  max resid 0.0001462722 
    ## ... Similar to previous best
    ## Run 190 stress 0.02693831 
    ## ... Procrustes: rmse 7.680625e-05  max resid 0.0001303736 
    ## ... Similar to previous best
    ## Run 191 stress 0.1862385 
    ## Run 192 stress 0.0737481 
    ## Run 193 stress 0.07374802 
    ## Run 194 stress 0.0269383 
    ## ... Procrustes: rmse 7.708855e-05  max resid 0.0001350213 
    ## ... Similar to previous best
    ## Run 195 stress 0.1756177 
    ## Run 196 stress 0.07547301 
    ## Run 197 stress 0.07374808 
    ## Run 198 stress 0.07374803 
    ## Run 199 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001595039  max resid 0.0002731493 
    ## ... Similar to previous best
    ## Run 200 stress 0.1756177 
    ## Run 201 stress 0.07374802 
    ## Run 202 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001410319  max resid 0.0002410058 
    ## ... Similar to previous best
    ## Run 203 stress 0.02693828 
    ## ... Procrustes: rmse 1.172573e-06  max resid 2.088967e-06 
    ## ... Similar to previous best
    ## Run 204 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001366284  max resid 0.0002340515 
    ## ... Similar to previous best
    ## Run 205 stress 0.073748 
    ## Run 206 stress 0.2650138 
    ## Run 207 stress 0.075473 
    ## Run 208 stress 0.179905 
    ## Run 209 stress 0.07547299 
    ## Run 210 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001501844  max resid 0.0002595347 
    ## ... Similar to previous best
    ## Run 211 stress 0.07547299 
    ## Run 212 stress 0.02693835 
    ## ... Procrustes: rmse 0.000142259  max resid 0.0002434546 
    ## ... Similar to previous best
    ## Run 213 stress 0.07547299 
    ## Run 214 stress 0.07547299 
    ## Run 215 stress 0.07547301 
    ## Run 216 stress 0.02693832 
    ## ... Procrustes: rmse 9.871937e-05  max resid 0.0001690615 
    ## ... Similar to previous best
    ## Run 217 stress 0.2273964 
    ## Run 218 stress 0.08382825 
    ## Run 219 stress 0.07547299 
    ## Run 220 stress 0.07374801 
    ## Run 221 stress 0.186239 
    ## Run 222 stress 0.0838283 
    ## Run 223 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001464184  max resid 0.0002482939 
    ## ... Similar to previous best
    ## Run 224 stress 0.08382832 
    ## Run 225 stress 0.07374801 
    ## Run 226 stress 0.07374802 
    ## Run 227 stress 0.075473 
    ## Run 228 stress 0.02693831 
    ## ... Procrustes: rmse 9.324507e-05  max resid 0.0001597369 
    ## ... Similar to previous best
    ## Run 229 stress 0.02693842 
    ## ... Procrustes: rmse 0.000199636  max resid 0.0003415141 
    ## ... Similar to previous best
    ## Run 230 stress 0.1889141 
    ## Run 231 stress 0.1756177 
    ## Run 232 stress 0.07374802 
    ## Run 233 stress 0.1889141 
    ## Run 234 stress 0.07374809 
    ## Run 235 stress 0.02693829 
    ## ... Procrustes: rmse 5.179237e-05  max resid 8.897975e-05 
    ## ... Similar to previous best
    ## Run 236 stress 0.0838283 
    ## Run 237 stress 0.08382844 
    ## Run 238 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001750174  max resid 0.0002999371 
    ## ... Similar to previous best
    ## Run 239 stress 0.0269383 
    ## ... Procrustes: rmse 5.54592e-05  max resid 9.453268e-05 
    ## ... Similar to previous best
    ## Run 240 stress 0.02693831 
    ## ... Procrustes: rmse 7.323841e-05  max resid 0.0001241552 
    ## ... Similar to previous best
    ## Run 241 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001065856  max resid 0.000182098 
    ## ... Similar to previous best
    ## Run 242 stress 0.07374801 
    ## Run 243 stress 0.07374808 
    ## Run 244 stress 0.07547301 
    ## Run 245 stress 0.07374806 
    ## Run 246 stress 0.1889142 
    ## Run 247 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001240592  max resid 0.0001908648 
    ## ... Similar to previous best
    ## Run 248 stress 0.1756177 
    ## Run 249 stress 0.179905 
    ## Run 250 stress 0.07547299 
    ## Run 251 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001232278  max resid 0.0002114673 
    ## ... Similar to previous best
    ## Run 252 stress 0.0269383 
    ## ... Procrustes: rmse 6.245071e-05  max resid 0.0001069401 
    ## ... Similar to previous best
    ## Run 253 stress 0.07374807 
    ## Run 254 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001012651  max resid 0.0001743163 
    ## ... Similar to previous best
    ## Run 255 stress 0.1756177 
    ## Run 256 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001314721  max resid 0.0002255257 
    ## ... Similar to previous best
    ## Run 257 stress 0.02693829 
    ## ... Procrustes: rmse 1.746597e-05  max resid 2.993409e-05 
    ## ... Similar to previous best
    ## Run 258 stress 0.02693829 
    ## ... Procrustes: rmse 5.518674e-05  max resid 9.447231e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.08382861 
    ## Run 260 stress 0.07374802 
    ## Run 261 stress 0.1889141 
    ## Run 262 stress 0.07374806 
    ## Run 263 stress 0.07374803 
    ## Run 264 stress 0.02693832 
    ## ... Procrustes: rmse 9.973301e-05  max resid 0.0001704159 
    ## ... Similar to previous best
    ## Run 265 stress 0.07547301 
    ## Run 266 stress 0.02693844 
    ## ... Procrustes: rmse 0.0001883806  max resid 0.0003246449 
    ## ... Similar to previous best
    ## Run 267 stress 0.08382829 
    ## Run 268 stress 0.1756177 
    ## Run 269 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001622571  max resid 0.0002782823 
    ## ... Similar to previous best
    ## Run 270 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001688623  max resid 0.0002890873 
    ## ... Similar to previous best
    ## Run 271 stress 0.02693832 
    ## ... Procrustes: rmse 9.128763e-05  max resid 0.0001549591 
    ## ... Similar to previous best
    ## Run 272 stress 0.02693841 
    ## ... Procrustes: rmse 0.0001886782  max resid 0.0003238938 
    ## ... Similar to previous best
    ## Run 273 stress 0.1889141 
    ## Run 274 stress 0.07374801 
    ## Run 275 stress 0.07374802 
    ## Run 276 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001266243  max resid 0.0002165399 
    ## ... Similar to previous best
    ## Run 277 stress 0.02693831 
    ## ... Procrustes: rmse 8.033033e-05  max resid 0.0001369437 
    ## ... Similar to previous best
    ## Run 278 stress 0.07374802 
    ## Run 279 stress 0.1799052 
    ## Run 280 stress 0.179905 
    ## Run 281 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001153731  max resid 0.0001974246 
    ## ... Similar to previous best
    ## Run 282 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001764538  max resid 0.0003016708 
    ## ... Similar to previous best
    ## Run 283 stress 0.07374805 
    ## Run 284 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001454775  max resid 0.0002497299 
    ## ... Similar to previous best
    ## Run 285 stress 0.179905 
    ## Run 286 stress 0.02693828 
    ## ... Procrustes: rmse 1.007622e-05  max resid 1.681175e-05 
    ## ... Similar to previous best
    ## Run 287 stress 0.07374803 
    ## Run 288 stress 0.02693832 
    ## ... Procrustes: rmse 6.496786e-05  max resid 0.0001115397 
    ## ... Similar to previous best
    ## Run 289 stress 0.07374802 
    ## Run 290 stress 0.1862389 
    ## Run 291 stress 0.02693842 
    ## ... Procrustes: rmse 0.0002004224  max resid 0.0003441307 
    ## ... Similar to previous best
    ## Run 292 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001814463  max resid 0.0003118237 
    ## ... Similar to previous best
    ## Run 293 stress 0.08382836 
    ## Run 294 stress 0.07374803 
    ## Run 295 stress 0.07547301 
    ## Run 296 stress 0.07547299 
    ## Run 297 stress 0.08382832 
    ## Run 298 stress 0.0269383 
    ## ... Procrustes: rmse 7.507357e-05  max resid 0.000128532 
    ## ... Similar to previous best
    ## Run 299 stress 0.07374802 
    ## Run 300 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001068527  max resid 0.0001832793 
    ## ... Similar to previous best
    ## Run 301 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001250582  max resid 0.0002195904 
    ## ... Similar to previous best
    ## Run 302 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001328084  max resid 0.0002317805 
    ## ... Similar to previous best
    ## Run 303 stress 0.07374801 
    ## Run 304 stress 0.08382846 
    ## Run 305 stress 0.07374801 
    ## Run 306 stress 0.1900871 
    ## Run 307 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001126449  max resid 0.0001940475 
    ## ... Similar to previous best
    ## Run 308 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001371658  max resid 0.0002348872 
    ## ... Similar to previous best
    ## Run 309 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001090163  max resid 0.0001861178 
    ## ... Similar to previous best
    ## Run 310 stress 0.07547299 
    ## Run 311 stress 0.07374803 
    ## Run 312 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001427722  max resid 0.0002445044 
    ## ... Similar to previous best
    ## Run 313 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001401494  max resid 0.0002405634 
    ## ... Similar to previous best
    ## Run 314 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001420947  max resid 0.0002446878 
    ## ... Similar to previous best
    ## Run 315 stress 0.07374801 
    ## Run 316 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001306497  max resid 0.0002236551 
    ## ... Similar to previous best
    ## Run 317 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001174716  max resid 0.000193226 
    ## ... Similar to previous best
    ## Run 318 stress 0.07547299 
    ## Run 319 stress 0.07374803 
    ## Run 320 stress 0.02693829 
    ## ... Procrustes: rmse 2.342632e-05  max resid 4.011081e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.02693831 
    ## ... Procrustes: rmse 8.61914e-05  max resid 0.0001511586 
    ## ... Similar to previous best
    ## Run 322 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001390342  max resid 0.0002379513 
    ## ... Similar to previous best
    ## Run 323 stress 0.1900873 
    ## Run 324 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001607197  max resid 0.0002761468 
    ## ... Similar to previous best
    ## Run 325 stress 0.07374803 
    ## Run 326 stress 0.02693842 
    ## ... Procrustes: rmse 0.000202899  max resid 0.0003483078 
    ## ... Similar to previous best
    ## Run 327 stress 0.08382867 
    ## Run 328 stress 0.07374803 
    ## Run 329 stress 0.07374803 
    ## Run 330 stress 0.07374801 
    ## Run 331 stress 0.07374804 
    ## Run 332 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001017959  max resid 0.0001744055 
    ## ... Similar to previous best
    ## Run 333 stress 0.07374809 
    ## Run 334 stress 0.02693831 
    ## ... Procrustes: rmse 8.438466e-05  max resid 0.0001448223 
    ## ... Similar to previous best
    ## Run 335 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001656923  max resid 0.0002836768 
    ## ... Similar to previous best
    ## Run 336 stress 0.179905 
    ## Run 337 stress 0.02693829 
    ## ... Procrustes: rmse 4.374067e-05  max resid 7.489939e-05 
    ## ... Similar to previous best
    ## Run 338 stress 0.02693837 
    ## ... Procrustes: rmse 0.000159803  max resid 0.0002729326 
    ## ... Similar to previous best
    ## Run 339 stress 0.07374801 
    ## Run 340 stress 0.073748 
    ## Run 341 stress 0.3083098 
    ## Run 342 stress 0.08382862 
    ## Run 343 stress 0.02693835 
    ## ... Procrustes: rmse 0.000143318  max resid 0.0002462496 
    ## ... Similar to previous best
    ## Run 344 stress 0.02693828 
    ## ... Procrustes: rmse 8.131456e-06  max resid 1.390263e-05 
    ## ... Similar to previous best
    ## Run 345 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001022103  max resid 0.0001749185 
    ## ... Similar to previous best
    ## Run 346 stress 0.075473 
    ## Run 347 stress 0.179905 
    ## Run 348 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001262104  max resid 0.0002167332 
    ## ... Similar to previous best
    ## Run 349 stress 0.0269383 
    ## ... Procrustes: rmse 6.965163e-05  max resid 0.0001193171 
    ## ... Similar to previous best
    ## Run 350 stress 0.08382832 
    ## Run 351 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001294062  max resid 0.0002220953 
    ## ... Similar to previous best
    ## Run 352 stress 0.02693829 
    ## ... Procrustes: rmse 1.733775e-05  max resid 2.306395e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.07374802 
    ## Run 354 stress 0.179905 
    ## Run 355 stress 0.0269383 
    ## ... Procrustes: rmse 7.538053e-05  max resid 0.0001095311 
    ## ... Similar to previous best
    ## Run 356 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001196831  max resid 0.0002045761 
    ## ... Similar to previous best
    ## Run 357 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001555417  max resid 0.0002642775 
    ## ... Similar to previous best
    ## Run 358 stress 0.02693831 
    ## ... Procrustes: rmse 8.343493e-05  max resid 0.0001421736 
    ## ... Similar to previous best
    ## Run 359 stress 0.02693829 
    ## ... Procrustes: rmse 2.613728e-05  max resid 4.454075e-05 
    ## ... Similar to previous best
    ## Run 360 stress 0.1799052 
    ## Run 361 stress 0.08382824 
    ## Run 362 stress 0.179905 
    ## Run 363 stress 0.08382842 
    ## Run 364 stress 0.1889141 
    ## Run 365 stress 0.1900871 
    ## Run 366 stress 0.075473 
    ## Run 367 stress 0.02693834 
    ## ... Procrustes: rmse 0.000128587  max resid 0.0002201932 
    ## ... Similar to previous best
    ## Run 368 stress 0.179905 
    ## Run 369 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001438704  max resid 0.0002453168 
    ## ... Similar to previous best
    ## Run 370 stress 0.07374805 
    ## Run 371 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001042163  max resid 0.0001784844 
    ## ... Similar to previous best
    ## Run 372 stress 0.0269383 
    ## ... Procrustes: rmse 6.933806e-05  max resid 0.000119101 
    ## ... Similar to previous best
    ## Run 373 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001432848  max resid 0.0002446605 
    ## ... Similar to previous best
    ## Run 374 stress 0.073748 
    ## Run 375 stress 0.07547299 
    ## Run 376 stress 0.1756177 
    ## Run 377 stress 0.1889142 
    ## Run 378 stress 0.07547299 
    ## Run 379 stress 0.179905 
    ## Run 380 stress 0.2398825 
    ## Run 381 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001841764  max resid 0.0003143577 
    ## ... Similar to previous best
    ## Run 382 stress 0.07374801 
    ## Run 383 stress 0.02693828 
    ## ... Procrustes: rmse 3.999175e-06  max resid 6.570145e-06 
    ## ... Similar to previous best
    ## Run 384 stress 0.07547299 
    ## Run 385 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001136276  max resid 0.000193977 
    ## ... Similar to previous best
    ## Run 386 stress 0.07374808 
    ## Run 387 stress 0.02693829 
    ## ... Procrustes: rmse 3.412727e-05  max resid 5.841717e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.07374807 
    ## Run 389 stress 0.1756177 
    ## Run 390 stress 0.02693831 
    ## ... Procrustes: rmse 7.744465e-05  max resid 0.000132986 
    ## ... Similar to previous best
    ## Run 391 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001713075  max resid 0.0002926857 
    ## ... Similar to previous best
    ## Run 392 stress 0.07374802 
    ## Run 393 stress 0.1862389 
    ## Run 394 stress 0.09391913 
    ## Run 395 stress 0.07547299 
    ## Run 396 stress 0.02693841 
    ## ... Procrustes: rmse 0.0001769306  max resid 0.0003008322 
    ## ... Similar to previous best
    ## Run 397 stress 0.08382837 
    ## Run 398 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001311851  max resid 0.0002251031 
    ## ... Similar to previous best
    ## Run 399 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001282158  max resid 0.000220148 
    ## ... Similar to previous best
    ## Run 400 stress 0.0269383 
    ## ... Procrustes: rmse 6.193685e-05  max resid 0.0001060309 
    ## ... Similar to previous best
    ## Run 401 stress 0.1889142 
    ## Run 402 stress 0.07547299 
    ## Run 403 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001244539  max resid 0.000213543 
    ## ... Similar to previous best
    ## Run 404 stress 0.0269383 
    ## ... Procrustes: rmse 6.046477e-05  max resid 0.0001039721 
    ## ... Similar to previous best
    ## Run 405 stress 0.07374801 
    ## Run 406 stress 0.07374805 
    ## Run 407 stress 0.02693834 
    ## ... Procrustes: rmse 0.000134856  max resid 0.0002305942 
    ## ... Similar to previous best
    ## Run 408 stress 0.07547299 
    ## Run 409 stress 0.02693833 
    ## ... Procrustes: rmse 7.738889e-05  max resid 0.000133704 
    ## ... Similar to previous best
    ## Run 410 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001796042  max resid 0.000306907 
    ## ... Similar to previous best
    ## Run 411 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001602518  max resid 0.0002739117 
    ## ... Similar to previous best
    ## Run 412 stress 0.07547302 
    ## Run 413 stress 0.1799051 
    ## Run 414 stress 0.179905 
    ## Run 415 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001831991  max resid 0.0003128326 
    ## ... Similar to previous best
    ## Run 416 stress 0.07547299 
    ## Run 417 stress 0.02693829 
    ## ... Procrustes: rmse 1.943574e-05  max resid 3.336584e-05 
    ## ... Similar to previous best
    ## Run 418 stress 0.07374803 
    ## Run 419 stress 0.179905 
    ## Run 420 stress 0.1889141 
    ## Run 421 stress 0.07374802 
    ## Run 422 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001558952  max resid 0.0002676892 
    ## ... Similar to previous best
    ## Run 423 stress 0.073748 
    ## Run 424 stress 0.08382839 
    ## Run 425 stress 0.07374808 
    ## Run 426 stress 0.0269383 
    ## ... Procrustes: rmse 4.49742e-05  max resid 7.762756e-05 
    ## ... Similar to previous best
    ## Run 427 stress 0.1862383 
    ## Run 428 stress 0.02693844 
    ## ... Procrustes: rmse 0.0001953132  max resid 0.0003319367 
    ## ... Similar to previous best
    ## Run 429 stress 0.0269383 
    ## ... Procrustes: rmse 6.596131e-05  max resid 0.0001125783 
    ## ... Similar to previous best
    ## Run 430 stress 0.1756177 
    ## Run 431 stress 0.02693841 
    ## ... Procrustes: rmse 0.0001927523  max resid 0.0003293134 
    ## ... Similar to previous best
    ## Run 432 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001662766  max resid 0.0002816253 
    ## ... Similar to previous best
    ## Run 433 stress 0.07374809 
    ## Run 434 stress 0.08382867 
    ## Run 435 stress 0.08382855 
    ## Run 436 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001133835  max resid 0.0001947389 
    ## ... Similar to previous best
    ## Run 437 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001103393  max resid 0.0001826377 
    ## ... Similar to previous best
    ## Run 438 stress 0.08382854 
    ## Run 439 stress 0.075473 
    ## Run 440 stress 0.07374811 
    ## Run 441 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001790667  max resid 0.0003073636 
    ## ... Similar to previous best
    ## Run 442 stress 0.07547299 
    ## Run 443 stress 0.02693829 
    ## ... Procrustes: rmse 5.026827e-05  max resid 8.633224e-05 
    ## ... Similar to previous best
    ## Run 444 stress 0.07547301 
    ## Run 445 stress 0.07547299 
    ## Run 446 stress 0.02693836 
    ## ... Procrustes: rmse 0.00011227  max resid 0.000193689 
    ## ... Similar to previous best
    ## Run 447 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001071781  max resid 0.0001835867 
    ## ... Similar to previous best
    ## Run 448 stress 0.0939191 
    ## Run 449 stress 0.073748 
    ## Run 450 stress 0.07374801 
    ## Run 451 stress 0.08382821 
    ## Run 452 stress 0.1756177 
    ## Run 453 stress 0.02693828 
    ## ... Procrustes: rmse 8.02565e-06  max resid 1.383824e-05 
    ## ... Similar to previous best
    ## Run 454 stress 0.07374803 
    ## Run 455 stress 0.07374806 
    ## Run 456 stress 0.073748 
    ## Run 457 stress 0.07374801 
    ## Run 458 stress 0.07547299 
    ## Run 459 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001018715  max resid 0.000175239 
    ## ... Similar to previous best
    ## Run 460 stress 0.07547299 
    ## Run 461 stress 0.02693838 
    ## ... Procrustes: rmse 0.000171896  max resid 0.0002947687 
    ## ... Similar to previous best
    ## Run 462 stress 0.02693831 
    ## ... Procrustes: rmse 7.978794e-05  max resid 0.0001427103 
    ## ... Similar to previous best
    ## Run 463 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001276471  max resid 0.0002192945 
    ## ... Similar to previous best
    ## Run 464 stress 0.02693829 
    ## ... Procrustes: rmse 5.045065e-05  max resid 8.61643e-05 
    ## ... Similar to previous best
    ## Run 465 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001588185  max resid 0.0002719038 
    ## ... Similar to previous best
    ## Run 466 stress 0.02693829 
    ## ... Procrustes: rmse 1.088106e-05  max resid 1.601135e-05 
    ## ... Similar to previous best
    ## Run 467 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001556425  max resid 0.0002661915 
    ## ... Similar to previous best
    ## Run 468 stress 0.1862389 
    ## Run 469 stress 0.02693835 
    ## ... Procrustes: rmse 0.000141346  max resid 0.0002414919 
    ## ... Similar to previous best
    ## Run 470 stress 0.02693842 
    ## ... Procrustes: rmse 0.0002023736  max resid 0.000346006 
    ## ... Similar to previous best
    ## Run 471 stress 0.08382824 
    ## Run 472 stress 0.02693831 
    ## ... Procrustes: rmse 7.323263e-05  max resid 0.0001237945 
    ## ... Similar to previous best
    ## Run 473 stress 0.02693831 
    ## ... Procrustes: rmse 9.667126e-05  max resid 0.0001657438 
    ## ... Similar to previous best
    ## Run 474 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001107571  max resid 0.0001890925 
    ## ... Similar to previous best
    ## Run 475 stress 0.1862382 
    ## Run 476 stress 0.073748 
    ## Run 477 stress 0.07374802 
    ## Run 478 stress 0.07374804 
    ## Run 479 stress 0.179905 
    ## Run 480 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001331585  max resid 0.0002279042 
    ## ... Similar to previous best
    ## Run 481 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001258633  max resid 0.0002159418 
    ## ... Similar to previous best
    ## Run 482 stress 0.073748 
    ## Run 483 stress 0.02693831 
    ## ... Procrustes: rmse 9.14179e-05  max resid 0.0001559266 
    ## ... Similar to previous best
    ## Run 484 stress 0.07374803 
    ## Run 485 stress 0.075473 
    ## Run 486 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001200291  max resid 0.0002060047 
    ## ... Similar to previous best
    ## Run 487 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001168581  max resid 0.0001972732 
    ## ... Similar to previous best
    ## Run 488 stress 0.07547301 
    ## Run 489 stress 0.075473 
    ## Run 490 stress 0.07374801 
    ## Run 491 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001463185  max resid 0.0002490749 
    ## ... Similar to previous best
    ## Run 492 stress 0.075473 
    ## Run 493 stress 0.08382842 
    ## Run 494 stress 0.07547299 
    ## Run 495 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001109345  max resid 0.0001901051 
    ## ... Similar to previous best
    ## Run 496 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001550141  max resid 0.0002655794 
    ## ... Similar to previous best
    ## Run 497 stress 0.08382848 
    ## Run 498 stress 0.08382853 
    ## Run 499 stress 0.0269383 
    ## ... Procrustes: rmse 7.440991e-05  max resid 0.0001273616 
    ## ... Similar to previous best
    ## Run 500 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001278818  max resid 0.0002191343 
    ## ... Similar to previous best
    ## *** Best solution repeated 157 times

``` r
### Environmental
# Surveyed sites 
PD_beta_env_NMDS <- metaMDS(PD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06330868 
    ## Run 1 stress 0.06402624 
    ## Run 2 stress 0.09290987 
    ## Run 3 stress 0.06402619 
    ## Run 4 stress 0.0633087 
    ## ... Procrustes: rmse 6.113799e-05  max resid 0.000106911 
    ## ... Similar to previous best
    ## Run 5 stress 0.06402636 
    ## Run 6 stress 0.09245946 
    ## Run 7 stress 0.06330871 
    ## ... Procrustes: rmse 7.472764e-05  max resid 0.0001659865 
    ## ... Similar to previous best
    ## Run 8 stress 0.06402617 
    ## Run 9 stress 0.0633087 
    ## ... Procrustes: rmse 6.730122e-05  max resid 0.0001505432 
    ## ... Similar to previous best
    ## Run 10 stress 0.06402616 
    ## Run 11 stress 0.06330871 
    ## ... Procrustes: rmse 8.978626e-05  max resid 0.0002009966 
    ## ... Similar to previous best
    ## Run 12 stress 0.06463315 
    ## Run 13 stress 0.09232143 
    ## Run 14 stress 0.09217179 
    ## Run 15 stress 0.06330868 
    ## ... Procrustes: rmse 1.766807e-05  max resid 3.089886e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.06402628 
    ## Run 17 stress 0.06402631 
    ## Run 18 stress 0.06463316 
    ## Run 19 stress 0.06463313 
    ## Run 20 stress 0.06463315 
    ## Run 21 stress 0.06330869 
    ## ... Procrustes: rmse 1.473398e-05  max resid 5.102477e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.06330869 
    ## ... Procrustes: rmse 2.886783e-05  max resid 6.55396e-05 
    ## ... Similar to previous best
    ## Run 23 stress 0.06330879 
    ## ... Procrustes: rmse 0.0001546765  max resid 0.0003440786 
    ## ... Similar to previous best
    ## Run 24 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 1.268722e-05  max resid 2.814209e-05 
    ## ... Similar to previous best
    ## Run 25 stress 0.06330877 
    ## ... Procrustes: rmse 0.0001153218  max resid 0.0002568696 
    ## ... Similar to previous best
    ## Run 26 stress 0.06330868 
    ## ... Procrustes: rmse 1.567815e-05  max resid 5.087582e-05 
    ## ... Similar to previous best
    ## Run 27 stress 0.06463314 
    ## Run 28 stress 0.06402625 
    ## Run 29 stress 0.06402626 
    ## Run 30 stress 0.06402628 
    ## Run 31 stress 0.06402617 
    ## Run 32 stress 0.06330868 
    ## ... Procrustes: rmse 1.400704e-05  max resid 3.171552e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.09533676 
    ## Run 34 stress 0.06330869 
    ## ... Procrustes: rmse 4.43263e-05  max resid 9.838243e-05 
    ## ... Similar to previous best
    ## Run 35 stress 0.08089064 
    ## Run 36 stress 0.06330869 
    ## ... Procrustes: rmse 3.603801e-05  max resid 9.600761e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.09428577 
    ## Run 38 stress 0.0640263 
    ## Run 39 stress 0.06463314 
    ## Run 40 stress 0.06402623 
    ## Run 41 stress 0.09485822 
    ## Run 42 stress 0.100585 
    ## Run 43 stress 0.06330874 
    ## ... Procrustes: rmse 0.0001050863  max resid 0.000235895 
    ## ... Similar to previous best
    ## Run 44 stress 0.09106974 
    ## Run 45 stress 0.06330868 
    ## ... Procrustes: rmse 1.426717e-05  max resid 5.203642e-05 
    ## ... Similar to previous best
    ## Run 46 stress 0.06330868 
    ## ... Procrustes: rmse 2.002053e-05  max resid 4.463628e-05 
    ## ... Similar to previous best
    ## Run 47 stress 0.09233456 
    ## Run 48 stress 0.06402619 
    ## Run 49 stress 0.06402618 
    ## Run 50 stress 0.0640262 
    ## Run 51 stress 0.06330868 
    ## ... Procrustes: rmse 2.607379e-05  max resid 5.79467e-05 
    ## ... Similar to previous best
    ## Run 52 stress 0.06330868 
    ## ... Procrustes: rmse 1.70794e-06  max resid 4.541871e-06 
    ## ... Similar to previous best
    ## Run 53 stress 0.08233733 
    ## Run 54 stress 0.06463313 
    ## Run 55 stress 0.06402615 
    ## Run 56 stress 0.09270094 
    ## Run 57 stress 0.0633087 
    ## ... Procrustes: rmse 4.944575e-05  max resid 0.0001007125 
    ## ... Similar to previous best
    ## Run 58 stress 0.06330868 
    ## ... Procrustes: rmse 1.661116e-05  max resid 3.831109e-05 
    ## ... Similar to previous best
    ## Run 59 stress 0.06463313 
    ## Run 60 stress 0.06402624 
    ## Run 61 stress 0.09995795 
    ## Run 62 stress 0.06330868 
    ## ... Procrustes: rmse 2.178766e-05  max resid 4.036652e-05 
    ## ... Similar to previous best
    ## Run 63 stress 0.08126162 
    ## Run 64 stress 0.08126172 
    ## Run 65 stress 0.06463313 
    ## Run 66 stress 0.06463316 
    ## Run 67 stress 0.06330869 
    ## ... Procrustes: rmse 2.435695e-05  max resid 7.500453e-05 
    ## ... Similar to previous best
    ## Run 68 stress 0.0633087 
    ## ... Procrustes: rmse 6.630908e-05  max resid 0.0001485118 
    ## ... Similar to previous best
    ## Run 69 stress 0.06402617 
    ## Run 70 stress 0.08126192 
    ## Run 71 stress 0.06330876 
    ## ... Procrustes: rmse 0.0001053285  max resid 0.0002326576 
    ## ... Similar to previous best
    ## Run 72 stress 0.06463318 
    ## Run 73 stress 0.06463314 
    ## Run 74 stress 0.09773519 
    ## Run 75 stress 0.08089053 
    ## Run 76 stress 0.06463315 
    ## Run 77 stress 0.06402616 
    ## Run 78 stress 0.06330874 
    ## ... Procrustes: rmse 0.0001082996  max resid 0.0002424132 
    ## ... Similar to previous best
    ## Run 79 stress 0.06330869 
    ## ... Procrustes: rmse 3.346093e-05  max resid 0.0001170366 
    ## ... Similar to previous best
    ## Run 80 stress 0.0640262 
    ## Run 81 stress 0.08357056 
    ## Run 82 stress 0.06463313 
    ## Run 83 stress 0.06402628 
    ## Run 84 stress 0.0633087 
    ## ... Procrustes: rmse 5.391566e-05  max resid 0.0001239805 
    ## ... Similar to previous best
    ## Run 85 stress 0.06463318 
    ## Run 86 stress 0.06402621 
    ## Run 87 stress 0.08089055 
    ## Run 88 stress 0.06402616 
    ## Run 89 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 5.125893e-06  max resid 1.121767e-05 
    ## ... Similar to previous best
    ## Run 90 stress 0.06330868 
    ## ... Procrustes: rmse 2.252494e-06  max resid 4.515765e-06 
    ## ... Similar to previous best
    ## Run 91 stress 0.06402618 
    ## Run 92 stress 0.06402616 
    ## Run 93 stress 0.0808907 
    ## Run 94 stress 0.06330868 
    ## ... Procrustes: rmse 1.527463e-05  max resid 3.36373e-05 
    ## ... Similar to previous best
    ## Run 95 stress 0.06330875 
    ## ... Procrustes: rmse 0.000110721  max resid 0.000239365 
    ## ... Similar to previous best
    ## Run 96 stress 0.06402627 
    ## Run 97 stress 0.06330868 
    ## ... Procrustes: rmse 2.093766e-05  max resid 4.661733e-05 
    ## ... Similar to previous best
    ## Run 98 stress 0.06330868 
    ## ... Procrustes: rmse 2.153759e-05  max resid 4.578124e-05 
    ## ... Similar to previous best
    ## Run 99 stress 0.06402627 
    ## Run 100 stress 0.06463313 
    ## Run 101 stress 0.06330869 
    ## ... Procrustes: rmse 3.719137e-05  max resid 0.0001003927 
    ## ... Similar to previous best
    ## Run 102 stress 0.06402624 
    ## Run 103 stress 0.06402637 
    ## Run 104 stress 0.06330868 
    ## ... Procrustes: rmse 2.140974e-05  max resid 4.956617e-05 
    ## ... Similar to previous best
    ## Run 105 stress 0.08301932 
    ## Run 106 stress 0.06402616 
    ## Run 107 stress 0.06402624 
    ## Run 108 stress 0.06402622 
    ## Run 109 stress 0.08315604 
    ## Run 110 stress 0.06330869 
    ## ... Procrustes: rmse 2.521612e-05  max resid 4.260481e-05 
    ## ... Similar to previous best
    ## Run 111 stress 0.06463313 
    ## Run 112 stress 0.06330872 
    ## ... Procrustes: rmse 7.965851e-05  max resid 0.0001739639 
    ## ... Similar to previous best
    ## Run 113 stress 0.08315607 
    ## Run 114 stress 0.0633087 
    ## ... Procrustes: rmse 5.814818e-05  max resid 0.0001283323 
    ## ... Similar to previous best
    ## Run 115 stress 0.06463313 
    ## Run 116 stress 0.06330868 
    ## ... Procrustes: rmse 2.342584e-05  max resid 8.364421e-05 
    ## ... Similar to previous best
    ## Run 117 stress 0.06463313 
    ## Run 118 stress 0.06463315 
    ## Run 119 stress 0.06330869 
    ## ... Procrustes: rmse 3.607099e-05  max resid 8.241706e-05 
    ## ... Similar to previous best
    ## Run 120 stress 0.06463313 
    ## Run 121 stress 0.0633087 
    ## ... Procrustes: rmse 1.659984e-05  max resid 3.041424e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.06330875 
    ## ... Procrustes: rmse 9.444468e-05  max resid 0.0002062373 
    ## ... Similar to previous best
    ## Run 123 stress 0.0633087 
    ## ... Procrustes: rmse 4.787595e-05  max resid 0.0001407116 
    ## ... Similar to previous best
    ## Run 124 stress 0.06330873 
    ## ... Procrustes: rmse 9.576767e-05  max resid 0.0002134215 
    ## ... Similar to previous best
    ## Run 125 stress 0.08301934 
    ## Run 126 stress 0.06402633 
    ## Run 127 stress 0.06402618 
    ## Run 128 stress 0.06330868 
    ## ... Procrustes: rmse 1.884581e-05  max resid 4.124198e-05 
    ## ... Similar to previous best
    ## Run 129 stress 0.06402634 
    ## Run 130 stress 0.06402622 
    ## Run 131 stress 0.0812621 
    ## Run 132 stress 0.06463313 
    ## Run 133 stress 0.06402622 
    ## Run 134 stress 0.06463313 
    ## Run 135 stress 0.06330868 
    ## ... Procrustes: rmse 1.101122e-05  max resid 2.897934e-05 
    ## ... Similar to previous best
    ## Run 136 stress 0.08301931 
    ## Run 137 stress 0.08126171 
    ## Run 138 stress 0.06330868 
    ## ... Procrustes: rmse 1.345695e-05  max resid 2.564019e-05 
    ## ... Similar to previous best
    ## Run 139 stress 0.0633087 
    ## ... Procrustes: rmse 1.339111e-05  max resid 3.692261e-05 
    ## ... Similar to previous best
    ## Run 140 stress 0.06330869 
    ## ... Procrustes: rmse 2.998841e-05  max resid 0.000103617 
    ## ... Similar to previous best
    ## Run 141 stress 0.0941671 
    ## Run 142 stress 0.09291 
    ## Run 143 stress 0.06330868 
    ## ... Procrustes: rmse 1.028086e-05  max resid 2.359946e-05 
    ## ... Similar to previous best
    ## Run 144 stress 0.06330868 
    ## ... Procrustes: rmse 1.313806e-05  max resid 2.770146e-05 
    ## ... Similar to previous best
    ## Run 145 stress 0.08344979 
    ## Run 146 stress 0.09290984 
    ## Run 147 stress 0.06402619 
    ## Run 148 stress 0.06330868 
    ## ... Procrustes: rmse 7.695648e-06  max resid 2.621462e-05 
    ## ... Similar to previous best
    ## Run 149 stress 0.08351229 
    ## Run 150 stress 0.06463313 
    ## Run 151 stress 0.06330871 
    ## ... Procrustes: rmse 7.489623e-05  max resid 0.0001630323 
    ## ... Similar to previous best
    ## Run 152 stress 0.06463316 
    ## Run 153 stress 0.06463314 
    ## Run 154 stress 0.06402617 
    ## Run 155 stress 0.06463313 
    ## Run 156 stress 0.08089061 
    ## Run 157 stress 0.06463313 
    ## Run 158 stress 0.06463314 
    ## Run 159 stress 0.06402617 
    ## Run 160 stress 0.06463314 
    ## Run 161 stress 0.06330873 
    ## ... Procrustes: rmse 9.391777e-05  max resid 0.0002016148 
    ## ... Similar to previous best
    ## Run 162 stress 0.06463313 
    ## Run 163 stress 0.06330871 
    ## ... Procrustes: rmse 6.45518e-05  max resid 0.0001438094 
    ## ... Similar to previous best
    ## Run 164 stress 0.0633087 
    ## ... Procrustes: rmse 4.067762e-05  max resid 8.582287e-05 
    ## ... Similar to previous best
    ## Run 165 stress 0.08315615 
    ## Run 166 stress 0.06463313 
    ## Run 167 stress 0.08126216 
    ## Run 168 stress 0.06463313 
    ## Run 169 stress 0.09429107 
    ## Run 170 stress 0.06463316 
    ## Run 171 stress 0.06330868 
    ## ... Procrustes: rmse 1.536244e-05  max resid 3.508365e-05 
    ## ... Similar to previous best
    ## Run 172 stress 0.08126166 
    ## Run 173 stress 0.06330868 
    ## ... Procrustes: rmse 5.162752e-06  max resid 1.19343e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.06463313 
    ## Run 175 stress 0.06463313 
    ## Run 176 stress 0.06463313 
    ## Run 177 stress 0.06330868 
    ## ... Procrustes: rmse 7.962905e-06  max resid 1.588363e-05 
    ## ... Similar to previous best
    ## Run 178 stress 0.06402617 
    ## Run 179 stress 0.0640263 
    ## Run 180 stress 0.06463314 
    ## Run 181 stress 0.06330869 
    ## ... Procrustes: rmse 4.399484e-05  max resid 9.88369e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.06330868 
    ## ... Procrustes: rmse 2.902723e-06  max resid 8.63012e-06 
    ## ... Similar to previous best
    ## Run 183 stress 0.06330871 
    ## ... Procrustes: rmse 6.673309e-05  max resid 0.0001499187 
    ## ... Similar to previous best
    ## Run 184 stress 0.09283196 
    ## Run 185 stress 0.06330878 
    ## ... Procrustes: rmse 0.0001285929  max resid 0.0002878822 
    ## ... Similar to previous best
    ## Run 186 stress 0.06463313 
    ## Run 187 stress 0.08089054 
    ## Run 188 stress 0.0808908 
    ## Run 189 stress 0.06402617 
    ## Run 190 stress 0.06463316 
    ## Run 191 stress 0.06463315 
    ## Run 192 stress 0.06463314 
    ## Run 193 stress 0.06463316 
    ## Run 194 stress 0.0640262 
    ## Run 195 stress 0.0633087 
    ## ... Procrustes: rmse 5.338143e-05  max resid 0.0001196209 
    ## ... Similar to previous best
    ## Run 196 stress 0.06463314 
    ## Run 197 stress 0.06463313 
    ## Run 198 stress 0.06402617 
    ## Run 199 stress 0.06402616 
    ## Run 200 stress 0.08403559 
    ## Run 201 stress 0.06330868 
    ## ... Procrustes: rmse 2.386926e-05  max resid 8.674162e-05 
    ## ... Similar to previous best
    ## Run 202 stress 0.09246001 
    ## Run 203 stress 0.06330868 
    ## ... Procrustes: rmse 1.803104e-05  max resid 6.002452e-05 
    ## ... Similar to previous best
    ## Run 204 stress 0.06463313 
    ## Run 205 stress 0.06402622 
    ## Run 206 stress 0.06463313 
    ## Run 207 stress 0.06402621 
    ## Run 208 stress 0.06463316 
    ## Run 209 stress 0.06463315 
    ## Run 210 stress 0.06402618 
    ## Run 211 stress 0.0633087 
    ## ... Procrustes: rmse 5.613386e-05  max resid 0.0001534578 
    ## ... Similar to previous best
    ## Run 212 stress 0.0633087 
    ## ... Procrustes: rmse 5.516839e-05  max resid 0.0001215317 
    ## ... Similar to previous best
    ## Run 213 stress 0.06463314 
    ## Run 214 stress 0.06463315 
    ## Run 215 stress 0.06463313 
    ## Run 216 stress 0.06463313 
    ## Run 217 stress 0.08233751 
    ## Run 218 stress 0.06402617 
    ## Run 219 stress 0.06330871 
    ## ... Procrustes: rmse 5.537497e-05  max resid 9.681577e-05 
    ## ... Similar to previous best
    ## Run 220 stress 0.06330869 
    ## ... Procrustes: rmse 3.018166e-05  max resid 6.534512e-05 
    ## ... Similar to previous best
    ## Run 221 stress 0.06463314 
    ## Run 222 stress 0.06330869 
    ## ... Procrustes: rmse 3.956392e-05  max resid 8.771519e-05 
    ## ... Similar to previous best
    ## Run 223 stress 0.06402621 
    ## Run 224 stress 0.09106917 
    ## Run 225 stress 0.06402631 
    ## Run 226 stress 0.06330875 
    ## ... Procrustes: rmse 0.000107692  max resid 0.0002418327 
    ## ... Similar to previous best
    ## Run 227 stress 0.06402632 
    ## Run 228 stress 0.0633087 
    ## ... Procrustes: rmse 5.265535e-05  max resid 0.0001168544 
    ## ... Similar to previous best
    ## Run 229 stress 0.06463315 
    ## Run 230 stress 0.06463313 
    ## Run 231 stress 0.06402619 
    ## Run 232 stress 0.08126183 
    ## Run 233 stress 0.06463314 
    ## Run 234 stress 0.08233734 
    ## Run 235 stress 0.06463315 
    ## Run 236 stress 0.06402624 
    ## Run 237 stress 0.06402616 
    ## Run 238 stress 0.06330868 
    ## ... Procrustes: rmse 1.228425e-05  max resid 2.539091e-05 
    ## ... Similar to previous best
    ## Run 239 stress 0.3698521 
    ## Run 240 stress 0.06330868 
    ## ... Procrustes: rmse 1.626752e-05  max resid 3.668596e-05 
    ## ... Similar to previous best
    ## Run 241 stress 0.06330868 
    ## ... Procrustes: rmse 2.865341e-06  max resid 7.336598e-06 
    ## ... Similar to previous best
    ## Run 242 stress 0.06402617 
    ## Run 243 stress 0.06402641 
    ## Run 244 stress 0.06463313 
    ## Run 245 stress 0.06330868 
    ## ... Procrustes: rmse 3.447055e-06  max resid 7.514447e-06 
    ## ... Similar to previous best
    ## Run 246 stress 0.06402621 
    ## Run 247 stress 0.06330868 
    ## ... Procrustes: rmse 1.000811e-05  max resid 2.273244e-05 
    ## ... Similar to previous best
    ## Run 248 stress 0.0977349 
    ## Run 249 stress 0.09233436 
    ## Run 250 stress 0.08379068 
    ## Run 251 stress 0.06402619 
    ## Run 252 stress 0.09472813 
    ## Run 253 stress 0.0633087 
    ## ... Procrustes: rmse 5.092634e-05  max resid 0.0001145404 
    ## ... Similar to previous best
    ## Run 254 stress 0.09106977 
    ## Run 255 stress 0.06402617 
    ## Run 256 stress 0.06330873 
    ## ... Procrustes: rmse 8.834751e-05  max resid 0.0001988402 
    ## ... Similar to previous best
    ## Run 257 stress 0.0633087 
    ## ... Procrustes: rmse 5.943851e-05  max resid 0.00013066 
    ## ... Similar to previous best
    ## Run 258 stress 0.08403568 
    ## Run 259 stress 0.06402616 
    ## Run 260 stress 0.06402623 
    ## Run 261 stress 0.08315607 
    ## Run 262 stress 0.06330868 
    ## ... Procrustes: rmse 1.874406e-05  max resid 4.108115e-05 
    ## ... Similar to previous best
    ## Run 263 stress 0.06463314 
    ## Run 264 stress 0.06463314 
    ## Run 265 stress 0.06463315 
    ## Run 266 stress 0.06402617 
    ## Run 267 stress 0.06330871 
    ## ... Procrustes: rmse 5.980398e-05  max resid 0.0001344196 
    ## ... Similar to previous best
    ## Run 268 stress 0.0921719 
    ## Run 269 stress 0.06402627 
    ## Run 270 stress 0.06402616 
    ## Run 271 stress 0.06463313 
    ## Run 272 stress 0.06463313 
    ## Run 273 stress 0.06402627 
    ## Run 274 stress 0.06330869 
    ## ... Procrustes: rmse 1.512902e-05  max resid 4.108921e-05 
    ## ... Similar to previous best
    ## Run 275 stress 0.06330868 
    ## ... Procrustes: rmse 5.264892e-06  max resid 9.667559e-06 
    ## ... Similar to previous best
    ## Run 276 stress 0.0941668 
    ## Run 277 stress 0.0640262 
    ## Run 278 stress 0.08089075 
    ## Run 279 stress 0.06402635 
    ## Run 280 stress 0.06463313 
    ## Run 281 stress 0.06463313 
    ## Run 282 stress 0.06330869 
    ## ... Procrustes: rmse 3.614545e-05  max resid 8.168084e-05 
    ## ... Similar to previous best
    ## Run 283 stress 0.06463313 
    ## Run 284 stress 0.06463314 
    ## Run 285 stress 0.06463313 
    ## Run 286 stress 0.0640263 
    ## Run 287 stress 0.06402616 
    ## Run 288 stress 0.06463313 
    ## Run 289 stress 0.06402622 
    ## Run 290 stress 0.06330868 
    ## ... Procrustes: rmse 1.909708e-05  max resid 5.205621e-05 
    ## ... Similar to previous best
    ## Run 291 stress 0.08126177 
    ## Run 292 stress 0.06463313 
    ## Run 293 stress 0.06330868 
    ## ... Procrustes: rmse 2.696683e-06  max resid 6.400435e-06 
    ## ... Similar to previous best
    ## Run 294 stress 0.08089017 
    ## Run 295 stress 0.0633087 
    ## ... Procrustes: rmse 5.950403e-05  max resid 0.0001294176 
    ## ... Similar to previous best
    ## Run 296 stress 0.06402629 
    ## Run 297 stress 0.06463313 
    ## Run 298 stress 0.06463313 
    ## Run 299 stress 0.06463314 
    ## Run 300 stress 0.08357057 
    ## Run 301 stress 0.09472685 
    ## Run 302 stress 0.08431267 
    ## Run 303 stress 0.09416691 
    ## Run 304 stress 0.06402636 
    ## Run 305 stress 0.06463313 
    ## Run 306 stress 0.06330868 
    ## ... Procrustes: rmse 1.400651e-05  max resid 4.866648e-05 
    ## ... Similar to previous best
    ## Run 307 stress 0.06402626 
    ## Run 308 stress 0.06402639 
    ## Run 309 stress 0.09547957 
    ## Run 310 stress 0.06463313 
    ## Run 311 stress 0.08301935 
    ## Run 312 stress 0.09232066 
    ## Run 313 stress 0.08327919 
    ## Run 314 stress 0.06463313 
    ## Run 315 stress 0.06330868 
    ## ... Procrustes: rmse 9.529071e-06  max resid 3.170523e-05 
    ## ... Similar to previous best
    ## Run 316 stress 0.06463318 
    ## Run 317 stress 0.08351253 
    ## Run 318 stress 0.3699023 
    ## Run 319 stress 0.06463313 
    ## Run 320 stress 0.06402619 
    ## Run 321 stress 0.06402618 
    ## Run 322 stress 0.06330868 
    ## ... Procrustes: rmse 6.362271e-06  max resid 1.67512e-05 
    ## ... Similar to previous best
    ## Run 323 stress 0.09533683 
    ## Run 324 stress 0.09712295 
    ## Run 325 stress 0.06463315 
    ## Run 326 stress 0.09563707 
    ## Run 327 stress 0.06463313 
    ## Run 328 stress 0.06463315 
    ## Run 329 stress 0.06402616 
    ## Run 330 stress 0.06463313 
    ## Run 331 stress 0.06330869 
    ## ... Procrustes: rmse 4.783403e-05  max resid 0.0001063631 
    ## ... Similar to previous best
    ## Run 332 stress 0.06463314 
    ## Run 333 stress 0.06330868 
    ## ... Procrustes: rmse 2.250232e-05  max resid 5.072445e-05 
    ## ... Similar to previous best
    ## Run 334 stress 0.06463313 
    ## Run 335 stress 0.06402624 
    ## Run 336 stress 0.06463314 
    ## Run 337 stress 0.06402616 
    ## Run 338 stress 0.06463313 
    ## Run 339 stress 0.0633088 
    ## ... Procrustes: rmse 0.0001389867  max resid 0.0003102229 
    ## ... Similar to previous best
    ## Run 340 stress 0.06402619 
    ## Run 341 stress 0.06402639 
    ## Run 342 stress 0.06402619 
    ## Run 343 stress 0.06463313 
    ## Run 344 stress 0.06402618 
    ## Run 345 stress 0.06330869 
    ## ... Procrustes: rmse 4.174649e-05  max resid 9.236822e-05 
    ## ... Similar to previous best
    ## Run 346 stress 0.09504606 
    ## Run 347 stress 0.09429109 
    ## Run 348 stress 0.06402637 
    ## Run 349 stress 0.09233441 
    ## Run 350 stress 0.09106925 
    ## Run 351 stress 0.3622739 
    ## Run 352 stress 0.06330868 
    ## ... Procrustes: rmse 2.357582e-05  max resid 7.668318e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.06402629 
    ## Run 354 stress 0.09429106 
    ## Run 355 stress 0.06330868 
    ## ... Procrustes: rmse 8.526831e-06  max resid 1.672902e-05 
    ## ... Similar to previous best
    ## Run 356 stress 0.08301936 
    ## Run 357 stress 0.08089028 
    ## Run 358 stress 0.06463314 
    ## Run 359 stress 0.0633087 
    ## ... Procrustes: rmse 4.191306e-05  max resid 9.315013e-05 
    ## ... Similar to previous best
    ## Run 360 stress 0.06330869 
    ## ... Procrustes: rmse 3.100738e-05  max resid 9.714912e-05 
    ## ... Similar to previous best
    ## Run 361 stress 0.06463314 
    ## Run 362 stress 0.06402624 
    ## Run 363 stress 0.06463313 
    ## Run 364 stress 0.06330869 
    ## ... Procrustes: rmse 3.069077e-05  max resid 0.0001152266 
    ## ... Similar to previous best
    ## Run 365 stress 0.06463317 
    ## Run 366 stress 0.06330868 
    ## ... Procrustes: rmse 2.467589e-05  max resid 5.624826e-05 
    ## ... Similar to previous best
    ## Run 367 stress 0.06463314 
    ## Run 368 stress 0.08315619 
    ## Run 369 stress 0.06402618 
    ## Run 370 stress 0.09206508 
    ## Run 371 stress 0.06330868 
    ## ... Procrustes: rmse 1.247528e-05  max resid 2.756072e-05 
    ## ... Similar to previous best
    ## Run 372 stress 0.06463313 
    ## Run 373 stress 0.09581765 
    ## Run 374 stress 0.06463314 
    ## Run 375 stress 0.0923344 
    ## Run 376 stress 0.06402617 
    ## Run 377 stress 0.06402627 
    ## Run 378 stress 0.06330869 
    ## ... Procrustes: rmse 3.285395e-05  max resid 6.460283e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.06463313 
    ## Run 380 stress 0.0959145 
    ## Run 381 stress 0.0633087 
    ## ... Procrustes: rmse 4.783597e-05  max resid 8.513117e-05 
    ## ... Similar to previous best
    ## Run 382 stress 0.06330868 
    ## ... Procrustes: rmse 1.424937e-05  max resid 4.75553e-05 
    ## ... Similar to previous best
    ## Run 383 stress 0.06463314 
    ## Run 384 stress 0.06330871 
    ## ... Procrustes: rmse 6.15766e-05  max resid 0.0001346233 
    ## ... Similar to previous best
    ## Run 385 stress 0.06402626 
    ## Run 386 stress 0.06463313 
    ## Run 387 stress 0.09290995 
    ## Run 388 stress 0.06463317 
    ## Run 389 stress 0.0633087 
    ## ... Procrustes: rmse 5.783845e-05  max resid 0.0001269377 
    ## ... Similar to previous best
    ## Run 390 stress 0.0633087 
    ## ... Procrustes: rmse 6.173603e-05  max resid 0.0001369903 
    ## ... Similar to previous best
    ## Run 391 stress 0.06330868 
    ## ... Procrustes: rmse 3.215011e-06  max resid 7.136895e-06 
    ## ... Similar to previous best
    ## Run 392 stress 0.0838997 
    ## Run 393 stress 0.06402617 
    ## Run 394 stress 0.06330868 
    ## ... Procrustes: rmse 9.705363e-06  max resid 3.42887e-05 
    ## ... Similar to previous best
    ## Run 395 stress 0.06402623 
    ## Run 396 stress 0.06402619 
    ## Run 397 stress 0.06330869 
    ## ... Procrustes: rmse 2.820955e-05  max resid 0.0001053053 
    ## ... Similar to previous best
    ## Run 398 stress 0.092701 
    ## Run 399 stress 0.06330871 
    ## ... Procrustes: rmse 6.458182e-05  max resid 0.0001440628 
    ## ... Similar to previous best
    ## Run 400 stress 0.08089071 
    ## Run 401 stress 0.09463929 
    ## Run 402 stress 0.06402626 
    ## Run 403 stress 0.06330868 
    ## ... Procrustes: rmse 3.313681e-06  max resid 5.630514e-06 
    ## ... Similar to previous best
    ## Run 404 stress 0.06402637 
    ## Run 405 stress 0.08233736 
    ## Run 406 stress 0.06463313 
    ## Run 407 stress 0.06330872 
    ## ... Procrustes: rmse 7.567923e-05  max resid 0.0001617963 
    ## ... Similar to previous best
    ## Run 408 stress 0.0640262 
    ## Run 409 stress 0.08089052 
    ## Run 410 stress 0.06463313 
    ## Run 411 stress 0.06402632 
    ## Run 412 stress 0.06463313 
    ## Run 413 stress 0.09625615 
    ## Run 414 stress 0.0834497 
    ## Run 415 stress 0.06330868 
    ## ... Procrustes: rmse 1.565032e-05  max resid 3.426612e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.06463314 
    ## Run 417 stress 0.06330871 
    ## ... Procrustes: rmse 6.932625e-05  max resid 0.0001528834 
    ## ... Similar to previous best
    ## Run 418 stress 0.06402633 
    ## Run 419 stress 0.0640262 
    ## Run 420 stress 0.06402628 
    ## Run 421 stress 0.06402619 
    ## Run 422 stress 0.06330871 
    ## ... Procrustes: rmse 5.528329e-05  max resid 0.0001094951 
    ## ... Similar to previous best
    ## Run 423 stress 0.08126186 
    ## Run 424 stress 0.06463313 
    ## Run 425 stress 0.08233735 
    ## Run 426 stress 0.06402619 
    ## Run 427 stress 0.06463314 
    ## Run 428 stress 0.06463313 
    ## Run 429 stress 0.06463314 
    ## Run 430 stress 0.09245996 
    ## Run 431 stress 0.06402623 
    ## Run 432 stress 0.06463313 
    ## Run 433 stress 0.06402621 
    ## Run 434 stress 0.06330869 
    ## ... Procrustes: rmse 4.015172e-05  max resid 9.093983e-05 
    ## ... Similar to previous best
    ## Run 435 stress 0.0823373 
    ## Run 436 stress 0.06463313 
    ## Run 437 stress 0.09504628 
    ## Run 438 stress 0.09245933 
    ## Run 439 stress 0.06330876 
    ## ... Procrustes: rmse 0.0001125537  max resid 0.0002520886 
    ## ... Similar to previous best
    ## Run 440 stress 0.09563686 
    ## Run 441 stress 0.0640262 
    ## Run 442 stress 0.06402619 
    ## Run 443 stress 0.06330868 
    ## ... Procrustes: rmse 6.875274e-06  max resid 1.240942e-05 
    ## ... Similar to previous best
    ## Run 444 stress 0.06463314 
    ## Run 445 stress 0.06330868 
    ## ... Procrustes: rmse 1.902016e-05  max resid 5.29045e-05 
    ## ... Similar to previous best
    ## Run 446 stress 0.06330868 
    ## ... Procrustes: rmse 2.675568e-05  max resid 5.94575e-05 
    ## ... Similar to previous best
    ## Run 447 stress 0.06402628 
    ## Run 448 stress 0.06330868 
    ## ... Procrustes: rmse 2.353752e-05  max resid 7.296287e-05 
    ## ... Similar to previous best
    ## Run 449 stress 0.06402628 
    ## Run 450 stress 0.06402617 
    ## Run 451 stress 0.06330868 
    ## ... Procrustes: rmse 3.276262e-06  max resid 1.099496e-05 
    ## ... Similar to previous best
    ## Run 452 stress 0.09876085 
    ## Run 453 stress 0.06463313 
    ## Run 454 stress 0.0840356 
    ## Run 455 stress 0.0640262 
    ## Run 456 stress 0.06402639 
    ## Run 457 stress 0.06402617 
    ## Run 458 stress 0.08126187 
    ## Run 459 stress 0.06402629 
    ## Run 460 stress 0.06402616 
    ## Run 461 stress 0.06402616 
    ## Run 462 stress 0.08089045 
    ## Run 463 stress 0.09233437 
    ## Run 464 stress 0.06330868 
    ## ... Procrustes: rmse 2.0585e-05  max resid 6.023227e-05 
    ## ... Similar to previous best
    ## Run 465 stress 0.08126193 
    ## Run 466 stress 0.06330869 
    ## ... Procrustes: rmse 2.546783e-05  max resid 4.485936e-05 
    ## ... Similar to previous best
    ## Run 467 stress 0.06402616 
    ## Run 468 stress 0.06402616 
    ## Run 469 stress 0.06463313 
    ## Run 470 stress 0.06402622 
    ## Run 471 stress 0.0946204 
    ## Run 472 stress 0.06330868 
    ## ... Procrustes: rmse 5.961835e-06  max resid 1.073175e-05 
    ## ... Similar to previous best
    ## Run 473 stress 0.08357049 
    ## Run 474 stress 0.06463314 
    ## Run 475 stress 0.06402636 
    ## Run 476 stress 0.08357052 
    ## Run 477 stress 0.09291006 
    ## Run 478 stress 0.06402642 
    ## Run 479 stress 0.06402623 
    ## Run 480 stress 0.06402636 
    ## Run 481 stress 0.06463315 
    ## Run 482 stress 0.06402621 
    ## Run 483 stress 0.06463313 
    ## Run 484 stress 0.09485781 
    ## Run 485 stress 0.06402616 
    ## Run 486 stress 0.09245981 
    ## Run 487 stress 0.06402634 
    ## Run 488 stress 0.06463314 
    ## Run 489 stress 0.06402628 
    ## Run 490 stress 0.06330869 
    ## ... Procrustes: rmse 4.399442e-05  max resid 0.0001004614 
    ## ... Similar to previous best
    ## Run 491 stress 0.08246345 
    ## Run 492 stress 0.0640262 
    ## Run 493 stress 0.08089074 
    ## Run 494 stress 0.06463316 
    ## Run 495 stress 0.06463315 
    ## Run 496 stress 0.06402616 
    ## Run 497 stress 0.08327924 
    ## Run 498 stress 0.06463313 
    ## Run 499 stress 0.0923345 
    ## Run 500 stress 0.06402621 
    ## *** Best solution repeated 102 times

``` r
# Mixed and stratified lakes
PD_beta_env_MS_NMDS <- metaMDS(PD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.05051299 
    ## ... Procrustes: rmse 0.03575975  max resid 0.1220672 
    ## Run 2 stress 0.06231304 
    ## Run 3 stress 0.06478523 
    ## Run 4 stress 0.05466711 
    ## Run 5 stress 0.05466706 
    ## Run 6 stress 0.05171058 
    ## Run 7 stress 0.05036793 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002254124  max resid 0.0005714592 
    ## ... Similar to previous best
    ## Run 8 stress 0.05979226 
    ## Run 9 stress 0.05171059 
    ## Run 10 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 3.68125e-05  max resid 8.683768e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.05171059 
    ## Run 12 stress 0.05171058 
    ## Run 13 stress 0.06904471 
    ## Run 14 stress 0.0505134 
    ## ... Procrustes: rmse 0.03569734  max resid 0.1217911 
    ## Run 15 stress 0.05036797 
    ## ... Procrustes: rmse 8.455938e-05  max resid 0.0001821458 
    ## ... Similar to previous best
    ## Run 16 stress 0.06943856 
    ## Run 17 stress 0.05036791 
    ## ... Procrustes: rmse 3.766421e-06  max resid 7.823995e-06 
    ## ... Similar to previous best
    ## Run 18 stress 0.05036798 
    ## ... Procrustes: rmse 9.576877e-05  max resid 0.0002334579 
    ## ... Similar to previous best
    ## Run 19 stress 0.05051312 
    ## ... Procrustes: rmse 0.03574086  max resid 0.1219548 
    ## Run 20 stress 0.05171062 
    ## Run 21 stress 0.05466703 
    ## Run 22 stress 0.05403579 
    ## Run 23 stress 0.06493686 
    ## Run 24 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001876382  max resid 0.0004905277 
    ## ... Similar to previous best
    ## Run 25 stress 0.05051303 
    ## ... Procrustes: rmse 0.03562292  max resid 0.1217285 
    ## Run 26 stress 0.05051342 
    ## ... Procrustes: rmse 0.03571323  max resid 0.1220191 
    ## Run 27 stress 0.05036796 
    ## ... Procrustes: rmse 6.560699e-05  max resid 0.00016143 
    ## ... Similar to previous best
    ## Run 28 stress 0.05171058 
    ## Run 29 stress 0.050368 
    ## ... Procrustes: rmse 0.0001745194  max resid 0.0004539177 
    ## ... Similar to previous best
    ## Run 30 stress 0.05466702 
    ## Run 31 stress 0.06630524 
    ## Run 32 stress 0.05928215 
    ## Run 33 stress 0.06409274 
    ## Run 34 stress 0.05466702 
    ## Run 35 stress 0.06640366 
    ## Run 36 stress 0.06231309 
    ## Run 37 stress 0.050368 
    ## ... Procrustes: rmse 0.0001174736  max resid 0.0002895139 
    ## ... Similar to previous best
    ## Run 38 stress 0.05036796 
    ## ... Procrustes: rmse 8.398573e-05  max resid 0.0001650469 
    ## ... Similar to previous best
    ## Run 39 stress 0.05036807 
    ## ... Procrustes: rmse 0.0001605374  max resid 0.0003953327 
    ## ... Similar to previous best
    ## Run 40 stress 0.06309851 
    ## Run 41 stress 0.05051295 
    ## ... Procrustes: rmse 0.0356755  max resid 0.1218828 
    ## Run 42 stress 0.05171061 
    ## Run 43 stress 0.05602043 
    ## Run 44 stress 0.05051316 
    ## ... Procrustes: rmse 0.03569993  max resid 0.1219688 
    ## Run 45 stress 0.05466702 
    ## Run 46 stress 0.05036795 
    ## ... Procrustes: rmse 6.433181e-05  max resid 0.0001455599 
    ## ... Similar to previous best
    ## Run 47 stress 0.05051278 
    ## ... Procrustes: rmse 0.03567998  max resid 0.1218189 
    ## Run 48 stress 0.06136682 
    ## Run 49 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001228487  max resid 0.0003145482 
    ## ... Similar to previous best
    ## Run 50 stress 0.06499242 
    ## Run 51 stress 0.05051318 
    ## ... Procrustes: rmse 0.03568894  max resid 0.12179 
    ## Run 52 stress 0.05171059 
    ## Run 53 stress 0.05968626 
    ## Run 54 stress 0.05466706 
    ## Run 55 stress 0.05171059 
    ## Run 56 stress 0.05337582 
    ## Run 57 stress 0.05051286 
    ## ... Procrustes: rmse 0.03567398  max resid 0.1217869 
    ## Run 58 stress 0.05036806 
    ## ... Procrustes: rmse 0.0001597661  max resid 0.0004043104 
    ## ... Similar to previous best
    ## Run 59 stress 0.06730247 
    ## Run 60 stress 0.05403625 
    ## Run 61 stress 0.05928238 
    ## Run 62 stress 0.06309843 
    ## Run 63 stress 0.05051338 
    ## ... Procrustes: rmse 0.03566682  max resid 0.1217004 
    ## Run 64 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001217588  max resid 0.0003174217 
    ## ... Similar to previous best
    ## Run 65 stress 0.05466705 
    ## Run 66 stress 0.05602079 
    ## Run 67 stress 0.05036792 
    ## ... Procrustes: rmse 8.057569e-05  max resid 0.0001983093 
    ## ... Similar to previous best
    ## Run 68 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001209218  max resid 0.0003383524 
    ## ... Similar to previous best
    ## Run 69 stress 0.0517106 
    ## Run 70 stress 0.05829948 
    ## Run 71 stress 0.06826695 
    ## Run 72 stress 0.06478534 
    ## Run 73 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001017346  max resid 0.0002526461 
    ## ... Similar to previous best
    ## Run 74 stress 0.05036799 
    ## ... Procrustes: rmse 9.979563e-05  max resid 0.0002471946 
    ## ... Similar to previous best
    ## Run 75 stress 0.05051314 
    ## ... Procrustes: rmse 0.03569509  max resid 0.1219534 
    ## Run 76 stress 0.06203043 
    ## Run 77 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001223213  max resid 0.0003084835 
    ## ... Similar to previous best
    ## Run 78 stress 0.05036797 
    ## ... Procrustes: rmse 8.785448e-05  max resid 0.0002450734 
    ## ... Similar to previous best
    ## Run 79 stress 0.05171058 
    ## Run 80 stress 0.05036797 
    ## ... Procrustes: rmse 0.000151515  max resid 0.0003834592 
    ## ... Similar to previous best
    ## Run 81 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001380182  max resid 0.0003491475 
    ## ... Similar to previous best
    ## Run 82 stress 0.05051291 
    ## ... Procrustes: rmse 0.03558899  max resid 0.1216159 
    ## Run 83 stress 0.05928211 
    ## Run 84 stress 0.05051326 
    ## ... Procrustes: rmse 0.03568038  max resid 0.1217541 
    ## Run 85 stress 0.05051291 
    ## ... Procrustes: rmse 0.03567387  max resid 0.1218748 
    ## Run 86 stress 0.05602058 
    ## Run 87 stress 0.050368 
    ## ... Procrustes: rmse 0.0001116792  max resid 0.0002688697 
    ## ... Similar to previous best
    ## Run 88 stress 0.07026205 
    ## Run 89 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001255034  max resid 0.0002771332 
    ## ... Similar to previous best
    ## Run 90 stress 0.05036792 
    ## ... Procrustes: rmse 2.600756e-05  max resid 6.442505e-05 
    ## ... Similar to previous best
    ## Run 91 stress 0.05036806 
    ## ... Procrustes: rmse 0.0001540608  max resid 0.0003852677 
    ## ... Similar to previous best
    ## Run 92 stress 0.05602049 
    ## Run 93 stress 0.05171061 
    ## Run 94 stress 0.06730222 
    ## Run 95 stress 0.05051288 
    ## ... Procrustes: rmse 0.03580914  max resid 0.1222392 
    ## Run 96 stress 0.05036792 
    ## ... Procrustes: rmse 1.480553e-05  max resid 3.358194e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.05051342 
    ## ... Procrustes: rmse 0.03565619  max resid 0.1216731 
    ## Run 98 stress 0.05036795 
    ## ... Procrustes: rmse 7.115305e-05  max resid 0.000162505 
    ## ... Similar to previous best
    ## Run 99 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001253132  max resid 0.0003392603 
    ## ... Similar to previous best
    ## Run 100 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001276215  max resid 0.0002894829 
    ## ... Similar to previous best
    ## Run 101 stress 0.05928228 
    ## Run 102 stress 0.06478532 
    ## Run 103 stress 0.05466702 
    ## Run 104 stress 0.06111764 
    ## Run 105 stress 0.05036794 
    ## ... Procrustes: rmse 5.478253e-05  max resid 0.0001351756 
    ## ... Similar to previous best
    ## Run 106 stress 0.05466708 
    ## Run 107 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001222841  max resid 0.0002993942 
    ## ... Similar to previous best
    ## Run 108 stress 0.05403582 
    ## Run 109 stress 0.05171061 
    ## Run 110 stress 0.05051315 
    ## ... Procrustes: rmse 0.03567845  max resid 0.1217594 
    ## Run 111 stress 0.06537607 
    ## Run 112 stress 0.05051291 
    ## ... Procrustes: rmse 0.03557339  max resid 0.1215613 
    ## Run 113 stress 0.05171061 
    ## Run 114 stress 0.0724227 
    ## Run 115 stress 0.06309863 
    ## Run 116 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001690828  max resid 0.0004329663 
    ## ... Similar to previous best
    ## Run 117 stress 0.06609053 
    ## Run 118 stress 0.05051335 
    ## ... Procrustes: rmse 0.0357169  max resid 0.1220278 
    ## Run 119 stress 0.06478527 
    ## Run 120 stress 0.0592823 
    ## Run 121 stress 0.05337594 
    ## Run 122 stress 0.05051314 
    ## ... Procrustes: rmse 0.03569529  max resid 0.1219539 
    ## Run 123 stress 0.0505134 
    ## ... Procrustes: rmse 0.03569908  max resid 0.1217984 
    ## Run 124 stress 0.05928212 
    ## Run 125 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001793295  max resid 0.0004830736 
    ## ... Similar to previous best
    ## Run 126 stress 0.06309865 
    ## Run 127 stress 0.05466712 
    ## Run 128 stress 0.05171058 
    ## Run 129 stress 0.05051349 
    ## ... Procrustes: rmse 0.03572173  max resid 0.1220424 
    ## Run 130 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001544264  max resid 0.0003885693 
    ## ... Similar to previous best
    ## Run 131 stress 0.05829938 
    ## Run 132 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001882508  max resid 0.0005057213 
    ## ... Similar to previous best
    ## Run 133 stress 0.05051313 
    ## ... Procrustes: rmse 0.03567881  max resid 0.1219039 
    ## Run 134 stress 0.054667 
    ## Run 135 stress 0.05171058 
    ## Run 136 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001480248  max resid 0.0003670419 
    ## ... Similar to previous best
    ## Run 137 stress 0.05036797 
    ## ... Procrustes: rmse 8.796051e-05  max resid 0.0002257243 
    ## ... Similar to previous best
    ## Run 138 stress 0.06309864 
    ## Run 139 stress 0.05928225 
    ## Run 140 stress 0.05036796 
    ## ... Procrustes: rmse 8.015836e-05  max resid 0.0001984972 
    ## ... Similar to previous best
    ## Run 141 stress 0.05171059 
    ## Run 142 stress 0.05036807 
    ## ... Procrustes: rmse 8.669617e-05  max resid 0.0001934988 
    ## ... Similar to previous best
    ## Run 143 stress 0.05171061 
    ## Run 144 stress 0.05403584 
    ## Run 145 stress 0.05403573 
    ## Run 146 stress 0.05171058 
    ## Run 147 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.037102e-05  max resid 2.461963e-05 
    ## ... Similar to previous best
    ## Run 148 stress 0.054667 
    ## Run 149 stress 0.06832188 
    ## Run 150 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001427652  max resid 0.0003804361 
    ## ... Similar to previous best
    ## Run 151 stress 0.05928214 
    ## Run 152 stress 0.05036796 
    ## ... Procrustes: rmse 8.458206e-05  max resid 0.0002010171 
    ## ... Similar to previous best
    ## Run 153 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001914497  max resid 0.0004708906 
    ## ... Similar to previous best
    ## Run 154 stress 0.05036791 
    ## ... Procrustes: rmse 1.813437e-05  max resid 3.831018e-05 
    ## ... Similar to previous best
    ## Run 155 stress 0.05979227 
    ## Run 156 stress 0.05337583 
    ## Run 157 stress 0.06885064 
    ## Run 158 stress 0.05829936 
    ## Run 159 stress 0.05051288 
    ## ... Procrustes: rmse 0.03554312  max resid 0.1214285 
    ## Run 160 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001658391  max resid 0.0004730277 
    ## ... Similar to previous best
    ## Run 161 stress 0.06231309 
    ## Run 162 stress 0.05036794 
    ## ... Procrustes: rmse 6.427001e-05  max resid 0.0001563032 
    ## ... Similar to previous best
    ## Run 163 stress 0.05171063 
    ## Run 164 stress 0.05829937 
    ## Run 165 stress 0.06730253 
    ## Run 166 stress 0.05968642 
    ## Run 167 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001502436  max resid 0.0003690603 
    ## ... Similar to previous best
    ## Run 168 stress 0.0505134 
    ## ... Procrustes: rmse 0.03567406  max resid 0.1217228 
    ## Run 169 stress 0.06027432 
    ## Run 170 stress 0.05829932 
    ## Run 171 stress 0.05036796 
    ## ... Procrustes: rmse 8.593065e-05  max resid 0.0002051784 
    ## ... Similar to previous best
    ## Run 172 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001491708  max resid 0.0003707373 
    ## ... Similar to previous best
    ## Run 173 stress 0.05036794 
    ## ... Procrustes: rmse 7.932546e-05  max resid 0.0001596462 
    ## ... Similar to previous best
    ## Run 174 stress 0.06478521 
    ## Run 175 stress 0.06640346 
    ## Run 176 stress 0.05051338 
    ## ... Procrustes: rmse 0.03571414  max resid 0.1220193 
    ## Run 177 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001375456  max resid 0.0003142522 
    ## ... Similar to previous best
    ## Run 178 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001731735  max resid 0.0004771933 
    ## ... Similar to previous best
    ## Run 179 stress 0.05403603 
    ## Run 180 stress 0.05051336 
    ## ... Procrustes: rmse 0.03572186  max resid 0.121872 
    ## Run 181 stress 0.05036792 
    ## ... Procrustes: rmse 1.963864e-05  max resid 4.604908e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.05928211 
    ## Run 183 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001779926  max resid 0.0004590214 
    ## ... Similar to previous best
    ## Run 184 stress 0.0723983 
    ## Run 185 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001238172  max resid 0.0003053177 
    ## ... Similar to previous best
    ## Run 186 stress 0.05171062 
    ## Run 187 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001833861  max resid 0.0005325188 
    ## ... Similar to previous best
    ## Run 188 stress 0.07242188 
    ## Run 189 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 3.306468e-05  max resid 7.014402e-05 
    ## ... Similar to previous best
    ## Run 190 stress 0.05051354 
    ## ... Procrustes: rmse 0.03570989  max resid 0.1220196 
    ## Run 191 stress 0.05051291 
    ## ... Procrustes: rmse 0.03569295  max resid 0.1218402 
    ## Run 192 stress 0.05337585 
    ## Run 193 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001424673  max resid 0.000341677 
    ## ... Similar to previous best
    ## Run 194 stress 0.05036793 
    ## ... Procrustes: rmse 7.500296e-05  max resid 0.0002150307 
    ## ... Similar to previous best
    ## Run 195 stress 0.0505133 
    ## ... Procrustes: rmse 0.03572148  max resid 0.1220451 
    ## Run 196 stress 0.05466702 
    ## Run 197 stress 0.05968632 
    ## Run 198 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001394088  max resid 0.0003823805 
    ## ... Similar to previous best
    ## Run 199 stress 0.06231286 
    ## Run 200 stress 0.05051318 
    ## ... Procrustes: rmse 0.03574882  max resid 0.1219787 
    ## Run 201 stress 0.05036793 
    ## ... Procrustes: rmse 5.232043e-05  max resid 0.0001394452 
    ## ... Similar to previous best
    ## Run 202 stress 0.05466701 
    ## Run 203 stress 0.06027432 
    ## Run 204 stress 0.05051364 
    ## ... Procrustes: rmse 0.03571678  max resid 0.1220412 
    ## Run 205 stress 0.05051346 
    ## ... Procrustes: rmse 0.03574365  max resid 0.1221173 
    ## Run 206 stress 0.05829938 
    ## Run 207 stress 0.05466708 
    ## Run 208 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001409605  max resid 0.0003400126 
    ## ... Similar to previous best
    ## Run 209 stress 0.06027436 
    ## Run 210 stress 0.05602058 
    ## Run 211 stress 0.05337593 
    ## Run 212 stress 0.05171058 
    ## Run 213 stress 0.05051309 
    ## ... Procrustes: rmse 0.035715  max resid 0.122017 
    ## Run 214 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001428532  max resid 0.0004179129 
    ## ... Similar to previous best
    ## Run 215 stress 0.05036794 
    ## ... Procrustes: rmse 7.913042e-05  max resid 0.0002217147 
    ## ... Similar to previous best
    ## Run 216 stress 0.0533758 
    ## Run 217 stress 0.05051307 
    ## ... Procrustes: rmse 0.03569999  max resid 0.1219718 
    ## Run 218 stress 0.05928246 
    ## Run 219 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001106395  max resid 0.0002635434 
    ## ... Similar to previous best
    ## Run 220 stress 0.05968641 
    ## Run 221 stress 0.05036801 
    ## ... Procrustes: rmse 0.000165328  max resid 0.0003925495 
    ## ... Similar to previous best
    ## Run 222 stress 0.05466702 
    ## Run 223 stress 0.05337587 
    ## Run 224 stress 0.05051318 
    ## ... Procrustes: rmse 0.03570843  max resid 0.1218554 
    ## Run 225 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001312249  max resid 0.0003497848 
    ## ... Similar to previous best
    ## Run 226 stress 0.05968626 
    ## Run 227 stress 0.05968632 
    ## Run 228 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001720767  max resid 0.0004462049 
    ## ... Similar to previous best
    ## Run 229 stress 0.05051367 
    ## ... Procrustes: rmse 0.03571135  max resid 0.1220247 
    ## Run 230 stress 0.05051356 
    ## ... Procrustes: rmse 0.03575434  max resid 0.1221527 
    ## Run 231 stress 0.0630984 
    ## Run 232 stress 0.06027446 
    ## Run 233 stress 0.06716866 
    ## Run 234 stress 0.05403606 
    ## Run 235 stress 0.0582994 
    ## Run 236 stress 0.05602053 
    ## Run 237 stress 0.05051338 
    ## ... Procrustes: rmse 0.03562378  max resid 0.1215869 
    ## Run 238 stress 0.06609076 
    ## Run 239 stress 0.06851912 
    ## Run 240 stress 0.06027436 
    ## Run 241 stress 0.05171059 
    ## Run 242 stress 0.05829934 
    ## Run 243 stress 0.05171059 
    ## Run 244 stress 0.05829933 
    ## Run 245 stress 0.05602077 
    ## Run 246 stress 0.05979235 
    ## Run 247 stress 0.05051303 
    ## ... Procrustes: rmse 0.03572111  max resid 0.12191 
    ## Run 248 stress 0.05171061 
    ## Run 249 stress 0.05171063 
    ## Run 250 stress 0.05602052 
    ## Run 251 stress 0.05036792 
    ## ... Procrustes: rmse 7.350894e-05  max resid 0.0001697856 
    ## ... Similar to previous best
    ## Run 252 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 8.470325e-06  max resid 1.55511e-05 
    ## ... Similar to previous best
    ## Run 253 stress 0.06938756 
    ## Run 254 stress 0.05171059 
    ## Run 255 stress 0.05403601 
    ## Run 256 stress 0.05403565 
    ## Run 257 stress 0.0517106 
    ## Run 258 stress 0.0505134 
    ## ... Procrustes: rmse 0.03574252  max resid 0.1221114 
    ## Run 259 stress 0.06478554 
    ## Run 260 stress 0.05466705 
    ## Run 261 stress 0.05602049 
    ## Run 262 stress 0.05171059 
    ## Run 263 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001738372  max resid 0.0004290147 
    ## ... Similar to previous best
    ## Run 264 stress 0.06111769 
    ## Run 265 stress 0.05829945 
    ## Run 266 stress 0.05466703 
    ## Run 267 stress 0.05051272 
    ## ... Procrustes: rmse 0.03566379  max resid 0.1218124 
    ## Run 268 stress 0.05051343 
    ## ... Procrustes: rmse 0.03570852  max resid 0.1218325 
    ## Run 269 stress 0.06640353 
    ## Run 270 stress 0.05171058 
    ## Run 271 stress 0.05337577 
    ## Run 272 stress 0.05036794 
    ## ... Procrustes: rmse 9.04125e-05  max resid 0.0002020445 
    ## ... Similar to previous best
    ## Run 273 stress 0.05036791 
    ## ... Procrustes: rmse 1.698152e-05  max resid 3.975261e-05 
    ## ... Similar to previous best
    ## Run 274 stress 0.05337577 
    ## Run 275 stress 0.0596863 
    ## Run 276 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001065497  max resid 0.0002817195 
    ## ... Similar to previous best
    ## Run 277 stress 0.06478532 
    ## Run 278 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001749826  max resid 0.0004341157 
    ## ... Similar to previous best
    ## Run 279 stress 0.06309865 
    ## Run 280 stress 0.05968632 
    ## Run 281 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001141824  max resid 0.0002730609 
    ## ... Similar to previous best
    ## Run 282 stress 0.0630984 
    ## Run 283 stress 0.05171062 
    ## Run 284 stress 0.0505134 
    ## ... Procrustes: rmse 0.03566046  max resid 0.1216893 
    ## Run 285 stress 0.05829931 
    ## Run 286 stress 0.06231288 
    ## Run 287 stress 0.05051305 
    ## ... Procrustes: rmse 0.0356589  max resid 0.1218454 
    ## Run 288 stress 0.05602048 
    ## Run 289 stress 0.05602053 
    ## Run 290 stress 0.05051338 
    ## ... Procrustes: rmse 0.03573595  max resid 0.1220915 
    ## Run 291 stress 0.06478544 
    ## Run 292 stress 0.05171064 
    ## Run 293 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001546696  max resid 0.0004267224 
    ## ... Similar to previous best
    ## Run 294 stress 0.05968633 
    ## Run 295 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001231981  max resid 0.0003395506 
    ## ... Similar to previous best
    ## Run 296 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001622017  max resid 0.0003921382 
    ## ... Similar to previous best
    ## Run 297 stress 0.05466701 
    ## Run 298 stress 0.06038943 
    ## Run 299 stress 0.05829937 
    ## Run 300 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001745259  max resid 0.0004207265 
    ## ... Similar to previous best
    ## Run 301 stress 0.05036794 
    ## ... Procrustes: rmse 7.12691e-05  max resid 0.0001439813 
    ## ... Similar to previous best
    ## Run 302 stress 0.05036794 
    ## ... Procrustes: rmse 8.664911e-05  max resid 0.0002067979 
    ## ... Similar to previous best
    ## Run 303 stress 0.050368 
    ## ... Procrustes: rmse 0.0001503766  max resid 0.0003585462 
    ## ... Similar to previous best
    ## Run 304 stress 0.05968625 
    ## Run 305 stress 0.07221087 
    ## Run 306 stress 0.06674551 
    ## Run 307 stress 0.05051323 
    ## ... Procrustes: rmse 0.03569258  max resid 0.121804 
    ## Run 308 stress 0.05171059 
    ## Run 309 stress 0.05403558 
    ## Run 310 stress 0.0647853 
    ## Run 311 stress 0.05171059 
    ## Run 312 stress 0.05036792 
    ## ... Procrustes: rmse 3.015781e-05  max resid 7.985613e-05 
    ## ... Similar to previous best
    ## Run 313 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001321191  max resid 0.0002731867 
    ## ... Similar to previous best
    ## Run 314 stress 0.05036793 
    ## ... Procrustes: rmse 7.053681e-05  max resid 0.0001618959 
    ## ... Similar to previous best
    ## Run 315 stress 0.05337583 
    ## Run 316 stress 0.06203039 
    ## Run 317 stress 0.0647852 
    ## Run 318 stress 0.05051334 
    ## ... Procrustes: rmse 0.03576872  max resid 0.1220336 
    ## Run 319 stress 0.05829933 
    ## Run 320 stress 0.05602046 
    ## Run 321 stress 0.05171059 
    ## Run 322 stress 0.06309865 
    ## Run 323 stress 0.06478536 
    ## Run 324 stress 0.05403575 
    ## Run 325 stress 0.06309842 
    ## Run 326 stress 0.06478535 
    ## Run 327 stress 0.06896874 
    ## Run 328 stress 0.05171059 
    ## Run 329 stress 0.06478543 
    ## Run 330 stress 0.05171064 
    ## Run 331 stress 0.0623129 
    ## Run 332 stress 0.05171059 
    ## Run 333 stress 0.05403555 
    ## Run 334 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001873794  max resid 0.0004638087 
    ## ... Similar to previous best
    ## Run 335 stress 0.05602074 
    ## Run 336 stress 0.05403581 
    ## Run 337 stress 0.05171058 
    ## Run 338 stress 0.06478519 
    ## Run 339 stress 0.05036795 
    ## ... Procrustes: rmse 9.32683e-05  max resid 0.0002627506 
    ## ... Similar to previous best
    ## Run 340 stress 0.05403588 
    ## Run 341 stress 0.06504221 
    ## Run 342 stress 0.05403567 
    ## Run 343 stress 0.05466703 
    ## Run 344 stress 0.05979231 
    ## Run 345 stress 0.05829932 
    ## Run 346 stress 0.05979231 
    ## Run 347 stress 0.05051314 
    ## ... Procrustes: rmse 0.03578432  max resid 0.1220952 
    ## Run 348 stress 0.05051358 
    ## ... Procrustes: rmse 0.03565388  max resid 0.1216572 
    ## Run 349 stress 0.05051356 
    ## ... Procrustes: rmse 0.03571793  max resid 0.1220438 
    ## Run 350 stress 0.06231295 
    ## Run 351 stress 0.06231303 
    ## Run 352 stress 0.05036792 
    ## ... Procrustes: rmse 4.077292e-05  max resid 9.215529e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.05171059 
    ## Run 354 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001222906  max resid 0.0002984627 
    ## ... Similar to previous best
    ## Run 355 stress 0.05036794 
    ## ... Procrustes: rmse 8.752849e-05  max resid 0.0002015319 
    ## ... Similar to previous best
    ## Run 356 stress 0.05051312 
    ## ... Procrustes: rmse 0.03571129  max resid 0.1220062 
    ## Run 357 stress 0.0673025 
    ## Run 358 stress 0.05051329 
    ## ... Procrustes: rmse 0.03566661  max resid 0.121718 
    ## Run 359 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001572667  max resid 0.000414986 
    ## ... Similar to previous best
    ## Run 360 stress 0.06362799 
    ## Run 361 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001275688  max resid 0.0003743407 
    ## ... Similar to previous best
    ## Run 362 stress 0.05051349 
    ## ... Procrustes: rmse 0.03573557  max resid 0.1220944 
    ## Run 363 stress 0.05171061 
    ## Run 364 stress 0.05403581 
    ## Run 365 stress 0.0505132 
    ## ... Procrustes: rmse 0.03572392  max resid 0.1220473 
    ## Run 366 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001417154  max resid 0.0004099207 
    ## ... Similar to previous best
    ## Run 367 stress 0.0517106 
    ## Run 368 stress 0.05403592 
    ## Run 369 stress 0.05036806 
    ## ... Procrustes: rmse 0.0001937385  max resid 0.0004808875 
    ## ... Similar to previous best
    ## Run 370 stress 0.06136686 
    ## Run 371 stress 0.05036795 
    ## ... Procrustes: rmse 9.903645e-05  max resid 0.0002601917 
    ## ... Similar to previous best
    ## Run 372 stress 0.05171059 
    ## Run 373 stress 0.05051349 
    ## ... Procrustes: rmse 0.03567618  max resid 0.1217375 
    ## Run 374 stress 0.06640363 
    ## Run 375 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001496864  max resid 0.0004161952 
    ## ... Similar to previous best
    ## Run 376 stress 0.06476222 
    ## Run 377 stress 0.0636279 
    ## Run 378 stress 0.05051303 
    ## ... Procrustes: rmse 0.03574835  max resid 0.1221088 
    ## Run 379 stress 0.05829935 
    ## Run 380 stress 0.05171058 
    ## Run 381 stress 0.06478522 
    ## Run 382 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001270989  max resid 0.0003120813 
    ## ... Similar to previous best
    ## Run 383 stress 0.05036793 
    ## ... Procrustes: rmse 6.782238e-05  max resid 0.0001986959 
    ## ... Similar to previous best
    ## Run 384 stress 0.05979225 
    ## Run 385 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001043156  max resid 0.0002940048 
    ## ... Similar to previous best
    ## Run 386 stress 0.06038943 
    ## Run 387 stress 0.05036811 
    ## ... Procrustes: rmse 0.0002059659  max resid 0.0005079731 
    ## ... Similar to previous best
    ## Run 388 stress 0.06231298 
    ## Run 389 stress 0.05968628 
    ## Run 390 stress 0.05171059 
    ## Run 391 stress 0.0517106 
    ## Run 392 stress 0.0613668 
    ## Run 393 stress 0.05036806 
    ## ... Procrustes: rmse 0.0001651848  max resid 0.0004048996 
    ## ... Similar to previous best
    ## Run 394 stress 0.05051322 
    ## ... Procrustes: rmse 0.03573932  max resid 0.1220952 
    ## Run 395 stress 0.05466709 
    ## Run 396 stress 0.05968627 
    ## Run 397 stress 0.05051318 
    ## ... Procrustes: rmse 0.03571802  max resid 0.12203 
    ## Run 398 stress 0.05602051 
    ## Run 399 stress 0.05051321 
    ## ... Procrustes: rmse 0.03569927  max resid 0.1219754 
    ## Run 400 stress 0.06730251 
    ## Run 401 stress 0.05171058 
    ## Run 402 stress 0.05051319 
    ## ... Procrustes: rmse 0.03573746  max resid 0.1220879 
    ## Run 403 stress 0.06027438 
    ## Run 404 stress 0.05051296 
    ## ... Procrustes: rmse 0.0356853  max resid 0.1218106 
    ## Run 405 stress 0.06493712 
    ## Run 406 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001289606  max resid 0.000368486 
    ## ... Similar to previous best
    ## Run 407 stress 0.05403591 
    ## Run 408 stress 0.06478543 
    ## Run 409 stress 0.05051342 
    ## ... Procrustes: rmse 0.03568092  max resid 0.121927 
    ## Run 410 stress 0.05051298 
    ## ... Procrustes: rmse 0.03569266  max resid 0.1219424 
    ## Run 411 stress 0.05602061 
    ## Run 412 stress 0.05968639 
    ## Run 413 stress 0.06647204 
    ## Run 414 stress 0.05337577 
    ## Run 415 stress 0.05036795 
    ## ... Procrustes: rmse 9.955775e-05  max resid 0.0002085713 
    ## ... Similar to previous best
    ## Run 416 stress 0.05928218 
    ## Run 417 stress 0.05036795 
    ## ... Procrustes: rmse 9.101895e-05  max resid 0.0002133334 
    ## ... Similar to previous best
    ## Run 418 stress 0.05051337 
    ## ... Procrustes: rmse 0.03570047  max resid 0.1219846 
    ## Run 419 stress 0.05051339 
    ## ... Procrustes: rmse 0.03564724  max resid 0.1216528 
    ## Run 420 stress 0.0533758 
    ## Run 421 stress 0.05036792 
    ## ... Procrustes: rmse 4.200241e-05  max resid 9.550762e-05 
    ## ... Similar to previous best
    ## Run 422 stress 0.06851886 
    ## Run 423 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001667168  max resid 0.0004058611 
    ## ... Similar to previous best
    ## Run 424 stress 0.05036795 
    ## ... Procrustes: rmse 9.227661e-05  max resid 0.0002484695 
    ## ... Similar to previous best
    ## Run 425 stress 0.05051302 
    ## ... Procrustes: rmse 0.03567603  max resid 0.1217769 
    ## Run 426 stress 0.05036807 
    ## ... Procrustes: rmse 0.0001928667  max resid 0.0004076985 
    ## ... Similar to previous best
    ## Run 427 stress 0.05829939 
    ## Run 428 stress 0.05337577 
    ## Run 429 stress 0.06797891 
    ## Run 430 stress 0.05968626 
    ## Run 431 stress 0.05602056 
    ## Run 432 stress 0.05968635 
    ## Run 433 stress 0.05337581 
    ## Run 434 stress 0.05968637 
    ## Run 435 stress 0.05036793 
    ## ... Procrustes: rmse 5.898942e-05  max resid 0.0001695874 
    ## ... Similar to previous best
    ## Run 436 stress 0.05171064 
    ## Run 437 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001864998  max resid 0.0004602615 
    ## ... Similar to previous best
    ## Run 438 stress 0.05466703 
    ## Run 439 stress 0.05171059 
    ## Run 440 stress 0.05829936 
    ## Run 441 stress 0.3592579 
    ## Run 442 stress 0.05051301 
    ## ... Procrustes: rmse 0.0356439  max resid 0.121794 
    ## Run 443 stress 0.0517106 
    ## Run 444 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001024366  max resid 0.0002412827 
    ## ... Similar to previous best
    ## Run 445 stress 0.06027433 
    ## Run 446 stress 0.05466705 
    ## Run 447 stress 0.05051325 
    ## ... Procrustes: rmse 0.03570144  max resid 0.1219839 
    ## Run 448 stress 0.06231288 
    ## Run 449 stress 0.05403591 
    ## Run 450 stress 0.05171059 
    ## Run 451 stress 0.0613668 
    ## Run 452 stress 0.05051313 
    ## ... Procrustes: rmse 0.03569531  max resid 0.1219597 
    ## Run 453 stress 0.05051326 
    ## ... Procrustes: rmse 0.03572441  max resid 0.1220527 
    ## Run 454 stress 0.05051285 
    ## ... Procrustes: rmse 0.03557516  max resid 0.1215603 
    ## Run 455 stress 0.05900148 
    ## Run 456 stress 0.05337581 
    ## Run 457 stress 0.05337581 
    ## Run 458 stress 0.06640355 
    ## Run 459 stress 0.07026225 
    ## Run 460 stress 0.05403627 
    ## Run 461 stress 0.05928217 
    ## Run 462 stress 0.05602051 
    ## Run 463 stress 0.05171059 
    ## Run 464 stress 0.05829938 
    ## Run 465 stress 0.05051328 
    ## ... Procrustes: rmse 0.03566679  max resid 0.1217194 
    ## Run 466 stress 0.06896907 
    ## Run 467 stress 0.05968633 
    ## Run 468 stress 0.05171059 
    ## Run 469 stress 0.3385138 
    ## Run 470 stress 0.05036793 
    ## ... Procrustes: rmse 6.575829e-05  max resid 0.0001636983 
    ## ... Similar to previous best
    ## Run 471 stress 0.05051346 
    ## ... Procrustes: rmse 0.03571961  max resid 0.1218627 
    ## Run 472 stress 0.0664034 
    ## Run 473 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001075532  max resid 0.0002783815 
    ## ... Similar to previous best
    ## Run 474 stress 0.05968629 
    ## Run 475 stress 0.05602048 
    ## Run 476 stress 0.05968626 
    ## Run 477 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001253176  max resid 0.000326482 
    ## ... Similar to previous best
    ## Run 478 stress 0.05171061 
    ## Run 479 stress 0.05928255 
    ## Run 480 stress 0.05968628 
    ## Run 481 stress 0.06609065 
    ## Run 482 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001663528  max resid 0.0004035775 
    ## ... Similar to previous best
    ## Run 483 stress 0.05466705 
    ## Run 484 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001194043  max resid 0.000287922 
    ## ... Similar to previous best
    ## Run 485 stress 0.05036793 
    ## ... Procrustes: rmse 6.592213e-05  max resid 0.0001736586 
    ## ... Similar to previous best
    ## Run 486 stress 0.05171064 
    ## Run 487 stress 0.05036794 
    ## ... Procrustes: rmse 8.691658e-05  max resid 0.0002064089 
    ## ... Similar to previous best
    ## Run 488 stress 0.06027443 
    ## Run 489 stress 0.05051308 
    ## ... Procrustes: rmse 0.03570874  max resid 0.1219974 
    ## Run 490 stress 0.05337582 
    ## Run 491 stress 0.06826681 
    ## Run 492 stress 0.06027438 
    ## Run 493 stress 0.05403595 
    ## Run 494 stress 0.06478541 
    ## Run 495 stress 0.05051306 
    ## ... Procrustes: rmse 0.03571103  max resid 0.122003 
    ## Run 496 stress 0.05051334 
    ## ... Procrustes: rmse 0.03569308  max resid 0.1217962 
    ## Run 497 stress 0.05036793 
    ## ... Procrustes: rmse 7.058323e-05  max resid 0.0001651415 
    ## ... Similar to previous best
    ## Run 498 stress 0.05466705 
    ## Run 499 stress 0.05051282 
    ## ... Procrustes: rmse 0.03571184  max resid 0.1219125 
    ## Run 500 stress 0.05051319 
    ## ... Procrustes: rmse 0.03567577  max resid 0.1217574 
    ## *** Best solution repeated 51 times

``` r
# Ocean sites and mixed lakes
PD_beta_env_OM_NMDS <- metaMDS(PD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08352634 
    ## Run 1 stress 0.1279309 
    ## Run 2 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 4.755506e-06  max resid 1.052055e-05 
    ## ... Similar to previous best
    ## Run 3 stress 0.1134014 
    ## Run 4 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.482606e-06  max resid 3.100743e-06 
    ## ... Similar to previous best
    ## Run 5 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.453976e-06  max resid 2.557635e-06 
    ## ... Similar to previous best
    ## Run 6 stress 0.2610955 
    ## Run 7 stress 0.1279309 
    ## Run 8 stress 0.1134015 
    ## Run 9 stress 0.08352634 
    ## ... Procrustes: rmse 2.955796e-06  max resid 6.988218e-06 
    ## ... Similar to previous best
    ## Run 10 stress 0.1134014 
    ## Run 11 stress 0.08352634 
    ## ... Procrustes: rmse 2.442014e-06  max resid 5.974019e-06 
    ## ... Similar to previous best
    ## Run 12 stress 0.08352634 
    ## ... Procrustes: rmse 1.503051e-06  max resid 2.636763e-06 
    ## ... Similar to previous best
    ## Run 13 stress 0.08352634 
    ## ... Procrustes: rmse 1.175381e-06  max resid 2.301604e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.1279309 
    ## Run 15 stress 0.1279309 
    ## Run 16 stress 0.1279309 
    ## Run 17 stress 0.1279309 
    ## Run 18 stress 0.1134015 
    ## Run 19 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.578539e-06  max resid 2.896841e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.1279309 
    ## Run 21 stress 0.08352634 
    ## ... Procrustes: rmse 3.770545e-06  max resid 8.395498e-06 
    ## ... Similar to previous best
    ## Run 22 stress 0.08352634 
    ## ... Procrustes: rmse 1.086606e-06  max resid 2.29256e-06 
    ## ... Similar to previous best
    ## Run 23 stress 0.127931 
    ## Run 24 stress 0.1134014 
    ## Run 25 stress 0.08352634 
    ## ... Procrustes: rmse 1.345387e-06  max resid 2.794846e-06 
    ## ... Similar to previous best
    ## Run 26 stress 0.127931 
    ## Run 27 stress 0.08352634 
    ## ... Procrustes: rmse 1.352563e-06  max resid 2.151648e-06 
    ## ... Similar to previous best
    ## Run 28 stress 0.178315 
    ## Run 29 stress 0.08352634 
    ## ... Procrustes: rmse 2.573798e-06  max resid 5.634763e-06 
    ## ... Similar to previous best
    ## Run 30 stress 0.127931 
    ## Run 31 stress 0.08352634 
    ## ... Procrustes: rmse 6.20439e-07  max resid 1.097319e-06 
    ## ... Similar to previous best
    ## Run 32 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 5.518553e-07  max resid 9.486661e-07 
    ## ... Similar to previous best
    ## Run 33 stress 0.08352634 
    ## ... Procrustes: rmse 3.792793e-06  max resid 8.583884e-06 
    ## ... Similar to previous best
    ## Run 34 stress 0.1279309 
    ## Run 35 stress 0.08352634 
    ## ... Procrustes: rmse 2.70268e-06  max resid 6.086028e-06 
    ## ... Similar to previous best
    ## Run 36 stress 0.08352634 
    ## ... Procrustes: rmse 2.257375e-06  max resid 4.413045e-06 
    ## ... Similar to previous best
    ## Run 37 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.082893e-06  max resid 2.349966e-06 
    ## ... Similar to previous best
    ## Run 38 stress 0.1134014 
    ## Run 39 stress 0.08352634 
    ## ... Procrustes: rmse 1.56063e-06  max resid 3.132391e-06 
    ## ... Similar to previous best
    ## Run 40 stress 0.08352634 
    ## ... Procrustes: rmse 1.675623e-06  max resid 3.661379e-06 
    ## ... Similar to previous best
    ## Run 41 stress 0.127931 
    ## Run 42 stress 0.08352634 
    ## ... Procrustes: rmse 1.698406e-06  max resid 3.7679e-06 
    ## ... Similar to previous best
    ## Run 43 stress 0.1279309 
    ## Run 44 stress 0.08352634 
    ## ... Procrustes: rmse 1.25927e-06  max resid 2.7783e-06 
    ## ... Similar to previous best
    ## Run 45 stress 0.08352634 
    ## ... Procrustes: rmse 1.114322e-06  max resid 2.344939e-06 
    ## ... Similar to previous best
    ## Run 46 stress 0.08352634 
    ## ... Procrustes: rmse 1.85469e-06  max resid 4.223507e-06 
    ## ... Similar to previous best
    ## Run 47 stress 0.3131414 
    ## Run 48 stress 0.08352634 
    ## ... Procrustes: rmse 1.73455e-06  max resid 3.584661e-06 
    ## ... Similar to previous best
    ## Run 49 stress 0.1134015 
    ## Run 50 stress 0.08352634 
    ## ... Procrustes: rmse 1.267881e-06  max resid 2.676771e-06 
    ## ... Similar to previous best
    ## Run 51 stress 0.08352634 
    ## ... Procrustes: rmse 1.66838e-06  max resid 3.70231e-06 
    ## ... Similar to previous best
    ## Run 52 stress 0.1134014 
    ## Run 53 stress 0.08352634 
    ## ... Procrustes: rmse 1.854489e-06  max resid 2.849592e-06 
    ## ... Similar to previous best
    ## Run 54 stress 0.08352634 
    ## ... Procrustes: rmse 5.323623e-07  max resid 1.163236e-06 
    ## ... Similar to previous best
    ## Run 55 stress 0.08352634 
    ## ... Procrustes: rmse 3.434084e-06  max resid 6.835152e-06 
    ## ... Similar to previous best
    ## Run 56 stress 0.08352634 
    ## ... Procrustes: rmse 2.52687e-06  max resid 5.62473e-06 
    ## ... Similar to previous best
    ## Run 57 stress 0.1279309 
    ## Run 58 stress 0.1134014 
    ## Run 59 stress 0.08352634 
    ## ... Procrustes: rmse 2.155486e-06  max resid 4.826017e-06 
    ## ... Similar to previous best
    ## Run 60 stress 0.1134015 
    ## Run 61 stress 0.08352634 
    ## ... Procrustes: rmse 3.243279e-07  max resid 6.823714e-07 
    ## ... Similar to previous best
    ## Run 62 stress 0.08352634 
    ## ... Procrustes: rmse 2.879214e-06  max resid 6.289864e-06 
    ## ... Similar to previous best
    ## Run 63 stress 0.08352634 
    ## ... Procrustes: rmse 1.660071e-06  max resid 3.593811e-06 
    ## ... Similar to previous best
    ## Run 64 stress 0.1134014 
    ## Run 65 stress 0.08352634 
    ## ... Procrustes: rmse 9.415106e-06  max resid 1.47654e-05 
    ## ... Similar to previous best
    ## Run 66 stress 0.1134014 
    ## Run 67 stress 0.3065578 
    ## Run 68 stress 0.08352634 
    ## ... Procrustes: rmse 4.283408e-06  max resid 9.483448e-06 
    ## ... Similar to previous best
    ## Run 69 stress 0.08352634 
    ## ... Procrustes: rmse 2.92912e-06  max resid 6.476508e-06 
    ## ... Similar to previous best
    ## Run 70 stress 0.2938011 
    ## Run 71 stress 0.08352634 
    ## ... Procrustes: rmse 7.89442e-07  max resid 1.498581e-06 
    ## ... Similar to previous best
    ## Run 72 stress 0.08352634 
    ## ... Procrustes: rmse 1.369264e-06  max resid 3.033722e-06 
    ## ... Similar to previous best
    ## Run 73 stress 0.08352634 
    ## ... Procrustes: rmse 3.134172e-06  max resid 7.005696e-06 
    ## ... Similar to previous best
    ## Run 74 stress 0.08352634 
    ## ... Procrustes: rmse 2.382839e-06  max resid 5.397584e-06 
    ## ... Similar to previous best
    ## Run 75 stress 0.08352634 
    ## ... Procrustes: rmse 3.851091e-06  max resid 8.442662e-06 
    ## ... Similar to previous best
    ## Run 76 stress 0.08352634 
    ## ... Procrustes: rmse 3.863426e-06  max resid 8.783378e-06 
    ## ... Similar to previous best
    ## Run 77 stress 0.1134014 
    ## Run 78 stress 0.1279309 
    ## Run 79 stress 0.08352634 
    ## ... Procrustes: rmse 1.056589e-06  max resid 2.217373e-06 
    ## ... Similar to previous best
    ## Run 80 stress 0.1134014 
    ## Run 81 stress 0.1134015 
    ## Run 82 stress 0.08352634 
    ## ... Procrustes: rmse 3.485149e-06  max resid 7.194148e-06 
    ## ... Similar to previous best
    ## Run 83 stress 0.1989489 
    ## Run 84 stress 0.1279309 
    ## Run 85 stress 0.08352634 
    ## ... Procrustes: rmse 1.680621e-06  max resid 3.839582e-06 
    ## ... Similar to previous best
    ## Run 86 stress 0.1279309 
    ## Run 87 stress 0.08352634 
    ## ... Procrustes: rmse 6.562666e-07  max resid 1.354902e-06 
    ## ... Similar to previous best
    ## Run 88 stress 0.1134014 
    ## Run 89 stress 0.08352634 
    ## ... Procrustes: rmse 2.473703e-06  max resid 5.502723e-06 
    ## ... Similar to previous best
    ## Run 90 stress 0.1279309 
    ## Run 91 stress 0.1134014 
    ## Run 92 stress 0.1134014 
    ## Run 93 stress 0.08352634 
    ## ... Procrustes: rmse 1.722006e-06  max resid 3.853821e-06 
    ## ... Similar to previous best
    ## Run 94 stress 0.08352634 
    ## ... Procrustes: rmse 4.247042e-06  max resid 9.313664e-06 
    ## ... Similar to previous best
    ## Run 95 stress 0.08352634 
    ## ... Procrustes: rmse 9.784867e-07  max resid 2.025662e-06 
    ## ... Similar to previous best
    ## Run 96 stress 0.1989158 
    ## Run 97 stress 0.08352634 
    ## ... Procrustes: rmse 6.934755e-07  max resid 1.340374e-06 
    ## ... Similar to previous best
    ## Run 98 stress 0.08352634 
    ## ... Procrustes: rmse 3.705007e-06  max resid 8.15567e-06 
    ## ... Similar to previous best
    ## Run 99 stress 0.1279309 
    ## Run 100 stress 0.08352634 
    ## ... Procrustes: rmse 9.226367e-07  max resid 1.924028e-06 
    ## ... Similar to previous best
    ## Run 101 stress 0.08352634 
    ## ... Procrustes: rmse 2.189042e-06  max resid 3.77681e-06 
    ## ... Similar to previous best
    ## Run 102 stress 0.1134015 
    ## Run 103 stress 0.08352634 
    ## ... Procrustes: rmse 2.649475e-06  max resid 4.986638e-06 
    ## ... Similar to previous best
    ## Run 104 stress 0.08352634 
    ## ... Procrustes: rmse 1.459457e-06  max resid 3.027443e-06 
    ## ... Similar to previous best
    ## Run 105 stress 0.1279309 
    ## Run 106 stress 0.08352634 
    ## ... Procrustes: rmse 3.036196e-06  max resid 6.770638e-06 
    ## ... Similar to previous best
    ## Run 107 stress 0.08352634 
    ## ... Procrustes: rmse 7.370771e-07  max resid 1.229966e-06 
    ## ... Similar to previous best
    ## Run 108 stress 0.1134014 
    ## Run 109 stress 0.1134014 
    ## Run 110 stress 0.1134014 
    ## Run 111 stress 0.08352634 
    ## ... Procrustes: rmse 2.274395e-06  max resid 3.488254e-06 
    ## ... Similar to previous best
    ## Run 112 stress 0.08352634 
    ## ... Procrustes: rmse 3.458724e-06  max resid 7.705401e-06 
    ## ... Similar to previous best
    ## Run 113 stress 0.08352634 
    ## ... Procrustes: rmse 6.48002e-07  max resid 1.309372e-06 
    ## ... Similar to previous best
    ## Run 114 stress 0.08352634 
    ## ... Procrustes: rmse 1.824841e-06  max resid 3.986719e-06 
    ## ... Similar to previous best
    ## Run 115 stress 0.127931 
    ## Run 116 stress 0.08352634 
    ## ... Procrustes: rmse 3.741894e-06  max resid 7.859652e-06 
    ## ... Similar to previous best
    ## Run 117 stress 0.08352634 
    ## ... Procrustes: rmse 2.234013e-06  max resid 4.914901e-06 
    ## ... Similar to previous best
    ## Run 118 stress 0.1134015 
    ## Run 119 stress 0.08352634 
    ## ... Procrustes: rmse 1.644396e-06  max resid 3.605733e-06 
    ## ... Similar to previous best
    ## Run 120 stress 0.1134015 
    ## Run 121 stress 0.08352634 
    ## ... Procrustes: rmse 5.531341e-07  max resid 7.893899e-07 
    ## ... Similar to previous best
    ## Run 122 stress 0.08352634 
    ## ... Procrustes: rmse 2.226963e-06  max resid 4.872237e-06 
    ## ... Similar to previous best
    ## Run 123 stress 0.1279309 
    ## Run 124 stress 0.127931 
    ## Run 125 stress 0.08352634 
    ## ... Procrustes: rmse 3.054646e-06  max resid 6.498966e-06 
    ## ... Similar to previous best
    ## Run 126 stress 0.08352634 
    ## ... Procrustes: rmse 4.321166e-06  max resid 9.856269e-06 
    ## ... Similar to previous best
    ## Run 127 stress 0.1134015 
    ## Run 128 stress 0.1134014 
    ## Run 129 stress 0.08352634 
    ## ... Procrustes: rmse 7.002272e-06  max resid 1.548302e-05 
    ## ... Similar to previous best
    ## Run 130 stress 0.08352634 
    ## ... Procrustes: rmse 6.532762e-07  max resid 1.109467e-06 
    ## ... Similar to previous best
    ## Run 131 stress 0.08352634 
    ## ... Procrustes: rmse 5.492889e-07  max resid 9.482366e-07 
    ## ... Similar to previous best
    ## Run 132 stress 0.1279309 
    ## Run 133 stress 0.08352634 
    ## ... Procrustes: rmse 4.770488e-07  max resid 9.469454e-07 
    ## ... Similar to previous best
    ## Run 134 stress 0.08352634 
    ## ... Procrustes: rmse 2.116227e-06  max resid 4.715304e-06 
    ## ... Similar to previous best
    ## Run 135 stress 0.08352634 
    ## ... Procrustes: rmse 4.476932e-06  max resid 9.871793e-06 
    ## ... Similar to previous best
    ## Run 136 stress 0.08352634 
    ## ... Procrustes: rmse 3.158538e-06  max resid 7.107259e-06 
    ## ... Similar to previous best
    ## Run 137 stress 0.1279309 
    ## Run 138 stress 0.1134015 
    ## Run 139 stress 0.1989488 
    ## Run 140 stress 0.1279309 
    ## Run 141 stress 0.08352634 
    ## ... Procrustes: rmse 1.959736e-06  max resid 4.025458e-06 
    ## ... Similar to previous best
    ## Run 142 stress 0.1134015 
    ## Run 143 stress 0.1134015 
    ## Run 144 stress 0.1279309 
    ## Run 145 stress 0.08352634 
    ## ... Procrustes: rmse 5.028603e-07  max resid 9.916292e-07 
    ## ... Similar to previous best
    ## Run 146 stress 0.08352634 
    ## ... Procrustes: rmse 1.957499e-06  max resid 3.708162e-06 
    ## ... Similar to previous best
    ## Run 147 stress 0.3275202 
    ## Run 148 stress 0.1279309 
    ## Run 149 stress 0.08352634 
    ## ... Procrustes: rmse 1.477608e-06  max resid 3.068615e-06 
    ## ... Similar to previous best
    ## Run 150 stress 0.08352634 
    ## ... Procrustes: rmse 2.264977e-06  max resid 5.029509e-06 
    ## ... Similar to previous best
    ## Run 151 stress 0.1134014 
    ## Run 152 stress 0.1279309 
    ## Run 153 stress 0.1134014 
    ## Run 154 stress 0.1134014 
    ## Run 155 stress 0.08352634 
    ## ... Procrustes: rmse 2.678334e-06  max resid 5.802545e-06 
    ## ... Similar to previous best
    ## Run 156 stress 0.1279309 
    ## Run 157 stress 0.08352634 
    ## ... Procrustes: rmse 9.089264e-07  max resid 2.009843e-06 
    ## ... Similar to previous best
    ## Run 158 stress 0.1134014 
    ## Run 159 stress 0.08352634 
    ## ... Procrustes: rmse 1.426803e-06  max resid 2.914598e-06 
    ## ... Similar to previous best
    ## Run 160 stress 0.1279309 
    ## Run 161 stress 0.08352634 
    ## ... Procrustes: rmse 1.483348e-06  max resid 3.299467e-06 
    ## ... Similar to previous best
    ## Run 162 stress 0.08352634 
    ## ... Procrustes: rmse 3.480602e-06  max resid 7.757276e-06 
    ## ... Similar to previous best
    ## Run 163 stress 0.08352634 
    ## ... Procrustes: rmse 3.964177e-07  max resid 8.578055e-07 
    ## ... Similar to previous best
    ## Run 164 stress 0.1279309 
    ## Run 165 stress 0.08352634 
    ## ... Procrustes: rmse 7.831523e-07  max resid 1.229721e-06 
    ## ... Similar to previous best
    ## Run 166 stress 0.1134015 
    ## Run 167 stress 0.08352634 
    ## ... Procrustes: rmse 2.431033e-06  max resid 5.062327e-06 
    ## ... Similar to previous best
    ## Run 168 stress 0.1134014 
    ## Run 169 stress 0.08352634 
    ## ... Procrustes: rmse 2.933506e-06  max resid 6.490659e-06 
    ## ... Similar to previous best
    ## Run 170 stress 0.1279309 
    ## Run 171 stress 0.08352634 
    ## ... Procrustes: rmse 2.821529e-06  max resid 6.332724e-06 
    ## ... Similar to previous best
    ## Run 172 stress 0.08352634 
    ## ... Procrustes: rmse 2.489195e-06  max resid 5.573405e-06 
    ## ... Similar to previous best
    ## Run 173 stress 0.1134015 
    ## Run 174 stress 0.08352634 
    ## ... Procrustes: rmse 2.313841e-06  max resid 5.189288e-06 
    ## ... Similar to previous best
    ## Run 175 stress 0.1134014 
    ## Run 176 stress 0.08352634 
    ## ... Procrustes: rmse 9.738236e-07  max resid 1.618831e-06 
    ## ... Similar to previous best
    ## Run 177 stress 0.08352634 
    ## ... Procrustes: rmse 3.550222e-06  max resid 7.954966e-06 
    ## ... Similar to previous best
    ## Run 178 stress 0.1134015 
    ## Run 179 stress 0.08352634 
    ## ... Procrustes: rmse 9.953302e-07  max resid 1.870178e-06 
    ## ... Similar to previous best
    ## Run 180 stress 0.08352634 
    ## ... Procrustes: rmse 9.008063e-07  max resid 1.553134e-06 
    ## ... Similar to previous best
    ## Run 181 stress 0.08352634 
    ## ... Procrustes: rmse 4.659482e-06  max resid 1.035792e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.08352634 
    ## ... Procrustes: rmse 4.547483e-06  max resid 1.036972e-05 
    ## ... Similar to previous best
    ## Run 183 stress 0.08352634 
    ## ... Procrustes: rmse 2.332509e-06  max resid 5.18734e-06 
    ## ... Similar to previous best
    ## Run 184 stress 0.08352634 
    ## ... Procrustes: rmse 2.204919e-06  max resid 4.897308e-06 
    ## ... Similar to previous best
    ## Run 185 stress 0.08352634 
    ## ... Procrustes: rmse 2.154451e-06  max resid 4.939798e-06 
    ## ... Similar to previous best
    ## Run 186 stress 0.08352634 
    ## ... Procrustes: rmse 2.49453e-06  max resid 5.781078e-06 
    ## ... Similar to previous best
    ## Run 187 stress 0.1134014 
    ## Run 188 stress 0.08352634 
    ## ... Procrustes: rmse 3.497681e-06  max resid 7.676088e-06 
    ## ... Similar to previous best
    ## Run 189 stress 0.1134014 
    ## Run 190 stress 0.08352634 
    ## ... Procrustes: rmse 2.468589e-06  max resid 5.627053e-06 
    ## ... Similar to previous best
    ## Run 191 stress 0.08352634 
    ## ... Procrustes: rmse 9.24942e-07  max resid 2.032774e-06 
    ## ... Similar to previous best
    ## Run 192 stress 0.08352634 
    ## ... Procrustes: rmse 3.177628e-06  max resid 7.354868e-06 
    ## ... Similar to previous best
    ## Run 193 stress 0.08352634 
    ## ... Procrustes: rmse 5.126618e-06  max resid 1.068466e-05 
    ## ... Similar to previous best
    ## Run 194 stress 0.1134015 
    ## Run 195 stress 0.127931 
    ## Run 196 stress 0.08352634 
    ## ... Procrustes: rmse 1.575563e-06  max resid 3.498521e-06 
    ## ... Similar to previous best
    ## Run 197 stress 0.1989157 
    ## Run 198 stress 0.1279309 
    ## Run 199 stress 0.08352634 
    ## ... Procrustes: rmse 7.276923e-07  max resid 1.565006e-06 
    ## ... Similar to previous best
    ## Run 200 stress 0.08352634 
    ## ... Procrustes: rmse 2.197261e-06  max resid 4.834045e-06 
    ## ... Similar to previous best
    ## Run 201 stress 0.1134014 
    ## Run 202 stress 0.1134014 
    ## Run 203 stress 0.1134014 
    ## Run 204 stress 0.2610955 
    ## Run 205 stress 0.1989156 
    ## Run 206 stress 0.08352634 
    ## ... Procrustes: rmse 2.263407e-06  max resid 4.947588e-06 
    ## ... Similar to previous best
    ## Run 207 stress 0.08352634 
    ## ... Procrustes: rmse 2.904934e-06  max resid 6.532119e-06 
    ## ... Similar to previous best
    ## Run 208 stress 0.1989156 
    ## Run 209 stress 0.1134014 
    ## Run 210 stress 0.08352634 
    ## ... Procrustes: rmse 1.910843e-06  max resid 4.302457e-06 
    ## ... Similar to previous best
    ## Run 211 stress 0.127931 
    ## Run 212 stress 0.08352634 
    ## ... Procrustes: rmse 1.767897e-06  max resid 3.918942e-06 
    ## ... Similar to previous best
    ## Run 213 stress 0.08352634 
    ## ... Procrustes: rmse 8.567019e-07  max resid 1.350163e-06 
    ## ... Similar to previous best
    ## Run 214 stress 0.08352634 
    ## ... Procrustes: rmse 3.332655e-06  max resid 7.50912e-06 
    ## ... Similar to previous best
    ## Run 215 stress 0.08352634 
    ## ... Procrustes: rmse 1.68202e-06  max resid 3.742994e-06 
    ## ... Similar to previous best
    ## Run 216 stress 0.1783132 
    ## Run 217 stress 0.1279309 
    ## Run 218 stress 0.08352634 
    ## ... Procrustes: rmse 5.405006e-06  max resid 1.053616e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.08352634 
    ## ... Procrustes: rmse 7.855201e-07  max resid 1.174968e-06 
    ## ... Similar to previous best
    ## Run 220 stress 0.08352634 
    ## ... Procrustes: rmse 6.616122e-07  max resid 1.059301e-06 
    ## ... Similar to previous best
    ## Run 221 stress 0.08352634 
    ## ... Procrustes: rmse 2.563626e-06  max resid 5.662982e-06 
    ## ... Similar to previous best
    ## Run 222 stress 0.1134015 
    ## Run 223 stress 0.1279309 
    ## Run 224 stress 0.08352634 
    ## ... Procrustes: rmse 3.697825e-07  max resid 6.322856e-07 
    ## ... Similar to previous best
    ## Run 225 stress 0.08352634 
    ## ... Procrustes: rmse 3.007448e-06  max resid 6.565669e-06 
    ## ... Similar to previous best
    ## Run 226 stress 0.08352634 
    ## ... Procrustes: rmse 1.600446e-06  max resid 3.364281e-06 
    ## ... Similar to previous best
    ## Run 227 stress 0.1279309 
    ## Run 228 stress 0.08352634 
    ## ... Procrustes: rmse 2.312754e-06  max resid 5.206719e-06 
    ## ... Similar to previous best
    ## Run 229 stress 0.08352634 
    ## ... Procrustes: rmse 1.215244e-06  max resid 2.600711e-06 
    ## ... Similar to previous best
    ## Run 230 stress 0.1134014 
    ## Run 231 stress 0.08352634 
    ## ... Procrustes: rmse 1.36405e-06  max resid 2.928079e-06 
    ## ... Similar to previous best
    ## Run 232 stress 0.1279309 
    ## Run 233 stress 0.08352634 
    ## ... Procrustes: rmse 1.455301e-06  max resid 3.420064e-06 
    ## ... Similar to previous best
    ## Run 234 stress 0.1134015 
    ## Run 235 stress 0.1134015 
    ## Run 236 stress 0.08352634 
    ## ... Procrustes: rmse 8.556982e-06  max resid 1.804347e-05 
    ## ... Similar to previous best
    ## Run 237 stress 0.08352634 
    ## ... Procrustes: rmse 3.636774e-06  max resid 7.5079e-06 
    ## ... Similar to previous best
    ## Run 238 stress 0.08352634 
    ## ... Procrustes: rmse 1.984369e-06  max resid 4.330357e-06 
    ## ... Similar to previous best
    ## Run 239 stress 0.08352634 
    ## ... Procrustes: rmse 2.421923e-06  max resid 5.325056e-06 
    ## ... Similar to previous best
    ## Run 240 stress 0.08352634 
    ## ... Procrustes: rmse 2.991448e-06  max resid 6.682966e-06 
    ## ... Similar to previous best
    ## Run 241 stress 0.1279309 
    ## Run 242 stress 0.08352634 
    ## ... Procrustes: rmse 3.370032e-06  max resid 7.433344e-06 
    ## ... Similar to previous best
    ## Run 243 stress 0.08352634 
    ## ... Procrustes: rmse 1.679161e-06  max resid 3.869546e-06 
    ## ... Similar to previous best
    ## Run 244 stress 0.08352634 
    ## ... Procrustes: rmse 1.0685e-06  max resid 1.98689e-06 
    ## ... Similar to previous best
    ## Run 245 stress 0.1134014 
    ## Run 246 stress 0.08352634 
    ## ... Procrustes: rmse 2.150461e-06  max resid 4.739456e-06 
    ## ... Similar to previous best
    ## Run 247 stress 0.08352634 
    ## ... Procrustes: rmse 2.723865e-06  max resid 5.08745e-06 
    ## ... Similar to previous best
    ## Run 248 stress 0.1989159 
    ## Run 249 stress 0.1134014 
    ## Run 250 stress 0.1134015 
    ## Run 251 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 5.536812e-07  max resid 9.822887e-07 
    ## ... Similar to previous best
    ## Run 252 stress 0.08352634 
    ## ... Procrustes: rmse 5.997265e-07  max resid 1.004031e-06 
    ## ... Similar to previous best
    ## Run 253 stress 0.08352634 
    ## ... Procrustes: rmse 8.014201e-07  max resid 1.575461e-06 
    ## ... Similar to previous best
    ## Run 254 stress 0.08352634 
    ## ... Procrustes: rmse 3.526301e-06  max resid 5.655398e-06 
    ## ... Similar to previous best
    ## Run 255 stress 0.1134015 
    ## Run 256 stress 0.08352634 
    ## ... Procrustes: rmse 6.647625e-07  max resid 1.321104e-06 
    ## ... Similar to previous best
    ## Run 257 stress 0.3266727 
    ## Run 258 stress 0.08352634 
    ## ... Procrustes: rmse 8.75049e-07  max resid 1.841808e-06 
    ## ... Similar to previous best
    ## Run 259 stress 0.08352634 
    ## ... Procrustes: rmse 3.434216e-06  max resid 7.21336e-06 
    ## ... Similar to previous best
    ## Run 260 stress 0.205005 
    ## Run 261 stress 0.1134014 
    ## Run 262 stress 0.1279309 
    ## Run 263 stress 0.08352634 
    ## ... Procrustes: rmse 3.591126e-06  max resid 8.048565e-06 
    ## ... Similar to previous best
    ## Run 264 stress 0.08352634 
    ## ... Procrustes: rmse 1.848317e-06  max resid 4.435873e-06 
    ## ... Similar to previous best
    ## Run 265 stress 0.1279309 
    ## Run 266 stress 0.1783136 
    ## Run 267 stress 0.1279309 
    ## Run 268 stress 0.1279309 
    ## Run 269 stress 0.1134014 
    ## Run 270 stress 0.3072351 
    ## Run 271 stress 0.1989157 
    ## Run 272 stress 0.1134014 
    ## Run 273 stress 0.1279309 
    ## Run 274 stress 0.08352634 
    ## ... Procrustes: rmse 1.664776e-06  max resid 3.734155e-06 
    ## ... Similar to previous best
    ## Run 275 stress 0.08352634 
    ## ... Procrustes: rmse 2.714782e-06  max resid 6.000879e-06 
    ## ... Similar to previous best
    ## Run 276 stress 0.1134015 
    ## Run 277 stress 0.08352634 
    ## ... Procrustes: rmse 1.804417e-06  max resid 4.173304e-06 
    ## ... Similar to previous best
    ## Run 278 stress 0.1279309 
    ## Run 279 stress 0.1134014 
    ## Run 280 stress 0.1279309 
    ## Run 281 stress 0.08352634 
    ## ... Procrustes: rmse 5.950466e-07  max resid 9.803857e-07 
    ## ... Similar to previous best
    ## Run 282 stress 0.1134015 
    ## Run 283 stress 0.1134014 
    ## Run 284 stress 0.08352634 
    ## ... Procrustes: rmse 2.42743e-06  max resid 5.398847e-06 
    ## ... Similar to previous best
    ## Run 285 stress 0.08352634 
    ## ... Procrustes: rmse 1.293065e-06  max resid 2.808838e-06 
    ## ... Similar to previous best
    ## Run 286 stress 0.1134014 
    ## Run 287 stress 0.08352634 
    ## ... Procrustes: rmse 4.009618e-06  max resid 8.953739e-06 
    ## ... Similar to previous best
    ## Run 288 stress 0.08352634 
    ## ... Procrustes: rmse 2.26023e-06  max resid 4.981687e-06 
    ## ... Similar to previous best
    ## Run 289 stress 0.08352634 
    ## ... Procrustes: rmse 6.75516e-07  max resid 1.347711e-06 
    ## ... Similar to previous best
    ## Run 290 stress 0.1134014 
    ## Run 291 stress 0.3061932 
    ## Run 292 stress 0.08352634 
    ## ... Procrustes: rmse 2.689743e-06  max resid 5.50377e-06 
    ## ... Similar to previous best
    ## Run 293 stress 0.08352634 
    ## ... Procrustes: rmse 3.671482e-06  max resid 8.101933e-06 
    ## ... Similar to previous best
    ## Run 294 stress 0.08352634 
    ## ... Procrustes: rmse 3.393598e-06  max resid 6.881368e-06 
    ## ... Similar to previous best
    ## Run 295 stress 0.08352634 
    ## ... Procrustes: rmse 9.844306e-07  max resid 1.792067e-06 
    ## ... Similar to previous best
    ## Run 296 stress 0.127931 
    ## Run 297 stress 0.1134015 
    ## Run 298 stress 0.127931 
    ## Run 299 stress 0.08352634 
    ## ... Procrustes: rmse 1.163861e-06  max resid 2.251678e-06 
    ## ... Similar to previous best
    ## Run 300 stress 0.08352634 
    ## ... Procrustes: rmse 2.262413e-06  max resid 4.704332e-06 
    ## ... Similar to previous best
    ## Run 301 stress 0.08352634 
    ## ... Procrustes: rmse 1.597814e-06  max resid 3.319481e-06 
    ## ... Similar to previous best
    ## Run 302 stress 0.08352634 
    ## ... Procrustes: rmse 4.349552e-06  max resid 9.99012e-06 
    ## ... Similar to previous best
    ## Run 303 stress 0.3154651 
    ## Run 304 stress 0.1134015 
    ## Run 305 stress 0.08352634 
    ## ... Procrustes: rmse 1.300699e-06  max resid 2.802618e-06 
    ## ... Similar to previous best
    ## Run 306 stress 0.08352634 
    ## ... Procrustes: rmse 1.822034e-06  max resid 3.926368e-06 
    ## ... Similar to previous best
    ## Run 307 stress 0.1279309 
    ## Run 308 stress 0.1279309 
    ## Run 309 stress 0.1783108 
    ## Run 310 stress 0.3266727 
    ## Run 311 stress 0.08352634 
    ## ... Procrustes: rmse 1.568642e-06  max resid 3.728372e-06 
    ## ... Similar to previous best
    ## Run 312 stress 0.08352634 
    ## ... Procrustes: rmse 7.980504e-07  max resid 1.665811e-06 
    ## ... Similar to previous best
    ## Run 313 stress 0.08352634 
    ## ... Procrustes: rmse 1.008747e-06  max resid 2.27047e-06 
    ## ... Similar to previous best
    ## Run 314 stress 0.08352634 
    ## ... Procrustes: rmse 8.759367e-07  max resid 1.783361e-06 
    ## ... Similar to previous best
    ## Run 315 stress 0.08352634 
    ## ... Procrustes: rmse 3.58204e-06  max resid 7.869741e-06 
    ## ... Similar to previous best
    ## Run 316 stress 0.1279309 
    ## Run 317 stress 0.1279309 
    ## Run 318 stress 0.08352634 
    ## ... Procrustes: rmse 3.15616e-06  max resid 6.865471e-06 
    ## ... Similar to previous best
    ## Run 319 stress 0.08352634 
    ## ... Procrustes: rmse 1.916465e-06  max resid 4.299233e-06 
    ## ... Similar to previous best
    ## Run 320 stress 0.08352634 
    ## ... Procrustes: rmse 5.56235e-06  max resid 1.271181e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.08352634 
    ## ... Procrustes: rmse 1.808825e-06  max resid 3.868927e-06 
    ## ... Similar to previous best
    ## Run 322 stress 0.08352634 
    ## ... Procrustes: rmse 3.897051e-06  max resid 8.517005e-06 
    ## ... Similar to previous best
    ## Run 323 stress 0.1279309 
    ## Run 324 stress 0.1134014 
    ## Run 325 stress 0.1279309 
    ## Run 326 stress 0.1989484 
    ## Run 327 stress 0.08352634 
    ## ... Procrustes: rmse 4.312197e-06  max resid 9.49996e-06 
    ## ... Similar to previous best
    ## Run 328 stress 0.1134014 
    ## Run 329 stress 0.3216389 
    ## Run 330 stress 0.08352634 
    ## ... Procrustes: rmse 2.263733e-06  max resid 5.044781e-06 
    ## ... Similar to previous best
    ## Run 331 stress 0.08352634 
    ## ... Procrustes: rmse 2.718676e-06  max resid 4.888123e-06 
    ## ... Similar to previous best
    ## Run 332 stress 0.1134014 
    ## Run 333 stress 0.1279309 
    ## Run 334 stress 0.1279309 
    ## Run 335 stress 0.1134014 
    ## Run 336 stress 0.1134014 
    ## Run 337 stress 0.08352634 
    ## ... Procrustes: rmse 1.953794e-06  max resid 4.516288e-06 
    ## ... Similar to previous best
    ## Run 338 stress 0.1279309 
    ## Run 339 stress 0.1134014 
    ## Run 340 stress 0.1989154 
    ## Run 341 stress 0.1134015 
    ## Run 342 stress 0.08352634 
    ## ... Procrustes: rmse 2.67067e-06  max resid 5.901137e-06 
    ## ... Similar to previous best
    ## Run 343 stress 0.2031893 
    ## Run 344 stress 0.08352634 
    ## ... Procrustes: rmse 9.203037e-07  max resid 1.661687e-06 
    ## ... Similar to previous best
    ## Run 345 stress 0.08352634 
    ## ... Procrustes: rmse 9.633931e-07  max resid 1.597986e-06 
    ## ... Similar to previous best
    ## Run 346 stress 0.08352634 
    ## ... Procrustes: rmse 1.270712e-06  max resid 2.248937e-06 
    ## ... Similar to previous best
    ## Run 347 stress 0.08352634 
    ## ... Procrustes: rmse 2.243346e-06  max resid 5.135481e-06 
    ## ... Similar to previous best
    ## Run 348 stress 0.08352634 
    ## ... Procrustes: rmse 3.767637e-06  max resid 8.223299e-06 
    ## ... Similar to previous best
    ## Run 349 stress 0.08352634 
    ## ... Procrustes: rmse 2.114805e-06  max resid 4.246984e-06 
    ## ... Similar to previous best
    ## Run 350 stress 0.1989486 
    ## Run 351 stress 0.1279309 
    ## Run 352 stress 0.1279309 
    ## Run 353 stress 0.1134014 
    ## Run 354 stress 0.1279309 
    ## Run 355 stress 0.08352634 
    ## ... Procrustes: rmse 2.77017e-06  max resid 5.546697e-06 
    ## ... Similar to previous best
    ## Run 356 stress 0.08352634 
    ## ... Procrustes: rmse 1.740736e-06  max resid 3.845103e-06 
    ## ... Similar to previous best
    ## Run 357 stress 0.1989485 
    ## Run 358 stress 0.08352634 
    ## ... Procrustes: rmse 2.09452e-06  max resid 4.631893e-06 
    ## ... Similar to previous best
    ## Run 359 stress 0.1279309 
    ## Run 360 stress 0.1134015 
    ## Run 361 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 2.514803e-07  max resid 3.634426e-07 
    ## ... Similar to previous best
    ## Run 362 stress 0.1279309 
    ## Run 363 stress 0.3242672 
    ## Run 364 stress 0.08352634 
    ## ... Procrustes: rmse 1.758067e-06  max resid 3.86083e-06 
    ## ... Similar to previous best
    ## Run 365 stress 0.08352634 
    ## ... Procrustes: rmse 1.830859e-06  max resid 3.785023e-06 
    ## ... Similar to previous best
    ## Run 366 stress 0.1279309 
    ## Run 367 stress 0.1134014 
    ## Run 368 stress 0.08352634 
    ## ... Procrustes: rmse 3.739258e-06  max resid 8.391351e-06 
    ## ... Similar to previous best
    ## Run 369 stress 0.08352634 
    ## ... Procrustes: rmse 1.688207e-06  max resid 3.725501e-06 
    ## ... Similar to previous best
    ## Run 370 stress 0.08352634 
    ## ... Procrustes: rmse 3.212942e-06  max resid 6.959694e-06 
    ## ... Similar to previous best
    ## Run 371 stress 0.1134014 
    ## Run 372 stress 0.08352634 
    ## ... Procrustes: rmse 5.489764e-07  max resid 8.419208e-07 
    ## ... Similar to previous best
    ## Run 373 stress 0.08352634 
    ## ... Procrustes: rmse 2.152056e-06  max resid 4.615872e-06 
    ## ... Similar to previous best
    ## Run 374 stress 0.1134015 
    ## Run 375 stress 0.3311494 
    ## Run 376 stress 0.08352634 
    ## ... Procrustes: rmse 1.322812e-06  max resid 2.933045e-06 
    ## ... Similar to previous best
    ## Run 377 stress 0.08352634 
    ## ... Procrustes: rmse 3.46411e-06  max resid 7.85596e-06 
    ## ... Similar to previous best
    ## Run 378 stress 0.08352634 
    ## ... Procrustes: rmse 1.290979e-06  max resid 2.699634e-06 
    ## ... Similar to previous best
    ## Run 379 stress 0.08352634 
    ## ... Procrustes: rmse 5.538999e-07  max resid 1.212782e-06 
    ## ... Similar to previous best
    ## Run 380 stress 0.1134015 
    ## Run 381 stress 0.08352634 
    ## ... Procrustes: rmse 1.447011e-06  max resid 3.054819e-06 
    ## ... Similar to previous best
    ## Run 382 stress 0.08352634 
    ## ... Procrustes: rmse 1.936812e-06  max resid 3.26168e-06 
    ## ... Similar to previous best
    ## Run 383 stress 0.1279309 
    ## Run 384 stress 0.08352634 
    ## ... Procrustes: rmse 1.613453e-06  max resid 3.485004e-06 
    ## ... Similar to previous best
    ## Run 385 stress 0.08352634 
    ## ... Procrustes: rmse 1.75162e-06  max resid 3.936206e-06 
    ## ... Similar to previous best
    ## Run 386 stress 0.08352634 
    ## ... Procrustes: rmse 8.39782e-07  max resid 1.710619e-06 
    ## ... Similar to previous best
    ## Run 387 stress 0.08352634 
    ## ... Procrustes: rmse 2.512281e-06  max resid 5.53085e-06 
    ## ... Similar to previous best
    ## Run 388 stress 0.08352634 
    ## ... Procrustes: rmse 1.833003e-06  max resid 3.976252e-06 
    ## ... Similar to previous best
    ## Run 389 stress 0.1279309 
    ## Run 390 stress 0.1134015 
    ## Run 391 stress 0.08352634 
    ## ... Procrustes: rmse 1.627293e-06  max resid 3.514449e-06 
    ## ... Similar to previous best
    ## Run 392 stress 0.1134014 
    ## Run 393 stress 0.08352634 
    ## ... Procrustes: rmse 2.328152e-06  max resid 5.066184e-06 
    ## ... Similar to previous best
    ## Run 394 stress 0.08352634 
    ## ... Procrustes: rmse 1.901726e-06  max resid 4.148794e-06 
    ## ... Similar to previous best
    ## Run 395 stress 0.08352634 
    ## ... Procrustes: rmse 1.172828e-06  max resid 2.504391e-06 
    ## ... Similar to previous best
    ## Run 396 stress 0.08352634 
    ## ... Procrustes: rmse 8.187014e-07  max resid 1.924607e-06 
    ## ... Similar to previous best
    ## Run 397 stress 0.08352634 
    ## ... Procrustes: rmse 1.641288e-06  max resid 3.695294e-06 
    ## ... Similar to previous best
    ## Run 398 stress 0.1134014 
    ## Run 399 stress 0.1134015 
    ## Run 400 stress 0.08352634 
    ## ... Procrustes: rmse 3.533509e-06  max resid 7.860434e-06 
    ## ... Similar to previous best
    ## Run 401 stress 0.08352634 
    ## ... Procrustes: rmse 2.103267e-06  max resid 4.476697e-06 
    ## ... Similar to previous best
    ## Run 402 stress 0.1134014 
    ## Run 403 stress 0.08352634 
    ## ... Procrustes: rmse 2.686537e-06  max resid 6.013759e-06 
    ## ... Similar to previous best
    ## Run 404 stress 0.08352634 
    ## ... Procrustes: rmse 2.105685e-06  max resid 4.67498e-06 
    ## ... Similar to previous best
    ## Run 405 stress 0.08352634 
    ## ... Procrustes: rmse 2.112814e-06  max resid 4.501223e-06 
    ## ... Similar to previous best
    ## Run 406 stress 0.08352634 
    ## ... Procrustes: rmse 1.534274e-06  max resid 3.43038e-06 
    ## ... Similar to previous best
    ## Run 407 stress 0.1783151 
    ## Run 408 stress 0.08352634 
    ## ... Procrustes: rmse 1.89515e-06  max resid 3.828272e-06 
    ## ... Similar to previous best
    ## Run 409 stress 0.2031893 
    ## Run 410 stress 0.1134014 
    ## Run 411 stress 0.1134014 
    ## Run 412 stress 0.08352634 
    ## ... Procrustes: rmse 9.617469e-07  max resid 2.054744e-06 
    ## ... Similar to previous best
    ## Run 413 stress 0.08352634 
    ## ... Procrustes: rmse 1.090435e-06  max resid 2.358046e-06 
    ## ... Similar to previous best
    ## Run 414 stress 0.08352634 
    ## ... Procrustes: rmse 2.917865e-06  max resid 6.473309e-06 
    ## ... Similar to previous best
    ## Run 415 stress 0.08352634 
    ## ... Procrustes: rmse 2.840129e-06  max resid 6.202993e-06 
    ## ... Similar to previous best
    ## Run 416 stress 0.08352634 
    ## ... Procrustes: rmse 3.590781e-06  max resid 7.918084e-06 
    ## ... Similar to previous best
    ## Run 417 stress 0.08352634 
    ## ... Procrustes: rmse 3.725212e-06  max resid 8.231621e-06 
    ## ... Similar to previous best
    ## Run 418 stress 0.1279309 
    ## Run 419 stress 0.1134015 
    ## Run 420 stress 0.1279309 
    ## Run 421 stress 0.08352634 
    ## ... Procrustes: rmse 1.673496e-06  max resid 3.639935e-06 
    ## ... Similar to previous best
    ## Run 422 stress 0.1134014 
    ## Run 423 stress 0.08352634 
    ## ... Procrustes: rmse 1.625891e-06  max resid 3.228438e-06 
    ## ... Similar to previous best
    ## Run 424 stress 0.1134014 
    ## Run 425 stress 0.1134014 
    ## Run 426 stress 0.1134015 
    ## Run 427 stress 0.08352634 
    ## ... Procrustes: rmse 9.560209e-07  max resid 2.098377e-06 
    ## ... Similar to previous best
    ## Run 428 stress 0.1134014 
    ## Run 429 stress 0.1783092 
    ## Run 430 stress 0.08352634 
    ## ... Procrustes: rmse 1.995446e-06  max resid 4.339133e-06 
    ## ... Similar to previous best
    ## Run 431 stress 0.08352635 
    ## ... Procrustes: rmse 1.028564e-05  max resid 1.750873e-05 
    ## ... Similar to previous best
    ## Run 432 stress 0.1134014 
    ## Run 433 stress 0.1279309 
    ## Run 434 stress 0.1134014 
    ## Run 435 stress 0.1134015 
    ## Run 436 stress 0.08352634 
    ## ... Procrustes: rmse 1.273251e-06  max resid 2.820478e-06 
    ## ... Similar to previous best
    ## Run 437 stress 0.1279309 
    ## Run 438 stress 0.08352634 
    ## ... Procrustes: rmse 2.002837e-06  max resid 4.320941e-06 
    ## ... Similar to previous best
    ## Run 439 stress 0.1134015 
    ## Run 440 stress 0.08352634 
    ## ... Procrustes: rmse 2.591494e-07  max resid 4.411542e-07 
    ## ... Similar to previous best
    ## Run 441 stress 0.1279309 
    ## Run 442 stress 0.1279309 
    ## Run 443 stress 0.1279309 
    ## Run 444 stress 0.08352634 
    ## ... Procrustes: rmse 3.945727e-06  max resid 8.820108e-06 
    ## ... Similar to previous best
    ## Run 445 stress 0.08352634 
    ## ... Procrustes: rmse 2.890109e-06  max resid 6.447102e-06 
    ## ... Similar to previous best
    ## Run 446 stress 0.08352634 
    ## ... Procrustes: rmse 2.886566e-06  max resid 6.388395e-06 
    ## ... Similar to previous best
    ## Run 447 stress 0.3154666 
    ## Run 448 stress 0.08352634 
    ## ... Procrustes: rmse 2.193965e-06  max resid 4.902506e-06 
    ## ... Similar to previous best
    ## Run 449 stress 0.08352634 
    ## ... Procrustes: rmse 3.163159e-06  max resid 7.068606e-06 
    ## ... Similar to previous best
    ## Run 450 stress 0.1279309 
    ## Run 451 stress 0.08352634 
    ## ... Procrustes: rmse 2.456589e-06  max resid 5.319688e-06 
    ## ... Similar to previous best
    ## Run 452 stress 0.08352634 
    ## ... Procrustes: rmse 1.643881e-06  max resid 3.455697e-06 
    ## ... Similar to previous best
    ## Run 453 stress 0.08352634 
    ## ... Procrustes: rmse 1.284033e-06  max resid 2.860968e-06 
    ## ... Similar to previous best
    ## Run 454 stress 0.08352634 
    ## ... Procrustes: rmse 3.055445e-06  max resid 6.791187e-06 
    ## ... Similar to previous best
    ## Run 455 stress 0.08352634 
    ## ... Procrustes: rmse 1.246827e-05  max resid 2.783555e-05 
    ## ... Similar to previous best
    ## Run 456 stress 0.08352634 
    ## ... Procrustes: rmse 1.045874e-06  max resid 2.350963e-06 
    ## ... Similar to previous best
    ## Run 457 stress 0.08352634 
    ## ... Procrustes: rmse 1.16264e-05  max resid 2.376651e-05 
    ## ... Similar to previous best
    ## Run 458 stress 0.08352634 
    ## ... Procrustes: rmse 3.000648e-06  max resid 6.537552e-06 
    ## ... Similar to previous best
    ## Run 459 stress 0.1279309 
    ## Run 460 stress 0.08352634 
    ## ... Procrustes: rmse 4.717549e-06  max resid 1.047318e-05 
    ## ... Similar to previous best
    ## Run 461 stress 0.3255118 
    ## Run 462 stress 0.1134015 
    ## Run 463 stress 0.08352634 
    ## ... Procrustes: rmse 5.20766e-07  max resid 1.155052e-06 
    ## ... Similar to previous best
    ## Run 464 stress 0.08352634 
    ## ... Procrustes: rmse 1.496098e-06  max resid 2.896082e-06 
    ## ... Similar to previous best
    ## Run 465 stress 0.1134015 
    ## Run 466 stress 0.1279309 
    ## Run 467 stress 0.1279309 
    ## Run 468 stress 0.08352634 
    ## ... Procrustes: rmse 1.598742e-06  max resid 3.521245e-06 
    ## ... Similar to previous best
    ## Run 469 stress 0.08352634 
    ## ... Procrustes: rmse 3.849529e-06  max resid 8.561942e-06 
    ## ... Similar to previous best
    ## Run 470 stress 0.1134014 
    ## Run 471 stress 0.08352634 
    ## ... Procrustes: rmse 9.017917e-07  max resid 1.236932e-06 
    ## ... Similar to previous best
    ## Run 472 stress 0.1989156 
    ## Run 473 stress 0.1989488 
    ## Run 474 stress 0.127931 
    ## Run 475 stress 0.08352634 
    ## ... Procrustes: rmse 1.055757e-06  max resid 1.735233e-06 
    ## ... Similar to previous best
    ## Run 476 stress 0.1134014 
    ## Run 477 stress 0.08352634 
    ## ... Procrustes: rmse 2.566553e-06  max resid 5.650882e-06 
    ## ... Similar to previous best
    ## Run 478 stress 0.08352634 
    ## ... Procrustes: rmse 2.228649e-06  max resid 4.926571e-06 
    ## ... Similar to previous best
    ## Run 479 stress 0.1279309 
    ## Run 480 stress 0.08352634 
    ## ... Procrustes: rmse 5.15761e-07  max resid 1.04693e-06 
    ## ... Similar to previous best
    ## Run 481 stress 0.08352634 
    ## ... Procrustes: rmse 2.756026e-06  max resid 5.791896e-06 
    ## ... Similar to previous best
    ## Run 482 stress 0.08352634 
    ## ... Procrustes: rmse 1.489052e-06  max resid 3.236923e-06 
    ## ... Similar to previous best
    ## Run 483 stress 0.3070022 
    ## Run 484 stress 0.1134015 
    ## Run 485 stress 0.08352634 
    ## ... Procrustes: rmse 1.70046e-06  max resid 3.714225e-06 
    ## ... Similar to previous best
    ## Run 486 stress 0.1279309 
    ## Run 487 stress 0.1279309 
    ## Run 488 stress 0.08352634 
    ## ... Procrustes: rmse 5.875118e-06  max resid 1.323305e-05 
    ## ... Similar to previous best
    ## Run 489 stress 0.1134014 
    ## Run 490 stress 0.1134014 
    ## Run 491 stress 0.08352634 
    ## ... Procrustes: rmse 1.725911e-06  max resid 3.820091e-06 
    ## ... Similar to previous best
    ## Run 492 stress 0.1134014 
    ## Run 493 stress 0.1279309 
    ## Run 494 stress 0.1989159 
    ## Run 495 stress 0.1134015 
    ## Run 496 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 2.594219e-07  max resid 5.202307e-07 
    ## ... Similar to previous best
    ## Run 497 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 3.711991e-07  max resid 7.973374e-07 
    ## ... Similar to previous best
    ## Run 498 stress 0.08352634 
    ## ... Procrustes: rmse 1.088786e-06  max resid 1.849389e-06 
    ## ... Similar to previous best
    ## Run 499 stress 0.08352634 
    ## ... Procrustes: rmse 4.020671e-06  max resid 8.805073e-06 
    ## ... Similar to previous best
    ## Run 500 stress 0.1989486 
    ## *** Best solution repeated 3 times

``` r
# Stratified lakes and ocean sites
PD_beta_env_SO_NMDS <- metaMDS(PD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.92006e-05 
    ## Run 1 stress 0.0008142653 
    ## Run 2 stress 9.635479e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003312611  max resid 0.0008957562 
    ## ... Similar to previous best
    ## Run 3 stress 9.629919e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000247327  max resid 0.0005974202 
    ## ... Similar to previous best
    ## Run 4 stress 9.90572e-05 
    ## ... Procrustes: rmse 0.0006504584  max resid 0.001168456 
    ## ... Similar to previous best
    ## Run 5 stress 9.742021e-05 
    ## ... Procrustes: rmse 0.0003101707  max resid 0.0007147551 
    ## ... Similar to previous best
    ## Run 6 stress 0.0006689901 
    ## Run 7 stress 0.001255429 
    ## Run 8 stress 9.889276e-05 
    ## ... Procrustes: rmse 0.0003303004  max resid 0.0006878053 
    ## ... Similar to previous best
    ## Run 9 stress 0.0007286044 
    ## Run 10 stress 9.994567e-05 
    ## ... Procrustes: rmse 1.436502e-05  max resid 2.099087e-05 
    ## ... Similar to previous best
    ## Run 11 stress 9.992995e-05 
    ## ... Procrustes: rmse 0.0002764506  max resid 0.0007079411 
    ## ... Similar to previous best
    ## Run 12 stress 9.382032e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 1.612811e-05  max resid 3.007718e-05 
    ## ... Similar to previous best
    ## Run 13 stress 0.0006393204 
    ## Run 14 stress 9.972734e-05 
    ## ... Procrustes: rmse 2.655574e-05  max resid 4.249061e-05 
    ## ... Similar to previous best
    ## Run 15 stress 9.972239e-05 
    ## ... Procrustes: rmse 0.0002713065  max resid 0.0005335629 
    ## ... Similar to previous best
    ## Run 16 stress 9.95809e-05 
    ## ... Procrustes: rmse 0.0002702266  max resid 0.0007088964 
    ## ... Similar to previous best
    ## Run 17 stress 9.924086e-05 
    ## ... Procrustes: rmse 0.0003259819  max resid 0.000673713 
    ## ... Similar to previous best
    ## Run 18 stress 9.96429e-05 
    ## ... Procrustes: rmse 0.0002702118  max resid 0.0005334011 
    ## ... Similar to previous best
    ## Run 19 stress 0.0006886604 
    ## Run 20 stress 9.66119e-05 
    ## ... Procrustes: rmse 0.0002999871  max resid 0.000689254 
    ## ... Similar to previous best
    ## Run 21 stress 9.93007e-05 
    ## ... Procrustes: rmse 0.0002536191  max resid 0.0007041043 
    ## ... Similar to previous best
    ## Run 22 stress 0.001258951 
    ## Run 23 stress 9.907667e-05 
    ## ... Procrustes: rmse 0.0002708213  max resid 0.0006891694 
    ## ... Similar to previous best
    ## Run 24 stress 9.967341e-05 
    ## ... Procrustes: rmse 0.0004470991  max resid 0.000818589 
    ## ... Similar to previous best
    ## Run 25 stress 0.0007537276 
    ## Run 26 stress 9.850537e-05 
    ## ... Procrustes: rmse 0.0003245722  max resid 0.0006720004 
    ## ... Similar to previous best
    ## Run 27 stress 0.0008194493 
    ## Run 28 stress 9.842877e-05 
    ## ... Procrustes: rmse 0.0003067099  max resid 0.0007061766 
    ## ... Similar to previous best
    ## Run 29 stress 0.00118797 
    ## Run 30 stress 0.0007649 
    ## Run 31 stress 0.0009990451 
    ## Run 32 stress 9.989985e-05 
    ## ... Procrustes: rmse 0.0003101697  max resid 0.0007135246 
    ## ... Similar to previous best
    ## Run 33 stress 9.828499e-05 
    ## ... Procrustes: rmse 0.00025021  max resid 0.0006060109 
    ## ... Similar to previous best
    ## Run 34 stress 9.944812e-05 
    ## ... Procrustes: rmse 0.0002095269  max resid 0.000595143 
    ## ... Similar to previous best
    ## Run 35 stress 9.910493e-05 
    ## ... Procrustes: rmse 0.0003018697  max resid 0.0006922555 
    ## ... Similar to previous best
    ## Run 36 stress 9.888735e-05 
    ## ... Procrustes: rmse 0.0002688791  max resid 0.0005297902 
    ## ... Similar to previous best
    ## Run 37 stress 9.716791e-05 
    ## ... Procrustes: rmse 0.0003230666  max resid 0.000666947 
    ## ... Similar to previous best
    ## Run 38 stress 0.0005365073 
    ## ... Procrustes: rmse 0.004473103  max resid 0.007045857 
    ## ... Similar to previous best
    ## Run 39 stress 9.826228e-05 
    ## ... Procrustes: rmse 0.0003069773  max resid 0.0007059366 
    ## ... Similar to previous best
    ## Run 40 stress 0.0006408801 
    ## Run 41 stress 9.726882e-05 
    ## ... Procrustes: rmse 0.0004140441  max resid 0.0006393759 
    ## ... Similar to previous best
    ## Run 42 stress 0.0003451745 
    ## ... Procrustes: rmse 0.002671473  max resid 0.00409764 
    ## ... Similar to previous best
    ## Run 43 stress 9.999715e-05 
    ## ... Procrustes: rmse 0.0003286432  max resid 0.0006783333 
    ## ... Similar to previous best
    ## Run 44 stress 0.001253315 
    ## Run 45 stress 9.991933e-05 
    ## ... Procrustes: rmse 2.477152e-05  max resid 3.547883e-05 
    ## ... Similar to previous best
    ## Run 46 stress 0.0008639419 
    ## Run 47 stress 0.0008964583 
    ## Run 48 stress 0.0003732831 
    ## ... Procrustes: rmse 0.002920677  max resid 0.00447923 
    ## ... Similar to previous best
    ## Run 49 stress 0.001019774 
    ## Run 50 stress 0.001092885 
    ## Run 51 stress 9.788897e-05 
    ## ... Procrustes: rmse 0.0003062895  max resid 0.000703945 
    ## ... Similar to previous best
    ## Run 52 stress 0.001110009 
    ## Run 53 stress 0.0009068284 
    ## Run 54 stress 9.961369e-05 
    ## ... Procrustes: rmse 0.0003277944  max resid 0.0006768056 
    ## ... Similar to previous best
    ## Run 55 stress 0.001291488 
    ## Run 56 stress 0.001002377 
    ## Run 57 stress 9.819947e-05 
    ## ... Procrustes: rmse 0.0002496724  max resid 0.0006065788 
    ## ... Similar to previous best
    ## Run 58 stress 0.0008532952 
    ## Run 59 stress 0.001106835 
    ## Run 60 stress 0.0004453559 
    ## ... Procrustes: rmse 0.003662884  max resid 0.005787679 
    ## ... Similar to previous best
    ## Run 61 stress 9.910127e-05 
    ## ... Procrustes: rmse 0.0003267464  max resid 0.0006750098 
    ## ... Similar to previous best
    ## Run 62 stress 9.937873e-05 
    ## ... Procrustes: rmse 0.0003272844  max resid 0.0006768197 
    ## ... Similar to previous best
    ## Run 63 stress 9.70297e-05 
    ## ... Procrustes: rmse 0.0002483139  max resid 0.0006014749 
    ## ... Similar to previous best
    ## Run 64 stress 9.819882e-05 
    ## ... Procrustes: rmse 0.0003241351  max resid 0.0006712321 
    ## ... Similar to previous best
    ## Run 65 stress 9.988496e-05 
    ## ... Procrustes: rmse 2.91715e-05  max resid 4.065767e-05 
    ## ... Similar to previous best
    ## Run 66 stress 9.733717e-05 
    ## ... Procrustes: rmse 0.0003236193  max resid 0.0006710249 
    ## ... Similar to previous best
    ## Run 67 stress 9.934981e-05 
    ## ... Procrustes: rmse 0.000253161  max resid 0.000700025 
    ## ... Similar to previous best
    ## Run 68 stress 0.0007521014 
    ## Run 69 stress 9.612008e-05 
    ## ... Procrustes: rmse 0.0002592438  max resid 0.0006638699 
    ## ... Similar to previous best
    ## Run 70 stress 0.001251028 
    ## Run 71 stress 0.000787581 
    ## Run 72 stress 9.441366e-05 
    ## ... Procrustes: rmse 0.0002415912  max resid 0.0005824931 
    ## ... Similar to previous best
    ## Run 73 stress 0.0006555849 
    ## Run 74 stress 9.965051e-05 
    ## ... Procrustes: rmse 0.0003092818  max resid 0.00071054 
    ## ... Similar to previous best
    ## Run 75 stress 0.001074328 
    ## Run 76 stress 9.84089e-05 
    ## ... Procrustes: rmse 0.0003072443  max resid 0.0007065563 
    ## ... Similar to previous best
    ## Run 77 stress 0.3055763 
    ## Run 78 stress 9.887228e-05 
    ## ... Procrustes: rmse 0.0002668603  max resid 0.0005282598 
    ## ... Similar to previous best
    ## Run 79 stress 9.720087e-05 
    ## ... Procrustes: rmse 0.000304994  max resid 0.0007006701 
    ## ... Similar to previous best
    ## Run 80 stress 0.0007953169 
    ## Run 81 stress 0.3192024 
    ## Run 82 stress 0.0005817948 
    ## ... Procrustes: rmse 0.004740973  max resid 0.007290395 
    ## ... Similar to previous best
    ## Run 83 stress 0.0009688427 
    ## Run 84 stress 0.00109733 
    ## Run 85 stress 9.703634e-05 
    ## ... Procrustes: rmse 0.0002465856  max resid 0.0005979922 
    ## ... Similar to previous best
    ## Run 86 stress 9.854419e-05 
    ## ... Procrustes: rmse 0.0003075274  max resid 0.0007064976 
    ## ... Similar to previous best
    ## Run 87 stress 0.001088169 
    ## Run 88 stress 9.857482e-05 
    ## ... Procrustes: rmse 2.5365e-05  max resid 3.88134e-05 
    ## ... Similar to previous best
    ## Run 89 stress 9.882403e-05 
    ## ... Procrustes: rmse 0.000268867  max resid 0.0005288311 
    ## ... Similar to previous best
    ## Run 90 stress 0.0007131585 
    ## Run 91 stress 0.001136442 
    ## Run 92 stress 9.897102e-05 
    ## ... Procrustes: rmse 0.0003085811  max resid 0.0007100113 
    ## ... Similar to previous best
    ## Run 93 stress 0.0006222768 
    ## Run 94 stress 9.996151e-05 
    ## ... Procrustes: rmse 0.0002735024  max resid 0.0006948235 
    ## ... Similar to previous best
    ## Run 95 stress 9.953974e-05 
    ## ... Procrustes: rmse 0.0003061908  max resid 0.0007031295 
    ## ... Similar to previous best
    ## Run 96 stress 9.842781e-05 
    ## ... Procrustes: rmse 0.0002503914  max resid 0.0006063025 
    ## ... Similar to previous best
    ## Run 97 stress 9.914988e-05 
    ## ... Procrustes: rmse 0.0002512652  max resid 0.0006074808 
    ## ... Similar to previous best
    ## Run 98 stress 9.772563e-05 
    ## ... Procrustes: rmse 0.0002469477  max resid 0.0005983524 
    ## ... Similar to previous best
    ## Run 99 stress 0.0008617809 
    ## Run 100 stress 0.000865056 
    ## Run 101 stress 9.848782e-05 
    ## ... Procrustes: rmse 0.0002682099  max resid 0.0006822697 
    ## ... Similar to previous best
    ## Run 102 stress 9.512928e-05 
    ## ... Procrustes: rmse 0.000262096  max resid 0.0006722001 
    ## ... Similar to previous best
    ## Run 103 stress 0.0009065865 
    ## Run 104 stress 0.2570691 
    ## Run 105 stress 9.917584e-05 
    ## ... Procrustes: rmse 0.0002508922  max resid 0.0006060928 
    ## ... Similar to previous best
    ## Run 106 stress 0.001154337 
    ## Run 107 stress 0.0009090095 
    ## Run 108 stress 0.001256255 
    ## Run 109 stress 9.963718e-05 
    ## ... Procrustes: rmse 0.0002884437  max resid 0.000757666 
    ## ... Similar to previous best
    ## Run 110 stress 9.748728e-05 
    ## ... Procrustes: rmse 0.0003005347  max resid 0.0006906318 
    ## ... Similar to previous best
    ## Run 111 stress 9.874915e-05 
    ## ... Procrustes: rmse 0.0003258526  max resid 0.0006732392 
    ## ... Similar to previous best
    ## Run 112 stress 9.93359e-05 
    ## ... Procrustes: rmse 0.0002515693  max resid 0.0006105013 
    ## ... Similar to previous best
    ## Run 113 stress 9.87134e-05 
    ## ... Procrustes: rmse 0.0002677781  max resid 0.0005256584 
    ## ... Similar to previous best
    ## Run 114 stress 9.905258e-05 
    ## ... Procrustes: rmse 0.0002693814  max resid 0.0005292349 
    ## ... Similar to previous best
    ## Run 115 stress 0.0005760267 
    ## ... Procrustes: rmse 0.004686086  max resid 0.007199166 
    ## ... Similar to previous best
    ## Run 116 stress 9.978527e-05 
    ## ... Procrustes: rmse 3.026388e-05  max resid 4.577091e-05 
    ## ... Similar to previous best
    ## Run 117 stress 9.710326e-05 
    ## ... Procrustes: rmse 0.0002657883  max resid 0.000698649 
    ## ... Similar to previous best
    ## Run 118 stress 0.0004388668 
    ## ... Procrustes: rmse 0.003489649  max resid 0.005349042 
    ## ... Similar to previous best
    ## Run 119 stress 9.816473e-05 
    ## ... Procrustes: rmse 0.0002703463  max resid 0.0006884453 
    ## ... Similar to previous best
    ## Run 120 stress 9.904827e-05 
    ## ... Procrustes: rmse 0.0002542238  max resid 0.0006311803 
    ## ... Similar to previous best
    ## Run 121 stress 9.885494e-05 
    ## ... Procrustes: rmse 0.0003081619  max resid 0.0007086583 
    ## ... Similar to previous best
    ## Run 122 stress 9.907821e-05 
    ## ... Procrustes: rmse 0.0002085123  max resid 0.0005927398 
    ## ... Similar to previous best
    ## Run 123 stress 0.3407337 
    ## Run 124 stress 9.732169e-05 
    ## ... Procrustes: rmse 0.0003226  max resid 0.0006674066 
    ## ... Similar to previous best
    ## Run 125 stress 0.000757079 
    ## Run 126 stress 9.611574e-05 
    ## ... Procrustes: rmse 0.0002659871  max resid 0.0006793668 
    ## ... Similar to previous best
    ## Run 127 stress 9.893476e-05 
    ## ... Procrustes: rmse 0.000326524  max resid 0.0006744988 
    ## ... Similar to previous best
    ## Run 128 stress 0.0006905599 
    ## Run 129 stress 0.0009009793 
    ## Run 130 stress 9.891079e-05 
    ## ... Procrustes: rmse 0.0002506218  max resid 0.0006070803 
    ## ... Similar to previous best
    ## Run 131 stress 9.936999e-05 
    ## ... Procrustes: rmse 0.00030902  max resid 0.0007102702 
    ## ... Similar to previous best
    ## Run 132 stress 0.0004352755 
    ## ... Procrustes: rmse 0.003196408  max resid 0.004947044 
    ## ... Similar to previous best
    ## Run 133 stress 0.0005994889 
    ## Run 134 stress 0.000778985 
    ## Run 135 stress 9.862687e-05 
    ## ... Procrustes: rmse 0.0002078219  max resid 0.0005917854 
    ## ... Similar to previous best
    ## Run 136 stress 9.821753e-05 
    ## ... Procrustes: rmse 0.0003247758  max resid 0.00067241 
    ## ... Similar to previous best
    ## Run 137 stress 0.001314472 
    ## Run 138 stress 9.824411e-05 
    ## ... Procrustes: rmse 0.0002502353  max resid 0.0006063197 
    ## ... Similar to previous best
    ## Run 139 stress 0.0008175954 
    ## Run 140 stress 9.997645e-05 
    ## ... Procrustes: rmse 0.0003082548  max resid 0.0007084933 
    ## ... Similar to previous best
    ## Run 141 stress 9.625314e-05 
    ## ... Procrustes: rmse 0.0003192567  max resid 0.0006610065 
    ## ... Similar to previous best
    ## Run 142 stress 0.001038453 
    ## Run 143 stress 0.0008558314 
    ## Run 144 stress 0.001018419 
    ## Run 145 stress 0.001194024 
    ## Run 146 stress 0.0007961412 
    ## Run 147 stress 0.0008284238 
    ## Run 148 stress 0.246707 
    ## Run 149 stress 0.0007525372 
    ## Run 150 stress 0.0006324245 
    ## Run 151 stress 0.0008635677 
    ## Run 152 stress 9.713926e-05 
    ## ... Procrustes: rmse 0.0003013283  max resid 0.0006924311 
    ## ... Similar to previous best
    ## Run 153 stress 9.778173e-05 
    ## ... Procrustes: rmse 0.0002473668  max resid 0.0005993421 
    ## ... Similar to previous best
    ## Run 154 stress 9.935642e-05 
    ## ... Procrustes: rmse 0.0002697264  max resid 0.0005289471 
    ## ... Similar to previous best
    ## Run 155 stress 0.0007939054 
    ## Run 156 stress 8.743362e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004224829  max resid 0.0006671482 
    ## ... Similar to previous best
    ## Run 157 stress 9.922827e-05 
    ## ... Procrustes: rmse 0.0004477488  max resid 0.0007016603 
    ## ... Similar to previous best
    ## Run 158 stress 9.623314e-05 
    ## ... Procrustes: rmse 0.000405143  max resid 0.0005976648 
    ## ... Similar to previous best
    ## Run 159 stress 9.860911e-05 
    ## ... Procrustes: rmse 0.0004255082  max resid 0.0006207412 
    ## ... Similar to previous best
    ## Run 160 stress 9.972281e-05 
    ## ... Procrustes: rmse 0.0004515836  max resid 0.0007033527 
    ## ... Similar to previous best
    ## Run 161 stress 0.0007037912 
    ## Run 162 stress 9.584532e-05 
    ## ... Procrustes: rmse 0.000447411  max resid 0.0007217594 
    ## ... Similar to previous best
    ## Run 163 stress 9.597929e-05 
    ## ... Procrustes: rmse 0.0004305584  max resid 0.0006861838 
    ## ... Similar to previous best
    ## Run 164 stress 0.001102544 
    ## Run 165 stress 0.0005473491 
    ## ... Procrustes: rmse 0.00475195  max resid 0.00725116 
    ## ... Similar to previous best
    ## Run 166 stress 9.894202e-05 
    ## ... Procrustes: rmse 0.0004924841  max resid 0.0007765651 
    ## ... Similar to previous best
    ## Run 167 stress 0.00065428 
    ## Run 168 stress 0.0009494332 
    ## Run 169 stress 9.964907e-05 
    ## ... Procrustes: rmse 0.0004836706  max resid 0.00077571 
    ## ... Similar to previous best
    ## Run 170 stress 0.0003126776 
    ## ... Procrustes: rmse 0.002742333  max resid 0.00415647 
    ## ... Similar to previous best
    ## Run 171 stress 0.0006541843 
    ## Run 172 stress 9.965423e-05 
    ## ... Procrustes: rmse 0.0004452367  max resid 0.000690298 
    ## ... Similar to previous best
    ## Run 173 stress 0.0005454041 
    ## ... Procrustes: rmse 0.004791694  max resid 0.007284561 
    ## ... Similar to previous best
    ## Run 174 stress 9.506494e-05 
    ## ... Procrustes: rmse 0.0004135257  max resid 0.0006847985 
    ## ... Similar to previous best
    ## Run 175 stress 0.0007557775 
    ## Run 176 stress 0.0007264958 
    ## Run 177 stress 9.996767e-05 
    ## ... Procrustes: rmse 0.0004288709  max resid 0.0006834113 
    ## ... Similar to previous best
    ## Run 178 stress 9.909485e-05 
    ## ... Procrustes: rmse 0.0004409478  max resid 0.0007073997 
    ## ... Similar to previous best
    ## Run 179 stress 0.3275199 
    ## Run 180 stress 9.932748e-05 
    ## ... Procrustes: rmse 0.00044441  max resid 0.0007118594 
    ## ... Similar to previous best
    ## Run 181 stress 0.0008237142 
    ## Run 182 stress 0.00122066 
    ## Run 183 stress 9.75604e-05 
    ## ... Procrustes: rmse 0.0004382541  max resid 0.0007022846 
    ## ... Similar to previous best
    ## Run 184 stress 9.916138e-05 
    ## ... Procrustes: rmse 0.0004471892  max resid 0.0007166692 
    ## ... Similar to previous best
    ## Run 185 stress 9.730438e-05 
    ## ... Procrustes: rmse 0.0004412838  max resid 0.0006896976 
    ## ... Similar to previous best
    ## Run 186 stress 9.978718e-05 
    ## ... Procrustes: rmse 0.0004486307  max resid 0.0007085888 
    ## ... Similar to previous best
    ## Run 187 stress 9.743097e-05 
    ## ... Procrustes: rmse 0.0004347951  max resid 0.0006742301 
    ## ... Similar to previous best
    ## Run 188 stress 9.397159e-05 
    ## ... Procrustes: rmse 0.0004206178  max resid 0.0006527539 
    ## ... Similar to previous best
    ## Run 189 stress 9.576555e-05 
    ## ... Procrustes: rmse 0.0003987  max resid 0.00057875 
    ## ... Similar to previous best
    ## Run 190 stress 9.8494e-05 
    ## ... Procrustes: rmse 0.0004252922  max resid 0.0006181621 
    ## ... Similar to previous best
    ## Run 191 stress 0.0006935427 
    ## Run 192 stress 0.001257607 
    ## Run 193 stress 0.0007539957 
    ## Run 194 stress 9.720102e-05 
    ## ... Procrustes: rmse 0.000451021  max resid 0.0007290405 
    ## ... Similar to previous best
    ## Run 195 stress 9.92617e-05 
    ## ... Procrustes: rmse 0.0004940986  max resid 0.0007784566 
    ## ... Similar to previous best
    ## Run 196 stress 0.0005932846 
    ## Run 197 stress 9.962599e-05 
    ## ... Procrustes: rmse 0.0004512535  max resid 0.0007035996 
    ## ... Similar to previous best
    ## Run 198 stress 9.929874e-05 
    ## ... Procrustes: rmse 0.0004215805  max resid 0.0006996998 
    ## ... Similar to previous best
    ## Run 199 stress 9.939426e-05 
    ## ... Procrustes: rmse 0.0004389987  max resid 0.0006871477 
    ## ... Similar to previous best
    ## Run 200 stress 0.001103336 
    ## Run 201 stress 0.000802315 
    ## Run 202 stress 9.898822e-05 
    ## ... Procrustes: rmse 0.0004263004  max resid 0.0006178472 
    ## ... Similar to previous best
    ## Run 203 stress 9.716518e-05 
    ## ... Procrustes: rmse 0.0004253508  max resid 0.000702799 
    ## ... Similar to previous best
    ## Run 204 stress 0.0005091537 
    ## ... Procrustes: rmse 0.004464741  max resid 0.006782962 
    ## ... Similar to previous best
    ## Run 205 stress 0.0007834355 
    ## Run 206 stress 0.0006853671 
    ## Run 207 stress 9.812029e-05 
    ## ... Procrustes: rmse 0.0004415064  max resid 0.0007089704 
    ## ... Similar to previous best
    ## Run 208 stress 9.858209e-05 
    ## ... Procrustes: rmse 0.0004311008  max resid 0.0007126884 
    ## ... Similar to previous best
    ## Run 209 stress 9.602411e-05 
    ## ... Procrustes: rmse 0.0004654104  max resid 0.0007476909 
    ## ... Similar to previous best
    ## Run 210 stress 0.0002629022 
    ## ... Procrustes: rmse 0.002302871  max resid 0.003479358 
    ## ... Similar to previous best
    ## Run 211 stress 9.960853e-05 
    ## ... Procrustes: rmse 0.0004606863  max resid 0.0007406859 
    ## ... Similar to previous best
    ## Run 212 stress 0.0007584937 
    ## Run 213 stress 0.0005259928 
    ## ... Procrustes: rmse 0.004630462  max resid 0.007075822 
    ## ... Similar to previous best
    ## Run 214 stress 0.0005818366 
    ## ... Procrustes: rmse 0.005108859  max resid 0.007778061 
    ## Run 215 stress 0.0006497719 
    ## Run 216 stress 9.872846e-05 
    ## ... Procrustes: rmse 0.0004172356  max resid 0.0006148609 
    ## ... Similar to previous best
    ## Run 217 stress 9.872547e-05 
    ## ... Procrustes: rmse 0.0004191642  max resid 0.0006049047 
    ## ... Similar to previous best
    ## Run 218 stress 0.0006411229 
    ## Run 219 stress 9.872015e-05 
    ## ... Procrustes: rmse 0.0004601034  max resid 0.0007398911 
    ## ... Similar to previous best
    ## Run 220 stress 9.750474e-05 
    ## ... Procrustes: rmse 0.0004403626  max resid 0.0007057064 
    ## ... Similar to previous best
    ## Run 221 stress 9.969729e-05 
    ## ... Procrustes: rmse 0.0004271293  max resid 0.0006266184 
    ## ... Similar to previous best
    ## Run 222 stress 9.952825e-05 
    ## ... Procrustes: rmse 0.000448642  max resid 0.0007176648 
    ## ... Similar to previous best
    ## Run 223 stress 9.878026e-05 
    ## ... Procrustes: rmse 0.0004607863  max resid 0.0007411798 
    ## ... Similar to previous best
    ## Run 224 stress 9.988219e-05 
    ## ... Procrustes: rmse 0.0004364512  max resid 0.0007198593 
    ## ... Similar to previous best
    ## Run 225 stress 0.0003101696 
    ## ... Procrustes: rmse 0.002734649  max resid 0.004106258 
    ## ... Similar to previous best
    ## Run 226 stress 9.382324e-05 
    ## ... Procrustes: rmse 0.0004241334  max resid 0.0006843741 
    ## ... Similar to previous best
    ## Run 227 stress 0.001178782 
    ## Run 228 stress 0.0007447762 
    ## Run 229 stress 9.687635e-05 
    ## ... Procrustes: rmse 0.0004375649  max resid 0.0006999789 
    ## ... Similar to previous best
    ## Run 230 stress 0.0003161958 
    ## ... Procrustes: rmse 0.00277049  max resid 0.004208921 
    ## ... Similar to previous best
    ## Run 231 stress 9.930666e-05 
    ## ... Procrustes: rmse 0.0004279471  max resid 0.0006193784 
    ## ... Similar to previous best
    ## Run 232 stress 9.615632e-05 
    ## ... Procrustes: rmse 0.0004685036  max resid 0.0007531686 
    ## ... Similar to previous best
    ## Run 233 stress 0.0003602178 
    ## ... Procrustes: rmse 0.003170936  max resid 0.00477283 
    ## ... Similar to previous best
    ## Run 234 stress 9.91903e-05 
    ## ... Procrustes: rmse 0.0004335549  max resid 0.0007151702 
    ## ... Similar to previous best
    ## Run 235 stress 0.0007458253 
    ## Run 236 stress 0.0005855567 
    ## ... Procrustes: rmse 0.005109596  max resid 0.007753787 
    ## Run 237 stress 0.001221022 
    ## Run 238 stress 9.900425e-05 
    ## ... Procrustes: rmse 0.0004274823  max resid 0.0006219981 
    ## ... Similar to previous best
    ## Run 239 stress 9.850865e-05 
    ## ... Procrustes: rmse 0.0004756401  max resid 0.0007679667 
    ## ... Similar to previous best
    ## Run 240 stress 0.0006214969 
    ## Run 241 stress 9.930752e-05 
    ## ... Procrustes: rmse 0.0004941979  max resid 0.0007799509 
    ## ... Similar to previous best
    ## Run 242 stress 0.001097404 
    ## Run 243 stress 9.915749e-05 
    ## ... Procrustes: rmse 0.0004334576  max resid 0.0007156344 
    ## ... Similar to previous best
    ## Run 244 stress 0.001176628 
    ## Run 245 stress 9.673799e-05 
    ## ... Procrustes: rmse 0.000417646  max resid 0.0006042588 
    ## ... Similar to previous best
    ## Run 246 stress 9.874604e-05 
    ## ... Procrustes: rmse 0.0004278974  max resid 0.0007053109 
    ## ... Similar to previous best
    ## Run 247 stress 9.721827e-05 
    ## ... Procrustes: rmse 0.0004384141  max resid 0.0007064548 
    ## ... Similar to previous best
    ## Run 248 stress 9.979583e-05 
    ## ... Procrustes: rmse 0.0004291159  max resid 0.0006220512 
    ## ... Similar to previous best
    ## Run 249 stress 9.995482e-05 
    ## ... Procrustes: rmse 0.0004511968  max resid 0.0007073874 
    ## ... Similar to previous best
    ## Run 250 stress 9.924846e-05 
    ## ... Procrustes: rmse 0.0004816056  max resid 0.0007715568 
    ## ... Similar to previous best
    ## Run 251 stress 9.220742e-05 
    ## ... Procrustes: rmse 0.0001520315  max resid 0.0002925147 
    ## ... Similar to previous best
    ## Run 252 stress 0.0007561245 
    ## Run 253 stress 9.912201e-05 
    ## ... Procrustes: rmse 0.0004931343  max resid 0.0007747691 
    ## ... Similar to previous best
    ## Run 254 stress 0.0006324984 
    ## Run 255 stress 9.851407e-05 
    ## ... Procrustes: rmse 0.0004310808  max resid 0.0007120521 
    ## ... Similar to previous best
    ## Run 256 stress 9.832065e-05 
    ## ... Procrustes: rmse 0.000458563  max resid 0.0007397258 
    ## ... Similar to previous best
    ## Run 257 stress 9.416096e-05 
    ## ... Procrustes: rmse 0.000421596  max resid 0.0006607375 
    ## ... Similar to previous best
    ## Run 258 stress 9.612062e-05 
    ## ... Procrustes: rmse 0.0004195801  max resid 0.0006966763 
    ## ... Similar to previous best
    ## Run 259 stress 0.0008564292 
    ## Run 260 stress 9.926772e-05 
    ## ... Procrustes: rmse 0.0004272384  max resid 0.0006272551 
    ## ... Similar to previous best
    ## Run 261 stress 0.0007320361 
    ## Run 262 stress 9.730228e-05 
    ## ... Procrustes: rmse 0.0004372496  max resid 0.0007010029 
    ## ... Similar to previous best
    ## Run 263 stress 9.7847e-05 
    ## ... Procrustes: rmse 0.0004549207  max resid 0.0007302164 
    ## ... Similar to previous best
    ## Run 264 stress 9.989817e-05 
    ## ... Procrustes: rmse 0.00048159  max resid 0.0007787017 
    ## ... Similar to previous best
    ## Run 265 stress 9.951451e-05 
    ## ... Procrustes: rmse 0.0004504186  max resid 0.0007032874 
    ## ... Similar to previous best
    ## Run 266 stress 0.001193824 
    ## Run 267 stress 9.984989e-05 
    ## ... Procrustes: rmse 0.0004288052  max resid 0.0007059184 
    ## ... Similar to previous best
    ## Run 268 stress 0.3407337 
    ## Run 269 stress 0.001211179 
    ## Run 270 stress 9.996042e-05 
    ## ... Procrustes: rmse 0.0004526136  max resid 0.0007050208 
    ## ... Similar to previous best
    ## Run 271 stress 9.837188e-05 
    ## ... Procrustes: rmse 0.0004232539  max resid 0.0006147974 
    ## ... Similar to previous best
    ## Run 272 stress 9.879688e-05 
    ## ... Procrustes: rmse 0.0004039006  max resid 0.0006446053 
    ## ... Similar to previous best
    ## Run 273 stress 0.001052262 
    ## Run 274 stress 9.731879e-05 
    ## ... Procrustes: rmse 0.0004109299  max resid 0.0006774104 
    ## ... Similar to previous best
    ## Run 275 stress 0.0009973986 
    ## Run 276 stress 0.0003884866 
    ## ... Procrustes: rmse 0.003412419  max resid 0.005199222 
    ## ... Similar to previous best
    ## Run 277 stress 9.965507e-05 
    ## ... Procrustes: rmse 0.0004491995  max resid 0.000702077 
    ## ... Similar to previous best
    ## Run 278 stress 9.938949e-05 
    ## ... Procrustes: rmse 0.0004288875  max resid 0.0006245109 
    ## ... Similar to previous best
    ## Run 279 stress 9.966582e-05 
    ## ... Procrustes: rmse 0.0004482947  max resid 0.0007181314 
    ## ... Similar to previous best
    ## Run 280 stress 0.0008158468 
    ## Run 281 stress 9.555214e-05 
    ## ... Procrustes: rmse 0.0004360592  max resid 0.0007043554 
    ## ... Similar to previous best
    ## Run 282 stress 9.976224e-05 
    ## ... Procrustes: rmse 0.0004518376  max resid 0.0007042058 
    ## ... Similar to previous best
    ## Run 283 stress 0.0009930567 
    ## Run 284 stress 9.582452e-05 
    ## ... Procrustes: rmse 0.0004668564  max resid 0.0007508841 
    ## ... Similar to previous best
    ## Run 285 stress 9.878723e-05 
    ## ... Procrustes: rmse 0.0004472524  max resid 0.0006958576 
    ## ... Similar to previous best
    ## Run 286 stress 0.0008487014 
    ## Run 287 stress 9.892019e-05 
    ## ... Procrustes: rmse 0.000445005  max resid 0.000702791 
    ## ... Similar to previous best
    ## Run 288 stress 0.001121184 
    ## Run 289 stress 0.3131342 
    ## Run 290 stress 8.879429e-05 
    ## ... Procrustes: rmse 7.901742e-05  max resid 0.0001903838 
    ## ... Similar to previous best
    ## Run 291 stress 0.001100673 
    ## Run 292 stress 9.949593e-05 
    ## ... Procrustes: rmse 0.0004205286  max resid 0.0006921405 
    ## ... Similar to previous best
    ## Run 293 stress 9.964913e-05 
    ## ... Procrustes: rmse 0.0004949855  max resid 0.0007745487 
    ## ... Similar to previous best
    ## Run 294 stress 0.0006052659 
    ## Run 295 stress 0.001096392 
    ## Run 296 stress 0.000804559 
    ## Run 297 stress 9.979757e-05 
    ## ... Procrustes: rmse 0.0004498311  max resid 0.0007052542 
    ## ... Similar to previous best
    ## Run 298 stress 0.001115111 
    ## Run 299 stress 0.0003546381 
    ## ... Procrustes: rmse 0.003105748  max resid 0.004730236 
    ## ... Similar to previous best
    ## Run 300 stress 9.435181e-05 
    ## ... Procrustes: rmse 0.0004051059  max resid 0.0005981379 
    ## ... Similar to previous best
    ## Run 301 stress 0.001116681 
    ## Run 302 stress 0.001157466 
    ## Run 303 stress 9.963491e-05 
    ## ... Procrustes: rmse 0.0004931206  max resid 0.0007832008 
    ## ... Similar to previous best
    ## Run 304 stress 9.938886e-05 
    ## ... Procrustes: rmse 0.0004478197  max resid 0.0006986085 
    ## ... Similar to previous best
    ## Run 305 stress 9.933756e-05 
    ## ... Procrustes: rmse 0.000478564  max resid 0.0007741735 
    ## ... Similar to previous best
    ## Run 306 stress 0.00118757 
    ## Run 307 stress 9.672888e-05 
    ## ... Procrustes: rmse 0.0004129516  max resid 0.0006040883 
    ## ... Similar to previous best
    ## Run 308 stress 9.947932e-05 
    ## ... Procrustes: rmse 0.000448785  max resid 0.0007029888 
    ## ... Similar to previous best
    ## Run 309 stress 9.837479e-05 
    ## ... Procrustes: rmse 0.0004442796  max resid 0.0006906702 
    ## ... Similar to previous best
    ## Run 310 stress 9.770363e-05 
    ## ... Procrustes: rmse 0.0004716757  max resid 0.000761351 
    ## ... Similar to previous best
    ## Run 311 stress 0.0007481832 
    ## Run 312 stress 9.767029e-05 
    ## ... Procrustes: rmse 0.0004269012  max resid 0.0007028051 
    ## ... Similar to previous best
    ## Run 313 stress 9.856966e-05 
    ## ... Procrustes: rmse 0.0004598594  max resid 0.0007402893 
    ## ... Similar to previous best
    ## Run 314 stress 0.000652917 
    ## Run 315 stress 0.001240492 
    ## Run 316 stress 9.820797e-05 
    ## ... Procrustes: rmse 0.0004292447  max resid 0.0007097268 
    ## ... Similar to previous best
    ## Run 317 stress 0.000626689 
    ## Run 318 stress 9.983848e-05 
    ## ... Procrustes: rmse 0.0004811214  max resid 0.0007735775 
    ## ... Similar to previous best
    ## Run 319 stress 9.868754e-05 
    ## ... Procrustes: rmse 0.0004591658  max resid 0.0007382877 
    ## ... Similar to previous best
    ## Run 320 stress 9.99768e-05 
    ## ... Procrustes: rmse 0.0004959277  max resid 0.0007792053 
    ## ... Similar to previous best
    ## Run 321 stress 9.86318e-05 
    ## ... Procrustes: rmse 0.0004449578  max resid 0.0007123245 
    ## ... Similar to previous best
    ## Run 322 stress 9.827093e-05 
    ## ... Procrustes: rmse 0.0004410736  max resid 0.0007037064 
    ## ... Similar to previous best
    ## Run 323 stress 0.001040763 
    ## Run 324 stress 9.739796e-05 
    ## ... Procrustes: rmse 0.0004543338  max resid 0.0007316165 
    ## ... Similar to previous best
    ## Run 325 stress 9.977542e-05 
    ## ... Procrustes: rmse 0.0004637989  max resid 0.00076098 
    ## ... Similar to previous best
    ## Run 326 stress 9.966919e-05 
    ## ... Procrustes: rmse 0.0004303146  max resid 0.0006252571 
    ## ... Similar to previous best
    ## Run 327 stress 9.770727e-05 
    ## ... Procrustes: rmse 0.0004358855  max resid 0.0007017773 
    ## ... Similar to previous best
    ## Run 328 stress 0.001252517 
    ## Run 329 stress 9.917862e-05 
    ## ... Procrustes: rmse 0.0004793535  max resid 0.0007745236 
    ## ... Similar to previous best
    ## Run 330 stress 0.0005688528 
    ## ... Procrustes: rmse 0.00499617  max resid 0.007601396 
    ## ... Similar to previous best
    ## Run 331 stress 9.209647e-05 
    ## ... Procrustes: rmse 0.0002292857  max resid 0.0004317931 
    ## ... Similar to previous best
    ## Run 332 stress 0.0005649078 
    ## ... Procrustes: rmse 0.004960976  max resid 0.007546386 
    ## ... Similar to previous best
    ## Run 333 stress 9.926091e-05 
    ## ... Procrustes: rmse 0.0004470159  max resid 0.0007156399 
    ## ... Similar to previous best
    ## Run 334 stress 0.0009927822 
    ## Run 335 stress 9.850222e-05 
    ## ... Procrustes: rmse 0.0004223192  max resid 0.000612803 
    ## ... Similar to previous best
    ## Run 336 stress 0.0007905807 
    ## Run 337 stress 9.9468e-05 
    ## ... Procrustes: rmse 0.0005075991  max resid 0.0007686622 
    ## ... Similar to previous best
    ## Run 338 stress 0.001320287 
    ## Run 339 stress 0.001282267 
    ## Run 340 stress 0.001106063 
    ## Run 341 stress 9.859814e-05 
    ## ... Procrustes: rmse 0.0004456661  max resid 0.0006982071 
    ## ... Similar to previous best
    ## Run 342 stress 0.3407335 
    ## Run 343 stress 9.639243e-05 
    ## ... Procrustes: rmse 0.0004636446  max resid 0.0007497968 
    ## ... Similar to previous best
    ## Run 344 stress 9.805767e-05 
    ## ... Procrustes: rmse 0.0004571952  max resid 0.0007363939 
    ## ... Similar to previous best
    ## Run 345 stress 9.94668e-05 
    ## ... Procrustes: rmse 0.0004488746  max resid 0.0007021802 
    ## ... Similar to previous best
    ## Run 346 stress 9.920783e-05 
    ## ... Procrustes: rmse 0.0004468419  max resid 0.0006968852 
    ## ... Similar to previous best
    ## Run 347 stress 9.799987e-05 
    ## ... Procrustes: rmse 0.0004409903  max resid 0.0007052961 
    ## ... Similar to previous best
    ## Run 348 stress 9.892277e-05 
    ## ... Procrustes: rmse 0.0004325113  max resid 0.0007146254 
    ## ... Similar to previous best
    ## Run 349 stress 0.000750592 
    ## Run 350 stress 9.934164e-05 
    ## ... Procrustes: rmse 0.0005050225  max resid 0.0007643117 
    ## ... Similar to previous best
    ## Run 351 stress 9.95503e-05 
    ## ... Procrustes: rmse 0.0004957664  max resid 0.0007801847 
    ## ... Similar to previous best
    ## Run 352 stress 9.924368e-05 
    ## ... Procrustes: rmse 0.0004622408  max resid 0.0007439259 
    ## ... Similar to previous best
    ## Run 353 stress 9.9788e-05 
    ## ... Procrustes: rmse 0.0004473749  max resid 0.0007005557 
    ## ... Similar to previous best
    ## Run 354 stress 9.34395e-05 
    ## ... Procrustes: rmse 0.0004100445  max resid 0.0006784225 
    ## ... Similar to previous best
    ## Run 355 stress 9.971952e-05 
    ## ... Procrustes: rmse 0.0004398859  max resid 0.0006845337 
    ## ... Similar to previous best
    ## Run 356 stress 9.865641e-05 
    ## ... Procrustes: rmse 0.0004738215  max resid 0.0007652213 
    ## ... Similar to previous best
    ## Run 357 stress 9.979847e-05 
    ## ... Procrustes: rmse 0.0004830274  max resid 0.0007768855 
    ## ... Similar to previous best
    ## Run 358 stress 0.0006774983 
    ## Run 359 stress 0.0008564095 
    ## Run 360 stress 9.883901e-05 
    ## ... Procrustes: rmse 0.0004264302  max resid 0.0006189225 
    ## ... Similar to previous best
    ## Run 361 stress 9.94598e-05 
    ## ... Procrustes: rmse 0.0004437417  max resid 0.0007152647 
    ## ... Similar to previous best
    ## Run 362 stress 9.61323e-05 
    ## ... Procrustes: rmse 0.0004209762  max resid 0.0006973082 
    ## ... Similar to previous best
    ## Run 363 stress 9.693762e-05 
    ## ... Procrustes: rmse 0.0004171635  max resid 0.0006085673 
    ## ... Similar to previous best
    ## Run 364 stress 9.853306e-05 
    ## ... Procrustes: rmse 0.0004465828  max resid 0.0006961857 
    ## ... Similar to previous best
    ## Run 365 stress 9.985265e-05 
    ## ... Procrustes: rmse 0.0005098894  max resid 0.0007717799 
    ## ... Similar to previous best
    ## Run 366 stress 9.83111e-05 
    ## ... Procrustes: rmse 0.0004434765  max resid 0.0007095324 
    ## ... Similar to previous best
    ## Run 367 stress 0.0007450732 
    ## Run 368 stress 0.0005382542 
    ## ... Procrustes: rmse 0.004666157  max resid 0.007119412 
    ## ... Similar to previous best
    ## Run 369 stress 9.8331e-05 
    ## ... Procrustes: rmse 0.0004428538  max resid 0.0006880956 
    ## ... Similar to previous best
    ## Run 370 stress 9.883029e-05 
    ## ... Procrustes: rmse 0.0004455261  max resid 0.0007003076 
    ## ... Similar to previous best
    ## Run 371 stress 9.849874e-05 
    ## ... Procrustes: rmse 0.0004448018  max resid 0.0006982291 
    ## ... Similar to previous best
    ## Run 372 stress 0.0005291875 
    ## ... Procrustes: rmse 0.004649852  max resid 0.007065899 
    ## ... Similar to previous best
    ## Run 373 stress 0.001114422 
    ## Run 374 stress 0.0006312566 
    ## Run 375 stress 9.975968e-05 
    ## ... Procrustes: rmse 0.0008677915  max resid 0.001233086 
    ## ... Similar to previous best
    ## Run 376 stress 0.0006825982 
    ## Run 377 stress 0.001268883 
    ## Run 378 stress 9.972456e-05 
    ## ... Procrustes: rmse 0.0004287201  max resid 0.000627205 
    ## ... Similar to previous best
    ## Run 379 stress 9.985801e-05 
    ## ... Procrustes: rmse 0.0004521991  max resid 0.0007061082 
    ## ... Similar to previous best
    ## Run 380 stress 9.880575e-05 
    ## ... Procrustes: rmse 0.0004448034  max resid 0.0007111117 
    ## ... Similar to previous best
    ## Run 381 stress 9.685868e-05 
    ## ... Procrustes: rmse 0.0004520357  max resid 0.0007310612 
    ## ... Similar to previous best
    ## Run 382 stress 0.000721053 
    ## Run 383 stress 9.847106e-05 
    ## ... Procrustes: rmse 0.0004431369  max resid 0.000711062 
    ## ... Similar to previous best
    ## Run 384 stress 9.89138e-05 
    ## ... Procrustes: rmse 0.0004807402  max resid 0.0007727097 
    ## ... Similar to previous best
    ## Run 385 stress 0.001058526 
    ## Run 386 stress 0.0006496577 
    ## Run 387 stress 9.838778e-05 
    ## ... Procrustes: rmse 0.000427554  max resid 0.0007092617 
    ## ... Similar to previous best
    ## Run 388 stress 9.715422e-05 
    ## ... Procrustes: rmse 0.0004386687  max resid 0.0007045761 
    ## ... Similar to previous best
    ## Run 389 stress 9.953086e-05 
    ## ... Procrustes: rmse 0.000434383  max resid 0.0007161311 
    ## ... Similar to previous best
    ## Run 390 stress 9.969791e-05 
    ## ... Procrustes: rmse 0.000423433  max resid 0.0006110976 
    ## ... Similar to previous best
    ## Run 391 stress 9.877576e-05 
    ## ... Procrustes: rmse 0.000491432  max resid 0.0007750058 
    ## ... Similar to previous best
    ## Run 392 stress 9.806837e-05 
    ## ... Procrustes: rmse 0.0004289344  max resid 0.0007074024 
    ## ... Similar to previous best
    ## Run 393 stress 9.911322e-05 
    ## ... Procrustes: rmse 0.0004809865  max resid 0.0007727303 
    ## ... Similar to previous best
    ## Run 394 stress 9.898243e-05 
    ## ... Procrustes: rmse 0.00040481  max resid 0.0006778426 
    ## ... Similar to previous best
    ## Run 395 stress 0.0002849015 
    ## ... Procrustes: rmse 0.002495495  max resid 0.003780601 
    ## ... Similar to previous best
    ## Run 396 stress 0.0008493264 
    ## Run 397 stress 9.982078e-05 
    ## ... Procrustes: rmse 0.0004520462  max resid 0.0007051622 
    ## ... Similar to previous best
    ## Run 398 stress 9.983714e-05 
    ## ... Procrustes: rmse 0.0004290102  max resid 0.0006212367 
    ## ... Similar to previous best
    ## Run 399 stress 9.590505e-05 
    ## ... Procrustes: rmse 0.0004213225  max resid 0.0006668678 
    ## ... Similar to previous best
    ## Run 400 stress 9.992803e-05 
    ## ... Procrustes: rmse 0.0004494786  max resid 0.0007059379 
    ## ... Similar to previous best
    ## Run 401 stress 0.0001676118 
    ## ... Procrustes: rmse 0.001484927  max resid 0.002178292 
    ## ... Similar to previous best
    ## Run 402 stress 0.0007967549 
    ## Run 403 stress 9.644816e-05 
    ## ... Procrustes: rmse 0.0004223845  max resid 0.0006987828 
    ## ... Similar to previous best
    ## Run 404 stress 9.979167e-05 
    ## ... Procrustes: rmse 0.0004470835  max resid 0.0007192601 
    ## ... Similar to previous best
    ## Run 405 stress 9.782921e-05 
    ## ... Procrustes: rmse 0.0004423269  max resid 0.0006879587 
    ## ... Similar to previous best
    ## Run 406 stress 9.981763e-05 
    ## ... Procrustes: rmse 0.0004298012  max resid 0.0006289206 
    ## ... Similar to previous best
    ## Run 407 stress 0.0006187348 
    ## Run 408 stress 0.0006459967 
    ## Run 409 stress 9.835151e-05 
    ## ... Procrustes: rmse 0.0004147805  max resid 0.000605719 
    ## ... Similar to previous best
    ## Run 410 stress 0.001267955 
    ## Run 411 stress 0.001252718 
    ## Run 412 stress 9.880824e-05 
    ## ... Procrustes: rmse 0.0004254421  max resid 0.000621031 
    ## ... Similar to previous best
    ## Run 413 stress 9.870824e-05 
    ## ... Procrustes: rmse 0.000431759  max resid 0.0007126452 
    ## ... Similar to previous best
    ## Run 414 stress 0.001282537 
    ## Run 415 stress 9.767063e-05 
    ## ... Procrustes: rmse 0.0004264215  max resid 0.0007058424 
    ## ... Similar to previous best
    ## Run 416 stress 0.0006772882 
    ## Run 417 stress 9.996926e-05 
    ## ... Procrustes: rmse 0.0004811255  max resid 0.0007772892 
    ## ... Similar to previous best
    ## Run 418 stress 9.233538e-05 
    ## ... Procrustes: rmse 8.264673e-05  max resid 0.0002070578 
    ## ... Similar to previous best
    ## Run 419 stress 0.001170613 
    ## Run 420 stress 9.935023e-05 
    ## ... Procrustes: rmse 0.0004583199  max resid 0.0007363155 
    ## ... Similar to previous best
    ## Run 421 stress 9.711237e-05 
    ## ... Procrustes: rmse 0.0003690189  max resid 0.000648234 
    ## ... Similar to previous best
    ## Run 422 stress 9.671598e-05 
    ## ... Procrustes: rmse 0.0004382605  max resid 0.0006791306 
    ## ... Similar to previous best
    ## Run 423 stress 0.0008319199 
    ## Run 424 stress 0.0005982353 
    ## Run 425 stress 0.001023292 
    ## Run 426 stress 0.0007977554 
    ## Run 427 stress 9.943195e-05 
    ## ... Procrustes: rmse 0.0004451328  max resid 0.0007158508 
    ## ... Similar to previous best
    ## Run 428 stress 9.847658e-05 
    ## ... Procrustes: rmse 0.0004181983  max resid 0.0006987553 
    ## ... Similar to previous best
    ## Run 429 stress 9.947156e-05 
    ## ... Procrustes: rmse 0.0004950123  max resid 0.0007798658 
    ## ... Similar to previous best
    ## Run 430 stress 0.0005578463 
    ## ... Procrustes: rmse 0.004914299  max resid 0.00751774 
    ## ... Similar to previous best
    ## Run 431 stress 0.0007904166 
    ## Run 432 stress 9.958813e-05 
    ## ... Procrustes: rmse 0.0004474667  max resid 0.0007017499 
    ## ... Similar to previous best
    ## Run 433 stress 0.0005805997 
    ## ... Procrustes: rmse 0.00509398  max resid 0.007743434 
    ## Run 434 stress 0.0008149534 
    ## Run 435 stress 9.876643e-05 
    ## ... Procrustes: rmse 0.0004450073  max resid 0.0007116181 
    ## ... Similar to previous best
    ## Run 436 stress 0.0007517273 
    ## Run 437 stress 9.650219e-05 
    ## ... Procrustes: rmse 0.0004503598  max resid 0.000726282 
    ## ... Similar to previous best
    ## Run 438 stress 9.901983e-05 
    ## ... Procrustes: rmse 0.0004931838  max resid 0.0007776876 
    ## ... Similar to previous best
    ## Run 439 stress 9.821213e-05 
    ## ... Procrustes: rmse 0.0004261481  max resid 0.0006992434 
    ## ... Similar to previous best
    ## Run 440 stress 9.806718e-05 
    ## ... Procrustes: rmse 0.0004235962  max resid 0.0006160132 
    ## ... Similar to previous best
    ## Run 441 stress 9.698653e-05 
    ## ... Procrustes: rmse 0.0004370162  max resid 0.0006889748 
    ## ... Similar to previous best
    ## Run 442 stress 9.80422e-05 
    ## ... Procrustes: rmse 0.0004418583  max resid 0.0007082288 
    ## ... Similar to previous best
    ## Run 443 stress 9.778493e-05 
    ## ... Procrustes: rmse 0.00042165  max resid 0.0006117104 
    ## ... Similar to previous best
    ## Run 444 stress 9.854463e-05 
    ## ... Procrustes: rmse 0.0004310563  max resid 0.0007123877 
    ## ... Similar to previous best
    ## Run 445 stress 9.801437e-05 
    ## ... Procrustes: rmse 0.0004983615  max resid 0.0007520009 
    ## ... Similar to previous best
    ## Run 446 stress 0.0007086178 
    ## Run 447 stress 9.856361e-05 
    ## ... Procrustes: rmse 0.0004280423  max resid 0.0006778851 
    ## ... Similar to previous best
    ## Run 448 stress 9.961052e-05 
    ## ... Procrustes: rmse 0.0004959384  max resid 0.0007808258 
    ## ... Similar to previous best
    ## Run 449 stress 0.0004456689 
    ## ... Procrustes: rmse 0.00340137  max resid 0.005226944 
    ## ... Similar to previous best
    ## Run 450 stress 9.991444e-05 
    ## ... Procrustes: rmse 0.0004312786  max resid 0.000626988 
    ## ... Similar to previous best
    ## Run 451 stress 9.51995e-05 
    ## ... Procrustes: rmse 0.0004036378  max resid 0.0005894038 
    ## ... Similar to previous best
    ## Run 452 stress 9.752542e-05 
    ## ... Procrustes: rmse 0.0004422307  max resid 0.0006883027 
    ## ... Similar to previous best
    ## Run 453 stress 9.897313e-05 
    ## ... Procrustes: rmse 0.000445349  max resid 0.0007205346 
    ## ... Similar to previous best
    ## Run 454 stress 9.77917e-05 
    ## ... Procrustes: rmse 0.0004725502  max resid 0.0007633783 
    ## ... Similar to previous best
    ## Run 455 stress 9.950895e-05 
    ## ... Procrustes: rmse 0.0004143476  max resid 0.0006102681 
    ## ... Similar to previous best
    ## Run 456 stress 9.77649e-05 
    ## ... Procrustes: rmse 0.0004376904  max resid 0.0007056214 
    ## ... Similar to previous best
    ## Run 457 stress 0.0007919543 
    ## Run 458 stress 9.285447e-05 
    ## ... Procrustes: rmse 0.0004202433  max resid 0.0006755215 
    ## ... Similar to previous best
    ## Run 459 stress 9.970006e-05 
    ## ... Procrustes: rmse 0.000449492  max resid 0.0007191656 
    ## ... Similar to previous best
    ## Run 460 stress 0.0008764249 
    ## Run 461 stress 9.755553e-05 
    ## ... Procrustes: rmse 0.0004138822  max resid 0.0006052987 
    ## ... Similar to previous best
    ## Run 462 stress 0.2348699 
    ## Run 463 stress 0.0006816092 
    ## Run 464 stress 9.897518e-05 
    ## ... Procrustes: rmse 0.0004469014  max resid 0.0007011233 
    ## ... Similar to previous best
    ## Run 465 stress 9.982261e-05 
    ## ... Procrustes: rmse 0.0004824173  max resid 0.0007790018 
    ## ... Similar to previous best
    ## Run 466 stress 9.932423e-05 
    ## ... Procrustes: rmse 0.000449138  max resid 0.0007016276 
    ## ... Similar to previous best
    ## Run 467 stress 9.848745e-05 
    ## ... Procrustes: rmse 0.0004301213  max resid 0.0007121452 
    ## ... Similar to previous best
    ## Run 468 stress 0.0005549417 
    ## ... Procrustes: rmse 0.004871745  max resid 0.007396064 
    ## ... Similar to previous best
    ## Run 469 stress 0.0006560081 
    ## Run 470 stress 9.86103e-05 
    ## ... Procrustes: rmse 0.0004431853  max resid 0.0007098924 
    ## ... Similar to previous best
    ## Run 471 stress 9.851906e-05 
    ## ... Procrustes: rmse 0.0004524481  max resid 0.0007411824 
    ## ... Similar to previous best
    ## Run 472 stress 9.77859e-05 
    ## ... Procrustes: rmse 0.0004761074  max resid 0.0007653758 
    ## ... Similar to previous best
    ## Run 473 stress 0.001119814 
    ## Run 474 stress 9.948025e-05 
    ## ... Procrustes: rmse 0.000444122  max resid 0.0007255738 
    ## ... Similar to previous best
    ## Run 475 stress 9.33582e-05 
    ## ... Procrustes: rmse 0.0004219018  max resid 0.0006682272 
    ## ... Similar to previous best
    ## Run 476 stress 0.0003429503 
    ## ... Procrustes: rmse 0.003021608  max resid 0.004549481 
    ## ... Similar to previous best
    ## Run 477 stress 0.0004670331 
    ## ... Procrustes: rmse 0.004072184  max resid 0.006246817 
    ## ... Similar to previous best
    ## Run 478 stress 0.001174569 
    ## Run 479 stress 9.888634e-05 
    ## ... Procrustes: rmse 0.0004782124  max resid 0.0007728553 
    ## ... Similar to previous best
    ## Run 480 stress 9.878014e-05 
    ## ... Procrustes: rmse 0.0005038471  max resid 0.000763075 
    ## ... Similar to previous best
    ## Run 481 stress 9.937029e-05 
    ## ... Procrustes: rmse 0.000446135  max resid 0.0007003251 
    ## ... Similar to previous best
    ## Run 482 stress 0.0001625873 
    ## ... Procrustes: rmse 0.001369861  max resid 0.002020108 
    ## ... Similar to previous best
    ## Run 483 stress 9.895336e-05 
    ## ... Procrustes: rmse 0.0004484899  max resid 0.0006994135 
    ## ... Similar to previous best
    ## Run 484 stress 9.818843e-05 
    ## ... Procrustes: rmse 0.0004774977  max resid 0.0007678974 
    ## ... Similar to previous best
    ## Run 485 stress 9.81341e-05 
    ## ... Procrustes: rmse 0.000404191  max resid 0.0005938984 
    ## ... Similar to previous best
    ## Run 486 stress 0.0009534207 
    ## Run 487 stress 9.901855e-05 
    ## ... Procrustes: rmse 0.0004016561  max resid 0.000651916 
    ## ... Similar to previous best
    ## Run 488 stress 9.775194e-05 
    ## ... Procrustes: rmse 0.000443416  max resid 0.0006917011 
    ## ... Similar to previous best
    ## Run 489 stress 9.820826e-05 
    ## ... Procrustes: rmse 0.0004226615  max resid 0.0006178885 
    ## ... Similar to previous best
    ## Run 490 stress 0.0008579824 
    ## Run 491 stress 9.771612e-05 
    ## ... Procrustes: rmse 0.000431163  max resid 0.0006879349 
    ## ... Similar to previous best
    ## Run 492 stress 0.0007537004 
    ## Run 493 stress 9.856578e-05 
    ## ... Procrustes: rmse 0.0004906869  max resid 0.0007715956 
    ## ... Similar to previous best
    ## Run 494 stress 0.0007397935 
    ## Run 495 stress 9.736932e-05 
    ## ... Procrustes: rmse 0.000420319  max resid 0.0006134897 
    ## ... Similar to previous best
    ## Run 496 stress 9.907389e-05 
    ## ... Procrustes: rmse 0.0004268286  max resid 0.0006226049 
    ## ... Similar to previous best
    ## Run 497 stress 9.832559e-05 
    ## ... Procrustes: rmse 0.0004428863  max resid 0.0007081999 
    ## ... Similar to previous best
    ## Run 498 stress 9.855243e-05 
    ## ... Procrustes: rmse 0.0004438131  max resid 0.0007096113 
    ## ... Similar to previous best
    ## Run 499 stress 9.863545e-05 
    ## ... Procrustes: rmse 0.0004595144  max resid 0.0007391877 
    ## ... Similar to previous best
    ## Run 500 stress 9.93384e-05 
    ## ... Procrustes: rmse 0.0004630833  max resid 0.0007395115 
    ## ... Similar to previous best
    ## *** Best solution repeated 241 times

    ## Warning in metaMDS(PD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
PD_beta_env_M_NMDS <- metaMDS(PD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.019724e-05 
    ## Run 1 stress 9.537473e-05 
    ## ... Procrustes: rmse 1.59976e-05  max resid 2.434499e-05 
    ## ... Similar to previous best
    ## Run 2 stress 9.429252e-05 
    ## ... Procrustes: rmse 0.0001754409  max resid 0.00035747 
    ## ... Similar to previous best
    ## Run 3 stress 8.730639e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001375923  max resid 0.0002388392 
    ## ... Similar to previous best
    ## Run 4 stress 9.059559e-05 
    ## ... Procrustes: rmse 0.0001954061  max resid 0.0004041779 
    ## ... Similar to previous best
    ## Run 5 stress 9.507648e-05 
    ## ... Procrustes: rmse 0.0001407476  max resid 0.0002485481 
    ## ... Similar to previous best
    ## Run 6 stress 9.80991e-05 
    ## ... Procrustes: rmse 0.0001877254  max resid 0.0003965357 
    ## ... Similar to previous best
    ## Run 7 stress 9.679951e-05 
    ## ... Procrustes: rmse 6.232991e-05  max resid 0.0001007105 
    ## ... Similar to previous best
    ## Run 8 stress 9.62454e-05 
    ## ... Procrustes: rmse 0.0001468971  max resid 0.0002474242 
    ## ... Similar to previous best
    ## Run 9 stress 9.618224e-05 
    ## ... Procrustes: rmse 0.0001432851  max resid 0.0002872792 
    ## ... Similar to previous best
    ## Run 10 stress 0.2319199 
    ## Run 11 stress 9.864719e-05 
    ## ... Procrustes: rmse 0.0001753439  max resid 0.0002942867 
    ## ... Similar to previous best
    ## Run 12 stress 9.623399e-05 
    ## ... Procrustes: rmse 0.000199208  max resid 0.0003169571 
    ## ... Similar to previous best
    ## Run 13 stress 8.421719e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001952911  max resid 0.0003588211 
    ## ... Similar to previous best
    ## Run 14 stress 9.219751e-05 
    ## ... Procrustes: rmse 0.0001315507  max resid 0.0002476805 
    ## ... Similar to previous best
    ## Run 15 stress 9.450743e-05 
    ## ... Procrustes: rmse 0.000178495  max resid 0.0003184261 
    ## ... Similar to previous best
    ## Run 16 stress 8.661974e-05 
    ## ... Procrustes: rmse 0.0002189754  max resid 0.000520822 
    ## ... Similar to previous best
    ## Run 17 stress 8.813599e-05 
    ## ... Procrustes: rmse 0.000196222  max resid 0.0003555288 
    ## ... Similar to previous best
    ## Run 18 stress 9.543529e-05 
    ## ... Procrustes: rmse 0.0002488711  max resid 0.000584996 
    ## ... Similar to previous best
    ## Run 19 stress 9.848036e-05 
    ## ... Procrustes: rmse 0.0002543503  max resid 0.0005970434 
    ## ... Similar to previous best
    ## Run 20 stress 9.602306e-05 
    ## ... Procrustes: rmse 0.000171853  max resid 0.0002519266 
    ## ... Similar to previous best
    ## Run 21 stress 9.633488e-05 
    ## ... Procrustes: rmse 0.0002178194  max resid 0.0003884029 
    ## ... Similar to previous best
    ## Run 22 stress 0.2848515 
    ## Run 23 stress 9.803997e-05 
    ## ... Procrustes: rmse 0.0002093216  max resid 0.0003927776 
    ## ... Similar to previous best
    ## Run 24 stress 9.925465e-05 
    ## ... Procrustes: rmse 0.0001720225  max resid 0.0002487514 
    ## ... Similar to previous best
    ## Run 25 stress 9.807744e-05 
    ## ... Procrustes: rmse 0.0001720796  max resid 0.0002595041 
    ## ... Similar to previous best
    ## Run 26 stress 9.69793e-05 
    ## ... Procrustes: rmse 0.0002136231  max resid 0.0003874572 
    ## ... Similar to previous best
    ## Run 27 stress 9.759726e-05 
    ## ... Procrustes: rmse 0.0002163692  max resid 0.0003848677 
    ## ... Similar to previous best
    ## Run 28 stress 9.085766e-05 
    ## ... Procrustes: rmse 0.0001893675  max resid 0.0003075134 
    ## ... Similar to previous best
    ## Run 29 stress 7.883525e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 8.789541e-05  max resid 0.0001518281 
    ## ... Similar to previous best
    ## Run 30 stress 8.12251e-05 
    ## ... Procrustes: rmse 0.0001683882  max resid 0.0003042901 
    ## ... Similar to previous best
    ## Run 31 stress 9.55471e-05 
    ## ... Procrustes: rmse 0.0001987513  max resid 0.0003170013 
    ## ... Similar to previous best
    ## Run 32 stress 9.209125e-05 
    ## ... Procrustes: rmse 0.0001834295  max resid 0.0004213114 
    ## ... Similar to previous best
    ## Run 33 stress 9.341458e-05 
    ## ... Procrustes: rmse 0.0002257264  max resid 0.0004276272 
    ## ... Similar to previous best
    ## Run 34 stress 8.73525e-05 
    ## ... Procrustes: rmse 0.0002062786  max resid 0.000390548 
    ## ... Similar to previous best
    ## Run 35 stress 9.233351e-05 
    ## ... Procrustes: rmse 0.0002229202  max resid 0.0004231498 
    ## ... Similar to previous best
    ## Run 36 stress 9.20846e-05 
    ## ... Procrustes: rmse 0.0001798622  max resid 0.0002630103 
    ## ... Similar to previous best
    ## Run 37 stress 8.835082e-05 
    ## ... Procrustes: rmse 0.0001053458  max resid 0.0001936392 
    ## ... Similar to previous best
    ## Run 38 stress 8.463243e-05 
    ## ... Procrustes: rmse 0.0001430721  max resid 0.0002130607 
    ## ... Similar to previous best
    ## Run 39 stress 9.947519e-05 
    ## ... Procrustes: rmse 0.0002014365  max resid 0.0004596778 
    ## ... Similar to previous best
    ## Run 40 stress 9.09095e-05 
    ## ... Procrustes: rmse 0.0001542327  max resid 0.0002319511 
    ## ... Similar to previous best
    ## Run 41 stress 9.863931e-05 
    ## ... Procrustes: rmse 0.0001818933  max resid 0.0002404989 
    ## ... Similar to previous best
    ## Run 42 stress 8.737428e-05 
    ## ... Procrustes: rmse 0.0001438366  max resid 0.0002684167 
    ## ... Similar to previous best
    ## Run 43 stress 9.071982e-05 
    ## ... Procrustes: rmse 0.000220214  max resid 0.0004189858 
    ## ... Similar to previous best
    ## Run 44 stress 8.243408e-05 
    ## ... Procrustes: rmse 0.0001405468  max resid 0.0002829403 
    ## ... Similar to previous best
    ## Run 45 stress 9.445036e-05 
    ## ... Procrustes: rmse 0.0002241895  max resid 0.00042775 
    ## ... Similar to previous best
    ## Run 46 stress 8.829062e-05 
    ## ... Procrustes: rmse 0.000163429  max resid 0.000304367 
    ## ... Similar to previous best
    ## Run 47 stress 9.971344e-05 
    ## ... Procrustes: rmse 0.0002366365  max resid 0.0004502755 
    ## ... Similar to previous best
    ## Run 48 stress 9.8063e-05 
    ## ... Procrustes: rmse 0.0001661675  max resid 0.0002449324 
    ## ... Similar to previous best
    ## Run 49 stress 9.347389e-05 
    ## ... Procrustes: rmse 0.0001861789  max resid 0.0004309277 
    ## ... Similar to previous best
    ## Run 50 stress 9.587335e-05 
    ## ... Procrustes: rmse 0.000232645  max resid 0.0004409676 
    ## ... Similar to previous best
    ## Run 51 stress 7.33482e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001252979  max resid 0.0002449237 
    ## ... Similar to previous best
    ## Run 52 stress 7.766039e-05 
    ## ... Procrustes: rmse 0.0001163292  max resid 0.0001817201 
    ## ... Similar to previous best
    ## Run 53 stress 9.763764e-05 
    ## ... Procrustes: rmse 0.000141961  max resid 0.0002418325 
    ## ... Similar to previous best
    ## Run 54 stress 9.345758e-05 
    ## ... Procrustes: rmse 0.0001977198  max resid 0.0003122101 
    ## ... Similar to previous best
    ## Run 55 stress 9.792951e-05 
    ## ... Procrustes: rmse 0.0001397998  max resid 0.0002427438 
    ## ... Similar to previous best
    ## Run 56 stress 9.51638e-05 
    ## ... Procrustes: rmse 5.439726e-05  max resid 7.535405e-05 
    ## ... Similar to previous best
    ## Run 57 stress 8.845027e-05 
    ## ... Procrustes: rmse 0.0001012977  max resid 0.000190775 
    ## ... Similar to previous best
    ## Run 58 stress 9.675017e-05 
    ## ... Procrustes: rmse 0.00013978  max resid 0.0002362164 
    ## ... Similar to previous best
    ## Run 59 stress 9.773055e-05 
    ## ... Procrustes: rmse 0.000187126  max resid 0.0002991122 
    ## ... Similar to previous best
    ## Run 60 stress 9.531443e-05 
    ## ... Procrustes: rmse 0.0001995647  max resid 0.0003108335 
    ## ... Similar to previous best
    ## Run 61 stress 9.468745e-05 
    ## ... Procrustes: rmse 0.0001322387  max resid 0.0002378491 
    ## ... Similar to previous best
    ## Run 62 stress 9.964288e-05 
    ## ... Procrustes: rmse 0.000210751  max resid 0.000324087 
    ## ... Similar to previous best
    ## Run 63 stress 8.863867e-05 
    ## ... Procrustes: rmse 0.0001624297  max resid 0.0002649731 
    ## ... Similar to previous best
    ## Run 64 stress 9.999407e-05 
    ## ... Procrustes: rmse 0.0001461409  max resid 0.0002979317 
    ## ... Similar to previous best
    ## Run 65 stress 9.177916e-05 
    ## ... Procrustes: rmse 0.0001456185  max resid 0.0002895328 
    ## ... Similar to previous best
    ## Run 66 stress 9.687357e-05 
    ## ... Procrustes: rmse 0.0001851166  max resid 0.0003003867 
    ## ... Similar to previous best
    ## Run 67 stress 9.642195e-05 
    ## ... Procrustes: rmse 0.0002019095  max resid 0.0003138358 
    ## ... Similar to previous best
    ## Run 68 stress 9.74805e-05 
    ## ... Procrustes: rmse 0.0001843908  max resid 0.0002983471 
    ## ... Similar to previous best
    ## Run 69 stress 9.318159e-05 
    ## ... Procrustes: rmse 0.000107792  max resid 0.0002503839 
    ## ... Similar to previous best
    ## Run 70 stress 9.886384e-05 
    ## ... Procrustes: rmse 0.0001449348  max resid 0.0002416085 
    ## ... Similar to previous best
    ## Run 71 stress 8.841788e-05 
    ## ... Procrustes: rmse 0.0001166354  max resid 0.0002078972 
    ## ... Similar to previous best
    ## Run 72 stress 8.645601e-05 
    ## ... Procrustes: rmse 0.0001668434  max resid 0.0002829884 
    ## ... Similar to previous best
    ## Run 73 stress 9.499977e-05 
    ## ... Procrustes: rmse 0.0001794478  max resid 0.000288541 
    ## ... Similar to previous best
    ## Run 74 stress 9.630141e-05 
    ## ... Procrustes: rmse 0.0001942896  max resid 0.0002865689 
    ## ... Similar to previous best
    ## Run 75 stress 9.083987e-05 
    ## ... Procrustes: rmse 0.000185554  max resid 0.0003019157 
    ## ... Similar to previous best
    ## Run 76 stress 8.852533e-05 
    ## ... Procrustes: rmse 0.0001882395  max resid 0.0003032139 
    ## ... Similar to previous best
    ## Run 77 stress 9.342675e-05 
    ## ... Procrustes: rmse 0.0001657766  max resid 0.0002704307 
    ## ... Similar to previous best
    ## Run 78 stress 9.572152e-05 
    ## ... Procrustes: rmse 0.0002002254  max resid 0.0003099373 
    ## ... Similar to previous best
    ## Run 79 stress 0.3083098 
    ## Run 80 stress 9.394992e-05 
    ## ... Procrustes: rmse 0.0001343366  max resid 0.0002329409 
    ## ... Similar to previous best
    ## Run 81 stress 9.996941e-05 
    ## ... Procrustes: rmse 0.0001880103  max resid 0.0003049348 
    ## ... Similar to previous best
    ## Run 82 stress 9.714647e-05 
    ## ... Procrustes: rmse 0.0002017011  max resid 0.0003122276 
    ## ... Similar to previous best
    ## Run 83 stress 9.556555e-05 
    ## ... Procrustes: rmse 0.000169377  max resid 0.0002674262 
    ## ... Similar to previous best
    ## Run 84 stress 7.799336e-05 
    ## ... Procrustes: rmse 9.385998e-05  max resid 0.0001870609 
    ## ... Similar to previous best
    ## Run 85 stress 9.961457e-05 
    ## ... Procrustes: rmse 0.0001890185  max resid 0.0002875338 
    ## ... Similar to previous best
    ## Run 86 stress 9.028003e-05 
    ## ... Procrustes: rmse 0.0001905871  max resid 0.0003061865 
    ## ... Similar to previous best
    ## Run 87 stress 9.709879e-05 
    ## ... Procrustes: rmse 0.0001999399  max resid 0.0003163344 
    ## ... Similar to previous best
    ## Run 88 stress 9.994507e-05 
    ## ... Procrustes: rmse 0.0001881881  max resid 0.0003052287 
    ## ... Similar to previous best
    ## Run 89 stress 9.345385e-05 
    ## ... Procrustes: rmse 0.0001728256  max resid 0.0002819205 
    ## ... Similar to previous best
    ## Run 90 stress 9.168055e-05 
    ## ... Procrustes: rmse 0.0001857024  max resid 0.0003001483 
    ## ... Similar to previous best
    ## Run 91 stress 9.50696e-05 
    ## ... Procrustes: rmse 0.000197943  max resid 0.0003076514 
    ## ... Similar to previous best
    ## Run 92 stress 9.97206e-05 
    ## ... Procrustes: rmse 0.0001530144  max resid 0.000316471 
    ## ... Similar to previous best
    ## Run 93 stress 8.845263e-05 
    ## ... Procrustes: rmse 0.0001400083  max resid 0.000249471 
    ## ... Similar to previous best
    ## Run 94 stress 9.392052e-05 
    ## ... Procrustes: rmse 0.0001676881  max resid 0.0002483297 
    ## ... Similar to previous best
    ## Run 95 stress 9.504562e-05 
    ## ... Procrustes: rmse 0.0001979659  max resid 0.0003086277 
    ## ... Similar to previous best
    ## Run 96 stress 9.696962e-05 
    ## ... Procrustes: rmse 0.000201513  max resid 0.00031595 
    ## ... Similar to previous best
    ## Run 97 stress 9.776136e-05 
    ## ... Procrustes: rmse 0.0002042215  max resid 0.0003183761 
    ## ... Similar to previous best
    ## Run 98 stress 8.213098e-05 
    ## ... Procrustes: rmse 0.0001703125  max resid 0.0002737673 
    ## ... Similar to previous best
    ## Run 99 stress 9.217984e-05 
    ## ... Procrustes: rmse 0.0001904481  max resid 0.0002992213 
    ## ... Similar to previous best
    ## Run 100 stress 9.272556e-05 
    ## ... Procrustes: rmse 0.0001779057  max resid 0.0002889272 
    ## ... Similar to previous best
    ## Run 101 stress 8.826017e-05 
    ## ... Procrustes: rmse 5.69089e-05  max resid 0.0001020747 
    ## ... Similar to previous best
    ## Run 102 stress 9.273208e-05 
    ## ... Procrustes: rmse 0.0001938959  max resid 0.000308155 
    ## ... Similar to previous best
    ## Run 103 stress 8.05192e-05 
    ## ... Procrustes: rmse 0.0001560514  max resid 0.0002544793 
    ## ... Similar to previous best
    ## Run 104 stress 8.984598e-05 
    ## ... Procrustes: rmse 9.758686e-05  max resid 0.000226902 
    ## ... Similar to previous best
    ## Run 105 stress 9.785323e-05 
    ## ... Procrustes: rmse 0.0001370641  max resid 0.0002279276 
    ## ... Similar to previous best
    ## Run 106 stress 8.729983e-05 
    ## ... Procrustes: rmse 0.000181279  max resid 0.0002896576 
    ## ... Similar to previous best
    ## Run 107 stress 9.253419e-05 
    ## ... Procrustes: rmse 0.000124421  max resid 0.0002784995 
    ## ... Similar to previous best
    ## Run 108 stress 9.008965e-05 
    ## ... Procrustes: rmse 0.0001703979  max resid 0.0002862719 
    ## ... Similar to previous best
    ## Run 109 stress 7.078721e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 9.726571e-05  max resid 0.000234192 
    ## ... Similar to previous best
    ## Run 110 stress 8.660696e-05 
    ## ... Procrustes: rmse 0.0002267303  max resid 0.000526074 
    ## ... Similar to previous best
    ## Run 111 stress 9.957077e-05 
    ## ... Procrustes: rmse 0.000221079  max resid 0.0003939175 
    ## ... Similar to previous best
    ## Run 112 stress 9.951667e-05 
    ## ... Procrustes: rmse 0.0002558015  max resid 0.0005880265 
    ## ... Similar to previous best
    ## Run 113 stress 7.044147e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001304198  max resid 0.0002116952 
    ## ... Similar to previous best
    ## Run 114 stress 8.88256e-05 
    ## ... Procrustes: rmse 0.0001459038  max resid 0.0002434836 
    ## ... Similar to previous best
    ## Run 115 stress 9.510922e-05 
    ## ... Procrustes: rmse 0.0002085109  max resid 0.0003515792 
    ## ... Similar to previous best
    ## Run 116 stress 9.805532e-05 
    ## ... Procrustes: rmse 0.0002057065  max resid 0.0003338808 
    ## ... Similar to previous best
    ## Run 117 stress 9.609147e-05 
    ## ... Procrustes: rmse 0.0002088065  max resid 0.0003444379 
    ## ... Similar to previous best
    ## Run 118 stress 8.754062e-05 
    ## ... Procrustes: rmse 0.0001448279  max resid 0.0002415679 
    ## ... Similar to previous best
    ## Run 119 stress 9.85077e-05 
    ## ... Procrustes: rmse 0.0002157256  max resid 0.0003570685 
    ## ... Similar to previous best
    ## Run 120 stress 0.3073311 
    ## Run 121 stress 9.800264e-05 
    ## ... Procrustes: rmse 0.000211233  max resid 0.0003494987 
    ## ... Similar to previous best
    ## Run 122 stress 9.724389e-05 
    ## ... Procrustes: rmse 0.0002125973  max resid 0.0003562236 
    ## ... Similar to previous best
    ## Run 123 stress 9.244628e-05 
    ## ... Procrustes: rmse 0.0002040506  max resid 0.0003402683 
    ## ... Similar to previous best
    ## Run 124 stress 9.825863e-05 
    ## ... Procrustes: rmse 0.0002141375  max resid 0.0003590948 
    ## ... Similar to previous best
    ## Run 125 stress 0.3082926 
    ## Run 126 stress 9.368206e-05 
    ## ... Procrustes: rmse 8.835491e-05  max resid 0.000137141 
    ## ... Similar to previous best
    ## Run 127 stress 8.968474e-05 
    ## ... Procrustes: rmse 0.0001954707  max resid 0.0002937339 
    ## ... Similar to previous best
    ## Run 128 stress 9.433665e-05 
    ## ... Procrustes: rmse 0.0001305287  max resid 0.0002278407 
    ## ... Similar to previous best
    ## Run 129 stress 9.325522e-05 
    ## ... Procrustes: rmse 0.0002055651  max resid 0.0003454312 
    ## ... Similar to previous best
    ## Run 130 stress 8.860211e-05 
    ## ... Procrustes: rmse 8.364028e-05  max resid 0.0001341398 
    ## ... Similar to previous best
    ## Run 131 stress 8.632379e-05 
    ## ... Procrustes: rmse 7.843951e-05  max resid 0.000129147 
    ## ... Similar to previous best
    ## Run 132 stress 9.882933e-05 
    ## ... Procrustes: rmse 0.0002118118  max resid 0.0003588239 
    ## ... Similar to previous best
    ## Run 133 stress 9.804694e-05 
    ## ... Procrustes: rmse 0.0002345674  max resid 0.0003470905 
    ## ... Similar to previous best
    ## Run 134 stress 9.943119e-05 
    ## ... Procrustes: rmse 0.0001519485  max resid 0.00025766 
    ## ... Similar to previous best
    ## Run 135 stress 9.745191e-05 
    ## ... Procrustes: rmse 0.0001834508  max resid 0.0002456581 
    ## ... Similar to previous best
    ## Run 136 stress 9.783487e-05 
    ## ... Procrustes: rmse 0.0002297061  max resid 0.000341253 
    ## ... Similar to previous best
    ## Run 137 stress 9.787778e-05 
    ## ... Procrustes: rmse 0.0002071438  max resid 0.0003613167 
    ## ... Similar to previous best
    ## Run 138 stress 9.532892e-05 
    ## ... Procrustes: rmse 0.0002081211  max resid 0.0003481699 
    ## ... Similar to previous best
    ## Run 139 stress 9.877664e-05 
    ## ... Procrustes: rmse 0.0002133252  max resid 0.0003511583 
    ## ... Similar to previous best
    ## Run 140 stress 6.810131e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001126208  max resid 0.0001697202 
    ## ... Similar to previous best
    ## Run 141 stress 8.297312e-05 
    ## ... Procrustes: rmse 4.229814e-05  max resid 7.053545e-05 
    ## ... Similar to previous best
    ## Run 142 stress 9.950244e-05 
    ## ... Procrustes: rmse 0.0002520986  max resid 0.0003743105 
    ## ... Similar to previous best
    ## Run 143 stress 9.042585e-05 
    ## ... Procrustes: rmse 0.0002260355  max resid 0.0003319452 
    ## ... Similar to previous best
    ## Run 144 stress 9.886775e-05 
    ## ... Procrustes: rmse 0.0001855041  max resid 0.0002951011 
    ## ... Similar to previous best
    ## Run 145 stress 9.453305e-05 
    ## ... Procrustes: rmse 0.0001978049  max resid 0.0003443535 
    ## ... Similar to previous best
    ## Run 146 stress 7.636633e-05 
    ## ... Procrustes: rmse 3.548076e-05  max resid 5.959103e-05 
    ## ... Similar to previous best
    ## Run 147 stress 9.199003e-05 
    ## ... Procrustes: rmse 0.0001330605  max resid 0.0002715784 
    ## ... Similar to previous best
    ## Run 148 stress 9.961465e-05 
    ## ... Procrustes: rmse 0.0001542186  max resid 0.0002407497 
    ## ... Similar to previous best
    ## Run 149 stress 9.355464e-05 
    ## ... Procrustes: rmse 0.00023099  max resid 0.0003319516 
    ## ... Similar to previous best
    ## Run 150 stress 9.90492e-05 
    ## ... Procrustes: rmse 0.0002146632  max resid 0.000376647 
    ## ... Similar to previous best
    ## Run 151 stress 8.631327e-05 
    ## ... Procrustes: rmse 0.0001182535  max resid 0.0002346722 
    ## ... Similar to previous best
    ## Run 152 stress 9.133836e-05 
    ## ... Procrustes: rmse 0.0001164755  max resid 0.0002517668 
    ## ... Similar to previous best
    ## Run 153 stress 9.589903e-05 
    ## ... Procrustes: rmse 0.0002389911  max resid 0.0003605934 
    ## ... Similar to previous best
    ## Run 154 stress 9.194819e-05 
    ## ... Procrustes: rmse 0.0001726073  max resid 0.0002725181 
    ## ... Similar to previous best
    ## Run 155 stress 8.951065e-05 
    ## ... Procrustes: rmse 0.0001072512  max resid 0.0001824707 
    ## ... Similar to previous best
    ## Run 156 stress 9.929341e-05 
    ## ... Procrustes: rmse 0.0002142722  max resid 0.0003791456 
    ## ... Similar to previous best
    ## Run 157 stress 0.2680462 
    ## Run 158 stress 9.086351e-05 
    ## ... Procrustes: rmse 0.0001744046  max resid 0.0002874847 
    ## ... Similar to previous best
    ## Run 159 stress 9.232308e-05 
    ## ... Procrustes: rmse 0.0002347074  max resid 0.000346975 
    ## ... Similar to previous best
    ## Run 160 stress 9.267791e-05 
    ## ... Procrustes: rmse 0.000182762  max resid 0.0003306577 
    ## ... Similar to previous best
    ## Run 161 stress 9.796486e-05 
    ## ... Procrustes: rmse 0.0002106532  max resid 0.0003650057 
    ## ... Similar to previous best
    ## Run 162 stress 9.544272e-05 
    ## ... Procrustes: rmse 0.0002062531  max resid 0.0003624562 
    ## ... Similar to previous best
    ## Run 163 stress 8.924096e-05 
    ## ... Procrustes: rmse 0.0001632443  max resid 0.0002938593 
    ## ... Similar to previous best
    ## Run 164 stress 0.3079292 
    ## Run 165 stress 9.963208e-05 
    ## ... Procrustes: rmse 0.0002527395  max resid 0.0003744516 
    ## ... Similar to previous best
    ## Run 166 stress 9.831703e-05 
    ## ... Procrustes: rmse 0.0002483531  max resid 0.0003703496 
    ## ... Similar to previous best
    ## Run 167 stress 9.535414e-05 
    ## ... Procrustes: rmse 0.000146778  max resid 0.0002274312 
    ## ... Similar to previous best
    ## Run 168 stress 9.422561e-05 
    ## ... Procrustes: rmse 0.0002309101  max resid 0.0003322246 
    ## ... Similar to previous best
    ## Run 169 stress 8.755276e-05 
    ## ... Procrustes: rmse 0.0002158916  max resid 0.000333119 
    ## ... Similar to previous best
    ## Run 170 stress 0.3083099 
    ## Run 171 stress 9.111378e-05 
    ## ... Procrustes: rmse 5.461965e-05  max resid 7.772674e-05 
    ## ... Similar to previous best
    ## Run 172 stress 9.404138e-05 
    ## ... Procrustes: rmse 0.0002288277  max resid 0.0003360758 
    ## ... Similar to previous best
    ## Run 173 stress 8.365961e-05 
    ## ... Procrustes: rmse 0.0001132467  max resid 0.0002238182 
    ## ... Similar to previous best
    ## Run 174 stress 8.557394e-05 
    ## ... Procrustes: rmse 0.0001027504  max resid 0.0001752763 
    ## ... Similar to previous best
    ## Run 175 stress 9.451299e-05 
    ## ... Procrustes: rmse 0.0001930414  max resid 0.0003317481 
    ## ... Similar to previous best
    ## Run 176 stress 9.52556e-05 
    ## ... Procrustes: rmse 0.0001131633  max resid 0.0002139006 
    ## ... Similar to previous best
    ## Run 177 stress 8.43968e-05 
    ## ... Procrustes: rmse 0.0001168027  max resid 0.000246977 
    ## ... Similar to previous best
    ## Run 178 stress 9.167831e-05 
    ## ... Procrustes: rmse 0.0001812035  max resid 0.0003224245 
    ## ... Similar to previous best
    ## Run 179 stress 9.997432e-05 
    ## ... Procrustes: rmse 0.0001950502  max resid 0.0003556622 
    ## ... Similar to previous best
    ## Run 180 stress 9.978682e-05 
    ## ... Procrustes: rmse 0.0002451231  max resid 0.0003543088 
    ## ... Similar to previous best
    ## Run 181 stress 9.843922e-05 
    ## ... Procrustes: rmse 0.0002423708  max resid 0.0003480353 
    ## ... Similar to previous best
    ## Run 182 stress 9.83718e-05 
    ## ... Procrustes: rmse 0.000245374  max resid 0.000370561 
    ## ... Similar to previous best
    ## Run 183 stress 9.300307e-05 
    ## ... Procrustes: rmse 0.000110999  max resid 0.0001886454 
    ## ... Similar to previous best
    ## Run 184 stress 9.323214e-05 
    ## ... Procrustes: rmse 0.0001348446  max resid 0.0001976379 
    ## ... Similar to previous best
    ## Run 185 stress 9.698175e-05 
    ## ... Procrustes: rmse 0.0001902831  max resid 0.0003468602 
    ## ... Similar to previous best
    ## Run 186 stress 9.733671e-05 
    ## ... Procrustes: rmse 6.091523e-05  max resid 8.244705e-05 
    ## ... Similar to previous best
    ## Run 187 stress 9.106015e-05 
    ## ... Procrustes: rmse 0.000181103  max resid 0.0003211528 
    ## ... Similar to previous best
    ## Run 188 stress 9.965334e-05 
    ## ... Procrustes: rmse 0.0002461273  max resid 0.0003482756 
    ## ... Similar to previous best
    ## Run 189 stress 9.555286e-05 
    ## ... Procrustes: rmse 0.0001813815  max resid 0.0003183462 
    ## ... Similar to previous best
    ## Run 190 stress 9.556664e-05 
    ## ... Procrustes: rmse 0.0002359977  max resid 0.0003649188 
    ## ... Similar to previous best
    ## Run 191 stress 9.933479e-05 
    ## ... Procrustes: rmse 0.0001877065  max resid 0.000304682 
    ## ... Similar to previous best
    ## Run 192 stress 9.641137e-05 
    ## ... Procrustes: rmse 0.0002002563  max resid 0.0003445157 
    ## ... Similar to previous best
    ## Run 193 stress 9.281989e-05 
    ## ... Procrustes: rmse 0.0001186263  max resid 0.0002196727 
    ## ... Similar to previous best
    ## Run 194 stress 9.476856e-05 
    ## ... Procrustes: rmse 0.0002402783  max resid 0.0003593445 
    ## ... Similar to previous best
    ## Run 195 stress 9.761433e-05 
    ## ... Procrustes: rmse 0.0002114458  max resid 0.0003650507 
    ## ... Similar to previous best
    ## Run 196 stress 0.3083098 
    ## Run 197 stress 0.3083098 
    ## Run 198 stress 9.616249e-05 
    ## ... Procrustes: rmse 0.0002032819  max resid 0.0003608198 
    ## ... Similar to previous best
    ## Run 199 stress 9.47341e-05 
    ## ... Procrustes: rmse 0.0001792041  max resid 0.0002867158 
    ## ... Similar to previous best
    ## Run 200 stress 0.2319186 
    ## Run 201 stress 0.3083098 
    ## Run 202 stress 8.493906e-05 
    ## ... Procrustes: rmse 0.0001339097  max resid 0.0001935826 
    ## ... Similar to previous best
    ## Run 203 stress 8.609266e-05 
    ## ... Procrustes: rmse 0.0001075983  max resid 0.0001777644 
    ## ... Similar to previous best
    ## Run 204 stress 9.408851e-05 
    ## ... Procrustes: rmse 0.0001372783  max resid 0.0003000958 
    ## ... Similar to previous best
    ## Run 205 stress 0.3083098 
    ## Run 206 stress 0.2944291 
    ## Run 207 stress 0.3083098 
    ## Run 208 stress 9.432459e-05 
    ## ... Procrustes: rmse 0.0002395843  max resid 0.0003542418 
    ## ... Similar to previous best
    ## Run 209 stress 9.920674e-05 
    ## ... Procrustes: rmse 0.0001887495  max resid 0.0003039522 
    ## ... Similar to previous best
    ## Run 210 stress 9.766008e-05 
    ## ... Procrustes: rmse 0.0001922414  max resid 0.0003436448 
    ## ... Similar to previous best
    ## Run 211 stress 9.244269e-05 
    ## ... Procrustes: rmse 5.654554e-05  max resid 8.544351e-05 
    ## ... Similar to previous best
    ## Run 212 stress 8.94595e-05 
    ## ... Procrustes: rmse 0.0001138091  max resid 0.0002109829 
    ## ... Similar to previous best
    ## Run 213 stress 9.968071e-05 
    ## ... Procrustes: rmse 0.0002171692  max resid 0.0002860006 
    ## ... Similar to previous best
    ## Run 214 stress 8.676723e-05 
    ## ... Procrustes: rmse 0.0001354151  max resid 0.0001956752 
    ## ... Similar to previous best
    ## Run 215 stress 9.125494e-05 
    ## ... Procrustes: rmse 0.0001961423  max resid 0.0003454168 
    ## ... Similar to previous best
    ## Run 216 stress 9.413918e-05 
    ## ... Procrustes: rmse 0.0002382026  max resid 0.000344485 
    ## ... Similar to previous best
    ## Run 217 stress 9.893092e-05 
    ## ... Procrustes: rmse 4.379585e-05  max resid 5.724962e-05 
    ## ... Similar to previous best
    ## Run 218 stress 8.949398e-05 
    ## ... Procrustes: rmse 0.0001144084  max resid 0.000202121 
    ## ... Similar to previous best
    ## Run 219 stress 9.795273e-05 
    ## ... Procrustes: rmse 0.0001998693  max resid 0.0003424456 
    ## ... Similar to previous best
    ## Run 220 stress 9.144069e-05 
    ## ... Procrustes: rmse 0.0002250989  max resid 0.0003405902 
    ## ... Similar to previous best
    ## Run 221 stress 8.635224e-05 
    ## ... Procrustes: rmse 0.0001275408  max resid 0.0002684662 
    ## ... Similar to previous best
    ## Run 222 stress 7.874918e-05 
    ## ... Procrustes: rmse 3.971243e-05  max resid 5.790418e-05 
    ## ... Similar to previous best
    ## Run 223 stress 9.870473e-05 
    ## ... Procrustes: rmse 0.0001909147  max resid 0.0003163038 
    ## ... Similar to previous best
    ## Run 224 stress 9.650798e-05 
    ## ... Procrustes: rmse 0.0001424606  max resid 0.0002170279 
    ## ... Similar to previous best
    ## Run 225 stress 9.193233e-05 
    ## ... Procrustes: rmse 0.0001176043  max resid 0.0002186532 
    ## ... Similar to previous best
    ## Run 226 stress 9.312311e-05 
    ## ... Procrustes: rmse 0.0001388251  max resid 0.0002130118 
    ## ... Similar to previous best
    ## Run 227 stress 9.085401e-05 
    ## ... Procrustes: rmse 0.000180823  max resid 0.0003206123 
    ## ... Similar to previous best
    ## Run 228 stress 9.120776e-05 
    ## ... Procrustes: rmse 0.0001935663  max resid 0.000339911 
    ## ... Similar to previous best
    ## Run 229 stress 9.367271e-05 
    ## ... Procrustes: rmse 0.0002345416  max resid 0.0003400807 
    ## ... Similar to previous best
    ## Run 230 stress 6.784207e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 2.446798e-05  max resid 3.818177e-05 
    ## ... Similar to previous best
    ## Run 231 stress 9.65011e-05 
    ## ... Procrustes: rmse 0.0001504972  max resid 0.0003158746 
    ## ... Similar to previous best
    ## Run 232 stress 9.207844e-05 
    ## ... Procrustes: rmse 6.937527e-05  max resid 0.0001175229 
    ## ... Similar to previous best
    ## Run 233 stress 9.896765e-05 
    ## ... Procrustes: rmse 0.000181855  max resid 0.0002818014 
    ## ... Similar to previous best
    ## Run 234 stress 8.699309e-05 
    ## ... Procrustes: rmse 6.030003e-05  max resid 0.0001072716 
    ## ... Similar to previous best
    ## Run 235 stress 9.879749e-05 
    ## ... Procrustes: rmse 0.0002144583  max resid 0.0003636592 
    ## ... Similar to previous best
    ## Run 236 stress 9.860991e-05 
    ## ... Procrustes: rmse 0.000243174  max resid 0.0003594361 
    ## ... Similar to previous best
    ## Run 237 stress 8.203278e-05 
    ## ... Procrustes: rmse 0.0001695728  max resid 0.0002777345 
    ## ... Similar to previous best
    ## Run 238 stress 8.842648e-05 
    ## ... Procrustes: rmse 0.0001185225  max resid 0.000232932 
    ## ... Similar to previous best
    ## Run 239 stress 9.408401e-05 
    ## ... Procrustes: rmse 0.0002351595  max resid 0.0003412547 
    ## ... Similar to previous best
    ## Run 240 stress 9.663357e-05 
    ## ... Procrustes: rmse 0.0002420733  max resid 0.0003585176 
    ## ... Similar to previous best
    ## Run 241 stress 9.865423e-05 
    ## ... Procrustes: rmse 0.000218283  max resid 0.0003662698 
    ## ... Similar to previous best
    ## Run 242 stress 9.977333e-05 
    ## ... Procrustes: rmse 0.0002207992  max resid 0.0003635221 
    ## ... Similar to previous best
    ## Run 243 stress 8.851811e-05 
    ## ... Procrustes: rmse 0.0002217553  max resid 0.0003239841 
    ## ... Similar to previous best
    ## Run 244 stress 9.985803e-05 
    ## ... Procrustes: rmse 0.0002433793  max resid 0.0003485548 
    ## ... Similar to previous best
    ## Run 245 stress 7.926785e-05 
    ## ... Procrustes: rmse 0.0001019143  max resid 0.0001740224 
    ## ... Similar to previous best
    ## Run 246 stress 9.707114e-05 
    ## ... Procrustes: rmse 0.0002414795  max resid 0.0003545963 
    ## ... Similar to previous best
    ## Run 247 stress 8.942488e-05 
    ## ... Procrustes: rmse 0.0001195426  max resid 0.0002338111 
    ## ... Similar to previous best
    ## Run 248 stress 8.580427e-05 
    ## ... Procrustes: rmse 0.0001207085  max resid 0.0002160978 
    ## ... Similar to previous best
    ## Run 249 stress 9.740901e-05 
    ## ... Procrustes: rmse 0.0002383127  max resid 0.0003380887 
    ## ... Similar to previous best
    ## Run 250 stress 9.443812e-05 
    ## ... Procrustes: rmse 0.0001394554  max resid 0.000283048 
    ## ... Similar to previous best
    ## Run 251 stress 8.680182e-05 
    ## ... Procrustes: rmse 0.000121376  max resid 0.0002304817 
    ## ... Similar to previous best
    ## Run 252 stress 9.551634e-05 
    ## ... Procrustes: rmse 0.0002088761  max resid 0.0003434841 
    ## ... Similar to previous best
    ## Run 253 stress 0.2680462 
    ## Run 254 stress 8.125624e-05 
    ## ... Procrustes: rmse 0.0001997312  max resid 0.0003000001 
    ## ... Similar to previous best
    ## Run 255 stress 0.2319195 
    ## Run 256 stress 8.808025e-05 
    ## ... Procrustes: rmse 0.0001379922  max resid 0.0002284532 
    ## ... Similar to previous best
    ## Run 257 stress 9.904845e-05 
    ## ... Procrustes: rmse 0.0001369466  max resid 0.0002692866 
    ## ... Similar to previous best
    ## Run 258 stress 9.239769e-05 
    ## ... Procrustes: rmse 0.0001645184  max resid 0.0002560676 
    ## ... Similar to previous best
    ## Run 259 stress 8.463532e-05 
    ## ... Procrustes: rmse 0.0001140886  max resid 0.0002242292 
    ## ... Similar to previous best
    ## Run 260 stress 8.92611e-05 
    ## ... Procrustes: rmse 0.000224079  max resid 0.0003294385 
    ## ... Similar to previous best
    ## Run 261 stress 8.88058e-05 
    ## ... Procrustes: rmse 0.0001042019  max resid 0.0001957611 
    ## ... Similar to previous best
    ## Run 262 stress 9.790977e-05 
    ## ... Procrustes: rmse 0.0001309323  max resid 0.000253409 
    ## ... Similar to previous best
    ## Run 263 stress 9.500785e-05 
    ## ... Procrustes: rmse 0.0002358363  max resid 0.0003496595 
    ## ... Similar to previous best
    ## Run 264 stress 9.177095e-05 
    ## ... Procrustes: rmse 0.0001251359  max resid 0.0002449517 
    ## ... Similar to previous best
    ## Run 265 stress 9.380433e-05 
    ## ... Procrustes: rmse 0.0002313374  max resid 0.0003395015 
    ## ... Similar to previous best
    ## Run 266 stress 8.952283e-05 
    ## ... Procrustes: rmse 0.0001101618  max resid 0.000204447 
    ## ... Similar to previous best
    ## Run 267 stress 7.377685e-05 
    ## ... Procrustes: rmse 4.014988e-05  max resid 7.113127e-05 
    ## ... Similar to previous best
    ## Run 268 stress 9.696355e-05 
    ## ... Procrustes: rmse 0.0002413743  max resid 0.0003596266 
    ## ... Similar to previous best
    ## Run 269 stress 9.460635e-05 
    ## ... Procrustes: rmse 0.000209345  max resid 0.0003519193 
    ## ... Similar to previous best
    ## Run 270 stress 9.222594e-05 
    ## ... Procrustes: rmse 0.0001177726  max resid 0.0002089777 
    ## ... Similar to previous best
    ## Run 271 stress 9.266055e-05 
    ## ... Procrustes: rmse 0.0001242472  max resid 0.0002458974 
    ## ... Similar to previous best
    ## Run 272 stress 9.994813e-05 
    ## ... Procrustes: rmse 7.948991e-05  max resid 0.0001185277 
    ## ... Similar to previous best
    ## Run 273 stress 9.849073e-05 
    ## ... Procrustes: rmse 0.0002147284  max resid 0.0003638992 
    ## ... Similar to previous best
    ## Run 274 stress 9.927696e-05 
    ## ... Procrustes: rmse 0.0001967653  max resid 0.0003568815 
    ## ... Similar to previous best
    ## Run 275 stress 9.699183e-05 
    ## ... Procrustes: rmse 0.0001949035  max resid 0.0003548602 
    ## ... Similar to previous best
    ## Run 276 stress 9.882132e-05 
    ## ... Procrustes: rmse 0.0001651884  max resid 0.0002799055 
    ## ... Similar to previous best
    ## Run 277 stress 9.342891e-05 
    ## ... Procrustes: rmse 0.0001826278  max resid 0.0003226357 
    ## ... Similar to previous best
    ## Run 278 stress 8.926523e-05 
    ## ... Procrustes: rmse 0.0001793431  max resid 0.0003067955 
    ## ... Similar to previous best
    ## Run 279 stress 9.760888e-05 
    ## ... Procrustes: rmse 0.0002164256  max resid 0.0003551035 
    ## ... Similar to previous best
    ## Run 280 stress 9.711794e-05 
    ## ... Procrustes: rmse 0.0001931932  max resid 0.0003501516 
    ## ... Similar to previous best
    ## Run 281 stress 9.663515e-05 
    ## ... Procrustes: rmse 0.0002379144  max resid 0.0003345302 
    ## ... Similar to previous best
    ## Run 282 stress 9.001643e-05 
    ## ... Procrustes: rmse 0.0001821258  max resid 0.0003226111 
    ## ... Similar to previous best
    ## Run 283 stress 9.38665e-05 
    ## ... Procrustes: rmse 0.0002362557  max resid 0.000344348 
    ## ... Similar to previous best
    ## Run 284 stress 9.640312e-05 
    ## ... Procrustes: rmse 0.0002414668  max resid 0.000354613 
    ## ... Similar to previous best
    ## Run 285 stress 8.673899e-05 
    ## ... Procrustes: rmse 0.0002125404  max resid 0.0003065402 
    ## ... Similar to previous best
    ## Run 286 stress 9.934167e-05 
    ## ... Procrustes: rmse 0.0002416867  max resid 0.0003641755 
    ## ... Similar to previous best
    ## Run 287 stress 6.658806e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 8.070487e-05  max resid 0.0001289153 
    ## ... Similar to previous best
    ## Run 288 stress 9.576635e-05 
    ## ... Procrustes: rmse 0.000176121  max resid 0.0002527543 
    ## ... Similar to previous best
    ## Run 289 stress 8.717815e-05 
    ## ... Procrustes: rmse 0.0001298104  max resid 0.0001863607 
    ## ... Similar to previous best
    ## Run 290 stress 9.929212e-05 
    ## ... Procrustes: rmse 0.0001612963  max resid 0.0002865555 
    ## ... Similar to previous best
    ## Run 291 stress 9.738198e-05 
    ## ... Procrustes: rmse 0.0002091159  max resid 0.0003590564 
    ## ... Similar to previous best
    ## Run 292 stress 9.973456e-05 
    ## ... Procrustes: rmse 0.000237537  max resid 0.0003782428 
    ## ... Similar to previous best
    ## Run 293 stress 9.828933e-05 
    ## ... Procrustes: rmse 0.0001284247  max resid 0.0001911511 
    ## ... Similar to previous best
    ## Run 294 stress 9.154274e-05 
    ## ... Procrustes: rmse 0.0001446696  max resid 0.0002551925 
    ## ... Similar to previous best
    ## Run 295 stress 9.821165e-05 
    ## ... Procrustes: rmse 0.0001427637  max resid 0.0002039935 
    ## ... Similar to previous best
    ## Run 296 stress 8.751544e-05 
    ## ... Procrustes: rmse 0.0001501784  max resid 0.0002873691 
    ## ... Similar to previous best
    ## Run 297 stress 9.613988e-05 
    ## ... Procrustes: rmse 0.0002072727  max resid 0.0003567054 
    ## ... Similar to previous best
    ## Run 298 stress 9.695816e-05 
    ## ... Procrustes: rmse 0.0001558986  max resid 0.0002482218 
    ## ... Similar to previous best
    ## Run 299 stress 8.938626e-05 
    ## ... Procrustes: rmse 0.0001543652  max resid 0.0003278732 
    ## ... Similar to previous best
    ## Run 300 stress 9.647022e-05 
    ## ... Procrustes: rmse 0.0001582744  max resid 0.0002776352 
    ## ... Similar to previous best
    ## Run 301 stress 0.3074301 
    ## Run 302 stress 8.071729e-05 
    ## ... Procrustes: rmse 0.0001148791  max resid 0.0001811861 
    ## ... Similar to previous best
    ## Run 303 stress 9.873296e-05 
    ## ... Procrustes: rmse 0.0002222141  max resid 0.0003776041 
    ## ... Similar to previous best
    ## Run 304 stress 9.543638e-05 
    ## ... Procrustes: rmse 0.0001786678  max resid 0.0002606759 
    ## ... Similar to previous best
    ## Run 305 stress 9.886105e-05 
    ## ... Procrustes: rmse 0.0002190466  max resid 0.0003655591 
    ## ... Similar to previous best
    ## Run 306 stress 9.616193e-05 
    ## ... Procrustes: rmse 0.0001512872  max resid 0.0002668065 
    ## ... Similar to previous best
    ## Run 307 stress 9.505279e-05 
    ## ... Procrustes: rmse 0.0001505593  max resid 0.0002660452 
    ## ... Similar to previous best
    ## Run 308 stress 9.488718e-05 
    ## ... Procrustes: rmse 0.0002056851  max resid 0.0003519122 
    ## ... Similar to previous best
    ## Run 309 stress 9.844174e-05 
    ## ... Procrustes: rmse 0.000215514  max resid 0.0003721686 
    ## ... Similar to previous best
    ## Run 310 stress 9.175658e-05 
    ## ... Procrustes: rmse 0.0002001335  max resid 0.000332716 
    ## ... Similar to previous best
    ## Run 311 stress 9.560363e-05 
    ## ... Procrustes: rmse 0.0002121182  max resid 0.0003625813 
    ## ... Similar to previous best
    ## Run 312 stress 9.625893e-05 
    ## ... Procrustes: rmse 0.0002292458  max resid 0.0003640628 
    ## ... Similar to previous best
    ## Run 313 stress 9.341552e-05 
    ## ... Procrustes: rmse 0.0001515783  max resid 0.0002444548 
    ## ... Similar to previous best
    ## Run 314 stress 9.871198e-05 
    ## ... Procrustes: rmse 0.0002209618  max resid 0.0003733595 
    ## ... Similar to previous best
    ## Run 315 stress 9.747554e-05 
    ## ... Procrustes: rmse 0.0002176561  max resid 0.0003710004 
    ## ... Similar to previous best
    ## Run 316 stress 8.968331e-05 
    ## ... Procrustes: rmse 0.0001648544  max resid 0.0002480801 
    ## ... Similar to previous best
    ## Run 317 stress 9.197159e-05 
    ## ... Procrustes: rmse 0.0002022993  max resid 0.0003437427 
    ## ... Similar to previous best
    ## Run 318 stress 9.546077e-05 
    ## ... Procrustes: rmse 0.0001172018  max resid 0.0001817531 
    ## ... Similar to previous best
    ## Run 319 stress 9.501726e-05 
    ## ... Procrustes: rmse 0.0002138894  max resid 0.0003487714 
    ## ... Similar to previous best
    ## Run 320 stress 9.425721e-05 
    ## ... Procrustes: rmse 0.000210529  max resid 0.0003577618 
    ## ... Similar to previous best
    ## Run 321 stress 9.098225e-05 
    ## ... Procrustes: rmse 0.0001428664  max resid 0.0002522811 
    ## ... Similar to previous best
    ## Run 322 stress 9.970939e-05 
    ## ... Procrustes: rmse 0.0001806839  max resid 0.0002826483 
    ## ... Similar to previous best
    ## Run 323 stress 9.459692e-05 
    ## ... Procrustes: rmse 0.0001479219  max resid 0.0002605261 
    ## ... Similar to previous best
    ## Run 324 stress 9.488229e-05 
    ## ... Procrustes: rmse 0.0002263778  max resid 0.0003571485 
    ## ... Similar to previous best
    ## Run 325 stress 9.460563e-05 
    ## ... Procrustes: rmse 0.0002123521  max resid 0.0003569837 
    ## ... Similar to previous best
    ## Run 326 stress 9.843808e-05 
    ## ... Procrustes: rmse 0.0001546709  max resid 0.0002728391 
    ## ... Similar to previous best
    ## Run 327 stress 9.428662e-05 
    ## ... Procrustes: rmse 0.0001517863  max resid 0.0002259193 
    ## ... Similar to previous best
    ## Run 328 stress 0.2848523 
    ## Run 329 stress 9.942811e-05 
    ## ... Procrustes: rmse 0.0002142289  max resid 0.0003681205 
    ## ... Similar to previous best
    ## Run 330 stress 8.504836e-05 
    ## ... Procrustes: rmse 0.0001869814  max resid 0.0003104455 
    ## ... Similar to previous best
    ## Run 331 stress 9.859551e-05 
    ## ... Procrustes: rmse 0.0001433815  max resid 0.0002095367 
    ## ... Similar to previous best
    ## Run 332 stress 9.384821e-05 
    ## ... Procrustes: rmse 0.0001759546  max resid 0.0002573152 
    ## ... Similar to previous best
    ## Run 333 stress 5.991221e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 8.356727e-05  max resid 0.0001641437 
    ## ... Similar to previous best
    ## Run 334 stress 8.305374e-05 
    ## ... Procrustes: rmse 8.348235e-05  max resid 0.0001396491 
    ## ... Similar to previous best
    ## Run 335 stress 8.970249e-05 
    ## ... Procrustes: rmse 0.0001306765  max resid 0.000270288 
    ## ... Similar to previous best
    ## Run 336 stress 9.7038e-05 
    ## ... Procrustes: rmse 0.0001039316  max resid 0.0001924558 
    ## ... Similar to previous best
    ## Run 337 stress 9.24986e-05 
    ## ... Procrustes: rmse 0.0001039469  max resid 0.0001864174 
    ## ... Similar to previous best
    ## Run 338 stress 0.3083098 
    ## Run 339 stress 9.404718e-05 
    ## ... Procrustes: rmse 0.0001415887  max resid 0.000294345 
    ## ... Similar to previous best
    ## Run 340 stress 9.570517e-05 
    ## ... Procrustes: rmse 0.0001928965  max resid 0.000317744 
    ## ... Similar to previous best
    ## Run 341 stress 9.830571e-05 
    ## ... Procrustes: rmse 0.0001856642  max resid 0.0002854573 
    ## ... Similar to previous best
    ## Run 342 stress 0.285215 
    ## Run 343 stress 9.817939e-05 
    ## ... Procrustes: rmse 0.0001479733  max resid 0.0003093392 
    ## ... Similar to previous best
    ## Run 344 stress 9.245173e-05 
    ## ... Procrustes: rmse 0.0001131148  max resid 0.000203614 
    ## ... Similar to previous best
    ## Run 345 stress 9.38414e-05 
    ## ... Procrustes: rmse 0.0002090403  max resid 0.0003041685 
    ## ... Similar to previous best
    ## Run 346 stress 8.772049e-05 
    ## ... Procrustes: rmse 0.0001250759  max resid 0.0002557504 
    ## ... Similar to previous best
    ## Run 347 stress 9.710803e-05 
    ## ... Procrustes: rmse 0.0002211073  max resid 0.00032828 
    ## ... Similar to previous best
    ## Run 348 stress 9.657778e-05 
    ## ... Procrustes: rmse 0.0001911565  max resid 0.0003179517 
    ## ... Similar to previous best
    ## Run 349 stress 9.976667e-05 
    ## ... Procrustes: rmse 0.0001255998  max resid 0.0002135558 
    ## ... Similar to previous best
    ## Run 350 stress 9.1654e-05 
    ## ... Procrustes: rmse 0.0001325392  max resid 0.0002670074 
    ## ... Similar to previous best
    ## Run 351 stress 9.111028e-05 
    ## ... Procrustes: rmse 0.000203879  max resid 0.0002951942 
    ## ... Similar to previous best
    ## Run 352 stress 0.2852145 
    ## Run 353 stress 9.858899e-05 
    ## ... Procrustes: rmse 0.000221545  max resid 0.0003240996 
    ## ... Similar to previous best
    ## Run 354 stress 9.678327e-05 
    ## ... Procrustes: rmse 0.0002191651  max resid 0.0003193868 
    ## ... Similar to previous best
    ## Run 355 stress 9.603865e-05 
    ## ... Procrustes: rmse 0.0001463596  max resid 0.0002567889 
    ## ... Similar to previous best
    ## Run 356 stress 9.965928e-05 
    ## ... Procrustes: rmse 0.0001992297  max resid 0.0003276042 
    ## ... Similar to previous best
    ## Run 357 stress 8.992844e-05 
    ## ... Procrustes: rmse 0.0001783539  max resid 0.000295031 
    ## ... Similar to previous best
    ## Run 358 stress 9.92687e-05 
    ## ... Procrustes: rmse 0.0002269148  max resid 0.0003352705 
    ## ... Similar to previous best
    ## Run 359 stress 9.71231e-05 
    ## ... Procrustes: rmse 0.0001086573  max resid 0.0001813192 
    ## ... Similar to previous best
    ## Run 360 stress 9.706519e-05 
    ## ... Procrustes: rmse 0.0001908702  max resid 0.0003184148 
    ## ... Similar to previous best
    ## Run 361 stress 9.918446e-05 
    ## ... Procrustes: rmse 0.0001894393  max resid 0.0002915912 
    ## ... Similar to previous best
    ## Run 362 stress 8.860351e-05 
    ## ... Procrustes: rmse 0.0001460981  max resid 0.0002714159 
    ## ... Similar to previous best
    ## Run 363 stress 8.661185e-05 
    ## ... Procrustes: rmse 0.0001356729  max resid 0.0002610361 
    ## ... Similar to previous best
    ## Run 364 stress 8.83345e-05 
    ## ... Procrustes: rmse 0.0001287453  max resid 0.0002659482 
    ## ... Similar to previous best
    ## Run 365 stress 9.416065e-05 
    ## ... Procrustes: rmse 0.0002119584  max resid 0.000312768 
    ## ... Similar to previous best
    ## Run 366 stress 8.973815e-05 
    ## ... Procrustes: rmse 0.0001345298  max resid 0.0002757632 
    ## ... Similar to previous best
    ## Run 367 stress 0.2319135 
    ## Run 368 stress 9.368647e-05 
    ## ... Procrustes: rmse 0.0002115257  max resid 0.0003105372 
    ## ... Similar to previous best
    ## Run 369 stress 9.790892e-05 
    ## ... Procrustes: rmse 0.0001848738  max resid 0.000288797 
    ## ... Similar to previous best
    ## Run 370 stress 9.464337e-05 
    ## ... Procrustes: rmse 0.0002136005  max resid 0.0003103043 
    ## ... Similar to previous best
    ## Run 371 stress 9.400926e-05 
    ## ... Procrustes: rmse 0.0001380228  max resid 0.0002861893 
    ## ... Similar to previous best
    ## Run 372 stress 8.612211e-05 
    ## ... Procrustes: rmse 0.0001076255  max resid 0.0001556143 
    ## ... Similar to previous best
    ## Run 373 stress 8.606851e-05 
    ## ... Procrustes: rmse 9.755909e-05  max resid 0.0001801043 
    ## ... Similar to previous best
    ## Run 374 stress 9.096064e-05 
    ## ... Procrustes: rmse 0.0001678815  max resid 0.0002877315 
    ## ... Similar to previous best
    ## Run 375 stress 8.486296e-05 
    ## ... Procrustes: rmse 0.0001199342  max resid 0.0002401242 
    ## ... Similar to previous best
    ## Run 376 stress 9.171819e-05 
    ## ... Procrustes: rmse 0.0001779575  max resid 0.0002946255 
    ## ... Similar to previous best
    ## Run 377 stress 9.412533e-05 
    ## ... Procrustes: rmse 0.0001597122  max resid 0.000273145 
    ## ... Similar to previous best
    ## Run 378 stress 0.2852159 
    ## Run 379 stress 8.647488e-05 
    ## ... Procrustes: rmse 0.0001904247  max resid 0.0002681663 
    ## ... Similar to previous best
    ## Run 380 stress 8.604628e-05 
    ## ... Procrustes: rmse 0.0001334163  max resid 0.0002779463 
    ## ... Similar to previous best
    ## Run 381 stress 9.412859e-05 
    ## ... Procrustes: rmse 0.0001504417  max resid 0.0002832786 
    ## ... Similar to previous best
    ## Run 382 stress 9.542794e-05 
    ## ... Procrustes: rmse 0.0001917506  max resid 0.0003192713 
    ## ... Similar to previous best
    ## Run 383 stress 9.205543e-05 
    ## ... Procrustes: rmse 0.0002042221  max resid 0.0003016306 
    ## ... Similar to previous best
    ## Run 384 stress 8.876281e-05 
    ## ... Procrustes: rmse 0.0001494006  max resid 0.0002775017 
    ## ... Similar to previous best
    ## Run 385 stress 8.136777e-05 
    ## ... Procrustes: rmse 0.0001007453  max resid 0.000152076 
    ## ... Similar to previous best
    ## Run 386 stress 9.789032e-05 
    ## ... Procrustes: rmse 0.0001830309  max resid 0.0002804196 
    ## ... Similar to previous best
    ## Run 387 stress 9.310053e-05 
    ## ... Procrustes: rmse 0.0002056791  max resid 0.0002987797 
    ## ... Similar to previous best
    ## Run 388 stress 9.064711e-05 
    ## ... Procrustes: rmse 7.97353e-05  max resid 0.0001526285 
    ## ... Similar to previous best
    ## Run 389 stress 9.377641e-05 
    ## ... Procrustes: rmse 0.0001857951  max resid 0.0003036876 
    ## ... Similar to previous best
    ## Run 390 stress 9.238011e-05 
    ## ... Procrustes: rmse 0.0001755176  max resid 0.0002701507 
    ## ... Similar to previous best
    ## Run 391 stress 0.2985159 
    ## Run 392 stress 9.565567e-05 
    ## ... Procrustes: rmse 0.0002064502  max resid 0.0002939944 
    ## ... Similar to previous best
    ## Run 393 stress 0.3083098 
    ## Run 394 stress 9.045754e-05 
    ## ... Procrustes: rmse 0.0001795784  max resid 0.0002943954 
    ## ... Similar to previous best
    ## Run 395 stress 9.368439e-05 
    ## ... Procrustes: rmse 0.0001777364  max resid 0.0002735088 
    ## ... Similar to previous best
    ## Run 396 stress 9.963447e-05 
    ## ... Procrustes: rmse 0.0001544503  max resid 0.0003253305 
    ## ... Similar to previous best
    ## Run 397 stress 9.742215e-05 
    ## ... Procrustes: rmse 0.0002079719  max resid 0.000296044 
    ## ... Similar to previous best
    ## Run 398 stress 9.038793e-05 
    ## ... Procrustes: rmse 0.00015042  max resid 0.0002803856 
    ## ... Similar to previous best
    ## Run 399 stress 9.030363e-05 
    ## ... Procrustes: rmse 0.0002004196  max resid 0.0002857622 
    ## ... Similar to previous best
    ## Run 400 stress 9.541239e-05 
    ## ... Procrustes: rmse 0.0001891832  max resid 0.0003146474 
    ## ... Similar to previous best
    ## Run 401 stress 9.189625e-05 
    ## ... Procrustes: rmse 0.0001743047  max resid 0.0002681061 
    ## ... Similar to previous best
    ## Run 402 stress 7.739716e-05 
    ## ... Procrustes: rmse 0.0001478341  max resid 0.0002345884 
    ## ... Similar to previous best
    ## Run 403 stress 9.862609e-05 
    ## ... Procrustes: rmse 0.0002241219  max resid 0.0003302126 
    ## ... Similar to previous best
    ## Run 404 stress 8.510964e-05 
    ## ... Procrustes: rmse 0.0001239998  max resid 0.0001849363 
    ## ... Similar to previous best
    ## Run 405 stress 9.119678e-05 
    ## ... Procrustes: rmse 0.0001739285  max resid 0.0002805492 
    ## ... Similar to previous best
    ## Run 406 stress 9.357335e-05 
    ## ... Procrustes: rmse 0.0001803888  max resid 0.000283446 
    ## ... Similar to previous best
    ## Run 407 stress 8.916736e-05 
    ## ... Procrustes: rmse 0.0001358245  max resid 0.000281111 
    ## ... Similar to previous best
    ## Run 408 stress 9.521966e-05 
    ## ... Procrustes: rmse 0.0001388419  max resid 0.0001980554 
    ## ... Similar to previous best
    ## Run 409 stress 9.209605e-05 
    ## ... Procrustes: rmse 9.349217e-05  max resid 0.000113153 
    ## ... Similar to previous best
    ## Run 410 stress 9.134475e-05 
    ## ... Procrustes: rmse 0.0001208073  max resid 0.0001908362 
    ## ... Similar to previous best
    ## Run 411 stress 9.760801e-05 
    ## ... Procrustes: rmse 0.0002029245  max resid 0.0002876033 
    ## ... Similar to previous best
    ## Run 412 stress 9.361411e-05 
    ## ... Procrustes: rmse 0.0002085688  max resid 0.0002997723 
    ## ... Similar to previous best
    ## Run 413 stress 9.173204e-05 
    ## ... Procrustes: rmse 0.000183813  max resid 0.0003016259 
    ## ... Similar to previous best
    ## Run 414 stress 9.530491e-05 
    ## ... Procrustes: rmse 0.0001661935  max resid 0.0002750554 
    ## ... Similar to previous best
    ## Run 415 stress 9.823338e-05 
    ## ... Procrustes: rmse 0.0001498031  max resid 0.0002570796 
    ## ... Similar to previous best
    ## Run 416 stress 8.934073e-05 
    ## ... Procrustes: rmse 0.0001566341  max resid 0.0002637584 
    ## ... Similar to previous best
    ## Run 417 stress 9.511381e-05 
    ## ... Procrustes: rmse 0.0001211873  max resid 0.0002033651 
    ## ... Similar to previous best
    ## Run 418 stress 8.740437e-05 
    ## ... Procrustes: rmse 0.000133401  max resid 0.0002775109 
    ## ... Similar to previous best
    ## Run 419 stress 9.084704e-05 
    ## ... Procrustes: rmse 0.000169271  max resid 0.0002813626 
    ## ... Similar to previous best
    ## Run 420 stress 9.904939e-05 
    ## ... Procrustes: rmse 0.0001030671  max resid 0.0001362583 
    ## ... Similar to previous best
    ## Run 421 stress 9.894782e-05 
    ## ... Procrustes: rmse 0.0001810309  max resid 0.0002950394 
    ## ... Similar to previous best
    ## Run 422 stress 8.469738e-05 
    ## ... Procrustes: rmse 9.840449e-05  max resid 0.0001747833 
    ## ... Similar to previous best
    ## Run 423 stress 9.484929e-05 
    ## ... Procrustes: rmse 0.0001164571  max resid 0.0002120461 
    ## ... Similar to previous best
    ## Run 424 stress 9.735059e-05 
    ## ... Procrustes: rmse 9.521436e-05  max resid 0.0001667088 
    ## ... Similar to previous best
    ## Run 425 stress 9.784555e-05 
    ## ... Procrustes: rmse 0.0001777989  max resid 0.0002689043 
    ## ... Similar to previous best
    ## Run 426 stress 8.871845e-05 
    ## ... Procrustes: rmse 0.000128425  max resid 0.0001900698 
    ## ... Similar to previous best
    ## Run 427 stress 6.050019e-05 
    ## ... Procrustes: rmse 0.0001046195  max resid 0.000170179 
    ## ... Similar to previous best
    ## Run 428 stress 9.633046e-05 
    ## ... Procrustes: rmse 0.0001950014  max resid 0.0003221079 
    ## ... Similar to previous best
    ## Run 429 stress 9.662086e-05 
    ## ... Procrustes: rmse 0.0001835147  max resid 0.0002865636 
    ## ... Similar to previous best
    ## Run 430 stress 7.959642e-05 
    ## ... Procrustes: rmse 8.017994e-05  max resid 0.0001219987 
    ## ... Similar to previous best
    ## Run 431 stress 0.3083098 
    ## Run 432 stress 8.395288e-05 
    ## ... Procrustes: rmse 0.0001324128  max resid 0.0002730571 
    ## ... Similar to previous best
    ## Run 433 stress 9.983191e-05 
    ## ... Procrustes: rmse 0.0002027492  max resid 0.000336445 
    ## ... Similar to previous best
    ## Run 434 stress 9.365989e-05 
    ## ... Procrustes: rmse 0.0001148634  max resid 0.0002060847 
    ## ... Similar to previous best
    ## Run 435 stress 8.327058e-05 
    ## ... Procrustes: rmse 0.0001764032  max resid 0.0002491339 
    ## ... Similar to previous best
    ## Run 436 stress 0.2842805 
    ## Run 437 stress 9.945486e-05 
    ## ... Procrustes: rmse 0.000198397  max resid 0.0003272858 
    ## ... Similar to previous best
    ## Run 438 stress 9.6518e-05 
    ## ... Procrustes: rmse 0.0001920626  max resid 0.0003149569 
    ## ... Similar to previous best
    ## Run 439 stress 9.584322e-05 
    ## ... Procrustes: rmse 0.0002170388  max resid 0.0003206735 
    ## ... Similar to previous best
    ## Run 440 stress 8.347728e-05 
    ## ... Procrustes: rmse 8.856715e-05  max resid 0.0001418843 
    ## ... Similar to previous best
    ## Run 441 stress 9.390631e-05 
    ## ... Procrustes: rmse 0.000115966  max resid 0.0002023947 
    ## ... Similar to previous best
    ## Run 442 stress 8.680474e-05 
    ## ... Procrustes: rmse 0.0001919669  max resid 0.0002662546 
    ## ... Similar to previous best
    ## Run 443 stress 9.622133e-05 
    ## ... Procrustes: rmse 0.0001915937  max resid 0.000315153 
    ## ... Similar to previous best
    ## Run 444 stress 9.778046e-05 
    ## ... Procrustes: rmse 0.0001854243  max resid 0.0002853852 
    ## ... Similar to previous best
    ## Run 445 stress 9.160732e-05 
    ## ... Procrustes: rmse 0.0001747396  max resid 0.0002692167 
    ## ... Similar to previous best
    ## Run 446 stress 9.96312e-05 
    ## ... Procrustes: rmse 0.0002255628  max resid 0.0003348156 
    ## ... Similar to previous best
    ## Run 447 stress 8.660348e-05 
    ## ... Procrustes: rmse 0.0001000192  max resid 0.0001874705 
    ## ... Similar to previous best
    ## Run 448 stress 8.422559e-05 
    ## ... Procrustes: rmse 8.348183e-05  max resid 0.0001334998 
    ## ... Similar to previous best
    ## Run 449 stress 9.505576e-05 
    ## ... Procrustes: rmse 0.0001838191  max resid 0.0003063303 
    ## ... Similar to previous best
    ## Run 450 stress 9.388033e-05 
    ## ... Procrustes: rmse 0.0002108214  max resid 0.0003124842 
    ## ... Similar to previous best
    ## Run 451 stress 9.895426e-05 
    ## ... Procrustes: rmse 0.0001188999  max resid 0.0002181517 
    ## ... Similar to previous best
    ## Run 452 stress 9.266752e-05 
    ## ... Procrustes: rmse 0.0001538666  max resid 0.0002646121 
    ## ... Similar to previous best
    ## Run 453 stress 0.2319185 
    ## Run 454 stress 0.284852 
    ## Run 455 stress 9.533237e-05 
    ## ... Procrustes: rmse 0.0001720595  max resid 0.0002790148 
    ## ... Similar to previous best
    ## Run 456 stress 8.946853e-05 
    ## ... Procrustes: rmse 0.0001571047  max resid 0.0002610202 
    ## ... Similar to previous best
    ## Run 457 stress 9.33089e-05 
    ## ... Procrustes: rmse 0.0001504773  max resid 0.0002744153 
    ## ... Similar to previous best
    ## Run 458 stress 9.453926e-05 
    ## ... Procrustes: rmse 0.0001789998  max resid 0.0002799317 
    ## ... Similar to previous best
    ## Run 459 stress 8.715522e-05 
    ## ... Procrustes: rmse 0.0001269397  max resid 0.0002549157 
    ## ... Similar to previous best
    ## Run 460 stress 9.881421e-05 
    ## ... Procrustes: rmse 0.0002153647  max resid 0.0003188625 
    ## ... Similar to previous best
    ## Run 461 stress 9.374983e-05 
    ## ... Procrustes: rmse 0.0001098971  max resid 0.0002047557 
    ## ... Similar to previous best
    ## Run 462 stress 8.884148e-05 
    ## ... Procrustes: rmse 0.00013913  max resid 0.0002900216 
    ## ... Similar to previous best
    ## Run 463 stress 9.39144e-05 
    ## ... Procrustes: rmse 0.0001061134  max resid 0.0001961262 
    ## ... Similar to previous best
    ## Run 464 stress 9.855941e-05 
    ## ... Procrustes: rmse 0.0002211488  max resid 0.0003231066 
    ## ... Similar to previous best
    ## Run 465 stress 8.92301e-05 
    ## ... Procrustes: rmse 0.0001278762  max resid 0.000190023 
    ## ... Similar to previous best
    ## Run 466 stress 9.877822e-05 
    ## ... Procrustes: rmse 0.0001254002  max resid 0.0002134308 
    ## ... Similar to previous best
    ## Run 467 stress 9.943219e-05 
    ## ... Procrustes: rmse 9.32734e-05  max resid 0.0001288891 
    ## ... Similar to previous best
    ## Run 468 stress 9.791328e-05 
    ## ... Procrustes: rmse 0.0001689944  max resid 0.0002789541 
    ## ... Similar to previous best
    ## Run 469 stress 9.408852e-05 
    ## ... Procrustes: rmse 0.0001001355  max resid 0.0001312821 
    ## ... Similar to previous best
    ## Run 470 stress 9.19494e-05 
    ## ... Procrustes: rmse 0.0001139091  max resid 0.0001950482 
    ## ... Similar to previous best
    ## Run 471 stress 9.127909e-05 
    ## ... Procrustes: rmse 0.0002039096  max resid 0.0002993706 
    ## ... Similar to previous best
    ## Run 472 stress 9.00792e-05 
    ## ... Procrustes: rmse 0.0001686401  max resid 0.0002683335 
    ## ... Similar to previous best
    ## Run 473 stress 9.493102e-05 
    ## ... Procrustes: rmse 0.0002152975  max resid 0.0003168699 
    ## ... Similar to previous best
    ## Run 474 stress 8.228428e-05 
    ## ... Procrustes: rmse 0.0001068575  max resid 0.0002070688 
    ## ... Similar to previous best
    ## Run 475 stress 9.410959e-05 
    ## ... Procrustes: rmse 0.0001877046  max resid 0.00031081 
    ## ... Similar to previous best
    ## Run 476 stress 9.001698e-05 
    ## ... Procrustes: rmse 0.0001552247  max resid 0.0002896776 
    ## ... Similar to previous best
    ## Run 477 stress 9.848355e-05 
    ## ... Procrustes: rmse 0.0002118425  max resid 0.0003031233 
    ## ... Similar to previous best
    ## Run 478 stress 8.540172e-05 
    ## ... Procrustes: rmse 8.855187e-05  max resid 0.0001357936 
    ## ... Similar to previous best
    ## Run 479 stress 9.155031e-05 
    ## ... Procrustes: rmse 7.559552e-05  max resid 0.0001413918 
    ## ... Similar to previous best
    ## Run 480 stress 9.320455e-05 
    ## ... Procrustes: rmse 0.0002121713  max resid 0.0003076404 
    ## ... Similar to previous best
    ## Run 481 stress 0.2848524 
    ## Run 482 stress 9.769264e-05 
    ## ... Procrustes: rmse 0.0002161969  max resid 0.0003217076 
    ## ... Similar to previous best
    ## Run 483 stress 9.754275e-05 
    ## ... Procrustes: rmse 0.0001932221  max resid 0.0003164211 
    ## ... Similar to previous best
    ## Run 484 stress 9.466059e-05 
    ## ... Procrustes: rmse 0.0002095259  max resid 0.0003114556 
    ## ... Similar to previous best
    ## Run 485 stress 8.59515e-05 
    ## ... Procrustes: rmse 9.523203e-05  max resid 0.0001788649 
    ## ... Similar to previous best
    ## Run 486 stress 9.883334e-05 
    ## ... Procrustes: rmse 0.0002221761  max resid 0.0003305836 
    ## ... Similar to previous best
    ## Run 487 stress 9.445561e-05 
    ## ... Procrustes: rmse 0.0001691802  max resid 0.0002664513 
    ## ... Similar to previous best
    ## Run 488 stress 9.44673e-05 
    ## ... Procrustes: rmse 0.0001825024  max resid 0.0003002598 
    ## ... Similar to previous best
    ## Run 489 stress 9.700394e-05 
    ## ... Procrustes: rmse 0.0002205019  max resid 0.0003266634 
    ## ... Similar to previous best
    ## Run 490 stress 9.722613e-05 
    ## ... Procrustes: rmse 0.0001829287  max resid 0.0002866725 
    ## ... Similar to previous best
    ## Run 491 stress 8.621098e-05 
    ## ... Procrustes: rmse 0.0001902111  max resid 0.0002689538 
    ## ... Similar to previous best
    ## Run 492 stress 9.523475e-05 
    ## ... Procrustes: rmse 0.0001701654  max resid 0.000296016 
    ## ... Similar to previous best
    ## Run 493 stress 9.472535e-05 
    ## ... Procrustes: rmse 0.0001904423  max resid 0.0003172625 
    ## ... Similar to previous best
    ## Run 494 stress 9.293096e-05 
    ## ... Procrustes: rmse 0.0001864239  max resid 0.0003059728 
    ## ... Similar to previous best
    ## Run 495 stress 8.963936e-05 
    ## ... Procrustes: rmse 0.0001006168  max resid 0.0001860619 
    ## ... Similar to previous best
    ## Run 496 stress 9.56881e-05 
    ## ... Procrustes: rmse 0.0002152277  max resid 0.0003144129 
    ## ... Similar to previous best
    ## Run 497 stress 9.904391e-05 
    ## ... Procrustes: rmse 0.0002162884  max resid 0.0003122714 
    ## ... Similar to previous best
    ## Run 498 stress 9.517155e-05 
    ## ... Procrustes: rmse 0.0001366478  max resid 0.0002038204 
    ## ... Similar to previous best
    ## Run 499 stress 9.305955e-05 
    ## ... Procrustes: rmse 0.0001853289  max resid 0.0003074421 
    ## ... Similar to previous best
    ## Run 500 stress 6.956137e-05 
    ## ... Procrustes: rmse 6.717576e-05  max resid 0.0001212229 
    ## ... Similar to previous best
    ## *** Best solution repeated 156 times

    ## Warning in metaMDS(PD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
PD_beta_geo_NMDS <- metaMDS(PD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07927535 
    ## Run 1 stress 0.093821 
    ## Run 2 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 3.134224e-05  max resid 7.618093e-05 
    ## ... Similar to previous best
    ## Run 3 stress 0.09252215 
    ## Run 4 stress 0.07927534 
    ## ... Procrustes: rmse 8.091947e-06  max resid 2.015137e-05 
    ## ... Similar to previous best
    ## Run 5 stress 0.07951445 
    ## ... Procrustes: rmse 0.01737005  max resid 0.05658854 
    ## Run 6 stress 0.0948569 
    ## Run 7 stress 0.07927534 
    ## ... Procrustes: rmse 9.737443e-06  max resid 3.672448e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.08985605 
    ## Run 9 stress 0.07927537 
    ## ... Procrustes: rmse 6.787324e-05  max resid 0.0001887379 
    ## ... Similar to previous best
    ## Run 10 stress 0.09106027 
    ## Run 11 stress 0.09106027 
    ## Run 12 stress 0.09344087 
    ## Run 13 stress 0.09072901 
    ## Run 14 stress 0.07927534 
    ## ... Procrustes: rmse 1.160187e-05  max resid 2.877141e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.07927535 
    ## ... Procrustes: rmse 3.629015e-05  max resid 8.832669e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.0898561 
    ## Run 17 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733636  max resid 0.05632225 
    ## Run 18 stress 0.09022796 
    ## Run 19 stress 0.07936951 
    ## ... Procrustes: rmse 0.009049803  max resid 0.03248676 
    ## Run 20 stress 0.0793695 
    ## ... Procrustes: rmse 0.009048346  max resid 0.03248772 
    ## Run 21 stress 0.07951443 
    ## ... Procrustes: rmse 0.017361  max resid 0.05652939 
    ## Run 22 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237483  max resid 0.0497106 
    ## Run 23 stress 0.09072913 
    ## Run 24 stress 0.09296239 
    ## Run 25 stress 0.1061419 
    ## Run 26 stress 0.09388349 
    ## Run 27 stress 0.09072899 
    ## Run 28 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734918  max resid 0.05638677 
    ## Run 29 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734048  max resid 0.0563424 
    ## Run 30 stress 0.08985612 
    ## Run 31 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237813  max resid 0.04968653 
    ## Run 32 stress 0.09485764 
    ## Run 33 stress 0.09382063 
    ## Run 34 stress 0.0793695 
    ## ... Procrustes: rmse 0.009035998  max resid 0.03241055 
    ## Run 35 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047853  max resid 0.03248293 
    ## Run 36 stress 0.0796926 
    ## ... Procrustes: rmse 0.01239519  max resid 0.04982793 
    ## Run 37 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173597  max resid 0.05646359 
    ## Run 38 stress 0.07951445 
    ## ... Procrustes: rmse 0.01732557  max resid 0.05622725 
    ## Run 39 stress 0.07969272 
    ## ... Procrustes: rmse 0.01241928  max resid 0.05000918 
    ## Run 40 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732698  max resid 0.05624697 
    ## Run 41 stress 0.07951446 
    ## ... Procrustes: rmse 0.01732144  max resid 0.05619569 
    ## Run 42 stress 0.1135254 
    ## Run 43 stress 0.07927534 
    ## ... Procrustes: rmse 1.177242e-05  max resid 4.028153e-05 
    ## ... Similar to previous best
    ## Run 44 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041272  max resid 0.0324551 
    ## Run 45 stress 0.09252223 
    ## Run 46 stress 0.07927536 
    ## ... Procrustes: rmse 5.503878e-05  max resid 0.0001263571 
    ## ... Similar to previous best
    ## Run 47 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734751  max resid 0.05638803 
    ## Run 48 stress 0.09106033 
    ## Run 49 stress 0.07927538 
    ## ... Procrustes: rmse 7.387985e-05  max resid 0.0001587886 
    ## ... Similar to previous best
    ## Run 50 stress 0.08985611 
    ## Run 51 stress 0.07927534 
    ## ... Procrustes: rmse 4.003521e-06  max resid 1.619394e-05 
    ## ... Similar to previous best
    ## Run 52 stress 0.09106029 
    ## Run 53 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237259  max resid 0.04966391 
    ## Run 54 stress 0.0910604 
    ## Run 55 stress 0.09565373 
    ## Run 56 stress 0.09106028 
    ## Run 57 stress 0.07936952 
    ## ... Procrustes: rmse 0.009033532  max resid 0.03239656 
    ## Run 58 stress 0.09022776 
    ## Run 59 stress 0.07927534 
    ## ... Procrustes: rmse 1.676115e-05  max resid 3.973777e-05 
    ## ... Similar to previous best
    ## Run 60 stress 0.07969267 
    ## ... Procrustes: rmse 0.01236689  max resid 0.04955018 
    ## Run 61 stress 0.09072893 
    ## Run 62 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734231  max resid 0.0563906 
    ## Run 63 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735064  max resid 0.05633806 
    ## Run 64 stress 0.0796926 
    ## ... Procrustes: rmse 0.01239641  max resid 0.04983268 
    ## Run 65 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 7.211344e-06  max resid 2.033751e-05 
    ## ... Similar to previous best
    ## Run 66 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734909  max resid 0.05634067 
    ## Run 67 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173475  max resid 0.05643408 
    ## Run 68 stress 0.08985607 
    ## Run 69 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734157  max resid 0.05634812 
    ## Run 70 stress 0.07951446 
    ## ... Procrustes: rmse 0.01733553  max resid 0.05630083 
    ## Run 71 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045601  max resid 0.03247615 
    ## Run 72 stress 0.091728 
    ## Run 73 stress 0.07951444 
    ## ... Procrustes: rmse 0.01734782  max resid 0.05638935 
    ## Run 74 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735933  max resid 0.05645835 
    ## Run 75 stress 0.09296187 
    ## Run 76 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173495  max resid 0.0563986 
    ## Run 77 stress 0.08985603 
    ## Run 78 stress 0.07927536 
    ## ... Procrustes: rmse 4.694178e-05  max resid 0.0001145426 
    ## ... Similar to previous best
    ## Run 79 stress 0.09072896 
    ## Run 80 stress 0.0929071 
    ## Run 81 stress 0.09106036 
    ## Run 82 stress 0.09344088 
    ## Run 83 stress 0.09544362 
    ## Run 84 stress 0.08985605 
    ## Run 85 stress 0.0938208 
    ## Run 86 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043383  max resid 0.03244973 
    ## Run 87 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735471  max resid 0.05643194 
    ## Run 88 stress 0.09106058 
    ## Run 89 stress 0.09474087 
    ## Run 90 stress 0.09172796 
    ## Run 91 stress 0.09072902 
    ## Run 92 stress 0.07927535 
    ## ... Procrustes: rmse 3.413204e-05  max resid 8.422672e-05 
    ## ... Similar to previous best
    ## Run 93 stress 0.09388344 
    ## Run 94 stress 0.09072893 
    ## Run 95 stress 0.07969267 
    ## ... Procrustes: rmse 0.01235516  max resid 0.04949584 
    ## Run 96 stress 0.09565394 
    ## Run 97 stress 0.0793695 
    ## ... Procrustes: rmse 0.009050162  max resid 0.03249832 
    ## Run 98 stress 0.07951441 
    ## ... Procrustes: rmse 0.017362  max resid 0.05647786 
    ## Run 99 stress 0.09172795 
    ## Run 100 stress 0.1061856 
    ## Run 101 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734766  max resid 0.05638841 
    ## Run 102 stress 0.07969259 
    ## ... Procrustes: rmse 0.01236944  max resid 0.04969303 
    ## Run 103 stress 0.07951443 
    ## ... Procrustes: rmse 0.01734915  max resid 0.05639563 
    ## Run 104 stress 0.09022784 
    ## Run 105 stress 0.07969258 
    ## ... Procrustes: rmse 0.01236917  max resid 0.04966454 
    ## Run 106 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734826  max resid 0.05638771 
    ## Run 107 stress 0.09252213 
    ## Run 108 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046113  max resid 0.03247368 
    ## Run 109 stress 0.07927534 
    ## ... Procrustes: rmse 2.092462e-05  max resid 5.322174e-05 
    ## ... Similar to previous best
    ## Run 110 stress 0.07969264 
    ## ... Procrustes: rmse 0.01236719  max resid 0.04957847 
    ## Run 111 stress 0.08985606 
    ## Run 112 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734734  max resid 0.05638313 
    ## Run 113 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735407  max resid 0.05642661 
    ## Run 114 stress 0.07936952 
    ## ... Procrustes: rmse 0.009033721  max resid 0.03239381 
    ## Run 115 stress 0.09106033 
    ## Run 116 stress 0.08985604 
    ## Run 117 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735615  max resid 0.05644213 
    ## Run 118 stress 0.07927534 
    ## ... Procrustes: rmse 1.219957e-06  max resid 2.49096e-06 
    ## ... Similar to previous best
    ## Run 119 stress 0.09022778 
    ## Run 120 stress 0.07927534 
    ## ... Procrustes: rmse 2.638447e-06  max resid 6.416471e-06 
    ## ... Similar to previous best
    ## Run 121 stress 0.07936954 
    ## ... Procrustes: rmse 0.009026805  max resid 0.03235017 
    ## Run 122 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732574  max resid 0.05623305 
    ## Run 123 stress 0.1041898 
    ## Run 124 stress 0.09043747 
    ## Run 125 stress 0.0795145 
    ## ... Procrustes: rmse 0.01731923  max resid 0.05617529 
    ## Run 126 stress 0.08985607 
    ## Run 127 stress 0.07927535 
    ## ... Procrustes: rmse 3.7853e-05  max resid 0.0001029319 
    ## ... Similar to previous best
    ## Run 128 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237456  max resid 0.04971985 
    ## Run 129 stress 0.09388336 
    ## Run 130 stress 0.07951444 
    ## ... Procrustes: rmse 0.01737948  max resid 0.05659482 
    ## Run 131 stress 0.09388327 
    ## Run 132 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732528  max resid 0.05626597 
    ## Run 133 stress 0.09252225 
    ## Run 134 stress 0.09043742 
    ## Run 135 stress 0.09290693 
    ## Run 136 stress 0.09296234 
    ## Run 137 stress 0.0793695 
    ## ... Procrustes: rmse 0.0090399  max resid 0.03243263 
    ## Run 138 stress 0.0793695 
    ## ... Procrustes: rmse 0.009042415  max resid 0.03244698 
    ## Run 139 stress 0.1083562 
    ## Run 140 stress 0.09043748 
    ## Run 141 stress 0.09247102 
    ## Run 142 stress 0.08985606 
    ## Run 143 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046011  max resid 0.03247299 
    ## Run 144 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736044  max resid 0.05646343 
    ## Run 145 stress 0.07936952 
    ## ... Procrustes: rmse 0.009043119  max resid 0.03245063 
    ## Run 146 stress 0.09072897 
    ## Run 147 stress 0.07951443 
    ## ... Procrustes: rmse 0.01735078  max resid 0.05640894 
    ## Run 148 stress 0.09474018 
    ## Run 149 stress 0.09290688 
    ## Run 150 stress 0.0904375 
    ## Run 151 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732539  max resid 0.0562319 
    ## Run 152 stress 0.07927535 
    ## ... Procrustes: rmse 3.794376e-05  max resid 9.309531e-05 
    ## ... Similar to previous best
    ## Run 153 stress 0.09382098 
    ## Run 154 stress 0.07951442 
    ## ... Procrustes: rmse 0.01736787  max resid 0.05652909 
    ## Run 155 stress 0.09388341 
    ## Run 156 stress 0.09290695 
    ## Run 157 stress 0.07936951 
    ## ... Procrustes: rmse 0.009033291  max resid 0.03239102 
    ## Run 158 stress 0.09344092 
    ## Run 159 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047489  max resid 0.03248229 
    ## Run 160 stress 0.07936951 
    ## ... Procrustes: rmse 0.009036465  max resid 0.03240826 
    ## Run 161 stress 0.07936951 
    ## ... Procrustes: rmse 0.009034869  max resid 0.03240064 
    ## Run 162 stress 0.07927536 
    ## ... Procrustes: rmse 4.401583e-05  max resid 0.0001206719 
    ## ... Similar to previous best
    ## Run 163 stress 0.07927538 
    ## ... Procrustes: rmse 7.322765e-05  max resid 0.0002071305 
    ## ... Similar to previous best
    ## Run 164 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733946  max resid 0.05632846 
    ## Run 165 stress 0.07927539 
    ## ... Procrustes: rmse 8.062843e-05  max resid 0.0002291873 
    ## ... Similar to previous best
    ## Run 166 stress 0.07927535 
    ## ... Procrustes: rmse 3.040186e-05  max resid 8.106491e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.0907291 
    ## Run 168 stress 0.07927535 
    ## ... Procrustes: rmse 1.972877e-05  max resid 3.654146e-05 
    ## ... Similar to previous best
    ## Run 169 stress 0.09290712 
    ## Run 170 stress 0.08985604 
    ## Run 171 stress 0.07927538 
    ## ... Procrustes: rmse 7.037801e-05  max resid 0.0001973263 
    ## ... Similar to previous best
    ## Run 172 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734698  max resid 0.05638271 
    ## Run 173 stress 0.09072897 
    ## Run 174 stress 0.08985612 
    ## Run 175 stress 0.07936952 
    ## ... Procrustes: rmse 0.0090617  max resid 0.03255352 
    ## Run 176 stress 0.0898562 
    ## Run 177 stress 0.09043788 
    ## Run 178 stress 0.09043756 
    ## Run 179 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237179  max resid 0.04962692 
    ## Run 180 stress 0.1060928 
    ## Run 181 stress 0.07927534 
    ## ... Procrustes: rmse 3.057976e-06  max resid 5.935719e-06 
    ## ... Similar to previous best
    ## Run 182 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733019  max resid 0.05632044 
    ## Run 183 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173465  max resid 0.05637981 
    ## Run 184 stress 0.09388365 
    ## Run 185 stress 0.07927536 
    ## ... Procrustes: rmse 3.701501e-05  max resid 7.442801e-05 
    ## ... Similar to previous best
    ## Run 186 stress 0.08985605 
    ## Run 187 stress 0.07969259 
    ## ... Procrustes: rmse 0.01236764  max resid 0.04963979 
    ## Run 188 stress 0.09247102 
    ## Run 189 stress 0.09106042 
    ## Run 190 stress 0.09022785 
    ## Run 191 stress 0.07936956 
    ## ... Procrustes: rmse 0.009022362  max resid 0.03232313 
    ## Run 192 stress 0.07927536 
    ## ... Procrustes: rmse 4.889816e-05  max resid 0.000111223 
    ## ... Similar to previous best
    ## Run 193 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173498  max resid 0.05639962 
    ## Run 194 stress 0.09106029 
    ## Run 195 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046747  max resid 0.03247968 
    ## Run 196 stress 0.09106041 
    ## Run 197 stress 0.09344085 
    ## Run 198 stress 0.0904375 
    ## Run 199 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173515  max resid 0.0564404 
    ## Run 200 stress 0.08985605 
    ## Run 201 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732828  max resid 0.05625315 
    ## Run 202 stress 0.07927534 
    ## ... Procrustes: rmse 3.995461e-06  max resid 1.000068e-05 
    ## ... Similar to previous best
    ## Run 203 stress 0.09072892 
    ## Run 204 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238351  max resid 0.04975028 
    ## Run 205 stress 0.08985603 
    ## Run 206 stress 0.07927535 
    ## ... Procrustes: rmse 3.869179e-05  max resid 0.0001079869 
    ## ... Similar to previous best
    ## Run 207 stress 0.09022788 
    ## Run 208 stress 0.09022791 
    ## Run 209 stress 0.09072891 
    ## Run 210 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734285  max resid 0.05635868 
    ## Run 211 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237313  max resid 0.04964722 
    ## Run 212 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734029  max resid 0.05634424 
    ## Run 213 stress 0.07936956 
    ## ... Procrustes: rmse 0.009023071  max resid 0.03232724 
    ## Run 214 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734883  max resid 0.0564066 
    ## Run 215 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237144  max resid 0.04964531 
    ## Run 216 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904488  max resid 0.0324584 
    ## Run 217 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237497  max resid 0.0497177 
    ## Run 218 stress 0.09296231 
    ## Run 219 stress 0.07951445 
    ## ... Procrustes: rmse 0.01738538  max resid 0.05662403 
    ## Run 220 stress 0.09252225 
    ## Run 221 stress 0.09565377 
    ## Run 222 stress 0.07969261 
    ## ... Procrustes: rmse 0.01236821  max resid 0.04960253 
    ## Run 223 stress 0.07927536 
    ## ... Procrustes: rmse 4.519956e-05  max resid 0.0001091989 
    ## ... Similar to previous best
    ## Run 224 stress 0.08985608 
    ## Run 225 stress 0.09072893 
    ## Run 226 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734604  max resid 0.05637397 
    ## Run 227 stress 0.108875 
    ## Run 228 stress 0.09022778 
    ## Run 229 stress 0.09043769 
    ## Run 230 stress 0.09252214 
    ## Run 231 stress 0.07927536 
    ## ... Procrustes: rmse 4.933023e-05  max resid 0.0001379682 
    ## ... Similar to previous best
    ## Run 232 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734425  max resid 0.05636809 
    ## Run 233 stress 0.07927535 
    ## ... Procrustes: rmse 1.090985e-05  max resid 2.707091e-05 
    ## ... Similar to previous best
    ## Run 234 stress 0.08985605 
    ## Run 235 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039168  max resid 0.03242722 
    ## Run 236 stress 0.09544393 
    ## Run 237 stress 0.09296265 
    ## Run 238 stress 0.08985609 
    ## Run 239 stress 0.09072906 
    ## Run 240 stress 0.09022788 
    ## Run 241 stress 0.07927535 
    ## ... Procrustes: rmse 1.438103e-05  max resid 3.060021e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.07969265 
    ## ... Procrustes: rmse 0.01242018  max resid 0.04997692 
    ## Run 243 stress 0.07927538 
    ## ... Procrustes: rmse 7.365104e-05  max resid 0.0002084123 
    ## ... Similar to previous best
    ## Run 244 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734454  max resid 0.05637004 
    ## Run 245 stress 0.07927536 
    ## ... Procrustes: rmse 4.127185e-05  max resid 0.0001077873 
    ## ... Similar to previous best
    ## Run 246 stress 0.09296231 
    ## Run 247 stress 0.09072891 
    ## Run 248 stress 0.09565395 
    ## Run 249 stress 0.09252228 
    ## Run 250 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238105  max resid 0.04973543 
    ## Run 251 stress 0.07927534 
    ## ... Procrustes: rmse 1.561195e-05  max resid 4.002591e-05 
    ## ... Similar to previous best
    ## Run 252 stress 0.07936951 
    ## ... Procrustes: rmse 0.0090521  max resid 0.03250172 
    ## Run 253 stress 0.09290717 
    ## Run 254 stress 0.0904377 
    ## Run 255 stress 0.07927534 
    ## ... Procrustes: rmse 3.751484e-06  max resid 7.095435e-06 
    ## ... Similar to previous best
    ## Run 256 stress 0.07927534 
    ## ... Procrustes: rmse 7.477023e-06  max resid 1.965332e-05 
    ## ... Similar to previous best
    ## Run 257 stress 0.09296217 
    ## Run 258 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173572  max resid 0.0564473 
    ## Run 259 stress 0.09106043 
    ## Run 260 stress 0.07969259 
    ## ... Procrustes: rmse 0.01236794  max resid 0.04962878 
    ## Run 261 stress 0.0925223 
    ## Run 262 stress 0.07927535 
    ## ... Procrustes: rmse 1.972833e-05  max resid 5.396618e-05 
    ## ... Similar to previous best
    ## Run 263 stress 0.09247102 
    ## Run 264 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039961  max resid 0.03243127 
    ## Run 265 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736626  max resid 0.05650506 
    ## Run 266 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237254  max resid 0.04968093 
    ## Run 267 stress 0.09022786 
    ## Run 268 stress 0.07969261 
    ## ... Procrustes: rmse 0.01238087  max resid 0.04975918 
    ## Run 269 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733224  max resid 0.0562785 
    ## Run 270 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735716  max resid 0.05637622 
    ## Run 271 stress 0.09290692 
    ## Run 272 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046107  max resid 0.03246758 
    ## Run 273 stress 0.07927535 
    ## ... Procrustes: rmse 2.856875e-05  max resid 8.356893e-05 
    ## ... Similar to previous best
    ## Run 274 stress 0.08985611 
    ## Run 275 stress 0.07927535 
    ## ... Procrustes: rmse 3.217626e-05  max resid 8.687329e-05 
    ## ... Similar to previous best
    ## Run 276 stress 0.07927534 
    ## ... Procrustes: rmse 1.548982e-05  max resid 4.249522e-05 
    ## ... Similar to previous best
    ## Run 277 stress 0.09022775 
    ## Run 278 stress 0.07927535 
    ## ... Procrustes: rmse 2.720959e-05  max resid 7.189689e-05 
    ## ... Similar to previous best
    ## Run 279 stress 0.08985614 
    ## Run 280 stress 0.07927534 
    ## ... Procrustes: rmse 1.353087e-05  max resid 3.633683e-05 
    ## ... Similar to previous best
    ## Run 281 stress 0.07927538 
    ## ... Procrustes: rmse 7.676072e-05  max resid 0.0002171191 
    ## ... Similar to previous best
    ## Run 282 stress 0.07969259 
    ## ... Procrustes: rmse 0.01235843  max resid 0.04961136 
    ## Run 283 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173512  max resid 0.05641897 
    ## Run 284 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735161  max resid 0.05642124 
    ## Run 285 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041771  max resid 0.03244776 
    ## Run 286 stress 0.0793695 
    ## ... Procrustes: rmse 0.009048915  max resid 0.03248112 
    ## Run 287 stress 0.09043735 
    ## Run 288 stress 0.0793695 
    ## ... Procrustes: rmse 0.009042265  max resid 0.03246579 
    ## Run 289 stress 0.07927534 
    ## ... Procrustes: rmse 1.469409e-05  max resid 4.011295e-05 
    ## ... Similar to previous best
    ## Run 290 stress 0.07927535 
    ## ... Procrustes: rmse 2.875937e-05  max resid 7.801495e-05 
    ## ... Similar to previous best
    ## Run 291 stress 0.07951444 
    ## ... Procrustes: rmse 0.01737138  max resid 0.05654916 
    ## Run 292 stress 0.07969264 
    ## ... Procrustes: rmse 0.01235728  max resid 0.04952338 
    ## Run 293 stress 0.09072895 
    ## Run 294 stress 0.07936952 
    ## ... Procrustes: rmse 0.009031843  max resid 0.03237977 
    ## Run 295 stress 0.09290691 
    ## Run 296 stress 0.08985605 
    ## Run 297 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734968  max resid 0.05638777 
    ## Run 298 stress 0.07927534 
    ## ... Procrustes: rmse 5.733382e-06  max resid 1.463139e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.09388337 
    ## Run 300 stress 0.09072895 
    ## Run 301 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734717  max resid 0.05638168 
    ## Run 302 stress 0.07927538 
    ## ... Procrustes: rmse 7.129404e-05  max resid 0.0002011068 
    ## ... Similar to previous best
    ## Run 303 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041597  max resid 0.03244313 
    ## Run 304 stress 0.07927535 
    ## ... Procrustes: rmse 2.80623e-05  max resid 7.898943e-05 
    ## ... Similar to previous best
    ## Run 305 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735398  max resid 0.05643695 
    ## Run 306 stress 0.07969262 
    ## ... Procrustes: rmse 0.01236996  max resid 0.04961152 
    ## Run 307 stress 0.07936951 
    ## ... Procrustes: rmse 0.009036218  max resid 0.0324068 
    ## Run 308 stress 0.09106029 
    ## Run 309 stress 0.09382066 
    ## Run 310 stress 0.07927534 
    ## ... Procrustes: rmse 2.361469e-05  max resid 6.527357e-05 
    ## ... Similar to previous best
    ## Run 311 stress 0.09388349 
    ## Run 312 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735601  max resid 0.05643782 
    ## Run 313 stress 0.07951443 
    ## ... Procrustes: rmse 0.01735799  max resid 0.05651515 
    ## Run 314 stress 0.09296223 
    ## Run 315 stress 0.09072909 
    ## Run 316 stress 0.09022802 
    ## Run 317 stress 0.09072896 
    ## Run 318 stress 0.07936952 
    ## ... Procrustes: rmse 0.009032206  max resid 0.03238429 
    ## Run 319 stress 0.09043772 
    ## Run 320 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046351  max resid 0.03247654 
    ## Run 321 stress 0.08985604 
    ## Run 322 stress 0.07951445 
    ## ... Procrustes: rmse 0.01731982  max resid 0.05619387 
    ## Run 323 stress 0.1079616 
    ## Run 324 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238965  max resid 0.04979789 
    ## Run 325 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046031  max resid 0.03247242 
    ## Run 326 stress 0.07927534 
    ## ... Procrustes: rmse 1.641088e-05  max resid 3.420242e-05 
    ## ... Similar to previous best
    ## Run 327 stress 0.07969267 
    ## ... Procrustes: rmse 0.01236907  max resid 0.0496005 
    ## Run 328 stress 0.07927538 
    ## ... Procrustes: rmse 7.933724e-05  max resid 0.0002249493 
    ## ... Similar to previous best
    ## Run 329 stress 0.09106052 
    ## Run 330 stress 0.09072901 
    ## Run 331 stress 0.07927535 
    ## ... Procrustes: rmse 2.147112e-05  max resid 5.63177e-05 
    ## ... Similar to previous best
    ## Run 332 stress 0.08985605 
    ## Run 333 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733465  max resid 0.05632058 
    ## Run 334 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238277  max resid 0.04975694 
    ## Run 335 stress 0.0793695 
    ## ... Procrustes: rmse 0.009053981  max resid 0.03252123 
    ## Run 336 stress 0.1088748 
    ## Run 337 stress 0.08985604 
    ## Run 338 stress 0.07927535 
    ## ... Procrustes: rmse 2.463277e-05  max resid 6.179787e-05 
    ## ... Similar to previous best
    ## Run 339 stress 0.0793695 
    ## ... Procrustes: rmse 0.009037461  max resid 0.03241741 
    ## Run 340 stress 0.07951443 
    ## ... Procrustes: rmse 0.01737076  max resid 0.05653543 
    ## Run 341 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238395  max resid 0.04973031 
    ## Run 342 stress 0.07927535 
    ## ... Procrustes: rmse 3.724283e-05  max resid 7.950341e-05 
    ## ... Similar to previous best
    ## Run 343 stress 0.07936954 
    ## ... Procrustes: rmse 0.009028018  max resid 0.03235042 
    ## Run 344 stress 0.09290711 
    ## Run 345 stress 0.104507 
    ## Run 346 stress 0.07969263 
    ## ... Procrustes: rmse 0.01237076  max resid 0.0496143 
    ## Run 347 stress 0.08985604 
    ## Run 348 stress 0.09382081 
    ## Run 349 stress 0.09043778 
    ## Run 350 stress 0.09252225 
    ## Run 351 stress 0.07969258 
    ## ... Procrustes: rmse 0.01236666  max resid 0.04965957 
    ## Run 352 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735866  max resid 0.05645765 
    ## Run 353 stress 0.07951443 
    ## ... Procrustes: rmse 0.01737126  max resid 0.05654379 
    ## Run 354 stress 0.07951447 
    ## ... Procrustes: rmse 0.0173938  max resid 0.05667592 
    ## Run 355 stress 0.09388327 
    ## Run 356 stress 0.09565389 
    ## Run 357 stress 0.07936952 
    ## ... Procrustes: rmse 0.009039874  max resid 0.03242647 
    ## Run 358 stress 0.07927534 
    ## ... Procrustes: rmse 1.059427e-05  max resid 2.713502e-05 
    ## ... Similar to previous best
    ## Run 359 stress 0.09106028 
    ## Run 360 stress 0.0904375 
    ## Run 361 stress 0.0898561 
    ## Run 362 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047058  max resid 0.03247994 
    ## Run 363 stress 0.0902278 
    ## Run 364 stress 0.09388326 
    ## Run 365 stress 0.09072896 
    ## Run 366 stress 0.09474016 
    ## Run 367 stress 0.07927535 
    ## ... Procrustes: rmse 4.425065e-05  max resid 0.0001219979 
    ## ... Similar to previous best
    ## Run 368 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045092  max resid 0.03246448 
    ## Run 369 stress 0.09043743 
    ## Run 370 stress 0.07936953 
    ## ... Procrustes: rmse 0.009028607  max resid 0.03236159 
    ## Run 371 stress 0.09106058 
    ## Run 372 stress 0.0904374 
    ## Run 373 stress 0.0917278 
    ## Run 374 stress 0.07951447 
    ## ... Procrustes: rmse 0.01732007  max resid 0.05618305 
    ## Run 375 stress 0.09544396 
    ## Run 376 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735026  max resid 0.05640456 
    ## Run 377 stress 0.07969262 
    ## ... Procrustes: rmse 0.01237134  max resid 0.04962768 
    ## Run 378 stress 0.08985603 
    ## Run 379 stress 0.09106039 
    ## Run 380 stress 0.09544365 
    ## Run 381 stress 0.07936959 
    ## ... Procrustes: rmse 0.00901706  max resid 0.03229035 
    ## Run 382 stress 0.07951442 
    ## ... Procrustes: rmse 0.01737067  max resid 0.05652285 
    ## Run 383 stress 0.106619 
    ## Run 384 stress 0.07927535 
    ## ... Procrustes: rmse 4.034437e-05  max resid 0.0001124275 
    ## ... Similar to previous best
    ## Run 385 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734927  max resid 0.05642294 
    ## Run 386 stress 0.0796927 
    ## ... Procrustes: rmse 0.01236716  max resid 0.04955584 
    ## Run 387 stress 0.09106031 
    ## Run 388 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734684  max resid 0.05638413 
    ## Run 389 stress 0.09022778 
    ## Run 390 stress 0.07927534 
    ## ... Procrustes: rmse 1.357055e-05  max resid 3.406014e-05 
    ## ... Similar to previous best
    ## Run 391 stress 0.07927534 
    ## ... Procrustes: rmse 6.799886e-06  max resid 1.89492e-05 
    ## ... Similar to previous best
    ## Run 392 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239901  max resid 0.04984758 
    ## Run 393 stress 0.09043771 
    ## Run 394 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237545  max resid 0.04969809 
    ## Run 395 stress 0.07927534 
    ## ... Procrustes: rmse 5.096542e-06  max resid 1.313612e-05 
    ## ... Similar to previous best
    ## Run 396 stress 0.07927534 
    ## ... Procrustes: rmse 1.401404e-05  max resid 3.701658e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.08985614 
    ## Run 398 stress 0.09043774 
    ## Run 399 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237786  max resid 0.04972865 
    ## Run 400 stress 0.07927534 
    ## ... Procrustes: rmse 2.904835e-06  max resid 7.371254e-06 
    ## ... Similar to previous best
    ## Run 401 stress 0.07951443 
    ## ... Procrustes: rmse 0.01737298  max resid 0.05657495 
    ## Run 402 stress 0.07927534 
    ## ... Procrustes: rmse 1.147102e-05  max resid 3.052755e-05 
    ## ... Similar to previous best
    ## Run 403 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736662  max resid 0.05650318 
    ## Run 404 stress 0.07936953 
    ## ... Procrustes: rmse 0.009034685  max resid 0.03239872 
    ## Run 405 stress 0.07969259 
    ## ... Procrustes: rmse 0.01235876  max resid 0.04960361 
    ## Run 406 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734478  max resid 0.05636729 
    ## Run 407 stress 0.09022775 
    ## Run 408 stress 0.07936956 
    ## ... Procrustes: rmse 0.009022591  max resid 0.03232319 
    ## Run 409 stress 0.09290713 
    ## Run 410 stress 0.09072905 
    ## Run 411 stress 0.07927535 
    ## ... Procrustes: rmse 3.486786e-05  max resid 9.513859e-05 
    ## ... Similar to previous best
    ## Run 412 stress 0.09290707 
    ## Run 413 stress 0.09474087 
    ## Run 414 stress 0.09296211 
    ## Run 415 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735417  max resid 0.0564296 
    ## Run 416 stress 0.08985609 
    ## Run 417 stress 0.08985605 
    ## Run 418 stress 0.08985605 
    ## Run 419 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733398  max resid 0.05631172 
    ## Run 420 stress 0.07927534 
    ## ... Procrustes: rmse 1.055115e-05  max resid 2.749458e-05 
    ## ... Similar to previous best
    ## Run 421 stress 0.09296199 
    ## Run 422 stress 0.0793695 
    ## ... Procrustes: rmse 0.009037477  max resid 0.0324163 
    ## Run 423 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237225  max resid 0.04963129 
    ## Run 424 stress 0.09106044 
    ## Run 425 stress 0.090729 
    ## Run 426 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237196  max resid 0.04965542 
    ## Run 427 stress 0.08985605 
    ## Run 428 stress 0.07927535 
    ## ... Procrustes: rmse 4.046698e-05  max resid 9.933352e-05 
    ## ... Similar to previous best
    ## Run 429 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733315  max resid 0.05629184 
    ## Run 430 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736368  max resid 0.0564904 
    ## Run 431 stress 0.08985627 
    ## Run 432 stress 0.09072901 
    ## Run 433 stress 0.09252218 
    ## Run 434 stress 0.09544365 
    ## Run 435 stress 0.0793695 
    ## ... Procrustes: rmse 0.009051954  max resid 0.03250876 
    ## Run 436 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237124  max resid 0.04963109 
    ## Run 437 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238257  max resid 0.04975301 
    ## Run 438 stress 0.08985605 
    ## Run 439 stress 0.07951445 
    ## ... Procrustes: rmse 0.01732379  max resid 0.0562152 
    ## Run 440 stress 0.07969259 
    ## ... Procrustes: rmse 0.0123794  max resid 0.04970369 
    ## Run 441 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237396  max resid 0.04970054 
    ## Run 442 stress 0.08985603 
    ## Run 443 stress 0.07969265 
    ## ... Procrustes: rmse 0.01236893  max resid 0.04958251 
    ## Run 444 stress 0.09544381 
    ## Run 445 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734571  max resid 0.05637257 
    ## Run 446 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043684  max resid 0.03245818 
    ## Run 447 stress 0.07936951 
    ## ... Procrustes: rmse 0.009039406  max resid 0.03242972 
    ## Run 448 stress 0.09290713 
    ## Run 449 stress 0.07951446 
    ## ... Procrustes: rmse 0.01732068  max resid 0.05618904 
    ## Run 450 stress 0.07936951 
    ## ... Procrustes: rmse 0.009034355  max resid 0.03239768 
    ## Run 451 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735403  max resid 0.05642269 
    ## Run 452 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735373  max resid 0.0564273 
    ## Run 453 stress 0.08985621 
    ## Run 454 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734985  max resid 0.05639472 
    ## Run 455 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736567  max resid 0.05650741 
    ## Run 456 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046603  max resid 0.03247804 
    ## Run 457 stress 0.3774314 
    ## Run 458 stress 0.092907 
    ## Run 459 stress 0.07927535 
    ## ... Procrustes: rmse 1.233407e-05  max resid 3.360619e-05 
    ## ... Similar to previous best
    ## Run 460 stress 0.07951447 
    ## ... Procrustes: rmse 0.0173198  max resid 0.05618245 
    ## Run 461 stress 0.09106044 
    ## Run 462 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735583  max resid 0.05644207 
    ## Run 463 stress 0.09072899 
    ## Run 464 stress 0.07936957 
    ## ... Procrustes: rmse 0.009021434  max resid 0.03231724 
    ## Run 465 stress 0.07927535 
    ## ... Procrustes: rmse 1.082268e-05  max resid 2.323133e-05 
    ## ... Similar to previous best
    ## Run 466 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732512  max resid 0.05625544 
    ## Run 467 stress 0.07951443 
    ## ... Procrustes: rmse 0.01736085  max resid 0.05649057 
    ## Run 468 stress 0.07927534 
    ## ... Procrustes: rmse 2.69238e-06  max resid 6.453699e-06 
    ## ... Similar to previous best
    ## Run 469 stress 0.09106033 
    ## Run 470 stress 0.07951442 
    ## ... Procrustes: rmse 0.01737475  max resid 0.0565527 
    ## Run 471 stress 0.09072891 
    ## Run 472 stress 0.380794 
    ## Run 473 stress 0.07969261 
    ## ... Procrustes: rmse 0.01236837  max resid 0.0496004 
    ## Run 474 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733249  max resid 0.05628249 
    ## Run 475 stress 0.07927534 
    ## ... Procrustes: rmse 1.684953e-05  max resid 4.960752e-05 
    ## ... Similar to previous best
    ## Run 476 stress 0.1046273 
    ## Run 477 stress 0.08985614 
    ## Run 478 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044458  max resid 0.03246251 
    ## Run 479 stress 0.07927534 
    ## ... Procrustes: rmse 3.26078e-06  max resid 7.512959e-06 
    ## ... Similar to previous best
    ## Run 480 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735758  max resid 0.05645329 
    ## Run 481 stress 0.09043759 
    ## Run 482 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046782  max resid 0.03247858 
    ## Run 483 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044526  max resid 0.03245667 
    ## Run 484 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734812  max resid 0.05639282 
    ## Run 485 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237009  max resid 0.04962267 
    ## Run 486 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173478  max resid 0.05638406 
    ## Run 487 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237314  max resid 0.04967745 
    ## Run 488 stress 0.07927534 
    ## ... Procrustes: rmse 1.692233e-05  max resid 4.630434e-05 
    ## ... Similar to previous best
    ## Run 489 stress 0.07936952 
    ## ... Procrustes: rmse 0.00903081  max resid 0.03237476 
    ## Run 490 stress 0.08985604 
    ## Run 491 stress 0.09022814 
    ## Run 492 stress 0.07969264 
    ## ... Procrustes: rmse 0.01241136  max resid 0.04992563 
    ## Run 493 stress 0.09388354 
    ## Run 494 stress 0.1041898 
    ## Run 495 stress 0.09344087 
    ## Run 496 stress 0.1039004 
    ## Run 497 stress 0.09474019 
    ## Run 498 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733463  max resid 0.05628793 
    ## Run 499 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237279  max resid 0.04965609 
    ## Run 500 stress 0.07969265 
    ## ... Procrustes: rmse 0.0123742  max resid 0.049599 
    ## *** Best solution repeated 64 times

``` r
# Mixed and stratified lakes
PD_beta_geo_MS_NMDS <- metaMDS(PD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.05051321 
    ## ... Procrustes: rmse 0.0357994  max resid 0.1221649 
    ## Run 2 stress 0.06493696 
    ## Run 3 stress 0.0533758 
    ## Run 4 stress 0.05036796 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002604847  max resid 0.0006578371 
    ## ... Similar to previous best
    ## Run 5 stress 0.06826676 
    ## Run 6 stress 0.05036795 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001986464  max resid 0.0004984389 
    ## ... Similar to previous best
    ## Run 7 stress 0.05036801 
    ## ... Procrustes: rmse 6.263595e-05  max resid 0.0001774486 
    ## ... Similar to previous best
    ## Run 8 stress 0.06027438 
    ## Run 9 stress 0.05051328 
    ## ... Procrustes: rmse 0.03573426  max resid 0.1219449 
    ## Run 10 stress 0.05337596 
    ## Run 11 stress 0.06743857 
    ## Run 12 stress 0.06136678 
    ## Run 13 stress 0.06478521 
    ## Run 14 stress 0.06136695 
    ## Run 15 stress 0.05403599 
    ## Run 16 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 8.034798e-05  max resid 0.0002122258 
    ## ... Similar to previous best
    ## Run 17 stress 0.05602055 
    ## Run 18 stress 0.05829932 
    ## Run 19 stress 0.05051297 
    ## ... Procrustes: rmse 0.03571097  max resid 0.1218904 
    ## Run 20 stress 0.05968629 
    ## Run 21 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001146925  max resid 0.0002852254 
    ## ... Similar to previous best
    ## Run 22 stress 0.06743854 
    ## Run 23 stress 0.05171061 
    ## Run 24 stress 0.06609053 
    ## Run 25 stress 0.05829945 
    ## Run 26 stress 0.050368 
    ## ... Procrustes: rmse 0.0001312841  max resid 0.0003455129 
    ## ... Similar to previous best
    ## Run 27 stress 0.05036796 
    ## ... Procrustes: rmse 8.862773e-05  max resid 0.0002385376 
    ## ... Similar to previous best
    ## Run 28 stress 0.05171059 
    ## Run 29 stress 0.06476227 
    ## Run 30 stress 0.05036797 
    ## ... Procrustes: rmse 0.000101915  max resid 0.0002615058 
    ## ... Similar to previous best
    ## Run 31 stress 0.05829939 
    ## Run 32 stress 0.05171059 
    ## Run 33 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001726116  max resid 0.0004411559 
    ## ... Similar to previous best
    ## Run 34 stress 0.05171058 
    ## Run 35 stress 0.05051318 
    ## ... Procrustes: rmse 0.03569  max resid 0.1218016 
    ## Run 36 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.698645e-05  max resid 2.93728e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.05036794 
    ## ... Procrustes: rmse 7.938862e-05  max resid 0.0001882685 
    ## ... Similar to previous best
    ## Run 38 stress 0.05829935 
    ## Run 39 stress 0.06136682 
    ## Run 40 stress 0.05403573 
    ## Run 41 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001578566  max resid 0.0004006105 
    ## ... Similar to previous best
    ## Run 42 stress 0.06716757 
    ## Run 43 stress 0.05051296 
    ## ... Procrustes: rmse 0.03569997  max resid 0.121961 
    ## Run 44 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001234736  max resid 0.0002963393 
    ## ... Similar to previous best
    ## Run 45 stress 0.05036795 
    ## ... Procrustes: rmse 9.426588e-05  max resid 0.0002300028 
    ## ... Similar to previous best
    ## Run 46 stress 0.05036799 
    ## ... Procrustes: rmse 0.000125865  max resid 0.0003628895 
    ## ... Similar to previous best
    ## Run 47 stress 0.06231298 
    ## Run 48 stress 0.06203054 
    ## Run 49 stress 0.05829931 
    ## Run 50 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001048567  max resid 0.0002707118 
    ## ... Similar to previous best
    ## Run 51 stress 0.05171058 
    ## Run 52 stress 0.05036795 
    ## ... Procrustes: rmse 9.829221e-05  max resid 0.0002693912 
    ## ... Similar to previous best
    ## Run 53 stress 0.05337581 
    ## Run 54 stress 0.050368 
    ## ... Procrustes: rmse 0.0001492275  max resid 0.000370611 
    ## ... Similar to previous best
    ## Run 55 stress 0.05051332 
    ## ... Procrustes: rmse 0.03573158  max resid 0.1220736 
    ## Run 56 stress 0.05829929 
    ## Run 57 stress 0.05928229 
    ## Run 58 stress 0.05036794 
    ## ... Procrustes: rmse 8.776242e-05  max resid 0.0002111512 
    ## ... Similar to previous best
    ## Run 59 stress 0.05171059 
    ## Run 60 stress 0.06309856 
    ## Run 61 stress 0.0636279 
    ## Run 62 stress 0.0582994 
    ## Run 63 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001041345  max resid 0.0002542537 
    ## ... Similar to previous best
    ## Run 64 stress 0.0636278 
    ## Run 65 stress 0.05036795 
    ## ... Procrustes: rmse 8.945656e-05  max resid 0.0002286516 
    ## ... Similar to previous best
    ## Run 66 stress 0.05051337 
    ## ... Procrustes: rmse 0.03569506  max resid 0.1217938 
    ## Run 67 stress 0.05337585 
    ## Run 68 stress 0.05171059 
    ## Run 69 stress 0.05466708 
    ## Run 70 stress 0.05602053 
    ## Run 71 stress 0.0533758 
    ## Run 72 stress 0.07217917 
    ## Run 73 stress 0.05051329 
    ## ... Procrustes: rmse 0.03572417  max resid 0.1220508 
    ## Run 74 stress 0.06309845 
    ## Run 75 stress 0.05051346 
    ## ... Procrustes: rmse 0.03573856  max resid 0.1220995 
    ## Run 76 stress 0.05051332 
    ## ... Procrustes: rmse 0.03574044  max resid 0.1221001 
    ## Run 77 stress 0.06231313 
    ## Run 78 stress 0.05979229 
    ## Run 79 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001493935  max resid 0.0003904213 
    ## ... Similar to previous best
    ## Run 80 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001089341  max resid 0.0002695009 
    ## ... Similar to previous best
    ## Run 81 stress 0.06478547 
    ## Run 82 stress 0.05051313 
    ## ... Procrustes: rmse 0.0357855  max resid 0.1222228 
    ## Run 83 stress 0.06111773 
    ## Run 84 stress 0.0517106 
    ## Run 85 stress 0.05403614 
    ## Run 86 stress 0.05171058 
    ## Run 87 stress 0.06231297 
    ## Run 88 stress 0.05051335 
    ## ... Procrustes: rmse 0.03570059  max resid 0.1219804 
    ## Run 89 stress 0.06362807 
    ## Run 90 stress 0.05051321 
    ## ... Procrustes: rmse 0.03576067  max resid 0.1220104 
    ## Run 91 stress 0.06309847 
    ## Run 92 stress 0.05968631 
    ## Run 93 stress 0.050368 
    ## ... Procrustes: rmse 0.0001537842  max resid 0.00038027 
    ## ... Similar to previous best
    ## Run 94 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001032295  max resid 0.000262653 
    ## ... Similar to previous best
    ## Run 95 stress 0.05466707 
    ## Run 96 stress 0.05171063 
    ## Run 97 stress 0.05968627 
    ## Run 98 stress 0.050368 
    ## ... Procrustes: rmse 0.0001498865  max resid 0.0003723768 
    ## ... Similar to previous best
    ## Run 99 stress 0.05968632 
    ## Run 100 stress 0.05928219 
    ## Run 101 stress 0.05337577 
    ## Run 102 stress 0.05403579 
    ## Run 103 stress 0.05051322 
    ## ... Procrustes: rmse 0.03570176  max resid 0.1219809 
    ## Run 104 stress 0.05051306 
    ## ... Procrustes: rmse 0.03570904  max resid 0.1219945 
    ## Run 105 stress 0.0540356 
    ## Run 106 stress 0.05337577 
    ## Run 107 stress 0.06499236 
    ## Run 108 stress 0.0596864 
    ## Run 109 stress 0.05466708 
    ## Run 110 stress 0.05171058 
    ## Run 111 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001173533  max resid 0.0003253023 
    ## ... Similar to previous best
    ## Run 112 stress 0.05036794 
    ## ... Procrustes: rmse 8.794081e-05  max resid 0.0002179829 
    ## ... Similar to previous best
    ## Run 113 stress 0.05829942 
    ## Run 114 stress 0.05337596 
    ## Run 115 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001177336  max resid 0.0002700121 
    ## ... Similar to previous best
    ## Run 116 stress 0.06789744 
    ## Run 117 stress 0.05036799 
    ## ... Procrustes: rmse 0.000141713  max resid 0.0003549424 
    ## ... Similar to previous best
    ## Run 118 stress 0.05051349 
    ## ... Procrustes: rmse 0.03570574  max resid 0.1218262 
    ## Run 119 stress 0.05051345 
    ## ... Procrustes: rmse 0.03574751  max resid 0.1221254 
    ## Run 120 stress 0.05051285 
    ## ... Procrustes: rmse 0.03569367  max resid 0.121932 
    ## Run 121 stress 0.05829931 
    ## Run 122 stress 0.06309849 
    ## Run 123 stress 0.05829933 
    ## Run 124 stress 0.050513 
    ## ... Procrustes: rmse 0.03569842  max resid 0.1219561 
    ## Run 125 stress 0.05051314 
    ## ... Procrustes: rmse 0.03569646  max resid 0.1219613 
    ## Run 126 stress 0.06476213 
    ## Run 127 stress 0.05051316 
    ## ... Procrustes: rmse 0.03557606  max resid 0.1215951 
    ## Run 128 stress 0.0596863 
    ## Run 129 stress 0.05051352 
    ## ... Procrustes: rmse 0.03574634  max resid 0.1221242 
    ## Run 130 stress 0.05171059 
    ## Run 131 stress 0.05036793 
    ## ... Procrustes: rmse 4.759506e-05  max resid 0.0001120771 
    ## ... Similar to previous best
    ## Run 132 stress 0.06136681 
    ## Run 133 stress 0.06136674 
    ## Run 134 stress 0.0582993 
    ## Run 135 stress 0.06943833 
    ## Run 136 stress 0.05051326 
    ## ... Procrustes: rmse 0.03567812  max resid 0.1217522 
    ## Run 137 stress 0.05602065 
    ## Run 138 stress 0.05968638 
    ## Run 139 stress 0.05036793 
    ## ... Procrustes: rmse 8.311038e-05  max resid 0.0002000574 
    ## ... Similar to previous best
    ## Run 140 stress 0.05051354 
    ## ... Procrustes: rmse 0.03571103  max resid 0.1220197 
    ## Run 141 stress 0.05051312 
    ## ... Procrustes: rmse 0.03572453  max resid 0.1219079 
    ## Run 142 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001889756  max resid 0.0004515902 
    ## ... Similar to previous best
    ## Run 143 stress 0.06704996 
    ## Run 144 stress 0.05051299 
    ## ... Procrustes: rmse 0.03570765  max resid 0.1219863 
    ## Run 145 stress 0.05403611 
    ## Run 146 stress 0.05051291 
    ## ... Procrustes: rmse 0.03569781  max resid 0.1219506 
    ## Run 147 stress 0.06499234 
    ## Run 148 stress 0.05602081 
    ## Run 149 stress 0.05337577 
    ## Run 150 stress 0.05051314 
    ## ... Procrustes: rmse 0.03571316  max resid 0.122011 
    ## Run 151 stress 0.050368 
    ## ... Procrustes: rmse 0.0001425131  max resid 0.0003662918 
    ## ... Similar to previous best
    ## Run 152 stress 0.05171058 
    ## Run 153 stress 0.05602045 
    ## Run 154 stress 0.05403586 
    ## Run 155 stress 0.06203028 
    ## Run 156 stress 0.05928225 
    ## Run 157 stress 0.0505133 
    ## ... Procrustes: rmse 0.03566527  max resid 0.1217109 
    ## Run 158 stress 0.054036 
    ## Run 159 stress 0.05036795 
    ## ... Procrustes: rmse 4.306172e-05  max resid 9.56895e-05 
    ## ... Similar to previous best
    ## Run 160 stress 0.06409264 
    ## Run 161 stress 0.05036794 
    ## ... Procrustes: rmse 8.06799e-05  max resid 0.0001913586 
    ## ... Similar to previous best
    ## Run 162 stress 0.05171058 
    ## Run 163 stress 0.05171058 
    ## Run 164 stress 0.0623129 
    ## Run 165 stress 0.05466702 
    ## Run 166 stress 0.05968637 
    ## Run 167 stress 0.05036799 
    ## ... Procrustes: rmse 0.000146846  max resid 0.0003598551 
    ## ... Similar to previous best
    ## Run 168 stress 0.06136684 
    ## Run 169 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001660091  max resid 0.0004144433 
    ## ... Similar to previous best
    ## Run 170 stress 0.05036794 
    ## ... Procrustes: rmse 8.418943e-05  max resid 0.0002043156 
    ## ... Similar to previous best
    ## Run 171 stress 0.06885084 
    ## Run 172 stress 0.05968627 
    ## Run 173 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001451453  max resid 0.0003928153 
    ## ... Similar to previous best
    ## Run 174 stress 0.05466718 
    ## Run 175 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001358983  max resid 0.0003614149 
    ## ... Similar to previous best
    ## Run 176 stress 0.06136674 
    ## Run 177 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001309796  max resid 0.0003221626 
    ## ... Similar to previous best
    ## Run 178 stress 0.05968641 
    ## Run 179 stress 0.06111766 
    ## Run 180 stress 0.05928224 
    ## Run 181 stress 0.05051341 
    ## ... Procrustes: rmse 0.03573607  max resid 0.1220884 
    ## Run 182 stress 0.05036808 
    ## ... Procrustes: rmse 0.000195996  max resid 0.0005107855 
    ## ... Similar to previous best
    ## Run 183 stress 0.06499238 
    ## Run 184 stress 0.06362786 
    ## Run 185 stress 0.05036797 
    ## ... Procrustes: rmse 9.325528e-05  max resid 0.0002273193 
    ## ... Similar to previous best
    ## Run 186 stress 0.05403606 
    ## Run 187 stress 0.05051334 
    ## ... Procrustes: rmse 0.03575745  max resid 0.121987 
    ## Run 188 stress 0.0582994 
    ## Run 189 stress 0.05036792 
    ## ... Procrustes: rmse 6.183228e-05  max resid 0.0001540953 
    ## ... Similar to previous best
    ## Run 190 stress 0.05403601 
    ## Run 191 stress 0.05968628 
    ## Run 192 stress 0.05928231 
    ## Run 193 stress 0.06362784 
    ## Run 194 stress 0.05602043 
    ## Run 195 stress 0.05466705 
    ## Run 196 stress 0.05036796 
    ## ... Procrustes: rmse 0.000120601  max resid 0.0003003581 
    ## ... Similar to previous best
    ## Run 197 stress 0.05051312 
    ## ... Procrustes: rmse 0.03567291  max resid 0.1217523 
    ## Run 198 stress 0.05171059 
    ## Run 199 stress 0.06730224 
    ## Run 200 stress 0.06499241 
    ## Run 201 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001729457  max resid 0.0004230206 
    ## ... Similar to previous best
    ## Run 202 stress 0.05171065 
    ## Run 203 stress 0.05829929 
    ## Run 204 stress 0.05968629 
    ## Run 205 stress 0.05968638 
    ## Run 206 stress 0.06309848 
    ## Run 207 stress 0.05466702 
    ## Run 208 stress 0.06231308 
    ## Run 209 stress 0.05979228 
    ## Run 210 stress 0.05968629 
    ## Run 211 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001780837  max resid 0.0004504782 
    ## ... Similar to previous best
    ## Run 212 stress 0.05602057 
    ## Run 213 stress 0.06027426 
    ## Run 214 stress 0.05051289 
    ## ... Procrustes: rmse 0.03562442  max resid 0.1217251 
    ## Run 215 stress 0.05466703 
    ## Run 216 stress 0.06362799 
    ## Run 217 stress 0.05403614 
    ## Run 218 stress 0.06362778 
    ## Run 219 stress 0.06231288 
    ## Run 220 stress 0.05829934 
    ## Run 221 stress 0.05036791 
    ## ... Procrustes: rmse 2.947114e-05  max resid 7.713953e-05 
    ## ... Similar to previous best
    ## Run 222 stress 0.05036794 
    ## ... Procrustes: rmse 8.712923e-05  max resid 0.0002100711 
    ## ... Similar to previous best
    ## Run 223 stress 0.05979227 
    ## Run 224 stress 0.05829942 
    ## Run 225 stress 0.05968629 
    ## Run 226 stress 0.05036792 
    ## ... Procrustes: rmse 4.603179e-05  max resid 0.0001136417 
    ## ... Similar to previous best
    ## Run 227 stress 0.05036806 
    ## ... Procrustes: rmse 0.0001738951  max resid 0.0004689528 
    ## ... Similar to previous best
    ## Run 228 stress 0.05829932 
    ## Run 229 stress 0.05036791 
    ## ... Procrustes: rmse 2.662485e-05  max resid 6.18797e-05 
    ## ... Similar to previous best
    ## Run 230 stress 0.06904521 
    ## Run 231 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001126756  max resid 0.0002856579 
    ## ... Similar to previous best
    ## Run 232 stress 0.3592581 
    ## Run 233 stress 0.05036794 
    ## ... Procrustes: rmse 7.276828e-05  max resid 0.0001662277 
    ## ... Similar to previous best
    ## Run 234 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001045549  max resid 0.0002601661 
    ## ... Similar to previous best
    ## Run 235 stress 0.05466702 
    ## Run 236 stress 0.06885112 
    ## Run 237 stress 0.06231299 
    ## Run 238 stress 0.05337577 
    ## Run 239 stress 0.05051278 
    ## ... Procrustes: rmse 0.0356926  max resid 0.1219208 
    ## Run 240 stress 0.05051334 
    ## ... Procrustes: rmse 0.03566187  max resid 0.1216974 
    ## Run 241 stress 0.05928249 
    ## Run 242 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001634753  max resid 0.0004097971 
    ## ... Similar to previous best
    ## Run 243 stress 0.05602053 
    ## Run 244 stress 0.0505129 
    ## ... Procrustes: rmse 0.03560621  max resid 0.1216666 
    ## Run 245 stress 0.050368 
    ## ... Procrustes: rmse 0.0001427278  max resid 0.0003646753 
    ## ... Similar to previous best
    ## Run 246 stress 0.05337577 
    ## Run 247 stress 0.05051329 
    ## ... Procrustes: rmse 0.03570634  max resid 0.1218341 
    ## Run 248 stress 0.06537607 
    ## Run 249 stress 0.05337579 
    ## Run 250 stress 0.06504228 
    ## Run 251 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001285365  max resid 0.0003218122 
    ## ... Similar to previous best
    ## Run 252 stress 0.0560205 
    ## Run 253 stress 0.0533758 
    ## Run 254 stress 0.05171059 
    ## Run 255 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001114501  max resid 0.0002784782 
    ## ... Similar to previous best
    ## Run 256 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001025307  max resid 0.0002418692 
    ## ... Similar to previous best
    ## Run 257 stress 0.06826709 
    ## Run 258 stress 0.05466701 
    ## Run 259 stress 0.05051317 
    ## ... Procrustes: rmse 0.0356959  max resid 0.1218153 
    ## Run 260 stress 0.05171063 
    ## Run 261 stress 0.0582994 
    ## Run 262 stress 0.05466712 
    ## Run 263 stress 0.05968639 
    ## Run 264 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001048214  max resid 0.000256094 
    ## ... Similar to previous best
    ## Run 265 stress 0.06943851 
    ## Run 266 stress 0.05602063 
    ## Run 267 stress 0.05337578 
    ## Run 268 stress 0.050368 
    ## ... Procrustes: rmse 0.0001385462  max resid 0.0003855225 
    ## ... Similar to previous best
    ## Run 269 stress 0.05928239 
    ## Run 270 stress 0.05051305 
    ## ... Procrustes: rmse 0.03567484  max resid 0.1217652 
    ## Run 271 stress 0.05466704 
    ## Run 272 stress 0.05036801 
    ## ... Procrustes: rmse 9.63519e-05  max resid 0.0002235646 
    ## ... Similar to previous best
    ## Run 273 stress 0.05171059 
    ## Run 274 stress 0.05928214 
    ## Run 275 stress 0.0667457 
    ## Run 276 stress 0.05036792 
    ## ... Procrustes: rmse 4.562232e-05  max resid 9.90563e-05 
    ## ... Similar to previous best
    ## Run 277 stress 0.05036792 
    ## ... Procrustes: rmse 4.981473e-05  max resid 0.0001036255 
    ## ... Similar to previous best
    ## Run 278 stress 0.0517106 
    ## Run 279 stress 0.05337584 
    ## Run 280 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001089556  max resid 0.0002736885 
    ## ... Similar to previous best
    ## Run 281 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001083186  max resid 0.0002738616 
    ## ... Similar to previous best
    ## Run 282 stress 0.06885089 
    ## Run 283 stress 0.05171064 
    ## Run 284 stress 0.0560206 
    ## Run 285 stress 0.05171058 
    ## Run 286 stress 0.05968626 
    ## Run 287 stress 0.05051325 
    ## ... Procrustes: rmse 0.0356975  max resid 0.1219673 
    ## Run 288 stress 0.05051283 
    ## ... Procrustes: rmse 0.03569311  max resid 0.1219288 
    ## Run 289 stress 0.06499252 
    ## Run 290 stress 0.05829929 
    ## Run 291 stress 0.05051363 
    ## ... Procrustes: rmse 0.03570647  max resid 0.1220085 
    ## Run 292 stress 0.05829936 
    ## Run 293 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001587422  max resid 0.0003975554 
    ## ... Similar to previous best
    ## Run 294 stress 0.06231308 
    ## Run 295 stress 0.05171059 
    ## Run 296 stress 0.05171058 
    ## Run 297 stress 0.05900148 
    ## Run 298 stress 0.07211128 
    ## Run 299 stress 0.05968631 
    ## Run 300 stress 0.06309836 
    ## Run 301 stress 0.0517106 
    ## Run 302 stress 0.05928256 
    ## Run 303 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001013723  max resid 0.0002478157 
    ## ... Similar to previous best
    ## Run 304 stress 0.05171059 
    ## Run 305 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001622718  max resid 0.0004069083 
    ## ... Similar to previous best
    ## Run 306 stress 0.06499239 
    ## Run 307 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001439305  max resid 0.000358979 
    ## ... Similar to previous best
    ## Run 308 stress 0.05051339 
    ## ... Procrustes: rmse 0.03572911  max resid 0.1220655 
    ## Run 309 stress 0.05036794 
    ## ... Procrustes: rmse 9.362287e-05  max resid 0.000228941 
    ## ... Similar to previous best
    ## Run 310 stress 0.05036799 
    ## ... Procrustes: rmse 0.000141532  max resid 0.0003469816 
    ## ... Similar to previous best
    ## Run 311 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001309828  max resid 0.0003347838 
    ## ... Similar to previous best
    ## Run 312 stress 0.05602059 
    ## Run 313 stress 0.05051316 
    ## ... Procrustes: rmse 0.0356936  max resid 0.1219532 
    ## Run 314 stress 0.05051316 
    ## ... Procrustes: rmse 0.03570237  max resid 0.1218366 
    ## Run 315 stress 0.05829931 
    ## Run 316 stress 0.05602047 
    ## Run 317 stress 0.05051322 
    ## ... Procrustes: rmse 0.03567821  max resid 0.12176 
    ## Run 318 stress 0.06938793 
    ## Run 319 stress 0.05466705 
    ## Run 320 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001474819  max resid 0.0003797642 
    ## ... Similar to previous best
    ## Run 321 stress 0.05036794 
    ## ... Procrustes: rmse 7.440551e-05  max resid 0.000188463 
    ## ... Similar to previous best
    ## Run 322 stress 0.05036795 
    ## ... Procrustes: rmse 8.463377e-05  max resid 0.0002440863 
    ## ... Similar to previous best
    ## Run 323 stress 0.06826689 
    ## Run 324 stress 0.050368 
    ## ... Procrustes: rmse 0.0001343708  max resid 0.0003528814 
    ## ... Similar to previous best
    ## Run 325 stress 0.05051352 
    ## ... Procrustes: rmse 0.03574476  max resid 0.1221187 
    ## Run 326 stress 0.06826694 
    ## Run 327 stress 0.05171058 
    ## Run 328 stress 0.05051304 
    ## ... Procrustes: rmse 0.03567815  max resid 0.1217784 
    ## Run 329 stress 0.07242319 
    ## Run 330 stress 0.05171059 
    ## Run 331 stress 0.0623132 
    ## Run 332 stress 0.05829934 
    ## Run 333 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001021025  max resid 0.0002624654 
    ## ... Similar to previous best
    ## Run 334 stress 0.05051336 
    ## ... Procrustes: rmse 0.0357257  max resid 0.122055 
    ## Run 335 stress 0.05051297 
    ## ... Procrustes: rmse 0.03569086  max resid 0.121824 
    ## Run 336 stress 0.05036792 
    ## ... Procrustes: rmse 4.076724e-05  max resid 9.969835e-05 
    ## ... Similar to previous best
    ## Run 337 stress 0.05602044 
    ## Run 338 stress 0.05337583 
    ## Run 339 stress 0.05051329 
    ## ... Procrustes: rmse 0.03573098  max resid 0.1220711 
    ## Run 340 stress 0.05171058 
    ## Run 341 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001622592  max resid 0.0004004537 
    ## ... Similar to previous best
    ## Run 342 stress 0.05036792 
    ## ... Procrustes: rmse 5.546616e-05  max resid 0.0001314233 
    ## ... Similar to previous best
    ## Run 343 stress 0.05171059 
    ## Run 344 stress 0.06478538 
    ## Run 345 stress 0.05036794 
    ## ... Procrustes: rmse 9.366337e-05  max resid 0.0002307988 
    ## ... Similar to previous best
    ## Run 346 stress 0.05051347 
    ## ... Procrustes: rmse 0.03574284  max resid 0.1219355 
    ## Run 347 stress 0.05036792 
    ## ... Procrustes: rmse 3.277194e-05  max resid 5.872792e-05 
    ## ... Similar to previous best
    ## Run 348 stress 0.06027434 
    ## Run 349 stress 0.06499245 
    ## Run 350 stress 0.06716828 
    ## Run 351 stress 0.05051308 
    ## ... Procrustes: rmse 0.03570612  max resid 0.1219867 
    ## Run 352 stress 0.05051298 
    ## ... Procrustes: rmse 0.03567528  max resid 0.1217763 
    ## Run 353 stress 0.05968636 
    ## Run 354 stress 0.06362782 
    ## Run 355 stress 0.05337589 
    ## Run 356 stress 0.0597923 
    ## Run 357 stress 0.05036793 
    ## ... Procrustes: rmse 7.694243e-05  max resid 0.0001828986 
    ## ... Similar to previous best
    ## Run 358 stress 0.06111766 
    ## Run 359 stress 0.05036804 
    ## ... Procrustes: rmse 0.000180038  max resid 0.0004488146 
    ## ... Similar to previous best
    ## Run 360 stress 0.05466704 
    ## Run 361 stress 0.0546671 
    ## Run 362 stress 0.05036795 
    ## ... Procrustes: rmse 9.35152e-05  max resid 0.0001901469 
    ## ... Similar to previous best
    ## Run 363 stress 0.05036793 
    ## ... Procrustes: rmse 5.661501e-05  max resid 0.0001360646 
    ## ... Similar to previous best
    ## Run 364 stress 0.0517106 
    ## Run 365 stress 0.06038947 
    ## Run 366 stress 0.05036791 
    ## ... Procrustes: rmse 1.305077e-05  max resid 3.080644e-05 
    ## ... Similar to previous best
    ## Run 367 stress 0.0684531 
    ## Run 368 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001158801  max resid 0.000257524 
    ## ... Similar to previous best
    ## Run 369 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001446079  max resid 0.0003661712 
    ## ... Similar to previous best
    ## Run 370 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001090011  max resid 0.0002698917 
    ## ... Similar to previous best
    ## Run 371 stress 0.06231315 
    ## Run 372 stress 0.05171059 
    ## Run 373 stress 0.06231289 
    ## Run 374 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001141235  max resid 0.0002845304 
    ## ... Similar to previous best
    ## Run 375 stress 0.06027428 
    ## Run 376 stress 0.05968632 
    ## Run 377 stress 0.05051321 
    ## ... Procrustes: rmse 0.0356724  max resid 0.1217417 
    ## Run 378 stress 0.0517106 
    ## Run 379 stress 0.06478547 
    ## Run 380 stress 0.05829929 
    ## Run 381 stress 0.06136682 
    ## Run 382 stress 0.05403591 
    ## Run 383 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001051212  max resid 0.0002520583 
    ## ... Similar to previous best
    ## Run 384 stress 0.05337581 
    ## Run 385 stress 0.06231307 
    ## Run 386 stress 0.05171061 
    ## Run 387 stress 0.0517106 
    ## Run 388 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001590964  max resid 0.0003382433 
    ## ... Similar to previous best
    ## Run 389 stress 0.05337591 
    ## Run 390 stress 0.05403594 
    ## Run 391 stress 0.0590016 
    ## Run 392 stress 0.05051311 
    ## ... Procrustes: rmse 0.03568749  max resid 0.1217971 
    ## Run 393 stress 0.05466705 
    ## Run 394 stress 0.3296764 
    ## Run 395 stress 0.06478537 
    ## Run 396 stress 0.05036794 
    ## ... Procrustes: rmse 7.173148e-05  max resid 0.0001761166 
    ## ... Similar to previous best
    ## Run 397 stress 0.05171059 
    ## Run 398 stress 0.06938803 
    ## Run 399 stress 0.06231286 
    ## Run 400 stress 0.05968629 
    ## Run 401 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001429591  max resid 0.0003608147 
    ## ... Similar to previous best
    ## Run 402 stress 0.05829932 
    ## Run 403 stress 0.05036806 
    ## ... Procrustes: rmse 0.0001733695  max resid 0.0004566707 
    ## ... Similar to previous best
    ## Run 404 stress 0.05171058 
    ## Run 405 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001099224  max resid 0.0002753989 
    ## ... Similar to previous best
    ## Run 406 stress 0.05051354 
    ## ... Procrustes: rmse 0.03564079  max resid 0.1216171 
    ## Run 407 stress 0.05036795 
    ## ... Procrustes: rmse 8.79505e-05  max resid 0.0002445334 
    ## ... Similar to previous best
    ## Run 408 stress 0.05036795 
    ## ... Procrustes: rmse 9.2891e-05  max resid 0.000229749 
    ## ... Similar to previous best
    ## Run 409 stress 0.05036794 
    ## ... Procrustes: rmse 6.345868e-05  max resid 0.0001506188 
    ## ... Similar to previous best
    ## Run 410 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001770627  max resid 0.0004484565 
    ## ... Similar to previous best
    ## Run 411 stress 0.05051317 
    ## ... Procrustes: rmse 0.03572313  max resid 0.1218971 
    ## Run 412 stress 0.05602051 
    ## Run 413 stress 0.05036793 
    ## ... Procrustes: rmse 6.991376e-05  max resid 0.0001701937 
    ## ... Similar to previous best
    ## Run 414 stress 0.05403586 
    ## Run 415 stress 0.06478523 
    ## Run 416 stress 0.06111766 
    ## Run 417 stress 0.05036794 
    ## ... Procrustes: rmse 9.210881e-05  max resid 0.0002289422 
    ## ... Similar to previous best
    ## Run 418 stress 0.06730212 
    ## Run 419 stress 0.0582993 
    ## Run 420 stress 0.05051332 
    ## ... Procrustes: rmse 0.03567312  max resid 0.1217324 
    ## Run 421 stress 0.05337578 
    ## Run 422 stress 0.0592822 
    ## Run 423 stress 0.05051319 
    ## ... Procrustes: rmse 0.03573874  max resid 0.1220895 
    ## Run 424 stress 0.0560208 
    ## Run 425 stress 0.05171058 
    ## Run 426 stress 0.05829937 
    ## Run 427 stress 0.06309838 
    ## Run 428 stress 0.07272446 
    ## Run 429 stress 0.05036794 
    ## ... Procrustes: rmse 7.673746e-05  max resid 0.0001860176 
    ## ... Similar to previous best
    ## Run 430 stress 0.06716778 
    ## Run 431 stress 0.05928212 
    ## Run 432 stress 0.05602046 
    ## Run 433 stress 0.05403609 
    ## Run 434 stress 0.05337587 
    ## Run 435 stress 0.05928218 
    ## Run 436 stress 0.06111765 
    ## Run 437 stress 0.0517106 
    ## Run 438 stress 0.05337587 
    ## Run 439 stress 0.05171061 
    ## Run 440 stress 0.0517106 
    ## Run 441 stress 0.05829942 
    ## Run 442 stress 0.05928226 
    ## Run 443 stress 0.0596863 
    ## Run 444 stress 0.06231319 
    ## Run 445 stress 0.05171062 
    ## Run 446 stress 0.06674566 
    ## Run 447 stress 0.0582993 
    ## Run 448 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001115218  max resid 0.0002836929 
    ## ... Similar to previous best
    ## Run 449 stress 0.05968628 
    ## Run 450 stress 0.05171061 
    ## Run 451 stress 0.06111764 
    ## Run 452 stress 0.05968626 
    ## Run 453 stress 0.3550249 
    ## Run 454 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001735976  max resid 0.0004382146 
    ## ... Similar to previous best
    ## Run 455 stress 0.05829929 
    ## Run 456 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001682989  max resid 0.0003492038 
    ## ... Similar to previous best
    ## Run 457 stress 0.06231287 
    ## Run 458 stress 0.05602049 
    ## Run 459 stress 0.05602065 
    ## Run 460 stress 0.05036793 
    ## ... Procrustes: rmse 6.116219e-05  max resid 0.0001126953 
    ## ... Similar to previous best
    ## Run 461 stress 0.05036793 
    ## ... Procrustes: rmse 4.920621e-05  max resid 0.0001162319 
    ## ... Similar to previous best
    ## Run 462 stress 0.05171059 
    ## Run 463 stress 0.05036797 
    ## ... Procrustes: rmse 9.742162e-05  max resid 0.0002336736 
    ## ... Similar to previous best
    ## Run 464 stress 0.0582993 
    ## Run 465 stress 0.05829946 
    ## Run 466 stress 0.05051314 
    ## ... Procrustes: rmse 0.03574848  max resid 0.121979 
    ## Run 467 stress 0.05171058 
    ## Run 468 stress 0.05051306 
    ## ... Procrustes: rmse 0.03569984  max resid 0.1218403 
    ## Run 469 stress 0.05337578 
    ## Run 470 stress 0.0590015 
    ## Run 471 stress 0.05051293 
    ## ... Procrustes: rmse 0.03565854  max resid 0.1218329 
    ## Run 472 stress 0.050368 
    ## ... Procrustes: rmse 0.0001388809  max resid 0.0003605077 
    ## ... Similar to previous best
    ## Run 473 stress 0.05466708 
    ## Run 474 stress 0.06231293 
    ## Run 475 stress 0.05051332 
    ## ... Procrustes: rmse 0.03567475  max resid 0.1217425 
    ## Run 476 stress 0.05036793 
    ## ... Procrustes: rmse 8.267132e-05  max resid 0.0002016719 
    ## ... Similar to previous best
    ## Run 477 stress 0.05171058 
    ## Run 478 stress 0.05051323 
    ## ... Procrustes: rmse 0.0357255  max resid 0.1220522 
    ## Run 479 stress 0.0505133 
    ## ... Procrustes: rmse 0.03577866  max resid 0.1220572 
    ## Run 480 stress 0.0505132 
    ## ... Procrustes: rmse 0.03569243  max resid 0.1218022 
    ## Run 481 stress 0.05968629 
    ## Run 482 stress 0.05051325 
    ## ... Procrustes: rmse 0.03572515  max resid 0.1220517 
    ## Run 483 stress 0.05051316 
    ## ... Procrustes: rmse 0.0356953  max resid 0.121958 
    ## Run 484 stress 0.05051302 
    ## ... Procrustes: rmse 0.03570038  max resid 0.1218469 
    ## Run 485 stress 0.05900148 
    ## Run 486 stress 0.05602055 
    ## Run 487 stress 0.05051328 
    ## ... Procrustes: rmse 0.03569881  max resid 0.1219746 
    ## Run 488 stress 0.05829931 
    ## Run 489 stress 0.06478518 
    ## Run 490 stress 0.05337577 
    ## Run 491 stress 0.05829932 
    ## Run 492 stress 0.05602053 
    ## Run 493 stress 0.06027442 
    ## Run 494 stress 0.05337596 
    ## Run 495 stress 0.07026147 
    ## Run 496 stress 0.05036793 
    ## ... Procrustes: rmse 5.893925e-05  max resid 0.0001406195 
    ## ... Similar to previous best
    ## Run 497 stress 0.06231311 
    ## Run 498 stress 0.05928254 
    ## Run 499 stress 0.05466703 
    ## Run 500 stress 0.05171058 
    ## *** Best solution repeated 107 times

``` r
# Ocean sites and mixed lakes
PD_beta_geo_OM_NMDS <- metaMDS(PD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1016168 
    ## Run 1 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0005560119  max resid 0.001554321 
    ## ... Similar to previous best
    ## Run 2 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001582932  max resid 0.0004448868 
    ## ... Similar to previous best
    ## Run 3 stress 0.1356145 
    ## Run 4 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 4.712755e-06  max resid 1.315949e-05 
    ## ... Similar to previous best
    ## Run 5 stress 0.1016164 
    ## ... Procrustes: rmse 4.84425e-05  max resid 0.0001365986 
    ## ... Similar to previous best
    ## Run 6 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004705309  max resid 0.001315745 
    ## ... Similar to previous best
    ## Run 7 stress 0.1356147 
    ## Run 8 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001046295  max resid 0.0002951332 
    ## ... Similar to previous best
    ## Run 9 stress 0.1016164 
    ## ... Procrustes: rmse 0.000109236  max resid 0.0003082009 
    ## ... Similar to previous best
    ## Run 10 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004094124  max resid 0.001141603 
    ## ... Similar to previous best
    ## Run 11 stress 0.101617 
    ## ... Procrustes: rmse 0.0006934422  max resid 0.001937071 
    ## ... Similar to previous best
    ## Run 12 stress 0.1016164 
    ## ... Procrustes: rmse 9.170659e-07  max resid 1.937471e-06 
    ## ... Similar to previous best
    ## Run 13 stress 0.1341868 
    ## Run 14 stress 0.1016164 
    ## ... Procrustes: rmse 3.496443e-05  max resid 9.801594e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.1016164 
    ## ... Procrustes: rmse 2.378269e-05  max resid 6.647766e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003392176  max resid 0.0009497301 
    ## ... Similar to previous best
    ## Run 17 stress 0.1016167 
    ## ... Procrustes: rmse 0.000415595  max resid 0.001163382 
    ## ... Similar to previous best
    ## Run 18 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003718307  max resid 0.001039178 
    ## ... Similar to previous best
    ## Run 19 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004232918  max resid 0.001184553 
    ## ... Similar to previous best
    ## Run 20 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001565313  max resid 0.0004400155 
    ## ... Similar to previous best
    ## Run 21 stress 0.135615 
    ## Run 22 stress 0.1016164 
    ## ... Procrustes: rmse 4.361417e-06  max resid 1.210612e-05 
    ## ... Similar to previous best
    ## Run 23 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001832879  max resid 0.0005151551 
    ## ... Similar to previous best
    ## Run 24 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003289865  max resid 0.0009194272 
    ## ... Similar to previous best
    ## Run 25 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005193075  max resid 0.001453517 
    ## ... Similar to previous best
    ## Run 26 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004351924  max resid 0.00121734 
    ## ... Similar to previous best
    ## Run 27 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002124019  max resid 0.0005942092 
    ## ... Similar to previous best
    ## Run 28 stress 0.1016164 
    ## ... Procrustes: rmse 3.913141e-05  max resid 0.0001095766 
    ## ... Similar to previous best
    ## Run 29 stress 0.1519033 
    ## Run 30 stress 0.1016165 
    ## ... Procrustes: rmse 0.000259707  max resid 0.0007281561 
    ## ... Similar to previous best
    ## Run 31 stress 0.1016164 
    ## ... Procrustes: rmse 2.858678e-05  max resid 7.631951e-05 
    ## ... Similar to previous best
    ## Run 32 stress 0.1341868 
    ## Run 33 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005873957  max resid 0.001635956 
    ## ... Similar to previous best
    ## Run 34 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004392998  max resid 0.001228811 
    ## ... Similar to previous best
    ## Run 35 stress 0.1341868 
    ## Run 36 stress 0.1356148 
    ## Run 37 stress 0.1016164 
    ## ... Procrustes: rmse 2.02217e-05  max resid 5.710252e-05 
    ## ... Similar to previous best
    ## Run 38 stress 0.1016165 
    ## ... Procrustes: rmse 9.292979e-05  max resid 0.0002335581 
    ## ... Similar to previous best
    ## Run 39 stress 0.1341868 
    ## Run 40 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004159286  max resid 0.001164244 
    ## ... Similar to previous best
    ## Run 41 stress 0.1341868 
    ## Run 42 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004701795  max resid 0.001316251 
    ## ... Similar to previous best
    ## Run 43 stress 0.1016166 
    ## ... Procrustes: rmse 0.000329263  max resid 0.0009224331 
    ## ... Similar to previous best
    ## Run 44 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003809034  max resid 0.001065661 
    ## ... Similar to previous best
    ## Run 45 stress 0.1016164 
    ## ... Procrustes: rmse 2.029677e-05  max resid 4.118391e-05 
    ## ... Similar to previous best
    ## Run 46 stress 0.1356148 
    ## Run 47 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002767541  max resid 0.0007744567 
    ## ... Similar to previous best
    ## Run 48 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002412675  max resid 0.0006772944 
    ## ... Similar to previous best
    ## Run 49 stress 0.1530426 
    ## Run 50 stress 0.1530427 
    ## Run 51 stress 0.1016164 
    ## ... Procrustes: rmse 4.304514e-05  max resid 0.0001178136 
    ## ... Similar to previous best
    ## Run 52 stress 0.1016164 
    ## ... Procrustes: rmse 9.874497e-06  max resid 2.506475e-05 
    ## ... Similar to previous best
    ## Run 53 stress 0.3364383 
    ## Run 54 stress 0.1519033 
    ## Run 55 stress 0.1016164 
    ## ... Procrustes: rmse 5.870757e-05  max resid 0.0001653671 
    ## ... Similar to previous best
    ## Run 56 stress 0.1016164 
    ## ... Procrustes: rmse 3.275583e-05  max resid 9.239337e-05 
    ## ... Similar to previous best
    ## Run 57 stress 0.1356148 
    ## Run 58 stress 0.1016164 
    ## ... Procrustes: rmse 4.183782e-05  max resid 0.0001179243 
    ## ... Similar to previous best
    ## Run 59 stress 0.1356149 
    ## Run 60 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005269823  max resid 0.001472759 
    ## ... Similar to previous best
    ## Run 61 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001040354  max resid 0.0002933189 
    ## ... Similar to previous best
    ## Run 62 stress 0.1341868 
    ## Run 63 stress 0.1341868 
    ## Run 64 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004806724  max resid 0.001343365 
    ## ... Similar to previous best
    ## Run 65 stress 0.1016164 
    ## ... Procrustes: rmse 9.965147e-06  max resid 2.833117e-05 
    ## ... Similar to previous best
    ## Run 66 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004535169  max resid 0.00126904 
    ## ... Similar to previous best
    ## Run 67 stress 0.1016164 
    ## ... Procrustes: rmse 2.526679e-05  max resid 7.12737e-05 
    ## ... Similar to previous best
    ## Run 68 stress 0.1016164 
    ## ... Procrustes: rmse 2.317338e-06  max resid 6.212824e-06 
    ## ... Similar to previous best
    ## Run 69 stress 0.1016164 
    ## ... Procrustes: rmse 7.714166e-05  max resid 0.000217722 
    ## ... Similar to previous best
    ## Run 70 stress 0.1016164 
    ## ... Procrustes: rmse 7.834935e-05  max resid 0.0002204359 
    ## ... Similar to previous best
    ## Run 71 stress 0.101617 
    ## ... Procrustes: rmse 0.0006438947  max resid 0.001799312 
    ## ... Similar to previous best
    ## Run 72 stress 0.1016164 
    ## ... Procrustes: rmse 7.435182e-06  max resid 2.094326e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.1356147 
    ## Run 74 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005184604  max resid 0.001450399 
    ## ... Similar to previous best
    ## Run 75 stress 0.1016164 
    ## ... Procrustes: rmse 1.330922e-05  max resid 3.727698e-05 
    ## ... Similar to previous best
    ## Run 76 stress 0.1016164 
    ## ... Procrustes: rmse 2.601951e-05  max resid 7.342989e-05 
    ## ... Similar to previous best
    ## Run 77 stress 0.1356143 
    ## Run 78 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001843484  max resid 0.000517005 
    ## ... Similar to previous best
    ## Run 79 stress 0.1519033 
    ## Run 80 stress 0.1016164 
    ## ... Procrustes: rmse 1.356828e-05  max resid 3.819503e-05 
    ## ... Similar to previous best
    ## Run 81 stress 0.1016164 
    ## ... Procrustes: rmse 1.390242e-05  max resid 3.627972e-05 
    ## ... Similar to previous best
    ## Run 82 stress 0.1341868 
    ## Run 83 stress 0.1016164 
    ## ... Procrustes: rmse 2.657084e-05  max resid 6.701822e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.1016164 
    ## ... Procrustes: rmse 6.539855e-06  max resid 1.845598e-05 
    ## ... Similar to previous best
    ## Run 85 stress 0.1016164 
    ## ... Procrustes: rmse 4.092947e-05  max resid 9.900069e-05 
    ## ... Similar to previous best
    ## Run 86 stress 0.1016164 
    ## ... Procrustes: rmse 3.707115e-05  max resid 0.0001041852 
    ## ... Similar to previous best
    ## Run 87 stress 0.1016164 
    ## ... Procrustes: rmse 3.519773e-06  max resid 9.856168e-06 
    ## ... Similar to previous best
    ## Run 88 stress 0.1016164 
    ## ... Procrustes: rmse 5.780893e-05  max resid 0.0001632897 
    ## ... Similar to previous best
    ## Run 89 stress 0.1341868 
    ## Run 90 stress 0.1341868 
    ## Run 91 stress 0.1356145 
    ## Run 92 stress 0.1341868 
    ## Run 93 stress 0.1530426 
    ## Run 94 stress 0.1016164 
    ## ... Procrustes: rmse 6.094352e-05  max resid 0.0001697577 
    ## ... Similar to previous best
    ## Run 95 stress 0.1341868 
    ## Run 96 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002447638  max resid 0.00068604 
    ## ... Similar to previous best
    ## Run 97 stress 0.1341868 
    ## Run 98 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002804097  max resid 0.0007865753 
    ## ... Similar to previous best
    ## Run 99 stress 0.1356145 
    ## Run 100 stress 0.1016164 
    ## ... Procrustes: rmse 5.364352e-05  max resid 0.0001513523 
    ## ... Similar to previous best
    ## Run 101 stress 0.1341868 
    ## Run 102 stress 0.1341868 
    ## Run 103 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006000721  max resid 0.001674767 
    ## ... Similar to previous best
    ## Run 104 stress 0.1016164 
    ## ... Procrustes: rmse 5.404005e-05  max resid 0.0001522079 
    ## ... Similar to previous best
    ## Run 105 stress 0.1341868 
    ## Run 106 stress 0.1016164 
    ## ... Procrustes: rmse 7.544008e-05  max resid 0.000209645 
    ## ... Similar to previous best
    ## Run 107 stress 0.1519033 
    ## Run 108 stress 0.1341868 
    ## Run 109 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003728649  max resid 0.001044589 
    ## ... Similar to previous best
    ## Run 110 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004860747  max resid 0.00135889 
    ## ... Similar to previous best
    ## Run 111 stress 0.1016164 
    ## ... Procrustes: rmse 5.775006e-05  max resid 0.0001631602 
    ## ... Similar to previous best
    ## Run 112 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003644566  max resid 0.001021438 
    ## ... Similar to previous best
    ## Run 113 stress 0.1016164 
    ## ... Procrustes: rmse 5.796366e-05  max resid 0.0001630012 
    ## ... Similar to previous best
    ## Run 114 stress 0.1016164 
    ## ... Procrustes: rmse 5.811416e-05  max resid 0.0001637078 
    ## ... Similar to previous best
    ## Run 115 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003852642  max resid 0.001078862 
    ## ... Similar to previous best
    ## Run 116 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005538681  max resid 0.00154804 
    ## ... Similar to previous best
    ## Run 117 stress 0.1016164 
    ## ... Procrustes: rmse 1.823249e-05  max resid 5.143922e-05 
    ## ... Similar to previous best
    ## Run 118 stress 0.1016164 
    ## ... Procrustes: rmse 1.46373e-05  max resid 4.142248e-05 
    ## ... Similar to previous best
    ## Run 119 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005673932  max resid 0.001586179 
    ## ... Similar to previous best
    ## Run 120 stress 0.1016164 
    ## ... Procrustes: rmse 4.144145e-05  max resid 0.0001143851 
    ## ... Similar to previous best
    ## Run 121 stress 0.1341868 
    ## Run 122 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004215973  max resid 0.00117943 
    ## ... Similar to previous best
    ## Run 123 stress 0.1016164 
    ## ... Procrustes: rmse 7.414036e-06  max resid 2.065743e-05 
    ## ... Similar to previous best
    ## Run 124 stress 0.1356152 
    ## Run 125 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005245403  max resid 0.001465831 
    ## ... Similar to previous best
    ## Run 126 stress 0.1016164 
    ## ... Procrustes: rmse 2.340146e-05  max resid 6.603489e-05 
    ## ... Similar to previous best
    ## Run 127 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004939551  max resid 0.001381446 
    ## ... Similar to previous best
    ## Run 128 stress 0.1537555 
    ## Run 129 stress 0.1016164 
    ## ... Procrustes: rmse 3.585982e-05  max resid 0.0001009505 
    ## ... Similar to previous best
    ## Run 130 stress 0.1356146 
    ## Run 131 stress 0.1016164 
    ## ... Procrustes: rmse 4.291354e-05  max resid 0.0001211439 
    ## ... Similar to previous best
    ## Run 132 stress 0.1016164 
    ## ... Procrustes: rmse 2.914871e-05  max resid 8.023309e-05 
    ## ... Similar to previous best
    ## Run 133 stress 0.1341868 
    ## Run 134 stress 0.1016164 
    ## ... Procrustes: rmse 8.532702e-05  max resid 0.0002399775 
    ## ... Similar to previous best
    ## Run 135 stress 0.1356147 
    ## Run 136 stress 0.1537555 
    ## Run 137 stress 0.1530425 
    ## Run 138 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002158782  max resid 0.0006040811 
    ## ... Similar to previous best
    ## Run 139 stress 0.1016164 
    ## ... Procrustes: rmse 2.039855e-05  max resid 5.753431e-05 
    ## ... Similar to previous best
    ## Run 140 stress 0.1341868 
    ## Run 141 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004483133  max resid 0.001253677 
    ## ... Similar to previous best
    ## Run 142 stress 0.1016164 
    ## ... Procrustes: rmse 2.709159e-05  max resid 7.65489e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.1016164 
    ## ... Procrustes: rmse 1.35017e-05  max resid 3.792281e-05 
    ## ... Similar to previous best
    ## Run 144 stress 0.101617 
    ## ... Procrustes: rmse 0.000692316  max resid 0.001934036 
    ## ... Similar to previous best
    ## Run 145 stress 0.101617 
    ## ... Procrustes: rmse 0.0006114065  max resid 0.00170482 
    ## ... Similar to previous best
    ## Run 146 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005318351  max resid 0.001485982 
    ## ... Similar to previous best
    ## Run 147 stress 0.1341868 
    ## Run 148 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001012977  max resid 0.0002857278 
    ## ... Similar to previous best
    ## Run 149 stress 0.1519033 
    ## Run 150 stress 0.1356148 
    ## Run 151 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003108896  max resid 0.0008712112 
    ## ... Similar to previous best
    ## Run 152 stress 0.135615 
    ## Run 153 stress 0.1016164 
    ## ... Procrustes: rmse 9.496397e-05  max resid 0.0002656621 
    ## ... Similar to previous best
    ## Run 154 stress 0.1341868 
    ## Run 155 stress 0.1341868 
    ## Run 156 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001029147  max resid 0.0002890405 
    ## ... Similar to previous best
    ## Run 157 stress 0.1519033 
    ## Run 158 stress 0.1356149 
    ## Run 159 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005398601  max resid 0.001508515 
    ## ... Similar to previous best
    ## Run 160 stress 0.1356151 
    ## Run 161 stress 0.349572 
    ## Run 162 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007200079  max resid 0.002011923 
    ## ... Similar to previous best
    ## Run 163 stress 0.1530422 
    ## Run 164 stress 0.1341868 
    ## Run 165 stress 0.1016171 
    ## ... Procrustes: rmse 0.0006893233  max resid 0.001923016 
    ## ... Similar to previous best
    ## Run 166 stress 0.1016164 
    ## ... Procrustes: rmse 1.243019e-05  max resid 3.319575e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.1016164 
    ## ... Procrustes: rmse 7.216357e-07  max resid 1.138792e-06 
    ## ... Similar to previous best
    ## Run 168 stress 0.1016164 
    ## ... Procrustes: rmse 4.860156e-05  max resid 0.0001360759 
    ## ... Similar to previous best
    ## Run 169 stress 0.1356148 
    ## Run 170 stress 0.1356148 
    ## Run 171 stress 0.1519033 
    ## Run 172 stress 0.1341868 
    ## Run 173 stress 0.1016164 
    ## ... Procrustes: rmse 1.111034e-05  max resid 3.133157e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.1356144 
    ## Run 175 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001273212  max resid 0.0003555351 
    ## ... Similar to previous best
    ## Run 176 stress 0.1016164 
    ## ... Procrustes: rmse 2.354667e-05  max resid 6.58975e-05 
    ## ... Similar to previous best
    ## Run 177 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002560981  max resid 0.0007179844 
    ## ... Similar to previous best
    ## Run 178 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001369224  max resid 0.0003851224 
    ## ... Similar to previous best
    ## Run 179 stress 0.1016164 
    ## ... Procrustes: rmse 2.45445e-05  max resid 6.795016e-05 
    ## ... Similar to previous best
    ## Run 180 stress 0.1016164 
    ## ... Procrustes: rmse 1.956206e-05  max resid 5.447216e-05 
    ## ... Similar to previous best
    ## Run 181 stress 0.1341868 
    ## Run 182 stress 0.1341868 
    ## Run 183 stress 0.1016164 
    ## ... Procrustes: rmse 1.223687e-05  max resid 3.429519e-05 
    ## ... Similar to previous best
    ## Run 184 stress 0.1341868 
    ## Run 185 stress 0.1356148 
    ## Run 186 stress 0.1341868 
    ## Run 187 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001135701  max resid 0.0003167424 
    ## ... Similar to previous best
    ## Run 188 stress 0.1016164 
    ## ... Procrustes: rmse 3.935424e-06  max resid 1.113337e-05 
    ## ... Similar to previous best
    ## Run 189 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001048935  max resid 0.0002951522 
    ## ... Similar to previous best
    ## Run 190 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002213069  max resid 0.0006040976 
    ## ... Similar to previous best
    ## Run 191 stress 0.1016164 
    ## ... Procrustes: rmse 6.52064e-06  max resid 1.803178e-05 
    ## ... Similar to previous best
    ## Run 192 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005181447  max resid 0.001447349 
    ## ... Similar to previous best
    ## Run 193 stress 0.1356148 
    ## Run 194 stress 0.1016164 
    ## ... Procrustes: rmse 3.770002e-05  max resid 0.0001064206 
    ## ... Similar to previous best
    ## Run 195 stress 0.1016164 
    ## ... Procrustes: rmse 4.495652e-05  max resid 0.0001260208 
    ## ... Similar to previous best
    ## Run 196 stress 0.1016164 
    ## ... Procrustes: rmse 5.978793e-05  max resid 0.0001683997 
    ## ... Similar to previous best
    ## Run 197 stress 0.1537555 
    ## Run 198 stress 0.1341868 
    ## Run 199 stress 0.1537555 
    ## Run 200 stress 0.1016164 
    ## ... Procrustes: rmse 6.205997e-06  max resid 1.753094e-05 
    ## ... Similar to previous best
    ## Run 201 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004206139  max resid 0.001177852 
    ## ... Similar to previous best
    ## Run 202 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003288624  max resid 0.0009212745 
    ## ... Similar to previous best
    ## Run 203 stress 0.3495558 
    ## Run 204 stress 0.1016164 
    ## ... Procrustes: rmse 1.125367e-05  max resid 3.13679e-05 
    ## ... Similar to previous best
    ## Run 205 stress 0.1016164 
    ## ... Procrustes: rmse 1.161412e-05  max resid 3.276885e-05 
    ## ... Similar to previous best
    ## Run 206 stress 0.1016164 
    ## ... Procrustes: rmse 1.683075e-05  max resid 3.318071e-05 
    ## ... Similar to previous best
    ## Run 207 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005016218  max resid 0.001402956 
    ## ... Similar to previous best
    ## Run 208 stress 0.1519033 
    ## Run 209 stress 0.1016164 
    ## ... Procrustes: rmse 9.335981e-05  max resid 0.0002627446 
    ## ... Similar to previous best
    ## Run 210 stress 0.1016164 
    ## ... Procrustes: rmse 6.658664e-05  max resid 0.0001878176 
    ## ... Similar to previous best
    ## Run 211 stress 0.1341868 
    ## Run 212 stress 0.1530427 
    ## Run 213 stress 0.1356146 
    ## Run 214 stress 0.1016164 
    ## ... Procrustes: rmse 6.577938e-05  max resid 0.0001857799 
    ## ... Similar to previous best
    ## Run 215 stress 0.1016164 
    ## ... Procrustes: rmse 1.251598e-05  max resid 3.520169e-05 
    ## ... Similar to previous best
    ## Run 216 stress 0.1016164 
    ## ... Procrustes: rmse 4.745041e-05  max resid 0.0001337567 
    ## ... Similar to previous best
    ## Run 217 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001934708  max resid 0.000543405 
    ## ... Similar to previous best
    ## Run 218 stress 0.1016164 
    ## ... Procrustes: rmse 4.824335e-05  max resid 0.000135831 
    ## ... Similar to previous best
    ## Run 219 stress 0.1341868 
    ## Run 220 stress 0.1016164 
    ## ... Procrustes: rmse 4.13841e-05  max resid 0.0001164388 
    ## ... Similar to previous best
    ## Run 221 stress 0.1519033 
    ## Run 222 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002812807  max resid 0.0007882694 
    ## ... Similar to previous best
    ## Run 223 stress 0.1519033 
    ## Run 224 stress 0.1341868 
    ## Run 225 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002052952  max resid 0.0005765157 
    ## ... Similar to previous best
    ## Run 226 stress 0.1016164 
    ## ... Procrustes: rmse 5.27127e-05  max resid 0.0001450885 
    ## ... Similar to previous best
    ## Run 227 stress 0.1016164 
    ## ... Procrustes: rmse 2.680015e-05  max resid 5.846587e-05 
    ## ... Similar to previous best
    ## Run 228 stress 0.1016164 
    ## ... Procrustes: rmse 4.269686e-05  max resid 0.0001183465 
    ## ... Similar to previous best
    ## Run 229 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005024594  max resid 0.001402723 
    ## ... Similar to previous best
    ## Run 230 stress 0.3242709 
    ## Run 231 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006246753  max resid 0.001745487 
    ## ... Similar to previous best
    ## Run 232 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003091663  max resid 0.0008656902 
    ## ... Similar to previous best
    ## Run 233 stress 0.1341868 
    ## Run 234 stress 0.1341868 
    ## Run 235 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003648237  max resid 0.001022071 
    ## ... Similar to previous best
    ## Run 236 stress 0.1341868 
    ## Run 237 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003141294  max resid 0.0008800379 
    ## ... Similar to previous best
    ## Run 238 stress 0.1016164 
    ## ... Procrustes: rmse 2.971994e-05  max resid 7.722782e-05 
    ## ... Similar to previous best
    ## Run 239 stress 0.1016164 
    ## ... Procrustes: rmse 2.298667e-06  max resid 6.290975e-06 
    ## ... Similar to previous best
    ## Run 240 stress 0.1341868 
    ## Run 241 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005772222  max resid 0.00161132 
    ## ... Similar to previous best
    ## Run 242 stress 0.1016164 
    ## ... Procrustes: rmse 7.423431e-05  max resid 0.0002095112 
    ## ... Similar to previous best
    ## Run 243 stress 0.1356151 
    ## Run 244 stress 0.1016164 
    ## ... Procrustes: rmse 1.500153e-05  max resid 3.709157e-05 
    ## ... Similar to previous best
    ## Run 245 stress 0.1016164 
    ## ... Procrustes: rmse 2.51001e-06  max resid 6.795371e-06 
    ## ... Similar to previous best
    ## Run 246 stress 0.1016164 
    ## ... Procrustes: rmse 6.465221e-05  max resid 0.0001790223 
    ## ... Similar to previous best
    ## Run 247 stress 0.1016164 
    ## ... Procrustes: rmse 3.078167e-05  max resid 8.680114e-05 
    ## ... Similar to previous best
    ## Run 248 stress 0.1341868 
    ## Run 249 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004079885  max resid 0.001141554 
    ## ... Similar to previous best
    ## Run 250 stress 0.1016164 
    ## ... Procrustes: rmse 1.986331e-05  max resid 4.705549e-05 
    ## ... Similar to previous best
    ## Run 251 stress 0.3411839 
    ## Run 252 stress 0.1356143 
    ## Run 253 stress 0.1016164 
    ## ... Procrustes: rmse 4.198196e-05  max resid 0.0001170728 
    ## ... Similar to previous best
    ## Run 254 stress 0.1016164 
    ## ... Procrustes: rmse 1.208859e-06  max resid 2.573961e-06 
    ## ... Similar to previous best
    ## Run 255 stress 0.1016164 
    ## ... Procrustes: rmse 7.067368e-05  max resid 0.0001980635 
    ## ... Similar to previous best
    ## Run 256 stress 0.1341868 
    ## Run 257 stress 0.1341868 
    ## Run 258 stress 0.1016164 
    ## ... Procrustes: rmse 3.398138e-05  max resid 9.585038e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.1016164 
    ## ... Procrustes: rmse 5.635458e-05  max resid 0.0001585012 
    ## ... Similar to previous best
    ## Run 260 stress 0.1016164 
    ## ... Procrustes: rmse 9.183958e-05  max resid 0.000258612 
    ## ... Similar to previous best
    ## Run 261 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005333405  max resid 0.001492132 
    ## ... Similar to previous best
    ## Run 262 stress 0.1356145 
    ## Run 263 stress 0.1016164 
    ## ... Procrustes: rmse 3.664881e-05  max resid 0.0001031401 
    ## ... Similar to previous best
    ## Run 264 stress 0.1341868 
    ## Run 265 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001074832  max resid 0.0003010191 
    ## ... Similar to previous best
    ## Run 266 stress 0.1341868 
    ## Run 267 stress 0.1341868 
    ## Run 268 stress 0.1341868 
    ## Run 269 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005762239  max resid 0.001610437 
    ## ... Similar to previous best
    ## Run 270 stress 0.1016164 
    ## ... Procrustes: rmse 1.284032e-05  max resid 2.823697e-05 
    ## ... Similar to previous best
    ## Run 271 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002689773  max resid 0.0007540116 
    ## ... Similar to previous best
    ## Run 272 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007347275  max resid 0.002051891 
    ## ... Similar to previous best
    ## Run 273 stress 0.1016164 
    ## ... Procrustes: rmse 4.447915e-06  max resid 1.012741e-05 
    ## ... Similar to previous best
    ## Run 274 stress 0.1341868 
    ## Run 275 stress 0.1519033 
    ## Run 276 stress 0.1341868 
    ## Run 277 stress 0.1016164 
    ## ... Procrustes: rmse 1.298934e-05  max resid 3.667397e-05 
    ## ... Similar to previous best
    ## Run 278 stress 0.1016164 
    ## ... Procrustes: rmse 7.272154e-05  max resid 0.0002051747 
    ## ... Similar to previous best
    ## Run 279 stress 0.1537555 
    ## Run 280 stress 0.1530426 
    ## Run 281 stress 0.1016164 
    ## ... Procrustes: rmse 4.587968e-05  max resid 0.0001296093 
    ## ... Similar to previous best
    ## Run 282 stress 0.1016165 
    ## ... Procrustes: rmse 0.000158502  max resid 0.0004410803 
    ## ... Similar to previous best
    ## Run 283 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005445025  max resid 0.001522221 
    ## ... Similar to previous best
    ## Run 284 stress 0.1016166 
    ## ... Procrustes: rmse 0.00034191  max resid 0.0009571488 
    ## ... Similar to previous best
    ## Run 285 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003313761  max resid 0.000928788 
    ## ... Similar to previous best
    ## Run 286 stress 0.1016164 
    ## ... Procrustes: rmse 1.109452e-05  max resid 3.059751e-05 
    ## ... Similar to previous best
    ## Run 287 stress 0.1341868 
    ## Run 288 stress 0.1341868 
    ## Run 289 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003090683  max resid 0.0008658479 
    ## ... Similar to previous best
    ## Run 290 stress 0.1016164 
    ## ... Procrustes: rmse 5.905438e-05  max resid 0.0001648851 
    ## ... Similar to previous best
    ## Run 291 stress 0.1356148 
    ## Run 292 stress 0.1016164 
    ## ... Procrustes: rmse 7.244591e-05  max resid 0.0002029095 
    ## ... Similar to previous best
    ## Run 293 stress 0.1356147 
    ## Run 294 stress 0.1016164 
    ## ... Procrustes: rmse 7.651717e-05  max resid 0.0002146277 
    ## ... Similar to previous best
    ## Run 295 stress 0.1341868 
    ## Run 296 stress 0.1356146 
    ## Run 297 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003378589  max resid 0.0009488007 
    ## ... Similar to previous best
    ## Run 298 stress 0.1341868 
    ## Run 299 stress 0.1016164 
    ## ... Procrustes: rmse 5.165294e-05  max resid 0.000144187 
    ## ... Similar to previous best
    ## Run 300 stress 0.1016164 
    ## ... Procrustes: rmse 3.033654e-05  max resid 8.345688e-05 
    ## ... Similar to previous best
    ## Run 301 stress 0.1341868 
    ## Run 302 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003763214  max resid 0.001052129 
    ## ... Similar to previous best
    ## Run 303 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005625463  max resid 0.001571953 
    ## ... Similar to previous best
    ## Run 304 stress 0.1519033 
    ## Run 305 stress 0.1016164 
    ## ... Procrustes: rmse 3.15332e-05  max resid 8.698545e-05 
    ## ... Similar to previous best
    ## Run 306 stress 0.1356151 
    ## Run 307 stress 0.101617 
    ## ... Procrustes: rmse 0.0006973881  max resid 0.001949581 
    ## ... Similar to previous best
    ## Run 308 stress 0.1341868 
    ## Run 309 stress 0.1341868 
    ## Run 310 stress 0.1341868 
    ## Run 311 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004616505  max resid 0.001291531 
    ## ... Similar to previous best
    ## Run 312 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001277237  max resid 0.0003598824 
    ## ... Similar to previous best
    ## Run 313 stress 0.1016164 
    ## ... Procrustes: rmse 1.953489e-05  max resid 5.454447e-05 
    ## ... Similar to previous best
    ## Run 314 stress 0.1341868 
    ## Run 315 stress 0.1341868 
    ## Run 316 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004688736  max resid 0.001307431 
    ## ... Similar to previous best
    ## Run 317 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002466072  max resid 0.000691838 
    ## ... Similar to previous best
    ## Run 318 stress 0.1016171 
    ## ... Procrustes: rmse 0.0006943542  max resid 0.001936205 
    ## ... Similar to previous best
    ## Run 319 stress 0.1356144 
    ## Run 320 stress 0.1356148 
    ## Run 321 stress 0.1356148 
    ## Run 322 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005244198  max resid 0.001464964 
    ## ... Similar to previous best
    ## Run 323 stress 0.1341868 
    ## Run 324 stress 0.1341868 
    ## Run 325 stress 0.101617 
    ## ... Procrustes: rmse 0.000647217  max resid 0.001807297 
    ## ... Similar to previous best
    ## Run 326 stress 0.1016164 
    ## ... Procrustes: rmse 1.452043e-05  max resid 4.106603e-05 
    ## ... Similar to previous best
    ## Run 327 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006231546  max resid 0.001741015 
    ## ... Similar to previous best
    ## Run 328 stress 0.1530423 
    ## Run 329 stress 0.1016164 
    ## ... Procrustes: rmse 8.745912e-05  max resid 0.0002470471 
    ## ... Similar to previous best
    ## Run 330 stress 0.1016164 
    ## ... Procrustes: rmse 2.351468e-05  max resid 6.620545e-05 
    ## ... Similar to previous best
    ## Run 331 stress 0.1016164 
    ## ... Procrustes: rmse 1.472931e-05  max resid 4.059549e-05 
    ## ... Similar to previous best
    ## Run 332 stress 0.1016164 
    ## ... Procrustes: rmse 3.628207e-05  max resid 0.0001010208 
    ## ... Similar to previous best
    ## Run 333 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005309761  max resid 0.001483441 
    ## ... Similar to previous best
    ## Run 334 stress 0.1016164 
    ## ... Procrustes: rmse 4.018034e-05  max resid 0.0001133289 
    ## ... Similar to previous best
    ## Run 335 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003334711  max resid 0.0009338175 
    ## ... Similar to previous best
    ## Run 336 stress 0.1356148 
    ## Run 337 stress 0.1356149 
    ## Run 338 stress 0.1356146 
    ## Run 339 stress 0.1016164 
    ## ... Procrustes: rmse 2.233436e-05  max resid 5.134908e-05 
    ## ... Similar to previous best
    ## Run 340 stress 0.1356147 
    ## Run 341 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005692516  max resid 0.001591485 
    ## ... Similar to previous best
    ## Run 342 stress 0.1016164 
    ## ... Procrustes: rmse 8.475907e-05  max resid 0.0002391589 
    ## ... Similar to previous best
    ## Run 343 stress 0.1356145 
    ## Run 344 stress 0.1341868 
    ## Run 345 stress 0.1519033 
    ## Run 346 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002101015  max resid 0.0005895008 
    ## ... Similar to previous best
    ## Run 347 stress 0.1530423 
    ## Run 348 stress 0.1016164 
    ## ... Procrustes: rmse 2.541865e-05  max resid 5.784908e-05 
    ## ... Similar to previous best
    ## Run 349 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001953203  max resid 0.0005466441 
    ## ... Similar to previous best
    ## Run 350 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003241042  max resid 0.0009083638 
    ## ... Similar to previous best
    ## Run 351 stress 0.1016164 
    ## ... Procrustes: rmse 3.49894e-05  max resid 9.871692e-05 
    ## ... Similar to previous best
    ## Run 352 stress 0.1016164 
    ## ... Procrustes: rmse 9.792671e-06  max resid 2.225263e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.1519033 
    ## Run 354 stress 0.1341868 
    ## Run 355 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001609798  max resid 0.0004514193 
    ## ... Similar to previous best
    ## Run 356 stress 0.1016164 
    ## ... Procrustes: rmse 1.125493e-05  max resid 3.180689e-05 
    ## ... Similar to previous best
    ## Run 357 stress 0.1530424 
    ## Run 358 stress 0.1016164 
    ## ... Procrustes: rmse 1.032889e-05  max resid 2.913615e-05 
    ## ... Similar to previous best
    ## Run 359 stress 0.1341868 
    ## Run 360 stress 0.1016164 
    ## ... Procrustes: rmse 6.532888e-05  max resid 0.0001826776 
    ## ... Similar to previous best
    ## Run 361 stress 0.1016164 
    ## ... Procrustes: rmse 9.620762e-05  max resid 0.0002708559 
    ## ... Similar to previous best
    ## Run 362 stress 0.1016164 
    ## ... Procrustes: rmse 3.761413e-05  max resid 0.0001061508 
    ## ... Similar to previous best
    ## Run 363 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002075143  max resid 0.0005821942 
    ## ... Similar to previous best
    ## Run 364 stress 0.1356148 
    ## Run 365 stress 0.1341868 
    ## Run 366 stress 0.1016164 
    ## ... Procrustes: rmse 3.322616e-05  max resid 9.318561e-05 
    ## ... Similar to previous best
    ## Run 367 stress 0.1016169 
    ## ... Procrustes: rmse 0.000577608  max resid 0.001615456 
    ## ... Similar to previous best
    ## Run 368 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002830381  max resid 0.0007939744 
    ## ... Similar to previous best
    ## Run 369 stress 0.1016164 
    ## ... Procrustes: rmse 3.905991e-05  max resid 0.0001099314 
    ## ... Similar to previous best
    ## Run 370 stress 0.1016164 
    ## ... Procrustes: rmse 8.657306e-05  max resid 0.0002437455 
    ## ... Similar to previous best
    ## Run 371 stress 0.1016168 
    ## ... Procrustes: rmse 0.000518082  max resid 0.001445872 
    ## ... Similar to previous best
    ## Run 372 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003973977  max resid 0.00111195 
    ## ... Similar to previous best
    ## Run 373 stress 0.1356146 
    ## Run 374 stress 0.1016164 
    ## ... Procrustes: rmse 5.051573e-05  max resid 0.0001421175 
    ## ... Similar to previous best
    ## Run 375 stress 0.1016164 
    ## ... Procrustes: rmse 3.429386e-05  max resid 9.661051e-05 
    ## ... Similar to previous best
    ## Run 376 stress 0.1530425 
    ## Run 377 stress 0.1016164 
    ## ... Procrustes: rmse 5.166074e-05  max resid 0.0001456294 
    ## ... Similar to previous best
    ## Run 378 stress 0.1016164 
    ## ... Procrustes: rmse 2.804312e-05  max resid 5.199776e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.1341868 
    ## Run 380 stress 0.1356149 
    ## Run 381 stress 0.1356148 
    ## Run 382 stress 0.1016164 
    ## ... Procrustes: rmse 6.819592e-06  max resid 1.926879e-05 
    ## ... Similar to previous best
    ## Run 383 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002871119  max resid 0.0008045263 
    ## ... Similar to previous best
    ## Run 384 stress 0.1519033 
    ## Run 385 stress 0.1341868 
    ## Run 386 stress 0.1016164 
    ## ... Procrustes: rmse 1.738558e-05  max resid 3.069558e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002043118  max resid 0.0005698956 
    ## ... Similar to previous best
    ## Run 388 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004768919  max resid 0.001334648 
    ## ... Similar to previous best
    ## Run 389 stress 0.1016164 
    ## ... Procrustes: rmse 1.419473e-05  max resid 3.996271e-05 
    ## ... Similar to previous best
    ## Run 390 stress 0.1356147 
    ## Run 391 stress 0.1530428 
    ## Run 392 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 1.385075e-06  max resid 4.133103e-06 
    ## ... Similar to previous best
    ## Run 393 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002656593  max resid 0.0007441134 
    ## ... Similar to previous best
    ## Run 394 stress 0.1341868 
    ## Run 395 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003378347  max resid 0.0009453805 
    ## ... Similar to previous best
    ## Run 396 stress 0.1016164 
    ## ... Procrustes: rmse 7.074741e-06  max resid 1.884195e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.1341868 
    ## Run 398 stress 0.1016164 
    ## ... Procrustes: rmse 7.998763e-05  max resid 0.0002218655 
    ## ... Similar to previous best
    ## Run 399 stress 0.1341868 
    ## Run 400 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007188101  max resid 0.002008453 
    ## ... Similar to previous best
    ## Run 401 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002497003  max resid 0.0006955684 
    ## ... Similar to previous best
    ## Run 402 stress 0.1016164 
    ## ... Procrustes: rmse 9.497813e-06  max resid 2.717061e-05 
    ## ... Similar to previous best
    ## Run 403 stress 0.1341868 
    ## Run 404 stress 0.1016164 
    ## ... Procrustes: rmse 3.940632e-05  max resid 0.0001106645 
    ## ... Similar to previous best
    ## Run 405 stress 0.1341868 
    ## Run 406 stress 0.1016164 
    ## ... Procrustes: rmse 5.636828e-05  max resid 0.0001585817 
    ## ... Similar to previous best
    ## Run 407 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002785712  max resid 0.0007774428 
    ## ... Similar to previous best
    ## Run 408 stress 0.1537554 
    ## Run 409 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001119748  max resid 0.0003152278 
    ## ... Similar to previous best
    ## Run 410 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004905365  max resid 0.001364944 
    ## ... Similar to previous best
    ## Run 411 stress 0.1016164 
    ## ... Procrustes: rmse 7.320592e-05  max resid 0.0002064961 
    ## ... Similar to previous best
    ## Run 412 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002826071  max resid 0.000790668 
    ## ... Similar to previous best
    ## Run 413 stress 0.1016164 
    ## ... Procrustes: rmse 5.888488e-06  max resid 1.633451e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.1530429 
    ## Run 415 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006344315  max resid 0.001770948 
    ## ... Similar to previous best
    ## Run 416 stress 0.1341868 
    ## Run 417 stress 0.1341868 
    ## Run 418 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006262579  max resid 0.001749644 
    ## ... Similar to previous best
    ## Run 419 stress 0.1341868 
    ## Run 420 stress 0.1356149 
    ## Run 421 stress 0.1537555 
    ## Run 422 stress 0.1356149 
    ## Run 423 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007310792  max resid 0.002041736 
    ## ... Similar to previous best
    ## Run 424 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004228216  max resid 0.001178981 
    ## ... Similar to previous best
    ## Run 425 stress 0.1341868 
    ## Run 426 stress 0.1016164 
    ## ... Procrustes: rmse 6.523521e-06  max resid 1.744745e-05 
    ## ... Similar to previous best
    ## Run 427 stress 0.2074952 
    ## Run 428 stress 0.1016164 
    ## ... Procrustes: rmse 3.633228e-05  max resid 0.0001023513 
    ## ... Similar to previous best
    ## Run 429 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005395002  max resid 0.001507789 
    ## ... Similar to previous best
    ## Run 430 stress 0.1341868 
    ## Run 431 stress 0.1519033 
    ## Run 432 stress 0.1016164 
    ## ... Procrustes: rmse 6.734101e-06  max resid 1.891824e-05 
    ## ... Similar to previous best
    ## Run 433 stress 0.1016164 
    ## ... Procrustes: rmse 4.312854e-05  max resid 0.0001217633 
    ## ... Similar to previous best
    ## Run 434 stress 0.1016164 
    ## ... Procrustes: rmse 5.191017e-05  max resid 0.0001466965 
    ## ... Similar to previous best
    ## Run 435 stress 0.1016164 
    ## ... Procrustes: rmse 3.320725e-05  max resid 9.363502e-05 
    ## ... Similar to previous best
    ## Run 436 stress 0.1356147 
    ## Run 437 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004679588  max resid 0.001308343 
    ## ... Similar to previous best
    ## Run 438 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001291841  max resid 0.0003624757 
    ## ... Similar to previous best
    ## Run 439 stress 0.1016164 
    ## ... Procrustes: rmse 1.180394e-05  max resid 2.489052e-05 
    ## ... Similar to previous best
    ## Run 440 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003570609  max resid 0.0009989279 
    ## ... Similar to previous best
    ## Run 441 stress 0.1530423 
    ## Run 442 stress 0.1356144 
    ## Run 443 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003200986  max resid 0.0008956442 
    ## ... Similar to previous best
    ## Run 444 stress 0.1016164 
    ## ... Procrustes: rmse 4.526709e-05  max resid 0.0001276135 
    ## ... Similar to previous best
    ## Run 445 stress 0.1356146 
    ## Run 446 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002080369  max resid 0.0005783408 
    ## ... Similar to previous best
    ## Run 447 stress 0.1016164 
    ## ... Procrustes: rmse 3.624592e-05  max resid 9.995395e-05 
    ## ... Similar to previous best
    ## Run 448 stress 0.1016164 
    ## ... Procrustes: rmse 5.794216e-05  max resid 0.0001626132 
    ## ... Similar to previous best
    ## Run 449 stress 0.1016164 
    ## ... Procrustes: rmse 3.092943e-05  max resid 8.315781e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.1016164 
    ## ... Procrustes: rmse 5.136195e-05  max resid 0.0001422981 
    ## ... Similar to previous best
    ## Run 451 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003183634  max resid 0.0008914926 
    ## ... Similar to previous best
    ## Run 452 stress 0.1016164 
    ## ... Procrustes: rmse 4.014128e-06  max resid 5.940116e-06 
    ## ... Similar to previous best
    ## Run 453 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005838945  max resid 0.001631411 
    ## ... Similar to previous best
    ## Run 454 stress 0.1016164 
    ## ... Procrustes: rmse 5.842302e-05  max resid 0.0001646564 
    ## ... Similar to previous best
    ## Run 455 stress 0.1016164 
    ## ... Procrustes: rmse 1.43454e-05  max resid 4.058553e-05 
    ## ... Similar to previous best
    ## Run 456 stress 0.1016164 
    ## ... Procrustes: rmse 5.306174e-06  max resid 1.533189e-05 
    ## ... Similar to previous best
    ## Run 457 stress 0.1016164 
    ## ... Procrustes: rmse 7.295858e-06  max resid 2.027417e-05 
    ## ... Similar to previous best
    ## Run 458 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003794857  max resid 0.001062006 
    ## ... Similar to previous best
    ## Run 459 stress 0.1016164 
    ## ... Procrustes: rmse 4.051575e-05  max resid 0.0001138776 
    ## ... Similar to previous best
    ## Run 460 stress 0.1537555 
    ## Run 461 stress 0.1341868 
    ## Run 462 stress 0.1341868 
    ## Run 463 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005040031  max resid 0.001402531 
    ## ... Similar to previous best
    ## Run 464 stress 0.135615 
    ## Run 465 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005930574  max resid 0.001655945 
    ## ... Similar to previous best
    ## Run 466 stress 0.1530422 
    ## Run 467 stress 0.1016164 
    ## ... Procrustes: rmse 1.720789e-05  max resid 4.797643e-05 
    ## ... Similar to previous best
    ## Run 468 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002657109  max resid 0.0007452854 
    ## ... Similar to previous best
    ## Run 469 stress 0.1016164 
    ## ... Procrustes: rmse 9.144529e-06  max resid 2.55845e-05 
    ## ... Similar to previous best
    ## Run 470 stress 0.1341868 
    ## Run 471 stress 0.1341868 
    ## Run 472 stress 0.1016164 
    ## ... Procrustes: rmse 8.605908e-05  max resid 0.0002424255 
    ## ... Similar to previous best
    ## Run 473 stress 0.1016164 
    ## ... Procrustes: rmse 5.963396e-05  max resid 0.0001683855 
    ## ... Similar to previous best
    ## Run 474 stress 0.1016164 
    ## ... Procrustes: rmse 6.339503e-05  max resid 0.0001787152 
    ## ... Similar to previous best
    ## Run 475 stress 0.1016164 
    ## ... Procrustes: rmse 2.349445e-05  max resid 6.629874e-05 
    ## ... Similar to previous best
    ## Run 476 stress 0.1537557 
    ## Run 477 stress 0.1016164 
    ## ... Procrustes: rmse 2.334122e-05  max resid 6.619193e-05 
    ## ... Similar to previous best
    ## Run 478 stress 0.1341868 
    ## Run 479 stress 0.1341868 
    ## Run 480 stress 0.1016169 
    ## ... Procrustes: rmse 0.000563543  max resid 0.001579566 
    ## ... Similar to previous best
    ## Run 481 stress 0.101617 
    ## ... Procrustes: rmse 0.000697967  max resid 0.001950516 
    ## ... Similar to previous best
    ## Run 482 stress 0.1016164 
    ## ... Procrustes: rmse 1.952549e-05  max resid 4.86211e-05 
    ## ... Similar to previous best
    ## Run 483 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004353365  max resid 0.001218048 
    ## ... Similar to previous best
    ## Run 484 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001023646  max resid 0.0002847747 
    ## ... Similar to previous best
    ## Run 485 stress 0.1530424 
    ## Run 486 stress 0.1016164 
    ## ... Procrustes: rmse 1.267537e-05  max resid 3.596445e-05 
    ## ... Similar to previous best
    ## Run 487 stress 0.1537555 
    ## Run 488 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003411879  max resid 0.0009555489 
    ## ... Similar to previous best
    ## Run 489 stress 0.135615 
    ## Run 490 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001079165  max resid 0.0003028165 
    ## ... Similar to previous best
    ## Run 491 stress 0.1341868 
    ## Run 492 stress 0.1016164 
    ## ... Procrustes: rmse 3.231638e-06  max resid 8.758951e-06 
    ## ... Similar to previous best
    ## Run 493 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001668465  max resid 0.0004659738 
    ## ... Similar to previous best
    ## Run 494 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004988608  max resid 0.001394904 
    ## ... Similar to previous best
    ## Run 495 stress 0.1341868 
    ## Run 496 stress 0.1341868 
    ## Run 497 stress 0.1016164 
    ## ... Procrustes: rmse 6.333451e-06  max resid 1.818357e-05 
    ## ... Similar to previous best
    ## Run 498 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005700009  max resid 0.001592388 
    ## ... Similar to previous best
    ## Run 499 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001392859  max resid 0.0003907901 
    ## ... Similar to previous best
    ## Run 500 stress 0.1530424 
    ## *** Best solution repeated 71 times

``` r
# Stratified lakes and ocean sites
PD_beta_geo_SO_NMDS <- metaMDS(PD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.04851417 
    ## Run 1 stress 0.0493862 
    ## Run 2 stress 0.06291425 
    ## Run 3 stress 0.06088971 
    ## Run 4 stress 0.05572688 
    ## Run 5 stress 0.06205755 
    ## Run 6 stress 0.05728066 
    ## Run 7 stress 0.05643263 
    ## Run 8 stress 0.04938615 
    ## Run 9 stress 0.05695873 
    ## Run 10 stress 0.05728089 
    ## Run 11 stress 0.05643263 
    ## Run 12 stress 0.05474813 
    ## Run 13 stress 0.06291427 
    ## Run 14 stress 0.05728059 
    ## Run 15 stress 0.04851419 
    ## ... Procrustes: rmse 2.9932e-05  max resid 7.806588e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.05959895 
    ## Run 17 stress 0.04938615 
    ## Run 18 stress 0.05439281 
    ## Run 19 stress 0.04938617 
    ## Run 20 stress 0.05439282 
    ## Run 21 stress 0.04612763 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04983951  max resid 0.139599 
    ## Run 22 stress 0.04851417 
    ## Run 23 stress 0.06291433 
    ## Run 24 stress 0.05474822 
    ## Run 25 stress 0.04938615 
    ## Run 26 stress 0.04851417 
    ## Run 27 stress 0.05612239 
    ## Run 28 stress 0.04938617 
    ## Run 29 stress 0.06205762 
    ## Run 30 stress 0.05451493 
    ## Run 31 stress 0.05260055 
    ## Run 32 stress 0.05137968 
    ## Run 33 stress 0.05607344 
    ## Run 34 stress 0.06205755 
    ## Run 35 stress 0.04612763 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002494396  max resid 0.0005082335 
    ## ... Similar to previous best
    ## Run 36 stress 0.04612759 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002071223  max resid 0.0004213966 
    ## ... Similar to previous best
    ## Run 37 stress 0.04851418 
    ## Run 38 stress 0.05276412 
    ## Run 39 stress 0.04938615 
    ## Run 40 stress 0.05612253 
    ## Run 41 stress 0.05276411 
    ## Run 42 stress 0.04938617 
    ## Run 43 stress 0.05593457 
    ## Run 44 stress 0.05643261 
    ## Run 45 stress 0.04938615 
    ## Run 46 stress 0.0543895 
    ## Run 47 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 8.457368e-05  max resid 0.000171778 
    ## ... Similar to previous best
    ## Run 48 stress 0.05216142 
    ## Run 49 stress 0.06291425 
    ## Run 50 stress 0.04938615 
    ## Run 51 stress 0.04935946 
    ## Run 52 stress 0.05438952 
    ## Run 53 stress 0.04860422 
    ## Run 54 stress 0.04851419 
    ## Run 55 stress 0.0527641 
    ## Run 56 stress 0.06400728 
    ## Run 57 stress 0.04938616 
    ## Run 58 stress 0.05728087 
    ## Run 59 stress 0.05216136 
    ## Run 60 stress 0.04938615 
    ## Run 61 stress 0.04938616 
    ## Run 62 stress 0.05433173 
    ## Run 63 stress 0.04851418 
    ## Run 64 stress 0.04938616 
    ## Run 65 stress 0.04938616 
    ## Run 66 stress 0.05605155 
    ## Run 67 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 1.188976e-05  max resid 1.920725e-05 
    ## ... Similar to previous best
    ## Run 68 stress 0.04851417 
    ## Run 69 stress 0.05607344 
    ## Run 70 stress 0.04851417 
    ## Run 71 stress 0.05260047 
    ## Run 72 stress 0.05216163 
    ## Run 73 stress 0.05488111 
    ## Run 74 stress 0.04938616 
    ## Run 75 stress 0.04612758 
    ## ... Procrustes: rmse 6.031579e-05  max resid 0.000105529 
    ## ... Similar to previous best
    ## Run 76 stress 0.05438936 
    ## Run 77 stress 0.04612758 
    ## ... Procrustes: rmse 4.785455e-05  max resid 9.572807e-05 
    ## ... Similar to previous best
    ## Run 78 stress 0.04851417 
    ## Run 79 stress 0.05643261 
    ## Run 80 stress 0.04851417 
    ## Run 81 stress 0.05728067 
    ## Run 82 stress 0.05488106 
    ## Run 83 stress 0.04851418 
    ## Run 84 stress 0.05607348 
    ## Run 85 stress 0.04851417 
    ## Run 86 stress 0.05260058 
    ## Run 87 stress 0.05855032 
    ## Run 88 stress 0.05433174 
    ## Run 89 stress 0.0620576 
    ## Run 90 stress 0.05959897 
    ## Run 91 stress 0.05607346 
    ## Run 92 stress 0.06291432 
    ## Run 93 stress 0.05643262 
    ## Run 94 stress 0.04851418 
    ## Run 95 stress 0.04860423 
    ## Run 96 stress 0.05643261 
    ## Run 97 stress 0.06933824 
    ## Run 98 stress 0.04612759 
    ## ... Procrustes: rmse 8.42879e-05  max resid 0.0001705964 
    ## ... Similar to previous best
    ## Run 99 stress 0.04851418 
    ## Run 100 stress 0.05607344 
    ## Run 101 stress 0.04612757 
    ## ... Procrustes: rmse 2.462183e-05  max resid 4.568227e-05 
    ## ... Similar to previous best
    ## Run 102 stress 0.04935943 
    ## Run 103 stress 0.06291424 
    ## Run 104 stress 0.05474814 
    ## Run 105 stress 0.05605147 
    ## Run 106 stress 0.0608897 
    ## Run 107 stress 0.05643261 
    ## Run 108 stress 0.05438939 
    ## Run 109 stress 0.06205753 
    ## Run 110 stress 0.05216148 
    ## Run 111 stress 0.05605153 
    ## Run 112 stress 0.06291427 
    ## Run 113 stress 0.05593449 
    ## Run 114 stress 0.05438948 
    ## Run 115 stress 0.05605149 
    ## Run 116 stress 0.05474814 
    ## Run 117 stress 0.04938617 
    ## Run 118 stress 0.05605149 
    ## Run 119 stress 0.05216149 
    ## Run 120 stress 0.05474819 
    ## Run 121 stress 0.0527641 
    ## Run 122 stress 0.05593466 
    ## Run 123 stress 0.0572806 
    ## Run 124 stress 0.06291425 
    ## Run 125 stress 0.05612245 
    ## Run 126 stress 0.05605146 
    ## Run 127 stress 0.05260048 
    ## Run 128 stress 0.05959902 
    ## Run 129 stress 0.05216131 
    ## Run 130 stress 0.04935947 
    ## Run 131 stress 0.05728068 
    ## Run 132 stress 0.05433176 
    ## Run 133 stress 0.05607345 
    ## Run 134 stress 0.04612759 
    ## ... Procrustes: rmse 7.258783e-05  max resid 0.0001489718 
    ## ... Similar to previous best
    ## Run 135 stress 0.05612249 
    ## Run 136 stress 0.05605154 
    ## Run 137 stress 0.05605151 
    ## Run 138 stress 0.04860422 
    ## Run 139 stress 0.06400727 
    ## Run 140 stress 0.05695854 
    ## Run 141 stress 0.05488106 
    ## Run 142 stress 0.04935945 
    ## Run 143 stress 0.06205752 
    ## Run 144 stress 0.05855034 
    ## Run 145 stress 0.05728069 
    ## Run 146 stress 0.04851418 
    ## Run 147 stress 0.05842786 
    ## Run 148 stress 0.06088971 
    ## Run 149 stress 0.06088973 
    ## Run 150 stress 0.04938616 
    ## Run 151 stress 0.05607344 
    ## Run 152 stress 0.04938616 
    ## Run 153 stress 0.07046331 
    ## Run 154 stress 0.04935942 
    ## Run 155 stress 0.04935948 
    ## Run 156 stress 0.04860423 
    ## Run 157 stress 0.0572807 
    ## Run 158 stress 0.04860422 
    ## Run 159 stress 0.06291425 
    ## Run 160 stress 0.05607345 
    ## Run 161 stress 0.05572688 
    ## Run 162 stress 0.05643261 
    ## Run 163 stress 0.0527641 
    ## Run 164 stress 0.05216134 
    ## Run 165 stress 0.05605149 
    ## Run 166 stress 0.04851417 
    ## Run 167 stress 0.04851417 
    ## Run 168 stress 0.05593452 
    ## Run 169 stress 0.06205754 
    ## Run 170 stress 0.05474829 
    ## Run 171 stress 0.04851417 
    ## Run 172 stress 0.05612239 
    ## Run 173 stress 0.04612758 
    ## ... Procrustes: rmse 4.665939e-05  max resid 9.621438e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.04851417 
    ## Run 175 stress 0.04935945 
    ## Run 176 stress 0.04935944 
    ## Run 177 stress 0.05695871 
    ## Run 178 stress 0.06400724 
    ## Run 179 stress 0.04938616 
    ## Run 180 stress 0.05612242 
    ## Run 181 stress 0.05216136 
    ## Run 182 stress 0.05842809 
    ## Run 183 stress 0.05605145 
    ## Run 184 stress 0.04935946 
    ## Run 185 stress 0.0543928 
    ## Run 186 stress 0.06205753 
    ## Run 187 stress 0.04935952 
    ## Run 188 stress 0.05260053 
    ## Run 189 stress 0.0629143 
    ## Run 190 stress 0.05643261 
    ## Run 191 stress 0.05276412 
    ## Run 192 stress 0.05959903 
    ## Run 193 stress 0.05593456 
    ## Run 194 stress 0.0585503 
    ## Run 195 stress 0.05438943 
    ## Run 196 stress 0.05607345 
    ## Run 197 stress 0.05260048 
    ## Run 198 stress 0.05451494 
    ## Run 199 stress 0.04612758 
    ## ... Procrustes: rmse 5.353168e-05  max resid 0.0001043338 
    ## ... Similar to previous best
    ## Run 200 stress 0.04938616 
    ## Run 201 stress 0.04860422 
    ## Run 202 stress 0.05474814 
    ## Run 203 stress 0.04612757 
    ## ... Procrustes: rmse 1.158344e-05  max resid 1.969207e-05 
    ## ... Similar to previous best
    ## Run 204 stress 0.04938615 
    ## Run 205 stress 0.05855031 
    ## Run 206 stress 0.05438936 
    ## Run 207 stress 0.04851417 
    ## Run 208 stress 0.05260056 
    ## Run 209 stress 0.0543928 
    ## Run 210 stress 0.0485142 
    ## Run 211 stress 0.05605155 
    ## Run 212 stress 0.04612762 
    ## ... Procrustes: rmse 0.000102778  max resid 0.0002094138 
    ## ... Similar to previous best
    ## Run 213 stress 0.05855031 
    ## Run 214 stress 0.05607344 
    ## Run 215 stress 0.06205754 
    ## Run 216 stress 0.05216163 
    ## Run 217 stress 0.05474817 
    ## Run 218 stress 0.04851417 
    ## Run 219 stress 0.05593457 
    ## Run 220 stress 0.04860422 
    ## Run 221 stress 0.05695858 
    ## Run 222 stress 0.05605154 
    ## Run 223 stress 0.04851417 
    ## Run 224 stress 0.06291425 
    ## Run 225 stress 0.05216171 
    ## Run 226 stress 0.0561224 
    ## Run 227 stress 0.05137964 
    ## Run 228 stress 0.05607345 
    ## Run 229 stress 0.05216125 
    ## Run 230 stress 0.0608897 
    ## Run 231 stress 0.05728053 
    ## Run 232 stress 0.3357174 
    ## Run 233 stress 0.04851417 
    ## Run 234 stress 0.05643267 
    ## Run 235 stress 0.0557269 
    ## Run 236 stress 0.04860423 
    ## Run 237 stress 0.04938618 
    ## Run 238 stress 0.05593449 
    ## Run 239 stress 0.04938615 
    ## Run 240 stress 0.05433178 
    ## Run 241 stress 0.05959896 
    ## Run 242 stress 0.05260047 
    ## Run 243 stress 0.04860423 
    ## Run 244 stress 0.3578482 
    ## Run 245 stress 0.05643261 
    ## Run 246 stress 0.0527641 
    ## Run 247 stress 0.05728058 
    ## Run 248 stress 0.05605148 
    ## Run 249 stress 0.05593451 
    ## Run 250 stress 0.04860422 
    ## Run 251 stress 0.06400731 
    ## Run 252 stress 0.05439279 
    ## Run 253 stress 0.04851417 
    ## Run 254 stress 0.06291425 
    ## Run 255 stress 0.04612758 
    ## ... Procrustes: rmse 5.951748e-05  max resid 0.0001196091 
    ## ... Similar to previous best
    ## Run 256 stress 0.05438945 
    ## Run 257 stress 0.04851418 
    ## Run 258 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001198862  max resid 0.0002461768 
    ## ... Similar to previous best
    ## Run 259 stress 0.0572807 
    ## Run 260 stress 0.06088975 
    ## Run 261 stress 0.05607344 
    ## Run 262 stress 0.05438937 
    ## Run 263 stress 0.05216158 
    ## Run 264 stress 0.04612765 
    ## ... Procrustes: rmse 0.0001473046  max resid 0.0002972252 
    ## ... Similar to previous best
    ## Run 265 stress 0.05572693 
    ## Run 266 stress 0.05728052 
    ## Run 267 stress 0.0547482 
    ## Run 268 stress 0.05728073 
    ## Run 269 stress 0.05643261 
    ## Run 270 stress 0.05572689 
    ## Run 271 stress 0.05605151 
    ## Run 272 stress 0.05439283 
    ## Run 273 stress 0.05488105 
    ## Run 274 stress 0.05695886 
    ## Run 275 stress 0.04860423 
    ## Run 276 stress 0.05137965 
    ## Run 277 stress 0.04938616 
    ## Run 278 stress 0.05474817 
    ## Run 279 stress 0.05607344 
    ## Run 280 stress 0.0559345 
    ## Run 281 stress 0.05572692 
    ## Run 282 stress 0.05260048 
    ## Run 283 stress 0.04612759 
    ## ... Procrustes: rmse 6.713847e-05  max resid 0.0001363362 
    ## ... Similar to previous best
    ## Run 284 stress 0.04851417 
    ## Run 285 stress 0.04612773 
    ## ... Procrustes: rmse 0.0001629602  max resid 0.0003377798 
    ## ... Similar to previous best
    ## Run 286 stress 0.05216127 
    ## Run 287 stress 0.05959898 
    ## Run 288 stress 0.04851417 
    ## Run 289 stress 0.05605152 
    ## Run 290 stress 0.05474814 
    ## Run 291 stress 0.05728057 
    ## Run 292 stress 0.05216146 
    ## Run 293 stress 0.05695874 
    ## Run 294 stress 0.04935953 
    ## Run 295 stress 0.06205753 
    ## Run 296 stress 0.05855031 
    ## Run 297 stress 0.05842813 
    ## Run 298 stress 0.04938619 
    ## Run 299 stress 0.05276413 
    ## Run 300 stress 0.05643265 
    ## Run 301 stress 0.06088976 
    ## Run 302 stress 0.04851418 
    ## Run 303 stress 0.05855036 
    ## Run 304 stress 0.05137972 
    ## Run 305 stress 0.05474813 
    ## Run 306 stress 0.04851417 
    ## Run 307 stress 0.05842815 
    ## Run 308 stress 0.05439282 
    ## Run 309 stress 0.05605149 
    ## Run 310 stress 0.05728073 
    ## Run 311 stress 0.05474829 
    ## Run 312 stress 0.04938615 
    ## Run 313 stress 0.05643261 
    ## Run 314 stress 0.06291429 
    ## Run 315 stress 0.05643261 
    ## Run 316 stress 0.04612757 
    ## ... Procrustes: rmse 1.605392e-05  max resid 3.214307e-05 
    ## ... Similar to previous best
    ## Run 317 stress 0.05572688 
    ## Run 318 stress 0.05643261 
    ## Run 319 stress 0.0585503 
    ## Run 320 stress 0.05593453 
    ## Run 321 stress 0.04938616 
    ## Run 322 stress 0.05137964 
    ## Run 323 stress 0.0493862 
    ## Run 324 stress 0.04860422 
    ## Run 325 stress 0.0557269 
    ## Run 326 stress 0.04851417 
    ## Run 327 stress 0.05438936 
    ## Run 328 stress 0.0557269 
    ## Run 329 stress 0.05728064 
    ## Run 330 stress 0.05607347 
    ## Run 331 stress 0.05695879 
    ## Run 332 stress 0.05137967 
    ## Run 333 stress 0.05593456 
    ## Run 334 stress 0.06400729 
    ## Run 335 stress 0.05216126 
    ## Run 336 stress 0.0461276 
    ## ... Procrustes: rmse 0.0001010159  max resid 0.000204074 
    ## ... Similar to previous best
    ## Run 337 stress 0.06205752 
    ## Run 338 stress 0.05643265 
    ## Run 339 stress 0.05612243 
    ## Run 340 stress 0.05607344 
    ## Run 341 stress 0.05451493 
    ## Run 342 stress 0.05607344 
    ## Run 343 stress 0.05474817 
    ## Run 344 stress 0.05855031 
    ## Run 345 stress 0.04851417 
    ## Run 346 stress 0.05605149 
    ## Run 347 stress 0.05842789 
    ## Run 348 stress 0.05643264 
    ## Run 349 stress 0.05842781 
    ## Run 350 stress 0.06291426 
    ## Run 351 stress 0.05260051 
    ## Run 352 stress 0.05842791 
    ## Run 353 stress 0.05216144 
    ## Run 354 stress 0.05451493 
    ## Run 355 stress 0.05695855 
    ## Run 356 stress 0.06205752 
    ## Run 357 stress 0.04935953 
    ## Run 358 stress 0.0557269 
    ## Run 359 stress 0.04938617 
    ## Run 360 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001182054  max resid 0.0002426407 
    ## ... Similar to previous best
    ## Run 361 stress 0.05607344 
    ## Run 362 stress 0.05605158 
    ## Run 363 stress 0.04935944 
    ## Run 364 stress 0.04938616 
    ## Run 365 stress 0.05216132 
    ## Run 366 stress 0.04851417 
    ## Run 367 stress 0.06400739 
    ## Run 368 stress 0.06205755 
    ## Run 369 stress 0.05607344 
    ## Run 370 stress 0.0640074 
    ## Run 371 stress 0.05728056 
    ## Run 372 stress 0.05605156 
    ## Run 373 stress 0.05842794 
    ## Run 374 stress 0.05728073 
    ## Run 375 stress 0.04935951 
    ## Run 376 stress 0.07046385 
    ## Run 377 stress 0.05607345 
    ## Run 378 stress 0.05438938 
    ## Run 379 stress 0.04938616 
    ## Run 380 stress 0.05643261 
    ## Run 381 stress 0.05216131 
    ## Run 382 stress 0.04938615 
    ## Run 383 stress 0.04860422 
    ## Run 384 stress 0.04612758 
    ## ... Procrustes: rmse 4.626943e-05  max resid 9.01669e-05 
    ## ... Similar to previous best
    ## Run 385 stress 0.05842783 
    ## Run 386 stress 0.04938616 
    ## Run 387 stress 0.0585503 
    ## Run 388 stress 0.0557269 
    ## Run 389 stress 0.05605149 
    ## Run 390 stress 0.05216174 
    ## Run 391 stress 0.05643263 
    ## Run 392 stress 0.06291425 
    ## Run 393 stress 0.04860422 
    ## Run 394 stress 0.05474814 
    ## Run 395 stress 0.05216144 
    ## Run 396 stress 0.04851418 
    ## Run 397 stress 0.05612247 
    ## Run 398 stress 0.05260049 
    ## Run 399 stress 0.05474814 
    ## Run 400 stress 0.3442052 
    ## Run 401 stress 0.05572688 
    ## Run 402 stress 0.06291424 
    ## Run 403 stress 0.05643263 
    ## Run 404 stress 0.0526005 
    ## Run 405 stress 0.04860422 
    ## Run 406 stress 0.04612761 
    ## ... Procrustes: rmse 9.032974e-05  max resid 0.0001814971 
    ## ... Similar to previous best
    ## Run 407 stress 0.04938616 
    ## Run 408 stress 0.0585503 
    ## Run 409 stress 0.04938619 
    ## Run 410 stress 0.04860422 
    ## Run 411 stress 0.0561224 
    ## Run 412 stress 0.04938618 
    ## Run 413 stress 0.05438947 
    ## Run 414 stress 0.05695864 
    ## Run 415 stress 0.05855032 
    ## Run 416 stress 0.04938617 
    ## Run 417 stress 0.0569588 
    ## Run 418 stress 0.04612757 
    ## ... Procrustes: rmse 3.668199e-05  max resid 7.489088e-05 
    ## ... Similar to previous best
    ## Run 419 stress 0.05451496 
    ## Run 420 stress 0.04860422 
    ## Run 421 stress 0.05605159 
    ## Run 422 stress 0.04938616 
    ## Run 423 stress 0.05439283 
    ## Run 424 stress 0.05643261 
    ## Run 425 stress 0.04935947 
    ## Run 426 stress 0.05605148 
    ## Run 427 stress 0.05572688 
    ## Run 428 stress 0.06205752 
    ## Run 429 stress 0.0573151 
    ## Run 430 stress 0.05605147 
    ## Run 431 stress 0.05438937 
    ## Run 432 stress 0.04851418 
    ## Run 433 stress 0.05607345 
    ## Run 434 stress 0.05855029 
    ## Run 435 stress 0.06291425 
    ## Run 436 stress 0.04851417 
    ## Run 437 stress 0.05137966 
    ## Run 438 stress 0.05216138 
    ## Run 439 stress 0.0526006 
    ## Run 440 stress 0.0461276 
    ## ... Procrustes: rmse 7.845887e-05  max resid 0.0001580725 
    ## ... Similar to previous best
    ## Run 441 stress 0.05855029 
    ## Run 442 stress 0.04938616 
    ## Run 443 stress 0.04851419 
    ## Run 444 stress 0.04851417 
    ## Run 445 stress 0.05728054 
    ## Run 446 stress 0.06088976 
    ## Run 447 stress 0.04938616 
    ## Run 448 stress 0.0584282 
    ## Run 449 stress 0.05959897 
    ## Run 450 stress 0.05438941 
    ## Run 451 stress 0.04851417 
    ## Run 452 stress 0.05474814 
    ## Run 453 stress 0.04612766 
    ## ... Procrustes: rmse 8.650671e-05  max resid 0.0001785735 
    ## ... Similar to previous best
    ## Run 454 stress 0.04860422 
    ## Run 455 stress 0.05855029 
    ## Run 456 stress 0.06205753 
    ## Run 457 stress 0.05643261 
    ## Run 458 stress 0.05728065 
    ## Run 459 stress 0.0572809 
    ## Run 460 stress 0.05572688 
    ## Run 461 stress 0.06291425 
    ## Run 462 stress 0.04851417 
    ## Run 463 stress 0.05474814 
    ## Run 464 stress 0.04938616 
    ## Run 465 stress 0.05855031 
    ## Run 466 stress 0.3303068 
    ## Run 467 stress 0.05451494 
    ## Run 468 stress 0.05842783 
    ## Run 469 stress 0.04935943 
    ## Run 470 stress 0.05216138 
    ## Run 471 stress 0.05605154 
    ## Run 472 stress 0.05643266 
    ## Run 473 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001250036  max resid 0.0002545894 
    ## ... Similar to previous best
    ## Run 474 stress 0.04851418 
    ## Run 475 stress 0.06291425 
    ## Run 476 stress 0.05474814 
    ## Run 477 stress 0.05731514 
    ## Run 478 stress 0.06400762 
    ## Run 479 stress 0.05216165 
    ## Run 480 stress 0.05731511 
    ## Run 481 stress 0.06185043 
    ## Run 482 stress 0.05451493 
    ## Run 483 stress 0.05607344 
    ## Run 484 stress 0.05612255 
    ## Run 485 stress 0.05260051 
    ## Run 486 stress 0.05959898 
    ## Run 487 stress 0.05842815 
    ## Run 488 stress 0.05438944 
    ## Run 489 stress 0.06088971 
    ## Run 490 stress 0.05605153 
    ## Run 491 stress 0.05488105 
    ## Run 492 stress 0.04935958 
    ## Run 493 stress 0.06291432 
    ## Run 494 stress 0.05643261 
    ## Run 495 stress 0.0527641 
    ## Run 496 stress 0.0561224 
    ## Run 497 stress 0.05593449 
    ## Run 498 stress 0.05474818 
    ## Run 499 stress 0.05607344 
    ## Run 500 stress 0.05695883 
    ## *** Best solution repeated 24 times

``` r
# Mixed lakes
PD_beta_geo_M_NMDS <- metaMDS(PD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.019724e-05 
    ## Run 1 stress 9.853594e-05 
    ## ... Procrustes: rmse 0.0001907711  max resid 0.0003675335 
    ## ... Similar to previous best
    ## Run 2 stress 9.940126e-05 
    ## ... Procrustes: rmse 1.95906e-05  max resid 3.7717e-05 
    ## ... Similar to previous best
    ## Run 3 stress 9.883443e-05 
    ## ... Procrustes: rmse 0.0001536755  max resid 0.0002867923 
    ## ... Similar to previous best
    ## Run 4 stress 9.184415e-05 
    ## ... Procrustes: rmse 0.0001523665  max resid 0.0002757072 
    ## ... Similar to previous best
    ## Run 5 stress 9.480423e-05 
    ## ... Procrustes: rmse 2.733063e-05  max resid 4.316826e-05 
    ## ... Similar to previous best
    ## Run 6 stress 8.684462e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001667676  max resid 0.0003376403 
    ## ... Similar to previous best
    ## Run 7 stress 8.868223e-05 
    ## ... Procrustes: rmse 0.00018257  max resid 0.0003181817 
    ## ... Similar to previous best
    ## Run 8 stress 7.864028e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001902548  max resid 0.0003396059 
    ## ... Similar to previous best
    ## Run 9 stress 9.663981e-05 
    ## ... Procrustes: rmse 0.000256998  max resid 0.000598132 
    ## ... Similar to previous best
    ## Run 10 stress 9.90381e-05 
    ## ... Procrustes: rmse 0.000262777  max resid 0.0006107803 
    ## ... Similar to previous best
    ## Run 11 stress 9.556798e-05 
    ## ... Procrustes: rmse 0.0001912157  max resid 0.0003265249 
    ## ... Similar to previous best
    ## Run 12 stress 9.435879e-05 
    ## ... Procrustes: rmse 0.0002515429  max resid 0.0005914032 
    ## ... Similar to previous best
    ## Run 13 stress 8.860957e-05 
    ## ... Procrustes: rmse 0.0001114635  max resid 0.000222773 
    ## ... Similar to previous best
    ## Run 14 stress 9.594064e-05 
    ## ... Procrustes: rmse 0.0001356244  max resid 0.0003135672 
    ## ... Similar to previous best
    ## Run 15 stress 9.068566e-05 
    ## ... Procrustes: rmse 0.0001315598  max resid 0.000257835 
    ## ... Similar to previous best
    ## Run 16 stress 9.409568e-05 
    ## ... Procrustes: rmse 0.0001232121  max resid 0.0002936181 
    ## ... Similar to previous best
    ## Run 17 stress 9.373983e-05 
    ## ... Procrustes: rmse 0.0001176017  max resid 0.0002209019 
    ## ... Similar to previous best
    ## Run 18 stress 9.535443e-05 
    ## ... Procrustes: rmse 0.000171757  max resid 0.0002549254 
    ## ... Similar to previous best
    ## Run 19 stress 8.775992e-05 
    ## ... Procrustes: rmse 0.0002320654  max resid 0.0005492001 
    ## ... Similar to previous best
    ## Run 20 stress 9.051446e-05 
    ## ... Procrustes: rmse 0.0002435517  max resid 0.0005706585 
    ## ... Similar to previous best
    ## Run 21 stress 9.479215e-05 
    ## ... Procrustes: rmse 0.0002474436  max resid 0.0005754278 
    ## ... Similar to previous best
    ## Run 22 stress 9.577087e-05 
    ## ... Procrustes: rmse 0.0002551643  max resid 0.0005938352 
    ## ... Similar to previous best
    ## Run 23 stress 9.223668e-05 
    ## ... Procrustes: rmse 0.0001250466  max resid 0.0002449158 
    ## ... Similar to previous best
    ## Run 24 stress 9.774446e-05 
    ## ... Procrustes: rmse 0.0002596344  max resid 0.0006024004 
    ## ... Similar to previous best
    ## Run 25 stress 9.497797e-05 
    ## ... Procrustes: rmse 0.0002091994  max resid 0.0003740708 
    ## ... Similar to previous best
    ## Run 26 stress 9.511264e-05 
    ## ... Procrustes: rmse 0.0002459995  max resid 0.000572236 
    ## ... Similar to previous best
    ## Run 27 stress 9.457172e-05 
    ## ... Procrustes: rmse 0.0002091504  max resid 0.0003727852 
    ## ... Similar to previous best
    ## Run 28 stress 9.12117e-05 
    ## ... Procrustes: rmse 0.0002391328  max resid 0.0005580643 
    ## ... Similar to previous best
    ## Run 29 stress 9.591735e-05 
    ## ... Procrustes: rmse 0.0002560766  max resid 0.0005972124 
    ## ... Similar to previous best
    ## Run 30 stress 9.379823e-05 
    ## ... Procrustes: rmse 0.0001958734  max resid 0.0003632537 
    ## ... Similar to previous best
    ## Run 31 stress 9.903319e-05 
    ## ... Procrustes: rmse 0.0002076024  max resid 0.0003852198 
    ## ... Similar to previous best
    ## Run 32 stress 0.3083098 
    ## Run 33 stress 8.666757e-05 
    ## ... Procrustes: rmse 0.0002272672  max resid 0.0005425345 
    ## ... Similar to previous best
    ## Run 34 stress 9.609838e-05 
    ## ... Procrustes: rmse 0.000209157  max resid 0.0003770595 
    ## ... Similar to previous best
    ## Run 35 stress 8.882475e-05 
    ## ... Procrustes: rmse 0.0002386694  max resid 0.0005630181 
    ## ... Similar to previous best
    ## Run 36 stress 9.016436e-05 
    ## ... Procrustes: rmse 0.0001984753  max resid 0.0003515549 
    ## ... Similar to previous best
    ## Run 37 stress 9.112292e-05 
    ## ... Procrustes: rmse 0.0001831294  max resid 0.000309288 
    ## ... Similar to previous best
    ## Run 38 stress 9.703939e-05 
    ## ... Procrustes: rmse 0.0001373363  max resid 0.0002469383 
    ## ... Similar to previous best
    ## Run 39 stress 8.332785e-05 
    ## ... Procrustes: rmse 1.655728e-05  max resid 3.057254e-05 
    ## ... Similar to previous best
    ## Run 40 stress 8.686022e-05 
    ## ... Procrustes: rmse 0.0001911215  max resid 0.0003473347 
    ## ... Similar to previous best
    ## Run 41 stress 8.967215e-05 
    ## ... Procrustes: rmse 0.0001185247  max resid 0.0002316239 
    ## ... Similar to previous best
    ## Run 42 stress 9.88404e-05 
    ## ... Procrustes: rmse 0.000211154  max resid 0.0003845295 
    ## ... Similar to previous best
    ## Run 43 stress 9.539395e-05 
    ## ... Procrustes: rmse 0.0002518868  max resid 0.0005869301 
    ## ... Similar to previous best
    ## Run 44 stress 0.3082558 
    ## Run 45 stress 9.128004e-05 
    ## ... Procrustes: rmse 0.000132651  max resid 0.0003029681 
    ## ... Similar to previous best
    ## Run 46 stress 9.618556e-05 
    ## ... Procrustes: rmse 0.0002532102  max resid 0.0005878955 
    ## ... Similar to previous best
    ## Run 47 stress 0.2319132 
    ## Run 48 stress 8.529929e-05 
    ## ... Procrustes: rmse 0.0001253895  max resid 0.0002559782 
    ## ... Similar to previous best
    ## Run 49 stress 9.843435e-05 
    ## ... Procrustes: rmse 0.0001961227  max resid 0.0003374585 
    ## ... Similar to previous best
    ## Run 50 stress 9.942779e-05 
    ## ... Procrustes: rmse 0.0001960916  max resid 0.0003381394 
    ## ... Similar to previous best
    ## Run 51 stress 9.470782e-05 
    ## ... Procrustes: rmse 0.0002091268  max resid 0.0003698549 
    ## ... Similar to previous best
    ## Run 52 stress 9.268876e-05 
    ## ... Procrustes: rmse 0.000186882  max resid 0.0003215569 
    ## ... Similar to previous best
    ## Run 53 stress 9.503171e-05 
    ## ... Procrustes: rmse 0.0002474518  max resid 0.0005812193 
    ## ... Similar to previous best
    ## Run 54 stress 9.185878e-05 
    ## ... Procrustes: rmse 0.0002111702  max resid 0.0004507561 
    ## ... Similar to previous best
    ## Run 55 stress 9.65928e-05 
    ## ... Procrustes: rmse 0.0001325526  max resid 0.0002541266 
    ## ... Similar to previous best
    ## Run 56 stress 9.63649e-05 
    ## ... Procrustes: rmse 0.0001943319  max resid 0.0003252613 
    ## ... Similar to previous best
    ## Run 57 stress 9.786367e-05 
    ## ... Procrustes: rmse 0.0002584618  max resid 0.0005981314 
    ## ... Similar to previous best
    ## Run 58 stress 9.803095e-05 
    ## ... Procrustes: rmse 0.0002593263  max resid 0.000605297 
    ## ... Similar to previous best
    ## Run 59 stress 9.804786e-05 
    ## ... Procrustes: rmse 0.0002020731  max resid 0.0003830557 
    ## ... Similar to previous best
    ## Run 60 stress 9.345276e-05 
    ## ... Procrustes: rmse 0.0001882214  max resid 0.0003202692 
    ## ... Similar to previous best
    ## Run 61 stress 9.70225e-05 
    ## ... Procrustes: rmse 0.0002539713  max resid 0.0005951986 
    ## ... Similar to previous best
    ## Run 62 stress 8.886651e-05 
    ## ... Procrustes: rmse 0.0001299086  max resid 0.0002548588 
    ## ... Similar to previous best
    ## Run 63 stress 9.937778e-05 
    ## ... Procrustes: rmse 0.0002636262  max resid 0.0006144134 
    ## ... Similar to previous best
    ## Run 64 stress 9.01555e-05 
    ## ... Procrustes: rmse 0.0002433326  max resid 0.0005712666 
    ## ... Similar to previous best
    ## Run 65 stress 9.462705e-05 
    ## ... Procrustes: rmse 0.0001858322  max resid 0.0003230697 
    ## ... Similar to previous best
    ## Run 66 stress 8.907032e-05 
    ## ... Procrustes: rmse 3.668437e-05  max resid 6.707472e-05 
    ## ... Similar to previous best
    ## Run 67 stress 9.794518e-05 
    ## ... Procrustes: rmse 2.462253e-05  max resid 3.837444e-05 
    ## ... Similar to previous best
    ## Run 68 stress 9.469049e-05 
    ## ... Procrustes: rmse 0.0001363305  max resid 0.0002269717 
    ## ... Similar to previous best
    ## Run 69 stress 9.402464e-05 
    ## ... Procrustes: rmse 0.0001344496  max resid 0.000270395 
    ## ... Similar to previous best
    ## Run 70 stress 9.860007e-05 
    ## ... Procrustes: rmse 0.0001931054  max resid 0.000339576 
    ## ... Similar to previous best
    ## Run 71 stress 9.254551e-05 
    ## ... Procrustes: rmse 7.79452e-05  max resid 0.000106854 
    ## ... Similar to previous best
    ## Run 72 stress 6.486534e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001861223  max resid 0.000452094 
    ## ... Similar to previous best
    ## Run 73 stress 9.167879e-05 
    ## ... Procrustes: rmse 0.0001585108  max resid 0.000291732 
    ## ... Similar to previous best
    ## Run 74 stress 9.629848e-05 
    ## ... Procrustes: rmse 0.0001899392  max resid 0.0004628337 
    ## ... Similar to previous best
    ## Run 75 stress 9.766015e-05 
    ## ... Procrustes: rmse 8.870494e-05  max resid 0.0001526582 
    ## ... Similar to previous best
    ## Run 76 stress 9.734572e-05 
    ## ... Procrustes: rmse 8.5751e-05  max resid 0.000143858 
    ## ... Similar to previous best
    ## Run 77 stress 9.446687e-05 
    ## ... Procrustes: rmse 7.951954e-05  max resid 0.0001316329 
    ## ... Similar to previous best
    ## Run 78 stress 8.902735e-05 
    ## ... Procrustes: rmse 0.0001401894  max resid 0.0003046164 
    ## ... Similar to previous best
    ## Run 79 stress 9.99786e-05 
    ## ... Procrustes: rmse 0.0001690888  max resid 0.0003020425 
    ## ... Similar to previous best
    ## Run 80 stress 8.767971e-05 
    ## ... Procrustes: rmse 6.576821e-05  max resid 0.0001199858 
    ## ... Similar to previous best
    ## Run 81 stress 9.360566e-05 
    ## ... Procrustes: rmse 0.0001643362  max resid 0.0002979971 
    ## ... Similar to previous best
    ## Run 82 stress 9.97348e-05 
    ## ... Procrustes: rmse 8.803257e-05  max resid 0.0001527274 
    ## ... Similar to previous best
    ## Run 83 stress 9.756837e-05 
    ## ... Procrustes: rmse 0.0001529824  max resid 0.0003226605 
    ## ... Similar to previous best
    ## Run 84 stress 9.236963e-05 
    ## ... Procrustes: rmse 0.0001210698  max resid 0.0002589685 
    ## ... Similar to previous best
    ## Run 85 stress 9.415554e-05 
    ## ... Procrustes: rmse 0.0001816238  max resid 0.0004427886 
    ## ... Similar to previous best
    ## Run 86 stress 9.962139e-05 
    ## ... Procrustes: rmse 9.38699e-05  max resid 0.0001558641 
    ## ... Similar to previous best
    ## Run 87 stress 9.115348e-05 
    ## ... Procrustes: rmse 0.0001816387  max resid 0.000444931 
    ## ... Similar to previous best
    ## Run 88 stress 9.342618e-05 
    ## ... Procrustes: rmse 7.790016e-05  max resid 0.0001327267 
    ## ... Similar to previous best
    ## Run 89 stress 9.230526e-05 
    ## ... Procrustes: rmse 0.0001672316  max resid 0.0003163662 
    ## ... Similar to previous best
    ## Run 90 stress 9.20757e-05 
    ## ... Procrustes: rmse 0.0001647349  max resid 0.000299969 
    ## ... Similar to previous best
    ## Run 91 stress 0.2848526 
    ## Run 92 stress 9.683472e-05 
    ## ... Procrustes: rmse 0.0001854021  max resid 0.0002302472 
    ## ... Similar to previous best
    ## Run 93 stress 9.837304e-05 
    ## ... Procrustes: rmse 0.0001578657  max resid 0.0002941432 
    ## ... Similar to previous best
    ## Run 94 stress 9.663929e-05 
    ## ... Procrustes: rmse 8.566703e-05  max resid 0.000147241 
    ## ... Similar to previous best
    ## Run 95 stress 9.756889e-05 
    ## ... Procrustes: rmse 0.0001272461  max resid 0.0002753463 
    ## ... Similar to previous best
    ## Run 96 stress 9.574469e-05 
    ## ... Procrustes: rmse 8.284329e-05  max resid 0.0001376881 
    ## ... Similar to previous best
    ## Run 97 stress 9.280658e-05 
    ## ... Procrustes: rmse 0.0001214123  max resid 0.0002598867 
    ## ... Similar to previous best
    ## Run 98 stress 8.730249e-05 
    ## ... Procrustes: rmse 0.0001220694  max resid 0.0002409896 
    ## ... Similar to previous best
    ## Run 99 stress 9.939737e-05 
    ## ... Procrustes: rmse 8.881941e-05  max resid 0.0001544525 
    ## ... Similar to previous best
    ## Run 100 stress 8.580205e-05 
    ## ... Procrustes: rmse 5.70032e-05  max resid 9.708132e-05 
    ## ... Similar to previous best
    ## Run 101 stress 8.984388e-05 
    ## ... Procrustes: rmse 0.0001604259  max resid 0.0002777172 
    ## ... Similar to previous best
    ## Run 102 stress 9.549829e-05 
    ## ... Procrustes: rmse 0.00017098  max resid 0.0003056992 
    ## ... Similar to previous best
    ## Run 103 stress 9.769966e-05 
    ## ... Procrustes: rmse 8.572088e-05  max resid 0.0001446926 
    ## ... Similar to previous best
    ## Run 104 stress 9.898861e-05 
    ## ... Procrustes: rmse 8.449194e-05  max resid 0.0001472777 
    ## ... Similar to previous best
    ## Run 105 stress 9.867602e-05 
    ## ... Procrustes: rmse 9.00478e-05  max resid 0.0001497379 
    ## ... Similar to previous best
    ## Run 106 stress 9.584129e-05 
    ## ... Procrustes: rmse 0.0001748401  max resid 0.0003218241 
    ## ... Similar to previous best
    ## Run 107 stress 9.911841e-05 
    ## ... Procrustes: rmse 8.561348e-05  max resid 0.000158499 
    ## ... Similar to previous best
    ## Run 108 stress 9.70447e-05 
    ## ... Procrustes: rmse 0.0001253607  max resid 0.0002699383 
    ## ... Similar to previous best
    ## Run 109 stress 9.755837e-05 
    ## ... Procrustes: rmse 0.0001208721  max resid 0.0002510311 
    ## ... Similar to previous best
    ## Run 110 stress 9.872115e-05 
    ## ... Procrustes: rmse 9.092844e-05  max resid 0.0001563159 
    ## ... Similar to previous best
    ## Run 111 stress 9.811175e-05 
    ## ... Procrustes: rmse 0.0001732401  max resid 0.000306395 
    ## ... Similar to previous best
    ## Run 112 stress 8.156454e-05 
    ## ... Procrustes: rmse 0.0001551922  max resid 0.0002927038 
    ## ... Similar to previous best
    ## Run 113 stress 9.839166e-05 
    ## ... Procrustes: rmse 0.0001268221  max resid 0.0002735577 
    ## ... Similar to previous best
    ## Run 114 stress 9.255844e-05 
    ## ... Procrustes: rmse 0.0001232117  max resid 0.0002497762 
    ## ... Similar to previous best
    ## Run 115 stress 9.159363e-05 
    ## ... Procrustes: rmse 0.0001221712  max resid 0.0002619382 
    ## ... Similar to previous best
    ## Run 116 stress 9.490127e-05 
    ## ... Procrustes: rmse 0.0001159453  max resid 0.0002433686 
    ## ... Similar to previous best
    ## Run 117 stress 8.527421e-05 
    ## ... Procrustes: rmse 0.0001749624  max resid 0.0003171413 
    ## ... Similar to previous best
    ## Run 118 stress 7.927014e-05 
    ## ... Procrustes: rmse 0.0001621368  max resid 0.0002751228 
    ## ... Similar to previous best
    ## Run 119 stress 9.974131e-05 
    ## ... Procrustes: rmse 0.000129578  max resid 0.0002809028 
    ## ... Similar to previous best
    ## Run 120 stress 9.134372e-05 
    ## ... Procrustes: rmse 6.762545e-05  max resid 0.0001237944 
    ## ... Similar to previous best
    ## Run 121 stress 9.966679e-05 
    ## ... Procrustes: rmse 0.0001792156  max resid 0.0003144962 
    ## ... Similar to previous best
    ## Run 122 stress 9.917239e-05 
    ## ... Procrustes: rmse 0.0001303761  max resid 0.0002820386 
    ## ... Similar to previous best
    ## Run 123 stress 9.808147e-05 
    ## ... Procrustes: rmse 9.051504e-05  max resid 0.0001561682 
    ## ... Similar to previous best
    ## Run 124 stress 9.97598e-05 
    ## ... Procrustes: rmse 0.0001772696  max resid 0.0003088951 
    ## ... Similar to previous best
    ## Run 125 stress 9.14571e-05 
    ## ... Procrustes: rmse 0.0001212444  max resid 0.0002586024 
    ## ... Similar to previous best
    ## Run 126 stress 9.624701e-05 
    ## ... Procrustes: rmse 0.0001263053  max resid 0.0002725996 
    ## ... Similar to previous best
    ## Run 127 stress 9.521049e-05 
    ## ... Procrustes: rmse 7.296597e-05  max resid 0.0001277235 
    ## ... Similar to previous best
    ## Run 128 stress 9.848595e-05 
    ## ... Procrustes: rmse 8.638705e-05  max resid 0.0001450478 
    ## ... Similar to previous best
    ## Run 129 stress 9.49205e-05 
    ## ... Procrustes: rmse 0.0001989246  max resid 0.0002684196 
    ## ... Similar to previous best
    ## Run 130 stress 9.103974e-05 
    ## ... Procrustes: rmse 0.0001632041  max resid 0.0002793962 
    ## ... Similar to previous best
    ## Run 131 stress 8.924318e-05 
    ## ... Procrustes: rmse 0.0001763585  max resid 0.0003211986 
    ## ... Similar to previous best
    ## Run 132 stress 9.32342e-05 
    ## ... Procrustes: rmse 0.0001217717  max resid 0.0002620969 
    ## ... Similar to previous best
    ## Run 133 stress 9.396968e-05 
    ## ... Procrustes: rmse 7.880193e-05  max resid 0.0001296983 
    ## ... Similar to previous best
    ## Run 134 stress 9.254897e-05 
    ## ... Procrustes: rmse 0.0001225047  max resid 0.0002631247 
    ## ... Similar to previous best
    ## Run 135 stress 9.593207e-05 
    ## ... Procrustes: rmse 0.0001256314  max resid 0.0002713442 
    ## ... Similar to previous best
    ## Run 136 stress 9.389223e-05 
    ## ... Procrustes: rmse 0.0001198032  max resid 0.0002672231 
    ## ... Similar to previous best
    ## Run 137 stress 9.994185e-05 
    ## ... Procrustes: rmse 0.0001798102  max resid 0.0003175194 
    ## ... Similar to previous best
    ## Run 138 stress 9.665953e-05 
    ## ... Procrustes: rmse 0.0001881604  max resid 0.0003005516 
    ## ... Similar to previous best
    ## Run 139 stress 9.561336e-05 
    ## ... Procrustes: rmse 0.0001664311  max resid 0.0003122229 
    ## ... Similar to previous best
    ## Run 140 stress 9.806231e-05 
    ## ... Procrustes: rmse 8.404337e-05  max resid 0.0001457401 
    ## ... Similar to previous best
    ## Run 141 stress 9.895358e-05 
    ## ... Procrustes: rmse 0.0001756734  max resid 0.0003151411 
    ## ... Similar to previous best
    ## Run 142 stress 9.431171e-05 
    ## ... Procrustes: rmse 0.0001658683  max resid 0.0003077655 
    ## ... Similar to previous best
    ## Run 143 stress 9.587473e-05 
    ## ... Procrustes: rmse 8.307797e-05  max resid 0.0001380361 
    ## ... Similar to previous best
    ## Run 144 stress 9.956177e-05 
    ## ... Procrustes: rmse 0.0001533295  max resid 0.000245901 
    ## ... Similar to previous best
    ## Run 145 stress 9.861041e-05 
    ## ... Procrustes: rmse 8.695478e-05  max resid 0.000150754 
    ## ... Similar to previous best
    ## Run 146 stress 9.283008e-05 
    ## ... Procrustes: rmse 6.807985e-05  max resid 0.0001117343 
    ## ... Similar to previous best
    ## Run 147 stress 9.867615e-05 
    ## ... Procrustes: rmse 0.0001771473  max resid 0.0003136554 
    ## ... Similar to previous best
    ## Run 148 stress 8.694198e-05 
    ## ... Procrustes: rmse 0.0001667691  max resid 0.0002880411 
    ## ... Similar to previous best
    ## Run 149 stress 8.02555e-05 
    ## ... Procrustes: rmse 0.0001704887  max resid 0.0002865245 
    ## ... Similar to previous best
    ## Run 150 stress 9.39159e-05 
    ## ... Procrustes: rmse 0.0001507759  max resid 0.0002929072 
    ## ... Similar to previous best
    ## Run 151 stress 9.955994e-05 
    ## ... Procrustes: rmse 9.393023e-05  max resid 0.0001584832 
    ## ... Similar to previous best
    ## Run 152 stress 9.304282e-05 
    ## ... Procrustes: rmse 0.0001636358  max resid 0.0002966716 
    ## ... Similar to previous best
    ## Run 153 stress 8.948519e-05 
    ## ... Procrustes: rmse 8.941583e-05  max resid 0.0002023116 
    ## ... Similar to previous best
    ## Run 154 stress 8.422457e-05 
    ## ... Procrustes: rmse 0.0001423566  max resid 0.0002739312 
    ## ... Similar to previous best
    ## Run 155 stress 9.44775e-05 
    ## ... Procrustes: rmse 7.834839e-05  max resid 0.0001289389 
    ## ... Similar to previous best
    ## Run 156 stress 9.490856e-05 
    ## ... Procrustes: rmse 7.999879e-05  max resid 0.0001331567 
    ## ... Similar to previous best
    ## Run 157 stress 9.435232e-05 
    ## ... Procrustes: rmse 0.0001673618  max resid 0.0003011564 
    ## ... Similar to previous best
    ## Run 158 stress 8.815687e-05 
    ## ... Procrustes: rmse 6.357454e-05  max resid 0.0001082566 
    ## ... Similar to previous best
    ## Run 159 stress 9.10664e-05 
    ## ... Procrustes: rmse 6.360928e-05  max resid 0.0001084867 
    ## ... Similar to previous best
    ## Run 160 stress 9.363745e-05 
    ## ... Procrustes: rmse 7.748416e-05  max resid 0.0001344646 
    ## ... Similar to previous best
    ## Run 161 stress 9.958553e-05 
    ## ... Procrustes: rmse 0.0001782095  max resid 0.0003166499 
    ## ... Similar to previous best
    ## Run 162 stress 9.159747e-05 
    ## ... Procrustes: rmse 7.310015e-05  max resid 0.0001268307 
    ## ... Similar to previous best
    ## Run 163 stress 9.806174e-05 
    ## ... Procrustes: rmse 0.0001724443  max resid 0.0003065163 
    ## ... Similar to previous best
    ## Run 164 stress 0.3073311 
    ## Run 165 stress 9.875288e-05 
    ## ... Procrustes: rmse 0.0001285781  max resid 0.000279411 
    ## ... Similar to previous best
    ## Run 166 stress 9.064797e-05 
    ## ... Procrustes: rmse 0.0001223585  max resid 0.0002607587 
    ## ... Similar to previous best
    ## Run 167 stress 9.432553e-05 
    ## ... Procrustes: rmse 7.090036e-05  max resid 0.0001231082 
    ## ... Similar to previous best
    ## Run 168 stress 0.2848513 
    ## Run 169 stress 8.557512e-05 
    ## ... Procrustes: rmse 5.646925e-05  max resid 9.578243e-05 
    ## ... Similar to previous best
    ## Run 170 stress 9.743843e-05 
    ## ... Procrustes: rmse 0.0001282988  max resid 0.0002784669 
    ## ... Similar to previous best
    ## Run 171 stress 8.78302e-05 
    ## ... Procrustes: rmse 0.0001547106  max resid 0.0002578785 
    ## ... Similar to previous best
    ## Run 172 stress 8.15011e-05 
    ## ... Procrustes: rmse 3.106047e-05  max resid 5.019067e-05 
    ## ... Similar to previous best
    ## Run 173 stress 9.410202e-05 
    ## ... Procrustes: rmse 0.0001161727  max resid 0.0002601549 
    ## ... Similar to previous best
    ## Run 174 stress 9.551108e-05 
    ## ... Procrustes: rmse 8.315166e-05  max resid 0.0001383338 
    ## ... Similar to previous best
    ## Run 175 stress 9.742609e-05 
    ## ... Procrustes: rmse 8.554994e-05  max resid 0.0001421447 
    ## ... Similar to previous best
    ## Run 176 stress 8.935614e-05 
    ## ... Procrustes: rmse 6.488211e-05  max resid 0.0001126663 
    ## ... Similar to previous best
    ## Run 177 stress 8.912377e-05 
    ## ... Procrustes: rmse 6.411967e-05  max resid 0.0001065022 
    ## ... Similar to previous best
    ## Run 178 stress 9.998513e-05 
    ## ... Procrustes: rmse 0.0001795771  max resid 0.0003170049 
    ## ... Similar to previous best
    ## Run 179 stress 9.721411e-05 
    ## ... Procrustes: rmse 8.567526e-05  max resid 0.0001421652 
    ## ... Similar to previous best
    ## Run 180 stress 9.898102e-05 
    ## ... Procrustes: rmse 0.0001966886  max resid 0.0004777078 
    ## ... Similar to previous best
    ## Run 181 stress 9.12584e-05 
    ## ... Procrustes: rmse 0.0001628449  max resid 0.0002966702 
    ## ... Similar to previous best
    ## Run 182 stress 8.670686e-05 
    ## ... Procrustes: rmse 0.0001571681  max resid 0.0003222466 
    ## ... Similar to previous best
    ## Run 183 stress 9.941987e-05 
    ## ... Procrustes: rmse 8.494703e-05  max resid 0.0001488107 
    ## ... Similar to previous best
    ## Run 184 stress 0.3083098 
    ## Run 185 stress 9.618097e-05 
    ## ... Procrustes: rmse 0.000163558  max resid 0.0003028136 
    ## ... Similar to previous best
    ## Run 186 stress 9.484071e-05 
    ## ... Procrustes: rmse 0.000124586  max resid 0.0002686843 
    ## ... Similar to previous best
    ## Run 187 stress 9.600621e-05 
    ## ... Procrustes: rmse 0.0001705508  max resid 0.0003238164 
    ## ... Similar to previous best
    ## Run 188 stress 9.633862e-05 
    ## ... Procrustes: rmse 8.525802e-05  max resid 0.0001471509 
    ## ... Similar to previous best
    ## Run 189 stress 9.493323e-05 
    ## ... Procrustes: rmse 0.0001695972  max resid 0.0003076658 
    ## ... Similar to previous best
    ## Run 190 stress 8.895186e-05 
    ## ... Procrustes: rmse 0.000145214  max resid 0.0003230677 
    ## ... Similar to previous best
    ## Run 191 stress 9.869345e-05 
    ## ... Procrustes: rmse 8.799227e-05  max resid 0.0001504029 
    ## ... Similar to previous best
    ## Run 192 stress 9.626277e-05 
    ## ... Procrustes: rmse 0.0001187  max resid 0.0002447925 
    ## ... Similar to previous best
    ## Run 193 stress 8.603451e-05 
    ## ... Procrustes: rmse 0.0001751584  max resid 0.0003145809 
    ## ... Similar to previous best
    ## Run 194 stress 9.634279e-05 
    ## ... Procrustes: rmse 0.0001695845  max resid 0.0003046702 
    ## ... Similar to previous best
    ## Run 195 stress 8.149793e-05 
    ## ... Procrustes: rmse 0.0001719379  max resid 0.0002963828 
    ## ... Similar to previous best
    ## Run 196 stress 9.755209e-05 
    ## ... Procrustes: rmse 9.879362e-05  max resid 0.0001711232 
    ## ... Similar to previous best
    ## Run 197 stress 9.777645e-05 
    ## ... Procrustes: rmse 0.0001257874  max resid 0.0002707505 
    ## ... Similar to previous best
    ## Run 198 stress 7.838983e-05 
    ## ... Procrustes: rmse 3.485232e-05  max resid 4.673263e-05 
    ## ... Similar to previous best
    ## Run 199 stress 9.662939e-05 
    ## ... Procrustes: rmse 0.0001639022  max resid 0.0002874623 
    ## ... Similar to previous best
    ## Run 200 stress 9.845149e-05 
    ## ... Procrustes: rmse 9.000633e-05  max resid 0.0001543003 
    ## ... Similar to previous best
    ## Run 201 stress 9.142717e-05 
    ## ... Procrustes: rmse 0.0001937536  max resid 0.0002798745 
    ## ... Similar to previous best
    ## Run 202 stress 8.525764e-05 
    ## ... Procrustes: rmse 0.0001164664  max resid 0.0002479048 
    ## ... Similar to previous best
    ## Run 203 stress 9.198966e-05 
    ## ... Procrustes: rmse 0.0001229192  max resid 0.0002644898 
    ## ... Similar to previous best
    ## Run 204 stress 9.221164e-05 
    ## ... Procrustes: rmse 0.0001222922  max resid 0.0002626022 
    ## ... Similar to previous best
    ## Run 205 stress 9.568877e-05 
    ## ... Procrustes: rmse 8.122207e-05  max resid 0.0001398058 
    ## ... Similar to previous best
    ## Run 206 stress 0.3083098 
    ## Run 207 stress 8.789791e-05 
    ## ... Procrustes: rmse 3.976954e-05  max resid 6.537408e-05 
    ## ... Similar to previous best
    ## Run 208 stress 9.818203e-05 
    ## ... Procrustes: rmse 9.008717e-05  max resid 0.0001508409 
    ## ... Similar to previous best
    ## Run 209 stress 9.619034e-05 
    ## ... Procrustes: rmse 0.0001715308  max resid 0.0003053581 
    ## ... Similar to previous best
    ## Run 210 stress 9.131961e-05 
    ## ... Procrustes: rmse 0.0001611966  max resid 0.0002985241 
    ## ... Similar to previous best
    ## Run 211 stress 9.507697e-05 
    ## ... Procrustes: rmse 0.0001386715  max resid 0.0002744409 
    ## ... Similar to previous best
    ## Run 212 stress 8.643259e-05 
    ## ... Procrustes: rmse 0.0001446216  max resid 0.000244252 
    ## ... Similar to previous best
    ## Run 213 stress 9.571049e-05 
    ## ... Procrustes: rmse 0.0001691099  max resid 0.0003024004 
    ## ... Similar to previous best
    ## Run 214 stress 9.214468e-05 
    ## ... Procrustes: rmse 0.0001721111  max resid 0.0002932413 
    ## ... Similar to previous best
    ## Run 215 stress 9.145121e-05 
    ## ... Procrustes: rmse 0.0001016288  max resid 0.0001703867 
    ## ... Similar to previous best
    ## Run 216 stress 9.102555e-05 
    ## ... Procrustes: rmse 0.0001729863  max resid 0.0002940549 
    ## ... Similar to previous best
    ## Run 217 stress 9.860767e-05 
    ## ... Procrustes: rmse 9.01875e-05  max resid 0.0001534114 
    ## ... Similar to previous best
    ## Run 218 stress 9.674963e-05 
    ## ... Procrustes: rmse 0.0001260329  max resid 0.0002736384 
    ## ... Similar to previous best
    ## Run 219 stress 9.77682e-05 
    ## ... Procrustes: rmse 0.0001714353  max resid 0.0003625319 
    ## ... Similar to previous best
    ## Run 220 stress 9.992031e-05 
    ## ... Procrustes: rmse 7.238571e-05  max resid 0.0001247561 
    ## ... Similar to previous best
    ## Run 221 stress 9.415891e-05 
    ## ... Procrustes: rmse 9.353092e-05  max resid 0.0001881723 
    ## ... Similar to previous best
    ## Run 222 stress 9.021079e-05 
    ## ... Procrustes: rmse 8.098659e-05  max resid 0.0001918325 
    ## ... Similar to previous best
    ## Run 223 stress 9.590415e-05 
    ## ... Procrustes: rmse 0.0001820402  max resid 0.0003036589 
    ## ... Similar to previous best
    ## Run 224 stress 8.551008e-05 
    ## ... Procrustes: rmse 0.0001848648  max resid 0.000244438 
    ## ... Similar to previous best
    ## Run 225 stress 9.805338e-05 
    ## ... Procrustes: rmse 9.008284e-05  max resid 0.0001571147 
    ## ... Similar to previous best
    ## Run 226 stress 8.507034e-05 
    ## ... Procrustes: rmse 5.418746e-05  max resid 9.395235e-05 
    ## ... Similar to previous best
    ## Run 227 stress 9.285148e-05 
    ## ... Procrustes: rmse 0.0001765962  max resid 0.0004137333 
    ## ... Similar to previous best
    ## Run 228 stress 9.372137e-05 
    ## ... Procrustes: rmse 7.888952e-05  max resid 0.0001298835 
    ## ... Similar to previous best
    ## Run 229 stress 9.590355e-05 
    ## ... Procrustes: rmse 0.0001764023  max resid 0.0003265524 
    ## ... Similar to previous best
    ## Run 230 stress 9.936384e-05 
    ## ... Procrustes: rmse 0.0001526925  max resid 0.0003522635 
    ## ... Similar to previous best
    ## Run 231 stress 9.270036e-05 
    ## ... Procrustes: rmse 0.0001510674  max resid 0.000323897 
    ## ... Similar to previous best
    ## Run 232 stress 9.583959e-05 
    ## ... Procrustes: rmse 8.440956e-05  max resid 0.0001452326 
    ## ... Similar to previous best
    ## Run 233 stress 9.325662e-05 
    ## ... Procrustes: rmse 0.0001237632  max resid 0.0002667202 
    ## ... Similar to previous best
    ## Run 234 stress 9.13793e-05 
    ## ... Procrustes: rmse 0.0001618544  max resid 0.0002967525 
    ## ... Similar to previous best
    ## Run 235 stress 8.161031e-05 
    ## ... Procrustes: rmse 4.624011e-05  max resid 7.84644e-05 
    ## ... Similar to previous best
    ## Run 236 stress 9.39285e-05 
    ## ... Procrustes: rmse 7.771729e-05  max resid 0.0001347169 
    ## ... Similar to previous best
    ## Run 237 stress 0.3083098 
    ## Run 238 stress 9.570478e-05 
    ## ... Procrustes: rmse 8.194816e-05  max resid 0.0001413219 
    ## ... Similar to previous best
    ## Run 239 stress 9.219655e-05 
    ## ... Procrustes: rmse 0.0001717453  max resid 0.0002978975 
    ## ... Similar to previous best
    ## Run 240 stress 8.532091e-05 
    ## ... Procrustes: rmse 0.0001677153  max resid 0.0002984982 
    ## ... Similar to previous best
    ## Run 241 stress 9.617594e-05 
    ## ... Procrustes: rmse 0.0001224453  max resid 0.0002624721 
    ## ... Similar to previous best
    ## Run 242 stress 9.627752e-05 
    ## ... Procrustes: rmse 0.0001891505  max resid 0.0004601397 
    ## ... Similar to previous best
    ## Run 243 stress 9.305162e-05 
    ## ... Procrustes: rmse 0.000177058  max resid 0.0003163428 
    ## ... Similar to previous best
    ## Run 244 stress 9.645359e-05 
    ## ... Procrustes: rmse 7.117764e-05  max resid 0.000128481 
    ## ... Similar to previous best
    ## Run 245 stress 0.284852 
    ## Run 246 stress 9.098321e-05 
    ## ... Procrustes: rmse 9.608384e-05  max resid 0.0001503923 
    ## ... Similar to previous best
    ## Run 247 stress 9.674795e-05 
    ## ... Procrustes: rmse 0.0001739255  max resid 0.0003263494 
    ## ... Similar to previous best
    ## Run 248 stress 0.2848518 
    ## Run 249 stress 9.253723e-05 
    ## ... Procrustes: rmse 0.0001176887  max resid 0.0002633016 
    ## ... Similar to previous best
    ## Run 250 stress 9.545518e-05 
    ## ... Procrustes: rmse 0.0001705226  max resid 0.0003080109 
    ## ... Similar to previous best
    ## Run 251 stress 9.512348e-05 
    ## ... Procrustes: rmse 7.145197e-05  max resid 0.0001227842 
    ## ... Similar to previous best
    ## Run 252 stress 9.365968e-05 
    ## ... Procrustes: rmse 0.0001798431  max resid 0.0003051415 
    ## ... Similar to previous best
    ## Run 253 stress 9.615065e-05 
    ## ... Procrustes: rmse 8.167193e-05  max resid 0.0001343097 
    ## ... Similar to previous best
    ## Run 254 stress 9.605982e-05 
    ## ... Procrustes: rmse 0.0001136018  max resid 0.0002388548 
    ## ... Similar to previous best
    ## Run 255 stress 9.967455e-05 
    ## ... Procrustes: rmse 9.30359e-05  max resid 0.0001568045 
    ## ... Similar to previous best
    ## Run 256 stress 7.333223e-05 
    ## ... Procrustes: rmse 3.074305e-05  max resid 5.333094e-05 
    ## ... Similar to previous best
    ## Run 257 stress 9.319523e-05 
    ## ... Procrustes: rmse 0.0001193419  max resid 0.0002546749 
    ## ... Similar to previous best
    ## Run 258 stress 9.880589e-05 
    ## ... Procrustes: rmse 0.0001713137  max resid 0.0003071626 
    ## ... Similar to previous best
    ## Run 259 stress 9.032049e-05 
    ## ... Procrustes: rmse 0.0001634499  max resid 0.0002776467 
    ## ... Similar to previous best
    ## Run 260 stress 9.176704e-05 
    ## ... Procrustes: rmse 0.0001921002  max resid 0.0002820084 
    ## ... Similar to previous best
    ## Run 261 stress 8.813693e-05 
    ## ... Procrustes: rmse 9.291348e-05  max resid 0.0001989237 
    ## ... Similar to previous best
    ## Run 262 stress 9.750567e-05 
    ## ... Procrustes: rmse 0.0001712061  max resid 0.0003113101 
    ## ... Similar to previous best
    ## Run 263 stress 9.386544e-05 
    ## ... Procrustes: rmse 7.549284e-05  max resid 0.0001293212 
    ## ... Similar to previous best
    ## Run 264 stress 9.079331e-05 
    ## ... Procrustes: rmse 0.0001210313  max resid 0.0002590863 
    ## ... Similar to previous best
    ## Run 265 stress 9.767942e-05 
    ## ... Procrustes: rmse 8.581394e-05  max resid 0.0001484379 
    ## ... Similar to previous best
    ## Run 266 stress 9.113732e-05 
    ## ... Procrustes: rmse 0.0001765869  max resid 0.0003274706 
    ## ... Similar to previous best
    ## Run 267 stress 9.351193e-05 
    ## ... Procrustes: rmse 0.0001810553  max resid 0.0003358636 
    ## ... Similar to previous best
    ## Run 268 stress 9.462041e-05 
    ## ... Procrustes: rmse 0.0001231096  max resid 0.0002643739 
    ## ... Similar to previous best
    ## Run 269 stress 9.682286e-05 
    ## ... Procrustes: rmse 0.0001784849  max resid 0.0003134179 
    ## ... Similar to previous best
    ## Run 270 stress 9.653532e-05 
    ## ... Procrustes: rmse 0.0001721191  max resid 0.0003088118 
    ## ... Similar to previous best
    ## Run 271 stress 9.660295e-05 
    ## ... Procrustes: rmse 0.0001246403  max resid 0.0002679354 
    ## ... Similar to previous best
    ## Run 272 stress 9.457113e-05 
    ## ... Procrustes: rmse 8.00017e-05  max resid 0.0001384069 
    ## ... Similar to previous best
    ## Run 273 stress 9.269057e-05 
    ## ... Procrustes: rmse 0.0001747718  max resid 0.0002969975 
    ## ... Similar to previous best
    ## Run 274 stress 9.915306e-05 
    ## ... Procrustes: rmse 8.586966e-05  max resid 0.0001439689 
    ## ... Similar to previous best
    ## Run 275 stress 9.767682e-05 
    ## ... Procrustes: rmse 7.978971e-05  max resid 0.0001287115 
    ## ... Similar to previous best
    ## Run 276 stress 9.494115e-05 
    ## ... Procrustes: rmse 0.0001771093  max resid 0.0002682088 
    ## ... Similar to previous best
    ## Run 277 stress 9.55414e-05 
    ## ... Procrustes: rmse 0.0001707718  max resid 0.0003072129 
    ## ... Similar to previous best
    ## Run 278 stress 9.481298e-05 
    ## ... Procrustes: rmse 8.184086e-05  max resid 0.0001368786 
    ## ... Similar to previous best
    ## Run 279 stress 9.655055e-05 
    ## ... Procrustes: rmse 8.543785e-05  max resid 0.0001435204 
    ## ... Similar to previous best
    ## Run 280 stress 9.295661e-05 
    ## ... Procrustes: rmse 7.215584e-05  max resid 0.0001243783 
    ## ... Similar to previous best
    ## Run 281 stress 9.203506e-05 
    ## ... Procrustes: rmse 7.147695e-05  max resid 0.0001222133 
    ## ... Similar to previous best
    ## Run 282 stress 9.143089e-05 
    ## ... Procrustes: rmse 9.782233e-05  max resid 0.0001822645 
    ## ... Similar to previous best
    ## Run 283 stress 9.664441e-05 
    ## ... Procrustes: rmse 0.0001719903  max resid 0.0003105389 
    ## ... Similar to previous best
    ## Run 284 stress 9.448815e-05 
    ## ... Procrustes: rmse 7.749804e-05  max resid 0.0001329514 
    ## ... Similar to previous best
    ## Run 285 stress 9.941325e-05 
    ## ... Procrustes: rmse 7.961813e-05  max resid 0.0001533853 
    ## ... Similar to previous best
    ## Run 286 stress 9.485017e-05 
    ## ... Procrustes: rmse 0.0001247337  max resid 0.0002690144 
    ## ... Similar to previous best
    ## Run 287 stress 8.745697e-05 
    ## ... Procrustes: rmse 0.0001529341  max resid 0.0002876501 
    ## ... Similar to previous best
    ## Run 288 stress 9.517755e-05 
    ## ... Procrustes: rmse 9.648054e-05  max resid 0.0001915565 
    ## ... Similar to previous best
    ## Run 289 stress 9.963325e-05 
    ## ... Procrustes: rmse 0.0001743262  max resid 0.0003134089 
    ## ... Similar to previous best
    ## Run 290 stress 9.429009e-05 
    ## ... Procrustes: rmse 7.914584e-05  max resid 0.0001291303 
    ## ... Similar to previous best
    ## Run 291 stress 8.272954e-05 
    ## ... Procrustes: rmse 0.0001155963  max resid 0.0002417663 
    ## ... Similar to previous best
    ## Run 292 stress 9.727918e-05 
    ## ... Procrustes: rmse 0.0001737229  max resid 0.0003075216 
    ## ... Similar to previous best
    ## Run 293 stress 8.984515e-05 
    ## ... Procrustes: rmse 0.0001872917  max resid 0.0004570606 
    ## ... Similar to previous best
    ## Run 294 stress 9.78975e-05 
    ## ... Procrustes: rmse 0.0001583136  max resid 0.0002641938 
    ## ... Similar to previous best
    ## Run 295 stress 9.880414e-05 
    ## ... Procrustes: rmse 0.0001737105  max resid 0.0003069018 
    ## ... Similar to previous best
    ## Run 296 stress 0.3080294 
    ## Run 297 stress 9.972553e-05 
    ## ... Procrustes: rmse 9.333509e-05  max resid 0.0001600961 
    ## ... Similar to previous best
    ## Run 298 stress 9.535398e-05 
    ## ... Procrustes: rmse 0.0001688231  max resid 0.0003080759 
    ## ... Similar to previous best
    ## Run 299 stress 9.314498e-05 
    ## ... Procrustes: rmse 0.0001149013  max resid 0.0002488082 
    ## ... Similar to previous best
    ## Run 300 stress 9.876241e-05 
    ## ... Procrustes: rmse 0.0001895287  max resid 0.0004603734 
    ## ... Similar to previous best
    ## Run 301 stress 6.592555e-05 
    ## ... Procrustes: rmse 0.0001456186  max resid 0.0002432594 
    ## ... Similar to previous best
    ## Run 302 stress 9.555948e-05 
    ## ... Procrustes: rmse 7.819511e-05  max resid 0.0001354368 
    ## ... Similar to previous best
    ## Run 303 stress 8.751159e-05 
    ## ... Procrustes: rmse 0.0001457299  max resid 0.0001947292 
    ## ... Similar to previous best
    ## Run 304 stress 9.157953e-05 
    ## ... Procrustes: rmse 7.187723e-05  max resid 0.0001184741 
    ## ... Similar to previous best
    ## Run 305 stress 0.3083098 
    ## Run 306 stress 9.270655e-05 
    ## ... Procrustes: rmse 0.000164862  max resid 0.0003443999 
    ## ... Similar to previous best
    ## Run 307 stress 9.258033e-05 
    ## ... Procrustes: rmse 0.0001852708  max resid 0.0004512347 
    ## ... Similar to previous best
    ## Run 308 stress 9.970846e-05 
    ## ... Procrustes: rmse 6.562116e-05  max resid 0.0001155324 
    ## ... Similar to previous best
    ## Run 309 stress 9.259372e-05 
    ## ... Procrustes: rmse 0.0001933343  max resid 0.0002582694 
    ## ... Similar to previous best
    ## Run 310 stress 9.030306e-05 
    ## ... Procrustes: rmse 0.0001184731  max resid 0.0002529513 
    ## ... Similar to previous best
    ## Run 311 stress 9.622042e-05 
    ## ... Procrustes: rmse 0.0001658532  max resid 0.000305761 
    ## ... Similar to previous best
    ## Run 312 stress 9.629192e-05 
    ## ... Procrustes: rmse 0.0001714416  max resid 0.000309308 
    ## ... Similar to previous best
    ## Run 313 stress 8.475935e-05 
    ## ... Procrustes: rmse 5.446353e-05  max resid 9.181938e-05 
    ## ... Similar to previous best
    ## Run 314 stress 9.940588e-05 
    ## ... Procrustes: rmse 0.0001612133  max resid 0.0003644942 
    ## ... Similar to previous best
    ## Run 315 stress 9.240682e-05 
    ## ... Procrustes: rmse 0.0001154463  max resid 0.0002589207 
    ## ... Similar to previous best
    ## Run 316 stress 9.562062e-05 
    ## ... Procrustes: rmse 0.0001680049  max resid 0.0003013497 
    ## ... Similar to previous best
    ## Run 317 stress 9.400573e-05 
    ## ... Procrustes: rmse 0.0001155015  max resid 0.0002575159 
    ## ... Similar to previous best
    ## Run 318 stress 9.386405e-05 
    ## ... Procrustes: rmse 7.811801e-05  max resid 0.0001310229 
    ## ... Similar to previous best
    ## Run 319 stress 0.285215 
    ## Run 320 stress 9.779613e-05 
    ## ... Procrustes: rmse 8.780294e-05  max resid 0.0001460914 
    ## ... Similar to previous best
    ## Run 321 stress 9.411632e-05 
    ## ... Procrustes: rmse 0.0001694762  max resid 0.0003158577 
    ## ... Similar to previous best
    ## Run 322 stress 9.514548e-05 
    ## ... Procrustes: rmse 0.0001704435  max resid 0.0003059746 
    ## ... Similar to previous best
    ## Run 323 stress 9.615782e-05 
    ## ... Procrustes: rmse 0.0001238631  max resid 0.0002659751 
    ## ... Similar to previous best
    ## Run 324 stress 7.459217e-05 
    ## ... Procrustes: rmse 8.969811e-05  max resid 0.0001735708 
    ## ... Similar to previous best
    ## Run 325 stress 0.3075827 
    ## Run 326 stress 9.694107e-05 
    ## ... Procrustes: rmse 0.0001945107  max resid 0.0002778309 
    ## ... Similar to previous best
    ## Run 327 stress 8.760419e-05 
    ## ... Procrustes: rmse 0.0001650323  max resid 0.0003542068 
    ## ... Similar to previous best
    ## Run 328 stress 9.380974e-05 
    ## ... Procrustes: rmse 0.0001220959  max resid 0.0002615021 
    ## ... Similar to previous best
    ## Run 329 stress 7.651054e-05 
    ## ... Procrustes: rmse 0.0001628229  max resid 0.0003620158 
    ## ... Similar to previous best
    ## Run 330 stress 9.90332e-05 
    ## ... Procrustes: rmse 9.183091e-05  max resid 0.0001534571 
    ## ... Similar to previous best
    ## Run 331 stress 8.153229e-05 
    ## ... Procrustes: rmse 4.604534e-05  max resid 7.894487e-05 
    ## ... Similar to previous best
    ## Run 332 stress 8.719589e-05 
    ## ... Procrustes: rmse 0.000159621  max resid 0.0002229427 
    ## ... Similar to previous best
    ## Run 333 stress 0.3083098 
    ## Run 334 stress 9.31929e-05 
    ## ... Procrustes: rmse 7.706342e-05  max resid 0.0001312369 
    ## ... Similar to previous best
    ## Run 335 stress 9.838516e-05 
    ## ... Procrustes: rmse 0.0001676507  max resid 0.0003480661 
    ## ... Similar to previous best
    ## Run 336 stress 9.158056e-05 
    ## ... Procrustes: rmse 0.0001209617  max resid 0.0002583602 
    ## ... Similar to previous best
    ## Run 337 stress 0.3083098 
    ## Run 338 stress 9.693503e-05 
    ## ... Procrustes: rmse 0.0001249614  max resid 0.0002687302 
    ## ... Similar to previous best
    ## Run 339 stress 9.217593e-05 
    ## ... Procrustes: rmse 0.0001594925  max resid 0.0003457254 
    ## ... Similar to previous best
    ## Run 340 stress 9.301465e-05 
    ## ... Procrustes: rmse 0.0001625376  max resid 0.0002972645 
    ## ... Similar to previous best
    ## Run 341 stress 9.991008e-05 
    ## ... Procrustes: rmse 0.0001802513  max resid 0.0003401908 
    ## ... Similar to previous best
    ## Run 342 stress 9.9404e-05 
    ## ... Procrustes: rmse 0.0001543089  max resid 0.0002717968 
    ## ... Similar to previous best
    ## Run 343 stress 7.737246e-05 
    ## ... Procrustes: rmse 0.0001381631  max resid 0.0002350107 
    ## ... Similar to previous best
    ## Run 344 stress 9.633002e-05 
    ## ... Procrustes: rmse 0.0001235634  max resid 0.0002650629 
    ## ... Similar to previous best
    ## Run 345 stress 9.122965e-05 
    ## ... Procrustes: rmse 8.694665e-05  max resid 0.0001877402 
    ## ... Similar to previous best
    ## Run 346 stress 9.37784e-05 
    ## ... Procrustes: rmse 7.340565e-05  max resid 0.0001261638 
    ## ... Similar to previous best
    ## Run 347 stress 9.383516e-05 
    ## ... Procrustes: rmse 7.454414e-05  max resid 0.0001286229 
    ## ... Similar to previous best
    ## Run 348 stress 8.911767e-05 
    ## ... Procrustes: rmse 8.841959e-05  max resid 0.0001960756 
    ## ... Similar to previous best
    ## Run 349 stress 9.5852e-05 
    ## ... Procrustes: rmse 8.138963e-05  max resid 0.0001392473 
    ## ... Similar to previous best
    ## Run 350 stress 9.087201e-05 
    ## ... Procrustes: rmse 7.205079e-05  max resid 0.0001235495 
    ## ... Similar to previous best
    ## Run 351 stress 9.614166e-05 
    ## ... Procrustes: rmse 7.689733e-05  max resid 0.0001328696 
    ## ... Similar to previous best
    ## Run 352 stress 8.517504e-05 
    ## ... Procrustes: rmse 0.0001712738  max resid 0.0003008083 
    ## ... Similar to previous best
    ## Run 353 stress 9.187839e-05 
    ## ... Procrustes: rmse 0.0001206121  max resid 0.0002577422 
    ## ... Similar to previous best
    ## Run 354 stress 9.50822e-05 
    ## ... Procrustes: rmse 8.189581e-05  max resid 0.0001408689 
    ## ... Similar to previous best
    ## Run 355 stress 9.93871e-05 
    ## ... Procrustes: rmse 0.0001763764  max resid 0.0003087871 
    ## ... Similar to previous best
    ## Run 356 stress 9.686946e-05 
    ## ... Procrustes: rmse 8.516011e-05  max resid 0.0001486353 
    ## ... Similar to previous best
    ## Run 357 stress 9.679e-05 
    ## ... Procrustes: rmse 8.323854e-05  max resid 0.0001404734 
    ## ... Similar to previous best
    ## Run 358 stress 7.355093e-05 
    ## ... Procrustes: rmse 0.0001791155  max resid 0.0004345816 
    ## ... Similar to previous best
    ## Run 359 stress 0.3079616 
    ## Run 360 stress 9.061094e-05 
    ## ... Procrustes: rmse 0.0001755262  max resid 0.0003240661 
    ## ... Similar to previous best
    ## Run 361 stress 9.886371e-05 
    ## ... Procrustes: rmse 0.000174468  max resid 0.0003082057 
    ## ... Similar to previous best
    ## Run 362 stress 9.28179e-05 
    ## ... Procrustes: rmse 7.113551e-05  max resid 0.0001195465 
    ## ... Similar to previous best
    ## Run 363 stress 9.497156e-05 
    ## ... Procrustes: rmse 0.0001689423  max resid 0.0003028042 
    ## ... Similar to previous best
    ## Run 364 stress 9.88498e-05 
    ## ... Procrustes: rmse 8.978542e-05  max resid 0.0001487799 
    ## ... Similar to previous best
    ## Run 365 stress 9.942406e-05 
    ## ... Procrustes: rmse 0.0001691074  max resid 0.0002927035 
    ## ... Similar to previous best
    ## Run 366 stress 9.627935e-05 
    ## ... Procrustes: rmse 8.5466e-05  max resid 0.000142984 
    ## ... Similar to previous best
    ## Run 367 stress 9.686199e-05 
    ## ... Procrustes: rmse 8.641856e-05  max resid 0.0001488424 
    ## ... Similar to previous best
    ## Run 368 stress 9.591649e-05 
    ## ... Procrustes: rmse 0.0001591405  max resid 0.0002918179 
    ## ... Similar to previous best
    ## Run 369 stress 8.698687e-05 
    ## ... Procrustes: rmse 0.000176035  max resid 0.0003226143 
    ## ... Similar to previous best
    ## Run 370 stress 7.012606e-05 
    ## ... Procrustes: rmse 0.0001395174  max resid 0.0002705642 
    ## ... Similar to previous best
    ## Run 371 stress 9.400948e-05 
    ## ... Procrustes: rmse 0.0001684371  max resid 0.0003047258 
    ## ... Similar to previous best
    ## Run 372 stress 9.454318e-05 
    ## ... Procrustes: rmse 0.0001824663  max resid 0.0003423588 
    ## ... Similar to previous best
    ## Run 373 stress 8.548718e-05 
    ## ... Procrustes: rmse 5.333338e-05  max resid 9.269778e-05 
    ## ... Similar to previous best
    ## Run 374 stress 9.220066e-05 
    ## ... Procrustes: rmse 0.0001679164  max resid 0.0003160938 
    ## ... Similar to previous best
    ## Run 375 stress 9.338913e-05 
    ## ... Procrustes: rmse 0.0001247055  max resid 0.0002689405 
    ## ... Similar to previous best
    ## Run 376 stress 9.241746e-05 
    ## ... Procrustes: rmse 6.56172e-05  max resid 0.0001141526 
    ## ... Similar to previous best
    ## Run 377 stress 9.607555e-05 
    ## ... Procrustes: rmse 0.0001714789  max resid 0.0003219034 
    ## ... Similar to previous best
    ## Run 378 stress 9.932049e-05 
    ## ... Procrustes: rmse 0.000177507  max resid 0.000312145 
    ## ... Similar to previous best
    ## Run 379 stress 8.731456e-05 
    ## ... Procrustes: rmse 0.0001942452  max resid 0.0004733992 
    ## ... Similar to previous best
    ## Run 380 stress 9.941538e-05 
    ## ... Procrustes: rmse 0.0001761113  max resid 0.0003099076 
    ## ... Similar to previous best
    ## Run 381 stress 9.293389e-05 
    ## ... Procrustes: rmse 0.0001231125  max resid 0.0002644385 
    ## ... Similar to previous best
    ## Run 382 stress 9.070308e-05 
    ## ... Procrustes: rmse 0.0001903326  max resid 0.0002539198 
    ## ... Similar to previous best
    ## Run 383 stress 9.897047e-05 
    ## ... Procrustes: rmse 9.197533e-05  max resid 0.000154129 
    ## ... Similar to previous best
    ## Run 384 stress 8.69783e-05 
    ## ... Procrustes: rmse 0.0001761317  max resid 0.0003216916 
    ## ... Similar to previous best
    ## Run 385 stress 8.122612e-05 
    ## ... Procrustes: rmse 0.0001012802  max resid 0.0002074323 
    ## ... Similar to previous best
    ## Run 386 stress 8.801747e-05 
    ## ... Procrustes: rmse 0.0001825471  max resid 0.0004448109 
    ## ... Similar to previous best
    ## Run 387 stress 9.472226e-05 
    ## ... Procrustes: rmse 7.978617e-05  max resid 0.0001391557 
    ## ... Similar to previous best
    ## Run 388 stress 9.511745e-05 
    ## ... Procrustes: rmse 0.0001682905  max resid 0.0003020769 
    ## ... Similar to previous best
    ## Run 389 stress 9.189522e-05 
    ## ... Procrustes: rmse 0.0001219048  max resid 0.0002615424 
    ## ... Similar to previous best
    ## Run 390 stress 9.843923e-05 
    ## ... Procrustes: rmse 9.975892e-05  max resid 0.0002164723 
    ## ... Similar to previous best
    ## Run 391 stress 9.323696e-05 
    ## ... Procrustes: rmse 7.819608e-05  max resid 0.0001341211 
    ## ... Similar to previous best
    ## Run 392 stress 9.914147e-05 
    ## ... Procrustes: rmse 0.0001874025  max resid 0.0002961873 
    ## ... Similar to previous best
    ## Run 393 stress 9.617522e-05 
    ## ... Procrustes: rmse 0.0001630686  max resid 0.0003118473 
    ## ... Similar to previous best
    ## Run 394 stress 9.664129e-05 
    ## ... Procrustes: rmse 8.202165e-05  max resid 0.0001413767 
    ## ... Similar to previous best
    ## Run 395 stress 9.896568e-05 
    ## ... Procrustes: rmse 0.0001759966  max resid 0.0003107956 
    ## ... Similar to previous best
    ## Run 396 stress 9.982408e-05 
    ## ... Procrustes: rmse 9.380823e-05  max resid 0.0001629297 
    ## ... Similar to previous best
    ## Run 397 stress 9.681023e-05 
    ## ... Procrustes: rmse 8.562459e-05  max resid 0.000147508 
    ## ... Similar to previous best
    ## Run 398 stress 9.739886e-05 
    ## ... Procrustes: rmse 8.66158e-05  max resid 0.0001415316 
    ## ... Similar to previous best
    ## Run 399 stress 8.836452e-05 
    ## ... Procrustes: rmse 0.0001952917  max resid 0.0004756326 
    ## ... Similar to previous best
    ## Run 400 stress 9.292785e-05 
    ## ... Procrustes: rmse 7.006613e-05  max resid 0.0001222313 
    ## ... Similar to previous best
    ## Run 401 stress 9.57615e-05 
    ## ... Procrustes: rmse 0.0001011543  max resid 0.0002151641 
    ## ... Similar to previous best
    ## Run 402 stress 9.257287e-05 
    ## ... Procrustes: rmse 0.000158377  max resid 0.0002778158 
    ## ... Similar to previous best
    ## Run 403 stress 0.307194 
    ## Run 404 stress 9.863574e-05 
    ## ... Procrustes: rmse 9.016994e-05  max resid 0.0001512958 
    ## ... Similar to previous best
    ## Run 405 stress 9.113866e-05 
    ## ... Procrustes: rmse 0.0001200488  max resid 0.0002598591 
    ## ... Similar to previous best
    ## Run 406 stress 9.978089e-05 
    ## ... Procrustes: rmse 6.682927e-05  max resid 0.0001147773 
    ## ... Similar to previous best
    ## Run 407 stress 7.759542e-05 
    ## ... Procrustes: rmse 0.0001829068  max resid 0.0002450977 
    ## ... Similar to previous best
    ## Run 408 stress 9.902496e-05 
    ## ... Procrustes: rmse 0.000176211  max resid 0.0003259877 
    ## ... Similar to previous best
    ## Run 409 stress 9.157873e-05 
    ## ... Procrustes: rmse 6.741411e-05  max resid 0.00011502 
    ## ... Similar to previous best
    ## Run 410 stress 9.85611e-05 
    ## ... Procrustes: rmse 0.0001749796  max resid 0.0003081938 
    ## ... Similar to previous best
    ## Run 411 stress 9.872561e-05 
    ## ... Procrustes: rmse 9.038348e-05  max resid 0.000150542 
    ## ... Similar to previous best
    ## Run 412 stress 8.631885e-05 
    ## ... Procrustes: rmse 0.0001930309  max resid 0.0002546583 
    ## ... Similar to previous best
    ## Run 413 stress 9.828423e-05 
    ## ... Procrustes: rmse 0.0001719647  max resid 0.0002896394 
    ## ... Similar to previous best
    ## Run 414 stress 8.90029e-05 
    ## ... Procrustes: rmse 0.0001527418  max resid 0.0002590911 
    ## ... Similar to previous best
    ## Run 415 stress 9.380584e-05 
    ## ... Procrustes: rmse 0.0001695817  max resid 0.0003173055 
    ## ... Similar to previous best
    ## Run 416 stress 9.488867e-05 
    ## ... Procrustes: rmse 8.06404e-05  max resid 0.0001360556 
    ## ... Similar to previous best
    ## Run 417 stress 9.756673e-05 
    ## ... Procrustes: rmse 7.305202e-05  max resid 0.0001256796 
    ## ... Similar to previous best
    ## Run 418 stress 8.747751e-05 
    ## ... Procrustes: rmse 5.791707e-05  max resid 0.0001011238 
    ## ... Similar to previous best
    ## Run 419 stress 9.848698e-05 
    ## ... Procrustes: rmse 0.000125756  max resid 0.0002706814 
    ## ... Similar to previous best
    ## Run 420 stress 9.825234e-05 
    ## ... Procrustes: rmse 7.45114e-05  max resid 0.00012341 
    ## ... Similar to previous best
    ## Run 421 stress 9.365402e-05 
    ## ... Procrustes: rmse 0.0001868895  max resid 0.000454731 
    ## ... Similar to previous best
    ## Run 422 stress 7.446587e-05 
    ## ... Procrustes: rmse 2.157611e-05  max resid 3.412787e-05 
    ## ... Similar to previous best
    ## Run 423 stress 9.606057e-05 
    ## ... Procrustes: rmse 0.0001637274  max resid 0.000295356 
    ## ... Similar to previous best
    ## Run 424 stress 9.724342e-05 
    ## ... Procrustes: rmse 8.586385e-05  max resid 0.0001506855 
    ## ... Similar to previous best
    ## Run 425 stress 8.371397e-05 
    ## ... Procrustes: rmse 0.0001628578  max resid 0.0003513742 
    ## ... Similar to previous best
    ## Run 426 stress 9.605722e-05 
    ## ... Procrustes: rmse 0.0001681153  max resid 0.0002892489 
    ## ... Similar to previous best
    ## Run 427 stress 9.984802e-05 
    ## ... Procrustes: rmse 9.451719e-05  max resid 0.0001590527 
    ## ... Similar to previous best
    ## Run 428 stress 9.655492e-05 
    ## ... Procrustes: rmse 8.175684e-05  max resid 0.0001409668 
    ## ... Similar to previous best
    ## Run 429 stress 9.535908e-05 
    ## ... Procrustes: rmse 0.0001025728  max resid 0.000221267 
    ## ... Similar to previous best
    ## Run 430 stress 9.297601e-05 
    ## ... Procrustes: rmse 0.0001656718  max resid 0.0002984187 
    ## ... Similar to previous best
    ## Run 431 stress 9.861911e-05 
    ## ... Procrustes: rmse 0.0001295909  max resid 0.0002808373 
    ## ... Similar to previous best
    ## Run 432 stress 9.109844e-05 
    ## ... Procrustes: rmse 0.0001793829  max resid 0.0003293943 
    ## ... Similar to previous best
    ## Run 433 stress 9.930006e-05 
    ## ... Procrustes: rmse 0.0001273459  max resid 0.0002744936 
    ## ... Similar to previous best
    ## Run 434 stress 9.965106e-05 
    ## ... Procrustes: rmse 0.0001281305  max resid 0.0002766411 
    ## ... Similar to previous best
    ## Run 435 stress 9.533863e-05 
    ## ... Procrustes: rmse 7.81173e-05  max resid 0.0001407565 
    ## ... Similar to previous best
    ## Run 436 stress 9.499171e-05 
    ## ... Procrustes: rmse 7.889667e-05  max resid 0.0001379305 
    ## ... Similar to previous best
    ## Run 437 stress 9.393129e-05 
    ## ... Procrustes: rmse 7.791151e-05  max resid 0.0001345012 
    ## ... Similar to previous best
    ## Run 438 stress 9.709784e-05 
    ## ... Procrustes: rmse 8.666974e-05  max resid 0.0001446059 
    ## ... Similar to previous best
    ## Run 439 stress 9.477504e-05 
    ## ... Procrustes: rmse 7.692988e-05  max resid 0.0001257954 
    ## ... Similar to previous best
    ## Run 440 stress 8.751966e-05 
    ## ... Procrustes: rmse 0.0001866885  max resid 0.0002576763 
    ## ... Similar to previous best
    ## Run 441 stress 9.804833e-05 
    ## ... Procrustes: rmse 0.0002036943  max resid 0.000273424 
    ## ... Similar to previous best
    ## Run 442 stress 9.755203e-05 
    ## ... Procrustes: rmse 8.830429e-05  max resid 0.0001524919 
    ## ... Similar to previous best
    ## Run 443 stress 9.828551e-05 
    ## ... Procrustes: rmse 9.01596e-05  max resid 0.0001545989 
    ## ... Similar to previous best
    ## Run 444 stress 9.519177e-05 
    ## ... Procrustes: rmse 8.282004e-05  max resid 0.0001385738 
    ## ... Similar to previous best
    ## Run 445 stress 9.137101e-05 
    ## ... Procrustes: rmse 0.0001207394  max resid 0.0002600549 
    ## ... Similar to previous best
    ## Run 446 stress 9.954744e-05 
    ## ... Procrustes: rmse 8.80147e-05  max resid 0.0001455564 
    ## ... Similar to previous best
    ## Run 447 stress 9.179599e-05 
    ## ... Procrustes: rmse 0.0001203233  max resid 0.0002567293 
    ## ... Similar to previous best
    ## Run 448 stress 0.3080897 
    ## Run 449 stress 8.919683e-05 
    ## ... Procrustes: rmse 0.0001598333  max resid 0.0003002631 
    ## ... Similar to previous best
    ## Run 450 stress 9.611504e-05 
    ## ... Procrustes: rmse 0.000164105  max resid 0.0002959331 
    ## ... Similar to previous best
    ## Run 451 stress 9.380901e-05 
    ## ... Procrustes: rmse 5.783917e-05  max resid 9.098804e-05 
    ## ... Similar to previous best
    ## Run 452 stress 8.940515e-05 
    ## ... Procrustes: rmse 0.000163539  max resid 0.0003067604 
    ## ... Similar to previous best
    ## Run 453 stress 7.405641e-05 
    ## ... Procrustes: rmse 0.0001112317  max resid 0.0002294331 
    ## ... Similar to previous best
    ## Run 454 stress 9.966734e-05 
    ## ... Procrustes: rmse 0.0001919994  max resid 0.0004660022 
    ## ... Similar to previous best
    ## Run 455 stress 9.341296e-05 
    ## ... Procrustes: rmse 0.0001218285  max resid 0.0002608759 
    ## ... Similar to previous best
    ## Run 456 stress 9.357668e-05 
    ## ... Procrustes: rmse 0.0001234118  max resid 0.0002656123 
    ## ... Similar to previous best
    ## Run 457 stress 9.918398e-05 
    ## ... Procrustes: rmse 9.375371e-05  max resid 0.0001598001 
    ## ... Similar to previous best
    ## Run 458 stress 9.542731e-05 
    ## ... Procrustes: rmse 8.313164e-05  max resid 0.0001387316 
    ## ... Similar to previous best
    ## Run 459 stress 8.353206e-05 
    ## ... Procrustes: rmse 0.0001469103  max resid 0.0002608088 
    ## ... Similar to previous best
    ## Run 460 stress 0.3083098 
    ## Run 461 stress 9.939411e-05 
    ## ... Procrustes: rmse 8.580608e-05  max resid 0.0001497928 
    ## ... Similar to previous best
    ## Run 462 stress 8.739583e-05 
    ## ... Procrustes: rmse 0.0001601251  max resid 0.0002743127 
    ## ... Similar to previous best
    ## Run 463 stress 8.496826e-05 
    ## ... Procrustes: rmse 5.17374e-05  max resid 8.665148e-05 
    ## ... Similar to previous best
    ## Run 464 stress 9.346687e-05 
    ## ... Procrustes: rmse 0.0001665077  max resid 0.0003038296 
    ## ... Similar to previous best
    ## Run 465 stress 9.614609e-05 
    ## ... Procrustes: rmse 0.0001722332  max resid 0.0003094656 
    ## ... Similar to previous best
    ## Run 466 stress 9.652052e-05 
    ## ... Procrustes: rmse 8.130424e-05  max resid 0.0001332839 
    ## ... Similar to previous best
    ## Run 467 stress 9.932653e-05 
    ## ... Procrustes: rmse 8.712807e-05  max resid 0.0001509692 
    ## ... Similar to previous best
    ## Run 468 stress 9.781258e-05 
    ## ... Procrustes: rmse 0.0001839359  max resid 0.0003121793 
    ## ... Similar to previous best
    ## Run 469 stress 8.427373e-05 
    ## ... Procrustes: rmse 0.0001889805  max resid 0.0002492108 
    ## ... Similar to previous best
    ## Run 470 stress 8.973116e-05 
    ## ... Procrustes: rmse 0.0001557412  max resid 0.0002730229 
    ## ... Similar to previous best
    ## Run 471 stress 9.817456e-05 
    ## ... Procrustes: rmse 0.0001758697  max resid 0.0003130839 
    ## ... Similar to previous best
    ## Run 472 stress 9.807115e-05 
    ## ... Procrustes: rmse 0.0001856587  max resid 0.000430968 
    ## ... Similar to previous best
    ## Run 473 stress 9.68448e-05 
    ## ... Procrustes: rmse 0.0001936315  max resid 0.0002610127 
    ## ... Similar to previous best
    ## Run 474 stress 9.418207e-05 
    ## ... Procrustes: rmse 0.0001683562  max resid 0.0003049829 
    ## ... Similar to previous best
    ## Run 475 stress 9.400571e-05 
    ## ... Procrustes: rmse 0.0001239253  max resid 0.0002669587 
    ## ... Similar to previous best
    ## Run 476 stress 9.93546e-05 
    ## ... Procrustes: rmse 0.0001009863  max resid 0.0001886918 
    ## ... Similar to previous best
    ## Run 477 stress 9.78308e-05 
    ## ... Procrustes: rmse 0.0001750078  max resid 0.0003092405 
    ## ... Similar to previous best
    ## Run 478 stress 9.285121e-05 
    ## ... Procrustes: rmse 0.0001589223  max resid 0.0002898202 
    ## ... Similar to previous best
    ## Run 479 stress 8.752904e-05 
    ## ... Procrustes: rmse 6.303698e-05  max resid 0.0001072156 
    ## ... Similar to previous best
    ## Run 480 stress 9.610681e-05 
    ## ... Procrustes: rmse 6.937443e-05  max resid 0.0001206673 
    ## ... Similar to previous best
    ## Run 481 stress 9.439911e-05 
    ## ... Procrustes: rmse 0.0001237346  max resid 0.0002677856 
    ## ... Similar to previous best
    ## Run 482 stress 4.917448e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001671476  max resid 0.0002345768 
    ## ... Similar to previous best
    ## Run 483 stress 9.734539e-05 
    ## ... Procrustes: rmse 0.0002478525  max resid 0.0003719133 
    ## ... Similar to previous best
    ## Run 484 stress 0.2775741 
    ## Run 485 stress 8.706907e-05 
    ## ... Procrustes: rmse 0.0001665622  max resid 0.0002748563 
    ## ... Similar to previous best
    ## Run 486 stress 9.24769e-05 
    ## ... Procrustes: rmse 0.000216704  max resid 0.0003463352 
    ## ... Similar to previous best
    ## Run 487 stress 9.785534e-05 
    ## ... Procrustes: rmse 0.0002479837  max resid 0.0003711631 
    ## ... Similar to previous best
    ## Run 488 stress 8.385708e-05 
    ## ... Procrustes: rmse 0.0001099366  max resid 0.0001901917 
    ## ... Similar to previous best
    ## Run 489 stress 9.15516e-05 
    ## ... Procrustes: rmse 0.0002353393  max resid 0.0003441006 
    ## ... Similar to previous best
    ## Run 490 stress 9.948617e-05 
    ## ... Procrustes: rmse 0.0002543234  max resid 0.0003773411 
    ## ... Similar to previous best
    ## Run 491 stress 9.79822e-05 
    ## ... Procrustes: rmse 0.0002389282  max resid 0.0003447974 
    ## ... Similar to previous best
    ## Run 492 stress 9.69406e-05 
    ## ... Procrustes: rmse 0.0002185838  max resid 0.0003483493 
    ## ... Similar to previous best
    ## Run 493 stress 9.906968e-05 
    ## ... Procrustes: rmse 0.0002485389  max resid 0.0003724972 
    ## ... Similar to previous best
    ## Run 494 stress 9.412541e-05 
    ## ... Procrustes: rmse 0.0002207998  max resid 0.0003545279 
    ## ... Similar to previous best
    ## Run 495 stress 9.857836e-05 
    ## ... Procrustes: rmse 0.0002328416  max resid 0.0003771447 
    ## ... Similar to previous best
    ## Run 496 stress 9.881769e-05 
    ## ... Procrustes: rmse 0.0001809214  max resid 0.0003494159 
    ## ... Similar to previous best
    ## Run 497 stress 9.441824e-05 
    ## ... Procrustes: rmse 0.000171241  max resid 0.0003299917 
    ## ... Similar to previous best
    ## Run 498 stress 8.924747e-05 
    ## ... Procrustes: rmse 8.384263e-05  max resid 0.0001592128 
    ## ... Similar to previous best
    ## Run 499 stress 9.006759e-05 
    ## ... Procrustes: rmse 0.000121063  max resid 0.0001994826 
    ## ... Similar to previous best
    ## Run 500 stress 9.575065e-05 
    ## ... Procrustes: rmse 0.0002468521  max resid 0.0003628567 
    ## ... Similar to previous best
    ## *** Best solution repeated 18 times

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
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31833 0.699  0.014 *
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
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31833 0.699  0.014 *
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
    ## temperature_median -0.75817 -0.65206 0.0724  0.681  
    ## salinity_median     0.94798  0.31833 0.6990  0.016 *
    ## oxygen_median       0.57184 -0.82036 0.6836  0.123  
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
    ## temperature_median -0.75817 -0.65206 0.0724  1.000  
    ## salinity_median     0.94798  0.31833 0.6990  0.048 *
    ## oxygen_median       0.57184 -0.82036 0.6836  0.369  
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
    ## temperature_median -0.92664  0.37594 0.0954  0.745  
    ## salinity_median     0.90472  0.42601 0.6840  0.026 *
    ## oxygen_median       0.44477 -0.89565 0.5837  0.198  
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
    ## temperature_median -0.92664  0.37594 0.0954  1.000  
    ## salinity_median     0.90472  0.42601 0.6840  0.078 .
    ## oxygen_median       0.44477 -0.89565 0.5837  0.594  
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
    ## temperature_median -0.52520 -0.85098 0.0223  0.896  
    ## salinity_median    -0.79785 -0.60286 0.7269  0.017 *
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.013 *
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
    ## temperature_median -0.52520 -0.85098 0.0223  1.000  
    ## salinity_median    -0.79785 -0.60286 0.7269  0.051 .
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.039 *
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
    ##                          NMDS1       NMDS2     r2 Pr(>r)  
    ## temperature_median  0.00063832  1.00000000 0.0377  0.760  
    ## salinity_median    -0.00064573  1.00000000 0.6087  0.077 .
    ## oxygen_median      -0.00233995 -1.00000000 0.7257  0.448  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ##                          NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median  0.00063832  1.00000000 0.0377  1.000
    ## salinity_median    -0.00064573  1.00000000 0.6087  0.231
    ## oxygen_median      -0.00233995 -1.00000000 0.7257  1.000
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
    ##                         NMDS1      NMDS2     r2 Pr(>r)  
    ## temperature_median -0.0070403  0.9999800 0.0200  0.927  
    ## salinity_median     0.0022377  1.0000000 0.3815  0.306  
    ## oxygen_median       0.0004436 -1.0000000 0.6951  0.069 .
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
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median -0.0070403  0.9999800 0.0200  1.000
    ## salinity_median     0.0022377  1.0000000 0.3815  0.918
    ## oxygen_median       0.0004436 -1.0000000 0.6951  0.207
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
    ## temperature_median       -0.39652 -0.91803 0.0381  0.913  
    ## salinity_median          -0.96164  0.27433 0.5830  0.117  
    ## oxygen_median             0.83740  0.54660 0.7450  0.032 *
    ## distance_to_ocean_mean_m -0.09987  0.99500 0.0319  0.920  
    ## max_depth                -0.87447  0.48507 0.1022  0.759  
    ## logArea                  -0.40925  0.91242 0.2148  0.547  
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
    ## temperature_median       -0.39652 -0.91803 0.0381  1.000
    ## salinity_median          -0.96164  0.27433 0.5830  0.702
    ## oxygen_median             0.83740  0.54660 0.7450  0.192
    ## distance_to_ocean_mean_m -0.09987  0.99500 0.0319  1.000
    ## max_depth                -0.87447  0.48507 0.1022  1.000
    ## logArea                  -0.40925  0.91242 0.2148  1.000
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
    ## distance_to_ocean_min_m -0.82790  0.56088 0.5558  0.063 .
    ## max_depth               -0.20869 -0.97798 0.0858  0.892  
    ## logArea                  0.25301 -0.96746 0.2012  0.290  
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
    ## distance_to_ocean_min_m -0.82790  0.56088 0.5558  0.189
    ## max_depth               -0.20869 -0.97798 0.0858  1.000
    ## logArea                  0.25301 -0.96746 0.2012  0.870
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
    ## distance_to_ocean_min_m -0.82790  0.56088 0.5558  0.075 .
    ## max_depth               -0.20869 -0.97798 0.0858  0.889  
    ## logArea                  0.25301 -0.96746 0.2012  0.296  
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
    ## distance_to_ocean_min_m -0.82790  0.56088 0.5558  0.225
    ## max_depth               -0.20869 -0.97798 0.0858  1.000
    ## logArea                  0.25301 -0.96746 0.2012  0.888
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
    ## distance_to_ocean_min_m -0.76293  0.64649 0.4672  0.146
    ## max_depth               -0.34720 -0.93779 0.0777  0.962
    ## logArea                  0.55379  0.83265 0.0074  0.968
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
    ## distance_to_ocean_min_m -0.76293  0.64649 0.4672  0.438
    ## max_depth               -0.34720 -0.93779 0.0777  1.000
    ## logArea                  0.55379  0.83265 0.0074  1.000
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
    ## distance_to_ocean_min_m  0.45616 -0.88990 0.4273  0.179   
    ## max_depth               -0.94119 -0.33789 0.6682  0.003 **
    ## logArea                 -0.86397  0.50354 0.2339  0.460   
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
    ## distance_to_ocean_min_m  0.45616 -0.88990 0.4273  0.537   
    ## max_depth               -0.94119 -0.33789 0.6682  0.009 **
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
    ## distance_to_ocean_min_m  0.99474  0.10244 0.6338  0.212
    ## max_depth                0.73014  0.68330 0.0374  0.972
    ## logArea                 -0.39136  0.92024 0.3079  0.311
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
    ## distance_to_ocean_min_m  0.99474  0.10244 0.6338  0.636
    ## max_depth                0.73014  0.68330 0.0374  1.000
    ## logArea                 -0.39136  0.92024 0.3079  0.933
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
    ## distance_to_ocean_min_m -0.00014803 -1.00000000 0.6528  0.074 .
    ## max_depth                0.00085537  1.00000000 0.7262  0.070 .
    ## logArea                  0.00033856 -1.00000000 0.7875  0.012 *
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
    ## distance_to_ocean_min_m -0.00014803 -1.00000000 0.6528  0.222  
    ## max_depth                0.00085537  1.00000000 0.7262  0.210  
    ## logArea                  0.00033856 -1.00000000 0.7875  0.036 *
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
    ##       Significance: 0.212 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.210 0.241 0.258 0.298 
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
    ##       Significance: 0.006 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.437 0.465 0.488 0.501 
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
    ##       Significance: 0.258 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.522 0.558 0.588 0.616 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_mant_pv <- rbind(PD_beta_env_mant_t$signif, PD_beta_env_mant_s$signif, PD_beta_env_mant_o$signif)
PD_beta_env_mant_pv <- PD_beta_env_mant_pv[,1]
(PD_beta_env_mant_pv <- p.adjust(PD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.636 0.018 0.774

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
    ##       Significance: 0.233 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.211 0.256 0.284 0.322 
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
    ##       Significance: 0.011 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.378 0.413 0.437 0.463 
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
    ##       Significance: 0.324 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.418 0.462 0.501 0.548 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_MS_mant_pv <- rbind(PD_beta_env_MS_mant_t$signif, PD_beta_env_MS_mant_s$signif, PD_beta_env_MS_mant_o$signif)
PD_beta_env_MS_mant_pv <- PD_beta_env_MS_mant_pv[,1]
(PD_beta_env_MS_mant_pv <- p.adjust(PD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.699 0.033 0.972

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
    ##       Significance: 0.54 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.130 0.184 0.242 0.351 
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
    ##       Significance: 0.022 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.298 0.380 0.440 0.483 
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
    ##       Significance: 0.014 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.322 0.400 0.498 0.588 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_OM_mant_pv <- rbind(PD_beta_env_OM_mant_t$signif, PD_beta_env_OM_mant_s$signif, PD_beta_env_OM_mant_o$signif)
PD_beta_env_OM_mant_pv <- PD_beta_env_OM_mant_pv[,1]
(PD_beta_env_OM_mant_pv <- p.adjust(PD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.066 0.042

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
    ##       Significance: 0.446 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0429 0.0858 0.1120 0.1340 
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
    ##       Significance: 0.027 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.392 0.441 0.478 0.507 
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
    ##       Significance: 0.233 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.693 0.718 0.735 0.757 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_SO_mant_pv <- rbind(PD_beta_env_SO_mant_t$signif, PD_beta_env_SO_mant_s$signif, PD_beta_env_SO_mant_o$signif)
PD_beta_env_SO_mant_pv <- PD_beta_env_SO_mant_pv[,1]
(PD_beta_env_SO_mant_pv <- p.adjust(PD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.081 0.699

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
    ##       Significance: 0.76 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.289 0.387 0.500 0.731 
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
    ##       Significance: 0.244 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.209 0.264 0.298 0.349 
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
    ##       Significance: 0.094 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.261 0.367 0.450 0.533 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_M_mant_pv <- rbind(PD_beta_env_M_mant_t$signif, PD_beta_env_M_mant_s$signif, PD_beta_env_M_mant_o$signif)
PD_beta_env_M_mant_pv <- PD_beta_env_M_mant_pv[,1]
(PD_beta_env_M_mant_pv <- p.adjust(PD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.732 0.282

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
    ##       Significance: 0.259 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.349 0.439 0.492 0.573 
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
    ##       Significance: 0.569 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.456 0.479 0.504 0.526 
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
    ##       Significance: 0.779 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.161 0.190 0.214 0.237 
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
    ##       Significance: 0.731 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.111 0.133 0.146 0.155 
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
    ##       Significance: 0.642 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.308 0.344 0.380 0.404 
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
    ##       Significance: 0.908 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.169 0.201 0.239 0.272 
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
    ##       Significance: 0.882 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.130 0.156 0.177 0.192 
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
    ##       Significance: 0.425 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.174 0.206 0.231 0.276 
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
    ##       Significance: 0.017 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.165 0.238 0.285 0.350 
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
    ##       Significance: 0.054 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.167 0.197 0.217 0.270 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_OM_mant_pv <- rbind( PD_beta_geo_OM_mant_dmean$signif, PD_beta_geo_OM_mant_md$signif, PD_beta_geo_OM_mant_la$signif)
PD_beta_geo_OM_mant_pv <- PD_beta_geo_OM_mant_pv[,1]
(PD_beta_geo_OM_mant_pv <- p.adjust(PD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.051 0.162

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
    ##       Significance: 0.555 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.623 0.649 0.671 0.689 
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
    ##       Significance: 0.795 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.115 0.142 0.169 0.198 
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
    ##       Significance: 0.213 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.372 0.395 0.413 0.427 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_SO_mant_pv <- rbind( PD_beta_geo_SO_mant_dmean$signif, PD_beta_geo_SO_mant_md$signif, PD_beta_geo_SO_mant_la$signif)
PD_beta_geo_SO_mant_pv <- PD_beta_geo_SO_mant_pv[,1]
(PD_beta_geo_SO_mant_pv <- p.adjust(PD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.639

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
    ##       Significance: 0.748 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.258 0.336 0.408 0.664 
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
    ##       Significance: 0.015 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.274 0.410 0.541 0.630 
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
    ##       Significance: 0.013 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.202 0.267 0.333 0.434 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_M_mant_pv <- rbind( PD_beta_geo_M_mant_dmean$signif, PD_beta_geo_M_mant_md$signif, PD_beta_geo_M_mant_la$signif)
PD_beta_geo_M_mant_pv <- PD_beta_geo_M_mant_pv[,1]
(PD_beta_geo_M_mant_pv <- p.adjust(PD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.045 0.039

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
        legend.position = "right", 
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
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_beta_ref_NMDS.jpg", PD_beta_ref_plot, width = 6.26, height = 6, units = "in")
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
