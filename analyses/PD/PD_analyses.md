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

    ##   |                  |          |   0%  |                  |          |   3%                                                                          |                  |.         |   7% [Bringing everything together load modifying files packages]             |                  |.         |  10%                                                                          |                  |.         |  13% [Bringing everything together load in modifying files]                   |                  |..        |  17%                                                                          |                  |..        |  20% [Check species names across files]                                       |                  |..        |  23%                                                                          |                  |...       |  27% [Modify environment data]                                                |                  |...       |  30%                                                                          |                  |...       |  33% [Modify incidence matrices]                                              |                  |....      |  37%                                                                          |                  |....      |  40% [Modify phylogeny]                                                       |                  |....      |  43%                                                                          |                  |.....     |  47% [Modify trait data]                                                      |                  |.....     |  50%                                                                          |                  |.....     |  53% [Modify Site_type data frames]                                           |                  |......    |  57%                                                                          |                  |......    |  60% [Modify Site_type trait data]                                            |                  |......    |  63%                                                                          |                  |.......   |  67% [Site_type trait data tests]                                             |                  |.......   |  70%                                                                          |                  |.......   |  73% [Modify site trait data frames]                                          |                  |........  |  77%                                                                          |                  |........  |  80% [Site trait data tests]                                                  |                  |........  |  83%                                                                          |                  |......... |  87% [unnamed-chunk-5]                                                        |                  |......... |  90%                                                                          |                  |......... |  93% [unnamed-chunk-6]                                                        |                  |..........|  97%                                                                          |                  |..........| 100% [Bringing everything together load out modified files and session info]

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

Site_type_group_ref <- env[,19]
Site_type_group <- env[surveyed_sites,19]
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
tree_disp_env$Site_type <- factor(tree_disp_env$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
row.names(tree_disp_env) <- tree_disp_env$X

# Mean pairwise distances
tree_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_sesmpd.csv")
tree_sesmpd_env <- merge(tree_sesmpd, env, by = "X", sort = F)
tree_sesmpd_env$measure <- "tree_sesmpd"
tree_sesmpd_env$Site_type <- factor(tree_sesmpd_env$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
row.names(tree_sesmpd_env) <- tree_sesmpd_env$X

# Mean nearest taxon distances
tree_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_sesmntd.csv")
tree_sesmntd_env <- merge(tree_sesmntd, env, by = "X", sort = F)
tree_sesmntd_env$measure <- "tree_sesmntd"
tree_sesmntd_env$Site_type <- factor(tree_sesmntd_env$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
row.names(tree_sesmntd_env) <- tree_sesmntd_env$X

#Faith's PD
tree_sespd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/tree_sespd.csv")
tree_sespd_env <- merge(tree_sespd, env, by = "X", sort = F)
tree_sespd_env$measure <- "tree_sespd"
tree_sespd_env$Site_type <- factor(tree_sespd_env$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
stree_sespd_env <- merge(tree_sespd_env, tree_disp, by = "X", sort = F)
row.names(tree_sespd_env) <- tree_sespd_env$X


## Surveyed sites
# Dispersion 
stree_disp <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_disp.csv")
stree_disp_env <- merge(stree_disp, env, by = "X", sort = F)
stree_disp_env$measure <- "stree_disp"
stree_disp_env$Site_type <- factor(stree_disp_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))
row.names(stree_disp_env) <- stree_disp_env$X

# Mean pairwise distances
stree_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_sesmpd.csv")
stree_sesmpd_env <- merge(stree_sesmpd, env, by = "X", sort = F)
stree_sesmpd_env$measure <- "stree_sesmpd"
stree_sesmpd_env$Site_type <- factor(stree_sesmpd_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))
row.names(stree_sesmpd_env) <- stree_sesmpd_env$X

# Mean nearest taxon distances
stree_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_sesmntd.csv")
stree_sesmntd_env <- merge(stree_sesmntd, env, by = "X", sort = F)
stree_sesmntd_env$measure <- "stree_sesmntd"
stree_sesmntd_env$Site_type <- factor(stree_sesmntd_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))
row.names(stree_sesmntd_env) <- stree_sesmntd_env$X

#Faith's PD
stree_sespd <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/PD/stree_sespd.csv")
stree_sespd_env <- merge(stree_sespd, env, by = "X", sort = F)
stree_sespd_env$measure <- "stree_sespd"
stree_sespd_env$Site_type <- factor(stree_sespd_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))
stree_sespd_env <- merge(stree_sespd_env, stree_disp, by = "X", sort = F)
row.names(stree_sespd_env) <- stree_sespd_env$X
```

### PD alpha outliers for each site type

``` r
outlier_PD_alpha <- stree_sespd_env %>%
  group_by(Site_type) %>%
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
(outlier_PD_alpha_plot <- ggplot(outlier_PD_alpha, aes(x = Site_type, y = pd.obs.z, fill = Site_type)) +
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
  labs(y = "PD alpha z-scores", x = "Site type", color = "Outlier:", tag = "b")) +
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

mod <-  lm(pd.obs ~ salinity_median + oxygen_median + temperature_median + pH_median, stree_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ salinity_median + oxygen_median + temperature_median + 
    ##     pH_median, data = stree_sespd_env)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -960.6 -344.4    5.4  266.2 1158.6 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -23046.44    7912.77  -2.913   0.0114 *
    ## salinity_median        29.70      77.33   0.384   0.7067  
    ## oxygen_median         187.18     209.39   0.894   0.3865  
    ## temperature_median   -226.09     151.83  -1.489   0.1586  
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
    ## salinity_median     1 19125479 19125479 46.9646 7.895e-06 ***
    ## oxygen_median       1  5747913  5747913 14.1146  0.002124 ** 
    ## temperature_median  1    40929    40929  0.1005  0.755903    
    ## pH_median           1  3432036  3432036  8.4277  0.011572 *  
    ## Residuals          14  5701242   407232                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
rms::vif(mod)
```

    ##    salinity_median      oxygen_median temperature_median          pH_median 
    ##           3.435415           2.209454           1.385456           4.821892

``` r
car::vif(mod)
```

    ##    salinity_median      oxygen_median temperature_median          pH_median 
    ##           3.435415           2.209454           1.385456           4.821892

``` r
mod <-  lm(pd.obs ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1060.5  -388.2  -180.5   295.4  1777.7 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        -5063.68    6020.36  -0.841  0.41350   
    ## salinity_median      203.57      59.81   3.403  0.00393 **
    ## oxygen_median        588.26     192.40   3.057  0.00798 **
    ## temperature_median   -43.83     169.04  -0.259  0.79895   
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

    ##    salinity_median      oxygen_median temperature_median 
    ##           1.374749           1.247588           1.148556

``` r
car::vif(mod)
```

    ##    salinity_median      oxygen_median temperature_median 
    ##           1.374749           1.247588           1.148556

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
### Site_type
# Run ANOVA on the original data
anova_PRic_result <- aov(pd.obs ~ Site_type, data = stree_sespd_env[surveyed_sites,])
summary(anova_PRic_result)
```

    ##             Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type    2 23597133 11798567   17.19 5.48e-05 ***
    ## Residuals   19 13044225   686538                     
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
    ## Fit: aov(formula = pd.obs ~ Site_type, data = stree_sespd_env[surveyed_sites, ])
    ## 
    ## $Site_type
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
Site_type <- tukey_PRic_result$Site_type

# q-values
(OM_q <- (abs(Site_type[1,1]))/USE)
```

    ## [1] 0.6915491

``` r
(MS_q <- (abs(Site_type[3,1]))/BSE)
```

    ## [1] 7.643793

``` r
(SO_q <- (abs(Site_type[2,1]))/USE)
```

    ## [1] 6.385228

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx

# Run ANOVA on the original data
anova_zscore_result <- aov(pd.obs.z ~ Site_type, data = stree_sespd_env[surveyed_sites,])
summary(anova_zscore_result)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2  1.838   0.919   0.888  0.428
    ## Residuals   19 19.668   1.035

``` r
# Perform Tukey's HSD test
tukey_zscore_result <- TukeyHSD(anova_zscore_result)
print(tukey_zscore_result)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = pd.obs.z ~ Site_type, data = stree_sespd_env[surveyed_sites, ])
    ## 
    ## $Site_type
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
Site_type <- tukey_zscore_result$Site_type

# q-values
(OM_q <- (abs(Site_type[1,1]))/USE)
```

    ## [1] 1.670047

``` r
(MS_q <- (abs(Site_type[3,1]))/BSE)
```

    ## [1] 0.0007977355

``` r
(SO_q <- (abs(Site_type[2,1]))/USE)
```

    ## [1] 1.670785

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx

# Run ANOVA on the original data
anova_PDisp_result <- aov(Dispersion ~ Site_type, data = stree_sespd_env[surveyed_sites,])
summary(anova_PDisp_result)
```

    ##             Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type    2 0.000707 0.0003537   0.703  0.507
    ## Residuals   19 0.009556 0.0005030

``` r
# Perform Tukey's HSD test
tukey_PDisp_result <- TukeyHSD(anova_PDisp_result)
print(tukey_PDisp_result)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = Dispersion ~ Site_type, data = stree_sespd_env[surveyed_sites, ])
    ## 
    ## $Site_type
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
Site_type <- tukey_PDisp_result$Site_type

# q-values
(OM_q <- (abs(Site_type[1,1]))/USE)
```

    ## [1] 0.825711

``` r
(MS_q <- (abs(Site_type[3,1]))/BSE)
```

    ## [1] 1.676295

``` r
(SO_q <- (abs(Site_type[2,1]))/USE)
```

    ## [1] 0.7262367

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx
```

### PD alpha with env and geo linear models & ANOVAs

``` r
### Environmental
# Surveyed sites
PD_PRic_env_am <- aov(log(pd.obs) ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_am)
```

    ##                              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median               1  6.080   6.080  56.580 6.79e-05 ***
    ## oxygen_median                 1  0.984   0.984   9.161   0.0164 *  
    ## temperature_median            1  0.002   0.002   0.020   0.8918    
    ## Site_type                     2  0.136   0.068   0.631   0.5565    
    ## salinity_median:Site_type     2  0.117   0.059   0.545   0.6000    
    ## oxygen_median:Site_type       2  0.083   0.041   0.385   0.6925    
    ## temperature_median:Site_type  1  0.016   0.016   0.146   0.7127    
    ## Residuals                     8  0.860   0.107                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_env_lm <- lm(log(pd.obs) ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = stree_sespd_env[surveyed_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.34125 -0.17839 -0.09832  0.17872  0.51241 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        2.210239   2.191942   1.008  0.32928    
    ## salinity_median    0.130178   0.021778   5.978 2.53e-05 ***
    ## oxygen_median      0.244896   0.070051   3.496  0.00325 ** 
    ## temperature_median 0.009972   0.061544   0.162  0.87345    
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

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##       1.0000000000       0.0001013439       0.0130007650       1.0000000000

``` r
PD_z_env_am <- aov(pd.obs.z ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_am)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median               1  1.819  1.8189   2.643 0.1427  
    ## oxygen_median                 1  0.693  0.6934   1.007 0.3449  
    ## temperature_median            1  0.055  0.0548   0.080 0.7850  
    ## Site_type                     2  4.877  2.4384   3.543 0.0791 .
    ## salinity_median:Site_type     2  1.536  0.7682   1.116 0.3736  
    ## oxygen_median:Site_type       2  4.134  2.0669   3.003 0.1064  
    ## temperature_median:Site_type  1  0.030  0.0303   0.044 0.8390  
    ## Residuals                     8  5.506  0.6882                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_env_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_env_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[surveyed_sites_env, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.0491 -0.6604  0.2329  0.6865  1.6970 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.05264    7.98903  -0.007    0.995
    ## salinity_median    -0.05571    0.07937  -0.702    0.493
    ## oxygen_median      -0.20248    0.25532  -0.793    0.440
    ## temperature_median  0.05069    0.22431   0.226    0.824
    ## 
    ## Residual standard error: 1.035 on 15 degrees of freedom
    ## Multiple R-squared:  0.1376, Adjusted R-squared:  -0.03483 
    ## F-statistic: 0.798 on 3 and 15 DF,  p-value: 0.514

``` r
p_values <- summary(PD_z_env_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_am <- aov(Dispersion ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_am)
```

    ##                              Df   Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median               1 0.000811 0.0008107   1.233  0.299
    ## oxygen_median                 1 0.000009 0.0000086   0.013  0.912
    ## temperature_median            1 0.002013 0.0020129   3.060  0.118
    ## Site_type                     2 0.000568 0.0002840   0.432  0.664
    ## salinity_median:Site_type     2 0.000028 0.0000138   0.021  0.979
    ## oxygen_median:Site_type       2 0.001112 0.0005561   0.846  0.464
    ## temperature_median:Site_type  1 0.000152 0.0001523   0.232  0.643
    ## Residuals                     8 0.005262 0.0006578

``` r
p_values <- summary(PD_PDisp_env_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_env_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[surveyed_sites_env, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.032704 -0.010212 -0.000927  0.012572  0.048834 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -0.2235120  0.1681181  -1.329   0.2036  
    ## salinity_median     0.0030130  0.0016703   1.804   0.0914 .
    ## oxygen_median      -0.0002191  0.0053728  -0.041   0.9680  
    ## temperature_median  0.0097192  0.0047203   2.059   0.0573 .
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

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          0.8142048          0.3655132          1.0000000          0.2291871

``` r
# Mixed and stratified lakes
PD_PRic_env_MS_am <- aov(log(pd.obs) ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_MS_am)
```

    ##                              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median               1  6.080   6.080  56.580 6.79e-05 ***
    ## oxygen_median                 1  0.984   0.984   9.161   0.0164 *  
    ## temperature_median            1  0.002   0.002   0.020   0.8918    
    ## Site_type                     2  0.136   0.068   0.631   0.5565    
    ## salinity_median:Site_type     2  0.117   0.059   0.545   0.6000    
    ## oxygen_median:Site_type       2  0.083   0.041   0.385   0.6925    
    ## temperature_median:Site_type  1  0.016   0.016   0.146   0.7127    
    ## Residuals                     8  0.860   0.107                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_env_MS_lm <- lm(log(pd.obs) ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_env_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = stree_sespd_env[mixed_stratified_lakes, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.40218 -0.20251 -0.06598  0.23118  0.45102 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         1.72817    2.52835   0.684 0.507264    
    ## salinity_median     0.13399    0.02390   5.606 0.000115 ***
    ## oxygen_median       0.28179    0.08638   3.262 0.006801 ** 
    ## temperature_median  0.01795    0.06991   0.257 0.801749    
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

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##       1.0000000000       0.0004600793       0.0272028786       1.0000000000

``` r
PD_z_env_MS_am <- aov(pd.obs.z ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_MS_am)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median               1  1.819  1.8189   2.643 0.1427  
    ## oxygen_median                 1  0.693  0.6934   1.007 0.3449  
    ## temperature_median            1  0.055  0.0548   0.080 0.7850  
    ## Site_type                     2  4.877  2.4384   3.543 0.0791 .
    ## salinity_median:Site_type     2  1.536  0.7682   1.116 0.3736  
    ## oxygen_median:Site_type       2  4.134  2.0669   3.003 0.1064  
    ## temperature_median:Site_type  1  0.030  0.0303   0.044 0.8390  
    ## Residuals                     8  5.506  0.6882                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_env_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_env_MS_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_z_env_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[mixed_stratified_lakes, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.40930 -0.48650  0.04376  0.57984  1.20044 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -7.01848    6.96162  -1.008    0.333
    ## salinity_median    -0.02443    0.06581  -0.371    0.717
    ## oxygen_median       0.13066    0.23784   0.549    0.593
    ## temperature_median  0.20879    0.19249   1.085    0.299
    ## 
    ## Residual standard error: 0.8432 on 12 degrees of freedom
    ## Multiple R-squared:  0.1298, Adjusted R-squared:  -0.08774 
    ## F-statistic: 0.5967 on 3 and 12 DF,  p-value: 0.6292

``` r
p_values <- summary(PD_z_env_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_MS_am <- aov(Dispersion ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_MS_am)
```

    ##                              Df   Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median               1 0.000811 0.0008107   1.233  0.299
    ## oxygen_median                 1 0.000009 0.0000086   0.013  0.912
    ## temperature_median            1 0.002013 0.0020129   3.060  0.118
    ## Site_type                     2 0.000568 0.0002840   0.432  0.664
    ## salinity_median:Site_type     2 0.000028 0.0000138   0.021  0.979
    ## oxygen_median:Site_type       2 0.001112 0.0005561   0.846  0.464
    ## temperature_median:Site_type  1 0.000152 0.0001523   0.232  0.643
    ## Residuals                     8 0.005262 0.0006578

``` r
p_values <- summary(PD_PDisp_env_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_env_MS_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_env_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[mixed_stratified_lakes, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.026782 -0.011687 -0.000974  0.010088  0.046961 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -0.325673   0.182766  -1.782   0.1001  
    ## salinity_median     0.003451   0.001728   1.997   0.0690 .
    ## oxygen_median       0.004063   0.006244   0.651   0.5275  
    ## temperature_median  0.012123   0.005054   2.399   0.0336 *
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

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          0.4002584          0.2759753          1.0000000          0.1343561

``` r
# Ocean sites and mixed lakes
PD_PRic_env_OM_am <- aov(log(pd.obs) ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_OM_am)
```

    ##                              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median               1  6.080   6.080  56.580 6.79e-05 ***
    ## oxygen_median                 1  0.984   0.984   9.161   0.0164 *  
    ## temperature_median            1  0.002   0.002   0.020   0.8918    
    ## Site_type                     2  0.136   0.068   0.631   0.5565    
    ## salinity_median:Site_type     2  0.117   0.059   0.545   0.6000    
    ## oxygen_median:Site_type       2  0.083   0.041   0.385   0.6925    
    ## temperature_median:Site_type  1  0.016   0.016   0.146   0.7127    
    ## Residuals                     8  0.860   0.107                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_env_OM_lm <- lm(log(pd.obs) ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PRic_env_OM_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = stree_sespd_env[ocean_mixed_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.34221 -0.18181  0.00482  0.10561  0.49457 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         3.67640   12.32519   0.298    0.774
    ## salinity_median     0.16910    0.31857   0.531    0.612
    ## oxygen_median       0.20225    0.25697   0.787    0.457
    ## temperature_median -0.07279    0.17016  -0.428    0.682
    ## 
    ## Residual standard error: 0.3262 on 7 degrees of freedom
    ## Multiple R-squared:  0.2849, Adjusted R-squared:  -0.02157 
    ## F-statistic: 0.9296 on 3 and 7 DF,  p-value: 0.4751

``` r
p_values <- summary(PD_PRic_env_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
PD_z_env_OM_am <- aov(pd.obs.z ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_OM_am)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median               1  1.819  1.8189   2.643 0.1427  
    ## oxygen_median                 1  0.693  0.6934   1.007 0.3449  
    ## temperature_median            1  0.055  0.0548   0.080 0.7850  
    ## Site_type                     2  4.877  2.4384   3.543 0.0791 .
    ## salinity_median:Site_type     2  1.536  0.7682   1.116 0.3736  
    ## oxygen_median:Site_type       2  4.134  2.0669   3.003 0.1064  
    ## temperature_median:Site_type  1  0.030  0.0303   0.044 0.8390  
    ## Residuals                     8  5.506  0.6882                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_env_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_env_OM_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_z_env_OM_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[ocean_mixed_sites_env, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.8318 -0.4833 -0.2593  0.5477  1.0164 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        31.53378   30.82399   1.023    0.340
    ## salinity_median    -0.77195    0.79671  -0.969    0.365
    ## oxygen_median      -1.15354    0.64267  -1.795    0.116
    ## temperature_median -0.05226    0.42556  -0.123    0.906
    ## 
    ## Residual standard error: 0.8159 on 7 degrees of freedom
    ## Multiple R-squared:  0.641,  Adjusted R-squared:  0.4872 
    ## F-statistic: 4.167 on 3 and 7 DF,  p-value: 0.05472

``` r
p_values <- summary(PD_z_env_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          1.0000000          0.4629584          1.0000000

``` r
PD_PDisp_env_OM_am <- aov(Dispersion ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_OM_am)
```

    ##                              Df   Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median               1 0.000811 0.0008107   1.233  0.299
    ## oxygen_median                 1 0.000009 0.0000086   0.013  0.912
    ## temperature_median            1 0.002013 0.0020129   3.060  0.118
    ## Site_type                     2 0.000568 0.0002840   0.432  0.664
    ## salinity_median:Site_type     2 0.000028 0.0000138   0.021  0.979
    ## oxygen_median:Site_type       2 0.001112 0.0005561   0.846  0.464
    ## temperature_median:Site_type  1 0.000152 0.0001523   0.232  0.643
    ## Residuals                     8 0.005262 0.0006578

``` r
p_values <- summary(PD_PDisp_env_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_env_OM_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PDisp_env_OM_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[ocean_mixed_sites_env, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.009814 -0.004982 -0.001190  0.001979  0.019442 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         2.457e-01  3.643e-01   0.674    0.522
    ## salinity_median    -8.417e-04  9.417e-03  -0.089    0.931
    ## oxygen_median      -9.043e-03  7.596e-03  -1.190    0.273
    ## temperature_median -5.727e-05  5.030e-03  -0.011    0.991
    ## 
    ## Residual standard error: 0.009644 on 7 degrees of freedom
    ## Multiple R-squared:  0.3056, Adjusted R-squared:  0.007994 
    ## F-statistic: 1.027 on 3 and 7 DF,  p-value: 0.4369

``` r
p_values <- summary(PD_PDisp_env_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
# Stratified lakes and ocean sites
PD_PRic_env_SO_am <- aov(log(pd.obs) ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_env_SO_am)
```

    ##                              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median               1  6.080   6.080  56.580 6.79e-05 ***
    ## oxygen_median                 1  0.984   0.984   9.161   0.0164 *  
    ## temperature_median            1  0.002   0.002   0.020   0.8918    
    ## Site_type                     2  0.136   0.068   0.631   0.5565    
    ## salinity_median:Site_type     2  0.117   0.059   0.545   0.6000    
    ## oxygen_median:Site_type       2  0.083   0.041   0.385   0.6925    
    ## temperature_median:Site_type  1  0.016   0.016   0.146   0.7127    
    ## Residuals                     8  0.860   0.107                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_env_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PRic_env_SO_lm <- lm(log(pd.obs) ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PRic_env_SO_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = stree_sespd_env[ocean_stratified_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.33792 -0.07625  0.00662  0.06006  0.33720 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         1.68269    1.85239   0.908 0.393871    
    ## salinity_median     0.11951    0.01902   6.283 0.000411 ***
    ## oxygen_median       0.21178    0.06018   3.519 0.009743 ** 
    ## temperature_median  0.03913    0.05299   0.738 0.484238    
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

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##        1.000000000        0.001642753        0.038971567        1.000000000

``` r
PD_z_env_SO_am <- aov(pd.obs.z ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_z_env_SO_am)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median               1  1.819  1.8189   2.643 0.1427  
    ## oxygen_median                 1  0.693  0.6934   1.007 0.3449  
    ## temperature_median            1  0.055  0.0548   0.080 0.7850  
    ## Site_type                     2  4.877  2.4384   3.543 0.0791 .
    ## salinity_median:Site_type     2  1.536  0.7682   1.116 0.3736  
    ## oxygen_median:Site_type       2  4.134  2.0669   3.003 0.1064  
    ## temperature_median:Site_type  1  0.030  0.0303   0.044 0.8390  
    ## Residuals                     8  5.506  0.6882                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_env_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_z_env_SO_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_z_env_SO_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[ocean_stratified_sites_env, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.2602 -0.6625 -0.1193  0.7336  1.1685 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         1.47357    8.55414   0.172    0.868
    ## salinity_median    -0.13008    0.08783  -1.481    0.182
    ## oxygen_median      -0.26149    0.27793  -0.941    0.378
    ## temperature_median  0.06716    0.24472   0.274    0.792
    ## 
    ## Residual standard error: 1.048 on 7 degrees of freedom
    ## Multiple R-squared:  0.4059, Adjusted R-squared:  0.1514 
    ## F-statistic: 1.594 on 3 and 7 DF,  p-value: 0.2747

``` r
p_values <- summary(PD_z_env_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          0.7285996          1.0000000          1.0000000

``` r
PD_PDisp_env_SO_am <- aov(Dispersion ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_env_SO_am)
```

    ##                              Df   Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median               1 0.000811 0.0008107   1.233  0.299
    ## oxygen_median                 1 0.000009 0.0000086   0.013  0.912
    ## temperature_median            1 0.002013 0.0020129   3.060  0.118
    ## Site_type                     2 0.000568 0.0002840   0.432  0.664
    ## salinity_median:Site_type     2 0.000028 0.0000138   0.021  0.979
    ## oxygen_median:Site_type       2 0.001112 0.0005561   0.846  0.464
    ## temperature_median:Site_type  1 0.000152 0.0001523   0.232  0.643
    ## Residuals                     8 0.005262 0.0006578

``` r
p_values <- summary(PD_PDisp_env_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
PD_PDisp_env_SO_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PDisp_env_SO_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[ocean_stratified_sites_env, ])
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.027570 -0.022235  0.006409  0.012965  0.050260 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.249416   0.239944  -1.039    0.333
    ## salinity_median     0.002342   0.002464   0.951    0.373
    ## oxygen_median      -0.001076   0.007796  -0.138    0.894
    ## temperature_median  0.011174   0.006864   1.628    0.148
    ## 
    ## Residual standard error: 0.02941 on 7 degrees of freedom
    ## Multiple R-squared:  0.3004, Adjusted R-squared:  0.0005571 
    ## F-statistic: 1.002 on 3 and 7 DF,  p-value: 0.4464

``` r
p_values <- summary(PD_PDisp_env_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          1.0000000          1.0000000          0.5902869

``` r
# Ocean sites
PD_PRic_env_O_lm <- lm(log(pd.obs) ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_sites,])
summary(PD_PRic_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         -7.9322        NaN     NaN      NaN
    ## salinity_median      0.4257        NaN     NaN      NaN
    ## oxygen_median        0.3226        NaN     NaN      NaN
    ## temperature_median       NA         NA      NA       NA
    ## 
    ## Residual standard error: NaN on 0 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
p_values <- summary(PD_PRic_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
PD_z_env_O_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_sites,])
summary(PD_z_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -123.440        NaN     NaN      NaN
    ## salinity_median       4.168        NaN     NaN      NaN
    ## oxygen_median        -3.395        NaN     NaN      NaN
    ## temperature_median       NA         NA      NA       NA
    ## 
    ## Residual standard error: NaN on 0 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
p_values <- summary(PD_z_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
PD_PDisp_env_O_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_sites,])
summary(PD_PDisp_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.014774        NaN     NaN      NaN
    ## salinity_median     0.009109        NaN     NaN      NaN
    ## oxygen_median      -0.023091        NaN     NaN      NaN
    ## temperature_median        NA         NA      NA       NA
    ## 
    ## Residual standard error: NaN on 0 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
p_values <- summary(PD_PDisp_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
# Mixed lakes
PD_PRic_env_M_lm <- lm(log(pd.obs) ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[mixed_lakes,])
summary(PD_PRic_env_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##      FLK      HLO      LLN      MLN      NLN      NLU      OLO      ULN 
    ##  0.04358 -0.23393  0.33491 -0.15413  0.48318 -0.31577 -0.28003  0.12220 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.72550   18.55796  -0.039    0.971
    ## salinity_median     0.28019    0.43143   0.649    0.551
    ## oxygen_median       0.29681    0.34898   0.850    0.443
    ## temperature_median -0.06134    0.26779  -0.229    0.830
    ## 
    ## Residual standard error: 0.3934 on 4 degrees of freedom
    ## Multiple R-squared:  0.3807, Adjusted R-squared:  -0.08376 
    ## F-statistic: 0.8197 on 3 and 4 DF,  p-value: 0.5469

``` r
p_values <- summary(PD_PRic_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
PD_z_env_M_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[mixed_lakes,])
summary(PD_z_env_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##     FLK     HLO     LLN     MLN     NLN     NLU     OLO     ULN 
    ## -0.6311 -0.9440 -0.1008  0.4796  0.7484  0.1079  0.1433  0.1967 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         18.8693    34.6682   0.544    0.615
    ## salinity_median     -0.7909     0.8060  -0.981    0.382
    ## oxygen_median       -0.6458     0.6519  -0.991    0.378
    ## temperature_median   0.3083     0.5003   0.616    0.571
    ## 
    ## Residual standard error: 0.7349 on 4 degrees of freedom
    ## Multiple R-squared:  0.5607, Adjusted R-squared:  0.2313 
    ## F-statistic: 1.702 on 3 and 4 DF,  p-value: 0.3034

``` r
p_values <- summary(PD_z_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_M_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[mixed_lakes,])
summary(PD_PDisp_env_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[mixed_lakes, ])
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
    ## salinity_median     0.002243   0.012316   0.182    0.864
    ## oxygen_median      -0.004627   0.009963  -0.464    0.666
    ## temperature_median  0.004574   0.007645   0.598    0.582
    ## 
    ## Residual standard error: 0.01123 on 4 degrees of freedom
    ## Multiple R-squared:  0.1322, Adjusted R-squared:  -0.5187 
    ## F-statistic: 0.2031 on 3 and 4 DF,  p-value: 0.8894

``` r
p_values <- summary(PD_PDisp_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
# Stratified lakes
PD_PRic_env_S_lm <- lm(log(pd.obs) ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[stratified_lakes,])
summary(PD_PRic_env_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = stree_sespd_env[stratified_lakes, 
    ##     ])
    ## 
    ## Residuals:
    ##       BCM       CLM       GLK       HLM       NLK       OTM       SLN       TLN 
    ##  0.222256  0.002541  0.093548  0.144056 -0.223573 -0.278494 -0.109216  0.148881 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         4.17718    2.71151   1.541    0.198
    ## salinity_median     0.06782    0.04211   1.610    0.183
    ## oxygen_median      -0.03224    0.18566  -0.174    0.871
    ## temperature_median  0.02913    0.06006   0.485    0.653
    ## 
    ## Residual standard error: 0.2452 on 4 degrees of freedom
    ## Multiple R-squared:  0.6503, Adjusted R-squared:  0.388 
    ## F-statistic: 2.479 on 3 and 4 DF,  p-value: 0.2006

``` r
p_values <- summary(PD_PRic_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          0.7931036          0.7303422          1.0000000          1.0000000

``` r
PD_z_env_S_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[stratified_lakes,])
summary(PD_z_env_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##     BCM     CLM     GLK     HLM     NLK     OTM     SLN     TLN 
    ##  0.6329 -0.5528  0.2931  0.0925 -0.7839  1.2305 -0.4533 -0.4590 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -9.36668   10.11196  -0.926    0.407
    ## salinity_median     0.03683    0.15704   0.235    0.826
    ## oxygen_median       0.57919    0.69238   0.837    0.450
    ## temperature_median  0.18240    0.22397   0.814    0.461
    ## 
    ## Residual standard error: 0.9145 on 4 degrees of freedom
    ## Multiple R-squared:  0.3152, Adjusted R-squared:  -0.1983 
    ## F-statistic: 0.6138 on 3 and 4 DF,  p-value: 0.6412

``` r
p_values <- summary(PD_z_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
PD_PDisp_env_S_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[stratified_lakes,])
summary(PD_PDisp_env_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##       BCM       CLM       GLK       HLM       NLK       OTM       SLN       TLN 
    ##  0.018547 -0.026500  0.037634 -0.004772 -0.034031  0.026157 -0.020451  0.003415 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -0.491926   0.381330  -1.290    0.267
    ## salinity_median     0.006309   0.005922   1.065    0.347
    ## oxygen_median       0.018246   0.026110   0.699    0.523
    ## temperature_median  0.013502   0.008446   1.599    0.185
    ## 
    ## Residual standard error: 0.03449 on 4 degrees of freedom
    ## Multiple R-squared:  0.4408, Adjusted R-squared:  0.0214 
    ## F-statistic: 1.051 on 3 and 4 DF,  p-value: 0.4619

``` r
p_values <- summary(PD_PDisp_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          1.0000000          1.0000000          0.7405587

``` r
### Geographical
# Surveyed sites
PD_PRic_geo_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PRic_geo_am)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1  4.802   4.802  56.041  2.1e-05 ***
    ## max_depth                          1  0.000   0.000   0.000 0.997311    
    ## logArea                            1  0.078   0.078   0.909 0.362755    
    ## Site_type                          2  2.571   1.286  15.005 0.000975 ***
    ## distance_to_ocean_min_m:Site_type  2  0.034   0.017   0.197 0.824545    
    ## max_depth:Site_type                2  0.138   0.069   0.808 0.472782    
    ## logArea:Site_type                  2  0.346   0.173   2.022 0.183079    
    ## Residuals                         10  0.857   0.086                     
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
PD_z_geo_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_z_geo_am)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1  4.802   4.802  56.041  2.1e-05 ***
    ## max_depth                          1  0.000   0.000   0.000 0.997311    
    ## logArea                            1  0.078   0.078   0.909 0.362755    
    ## Site_type                          2  2.571   1.286  15.005 0.000975 ***
    ## distance_to_ocean_min_m:Site_type  2  0.034   0.017   0.197 0.824545    
    ## max_depth:Site_type                2  0.138   0.069   0.808 0.472782    
    ## logArea:Site_type                  2  0.346   0.173   2.022 0.183079    
    ## Residuals                         10  0.857   0.086                     
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
PD_PDisp_geo_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PDisp_geo_am)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1  4.802   4.802  56.041  2.1e-05 ***
    ## max_depth                          1  0.000   0.000   0.000 0.997311    
    ## logArea                            1  0.078   0.078   0.909 0.362755    
    ## Site_type                          2  2.571   1.286  15.005 0.000975 ***
    ## distance_to_ocean_min_m:Site_type  2  0.034   0.017   0.197 0.824545    
    ## max_depth:Site_type                2  0.138   0.069   0.808 0.472782    
    ## logArea:Site_type                  2  0.346   0.173   2.022 0.183079    
    ## Residuals                         10  0.857   0.086                     
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
PD_PRic_geo_MS_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_geo_MS_am)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1  3.604   3.604  48.493 0.000117 ***
    ## max_depth                          1  0.023   0.023   0.308 0.593953    
    ## logArea                            1  1.285   1.285  17.295 0.003171 ** 
    ## Site_type                          1  1.610   1.610  21.671 0.001634 ** 
    ## distance_to_ocean_min_m:Site_type  1  0.003   0.003   0.039 0.848475    
    ## max_depth:Site_type                1  0.048   0.048   0.648 0.444062    
    ## logArea:Site_type                  1  0.004   0.004   0.050 0.829172    
    ## Residuals                          8  0.594   0.074                     
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
PD_z_geo_MS_am <- aov(pd.obs.z ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_z_geo_MS_am)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m            1  0.201  0.2007   0.313 0.5909  
    ## max_depth                          1  0.040  0.0401   0.063 0.8086  
    ## logArea                            1  0.475  0.4747   0.741 0.4143  
    ## Site_type                          1  0.014  0.0144   0.022 0.8847  
    ## distance_to_ocean_min_m:Site_type  1  0.054  0.0545   0.085 0.7780  
    ## max_depth:Site_type                1  2.287  2.2872   3.572 0.0954 .
    ## logArea:Site_type                  1  1.610  1.6104   2.515 0.1514  
    ## Residuals                          8  5.122  0.6402                 
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
PD_PDisp_geo_MS_am <- aov(Dispersion ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_geo_MS_am)
```

    ##                                   Df    Sum Sq   Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m            1 0.0019693 0.0019693   6.241 0.0370 *
    ## max_depth                          1 0.0006945 0.0006945   2.201 0.1762  
    ## logArea                            1 0.0009206 0.0009206   2.917 0.1260  
    ## Site_type                          1 0.0000004 0.0000004   0.001 0.9731  
    ## distance_to_ocean_min_m:Site_type  1 0.0000779 0.0000779   0.247 0.6325  
    ## max_depth:Site_type                1 0.0013180 0.0013180   4.177 0.0752 .
    ## logArea:Site_type                  1 0.0022907 0.0022907   7.260 0.0273 *
    ## Residuals                          8 0.0025243 0.0003155                 
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
PD_PRic_geo_OM_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PRic_geo_OM_am)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m            1 0.0038  0.0038   0.045 0.8385  
    ## max_depth                          1 0.3603  0.3603   4.268 0.0844 .
    ## logArea                            1 0.0498  0.0498   0.590 0.4715  
    ## Site_type                          1 0.0243  0.0243   0.288 0.6107  
    ## distance_to_ocean_min_m:Site_type  1 0.0007  0.0007   0.008 0.9300  
    ## max_depth:Site_type                1 0.0842  0.0842   0.997 0.3566  
    ## logArea:Site_type                  1 0.3205  0.3205   3.796 0.0993 .
    ## Residuals                          6 0.5065  0.0844                 
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
PD_z_geo_OM_am <- aov(pd.obs.z ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_z_geo_OM_am)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m            1  1.387   1.387   1.871 0.2204  
    ## max_depth                          1  8.062   8.062  10.874 0.0165 *
    ## logArea                            1  0.094   0.094   0.126 0.7346  
    ## Site_type                          1  1.077   1.077   1.453 0.2734  
    ## distance_to_ocean_min_m:Site_type  1  0.077   0.077   0.104 0.7583  
    ## max_depth:Site_type                1  0.858   0.858   1.157 0.3234  
    ## logArea:Site_type                  1  0.223   0.223   0.301 0.6032  
    ## Residuals                          6  4.448   0.741                 
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
PD_PDisp_geo_OM_am <- aov(Dispersion ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PDisp_geo_OM_am)
```

    ##                                   Df    Sum Sq   Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1 0.0000369 3.687e-05   0.434  0.534
    ## max_depth                          1 0.0001463 1.463e-04   1.723  0.237
    ## logArea                            1 0.0000945 9.450e-05   1.113  0.332
    ## Site_type                          1 0.0002191 2.191e-04   2.579  0.159
    ## distance_to_ocean_min_m:Site_type  1 0.0000172 1.719e-05   0.202  0.669
    ## max_depth:Site_type                1 0.0000351 3.512e-05   0.413  0.544
    ## logArea:Site_type                  1 0.0001612 1.612e-04   1.898  0.217
    ## Residuals                          6 0.0005096 8.494e-05

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
PD_PRic_geo_SO_am <- aov(log(pd.obs) ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PRic_geo_SO_am)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1  4.012   4.012  39.287 0.000766 ***
    ## max_depth                          1  0.043   0.043   0.419 0.541244    
    ## logArea                            1  0.077   0.077   0.753 0.418950    
    ## Site_type                          1  0.624   0.624   6.109 0.048356 *  
    ## distance_to_ocean_min_m:Site_type  1  0.001   0.001   0.011 0.919385    
    ## max_depth:Site_type                1  0.008   0.008   0.081 0.785905    
    ## logArea:Site_type                  1  0.060   0.060   0.586 0.473152    
    ## Residuals                          6  0.613   0.102                     
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
PD_z_geo_SO_am <- aov(pd.obs.z ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_z_geo_SO_am)
```

    ##                                   Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m            1  1.287   1.287   2.746 0.14855   
    ## max_depth                          1  1.433   1.433   3.059 0.13088   
    ## logArea                            1  0.011   0.011   0.023 0.88351   
    ## Site_type                          1  0.520   0.520   1.111 0.33254   
    ## distance_to_ocean_min_m:Site_type  1  0.582   0.582   1.242 0.30765   
    ## max_depth:Site_type                1  6.959   6.959  14.851 0.00842 **
    ## logArea:Site_type                  1  2.590   2.590   5.527 0.05698 . 
    ## Residuals                          6  2.812   0.469                   
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
PD_PDisp_geo_SO_am <- aov(Dispersion ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PDisp_geo_SO_am)
```

    ##                                   Df   Sum Sq  Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m            1 0.000885 0.000885   2.350 0.1761  
    ## max_depth                          1 0.000625 0.000625   1.660 0.2451  
    ## logArea                            1 0.000006 0.000006   0.016 0.9033  
    ## Site_type                          1 0.000385 0.000385   1.022 0.3511  
    ## distance_to_ocean_min_m:Site_type  1 0.000011 0.000011   0.031 0.8671  
    ## max_depth:Site_type                1 0.001583 0.001583   4.205 0.0862 .
    ## logArea:Site_type                  1 0.003353 0.003353   8.907 0.0245 *
    ## Residuals                          6 0.002259 0.000376                 
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
PD_alpha_LA_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.z, x = logArea, color = Site_type, fill = Site_type)) + 
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
  labs(x="Log Area (m"^"2"~")", y="PD alpha z-score", colour = "Site type:", fill = "Site type:", tag = "c")
(PD_alpha_LA_plot <- PD_alpha_LA_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_alpha_LA_plot.jpg", plot = PD_alpha_LA_plot, width = 6.26, height = 6, units = "in")

# Distance from the ocean mean
PD_alpha_D_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.z, x = distance_to_ocean_min_m, color = Site_type, fill = Site_type)) + 
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
  labs(x="Isolation (m)", y="PD alpha z-score", colour = "Site type:", fill = "Site type:", tag = "b")
(PD_alpha_D_plot <- PD_alpha_D_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_alpha_D_plot.jpg", plot = PD_alpha_D_plot, width = 6.26, height = 6, units = "in")

# Max depth
PD_alpha_MD_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.z, x = max_depth, color = Site_type, fill = Site_type)) + 
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
  labs(x="Age (m)", y="PD alpha z-score", colour = "Site type:", fill = "Site type:", tag = "a")
(PD_alpha_MD_plot <- PD_alpha_MD_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/PD_alpha_MD_plot.jpg", plot = PD_alpha_MD_plot, width = 6.26, height = 6, units = "in")
```

### PD alpha z-scores & dispersion

``` r
PD_alpha_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.z, x = Dispersion, color = Site_type, fill = Site_type)) +
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
  labs(x="PD alpha Dispersion", y="PD alpha z-score", colour = "Site type:", fill = "Site type:", tag = "a")
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
(PD_beta_ref_BD <- betadisper(PD_beta_ref_dist$Btotal, Site_type_group_ref))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = PD_beta_ref_dist$Btotal, group =
    ## Site_type_group_ref)
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
    ## 1     Stratified vs Mixed  1 1.2940870 9.415922 0.4021162   0.001      0.006
    ## 2     Stratified vs Ocean  1 1.1858559 8.357071 0.4105242   0.002      0.012
    ## 3 Stratified vs Reference  1 0.7729479 7.110647 0.5039207   0.100      0.600
    ## 4          Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.101      0.606
    ## 5      Mixed vs Reference  1 0.6592270 3.967203 0.3617334   0.125      0.750
    ## 6      Ocean vs Reference  1 0.6319662 3.354877 0.4015472   0.146      0.876
    ##   sig
    ## 1   *
    ## 2   .
    ## 3    
    ## 4    
    ## 5    
    ## 6

``` r
# Surveyed sites
PD_beta_dist <- BAT::beta(surveyed_sites_lake, stree, abund = F)
(PD_beta_BD <- betadisper(PD_beta_dist$Btotal, Site_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = PD_beta_dist$Btotal, group = Site_type_group)
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
    ## 3      Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.090      0.270

``` r
# Without LCN
PD_beta_wo_LCN_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_LCN,], stree, abund = F)
(PD_beta_wo_LCN_BD <- betadisper(PD_beta_wo_LCN_dist$Btotal, Site_type_group[-8]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = PD_beta_wo_LCN_dist$Btotal, group =
    ## Site_type_group[-8])
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
    ## 2 Stratified vs Ocean  1 1.2083230 9.439158 0.4618174   0.004      0.012   .
    ## 3      Mixed vs Ocean  1 0.3347638 2.034034 0.1560556   0.014      0.042   .

``` r
# Without TLN and HLM
PD_beta_wo_TLN_HLM_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_TLN_HLM,], stree, abund = F)
(PD_beta_wo_TLN_HLM_BD <- betadisper(PD_beta_wo_TLN_HLM_dist$Btotal, Site_type_group[-c(5,21)]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = PD_beta_wo_TLN_HLM_dist$Btotal, group =
    ## Site_type_group[-c(5, 21)])
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
    ## 2 Stratified vs Ocean  1 1.2146330  9.192717 0.4789690   0.002      0.006   *
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

PD_beta_BD_dist_env$Site_type <- factor(PD_beta_BD_dist_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

(PD_beta_BD_dist_env_plot <- ggplot(PD_beta_BD_dist_env, aes(x = Site_type, y = PD_beta_BD_dist, color = Site_type, fill = Site_type)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Site_type, fill = Site_type)) +
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
  labs(color = "Site_type", tag = "b"))
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
    ## Run 1 stress 6.798418e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002386498  max resid 0.00039228 
    ## ... Similar to previous best
    ## Run 2 stress 9.982179e-05 
    ## ... Procrustes: rmse 5.675929e-05  max resid 0.0001309431 
    ## ... Similar to previous best
    ## Run 3 stress 9.244747e-05 
    ## ... Procrustes: rmse 0.0001962285  max resid 0.0003252607 
    ## ... Similar to previous best
    ## Run 4 stress 9.989435e-05 
    ## ... Procrustes: rmse 0.00022641  max resid 0.000376785 
    ## ... Similar to previous best
    ## Run 5 stress 9.996774e-05 
    ## ... Procrustes: rmse 0.0002564658  max resid 0.0004297573 
    ## ... Similar to previous best
    ## Run 6 stress 9.97744e-05 
    ## ... Procrustes: rmse 0.0002607467  max resid 0.0004348562 
    ## ... Similar to previous best
    ## Run 7 stress 7.250394e-05 
    ## ... Procrustes: rmse 4.685126e-05  max resid 9.638365e-05 
    ## ... Similar to previous best
    ## Run 8 stress 9.928659e-05 
    ## ... Procrustes: rmse 5.31893e-05  max resid 8.194523e-05 
    ## ... Similar to previous best
    ## Run 9 stress 9.715532e-05 
    ## ... Procrustes: rmse 0.0002168378  max resid 0.0003662847 
    ## ... Similar to previous best
    ## Run 10 stress 8.893505e-05 
    ## ... Procrustes: rmse 0.0001816228  max resid 0.0003094947 
    ## ... Similar to previous best
    ## Run 11 stress 9.804733e-05 
    ## ... Procrustes: rmse 0.0002175309  max resid 0.0003712279 
    ## ... Similar to previous best
    ## Run 12 stress 9.531072e-05 
    ## ... Procrustes: rmse 0.0002131916  max resid 0.0003551992 
    ## ... Similar to previous best
    ## Run 13 stress 9.971624e-05 
    ## ... Procrustes: rmse 0.0002474813  max resid 0.0004134496 
    ## ... Similar to previous best
    ## Run 14 stress 9.917956e-05 
    ## ... Procrustes: rmse 0.0002149208  max resid 0.0003642387 
    ## ... Similar to previous best
    ## Run 15 stress 9.425879e-05 
    ## ... Procrustes: rmse 0.0002059148  max resid 0.0003517668 
    ## ... Similar to previous best
    ## Run 16 stress 9.651611e-05 
    ## ... Procrustes: rmse 0.0002496925  max resid 0.0004237667 
    ## ... Similar to previous best
    ## Run 17 stress 4.459683e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 4.525342e-05  max resid 9.71169e-05 
    ## ... Similar to previous best
    ## Run 18 stress 9.93032e-05 
    ## ... Procrustes: rmse 0.0002264949  max resid 0.0003842432 
    ## ... Similar to previous best
    ## Run 19 stress 9.642329e-05 
    ## ... Procrustes: rmse 0.0002024447  max resid 0.0003369485 
    ## ... Similar to previous best
    ## Run 20 stress 9.994952e-05 
    ## ... Procrustes: rmse 0.0002585253  max resid 0.0004278014 
    ## ... Similar to previous best
    ## Run 21 stress 8.04138e-05 
    ## ... Procrustes: rmse 4.86557e-05  max resid 0.0001212015 
    ## ... Similar to previous best
    ## Run 22 stress 9.854248e-05 
    ## ... Procrustes: rmse 0.0002203985  max resid 0.0003688758 
    ## ... Similar to previous best
    ## Run 23 stress 9.538415e-05 
    ## ... Procrustes: rmse 0.000188791  max resid 0.0003337818 
    ## ... Similar to previous best
    ## Run 24 stress 4.343281e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 3.865866e-05  max resid 8.348451e-05 
    ## ... Similar to previous best
    ## Run 25 stress 9.925986e-05 
    ## ... Procrustes: rmse 0.0002267084  max resid 0.0003572137 
    ## ... Similar to previous best
    ## Run 26 stress 9.883937e-05 
    ## ... Procrustes: rmse 0.0002263359  max resid 0.0003634272 
    ## ... Similar to previous best
    ## Run 27 stress 8.054055e-05 
    ## ... Procrustes: rmse 8.208593e-05  max resid 0.0001367194 
    ## ... Similar to previous best
    ## Run 28 stress 4.052575e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 4.064344e-05  max resid 7.035911e-05 
    ## ... Similar to previous best
    ## Run 29 stress 9.96302e-05 
    ## ... Procrustes: rmse 0.000227878  max resid 0.0003661531 
    ## ... Similar to previous best
    ## Run 30 stress 6.769681e-05 
    ## ... Procrustes: rmse 5.045491e-05  max resid 8.199511e-05 
    ## ... Similar to previous best
    ## Run 31 stress 8.813701e-05 
    ## ... Procrustes: rmse 0.0001469865  max resid 0.0002304239 
    ## ... Similar to previous best
    ## Run 32 stress 9.663168e-05 
    ## ... Procrustes: rmse 0.0001297956  max resid 0.0002255483 
    ## ... Similar to previous best
    ## Run 33 stress 9.84988e-05 
    ## ... Procrustes: rmse 0.0002224396  max resid 0.0003408144 
    ## ... Similar to previous best
    ## Run 34 stress 9.827808e-05 
    ## ... Procrustes: rmse 0.000255684  max resid 0.000400775 
    ## ... Similar to previous best
    ## Run 35 stress 6.916033e-05 
    ## ... Procrustes: rmse 3.756156e-05  max resid 7.18954e-05 
    ## ... Similar to previous best
    ## Run 36 stress 9.758903e-05 
    ## ... Procrustes: rmse 0.0002183938  max resid 0.000339049 
    ## ... Similar to previous best
    ## Run 37 stress 9.845996e-05 
    ## ... Procrustes: rmse 0.0002244066  max resid 0.0003490212 
    ## ... Similar to previous best
    ## Run 38 stress 9.76999e-05 
    ## ... Procrustes: rmse 0.0002202412  max resid 0.0003654395 
    ## ... Similar to previous best
    ## Run 39 stress 8.992296e-05 
    ## ... Procrustes: rmse 5.041619e-05  max resid 9.53901e-05 
    ## ... Similar to previous best
    ## Run 40 stress 9.644139e-05 
    ## ... Procrustes: rmse 0.0002531202  max resid 0.0003986462 
    ## ... Similar to previous best
    ## Run 41 stress 9.914556e-05 
    ## ... Procrustes: rmse 0.0002284597  max resid 0.0003677264 
    ## ... Similar to previous best
    ## Run 42 stress 7.667185e-05 
    ## ... Procrustes: rmse 4.082971e-05  max resid 7.881298e-05 
    ## ... Similar to previous best
    ## Run 43 stress 9.70236e-05 
    ## ... Procrustes: rmse 0.0002536469  max resid 0.0004022015 
    ## ... Similar to previous best
    ## Run 44 stress 9.83963e-05 
    ## ... Procrustes: rmse 0.0002452999  max resid 0.0004055593 
    ## ... Similar to previous best
    ## Run 45 stress 9.939428e-05 
    ## ... Procrustes: rmse 0.0002282383  max resid 0.0003672964 
    ## ... Similar to previous best
    ## Run 46 stress 9.804567e-05 
    ## ... Procrustes: rmse 0.0002213839  max resid 0.0003531274 
    ## ... Similar to previous best
    ## Run 47 stress 6.847729e-05 
    ## ... Procrustes: rmse 4.197944e-05  max resid 7.717474e-05 
    ## ... Similar to previous best
    ## Run 48 stress 9.309917e-05 
    ## ... Procrustes: rmse 0.0002138851  max resid 0.0003450343 
    ## ... Similar to previous best
    ## Run 49 stress 9.343456e-05 
    ## ... Procrustes: rmse 5.008606e-05  max resid 0.0001069461 
    ## ... Similar to previous best
    ## Run 50 stress 0.3793399 
    ## Run 51 stress 9.866183e-05 
    ## ... Procrustes: rmse 0.0002172871  max resid 0.0003471107 
    ## ... Similar to previous best
    ## Run 52 stress 9.993799e-05 
    ## ... Procrustes: rmse 0.0002535174  max resid 0.0004143409 
    ## ... Similar to previous best
    ## Run 53 stress 8.65951e-05 
    ## ... Procrustes: rmse 4.373669e-05  max resid 8.638337e-05 
    ## ... Similar to previous best
    ## Run 54 stress 9.966144e-05 
    ## ... Procrustes: rmse 0.0001932255  max resid 0.0003182842 
    ## ... Similar to previous best
    ## Run 55 stress 9.91576e-05 
    ## ... Procrustes: rmse 0.0002569276  max resid 0.0004029074 
    ## ... Similar to previous best
    ## Run 56 stress 9.556775e-05 
    ## ... Procrustes: rmse 0.0002387939  max resid 0.0003640581 
    ## ... Similar to previous best
    ## Run 57 stress 7.463744e-05 
    ## ... Procrustes: rmse 5.450221e-05  max resid 7.404054e-05 
    ## ... Similar to previous best
    ## Run 58 stress 9.762018e-05 
    ## ... Procrustes: rmse 0.0002262894  max resid 0.0003657797 
    ## ... Similar to previous best
    ## Run 59 stress 8.837095e-05 
    ## ... Procrustes: rmse 0.0001033075  max resid 0.0001859139 
    ## ... Similar to previous best
    ## Run 60 stress 9.430145e-05 
    ## ... Procrustes: rmse 0.0001519051  max resid 0.0002367105 
    ## ... Similar to previous best
    ## Run 61 stress 9.794921e-05 
    ## ... Procrustes: rmse 0.0002220704  max resid 0.00034778 
    ## ... Similar to previous best
    ## Run 62 stress 9.941611e-05 
    ## ... Procrustes: rmse 0.0002572812  max resid 0.0004021169 
    ## ... Similar to previous best
    ## Run 63 stress 9.756652e-05 
    ## ... Procrustes: rmse 0.0002214382  max resid 0.0003720336 
    ## ... Similar to previous best
    ## Run 64 stress 9.62703e-05 
    ## ... Procrustes: rmse 0.0002106785  max resid 0.0003661147 
    ## ... Similar to previous best
    ## Run 65 stress 9.970925e-05 
    ## ... Procrustes: rmse 0.0002274654  max resid 0.0003587848 
    ## ... Similar to previous best
    ## Run 66 stress 9.606922e-05 
    ## ... Procrustes: rmse 0.0002519828  max resid 0.0004022942 
    ## ... Similar to previous best
    ## Run 67 stress 9.954194e-05 
    ## ... Procrustes: rmse 0.0002232783  max resid 0.0003673914 
    ## ... Similar to previous best
    ## Run 68 stress 9.856903e-05 
    ## ... Procrustes: rmse 0.0002556197  max resid 0.0004113106 
    ## ... Similar to previous best
    ## Run 69 stress 9.956144e-05 
    ## ... Procrustes: rmse 0.0001818354  max resid 0.0003071425 
    ## ... Similar to previous best
    ## Run 70 stress 9.952174e-05 
    ## ... Procrustes: rmse 0.0002542805  max resid 0.0004073152 
    ## ... Similar to previous best
    ## Run 71 stress 9.587563e-05 
    ## ... Procrustes: rmse 0.0002160377  max resid 0.0003699178 
    ## ... Similar to previous best
    ## Run 72 stress 9.948061e-05 
    ## ... Procrustes: rmse 0.0002540558  max resid 0.000397779 
    ## ... Similar to previous best
    ## Run 73 stress 9.952532e-05 
    ## ... Procrustes: rmse 0.0002533702  max resid 0.0003961973 
    ## ... Similar to previous best
    ## Run 74 stress 7.143888e-05 
    ## ... Procrustes: rmse 4.220668e-05  max resid 7.154232e-05 
    ## ... Similar to previous best
    ## Run 75 stress 9.581487e-05 
    ## ... Procrustes: rmse 0.0002044526  max resid 0.0003606572 
    ## ... Similar to previous best
    ## Run 76 stress 9.864963e-05 
    ## ... Procrustes: rmse 0.000226531  max resid 0.0003561057 
    ## ... Similar to previous best
    ## Run 77 stress 9.660367e-05 
    ## ... Procrustes: rmse 0.0002491148  max resid 0.0003990517 
    ## ... Similar to previous best
    ## Run 78 stress 9.068687e-05 
    ## ... Procrustes: rmse 0.0001916911  max resid 0.000341111 
    ## ... Similar to previous best
    ## Run 79 stress 9.939864e-05 
    ## ... Procrustes: rmse 0.000112449  max resid 0.0001900292 
    ## ... Similar to previous best
    ## Run 80 stress 9.974945e-05 
    ## ... Procrustes: rmse 5.011235e-05  max resid 9.739037e-05 
    ## ... Similar to previous best
    ## Run 81 stress 9.772355e-05 
    ## ... Procrustes: rmse 0.000219515  max resid 0.0003453759 
    ## ... Similar to previous best
    ## Run 82 stress 9.74025e-05 
    ## ... Procrustes: rmse 0.0001795086  max resid 0.0003230335 
    ## ... Similar to previous best
    ## Run 83 stress 6.858265e-05 
    ## ... Procrustes: rmse 3.998354e-05  max resid 6.429005e-05 
    ## ... Similar to previous best
    ## Run 84 stress 9.992483e-05 
    ## ... Procrustes: rmse 0.0002317071  max resid 0.0003670092 
    ## ... Similar to previous best
    ## Run 85 stress 5.677672e-05 
    ## ... Procrustes: rmse 7.295398e-05  max resid 0.0001253404 
    ## ... Similar to previous best
    ## Run 86 stress 9.585423e-05 
    ## ... Procrustes: rmse 0.0002060795  max resid 0.0003553811 
    ## ... Similar to previous best
    ## Run 87 stress 8.375758e-05 
    ## ... Procrustes: rmse 5.228943e-05  max resid 9.114085e-05 
    ## ... Similar to previous best
    ## Run 88 stress 7.517466e-05 
    ## ... Procrustes: rmse 4.514627e-05  max resid 8.427043e-05 
    ## ... Similar to previous best
    ## Run 89 stress 9.961893e-05 
    ## ... Procrustes: rmse 0.0002468234  max resid 0.0003972754 
    ## ... Similar to previous best
    ## Run 90 stress 9.508026e-05 
    ## ... Procrustes: rmse 0.0001466115  max resid 0.0002773064 
    ## ... Similar to previous best
    ## Run 91 stress 9.507736e-05 
    ## ... Procrustes: rmse 0.0002211928  max resid 0.0003551582 
    ## ... Similar to previous best
    ## Run 92 stress 8.762118e-05 
    ## ... Procrustes: rmse 4.821902e-05  max resid 8.945719e-05 
    ## ... Similar to previous best
    ## Run 93 stress 8.774173e-05 
    ## ... Procrustes: rmse 0.000141344  max resid 0.0002415649 
    ## ... Similar to previous best
    ## Run 94 stress 9.713172e-05 
    ## ... Procrustes: rmse 0.0002176263  max resid 0.0003540665 
    ## ... Similar to previous best
    ## Run 95 stress 9.285464e-05 
    ## ... Procrustes: rmse 5.997755e-05  max resid 9.639376e-05 
    ## ... Similar to previous best
    ## Run 96 stress 9.334483e-05 
    ## ... Procrustes: rmse 0.0002135344  max resid 0.0003472416 
    ## ... Similar to previous best
    ## Run 97 stress 9.928819e-05 
    ## ... Procrustes: rmse 0.0002496917  max resid 0.000410431 
    ## ... Similar to previous best
    ## Run 98 stress 9.918069e-05 
    ## ... Procrustes: rmse 0.0002549792  max resid 0.0003998317 
    ## ... Similar to previous best
    ## Run 99 stress 9.73511e-05 
    ## ... Procrustes: rmse 0.0002543128  max resid 0.000404537 
    ## ... Similar to previous best
    ## Run 100 stress 9.889004e-05 
    ## ... Procrustes: rmse 0.0002225172  max resid 0.0003679905 
    ## ... Similar to previous best
    ## Run 101 stress 7.796319e-05 
    ## ... Procrustes: rmse 0.000126933  max resid 0.0002314761 
    ## ... Similar to previous best
    ## Run 102 stress 9.837043e-05 
    ## ... Procrustes: rmse 9.087459e-05  max resid 0.0002023273 
    ## ... Similar to previous best
    ## Run 103 stress 9.845763e-05 
    ## ... Procrustes: rmse 0.0001956427  max resid 0.0003413083 
    ## ... Similar to previous best
    ## Run 104 stress 8.99801e-05 
    ## ... Procrustes: rmse 0.0001856225  max resid 0.0003360935 
    ## ... Similar to previous best
    ## Run 105 stress 9.99631e-05 
    ## ... Procrustes: rmse 5.159205e-05  max resid 8.809554e-05 
    ## ... Similar to previous best
    ## Run 106 stress 5.717635e-05 
    ## ... Procrustes: rmse 2.995492e-05  max resid 8.016742e-05 
    ## ... Similar to previous best
    ## Run 107 stress 6.147395e-05 
    ## ... Procrustes: rmse 3.475825e-05  max resid 6.678524e-05 
    ## ... Similar to previous best
    ## Run 108 stress 4.99772e-05 
    ## ... Procrustes: rmse 3.269133e-05  max resid 5.634277e-05 
    ## ... Similar to previous best
    ## Run 109 stress 9.698439e-05 
    ## ... Procrustes: rmse 0.0002519236  max resid 0.0003939405 
    ## ... Similar to previous best
    ## Run 110 stress 9.885377e-05 
    ## ... Procrustes: rmse 0.0001837462  max resid 0.0003287996 
    ## ... Similar to previous best
    ## Run 111 stress 5.990195e-05 
    ## ... Procrustes: rmse 3.370195e-05  max resid 5.983174e-05 
    ## ... Similar to previous best
    ## Run 112 stress 7.361923e-05 
    ## ... Procrustes: rmse 4.006615e-05  max resid 6.313491e-05 
    ## ... Similar to previous best
    ## Run 113 stress 9.094389e-05 
    ## ... Procrustes: rmse 8.128456e-05  max resid 0.0001836261 
    ## ... Similar to previous best
    ## Run 114 stress 9.985546e-05 
    ## ... Procrustes: rmse 0.0002290592  max resid 0.0003664498 
    ## ... Similar to previous best
    ## Run 115 stress 9.952914e-05 
    ## ... Procrustes: rmse 0.0002124405  max resid 0.0003835448 
    ## ... Similar to previous best
    ## Run 116 stress 9.061727e-05 
    ## ... Procrustes: rmse 5.972709e-05  max resid 0.0001047721 
    ## ... Similar to previous best
    ## Run 117 stress 9.734075e-05 
    ## ... Procrustes: rmse 0.0002026319  max resid 0.0003613003 
    ## ... Similar to previous best
    ## Run 118 stress 9.869102e-05 
    ## ... Procrustes: rmse 0.0002157808  max resid 0.0003726606 
    ## ... Similar to previous best
    ## Run 119 stress 9.964492e-05 
    ## ... Procrustes: rmse 0.0002489122  max resid 0.0003934417 
    ## ... Similar to previous best
    ## Run 120 stress 9.816379e-05 
    ## ... Procrustes: rmse 0.0002110859  max resid 0.0003736861 
    ## ... Similar to previous best
    ## Run 121 stress 9.987422e-05 
    ## ... Procrustes: rmse 0.00025597  max resid 0.000397093 
    ## ... Similar to previous best
    ## Run 122 stress 9.846258e-05 
    ## ... Procrustes: rmse 0.0002214564  max resid 0.0003595226 
    ## ... Similar to previous best
    ## Run 123 stress 9.485782e-05 
    ## ... Procrustes: rmse 0.0001982054  max resid 0.0003531042 
    ## ... Similar to previous best
    ## Run 124 stress 9.788346e-05 
    ## ... Procrustes: rmse 0.0002570383  max resid 0.0004070072 
    ## ... Similar to previous best
    ## Run 125 stress 7.177413e-05 
    ## ... Procrustes: rmse 4.182682e-05  max resid 8.339571e-05 
    ## ... Similar to previous best
    ## Run 126 stress 9.780717e-05 
    ## ... Procrustes: rmse 0.000220175  max resid 0.0003505006 
    ## ... Similar to previous best
    ## Run 127 stress 5.83656e-05 
    ## ... Procrustes: rmse 4.160733e-05  max resid 7.495805e-05 
    ## ... Similar to previous best
    ## Run 128 stress 9.775707e-05 
    ## ... Procrustes: rmse 0.0002190113  max resid 0.0003683265 
    ## ... Similar to previous best
    ## Run 129 stress 9.566809e-05 
    ## ... Procrustes: rmse 0.0002001816  max resid 0.0003574777 
    ## ... Similar to previous best
    ## Run 130 stress 9.852287e-05 
    ## ... Procrustes: rmse 0.0002274007  max resid 0.0003746938 
    ## ... Similar to previous best
    ## Run 131 stress 9.89664e-05 
    ## ... Procrustes: rmse 0.0002559836  max resid 0.0003994339 
    ## ... Similar to previous best
    ## Run 132 stress 9.981936e-05 
    ## ... Procrustes: rmse 0.0001765631  max resid 0.0002970144 
    ## ... Similar to previous best
    ## Run 133 stress 9.551349e-05 
    ## ... Procrustes: rmse 0.0002500763  max resid 0.0003968634 
    ## ... Similar to previous best
    ## Run 134 stress 9.900643e-05 
    ## ... Procrustes: rmse 0.0002579547  max resid 0.0004044245 
    ## ... Similar to previous best
    ## Run 135 stress 9.707792e-05 
    ## ... Procrustes: rmse 0.0002137161  max resid 0.0003733581 
    ## ... Similar to previous best
    ## Run 136 stress 9.872633e-05 
    ## ... Procrustes: rmse 0.0002551095  max resid 0.0003981224 
    ## ... Similar to previous best
    ## Run 137 stress 6.436084e-05 
    ## ... Procrustes: rmse 3.54824e-05  max resid 6.081529e-05 
    ## ... Similar to previous best
    ## Run 138 stress 9.981553e-05 
    ## ... Procrustes: rmse 0.0002291846  max resid 0.0003600043 
    ## ... Similar to previous best
    ## Run 139 stress 9.720837e-05 
    ## ... Procrustes: rmse 0.0002478208  max resid 0.0003859868 
    ## ... Similar to previous best
    ## Run 140 stress 4.442366e-05 
    ## ... Procrustes: rmse 3.0027e-05  max resid 6.002334e-05 
    ## ... Similar to previous best
    ## Run 141 stress 9.914614e-05 
    ## ... Procrustes: rmse 0.0002496708  max resid 0.0004022217 
    ## ... Similar to previous best
    ## Run 142 stress 9.782026e-05 
    ## ... Procrustes: rmse 0.000256085  max resid 0.000403925 
    ## ... Similar to previous best
    ## Run 143 stress 5.088251e-05 
    ## ... Procrustes: rmse 3.349415e-05  max resid 5.447541e-05 
    ## ... Similar to previous best
    ## Run 144 stress 9.945902e-05 
    ## ... Procrustes: rmse 0.0002471306  max resid 0.0003899717 
    ## ... Similar to previous best
    ## Run 145 stress 9.682037e-05 
    ## ... Procrustes: rmse 0.0002053554  max resid 0.0003221993 
    ## ... Similar to previous best
    ## Run 146 stress 9.77682e-05 
    ## ... Procrustes: rmse 0.0002520873  max resid 0.0004039426 
    ## ... Similar to previous best
    ## Run 147 stress 9.592284e-05 
    ## ... Procrustes: rmse 0.0001854707  max resid 0.0002915748 
    ## ... Similar to previous best
    ## Run 148 stress 8.853632e-05 
    ## ... Procrustes: rmse 5.113554e-05  max resid 7.919961e-05 
    ## ... Similar to previous best
    ## Run 149 stress 9.422076e-05 
    ## ... Procrustes: rmse 0.0001882717  max resid 0.0002955465 
    ## ... Similar to previous best
    ## Run 150 stress 9.220729e-05 
    ## ... Procrustes: rmse 5.227618e-05  max resid 9.608915e-05 
    ## ... Similar to previous best
    ## Run 151 stress 9.687328e-05 
    ## ... Procrustes: rmse 0.0002021666  max resid 0.0003572186 
    ## ... Similar to previous best
    ## Run 152 stress 9.722268e-05 
    ## ... Procrustes: rmse 0.0002195738  max resid 0.0003457042 
    ## ... Similar to previous best
    ## Run 153 stress 9.645019e-05 
    ## ... Procrustes: rmse 0.0001835663  max resid 0.0003107239 
    ## ... Similar to previous best
    ## Run 154 stress 9.686025e-05 
    ## ... Procrustes: rmse 0.0002449734  max resid 0.0004010733 
    ## ... Similar to previous best
    ## Run 155 stress 9.848971e-05 
    ## ... Procrustes: rmse 0.000219206  max resid 0.0003667768 
    ## ... Similar to previous best
    ## Run 156 stress 9.973739e-05 
    ## ... Procrustes: rmse 0.0002224873  max resid 0.0003446283 
    ## ... Similar to previous best
    ## Run 157 stress 9.861501e-05 
    ## ... Procrustes: rmse 0.0002463442  max resid 0.0003804175 
    ## ... Similar to previous best
    ## Run 158 stress 9.00638e-05 
    ## ... Procrustes: rmse 4.496673e-05  max resid 7.874757e-05 
    ## ... Similar to previous best
    ## Run 159 stress 8.685324e-05 
    ## ... Procrustes: rmse 4.284051e-05  max resid 9.02057e-05 
    ## ... Similar to previous best
    ## Run 160 stress 9.702586e-05 
    ## ... Procrustes: rmse 0.0002395813  max resid 0.0004130088 
    ## ... Similar to previous best
    ## Run 161 stress 9.996646e-05 
    ## ... Procrustes: rmse 0.0002502921  max resid 0.0004083193 
    ## ... Similar to previous best
    ## Run 162 stress 9.979541e-05 
    ## ... Procrustes: rmse 0.0002297882  max resid 0.0003686285 
    ## ... Similar to previous best
    ## Run 163 stress 7.159813e-05 
    ## ... Procrustes: rmse 0.0001097798  max resid 0.0001669209 
    ## ... Similar to previous best
    ## Run 164 stress 9.825851e-05 
    ## ... Procrustes: rmse 0.0002568115  max resid 0.0004078235 
    ## ... Similar to previous best
    ## Run 165 stress 9.986642e-05 
    ## ... Procrustes: rmse 0.0002581466  max resid 0.0004051546 
    ## ... Similar to previous best
    ## Run 166 stress 9.860064e-05 
    ## ... Procrustes: rmse 0.0002286343  max resid 0.0003615306 
    ## ... Similar to previous best
    ## Run 167 stress 4.804026e-05 
    ## ... Procrustes: rmse 2.543375e-05  max resid 7.956855e-05 
    ## ... Similar to previous best
    ## Run 168 stress 9.904998e-05 
    ## ... Procrustes: rmse 0.0002256332  max resid 0.0003522196 
    ## ... Similar to previous best
    ## Run 169 stress 9.971676e-05 
    ## ... Procrustes: rmse 0.0002137183  max resid 0.0003454529 
    ## ... Similar to previous best
    ## Run 170 stress 9.635689e-05 
    ## ... Procrustes: rmse 0.0002191682  max resid 0.000350302 
    ## ... Similar to previous best
    ## Run 171 stress 9.902474e-05 
    ## ... Procrustes: rmse 0.0002423813  max resid 0.0004158194 
    ## ... Similar to previous best
    ## Run 172 stress 9.752964e-05 
    ## ... Procrustes: rmse 0.0001995345  max resid 0.0003526993 
    ## ... Similar to previous best
    ## Run 173 stress 6.654264e-05 
    ## ... Procrustes: rmse 3.400667e-05  max resid 7.186758e-05 
    ## ... Similar to previous best
    ## Run 174 stress 7.093004e-05 
    ## ... Procrustes: rmse 4.430983e-05  max resid 7.16206e-05 
    ## ... Similar to previous best
    ## Run 175 stress 9.912079e-05 
    ## ... Procrustes: rmse 0.0002248075  max resid 0.0003655642 
    ## ... Similar to previous best
    ## Run 176 stress 9.816225e-05 
    ## ... Procrustes: rmse 0.0002522572  max resid 0.0004101951 
    ## ... Similar to previous best
    ## Run 177 stress 6.859958e-05 
    ## ... Procrustes: rmse 8.341892e-05  max resid 0.0001434534 
    ## ... Similar to previous best
    ## Run 178 stress 9.033217e-05 
    ## ... Procrustes: rmse 0.0001080048  max resid 0.0002013524 
    ## ... Similar to previous best
    ## Run 179 stress 8.82356e-05 
    ## ... Procrustes: rmse 6.215077e-05  max resid 0.0001015865 
    ## ... Similar to previous best
    ## Run 180 stress 9.223282e-05 
    ## ... Procrustes: rmse 0.0001825403  max resid 0.0002954839 
    ## ... Similar to previous best
    ## Run 181 stress 4.475769e-05 
    ## ... Procrustes: rmse 3.009543e-05  max resid 5.706678e-05 
    ## ... Similar to previous best
    ## Run 182 stress 6.609227e-05 
    ## ... Procrustes: rmse 3.685881e-05  max resid 7.3185e-05 
    ## ... Similar to previous best
    ## Run 183 stress 9.351649e-05 
    ## ... Procrustes: rmse 0.0002167753  max resid 0.0003464086 
    ## ... Similar to previous best
    ## Run 184 stress 9.99587e-05 
    ## ... Procrustes: rmse 0.0002602224  max resid 0.0004072401 
    ## ... Similar to previous best
    ## Run 185 stress 9.826425e-05 
    ## ... Procrustes: rmse 0.0002279894  max resid 0.0004010666 
    ## ... Similar to previous best
    ## Run 186 stress 9.688751e-05 
    ## ... Procrustes: rmse 7.737848e-05  max resid 0.0001489508 
    ## ... Similar to previous best
    ## Run 187 stress 9.221781e-05 
    ## ... Procrustes: rmse 0.0001269211  max resid 0.0002309185 
    ## ... Similar to previous best
    ## Run 188 stress 6.234865e-05 
    ## ... Procrustes: rmse 8.240872e-05  max resid 0.0001622458 
    ## ... Similar to previous best
    ## Run 189 stress 9.735325e-05 
    ## ... Procrustes: rmse 0.0002495424  max resid 0.0003891836 
    ## ... Similar to previous best
    ## Run 190 stress 9.775356e-05 
    ## ... Procrustes: rmse 0.0002494481  max resid 0.0003907775 
    ## ... Similar to previous best
    ## Run 191 stress 9.791183e-05 
    ## ... Procrustes: rmse 0.0002553194  max resid 0.0004123302 
    ## ... Similar to previous best
    ## Run 192 stress 9.898097e-05 
    ## ... Procrustes: rmse 0.0002527963  max resid 0.0003955256 
    ## ... Similar to previous best
    ## Run 193 stress 9.62657e-05 
    ## ... Procrustes: rmse 0.0001622699  max resid 0.0002785438 
    ## ... Similar to previous best
    ## Run 194 stress 9.759601e-05 
    ## ... Procrustes: rmse 0.0002526486  max resid 0.0003971982 
    ## ... Similar to previous best
    ## Run 195 stress 8.530112e-05 
    ## ... Procrustes: rmse 8.731405e-05  max resid 0.0001668058 
    ## ... Similar to previous best
    ## Run 196 stress 9.839733e-05 
    ## ... Procrustes: rmse 0.0002299378  max resid 0.0003722195 
    ## ... Similar to previous best
    ## Run 197 stress 9.864658e-05 
    ## ... Procrustes: rmse 0.0002541359  max resid 0.0003970152 
    ## ... Similar to previous best
    ## Run 198 stress 9.953049e-05 
    ## ... Procrustes: rmse 0.0001704882  max resid 0.0003041998 
    ## ... Similar to previous best
    ## Run 199 stress 9.973283e-05 
    ## ... Procrustes: rmse 0.0002303684  max resid 0.0003733735 
    ## ... Similar to previous best
    ## Run 200 stress 9.823164e-05 
    ## ... Procrustes: rmse 0.0002202263  max resid 0.0003577189 
    ## ... Similar to previous best
    ## Run 201 stress 5.971877e-05 
    ## ... Procrustes: rmse 6.829204e-05  max resid 0.0001222369 
    ## ... Similar to previous best
    ## Run 202 stress 6.545558e-05 
    ## ... Procrustes: rmse 3.558939e-05  max resid 6.459778e-05 
    ## ... Similar to previous best
    ## Run 203 stress 9.898542e-05 
    ## ... Procrustes: rmse 0.0002327613  max resid 0.0003822177 
    ## ... Similar to previous best
    ## Run 204 stress 9.999544e-05 
    ## ... Procrustes: rmse 0.0002625219  max resid 0.0004158072 
    ## ... Similar to previous best
    ## Run 205 stress 9.901263e-05 
    ## ... Procrustes: rmse 0.0002587836  max resid 0.0004066676 
    ## ... Similar to previous best
    ## Run 206 stress 9.830337e-05 
    ## ... Procrustes: rmse 0.0002080911  max resid 0.0003633076 
    ## ... Similar to previous best
    ## Run 207 stress 8.890041e-05 
    ## ... Procrustes: rmse 0.0001953771  max resid 0.000326481 
    ## ... Similar to previous best
    ## Run 208 stress 9.743327e-05 
    ## ... Procrustes: rmse 0.0002517526  max resid 0.0004026241 
    ## ... Similar to previous best
    ## Run 209 stress 9.909516e-05 
    ## ... Procrustes: rmse 0.0002592087  max resid 0.0004164147 
    ## ... Similar to previous best
    ## Run 210 stress 8.791255e-05 
    ## ... Procrustes: rmse 0.0001996796  max resid 0.000315772 
    ## ... Similar to previous best
    ## Run 211 stress 7.73824e-05 
    ## ... Procrustes: rmse 4.279703e-05  max resid 7.819701e-05 
    ## ... Similar to previous best
    ## Run 212 stress 9.824127e-05 
    ## ... Procrustes: rmse 0.0002581512  max resid 0.0004103849 
    ## ... Similar to previous best
    ## Run 213 stress 9.073368e-05 
    ## ... Procrustes: rmse 0.0001905626  max resid 0.0003427926 
    ## ... Similar to previous best
    ## Run 214 stress 9.857041e-05 
    ## ... Procrustes: rmse 0.0002534954  max resid 0.0003970891 
    ## ... Similar to previous best
    ## Run 215 stress 9.929418e-05 
    ## ... Procrustes: rmse 0.0001979219  max resid 0.0003284829 
    ## ... Similar to previous best
    ## Run 216 stress 8.50574e-05 
    ## ... Procrustes: rmse 3.99907e-05  max resid 7.726637e-05 
    ## ... Similar to previous best
    ## Run 217 stress 8.976532e-05 
    ## ... Procrustes: rmse 4.731241e-05  max resid 8.732964e-05 
    ## ... Similar to previous best
    ## Run 218 stress 9.846492e-05 
    ## ... Procrustes: rmse 0.0002180252  max resid 0.0003678321 
    ## ... Similar to previous best
    ## Run 219 stress 9.934882e-05 
    ## ... Procrustes: rmse 0.0002271554  max resid 0.0003661546 
    ## ... Similar to previous best
    ## Run 220 stress 9.495546e-05 
    ## ... Procrustes: rmse 0.0002143732  max resid 0.0003382822 
    ## ... Similar to previous best
    ## Run 221 stress 9.968801e-05 
    ## ... Procrustes: rmse 0.0002600302  max resid 0.0004079042 
    ## ... Similar to previous best
    ## Run 222 stress 9.506256e-05 
    ## ... Procrustes: rmse 0.0002184272  max resid 0.0003479299 
    ## ... Similar to previous best
    ## Run 223 stress 9.981809e-05 
    ## ... Procrustes: rmse 0.000259308  max resid 0.0004093478 
    ## ... Similar to previous best
    ## Run 224 stress 7.781492e-05 
    ## ... Procrustes: rmse 4.411904e-05  max resid 8.147811e-05 
    ## ... Similar to previous best
    ## Run 225 stress 6.410869e-05 
    ## ... Procrustes: rmse 4.312332e-05  max resid 6.429781e-05 
    ## ... Similar to previous best
    ## Run 226 stress 9.983539e-05 
    ## ... Procrustes: rmse 0.0002541287  max resid 0.0003978639 
    ## ... Similar to previous best
    ## Run 227 stress 9.931894e-05 
    ## ... Procrustes: rmse 0.0002278598  max resid 0.0003633105 
    ## ... Similar to previous best
    ## Run 228 stress 9.968993e-05 
    ## ... Procrustes: rmse 0.0002214438  max resid 0.0003741664 
    ## ... Similar to previous best
    ## Run 229 stress 9.619244e-05 
    ## ... Procrustes: rmse 0.0001290242  max resid 0.0002789849 
    ## ... Similar to previous best
    ## Run 230 stress 9.808848e-05 
    ## ... Procrustes: rmse 0.0002031702  max resid 0.0003578444 
    ## ... Similar to previous best
    ## Run 231 stress 9.821931e-05 
    ## ... Procrustes: rmse 0.0002511759  max resid 0.0004086028 
    ## ... Similar to previous best
    ## Run 232 stress 9.951429e-05 
    ## ... Procrustes: rmse 0.0002596157  max resid 0.0004058327 
    ## ... Similar to previous best
    ## Run 233 stress 9.693323e-05 
    ## ... Procrustes: rmse 0.0001859194  max resid 0.0003019627 
    ## ... Similar to previous best
    ## Run 234 stress 9.990248e-05 
    ## ... Procrustes: rmse 0.0002306421  max resid 0.000376434 
    ## ... Similar to previous best
    ## Run 235 stress 9.973726e-05 
    ## ... Procrustes: rmse 0.0002496054  max resid 0.0004245176 
    ## ... Similar to previous best
    ## Run 236 stress 9.704271e-05 
    ## ... Procrustes: rmse 0.0002209923  max resid 0.000344993 
    ## ... Similar to previous best
    ## Run 237 stress 9.845931e-05 
    ## ... Procrustes: rmse 0.0002553161  max resid 0.0003988634 
    ## ... Similar to previous best
    ## Run 238 stress 9.63519e-05 
    ## ... Procrustes: rmse 0.0002428222  max resid 0.0003980219 
    ## ... Similar to previous best
    ## Run 239 stress 9.775959e-05 
    ## ... Procrustes: rmse 0.0002342981  max resid 0.0003919785 
    ## ... Similar to previous best
    ## Run 240 stress 9.310961e-05 
    ## ... Procrustes: rmse 0.0001947458  max resid 0.0003437659 
    ## ... Similar to previous best
    ## Run 241 stress 9.135578e-05 
    ## ... Procrustes: rmse 4.877951e-05  max resid 7.940018e-05 
    ## ... Similar to previous best
    ## Run 242 stress 9.895835e-05 
    ## ... Procrustes: rmse 0.0002565907  max resid 0.0004103433 
    ## ... Similar to previous best
    ## Run 243 stress 9.406866e-05 
    ## ... Procrustes: rmse 0.0001818652  max resid 0.0003235534 
    ## ... Similar to previous best
    ## Run 244 stress 9.766408e-05 
    ## ... Procrustes: rmse 0.000256624  max resid 0.0004085596 
    ## ... Similar to previous best
    ## Run 245 stress 9.632377e-05 
    ## ... Procrustes: rmse 0.0002492787  max resid 0.0003937942 
    ## ... Similar to previous best
    ## Run 246 stress 9.810484e-05 
    ## ... Procrustes: rmse 0.0002399288  max resid 0.0003742881 
    ## ... Similar to previous best
    ## Run 247 stress 9.752917e-05 
    ## ... Procrustes: rmse 0.0002126415  max resid 0.000371984 
    ## ... Similar to previous best
    ## Run 248 stress 9.938052e-05 
    ## ... Procrustes: rmse 0.0002597045  max resid 0.0004149922 
    ## ... Similar to previous best
    ## Run 249 stress 7.216972e-05 
    ## ... Procrustes: rmse 4.310102e-05  max resid 0.0001059448 
    ## ... Similar to previous best
    ## Run 250 stress 9.897518e-05 
    ## ... Procrustes: rmse 0.0002495635  max resid 0.0003894618 
    ## ... Similar to previous best
    ## Run 251 stress 6.985701e-05 
    ## ... Procrustes: rmse 3.777764e-05  max resid 8.084372e-05 
    ## ... Similar to previous best
    ## Run 252 stress 9.672054e-05 
    ## ... Procrustes: rmse 0.0001842517  max resid 0.0003376347 
    ## ... Similar to previous best
    ## Run 253 stress 9.706125e-05 
    ## ... Procrustes: rmse 0.0002148218  max resid 0.0003713098 
    ## ... Similar to previous best
    ## Run 254 stress 9.883724e-05 
    ## ... Procrustes: rmse 0.0002159007  max resid 0.0003698296 
    ## ... Similar to previous best
    ## Run 255 stress 8.494639e-05 
    ## ... Procrustes: rmse 6.886336e-05  max resid 0.0001165177 
    ## ... Similar to previous best
    ## Run 256 stress 9.652879e-05 
    ## ... Procrustes: rmse 0.0002009725  max resid 0.0003320783 
    ## ... Similar to previous best
    ## Run 257 stress 9.100415e-05 
    ## ... Procrustes: rmse 5.055186e-05  max resid 8.679131e-05 
    ## ... Similar to previous best
    ## Run 258 stress 8.12963e-05 
    ## ... Procrustes: rmse 4.505272e-05  max resid 7.032377e-05 
    ## ... Similar to previous best
    ## Run 259 stress 9.86915e-05 
    ## ... Procrustes: rmse 0.0002591758  max resid 0.0004106849 
    ## ... Similar to previous best
    ## Run 260 stress 9.809035e-05 
    ## ... Procrustes: rmse 0.0002276918  max resid 0.0003857547 
    ## ... Similar to previous best
    ## Run 261 stress 9.776294e-05 
    ## ... Procrustes: rmse 0.0001976249  max resid 0.0003468461 
    ## ... Similar to previous best
    ## Run 262 stress 9.987961e-05 
    ## ... Procrustes: rmse 0.0002549088  max resid 0.0004067623 
    ## ... Similar to previous best
    ## Run 263 stress 9.887045e-05 
    ## ... Procrustes: rmse 0.0002283928  max resid 0.0003768978 
    ## ... Similar to previous best
    ## Run 264 stress 8.455259e-05 
    ## ... Procrustes: rmse 4.816382e-05  max resid 0.0001017478 
    ## ... Similar to previous best
    ## Run 265 stress 8.838205e-05 
    ## ... Procrustes: rmse 4.266513e-05  max resid 9.575717e-05 
    ## ... Similar to previous best
    ## Run 266 stress 9.914479e-05 
    ## ... Procrustes: rmse 0.0002087559  max resid 0.0003715189 
    ## ... Similar to previous best
    ## Run 267 stress 9.378768e-05 
    ## ... Procrustes: rmse 0.0001635584  max resid 0.0002995387 
    ## ... Similar to previous best
    ## Run 268 stress 9.584004e-05 
    ## ... Procrustes: rmse 0.0001910578  max resid 0.0003271996 
    ## ... Similar to previous best
    ## Run 269 stress 8.645279e-05 
    ## ... Procrustes: rmse 5.45469e-05  max resid 7.800422e-05 
    ## ... Similar to previous best
    ## Run 270 stress 9.910645e-05 
    ## ... Procrustes: rmse 6.219354e-05  max resid 9.714718e-05 
    ## ... Similar to previous best
    ## Run 271 stress 9.843645e-05 
    ## ... Procrustes: rmse 0.0002561664  max resid 0.000410969 
    ## ... Similar to previous best
    ## Run 272 stress 7.832986e-05 
    ## ... Procrustes: rmse 3.851804e-05  max resid 6.59475e-05 
    ## ... Similar to previous best
    ## Run 273 stress 9.968374e-05 
    ## ... Procrustes: rmse 0.0002617922  max resid 0.0004170064 
    ## ... Similar to previous best
    ## Run 274 stress 9.791723e-05 
    ## ... Procrustes: rmse 0.0002291278  max resid 0.0003684419 
    ## ... Similar to previous best
    ## Run 275 stress 9.895323e-05 
    ## ... Procrustes: rmse 0.0002256872  max resid 0.0003880519 
    ## ... Similar to previous best
    ## Run 276 stress 9.807555e-05 
    ## ... Procrustes: rmse 0.0002538268  max resid 0.0003963096 
    ## ... Similar to previous best
    ## Run 277 stress 9.194916e-05 
    ## ... Procrustes: rmse 8.909138e-05  max resid 0.0001931996 
    ## ... Similar to previous best
    ## Run 278 stress 9.838688e-05 
    ## ... Procrustes: rmse 0.0002235698  max resid 0.0003836724 
    ## ... Similar to previous best
    ## Run 279 stress 9.587538e-05 
    ## ... Procrustes: rmse 0.0002482919  max resid 0.0003906365 
    ## ... Similar to previous best
    ## Run 280 stress 9.893121e-05 
    ## ... Procrustes: rmse 0.0002375193  max resid 0.0004151427 
    ## ... Similar to previous best
    ## Run 281 stress 9.556024e-05 
    ## ... Procrustes: rmse 0.0002036455  max resid 0.0003598435 
    ## ... Similar to previous best
    ## Run 282 stress 8.931964e-05 
    ## ... Procrustes: rmse 5.458771e-05  max resid 9.296608e-05 
    ## ... Similar to previous best
    ## Run 283 stress 9.555151e-05 
    ## ... Procrustes: rmse 0.0002440692  max resid 0.0003938616 
    ## ... Similar to previous best
    ## Run 284 stress 7.603978e-05 
    ## ... Procrustes: rmse 4.128816e-05  max resid 7.690653e-05 
    ## ... Similar to previous best
    ## Run 285 stress 9.786305e-05 
    ## ... Procrustes: rmse 0.0002378873  max resid 0.0003929922 
    ## ... Similar to previous best
    ## Run 286 stress 9.91605e-05 
    ## ... Procrustes: rmse 0.0002375896  max resid 0.0003848671 
    ## ... Similar to previous best
    ## Run 287 stress 9.688736e-05 
    ## ... Procrustes: rmse 0.0002238396  max resid 0.000353998 
    ## ... Similar to previous best
    ## Run 288 stress 9.977581e-05 
    ## ... Procrustes: rmse 0.0002621147  max resid 0.0004148341 
    ## ... Similar to previous best
    ## Run 289 stress 6.412084e-05 
    ## ... Procrustes: rmse 3.576786e-05  max resid 6.930417e-05 
    ## ... Similar to previous best
    ## Run 290 stress 9.821421e-05 
    ## ... Procrustes: rmse 0.0002556034  max resid 0.0004105271 
    ## ... Similar to previous best
    ## Run 291 stress 9.900336e-05 
    ## ... Procrustes: rmse 0.0002590359  max resid 0.0004157363 
    ## ... Similar to previous best
    ## Run 292 stress 9.885207e-05 
    ## ... Procrustes: rmse 0.0002536414  max resid 0.0004045532 
    ## ... Similar to previous best
    ## Run 293 stress 9.913707e-05 
    ## ... Procrustes: rmse 0.0002171288  max resid 0.0003495628 
    ## ... Similar to previous best
    ## Run 294 stress 9.503851e-05 
    ## ... Procrustes: rmse 0.0002486917  max resid 0.0003938096 
    ## ... Similar to previous best
    ## Run 295 stress 8.269657e-05 
    ## ... Procrustes: rmse 4.566351e-05  max resid 7.552136e-05 
    ## ... Similar to previous best
    ## Run 296 stress 7.636706e-05 
    ## ... Procrustes: rmse 4.762436e-05  max resid 9.632174e-05 
    ## ... Similar to previous best
    ## Run 297 stress 8.790118e-05 
    ## ... Procrustes: rmse 0.000141354  max resid 0.0002523989 
    ## ... Similar to previous best
    ## Run 298 stress 9.434842e-05 
    ## ... Procrustes: rmse 0.0002171319  max resid 0.0003410806 
    ## ... Similar to previous best
    ## Run 299 stress 9.944754e-05 
    ## ... Procrustes: rmse 0.000210033  max resid 0.000369925 
    ## ... Similar to previous best
    ## Run 300 stress 9.594373e-05 
    ## ... Procrustes: rmse 0.0001732165  max resid 0.0003210352 
    ## ... Similar to previous best
    ## Run 301 stress 9.526613e-05 
    ## ... Procrustes: rmse 0.000246272  max resid 0.00038664 
    ## ... Similar to previous best
    ## Run 302 stress 9.194971e-05 
    ## ... Procrustes: rmse 5.119166e-05  max resid 7.841509e-05 
    ## ... Similar to previous best
    ## Run 303 stress 9.843666e-05 
    ## ... Procrustes: rmse 0.0002445231  max resid 0.0003860587 
    ## ... Similar to previous best
    ## Run 304 stress 9.969238e-05 
    ## ... Procrustes: rmse 0.0002498332  max resid 0.0004099204 
    ## ... Similar to previous best
    ## Run 305 stress 9.440204e-05 
    ## ... Procrustes: rmse 4.381136e-05  max resid 8.854199e-05 
    ## ... Similar to previous best
    ## Run 306 stress 9.806624e-05 
    ## ... Procrustes: rmse 0.0002269182  max resid 0.0003654662 
    ## ... Similar to previous best
    ## Run 307 stress 9.819308e-05 
    ## ... Procrustes: rmse 0.0002501485  max resid 0.0003865005 
    ## ... Similar to previous best
    ## Run 308 stress 8.132226e-05 
    ## ... Procrustes: rmse 6.435443e-05  max resid 0.0001091431 
    ## ... Similar to previous best
    ## Run 309 stress 9.839274e-05 
    ## ... Procrustes: rmse 0.0002277508  max resid 0.0003636515 
    ## ... Similar to previous best
    ## Run 310 stress 9.493718e-05 
    ## ... Procrustes: rmse 0.0001120514  max resid 0.000227454 
    ## ... Similar to previous best
    ## Run 311 stress 9.868068e-05 
    ## ... Procrustes: rmse 0.0002450639  max resid 0.0003876691 
    ## ... Similar to previous best
    ## Run 312 stress 9.584973e-05 
    ## ... Procrustes: rmse 0.0002399429  max resid 0.0003769411 
    ## ... Similar to previous best
    ## Run 313 stress 9.815319e-05 
    ## ... Procrustes: rmse 0.0002571799  max resid 0.0004096779 
    ## ... Similar to previous best
    ## Run 314 stress 9.827629e-05 
    ## ... Procrustes: rmse 0.0002513841  max resid 0.0003909059 
    ## ... Similar to previous best
    ## Run 315 stress 9.378236e-05 
    ## ... Procrustes: rmse 0.0001998332  max resid 0.0003231298 
    ## ... Similar to previous best
    ## Run 316 stress 9.853213e-05 
    ## ... Procrustes: rmse 0.0002530859  max resid 0.0003948957 
    ## ... Similar to previous best
    ## Run 317 stress 9.863706e-05 
    ## ... Procrustes: rmse 0.0002413717  max resid 0.0004141781 
    ## ... Similar to previous best
    ## Run 318 stress 9.822736e-05 
    ## ... Procrustes: rmse 0.0002570174  max resid 0.0004098 
    ## ... Similar to previous best
    ## Run 319 stress 9.432937e-05 
    ## ... Procrustes: rmse 4.113924e-05  max resid 8.435654e-05 
    ## ... Similar to previous best
    ## Run 320 stress 9.685136e-05 
    ## ... Procrustes: rmse 0.0002223639  max resid 0.0003526838 
    ## ... Similar to previous best
    ## Run 321 stress 9.918327e-05 
    ## ... Procrustes: rmse 0.000211695  max resid 0.0003726272 
    ## ... Similar to previous best
    ## Run 322 stress 9.050156e-05 
    ## ... Procrustes: rmse 4.730566e-05  max resid 9.800077e-05 
    ## ... Similar to previous best
    ## Run 323 stress 9.539699e-05 
    ## ... Procrustes: rmse 0.0002059969  max resid 0.0003556039 
    ## ... Similar to previous best
    ## Run 324 stress 9.57909e-05 
    ## ... Procrustes: rmse 0.0002061289  max resid 0.0003595002 
    ## ... Similar to previous best
    ## Run 325 stress 9.830831e-05 
    ## ... Procrustes: rmse 0.0002562502  max resid 0.0004067326 
    ## ... Similar to previous best
    ## Run 326 stress 9.681402e-05 
    ## ... Procrustes: rmse 0.0002462261  max resid 0.0003862232 
    ## ... Similar to previous best
    ## Run 327 stress 8.118912e-05 
    ## ... Procrustes: rmse 4.630888e-05  max resid 9.263849e-05 
    ## ... Similar to previous best
    ## Run 328 stress 9.35895e-05 
    ## ... Procrustes: rmse 0.0002069424  max resid 0.0003297229 
    ## ... Similar to previous best
    ## Run 329 stress 9.531973e-05 
    ## ... Procrustes: rmse 0.0001243876  max resid 0.00020953 
    ## ... Similar to previous best
    ## Run 330 stress 9.98534e-05 
    ## ... Procrustes: rmse 0.0002147787  max resid 0.000373469 
    ## ... Similar to previous best
    ## Run 331 stress 9.709711e-05 
    ## ... Procrustes: rmse 0.0001904451  max resid 0.0003374799 
    ## ... Similar to previous best
    ## Run 332 stress 9.763831e-05 
    ## ... Procrustes: rmse 0.0002161717  max resid 0.0003485874 
    ## ... Similar to previous best
    ## Run 333 stress 9.377873e-05 
    ## ... Procrustes: rmse 0.0002195104  max resid 0.0003510572 
    ## ... Similar to previous best
    ## Run 334 stress 9.559118e-05 
    ## ... Procrustes: rmse 0.0002016504  max resid 0.0003217799 
    ## ... Similar to previous best
    ## Run 335 stress 9.995939e-05 
    ## ... Procrustes: rmse 0.0002586978  max resid 0.000404355 
    ## ... Similar to previous best
    ## Run 336 stress 9.974029e-05 
    ## ... Procrustes: rmse 0.0002197847  max resid 0.0003474101 
    ## ... Similar to previous best
    ## Run 337 stress 9.794576e-05 
    ## ... Procrustes: rmse 0.0002216025  max resid 0.0003596393 
    ## ... Similar to previous best
    ## Run 338 stress 9.662938e-05 
    ## ... Procrustes: rmse 0.0001288098  max resid 0.0002452351 
    ## ... Similar to previous best
    ## Run 339 stress 6.635974e-05 
    ## ... Procrustes: rmse 3.833299e-05  max resid 6.679557e-05 
    ## ... Similar to previous best
    ## Run 340 stress 9.414683e-05 
    ## ... Procrustes: rmse 0.0002068201  max resid 0.0003564177 
    ## ... Similar to previous best
    ## Run 341 stress 9.796039e-05 
    ## ... Procrustes: rmse 0.0002436547  max resid 0.0004125616 
    ## ... Similar to previous best
    ## Run 342 stress 6.332738e-05 
    ## ... Procrustes: rmse 3.422035e-05  max resid 6.365661e-05 
    ## ... Similar to previous best
    ## Run 343 stress 9.354459e-05 
    ## ... Procrustes: rmse 7.760576e-05  max resid 0.0001152975 
    ## ... Similar to previous best
    ## Run 344 stress 8.755254e-05 
    ## ... Procrustes: rmse 5.481905e-05  max resid 8.483081e-05 
    ## ... Similar to previous best
    ## Run 345 stress 8.821303e-05 
    ## ... Procrustes: rmse 7.851309e-05  max resid 0.0001040934 
    ## ... Similar to previous best
    ## Run 346 stress 9.357969e-05 
    ## ... Procrustes: rmse 0.00018427  max resid 0.0003113932 
    ## ... Similar to previous best
    ## Run 347 stress 9.863349e-05 
    ## ... Procrustes: rmse 0.0002584758  max resid 0.0004128924 
    ## ... Similar to previous best
    ## Run 348 stress 9.927807e-05 
    ## ... Procrustes: rmse 0.0001856196  max resid 0.0003508291 
    ## ... Similar to previous best
    ## Run 349 stress 6.21068e-05 
    ## ... Procrustes: rmse 4.092413e-05  max resid 7.088835e-05 
    ## ... Similar to previous best
    ## Run 350 stress 7.94589e-05 
    ## ... Procrustes: rmse 4.367325e-05  max resid 8.522421e-05 
    ## ... Similar to previous best
    ## Run 351 stress 9.95519e-05 
    ## ... Procrustes: rmse 0.0002604286  max resid 0.000412292 
    ## ... Similar to previous best
    ## Run 352 stress 8.401085e-05 
    ## ... Procrustes: rmse 5.125279e-05  max resid 0.0001144213 
    ## ... Similar to previous best
    ## Run 353 stress 9.85631e-05 
    ## ... Procrustes: rmse 0.0002515042  max resid 0.0003878317 
    ## ... Similar to previous best
    ## Run 354 stress 9.732976e-05 
    ## ... Procrustes: rmse 0.0002553991  max resid 0.0004065447 
    ## ... Similar to previous best
    ## Run 355 stress 6.286469e-05 
    ## ... Procrustes: rmse 3.50828e-05  max resid 6.741769e-05 
    ## ... Similar to previous best
    ## Run 356 stress 9.825206e-05 
    ## ... Procrustes: rmse 0.000208962  max resid 0.0003436237 
    ## ... Similar to previous best
    ## Run 357 stress 9.836492e-05 
    ## ... Procrustes: rmse 0.0002565696  max resid 0.0004104601 
    ## ... Similar to previous best
    ## Run 358 stress 3.898412e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 2.973122e-05  max resid 5.678445e-05 
    ## ... Similar to previous best
    ## Run 359 stress 6.983051e-05 
    ## ... Procrustes: rmse 9.87237e-05  max resid 0.0001608677 
    ## ... Similar to previous best
    ## Run 360 stress 9.820795e-05 
    ## ... Procrustes: rmse 0.0001935836  max resid 0.0003327237 
    ## ... Similar to previous best
    ## Run 361 stress 9.40573e-05 
    ## ... Procrustes: rmse 0.0001824819  max resid 0.000334103 
    ## ... Similar to previous best
    ## Run 362 stress 9.875959e-05 
    ## ... Procrustes: rmse 0.0001999475  max resid 0.0003134583 
    ## ... Similar to previous best
    ## Run 363 stress 9.72858e-05 
    ## ... Procrustes: rmse 0.0002073082  max resid 0.0003193793 
    ## ... Similar to previous best
    ## Run 364 stress 9.601673e-05 
    ## ... Procrustes: rmse 0.0002195362  max resid 0.0003795172 
    ## ... Similar to previous best
    ## Run 365 stress 4.271427e-05 
    ## ... Procrustes: rmse 2.626047e-05  max resid 5.381251e-05 
    ## ... Similar to previous best
    ## Run 366 stress 9.995014e-05 
    ## ... Procrustes: rmse 0.0001802465  max resid 0.0003294335 
    ## ... Similar to previous best
    ## Run 367 stress 8.674334e-05 
    ## ... Procrustes: rmse 6.625669e-05  max resid 0.0001094609 
    ## ... Similar to previous best
    ## Run 368 stress 9.752298e-05 
    ## ... Procrustes: rmse 0.0002191696  max resid 0.0003826506 
    ## ... Similar to previous best
    ## Run 369 stress 4.565137e-05 
    ## ... Procrustes: rmse 2.841588e-05  max resid 8.6863e-05 
    ## ... Similar to previous best
    ## Run 370 stress 9.959694e-05 
    ## ... Procrustes: rmse 0.0002061633  max resid 0.0003354704 
    ## ... Similar to previous best
    ## Run 371 stress 9.93171e-05 
    ## ... Procrustes: rmse 0.0002043915  max resid 0.0003579105 
    ## ... Similar to previous best
    ## Run 372 stress 9.96707e-05 
    ## ... Procrustes: rmse 0.000212731  max resid 0.0003534167 
    ## ... Similar to previous best
    ## Run 373 stress 9.967107e-05 
    ## ... Procrustes: rmse 0.0002090156  max resid 0.0003578787 
    ## ... Similar to previous best
    ## Run 374 stress 9.974312e-05 
    ## ... Procrustes: rmse 0.0002068548  max resid 0.0003413558 
    ## ... Similar to previous best
    ## Run 375 stress 9.603204e-05 
    ## ... Procrustes: rmse 9.848194e-05  max resid 0.0001703346 
    ## ... Similar to previous best
    ## Run 376 stress 7.163512e-05 
    ## ... Procrustes: rmse 4.716338e-05  max resid 7.868327e-05 
    ## ... Similar to previous best
    ## Run 377 stress 8.109704e-05 
    ## ... Procrustes: rmse 3.957789e-05  max resid 7.466451e-05 
    ## ... Similar to previous best
    ## Run 378 stress 9.715757e-05 
    ## ... Procrustes: rmse 4.120863e-05  max resid 7.266481e-05 
    ## ... Similar to previous best
    ## Run 379 stress 9.929991e-05 
    ## ... Procrustes: rmse 0.0001804958  max resid 0.000283589 
    ## ... Similar to previous best
    ## Run 380 stress 9.804596e-05 
    ## ... Procrustes: rmse 0.000237858  max resid 0.0003843869 
    ## ... Similar to previous best
    ## Run 381 stress 6.12543e-05 
    ## ... Procrustes: rmse 2.835646e-05  max resid 7.249092e-05 
    ## ... Similar to previous best
    ## Run 382 stress 9.321683e-05 
    ## ... Procrustes: rmse 0.0001942784  max resid 0.0003361352 
    ## ... Similar to previous best
    ## Run 383 stress 9.874209e-05 
    ## ... Procrustes: rmse 0.0002393066  max resid 0.0003824598 
    ## ... Similar to previous best
    ## Run 384 stress 9.766672e-05 
    ## ... Procrustes: rmse 0.0002363463  max resid 0.0003869235 
    ## ... Similar to previous best
    ## Run 385 stress 8.948667e-05 
    ## ... Procrustes: rmse 4.802595e-05  max resid 7.910362e-05 
    ## ... Similar to previous best
    ## Run 386 stress 9.944108e-05 
    ## ... Procrustes: rmse 0.0002257895  max resid 0.0003756918 
    ## ... Similar to previous best
    ## Run 387 stress 6.328329e-05 
    ## ... Procrustes: rmse 5.029901e-05  max resid 0.0001129264 
    ## ... Similar to previous best
    ## Run 388 stress 9.931224e-05 
    ## ... Procrustes: rmse 0.0002044267  max resid 0.0003298569 
    ## ... Similar to previous best
    ## Run 389 stress 2.975718e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 2.489187e-05  max resid 6.203058e-05 
    ## ... Similar to previous best
    ## Run 390 stress 9.387132e-05 
    ## ... Procrustes: rmse 0.0001851174  max resid 0.000305304 
    ## ... Similar to previous best
    ## Run 391 stress 9.443006e-05 
    ## ... Procrustes: rmse 0.0001872786  max resid 0.0003149923 
    ## ... Similar to previous best
    ## Run 392 stress 9.749517e-05 
    ## ... Procrustes: rmse 0.000246296  max resid 0.0004153733 
    ## ... Similar to previous best
    ## Run 393 stress 5.709178e-05 
    ## ... Procrustes: rmse 3.188925e-05  max resid 6.336118e-05 
    ## ... Similar to previous best
    ## Run 394 stress 9.716347e-05 
    ## ... Procrustes: rmse 0.0002470169  max resid 0.0004169546 
    ## ... Similar to previous best
    ## Run 395 stress 6.722047e-05 
    ## ... Procrustes: rmse 6.47397e-05  max resid 9.118088e-05 
    ## ... Similar to previous best
    ## Run 396 stress 9.712817e-05 
    ## ... Procrustes: rmse 0.0002018999  max resid 0.0003377425 
    ## ... Similar to previous best
    ## Run 397 stress 9.66297e-05 
    ## ... Procrustes: rmse 0.0002170374  max resid 0.0003710363 
    ## ... Similar to previous best
    ## Run 398 stress 9.775064e-05 
    ## ... Procrustes: rmse 0.000245605  max resid 0.0004123119 
    ## ... Similar to previous best
    ## Run 399 stress 9.948849e-05 
    ## ... Procrustes: rmse 0.0002496691  max resid 0.0004211923 
    ## ... Similar to previous best
    ## Run 400 stress 9.324026e-05 
    ## ... Procrustes: rmse 0.000179638  max resid 0.0003167544 
    ## ... Similar to previous best
    ## Run 401 stress 9.996766e-05 
    ## ... Procrustes: rmse 0.0002039945  max resid 0.0003475346 
    ## ... Similar to previous best
    ## Run 402 stress 9.930512e-05 
    ## ... Procrustes: rmse 0.0002166605  max resid 0.000365244 
    ## ... Similar to previous best
    ## Run 403 stress 9.166732e-05 
    ## ... Procrustes: rmse 5.222422e-05  max resid 0.0001041485 
    ## ... Similar to previous best
    ## Run 404 stress 9.81818e-05 
    ## ... Procrustes: rmse 0.0002509193  max resid 0.0004229662 
    ## ... Similar to previous best
    ## Run 405 stress 9.952822e-05 
    ## ... Procrustes: rmse 0.0002504433  max resid 0.000412033 
    ## ... Similar to previous best
    ## Run 406 stress 9.911254e-05 
    ## ... Procrustes: rmse 0.000221127  max resid 0.000379398 
    ## ... Similar to previous best
    ## Run 407 stress 9.854961e-05 
    ## ... Procrustes: rmse 4.775258e-05  max resid 7.736943e-05 
    ## ... Similar to previous best
    ## Run 408 stress 8.129409e-05 
    ## ... Procrustes: rmse 5.319448e-05  max resid 0.0001103499 
    ## ... Similar to previous best
    ## Run 409 stress 9.885116e-05 
    ## ... Procrustes: rmse 0.0001928626  max resid 0.0003295514 
    ## ... Similar to previous best
    ## Run 410 stress 9.828978e-05 
    ## ... Procrustes: rmse 0.0002201262  max resid 0.0003744225 
    ## ... Similar to previous best
    ## Run 411 stress 9.840705e-05 
    ## ... Procrustes: rmse 0.0002515135  max resid 0.0004214365 
    ## ... Similar to previous best
    ## Run 412 stress 8.145092e-05 
    ## ... Procrustes: rmse 3.779898e-05  max resid 7.738884e-05 
    ## ... Similar to previous best
    ## Run 413 stress 6.244156e-05 
    ## ... Procrustes: rmse 3.052932e-05  max resid 5.970629e-05 
    ## ... Similar to previous best
    ## Run 414 stress 9.960645e-05 
    ## ... Procrustes: rmse 0.000252587  max resid 0.0004278705 
    ## ... Similar to previous best
    ## Run 415 stress 8.955105e-05 
    ## ... Procrustes: rmse 7.272258e-05  max resid 0.0001096382 
    ## ... Similar to previous best
    ## Run 416 stress 9.921214e-05 
    ## ... Procrustes: rmse 0.0001687048  max resid 0.0002741119 
    ## ... Similar to previous best
    ## Run 417 stress 7.353601e-05 
    ## ... Procrustes: rmse 0.0001304885  max resid 0.0002186399 
    ## ... Similar to previous best
    ## Run 418 stress 9.884193e-05 
    ## ... Procrustes: rmse 0.0002212771  max resid 0.0003772641 
    ## ... Similar to previous best
    ## Run 419 stress 9.923995e-05 
    ## ... Procrustes: rmse 0.0002206094  max resid 0.0003753817 
    ## ... Similar to previous best
    ## Run 420 stress 9.486044e-05 
    ## ... Procrustes: rmse 4.878815e-05  max resid 7.099748e-05 
    ## ... Similar to previous best
    ## Run 421 stress 9.451283e-05 
    ## ... Procrustes: rmse 0.0002403647  max resid 0.0004070855 
    ## ... Similar to previous best
    ## Run 422 stress 9.928186e-05 
    ## ... Procrustes: rmse 0.0002541223  max resid 0.0004252375 
    ## ... Similar to previous best
    ## Run 423 stress 9.337101e-05 
    ## ... Procrustes: rmse 0.0002378144  max resid 0.0003995249 
    ## ... Similar to previous best
    ## Run 424 stress 9.953691e-05 
    ## ... Procrustes: rmse 0.0002052937  max resid 0.0003440173 
    ## ... Similar to previous best
    ## Run 425 stress 9.800431e-05 
    ## ... Procrustes: rmse 0.0002480368  max resid 0.0004155767 
    ## ... Similar to previous best
    ## Run 426 stress 7.570459e-05 
    ## ... Procrustes: rmse 4.402449e-05  max resid 7.499852e-05 
    ## ... Similar to previous best
    ## Run 427 stress 9.263616e-05 
    ## ... Procrustes: rmse 0.0001982374  max resid 0.0003339762 
    ## ... Similar to previous best
    ## Run 428 stress 9.923797e-05 
    ## ... Procrustes: rmse 0.000211155  max resid 0.0003565369 
    ## ... Similar to previous best
    ## Run 429 stress 9.695991e-05 
    ## ... Procrustes: rmse 0.0002483057  max resid 0.0004148675 
    ## ... Similar to previous best
    ## Run 430 stress 9.96447e-05 
    ## ... Procrustes: rmse 0.000252091  max resid 0.0004232013 
    ## ... Similar to previous best
    ## Run 431 stress 3.95404e-05 
    ## ... Procrustes: rmse 2.194811e-05  max resid 4.33618e-05 
    ## ... Similar to previous best
    ## Run 432 stress 9.793267e-05 
    ## ... Procrustes: rmse 0.0001962536  max resid 0.0003324569 
    ## ... Similar to previous best
    ## Run 433 stress 5.380587e-05 
    ## ... Procrustes: rmse 3.916424e-05  max resid 6.910213e-05 
    ## ... Similar to previous best
    ## Run 434 stress 6.416866e-05 
    ## ... Procrustes: rmse 4.587228e-05  max resid 7.761023e-05 
    ## ... Similar to previous best
    ## Run 435 stress 9.745229e-05 
    ## ... Procrustes: rmse 0.0002119879  max resid 0.0003550461 
    ## ... Similar to previous best
    ## Run 436 stress 9.659655e-05 
    ## ... Procrustes: rmse 0.0002468339  max resid 0.0004137942 
    ## ... Similar to previous best
    ## Run 437 stress 9.896646e-05 
    ## ... Procrustes: rmse 0.0002253871  max resid 0.0003826286 
    ## ... Similar to previous best
    ## Run 438 stress 9.962905e-05 
    ## ... Procrustes: rmse 0.0002137665  max resid 0.000352467 
    ## ... Similar to previous best
    ## Run 439 stress 9.982923e-05 
    ## ... Procrustes: rmse 0.000245778  max resid 0.0004133652 
    ## ... Similar to previous best
    ## Run 440 stress 7.840087e-05 
    ## ... Procrustes: rmse 7.22349e-05  max resid 0.0001196368 
    ## ... Similar to previous best
    ## Run 441 stress 9.439802e-05 
    ## ... Procrustes: rmse 0.0002372916  max resid 0.0003961939 
    ## ... Similar to previous best
    ## Run 442 stress 9.987732e-05 
    ## ... Procrustes: rmse 0.0002477143  max resid 0.0004114352 
    ## ... Similar to previous best
    ## Run 443 stress 9.922931e-05 
    ## ... Procrustes: rmse 0.0002224039  max resid 0.0003843065 
    ## ... Similar to previous best
    ## Run 444 stress 9.245192e-05 
    ## ... Procrustes: rmse 0.0002288107  max resid 0.0003831253 
    ## ... Similar to previous best
    ## Run 445 stress 9.764825e-05 
    ## ... Procrustes: rmse 0.0002037413  max resid 0.0003324878 
    ## ... Similar to previous best
    ## Run 446 stress 5.797676e-05 
    ## ... Procrustes: rmse 2.959617e-05  max resid 5.494554e-05 
    ## ... Similar to previous best
    ## Run 447 stress 9.76763e-05 
    ## ... Procrustes: rmse 0.0002139285  max resid 0.0003620727 
    ## ... Similar to previous best
    ## Run 448 stress 9.681329e-05 
    ## ... Procrustes: rmse 0.0002488459  max resid 0.0004162669 
    ## ... Similar to previous best
    ## Run 449 stress 6.402985e-05 
    ## ... Procrustes: rmse 5.715269e-05  max resid 0.0001032007 
    ## ... Similar to previous best
    ## Run 450 stress 6.843393e-05 
    ## ... Procrustes: rmse 3.732714e-05  max resid 5.555581e-05 
    ## ... Similar to previous best
    ## Run 451 stress 9.985786e-05 
    ## ... Procrustes: rmse 0.0002247112  max resid 0.0003855189 
    ## ... Similar to previous best
    ## Run 452 stress 9.336598e-05 
    ## ... Procrustes: rmse 0.0002026222  max resid 0.0003376205 
    ## ... Similar to previous best
    ## Run 453 stress 9.706198e-05 
    ## ... Procrustes: rmse 0.0002161214  max resid 0.0003685976 
    ## ... Similar to previous best
    ## Run 454 stress 9.970956e-05 
    ## ... Procrustes: rmse 0.0002285646  max resid 0.0003791482 
    ## ... Similar to previous best
    ## Run 455 stress 9.889198e-05 
    ## ... Procrustes: rmse 0.0002164948  max resid 0.0003669537 
    ## ... Similar to previous best
    ## Run 456 stress 9.026447e-05 
    ## ... Procrustes: rmse 0.0001045526  max resid 0.0001506862 
    ## ... Similar to previous best
    ## Run 457 stress 9.110549e-05 
    ## ... Procrustes: rmse 4.936542e-05  max resid 0.0001012212 
    ## ... Similar to previous best
    ## Run 458 stress 9.943514e-05 
    ## ... Procrustes: rmse 0.00024156  max resid 0.0003987503 
    ## ... Similar to previous best
    ## Run 459 stress 8.892058e-05 
    ## ... Procrustes: rmse 0.0001179275  max resid 0.0001824656 
    ## ... Similar to previous best
    ## Run 460 stress 9.541404e-05 
    ## ... Procrustes: rmse 0.0001609789  max resid 0.0002736243 
    ## ... Similar to previous best
    ## Run 461 stress 9.938121e-05 
    ## ... Procrustes: rmse 0.0002117903  max resid 0.0003549756 
    ## ... Similar to previous best
    ## Run 462 stress 9.680013e-05 
    ## ... Procrustes: rmse 0.0002051971  max resid 0.0003488591 
    ## ... Similar to previous best
    ## Run 463 stress 9.97619e-05 
    ## ... Procrustes: rmse 0.0002549688  max resid 0.0004281893 
    ## ... Similar to previous best
    ## Run 464 stress 8.754873e-05 
    ## ... Procrustes: rmse 4.389708e-05  max resid 9.21111e-05 
    ## ... Similar to previous best
    ## Run 465 stress 8.589335e-05 
    ## ... Procrustes: rmse 5.324823e-05  max resid 0.0001040683 
    ## ... Similar to previous best
    ## Run 466 stress 6.393184e-05 
    ## ... Procrustes: rmse 7.848365e-05  max resid 0.0001387359 
    ## ... Similar to previous best
    ## Run 467 stress 9.858293e-05 
    ## ... Procrustes: rmse 0.0002514706  max resid 0.0004223967 
    ## ... Similar to previous best
    ## Run 468 stress 9.987723e-05 
    ## ... Procrustes: rmse 0.0002050855  max resid 0.00034342 
    ## ... Similar to previous best
    ## Run 469 stress 9.849565e-05 
    ## ... Procrustes: rmse 0.0002470361  max resid 0.0004146472 
    ## ... Similar to previous best
    ## Run 470 stress 9.725562e-05 
    ## ... Procrustes: rmse 0.0001917675  max resid 0.0003249939 
    ## ... Similar to previous best
    ## Run 471 stress 9.83286e-05 
    ## ... Procrustes: rmse 0.0002057011  max resid 0.0003442611 
    ## ... Similar to previous best
    ## Run 472 stress 9.849389e-05 
    ## ... Procrustes: rmse 0.0002534298  max resid 0.0004247739 
    ## ... Similar to previous best
    ## Run 473 stress 6.606106e-05 
    ## ... Procrustes: rmse 3.047558e-05  max resid 5.386715e-05 
    ## ... Similar to previous best
    ## Run 474 stress 9.779622e-05 
    ## ... Procrustes: rmse 7.236378e-05  max resid 0.0001123865 
    ## ... Similar to previous best
    ## Run 475 stress 9.285761e-05 
    ## ... Procrustes: rmse 0.0001941038  max resid 0.0003302057 
    ## ... Similar to previous best
    ## Run 476 stress 9.875271e-05 
    ## ... Procrustes: rmse 0.0002514884  max resid 0.0004226048 
    ## ... Similar to previous best
    ## Run 477 stress 8.051405e-05 
    ## ... Procrustes: rmse 0.0001237164  max resid 0.0002046147 
    ## ... Similar to previous best
    ## Run 478 stress 7.666153e-05 
    ## ... Procrustes: rmse 9.371349e-05  max resid 0.0001687111 
    ## ... Similar to previous best
    ## Run 479 stress 9.675059e-05 
    ## ... Procrustes: rmse 0.0002480478  max resid 0.0004136335 
    ## ... Similar to previous best
    ## Run 480 stress 9.975949e-05 
    ## ... Procrustes: rmse 0.0002072841  max resid 0.0003507996 
    ## ... Similar to previous best
    ## Run 481 stress 8.064642e-05 
    ## ... Procrustes: rmse 5.870499e-05  max resid 0.0001089232 
    ## ... Similar to previous best
    ## Run 482 stress 9.953101e-05 
    ## ... Procrustes: rmse 0.0002494991  max resid 0.0004177111 
    ## ... Similar to previous best
    ## Run 483 stress 9.904152e-05 
    ## ... Procrustes: rmse 5.172552e-05  max resid 0.0001020164 
    ## ... Similar to previous best
    ## Run 484 stress 9.864017e-05 
    ## ... Procrustes: rmse 0.000212814  max resid 0.0003588876 
    ## ... Similar to previous best
    ## Run 485 stress 9.877935e-05 
    ## ... Procrustes: rmse 0.0002122975  max resid 0.0003497095 
    ## ... Similar to previous best
    ## Run 486 stress 9.541756e-05 
    ## ... Procrustes: rmse 4.766414e-05  max resid 9.825508e-05 
    ## ... Similar to previous best
    ## Run 487 stress 8.752218e-05 
    ## ... Procrustes: rmse 4.22785e-05  max resid 8.184347e-05 
    ## ... Similar to previous best
    ## Run 488 stress 9.890379e-05 
    ## ... Procrustes: rmse 0.0002200443  max resid 0.0003626134 
    ## ... Similar to previous best
    ## Run 489 stress 9.830754e-05 
    ## ... Procrustes: rmse 0.0002382614  max resid 0.0003938003 
    ## ... Similar to previous best
    ## Run 490 stress 7.307887e-05 
    ## ... Procrustes: rmse 9.416633e-05  max resid 0.0001707495 
    ## ... Similar to previous best
    ## Run 491 stress 9.659093e-05 
    ## ... Procrustes: rmse 0.0002196565  max resid 0.0003637991 
    ## ... Similar to previous best
    ## Run 492 stress 8.964109e-05 
    ## ... Procrustes: rmse 4.808695e-05  max resid 7.395883e-05 
    ## ... Similar to previous best
    ## Run 493 stress 9.935193e-05 
    ## ... Procrustes: rmse 0.0001682401  max resid 0.0002780466 
    ## ... Similar to previous best
    ## Run 494 stress 9.731643e-05 
    ## ... Procrustes: rmse 0.0002132475  max resid 0.0003656811 
    ## ... Similar to previous best
    ## Run 495 stress 7.237203e-05 
    ## ... Procrustes: rmse 6.894397e-05  max resid 0.0001183757 
    ## ... Similar to previous best
    ## Run 496 stress 9.403732e-05 
    ## ... Procrustes: rmse 0.0002064041  max resid 0.0003490438 
    ## ... Similar to previous best
    ## Run 497 stress 9.814009e-05 
    ## ... Procrustes: rmse 0.0002046796  max resid 0.0003390623 
    ## ... Similar to previous best
    ## Run 498 stress 9.821707e-05 
    ## ... Procrustes: rmse 0.0002191982  max resid 0.0003739357 
    ## ... Similar to previous best
    ## Run 499 stress 9.775248e-05 
    ## ... Procrustes: rmse 0.0001772656  max resid 0.0002997919 
    ## ... Similar to previous best
    ## Run 500 stress 9.797504e-05 
    ## ... Procrustes: rmse 0.000218536  max resid 0.0003713855 
    ## ... Similar to previous best
    ## *** Best solution repeated 112 times

    ## Warning in metaMDS(PD_beta_ref_dist$Btotal, try = 1000, parallel = 4, trymax =
    ## 500, : stress is (nearly) zero: you may have insufficient data

``` r
# Surveyed sites 
PD_beta_NMDS <- metaMDS(PD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07927535 
    ## Run 1 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237959  max resid 0.04977607 
    ## Run 2 stress 0.09106036 
    ## Run 3 stress 0.0910604 
    ## Run 4 stress 0.09544379 
    ## Run 5 stress 0.07927536 
    ## ... Procrustes: rmse 1.885297e-05  max resid 5.419731e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.07936952 
    ## ... Procrustes: rmse 0.009045158  max resid 0.03248274 
    ## Run 7 stress 0.09544373 
    ## Run 8 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044085  max resid 0.03246691 
    ## Run 9 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733521  max resid 0.05639641 
    ## Run 10 stress 0.09247106 
    ## Run 11 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 3.199107e-05  max resid 7.886017e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.09072905 
    ## Run 13 stress 0.07936951 
    ## ... Procrustes: rmse 0.009036365  max resid 0.03240702 
    ## Run 14 stress 0.104507 
    ## Run 15 stress 0.09072908 
    ## Run 16 stress 0.09106032 
    ## Run 17 stress 0.0793695 
    ## ... Procrustes: rmse 0.009050272  max resid 0.03249675 
    ## Run 18 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173385  max resid 0.05632909 
    ## Run 19 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237996  max resid 0.04972887 
    ## Run 20 stress 0.08985608 
    ## Run 21 stress 0.09474065 
    ## Run 22 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237501  max resid 0.04968757 
    ## Run 23 stress 0.07927534 
    ## ... Procrustes: rmse 1.955085e-05  max resid 5.019943e-05 
    ## ... Similar to previous best
    ## Run 24 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 8.488353e-06  max resid 1.854055e-05 
    ## ... Similar to previous best
    ## Run 25 stress 0.09247103 
    ## Run 26 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734847  max resid 0.05639804 
    ## Run 27 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732544  max resid 0.05623741 
    ## Run 28 stress 0.09296227 
    ## Run 29 stress 0.07927534 
    ## ... Procrustes: rmse 2.760286e-05  max resid 7.405225e-05 
    ## ... Similar to previous best
    ## Run 30 stress 0.07927535 
    ## ... Procrustes: rmse 4.354485e-05  max resid 9.95942e-05 
    ## ... Similar to previous best
    ## Run 31 stress 0.07951445 
    ## ... Procrustes: rmse 0.01737968  max resid 0.05660388 
    ## Run 32 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 5.895311e-06  max resid 1.797396e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.09252224 
    ## Run 34 stress 0.07927538 
    ## ... Procrustes: rmse 4.782199e-05  max resid 0.0001266252 
    ## ... Similar to previous best
    ## Run 35 stress 0.07927537 
    ## ... Procrustes: rmse 4.865013e-05  max resid 0.0001333271 
    ## ... Similar to previous best
    ## Run 36 stress 0.0793695 
    ## ... Procrustes: rmse 0.0090488  max resid 0.03248394 
    ## Run 37 stress 0.09247106 
    ## Run 38 stress 0.07951446 
    ## ... Procrustes: rmse 0.0173233  max resid 0.05620164 
    ## Run 39 stress 0.08985604 
    ## Run 40 stress 0.0954437 
    ## Run 41 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734114  max resid 0.05634499 
    ## Run 42 stress 0.09022792 
    ## Run 43 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237082  max resid 0.04963155 
    ## Run 44 stress 0.07936952 
    ## ... Procrustes: rmse 0.009048128  max resid 0.03248825 
    ## Run 45 stress 0.07936955 
    ## ... Procrustes: rmse 0.009025463  max resid 0.03234235 
    ## Run 46 stress 0.107482 
    ## Run 47 stress 0.07969259 
    ## ... Procrustes: rmse 0.0123883  max resid 0.04978108 
    ## Run 48 stress 0.09022801 
    ## Run 49 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237221  max resid 0.04967661 
    ## Run 50 stress 0.09290699 
    ## Run 51 stress 0.09072904 
    ## Run 52 stress 0.09043759 
    ## Run 53 stress 0.0954437 
    ## Run 54 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734381  max resid 0.05636736 
    ## Run 55 stress 0.07927534 
    ## ... Procrustes: rmse 6.470058e-06  max resid 1.554745e-05 
    ## ... Similar to previous best
    ## Run 56 stress 0.08985603 
    ## Run 57 stress 0.07936952 
    ## ... Procrustes: rmse 0.009034601  max resid 0.03239917 
    ## Run 58 stress 0.09388344 
    ## Run 59 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734366  max resid 0.05637546 
    ## Run 60 stress 0.07927534 
    ## ... Procrustes: rmse 1.225935e-05  max resid 3.329631e-05 
    ## ... Similar to previous best
    ## Run 61 stress 0.1033494 
    ## Run 62 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904975  max resid 0.03249383 
    ## Run 63 stress 0.09247105 
    ## Run 64 stress 0.07936953 
    ## ... Procrustes: rmse 0.009032802  max resid 0.03238691 
    ## Run 65 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735226  max resid 0.05641259 
    ## Run 66 stress 0.09043737 
    ## Run 67 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734861  max resid 0.05639089 
    ## Run 68 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045744  max resid 0.03247433 
    ## Run 69 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237318  max resid 0.04968131 
    ## Run 70 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239889  max resid 0.04985702 
    ## Run 71 stress 0.09106042 
    ## Run 72 stress 0.09022776 
    ## Run 73 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735738  max resid 0.05644969 
    ## Run 74 stress 0.07969263 
    ## ... Procrustes: rmse 0.01236896  max resid 0.0495945 
    ## Run 75 stress 0.09072905 
    ## Run 76 stress 0.0898561 
    ## Run 77 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237366  max resid 0.04968158 
    ## Run 78 stress 0.08985618 
    ## Run 79 stress 0.0947401 
    ## Run 80 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734213  max resid 0.05634827 
    ## Run 81 stress 0.07927536 
    ## ... Procrustes: rmse 4.217155e-05  max resid 0.0001040375 
    ## ... Similar to previous best
    ## Run 82 stress 0.09290698 
    ## Run 83 stress 0.07936951 
    ## ... Procrustes: rmse 0.009036779  max resid 0.0324117 
    ## Run 84 stress 0.07927534 
    ## ... Procrustes: rmse 1.894326e-05  max resid 5.136993e-05 
    ## ... Similar to previous best
    ## Run 85 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734269  max resid 0.05640671 
    ## Run 86 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735112  max resid 0.05640159 
    ## Run 87 stress 0.07927534 
    ## ... Procrustes: rmse 6.28337e-06  max resid 1.92554e-05 
    ## ... Similar to previous best
    ## Run 88 stress 0.09022783 
    ## Run 89 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237859  max resid 0.04973951 
    ## Run 90 stress 0.09072895 
    ## Run 91 stress 0.1112806 
    ## Run 92 stress 0.07927534 
    ## ... Procrustes: rmse 1.004841e-05  max resid 2.977543e-05 
    ## ... Similar to previous best
    ## Run 93 stress 0.1047632 
    ## Run 94 stress 0.09043754 
    ## Run 95 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237313  max resid 0.04964536 
    ## Run 96 stress 0.0910606 
    ## Run 97 stress 0.07936956 
    ## ... Procrustes: rmse 0.00902364  max resid 0.03233103 
    ## Run 98 stress 0.09072897 
    ## Run 99 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046874  max resid 0.03247956 
    ## Run 100 stress 0.07927534 
    ## ... Procrustes: rmse 6.992686e-06  max resid 1.796188e-05 
    ## ... Similar to previous best
    ## Run 101 stress 0.07969262 
    ## ... Procrustes: rmse 0.01239952  max resid 0.04985203 
    ## Run 102 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045988  max resid 0.03246862 
    ## Run 103 stress 0.09022789 
    ## Run 104 stress 0.07927537 
    ## ... Procrustes: rmse 6.904247e-05  max resid 0.0001945992 
    ## ... Similar to previous best
    ## Run 105 stress 0.09290689 
    ## Run 106 stress 0.08985605 
    ## Run 107 stress 0.09544394 
    ## Run 108 stress 0.09043765 
    ## Run 109 stress 0.0910603 
    ## Run 110 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237978  max resid 0.04968474 
    ## Run 111 stress 0.09382064 
    ## Run 112 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046101  max resid 0.03247315 
    ## Run 113 stress 0.09022791 
    ## Run 114 stress 0.09252214 
    ## Run 115 stress 0.08985604 
    ## Run 116 stress 0.08985606 
    ## Run 117 stress 0.1077586 
    ## Run 118 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735035  max resid 0.05641944 
    ## Run 119 stress 0.106642 
    ## Run 120 stress 0.07936951 
    ## ... Procrustes: rmse 0.009040828  max resid 0.03243815 
    ## Run 121 stress 0.09388331 
    ## Run 122 stress 0.0917278 
    ## Run 123 stress 0.07936951 
    ## ... Procrustes: rmse 0.009038605  max resid 0.03242979 
    ## Run 124 stress 0.09172781 
    ## Run 125 stress 0.08985606 
    ## Run 126 stress 0.07927534 
    ## ... Procrustes: rmse 3.697304e-06  max resid 7.709698e-06 
    ## ... Similar to previous best
    ## Run 127 stress 0.07927535 
    ## ... Procrustes: rmse 3.444335e-05  max resid 9.513076e-05 
    ## ... Similar to previous best
    ## Run 128 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735802  max resid 0.05649823 
    ## Run 129 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735254  max resid 0.05641764 
    ## Run 130 stress 0.07927536 
    ## ... Procrustes: rmse 5.836563e-05  max resid 0.0001381584 
    ## ... Similar to previous best
    ## Run 131 stress 0.07936951 
    ## ... Procrustes: rmse 0.009037858  max resid 0.03241963 
    ## Run 132 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239239  max resid 0.04980373 
    ## Run 133 stress 0.07927534 
    ## ... Procrustes: rmse 7.57263e-06  max resid 1.610346e-05 
    ## ... Similar to previous best
    ## Run 134 stress 0.08985605 
    ## Run 135 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735267  max resid 0.05641354 
    ## Run 136 stress 0.09565372 
    ## Run 137 stress 0.07951445 
    ## ... Procrustes: rmse 0.01731515  max resid 0.0561722 
    ## Run 138 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043595  max resid 0.0324589 
    ## Run 139 stress 0.09252224 
    ## Run 140 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735673  max resid 0.05644477 
    ## Run 141 stress 0.09247104 
    ## Run 142 stress 0.07927534 
    ## ... Procrustes: rmse 3.789608e-06  max resid 1.073863e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.09544369 
    ## Run 144 stress 0.1033494 
    ## Run 145 stress 0.07936956 
    ## ... Procrustes: rmse 0.009022852  max resid 0.03232656 
    ## Run 146 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237714  max resid 0.049707 
    ## Run 147 stress 0.07969267 
    ## ... Procrustes: rmse 0.0123614  max resid 0.04952543 
    ## Run 148 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734472  max resid 0.05637311 
    ## Run 149 stress 0.0793695 
    ## ... Procrustes: rmse 0.009051065  max resid 0.03250201 
    ## Run 150 stress 0.07936949 
    ## ... Procrustes: rmse 0.009049561  max resid 0.03249488 
    ## Run 151 stress 0.0898561 
    ## Run 152 stress 0.09344086 
    ## Run 153 stress 0.07951443 
    ## ... Procrustes: rmse 0.01736867  max resid 0.05653616 
    ## Run 154 stress 0.07951442 
    ## ... Procrustes: rmse 0.01736006  max resid 0.05650931 
    ## Run 155 stress 0.09252217 
    ## Run 156 stress 0.09106043 
    ## Run 157 stress 0.08985608 
    ## Run 158 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045378  max resid 0.03246612 
    ## Run 159 stress 0.07927535 
    ## ... Procrustes: rmse 3.856309e-05  max resid 0.0001057389 
    ## ... Similar to previous best
    ## Run 160 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735646  max resid 0.05645372 
    ## Run 161 stress 0.09072897 
    ## Run 162 stress 0.09072894 
    ## Run 163 stress 0.07927536 
    ## ... Procrustes: rmse 4.811663e-05  max resid 0.0001341772 
    ## ... Similar to previous best
    ## Run 164 stress 0.09296255 
    ## Run 165 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 2.541211e-06  max resid 5.271354e-06 
    ## ... Similar to previous best
    ## Run 166 stress 0.09290689 
    ## Run 167 stress 0.09043762 
    ## Run 168 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039822  max resid 0.03243169 
    ## Run 169 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238937  max resid 0.04978298 
    ## Run 170 stress 0.08985603 
    ## Run 171 stress 0.07927534 
    ## ... Procrustes: rmse 6.708607e-06  max resid 1.491005e-05 
    ## ... Similar to previous best
    ## Run 172 stress 0.07969263 
    ## ... Procrustes: rmse 0.01241203  max resid 0.04992456 
    ## Run 173 stress 0.07969263 
    ## ... Procrustes: rmse 0.01241149  max resid 0.04992495 
    ## Run 174 stress 0.1039005 
    ## Run 175 stress 0.09072907 
    ## Run 176 stress 0.07951444 
    ## ... Procrustes: rmse 0.01733772  max resid 0.05631539 
    ## Run 177 stress 0.08985603 
    ## Run 178 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239531  max resid 0.04980311 
    ## Run 179 stress 0.08985608 
    ## Run 180 stress 0.09106039 
    ## Run 181 stress 0.07927534 
    ## ... Procrustes: rmse 8.989007e-06  max resid 2.098893e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.09485705 
    ## Run 183 stress 0.09022785 
    ## Run 184 stress 0.07927534 
    ## ... Procrustes: rmse 5.601321e-06  max resid 1.462056e-05 
    ## ... Similar to previous best
    ## Run 185 stress 0.08985614 
    ## Run 186 stress 0.08985603 
    ## Run 187 stress 0.07936955 
    ## ... Procrustes: rmse 0.009026022  max resid 0.03234544 
    ## Run 188 stress 0.09382085 
    ## Run 189 stress 0.09252216 
    ## Run 190 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904374  max resid 0.03245872 
    ## Run 191 stress 0.07927535 
    ## ... Procrustes: rmse 2.549988e-05  max resid 6.716661e-05 
    ## ... Similar to previous best
    ## Run 192 stress 0.08985603 
    ## Run 193 stress 0.07927536 
    ## ... Procrustes: rmse 1.58543e-05  max resid 3.063909e-05 
    ## ... Similar to previous best
    ## Run 194 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047379  max resid 0.03247808 
    ## Run 195 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041259  max resid 0.03245348 
    ## Run 196 stress 0.08985603 
    ## Run 197 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734337  max resid 0.05635687 
    ## Run 198 stress 0.07969261 
    ## ... Procrustes: rmse 0.01239238  max resid 0.04980483 
    ## Run 199 stress 0.09565383 
    ## Run 200 stress 0.07951445 
    ## ... Procrustes: rmse 0.01732199  max resid 0.05620133 
    ## Run 201 stress 0.07927534 
    ## ... Procrustes: rmse 1.489283e-06  max resid 3.649536e-06 
    ## ... Similar to previous best
    ## Run 202 stress 0.07936951 
    ## ... Procrustes: rmse 0.009039589  max resid 0.03241793 
    ## Run 203 stress 0.1083561 
    ## Run 204 stress 0.08985605 
    ## Run 205 stress 0.07927537 
    ## ... Procrustes: rmse 6.055612e-05  max resid 0.0001492503 
    ## ... Similar to previous best
    ## Run 206 stress 0.09043769 
    ## Run 207 stress 0.1033495 
    ## Run 208 stress 0.09043785 
    ## Run 209 stress 0.09072909 
    ## Run 210 stress 0.07927534 
    ## ... Procrustes: rmse 1.711081e-05  max resid 4.652345e-05 
    ## ... Similar to previous best
    ## Run 211 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735082  max resid 0.05643686 
    ## Run 212 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735341  max resid 0.05643813 
    ## Run 213 stress 0.07951444 
    ## ... Procrustes: rmse 0.01735949  max resid 0.05649277 
    ## Run 214 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734669  max resid 0.05638043 
    ## Run 215 stress 0.09106059 
    ## Run 216 stress 0.07951443 
    ## ... Procrustes: rmse 0.01733877  max resid 0.05635175 
    ## Run 217 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735926  max resid 0.05649065 
    ## Run 218 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237938  max resid 0.04972628 
    ## Run 219 stress 0.08985605 
    ## Run 220 stress 0.07927534 
    ## ... Procrustes: rmse 7.416089e-06  max resid 2.120495e-05 
    ## ... Similar to previous best
    ## Run 221 stress 0.09356315 
    ## Run 222 stress 0.1039003 
    ## Run 223 stress 0.09072891 
    ## Run 224 stress 0.07969259 
    ## ... Procrustes: rmse 0.01235976  max resid 0.04962809 
    ## Run 225 stress 0.09344097 
    ## Run 226 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735497  max resid 0.05646267 
    ## Run 227 stress 0.09252214 
    ## Run 228 stress 0.09252214 
    ## Run 229 stress 0.07951445 
    ## ... Procrustes: rmse 0.01737815  max resid 0.0565827 
    ## Run 230 stress 0.09388339 
    ## Run 231 stress 0.09043759 
    ## Run 232 stress 0.0792754 
    ## ... Procrustes: rmse 9.39268e-05  max resid 0.0002676755 
    ## ... Similar to previous best
    ## Run 233 stress 0.09252214 
    ## Run 234 stress 0.09382075 
    ## Run 235 stress 0.07969266 
    ## ... Procrustes: rmse 0.01237202  max resid 0.04958161 
    ## Run 236 stress 0.07936957 
    ## ... Procrustes: rmse 0.00905021  max resid 0.03249346 
    ## Run 237 stress 0.0910604 
    ## Run 238 stress 0.08985614 
    ## Run 239 stress 0.07927538 
    ## ... Procrustes: rmse 5.943807e-05  max resid 0.0001698269 
    ## ... Similar to previous best
    ## Run 240 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733661  max resid 0.05631058 
    ## Run 241 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237469  max resid 0.04968069 
    ## Run 242 stress 0.09072893 
    ## Run 243 stress 0.07927534 
    ## ... Procrustes: rmse 1.115614e-05  max resid 3.12013e-05 
    ## ... Similar to previous best
    ## Run 244 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735116  max resid 0.05639371 
    ## Run 245 stress 0.09106032 
    ## Run 246 stress 0.07927536 
    ## ... Procrustes: rmse 4.679005e-05  max resid 0.0001337924 
    ## ... Similar to previous best
    ## Run 247 stress 0.07951442 
    ## ... Procrustes: rmse 0.01736278  max resid 0.05651416 
    ## Run 248 stress 0.09106054 
    ## Run 249 stress 0.07927536 
    ## ... Procrustes: rmse 4.699011e-05  max resid 0.0001186193 
    ## ... Similar to previous best
    ## Run 250 stress 0.09043756 
    ## Run 251 stress 0.09022776 
    ## Run 252 stress 0.09388338 
    ## Run 253 stress 0.07927534 
    ## ... Procrustes: rmse 6.166506e-06  max resid 1.508197e-05 
    ## ... Similar to previous best
    ## Run 254 stress 0.07927535 
    ## ... Procrustes: rmse 3.804133e-05  max resid 0.0001025869 
    ## ... Similar to previous best
    ## Run 255 stress 0.07936951 
    ## ... Procrustes: rmse 0.00903484  max resid 0.03240024 
    ## Run 256 stress 0.07969259 
    ## ... Procrustes: rmse 0.0123811  max resid 0.04974124 
    ## Run 257 stress 0.09247104 
    ## Run 258 stress 0.09072899 
    ## Run 259 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735059  max resid 0.05640877 
    ## Run 260 stress 0.09072892 
    ## Run 261 stress 0.07936951 
    ## ... Procrustes: rmse 0.009034175  max resid 0.03239404 
    ## Run 262 stress 0.09296263 
    ## Run 263 stress 0.09172787 
    ## Run 264 stress 0.09072892 
    ## Run 265 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237652  max resid 0.04970498 
    ## Run 266 stress 0.09022778 
    ## Run 267 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238571  max resid 0.04976313 
    ## Run 268 stress 0.0929071 
    ## Run 269 stress 0.0910605 
    ## Run 270 stress 0.08985607 
    ## Run 271 stress 0.09043754 
    ## Run 272 stress 0.08985606 
    ## Run 273 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734121  max resid 0.05634448 
    ## Run 274 stress 0.09290708 
    ## Run 275 stress 0.08985603 
    ## Run 276 stress 0.0793695 
    ## ... Procrustes: rmse 0.009053564  max resid 0.03250699 
    ## Run 277 stress 0.07927535 
    ## ... Procrustes: rmse 1.220222e-05  max resid 3.491872e-05 
    ## ... Similar to previous best
    ## Run 278 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237684  max resid 0.0497158 
    ## Run 279 stress 0.09290709 
    ## Run 280 stress 0.09106044 
    ## Run 281 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735483  max resid 0.05646908 
    ## Run 282 stress 0.08985606 
    ## Run 283 stress 0.09172786 
    ## Run 284 stress 0.08985616 
    ## Run 285 stress 0.0917279 
    ## Run 286 stress 0.09344097 
    ## Run 287 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044869  max resid 0.03246365 
    ## Run 288 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240135  max resid 0.04986334 
    ## Run 289 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237406  max resid 0.0496812 
    ## Run 290 stress 0.07927538 
    ## ... Procrustes: rmse 6.790586e-05  max resid 0.0001692395 
    ## ... Similar to previous best
    ## Run 291 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238284  max resid 0.0497725 
    ## Run 292 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237138  max resid 0.04964114 
    ## Run 293 stress 0.07951443 
    ## ... Procrustes: rmse 0.01736664  max resid 0.05655885 
    ## Run 294 stress 0.08985603 
    ## Run 295 stress 0.07927534 
    ## ... Procrustes: rmse 1.374131e-05  max resid 3.866237e-05 
    ## ... Similar to previous best
    ## Run 296 stress 0.09022786 
    ## Run 297 stress 0.1039003 
    ## Run 298 stress 0.09172781 
    ## Run 299 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735275  max resid 0.05647367 
    ## Run 300 stress 0.09296237 
    ## Run 301 stress 0.09565384 
    ## Run 302 stress 0.09290702 
    ## Run 303 stress 0.07969266 
    ## ... Procrustes: rmse 0.01236187  max resid 0.04953435 
    ## Run 304 stress 0.07927534 
    ## ... Procrustes: rmse 4.482028e-06  max resid 1.337804e-05 
    ## ... Similar to previous best
    ## Run 305 stress 0.08985603 
    ## Run 306 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734885  max resid 0.05640738 
    ## Run 307 stress 0.09296244 
    ## Run 308 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734  max resid 0.05633422 
    ## Run 309 stress 0.09474042 
    ## Run 310 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735168  max resid 0.05640262 
    ## Run 311 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173521  max resid 0.05644098 
    ## Run 312 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046529  max resid 0.03247729 
    ## Run 313 stress 0.08985604 
    ## Run 314 stress 0.07969272 
    ## ... Procrustes: rmse 0.01236883  max resid 0.04953102 
    ## Run 315 stress 0.09106038 
    ## Run 316 stress 0.09043744 
    ## Run 317 stress 0.07969261 
    ## ... Procrustes: rmse 0.01238192  max resid 0.04973599 
    ## Run 318 stress 0.0898561 
    ## Run 319 stress 0.07927534 
    ## ... Procrustes: rmse 2.192995e-05  max resid 6.112545e-05 
    ## ... Similar to previous best
    ## Run 320 stress 0.08985606 
    ## Run 321 stress 0.0796926 
    ## ... Procrustes: rmse 0.01236466  max resid 0.04963355 
    ## Run 322 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732778  max resid 0.05623416 
    ## Run 323 stress 0.09344092 
    ## Run 324 stress 0.0956538 
    ## Run 325 stress 0.09072898 
    ## Run 326 stress 0.1097218 
    ## Run 327 stress 0.07927534 
    ## ... Procrustes: rmse 8.696458e-06  max resid 2.197648e-05 
    ## ... Similar to previous best
    ## Run 328 stress 0.09485757 
    ## Run 329 stress 0.07951445 
    ## ... Procrustes: rmse 0.01737178  max resid 0.05656213 
    ## Run 330 stress 0.07927535 
    ## ... Procrustes: rmse 3.856655e-05  max resid 0.0001030683 
    ## ... Similar to previous best
    ## Run 331 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045671  max resid 0.0324754 
    ## Run 332 stress 0.07951442 
    ## ... Procrustes: rmse 0.01738015  max resid 0.0565727 
    ## Run 333 stress 0.09106027 
    ## Run 334 stress 0.09290712 
    ## Run 335 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045882  max resid 0.03247374 
    ## Run 336 stress 0.09382092 
    ## Run 337 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734862  max resid 0.0563908 
    ## Run 338 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735572  max resid 0.05644782 
    ## Run 339 stress 0.09043761 
    ## Run 340 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734533  max resid 0.05639344 
    ## Run 341 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736472  max resid 0.05649229 
    ## Run 342 stress 0.1033494 
    ## Run 343 stress 0.07951442 
    ## ... Procrustes: rmse 0.0173562  max resid 0.0564785 
    ## Run 344 stress 0.09022787 
    ## Run 345 stress 0.07969266 
    ## ... Procrustes: rmse 0.01236804  max resid 0.04956814 
    ## Run 346 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045859  max resid 0.03247085 
    ## Run 347 stress 0.0796926 
    ## ... Procrustes: rmse 0.01239876  max resid 0.04982813 
    ## Run 348 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045446  max resid 0.03246913 
    ## Run 349 stress 0.1079617 
    ## Run 350 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237116  max resid 0.04965763 
    ## Run 351 stress 0.09072899 
    ## Run 352 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734045  max resid 0.05635257 
    ## Run 353 stress 0.09022786 
    ## Run 354 stress 0.08985603 
    ## Run 355 stress 0.1046273 
    ## Run 356 stress 0.07927534 
    ## ... Procrustes: rmse 4.982593e-06  max resid 1.071151e-05 
    ## ... Similar to previous best
    ## Run 357 stress 0.0929622 
    ## Run 358 stress 0.08985606 
    ## Run 359 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735186  max resid 0.05643382 
    ## Run 360 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043441  max resid 0.03246028 
    ## Run 361 stress 0.07951444 
    ## ... Procrustes: rmse 0.01735513  max resid 0.05652386 
    ## Run 362 stress 0.09106029 
    ## Run 363 stress 0.09485682 
    ## Run 364 stress 0.0793695 
    ## ... Procrustes: rmse 0.009040413  max resid 0.03242547 
    ## Run 365 stress 0.09072899 
    ## Run 366 stress 0.07936952 
    ## ... Procrustes: rmse 0.009050966  max resid 0.03249736 
    ## Run 367 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734898  max resid 0.05638554 
    ## Run 368 stress 0.08985611 
    ## Run 369 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735466  max resid 0.05646819 
    ## Run 370 stress 0.0793695 
    ## ... Procrustes: rmse 0.009051437  max resid 0.03249954 
    ## Run 371 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904526  max resid 0.03246933 
    ## Run 372 stress 0.08985606 
    ## Run 373 stress 0.07927534 
    ## ... Procrustes: rmse 1.338125e-05  max resid 3.474264e-05 
    ## ... Similar to previous best
    ## Run 374 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733818  max resid 0.05632233 
    ## Run 375 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240318  max resid 0.04987096 
    ## Run 376 stress 0.09106028 
    ## Run 377 stress 0.07969262 
    ## ... Procrustes: rmse 0.0124053  max resid 0.04989157 
    ## Run 378 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238757  max resid 0.04978987 
    ## Run 379 stress 0.07927534 
    ## ... Procrustes: rmse 1.200002e-05  max resid 3.301574e-05 
    ## ... Similar to previous best
    ## Run 380 stress 0.07936953 
    ## ... Procrustes: rmse 0.009029121  max resid 0.0323644 
    ## Run 381 stress 0.09043744 
    ## Run 382 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237456  max resid 0.04969443 
    ## Run 383 stress 0.09388343 
    ## Run 384 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735326  max resid 0.05642064 
    ## Run 385 stress 0.09290698 
    ## Run 386 stress 0.07927535 
    ## ... Procrustes: rmse 1.034938e-05  max resid 2.975401e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734751  max resid 0.0563892 
    ## Run 388 stress 0.108875 
    ## Run 389 stress 0.09388326 
    ## Run 390 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043843  max resid 0.03245617 
    ## Run 391 stress 0.08985609 
    ## Run 392 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735914  max resid 0.05646027 
    ## Run 393 stress 0.07927537 
    ## ... Procrustes: rmse 6.428297e-05  max resid 0.0001749403 
    ## ... Similar to previous best
    ## Run 394 stress 0.09043745 
    ## Run 395 stress 0.07927535 
    ## ... Procrustes: rmse 3.006077e-05  max resid 7.806415e-05 
    ## ... Similar to previous best
    ## Run 396 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733771  max resid 0.05632008 
    ## Run 397 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044327  max resid 0.0324614 
    ## Run 398 stress 0.105848 
    ## Run 399 stress 0.09388342 
    ## Run 400 stress 0.08985605 
    ## Run 401 stress 0.1033494 
    ## Run 402 stress 0.07927534 
    ## ... Procrustes: rmse 8.66705e-06  max resid 2.153059e-05 
    ## ... Similar to previous best
    ## Run 403 stress 0.09388334 
    ## Run 404 stress 0.07927534 
    ## ... Procrustes: rmse 8.256824e-06  max resid 2.292055e-05 
    ## ... Similar to previous best
    ## Run 405 stress 0.0793695 
    ## ... Procrustes: rmse 0.009050599  max resid 0.03249599 
    ## Run 406 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732232  max resid 0.05626874 
    ## Run 407 stress 0.08985605 
    ## Run 408 stress 0.0898562 
    ## Run 409 stress 0.07951443 
    ## ... Procrustes: rmse 0.01736638  max resid 0.05650223 
    ## Run 410 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046391  max resid 0.03247929 
    ## Run 411 stress 0.09474051 
    ## Run 412 stress 0.0924711 
    ## Run 413 stress 0.07927535 
    ## ... Procrustes: rmse 3.175958e-05  max resid 7.771287e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.07927535 
    ## ... Procrustes: rmse 1.469742e-05  max resid 4.110869e-05 
    ## ... Similar to previous best
    ## Run 415 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238259  max resid 0.04970064 
    ## Run 416 stress 0.08985603 
    ## Run 417 stress 0.09022776 
    ## Run 418 stress 0.106491 
    ## Run 419 stress 0.09072897 
    ## Run 420 stress 0.08985615 
    ## Run 421 stress 0.09072908 
    ## Run 422 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735324  max resid 0.05640827 
    ## Run 423 stress 0.08985619 
    ## Run 424 stress 0.09474044 
    ## Run 425 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735116  max resid 0.05640679 
    ## Run 426 stress 0.09382055 
    ## Run 427 stress 0.09565389 
    ## Run 428 stress 0.07951444 
    ## ... Procrustes: rmse 0.01737967  max resid 0.05660835 
    ## Run 429 stress 0.09043739 
    ## Run 430 stress 0.07927534 
    ## ... Procrustes: rmse 2.367333e-06  max resid 7.75767e-06 
    ## ... Similar to previous best
    ## Run 431 stress 0.09252228 
    ## Run 432 stress 0.07936959 
    ## ... Procrustes: rmse 0.009046406  max resid 0.0324695 
    ## Run 433 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735597  max resid 0.05647381 
    ## Run 434 stress 0.0793695 
    ## ... Procrustes: rmse 0.009051655  max resid 0.03250831 
    ## Run 435 stress 0.1079617 
    ## Run 436 stress 0.09043754 
    ## Run 437 stress 0.07927537 
    ## ... Procrustes: rmse 6.498051e-05  max resid 0.0001828114 
    ## ... Similar to previous best
    ## Run 438 stress 0.08985624 
    ## Run 439 stress 0.07927535 
    ## ... Procrustes: rmse 2.54147e-05  max resid 6.112817e-05 
    ## ... Similar to previous best
    ## Run 440 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237455  max resid 0.04970112 
    ## Run 441 stress 0.08985606 
    ## Run 442 stress 0.07927534 
    ## ... Procrustes: rmse 3.878394e-06  max resid 1.15445e-05 
    ## ... Similar to previous best
    ## Run 443 stress 0.08985613 
    ## Run 444 stress 0.1109895 
    ## Run 445 stress 0.07936951 
    ## ... Procrustes: rmse 0.009043846  max resid 0.03245767 
    ## Run 446 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237272  max resid 0.04966458 
    ## Run 447 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173574  max resid 0.05645072 
    ## Run 448 stress 0.07951442 
    ## ... Procrustes: rmse 0.01734933  max resid 0.05645361 
    ## Run 449 stress 0.09106027 
    ## Run 450 stress 0.07969263 
    ## ... Procrustes: rmse 0.01240854  max resid 0.04990975 
    ## Run 451 stress 0.094857 
    ## Run 452 stress 0.09022787 
    ## Run 453 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735244  max resid 0.05641574 
    ## Run 454 stress 0.1097216 
    ## Run 455 stress 0.07927534 
    ## ... Procrustes: rmse 1.395221e-05  max resid 3.624963e-05 
    ## ... Similar to previous best
    ## Run 456 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173562  max resid 0.05646762 
    ## Run 457 stress 0.09043736 
    ## Run 458 stress 0.07927534 
    ## ... Procrustes: rmse 1.75346e-05  max resid 4.866029e-05 
    ## ... Similar to previous best
    ## Run 459 stress 0.07927535 
    ## ... Procrustes: rmse 4.336441e-05  max resid 0.0001213754 
    ## ... Similar to previous best
    ## Run 460 stress 0.08985615 
    ## Run 461 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733064  max resid 0.05626539 
    ## Run 462 stress 0.07927535 
    ## ... Procrustes: rmse 1.657268e-05  max resid 3.646347e-05 
    ## ... Similar to previous best
    ## Run 463 stress 0.09247102 
    ## Run 464 stress 0.0793695 
    ## ... Procrustes: rmse 0.009051603  max resid 0.03250228 
    ## Run 465 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736644  max resid 0.05650061 
    ## Run 466 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734863  max resid 0.05639424 
    ## Run 467 stress 0.09072905 
    ## Run 468 stress 0.07936949 
    ## ... Procrustes: rmse 0.009044284  max resid 0.03246051 
    ## Run 469 stress 0.09252219 
    ## Run 470 stress 0.08985606 
    ## Run 471 stress 0.07927538 
    ## ... Procrustes: rmse 7.123139e-05  max resid 0.0001625134 
    ## ... Similar to previous best
    ## Run 472 stress 0.09106032 
    ## Run 473 stress 0.08985603 
    ## Run 474 stress 0.08985606 
    ## Run 475 stress 0.07936956 
    ## ... Procrustes: rmse 0.009022172  max resid 0.03232075 
    ## Run 476 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046695  max resid 0.03247642 
    ## Run 477 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734929  max resid 0.05637504 
    ## Run 478 stress 0.07927535 
    ## ... Procrustes: rmse 3.357358e-05  max resid 9.375331e-05 
    ## ... Similar to previous best
    ## Run 479 stress 0.07969262 
    ## ... Procrustes: rmse 0.01237055  max resid 0.04960748 
    ## Run 480 stress 0.09043737 
    ## Run 481 stress 0.09485748 
    ## Run 482 stress 0.07936949 
    ## ... Procrustes: rmse 0.009043707  max resid 0.03245829 
    ## Run 483 stress 0.07951444 
    ## ... Procrustes: rmse 0.01735563  max resid 0.05652557 
    ## Run 484 stress 0.07927535 
    ## ... Procrustes: rmse 3.904102e-05  max resid 0.0001091105 
    ## ... Similar to previous best
    ## Run 485 stress 0.08985608 
    ## Run 486 stress 0.07936951 
    ## ... Procrustes: rmse 0.009033948  max resid 0.03239432 
    ## Run 487 stress 0.0902279 
    ## Run 488 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735068  max resid 0.05641419 
    ## Run 489 stress 0.09072891 
    ## Run 490 stress 0.09106034 
    ## Run 491 stress 0.07951446 
    ## ... Procrustes: rmse 0.01739447  max resid 0.05665694 
    ## Run 492 stress 0.09043771 
    ## Run 493 stress 0.0793695 
    ## ... Procrustes: rmse 0.009038937  max resid 0.0324256 
    ## Run 494 stress 0.1058826 
    ## Run 495 stress 0.1138774 
    ## Run 496 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734515  max resid 0.05636203 
    ## Run 497 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237796  max resid 0.04969901 
    ## Run 498 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734389  max resid 0.05635766 
    ## Run 499 stress 0.07927534 
    ## ... Procrustes: rmse 7.203148e-06  max resid 2.106546e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.08985619 
    ## *** Best solution repeated 46 times

``` r
round(PD_beta_NMDS$stress, digits = 2)
```

    ## [1] 0.08

``` r
PD_beta_rep_NMDS <- metaMDS(PD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.3173923 
    ## Run 1 stress 0.3224113 
    ## Run 2 stress 0.3270627 
    ## Run 3 stress 0.3318741 
    ## Run 4 stress 0.3362625 
    ## Run 5 stress 0.3227267 
    ## Run 6 stress 0.3273855 
    ## Run 7 stress 0.3273489 
    ## Run 8 stress 0.321514 
    ## Run 9 stress 0.331593 
    ## Run 10 stress 0.3331967 
    ## Run 11 stress 0.3423502 
    ## Run 12 stress 0.3174232 
    ## ... Procrustes: rmse 0.17039  max resid 0.2611677 
    ## Run 13 stress 0.3115791 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1694969  max resid 0.3670726 
    ## Run 14 stress 0.3268008 
    ## Run 15 stress 0.3291027 
    ## Run 16 stress 0.3170161 
    ## Run 17 stress 0.3265977 
    ## Run 18 stress 0.3192591 
    ## Run 19 stress 0.3293816 
    ## Run 20 stress 0.320421 
    ## Run 21 stress 0.3270724 
    ## Run 22 stress 0.3316075 
    ## Run 23 stress 0.3260771 
    ## Run 24 stress 0.3256879 
    ## Run 25 stress 0.3347457 
    ## Run 26 stress 0.3332456 
    ## Run 27 stress 0.3171469 
    ## Run 28 stress 0.3207882 
    ## Run 29 stress 0.3362098 
    ## Run 30 stress 0.3309064 
    ## Run 31 stress 0.3311809 
    ## Run 32 stress 0.3807936 
    ## Run 33 stress 0.3189282 
    ## Run 34 stress 0.3234334 
    ## Run 35 stress 0.3168196 
    ## Run 36 stress 0.3269632 
    ## Run 37 stress 0.3321477 
    ## Run 38 stress 0.3285483 
    ## Run 39 stress 0.31652 
    ## Run 40 stress 0.3313573 
    ## Run 41 stress 0.3145287 
    ## Run 42 stress 0.32193 
    ## Run 43 stress 0.3294319 
    ## Run 44 stress 0.3273077 
    ## Run 45 stress 0.3152818 
    ## Run 46 stress 0.3351073 
    ## Run 47 stress 0.3183876 
    ## Run 48 stress 0.3279048 
    ## Run 49 stress 0.3341867 
    ## Run 50 stress 0.3327315 
    ## Run 51 stress 0.3161354 
    ## Run 52 stress 0.3172839 
    ## Run 53 stress 0.3340467 
    ## Run 54 stress 0.3299675 
    ## Run 55 stress 0.3380998 
    ## Run 56 stress 0.3229498 
    ## Run 57 stress 0.3133285 
    ## Run 58 stress 0.3380046 
    ## Run 59 stress 0.328147 
    ## Run 60 stress 0.3189354 
    ## Run 61 stress 0.3270451 
    ## Run 62 stress 0.3430095 
    ## Run 63 stress 0.3330908 
    ## Run 64 stress 0.3330297 
    ## Run 65 stress 0.3299218 
    ## Run 66 stress 0.3305431 
    ## Run 67 stress 0.3150449 
    ## Run 68 stress 0.3252018 
    ## Run 69 stress 0.3158019 
    ## Run 70 stress 0.3298335 
    ## Run 71 stress 0.3353778 
    ## Run 72 stress 0.326942 
    ## Run 73 stress 0.3213834 
    ## Run 74 stress 0.3269363 
    ## Run 75 stress 0.3367753 
    ## Run 76 stress 0.3328821 
    ## Run 77 stress 0.3337466 
    ## Run 78 stress 0.3178324 
    ## Run 79 stress 0.3283289 
    ## Run 80 stress 0.3246975 
    ## Run 81 stress 0.3189671 
    ## Run 82 stress 0.3283416 
    ## Run 83 stress 0.3166002 
    ## Run 84 stress 0.3347242 
    ## Run 85 stress 0.3338018 
    ## Run 86 stress 0.3223408 
    ## Run 87 stress 0.3333845 
    ## Run 88 stress 0.3450533 
    ## Run 89 stress 0.3280074 
    ## Run 90 stress 0.325956 
    ## Run 91 stress 0.3297832 
    ## Run 92 stress 0.3148609 
    ## Run 93 stress 0.3269097 
    ## Run 94 stress 0.3324304 
    ## Run 95 stress 0.3269064 
    ## Run 96 stress 0.3181101 
    ## Run 97 stress 0.3174334 
    ## Run 98 stress 0.3773854 
    ## Run 99 stress 0.3317544 
    ## Run 100 stress 0.3239248 
    ## Run 101 stress 0.3299957 
    ## Run 102 stress 0.3379926 
    ## Run 103 stress 0.3335524 
    ## Run 104 stress 0.3219398 
    ## Run 105 stress 0.3357625 
    ## Run 106 stress 0.3276201 
    ## Run 107 stress 0.3393697 
    ## Run 108 stress 0.3382774 
    ## Run 109 stress 0.3347501 
    ## Run 110 stress 0.3368598 
    ## Run 111 stress 0.3237879 
    ## Run 112 stress 0.3330069 
    ## Run 113 stress 0.3246251 
    ## Run 114 stress 0.3423346 
    ## Run 115 stress 0.3173625 
    ## Run 116 stress 0.3310092 
    ## Run 117 stress 0.3325979 
    ## Run 118 stress 0.3342923 
    ## Run 119 stress 0.3360555 
    ## Run 120 stress 0.3179288 
    ## Run 121 stress 0.3281531 
    ## Run 122 stress 0.3398909 
    ## Run 123 stress 0.3241861 
    ## Run 124 stress 0.3326355 
    ## Run 125 stress 0.3129897 
    ## Run 126 stress 0.3302799 
    ## Run 127 stress 0.3269778 
    ## Run 128 stress 0.3260227 
    ## Run 129 stress 0.3152663 
    ## Run 130 stress 0.3310506 
    ## Run 131 stress 0.3152675 
    ## Run 132 stress 0.3327384 
    ## Run 133 stress 0.3297779 
    ## Run 134 stress 0.3238164 
    ## Run 135 stress 0.3171366 
    ## Run 136 stress 0.331381 
    ## Run 137 stress 0.3269813 
    ## Run 138 stress 0.3135526 
    ## Run 139 stress 0.332038 
    ## Run 140 stress 0.3263944 
    ## Run 141 stress 0.3321222 
    ## Run 142 stress 0.3297219 
    ## Run 143 stress 0.31364 
    ## Run 144 stress 0.3201695 
    ## Run 145 stress 0.3221556 
    ## Run 146 stress 0.3303591 
    ## Run 147 stress 0.3289448 
    ## Run 148 stress 0.3320731 
    ## Run 149 stress 0.3329263 
    ## Run 150 stress 0.3190501 
    ## Run 151 stress 0.3237971 
    ## Run 152 stress 0.3234538 
    ## Run 153 stress 0.312114 
    ## Run 154 stress 0.3403999 
    ## Run 155 stress 0.3416015 
    ## Run 156 stress 0.3237227 
    ## Run 157 stress 0.3171499 
    ## Run 158 stress 0.3207932 
    ## Run 159 stress 0.3246922 
    ## Run 160 stress 0.3302342 
    ## Run 161 stress 0.3347101 
    ## Run 162 stress 0.3126141 
    ## Run 163 stress 0.3270947 
    ## Run 164 stress 0.3329629 
    ## Run 165 stress 0.327207 
    ## Run 166 stress 0.3375525 
    ## Run 167 stress 0.3154074 
    ## Run 168 stress 0.3241816 
    ## Run 169 stress 0.3331669 
    ## Run 170 stress 0.3377613 
    ## Run 171 stress 0.3329732 
    ## Run 172 stress 0.3145721 
    ## Run 173 stress 0.333891 
    ## Run 174 stress 0.3303532 
    ## Run 175 stress 0.330498 
    ## Run 176 stress 0.3218518 
    ## Run 177 stress 0.3278267 
    ## Run 178 stress 0.3149295 
    ## Run 179 stress 0.3357241 
    ## Run 180 stress 0.3305324 
    ## Run 181 stress 0.3312979 
    ## Run 182 stress 0.3329201 
    ## Run 183 stress 0.3337223 
    ## Run 184 stress 0.3361927 
    ## Run 185 stress 0.3192929 
    ## Run 186 stress 0.3362634 
    ## Run 187 stress 0.3232054 
    ## Run 188 stress 0.3211222 
    ## Run 189 stress 0.3342659 
    ## Run 190 stress 0.3359537 
    ## Run 191 stress 0.3310531 
    ## Run 192 stress 0.3307298 
    ## Run 193 stress 0.3276612 
    ## Run 194 stress 0.3382113 
    ## Run 195 stress 0.3217314 
    ## Run 196 stress 0.321641 
    ## Run 197 stress 0.3179069 
    ## Run 198 stress 0.3294986 
    ## Run 199 stress 0.3132812 
    ## Run 200 stress 0.3344655 
    ## Run 201 stress 0.3280095 
    ## Run 202 stress 0.3328586 
    ## Run 203 stress 0.3253563 
    ## Run 204 stress 0.3142881 
    ## Run 205 stress 0.3205993 
    ## Run 206 stress 0.3183361 
    ## Run 207 stress 0.324513 
    ## Run 208 stress 0.315286 
    ## Run 209 stress 0.332534 
    ## Run 210 stress 0.316357 
    ## Run 211 stress 0.3184634 
    ## Run 212 stress 0.3149416 
    ## Run 213 stress 0.3165079 
    ## Run 214 stress 0.3170033 
    ## Run 215 stress 0.3210856 
    ## Run 216 stress 0.3327529 
    ## Run 217 stress 0.322662 
    ## Run 218 stress 0.3186871 
    ## Run 219 stress 0.317078 
    ## Run 220 stress 0.3161624 
    ## Run 221 stress 0.3199391 
    ## Run 222 stress 0.3222802 
    ## Run 223 stress 0.3336527 
    ## Run 224 stress 0.332642 
    ## Run 225 stress 0.3249779 
    ## Run 226 stress 0.3333827 
    ## Run 227 stress 0.3276856 
    ## Run 228 stress 0.3243195 
    ## Run 229 stress 0.3309055 
    ## Run 230 stress 0.3380628 
    ## Run 231 stress 0.3328971 
    ## Run 232 stress 0.3205105 
    ## Run 233 stress 0.3292426 
    ## Run 234 stress 0.3282515 
    ## Run 235 stress 0.3254299 
    ## Run 236 stress 0.3385376 
    ## Run 237 stress 0.3334466 
    ## Run 238 stress 0.3338742 
    ## Run 239 stress 0.3313906 
    ## Run 240 stress 0.3385612 
    ## Run 241 stress 0.3198444 
    ## Run 242 stress 0.3235555 
    ## Run 243 stress 0.3228281 
    ## Run 244 stress 0.3162659 
    ## Run 245 stress 0.3226315 
    ## Run 246 stress 0.3246502 
    ## Run 247 stress 0.319987 
    ## Run 248 stress 0.3322102 
    ## Run 249 stress 0.3337514 
    ## Run 250 stress 0.3301647 
    ## Run 251 stress 0.331054 
    ## Run 252 stress 0.3298983 
    ## Run 253 stress 0.3309545 
    ## Run 254 stress 0.3166476 
    ## Run 255 stress 0.334179 
    ## Run 256 stress 0.3355726 
    ## Run 257 stress 0.3198595 
    ## Run 258 stress 0.3338823 
    ## Run 259 stress 0.3343364 
    ## Run 260 stress 0.327241 
    ## Run 261 stress 0.3302521 
    ## Run 262 stress 0.3289139 
    ## Run 263 stress 0.3274346 
    ## Run 264 stress 0.3314699 
    ## Run 265 stress 0.3179431 
    ## Run 266 stress 0.3319381 
    ## Run 267 stress 0.3144882 
    ## Run 268 stress 0.3283997 
    ## Run 269 stress 0.320357 
    ## Run 270 stress 0.333358 
    ## Run 271 stress 0.3181546 
    ## Run 272 stress 0.3282438 
    ## Run 273 stress 0.333969 
    ## Run 274 stress 0.3356215 
    ## Run 275 stress 0.3314016 
    ## Run 276 stress 0.3189345 
    ## Run 277 stress 0.3374007 
    ## Run 278 stress 0.3356574 
    ## Run 279 stress 0.3291019 
    ## Run 280 stress 0.31986 
    ## Run 281 stress 0.3154972 
    ## Run 282 stress 0.311165 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1733616  max resid 0.300518 
    ## Run 283 stress 0.3152963 
    ## Run 284 stress 0.3385757 
    ## Run 285 stress 0.326155 
    ## Run 286 stress 0.3183714 
    ## Run 287 stress 0.3184768 
    ## Run 288 stress 0.3327251 
    ## Run 289 stress 0.3383209 
    ## Run 290 stress 0.3123206 
    ## Run 291 stress 0.3190922 
    ## Run 292 stress 0.3374124 
    ## Run 293 stress 0.3267759 
    ## Run 294 stress 0.3326873 
    ## Run 295 stress 0.3283714 
    ## Run 296 stress 0.3259399 
    ## Run 297 stress 0.3269965 
    ## Run 298 stress 0.3204639 
    ## Run 299 stress 0.339345 
    ## Run 300 stress 0.3232987 
    ## Run 301 stress 0.3328108 
    ## Run 302 stress 0.3335141 
    ## Run 303 stress 0.3243225 
    ## Run 304 stress 0.3197125 
    ## Run 305 stress 0.32413 
    ## Run 306 stress 0.3398004 
    ## Run 307 stress 0.316531 
    ## Run 308 stress 0.3284069 
    ## Run 309 stress 0.3217767 
    ## Run 310 stress 0.3140358 
    ## Run 311 stress 0.335558 
    ## Run 312 stress 0.3269583 
    ## Run 313 stress 0.3169854 
    ## Run 314 stress 0.337052 
    ## Run 315 stress 0.3146542 
    ## Run 316 stress 0.3387096 
    ## Run 317 stress 0.3324355 
    ## Run 318 stress 0.3245466 
    ## Run 319 stress 0.3321154 
    ## Run 320 stress 0.364006 
    ## Run 321 stress 0.3193684 
    ## Run 322 stress 0.3284913 
    ## Run 323 stress 0.334538 
    ## Run 324 stress 0.3406407 
    ## Run 325 stress 0.3183165 
    ## Run 326 stress 0.3273703 
    ## Run 327 stress 0.3217114 
    ## Run 328 stress 0.311841 
    ## Run 329 stress 0.3243234 
    ## Run 330 stress 0.3289447 
    ## Run 331 stress 0.3263899 
    ## Run 332 stress 0.3291557 
    ## Run 333 stress 0.3293888 
    ## Run 334 stress 0.3332008 
    ## Run 335 stress 0.3160033 
    ## Run 336 stress 0.3195683 
    ## Run 337 stress 0.3296277 
    ## Run 338 stress 0.319181 
    ## Run 339 stress 0.3257481 
    ## Run 340 stress 0.3206681 
    ## Run 341 stress 0.3373832 
    ## Run 342 stress 0.3239695 
    ## Run 343 stress 0.3356876 
    ## Run 344 stress 0.3339534 
    ## Run 345 stress 0.3322911 
    ## Run 346 stress 0.3214581 
    ## Run 347 stress 0.3312958 
    ## Run 348 stress 0.3298686 
    ## Run 349 stress 0.3295338 
    ## Run 350 stress 0.3773784 
    ## Run 351 stress 0.3249179 
    ## Run 352 stress 0.335809 
    ## Run 353 stress 0.3185015 
    ## Run 354 stress 0.3181809 
    ## Run 355 stress 0.3269839 
    ## Run 356 stress 0.3335187 
    ## Run 357 stress 0.3297266 
    ## Run 358 stress 0.3269044 
    ## Run 359 stress 0.3350086 
    ## Run 360 stress 0.3216461 
    ## Run 361 stress 0.3322688 
    ## Run 362 stress 0.3309835 
    ## Run 363 stress 0.3234837 
    ## Run 364 stress 0.3308252 
    ## Run 365 stress 0.3328359 
    ## Run 366 stress 0.3211273 
    ## Run 367 stress 0.3260929 
    ## Run 368 stress 0.3185861 
    ## Run 369 stress 0.3349288 
    ## Run 370 stress 0.3159188 
    ## Run 371 stress 0.3331084 
    ## Run 372 stress 0.3222582 
    ## Run 373 stress 0.3122743 
    ## Run 374 stress 0.3316514 
    ## Run 375 stress 0.3292251 
    ## Run 376 stress 0.3179447 
    ## Run 377 stress 0.3269461 
    ## Run 378 stress 0.3242886 
    ## Run 379 stress 0.3359235 
    ## Run 380 stress 0.3188416 
    ## Run 381 stress 0.3153885 
    ## Run 382 stress 0.3241327 
    ## Run 383 stress 0.3193533 
    ## Run 384 stress 0.3357097 
    ## Run 385 stress 0.3206922 
    ## Run 386 stress 0.3344469 
    ## Run 387 stress 0.3282106 
    ## Run 388 stress 0.3270823 
    ## Run 389 stress 0.3341621 
    ## Run 390 stress 0.3290088 
    ## Run 391 stress 0.3206336 
    ## Run 392 stress 0.3228108 
    ## Run 393 stress 0.3358849 
    ## Run 394 stress 0.3311793 
    ## Run 395 stress 0.3309655 
    ## Run 396 stress 0.3382925 
    ## Run 397 stress 0.330774 
    ## Run 398 stress 0.3158687 
    ## Run 399 stress 0.3298597 
    ## Run 400 stress 0.3333187 
    ## Run 401 stress 0.3168205 
    ## Run 402 stress 0.3383268 
    ## Run 403 stress 0.3342566 
    ## Run 404 stress 0.3344712 
    ## Run 405 stress 0.3223147 
    ## Run 406 stress 0.3210929 
    ## Run 407 stress 0.3237408 
    ## Run 408 stress 0.3171416 
    ## Run 409 stress 0.3178408 
    ## Run 410 stress 0.3321412 
    ## Run 411 stress 0.333407 
    ## Run 412 stress 0.332185 
    ## Run 413 stress 0.3405004 
    ## Run 414 stress 0.3206879 
    ## Run 415 stress 0.3408535 
    ## Run 416 stress 0.3196347 
    ## Run 417 stress 0.3216378 
    ## Run 418 stress 0.3216535 
    ## Run 419 stress 0.3334611 
    ## Run 420 stress 0.3255696 
    ## Run 421 stress 0.3344613 
    ## Run 422 stress 0.3305294 
    ## Run 423 stress 0.3232165 
    ## Run 424 stress 0.3334418 
    ## Run 425 stress 0.3159579 
    ## Run 426 stress 0.3345099 
    ## Run 427 stress 0.3386048 
    ## Run 428 stress 0.3299404 
    ## Run 429 stress 0.3345207 
    ## Run 430 stress 0.3160723 
    ## Run 431 stress 0.3337416 
    ## Run 432 stress 0.32004 
    ## Run 433 stress 0.3256209 
    ## Run 434 stress 0.3209409 
    ## Run 435 stress 0.3336975 
    ## Run 436 stress 0.3294868 
    ## Run 437 stress 0.3233572 
    ## Run 438 stress 0.3167214 
    ## Run 439 stress 0.3378555 
    ## Run 440 stress 0.3348793 
    ## Run 441 stress 0.3358684 
    ## Run 442 stress 0.3357419 
    ## Run 443 stress 0.3227542 
    ## Run 444 stress 0.3170364 
    ## Run 445 stress 0.3156859 
    ## Run 446 stress 0.3144129 
    ## Run 447 stress 0.3186663 
    ## Run 448 stress 0.3301857 
    ## Run 449 stress 0.3374079 
    ## Run 450 stress 0.3151913 
    ## Run 451 stress 0.3209163 
    ## Run 452 stress 0.3256515 
    ## Run 453 stress 0.319927 
    ## Run 454 stress 0.3169633 
    ## Run 455 stress 0.324192 
    ## Run 456 stress 0.3283091 
    ## Run 457 stress 0.3334767 
    ## Run 458 stress 0.3144921 
    ## Run 459 stress 0.3365061 
    ## Run 460 stress 0.326675 
    ## Run 461 stress 0.3136821 
    ## Run 462 stress 0.3226996 
    ## Run 463 stress 0.3309902 
    ## Run 464 stress 0.3264395 
    ## Run 465 stress 0.3350768 
    ## Run 466 stress 0.3232092 
    ## Run 467 stress 0.3183955 
    ## Run 468 stress 0.3369049 
    ## Run 469 stress 0.3201677 
    ## Run 470 stress 0.3369896 
    ## Run 471 stress 0.3180423 
    ## Run 472 stress 0.3237118 
    ## Run 473 stress 0.3238293 
    ## Run 474 stress 0.3346688 
    ## Run 475 stress 0.33854 
    ## Run 476 stress 0.3248679 
    ## Run 477 stress 0.3290895 
    ## Run 478 stress 0.3199568 
    ## Run 479 stress 0.3381336 
    ## Run 480 stress 0.3168592 
    ## Run 481 stress 0.317845 
    ## Run 482 stress 0.3202747 
    ## Run 483 stress 0.3398344 
    ## Run 484 stress 0.3262227 
    ## Run 485 stress 0.3342301 
    ## Run 486 stress 0.340145 
    ## Run 487 stress 0.325834 
    ## Run 488 stress 0.3246395 
    ## Run 489 stress 0.3147187 
    ## Run 490 stress 0.3308985 
    ## Run 491 stress 0.3257015 
    ## Run 492 stress 0.3226379 
    ## Run 493 stress 0.3292268 
    ## Run 494 stress 0.3148602 
    ## Run 495 stress 0.3283894 
    ## Run 496 stress 0.3158465 
    ## Run 497 stress 0.32604 
    ## Run 498 stress 0.3153034 
    ## Run 499 stress 0.3246507 
    ## Run 500 stress 0.3303115 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##    500: stress ratio > sratmax

``` r
PD_beta_ric_NMDS <- metaMDS(PD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.01115177 
    ## Run 1 stress 0.01115187 
    ## ... Procrustes: rmse 3.471152e-05  max resid 7.833082e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.01558312 
    ## Run 3 stress 0.01123235 
    ## ... Procrustes: rmse 0.007469264  max resid 0.01681186 
    ## Run 4 stress 0.0111879 
    ## ... Procrustes: rmse 0.005308882  max resid 0.0119554 
    ## Run 5 stress 0.01700767 
    ## Run 6 stress 0.01552945 
    ## Run 7 stress 0.01549946 
    ## Run 8 stress 0.01528201 
    ## Run 9 stress 0.0111517 
    ## ... New best solution
    ## ... Procrustes: rmse 2.950173e-05  max resid 6.657583e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.01549953 
    ## Run 11 stress 0.01115182 
    ## ... Procrustes: rmse 5.202926e-05  max resid 0.000117613 
    ## ... Similar to previous best
    ## Run 12 stress 0.01498594 
    ## Run 13 stress 0.0154551 
    ## Run 14 stress 0.01115186 
    ## ... Procrustes: rmse 6.834397e-05  max resid 0.0001546063 
    ## ... Similar to previous best
    ## Run 15 stress 0.01116188 
    ## ... Procrustes: rmse 0.003204739  max resid 0.007223758 
    ## ... Similar to previous best
    ## Run 16 stress 0.01526228 
    ## Run 17 stress 0.01515712 
    ## Run 18 stress 0.01558243 
    ## Run 19 stress 0.01521389 
    ## Run 20 stress 0.01527185 
    ## Run 21 stress 0.01547486 
    ## Run 22 stress 0.01574937 
    ## Run 23 stress 0.01649862 
    ## Run 24 stress 0.01556927 
    ## Run 25 stress 0.01498606 
    ## Run 26 stress 0.01515726 
    ## Run 27 stress 0.01529419 
    ## Run 28 stress 0.01556151 
    ## Run 29 stress 0.01536086 
    ## Run 30 stress 0.01556101 
    ## Run 31 stress 0.01498596 
    ## Run 32 stress 0.01498613 
    ## Run 33 stress 0.0157627 
    ## Run 34 stress 0.01553036 
    ## Run 35 stress 0.01546071 
    ## Run 36 stress 0.01562164 
    ## Run 37 stress 0.01115172 
    ## ... Procrustes: rmse 7.008021e-06  max resid 1.602276e-05 
    ## ... Similar to previous best
    ## Run 38 stress 0.01515705 
    ## Run 39 stress 0.01649867 
    ## Run 40 stress 0.01529073 
    ## Run 41 stress 0.01115184 
    ## ... Procrustes: rmse 5.76885e-05  max resid 0.0001302915 
    ## ... Similar to previous best
    ## Run 42 stress 0.01549935 
    ## Run 43 stress 0.01115195 
    ## ... Procrustes: rmse 9.998296e-05  max resid 0.0002258255 
    ## ... Similar to previous best
    ## Run 44 stress 0.01552929 
    ## Run 45 stress 0.01515744 
    ## Run 46 stress 0.0168542 
    ## Run 47 stress 0.01527929 
    ## Run 48 stress 0.01515671 
    ## Run 49 stress 0.01516388 
    ## Run 50 stress 0.01668081 
    ## Run 51 stress 0.01515693 
    ## Run 52 stress 0.01668205 
    ## Run 53 stress 0.01527577 
    ## Run 54 stress 0.01655699 
    ## Run 55 stress 0.01493173 
    ## Run 56 stress 0.01283353 
    ## Run 57 stress 0.01583467 
    ## Run 58 stress 0.01564567 
    ## Run 59 stress 0.01148852 
    ## ... Procrustes: rmse 0.01332874  max resid 0.02905266 
    ## Run 60 stress 0.01118036 
    ## ... Procrustes: rmse 0.004836385  max resid 0.01089729 
    ## Run 61 stress 0.01543904 
    ## Run 62 stress 0.01543179 
    ## Run 63 stress 0.01552976 
    ## Run 64 stress 0.01115194 
    ## ... Procrustes: rmse 9.878091e-05  max resid 0.0002232469 
    ## ... Similar to previous best
    ## Run 65 stress 0.01528212 
    ## Run 66 stress 0.01115176 
    ## ... Procrustes: rmse 0.001341155  max resid 0.003023628 
    ## ... Similar to previous best
    ## Run 67 stress 0.01499743 
    ## Run 68 stress 0.01583038 
    ## Run 69 stress 0.01499111 
    ## Run 70 stress 0.01582937 
    ## Run 71 stress 0.01115203 
    ## ... Procrustes: rmse 0.001445304  max resid 0.003258207 
    ## ... Similar to previous best
    ## Run 72 stress 0.0151567 
    ## Run 73 stress 0.01525807 
    ## Run 74 stress 0.01499274 
    ## Run 75 stress 0.01493159 
    ## Run 76 stress 0.01556271 
    ## Run 77 stress 0.01547206 
    ## Run 78 stress 0.0167181 
    ## Run 79 stress 0.01576366 
    ## Run 80 stress 0.01525788 
    ## Run 81 stress 0.01527565 
    ## Run 82 stress 0.01115178 
    ## ... Procrustes: rmse 3.585068e-05  max resid 8.095083e-05 
    ## ... Similar to previous best
    ## Run 83 stress 0.01638348 
    ## Run 84 stress 0.01673573 
    ## Run 85 stress 0.01528164 
    ## Run 86 stress 0.01581979 
    ## Run 87 stress 0.01663496 
    ## Run 88 stress 0.0155738 
    ## Run 89 stress 0.0149859 
    ## Run 90 stress 0.01558318 
    ## Run 91 stress 0.01583074 
    ## Run 92 stress 0.01498589 
    ## Run 93 stress 0.01648667 
    ## Run 94 stress 0.01546647 
    ## Run 95 stress 0.01555998 
    ## Run 96 stress 0.01157654 
    ## ... Procrustes: rmse 0.01665466  max resid 0.03682994 
    ## Run 97 stress 0.01499247 
    ## Run 98 stress 0.01558207 
    ## Run 99 stress 0.01115343 
    ## ... Procrustes: rmse 0.001870748  max resid 0.004217891 
    ## ... Similar to previous best
    ## Run 100 stress 0.01658193 
    ## Run 101 stress 0.01129312 
    ## ... Procrustes: rmse 0.00985043  max resid 0.0221553 
    ## Run 102 stress 0.01115173 
    ## ... Procrustes: rmse 0.001336063  max resid 0.003012437 
    ## ... Similar to previous best
    ## Run 103 stress 0.01248106 
    ## Run 104 stress 0.01115969 
    ## ... Procrustes: rmse 0.00293721  max resid 0.006620988 
    ## ... Similar to previous best
    ## Run 105 stress 0.01136018 
    ## ... Procrustes: rmse 0.01106236  max resid 0.02472366 
    ## Run 106 stress 0.0153378 
    ## Run 107 stress 0.01115196 
    ## ... Procrustes: rmse 0.001423547  max resid 0.003210171 
    ## ... Similar to previous best
    ## Run 108 stress 0.01583018 
    ## Run 109 stress 0.01115204 
    ## ... Procrustes: rmse 0.001453242  max resid 0.003277165 
    ## ... Similar to previous best
    ## Run 110 stress 0.01115174 
    ## ... Procrustes: rmse 0.001338206  max resid 0.003017157 
    ## ... Similar to previous best
    ## Run 111 stress 0.01119215 
    ## ... Procrustes: rmse 0.00559904  max resid 0.01261092 
    ## Run 112 stress 0.0154741 
    ## Run 113 stress 0.01499332 
    ## Run 114 stress 0.01641894 
    ## Run 115 stress 0.01545361 
    ## Run 116 stress 0.01499285 
    ## Run 117 stress 0.01556991 
    ## Run 118 stress 0.01545268 
    ## Run 119 stress 0.01659998 
    ## Run 120 stress 0.01115172 
    ## ... Procrustes: rmse 0.001323403  max resid 0.002983557 
    ## ... Similar to previous best
    ## Run 121 stress 0.01583007 
    ## Run 122 stress 0.01549954 
    ## Run 123 stress 0.01533816 
    ## Run 124 stress 0.01115168 
    ## ... New best solution
    ## ... Procrustes: rmse 1.224193e-05  max resid 2.55178e-05 
    ## ... Similar to previous best
    ## Run 125 stress 0.01545449 
    ## Run 126 stress 0.01150545 
    ## ... Procrustes: rmse 0.01381474  max resid 0.03001973 
    ## Run 127 stress 0.01123129 
    ## ... Procrustes: rmse 0.007555919  max resid 0.01700142 
    ## Run 128 stress 0.01115133 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001974267  max resid 0.0004454141 
    ## ... Similar to previous best
    ## Run 129 stress 0.0111517 
    ## ... Procrustes: rmse 0.000206837  max resid 0.0004667331 
    ## ... Similar to previous best
    ## Run 130 stress 0.01493188 
    ## Run 131 stress 0.01498601 
    ## Run 132 stress 0.01544876 
    ## Run 133 stress 0.01529113 
    ## Run 134 stress 0.01528158 
    ## Run 135 stress 0.01115536 
    ## ... Procrustes: rmse 0.002067692  max resid 0.004661677 
    ## ... Similar to previous best
    ## Run 136 stress 0.01550062 
    ## Run 137 stress 0.01136802 
    ## ... Procrustes: rmse 0.0109571  max resid 0.0244695 
    ## Run 138 stress 0.01115183 
    ## ... Procrustes: rmse 0.0002599272  max resid 0.0005870221 
    ## ... Similar to previous best
    ## Run 139 stress 0.01577209 
    ## Run 140 stress 0.01515732 
    ## Run 141 stress 0.01542735 
    ## Run 142 stress 0.01566991 
    ## Run 143 stress 0.01115166 
    ## ... Procrustes: rmse 0.001089186  max resid 0.00245543 
    ## ... Similar to previous best
    ## Run 144 stress 0.01556092 
    ## Run 145 stress 0.01115189 
    ## ... Procrustes: rmse 0.0002827863  max resid 0.0006384096 
    ## ... Similar to previous best
    ## Run 146 stress 0.01554727 
    ## Run 147 stress 0.01659995 
    ## Run 148 stress 0.01564595 
    ## Run 149 stress 0.01515716 
    ## Run 150 stress 0.01528207 
    ## Run 151 stress 0.01528178 
    ## Run 152 stress 0.01515202 
    ## Run 153 stress 0.01549927 
    ## Run 154 stress 0.01116913 
    ## ... Procrustes: rmse 0.003735732  max resid 0.008419394 
    ## ... Similar to previous best
    ## Run 155 stress 0.01542716 
    ## Run 156 stress 0.01154079 
    ## ... Procrustes: rmse 0.01467336  max resid 0.03175127 
    ## Run 157 stress 0.01549914 
    ## Run 158 stress 0.0151799 
    ## Run 159 stress 0.015157 
    ## Run 160 stress 0.01567245 
    ## Run 161 stress 0.0155696 
    ## Run 162 stress 0.01127408 
    ## ... Procrustes: rmse 0.008955529  max resid 0.02014174 
    ## Run 163 stress 0.0152818 
    ## Run 164 stress 0.01543874 
    ## Run 165 stress 0.01115201 
    ## ... Procrustes: rmse 0.0003269202  max resid 0.000738324 
    ## ... Similar to previous best
    ## Run 166 stress 0.01674199 
    ## Run 167 stress 0.01558321 
    ## Run 168 stress 0.01527586 
    ## Run 169 stress 0.01655603 
    ## Run 170 stress 0.01558617 
    ## Run 171 stress 0.01545708 
    ## Run 172 stress 0.0111952 
    ## ... Procrustes: rmse 0.005554718  max resid 0.01250952 
    ## Run 173 stress 0.01658214 
    ## Run 174 stress 0.01564595 
    ## Run 175 stress 0.01558331 
    ## Run 176 stress 0.0152817 
    ## Run 177 stress 0.01562184 
    ## Run 178 stress 0.01152921 
    ## ... Procrustes: rmse 0.01431467  max resid 0.03100112 
    ## Run 179 stress 0.01525953 
    ## Run 180 stress 0.01533784 
    ## Run 181 stress 0.01659958 
    ## Run 182 stress 0.01569098 
    ## Run 183 stress 0.01545811 
    ## Run 184 stress 0.01498596 
    ## Run 185 stress 0.0149919 
    ## Run 186 stress 0.01533729 
    ## Run 187 stress 0.015258 
    ## Run 188 stress 0.01515726 
    ## Run 189 stress 0.01562201 
    ## Run 190 stress 0.01514929 
    ## Run 191 stress 0.01564572 
    ## Run 192 stress 0.015579 
    ## Run 193 stress 0.01557867 
    ## Run 194 stress 0.01571559 
    ## Run 195 stress 0.01556964 
    ## Run 196 stress 0.0158196 
    ## Run 197 stress 0.01499197 
    ## Run 198 stress 0.01558358 
    ## Run 199 stress 0.01548131 
    ## Run 200 stress 0.01581904 
    ## Run 201 stress 0.01581392 
    ## Run 202 stress 0.01498578 
    ## Run 203 stress 0.01547387 
    ## Run 204 stress 0.01650625 
    ## Run 205 stress 0.01493172 
    ## Run 206 stress 0.0156215 
    ## Run 207 stress 0.01493149 
    ## Run 208 stress 0.01536149 
    ## Run 209 stress 0.01116132 
    ## ... Procrustes: rmse 0.002829166  max resid 0.006374331 
    ## ... Similar to previous best
    ## Run 210 stress 0.01115177 
    ## ... Procrustes: rmse 0.001143458  max resid 0.002578229 
    ## ... Similar to previous best
    ## Run 211 stress 0.01498598 
    ## Run 212 stress 0.0154815 
    ## Run 213 stress 0.01499249 
    ## Run 214 stress 0.01557897 
    ## Run 215 stress 0.01548002 
    ## Run 216 stress 0.01610254 
    ## Run 217 stress 0.01127658 
    ## ... Procrustes: rmse 0.009090685  max resid 0.0204148 
    ## Run 218 stress 0.01528199 
    ## Run 219 stress 0.01115175 
    ## ... Procrustes: rmse 0.0002310497  max resid 0.0005215106 
    ## ... Similar to previous best
    ## Run 220 stress 0.01115548 
    ## ... Procrustes: rmse 0.002085737  max resid 0.004701915 
    ## ... Similar to previous best
    ## Run 221 stress 0.01115241 
    ## ... Procrustes: rmse 0.001379337  max resid 0.003109974 
    ## ... Similar to previous best
    ## Run 222 stress 0.01555846 
    ## Run 223 stress 0.01499232 
    ## Run 224 stress 0.01560083 
    ## Run 225 stress 0.01582971 
    ## Run 226 stress 0.01498605 
    ## Run 227 stress 0.01556876 
    ## Run 228 stress 0.01527555 
    ## Run 229 stress 0.01576788 
    ## Run 230 stress 0.01649852 
    ## Run 231 stress 0.01515012 
    ## Run 232 stress 0.01552999 
    ## Run 233 stress 0.01538082 
    ## Run 234 stress 0.01525789 
    ## Run 235 stress 0.01538056 
    ## Run 236 stress 0.01700829 
    ## Run 237 stress 0.01557932 
    ## Run 238 stress 0.01522322 
    ## Run 239 stress 0.01688316 
    ## Run 240 stress 0.01557873 
    ## Run 241 stress 0.01525785 
    ## Run 242 stress 0.01528192 
    ## Run 243 stress 0.01641871 
    ## Run 244 stress 0.01574933 
    ## Run 245 stress 0.014986 
    ## Run 246 stress 0.01124043 
    ## ... Procrustes: rmse 0.007745886  max resid 0.01743036 
    ## Run 247 stress 0.0155792 
    ## Run 248 stress 0.01515091 
    ## Run 249 stress 0.01555896 
    ## Run 250 stress 0.01551821 
    ## Run 251 stress 0.014986 
    ## Run 252 stress 0.01547287 
    ## Run 253 stress 0.01558363 
    ## Run 254 stress 0.01564576 
    ## Run 255 stress 0.01129198 
    ## ... Procrustes: rmse 0.009617095  max resid 0.0216288 
    ## Run 256 stress 0.01493169 
    ## Run 257 stress 0.01529481 
    ## Run 258 stress 0.01115187 
    ## ... Procrustes: rmse 0.0002809005  max resid 0.0006343988 
    ## ... Similar to previous best
    ## Run 259 stress 0.01648685 
    ## Run 260 stress 0.01115364 
    ## ... Procrustes: rmse 0.001712353  max resid 0.003861008 
    ## ... Similar to previous best
    ## Run 261 stress 0.01562155 
    ## Run 262 stress 0.01546488 
    ## Run 263 stress 0.01553735 
    ## Run 264 stress 0.01549944 
    ## Run 265 stress 0.01561906 
    ## Run 266 stress 0.01542695 
    ## Run 267 stress 0.01543988 
    ## Run 268 stress 0.01498583 
    ## Run 269 stress 0.3596189 
    ## Run 270 stress 0.01115179 
    ## ... Procrustes: rmse 0.0002461162  max resid 0.0005556415 
    ## ... Similar to previous best
    ## Run 271 stress 0.0153375 
    ## Run 272 stress 0.01528172 
    ## Run 273 stress 0.01560102 
    ## Run 274 stress 0.0153612 
    ## Run 275 stress 0.01115192 
    ## ... Procrustes: rmse 0.00119552  max resid 0.002695107 
    ## ... Similar to previous best
    ## Run 276 stress 0.01493173 
    ## Run 277 stress 0.0155609 
    ## Run 278 stress 0.01549962 
    ## Run 279 stress 0.01529112 
    ## Run 280 stress 0.01547421 
    ## Run 281 stress 0.01544948 
    ## Run 282 stress 0.01663507 
    ## Run 283 stress 0.01547958 
    ## Run 284 stress 0.01558322 
    ## Run 285 stress 0.01528191 
    ## Run 286 stress 0.01685543 
    ## Run 287 stress 0.01603847 
    ## Run 288 stress 0.01115188 
    ## ... Procrustes: rmse 0.001190485  max resid 0.00268433 
    ## ... Similar to previous best
    ## Run 289 stress 0.01555562 
    ## Run 290 stress 0.01498578 
    ## Run 291 stress 0.01528186 
    ## Run 292 stress 0.01115173 
    ## ... Procrustes: rmse 0.001120453  max resid 0.002526013 
    ## ... Similar to previous best
    ## Run 293 stress 0.01701616 
    ## Run 294 stress 0.0111519 
    ## ... Procrustes: rmse 0.0002927705  max resid 0.0006613369 
    ## ... Similar to previous best
    ## Run 295 stress 0.01516727 
    ## Run 296 stress 0.01547397 
    ## Run 297 stress 0.01115171 
    ## ... Procrustes: rmse 0.0002093685  max resid 0.000472295 
    ## ... Similar to previous best
    ## Run 298 stress 0.01700746 
    ## Run 299 stress 0.0153379 
    ## Run 300 stress 0.01536087 
    ## Run 301 stress 0.01115182 
    ## ... Procrustes: rmse 0.0002608771  max resid 0.000588965 
    ## ... Similar to previous best
    ## Run 302 stress 0.01115178 
    ## ... Procrustes: rmse 0.0002417358  max resid 0.0005457025 
    ## ... Similar to previous best
    ## Run 303 stress 0.01499204 
    ## Run 304 stress 0.01641898 
    ## Run 305 stress 0.01641865 
    ## Run 306 stress 0.01558383 
    ## Run 307 stress 0.01637152 
    ## Run 308 stress 0.01115451 
    ## ... Procrustes: rmse 0.001902585  max resid 0.004289759 
    ## ... Similar to previous best
    ## Run 309 stress 0.015427 
    ## Run 310 stress 0.01498606 
    ## Run 311 stress 0.01565346 
    ## Run 312 stress 0.01562278 
    ## Run 313 stress 0.01685466 
    ## Run 314 stress 0.01574889 
    ## Run 315 stress 0.01556498 
    ## Run 316 stress 0.01549969 
    ## Run 317 stress 0.01493187 
    ## Run 318 stress 0.0156458 
    ## Run 319 stress 0.01536176 
    ## Run 320 stress 0.01528182 
    ## Run 321 stress 0.01576319 
    ## Run 322 stress 0.01543169 
    ## Run 323 stress 0.01169169 
    ## Run 324 stress 0.01557886 
    ## Run 325 stress 0.01556146 
    ## Run 326 stress 0.015282 
    ## Run 327 stress 0.01115199 
    ## ... Procrustes: rmse 0.001222427  max resid 0.002756801 
    ## ... Similar to previous best
    ## Run 328 stress 0.01582311 
    ## Run 329 stress 0.01115167 
    ## ... Procrustes: rmse 0.000191818  max resid 0.0004328184 
    ## ... Similar to previous best
    ## Run 330 stress 0.01498583 
    ## Run 331 stress 0.01127136 
    ## ... Procrustes: rmse 0.008689251  max resid 0.0195492 
    ## Run 332 stress 0.01115193 
    ## ... Procrustes: rmse 0.001203738  max resid 0.002713682 
    ## ... Similar to previous best
    ## Run 333 stress 0.01167486 
    ## Run 334 stress 0.0111847 
    ## ... Procrustes: rmse 0.004925687  max resid 0.0110986 
    ## Run 335 stress 0.01533731 
    ## Run 336 stress 0.01557909 
    ## Run 337 stress 0.01515734 
    ## Run 338 stress 0.0111598 
    ## ... Procrustes: rmse 0.002742968  max resid 0.006183616 
    ## ... Similar to previous best
    ## Run 339 stress 0.01542713 
    ## Run 340 stress 0.01115179 
    ## ... Procrustes: rmse 0.0002487674  max resid 0.000561715 
    ## ... Similar to previous best
    ## Run 341 stress 0.01281553 
    ## Run 342 stress 0.01542775 
    ## Run 343 stress 0.0154269 
    ## Run 344 stress 0.01521574 
    ## Run 345 stress 0.01576288 
    ## Run 346 stress 0.01565117 
    ## Run 347 stress 0.01685399 
    ## Run 348 stress 0.01129287 
    ## ... Procrustes: rmse 0.009486446  max resid 0.02133513 
    ## Run 349 stress 0.01564597 
    ## Run 350 stress 0.01557908 
    ## Run 351 stress 0.01542815 
    ## Run 352 stress 0.01115181 
    ## ... Procrustes: rmse 0.0002481781  max resid 0.0005599547 
    ## ... Similar to previous best
    ## Run 353 stress 0.01493149 
    ## Run 354 stress 0.01571579 
    ## Run 355 stress 0.01529066 
    ## Run 356 stress 0.01323099 
    ## Run 357 stress 0.01558342 
    ## Run 358 stress 0.01121081 
    ## ... Procrustes: rmse 0.006407885  max resid 0.01442465 
    ## Run 359 stress 0.0157468 
    ## Run 360 stress 0.0156741 
    ## Run 361 stress 0.0155317 
    ## Run 362 stress 0.01498607 
    ## Run 363 stress 0.01528184 
    ## Run 364 stress 0.01543168 
    ## Run 365 stress 0.01547113 
    ## Run 366 stress 0.01582054 
    ## Run 367 stress 0.01515716 
    ## Run 368 stress 0.01546679 
    ## Run 369 stress 0.01558439 
    ## Run 370 stress 0.01533784 
    ## Run 371 stress 0.01533917 
    ## Run 372 stress 0.01131057 
    ## ... Procrustes: rmse 0.0101486  max resid 0.02281771 
    ## Run 373 stress 0.01545671 
    ## Run 374 stress 0.0169457 
    ## Run 375 stress 0.01115139 
    ## ... Procrustes: rmse 4.92216e-05  max resid 0.000110336 
    ## ... Similar to previous best
    ## Run 376 stress 0.0163711 
    ## Run 377 stress 0.01582695 
    ## Run 378 stress 0.01498581 
    ## Run 379 stress 0.01528174 
    ## Run 380 stress 0.01498584 
    ## Run 381 stress 0.01498607 
    ## Run 382 stress 0.01131737 
    ## ... Procrustes: rmse 0.01025812  max resid 0.02304453 
    ## Run 383 stress 0.01545279 
    ## Run 384 stress 0.01547942 
    ## Run 385 stress 0.01557888 
    ## Run 386 stress 0.01116369 
    ## ... Procrustes: rmse 0.00313798  max resid 0.007075381 
    ## ... Similar to previous best
    ## Run 387 stress 0.0149859 
    ## Run 388 stress 0.01115167 
    ## ... Procrustes: rmse 0.001100527  max resid 0.002481481 
    ## ... Similar to previous best
    ## Run 389 stress 0.01574915 
    ## Run 390 stress 0.01150849 
    ## ... Procrustes: rmse 0.013708  max resid 0.0297556 
    ## Run 391 stress 0.0149317 
    ## Run 392 stress 0.01567222 
    ## Run 393 stress 0.01498598 
    ## Run 394 stress 0.01551001 
    ## Run 395 stress 0.01567515 
    ## Run 396 stress 0.01498573 
    ## Run 397 stress 0.01556911 
    ## Run 398 stress 0.01675224 
    ## Run 399 stress 0.01669122 
    ## Run 400 stress 0.01515729 
    ## Run 401 stress 0.01549911 
    ## Run 402 stress 0.01147843 
    ## ... Procrustes: rmse 0.01284297  max resid 0.0280455 
    ## Run 403 stress 0.01156638 
    ## ... Procrustes: rmse 0.01527343  max resid 0.03301813 
    ## Run 404 stress 0.01122002 
    ## ... Procrustes: rmse 0.006854458  max resid 0.01542593 
    ## Run 405 stress 0.01119718 
    ## ... Procrustes: rmse 0.005687605  max resid 0.01280922 
    ## Run 406 stress 0.0112283 
    ## ... Procrustes: rmse 0.007228676  max resid 0.01626547 
    ## Run 407 stress 0.01499859 
    ## Run 408 stress 0.01581902 
    ## Run 409 stress 0.01557908 
    ## Run 410 stress 0.01120182 
    ## ... Procrustes: rmse 0.005611171  max resid 0.01264591 
    ## Run 411 stress 0.01515659 
    ## Run 412 stress 0.01557893 
    ## Run 413 stress 0.01547631 
    ## Run 414 stress 0.01557886 
    ## Run 415 stress 0.01549958 
    ## Run 416 stress 0.01118504 
    ## ... Procrustes: rmse 0.004910317  max resid 0.01106642 
    ## Run 417 stress 0.01567262 
    ## Run 418 stress 0.01528196 
    ## Run 419 stress 0.01525791 
    ## Run 420 stress 0.01553496 
    ## Run 421 stress 0.01493172 
    ## Run 422 stress 0.01116241 
    ## ... Procrustes: rmse 0.003059454  max resid 0.006896647 
    ## ... Similar to previous best
    ## Run 423 stress 0.01700867 
    ## Run 424 stress 0.01121706 
    ## ... Procrustes: rmse 0.006715212  max resid 0.01511355 
    ## Run 425 stress 0.01552955 
    ## Run 426 stress 0.0151612 
    ## Run 427 stress 0.0111517 
    ## ... Procrustes: rmse 0.0002108269  max resid 0.0004757402 
    ## ... Similar to previous best
    ## Run 428 stress 0.01544558 
    ## Run 429 stress 0.01499154 
    ## Run 430 stress 0.01689046 
    ## Run 431 stress 0.01581908 
    ## Run 432 stress 0.01544696 
    ## Run 433 stress 0.01268837 
    ## Run 434 stress 0.01668202 
    ## Run 435 stress 0.01533762 
    ## Run 436 stress 0.01545214 
    ## Run 437 stress 0.01498582 
    ## Run 438 stress 0.01336293 
    ## Run 439 stress 0.01655634 
    ## Run 440 stress 0.01637117 
    ## Run 441 stress 0.01115191 
    ## ... Procrustes: rmse 0.0002840712  max resid 0.0006417182 
    ## ... Similar to previous best
    ## Run 442 stress 0.01558503 
    ## Run 443 stress 0.01515723 
    ## Run 444 stress 0.01545332 
    ## Run 445 stress 0.01162297 
    ## ... Procrustes: rmse 0.01704769  max resid 0.03682631 
    ## Run 446 stress 0.01150189 
    ## ... Procrustes: rmse 0.01351497  max resid 0.0293666 
    ## Run 447 stress 0.01498581 
    ## Run 448 stress 0.01498572 
    ## Run 449 stress 0.01314065 
    ## Run 450 stress 0.0154723 
    ## Run 451 stress 0.01498655 
    ## Run 452 stress 0.01500106 
    ## Run 453 stress 0.01529118 
    ## Run 454 stress 0.01545394 
    ## Run 455 stress 0.01295342 
    ## Run 456 stress 0.01545337 
    ## Run 457 stress 0.01545285 
    ## Run 458 stress 0.01115169 
    ## ... Procrustes: rmse 0.001109138  max resid 0.002500879 
    ## ... Similar to previous best
    ## Run 459 stress 0.01558673 
    ## Run 460 stress 0.01528207 
    ## Run 461 stress 0.01546286 
    ## Run 462 stress 0.01145771 
    ## ... Procrustes: rmse 0.01230885  max resid 0.02701991 
    ## Run 463 stress 0.01493167 
    ## Run 464 stress 0.01556088 
    ## Run 465 stress 0.01115184 
    ## ... Procrustes: rmse 0.0002657269  max resid 0.0006001629 
    ## ... Similar to previous best
    ## Run 466 stress 0.01115166 
    ## ... Procrustes: rmse 0.000190634  max resid 0.0004300297 
    ## ... Similar to previous best
    ## Run 467 stress 0.01567528 
    ## Run 468 stress 0.01515091 
    ## Run 469 stress 0.01127861 
    ## ... Procrustes: rmse 0.009172906  max resid 0.02062744 
    ## Run 470 stress 0.2893968 
    ## Run 471 stress 0.0154475 
    ## Run 472 stress 0.01120963 
    ## ... Procrustes: rmse 0.006322598  max resid 0.01423493 
    ## Run 473 stress 0.37386 
    ## Run 474 stress 0.01694553 
    ## Run 475 stress 0.01567114 
    ## Run 476 stress 0.01115171 
    ## ... Procrustes: rmse 0.0002108148  max resid 0.0004754693 
    ## ... Similar to previous best
    ## Run 477 stress 0.01169622 
    ## Run 478 stress 0.01527945 
    ## Run 479 stress 0.015748 
    ## Run 480 stress 0.0111518 
    ## ... Procrustes: rmse 0.0002532976  max resid 0.0005719316 
    ## ... Similar to previous best
    ## Run 481 stress 0.0149922 
    ## Run 482 stress 0.01577212 
    ## Run 483 stress 0.01499257 
    ## Run 484 stress 0.0154991 
    ## Run 485 stress 0.0152635 
    ## Run 486 stress 0.0154278 
    ## Run 487 stress 0.01663487 
    ## Run 488 stress 0.0156522 
    ## Run 489 stress 0.01552969 
    ## Run 490 stress 0.01668205 
    ## Run 491 stress 0.01115172 
    ## ... Procrustes: rmse 0.001121873  max resid 0.002529473 
    ## ... Similar to previous best
    ## Run 492 stress 0.01117061 
    ## ... Procrustes: rmse 0.003866979  max resid 0.00871533 
    ## ... Similar to previous best
    ## Run 493 stress 0.01122923 
    ## ... Procrustes: rmse 0.007264421  max resid 0.01634769 
    ## Run 494 stress 0.01115189 
    ## ... Procrustes: rmse 0.0002823547  max resid 0.0006379482 
    ## ... Similar to previous best
    ## Run 495 stress 0.01544718 
    ## Run 496 stress 0.01514993 
    ## Run 497 stress 0.0152817 
    ## Run 498 stress 0.01527269 
    ## Run 499 stress 0.01555194 
    ## Run 500 stress 0.01515704 
    ## *** Best solution repeated 44 times

``` r
# Mixed and stratified lakes
PD_beta_MS_NMDS <- metaMDS(PD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.06111765 
    ## Run 2 stress 0.05171062 
    ## Run 3 stress 0.06499235 
    ## Run 4 stress 0.05051307 
    ## ... Procrustes: rmse 0.03573358  max resid 0.1219801 
    ## Run 5 stress 0.05337579 
    ## Run 6 stress 0.06493714 
    ## Run 7 stress 0.06231288 
    ## Run 8 stress 0.06851823 
    ## Run 9 stress 0.05036798 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000285822  max resid 0.0007239979 
    ## ... Similar to previous best
    ## Run 10 stress 0.05171058 
    ## Run 11 stress 0.05036796 
    ## ... New best solution
    ## ... Procrustes: rmse 1.666657e-05  max resid 4.51422e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.05900148 
    ## Run 13 stress 0.05051351 
    ## ... Procrustes: rmse 0.03564503  max resid 0.1217999 
    ## Run 14 stress 0.05051327 
    ## ... Procrustes: rmse 0.03564209  max resid 0.1216207 
    ## Run 15 stress 0.05036809 
    ## ... Procrustes: rmse 0.0003108198  max resid 0.0008099815 
    ## ... Similar to previous best
    ## Run 16 stress 0.05036804 
    ## ... Procrustes: rmse 6.299391e-05  max resid 0.0001654232 
    ## ... Similar to previous best
    ## Run 17 stress 0.06478527 
    ## Run 18 stress 0.06136691 
    ## Run 19 stress 0.05968631 
    ## Run 20 stress 0.05829935 
    ## Run 21 stress 0.0689694 
    ## Run 22 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001803967  max resid 0.0005263452 
    ## ... Similar to previous best
    ## Run 23 stress 0.05051316 
    ## ... Procrustes: rmse 0.03569116  max resid 0.1217792 
    ## Run 24 stress 0.05968643 
    ## Run 25 stress 0.0517106 
    ## Run 26 stress 0.05171059 
    ## Run 27 stress 0.05036806 
    ## ... Procrustes: rmse 8.057399e-05  max resid 0.0001982474 
    ## ... Similar to previous best
    ## Run 28 stress 0.05051312 
    ## ... Procrustes: rmse 0.0356467  max resid 0.1216546 
    ## Run 29 stress 0.06478533 
    ## Run 30 stress 0.05036794 
    ## ... New best solution
    ## ... Procrustes: rmse 2.489239e-05  max resid 5.932836e-05 
    ## ... Similar to previous best
    ## Run 31 stress 0.05337583 
    ## Run 32 stress 0.05051337 
    ## ... Procrustes: rmse 0.03569583  max resid 0.121955 
    ## Run 33 stress 0.050368 
    ## ... Procrustes: rmse 0.0001927692  max resid 0.0005640433 
    ## ... Similar to previous best
    ## Run 34 stress 0.05036793 
    ## ... New best solution
    ## ... Procrustes: rmse 2.333587e-05  max resid 6.740214e-05 
    ## ... Similar to previous best
    ## Run 35 stress 0.05051318 
    ## ... Procrustes: rmse 0.03565791  max resid 0.1218368 
    ## Run 36 stress 0.06476226 
    ## Run 37 stress 0.0596864 
    ## Run 38 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001509024  max resid 0.0004052783 
    ## ... Similar to previous best
    ## Run 39 stress 0.05171059 
    ## Run 40 stress 0.06136689 
    ## Run 41 stress 0.05602057 
    ## Run 42 stress 0.05036801 
    ## ... Procrustes: rmse 7.668078e-05  max resid 0.0002009239 
    ## ... Similar to previous best
    ## Run 43 stress 0.05968642 
    ## Run 44 stress 0.05337584 
    ## Run 45 stress 0.05036795 
    ## ... Procrustes: rmse 9.465495e-05  max resid 0.0002675481 
    ## ... Similar to previous best
    ## Run 46 stress 0.06038949 
    ## Run 47 stress 0.05036793 
    ## ... New best solution
    ## ... Procrustes: rmse 3.495555e-05  max resid 6.710758e-05 
    ## ... Similar to previous best
    ## Run 48 stress 0.06938796 
    ## Run 49 stress 0.06609067 
    ## Run 50 stress 0.0533758 
    ## Run 51 stress 0.06499241 
    ## Run 52 stress 0.05602071 
    ## Run 53 stress 0.0613668 
    ## Run 54 stress 0.0505135 
    ## ... Procrustes: rmse 0.03568051  max resid 0.1219204 
    ## Run 55 stress 0.05036799 
    ## ... Procrustes: rmse 9.493074e-05  max resid 0.0002389863 
    ## ... Similar to previous best
    ## Run 56 stress 0.05337577 
    ## Run 57 stress 0.06499234 
    ## Run 58 stress 0.06872752 
    ## Run 59 stress 0.05602066 
    ## Run 60 stress 0.05036792 
    ## ... New best solution
    ## ... Procrustes: rmse 7.645765e-05  max resid 0.0001809517 
    ## ... Similar to previous best
    ## Run 61 stress 0.05466702 
    ## Run 62 stress 0.06904456 
    ## Run 63 stress 0.05171059 
    ## Run 64 stress 0.0630985 
    ## Run 65 stress 0.05036792 
    ## ... Procrustes: rmse 8.95141e-05  max resid 0.0002132005 
    ## ... Similar to previous best
    ## Run 66 stress 0.05036793 
    ## ... Procrustes: rmse 0.0001060521  max resid 0.0002556552 
    ## ... Similar to previous best
    ## Run 67 stress 0.0590015 
    ## Run 68 stress 0.05051354 
    ## ... Procrustes: rmse 0.03567128  max resid 0.1217184 
    ## Run 69 stress 0.0636279 
    ## Run 70 stress 0.05036794 
    ## ... Procrustes: rmse 5.048298e-05  max resid 0.0001284912 
    ## ... Similar to previous best
    ## Run 71 stress 0.3527961 
    ## Run 72 stress 0.06476218 
    ## Run 73 stress 0.06499237 
    ## Run 74 stress 0.05171059 
    ## Run 75 stress 0.05051342 
    ## ... Procrustes: rmse 0.03574341  max resid 0.1219456 
    ## Run 76 stress 0.05829939 
    ## Run 77 stress 0.0620301 
    ## Run 78 stress 0.05403611 
    ## Run 79 stress 0.06504218 
    ## Run 80 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001375417  max resid 0.0003362724 
    ## ... Similar to previous best
    ## Run 81 stress 0.3406432 
    ## Run 82 stress 0.05051357 
    ## ... Procrustes: rmse 0.03572176  max resid 0.1220594 
    ## Run 83 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001744903  max resid 0.0004236938 
    ## ... Similar to previous best
    ## Run 84 stress 0.05466712 
    ## Run 85 stress 0.05171058 
    ## Run 86 stress 0.05829935 
    ## Run 87 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001456245  max resid 0.0003444097 
    ## ... Similar to previous best
    ## Run 88 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001233941  max resid 0.0003399768 
    ## ... Similar to previous best
    ## Run 89 stress 0.05051306 
    ## ... Procrustes: rmse 0.03570843  max resid 0.1220008 
    ## Run 90 stress 0.050368 
    ## ... Procrustes: rmse 0.0001817647  max resid 0.0004402977 
    ## ... Similar to previous best
    ## Run 91 stress 0.05036811 
    ## ... Procrustes: rmse 0.0002378672  max resid 0.0005832991 
    ## ... Similar to previous best
    ## Run 92 stress 0.05036796 
    ## ... Procrustes: rmse 7.176296e-05  max resid 0.0002085941 
    ## ... Similar to previous best
    ## Run 93 stress 0.06309854 
    ## Run 94 stress 0.06309872 
    ## Run 95 stress 0.06640334 
    ## Run 96 stress 0.05051353 
    ## ... Procrustes: rmse 0.03575072  max resid 0.1221477 
    ## Run 97 stress 0.05171059 
    ## Run 98 stress 0.0650424 
    ## Run 99 stress 0.05928214 
    ## Run 100 stress 0.05051322 
    ## ... Procrustes: rmse 0.03573866  max resid 0.1220995 
    ## Run 101 stress 0.05337582 
    ## Run 102 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001439792  max resid 0.0003409784 
    ## ... Similar to previous best
    ## Run 103 stress 0.06478537 
    ## Run 104 stress 0.05337585 
    ## Run 105 stress 0.05968629 
    ## Run 106 stress 0.05036794 
    ## ... Procrustes: rmse 0.0001174392  max resid 0.0002447204 
    ## ... Similar to previous best
    ## Run 107 stress 0.05171059 
    ## Run 108 stress 0.05036794 
    ## ... Procrustes: rmse 0.0001188174  max resid 0.0002864227 
    ## ... Similar to previous best
    ## Run 109 stress 0.05171058 
    ## Run 110 stress 0.05051358 
    ## ... Procrustes: rmse 0.03574747  max resid 0.1221374 
    ## Run 111 stress 0.0505134 
    ## ... Procrustes: rmse 0.03568278  max resid 0.1217652 
    ## Run 112 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001510169  max resid 0.0003693829 
    ## ... Similar to previous best
    ## Run 113 stress 0.05968625 
    ## Run 114 stress 0.05900154 
    ## Run 115 stress 0.06136679 
    ## Run 116 stress 0.06136677 
    ## Run 117 stress 0.05466701 
    ## Run 118 stress 0.06537616 
    ## Run 119 stress 0.05403554 
    ## Run 120 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001374304  max resid 0.0003020921 
    ## ... Similar to previous best
    ## Run 121 stress 0.0582993 
    ## Run 122 stress 0.05968636 
    ## Run 123 stress 0.05829936 
    ## Run 124 stress 0.06027436 
    ## Run 125 stress 0.05337583 
    ## Run 126 stress 0.050368 
    ## ... Procrustes: rmse 0.0001723607  max resid 0.0003928413 
    ## ... Similar to previous best
    ## Run 127 stress 0.05051304 
    ## ... Procrustes: rmse 0.03572115  max resid 0.1220381 
    ## Run 128 stress 0.05829943 
    ## Run 129 stress 0.05051284 
    ## ... Procrustes: rmse 0.03563834  max resid 0.1217709 
    ## Run 130 stress 0.05968638 
    ## Run 131 stress 0.05051303 
    ## ... Procrustes: rmse 0.03573144  max resid 0.1220676 
    ## Run 132 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001930058  max resid 0.000445772 
    ## ... Similar to previous best
    ## Run 133 stress 0.05036795 
    ## ... Procrustes: rmse 0.000136562  max resid 0.0003322182 
    ## ... Similar to previous best
    ## Run 134 stress 0.05337591 
    ## Run 135 stress 0.06478531 
    ## Run 136 stress 0.06476216 
    ## Run 137 stress 0.06027432 
    ## Run 138 stress 0.0596863 
    ## Run 139 stress 0.0636279 
    ## Run 140 stress 0.05968625 
    ## Run 141 stress 0.3592593 
    ## Run 142 stress 0.053376 
    ## Run 143 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001591977  max resid 0.0003875727 
    ## ... Similar to previous best
    ## Run 144 stress 0.05036797 
    ## ... Procrustes: rmse 8.679497e-05  max resid 0.0002487873 
    ## ... Similar to previous best
    ## Run 145 stress 0.05051337 
    ## ... Procrustes: rmse 0.03574201  max resid 0.1221134 
    ## Run 146 stress 0.05466702 
    ## Run 147 stress 0.06038952 
    ## Run 148 stress 0.05928211 
    ## Run 149 stress 0.05466708 
    ## Run 150 stress 0.06136687 
    ## Run 151 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 6.492572e-05  max resid 0.000155574 
    ## ... Similar to previous best
    ## Run 152 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001451848  max resid 0.0003666491 
    ## ... Similar to previous best
    ## Run 153 stress 0.05051311 
    ## ... Procrustes: rmse 0.03567637  max resid 0.1218919 
    ## Run 154 stress 0.05036797 
    ## ... Procrustes: rmse 8.875923e-05  max resid 0.0002214856 
    ## ... Similar to previous best
    ## Run 155 stress 0.06231287 
    ## Run 156 stress 0.05036796 
    ## ... Procrustes: rmse 8.892814e-05  max resid 0.0001785585 
    ## ... Similar to previous best
    ## Run 157 stress 0.07242246 
    ## Run 158 stress 0.05466704 
    ## Run 159 stress 0.05171058 
    ## Run 160 stress 0.05036793 
    ## ... Procrustes: rmse 2.137656e-05  max resid 4.062657e-05 
    ## ... Similar to previous best
    ## Run 161 stress 0.05337593 
    ## Run 162 stress 0.05036792 
    ## ... Procrustes: rmse 6.826731e-05  max resid 0.0001678281 
    ## ... Similar to previous best
    ## Run 163 stress 0.05337577 
    ## Run 164 stress 0.05171059 
    ## Run 165 stress 0.05051302 
    ## ... Procrustes: rmse 0.03568796  max resid 0.1218036 
    ## Run 166 stress 0.05466708 
    ## Run 167 stress 0.05403585 
    ## Run 168 stress 0.05036794 
    ## ... Procrustes: rmse 5.242942e-05  max resid 0.0001281966 
    ## ... Similar to previous best
    ## Run 169 stress 0.05036795 
    ## ... Procrustes: rmse 2.435339e-05  max resid 5.060753e-05 
    ## ... Similar to previous best
    ## Run 170 stress 0.05171059 
    ## Run 171 stress 0.05051286 
    ## ... Procrustes: rmse 0.03568334  max resid 0.1218124 
    ## Run 172 stress 0.05602068 
    ## Run 173 stress 0.05036792 
    ## ... Procrustes: rmse 1.272502e-05  max resid 2.428605e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.05403606 
    ## Run 175 stress 0.05036798 
    ## ... Procrustes: rmse 6.447036e-05  max resid 0.0001578902 
    ## ... Similar to previous best
    ## Run 176 stress 0.05602052 
    ## Run 177 stress 0.05602043 
    ## Run 178 stress 0.05036796 
    ## ... Procrustes: rmse 8.652175e-05  max resid 0.0002102771 
    ## ... Similar to previous best
    ## Run 179 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001413883  max resid 0.0003557874 
    ## ... Similar to previous best
    ## Run 180 stress 0.05337595 
    ## Run 181 stress 0.05051291 
    ## ... Procrustes: rmse 0.03567487  max resid 0.1218787 
    ## Run 182 stress 0.05829933 
    ## Run 183 stress 0.0590015 
    ## Run 184 stress 0.05051295 
    ## ... Procrustes: rmse 0.03568361  max resid 0.1219074 
    ## Run 185 stress 0.05051309 
    ## ... Procrustes: rmse 0.03569488  max resid 0.1219484 
    ## Run 186 stress 0.05051322 
    ## ... Procrustes: rmse 0.03571729  max resid 0.1220239 
    ## Run 187 stress 0.05602045 
    ## Run 188 stress 0.05036795 
    ## ... Procrustes: rmse 6.930302e-05  max resid 0.0001634947 
    ## ... Similar to previous best
    ## Run 189 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001306989  max resid 0.0003434437 
    ## ... Similar to previous best
    ## Run 190 stress 0.05829935 
    ## Run 191 stress 0.3525959 
    ## Run 192 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001313947  max resid 0.0003306111 
    ## ... Similar to previous best
    ## Run 193 stress 0.05602057 
    ## Run 194 stress 0.05900152 
    ## Run 195 stress 0.0505132 
    ## ... Procrustes: rmse 0.03570386  max resid 0.1219832 
    ## Run 196 stress 0.05403623 
    ## Run 197 stress 0.05036793 
    ## ... Procrustes: rmse 7.892279e-05  max resid 0.0001943828 
    ## ... Similar to previous best
    ## Run 198 stress 0.05968627 
    ## Run 199 stress 0.05051297 
    ## ... Procrustes: rmse 0.03574784  max resid 0.1220978 
    ## Run 200 stress 0.06136674 
    ## Run 201 stress 0.05036795 
    ## ... Procrustes: rmse 6.638797e-05  max resid 0.0001595787 
    ## ... Similar to previous best
    ## Run 202 stress 0.05602061 
    ## Run 203 stress 0.05051332 
    ## ... Procrustes: rmse 0.03569076  max resid 0.121949 
    ## Run 204 stress 0.06478512 
    ## Run 205 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001327816  max resid 0.000334404 
    ## ... Similar to previous best
    ## Run 206 stress 0.05036793 
    ## ... Procrustes: rmse 9.162872e-05  max resid 0.0002261124 
    ## ... Similar to previous best
    ## Run 207 stress 0.05466714 
    ## Run 208 stress 0.05928215 
    ## Run 209 stress 0.05337584 
    ## Run 210 stress 0.05403594 
    ## Run 211 stress 0.05829933 
    ## Run 212 stress 0.06203071 
    ## Run 213 stress 0.05337587 
    ## Run 214 stress 0.05829935 
    ## Run 215 stress 0.05171058 
    ## Run 216 stress 0.06832165 
    ## Run 217 stress 0.06647211 
    ## Run 218 stress 0.05171061 
    ## Run 219 stress 0.05403615 
    ## Run 220 stress 0.05171058 
    ## Run 221 stress 0.05171059 
    ## Run 222 stress 0.05036795 
    ## ... Procrustes: rmse 4.632539e-05  max resid 0.0001079452 
    ## ... Similar to previous best
    ## Run 223 stress 0.05051348 
    ## ... Procrustes: rmse 0.03572424  max resid 0.1220533 
    ## Run 224 stress 0.05051341 
    ## ... Procrustes: rmse 0.03572696  max resid 0.1220605 
    ## Run 225 stress 0.05051302 
    ## ... Procrustes: rmse 0.03567924  max resid 0.1218966 
    ## Run 226 stress 0.05036798 
    ## ... Procrustes: rmse 9.741131e-05  max resid 0.0002135677 
    ## ... Similar to previous best
    ## Run 227 stress 0.0560207 
    ## Run 228 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001492695  max resid 0.0003702088 
    ## ... Similar to previous best
    ## Run 229 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001357878  max resid 0.0003866199 
    ## ... Similar to previous best
    ## Run 230 stress 0.0582993 
    ## Run 231 stress 0.05051314 
    ## ... Procrustes: rmse 0.03567147  max resid 0.1217417 
    ## Run 232 stress 0.05051359 
    ## ... Procrustes: rmse 0.03564277  max resid 0.1216174 
    ## Run 233 stress 0.05036794 
    ## ... Procrustes: rmse 5.443482e-05  max resid 0.00013325 
    ## ... Similar to previous best
    ## Run 234 stress 0.0505131 
    ## ... Procrustes: rmse 0.03565711  max resid 0.1217011 
    ## Run 235 stress 0.05036793 
    ## ... Procrustes: rmse 3.374604e-05  max resid 7.596949e-05 
    ## ... Similar to previous best
    ## Run 236 stress 0.05829929 
    ## Run 237 stress 0.06362789 
    ## Run 238 stress 0.06111771 
    ## Run 239 stress 0.05466705 
    ## Run 240 stress 0.050368 
    ## ... Procrustes: rmse 0.0001696631  max resid 0.0004760707 
    ## ... Similar to previous best
    ## Run 241 stress 0.05036793 
    ## ... Procrustes: rmse 9.108525e-05  max resid 0.0002519182 
    ## ... Similar to previous best
    ## Run 242 stress 0.05051304 
    ## ... Procrustes: rmse 0.03569016  max resid 0.1219331 
    ## Run 243 stress 0.05036794 
    ## ... Procrustes: rmse 0.0001045882  max resid 0.000237501 
    ## ... Similar to previous best
    ## Run 244 stress 0.05171064 
    ## Run 245 stress 0.06478532 
    ## Run 246 stress 0.05403589 
    ## Run 247 stress 0.05829929 
    ## Run 248 stress 0.05968627 
    ## Run 249 stress 0.05337581 
    ## Run 250 stress 0.05036793 
    ## ... Procrustes: rmse 0.0001043194  max resid 0.0002597947 
    ## ... Similar to previous best
    ## Run 251 stress 0.05337593 
    ## Run 252 stress 0.06938771 
    ## Run 253 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001585923  max resid 0.00042163 
    ## ... Similar to previous best
    ## Run 254 stress 0.05337577 
    ## Run 255 stress 0.06647208 
    ## Run 256 stress 0.05051321 
    ## ... Procrustes: rmse 0.03565244  max resid 0.1216756 
    ## Run 257 stress 0.0517106 
    ## Run 258 stress 0.06136684 
    ## Run 259 stress 0.05466702 
    ## Run 260 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001591371  max resid 0.0003872917 
    ## ... Similar to previous best
    ## Run 261 stress 0.05403576 
    ## Run 262 stress 0.0667456 
    ## Run 263 stress 0.05829929 
    ## Run 264 stress 0.05051347 
    ## ... Procrustes: rmse 0.035694  max resid 0.1219642 
    ## Run 265 stress 0.05051317 
    ## ... Procrustes: rmse 0.0357858  max resid 0.122099 
    ## Run 266 stress 0.06826688 
    ## Run 267 stress 0.05051339 
    ## ... Procrustes: rmse 0.03562732  max resid 0.1217537 
    ## Run 268 stress 0.05036805 
    ## ... Procrustes: rmse 0.0002032471  max resid 0.0005375438 
    ## ... Similar to previous best
    ## Run 269 stress 0.05403645 
    ## Run 270 stress 0.06844316 
    ## Run 271 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001051542  max resid 0.0002655678 
    ## ... Similar to previous best
    ## Run 272 stress 0.05036806 
    ## ... Procrustes: rmse 0.0002107881  max resid 0.0005674237 
    ## ... Similar to previous best
    ## Run 273 stress 0.06478525 
    ## Run 274 stress 0.05051304 
    ## ... Procrustes: rmse 0.0356908  max resid 0.1219357 
    ## Run 275 stress 0.05968629 
    ## Run 276 stress 0.050368 
    ## ... Procrustes: rmse 0.000104704  max resid 0.0002311075 
    ## ... Similar to previous best
    ## Run 277 stress 0.05036794 
    ## ... Procrustes: rmse 6.648046e-05  max resid 0.0001649557 
    ## ... Similar to previous best
    ## Run 278 stress 0.050368 
    ## ... Procrustes: rmse 0.0001614755  max resid 0.0004629189 
    ## ... Similar to previous best
    ## Run 279 stress 0.05900156 
    ## Run 280 stress 0.05829933 
    ## Run 281 stress 0.06111767 
    ## Run 282 stress 0.05051287 
    ## ... Procrustes: rmse 0.03567043  max resid 0.1218588 
    ## Run 283 stress 0.06478527 
    ## Run 284 stress 0.05403578 
    ## Run 285 stress 0.05337582 
    ## Run 286 stress 0.05051296 
    ## ... Procrustes: rmse 0.03554754  max resid 0.1214124 
    ## Run 287 stress 0.05036794 
    ## ... Procrustes: rmse 0.0001172725  max resid 0.000283432 
    ## ... Similar to previous best
    ## Run 288 stress 0.05829944 
    ## Run 289 stress 0.05051315 
    ## ... Procrustes: rmse 0.03566197  max resid 0.1217109 
    ## Run 290 stress 0.06609064 
    ## Run 291 stress 0.05171058 
    ## Run 292 stress 0.05466702 
    ## Run 293 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001590302  max resid 0.0004090384 
    ## ... Similar to previous best
    ## Run 294 stress 0.05051316 
    ## ... Procrustes: rmse 0.03562506  max resid 0.121746 
    ## Run 295 stress 0.06938751 
    ## Run 296 stress 0.05337576 
    ## Run 297 stress 0.05171059 
    ## Run 298 stress 0.05036794 
    ## ... Procrustes: rmse 5.602338e-05  max resid 0.0001368728 
    ## ... Similar to previous best
    ## Run 299 stress 0.05171058 
    ## Run 300 stress 0.06027435 
    ## Run 301 stress 0.05036797 
    ## ... Procrustes: rmse 5.284471e-05  max resid 0.0001221675 
    ## ... Similar to previous best
    ## Run 302 stress 0.050368 
    ## ... Procrustes: rmse 0.0001245989  max resid 0.0003103697 
    ## ... Similar to previous best
    ## Run 303 stress 0.06027426 
    ## Run 304 stress 0.05051313 
    ## ... Procrustes: rmse 0.03568483  max resid 0.1217854 
    ## Run 305 stress 0.06904496 
    ## Run 306 stress 0.07033252 
    ## Run 307 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001374254  max resid 0.000342154 
    ## ... Similar to previous best
    ## Run 308 stress 0.0596863 
    ## Run 309 stress 0.06136678 
    ## Run 310 stress 0.05900153 
    ## Run 311 stress 0.06943828 
    ## Run 312 stress 0.05036807 
    ## ... Procrustes: rmse 0.0001647503  max resid 0.0004135028 
    ## ... Similar to previous best
    ## Run 313 stress 0.05036794 
    ## ... Procrustes: rmse 0.0001076068  max resid 0.0002898213 
    ## ... Similar to previous best
    ## Run 314 stress 0.05602053 
    ## Run 315 stress 0.05051328 
    ## ... Procrustes: rmse 0.03568832  max resid 0.1219402 
    ## Run 316 stress 0.06951046 
    ## Run 317 stress 0.050368 
    ## ... Procrustes: rmse 0.0001688071  max resid 0.0004462034 
    ## ... Similar to previous best
    ## Run 318 stress 0.05171062 
    ## Run 319 stress 0.05051287 
    ## ... Procrustes: rmse 0.03567827  max resid 0.1218852 
    ## Run 320 stress 0.05602048 
    ## Run 321 stress 0.05900147 
    ## Run 322 stress 0.05051339 
    ## ... Procrustes: rmse 0.03572087  max resid 0.1220416 
    ## Run 323 stress 0.05036796 
    ## ... Procrustes: rmse 7.555398e-05  max resid 0.0001864849 
    ## ... Similar to previous best
    ## Run 324 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 2.951638e-05  max resid 8.617958e-05 
    ## ... Similar to previous best
    ## Run 325 stress 0.06943801 
    ## Run 326 stress 0.050368 
    ## ... Procrustes: rmse 0.0001464738  max resid 0.0003750792 
    ## ... Similar to previous best
    ## Run 327 stress 0.05036793 
    ## ... Procrustes: rmse 7.018044e-05  max resid 0.0001452261 
    ## ... Similar to previous best
    ## Run 328 stress 0.05051305 
    ## ... Procrustes: rmse 0.03567864  max resid 0.1217794 
    ## Run 329 stress 0.0505133 
    ## ... Procrustes: rmse 0.03568413  max resid 0.1217699 
    ## Run 330 stress 0.05051333 
    ## ... Procrustes: rmse 0.035731  max resid 0.1220742 
    ## Run 331 stress 0.05928233 
    ## Run 332 stress 0.05051289 
    ## ... Procrustes: rmse 0.03567873  max resid 0.1218004 
    ## Run 333 stress 0.05171063 
    ## Run 334 stress 0.05036792 
    ## ... Procrustes: rmse 4.456606e-05  max resid 0.0001261403 
    ## ... Similar to previous best
    ## Run 335 stress 0.05036793 
    ## ... Procrustes: rmse 6.161923e-05  max resid 0.0001680861 
    ## ... Similar to previous best
    ## Run 336 stress 0.05403626 
    ## Run 337 stress 0.0582994 
    ## Run 338 stress 0.05036797 
    ## ... Procrustes: rmse 0.000117114  max resid 0.0002933889 
    ## ... Similar to previous best
    ## Run 339 stress 0.05171062 
    ## Run 340 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001711752  max resid 0.0004431036 
    ## ... Similar to previous best
    ## Run 341 stress 0.05036794 
    ## ... Procrustes: rmse 7.233796e-05  max resid 0.0001571117 
    ## ... Similar to previous best
    ## Run 342 stress 0.05036795 
    ## ... Procrustes: rmse 9.692506e-05  max resid 0.0002535393 
    ## ... Similar to previous best
    ## Run 343 stress 0.0505129 
    ## ... Procrustes: rmse 0.0356808  max resid 0.1219008 
    ## Run 344 stress 0.05466702 
    ## Run 345 stress 0.0582993 
    ## Run 346 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001693565  max resid 0.0004485246 
    ## ... Similar to previous best
    ## Run 347 stress 0.05051321 
    ## ... Procrustes: rmse 0.03565399  max resid 0.1216927 
    ## Run 348 stress 0.05171066 
    ## Run 349 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001170364  max resid 0.000325739 
    ## ... Similar to previous best
    ## Run 350 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001681734  max resid 0.0004429544 
    ## ... Similar to previous best
    ## Run 351 stress 0.06038962 
    ## Run 352 stress 0.365741 
    ## Run 353 stress 0.05036793 
    ## ... Procrustes: rmse 6.649545e-05  max resid 0.0001443795 
    ## ... Similar to previous best
    ## Run 354 stress 0.05051325 
    ## ... Procrustes: rmse 0.03576036  max resid 0.1220111 
    ## Run 355 stress 0.05171059 
    ## Run 356 stress 0.06640348 
    ## Run 357 stress 0.05036794 
    ## ... Procrustes: rmse 7.231136e-05  max resid 0.000194617 
    ## ... Similar to previous best
    ## Run 358 stress 0.05036807 
    ## ... Procrustes: rmse 0.0001740488  max resid 0.0004041995 
    ## ... Similar to previous best
    ## Run 359 stress 0.05466704 
    ## Run 360 stress 0.05403616 
    ## Run 361 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001484833  max resid 0.000358065 
    ## ... Similar to previous best
    ## Run 362 stress 0.06111768 
    ## Run 363 stress 0.06845265 
    ## Run 364 stress 0.05466702 
    ## Run 365 stress 0.0517106 
    ## Run 366 stress 0.05036796 
    ## ... Procrustes: rmse 9.818754e-05  max resid 0.0002591997 
    ## ... Similar to previous best
    ## Run 367 stress 0.06674546 
    ## Run 368 stress 0.06309848 
    ## Run 369 stress 0.06674553 
    ## Run 370 stress 0.06136696 
    ## Run 371 stress 0.05602061 
    ## Run 372 stress 0.05602073 
    ## Run 373 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001228516  max resid 0.0002829901 
    ## ... Similar to previous best
    ## Run 374 stress 0.05036792 
    ## ... Procrustes: rmse 3.079471e-05  max resid 8.625444e-05 
    ## ... Similar to previous best
    ## Run 375 stress 0.05036795 
    ## ... Procrustes: rmse 9.77342e-05  max resid 0.0002610011 
    ## ... Similar to previous best
    ## Run 376 stress 0.05036794 
    ## ... Procrustes: rmse 7.456838e-05  max resid 0.0001476036 
    ## ... Similar to previous best
    ## Run 377 stress 0.06231302 
    ## Run 378 stress 0.06499239 
    ## Run 379 stress 0.06478544 
    ## Run 380 stress 0.06943741 
    ## Run 381 stress 0.05829938 
    ## Run 382 stress 0.05829931 
    ## Run 383 stress 0.05900156 
    ## Run 384 stress 0.05466705 
    ## Run 385 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001103882  max resid 0.0002916785 
    ## ... Similar to previous best
    ## Run 386 stress 0.06904502 
    ## Run 387 stress 0.06896893 
    ## Run 388 stress 0.05602049 
    ## Run 389 stress 0.05036794 
    ## ... Procrustes: rmse 7.986845e-05  max resid 0.00021284 
    ## ... Similar to previous best
    ## Run 390 stress 0.05968626 
    ## Run 391 stress 0.07461379 
    ## Run 392 stress 0.05171059 
    ## Run 393 stress 0.05171058 
    ## Run 394 stress 0.05968633 
    ## Run 395 stress 0.05968636 
    ## Run 396 stress 0.06231291 
    ## Run 397 stress 0.06640348 
    ## Run 398 stress 0.06885085 
    ## Run 399 stress 0.05979229 
    ## Run 400 stress 0.06789768 
    ## Run 401 stress 0.05036797 
    ## ... Procrustes: rmse 8.274872e-05  max resid 0.0001865599 
    ## ... Similar to previous best
    ## Run 402 stress 0.05036794 
    ## ... Procrustes: rmse 8.029446e-05  max resid 0.0002076565 
    ## ... Similar to previous best
    ## Run 403 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001535599  max resid 0.0003999146 
    ## ... Similar to previous best
    ## Run 404 stress 0.05979234 
    ## Run 405 stress 0.05036792 
    ## ... Procrustes: rmse 3.72056e-05  max resid 8.092284e-05 
    ## ... Similar to previous best
    ## Run 406 stress 0.0596864 
    ## Run 407 stress 0.05051317 
    ## ... Procrustes: rmse 0.03568739  max resid 0.121937 
    ## Run 408 stress 0.05900156 
    ## Run 409 stress 0.0597923 
    ## Run 410 stress 0.06943848 
    ## Run 411 stress 0.05602057 
    ## Run 412 stress 0.05036791 
    ## ... Procrustes: rmse 1.790785e-05  max resid 3.624005e-05 
    ## ... Similar to previous best
    ## Run 413 stress 0.06038942 
    ## Run 414 stress 0.0592824 
    ## Run 415 stress 0.0592824 
    ## Run 416 stress 0.06027441 
    ## Run 417 stress 0.05337579 
    ## Run 418 stress 0.05051328 
    ## ... Procrustes: rmse 0.03572563  max resid 0.1220563 
    ## Run 419 stress 0.06499232 
    ## Run 420 stress 0.0517106 
    ## Run 421 stress 0.06136676 
    ## Run 422 stress 0.0611177 
    ## Run 423 stress 0.06630526 
    ## Run 424 stress 0.0582993 
    ## Run 425 stress 0.05968627 
    ## Run 426 stress 0.05036792 
    ## ... Procrustes: rmse 5.304308e-05  max resid 0.0001126277 
    ## ... Similar to previous best
    ## Run 427 stress 0.05900149 
    ## Run 428 stress 0.0647854 
    ## Run 429 stress 0.05051301 
    ## ... Procrustes: rmse 0.03568018  max resid 0.1217887 
    ## Run 430 stress 0.05051344 
    ## ... Procrustes: rmse 0.03572922  max resid 0.1220728 
    ## Run 431 stress 0.05051281 
    ## ... Procrustes: rmse 0.0355717  max resid 0.1215314 
    ## Run 432 stress 0.05900151 
    ## Run 433 stress 0.05466702 
    ## Run 434 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001244639  max resid 0.000346496 
    ## ... Similar to previous best
    ## Run 435 stress 0.05171058 
    ## Run 436 stress 0.06136676 
    ## Run 437 stress 0.06309855 
    ## Run 438 stress 0.05928224 
    ## Run 439 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001102028  max resid 0.0002962418 
    ## ... Similar to previous best
    ## Run 440 stress 0.05051302 
    ## ... Procrustes: rmse 0.03569385  max resid 0.1219491 
    ## Run 441 stress 0.06499234 
    ## Run 442 stress 0.05829932 
    ## Run 443 stress 0.05979228 
    ## Run 444 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001112742  max resid 0.0002913827 
    ## ... Similar to previous best
    ## Run 445 stress 0.0664036 
    ## Run 446 stress 0.05829944 
    ## Run 447 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001771706  max resid 0.0004692394 
    ## ... Similar to previous best
    ## Run 448 stress 0.05036794 
    ## ... Procrustes: rmse 8.904042e-05  max resid 0.0002369788 
    ## ... Similar to previous best
    ## Run 449 stress 0.0596864 
    ## Run 450 stress 0.05466702 
    ## Run 451 stress 0.05602054 
    ## Run 452 stress 0.05036795 
    ## ... Procrustes: rmse 9.014261e-05  max resid 0.0002412118 
    ## ... Similar to previous best
    ## Run 453 stress 0.054667 
    ## Run 454 stress 0.06938799 
    ## Run 455 stress 0.06111771 
    ## Run 456 stress 0.06872635 
    ## Run 457 stress 0.05036796 
    ## ... Procrustes: rmse 5.69247e-05  max resid 0.0001487985 
    ## ... Similar to previous best
    ## Run 458 stress 0.06476216 
    ## Run 459 stress 0.05051319 
    ## ... Procrustes: rmse 0.03572302  max resid 0.122045 
    ## Run 460 stress 0.06027428 
    ## Run 461 stress 0.06231291 
    ## Run 462 stress 0.06409295 
    ## Run 463 stress 0.05406615 
    ## Run 464 stress 0.05171059 
    ## Run 465 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001255994  max resid 0.0003321264 
    ## ... Similar to previous best
    ## Run 466 stress 0.05968636 
    ## Run 467 stress 0.05403586 
    ## Run 468 stress 0.05171058 
    ## Run 469 stress 0.06622012 
    ## Run 470 stress 0.06231309 
    ## Run 471 stress 0.07221109 
    ## Run 472 stress 0.06027427 
    ## Run 473 stress 0.05051307 
    ## ... Procrustes: rmse 0.0356709  max resid 0.1217538 
    ## Run 474 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001299498  max resid 0.0003413416 
    ## ... Similar to previous best
    ## Run 475 stress 0.05337578 
    ## Run 476 stress 0.05051301 
    ## ... Procrustes: rmse 0.03566458  max resid 0.1217432 
    ## Run 477 stress 0.050368 
    ## ... Procrustes: rmse 0.0001519183  max resid 0.0003989573 
    ## ... Similar to previous best
    ## Run 478 stress 0.05051293 
    ## ... Procrustes: rmse 0.03570217  max resid 0.1218682 
    ## Run 479 stress 0.06231301 
    ## Run 480 stress 0.05171062 
    ## Run 481 stress 0.05036792 
    ## ... Procrustes: rmse 4.357035e-05  max resid 9.569311e-05 
    ## ... Similar to previous best
    ## Run 482 stress 0.05051297 
    ## ... Procrustes: rmse 0.0356794  max resid 0.1219019 
    ## Run 483 stress 0.3594406 
    ## Run 484 stress 0.05036809 
    ## ... Procrustes: rmse 0.0002047812  max resid 0.0005223316 
    ## ... Similar to previous best
    ## Run 485 stress 0.05171058 
    ## Run 486 stress 0.05602049 
    ## Run 487 stress 0.05466709 
    ## Run 488 stress 0.05051334 
    ## ... Procrustes: rmse 0.03569982  max resid 0.1218134 
    ## Run 489 stress 0.05036794 
    ## ... Procrustes: rmse 7.810201e-05  max resid 0.0001980528 
    ## ... Similar to previous best
    ## Run 490 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001672378  max resid 0.0004266318 
    ## ... Similar to previous best
    ## Run 491 stress 0.050368 
    ## ... Procrustes: rmse 0.000127764  max resid 0.0003158614 
    ## ... Similar to previous best
    ## Run 492 stress 0.05403583 
    ## Run 493 stress 0.0582994 
    ## Run 494 stress 0.06231289 
    ## Run 495 stress 0.05171059 
    ## Run 496 stress 0.06111766 
    ## Run 497 stress 0.05337578 
    ## Run 498 stress 0.0517106 
    ## Run 499 stress 0.05051321 
    ## ... Procrustes: rmse 0.03565338  max resid 0.1216869 
    ## Run 500 stress 0.05466703 
    ## *** Best solution repeated 44 times

``` r
# Ocean sites and mixed lakes
PD_beta_OM_NMDS <- metaMDS(PD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1016168 
    ## Run 1 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0005919866  max resid 0.001657047 
    ## ... Similar to previous best
    ## Run 2 stress 0.349556 
    ## Run 3 stress 0.1341868 
    ## Run 4 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 2.743778e-05  max resid 7.715932e-05 
    ## ... Similar to previous best
    ## Run 5 stress 0.1341868 
    ## Run 6 stress 0.1016164 
    ## ... Procrustes: rmse 3.626639e-05  max resid 0.0001019064 
    ## ... Similar to previous best
    ## Run 7 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 9.525895e-06  max resid 2.673775e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.1530426 
    ## Run 9 stress 0.1016164 
    ## ... Procrustes: rmse 3.34086e-05  max resid 9.412231e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.1016164 
    ## ... Procrustes: rmse 9.273593e-06  max resid 2.591155e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.3147478 
    ## Run 12 stress 0.1016164 
    ## ... Procrustes: rmse 6.190352e-05  max resid 0.000174692 
    ## ... Similar to previous best
    ## Run 13 stress 0.1356144 
    ## Run 14 stress 0.1016164 
    ## ... Procrustes: rmse 5.60091e-06  max resid 1.574945e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 6.224777e-06  max resid 1.723838e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.1341868 
    ## Run 17 stress 0.1016164 
    ## ... Procrustes: rmse 6.452495e-05  max resid 0.0001816995 
    ## ... Similar to previous best
    ## Run 18 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004655055  max resid 0.001297499 
    ## ... Similar to previous best
    ## Run 19 stress 0.1016164 
    ## ... Procrustes: rmse 8.105114e-06  max resid 2.210008e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.1016164 
    ## ... Procrustes: rmse 2.054841e-05  max resid 5.760673e-05 
    ## ... Similar to previous best
    ## Run 21 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001123245  max resid 0.000316614 
    ## ... Similar to previous best
    ## Run 22 stress 0.1016164 
    ## ... Procrustes: rmse 9.916087e-05  max resid 0.0002778681 
    ## ... Similar to previous best
    ## Run 23 stress 0.1341868 
    ## Run 24 stress 0.1341868 
    ## Run 25 stress 0.1356147 
    ## Run 26 stress 0.3578479 
    ## Run 27 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006089425  max resid 0.001701843 
    ## ... Similar to previous best
    ## Run 28 stress 0.1341868 
    ## Run 29 stress 0.1016164 
    ## ... Procrustes: rmse 1.728507e-05  max resid 4.822165e-05 
    ## ... Similar to previous best
    ## Run 30 stress 0.1341868 
    ## Run 31 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005946874  max resid 0.001663133 
    ## ... Similar to previous best
    ## Run 32 stress 0.1016164 
    ## ... Procrustes: rmse 1.877259e-05  max resid 5.236171e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.1016164 
    ## ... Procrustes: rmse 7.873672e-05  max resid 0.0002215355 
    ## ... Similar to previous best
    ## Run 34 stress 0.1016164 
    ## ... Procrustes: rmse 4.859653e-05  max resid 0.0001361869 
    ## ... Similar to previous best
    ## Run 35 stress 0.1530423 
    ## Run 36 stress 0.1016164 
    ## ... Procrustes: rmse 5.975215e-05  max resid 0.0001686772 
    ## ... Similar to previous best
    ## Run 37 stress 0.1341868 
    ## Run 38 stress 0.1519033 
    ## Run 39 stress 0.1016164 
    ## ... Procrustes: rmse 4.347456e-05  max resid 0.0001226148 
    ## ... Similar to previous best
    ## Run 40 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001120243  max resid 0.0003162225 
    ## ... Similar to previous best
    ## Run 41 stress 0.1341868 
    ## Run 42 stress 0.1341868 
    ## Run 43 stress 0.1016164 
    ## ... Procrustes: rmse 5.201264e-05  max resid 0.0001464967 
    ## ... Similar to previous best
    ## Run 44 stress 0.1016164 
    ## ... Procrustes: rmse 8.404771e-05  max resid 0.0002358921 
    ## ... Similar to previous best
    ## Run 45 stress 0.1356143 
    ## Run 46 stress 0.1016164 
    ## ... Procrustes: rmse 8.624729e-05  max resid 0.0002430363 
    ## ... Similar to previous best
    ## Run 47 stress 0.1341868 
    ## Run 48 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001314474  max resid 0.0003695939 
    ## ... Similar to previous best
    ## Run 49 stress 0.1016164 
    ## ... Procrustes: rmse 4.547026e-06  max resid 1.267217e-05 
    ## ... Similar to previous best
    ## Run 50 stress 0.1016164 
    ## ... Procrustes: rmse 5.021643e-05  max resid 0.0001416374 
    ## ... Similar to previous best
    ## Run 51 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005557902  max resid 0.001553568 
    ## ... Similar to previous best
    ## Run 52 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004523328  max resid 0.001265506 
    ## ... Similar to previous best
    ## Run 53 stress 0.1530424 
    ## Run 54 stress 0.1356147 
    ## Run 55 stress 0.1016164 
    ## ... Procrustes: rmse 0.000125547  max resid 0.0003537647 
    ## ... Similar to previous best
    ## Run 56 stress 0.1016164 
    ## ... Procrustes: rmse 4.206937e-05  max resid 0.0001185168 
    ## ... Similar to previous best
    ## Run 57 stress 0.1016164 
    ## ... Procrustes: rmse 5.580746e-05  max resid 0.0001575786 
    ## ... Similar to previous best
    ## Run 58 stress 0.1016164 
    ## ... Procrustes: rmse 1.57082e-05  max resid 4.357219e-05 
    ## ... Similar to previous best
    ## Run 59 stress 0.1341868 
    ## Run 60 stress 0.1530425 
    ## Run 61 stress 0.1341868 
    ## Run 62 stress 0.1530422 
    ## Run 63 stress 0.1016164 
    ## ... Procrustes: rmse 2.23556e-05  max resid 6.296809e-05 
    ## ... Similar to previous best
    ## Run 64 stress 0.1016164 
    ## ... Procrustes: rmse 1.603891e-05  max resid 4.502617e-05 
    ## ... Similar to previous best
    ## Run 65 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003155246  max resid 0.0008837439 
    ## ... Similar to previous best
    ## Run 66 stress 0.1016164 
    ## ... Procrustes: rmse 1.834506e-05  max resid 5.08453e-05 
    ## ... Similar to previous best
    ## Run 67 stress 0.1356146 
    ## Run 68 stress 0.2109397 
    ## Run 69 stress 0.1341868 
    ## Run 70 stress 0.1356147 
    ## Run 71 stress 0.101617 
    ## ... Procrustes: rmse 0.0006594253  max resid 0.001843441 
    ## ... Similar to previous best
    ## Run 72 stress 0.1016164 
    ## ... Procrustes: rmse 2.011566e-05  max resid 5.497221e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.1016164 
    ## ... Procrustes: rmse 5.469799e-05  max resid 0.0001533018 
    ## ... Similar to previous best
    ## Run 74 stress 0.1356147 
    ## Run 75 stress 0.1016164 
    ## ... Procrustes: rmse 1.900481e-05  max resid 5.329247e-05 
    ## ... Similar to previous best
    ## Run 76 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001280249  max resid 0.0003569477 
    ## ... Similar to previous best
    ## Run 77 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004616024  max resid 0.001291336 
    ## ... Similar to previous best
    ## Run 78 stress 0.1016164 
    ## ... Procrustes: rmse 7.534404e-05  max resid 0.0002123227 
    ## ... Similar to previous best
    ## Run 79 stress 0.1016164 
    ## ... Procrustes: rmse 1.365335e-05  max resid 3.83641e-05 
    ## ... Similar to previous best
    ## Run 80 stress 0.3493327 
    ## Run 81 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001706789  max resid 0.0004798121 
    ## ... Similar to previous best
    ## Run 82 stress 0.1356143 
    ## Run 83 stress 0.1341868 
    ## Run 84 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 3.599674e-06  max resid 8.422108e-06 
    ## ... Similar to previous best
    ## Run 85 stress 0.1519033 
    ## Run 86 stress 0.1016164 
    ## ... Procrustes: rmse 2.515173e-05  max resid 7.151095e-05 
    ## ... Similar to previous best
    ## Run 87 stress 0.1341868 
    ## Run 88 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004770435  max resid 0.001333529 
    ## ... Similar to previous best
    ## Run 89 stress 0.1016164 
    ## ... Procrustes: rmse 4.244679e-05  max resid 0.0001202821 
    ## ... Similar to previous best
    ## Run 90 stress 0.1530428 
    ## Run 91 stress 0.1016164 
    ## ... Procrustes: rmse 1.277707e-05  max resid 3.521403e-05 
    ## ... Similar to previous best
    ## Run 92 stress 0.1016164 
    ## ... Procrustes: rmse 7.571872e-05  max resid 0.0002144842 
    ## ... Similar to previous best
    ## Run 93 stress 0.1016164 
    ## ... Procrustes: rmse 3.578089e-05  max resid 0.0001013711 
    ## ... Similar to previous best
    ## Run 94 stress 0.1356147 
    ## Run 95 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002106944  max resid 0.0005890038 
    ## ... Similar to previous best
    ## Run 96 stress 0.1016164 
    ## ... Procrustes: rmse 7.962594e-05  max resid 0.0002250729 
    ## ... Similar to previous best
    ## Run 97 stress 0.1016164 
    ## ... Procrustes: rmse 4.72057e-05  max resid 0.0001331537 
    ## ... Similar to previous best
    ## Run 98 stress 0.1341868 
    ## Run 99 stress 0.1341868 
    ## Run 100 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001054783  max resid 0.0002982588 
    ## ... Similar to previous best
    ## Run 101 stress 0.1016164 
    ## ... Procrustes: rmse 1.820621e-05  max resid 5.034432e-05 
    ## ... Similar to previous best
    ## Run 102 stress 0.1016164 
    ## ... Procrustes: rmse 4.06571e-05  max resid 0.0001137405 
    ## ... Similar to previous best
    ## Run 103 stress 0.1016164 
    ## ... Procrustes: rmse 2.299781e-06  max resid 5.188627e-06 
    ## ... Similar to previous best
    ## Run 104 stress 0.1016164 
    ## ... Procrustes: rmse 5.043158e-05  max resid 0.0001427738 
    ## ... Similar to previous best
    ## Run 105 stress 0.1016164 
    ## ... Procrustes: rmse 5.734174e-06  max resid 1.658918e-05 
    ## ... Similar to previous best
    ## Run 106 stress 0.336051 
    ## Run 107 stress 0.1016164 
    ## ... Procrustes: rmse 1.714005e-05  max resid 4.88678e-05 
    ## ... Similar to previous best
    ## Run 108 stress 0.1016164 
    ## ... Procrustes: rmse 1.120329e-05  max resid 2.61136e-05 
    ## ... Similar to previous best
    ## Run 109 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 1.333865e-06  max resid 2.293747e-06 
    ## ... Similar to previous best
    ## Run 110 stress 0.1341868 
    ## Run 111 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003755958  max resid 0.00105291 
    ## ... Similar to previous best
    ## Run 112 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002335206  max resid 0.0006545048 
    ## ... Similar to previous best
    ## Run 113 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004273081  max resid 0.001194654 
    ## ... Similar to previous best
    ## Run 114 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004359382  max resid 0.001217091 
    ## ... Similar to previous best
    ## Run 115 stress 0.1341868 
    ## Run 116 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002344107  max resid 0.0006568423 
    ## ... Similar to previous best
    ## Run 117 stress 0.320782 
    ## Run 118 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001593833  max resid 0.0004465564 
    ## ... Similar to previous best
    ## Run 119 stress 0.1341868 
    ## Run 120 stress 0.1341868 
    ## Run 121 stress 0.1016164 
    ## ... Procrustes: rmse 1.403371e-05  max resid 3.648595e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.1341868 
    ## Run 123 stress 0.1016164 
    ## ... Procrustes: rmse 1.261161e-05  max resid 2.530469e-05 
    ## ... Similar to previous best
    ## Run 124 stress 0.1016164 
    ## ... Procrustes: rmse 2.305973e-05  max resid 6.206449e-05 
    ## ... Similar to previous best
    ## Run 125 stress 0.1016164 
    ## ... Procrustes: rmse 5.730929e-05  max resid 0.0001613865 
    ## ... Similar to previous best
    ## Run 126 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004823537  max resid 0.001344851 
    ## ... Similar to previous best
    ## Run 127 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003598379  max resid 0.001007764 
    ## ... Similar to previous best
    ## Run 128 stress 0.1016164 
    ## ... Procrustes: rmse 6.401739e-05  max resid 0.000180952 
    ## ... Similar to previous best
    ## Run 129 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002339631  max resid 0.0006538045 
    ## ... Similar to previous best
    ## Run 130 stress 0.1356146 
    ## Run 131 stress 0.1016164 
    ## ... Procrustes: rmse 1.019461e-05  max resid 2.389742e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.1016164 
    ## ... Procrustes: rmse 2.130268e-05  max resid 6.046249e-05 
    ## ... Similar to previous best
    ## Run 133 stress 0.1341868 
    ## Run 134 stress 0.1016164 
    ## ... Procrustes: rmse 4.735123e-05  max resid 0.0001340406 
    ## ... Similar to previous best
    ## Run 135 stress 0.1016164 
    ## ... Procrustes: rmse 1.834298e-05  max resid 5.066961e-05 
    ## ... Similar to previous best
    ## Run 136 stress 0.1016164 
    ## ... Procrustes: rmse 4.978414e-06  max resid 1.243349e-05 
    ## ... Similar to previous best
    ## Run 137 stress 0.1341868 
    ## Run 138 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004098024  max resid 0.001146354 
    ## ... Similar to previous best
    ## Run 139 stress 0.1016164 
    ## ... Procrustes: rmse 8.578094e-05  max resid 0.0002418121 
    ## ... Similar to previous best
    ## Run 140 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001756825  max resid 0.0004911054 
    ## ... Similar to previous best
    ## Run 141 stress 0.1341868 
    ## Run 142 stress 0.1016164 
    ## ... Procrustes: rmse 6.83591e-05  max resid 0.0001881247 
    ## ... Similar to previous best
    ## Run 143 stress 0.1530429 
    ## Run 144 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005240597  max resid 0.001465407 
    ## ... Similar to previous best
    ## Run 145 stress 0.1016164 
    ## ... Procrustes: rmse 6.587301e-05  max resid 0.0001862122 
    ## ... Similar to previous best
    ## Run 146 stress 0.1356144 
    ## Run 147 stress 0.101617 
    ## ... Procrustes: rmse 0.0006427816  max resid 0.001792186 
    ## ... Similar to previous best
    ## Run 148 stress 0.1519033 
    ## Run 149 stress 0.1016164 
    ## ... Procrustes: rmse 2.570764e-05  max resid 7.156561e-05 
    ## ... Similar to previous best
    ## Run 150 stress 0.1016164 
    ## ... Procrustes: rmse 8.918502e-06  max resid 2.525331e-05 
    ## ... Similar to previous best
    ## Run 151 stress 0.1016164 
    ## ... Procrustes: rmse 1.80455e-05  max resid 5.068191e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.1016164 
    ## ... Procrustes: rmse 5.464499e-06  max resid 1.526751e-05 
    ## ... Similar to previous best
    ## Run 153 stress 0.1016164 
    ## ... Procrustes: rmse 8.134933e-05  max resid 0.0002295296 
    ## ... Similar to previous best
    ## Run 154 stress 0.1356148 
    ## Run 155 stress 0.3056532 
    ## Run 156 stress 0.1530428 
    ## Run 157 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005574884  max resid 0.001559514 
    ## ... Similar to previous best
    ## Run 158 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003056842  max resid 0.0008552531 
    ## ... Similar to previous best
    ## Run 159 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001547097  max resid 0.0004344098 
    ## ... Similar to previous best
    ## Run 160 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003504026  max resid 0.0009805544 
    ## ... Similar to previous best
    ## Run 161 stress 0.1016164 
    ## ... Procrustes: rmse 1.634004e-05  max resid 3.710563e-05 
    ## ... Similar to previous best
    ## Run 162 stress 0.1016164 
    ## ... Procrustes: rmse 4.472443e-05  max resid 0.0001266076 
    ## ... Similar to previous best
    ## Run 163 stress 0.1016164 
    ## ... Procrustes: rmse 1.228739e-05  max resid 3.391198e-05 
    ## ... Similar to previous best
    ## Run 164 stress 0.1341868 
    ## Run 165 stress 0.1016164 
    ## ... Procrustes: rmse 5.54357e-06  max resid 1.569045e-05 
    ## ... Similar to previous best
    ## Run 166 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002670102  max resid 0.0007479248 
    ## ... Similar to previous best
    ## Run 167 stress 0.1016164 
    ## ... Procrustes: rmse 3.925769e-05  max resid 7.028292e-05 
    ## ... Similar to previous best
    ## Run 168 stress 0.1356149 
    ## Run 169 stress 0.1016164 
    ## ... Procrustes: rmse 1.358716e-05  max resid 3.857121e-05 
    ## ... Similar to previous best
    ## Run 170 stress 0.1016164 
    ## ... Procrustes: rmse 8.726861e-06  max resid 2.475504e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.101617 
    ## ... Procrustes: rmse 0.0007090617  max resid 0.001980658 
    ## ... Similar to previous best
    ## Run 172 stress 0.1016164 
    ## ... Procrustes: rmse 3.864608e-05  max resid 0.0001093827 
    ## ... Similar to previous best
    ## Run 173 stress 0.1341868 
    ## Run 174 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004166197  max resid 0.001165861 
    ## ... Similar to previous best
    ## Run 175 stress 0.1016164 
    ## ... Procrustes: rmse 4.276116e-05  max resid 0.0001206696 
    ## ... Similar to previous best
    ## Run 176 stress 0.1016164 
    ## ... Procrustes: rmse 1.150427e-05  max resid 2.095034e-05 
    ## ... Similar to previous best
    ## Run 177 stress 0.1530429 
    ## Run 178 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004951289  max resid 0.001384029 
    ## ... Similar to previous best
    ## Run 179 stress 0.1356145 
    ## Run 180 stress 0.1016164 
    ## ... Procrustes: rmse 8.964757e-05  max resid 0.000251685 
    ## ... Similar to previous best
    ## Run 181 stress 0.1016168 
    ## ... Procrustes: rmse 0.000499449  max resid 0.00139656 
    ## ... Similar to previous best
    ## Run 182 stress 0.1016164 
    ## ... Procrustes: rmse 2.53551e-05  max resid 6.866232e-05 
    ## ... Similar to previous best
    ## Run 183 stress 0.1016164 
    ## ... Procrustes: rmse 8.50865e-05  max resid 0.0002402148 
    ## ... Similar to previous best
    ## Run 184 stress 0.1016164 
    ## ... Procrustes: rmse 8.378088e-05  max resid 0.0002368133 
    ## ... Similar to previous best
    ## Run 185 stress 0.1356151 
    ## Run 186 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 2.20571e-06  max resid 5.579326e-06 
    ## ... Similar to previous best
    ## Run 187 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005124104  max resid 0.001430711 
    ## ... Similar to previous best
    ## Run 188 stress 0.1016164 
    ## ... Procrustes: rmse 9.237143e-06  max resid 2.557241e-05 
    ## ... Similar to previous best
    ## Run 189 stress 0.1016164 
    ## ... Procrustes: rmse 1.348292e-05  max resid 3.762313e-05 
    ## ... Similar to previous best
    ## Run 190 stress 0.1341868 
    ## Run 191 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007011036  max resid 0.001957655 
    ## ... Similar to previous best
    ## Run 192 stress 0.1341868 
    ## Run 193 stress 0.1016164 
    ## ... Procrustes: rmse 1.859552e-05  max resid 5.223975e-05 
    ## ... Similar to previous best
    ## Run 194 stress 0.1530429 
    ## Run 195 stress 0.1016164 
    ## ... Procrustes: rmse 4.08119e-05  max resid 0.000113292 
    ## ... Similar to previous best
    ## Run 196 stress 0.1537555 
    ## Run 197 stress 0.1016164 
    ## ... Procrustes: rmse 2.195077e-05  max resid 6.158892e-05 
    ## ... Similar to previous best
    ## Run 198 stress 0.1530424 
    ## Run 199 stress 0.101617 
    ## ... Procrustes: rmse 0.000698667  max resid 0.001952062 
    ## ... Similar to previous best
    ## Run 200 stress 0.1016164 
    ## ... Procrustes: rmse 1.97769e-05  max resid 5.240465e-05 
    ## ... Similar to previous best
    ## Run 201 stress 0.1341868 
    ## Run 202 stress 0.1356149 
    ## Run 203 stress 0.1341868 
    ## Run 204 stress 0.1016164 
    ## ... Procrustes: rmse 3.18987e-06  max resid 8.466961e-06 
    ## ... Similar to previous best
    ## Run 205 stress 0.1016164 
    ## ... Procrustes: rmse 5.946994e-06  max resid 1.658626e-05 
    ## ... Similar to previous best
    ## Run 206 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004951317  max resid 0.001383541 
    ## ... Similar to previous best
    ## Run 207 stress 0.1016164 
    ## ... Procrustes: rmse 7.502112e-05  max resid 0.0002093501 
    ## ... Similar to previous best
    ## Run 208 stress 0.1016164 
    ## ... Procrustes: rmse 5.212161e-06  max resid 1.407782e-05 
    ## ... Similar to previous best
    ## Run 209 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003822687  max resid 0.001068092 
    ## ... Similar to previous best
    ## Run 210 stress 0.1016164 
    ## ... Procrustes: rmse 6.541926e-05  max resid 0.0001842811 
    ## ... Similar to previous best
    ## Run 211 stress 0.3578481 
    ## Run 212 stress 0.1530425 
    ## Run 213 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004136552  max resid 0.001158557 
    ## ... Similar to previous best
    ## Run 214 stress 0.1016164 
    ## ... Procrustes: rmse 4.979591e-05  max resid 0.0001377878 
    ## ... Similar to previous best
    ## Run 215 stress 0.1341868 
    ## Run 216 stress 0.101617 
    ## ... Procrustes: rmse 0.0006807323  max resid 0.0019021 
    ## ... Similar to previous best
    ## Run 217 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003587389  max resid 0.001009509 
    ## ... Similar to previous best
    ## Run 218 stress 0.1016164 
    ## ... Procrustes: rmse 1.903122e-05  max resid 5.312262e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.1016164 
    ## ... Procrustes: rmse 1.291013e-06  max resid 3.007533e-06 
    ## ... Similar to previous best
    ## Run 220 stress 0.1016164 
    ## ... Procrustes: rmse 1.369048e-05  max resid 3.812145e-05 
    ## ... Similar to previous best
    ## Run 221 stress 0.1341868 
    ## Run 222 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005205894  max resid 0.001455242 
    ## ... Similar to previous best
    ## Run 223 stress 0.1341868 
    ## Run 224 stress 0.1530423 
    ## Run 225 stress 0.1341868 
    ## Run 226 stress 0.1016164 
    ## ... Procrustes: rmse 3.853234e-06  max resid 1.009636e-05 
    ## ... Similar to previous best
    ## Run 227 stress 0.1356149 
    ## Run 228 stress 0.1016164 
    ## ... Procrustes: rmse 1.422103e-05  max resid 3.711581e-05 
    ## ... Similar to previous best
    ## Run 229 stress 0.1356143 
    ## Run 230 stress 0.1016164 
    ## ... Procrustes: rmse 4.914231e-05  max resid 0.0001368324 
    ## ... Similar to previous best
    ## Run 231 stress 0.1356149 
    ## Run 232 stress 0.101617 
    ## ... Procrustes: rmse 0.0006956501  max resid 0.001943603 
    ## ... Similar to previous best
    ## Run 233 stress 0.1356149 
    ## Run 234 stress 0.1016164 
    ## ... Procrustes: rmse 7.698076e-05  max resid 0.0002169943 
    ## ... Similar to previous best
    ## Run 235 stress 0.1016164 
    ## ... Procrustes: rmse 5.743924e-05  max resid 0.0001619262 
    ## ... Similar to previous best
    ## Run 236 stress 0.1341868 
    ## Run 237 stress 0.1356152 
    ## Run 238 stress 0.1341868 
    ## Run 239 stress 0.1341868 
    ## Run 240 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003541979  max resid 0.0009931343 
    ## ... Similar to previous best
    ## Run 241 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003027539  max resid 0.0008475827 
    ## ... Similar to previous best
    ## Run 242 stress 0.1016164 
    ## ... Procrustes: rmse 4.386896e-05  max resid 0.0001228879 
    ## ... Similar to previous best
    ## Run 243 stress 0.1356144 
    ## Run 244 stress 0.1356148 
    ## Run 245 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004536783  max resid 0.001270133 
    ## ... Similar to previous best
    ## Run 246 stress 0.1341868 
    ## Run 247 stress 0.1016164 
    ## ... Procrustes: rmse 2.086486e-06  max resid 5.911602e-06 
    ## ... Similar to previous best
    ## Run 248 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002985828  max resid 0.0008349924 
    ## ... Similar to previous best
    ## Run 249 stress 0.101617 
    ## ... Procrustes: rmse 0.0006223005  max resid 0.001736684 
    ## ... Similar to previous best
    ## Run 250 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006124134  max resid 0.001710588 
    ## ... Similar to previous best
    ## Run 251 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005603302  max resid 0.001566929 
    ## ... Similar to previous best
    ## Run 252 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004433061  max resid 0.001241005 
    ## ... Similar to previous best
    ## Run 253 stress 0.1356144 
    ## Run 254 stress 0.1016164 
    ## ... Procrustes: rmse 7.435818e-05  max resid 0.0002085827 
    ## ... Similar to previous best
    ## Run 255 stress 0.1016164 
    ## ... Procrustes: rmse 4.217489e-06  max resid 1.207575e-05 
    ## ... Similar to previous best
    ## Run 256 stress 0.1016167 
    ## ... Procrustes: rmse 0.000390286  max resid 0.001091685 
    ## ... Similar to previous best
    ## Run 257 stress 0.1356148 
    ## Run 258 stress 0.1016164 
    ## ... Procrustes: rmse 6.98004e-06  max resid 2.003441e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.1016164 
    ## ... Procrustes: rmse 5.009447e-05  max resid 0.0001405582 
    ## ... Similar to previous best
    ## Run 260 stress 0.1356145 
    ## Run 261 stress 0.1537555 
    ## Run 262 stress 0.1016167 
    ## ... Procrustes: rmse 0.000470442  max resid 0.001315808 
    ## ... Similar to previous best
    ## Run 263 stress 0.1016164 
    ## ... Procrustes: rmse 6.891616e-05  max resid 0.00019336 
    ## ... Similar to previous best
    ## Run 264 stress 0.1341868 
    ## Run 265 stress 0.1016164 
    ## ... Procrustes: rmse 6.18343e-05  max resid 0.0001740333 
    ## ... Similar to previous best
    ## Run 266 stress 0.1016164 
    ## ... Procrustes: rmse 2.80821e-06  max resid 7.411467e-06 
    ## ... Similar to previous best
    ## Run 267 stress 0.1519033 
    ## Run 268 stress 0.1341868 
    ## Run 269 stress 0.1356152 
    ## Run 270 stress 0.3397632 
    ## Run 271 stress 0.1016164 
    ## ... Procrustes: rmse 3.819034e-05  max resid 0.000106811 
    ## ... Similar to previous best
    ## Run 272 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001075584  max resid 0.0003028113 
    ## ... Similar to previous best
    ## Run 273 stress 0.1356152 
    ## Run 274 stress 0.1016164 
    ## ... Procrustes: rmse 2.98805e-05  max resid 8.377614e-05 
    ## ... Similar to previous best
    ## Run 275 stress 0.1530429 
    ## Run 276 stress 0.1341868 
    ## Run 277 stress 0.1016164 
    ## ... Procrustes: rmse 2.972674e-05  max resid 8.327472e-05 
    ## ... Similar to previous best
    ## Run 278 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001450719  max resid 0.0004080101 
    ## ... Similar to previous best
    ## Run 279 stress 0.1016164 
    ## ... Procrustes: rmse 5.724648e-05  max resid 0.0001608532 
    ## ... Similar to previous best
    ## Run 280 stress 0.1016164 
    ## ... Procrustes: rmse 5.55068e-06  max resid 1.524302e-05 
    ## ... Similar to previous best
    ## Run 281 stress 0.1356148 
    ## Run 282 stress 0.1341868 
    ## Run 283 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003374044  max resid 0.0009442737 
    ## ... Similar to previous best
    ## Run 284 stress 0.1016164 
    ## ... Procrustes: rmse 5.157392e-05  max resid 0.0001412111 
    ## ... Similar to previous best
    ## Run 285 stress 0.1016164 
    ## ... Procrustes: rmse 1.123022e-05  max resid 2.909191e-05 
    ## ... Similar to previous best
    ## Run 286 stress 0.1016164 
    ## ... Procrustes: rmse 1.604426e-05  max resid 4.478581e-05 
    ## ... Similar to previous best
    ## Run 287 stress 0.1341868 
    ## Run 288 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001458954  max resid 0.0004108849 
    ## ... Similar to previous best
    ## Run 289 stress 0.1016164 
    ## ... Procrustes: rmse 1.318331e-05  max resid 3.720274e-05 
    ## ... Similar to previous best
    ## Run 290 stress 0.1356152 
    ## Run 291 stress 0.1356145 
    ## Run 292 stress 0.1016164 
    ## ... Procrustes: rmse 6.095932e-06  max resid 1.674217e-05 
    ## ... Similar to previous best
    ## Run 293 stress 0.101617 
    ## ... Procrustes: rmse 0.0006855107  max resid 0.001915686 
    ## ... Similar to previous best
    ## Run 294 stress 0.1016164 
    ## ... Procrustes: rmse 3.257373e-05  max resid 9.091479e-05 
    ## ... Similar to previous best
    ## Run 295 stress 0.1016164 
    ## ... Procrustes: rmse 4.531102e-05  max resid 0.0001273584 
    ## ... Similar to previous best
    ## Run 296 stress 0.1341868 
    ## Run 297 stress 0.135615 
    ## Run 298 stress 0.1356144 
    ## Run 299 stress 0.1016164 
    ## ... Procrustes: rmse 4.675835e-05  max resid 0.0001311982 
    ## ... Similar to previous best
    ## Run 300 stress 0.1341868 
    ## Run 301 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005395037  max resid 0.001507743 
    ## ... Similar to previous best
    ## Run 302 stress 0.1356146 
    ## Run 303 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002225903  max resid 0.0006250097 
    ## ... Similar to previous best
    ## Run 304 stress 0.1016164 
    ## ... Procrustes: rmse 4.67813e-05  max resid 0.0001318117 
    ## ... Similar to previous best
    ## Run 305 stress 0.1016164 
    ## ... Procrustes: rmse 2.853321e-06  max resid 7.988891e-06 
    ## ... Similar to previous best
    ## Run 306 stress 0.1016164 
    ## ... Procrustes: rmse 6.469204e-06  max resid 1.759598e-05 
    ## ... Similar to previous best
    ## Run 307 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005204975  max resid 0.001454887 
    ## ... Similar to previous best
    ## Run 308 stress 0.1341868 
    ## Run 309 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001743193  max resid 0.0004882103 
    ## ... Similar to previous best
    ## Run 310 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004183058  max resid 0.001170486 
    ## ... Similar to previous best
    ## Run 311 stress 0.1016164 
    ## ... Procrustes: rmse 1.229285e-05  max resid 3.503881e-05 
    ## ... Similar to previous best
    ## Run 312 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006126754  max resid 0.001712322 
    ## ... Similar to previous best
    ## Run 313 stress 0.1341868 
    ## Run 314 stress 0.1016164 
    ## ... Procrustes: rmse 3.904254e-05  max resid 0.0001097105 
    ## ... Similar to previous best
    ## Run 315 stress 0.1016164 
    ## ... Procrustes: rmse 1.875184e-06  max resid 4.555199e-06 
    ## ... Similar to previous best
    ## Run 316 stress 0.1016164 
    ## ... Procrustes: rmse 1.616355e-05  max resid 4.445098e-05 
    ## ... Similar to previous best
    ## Run 317 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001845715  max resid 0.0005140302 
    ## ... Similar to previous best
    ## Run 318 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004081201  max resid 0.001142924 
    ## ... Similar to previous best
    ## Run 319 stress 0.1341868 
    ## Run 320 stress 0.1016164 
    ## ... Procrustes: rmse 2.052426e-05  max resid 4.683066e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.1016164 
    ## ... Procrustes: rmse 3.490626e-05  max resid 9.78861e-05 
    ## ... Similar to previous best
    ## Run 322 stress 0.1016164 
    ## ... Procrustes: rmse 4.560042e-06  max resid 1.244116e-05 
    ## ... Similar to previous best
    ## Run 323 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004018398  max resid 0.00112413 
    ## ... Similar to previous best
    ## Run 324 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005727162  max resid 0.001598581 
    ## ... Similar to previous best
    ## Run 325 stress 0.1016164 
    ## ... Procrustes: rmse 9.497861e-05  max resid 0.0002676238 
    ## ... Similar to previous best
    ## Run 326 stress 0.1356147 
    ## Run 327 stress 0.1016164 
    ## ... Procrustes: rmse 4.43193e-05  max resid 0.0001245245 
    ## ... Similar to previous best
    ## Run 328 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004422358  max resid 0.001237494 
    ## ... Similar to previous best
    ## Run 329 stress 0.1341868 
    ## Run 330 stress 0.1537555 
    ## Run 331 stress 0.1341868 
    ## Run 332 stress 0.1016164 
    ## ... Procrustes: rmse 1.80187e-06  max resid 3.884862e-06 
    ## ... Similar to previous best
    ## Run 333 stress 0.1016167 
    ## ... Procrustes: rmse 0.000353536  max resid 0.0009843909 
    ## ... Similar to previous best
    ## Run 334 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003891704  max resid 0.001088005 
    ## ... Similar to previous best
    ## Run 335 stress 0.1530424 
    ## Run 336 stress 0.2109396 
    ## Run 337 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004467031  max resid 0.001249903 
    ## ... Similar to previous best
    ## Run 338 stress 0.1016164 
    ## ... Procrustes: rmse 8.227193e-06  max resid 2.370423e-05 
    ## ... Similar to previous best
    ## Run 339 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002113998  max resid 0.0005885414 
    ## ... Similar to previous best
    ## Run 340 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005825824  max resid 0.00162785 
    ## ... Similar to previous best
    ## Run 341 stress 0.101617 
    ## ... Procrustes: rmse 0.0006192962  max resid 0.001726857 
    ## ... Similar to previous best
    ## Run 342 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001873165  max resid 0.0005267928 
    ## ... Similar to previous best
    ## Run 343 stress 0.1341868 
    ## Run 344 stress 0.1016171 
    ## ... Procrustes: rmse 0.0006526215  max resid 0.001816423 
    ## ... Similar to previous best
    ## Run 345 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003307633  max resid 0.0009237087 
    ## ... Similar to previous best
    ## Run 346 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004686681  max resid 0.001309564 
    ## ... Similar to previous best
    ## Run 347 stress 0.1016164 
    ## ... Procrustes: rmse 4.263526e-05  max resid 0.0001191185 
    ## ... Similar to previous best
    ## Run 348 stress 0.1341868 
    ## Run 349 stress 0.1016164 
    ## ... Procrustes: rmse 8.647845e-05  max resid 0.0002436617 
    ## ... Similar to previous best
    ## Run 350 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003936046  max resid 0.001100502 
    ## ... Similar to previous best
    ## Run 351 stress 0.1016164 
    ## ... Procrustes: rmse 3.988254e-05  max resid 0.000111946 
    ## ... Similar to previous best
    ## Run 352 stress 0.1016164 
    ## ... Procrustes: rmse 3.429988e-05  max resid 7.4905e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003370383  max resid 0.0009393581 
    ## ... Similar to previous best
    ## Run 354 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004453802  max resid 0.001245908 
    ## ... Similar to previous best
    ## Run 355 stress 0.1016164 
    ## ... Procrustes: rmse 6.41545e-05  max resid 0.0001794411 
    ## ... Similar to previous best
    ## Run 356 stress 0.1016164 
    ## ... Procrustes: rmse 4.723925e-06  max resid 1.063626e-05 
    ## ... Similar to previous best
    ## Run 357 stress 0.1530426 
    ## Run 358 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003749364  max resid 0.001045917 
    ## ... Similar to previous best
    ## Run 359 stress 0.101617 
    ## ... Procrustes: rmse 0.0006297299  max resid 0.001760619 
    ## ... Similar to previous best
    ## Run 360 stress 0.1016164 
    ## ... Procrustes: rmse 1.788789e-05  max resid 4.984877e-05 
    ## ... Similar to previous best
    ## Run 361 stress 0.1341868 
    ## Run 362 stress 0.1016164 
    ## ... Procrustes: rmse 9.67604e-06  max resid 2.673615e-05 
    ## ... Similar to previous best
    ## Run 363 stress 0.1341868 
    ## Run 364 stress 0.1016164 
    ## ... Procrustes: rmse 2.77419e-05  max resid 5.58338e-05 
    ## ... Similar to previous best
    ## Run 365 stress 0.1356152 
    ## Run 366 stress 0.1016164 
    ## ... Procrustes: rmse 4.513023e-06  max resid 1.192062e-05 
    ## ... Similar to previous best
    ## Run 367 stress 0.1356144 
    ## Run 368 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004222761  max resid 0.001181978 
    ## ... Similar to previous best
    ## Run 369 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003355095  max resid 0.0009364785 
    ## ... Similar to previous best
    ## Run 370 stress 0.101617 
    ## ... Procrustes: rmse 0.0006055055  max resid 0.001686237 
    ## ... Similar to previous best
    ## Run 371 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004621558  max resid 0.001292111 
    ## ... Similar to previous best
    ## Run 372 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004624327  max resid 0.001294519 
    ## ... Similar to previous best
    ## Run 373 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005883586  max resid 0.00164475 
    ## ... Similar to previous best
    ## Run 374 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005505017  max resid 0.001539218 
    ## ... Similar to previous best
    ## Run 375 stress 0.1016164 
    ## ... Procrustes: rmse 7.120904e-05  max resid 0.0002003931 
    ## ... Similar to previous best
    ## Run 376 stress 0.1016166 
    ## ... Procrustes: rmse 0.000379803  max resid 0.001062501 
    ## ... Similar to previous best
    ## Run 377 stress 0.1016171 
    ## ... Procrustes: rmse 0.0006365077  max resid 0.001780648 
    ## ... Similar to previous best
    ## Run 378 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001815472  max resid 0.0005089208 
    ## ... Similar to previous best
    ## Run 379 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004214245  max resid 0.001178581 
    ## ... Similar to previous best
    ## Run 380 stress 0.1016164 
    ## ... Procrustes: rmse 5.797684e-05  max resid 0.0001630493 
    ## ... Similar to previous best
    ## Run 381 stress 0.1530427 
    ## Run 382 stress 0.1016164 
    ## ... Procrustes: rmse 2.056279e-06  max resid 5.208012e-06 
    ## ... Similar to previous best
    ## Run 383 stress 0.1341868 
    ## Run 384 stress 0.1356152 
    ## Run 385 stress 0.1016164 
    ## ... Procrustes: rmse 3.355699e-05  max resid 9.445516e-05 
    ## ... Similar to previous best
    ## Run 386 stress 0.1016164 
    ## ... Procrustes: rmse 1.943784e-06  max resid 5.515313e-06 
    ## ... Similar to previous best
    ## Run 387 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001515288  max resid 0.0004263906 
    ## ... Similar to previous best
    ## Run 388 stress 0.1016164 
    ## ... Procrustes: rmse 1.759304e-05  max resid 4.725649e-05 
    ## ... Similar to previous best
    ## Run 389 stress 0.1016164 
    ## ... Procrustes: rmse 4.29975e-05  max resid 0.0001206386 
    ## ... Similar to previous best
    ## Run 390 stress 0.1341868 
    ## Run 391 stress 0.1016164 
    ## ... Procrustes: rmse 2.24318e-05  max resid 6.140062e-05 
    ## ... Similar to previous best
    ## Run 392 stress 0.1016164 
    ## ... Procrustes: rmse 2.640147e-05  max resid 7.423348e-05 
    ## ... Similar to previous best
    ## Run 393 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001343388  max resid 0.0003771336 
    ## ... Similar to previous best
    ## Run 394 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005198475  max resid 0.001454305 
    ## ... Similar to previous best
    ## Run 395 stress 0.1016164 
    ## ... Procrustes: rmse 7.486257e-05  max resid 0.0002110534 
    ## ... Similar to previous best
    ## Run 396 stress 0.1016164 
    ## ... Procrustes: rmse 8.654113e-05  max resid 0.0002424496 
    ## ... Similar to previous best
    ## Run 397 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001106107  max resid 0.000311494 
    ## ... Similar to previous best
    ## Run 398 stress 0.1356147 
    ## Run 399 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004761186  max resid 0.001328168 
    ## ... Similar to previous best
    ## Run 400 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004445093  max resid 0.001244371 
    ## ... Similar to previous best
    ## Run 401 stress 0.1530424 
    ## Run 402 stress 0.1356146 
    ## Run 403 stress 0.1356148 
    ## Run 404 stress 0.1341868 
    ## Run 405 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004133533  max resid 0.001157628 
    ## ... Similar to previous best
    ## Run 406 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005045374  max resid 0.001410786 
    ## ... Similar to previous best
    ## Run 407 stress 0.1341868 
    ## Run 408 stress 0.1016164 
    ## ... Procrustes: rmse 1.869356e-06  max resid 3.537261e-06 
    ## ... Similar to previous best
    ## Run 409 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003562357  max resid 0.000999183 
    ## ... Similar to previous best
    ## Run 410 stress 0.1519033 
    ## Run 411 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005834209  max resid 0.001632855 
    ## ... Similar to previous best
    ## Run 412 stress 0.1016164 
    ## ... Procrustes: rmse 3.854009e-05  max resid 0.0001078094 
    ## ... Similar to previous best
    ## Run 413 stress 0.1530428 
    ## Run 414 stress 0.1016168 
    ## ... Procrustes: rmse 0.000500487  max resid 0.001400664 
    ## ... Similar to previous best
    ## Run 415 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004676219  max resid 0.001308496 
    ## ... Similar to previous best
    ## Run 416 stress 0.1356145 
    ## Run 417 stress 0.1016164 
    ## ... Procrustes: rmse 4.615564e-05  max resid 0.0001289662 
    ## ... Similar to previous best
    ## Run 418 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002783721  max resid 0.0007782404 
    ## ... Similar to previous best
    ## Run 419 stress 0.1016164 
    ## ... Procrustes: rmse 5.720257e-05  max resid 0.0001609769 
    ## ... Similar to previous best
    ## Run 420 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001435705  max resid 0.0004031974 
    ## ... Similar to previous best
    ## Run 421 stress 0.1016164 
    ## ... Procrustes: rmse 4.81664e-05  max resid 0.0001356173 
    ## ... Similar to previous best
    ## Run 422 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003172495  max resid 0.0008873278 
    ## ... Similar to previous best
    ## Run 423 stress 0.1356152 
    ## Run 424 stress 0.1016169 
    ## ... Procrustes: rmse 0.000598548  max resid 0.001674217 
    ## ... Similar to previous best
    ## Run 425 stress 0.1016167 
    ## ... Procrustes: rmse 0.000387209  max resid 0.001082947 
    ## ... Similar to previous best
    ## Run 426 stress 0.1530428 
    ## Run 427 stress 0.1356144 
    ## Run 428 stress 0.1016164 
    ## ... Procrustes: rmse 6.950689e-05  max resid 0.0001941975 
    ## ... Similar to previous best
    ## Run 429 stress 0.1341868 
    ## Run 430 stress 0.1016164 
    ## ... Procrustes: rmse 8.810461e-05  max resid 0.0002480128 
    ## ... Similar to previous best
    ## Run 431 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005498213  max resid 0.001537078 
    ## ... Similar to previous best
    ## Run 432 stress 0.1016164 
    ## ... Procrustes: rmse 5.104493e-05  max resid 0.0001434844 
    ## ... Similar to previous best
    ## Run 433 stress 0.1016164 
    ## ... Procrustes: rmse 7.788884e-05  max resid 0.0002193216 
    ## ... Similar to previous best
    ## Run 434 stress 0.1016164 
    ## ... Procrustes: rmse 2.725488e-05  max resid 7.623416e-05 
    ## ... Similar to previous best
    ## Run 435 stress 0.1016164 
    ## ... Procrustes: rmse 1.412865e-05  max resid 3.937682e-05 
    ## ... Similar to previous best
    ## Run 436 stress 0.1016164 
    ## ... Procrustes: rmse 5.545542e-06  max resid 1.461757e-05 
    ## ... Similar to previous best
    ## Run 437 stress 0.1356149 
    ## Run 438 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003221792  max resid 0.0009023399 
    ## ... Similar to previous best
    ## Run 439 stress 0.1016164 
    ## ... Procrustes: rmse 2.136539e-05  max resid 6.060361e-05 
    ## ... Similar to previous best
    ## Run 440 stress 0.1016164 
    ## ... Procrustes: rmse 6.227346e-06  max resid 1.720196e-05 
    ## ... Similar to previous best
    ## Run 441 stress 0.1356146 
    ## Run 442 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004271501  max resid 0.001194848 
    ## ... Similar to previous best
    ## Run 443 stress 0.1341868 
    ## Run 444 stress 0.1341868 
    ## Run 445 stress 0.1341868 
    ## Run 446 stress 0.1356149 
    ## Run 447 stress 0.101617 
    ## ... Procrustes: rmse 0.0006527498  max resid 0.001825402 
    ## ... Similar to previous best
    ## Run 448 stress 0.1016164 
    ## ... Procrustes: rmse 6.011723e-05  max resid 0.0001694524 
    ## ... Similar to previous best
    ## Run 449 stress 0.1341868 
    ## Run 450 stress 0.1016164 
    ## ... Procrustes: rmse 4.162623e-05  max resid 0.0001167464 
    ## ... Similar to previous best
    ## Run 451 stress 0.1016164 
    ## ... Procrustes: rmse 1.192024e-05  max resid 2.337994e-05 
    ## ... Similar to previous best
    ## Run 452 stress 0.1356148 
    ## Run 453 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004967766  max resid 0.001388668 
    ## ... Similar to previous best
    ## Run 454 stress 0.101617 
    ## ... Procrustes: rmse 0.0006591877  max resid 0.00184117 
    ## ... Similar to previous best
    ## Run 455 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006308543  max resid 0.001763208 
    ## ... Similar to previous best
    ## Run 456 stress 0.1530426 
    ## Run 457 stress 0.1016164 
    ## ... Procrustes: rmse 7.783067e-06  max resid 2.161455e-05 
    ## ... Similar to previous best
    ## Run 458 stress 0.1016164 
    ## ... Procrustes: rmse 8.312776e-05  max resid 0.0002335328 
    ## ... Similar to previous best
    ## Run 459 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004500751  max resid 0.001259126 
    ## ... Similar to previous best
    ## Run 460 stress 0.1341868 
    ## Run 461 stress 0.1356147 
    ## Run 462 stress 0.1341868 
    ## Run 463 stress 0.1356146 
    ## Run 464 stress 0.1016164 
    ## ... Procrustes: rmse 9.055183e-06  max resid 2.508331e-05 
    ## ... Similar to previous best
    ## Run 465 stress 0.1356144 
    ## Run 466 stress 0.1016164 
    ## ... Procrustes: rmse 4.121188e-05  max resid 0.0001159058 
    ## ... Similar to previous best
    ## Run 467 stress 0.1341868 
    ## Run 468 stress 0.1016164 
    ## ... Procrustes: rmse 1.787427e-05  max resid 5.077911e-05 
    ## ... Similar to previous best
    ## Run 469 stress 0.1016164 
    ## ... Procrustes: rmse 4.347378e-05  max resid 0.000122158 
    ## ... Similar to previous best
    ## Run 470 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001692109  max resid 0.0004769524 
    ## ... Similar to previous best
    ## Run 471 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004628869  max resid 0.001293992 
    ## ... Similar to previous best
    ## Run 472 stress 0.1356149 
    ## Run 473 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003468456  max resid 0.0009713086 
    ## ... Similar to previous best
    ## Run 474 stress 0.1341868 
    ## Run 475 stress 0.1530424 
    ## Run 476 stress 0.1341868 
    ## Run 477 stress 0.1016164 
    ## ... Procrustes: rmse 3.760472e-06  max resid 1.057913e-05 
    ## ... Similar to previous best
    ## Run 478 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002315985  max resid 0.0006506707 
    ## ... Similar to previous best
    ## Run 479 stress 0.1016164 
    ## ... Procrustes: rmse 2.838624e-05  max resid 7.663675e-05 
    ## ... Similar to previous best
    ## Run 480 stress 0.1016164 
    ## ... Procrustes: rmse 0.000103658  max resid 0.0002914541 
    ## ... Similar to previous best
    ## Run 481 stress 0.1341868 
    ## Run 482 stress 0.1356148 
    ## Run 483 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003732504  max resid 0.001045523 
    ## ... Similar to previous best
    ## Run 484 stress 0.1341868 
    ## Run 485 stress 0.1016164 
    ## ... Procrustes: rmse 5.858389e-05  max resid 0.0001645569 
    ## ... Similar to previous best
    ## Run 486 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003182272  max resid 0.000890216 
    ## ... Similar to previous best
    ## Run 487 stress 0.1341868 
    ## Run 488 stress 0.1356145 
    ## Run 489 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003702733  max resid 0.001036801 
    ## ... Similar to previous best
    ## Run 490 stress 0.1530425 
    ## Run 491 stress 0.1341868 
    ## Run 492 stress 0.1356145 
    ## Run 493 stress 0.1016164 
    ## ... Procrustes: rmse 1.522894e-05  max resid 4.203433e-05 
    ## ... Similar to previous best
    ## Run 494 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004539545  max resid 0.001267989 
    ## ... Similar to previous best
    ## Run 495 stress 0.1341868 
    ## Run 496 stress 0.1016164 
    ## ... Procrustes: rmse 5.76035e-05  max resid 0.0001621256 
    ## ... Similar to previous best
    ## Run 497 stress 0.1016164 
    ## ... Procrustes: rmse 3.685619e-06  max resid 9.628635e-06 
    ## ... Similar to previous best
    ## Run 498 stress 0.1356147 
    ## Run 499 stress 0.1016164 
    ## ... Procrustes: rmse 8.028646e-06  max resid 2.233043e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003638197  max resid 0.001017344 
    ## ... Similar to previous best
    ## *** Best solution repeated 205 times

``` r
# Stratified lakes and ocean sites
PD_beta_SO_NMDS <- metaMDS(PD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.04851417 
    ## Run 1 stress 0.05728088 
    ## Run 2 stress 0.05612251 
    ## Run 3 stress 0.04938617 
    ## Run 4 stress 0.05607344 
    ## Run 5 stress 0.05695888 
    ## Run 6 stress 0.05605151 
    ## Run 7 stress 0.05643261 
    ## Run 8 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04980871  max resid 0.1394444 
    ## Run 9 stress 0.05216146 
    ## Run 10 stress 0.05695857 
    ## Run 11 stress 0.05438946 
    ## Run 12 stress 0.06185052 
    ## Run 13 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001513324  max resid 0.0003051093 
    ## ... Similar to previous best
    ## Run 14 stress 0.04938617 
    ## Run 15 stress 0.05695855 
    ## Run 16 stress 0.06185047 
    ## Run 17 stress 0.06400754 
    ## Run 18 stress 0.05842801 
    ## Run 19 stress 0.0704634 
    ## Run 20 stress 0.05728091 
    ## Run 21 stress 0.05451493 
    ## Run 22 stress 0.04860422 
    ## Run 23 stress 0.05607345 
    ## Run 24 stress 0.05605152 
    ## Run 25 stress 0.0585503 
    ## Run 26 stress 0.05605147 
    ## Run 27 stress 0.05607346 
    ## Run 28 stress 0.04851418 
    ## Run 29 stress 0.05438955 
    ## Run 30 stress 0.05728062 
    ## Run 31 stress 0.05855033 
    ## Run 32 stress 0.0493862 
    ## Run 33 stress 0.05695856 
    ## Run 34 stress 0.04938616 
    ## Run 35 stress 0.06291425 
    ## Run 36 stress 0.05488105 
    ## Run 37 stress 0.06291425 
    ## Run 38 stress 0.04938615 
    ## Run 39 stress 0.04938616 
    ## Run 40 stress 0.05474816 
    ## Run 41 stress 0.05260047 
    ## Run 42 stress 0.04935944 
    ## Run 43 stress 0.04851418 
    ## Run 44 stress 0.05959899 
    ## Run 45 stress 0.06205755 
    ## Run 46 stress 0.05276412 
    ## Run 47 stress 0.0527641 
    ## Run 48 stress 0.05216135 
    ## Run 49 stress 0.05451494 
    ## Run 50 stress 0.05695882 
    ## Run 51 stress 0.04851418 
    ## Run 52 stress 0.05728078 
    ## Run 53 stress 0.0526006 
    ## Run 54 stress 0.05260056 
    ## Run 55 stress 0.0560515 
    ## Run 56 stress 0.05728065 
    ## Run 57 stress 0.05216158 
    ## Run 58 stress 0.05474817 
    ## Run 59 stress 0.05433176 
    ## Run 60 stress 0.0608898 
    ## Run 61 stress 0.04612759 
    ## ... Procrustes: rmse 9.290607e-05  max resid 0.0001880687 
    ## ... Similar to previous best
    ## Run 62 stress 0.05612257 
    ## Run 63 stress 0.06291434 
    ## Run 64 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001080242  max resid 0.0002102222 
    ## ... Similar to previous best
    ## Run 65 stress 0.04938615 
    ## Run 66 stress 0.05216128 
    ## Run 67 stress 0.04938621 
    ## Run 68 stress 0.05260059 
    ## Run 69 stress 0.05695867 
    ## Run 70 stress 0.05605156 
    ## Run 71 stress 0.05605151 
    ## Run 72 stress 0.0629143 
    ## Run 73 stress 0.05572688 
    ## Run 74 stress 0.05728062 
    ## Run 75 stress 0.05593463 
    ## Run 76 stress 0.05260047 
    ## Run 77 stress 0.05137965 
    ## Run 78 stress 0.05959895 
    ## Run 79 stress 0.05593459 
    ## Run 80 stress 0.0543894 
    ## Run 81 stress 0.0543928 
    ## Run 82 stress 0.0561224 
    ## Run 83 stress 0.05451493 
    ## Run 84 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001138528  max resid 0.0002297924 
    ## ... Similar to previous best
    ## Run 85 stress 0.04938615 
    ## Run 86 stress 0.05728055 
    ## Run 87 stress 0.05855032 
    ## Run 88 stress 0.05728071 
    ## Run 89 stress 0.05260056 
    ## Run 90 stress 0.05607345 
    ## Run 91 stress 0.04851417 
    ## Run 92 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001557921  max resid 0.0003112403 
    ## ... Similar to previous best
    ## Run 93 stress 0.06291431 
    ## Run 94 stress 0.05474814 
    ## Run 95 stress 0.05695876 
    ## Run 96 stress 0.0527641 
    ## Run 97 stress 0.05474814 
    ## Run 98 stress 0.04860422 
    ## Run 99 stress 0.04938617 
    ## Run 100 stress 0.05216149 
    ## Run 101 stress 0.06185055 
    ## Run 102 stress 0.05474814 
    ## Run 103 stress 0.06185047 
    ## Run 104 stress 0.05488112 
    ## Run 105 stress 0.05607345 
    ## Run 106 stress 0.05474815 
    ## Run 107 stress 0.05593452 
    ## Run 108 stress 0.05260048 
    ## Run 109 stress 0.05474816 
    ## Run 110 stress 0.05607345 
    ## Run 111 stress 0.05605155 
    ## Run 112 stress 0.05607344 
    ## Run 113 stress 0.05216143 
    ## Run 114 stress 0.05593466 
    ## Run 115 stress 0.04612757 
    ## ... Procrustes: rmse 1.850379e-05  max resid 3.005215e-05 
    ## ... Similar to previous best
    ## Run 116 stress 0.04851417 
    ## Run 117 stress 0.05728065 
    ## Run 118 stress 0.04851418 
    ## Run 119 stress 0.05137963 
    ## Run 120 stress 0.05855032 
    ## Run 121 stress 0.05451493 
    ## Run 122 stress 0.05474817 
    ## Run 123 stress 0.0547482 
    ## Run 124 stress 0.04935945 
    ## Run 125 stress 0.04938615 
    ## Run 126 stress 0.04935947 
    ## Run 127 stress 0.06400738 
    ## Run 128 stress 0.04851419 
    ## Run 129 stress 0.04935943 
    ## Run 130 stress 0.05607344 
    ## Run 131 stress 0.05260052 
    ## Run 132 stress 0.05607345 
    ## Run 133 stress 0.05612244 
    ## Run 134 stress 0.05451494 
    ## Run 135 stress 0.05593457 
    ## Run 136 stress 0.04851418 
    ## Run 137 stress 0.05607348 
    ## Run 138 stress 0.05842797 
    ## Run 139 stress 0.05728064 
    ## Run 140 stress 0.05474813 
    ## Run 141 stress 0.0585503 
    ## Run 142 stress 0.05605147 
    ## Run 143 stress 0.05216137 
    ## Run 144 stress 0.05474814 
    ## Run 145 stress 0.05643261 
    ## Run 146 stress 0.05260051 
    ## Run 147 stress 0.0561225 
    ## Run 148 stress 0.05474818 
    ## Run 149 stress 0.05959898 
    ## Run 150 stress 0.05607347 
    ## Run 151 stress 0.05959896 
    ## Run 152 stress 0.05612243 
    ## Run 153 stress 0.05855035 
    ## Run 154 stress 0.04938615 
    ## Run 155 stress 0.04851417 
    ## Run 156 stress 0.05488107 
    ## Run 157 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001062718  max resid 0.0002134408 
    ## ... Similar to previous best
    ## Run 158 stress 0.05695883 
    ## Run 159 stress 0.0485142 
    ## Run 160 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 1.797342e-05  max resid 2.866014e-05 
    ## ... Similar to previous best
    ## Run 161 stress 0.05643265 
    ## Run 162 stress 0.05216133 
    ## Run 163 stress 0.0560735 
    ## Run 164 stress 0.05728057 
    ## Run 165 stress 0.05728052 
    ## Run 166 stress 0.0461276 
    ## ... Procrustes: rmse 8.415418e-05  max resid 0.000169775 
    ## ... Similar to previous best
    ## Run 167 stress 0.04935956 
    ## Run 168 stress 0.04851417 
    ## Run 169 stress 0.06205752 
    ## Run 170 stress 0.05605151 
    ## Run 171 stress 0.04938616 
    ## Run 172 stress 0.05438939 
    ## Run 173 stress 0.05137964 
    ## Run 174 stress 0.05276411 
    ## Run 175 stress 0.04935944 
    ## Run 176 stress 0.0560515 
    ## Run 177 stress 0.05451494 
    ## Run 178 stress 0.05728056 
    ## Run 179 stress 0.05216152 
    ## Run 180 stress 0.05959906 
    ## Run 181 stress 0.04851417 
    ## Run 182 stress 0.05572689 
    ## Run 183 stress 0.04938617 
    ## Run 184 stress 0.05451498 
    ## Run 185 stress 0.05855029 
    ## Run 186 stress 0.05137964 
    ## Run 187 stress 0.06205754 
    ## Run 188 stress 0.05216141 
    ## Run 189 stress 0.05572689 
    ## Run 190 stress 0.05451494 
    ## Run 191 stress 0.05607344 
    ## Run 192 stress 0.05488106 
    ## Run 193 stress 0.05605151 
    ## Run 194 stress 0.05433172 
    ## Run 195 stress 0.05607345 
    ## Run 196 stress 0.05607344 
    ## Run 197 stress 0.04851418 
    ## Run 198 stress 0.0461276 
    ## ... Procrustes: rmse 7.436707e-05  max resid 0.0001516457 
    ## ... Similar to previous best
    ## Run 199 stress 0.05474817 
    ## Run 200 stress 0.0493862 
    ## Run 201 stress 0.05959897 
    ## Run 202 stress 0.05643262 
    ## Run 203 stress 0.04851417 
    ## Run 204 stress 0.05607344 
    ## Run 205 stress 0.05728063 
    ## Run 206 stress 0.0521613 
    ## Run 207 stress 0.0569588 
    ## Run 208 stress 0.04851418 
    ## Run 209 stress 0.05572688 
    ## Run 210 stress 0.04851418 
    ## Run 211 stress 0.05137966 
    ## Run 212 stress 0.05728076 
    ## Run 213 stress 0.05605151 
    ## Run 214 stress 0.05488103 
    ## Run 215 stress 0.05959902 
    ## Run 216 stress 0.05695869 
    ## Run 217 stress 0.05433176 
    ## Run 218 stress 0.04851417 
    ## Run 219 stress 0.06205753 
    ## Run 220 stress 0.04938616 
    ## Run 221 stress 0.05605153 
    ## Run 222 stress 0.05605148 
    ## Run 223 stress 0.04938616 
    ## Run 224 stress 0.05474817 
    ## Run 225 stress 0.04935942 
    ## Run 226 stress 0.05607344 
    ## Run 227 stress 0.05855042 
    ## Run 228 stress 0.3584402 
    ## Run 229 stress 0.05451499 
    ## Run 230 stress 0.04860423 
    ## Run 231 stress 0.05643261 
    ## Run 232 stress 0.05276411 
    ## Run 233 stress 0.05612253 
    ## Run 234 stress 0.05728053 
    ## Run 235 stress 0.05695863 
    ## Run 236 stress 0.04938615 
    ## Run 237 stress 0.05488104 
    ## Run 238 stress 0.05728074 
    ## Run 239 stress 0.04851417 
    ## Run 240 stress 0.04851417 
    ## Run 241 stress 0.05607346 
    ## Run 242 stress 0.04612767 
    ## ... Procrustes: rmse 0.000152296  max resid 0.0003137266 
    ## ... Similar to previous best
    ## Run 243 stress 0.05643264 
    ## Run 244 stress 0.04851417 
    ## Run 245 stress 0.05474822 
    ## Run 246 stress 0.04860423 
    ## Run 247 stress 0.04935954 
    ## Run 248 stress 0.0560515 
    ## Run 249 stress 0.05593455 
    ## Run 250 stress 0.0585503 
    ## Run 251 stress 0.05593452 
    ## Run 252 stress 0.04851417 
    ## Run 253 stress 0.0543928 
    ## Run 254 stress 0.05612243 
    ## Run 255 stress 0.04860424 
    ## Run 256 stress 0.05438939 
    ## Run 257 stress 0.04938615 
    ## Run 258 stress 0.04612758 
    ## ... Procrustes: rmse 4.218418e-05  max resid 8.572952e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.06205754 
    ## Run 260 stress 0.05607345 
    ## Run 261 stress 0.04851417 
    ## Run 262 stress 0.05572688 
    ## Run 263 stress 0.05276411 
    ## Run 264 stress 0.04851417 
    ## Run 265 stress 0.05728055 
    ## Run 266 stress 0.05593456 
    ## Run 267 stress 0.05612255 
    ## Run 268 stress 0.0543928 
    ## Run 269 stress 0.05607344 
    ## Run 270 stress 0.04938616 
    ## Run 271 stress 0.0585503 
    ## Run 272 stress 0.06185052 
    ## Run 273 stress 0.05216129 
    ## Run 274 stress 0.04612761 
    ## ... Procrustes: rmse 9.5029e-05  max resid 0.0001919689 
    ## ... Similar to previous best
    ## Run 275 stress 0.05593459 
    ## Run 276 stress 0.04938616 
    ## Run 277 stress 0.04612761 
    ## ... Procrustes: rmse 8.792468e-05  max resid 0.0001751042 
    ## ... Similar to previous best
    ## Run 278 stress 0.06291428 
    ## Run 279 stress 0.05438945 
    ## Run 280 stress 0.0461276 
    ## ... Procrustes: rmse 7.394548e-05  max resid 0.0001502113 
    ## ... Similar to previous best
    ## Run 281 stress 0.05612239 
    ## Run 282 stress 0.05643262 
    ## Run 283 stress 0.05695864 
    ## Run 284 stress 0.05137968 
    ## Run 285 stress 0.04851417 
    ## Run 286 stress 0.04938616 
    ## Run 287 stress 0.04938617 
    ## Run 288 stress 0.05612239 
    ## Run 289 stress 0.05728063 
    ## Run 290 stress 0.3578481 
    ## Run 291 stress 0.05612246 
    ## Run 292 stress 0.05959897 
    ## Run 293 stress 0.04938616 
    ## Run 294 stress 0.04938616 
    ## Run 295 stress 0.04851417 
    ## Run 296 stress 0.05216129 
    ## Run 297 stress 0.05728052 
    ## Run 298 stress 0.04851418 
    ## Run 299 stress 0.05695873 
    ## Run 300 stress 0.05605151 
    ## Run 301 stress 0.05438945 
    ## Run 302 stress 0.04851417 
    ## Run 303 stress 0.05474813 
    ## Run 304 stress 0.05433175 
    ## Run 305 stress 0.0608898 
    ## Run 306 stress 0.05607344 
    ## Run 307 stress 0.05607344 
    ## Run 308 stress 0.05216136 
    ## Run 309 stress 0.06088975 
    ## Run 310 stress 0.06205753 
    ## Run 311 stress 0.05728058 
    ## Run 312 stress 0.05137969 
    ## Run 313 stress 0.05959896 
    ## Run 314 stress 0.04860424 
    ## Run 315 stress 0.05643262 
    ## Run 316 stress 0.05451495 
    ## Run 317 stress 0.05605153 
    ## Run 318 stress 0.05855041 
    ## Run 319 stress 0.05855033 
    ## Run 320 stress 0.05216138 
    ## Run 321 stress 0.05474815 
    ## Run 322 stress 0.05959895 
    ## Run 323 stress 0.04938616 
    ## Run 324 stress 0.05605148 
    ## Run 325 stress 0.05728056 
    ## Run 326 stress 0.05439287 
    ## Run 327 stress 0.05842785 
    ## Run 328 stress 0.05572691 
    ## Run 329 stress 0.06400735 
    ## Run 330 stress 0.0543929 
    ## Run 331 stress 0.05842794 
    ## Run 332 stress 0.0461276 
    ## ... Procrustes: rmse 8.386931e-05  max resid 0.0001707015 
    ## ... Similar to previous best
    ## Run 333 stress 0.0561224 
    ## Run 334 stress 0.04938616 
    ## Run 335 stress 0.06205754 
    ## Run 336 stress 0.05695888 
    ## Run 337 stress 0.05260049 
    ## Run 338 stress 0.05216157 
    ## Run 339 stress 0.05728057 
    ## Run 340 stress 0.05605145 
    ## Run 341 stress 0.05728053 
    ## Run 342 stress 0.05607345 
    ## Run 343 stress 0.04612762 
    ## ... Procrustes: rmse 0.0001263321  max resid 0.0002559677 
    ## ... Similar to previous best
    ## Run 344 stress 0.04851417 
    ## Run 345 stress 0.05607346 
    ## Run 346 stress 0.05605145 
    ## Run 347 stress 0.05572688 
    ## Run 348 stress 0.05451494 
    ## Run 349 stress 0.06400727 
    ## Run 350 stress 0.05474817 
    ## Run 351 stress 0.04938615 
    ## Run 352 stress 0.05451493 
    ## Run 353 stress 0.05695857 
    ## Run 354 stress 0.0543928 
    ## Run 355 stress 0.06205758 
    ## Run 356 stress 0.05695869 
    ## Run 357 stress 0.05216141 
    ## Run 358 stress 0.05474819 
    ## Run 359 stress 0.04851417 
    ## Run 360 stress 0.0608897 
    ## Run 361 stress 0.05593461 
    ## Run 362 stress 0.05216166 
    ## Run 363 stress 0.05607345 
    ## Run 364 stress 0.05728068 
    ## Run 365 stress 0.05216133 
    ## Run 366 stress 0.0461276 
    ## ... Procrustes: rmse 5.957477e-05  max resid 0.0001148574 
    ## ... Similar to previous best
    ## Run 367 stress 0.04860422 
    ## Run 368 stress 0.05643263 
    ## Run 369 stress 0.05438954 
    ## Run 370 stress 0.05438937 
    ## Run 371 stress 0.05959899 
    ## Run 372 stress 0.04938616 
    ## Run 373 stress 0.05276412 
    ## Run 374 stress 0.05488103 
    ## Run 375 stress 0.05607347 
    ## Run 376 stress 0.04938618 
    ## Run 377 stress 0.05643261 
    ## Run 378 stress 0.04851417 
    ## Run 379 stress 0.0559345 
    ## Run 380 stress 0.05643261 
    ## Run 381 stress 0.05216156 
    ## Run 382 stress 0.06205753 
    ## Run 383 stress 0.04938617 
    ## Run 384 stress 0.04851418 
    ## Run 385 stress 0.05488103 
    ## Run 386 stress 0.05447825 
    ## Run 387 stress 0.05260051 
    ## Run 388 stress 0.05439283 
    ## Run 389 stress 0.05643261 
    ## Run 390 stress 0.04938616 
    ## Run 391 stress 0.05728052 
    ## Run 392 stress 0.04851419 
    ## Run 393 stress 0.05260047 
    ## Run 394 stress 0.05276414 
    ## Run 395 stress 0.06088976 
    ## Run 396 stress 0.05643262 
    ## Run 397 stress 0.05842789 
    ## Run 398 stress 0.05474814 
    ## Run 399 stress 0.05216125 
    ## Run 400 stress 0.05695869 
    ## Run 401 stress 0.04938618 
    ## Run 402 stress 0.0527641 
    ## Run 403 stress 0.0585503 
    ## Run 404 stress 0.0526005 
    ## Run 405 stress 0.05276413 
    ## Run 406 stress 0.05137964 
    ## Run 407 stress 0.0527641 
    ## Run 408 stress 0.04938616 
    ## Run 409 stress 0.06291425 
    ## Run 410 stress 0.05439287 
    ## Run 411 stress 0.05438945 
    ## Run 412 stress 0.05643261 
    ## Run 413 stress 0.04612757 
    ## ... Procrustes: rmse 1.2421e-05  max resid 2.152309e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.05474817 
    ## Run 415 stress 0.05260052 
    ## Run 416 stress 0.0585503 
    ## Run 417 stress 0.05607351 
    ## Run 418 stress 0.0543928 
    ## Run 419 stress 0.05643261 
    ## Run 420 stress 0.04851417 
    ## Run 421 stress 0.05643261 
    ## Run 422 stress 0.04938619 
    ## Run 423 stress 0.05607344 
    ## Run 424 stress 0.0547482 
    ## Run 425 stress 0.05593453 
    ## Run 426 stress 0.05137964 
    ## Run 427 stress 0.05605149 
    ## Run 428 stress 0.05643263 
    ## Run 429 stress 0.04935948 
    ## Run 430 stress 0.04935943 
    ## Run 431 stress 0.0543895 
    ## Run 432 stress 0.05572688 
    ## Run 433 stress 0.05451494 
    ## Run 434 stress 0.06291428 
    ## Run 435 stress 0.05612252 
    ## Run 436 stress 0.05451498 
    ## Run 437 stress 0.04938616 
    ## Run 438 stress 0.04935947 
    ## Run 439 stress 0.04612766 
    ## ... Procrustes: rmse 0.0001337699  max resid 0.0002730811 
    ## ... Similar to previous best
    ## Run 440 stress 0.04938616 
    ## Run 441 stress 0.04938616 
    ## Run 442 stress 0.05572689 
    ## Run 443 stress 0.06291426 
    ## Run 444 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001318814  max resid 0.0002649757 
    ## ... Similar to previous best
    ## Run 445 stress 0.04860422 
    ## Run 446 stress 0.05842788 
    ## Run 447 stress 0.3229249 
    ## Run 448 stress 0.04860422 
    ## Run 449 stress 0.06205756 
    ## Run 450 stress 0.05607348 
    ## Run 451 stress 0.05855034 
    ## Run 452 stress 0.04851417 
    ## Run 453 stress 0.05438942 
    ## Run 454 stress 0.04938616 
    ## Run 455 stress 0.04612758 
    ## ... Procrustes: rmse 5.410066e-05  max resid 0.0001112449 
    ## ... Similar to previous best
    ## Run 456 stress 0.06291425 
    ## Run 457 stress 0.05260057 
    ## Run 458 stress 0.0560515 
    ## Run 459 stress 0.05728074 
    ## Run 460 stress 0.04935942 
    ## Run 461 stress 0.05572688 
    ## Run 462 stress 0.06291427 
    ## Run 463 stress 0.05643262 
    ## Run 464 stress 0.04860423 
    ## Run 465 stress 0.04935948 
    ## Run 466 stress 0.05447811 
    ## Run 467 stress 0.04851418 
    ## Run 468 stress 0.05731516 
    ## Run 469 stress 0.04851417 
    ## Run 470 stress 0.05607344 
    ## Run 471 stress 0.05439284 
    ## Run 472 stress 0.06205753 
    ## Run 473 stress 0.04935945 
    ## Run 474 stress 0.06291428 
    ## Run 475 stress 0.04612757 
    ## ... Procrustes: rmse 1.505438e-05  max resid 3.028185e-05 
    ## ... Similar to previous best
    ## Run 476 stress 0.05695869 
    ## Run 477 stress 0.05137963 
    ## Run 478 stress 0.3495556 
    ## Run 479 stress 0.0572806 
    ## Run 480 stress 0.04851417 
    ## Run 481 stress 0.04938616 
    ## Run 482 stress 0.05695873 
    ## Run 483 stress 0.05959901 
    ## Run 484 stress 0.0493595 
    ## Run 485 stress 0.05959895 
    ## Run 486 stress 0.04935947 
    ## Run 487 stress 0.05216161 
    ## Run 488 stress 0.06185061 
    ## Run 489 stress 0.05728074 
    ## Run 490 stress 0.06400736 
    ## Run 491 stress 0.05572689 
    ## Run 492 stress 0.05216146 
    ## Run 493 stress 0.05572688 
    ## Run 494 stress 0.05643263 
    ## Run 495 stress 0.06400723 
    ## Run 496 stress 0.04612758 
    ## ... Procrustes: rmse 3.963059e-05  max resid 7.974688e-05 
    ## ... Similar to previous best
    ## Run 497 stress 0.04860424 
    ## Run 498 stress 0.05451493 
    ## Run 499 stress 0.05695873 
    ## Run 500 stress 0.05607344 
    ## *** Best solution repeated 17 times

``` r
# Stratified lakes
PD_beta_S_NMDS <- metaMDS(PD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.02693829 
    ## Run 1 stress 0.02693829 
    ## ... New best solution
    ## ... Procrustes: rmse 5.842146e-05  max resid 0.0001087134 
    ## ... Similar to previous best
    ## Run 2 stress 0.02693831 
    ## ... Procrustes: rmse 0.0001230457  max resid 0.0002204617 
    ## ... Similar to previous best
    ## Run 3 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001779398  max resid 0.0003101397 
    ## ... Similar to previous best
    ## Run 4 stress 0.0269383 
    ## ... Procrustes: rmse 3.714232e-05  max resid 5.047352e-05 
    ## ... Similar to previous best
    ## Run 5 stress 0.1889141 
    ## Run 6 stress 0.07547299 
    ## Run 7 stress 0.1862384 
    ## Run 8 stress 0.1862387 
    ## Run 9 stress 0.1889141 
    ## Run 10 stress 0.02693831 
    ## ... Procrustes: rmse 5.465495e-05  max resid 8.137694e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.07547299 
    ## Run 12 stress 0.179905 
    ## Run 13 stress 0.075473 
    ## Run 14 stress 0.07547299 
    ## Run 15 stress 0.08382871 
    ## Run 16 stress 0.07374801 
    ## Run 17 stress 0.07374813 
    ## Run 18 stress 0.07374806 
    ## Run 19 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001415234  max resid 0.0002499792 
    ## ... Similar to previous best
    ## Run 20 stress 0.07547302 
    ## Run 21 stress 0.07374802 
    ## Run 22 stress 0.07374808 
    ## Run 23 stress 0.02693839 
    ## ... Procrustes: rmse 0.00020114  max resid 0.0003540468 
    ## ... Similar to previous best
    ## Run 24 stress 0.07374804 
    ## Run 25 stress 0.07374804 
    ## Run 26 stress 0.02693839 
    ## ... Procrustes: rmse 0.0002072434  max resid 0.0003641142 
    ## ... Similar to previous best
    ## Run 27 stress 0.179905 
    ## Run 28 stress 0.08382857 
    ## Run 29 stress 0.075473 
    ## Run 30 stress 0.02693829 
    ## ... Procrustes: rmse 6.613887e-05  max resid 0.0001212806 
    ## ... Similar to previous best
    ## Run 31 stress 0.02693829 
    ## ... Procrustes: rmse 1.249393e-05  max resid 2.310638e-05 
    ## ... Similar to previous best
    ## Run 32 stress 0.02693832 
    ## ... Procrustes: rmse 0.000125387  max resid 0.0002244113 
    ## ... Similar to previous best
    ## Run 33 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001589613  max resid 0.00028373 
    ## ... Similar to previous best
    ## Run 34 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001925093  max resid 0.0003396082 
    ## ... Similar to previous best
    ## Run 35 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001019081  max resid 0.0001629837 
    ## ... Similar to previous best
    ## Run 36 stress 0.02693831 
    ## ... Procrustes: rmse 6.222291e-05  max resid 9.463134e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.09391893 
    ## Run 38 stress 0.07374805 
    ## Run 39 stress 0.02693829 
    ## ... New best solution
    ## ... Procrustes: rmse 4.717169e-05  max resid 8.936076e-05 
    ## ... Similar to previous best
    ## Run 40 stress 0.1862385 
    ## Run 41 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001082192  max resid 0.0001850764 
    ## ... Similar to previous best
    ## Run 42 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001680607  max resid 0.0002883438 
    ## ... Similar to previous best
    ## Run 43 stress 0.07374802 
    ## Run 44 stress 0.02693829 
    ## ... Procrustes: rmse 2.037967e-05  max resid 3.350062e-05 
    ## ... Similar to previous best
    ## Run 45 stress 0.1799052 
    ## Run 46 stress 0.07374803 
    ## Run 47 stress 0.075473 
    ## Run 48 stress 0.02693831 
    ## ... Procrustes: rmse 6.817093e-05  max resid 0.0001082489 
    ## ... Similar to previous best
    ## Run 49 stress 0.02693829 
    ## ... Procrustes: rmse 5.390076e-05  max resid 9.555174e-05 
    ## ... Similar to previous best
    ## Run 50 stress 0.1799051 
    ## Run 51 stress 0.02693829 
    ## ... Procrustes: rmse 6.200575e-05  max resid 0.0001068889 
    ## ... Similar to previous best
    ## Run 52 stress 0.02693829 
    ## ... Procrustes: rmse 3.276922e-05  max resid 5.478006e-05 
    ## ... Similar to previous best
    ## Run 53 stress 0.075473 
    ## Run 54 stress 0.075473 
    ## Run 55 stress 0.07374809 
    ## Run 56 stress 0.08382874 
    ## Run 57 stress 0.179905 
    ## Run 58 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001656359  max resid 0.0002826999 
    ## ... Similar to previous best
    ## Run 59 stress 0.179905 
    ## Run 60 stress 0.08382837 
    ## Run 61 stress 0.08382874 
    ## Run 62 stress 0.179905 
    ## Run 63 stress 0.08382841 
    ## Run 64 stress 0.1889141 
    ## Run 65 stress 0.07374802 
    ## Run 66 stress 0.07374803 
    ## Run 67 stress 0.07374812 
    ## Run 68 stress 0.02693833 
    ## ... Procrustes: rmse 0.000100954  max resid 0.0001724091 
    ## ... Similar to previous best
    ## Run 69 stress 0.08382852 
    ## Run 70 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001423603  max resid 0.0002447867 
    ## ... Similar to previous best
    ## Run 71 stress 0.1889141 
    ## Run 72 stress 0.08382858 
    ## Run 73 stress 0.07547299 
    ## Run 74 stress 0.0269383 
    ## ... Procrustes: rmse 3.093368e-05  max resid 5.201973e-05 
    ## ... Similar to previous best
    ## Run 75 stress 0.07547299 
    ## Run 76 stress 0.02693833 
    ## ... Procrustes: rmse 8.415892e-05  max resid 0.0001440275 
    ## ... Similar to previous best
    ## Run 77 stress 0.07374802 
    ## Run 78 stress 0.07374802 
    ## Run 79 stress 0.02693829 
    ## ... Procrustes: rmse 2.703544e-05  max resid 4.515526e-05 
    ## ... Similar to previous best
    ## Run 80 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001061444  max resid 0.000181066 
    ## ... Similar to previous best
    ## Run 81 stress 0.08382852 
    ## Run 82 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001263416  max resid 0.0002185902 
    ## ... Similar to previous best
    ## Run 83 stress 0.07374803 
    ## Run 84 stress 0.179905 
    ## Run 85 stress 0.08382873 
    ## Run 86 stress 0.02693829 
    ## ... New best solution
    ## ... Procrustes: rmse 2.773688e-06  max resid 4.289099e-06 
    ## ... Similar to previous best
    ## Run 87 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001824791  max resid 0.0003127793 
    ## ... Similar to previous best
    ## Run 88 stress 0.07547299 
    ## Run 89 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001285611  max resid 0.0002205662 
    ## ... Similar to previous best
    ## Run 90 stress 0.1889141 
    ## Run 91 stress 0.09391893 
    ## Run 92 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001448728  max resid 0.0002483174 
    ## ... Similar to previous best
    ## Run 93 stress 0.07547301 
    ## Run 94 stress 0.02693837 
    ## ... Procrustes: rmse 0.000126454  max resid 0.0002177269 
    ## ... Similar to previous best
    ## Run 95 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001889133  max resid 0.0003224497 
    ## ... Similar to previous best
    ## Run 96 stress 0.02693829 
    ## ... Procrustes: rmse 1.153082e-05  max resid 1.882638e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.02693829 
    ## ... Procrustes: rmse 3.614236e-05  max resid 6.188344e-05 
    ## ... Similar to previous best
    ## Run 98 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001787731  max resid 0.0003064705 
    ## ... Similar to previous best
    ## Run 99 stress 0.02693837 
    ## ... Procrustes: rmse 0.000128106  max resid 0.0002207116 
    ## ... Similar to previous best
    ## Run 100 stress 0.1756177 
    ## Run 101 stress 0.075473 
    ## Run 102 stress 0.07547299 
    ## Run 103 stress 0.09391912 
    ## Run 104 stress 0.2273964 
    ## Run 105 stress 0.1889143 
    ## Run 106 stress 0.08382839 
    ## Run 107 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001121277  max resid 0.0001931554 
    ## ... Similar to previous best
    ## Run 108 stress 0.0269384 
    ## ... Procrustes: rmse 0.000172178  max resid 0.0002935666 
    ## ... Similar to previous best
    ## Run 109 stress 0.08382873 
    ## Run 110 stress 0.02693833 
    ## ... Procrustes: rmse 9.192292e-05  max resid 0.000158198 
    ## ... Similar to previous best
    ## Run 111 stress 0.073748 
    ## Run 112 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001354828  max resid 0.0002318739 
    ## ... Similar to previous best
    ## Run 113 stress 0.0269383 
    ## ... Procrustes: rmse 4.866458e-05  max resid 8.251822e-05 
    ## ... Similar to previous best
    ## Run 114 stress 0.075473 
    ## Run 115 stress 0.07374814 
    ## Run 116 stress 0.07374805 
    ## Run 117 stress 0.02693828 
    ## ... New best solution
    ## ... Procrustes: rmse 2.749591e-05  max resid 4.735663e-05 
    ## ... Similar to previous best
    ## Run 118 stress 0.07374808 
    ## Run 119 stress 0.1889143 
    ## Run 120 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001332712  max resid 0.0002285139 
    ## ... Similar to previous best
    ## Run 121 stress 0.283614 
    ## Run 122 stress 0.02693831 
    ## ... Procrustes: rmse 8.910572e-05  max resid 0.0001532623 
    ## ... Similar to previous best
    ## Run 123 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001634087  max resid 0.0002802605 
    ## ... Similar to previous best
    ## Run 124 stress 0.08382844 
    ## Run 125 stress 0.07374801 
    ## Run 126 stress 0.07374814 
    ## Run 127 stress 0.02693833 
    ## ... Procrustes: rmse 9.013987e-05  max resid 0.0001570573 
    ## ... Similar to previous best
    ## Run 128 stress 0.07374802 
    ## Run 129 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001857141  max resid 0.0003166476 
    ## ... Similar to previous best
    ## Run 130 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001188718  max resid 0.0002043303 
    ## ... Similar to previous best
    ## Run 131 stress 0.07374801 
    ## Run 132 stress 0.1756177 
    ## Run 133 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001330995  max resid 0.0002264922 
    ## ... Similar to previous best
    ## Run 134 stress 0.073748 
    ## Run 135 stress 0.1862389 
    ## Run 136 stress 0.07374806 
    ## Run 137 stress 0.0737481 
    ## Run 138 stress 0.1799051 
    ## Run 139 stress 0.179905 
    ## Run 140 stress 0.179905 
    ## Run 141 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001103021  max resid 0.000188442 
    ## ... Similar to previous best
    ## Run 142 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001209208  max resid 0.0002067572 
    ## ... Similar to previous best
    ## Run 143 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001394912  max resid 0.0002375594 
    ## ... Similar to previous best
    ## Run 144 stress 0.179905 
    ## Run 145 stress 0.07547301 
    ## Run 146 stress 0.07374801 
    ## Run 147 stress 0.2273964 
    ## Run 148 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001877158  max resid 0.0003207298 
    ## ... Similar to previous best
    ## Run 149 stress 0.02693829 
    ## ... Procrustes: rmse 2.739327e-05  max resid 4.716214e-05 
    ## ... Similar to previous best
    ## Run 150 stress 0.02693832 
    ## ... Procrustes: rmse 9.153897e-05  max resid 0.0001563218 
    ## ... Similar to previous best
    ## Run 151 stress 0.07547299 
    ## Run 152 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001495655  max resid 0.0002567723 
    ## ... Similar to previous best
    ## Run 153 stress 0.08382823 
    ## Run 154 stress 0.07374805 
    ## Run 155 stress 0.02693829 
    ## ... Procrustes: rmse 2.226468e-05  max resid 4.333513e-05 
    ## ... Similar to previous best
    ## Run 156 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001662322  max resid 0.000285184 
    ## ... Similar to previous best
    ## Run 157 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001137164  max resid 0.0001954351 
    ## ... Similar to previous best
    ## Run 158 stress 0.179905 
    ## Run 159 stress 0.07374806 
    ## Run 160 stress 0.07547299 
    ## Run 161 stress 0.179905 
    ## Run 162 stress 0.07547299 
    ## Run 163 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001698  max resid 0.0002904341 
    ## ... Similar to previous best
    ## Run 164 stress 0.08382851 
    ## Run 165 stress 0.02693829 
    ## ... Procrustes: rmse 3.804725e-05  max resid 6.537679e-05 
    ## ... Similar to previous best
    ## Run 166 stress 0.075473 
    ## Run 167 stress 0.07374811 
    ## Run 168 stress 0.02693832 
    ## ... Procrustes: rmse 8.685413e-05  max resid 0.0001485524 
    ## ... Similar to previous best
    ## Run 169 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001307714  max resid 0.0002246542 
    ## ... Similar to previous best
    ## Run 170 stress 0.02693831 
    ## ... Procrustes: rmse 7.354922e-05  max resid 0.000125058 
    ## ... Similar to previous best
    ## Run 171 stress 0.1889142 
    ## Run 172 stress 0.07374802 
    ## Run 173 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001341412  max resid 0.0002286109 
    ## ... Similar to previous best
    ## Run 174 stress 0.02693829 
    ## ... Procrustes: rmse 4.365813e-05  max resid 7.413747e-05 
    ## ... Similar to previous best
    ## Run 175 stress 0.075473 
    ## Run 176 stress 0.07374807 
    ## Run 177 stress 0.073748 
    ## Run 178 stress 0.1862385 
    ## Run 179 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001538143  max resid 0.0002635742 
    ## ... Similar to previous best
    ## Run 180 stress 0.08382821 
    ## Run 181 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001561769  max resid 0.0002687725 
    ## ... Similar to previous best
    ## Run 182 stress 0.179905 
    ## Run 183 stress 0.0269383 
    ## ... Procrustes: rmse 8.576012e-05  max resid 0.0001440519 
    ## ... Similar to previous best
    ## Run 184 stress 0.02693832 
    ## ... Procrustes: rmse 9.293129e-05  max resid 0.0001590102 
    ## ... Similar to previous best
    ## Run 185 stress 0.1756177 
    ## Run 186 stress 0.07374802 
    ## Run 187 stress 0.07374803 
    ## Run 188 stress 0.02693829 
    ## ... Procrustes: rmse 2.653977e-05  max resid 4.511025e-05 
    ## ... Similar to previous best
    ## Run 189 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001537841  max resid 0.0002641448 
    ## ... Similar to previous best
    ## Run 190 stress 0.02693832 
    ## ... Procrustes: rmse 8.827802e-05  max resid 0.0001520782 
    ## ... Similar to previous best
    ## Run 191 stress 0.09391923 
    ## Run 192 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001058386  max resid 0.0001817592 
    ## ... Similar to previous best
    ## Run 193 stress 0.07374809 
    ## Run 194 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001176556  max resid 0.0002078629 
    ## ... Similar to previous best
    ## Run 195 stress 0.07374801 
    ## Run 196 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001510473  max resid 0.0002582572 
    ## ... Similar to previous best
    ## Run 197 stress 0.073748 
    ## Run 198 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001292716  max resid 0.0002213883 
    ## ... Similar to previous best
    ## Run 199 stress 0.02693835 
    ## ... Procrustes: rmse 0.000128569  max resid 0.000220509 
    ## ... Similar to previous best
    ## Run 200 stress 0.073748 
    ## Run 201 stress 0.07374801 
    ## Run 202 stress 0.1756177 
    ## Run 203 stress 0.07374803 
    ## Run 204 stress 0.02693829 
    ## ... Procrustes: rmse 4.438299e-05  max resid 7.577934e-05 
    ## ... Similar to previous best
    ## Run 205 stress 0.179905 
    ## Run 206 stress 0.02693831 
    ## ... Procrustes: rmse 7.50908e-05  max resid 0.0001285832 
    ## ... Similar to previous best
    ## Run 207 stress 0.08382876 
    ## Run 208 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001410373  max resid 0.0002411822 
    ## ... Similar to previous best
    ## Run 209 stress 0.179905 
    ## Run 210 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001724366  max resid 0.0002952131 
    ## ... Similar to previous best
    ## Run 211 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001689847  max resid 0.0002899737 
    ## ... Similar to previous best
    ## Run 212 stress 0.07547299 
    ## Run 213 stress 0.07547302 
    ## Run 214 stress 0.08382847 
    ## Run 215 stress 0.073748 
    ## Run 216 stress 0.179905 
    ## Run 217 stress 0.0269383 
    ## ... Procrustes: rmse 6.83298e-05  max resid 0.0001171399 
    ## ... Similar to previous best
    ## Run 218 stress 0.1889141 
    ## Run 219 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001511445  max resid 0.0002595428 
    ## ... Similar to previous best
    ## Run 220 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001204198  max resid 0.0002038754 
    ## ... Similar to previous best
    ## Run 221 stress 0.07374802 
    ## Run 222 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001055632  max resid 0.0001806879 
    ## ... Similar to previous best
    ## Run 223 stress 0.07374801 
    ## Run 224 stress 0.179905 
    ## Run 225 stress 0.07374803 
    ## Run 226 stress 0.02693833 
    ## ... Procrustes: rmse 9.203431e-05  max resid 0.0001585156 
    ## ... Similar to previous best
    ## Run 227 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001814801  max resid 0.0003119332 
    ## ... Similar to previous best
    ## Run 228 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001530633  max resid 0.0002600586 
    ## ... Similar to previous best
    ## Run 229 stress 0.08382822 
    ## Run 230 stress 0.02693831 
    ## ... Procrustes: rmse 9.482876e-05  max resid 0.0001625032 
    ## ... Similar to previous best
    ## Run 231 stress 0.1756177 
    ## Run 232 stress 0.08382851 
    ## Run 233 stress 0.1799051 
    ## Run 234 stress 0.1889141 
    ## Run 235 stress 0.07547299 
    ## Run 236 stress 0.1756177 
    ## Run 237 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001448862  max resid 0.0002480361 
    ## ... Similar to previous best
    ## Run 238 stress 0.02693831 
    ## ... Procrustes: rmse 9.610198e-05  max resid 0.0001645426 
    ## ... Similar to previous best
    ## Run 239 stress 0.07374801 
    ## Run 240 stress 0.179905 
    ## Run 241 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001231511  max resid 0.000210489 
    ## ... Similar to previous best
    ## Run 242 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001010969  max resid 0.000172403 
    ## ... Similar to previous best
    ## Run 243 stress 0.0737481 
    ## Run 244 stress 0.07374808 
    ## Run 245 stress 0.07374804 
    ## Run 246 stress 0.07374806 
    ## Run 247 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001726899  max resid 0.0002958512 
    ## ... Similar to previous best
    ## Run 248 stress 0.07547301 
    ## Run 249 stress 0.02693832 
    ## ... Procrustes: rmse 8.770723e-05  max resid 0.0001507191 
    ## ... Similar to previous best
    ## Run 250 stress 0.1756177 
    ## Run 251 stress 0.08382832 
    ## Run 252 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001210268  max resid 0.000207036 
    ## ... Similar to previous best
    ## Run 253 stress 0.075473 
    ## Run 254 stress 0.2273964 
    ## Run 255 stress 0.1756177 
    ## Run 256 stress 0.075473 
    ## Run 257 stress 0.07374811 
    ## Run 258 stress 0.02693832 
    ## ... Procrustes: rmse 8.998882e-05  max resid 0.0001606943 
    ## ... Similar to previous best
    ## Run 259 stress 0.1889141 
    ## Run 260 stress 0.2714151 
    ## Run 261 stress 0.179905 
    ## Run 262 stress 0.07374801 
    ## Run 263 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001343054  max resid 0.0002302305 
    ## ... Similar to previous best
    ## Run 264 stress 0.3083098 
    ## Run 265 stress 0.07374801 
    ## Run 266 stress 0.073748 
    ## Run 267 stress 0.075473 
    ## Run 268 stress 0.07547299 
    ## Run 269 stress 0.02693834 
    ## ... Procrustes: rmse 0.000102485  max resid 0.0001756404 
    ## ... Similar to previous best
    ## Run 270 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001715232  max resid 0.0002921094 
    ## ... Similar to previous best
    ## Run 271 stress 0.07374807 
    ## Run 272 stress 0.1756177 
    ## Run 273 stress 0.07547301 
    ## Run 274 stress 0.179905 
    ## Run 275 stress 0.07547299 
    ## Run 276 stress 0.07374802 
    ## Run 277 stress 0.07374803 
    ## Run 278 stress 0.02693831 
    ## ... Procrustes: rmse 7.091119e-05  max resid 0.0001211483 
    ## ... Similar to previous best
    ## Run 279 stress 0.075473 
    ## Run 280 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001296696  max resid 0.0002235154 
    ## ... Similar to previous best
    ## Run 281 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001375394  max resid 0.0002350644 
    ## ... Similar to previous best
    ## Run 282 stress 0.1799051 
    ## Run 283 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001203548  max resid 0.0002048541 
    ## ... Similar to previous best
    ## Run 284 stress 0.07374801 
    ## Run 285 stress 0.08382873 
    ## Run 286 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001606309  max resid 0.0002748032 
    ## ... Similar to previous best
    ## Run 287 stress 0.02693829 
    ## ... Procrustes: rmse 2.987598e-05  max resid 5.13505e-05 
    ## ... Similar to previous best
    ## Run 288 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001219719  max resid 0.0002081144 
    ## ... Similar to previous best
    ## Run 289 stress 0.07547302 
    ## Run 290 stress 0.08382848 
    ## Run 291 stress 0.09391904 
    ## Run 292 stress 0.02693841 
    ## ... Procrustes: rmse 0.0002028375  max resid 0.0003414052 
    ## ... Similar to previous best
    ## Run 293 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001131099  max resid 0.0001940547 
    ## ... Similar to previous best
    ## Run 294 stress 0.0838283 
    ## Run 295 stress 0.07374804 
    ## Run 296 stress 0.07547299 
    ## Run 297 stress 0.179905 
    ## Run 298 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001046307  max resid 0.0001798286 
    ## ... Similar to previous best
    ## Run 299 stress 0.1900877 
    ## Run 300 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001642476  max resid 0.0002803813 
    ## ... Similar to previous best
    ## Run 301 stress 0.179905 
    ## Run 302 stress 0.1756177 
    ## Run 303 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001683532  max resid 0.0002879296 
    ## ... Similar to previous best
    ## Run 304 stress 0.02693831 
    ## ... Procrustes: rmse 9.939683e-05  max resid 0.000171114 
    ## ... Similar to previous best
    ## Run 305 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001110615  max resid 0.0001892159 
    ## ... Similar to previous best
    ## Run 306 stress 0.07547299 
    ## Run 307 stress 0.073748 
    ## Run 308 stress 0.07374806 
    ## Run 309 stress 0.07374801 
    ## Run 310 stress 0.07374802 
    ## Run 311 stress 0.075473 
    ## Run 312 stress 0.02693829 
    ## ... Procrustes: rmse 6.188778e-05  max resid 0.0001059363 
    ## ... Similar to previous best
    ## Run 313 stress 0.1756177 
    ## Run 314 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001447239  max resid 0.0002487245 
    ## ... Similar to previous best
    ## Run 315 stress 0.0269383 
    ## ... Procrustes: rmse 7.650443e-05  max resid 0.0001317755 
    ## ... Similar to previous best
    ## Run 316 stress 0.08382823 
    ## Run 317 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001445919  max resid 0.000246998 
    ## ... Similar to previous best
    ## Run 318 stress 0.02693838 
    ## ... Procrustes: rmse 0.0001548731  max resid 0.0002846818 
    ## ... Similar to previous best
    ## Run 319 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001614339  max resid 0.0002771823 
    ## ... Similar to previous best
    ## Run 320 stress 0.02693833 
    ## ... Procrustes: rmse 0.000128833  max resid 0.0002207238 
    ## ... Similar to previous best
    ## Run 321 stress 0.07547301 
    ## Run 322 stress 0.3083098 
    ## Run 323 stress 0.07374803 
    ## Run 324 stress 0.075473 
    ## Run 325 stress 0.02693828 
    ## ... New best solution
    ## ... Procrustes: rmse 7.869097e-07  max resid 1.479961e-06 
    ## ... Similar to previous best
    ## Run 326 stress 0.0269384 
    ## ... Procrustes: rmse 0.000176905  max resid 0.0003032912 
    ## ... Similar to previous best
    ## Run 327 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001392444  max resid 0.0002379539 
    ## ... Similar to previous best
    ## Run 328 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001639459  max resid 0.0002807565 
    ## ... Similar to previous best
    ## Run 329 stress 0.07374807 
    ## Run 330 stress 0.1756177 
    ## Run 331 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001573921  max resid 0.0002703614 
    ## ... Similar to previous best
    ## Run 332 stress 0.07374804 
    ## Run 333 stress 0.07374811 
    ## Run 334 stress 0.07547299 
    ## Run 335 stress 0.0269384 
    ## ... Procrustes: rmse 0.0001746961  max resid 0.0002985601 
    ## ... Similar to previous best
    ## Run 336 stress 0.07547299 
    ## Run 337 stress 0.073748 
    ## Run 338 stress 0.1862386 
    ## Run 339 stress 0.08382855 
    ## Run 340 stress 0.1756177 
    ## Run 341 stress 0.08382858 
    ## Run 342 stress 0.07374802 
    ## Run 343 stress 0.073748 
    ## Run 344 stress 0.07374808 
    ## Run 345 stress 0.02693842 
    ## ... Procrustes: rmse 0.0001805468  max resid 0.0003103875 
    ## ... Similar to previous best
    ## Run 346 stress 0.07547299 
    ## Run 347 stress 0.1756177 
    ## Run 348 stress 0.07547299 
    ## Run 349 stress 0.07374802 
    ## Run 350 stress 0.02693828 
    ## ... Procrustes: rmse 2.35507e-05  max resid 4.140269e-05 
    ## ... Similar to previous best
    ## Run 351 stress 0.02693831 
    ## ... Procrustes: rmse 8.914181e-05  max resid 0.0001515533 
    ## ... Similar to previous best
    ## Run 352 stress 0.02693839 
    ## ... Procrustes: rmse 0.0001768816  max resid 0.0003029448 
    ## ... Similar to previous best
    ## Run 353 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001595788  max resid 0.0002728463 
    ## ... Similar to previous best
    ## Run 354 stress 0.07374801 
    ## Run 355 stress 0.07374806 
    ## Run 356 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001320118  max resid 0.0002255237 
    ## ... Similar to previous best
    ## Run 357 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001723301  max resid 0.0002948132 
    ## ... Similar to previous best
    ## Run 358 stress 0.02693829 
    ## ... Procrustes: rmse 4.301982e-05  max resid 7.383845e-05 
    ## ... Similar to previous best
    ## Run 359 stress 0.073748 
    ## Run 360 stress 0.02693833 
    ## ... Procrustes: rmse 0.000125797  max resid 0.0002155156 
    ## ... Similar to previous best
    ## Run 361 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001465592  max resid 0.0002507315 
    ## ... Similar to previous best
    ## Run 362 stress 0.07547301 
    ## Run 363 stress 0.0269383 
    ## ... Procrustes: rmse 8.012945e-05  max resid 0.0001373991 
    ## ... Similar to previous best
    ## Run 364 stress 0.07374805 
    ## Run 365 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001505229  max resid 0.0002569887 
    ## ... Similar to previous best
    ## Run 366 stress 0.1889143 
    ## Run 367 stress 0.08382874 
    ## Run 368 stress 0.0269383 
    ## ... Procrustes: rmse 7.63244e-05  max resid 0.0001304091 
    ## ... Similar to previous best
    ## Run 369 stress 0.0269383 
    ## ... Procrustes: rmse 7.920693e-05  max resid 0.0001355853 
    ## ... Similar to previous best
    ## Run 370 stress 0.2842805 
    ## Run 371 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001051222  max resid 0.000179577 
    ## ... Similar to previous best
    ## Run 372 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001516428  max resid 0.000259387 
    ## ... Similar to previous best
    ## Run 373 stress 0.07374806 
    ## Run 374 stress 0.0269384 
    ## ... Procrustes: rmse 0.000186981  max resid 0.0003191126 
    ## ... Similar to previous best
    ## Run 375 stress 0.08382821 
    ## Run 376 stress 0.07547299 
    ## Run 377 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001452803  max resid 0.0002472186 
    ## ... Similar to previous best
    ## Run 378 stress 0.07374812 
    ## Run 379 stress 0.08382846 
    ## Run 380 stress 0.02693829 
    ## ... Procrustes: rmse 3.901959e-05  max resid 6.623053e-05 
    ## ... Similar to previous best
    ## Run 381 stress 0.1889143 
    ## Run 382 stress 0.07547299 
    ## Run 383 stress 0.08382833 
    ## Run 384 stress 0.07374802 
    ## Run 385 stress 0.07374802 
    ## Run 386 stress 0.02693832 
    ## ... Procrustes: rmse 8.89343e-05  max resid 0.0001529134 
    ## ... Similar to previous best
    ## Run 387 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001329973  max resid 0.0002278906 
    ## ... Similar to previous best
    ## Run 388 stress 0.07374802 
    ## Run 389 stress 0.02693831 
    ## ... Procrustes: rmse 8.706398e-05  max resid 0.0001482631 
    ## ... Similar to previous best
    ## Run 390 stress 0.02693829 
    ## ... Procrustes: rmse 1.123104e-05  max resid 1.582639e-05 
    ## ... Similar to previous best
    ## Run 391 stress 0.08382851 
    ## Run 392 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001556573  max resid 0.0002671562 
    ## ... Similar to previous best
    ## Run 393 stress 0.02693831 
    ## ... Procrustes: rmse 6.652441e-05  max resid 0.0001148365 
    ## ... Similar to previous best
    ## Run 394 stress 0.07374807 
    ## Run 395 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001619943  max resid 0.0002774441 
    ## ... Similar to previous best
    ## Run 396 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001104063  max resid 0.0001896958 
    ## ... Similar to previous best
    ## Run 397 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001011388  max resid 0.0001741434 
    ## ... Similar to previous best
    ## Run 398 stress 0.1862384 
    ## Run 399 stress 0.02693831 
    ## ... Procrustes: rmse 0.0001038869  max resid 0.0001774456 
    ## ... Similar to previous best
    ## Run 400 stress 0.07374801 
    ## Run 401 stress 0.02693837 
    ## ... Procrustes: rmse 0.000159389  max resid 0.0002740024 
    ## ... Similar to previous best
    ## Run 402 stress 0.02693832 
    ## ... Procrustes: rmse 9.296069e-05  max resid 0.0001588027 
    ## ... Similar to previous best
    ## Run 403 stress 0.1889141 
    ## Run 404 stress 0.02693832 
    ## ... Procrustes: rmse 0.000111089  max resid 0.0001894755 
    ## ... Similar to previous best
    ## Run 405 stress 0.073748 
    ## Run 406 stress 0.02693829 
    ## ... Procrustes: rmse 5.388906e-05  max resid 9.177159e-05 
    ## ... Similar to previous best
    ## Run 407 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001031267  max resid 0.0001727386 
    ## ... Similar to previous best
    ## Run 408 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001623151  max resid 0.0002782264 
    ## ... Similar to previous best
    ## Run 409 stress 0.07374801 
    ## Run 410 stress 0.0838284 
    ## Run 411 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001222632  max resid 0.0002092057 
    ## ... Similar to previous best
    ## Run 412 stress 0.08382859 
    ## Run 413 stress 0.0269383 
    ## ... Procrustes: rmse 7.333565e-05  max resid 0.0001252857 
    ## ... Similar to previous best
    ## Run 414 stress 0.02693835 
    ## ... Procrustes: rmse 0.000146479  max resid 0.0002511552 
    ## ... Similar to previous best
    ## Run 415 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001156438  max resid 0.0001988293 
    ## ... Similar to previous best
    ## Run 416 stress 0.07374801 
    ## Run 417 stress 0.02693829 
    ## ... Procrustes: rmse 2.152368e-05  max resid 3.695302e-05 
    ## ... Similar to previous best
    ## Run 418 stress 0.07374802 
    ## Run 419 stress 0.07374803 
    ## Run 420 stress 0.08382837 
    ## Run 421 stress 0.07547299 
    ## Run 422 stress 0.07374813 
    ## Run 423 stress 0.07374804 
    ## Run 424 stress 0.07374808 
    ## Run 425 stress 0.07374802 
    ## Run 426 stress 0.08382863 
    ## Run 427 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001422637  max resid 0.0002444315 
    ## ... Similar to previous best
    ## Run 428 stress 0.07374804 
    ## Run 429 stress 0.1756177 
    ## Run 430 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001603896  max resid 0.0002752543 
    ## ... Similar to previous best
    ## Run 431 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001267741  max resid 0.0002155773 
    ## ... Similar to previous best
    ## Run 432 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001017245  max resid 0.0001735495 
    ## ... Similar to previous best
    ## Run 433 stress 0.02693839 
    ## ... Procrustes: rmse 0.000150379  max resid 0.0002590016 
    ## ... Similar to previous best
    ## Run 434 stress 0.075473 
    ## Run 435 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001308678  max resid 0.0002253613 
    ## ... Similar to previous best
    ## Run 436 stress 0.02693829 
    ## ... Procrustes: rmse 3.42087e-05  max resid 5.887561e-05 
    ## ... Similar to previous best
    ## Run 437 stress 0.1889141 
    ## Run 438 stress 0.02693837 
    ## ... Procrustes: rmse 0.0001414498  max resid 0.0002405435 
    ## ... Similar to previous best
    ## Run 439 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001233712  max resid 0.0002111185 
    ## ... Similar to previous best
    ## Run 440 stress 0.0269383 
    ## ... Procrustes: rmse 6.121886e-05  max resid 0.0001042673 
    ## ... Similar to previous best
    ## Run 441 stress 0.02693831 
    ## ... Procrustes: rmse 9.51018e-05  max resid 0.0001633444 
    ## ... Similar to previous best
    ## Run 442 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001424259  max resid 0.0002424436 
    ## ... Similar to previous best
    ## Run 443 stress 0.07374809 
    ## Run 444 stress 0.07374813 
    ## Run 445 stress 0.02693832 
    ## ... Procrustes: rmse 0.000115819  max resid 0.0001979343 
    ## ... Similar to previous best
    ## Run 446 stress 0.1889144 
    ## Run 447 stress 0.08382829 
    ## Run 448 stress 0.0269383 
    ## ... Procrustes: rmse 6.041076e-05  max resid 0.0001037804 
    ## ... Similar to previous best
    ## Run 449 stress 0.07374802 
    ## Run 450 stress 0.07374803 
    ## Run 451 stress 0.07547299 
    ## Run 452 stress 0.08382853 
    ## Run 453 stress 0.07547301 
    ## Run 454 stress 0.02693832 
    ## ... Procrustes: rmse 0.0001108815  max resid 0.0001901441 
    ## ... Similar to previous best
    ## Run 455 stress 0.0269383 
    ## ... Procrustes: rmse 6.329556e-05  max resid 0.0001084476 
    ## ... Similar to previous best
    ## Run 456 stress 0.07374801 
    ## Run 457 stress 0.07374809 
    ## Run 458 stress 0.07374813 
    ## Run 459 stress 0.07547301 
    ## Run 460 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001128266  max resid 0.0001932858 
    ## ... Similar to previous best
    ## Run 461 stress 0.08382856 
    ## Run 462 stress 0.07547299 
    ## Run 463 stress 0.07374801 
    ## Run 464 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001202569  max resid 0.0002061276 
    ## ... Similar to previous best
    ## Run 465 stress 0.1756177 
    ## Run 466 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001585457  max resid 0.0002712197 
    ## ... Similar to previous best
    ## Run 467 stress 0.07374807 
    ## Run 468 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001280417  max resid 0.0002186658 
    ## ... Similar to previous best
    ## Run 469 stress 0.1889141 
    ## Run 470 stress 0.02693836 
    ## ... Procrustes: rmse 0.000156096  max resid 0.0002495772 
    ## ... Similar to previous best
    ## Run 471 stress 0.07374802 
    ## Run 472 stress 0.0269383 
    ## ... Procrustes: rmse 4.323688e-05  max resid 7.320157e-05 
    ## ... Similar to previous best
    ## Run 473 stress 0.07374802 
    ## Run 474 stress 0.07374803 
    ## Run 475 stress 0.0269383 
    ## ... Procrustes: rmse 6.284562e-05  max resid 0.0001069285 
    ## ... Similar to previous best
    ## Run 476 stress 0.08382834 
    ## Run 477 stress 0.07547304 
    ## Run 478 stress 0.07374803 
    ## Run 479 stress 0.07374803 
    ## Run 480 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001433685  max resid 0.0002451879 
    ## ... Similar to previous best
    ## Run 481 stress 0.09391915 
    ## Run 482 stress 0.08382875 
    ## Run 483 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001033157  max resid 0.0001767638 
    ## ... Similar to previous best
    ## Run 484 stress 0.1889141 
    ## Run 485 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001118639  max resid 0.0001920572 
    ## ... Similar to previous best
    ## Run 486 stress 0.02693835 
    ## ... Procrustes: rmse 0.0001243768  max resid 0.0002127713 
    ## ... Similar to previous best
    ## Run 487 stress 0.02693833 
    ## ... Procrustes: rmse 0.0001030083  max resid 0.0001754845 
    ## ... Similar to previous best
    ## Run 488 stress 0.02693828 
    ## ... New best solution
    ## ... Procrustes: rmse 6.817654e-06  max resid 1.14754e-05 
    ## ... Similar to previous best
    ## Run 489 stress 0.179905 
    ## Run 490 stress 0.02693831 
    ## ... Procrustes: rmse 8.313747e-05  max resid 0.0001426792 
    ## ... Similar to previous best
    ## Run 491 stress 0.07547299 
    ## Run 492 stress 0.02693836 
    ## ... Procrustes: rmse 0.0001467377  max resid 0.000251006 
    ## ... Similar to previous best
    ## Run 493 stress 0.073748 
    ## Run 494 stress 0.02693834 
    ## ... Procrustes: rmse 0.0001214938  max resid 0.0002084831 
    ## ... Similar to previous best
    ## Run 495 stress 0.07374805 
    ## Run 496 stress 0.1862385 
    ## Run 497 stress 0.07374804 
    ## Run 498 stress 0.07374805 
    ## Run 499 stress 0.02693829 
    ## ... Procrustes: rmse 3.28129e-05  max resid 5.602584e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.07374802 
    ## *** Best solution repeated 5 times

``` r
### Environmental
# Surveyed sites 
PD_beta_env_NMDS <- metaMDS(PD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06330868 
    ## Run 1 stress 0.06402617 
    ## Run 2 stress 0.06463313 
    ## Run 3 stress 0.06402621 
    ## Run 4 stress 0.06463313 
    ## Run 5 stress 0.08327913 
    ## Run 6 stress 0.06330868 
    ## ... Procrustes: rmse 1.010558e-05  max resid 2.388059e-05 
    ## ... Similar to previous best
    ## Run 7 stress 0.06330869 
    ## ... Procrustes: rmse 3.091477e-05  max resid 0.0001155194 
    ## ... Similar to previous best
    ## Run 8 stress 0.0640263 
    ## Run 9 stress 0.06463314 
    ## Run 10 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 9.821587e-06  max resid 2.177912e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.06402629 
    ## Run 12 stress 0.06402627 
    ## Run 13 stress 0.09999154 
    ## Run 14 stress 0.06330868 
    ## ... Procrustes: rmse 8.391589e-06  max resid 2.568518e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.06463313 
    ## Run 16 stress 0.06463314 
    ## Run 17 stress 0.06402622 
    ## Run 18 stress 0.06463313 
    ## Run 19 stress 0.06330869 
    ## ... Procrustes: rmse 3.740374e-05  max resid 8.667861e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.06330869 
    ## ... Procrustes: rmse 2.042598e-05  max resid 7.148754e-05 
    ## ... Similar to previous best
    ## Run 21 stress 0.06330868 
    ## ... Procrustes: rmse 1.645613e-05  max resid 5.080902e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.0633087 
    ## ... Procrustes: rmse 6.283126e-05  max resid 0.0001383048 
    ## ... Similar to previous best
    ## Run 23 stress 0.06402617 
    ## Run 24 stress 0.08301941 
    ## Run 25 stress 0.06402626 
    ## Run 26 stress 0.06463316 
    ## Run 27 stress 0.0646332 
    ## Run 28 stress 0.06330868 
    ## ... Procrustes: rmse 7.056695e-06  max resid 1.456546e-05 
    ## ... Similar to previous best
    ## Run 29 stress 0.06402633 
    ## Run 30 stress 0.06330869 
    ## ... Procrustes: rmse 2.927378e-05  max resid 6.181713e-05 
    ## ... Similar to previous best
    ## Run 31 stress 0.06463313 
    ## Run 32 stress 0.06463313 
    ## Run 33 stress 0.06402623 
    ## Run 34 stress 0.09336659 
    ## Run 35 stress 0.06330869 
    ## ... Procrustes: rmse 3.758199e-05  max resid 6.787684e-05 
    ## ... Similar to previous best
    ## Run 36 stress 0.06463315 
    ## Run 37 stress 0.06463314 
    ## Run 38 stress 0.06463313 
    ## Run 39 stress 0.06330869 
    ## ... Procrustes: rmse 5.229681e-05  max resid 0.000114473 
    ## ... Similar to previous best
    ## Run 40 stress 0.06402622 
    ## Run 41 stress 0.06402629 
    ## Run 42 stress 0.06463313 
    ## Run 43 stress 0.06402625 
    ## Run 44 stress 0.0633087 
    ## ... Procrustes: rmse 3.630362e-05  max resid 6.764138e-05 
    ## ... Similar to previous best
    ## Run 45 stress 0.09462076 
    ## Run 46 stress 0.0633087 
    ## ... Procrustes: rmse 6.436226e-05  max resid 0.0001415645 
    ## ... Similar to previous best
    ## Run 47 stress 0.0831563 
    ## Run 48 stress 0.06402616 
    ## Run 49 stress 0.06463313 
    ## Run 50 stress 0.06463314 
    ## Run 51 stress 0.06330868 
    ## ... New best solution
    ## ... Procrustes: rmse 9.460554e-06  max resid 3.163921e-05 
    ## ... Similar to previous best
    ## Run 52 stress 0.06330868 
    ## ... Procrustes: rmse 2.403156e-05  max resid 5.268242e-05 
    ## ... Similar to previous best
    ## Run 53 stress 0.06402623 
    ## Run 54 stress 0.06463313 
    ## Run 55 stress 0.08357051 
    ## Run 56 stress 0.09586232 
    ## Run 57 stress 0.0640262 
    ## Run 58 stress 0.06463313 
    ## Run 59 stress 0.08357049 
    ## Run 60 stress 0.06402629 
    ## Run 61 stress 0.06330869 
    ## ... Procrustes: rmse 2.674138e-05  max resid 6.162299e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.06463315 
    ## Run 63 stress 0.06330872 
    ## ... Procrustes: rmse 7.656354e-05  max resid 0.000172016 
    ## ... Similar to previous best
    ## Run 64 stress 0.06402624 
    ## Run 65 stress 0.06330868 
    ## ... Procrustes: rmse 6.877314e-06  max resid 2.554667e-05 
    ## ... Similar to previous best
    ## Run 66 stress 0.06330872 
    ## ... Procrustes: rmse 8.271698e-05  max resid 0.0001832738 
    ## ... Similar to previous best
    ## Run 67 stress 0.06330868 
    ## ... Procrustes: rmse 2.34307e-05  max resid 5.047152e-05 
    ## ... Similar to previous best
    ## Run 68 stress 0.0640263 
    ## Run 69 stress 0.08431269 
    ## Run 70 stress 0.06402618 
    ## Run 71 stress 0.06330871 
    ## ... Procrustes: rmse 1.740912e-05  max resid 5.245817e-05 
    ## ... Similar to previous best
    ## Run 72 stress 0.06463317 
    ## Run 73 stress 0.06402629 
    ## Run 74 stress 0.06463314 
    ## Run 75 stress 0.09290986 
    ## Run 76 stress 0.09846453 
    ## Run 77 stress 0.06330868 
    ## ... Procrustes: rmse 2.382594e-05  max resid 8.771328e-05 
    ## ... Similar to previous best
    ## Run 78 stress 0.09290996 
    ## Run 79 stress 0.08390026 
    ## Run 80 stress 0.0950462 
    ## Run 81 stress 0.06330869 
    ## ... Procrustes: rmse 4.120827e-05  max resid 9.082139e-05 
    ## ... Similar to previous best
    ## Run 82 stress 0.06463315 
    ## Run 83 stress 0.09233467 
    ## Run 84 stress 0.06463313 
    ## Run 85 stress 0.09378692 
    ## Run 86 stress 0.06463314 
    ## Run 87 stress 0.06330868 
    ## ... Procrustes: rmse 1.658506e-05  max resid 3.721064e-05 
    ## ... Similar to previous best
    ## Run 88 stress 0.06402629 
    ## Run 89 stress 0.06402624 
    ## Run 90 stress 0.06463321 
    ## Run 91 stress 0.0633087 
    ## ... Procrustes: rmse 4.363275e-05  max resid 0.0001376906 
    ## ... Similar to previous best
    ## Run 92 stress 0.06463315 
    ## Run 93 stress 0.06330869 
    ## ... Procrustes: rmse 2.898486e-05  max resid 7.430954e-05 
    ## ... Similar to previous best
    ## Run 94 stress 0.09283194 
    ## Run 95 stress 0.08089077 
    ## Run 96 stress 0.06402625 
    ## Run 97 stress 0.06463314 
    ## Run 98 stress 0.09157938 
    ## Run 99 stress 0.06463317 
    ## Run 100 stress 0.06402628 
    ## Run 101 stress 0.06463314 
    ## Run 102 stress 0.08126172 
    ## Run 103 stress 0.0835705 
    ## Run 104 stress 0.06402619 
    ## Run 105 stress 0.06330868 
    ## ... Procrustes: rmse 4.915667e-06  max resid 1.040049e-05 
    ## ... Similar to previous best
    ## Run 106 stress 0.06330868 
    ## ... Procrustes: rmse 1.174118e-05  max resid 2.472941e-05 
    ## ... Similar to previous best
    ## Run 107 stress 0.0633087 
    ## ... Procrustes: rmse 5.015094e-05  max resid 0.0001181291 
    ## ... Similar to previous best
    ## Run 108 stress 0.08089059 
    ## Run 109 stress 0.06402634 
    ## Run 110 stress 0.06463313 
    ## Run 111 stress 0.08233728 
    ## Run 112 stress 0.06330876 
    ## ... Procrustes: rmse 0.0001171245  max resid 0.0002640132 
    ## ... Similar to previous best
    ## Run 113 stress 0.06463313 
    ## Run 114 stress 0.06330868 
    ## ... Procrustes: rmse 4.163089e-06  max resid 9.557367e-06 
    ## ... Similar to previous best
    ## Run 115 stress 0.06402619 
    ## Run 116 stress 0.06402619 
    ## Run 117 stress 0.06330868 
    ## ... Procrustes: rmse 1.693297e-05  max resid 6.209196e-05 
    ## ... Similar to previous best
    ## Run 118 stress 0.06330869 
    ## ... Procrustes: rmse 2.920577e-05  max resid 7.833404e-05 
    ## ... Similar to previous best
    ## Run 119 stress 0.08246351 
    ## Run 120 stress 0.08089071 
    ## Run 121 stress 0.06463315 
    ## Run 122 stress 0.06463315 
    ## Run 123 stress 0.06330877 
    ## ... Procrustes: rmse 9.865965e-05  max resid 0.000220479 
    ## ... Similar to previous best
    ## Run 124 stress 0.06463313 
    ## Run 125 stress 0.06402636 
    ## Run 126 stress 0.08126195 
    ## Run 127 stress 0.08089059 
    ## Run 128 stress 0.0646332 
    ## Run 129 stress 0.09157936 
    ## Run 130 stress 0.06402615 
    ## Run 131 stress 0.06402632 
    ## Run 132 stress 0.08357054 
    ## Run 133 stress 0.06402631 
    ## Run 134 stress 0.09472712 
    ## Run 135 stress 0.06402624 
    ## Run 136 stress 0.09602322 
    ## Run 137 stress 0.06463313 
    ## Run 138 stress 0.0633087 
    ## ... Procrustes: rmse 5.696487e-05  max resid 0.0001269494 
    ## ... Similar to previous best
    ## Run 139 stress 0.06402615 
    ## Run 140 stress 0.06402638 
    ## Run 141 stress 0.0640262 
    ## Run 142 stress 0.09563688 
    ## Run 143 stress 0.06402619 
    ## Run 144 stress 0.06330868 
    ## ... Procrustes: rmse 4.61355e-06  max resid 1.342349e-05 
    ## ... Similar to previous best
    ## Run 145 stress 0.06402624 
    ## Run 146 stress 0.06330874 
    ## ... Procrustes: rmse 0.0001026815  max resid 0.0002311324 
    ## ... Similar to previous best
    ## Run 147 stress 0.06402633 
    ## Run 148 stress 0.0831561 
    ## Run 149 stress 0.08233729 
    ## Run 150 stress 0.08126181 
    ## Run 151 stress 0.06330868 
    ## ... Procrustes: rmse 2.263552e-05  max resid 4.378061e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.06463313 
    ## Run 153 stress 0.06463314 
    ## Run 154 stress 0.08315603 
    ## Run 155 stress 0.06330875 
    ## ... Procrustes: rmse 0.0001042408  max resid 0.0002780299 
    ## ... Similar to previous best
    ## Run 156 stress 0.06330871 
    ## ... Procrustes: rmse 6.932707e-05  max resid 0.0001534149 
    ## ... Similar to previous best
    ## Run 157 stress 0.06402618 
    ## Run 158 stress 0.06402639 
    ## Run 159 stress 0.06330869 
    ## ... Procrustes: rmse 3.199586e-05  max resid 7.1789e-05 
    ## ... Similar to previous best
    ## Run 160 stress 0.06463313 
    ## Run 161 stress 0.06463313 
    ## Run 162 stress 0.06330868 
    ## ... Procrustes: rmse 4.100244e-06  max resid 8.104517e-06 
    ## ... Similar to previous best
    ## Run 163 stress 0.06402616 
    ## Run 164 stress 0.06402615 
    ## Run 165 stress 0.06402626 
    ## Run 166 stress 0.06330875 
    ## ... Procrustes: rmse 0.000108194  max resid 0.0002391268 
    ## ... Similar to previous best
    ## Run 167 stress 0.06330871 
    ## ... Procrustes: rmse 6.992182e-05  max resid 0.0001555528 
    ## ... Similar to previous best
    ## Run 168 stress 0.06402632 
    ## Run 169 stress 0.06463313 
    ## Run 170 stress 0.06402616 
    ## Run 171 stress 0.06463314 
    ## Run 172 stress 0.06330868 
    ## ... Procrustes: rmse 1.264784e-05  max resid 4.67249e-05 
    ## ... Similar to previous best
    ## Run 173 stress 0.06330868 
    ## ... Procrustes: rmse 1.908367e-05  max resid 3.925945e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.06463314 
    ## Run 175 stress 0.08315605 
    ## Run 176 stress 0.08431258 
    ## Run 177 stress 0.06463313 
    ## Run 178 stress 0.06402619 
    ## Run 179 stress 0.09472738 
    ## Run 180 stress 0.06463313 
    ## Run 181 stress 0.06463313 
    ## Run 182 stress 0.06330869 
    ## ... Procrustes: rmse 4.12874e-05  max resid 9.276242e-05 
    ## ... Similar to previous best
    ## Run 183 stress 0.06463313 
    ## Run 184 stress 0.06402626 
    ## Run 185 stress 0.06330868 
    ## ... Procrustes: rmse 1.242986e-05  max resid 2.168277e-05 
    ## ... Similar to previous best
    ## Run 186 stress 0.0633087 
    ## ... Procrustes: rmse 1.189591e-05  max resid 3.107904e-05 
    ## ... Similar to previous best
    ## Run 187 stress 0.06330868 
    ## ... Procrustes: rmse 2.412919e-05  max resid 5.236595e-05 
    ## ... Similar to previous best
    ## Run 188 stress 0.06402631 
    ## Run 189 stress 0.06330868 
    ## ... Procrustes: rmse 5.63301e-06  max resid 9.472366e-06 
    ## ... Similar to previous best
    ## Run 190 stress 0.06463316 
    ## Run 191 stress 0.06330874 
    ## ... Procrustes: rmse 8.703415e-05  max resid 0.000196045 
    ## ... Similar to previous best
    ## Run 192 stress 0.06330868 
    ## ... Procrustes: rmse 1.353193e-05  max resid 2.893164e-05 
    ## ... Similar to previous best
    ## Run 193 stress 0.06330868 
    ## ... Procrustes: rmse 1.046994e-05  max resid 2.287458e-05 
    ## ... Similar to previous best
    ## Run 194 stress 0.06330875 
    ## ... Procrustes: rmse 0.0001041993  max resid 0.0002334922 
    ## ... Similar to previous best
    ## Run 195 stress 0.06463313 
    ## Run 196 stress 0.06330868 
    ## ... Procrustes: rmse 9.90522e-06  max resid 3.102219e-05 
    ## ... Similar to previous best
    ## Run 197 stress 0.06463313 
    ## Run 198 stress 0.08357055 
    ## Run 199 stress 0.06463317 
    ## Run 200 stress 0.06463313 
    ## Run 201 stress 0.06463313 
    ## Run 202 stress 0.06463315 
    ## Run 203 stress 0.06330868 
    ## ... Procrustes: rmse 1.15239e-05  max resid 4.126068e-05 
    ## ... Similar to previous best
    ## Run 204 stress 0.06402633 
    ## Run 205 stress 0.09157936 
    ## Run 206 stress 0.06402621 
    ## Run 207 stress 0.08126169 
    ## Run 208 stress 0.06463315 
    ## Run 209 stress 0.0640263 
    ## Run 210 stress 0.0633087 
    ## ... Procrustes: rmse 5.259387e-05  max resid 0.0001665794 
    ## ... Similar to previous best
    ## Run 211 stress 0.0640262 
    ## Run 212 stress 0.08357066 
    ## Run 213 stress 0.06463313 
    ## Run 214 stress 0.08357056 
    ## Run 215 stress 0.06402618 
    ## Run 216 stress 0.06402627 
    ## Run 217 stress 0.0633087 
    ## ... Procrustes: rmse 4.933445e-05  max resid 0.000108603 
    ## ... Similar to previous best
    ## Run 218 stress 0.06463315 
    ## Run 219 stress 0.09106927 
    ## Run 220 stress 0.06330877 
    ## ... Procrustes: rmse 0.0001201579  max resid 0.0002695643 
    ## ... Similar to previous best
    ## Run 221 stress 0.06463315 
    ## Run 222 stress 0.06463317 
    ## Run 223 stress 0.06330871 
    ## ... Procrustes: rmse 6.37556e-05  max resid 0.0001091791 
    ## ... Similar to previous best
    ## Run 224 stress 0.06330874 
    ## ... Procrustes: rmse 9.90397e-05  max resid 0.0002161295 
    ## ... Similar to previous best
    ## Run 225 stress 0.06330868 
    ## ... Procrustes: rmse 1.057154e-05  max resid 2.317275e-05 
    ## ... Similar to previous best
    ## Run 226 stress 0.06402621 
    ## Run 227 stress 0.06330869 
    ## ... Procrustes: rmse 2.730841e-05  max resid 6.380842e-05 
    ## ... Similar to previous best
    ## Run 228 stress 0.0633087 
    ## ... Procrustes: rmse 5.098214e-05  max resid 0.000111619 
    ## ... Similar to previous best
    ## Run 229 stress 0.06330868 
    ## ... Procrustes: rmse 1.046133e-05  max resid 2.371143e-05 
    ## ... Similar to previous best
    ## Run 230 stress 0.06330873 
    ## ... Procrustes: rmse 8.643007e-05  max resid 0.0001899605 
    ## ... Similar to previous best
    ## Run 231 stress 0.08246376 
    ## Run 232 stress 0.08315619 
    ## Run 233 stress 0.08344966 
    ## Run 234 stress 0.0633087 
    ## ... Procrustes: rmse 5.459794e-05  max resid 0.0001939508 
    ## ... Similar to previous best
    ## Run 235 stress 0.06402639 
    ## Run 236 stress 0.0640262 
    ## Run 237 stress 0.06402627 
    ## Run 238 stress 0.0633087 
    ## ... Procrustes: rmse 4.38209e-05  max resid 9.724061e-05 
    ## ... Similar to previous best
    ## Run 239 stress 0.06402618 
    ## Run 240 stress 0.06463315 
    ## Run 241 stress 0.06330869 
    ## ... Procrustes: rmse 4.024452e-05  max resid 0.0001445002 
    ## ... Similar to previous best
    ## Run 242 stress 0.09270079 
    ## Run 243 stress 0.08089076 
    ## Run 244 stress 0.09504629 
    ## Run 245 stress 0.06402636 
    ## Run 246 stress 0.06463316 
    ## Run 247 stress 0.06463314 
    ## Run 248 stress 0.06402629 
    ## Run 249 stress 0.06463319 
    ## Run 250 stress 0.06402634 
    ## Run 251 stress 0.08089038 
    ## Run 252 stress 0.06330869 
    ## ... Procrustes: rmse 3.189797e-05  max resid 5.861923e-05 
    ## ... Similar to previous best
    ## Run 253 stress 0.06463313 
    ## Run 254 stress 0.06402626 
    ## Run 255 stress 0.08431267 
    ## Run 256 stress 0.06330868 
    ## ... Procrustes: rmse 4.458485e-06  max resid 8.295033e-06 
    ## ... Similar to previous best
    ## Run 257 stress 0.06330868 
    ## ... Procrustes: rmse 5.977606e-06  max resid 1.165187e-05 
    ## ... Similar to previous best
    ## Run 258 stress 0.06330868 
    ## ... Procrustes: rmse 9.876881e-06  max resid 2.264244e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.09472844 
    ## Run 260 stress 0.06402619 
    ## Run 261 stress 0.06463313 
    ## Run 262 stress 0.06330868 
    ## ... Procrustes: rmse 2.187416e-05  max resid 3.796212e-05 
    ## ... Similar to previous best
    ## Run 263 stress 0.06402617 
    ## Run 264 stress 0.06402619 
    ## Run 265 stress 0.06330868 
    ## ... Procrustes: rmse 3.183277e-06  max resid 5.287031e-06 
    ## ... Similar to previous best
    ## Run 266 stress 0.06330868 
    ## ... Procrustes: rmse 5.170431e-06  max resid 1.181778e-05 
    ## ... Similar to previous best
    ## Run 267 stress 0.08126211 
    ## Run 268 stress 0.06463315 
    ## Run 269 stress 0.06463314 
    ## Run 270 stress 0.06402618 
    ## Run 271 stress 0.06402625 
    ## Run 272 stress 0.06463313 
    ## Run 273 stress 0.09542741 
    ## Run 274 stress 0.06402619 
    ## Run 275 stress 0.06402636 
    ## Run 276 stress 0.06463314 
    ## Run 277 stress 0.06330868 
    ## ... Procrustes: rmse 5.790258e-06  max resid 1.836212e-05 
    ## ... Similar to previous best
    ## Run 278 stress 0.08089056 
    ## Run 279 stress 0.06402627 
    ## Run 280 stress 0.09602316 
    ## Run 281 stress 0.06330868 
    ## ... Procrustes: rmse 1.611491e-05  max resid 5.82196e-05 
    ## ... Similar to previous best
    ## Run 282 stress 0.09547937 
    ## Run 283 stress 0.06330869 
    ## ... Procrustes: rmse 1.326777e-05  max resid 3.321878e-05 
    ## ... Similar to previous best
    ## Run 284 stress 0.0633087 
    ## ... Procrustes: rmse 6.276402e-05  max resid 0.0001363841 
    ## ... Similar to previous best
    ## Run 285 stress 0.09550955 
    ## Run 286 stress 0.06463314 
    ## Run 287 stress 0.08089065 
    ## Run 288 stress 0.06463314 
    ## Run 289 stress 0.06402619 
    ## Run 290 stress 0.06463313 
    ## Run 291 stress 0.06402615 
    ## Run 292 stress 0.06330872 
    ## ... Procrustes: rmse 8.582766e-05  max resid 0.0001923444 
    ## ... Similar to previous best
    ## Run 293 stress 0.08126186 
    ## Run 294 stress 0.06463314 
    ## Run 295 stress 0.06463313 
    ## Run 296 stress 0.06330869 
    ## ... Procrustes: rmse 3.467895e-05  max resid 7.66129e-05 
    ## ... Similar to previous best
    ## Run 297 stress 0.0633087 
    ## ... Procrustes: rmse 5.371774e-05  max resid 0.0001207391 
    ## ... Similar to previous best
    ## Run 298 stress 0.06330868 
    ## ... Procrustes: rmse 4.58003e-06  max resid 8.310217e-06 
    ## ... Similar to previous best
    ## Run 299 stress 0.06330868 
    ## ... Procrustes: rmse 8.656537e-06  max resid 1.769819e-05 
    ## ... Similar to previous best
    ## Run 300 stress 0.0646332 
    ## Run 301 stress 0.06402616 
    ## Run 302 stress 0.06463313 
    ## Run 303 stress 0.06330871 
    ## ... Procrustes: rmse 7.641186e-05  max resid 0.0001692653 
    ## ... Similar to previous best
    ## Run 304 stress 0.06330874 
    ## ... Procrustes: rmse 0.0001007956  max resid 0.0002205544 
    ## ... Similar to previous best
    ## Run 305 stress 0.06463315 
    ## Run 306 stress 0.06402627 
    ## Run 307 stress 0.0633087 
    ## ... Procrustes: rmse 5.38122e-05  max resid 0.0001203093 
    ## ... Similar to previous best
    ## Run 308 stress 0.0633087 
    ## ... Procrustes: rmse 5.665959e-05  max resid 0.0001269748 
    ## ... Similar to previous best
    ## Run 309 stress 0.0812619 
    ## Run 310 stress 0.0633087 
    ## ... Procrustes: rmse 5.474168e-05  max resid 0.0001231374 
    ## ... Similar to previous best
    ## Run 311 stress 0.06402629 
    ## Run 312 stress 0.06330871 
    ## ... Procrustes: rmse 2.41349e-05  max resid 4.787644e-05 
    ## ... Similar to previous best
    ## Run 313 stress 0.06463313 
    ## Run 314 stress 0.06463313 
    ## Run 315 stress 0.09485808 
    ## Run 316 stress 0.06463313 
    ## Run 317 stress 0.06463318 
    ## Run 318 stress 0.06402635 
    ## Run 319 stress 0.09523979 
    ## Run 320 stress 0.0633088 
    ## ... Procrustes: rmse 0.0001421951  max resid 0.0003196495 
    ## ... Similar to previous best
    ## Run 321 stress 0.06330877 
    ## ... Procrustes: rmse 0.0001162953  max resid 0.0002581425 
    ## ... Similar to previous best
    ## Run 322 stress 0.06402625 
    ## Run 323 stress 0.06330872 
    ## ... Procrustes: rmse 3.438142e-05  max resid 6.500856e-05 
    ## ... Similar to previous best
    ## Run 324 stress 0.06330868 
    ## ... Procrustes: rmse 8.542752e-06  max resid 2.398441e-05 
    ## ... Similar to previous best
    ## Run 325 stress 0.06463317 
    ## Run 326 stress 0.09233437 
    ## Run 327 stress 0.06330876 
    ## ... Procrustes: rmse 0.0001177423  max resid 0.0002664555 
    ## ... Similar to previous best
    ## Run 328 stress 0.06330879 
    ## ... Procrustes: rmse 0.0001404991  max resid 0.0003173534 
    ## ... Similar to previous best
    ## Run 329 stress 0.06402616 
    ## Run 330 stress 0.08389954 
    ## Run 331 stress 0.06402629 
    ## Run 332 stress 0.08089039 
    ## Run 333 stress 0.06330876 
    ## ... Procrustes: rmse 0.0001135083  max resid 0.0002497478 
    ## ... Similar to previous best
    ## Run 334 stress 0.06402615 
    ## Run 335 stress 0.06330868 
    ## ... Procrustes: rmse 1.37862e-05  max resid 4.985483e-05 
    ## ... Similar to previous best
    ## Run 336 stress 0.06402631 
    ## Run 337 stress 0.09157919 
    ## Run 338 stress 0.06402641 
    ## Run 339 stress 0.09485786 
    ## Run 340 stress 0.06330872 
    ## ... Procrustes: rmse 8.396479e-05  max resid 0.0001895114 
    ## ... Similar to previous best
    ## Run 341 stress 0.06330869 
    ## ... Procrustes: rmse 2.743344e-05  max resid 6.168195e-05 
    ## ... Similar to previous best
    ## Run 342 stress 0.0835706 
    ## Run 343 stress 0.06402621 
    ## Run 344 stress 0.08126203 
    ## Run 345 stress 0.08126206 
    ## Run 346 stress 0.06463313 
    ## Run 347 stress 0.06402621 
    ## Run 348 stress 0.06463315 
    ## Run 349 stress 0.06402617 
    ## Run 350 stress 0.06330872 
    ## ... Procrustes: rmse 8.704761e-05  max resid 0.0001988149 
    ## ... Similar to previous best
    ## Run 351 stress 0.06402618 
    ## Run 352 stress 0.06463313 
    ## Run 353 stress 0.06402632 
    ## Run 354 stress 0.06463314 
    ## Run 355 stress 0.06330868 
    ## ... Procrustes: rmse 4.286022e-06  max resid 8.771346e-06 
    ## ... Similar to previous best
    ## Run 356 stress 0.0640262 
    ## Run 357 stress 0.09528427 
    ## Run 358 stress 0.06330877 
    ## ... Procrustes: rmse 0.0001231042  max resid 0.000277029 
    ## ... Similar to previous best
    ## Run 359 stress 0.09106927 
    ## Run 360 stress 0.06330869 
    ## ... Procrustes: rmse 3.088178e-05  max resid 8.60778e-05 
    ## ... Similar to previous best
    ## Run 361 stress 0.06463313 
    ## Run 362 stress 0.06330869 
    ## ... Procrustes: rmse 3.99969e-05  max resid 8.722616e-05 
    ## ... Similar to previous best
    ## Run 363 stress 0.06330871 
    ## ... Procrustes: rmse 7.272244e-05  max resid 0.0001634352 
    ## ... Similar to previous best
    ## Run 364 stress 0.08459962 
    ## Run 365 stress 0.06330868 
    ## ... Procrustes: rmse 9.52912e-06  max resid 2.389828e-05 
    ## ... Similar to previous best
    ## Run 366 stress 0.06463313 
    ## Run 367 stress 0.06330868 
    ## ... Procrustes: rmse 2.23532e-05  max resid 3.806113e-05 
    ## ... Similar to previous best
    ## Run 368 stress 0.06463319 
    ## Run 369 stress 0.0633087 
    ## ... Procrustes: rmse 4.14397e-05  max resid 9.004203e-05 
    ## ... Similar to previous best
    ## Run 370 stress 0.0633087 
    ## ... Procrustes: rmse 4.779669e-05  max resid 8.325044e-05 
    ## ... Similar to previous best
    ## Run 371 stress 0.06402633 
    ## Run 372 stress 0.06463314 
    ## Run 373 stress 0.06463313 
    ## Run 374 stress 0.08301933 
    ## Run 375 stress 0.06402617 
    ## Run 376 stress 0.0808905 
    ## Run 377 stress 0.08315608 
    ## Run 378 stress 0.06330868 
    ## ... Procrustes: rmse 5.232496e-06  max resid 1.447794e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.06402621 
    ## Run 380 stress 0.06463313 
    ## Run 381 stress 0.0640263 
    ## Run 382 stress 0.06463313 
    ## Run 383 stress 0.06402616 
    ## Run 384 stress 0.06402635 
    ## Run 385 stress 0.06402616 
    ## Run 386 stress 0.06463317 
    ## Run 387 stress 0.06330868 
    ## ... Procrustes: rmse 1.759737e-05  max resid 3.825126e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.06463313 
    ## Run 389 stress 0.06330872 
    ## ... Procrustes: rmse 8.709951e-05  max resid 0.0001904589 
    ## ... Similar to previous best
    ## Run 390 stress 0.06330869 
    ## ... Procrustes: rmse 3.095344e-05  max resid 6.378352e-05 
    ## ... Similar to previous best
    ## Run 391 stress 0.06330871 
    ## ... Procrustes: rmse 4.750558e-05  max resid 0.0001025166 
    ## ... Similar to previous best
    ## Run 392 stress 0.06330868 
    ## ... Procrustes: rmse 1.406122e-05  max resid 3.929042e-05 
    ## ... Similar to previous best
    ## Run 393 stress 0.06330874 
    ## ... Procrustes: rmse 0.000100104  max resid 0.0002173843 
    ## ... Similar to previous best
    ## Run 394 stress 0.06463315 
    ## Run 395 stress 0.09596594 
    ## Run 396 stress 0.06330868 
    ## ... Procrustes: rmse 1.076512e-05  max resid 2.367795e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.06330869 
    ## ... Procrustes: rmse 2.968044e-05  max resid 7.547781e-05 
    ## ... Similar to previous best
    ## Run 398 stress 0.09217183 
    ## Run 399 stress 0.06330868 
    ## ... Procrustes: rmse 1.091441e-05  max resid 2.477921e-05 
    ## ... Similar to previous best
    ## Run 400 stress 0.06330868 
    ## ... Procrustes: rmse 7.270838e-06  max resid 1.79327e-05 
    ## ... Similar to previous best
    ## Run 401 stress 0.06463315 
    ## Run 402 stress 0.06402619 
    ## Run 403 stress 0.06463313 
    ## Run 404 stress 0.06463315 
    ## Run 405 stress 0.06402632 
    ## Run 406 stress 0.06402634 
    ## Run 407 stress 0.06402617 
    ## Run 408 stress 0.08301932 
    ## Run 409 stress 0.06463313 
    ## Run 410 stress 0.06330869 
    ## ... Procrustes: rmse 2.943857e-05  max resid 6.652661e-05 
    ## ... Similar to previous best
    ## Run 411 stress 0.06402616 
    ## Run 412 stress 0.06463315 
    ## Run 413 stress 0.0633087 
    ## ... Procrustes: rmse 5.905577e-05  max resid 0.0001282726 
    ## ... Similar to previous best
    ## Run 414 stress 0.06463318 
    ## Run 415 stress 0.08089078 
    ## Run 416 stress 0.06330871 
    ## ... Procrustes: rmse 6.651823e-05  max resid 0.0001493458 
    ## ... Similar to previous best
    ## Run 417 stress 0.06463318 
    ## Run 418 stress 0.0633087 
    ## ... Procrustes: rmse 5.191157e-05  max resid 0.0001148143 
    ## ... Similar to previous best
    ## Run 419 stress 0.08390003 
    ## Run 420 stress 0.06463313 
    ## Run 421 stress 0.06330868 
    ## ... Procrustes: rmse 7.123054e-06  max resid 2.556925e-05 
    ## ... Similar to previous best
    ## Run 422 stress 0.09489759 
    ## Run 423 stress 0.0640263 
    ## Run 424 stress 0.06402618 
    ## Run 425 stress 0.06330871 
    ## ... Procrustes: rmse 6.6376e-05  max resid 0.000150176 
    ## ... Similar to previous best
    ## Run 426 stress 0.06330869 
    ## ... Procrustes: rmse 2.514981e-05  max resid 5.214501e-05 
    ## ... Similar to previous best
    ## Run 427 stress 0.06402616 
    ## Run 428 stress 0.06402618 
    ## Run 429 stress 0.06463313 
    ## Run 430 stress 0.06330869 
    ## ... Procrustes: rmse 3.465855e-05  max resid 7.495364e-05 
    ## ... Similar to previous best
    ## Run 431 stress 0.08390015 
    ## Run 432 stress 0.06463314 
    ## Run 433 stress 0.06402624 
    ## Run 434 stress 0.06463313 
    ## Run 435 stress 0.06402628 
    ## Run 436 stress 0.09283161 
    ## Run 437 stress 0.06330872 
    ## ... Procrustes: rmse 8.107019e-05  max resid 0.0001839207 
    ## ... Similar to previous best
    ## Run 438 stress 0.09283145 
    ## Run 439 stress 0.06463313 
    ## Run 440 stress 0.06463313 
    ## Run 441 stress 0.06402616 
    ## Run 442 stress 0.06463314 
    ## Run 443 stress 0.06402616 
    ## Run 444 stress 0.06402628 
    ## Run 445 stress 0.06402634 
    ## Run 446 stress 0.06402639 
    ## Run 447 stress 0.08089071 
    ## Run 448 stress 0.09954592 
    ## Run 449 stress 0.06402615 
    ## Run 450 stress 0.06330871 
    ## ... Procrustes: rmse 6.375277e-05  max resid 0.0002183839 
    ## ... Similar to previous best
    ## Run 451 stress 0.0633087 
    ## ... Procrustes: rmse 2.239672e-05  max resid 5.354305e-05 
    ## ... Similar to previous best
    ## Run 452 stress 0.06402621 
    ## Run 453 stress 0.09503234 
    ## Run 454 stress 0.06463315 
    ## Run 455 stress 0.08344973 
    ## Run 456 stress 0.06463319 
    ## Run 457 stress 0.09602274 
    ## Run 458 stress 0.06463313 
    ## Run 459 stress 0.09233453 
    ## Run 460 stress 0.06402623 
    ## Run 461 stress 0.06402617 
    ## Run 462 stress 0.08315612 
    ## Run 463 stress 0.06402617 
    ## Run 464 stress 0.08089054 
    ## Run 465 stress 0.09682934 
    ## Run 466 stress 0.09245945 
    ## Run 467 stress 0.09935143 
    ## Run 468 stress 0.06463313 
    ## Run 469 stress 0.09875994 
    ## Run 470 stress 0.06463314 
    ## Run 471 stress 0.0640263 
    ## Run 472 stress 0.06330869 
    ## ... Procrustes: rmse 4.910002e-05  max resid 0.0001133442 
    ## ... Similar to previous best
    ## Run 473 stress 0.06463313 
    ## Run 474 stress 0.06402634 
    ## Run 475 stress 0.06330868 
    ## ... Procrustes: rmse 6.240848e-06  max resid 1.459028e-05 
    ## ... Similar to previous best
    ## Run 476 stress 0.06330868 
    ## ... Procrustes: rmse 5.101982e-06  max resid 1.055187e-05 
    ## ... Similar to previous best
    ## Run 477 stress 0.06463316 
    ## Run 478 stress 0.06402621 
    ## Run 479 stress 0.06330871 
    ## ... Procrustes: rmse 6.256483e-05  max resid 0.0001683541 
    ## ... Similar to previous best
    ## Run 480 stress 0.06402625 
    ## Run 481 stress 0.08315608 
    ## Run 482 stress 0.06463314 
    ## Run 483 stress 0.06330868 
    ## ... Procrustes: rmse 7.504801e-06  max resid 2.287153e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.08315607 
    ## Run 485 stress 0.06463315 
    ## Run 486 stress 0.06463316 
    ## Run 487 stress 0.06463315 
    ## Run 488 stress 0.06463314 
    ## Run 489 stress 0.06330881 
    ## ... Procrustes: rmse 0.0001293382  max resid 0.0002794448 
    ## ... Similar to previous best
    ## Run 490 stress 0.06402633 
    ## Run 491 stress 0.06402617 
    ## Run 492 stress 0.06330869 
    ## ... Procrustes: rmse 4.922496e-05  max resid 0.0001118101 
    ## ... Similar to previous best
    ## Run 493 stress 0.06402616 
    ## Run 494 stress 0.06463321 
    ## Run 495 stress 0.06463313 
    ## Run 496 stress 0.0633087 
    ## ... Procrustes: rmse 5.239171e-05  max resid 0.0001164754 
    ## ... Similar to previous best
    ## Run 497 stress 0.09157928 
    ## Run 498 stress 0.09232124 
    ## Run 499 stress 0.06402616 
    ## Run 500 stress 0.06330868 
    ## ... Procrustes: rmse 5.638259e-06  max resid 1.199751e-05 
    ## ... Similar to previous best
    ## *** Best solution repeated 130 times

``` r
# Mixed and stratified lakes
PD_beta_env_MS_NMDS <- metaMDS(PD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.06904483 
    ## Run 2 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001327767  max resid 0.0003494306 
    ## ... Similar to previous best
    ## Run 3 stress 0.3657354 
    ## Run 4 stress 0.05036799 
    ## ... Procrustes: rmse 0.000161441  max resid 0.0004464304 
    ## ... Similar to previous best
    ## Run 5 stress 0.05928227 
    ## Run 6 stress 0.06136676 
    ## Run 7 stress 0.3292803 
    ## Run 8 stress 0.05829934 
    ## Run 9 stress 0.05403582 
    ## Run 10 stress 0.06885089 
    ## Run 11 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001933439  max resid 0.0004595577 
    ## ... Similar to previous best
    ## Run 12 stress 0.05337577 
    ## Run 13 stress 0.07024642 
    ## Run 14 stress 0.05968626 
    ## Run 15 stress 0.05968637 
    ## Run 16 stress 0.050513 
    ## ... Procrustes: rmse 0.03558898  max resid 0.1216294 
    ## Run 17 stress 0.05337589 
    ## Run 18 stress 0.05051299 
    ## ... Procrustes: rmse 0.03569201  max resid 0.121835 
    ## Run 19 stress 0.06038956 
    ## Run 20 stress 0.05900164 
    ## Run 21 stress 0.05171059 
    ## Run 22 stress 0.05829931 
    ## Run 23 stress 0.05171059 
    ## Run 24 stress 0.05900149 
    ## Run 25 stress 0.05171064 
    ## Run 26 stress 0.05900152 
    ## Run 27 stress 0.06231291 
    ## Run 28 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001915347  max resid 0.0004773274 
    ## ... Similar to previous best
    ## Run 29 stress 0.06640376 
    ## Run 30 stress 0.05051343 
    ## ... Procrustes: rmse 0.03574724  max resid 0.1221289 
    ## Run 31 stress 0.06231287 
    ## Run 32 stress 0.05337577 
    ## Run 33 stress 0.0613668 
    ## Run 34 stress 0.05171059 
    ## Run 35 stress 0.05602064 
    ## Run 36 stress 0.05171059 
    ## Run 37 stress 0.06478516 
    ## Run 38 stress 0.05036794 
    ## ... Procrustes: rmse 2.850875e-05  max resid 7.317066e-05 
    ## ... Similar to previous best
    ## Run 39 stress 0.05051315 
    ## ... Procrustes: rmse 0.03573844  max resid 0.1219525 
    ## Run 40 stress 0.05051314 
    ## ... Procrustes: rmse 0.0356803  max resid 0.1217817 
    ## Run 41 stress 0.050368 
    ## ... Procrustes: rmse 0.000123011  max resid 0.0003581064 
    ## ... Similar to previous best
    ## Run 42 stress 0.050368 
    ## ... Procrustes: rmse 0.0001694786  max resid 0.000415179 
    ## ... Similar to previous best
    ## Run 43 stress 0.06674553 
    ## Run 44 stress 0.05036793 
    ## ... Procrustes: rmse 3.98807e-05  max resid 0.0001096897 
    ## ... Similar to previous best
    ## Run 45 stress 0.05036806 
    ## ... Procrustes: rmse 0.0002152487  max resid 0.0005352561 
    ## ... Similar to previous best
    ## Run 46 stress 0.05036794 
    ## ... Procrustes: rmse 4.770415e-05  max resid 0.0001297755 
    ## ... Similar to previous best
    ## Run 47 stress 0.05171061 
    ## Run 48 stress 0.05036798 
    ## ... Procrustes: rmse 0.000102452  max resid 0.0002993285 
    ## ... Similar to previous best
    ## Run 49 stress 0.05036794 
    ## ... Procrustes: rmse 0.0001146621  max resid 0.0002673194 
    ## ... Similar to previous best
    ## Run 50 stress 0.06111765 
    ## Run 51 stress 0.05403612 
    ## Run 52 stress 0.05036793 
    ## ... Procrustes: rmse 4.855525e-05  max resid 0.0001368078 
    ## ... Similar to previous best
    ## Run 53 stress 0.05829935 
    ## Run 54 stress 0.06136686 
    ## Run 55 stress 0.06231304 
    ## Run 56 stress 0.05051282 
    ## ... Procrustes: rmse 0.03569835  max resid 0.1219473 
    ## Run 57 stress 0.05051284 
    ## ... Procrustes: rmse 0.03563163  max resid 0.1217453 
    ## Run 58 stress 0.06309846 
    ## Run 59 stress 0.05036792 
    ## ... Procrustes: rmse 7.85131e-05  max resid 0.0001894709 
    ## ... Similar to previous best
    ## Run 60 stress 0.05466701 
    ## Run 61 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.228301e-05  max resid 2.070764e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.05602052 
    ## Run 63 stress 0.05602054 
    ## Run 64 stress 0.05171063 
    ## Run 65 stress 0.06826693 
    ## Run 66 stress 0.05928232 
    ## Run 67 stress 0.05829934 
    ## Run 68 stress 0.05337583 
    ## Run 69 stress 0.05171064 
    ## Run 70 stress 0.05829931 
    ## Run 71 stress 0.05051319 
    ## ... Procrustes: rmse 0.03570808  max resid 0.1218569 
    ## Run 72 stress 0.06136692 
    ## Run 73 stress 0.05403603 
    ## Run 74 stress 0.05337579 
    ## Run 75 stress 0.06309838 
    ## Run 76 stress 0.05403622 
    ## Run 77 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001199079  max resid 0.0002916919 
    ## ... Similar to previous best
    ## Run 78 stress 0.05051337 
    ## ... Procrustes: rmse 0.03574798  max resid 0.122129 
    ## Run 79 stress 0.06789734 
    ## Run 80 stress 0.05051332 
    ## ... Procrustes: rmse 0.03567784  max resid 0.1217544 
    ## Run 81 stress 0.05337584 
    ## Run 82 stress 0.05051294 
    ## ... Procrustes: rmse 0.03569435  max resid 0.121947 
    ## Run 83 stress 0.05968639 
    ## Run 84 stress 0.05051289 
    ## ... Procrustes: rmse 0.03562697  max resid 0.1217373 
    ## Run 85 stress 0.05900152 
    ## Run 86 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 3.041944e-05  max resid 5.91591e-05 
    ## ... Similar to previous best
    ## Run 87 stress 0.05036794 
    ## ... Procrustes: rmse 7.48451e-05  max resid 0.0001972302 
    ## ... Similar to previous best
    ## Run 88 stress 0.06231317 
    ## Run 89 stress 0.06111772 
    ## Run 90 stress 0.050513 
    ## ... Procrustes: rmse 0.0357327  max resid 0.121948 
    ## Run 91 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001249214  max resid 0.0003114404 
    ## ... Similar to previous best
    ## Run 92 stress 0.05968633 
    ## Run 93 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001278394  max resid 0.0003362206 
    ## ... Similar to previous best
    ## Run 94 stress 0.06476223 
    ## Run 95 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001274035  max resid 0.0003122609 
    ## ... Similar to previous best
    ## Run 96 stress 0.05829934 
    ## Run 97 stress 0.0663053 
    ## Run 98 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001201061  max resid 0.0003091627 
    ## ... Similar to previous best
    ## Run 99 stress 0.05403588 
    ## Run 100 stress 0.06478538 
    ## Run 101 stress 0.06499238 
    ## Run 102 stress 0.05036807 
    ## ... Procrustes: rmse 0.0002057535  max resid 0.0005387945 
    ## ... Similar to previous best
    ## Run 103 stress 0.05968626 
    ## Run 104 stress 0.06038944 
    ## Run 105 stress 0.0560205 
    ## Run 106 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001155686  max resid 0.0002912616 
    ## ... Similar to previous best
    ## Run 107 stress 0.05900157 
    ## Run 108 stress 0.05979229 
    ## Run 109 stress 0.05051348 
    ## ... Procrustes: rmse 0.03572386  max resid 0.1220559 
    ## Run 110 stress 0.06478545 
    ## Run 111 stress 0.05602057 
    ## Run 112 stress 0.05337578 
    ## Run 113 stress 0.05036794 
    ## ... Procrustes: rmse 9.853093e-05  max resid 0.0002416662 
    ## ... Similar to previous best
    ## Run 114 stress 0.05171059 
    ## Run 115 stress 0.05171059 
    ## Run 116 stress 0.06309849 
    ## Run 117 stress 0.06894592 
    ## Run 118 stress 0.0689686 
    ## Run 119 stress 0.0505136 
    ## ... Procrustes: rmse 0.03569908  max resid 0.1219857 
    ## Run 120 stress 0.06476232 
    ## Run 121 stress 0.05466704 
    ## Run 122 stress 0.06027437 
    ## Run 123 stress 0.05051336 
    ## ... Procrustes: rmse 0.03568687  max resid 0.1219417 
    ## Run 124 stress 0.05171058 
    ## Run 125 stress 0.05968629 
    ## Run 126 stress 0.05036793 
    ## ... Procrustes: rmse 8.596151e-05  max resid 0.0002030618 
    ## ... Similar to previous best
    ## Run 127 stress 0.05051313 
    ## ... Procrustes: rmse 0.0357131  max resid 0.1218716 
    ## Run 128 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 2.560223e-05  max resid 5.457808e-05 
    ## ... Similar to previous best
    ## Run 129 stress 0.05466711 
    ## Run 130 stress 0.05466702 
    ## Run 131 stress 0.06231306 
    ## Run 132 stress 0.05036791 
    ## ... Procrustes: rmse 2.069861e-05  max resid 4.352411e-05 
    ## ... Similar to previous best
    ## Run 133 stress 0.347363 
    ## Run 134 stress 0.06038952 
    ## Run 135 stress 0.063628 
    ## Run 136 stress 0.05466701 
    ## Run 137 stress 0.05337582 
    ## Run 138 stress 0.05403578 
    ## Run 139 stress 0.05968625 
    ## Run 140 stress 0.06136684 
    ## Run 141 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001170491  max resid 0.0003008902 
    ## ... Similar to previous best
    ## Run 142 stress 0.05171062 
    ## Run 143 stress 0.0582994 
    ## Run 144 stress 0.06362798 
    ## Run 145 stress 0.05466703 
    ## Run 146 stress 0.06136686 
    ## Run 147 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001009966  max resid 0.000242698 
    ## ... Similar to previous best
    ## Run 148 stress 0.05337577 
    ## Run 149 stress 0.05968625 
    ## Run 150 stress 0.050513 
    ## ... Procrustes: rmse 0.03570283  max resid 0.1218594 
    ## Run 151 stress 0.05900154 
    ## Run 152 stress 0.05036795 
    ## ... Procrustes: rmse 8.43141e-05  max resid 0.0001881495 
    ## ... Similar to previous best
    ## Run 153 stress 0.05171058 
    ## Run 154 stress 0.06647196 
    ## Run 155 stress 0.05337579 
    ## Run 156 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001134987  max resid 0.0003289512 
    ## ... Similar to previous best
    ## Run 157 stress 0.06789799 
    ## Run 158 stress 0.05171059 
    ## Run 159 stress 0.05051341 
    ## ... Procrustes: rmse 0.0357293  max resid 0.1220738 
    ## Run 160 stress 0.05051309 
    ## ... Procrustes: rmse 0.03570142  max resid 0.1218507 
    ## Run 161 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001485169  max resid 0.0003262507 
    ## ... Similar to previous best
    ## Run 162 stress 0.05602053 
    ## Run 163 stress 0.05051304 
    ## ... Procrustes: rmse 0.03569502  max resid 0.1219535 
    ## Run 164 stress 0.05051342 
    ## ... Procrustes: rmse 0.03573214  max resid 0.1220818 
    ## Run 165 stress 0.06478522 
    ## Run 166 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001665038  max resid 0.0004069251 
    ## ... Similar to previous best
    ## Run 167 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001591382  max resid 0.0003900304 
    ## ... Similar to previous best
    ## Run 168 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001035122  max resid 0.0002945333 
    ## ... Similar to previous best
    ## Run 169 stress 0.05051352 
    ## ... Procrustes: rmse 0.03569658  max resid 0.1217876 
    ## Run 170 stress 0.05337589 
    ## Run 171 stress 0.05403557 
    ## Run 172 stress 0.05337577 
    ## Run 173 stress 0.05337577 
    ## Run 174 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001372336  max resid 0.0003182223 
    ## ... Similar to previous best
    ## Run 175 stress 0.05466704 
    ## Run 176 stress 0.05051307 
    ## ... Procrustes: rmse 0.0356998  max resid 0.1218421 
    ## Run 177 stress 0.05829929 
    ## Run 178 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001022855  max resid 0.0002978372 
    ## ... Similar to previous best
    ## Run 179 stress 0.05403595 
    ## Run 180 stress 0.05051323 
    ## ... Procrustes: rmse 0.03569594  max resid 0.1218138 
    ## Run 181 stress 0.06309836 
    ## Run 182 stress 0.06789759 
    ## Run 183 stress 0.05051343 
    ## ... Procrustes: rmse 0.03570506  max resid 0.1220021 
    ## Run 184 stress 0.06826667 
    ## Run 185 stress 0.06362789 
    ## Run 186 stress 0.05171059 
    ## Run 187 stress 0.05171058 
    ## Run 188 stress 0.06309838 
    ## Run 189 stress 0.05051315 
    ## ... Procrustes: rmse 0.03566572  max resid 0.1217314 
    ## Run 190 stress 0.05968628 
    ## Run 191 stress 0.05602052 
    ## Run 192 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001322237  max resid 0.0003136072 
    ## ... Similar to previous best
    ## Run 193 stress 0.05900151 
    ## Run 194 stress 0.05403609 
    ## Run 195 stress 0.05928248 
    ## Run 196 stress 0.05036794 
    ## ... Procrustes: rmse 9.590492e-05  max resid 0.0002286689 
    ## ... Similar to previous best
    ## Run 197 stress 0.06203028 
    ## Run 198 stress 0.06609068 
    ## Run 199 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001239618  max resid 0.0002992933 
    ## ... Similar to previous best
    ## Run 200 stress 0.05051311 
    ## ... Procrustes: rmse 0.03568327  max resid 0.1217935 
    ## Run 201 stress 0.06362787 
    ## Run 202 stress 0.06136683 
    ## Run 203 stress 0.05829931 
    ## Run 204 stress 0.05968628 
    ## Run 205 stress 0.06136682 
    ## Run 206 stress 0.05928254 
    ## Run 207 stress 0.05051293 
    ## ... Procrustes: rmse 0.0357141  max resid 0.1219036 
    ## Run 208 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001810893  max resid 0.0004472517 
    ## ... Similar to previous best
    ## Run 209 stress 0.05036794 
    ## ... Procrustes: rmse 7.378261e-05  max resid 0.0001926208 
    ## ... Similar to previous best
    ## Run 210 stress 0.06478526 
    ## Run 211 stress 0.05036803 
    ## ... Procrustes: rmse 0.000162168  max resid 0.0004381962 
    ## ... Similar to previous best
    ## Run 212 stress 0.0533758 
    ## Run 213 stress 0.05171063 
    ## Run 214 stress 0.05829934 
    ## Run 215 stress 0.05051343 
    ## ... Procrustes: rmse 0.0357259  max resid 0.1220645 
    ## Run 216 stress 0.05337586 
    ## Run 217 stress 0.06111764 
    ## Run 218 stress 0.06231305 
    ## Run 219 stress 0.0560205 
    ## Run 220 stress 0.05036792 
    ## ... Procrustes: rmse 5.812809e-05  max resid 0.0001332369 
    ## ... Similar to previous best
    ## Run 221 stress 0.05403619 
    ## Run 222 stress 0.06309858 
    ## Run 223 stress 0.0596863 
    ## Run 224 stress 0.06231322 
    ## Run 225 stress 0.06038953 
    ## Run 226 stress 0.05036808 
    ## ... Procrustes: rmse 0.00019021  max resid 0.0005189026 
    ## ... Similar to previous best
    ## Run 227 stress 0.05602054 
    ## Run 228 stress 0.05928239 
    ## Run 229 stress 0.05403575 
    ## Run 230 stress 0.05171062 
    ## Run 231 stress 0.05900153 
    ## Run 232 stress 0.06730254 
    ## Run 233 stress 0.06674552 
    ## Run 234 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001851631  max resid 0.0004447059 
    ## ... Similar to previous best
    ## Run 235 stress 0.05829945 
    ## Run 236 stress 0.05466701 
    ## Run 237 stress 0.05171058 
    ## Run 238 stress 0.05171064 
    ## Run 239 stress 0.05171058 
    ## Run 240 stress 0.05829938 
    ## Run 241 stress 0.05051283 
    ## ... Procrustes: rmse 0.03569391  max resid 0.121857 
    ## Run 242 stress 0.05036793 
    ## ... Procrustes: rmse 6.440105e-05  max resid 0.000188003 
    ## ... Similar to previous best
    ## Run 243 stress 0.06231295 
    ## Run 244 stress 0.06478522 
    ## Run 245 stress 0.06478528 
    ## Run 246 stress 0.05968627 
    ## Run 247 stress 0.06136677 
    ## Run 248 stress 0.05171059 
    ## Run 249 stress 0.05602045 
    ## Run 250 stress 0.05466703 
    ## Run 251 stress 0.06027441 
    ## Run 252 stress 0.05051311 
    ## ... Procrustes: rmse 0.03566727  max resid 0.12174 
    ## Run 253 stress 0.06231288 
    ## Run 254 stress 0.05036792 
    ## ... Procrustes: rmse 5.861943e-05  max resid 0.000130796 
    ## ... Similar to previous best
    ## Run 255 stress 0.05602048 
    ## Run 256 stress 0.05829934 
    ## Run 257 stress 0.05829942 
    ## Run 258 stress 0.05051299 
    ## ... Procrustes: rmse 0.03568112  max resid 0.1217949 
    ## Run 259 stress 0.05051298 
    ## ... Procrustes: rmse 0.03564752  max resid 0.1218076 
    ## Run 260 stress 0.05337584 
    ## Run 261 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001421187  max resid 0.0004050222 
    ## ... Similar to previous best
    ## Run 262 stress 0.0623129 
    ## Run 263 stress 0.05979229 
    ## Run 264 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001299449  max resid 0.0003135863 
    ## ... Similar to previous best
    ## Run 265 stress 0.05051358 
    ## ... Procrustes: rmse 0.03567928  max resid 0.1217339 
    ## Run 266 stress 0.05051357 
    ## ... Procrustes: rmse 0.03573547  max resid 0.1220976 
    ## Run 267 stress 0.05171059 
    ## Run 268 stress 0.05051309 
    ## ... Procrustes: rmse 0.03571297  max resid 0.1220118 
    ## Run 269 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001211616  max resid 0.0002904386 
    ## ... Similar to previous best
    ## Run 270 stress 0.05171059 
    ## Run 271 stress 0.05051299 
    ## ... Procrustes: rmse 0.03569498  max resid 0.1219517 
    ## Run 272 stress 0.05602045 
    ## Run 273 stress 0.05968639 
    ## Run 274 stress 0.05968638 
    ## Run 275 stress 0.06136677 
    ## Run 276 stress 0.05466704 
    ## Run 277 stress 0.06537715 
    ## Run 278 stress 0.3545774 
    ## Run 279 stress 0.050368 
    ## ... Procrustes: rmse 0.000156368  max resid 0.0003821368 
    ## ... Similar to previous best
    ## Run 280 stress 0.05051271 
    ## ... Procrustes: rmse 0.03570045  max resid 0.1219119 
    ## Run 281 stress 0.0517106 
    ## Run 282 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001524723  max resid 0.000367824 
    ## ... Similar to previous best
    ## Run 283 stress 0.06111772 
    ## Run 284 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001207772  max resid 0.0002902431 
    ## ... Similar to previous best
    ## Run 285 stress 0.05968644 
    ## Run 286 stress 0.05337579 
    ## Run 287 stress 0.05036792 
    ## ... Procrustes: rmse 3.239412e-05  max resid 8.929067e-05 
    ## ... Similar to previous best
    ## Run 288 stress 0.05602046 
    ## Run 289 stress 0.06837501 
    ## Run 290 stress 0.05602071 
    ## Run 291 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001508139  max resid 0.0003561962 
    ## ... Similar to previous best
    ## Run 292 stress 0.05051291 
    ## ... Procrustes: rmse 0.03568978  max resid 0.1218313 
    ## Run 293 stress 0.06111767 
    ## Run 294 stress 0.05036792 
    ## ... Procrustes: rmse 3.668231e-05  max resid 9.626137e-05 
    ## ... Similar to previous best
    ## Run 295 stress 0.05051336 
    ## ... Procrustes: rmse 0.03570697  max resid 0.1220041 
    ## Run 296 stress 0.06038945 
    ## Run 297 stress 0.05051322 
    ## ... Procrustes: rmse 0.03572736  max resid 0.1220595 
    ## Run 298 stress 0.07242235 
    ## Run 299 stress 0.06136689 
    ## Run 300 stress 0.05337581 
    ## Run 301 stress 0.05602042 
    ## Run 302 stress 0.05829939 
    ## Run 303 stress 0.05337587 
    ## Run 304 stress 0.050368 
    ## ... Procrustes: rmse 0.0001361822  max resid 0.0003930423 
    ## ... Similar to previous best
    ## Run 305 stress 0.05036794 
    ## ... Procrustes: rmse 9.752503e-05  max resid 0.0002332228 
    ## ... Similar to previous best
    ## Run 306 stress 0.06478542 
    ## Run 307 stress 0.05051321 
    ## ... Procrustes: rmse 0.03572869  max resid 0.1220645 
    ## Run 308 stress 0.05466703 
    ## Run 309 stress 0.05051321 
    ## ... Procrustes: rmse 0.03574616  max resid 0.1219681 
    ## Run 310 stress 0.05337576 
    ## Run 311 stress 0.05051313 
    ## ... Procrustes: rmse 0.03571328  max resid 0.1220148 
    ## Run 312 stress 0.05171059 
    ## Run 313 stress 0.06655687 
    ## Run 314 stress 0.05051294 
    ## ... Procrustes: rmse 0.03571169  max resid 0.1218941 
    ## Run 315 stress 0.06111764 
    ## Run 316 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001533239  max resid 0.0003790917 
    ## ... Similar to previous best
    ## Run 317 stress 0.05602052 
    ## Run 318 stress 0.05602044 
    ## Run 319 stress 0.05979232 
    ## Run 320 stress 0.05171058 
    ## Run 321 stress 0.05171059 
    ## Run 322 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001465299  max resid 0.0003937436 
    ## ... Similar to previous best
    ## Run 323 stress 0.05968629 
    ## Run 324 stress 0.0505133 
    ## ... Procrustes: rmse 0.03574537  max resid 0.121959 
    ## Run 325 stress 0.06730204 
    ## Run 326 stress 0.05602064 
    ## Run 327 stress 0.05036792 
    ## ... Procrustes: rmse 2.260344e-05  max resid 4.106642e-05 
    ## ... Similar to previous best
    ## Run 328 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001038216  max resid 0.0002531725 
    ## ... Similar to previous best
    ## Run 329 stress 0.05968637 
    ## Run 330 stress 0.05928222 
    ## Run 331 stress 0.05602048 
    ## Run 332 stress 0.05051284 
    ## ... Procrustes: rmse 0.03569247  max resid 0.1219307 
    ## Run 333 stress 0.05171059 
    ## Run 334 stress 0.06609064 
    ## Run 335 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001802147  max resid 0.0004440069 
    ## ... Similar to previous best
    ## Run 336 stress 0.05466701 
    ## Run 337 stress 0.05337594 
    ## Run 338 stress 0.06655654 
    ## Run 339 stress 0.05036796 
    ## ... Procrustes: rmse 0.000119066  max resid 0.0002773162 
    ## ... Similar to previous best
    ## Run 340 stress 0.06038951 
    ## Run 341 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001362086  max resid 0.0003436688 
    ## ... Similar to previous best
    ## Run 342 stress 0.05602063 
    ## Run 343 stress 0.06537635 
    ## Run 344 stress 0.05337578 
    ## Run 345 stress 0.06938807 
    ## Run 346 stress 0.05051327 
    ## ... Procrustes: rmse 0.03570317  max resid 0.1219908 
    ## Run 347 stress 0.06231287 
    ## Run 348 stress 0.05968634 
    ## Run 349 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001301722  max resid 0.0003688464 
    ## ... Similar to previous best
    ## Run 350 stress 0.05051333 
    ## ... Procrustes: rmse 0.03571757  max resid 0.1218706 
    ## Run 351 stress 0.05171066 
    ## Run 352 stress 0.05403615 
    ## Run 353 stress 0.06203074 
    ## Run 354 stress 0.050368 
    ## ... Procrustes: rmse 0.0001592173  max resid 0.0003862597 
    ## ... Similar to previous best
    ## Run 355 stress 0.06111767 
    ## Run 356 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001038484  max resid 0.0002461471 
    ## ... Similar to previous best
    ## Run 357 stress 0.06231304 
    ## Run 358 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001532593  max resid 0.0004125904 
    ## ... Similar to previous best
    ## Run 359 stress 0.06478536 
    ## Run 360 stress 0.05171058 
    ## Run 361 stress 0.05979226 
    ## Run 362 stress 0.05979236 
    ## Run 363 stress 0.05968627 
    ## Run 364 stress 0.05051347 
    ## ... Procrustes: rmse 0.03570613  max resid 0.1220069 
    ## Run 365 stress 0.05051329 
    ## ... Procrustes: rmse 0.03563648  max resid 0.1217917 
    ## Run 366 stress 0.06362804 
    ## Run 367 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001448072  max resid 0.0003508807 
    ## ... Similar to previous best
    ## Run 368 stress 0.0613668 
    ## Run 369 stress 0.0517106 
    ## Run 370 stress 0.05051296 
    ## ... Procrustes: rmse 0.0357018  max resid 0.1218636 
    ## Run 371 stress 0.06136681 
    ## Run 372 stress 0.06630525 
    ## Run 373 stress 0.05051315 
    ## ... Procrustes: rmse 0.03564003  max resid 0.1217945 
    ## Run 374 stress 0.06609065 
    ## Run 375 stress 0.06621918 
    ## Run 376 stress 0.05466717 
    ## Run 377 stress 0.06478527 
    ## Run 378 stress 0.05051297 
    ## ... Procrustes: rmse 0.03569485  max resid 0.1219497 
    ## Run 379 stress 0.06309851 
    ## Run 380 stress 0.06478535 
    ## Run 381 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001916817  max resid 0.0004899141 
    ## ... Similar to previous best
    ## Run 382 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001470696  max resid 0.0004065974 
    ## ... Similar to previous best
    ## Run 383 stress 0.05602053 
    ## Run 384 stress 0.05602045 
    ## Run 385 stress 0.05968639 
    ## Run 386 stress 0.06309864 
    ## Run 387 stress 0.05051295 
    ## ... Procrustes: rmse 0.03570116  max resid 0.1219673 
    ## Run 388 stress 0.05051292 
    ## ... Procrustes: rmse 0.03583006  max resid 0.1222858 
    ## Run 389 stress 0.05036791 
    ## ... Procrustes: rmse 1.068378e-05  max resid 2.99376e-05 
    ## ... Similar to previous best
    ## Run 390 stress 0.06674561 
    ## Run 391 stress 0.05403567 
    ## Run 392 stress 0.0517106 
    ## Run 393 stress 0.05829928 
    ## Run 394 stress 0.06111767 
    ## Run 395 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001853378  max resid 0.0004454467 
    ## ... Similar to previous best
    ## Run 396 stress 0.05051315 
    ## ... Procrustes: rmse 0.03569928  max resid 0.1219737 
    ## Run 397 stress 0.0505136 
    ## ... Procrustes: rmse 0.03574133  max resid 0.1221141 
    ## Run 398 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001502928  max resid 0.0003705433 
    ## ... Similar to previous best
    ## Run 399 stress 0.0505134 
    ## ... Procrustes: rmse 0.03576193  max resid 0.1221702 
    ## Run 400 stress 0.06504214 
    ## Run 401 stress 0.06362793 
    ## Run 402 stress 0.05051306 
    ## ... Procrustes: rmse 0.03569206  max resid 0.1218188 
    ## Run 403 stress 0.0597923 
    ## Run 404 stress 0.05928245 
    ## Run 405 stress 0.06730239 
    ## Run 406 stress 0.06038956 
    ## Run 407 stress 0.0663053 
    ## Run 408 stress 0.06362779 
    ## Run 409 stress 0.05968625 
    ## Run 410 stress 0.050368 
    ## ... Procrustes: rmse 0.0001348604  max resid 0.0003584918 
    ## ... Similar to previous best
    ## Run 411 stress 0.05403614 
    ## Run 412 stress 0.05051326 
    ## ... Procrustes: rmse 0.03574483  max resid 0.1221119 
    ## Run 413 stress 0.05051311 
    ## ... Procrustes: rmse 0.03569798  max resid 0.1218318 
    ## Run 414 stress 0.06309857 
    ## Run 415 stress 0.05829936 
    ## Run 416 stress 0.05051333 
    ## ... Procrustes: rmse 0.0356647  max resid 0.1217102 
    ## Run 417 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001091203  max resid 0.0002973001 
    ## ... Similar to previous best
    ## Run 418 stress 0.0517106 
    ## Run 419 stress 0.05466709 
    ## Run 420 stress 0.05829933 
    ## Run 421 stress 0.05036795 
    ## ... Procrustes: rmse 8.972938e-05  max resid 0.0002371181 
    ## ... Similar to previous best
    ## Run 422 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001162937  max resid 0.000335924 
    ## ... Similar to previous best
    ## Run 423 stress 0.05466702 
    ## Run 424 stress 0.05036791 
    ## ... Procrustes: rmse 2.854474e-05  max resid 6.069519e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.05403592 
    ## Run 426 stress 0.05928212 
    ## Run 427 stress 0.06845252 
    ## Run 428 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001482471  max resid 0.000362816 
    ## ... Similar to previous best
    ## Run 429 stress 0.05051332 
    ## ... Procrustes: rmse 0.03571882  max resid 0.1220369 
    ## Run 430 stress 0.06362805 
    ## Run 431 stress 0.06203055 
    ## Run 432 stress 0.0590016 
    ## Run 433 stress 0.05051318 
    ## ... Procrustes: rmse 0.03566418  max resid 0.121726 
    ## Run 434 stress 0.05968634 
    ## Run 435 stress 0.05036809 
    ## ... Procrustes: rmse 0.0002086466  max resid 0.0005148435 
    ## ... Similar to previous best
    ## Run 436 stress 0.05171058 
    ## Run 437 stress 0.06730239 
    ## Run 438 stress 0.06885107 
    ## Run 439 stress 0.0582993 
    ## Run 440 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001897637  max resid 0.0004660737 
    ## ... Similar to previous best
    ## Run 441 stress 0.05051305 
    ## ... Procrustes: rmse 0.03570518  max resid 0.1219861 
    ## Run 442 stress 0.05051303 
    ## ... Procrustes: rmse 0.03569848  max resid 0.1218424 
    ## Run 443 stress 0.05829931 
    ## Run 444 stress 0.05051359 
    ## ... Procrustes: rmse 0.03574343  max resid 0.1221181 
    ## Run 445 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001871969  max resid 0.0004459169 
    ## ... Similar to previous best
    ## Run 446 stress 0.0517106 
    ## Run 447 stress 0.06136679 
    ## Run 448 stress 0.05051296 
    ## ... Procrustes: rmse 0.03569355  max resid 0.1219439 
    ## Run 449 stress 0.05900155 
    ## Run 450 stress 0.0533758 
    ## Run 451 stress 0.05968638 
    ## Run 452 stress 0.05602055 
    ## Run 453 stress 0.05036802 
    ## ... Procrustes: rmse 0.000172995  max resid 0.0004097808 
    ## ... Similar to previous best
    ## Run 454 stress 0.05051351 
    ## ... Procrustes: rmse 0.03567897  max resid 0.1217406 
    ## Run 455 stress 0.05968638 
    ## Run 456 stress 0.06362811 
    ## Run 457 stress 0.05337577 
    ## Run 458 stress 0.05979235 
    ## Run 459 stress 0.06136691 
    ## Run 460 stress 0.0582993 
    ## Run 461 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001457303  max resid 0.0003559689 
    ## ... Similar to previous best
    ## Run 462 stress 0.05602071 
    ## Run 463 stress 0.06844318 
    ## Run 464 stress 0.05036795 
    ## ... Procrustes: rmse 9.676426e-05  max resid 0.0002296594 
    ## ... Similar to previous best
    ## Run 465 stress 0.05602053 
    ## Run 466 stress 0.05036792 
    ## ... Procrustes: rmse 4.772319e-05  max resid 0.0001089309 
    ## ... Similar to previous best
    ## Run 467 stress 0.06743844 
    ## Run 468 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001202262  max resid 0.0002901585 
    ## ... Similar to previous best
    ## Run 469 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001307553  max resid 0.0003127274 
    ## ... Similar to previous best
    ## Run 470 stress 0.05051285 
    ## ... Procrustes: rmse 0.03570671  max resid 0.1218927 
    ## Run 471 stress 0.0630986 
    ## Run 472 stress 0.05171062 
    ## Run 473 stress 0.06630526 
    ## Run 474 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001201798  max resid 0.0002885775 
    ## ... Similar to previous best
    ## Run 475 stress 0.05171059 
    ## Run 476 stress 0.05171058 
    ## Run 477 stress 0.05466709 
    ## Run 478 stress 0.05036792 
    ## ... Procrustes: rmse 5.670042e-05  max resid 0.0001230838 
    ## ... Similar to previous best
    ## Run 479 stress 0.05036797 
    ## ... Procrustes: rmse 9.37072e-05  max resid 0.0002429556 
    ## ... Similar to previous best
    ## Run 480 stress 0.06640339 
    ## Run 481 stress 0.05829933 
    ## Run 482 stress 0.05337582 
    ## Run 483 stress 0.06309845 
    ## Run 484 stress 0.06609066 
    ## Run 485 stress 0.05171058 
    ## Run 486 stress 0.05036792 
    ## ... Procrustes: rmse 2.599975e-05  max resid 7.139208e-05 
    ## ... Similar to previous best
    ## Run 487 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001071293  max resid 0.0002593547 
    ## ... Similar to previous best
    ## Run 488 stress 0.0505134 
    ## ... Procrustes: rmse 0.03575086  max resid 0.1219699 
    ## Run 489 stress 0.05036793 
    ## ... Procrustes: rmse 7.428193e-05  max resid 0.0001765837 
    ## ... Similar to previous best
    ## Run 490 stress 0.05466704 
    ## Run 491 stress 0.05466703 
    ## Run 492 stress 0.05036795 
    ## ... Procrustes: rmse 7.002429e-05  max resid 0.0001966852 
    ## ... Similar to previous best
    ## Run 493 stress 0.05036794 
    ## ... Procrustes: rmse 7.915559e-05  max resid 0.0002105113 
    ## ... Similar to previous best
    ## Run 494 stress 0.05337578 
    ## Run 495 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001894904  max resid 0.0004685112 
    ## ... Similar to previous best
    ## Run 496 stress 0.05051323 
    ## ... Procrustes: rmse 0.03570497  max resid 0.1219947 
    ## Run 497 stress 0.05829946 
    ## Run 498 stress 0.05051302 
    ## ... Procrustes: rmse 0.03570387  max resid 0.1219789 
    ## Run 499 stress 0.0603895 
    ## Run 500 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001509179  max resid 0.0003993863 
    ## ... Similar to previous best
    ## *** Best solution repeated 76 times

``` r
# Ocean sites and mixed lakes
PD_beta_env_OM_NMDS <- metaMDS(PD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08352634 
    ## Run 1 stress 0.1134014 
    ## Run 2 stress 0.08352634 
    ## ... Procrustes: rmse 7.89325e-06  max resid 1.761608e-05 
    ## ... Similar to previous best
    ## Run 3 stress 0.08352634 
    ## ... Procrustes: rmse 4.056456e-06  max resid 8.472505e-06 
    ## ... Similar to previous best
    ## Run 4 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 4.836882e-07  max resid 8.893234e-07 
    ## ... Similar to previous best
    ## Run 5 stress 0.08352634 
    ## ... Procrustes: rmse 6.825662e-06  max resid 1.519864e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.1134015 
    ## Run 7 stress 0.1134014 
    ## Run 8 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.759966e-06  max resid 3.940434e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.08352634 
    ## ... Procrustes: rmse 6.722209e-07  max resid 1.418431e-06 
    ## ... Similar to previous best
    ## Run 10 stress 0.08352634 
    ## ... Procrustes: rmse 1.00068e-06  max resid 2.002053e-06 
    ## ... Similar to previous best
    ## Run 11 stress 0.08352634 
    ## ... Procrustes: rmse 2.70243e-06  max resid 5.851172e-06 
    ## ... Similar to previous best
    ## Run 12 stress 0.08352634 
    ## ... Procrustes: rmse 3.720454e-06  max resid 8.452876e-06 
    ## ... Similar to previous best
    ## Run 13 stress 0.1134014 
    ## Run 14 stress 0.205005 
    ## Run 15 stress 0.08352634 
    ## ... Procrustes: rmse 4.008697e-06  max resid 9.006015e-06 
    ## ... Similar to previous best
    ## Run 16 stress 0.1279309 
    ## Run 17 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.494235e-06  max resid 2.968975e-06 
    ## ... Similar to previous best
    ## Run 18 stress 0.08352634 
    ## ... Procrustes: rmse 3.148229e-06  max resid 6.858554e-06 
    ## ... Similar to previous best
    ## Run 19 stress 0.1134014 
    ## Run 20 stress 0.08352634 
    ## ... Procrustes: rmse 3.048338e-06  max resid 6.358481e-06 
    ## ... Similar to previous best
    ## Run 21 stress 0.1279309 
    ## Run 22 stress 0.1134015 
    ## Run 23 stress 0.08352634 
    ## ... Procrustes: rmse 1.987535e-06  max resid 4.064243e-06 
    ## ... Similar to previous best
    ## Run 24 stress 0.08352634 
    ## ... Procrustes: rmse 2.74291e-06  max resid 5.972752e-06 
    ## ... Similar to previous best
    ## Run 25 stress 0.1279309 
    ## Run 26 stress 0.1134014 
    ## Run 27 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.339553e-06  max resid 2.395029e-06 
    ## ... Similar to previous best
    ## Run 28 stress 0.08352634 
    ## ... Procrustes: rmse 3.966702e-06  max resid 6.021563e-06 
    ## ... Similar to previous best
    ## Run 29 stress 0.1279309 
    ## Run 30 stress 0.08352634 
    ## ... Procrustes: rmse 3.87637e-06  max resid 8.107853e-06 
    ## ... Similar to previous best
    ## Run 31 stress 0.127931 
    ## Run 32 stress 0.08352634 
    ## ... Procrustes: rmse 1.766915e-06  max resid 3.28084e-06 
    ## ... Similar to previous best
    ## Run 33 stress 0.1279309 
    ## Run 34 stress 0.08352634 
    ## ... Procrustes: rmse 3.730727e-06  max resid 8.097995e-06 
    ## ... Similar to previous best
    ## Run 35 stress 0.2031893 
    ## Run 36 stress 0.1134014 
    ## Run 37 stress 0.08352634 
    ## ... Procrustes: rmse 2.171228e-06  max resid 4.927112e-06 
    ## ... Similar to previous best
    ## Run 38 stress 0.08352634 
    ## ... Procrustes: rmse 2.306925e-06  max resid 4.988219e-06 
    ## ... Similar to previous best
    ## Run 39 stress 0.1989485 
    ## Run 40 stress 0.08352634 
    ## ... Procrustes: rmse 2.05446e-06  max resid 4.222141e-06 
    ## ... Similar to previous best
    ## Run 41 stress 0.1134014 
    ## Run 42 stress 0.08352634 
    ## ... Procrustes: rmse 9.101803e-07  max resid 1.909852e-06 
    ## ... Similar to previous best
    ## Run 43 stress 0.08352634 
    ## ... Procrustes: rmse 3.32881e-06  max resid 7.315749e-06 
    ## ... Similar to previous best
    ## Run 44 stress 0.08352634 
    ## ... Procrustes: rmse 7.137018e-07  max resid 1.41335e-06 
    ## ... Similar to previous best
    ## Run 45 stress 0.08352634 
    ## ... Procrustes: rmse 1.629649e-06  max resid 3.091766e-06 
    ## ... Similar to previous best
    ## Run 46 stress 0.1134014 
    ## Run 47 stress 0.1134015 
    ## Run 48 stress 0.08352634 
    ## ... Procrustes: rmse 2.267453e-06  max resid 5.089765e-06 
    ## ... Similar to previous best
    ## Run 49 stress 0.1279309 
    ## Run 50 stress 0.1279309 
    ## Run 51 stress 0.1279309 
    ## Run 52 stress 0.08352634 
    ## ... Procrustes: rmse 3.491997e-06  max resid 7.636338e-06 
    ## ... Similar to previous best
    ## Run 53 stress 0.1279309 
    ## Run 54 stress 0.1134014 
    ## Run 55 stress 0.08352634 
    ## ... Procrustes: rmse 2.033464e-06  max resid 4.633796e-06 
    ## ... Similar to previous best
    ## Run 56 stress 0.08352634 
    ## ... Procrustes: rmse 1.649275e-06  max resid 3.078343e-06 
    ## ... Similar to previous best
    ## Run 57 stress 0.08352634 
    ## ... Procrustes: rmse 1.547218e-06  max resid 3.396065e-06 
    ## ... Similar to previous best
    ## Run 58 stress 0.1989155 
    ## Run 59 stress 0.1279309 
    ## Run 60 stress 0.08352634 
    ## ... Procrustes: rmse 3.33753e-06  max resid 7.443061e-06 
    ## ... Similar to previous best
    ## Run 61 stress 0.08352634 
    ## ... Procrustes: rmse 1.973618e-06  max resid 4.457846e-06 
    ## ... Similar to previous best
    ## Run 62 stress 0.1134014 
    ## Run 63 stress 0.08352634 
    ## ... Procrustes: rmse 1.414861e-06  max resid 2.599013e-06 
    ## ... Similar to previous best
    ## Run 64 stress 0.08352634 
    ## ... Procrustes: rmse 2.933451e-06  max resid 6.637028e-06 
    ## ... Similar to previous best
    ## Run 65 stress 0.1989157 
    ## Run 66 stress 0.1134014 
    ## Run 67 stress 0.1134015 
    ## Run 68 stress 0.1134014 
    ## Run 69 stress 0.08352634 
    ## ... Procrustes: rmse 2.760242e-06  max resid 6.239813e-06 
    ## ... Similar to previous best
    ## Run 70 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.839527e-06  max resid 3.963131e-06 
    ## ... Similar to previous best
    ## Run 71 stress 0.08352634 
    ## ... Procrustes: rmse 2.677546e-06  max resid 5.882765e-06 
    ## ... Similar to previous best
    ## Run 72 stress 0.08352634 
    ## ... Procrustes: rmse 4.002804e-06  max resid 8.839376e-06 
    ## ... Similar to previous best
    ## Run 73 stress 0.08352634 
    ## ... Procrustes: rmse 2.671376e-06  max resid 5.998724e-06 
    ## ... Similar to previous best
    ## Run 74 stress 0.1134014 
    ## Run 75 stress 0.1134014 
    ## Run 76 stress 0.08352634 
    ## ... Procrustes: rmse 1.136971e-06  max resid 1.739419e-06 
    ## ... Similar to previous best
    ## Run 77 stress 0.1134014 
    ## Run 78 stress 0.3154697 
    ## Run 79 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 1.340186e-06  max resid 2.724722e-06 
    ## ... Similar to previous best
    ## Run 80 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 5.355945e-07  max resid 1.053937e-06 
    ## ... Similar to previous best
    ## Run 81 stress 0.08352634 
    ## ... Procrustes: rmse 3.149716e-06  max resid 6.894619e-06 
    ## ... Similar to previous best
    ## Run 82 stress 0.1134014 
    ## Run 83 stress 0.1279309 
    ## Run 84 stress 0.1279309 
    ## Run 85 stress 0.08352634 
    ## ... Procrustes: rmse 2.39429e-06  max resid 5.481571e-06 
    ## ... Similar to previous best
    ## Run 86 stress 0.1134015 
    ## Run 87 stress 0.08352634 
    ## ... Procrustes: rmse 2.359757e-06  max resid 5.002201e-06 
    ## ... Similar to previous best
    ## Run 88 stress 0.1279309 
    ## Run 89 stress 0.08352634 
    ## ... Procrustes: rmse 5.883419e-06  max resid 1.313117e-05 
    ## ... Similar to previous best
    ## Run 90 stress 0.08352634 
    ## ... Procrustes: rmse 4.781899e-07  max resid 7.809598e-07 
    ## ... Similar to previous best
    ## Run 91 stress 0.08352634 
    ## ... Procrustes: rmse 3.009887e-07  max resid 4.31607e-07 
    ## ... Similar to previous best
    ## Run 92 stress 0.1134014 
    ## Run 93 stress 0.08352634 
    ## ... Procrustes: rmse 1.777803e-06  max resid 3.202743e-06 
    ## ... Similar to previous best
    ## Run 94 stress 0.08352634 
    ## ... Procrustes: rmse 2.049083e-06  max resid 4.725316e-06 
    ## ... Similar to previous best
    ## Run 95 stress 0.1134015 
    ## Run 96 stress 0.1134014 
    ## Run 97 stress 0.08352634 
    ## ... Procrustes: rmse 3.874245e-06  max resid 8.463441e-06 
    ## ... Similar to previous best
    ## Run 98 stress 0.1279309 
    ## Run 99 stress 0.1134014 
    ## Run 100 stress 0.08352634 
    ## ... Procrustes: rmse 2.426241e-06  max resid 5.245678e-06 
    ## ... Similar to previous best
    ## Run 101 stress 0.08352634 
    ## ... Procrustes: rmse 1.477273e-06  max resid 3.023897e-06 
    ## ... Similar to previous best
    ## Run 102 stress 0.08352634 
    ## ... Procrustes: rmse 7.346906e-07  max resid 1.534384e-06 
    ## ... Similar to previous best
    ## Run 103 stress 0.08352634 
    ## ... Procrustes: rmse 3.005639e-06  max resid 6.572393e-06 
    ## ... Similar to previous best
    ## Run 104 stress 0.1279309 
    ## Run 105 stress 0.08352634 
    ## ... Procrustes: rmse 9.301433e-07  max resid 1.850952e-06 
    ## ... Similar to previous best
    ## Run 106 stress 0.1134015 
    ## Run 107 stress 0.08352634 
    ## ... Procrustes: rmse 1.829387e-06  max resid 3.91701e-06 
    ## ... Similar to previous best
    ## Run 108 stress 0.08352634 
    ## ... Procrustes: rmse 2.969205e-06  max resid 6.72852e-06 
    ## ... Similar to previous best
    ## Run 109 stress 0.08352634 
    ## ... Procrustes: rmse 1.738189e-06  max resid 3.431963e-06 
    ## ... Similar to previous best
    ## Run 110 stress 0.08352634 
    ## ... Procrustes: rmse 2.610485e-06  max resid 5.906714e-06 
    ## ... Similar to previous best
    ## Run 111 stress 0.1134015 
    ## Run 112 stress 0.1134014 
    ## Run 113 stress 0.08352634 
    ## ... Procrustes: rmse 1.226018e-06  max resid 2.631628e-06 
    ## ... Similar to previous best
    ## Run 114 stress 0.1279309 
    ## Run 115 stress 0.08352634 
    ## ... Procrustes: rmse 1.985187e-06  max resid 4.251779e-06 
    ## ... Similar to previous best
    ## Run 116 stress 0.08352634 
    ## ... Procrustes: rmse 3.239387e-06  max resid 7.025321e-06 
    ## ... Similar to previous best
    ## Run 117 stress 0.1279309 
    ## Run 118 stress 0.1134015 
    ## Run 119 stress 0.08352634 
    ## ... Procrustes: rmse 2.395875e-06  max resid 5.223428e-06 
    ## ... Similar to previous best
    ## Run 120 stress 0.2034292 
    ## Run 121 stress 0.1989154 
    ## Run 122 stress 0.08352634 
    ## ... Procrustes: rmse 1.946216e-06  max resid 4.361439e-06 
    ## ... Similar to previous best
    ## Run 123 stress 0.08352634 
    ## ... Procrustes: rmse 1.629023e-06  max resid 3.213978e-06 
    ## ... Similar to previous best
    ## Run 124 stress 0.1279309 
    ## Run 125 stress 0.08352634 
    ## ... Procrustes: rmse 3.204389e-06  max resid 6.626473e-06 
    ## ... Similar to previous best
    ## Run 126 stress 0.1279309 
    ## Run 127 stress 0.08352634 
    ## ... Procrustes: rmse 3.359411e-06  max resid 7.561135e-06 
    ## ... Similar to previous best
    ## Run 128 stress 0.1279309 
    ## Run 129 stress 0.08352634 
    ## ... Procrustes: rmse 1.619977e-06  max resid 3.71565e-06 
    ## ... Similar to previous best
    ## Run 130 stress 0.1134015 
    ## Run 131 stress 0.08352634 
    ## ... Procrustes: rmse 2.403446e-06  max resid 5.112506e-06 
    ## ... Similar to previous best
    ## Run 132 stress 0.08352634 
    ## ... Procrustes: rmse 3.707344e-06  max resid 8.089734e-06 
    ## ... Similar to previous best
    ## Run 133 stress 0.08352634 
    ## ... Procrustes: rmse 6.542642e-07  max resid 9.929911e-07 
    ## ... Similar to previous best
    ## Run 134 stress 0.08352634 
    ## ... Procrustes: rmse 1.613384e-06  max resid 3.241327e-06 
    ## ... Similar to previous best
    ## Run 135 stress 0.08352634 
    ## ... Procrustes: rmse 3.400755e-06  max resid 7.717996e-06 
    ## ... Similar to previous best
    ## Run 136 stress 0.1279309 
    ## Run 137 stress 0.08352634 
    ## ... Procrustes: rmse 1.937516e-06  max resid 4.464752e-06 
    ## ... Similar to previous best
    ## Run 138 stress 0.08352634 
    ## ... Procrustes: rmse 1.476029e-06  max resid 3.391085e-06 
    ## ... Similar to previous best
    ## Run 139 stress 0.1134015 
    ## Run 140 stress 0.08352634 
    ## ... Procrustes: rmse 2.029275e-06  max resid 4.308444e-06 
    ## ... Similar to previous best
    ## Run 141 stress 0.08352634 
    ## ... Procrustes: rmse 3.29836e-06  max resid 7.163362e-06 
    ## ... Similar to previous best
    ## Run 142 stress 0.08352634 
    ## ... Procrustes: rmse 5.715283e-06  max resid 1.253645e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.08352634 
    ## ... Procrustes: rmse 2.012876e-06  max resid 3.249803e-06 
    ## ... Similar to previous best
    ## Run 144 stress 0.1134014 
    ## Run 145 stress 0.1279309 
    ## Run 146 stress 0.08352634 
    ## ... Procrustes: rmse 2.634822e-06  max resid 5.761133e-06 
    ## ... Similar to previous best
    ## Run 147 stress 0.1279309 
    ## Run 148 stress 0.2031894 
    ## Run 149 stress 0.08352634 
    ## ... Procrustes: rmse 3.373238e-06  max resid 7.583399e-06 
    ## ... Similar to previous best
    ## Run 150 stress 0.1279309 
    ## Run 151 stress 0.08352634 
    ## ... Procrustes: rmse 1.423325e-06  max resid 2.171719e-06 
    ## ... Similar to previous best
    ## Run 152 stress 0.1134014 
    ## Run 153 stress 0.127931 
    ## Run 154 stress 0.08352634 
    ## ... Procrustes: rmse 6.709134e-07  max resid 1.193447e-06 
    ## ... Similar to previous best
    ## Run 155 stress 0.1134014 
    ## Run 156 stress 0.1134015 
    ## Run 157 stress 0.08352634 
    ## ... Procrustes: rmse 1.334995e-06  max resid 2.442041e-06 
    ## ... Similar to previous best
    ## Run 158 stress 0.1134014 
    ## Run 159 stress 0.08352634 
    ## ... Procrustes: rmse 8.605457e-07  max resid 1.828113e-06 
    ## ... Similar to previous best
    ## Run 160 stress 0.1134015 
    ## Run 161 stress 0.1134015 
    ## Run 162 stress 0.1134015 
    ## Run 163 stress 0.1134014 
    ## Run 164 stress 0.1279309 
    ## Run 165 stress 0.08352634 
    ## ... Procrustes: rmse 9.057339e-07  max resid 1.674324e-06 
    ## ... Similar to previous best
    ## Run 166 stress 0.1134015 
    ## Run 167 stress 0.08352634 
    ## ... Procrustes: rmse 1.112173e-06  max resid 2.357482e-06 
    ## ... Similar to previous best
    ## Run 168 stress 0.08352634 
    ## ... Procrustes: rmse 5.250701e-07  max resid 1.052958e-06 
    ## ... Similar to previous best
    ## Run 169 stress 0.1134015 
    ## Run 170 stress 0.1279309 
    ## Run 171 stress 0.1134014 
    ## Run 172 stress 0.1134015 
    ## Run 173 stress 0.08352634 
    ## ... Procrustes: rmse 4.175096e-06  max resid 9.152657e-06 
    ## ... Similar to previous best
    ## Run 174 stress 0.1989156 
    ## Run 175 stress 0.1134014 
    ## Run 176 stress 0.08352634 
    ## ... Procrustes: rmse 1.135861e-06  max resid 1.900391e-06 
    ## ... Similar to previous best
    ## Run 177 stress 0.08352634 
    ## ... Procrustes: rmse 2.237575e-06  max resid 4.933766e-06 
    ## ... Similar to previous best
    ## Run 178 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 3.776149e-07  max resid 5.775079e-07 
    ## ... Similar to previous best
    ## Run 179 stress 0.1279309 
    ## Run 180 stress 0.1134014 
    ## Run 181 stress 0.08352634 
    ## ... Procrustes: rmse 5.662441e-07  max resid 1.357347e-06 
    ## ... Similar to previous best
    ## Run 182 stress 0.08352634 
    ## ... Procrustes: rmse 2.56339e-06  max resid 5.915787e-06 
    ## ... Similar to previous best
    ## Run 183 stress 0.3116494 
    ## Run 184 stress 0.1134015 
    ## Run 185 stress 0.08352634 
    ## ... Procrustes: rmse 6.792856e-07  max resid 1.195193e-06 
    ## ... Similar to previous best
    ## Run 186 stress 0.1279309 
    ## Run 187 stress 0.1989156 
    ## Run 188 stress 0.08352634 
    ## ... Procrustes: rmse 1.706977e-06  max resid 2.845298e-06 
    ## ... Similar to previous best
    ## Run 189 stress 0.08352634 
    ## ... Procrustes: rmse 1.515503e-06  max resid 3.537901e-06 
    ## ... Similar to previous best
    ## Run 190 stress 0.08352634 
    ## ... Procrustes: rmse 2.635951e-06  max resid 5.619932e-06 
    ## ... Similar to previous best
    ## Run 191 stress 0.1989157 
    ## Run 192 stress 0.08352634 
    ## ... Procrustes: rmse 1.126949e-06  max resid 2.79331e-06 
    ## ... Similar to previous best
    ## Run 193 stress 0.1989486 
    ## Run 194 stress 0.1134014 
    ## Run 195 stress 0.1134014 
    ## Run 196 stress 0.08352634 
    ## ... Procrustes: rmse 1.898922e-06  max resid 3.96592e-06 
    ## ... Similar to previous best
    ## Run 197 stress 0.08352634 
    ## ... Procrustes: rmse 1.322599e-06  max resid 3.00974e-06 
    ## ... Similar to previous best
    ## Run 198 stress 0.1134014 
    ## Run 199 stress 0.08352634 
    ## ... Procrustes: rmse 1.568803e-06  max resid 3.524411e-06 
    ## ... Similar to previous best
    ## Run 200 stress 0.08352634 
    ## ... Procrustes: rmse 2.018304e-06  max resid 4.751041e-06 
    ## ... Similar to previous best
    ## Run 201 stress 0.08352634 
    ## ... Procrustes: rmse 2.64917e-06  max resid 6.134534e-06 
    ## ... Similar to previous best
    ## Run 202 stress 0.08352634 
    ## ... Procrustes: rmse 5.14252e-07  max resid 1.272729e-06 
    ## ... Similar to previous best
    ## Run 203 stress 0.08352634 
    ## ... Procrustes: rmse 7.728649e-07  max resid 1.163384e-06 
    ## ... Similar to previous best
    ## Run 204 stress 0.08352634 
    ## ... Procrustes: rmse 2.859474e-06  max resid 6.378334e-06 
    ## ... Similar to previous best
    ## Run 205 stress 0.08352634 
    ## ... Procrustes: rmse 6.292322e-06  max resid 1.450149e-05 
    ## ... Similar to previous best
    ## Run 206 stress 0.08352634 
    ## ... Procrustes: rmse 2.905105e-06  max resid 5.990122e-06 
    ## ... Similar to previous best
    ## Run 207 stress 0.1134014 
    ## Run 208 stress 0.08352634 
    ## ... Procrustes: rmse 3.827632e-06  max resid 8.190675e-06 
    ## ... Similar to previous best
    ## Run 209 stress 0.1134014 
    ## Run 210 stress 0.08352634 
    ## ... Procrustes: rmse 1.907523e-06  max resid 3.87379e-06 
    ## ... Similar to previous best
    ## Run 211 stress 0.1279309 
    ## Run 212 stress 0.08352634 
    ## ... Procrustes: rmse 1.86135e-06  max resid 4.431857e-06 
    ## ... Similar to previous best
    ## Run 213 stress 0.08352634 
    ## ... Procrustes: rmse 3.398965e-06  max resid 7.781438e-06 
    ## ... Similar to previous best
    ## Run 214 stress 0.331149 
    ## Run 215 stress 0.1279309 
    ## Run 216 stress 0.08352634 
    ## ... Procrustes: rmse 1.005108e-06  max resid 2.211842e-06 
    ## ... Similar to previous best
    ## Run 217 stress 0.1134014 
    ## Run 218 stress 0.1279309 
    ## Run 219 stress 0.08352634 
    ## ... Procrustes: rmse 1.455031e-06  max resid 3.449287e-06 
    ## ... Similar to previous best
    ## Run 220 stress 0.08352634 
    ## ... Procrustes: rmse 3.04423e-07  max resid 4.647275e-07 
    ## ... Similar to previous best
    ## Run 221 stress 0.08352634 
    ## ... Procrustes: rmse 3.477516e-06  max resid 7.977424e-06 
    ## ... Similar to previous best
    ## Run 222 stress 0.1134014 
    ## Run 223 stress 0.08352634 
    ## ... Procrustes: rmse 1.213815e-06  max resid 2.160329e-06 
    ## ... Similar to previous best
    ## Run 224 stress 0.08352634 
    ## ... Procrustes: rmse 2.790853e-06  max resid 6.46851e-06 
    ## ... Similar to previous best
    ## Run 225 stress 0.205005 
    ## Run 226 stress 0.08352634 
    ## ... Procrustes: rmse 1.641048e-06  max resid 3.846458e-06 
    ## ... Similar to previous best
    ## Run 227 stress 0.08352634 
    ## ... Procrustes: rmse 9.46621e-07  max resid 1.94729e-06 
    ## ... Similar to previous best
    ## Run 228 stress 0.1134014 
    ## Run 229 stress 0.08352634 
    ## ... Procrustes: rmse 1.717176e-06  max resid 4.071648e-06 
    ## ... Similar to previous best
    ## Run 230 stress 0.08352634 
    ## ... Procrustes: rmse 1.08836e-06  max resid 1.864694e-06 
    ## ... Similar to previous best
    ## Run 231 stress 0.1279309 
    ## Run 232 stress 0.08352634 
    ## ... New best solution
    ## ... Procrustes: rmse 5.526229e-07  max resid 7.70602e-07 
    ## ... Similar to previous best
    ## Run 233 stress 0.08352634 
    ## ... Procrustes: rmse 7.827527e-07  max resid 1.521507e-06 
    ## ... Similar to previous best
    ## Run 234 stress 0.127931 
    ## Run 235 stress 0.08352634 
    ## ... Procrustes: rmse 3.364043e-06  max resid 7.59811e-06 
    ## ... Similar to previous best
    ## Run 236 stress 0.1989159 
    ## Run 237 stress 0.08352634 
    ## ... Procrustes: rmse 5.731815e-06  max resid 1.252137e-05 
    ## ... Similar to previous best
    ## Run 238 stress 0.08352634 
    ## ... Procrustes: rmse 3.107549e-06  max resid 6.951016e-06 
    ## ... Similar to previous best
    ## Run 239 stress 0.1134014 
    ## Run 240 stress 0.08352634 
    ## ... Procrustes: rmse 3.758013e-06  max resid 8.192708e-06 
    ## ... Similar to previous best
    ## Run 241 stress 0.1279309 
    ## Run 242 stress 0.1134014 
    ## Run 243 stress 0.1134014 
    ## Run 244 stress 0.08352634 
    ## ... Procrustes: rmse 3.752764e-06  max resid 8.299255e-06 
    ## ... Similar to previous best
    ## Run 245 stress 0.08352634 
    ## ... Procrustes: rmse 4.34881e-07  max resid 8.781831e-07 
    ## ... Similar to previous best
    ## Run 246 stress 0.1134014 
    ## Run 247 stress 0.08352634 
    ## ... Procrustes: rmse 1.31628e-06  max resid 2.617413e-06 
    ## ... Similar to previous best
    ## Run 248 stress 0.08352634 
    ## ... Procrustes: rmse 3.471651e-06  max resid 7.637048e-06 
    ## ... Similar to previous best
    ## Run 249 stress 0.1279309 
    ## Run 250 stress 0.1134015 
    ## Run 251 stress 0.08352634 
    ## ... Procrustes: rmse 1.390498e-06  max resid 2.844647e-06 
    ## ... Similar to previous best
    ## Run 252 stress 0.08352634 
    ## ... Procrustes: rmse 3.814716e-06  max resid 8.520752e-06 
    ## ... Similar to previous best
    ## Run 253 stress 0.08352634 
    ## ... Procrustes: rmse 1.976919e-06  max resid 4.131062e-06 
    ## ... Similar to previous best
    ## Run 254 stress 0.1134014 
    ## Run 255 stress 0.1134015 
    ## Run 256 stress 0.08352634 
    ## ... Procrustes: rmse 4.512145e-06  max resid 9.835622e-06 
    ## ... Similar to previous best
    ## Run 257 stress 0.08352634 
    ## ... Procrustes: rmse 1.375207e-06  max resid 2.428491e-06 
    ## ... Similar to previous best
    ## Run 258 stress 0.08352634 
    ## ... Procrustes: rmse 2.481939e-06  max resid 5.448463e-06 
    ## ... Similar to previous best
    ## Run 259 stress 0.08352634 
    ## ... Procrustes: rmse 7.019718e-06  max resid 1.588909e-05 
    ## ... Similar to previous best
    ## Run 260 stress 0.2938009 
    ## Run 261 stress 0.08352634 
    ## ... Procrustes: rmse 4.024818e-07  max resid 7.677695e-07 
    ## ... Similar to previous best
    ## Run 262 stress 0.08352634 
    ## ... Procrustes: rmse 3.099324e-06  max resid 5.55083e-06 
    ## ... Similar to previous best
    ## Run 263 stress 0.1134014 
    ## Run 264 stress 0.08352634 
    ## ... Procrustes: rmse 2.269757e-06  max resid 4.782373e-06 
    ## ... Similar to previous best
    ## Run 265 stress 0.08352634 
    ## ... Procrustes: rmse 5.6317e-07  max resid 1.176862e-06 
    ## ... Similar to previous best
    ## Run 266 stress 0.08352634 
    ## ... Procrustes: rmse 2.259498e-06  max resid 5.008732e-06 
    ## ... Similar to previous best
    ## Run 267 stress 0.08352634 
    ## ... Procrustes: rmse 3.313725e-06  max resid 7.544393e-06 
    ## ... Similar to previous best
    ## Run 268 stress 0.08352634 
    ## ... Procrustes: rmse 6.471447e-07  max resid 1.398733e-06 
    ## ... Similar to previous best
    ## Run 269 stress 0.08352634 
    ## ... Procrustes: rmse 4.051664e-06  max resid 9.096832e-06 
    ## ... Similar to previous best
    ## Run 270 stress 0.1134014 
    ## Run 271 stress 0.1279309 
    ## Run 272 stress 0.1134015 
    ## Run 273 stress 0.1134015 
    ## Run 274 stress 0.1279309 
    ## Run 275 stress 0.08352634 
    ## ... Procrustes: rmse 4.967232e-06  max resid 1.065575e-05 
    ## ... Similar to previous best
    ## Run 276 stress 0.08352634 
    ## ... Procrustes: rmse 3.086117e-06  max resid 6.671784e-06 
    ## ... Similar to previous best
    ## Run 277 stress 0.08352634 
    ## ... Procrustes: rmse 1.237049e-06  max resid 2.709951e-06 
    ## ... Similar to previous best
    ## Run 278 stress 0.08352634 
    ## ... Procrustes: rmse 1.284246e-06  max resid 2.836938e-06 
    ## ... Similar to previous best
    ## Run 279 stress 0.08352634 
    ## ... Procrustes: rmse 8.856652e-07  max resid 1.979809e-06 
    ## ... Similar to previous best
    ## Run 280 stress 0.1279309 
    ## Run 281 stress 0.08352634 
    ## ... Procrustes: rmse 3.333373e-06  max resid 7.458347e-06 
    ## ... Similar to previous best
    ## Run 282 stress 0.08352634 
    ## ... Procrustes: rmse 2.02608e-06  max resid 4.624011e-06 
    ## ... Similar to previous best
    ## Run 283 stress 0.1134014 
    ## Run 284 stress 0.1134015 
    ## Run 285 stress 0.1279309 
    ## Run 286 stress 0.1279309 
    ## Run 287 stress 0.08352634 
    ## ... Procrustes: rmse 2.409611e-06  max resid 5.330679e-06 
    ## ... Similar to previous best
    ## Run 288 stress 0.08352634 
    ## ... Procrustes: rmse 3.061359e-06  max resid 6.681138e-06 
    ## ... Similar to previous best
    ## Run 289 stress 0.1279309 
    ## Run 290 stress 0.08352634 
    ## ... Procrustes: rmse 6.501633e-07  max resid 1.546642e-06 
    ## ... Similar to previous best
    ## Run 291 stress 0.1989486 
    ## Run 292 stress 0.08352634 
    ## ... Procrustes: rmse 1.993508e-06  max resid 4.469464e-06 
    ## ... Similar to previous best
    ## Run 293 stress 0.08352634 
    ## ... Procrustes: rmse 2.439984e-06  max resid 5.136172e-06 
    ## ... Similar to previous best
    ## Run 294 stress 0.1134015 
    ## Run 295 stress 0.1134014 
    ## Run 296 stress 0.08352634 
    ## ... Procrustes: rmse 4.988704e-07  max resid 8.01223e-07 
    ## ... Similar to previous best
    ## Run 297 stress 0.08352634 
    ## ... Procrustes: rmse 3.646469e-06  max resid 8.219344e-06 
    ## ... Similar to previous best
    ## Run 298 stress 0.08352634 
    ## ... Procrustes: rmse 7.031161e-07  max resid 1.397105e-06 
    ## ... Similar to previous best
    ## Run 299 stress 0.08352634 
    ## ... Procrustes: rmse 5.543332e-06  max resid 1.248863e-05 
    ## ... Similar to previous best
    ## Run 300 stress 0.1134015 
    ## Run 301 stress 0.1279309 
    ## Run 302 stress 0.1279309 
    ## Run 303 stress 0.08352634 
    ## ... Procrustes: rmse 2.229102e-06  max resid 5.284174e-06 
    ## ... Similar to previous best
    ## Run 304 stress 0.08352634 
    ## ... Procrustes: rmse 1.402948e-06  max resid 3.337359e-06 
    ## ... Similar to previous best
    ## Run 305 stress 0.08352634 
    ## ... Procrustes: rmse 2.422773e-06  max resid 5.433786e-06 
    ## ... Similar to previous best
    ## Run 306 stress 0.08352634 
    ## ... Procrustes: rmse 8.083385e-06  max resid 1.818214e-05 
    ## ... Similar to previous best
    ## Run 307 stress 0.08352634 
    ## ... Procrustes: rmse 5.580209e-07  max resid 1.127357e-06 
    ## ... Similar to previous best
    ## Run 308 stress 0.1279309 
    ## Run 309 stress 0.08352634 
    ## ... Procrustes: rmse 2.284733e-06  max resid 4.985427e-06 
    ## ... Similar to previous best
    ## Run 310 stress 0.08352634 
    ## ... Procrustes: rmse 1.726461e-06  max resid 3.46222e-06 
    ## ... Similar to previous best
    ## Run 311 stress 0.08352634 
    ## ... Procrustes: rmse 2.109324e-06  max resid 3.226595e-06 
    ## ... Similar to previous best
    ## Run 312 stress 0.08352634 
    ## ... Procrustes: rmse 1.961769e-06  max resid 4.27806e-06 
    ## ... Similar to previous best
    ## Run 313 stress 0.08352634 
    ## ... Procrustes: rmse 2.349633e-06  max resid 4.723061e-06 
    ## ... Similar to previous best
    ## Run 314 stress 0.1989156 
    ## Run 315 stress 0.08352634 
    ## ... Procrustes: rmse 4.713018e-07  max resid 8.961236e-07 
    ## ... Similar to previous best
    ## Run 316 stress 0.08352634 
    ## ... Procrustes: rmse 1.432352e-06  max resid 2.018052e-06 
    ## ... Similar to previous best
    ## Run 317 stress 0.1134014 
    ## Run 318 stress 0.08352634 
    ## ... Procrustes: rmse 5.367239e-07  max resid 1.059547e-06 
    ## ... Similar to previous best
    ## Run 319 stress 0.1989154 
    ## Run 320 stress 0.08352635 
    ## ... Procrustes: rmse 8.588804e-06  max resid 1.346241e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.08352634 
    ## ... Procrustes: rmse 2.357018e-06  max resid 5.257375e-06 
    ## ... Similar to previous best
    ## Run 322 stress 0.1279309 
    ## Run 323 stress 0.1279309 
    ## Run 324 stress 0.1134014 
    ## Run 325 stress 0.08352634 
    ## ... Procrustes: rmse 2.041059e-06  max resid 4.061512e-06 
    ## ... Similar to previous best
    ## Run 326 stress 0.1134015 
    ## Run 327 stress 0.08352634 
    ## ... Procrustes: rmse 1.70033e-06  max resid 3.697648e-06 
    ## ... Similar to previous best
    ## Run 328 stress 0.08352634 
    ## ... Procrustes: rmse 6.192253e-06  max resid 1.411945e-05 
    ## ... Similar to previous best
    ## Run 329 stress 0.08352634 
    ## ... Procrustes: rmse 6.554585e-06  max resid 1.312625e-05 
    ## ... Similar to previous best
    ## Run 330 stress 0.08352634 
    ## ... Procrustes: rmse 2.106106e-06  max resid 4.669568e-06 
    ## ... Similar to previous best
    ## Run 331 stress 0.1279309 
    ## Run 332 stress 0.1279309 
    ## Run 333 stress 0.1279309 
    ## Run 334 stress 0.08352634 
    ## ... Procrustes: rmse 2.628136e-06  max resid 5.842636e-06 
    ## ... Similar to previous best
    ## Run 335 stress 0.08352634 
    ## ... Procrustes: rmse 2.076051e-06  max resid 4.667347e-06 
    ## ... Similar to previous best
    ## Run 336 stress 0.08352634 
    ## ... Procrustes: rmse 4.856569e-06  max resid 1.097017e-05 
    ## ... Similar to previous best
    ## Run 337 stress 0.1279309 
    ## Run 338 stress 0.08352634 
    ## ... Procrustes: rmse 1.911382e-06  max resid 4.476434e-06 
    ## ... Similar to previous best
    ## Run 339 stress 0.3154684 
    ## Run 340 stress 0.08352634 
    ## ... Procrustes: rmse 1.083778e-06  max resid 2.143205e-06 
    ## ... Similar to previous best
    ## Run 341 stress 0.08352634 
    ## ... Procrustes: rmse 3.010648e-06  max resid 6.788927e-06 
    ## ... Similar to previous best
    ## Run 342 stress 0.08352634 
    ## ... Procrustes: rmse 1.603463e-06  max resid 3.491404e-06 
    ## ... Similar to previous best
    ## Run 343 stress 0.08352634 
    ## ... Procrustes: rmse 1.793927e-06  max resid 4.068951e-06 
    ## ... Similar to previous best
    ## Run 344 stress 0.08352634 
    ## ... Procrustes: rmse 1.119888e-06  max resid 2.500927e-06 
    ## ... Similar to previous best
    ## Run 345 stress 0.08352634 
    ## ... Procrustes: rmse 1.431435e-06  max resid 3.166745e-06 
    ## ... Similar to previous best
    ## Run 346 stress 0.1989155 
    ## Run 347 stress 0.1134015 
    ## Run 348 stress 0.08352634 
    ## ... Procrustes: rmse 1.073758e-06  max resid 1.933885e-06 
    ## ... Similar to previous best
    ## Run 349 stress 0.08352634 
    ## ... Procrustes: rmse 1.305463e-06  max resid 2.210796e-06 
    ## ... Similar to previous best
    ## Run 350 stress 0.08352634 
    ## ... Procrustes: rmse 2.133337e-06  max resid 4.599582e-06 
    ## ... Similar to previous best
    ## Run 351 stress 0.08352634 
    ## ... Procrustes: rmse 7.680704e-07  max resid 1.297711e-06 
    ## ... Similar to previous best
    ## Run 352 stress 0.08352634 
    ## ... Procrustes: rmse 2.776506e-06  max resid 6.123803e-06 
    ## ... Similar to previous best
    ## Run 353 stress 0.08352634 
    ## ... Procrustes: rmse 3.064456e-06  max resid 6.746683e-06 
    ## ... Similar to previous best
    ## Run 354 stress 0.08352634 
    ## ... Procrustes: rmse 6.700102e-07  max resid 1.254855e-06 
    ## ... Similar to previous best
    ## Run 355 stress 0.1989487 
    ## Run 356 stress 0.08352634 
    ## ... Procrustes: rmse 3.561434e-06  max resid 7.871252e-06 
    ## ... Similar to previous best
    ## Run 357 stress 0.08352634 
    ## ... Procrustes: rmse 2.134032e-06  max resid 4.575838e-06 
    ## ... Similar to previous best
    ## Run 358 stress 0.08352634 
    ## ... Procrustes: rmse 5.053429e-07  max resid 1.068198e-06 
    ## ... Similar to previous best
    ## Run 359 stress 0.3154694 
    ## Run 360 stress 0.1134015 
    ## Run 361 stress 0.08352634 
    ## ... Procrustes: rmse 2.445949e-06  max resid 3.804562e-06 
    ## ... Similar to previous best
    ## Run 362 stress 0.08352634 
    ## ... Procrustes: rmse 3.142707e-06  max resid 7.049961e-06 
    ## ... Similar to previous best
    ## Run 363 stress 0.1134014 
    ## Run 364 stress 0.08352634 
    ## ... Procrustes: rmse 2.329335e-06  max resid 5.108515e-06 
    ## ... Similar to previous best
    ## Run 365 stress 0.08352634 
    ## ... Procrustes: rmse 4.042916e-06  max resid 9.11392e-06 
    ## ... Similar to previous best
    ## Run 366 stress 0.1134015 
    ## Run 367 stress 0.08352634 
    ## ... Procrustes: rmse 1.514217e-06  max resid 3.218558e-06 
    ## ... Similar to previous best
    ## Run 368 stress 0.08352634 
    ## ... Procrustes: rmse 2.513675e-06  max resid 5.712808e-06 
    ## ... Similar to previous best
    ## Run 369 stress 0.1134015 
    ## Run 370 stress 0.08352634 
    ## ... Procrustes: rmse 1.761829e-06  max resid 3.794562e-06 
    ## ... Similar to previous best
    ## Run 371 stress 0.08352634 
    ## ... Procrustes: rmse 5.826229e-06  max resid 1.299957e-05 
    ## ... Similar to previous best
    ## Run 372 stress 0.08352634 
    ## ... Procrustes: rmse 5.340659e-07  max resid 1.058058e-06 
    ## ... Similar to previous best
    ## Run 373 stress 0.08352634 
    ## ... Procrustes: rmse 2.997808e-06  max resid 6.712074e-06 
    ## ... Similar to previous best
    ## Run 374 stress 0.08352634 
    ## ... Procrustes: rmse 1.752582e-06  max resid 3.863647e-06 
    ## ... Similar to previous best
    ## Run 375 stress 0.08352634 
    ## ... Procrustes: rmse 1.895394e-06  max resid 4.159846e-06 
    ## ... Similar to previous best
    ## Run 376 stress 0.08352634 
    ## ... Procrustes: rmse 2.032228e-06  max resid 4.279575e-06 
    ## ... Similar to previous best
    ## Run 377 stress 0.08352634 
    ## ... Procrustes: rmse 7.972406e-07  max resid 1.702155e-06 
    ## ... Similar to previous best
    ## Run 378 stress 0.08352634 
    ## ... Procrustes: rmse 2.464993e-06  max resid 4.92209e-06 
    ## ... Similar to previous best
    ## Run 379 stress 0.1134014 
    ## Run 380 stress 0.08352634 
    ## ... Procrustes: rmse 3.32051e-06  max resid 7.035115e-06 
    ## ... Similar to previous best
    ## Run 381 stress 0.1279309 
    ## Run 382 stress 0.1134015 
    ## Run 383 stress 0.1279309 
    ## Run 384 stress 0.08352634 
    ## ... Procrustes: rmse 1.048929e-06  max resid 2.171152e-06 
    ## ... Similar to previous best
    ## Run 385 stress 0.08352634 
    ## ... Procrustes: rmse 9.540235e-07  max resid 2.278195e-06 
    ## ... Similar to previous best
    ## Run 386 stress 0.1279309 
    ## Run 387 stress 0.08352634 
    ## ... Procrustes: rmse 1.204332e-06  max resid 2.482658e-06 
    ## ... Similar to previous best
    ## Run 388 stress 0.08352634 
    ## ... Procrustes: rmse 1.718122e-06  max resid 3.784203e-06 
    ## ... Similar to previous best
    ## Run 389 stress 0.08352634 
    ## ... Procrustes: rmse 1.906136e-06  max resid 4.272063e-06 
    ## ... Similar to previous best
    ## Run 390 stress 0.1134014 
    ## Run 391 stress 0.08352634 
    ## ... Procrustes: rmse 5.586395e-06  max resid 1.202812e-05 
    ## ... Similar to previous best
    ## Run 392 stress 0.08352634 
    ## ... Procrustes: rmse 5.307864e-07  max resid 8.180323e-07 
    ## ... Similar to previous best
    ## Run 393 stress 0.1134015 
    ## Run 394 stress 0.08352634 
    ## ... Procrustes: rmse 2.608593e-06  max resid 5.684859e-06 
    ## ... Similar to previous best
    ## Run 395 stress 0.08352634 
    ## ... Procrustes: rmse 1.334859e-06  max resid 2.353386e-06 
    ## ... Similar to previous best
    ## Run 396 stress 0.08352634 
    ## ... Procrustes: rmse 7.64765e-07  max resid 1.564795e-06 
    ## ... Similar to previous best
    ## Run 397 stress 0.08352634 
    ## ... Procrustes: rmse 3.308421e-06  max resid 7.313717e-06 
    ## ... Similar to previous best
    ## Run 398 stress 0.1134014 
    ## Run 399 stress 0.1134014 
    ## Run 400 stress 0.3407343 
    ## Run 401 stress 0.08352634 
    ## ... Procrustes: rmse 8.928497e-07  max resid 1.455766e-06 
    ## ... Similar to previous best
    ## Run 402 stress 0.1134015 
    ## Run 403 stress 0.08352634 
    ## ... Procrustes: rmse 5.832429e-07  max resid 1.229629e-06 
    ## ... Similar to previous best
    ## Run 404 stress 0.08352634 
    ## ... Procrustes: rmse 4.908428e-06  max resid 1.111545e-05 
    ## ... Similar to previous best
    ## Run 405 stress 0.3407337 
    ## Run 406 stress 0.1279309 
    ## Run 407 stress 0.1134014 
    ## Run 408 stress 0.08352634 
    ## ... Procrustes: rmse 2.88084e-06  max resid 6.491602e-06 
    ## ... Similar to previous best
    ## Run 409 stress 0.127931 
    ## Run 410 stress 0.1134015 
    ## Run 411 stress 0.08352634 
    ## ... Procrustes: rmse 1.219473e-06  max resid 2.684356e-06 
    ## ... Similar to previous best
    ## Run 412 stress 0.08352634 
    ## ... Procrustes: rmse 3.048394e-06  max resid 6.782907e-06 
    ## ... Similar to previous best
    ## Run 413 stress 0.08352634 
    ## ... Procrustes: rmse 2.966161e-06  max resid 6.543286e-06 
    ## ... Similar to previous best
    ## Run 414 stress 0.1134014 
    ## Run 415 stress 0.08352634 
    ## ... Procrustes: rmse 2.029643e-06  max resid 4.495591e-06 
    ## ... Similar to previous best
    ## Run 416 stress 0.1134014 
    ## Run 417 stress 0.08352634 
    ## ... Procrustes: rmse 3.859286e-06  max resid 8.463586e-06 
    ## ... Similar to previous best
    ## Run 418 stress 0.08352634 
    ## ... Procrustes: rmse 3.027477e-06  max resid 6.730466e-06 
    ## ... Similar to previous best
    ## Run 419 stress 0.1134014 
    ## Run 420 stress 0.1134014 
    ## Run 421 stress 0.08352634 
    ## ... Procrustes: rmse 2.04229e-06  max resid 4.603421e-06 
    ## ... Similar to previous best
    ## Run 422 stress 0.08352634 
    ## ... Procrustes: rmse 1.914764e-06  max resid 4.077878e-06 
    ## ... Similar to previous best
    ## Run 423 stress 0.127931 
    ## Run 424 stress 0.08352634 
    ## ... Procrustes: rmse 3.527117e-07  max resid 6.805092e-07 
    ## ... Similar to previous best
    ## Run 425 stress 0.1134015 
    ## Run 426 stress 0.1279309 
    ## Run 427 stress 0.1134015 
    ## Run 428 stress 0.08352634 
    ## ... Procrustes: rmse 4.96864e-06  max resid 1.129059e-05 
    ## ... Similar to previous best
    ## Run 429 stress 0.08352634 
    ## ... Procrustes: rmse 1.159421e-06  max resid 2.539901e-06 
    ## ... Similar to previous best
    ## Run 430 stress 0.08352634 
    ## ... Procrustes: rmse 1.60972e-06  max resid 3.89629e-06 
    ## ... Similar to previous best
    ## Run 431 stress 0.1134015 
    ## Run 432 stress 0.1134015 
    ## Run 433 stress 0.08352634 
    ## ... Procrustes: rmse 1.222696e-06  max resid 2.639995e-06 
    ## ... Similar to previous best
    ## Run 434 stress 0.1134015 
    ## Run 435 stress 0.3275377 
    ## Run 436 stress 0.1989484 
    ## Run 437 stress 0.08352634 
    ## ... Procrustes: rmse 4.341841e-07  max resid 9.683561e-07 
    ## ... Similar to previous best
    ## Run 438 stress 0.08352634 
    ## ... Procrustes: rmse 1.631889e-06  max resid 3.636591e-06 
    ## ... Similar to previous best
    ## Run 439 stress 0.08352634 
    ## ... Procrustes: rmse 3.092661e-06  max resid 7.012954e-06 
    ## ... Similar to previous best
    ## Run 440 stress 0.08352634 
    ## ... Procrustes: rmse 1.565705e-06  max resid 3.333633e-06 
    ## ... Similar to previous best
    ## Run 441 stress 0.08352634 
    ## ... Procrustes: rmse 2.831146e-06  max resid 6.186845e-06 
    ## ... Similar to previous best
    ## Run 442 stress 0.08352634 
    ## ... Procrustes: rmse 3.380827e-06  max resid 5.344787e-06 
    ## ... Similar to previous best
    ## Run 443 stress 0.3154675 
    ## Run 444 stress 0.1279309 
    ## Run 445 stress 0.08352634 
    ## ... Procrustes: rmse 1.371003e-06  max resid 2.887821e-06 
    ## ... Similar to previous best
    ## Run 446 stress 0.08352634 
    ## ... Procrustes: rmse 4.680845e-06  max resid 1.049483e-05 
    ## ... Similar to previous best
    ## Run 447 stress 0.08352634 
    ## ... Procrustes: rmse 3.16261e-07  max resid 5.116204e-07 
    ## ... Similar to previous best
    ## Run 448 stress 0.1279309 
    ## Run 449 stress 0.08352634 
    ## ... Procrustes: rmse 6.482005e-06  max resid 1.425789e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.1134014 
    ## Run 451 stress 0.08352634 
    ## ... Procrustes: rmse 2.532464e-06  max resid 4.713047e-06 
    ## ... Similar to previous best
    ## Run 452 stress 0.1134014 
    ## Run 453 stress 0.1134014 
    ## Run 454 stress 0.1134015 
    ## Run 455 stress 0.1134014 
    ## Run 456 stress 0.08352634 
    ## ... Procrustes: rmse 6.564357e-07  max resid 1.088657e-06 
    ## ... Similar to previous best
    ## Run 457 stress 0.08352634 
    ## ... Procrustes: rmse 2.095845e-06  max resid 4.585799e-06 
    ## ... Similar to previous best
    ## Run 458 stress 0.08352634 
    ## ... Procrustes: rmse 2.112303e-06  max resid 4.656272e-06 
    ## ... Similar to previous best
    ## Run 459 stress 0.1134014 
    ## Run 460 stress 0.08352634 
    ## ... Procrustes: rmse 1.236201e-06  max resid 2.185809e-06 
    ## ... Similar to previous best
    ## Run 461 stress 0.08352634 
    ## ... Procrustes: rmse 3.881329e-07  max resid 5.837917e-07 
    ## ... Similar to previous best
    ## Run 462 stress 0.08352634 
    ## ... Procrustes: rmse 2.923339e-06  max resid 6.399791e-06 
    ## ... Similar to previous best
    ## Run 463 stress 0.08352634 
    ## ... Procrustes: rmse 1.545632e-06  max resid 3.453603e-06 
    ## ... Similar to previous best
    ## Run 464 stress 0.08352634 
    ## ... Procrustes: rmse 3.710631e-07  max resid 7.26174e-07 
    ## ... Similar to previous best
    ## Run 465 stress 0.08352634 
    ## ... Procrustes: rmse 8.245976e-07  max resid 1.833924e-06 
    ## ... Similar to previous best
    ## Run 466 stress 0.08352634 
    ## ... Procrustes: rmse 1.634392e-06  max resid 3.59028e-06 
    ## ... Similar to previous best
    ## Run 467 stress 0.1134015 
    ## Run 468 stress 0.08352634 
    ## ... Procrustes: rmse 5.241917e-06  max resid 8.450316e-06 
    ## ... Similar to previous best
    ## Run 469 stress 0.08352634 
    ## ... Procrustes: rmse 1.676412e-06  max resid 3.582077e-06 
    ## ... Similar to previous best
    ## Run 470 stress 0.08352634 
    ## ... Procrustes: rmse 2.089942e-06  max resid 4.770391e-06 
    ## ... Similar to previous best
    ## Run 471 stress 0.08352634 
    ## ... Procrustes: rmse 1.968174e-06  max resid 4.280162e-06 
    ## ... Similar to previous best
    ## Run 472 stress 0.1134014 
    ## Run 473 stress 0.1134014 
    ## Run 474 stress 0.1279309 
    ## Run 475 stress 0.08352634 
    ## ... Procrustes: rmse 1.087757e-06  max resid 2.231929e-06 
    ## ... Similar to previous best
    ## Run 476 stress 0.1279309 
    ## Run 477 stress 0.1134014 
    ## Run 478 stress 0.08352634 
    ## ... Procrustes: rmse 7.171863e-07  max resid 1.598236e-06 
    ## ... Similar to previous best
    ## Run 479 stress 0.1134015 
    ## Run 480 stress 0.1134014 
    ## Run 481 stress 0.1989486 
    ## Run 482 stress 0.1134014 
    ## Run 483 stress 0.08352634 
    ## ... Procrustes: rmse 2.96339e-06  max resid 6.453339e-06 
    ## ... Similar to previous best
    ## Run 484 stress 0.1134014 
    ## Run 485 stress 0.1989487 
    ## Run 486 stress 0.08352634 
    ## ... Procrustes: rmse 1.070859e-06  max resid 1.982108e-06 
    ## ... Similar to previous best
    ## Run 487 stress 0.08352634 
    ## ... Procrustes: rmse 1.56606e-06  max resid 3.47947e-06 
    ## ... Similar to previous best
    ## Run 488 stress 0.1279309 
    ## Run 489 stress 0.08352634 
    ## ... Procrustes: rmse 1.718022e-06  max resid 3.609767e-06 
    ## ... Similar to previous best
    ## Run 490 stress 0.08352634 
    ## ... Procrustes: rmse 1.96528e-06  max resid 4.186424e-06 
    ## ... Similar to previous best
    ## Run 491 stress 0.1134015 
    ## Run 492 stress 0.08352634 
    ## ... Procrustes: rmse 1.51333e-06  max resid 3.158249e-06 
    ## ... Similar to previous best
    ## Run 493 stress 0.08352634 
    ## ... Procrustes: rmse 2.028932e-06  max resid 3.740094e-06 
    ## ... Similar to previous best
    ## Run 494 stress 0.3266727 
    ## Run 495 stress 0.1134015 
    ## Run 496 stress 0.08352634 
    ## ... Procrustes: rmse 2.660602e-06  max resid 5.961753e-06 
    ## ... Similar to previous best
    ## Run 497 stress 0.08352634 
    ## ... Procrustes: rmse 5.725463e-06  max resid 1.287551e-05 
    ## ... Similar to previous best
    ## Run 498 stress 0.08352634 
    ## ... Procrustes: rmse 1.799914e-06  max resid 4.037109e-06 
    ## ... Similar to previous best
    ## Run 499 stress 0.1134014 
    ## Run 500 stress 0.08352634 
    ## ... Procrustes: rmse 2.951753e-06  max resid 6.522472e-06 
    ## ... Similar to previous best
    ## *** Best solution repeated 163 times

``` r
# Stratified lakes and ocean sites
PD_beta_env_SO_NMDS <- metaMDS(PD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.92006e-05 
    ## Run 1 stress 0.0002997702 
    ## ... Procrustes: rmse 0.002266303  max resid 0.003413405 
    ## ... Similar to previous best
    ## Run 2 stress 9.867324e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002284394  max resid 0.0004624251 
    ## ... Similar to previous best
    ## Run 3 stress 9.85232e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002008334  max resid 0.0004663151 
    ## ... Similar to previous best
    ## Run 4 stress 9.551842e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002490729  max resid 0.0006060279 
    ## ... Similar to previous best
    ## Run 5 stress 9.928493e-05 
    ## ... Procrustes: rmse 0.0002076655  max resid 0.0005946816 
    ## ... Similar to previous best
    ## Run 6 stress 9.649181e-05 
    ## ... Procrustes: rmse 0.000246359  max resid 0.000598751 
    ## ... Similar to previous best
    ## Run 7 stress 9.832155e-05 
    ## ... Procrustes: rmse 0.0002492537  max resid 0.0006991443 
    ## ... Similar to previous best
    ## Run 8 stress 9.800814e-05 
    ## ... Procrustes: rmse 0.0003052213  max resid 0.0007042117 
    ## ... Similar to previous best
    ## Run 9 stress 9.85377e-05 
    ## ... Procrustes: rmse 0.000306202  max resid 0.0007065398 
    ## ... Similar to previous best
    ## Run 10 stress 0.0004378958 
    ## ... Procrustes: rmse 0.003568134  max resid 0.005653862 
    ## ... Similar to previous best
    ## Run 11 stress 9.569368e-05 
    ## ... Procrustes: rmse 6.833533e-06  max resid 1.54846e-05 
    ## ... Similar to previous best
    ## Run 12 stress 9.916094e-05 
    ## ... Procrustes: rmse 0.0002502533  max resid 0.0006092527 
    ## ... Similar to previous best
    ## Run 13 stress 0.3192024 
    ## Run 14 stress 9.762884e-05 
    ## ... Procrustes: rmse 0.0002664144  max resid 0.0006815758 
    ## ... Similar to previous best
    ## Run 15 stress 9.926169e-05 
    ## ... Procrustes: rmse 0.0003180535  max resid 0.0006626717 
    ## ... Similar to previous best
    ## Run 16 stress 0.0004152576 
    ## ... Procrustes: rmse 0.003370449  max resid 0.00534194 
    ## ... Similar to previous best
    ## Run 17 stress 0.0008641716 
    ## Run 18 stress 9.923569e-05 
    ## ... Procrustes: rmse 2.124816e-05  max resid 3.287656e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.3403702 
    ## Run 20 stress 9.956811e-05 
    ## ... Procrustes: rmse 0.0003072591  max resid 0.000709256 
    ## ... Similar to previous best
    ## Run 21 stress 0.001042378 
    ## Run 22 stress 0.3407339 
    ## Run 23 stress 9.908002e-05 
    ## ... Procrustes: rmse 0.0002487699  max resid 0.0006063597 
    ## ... Similar to previous best
    ## Run 24 stress 0.001201385 
    ## Run 25 stress 0.001572228 
    ## Run 26 stress 9.90808e-05 
    ## ... Procrustes: rmse 0.0002509434  max resid 0.00070176 
    ## ... Similar to previous best
    ## Run 27 stress 9.883334e-05 
    ## ... Procrustes: rmse 1.963292e-05  max resid 3.227644e-05 
    ## ... Similar to previous best
    ## Run 28 stress 0.3275198 
    ## Run 29 stress 9.90718e-05 
    ## ... Procrustes: rmse 0.0003256309  max resid 0.0006735184 
    ## ... Similar to previous best
    ## Run 30 stress 0.0007470492 
    ## Run 31 stress 9.748551e-05 
    ## ... Procrustes: rmse 0.0003221195  max resid 0.0006675132 
    ## ... Similar to previous best
    ## Run 32 stress 9.988684e-05 
    ## ... Procrustes: rmse 2.360342e-05  max resid 3.988482e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.0004666914 
    ## ... Procrustes: rmse 0.003817879  max resid 0.006040565 
    ## ... Similar to previous best
    ## Run 34 stress 9.875288e-05 
    ## ... Procrustes: rmse 0.000325076  max resid 0.0006725388 
    ## ... Similar to previous best
    ## Run 35 stress 0.0007690934 
    ## Run 36 stress 9.972704e-05 
    ## ... Procrustes: rmse 0.000269879  max resid 0.0007103339 
    ## ... Similar to previous best
    ## Run 37 stress 9.920188e-05 
    ## ... Procrustes: rmse 0.0002520095  max resid 0.0007042857 
    ## ... Similar to previous best
    ## Run 38 stress 0.001100242 
    ## Run 39 stress 0.0009998608 
    ## Run 40 stress 0.001118707 
    ## Run 41 stress 9.859784e-05 
    ## ... Procrustes: rmse 0.0002594821  max resid 0.0006691746 
    ## ... Similar to previous best
    ## Run 42 stress 9.88305e-05 
    ## ... Procrustes: rmse 0.0003243714  max resid 0.0006688841 
    ## ... Similar to previous best
    ## Run 43 stress 9.857723e-05 
    ## ... Procrustes: rmse 0.0001993043  max resid 0.0005763753 
    ## ... Similar to previous best
    ## Run 44 stress 9.271918e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004087714  max resid 0.0005999292 
    ## ... Similar to previous best
    ## Run 45 stress 9.52944e-05 
    ## ... Procrustes: rmse 0.0004290527  max resid 0.0006894616 
    ## ... Similar to previous best
    ## Run 46 stress 9.800715e-05 
    ## ... Procrustes: rmse 0.0004575853  max resid 0.0007420133 
    ## ... Similar to previous best
    ## Run 47 stress 9.950506e-05 
    ## ... Procrustes: rmse 0.0004275103  max resid 0.0006194071 
    ## ... Similar to previous best
    ## Run 48 stress 0.0001477598 
    ## ... Procrustes: rmse 0.001290662  max resid 0.001888092 
    ## ... Similar to previous best
    ## Run 49 stress 0.0004160645 
    ## ... Procrustes: rmse 0.003659719  max resid 0.00553061 
    ## ... Similar to previous best
    ## Run 50 stress 9.972549e-05 
    ## ... Procrustes: rmse 0.0004820847  max resid 0.0007836591 
    ## ... Similar to previous best
    ## Run 51 stress 0.0007902969 
    ## Run 52 stress 0.0006208539 
    ## Run 53 stress 9.856809e-05 
    ## ... Procrustes: rmse 0.0004464581  max resid 0.0007194622 
    ## ... Similar to previous best
    ## Run 54 stress 9.63754e-05 
    ## ... Procrustes: rmse 0.0004156088  max resid 0.0006038928 
    ## ... Similar to previous best
    ## Run 55 stress 9.796795e-05 
    ## ... Procrustes: rmse 0.0004446059  max resid 0.0007167156 
    ## ... Similar to previous best
    ## Run 56 stress 9.839294e-05 
    ## ... Procrustes: rmse 0.0004474202  max resid 0.0006972975 
    ## ... Similar to previous best
    ## Run 57 stress 9.646947e-05 
    ## ... Procrustes: rmse 0.000494836  max resid 0.0007545682 
    ## ... Similar to previous best
    ## Run 58 stress 9.756396e-05 
    ## ... Procrustes: rmse 0.0004213075  max resid 0.0006986806 
    ## ... Similar to previous best
    ## Run 59 stress 9.921291e-05 
    ## ... Procrustes: rmse 0.0004505098  max resid 0.0007062538 
    ## ... Similar to previous best
    ## Run 60 stress 0.000831307 
    ## Run 61 stress 9.900557e-05 
    ## ... Procrustes: rmse 0.0004266653  max resid 0.0006181193 
    ## ... Similar to previous best
    ## Run 62 stress 0.0005770542 
    ## ... Procrustes: rmse 0.005060644  max resid 0.007707869 
    ## Run 63 stress 0.001164626 
    ## Run 64 stress 0.001095486 
    ## Run 65 stress 9.870447e-05 
    ## ... Procrustes: rmse 0.0004467813  max resid 0.0007192944 
    ## ... Similar to previous best
    ## Run 66 stress 0.0005205793 
    ## ... Procrustes: rmse 0.004573921  max resid 0.006946367 
    ## ... Similar to previous best
    ## Run 67 stress 9.844031e-05 
    ## ... Procrustes: rmse 0.0004318745  max resid 0.0007176629 
    ## ... Similar to previous best
    ## Run 68 stress 9.445675e-05 
    ## ... Procrustes: rmse 1.452747e-05  max resid 2.51811e-05 
    ## ... Similar to previous best
    ## Run 69 stress 0.001069139 
    ## Run 70 stress 9.941929e-05 
    ## ... Procrustes: rmse 0.0004268441  max resid 0.0006217544 
    ## ... Similar to previous best
    ## Run 71 stress 9.969551e-05 
    ## ... Procrustes: rmse 0.0004364712  max resid 0.0007276425 
    ## ... Similar to previous best
    ## Run 72 stress 0.0010412 
    ## Run 73 stress 9.968994e-05 
    ## ... Procrustes: rmse 0.000453471  max resid 0.0007092513 
    ## ... Similar to previous best
    ## Run 74 stress 9.901818e-05 
    ## ... Procrustes: rmse 0.0004276848  max resid 0.0007136764 
    ## ... Similar to previous best
    ## Run 75 stress 9.865538e-05 
    ## ... Procrustes: rmse 0.0004311039  max resid 0.0007191387 
    ## ... Similar to previous best
    ## Run 76 stress 9.935097e-05 
    ## ... Procrustes: rmse 0.000435793  max resid 0.0007251169 
    ## ... Similar to previous best
    ## Run 77 stress 0.0006205841 
    ## Run 78 stress 9.75775e-05 
    ## ... Procrustes: rmse 0.0004808816  max resid 0.0007627042 
    ## ... Similar to previous best
    ## Run 79 stress 0.001107935 
    ## Run 80 stress 9.547142e-05 
    ## ... Procrustes: rmse 0.0004340215  max resid 0.0007003543 
    ## ... Similar to previous best
    ## Run 81 stress 0.001124434 
    ## Run 82 stress 0.0007819547 
    ## Run 83 stress 0.0009402724 
    ## Run 84 stress 9.988541e-05 
    ## ... Procrustes: rmse 0.0005115627  max resid 0.0007790787 
    ## ... Similar to previous best
    ## Run 85 stress 9.926679e-05 
    ## ... Procrustes: rmse 0.0004491648  max resid 0.0007094936 
    ## ... Similar to previous best
    ## Run 86 stress 9.958493e-05 
    ## ... Procrustes: rmse 0.0004428906  max resid 0.0007225135 
    ## ... Similar to previous best
    ## Run 87 stress 0.00110996 
    ## Run 88 stress 0.2629798 
    ## Run 89 stress 0.001089056 
    ## Run 90 stress 9.731584e-05 
    ## ... Procrustes: rmse 0.0002442662  max resid 0.0004542951 
    ## ... Similar to previous best
    ## Run 91 stress 0.001182482 
    ## Run 92 stress 9.906363e-05 
    ## ... Procrustes: rmse 0.0004764742  max resid 0.0007750566 
    ## ... Similar to previous best
    ## Run 93 stress 9.905832e-05 
    ## ... Procrustes: rmse 0.000460368  max resid 0.0007451463 
    ## ... Similar to previous best
    ## Run 94 stress 9.979544e-05 
    ## ... Procrustes: rmse 0.0004530656  max resid 0.0007099153 
    ## ... Similar to previous best
    ## Run 95 stress 0.0002802414 
    ## ... Procrustes: rmse 0.002469437  max resid 0.003688177 
    ## ... Similar to previous best
    ## Run 96 stress 9.950074e-05 
    ## ... Procrustes: rmse 0.0004287104  max resid 0.0006209816 
    ## ... Similar to previous best
    ## Run 97 stress 0.0007949823 
    ## Run 98 stress 8.152826e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 9.324015e-05  max resid 0.0001439385 
    ## ... Similar to previous best
    ## Run 99 stress 9.870508e-05 
    ## ... Procrustes: rmse 0.0004174271  max resid 0.0006630914 
    ## ... Similar to previous best
    ## Run 100 stress 9.827126e-05 
    ## ... Procrustes: rmse 0.0003774198  max resid 0.0006023803 
    ## ... Similar to previous best
    ## Run 101 stress 9.912815e-05 
    ## ... Procrustes: rmse 0.000443798  max resid 0.0006735502 
    ## ... Similar to previous best
    ## Run 102 stress 9.963284e-05 
    ## ... Procrustes: rmse 0.0004666367  max resid 0.0007339252 
    ## ... Similar to previous best
    ## Run 103 stress 0.0002320143 
    ## ... Procrustes: rmse 0.002046352  max resid 0.003172888 
    ## ... Similar to previous best
    ## Run 104 stress 9.233726e-05 
    ## ... Procrustes: rmse 0.0001138401  max resid 0.0002062191 
    ## ... Similar to previous best
    ## Run 105 stress 9.778097e-05 
    ## ... Procrustes: rmse 0.0004833964  max resid 0.0007113066 
    ## ... Similar to previous best
    ## Run 106 stress 9.904482e-05 
    ## ... Procrustes: rmse 0.0004405706  max resid 0.0006533977 
    ## ... Similar to previous best
    ## Run 107 stress 0.0007504212 
    ## Run 108 stress 9.905177e-05 
    ## ... Procrustes: rmse 0.0004401764  max resid 0.0006526807 
    ## ... Similar to previous best
    ## Run 109 stress 9.658537e-05 
    ## ... Procrustes: rmse 0.0004166059  max resid 0.0006371881 
    ## ... Similar to previous best
    ## Run 110 stress 0.001132278 
    ## Run 111 stress 9.983429e-05 
    ## ... Procrustes: rmse 0.0004674621  max resid 0.0007340575 
    ## ... Similar to previous best
    ## Run 112 stress 0.0008464568 
    ## Run 113 stress 9.870753e-05 
    ## ... Procrustes: rmse 0.0004390109  max resid 0.0006518987 
    ## ... Similar to previous best
    ## Run 114 stress 0.001112871 
    ## Run 115 stress 9.640617e-05 
    ## ... Procrustes: rmse 0.0004145796  max resid 0.000614535 
    ## ... Similar to previous best
    ## Run 116 stress 0.0006799169 
    ## Run 117 stress 0.0007717004 
    ## Run 118 stress 9.970431e-05 
    ## ... Procrustes: rmse 0.0004223747  max resid 0.0006724825 
    ## ... Similar to previous best
    ## Run 119 stress 0.001117314 
    ## Run 120 stress 9.830231e-05 
    ## ... Procrustes: rmse 0.0004616487  max resid 0.0007072858 
    ## ... Similar to previous best
    ## Run 121 stress 0.0007691359 
    ## Run 122 stress 9.958143e-05 
    ## ... Procrustes: rmse 0.0004662293  max resid 0.0007324449 
    ## ... Similar to previous best
    ## Run 123 stress 9.927631e-05 
    ## ... Procrustes: rmse 0.0004417867  max resid 0.0006561578 
    ## ... Similar to previous best
    ## Run 124 stress 9.941274e-05 
    ## ... Procrustes: rmse 0.0004300457  max resid 0.0006406112 
    ## ... Similar to previous best
    ## Run 125 stress 9.727714e-05 
    ## ... Procrustes: rmse 0.0004365666  max resid 0.0006656552 
    ## ... Similar to previous best
    ## Run 126 stress 0.0008114812 
    ## Run 127 stress 9.988514e-05 
    ## ... Procrustes: rmse 0.0004276498  max resid 0.0006877063 
    ## ... Similar to previous best
    ## Run 128 stress 0.0006557558 
    ## Run 129 stress 9.919068e-05 
    ## ... Procrustes: rmse 0.000462488  max resid 0.0007074218 
    ## ... Similar to previous best
    ## Run 130 stress 0.001081127 
    ## Run 131 stress 0.00117537 
    ## Run 132 stress 9.395258e-05 
    ## ... Procrustes: rmse 0.0003979008  max resid 0.000636656 
    ## ... Similar to previous best
    ## Run 133 stress 9.866406e-05 
    ## ... Procrustes: rmse 0.0004244812  max resid 0.0006284125 
    ## ... Similar to previous best
    ## Run 134 stress 9.858436e-05 
    ## ... Procrustes: rmse 0.000413072  max resid 0.0006561225 
    ## ... Similar to previous best
    ## Run 135 stress 9.923064e-05 
    ## ... Procrustes: rmse 0.0004275843  max resid 0.0006497972 
    ## ... Similar to previous best
    ## Run 136 stress 9.925664e-05 
    ## ... Procrustes: rmse 0.0004210126  max resid 0.0006700751 
    ## ... Similar to previous best
    ## Run 137 stress 9.974754e-05 
    ## ... Procrustes: rmse 0.0004216972  max resid 0.0006723338 
    ## ... Similar to previous best
    ## Run 138 stress 9.982388e-05 
    ## ... Procrustes: rmse 0.0004286901  max resid 0.00068446 
    ## ... Similar to previous best
    ## Run 139 stress 9.970743e-05 
    ## ... Procrustes: rmse 0.0004278167  max resid 0.0006348761 
    ## ... Similar to previous best
    ## Run 140 stress 9.895623e-05 
    ## ... Procrustes: rmse 0.0004238423  max resid 0.0006731986 
    ## ... Similar to previous best
    ## Run 141 stress 9.882244e-05 
    ## ... Procrustes: rmse 0.0003656245  max resid 0.0006132555 
    ## ... Similar to previous best
    ## Run 142 stress 9.60984e-05 
    ## ... Procrustes: rmse 0.0004120315  max resid 0.0006693339 
    ## ... Similar to previous best
    ## Run 143 stress 9.907981e-05 
    ## ... Procrustes: rmse 0.0004277247  max resid 0.0006530135 
    ## ... Similar to previous best
    ## Run 144 stress 9.688056e-05 
    ## ... Procrustes: rmse 0.0004174103  max resid 0.0006350999 
    ## ... Similar to previous best
    ## Run 145 stress 9.124001e-05 
    ## ... Procrustes: rmse 0.0001075203  max resid 0.0001732126 
    ## ... Similar to previous best
    ## Run 146 stress 9.95804e-05 
    ## ... Procrustes: rmse 0.0004282899  max resid 0.0006344033 
    ## ... Similar to previous best
    ## Run 147 stress 9.873284e-05 
    ## ... Procrustes: rmse 0.000425279  max resid 0.0006474553 
    ## ... Similar to previous best
    ## Run 148 stress 9.72392e-05 
    ## ... Procrustes: rmse 0.0002110999  max resid 0.0003833317 
    ## ... Similar to previous best
    ## Run 149 stress 0.001195176 
    ## Run 150 stress 9.893008e-05 
    ## ... Procrustes: rmse 0.0004861163  max resid 0.0007359766 
    ## ... Similar to previous best
    ## Run 151 stress 9.678779e-05 
    ## ... Procrustes: rmse 0.0004099831  max resid 0.0006547867 
    ## ... Similar to previous best
    ## Run 152 stress 9.589649e-05 
    ## ... Procrustes: rmse 0.0004319189  max resid 0.0006581072 
    ## ... Similar to previous best
    ## Run 153 stress 0.001148894 
    ## Run 154 stress 0.0005088039 
    ## ... Procrustes: rmse 0.004493652  max resid 0.006948657 
    ## ... Similar to previous best
    ## Run 155 stress 9.748969e-05 
    ## ... Procrustes: rmse 0.0004187222  max resid 0.0006359868 
    ## ... Similar to previous best
    ## Run 156 stress 9.196417e-05 
    ## ... Procrustes: rmse 0.0003976942  max resid 0.0006103367 
    ## ... Similar to previous best
    ## Run 157 stress 0.0007119992 
    ## Run 158 stress 9.954107e-05 
    ## ... Procrustes: rmse 0.0004278501  max resid 0.0006341949 
    ## ... Similar to previous best
    ## Run 159 stress 0.00122895 
    ## Run 160 stress 9.652388e-05 
    ## ... Procrustes: rmse 0.0004349238  max resid 0.0006633925 
    ## ... Similar to previous best
    ## Run 161 stress 9.613872e-05 
    ## ... Procrustes: rmse 0.0004134229  max resid 0.0006117755 
    ## ... Similar to previous best
    ## Run 162 stress 9.98608e-05 
    ## ... Procrustes: rmse 0.0004218331  max resid 0.0006722032 
    ## ... Similar to previous best
    ## Run 163 stress 9.851915e-05 
    ## ... Procrustes: rmse 0.0004111388  max resid 0.0006491021 
    ## ... Similar to previous best
    ## Run 164 stress 0.001316003 
    ## Run 165 stress 0.0008343638 
    ## Run 166 stress 9.998905e-05 
    ## ... Procrustes: rmse 0.0004942858  max resid 0.0007249672 
    ## ... Similar to previous best
    ## Run 167 stress 9.826019e-05 
    ## ... Procrustes: rmse 0.0004162097  max resid 0.0006610254 
    ## ... Similar to previous best
    ## Run 168 stress 9.908103e-05 
    ## ... Procrustes: rmse 0.0004251415  max resid 0.0006858691 
    ## ... Similar to previous best
    ## Run 169 stress 0.0005778943 
    ## ... Procrustes: rmse 0.005072273  max resid 0.00785583 
    ## Run 170 stress 0.0008093383 
    ## Run 171 stress 0.001160369 
    ## Run 172 stress 9.893544e-05 
    ## ... Procrustes: rmse 0.0004404433  max resid 0.0006710391 
    ## ... Similar to previous best
    ## Run 173 stress 9.647019e-05 
    ## ... Procrustes: rmse 0.000408859  max resid 0.0006558259 
    ## ... Similar to previous best
    ## Run 174 stress 0.0007220897 
    ## Run 175 stress 9.971541e-05 
    ## ... Procrustes: rmse 0.0004435974  max resid 0.0006569314 
    ## ... Similar to previous best
    ## Run 176 stress 9.851349e-05 
    ## ... Procrustes: rmse 0.0004628318  max resid 0.0007104993 
    ## ... Similar to previous best
    ## Run 177 stress 0.0005860083 
    ## Run 178 stress 0.0007189624 
    ## Run 179 stress 0.000525528 
    ## ... Procrustes: rmse 0.004614032  max resid 0.007134505 
    ## ... Similar to previous best
    ## Run 180 stress 0.0007035304 
    ## Run 181 stress 9.933505e-05 
    ## ... Procrustes: rmse 0.0004227674  max resid 0.0006279737 
    ## ... Similar to previous best
    ## Run 182 stress 0.001068784 
    ## Run 183 stress 9.76365e-05 
    ## ... Procrustes: rmse 0.0004083682  max resid 0.0006481721 
    ## ... Similar to previous best
    ## Run 184 stress 9.976154e-05 
    ## ... Procrustes: rmse 0.0004274271  max resid 0.0006296934 
    ## ... Similar to previous best
    ## Run 185 stress 9.881544e-05 
    ## ... Procrustes: rmse 0.0004397584  max resid 0.0006515649 
    ## ... Similar to previous best
    ## Run 186 stress 0.000814061 
    ## Run 187 stress 0.0006601086 
    ## Run 188 stress 0.0007809712 
    ## Run 189 stress 0.001040095 
    ## Run 190 stress 0.001203329 
    ## Run 191 stress 0.0006963512 
    ## Run 192 stress 9.881232e-05 
    ## ... Procrustes: rmse 0.0004634398  max resid 0.0007108888 
    ## ... Similar to previous best
    ## Run 193 stress 9.746876e-05 
    ## ... Procrustes: rmse 0.0004802114  max resid 0.0007057417 
    ## ... Similar to previous best
    ## Run 194 stress 0.000334439 
    ## ... Procrustes: rmse 0.002934425  max resid 0.004555156 
    ## ... Similar to previous best
    ## Run 195 stress 0.0006264581 
    ## Run 196 stress 9.874619e-05 
    ## ... Procrustes: rmse 0.0004185192  max resid 0.0006646639 
    ## ... Similar to previous best
    ## Run 197 stress 0.0002577073 
    ## ... Procrustes: rmse 0.002268135  max resid 0.003507212 
    ## ... Similar to previous best
    ## Run 198 stress 0.001278871 
    ## Run 199 stress 9.711786e-05 
    ## ... Procrustes: rmse 0.000470916  max resid 0.0007194774 
    ## ... Similar to previous best
    ## Run 200 stress 0.0005814343 
    ## ... Procrustes: rmse 0.005104001  max resid 0.007895433 
    ## Run 201 stress 0.0009998002 
    ## Run 202 stress 0.001045204 
    ## Run 203 stress 9.888719e-05 
    ## ... Procrustes: rmse 0.0004594174  max resid 0.0007048281 
    ## ... Similar to previous best
    ## Run 204 stress 0.001174926 
    ## Run 205 stress 9.214127e-05 
    ## ... Procrustes: rmse 0.0003984817  max resid 0.0006088615 
    ## ... Similar to previous best
    ## Run 206 stress 9.69704e-05 
    ## ... Procrustes: rmse 0.0004169  max resid 0.0006183522 
    ## ... Similar to previous best
    ## Run 207 stress 9.623402e-05 
    ## ... Procrustes: rmse 0.0004076051  max resid 0.0006071306 
    ## ... Similar to previous best
    ## Run 208 stress 0.0002480904 
    ## ... Procrustes: rmse 0.002188426  max resid 0.003397485 
    ## ... Similar to previous best
    ## Run 209 stress 9.89098e-05 
    ## ... Procrustes: rmse 0.0004456398  max resid 0.0006789074 
    ## ... Similar to previous best
    ## Run 210 stress 9.987631e-05 
    ## ... Procrustes: rmse 0.0004292406  max resid 0.0006884338 
    ## ... Similar to previous best
    ## Run 211 stress 0.0009592467 
    ## Run 212 stress 9.746743e-05 
    ## ... Procrustes: rmse 0.0004581949  max resid 0.0007042902 
    ## ... Similar to previous best
    ## Run 213 stress 9.597838e-05 
    ## ... Procrustes: rmse 0.00042564  max resid 0.000631891 
    ## ... Similar to previous best
    ## Run 214 stress 0.0007969716 
    ## Run 215 stress 9.871009e-05 
    ## ... Procrustes: rmse 0.0003735798  max resid 0.0005971235 
    ## ... Similar to previous best
    ## Run 216 stress 9.970805e-05 
    ## ... Procrustes: rmse 0.0004607268  max resid 0.0007183854 
    ## ... Similar to previous best
    ## Run 217 stress 9.969357e-05 
    ## ... Procrustes: rmse 0.0004449897  max resid 0.0006703179 
    ## ... Similar to previous best
    ## Run 218 stress 0.0008264818 
    ## Run 219 stress 9.911183e-05 
    ## ... Procrustes: rmse 0.0004197176  max resid 0.0006690432 
    ## ... Similar to previous best
    ## Run 220 stress 9.988952e-05 
    ## ... Procrustes: rmse 0.0004692735  max resid 0.0007193978 
    ## ... Similar to previous best
    ## Run 221 stress 9.68764e-05 
    ## ... Procrustes: rmse 0.0004250294  max resid 0.0006374562 
    ## ... Similar to previous best
    ## Run 222 stress 0.001149008 
    ## Run 223 stress 9.993154e-05 
    ## ... Procrustes: rmse 0.0004284941  max resid 0.0006840023 
    ## ... Similar to previous best
    ## Run 224 stress 9.896459e-05 
    ## ... Procrustes: rmse 0.0004250109  max resid 0.0006458321 
    ## ... Similar to previous best
    ## Run 225 stress 9.962642e-05 
    ## ... Procrustes: rmse 0.0004680823  max resid 0.0007170936 
    ## ... Similar to previous best
    ## Run 226 stress 9.872459e-05 
    ## ... Procrustes: rmse 0.0004633625  max resid 0.0007100603 
    ## ... Similar to previous best
    ## Run 227 stress 0.0005852213 
    ## Run 228 stress 9.68334e-05 
    ## ... Procrustes: rmse 0.0004019408  max resid 0.000585716 
    ## ... Similar to previous best
    ## Run 229 stress 9.991377e-05 
    ## ... Procrustes: rmse 0.0004298357  max resid 0.0006362805 
    ## ... Similar to previous best
    ## Run 230 stress 0.001068119 
    ## Run 231 stress 9.910843e-05 
    ## ... Procrustes: rmse 0.0004639708  max resid 0.0007267555 
    ## ... Similar to previous best
    ## Run 232 stress 0.000751783 
    ## Run 233 stress 0.000662803 
    ## Run 234 stress 0.0007074821 
    ## Run 235 stress 0.0007500199 
    ## Run 236 stress 0.0005511255 
    ## ... Procrustes: rmse 0.004865509  max resid 0.007519357 
    ## ... Similar to previous best
    ## Run 237 stress 0.001160339 
    ## Run 238 stress 9.859361e-05 
    ## ... Procrustes: rmse 0.0004179353  max resid 0.000665376 
    ## ... Similar to previous best
    ## Run 239 stress 9.767274e-05 
    ## ... Procrustes: rmse 0.0004183495  max resid 0.0006169005 
    ## ... Similar to previous best
    ## Run 240 stress 0.0006037056 
    ## Run 241 stress 9.027856e-05 
    ## ... Procrustes: rmse 0.000115227  max resid 0.0001677243 
    ## ... Similar to previous best
    ## Run 242 stress 9.94582e-05 
    ## ... Procrustes: rmse 0.0004662867  max resid 0.0007116698 
    ## ... Similar to previous best
    ## Run 243 stress 0.0007692459 
    ## Run 244 stress 9.921627e-05 
    ## ... Procrustes: rmse 0.0004235378  max resid 0.0006500959 
    ## ... Similar to previous best
    ## Run 245 stress 9.946268e-05 
    ## ... Procrustes: rmse 0.0004238953  max resid 0.0006310115 
    ## ... Similar to previous best
    ## Run 246 stress 0.0005957644 
    ## Run 247 stress 0.0009884572 
    ## Run 248 stress 0.0008527212 
    ## Run 249 stress 0.001100676 
    ## Run 250 stress 9.880545e-05 
    ## ... Procrustes: rmse 0.0004192783  max resid 0.0006669992 
    ## ... Similar to previous best
    ## Run 251 stress 9.529208e-05 
    ## ... Procrustes: rmse 0.0004096141  max resid 0.0006549329 
    ## ... Similar to previous best
    ## Run 252 stress 0.2348703 
    ## Run 253 stress 0.0006184911 
    ## Run 254 stress 9.89128e-05 
    ## ... Procrustes: rmse 0.0004588889  max resid 0.0007211412 
    ## ... Similar to previous best
    ## Run 255 stress 9.881034e-05 
    ## ... Procrustes: rmse 0.0004636539  max resid 0.0007087534 
    ## ... Similar to previous best
    ## Run 256 stress 0.0004433281 
    ## ... Procrustes: rmse 0.003911866  max resid 0.006040523 
    ## ... Similar to previous best
    ## Run 257 stress 9.959093e-05 
    ## ... Procrustes: rmse 0.0004651052  max resid 0.0007321078 
    ## ... Similar to previous best
    ## Run 258 stress 9.88102e-05 
    ## ... Procrustes: rmse 0.0004397803  max resid 0.0006528028 
    ## ... Similar to previous best
    ## Run 259 stress 0.0004305467 
    ## ... Procrustes: rmse 0.003785278  max resid 0.005854361 
    ## ... Similar to previous best
    ## Run 260 stress 0.0006307127 
    ## Run 261 stress 9.877494e-05 
    ## ... Procrustes: rmse 0.0004447816  max resid 0.0006776863 
    ## ... Similar to previous best
    ## Run 262 stress 0.0002767346 
    ## ... Procrustes: rmse 0.002438304  max resid 0.003765453 
    ## ... Similar to previous best
    ## Run 263 stress 9.656884e-05 
    ## ... Procrustes: rmse 0.0004155029  max resid 0.0006646457 
    ## ... Similar to previous best
    ## Run 264 stress 0.001213823 
    ## Run 265 stress 9.886925e-05 
    ## ... Procrustes: rmse 0.0004397661  max resid 0.0006533766 
    ## ... Similar to previous best
    ## Run 266 stress 9.725682e-05 
    ## ... Procrustes: rmse 0.0004163596  max resid 0.0006187524 
    ## ... Similar to previous best
    ## Run 267 stress 9.874444e-05 
    ## ... Procrustes: rmse 0.0004263389  max resid 0.0006504665 
    ## ... Similar to previous best
    ## Run 268 stress 9.775131e-05 
    ## ... Procrustes: rmse 0.0004394964  max resid 0.0006681818 
    ## ... Similar to previous best
    ## Run 269 stress 0.000658288 
    ## Run 270 stress 9.776433e-05 
    ## ... Procrustes: rmse 0.0004689698  max resid 0.0006978272 
    ## ... Similar to previous best
    ## Run 271 stress 9.782588e-05 
    ## ... Procrustes: rmse 0.0004157704  max resid 0.0006584308 
    ## ... Similar to previous best
    ## Run 272 stress 9.821751e-05 
    ## ... Procrustes: rmse 0.000483025  max resid 0.0007121206 
    ## ... Similar to previous best
    ## Run 273 stress 9.897967e-05 
    ## ... Procrustes: rmse 0.0004405931  max resid 0.000655053 
    ## ... Similar to previous best
    ## Run 274 stress 9.831404e-05 
    ## ... Procrustes: rmse 0.0004209096  max resid 0.0006199272 
    ## ... Similar to previous best
    ## Run 275 stress 9.986463e-05 
    ## ... Procrustes: rmse 0.0004353514  max resid 0.0006624694 
    ## ... Similar to previous best
    ## Run 276 stress 9.654711e-05 
    ## ... Procrustes: rmse 0.0004098022  max resid 0.0006519693 
    ## ... Similar to previous best
    ## Run 277 stress 8.967613e-05 
    ## ... Procrustes: rmse 0.0001307906  max resid 0.0001943599 
    ## ... Similar to previous best
    ## Run 278 stress 0.0006470484 
    ## Run 279 stress 9.545656e-05 
    ## ... Procrustes: rmse 0.0004105194  max resid 0.0006561301 
    ## ... Similar to previous best
    ## Run 280 stress 9.904892e-05 
    ## ... Procrustes: rmse 0.000438315  max resid 0.000650589 
    ## ... Similar to previous best
    ## Run 281 stress 0.001031884 
    ## Run 282 stress 9.959975e-05 
    ## ... Procrustes: rmse 0.000448614  max resid 0.0006823976 
    ## ... Similar to previous best
    ## Run 283 stress 0.0006905779 
    ## Run 284 stress 9.923235e-05 
    ## ... Procrustes: rmse 0.0004259874  max resid 0.0006291659 
    ## ... Similar to previous best
    ## Run 285 stress 9.838246e-05 
    ## ... Procrustes: rmse 0.0004407602  max resid 0.0006709441 
    ## ... Similar to previous best
    ## Run 286 stress 0.0008865352 
    ## Run 287 stress 0.0007938325 
    ## Run 288 stress 9.922e-05 
    ## ... Procrustes: rmse 0.0004230367  max resid 0.000629791 
    ## ... Similar to previous best
    ## Run 289 stress 9.956361e-05 
    ## ... Procrustes: rmse 0.0004663303  max resid 0.000732891 
    ## ... Similar to previous best
    ## Run 290 stress 9.576928e-05 
    ## ... Procrustes: rmse 0.000412032  max resid 0.0006087556 
    ## ... Similar to previous best
    ## Run 291 stress 0.0007847406 
    ## Run 292 stress 9.84372e-05 
    ## ... Procrustes: rmse 0.0004233888  max resid 0.0006253907 
    ## ... Similar to previous best
    ## Run 293 stress 0.001288958 
    ## Run 294 stress 9.93823e-05 
    ## ... Procrustes: rmse 0.0004209433  max resid 0.0006701598 
    ## ... Similar to previous best
    ## Run 295 stress 9.999229e-05 
    ## ... Procrustes: rmse 0.0004922347  max resid 0.0007200924 
    ## ... Similar to previous best
    ## Run 296 stress 9.35587e-05 
    ## ... Procrustes: rmse 0.0001020677  max resid 0.0001552632 
    ## ... Similar to previous best
    ## Run 297 stress 9.925084e-05 
    ## ... Procrustes: rmse 0.000421053  max resid 0.0006700278 
    ## ... Similar to previous best
    ## Run 298 stress 0.0004918799 
    ## ... Procrustes: rmse 0.004321278  max resid 0.006682631 
    ## ... Similar to previous best
    ## Run 299 stress 9.856189e-05 
    ## ... Procrustes: rmse 0.0004378838  max resid 0.0006685156 
    ## ... Similar to previous best
    ## Run 300 stress 0.3407338 
    ## Run 301 stress 0.2711872 
    ## Run 302 stress 9.851593e-05 
    ## ... Procrustes: rmse 0.0004110321  max resid 0.0006533148 
    ## ... Similar to previous best
    ## Run 303 stress 9.461529e-05 
    ## ... Procrustes: rmse 0.0004061616  max resid 0.0006261807 
    ## ... Similar to previous best
    ## Run 304 stress 9.960449e-05 
    ## ... Procrustes: rmse 0.0004925082  max resid 0.0007229617 
    ## ... Similar to previous best
    ## Run 305 stress 9.92473e-05 
    ## ... Procrustes: rmse 0.0004644569  max resid 0.0007283071 
    ## ... Similar to previous best
    ## Run 306 stress 9.981546e-05 
    ## ... Procrustes: rmse 0.0004366203  max resid 0.0006499614 
    ## ... Similar to previous best
    ## Run 307 stress 0.0002269998 
    ## ... Procrustes: rmse 0.001999006  max resid 0.003090381 
    ## ... Similar to previous best
    ## Run 308 stress 9.80689e-05 
    ## ... Procrustes: rmse 0.0004419574  max resid 0.0006732209 
    ## ... Similar to previous best
    ## Run 309 stress 0.0006932854 
    ## Run 310 stress 9.779378e-05 
    ## ... Procrustes: rmse 0.0004204789  max resid 0.0006395125 
    ## ... Similar to previous best
    ## Run 311 stress 9.929529e-05 
    ## ... Procrustes: rmse 0.0004280204  max resid 0.0006541437 
    ## ... Similar to previous best
    ## Run 312 stress 0.0006461156 
    ## Run 313 stress 0.0008229766 
    ## Run 314 stress 0.0004240766 
    ## ... Procrustes: rmse 0.003728813  max resid 0.005765028 
    ## ... Similar to previous best
    ## Run 315 stress 9.947972e-05 
    ## ... Procrustes: rmse 0.0004655046  max resid 0.0007316629 
    ## ... Similar to previous best
    ## Run 316 stress 0.0008997435 
    ## Run 317 stress 0.0008857586 
    ## Run 318 stress 9.865039e-05 
    ## ... Procrustes: rmse 0.0004361737  max resid 0.0006517164 
    ## ... Similar to previous best
    ## Run 319 stress 0.001089737 
    ## Run 320 stress 9.784903e-05 
    ## ... Procrustes: rmse 0.0004126933  max resid 0.0006647407 
    ## ... Similar to previous best
    ## Run 321 stress 9.964033e-05 
    ## ... Procrustes: rmse 0.0004925093  max resid 0.0007225318 
    ## ... Similar to previous best
    ## Run 322 stress 9.900405e-05 
    ## ... Procrustes: rmse 0.0004259254  max resid 0.0006309504 
    ## ... Similar to previous best
    ## Run 323 stress 9.787654e-05 
    ## ... Procrustes: rmse 0.0004206531  max resid 0.0006225107 
    ## ... Similar to previous best
    ## Run 324 stress 0.3407335 
    ## Run 325 stress 9.938748e-05 
    ## ... Procrustes: rmse 0.0004254569  max resid 0.0006772386 
    ## ... Similar to previous best
    ## Run 326 stress 0.0006850795 
    ## Run 327 stress 0.0007929824 
    ## Run 328 stress 9.834703e-05 
    ## ... Procrustes: rmse 0.0004225953  max resid 0.0006242863 
    ## ... Similar to previous best
    ## Run 329 stress 9.992282e-05 
    ## ... Procrustes: rmse 0.0004243779  max resid 0.0006709493 
    ## ... Similar to previous best
    ## Run 330 stress 9.961393e-05 
    ## ... Procrustes: rmse 0.0004468215  max resid 0.000680306 
    ## ... Similar to previous best
    ## Run 331 stress 9.648719e-05 
    ## ... Procrustes: rmse 0.0004068008  max resid 0.0006266307 
    ## ... Similar to previous best
    ## Run 332 stress 9.874507e-05 
    ## ... Procrustes: rmse 0.0004380651  max resid 0.0006536655 
    ## ... Similar to previous best
    ## Run 333 stress 9.877339e-05 
    ## ... Procrustes: rmse 0.0004248378  max resid 0.0006797028 
    ## ... Similar to previous best
    ## Run 334 stress 9.286382e-05 
    ## ... Procrustes: rmse 0.0003928655  max resid 0.0006249241 
    ## ... Similar to previous best
    ## Run 335 stress 0.001159483 
    ## Run 336 stress 9.931483e-05 
    ## ... Procrustes: rmse 0.0004644347  max resid 0.0007299726 
    ## ... Similar to previous best
    ## Run 337 stress 9.963309e-05 
    ## ... Procrustes: rmse 0.0004925537  max resid 0.0007226393 
    ## ... Similar to previous best
    ## Run 338 stress 9.963501e-05 
    ## ... Procrustes: rmse 0.0004051996  max resid 0.0005898044 
    ## ... Similar to previous best
    ## Run 339 stress 0.0006687146 
    ## Run 340 stress 9.984134e-05 
    ## ... Procrustes: rmse 0.0004899266  max resid 0.0007365591 
    ## ... Similar to previous best
    ## Run 341 stress 0.0007735283 
    ## Run 342 stress 9.937586e-05 
    ## ... Procrustes: rmse 0.0004245941  max resid 0.0006314264 
    ## ... Similar to previous best
    ## Run 343 stress 0.0008067227 
    ## Run 344 stress 9.730446e-05 
    ## ... Procrustes: rmse 0.0004180358  max resid 0.0006202627 
    ## ... Similar to previous best
    ## Run 345 stress 0.0005081332 
    ## ... Procrustes: rmse 0.004462807  max resid 0.006911983 
    ## ... Similar to previous best
    ## Run 346 stress 9.957507e-05 
    ## ... Procrustes: rmse 0.000427156  max resid 0.0006861799 
    ## ... Similar to previous best
    ## Run 347 stress 9.877522e-05 
    ## ... Procrustes: rmse 0.0004241666  max resid 0.0006272423 
    ## ... Similar to previous best
    ## Run 348 stress 9.937504e-05 
    ## ... Procrustes: rmse 0.0004289651  max resid 0.0006547956 
    ## ... Similar to previous best
    ## Run 349 stress 9.747837e-05 
    ## ... Procrustes: rmse 0.0004193824  max resid 0.0006212155 
    ## ... Similar to previous best
    ## Run 350 stress 0.0006430016 
    ## Run 351 stress 9.902181e-05 
    ## ... Procrustes: rmse 0.0004896994  max resid 0.0007185775 
    ## ... Similar to previous best
    ## Run 352 stress 9.946324e-05 
    ## ... Procrustes: rmse 0.000435826  max resid 0.000662621 
    ## ... Similar to previous best
    ## Run 353 stress 0.0005125201 
    ## ... Procrustes: rmse 0.004489726  max resid 0.006962434 
    ## ... Similar to previous best
    ## Run 354 stress 9.900832e-05 
    ## ... Procrustes: rmse 0.0004244409  max resid 0.0006287957 
    ## ... Similar to previous best
    ## Run 355 stress 9.58878e-05 
    ## ... Procrustes: rmse 0.0004099462  max resid 0.0006548836 
    ## ... Similar to previous best
    ## Run 356 stress 9.871815e-05 
    ## ... Procrustes: rmse 0.0004243992  max resid 0.0006444722 
    ## ... Similar to previous best
    ## Run 357 stress 0.0005942287 
    ## Run 358 stress 9.861227e-05 
    ## ... Procrustes: rmse 0.0004381322  max resid 0.0006544296 
    ## ... Similar to previous best
    ## Run 359 stress 0.0009491706 
    ## Run 360 stress 0.3407343 
    ## Run 361 stress 0.0006385323 
    ## Run 362 stress 9.976479e-05 
    ## ... Procrustes: rmse 0.0004190408  max resid 0.0006620596 
    ## ... Similar to previous best
    ## Run 363 stress 9.601048e-05 
    ## ... Procrustes: rmse 0.000425345  max resid 0.0006358925 
    ## ... Similar to previous best
    ## Run 364 stress 9.877745e-05 
    ## ... Procrustes: rmse 0.0004190548  max resid 0.0006665637 
    ## ... Similar to previous best
    ## Run 365 stress 9.8497e-05 
    ## ... Procrustes: rmse 0.0004842968  max resid 0.0007120228 
    ## ... Similar to previous best
    ## Run 366 stress 9.867124e-05 
    ## ... Procrustes: rmse 0.0004440469  max resid 0.0006762622 
    ## ... Similar to previous best
    ## Run 367 stress 9.739041e-05 
    ## ... Procrustes: rmse 0.0004194994  max resid 0.0006730052 
    ## ... Similar to previous best
    ## Run 368 stress 9.99945e-05 
    ## ... Procrustes: rmse 0.0004495844  max resid 0.0006843929 
    ## ... Similar to previous best
    ## Run 369 stress 9.953052e-05 
    ## ... Procrustes: rmse 0.0004661901  max resid 0.0007329721 
    ## ... Similar to previous best
    ## Run 370 stress 9.839902e-05 
    ## ... Procrustes: rmse 0.000417604  max resid 0.0006641659 
    ## ... Similar to previous best
    ## Run 371 stress 9.86449e-05 
    ## ... Procrustes: rmse 0.0004241514  max resid 0.0006264459 
    ## ... Similar to previous best
    ## Run 372 stress 9.835144e-05 
    ## ... Procrustes: rmse 0.0004240532  max resid 0.0006483296 
    ## ... Similar to previous best
    ## Run 373 stress 9.266717e-05 
    ## ... Procrustes: rmse 0.0003940314  max resid 0.0005997486 
    ## ... Similar to previous best
    ## Run 374 stress 9.864919e-05 
    ## ... Procrustes: rmse 0.0004254109  max resid 0.000646892 
    ## ... Similar to previous best
    ## Run 375 stress 9.79112e-05 
    ## ... Procrustes: rmse 0.000412911  max resid 0.0006565308 
    ## ... Similar to previous best
    ## Run 376 stress 9.69338e-05 
    ## ... Procrustes: rmse 0.0004794069  max resid 0.0007045024 
    ## ... Similar to previous best
    ## Run 377 stress 9.574903e-05 
    ## ... Procrustes: rmse 0.0004132782  max resid 0.0006323127 
    ## ... Similar to previous best
    ## Run 378 stress 0.0005604971 
    ## ... Procrustes: rmse 0.004916113  max resid 0.00759928 
    ## ... Similar to previous best
    ## Run 379 stress 9.923394e-05 
    ## ... Procrustes: rmse 0.000414818  max resid 0.000662417 
    ## ... Similar to previous best
    ## Run 380 stress 9.61571e-05 
    ## ... Procrustes: rmse 0.000474947  max resid 0.0006972648 
    ## ... Similar to previous best
    ## Run 381 stress 9.990353e-05 
    ## ... Procrustes: rmse 0.0004296202  max resid 0.0006347217 
    ## ... Similar to previous best
    ## Run 382 stress 9.697011e-05 
    ## ... Procrustes: rmse 0.0004304379  max resid 0.0006416522 
    ## ... Similar to previous best
    ## Run 383 stress 9.756226e-05 
    ## ... Procrustes: rmse 0.0004271392  max resid 0.0006553181 
    ## ... Similar to previous best
    ## Run 384 stress 9.456669e-05 
    ## ... Procrustes: rmse 0.0004070962  max resid 0.0006045632 
    ## ... Similar to previous best
    ## Run 385 stress 9.847082e-05 
    ## ... Procrustes: rmse 0.0004169365  max resid 0.0006655037 
    ## ... Similar to previous best
    ## Run 386 stress 9.943968e-05 
    ## ... Procrustes: rmse 0.0004216255  max resid 0.0006648708 
    ## ... Similar to previous best
    ## Run 387 stress 0.0006730592 
    ## Run 388 stress 9.967756e-05 
    ## ... Procrustes: rmse 0.0004484053  max resid 0.0006831126 
    ## ... Similar to previous best
    ## Run 389 stress 9.712647e-05 
    ## ... Procrustes: rmse 0.0004180449  max resid 0.0006684651 
    ## ... Similar to previous best
    ## Run 390 stress 0.0007291634 
    ## Run 391 stress 0.0005934524 
    ## Run 392 stress 9.981911e-05 
    ## ... Procrustes: rmse 0.0004636453  max resid 0.0007116322 
    ## ... Similar to previous best
    ## Run 393 stress 9.755292e-05 
    ## ... Procrustes: rmse 0.0004101926  max resid 0.0006564387 
    ## ... Similar to previous best
    ## Run 394 stress 0.0006394944 
    ## Run 395 stress 9.936978e-05 
    ## ... Procrustes: rmse 0.0004914866  max resid 0.0007218607 
    ## ... Similar to previous best
    ## Run 396 stress 9.849696e-05 
    ## ... Procrustes: rmse 0.0004255977  max resid 0.0006481989 
    ## ... Similar to previous best
    ## Run 397 stress 9.980452e-05 
    ## ... Procrustes: rmse 0.0004295141  max resid 0.0006879764 
    ## ... Similar to previous best
    ## Run 398 stress 0.0005224216 
    ## ... Procrustes: rmse 0.004589448  max resid 0.007105638 
    ## ... Similar to previous best
    ## Run 399 stress 9.698295e-05 
    ## ... Procrustes: rmse 0.0004090962  max resid 0.0006494947 
    ## ... Similar to previous best
    ## Run 400 stress 0.001222913 
    ## Run 401 stress 9.780829e-05 
    ## ... Procrustes: rmse 0.0003159858  max resid 0.0005554929 
    ## ... Similar to previous best
    ## Run 402 stress 0.0005117356 
    ## ... Procrustes: rmse 0.004495727  max resid 0.006954255 
    ## ... Similar to previous best
    ## Run 403 stress 9.430555e-05 
    ## ... Procrustes: rmse 0.0003998699  max resid 0.0006367376 
    ## ... Similar to previous best
    ## Run 404 stress 9.761042e-05 
    ## ... Procrustes: rmse 0.0004185435  max resid 0.0006755084 
    ## ... Similar to previous best
    ## Run 405 stress 0.0006637917 
    ## Run 406 stress 0.0004104339 
    ## ... Procrustes: rmse 0.003622179  max resid 0.005597265 
    ## ... Similar to previous best
    ## Run 407 stress 9.966439e-05 
    ## ... Procrustes: rmse 0.0004594153  max resid 0.0007035392 
    ## ... Similar to previous best
    ## Run 408 stress 0.001013225 
    ## Run 409 stress 9.753534e-05 
    ## ... Procrustes: rmse 0.0004342097  max resid 0.0006441246 
    ## ... Similar to previous best
    ## Run 410 stress 9.840691e-05 
    ## ... Procrustes: rmse 0.0004432921  max resid 0.0006753715 
    ## ... Similar to previous best
    ## Run 411 stress 9.995714e-05 
    ## ... Procrustes: rmse 0.0004681446  max resid 0.0007354338 
    ## ... Similar to previous best
    ## Run 412 stress 0.001453294 
    ## Run 413 stress 0.001120418 
    ## Run 414 stress 0.001017588 
    ## Run 415 stress 0.001140978 
    ## Run 416 stress 9.965326e-05 
    ## ... Procrustes: rmse 0.000426979  max resid 0.0006469387 
    ## ... Similar to previous best
    ## Run 417 stress 0.0002700853 
    ## ... Procrustes: rmse 0.002378268  max resid 0.003678965 
    ## ... Similar to previous best
    ## Run 418 stress 9.923972e-05 
    ## ... Procrustes: rmse 0.0004261078  max resid 0.000628789 
    ## ... Similar to previous best
    ## Run 419 stress 9.814454e-05 
    ## ... Procrustes: rmse 0.0004236775  max resid 0.0006446171 
    ## ... Similar to previous best
    ## Run 420 stress 0.0007350944 
    ## Run 421 stress 0.0006858117 
    ## Run 422 stress 9.989415e-05 
    ## ... Procrustes: rmse 0.000416611  max resid 0.0006653694 
    ## ... Similar to previous best
    ## Run 423 stress 0.0005377225 
    ## ... Procrustes: rmse 0.004716391  max resid 0.00730972 
    ## ... Similar to previous best
    ## Run 424 stress 9.9231e-05 
    ## ... Procrustes: rmse 0.000420001  max resid 0.0006671733 
    ## ... Similar to previous best
    ## Run 425 stress 9.919817e-05 
    ## ... Procrustes: rmse 0.0004869306  max resid 0.0007128391 
    ## ... Similar to previous best
    ## Run 426 stress 9.898277e-05 
    ## ... Procrustes: rmse 0.0004808181  max resid 0.0007320434 
    ## ... Similar to previous best
    ## Run 427 stress 0.0009039244 
    ## Run 428 stress 9.90178e-05 
    ## ... Procrustes: rmse 0.0004207078  max resid 0.000626944 
    ## ... Similar to previous best
    ## Run 429 stress 9.825895e-05 
    ## ... Procrustes: rmse 0.0004229271  max resid 0.0006767745 
    ## ... Similar to previous best
    ## Run 430 stress 9.897866e-05 
    ## ... Procrustes: rmse 0.0004651205  max resid 0.0007126922 
    ## ... Similar to previous best
    ## Run 431 stress 9.647927e-05 
    ## ... Procrustes: rmse 0.0004154536  max resid 0.0006644998 
    ## ... Similar to previous best
    ## Run 432 stress 9.704892e-05 
    ## ... Procrustes: rmse 0.0004083161  max resid 0.0006104397 
    ## ... Similar to previous best
    ## Run 433 stress 9.850697e-05 
    ## ... Procrustes: rmse 0.0004238557  max resid 0.0006274887 
    ## ... Similar to previous best
    ## Run 434 stress 9.5286e-05 
    ## ... Procrustes: rmse 0.0004025778  max resid 0.0006428523 
    ## ... Similar to previous best
    ## Run 435 stress 9.453342e-05 
    ## ... Procrustes: rmse 0.0003969722  max resid 0.0005831058 
    ## ... Similar to previous best
    ## Run 436 stress 9.910055e-05 
    ## ... Procrustes: rmse 0.000425691  max resid 0.0006847296 
    ## ... Similar to previous best
    ## Run 437 stress 0.257069 
    ## Run 438 stress 8.889828e-05 
    ## ... Procrustes: rmse 0.0001036409  max resid 0.0001713528 
    ## ... Similar to previous best
    ## Run 439 stress 9.657117e-05 
    ## ... Procrustes: rmse 0.0004153899  max resid 0.0006157788 
    ## ... Similar to previous best
    ## Run 440 stress 0.3407336 
    ## Run 441 stress 9.797007e-05 
    ## ... Procrustes: rmse 0.0004828239  max resid 0.0007288139 
    ## ... Similar to previous best
    ## Run 442 stress 0.001128964 
    ## Run 443 stress 8.877689e-05 
    ## ... Procrustes: rmse 9.403579e-05  max resid 0.000151532 
    ## ... Similar to previous best
    ## Run 444 stress 0.001260311 
    ## Run 445 stress 9.906872e-05 
    ## ... Procrustes: rmse 0.0004254934  max resid 0.0006277407 
    ## ... Similar to previous best
    ## Run 446 stress 0.001076847 
    ## Run 447 stress 9.997814e-05 
    ## ... Procrustes: rmse 0.0004236609  max resid 0.000672894 
    ## ... Similar to previous best
    ## Run 448 stress 0.00102513 
    ## Run 449 stress 9.263448e-05 
    ## ... Procrustes: rmse 0.0004121503  max resid 0.0006180888 
    ## ... Similar to previous best
    ## Run 450 stress 9.934746e-05 
    ## ... Procrustes: rmse 0.0004183119  max resid 0.0006597069 
    ## ... Similar to previous best
    ## Run 451 stress 9.966262e-05 
    ## ... Procrustes: rmse 0.0004928552  max resid 0.0007231447 
    ## ... Similar to previous best
    ## Run 452 stress 9.762291e-05 
    ## ... Procrustes: rmse 0.0004088406  max resid 0.0006525088 
    ## ... Similar to previous best
    ## Run 453 stress 9.978528e-05 
    ## ... Procrustes: rmse 0.000493319  max resid 0.0007237939 
    ## ... Similar to previous best
    ## Run 454 stress 9.930606e-05 
    ## ... Procrustes: rmse 0.0003781728  max resid 0.0006064651 
    ## ... Similar to previous best
    ## Run 455 stress 9.968496e-05 
    ## ... Procrustes: rmse 0.0004625861  max resid 0.0007304549 
    ## ... Similar to previous best
    ## Run 456 stress 9.799753e-05 
    ## ... Procrustes: rmse 0.0004223889  max resid 0.0006486213 
    ## ... Similar to previous best
    ## Run 457 stress 9.808178e-05 
    ## ... Procrustes: rmse 0.0004216265  max resid 0.0006746348 
    ## ... Similar to previous best
    ## Run 458 stress 9.529286e-05 
    ## ... Procrustes: rmse 0.0004022462  max resid 0.0006421416 
    ## ... Similar to previous best
    ## Run 459 stress 9.973932e-05 
    ## ... Procrustes: rmse 0.0004407357  max resid 0.000671304 
    ## ... Similar to previous best
    ## Run 460 stress 0.0008185809 
    ## Run 461 stress 0.0008591836 
    ## Run 462 stress 0.0006465137 
    ## Run 463 stress 0.2467071 
    ## Run 464 stress 0.0007970348 
    ## Run 465 stress 0.0009556731 
    ## Run 466 stress 0.3192024 
    ## Run 467 stress 9.895975e-05 
    ## ... Procrustes: rmse 0.0004615672  max resid 0.0007258861 
    ## ... Similar to previous best
    ## Run 468 stress 0.0008025419 
    ## Run 469 stress 9.891949e-05 
    ## ... Procrustes: rmse 0.0004403528  max resid 0.000655091 
    ## ... Similar to previous best
    ## Run 470 stress 0.001097662 
    ## Run 471 stress 9.190695e-05 
    ## ... Procrustes: rmse 0.0003889302  max resid 0.0006274819 
    ## ... Similar to previous best
    ## Run 472 stress 9.951804e-05 
    ## ... Procrustes: rmse 0.0004655624  max resid 0.0007326494 
    ## ... Similar to previous best
    ## Run 473 stress 9.902463e-05 
    ## ... Procrustes: rmse 0.0004434066  max resid 0.0006758302 
    ## ... Similar to previous best
    ## Run 474 stress 9.783222e-05 
    ## ... Procrustes: rmse 0.0004823681  max resid 0.0007276367 
    ## ... Similar to previous best
    ## Run 475 stress 9.890308e-05 
    ## ... Procrustes: rmse 0.0004372819  max resid 0.0006514312 
    ## ... Similar to previous best
    ## Run 476 stress 9.899903e-05 
    ## ... Procrustes: rmse 0.0004248302  max resid 0.0006772117 
    ## ... Similar to previous best
    ## Run 477 stress 9.928514e-05 
    ## ... Procrustes: rmse 0.0003039268  max resid 0.0005185012 
    ## ... Similar to previous best
    ## Run 478 stress 9.819432e-05 
    ## ... Procrustes: rmse 0.0004148053  max resid 0.0006592727 
    ## ... Similar to previous best
    ## Run 479 stress 9.849976e-05 
    ## ... Procrustes: rmse 0.0004178108  max resid 0.0006652244 
    ## ... Similar to previous best
    ## Run 480 stress 9.747694e-05 
    ## ... Procrustes: rmse 0.0004193325  max resid 0.0006739633 
    ## ... Similar to previous best
    ## Run 481 stress 0.0007730564 
    ## Run 482 stress 9.91387e-05 
    ## ... Procrustes: rmse 0.0004266344  max resid 0.0006829824 
    ## ... Similar to previous best
    ## Run 483 stress 9.820805e-05 
    ## ... Procrustes: rmse 0.0004204412  max resid 0.0006415869 
    ## ... Similar to previous best
    ## Run 484 stress 9.698674e-05 
    ## ... Procrustes: rmse 0.0004129708  max resid 0.0006121534 
    ## ... Similar to previous best
    ## Run 485 stress 9.869706e-05 
    ## ... Procrustes: rmse 0.0004176553  max resid 0.0006637271 
    ## ... Similar to previous best
    ## Run 486 stress 9.95492e-05 
    ## ... Procrustes: rmse 0.0004144373  max resid 0.0006619833 
    ## ... Similar to previous best
    ## Run 487 stress 0.0006470739 
    ## Run 488 stress 9.908247e-05 
    ## ... Procrustes: rmse 0.000427523  max resid 0.0006511676 
    ## ... Similar to previous best
    ## Run 489 stress 0.0001329014 
    ## ... Procrustes: rmse 0.001167845  max resid 0.001808616 
    ## ... Similar to previous best
    ## Run 490 stress 0.001343603 
    ## Run 491 stress 0.0009775291 
    ## Run 492 stress 9.809515e-05 
    ## ... Procrustes: rmse 0.0004198895  max resid 0.0006724158 
    ## ... Similar to previous best
    ## Run 493 stress 9.907947e-05 
    ## ... Procrustes: rmse 0.0004263885  max resid 0.0006860499 
    ## ... Similar to previous best
    ## Run 494 stress 9.990213e-05 
    ## ... Procrustes: rmse 0.0004886828  max resid 0.000720139 
    ## ... Similar to previous best
    ## Run 495 stress 0.001281701 
    ## Run 496 stress 0.0009832615 
    ## Run 497 stress 0.001042749 
    ## Run 498 stress 9.752629e-05 
    ## ... Procrustes: rmse 0.0004056767  max resid 0.0006140011 
    ## ... Similar to previous best
    ## Run 499 stress 0.00111759 
    ## Run 500 stress 9.988914e-05 
    ## ... Procrustes: rmse 0.000867062  max resid 0.001353675 
    ## ... Similar to previous best
    ## *** Best solution repeated 279 times

    ## Warning in metaMDS(PD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
PD_beta_env_M_NMDS <- metaMDS(PD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.019724e-05 
    ## Run 1 stress 9.471983e-05 
    ## ... Procrustes: rmse 0.00025125  max resid 0.0003890997 
    ## ... Similar to previous best
    ## Run 2 stress 9.555328e-05 
    ## ... Procrustes: rmse 1.525434e-05  max resid 1.984421e-05 
    ## ... Similar to previous best
    ## Run 3 stress 9.189121e-05 
    ## ... Procrustes: rmse 0.0001926966  max resid 0.0003584127 
    ## ... Similar to previous best
    ## Run 4 stress 8.680582e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.00015379  max resid 0.0002709072 
    ## ... Similar to previous best
    ## Run 5 stress 9.672093e-05 
    ## ... Procrustes: rmse 0.0001615013  max resid 0.0002871716 
    ## ... Similar to previous best
    ## Run 6 stress 9.61254e-05 
    ## ... Procrustes: rmse 0.0001025013  max resid 0.0001650299 
    ## ... Similar to previous best
    ## Run 7 stress 0.3083098 
    ## Run 8 stress 9.570114e-05 
    ## ... Procrustes: rmse 0.00021701  max resid 0.0003373005 
    ## ... Similar to previous best
    ## Run 9 stress 9.508425e-05 
    ## ... Procrustes: rmse 0.0002020014  max resid 0.0004079274 
    ## ... Similar to previous best
    ## Run 10 stress 9.351355e-05 
    ## ... Procrustes: rmse 0.0001599641  max resid 0.0002855134 
    ## ... Similar to previous best
    ## Run 11 stress 8.98445e-05 
    ## ... Procrustes: rmse 1.270498e-05  max resid 1.832591e-05 
    ## ... Similar to previous best
    ## Run 12 stress 9.483604e-05 
    ## ... Procrustes: rmse 0.0001521226  max resid 0.0003102182 
    ## ... Similar to previous best
    ## Run 13 stress 9.844867e-05 
    ## ... Procrustes: rmse 0.0001652485  max resid 0.000336591 
    ## ... Similar to previous best
    ## Run 14 stress 9.065852e-05 
    ## ... Procrustes: rmse 0.0001535507  max resid 0.0002716456 
    ## ... Similar to previous best
    ## Run 15 stress 9.508612e-05 
    ## ... Procrustes: rmse 0.000163031  max resid 0.000292034 
    ## ... Similar to previous best
    ## Run 16 stress 9.719002e-05 
    ## ... Procrustes: rmse 0.0001677797  max resid 0.0003001786 
    ## ... Similar to previous best
    ## Run 17 stress 9.367833e-05 
    ## ... Procrustes: rmse 0.0001573041  max resid 0.0002772657 
    ## ... Similar to previous best
    ## Run 18 stress 9.63494e-05 
    ## ... Procrustes: rmse 0.0001647689  max resid 0.0002929873 
    ## ... Similar to previous best
    ## Run 19 stress 9.471501e-05 
    ## ... Procrustes: rmse 1.614177e-05  max resid 2.839873e-05 
    ## ... Similar to previous best
    ## Run 20 stress 8.760169e-05 
    ## ... Procrustes: rmse 0.0001924849  max resid 0.0003845442 
    ## ... Similar to previous best
    ## Run 21 stress 9.839706e-05 
    ## ... Procrustes: rmse 2.811275e-05  max resid 4.062371e-05 
    ## ... Similar to previous best
    ## Run 22 stress 9.256564e-05 
    ## ... Procrustes: rmse 1.164013e-05  max resid 2.009411e-05 
    ## ... Similar to previous best
    ## Run 23 stress 9.929509e-05 
    ## ... Procrustes: rmse 0.0002143965  max resid 0.0003339865 
    ## ... Similar to previous best
    ## Run 24 stress 9.918564e-05 
    ## ... Procrustes: rmse 0.0002256974  max resid 0.0003457763 
    ## ... Similar to previous best
    ## Run 25 stress 9.665201e-05 
    ## ... Procrustes: rmse 0.0002228851  max resid 0.0003429309 
    ## ... Similar to previous best
    ## Run 26 stress 6.614467e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001686458  max resid 0.0002933697 
    ## ... Similar to previous best
    ## Run 27 stress 9.460154e-05 
    ## ... Procrustes: rmse 0.0001825507  max resid 0.0003203531 
    ## ... Similar to previous best
    ## Run 28 stress 9.226638e-05 
    ## ... Procrustes: rmse 0.0002391557  max resid 0.0005461874 
    ## ... Similar to previous best
    ## Run 29 stress 9.669226e-05 
    ## ... Procrustes: rmse 0.0002499333  max resid 0.0005666918 
    ## ... Similar to previous best
    ## Run 30 stress 9.699511e-05 
    ## ... Procrustes: rmse 0.000249283  max resid 0.0005627644 
    ## ... Similar to previous best
    ## Run 31 stress 9.313351e-05 
    ## ... Procrustes: rmse 0.000240841  max resid 0.0005484656 
    ## ... Similar to previous best
    ## Run 32 stress 0.2854308 
    ## Run 33 stress 0.3069398 
    ## Run 34 stress 8.663192e-05 
    ## ... Procrustes: rmse 0.0001193906  max resid 0.0002087797 
    ## ... Similar to previous best
    ## Run 35 stress 7.453168e-05 
    ## ... Procrustes: rmse 0.0001095558  max resid 0.0002051907 
    ## ... Similar to previous best
    ## Run 36 stress 9.221192e-05 
    ## ... Procrustes: rmse 0.0001178591  max resid 0.0002320994 
    ## ... Similar to previous best
    ## Run 37 stress 9.855191e-05 
    ## ... Procrustes: rmse 0.000245343  max resid 0.0005608306 
    ## ... Similar to previous best
    ## Run 38 stress 9.398737e-05 
    ## ... Procrustes: rmse 0.0002082826  max resid 0.0003655901 
    ## ... Similar to previous best
    ## Run 39 stress 9.840553e-05 
    ## ... Procrustes: rmse 0.000253412  max resid 0.00057134 
    ## ... Similar to previous best
    ## Run 40 stress 9.350063e-05 
    ## ... Procrustes: rmse 0.0001555361  max resid 0.0002222288 
    ## ... Similar to previous best
    ## Run 41 stress 9.932523e-05 
    ## ... Procrustes: rmse 0.0002463325  max resid 0.0005631501 
    ## ... Similar to previous best
    ## Run 42 stress 9.999242e-05 
    ## ... Procrustes: rmse 0.0002563609  max resid 0.0005762984 
    ## ... Similar to previous best
    ## Run 43 stress 9.779813e-05 
    ## ... Procrustes: rmse 0.0002179989  max resid 0.0003943175 
    ## ... Similar to previous best
    ## Run 44 stress 9.747894e-05 
    ## ... Procrustes: rmse 0.0002187195  max resid 0.0003841217 
    ## ... Similar to previous best
    ## Run 45 stress 9.684318e-05 
    ## ... Procrustes: rmse 0.0002223303  max resid 0.0004882833 
    ## ... Similar to previous best
    ## Run 46 stress 8.299141e-05 
    ## ... Procrustes: rmse 0.0001656101  max resid 0.0002815142 
    ## ... Similar to previous best
    ## Run 47 stress 6.999892e-05 
    ## ... Procrustes: rmse 0.0001237694  max resid 0.0001945864 
    ## ... Similar to previous best
    ## Run 48 stress 8.688409e-05 
    ## ... Procrustes: rmse 3.505195e-05  max resid 5.749281e-05 
    ## ... Similar to previous best
    ## Run 49 stress 8.463937e-05 
    ## ... Procrustes: rmse 0.0001635786  max resid 0.0002800568 
    ## ... Similar to previous best
    ## Run 50 stress 8.718807e-05 
    ## ... Procrustes: rmse 0.0001828297  max resid 0.0003418188 
    ## ... Similar to previous best
    ## Run 51 stress 9.965839e-05 
    ## ... Procrustes: rmse 4.991009e-05  max resid 8.928192e-05 
    ## ... Similar to previous best
    ## Run 52 stress 9.708858e-05 
    ## ... Procrustes: rmse 3.908854e-05  max resid 5.921048e-05 
    ## ... Similar to previous best
    ## Run 53 stress 9.89023e-05 
    ## ... Procrustes: rmse 7.236549e-05  max resid 0.0001481328 
    ## ... Similar to previous best
    ## Run 54 stress 9.179179e-05 
    ## ... Procrustes: rmse 0.0001805387  max resid 0.000315185 
    ## ... Similar to previous best
    ## Run 55 stress 8.23634e-05 
    ## ... Procrustes: rmse 0.0001617103  max resid 0.0003043954 
    ## ... Similar to previous best
    ## Run 56 stress 9.001441e-05 
    ## ... Procrustes: rmse 0.0001183214  max resid 0.0002055111 
    ## ... Similar to previous best
    ## Run 57 stress 0.2848523 
    ## Run 58 stress 0.2987885 
    ## Run 59 stress 9.946601e-05 
    ## ... Procrustes: rmse 0.000212425  max resid 0.0003759623 
    ## ... Similar to previous best
    ## Run 60 stress 9.355817e-05 
    ## ... Procrustes: rmse 0.0001794862  max resid 0.000316513 
    ## ... Similar to previous best
    ## Run 61 stress 8.473844e-05 
    ## ... Procrustes: rmse 0.0001678025  max resid 0.0002988729 
    ## ... Similar to previous best
    ## Run 62 stress 9.688322e-05 
    ## ... Procrustes: rmse 0.0001634419  max resid 0.0003559298 
    ## ... Similar to previous best
    ## Run 63 stress 8.885868e-05 
    ## ... Procrustes: rmse 0.0001208085  max resid 0.0002169668 
    ## ... Similar to previous best
    ## Run 64 stress 9.009784e-05 
    ## ... Procrustes: rmse 0.0001536238  max resid 0.000331269 
    ## ... Similar to previous best
    ## Run 65 stress 9.243035e-05 
    ## ... Procrustes: rmse 0.00017825  max resid 0.000307417 
    ## ... Similar to previous best
    ## Run 66 stress 9.21972e-05 
    ## ... Procrustes: rmse 0.0001388345  max resid 0.0003130726 
    ## ... Similar to previous best
    ## Run 67 stress 5.746264e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 9.036349e-05  max resid 0.00013923 
    ## ... Similar to previous best
    ## Run 68 stress 9.478604e-05 
    ## ... Procrustes: rmse 0.0002248177  max resid 0.0004232852 
    ## ... Similar to previous best
    ## Run 69 stress 8.707407e-05 
    ## ... Procrustes: rmse 0.0001058633  max resid 0.0001450679 
    ## ... Similar to previous best
    ## Run 70 stress 9.961753e-05 
    ## ... Procrustes: rmse 0.0001798254  max resid 0.0002500629 
    ## ... Similar to previous best
    ## Run 71 stress 9.852223e-05 
    ## ... Procrustes: rmse 0.0002268205  max resid 0.0004458996 
    ## ... Similar to previous best
    ## Run 72 stress 0.2854383 
    ## Run 73 stress 9.013495e-05 
    ## ... Procrustes: rmse 0.000133919  max resid 0.0002515856 
    ## ... Similar to previous best
    ## Run 74 stress 9.710498e-05 
    ## ... Procrustes: rmse 0.0002296833  max resid 0.0004559043 
    ## ... Similar to previous best
    ## Run 75 stress 9.689965e-05 
    ## ... Procrustes: rmse 0.0001509286  max resid 0.0002889027 
    ## ... Similar to previous best
    ## Run 76 stress 8.565992e-05 
    ## ... Procrustes: rmse 0.0002017171  max resid 0.0004117676 
    ## ... Similar to previous best
    ## Run 77 stress 9.123432e-05 
    ## ... Procrustes: rmse 0.0002194849  max resid 0.0004420178 
    ## ... Similar to previous best
    ## Run 78 stress 0.2852173 
    ## Run 79 stress 9.212499e-05 
    ## ... Procrustes: rmse 0.0001646756  max resid 0.0002344255 
    ## ... Similar to previous best
    ## Run 80 stress 9.370836e-05 
    ## ... Procrustes: rmse 0.0002225368  max resid 0.0004490338 
    ## ... Similar to previous best
    ## Run 81 stress 9.770205e-05 
    ## ... Procrustes: rmse 0.0002352021  max resid 0.0004675981 
    ## ... Similar to previous best
    ## Run 82 stress 9.882224e-05 
    ## ... Procrustes: rmse 0.0002282147  max resid 0.0004625853 
    ## ... Similar to previous best
    ## Run 83 stress 9.76505e-05 
    ## ... Procrustes: rmse 0.0002355025  max resid 0.0004708221 
    ## ... Similar to previous best
    ## Run 84 stress 9.355249e-05 
    ## ... Procrustes: rmse 0.0001220171  max resid 0.0001717685 
    ## ... Similar to previous best
    ## Run 85 stress 9.234603e-05 
    ## ... Procrustes: rmse 0.000231359  max resid 0.0004532074 
    ## ... Similar to previous best
    ## Run 86 stress 8.899791e-05 
    ## ... Procrustes: rmse 0.0001109554  max resid 0.0001583768 
    ## ... Similar to previous best
    ## Run 87 stress 8.991896e-05 
    ## ... Procrustes: rmse 0.0002160713  max resid 0.0004353958 
    ## ... Similar to previous best
    ## Run 88 stress 9.241804e-05 
    ## ... Procrustes: rmse 0.0002297702  max resid 0.0004403137 
    ## ... Similar to previous best
    ## Run 89 stress 9.918682e-05 
    ## ... Procrustes: rmse 0.000179789  max resid 0.0002502695 
    ## ... Similar to previous best
    ## Run 90 stress 9.721578e-05 
    ## ... Procrustes: rmse 0.0001760487  max resid 0.0002447989 
    ## ... Similar to previous best
    ## Run 91 stress 9.379685e-05 
    ## ... Procrustes: rmse 0.0001439417  max resid 0.0002719035 
    ## ... Similar to previous best
    ## Run 92 stress 9.026374e-05 
    ## ... Procrustes: rmse 0.0001608073  max resid 0.0002299161 
    ## ... Similar to previous best
    ## Run 93 stress 9.420307e-05 
    ## ... Procrustes: rmse 0.0002191534  max resid 0.0004418128 
    ## ... Similar to previous best
    ## Run 94 stress 9.555284e-05 
    ## ... Procrustes: rmse 0.0002390731  max resid 0.0004550249 
    ## ... Similar to previous best
    ## Run 95 stress 9.974104e-05 
    ## ... Procrustes: rmse 0.000180787  max resid 0.0002517003 
    ## ... Similar to previous best
    ## Run 96 stress 9.181284e-05 
    ## ... Procrustes: rmse 0.000165876  max resid 0.0002300014 
    ## ... Similar to previous best
    ## Run 97 stress 9.827131e-05 
    ## ... Procrustes: rmse 0.0001782224  max resid 0.0002480253 
    ## ... Similar to previous best
    ## Run 98 stress 8.993549e-05 
    ## ... Procrustes: rmse 0.000136169  max resid 0.0002569674 
    ## ... Similar to previous best
    ## Run 99 stress 9.649788e-05 
    ## ... Procrustes: rmse 0.0001748532  max resid 0.0002486088 
    ## ... Similar to previous best
    ## Run 100 stress 8.500751e-05 
    ## ... Procrustes: rmse 0.0002020851  max resid 0.0004054272 
    ## ... Similar to previous best
    ## Run 101 stress 9.598105e-05 
    ## ... Procrustes: rmse 0.0002300903  max resid 0.0004576036 
    ## ... Similar to previous best
    ## Run 102 stress 9.875834e-05 
    ## ... Procrustes: rmse 0.0001197042  max resid 0.0001876873 
    ## ... Similar to previous best
    ## Run 103 stress 9.645644e-05 
    ## ... Procrustes: rmse 0.0002313722  max resid 0.0004633097 
    ## ... Similar to previous best
    ## Run 104 stress 8.890745e-05 
    ## ... Procrustes: rmse 0.0001683663  max resid 0.0003273673 
    ## ... Similar to previous best
    ## Run 105 stress 9.550911e-05 
    ## ... Procrustes: rmse 0.0002293204  max resid 0.0004593518 
    ## ... Similar to previous best
    ## Run 106 stress 9.854571e-05 
    ## ... Procrustes: rmse 0.0002339872  max resid 0.0004640088 
    ## ... Similar to previous best
    ## Run 107 stress 9.413428e-05 
    ## ... Procrustes: rmse 0.0002209919  max resid 0.000439976 
    ## ... Similar to previous best
    ## Run 108 stress 9.76921e-05 
    ## ... Procrustes: rmse 0.0002322161  max resid 0.0004612434 
    ## ... Similar to previous best
    ## Run 109 stress 9.375891e-05 
    ## ... Procrustes: rmse 0.000220029  max resid 0.0004437906 
    ## ... Similar to previous best
    ## Run 110 stress 9.819031e-05 
    ## ... Procrustes: rmse 0.000243982  max resid 0.0004697405 
    ## ... Similar to previous best
    ## Run 111 stress 9.724068e-05 
    ## ... Procrustes: rmse 0.0001783083  max resid 0.0003426867 
    ## ... Similar to previous best
    ## Run 112 stress 8.66534e-05 
    ## ... Procrustes: rmse 0.0001351657  max resid 0.0002550435 
    ## ... Similar to previous best
    ## Run 113 stress 0.2680462 
    ## Run 114 stress 9.625029e-05 
    ## ... Procrustes: rmse 0.0002307042  max resid 0.0004629865 
    ## ... Similar to previous best
    ## Run 115 stress 9.956616e-05 
    ## ... Procrustes: rmse 0.0001793656  max resid 0.0002549127 
    ## ... Similar to previous best
    ## Run 116 stress 9.334986e-05 
    ## ... Procrustes: rmse 0.0002322359  max resid 0.0004531922 
    ## ... Similar to previous best
    ## Run 117 stress 9.568058e-05 
    ## ... Procrustes: rmse 0.0001728846  max resid 0.0002401169 
    ## ... Similar to previous best
    ## Run 118 stress 8.893345e-05 
    ## ... Procrustes: rmse 0.0001143459  max resid 0.0001737098 
    ## ... Similar to previous best
    ## Run 119 stress 9.575974e-05 
    ## ... Procrustes: rmse 0.0002267981  max resid 0.0004495966 
    ## ... Similar to previous best
    ## Run 120 stress 0.231917 
    ## Run 121 stress 9.176811e-05 
    ## ... Procrustes: rmse 0.0001440056  max resid 0.0002930766 
    ## ... Similar to previous best
    ## Run 122 stress 9.24235e-05 
    ## ... Procrustes: rmse 0.0001699716  max resid 0.0002401485 
    ## ... Similar to previous best
    ## Run 123 stress 0.2848526 
    ## Run 124 stress 0.3028823 
    ## Run 125 stress 8.371184e-05 
    ## ... Procrustes: rmse 0.0001539241  max resid 0.0002131014 
    ## ... Similar to previous best
    ## Run 126 stress 9.54848e-05 
    ## ... Procrustes: rmse 0.0001405626  max resid 0.0002651467 
    ## ... Similar to previous best
    ## Run 127 stress 9.588493e-05 
    ## ... Procrustes: rmse 0.0001304741  max resid 0.0002220977 
    ## ... Similar to previous best
    ## Run 128 stress 9.309413e-05 
    ## ... Procrustes: rmse 0.0001389027  max resid 0.0002522703 
    ## ... Similar to previous best
    ## Run 129 stress 7.517549e-05 
    ## ... Procrustes: rmse 0.0001288017  max resid 0.0002765842 
    ## ... Similar to previous best
    ## Run 130 stress 9.674656e-05 
    ## ... Procrustes: rmse 0.0001175807  max resid 0.0001863815 
    ## ... Similar to previous best
    ## Run 131 stress 9.545103e-05 
    ## ... Procrustes: rmse 0.0002301505  max resid 0.0004301907 
    ## ... Similar to previous best
    ## Run 132 stress 0.2848517 
    ## Run 133 stress 5.749241e-05 
    ## ... Procrustes: rmse 7.967796e-05  max resid 0.0001189434 
    ## ... Similar to previous best
    ## Run 134 stress 0.3075827 
    ## Run 135 stress 9.353704e-05 
    ## ... Procrustes: rmse 0.0001700932  max resid 0.0002362295 
    ## ... Similar to previous best
    ## Run 136 stress 9.634714e-05 
    ## ... Procrustes: rmse 0.0002421572  max resid 0.0004668405 
    ## ... Similar to previous best
    ## Run 137 stress 9.886711e-05 
    ## ... Procrustes: rmse 0.000236292  max resid 0.000469645 
    ## ... Similar to previous best
    ## Run 138 stress 8.287518e-05 
    ## ... Procrustes: rmse 0.000131666  max resid 0.0002340495 
    ## ... Similar to previous best
    ## Run 139 stress 9.43561e-05 
    ## ... Procrustes: rmse 0.0002344906  max resid 0.0004484366 
    ## ... Similar to previous best
    ## Run 140 stress 8.417335e-05 
    ## ... Procrustes: rmse 0.0001116298  max resid 0.0001628736 
    ## ... Similar to previous best
    ## Run 141 stress 9.827438e-05 
    ## ... Procrustes: rmse 0.0002295975  max resid 0.000455814 
    ## ... Similar to previous best
    ## Run 142 stress 9.9479e-05 
    ## ... Procrustes: rmse 0.0001791851  max resid 0.0002545637 
    ## ... Similar to previous best
    ## Run 143 stress 9.438415e-05 
    ## ... Procrustes: rmse 0.0001735469  max resid 0.0003154183 
    ## ... Similar to previous best
    ## Run 144 stress 9.33522e-05 
    ## ... Procrustes: rmse 0.0002283696  max resid 0.0004343056 
    ## ... Similar to previous best
    ## Run 145 stress 8.758915e-05 
    ## ... Procrustes: rmse 0.0002102012  max resid 0.000424674 
    ## ... Similar to previous best
    ## Run 146 stress 9.87883e-05 
    ## ... Procrustes: rmse 0.0002347681  max resid 0.0004656945 
    ## ... Similar to previous best
    ## Run 147 stress 0.2319145 
    ## Run 148 stress 9.821718e-05 
    ## ... Procrustes: rmse 0.0002216887  max resid 0.000452798 
    ## ... Similar to previous best
    ## Run 149 stress 9.473763e-05 
    ## ... Procrustes: rmse 0.0001829639  max resid 0.0003117352 
    ## ... Similar to previous best
    ## Run 150 stress 8.589123e-05 
    ## ... Procrustes: rmse 8.171667e-05  max resid 0.0001396158 
    ## ... Similar to previous best
    ## Run 151 stress 9.925656e-05 
    ## ... Procrustes: rmse 0.0002402497  max resid 0.0004498435 
    ## ... Similar to previous best
    ## Run 152 stress 9.817988e-05 
    ## ... Procrustes: rmse 0.0002425555  max resid 0.0004685155 
    ## ... Similar to previous best
    ## Run 153 stress 0.2848518 
    ## Run 154 stress 8.73921e-05 
    ## ... Procrustes: rmse 0.0002040902  max resid 0.000416377 
    ## ... Similar to previous best
    ## Run 155 stress 9.921127e-05 
    ## ... Procrustes: rmse 0.0001825589  max resid 0.000236688 
    ## ... Similar to previous best
    ## Run 156 stress 9.609125e-05 
    ## ... Procrustes: rmse 0.000239854  max resid 0.0004555269 
    ## ... Similar to previous best
    ## Run 157 stress 9.628309e-05 
    ## ... Procrustes: rmse 0.0001211991  max resid 0.0001762475 
    ## ... Similar to previous best
    ## Run 158 stress 9.368047e-05 
    ## ... Procrustes: rmse 0.0001736982  max resid 0.0003008931 
    ## ... Similar to previous best
    ## Run 159 stress 9.979078e-05 
    ## ... Procrustes: rmse 0.000234491  max resid 0.0004639483 
    ## ... Similar to previous best
    ## Run 160 stress 9.511038e-05 
    ## ... Procrustes: rmse 0.0002273928  max resid 0.0004603785 
    ## ... Similar to previous best
    ## Run 161 stress 9.919182e-05 
    ## ... Procrustes: rmse 0.0001776717  max resid 0.0002470415 
    ## ... Similar to previous best
    ## Run 162 stress 9.320271e-05 
    ## ... Procrustes: rmse 0.0002151749  max resid 0.0004234191 
    ## ... Similar to previous best
    ## Run 163 stress 9.021578e-05 
    ## ... Procrustes: rmse 0.0001454242  max resid 0.0002458755 
    ## ... Similar to previous best
    ## Run 164 stress 9.906811e-05 
    ## ... Procrustes: rmse 0.0001924783  max resid 0.0003366402 
    ## ... Similar to previous best
    ## Run 165 stress 9.497045e-05 
    ## ... Procrustes: rmse 0.0001660977  max resid 0.0003381022 
    ## ... Similar to previous best
    ## Run 166 stress 9.013892e-05 
    ## ... Procrustes: rmse 0.0001426007  max resid 0.0002690787 
    ## ... Similar to previous best
    ## Run 167 stress 9.305009e-05 
    ## ... Procrustes: rmse 0.0001450531  max resid 0.0002679641 
    ## ... Similar to previous best
    ## Run 168 stress 9.421809e-05 
    ## ... Procrustes: rmse 0.0001695819  max resid 0.0002420501 
    ## ... Similar to previous best
    ## Run 169 stress 9.267681e-05 
    ## ... Procrustes: rmse 0.0001353196  max resid 0.0002578954 
    ## ... Similar to previous best
    ## Run 170 stress 9.540707e-05 
    ## ... Procrustes: rmse 0.000151447  max resid 0.0002893528 
    ## ... Similar to previous best
    ## Run 171 stress 9.815595e-05 
    ## ... Procrustes: rmse 0.0001779192  max resid 0.0002475788 
    ## ... Similar to previous best
    ## Run 172 stress 9.122745e-05 
    ## ... Procrustes: rmse 0.000171064  max resid 0.0002353112 
    ## ... Similar to previous best
    ## Run 173 stress 0.3075827 
    ## Run 174 stress 9.529243e-05 
    ## ... Procrustes: rmse 0.0001813081  max resid 0.0002461351 
    ## ... Similar to previous best
    ## Run 175 stress 9.390948e-05 
    ## ... Procrustes: rmse 0.0002335026  max resid 0.0004547588 
    ## ... Similar to previous best
    ## Run 176 stress 9.621178e-05 
    ## ... Procrustes: rmse 0.0002304673  max resid 0.0004645889 
    ## ... Similar to previous best
    ## Run 177 stress 9.521891e-05 
    ## ... Procrustes: rmse 0.0001794479  max resid 0.0003329728 
    ## ... Similar to previous best
    ## Run 178 stress 9.021973e-05 
    ## ... Procrustes: rmse 0.0002147498  max resid 0.0004302309 
    ## ... Similar to previous best
    ## Run 179 stress 9.10051e-05 
    ## ... Procrustes: rmse 0.0001433446  max resid 0.0002719968 
    ## ... Similar to previous best
    ## Run 180 stress 9.67888e-05 
    ## ... Procrustes: rmse 0.0001742179  max resid 0.0002470566 
    ## ... Similar to previous best
    ## Run 181 stress 9.077046e-05 
    ## ... Procrustes: rmse 0.0002132601  max resid 0.0004276444 
    ## ... Similar to previous best
    ## Run 182 stress 9.160351e-05 
    ## ... Procrustes: rmse 0.0001680342  max resid 0.0002336775 
    ## ... Similar to previous best
    ## Run 183 stress 8.670072e-05 
    ## ... Procrustes: rmse 7.798076e-05  max resid 0.0001088443 
    ## ... Similar to previous best
    ## Run 184 stress 9.59916e-05 
    ## ... Procrustes: rmse 0.0002245074  max resid 0.0004489502 
    ## ... Similar to previous best
    ## Run 185 stress 8.408437e-05 
    ## ... Procrustes: rmse 0.0001102476  max resid 0.0001557965 
    ## ... Similar to previous best
    ## Run 186 stress 8.962002e-05 
    ## ... Procrustes: rmse 0.0002005211  max resid 0.0004060474 
    ## ... Similar to previous best
    ## Run 187 stress 9.716399e-05 
    ## ... Procrustes: rmse 0.000226424  max resid 0.0004493631 
    ## ... Similar to previous best
    ## Run 188 stress 9.814633e-05 
    ## ... Procrustes: rmse 0.0002351711  max resid 0.0004670788 
    ## ... Similar to previous best
    ## Run 189 stress 9.860642e-05 
    ## ... Procrustes: rmse 0.0001767088  max resid 0.0002470857 
    ## ... Similar to previous best
    ## Run 190 stress 9.086589e-05 
    ## ... Procrustes: rmse 0.0001178935  max resid 0.0002037444 
    ## ... Similar to previous best
    ## Run 191 stress 8.62619e-05 
    ## ... Procrustes: rmse 0.0001630302  max resid 0.0002205955 
    ## ... Similar to previous best
    ## Run 192 stress 9.277154e-05 
    ## ... Procrustes: rmse 0.0001830369  max resid 0.0003517498 
    ## ... Similar to previous best
    ## Run 193 stress 8.77926e-05 
    ## ... Procrustes: rmse 9.323429e-05  max resid 0.0001492507 
    ## ... Similar to previous best
    ## Run 194 stress 8.910312e-05 
    ## ... Procrustes: rmse 0.0002228837  max resid 0.0004303985 
    ## ... Similar to previous best
    ## Run 195 stress 9.485833e-05 
    ## ... Procrustes: rmse 0.0002268357  max resid 0.0004519855 
    ## ... Similar to previous best
    ## Run 196 stress 9.611132e-05 
    ## ... Procrustes: rmse 0.0001279303  max resid 0.0001911455 
    ## ... Similar to previous best
    ## Run 197 stress 9.86773e-05 
    ## ... Procrustes: rmse 0.00012791  max resid 0.0001793027 
    ## ... Similar to previous best
    ## Run 198 stress 8.01704e-05 
    ## ... Procrustes: rmse 0.0001086329  max resid 0.0001698366 
    ## ... Similar to previous best
    ## Run 199 stress 9.34942e-05 
    ## ... Procrustes: rmse 0.0002138877  max resid 0.000424055 
    ## ... Similar to previous best
    ## Run 200 stress 9.364243e-05 
    ## ... Procrustes: rmse 0.0001678718  max resid 0.0002365778 
    ## ... Similar to previous best
    ## Run 201 stress 8.83036e-05 
    ## ... Procrustes: rmse 0.0001164013  max resid 0.0001666685 
    ## ... Similar to previous best
    ## Run 202 stress 6.793447e-05 
    ## ... Procrustes: rmse 0.0001109194  max resid 0.0001854358 
    ## ... Similar to previous best
    ## Run 203 stress 9.926141e-05 
    ## ... Procrustes: rmse 0.000179923  max resid 0.0002504607 
    ## ... Similar to previous best
    ## Run 204 stress 9.778823e-05 
    ## ... Procrustes: rmse 0.0002433513  max resid 0.0004705934 
    ## ... Similar to previous best
    ## Run 205 stress 9.734934e-05 
    ## ... Procrustes: rmse 0.0001802979  max resid 0.0002492685 
    ## ... Similar to previous best
    ## Run 206 stress 9.449453e-05 
    ## ... Procrustes: rmse 0.0002304218  max resid 0.0004530042 
    ## ... Similar to previous best
    ## Run 207 stress 9.572982e-05 
    ## ... Procrustes: rmse 0.0002253003  max resid 0.0004511164 
    ## ... Similar to previous best
    ## Run 208 stress 8.569355e-05 
    ## ... Procrustes: rmse 0.0001439515  max resid 0.0002431802 
    ## ... Similar to previous best
    ## Run 209 stress 9.541293e-05 
    ## ... Procrustes: rmse 0.0002005106  max resid 0.000411073 
    ## ... Similar to previous best
    ## Run 210 stress 9.922277e-05 
    ## ... Procrustes: rmse 0.0001794138  max resid 0.0002495731 
    ## ... Similar to previous best
    ## Run 211 stress 9.673991e-05 
    ## ... Procrustes: rmse 0.0002360525  max resid 0.0004422497 
    ## ... Similar to previous best
    ## Run 212 stress 9.741607e-05 
    ## ... Procrustes: rmse 0.0002316758  max resid 0.0004582478 
    ## ... Similar to previous best
    ## Run 213 stress 9.367655e-05 
    ## ... Procrustes: rmse 0.0002171507  max resid 0.0004393603 
    ## ... Similar to previous best
    ## Run 214 stress 9.421055e-05 
    ## ... Procrustes: rmse 0.000147016  max resid 0.0002805813 
    ## ... Similar to previous best
    ## Run 215 stress 0.2848514 
    ## Run 216 stress 8.701534e-05 
    ## ... Procrustes: rmse 0.0002027271  max resid 0.0004074061 
    ## ... Similar to previous best
    ## Run 217 stress 9.330126e-05 
    ## ... Procrustes: rmse 0.0001521786  max resid 0.000251627 
    ## ... Similar to previous best
    ## Run 218 stress 9.612959e-05 
    ## ... Procrustes: rmse 0.0002411148  max resid 0.0004635716 
    ## ... Similar to previous best
    ## Run 219 stress 9.970782e-05 
    ## ... Procrustes: rmse 0.0001579739  max resid 0.0002969528 
    ## ... Similar to previous best
    ## Run 220 stress 9.715149e-05 
    ## ... Procrustes: rmse 0.0001834008  max resid 0.0003209493 
    ## ... Similar to previous best
    ## Run 221 stress 9.407025e-05 
    ## ... Procrustes: rmse 0.0002292573  max resid 0.000439375 
    ## ... Similar to previous best
    ## Run 222 stress 9.206094e-05 
    ## ... Procrustes: rmse 8.499574e-05  max resid 0.0001163495 
    ## ... Similar to previous best
    ## Run 223 stress 9.945099e-05 
    ## ... Procrustes: rmse 0.0002270785  max resid 0.0004493329 
    ## ... Similar to previous best
    ## Run 224 stress 8.181533e-05 
    ## ... Procrustes: rmse 0.00019039  max resid 0.000388768 
    ## ... Similar to previous best
    ## Run 225 stress 8.963064e-05 
    ## ... Procrustes: rmse 0.0001116312  max resid 0.0001790318 
    ## ... Similar to previous best
    ## Run 226 stress 9.889891e-05 
    ## ... Procrustes: rmse 0.0001817257  max resid 0.000254224 
    ## ... Similar to previous best
    ## Run 227 stress 9.625375e-05 
    ## ... Procrustes: rmse 0.0002134151  max resid 0.0004271039 
    ## ... Similar to previous best
    ## Run 228 stress 9.859316e-05 
    ## ... Procrustes: rmse 0.0002362923  max resid 0.0004705446 
    ## ... Similar to previous best
    ## Run 229 stress 9.752407e-05 
    ## ... Procrustes: rmse 0.0001549616  max resid 0.00025526 
    ## ... Similar to previous best
    ## Run 230 stress 8.522342e-05 
    ## ... Procrustes: rmse 0.0001862096  max resid 0.0003731163 
    ## ... Similar to previous best
    ## Run 231 stress 9.531833e-05 
    ## ... Procrustes: rmse 0.0002284951  max resid 0.000458584 
    ## ... Similar to previous best
    ## Run 232 stress 9.430122e-05 
    ## ... Procrustes: rmse 0.0002190318  max resid 0.0004381096 
    ## ... Similar to previous best
    ## Run 233 stress 9.961393e-05 
    ## ... Procrustes: rmse 0.0002456929  max resid 0.0004697769 
    ## ... Similar to previous best
    ## Run 234 stress 0.2816758 
    ## Run 235 stress 9.595081e-05 
    ## ... Procrustes: rmse 0.0002279195  max resid 0.0004532785 
    ## ... Similar to previous best
    ## Run 236 stress 9.076077e-05 
    ## ... Procrustes: rmse 0.0001284342  max resid 0.0002387884 
    ## ... Similar to previous best
    ## Run 237 stress 9.675846e-05 
    ## ... Procrustes: rmse 0.000230312  max resid 0.0004618467 
    ## ... Similar to previous best
    ## Run 238 stress 8.994576e-05 
    ## ... Procrustes: rmse 6.723993e-05  max resid 9.553864e-05 
    ## ... Similar to previous best
    ## Run 239 stress 9.878093e-05 
    ## ... Procrustes: rmse 0.0002475887  max resid 0.0004739258 
    ## ... Similar to previous best
    ## Run 240 stress 9.368232e-05 
    ## ... Procrustes: rmse 0.0001699578  max resid 0.0002745353 
    ## ... Similar to previous best
    ## Run 241 stress 8.597687e-05 
    ## ... Procrustes: rmse 0.0001079535  max resid 0.0001609955 
    ## ... Similar to previous best
    ## Run 242 stress 8.944002e-05 
    ## ... Procrustes: rmse 0.0002199346  max resid 0.0004284113 
    ## ... Similar to previous best
    ## Run 243 stress 9.412085e-05 
    ## ... Procrustes: rmse 0.0002354321  max resid 0.0004503659 
    ## ... Similar to previous best
    ## Run 244 stress 9.493481e-05 
    ## ... Procrustes: rmse 0.0001491546  max resid 0.0002695042 
    ## ... Similar to previous best
    ## Run 245 stress 9.664682e-05 
    ## ... Procrustes: rmse 0.000168228  max resid 0.0002486951 
    ## ... Similar to previous best
    ## Run 246 stress 9.017138e-05 
    ## ... Procrustes: rmse 0.0002159122  max resid 0.0004327125 
    ## ... Similar to previous best
    ## Run 247 stress 9.756084e-05 
    ## ... Procrustes: rmse 0.0001286514  max resid 0.0001913322 
    ## ... Similar to previous best
    ## Run 248 stress 9.025043e-05 
    ## ... Procrustes: rmse 0.0002083433  max resid 0.0004232229 
    ## ... Similar to previous best
    ## Run 249 stress 9.049335e-05 
    ## ... Procrustes: rmse 0.0001795262  max resid 0.0003359335 
    ## ... Similar to previous best
    ## Run 250 stress 9.001207e-05 
    ## ... Procrustes: rmse 0.0001602202  max resid 0.0002267468 
    ## ... Similar to previous best
    ## Run 251 stress 9.962828e-05 
    ## ... Procrustes: rmse 0.0001774921  max resid 0.0002448701 
    ## ... Similar to previous best
    ## Run 252 stress 8.757434e-05 
    ## ... Procrustes: rmse 0.0001351653  max resid 0.0002288314 
    ## ... Similar to previous best
    ## Run 253 stress 9.132256e-05 
    ## ... Procrustes: rmse 0.0001006596  max resid 0.0001576982 
    ## ... Similar to previous best
    ## Run 254 stress 9.13371e-05 
    ## ... Procrustes: rmse 0.0002218941  max resid 0.0004237832 
    ## ... Similar to previous best
    ## Run 255 stress 9.762264e-05 
    ## ... Procrustes: rmse 0.0002346101  max resid 0.0004706993 
    ## ... Similar to previous best
    ## Run 256 stress 9.476175e-05 
    ## ... Procrustes: rmse 0.0001649284  max resid 0.0003499866 
    ## ... Similar to previous best
    ## Run 257 stress 9.5182e-05 
    ## ... Procrustes: rmse 0.000152157  max resid 0.0002906407 
    ## ... Similar to previous best
    ## Run 258 stress 9.53467e-05 
    ## ... Procrustes: rmse 0.0002283027  max resid 0.0004585668 
    ## ... Similar to previous best
    ## Run 259 stress 9.439437e-05 
    ## ... Procrustes: rmse 0.0001720902  max resid 0.0002437038 
    ## ... Similar to previous best
    ## Run 260 stress 9.796094e-05 
    ## ... Procrustes: rmse 0.0002285146  max resid 0.0004615023 
    ## ... Similar to previous best
    ## Run 261 stress 0.3083098 
    ## Run 262 stress 9.953198e-05 
    ## ... Procrustes: rmse 0.0002385608  max resid 0.0004731532 
    ## ... Similar to previous best
    ## Run 263 stress 9.918374e-05 
    ## ... Procrustes: rmse 0.0001808966  max resid 0.0002536167 
    ## ... Similar to previous best
    ## Run 264 stress 9.507748e-05 
    ## ... Procrustes: rmse 0.0002257121  max resid 0.0004572281 
    ## ... Similar to previous best
    ## Run 265 stress 9.841195e-05 
    ## ... Procrustes: rmse 0.00016648  max resid 0.0002826388 
    ## ... Similar to previous best
    ## Run 266 stress 8.338092e-05 
    ## ... Procrustes: rmse 0.0001336846  max resid 0.0002482693 
    ## ... Similar to previous best
    ## Run 267 stress 8.645795e-05 
    ## ... Procrustes: rmse 0.0001345674  max resid 0.0002426457 
    ## ... Similar to previous best
    ## Run 268 stress 9.568128e-05 
    ## ... Procrustes: rmse 0.0002258756  max resid 0.0004485442 
    ## ... Similar to previous best
    ## Run 269 stress 9.845376e-05 
    ## ... Procrustes: rmse 0.0002354572  max resid 0.0004674536 
    ## ... Similar to previous best
    ## Run 270 stress 9.486581e-05 
    ## ... Procrustes: rmse 0.0001719049  max resid 0.0002388843 
    ## ... Similar to previous best
    ## Run 271 stress 9.44594e-05 
    ## ... Procrustes: rmse 0.0001752266  max resid 0.0003414527 
    ## ... Similar to previous best
    ## Run 272 stress 9.04928e-05 
    ## ... Procrustes: rmse 0.0002228306  max resid 0.00043181 
    ## ... Similar to previous best
    ## Run 273 stress 9.897642e-05 
    ## ... Procrustes: rmse 0.0002359856  max resid 0.000467771 
    ## ... Similar to previous best
    ## Run 274 stress 9.845156e-05 
    ## ... Procrustes: rmse 0.0002371013  max resid 0.000473575 
    ## ... Similar to previous best
    ## Run 275 stress 9.41142e-05 
    ## ... Procrustes: rmse 0.0001240743  max resid 0.0001865952 
    ## ... Similar to previous best
    ## Run 276 stress 8.850072e-05 
    ## ... Procrustes: rmse 0.0001654964  max resid 0.0002314165 
    ## ... Similar to previous best
    ## Run 277 stress 9.442948e-05 
    ## ... Procrustes: rmse 0.0002206552  max resid 0.0004447477 
    ## ... Similar to previous best
    ## Run 278 stress 9.930042e-05 
    ## ... Procrustes: rmse 0.000180329  max resid 0.00025603 
    ## ... Similar to previous best
    ## Run 279 stress 0.3082558 
    ## Run 280 stress 9.217554e-05 
    ## ... Procrustes: rmse 0.0001730991  max resid 0.0003116213 
    ## ... Similar to previous best
    ## Run 281 stress 9.183104e-05 
    ## ... Procrustes: rmse 0.000230499  max resid 0.0004441385 
    ## ... Similar to previous best
    ## Run 282 stress 9.443749e-05 
    ## ... Procrustes: rmse 0.000226356  max resid 0.0004562126 
    ## ... Similar to previous best
    ## Run 283 stress 0.2848512 
    ## Run 284 stress 9.792373e-05 
    ## ... Procrustes: rmse 0.0001210566  max resid 0.0001866073 
    ## ... Similar to previous best
    ## Run 285 stress 9.764702e-05 
    ## ... Procrustes: rmse 0.0001492564  max resid 0.0002746126 
    ## ... Similar to previous best
    ## Run 286 stress 9.412038e-05 
    ## ... Procrustes: rmse 0.0001699809  max resid 0.0002359088 
    ## ... Similar to previous best
    ## Run 287 stress 9.80908e-05 
    ## ... Procrustes: rmse 0.0002360224  max resid 0.0004731793 
    ## ... Similar to previous best
    ## Run 288 stress 9.153506e-05 
    ## ... Procrustes: rmse 0.0001656102  max resid 0.0002296286 
    ## ... Similar to previous best
    ## Run 289 stress 9.938769e-05 
    ## ... Procrustes: rmse 0.0001866563  max resid 0.0003312697 
    ## ... Similar to previous best
    ## Run 290 stress 8.887679e-05 
    ## ... Procrustes: rmse 0.0001952281  max resid 0.0004024371 
    ## ... Similar to previous best
    ## Run 291 stress 9.944615e-05 
    ## ... Procrustes: rmse 0.0002492967  max resid 0.0004739232 
    ## ... Similar to previous best
    ## Run 292 stress 9.283161e-05 
    ## ... Procrustes: rmse 0.0001258499  max resid 0.0001833543 
    ## ... Similar to previous best
    ## Run 293 stress 9.358782e-05 
    ## ... Procrustes: rmse 0.0001535996  max resid 0.000269019 
    ## ... Similar to previous best
    ## Run 294 stress 8.87486e-05 
    ## ... Procrustes: rmse 0.0001729786  max resid 0.0002399379 
    ## ... Similar to previous best
    ## Run 295 stress 9.627147e-05 
    ## ... Procrustes: rmse 0.0002407107  max resid 0.0004643107 
    ## ... Similar to previous best
    ## Run 296 stress 8.734364e-05 
    ## ... Procrustes: rmse 0.0002181366  max resid 0.0004259949 
    ## ... Similar to previous best
    ## Run 297 stress 9.864412e-05 
    ## ... Procrustes: rmse 0.0002348686  max resid 0.0004659055 
    ## ... Similar to previous best
    ## Run 298 stress 9.222511e-05 
    ## ... Procrustes: rmse 0.000140182  max resid 0.000265714 
    ## ... Similar to previous best
    ## Run 299 stress 9.83588e-05 
    ## ... Procrustes: rmse 0.0002311403  max resid 0.000463738 
    ## ... Similar to previous best
    ## Run 300 stress 9.526266e-05 
    ## ... Procrustes: rmse 0.0002181  max resid 0.000432414 
    ## ... Similar to previous best
    ## Run 301 stress 9.076432e-05 
    ## ... Procrustes: rmse 0.00016226  max resid 0.0002230484 
    ## ... Similar to previous best
    ## Run 302 stress 8.747278e-05 
    ## ... Procrustes: rmse 0.0001471437  max resid 0.000247363 
    ## ... Similar to previous best
    ## Run 303 stress 7.874647e-05 
    ## ... Procrustes: rmse 0.0001119816  max resid 0.0001844268 
    ## ... Similar to previous best
    ## Run 304 stress 9.494587e-05 
    ## ... Procrustes: rmse 0.0002241615  max resid 0.0004464296 
    ## ... Similar to previous best
    ## Run 305 stress 9.20549e-05 
    ## ... Procrustes: rmse 0.0001668539  max resid 0.0002387839 
    ## ... Similar to previous best
    ## Run 306 stress 9.133546e-05 
    ## ... Procrustes: rmse 0.0001658525  max resid 0.000230913 
    ## ... Similar to previous best
    ## Run 307 stress 9.855599e-05 
    ## ... Procrustes: rmse 0.0001856474  max resid 0.0003337557 
    ## ... Similar to previous best
    ## Run 308 stress 9.844534e-05 
    ## ... Procrustes: rmse 0.000199234  max resid 0.0003817724 
    ## ... Similar to previous best
    ## Run 309 stress 9.616023e-05 
    ## ... Procrustes: rmse 0.0001877624  max resid 0.0003515318 
    ## ... Similar to previous best
    ## Run 310 stress 9.799954e-05 
    ## ... Procrustes: rmse 0.0001469935  max resid 0.0002661504 
    ## ... Similar to previous best
    ## Run 311 stress 9.263132e-05 
    ## ... Procrustes: rmse 5.722562e-05  max resid 9.082992e-05 
    ## ... Similar to previous best
    ## Run 312 stress 9.206711e-05 
    ## ... Procrustes: rmse 0.0002124879  max resid 0.0004265099 
    ## ... Similar to previous best
    ## Run 313 stress 9.430327e-05 
    ## ... Procrustes: rmse 0.0002303358  max resid 0.0004379902 
    ## ... Similar to previous best
    ## Run 314 stress 9.606746e-05 
    ## ... Procrustes: rmse 0.0002249504  max resid 0.0004560841 
    ## ... Similar to previous best
    ## Run 315 stress 8.91846e-05 
    ## ... Procrustes: rmse 0.00012126  max resid 0.0001979583 
    ## ... Similar to previous best
    ## Run 316 stress 8.523165e-05 
    ## ... Procrustes: rmse 0.0001601093  max resid 0.0002842995 
    ## ... Similar to previous best
    ## Run 317 stress 9.494175e-05 
    ## ... Procrustes: rmse 0.0001542027  max resid 0.0002210758 
    ## ... Similar to previous best
    ## Run 318 stress 0.2319192 
    ## Run 319 stress 9.339557e-05 
    ## ... Procrustes: rmse 0.0002348076  max resid 0.000451297 
    ## ... Similar to previous best
    ## Run 320 stress 9.493329e-05 
    ## ... Procrustes: rmse 0.0002264712  max resid 0.0004518851 
    ## ... Similar to previous best
    ## Run 321 stress 8.019117e-05 
    ## ... Procrustes: rmse 0.000188217  max resid 0.0003865244 
    ## ... Similar to previous best
    ## Run 322 stress 0.2971485 
    ## Run 323 stress 9.365475e-05 
    ## ... Procrustes: rmse 0.0002218354  max resid 0.0004456556 
    ## ... Similar to previous best
    ## Run 324 stress 9.215466e-05 
    ## ... Procrustes: rmse 0.0001685244  max resid 0.0002243069 
    ## ... Similar to previous best
    ## Run 325 stress 9.575732e-05 
    ## ... Procrustes: rmse 0.0001514869  max resid 0.000263683 
    ## ... Similar to previous best
    ## Run 326 stress 9.89125e-05 
    ## ... Procrustes: rmse 0.0002253727  max resid 0.0004523679 
    ## ... Similar to previous best
    ## Run 327 stress 9.31887e-05 
    ## ... Procrustes: rmse 0.0002191129  max resid 0.0004375881 
    ## ... Similar to previous best
    ## Run 328 stress 9.530236e-05 
    ## ... Procrustes: rmse 0.0002375323  max resid 0.0004605775 
    ## ... Similar to previous best
    ## Run 329 stress 9.742179e-05 
    ## ... Procrustes: rmse 0.0002232367  max resid 0.000463025 
    ## ... Similar to previous best
    ## Run 330 stress 9.34344e-05 
    ## ... Procrustes: rmse 0.0002180617  max resid 0.0004375203 
    ## ... Similar to previous best
    ## Run 331 stress 0.2319133 
    ## Run 332 stress 9.616977e-05 
    ## ... Procrustes: rmse 0.0002301962  max resid 0.0004440839 
    ## ... Similar to previous best
    ## Run 333 stress 8.552861e-05 
    ## ... Procrustes: rmse 7.656982e-05  max resid 0.0001074023 
    ## ... Similar to previous best
    ## Run 334 stress 9.092319e-05 
    ## ... Procrustes: rmse 0.0001215628  max resid 0.0002244079 
    ## ... Similar to previous best
    ## Run 335 stress 9.131768e-05 
    ## ... Procrustes: rmse 0.0002278903  max resid 0.0004380472 
    ## ... Similar to previous best
    ## Run 336 stress 9.424592e-05 
    ## ... Procrustes: rmse 0.0001834437  max resid 0.0003338056 
    ## ... Similar to previous best
    ## Run 337 stress 9.517842e-05 
    ## ... Procrustes: rmse 0.0002213573  max resid 0.0004460516 
    ## ... Similar to previous best
    ## Run 338 stress 9.130905e-05 
    ## ... Procrustes: rmse 0.0002025147  max resid 0.0004095159 
    ## ... Similar to previous best
    ## Run 339 stress 9.517524e-05 
    ## ... Procrustes: rmse 0.0002267183  max resid 0.0004505541 
    ## ... Similar to previous best
    ## Run 340 stress 9.692175e-05 
    ## ... Procrustes: rmse 0.0002406879  max resid 0.0004624703 
    ## ... Similar to previous best
    ## Run 341 stress 9.969002e-05 
    ## ... Procrustes: rmse 0.0001249584  max resid 0.0001972721 
    ## ... Similar to previous best
    ## Run 342 stress 9.604439e-05 
    ## ... Procrustes: rmse 0.0002274751  max resid 0.0004523343 
    ## ... Similar to previous best
    ## Run 343 stress 9.728234e-05 
    ## ... Procrustes: rmse 0.0002306152  max resid 0.0004636952 
    ## ... Similar to previous best
    ## Run 344 stress 9.567528e-05 
    ## ... Procrustes: rmse 0.0001829767  max resid 0.0003370651 
    ## ... Similar to previous best
    ## Run 345 stress 8.779107e-05 
    ## ... Procrustes: rmse 0.0001582504  max resid 0.00022432 
    ## ... Similar to previous best
    ## Run 346 stress 9.429537e-05 
    ## ... Procrustes: rmse 0.0002371952  max resid 0.000455665 
    ## ... Similar to previous best
    ## Run 347 stress 9.226538e-05 
    ## ... Procrustes: rmse 0.0001244798  max resid 0.0001700786 
    ## ... Similar to previous best
    ## Run 348 stress 9.645106e-05 
    ## ... Procrustes: rmse 0.0001911214  max resid 0.0003748243 
    ## ... Similar to previous best
    ## Run 349 stress 9.743252e-05 
    ## ... Procrustes: rmse 0.0001255342  max resid 0.0001739714 
    ## ... Similar to previous best
    ## Run 350 stress 8.636328e-05 
    ## ... Procrustes: rmse 0.0002116822  max resid 0.0004124524 
    ## ... Similar to previous best
    ## Run 351 stress 9.591818e-05 
    ## ... Procrustes: rmse 0.0001223127  max resid 0.0001729852 
    ## ... Similar to previous best
    ## Run 352 stress 9.408227e-05 
    ## ... Procrustes: rmse 0.0002240277  max resid 0.0004466513 
    ## ... Similar to previous best
    ## Run 353 stress 8.890178e-05 
    ## ... Procrustes: rmse 0.0002197788  max resid 0.0004245063 
    ## ... Similar to previous best
    ## Run 354 stress 9.054052e-05 
    ## ... Procrustes: rmse 0.0002157908  max resid 0.0004338191 
    ## ... Similar to previous best
    ## Run 355 stress 9.49623e-05 
    ## ... Procrustes: rmse 0.0002241384  max resid 0.0004453388 
    ## ... Similar to previous best
    ## Run 356 stress 8.404676e-05 
    ## ... Procrustes: rmse 0.0001317491  max resid 0.0002493165 
    ## ... Similar to previous best
    ## Run 357 stress 9.433228e-05 
    ## ... Procrustes: rmse 0.0002221677  max resid 0.0004424742 
    ## ... Similar to previous best
    ## Run 358 stress 9.512565e-05 
    ## ... Procrustes: rmse 0.000221614  max resid 0.0004197717 
    ## ... Similar to previous best
    ## Run 359 stress 9.697187e-05 
    ## ... Procrustes: rmse 0.000230572  max resid 0.0004626055 
    ## ... Similar to previous best
    ## Run 360 stress 9.66305e-05 
    ## ... Procrustes: rmse 0.0002247026  max resid 0.0004528874 
    ## ... Similar to previous best
    ## Run 361 stress 9.966219e-05 
    ## ... Procrustes: rmse 0.0001830905  max resid 0.0003618692 
    ## ... Similar to previous best
    ## Run 362 stress 9.252794e-05 
    ## ... Procrustes: rmse 0.0002152009  max resid 0.0004289276 
    ## ... Similar to previous best
    ## Run 363 stress 9.751125e-05 
    ## ... Procrustes: rmse 0.000234408  max resid 0.0004655987 
    ## ... Similar to previous best
    ## Run 364 stress 9.952375e-05 
    ## ... Procrustes: rmse 0.0002449715  max resid 0.000464491 
    ## ... Similar to previous best
    ## Run 365 stress 9.761618e-05 
    ## ... Procrustes: rmse 0.0002342749  max resid 0.0004665072 
    ## ... Similar to previous best
    ## Run 366 stress 9.695019e-05 
    ## ... Procrustes: rmse 0.0002372403  max resid 0.0004504174 
    ## ... Similar to previous best
    ## Run 367 stress 9.930557e-05 
    ## ... Procrustes: rmse 0.0002346196  max resid 0.000465038 
    ## ... Similar to previous best
    ## Run 368 stress 9.992376e-05 
    ## ... Procrustes: rmse 0.0002324743  max resid 0.0004594376 
    ## ... Similar to previous best
    ## Run 369 stress 9.430863e-05 
    ## ... Procrustes: rmse 0.0002353539  max resid 0.0004564472 
    ## ... Similar to previous best
    ## Run 370 stress 8.884912e-05 
    ## ... Procrustes: rmse 0.0001377418  max resid 0.0002520076 
    ## ... Similar to previous best
    ## Run 371 stress 9.269763e-05 
    ## ... Procrustes: rmse 0.0001494665  max resid 0.0002607817 
    ## ... Similar to previous best
    ## Run 372 stress 8.335828e-05 
    ## ... Procrustes: rmse 0.0001066468  max resid 0.0001882648 
    ## ... Similar to previous best
    ## Run 373 stress 9.371148e-05 
    ## ... Procrustes: rmse 0.0001791816  max resid 0.0003116904 
    ## ... Similar to previous best
    ## Run 374 stress 9.286588e-05 
    ## ... Procrustes: rmse 0.0002186646  max resid 0.0004381796 
    ## ... Similar to previous best
    ## Run 375 stress 9.62836e-05 
    ## ... Procrustes: rmse 0.0002333899  max resid 0.0004552947 
    ## ... Similar to previous best
    ## Run 376 stress 8.60595e-05 
    ## ... Procrustes: rmse 0.0001452321  max resid 0.0002766598 
    ## ... Similar to previous best
    ## Run 377 stress 9.663488e-05 
    ## ... Procrustes: rmse 0.0002295996  max resid 0.0004324284 
    ## ... Similar to previous best
    ## Run 378 stress 9.224086e-05 
    ## ... Procrustes: rmse 0.0002045764  max resid 0.0004146582 
    ## ... Similar to previous best
    ## Run 379 stress 9.58358e-05 
    ## ... Procrustes: rmse 0.0001208509  max resid 0.0001926715 
    ## ... Similar to previous best
    ## Run 380 stress 9.922456e-05 
    ## ... Procrustes: rmse 0.0001756848  max resid 0.0003602927 
    ## ... Similar to previous best
    ## Run 381 stress 9.364477e-05 
    ## ... Procrustes: rmse 0.0001691768  max resid 0.0002397256 
    ## ... Similar to previous best
    ## Run 382 stress 9.941836e-05 
    ## ... Procrustes: rmse 0.0001851955  max resid 0.0002455478 
    ## ... Similar to previous best
    ## Run 383 stress 9.622135e-05 
    ## ... Procrustes: rmse 0.0001833113  max resid 0.0003358439 
    ## ... Similar to previous best
    ## Run 384 stress 9.992662e-05 
    ## ... Procrustes: rmse 0.0001718636  max resid 0.0002547856 
    ## ... Similar to previous best
    ## Run 385 stress 9.933784e-05 
    ## ... Procrustes: rmse 0.0002443534  max resid 0.0004690655 
    ## ... Similar to previous best
    ## Run 386 stress 9.33999e-05 
    ## ... Procrustes: rmse 0.0001694382  max resid 0.0002351911 
    ## ... Similar to previous best
    ## Run 387 stress 9.999842e-05 
    ## ... Procrustes: rmse 0.0001949117  max resid 0.0003435208 
    ## ... Similar to previous best
    ## Run 388 stress 9.739836e-05 
    ## ... Procrustes: rmse 0.0002307224  max resid 0.0004582286 
    ## ... Similar to previous best
    ## Run 389 stress 9.638159e-05 
    ## ... Procrustes: rmse 0.0002262331  max resid 0.0004508714 
    ## ... Similar to previous best
    ## Run 390 stress 0.2319142 
    ## Run 391 stress 9.776486e-05 
    ## ... Procrustes: rmse 0.0002353526  max resid 0.0004685015 
    ## ... Similar to previous best
    ## Run 392 stress 9.296097e-05 
    ## ... Procrustes: rmse 0.0002162676  max resid 0.0004309219 
    ## ... Similar to previous best
    ## Run 393 stress 9.748026e-05 
    ## ... Procrustes: rmse 0.0001793037  max resid 0.0003458132 
    ## ... Similar to previous best
    ## Run 394 stress 0.2680462 
    ## Run 395 stress 9.378622e-05 
    ## ... Procrustes: rmse 0.0002259498  max resid 0.0004525324 
    ## ... Similar to previous best
    ## Run 396 stress 9.546665e-05 
    ## ... Procrustes: rmse 0.0001730709  max resid 0.0002405973 
    ## ... Similar to previous best
    ## Run 397 stress 9.453062e-05 
    ## ... Procrustes: rmse 0.0001571497  max resid 0.0002917031 
    ## ... Similar to previous best
    ## Run 398 stress 9.760381e-05 
    ## ... Procrustes: rmse 0.0002452345  max resid 0.0004671836 
    ## ... Similar to previous best
    ## Run 399 stress 9.051808e-05 
    ## ... Procrustes: rmse 0.0002142591  max resid 0.0004346474 
    ## ... Similar to previous best
    ## Run 400 stress 9.3788e-05 
    ## ... Procrustes: rmse 0.0002190588  max resid 0.0004360507 
    ## ... Similar to previous best
    ## Run 401 stress 9.043161e-05 
    ## ... Procrustes: rmse 0.0001574706  max resid 0.0002125571 
    ## ... Similar to previous best
    ## Run 402 stress 9.576514e-05 
    ## ... Procrustes: rmse 0.0002231044  max resid 0.0004430355 
    ## ... Similar to previous best
    ## Run 403 stress 9.840707e-05 
    ## ... Procrustes: rmse 0.0002366486  max resid 0.0004701866 
    ## ... Similar to previous best
    ## Run 404 stress 9.340108e-05 
    ## ... Procrustes: rmse 0.000231209  max resid 0.0004491733 
    ## ... Similar to previous best
    ## Run 405 stress 9.329824e-05 
    ## ... Procrustes: rmse 0.0002300919  max resid 0.000451243 
    ## ... Similar to previous best
    ## Run 406 stress 9.627932e-05 
    ## ... Procrustes: rmse 0.0001631134  max resid 0.0003028745 
    ## ... Similar to previous best
    ## Run 407 stress 9.843355e-05 
    ## ... Procrustes: rmse 0.0001920409  max resid 0.0003567379 
    ## ... Similar to previous best
    ## Run 408 stress 9.252793e-05 
    ## ... Procrustes: rmse 0.0001664232  max resid 0.0002358663 
    ## ... Similar to previous best
    ## Run 409 stress 8.810241e-05 
    ## ... Procrustes: rmse 9.390661e-05  max resid 0.0001674603 
    ## ... Similar to previous best
    ## Run 410 stress 9.666009e-05 
    ## ... Procrustes: rmse 0.0001603286  max resid 0.0002369135 
    ## ... Similar to previous best
    ## Run 411 stress 9.087618e-05 
    ## ... Procrustes: rmse 0.0001282575  max resid 0.0002010052 
    ## ... Similar to previous best
    ## Run 412 stress 9.981397e-05 
    ## ... Procrustes: rmse 0.0001795  max resid 0.0002552252 
    ## ... Similar to previous best
    ## Run 413 stress 9.291915e-05 
    ## ... Procrustes: rmse 0.0001686132  max resid 0.0002342242 
    ## ... Similar to previous best
    ## Run 414 stress 8.999087e-05 
    ## ... Procrustes: rmse 0.0001657635  max resid 0.000231197 
    ## ... Similar to previous best
    ## Run 415 stress 9.091369e-05 
    ## ... Procrustes: rmse 0.0001818813  max resid 0.0003375697 
    ## ... Similar to previous best
    ## Run 416 stress 9.151815e-05 
    ## ... Procrustes: rmse 0.0001422846  max resid 0.0002614452 
    ## ... Similar to previous best
    ## Run 417 stress 0.3083099 
    ## Run 418 stress 9.227248e-05 
    ## ... Procrustes: rmse 0.0002248022  max resid 0.0004282002 
    ## ... Similar to previous best
    ## Run 419 stress 9.5857e-05 
    ## ... Procrustes: rmse 0.0001291552  max resid 0.0001887461 
    ## ... Similar to previous best
    ## Run 420 stress 9.56131e-05 
    ## ... Procrustes: rmse 0.0001771398  max resid 0.0003445292 
    ## ... Similar to previous best
    ## Run 421 stress 9.801629e-05 
    ## ... Procrustes: rmse 0.0002453869  max resid 0.0004755682 
    ## ... Similar to previous best
    ## Run 422 stress 9.672943e-05 
    ## ... Procrustes: rmse 0.0001612575  max resid 0.0002350461 
    ## ... Similar to previous best
    ## Run 423 stress 0.2680462 
    ## Run 424 stress 8.341616e-05 
    ## ... Procrustes: rmse 0.0001660243  max resid 0.0003110593 
    ## ... Similar to previous best
    ## Run 425 stress 9.312058e-05 
    ## ... Procrustes: rmse 0.0001766556  max resid 0.0003234971 
    ## ... Similar to previous best
    ## Run 426 stress 9.382119e-05 
    ## ... Procrustes: rmse 0.0002202528  max resid 0.0004469754 
    ## ... Similar to previous best
    ## Run 427 stress 9.229735e-05 
    ## ... Procrustes: rmse 0.0001674156  max resid 0.0002324706 
    ## ... Similar to previous best
    ## Run 428 stress 8.556303e-05 
    ## ... Procrustes: rmse 0.0001082124  max resid 0.0001657637 
    ## ... Similar to previous best
    ## Run 429 stress 9.373421e-05 
    ## ... Procrustes: rmse 0.0001345565  max resid 0.0002514753 
    ## ... Similar to previous best
    ## Run 430 stress 9.60384e-05 
    ## ... Procrustes: rmse 0.0002250401  max resid 0.0004467557 
    ## ... Similar to previous best
    ## Run 431 stress 9.626883e-05 
    ## ... Procrustes: rmse 0.0002297542  max resid 0.0004545506 
    ## ... Similar to previous best
    ## Run 432 stress 9.228831e-05 
    ## ... Procrustes: rmse 0.0001447473  max resid 0.0002694026 
    ## ... Similar to previous best
    ## Run 433 stress 9.269027e-05 
    ## ... Procrustes: rmse 0.0002180018  max resid 0.0004421127 
    ## ... Similar to previous best
    ## Run 434 stress 0.23192 
    ## Run 435 stress 9.04138e-05 
    ## ... Procrustes: rmse 0.000138889  max resid 0.0002526879 
    ## ... Similar to previous best
    ## Run 436 stress 8.959298e-05 
    ## ... Procrustes: rmse 0.0001686007  max resid 0.0003034297 
    ## ... Similar to previous best
    ## Run 437 stress 0.2848512 
    ## Run 438 stress 9.883236e-05 
    ## ... Procrustes: rmse 0.0001223556  max resid 0.0001674965 
    ## ... Similar to previous best
    ## Run 439 stress 9.374939e-05 
    ## ... Procrustes: rmse 0.0002249379  max resid 0.0004515821 
    ## ... Similar to previous best
    ## Run 440 stress 8.477108e-05 
    ## ... Procrustes: rmse 0.0001341875  max resid 0.000248874 
    ## ... Similar to previous best
    ## Run 441 stress 9.605446e-05 
    ## ... Procrustes: rmse 0.0001285134  max resid 0.0001918315 
    ## ... Similar to previous best
    ## Run 442 stress 9.567797e-05 
    ## ... Procrustes: rmse 0.0002261742  max resid 0.0004562135 
    ## ... Similar to previous best
    ## Run 443 stress 9.660669e-05 
    ## ... Procrustes: rmse 0.000150024  max resid 0.0002791446 
    ## ... Similar to previous best
    ## Run 444 stress 8.908613e-05 
    ## ... Procrustes: rmse 0.0001282079  max resid 0.0001836609 
    ## ... Similar to previous best
    ## Run 445 stress 9.43212e-05 
    ## ... Procrustes: rmse 0.0002268275  max resid 0.0004525532 
    ## ... Similar to previous best
    ## Run 446 stress 9.52721e-05 
    ## ... Procrustes: rmse 8.968253e-05  max resid 0.0001214884 
    ## ... Similar to previous best
    ## Run 447 stress 9.608803e-05 
    ## ... Procrustes: rmse 0.0002395728  max resid 0.0004640543 
    ## ... Similar to previous best
    ## Run 448 stress 9.791464e-05 
    ## ... Procrustes: rmse 0.0001980781  max resid 0.0004027811 
    ## ... Similar to previous best
    ## Run 449 stress 9.400025e-05 
    ## ... Procrustes: rmse 0.0002265708  max resid 0.0004287612 
    ## ... Similar to previous best
    ## Run 450 stress 9.140742e-05 
    ## ... Procrustes: rmse 0.0001437028  max resid 0.0002657854 
    ## ... Similar to previous best
    ## Run 451 stress 9.643188e-05 
    ## ... Procrustes: rmse 0.0001753953  max resid 0.0002437008 
    ## ... Similar to previous best
    ## Run 452 stress 9.528664e-05 
    ## ... Procrustes: rmse 0.0002343762  max resid 0.0004596423 
    ## ... Similar to previous best
    ## Run 453 stress 9.789813e-05 
    ## ... Procrustes: rmse 0.0001528846  max resid 0.0002928418 
    ## ... Similar to previous best
    ## Run 454 stress 9.12461e-05 
    ## ... Procrustes: rmse 0.0002273164  max resid 0.000445989 
    ## ... Similar to previous best
    ## Run 455 stress 8.96662e-05 
    ## ... Procrustes: rmse 0.0001107094  max resid 0.0001639784 
    ## ... Similar to previous best
    ## Run 456 stress 9.168463e-05 
    ## ... Procrustes: rmse 0.0002206293  max resid 0.0004442576 
    ## ... Similar to previous best
    ## Run 457 stress 7.677613e-05 
    ## ... Procrustes: rmse 0.0001528201  max resid 0.0003149453 
    ## ... Similar to previous best
    ## Run 458 stress 9.597177e-05 
    ## ... Procrustes: rmse 0.0001710127  max resid 0.0002431049 
    ## ... Similar to previous best
    ## Run 459 stress 9.390084e-05 
    ## ... Procrustes: rmse 0.0002285846  max resid 0.0004523866 
    ## ... Similar to previous best
    ## Run 460 stress 8.859912e-05 
    ## ... Procrustes: rmse 0.0002025876  max resid 0.0004042667 
    ## ... Similar to previous best
    ## Run 461 stress 9.187714e-05 
    ## ... Procrustes: rmse 0.0002199589  max resid 0.0004419113 
    ## ... Similar to previous best
    ## Run 462 stress 9.730406e-05 
    ## ... Procrustes: rmse 0.0002429143  max resid 0.0004598979 
    ## ... Similar to previous best
    ## Run 463 stress 9.95933e-05 
    ## ... Procrustes: rmse 0.0001710157  max resid 0.0002896218 
    ## ... Similar to previous best
    ## Run 464 stress 9.833634e-05 
    ## ... Procrustes: rmse 0.0002456594  max resid 0.000464936 
    ## ... Similar to previous best
    ## Run 465 stress 9.564833e-05 
    ## ... Procrustes: rmse 0.0002346114  max resid 0.0004478044 
    ## ... Similar to previous best
    ## Run 466 stress 9.676371e-05 
    ## ... Procrustes: rmse 0.0002379869  max resid 0.0004494548 
    ## ... Similar to previous best
    ## Run 467 stress 9.697679e-05 
    ## ... Procrustes: rmse 0.0001623708  max resid 0.0002366417 
    ## ... Similar to previous best
    ## Run 468 stress 9.187724e-05 
    ## ... Procrustes: rmse 0.0001343306  max resid 0.0002241144 
    ## ... Similar to previous best
    ## Run 469 stress 9.23777e-05 
    ## ... Procrustes: rmse 0.0002296375  max resid 0.000439294 
    ## ... Similar to previous best
    ## Run 470 stress 9.397916e-05 
    ## ... Procrustes: rmse 0.000224349  max resid 0.0004466522 
    ## ... Similar to previous best
    ## Run 471 stress 9.803283e-05 
    ## ... Procrustes: rmse 0.0002346629  max resid 0.0004657659 
    ## ... Similar to previous best
    ## Run 472 stress 9.374462e-05 
    ## ... Procrustes: rmse 0.0001414453  max resid 0.000252648 
    ## ... Similar to previous best
    ## Run 473 stress 9.266312e-05 
    ## ... Procrustes: rmse 0.0001678172  max resid 0.0002326627 
    ## ... Similar to previous best
    ## Run 474 stress 8.744475e-05 
    ## ... Procrustes: rmse 0.000113644  max resid 0.0001651541 
    ## ... Similar to previous best
    ## Run 475 stress 9.765904e-05 
    ## ... Procrustes: rmse 0.0001765666  max resid 0.0002454706 
    ## ... Similar to previous best
    ## Run 476 stress 9.227111e-05 
    ## ... Procrustes: rmse 0.0001662652  max resid 0.0002355392 
    ## ... Similar to previous best
    ## Run 477 stress 8.737113e-05 
    ## ... Procrustes: rmse 0.0001136588  max resid 0.0002117931 
    ## ... Similar to previous best
    ## Run 478 stress 9.494995e-05 
    ## ... Procrustes: rmse 0.0001479794  max resid 0.0002170746 
    ## ... Similar to previous best
    ## Run 479 stress 9.357399e-05 
    ## ... Procrustes: rmse 0.0002249886  max resid 0.0004519421 
    ## ... Similar to previous best
    ## Run 480 stress 9.210339e-05 
    ## ... Procrustes: rmse 0.0001678886  max resid 0.0003274517 
    ## ... Similar to previous best
    ## Run 481 stress 9.229466e-05 
    ## ... Procrustes: rmse 0.0001673938  max resid 0.000232421 
    ## ... Similar to previous best
    ## Run 482 stress 8.790394e-05 
    ## ... Procrustes: rmse 0.0001612715  max resid 0.0002970989 
    ## ... Similar to previous best
    ## Run 483 stress 7.844648e-05 
    ## ... Procrustes: rmse 0.000132619  max resid 0.0002181775 
    ## ... Similar to previous best
    ## Run 484 stress 9.361311e-05 
    ## ... Procrustes: rmse 0.0001227552  max resid 0.0001732964 
    ## ... Similar to previous best
    ## Run 485 stress 8.69345e-05 
    ## ... Procrustes: rmse 0.0001328369  max resid 0.0002407159 
    ## ... Similar to previous best
    ## Run 486 stress 9.892271e-05 
    ## ... Procrustes: rmse 0.0001575671  max resid 0.0002778138 
    ## ... Similar to previous best
    ## Run 487 stress 9.600757e-05 
    ## ... Procrustes: rmse 0.0002251755  max resid 0.00045239 
    ## ... Similar to previous best
    ## Run 488 stress 8.1476e-05 
    ## ... Procrustes: rmse 0.0001286824  max resid 0.0002055054 
    ## ... Similar to previous best
    ## Run 489 stress 8.681114e-05 
    ## ... Procrustes: rmse 0.0001470947  max resid 0.0002554189 
    ## ... Similar to previous best
    ## Run 490 stress 9.00238e-05 
    ## ... Procrustes: rmse 0.0001778229  max resid 0.0003102603 
    ## ... Similar to previous best
    ## Run 491 stress 9.955425e-05 
    ## ... Procrustes: rmse 0.000179146  max resid 0.0002545276 
    ## ... Similar to previous best
    ## Run 492 stress 8.762829e-05 
    ## ... Procrustes: rmse 0.000201037  max resid 0.000408 
    ## ... Similar to previous best
    ## Run 493 stress 9.100849e-05 
    ## ... Procrustes: rmse 0.0002166266  max resid 0.0004355534 
    ## ... Similar to previous best
    ## Run 494 stress 8.864114e-05 
    ## ... Procrustes: rmse 0.0001654528  max resid 0.0002928994 
    ## ... Similar to previous best
    ## Run 495 stress 9.965689e-05 
    ## ... Procrustes: rmse 0.0001207044  max resid 0.0001749905 
    ## ... Similar to previous best
    ## Run 496 stress 9.39902e-05 
    ## ... Procrustes: rmse 0.0002264768  max resid 0.0004543918 
    ## ... Similar to previous best
    ## Run 497 stress 8.813126e-05 
    ## ... Procrustes: rmse 0.0002074865  max resid 0.0004190101 
    ## ... Similar to previous best
    ## Run 498 stress 9.153186e-05 
    ## ... Procrustes: rmse 0.0001441795  max resid 0.0002626088 
    ## ... Similar to previous best
    ## Run 499 stress 9.6982e-05 
    ## ... Procrustes: rmse 0.0001574458  max resid 0.0002657335 
    ## ... Similar to previous best
    ## Run 500 stress 9.304373e-05 
    ## ... Procrustes: rmse 0.000179178  max resid 0.0003420575 
    ## ... Similar to previous best
    ## *** Best solution repeated 409 times

    ## Warning in metaMDS(PD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
PD_beta_geo_NMDS <- metaMDS(PD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07927535 
    ## Run 1 stress 0.1088749 
    ## Run 2 stress 0.08985618 
    ## Run 3 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732693  max resid 0.05627204 
    ## Run 4 stress 0.08985605 
    ## Run 5 stress 0.08985604 
    ## Run 6 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735029  max resid 0.0564776 
    ## Run 7 stress 0.08985603 
    ## Run 8 stress 0.09565364 
    ## Run 9 stress 0.1038073 
    ## Run 10 stress 0.07936952 
    ## ... Procrustes: rmse 0.009032617  max resid 0.03239468 
    ## Run 11 stress 0.09106043 
    ## Run 12 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173392  max resid 0.05638186 
    ## Run 13 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047487  max resid 0.03249374 
    ## Run 14 stress 0.08985603 
    ## Run 15 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238657  max resid 0.04981898 
    ## Run 16 stress 0.09252216 
    ## Run 17 stress 0.1058826 
    ## Run 18 stress 0.07936949 
    ## ... Procrustes: rmse 0.009043462  max resid 0.03246848 
    ## Run 19 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 3.503876e-05  max resid 9.220007e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.3655835 
    ## Run 21 stress 0.07927534 
    ## ... Procrustes: rmse 3.794418e-06  max resid 9.571286e-06 
    ## ... Similar to previous best
    ## Run 22 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237104  max resid 0.04966132 
    ## Run 23 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735597  max resid 0.05644165 
    ## Run 24 stress 0.07927535 
    ## ... Procrustes: rmse 4.822771e-05  max resid 0.0001286142 
    ## ... Similar to previous best
    ## Run 25 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237438  max resid 0.0496661 
    ## Run 26 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734728  max resid 0.05636916 
    ## Run 27 stress 0.07969263 
    ## ... Procrustes: rmse 0.01236982  max resid 0.04958003 
    ## Run 28 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733172  max resid 0.05626043 
    ## Run 29 stress 0.08985606 
    ## Run 30 stress 0.09544357 
    ## Run 31 stress 0.07936952 
    ## ... Procrustes: rmse 0.009041336  max resid 0.03244752 
    ## Run 32 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735099  max resid 0.05639522 
    ## Run 33 stress 0.07927535 
    ## ... Procrustes: rmse 4.588014e-05  max resid 0.0001147966 
    ## ... Similar to previous best
    ## Run 34 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735218  max resid 0.05641467 
    ## Run 35 stress 0.07969261 
    ## ... Procrustes: rmse 0.01238155  max resid 0.04975096 
    ## Run 36 stress 0.07969262 
    ## ... Procrustes: rmse 0.01240627  max resid 0.04988272 
    ## Run 37 stress 0.09296238 
    ## Run 38 stress 0.08985621 
    ## Run 39 stress 0.09106032 
    ## Run 40 stress 0.09043739 
    ## Run 41 stress 0.07951443 
    ## ... Procrustes: rmse 0.01733419  max resid 0.05627556 
    ## Run 42 stress 0.07936953 
    ## ... Procrustes: rmse 0.009032398  max resid 0.03238118 
    ## Run 43 stress 0.1090097 
    ## Run 44 stress 0.08985606 
    ## Run 45 stress 0.09356334 
    ## Run 46 stress 0.07936954 
    ## ... Procrustes: rmse 0.009030708  max resid 0.03237159 
    ## Run 47 stress 0.07969259 
    ## ... Procrustes: rmse 0.01236559  max resid 0.04962008 
    ## Run 48 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733311  max resid 0.05627548 
    ## Run 49 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047398  max resid 0.03248204 
    ## Run 50 stress 0.0904377 
    ## Run 51 stress 0.09388356 
    ## Run 52 stress 0.1047993 
    ## Run 53 stress 0.0796926 
    ## ... Procrustes: rmse 0.0123728  max resid 0.04963377 
    ## Run 54 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735191  max resid 0.05643485 
    ## Run 55 stress 0.07927534 
    ## ... Procrustes: rmse 8.36305e-06  max resid 1.704296e-05 
    ## ... Similar to previous best
    ## Run 56 stress 0.09022785 
    ## Run 57 stress 0.08985612 
    ## Run 58 stress 0.1125351 
    ## Run 59 stress 0.1047632 
    ## Run 60 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735084  max resid 0.05640359 
    ## Run 61 stress 0.07927535 
    ## ... Procrustes: rmse 3.139503e-05  max resid 7.111954e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.09388327 
    ## Run 63 stress 0.07936951 
    ## ... Procrustes: rmse 0.009035336  max resid 0.03239787 
    ## Run 64 stress 0.1041925 
    ## Run 65 stress 0.07927538 
    ## ... Procrustes: rmse 6.066426e-05  max resid 0.0001777114 
    ## ... Similar to previous best
    ## Run 66 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237975  max resid 0.0497046 
    ## Run 67 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049038  max resid 0.03246916 
    ## Run 68 stress 0.07927535 
    ## ... Procrustes: rmse 2.751569e-05  max resid 7.224805e-05 
    ## ... Similar to previous best
    ## Run 69 stress 0.07936951 
    ## ... Procrustes: rmse 0.009036338  max resid 0.03240555 
    ## Run 70 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237652  max resid 0.04969233 
    ## Run 71 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237513  max resid 0.04966347 
    ## Run 72 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237759  max resid 0.04970767 
    ## Run 73 stress 0.09072893 
    ## Run 74 stress 0.08985609 
    ## Run 75 stress 0.09252216 
    ## Run 76 stress 0.08985604 
    ## Run 77 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045191  max resid 0.03246268 
    ## Run 78 stress 0.07927534 
    ## ... Procrustes: rmse 8.179739e-06  max resid 2.344242e-05 
    ## ... Similar to previous best
    ## Run 79 stress 0.09172785 
    ## Run 80 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041948  max resid 0.03243901 
    ## Run 81 stress 0.07927538 
    ## ... Procrustes: rmse 7.350694e-05  max resid 0.0002073887 
    ## ... Similar to previous best
    ## Run 82 stress 0.09106048 
    ## Run 83 stress 0.09043758 
    ## Run 84 stress 0.07951444 
    ## ... Procrustes: rmse 0.01735103  max resid 0.05639566 
    ## Run 85 stress 0.0793695 
    ## ... Procrustes: rmse 0.009039672  max resid 0.03242615 
    ## Run 86 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238213  max resid 0.04968336 
    ## Run 87 stress 0.09247106 
    ## Run 88 stress 0.09043749 
    ## Run 89 stress 0.09485733 
    ## Run 90 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046218  max resid 0.03247099 
    ## Run 91 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733707  max resid 0.0563041 
    ## Run 92 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238595  max resid 0.04974916 
    ## Run 93 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735165  max resid 0.0563708 
    ## Run 94 stress 0.08985605 
    ## Run 95 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734946  max resid 0.05638392 
    ## Run 96 stress 0.09252221 
    ## Run 97 stress 0.09565362 
    ## Run 98 stress 0.09106032 
    ## Run 99 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734287  max resid 0.05634075 
    ## Run 100 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047744  max resid 0.03247752 
    ## Run 101 stress 0.0917278 
    ## Run 102 stress 0.09072897 
    ## Run 103 stress 0.09022787 
    ## Run 104 stress 0.08985605 
    ## Run 105 stress 0.07936951 
    ## ... Procrustes: rmse 0.009045329  max resid 0.03245821 
    ## Run 106 stress 0.0793695 
    ## ... Procrustes: rmse 0.009051679  max resid 0.03249903 
    ## Run 107 stress 0.0793695 
    ## ... Procrustes: rmse 0.009047331  max resid 0.03247316 
    ## Run 108 stress 0.08985605 
    ## Run 109 stress 0.09247112 
    ## Run 110 stress 0.09072893 
    ## Run 111 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237619  max resid 0.04971129 
    ## Run 112 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237973  max resid 0.04972059 
    ## Run 113 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735593  max resid 0.05642591 
    ## Run 114 stress 0.09296248 
    ## Run 115 stress 0.1046274 
    ## Run 116 stress 0.07927537 
    ## ... Procrustes: rmse 6.82883e-05  max resid 0.0001457889 
    ## ... Similar to previous best
    ## Run 117 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735287  max resid 0.05640576 
    ## Run 118 stress 0.07927536 
    ## ... Procrustes: rmse 3.557938e-05  max resid 8.87268e-05 
    ## ... Similar to previous best
    ## Run 119 stress 0.09388336 
    ## Run 120 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237759  max resid 0.04969803 
    ## Run 121 stress 0.07927537 
    ## ... Procrustes: rmse 7.320809e-05  max resid 0.0001780422 
    ## ... Similar to previous best
    ## Run 122 stress 0.07969261 
    ## ... Procrustes: rmse 0.01240247  max resid 0.04980935 
    ## Run 123 stress 0.07969265 
    ## ... Procrustes: rmse 0.01242748  max resid 0.04996307 
    ## Run 124 stress 0.07927535 
    ## ... Procrustes: rmse 2.447906e-05  max resid 6.253651e-05 
    ## ... Similar to previous best
    ## Run 125 stress 0.09022776 
    ## Run 126 stress 0.09565372 
    ## Run 127 stress 0.09382056 
    ## Run 128 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237572  max resid 0.04968362 
    ## Run 129 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237908  max resid 0.04970785 
    ## Run 130 stress 0.07927534 
    ## ... Procrustes: rmse 6.768785e-06  max resid 2.447225e-05 
    ## ... Similar to previous best
    ## Run 131 stress 0.106141 
    ## Run 132 stress 0.1033494 
    ## Run 133 stress 0.09072898 
    ## Run 134 stress 0.07969271 
    ## ... Procrustes: rmse 0.01236678  max resid 0.04953324 
    ## Run 135 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049531  max resid 0.03248875 
    ## Run 136 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735661  max resid 0.05641554 
    ## Run 137 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904335  max resid 0.03245133 
    ## Run 138 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237219  max resid 0.0496645 
    ## Run 139 stress 0.07927534 
    ## ... Procrustes: rmse 7.571094e-06  max resid 2.15252e-05 
    ## ... Similar to previous best
    ## Run 140 stress 0.08985608 
    ## Run 141 stress 0.07927535 
    ## ... Procrustes: rmse 2.454212e-05  max resid 6.47562e-05 
    ## ... Similar to previous best
    ## Run 142 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735338  max resid 0.05642647 
    ## Run 143 stress 0.07927534 
    ## ... Procrustes: rmse 1.23124e-05  max resid 2.943712e-05 
    ## ... Similar to previous best
    ## Run 144 stress 0.1084765 
    ## Run 145 stress 0.09106036 
    ## Run 146 stress 0.07927534 
    ## ... Procrustes: rmse 1.19016e-05  max resid 3.30432e-05 
    ## ... Similar to previous best
    ## Run 147 stress 0.07927536 
    ## ... Procrustes: rmse 3.863542e-05  max resid 0.0001054128 
    ## ... Similar to previous best
    ## Run 148 stress 0.07927536 
    ## ... Procrustes: rmse 1.599012e-05  max resid 4.035223e-05 
    ## ... Similar to previous best
    ## Run 149 stress 0.09382072 
    ## Run 150 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734984  max resid 0.05639218 
    ## Run 151 stress 0.09565379 
    ## Run 152 stress 0.09485729 
    ## Run 153 stress 0.09388333 
    ## Run 154 stress 0.09565383 
    ## Run 155 stress 0.07936951 
    ## ... Procrustes: rmse 0.009052654  max resid 0.03250057 
    ## Run 156 stress 0.07969268 
    ## ... Procrustes: rmse 0.01237218  max resid 0.04956454 
    ## Run 157 stress 0.07951448 
    ## ... Procrustes: rmse 0.01731737  max resid 0.05615401 
    ## Run 158 stress 0.0796926 
    ## ... Procrustes: rmse 0.0123807  max resid 0.04966049 
    ## Run 159 stress 0.0904374 
    ## Run 160 stress 0.08985617 
    ## Run 161 stress 0.08985603 
    ## Run 162 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733529  max resid 0.05628866 
    ## Run 163 stress 0.07969263 
    ## ... Procrustes: rmse 0.01239738  max resid 0.04985134 
    ## Run 164 stress 0.1089997 
    ## Run 165 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734434  max resid 0.05635141 
    ## Run 166 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237729  max resid 0.0497174 
    ## Run 167 stress 0.08985605 
    ## Run 168 stress 0.07951443 
    ## ... Procrustes: rmse 0.01733168  max resid 0.05626108 
    ## Run 169 stress 0.09290698 
    ## Run 170 stress 0.08985607 
    ## Run 171 stress 0.07969262 
    ## ... Procrustes: rmse 0.01238232  max resid 0.04964978 
    ## Run 172 stress 0.07927536 
    ## ... Procrustes: rmse 4.065309e-05  max resid 0.0001090768 
    ## ... Similar to previous best
    ## Run 173 stress 0.09356327 
    ## Run 174 stress 0.09296194 
    ## Run 175 stress 0.09043757 
    ## Run 176 stress 0.07969262 
    ## ... Procrustes: rmse 0.01238739  max resid 0.04978618 
    ## Run 177 stress 0.09072899 
    ## Run 178 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046859  max resid 0.03247654 
    ## Run 179 stress 0.0910603 
    ## Run 180 stress 0.0793695 
    ## ... Procrustes: rmse 0.009045069  max resid 0.03246266 
    ## Run 181 stress 0.09072893 
    ## Run 182 stress 0.0793695 
    ## ... Procrustes: rmse 0.009042895  max resid 0.03244874 
    ## Run 183 stress 0.08985614 
    ## Run 184 stress 0.09544368 
    ## Run 185 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237104  max resid 0.04965255 
    ## Run 186 stress 0.09388331 
    ## Run 187 stress 0.08985612 
    ## Run 188 stress 0.09290716 
    ## Run 189 stress 0.09474016 
    ## Run 190 stress 0.09072893 
    ## Run 191 stress 0.07927534 
    ## ... Procrustes: rmse 8.557409e-06  max resid 2.105442e-05 
    ## ... Similar to previous best
    ## Run 192 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735049  max resid 0.05637721 
    ## Run 193 stress 0.07951445 
    ## ... Procrustes: rmse 0.01738604  max resid 0.05662956 
    ## Run 194 stress 0.09072914 
    ## Run 195 stress 0.09106049 
    ## Run 196 stress 0.09072901 
    ## Run 197 stress 0.09106031 
    ## Run 198 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735251  max resid 0.05640382 
    ## Run 199 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237292  max resid 0.0496712 
    ## Run 200 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044076  max resid 0.0324556 
    ## Run 201 stress 0.0795144 
    ## ... Procrustes: rmse 0.0173461  max resid 0.0563694 
    ## Run 202 stress 0.0793695 
    ## ... Procrustes: rmse 0.009056017  max resid 0.03252517 
    ## Run 203 stress 0.090729 
    ## Run 204 stress 0.07927536 
    ## ... Procrustes: rmse 4.798297e-05  max resid 0.0001225027 
    ## ... Similar to previous best
    ## Run 205 stress 0.08985604 
    ## Run 206 stress 0.07951443 
    ## ... Procrustes: rmse 0.01736231  max resid 0.05653369 
    ## Run 207 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047446  max resid 0.03247819 
    ## Run 208 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735621  max resid 0.05642715 
    ## Run 209 stress 0.07927534 
    ## ... New best solution
    ## ... Procrustes: rmse 5.555316e-06  max resid 1.546982e-05 
    ## ... Similar to previous best
    ## Run 210 stress 0.09106031 
    ## Run 211 stress 0.08985615 
    ## Run 212 stress 0.1038073 
    ## Run 213 stress 0.07969266 
    ## ... Procrustes: rmse 0.01236654  max resid 0.04954745 
    ## Run 214 stress 0.07951441 
    ## ... Procrustes: rmse 0.01733353  max resid 0.05628569 
    ## Run 215 stress 0.08985612 
    ## Run 216 stress 0.09290693 
    ## Run 217 stress 0.07936954 
    ## ... Procrustes: rmse 0.009027267  max resid 0.03235326 
    ## Run 218 stress 0.09290713 
    ## Run 219 stress 0.08985604 
    ## Run 220 stress 0.07927535 
    ## ... Procrustes: rmse 2.661235e-05  max resid 6.491305e-05 
    ## ... Similar to previous best
    ## Run 221 stress 0.105681 
    ## Run 222 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237401  max resid 0.0496806 
    ## Run 223 stress 0.0795144 
    ## ... Procrustes: rmse 0.017345  max resid 0.05636454 
    ## Run 224 stress 0.07927535 
    ## ... Procrustes: rmse 2.734691e-05  max resid 7.206204e-05 
    ## ... Similar to previous best
    ## Run 225 stress 0.07927536 
    ## ... Procrustes: rmse 4.835514e-05  max resid 0.0001255686 
    ## ... Similar to previous best
    ## Run 226 stress 0.09565349 
    ## Run 227 stress 0.07951448 
    ## ... Procrustes: rmse 0.01732299  max resid 0.05620086 
    ## Run 228 stress 0.07936951 
    ## ... Procrustes: rmse 0.00904958  max resid 0.03249065 
    ## Run 229 stress 0.104507 
    ## Run 230 stress 0.09290697 
    ## Run 231 stress 0.08985608 
    ## Run 232 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045297  max resid 0.03247154 
    ## Run 233 stress 0.07951443 
    ## ... Procrustes: rmse 0.0173602  max resid 0.05651295 
    ## Run 234 stress 0.07951448 
    ## ... Procrustes: rmse 0.01739308  max resid 0.05669116 
    ## Run 235 stress 0.09544368 
    ## Run 236 stress 0.09544375 
    ## Run 237 stress 0.09247112 
    ## Run 238 stress 0.07969261 
    ## ... Procrustes: rmse 0.01236362  max resid 0.04958054 
    ## Run 239 stress 0.07927534 
    ## ... Procrustes: rmse 1.854416e-05  max resid 4.948729e-05 
    ## ... Similar to previous best
    ## Run 240 stress 0.09043776 
    ## Run 241 stress 0.09072895 
    ## Run 242 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045894  max resid 0.03247137 
    ## Run 243 stress 0.07927535 
    ## ... Procrustes: rmse 3.668485e-05  max resid 8.845565e-05 
    ## ... Similar to previous best
    ## Run 244 stress 0.07927535 
    ## ... Procrustes: rmse 2.792822e-05  max resid 7.608622e-05 
    ## ... Similar to previous best
    ## Run 245 stress 0.07927535 
    ## ... Procrustes: rmse 3.320562e-05  max resid 9.117807e-05 
    ## ... Similar to previous best
    ## Run 246 stress 0.09172781 
    ## Run 247 stress 0.07927537 
    ## ... Procrustes: rmse 6.207658e-05  max resid 0.0001736232 
    ## ... Similar to previous best
    ## Run 248 stress 0.09247102 
    ## Run 249 stress 0.07927534 
    ## ... Procrustes: rmse 2.162452e-05  max resid 5.910458e-05 
    ## ... Similar to previous best
    ## Run 250 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734599  max resid 0.0563664 
    ## Run 251 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237027  max resid 0.04960366 
    ## Run 252 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735512  max resid 0.05643788 
    ## Run 253 stress 0.08985618 
    ## Run 254 stress 0.07927534 
    ## ... Procrustes: rmse 8.21551e-06  max resid 1.903606e-05 
    ## ... Similar to previous best
    ## Run 255 stress 0.08985607 
    ## Run 256 stress 0.1038074 
    ## Run 257 stress 0.07951442 
    ## ... Procrustes: rmse 0.01736243  max resid 0.05650716 
    ## Run 258 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173599  max resid 0.05649423 
    ## Run 259 stress 0.08985605 
    ## Run 260 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046609  max resid 0.03247681 
    ## Run 261 stress 0.09022775 
    ## Run 262 stress 0.08985603 
    ## Run 263 stress 0.07927534 
    ## ... Procrustes: rmse 4.953546e-06  max resid 1.185436e-05 
    ## ... Similar to previous best
    ## Run 264 stress 0.09106039 
    ## Run 265 stress 0.09106056 
    ## Run 266 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734488  max resid 0.05639821 
    ## Run 267 stress 0.09072898 
    ## Run 268 stress 0.09252228 
    ## Run 269 stress 0.0910603 
    ## Run 270 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049444  max resid 0.03248603 
    ## Run 271 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734102  max resid 0.05634539 
    ## Run 272 stress 0.09382092 
    ## Run 273 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044982  max resid 0.03247434 
    ## Run 274 stress 0.07969266 
    ## ... Procrustes: rmse 0.01242308  max resid 0.04999083 
    ## Run 275 stress 0.09344108 
    ## Run 276 stress 0.07927534 
    ## ... Procrustes: rmse 1.724136e-05  max resid 4.434072e-05 
    ## ... Similar to previous best
    ## Run 277 stress 0.09072896 
    ## Run 278 stress 0.07927536 
    ## ... Procrustes: rmse 5.919096e-05  max resid 0.0001404016 
    ## ... Similar to previous best
    ## Run 279 stress 0.07951445 
    ## ... Procrustes: rmse 0.01738632  max resid 0.05663923 
    ## Run 280 stress 0.07936951 
    ## ... Procrustes: rmse 0.009036439  max resid 0.0324104 
    ## Run 281 stress 0.09106027 
    ## Run 282 stress 0.07936951 
    ## ... Procrustes: rmse 0.009051015  max resid 0.03249006 
    ## Run 283 stress 0.09485675 
    ## Run 284 stress 0.07927534 
    ## ... Procrustes: rmse 1.67063e-05  max resid 4.279848e-05 
    ## ... Similar to previous best
    ## Run 285 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043628  max resid 0.03245655 
    ## Run 286 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734279  max resid 0.05634932 
    ## Run 287 stress 0.07969259 
    ## ... Procrustes: rmse 0.01239286  max resid 0.04978076 
    ## Run 288 stress 0.07969261 
    ## ... Procrustes: rmse 0.01237315  max resid 0.04962805 
    ## Run 289 stress 0.09072894 
    ## Run 290 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735769  max resid 0.05649967 
    ## Run 291 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237668  max resid 0.04970592 
    ## Run 292 stress 0.07927535 
    ## ... Procrustes: rmse 2.535156e-05  max resid 7.068281e-05 
    ## ... Similar to previous best
    ## Run 293 stress 0.09565374 
    ## Run 294 stress 0.07969262 
    ## ... Procrustes: rmse 0.01239894  max resid 0.04986087 
    ## Run 295 stress 0.07936949 
    ## ... Procrustes: rmse 0.009048584  max resid 0.03248791 
    ## Run 296 stress 0.09022804 
    ## Run 297 stress 0.09382059 
    ## Run 298 stress 0.09252227 
    ## Run 299 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732788  max resid 0.05624734 
    ## Run 300 stress 0.09474082 
    ## Run 301 stress 0.09544389 
    ## Run 302 stress 0.07969263 
    ## ... Procrustes: rmse 0.01237224  max resid 0.04967024 
    ## Run 303 stress 0.09565356 
    ## Run 304 stress 0.07927534 
    ## ... Procrustes: rmse 1.006048e-05  max resid 2.799268e-05 
    ## ... Similar to previous best
    ## Run 305 stress 0.09290706 
    ## Run 306 stress 0.07951442 
    ## ... Procrustes: rmse 0.01737032  max resid 0.05651579 
    ## Run 307 stress 0.07951441 
    ## ... Procrustes: rmse 0.01736363  max resid 0.05648437 
    ## Run 308 stress 0.07951445 
    ## ... Procrustes: rmse 0.01732703  max resid 0.0562259 
    ## Run 309 stress 0.07969263 
    ## ... Procrustes: rmse 0.01236783  max resid 0.04958064 
    ## Run 310 stress 0.07927534 
    ## ... Procrustes: rmse 6.390047e-06  max resid 1.653676e-05 
    ## ... Similar to previous best
    ## Run 311 stress 0.07951447 
    ## ... Procrustes: rmse 0.01731864  max resid 0.05617048 
    ## Run 312 stress 0.09252219 
    ## Run 313 stress 0.07951446 
    ## ... Procrustes: rmse 0.0173224  max resid 0.0562315 
    ## Run 314 stress 0.09022777 
    ## Run 315 stress 0.0904374 
    ## Run 316 stress 0.09252216 
    ## Run 317 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049579  max resid 0.03249003 
    ## Run 318 stress 0.07936958 
    ## ... Procrustes: rmse 0.009019053  max resid 0.03230265 
    ## Run 319 stress 0.07936949 
    ## ... Procrustes: rmse 0.009046667  max resid 0.03247546 
    ## Run 320 stress 0.09106029 
    ## Run 321 stress 0.07951441 
    ## ... Procrustes: rmse 0.01734852  max resid 0.05638362 
    ## Run 322 stress 0.08985605 
    ## Run 323 stress 0.07969263 
    ## ... Procrustes: rmse 0.01237081  max resid 0.04959726 
    ## Run 324 stress 0.09106029 
    ## Run 325 stress 0.07969264 
    ## ... Procrustes: rmse 0.01237141  max resid 0.04959284 
    ## Run 326 stress 0.07927534 
    ## ... Procrustes: rmse 1.092018e-05  max resid 2.996391e-05 
    ## ... Similar to previous best
    ## Run 327 stress 0.09022784 
    ## Run 328 stress 0.0796926 
    ## ... Procrustes: rmse 0.01237151  max resid 0.04962694 
    ## Run 329 stress 0.07927536 
    ## ... Procrustes: rmse 5.363032e-05  max resid 0.0001244065 
    ## ... Similar to previous best
    ## Run 330 stress 0.0902279 
    ## Run 331 stress 0.07936952 
    ## ... Procrustes: rmse 0.009052085  max resid 0.03249099 
    ## Run 332 stress 0.0793695 
    ## ... Procrustes: rmse 0.00905014  max resid 0.03248887 
    ## Run 333 stress 0.0793695 
    ## ... Procrustes: rmse 0.00904756  max resid 0.03247885 
    ## Run 334 stress 0.07927535 
    ## ... Procrustes: rmse 3.804422e-05  max resid 0.0001060369 
    ## ... Similar to previous best
    ## Run 335 stress 0.090729 
    ## Run 336 stress 0.07936951 
    ## ... Procrustes: rmse 0.009037208  max resid 0.03241502 
    ## Run 337 stress 0.09022788 
    ## Run 338 stress 0.07951442 
    ## ... Procrustes: rmse 0.01736384  max resid 0.05652346 
    ## Run 339 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732852  max resid 0.05623653 
    ## Run 340 stress 0.08985603 
    ## Run 341 stress 0.09382074 
    ## Run 342 stress 0.07927538 
    ## ... Procrustes: rmse 2.489487e-05  max resid 4.682875e-05 
    ## ... Similar to previous best
    ## Run 343 stress 0.08985606 
    ## Run 344 stress 0.08985604 
    ## Run 345 stress 0.09106034 
    ## Run 346 stress 0.09252217 
    ## Run 347 stress 0.07927535 
    ## ... Procrustes: rmse 3.565403e-05  max resid 9.879245e-05 
    ## ... Similar to previous best
    ## Run 348 stress 0.0796926 
    ## ... Procrustes: rmse 0.0123801  max resid 0.0497545 
    ## Run 349 stress 0.09043754 
    ## Run 350 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049317  max resid 0.03249195 
    ## Run 351 stress 0.1033494 
    ## Run 352 stress 0.0793695 
    ## ... Procrustes: rmse 0.009043042  max resid 0.0324497 
    ## Run 353 stress 0.09106028 
    ## Run 354 stress 0.07969262 
    ## ... Procrustes: rmse 0.01237412  max resid 0.04960676 
    ## Run 355 stress 0.07927534 
    ## ... Procrustes: rmse 1.425062e-05  max resid 3.473471e-05 
    ## ... Similar to previous best
    ## Run 356 stress 0.07951442 
    ## ... Procrustes: rmse 0.01733981  max resid 0.05643851 
    ## Run 357 stress 0.07927538 
    ## ... Procrustes: rmse 6.956242e-05  max resid 0.0001966938 
    ## ... Similar to previous best
    ## Run 358 stress 0.0956537 
    ## Run 359 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734417  max resid 0.05635899 
    ## Run 360 stress 0.08985604 
    ## Run 361 stress 0.09382066 
    ## Run 362 stress 0.09382056 
    ## Run 363 stress 0.0929071 
    ## Run 364 stress 0.07951442 
    ## ... Procrustes: rmse 0.01737103  max resid 0.05652991 
    ## Run 365 stress 0.09106044 
    ## Run 366 stress 0.09296247 
    ## Run 367 stress 0.08985605 
    ## Run 368 stress 0.0796926 
    ## ... Procrustes: rmse 0.01238834  max resid 0.04978184 
    ## Run 369 stress 0.09072899 
    ## Run 370 stress 0.09565381 
    ## Run 371 stress 0.09388339 
    ## Run 372 stress 0.07927534 
    ## ... Procrustes: rmse 1.149571e-05  max resid 3.082163e-05 
    ## ... Similar to previous best
    ## Run 373 stress 0.08985611 
    ## Run 374 stress 0.09072892 
    ## Run 375 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237507  max resid 0.04967429 
    ## Run 376 stress 0.09474092 
    ## Run 377 stress 0.0938835 
    ## Run 378 stress 0.08985604 
    ## Run 379 stress 0.09388341 
    ## Run 380 stress 0.09290711 
    ## Run 381 stress 0.07969266 
    ## ... Procrustes: rmse 0.01239291  max resid 0.04985954 
    ## Run 382 stress 0.09344104 
    ## Run 383 stress 0.1108125 
    ## Run 384 stress 0.09043741 
    ## Run 385 stress 0.07927535 
    ## ... Procrustes: rmse 4.453754e-05  max resid 0.0001045329 
    ## ... Similar to previous best
    ## Run 386 stress 0.09247102 
    ## Run 387 stress 0.07936949 
    ## ... Procrustes: rmse 0.009045297  max resid 0.03247197 
    ## Run 388 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238684  max resid 0.04976447 
    ## Run 389 stress 0.1112806 
    ## Run 390 stress 0.07969258 
    ## ... Procrustes: rmse 0.01236915  max resid 0.04966071 
    ## Run 391 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735187  max resid 0.05641605 
    ## Run 392 stress 0.07936951 
    ## ... Procrustes: rmse 0.009045147  max resid 0.03246627 
    ## Run 393 stress 0.07927534 
    ## ... Procrustes: rmse 1.038894e-05  max resid 3.593861e-05 
    ## ... Similar to previous best
    ## Run 394 stress 0.07927534 
    ## ... Procrustes: rmse 2.372282e-05  max resid 5.205522e-05 
    ## ... Similar to previous best
    ## Run 395 stress 0.1067418 
    ## Run 396 stress 0.07927535 
    ## ... Procrustes: rmse 2.041867e-05  max resid 4.781746e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.0793695 
    ## ... Procrustes: rmse 0.009046262  max resid 0.03247587 
    ## Run 398 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238452  max resid 0.04975275 
    ## Run 399 stress 0.07969265 
    ## ... Procrustes: rmse 0.01237242  max resid 0.04958251 
    ## Run 400 stress 0.0956536 
    ## Run 401 stress 0.07927535 
    ## ... Procrustes: rmse 3.191841e-05  max resid 7.844547e-05 
    ## ... Similar to previous best
    ## Run 402 stress 0.08985606 
    ## Run 403 stress 0.08985605 
    ## Run 404 stress 0.09565368 
    ## Run 405 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237417  max resid 0.04970442 
    ## Run 406 stress 0.09290697 
    ## Run 407 stress 0.09565376 
    ## Run 408 stress 0.07969259 
    ## ... Procrustes: rmse 0.01237344  max resid 0.0496611 
    ## Run 409 stress 0.07927536 
    ## ... Procrustes: rmse 3.460838e-05  max resid 9.159175e-05 
    ## ... Similar to previous best
    ## Run 410 stress 0.09043736 
    ## Run 411 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735164  max resid 0.05638586 
    ## Run 412 stress 0.07927536 
    ## ... Procrustes: rmse 4.860192e-05  max resid 0.0001360104 
    ## ... Similar to previous best
    ## Run 413 stress 0.09072893 
    ## Run 414 stress 0.0793695 
    ## ... Procrustes: rmse 0.009040061  max resid 0.03243293 
    ## Run 415 stress 0.09043751 
    ## Run 416 stress 0.09290703 
    ## Run 417 stress 0.0954436 
    ## Run 418 stress 0.07936951 
    ## ... Procrustes: rmse 0.009036948  max resid 0.03241316 
    ## Run 419 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237824  max resid 0.04970052 
    ## Run 420 stress 0.09106037 
    ## Run 421 stress 0.09072893 
    ## Run 422 stress 0.08985613 
    ## Run 423 stress 0.09382057 
    ## Run 424 stress 0.09043761 
    ## Run 425 stress 0.07936953 
    ## ... Procrustes: rmse 0.009029125  max resid 0.03236465 
    ## Run 426 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238463  max resid 0.04976307 
    ## Run 427 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237059  max resid 0.04966628 
    ## Run 428 stress 0.09072898 
    ## Run 429 stress 0.07927534 
    ## ... Procrustes: rmse 5.890665e-06  max resid 1.552513e-05 
    ## ... Similar to previous best
    ## Run 430 stress 0.104507 
    ## Run 431 stress 0.1089997 
    ## Run 432 stress 0.09043763 
    ## Run 433 stress 0.07951441 
    ## ... Procrustes: rmse 0.0173361  max resid 0.05634224 
    ## Run 434 stress 0.0793695 
    ## ... Procrustes: rmse 0.0090435  max resid 0.03245563 
    ## Run 435 stress 0.09172791 
    ## Run 436 stress 0.07927534 
    ## ... Procrustes: rmse 7.311115e-06  max resid 1.527569e-05 
    ## ... Similar to previous best
    ## Run 437 stress 0.0793695 
    ## ... Procrustes: rmse 0.009049419  max resid 0.03249289 
    ## Run 438 stress 0.08985603 
    ## Run 439 stress 0.08985608 
    ## Run 440 stress 0.08985603 
    ## Run 441 stress 0.07927534 
    ## ... Procrustes: rmse 2.05549e-05  max resid 5.26943e-05 
    ## ... Similar to previous best
    ## Run 442 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735075  max resid 0.05640224 
    ## Run 443 stress 0.08985605 
    ## Run 444 stress 0.09022776 
    ## Run 445 stress 0.07951441 
    ## ... Procrustes: rmse 0.01735224  max resid 0.05639516 
    ## Run 446 stress 0.07927535 
    ## ... Procrustes: rmse 2.279733e-05  max resid 6.231376e-05 
    ## ... Similar to previous best
    ## Run 447 stress 0.09565352 
    ## Run 448 stress 0.07936959 
    ## ... Procrustes: rmse 0.009018024  max resid 0.03229646 
    ## Run 449 stress 0.07936951 
    ## ... Procrustes: rmse 0.009036954  max resid 0.03241168 
    ## Run 450 stress 0.09072902 
    ## Run 451 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044294  max resid 0.03246363 
    ## Run 452 stress 0.09043767 
    ## Run 453 stress 0.07936949 
    ## ... Procrustes: rmse 0.00904731  max resid 0.03248403 
    ## Run 454 stress 0.09474092 
    ## Run 455 stress 0.09474028 
    ## Run 456 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735271  max resid 0.05641881 
    ## Run 457 stress 0.1109853 
    ## Run 458 stress 0.09072898 
    ## Run 459 stress 0.07951444 
    ## ... Procrustes: rmse 0.01735751  max resid 0.05650426 
    ## Run 460 stress 0.0907291 
    ## Run 461 stress 0.09106046 
    ## Run 462 stress 0.0910603 
    ## Run 463 stress 0.07927538 
    ## ... Procrustes: rmse 2.145433e-05  max resid 6.200539e-05 
    ## ... Similar to previous best
    ## Run 464 stress 0.1047633 
    ## Run 465 stress 0.09106036 
    ## Run 466 stress 0.07969258 
    ## ... Procrustes: rmse 0.01237356  max resid 0.04967974 
    ## Run 467 stress 0.07927536 
    ## ... Procrustes: rmse 4.326205e-05  max resid 0.000119008 
    ## ... Similar to previous best
    ## Run 468 stress 0.07969259 
    ## ... Procrustes: rmse 0.01238156  max resid 0.04973314 
    ## Run 469 stress 0.1044027 
    ## Run 470 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044667  max resid 0.03246464 
    ## Run 471 stress 0.09382082 
    ## Run 472 stress 0.07951442 
    ## ... Procrustes: rmse 0.01735641  max resid 0.05643779 
    ## Run 473 stress 0.1038073 
    ## Run 474 stress 0.09106031 
    ## Run 475 stress 0.09106044 
    ## Run 476 stress 0.0793695 
    ## ... Procrustes: rmse 0.009041653  max resid 0.03244355 
    ## Run 477 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044599  max resid 0.0324637 
    ## Run 478 stress 0.07969258 
    ## ... Procrustes: rmse 0.0123731  max resid 0.0496751 
    ## Run 479 stress 0.09344096 
    ## Run 480 stress 0.09022791 
    ## Run 481 stress 0.07927535 
    ## ... Procrustes: rmse 3.134427e-05  max resid 8.577915e-05 
    ## ... Similar to previous best
    ## Run 482 stress 0.07969261 
    ## ... Procrustes: rmse 0.01238104  max resid 0.0497715 
    ## Run 483 stress 0.0795144 
    ## ... Procrustes: rmse 0.01735325  max resid 0.05641601 
    ## Run 484 stress 0.0793695 
    ## ... Procrustes: rmse 0.009044037  max resid 0.03245843 
    ## Run 485 stress 0.07927534 
    ## ... Procrustes: rmse 1.469562e-05  max resid 3.582926e-05 
    ## ... Similar to previous best
    ## Run 486 stress 0.09290693 
    ## Run 487 stress 0.0910605 
    ## Run 488 stress 0.0795144 
    ## ... Procrustes: rmse 0.01734972  max resid 0.05639538 
    ## Run 489 stress 0.08985605 
    ## Run 490 stress 0.07927535 
    ## ... Procrustes: rmse 1.682168e-05  max resid 4.279526e-05 
    ## ... Similar to previous best
    ## Run 491 stress 0.07927537 
    ## ... Procrustes: rmse 5.001595e-05  max resid 0.0001379018 
    ## ... Similar to previous best
    ## Run 492 stress 0.09474004 
    ## Run 493 stress 0.07951443 
    ## ... Procrustes: rmse 0.01732543  max resid 0.05636696 
    ## Run 494 stress 0.07927536 
    ## ... Procrustes: rmse 5.781322e-05  max resid 0.000137784 
    ## ... Similar to previous best
    ## Run 495 stress 0.07951444 
    ## ... Procrustes: rmse 0.01732419  max resid 0.05621551 
    ## Run 496 stress 0.07936949 
    ## ... Procrustes: rmse 0.009047836  max resid 0.03248405 
    ## Run 497 stress 0.07969261 
    ## ... Procrustes: rmse 0.0123994  max resid 0.04984132 
    ## Run 498 stress 0.07969266 
    ## ... Procrustes: rmse 0.0123714  max resid 0.04960961 
    ## Run 499 stress 0.09043743 
    ## Run 500 stress 0.09290723 
    ## *** Best solution repeated 44 times

``` r
# Mixed and stratified lakes
PD_beta_geo_MS_NMDS <- metaMDS(PD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.05036802 
    ## Run 1 stress 0.05900159 
    ## Run 2 stress 0.05051326 
    ## ... Procrustes: rmse 0.03578053  max resid 0.1222478 
    ## Run 3 stress 0.0596863 
    ## Run 4 stress 0.05051361 
    ## ... Procrustes: rmse 0.03580884  max resid 0.1223409 
    ## Run 5 stress 0.05829929 
    ## Run 6 stress 0.05337589 
    ## Run 7 stress 0.05337581 
    ## Run 8 stress 0.05337579 
    ## Run 9 stress 0.05051311 
    ## ... Procrustes: rmse 0.03574919  max resid 0.1220214 
    ## Run 10 stress 0.06938763 
    ## Run 11 stress 0.06309838 
    ## Run 12 stress 0.06478545 
    ## Run 13 stress 0.05928212 
    ## Run 14 stress 0.05051337 
    ## ... Procrustes: rmse 0.03579063  max resid 0.1222815 
    ## Run 15 stress 0.05171061 
    ## Run 16 stress 0.05403591 
    ## Run 17 stress 0.06730253 
    ## Run 18 stress 0.06609051 
    ## Run 19 stress 0.05968629 
    ## Run 20 stress 0.06038949 
    ## Run 21 stress 0.05036803 
    ## ... Procrustes: rmse 0.0003311936  max resid 0.000842938 
    ## ... Similar to previous best
    ## Run 22 stress 0.05036792 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001194093  max resid 0.0003070453 
    ## ... Similar to previous best
    ## Run 23 stress 0.0540358 
    ## Run 24 stress 0.05036794 
    ## ... Procrustes: rmse 0.0001230186  max resid 0.0002983243 
    ## ... Similar to previous best
    ## Run 25 stress 0.3592585 
    ## Run 26 stress 0.05051287 
    ## ... Procrustes: rmse 0.03559686  max resid 0.1216354 
    ## Run 27 stress 0.05171058 
    ## Run 28 stress 0.05036793 
    ## ... Procrustes: rmse 2.553388e-05  max resid 5.334723e-05 
    ## ... Similar to previous best
    ## Run 29 stress 0.05036803 
    ## ... Procrustes: rmse 0.000208167  max resid 0.0005220047 
    ## ... Similar to previous best
    ## Run 30 stress 0.06136682 
    ## Run 31 stress 0.06136675 
    ## Run 32 stress 0.05051332 
    ## ... Procrustes: rmse 0.0356919  max resid 0.1217992 
    ## Run 33 stress 0.05036799 
    ## ... Procrustes: rmse 7.88244e-05  max resid 0.0002260865 
    ## ... Similar to previous best
    ## Run 34 stress 0.05337587 
    ## Run 35 stress 0.05036793 
    ## ... Procrustes: rmse 0.0001163104  max resid 0.0002859512 
    ## ... Similar to previous best
    ## Run 36 stress 0.05171061 
    ## Run 37 stress 0.05602044 
    ## Run 38 stress 0.05602048 
    ## Run 39 stress 0.05829934 
    ## Run 40 stress 0.06231289 
    ## Run 41 stress 0.06478546 
    ## Run 42 stress 0.06111767 
    ## Run 43 stress 0.05051354 
    ## ... Procrustes: rmse 0.03575229  max resid 0.1221495 
    ## Run 44 stress 0.0611177 
    ## Run 45 stress 0.0689692 
    ## Run 46 stress 0.06111774 
    ## Run 47 stress 0.06640362 
    ## Run 48 stress 0.05403593 
    ## Run 49 stress 0.05036792 
    ## ... Procrustes: rmse 8.949302e-05  max resid 0.0002241695 
    ## ... Similar to previous best
    ## Run 50 stress 0.05928216 
    ## Run 51 stress 0.05968636 
    ## Run 52 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001529475  max resid 0.0003841767 
    ## ... Similar to previous best
    ## Run 53 stress 0.05968629 
    ## Run 54 stress 0.05466702 
    ## Run 55 stress 0.06499239 
    ## Run 56 stress 0.05036792 
    ## ... Procrustes: rmse 1.377913e-05  max resid 3.541149e-05 
    ## ... Similar to previous best
    ## Run 57 stress 0.05337578 
    ## Run 58 stress 0.06705055 
    ## Run 59 stress 0.05051338 
    ## ... Procrustes: rmse 0.03567559  max resid 0.1217459 
    ## Run 60 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001447221  max resid 0.0003498194 
    ## ... Similar to previous best
    ## Run 61 stress 0.05171058 
    ## Run 62 stress 0.06478524 
    ## Run 63 stress 0.06111774 
    ## Run 64 stress 0.05036793 
    ## ... Procrustes: rmse 0.0001052087  max resid 0.000266357 
    ## ... Similar to previous best
    ## Run 65 stress 0.05036803 
    ## ... Procrustes: rmse 0.0002018082  max resid 0.0005059888 
    ## ... Similar to previous best
    ## Run 66 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001926756  max resid 0.0004854003 
    ## ... Similar to previous best
    ## Run 67 stress 0.05051315 
    ## ... Procrustes: rmse 0.03572574  max resid 0.1220572 
    ## Run 68 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 5.358783e-05  max resid 0.0001398868 
    ## ... Similar to previous best
    ## Run 69 stress 0.06730237 
    ## Run 70 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001589964  max resid 0.0003646078 
    ## ... Similar to previous best
    ## Run 71 stress 0.05171058 
    ## Run 72 stress 0.05051292 
    ## ... Procrustes: rmse 0.03569722  max resid 0.1219484 
    ## Run 73 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001786222  max resid 0.0004788839 
    ## ... Similar to previous best
    ## Run 74 stress 0.05051335 
    ## ... Procrustes: rmse 0.03570135  max resid 0.1218127 
    ## Run 75 stress 0.0721793 
    ## Run 76 stress 0.05171062 
    ## Run 77 stress 0.05051333 
    ## ... Procrustes: rmse 0.03569989  max resid 0.1218108 
    ## Run 78 stress 0.05171061 
    ## Run 79 stress 0.0582993 
    ## Run 80 stress 0.0590015 
    ## Run 81 stress 0.05466706 
    ## Run 82 stress 0.05051308 
    ## ... Procrustes: rmse 0.03570128  max resid 0.1219714 
    ## Run 83 stress 0.05979229 
    ## Run 84 stress 0.06309843 
    ## Run 85 stress 0.05051292 
    ## ... Procrustes: rmse 0.03568993  max resid 0.121925 
    ## Run 86 stress 0.06478539 
    ## Run 87 stress 0.06478524 
    ## Run 88 stress 0.05051353 
    ## ... Procrustes: rmse 0.0356951  max resid 0.1219712 
    ## Run 89 stress 0.050368 
    ## ... Procrustes: rmse 0.0001561591  max resid 0.0004483987 
    ## ... Similar to previous best
    ## Run 90 stress 0.05829931 
    ## Run 91 stress 0.05051327 
    ## ... Procrustes: rmse 0.03568966  max resid 0.1219459 
    ## Run 92 stress 0.05036793 
    ## ... Procrustes: rmse 7.58072e-05  max resid 0.0001828573 
    ## ... Similar to previous best
    ## Run 93 stress 0.05051329 
    ## ... Procrustes: rmse 0.03571383  max resid 0.1220177 
    ## Run 94 stress 0.06136688 
    ## Run 95 stress 0.05036795 
    ## ... Procrustes: rmse 6.593747e-05  max resid 0.000156951 
    ## ... Similar to previous best
    ## Run 96 stress 0.05051347 
    ## ... Procrustes: rmse 0.03570548  max resid 0.1218133 
    ## Run 97 stress 0.05602057 
    ## Run 98 stress 0.064937 
    ## Run 99 stress 0.0505133 
    ## ... Procrustes: rmse 0.03571631  max resid 0.122025 
    ## Run 100 stress 0.05051339 
    ## ... Procrustes: rmse 0.03574892  max resid 0.1221261 
    ## Run 101 stress 0.05602044 
    ## Run 102 stress 0.05036794 
    ## ... Procrustes: rmse 6.23947e-05  max resid 0.0001415137 
    ## ... Similar to previous best
    ## Run 103 stress 0.05968638 
    ## Run 104 stress 0.05602048 
    ## Run 105 stress 0.06231287 
    ## Run 106 stress 0.05466709 
    ## Run 107 stress 0.0582993 
    ## Run 108 stress 0.05036793 
    ## ... Procrustes: rmse 7.162146e-05  max resid 0.0001763011 
    ## ... Similar to previous best
    ## Run 109 stress 0.06231287 
    ## Run 110 stress 0.06136681 
    ## Run 111 stress 0.05051329 
    ## ... Procrustes: rmse 0.03571774  max resid 0.1220303 
    ## Run 112 stress 0.05829929 
    ## Run 113 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001013608  max resid 0.0002474991 
    ## ... Similar to previous best
    ## Run 114 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001523104  max resid 0.0003796136 
    ## ... Similar to previous best
    ## Run 115 stress 0.05051359 
    ## ... Procrustes: rmse 0.03565505  max resid 0.1216579 
    ## Run 116 stress 0.05829935 
    ## Run 117 stress 0.05968632 
    ## Run 118 stress 0.05337577 
    ## Run 119 stress 0.05466703 
    ## Run 120 stress 0.05466705 
    ## Run 121 stress 0.05051299 
    ## ... Procrustes: rmse 0.03568226  max resid 0.1219091 
    ## Run 122 stress 0.06622045 
    ## Run 123 stress 0.06362796 
    ## Run 124 stress 0.06730256 
    ## Run 125 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001860072  max resid 0.0005010404 
    ## ... Similar to previous best
    ## Run 126 stress 0.05602051 
    ## Run 127 stress 0.05171058 
    ## Run 128 stress 0.05602046 
    ## Run 129 stress 0.05602059 
    ## Run 130 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001331261  max resid 0.0003286064 
    ## ... Similar to previous best
    ## Run 131 stress 0.05051302 
    ## ... Procrustes: rmse 0.03567105  max resid 0.1217558 
    ## Run 132 stress 0.06309862 
    ## Run 133 stress 0.05051303 
    ## ... Procrustes: rmse 0.03572063  max resid 0.1219055 
    ## Run 134 stress 0.05051298 
    ## ... Procrustes: rmse 0.03569016  max resid 0.1219315 
    ## Run 135 stress 0.06309838 
    ## Run 136 stress 0.05337586 
    ## Run 137 stress 0.05829937 
    ## Run 138 stress 0.05337578 
    ## Run 139 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001855427  max resid 0.0005013726 
    ## ... Similar to previous best
    ## Run 140 stress 0.06362795 
    ## Run 141 stress 0.050368 
    ## ... Procrustes: rmse 0.0001342533  max resid 0.0003442389 
    ## ... Similar to previous best
    ## Run 142 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001660182  max resid 0.000403775 
    ## ... Similar to previous best
    ## Run 143 stress 0.05829937 
    ## Run 144 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001458204  max resid 0.0003568798 
    ## ... Similar to previous best
    ## Run 145 stress 0.05051322 
    ## ... Procrustes: rmse 0.03573289  max resid 0.1219208 
    ## Run 146 stress 0.05171063 
    ## Run 147 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001205895  max resid 0.0002625246 
    ## ... Similar to previous best
    ## Run 148 stress 0.06730248 
    ## Run 149 stress 0.05968637 
    ## Run 150 stress 0.0582993 
    ## Run 151 stress 0.05171059 
    ## Run 152 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001350571  max resid 0.0003808112 
    ## ... Similar to previous best
    ## Run 153 stress 0.05337578 
    ## Run 154 stress 0.05466701 
    ## Run 155 stress 0.05928214 
    ## Run 156 stress 0.05051324 
    ## ... Procrustes: rmse 0.03571647  max resid 0.1220248 
    ## Run 157 stress 0.05051358 
    ## ... Procrustes: rmse 0.03569725  max resid 0.1219792 
    ## Run 158 stress 0.05466706 
    ## Run 159 stress 0.06231303 
    ## Run 160 stress 0.06647186 
    ## Run 161 stress 0.05602052 
    ## Run 162 stress 0.06136675 
    ## Run 163 stress 0.05171059 
    ## Run 164 stress 0.05337593 
    ## Run 165 stress 0.05051329 
    ## ... Procrustes: rmse 0.03574085  max resid 0.1220975 
    ## Run 166 stress 0.05051341 
    ## ... Procrustes: rmse 0.03569767  max resid 0.1217998 
    ## Run 167 stress 0.05337578 
    ## Run 168 stress 0.05466701 
    ## Run 169 stress 0.05051317 
    ## ... Procrustes: rmse 0.03566698  max resid 0.1217262 
    ## Run 170 stress 0.06476217 
    ## Run 171 stress 0.07242205 
    ## Run 172 stress 0.05051328 
    ## ... Procrustes: rmse 0.035728  max resid 0.1220603 
    ## Run 173 stress 0.06231316 
    ## Run 174 stress 0.05036792 
    ## ... Procrustes: rmse 4.268509e-05  max resid 9.975598e-05 
    ## ... Similar to previous best
    ## Run 175 stress 0.05403581 
    ## Run 176 stress 0.0533758 
    ## Run 177 stress 0.05979237 
    ## Run 178 stress 0.05036793 
    ## ... Procrustes: rmse 8.60551e-05  max resid 0.000229453 
    ## ... Similar to previous best
    ## Run 179 stress 0.05051355 
    ## ... Procrustes: rmse 0.03569846  max resid 0.1219818 
    ## Run 180 stress 0.05928249 
    ## Run 181 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001384477  max resid 0.0003406021 
    ## ... Similar to previous best
    ## Run 182 stress 0.05403588 
    ## Run 183 stress 0.0582994 
    ## Run 184 stress 0.05036792 
    ## ... Procrustes: rmse 5.808461e-05  max resid 0.0001432152 
    ## ... Similar to previous best
    ## Run 185 stress 0.05171058 
    ## Run 186 stress 0.05171058 
    ## Run 187 stress 0.05403615 
    ## Run 188 stress 0.05171058 
    ## Run 189 stress 0.05337583 
    ## Run 190 stress 0.0533758 
    ## Run 191 stress 0.05036796 
    ## ... Procrustes: rmse 0.0001217789  max resid 0.0003218936 
    ## ... Similar to previous best
    ## Run 192 stress 0.0517106 
    ## Run 193 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001142834  max resid 0.0002748854 
    ## ... Similar to previous best
    ## Run 194 stress 0.05171063 
    ## Run 195 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001270831  max resid 0.0003609959 
    ## ... Similar to previous best
    ## Run 196 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001493109  max resid 0.0003693154 
    ## ... Similar to previous best
    ## Run 197 stress 0.05036795 
    ## ... Procrustes: rmse 7.523185e-05  max resid 0.0001617793 
    ## ... Similar to previous best
    ## Run 198 stress 0.05051332 
    ## ... Procrustes: rmse 0.03569676  max resid 0.1218002 
    ## Run 199 stress 0.05171058 
    ## Run 200 stress 0.06476224 
    ## Run 201 stress 0.05051295 
    ## ... Procrustes: rmse 0.0357214  max resid 0.1219172 
    ## Run 202 stress 0.05171059 
    ## Run 203 stress 0.3451223 
    ## Run 204 stress 0.05051317 
    ## ... Procrustes: rmse 0.03571465  max resid 0.1220158 
    ## Run 205 stress 0.05051331 
    ## ... Procrustes: rmse 0.03569268  max resid 0.1217902 
    ## Run 206 stress 0.06789755 
    ## Run 207 stress 0.05403573 
    ## Run 208 stress 0.0505129 
    ## ... Procrustes: rmse 0.03567855  max resid 0.1217943 
    ## Run 209 stress 0.06478533 
    ## Run 210 stress 0.05466703 
    ## Run 211 stress 0.05051277 
    ## ... Procrustes: rmse 0.03559442  max resid 0.121594 
    ## Run 212 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001084341  max resid 0.0002653939 
    ## ... Similar to previous best
    ## Run 213 stress 0.05036797 
    ## ... Procrustes: rmse 9.815974e-05  max resid 0.0002357928 
    ## ... Similar to previous best
    ## Run 214 stress 0.05051346 
    ## ... Procrustes: rmse 0.03571851  max resid 0.122039 
    ## Run 215 stress 0.05829931 
    ## Run 216 stress 0.05466704 
    ## Run 217 stress 0.0517106 
    ## Run 218 stress 0.05829948 
    ## Run 219 stress 0.05602046 
    ## Run 220 stress 0.05171059 
    ## Run 221 stress 0.05036795 
    ## ... Procrustes: rmse 8.854747e-05  max resid 0.0002145445 
    ## ... Similar to previous best
    ## Run 222 stress 0.05051276 
    ## ... Procrustes: rmse 0.03561905  max resid 0.1216879 
    ## Run 223 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001289713  max resid 0.0003351157 
    ## ... Similar to previous best
    ## Run 224 stress 0.05602043 
    ## Run 225 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001572708  max resid 0.0004515387 
    ## ... Similar to previous best
    ## Run 226 stress 0.05403554 
    ## Run 227 stress 0.05036797 
    ## ... Procrustes: rmse 9.284041e-05  max resid 0.000228944 
    ## ... Similar to previous best
    ## Run 228 stress 0.05968625 
    ## Run 229 stress 0.06231291 
    ## Run 230 stress 0.05036791 
    ## ... Procrustes: rmse 5.004511e-05  max resid 0.0001298146 
    ## ... Similar to previous best
    ## Run 231 stress 0.05979235 
    ## Run 232 stress 0.06362783 
    ## Run 233 stress 0.05171058 
    ## Run 234 stress 0.05036797 
    ## ... Procrustes: rmse 0.00012518  max resid 0.0003633335 
    ## ... Similar to previous best
    ## Run 235 stress 0.05403554 
    ## Run 236 stress 0.05051326 
    ## ... Procrustes: rmse 0.03570499  max resid 0.121991 
    ## Run 237 stress 0.06309845 
    ## Run 238 stress 0.05968634 
    ## Run 239 stress 0.05829929 
    ## Run 240 stress 0.06231308 
    ## Run 241 stress 0.06231313 
    ## Run 242 stress 0.050368 
    ## ... Procrustes: rmse 0.0001627055  max resid 0.0004363138 
    ## ... Similar to previous best
    ## Run 243 stress 0.05036806 
    ## ... Procrustes: rmse 0.000176976  max resid 0.0004867605 
    ## ... Similar to previous best
    ## Run 244 stress 0.05036801 
    ## ... Procrustes: rmse 0.000141266  max resid 0.0003479966 
    ## ... Similar to previous best
    ## Run 245 stress 0.06366665 
    ## Run 246 stress 0.05968627 
    ## Run 247 stress 0.05036793 
    ## ... Procrustes: rmse 6.346518e-05  max resid 0.0001461476 
    ## ... Similar to previous best
    ## Run 248 stress 0.05036793 
    ## ... Procrustes: rmse 9.088164e-05  max resid 0.0002292265 
    ## ... Similar to previous best
    ## Run 249 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001081069  max resid 0.0002557277 
    ## ... Similar to previous best
    ## Run 250 stress 0.06621927 
    ## Run 251 stress 0.05051329 
    ## ... Procrustes: rmse 0.03569348  max resid 0.1219578 
    ## Run 252 stress 0.05337594 
    ## Run 253 stress 0.06630528 
    ## Run 254 stress 0.0517106 
    ## Run 255 stress 0.05829936 
    ## Run 256 stress 0.0517106 
    ## Run 257 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001111933  max resid 0.0002730463 
    ## ... Similar to previous best
    ## Run 258 stress 0.06038955 
    ## Run 259 stress 0.05036796 
    ## ... Procrustes: rmse 9.997291e-05  max resid 0.0002281777 
    ## ... Similar to previous best
    ## Run 260 stress 0.06362802 
    ## Run 261 stress 0.06231295 
    ## Run 262 stress 0.06136679 
    ## Run 263 stress 0.05337589 
    ## Run 264 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001198334  max resid 0.0003079095 
    ## ... Similar to previous best
    ## Run 265 stress 0.06309856 
    ## Run 266 stress 0.05337592 
    ## Run 267 stress 0.05602055 
    ## Run 268 stress 0.05036797 
    ## ... Procrustes: rmse 0.0001053644  max resid 0.0002484337 
    ## ... Similar to previous best
    ## Run 269 stress 0.06231297 
    ## Run 270 stress 0.05036804 
    ## ... Procrustes: rmse 9.028655e-05  max resid 0.0002149134 
    ## ... Similar to previous best
    ## Run 271 stress 0.0517106 
    ## Run 272 stress 0.05171059 
    ## Run 273 stress 0.05968633 
    ## Run 274 stress 0.05829929 
    ## Run 275 stress 0.05051325 
    ## ... Procrustes: rmse 0.03571423  max resid 0.1220186 
    ## Run 276 stress 0.05602069 
    ## Run 277 stress 0.05036794 
    ## ... Procrustes: rmse 8.583837e-05  max resid 0.0002213659 
    ## ... Similar to previous best
    ## Run 278 stress 0.05036792 
    ## ... Procrustes: rmse 5.350543e-05  max resid 0.0001314542 
    ## ... Similar to previous best
    ## Run 279 stress 0.05928255 
    ## Run 280 stress 0.05602063 
    ## Run 281 stress 0.05171059 
    ## Run 282 stress 0.06476221 
    ## Run 283 stress 0.05036794 
    ## ... Procrustes: rmse 9.907151e-05  max resid 0.0002815011 
    ## ... Similar to previous best
    ## Run 284 stress 0.05968636 
    ## Run 285 stress 0.05051325 
    ## ... Procrustes: rmse 0.03566531  max resid 0.1217262 
    ## Run 286 stress 0.05900149 
    ## Run 287 stress 0.0505133 
    ## ... Procrustes: rmse 0.03570212  max resid 0.1218189 
    ## Run 288 stress 0.0546671 
    ## Run 289 stress 0.05051358 
    ## ... Procrustes: rmse 0.03569703  max resid 0.121978 
    ## Run 290 stress 0.05036791 
    ## ... Procrustes: rmse 2.821617e-05  max resid 5.703412e-05 
    ## ... Similar to previous best
    ## Run 291 stress 0.05036801 
    ## ... Procrustes: rmse 0.0001697466  max resid 0.0003706727 
    ## ... Similar to previous best
    ## Run 292 stress 0.05051335 
    ## ... Procrustes: rmse 0.03564267  max resid 0.1216369 
    ## Run 293 stress 0.06136675 
    ## Run 294 stress 0.05337579 
    ## Run 295 stress 0.05051351 
    ## ... Procrustes: rmse 0.03572555  max resid 0.1220614 
    ## Run 296 stress 0.05036794 
    ## ... Procrustes: rmse 7.906694e-05  max resid 0.0002000072 
    ## ... Similar to previous best
    ## Run 297 stress 0.06885099 
    ## Run 298 stress 0.05036803 
    ## ... Procrustes: rmse 0.0001525557  max resid 0.0003781319 
    ## ... Similar to previous best
    ## Run 299 stress 0.05337577 
    ## Run 300 stress 0.05466708 
    ## Run 301 stress 0.05171059 
    ## Run 302 stress 0.0597923 
    ## Run 303 stress 0.05051283 
    ## ... Procrustes: rmse 0.03569615  max resid 0.1218577 
    ## Run 304 stress 0.06537577 
    ## Run 305 stress 0.05900148 
    ## Run 306 stress 0.05051331 
    ## ... Procrustes: rmse 0.03575072  max resid 0.1219752 
    ## Run 307 stress 0.05036792 
    ## ... Procrustes: rmse 5.251568e-05  max resid 0.0001300734 
    ## ... Similar to previous best
    ## Run 308 stress 0.05051316 
    ## ... Procrustes: rmse 0.03565187  max resid 0.1216824 
    ## Run 309 stress 0.05171059 
    ## Run 310 stress 0.06478517 
    ## Run 311 stress 0.05602064 
    ## Run 312 stress 0.05051358 
    ## ... Procrustes: rmse 0.03571808  max resid 0.1220408 
    ## Run 313 stress 0.06730202 
    ## Run 314 stress 0.05968626 
    ## Run 315 stress 0.06640329 
    ## Run 316 stress 0.05337596 
    ## Run 317 stress 0.06478529 
    ## Run 318 stress 0.05036793 
    ## ... Procrustes: rmse 5.434152e-05  max resid 0.0001262849 
    ## ... Similar to previous best
    ## Run 319 stress 0.05036806 
    ## ... Procrustes: rmse 0.0002021034  max resid 0.0005691672 
    ## ... Similar to previous best
    ## Run 320 stress 0.05051321 
    ## ... Procrustes: rmse 0.03568626  max resid 0.1217794 
    ## Run 321 stress 0.05036795 
    ## ... Procrustes: rmse 8.381106e-05  max resid 0.0002006656 
    ## ... Similar to previous best
    ## Run 322 stress 0.05171058 
    ## Run 323 stress 0.05900148 
    ## Run 324 stress 0.0703329 
    ## Run 325 stress 0.06136685 
    ## Run 326 stress 0.05051326 
    ## ... Procrustes: rmse 0.03568536  max resid 0.1217728 
    ## Run 327 stress 0.05928237 
    ## Run 328 stress 0.05337576 
    ## Run 329 stress 0.05051316 
    ## ... Procrustes: rmse 0.03569124  max resid 0.1218003 
    ## Run 330 stress 0.06136681 
    ## Run 331 stress 0.06231301 
    ## Run 332 stress 0.0647623 
    ## Run 333 stress 0.05051346 
    ## ... Procrustes: rmse 0.0357198  max resid 0.1220414 
    ## Run 334 stress 0.05051291 
    ## ... Procrustes: rmse 0.03565797  max resid 0.1218292 
    ## Run 335 stress 0.06896885 
    ## Run 336 stress 0.05036798 
    ## ... Procrustes: rmse 0.0001368795  max resid 0.0003477336 
    ## ... Similar to previous best
    ## Run 337 stress 0.07026164 
    ## Run 338 stress 0.06136678 
    ## Run 339 stress 0.06872618 
    ## Run 340 stress 0.05036796 
    ## ... Procrustes: rmse 8.023244e-05  max resid 0.000174503 
    ## ... Similar to previous best
    ## Run 341 stress 0.06409274 
    ## Run 342 stress 0.050513 
    ## ... Procrustes: rmse 0.03574467  max resid 0.1220942 
    ## Run 343 stress 0.06231309 
    ## Run 344 stress 0.0663053 
    ## Run 345 stress 0.05171058 
    ## Run 346 stress 0.06655663 
    ## Run 347 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001758897  max resid 0.0004895098 
    ## ... Similar to previous best
    ## Run 348 stress 0.06231287 
    ## Run 349 stress 0.06309843 
    ## Run 350 stress 0.05928231 
    ## Run 351 stress 0.05051323 
    ## ... Procrustes: rmse 0.03575294  max resid 0.1219821 
    ## Run 352 stress 0.05036792 
    ## ... Procrustes: rmse 2.807254e-05  max resid 6.437521e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.05928223 
    ## Run 354 stress 0.0517106 
    ## Run 355 stress 0.05829933 
    ## Run 356 stress 0.06231319 
    ## Run 357 stress 0.05968629 
    ## Run 358 stress 0.05051333 
    ## ... Procrustes: rmse 0.03567345  max resid 0.1217312 
    ## Run 359 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001271878  max resid 0.0003125549 
    ## ... Similar to previous best
    ## Run 360 stress 0.05928222 
    ## Run 361 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001665816  max resid 0.0003491914 
    ## ... Similar to previous best
    ## Run 362 stress 0.05051317 
    ## ... Procrustes: rmse 0.03568553  max resid 0.1217811 
    ## Run 363 stress 0.06938775 
    ## Run 364 stress 0.05171059 
    ## Run 365 stress 0.05051316 
    ## ... Procrustes: rmse 0.03569943  max resid 0.1219704 
    ## Run 366 stress 0.05051317 
    ## ... Procrustes: rmse 0.03567706  max resid 0.1219039 
    ## Run 367 stress 0.05036796 
    ## ... Procrustes: rmse 9.177014e-05  max resid 0.0002106159 
    ## ... Similar to previous best
    ## Run 368 stress 0.05051349 
    ## ... Procrustes: rmse 0.03570997  max resid 0.1218308 
    ## Run 369 stress 0.05051318 
    ## ... Procrustes: rmse 0.0357331  max resid 0.1219263 
    ## Run 370 stress 0.05036797 
    ## ... Procrustes: rmse 9.608774e-05  max resid 0.0002354135 
    ## ... Similar to previous best
    ## Run 371 stress 0.06111773 
    ## Run 372 stress 0.05829937 
    ## Run 373 stress 0.06504247 
    ## Run 374 stress 0.05968626 
    ## Run 375 stress 0.06038959 
    ## Run 376 stress 0.06231296 
    ## Run 377 stress 0.05051284 
    ## ... Procrustes: rmse 0.03573079  max resid 0.1220378 
    ## Run 378 stress 0.05051341 
    ## ... Procrustes: rmse 0.03570137  max resid 0.121807 
    ## Run 379 stress 0.05171059 
    ## Run 380 stress 0.06904478 
    ## Run 381 stress 0.05051342 
    ## ... Procrustes: rmse 0.03575402  max resid 0.1221394 
    ## Run 382 stress 0.05602049 
    ## Run 383 stress 0.05900147 
    ## Run 384 stress 0.06743844 
    ## Run 385 stress 0.05900155 
    ## Run 386 stress 0.06231302 
    ## Run 387 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001725711  max resid 0.0004206924 
    ## ... Similar to previous best
    ## Run 388 stress 0.05900152 
    ## Run 389 stress 0.06111769 
    ## Run 390 stress 0.05403566 
    ## Run 391 stress 0.0517106 
    ## Run 392 stress 0.05829932 
    ## Run 393 stress 0.05036793 
    ## ... Procrustes: rmse 6.113798e-05  max resid 0.0001473265 
    ## ... Similar to previous best
    ## Run 394 stress 0.05403591 
    ## Run 395 stress 0.05051329 
    ## ... Procrustes: rmse 0.03567739  max resid 0.121746 
    ## Run 396 stress 0.05051364 
    ## ... Procrustes: rmse 0.0357432  max resid 0.1221143 
    ## Run 397 stress 0.05928224 
    ## Run 398 stress 0.06231311 
    ## Run 399 stress 0.06885088 
    ## Run 400 stress 0.05337585 
    ## Run 401 stress 0.05171065 
    ## Run 402 stress 0.05968636 
    ## Run 403 stress 0.05829943 
    ## Run 404 stress 0.05036804 
    ## ... Procrustes: rmse 0.0001630528  max resid 0.0003951567 
    ## ... Similar to previous best
    ## Run 405 stress 0.05466702 
    ## Run 406 stress 0.05051308 
    ## ... Procrustes: rmse 0.03570028  max resid 0.1219684 
    ## Run 407 stress 0.05051335 
    ## ... Procrustes: rmse 0.03572682  max resid 0.1220599 
    ## Run 408 stress 0.05051319 
    ## ... Procrustes: rmse 0.03566958  max resid 0.1217337 
    ## Run 409 stress 0.06038961 
    ## Run 410 stress 0.05602047 
    ## Run 411 stress 0.05051318 
    ## ... Procrustes: rmse 0.03555893  max resid 0.1214168 
    ## Run 412 stress 0.05602055 
    ## Run 413 stress 0.0517106 
    ## Run 414 stress 0.05036791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.803951e-05  max resid 5.009625e-05 
    ## ... Similar to previous best
    ## Run 415 stress 0.06478549 
    ## Run 416 stress 0.050513 
    ## ... Procrustes: rmse 0.03568041  max resid 0.1217913 
    ## Run 417 stress 0.06499234 
    ## Run 418 stress 0.05036795 
    ## ... Procrustes: rmse 9.899373e-05  max resid 0.0002458001 
    ## ... Similar to previous best
    ## Run 419 stress 0.06231288 
    ## Run 420 stress 0.06136682 
    ## Run 421 stress 0.0613668 
    ## Run 422 stress 0.0560206 
    ## Run 423 stress 0.05051326 
    ## ... Procrustes: rmse 0.03572542  max resid 0.1220554 
    ## Run 424 stress 0.06674555 
    ## Run 425 stress 0.05051294 
    ## ... Procrustes: rmse 0.03565168  max resid 0.1218171 
    ## Run 426 stress 0.06904519 
    ## Run 427 stress 0.06038955 
    ## Run 428 stress 0.05979235 
    ## Run 429 stress 0.05403603 
    ## Run 430 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001374986  max resid 0.0003431444 
    ## ... Similar to previous best
    ## Run 431 stress 0.05979243 
    ## Run 432 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001828779  max resid 0.00046165 
    ## ... Similar to previous best
    ## Run 433 stress 0.0517106 
    ## Run 434 stress 0.05036804 
    ## ... Procrustes: rmse 0.000176003  max resid 0.0004418265 
    ## ... Similar to previous best
    ## Run 435 stress 0.05036792 
    ## ... Procrustes: rmse 1.461906e-05  max resid 2.738657e-05 
    ## ... Similar to previous best
    ## Run 436 stress 0.05036794 
    ## ... Procrustes: rmse 8.450232e-05  max resid 0.0002146811 
    ## ... Similar to previous best
    ## Run 437 stress 0.05036793 
    ## ... Procrustes: rmse 6.70132e-05  max resid 0.0001694403 
    ## ... Similar to previous best
    ## Run 438 stress 0.06674551 
    ## Run 439 stress 0.05968628 
    ## Run 440 stress 0.05928241 
    ## Run 441 stress 0.05051296 
    ## ... Procrustes: rmse 0.03570524  max resid 0.1219789 
    ## Run 442 stress 0.05036794 
    ## ... Procrustes: rmse 7.260836e-05  max resid 0.0001989468 
    ## ... Similar to previous best
    ## Run 443 stress 0.05051283 
    ## ... Procrustes: rmse 0.03568414  max resid 0.1218268 
    ## Run 444 stress 0.05051336 
    ## ... Procrustes: rmse 0.03572768  max resid 0.1220666 
    ## Run 445 stress 0.05051343 
    ## ... Procrustes: rmse 0.03564039  max resid 0.1216285 
    ## Run 446 stress 0.05466702 
    ## Run 447 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001647798  max resid 0.000415164 
    ## ... Similar to previous best
    ## Run 448 stress 0.06136681 
    ## Run 449 stress 0.05171059 
    ## Run 450 stress 0.06111766 
    ## Run 451 stress 0.05928223 
    ## Run 452 stress 0.06476229 
    ## Run 453 stress 0.05051373 
    ## ... Procrustes: rmse 0.03573291  max resid 0.1220905 
    ## Run 454 stress 0.06309836 
    ## Run 455 stress 0.05403628 
    ## Run 456 stress 0.06704945 
    ## Run 457 stress 0.05337595 
    ## Run 458 stress 0.05036793 
    ## ... Procrustes: rmse 6.781846e-05  max resid 0.0001700535 
    ## ... Similar to previous best
    ## Run 459 stress 0.05171058 
    ## Run 460 stress 0.05036794 
    ## ... Procrustes: rmse 9.641103e-05  max resid 0.0002385954 
    ## ... Similar to previous best
    ## Run 461 stress 0.05051342 
    ## ... Procrustes: rmse 0.0356527  max resid 0.1216651 
    ## Run 462 stress 0.05979228 
    ## Run 463 stress 0.05928234 
    ## Run 464 stress 0.05051339 
    ## ... Procrustes: rmse 0.03569683  max resid 0.1219753 
    ## Run 465 stress 0.05036792 
    ## ... Procrustes: rmse 3.804225e-05  max resid 8.941266e-05 
    ## ... Similar to previous best
    ## Run 466 stress 0.05036802 
    ## ... Procrustes: rmse 0.0001020435  max resid 0.0002435281 
    ## ... Similar to previous best
    ## Run 467 stress 0.05337581 
    ## Run 468 stress 0.05979226 
    ## Run 469 stress 0.05466701 
    ## Run 470 stress 0.05051319 
    ## ... Procrustes: rmse 0.03572962  max resid 0.1219217 
    ## Run 471 stress 0.0505131 
    ## ... Procrustes: rmse 0.03570179  max resid 0.121976 
    ## Run 472 stress 0.06896891 
    ## Run 473 stress 0.05968626 
    ## Run 474 stress 0.05036805 
    ## ... Procrustes: rmse 0.0001652868  max resid 0.0004028749 
    ## ... Similar to previous best
    ## Run 475 stress 0.05171062 
    ## Run 476 stress 0.0611177 
    ## Run 477 stress 0.05036796 
    ## ... Procrustes: rmse 0.000111426  max resid 0.0002757606 
    ## ... Similar to previous best
    ## Run 478 stress 0.0517106 
    ## Run 479 stress 0.05051331 
    ## ... Procrustes: rmse 0.03574912  max resid 0.121967 
    ## Run 480 stress 0.05979228 
    ## Run 481 stress 0.05466705 
    ## Run 482 stress 0.06027441 
    ## Run 483 stress 0.05171058 
    ## Run 484 stress 0.05466703 
    ## Run 485 stress 0.05968637 
    ## Run 486 stress 0.06789747 
    ## Run 487 stress 0.05036792 
    ## ... Procrustes: rmse 4.861587e-05  max resid 0.0001227112 
    ## ... Similar to previous best
    ## Run 488 stress 0.3657343 
    ## Run 489 stress 0.05051316 
    ## ... Procrustes: rmse 0.03571055  max resid 0.1220071 
    ## Run 490 stress 0.0590016 
    ## Run 491 stress 0.05829929 
    ## Run 492 stress 0.05968633 
    ## Run 493 stress 0.05602048 
    ## Run 494 stress 0.05829933 
    ## Run 495 stress 0.05051323 
    ## ... Procrustes: rmse 0.03568548  max resid 0.1217819 
    ## Run 496 stress 0.05968625 
    ## Run 497 stress 0.05036795 
    ## ... Procrustes: rmse 0.0001025929  max resid 0.0002578415 
    ## ... Similar to previous best
    ## Run 498 stress 0.05036799 
    ## ... Procrustes: rmse 0.0001435305  max resid 0.0003607904 
    ## ... Similar to previous best
    ## Run 499 stress 0.06136676 
    ## Run 500 stress 0.05968626 
    ## *** Best solution repeated 19 times

``` r
# Ocean sites and mixed lakes
PD_beta_geo_OM_NMDS <- metaMDS(PD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1016168 
    ## Run 1 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0006143089  max resid 0.001720355 
    ## ... Similar to previous best
    ## Run 2 stress 0.1537554 
    ## Run 3 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 1.327855e-05  max resid 3.306773e-05 
    ## ... Similar to previous best
    ## Run 4 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 1.031852e-05  max resid 2.77126e-05 
    ## ... Similar to previous best
    ## Run 5 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 1.272401e-05  max resid 3.294286e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.1341868 
    ## Run 7 stress 0.1016166 
    ## ... Procrustes: rmse 0.0004023251  max resid 0.001128559 
    ## ... Similar to previous best
    ## Run 8 stress 0.1356147 
    ## Run 9 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005677355  max resid 0.001587458 
    ## ... Similar to previous best
    ## Run 10 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006101691  max resid 0.001706188 
    ## ... Similar to previous best
    ## Run 11 stress 0.1016164 
    ## ... Procrustes: rmse 5.526452e-05  max resid 0.0001543201 
    ## ... Similar to previous best
    ## Run 12 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 2.733063e-05  max resid 7.684386e-05 
    ## ... Similar to previous best
    ## Run 13 stress 0.1016164 
    ## ... Procrustes: rmse 7.247347e-05  max resid 0.000198113 
    ## ... Similar to previous best
    ## Run 14 stress 0.1016164 
    ## ... Procrustes: rmse 3.85657e-05  max resid 0.0001086859 
    ## ... Similar to previous best
    ## Run 15 stress 0.3156808 
    ## Run 16 stress 0.1016164 
    ## ... Procrustes: rmse 3.646967e-05  max resid 0.0001024271 
    ## ... Similar to previous best
    ## Run 17 stress 0.1341868 
    ## Run 18 stress 0.3411811 
    ## Run 19 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 6.870525e-06  max resid 1.884275e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.1016164 
    ## ... Procrustes: rmse 3.577749e-05  max resid 0.0001008013 
    ## ... Similar to previous best
    ## Run 21 stress 0.1341868 
    ## Run 22 stress 0.1016164 
    ## ... Procrustes: rmse 4.188231e-05  max resid 0.0001190811 
    ## ... Similar to previous best
    ## Run 23 stress 0.1016167 
    ## ... Procrustes: rmse 0.0002788204  max resid 0.0007722716 
    ## ... Similar to previous best
    ## Run 24 stress 0.1016164 
    ## ... Procrustes: rmse 5.988306e-06  max resid 1.604191e-05 
    ## ... Similar to previous best
    ## Run 25 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001935811  max resid 0.0005412964 
    ## ... Similar to previous best
    ## Run 26 stress 0.1016164 
    ## ... Procrustes: rmse 8.933175e-05  max resid 0.0002519213 
    ## ... Similar to previous best
    ## Run 27 stress 0.1016164 
    ## ... Procrustes: rmse 9.726255e-05  max resid 0.0002732009 
    ## ... Similar to previous best
    ## Run 28 stress 0.1341868 
    ## Run 29 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005906097  max resid 0.001651127 
    ## ... Similar to previous best
    ## Run 30 stress 0.135615 
    ## Run 31 stress 0.1341868 
    ## Run 32 stress 0.1356148 
    ## Run 33 stress 0.1341868 
    ## Run 34 stress 0.1016164 
    ## ... Procrustes: rmse 8.717544e-05  max resid 0.0002461635 
    ## ... Similar to previous best
    ## Run 35 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001058321  max resid 0.0002828342 
    ## ... Similar to previous best
    ## Run 36 stress 0.1341868 
    ## Run 37 stress 0.1356149 
    ## Run 38 stress 0.1016164 
    ## ... Procrustes: rmse 1.944094e-05  max resid 5.453293e-05 
    ## ... Similar to previous best
    ## Run 39 stress 0.1341868 
    ## Run 40 stress 0.1356147 
    ## Run 41 stress 0.1341868 
    ## Run 42 stress 0.1016164 
    ## ... Procrustes: rmse 8.701077e-06  max resid 1.481939e-05 
    ## ... Similar to previous best
    ## Run 43 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002446474  max resid 0.000684391 
    ## ... Similar to previous best
    ## Run 44 stress 0.1016164 
    ## ... Procrustes: rmse 9.179313e-06  max resid 2.395617e-05 
    ## ... Similar to previous best
    ## Run 45 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004965438  max resid 0.001388493 
    ## ... Similar to previous best
    ## Run 46 stress 0.1356152 
    ## Run 47 stress 0.1016164 
    ## ... Procrustes: rmse 1.34076e-05  max resid 3.766681e-05 
    ## ... Similar to previous best
    ## Run 48 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005502135  max resid 0.001536647 
    ## ... Similar to previous best
    ## Run 49 stress 0.1530428 
    ## Run 50 stress 0.1341868 
    ## Run 51 stress 0.1016164 
    ## ... Procrustes: rmse 5.531501e-06  max resid 1.542646e-05 
    ## ... Similar to previous best
    ## Run 52 stress 0.1016172 
    ## ... Procrustes: rmse 0.0008041958  max resid 0.002247385 
    ## ... Similar to previous best
    ## Run 53 stress 0.135615 
    ## Run 54 stress 0.1341868 
    ## Run 55 stress 0.1016164 
    ## ... Procrustes: rmse 4.946643e-05  max resid 0.0001394092 
    ## ... Similar to previous best
    ## Run 56 stress 0.1016164 
    ## ... Procrustes: rmse 3.575207e-05  max resid 9.945923e-05 
    ## ... Similar to previous best
    ## Run 57 stress 0.1341868 
    ## Run 58 stress 0.1016164 
    ## ... Procrustes: rmse 3.067765e-05  max resid 8.359152e-05 
    ## ... Similar to previous best
    ## Run 59 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004290674  max resid 0.001201306 
    ## ... Similar to previous best
    ## Run 60 stress 0.3361428 
    ## Run 61 stress 0.1016164 
    ## ... Procrustes: rmse 7.443578e-06  max resid 2.048305e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005579746  max resid 0.001559981 
    ## ... Similar to previous best
    ## Run 63 stress 0.1016164 
    ## ... Procrustes: rmse 6.673507e-05  max resid 0.0001869223 
    ## ... Similar to previous best
    ## Run 64 stress 0.1016168 
    ## ... Procrustes: rmse 0.000538765  max resid 0.001504077 
    ## ... Similar to previous best
    ## Run 65 stress 0.1356147 
    ## Run 66 stress 0.1356148 
    ## Run 67 stress 0.1341868 
    ## Run 68 stress 0.1016167 
    ## ... Procrustes: rmse 0.00040062  max resid 0.001117799 
    ## ... Similar to previous best
    ## Run 69 stress 0.1016164 
    ## ... Procrustes: rmse 2.721415e-05  max resid 6.64685e-05 
    ## ... Similar to previous best
    ## Run 70 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003988272  max resid 0.001115722 
    ## ... Similar to previous best
    ## Run 71 stress 0.1016164 
    ## ... Procrustes: rmse 9.735624e-05  max resid 0.0002728725 
    ## ... Similar to previous best
    ## Run 72 stress 0.1016164 
    ## ... Procrustes: rmse 6.353468e-05  max resid 0.0001785298 
    ## ... Similar to previous best
    ## Run 73 stress 0.1016164 
    ## ... Procrustes: rmse 7.580851e-06  max resid 2.137417e-05 
    ## ... Similar to previous best
    ## Run 74 stress 0.1519033 
    ## Run 75 stress 0.1356148 
    ## Run 76 stress 0.1016168 
    ## ... Procrustes: rmse 0.000530401  max resid 0.001482068 
    ## ... Similar to previous best
    ## Run 77 stress 0.1356149 
    ## Run 78 stress 0.1341868 
    ## Run 79 stress 0.1341868 
    ## Run 80 stress 0.1016164 
    ## ... Procrustes: rmse 5.520201e-05  max resid 0.0001555153 
    ## ... Similar to previous best
    ## Run 81 stress 0.1016164 
    ## ... Procrustes: rmse 5.781599e-05  max resid 0.0001630282 
    ## ... Similar to previous best
    ## Run 82 stress 0.1356144 
    ## Run 83 stress 0.1016164 
    ## ... Procrustes: rmse 4.778751e-06  max resid 1.191499e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006195753  max resid 0.00173269 
    ## ... Similar to previous best
    ## Run 85 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007345395  max resid 0.002054115 
    ## ... Similar to previous best
    ## Run 86 stress 0.1016164 
    ## ... Procrustes: rmse 8.544615e-05  max resid 0.0002411869 
    ## ... Similar to previous best
    ## Run 87 stress 0.1016164 
    ## ... Procrustes: rmse 2.820194e-05  max resid 7.963365e-05 
    ## ... Similar to previous best
    ## Run 88 stress 0.1016164 
    ## ... Procrustes: rmse 1.556593e-05  max resid 3.919301e-05 
    ## ... Similar to previous best
    ## Run 89 stress 0.1016164 
    ## ... Procrustes: rmse 2.692716e-05  max resid 7.56869e-05 
    ## ... Similar to previous best
    ## Run 90 stress 0.1016164 
    ## ... Procrustes: rmse 2.294395e-05  max resid 6.368981e-05 
    ## ... Similar to previous best
    ## Run 91 stress 0.1016164 
    ## ... Procrustes: rmse 6.292541e-05  max resid 0.0001771459 
    ## ... Similar to previous best
    ## Run 92 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002084016  max resid 0.0005844852 
    ## ... Similar to previous best
    ## Run 93 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002729005  max resid 0.000763805 
    ## ... Similar to previous best
    ## Run 94 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001253953  max resid 0.0003537284 
    ## ... Similar to previous best
    ## Run 95 stress 0.1530422 
    ## Run 96 stress 0.1016164 
    ## ... Procrustes: rmse 8.458843e-05  max resid 0.0002384306 
    ## ... Similar to previous best
    ## Run 97 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002614193  max resid 0.0007333898 
    ## ... Similar to previous best
    ## Run 98 stress 0.1016164 
    ## ... Procrustes: rmse 2.713449e-05  max resid 7.593151e-05 
    ## ... Similar to previous best
    ## Run 99 stress 0.1016164 
    ## ... Procrustes: rmse 6.184851e-05  max resid 0.0001741218 
    ## ... Similar to previous best
    ## Run 100 stress 0.1016164 
    ## ... Procrustes: rmse 3.4773e-06  max resid 6.2429e-06 
    ## ... Similar to previous best
    ## Run 101 stress 0.1341868 
    ## Run 102 stress 0.1356149 
    ## Run 103 stress 0.1016164 
    ## ... Procrustes: rmse 5.676161e-06  max resid 1.542892e-05 
    ## ... Similar to previous best
    ## Run 104 stress 0.1341868 
    ## Run 105 stress 0.1016164 
    ## ... Procrustes: rmse 4.685689e-05  max resid 0.000130907 
    ## ... Similar to previous best
    ## Run 106 stress 0.1530423 
    ## Run 107 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001104913  max resid 0.0003102293 
    ## ... Similar to previous best
    ## Run 108 stress 0.1016164 
    ## ... Procrustes: rmse 4.508469e-06  max resid 1.238444e-05 
    ## ... Similar to previous best
    ## Run 109 stress 0.101617 
    ## ... Procrustes: rmse 0.0006497027  max resid 0.001813165 
    ## ... Similar to previous best
    ## Run 110 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004802131  max resid 0.001343246 
    ## ... Similar to previous best
    ## Run 111 stress 0.1341868 
    ## Run 112 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004946051  max resid 0.001382903 
    ## ... Similar to previous best
    ## Run 113 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007495577  max resid 0.002092662 
    ## ... Similar to previous best
    ## Run 114 stress 0.1341868 
    ## Run 115 stress 0.101617 
    ## ... Procrustes: rmse 0.0006319561  max resid 0.001766884 
    ## ... Similar to previous best
    ## Run 116 stress 0.1341868 
    ## Run 117 stress 0.1016164 
    ## ... Procrustes: rmse 2.768266e-05  max resid 7.845002e-05 
    ## ... Similar to previous best
    ## Run 118 stress 0.1519033 
    ## Run 119 stress 0.1016164 
    ## ... Procrustes: rmse 6.651372e-05  max resid 0.0001872692 
    ## ... Similar to previous best
    ## Run 120 stress 0.3455887 
    ## Run 121 stress 0.1341868 
    ## Run 122 stress 0.101617 
    ## ... Procrustes: rmse 0.0007059409  max resid 0.001971783 
    ## ... Similar to previous best
    ## Run 123 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004761126  max resid 0.001329484 
    ## ... Similar to previous best
    ## Run 124 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003382814  max resid 0.0009472952 
    ## ... Similar to previous best
    ## Run 125 stress 0.1016164 
    ## ... Procrustes: rmse 1.626466e-05  max resid 4.578888e-05 
    ## ... Similar to previous best
    ## Run 126 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004730235  max resid 0.001323418 
    ## ... Similar to previous best
    ## Run 127 stress 0.1341868 
    ## Run 128 stress 0.1016164 
    ## ... Procrustes: rmse 2.36209e-05  max resid 6.570307e-05 
    ## ... Similar to previous best
    ## Run 129 stress 0.101617 
    ## ... Procrustes: rmse 0.00064585  max resid 0.001802586 
    ## ... Similar to previous best
    ## Run 130 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 1.881044e-06  max resid 4.619852e-06 
    ## ... Similar to previous best
    ## Run 131 stress 0.1341868 
    ## Run 132 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002132609  max resid 0.0005980743 
    ## ... Similar to previous best
    ## Run 133 stress 0.1341868 
    ## Run 134 stress 0.1356146 
    ## Run 135 stress 0.1356146 
    ## Run 136 stress 0.1016164 
    ## ... Procrustes: rmse 1.383811e-05  max resid 3.88503e-05 
    ## ... Similar to previous best
    ## Run 137 stress 0.1016164 
    ## ... Procrustes: rmse 6.764614e-05  max resid 0.0001851143 
    ## ... Similar to previous best
    ## Run 138 stress 0.1016164 
    ## ... Procrustes: rmse 1.494356e-05  max resid 4.206514e-05 
    ## ... Similar to previous best
    ## Run 139 stress 0.1016164 
    ## ... Procrustes: rmse 9.246493e-06  max resid 2.624857e-05 
    ## ... Similar to previous best
    ## Run 140 stress 0.3014422 
    ## Run 141 stress 0.1016164 
    ## ... Procrustes: rmse 3.533066e-06  max resid 5.510098e-06 
    ## ... Similar to previous best
    ## Run 142 stress 0.1016164 
    ## ... Procrustes: rmse 2.739292e-05  max resid 7.68729e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.1356145 
    ## Run 144 stress 0.1341868 
    ## Run 145 stress 0.1016164 
    ## ... Procrustes: rmse 7.620809e-05  max resid 0.000213907 
    ## ... Similar to previous best
    ## Run 146 stress 0.1356145 
    ## Run 147 stress 0.1016164 
    ## ... Procrustes: rmse 4.417578e-05  max resid 0.0001244781 
    ## ... Similar to previous best
    ## Run 148 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001919615  max resid 0.0005385644 
    ## ... Similar to previous best
    ## Run 149 stress 0.1341868 
    ## Run 150 stress 0.1016164 
    ## ... Procrustes: rmse 6.088295e-05  max resid 0.0001711694 
    ## ... Similar to previous best
    ## Run 151 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 1.188428e-06  max resid 2.107937e-06 
    ## ... Similar to previous best
    ## Run 152 stress 0.1530425 
    ## Run 153 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004816894  max resid 0.001346125 
    ## ... Similar to previous best
    ## Run 154 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003203611  max resid 0.0008983163 
    ## ... Similar to previous best
    ## Run 155 stress 0.1530423 
    ## Run 156 stress 0.1530426 
    ## Run 157 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001207085  max resid 0.0003407948 
    ## ... Similar to previous best
    ## Run 158 stress 0.1016164 
    ## ... Procrustes: rmse 4.640053e-05  max resid 0.0001309955 
    ## ... Similar to previous best
    ## Run 159 stress 0.1016164 
    ## ... Procrustes: rmse 2.434026e-05  max resid 6.878993e-05 
    ## ... Similar to previous best
    ## Run 160 stress 0.1016164 
    ## ... Procrustes: rmse 4.109965e-05  max resid 0.0001140257 
    ## ... Similar to previous best
    ## Run 161 stress 0.1016164 
    ## ... Procrustes: rmse 2.359601e-05  max resid 6.650009e-05 
    ## ... Similar to previous best
    ## Run 162 stress 0.1016164 
    ## ... Procrustes: rmse 8.937596e-06  max resid 2.514098e-05 
    ## ... Similar to previous best
    ## Run 163 stress 0.1016164 
    ## ... Procrustes: rmse 4.067952e-05  max resid 0.0001142294 
    ## ... Similar to previous best
    ## Run 164 stress 0.1519033 
    ## Run 165 stress 0.1356148 
    ## Run 166 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002268729  max resid 0.0006363268 
    ## ... Similar to previous best
    ## Run 167 stress 0.1341868 
    ## Run 168 stress 0.1341868 
    ## Run 169 stress 0.1341868 
    ## Run 170 stress 0.1341868 
    ## Run 171 stress 0.1016164 
    ## ... Procrustes: rmse 2.991597e-06  max resid 7.295493e-06 
    ## ... Similar to previous best
    ## Run 172 stress 0.1356149 
    ## Run 173 stress 0.1530424 
    ## Run 174 stress 0.1016164 
    ## ... Procrustes: rmse 1.643317e-06  max resid 4.107556e-06 
    ## ... Similar to previous best
    ## Run 175 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003193455  max resid 0.0008943452 
    ## ... Similar to previous best
    ## Run 176 stress 0.1016164 
    ## ... Procrustes: rmse 4.864622e-05  max resid 0.0001372533 
    ## ... Similar to previous best
    ## Run 177 stress 0.1016164 
    ## ... Procrustes: rmse 6.429314e-05  max resid 0.0001806178 
    ## ... Similar to previous best
    ## Run 178 stress 0.1341868 
    ## Run 179 stress 0.101617 
    ## ... Procrustes: rmse 0.0005993067  max resid 0.001670009 
    ## ... Similar to previous best
    ## Run 180 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005530098  max resid 0.001544072 
    ## ... Similar to previous best
    ## Run 181 stress 0.1016164 
    ## ... Procrustes: rmse 8.739202e-06  max resid 1.950178e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005260168  max resid 0.001469762 
    ## ... Similar to previous best
    ## Run 183 stress 0.1341868 
    ## Run 184 stress 0.1341868 
    ## Run 185 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007147202  max resid 0.001997493 
    ## ... Similar to previous best
    ## Run 186 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004281064  max resid 0.001196293 
    ## ... Similar to previous best
    ## Run 187 stress 0.1016164 
    ## ... Procrustes: rmse 4.872951e-05  max resid 0.0001361949 
    ## ... Similar to previous best
    ## Run 188 stress 0.1016164 
    ## ... New best solution
    ## ... Procrustes: rmse 3.049864e-06  max resid 8.657225e-06 
    ## ... Similar to previous best
    ## Run 189 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005036375  max resid 0.001405669 
    ## ... Similar to previous best
    ## Run 190 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005126802  max resid 0.001432312 
    ## ... Similar to previous best
    ## Run 191 stress 0.1016164 
    ## ... Procrustes: rmse 3.784048e-05  max resid 0.0001067919 
    ## ... Similar to previous best
    ## Run 192 stress 0.1016166 
    ## ... Procrustes: rmse 0.000236992  max resid 0.0006637629 
    ## ... Similar to previous best
    ## Run 193 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005527297  max resid 0.001543641 
    ## ... Similar to previous best
    ## Run 194 stress 0.3578488 
    ## Run 195 stress 0.1016164 
    ## ... Procrustes: rmse 7.006053e-05  max resid 0.0001971059 
    ## ... Similar to previous best
    ## Run 196 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004675943  max resid 0.001307614 
    ## ... Similar to previous best
    ## Run 197 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006392581  max resid 0.001785774 
    ## ... Similar to previous best
    ## Run 198 stress 0.1341868 
    ## Run 199 stress 0.1341868 
    ## Run 200 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004267572  max resid 0.001194447 
    ## ... Similar to previous best
    ## Run 201 stress 0.325277 
    ## Run 202 stress 0.1356151 
    ## Run 203 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004331381  max resid 0.001211751 
    ## ... Similar to previous best
    ## Run 204 stress 0.101617 
    ## ... Procrustes: rmse 0.0006614256  max resid 0.001846252 
    ## ... Similar to previous best
    ## Run 205 stress 0.135615 
    ## Run 206 stress 0.1016164 
    ## ... Procrustes: rmse 3.300628e-05  max resid 9.2948e-05 
    ## ... Similar to previous best
    ## Run 207 stress 0.1016164 
    ## ... Procrustes: rmse 3.374023e-05  max resid 9.521879e-05 
    ## ... Similar to previous best
    ## Run 208 stress 0.1530423 
    ## Run 209 stress 0.1016171 
    ## ... Procrustes: rmse 0.0007343266  max resid 0.00205124 
    ## ... Similar to previous best
    ## Run 210 stress 0.1356147 
    ## Run 211 stress 0.1356148 
    ## Run 212 stress 0.1356149 
    ## Run 213 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004967837  max resid 0.001386477 
    ## ... Similar to previous best
    ## Run 214 stress 0.1016164 
    ## ... Procrustes: rmse 1.33216e-05  max resid 3.744676e-05 
    ## ... Similar to previous best
    ## Run 215 stress 0.1016164 
    ## ... Procrustes: rmse 3.386745e-05  max resid 9.516583e-05 
    ## ... Similar to previous best
    ## Run 216 stress 0.1356147 
    ## Run 217 stress 0.101617 
    ## ... Procrustes: rmse 0.0006704522  max resid 0.001873365 
    ## ... Similar to previous best
    ## Run 218 stress 0.1016164 
    ## ... Procrustes: rmse 5.016577e-05  max resid 0.0001415684 
    ## ... Similar to previous best
    ## Run 219 stress 0.1341868 
    ## Run 220 stress 0.1016167 
    ## ... Procrustes: rmse 0.000451822  max resid 0.001265023 
    ## ... Similar to previous best
    ## Run 221 stress 0.1016164 
    ## ... Procrustes: rmse 4.374651e-05  max resid 0.0001229781 
    ## ... Similar to previous best
    ## Run 222 stress 0.101617 
    ## ... Procrustes: rmse 0.0006143633  max resid 0.001714447 
    ## ... Similar to previous best
    ## Run 223 stress 0.1356151 
    ## Run 224 stress 0.101617 
    ## ... Procrustes: rmse 0.0006747621  max resid 0.001886369 
    ## ... Similar to previous best
    ## Run 225 stress 0.1341868 
    ## Run 226 stress 0.101617 
    ## ... Procrustes: rmse 0.0007045013  max resid 0.001967998 
    ## ... Similar to previous best
    ## Run 227 stress 0.1530424 
    ## Run 228 stress 0.1341868 
    ## Run 229 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001321271  max resid 0.000371332 
    ## ... Similar to previous best
    ## Run 230 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004173831  max resid 0.001166831 
    ## ... Similar to previous best
    ## Run 231 stress 0.135615 
    ## Run 232 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002468459  max resid 0.0006910028 
    ## ... Similar to previous best
    ## Run 233 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001846326  max resid 0.0005186857 
    ## ... Similar to previous best
    ## Run 234 stress 0.1016164 
    ## ... Procrustes: rmse 7.028738e-05  max resid 0.000198019 
    ## ... Similar to previous best
    ## Run 235 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005853806  max resid 0.001636612 
    ## ... Similar to previous best
    ## Run 236 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005892978  max resid 0.001647127 
    ## ... Similar to previous best
    ## Run 237 stress 0.1341868 
    ## Run 238 stress 0.1016164 
    ## ... Procrustes: rmse 7.32232e-05  max resid 0.0002064162 
    ## ... Similar to previous best
    ## Run 239 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004988427  max resid 0.001394665 
    ## ... Similar to previous best
    ## Run 240 stress 0.1356151 
    ## Run 241 stress 0.1016164 
    ## ... Procrustes: rmse 4.809239e-05  max resid 0.0001355872 
    ## ... Similar to previous best
    ## Run 242 stress 0.1341868 
    ## Run 243 stress 0.1341868 
    ## Run 244 stress 0.1356146 
    ## Run 245 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004652281  max resid 0.0013026 
    ## ... Similar to previous best
    ## Run 246 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005323032  max resid 0.001488325 
    ## ... Similar to previous best
    ## Run 247 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004818891  max resid 0.001347018 
    ## ... Similar to previous best
    ## Run 248 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005445421  max resid 0.001522302 
    ## ... Similar to previous best
    ## Run 249 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004640918  max resid 0.001298025 
    ## ... Similar to previous best
    ## Run 250 stress 0.1016164 
    ## ... Procrustes: rmse 6.009901e-05  max resid 0.0001621744 
    ## ... Similar to previous best
    ## Run 251 stress 0.1341868 
    ## Run 252 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003537999  max resid 0.0009893826 
    ## ... Similar to previous best
    ## Run 253 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004679701  max resid 0.001304936 
    ## ... Similar to previous best
    ## Run 254 stress 0.1537555 
    ## Run 255 stress 0.1016164 
    ## ... Procrustes: rmse 9.153132e-06  max resid 2.571863e-05 
    ## ... Similar to previous best
    ## Run 256 stress 0.1016167 
    ## ... Procrustes: rmse 0.000401567  max resid 0.001124214 
    ## ... Similar to previous best
    ## Run 257 stress 0.1341868 
    ## Run 258 stress 0.1016164 
    ## ... Procrustes: rmse 1.581186e-05  max resid 4.408628e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.1341868 
    ## Run 260 stress 0.1356145 
    ## Run 261 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002343004  max resid 0.0006473948 
    ## ... Similar to previous best
    ## Run 262 stress 0.1016164 
    ## ... Procrustes: rmse 6.56452e-05  max resid 0.0001846607 
    ## ... Similar to previous best
    ## Run 263 stress 0.1341868 
    ## Run 264 stress 0.1016164 
    ## ... Procrustes: rmse 5.250262e-07  max resid 1.289723e-06 
    ## ... Similar to previous best
    ## Run 265 stress 0.1356147 
    ## Run 266 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006367313  max resid 0.001780192 
    ## ... Similar to previous best
    ## Run 267 stress 0.1016164 
    ## ... Procrustes: rmse 5.456632e-05  max resid 0.0001524792 
    ## ... Similar to previous best
    ## Run 268 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004881628  max resid 0.001365884 
    ## ... Similar to previous best
    ## Run 269 stress 0.1341868 
    ## Run 270 stress 0.1016164 
    ## ... Procrustes: rmse 5.850565e-06  max resid 1.611504e-05 
    ## ... Similar to previous best
    ## Run 271 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003670309  max resid 0.001021477 
    ## ... Similar to previous best
    ## Run 272 stress 0.1016164 
    ## ... Procrustes: rmse 1.552672e-05  max resid 4.317932e-05 
    ## ... Similar to previous best
    ## Run 273 stress 0.1016164 
    ## ... Procrustes: rmse 2.081908e-05  max resid 5.855173e-05 
    ## ... Similar to previous best
    ## Run 274 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005324799  max resid 0.001487163 
    ## ... Similar to previous best
    ## Run 275 stress 0.1016166 
    ## ... Procrustes: rmse 0.000396131  max resid 0.001108633 
    ## ... Similar to previous best
    ## Run 276 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002653996  max resid 0.0007434839 
    ## ... Similar to previous best
    ## Run 277 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001652176  max resid 0.0004632624 
    ## ... Similar to previous best
    ## Run 278 stress 0.1530427 
    ## Run 279 stress 0.1341868 
    ## Run 280 stress 0.1341868 
    ## Run 281 stress 0.1356149 
    ## Run 282 stress 0.1016164 
    ## ... Procrustes: rmse 3.566286e-07  max resid 5.933926e-07 
    ## ... Similar to previous best
    ## Run 283 stress 0.1537555 
    ## Run 284 stress 0.1016164 
    ## ... Procrustes: rmse 4.372032e-06  max resid 1.153683e-05 
    ## ... Similar to previous best
    ## Run 285 stress 0.1356149 
    ## Run 286 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004332324  max resid 0.0012129 
    ## ... Similar to previous best
    ## Run 287 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003694197  max resid 0.001033072 
    ## ... Similar to previous best
    ## Run 288 stress 0.1341868 
    ## Run 289 stress 0.1341868 
    ## Run 290 stress 0.1016164 
    ## ... Procrustes: rmse 1.048449e-05  max resid 2.962588e-05 
    ## ... Similar to previous best
    ## Run 291 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004257611  max resid 0.00118977 
    ## ... Similar to previous best
    ## Run 292 stress 0.101617 
    ## ... Procrustes: rmse 0.0006806638  max resid 0.001901033 
    ## ... Similar to previous best
    ## Run 293 stress 0.1356145 
    ## Run 294 stress 0.1016164 
    ## ... Procrustes: rmse 6.409986e-06  max resid 1.81115e-05 
    ## ... Similar to previous best
    ## Run 295 stress 0.1016164 
    ## ... Procrustes: rmse 1.425746e-05  max resid 3.866965e-05 
    ## ... Similar to previous best
    ## Run 296 stress 0.1016166 
    ## ... Procrustes: rmse 0.000252746  max resid 0.0006994254 
    ## ... Similar to previous best
    ## Run 297 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005676526  max resid 0.00158531 
    ## ... Similar to previous best
    ## Run 298 stress 0.1341868 
    ## Run 299 stress 0.1016164 
    ## ... Procrustes: rmse 1.63064e-05  max resid 4.436564e-05 
    ## ... Similar to previous best
    ## Run 300 stress 0.1016164 
    ## ... Procrustes: rmse 1.602502e-05  max resid 3.447935e-05 
    ## ... Similar to previous best
    ## Run 301 stress 0.1016164 
    ## ... Procrustes: rmse 3.781838e-05  max resid 0.00010625 
    ## ... Similar to previous best
    ## Run 302 stress 0.1356147 
    ## Run 303 stress 0.1356147 
    ## Run 304 stress 0.1016164 
    ## ... Procrustes: rmse 2.01443e-05  max resid 5.68448e-05 
    ## ... Similar to previous best
    ## Run 305 stress 0.1341868 
    ## Run 306 stress 0.1530426 
    ## Run 307 stress 0.1356144 
    ## Run 308 stress 0.1016169 
    ## ... Procrustes: rmse 0.000618035  max resid 0.00172648 
    ## ... Similar to previous best
    ## Run 309 stress 0.1016164 
    ## ... Procrustes: rmse 8.699103e-07  max resid 2.324687e-06 
    ## ... Similar to previous best
    ## Run 310 stress 0.1016164 
    ## ... Procrustes: rmse 2.636169e-05  max resid 7.102314e-05 
    ## ... Similar to previous best
    ## Run 311 stress 0.1530425 
    ## Run 312 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002331504  max resid 0.0006542524 
    ## ... Similar to previous best
    ## Run 313 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001092926  max resid 0.0003087503 
    ## ... Similar to previous best
    ## Run 314 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002409744  max resid 0.0006766224 
    ## ... Similar to previous best
    ## Run 315 stress 0.1016164 
    ## ... Procrustes: rmse 5.590038e-05  max resid 0.0001574573 
    ## ... Similar to previous best
    ## Run 316 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004289506  max resid 0.001199547 
    ## ... Similar to previous best
    ## Run 317 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004784408  max resid 0.00133836 
    ## ... Similar to previous best
    ## Run 318 stress 0.1537554 
    ## Run 319 stress 0.1016164 
    ## ... Procrustes: rmse 2.886405e-05  max resid 8.118558e-05 
    ## ... Similar to previous best
    ## Run 320 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005200477  max resid 0.001454445 
    ## ... Similar to previous best
    ## Run 321 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003729183  max resid 0.001042617 
    ## ... Similar to previous best
    ## Run 322 stress 0.1356149 
    ## Run 323 stress 0.101617 
    ## ... Procrustes: rmse 0.0006616063  max resid 0.001849814 
    ## ... Similar to previous best
    ## Run 324 stress 0.1016164 
    ## ... Procrustes: rmse 8.165032e-05  max resid 0.0002298996 
    ## ... Similar to previous best
    ## Run 325 stress 0.1341868 
    ## Run 326 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006131034  max resid 0.001713031 
    ## ... Similar to previous best
    ## Run 327 stress 0.1356144 
    ## Run 328 stress 0.2109396 
    ## Run 329 stress 0.1356149 
    ## Run 330 stress 0.1341868 
    ## Run 331 stress 0.1519033 
    ## Run 332 stress 0.1016164 
    ## ... Procrustes: rmse 8.453325e-07  max resid 1.609366e-06 
    ## ... Similar to previous best
    ## Run 333 stress 0.1016164 
    ## ... Procrustes: rmse 4.285911e-05  max resid 0.0001201251 
    ## ... Similar to previous best
    ## Run 334 stress 0.1016164 
    ## ... Procrustes: rmse 7.133088e-06  max resid 2.013246e-05 
    ## ... Similar to previous best
    ## Run 335 stress 0.1016171 
    ## ... Procrustes: rmse 0.0006572523  max resid 0.001831011 
    ## ... Similar to previous best
    ## Run 336 stress 0.101617 
    ## ... Procrustes: rmse 0.00070589  max resid 0.001972667 
    ## ... Similar to previous best
    ## Run 337 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004740945  max resid 0.001326116 
    ## ... Similar to previous best
    ## Run 338 stress 0.1016166 
    ## ... Procrustes: rmse 0.0002969459  max resid 0.0008295887 
    ## ... Similar to previous best
    ## Run 339 stress 0.1341868 
    ## Run 340 stress 0.1016164 
    ## ... Procrustes: rmse 4.062705e-06  max resid 9.567839e-06 
    ## ... Similar to previous best
    ## Run 341 stress 0.101617 
    ## ... Procrustes: rmse 0.0006461026  max resid 0.001804783 
    ## ... Similar to previous best
    ## Run 342 stress 0.1016164 
    ## ... Procrustes: rmse 2.086804e-05  max resid 5.886199e-05 
    ## ... Similar to previous best
    ## Run 343 stress 0.1356148 
    ## Run 344 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003590417  max resid 0.001004831 
    ## ... Similar to previous best
    ## Run 345 stress 0.1016164 
    ## ... Procrustes: rmse 9.425167e-06  max resid 2.43361e-05 
    ## ... Similar to previous best
    ## Run 346 stress 0.1341868 
    ## Run 347 stress 0.1341868 
    ## Run 348 stress 0.1341868 
    ## Run 349 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003410591  max resid 0.0009542573 
    ## ... Similar to previous best
    ## Run 350 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004482376  max resid 0.001251356 
    ## ... Similar to previous best
    ## Run 351 stress 0.1016164 
    ## ... Procrustes: rmse 3.028363e-05  max resid 8.530252e-05 
    ## ... Similar to previous best
    ## Run 352 stress 0.1530426 
    ## Run 353 stress 0.1530423 
    ## Run 354 stress 0.1016164 
    ## ... Procrustes: rmse 1.795994e-05  max resid 5.078834e-05 
    ## ... Similar to previous best
    ## Run 355 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006008679  max resid 0.001680667 
    ## ... Similar to previous best
    ## Run 356 stress 0.1341868 
    ## Run 357 stress 0.1016164 
    ## ... Procrustes: rmse 2.290504e-05  max resid 6.473657e-05 
    ## ... Similar to previous best
    ## Run 358 stress 0.1341868 
    ## Run 359 stress 0.1016164 
    ## ... Procrustes: rmse 3.222484e-05  max resid 9.078595e-05 
    ## ... Similar to previous best
    ## Run 360 stress 0.1016164 
    ## ... Procrustes: rmse 9.302687e-05  max resid 0.0002620236 
    ## ... Similar to previous best
    ## Run 361 stress 0.1341868 
    ## Run 362 stress 0.1016164 
    ## ... Procrustes: rmse 9.364818e-05  max resid 0.0002645131 
    ## ... Similar to previous best
    ## Run 363 stress 0.1016164 
    ## ... Procrustes: rmse 8.879821e-05  max resid 0.0002492967 
    ## ... Similar to previous best
    ## Run 364 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005558238  max resid 0.001553182 
    ## ... Similar to previous best
    ## Run 365 stress 0.1356144 
    ## Run 366 stress 0.1016164 
    ## ... Procrustes: rmse 7.228359e-06  max resid 2.035593e-05 
    ## ... Similar to previous best
    ## Run 367 stress 0.1341868 
    ## Run 368 stress 0.1016164 
    ## ... Procrustes: rmse 3.068022e-05  max resid 8.208425e-05 
    ## ... Similar to previous best
    ## Run 369 stress 0.1519033 
    ## Run 370 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003882844  max resid 0.0010862 
    ## ... Similar to previous best
    ## Run 371 stress 0.1356151 
    ## Run 372 stress 0.1016164 
    ## ... Procrustes: rmse 9.092658e-05  max resid 0.0002563556 
    ## ... Similar to previous best
    ## Run 373 stress 0.1016164 
    ## ... Procrustes: rmse 2.669005e-05  max resid 7.465591e-05 
    ## ... Similar to previous best
    ## Run 374 stress 0.1341868 
    ## Run 375 stress 0.1016164 
    ## ... Procrustes: rmse 1.832124e-06  max resid 3.753037e-06 
    ## ... Similar to previous best
    ## Run 376 stress 0.1341868 
    ## Run 377 stress 0.1016164 
    ## ... Procrustes: rmse 3.36537e-05  max resid 9.366443e-05 
    ## ... Similar to previous best
    ## Run 378 stress 0.1016164 
    ## ... Procrustes: rmse 2.04354e-05  max resid 5.766866e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003804521  max resid 0.001065686 
    ## ... Similar to previous best
    ## Run 380 stress 0.1016164 
    ## ... Procrustes: rmse 1.472026e-05  max resid 4.131393e-05 
    ## ... Similar to previous best
    ## Run 381 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003108567  max resid 0.0008705379 
    ## ... Similar to previous best
    ## Run 382 stress 0.1356147 
    ## Run 383 stress 0.1356153 
    ## Run 384 stress 0.1016164 
    ## ... Procrustes: rmse 2.292969e-05  max resid 6.454194e-05 
    ## ... Similar to previous best
    ## Run 385 stress 0.1341868 
    ## Run 386 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005250091  max resid 0.001467767 
    ## ... Similar to previous best
    ## Run 387 stress 0.1016164 
    ## ... Procrustes: rmse 7.279371e-05  max resid 0.000204219 
    ## ... Similar to previous best
    ## Run 388 stress 0.1016164 
    ## ... Procrustes: rmse 3.31068e-05  max resid 9.330546e-05 
    ## ... Similar to previous best
    ## Run 389 stress 0.1356146 
    ## Run 390 stress 0.1530426 
    ## Run 391 stress 0.1016164 
    ## ... Procrustes: rmse 5.949637e-05  max resid 0.0001671947 
    ## ... Similar to previous best
    ## Run 392 stress 0.1341868 
    ## Run 393 stress 0.1016164 
    ## ... Procrustes: rmse 6.025089e-06  max resid 1.422741e-05 
    ## ... Similar to previous best
    ## Run 394 stress 0.1016168 
    ## ... Procrustes: rmse 0.0004989477  max resid 0.001394735 
    ## ... Similar to previous best
    ## Run 395 stress 0.1341868 
    ## Run 396 stress 0.1016164 
    ## ... Procrustes: rmse 3.249847e-05  max resid 9.165556e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004143907  max resid 0.001158118 
    ## ... Similar to previous best
    ## Run 398 stress 0.1356149 
    ## Run 399 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003681405  max resid 0.001026401 
    ## ... Similar to previous best
    ## Run 400 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002669472  max resid 0.0007479997 
    ## ... Similar to previous best
    ## Run 401 stress 0.1356143 
    ## Run 402 stress 0.1341868 
    ## Run 403 stress 0.1016164 
    ## ... Procrustes: rmse 3.291864e-05  max resid 9.282423e-05 
    ## ... Similar to previous best
    ## Run 404 stress 0.1341868 
    ## Run 405 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005103491  max resid 0.001426038 
    ## ... Similar to previous best
    ## Run 406 stress 0.1016164 
    ## ... Procrustes: rmse 9.727937e-06  max resid 2.683333e-05 
    ## ... Similar to previous best
    ## Run 407 stress 0.1537555 
    ## Run 408 stress 0.1341868 
    ## Run 409 stress 0.1016164 
    ## ... Procrustes: rmse 3.126756e-05  max resid 8.799366e-05 
    ## ... Similar to previous best
    ## Run 410 stress 0.1016164 
    ## ... Procrustes: rmse 7.025882e-05  max resid 0.0001972927 
    ## ... Similar to previous best
    ## Run 411 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005318044  max resid 0.001487399 
    ## ... Similar to previous best
    ## Run 412 stress 0.1016164 
    ## ... Procrustes: rmse 5.304266e-06  max resid 1.485118e-05 
    ## ... Similar to previous best
    ## Run 413 stress 0.1016164 
    ## ... Procrustes: rmse 1.14046e-05  max resid 3.214168e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006178383  max resid 0.001724587 
    ## ... Similar to previous best
    ## Run 415 stress 0.1016164 
    ## ... Procrustes: rmse 8.66809e-05  max resid 0.0002442243 
    ## ... Similar to previous best
    ## Run 416 stress 0.1016164 
    ## ... Procrustes: rmse 7.619462e-05  max resid 0.0002150652 
    ## ... Similar to previous best
    ## Run 417 stress 0.1341868 
    ## Run 418 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004576766  max resid 0.001279298 
    ## ... Similar to previous best
    ## Run 419 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005497698  max resid 0.001537408 
    ## ... Similar to previous best
    ## Run 420 stress 0.1341868 
    ## Run 421 stress 0.1016164 
    ## ... Procrustes: rmse 8.221869e-05  max resid 0.000230498 
    ## ... Similar to previous best
    ## Run 422 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002309263  max resid 0.000647327 
    ## ... Similar to previous best
    ## Run 423 stress 0.1530425 
    ## Run 424 stress 0.1016164 
    ## ... Procrustes: rmse 1.775671e-05  max resid 4.06229e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004501385  max resid 0.001258115 
    ## ... Similar to previous best
    ## Run 426 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003978608  max resid 0.001108848 
    ## ... Similar to previous best
    ## Run 427 stress 0.1016164 
    ## ... Procrustes: rmse 5.316315e-05  max resid 0.0001497435 
    ## ... Similar to previous best
    ## Run 428 stress 0.1016165 
    ## ... Procrustes: rmse 0.0001759764  max resid 0.0004888436 
    ## ... Similar to previous best
    ## Run 429 stress 0.1530422 
    ## Run 430 stress 0.2074952 
    ## Run 431 stress 0.1341868 
    ## Run 432 stress 0.1016164 
    ## ... Procrustes: rmse 3.645592e-05  max resid 0.0001026761 
    ## ... Similar to previous best
    ## Run 433 stress 0.1016164 
    ## ... Procrustes: rmse 8.09994e-05  max resid 0.0002279342 
    ## ... Similar to previous best
    ## Run 434 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001134811  max resid 0.000319901 
    ## ... Similar to previous best
    ## Run 435 stress 0.1341868 
    ## Run 436 stress 0.1356149 
    ## Run 437 stress 0.1016164 
    ## ... Procrustes: rmse 2.541144e-05  max resid 7.159457e-05 
    ## ... Similar to previous best
    ## Run 438 stress 0.1016168 
    ## ... Procrustes: rmse 0.0005588052  max resid 0.001561667 
    ## ... Similar to previous best
    ## Run 439 stress 0.1016164 
    ## ... Procrustes: rmse 1.142281e-05  max resid 3.104646e-05 
    ## ... Similar to previous best
    ## Run 440 stress 0.1016167 
    ## ... Procrustes: rmse 0.0003976743  max resid 0.001111644 
    ## ... Similar to previous best
    ## Run 441 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001078716  max resid 0.0003004393 
    ## ... Similar to previous best
    ## Run 442 stress 0.1016164 
    ## ... Procrustes: rmse 7.209831e-05  max resid 0.0002031096 
    ## ... Similar to previous best
    ## Run 443 stress 0.1016164 
    ## ... Procrustes: rmse 3.219477e-05  max resid 9.074409e-05 
    ## ... Similar to previous best
    ## Run 444 stress 0.1341868 
    ## Run 445 stress 0.1016164 
    ## ... Procrustes: rmse 4.301225e-05  max resid 0.0001211291 
    ## ... Similar to previous best
    ## Run 446 stress 0.1341868 
    ## Run 447 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003691957  max resid 0.001031799 
    ## ... Similar to previous best
    ## Run 448 stress 0.1016164 
    ## ... Procrustes: rmse 2.402289e-05  max resid 6.773508e-05 
    ## ... Similar to previous best
    ## Run 449 stress 0.1341868 
    ## Run 450 stress 0.1016164 
    ## ... Procrustes: rmse 2.132064e-06  max resid 5.719439e-06 
    ## ... Similar to previous best
    ## Run 451 stress 0.1016164 
    ## ... Procrustes: rmse 3.720401e-05  max resid 0.0001020478 
    ## ... Similar to previous best
    ## Run 452 stress 0.1530425 
    ## Run 453 stress 0.1016164 
    ## ... Procrustes: rmse 1.580013e-05  max resid 4.453653e-05 
    ## ... Similar to previous best
    ## Run 454 stress 0.1016164 
    ## ... Procrustes: rmse 1.854993e-05  max resid 5.223471e-05 
    ## ... Similar to previous best
    ## Run 455 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004451223  max resid 0.001245425 
    ## ... Similar to previous best
    ## Run 456 stress 0.1356148 
    ## Run 457 stress 0.1016164 
    ## ... Procrustes: rmse 3.275204e-05  max resid 9.079718e-05 
    ## ... Similar to previous best
    ## Run 458 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005806419  max resid 0.001622076 
    ## ... Similar to previous best
    ## Run 459 stress 0.1341868 
    ## Run 460 stress 0.1016164 
    ## ... Procrustes: rmse 4.512199e-05  max resid 0.0001274874 
    ## ... Similar to previous best
    ## Run 461 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004762212  max resid 0.00133102 
    ## ... Similar to previous best
    ## Run 462 stress 0.1016164 
    ## ... Procrustes: rmse 3.195473e-05  max resid 8.987391e-05 
    ## ... Similar to previous best
    ## Run 463 stress 0.1016164 
    ## ... Procrustes: rmse 1.280354e-05  max resid 3.563953e-05 
    ## ... Similar to previous best
    ## Run 464 stress 0.2109396 
    ## Run 465 stress 0.1356152 
    ## Run 466 stress 0.1356147 
    ## Run 467 stress 0.1016164 
    ## ... Procrustes: rmse 2.512025e-05  max resid 7.088106e-05 
    ## ... Similar to previous best
    ## Run 468 stress 0.1016164 
    ## ... Procrustes: rmse 3.617291e-05  max resid 0.0001018525 
    ## ... Similar to previous best
    ## Run 469 stress 0.1341868 
    ## Run 470 stress 0.1356145 
    ## Run 471 stress 0.1016164 
    ## ... Procrustes: rmse 3.641501e-05  max resid 0.0001028995 
    ## ... Similar to previous best
    ## Run 472 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004397423  max resid 0.001229405 
    ## ... Similar to previous best
    ## Run 473 stress 0.1016164 
    ## ... Procrustes: rmse 3.168877e-05  max resid 8.9236e-05 
    ## ... Similar to previous best
    ## Run 474 stress 0.1341868 
    ## Run 475 stress 0.1016164 
    ## ... Procrustes: rmse 1.850884e-05  max resid 5.220366e-05 
    ## ... Similar to previous best
    ## Run 476 stress 0.1016165 
    ## ... Procrustes: rmse 0.0002589118  max resid 0.0007247151 
    ## ... Similar to previous best
    ## Run 477 stress 0.1356147 
    ## Run 478 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006270762  max resid 0.001751746 
    ## ... Similar to previous best
    ## Run 479 stress 0.1016166 
    ## ... Procrustes: rmse 0.0003575802  max resid 0.001000482 
    ## ... Similar to previous best
    ## Run 480 stress 0.101617 
    ## ... Procrustes: rmse 0.0006410226  max resid 0.001789523 
    ## ... Similar to previous best
    ## Run 481 stress 0.1016164 
    ## ... Procrustes: rmse 5.925799e-05  max resid 0.0001669978 
    ## ... Similar to previous best
    ## Run 482 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005916717  max resid 0.001653009 
    ## ... Similar to previous best
    ## Run 483 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001048219  max resid 0.0002953009 
    ## ... Similar to previous best
    ## Run 484 stress 0.1016169 
    ## ... Procrustes: rmse 0.0006274082  max resid 0.001751208 
    ## ... Similar to previous best
    ## Run 485 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004548334  max resid 0.001271346 
    ## ... Similar to previous best
    ## Run 486 stress 0.1016164 
    ## ... Procrustes: rmse 0.0001191411  max resid 0.0003355783 
    ## ... Similar to previous best
    ## Run 487 stress 0.1016164 
    ## ... Procrustes: rmse 8.134161e-05  max resid 0.000227956 
    ## ... Similar to previous best
    ## Run 488 stress 0.1016167 
    ## ... Procrustes: rmse 0.0004372886  max resid 0.001220081 
    ## ... Similar to previous best
    ## Run 489 stress 0.1016164 
    ## ... Procrustes: rmse 4.754918e-05  max resid 0.0001342054 
    ## ... Similar to previous best
    ## Run 490 stress 0.1016165 
    ## ... Procrustes: rmse 0.000176194  max resid 0.0004940921 
    ## ... Similar to previous best
    ## Run 491 stress 0.101617 
    ## ... Procrustes: rmse 0.0006496009  max resid 0.001812834 
    ## ... Similar to previous best
    ## Run 492 stress 0.1016164 
    ## ... Procrustes: rmse 8.810983e-06  max resid 2.22225e-05 
    ## ... Similar to previous best
    ## Run 493 stress 0.1341868 
    ## Run 494 stress 0.1356144 
    ## Run 495 stress 0.1016168 
    ## ... Procrustes: rmse 0.000523371  max resid 0.001461257 
    ## ... Similar to previous best
    ## Run 496 stress 0.1341868 
    ## Run 497 stress 0.1016167 
    ## ... Procrustes: rmse 0.000475548  max resid 0.001328991 
    ## ... Similar to previous best
    ## Run 498 stress 0.1016169 
    ## ... Procrustes: rmse 0.0005953048  max resid 0.001663121 
    ## ... Similar to previous best
    ## Run 499 stress 0.1016164 
    ## ... Procrustes: rmse 1.36366e-05  max resid 3.856161e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.1016164 
    ## ... Procrustes: rmse 7.517398e-05  max resid 0.0002109032 
    ## ... Similar to previous best
    ## *** Best solution repeated 206 times

``` r
# Stratified lakes and ocean sites
PD_beta_geo_SO_NMDS <- metaMDS(PD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.04851417 
    ## Run 1 stress 0.0573152 
    ## Run 2 stress 0.05855029 
    ## Run 3 stress 0.06088977 
    ## Run 4 stress 0.06088983 
    ## Run 5 stress 0.05855029 
    ## Run 6 stress 0.05474815 
    ## Run 7 stress 0.04851417 
    ## ... New best solution
    ## ... Procrustes: rmse 2.574732e-05  max resid 7.176035e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.04935948 
    ## Run 9 stress 0.04938616 
    ## Run 10 stress 0.05855031 
    ## Run 11 stress 0.05488104 
    ## Run 12 stress 0.05488105 
    ## Run 13 stress 0.04851417 
    ## ... Procrustes: rmse 1.231672e-05  max resid 3.248736e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.05607344 
    ## Run 15 stress 0.05855034 
    ## Run 16 stress 0.05474814 
    ## Run 17 stress 0.06205753 
    ## Run 18 stress 0.05728053 
    ## Run 19 stress 0.04935943 
    ## Run 20 stress 0.05607344 
    ## Run 21 stress 0.04851417 
    ## ... Procrustes: rmse 1.462688e-05  max resid 3.952575e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.05612239 
    ## Run 23 stress 0.05474818 
    ## Run 24 stress 0.05607348 
    ## Run 25 stress 0.04935945 
    ## Run 26 stress 0.05842785 
    ## Run 27 stress 0.05612241 
    ## Run 28 stress 0.05959895 
    ## Run 29 stress 0.06205752 
    ## Run 30 stress 0.04612758 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04983686  max resid 0.1395558 
    ## Run 31 stress 0.05216157 
    ## Run 32 stress 0.05842797 
    ## Run 33 stress 0.06205755 
    ## Run 34 stress 0.04938615 
    ## Run 35 stress 0.04938616 
    ## Run 36 stress 0.04612758 
    ## ... Procrustes: rmse 0.0001028098  max resid 0.0002017361 
    ## ... Similar to previous best
    ## Run 37 stress 0.05216149 
    ## Run 38 stress 0.05439282 
    ## Run 39 stress 0.04851422 
    ## Run 40 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 8.712277e-05  max resid 0.0001747801 
    ## ... Similar to previous best
    ## Run 41 stress 0.06185054 
    ## Run 42 stress 0.05643264 
    ## Run 43 stress 0.05605153 
    ## Run 44 stress 0.05260047 
    ## Run 45 stress 0.05605158 
    ## Run 46 stress 0.05439282 
    ## Run 47 stress 0.05433181 
    ## Run 48 stress 0.04612759 
    ## ... Procrustes: rmse 4.496124e-05  max resid 9.222653e-05 
    ## ... Similar to previous best
    ## Run 49 stress 0.05842821 
    ## Run 50 stress 0.04935947 
    ## Run 51 stress 0.04851418 
    ## Run 52 stress 0.04860423 
    ## Run 53 stress 0.05695873 
    ## Run 54 stress 0.05695884 
    ## Run 55 stress 0.04860422 
    ## Run 56 stress 0.05605154 
    ## Run 57 stress 0.05216129 
    ## Run 58 stress 0.04935943 
    ## Run 59 stress 0.05451493 
    ## Run 60 stress 0.06185047 
    ## Run 61 stress 0.05605156 
    ## Run 62 stress 0.04938616 
    ## Run 63 stress 0.05572688 
    ## Run 64 stress 0.04935945 
    ## Run 65 stress 0.0543928 
    ## Run 66 stress 0.05433179 
    ## Run 67 stress 0.05276415 
    ## Run 68 stress 0.0561224 
    ## Run 69 stress 0.04612765 
    ## ... Procrustes: rmse 9.995499e-05  max resid 0.0002002432 
    ## ... Similar to previous best
    ## Run 70 stress 0.04938617 
    ## Run 71 stress 0.04612757 
    ## ... New best solution
    ## ... Procrustes: rmse 4.48663e-05  max resid 9.02045e-05 
    ## ... Similar to previous best
    ## Run 72 stress 0.05605151 
    ## Run 73 stress 0.04938615 
    ## Run 74 stress 0.04851417 
    ## Run 75 stress 0.05607346 
    ## Run 76 stress 0.05607344 
    ## Run 77 stress 0.04938619 
    ## Run 78 stress 0.04851417 
    ## Run 79 stress 0.05276412 
    ## Run 80 stress 0.05612242 
    ## Run 81 stress 0.04612757 
    ## ... Procrustes: rmse 2.020045e-05  max resid 3.677396e-05 
    ## ... Similar to previous best
    ## Run 82 stress 0.04935945 
    ## Run 83 stress 0.06400756 
    ## Run 84 stress 0.05605146 
    ## Run 85 stress 0.05612244 
    ## Run 86 stress 0.05137964 
    ## Run 87 stress 0.05728065 
    ## Run 88 stress 0.04938616 
    ## Run 89 stress 0.04938616 
    ## Run 90 stress 0.06205753 
    ## Run 91 stress 0.05612239 
    ## Run 92 stress 0.04851417 
    ## Run 93 stress 0.05607346 
    ## Run 94 stress 0.04851417 
    ## Run 95 stress 0.04612758 
    ## ... Procrustes: rmse 3.828065e-05  max resid 7.951891e-05 
    ## ... Similar to previous best
    ## Run 96 stress 0.04938618 
    ## Run 97 stress 0.05605151 
    ## Run 98 stress 0.04851417 
    ## Run 99 stress 0.04938616 
    ## Run 100 stress 0.05474814 
    ## Run 101 stress 0.06185053 
    ## Run 102 stress 0.06185056 
    ## Run 103 stress 0.3578487 
    ## Run 104 stress 0.05643262 
    ## Run 105 stress 0.05438944 
    ## Run 106 stress 0.05855033 
    ## Run 107 stress 0.05572691 
    ## Run 108 stress 0.04851418 
    ## Run 109 stress 0.05474817 
    ## Run 110 stress 0.04860422 
    ## Run 111 stress 0.05451495 
    ## Run 112 stress 0.05474815 
    ## Run 113 stress 0.0527641 
    ## Run 114 stress 0.05728091 
    ## Run 115 stress 0.05643261 
    ## Run 116 stress 0.05842791 
    ## Run 117 stress 0.04860423 
    ## Run 118 stress 0.05695878 
    ## Run 119 stress 0.04860422 
    ## Run 120 stress 0.06185047 
    ## Run 121 stress 0.0543928 
    ## Run 122 stress 0.04938619 
    ## Run 123 stress 0.05216163 
    ## Run 124 stress 0.06205752 
    ## Run 125 stress 0.04851417 
    ## Run 126 stress 0.04938615 
    ## Run 127 stress 0.05855032 
    ## Run 128 stress 0.05474818 
    ## Run 129 stress 0.06088976 
    ## Run 130 stress 0.05572689 
    ## Run 131 stress 0.05643261 
    ## Run 132 stress 0.05137964 
    ## Run 133 stress 0.06205754 
    ## Run 134 stress 0.05842814 
    ## Run 135 stress 0.04851417 
    ## Run 136 stress 0.05593458 
    ## Run 137 stress 0.05451494 
    ## Run 138 stress 0.05643261 
    ## Run 139 stress 0.05855036 
    ## Run 140 stress 0.06088971 
    ## Run 141 stress 0.06291425 
    ## Run 142 stress 0.04851417 
    ## Run 143 stress 0.0585503 
    ## Run 144 stress 0.06291428 
    ## Run 145 stress 0.05607345 
    ## Run 146 stress 0.04938617 
    ## Run 147 stress 0.04935944 
    ## Run 148 stress 0.05728052 
    ## Run 149 stress 0.05488105 
    ## Run 150 stress 0.05728059 
    ## Run 151 stress 0.06400725 
    ## Run 152 stress 0.05855037 
    ## Run 153 stress 0.04935949 
    ## Run 154 stress 0.05643261 
    ## Run 155 stress 0.05607344 
    ## Run 156 stress 0.05276411 
    ## Run 157 stress 0.05607348 
    ## Run 158 stress 0.05137969 
    ## Run 159 stress 0.05612242 
    ## Run 160 stress 0.05643267 
    ## Run 161 stress 0.04612758 
    ## ... Procrustes: rmse 3.798459e-05  max resid 7.259442e-05 
    ## ... Similar to previous best
    ## Run 162 stress 0.05137967 
    ## Run 163 stress 0.04851418 
    ## Run 164 stress 0.05474816 
    ## Run 165 stress 0.04851417 
    ## Run 166 stress 0.05728065 
    ## Run 167 stress 0.05643262 
    ## Run 168 stress 0.0493862 
    ## Run 169 stress 0.05695862 
    ## Run 170 stress 0.05216134 
    ## Run 171 stress 0.04612758 
    ## ... Procrustes: rmse 4.52771e-05  max resid 9.148097e-05 
    ## ... Similar to previous best
    ## Run 172 stress 0.06291436 
    ## Run 173 stress 0.04938618 
    ## Run 174 stress 0.06291426 
    ## Run 175 stress 0.05593451 
    ## Run 176 stress 0.05572688 
    ## Run 177 stress 0.05474821 
    ## Run 178 stress 0.05438953 
    ## Run 179 stress 0.0585504 
    ## Run 180 stress 0.04851419 
    ## Run 181 stress 0.0521613 
    ## Run 182 stress 0.05451497 
    ## Run 183 stress 0.05572689 
    ## Run 184 stress 0.05451494 
    ## Run 185 stress 0.0527641 
    ## Run 186 stress 0.05643261 
    ## Run 187 stress 0.05612251 
    ## Run 188 stress 0.04851418 
    ## Run 189 stress 0.04612772 
    ## ... Procrustes: rmse 0.0001949145  max resid 0.0003938589 
    ## ... Similar to previous best
    ## Run 190 stress 0.04612757 
    ## ... Procrustes: rmse 7.792063e-06  max resid 1.461154e-05 
    ## ... Similar to previous best
    ## Run 191 stress 0.06291425 
    ## Run 192 stress 0.0640074 
    ## Run 193 stress 0.06291429 
    ## Run 194 stress 0.06291425 
    ## Run 195 stress 0.06291426 
    ## Run 196 stress 0.05855031 
    ## Run 197 stress 0.05643264 
    ## Run 198 stress 0.04935947 
    ## Run 199 stress 0.05605147 
    ## Run 200 stress 0.04938616 
    ## Run 201 stress 0.05695882 
    ## Run 202 stress 0.0461276 
    ## ... Procrustes: rmse 6.55047e-05  max resid 0.000129054 
    ## ... Similar to previous best
    ## Run 203 stress 0.05260048 
    ## Run 204 stress 0.04851417 
    ## Run 205 stress 0.04851417 
    ## Run 206 stress 0.06205755 
    ## Run 207 stress 0.04851418 
    ## Run 208 stress 0.05643261 
    ## Run 209 stress 0.05643264 
    ## Run 210 stress 0.04851417 
    ## Run 211 stress 0.05276415 
    ## Run 212 stress 0.0527641 
    ## Run 213 stress 0.04938616 
    ## Run 214 stress 0.05607345 
    ## Run 215 stress 0.04851417 
    ## Run 216 stress 0.05447815 
    ## Run 217 stress 0.0527641 
    ## Run 218 stress 0.05216138 
    ## Run 219 stress 0.05855031 
    ## Run 220 stress 0.05439281 
    ## Run 221 stress 0.04935951 
    ## Run 222 stress 0.05612249 
    ## Run 223 stress 0.05728067 
    ## Run 224 stress 0.04612761 
    ## ... Procrustes: rmse 9.04851e-05  max resid 0.0001828309 
    ## ... Similar to previous best
    ## Run 225 stress 0.04938615 
    ## Run 226 stress 0.05855035 
    ## Run 227 stress 0.06933819 
    ## Run 228 stress 0.05855029 
    ## Run 229 stress 0.05451498 
    ## Run 230 stress 0.05728074 
    ## Run 231 stress 0.05643263 
    ## Run 232 stress 0.04860423 
    ## Run 233 stress 0.05695881 
    ## Run 234 stress 0.05612242 
    ## Run 235 stress 0.05593456 
    ## Run 236 stress 0.05260046 
    ## Run 237 stress 0.04851418 
    ## Run 238 stress 0.06205752 
    ## Run 239 stress 0.0493862 
    ## Run 240 stress 0.05276411 
    ## Run 241 stress 0.05438939 
    ## Run 242 stress 0.05959897 
    ## Run 243 stress 0.04612757 
    ## ... Procrustes: rmse 2.635573e-05  max resid 5.175277e-05 
    ## ... Similar to previous best
    ## Run 244 stress 0.05695863 
    ## Run 245 stress 0.0618505 
    ## Run 246 stress 0.04851417 
    ## Run 247 stress 0.05607349 
    ## Run 248 stress 0.04860423 
    ## Run 249 stress 0.05474813 
    ## Run 250 stress 0.0526005 
    ## Run 251 stress 0.05643261 
    ## Run 252 stress 0.04935944 
    ## Run 253 stress 0.05438941 
    ## Run 254 stress 0.04612757 
    ## ... Procrustes: rmse 1.147679e-05  max resid 1.936743e-05 
    ## ... Similar to previous best
    ## Run 255 stress 0.05643263 
    ## Run 256 stress 0.04851417 
    ## Run 257 stress 0.05605148 
    ## Run 258 stress 0.04851417 
    ## Run 259 stress 0.05728063 
    ## Run 260 stress 0.05605148 
    ## Run 261 stress 0.05216144 
    ## Run 262 stress 0.0585503 
    ## Run 263 stress 0.05605148 
    ## Run 264 stress 0.05474816 
    ## Run 265 stress 0.04938617 
    ## Run 266 stress 0.04938616 
    ## Run 267 stress 0.05728053 
    ## Run 268 stress 0.04851418 
    ## Run 269 stress 0.04935949 
    ## Run 270 stress 0.0543928 
    ## Run 271 stress 0.05959896 
    ## Run 272 stress 0.3495727 
    ## Run 273 stress 0.06205752 
    ## Run 274 stress 0.05474825 
    ## Run 275 stress 0.05474819 
    ## Run 276 stress 0.05605149 
    ## Run 277 stress 0.05695884 
    ## Run 278 stress 0.06291425 
    ## Run 279 stress 0.04938616 
    ## Run 280 stress 0.05605148 
    ## Run 281 stress 0.04851417 
    ## Run 282 stress 0.0569586 
    ## Run 283 stress 0.05474816 
    ## Run 284 stress 0.04851417 
    ## Run 285 stress 0.05842787 
    ## Run 286 stress 0.05643264 
    ## Run 287 stress 0.04851417 
    ## Run 288 stress 0.05216124 
    ## Run 289 stress 0.04938616 
    ## Run 290 stress 0.05643263 
    ## Run 291 stress 0.0543928 
    ## Run 292 stress 0.05607344 
    ## Run 293 stress 0.05728082 
    ## Run 294 stress 0.06205759 
    ## Run 295 stress 0.04851417 
    ## Run 296 stress 0.04851418 
    ## Run 297 stress 0.05216144 
    ## Run 298 stress 0.05643262 
    ## Run 299 stress 0.05593465 
    ## Run 300 stress 0.05607344 
    ## Run 301 stress 0.04612762 
    ## ... Procrustes: rmse 0.0001248112  max resid 0.0002538104 
    ## ... Similar to previous best
    ## Run 302 stress 0.05607345 
    ## Run 303 stress 0.0493862 
    ## Run 304 stress 0.05643262 
    ## Run 305 stress 0.04851417 
    ## Run 306 stress 0.05216171 
    ## Run 307 stress 0.05451493 
    ## Run 308 stress 0.06291425 
    ## Run 309 stress 0.05474814 
    ## Run 310 stress 0.04851417 
    ## Run 311 stress 0.04851417 
    ## Run 312 stress 0.04851417 
    ## Run 313 stress 0.04851417 
    ## Run 314 stress 0.3316845 
    ## Run 315 stress 0.0493862 
    ## Run 316 stress 0.04938616 
    ## Run 317 stress 0.05137966 
    ## Run 318 stress 0.04938617 
    ## Run 319 stress 0.05593452 
    ## Run 320 stress 0.05451495 
    ## Run 321 stress 0.05855029 
    ## Run 322 stress 0.04851418 
    ## Run 323 stress 0.05607344 
    ## Run 324 stress 0.04935945 
    ## Run 325 stress 0.05260047 
    ## Run 326 stress 0.06185044 
    ## Run 327 stress 0.04938615 
    ## Run 328 stress 0.05607345 
    ## Run 329 stress 0.05474825 
    ## Run 330 stress 0.0572806 
    ## Run 331 stress 0.05728062 
    ## Run 332 stress 0.05605148 
    ## Run 333 stress 0.04938616 
    ## Run 334 stress 0.05728054 
    ## Run 335 stress 0.05855031 
    ## Run 336 stress 0.05643264 
    ## Run 337 stress 0.05855033 
    ## Run 338 stress 0.04851418 
    ## Run 339 stress 0.05216138 
    ## Run 340 stress 0.05605152 
    ## Run 341 stress 0.04938617 
    ## Run 342 stress 0.05605146 
    ## Run 343 stress 0.05438945 
    ## Run 344 stress 0.05474817 
    ## Run 345 stress 0.06291426 
    ## Run 346 stress 0.05728057 
    ## Run 347 stress 0.05728053 
    ## Run 348 stress 0.06400731 
    ## Run 349 stress 0.04860422 
    ## Run 350 stress 0.05605156 
    ## Run 351 stress 0.06205755 
    ## Run 352 stress 0.05474816 
    ## Run 353 stress 0.05451493 
    ## Run 354 stress 0.05451494 
    ## Run 355 stress 0.04938618 
    ## Run 356 stress 0.04938618 
    ## Run 357 stress 0.05605151 
    ## Run 358 stress 0.0493862 
    ## Run 359 stress 0.05474822 
    ## Run 360 stress 0.05216135 
    ## Run 361 stress 0.05695868 
    ## Run 362 stress 0.06291425 
    ## Run 363 stress 0.05855033 
    ## Run 364 stress 0.05607345 
    ## Run 365 stress 0.05959897 
    ## Run 366 stress 0.05607344 
    ## Run 367 stress 0.05276411 
    ## Run 368 stress 0.05438946 
    ## Run 369 stress 0.04935943 
    ## Run 370 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001171927  max resid 0.0002370792 
    ## ... Similar to previous best
    ## Run 371 stress 0.04612758 
    ## ... Procrustes: rmse 4.122183e-05  max resid 8.370294e-05 
    ## ... Similar to previous best
    ## Run 372 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001228057  max resid 0.0002462822 
    ## ... Similar to previous best
    ## Run 373 stress 0.06291429 
    ## Run 374 stress 0.04938615 
    ## Run 375 stress 0.05959897 
    ## Run 376 stress 0.05137964 
    ## Run 377 stress 0.05728071 
    ## Run 378 stress 0.0493862 
    ## Run 379 stress 0.05607345 
    ## Run 380 stress 0.05842798 
    ## Run 381 stress 0.05438941 
    ## Run 382 stress 0.05488108 
    ## Run 383 stress 0.05728052 
    ## Run 384 stress 0.06088971 
    ## Run 385 stress 0.04860423 
    ## Run 386 stress 0.05855032 
    ## Run 387 stress 0.05695878 
    ## Run 388 stress 0.05474813 
    ## Run 389 stress 0.05643261 
    ## Run 390 stress 0.04860423 
    ## Run 391 stress 0.05439284 
    ## Run 392 stress 0.04851417 
    ## Run 393 stress 0.04612763 
    ## ... Procrustes: rmse 0.0001255068  max resid 0.0002522773 
    ## ... Similar to previous best
    ## Run 394 stress 0.0560515 
    ## Run 395 stress 0.04612757 
    ## ... Procrustes: rmse 3.987281e-05  max resid 7.990331e-05 
    ## ... Similar to previous best
    ## Run 396 stress 0.0560515 
    ## Run 397 stress 0.0543928 
    ## Run 398 stress 0.05593449 
    ## Run 399 stress 0.04612759 
    ## ... Procrustes: rmse 7.404532e-05  max resid 0.0001509945 
    ## ... Similar to previous best
    ## Run 400 stress 0.05643262 
    ## Run 401 stress 0.05451493 
    ## Run 402 stress 0.04612757 
    ## ... Procrustes: rmse 2.758332e-05  max resid 4.583293e-05 
    ## ... Similar to previous best
    ## Run 403 stress 0.05572689 
    ## Run 404 stress 0.05728057 
    ## Run 405 stress 0.04851417 
    ## Run 406 stress 0.05612245 
    ## Run 407 stress 0.04938616 
    ## Run 408 stress 0.04938616 
    ## Run 409 stress 0.04612757 
    ## ... Procrustes: rmse 2.231666e-05  max resid 3.640025e-05 
    ## ... Similar to previous best
    ## Run 410 stress 0.04851417 
    ## Run 411 stress 0.04612764 
    ## ... Procrustes: rmse 0.0001184204  max resid 0.0002392037 
    ## ... Similar to previous best
    ## Run 412 stress 0.04851417 
    ## Run 413 stress 0.05438939 
    ## Run 414 stress 0.06933824 
    ## Run 415 stress 0.05451493 
    ## Run 416 stress 0.05855036 
    ## Run 417 stress 0.05607346 
    ## Run 418 stress 0.05276417 
    ## Run 419 stress 0.04851417 
    ## Run 420 stress 0.05643261 
    ## Run 421 stress 0.05643261 
    ## Run 422 stress 0.04851417 
    ## Run 423 stress 0.04851418 
    ## Run 424 stress 0.05216155 
    ## Run 425 stress 0.05593454 
    ## Run 426 stress 0.05572689 
    ## Run 427 stress 0.04935947 
    ## Run 428 stress 0.04935948 
    ## Run 429 stress 0.06291426 
    ## Run 430 stress 0.04938616 
    ## Run 431 stress 0.05572692 
    ## Run 432 stress 0.04851417 
    ## Run 433 stress 0.05474813 
    ## Run 434 stress 0.04938615 
    ## Run 435 stress 0.06291424 
    ## Run 436 stress 0.05474817 
    ## Run 437 stress 0.0543894 
    ## Run 438 stress 0.05728053 
    ## Run 439 stress 0.05607344 
    ## Run 440 stress 0.05605156 
    ## Run 441 stress 0.05643261 
    ## Run 442 stress 0.04938615 
    ## Run 443 stress 0.06185068 
    ## Run 444 stress 0.05607347 
    ## Run 445 stress 0.04851419 
    ## Run 446 stress 0.04612761 
    ## ... Procrustes: rmse 8.841703e-05  max resid 0.0001780853 
    ## ... Similar to previous best
    ## Run 447 stress 0.05855031 
    ## Run 448 stress 0.05842784 
    ## Run 449 stress 0.05728061 
    ## Run 450 stress 0.04860422 
    ## Run 451 stress 0.05643261 
    ## Run 452 stress 0.0493862 
    ## Run 453 stress 0.05451493 
    ## Run 454 stress 0.05643262 
    ## Run 455 stress 0.05572688 
    ## Run 456 stress 0.04935944 
    ## Run 457 stress 0.05959898 
    ## Run 458 stress 0.05855029 
    ## Run 459 stress 0.05695863 
    ## Run 460 stress 0.04938615 
    ## Run 461 stress 0.05439281 
    ## Run 462 stress 0.05433176 
    ## Run 463 stress 0.04938615 
    ## Run 464 stress 0.05137968 
    ## Run 465 stress 0.04860423 
    ## Run 466 stress 0.05572688 
    ## Run 467 stress 0.05605148 
    ## Run 468 stress 0.06400726 
    ## Run 469 stress 0.04938617 
    ## Run 470 stress 0.05593465 
    ## Run 471 stress 0.05593464 
    ## Run 472 stress 0.05572689 
    ## Run 473 stress 0.05728065 
    ## Run 474 stress 0.05643261 
    ## Run 475 stress 0.05607344 
    ## Run 476 stress 0.05855032 
    ## Run 477 stress 0.05607347 
    ## Run 478 stress 0.04851418 
    ## Run 479 stress 0.04938615 
    ## Run 480 stress 0.04938616 
    ## Run 481 stress 0.05607344 
    ## Run 482 stress 0.05488103 
    ## Run 483 stress 0.05728058 
    ## Run 484 stress 0.05607345 
    ## Run 485 stress 0.04851417 
    ## Run 486 stress 0.05216154 
    ## Run 487 stress 0.05643262 
    ## Run 488 stress 0.05438955 
    ## Run 489 stress 0.04851417 
    ## Run 490 stress 0.06400746 
    ## Run 491 stress 0.05438948 
    ## Run 492 stress 0.04935948 
    ## Run 493 stress 0.05474813 
    ## Run 494 stress 0.04851417 
    ## Run 495 stress 0.04938617 
    ## Run 496 stress 0.0560735 
    ## Run 497 stress 0.05607345 
    ## Run 498 stress 0.05433186 
    ## Run 499 stress 0.06185044 
    ## Run 500 stress 0.04851417 
    ## *** Best solution repeated 22 times

``` r
# Mixed lakes
PD_beta_geo_M_NMDS <- metaMDS(PD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.019724e-05 
    ## Run 1 stress 9.924721e-05 
    ## ... Procrustes: rmse 2.487167e-05  max resid 4.042574e-05 
    ## ... Similar to previous best
    ## Run 2 stress 9.386387e-05 
    ## ... Procrustes: rmse 0.0001933354  max resid 0.0003616781 
    ## ... Similar to previous best
    ## Run 3 stress 9.676459e-05 
    ## ... Procrustes: rmse 0.0001768162  max resid 0.0003595335 
    ## ... Similar to previous best
    ## Run 4 stress 9.088985e-05 
    ## ... Procrustes: rmse 0.0001547113  max resid 0.0002771397 
    ## ... Similar to previous best
    ## Run 5 stress 9.927096e-05 
    ## ... Procrustes: rmse 0.0002603645  max resid 0.0006019473 
    ## ... Similar to previous best
    ## Run 6 stress 9.567094e-05 
    ## ... Procrustes: rmse 0.0001745846  max resid 0.0003549862 
    ## ... Similar to previous best
    ## Run 7 stress 9.91026e-05 
    ## ... Procrustes: rmse 2.334797e-05  max resid 4.6454e-05 
    ## ... Similar to previous best
    ## Run 8 stress 8.62003e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001164361  max resid 0.0002351157 
    ## ... Similar to previous best
    ## Run 9 stress 8.756065e-05 
    ## ... Procrustes: rmse 0.0001724686  max resid 0.0002695414 
    ## ... Similar to previous best
    ## Run 10 stress 0.3028823 
    ## Run 11 stress 9.145965e-05 
    ## ... Procrustes: rmse 0.000155881  max resid 0.0003423826 
    ## ... Similar to previous best
    ## Run 12 stress 9.0253e-05 
    ## ... Procrustes: rmse 0.0001969893  max resid 0.0003032619 
    ## ... Similar to previous best
    ## Run 13 stress 8.956641e-05 
    ## ... Procrustes: rmse 0.0001331926  max resid 0.0003092056 
    ## ... Similar to previous best
    ## Run 14 stress 9.76214e-05 
    ## ... Procrustes: rmse 9.869895e-05  max resid 0.0001942286 
    ## ... Similar to previous best
    ## Run 15 stress 9.741495e-05 
    ## ... Procrustes: rmse 2.621132e-05  max resid 3.812977e-05 
    ## ... Similar to previous best
    ## Run 16 stress 9.295741e-05 
    ## ... Procrustes: rmse 0.0001728214  max resid 0.0003696361 
    ## ... Similar to previous best
    ## Run 17 stress 8.892554e-05 
    ## ... Procrustes: rmse 0.0001290278  max resid 0.0002849422 
    ## ... Similar to previous best
    ## Run 18 stress 9.486294e-05 
    ## ... Procrustes: rmse 0.0001240416  max resid 0.0002413443 
    ## ... Similar to previous best
    ## Run 19 stress 9.555553e-05 
    ## ... Procrustes: rmse 0.0002043699  max resid 0.0003566543 
    ## ... Similar to previous best
    ## Run 20 stress 9.266658e-05 
    ## ... Procrustes: rmse 0.0001311014  max resid 0.0002886705 
    ## ... Similar to previous best
    ## Run 21 stress 9.380114e-05 
    ## ... Procrustes: rmse 3.096893e-05  max resid 5.624543e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.2775739 
    ## Run 23 stress 6.841475e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000133979  max resid 0.0002858029 
    ## ... Similar to previous best
    ## Run 24 stress 9.896597e-05 
    ## ... Procrustes: rmse 0.0002026565  max resid 0.0002947086 
    ## ... Similar to previous best
    ## Run 25 stress 8.907241e-05 
    ## ... Procrustes: rmse 0.0001426909  max resid 0.000293314 
    ## ... Similar to previous best
    ## Run 26 stress 9.213771e-05 
    ## ... Procrustes: rmse 0.0001312288  max resid 0.0002152248 
    ## ... Similar to previous best
    ## Run 27 stress 9.241229e-05 
    ## ... Procrustes: rmse 0.0001750354  max resid 0.0002870338 
    ## ... Similar to previous best
    ## Run 28 stress 8.648207e-05 
    ## ... Procrustes: rmse 0.0001145251  max resid 0.0001939957 
    ## ... Similar to previous best
    ## Run 29 stress 8.690452e-05 
    ## ... Procrustes: rmse 0.0001348772  max resid 0.0002746266 
    ## ... Similar to previous best
    ## Run 30 stress 9.195809e-05 
    ## ... Procrustes: rmse 0.0001659959  max resid 0.0003690332 
    ## ... Similar to previous best
    ## Run 31 stress 9.366656e-05 
    ## ... Procrustes: rmse 0.0001243717  max resid 0.0002826093 
    ## ... Similar to previous best
    ## Run 32 stress 9.577848e-05 
    ## ... Procrustes: rmse 0.0001293188  max resid 0.0002110884 
    ## ... Similar to previous best
    ## Run 33 stress 0.2680462 
    ## Run 34 stress 9.749379e-05 
    ## ... Procrustes: rmse 0.0001473655  max resid 0.000240172 
    ## ... Similar to previous best
    ## Run 35 stress 7.312205e-05 
    ## ... Procrustes: rmse 0.0001321102  max resid 0.0002316875 
    ## ... Similar to previous best
    ## Run 36 stress 9.594072e-05 
    ## ... Procrustes: rmse 0.0001510385  max resid 0.0003132104 
    ## ... Similar to previous best
    ## Run 37 stress 8.440472e-05 
    ## ... Procrustes: rmse 9.320585e-05  max resid 0.000211339 
    ## ... Similar to previous best
    ## Run 38 stress 9.528299e-05 
    ## ... Procrustes: rmse 0.0001784674  max resid 0.0002923737 
    ## ... Similar to previous best
    ## Run 39 stress 9.225384e-05 
    ## ... Procrustes: rmse 0.0001596647  max resid 0.0002619421 
    ## ... Similar to previous best
    ## Run 40 stress 9.823287e-05 
    ## ... Procrustes: rmse 0.0002059678  max resid 0.000302317 
    ## ... Similar to previous best
    ## Run 41 stress 9.168344e-05 
    ## ... Procrustes: rmse 0.0001513853  max resid 0.0002922995 
    ## ... Similar to previous best
    ## Run 42 stress 7.911015e-05 
    ## ... Procrustes: rmse 0.0001012987  max resid 0.0001869883 
    ## ... Similar to previous best
    ## Run 43 stress 8.283724e-05 
    ## ... Procrustes: rmse 0.0001231499  max resid 0.0002117492 
    ## ... Similar to previous best
    ## Run 44 stress 9.967678e-05 
    ## ... Procrustes: rmse 0.0002100051  max resid 0.0003060562 
    ## ... Similar to previous best
    ## Run 45 stress 9.664955e-05 
    ## ... Procrustes: rmse 0.0001399723  max resid 0.0002718333 
    ## ... Similar to previous best
    ## Run 46 stress 9.09752e-05 
    ## ... Procrustes: rmse 0.0001238288  max resid 0.0002640201 
    ## ... Similar to previous best
    ## Run 47 stress 8.557232e-05 
    ## ... Procrustes: rmse 0.0001512189  max resid 0.0002428153 
    ## ... Similar to previous best
    ## Run 48 stress 9.065753e-05 
    ## ... Procrustes: rmse 9.265229e-05  max resid 0.0002038532 
    ## ... Similar to previous best
    ## Run 49 stress 6.058916e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001098112  max resid 0.00019978 
    ## ... Similar to previous best
    ## Run 50 stress 9.225903e-05 
    ## ... Procrustes: rmse 0.0001633038  max resid 0.0003262103 
    ## ... Similar to previous best
    ## Run 51 stress 9.584771e-05 
    ## ... Procrustes: rmse 0.0001716682  max resid 0.0003432533 
    ## ... Similar to previous best
    ## Run 52 stress 9.111465e-05 
    ## ... Procrustes: rmse 0.0001594924  max resid 0.0003193571 
    ## ... Similar to previous best
    ## Run 53 stress 9.379482e-05 
    ## ... Procrustes: rmse 0.0002203353  max resid 0.0003469671 
    ## ... Similar to previous best
    ## Run 54 stress 9.750564e-05 
    ## ... Procrustes: rmse 0.0002246085  max resid 0.0003639386 
    ## ... Similar to previous best
    ## Run 55 stress 9.398177e-05 
    ## ... Procrustes: rmse 0.0002220031  max resid 0.0003556851 
    ## ... Similar to previous best
    ## Run 56 stress 9.236902e-05 
    ## ... Procrustes: rmse 0.0002339639  max resid 0.0003440787 
    ## ... Similar to previous best
    ## Run 57 stress 0.3075829 
    ## Run 58 stress 9.684066e-05 
    ## ... Procrustes: rmse 0.0002279586  max resid 0.0003671062 
    ## ... Similar to previous best
    ## Run 59 stress 9.973857e-05 
    ## ... Procrustes: rmse 0.0001128532  max resid 0.0001924967 
    ## ... Similar to previous best
    ## Run 60 stress 9.934315e-05 
    ## ... Procrustes: rmse 0.0001687835  max resid 0.0003076348 
    ## ... Similar to previous best
    ## Run 61 stress 8.842505e-05 
    ## ... Procrustes: rmse 0.0001016387  max resid 0.0001674791 
    ## ... Similar to previous best
    ## Run 62 stress 0.3083098 
    ## Run 63 stress 9.613088e-05 
    ## ... Procrustes: rmse 0.0002511796  max resid 0.0005151161 
    ## ... Similar to previous best
    ## Run 64 stress 9.964408e-05 
    ## ... Procrustes: rmse 0.0002323676  max resid 0.0003621417 
    ## ... Similar to previous best
    ## Run 65 stress 9.87118e-05 
    ## ... Procrustes: rmse 0.0001395505  max resid 0.0002221285 
    ## ... Similar to previous best
    ## Run 66 stress 9.811507e-05 
    ## ... Procrustes: rmse 0.0001075658  max resid 0.000148379 
    ## ... Similar to previous best
    ## Run 67 stress 9.557384e-05 
    ## ... Procrustes: rmse 0.0001736243  max resid 0.0003370479 
    ## ... Similar to previous best
    ## Run 68 stress 0.2842805 
    ## Run 69 stress 9.101249e-05 
    ## ... Procrustes: rmse 0.0001165947  max resid 0.0001907707 
    ## ... Similar to previous best
    ## Run 70 stress 9.598807e-05 
    ## ... Procrustes: rmse 0.0002507309  max resid 0.0003504294 
    ## ... Similar to previous best
    ## Run 71 stress 9.647795e-05 
    ## ... Procrustes: rmse 0.0002028936  max resid 0.0003865044 
    ## ... Similar to previous best
    ## Run 72 stress 9.509991e-05 
    ## ... Procrustes: rmse 0.0002177958  max resid 0.0003540962 
    ## ... Similar to previous best
    ## Run 73 stress 9.815316e-05 
    ## ... Procrustes: rmse 0.0002293913  max resid 0.000373334 
    ## ... Similar to previous best
    ## Run 74 stress 9.839368e-05 
    ## ... Procrustes: rmse 0.0002594412  max resid 0.0003694822 
    ## ... Similar to previous best
    ## Run 75 stress 9.489129e-05 
    ## ... Procrustes: rmse 0.0002500803  max resid 0.0005185104 
    ## ... Similar to previous best
    ## Run 76 stress 9.282476e-05 
    ## ... Procrustes: rmse 0.0001966156  max resid 0.0003633094 
    ## ... Similar to previous best
    ## Run 77 stress 8.523097e-05 
    ## ... Procrustes: rmse 0.0002201745  max resid 0.0004653516 
    ## ... Similar to previous best
    ## Run 78 stress 0.3083098 
    ## Run 79 stress 0.3079616 
    ## Run 80 stress 9.663075e-05 
    ## ... Procrustes: rmse 0.0001974635  max resid 0.0002987606 
    ## ... Similar to previous best
    ## Run 81 stress 9.847235e-05 
    ## ... Procrustes: rmse 0.0001640585  max resid 0.0003079748 
    ## ... Similar to previous best
    ## Run 82 stress 9.662238e-05 
    ## ... Procrustes: rmse 0.0002082064  max resid 0.0003159111 
    ## ... Similar to previous best
    ## Run 83 stress 9.95171e-05 
    ## ... Procrustes: rmse 0.0002249838  max resid 0.0003547042 
    ## ... Similar to previous best
    ## Run 84 stress 9.773775e-05 
    ## ... Procrustes: rmse 0.0002263848  max resid 0.0003556748 
    ## ... Similar to previous best
    ## Run 85 stress 9.704405e-05 
    ## ... Procrustes: rmse 9.229938e-05  max resid 0.0001951421 
    ## ... Similar to previous best
    ## Run 86 stress 9.209959e-05 
    ## ... Procrustes: rmse 0.0001272248  max resid 0.0002063168 
    ## ... Similar to previous best
    ## Run 87 stress 8.136696e-05 
    ## ... Procrustes: rmse 0.0002142474  max resid 0.0004497174 
    ## ... Similar to previous best
    ## Run 88 stress 9.231275e-05 
    ## ... Procrustes: rmse 0.0001681127  max resid 0.0003349377 
    ## ... Similar to previous best
    ## Run 89 stress 9.918641e-05 
    ## ... Procrustes: rmse 8.647398e-05  max resid 0.0001480877 
    ## ... Similar to previous best
    ## Run 90 stress 9.748098e-05 
    ## ... Procrustes: rmse 0.0002559506  max resid 0.000366542 
    ## ... Similar to previous best
    ## Run 91 stress 9.686368e-05 
    ## ... Procrustes: rmse 9.485453e-05  max resid 0.00019727 
    ## ... Similar to previous best
    ## Run 92 stress 9.880414e-05 
    ## ... Procrustes: rmse 0.0002566222  max resid 0.0005250389 
    ## ... Similar to previous best
    ## Run 93 stress 9.568681e-05 
    ## ... Procrustes: rmse 0.0002515734  max resid 0.0005188682 
    ## ... Similar to previous best
    ## Run 94 stress 8.74554e-05 
    ## ... Procrustes: rmse 0.0002314782  max resid 0.0003316006 
    ## ... Similar to previous best
    ## Run 95 stress 9.10416e-05 
    ## ... Procrustes: rmse 0.0001261093  max resid 0.0002205937 
    ## ... Similar to previous best
    ## Run 96 stress 9.62597e-05 
    ## ... Procrustes: rmse 0.0002463518  max resid 0.0005148676 
    ## ... Similar to previous best
    ## Run 97 stress 9.747821e-05 
    ## ... Procrustes: rmse 0.0002515957  max resid 0.0005140799 
    ## ... Similar to previous best
    ## Run 98 stress 9.522774e-05 
    ## ... Procrustes: rmse 0.0002489065  max resid 0.0003467766 
    ## ... Similar to previous best
    ## Run 99 stress 9.251823e-05 
    ## ... Procrustes: rmse 0.0002405575  max resid 0.0005004669 
    ## ... Similar to previous best
    ## Run 100 stress 9.948338e-05 
    ## ... Procrustes: rmse 0.0001653815  max resid 0.0003222832 
    ## ... Similar to previous best
    ## Run 101 stress 8.821154e-05 
    ## ... Procrustes: rmse 0.0001355162  max resid 0.0002477806 
    ## ... Similar to previous best
    ## Run 102 stress 9.94365e-05 
    ## ... Procrustes: rmse 0.0002279284  max resid 0.0003588859 
    ## ... Similar to previous best
    ## Run 103 stress 9.700959e-05 
    ## ... Procrustes: rmse 0.0002545152  max resid 0.0003697891 
    ## ... Similar to previous best
    ## Run 104 stress 8.825487e-05 
    ## ... Procrustes: rmse 0.000149948  max resid 0.0003063575 
    ## ... Similar to previous best
    ## Run 105 stress 9.772605e-05 
    ## ... Procrustes: rmse 0.0002572907  max resid 0.0005291649 
    ## ... Similar to previous best
    ## Run 106 stress 9.710061e-05 
    ## ... Procrustes: rmse 0.0002555559  max resid 0.0005285262 
    ## ... Similar to previous best
    ## Run 107 stress 8.770254e-05 
    ## ... Procrustes: rmse 0.0002301864  max resid 0.0004810279 
    ## ... Similar to previous best
    ## Run 108 stress 9.316404e-05 
    ## ... Procrustes: rmse 0.0001679573  max resid 0.0003354283 
    ## ... Similar to previous best
    ## Run 109 stress 7.838591e-05 
    ## ... Procrustes: rmse 0.0001066442  max resid 0.0001455584 
    ## ... Similar to previous best
    ## Run 110 stress 9.104351e-05 
    ## ... Procrustes: rmse 0.0001291125  max resid 0.0001743525 
    ## ... Similar to previous best
    ## Run 111 stress 9.944097e-05 
    ## ... Procrustes: rmse 0.0002316834  max resid 0.0003748672 
    ## ... Similar to previous best
    ## Run 112 stress 9.51998e-05 
    ## ... Procrustes: rmse 0.000222282  max resid 0.0003470828 
    ## ... Similar to previous best
    ## Run 113 stress 9.330786e-05 
    ## ... Procrustes: rmse 0.0002207239  max resid 0.0003313352 
    ## ... Similar to previous best
    ## Run 114 stress 8.922975e-05 
    ## ... Procrustes: rmse 0.0001061133  max resid 0.0001727138 
    ## ... Similar to previous best
    ## Run 115 stress 8.954914e-05 
    ## ... Procrustes: rmse 0.0001049112  max resid 0.0001708346 
    ## ... Similar to previous best
    ## Run 116 stress 7.821282e-05 
    ## ... Procrustes: rmse 9.105069e-05  max resid 0.0001355848 
    ## ... Similar to previous best
    ## Run 117 stress 9.507914e-05 
    ## ... Procrustes: rmse 0.0001738661  max resid 0.0003385777 
    ## ... Similar to previous best
    ## Run 118 stress 9.520168e-05 
    ## ... Procrustes: rmse 0.0002468849  max resid 0.0005070173 
    ## ... Similar to previous best
    ## Run 119 stress 9.677728e-05 
    ## ... Procrustes: rmse 0.0001157282  max resid 0.0001927617 
    ## ... Similar to previous best
    ## Run 120 stress 9.695206e-05 
    ## ... Procrustes: rmse 0.0002276065  max resid 0.0003669617 
    ## ... Similar to previous best
    ## Run 121 stress 9.95587e-05 
    ## ... Procrustes: rmse 0.0002321461  max resid 0.0003742473 
    ## ... Similar to previous best
    ## Run 122 stress 9.366214e-05 
    ## ... Procrustes: rmse 0.0002472342  max resid 0.0005109751 
    ## ... Similar to previous best
    ## Run 123 stress 9.961017e-05 
    ## ... Procrustes: rmse 0.0002604066  max resid 0.0003785504 
    ## ... Similar to previous best
    ## Run 124 stress 9.490424e-05 
    ## ... Procrustes: rmse 0.0001630637  max resid 0.0003268484 
    ## ... Similar to previous best
    ## Run 125 stress 9.932395e-05 
    ## ... Procrustes: rmse 0.0002608701  max resid 0.0003711228 
    ## ... Similar to previous best
    ## Run 126 stress 9.120497e-05 
    ## ... Procrustes: rmse 9.990647e-05  max resid 0.0001691532 
    ## ... Similar to previous best
    ## Run 127 stress 9.646666e-05 
    ## ... Procrustes: rmse 0.0001310594  max resid 0.0002225224 
    ## ... Similar to previous best
    ## Run 128 stress 0.3079294 
    ## Run 129 stress 9.379455e-05 
    ## ... Procrustes: rmse 0.0002454024  max resid 0.000504743 
    ## ... Similar to previous best
    ## Run 130 stress 9.936225e-05 
    ## ... Procrustes: rmse 0.0001786616  max resid 0.0003565171 
    ## ... Similar to previous best
    ## Run 131 stress 0.2842088 
    ## Run 132 stress 9.818383e-05 
    ## ... Procrustes: rmse 0.0001492156  max resid 0.0002648622 
    ## ... Similar to previous best
    ## Run 133 stress 9.415323e-05 
    ## ... Procrustes: rmse 0.0002455223  max resid 0.0003541253 
    ## ... Similar to previous best
    ## Run 134 stress 9.577317e-05 
    ## ... Procrustes: rmse 0.000174901  max resid 0.0003403826 
    ## ... Similar to previous best
    ## Run 135 stress 9.149624e-05 
    ## ... Procrustes: rmse 0.0001834202  max resid 0.000262054 
    ## ... Similar to previous best
    ## Run 136 stress 9.554484e-05 
    ## ... Procrustes: rmse 0.0001923434  max resid 0.0002999102 
    ## ... Similar to previous best
    ## Run 137 stress 8.932281e-05 
    ## ... Procrustes: rmse 0.0001301614  max resid 0.0002219363 
    ## ... Similar to previous best
    ## Run 138 stress 7.978924e-05 
    ## ... Procrustes: rmse 6.479617e-05  max resid 0.0001249115 
    ## ... Similar to previous best
    ## Run 139 stress 9.704427e-05 
    ## ... Procrustes: rmse 0.0002495309  max resid 0.0005196552 
    ## ... Similar to previous best
    ## Run 140 stress 9.755432e-05 
    ## ... Procrustes: rmse 0.0002487312  max resid 0.0005179295 
    ## ... Similar to previous best
    ## Run 141 stress 8.32573e-05 
    ## ... Procrustes: rmse 6.965864e-05  max resid 8.977718e-05 
    ## ... Similar to previous best
    ## Run 142 stress 9.438015e-05 
    ## ... Procrustes: rmse 0.0002454116  max resid 0.0005094279 
    ## ... Similar to previous best
    ## Run 143 stress 9.741045e-05 
    ## ... Procrustes: rmse 0.0002296248  max resid 0.0003617561 
    ## ... Similar to previous best
    ## Run 144 stress 9.917302e-05 
    ## ... Procrustes: rmse 0.0002335726  max resid 0.000366407 
    ## ... Similar to previous best
    ## Run 145 stress 9.981123e-05 
    ## ... Procrustes: rmse 0.0002155189  max resid 0.0003572504 
    ## ... Similar to previous best
    ## Run 146 stress 9.740835e-05 
    ## ... Procrustes: rmse 0.0001650297  max resid 0.0003152052 
    ## ... Similar to previous best
    ## Run 147 stress 0.2842805 
    ## Run 148 stress 9.919998e-05 
    ## ... Procrustes: rmse 0.0001807535  max resid 0.0003515288 
    ## ... Similar to previous best
    ## Run 149 stress 9.008898e-05 
    ## ... Procrustes: rmse 0.0001641613  max resid 0.0003263791 
    ## ... Similar to previous best
    ## Run 150 stress 9.449845e-05 
    ## ... Procrustes: rmse 0.0002400448  max resid 0.0003502009 
    ## ... Similar to previous best
    ## Run 151 stress 8.945454e-05 
    ## ... Procrustes: rmse 0.0001513877  max resid 0.0002869321 
    ## ... Similar to previous best
    ## Run 152 stress 9.735869e-05 
    ## ... Procrustes: rmse 0.000177626  max resid 0.0003455769 
    ## ... Similar to previous best
    ## Run 153 stress 7.443861e-05 
    ## ... Procrustes: rmse 8.359868e-05  max resid 0.0001296006 
    ## ... Similar to previous best
    ## Run 154 stress 9.224335e-05 
    ## ... Procrustes: rmse 0.0001705526  max resid 0.0003343013 
    ## ... Similar to previous best
    ## Run 155 stress 9.688876e-05 
    ## ... Procrustes: rmse 0.0002523593  max resid 0.0005167387 
    ## ... Similar to previous best
    ## Run 156 stress 9.409993e-05 
    ## ... Procrustes: rmse 0.0001901285  max resid 0.0002834652 
    ## ... Similar to previous best
    ## Run 157 stress 9.902327e-05 
    ## ... Procrustes: rmse 0.000258241  max resid 0.0003641033 
    ## ... Similar to previous best
    ## Run 158 stress 0.3081021 
    ## Run 159 stress 9.685077e-05 
    ## ... Procrustes: rmse 0.0001852196  max resid 0.0002895215 
    ## ... Similar to previous best
    ## Run 160 stress 9.825743e-05 
    ## ... Procrustes: rmse 0.0001674791  max resid 0.0003227289 
    ## ... Similar to previous best
    ## Run 161 stress 8.985483e-05 
    ## ... Procrustes: rmse 0.0001614817  max resid 0.0003229764 
    ## ... Similar to previous best
    ## Run 162 stress 0.3083098 
    ## Run 163 stress 0.2854403 
    ## Run 164 stress 9.503685e-05 
    ## ... Procrustes: rmse 0.0002460791  max resid 0.000361016 
    ## ... Similar to previous best
    ## Run 165 stress 9.538023e-05 
    ## ... Procrustes: rmse 0.0002486665  max resid 0.0005168917 
    ## ... Similar to previous best
    ## Run 166 stress 8.720681e-05 
    ## ... Procrustes: rmse 0.0001321668  max resid 0.0002182664 
    ## ... Similar to previous best
    ## Run 167 stress 9.855583e-05 
    ## ... Procrustes: rmse 0.0002574629  max resid 0.000526005 
    ## ... Similar to previous best
    ## Run 168 stress 9.054532e-05 
    ## ... Procrustes: rmse 0.0001252936  max resid 0.0002148077 
    ## ... Similar to previous best
    ## Run 169 stress 8.383843e-05 
    ## ... Procrustes: rmse 6.974504e-05  max resid 0.0001248174 
    ## ... Similar to previous best
    ## Run 170 stress 9.237307e-05 
    ## ... Procrustes: rmse 0.0001859936  max resid 0.0002829715 
    ## ... Similar to previous best
    ## Run 171 stress 9.698279e-05 
    ## ... Procrustes: rmse 0.000249087  max resid 0.0005096165 
    ## ... Similar to previous best
    ## Run 172 stress 9.252697e-05 
    ## ... Procrustes: rmse 0.0002430783  max resid 0.0005002588 
    ## ... Similar to previous best
    ## Run 173 stress 7.401892e-05 
    ## ... Procrustes: rmse 0.0001628693  max resid 0.0002650964 
    ## ... Similar to previous best
    ## Run 174 stress 9.038559e-05 
    ## ... Procrustes: rmse 0.0002255176  max resid 0.0003362212 
    ## ... Similar to previous best
    ## Run 175 stress 9.29379e-05 
    ## ... Procrustes: rmse 0.000167674  max resid 0.0003348347 
    ## ... Similar to previous best
    ## Run 176 stress 9.42942e-05 
    ## ... Procrustes: rmse 0.0001569852  max resid 0.0002894925 
    ## ... Similar to previous best
    ## Run 177 stress 9.610894e-05 
    ## ... Procrustes: rmse 0.0002083332  max resid 0.0003326802 
    ## ... Similar to previous best
    ## Run 178 stress 9.840409e-05 
    ## ... Procrustes: rmse 0.0002308831  max resid 0.0003715545 
    ## ... Similar to previous best
    ## Run 179 stress 8.555481e-05 
    ## ... Procrustes: rmse 9.854819e-05  max resid 0.0001591251 
    ## ... Similar to previous best
    ## Run 180 stress 8.811511e-05 
    ## ... Procrustes: rmse 0.0002336426  max resid 0.0003343229 
    ## ... Similar to previous best
    ## Run 181 stress 9.682391e-05 
    ## ... Procrustes: rmse 0.0001348032  max resid 0.0002158454 
    ## ... Similar to previous best
    ## Run 182 stress 9.131113e-05 
    ## ... Procrustes: rmse 0.0001667339  max resid 0.0003241423 
    ## ... Similar to previous best
    ## Run 183 stress 9.295652e-05 
    ## ... Procrustes: rmse 8.209177e-05  max resid 0.0001766711 
    ## ... Similar to previous best
    ## Run 184 stress 9.157403e-05 
    ## ... Procrustes: rmse 0.0001738782  max resid 0.0002531088 
    ## ... Similar to previous best
    ## Run 185 stress 9.52429e-05 
    ## ... Procrustes: rmse 0.0001309696  max resid 0.0001753016 
    ## ... Similar to previous best
    ## Run 186 stress 9.578633e-05 
    ## ... Procrustes: rmse 0.0002246604  max resid 0.0003620229 
    ## ... Similar to previous best
    ## Run 187 stress 9.43156e-05 
    ## ... Procrustes: rmse 0.0002473474  max resid 0.0003491209 
    ## ... Similar to previous best
    ## Run 188 stress 9.361685e-05 
    ## ... Procrustes: rmse 0.0002446243  max resid 0.000503171 
    ## ... Similar to previous best
    ## Run 189 stress 9.126684e-05 
    ## ... Procrustes: rmse 0.0002115714  max resid 0.0003448212 
    ## ... Similar to previous best
    ## Run 190 stress 8.229162e-05 
    ## ... Procrustes: rmse 0.0001187713  max resid 0.0002272604 
    ## ... Similar to previous best
    ## Run 191 stress 9.05844e-05 
    ## ... Procrustes: rmse 0.0002047916  max resid 0.0003177269 
    ## ... Similar to previous best
    ## Run 192 stress 0.3083098 
    ## Run 193 stress 9.447623e-05 
    ## ... Procrustes: rmse 0.0001329039  max resid 0.0001828597 
    ## ... Similar to previous best
    ## Run 194 stress 8.898452e-05 
    ## ... Procrustes: rmse 0.0002349486  max resid 0.0003382955 
    ## ... Similar to previous best
    ## Run 195 stress 8.720073e-05 
    ## ... Procrustes: rmse 7.516236e-05  max resid 0.000121729 
    ## ... Similar to previous best
    ## Run 196 stress 9.736579e-05 
    ## ... Procrustes: rmse 9.28146e-05  max resid 0.0001960622 
    ## ... Similar to previous best
    ## Run 197 stress 0.2848518 
    ## Run 198 stress 9.980304e-05 
    ## ... Procrustes: rmse 0.00026249  max resid 0.0005383368 
    ## ... Similar to previous best
    ## Run 199 stress 7.371613e-05 
    ## ... Procrustes: rmse 0.0001765557  max resid 0.0002468398 
    ## ... Similar to previous best
    ## Run 200 stress 8.217652e-05 
    ## ... Procrustes: rmse 0.0001516035  max resid 0.0003003172 
    ## ... Similar to previous best
    ## Run 201 stress 9.523662e-05 
    ## ... Procrustes: rmse 0.0002253741  max resid 0.0003547334 
    ## ... Similar to previous best
    ## Run 202 stress 9.933001e-05 
    ## ... Procrustes: rmse 0.0002329099  max resid 0.000367383 
    ## ... Similar to previous best
    ## Run 203 stress 8.679943e-05 
    ## ... Procrustes: rmse 7.436923e-05  max resid 0.0001059554 
    ## ... Similar to previous best
    ## Run 204 stress 9.549188e-05 
    ## ... Procrustes: rmse 0.0001953722  max resid 0.0003815245 
    ## ... Similar to previous best
    ## Run 205 stress 9.667842e-05 
    ## ... Procrustes: rmse 0.0002242486  max resid 0.0003507484 
    ## ... Similar to previous best
    ## Run 206 stress 9.881159e-05 
    ## ... Procrustes: rmse 0.0002338541  max resid 0.0003738466 
    ## ... Similar to previous best
    ## Run 207 stress 9.974326e-05 
    ## ... Procrustes: rmse 0.0002536513  max resid 0.0003713201 
    ## ... Similar to previous best
    ## Run 208 stress 0.2848523 
    ## Run 209 stress 0.3075828 
    ## Run 210 stress 9.141556e-05 
    ## ... Procrustes: rmse 0.0002382714  max resid 0.00033943 
    ## ... Similar to previous best
    ## Run 211 stress 0.3083098 
    ## Run 212 stress 9.901726e-05 
    ## ... Procrustes: rmse 0.0002524388  max resid 0.0003502018 
    ## ... Similar to previous best
    ## Run 213 stress 0.3083098 
    ## Run 214 stress 9.18836e-05 
    ## ... Procrustes: rmse 0.0002346036  max resid 0.0003380655 
    ## ... Similar to previous best
    ## Run 215 stress 9.595979e-05 
    ## ... Procrustes: rmse 8.726522e-05  max resid 0.0001875044 
    ## ... Similar to previous best
    ## Run 216 stress 9.003471e-05 
    ## ... Procrustes: rmse 0.0001236348  max resid 0.0002124916 
    ## ... Similar to previous best
    ## Run 217 stress 9.997297e-05 
    ## ... Procrustes: rmse 0.000177561  max resid 0.000341842 
    ## ... Similar to previous best
    ## Run 218 stress 9.138475e-05 
    ## ... Procrustes: rmse 0.0002321678  max resid 0.0003353807 
    ## ... Similar to previous best
    ## Run 219 stress 8.464108e-05 
    ## ... Procrustes: rmse 0.0001963155  max resid 0.0003130014 
    ## ... Similar to previous best
    ## Run 220 stress 9.699441e-05 
    ## ... Procrustes: rmse 0.0002235218  max resid 0.0003465229 
    ## ... Similar to previous best
    ## Run 221 stress 9.81275e-05 
    ## ... Procrustes: rmse 0.0002583875  max resid 0.0005315896 
    ## ... Similar to previous best
    ## Run 222 stress 9.17576e-05 
    ## ... Procrustes: rmse 0.0002365869  max resid 0.0004949771 
    ## ... Similar to previous best
    ## Run 223 stress 9.061962e-05 
    ## ... Procrustes: rmse 0.0002059958  max resid 0.0003268359 
    ## ... Similar to previous best
    ## Run 224 stress 8.379838e-05 
    ## ... Procrustes: rmse 0.0002166268  max resid 0.0003123397 
    ## ... Similar to previous best
    ## Run 225 stress 9.571683e-05 
    ## ... Procrustes: rmse 8.626561e-05  max resid 0.0001845612 
    ## ... Similar to previous best
    ## Run 226 stress 9.789557e-05 
    ## ... Procrustes: rmse 0.0001988968  max resid 0.0003085484 
    ## ... Similar to previous best
    ## Run 227 stress 9.478813e-05 
    ## ... Procrustes: rmse 0.0002498382  max resid 0.0005147827 
    ## ... Similar to previous best
    ## Run 228 stress 9.844674e-05 
    ## ... Procrustes: rmse 0.0002058642  max resid 0.000391923 
    ## ... Similar to previous best
    ## Run 229 stress 8.230382e-05 
    ## ... Procrustes: rmse 0.000187861  max resid 0.0002937523 
    ## ... Similar to previous best
    ## Run 230 stress 9.430762e-05 
    ## ... Procrustes: rmse 0.0001346572  max resid 0.0002196927 
    ## ... Similar to previous best
    ## Run 231 stress 9.67441e-05 
    ## ... Procrustes: rmse 0.0002528515  max resid 0.0003650276 
    ## ... Similar to previous best
    ## Run 232 stress 9.700939e-05 
    ## ... Procrustes: rmse 0.000171712  max resid 0.000345192 
    ## ... Similar to previous best
    ## Run 233 stress 9.742678e-05 
    ## ... Procrustes: rmse 0.0001765456  max resid 0.0003646403 
    ## ... Similar to previous best
    ## Run 234 stress 8.167978e-05 
    ## ... Procrustes: rmse 0.0001189501  max resid 0.0001941805 
    ## ... Similar to previous best
    ## Run 235 stress 8.947226e-05 
    ## ... Procrustes: rmse 0.0002352763  max resid 0.0003357354 
    ## ... Similar to previous best
    ## Run 236 stress 9.848833e-05 
    ## ... Procrustes: rmse 0.0002221146  max resid 0.0003448158 
    ## ... Similar to previous best
    ## Run 237 stress 9.776194e-05 
    ## ... Procrustes: rmse 0.0002206963  max resid 0.0003640633 
    ## ... Similar to previous best
    ## Run 238 stress 0.3079293 
    ## Run 239 stress 9.179723e-05 
    ## ... Procrustes: rmse 0.0001371873  max resid 0.0002231472 
    ## ... Similar to previous best
    ## Run 240 stress 9.477431e-05 
    ## ... Procrustes: rmse 0.0001674506  max resid 0.0003239585 
    ## ... Similar to previous best
    ## Run 241 stress 9.265055e-05 
    ## ... Procrustes: rmse 0.0002297774  max resid 0.0003283059 
    ## ... Similar to previous best
    ## Run 242 stress 9.944343e-05 
    ## ... Procrustes: rmse 0.0002311535  max resid 0.0003638785 
    ## ... Similar to previous best
    ## Run 243 stress 9.969539e-05 
    ## ... Procrustes: rmse 0.0001676483  max resid 0.0003141317 
    ## ... Similar to previous best
    ## Run 244 stress 9.569028e-05 
    ## ... Procrustes: rmse 0.0002458117  max resid 0.0005129854 
    ## ... Similar to previous best
    ## Run 245 stress 9.402691e-05 
    ## ... Procrustes: rmse 0.0001716776  max resid 0.0003419081 
    ## ... Similar to previous best
    ## Run 246 stress 9.529565e-05 
    ## ... Procrustes: rmse 0.0001926362  max resid 0.0002921096 
    ## ... Similar to previous best
    ## Run 247 stress 9.5536e-05 
    ## ... Procrustes: rmse 0.0001338538  max resid 0.0002224673 
    ## ... Similar to previous best
    ## Run 248 stress 8.770407e-05 
    ## ... Procrustes: rmse 7.968211e-05  max resid 0.0001740518 
    ## ... Similar to previous best
    ## Run 249 stress 9.944363e-05 
    ## ... Procrustes: rmse 0.0002267302  max resid 0.0003501503 
    ## ... Similar to previous best
    ## Run 250 stress 9.571561e-05 
    ## ... Procrustes: rmse 0.0001948345  max resid 0.0003037184 
    ## ... Similar to previous best
    ## Run 251 stress 9.229173e-05 
    ## ... Procrustes: rmse 0.0001858467  max resid 0.0003543851 
    ## ... Similar to previous best
    ## Run 252 stress 9.468238e-05 
    ## ... Procrustes: rmse 0.0001936969  max resid 0.0002993174 
    ## ... Similar to previous best
    ## Run 253 stress 0.3120129 
    ## Run 254 stress 9.190198e-05 
    ## ... Procrustes: rmse 0.000212832  max resid 0.0003315667 
    ## ... Similar to previous best
    ## Run 255 stress 9.463874e-05 
    ## ... Procrustes: rmse 0.0002177416  max resid 0.0003373752 
    ## ... Similar to previous best
    ## Run 256 stress 9.815052e-05 
    ## ... Procrustes: rmse 0.0002577798  max resid 0.0005312432 
    ## ... Similar to previous best
    ## Run 257 stress 7.254546e-05 
    ## ... Procrustes: rmse 0.00011275  max resid 0.0001984955 
    ## ... Similar to previous best
    ## Run 258 stress 9.633586e-05 
    ## ... Procrustes: rmse 0.0002479672  max resid 0.0003659974 
    ## ... Similar to previous best
    ## Run 259 stress 9.88145e-05 
    ## ... Procrustes: rmse 0.000231797  max resid 0.000370923 
    ## ... Similar to previous best
    ## Run 260 stress 9.513672e-05 
    ## ... Procrustes: rmse 0.0002471457  max resid 0.0003516413 
    ## ... Similar to previous best
    ## Run 261 stress 9.799588e-05 
    ## ... Procrustes: rmse 0.0002542991  max resid 0.000527295 
    ## ... Similar to previous best
    ## Run 262 stress 9.518952e-05 
    ## ... Procrustes: rmse 0.0002501734  max resid 0.0003523881 
    ## ... Similar to previous best
    ## Run 263 stress 9.746426e-05 
    ## ... Procrustes: rmse 0.0002494856  max resid 0.0005198587 
    ## ... Similar to previous best
    ## Run 264 stress 9.830567e-05 
    ## ... Procrustes: rmse 0.0001387284  max resid 0.0001937221 
    ## ... Similar to previous best
    ## Run 265 stress 9.947232e-05 
    ## ... Procrustes: rmse 0.0002345856  max resid 0.0003689611 
    ## ... Similar to previous best
    ## Run 266 stress 8.564722e-05 
    ## ... Procrustes: rmse 6.989768e-05  max resid 0.0001383291 
    ## ... Similar to previous best
    ## Run 267 stress 9.282366e-05 
    ## ... Procrustes: rmse 0.00024497  max resid 0.0005069469 
    ## ... Similar to previous best
    ## Run 268 stress 9.901312e-05 
    ## ... Procrustes: rmse 0.0001698773  max resid 0.0003418758 
    ## ... Similar to previous best
    ## Run 269 stress 9.252623e-05 
    ## ... Procrustes: rmse 0.0001005038  max resid 0.0001605381 
    ## ... Similar to previous best
    ## Run 270 stress 9.74451e-05 
    ## ... Procrustes: rmse 0.0002514335  max resid 0.0005200736 
    ## ... Similar to previous best
    ## Run 271 stress 9.144012e-05 
    ## ... Procrustes: rmse 7.986572e-05  max resid 0.0001569225 
    ## ... Similar to previous best
    ## Run 272 stress 9.509803e-05 
    ## ... Procrustes: rmse 9.579474e-05  max resid 0.0001938497 
    ## ... Similar to previous best
    ## Run 273 stress 9.864819e-05 
    ## ... Procrustes: rmse 0.0001628725  max resid 0.0003123349 
    ## ... Similar to previous best
    ## Run 274 stress 9.659782e-05 
    ## ... Procrustes: rmse 0.000251096  max resid 0.0005170105 
    ## ... Similar to previous best
    ## Run 275 stress 9.384439e-05 
    ## ... Procrustes: rmse 0.0002203973  max resid 0.0003496259 
    ## ... Similar to previous best
    ## Run 276 stress 9.89733e-05 
    ## ... Procrustes: rmse 0.000256935  max resid 0.0003698111 
    ## ... Similar to previous best
    ## Run 277 stress 8.698438e-05 
    ## ... Procrustes: rmse 0.0001431634  max resid 0.0002320862 
    ## ... Similar to previous best
    ## Run 278 stress 9.70293e-05 
    ## ... Procrustes: rmse 0.0001748832  max resid 0.000345339 
    ## ... Similar to previous best
    ## Run 279 stress 9.976817e-05 
    ## ... Procrustes: rmse 0.0002296989  max resid 0.0003551215 
    ## ... Similar to previous best
    ## Run 280 stress 9.685308e-05 
    ## ... Procrustes: rmse 0.0001423445  max resid 0.0002279541 
    ## ... Similar to previous best
    ## Run 281 stress 8.51848e-05 
    ## ... Procrustes: rmse 0.0001204478  max resid 0.0001821985 
    ## ... Similar to previous best
    ## Run 282 stress 8.850314e-05 
    ## ... Procrustes: rmse 0.0001708178  max resid 0.00025145 
    ## ... Similar to previous best
    ## Run 283 stress 9.520364e-05 
    ## ... Procrustes: rmse 0.0002499369  max resid 0.0005131708 
    ## ... Similar to previous best
    ## Run 284 stress 9.845535e-05 
    ## ... Procrustes: rmse 0.0002474742  max resid 0.0003622794 
    ## ... Similar to previous best
    ## Run 285 stress 9.608072e-05 
    ## ... Procrustes: rmse 0.0002481853  max resid 0.0005182677 
    ## ... Similar to previous best
    ## Run 286 stress 9.428783e-05 
    ## ... Procrustes: rmse 0.0001725224  max resid 0.0003359192 
    ## ... Similar to previous best
    ## Run 287 stress 9.989319e-05 
    ## ... Procrustes: rmse 0.0002595224  max resid 0.0005344469 
    ## ... Similar to previous best
    ## Run 288 stress 9.18614e-05 
    ## ... Procrustes: rmse 0.000187155  max resid 0.0002949538 
    ## ... Similar to previous best
    ## Run 289 stress 9.381683e-05 
    ## ... Procrustes: rmse 0.0001886991  max resid 0.0002954725 
    ## ... Similar to previous best
    ## Run 290 stress 9.334511e-05 
    ## ... Procrustes: rmse 0.0001682514  max resid 0.0003359887 
    ## ... Similar to previous best
    ## Run 291 stress 8.75291e-05 
    ## ... Procrustes: rmse 9.832618e-05  max resid 0.000155392 
    ## ... Similar to previous best
    ## Run 292 stress 9.702677e-05 
    ## ... Procrustes: rmse 0.0002204181  max resid 0.0003588528 
    ## ... Similar to previous best
    ## Run 293 stress 9.998928e-05 
    ## ... Procrustes: rmse 9.802165e-05  max resid 0.0001875366 
    ## ... Similar to previous best
    ## Run 294 stress 8.120611e-05 
    ## ... Procrustes: rmse 6.686756e-05  max resid 0.0001200961 
    ## ... Similar to previous best
    ## Run 295 stress 9.718799e-05 
    ## ... Procrustes: rmse 0.0001247076  max resid 0.0001799002 
    ## ... Similar to previous best
    ## Run 296 stress 9.916505e-05 
    ## ... Procrustes: rmse 0.0002591631  max resid 0.0005293882 
    ## ... Similar to previous best
    ## Run 297 stress 9.353247e-05 
    ## ... Procrustes: rmse 0.0001708583  max resid 0.0003324004 
    ## ... Similar to previous best
    ## Run 298 stress 9.444272e-05 
    ## ... Procrustes: rmse 0.0002230455  max resid 0.0003573075 
    ## ... Similar to previous best
    ## Run 299 stress 9.850311e-05 
    ## ... Procrustes: rmse 0.000229033  max resid 0.0003640756 
    ## ... Similar to previous best
    ## Run 300 stress 8.618053e-05 
    ## ... Procrustes: rmse 9.973295e-05  max resid 0.0001588949 
    ## ... Similar to previous best
    ## Run 301 stress 0.2854266 
    ## Run 302 stress 9.728111e-05 
    ## ... Procrustes: rmse 0.0002506863  max resid 0.0005192477 
    ## ... Similar to previous best
    ## Run 303 stress 9.602034e-05 
    ## ... Procrustes: rmse 0.000136196  max resid 0.0002244001 
    ## ... Similar to previous best
    ## Run 304 stress 8.372316e-05 
    ## ... Procrustes: rmse 0.0001195957  max resid 0.000207495 
    ## ... Similar to previous best
    ## Run 305 stress 9.090939e-05 
    ## ... Procrustes: rmse 0.0001665212  max resid 0.0003241405 
    ## ... Similar to previous best
    ## Run 306 stress 8.212209e-05 
    ## ... Procrustes: rmse 0.000179696  max resid 0.0003026386 
    ## ... Similar to previous best
    ## Run 307 stress 9.802475e-05 
    ## ... Procrustes: rmse 0.0001842156  max resid 0.0002858277 
    ## ... Similar to previous best
    ## Run 308 stress 9.800174e-05 
    ## ... Procrustes: rmse 0.0001788551  max resid 0.0003557909 
    ## ... Similar to previous best
    ## Run 309 stress 9.298782e-05 
    ## ... Procrustes: rmse 0.00018894  max resid 0.0002904838 
    ## ... Similar to previous best
    ## Run 310 stress 9.401207e-05 
    ## ... Procrustes: rmse 0.0002343812  max resid 0.0003583576 
    ## ... Similar to previous best
    ## Run 311 stress 9.66277e-05 
    ## ... Procrustes: rmse 0.00016164  max resid 0.0002949992 
    ## ... Similar to previous best
    ## Run 312 stress 9.34184e-05 
    ## ... Procrustes: rmse 0.000153523  max resid 0.0002851174 
    ## ... Similar to previous best
    ## Run 313 stress 9.275314e-05 
    ## ... Procrustes: rmse 0.0001976096  max resid 0.0003705392 
    ## ... Similar to previous best
    ## Run 314 stress 0.2848512 
    ## Run 315 stress 9.450927e-05 
    ## ... Procrustes: rmse 0.0002429756  max resid 0.0005078319 
    ## ... Similar to previous best
    ## Run 316 stress 9.894259e-05 
    ## ... Procrustes: rmse 0.0002602648  max resid 0.0003715865 
    ## ... Similar to previous best
    ## Run 317 stress 9.597058e-05 
    ## ... Procrustes: rmse 0.0001388117  max resid 0.0001959253 
    ## ... Similar to previous best
    ## Run 318 stress 7.7878e-05 
    ## ... Procrustes: rmse 8.337184e-05  max resid 0.0001300498 
    ## ... Similar to previous best
    ## Run 319 stress 9.682022e-05 
    ## ... Procrustes: rmse 0.0002535973  max resid 0.0003579202 
    ## ... Similar to previous best
    ## Run 320 stress 9.773071e-05 
    ## ... Procrustes: rmse 0.0001767487  max resid 0.000342883 
    ## ... Similar to previous best
    ## Run 321 stress 0.2775753 
    ## Run 322 stress 9.849305e-05 
    ## ... Procrustes: rmse 0.0002293729  max resid 0.0003566607 
    ## ... Similar to previous best
    ## Run 323 stress 8.6178e-05 
    ## ... Procrustes: rmse 0.0002279625  max resid 0.0003276441 
    ## ... Similar to previous best
    ## Run 324 stress 9.841597e-05 
    ## ... Procrustes: rmse 0.0002310846  max resid 0.0003709996 
    ## ... Similar to previous best
    ## Run 325 stress 9.667093e-05 
    ## ... Procrustes: rmse 0.0001742192  max resid 0.0003477206 
    ## ... Similar to previous best
    ## Run 326 stress 9.795177e-05 
    ## ... Procrustes: rmse 0.0001646645  max resid 0.0003011556 
    ## ... Similar to previous best
    ## Run 327 stress 0.2479299 
    ## Run 328 stress 9.208924e-05 
    ## ... Procrustes: rmse 0.0002384464  max resid 0.0004984018 
    ## ... Similar to previous best
    ## Run 329 stress 9.833983e-05 
    ## ... Procrustes: rmse 0.0002141132  max resid 0.0003530342 
    ## ... Similar to previous best
    ## Run 330 stress 0.2848516 
    ## Run 331 stress 9.676822e-05 
    ## ... Procrustes: rmse 9.329325e-05  max resid 0.0001802516 
    ## ... Similar to previous best
    ## Run 332 stress 9.726394e-05 
    ## ... Procrustes: rmse 0.0002285712  max resid 0.0003558725 
    ## ... Similar to previous best
    ## Run 333 stress 9.78724e-05 
    ## ... Procrustes: rmse 0.00022664  max resid 0.0003520288 
    ## ... Similar to previous best
    ## Run 334 stress 9.222471e-05 
    ## ... Procrustes: rmse 0.0001684227  max resid 0.0003275106 
    ## ... Similar to previous best
    ## Run 335 stress 9.540141e-05 
    ## ... Procrustes: rmse 0.0001312756  max resid 0.0002228667 
    ## ... Similar to previous best
    ## Run 336 stress 9.585871e-05 
    ## ... Procrustes: rmse 0.0002518664  max resid 0.0005198513 
    ## ... Similar to previous best
    ## Run 337 stress 9.430974e-05 
    ## ... Procrustes: rmse 0.0001909587  max resid 0.0002942971 
    ## ... Similar to previous best
    ## Run 338 stress 9.168921e-05 
    ## ... Procrustes: rmse 0.0002130045  max resid 0.0003346838 
    ## ... Similar to previous best
    ## Run 339 stress 9.707292e-05 
    ## ... Procrustes: rmse 0.0002471586  max resid 0.0003680484 
    ## ... Similar to previous best
    ## Run 340 stress 0.2680462 
    ## Run 341 stress 9.485613e-05 
    ## ... Procrustes: rmse 0.0001412137  max resid 0.0002474696 
    ## ... Similar to previous best
    ## Run 342 stress 9.530315e-05 
    ## ... Procrustes: rmse 0.0002451415  max resid 0.0005118119 
    ## ... Similar to previous best
    ## Run 343 stress 9.411545e-05 
    ## ... Procrustes: rmse 0.0002403894  max resid 0.000337987 
    ## ... Similar to previous best
    ## Run 344 stress 9.977254e-05 
    ## ... Procrustes: rmse 0.0002559732  max resid 0.0005232832 
    ## ... Similar to previous best
    ## Run 345 stress 9.449455e-05 
    ## ... Procrustes: rmse 0.0002445027  max resid 0.0003523308 
    ## ... Similar to previous best
    ## Run 346 stress 8.836235e-05 
    ## ... Procrustes: rmse 0.0001312539  max resid 0.0002222017 
    ## ... Similar to previous best
    ## Run 347 stress 9.909493e-05 
    ## ... Procrustes: rmse 0.0002013026  max resid 0.0003038077 
    ## ... Similar to previous best
    ## Run 348 stress 9.783066e-05 
    ## ... Procrustes: rmse 0.0001574515  max resid 0.0002813784 
    ## ... Similar to previous best
    ## Run 349 stress 9.293807e-05 
    ## ... Procrustes: rmse 9.339113e-05  max resid 0.000166649 
    ## ... Similar to previous best
    ## Run 350 stress 9.686292e-05 
    ## ... Procrustes: rmse 8.918339e-05  max resid 0.0001700716 
    ## ... Similar to previous best
    ## Run 351 stress 9.718123e-05 
    ## ... Procrustes: rmse 0.0001793806  max resid 0.0003353028 
    ## ... Similar to previous best
    ## Run 352 stress 9.128294e-05 
    ## ... Procrustes: rmse 0.0001661802  max resid 0.0003226749 
    ## ... Similar to previous best
    ## Run 353 stress 9.876992e-05 
    ## ... Procrustes: rmse 0.000110132  max resid 0.0001591785 
    ## ... Similar to previous best
    ## Run 354 stress 8.283643e-05 
    ## ... Procrustes: rmse 0.0001886658  max resid 0.0003146313 
    ## ... Similar to previous best
    ## Run 355 stress 9.534526e-05 
    ## ... Procrustes: rmse 0.0002189134  max resid 0.0003410682 
    ## ... Similar to previous best
    ## Run 356 stress 9.099298e-05 
    ## ... Procrustes: rmse 0.0001779107  max resid 0.0002625219 
    ## ... Similar to previous best
    ## Run 357 stress 9.756862e-05 
    ## ... Procrustes: rmse 0.0002228885  max resid 0.0003683625 
    ## ... Similar to previous best
    ## Run 358 stress 8.916735e-05 
    ## ... Procrustes: rmse 0.0002325841  max resid 0.0003296449 
    ## ... Similar to previous best
    ## Run 359 stress 8.845673e-05 
    ## ... Procrustes: rmse 0.000176495  max resid 0.0002734822 
    ## ... Similar to previous best
    ## Run 360 stress 9.593234e-05 
    ## ... Procrustes: rmse 0.0002434944  max resid 0.0005089959 
    ## ... Similar to previous best
    ## Run 361 stress 8.894996e-05 
    ## ... Procrustes: rmse 0.0002323787  max resid 0.0003357317 
    ## ... Similar to previous best
    ## Run 362 stress 9.542107e-05 
    ## ... Procrustes: rmse 0.0001357762  max resid 0.0001888002 
    ## ... Similar to previous best
    ## Run 363 stress 9.457578e-05 
    ## ... Procrustes: rmse 0.0002485327  max resid 0.0005158619 
    ## ... Similar to previous best
    ## Run 364 stress 9.44046e-05 
    ## ... Procrustes: rmse 0.0002471832  max resid 0.0005080681 
    ## ... Similar to previous best
    ## Run 365 stress 7.606089e-05 
    ## ... Procrustes: rmse 0.0001361208  max resid 0.0002177636 
    ## ... Similar to previous best
    ## Run 366 stress 7.61585e-05 
    ## ... Procrustes: rmse 0.0001405916  max resid 0.0002836031 
    ## ... Similar to previous best
    ## Run 367 stress 7.35045e-05 
    ## ... Procrustes: rmse 6.768276e-05  max resid 0.0001508407 
    ## ... Similar to previous best
    ## Run 368 stress 0.3083098 
    ## Run 369 stress 8.903064e-05 
    ## ... Procrustes: rmse 0.0002324956  max resid 0.0003322952 
    ## ... Similar to previous best
    ## Run 370 stress 9.458623e-05 
    ## ... Procrustes: rmse 0.0001727833  max resid 0.0003362199 
    ## ... Similar to previous best
    ## Run 371 stress 0.3083098 
    ## Run 372 stress 8.602668e-05 
    ## ... Procrustes: rmse 0.0001095776  max resid 0.0001665774 
    ## ... Similar to previous best
    ## Run 373 stress 8.754219e-05 
    ## ... Procrustes: rmse 8.165771e-05  max resid 0.0001770494 
    ## ... Similar to previous best
    ## Run 374 stress 8.848729e-05 
    ## ... Procrustes: rmse 0.0001027935  max resid 0.0001660578 
    ## ... Similar to previous best
    ## Run 375 stress 9.142653e-05 
    ## ... Procrustes: rmse 0.0002400791  max resid 0.0003464917 
    ## ... Similar to previous best
    ## Run 376 stress 9.036493e-05 
    ## ... Procrustes: rmse 0.0002387535  max resid 0.0003416129 
    ## ... Similar to previous best
    ## Run 377 stress 9.736914e-05 
    ## ... Procrustes: rmse 0.0002568525  max resid 0.0005284572 
    ## ... Similar to previous best
    ## Run 378 stress 8.843284e-05 
    ## ... Procrustes: rmse 8.498005e-05  max resid 0.0001719016 
    ## ... Similar to previous best
    ## Run 379 stress 8.765428e-05 
    ## ... Procrustes: rmse 0.0002257859  max resid 0.0003189037 
    ## ... Similar to previous best
    ## Run 380 stress 9.861683e-05 
    ## ... Procrustes: rmse 0.000230976  max resid 0.0003717321 
    ## ... Similar to previous best
    ## Run 381 stress 9.486131e-05 
    ## ... Procrustes: rmse 8.08341e-05  max resid 0.000138809 
    ## ... Similar to previous best
    ## Run 382 stress 8.492312e-05 
    ## ... Procrustes: rmse 7.966566e-05  max resid 0.0001658546 
    ## ... Similar to previous best
    ## Run 383 stress 9.818142e-05 
    ## ... Procrustes: rmse 0.0002321664  max resid 0.0003691625 
    ## ... Similar to previous best
    ## Run 384 stress 9.827407e-05 
    ## ... Procrustes: rmse 0.0002286866  max resid 0.0003549785 
    ## ... Similar to previous best
    ## Run 385 stress 9.288467e-05 
    ## ... Procrustes: rmse 0.0001716633  max resid 0.0003365277 
    ## ... Similar to previous best
    ## Run 386 stress 9.595073e-05 
    ## ... Procrustes: rmse 0.0001611256  max resid 0.0003058034 
    ## ... Similar to previous best
    ## Run 387 stress 9.745724e-05 
    ## ... Procrustes: rmse 0.000253918  max resid 0.0003744006 
    ## ... Similar to previous best
    ## Run 388 stress 9.185302e-05 
    ## ... Procrustes: rmse 0.0001830908  max resid 0.0002849156 
    ## ... Similar to previous best
    ## Run 389 stress 9.161301e-05 
    ## ... Procrustes: rmse 0.0001822291  max resid 0.0002777686 
    ## ... Similar to previous best
    ## Run 390 stress 5.341047e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 1.409115e-05  max resid 2.2694e-05 
    ## ... Similar to previous best
    ## Run 391 stress 9.23668e-05 
    ## ... Procrustes: rmse 8.230253e-05  max resid 0.0001554362 
    ## ... Similar to previous best
    ## Run 392 stress 9.993858e-05 
    ## ... Procrustes: rmse 0.000114204  max resid 0.0002024903 
    ## ... Similar to previous best
    ## Run 393 stress 9.778449e-05 
    ## ... Procrustes: rmse 0.0002386122  max resid 0.0003646989 
    ## ... Similar to previous best
    ## Run 394 stress 9.139501e-05 
    ## ... Procrustes: rmse 0.0002323374  max resid 0.0003293118 
    ## ... Similar to previous best
    ## Run 395 stress 9.865787e-05 
    ## ... Procrustes: rmse 0.000181689  max resid 0.0003514844 
    ## ... Similar to previous best
    ## Run 396 stress 9.906043e-05 
    ## ... Procrustes: rmse 0.0002292295  max resid 0.0003627546 
    ## ... Similar to previous best
    ## Run 397 stress 8.913816e-05 
    ## ... Procrustes: rmse 7.253243e-05  max resid 0.0001512487 
    ## ... Similar to previous best
    ## Run 398 stress 9.47312e-05 
    ## ... Procrustes: rmse 0.000238771  max resid 0.0003532445 
    ## ... Similar to previous best
    ## Run 399 stress 9.59622e-05 
    ## ... Procrustes: rmse 0.0002171553  max resid 0.0003461701 
    ## ... Similar to previous best
    ## Run 400 stress 9.151351e-05 
    ## ... Procrustes: rmse 0.0001345528  max resid 0.0002450361 
    ## ... Similar to previous best
    ## Run 401 stress 9.418114e-05 
    ## ... Procrustes: rmse 0.0002112697  max resid 0.0003310656 
    ## ... Similar to previous best
    ## Run 402 stress 9.937335e-05 
    ## ... Procrustes: rmse 0.0002181039  max resid 0.0004565201 
    ## ... Similar to previous best
    ## Run 403 stress 9.281285e-05 
    ## ... Procrustes: rmse 0.0001725283  max resid 0.0003404096 
    ## ... Similar to previous best
    ## Run 404 stress 9.28051e-05 
    ## ... Procrustes: rmse 0.0001693121  max resid 0.0003354754 
    ## ... Similar to previous best
    ## Run 405 stress 9.416637e-05 
    ## ... Procrustes: rmse 0.0002283473  max resid 0.0003415701 
    ## ... Similar to previous best
    ## Run 406 stress 9.207756e-05 
    ## ... Procrustes: rmse 0.0001256348  max resid 0.0002054194 
    ## ... Similar to previous best
    ## Run 407 stress 9.694923e-05 
    ## ... Procrustes: rmse 0.0002052867  max resid 0.0003758292 
    ## ... Similar to previous best
    ## Run 408 stress 9.571462e-05 
    ## ... Procrustes: rmse 0.0002183037  max resid 0.0003449915 
    ## ... Similar to previous best
    ## Run 409 stress 9.414454e-05 
    ## ... Procrustes: rmse 0.000172037  max resid 0.0003350886 
    ## ... Similar to previous best
    ## Run 410 stress 9.556362e-05 
    ## ... Procrustes: rmse 0.0002347701  max resid 0.0003348665 
    ## ... Similar to previous best
    ## Run 411 stress 8.952709e-05 
    ## ... Procrustes: rmse 0.0001618416  max resid 0.0003218675 
    ## ... Similar to previous best
    ## Run 412 stress 8.623039e-05 
    ## ... Procrustes: rmse 0.0001051759  max resid 0.0001796955 
    ## ... Similar to previous best
    ## Run 413 stress 9.971993e-05 
    ## ... Procrustes: rmse 0.0002528415  max resid 0.0003742663 
    ## ... Similar to previous best
    ## Run 414 stress 9.919125e-05 
    ## ... Procrustes: rmse 0.0002536955  max resid 0.0003663704 
    ## ... Similar to previous best
    ## Run 415 stress 9.735778e-05 
    ## ... Procrustes: rmse 0.0001939732  max resid 0.0003729323 
    ## ... Similar to previous best
    ## Run 416 stress 8.088891e-05 
    ## ... Procrustes: rmse 0.0001516183  max resid 0.0002968809 
    ## ... Similar to previous best
    ## Run 417 stress 9.937608e-05 
    ## ... Procrustes: rmse 0.0002254508  max resid 0.0003714774 
    ## ... Similar to previous best
    ## Run 418 stress 9.659363e-05 
    ## ... Procrustes: rmse 0.0002403005  max resid 0.000340721 
    ## ... Similar to previous best
    ## Run 419 stress 9.303726e-05 
    ## ... Procrustes: rmse 0.0001019407  max resid 0.0001749378 
    ## ... Similar to previous best
    ## Run 420 stress 9.698727e-05 
    ## ... Procrustes: rmse 0.0001795979  max resid 0.0003540906 
    ## ... Similar to previous best
    ## Run 421 stress 8.986661e-05 
    ## ... Procrustes: rmse 0.0001053195  max resid 0.0001849532 
    ## ... Similar to previous best
    ## Run 422 stress 9.898305e-05 
    ## ... Procrustes: rmse 0.00025049  max resid 0.0003748375 
    ## ... Similar to previous best
    ## Run 423 stress 9.152779e-05 
    ## ... Procrustes: rmse 0.000233466  max resid 0.0003439967 
    ## ... Similar to previous best
    ## Run 424 stress 8.438919e-05 
    ## ... Procrustes: rmse 0.0002107745  max resid 0.0003011505 
    ## ... Similar to previous best
    ## Run 425 stress 8.523591e-05 
    ## ... Procrustes: rmse 0.0001202154  max resid 0.0001878259 
    ## ... Similar to previous best
    ## Run 426 stress 8.566207e-05 
    ## ... Procrustes: rmse 0.0001239474  max resid 0.0001897213 
    ## ... Similar to previous best
    ## Run 427 stress 9.542501e-05 
    ## ... Procrustes: rmse 0.000212377  max resid 0.0003317591 
    ## ... Similar to previous best
    ## Run 428 stress 9.827051e-05 
    ## ... Procrustes: rmse 0.0001107284  max resid 0.0001946599 
    ## ... Similar to previous best
    ## Run 429 stress 9.540584e-05 
    ## ... Procrustes: rmse 0.00017102  max resid 0.0003289947 
    ## ... Similar to previous best
    ## Run 430 stress 9.607678e-05 
    ## ... Procrustes: rmse 0.0002205595  max resid 0.0003647449 
    ## ... Similar to previous best
    ## Run 431 stress 9.232654e-05 
    ## ... Procrustes: rmse 0.0001768896  max resid 0.0002899345 
    ## ... Similar to previous best
    ## Run 432 stress 7.972633e-05 
    ## ... Procrustes: rmse 9.481654e-05  max resid 0.0001670248 
    ## ... Similar to previous best
    ## Run 433 stress 9.76277e-05 
    ## ... Procrustes: rmse 0.0002469385  max resid 0.0003537021 
    ## ... Similar to previous best
    ## Run 434 stress 8.360195e-05 
    ## ... Procrustes: rmse 0.0001139718  max resid 0.0001933238 
    ## ... Similar to previous best
    ## Run 435 stress 9.741934e-05 
    ## ... Procrustes: rmse 0.0002493874  max resid 0.0003644666 
    ## ... Similar to previous best
    ## Run 436 stress 9.296795e-05 
    ## ... Procrustes: rmse 8.922696e-05  max resid 0.0001758738 
    ## ... Similar to previous best
    ## Run 437 stress 9.842124e-05 
    ## ... Procrustes: rmse 0.0001986602  max resid 0.0002531509 
    ## ... Similar to previous best
    ## Run 438 stress 8.829632e-05 
    ## ... Procrustes: rmse 0.0001193089  max resid 0.0001659208 
    ## ... Similar to previous best
    ## Run 439 stress 9.751153e-05 
    ## ... Procrustes: rmse 0.0001056292  max resid 0.0001819177 
    ## ... Similar to previous best
    ## Run 440 stress 8.138114e-05 
    ## ... Procrustes: rmse 6.398694e-05  max resid 9.30281e-05 
    ## ... Similar to previous best
    ## Run 441 stress 9.592335e-05 
    ## ... Procrustes: rmse 0.0001560266  max resid 0.000287174 
    ## ... Similar to previous best
    ## Run 442 stress 0.231921 
    ## Run 443 stress 9.564849e-05 
    ## ... Procrustes: rmse 0.0002171935  max resid 0.000355637 
    ## ... Similar to previous best
    ## Run 444 stress 9.489371e-05 
    ## ... Procrustes: rmse 0.0002416254  max resid 0.0003547494 
    ## ... Similar to previous best
    ## Run 445 stress 9.660607e-05 
    ## ... Procrustes: rmse 0.0002415603  max resid 0.000342485 
    ## ... Similar to previous best
    ## Run 446 stress 9.419349e-05 
    ## ... Procrustes: rmse 0.0001142758  max resid 0.000205644 
    ## ... Similar to previous best
    ## Run 447 stress 9.647149e-05 
    ## ... Procrustes: rmse 0.0002216672  max resid 0.0003578141 
    ## ... Similar to previous best
    ## Run 448 stress 8.837197e-05 
    ## ... Procrustes: rmse 6.686329e-05  max resid 0.0001187393 
    ## ... Similar to previous best
    ## Run 449 stress 9.466631e-05 
    ## ... Procrustes: rmse 0.0002181887  max resid 0.0003471098 
    ## ... Similar to previous best
    ## Run 450 stress 9.50474e-05 
    ## ... Procrustes: rmse 0.0002123135  max resid 0.0003365773 
    ## ... Similar to previous best
    ## Run 451 stress 9.2141e-05 
    ## ... Procrustes: rmse 0.000167368  max resid 0.0002866793 
    ## ... Similar to previous best
    ## Run 452 stress 9.772293e-05 
    ## ... Procrustes: rmse 0.0002489803  max resid 0.0003702118 
    ## ... Similar to previous best
    ## Run 453 stress 9.373081e-05 
    ## ... Procrustes: rmse 0.000111409  max resid 0.0001986901 
    ## ... Similar to previous best
    ## Run 454 stress 9.901043e-05 
    ## ... Procrustes: rmse 0.0002517413  max resid 0.0003720772 
    ## ... Similar to previous best
    ## Run 455 stress 8.538065e-05 
    ## ... Procrustes: rmse 0.0001494415  max resid 0.0002319977 
    ## ... Similar to previous best
    ## Run 456 stress 9.502905e-05 
    ## ... Procrustes: rmse 0.0002113867  max resid 0.0003501439 
    ## ... Similar to previous best
    ## Run 457 stress 9.707956e-05 
    ## ... Procrustes: rmse 0.0002296864  max resid 0.0003417211 
    ## ... Similar to previous best
    ## Run 458 stress 9.157187e-05 
    ## ... Procrustes: rmse 0.000169187  max resid 0.0003278891 
    ## ... Similar to previous best
    ## Run 459 stress 9.301526e-05 
    ## ... Procrustes: rmse 0.000198747  max resid 0.000366573 
    ## ... Similar to previous best
    ## Run 460 stress 9.034437e-05 
    ## ... Procrustes: rmse 0.0001054751  max resid 0.0001848045 
    ## ... Similar to previous best
    ## Run 461 stress 9.720201e-05 
    ## ... Procrustes: rmse 0.0002472922  max resid 0.0003650443 
    ## ... Similar to previous best
    ## Run 462 stress 9.03826e-05 
    ## ... Procrustes: rmse 0.000149072  max resid 0.0002829027 
    ## ... Similar to previous best
    ## Run 463 stress 8.790241e-05 
    ## ... Procrustes: rmse 5.623526e-05  max resid 8.013196e-05 
    ## ... Similar to previous best
    ## Run 464 stress 9.365446e-05 
    ## ... Procrustes: rmse 0.0002346447  max resid 0.0003496618 
    ## ... Similar to previous best
    ## Run 465 stress 9.506014e-05 
    ## ... Procrustes: rmse 0.0001351351  max resid 0.0002064091 
    ## ... Similar to previous best
    ## Run 466 stress 9.271106e-05 
    ## ... Procrustes: rmse 0.0002262399  max resid 0.0003166364 
    ## ... Similar to previous best
    ## Run 467 stress 9.618863e-05 
    ## ... Procrustes: rmse 0.0002425096  max resid 0.0003602378 
    ## ... Similar to previous best
    ## Run 468 stress 9.685788e-05 
    ## ... Procrustes: rmse 0.0002431641  max resid 0.000345745 
    ## ... Similar to previous best
    ## Run 469 stress 9.890545e-05 
    ## ... Procrustes: rmse 0.0001582952  max resid 0.0002844658 
    ## ... Similar to previous best
    ## Run 470 stress 9.807404e-05 
    ## ... Procrustes: rmse 0.0001110936  max resid 0.0002020352 
    ## ... Similar to previous best
    ## Run 471 stress 7.289631e-05 
    ## ... Procrustes: rmse 5.613961e-05  max resid 0.0001034168 
    ## ... Similar to previous best
    ## Run 472 stress 9.462338e-05 
    ## ... Procrustes: rmse 0.0001118003  max resid 0.0001989703 
    ## ... Similar to previous best
    ## Run 473 stress 8.216085e-05 
    ## ... Procrustes: rmse 8.429341e-05  max resid 0.0001231615 
    ## ... Similar to previous best
    ## Run 474 stress 9.752553e-05 
    ## ... Procrustes: rmse 0.0002451381  max resid 0.0003707837 
    ## ... Similar to previous best
    ## Run 475 stress 9.434602e-05 
    ## ... Procrustes: rmse 0.0002180391  max resid 0.0003559077 
    ## ... Similar to previous best
    ## Run 476 stress 8.85826e-05 
    ## ... Procrustes: rmse 5.793125e-05  max resid 0.0001033197 
    ## ... Similar to previous best
    ## Run 477 stress 9.447344e-05 
    ## ... Procrustes: rmse 0.0001271598  max resid 0.000210509 
    ## ... Similar to previous best
    ## Run 478 stress 9.123586e-05 
    ## ... Procrustes: rmse 0.0001769108  max resid 0.0002869407 
    ## ... Similar to previous best
    ## Run 479 stress 8.541358e-05 
    ## ... Procrustes: rmse 7.366486e-05  max resid 0.0001353926 
    ## ... Similar to previous best
    ## Run 480 stress 9.682291e-05 
    ## ... Procrustes: rmse 0.0002478453  max resid 0.000363189 
    ## ... Similar to previous best
    ## Run 481 stress 9.409299e-05 
    ## ... Procrustes: rmse 0.0002404203  max resid 0.0003492362 
    ## ... Similar to previous best
    ## Run 482 stress 9.737535e-05 
    ## ... Procrustes: rmse 0.0002249651  max resid 0.0003603725 
    ## ... Similar to previous best
    ## Run 483 stress 9.196481e-05 
    ## ... Procrustes: rmse 0.000108397  max resid 0.0001838873 
    ## ... Similar to previous best
    ## Run 484 stress 9.422961e-05 
    ## ... Procrustes: rmse 0.0002388371  max resid 0.0003414127 
    ## ... Similar to previous best
    ## Run 485 stress 9.319044e-05 
    ## ... Procrustes: rmse 0.0001519397  max resid 0.0002906591 
    ## ... Similar to previous best
    ## Run 486 stress 9.922765e-05 
    ## ... Procrustes: rmse 0.0002519776  max resid 0.0003727068 
    ## ... Similar to previous best
    ## Run 487 stress 8.758968e-05 
    ## ... Procrustes: rmse 0.0001213418  max resid 0.0001785129 
    ## ... Similar to previous best
    ## Run 488 stress 9.945462e-05 
    ## ... Procrustes: rmse 0.0001974372  max resid 0.0003182558 
    ## ... Similar to previous best
    ## Run 489 stress 8.258601e-05 
    ## ... Procrustes: rmse 0.0002056753  max resid 0.0002923871 
    ## ... Similar to previous best
    ## Run 490 stress 9.479609e-05 
    ## ... Procrustes: rmse 0.0002367868  max resid 0.0003347817 
    ## ... Similar to previous best
    ## Run 491 stress 9.01606e-05 
    ## ... Procrustes: rmse 6.606942e-05  max resid 0.000113854 
    ## ... Similar to previous best
    ## Run 492 stress 8.736956e-05 
    ## ... Procrustes: rmse 0.0001993663  max resid 0.0003255952 
    ## ... Similar to previous best
    ## Run 493 stress 8.816904e-05 
    ## ... Procrustes: rmse 0.0002189941  max resid 0.0003175195 
    ## ... Similar to previous best
    ## Run 494 stress 9.880423e-05 
    ## ... Procrustes: rmse 0.0001378127  max resid 0.0002059177 
    ## ... Similar to previous best
    ## Run 495 stress 9.538249e-05 
    ## ... Procrustes: rmse 8.800409e-05  max resid 0.0001742664 
    ## ... Similar to previous best
    ## Run 496 stress 9.348835e-05 
    ## ... Procrustes: rmse 0.0001964367  max resid 0.0003636817 
    ## ... Similar to previous best
    ## Run 497 stress 7.337882e-05 
    ## ... Procrustes: rmse 0.0001013678  max resid 0.0001356778 
    ## ... Similar to previous best
    ## Run 498 stress 8.889798e-05 
    ## ... Procrustes: rmse 7.0243e-05  max resid 0.000100437 
    ## ... Similar to previous best
    ## Run 499 stress 8.737808e-05 
    ## ... Procrustes: rmse 0.000118274  max resid 0.000191175 
    ## ... Similar to previous best
    ## Run 500 stress 9.532352e-05 
    ## ... Procrustes: rmse 0.0002126623  max resid 0.0003380585 
    ## ... Similar to previous best
    ## *** Best solution repeated 110 times

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
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31833 0.699  0.019 *
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
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31833 0.699  0.019 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by Site_type
(PD_beta_env_A_ef <- envfit(PD_beta_env_NMDS, env[surveyed_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.75819 -0.65203 0.0724  0.699  
    ## salinity_median     0.94798  0.31833 0.6990  0.023 *
    ## oxygen_median       0.57185 -0.82036 0.6836  0.115  
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
    ## temperature_median -0.75819 -0.65203 0.0724  1.000  
    ## salinity_median     0.94798  0.31833 0.6990  0.069 .
    ## oxygen_median       0.57185 -0.82036 0.6836  0.345  
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
    ## temperature_median -0.92661  0.37601 0.0954  0.747  
    ## salinity_median     0.90470  0.42604 0.6840  0.034 *
    ## oxygen_median       0.44477 -0.89564 0.5837  0.198  
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
    ## temperature_median -0.92661  0.37601 0.0954  1.000
    ## salinity_median     0.90470  0.42604 0.6840  0.102
    ## oxygen_median       0.44477 -0.89564 0.5837  0.594
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
    ## temperature_median -0.52519 -0.85098 0.0223  0.898  
    ## salinity_median    -0.79785 -0.60286 0.7269  0.016 *
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.011 *
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
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.033 *
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
    ## temperature_median  0.0001220  1.0000000 0.1090  0.362
    ## salinity_median    -0.0091253 -0.9999600 0.4823  0.944
    ## oxygen_median      -0.0025077  1.0000000 0.7159  0.706
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
    ## temperature_median  0.0001220  1.0000000 0.1090      1
    ## salinity_median    -0.0091253 -0.9999600 0.4823      1
    ## oxygen_median      -0.0025077  1.0000000 0.7159      1
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
    ## temperature_median -0.00006819 -1.00000000 0.1360  0.653  
    ## salinity_median     0.00019801  1.00000000 0.6386  0.069 .
    ## oxygen_median       0.00042353 -1.00000000 0.6682  0.058 .
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
    ## temperature_median -0.00006819 -1.00000000 0.1360  1.000
    ## salinity_median     0.00019801  1.00000000 0.6386  0.207
    ## oxygen_median       0.00042353 -1.00000000 0.6682  0.174
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
    ## temperature_median       -0.39652 -0.91802 0.0381  0.907  
    ## salinity_median          -0.96164  0.27433 0.5830  0.109  
    ## oxygen_median             0.83740  0.54660 0.7450  0.027 *
    ## distance_to_ocean_mean_m -0.09987  0.99500 0.0319  0.921  
    ## max_depth                -0.87446  0.48510 0.1022  0.771  
    ## logArea                  -0.40925  0.91242 0.2148  0.545  
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
    ## temperature_median       -0.39652 -0.91802 0.0381  1.000
    ## salinity_median          -0.96164  0.27433 0.5830  0.654
    ## oxygen_median             0.83740  0.54660 0.7450  0.162
    ## distance_to_ocean_mean_m -0.09987  0.99500 0.0319  1.000
    ## max_depth                -0.87446  0.48510 0.1022  1.000
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
    ## distance_to_ocean_min_m -0.82789  0.56090 0.5558  0.070 .
    ## max_depth               -0.20870 -0.97798 0.0858  0.877  
    ## logArea                  0.25302 -0.96746 0.2012  0.268  
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
    ## distance_to_ocean_min_m -0.82789  0.56090 0.5558  0.210
    ## max_depth               -0.20870 -0.97798 0.0858  1.000
    ## logArea                  0.25302 -0.96746 0.2012  0.804
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by Site_type
(PD_beta_geo_A_ef <- envfit(PD_beta_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.82789  0.56090 0.5558  0.088 .
    ## max_depth               -0.20870 -0.97798 0.0858  0.894  
    ## logArea                  0.25302 -0.96746 0.2012  0.271  
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
    ## distance_to_ocean_min_m -0.82789  0.56090 0.5558  0.264
    ## max_depth               -0.20870 -0.97798 0.0858  1.000
    ## logArea                  0.25302 -0.96746 0.2012  0.813
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
    ## distance_to_ocean_min_m -0.76303  0.64636 0.4671  0.121
    ## max_depth               -0.34716 -0.93781 0.0777  0.972
    ## logArea                  0.55425  0.83235 0.0074  0.967
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
    ## distance_to_ocean_min_m -0.76303  0.64636 0.4671  0.363
    ## max_depth               -0.34716 -0.93781 0.0777  1.000
    ## logArea                  0.55425  0.83235 0.0074  1.000
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
    ## distance_to_ocean_min_m  0.45615 -0.88990 0.4273  0.165   
    ## max_depth               -0.94119 -0.33789 0.6682  0.002 **
    ## logArea                 -0.86397  0.50355 0.2339  0.445   
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
    ## distance_to_ocean_min_m  0.45615 -0.88990 0.4273  0.495   
    ## max_depth               -0.94119 -0.33789 0.6682  0.006 **
    ## logArea                 -0.86397  0.50355 0.2339  1.000   
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
    ## distance_to_ocean_min_m  0.99475  0.10238 0.6338  0.238
    ## max_depth                0.73039  0.68303 0.0374  0.976
    ## logArea                 -0.39144  0.92020 0.3078  0.296
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
    ## distance_to_ocean_min_m  0.99475  0.10238 0.6338  0.714
    ## max_depth                0.73039  0.68303 0.0374  1.000
    ## logArea                 -0.39144  0.92020 0.3078  0.888
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
    ## distance_to_ocean_min_m -0.00016906 -1.00000000 0.6362  0.081 .
    ## max_depth                0.00121689  1.00000000 0.7172  0.055 .
    ## logArea                  0.00038183 -1.00000000 0.7842  0.017 *
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
    ## distance_to_ocean_min_m -0.00016906 -1.00000000 0.6362  0.243  
    ## max_depth                0.00121689  1.00000000 0.7172  0.165  
    ## logArea                  0.00038183 -1.00000000 0.7842  0.051 .
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
    ##       Significance: 0.205 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.212 0.240 0.260 0.294 
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
    ##       Significance: 0.007 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.436 0.463 0.480 0.507 
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
    ##       Significance: 0.249 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.519 0.555 0.583 0.606 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_mant_pv <- rbind(PD_beta_env_mant_t$signif, PD_beta_env_mant_s$signif, PD_beta_env_mant_o$signif)
PD_beta_env_mant_pv <- PD_beta_env_mant_pv[,1]
(PD_beta_env_mant_pv <- p.adjust(PD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.615 0.021 0.747

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
    ## 0.209 0.247 0.267 0.302 
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
    ##       Significance: 0.01 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.374 0.404 0.429 0.453 
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
    ##       Significance: 0.32 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.415 0.449 0.484 0.538 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_MS_mant_pv <- rbind(PD_beta_env_MS_mant_t$signif, PD_beta_env_MS_mant_s$signif, PD_beta_env_MS_mant_o$signif)
PD_beta_env_MS_mant_pv <- PD_beta_env_MS_mant_pv[,1]
(PD_beta_env_MS_mant_pv <- p.adjust(PD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.699 0.030 0.960

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
    ##       Significance: 0.512 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.134 0.182 0.225 0.334 
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
    ##       Significance: 0.024 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.317 0.368 0.437 0.494 
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
    ## 0.347 0.430 0.502 0.574 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_OM_mant_pv <- rbind(PD_beta_env_OM_mant_t$signif, PD_beta_env_OM_mant_s$signif, PD_beta_env_OM_mant_o$signif)
PD_beta_env_OM_mant_pv <- PD_beta_env_OM_mant_pv[,1]
(PD_beta_env_OM_mant_pv <- p.adjust(PD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.072 0.042

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
    ##       Significance: 0.454 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0367 0.0652 0.0950 0.1289 
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
    ##       Significance: 0.024 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.395 0.435 0.469 0.501 
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
    ##       Significance: 0.237 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.693 0.718 0.738 0.771 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_SO_mant_pv <- rbind(PD_beta_env_SO_mant_t$signif, PD_beta_env_SO_mant_s$signif, PD_beta_env_SO_mant_o$signif)
PD_beta_env_SO_mant_pv <- PD_beta_env_SO_mant_pv[,1]
(PD_beta_env_SO_mant_pv <- p.adjust(PD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.072 0.711

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
    ##       Significance: 0.782 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.235 0.359 0.465 0.700 
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
    ##       Significance: 0.229 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.202 0.261 0.306 0.371 
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
    ##       Significance: 0.078 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.238 0.347 0.423 0.587 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_M_mant_pv <- rbind(PD_beta_env_M_mant_t$signif, PD_beta_env_M_mant_s$signif, PD_beta_env_M_mant_o$signif)
PD_beta_env_M_mant_pv <- PD_beta_env_M_mant_pv[,1]
(PD_beta_env_M_mant_pv <- p.adjust(PD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.687 0.234

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
    ##       Significance: 0.279 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.344 0.452 0.535 0.591 
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
    ##       Significance: 0.622 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.457 0.485 0.507 0.524 
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
    ##       Significance: 0.79 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.164 0.190 0.214 0.238 
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
    ##       Significance: 0.721 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.112 0.132 0.146 0.163 
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
    ##       Significance: 0.654 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.310 0.342 0.371 0.407 
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
    ##       Significance: 0.915 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.156 0.192 0.212 0.254 
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
    ## 0.128 0.153 0.172 0.199 
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
    ##       Significance: 0.422 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.178 0.208 0.238 0.272 
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
    ##       Significance: 0.018 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.162 0.220 0.295 0.349 
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
    ##       Significance: 0.07 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.178 0.213 0.249 0.295 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_OM_mant_pv <- rbind( PD_beta_geo_OM_mant_dmean$signif, PD_beta_geo_OM_mant_md$signif, PD_beta_geo_OM_mant_la$signif)
PD_beta_geo_OM_mant_pv <- PD_beta_geo_OM_mant_pv[,1]
(PD_beta_geo_OM_mant_pv <- p.adjust(PD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.054 0.210

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
    ##       Significance: 0.548 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.621 0.649 0.667 0.697 
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
    ##       Significance: 0.794 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.115 0.143 0.163 0.205 
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
    ##       Significance: 0.212 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.376 0.398 0.417 0.430 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_SO_mant_pv <- rbind( PD_beta_geo_SO_mant_dmean$signif, PD_beta_geo_SO_mant_md$signif, PD_beta_geo_SO_mant_la$signif)
PD_beta_geo_SO_mant_pv <- PD_beta_geo_SO_mant_pv[,1]
(PD_beta_geo_SO_mant_pv <- p.adjust(PD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.636

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
    ##       Significance: 0.719 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.225 0.311 0.420 0.683 
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
    ##       Significance: 0.009 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.261 0.409 0.519 0.568 
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
    ##       Significance: 0.016 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.201 0.281 0.335 0.420 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_M_mant_pv <- rbind( PD_beta_geo_M_mant_dmean$signif, PD_beta_geo_M_mant_md$signif, PD_beta_geo_M_mant_la$signif)
PD_beta_geo_M_mant_pv <- PD_beta_geo_M_mant_pv[,1]
(PD_beta_geo_M_mant_pv <- p.adjust(PD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.027 0.048

### PD beta NMDS ordination plots

``` r
# PD beta total NMDS scores
PD_beta_NMDS_data.scores <- as.data.frame(scores(PD_beta_NMDS))
PD_beta_NMDS_data.scores$Site_type <- env[surveyed_sites,19]
PD_beta_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
PD_beta_NMDS_data.scores$Site_type <- factor(PD_beta_NMDS_data.scores$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

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
PD_beta_ef_plot <- ggplot(data = PD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = PD_beta_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
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
PD_beta_ef_plot_CI90 <- ggplot(data = PD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.90) +
  geom_point(data = PD_beta_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
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
PD_beta_rep_NMDS_data.scores$Site_type <- env[surveyed_sites,19]
PD_beta_rep_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
PD_beta_rep_NMDS_data.scores$Site_type <- factor(PD_beta_rep_NMDS_data.scores$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of PD beta replacement with CI = 0.95
PD_beta_rep_plot <- ggplot(data = PD_beta_rep_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = PD_beta_rep_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
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
PD_beta_ric_NMDS_data.scores$Site_type <- env[surveyed_sites,19]
PD_beta_ric_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
PD_beta_ric_NMDS_data.scores$Site_type <- factor(PD_beta_ric_NMDS_data.scores$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of PD beta richness with CI = 0.95
PD_beta_ric_plot <- ggplot(data = PD_beta_ric_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = PD_beta_ric_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
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
PD_beta_ref_NMDS_data.scores$Site_type <- env[,19]
PD_beta_ref_NMDS_data.scores$Lakes <- env[,1]
PD_beta_ref_NMDS_data.scores$Site_type <- factor(PD_beta_ref_NMDS_data.scores$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of PD beta ref with CI = 0.95
PD_beta_ref_plot <- ggplot(data = PD_beta_ref_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = PD_beta_ref_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
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
  group_by(Site_type) %>%
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
(outlier_PD_beta_NMDS1_plot <- ggplot(outlier_PD_beta_NMDS1, aes(x = Site_type, y = NMDS1, fill = Site_type)) +
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

![](PD_analyses_files/figure-gfm/PD%20beta%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD/outlier_PD_beta_NMDS1.jpg", outlier_PD_beta_NMDS1_plot, width = 6.26, height = 6, units = "in")


# Determine NMDS2 outliers
ordered_PD_beta_NMDS2_data.scores <- PD_beta_NMDS_data.scores[order(PD_beta_NMDS_data.scores$NMDS2), ]

outlier_PD_beta_NMDS2 <- ordered_PD_beta_NMDS2_data.scores %>%
  group_by(Site_type) %>%
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
(outlier_PD_beta_NMDS2_plot <- ggplot(outlier_PD_beta_NMDS2, aes(x = Site_type, y = NMDS2, fill = Site_type)) +
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
PD_beta_mean_dist <- aggregate(PD_beta_dist_matrix, by = list(Site_type_group), FUN = mean, na.rm = TRUE)

# Rename rows, delete column, and transpose matrix and make it a data frame
row.names(PD_beta_mean_dist) <- PD_beta_mean_dist$Group.1
PD_beta_mean_dist <- PD_beta_mean_dist[,-1] 
PD_beta_mean_dist <- as.data.frame(t(PD_beta_mean_dist))

# Add Site_type column
PD_beta_mean_dist$Site_type <- env[surveyed_sites,19]
PD_beta_mean_dist$Site_type <- factor(PD_beta_mean_dist$Site_type, levels = c("Ocean", "Mixed", "Stratified"))


# Determine ocean sites outliers based on distance from other ocean sites
outlier_PD_beta_mean_dist_ocean <- PD_beta_mean_dist[ocean_sites,] %>%
  group_by(Site_type) %>%
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
  group_by(Site_type) %>%
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
  group_by(Site_type) %>%
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

(outlier_PD_beta_mean_dist_plot <- ggplot(outlier_PD_beta_mean_dist, aes(x = Site_type, y = Distances, fill = Site_type)) +
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
stree_sesmpd_plot <- ggplot(data = stree_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Site_type, fill = Site_type)) + 
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
  labs(y= "p-value", x= "z-score", colour = "Site type:", fill = "Site type:", tag = "D") +
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
stree_sesmntd_plot <- ggplot(data = stree_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Site_type, fill = Site_type)) + 
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
  labs(y= "p-value", x= "z-score", colour = "Site type:", fill = "Site type:", tag = "E") +
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
stree_sespd_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Site_type, fill = Site_type)) + 
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
  labs(y= "p-value", x= "z-score", colour = "Site type:", fill = "Site type:", tag = "F") +
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
stree_sesmpd_plot <- ggplot(data = stree_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Site_type, fill = Site_type)) + 
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
stree_sesmntd_plot <- ggplot(data = stree_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Site_type, fill = Site_type)) + 
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
stree_sespd_plot <- ggplot(data = stree_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Site_type, fill = Site_type)) + 
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
  labs(y= "Observed PD", x= "z-score", colour = "Site type:", fill = "Site type:") +
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
