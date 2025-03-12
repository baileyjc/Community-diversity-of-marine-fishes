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
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                     Sum Sq Df F value  Pr(>F)  
    ## (Intercept)        3454556  1  8.4830 0.01136 *
    ## salinity_median      60080  1  0.1475 0.70668  
    ## oxygen_median       325406  1  0.7991 0.38648  
    ## temperature_median  902998  1  2.2174 0.15864  
    ## pH_median          3432036  1  8.4277 0.01157 *
    ## Residuals          5701242 14                  
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
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                     Sum Sq Df F value   Pr(>F)   
    ## (Intercept)         430749  1  0.7074 0.413504   
    ## salinity_median    7052231  1 11.5822 0.003931 **
    ## oxygen_median      5692040  1  9.3483 0.007980 **
    ## temperature_median   40929  1  0.0672 0.798954   
    ## Residuals          9133278 15                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

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
environment <- c("salinity_median", "oxygen_median", "temperature_median")

mod <- aov(salinity_median ~ Site_type, stree_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2 147.34   73.67   13.61 0.000353 ***
    ## Residuals   16  86.62    5.41                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 3 observations deleted due to missingness

``` r
mod <- aov(oxygen_median ~ Site_type, stree_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2 14.440    7.22      19 5.95e-05 ***
    ## Residuals   16  6.081    0.38                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 3 observations deleted due to missingness

``` r
mod <- aov(temperature_median ~ Site_type, stree_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2   2.94   1.470   1.092  0.359
    ## Residuals   16  21.54   1.346               
    ## 3 observations deleted due to missingness

``` r
SR_ss_env <- stree_sespd_env[surveyed_sites_env,]
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
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                             Sum Sq Df F value  Pr(>F)  
    ## (Intercept)                 334277  1  0.3577 0.56456  
    ## volume_m3_w_chemocline      679038  1  0.7265 0.41612  
    ## volume_m3                   454018  1  0.4858 0.50342  
    ## surface_area_m2            3299154  1  3.5299 0.09298 .
    ## distance_to_ocean_min_m      99395  1  0.1063 0.75180  
    ## distance_to_ocean_mean_m    245998  1  0.2632 0.62028  
    ## distance_to_ocean_median_m  316251  1  0.3384 0.57505  
    ## tidal_lag_time_minutes      342739  1  0.3667 0.55976  
    ## tidal_efficiency             19133  1  0.0205 0.88938  
    ## perimeter_fromSat           508096  1  0.5436 0.47971  
    ## max_depth                   385984  1  0.4130 0.53649  
    ## logArea                    2162579  1  2.3138 0.16255  
    ## Residuals                  8411739  9                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

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
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                           Sum Sq Df F value  Pr(>F)  
    ## (Intercept)              1851004  1  1.5508 0.22898  
    ## distance_to_ocean_min_m  9591272  1  8.0359 0.01099 *
    ## max_depth                  82010  1  0.0687 0.79620  
    ## logArea                   899626  1  0.7537 0.39672  
    ## Residuals               21484032 18                  
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

mod <- aov(distance_to_ocean_min_m ~ Site_type, stree_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2  80935   40468   16.49 7.05e-05 ***
    ## Residuals   19  46637    2455                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mod <- aov(max_depth ~ Site_type, stree_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2    535   267.5   2.235  0.134
    ## Residuals   19   2274   119.7

``` r
mod <- aov(logArea ~ Site_type, stree_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2  14.37   7.184   2.341  0.123
    ## Residuals   19  58.32   3.069

``` r
SR_ss_geo <- stree_sespd_env[surveyed_sites,]
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
mod <-  lm(pd.obs ~ salinity_median + oxygen_median + temperature_median + distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ salinity_median + oxygen_median + temperature_median + 
    ##     distance_to_ocean_min_m + max_depth + logArea, data = stree_sespd_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1338.44  -234.51   -29.79   376.69  1083.33 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -6532.6687  6930.2272  -0.943   0.3645  
    ## salinity_median           221.6224    90.9978   2.435   0.0314 *
    ## oxygen_median             422.0351   194.7223   2.167   0.0510 .
    ## temperature_median       -104.4125   186.7860  -0.559   0.5864  
    ## distance_to_ocean_min_m     0.9776     4.2929   0.228   0.8237  
    ## max_depth                 -39.9509    25.2787  -1.580   0.1400  
    ## logArea                   409.5474   173.0952   2.366   0.0357 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 717.9 on 12 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.8183, Adjusted R-squared:  0.7275 
    ## F-statistic:  9.01 on 6 and 12 DF,  p-value: 0.0007177

``` r
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                          Sum Sq Df F value  Pr(>F)  
    ## (Intercept)              457976  1  0.8886 0.36446  
    ## salinity_median         3057189  1  5.9315 0.03142 *
    ## oxygen_median           2421151  1  4.6975 0.05103 .
    ## temperature_median       161054  1  0.3125 0.58645  
    ## distance_to_ocean_min_m   26727  1  0.0519 0.82370  
    ## max_depth               1287359  1  2.4977 0.14000  
    ## logArea                 2885324  1  5.5981 0.03566 *
    ## Residuals               6184967 12                  
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
mod <-  lm(pd.obs ~ oxygen_median + temperature_median + distance_to_ocean_min_m + logArea, stree_sespd_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ oxygen_median + temperature_median + distance_to_ocean_min_m + 
    ##     logArea, data = stree_sespd_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1089.1  -515.5  -131.7   467.2  1618.6 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             6866.106   6017.217   1.141   0.2730  
    ## oxygen_median            529.669    224.507   2.359   0.0334 *
    ## temperature_median      -276.362    207.914  -1.329   0.2050  
    ## distance_to_ocean_min_m   -7.323      3.103  -2.360   0.0333 *
    ## logArea                  240.988    144.789   1.664   0.1182  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 856.9 on 14 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:  0.6981, Adjusted R-squared:  0.6118 
    ## F-statistic: 8.092 on 4 and 14 DF,  p-value: 0.001346

``` r
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                           Sum Sq Df F value  Pr(>F)  
    ## (Intercept)               956066  1  1.3021 0.27298  
    ## oxygen_median            4087034  1  5.5661 0.03337 *
    ## temperature_median       1297314  1  1.7668 0.20503  
    ## distance_to_ocean_min_m  4089851  1  5.5699 0.03332 *
    ## logArea                  2034131  1  2.7703 0.11824  
    ## Residuals               10279839 14                  
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
par(mfrow=c(2,2)) 
#### Environmental
### logSRic ANOVA
## Surveyed sites
# Temperature
PD_z_lm_temp <- lm(pd.obs.z ~  Site_type * temperature_median, data = stree_sespd_env[surveyed_sites_env,])
plot(PD_z_lm_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-1.png)<!-- -->

``` r
PD_z_am_temp <- aov(pd.obs.z ~  Site_type * temperature_median, data = stree_sespd_env[surveyed_sites_env,])
# Salinity
PD_z_lm_sal <- lm(pd.obs.z ~  Site_type * salinity_median, data = stree_sespd_env[surveyed_sites_env,])
plot(PD_z_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-2.png)<!-- -->

``` r
PD_z_am_sal <- aov(pd.obs.z ~  Site_type * salinity_median, data = stree_sespd_env[surveyed_sites_env,])
# Oxygen
PD_z_lm_oxy <- lm(pd.obs.z ~  Site_type * oxygen_median, data = stree_sespd_env[surveyed_sites_env,])
plot(PD_z_lm_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-3.png)<!-- -->

``` r
PD_z_am_oxy <- aov(pd.obs.z ~  Site_type * oxygen_median, data = stree_sespd_env[surveyed_sites_env,])
# Anova outputs
summary(PD_z_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                     2  5.732  2.8659   3.806  0.050 *
    ## temperature_median            1  0.507  0.5066   0.673  0.427  
    ## Site_type:temperature_median  2  2.624  1.3121   1.743  0.214  
    ## Residuals                    13  9.788  0.7529                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_z_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                  2  5.732  2.8659   3.758 0.0515 .
    ## salinity_median            1  0.751  0.7515   0.985 0.3390  
    ## Site_type:salinity_median  2  2.253  1.1263   1.477 0.2643  
    ## Residuals                 13  9.914  0.7626                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_z_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                2  5.732  2.8659   4.736 0.0285 *
    ## oxygen_median            1  0.014  0.0137   0.023 0.8826  
    ## Site_type:oxygen_median  2  5.038  2.5189   4.163 0.0401 *
    ## Residuals               13  7.867  0.6051                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_z_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_z_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_z_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.4497349 1.0000000 1.0000000 0.4637951 1.0000000 1.0000000 0.2565862
    ## [8] 1.0000000 0.3606397

``` r
## Mixed and stratified lakes
# Temperature
PD_z_MS_lm_temp <- lm(pd.obs.z ~  Site_type * temperature_median, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_z_MS_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-4.png)<!-- -->

``` r
PD_z_MS_am_temp <- aov(pd.obs.z ~  Site_type * temperature_median, data = stree_sespd_env[mixed_stratified_lakes,])
# Salinity
PD_z_MS_lm_sal <- lm(pd.obs.z ~  Site_type * salinity_median, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_z_MS_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-5.png)<!-- -->

``` r
PD_z_MS_am_sal <- aov(pd.obs.z ~  Site_type * salinity_median, data = stree_sespd_env[mixed_stratified_lakes,])
# Oxygen
PD_z_MS_lm_oxy <- lm(pd.obs.z ~  Site_type * oxygen_median, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_z_MS_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-6.png)<!-- -->

``` r
PD_z_MS_am_oxy <- aov(pd.obs.z ~  Site_type * oxygen_median, data = stree_sespd_env[mixed_stratified_lakes,])
# Anova outputs
summary(PD_z_MS_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type                     1  0.000  0.0000   0.000  0.999
    ## temperature_median            1  1.166  1.1663   1.694  0.218
    ## Site_type:temperature_median  1  0.373  0.3734   0.542  0.476
    ## Residuals                    12  8.264  0.6887

``` r
summary(PD_z_MS_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type                  1  0.000  0.0000   0.000  0.999
    ## salinity_median            1  0.775  0.7746   1.287  0.279
    ## Site_type:salinity_median  1  1.807  1.8070   3.002  0.109
    ## Residuals                 12  7.222  0.6019

``` r
summary(PD_z_MS_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                1  0.000  0.0000   0.000 0.9994  
    ## oxygen_median            1  0.039  0.0393   0.065 0.8032  
    ## Site_type:oxygen_median  1  2.500  2.4996   4.129 0.0649 .
    ## Residuals               12  7.265  0.6054                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_z_MS_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_z_MS_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_z_MS_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 0.9786888 1.0000000
    ## [8] 1.0000000 0.5841744

``` r
## Ocean sites and mixed lakes
# Temperature
PD_z_OM_lm_temp <- lm(pd.obs.z ~  Site_type * temperature_median, data = stree_sespd_env[ocean_mixed_sites_env,])
plot(PD_z_OM_lm_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-7.png)<!-- -->

``` r
PD_z_OM_am_temp <- aov(pd.obs.z ~  Site_type * temperature_median, data = stree_sespd_env[ocean_mixed_sites_env,])
# Salinity
PD_z_OM_lm_sal <- lm(pd.obs.z ~  Site_type * salinity_median, data = stree_sespd_env[ocean_mixed_sites_env,])
plot(PD_z_OM_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-8.png)<!-- -->

``` r
PD_z_OM_am_sal <- aov(pd.obs.z ~  Site_type * salinity_median, data = stree_sespd_env[ocean_mixed_sites_env,])
# Oxygen
PD_z_OM_lm_oxy <- lm(pd.obs.z ~  Site_type * oxygen_median, data = stree_sespd_env[ocean_mixed_sites_env,])
plot(PD_z_OM_lm_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-9.png)<!-- -->

``` r
PD_z_OM_am_oxy <- aov(pd.obs.z ~  Site_type * oxygen_median, data = stree_sespd_env[ocean_mixed_sites_env,])
# Anova outputs
summary(PD_z_OM_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                     1  4.949   4.949   6.278 0.0407 *
    ## temperature_median            1  0.000   0.000   0.000 0.9886  
    ## Site_type:temperature_median  1  2.515   2.515   3.190 0.1173  
    ## Residuals                     7  5.518   0.788                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_z_OM_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                  1  4.949   4.949   6.256 0.0409 *
    ## salinity_median            1  1.722   1.722   2.176 0.1837  
    ## Site_type:salinity_median  1  0.774   0.774   0.978 0.3556  
    ## Residuals                  7  5.538   0.791                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_z_OM_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                1  4.949   4.949   8.745 0.0212 *
    ## oxygen_median            1  3.012   3.012   5.322 0.0544 .
    ## Site_type:oxygen_median  1  1.059   1.059   1.871 0.2136  
    ## Residuals                7  3.962   0.566                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_z_OM_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_z_OM_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_z_OM_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.3659025 1.0000000 1.0000000 0.3682672 1.0000000 1.0000000 0.1907168
    ## [8] 0.4898398 1.0000000

``` r
## Stratified lakes and ocean sites
# Temperature
PD_z_SO_lm_temp <- lm(pd.obs.z ~  Site_type * temperature_median, data = stree_sespd_env[ocean_stratified_sites_env,])
plot(PD_z_SO_lm_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-10.png)<!-- -->

``` r
PD_z_SO_am_temp <- aov(pd.obs.z ~  Site_type * temperature_median, data = stree_sespd_env[ocean_stratified_sites_env,])
# Salinity
PD_z_SO_lm_sal <- lm(pd.obs.z ~  Site_type * salinity_median, data = stree_sespd_env[ocean_stratified_sites_env,])
plot(PD_z_SO_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-11.png)<!-- -->

``` r
PD_z_SO_am_sal <- aov(pd.obs.z ~  Site_type * salinity_median, data = stree_sespd_env[ocean_stratified_sites_env,])
# Oxygen
PD_z_SO_lm_oxy <- lm(pd.obs.z ~  Site_type * oxygen_median, data = stree_sespd_env[ocean_stratified_sites_env,])
plot(PD_z_SO_lm_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-12.png)<!-- -->

``` r
PD_z_SO_am_oxy <- aov(pd.obs.z ~  Site_type * oxygen_median, data = stree_sespd_env[ocean_stratified_sites_env,])
# Anova outputs
summary(PD_z_SO_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                     1  4.951   4.951   5.983 0.0444 *
    ## temperature_median            1  0.161   0.161   0.194 0.6727  
    ## Site_type:temperature_median  1  2.047   2.047   2.473 0.1598  
    ## Residuals                     7  5.793   0.828                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_z_SO_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                  1  4.951   4.951   4.903 0.0624 .
    ## salinity_median            1  0.490   0.490   0.485 0.5085  
    ## Site_type:salinity_median  1  0.441   0.441   0.437 0.5297  
    ## Residuals                  7  7.069   1.010                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_z_SO_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type                1  4.951   4.951   7.690 0.0276 *
    ## oxygen_median            1  0.362   0.362   0.563 0.4775  
    ## Site_type:oxygen_median  1  3.131   3.131   4.863 0.0633 .
    ## Residuals                7  4.507   0.644                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_z_SO_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_z_SO_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_z_SO_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.3992962 1.0000000 1.0000000 0.5616393 1.0000000 1.0000000 0.2481493
    ## [8] 1.0000000 0.5692635

``` r
## Ocean sites
PD_z_env_O_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_sites_env,])
summary(PD_z_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[ocean_sites_env, ])
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
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
p_values <- summary(PD_z_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
## Mixed lakes
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
Anova(PD_z_env_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                     Sum Sq Df F value Pr(>F)
    ## (Intercept)        0.16000  1  0.2962 0.6152
    ## salinity_median    0.52010  1  0.9630 0.3820
    ## oxygen_median      0.52992  1  0.9811 0.3780
    ## temperature_median 0.20518  1  0.3799 0.5710
    ## Residuals          2.16043  4

``` r
p_values <- summary(PD_z_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
## Stratified lakes
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
Anova(PD_z_env_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                    Sum Sq Df F value Pr(>F)
    ## (Intercept)        0.7176  1  0.8580 0.4067
    ## salinity_median    0.0460  1  0.0550 0.8261
    ## oxygen_median      0.5853  1  0.6998 0.4499
    ## temperature_median 0.5547  1  0.6633 0.4611
    ## Residuals          3.3454  4

``` r
p_values <- summary(PD_z_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
### PRic ANOVA
## Surveyed sites
# Temperature
PD_PRic_lm_temp <- lm(pd.obs ~ Site_type * temperature_median, data = stree_sespd_env[surveyed_sites_env,])
plot(PD_PRic_lm_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-13.png)<!-- -->

``` r
PD_PRic_temp_am <- aov(pd.obs ~ Site_type * temperature_median, data = stree_sespd_env[surveyed_sites_env,])
# Salinity
PD_PRic_lm_sal <- lm(pd.obs ~ Site_type * salinity_median, data = stree_sespd_env[surveyed_sites_env,])
plot(PD_PRic_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-14.png)<!-- -->

``` r
PD_PRic_sal_am <- aov(pd.obs ~ Site_type * salinity_median, data = stree_sespd_env[surveyed_sites_env,])
# Oxygen
PD_PRic_lm_oxy <- lm(pd.obs ~ Site_type * oxygen_median, data = stree_sespd_env[surveyed_sites_env,])
plot(PD_PRic_lm_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-15.png)<!-- -->

``` r
PD_PRic_oxy_am <- aov(pd.obs ~ Site_type * oxygen_median, data = stree_sespd_env[surveyed_sites_env,])
# Anova outputs
summary(PD_PRic_temp_am)
```

    ##                              Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                     2 23433194 11716597  16.621 0.000262 ***
    ## temperature_median            1    73004    73004   0.104 0.752712    
    ## Site_type:temperature_median  2  1377325   688662   0.977 0.402473    
    ## Residuals                    13  9164076   704929                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_sal_am)
```

    ##                           Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                  2 23433194 11716597  18.803 0.000146 ***
    ## salinity_median            1   781794   781794   1.255 0.282936    
    ## Site_type:salinity_median  2  1732017   866008   1.390 0.283800    
    ## Residuals                 13  8100594   623123                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_oxy_am)
```

    ##                         Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                2 23433194 11716597  19.867 0.000111 ***
    ## oxygen_median            1   151334   151334   0.257 0.620941    
    ## Site_type:oxygen_median  2  2796399  1398199   2.371 0.132483    
    ## Residuals               13  7666672   589744                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_temp_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_sal_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_oxy_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.002355813 1.000000000 1.000000000 0.001310863 1.000000000 1.000000000
    ## [7] 0.001002898 1.000000000 1.000000000

``` r
## Mixed and stratified lakes
# Temperature
PD_PRic_MS_lm_temp <- lm(pd.obs ~ Site_type * temperature_median, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-16.png)<!-- -->

``` r
PD_PRic_MS_temp_am <- aov(pd.obs ~ Site_type * temperature_median, data = stree_sespd_env[mixed_stratified_lakes,])
# Salinity
PD_PRic_MS_lm_sal <- lm(pd.obs ~ Site_type * salinity_median, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-17.png)<!-- -->

``` r
PD_PRic_MS_sal_am <- aov(pd.obs ~ Site_type * salinity_median, data = stree_sespd_env[mixed_stratified_lakes,])
# Oxygen
PD_PRic_MS_lm_oxy <- lm(pd.obs ~ Site_type * oxygen_median, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-18.png)<!-- -->

``` r
PD_PRic_MS_oxy_am <- aov(pd.obs ~ Site_type * oxygen_median, data = stree_sespd_env[mixed_stratified_lakes,])
# Anova outputs
summary(PD_PRic_MS_temp_am)
```

    ##                              Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                     1 20056377 20056377  26.274 0.000251 ***
    ## temperature_median            1   195004   195004   0.255 0.622416    
    ## Site_type:temperature_median  1   900909   900909   1.180 0.298659    
    ## Residuals                    12  9160361   763363                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_MS_sal_am)
```

    ##                           Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                  1 20056377 20056377  30.759 0.000127 ***
    ## salinity_median            1   771974   771974   1.184 0.297930    
    ## Site_type:salinity_median  1  1659724  1659724   2.545 0.136600    
    ## Residuals                 12  7824577   652048                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_MS_oxy_am)
```

    ##                         Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                1 20056377 20056377   31.64 0.000112 ***
    ## oxygen_median            1    82535    82535    0.13 0.724495    
    ## Site_type:oxygen_median  1  2566988  2566988    4.05 0.067183 .  
    ## Residuals               12  7606751   633896                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_MS_temp_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_MS_sal_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_MS_oxy_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.002257563 1.000000000 1.000000000 0.001139070 1.000000000 1.000000000
    ## [7] 0.001004657 1.000000000 0.604643861

``` r
## Ocean sites and mixed lakes
# Temperature
PD_PRic_OM_lm_temp <- lm(pd.obs ~ Site_type * temperature_median, data = stree_sespd_env[ocean_mixed_sites_env,])
plot(PD_PRic_OM_lm_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-19.png)<!-- -->

``` r
PD_PRic_OM_temp_am <- aov(pd.obs ~ Site_type * temperature_median, data = stree_sespd_env[ocean_mixed_sites_env,])
# Salinity
PD_PRic_OM_lm_sal <- lm(pd.obs ~ Site_type * salinity_median, data = stree_sespd_env[ocean_mixed_sites_env,])
plot(PD_PRic_OM_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-20.png)<!-- -->

``` r
PD_PRic_OM_sal_am <- aov(pd.obs ~ Site_type * salinity_median, data = stree_sespd_env[ocean_mixed_sites_env,])
# Oxygen
PD_PRic_OM_lm_oxy <- lm(pd.obs ~ Site_type * oxygen_median, data = stree_sespd_env[ocean_mixed_sites_env,])
plot(PD_PRic_OM_lm_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-21.png)<!-- -->

``` r
PD_PRic_OM_oxy_am <- aov(pd.obs ~ Site_type * oxygen_median, data = stree_sespd_env[ocean_mixed_sites_env,])
# Anova outputs
summary(PD_PRic_OM_temp_am)
```

    ##                              Df  Sum Sq Mean Sq F value Pr(>F)
    ## Site_type                     1    2912    2912   0.002  0.962
    ## temperature_median            1  229779  229779   0.192  0.674
    ## Site_type:temperature_median  1 1216000 1216000   1.016  0.347
    ## Residuals                     7 8379005 1197001

``` r
summary(PD_PRic_OM_sal_am)
```

    ##                           Df  Sum Sq Mean Sq F value Pr(>F)
    ## Site_type                  1    2912    2912   0.003  0.961
    ## salinity_median            1 1994704 1994704   1.784  0.223
    ## Site_type:salinity_median  1    2836    2836   0.003  0.961
    ## Residuals                  7 7827244 1118178

``` r
summary(PD_PRic_OM_oxy_am)
```

    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## Site_type                1    2912    2912   0.003  0.959
    ## oxygen_median            1 2602706 2602706   2.523  0.156
    ## Site_type:oxygen_median  1    1412    1412   0.001  0.972
    ## Residuals                7 7220665 1031524

``` r
# p-values
temp_p_values <- summary(PD_PRic_OM_temp_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_OM_sal_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_OM_oxy_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1 1 1 1 1 1 1 1 1

``` r
## Stratified lakes and ocean sites
# Temperature
PD_PRic_SO_lm_temp <- lm(pd.obs ~ Site_type * temperature_median, data = stree_sespd_env[ocean_stratified_sites_env,])
plot(PD_PRic_SO_lm_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-22.png)<!-- -->

``` r
PD_PRic_SO_temp_am <- aov(pd.obs ~ Site_type * temperature_median, data = stree_sespd_env[ocean_stratified_sites_env,])
# Salinity
PD_PRic_SO_lm_sal <- lm(pd.obs ~ Site_type * salinity_median, data = stree_sespd_env[ocean_stratified_sites_env,])
plot(PD_PRic_SO_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-23.png)<!-- -->

``` r
PD_PRic_SO_sal_am <- aov(pd.obs ~ Site_type * salinity_median, data = stree_sespd_env[ocean_stratified_sites_env,])
# Oxygen
PD_PRic_SO_lm_oxy <- lm(pd.obs ~ Site_type * oxygen_median, data = stree_sespd_env[ocean_stratified_sites_env,])
plot(PD_PRic_SO_lm_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-24.png)<!-- -->

``` r
PD_PRic_SO_oxy_am <- aov(pd.obs ~ Site_type * oxygen_median, data = stree_sespd_env[ocean_stratified_sites_env,])
# Anova outputs
summary(PD_PRic_SO_temp_am)
```

    ##                              Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                     1 11299693 11299693 100.278 2.12e-05 ***
    ## temperature_median            1    10349    10349   0.092    0.771    
    ## Site_type:temperature_median  1   348618   348618   3.094    0.122    
    ## Residuals                     7   788785   112684                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_SO_sal_am)
```

    ##                           Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                  1 11299693 11299693 143.980 6.36e-06 ***
    ## salinity_median            1   524411   524411   6.682   0.0362 *  
    ## Site_type:salinity_median  1    73973    73973   0.943   0.3640    
    ## Residuals                  7   549367    78481                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_SO_oxy_am)
```

    ##                         Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                1 11299693 11299693 156.342 4.82e-06 ***
    ## oxygen_median            1   198922   198922   2.752   0.1411    
    ## Site_type:oxygen_median  1   442904   442904   6.128   0.0425 *  
    ## Residuals                7   505927    72275                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_SO_temp_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_SO_sal_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_SO_oxy_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.907936e-04 1.000000e+00 1.000000e+00 5.725162e-05 3.258875e-01
    ## [6] 1.000000e+00 4.340523e-05 1.000000e+00 3.824133e-01

``` r
## Ocean sites
PD_PRic_env_O_lm <- lm(pd.obs ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_sites,])
summary(PD_PRic_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)          -46715        NaN     NaN      NaN
    ## salinity_median        1315        NaN     NaN      NaN
    ## oxygen_median          1087        NaN     NaN      NaN
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
## Mixed lakes
PD_PRic_env_M_lm <- lm(pd.obs ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[mixed_lakes,])
summary(PD_PRic_env_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##    FLK    HLO    LLN    MLN    NLN    NLU    OLO    ULN 
    ##  109.4 -791.4 1111.9 -404.1 1547.2 -968.6 -850.6  246.2 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         -6360.2    58511.4  -0.109    0.919
    ## salinity_median       533.5     1360.3   0.392    0.715
    ## oxygen_median         945.1     1100.3   0.859    0.439
    ## temperature_median   -407.3      844.3  -0.482    0.655
    ## 
    ## Residual standard error: 1240 on 4 degrees of freedom
    ## Multiple R-squared:  0.3499, Adjusted R-squared:  -0.1376 
    ## F-statistic: 0.7177 on 3 and 4 DF,  p-value: 0.5912

``` r
Anova(PD_PRic_env_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                     Sum Sq Df F value Pr(>F)
    ## (Intercept)          18179  1  0.0118 0.9187
    ## salinity_median     236683  1  0.1538 0.7149
    ## oxygen_median      1135172  1  0.7378 0.4388
    ## temperature_median  358075  1  0.2327 0.6547
    ## Residuals          6154001  4

``` r
p_values <- summary(PD_PRic_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
## Stratified lakes
PD_PRic_env_S_lm <- lm(pd.obs ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[stratified_lakes,])
summary(PD_PRic_env_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##     BCM     CLM     GLK     HLM     NLK     OTM     SLN     TLN 
    ##  212.22   56.93   61.76  181.23 -243.79 -297.72 -114.92  144.28 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -1322.38    2856.87  -0.463    0.668
    ## salinity_median       72.51      44.37   1.634    0.178
    ## oxygen_median        -37.31     195.61  -0.191    0.858
    ## temperature_median    13.92      63.28   0.220    0.837
    ## 
    ## Residual standard error: 258.4 on 4 degrees of freedom
    ## Multiple R-squared:  0.6618, Adjusted R-squared:  0.4082 
    ## F-statistic: 2.609 on 3 and 4 DF,  p-value: 0.1885

``` r
Anova(PD_PRic_env_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                    Sum Sq Df F value Pr(>F)
    ## (Intercept)         14303  1  0.2143 0.6675
    ## salinity_median    178300  1  2.6709 0.1775
    ## oxygen_median        2429  1  0.0364 0.8580
    ## temperature_median   3233  1  0.0484 0.8366
    ## Residuals          267028  4

``` r
p_values <- summary(PD_PRic_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          0.7101534          1.0000000          1.0000000

``` r
### PDisp ANOVA
## Surveyed sites
# Temperature
PD_PRic_lm_temp <- lm(Dispersion ~ Site_type * temperature_median, data = stree_sespd_env[surveyed_sites_env,])
plot(PD_PRic_lm_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-25.png)<!-- -->

``` r
PD_PRic_temp_am <- aov(Dispersion ~ Site_type * temperature_median, data = stree_sespd_env[surveyed_sites_env,])
# Salinity
PD_PRic_lm_sal <- lm(Dispersion ~ Site_type * salinity_median, data = stree_sespd_env[surveyed_sites_env,])
plot(PD_PRic_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-26.png)<!-- -->

``` r
PD_PRic_sal_am <- aov(Dispersion ~ Site_type * salinity_median, data = stree_sespd_env[surveyed_sites_env,])
# Oxygen
PD_PRic_lm_oxy <- lm(Dispersion ~ Site_type * oxygen_median, data = stree_sespd_env[surveyed_sites_env,])
plot(PD_PRic_lm_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-27.png)<!-- -->

``` r
PD_PRic_oxy_am <- aov(Dispersion ~ Site_type * oxygen_median, data = stree_sespd_env[surveyed_sites_env,])
# Anova outputs
summary(PD_PRic_temp_am)
```

    ##                              Df   Sum Sq   Mean Sq F value Pr(>F)  
    ## Site_type                     2 0.000740 0.0003701   0.717 0.5063  
    ## temperature_median            1 0.001896 0.0018956   3.675 0.0775 .
    ## Site_type:temperature_median  2 0.000613 0.0003063   0.594 0.5666  
    ## Residuals                    13 0.006706 0.0005158                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_sal_am)
```

    ##                           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                  2 0.000740 0.0003701   0.545  0.593
    ## salinity_median            1 0.000353 0.0003531   0.520  0.484
    ## Site_type:salinity_median  2 0.000030 0.0000148   0.022  0.979
    ## Residuals                 13 0.008832 0.0006793

``` r
summary(PD_PRic_oxy_am)
```

    ##                         Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                2 0.000740 0.0003701   0.532  0.600
    ## oxygen_median            1 0.000056 0.0000561   0.081  0.781
    ## Site_type:oxygen_median  2 0.000106 0.0000530   0.076  0.927
    ## Residuals               13 0.009052 0.0006963

``` r
# p-values
temp_p_values <- summary(PD_PRic_temp_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_sal_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_oxy_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 0.6973644 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [8] 1.0000000 1.0000000

``` r
## Mixed and stratified lakes
# Temperature
PD_PRic_MS_lm_temp <- lm(Dispersion ~ Site_type * temperature_median, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-28.png)<!-- -->

``` r
PD_PRic_MS_temp_am <- aov(Dispersion ~ Site_type * temperature_median, data = stree_sespd_env[mixed_stratified_lakes,])
# Salinity
PD_PRic_MS_lm_sal <- lm(Dispersion ~ Site_type * salinity_median, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-29.png)<!-- -->

``` r
PD_PRic_MS_sal_am <- aov(Dispersion ~ Site_type * salinity_median, data = stree_sespd_env[mixed_stratified_lakes,])
# Oxygen
PD_PRic_MS_lm_oxy <- lm(Dispersion ~ Site_type * oxygen_median, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-30.png)<!-- -->

``` r
PD_PRic_MS_oxy_am <- aov(Dispersion ~ Site_type * oxygen_median, data = stree_sespd_env[mixed_stratified_lakes,])
# Anova outputs
summary(PD_PRic_MS_temp_am)
```

    ##                              Df   Sum Sq   Mean Sq F value Pr(>F)  
    ## Site_type                     1 0.000707 0.0007066   1.269 0.2819  
    ## temperature_median            1 0.002287 0.0022868   4.108 0.0655 .
    ## Site_type:temperature_median  1 0.000122 0.0001220   0.219 0.6481  
    ## Residuals                    12 0.006680 0.0005567                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_MS_sal_am)
```

    ##                           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                  1 0.000707 0.0007066   0.974  0.343
    ## salinity_median            1 0.000353 0.0003527   0.486  0.499
    ## Site_type:salinity_median  1 0.000029 0.0000293   0.040  0.844
    ## Residuals                 12 0.008707 0.0007256

``` r
summary(PD_PRic_MS_oxy_am)
```

    ##                         Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                1 0.000707 0.0007066   0.937  0.352
    ## oxygen_median            1 0.000029 0.0000294   0.039  0.847
    ## Site_type:oxygen_median  1 0.000010 0.0000104   0.014  0.908
    ## Residuals               12 0.009049 0.0007541

``` r
# p-values
temp_p_values <- summary(PD_PRic_MS_temp_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_MS_sal_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_MS_oxy_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 0.5894609 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [8] 1.0000000 1.0000000

``` r
## Ocean sites and mixed lakes
# Temperature
PD_PRic_OM_lm_temp <- lm(Dispersion ~ Site_type * temperature_median, data = stree_sespd_env[ocean_mixed_sites_env,])
plot(PD_PRic_OM_lm_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-31.png)<!-- -->

``` r
PD_PRic_OM_temp_am <- aov(Dispersion ~ Site_type * temperature_median, data = stree_sespd_env[ocean_mixed_sites_env,])
# Salinity
PD_PRic_OM_lm_sal <- lm(Dispersion ~ Site_type * salinity_median, data = stree_sespd_env[ocean_mixed_sites_env,])
plot(PD_PRic_OM_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-32.png)<!-- -->

``` r
PD_PRic_OM_sal_am <- aov(Dispersion ~ Site_type * salinity_median, data = stree_sespd_env[ocean_mixed_sites_env,])
# Oxygen
PD_PRic_OM_lm_oxy <- lm(Dispersion ~ Site_type * oxygen_median, data = stree_sespd_env[ocean_mixed_sites_env,])
plot(PD_PRic_OM_lm_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-33.png)<!-- -->

``` r
PD_PRic_OM_oxy_am <- aov(Dispersion ~ Site_type * oxygen_median, data = stree_sespd_env[ocean_mixed_sites_env,])
# Anova outputs
summary(PD_PRic_OM_temp_am)
```

    ##                              Df    Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                     1 0.0002310 2.311e-04   2.901  0.132
    ## temperature_median            1 0.0000001 1.100e-07   0.001  0.971
    ## Site_type:temperature_median  1 0.0001488 1.489e-04   1.869  0.214
    ## Residuals                     7 0.0005576 7.966e-05

``` r
summary(PD_PRic_OM_sal_am)
```

    ##                           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                  1 2.31e-04 2.311e-04   2.324  0.171
    ## salinity_median            1 8.70e-06 8.730e-06   0.088  0.776
    ## Site_type:salinity_median  1 1.80e-06 1.790e-06   0.018  0.897
    ## Residuals                  7 6.96e-04 9.943e-05

``` r
summary(PD_PRic_OM_oxy_am)
```

    ##                         Df    Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                1 0.0002310 2.311e-04   2.924  0.131
    ## oxygen_median            1 0.0000843 8.432e-05   1.067  0.336
    ## Site_type:oxygen_median  1 0.0000691 6.905e-05   0.874  0.381
    ## Residuals                7 0.0005532 7.903e-05

``` r
# p-values
temp_p_values <- summary(PD_PRic_OM_temp_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_OM_sal_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_OM_oxy_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1 1 1 1 1 1 1 1 1

``` r
## Stratified lakes and ocean sites
# Temperature
PD_PRic_SO_lm_temp <- lm(Dispersion ~ Site_type * temperature_median, data = stree_sespd_env[ocean_stratified_sites_env,])
plot(PD_PRic_SO_lm_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-34.png)<!-- -->

``` r
PD_PRic_SO_temp_am <- aov(Dispersion ~ Site_type * temperature_median, data = stree_sespd_env[ocean_stratified_sites_env,])
# Salinity
PD_PRic_SO_lm_sal <- lm(Dispersion ~ Site_type * salinity_median, data = stree_sespd_env[ocean_stratified_sites_env,])
plot(PD_PRic_SO_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-35.png)<!-- -->

``` r
PD_PRic_SO_sal_am <- aov(Dispersion ~ Site_type * salinity_median, data = stree_sespd_env[ocean_stratified_sites_env,])
# Oxygen
PD_PRic_SO_lm_oxy <- lm(Dispersion ~ Site_type * oxygen_median, data = stree_sespd_env[ocean_stratified_sites_env,])
plot(PD_PRic_SO_lm_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-36.png)<!-- -->

``` r
PD_PRic_SO_oxy_am <- aov(Dispersion ~ Site_type * oxygen_median, data = stree_sespd_env[ocean_stratified_sites_env,])
# Anova outputs
summary(PD_PRic_SO_temp_am)
```

    ##                              Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                     1 0.000020 0.0000196   0.022  0.886
    ## temperature_median            1 0.001922 0.0019225   2.180  0.183
    ## Site_type:temperature_median  1 0.000536 0.0005361   0.608  0.461
    ## Residuals                     7 0.006174 0.0008820

``` r
summary(PD_PRic_SO_sal_am)
```

    ##                           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                  1 0.000020 0.0000196   0.017  0.901
    ## salinity_median            1 0.000373 0.0003727   0.316  0.592
    ## Site_type:salinity_median  1 0.000000 0.0000002   0.000  0.990
    ## Residuals                  7 0.008260 0.0011800

``` r
summary(PD_PRic_SO_oxy_am)
```

    ##                         Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                1 0.000020 0.0000196   0.016  0.902
    ## oxygen_median            1 0.000029 0.0000292   0.024  0.881
    ## Site_type:oxygen_median  1 0.000102 0.0001019   0.084  0.781
    ## Residuals                7 0.008502 0.0012145

``` r
# p-values
temp_p_values <- summary(PD_PRic_SO_temp_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_SO_sal_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_SO_oxy_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1 1 1 1 1 1 1 1 1

``` r
## Ocean sites
PD_PRic_env_O_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_sites,])
summary(PD_PRic_env_O_lm)
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
p_values <- summary(PD_PRic_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
## Mixed lakes
PD_PRic_env_M_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[mixed_lakes,])
summary(PD_PRic_env_M_lm)
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
Anova(PD_PRic_env_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                        Sum Sq Df F value Pr(>F)
    ## (Intercept)        0.00000012  1  0.0010 0.9766
    ## salinity_median    0.00000418  1  0.0332 0.8643
    ## oxygen_median      0.00002720  1  0.2157 0.6665
    ## temperature_median 0.00004515  1  0.3579 0.5819
    ## Residuals          0.00050452  4

``` r
p_values <- summary(PD_PRic_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
## Stratified lakes
PD_PRic_env_S_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[stratified_lakes,])
summary(PD_PRic_env_S_lm)
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
Anova(PD_PRic_env_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                       Sum Sq Df F value Pr(>F)
    ## (Intercept)        0.0019793  1  1.6642 0.2666
    ## salinity_median    0.0013499  1  1.1349 0.3468
    ## oxygen_median      0.0005808  1  0.4883 0.5232
    ## temperature_median 0.0030397  1  2.5557 0.1851
    ## Residuals          0.0047575  4

``` r
p_values <- summary(PD_PRic_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          1.0000000          1.0000000          0.7405587

``` r
#### Geographical
### logPRic ANOVA
## Surveyed sites
# Distance to the ocean
PD_z_lm_dist <- lm(pd.obs.z ~  Site_type * distance_to_ocean_min_m, data = stree_sespd_env[surveyed_sites,])
plot(PD_z_lm_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   19

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-37.png)<!-- -->

``` r
PD_z_am_dist <- aov(pd.obs.z ~  Site_type * distance_to_ocean_min_m, data = stree_sespd_env[surveyed_sites,])
# Max depth
PD_z_lm_mxd <- lm(pd.obs.z ~  Site_type * max_depth, data = stree_sespd_env[surveyed_sites,])
plot(PD_z_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-38.png)<!-- -->

``` r
PD_z_am_mxd <- aov(pd.obs.z ~  Site_type * max_depth, data = stree_sespd_env[surveyed_sites,])
# Log Area
PD_z_lm_lga <- lm(pd.obs.z ~  Site_type * logArea, data = stree_sespd_env[surveyed_sites,])
plot(PD_z_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-39.png)<!-- -->

``` r
PD_z_am_lga <- aov(pd.obs.z ~  Site_type * logArea, data = stree_sespd_env[surveyed_sites,])
# Anova outputs
summary(PD_z_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type                          2  1.838  0.9190   0.807  0.463
    ## distance_to_ocean_min_m            1  0.228  0.2284   0.201  0.660
    ## Site_type:distance_to_ocean_min_m  2  1.224  0.6118   0.537  0.594
    ## Residuals                         16 18.216  1.1385

``` r
summary(PD_z_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value  Pr(>F)   
    ## Site_type            2  1.838   0.919   1.636 0.22577   
    ## max_depth            1  2.719   2.719   4.839 0.04287 * 
    ## Site_type:max_depth  2  7.960   3.980   7.083 0.00626 **
    ## Residuals           16  8.990   0.562                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_z_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type          2  1.838  0.9190   0.916  0.420
    ## logArea            1  0.311  0.3114   0.310  0.585
    ## Site_type:logArea  2  3.297  1.6485   1.642  0.225
    ## Residuals         16 16.060  1.0037

``` r
# p-values
temp_p_values <- summary(PD_z_am_dist)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_z_am_mxd)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_z_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.00000000 1.00000000 1.00000000 1.00000000 0.38583104 0.05637211 1.00000000
    ## [8] 1.00000000 1.00000000

``` r
## Mixed and stratified lakes
# Distance to the ocean
PD_z_MS_lm_dist <- lm(pd.obs.z ~  Site_type * distance_to_ocean_min_m, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_z_MS_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-40.png)<!-- -->

``` r
PD_z_MS_am_dist <- aov(pd.obs.z ~  Site_type * distance_to_ocean_min_m, data = stree_sespd_env[mixed_stratified_lakes,])
# Max depth
PD_z_MS_lm_mxd <- lm(pd.obs.z ~  Site_type * max_depth, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_z_MS_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-41.png)<!-- -->

``` r
PD_z_MS_am_mxd <- aov(pd.obs.z ~  Site_type * max_depth, data = stree_sespd_env[mixed_stratified_lakes,])
# Log Area
PD_z_MS_lm_lga <- lm(pd.obs.z ~  Site_type * logArea, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_z_MS_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-42.png)<!-- -->

``` r
PD_z_MS_am_lga <- aov(pd.obs.z ~  Site_type * logArea, data = stree_sespd_env[mixed_stratified_lakes,])
# Anova outputs
summary(PD_z_MS_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type                          1  0.000  0.0000   0.000  0.999
    ## distance_to_ocean_min_m            1  0.331  0.3311   0.427  0.526
    ## Site_type:distance_to_ocean_min_m  1  0.171  0.1708   0.220  0.647
    ## Residuals                         12  9.302  0.7752

``` r
summary(PD_z_MS_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type            1  0.000  0.0000   0.000 0.9994  
    ## max_depth            1  0.004  0.0039   0.006 0.9388  
    ## Site_type:max_depth  1  2.139  2.1389   3.350 0.0921 .
    ## Residuals           12  7.661  0.6384                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_z_MS_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type          1  0.000  0.0000   0.000 0.9994  
    ## logArea            1  0.194  0.1940   0.312 0.5869  
    ## Site_type:logArea  1  2.141  2.1410   3.440 0.0884 .
    ## Residuals         12  7.469  0.6224                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_z_MS_am_dist)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_z_MS_am_mxd)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_z_MS_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 0.8291696 1.0000000
    ## [8] 1.0000000 0.7952892

``` r
## Ocean sites and mixed lakes
# Distance to the ocean
PD_z_OM_lm_dist <- lm(pd.obs.z ~  Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_mixed_sites,])
plot(PD_z_OM_lm_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   6

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-43.png)<!-- -->

``` r
PD_z_OM_am_dist <- aov(pd.obs.z ~  Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_mixed_sites,])
# Max depth
PD_z_OM_lm_mxd <- lm(pd.obs.z ~  Site_type * max_depth, data = stree_sespd_env[ocean_mixed_sites,])
plot(PD_z_OM_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-44.png)<!-- -->

``` r
PD_z_OM_am_mxd <- aov(pd.obs.z ~  Site_type * max_depth, data = stree_sespd_env[ocean_mixed_sites,])
# Log Area
PD_z_OM_lm_lga <- lm(pd.obs.z ~  Site_type * logArea, data = stree_sespd_env[ocean_mixed_sites,])
plot(PD_z_OM_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-45.png)<!-- -->

``` r
PD_z_OM_am_lga <- aov(pd.obs.z ~  Site_type * logArea, data = stree_sespd_env[ocean_mixed_sites,])
# Anova outputs
summary(PD_z_OM_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type                          1  1.444  1.4436   1.072  0.325
    ## distance_to_ocean_min_m            1  0.118  0.1182   0.088  0.773
    ## Site_type:distance_to_ocean_min_m  1  1.195  1.1950   0.887  0.368
    ## Residuals                         10 13.469  1.3469

``` r
summary(PD_z_OM_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value  Pr(>F)   
    ## Site_type            1  1.444   1.444   3.010 0.11339   
    ## max_depth            1  8.785   8.785  18.319 0.00161 **
    ## Site_type:max_depth  1  1.203   1.203   2.509 0.14431   
    ## Residuals           10  4.795   0.480                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_z_OM_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type          1  1.444  1.4436   1.095   0.32
    ## logArea            1  1.577  1.5768   1.196   0.30
    ## Site_type:logArea  1  0.026  0.0264   0.020   0.89
    ## Residuals         10 13.180  1.3180

``` r
# p-values
temp_p_values <- summary(PD_z_OM_am_dist)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_z_OM_am_mxd)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_z_OM_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000 1.0000000 1.0000000 0.0144999 1.0000000 1.0000000
    ## [8] 1.0000000 1.0000000

``` r
## Stratified lakes and ocean sites
# Distance to the ocean
PD_z_SO_lm_dist <- lm(pd.obs.z ~  Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_stratified_sites,])
plot(PD_z_SO_lm_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   6

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-46.png)<!-- -->

``` r
PD_z_SO_am_dist <- aov(pd.obs.z ~  Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_stratified_sites,])
# Max depth
PD_z_SO_lm_mxd <- lm(pd.obs.z ~  Site_type * max_depth, data = stree_sespd_env[ocean_stratified_sites,])
plot(PD_z_SO_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-47.png)<!-- -->

``` r
PD_z_SO_am_mxd <- aov(pd.obs.z ~  Site_type * max_depth, data = stree_sespd_env[ocean_stratified_sites,])
# Log Area
PD_z_SO_lm_lga <- lm(pd.obs.z ~  Site_type * logArea, data = stree_sespd_env[ocean_stratified_sites,])
plot(PD_z_SO_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-48.png)<!-- -->

``` r
PD_z_SO_am_lga <- aov(pd.obs.z ~  Site_type * logArea, data = stree_sespd_env[ocean_stratified_sites,])
# Anova outputs
summary(PD_z_SO_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type                          1  1.445  1.4449   1.058  0.328
    ## distance_to_ocean_min_m            1  0.072  0.0717   0.053  0.823
    ## Site_type:distance_to_ocean_min_m  1  1.017  1.0171   0.744  0.408
    ## Residuals                         10 13.661  1.3661

``` r
summary(PD_z_SO_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value  Pr(>F)   
    ## Site_type            1  1.445   1.445   2.616 0.13688   
    ## max_depth            1  1.453   1.453   2.631 0.13585   
    ## Site_type:max_depth  1  7.773   7.773  14.072 0.00378 **
    ## Residuals           10  5.524   0.552                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_z_SO_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type          1  1.445   1.445   1.260  0.288
    ## logArea            1  0.094   0.094   0.082  0.781
    ## Site_type:logArea  1  3.185   3.185   2.776  0.127
    ## Residuals         10 11.471   1.147

``` r
# p-values
temp_p_values <- summary(PD_z_SO_am_dist)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_z_SO_am_mxd)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_z_SO_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 0.03398504 1.00000000
    ## [8] 1.00000000 1.00000000

``` r
## Ocean sites
PD_z_geo_O_lm <- lm(pd.obs.z ~ max_depth + logArea + distance_to_ocean_min_m, stree_sespd_env[ocean_sites,])
summary(PD_z_geo_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ max_depth + logArea + distance_to_ocean_min_m, 
    ##     data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ##        IBK        LCN        NCN        OOO        OOM        RCA 
    ##  1.413e-01  7.442e-01 -4.173e-01  9.104e-02 -5.592e-01 -5.696e-17 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              1.21947    1.87239   0.651   0.5817  
    ## max_depth               -0.09940    0.03360  -2.959   0.0978 .
    ## logArea                 -0.09907    0.16662  -0.595   0.6124  
    ## distance_to_ocean_min_m -0.02868    0.04447  -0.645   0.5851  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7311 on 2 degrees of freedom
    ## Multiple R-squared:  0.8916, Adjusted R-squared:  0.7291 
    ## F-statistic: 5.485 on 3 and 2 DF,  p-value: 0.1581

``` r
Anova(PD_z_geo_O_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                         Sum Sq Df F value  Pr(>F)  
    ## (Intercept)             0.2267  1  0.4242 0.58169  
    ## max_depth               4.6790  1  8.7535 0.09777 .
    ## logArea                 0.1890  1  0.3535 0.61243  
    ## distance_to_ocean_min_m 0.2222  1  0.4158 0.58514  
    ## Residuals               1.0691  2                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_geo_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept)               max_depth                 logArea 
    ##               1.0000000               0.3910896               1.0000000 
    ## distance_to_ocean_min_m 
    ##               1.0000000

``` r
## Mixed lakes
PD_z_geo_M_lm <- lm(pd.obs.z ~ max_depth + logArea + distance_to_ocean_min_m, stree_sespd_env[mixed_lakes,])
summary(PD_z_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ max_depth + logArea + distance_to_ocean_min_m, 
    ##     data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##     FLK     HLO     LLN     MLN     NLN     NLU     OLO     ULN 
    ## -0.2123 -0.7466 -0.5777  0.4045  1.3559 -0.0301  0.3621 -0.5558 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)             -0.746784   2.562853  -0.291    0.785
    ## max_depth               -0.066912   0.066221  -1.010    0.369
    ## logArea                  0.106039   0.339921   0.312    0.771
    ## distance_to_ocean_min_m -0.003258   0.015157  -0.215    0.840
    ## 
    ## Residual standard error: 0.9192 on 4 degrees of freedom
    ## Multiple R-squared:  0.3129, Adjusted R-squared:  -0.2024 
    ## F-statistic: 0.6072 on 3 and 4 DF,  p-value: 0.6446

``` r
Anova(PD_z_geo_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                         Sum Sq Df F value Pr(>F)
    ## (Intercept)             0.0717  1  0.0849 0.7852
    ## max_depth               0.8626  1  1.0210 0.3694
    ## logArea                 0.0822  1  0.0973 0.7707
    ## distance_to_ocean_min_m 0.0390  1  0.0462 0.8403
    ## Residuals               3.3794  4

``` r
p_values <- summary(PD_z_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept)               max_depth                 logArea 
    ##                       1                       1                       1 
    ## distance_to_ocean_min_m 
    ##                       1

``` r
## Stratified lakes
PD_z_geo_S_lm <- lm(pd.obs.z ~ max_depth + logArea + distance_to_ocean_min_m, stree_sespd_env[stratified_lakes,])
summary(PD_z_geo_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ max_depth + logArea + distance_to_ocean_min_m, 
    ##     data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##      BCM      CLM      GLK      HLM      NLK      OTM      SLN      TLN 
    ## -0.12449 -0.24694  0.59937 -0.58482  0.07139  0.72921  0.18406 -0.62777 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -10.973624   4.064548  -2.700   0.0541 .
    ## max_depth                -0.079262   0.049056  -1.616   0.1815  
    ## logArea                   1.180681   0.502819   2.348   0.0787 .
    ## distance_to_ocean_min_m  -0.000481   0.003454  -0.139   0.8960  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.66 on 4 degrees of freedom
    ## Multiple R-squared:  0.6433, Adjusted R-squared:  0.3758 
    ## F-statistic: 2.405 on 3 and 4 DF,  p-value: 0.2079

``` r
Anova(PD_z_geo_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                         Sum Sq Df F value  Pr(>F)  
    ## (Intercept)             3.1755  1  7.2891 0.05410 .
    ## max_depth               1.1373  1  2.6106 0.18145  
    ## logArea                 2.4020  1  5.5137 0.07868 .
    ## distance_to_ocean_min_m 0.0084  1  0.0194 0.89597  
    ## Residuals               1.7426  4                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_z_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept)               max_depth                 logArea 
    ##               0.2164135               0.7258090               0.3147069 
    ## distance_to_ocean_min_m 
    ##               1.0000000

``` r
### PRic ANOVA
## Surveyed sites
# Distance to the ocean
PD_PRic_lm_dist <- lm(pd.obs ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[surveyed_sites,])
plot(PD_PRic_lm_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   19

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-49.png)<!-- -->

``` r
PD_PRic_dist_am <- aov(pd.obs ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[surveyed_sites,])
# Max depth
PD_PRic_lm_mxd <- lm(pd.obs ~ Site_type * max_depth, data = stree_sespd_env[surveyed_sites,])
plot(PD_PRic_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-50.png)<!-- -->

``` r
PD_PRic_mxd_am <- aov(pd.obs ~ Site_type * max_depth, data = stree_sespd_env[surveyed_sites,])
# Log Area
PD_PRic_lm_lga <- lm(pd.obs ~ Site_type * logArea, data = stree_sespd_env[surveyed_sites,])
plot(PD_PRic_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-51.png)<!-- -->

``` r
PD_PRic_lga_am <- aov(pd.obs ~ Site_type * logArea, data = stree_sespd_env[surveyed_sites,])
# Anova outputs
summary(PD_PRic_dist_am)
```

    ##                                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                          2 23597133 11798567  15.029 0.000212 ***
    ## distance_to_ocean_min_m            1   447734   447734   0.570 0.461104    
    ## Site_type:distance_to_ocean_min_m  2    35606    17803   0.023 0.977609    
    ## Residuals                         16 12560885   785055                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_mxd_am)
```

    ##                     Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type            2 23597133 11798567  19.671 4.88e-05 ***
    ## max_depth            1  1698378  1698378   2.832    0.112    
    ## Site_type:max_depth  2  1749289   874645   1.458    0.262    
    ## Residuals           16  9596558   599785                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_lga_am)
```

    ##                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type          2 23597133 11798567  30.727 3.32e-06 ***
    ## logArea            1  2405280  2405280   6.264   0.0235 *  
    ## Site_type:logArea  2  4495320  2247660   5.854   0.0124 *  
    ## Residuals         16  6143625   383977                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_dist_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_mxd_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_lga_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.908831e-03 1.000000e+00 1.000000e+00 4.392604e-04 1.000000e+00
    ## [6] 1.000000e+00 2.984171e-05 2.118605e-01 1.112888e-01

``` r
## Mixed and stratified lakes
# Distance to the ocean
PD_PRic_MS_lm_dist <- lm(pd.obs ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-52.png)<!-- -->

``` r
PD_PRic_MS_dist_am <- aov(pd.obs ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[mixed_stratified_lakes,])
# Max depth
PD_PRic_MS_lm_mxd <- lm(pd.obs ~ Site_type * max_depth, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-53.png)<!-- -->

``` r
PD_PRic_MS_mxd_am <- aov(pd.obs ~ Site_type * max_depth, data = stree_sespd_env[mixed_stratified_lakes,])
# Log Area
PD_PRic_MS_lm_lga <- lm(pd.obs ~ Site_type * logArea, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-54.png)<!-- -->

``` r
PD_PRic_MS_lga_am <- aov(pd.obs ~ Site_type * logArea, data = stree_sespd_env[mixed_stratified_lakes,])
# Anova outputs
summary(PD_PRic_MS_dist_am)
```

    ##                                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                          1 20056377 20056377  24.588 0.000332 ***
    ## distance_to_ocean_min_m            1   435891   435891   0.534 0.478801    
    ## Site_type:distance_to_ocean_min_m  1    32123    32123   0.039 0.846018    
    ## Residuals                         12  9788260   815688                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_MS_mxd_am)
```

    ##                     Df   Sum Sq  Mean Sq F value  Pr(>F)    
    ## Site_type            1 20056377 20056377  32.955 9.3e-05 ***
    ## max_depth            1  1203987  1203987   1.978   0.185    
    ## Site_type:max_depth  1  1749176  1749176   2.874   0.116    
    ## Residuals           12  7303111   608593                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_MS_lga_am)
```

    ##                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type          1 20056377 20056377  70.714 2.25e-06 ***
    ## logArea            1  4421808  4421808  15.590  0.00193 ** 
    ## Site_type:logArea  1  2430958  2430958   8.571  0.01266 *  
    ## Residuals         12  3403507   283626                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_MS_dist_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_MS_mxd_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_MS_lga_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 2.983691e-03 1.000000e+00 1.000000e+00 8.368673e-04 1.000000e+00
    ## [6] 1.000000e+00 2.023307e-05 1.739855e-02 1.139366e-01

``` r
## Ocean sites and mixed lakes
# Distance to the ocean
PD_PRic_OM_lm_dist <- lm(pd.obs ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_mixed_sites,])
plot(PD_PRic_OM_lm_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   6

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-55.png)<!-- -->

``` r
PD_PRic_OM_dist_am <- aov(pd.obs ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_mixed_sites,])
# Max depth
PD_PRic_OM_lm_mxd <- lm(pd.obs ~ Site_type * max_depth, data = stree_sespd_env[ocean_mixed_sites,])
plot(PD_PRic_OM_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-56.png)<!-- -->

``` r
PD_PRic_OM_mxd_am <- aov(pd.obs ~ Site_type * max_depth, data = stree_sespd_env[ocean_mixed_sites,])
# Log Area
PD_PRic_OM_lm_lga <- lm(pd.obs ~ Site_type * logArea, data = stree_sespd_env[ocean_mixed_sites,])
plot(PD_PRic_OM_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-57.png)<!-- -->

``` r
PD_PRic_OM_lga_am <- aov(pd.obs ~ Site_type * logArea, data = stree_sespd_env[ocean_mixed_sites,])
# Anova outputs
summary(PD_PRic_OM_dist_am)
```

    ##                                   Df   Sum Sq Mean Sq F value Pr(>F)
    ## Site_type                          1   164165  164165   0.136  0.720
    ## distance_to_ocean_min_m            1   190695  190695   0.158  0.699
    ## Site_type:distance_to_ocean_min_m  1      279     279   0.000  0.988
    ## Residuals                         10 12063630 1206363

``` r
summary(PD_PRic_OM_mxd_am)
```

    ##                     Df  Sum Sq Mean Sq F value Pr(>F)
    ## Site_type            1  164165  164165   0.186  0.675
    ## max_depth            1 2826617 2826617   3.207  0.104
    ## Site_type:max_depth  1  614819  614819   0.698  0.423
    ## Residuals           10 8813167  881317

``` r
summary(PD_PRic_OM_lga_am)
```

    ##                   Df  Sum Sq Mean Sq F value Pr(>F)  
    ## Site_type          1  164165  164165   0.306 0.5923  
    ## logArea            1 2742783 2742783   5.113 0.0473 *
    ## Site_type:logArea  1 4147050 4147050   7.730 0.0194 *
    ## Residuals         10 5364770  536477                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_OM_dist_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_OM_mxd_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_OM_lga_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000 1.0000000 1.0000000 0.9322017 1.0000000 1.0000000
    ## [8] 0.4255192 0.1749611

``` r
## Stratified lakes and ocean sites
# Distance to the ocean
PD_PRic_SO_lm_dist <- lm(pd.obs ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_stratified_sites,])
plot(PD_PRic_SO_lm_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   6

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-58.png)<!-- -->

``` r
PD_PRic_SO_dist_am <- aov(pd.obs ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_stratified_sites,])
# Max depth
PD_PRic_SO_lm_mxd <- lm(pd.obs ~ Site_type * max_depth, data = stree_sespd_env[ocean_stratified_sites,])
plot(PD_PRic_SO_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-59.png)<!-- -->

``` r
PD_PRic_SO_mxd_am <- aov(pd.obs ~ Site_type * max_depth, data = stree_sespd_env[ocean_stratified_sites,])
# Log Area
PD_PRic_SO_lm_lga <- lm(pd.obs ~ Site_type * logArea, data = stree_sespd_env[ocean_stratified_sites,])
plot(PD_PRic_SO_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-60.png)<!-- -->

``` r
PD_PRic_SO_lga_am <- aov(pd.obs ~ Site_type * logArea, data = stree_sespd_env[ocean_stratified_sites,])
# Anova outputs
summary(PD_PRic_SO_dist_am)
```

    ##                                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type                          1 13995470 13995470  42.801 6.54e-05 ***
    ## distance_to_ocean_min_m            1   303308   303308   0.928    0.358    
    ## Site_type:distance_to_ocean_min_m  1     4385     4385   0.013    0.910    
    ## Residuals                         10  3269879   326988                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_SO_mxd_am)
```

    ##                     Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type            1 13995470 13995470  45.487 5.08e-05 ***
    ## max_depth            1   256212   256212   0.833    0.383    
    ## Site_type:max_depth  1   244522   244522   0.795    0.394    
    ## Residuals           10  3076838   307684                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_SO_lga_am)
```

    ##                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## Site_type          1 13995470 13995470  39.771 8.83e-05 ***
    ## logArea            1    58391    58391   0.166    0.692    
    ## Site_type:logArea  1      209      209   0.001    0.981    
    ## Residuals         10  3518972   351897                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_SO_dist_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_SO_mxd_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_SO_lga_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.0005882776 1.0000000000 1.0000000000 0.0004569670 1.0000000000
    ## [6] 1.0000000000 0.0007950846 1.0000000000 1.0000000000

``` r
## Ocean sites
PD_PRic_geo_O_lm <- lm(pd.obs ~ max_depth + logArea + distance_to_ocean_min_m, stree_sespd_env[ocean_sites,])
summary(PD_PRic_geo_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ max_depth + logArea + distance_to_ocean_min_m, 
    ##     data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ##        IBK        LCN        NCN        OOO        OOM        RCA 
    ##  3.468e+01 -3.594e+02  5.753e+01 -8.603e+02  1.127e+03 -2.155e-13 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              3217.41    2652.24   1.213    0.349
    ## max_depth                  35.43      47.59   0.744    0.534
    ## logArea                   -53.45     236.02  -0.226    0.842
    ## distance_to_ocean_min_m   -23.43      63.00  -0.372    0.746
    ## 
    ## Residual standard error: 1036 on 2 degrees of freedom
    ## Multiple R-squared:  0.2306, Adjusted R-squared:  -0.9235 
    ## F-statistic: 0.1998 on 3 and 2 DF,  p-value: 0.8893

``` r
Anova(PD_PRic_geo_O_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                          Sum Sq Df F value Pr(>F)
    ## (Intercept)             1578309  1  1.4716 0.3489
    ## max_depth                594414  1  0.5542 0.5342
    ## logArea                   55014  1  0.0513 0.8419
    ## distance_to_ocean_min_m  148322  1  0.1383 0.7457
    ## Residuals               2145041  2

``` r
p_values <- summary(PD_PRic_geo_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept)               max_depth                 logArea 
    ##                       1                       1                       1 
    ## distance_to_ocean_min_m 
    ##                       1

``` r
## Mixed lakes
PD_PRic_geo_M_lm <- lm(pd.obs ~ max_depth + logArea + distance_to_ocean_min_m, stree_sespd_env[mixed_lakes,])
summary(PD_PRic_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ max_depth + logArea + distance_to_ocean_min_m, 
    ##     data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##      FLK      HLO      LLN      MLN      NLN      NLU      OLO      ULN 
    ##   209.26  -644.99   653.52    91.81   242.45   478.25 -1131.39   101.09 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -3042.901   2192.435  -1.388   0.2375  
    ## max_depth                 -11.498     56.650  -0.203   0.8491  
    ## logArea                   707.357    290.791   2.433   0.0718 .
    ## distance_to_ocean_min_m    -6.087     12.967  -0.469   0.6632  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 786.3 on 4 degrees of freedom
    ## Multiple R-squared:  0.7388, Adjusted R-squared:  0.5428 
    ## F-statistic:  3.77 on 3 and 4 DF,  p-value: 0.1162

``` r
Anova(PD_PRic_geo_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                          Sum Sq Df F value  Pr(>F)  
    ## (Intercept)             1190980  1  1.9263 0.23747  
    ## max_depth                 25470  1  0.0412 0.84907  
    ## logArea                 3658461  1  5.9172 0.07178 .
    ## distance_to_ocean_min_m  136260  1  0.2204 0.66319  
    ## Residuals               2473103  4                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept)               max_depth                 logArea 
    ##               0.9498828               1.0000000               0.2871387 
    ## distance_to_ocean_min_m 
    ##               1.0000000

``` r
## Stratified lakes
PD_PRic_geo_S_lm <- lm(pd.obs ~ max_depth + logArea + distance_to_ocean_min_m, stree_sespd_env[stratified_lakes,])
summary(PD_PRic_geo_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ max_depth + logArea + distance_to_ocean_min_m, 
    ##     data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##     BCM     CLM     GLK     HLM     NLK     OTM     SLN     TLN 
    ##  105.32  -31.80 -227.61  376.51   25.78 -424.77  -22.82  199.40 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              123.408   2012.127   0.061    0.954
    ## max_depth                 -7.210     24.285  -0.297    0.781
    ## logArea                  147.399    248.917   0.592    0.586
    ## distance_to_ocean_min_m   -3.098      1.710  -1.812    0.144
    ## 
    ## Residual standard error: 326.7 on 4 degrees of freedom
    ## Multiple R-squared:  0.4592, Adjusted R-squared:  0.05355 
    ## F-statistic: 1.132 on 3 and 4 DF,  p-value: 0.4364

``` r
Anova(PD_PRic_geo_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                         Sum Sq Df F value Pr(>F)
    ## (Intercept)                402  1  0.0038 0.9540
    ## max_depth                 9410  1  0.0881 0.7813
    ## logArea                  37437  1  0.3507 0.5856
    ## distance_to_ocean_min_m 350555  1  3.2835 0.1442
    ## Residuals               427048  4

``` r
p_values <- summary(PD_PRic_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept)               max_depth                 logArea 
    ##               1.0000000               1.0000000               1.0000000 
    ## distance_to_ocean_min_m 
    ##               0.5768139

``` r
### PDisp ANOVA
## Surveyed sites
# Distance to the ocean
PD_PRic_lm_dist <- lm(Dispersion ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[surveyed_sites,])
plot(PD_PRic_lm_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   19

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-61.png)<!-- -->

``` r
PD_PRic_dist_am <- aov(Dispersion ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[surveyed_sites,])
# Max depth
PD_PRic_lm_mxd <- lm(Dispersion ~ Site_type * max_depth, data = stree_sespd_env[surveyed_sites,])
plot(PD_PRic_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-62.png)<!-- -->

``` r
PD_PRic_mxd_am <- aov(Dispersion ~ Site_type * max_depth, data = stree_sespd_env[surveyed_sites,])
# Log Area
PD_PRic_lm_lga <- lm(Dispersion ~ Site_type * logArea, data = stree_sespd_env[surveyed_sites,])
plot(PD_PRic_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-63.png)<!-- -->

``` r
PD_PRic_lga_am <- aov(Dispersion ~ Site_type * logArea, data = stree_sespd_env[surveyed_sites,])
# Anova outputs
summary(PD_PRic_dist_am)
```

    ##                                   Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                          2 0.000707 0.0003537   0.686  0.518
    ## distance_to_ocean_min_m            1 0.001259 0.0012587   2.441  0.138
    ## Site_type:distance_to_ocean_min_m  2 0.000046 0.0000228   0.044  0.957
    ## Residuals                         16 0.008252 0.0005157

``` r
summary(PD_PRic_mxd_am)
```

    ##                     Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type            2 0.000707 0.0003537   0.701  0.511
    ## max_depth            1 0.000201 0.0002007   0.398  0.537
    ## Site_type:max_depth  2 0.001284 0.0006421   1.273  0.307
    ## Residuals           16 0.008071 0.0005044

``` r
summary(PD_PRic_lga_am)
```

    ##                   Df   Sum Sq   Mean Sq F value Pr(>F)  
    ## Site_type          2 0.000707 0.0003537   0.851 0.4455  
    ## logArea            1 0.000164 0.0001636   0.393 0.5394  
    ## Site_type:logArea  2 0.002741 0.0013705   3.297 0.0633 .
    ## Residuals         16 0.006652 0.0004157                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_dist_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_mxd_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_lga_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [8] 1.0000000 0.5693403

``` r
## Mixed and stratified lakes
# Distance to the ocean
PD_PRic_MS_lm_dist <- lm(Dispersion ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-64.png)<!-- -->

``` r
PD_PRic_MS_dist_am <- aov(Dispersion ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[mixed_stratified_lakes,])
# Max depth
PD_PRic_MS_lm_mxd <- lm(Dispersion ~ Site_type * max_depth, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-65.png)<!-- -->

``` r
PD_PRic_MS_mxd_am <- aov(Dispersion ~ Site_type * max_depth, data = stree_sespd_env[mixed_stratified_lakes,])
# Log Area
PD_PRic_MS_lm_lga <- lm(Dispersion ~ Site_type * logArea, data = stree_sespd_env[mixed_stratified_lakes,])
plot(PD_PRic_MS_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-66.png)<!-- -->

``` r
PD_PRic_MS_lga_am <- aov(Dispersion ~ Site_type * logArea, data = stree_sespd_env[mixed_stratified_lakes,])
# Anova outputs
summary(PD_PRic_MS_dist_am)
```

    ##                                   Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                          1 0.000707 0.0007066   1.089  0.317
    ## distance_to_ocean_min_m            1 0.001265 0.0012654   1.951  0.188
    ## Site_type:distance_to_ocean_min_m  1 0.000038 0.0000383   0.059  0.812
    ## Residuals                         12 0.007785 0.0006488

``` r
summary(PD_PRic_MS_mxd_am)
```

    ##                     Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type            1 0.000707 0.0007066   1.078  0.320
    ## max_depth            1 0.000758 0.0007581   1.156  0.303
    ## Site_type:max_depth  1 0.000462 0.0004621   0.705  0.418
    ## Residuals           12 0.007869 0.0006557

``` r
summary(PD_PRic_MS_lga_am)
```

    ##                   Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type          1 0.000707 0.0007066   1.332  0.271
    ## logArea            1 0.001209 0.0012085   2.278  0.157
    ## Site_type:logArea  1 0.001514 0.0015144   2.855  0.117
    ## Residuals         12 0.006366 0.0005305

``` r
# p-values
temp_p_values <- summary(PD_PRic_MS_dist_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_MS_mxd_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_MS_lga_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1 1 1 1 1 1 1 1 1

``` r
## Ocean sites and mixed lakes
# Distance to the ocean
PD_PRic_OM_lm_dist <- lm(Dispersion ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_mixed_sites,])
plot(PD_PRic_OM_lm_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   6

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-67.png)<!-- -->

``` r
PD_PRic_OM_dist_am <- aov(Dispersion ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_mixed_sites,])
# Max depth
PD_PRic_OM_lm_mxd <- lm(Dispersion ~ Site_type * max_depth, data = stree_sespd_env[ocean_mixed_sites,])
plot(PD_PRic_OM_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-68.png)<!-- -->

``` r
PD_PRic_OM_mxd_am <- aov(Dispersion ~ Site_type * max_depth, data = stree_sespd_env[ocean_mixed_sites,])
# Log Area
PD_PRic_OM_lm_lga <- lm(Dispersion ~ Site_type * logArea, data = stree_sespd_env[ocean_mixed_sites,])
plot(PD_PRic_OM_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-69.png)<!-- -->

``` r
PD_PRic_OM_lga_am <- aov(Dispersion ~ Site_type * logArea, data = stree_sespd_env[ocean_mixed_sites,])
# Anova outputs
summary(PD_PRic_OM_dist_am)
```

    ##                                   Df    Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                          1 0.0001715 1.715e-04   1.741  0.216
    ## distance_to_ocean_min_m            1 0.0000623 6.235e-05   0.633  0.445
    ## Site_type:distance_to_ocean_min_m  1 0.0000015 1.510e-06   0.015  0.904
    ## Residuals                         10 0.0009846 9.846e-05

``` r
summary(PD_PRic_OM_mxd_am)
```

    ##                     Df    Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type            1 0.0001715 1.715e-04   2.189  0.170
    ## max_depth            1 0.0001523 1.523e-04   1.944  0.193
    ## Site_type:max_depth  1 0.0001129 1.129e-04   1.441  0.258
    ## Residuals           10 0.0007834 7.834e-05

``` r
summary(PD_PRic_OM_lga_am)
```

    ##                   Df    Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type          1 0.0001715 1.715e-04   1.994  0.188
    ## logArea            1 0.0000924 9.242e-05   1.075  0.324
    ## Site_type:logArea  1 0.0000964 9.639e-05   1.121  0.315
    ## Residuals         10 0.0008597 8.597e-05

``` r
# p-values
temp_p_values <- summary(PD_PRic_OM_dist_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_OM_mxd_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_OM_lga_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1 1 1 1 1 1 1 1 1

``` r
## Stratified lakes and ocean sites
# Distance to the ocean
PD_PRic_SO_lm_dist <- lm(Dispersion ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_stratified_sites,])
plot(PD_PRic_SO_lm_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   6

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-70.png)<!-- -->

``` r
PD_PRic_SO_dist_am <- aov(Dispersion ~ Site_type * distance_to_ocean_min_m, data = stree_sespd_env[ocean_stratified_sites,])
# Max depth
PD_PRic_SO_lm_mxd <- lm(Dispersion ~ Site_type * max_depth, data = stree_sespd_env[ocean_stratified_sites,])
plot(PD_PRic_SO_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-71.png)<!-- -->

``` r
PD_PRic_SO_mxd_am <- aov(Dispersion ~ Site_type * max_depth, data = stree_sespd_env[ocean_stratified_sites,])
# Log Area
PD_PRic_SO_lm_lga <- lm(Dispersion ~ Site_type * logArea, data = stree_sespd_env[ocean_stratified_sites,])
plot(PD_PRic_SO_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-72.png)<!-- -->

``` r
PD_PRic_SO_lga_am <- aov(Dispersion ~ Site_type * logArea, data = stree_sespd_env[ocean_stratified_sites,])
# Anova outputs
summary(PD_PRic_SO_dist_am)
```

    ##                                   Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type                          1 0.000133 0.0001326   0.171  0.688
    ## distance_to_ocean_min_m            1 0.001232 0.0012323   1.593  0.235
    ## Site_type:distance_to_ocean_min_m  1 0.000009 0.0000087   0.011  0.918
    ## Residuals                         10 0.007734 0.0007734

``` r
summary(PD_PRic_SO_mxd_am)
```

    ##                     Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Site_type            1 0.000133 0.0001326   0.177  0.683
    ## max_depth            1 0.000280 0.0002803   0.374  0.554
    ## Site_type:max_depth  1 0.001204 0.0012044   1.608  0.234
    ## Residuals           10 0.007490 0.0007490

``` r
summary(PD_PRic_SO_lga_am)
```

    ##                   Df   Sum Sq   Mean Sq F value Pr(>F)  
    ## Site_type          1 0.000133 0.0001326   0.218 0.6504  
    ## logArea            1 0.000177 0.0001769   0.291 0.6014  
    ## Site_type:logArea  1 0.002720 0.0027205   4.476 0.0605 .
    ## Residuals         10 0.006077 0.0006077                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_SO_dist_am)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_SO_mxd_am)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_SO_lga_am)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    ## [8] 1.0000000 0.5441027

``` r
## Ocean sites
PD_PRic_geo_O_lm <- lm(Dispersion ~ max_depth + logArea + distance_to_ocean_min_m, stree_sespd_env[ocean_sites,])
summary(PD_PRic_geo_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ max_depth + logArea + distance_to_ocean_min_m, 
    ##     data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ##        IBK        LCN        NCN        OOO        OOM        RCA 
    ##  2.187e-03  6.768e-03 -5.057e-03 -6.328e-03  2.430e-03 -2.458e-18 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              0.1969825  0.0200118   9.843   0.0102 *
    ## max_depth               -0.0004631  0.0003591  -1.290   0.3262  
    ## logArea                 -0.0018941  0.0017808  -1.064   0.3989  
    ## distance_to_ocean_min_m -0.0001342  0.0004753  -0.282   0.8042  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.007814 on 2 degrees of freedom
    ## Multiple R-squared:  0.7386, Adjusted R-squared:  0.3465 
    ## F-statistic: 1.884 on 3 and 2 DF,  p-value: 0.3653

``` r
Anova(PD_PRic_geo_O_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                            Sum Sq Df F value  Pr(>F)  
    ## (Intercept)             0.0059161  1 96.8908 0.01016 *
    ## max_depth               0.0001015  1  1.6629 0.32621  
    ## logArea                 0.0000691  1  1.1313 0.39893  
    ## distance_to_ocean_min_m 0.0000049  1  0.0797 0.80422  
    ## Residuals               0.0001221  2                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept)               max_depth                 logArea 
    ##              0.04065527              1.00000000              1.00000000 
    ## distance_to_ocean_min_m 
    ##              1.00000000

``` r
## Mixed lakes
PD_PRic_geo_M_lm <- lm(Dispersion ~ max_depth + logArea + distance_to_ocean_min_m, stree_sespd_env[mixed_lakes,])
summary(PD_PRic_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ max_depth + logArea + distance_to_ocean_min_m, 
    ##     data = stree_sespd_env[mixed_lakes, ])
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
    ## max_depth               -0.0008034  0.0007091  -1.133  0.32054   
    ## logArea                  0.0036215  0.0036400   0.995  0.37609   
    ## distance_to_ocean_min_m -0.0002204  0.0001623  -1.358  0.24611   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.009843 on 4 degrees of freedom
    ## Multiple R-squared:  0.3334, Adjusted R-squared:  -0.1665 
    ## F-statistic: 0.6669 on 3 and 4 DF,  p-value: 0.615

``` r
Anova(PD_PRic_geo_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                            Sum Sq Df F value  Pr(>F)   
    ## (Intercept)             0.0035139  1 36.2705 0.00383 **
    ## max_depth               0.0001244  1  1.2836 0.32054   
    ## logArea                 0.0000959  1  0.9898 0.37609   
    ## distance_to_ocean_min_m 0.0001786  1  1.8433 0.24611   
    ## Residuals               0.0003875  4                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept)               max_depth                 logArea 
    ##              0.01531847              1.00000000              1.00000000 
    ## distance_to_ocean_min_m 
    ##              0.98443978

``` r
## Stratified lakes
PD_PRic_geo_S_lm <- lm(Dispersion ~ max_depth + logArea + distance_to_ocean_min_m, stree_sespd_env[stratified_lakes,])
summary(PD_PRic_geo_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ max_depth + logArea + distance_to_ocean_min_m, 
    ##     data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##       BCM       CLM       GLK       HLM       NLK       OTM       SLN       TLN 
    ##  0.002563 -0.008611  0.028940 -0.025547  0.020337 -0.002319 -0.011485 -0.003878 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -0.1881649  0.1423307  -1.322   0.2567  
    ## max_depth               -0.0025469  0.0017178  -1.483   0.2123  
    ## logArea                  0.0441519  0.0176075   2.508   0.0662 .
    ## distance_to_ocean_min_m -0.0002779  0.0001209  -2.298   0.0831 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02311 on 4 degrees of freedom
    ## Multiple R-squared:  0.7488, Adjusted R-squared:  0.5605 
    ## F-statistic: 3.975 on 3 and 4 DF,  p-value: 0.1079

``` r
Anova(PD_PRic_geo_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                            Sum Sq Df F value  Pr(>F)  
    ## (Intercept)             0.0009336  1  1.7478 0.25669  
    ## max_depth               0.0011743  1  2.1983 0.21231  
    ## logArea                 0.0033590  1  6.2879 0.06623 .
    ## distance_to_ocean_min_m 0.0028211  1  5.2810 0.08312 .
    ## Residuals               0.0021368  4                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept)               max_depth                 logArea 
    ##               1.0000000               0.8492304               0.2649182 
    ## distance_to_ocean_min_m 
    ##               0.3324746

``` r
par(mfrow=c(1,1))
```

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
    ## 2     Stratified vs Ocean  1 1.1858559 8.357071 0.4105242   0.001      0.006
    ## 3 Stratified vs Reference  1 0.7729479 7.110647 0.5039207   0.097      0.582
    ## 4          Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.098      0.588
    ## 5      Mixed vs Reference  1 0.6592270 3.967203 0.3617334   0.114      0.684
    ## 6      Ocean vs Reference  1 0.6319662 3.354877 0.4015472   0.144      0.864
    ##   sig
    ## 1   *
    ## 2   *
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
    ## 2 Stratified vs Ocean  1 1.1858559 8.357071 0.4105242   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.107      0.321

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
    ## 1 Stratified vs Mixed  1 1.2940870 9.415922 0.4021162   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.2083230 9.439158 0.4618174   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.3347638 2.034034 0.1560556   0.026      0.078

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
    ## 1 Stratified vs Mixed  1 1.3499033 10.500850 0.4666868   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.2146330  9.192717 0.4789690   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.2573592  1.467099 0.1089395   0.109      0.327

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
# Reference
PD_beta_ref_NMDS <- metaMDS(PD_beta_ref_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(PD_beta_ref_dist$Btotal, try = 1000, parallel = 4, trymax =
    ## 500, : stress is (nearly) zero: you may have insufficient data

``` r
# Surveyed sites 
PD_beta_NMDS <- metaMDS(PD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
round(PD_beta_NMDS$stress, digits = 2)
PD_beta_rep_NMDS <- metaMDS(PD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
PD_beta_ric_NMDS <- metaMDS(PD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
PD_beta_MS_NMDS <- metaMDS(PD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
PD_beta_OM_NMDS <- metaMDS(PD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
PD_beta_SO_NMDS <- metaMDS(PD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes
PD_beta_S_NMDS <- metaMDS(PD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)


### Environmental
# Surveyed sites 
PD_beta_env_NMDS <- metaMDS(PD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
PD_beta_env_MS_NMDS <- metaMDS(PD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
PD_beta_env_OM_NMDS <- metaMDS(PD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
PD_beta_env_SO_NMDS <- metaMDS(PD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(PD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
PD_beta_env_M_NMDS <- metaMDS(PD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(PD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
PD_beta_geo_NMDS <- metaMDS(PD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
PD_beta_geo_MS_NMDS <- metaMDS(PD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
PD_beta_geo_OM_NMDS <- metaMDS(PD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
PD_beta_geo_SO_NMDS <- metaMDS(PD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed lakes
PD_beta_geo_M_NMDS <- metaMDS(PD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

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
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31832 0.699  0.012 *
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
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31832 0.699  0.012 *
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
    ## salinity_median     0.94798  0.31832 0.6990  0.008 **
    ## oxygen_median       0.57185 -0.82036 0.6836  0.110   
    ## temperature_median -0.75815 -0.65208 0.0724  0.658   
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
    ## salinity_median     0.94798  0.31832 0.6990  0.024 *
    ## oxygen_median       0.57185 -0.82036 0.6836  0.330  
    ## temperature_median -0.75815 -0.65208 0.0724  1.000  
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
    ## salinity_median     0.90473  0.42600 0.6840  0.026 *
    ## oxygen_median       0.44477 -0.89564 0.5836  0.200  
    ## temperature_median -0.92696  0.37515 0.0954  0.737  
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
    ## salinity_median     0.90473  0.42600 0.6840  0.078 .
    ## oxygen_median       0.44477 -0.89564 0.5836  0.600  
    ## temperature_median -0.92696  0.37515 0.0954  1.000  
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
    ## salinity_median    -0.79785 -0.60286 0.7269  0.014 *
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.015 *
    ## temperature_median -0.52519 -0.85098 0.0223  0.913  
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
    ## salinity_median    -0.79785 -0.60286 0.7269  0.042 *
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.045 *
    ## temperature_median -0.52519 -0.85098 0.0223  1.000  
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
    ## salinity_median    -0.0052518  0.9999900 0.4834  0.863  
    ## oxygen_median      -0.0008316  1.0000000 0.7903  0.065 .
    ## temperature_median  0.0007131 -1.0000000 0.0342  0.844  
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
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## salinity_median    -0.0052518  0.9999900 0.4834  1.000
    ## oxygen_median      -0.0008316  1.0000000 0.7903  0.195
    ## temperature_median  0.0007131 -1.0000000 0.0342  1.000
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
    ## salinity_median     0.00025711 -1.00000000 0.4119  0.261  
    ## oxygen_median       0.00012696 -1.00000000 0.7912  0.022 *
    ## temperature_median -0.00017570  1.00000000 0.0238  0.925  
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
    ## salinity_median     0.00025711 -1.00000000 0.4119  0.783  
    ## oxygen_median       0.00012696 -1.00000000 0.7912  0.066 .
    ## temperature_median -0.00017570  1.00000000 0.0238  1.000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ## temperature_median       -0.39650 -0.91803 0.0381  0.886  
    ## salinity_median          -0.96164  0.27432 0.5830  0.116  
    ## oxygen_median             0.83740  0.54659 0.7450  0.030 *
    ## distance_to_ocean_mean_m -0.09986  0.99500 0.0319  0.930  
    ## max_depth                -0.87451  0.48501 0.1022  0.777  
    ## logArea                  -0.40927  0.91241 0.2148  0.579  
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
    ## temperature_median       -0.39650 -0.91803 0.0381  1.000
    ## salinity_median          -0.96164  0.27432 0.5830  0.696
    ## oxygen_median             0.83740  0.54659 0.7450  0.180
    ## distance_to_ocean_mean_m -0.09986  0.99500 0.0319  1.000
    ## max_depth                -0.87451  0.48501 0.1022  1.000
    ## logArea                  -0.40927  0.91241 0.2148  1.000
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
    ## distance_to_ocean_min_m -0.82790  0.56088 0.5558  0.053 .
    ## max_depth               -0.20870 -0.97798 0.0858  0.904  
    ## logArea                  0.25301 -0.96746 0.2012  0.317  
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
    ## distance_to_ocean_min_m -0.82790  0.56088 0.5558  0.159
    ## max_depth               -0.20870 -0.97798 0.0858  1.000
    ## logArea                  0.25301 -0.96746 0.2012  0.951
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
    ## distance_to_ocean_min_m -0.82790  0.56088 0.5558  0.074 .
    ## max_depth               -0.20870 -0.97798 0.0858  0.921  
    ## logArea                  0.25301 -0.96746 0.2012  0.293  
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
    ## distance_to_ocean_min_m -0.82790  0.56088 0.5558  0.222
    ## max_depth               -0.20870 -0.97798 0.0858  1.000
    ## logArea                  0.25301 -0.96746 0.2012  0.879
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
    ## distance_to_ocean_min_m -0.76294  0.64648 0.4672  0.128
    ## max_depth               -0.34723 -0.93778 0.0777  0.968
    ## logArea                  0.55358  0.83279 0.0074  0.969
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
    ## distance_to_ocean_min_m -0.76294  0.64648 0.4672  0.384
    ## max_depth               -0.34723 -0.93778 0.0777  1.000
    ## logArea                  0.55358  0.83279 0.0074  1.000
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
    ## max_depth               -0.94119 -0.33789 0.6682  0.007 **
    ## logArea                 -0.86397  0.50354 0.2339  0.443   
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
    ## max_depth               -0.94119 -0.33789 0.6682  0.021 *
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
    ## distance_to_ocean_min_m  0.99474  0.10244 0.6338  0.228
    ## max_depth                0.73004  0.68340 0.0374  0.983
    ## logArea                 -0.39131  0.92026 0.3079  0.290
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
    ## distance_to_ocean_min_m  0.99474  0.10244 0.6338  0.684
    ## max_depth                0.73004  0.68340 0.0374  1.000
    ## logArea                 -0.39131  0.92026 0.3079  0.870
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
    ## distance_to_ocean_min_m -0.00025919 -1.00000000 0.5355  0.144  
    ## max_depth                0.00094670  1.00000000 0.7411  0.045 *
    ## logArea                  0.00066133 -1.00000000 0.7230  0.030 *
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
    ## distance_to_ocean_min_m -0.00025919 -1.00000000 0.5355  0.432  
    ## max_depth                0.00094670  1.00000000 0.7411  0.135  
    ## logArea                  0.00066133 -1.00000000 0.7230  0.090 .
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
    ##       Significance: 0.217 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.214 0.245 0.268 0.294 
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
    ##       Significance: 0.008 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.437 0.469 0.489 0.510 
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
    ##       Significance: 0.262 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.525 0.565 0.600 0.648 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_mant_pv <- rbind(PD_beta_env_mant_t$signif, PD_beta_env_mant_s$signif, PD_beta_env_mant_o$signif)
PD_beta_env_mant_pv <- PD_beta_env_mant_pv[,1]
(PD_beta_env_mant_pv <- p.adjust(PD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.651 0.024 0.786

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
    ##       Significance: 0.24 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.209 0.235 0.280 0.299 
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
    ## 0.374 0.408 0.433 0.460 
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
    ##       Significance: 0.287 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.413 0.447 0.480 0.527 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_MS_mant_pv <- rbind(PD_beta_env_MS_mant_t$signif, PD_beta_env_MS_mant_s$signif, PD_beta_env_MS_mant_o$signif)
PD_beta_env_MS_mant_pv <- PD_beta_env_MS_mant_pv[,1]
(PD_beta_env_MS_mant_pv <- p.adjust(PD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.720 0.030 0.861

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
    ##       Significance: 0.503 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.129 0.182 0.223 0.356 
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
    ##       Significance: 0.025 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.323 0.375 0.453 0.503 
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
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.285 0.390 0.483 0.563 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_OM_mant_pv <- rbind(PD_beta_env_OM_mant_t$signif, PD_beta_env_OM_mant_s$signif, PD_beta_env_OM_mant_o$signif)
PD_beta_env_OM_mant_pv <- PD_beta_env_OM_mant_pv[,1]
(PD_beta_env_OM_mant_pv <- p.adjust(PD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.075 0.036

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
    ##       Significance: 0.467 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0455 0.0897 0.1245 0.1577 
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
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.381 0.422 0.443 0.481 
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
    ##       Significance: 0.217 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.685 0.712 0.730 0.767 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_SO_mant_pv <- rbind(PD_beta_env_SO_mant_t$signif, PD_beta_env_SO_mant_s$signif, PD_beta_env_SO_mant_o$signif)
PD_beta_env_SO_mant_pv <- PD_beta_env_SO_mant_pv[,1]
(PD_beta_env_SO_mant_pv <- p.adjust(PD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.036 0.651

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
    ##       Significance: 0.745 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.258 0.355 0.452 0.676 
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
    ##       Significance: 0.255 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.222 0.284 0.348 0.393 
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
    ##       Significance: 0.086 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.251 0.361 0.435 0.510 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_M_mant_pv <- rbind(PD_beta_env_M_mant_t$signif, PD_beta_env_M_mant_s$signif, PD_beta_env_M_mant_o$signif)
PD_beta_env_M_mant_pv <- PD_beta_env_M_mant_pv[,1]
(PD_beta_env_M_mant_pv <- p.adjust(PD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.765 0.258

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
    ##       Significance: 0.284 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.351 0.434 0.516 0.583 
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
    ##       Significance: 0.586 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.461 0.492 0.513 0.543 
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
    ##       Significance: 0.787 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.166 0.190 0.215 0.245 
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
    ##       Significance: 0.725 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.113 0.134 0.147 0.157 
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
    ##       Significance: 0.655 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.316 0.354 0.383 0.410 
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
    ##       Significance: 0.933 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.162 0.204 0.233 0.275 
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
    ##       Significance: 0.894 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.133 0.154 0.172 0.197 
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
    ##       Significance: 0.396 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.166 0.201 0.216 0.264 
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
    ##       Significance: 0.014 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.159 0.216 0.264 0.330 
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
    ##       Significance: 0.068 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.171 0.206 0.226 0.268 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_OM_mant_pv <- rbind( PD_beta_geo_OM_mant_dmean$signif, PD_beta_geo_OM_mant_md$signif, PD_beta_geo_OM_mant_la$signif)
PD_beta_geo_OM_mant_pv <- PD_beta_geo_OM_mant_pv[,1]
(PD_beta_geo_OM_mant_pv <- p.adjust(PD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.042 0.204

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
    ##       Significance: 0.537 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.620 0.640 0.668 0.686 
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
    ##       Significance: 0.802 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.115 0.148 0.179 0.203 
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
    ##       Significance: 0.218 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.372 0.400 0.418 0.437 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_SO_mant_pv <- rbind( PD_beta_geo_SO_mant_dmean$signif, PD_beta_geo_SO_mant_md$signif, PD_beta_geo_SO_mant_la$signif)
PD_beta_geo_SO_mant_pv <- PD_beta_geo_SO_mant_pv[,1]
(PD_beta_geo_SO_mant_pv <- p.adjust(PD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.654

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
    ##       Significance: 0.695 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.223 0.300 0.392 0.699 
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
    ## 0.280 0.403 0.541 0.588 
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
    ##       Significance: 0.022 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.203 0.289 0.383 0.436 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_M_mant_pv <- rbind( PD_beta_geo_M_mant_dmean$signif, PD_beta_geo_M_mant_md$signif, PD_beta_geo_M_mant_la$signif)
PD_beta_geo_M_mant_pv <- PD_beta_geo_M_mant_pv[,1]
(PD_beta_geo_M_mant_pv <- p.adjust(PD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.033 0.066

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
    ##  [1] MASS_7.3-65          pairwiseAdonis_0.4.1 cluster_2.1.8       
    ##  [4] BAT_2.9.6            caret_7.0-1          ggrepel_0.9.6       
    ##  [7] viridis_0.6.5        viridisLite_0.4.2    ggplot2_3.5.1       
    ## [10] picante_1.8.2        nlme_3.1-167         vegan_2.6-10        
    ## [13] lattice_0.22-6       permute_0.9-7        car_3.1-3           
    ## [16] carData_3.0-5        tidyr_1.3.1          phytools_2.4-4      
    ## [19] maps_3.4.2.1         ape_5.8-1            reshape2_1.4.4      
    ## [22] stringr_1.5.1        dplyr_1.1.4          knitr_1.49          
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
    ##  [52] quantreg_6.00           lava_1.8.1              scatterplot3d_0.3-44   
    ##  [55] ModelMetrics_1.2.2.2    tools_4.4.3             foreign_0.8-88         
    ##  [58] future.apply_1.11.3     nnet_7.3-20             glue_1.8.0             
    ##  [61] quadprog_1.5-8          grid_4.4.3              checkmate_2.3.2        
    ##  [64] generics_0.1.3          recipes_1.1.1           gtable_0.3.6           
    ##  [67] class_7.3-23            data.table_1.17.0       hms_1.1.3              
    ##  [70] foreach_1.5.2           pillar_1.10.1           splines_4.4.3          
    ##  [73] survival_3.8-3          SparseM_1.84-2          ks_1.14.3              
    ##  [76] tidyselect_1.2.1        rms_7.0-0               gridExtra_2.3          
    ##  [79] stats4_4.4.3            xfun_0.51               expm_1.0-0             
    ##  [82] hardhat_1.4.1           timeDate_4041.110       proto_1.0.0            
    ##  [85] stringi_1.8.4           yaml_2.3.10             evaluate_1.0.3         
    ##  [88] codetools_0.2-20        tibble_3.2.1            cli_3.6.4              
    ##  [91] rpart_4.1.24            nls2_0.3-4              geometry_0.5.2         
    ##  [94] systemfonts_1.2.1       munsell_0.5.1           Rcpp_1.0.14            
    ##  [97] globals_0.16.3          coda_0.19-4.1           fastcluster_1.2.6      
    ## [100] MatrixModels_0.5-3      gower_1.0.2             prettyunits_1.2.0      
    ## [103] mclust_6.1.1            listenv_0.9.1           phangorn_2.12.1        
    ## [106] mvtnorm_1.3-3           ipred_0.9-15            scales_1.3.0           
    ## [109] prodlim_2024.06.25      e1071_1.7-16            purrr_1.0.4            
    ## [112] crayon_1.5.3            combinat_0.0-8          rlang_1.1.5            
    ## [115] multcomp_1.4-28         fastmatch_1.1-6         mnormt_2.1.1           
    ## [118] hypervolume_3.1.5

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
