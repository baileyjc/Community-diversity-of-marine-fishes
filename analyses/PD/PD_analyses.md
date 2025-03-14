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
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                     Sum Sq Df F value  Pr(>F)  
    ## (Intercept)        3454556  1  8.4830 0.01136 *
    ## temperature_median  902998  1  2.2174 0.15864  
    ## salinity_median      60080  1  0.1475 0.70668  
    ## oxygen_median       325406  1  0.7991 0.38648  
    ## pH_median          3432036  1  8.4277 0.01157 *
    ## Residuals          5701242 14                  
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
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs
    ##                     Sum Sq Df F value   Pr(>F)   
    ## (Intercept)         430749  1  0.7074 0.413504   
    ## temperature_median   40929  1  0.0672 0.798954   
    ## salinity_median    7052231  1 11.5822 0.003931 **
    ## oxygen_median      5692040  1  9.3483 0.007980 **
    ## Residuals          9133278 15                    
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

mod <- aov(temperature_median ~ Site_type, stree_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2   2.94   1.470   1.092  0.359
    ## Residuals   16  21.54   1.346               
    ## 3 observations deleted due to missingness

``` r
mod <- aov(pd.obs ~ temperature_median + Site_type, stree_sespd_env)
car::vif(mod)
```

    ##                        GVIF Df GVIF^(1/(2*Df))
    ## temperature_median 1.136543  1        1.066088
    ## Site_type          1.136543  2        1.032515

``` r
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
mod <- aov(pd.obs ~ salinity_median + Site_type, stree_sespd_env)
car::vif(mod)
```

    ##                     GVIF Df GVIF^(1/(2*Df))
    ## salinity_median 2.701125  1        1.643510
    ## Site_type       2.701125  2        1.281995

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
mod <- aov(pd.obs ~ oxygen_median + Site_type, stree_sespd_env)
car::vif(mod)
```

    ##                   GVIF Df GVIF^(1/(2*Df))
    ## oxygen_median 3.374424  1        1.836961
    ## Site_type     3.374424  2        1.355345

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
mod <- aov(pd.obs ~ distance_to_ocean_min_m + Site_type, stree_sespd_env)
car::vif(mod)
```

    ##                             GVIF Df GVIF^(1/(2*Df))
    ## distance_to_ocean_min_m 2.735421  1        1.653911
    ## Site_type               2.735421  2        1.286045

``` r
mod <- aov(max_depth ~ Site_type, stree_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2    535   267.5   2.235  0.134
    ## Residuals   19   2274   119.7

``` r
mod <- aov(pd.obs ~ max_depth + Site_type, stree_sespd_env)
car::vif(mod)
```

    ##               GVIF Df GVIF^(1/(2*Df))
    ## max_depth 1.235315  1        1.111447
    ## Site_type 1.235315  2        1.054252

``` r
mod <- aov(logArea ~ Site_type, stree_sespd_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2  14.37   7.184   2.341  0.123
    ## Residuals   19  58.32   3.069

``` r
mod <- aov(pd.obs ~ logArea + Site_type, stree_sespd_env)
car::vif(mod)
```

    ##               GVIF Df GVIF^(1/(2*Df))
    ## logArea   1.246371  1        1.116410
    ## Site_type 1.246371  2        1.056603

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
anova_PRicz_result <- aov(pd.obs.z ~ Site_type, data = stree_sespd_env[surveyed_sites,])
summary(anova_PRicz_result)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2  1.838   0.919   0.888  0.428
    ## Residuals   19 19.668   1.035

``` r
# Perform Tukey's HSD test
tukey_PRicz_result <- TukeyHSD(anova_PRicz_result)
print(tukey_PRicz_result)
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
N <- length(anova_PRicz_result$residuals)
# Categories
k <- length(anova_PRicz_result$coefficients)
# Sum of squares
SS <- sum(resid(anova_PRicz_result)^2) # from anova, the residual SS
# Mean squares
MS <- SS/(N-k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

# Get differences
Site_type <- tukey_PRicz_result$Site_type

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
### PRic ANOVA
## Surveyed sites
# Temperature
# Interaction
PD_PRic_am_temp <- aov(pd.obs ~ temperature_median * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_am_temp)
```

    ##                              Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## temperature_median            1  3431415  3431415   4.868 0.045963 *  
    ## Site_type                     2 20074783 10037392  14.239 0.000531 ***
    ## temperature_median:Site_type  2  1377325   688662   0.977 0.402473    
    ## Residuals                    13  9164076   704929                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_lm_temp <- lm(pd.obs ~ temperature_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_lm_temp)
    ## W = 0.96976, p-value = 0.7716

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_lm_temp) ~ stree_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  2  7.5978 0.004789 **
    ##       16                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN  2.81112           0.013873      0.26359

``` r
# Plot residuals
plot(PD_PRic_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-1.png)<!-- -->

``` r
# Test relationship
PD_PRic_am_temp <- aov(pd.obs ~ temperature_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])

# Salinity
# Interaction
PD_PRic_am_sal <- aov(pd.obs ~ salinity_median * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_am_sal)
```

    ##                           Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## salinity_median            1 19125479 19125479  30.693 9.54e-05 ***
    ## Site_type                  2  5089509  2544754   4.084    0.042 *  
    ## salinity_median:Site_type  2  1732017   866008   1.390    0.284    
    ## Residuals                 13  8100594   623123                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_lm_sal <- lm(pd.obs ~ salinity_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_lm_sal)
    ## W = 0.95265, p-value = 0.4379

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_lm_sal) ~ stree_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  2  10.107 0.001452 **
    ##       16                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.973688           0.010063      0.19119

``` r
# Plot residuals
plot(PD_PRic_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-2.png)<!-- -->

``` r
# Test relationship
PD_PRic_am_sal <- aov(pd.obs ~ salinity_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])

# Oxygen
# Interaction
PD_PRic_am_oxy <- aov(pd.obs ~ oxygen_median * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRic_am_oxy)
```

    ##                         Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## oxygen_median            1 16717491 16717491  28.347 0.000138 ***
    ## Site_type                2  6867037  3433519   5.822 0.015649 *  
    ## oxygen_median:Site_type  2  2796399  1398199   2.371 0.132483    
    ## Residuals               13  7666672   589744                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_lm_oxy <- lm(pd.obs ~ oxygen_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_lm_oxy)
    ## W = 0.96468, p-value = 0.6672

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_lm_oxy) ~ stree_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  2  6.2702 0.009756 **
    ##       16                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.777901            0.01481       0.2814

``` r
# Plot residuals
plot(PD_PRic_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-3.png)<!-- -->

``` r
# Test relationship
PD_PRic_am_oxy <- aov(pd.obs ~ oxygen_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])

# Anova outputs
summary(PD_PRic_am_temp)
```

    ##                    Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## temperature_median  1  3431415  3431415   4.883 0.043090 *  
    ## Site_type           2 20074783 10037392  14.283 0.000337 ***
    ## Residuals          15 10541401   702760                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_am_sal)
```

    ##                 Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## salinity_median  1 19125479 19125479  29.177 7.35e-05 ***
    ## Site_type        2  5089509  2544754   3.882   0.0438 *  
    ## Residuals       15  9832611   655507                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_am_oxy)
```

    ##               Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## oxygen_median  1 16717491 16717491  23.966 0.000194 ***
    ## Site_type      2  6867037  3433519   4.922 0.022721 *  
    ## Residuals     15 10463070   697538                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.2585401581 0.0020195553           NA 0.0004408375 0.2626870911
    ## [6]           NA 0.0011641337 0.1363281099           NA

``` r
## Mixed and stratified lakes
# Temperature
# Interaction
PD_PRic_MS_am_temp <- aov(pd.obs ~ temperature_median * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_MS_am_temp)
```

    ##                              Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## temperature_median            1  4037355  4037355   5.289 0.040215 *  
    ## Site_type                     1 16214027 16214027  21.240 0.000602 ***
    ## temperature_median:Site_type  1   900909   900909   1.180 0.298659    
    ## Residuals                    12  9160361   763363                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_MS_lm_temp <- lm(pd.obs ~ temperature_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_MS_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_MS_lm_temp)
    ## W = 0.97644, p-value = 0.9291

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_MS_lm_temp) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  1  12.252 0.003533 **
    ##       14                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_MS_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN  2.61622           0.022542      0.36066

``` r
# Plot residuals
plot(PD_PRic_MS_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-4.png)<!-- -->

``` r
# Test relationship
PD_PRic_MS_am_temp <- aov(pd.obs ~ temperature_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Salinity
# Interaction
PD_PRic_MS_am_sal <- aov(pd.obs ~ salinity_median * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_MS_am_sal)
```

    ##                           Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## salinity_median            1 15851356 15851356  24.310 0.000348 ***
    ## Site_type                  1  4976995  4976995   7.633 0.017191 *  
    ## salinity_median:Site_type  1  1659724  1659724   2.545 0.136600    
    ## Residuals                 12  7824577   652048                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_MS_lm_sal <- lm(pd.obs ~ salinity_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_MS_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_MS_lm_sal)
    ## W = 0.95244, p-value = 0.5294

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_MS_lm_sal) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)    
    ## group  1  17.701 0.000878 ***
    ##       14                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_MS_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.836453           0.014995      0.23992

``` r
# Plot residuals
plot(PD_PRic_MS_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-5.png)<!-- -->

``` r
# Test relationship
PD_PRic_MS_am_sal <- aov(pd.obs ~ salinity_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Oxygen
# Interaction
PD_PRic_MS_am_oxy <- aov(pd.obs ~ oxygen_median * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_MS_am_oxy)
```

    ##                         Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## oxygen_median            1 13219742 13219742   20.86 0.000647 ***
    ## Site_type                1  6919170  6919170   10.91 0.006296 ** 
    ## oxygen_median:Site_type  1  2566988  2566988    4.05 0.067183 .  
    ## Residuals               12  7606751   633896                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_MS_lm_oxy <- lm(pd.obs ~ oxygen_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_MS_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_MS_lm_oxy)
    ## W = 0.96513, p-value = 0.7549

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_MS_lm_oxy) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  1   10.52 0.005889 **
    ##       14                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_MS_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.656788           0.020915      0.33464

``` r
# Plot residuals
plot(PD_PRic_MS_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-6.png)<!-- -->

``` r
# Test relationship
PD_PRic_MS_am_oxy <- aov(pd.obs ~ oxygen_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Anova outputs
summary(PD_PRic_MS_am_temp)
```

    ##                    Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## temperature_median  1  4037355  4037355   5.217 0.039824 *  
    ## Site_type           1 16214027 16214027  20.950 0.000519 ***
    ## Residuals          13 10061270   773944                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_MS_am_sal)
```

    ##                 Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## salinity_median  1 15851356 15851356  21.727 0.000445 ***
    ## Site_type        1  4976995  4976995   6.822 0.021515 *  
    ## Residuals       13  9484300   729562                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_MS_am_oxy)
```

    ##               Df   Sum Sq  Mean Sq F value  Pr(>F)   
    ## oxygen_median  1 13219742 13219742  16.892 0.00123 **
    ## Site_type      1  6919170  6919170   8.841 0.01078 * 
    ## Residuals     13 10173739   782595                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_MS_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_MS_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_MS_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.238942498 0.003113216          NA 0.002672124 0.129092974          NA
    ## [7] 0.007379200 0.064652705          NA

``` r
## Ocean sites and mixed lakes
# Temperature
# Interaction
PD_PRic_OM_am_temp <- aov(pd.obs ~ temperature_median * Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PRic_OM_am_temp)
```

    ##                              Df  Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median            1  204662  204662   0.171  0.692
    ## Site_type                     1   28028   28028   0.023  0.883
    ## temperature_median:Site_type  1 1216000 1216000   1.016  0.347
    ## Residuals                     7 8379005 1197001

``` r
# Linear model
PD_PRic_OM_lm_temp <- lm(pd.obs ~ temperature_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_OM_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_OM_lm_temp)
    ## W = 0.97522, p-value = 0.9339

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_OM_lm_temp) ~ stree_sespd_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.4294 0.1535
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRic_OM_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.218473           0.062017      0.68219

``` r
# Plot residuals
plot(PD_PRic_OM_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-7.png)<!-- -->

``` r
# Test relationship
PD_PRic_OM_am_temp <- aov(pd.obs ~ temperature_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])

# Salinity
# Interaction
PD_PRic_OM_am_sal <- aov(pd.obs ~ salinity_median * Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PRic_OM_am_sal)
```

    ##                           Df  Sum Sq Mean Sq F value Pr(>F)
    ## salinity_median            1 1365849 1365849   1.221  0.306
    ## Site_type                  1  631766  631766   0.565  0.477
    ## salinity_median:Site_type  1    2836    2836   0.003  0.961
    ## Residuals                  7 7827244 1118178

``` r
# Linear model
PD_PRic_OM_lm_sal <- lm(pd.obs ~ salinity_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_OM_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_OM_lm_sal)
    ## W = 0.91981, p-value = 0.3171

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_OM_lm_sal) ~ stree_sespd_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1   1.007 0.3418
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRic_OM_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN  2.27592           0.056978      0.62676

``` r
# Plot residuals
plot(PD_PRic_OM_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-8.png)<!-- -->

``` r
# Test relationship
PD_PRic_OM_am_sal <- aov(pd.obs ~ salinity_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])

# Oxygen
# Interaction
PD_PRic_OM_am_oxy <- aov(pd.obs ~ oxygen_median * Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PRic_OM_am_oxy)
```

    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## oxygen_median            1 1569097 1569097   1.521  0.257
    ## Site_type                1 1036520 1036520   1.005  0.350
    ## oxygen_median:Site_type  1    1412    1412   0.001  0.972
    ## Residuals                7 7220665 1031524

``` r
# Linear model
PD_PRic_OM_lm_oxy <- lm(pd.obs ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_OM_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_OM_lm_oxy)
    ## W = 0.97507, p-value = 0.9326

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_OM_lm_oxy) ~ stree_sespd_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  4.1076 0.07332 .
    ##        9                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_OM_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.020842           0.083027       0.9133

``` r
# Plot residuals
plot(PD_PRic_OM_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-9.png)<!-- -->

``` r
# Test relationship
PD_PRic_OM_am_oxy <- aov(pd.obs ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])

# Anova outputs
summary(PD_PRic_OM_am_temp)
```

    ##                    Df  Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1  204662  204662   0.171  0.690
    ## Site_type           1   28028   28028   0.023  0.882
    ## Residuals           8 9595005 1199376

``` r
summary(PD_PRic_OM_am_sal)
```

    ##                 Df  Sum Sq Mean Sq F value Pr(>F)
    ## salinity_median  1 1365849 1365849   1.395  0.271
    ## Site_type        1  631766  631766   0.645  0.445
    ## Residuals        8 7830080  978760

``` r
summary(PD_PRic_OM_am_oxy)
```

    ##               Df  Sum Sq Mean Sq F value Pr(>F)
    ## oxygen_median  1 1569097 1569097   1.738  0.224
    ## Site_type      1 1036520 1036520   1.148  0.315
    ## Residuals      8 7222077  902760

``` r
# p-values
temp_p_values <- summary(PD_PRic_OM_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_OM_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_OM_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1]  1  1 NA  1  1 NA  1  1 NA

``` r
## Stratified lakes and ocean sites
# Temperature
# Interaction
PD_PRic_SO_am_temp <- aov(pd.obs ~ temperature_median * Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PRic_SO_am_temp)
```

    ##                              Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## temperature_median            1   229718   229718   2.039    0.196    
    ## Site_type                     1 11080323 11080323  98.331 2.26e-05 ***
    ## temperature_median:Site_type  1   348618   348618   3.094    0.122    
    ## Residuals                     7   788785   112684                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_SO_lm_temp <- lm(pd.obs ~ temperature_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_SO_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_SO_lm_temp)
    ## W = 0.84999, p-value = 0.04267

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_SO_lm_temp) ~ stree_sespd_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0376 0.8505
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRic_SO_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLM 2.450886            0.04405      0.48455

``` r
# Plot residuals
plot(PD_PRic_SO_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-10.png)<!-- -->

``` r
# Test relationship
PD_PRic_SO_am_temp <- aov(pd.obs ~ temperature_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])

# Salinity
# Interaction
PD_PRic_SO_am_sal <- aov(pd.obs ~ salinity_median * Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PRic_SO_am_sal)
```

    ##                           Df  Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median            1 8152385 8152385 103.877 1.89e-05 ***
    ## Site_type                  1 3671719 3671719  46.785 0.000244 ***
    ## salinity_median:Site_type  1   73973   73973   0.943 0.363962    
    ## Residuals                  7  549367   78481                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_SO_lm_sal <- lm(pd.obs ~ salinity_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_SO_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_SO_lm_sal)
    ## W = 0.96881, p-value = 0.8742

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_SO_lm_sal) ~ stree_sespd_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.9872 0.3464
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRic_SO_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## IBK 2.697226           0.030762      0.33839

``` r
# Plot residuals
plot(PD_PRic_SO_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-11.png)<!-- -->

``` r
# Test relationship
PD_PRic_SO_am_sal <- aov(pd.obs ~ salinity_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])

# Oxygen
# Interaction
PD_PRic_SO_am_oxy <- aov(pd.obs ~ oxygen_median * Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PRic_SO_am_oxy)
```

    ##                         Df  Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median            1 6737095 6737095  93.214 2.70e-05 ***
    ## Site_type                1 4761519 4761519  65.880 8.31e-05 ***
    ## oxygen_median:Site_type  1  442904  442904   6.128   0.0425 *  
    ## Residuals                7  505927   72275                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_SO_lm_oxy <- lm(pd.obs ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_SO_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_SO_lm_oxy)
    ## W = 0.91805, p-value = 0.3028

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_SO_lm_oxy) ~ stree_sespd_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.6897 0.4277
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRic_SO_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## IBK 2.620448           0.034387      0.37825

``` r
# Plot residuals
plot(PD_PRic_SO_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-12.png)<!-- -->

``` r
# Test relationship
PD_PRic_SO_am_oxy <- aov(pd.obs ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])

# Anova outputs
summary(PD_PRic_SO_am_temp)
```

    ##                    Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## temperature_median  1   229718   229718   1.616    0.239    
    ## Site_type           1 11080323 11080323  77.934 2.14e-05 ***
    ## Residuals           8  1137403   142175                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_SO_am_sal)
```

    ##                 Df  Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median  1 8152385 8152385  104.63 7.17e-06 ***
    ## Site_type        1 3671719 3671719   47.12 0.000129 ***
    ## Residuals        8  623341   77918                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_SO_am_oxy)
```

    ##               Df  Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median  1 6737095 6737095   56.80 6.69e-05 ***
    ## Site_type      1 4761519 4761519   40.15 0.000224 ***
    ## Residuals      8  948830  118604                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRic_SO_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRic_SO_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRic_SO_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.000000e+00 1.281018e-04           NA 4.300386e-05 7.743699e-04
    ## [6]           NA 4.014569e-04 1.343510e-03           NA

``` r
## Ocean sites
PD_PRic_env_O_lm <- lm(pd.obs ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_sites_env,])
summary(PD_PRic_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[ocean_sites_env, ])
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
#### Geographical
### PRic ANOVA
## Surveyed sites
# Distance
# Interaction
PD_PRic_am_dist <- aov(pd.obs ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PRic_am_dist)
```

    ##                                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1 14204403 14204403  18.094 0.000607 ***
    ## Site_type                          2  9840464  4920232   6.267 0.009772 ** 
    ## distance_to_ocean_min_m:Site_type  2    35606    17803   0.023 0.977609    
    ## Residuals                         16 12560885   785055                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_lm_dist <- lm(pd.obs ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_lm_dist)
    ## W = 0.97769, p-value = 0.8763

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_lm_dist) ~ stree_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)   
    ## group  2  7.2171 0.00466 **
    ##       19                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.972697          0.0085371      0.18782

``` r
# Plot residuals
plot(PD_PRic_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-13.png)<!-- -->

``` r
# Test relationship
PD_PRic_am_dist <- aov(pd.obs ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[surveyed_sites,])

# Max depth
# Interaction
PD_PRic_am_mxd <- aov(pd.obs ~ max_depth * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PRic_am_mxd)
```

    ##                     Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## max_depth            1   897140   897140   1.496    0.239    
    ## Site_type            2 24398371 12199185  20.339 4.03e-05 ***
    ## max_depth:Site_type  2  1749289   874645   1.458    0.262    
    ## Residuals           16  9596558   599785                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_lm_mxd <- lm(pd.obs ~ max_depth + Site_type, data = stree_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_lm_mxd)
    ## W = 0.96712, p-value = 0.6445

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_lm_mxd) ~ stree_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  2.7833 0.08707 .
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.763332           0.013291       0.2924

``` r
# Plot residuals
plot(PD_PRic_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-14.png)<!-- -->

``` r
# Test relationship
PD_PRic_am_mxd <- aov(pd.obs ~ max_depth + Site_type, data = stree_sespd_env[surveyed_sites,])

# Log area
# Interaction
PD_PRic_am_lga <- aov(pd.obs ~ logArea * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PRic_am_lga)
```

    ##                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## logArea            1  2396779  2396779   6.242   0.0238 *  
    ## Site_type          2 23605634 11802817  30.738 3.31e-06 ***
    ## logArea:Site_type  2  4495320  2247660   5.854   0.0124 *  
    ## Residuals         16  6143625   383977                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_lm_lga <- lm(pd.obs ~ logArea + Site_type, data = stree_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_lm_lga)
    ## W = 0.97487, p-value = 0.8201

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_lm_lga) ~ stree_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  2.1613 0.1427
    ##       19

``` r
# Check for outliers
outlierTest(PD_PRic_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OOO -2.810007           0.012049      0.26508

``` r
# Plot residuals
plot(PD_PRic_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-15.png)<!-- -->

``` r
# Test relationship
PD_PRic_am_lga <- aov(pd.obs ~ logArea + Site_type, data = stree_sespd_env[surveyed_sites,])

# Anova outputs
summary(PD_PRic_am_dist)
```

    ##                         Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m  1 14204403 14204403  20.298 0.000274 ***
    ## Site_type                2  9840464  4920232   7.031 0.005541 ** 
    ## Residuals               18 12596491   699805                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_am_mxd)
```

    ##             Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## max_depth    1   897140   897140   1.423    0.248    
    ## Site_type    2 24398371 12199185  19.354 3.27e-05 ***
    ## Residuals   18 11345847   630325                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_am_lga)
```

    ##             Df   Sum Sq  Mean Sq F value  Pr(>F)    
    ## logArea      1  2396779  2396779   4.055  0.0592 .  
    ## Site_type    2 23605634 11802817  19.969 2.7e-05 ***
    ## Residuals   18 10638945   591053                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
dist_p_values <- summary(PD_PRic_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PRic_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PRic_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.0016427115 0.0332446855           NA 1.0000000000 0.0001962598
    ## [6]           NA 0.3554092488 0.0001617752           NA

``` r
## Mixed and stratified lakes
# Distance
# Interaction
PD_PRic_MS_am_dist <- aov(pd.obs ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_MS_am_dist)
```

    ##                                   Df   Sum Sq  Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m            1 11074223 11074223  13.577 0.00312 **
    ## Site_type                          1  9418046  9418046  11.546 0.00529 **
    ## distance_to_ocean_min_m:Site_type  1    32123    32123   0.039 0.84602   
    ## Residuals                         12  9788260   815688                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_MS_lm_dist <- lm(pd.obs ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_MS_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_MS_lm_dist)
    ## W = 0.9655, p-value = 0.7615

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_MS_lm_dist) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  1  14.581 0.001881 **
    ##       14                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_MS_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.060976          0.0098811       0.1581

``` r
# Plot residuals
plot(PD_PRic_MS_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-16.png)<!-- -->

``` r
# Test relationship
PD_PRic_MS_am_dist <- aov(pd.obs ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Max depth
# Interaction
PD_PRic_MS_am_mxd <- aov(pd.obs ~ max_depth * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_MS_am_mxd)
```

    ##                     Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## max_depth            1  1338297  1338297   2.199    0.164    
    ## Site_type            1 19922067 19922067  32.735 9.58e-05 ***
    ## max_depth:Site_type  1  1749176  1749176   2.874    0.116    
    ## Residuals           12  7303111   608593                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_MS_lm_mxd <- lm(pd.obs ~ max_depth + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_MS_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_MS_lm_mxd)
    ## W = 0.95538, p-value = 0.5794

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_MS_lm_mxd) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  4.9978 0.04218 *
    ##       14                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_MS_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.764629           0.017132      0.27411

``` r
# Plot residuals
plot(PD_PRic_MS_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-17.png)<!-- -->

``` r
# Test relationship
PD_PRic_MS_am_mxd <- aov(pd.obs ~ max_depth + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Log area
# Interaction
PD_PRic_MS_am_lga <- aov(pd.obs ~ logArea * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRic_MS_am_lga)
```

    ##                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## logArea            1  1001630  1001630   3.532   0.0847 .  
    ## Site_type          1 23476556 23476556  82.773 9.85e-07 ***
    ## logArea:Site_type  1  2430958  2430958   8.571   0.0127 *  
    ## Residuals         12  3403507   283626                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_MS_lm_lga <- lm(pd.obs ~ logArea + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_MS_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_MS_lm_lga)
    ## W = 0.95954, p-value = 0.6534

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_MS_lm_lga) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.2201 0.6462
    ##       14

``` r
# Check for outliers
outlierTest(PD_PRic_MS_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OLO -2.424884            0.03203      0.51249

``` r
# Plot residuals
plot(PD_PRic_MS_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-18.png)<!-- -->

``` r
# Test relationship
PD_PRic_MS_am_lga <- aov(pd.obs ~ logArea + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Anova outputs
summary(PD_PRic_MS_am_dist)
```

    ##                         Df   Sum Sq  Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m  1 11074223 11074223   14.66 0.00209 **
    ## Site_type                1  9418046  9418046   12.47 0.00369 **
    ## Residuals               13  9820383   755414                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_MS_am_mxd)
```

    ##             Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## max_depth    1  1338297  1338297   1.922 0.188962    
    ## Site_type    1 19922067 19922067  28.610 0.000132 ***
    ## Residuals   13  9052287   696330                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_MS_am_lga)
```

    ##             Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## logArea      1  1001630  1001630   2.232    0.159    
    ## Site_type    1 23476556 23476556  52.309 6.62e-06 ***
    ## Residuals   13  5834466   448805                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
dist_p_values <- summary(PD_PRic_MS_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PRic_MS_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PRic_MS_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.253904e-02 2.213284e-02           NA 1.000000e+00 7.935412e-04
    ## [6]           NA 9.544193e-01 3.974896e-05           NA

``` r
## Ocean sites and mixed lakes
# Distance
# Interaction
PD_PRic_OM_am_dist <- aov(pd.obs ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PRic_OM_am_dist)
```

    ##                                   Df   Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1     5800    5800   0.005  0.946
    ## Site_type                          1   349059  349059   0.289  0.602
    ## distance_to_ocean_min_m:Site_type  1      279     279   0.000  0.988
    ## Residuals                         10 12063630 1206363

``` r
# Linear model
PD_PRic_OM_lm_dist <- lm(pd.obs ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_OM_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_OM_lm_dist)
    ## W = 0.96068, p-value = 0.7343

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_OM_lm_dist) ~ stree_sespd_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.9063 0.1925
    ##       12

``` r
# Check for outliers
outlierTest(PD_PRic_OM_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.715308           0.021735      0.30429

``` r
# Plot residuals
plot(PD_PRic_OM_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-19.png)<!-- -->

``` r
# Test relationship
PD_PRic_OM_am_dist <- aov(pd.obs ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_mixed_sites,])

# Max depth
# Interaction
PD_PRic_OM_am_mxd <- aov(pd.obs ~ max_depth * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PRic_OM_am_mxd)
```

    ##                     Df  Sum Sq Mean Sq F value Pr(>F)
    ## max_depth            1 2723099 2723099   3.090  0.109
    ## Site_type            1  267683  267683   0.304  0.594
    ## max_depth:Site_type  1  614819  614819   0.698  0.423
    ## Residuals           10 8813167  881317

``` r
# Linear model
PD_PRic_OM_lm_mxd <- lm(pd.obs ~ max_depth + Site_type, data = stree_sespd_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_OM_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_OM_lm_mxd)
    ## W = 0.90557, p-value = 0.1359

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_OM_lm_mxd) ~ stree_sespd_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.7438 0.4054
    ##       12

``` r
# Check for outliers
outlierTest(PD_PRic_OM_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.310921           0.043439      0.60815

``` r
# Plot residuals
plot(PD_PRic_OM_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-20.png)<!-- -->

``` r
# Test relationship
PD_PRic_OM_am_mxd <- aov(pd.obs ~ max_depth + Site_type, data = stree_sespd_env[ocean_mixed_sites,])

# Log area
# Interaction
PD_PRic_OM_am_lga <- aov(pd.obs ~ logArea * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PRic_OM_am_lga)
```

    ##                   Df  Sum Sq Mean Sq F value Pr(>F)  
    ## logArea            1 1603123 1603123   2.988 0.1146  
    ## Site_type          1 1303825 1303825   2.430 0.1501  
    ## logArea:Site_type  1 4147050 4147050   7.730 0.0194 *
    ## Residuals         10 5364770  536477                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_OM_lm_lga <- lm(pd.obs ~ logArea + Site_type, data = stree_sespd_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_OM_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_OM_lm_lga)
    ## W = 0.96662, p-value = 0.8285

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_OM_lm_lga) ~ stree_sespd_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1596 0.6966
    ##       12

``` r
# Check for outliers
outlierTest(PD_PRic_OM_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OOO -2.628264           0.025237      0.35332

``` r
# Plot residuals
plot(PD_PRic_OM_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-21.png)<!-- -->

``` r
# Test relationship
PD_PRic_OM_am_lga <- aov(pd.obs ~ logArea + Site_type, data = stree_sespd_env[ocean_mixed_sites,])

# Anova outputs
summary(PD_PRic_OM_am_dist)
```

    ##                         Df   Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1     5800    5800   0.005  0.943
    ## Site_type                1   349059  349059   0.318  0.584
    ## Residuals               11 12063908 1096719

``` r
summary(PD_PRic_OM_am_mxd)
```

    ##             Df  Sum Sq Mean Sq F value Pr(>F)
    ## max_depth    1 2723099 2723099   3.177  0.102
    ## Site_type    1  267683  267683   0.312  0.587
    ## Residuals   11 9427986  857090

``` r
summary(PD_PRic_OM_am_lga)
```

    ##             Df  Sum Sq Mean Sq F value Pr(>F)
    ## logArea      1 1603123 1603123   1.854  0.201
    ## Site_type    1 1303825 1303825   1.508  0.245
    ## Residuals   11 9511820  864711

``` r
# p-values
dist_p_values <- summary(PD_PRic_OM_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PRic_OM_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PRic_OM_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000        NA 0.6135746 1.0000000        NA 1.0000000
    ## [8] 1.0000000        NA

``` r
## Stratified lakes and ocean sites
# Distance
# Interaction
PD_PRic_SO_am_dist <- aov(pd.obs ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PRic_SO_am_dist)
```

    ##                                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1 11333785 11333785  34.661 0.000154 ***
    ## Site_type                          1  2964993  2964993   9.068 0.013090 *  
    ## distance_to_ocean_min_m:Site_type  1     4385     4385   0.013 0.910100    
    ## Residuals                         10  3269879   326988                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_SO_lm_dist <- lm(pd.obs ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_SO_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_SO_lm_dist)
    ## W = 0.97663, p-value = 0.9505

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_SO_lm_dist) ~ stree_sespd_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  5.9308 0.03143 *
    ##       12                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_SO_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OOO -2.660798           0.023867      0.33413

``` r
# Plot residuals
plot(PD_PRic_SO_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-22.png)<!-- -->

``` r
# Test relationship
PD_PRic_SO_am_dist <- aov(pd.obs ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_stratified_sites,])

# Max depth
# Interaction
PD_PRic_SO_am_mxd <- aov(pd.obs ~ max_depth * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PRic_SO_am_mxd)
```

    ##                     Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## max_depth            1   997487   997487   3.242    0.102    
    ## Site_type            1 13254195 13254195  43.077 6.37e-05 ***
    ## max_depth:Site_type  1   244522   244522   0.795    0.394    
    ## Residuals           10  3076838   307684                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_SO_lm_mxd <- lm(pd.obs ~ max_depth + Site_type, data = stree_sespd_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_SO_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_SO_lm_mxd)
    ## W = 0.97603, p-value = 0.9452

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_SO_lm_mxd) ~ stree_sespd_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1   3.298 0.09441 .
    ##       12                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_SO_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OOO -2.407113           0.036865      0.51611

``` r
# Plot residuals
plot(PD_PRic_SO_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-23.png)<!-- -->

``` r
# Test relationship
PD_PRic_SO_am_mxd <- aov(pd.obs ~ max_depth + Site_type, data = stree_sespd_env[ocean_stratified_sites,])

# Log area
# Interaction
PD_PRic_SO_am_lga <- aov(pd.obs ~ logArea * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PRic_SO_am_lga)
```

    ##                   Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## logArea            1  2611857  2611857   7.422 0.021399 *  
    ## Site_type          1 11442004 11442004  32.515 0.000198 ***
    ## logArea:Site_type  1      209      209   0.001 0.981045    
    ## Residuals         10  3518972   351897                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRic_SO_lm_lga <- lm(pd.obs ~ logArea + Site_type, data = stree_sespd_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRic_SO_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRic_SO_lm_lga)
    ## W = 0.95582, p-value = 0.6541

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRic_SO_lm_lga) ~ stree_sespd_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)  
    ## group  1  3.6679 0.0796 .
    ##       12                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PRic_SO_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OOO -3.451756          0.0062079      0.08691

``` r
# Plot residuals
plot(PD_PRic_SO_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-24.png)<!-- -->

``` r
# Test relationship
PD_PRic_SO_am_lga <- aov(pd.obs ~ logArea + Site_type, data = stree_sespd_env[ocean_stratified_sites,])

# Anova outputs
summary(PD_PRic_SO_am_dist)
```

    ##                         Df   Sum Sq  Mean Sq F value  Pr(>F)    
    ## distance_to_ocean_min_m  1 11333785 11333785  38.076   7e-05 ***
    ## Site_type                1  2964993  2964993   9.961 0.00914 ** 
    ## Residuals               11  3274264   297660                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_SO_am_mxd)
```

    ##             Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## max_depth    1   997487   997487   3.304   0.0964 .  
    ## Site_type    1 13254195 13254195  43.897 3.73e-05 ***
    ## Residuals   11  3321360   301942                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRic_SO_am_lga)
```

    ##             Df   Sum Sq  Mean Sq F value   Pr(>F)    
    ## logArea      1  2611857  2611857   8.164   0.0156 *  
    ## Site_type    1 11442004 11442004  35.765 9.18e-05 ***
    ## Residuals   11  3519181   319926                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
dist_p_values <- summary(PD_PRic_SO_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PRic_SO_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PRic_SO_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.0004200309 0.0548531038           NA 0.5786528237 0.0002239774
    ## [6]           NA 0.0935568185 0.0005509114           NA

``` r
## Ocean sites
PD_PRic_geo_O_lm <- lm(pd.obs ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_sites,])
summary(PD_PRic_geo_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ distance_to_ocean_min_m + max_depth + logArea, 
    ##     data = stree_sespd_env[ocean_sites, ])
    ## 
    ## Residuals:
    ##        IBK        LCN        NCN        OOO        OOM        RCA 
    ##  3.468e+01 -3.594e+02  5.753e+01 -8.603e+02  1.127e+03  5.742e-14 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              3217.41    2652.24   1.213    0.349
    ## distance_to_ocean_min_m   -23.43      63.00  -0.372    0.746
    ## max_depth                  35.43      47.59   0.744    0.534
    ## logArea                   -53.45     236.02  -0.226    0.842
    ## 
    ## Residual standard error: 1036 on 2 degrees of freedom
    ## Multiple R-squared:  0.2306, Adjusted R-squared:  -0.9235 
    ## F-statistic: 0.1998 on 3 and 2 DF,  p-value: 0.8893

``` r
p_values <- summary(PD_PRic_geo_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                       1                       1                       1 
    ##                 logArea 
    ##                       1

``` r
## Mixed lakes
PD_PRic_geo_M_lm <- lm(pd.obs ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_lakes,])
summary(PD_PRic_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ distance_to_ocean_min_m + max_depth + logArea, 
    ##     data = stree_sespd_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##      FLK      HLO      LLN      MLN      NLN      NLU      OLO      ULN 
    ##   209.26  -644.99   653.52    91.81   242.45   478.25 -1131.39   101.09 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -3042.901   2192.435  -1.388   0.2375  
    ## distance_to_ocean_min_m    -6.087     12.967  -0.469   0.6632  
    ## max_depth                 -11.498     56.650  -0.203   0.8491  
    ## logArea                   707.357    290.791   2.433   0.0718 .
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
    ## distance_to_ocean_min_m  136260  1  0.2204 0.66319  
    ## max_depth                 25470  1  0.0412 0.84907  
    ## logArea                 3658461  1  5.9172 0.07178 .
    ## Residuals               2473103  4                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRic_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.9498828               1.0000000               1.0000000 
    ##                 logArea 
    ##               0.2871387

``` r
## Stratified lakes
PD_PRic_geo_S_lm <- lm(pd.obs ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[stratified_lakes,])
summary(PD_PRic_geo_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ distance_to_ocean_min_m + max_depth + logArea, 
    ##     data = stree_sespd_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##     BCM     CLM     GLK     HLM     NLK     OTM     SLN     TLN 
    ##  105.32  -31.80 -227.61  376.51   25.78 -424.77  -22.82  199.40 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              123.408   2012.127   0.061    0.954
    ## distance_to_ocean_min_m   -3.098      1.710  -1.812    0.144
    ## max_depth                 -7.210     24.285  -0.297    0.781
    ## logArea                  147.399    248.917   0.592    0.586
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
    ## distance_to_ocean_min_m 350555  1  3.2835 0.1442
    ## max_depth                 9410  1  0.0881 0.7813
    ## logArea                  37437  1  0.3507 0.5856
    ## Residuals               427048  4

``` r
p_values <- summary(PD_PRic_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               0.5768139               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
#### Environmental
### PRicz ANOVA
## Surveyed sites
# Temperature
# Interaction
PD_PRicz_am_temp <- aov(pd.obs.z ~ temperature_median * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRicz_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median            1  0.541  0.5408   0.718 0.4120  
    ## Site_type                     2  5.698  2.8488   3.784 0.0507 .
    ## temperature_median:Site_type  2  2.624  1.3121   1.743 0.2135  
    ## Residuals                    13  9.788  0.7529                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_lm_temp <- lm(pd.obs.z ~ temperature_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_lm_temp)
    ## W = 0.95365, p-value = 0.4548

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_lm_temp) ~ stree_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.3825 0.6882
    ##       16

``` r
# Check for outliers
outlierTest(PD_PRicz_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## IBK -2.137482           0.050684        0.963

``` r
# Plot residuals
plot(PD_PRicz_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-25.png)<!-- -->

``` r
# Test relationship
PD_PRicz_am_temp <- aov(pd.obs.z ~ temperature_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])

# Salinity
# Interaction
PD_PRicz_am_sal <- aov(pd.obs.z ~ salinity_median * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRicz_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median            1  1.819  1.8189   2.385 0.1465  
    ## Site_type                  2  4.664  2.3322   3.058 0.0816 .
    ## salinity_median:Site_type  2  2.253  1.1263   1.477 0.2643  
    ## Residuals                 13  9.914  0.7626                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_lm_sal <- lm(pd.obs.z ~ salinity_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_lm_sal)
    ## W = 0.95766, p-value = 0.5272

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_lm_sal) ~ stree_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.2363 0.7923
    ##       16

``` r
# Check for outliers
outlierTest(PD_PRicz_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NCN  1.87316           0.082076           NA

``` r
# Plot residuals
plot(PD_PRicz_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-26.png)<!-- -->

``` r
# Test relationship
PD_PRicz_am_sal <- aov(pd.obs.z ~ salinity_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])

# Oxygen
# Interaction
PD_PRicz_am_oxy <- aov(pd.obs.z ~ oxygen_median * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PRicz_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)  
    ## oxygen_median            1  1.808  1.8078   2.987 0.1076  
    ## Site_type                2  3.938  1.9688   3.253 0.0715 .
    ## oxygen_median:Site_type  2  5.038  2.5189   4.163 0.0401 *
    ## Residuals               13  7.867  0.6051                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_lm_oxy <- lm(pd.obs.z ~ oxygen_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_lm_oxy)
    ## W = 0.94023, p-value = 0.2663

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_lm_oxy) ~ stree_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.1531 0.8593
    ##       16

``` r
# Check for outliers
outlierTest(PD_PRicz_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## TLN -1.972283           0.068664           NA

``` r
# Plot residuals
plot(PD_PRicz_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-27.png)<!-- -->

``` r
# Test relationship
PD_PRicz_am_oxy <- aov(pd.obs.z ~ oxygen_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])

# Anova outputs
summary(PD_PRicz_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median  1  0.541  0.5408   0.654 0.4315  
    ## Site_type           2  5.698  2.8488   3.443 0.0588 .
    ## Residuals          15 12.412  0.8275                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRicz_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median  1  1.819  1.8189   2.242 0.1550  
    ## Site_type        2  4.664  2.3322   2.875 0.0877 .
    ## Residuals       15 12.167  0.8111                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRicz_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## oxygen_median  1  1.808  1.8078   2.101  0.168
    ## Site_type      2  3.938  1.9688   2.288  0.136
    ## Residuals     15 12.905  0.8603

``` r
# p-values
temp_p_values <- summary(PD_PRicz_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRicz_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRicz_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 0.3528936        NA 0.9301078 0.5261841        NA 1.0000000
    ## [8] 0.8142227        NA

``` r
## Mixed and stratified lakes
# Temperature
# Interaction
PD_PRicz_MS_am_temp <- aov(pd.obs.z ~ temperature_median * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRicz_MS_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median            1  1.018  1.0185   1.479  0.247
    ## Site_type                     1  0.148  0.1479   0.215  0.651
    ## temperature_median:Site_type  1  0.373  0.3734   0.542  0.476
    ## Residuals                    12  8.264  0.6887

``` r
# Linear model
PD_PRicz_MS_lm_temp <- lm(pd.obs.z ~ temperature_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_MS_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_MS_lm_temp)
    ## W = 0.95281, p-value = 0.5355

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_MS_lm_temp) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0446 0.8357
    ##       14

``` r
# Check for outliers
outlierTest(PD_PRicz_MS_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -2.082846           0.059323      0.94917

``` r
# Plot residuals
plot(PD_PRicz_MS_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-28.png)<!-- -->

``` r
# Test relationship
PD_PRicz_MS_am_temp <- aov(pd.obs.z ~ temperature_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Salinity
# Interaction
PD_PRicz_MS_am_sal <- aov(pd.obs.z ~ salinity_median * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRicz_MS_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)
    ## salinity_median            1  0.325  0.3254   0.541  0.476
    ## Site_type                  1  0.449  0.4492   0.746  0.405
    ## salinity_median:Site_type  1  1.807  1.8070   3.002  0.109
    ## Residuals                 12  7.222  0.6019

``` r
# Linear model
PD_PRicz_MS_lm_sal <- lm(pd.obs.z ~ salinity_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_MS_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_MS_lm_sal)
    ## W = 0.97155, p-value = 0.8631

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_MS_lm_sal) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0847 0.7753
    ##       14

``` r
# Check for outliers
outlierTest(PD_PRicz_MS_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -1.826192            0.09279           NA

``` r
# Plot residuals
plot(PD_PRicz_MS_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-29.png)<!-- -->

``` r
# Test relationship
PD_PRicz_MS_am_sal <- aov(pd.obs.z ~ salinity_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Oxygen
# Interaction
PD_PRicz_MS_am_oxy <- aov(pd.obs.z ~ oxygen_median * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRicz_MS_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)  
    ## oxygen_median            1  0.016  0.0158   0.026 0.8742  
    ## Site_type                1  0.023  0.0235   0.039 0.8472  
    ## oxygen_median:Site_type  1  2.500  2.4996   4.129 0.0649 .
    ## Residuals               12  7.265  0.6054                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_MS_lm_oxy <- lm(pd.obs.z ~ oxygen_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_MS_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_MS_lm_oxy)
    ## W = 0.96627, p-value = 0.7752

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_MS_lm_oxy) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0325 0.8596
    ##       14

``` r
# Check for outliers
outlierTest(PD_PRicz_MS_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -1.919228           0.079046           NA

``` r
# Plot residuals
plot(PD_PRicz_MS_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-30.png)<!-- -->

``` r
# Test relationship
PD_PRicz_MS_am_oxy <- aov(pd.obs.z ~ oxygen_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Anova outputs
summary(PD_PRicz_MS_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1  1.018  1.0185   1.533  0.238
    ## Site_type           1  0.148  0.1479   0.223  0.645
    ## Residuals          13  8.638  0.6644

``` r
summary(PD_PRicz_MS_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value Pr(>F)
    ## salinity_median  1  0.325  0.3254   0.468  0.506
    ## Site_type        1  0.449  0.4492   0.647  0.436
    ## Residuals       13  9.029  0.6946

``` r
summary(PD_PRicz_MS_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## oxygen_median  1  0.016  0.0158   0.021  0.887
    ## Site_type      1  0.023  0.0235   0.031  0.862
    ## Residuals     13  9.765  0.7511

``` r
# p-values
temp_p_values <- summary(PD_PRicz_MS_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRicz_MS_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRicz_MS_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1]  1  1 NA  1  1 NA  1  1 NA

``` r
## Ocean sites and mixed lakes
# Temperature
# Interaction
PD_PRicz_OM_am_temp <- aov(pd.obs.z ~ temperature_median * Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PRicz_OM_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median            1  0.272   0.272   0.345  0.576  
    ## Site_type                     1  4.678   4.678   5.934  0.045 *
    ## temperature_median:Site_type  1  2.515   2.515   3.190  0.117  
    ## Residuals                     7  5.518   0.788                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_OM_lm_temp <- lm(pd.obs.z ~ temperature_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_OM_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_OM_lm_temp)
    ## W = 0.95149, p-value = 0.663

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_OM_lm_temp) ~ stree_sespd_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.2594 0.6228
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRicz_OM_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## IBK -2.316377            0.05368      0.59048

``` r
# Plot residuals
plot(PD_PRicz_OM_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-31.png)<!-- -->

``` r
# Test relationship
PD_PRicz_OM_am_temp <- aov(pd.obs.z ~ temperature_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])

# Salinity
# Interaction
PD_PRicz_OM_am_sal <- aov(pd.obs.z ~ salinity_median * Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PRicz_OM_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median            1  5.647   5.647   7.138 0.0319 *
    ## Site_type                  1  1.024   1.024   1.294 0.2928  
    ## salinity_median:Site_type  1  0.774   0.774   0.978 0.3556  
    ## Residuals                  7  5.538   0.791                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_OM_lm_sal <- lm(pd.obs.z ~ salinity_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_OM_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_OM_lm_sal)
    ## W = 0.96289, p-value = 0.8072

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_OM_lm_sal) ~ stree_sespd_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.3639 0.2729
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRicz_OM_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NCN  2.51088           0.040345      0.44379

``` r
# Plot residuals
plot(PD_PRicz_OM_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-32.png)<!-- -->

``` r
# Test relationship
PD_PRicz_OM_am_sal <- aov(pd.obs.z ~ salinity_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])

# Oxygen
# Interaction
PD_PRicz_OM_am_oxy <- aov(pd.obs.z ~ oxygen_median * Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PRicz_OM_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value  Pr(>F)   
    ## oxygen_median            1  7.670   7.670  13.553 0.00785 **
    ## Site_type                1  0.291   0.291   0.514 0.49656   
    ## oxygen_median:Site_type  1  1.059   1.059   1.871 0.21365   
    ## Residuals                7  3.962   0.566                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_OM_lm_oxy <- lm(pd.obs.z ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_OM_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_OM_lm_oxy)
    ## W = 0.90199, p-value = 0.1955

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_OM_lm_oxy) ~ stree_sespd_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0037 0.9531
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRicz_OM_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NCN 1.705195            0.13193           NA

``` r
# Plot residuals
plot(PD_PRicz_OM_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-33.png)<!-- -->

``` r
# Test relationship
PD_PRicz_OM_am_oxy <- aov(pd.obs.z ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])

# Anova outputs
summary(PD_PRicz_OM_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median  1  0.272   0.272   0.271 0.6171  
    ## Site_type           1  4.678   4.678   4.659 0.0629 .
    ## Residuals           8  8.033   1.004                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRicz_OM_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median  1  5.647   5.647   7.158 0.0281 *
    ## Site_type        1  1.024   1.024   1.297 0.2876  
    ## Residuals        8  6.311   0.789                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRicz_OM_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value  Pr(>F)   
    ## oxygen_median  1  7.670   7.670  12.222 0.00813 **
    ## Site_type      1  0.291   0.291   0.464 0.51512   
    ## Residuals      8  5.021   0.628                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PRicz_OM_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRicz_OM_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRicz_OM_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.00000000 0.37767692         NA 0.16873602 1.00000000         NA 0.04875932
    ## [8] 1.00000000         NA

``` r
## Stratified lakes and ocean sites
# Temperature
# Interaction
PD_PRicz_SO_am_temp <- aov(pd.obs.z ~ temperature_median * Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PRicz_SO_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median            1  0.606   0.606   0.732 0.4205  
    ## Site_type                     1  4.506   4.506   5.445 0.0523 .
    ## temperature_median:Site_type  1  2.047   2.047   2.473 0.1598  
    ## Residuals                     7  5.793   0.828                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_SO_lm_temp <- lm(pd.obs.z ~ temperature_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_SO_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_SO_lm_temp)
    ## W = 0.94168, p-value = 0.5404

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_SO_lm_temp) ~ stree_sespd_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1   0.444 0.5219
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRicz_SO_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## IBK -1.985938           0.087413      0.96154

``` r
# Plot residuals
plot(PD_PRicz_SO_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-34.png)<!-- -->

``` r
# Test relationship
PD_PRicz_SO_am_temp <- aov(pd.obs.z ~ temperature_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])

# Salinity
# Interaction
PD_PRicz_SO_am_sal <- aov(pd.obs.z ~ salinity_median * Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PRicz_SO_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median            1  4.198   4.198   4.157 0.0809 .
    ## Site_type                  1  1.244   1.244   1.231 0.3038  
    ## salinity_median:Site_type  1  0.441   0.441   0.437 0.5297  
    ## Residuals                  7  7.069   1.010                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_SO_lm_sal <- lm(pd.obs.z ~ salinity_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_SO_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_SO_lm_sal)
    ## W = 0.94609, p-value = 0.5945

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_SO_lm_sal) ~ stree_sespd_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1   0.286 0.6058
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRicz_SO_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NCN 1.830114            0.10993           NA

``` r
# Plot residuals
plot(PD_PRicz_SO_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-35.png)<!-- -->

``` r
# Test relationship
PD_PRicz_SO_am_sal <- aov(pd.obs.z ~ salinity_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])

# Oxygen
# Interaction
PD_PRicz_SO_am_oxy <- aov(pd.obs.z ~ oxygen_median * Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PRicz_SO_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)  
    ## oxygen_median            1  2.413  2.4127   3.747 0.0941 .
    ## Site_type                1  2.901  2.9009   4.506 0.0714 .
    ## oxygen_median:Site_type  1  3.131  3.1307   4.863 0.0633 .
    ## Residuals                7  4.507  0.6438                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_SO_lm_oxy <- lm(pd.obs.z ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_SO_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_SO_lm_oxy)
    ## W = 0.96578, p-value = 0.841

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_SO_lm_oxy) ~ stree_sespd_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.6295 0.4479
    ##        9

``` r
# Check for outliers
outlierTest(PD_PRicz_SO_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## IBK -2.043276           0.080324      0.88356

``` r
# Plot residuals
plot(PD_PRicz_SO_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-36.png)<!-- -->

``` r
# Test relationship
PD_PRicz_SO_am_oxy <- aov(pd.obs.z ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])

# Anova outputs
summary(PD_PRicz_SO_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median  1  0.606   0.606   0.618 0.4543  
    ## Site_type           1  4.506   4.506   4.598 0.0643 .
    ## Residuals           8  7.839   0.980                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRicz_SO_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median  1  4.198   4.198   4.471 0.0674 .
    ## Site_type        1  1.244   1.244   1.325 0.2830  
    ## Residuals        8  7.510   0.939                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRicz_SO_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## oxygen_median  1  2.413  2.4127   2.527  0.151
    ## Site_type      1  2.901  2.9009   3.039  0.119
    ## Residuals      8  7.638  0.9547

``` r
# p-values
temp_p_values <- summary(PD_PRicz_SO_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PRicz_SO_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PRicz_SO_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 0.3860551        NA 0.4043595 1.0000000        NA 0.9034016
    ## [8] 0.7168322        NA

``` r
## Ocean sites
PD_PRicz_env_O_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_sites_env,])
summary(PD_PRicz_env_O_lm)
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
p_values <- summary(PD_PRicz_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
## Mixed lakes
PD_PRicz_env_M_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[mixed_lakes,])
summary(PD_PRicz_env_M_lm)
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
Anova(PD_PRicz_env_M_lm, type = 3)
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
p_values <- summary(PD_PRicz_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
## Stratified lakes
PD_PRicz_env_S_lm <- lm(pd.obs.z ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[stratified_lakes,])
summary(PD_PRicz_env_S_lm)
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
Anova(PD_PRicz_env_S_lm, type = 3)
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
p_values <- summary(PD_PRicz_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
#### Geographical
### PRicz ANOVA
## Surveyed sites
# Distance
# Interaction
PD_PRicz_am_dist <- aov(pd.obs.z ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PRicz_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1  1.310  1.3103   1.151  0.299
    ## Site_type                          2  0.756  0.3781   0.332  0.722
    ## distance_to_ocean_min_m:Site_type  2  1.224  0.6118   0.537  0.594
    ## Residuals                         16 18.216  1.1385

``` r
# Linear model
PD_PRicz_lm_dist <- lm(pd.obs.z ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_lm_dist)
    ## W = 0.99249, p-value = 0.9996

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_lm_dist) ~ stree_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2   0.546 0.5881
    ##       19

``` r
# Check for outliers
outlierTest(PD_PRicz_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LCN 2.532324           0.021476      0.47248

``` r
# Plot residuals
plot(PD_PRicz_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-37.png)<!-- -->

``` r
# Test relationship
PD_PRicz_am_dist <- aov(pd.obs.z ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[surveyed_sites,])

# Max depth
# Interaction
PD_PRicz_am_mxd <- aov(pd.obs.z ~ max_depth * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PRicz_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value  Pr(>F)   
    ## max_depth            1  1.600   1.600   2.848 0.11086   
    ## Site_type            2  2.956   1.478   2.631 0.10284   
    ## max_depth:Site_type  2  7.960   3.980   7.083 0.00626 **
    ## Residuals           16  8.990   0.562                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_lm_mxd <- lm(pd.obs.z ~ max_depth + Site_type, data = stree_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_lm_mxd)
    ## W = 0.96682, p-value = 0.6377

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_lm_mxd) ~ stree_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.3217 0.7288
    ##       19

``` r
# Check for outliers
outlierTest(PD_PRicz_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LCN 2.202306           0.041735      0.91817

``` r
# Plot residuals
plot(PD_PRicz_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-38.png)<!-- -->

``` r
# Test relationship
PD_PRicz_am_mxd <- aov(pd.obs.z ~ max_depth + Site_type, data = stree_sespd_env[surveyed_sites,])

# Log area
# Interaction
PD_PRicz_am_lga <- aov(pd.obs.z ~ logArea * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PRicz_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea            1  1.147  1.1473   1.143  0.301
    ## Site_type          2  1.002  0.5011   0.499  0.616
    ## logArea:Site_type  2  3.297  1.6485   1.642  0.225
    ## Residuals         16 16.060  1.0037

``` r
# Linear model
PD_PRicz_lm_lga <- lm(pd.obs.z ~ logArea + Site_type, data = stree_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_lm_lga)
    ## W = 0.98772, p-value = 0.9911

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_lm_lga) ~ stree_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.4736 0.6299
    ##       19

``` r
# Check for outliers
outlierTest(PD_PRicz_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## IBK -2.45073           0.025375      0.55826

``` r
# Plot residuals
plot(PD_PRicz_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-39.png)<!-- -->

``` r
# Test relationship
PD_PRicz_am_lga <- aov(pd.obs.z ~ logArea + Site_type, data = stree_sespd_env[surveyed_sites,])

# Anova outputs
summary(PD_PRicz_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1  1.310  1.3103   1.213  0.285
    ## Site_type                2  0.756  0.3781   0.350  0.709
    ## Residuals               18 19.440  1.0800

``` r
summary(PD_PRicz_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## max_depth    1  1.600  1.6005    1.70  0.209
    ## Site_type    2  2.956  1.4782    1.57  0.235
    ## Residuals   18 16.950  0.9416

``` r
summary(PD_PRicz_am_lga)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea      1  1.147  1.1473   1.067  0.315
    ## Site_type    2  1.002  0.5011   0.466  0.635
    ## Residuals   18 19.357  1.0754

``` r
# p-values
dist_p_values <- summary(PD_PRicz_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PRicz_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PRicz_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1]  1  1 NA  1  1 NA  1  1 NA

``` r
## Mixed and stratified lakes
# Distance
# Interaction
PD_PRicz_MS_am_dist <- aov(pd.obs.z ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRicz_MS_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1  0.201  0.2007   0.259  0.620
    ## Site_type                          1  0.130  0.1304   0.168  0.689
    ## distance_to_ocean_min_m:Site_type  1  0.171  0.1708   0.220  0.647
    ## Residuals                         12  9.302  0.7752

``` r
# Linear model
PD_PRicz_MS_lm_dist <- lm(pd.obs.z ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_MS_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_MS_lm_dist)
    ## W = 0.9547, p-value = 0.5675

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_MS_lm_dist) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1   0.053 0.8213
    ##       14

``` r
# Check for outliers
outlierTest(PD_PRicz_MS_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NLN   1.8388           0.090808           NA

``` r
# Plot residuals
plot(PD_PRicz_MS_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-40.png)<!-- -->

``` r
# Test relationship
PD_PRicz_MS_am_dist <- aov(pd.obs.z ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Max depth
# Interaction
PD_PRicz_MS_am_mxd <- aov(pd.obs.z ~ max_depth * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRicz_MS_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value Pr(>F)  
    ## max_depth            1  0.003  0.0030   0.005 0.9464  
    ## Site_type            1  0.001  0.0009   0.001 0.9705  
    ## max_depth:Site_type  1  2.139  2.1389   3.350 0.0921 .
    ## Residuals           12  7.661  0.6384                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_MS_lm_mxd <- lm(pd.obs.z ~ max_depth + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_MS_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_MS_lm_mxd)
    ## W = 0.9589, p-value = 0.6419

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_MS_lm_mxd) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0201 0.8893
    ##       14

``` r
# Check for outliers
outlierTest(PD_PRicz_MS_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -1.834797           0.091433           NA

``` r
# Plot residuals
plot(PD_PRicz_MS_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-41.png)<!-- -->

``` r
# Test relationship
PD_PRicz_MS_am_mxd <- aov(pd.obs.z ~ max_depth + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Log area
# Interaction
PD_PRicz_MS_am_lga <- aov(pd.obs.z ~ logArea * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PRicz_MS_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## logArea            1  0.184  0.1836   0.295 0.5970  
    ## Site_type          1  0.010  0.0104   0.017 0.8991  
    ## logArea:Site_type  1  2.141  2.1410   3.440 0.0884 .
    ## Residuals         12  7.469  0.6224                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_MS_lm_lga <- lm(pd.obs.z ~ logArea + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_MS_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_MS_lm_lga)
    ## W = 0.94029, p-value = 0.3526

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_MS_lm_lga) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.2073 0.6558
    ##       14

``` r
# Check for outliers
outlierTest(PD_PRicz_MS_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## HLO -1.828438           0.092434           NA

``` r
# Plot residuals
plot(PD_PRicz_MS_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-42.png)<!-- -->

``` r
# Test relationship
PD_PRicz_MS_am_lga <- aov(pd.obs.z ~ logArea + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Anova outputs
summary(PD_PRicz_MS_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1  0.201  0.2007   0.275  0.609
    ## Site_type                1  0.130  0.1304   0.179  0.679
    ## Residuals               13  9.473  0.7287

``` r
summary(PD_PRicz_MS_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## max_depth    1  0.003  0.0030   0.004  0.951
    ## Site_type    1  0.001  0.0009   0.001  0.973
    ## Residuals   13  9.800  0.7538

``` r
summary(PD_PRicz_MS_am_lga)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea      1  0.184  0.1836   0.248  0.627
    ## Site_type    1  0.010  0.0104   0.014  0.907
    ## Residuals   13  9.610  0.7392

``` r
# p-values
dist_p_values <- summary(PD_PRicz_MS_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PRicz_MS_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PRicz_MS_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1]  1  1 NA  1  1 NA  1  1 NA

``` r
## Ocean sites and mixed lakes
# Distance
# Interaction
PD_PRicz_OM_am_dist <- aov(pd.obs.z ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PRicz_OM_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1  1.387  1.3871   1.030  0.334
    ## Site_type                          1  0.175  0.1747   0.130  0.726
    ## distance_to_ocean_min_m:Site_type  1  1.195  1.1950   0.887  0.368
    ## Residuals                         10 13.469  1.3469

``` r
# Linear model
PD_PRicz_OM_lm_dist <- lm(pd.obs.z ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_OM_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_OM_lm_dist)
    ## W = 0.98458, p-value = 0.9932

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_OM_lm_dist) ~ stree_sespd_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.8484 0.3751
    ##       12

``` r
# Check for outliers
outlierTest(PD_PRicz_OM_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LCN 2.403997           0.037062      0.51886

``` r
# Plot residuals
plot(PD_PRicz_OM_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-43.png)<!-- -->

``` r
# Test relationship
PD_PRicz_OM_am_dist <- aov(pd.obs.z ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_mixed_sites,])

# Max depth
# Interaction
PD_PRicz_OM_am_mxd <- aov(pd.obs.z ~ max_depth * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PRicz_OM_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value  Pr(>F)   
    ## max_depth            1  9.229   9.229  19.247 0.00136 **
    ## Site_type            1  0.999   0.999   2.083 0.17957   
    ## max_depth:Site_type  1  1.203   1.203   2.509 0.14431   
    ## Residuals           10  4.795   0.480                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_OM_lm_mxd <- lm(pd.obs.z ~ max_depth + Site_type, data = stree_sespd_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_OM_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_OM_lm_mxd)
    ## W = 0.79399, p-value = 0.00421

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_OM_lm_mxd) ~ stree_sespd_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1   0.056 0.8169
    ##       12

``` r
# Check for outliers
outlierTest(PD_PRicz_OM_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NLN 2.963726           0.014197      0.19876

``` r
# Plot residuals
plot(PD_PRicz_OM_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-44.png)<!-- -->

``` r
# Test relationship
PD_PRicz_OM_am_mxd <- aov(pd.obs.z ~ max_depth + Site_type, data = stree_sespd_env[ocean_mixed_sites,])

# Log area
# Interaction
PD_PRicz_OM_am_lga <- aov(pd.obs.z ~ logArea * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PRicz_OM_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea            1  2.806  2.8064   2.129  0.175
    ## Site_type          1  0.214  0.2140   0.162  0.695
    ## logArea:Site_type  1  0.026  0.0264   0.020  0.890
    ## Residuals         10 13.180  1.3180

``` r
# Linear model
PD_PRicz_OM_lm_lga <- lm(pd.obs.z ~ logArea + Site_type, data = stree_sespd_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_OM_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_OM_lm_lga)
    ## W = 0.94771, p-value = 0.5257

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_OM_lm_lga) ~ stree_sespd_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.1471 0.1685
    ##       12

``` r
# Check for outliers
outlierTest(PD_PRicz_OM_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LCN  2.20534           0.051969      0.72757

``` r
# Plot residuals
plot(PD_PRicz_OM_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-45.png)<!-- -->

``` r
# Test relationship
PD_PRicz_OM_am_lga <- aov(pd.obs.z ~ logArea + Site_type, data = stree_sespd_env[ocean_mixed_sites,])

# Anova outputs
summary(PD_PRicz_OM_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1  1.387  1.3871   1.040  0.330
    ## Site_type                1  0.175  0.1747   0.131  0.724
    ## Residuals               11 14.665  1.3331

``` r
summary(PD_PRicz_OM_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## max_depth    1  9.229   9.229  16.926 0.00172 **
    ## Site_type    1  0.999   0.999   1.831 0.20311   
    ## Residuals   11  5.998   0.545                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PRicz_OM_am_lga)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea      1  2.806   2.806   2.338  0.155
    ## Site_type    1  0.214   0.214   0.178  0.681
    ## Residuals   11 13.206   1.200

``` r
# p-values
dist_p_values <- summary(PD_PRicz_OM_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PRicz_OM_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PRicz_OM_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.00000000 1.00000000         NA 0.01030662 1.00000000         NA 0.92711180
    ## [8] 1.00000000         NA

``` r
## Stratified lakes and ocean sites
# Distance
# Interaction
PD_PRicz_SO_am_dist <- aov(pd.obs.z ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PRicz_SO_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1  1.287  1.2869   0.942  0.355
    ## Site_type                          1  0.230  0.2297   0.168  0.690
    ## distance_to_ocean_min_m:Site_type  1  1.017  1.0171   0.744  0.408
    ## Residuals                         10 13.661  1.3661

``` r
# Linear model
PD_PRicz_SO_lm_dist <- lm(pd.obs.z ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_SO_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_SO_lm_dist)
    ## W = 0.98535, p-value = 0.9949

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_SO_lm_dist) ~ stree_sespd_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.6293  0.443
    ##       12

``` r
# Check for outliers
outlierTest(PD_PRicz_SO_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LCN 2.380078           0.038607       0.5405

``` r
# Plot residuals
plot(PD_PRicz_SO_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-46.png)<!-- -->

``` r
# Test relationship
PD_PRicz_SO_am_dist <- aov(pd.obs.z ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_stratified_sites,])

# Max depth
# Interaction
PD_PRicz_SO_am_mxd <- aov(pd.obs.z ~ max_depth * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PRicz_SO_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value  Pr(>F)   
    ## max_depth            1  0.408   0.408   0.739 0.41020   
    ## Site_type            1  2.490   2.490   4.508 0.05970 . 
    ## max_depth:Site_type  1  7.773   7.773  14.072 0.00378 **
    ## Residuals           10  5.524   0.552                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PRicz_SO_lm_mxd <- lm(pd.obs.z ~ max_depth + Site_type, data = stree_sespd_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_SO_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_SO_lm_mxd)
    ## W = 0.94225, p-value = 0.448

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_SO_lm_mxd) ~ stree_sespd_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0118 0.9151
    ##       12

``` r
# Check for outliers
outlierTest(PD_PRicz_SO_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LCN 2.121172           0.059905      0.83867

``` r
# Plot residuals
plot(PD_PRicz_SO_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-47.png)<!-- -->

``` r
# Test relationship
PD_PRicz_SO_am_mxd <- aov(pd.obs.z ~ max_depth + Site_type, data = stree_sespd_env[ocean_stratified_sites,])

# Log area
# Interaction
PD_PRicz_SO_am_lga <- aov(pd.obs.z ~ logArea * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PRicz_SO_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea            1  0.535   0.535   0.467  0.510
    ## Site_type          1  1.003   1.003   0.875  0.372
    ## logArea:Site_type  1  3.185   3.185   2.776  0.127
    ## Residuals         10 11.471   1.147

``` r
# Linear model
PD_PRicz_SO_lm_lga <- lm(pd.obs.z ~ logArea + Site_type, data = stree_sespd_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PRicz_SO_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PRicz_SO_lm_lga)
    ## W = 0.98747, p-value = 0.998

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PRicz_SO_lm_lga) ~ stree_sespd_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.3725  0.553
    ##       12

``` r
# Check for outliers
outlierTest(PD_PRicz_SO_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## IBK -2.564262           0.028167      0.39433

``` r
# Plot residuals
plot(PD_PRicz_SO_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-48.png)<!-- -->

``` r
# Test relationship
PD_PRicz_SO_am_lga <- aov(pd.obs.z ~ logArea + Site_type, data = stree_sespd_env[ocean_stratified_sites,])

# Anova outputs
summary(PD_PRicz_SO_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1  1.287  1.2869   0.964  0.347
    ## Site_type                1  0.230  0.2297   0.172  0.686
    ## Residuals               11 14.678  1.3344

``` r
summary(PD_PRicz_SO_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## max_depth    1  0.408  0.4081   0.338  0.573
    ## Site_type    1  2.490  2.4902   2.060  0.179
    ## Residuals   11 13.297  1.2088

``` r
summary(PD_PRicz_SO_am_lga)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea      1  0.535  0.5355   0.402  0.539
    ## Site_type    1  1.003  1.0034   0.753  0.404
    ## Residuals   11 14.656  1.3324

``` r
# p-values
dist_p_values <- summary(PD_PRicz_SO_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PRicz_SO_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PRicz_SO_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1]  1  1 NA  1  1 NA  1  1 NA

``` r
## Ocean sites
PD_PRicz_geo_O_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[ocean_sites,])
summary(PD_PRicz_geo_O_lm)
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
p_values <- summary(PD_PRicz_geo_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               1.0000000               0.3910896 
    ##                 logArea 
    ##               1.0000000

``` r
## Mixed lakes
PD_PRicz_geo_M_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[mixed_lakes,])
summary(PD_PRicz_geo_M_lm)
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
Anova(PD_PRicz_geo_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                         Sum Sq Df F value Pr(>F)
    ## (Intercept)             0.0717  1  0.0849 0.7852
    ## distance_to_ocean_min_m 0.0390  1  0.0462 0.8403
    ## max_depth               0.8626  1  1.0210 0.3694
    ## logArea                 0.0822  1  0.0973 0.7707
    ## Residuals               3.3794  4

``` r
p_values <- summary(PD_PRicz_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                       1                       1                       1 
    ##                 logArea 
    ##                       1

``` r
## Stratified lakes
PD_PRicz_geo_S_lm <- lm(pd.obs.z ~ distance_to_ocean_min_m + max_depth + logArea, stree_sespd_env[stratified_lakes,])
summary(PD_PRicz_geo_S_lm)
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
Anova(PD_PRicz_geo_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: pd.obs.z
    ##                         Sum Sq Df F value  Pr(>F)  
    ## (Intercept)             3.1755  1  7.2891 0.05410 .
    ## distance_to_ocean_min_m 0.0084  1  0.0194 0.89597  
    ## max_depth               1.1373  1  2.6106 0.18145  
    ## logArea                 2.4020  1  5.5137 0.07868 .
    ## Residuals               1.7426  4                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PRicz_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.2164135               1.0000000               0.7258090 
    ##                 logArea 
    ##               0.3147069

``` r
#### Environmental
### PDisp ANOVA
## Surveyed sites
# Temperature
# Interaction
PD_PDisp_am_temp <- aov(Dispersion ~ temperature_median * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_am_temp)
```

    ##                              Df   Sum Sq   Mean Sq F value Pr(>F)
    ## temperature_median            1 0.001013 0.0010125   1.963  0.185
    ## Site_type                     2 0.001623 0.0008117   1.573  0.244
    ## temperature_median:Site_type  2 0.000613 0.0003063   0.594  0.567
    ## Residuals                    13 0.006706 0.0005158

``` r
# Linear model
PD_PDisp_lm_temp <- lm(Dispersion ~ temperature_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_lm_temp)
    ## W = 0.96789, p-value = 0.7336

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_lm_temp) ~ stree_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2   4.485 0.02842 *
    ##       16                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_lm_temp)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -3.938722          0.0014843     0.028202

``` r
# Plot residuals
plot(PD_PDisp_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-49.png)<!-- -->

``` r
# Test relationship
PD_PDisp_am_temp <- aov(Dispersion ~ temperature_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])

# Salinity
# Interaction
PD_PDisp_am_sal <- aov(Dispersion ~ salinity_median * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_am_sal)
```

    ##                           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median            1 0.000811 0.0008107   1.193  0.294
    ## Site_type                  2 0.000283 0.0001413   0.208  0.815
    ## salinity_median:Site_type  2 0.000030 0.0000148   0.022  0.979
    ## Residuals                 13 0.008832 0.0006793

``` r
# Linear model
PD_PDisp_lm_sal <- lm(Dispersion ~ salinity_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_lm_sal)
    ## W = 0.84924, p-value = 0.006521

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_lm_sal) ~ stree_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  2.6351 0.1025
    ##       16

``` r
# Check for outliers
outlierTest(PD_PDisp_lm_sal)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -4.394898         0.00061078     0.011605

``` r
# Plot residuals
plot(PD_PDisp_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-50.png)<!-- -->

``` r
# Test relationship
PD_PDisp_am_sal <- aov(Dispersion ~ salinity_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])

# Oxygen
# Interaction
PD_PDisp_am_oxy <- aov(Dispersion ~ oxygen_median * Site_type, data = stree_sespd_env[surveyed_sites_env,])
summary(PD_PDisp_am_oxy)
```

    ##                         Df   Sum Sq   Mean Sq F value Pr(>F)
    ## oxygen_median            1 0.000100 0.0001001   0.144  0.711
    ## Site_type                2 0.000696 0.0003481   0.500  0.618
    ## oxygen_median:Site_type  2 0.000106 0.0000530   0.076  0.927
    ## Residuals               13 0.009052 0.0006963

``` r
# Linear model
PD_PDisp_lm_oxy <- lm(Dispersion ~ oxygen_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_lm_oxy)
    ## W = 0.85238, p-value = 0.007336

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_lm_oxy) ~ stree_sespd_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  3.2924 0.06345 .
    ##       16                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_lm_oxy)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -4.551226         0.00045264    0.0086001

``` r
# Plot residuals
plot(PD_PDisp_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-51.png)<!-- -->

``` r
# Test relationship
PD_PDisp_am_oxy <- aov(Dispersion ~ oxygen_median + Site_type, data = stree_sespd_env[surveyed_sites_env,])

# Anova outputs
summary(PD_PDisp_am_temp)
```

    ##                    Df   Sum Sq   Mean Sq F value Pr(>F)
    ## temperature_median  1 0.001013 0.0010125   2.075  0.170
    ## Site_type           2 0.001623 0.0008117   1.664  0.223
    ## Residuals          15 0.007319 0.0004879

``` r
summary(PD_PDisp_am_sal)
```

    ##                 Df   Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median  1 0.000811 0.0008107   1.372   0.26
    ## Site_type        2 0.000283 0.0001413   0.239   0.79
    ## Residuals       15 0.008861 0.0005907

``` r
summary(PD_PDisp_am_oxy)
```

    ##               Df   Sum Sq   Mean Sq F value Pr(>F)
    ## oxygen_median  1 0.000100 0.0001001   0.164  0.691
    ## Site_type      2 0.000696 0.0003481   0.570  0.577
    ## Residuals     15 0.009158 0.0006105

``` r
# p-values
temp_p_values <- summary(PD_PDisp_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PDisp_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PDisp_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1]  1  1 NA  1  1 NA  1  1 NA

``` r
## Mixed and stratified lakes
# Temperature
# Interaction
PD_PDisp_MS_am_temp <- aov(Dispersion ~ temperature_median * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_MS_am_temp)
```

    ##                              Df   Sum Sq   Mean Sq F value Pr(>F)
    ## temperature_median            1 0.001239 0.0012390   2.226  0.162
    ## Site_type                     1 0.001754 0.0017544   3.152  0.101
    ## temperature_median:Site_type  1 0.000122 0.0001220   0.219  0.648
    ## Residuals                    12 0.006680 0.0005567

``` r
# Linear model
PD_PDisp_MS_lm_temp <- lm(Dispersion ~ temperature_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_MS_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_MS_lm_temp)
    ## W = 0.97459, p-value = 0.9063

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_MS_lm_temp) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  1  9.2145 0.008902 **
    ##       14                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_MS_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -3.436283          0.0049282     0.078851

``` r
# Plot residuals
plot(PD_PDisp_MS_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-52.png)<!-- -->

``` r
# Test relationship
PD_PDisp_MS_am_temp <- aov(Dispersion ~ temperature_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Salinity
# Interaction
PD_PDisp_MS_am_sal <- aov(Dispersion ~ salinity_median * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_MS_am_sal)
```

    ##                           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median            1 0.001051 0.0010509   1.448  0.252
    ## Site_type                  1 0.000008 0.0000085   0.012  0.916
    ## salinity_median:Site_type  1 0.000029 0.0000293   0.040  0.844
    ## Residuals                 12 0.008707 0.0007256

``` r
# Linear model
PD_PDisp_MS_lm_sal <- lm(Dispersion ~ salinity_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_MS_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_MS_lm_sal)
    ## W = 0.85785, p-value = 0.01781

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_MS_lm_sal) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)  
    ## group  1  3.9465 0.0669 .
    ##       14                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_MS_lm_sal)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -4.139222          0.0013727     0.021964

``` r
# Plot residuals
plot(PD_PDisp_MS_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-53.png)<!-- -->

``` r
# Test relationship
PD_PDisp_MS_am_sal <- aov(Dispersion ~ salinity_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Oxygen
# Interaction
PD_PDisp_MS_am_oxy <- aov(Dispersion ~ oxygen_median * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_MS_am_oxy)
```

    ##                         Df   Sum Sq   Mean Sq F value Pr(>F)
    ## oxygen_median            1 0.000290 0.0002905   0.385  0.546
    ## Site_type                1 0.000446 0.0004456   0.591  0.457
    ## oxygen_median:Site_type  1 0.000010 0.0000104   0.014  0.908
    ## Residuals               12 0.009049 0.0007541

``` r
# Linear model
PD_PDisp_MS_lm_oxy <- lm(Dispersion ~ oxygen_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_MS_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_MS_lm_oxy)
    ## W = 0.86923, p-value = 0.02647

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_MS_lm_oxy) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  4.7725 0.04642 *
    ##       14                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_MS_lm_oxy)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -4.329353         0.00097982     0.015677

``` r
# Plot residuals
plot(PD_PDisp_MS_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-54.png)<!-- -->

``` r
# Test relationship
PD_PDisp_MS_am_oxy <- aov(Dispersion ~ oxygen_median + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Anova outputs
summary(PD_PDisp_MS_am_temp)
```

    ##                    Df   Sum Sq   Mean Sq F value Pr(>F)  
    ## temperature_median  1 0.001239 0.0012390   2.368 0.1478  
    ## Site_type           1 0.001754 0.0017544   3.353 0.0901 .
    ## Residuals          13 0.006802 0.0005232                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PDisp_MS_am_sal)
```

    ##                 Df   Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median  1 0.001051 0.0010509   1.564  0.233
    ## Site_type        1 0.000008 0.0000085   0.013  0.912
    ## Residuals       13 0.008736 0.0006720

``` r
summary(PD_PDisp_MS_am_oxy)
```

    ##               Df   Sum Sq   Mean Sq F value Pr(>F)
    ## oxygen_median  1 0.000290 0.0002905   0.417  0.530
    ## Site_type      1 0.000446 0.0004456   0.639  0.438
    ## Residuals     13 0.009060 0.0006969

``` r
# p-values
temp_p_values <- summary(PD_PDisp_MS_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PDisp_MS_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PDisp_MS_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.8869469 0.5405565        NA 1.0000000 1.0000000        NA 1.0000000
    ## [8] 1.0000000        NA

``` r
## Ocean sites and mixed lakes
# Temperature
# Interaction
PD_PDisp_OM_am_temp <- aov(Dispersion ~ temperature_median * Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PDisp_OM_am_temp)
```

    ##                              Df    Sum Sq   Mean Sq F value Pr(>F)
    ## temperature_median            1 0.0000158 1.581e-05   0.198  0.669
    ## Site_type                     1 0.0002154 2.154e-04   2.704  0.144
    ## temperature_median:Site_type  1 0.0001488 1.489e-04   1.869  0.214
    ## Residuals                     7 0.0005576 7.966e-05

``` r
# Linear model
PD_PDisp_OM_lm_temp <- lm(Dispersion ~ temperature_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_OM_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_OM_lm_temp)
    ## W = 0.866, p-value = 0.06885

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_OM_lm_temp) ~ stree_sespd_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0184  0.895
    ##        9

``` r
# Check for outliers
outlierTest(PD_PDisp_OM_lm_temp)
```

    ##     rstudent unadjusted p-value Bonferroni p
    ## NLN  5.46505          0.0009408     0.010349

``` r
# Plot residuals
plot(PD_PDisp_OM_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-55.png)<!-- -->

``` r
# Test relationship
PD_PDisp_OM_am_temp <- aov(Dispersion ~ temperature_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])

# Salinity
# Interaction
PD_PDisp_OM_am_sal <- aov(Dispersion ~ salinity_median * Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PDisp_OM_am_sal)
```

    ##                           Df    Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median            1 0.0001300 1.300e-04   1.307  0.291
    ## Site_type                  1 0.0001098 1.098e-04   1.104  0.328
    ## salinity_median:Site_type  1 0.0000018 1.790e-06   0.018  0.897
    ## Residuals                  7 0.0006960 9.943e-05

``` r
# Linear model
PD_PDisp_OM_lm_sal <- lm(Dispersion ~ salinity_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_OM_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_OM_lm_sal)
    ## W = 0.89544, p-value = 0.1625

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_OM_lm_sal) ~ stree_sespd_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0394 0.8472
    ##        9

``` r
# Check for outliers
outlierTest(PD_PDisp_OM_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NLN 3.684173          0.0078165     0.085981

``` r
# Plot residuals
plot(PD_PDisp_OM_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-56.png)<!-- -->

``` r
# Test relationship
PD_PDisp_OM_am_sal <- aov(Dispersion ~ salinity_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])

# Oxygen
# Interaction
PD_PDisp_OM_am_oxy <- aov(Dispersion ~ oxygen_median * Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
summary(PD_PDisp_OM_am_oxy)
```

    ##                         Df    Sum Sq   Mean Sq F value Pr(>F)  
    ## oxygen_median            1 0.0002858 2.857e-04   3.616  0.099 .
    ## Site_type                1 0.0000296 2.961e-05   0.375  0.560  
    ## oxygen_median:Site_type  1 0.0000691 6.905e-05   0.874  0.381  
    ## Residuals                7 0.0005532 7.903e-05                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PDisp_OM_lm_oxy <- lm(Dispersion ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_OM_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_OM_lm_oxy)
    ## W = 0.86393, p-value = 0.06474

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_OM_lm_oxy) ~ stree_sespd_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1368   0.72
    ##        9

``` r
# Check for outliers
outlierTest(PD_PDisp_OM_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NLN 3.899055          0.0059061     0.064967

``` r
# Plot residuals
plot(PD_PDisp_OM_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-57.png)<!-- -->

``` r
# Test relationship
PD_PDisp_OM_am_oxy <- aov(Dispersion ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_mixed_sites_env,])

# Anova outputs
summary(PD_PDisp_OM_am_temp)
```

    ##                    Df    Sum Sq   Mean Sq F value Pr(>F)
    ## temperature_median  1 0.0000158 1.581e-05   0.179  0.683
    ## Site_type           1 0.0002154 2.154e-04   2.439  0.157
    ## Residuals           8 0.0007064 8.830e-05

``` r
summary(PD_PDisp_OM_am_sal)
```

    ##                 Df    Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median  1 0.0001300 1.300e-04   1.490  0.257
    ## Site_type        1 0.0001098 1.098e-04   1.259  0.294
    ## Residuals        8 0.0006978 8.723e-05

``` r
summary(PD_PDisp_OM_am_oxy)
```

    ##               Df    Sum Sq   Mean Sq F value Pr(>F)  
    ## oxygen_median  1 0.0002858 2.857e-04   3.674 0.0916 .
    ## Site_type      1 0.0000296 2.961e-05   0.381 0.5543  
    ## Residuals      8 0.0006222 7.778e-05                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(PD_PDisp_OM_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PDisp_OM_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PDisp_OM_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 0.9419720        NA 1.0000000 1.0000000        NA 0.5494481
    ## [8] 1.0000000        NA

``` r
## Stratified lakes and ocean sites
# Temperature
# Interaction
PD_PDisp_SO_am_temp <- aov(Dispersion ~ temperature_median * Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PDisp_SO_am_temp)
```

    ##                              Df   Sum Sq   Mean Sq F value Pr(>F)
    ## temperature_median            1 0.001800 0.0017999   2.041  0.196
    ## Site_type                     1 0.000142 0.0001422   0.161  0.700
    ## temperature_median:Site_type  1 0.000536 0.0005361   0.608  0.461
    ## Residuals                     7 0.006174 0.0008820

``` r
# Linear model
PD_PDisp_SO_lm_temp <- lm(Dispersion ~ temperature_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_SO_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_SO_lm_temp)
    ## W = 0.95062, p-value = 0.6518

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_SO_lm_temp) ~ stree_sespd_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.1878 0.1732
    ##        9

``` r
# Check for outliers
outlierTest(PD_PDisp_SO_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -3.301368           0.013096      0.14406

``` r
# Plot residuals
plot(PD_PDisp_SO_lm_temp)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-58.png)<!-- -->

``` r
# Test relationship
PD_PDisp_SO_am_temp <- aov(Dispersion ~ temperature_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])

# Salinity
# Interaction
PD_PDisp_SO_am_sal <- aov(Dispersion ~ salinity_median * Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PDisp_SO_am_sal)
```

    ##                           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median            1 0.000288 0.0002880   0.244  0.636
    ## Site_type                  1 0.000104 0.0001043   0.088  0.775
    ## salinity_median:Site_type  1 0.000000 0.0000002   0.000  0.990
    ## Residuals                  7 0.008260 0.0011800

``` r
# Linear model
PD_PDisp_SO_lm_sal <- lm(Dispersion ~ salinity_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_SO_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_SO_lm_sal)
    ## W = 0.91654, p-value = 0.2908

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_SO_lm_sal) ~ stree_sespd_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.6659  0.229
    ##        9

``` r
# Check for outliers
outlierTest(PD_PDisp_SO_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -3.377305           0.011802      0.12982

``` r
# Plot residuals
plot(PD_PDisp_SO_lm_sal)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-59.png)<!-- -->

``` r
# Test relationship
PD_PDisp_SO_am_sal <- aov(Dispersion ~ salinity_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])

# Oxygen
# Interaction
PD_PDisp_SO_am_oxy <- aov(Dispersion ~ oxygen_median * Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
summary(PD_PDisp_SO_am_oxy)
```

    ##                         Df   Sum Sq   Mean Sq F value Pr(>F)
    ## oxygen_median            1 0.000001 0.0000007   0.001  0.982
    ## Site_type                1 0.000048 0.0000482   0.040  0.848
    ## oxygen_median:Site_type  1 0.000102 0.0001019   0.084  0.781
    ## Residuals                7 0.008502 0.0012145

``` r
# Linear model
PD_PDisp_SO_lm_oxy <- lm(Dispersion ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_SO_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_SO_lm_oxy)
    ## W = 0.9178, p-value = 0.3008

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_SO_lm_oxy) ~ stree_sespd_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.0511 0.1859
    ##        9

``` r
# Check for outliers
outlierTest(PD_PDisp_SO_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -3.574215          0.0090458     0.099504

``` r
# Plot residuals
plot(PD_PDisp_SO_lm_oxy)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-60.png)<!-- -->

``` r
# Test relationship
PD_PDisp_SO_am_oxy <- aov(Dispersion ~ oxygen_median + Site_type, data = stree_sespd_env[ocean_stratified_sites_env,])

# Anova outputs
summary(PD_PDisp_SO_am_temp)
```

    ##                    Df   Sum Sq   Mean Sq F value Pr(>F)
    ## temperature_median  1 0.001800 0.0017999   2.146  0.181
    ## Site_type           1 0.000142 0.0001422   0.170  0.691
    ## Residuals           8 0.006710 0.0008388

``` r
summary(PD_PDisp_SO_am_sal)
```

    ##                 Df   Sum Sq   Mean Sq F value Pr(>F)
    ## salinity_median  1 0.000288 0.0002880   0.279  0.612
    ## Site_type        1 0.000104 0.0001043   0.101  0.759
    ## Residuals        8 0.008260 0.0010325

``` r
summary(PD_PDisp_SO_am_oxy)
```

    ##               Df   Sum Sq   Mean Sq F value Pr(>F)
    ## oxygen_median  1 0.000001 0.0000007   0.001  0.980
    ## Site_type      1 0.000048 0.0000482   0.045  0.838
    ## Residuals      8 0.008604 0.0010755

``` r
# p-values
temp_p_values <- summary(PD_PDisp_SO_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(PD_PDisp_SO_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(PD_PDisp_SO_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1]  1  1 NA  1  1 NA  1  1 NA

``` r
## Ocean sites
PD_PDisp_env_O_lm <- lm(Dispersion ~ salinity_median + oxygen_median + temperature_median, stree_sespd_env[ocean_sites_env,])
summary(PD_PDisp_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = Dispersion ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = stree_sespd_env[ocean_sites_env, ])
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
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
p_values <- summary(PD_PDisp_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
## Mixed lakes
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
Anova(PD_PDisp_env_M_lm, type = 3)
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
p_values <- summary(PD_PDisp_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
## Stratified lakes
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
Anova(PD_PDisp_env_S_lm, type = 3)
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
p_values <- summary(PD_PDisp_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          1.0000000          1.0000000          0.7405587

``` r
#### Geographical
### PDisp ANOVA
## Surveyed sites
# Distance
# Interaction
PD_PDisp_am_dist <- aov(Dispersion ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PDisp_am_dist)
```

    ##                                   Df   Sum Sq   Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1 0.001143 0.0011432   2.217  0.156
    ## Site_type                          2 0.000823 0.0004115   0.798  0.467
    ## distance_to_ocean_min_m:Site_type  2 0.000046 0.0000228   0.044  0.957
    ## Residuals                         16 0.008252 0.0005157

``` r
# Linear model
PD_PDisp_lm_dist <- lm(Dispersion ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_lm_dist)
    ## W = 0.87644, p-value = 0.01034

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_lm_dist) ~ stree_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  4.0597 0.03404 *
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_lm_dist)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -4.693587         0.00020903    0.0045987

``` r
# Plot residuals
plot(PD_PDisp_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-61.png)<!-- -->

``` r
# Test relationship
PD_PDisp_am_dist <- aov(Dispersion ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[surveyed_sites,])

# Max depth
# Interaction
PD_PDisp_am_mxd <- aov(Dispersion ~ max_depth * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PDisp_am_mxd)
```

    ##                     Df   Sum Sq   Mean Sq F value Pr(>F)
    ## max_depth            1 0.000004 0.0000044   0.009  0.926
    ## Site_type            2 0.000904 0.0004519   0.896  0.428
    ## max_depth:Site_type  2 0.001284 0.0006421   1.273  0.307
    ## Residuals           16 0.008071 0.0005044

``` r
# Linear model
PD_PDisp_lm_mxd <- lm(Dispersion ~ max_depth + Site_type, data = stree_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_lm_mxd)
    ## W = 0.90144, p-value = 0.03177

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_lm_mxd) ~ stree_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)  
    ## group  2  2.8069 0.0855 .
    ##       19                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_lm_mxd)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -5.068094         9.5051e-05    0.0020911

``` r
# Plot residuals
plot(PD_PDisp_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-62.png)<!-- -->

``` r
# Test relationship
PD_PDisp_am_mxd <- aov(Dispersion ~ max_depth + Site_type, data = stree_sespd_env[surveyed_sites,])

# Log area
# Interaction
PD_PDisp_am_lga <- aov(Dispersion ~ logArea * Site_type, data = stree_sespd_env[surveyed_sites,])
summary(PD_PDisp_am_lga)
```

    ##                   Df   Sum Sq   Mean Sq F value Pr(>F)  
    ## logArea            1 0.000054 0.0000535   0.129 0.7244  
    ## Site_type          2 0.000817 0.0004087   0.983 0.3956  
    ## logArea:Site_type  2 0.002741 0.0013705   3.297 0.0633 .
    ## Residuals         16 0.006652 0.0004157                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PDisp_lm_lga <- lm(Dispersion ~ logArea + Site_type, data = stree_sespd_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_lm_lga)
    ## W = 0.89993, p-value = 0.02964

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_lm_lga) ~ stree_sespd_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  2.7222 0.09129 .
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_lm_lga)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -4.817792         0.00016071    0.0035356

``` r
# Plot residuals
plot(PD_PDisp_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-63.png)<!-- -->

``` r
# Test relationship
PD_PDisp_am_lga <- aov(Dispersion ~ logArea + Site_type, data = stree_sespd_env[surveyed_sites,])

# Anova outputs
summary(PD_PDisp_am_dist)
```

    ##                         Df   Sum Sq   Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.001143 0.0011432   2.480  0.133
    ## Site_type                2 0.000823 0.0004115   0.893  0.427
    ## Residuals               18 0.008297 0.0004610

``` r
summary(PD_PDisp_am_mxd)
```

    ##             Df   Sum Sq   Mean Sq F value Pr(>F)
    ## max_depth    1 0.000004 0.0000044   0.009  0.927
    ## Site_type    2 0.000904 0.0004519   0.869  0.436
    ## Residuals   18 0.009355 0.0005197

``` r
summary(PD_PDisp_am_lga)
```

    ##             Df   Sum Sq   Mean Sq F value Pr(>F)
    ## logArea      1 0.000054 0.0000535   0.103  0.752
    ## Site_type    2 0.000817 0.0004087   0.783  0.472
    ## Residuals   18 0.009393 0.0005218

``` r
# p-values
dist_p_values <- summary(PD_PDisp_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PDisp_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PDisp_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.796261 1.000000       NA 1.000000 1.000000       NA 1.000000 1.000000
    ## [9]       NA

``` r
## Mixed and stratified lakes
# Distance
# Interaction
PD_PDisp_MS_am_dist <- aov(Dispersion ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_MS_am_dist)
```

    ##                                   Df   Sum Sq   Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1 0.001969 0.0019693   3.035  0.107
    ## Site_type                          1 0.000003 0.0000028   0.004  0.949
    ## distance_to_ocean_min_m:Site_type  1 0.000038 0.0000383   0.059  0.812
    ## Residuals                         12 0.007785 0.0006488

``` r
# Linear model
PD_PDisp_MS_lm_dist <- lm(Dispersion ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_MS_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_MS_lm_dist)
    ## W = 0.88826, p-value = 0.05228

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_MS_lm_dist) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  5.2903 0.03734 *
    ##       14                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_MS_lm_dist)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -4.224355          0.0011798     0.018876

``` r
# Plot residuals
plot(PD_PDisp_MS_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-64.png)<!-- -->

``` r
# Test relationship
PD_PDisp_MS_am_dist <- aov(Dispersion ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Max depth
# Interaction
PD_PDisp_MS_am_mxd <- aov(Dispersion ~ max_depth * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_MS_am_mxd)
```

    ##                     Df   Sum Sq   Mean Sq F value Pr(>F)
    ## max_depth            1 0.000136 0.0001355   0.207  0.657
    ## Site_type            1 0.001329 0.0013292   2.027  0.180
    ## max_depth:Site_type  1 0.000462 0.0004621   0.705  0.418
    ## Residuals           12 0.007869 0.0006557

``` r
# Linear model
PD_PDisp_MS_lm_mxd <- lm(Dispersion ~ max_depth + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_MS_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_MS_lm_mxd)
    ## W = 0.949, p-value = 0.4741

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_MS_lm_mxd) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.8758 0.06911 .
    ##       14                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_MS_lm_mxd)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -4.088236          0.0015037     0.024059

``` r
# Plot residuals
plot(PD_PDisp_MS_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-65.png)<!-- -->

``` r
# Test relationship
PD_PDisp_MS_am_mxd <- aov(Dispersion ~ max_depth + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Log area
# Interaction
PD_PDisp_MS_am_lga <- aov(Dispersion ~ logArea * Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
summary(PD_PDisp_MS_am_lga)
```

    ##                   Df   Sum Sq   Mean Sq F value Pr(>F)
    ## logArea            1 0.000762 0.0007622   1.437  0.254
    ## Site_type          1 0.001153 0.0011529   2.173  0.166
    ## logArea:Site_type  1 0.001514 0.0015144   2.855  0.117
    ## Residuals         12 0.006366 0.0005305

``` r
# Linear model
PD_PDisp_MS_lm_lga <- lm(Dispersion ~ logArea + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_MS_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_MS_lm_lga)
    ## W = 0.9575, p-value = 0.6168

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_MS_lm_lga) ~ stree_sespd_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.2161 0.1588
    ##       14

``` r
# Check for outliers
outlierTest(PD_PDisp_MS_lm_lga)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -3.779688          0.0026251     0.042001

``` r
# Plot residuals
plot(PD_PDisp_MS_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-66.png)<!-- -->

``` r
# Test relationship
PD_PDisp_MS_am_lga <- aov(Dispersion ~ logArea + Site_type, data = stree_sespd_env[mixed_stratified_lakes,])

# Anova outputs
summary(PD_PDisp_MS_am_dist)
```

    ##                         Df   Sum Sq   Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m  1 0.001969 0.0019693   3.272 0.0936 .
    ## Site_type                1 0.000003 0.0000028   0.005 0.9466  
    ## Residuals               13 0.007824 0.0006018                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(PD_PDisp_MS_am_mxd)
```

    ##             Df   Sum Sq   Mean Sq F value Pr(>F)
    ## max_depth    1 0.000136 0.0001355   0.212  0.653
    ## Site_type    1 0.001329 0.0013292   2.074  0.173
    ## Residuals   13 0.008331 0.0006408

``` r
summary(PD_PDisp_MS_am_lga)
```

    ##             Df   Sum Sq   Mean Sq F value Pr(>F)
    ## logArea      1 0.000762 0.0007622   1.257  0.282
    ## Site_type    1 0.001153 0.0011529   1.902  0.191
    ## Residuals   13 0.007881 0.0006062

``` r
# p-values
dist_p_values <- summary(PD_PDisp_MS_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PDisp_MS_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PDisp_MS_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.5618147 1.0000000        NA 1.0000000 1.0000000        NA 1.0000000
    ## [8] 1.0000000        NA

``` r
## Ocean sites and mixed lakes
# Distance
# Interaction
PD_PDisp_OM_am_dist <- aov(Dispersion ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PDisp_OM_am_dist)
```

    ##                                   Df    Sum Sq   Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1 0.0000369 3.687e-05   0.374  0.554
    ## Site_type                          1 0.0001969 1.969e-04   2.000  0.188
    ## distance_to_ocean_min_m:Site_type  1 0.0000015 1.510e-06   0.015  0.904
    ## Residuals                         10 0.0009846 9.846e-05

``` r
# Linear model
PD_PDisp_OM_lm_dist <- lm(Dispersion ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_OM_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_OM_lm_dist)
    ## W = 0.93105, p-value = 0.3156

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_OM_lm_dist) ~ stree_sespd_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1364 0.7183
    ##       12

``` r
# Check for outliers
outlierTest(PD_PDisp_OM_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NLN 2.670446           0.023475      0.32864

``` r
# Plot residuals
plot(PD_PDisp_OM_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-67.png)<!-- -->

``` r
# Test relationship
PD_PDisp_OM_am_dist <- aov(Dispersion ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_mixed_sites,])

# Max depth
# Interaction
PD_PDisp_OM_am_mxd <- aov(Dispersion ~ max_depth * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PDisp_OM_am_mxd)
```

    ##                     Df    Sum Sq   Mean Sq F value Pr(>F)
    ## max_depth            1 0.0001741 1.741e-04   2.222  0.167
    ## Site_type            1 0.0001497 1.497e-04   1.911  0.197
    ## max_depth:Site_type  1 0.0001129 1.129e-04   1.441  0.258
    ## Residuals           10 0.0007834 7.834e-05

``` r
# Linear model
PD_PDisp_OM_lm_mxd <- lm(Dispersion ~ max_depth + Site_type, data = stree_sespd_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_OM_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_OM_lm_mxd)
    ## W = 0.89216, p-value = 0.08684

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_OM_lm_mxd) ~ stree_sespd_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1157 0.7396
    ##       12

``` r
# Check for outliers
outlierTest(PD_PDisp_OM_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NLN 3.699459          0.0041122      0.05757

``` r
# Plot residuals
plot(PD_PDisp_OM_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-68.png)<!-- -->

``` r
# Test relationship
PD_PDisp_OM_am_mxd <- aov(Dispersion ~ max_depth + Site_type, data = stree_sespd_env[ocean_mixed_sites,])

# Log area
# Interaction
PD_PDisp_OM_am_lga <- aov(Dispersion ~ logArea * Site_type, data = stree_sespd_env[ocean_mixed_sites,])
summary(PD_PDisp_OM_am_lga)
```

    ##                   Df    Sum Sq   Mean Sq F value Pr(>F)
    ## logArea            1 0.0002153 2.153e-04   2.504  0.145
    ## Site_type          1 0.0000486 4.861e-05   0.565  0.469
    ## logArea:Site_type  1 0.0000964 9.639e-05   1.121  0.315
    ## Residuals         10 0.0008597 8.597e-05

``` r
# Linear model
PD_PDisp_OM_lm_lga <- lm(Dispersion ~ logArea + Site_type, data = stree_sespd_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_OM_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_OM_lm_lga)
    ## W = 0.81747, p-value = 0.008277

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_OM_lm_lga) ~ stree_sespd_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0073 0.9333
    ##       12

``` r
# Check for outliers
outlierTest(PD_PDisp_OM_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## NLN 3.533486          0.0054152     0.075813

``` r
# Plot residuals
plot(PD_PDisp_OM_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-69.png)<!-- -->

``` r
# Test relationship
PD_PDisp_OM_am_lga <- aov(Dispersion ~ logArea + Site_type, data = stree_sespd_env[ocean_mixed_sites,])

# Anova outputs
summary(PD_PDisp_OM_am_dist)
```

    ##                         Df    Sum Sq   Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.0000369 3.687e-05   0.411  0.534
    ## Site_type                1 0.0001969 1.969e-04   2.197  0.166
    ## Residuals               11 0.0009862 8.965e-05

``` r
summary(PD_PDisp_OM_am_mxd)
```

    ##             Df    Sum Sq   Mean Sq F value Pr(>F)
    ## max_depth    1 0.0001741 1.741e-04   2.136  0.172
    ## Site_type    1 0.0001497 1.497e-04   1.837  0.202
    ## Residuals   11 0.0008962 8.147e-05

``` r
summary(PD_PDisp_OM_am_lga)
```

    ##             Df    Sum Sq   Mean Sq F value Pr(>F)
    ## logArea      1 0.0002153 2.153e-04   2.477  0.144
    ## Site_type    1 0.0000486 4.861e-05   0.559  0.470
    ## Residuals   11 0.0009561 8.692e-05

``` r
# p-values
dist_p_values <- summary(PD_PDisp_OM_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PDisp_OM_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PDisp_OM_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 0.9982671        NA 1.0000000 1.0000000        NA 0.8630313
    ## [8] 1.0000000        NA

``` r
## Stratified lakes and ocean sites
# Distance
# Interaction
PD_PDisp_SO_am_dist <- aov(Dispersion ~ distance_to_ocean_min_m * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PDisp_SO_am_dist)
```

    ##                                   Df   Sum Sq   Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1 0.000885 0.0008849   1.144  0.310
    ## Site_type                          1 0.000480 0.0004800   0.621  0.449
    ## distance_to_ocean_min_m:Site_type  1 0.000009 0.0000087   0.011  0.918
    ## Residuals                         10 0.007734 0.0007734

``` r
# Linear model
PD_PDisp_SO_lm_dist <- lm(Dispersion ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_SO_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_SO_lm_dist)
    ## W = 0.9079, p-value = 0.1469

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_SO_lm_dist) ~ stree_sespd_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.5969 0.08221 .
    ##       12                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(PD_PDisp_SO_lm_dist)
```

    ##     rstudent unadjusted p-value Bonferroni p
    ## SLN -3.87417          0.0030881     0.043234

``` r
# Plot residuals
plot(PD_PDisp_SO_lm_dist)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-70.png)<!-- -->

``` r
# Test relationship
PD_PDisp_SO_am_dist <- aov(Dispersion ~ distance_to_ocean_min_m + Site_type, data = stree_sespd_env[ocean_stratified_sites,])

# Max depth
# Interaction
PD_PDisp_SO_am_mxd <- aov(Dispersion ~ max_depth * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PDisp_SO_am_mxd)
```

    ##                     Df   Sum Sq   Mean Sq F value Pr(>F)
    ## max_depth            1 0.000119 0.0001188   0.159  0.699
    ## Site_type            1 0.000294 0.0002942   0.393  0.545
    ## max_depth:Site_type  1 0.001204 0.0012044   1.608  0.234
    ## Residuals           10 0.007490 0.0007490

``` r
# Linear model
PD_PDisp_SO_lm_mxd <- lm(Dispersion ~ max_depth + Site_type, data = stree_sespd_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_SO_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_SO_lm_mxd)
    ## W = 0.94428, p-value = 0.4759

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_SO_lm_mxd) ~ stree_sespd_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.7799 0.2069
    ##       12

``` r
# Check for outliers
outlierTest(PD_PDisp_SO_lm_mxd)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -4.346481          0.0014516     0.020322

``` r
# Plot residuals
plot(PD_PDisp_SO_lm_mxd)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-71.png)<!-- -->

``` r
# Test relationship
PD_PDisp_SO_am_mxd <- aov(Dispersion ~ max_depth + Site_type, data = stree_sespd_env[ocean_stratified_sites,])

# Log area
# Interaction
PD_PDisp_SO_am_lga <- aov(Dispersion ~ logArea * Site_type, data = stree_sespd_env[ocean_stratified_sites,])
summary(PD_PDisp_SO_am_lga)
```

    ##                   Df   Sum Sq   Mean Sq F value Pr(>F)  
    ## logArea            1 0.000277 0.0002765   0.455 0.5153  
    ## Site_type          1 0.000033 0.0000330   0.054 0.8205  
    ## logArea:Site_type  1 0.002720 0.0027205   4.476 0.0605 .
    ## Residuals         10 0.006077 0.0006077                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
PD_PDisp_SO_lm_lga <- lm(Dispersion ~ logArea + Site_type, data = stree_sespd_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(PD_PDisp_SO_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(PD_PDisp_SO_lm_lga)
    ## W = 0.93533, p-value = 0.3618

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(PD_PDisp_SO_lm_lga) ~ stree_sespd_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.6714 0.2204
    ##       12

``` r
# Check for outliers
outlierTest(PD_PDisp_SO_lm_lga)
```

    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -4.073984          0.0022358     0.031302

``` r
# Plot residuals
plot(PD_PDisp_SO_lm_lga)
```

![](PD_analyses_files/figure-gfm/PD%20alpha%20lm%20and%20ANOVA-72.png)<!-- -->

``` r
# Test relationship
PD_PDisp_SO_am_lga <- aov(Dispersion ~ logArea + Site_type, data = stree_sespd_env[ocean_stratified_sites,])

# Anova outputs
summary(PD_PDisp_SO_am_dist)
```

    ##                         Df   Sum Sq   Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.000885 0.0008849   1.257  0.286
    ## Site_type                1 0.000480 0.0004800   0.682  0.426
    ## Residuals               11 0.007743 0.0007039

``` r
summary(PD_PDisp_SO_am_mxd)
```

    ##             Df   Sum Sq   Mean Sq F value Pr(>F)
    ## max_depth    1 0.000119 0.0001188   0.150  0.706
    ## Site_type    1 0.000294 0.0002942   0.372  0.554
    ## Residuals   11 0.008695 0.0007904

``` r
summary(PD_PDisp_SO_am_lga)
```

    ##             Df   Sum Sq   Mean Sq F value Pr(>F)
    ## logArea      1 0.000277 0.0002765   0.346  0.568
    ## Site_type    1 0.000033 0.0000330   0.041  0.843
    ## Residuals   11 0.008798 0.0007998

``` r
# p-values
dist_p_values <- summary(PD_PDisp_SO_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(PD_PDisp_SO_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(PD_PDisp_SO_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1]  1  1 NA  1  1 NA  1  1 NA

``` r
## Ocean sites
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
## Mixed lakes
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
Anova(PD_PDisp_geo_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                            Sum Sq Df F value  Pr(>F)   
    ## (Intercept)             0.0035139  1 36.2705 0.00383 **
    ## distance_to_ocean_min_m 0.0001786  1  1.8433 0.24611   
    ## max_depth               0.0001244  1  1.2836 0.32054   
    ## logArea                 0.0000959  1  0.9898 0.37609   
    ## Residuals               0.0003875  4                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              0.01531847              0.98443978              1.00000000 
    ##                 logArea 
    ##              1.00000000

``` r
## Stratified lakes
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
Anova(PD_PDisp_geo_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: Dispersion
    ##                            Sum Sq Df F value  Pr(>F)  
    ## (Intercept)             0.0009336  1  1.7478 0.25669  
    ## distance_to_ocean_min_m 0.0028211  1  5.2810 0.08312 .
    ## max_depth               0.0011743  1  2.1983 0.21231  
    ## logArea                 0.0033590  1  6.2879 0.06623 .
    ## Residuals               0.0021368  4                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(PD_PDisp_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               0.3324746               0.8492304 
    ##                 logArea 
    ##               0.2649182

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
    ## 3 Stratified vs Reference  1 0.7729479 7.110647 0.5039207   0.098      0.588
    ## 4          Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.132      0.792
    ## 5      Mixed vs Reference  1 0.6592270 3.967203 0.3617334   0.113      0.678
    ## 6      Ocean vs Reference  1 0.6319662 3.354877 0.4015472   0.164      0.984
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
    ## 2 Stratified vs Ocean  1 1.1858559 8.357071 0.4105242   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.2573592 1.467099 0.1089395   0.095      0.285

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
    ## 3      Mixed vs Ocean  1 0.3347638 2.034034 0.1560556   0.018      0.054

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
    ## 1 Stratified vs Mixed  1 1.3499033 10.500850 0.4666868   0.003      0.009   *
    ## 2 Stratified vs Ocean  1 1.2146330  9.192717 0.4789690   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.2573592  1.467099 0.1089395   0.101      0.303

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
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31832 0.699  0.015 *
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
    ## env[surveyed_sites_env, c(34)] 0.94798 0.31832 0.699  0.015 *
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
    ## temperature_median -0.75816 -0.65206 0.0724  0.690   
    ## salinity_median     0.94798  0.31832 0.6990  0.009 **
    ## oxygen_median       0.57185 -0.82036 0.6836  0.122   
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
    ## temperature_median -0.75816 -0.65206 0.0724  1.000  
    ## salinity_median     0.94798  0.31832 0.6990  0.027 *
    ## oxygen_median       0.57185 -0.82036 0.6836  0.366  
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
    ## temperature_median -0.92674  0.37571 0.0954  0.714  
    ## salinity_median     0.90472  0.42600 0.6840  0.031 *
    ## oxygen_median       0.44478 -0.89564 0.5836  0.201  
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
    ## temperature_median -0.92674  0.37571 0.0954  1.000  
    ## salinity_median     0.90472  0.42600 0.6840  0.093 .
    ## oxygen_median       0.44478 -0.89564 0.5836  0.603  
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
    ## temperature_median -0.52519 -0.85098 0.0223  0.900  
    ## salinity_median    -0.79785 -0.60286 0.7269  0.019 *
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.016 *
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
    ## salinity_median    -0.79785 -0.60286 0.7269  0.057 .
    ## oxygen_median      -0.91393 -0.40587 0.8026  0.048 *
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
    ## temperature_median  0.0002545 -1.0000000 0.0757  0.504
    ## salinity_median    -0.0051338 -0.9999900 0.4840  0.847
    ## oxygen_median      -0.0011559  1.0000000 0.7647  0.144
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
    ## temperature_median  0.0002545 -1.0000000 0.0757  1.000
    ## salinity_median    -0.0051338 -0.9999900 0.4840  1.000
    ## oxygen_median      -0.0011559  1.0000000 0.7647  0.432
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
    ## temperature_median -0.00016277  1.00000000 0.1197  0.741  
    ## salinity_median     0.00052765 -1.00000000 0.5577  0.104  
    ## oxygen_median       0.00068652 -1.00000000 0.7434  0.033 *
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
    ## temperature_median -0.00016277  1.00000000 0.1197  1.000  
    ## salinity_median     0.00052765 -1.00000000 0.5577  0.312  
    ## oxygen_median       0.00068652 -1.00000000 0.7434  0.099 .
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
    ## temperature_median       -0.39653 -0.91802 0.0381  0.901  
    ## salinity_median          -0.96164  0.27433 0.5830  0.095 .
    ## oxygen_median             0.83739  0.54660 0.7450  0.026 *
    ## distance_to_ocean_mean_m -0.09987  0.99500 0.0319  0.918  
    ## max_depth                -0.87444  0.48513 0.1022  0.774  
    ## logArea                  -0.40923  0.91243 0.2149  0.552  
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
    ## temperature_median       -0.39653 -0.91802 0.0381  1.000
    ## salinity_median          -0.96164  0.27433 0.5830  0.570
    ## oxygen_median             0.83739  0.54660 0.7450  0.156
    ## distance_to_ocean_mean_m -0.09987  0.99500 0.0319  1.000
    ## max_depth                -0.87444  0.48513 0.1022  1.000
    ## logArea                  -0.40923  0.91243 0.2149  1.000
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
    ## distance_to_ocean_min_m -0.82791  0.56086 0.5558  0.070 .
    ## max_depth               -0.20870 -0.97798 0.0858  0.868  
    ## logArea                  0.25301 -0.96746 0.2012  0.300  
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
    ## distance_to_ocean_min_m -0.82791  0.56086 0.5558   0.21
    ## max_depth               -0.20870 -0.97798 0.0858   1.00
    ## logArea                  0.25301 -0.96746 0.2012   0.90
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
    ## distance_to_ocean_min_m -0.82791  0.56086 0.5558  0.068 .
    ## max_depth               -0.20870 -0.97798 0.0858  0.878  
    ## logArea                  0.25301 -0.96746 0.2012  0.301  
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
    ## max_depth               -0.20870 -0.97798 0.0858  1.000
    ## logArea                  0.25301 -0.96746 0.2012  0.903
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
    ## distance_to_ocean_min_m -0.76299  0.64641 0.4671  0.124
    ## max_depth               -0.34719 -0.93780 0.0777  0.966
    ## logArea                  0.55404  0.83249 0.0074  0.960
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
    ## distance_to_ocean_min_m -0.76299  0.64641 0.4671  0.372
    ## max_depth               -0.34719 -0.93780 0.0777  1.000
    ## logArea                  0.55404  0.83249 0.0074  1.000
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
    ## distance_to_ocean_min_m  0.45615 -0.88990 0.4273  0.172   
    ## max_depth               -0.94119 -0.33789 0.6682  0.002 **
    ## logArea                 -0.86397  0.50354 0.2339  0.455   
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
    ## distance_to_ocean_min_m  0.45615 -0.88990 0.4273  0.516   
    ## max_depth               -0.94119 -0.33789 0.6682  0.006 **
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
    ## distance_to_ocean_min_m  0.99474  0.10241 0.6338  0.228
    ## max_depth                0.73039  0.68303 0.0374  0.974
    ## logArea                 -0.39145  0.92020 0.3078  0.345
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
    ## distance_to_ocean_min_m  0.99474  0.10241 0.6338  0.684
    ## max_depth                0.73039  0.68303 0.0374  1.000
    ## logArea                 -0.39145  0.92020 0.3078  1.000
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
    ## distance_to_ocean_min_m -0.00022472 -1.00000000 0.5336  0.149  
    ## max_depth                0.00068729  1.00000000 0.7569  0.045 *
    ## logArea                  0.00053980 -1.00000000 0.7316  0.034 *
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
    ## distance_to_ocean_min_m -0.00022472 -1.00000000 0.5336  0.447
    ## max_depth                0.00068729  1.00000000 0.7569  0.135
    ## logArea                  0.00053980 -1.00000000 0.7316  0.102
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
    ##       Significance: 0.2 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.213 0.243 0.261 0.274 
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
    ## 0.427 0.458 0.482 0.500 
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
    ##       Significance: 0.288 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.524 0.564 0.589 0.611 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_mant_pv <- rbind(PD_beta_env_mant_t$signif, PD_beta_env_mant_s$signif, PD_beta_env_mant_o$signif)
PD_beta_env_mant_pv <- PD_beta_env_mant_pv[,1]
(PD_beta_env_mant_pv <- p.adjust(PD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.600 0.015 0.864

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
    ##       Significance: 0.23 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.215 0.246 0.284 0.315 
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
    ##       Significance: 0.008 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.375 0.403 0.424 0.455 
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
    ##       Significance: 0.321 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.414 0.455 0.492 0.525 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_MS_mant_pv <- rbind(PD_beta_env_MS_mant_t$signif, PD_beta_env_MS_mant_s$signif, PD_beta_env_MS_mant_o$signif)
PD_beta_env_MS_mant_pv <- PD_beta_env_MS_mant_pv[,1]
(PD_beta_env_MS_mant_pv <- p.adjust(PD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.690 0.024 0.963

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
    ##       Significance: 0.504 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.142 0.199 0.247 0.373 
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
    ##       Significance: 0.013 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.291 0.356 0.424 0.465 
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
    ## 0.338 0.421 0.473 0.611 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_OM_mant_pv <- rbind(PD_beta_env_OM_mant_t$signif, PD_beta_env_OM_mant_s$signif, PD_beta_env_OM_mant_o$signif)
PD_beta_env_OM_mant_pv <- PD_beta_env_OM_mant_pv[,1]
(PD_beta_env_OM_mant_pv <- p.adjust(PD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.039 0.042

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
    ## 0.0390 0.0735 0.1046 0.1247 
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
    ##       Significance: 0.014 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.391 0.427 0.451 0.501 
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
    ##       Significance: 0.218 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.690 0.718 0.739 0.776 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_SO_mant_pv <- rbind(PD_beta_env_SO_mant_t$signif, PD_beta_env_SO_mant_s$signif, PD_beta_env_SO_mant_o$signif)
PD_beta_env_SO_mant_pv <- PD_beta_env_SO_mant_pv[,1]
(PD_beta_env_SO_mant_pv <- p.adjust(PD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.042 0.654

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
    ##       Significance: 0.779 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.258 0.354 0.452 0.688 
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
    ## 0.214 0.286 0.325 0.369 
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
    ##       Significance: 0.084 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.241 0.353 0.461 0.564 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_env_M_mant_pv <- rbind(PD_beta_env_M_mant_t$signif, PD_beta_env_M_mant_s$signif, PD_beta_env_M_mant_o$signif)
PD_beta_env_M_mant_pv <- PD_beta_env_M_mant_pv[,1]
(PD_beta_env_M_mant_pv <- p.adjust(PD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.702 0.252

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
    ##       Significance: 0.268 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.346 0.439 0.545 0.637 
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
    ##       Significance: 0.583 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.460 0.483 0.509 0.535 
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
    ##       Significance: 0.769 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.164 0.193 0.214 0.245 
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
    ##       Significance: 0.704 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.111 0.128 0.143 0.167 
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
    ##       Significance: 0.625 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.296 0.340 0.374 0.424 
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
    ## 0.154 0.193 0.224 0.262 
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
    ##       Significance: 0.874 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.128 0.153 0.178 0.206 
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
    ##       Significance: 0.392 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.166 0.195 0.216 0.254 
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
    ##       Significance: 0.023 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.181 0.238 0.286 0.403 
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
    ##       Significance: 0.064 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.169 0.206 0.232 0.264 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_OM_mant_pv <- rbind( PD_beta_geo_OM_mant_dmean$signif, PD_beta_geo_OM_mant_md$signif, PD_beta_geo_OM_mant_la$signif)
PD_beta_geo_OM_mant_pv <- PD_beta_geo_OM_mant_pv[,1]
(PD_beta_geo_OM_mant_pv <- p.adjust(PD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.069 0.192

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
    ##       Significance: 0.521 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.622 0.657 0.681 0.696 
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
    ##       Significance: 0.785 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.116 0.140 0.164 0.192 
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
    ##       Significance: 0.198 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.370 0.399 0.421 0.441 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_SO_mant_pv <- rbind( PD_beta_geo_SO_mant_dmean$signif, PD_beta_geo_SO_mant_md$signif, PD_beta_geo_SO_mant_la$signif)
PD_beta_geo_SO_mant_pv <- PD_beta_geo_SO_mant_pv[,1]
(PD_beta_geo_SO_mant_pv <- p.adjust(PD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.594

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
    ##       Significance: 0.702 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.228 0.303 0.364 0.664 
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
    ##       Significance: 0.017 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.297 0.449 0.548 0.624 
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
    ## 0.203 0.292 0.344 0.418 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
PD_beta_geo_M_mant_pv <- rbind( PD_beta_geo_M_mant_dmean$signif, PD_beta_geo_M_mant_md$signif, PD_beta_geo_M_mant_la$signif)
PD_beta_geo_M_mant_pv <- PD_beta_geo_M_mant_pv[,1]
(PD_beta_geo_M_mant_pv <- p.adjust(PD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.051 0.051

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
    ##  [7] ggplot2_3.5.1        picante_1.8.2        nlme_3.1-167        
    ## [10] vegan_2.6-10         lattice_0.22-6       permute_0.9-7       
    ## [13] car_3.1-3            carData_3.0-5        tidyr_1.3.1         
    ## [16] phytools_2.4-4       maps_3.4.2.1         ape_5.8-1           
    ## [19] reshape2_1.4.4       stringr_1.5.1        dplyr_1.1.4         
    ## [22] knitr_1.49          
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
