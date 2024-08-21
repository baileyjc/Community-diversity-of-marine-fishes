Trait diversity analyses
================

#### R Markdown

## Load packages

``` r
# Load the knitr package if not already loaded
library(knitr)

# Source the R Markdown file
knit("/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.Rmd", output = "/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md")
```

    ## 
    ## 
    ## processing file: /Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.Rmd

    ##   |                  |          |   0%  |                  |          |   3%                                                                          |                  |.         |   6% [Bringing everything together load modifying files packages]             |                  |.         |   9%                                                                          |                  |.         |  12% [Bringing everything together load in modifying files]                   |                  |..        |  16%                                                                          |                  |..        |  19% [Check species names across files]                                       |                  |..        |  22%                                                                          |                  |..        |  25% [Modify environment data]                                                |                  |...       |  28%                                                                          |                  |...       |  31% [Modify incidence matrices]                                              |                  |...       |  34%                                                                          |                  |....      |  38% [Modify phylogeny]                                                       |                  |....      |  41%                                                                          |                  |....      |  44% [Modify trait data]                                                      |                  |.....     |  47%                                                                          |                  |.....     |  50% [Modify community data frames]                                           |                  |.....     |  53%                                                                          |                  |......    |  56% [Modify community trait data]                                            |                  |......    |  59%                                                                          |                  |......    |  62% [Community trait data tests]                                             |                  |.......   |  66%                                                                          |                  |.......   |  69% [Modify stratification data frames]                                      |                  |.......   |  72%                                                                          |                  |........  |  75% [Modify stratification trait data]                                       |                  |........  |  78%                                                                          |                  |........  |  81% [Stratification trait data tests]                                        |                  |........  |  84%                                                                          |                  |......... |  88% [Modify site trait data frames]                                          |                  |......... |  91%                                                                          |                  |......... |  94% [Site trait data tests]                                                  |                  |..........|  97%                                                                          |                  |..........| 100% [Bringing everything together load out modified files and session info]

    ## output file: /Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md

    ## [1] "/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md"

``` r
library(picante)
```

    ## Loading required package: vegan

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.6-4

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

# Define your custom colors
custom_colors <- c("Reference" = "black", "Ocean" = "#EE6363", "Mixed" = "#87CEFA", "Stratified" = "#6E8B3D")
```

## Functional Diversity MNTD, MPD, and PD

``` r
# Trait distance calculation
traits_dist <- gowdis(traits)
traits_clust <- hclust(traits_dist,"average")
traits_tree <- as.phylo(traits_clust)


#Mean pairwise differences
traits_sesmpd <- ses.mpd(presabs_lake, traits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(traits_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sesmpd.csv")


#Mean nearest taxon distance
traits_sesmntd <- ses.mntd(presabs_lake, traits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(traits_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sesmntd.csv")


#Faith's PD
traits_sespd <- ses.pd(presabs_lake, traits_tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)

write.csv(traits_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sespd.csv")


# Trait distance calculation
straits_dist <- gowdis(straits)
straits_clust <- hclust(straits_dist,"average")
straits_tree <- as.phylo(straits_clust)


#Mean pairwise differences
straits_sesmpd <- ses.mpd(surveyed_sites_lake, straits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(straits_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sesmpd.csv")


#Mean nearest taxon distance
straits_sesmntd <- ses.mntd(surveyed_sites_lake, straits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(straits_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sesmntd.csv")


#Faith's PD
straits_sespd <- ses.pd(surveyed_sites_lake, straits_tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)

write.csv(straits_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sespd.csv")
```

## Subset of traits Functional Diversity MNTD, MPD, and PD

``` r
# Trait distance calculation
keep <- c("BodyShapeI", "OperculumPresent", "Troph", "RepGuild1", "ParentalCare", "WaterPref", "DorsalSpinesMeann")
traits_subset <- traits[,keep]
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
straits_subset <- straits[,keep]
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

## Read in FD mpd, mntd, pd

- Read in Functional Diversity files and combine with env data

``` r
traits_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sesmpd.csv")

traits_sesmpd_env <- merge(traits_sesmpd, env, by = "X", sort = F)
traits_sesmpd_env$measure <- "traits_sesmpd"
traits_sesmpd_env$Stratification <- factor(traits_sesmpd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


traits_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sesmntd.csv")

traits_sesmntd_env <- merge(traits_sesmntd, env, by = "X", sort = F)
traits_sesmntd_env$measure <- "traits_sesmntd"
traits_sesmntd_env$Stratification <- factor(traits_sesmntd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


traits_sespd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sespd.csv")

traits_sespd_env <- merge(traits_sespd, env, by = "X", sort = F)
traits_sespd_env$measure <- "traits_sespd"
traits_sespd_env$Stratification <- factor(traits_sespd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


straits_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sesmpd.csv")

straits_sesmpd_env <- merge(straits_sesmpd, env, by = "X", sort = F)
straits_sesmpd_env$measure <- "straits_sesmpd"
straits_sesmpd_env$Stratification <- factor(straits_sesmpd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))


straits_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sesmntd.csv")

straits_sesmntd_env <- merge(straits_sesmntd, env, by = "X", sort = F)
straits_sesmntd_env$measure <- "straits_sesmntd"
straits_sesmntd_env$Stratification <- factor(straits_sesmntd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))


straits_sespd <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sespd.csv")

straits_sespd_env <- merge(straits_sespd, env, by = "X", sort = F)
straits_sespd_env$measure <- "straits_sespd"
straits_sespd_env$Stratification <- factor(straits_sespd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))
```

## FR complete island Linear Models

``` r
FR_model <- lm(log(pd.obs) ~ volume_m3 + distance_to_ocean_mean_m + tidal_lag_time_minutes + max_depth + logArea, data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(FR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ volume_m3 + distance_to_ocean_mean_m + 
    ##     tidal_lag_time_minutes + max_depth + logArea, data = traits_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##        1        2        4        5       12       17       21       22 
    ## -0.02894 -0.28497  0.13805  0.75170 -0.02437 -0.60545 -0.11837  0.17235 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              -8.364e-01  5.405e+00  -0.155    0.891
    ## volume_m3                 1.133e-06  1.831e-06   0.619    0.599
    ## distance_to_ocean_mean_m -4.152e-03  4.178e-03  -0.994    0.425
    ## tidal_lag_time_minutes    5.555e-03  1.102e-02   0.504    0.664
    ## max_depth                -5.577e-02  8.410e-02  -0.663    0.575
    ## logArea                   1.296e-01  5.665e-01   0.229    0.840
    ## 
    ## Residual standard error: 0.7338 on 2 degrees of freedom
    ## Multiple R-squared:  0.3904, Adjusted R-squared:  -1.134 
    ## F-statistic: 0.2561 on 5 and 2 DF,  p-value: 0.9048

``` r
p_values <- summary(FR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(FR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                          Df  Sum Sq Mean Sq F value Pr(>F)
    ## volume_m3                 1 0.00405 0.00405  0.0075 0.9388
    ## distance_to_ocean_mean_m  1 0.34723 0.34723  0.6448 0.5062
    ## tidal_lag_time_minutes    1 0.08157 0.08157  0.1515 0.7347
    ## max_depth                 1 0.22864 0.22864  0.4246 0.5815
    ## logArea                   1 0.02817 0.02817  0.0523 0.8403
    ## Residuals                 2 1.07703 0.53851

``` r
# Get R-squared value
r_squared <- summary(FR_model)$r.squared

# Get number of observations
n <- nrow(traits_sespd_env[c(1:2,4,5,12,17,21:22),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(FR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 0.2561359

``` r
FR_model <- lm(pd.obs ~ temperature_median + salinity_median + oxygen_median + pH_median, data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(FR_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = traits_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##        1        2        4        5       12       17       21       22 
    ##  0.18487 -0.17062  0.27164 -0.02230 -0.13209 -0.28022 -0.08017  0.22889 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -6.65529    5.52857  -1.204    0.315
    ## temperature_median -0.07557    0.09853  -0.767    0.499
    ## salinity_median     0.02337    0.06690   0.349    0.750
    ## oxygen_median      -0.26654    0.24478  -1.089    0.356
    ## pH_median           1.34567    1.01008   1.332    0.275
    ## 
    ## Residual standard error: 0.3122 on 3 degrees of freedom
    ## Multiple R-squared:  0.806,  Adjusted R-squared:  0.5474 
    ## F-statistic: 3.117 on 4 and 3 DF,  p-value: 0.1887

``` r
p_values <- summary(FR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(FR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##                    Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## temperature_median  1 0.01523 0.01523  0.1563 0.71903  
    ## salinity_median     1 0.96913 0.96913  9.9444 0.05113 .
    ## oxygen_median       1 0.05758 0.05758  0.5908 0.49807  
    ## pH_median           1 0.17297 0.17297  1.7749 0.27494  
    ## Residuals           3 0.29236 0.09745                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get R-squared value
r_squared <- summary(FR_model)$r.squared

# Get number of observations
n <- nrow(traits_sespd_env[c(1:2,4,5,12,17,21:22),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(FR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 3.116592

``` r
FR_model <- lm(pd.obs ~ salinity_median, data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(FR_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ salinity_median, data = traits_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.37381 -0.19949  0.01921  0.15600  0.39777 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)     -2.09412    0.88619  -2.363   0.0560 .
    ## salinity_median  0.10723    0.03205   3.346   0.0155 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2961 on 6 degrees of freedom
    ## Multiple R-squared:  0.6511, Adjusted R-squared:  0.5929 
    ## F-statistic: 11.19 on 1 and 6 DF,  p-value: 0.0155

``` r
anova(FR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##                 Df  Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median  1 0.98132 0.98132  11.195 0.0155 *
    ## Residuals        6 0.52595 0.08766                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(traits_sespd_env[c(1:2,4,5,12,17,21:22),3], traits_sespd_env[c(1:2,4,5,12,17,21:22),14])
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# Island biogeography common relationships
fr_model <- lm(log(pd.obs) ~ logArea, data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(fr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ logArea, data = traits_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.61822 -0.33452 -0.06349  0.19407  0.74465 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept) -0.50232    1.77165  -0.284    0.786
    ## logArea      0.02178    0.17166   0.127    0.903
    ## 
    ## Residual standard error: 0.5419 on 6 degrees of freedom
    ## Multiple R-squared:  0.002676,   Adjusted R-squared:  -0.1635 
    ## F-statistic: 0.0161 on 1 and 6 DF,  p-value: 0.9032

``` r
anova(fr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## logArea    1 0.00473 0.004728  0.0161 0.9032
    ## Residuals  6 1.76197 0.293661

``` r
frz_model <- lm(pd.obs.z ~ logArea, data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(frz_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ logArea, data = traits_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.09318 -0.58791  0.01385  0.43173  1.19272 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)  -1.9332     2.8311  -0.683    0.520
    ## logArea       0.1677     0.2743   0.612    0.563
    ## 
    ## Residual standard error: 0.866 on 6 degrees of freedom
    ## Multiple R-squared:  0.05867,    Adjusted R-squared:  -0.09822 
    ## F-statistic: 0.374 on 1 and 6 DF,  p-value: 0.5633

``` r
anova(frz_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##           Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea    1 0.2804 0.28044   0.374 0.5633
    ## Residuals  6 4.4993 0.74989

``` r
fr_model <- lm(log(pd.obs) ~ log(distance_to_ocean_mean_m), data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(fr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ log(distance_to_ocean_mean_m), data = traits_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5345 -0.3486 -0.1077  0.2601  0.8661 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                     2.0620     2.6516   0.778    0.466
    ## log(distance_to_ocean_mean_m)  -0.4234     0.4785  -0.885    0.410
    ## 
    ## Residual standard error: 0.5104 on 6 degrees of freedom
    ## Multiple R-squared:  0.1154, Adjusted R-squared:  -0.03199 
    ## F-statistic: 0.783 on 1 and 6 DF,  p-value: 0.4103

``` r
anova(fr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                               Df  Sum Sq Mean Sq F value Pr(>F)
    ## log(distance_to_ocean_mean_m)  1 0.20394 0.20394   0.783 0.4103
    ## Residuals                      6 1.56275 0.26046

``` r
fr_model <- lm(log(pd.obs) ~ log(max_depth), data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(fr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ log(max_depth), data = traits_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5749 -0.3222 -0.0520  0.1950  0.7971 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)     0.03507    1.23752   0.028    0.978
    ## log(max_depth) -0.10099    0.39336  -0.257    0.806
    ## 
    ## Residual standard error: 0.5397 on 6 degrees of freedom
    ## Multiple R-squared:  0.01087,    Adjusted R-squared:  -0.154 
    ## F-statistic: 0.06591 on 1 and 6 DF,  p-value: 0.806

``` r
anova(fr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                Df Sum Sq  Mean Sq F value Pr(>F)
    ## log(max_depth)  1 0.0192 0.019197  0.0659  0.806
    ## Residuals       6 1.7475 0.291250

## FR habitat island Linear Models

``` r
FR_model <- lm(log(pd.obs) ~ volume_m3_w_chemocline + distance_to_ocean_mean_m + tidal_lag_time_minutes + max_depth + logArea, data = traits_sespd_env[c(3,6,8,10,13:15,23),])
summary(FR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ volume_m3_w_chemocline + distance_to_ocean_mean_m + 
    ##     tidal_lag_time_minutes + max_depth + logArea, data = traits_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##         3         6        10        13        14        15        23 
    ##  0.120455  0.005801  0.053594 -0.098313 -0.053415 -0.089692  0.061570 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)               4.694e-01  1.061e+00   0.443    0.735
    ## volume_m3_w_chemocline    2.300e-06  9.835e-07   2.338    0.257
    ## distance_to_ocean_mean_m  1.798e-04  3.172e-03   0.057    0.964
    ## tidal_lag_time_minutes   -6.387e-04  3.548e-03  -0.180    0.887
    ## max_depth                -2.922e-02  2.261e-02  -1.292    0.419
    ## logArea                   4.152e-02  1.407e-01   0.295    0.817
    ## 
    ## Residual standard error: 0.2044 on 1 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.9223, Adjusted R-squared:  0.534 
    ## F-statistic: 2.375 on 5 and 1 DF,  p-value: 0.455

``` r
p_values <- summary(FR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(FR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                          Df  Sum Sq Mean Sq F value Pr(>F)
    ## volume_m3_w_chemocline    1 0.39843 0.39843  9.5388 0.1993
    ## distance_to_ocean_mean_m  1 0.00218 0.00218  0.0523 0.8569
    ## tidal_lag_time_minutes    1 0.02084 0.02084  0.4990 0.6085
    ## max_depth                 1 0.07091 0.07091  1.6976 0.4167
    ## logArea                   1 0.00364 0.00364  0.0871 0.8173
    ## Residuals                 1 0.04177 0.04177

``` r
# Get R-squared value
r_squared <- summary(FR_model)$r.squared

# Get number of observations
n <- nrow(traits_sespd_env[c(3,6,8,10,13:15,23),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(FR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 4.749898

``` r
FR_model <- lm(pd.obs ~ temperature_median + salinity_median + oxygen_median + pH_median, data = traits_sespd_env[c(3,6,8,10,13:15,23),])
summary(FR_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = traits_sespd_env[c(3, 6, 
    ##     8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##         3         6         8        10        13        14        15        23 
    ## -0.171803  0.149576  0.005746 -0.344991  0.421798  0.065097 -0.372012  0.246587 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -37.176048  21.141060  -1.758   0.1769  
    ## temperature_median  -0.359475   0.292248  -1.230   0.3063  
    ## salinity_median      0.002123   0.504160   0.004   0.9969  
    ## oxygen_median       -1.122026   0.569559  -1.970   0.1434  
    ## pH_median            7.158697   1.832729   3.906   0.0298 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4291 on 3 degrees of freedom
    ## Multiple R-squared:  0.9009, Adjusted R-squared:  0.7688 
    ## F-statistic: 6.819 on 4 and 3 DF,  p-value: 0.07334

``` r
p_values <- summary(FR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(FR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##                    Df  Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median  1 0.85265 0.85265  4.6315 0.1205  
    ## salinity_median     1 0.99925 0.99925  5.4279 0.1022  
    ## oxygen_median       1 0.36100 0.36100  1.9609 0.2559  
    ## pH_median           1 2.80878 2.80878 15.2571 0.0298 *
    ## Residuals           3 0.55229 0.18410                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get R-squared value
r_squared <- summary(FR_model)$r.squared

# Get number of observations
n <- nrow(traits_sespd_env[c(3,6,8,10,13:15,23),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(FR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 6.819347

``` r
FR_model <- lm(pd.obs ~ pH_median, data = traits_sespd_env[c(3,6,8,10,13:15,23),])
summary(FR_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ pH_median, data = traits_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7775 -0.3734  0.2283  0.3601  0.4659 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  -34.125      9.700  -3.518  0.01255 * 
    ## pH_median      4.717      1.238   3.810  0.00887 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5213 on 6 degrees of freedom
    ## Multiple R-squared:  0.7075, Adjusted R-squared:  0.6588 
    ## F-statistic: 14.51 on 1 and 6 DF,  p-value: 0.008865

``` r
anova(FR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## pH_median  1 3.9438  3.9438  14.515 0.008865 **
    ## Residuals  6 1.6302  0.2717                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(traits_sespd_env[c(3,6,8,10,13:15,23),3], traits_sespd_env[c(3,6,8,10,13:15,23),18])
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Island biogeography common relationships
fr_model <- lm(log(pd.obs) ~ logArea, data = traits_sespd_env[c(3,6,8,10,13:15,23),])
summary(fr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ logArea, data = traits_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.41701 -0.02270  0.08463  0.12394  0.15847 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept) -0.68082    0.54487  -1.250    0.258  
    ## logArea      0.17331    0.05585   3.103    0.021 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2191 on 6 degrees of freedom
    ## Multiple R-squared:  0.6161, Adjusted R-squared:  0.5521 
    ## F-statistic: 9.629 on 1 and 6 DF,  p-value: 0.02103

``` r
anova(fr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##           Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## logArea    1 0.46214 0.46214   9.629 0.02103 *
    ## Residuals  6 0.28797 0.04799                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
frz_model <- lm(pd.obs.z ~ logArea, data = traits_sespd_env[c(3,6,8,10,13:15,23),])
summary(frz_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ logArea, data = traits_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1524 -0.2958 -0.1765  0.7527  1.6566 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   5.3455     3.0304   1.764   0.1282  
    ## logArea      -0.7657     0.3106  -2.465   0.0488 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.218 on 6 degrees of freedom
    ## Multiple R-squared:  0.5031, Adjusted R-squared:  0.4203 
    ## F-statistic: 6.076 on 1 and 6 DF,  p-value: 0.04879

``` r
anova(frz_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## logArea    1 9.0201  9.0201  6.0758 0.04879 *
    ## Residuals  6 8.9076  1.4846                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
fr_model <- lm(log(pd.obs) ~ log(distance_to_ocean_mean_m), data = traits_sespd_env[c(3,6,8,10,13:15,23),])
summary(fr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ log(distance_to_ocean_mean_m), data = traits_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.54961 -0.20769  0.09709  0.19735  0.35720 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                    -0.3901     1.5297  -0.255    0.807
    ## log(distance_to_ocean_mean_m)   0.2871     0.3167   0.907    0.400
    ## 
    ## Residual standard error: 0.3316 on 6 degrees of freedom
    ## Multiple R-squared:  0.1205, Adjusted R-squared:  -0.02608 
    ## F-statistic: 0.8221 on 1 and 6 DF,  p-value: 0.3995

``` r
anova(fr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                               Df  Sum Sq  Mean Sq F value Pr(>F)
    ## log(distance_to_ocean_mean_m)  1 0.09039 0.090389  0.8221 0.3995
    ## Residuals                      6 0.65972 0.109953

``` r
fr_model <- lm(log(pd.obs) ~ log(max_depth), data = traits_sespd_env[c(3,6,8,10,13:15,23),])
summary(fr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ log(max_depth), data = traits_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.43030 -0.13217  0.04773  0.17077  0.33963 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)      0.4321     0.3482   1.241    0.261
    ## log(max_depth)   0.2352     0.1395   1.685    0.143
    ## 
    ## Residual standard error: 0.2913 on 6 degrees of freedom
    ## Multiple R-squared:  0.3213, Adjusted R-squared:  0.2082 
    ## F-statistic:  2.84 on 1 and 6 DF,  p-value: 0.1429

``` r
anova(fr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                Df  Sum Sq  Mean Sq F value Pr(>F)
    ## log(max_depth)  1 0.24099 0.240995  2.8402 0.1429
    ## Residuals       6 0.50911 0.084852

## Plot functional richness for all sites with more than 1 species

``` r
fr_plot <- ggplot(data = traits_sespd_env[-20,], mapping = aes(y = log(pd.obs), x = logArea, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = traits_sespd_env[-20,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Log Area (m"^"2"~")", y="FRic", colour = "Site type", fill = "Site type")
fr_plot <- fr_plot + guides(color = guide_legend(override.aes = list(label = "")))
fr_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/fr_plot_area.png", plot = fr_plot, width = 6, height = 4, units = "in")

fr_plot <- ggplot(data = traits_sespd_env[-20,], mapping = aes(y = log(pd.obs), x = log(distance_to_ocean_mean_m), color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = traits_sespd_env[-20,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Log Isolation", y="FRic", colour = "Site type", fill = "Site type")
fr_plot <- fr_plot + guides(color = guide_legend(override.aes = list(label = "")))
fr_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/fr_plot_dist.png", plot = fr_plot, width = 6, height = 4, units = "in")

fr_plot <- ggplot(data = traits_sespd_env[-20,], mapping = aes(y = log(pd.obs), x = log(max_depth), color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = traits_sespd_env[-20,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Log Age", y="FRic", colour = "Site type", fill = "Site type")
fr_plot <- fr_plot + guides(color = guide_legend(override.aes = list(label = "")))
fr_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/fr_plot_depth.png", plot = fr_plot, width = 6, height = 4, units = "in")
```

## Traits Z score vs p value plots

``` r
traits_sesmpd_plot <- ggplot(data = traits_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Stratification, fill = Stratification, label = X)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = traits_sesmpd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x= "z-score", y= "p-value", colour = "Site type", fill = "Site type", tag = "A") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sesmpd_plot <- traits_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sesmpd_plot
```

![](trait_div_analyses_files/figure-gfm/FD%20z%20score%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sesmpd_plot_z.png", traits_sesmpd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 5 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_traits_sesmpd_env <- subset(traits_sesmpd_env, Stratification == "Ocean")
ordered_subset_traits_sesmpd_env <- subset_traits_sesmpd_env[order(subset_traits_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_traits_sesmpd_env[2,7]
Q3 <- ordered_subset_traits_sesmpd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -6.731934

``` r
Q3 + 1.5*IQR1
```

    ## [1] -1.023946

``` r
ordered_subset_traits_sesmpd_env$mpd.obs.z
```

    ## [1] -6.90169201 -4.59143848 -4.24911374 -3.17193646 -3.16444130 -0.06529561

``` r
subset_traits_sesmpd_env <- subset(traits_sesmpd_env, Stratification == "Mixed")
ordered_subset_traits_sesmpd_env <- subset_traits_sesmpd_env[order(subset_traits_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_traits_sesmpd_env[3,7]
Q3 <- ordered_subset_traits_sesmpd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -5.145994

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.98926

``` r
ordered_subset_traits_sesmpd_env$mpd.obs.z
```

    ## [1] -5.2091747 -3.5014858 -2.4702739 -2.1486624 -1.8448545 -0.6864604 -0.4115270
    ## [8]  0.4941511

``` r
subset_traits_sesmpd_env <- subset(traits_sesmpd_env, Stratification == "Stratified")
ordered_subset_traits_sesmpd_env <- subset_traits_sesmpd_env[order(subset_traits_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_traits_sesmpd_env[3,7]
Q3 <- ordered_subset_traits_sesmpd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -2.244262

``` r
Q3 + 1.5*IQR1
```

    ## [1] 2.582766

``` r
ordered_subset_traits_sesmpd_env$mpd.obs.z
```

    ## [1] -0.9672644 -0.7940237 -0.4341267 -0.3979434  0.1640887  0.7726302  0.8829111
    ## [8]  0.9743155

``` r
traits_sesmntd_plot <- ggplot(data = traits_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = traits_sesmntd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x= "z-score", y= "p-value", colour = "Site type", fill = "Site type", tag = "B") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sesmntd_plot <- traits_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sesmntd_plot
```

    ## Warning: ggrepel: 7 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](trait_div_analyses_files/figure-gfm/FD%20z%20score%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sesmntd_plot_z.png", traits_sesmntd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_traits_sesmntd_env <- subset(traits_sesmntd_env, Stratification == "Ocean")
ordered_subset_traits_sesmntd_env <- subset_traits_sesmntd_env[order(subset_traits_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_traits_sesmntd_env[2,7]
Q3 <- ordered_subset_traits_sesmntd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -4.315246

``` r
Q3 + 1.5*IQR1
```

    ## [1] -1.22859

``` r
ordered_subset_traits_sesmntd_env$mntd.obs.z
```

    ## [1] -4.5637766 -3.1577499 -3.0996684 -2.4421858 -2.3860860 -0.6933401

``` r
subset_traits_sesmntd_env <- subset(traits_sesmntd_env, Stratification == "Mixed")
ordered_subset_traits_sesmntd_env <- subset_traits_sesmntd_env[order(subset_traits_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_traits_sesmntd_env[3,7]
Q3 <- ordered_subset_traits_sesmntd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -4.696753

``` r
Q3 + 1.5*IQR1
```

    ## [1] 0.7220788

``` r
ordered_subset_traits_sesmntd_env$mntd.obs.z
```

    ## [1] -4.0630971 -3.7344438 -2.6646909 -2.2821660 -1.8326524 -1.3099831 -0.7445251
    ## [8] -0.3851049

``` r
subset_traits_sesmntd_env <- subset(traits_sesmntd_env, Stratification == "Stratified")
ordered_subset_traits_sesmntd_env <- subset_traits_sesmntd_env[order(subset_traits_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_traits_sesmntd_env[3,7]
Q3 <- ordered_subset_traits_sesmntd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -1.228819

``` r
Q3 + 1.5*IQR1
```

    ## [1] 0.3183194

``` r
ordered_subset_traits_sesmntd_env$mntd.obs.z
```

    ## [1] -1.0914083 -0.9043121 -0.6486423 -0.3219626 -0.3096903 -0.2618576  0.1908168
    ## [8]  1.6365448

``` r
traits_sespd_plot <- ggplot(data = traits_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = traits_sespd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x= "z-score", y= "p-value", colour = "Site type", fill = "Site type", tag = "C") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) + 
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sespd_plot <- traits_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sespd_plot
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_text_repel()`).

    ## Warning: ggrepel: 1 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](trait_div_analyses_files/figure-gfm/FD%20z%20score%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sespd_plot_z.png", traits_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_text_repel()`).

    ## Warning: ggrepel: 6 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_traits_sespd_env <- subset(traits_sespd_env, Stratification == "Ocean")
ordered_subset_traits_sespd_env <- subset_traits_sespd_env[order(subset_traits_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_traits_sespd_env[2,7]
Q3 <- ordered_subset_traits_sespd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -6.140314

``` r
Q3 + 1.5*IQR1
```

    ## [1] -0.3158404

``` r
ordered_subset_traits_sespd_env$pd.obs.z
```

    ## [1] -5.501649 -3.956136 -3.581414 -3.220533 -2.500018 -0.451372

``` r
subset_traits_sespd_env <- subset(traits_sespd_env, Stratification == "Mixed")
ordered_subset_traits_sespd_env <- subset_traits_sespd_env[order(subset_traits_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_traits_sespd_env[3,7]
Q3 <- ordered_subset_traits_sespd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -5.080133

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.094516

``` r
ordered_subset_traits_sespd_env$pd.obs.z
```

    ## [1] -4.42104675 -4.00730873 -2.76463966 -2.00405916 -1.63330984 -1.22097740
    ## [7] -0.24841483 -0.08647715

``` r
subset_traits_sespd_env <- subset(traits_sespd_env, Stratification == "Stratified")
ordered_subset_traits_sespd_env <- subset_traits_sespd_env[order(subset_traits_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_traits_sespd_env[3,7]
Q3 <- ordered_subset_traits_sespd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -1.376311

``` r
Q3 + 1.5*IQR1
```

    ## [1] 0.8000252

``` r
ordered_subset_traits_sespd_env$pd.obs.z
```

    ## [1] -1.19943142 -0.93204586 -0.56018505 -0.50721662 -0.48323204 -0.01610094
    ## [7]  0.89133996  1.11034836

## Traits Z score vs obs value plots

``` r
traits_sesmpd_plot <- ggplot(data = traits_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed MPD", x= "z-score", colour = "Site type", fill = "Site type") +
  # ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sesmpd_plot <- traits_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sesmpd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sesmpd_plot_obs+z.png", traits_sesmpd_plot, width = 6, height = 4, units = "in")

traits_sesmntd_plot <- ggplot(data = traits_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed MNTD", x= "z-score", colour = "Site type", fill = "Site type") +
  # ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sesmntd_plot <- traits_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sesmntd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sesmntd_plot_obs+z.png", traits_sesmntd_plot, width = 6, height = 4, units = "in")

traits_sespd_plot <- ggplot(data = traits_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed FD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,6)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=6, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sespd_plot <- traits_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sespd_plot
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sespd_plot_obs+z.png", traits_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

## Straits Z score vs p value plots

``` r
straits_sesmpd_plot <- ggplot(data = straits_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = straits_sesmpd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "D") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmpd_plot <- straits_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmpd_plot
```

    ## Warning: ggrepel: 1 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sesmpd_plot_z.png", straits_sesmpd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 6 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_straits_sesmpd_env <- subset(straits_sesmpd_env, Stratification == "Ocean")
ordered_subset_straits_sesmpd_env <- subset_straits_sesmpd_env[order(subset_straits_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_straits_sesmpd_env[2,7]
Q3 <- ordered_subset_straits_sesmpd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -5.417994

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.250204

``` r
ordered_subset_straits_sesmpd_env$mpd.obs.z
```

    ## [1] -4.751094 -2.917420 -2.317507 -1.928690 -1.250370  1.252474

``` r
subset_straits_sesmpd_env <- subset(straits_sesmpd_env, Stratification == "Mixed")
ordered_subset_straits_sesmpd_env <- subset_straits_sesmpd_env[order(subset_straits_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_straits_sesmpd_env[3,7]
Q3 <- ordered_subset_straits_sesmpd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -2.643419

``` r
Q3 + 1.5*IQR1
```

    ## [1] 3.237917

``` r
ordered_subset_straits_sesmpd_env$mpd.obs.z
```

    ## [1] -2.6243314 -1.8281103 -0.4379180  0.0691394  0.9701496  1.0324162  1.2778028
    ## [8]  1.6597409

``` r
subset_straits_sesmpd_env <- subset(straits_sesmpd_env, Stratification == "Stratified")
ordered_subset_straits_sesmpd_env <- subset_straits_sesmpd_env[order(subset_straits_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_straits_sesmpd_env[3,7]
Q3 <- ordered_subset_straits_sesmpd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -0.2636995

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.795727

``` r
ordered_subset_straits_sesmpd_env$mpd.obs.z
```

    ## [1] -0.8125158 -0.7277003  0.5085854  0.5890844  0.7748363  1.0234420  1.1236771
    ## [8]  1.9051685

``` r
straits_sesmntd_plot <- ggplot(data = straits_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = straits_sesmntd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "E") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmntd_plot <- straits_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmntd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sesmntd_plot_z.png", straits_sesmntd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 5 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_straits_sesmntd_env <- subset(straits_sesmntd_env, Stratification == "Ocean")
ordered_subset_straits_sesmntd_env <- subset_straits_sesmntd_env[order(subset_straits_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_straits_sesmntd_env[2,7]
Q3 <- ordered_subset_straits_sesmntd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -4.542634

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.318892

``` r
ordered_subset_straits_sesmntd_env$mntd.obs.z
```

    ## [1] -3.1408629 -2.3445619 -1.4553954 -1.3413603 -0.8791803  0.4600560

``` r
subset_straits_sesmntd_env <- subset(straits_sesmntd_env, Stratification == "Mixed")
ordered_subset_straits_sesmntd_env <- subset_straits_sesmntd_env[order(subset_straits_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_straits_sesmntd_env[3,7]
Q3 <- ordered_subset_straits_sesmntd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -3.842263

``` r
Q3 + 1.5*IQR1
```

    ## [1] 2.961594

``` r
ordered_subset_straits_sesmntd_env$mntd.obs.z
```

    ## [1] -3.2827666 -1.6189741 -1.2908162 -0.4011833 -0.2256625  0.4101480  0.7688688
    ## [8]  0.8915420

``` r
subset_straits_sesmntd_env <- subset(straits_sesmntd_env, Stratification == "Stratified")
ordered_subset_straits_sesmntd_env <- subset_straits_sesmntd_env[order(subset_straits_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_straits_sesmntd_env[3,7]
Q3 <- ordered_subset_straits_sesmntd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -1.733817

``` r
Q3 + 1.5*IQR1
```

    ## [1] 2.23897

``` r
ordered_subset_straits_sesmntd_env$mntd.obs.z
```

    ## [1] -1.10300288 -0.76681746 -0.24402176  0.01115308  0.47165515  0.74917478
    ## [7]  1.07169951  1.73092352

``` r
straits_sespd_plot <- ggplot(data = straits_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = straits_sespd_env[,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "F") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sespd_plot <- straits_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sespd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sespd_plot_z.png", straits_sespd_plot, width = 6, height = 4, units = "in")

subset_straits_sespd_env <- subset(straits_sespd_env, Stratification == "Ocean")
ordered_subset_straits_sespd_env <- subset_straits_sespd_env[order(subset_straits_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_straits_sespd_env[2,7]
Q3 <- ordered_subset_straits_sespd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -5.123513

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.537295

``` r
ordered_subset_straits_sespd_env$pd.obs.z
```

    ## [1] -4.2146760 -2.6257099 -2.2496818 -1.0530260 -0.9605080  0.3023603

``` r
subset_straits_sespd_env <- subset(straits_sespd_env, Stratification == "Mixed")
ordered_subset_straits_sespd_env <- subset_straits_sespd_env[order(subset_straits_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_straits_sespd_env[3,7]
Q3 <- ordered_subset_straits_sespd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -3.234449

``` r
Q3 + 1.5*IQR1
```

    ## [1] 2.978314

``` r
ordered_subset_straits_sespd_env$pd.obs.z
```

    ## [1] -2.94642728 -1.64194190 -0.90466291 -0.09250091  0.04083578  0.64852803
    ## [7]  1.47529669  1.57400278

``` r
subset_straits_sespd_env <- subset(straits_sespd_env, Stratification == "Stratified")
ordered_subset_straits_sespd_env <- subset_straits_sespd_env[order(subset_straits_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_straits_sespd_env[3,7]
Q3 <- ordered_subset_straits_sespd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -2.045444

``` r
Q3 + 1.5*IQR1
```

    ## [1] 2.354706

``` r
ordered_subset_straits_sespd_env$pd.obs.z
```

    ## [1] -1.5377502 -1.0879618 -0.3953878  0.1421412  0.3918870  0.7046498  1.0192597
    ## [8]  1.0995074

## Straits Z score vs obs value plots

``` r
straits_sesmpd_plot <- ggplot(data = straits_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed MPD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmpd_plot <- straits_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmpd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sesmpd_plot_obs+z.png", straits_sesmpd_plot, width = 6, height = 4, units = "in")

straits_sesmntd_plot <- ggplot(data = straits_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed MNTD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmntd_plot <- straits_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmntd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sesmntd_plot_obs+z.png", straits_sesmntd_plot, width = 6, height = 4, units = "in")

straits_sespd_plot <- ggplot(data = straits_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22),legend.title =  element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed FD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,6)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=6, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sespd_plot <- straits_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sespd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sespd_plot_obs+z.png", straits_sespd_plot, width = 6, height = 4, units = "in")
```

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
    ##  [1] ggrepel_0.9.4     viridis_0.6.4     viridisLite_0.4.2 ggplot2_3.5.1    
    ##  [5] picante_1.8.2     nlme_3.1-164      vegan_2.6-4       lattice_0.22-5   
    ##  [9] permute_0.9-7     tidyr_1.3.0       phytools_2.0-3    maps_3.4.1.1     
    ## [13] ape_5.7-1         reshape2_1.4.4    stringr_1.5.1     dplyr_1.1.4      
    ## [17] knitr_1.45       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fastmatch_1.1-4         gtable_0.3.4            xfun_0.41              
    ##  [4] numDeriv_2016.8-1.1     quadprog_1.5-8          vctrs_0.6.5            
    ##  [7] tools_4.3.1             generics_0.1.3          tibble_3.2.1           
    ## [10] fansi_1.0.6             highr_0.10              cluster_2.1.6          
    ## [13] pkgconfig_2.0.3         Matrix_1.6-4            scatterplot3d_0.3-44   
    ## [16] lifecycle_1.0.4         farver_2.1.1            compiler_4.3.1         
    ## [19] textshaping_0.3.7       munsell_0.5.0           mnormt_2.1.1           
    ## [22] combinat_0.0-8          codetools_0.2-19        htmltools_0.5.7        
    ## [25] yaml_2.3.7              pillar_1.9.0            MASS_7.3-60            
    ## [28] clusterGeneration_1.3.8 iterators_1.0.14        foreach_1.5.2          
    ## [31] phangorn_2.11.1         tidyselect_1.2.0        digest_0.6.33          
    ## [34] stringi_1.8.2           purrr_1.0.2             labeling_0.4.3         
    ## [37] splines_4.3.1           fastmap_1.1.1           grid_4.3.1             
    ## [40] colorspace_2.1-0        expm_0.999-8            cli_3.6.1              
    ## [43] magrittr_2.0.3          optimParallel_1.0-2     utf8_1.2.4             
    ## [46] withr_2.5.2             scales_1.3.0            rmarkdown_2.25         
    ## [49] igraph_1.5.1            gridExtra_2.3           ragg_1.2.6             
    ## [52] coda_0.19-4             evaluate_0.23           doParallel_1.0.17      
    ## [55] mgcv_1.9-0              rlang_1.1.2             Rcpp_1.0.11            
    ## [58] glue_1.6.2              rstudioapi_0.15.0       R6_2.5.1               
    ## [61] plyr_1.8.9              systemfonts_1.0.5

## Regression of MPD

``` r
#MPD
man_traits_sesmpd_env_chem <- manova(cbind(temperature_median, salinity_median, oxygen_median, pH_median) ~ mpd.obs.z, data = traits_sesmpd_env)
summary.aov(man_traits_mpd_env_chem)

# traits_sesmpd_env_chem_model <- glmulti("mpd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = traits_sesmpd_env, method = "h")

man_traits_sesmpd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth) ~ mpd.obs.z, data = traits_sesmpd_env)
summary.aov(man_traits_mpd_env_phys)

# traits_sesmpd_env_phys_model <- glmulti("mpd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = traits_sesmpd_env, method = "h")
# 
# traits_sesmpd_env_model <- glmulti("mpd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "tidal_efficiency", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = traits_sesmpd_env, method = "h")

traits_sesmpd_env_model <- lm(mpd.obs.z ~ -1 + distance_to_ocean_min_m + tidal_efficiency + salinity_median + oxygen_median, data = traits_sesmpd_env)
summary(traits_mpd_env_model)
anova(traits_mpd_env_model)
```

## Plot MPD vs Oxygen Concentration

``` r
traits_sesmpd_oxygen_plot <- ggplot(data = traits_sesmpd_env, mapping = aes(x = oxygen_median, y = mpd.obs.z, color = Stratification, fill = Stratification, label = X)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "sesMPD") +
  ylim(c(-6,3))
traits_sesmpd_oxygen_plot <- traits_sesmpd_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = mpd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sesmpd_oxygen_plot.png", traits_sesmpd_oxygen_plot, width = 2, height = 2)
```

## Plot MPD vs Tidal Efficiency

``` r
traits_sesmpd_oxygen_plot <- ggplot(data = traits_sesmpd_env, mapping = aes(x = oxygen_median, y = mpd.obs.z, color = Stratification, fill = Stratification, label = X)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "sesMPD") +
  ylim(c(-6,3))
traits_sesmpd_oxygen_plot <- traits_sesmpd_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = mpd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sesmpd_oxygen_plot.png", traits_sesmpd_oxygen_plot, width = 2, height = 2)
```

## Regression of MNTD

``` r
#MNTD
man_traits_sesmntd_env_chem <- manova(cbind(temperature_median, conductivity_median, salinity_median, oxygen_median, pH_median) ~ mntd.obs.z, data = traits_sesmntd_env)
summary.aov(man_traits_mntd_env_chem)

man_traits_sesmntd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth) ~ mntd.obs.z, data = traits_sesmntd_env)
summary.aov(man_traits_mntd_env_phys)

traits_sesmntd_env_phys_model <- glmulti("mntd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = traits_sesmntd_env, method = "h")

traits_sesmntd_env_chem_model <- glmulti("mntd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = traits_sesmntd_env, method = "h")

traits_sesmntd_env_model <- glmulti("mntd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "tidal_efficiency", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = traits_sesmntd_env, method = "h")

traits_sesmntd_env_model <- lm(mntd.obs.z ~ -1 + distance_to_ocean_min_m + salinity_median + oxygen_median, data = traits_sesmntd_env)
summary(traits_mntd_env_model)
anova(traits_mntd_env_model)
```

## Plot MNTD vs Oxygen Concentration

``` r
traits_sesmntd_oxygen_plot <- ggplot(data = traits_sesmntd_env, mapping = aes(x = oxygen_median, y = mntd.obs.z, color = Stratification, fill = Stratification, label = X)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "sesMNTD") +
  ylim(c(-6,2))
traits_sesmntd_oxygen_plot <- traits_sesmntd_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = mntd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sesmntd_oxygen_plot.png", traits_sesmntd_oxygen_plot, width = 2, height = 2)
```

## Plot MNTD vs Tidal Efficiency

``` r
traits_sesmntd_oxygen_plot <- ggplot(data = traits_sesmntd_env, mapping = aes(x = oxygen_median, y = mntd.obs.z, color = Stratification, fill = Stratification, label = X)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "sesMNTD") +
  ylim(c(-6,2))
traits_sesmntd_oxygen_plot <- traits_sesmntd_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = mntd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sesmntd_oxygen_plot.png", traits_sesmntd_oxygen_plot, width = 2, height = 2)
```

## Regression of PD

``` r
#PD
man_traits_sespd_env_chem <- manova(cbind(temperature_median, salinity_median, oxygen_median, pH_median) ~ pd.obs.z, data = traits_sespd_env)
summary.aov(man_traits_sespd_env_chem)

lm_traits_sespd_env_chem <- lm(pd.obs.z ~ oxygen_median, data = traits_sespd_env)
summary(lm_traits_sespd_env_chem)

traits_sespd_env_chem_model <- glmulti("pd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = traits_sespd_env, method = "h")

man_traits_sespd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth, logArea) ~ pd.obs.z, data = traits_sespd_env)
summary.aov(man_traits_sespd_env_phys)

lm_traits_sespd_tidal <- lm(pd.obs.z ~ tidal_efficiency, data = traits_sespd_env)
summary(lm_traits_sespd_tidal)

traits_sespd_env_phys_model <- glmulti("pd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = traits_sespd_env, method = "h")

traits_sespd_env_model <- lm(pd.obs.z ~ -1 + distance_to_ocean_min_m + tidal_efficiency + oxygen_median + pH_median, data = traits_sespd_env)
summary(traits_pd_env_model)
anova(traits_pd_env_model)
```

## Plot PD vs Oxygen Concentration

``` r
traits_sespd_env_plot <- ggplot(data = traits_sespd_env, mapping = aes(x = oxygen_median, y = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 3,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.5, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 18),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Env Variable", y= "sesFPD")
traits_sespd_env_plot + geom_smooth(mapping = aes(x = env_variable, y = pd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sespd_env_plot.png", traits_sespd_env_plot, width = 6, height = 4, units = "in")
```

## Plot PD vs Tidal Efficiency

``` r
traits_sespd_env_plot <- ggplot(data = traits_sespd_env, mapping = aes(x = oxygen_median, y = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 3,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.5, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 18),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Env Variable", y= "sesFPD")
traits_sespd_env_plot + geom_smooth(mapping = aes(x = env_variable, y = pd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sespd_env_plot.png", traits_sespd_env_plot, width = 6, height = 4, units = "in")
```

## Plot z-scores and environmental variables

``` r
traits_phylo_ses_tidal_plot <- ggplot(data =  traits_phylo_ses_env, mapping = aes(x = tidal_efficiency, y = obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.5, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 24),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(x= "Tidal Efficiency", y= "z-score") +
  ylim(c(-10,10))
traits_phylo_ses_tidal_plot <- traits_phylo_ses_tidal_plot + geom_smooth(mapping = aes(x = tidal_efficiency, y = obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
traits_phylo_ses_tidal_plot <- traits_phylo_ses_tidal_plot + facet_wrap(ncol=3)
traits_phylo_ses_tidal_plot
ggsave("traits_phylo_ses_tidal_plot.png", traits_phylo_ses_tidal_plot, width = 2, height = 2)


traits_phylo_ses_oxygen_plot <- ggplot(data = traits_phylo_FD_all_env[-c(65:72,121:128,137:144,153:160),], mapping = aes(x = oxygen_median, y = obs.z, color = Environment)) + 
  geom_point(stat = 'identity',
    size = 8,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, end = .95, discrete = TRUE, option = "H") +
  theme_bw() +
  theme(text = element_text(size = 24),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "z-score") +
  ylim(c(-10,10))
traits_phylo_ses_oxygen_plot <- traits_phylo_ses_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
traits_phylo_ses_oxygen_plot <- traits_phylo_ses_oxygen_plot + facet_wrap(vars(Indice), ncol=2)
traits_phylo_ses_oxygen_plot
ggsave("traits_phylo_ses_oxygen_plot.png", traits_phylo_ses_oxygen_plot, width = 2, height = 2)
```

## Beta Diversity of Traits

``` r
traits_dist <- gowdis(main_traits)
traits_dist2 <- cailliez(traits_dist)
traits_dist2_pco <- dudi.pco(traits_dist2, scannf = FALSE, full = TRUE)
traits_pco_round <- round(traits_dist2_pco$li, .Machine$double.exponent)
traits_pco_complete <- traits_pco_round[, 1:2]
qual.FRic <- sum(traits_dist2_pco$eig[1:2]) / sum(traits_dist2_pco$eig)

# Do the following to order lakes by the circulation pattern, then delete env data
####https://pedrohbraga.github.io/StratificationPhylogenetics-Workshop/StratificationPhylogenetics-Workshop.html#between-assemblage-phylogenetic-structure####
presabs_env <- merge(pres_abs_by_lake, env_byLake, by = 'X')
presabs_env <- presabs_env[order(presabs_env$Environment),]
row.names(presabs_env) <- presabs_env$X
presabs_env_ordered <- presabs_env[,-c(1, 1692:1722)]

trait_beta <- functional.betapart.core.pairwise(presabs_env_ordered[-23,], traits_dist, parallel = TRUE)


trait_betap_jac <- functional.beta.pair(presabs_env_ordered[-23,], main_traits[,-c(1,5)], index.family="jaccard")

# Turnover matrix:
trait.jac.turn <- trait_betap_jac$beta.jtu %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("sites") + ylab("sites") + labs(fill = "Turnover") + # edit axis and legend titles
  scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels

# Nestedness matrix:
trait.jac.nest <- trait_betap_jac$beta.jne %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("sites") + ylab("sites") + labs(fill = "Nestedness") + # edit axis and legend titles
    scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels

# Nestedness matrix:
trait.jac <- trait_betap_jac$beta.jac %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("sites") + ylab("sites") + labs(fill = "Beta Diversity") + # edit axis and legend titles
    scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels
trait.jac

# plot both heatmaps next to each other
gridExtra::grid.arrange(trait.jac.turn, trait.jac.nest, ncol = 2)
```

## Phylogenetic Signal

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
