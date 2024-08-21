Phylogenetic diversity analyses
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

## Phylogenetic diversity MPD, MNTD, PD for all

``` r
#Mean pairwise differences
phylo_sesmpd <- ses.mpd(presabs_lake, cophenetic(tree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(phylo_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sesmpd.csv")


#Mean nearest taxon distance
phylo_sesmntd <- ses.mntd(presabs_lake, cophenetic(tree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(phylo_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sesmntd.csv")


#Faiths PD
phylo_sespd <- ses.pd(presabs_lake, tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)

write.csv(phylo_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sespd.csv")
```

## Phylogenetic diversity MPD, MNTD, PD for sites

``` r
#Mean pairwise differences
sphylo_sesmpd <- ses.mpd(surveyed_sites_lake, cophenetic(stree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(sphylo_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sesmpd.csv")


#Mean nearest taxon distance
sphylo_sesmntd <- ses.mntd(surveyed_sites_lake, cophenetic(stree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(sphylo_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sesmntd.csv")


#Faiths PD
sphylo_sespd <- ses.pd(surveyed_sites_lake, stree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)

write.csv(sphylo_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sespd.csv")
```

## Read in PD mpd, mntd, pd

- Read in Phylogenetic Diversity files and combine with env data

``` r
phylo_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sesmpd.csv")

phylo_sesmpd_env <- merge(phylo_sesmpd, env, by = "X", sort = F)
phylo_sesmpd_env$measure <- "phylo_sesmpd"
phylo_sesmpd_env$Stratification <- factor(phylo_sesmpd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


phylo_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sesmntd.csv")

phylo_sesmntd_env <- merge(phylo_sesmntd, env, by = "X", sort = F)
phylo_sesmntd_env$measure <- "phylo_sesmntd"
phylo_sesmntd_env$Stratification <- factor(phylo_sesmntd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


phylo_sespd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sespd.csv")

phylo_sespd_env <- merge(phylo_sespd, env, by = "X", sort = F)
phylo_sespd_env$measure <- "phylo_sespd"
phylo_sespd_env$Stratification <- factor(phylo_sespd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


sphylo_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sesmpd.csv")

sphylo_sesmpd_env <- merge(sphylo_sesmpd, env[-20,], by = "X", sort = F)
sphylo_sesmpd_env$measure <- "sphylo_sesmpd"
sphylo_sesmpd_env$Stratification <- factor(sphylo_sesmpd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))


sphylo_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sesmntd.csv")

sphylo_sesmntd_env <- merge(sphylo_sesmntd, env[-20,], by = "X", sort = F)
sphylo_sesmntd_env$measure <- "sphylo_sesmntd"
sphylo_sesmntd_env$Stratification <- factor(sphylo_sesmntd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))


sphylo_sespd <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sespd.csv")

sphylo_sespd_env <- merge(sphylo_sespd, env[-20,], by = "X", sort = F)
sphylo_sespd_env$measure <- "sphylo_sespd"
sphylo_sespd_env$Stratification <- factor(sphylo_sespd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))
```

## PR complete island Linear Models

``` r
PR_model <- lm(log(pd.obs) ~ volume_m3_w_chemocline + surface_area_m2 + distance_to_ocean_mean_m + tidal_lag_time_minutes + max_depth + logArea, data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(PR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ volume_m3_w_chemocline + surface_area_m2 + 
    ##     distance_to_ocean_mean_m + tidal_lag_time_minutes + max_depth + 
    ##     logArea, data = phylo_sespd_env[c(1:2, 4, 5, 12, 17, 21:22), 
    ##     ])
    ## 
    ## Residuals:
    ##         1         2         4         5        12        17        21        22 
    ##  0.033811 -0.041902  0.018510 -0.006714  0.042853 -0.036917  0.003337 -0.012978 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)               3.796e+00  9.548e-01   3.976   0.1569  
    ## volume_m3_w_chemocline    1.170e-06  1.502e-07   7.793   0.0812 .
    ## surface_area_m2          -1.121e-05  5.326e-06  -2.104   0.2825  
    ## distance_to_ocean_mean_m -1.680e-03  5.457e-04  -3.079   0.1999  
    ## tidal_lag_time_minutes    5.481e-03  1.279e-03   4.284   0.1460  
    ## max_depth                -4.324e-02  7.105e-03  -6.086   0.1037  
    ## logArea                   3.586e-01  1.133e-01   3.165   0.1948  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08164 on 1 degrees of freedom
    ## Multiple R-squared:  0.9903, Adjusted R-squared:  0.9322 
    ## F-statistic: 17.03 on 6 and 1 DF,  p-value: 0.1834

``` r
p_values <- summary(PR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(PR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                          Df   Sum Sq  Mean Sq F value Pr(>F)
    ## volume_m3_w_chemocline    1 0.166446 0.166446 24.9715 0.1257
    ## surface_area_m2           1 0.123657 0.123657 18.5519 0.1452
    ## distance_to_ocean_mean_m  1 0.089811 0.089811 13.4741 0.1693
    ## tidal_lag_time_minutes    1 0.053402 0.053402  8.0118 0.2162
    ## max_depth                 1 0.181034 0.181034 27.1601 0.1207
    ## logArea                   1 0.066774 0.066774 10.0179 0.1948
    ## Residuals                 1 0.006665 0.006665

``` r
# Get R-squared value
r_squared <- summary(PR_model)$r.squared

# Get number of observations
n <- nrow(phylo_sespd_env[c(1:2,4,5,12,17,21:22),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(PR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 17.03121

``` r
PR_model <- lm(pd.obs ~ temperature_median + salinity_median + oxygen_median + pH_median, data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(PR_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = phylo_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##        1        2        4        5       12       17       21       22 
    ##  128.442 -115.678  201.540   -8.766 -121.643 -162.365  -71.738  150.207 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -6263.70    3818.09  -1.641    0.199
    ## temperature_median   -57.19      68.04  -0.840    0.462
    ## salinity_median       26.72      46.20   0.578    0.604
    ## oxygen_median       -110.20     169.05  -0.652    0.561
    ## pH_median           1155.73     697.57   1.657    0.196
    ## 
    ## Residual standard error: 215.6 on 3 degrees of freedom
    ## Multiple R-squared:  0.8234, Adjusted R-squared:  0.5879 
    ## F-statistic: 3.497 on 4 and 3 DF,  p-value: 0.1659

``` r
p_values <- summary(PR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(PR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##                    Df Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median  1   4551    4551  0.0979 0.7748  
    ## salinity_median     1 515614  515614 11.0931 0.0447 *
    ## oxygen_median       1   2429    2429  0.0523 0.8339  
    ## pH_median           1 127586  127586  2.7449 0.1961  
    ## Residuals           3 139442   46481                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get R-squared value
r_squared <- summary(PR_model)$r.squared

# Get number of observations
n <- nrow(phylo_sespd_env[c(1:2,4,5,12,17,21:22),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(PR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 3.497052

``` r
PR_model <- lm(pd.obs ~ salinity_median, data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(PR_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ salinity_median, data = phylo_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -281.8 -179.5   54.8  161.3  225.5 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)     -1150.75     638.87  -1.801   0.1217  
    ## salinity_median    77.78      23.10   3.366   0.0151 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 213.4 on 6 degrees of freedom
    ## Multiple R-squared:  0.6538, Adjusted R-squared:  0.5961 
    ## F-statistic: 11.33 on 1 and 6 DF,  p-value: 0.01511

``` r
anova(PR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##                 Df Sum Sq Mean Sq F value  Pr(>F)  
    ## salinity_median  1 516271  516271  11.332 0.01511 *
    ## Residuals        6 273350   45558                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(phylo_sespd_env[c(1:2,4,5,12,17,21:22),3], phylo_sespd_env[c(1:2,4,5,12,17,21:22),14])
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# Island biogeography common relationships
pr_model <- lm(log(pd.obs) ~ logArea, data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(pr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ logArea, data = phylo_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.31149 -0.23963 -0.05928  0.10755  0.49611 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  6.45326    1.09502   5.893  0.00106 **
    ## logArea      0.03838    0.10610   0.362  0.72993   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3349 on 6 degrees of freedom
    ## Multiple R-squared:  0.02135,    Adjusted R-squared:  -0.1418 
    ## F-statistic: 0.1309 on 1 and 6 DF,  p-value: 0.7299

``` r
anova(pr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## logArea    1 0.01468 0.014682  0.1309 0.7299
    ## Residuals  6 0.67311 0.112184

``` r
prz_model <- lm(pd.obs.z ~ logArea, data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(prz_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ logArea, data = phylo_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.0090 -0.4901  0.2694  0.5057  0.6676 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)  -3.9468     2.4216  -1.630    0.154
    ## logArea       0.2777     0.2346   1.183    0.281
    ## 
    ## Residual standard error: 0.7407 on 6 degrees of freedom
    ## Multiple R-squared:  0.1892, Adjusted R-squared:  0.05412 
    ## F-statistic:   1.4 on 1 and 6 DF,  p-value: 0.2814

``` r
anova(prz_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##           Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea    1 0.7684 0.76840  1.4005 0.2814
    ## Residuals  6 3.2920 0.54867

``` r
pr_model <- lm(log(pd.obs) ~ log(distance_to_ocean_mean_m), data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(pr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ log(distance_to_ocean_mean_m), data = phylo_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.37557 -0.20516 -0.06024  0.12194  0.58454 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)                     7.8635     1.7092   4.601  0.00369 **
    ## log(distance_to_ocean_mean_m)  -0.1838     0.3084  -0.596  0.57293   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.329 on 6 degrees of freedom
    ## Multiple R-squared:  0.0559, Adjusted R-squared:  -0.1014 
    ## F-statistic: 0.3553 on 1 and 6 DF,  p-value: 0.5729

``` r
anova(pr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                               Df  Sum Sq  Mean Sq F value Pr(>F)
    ## log(distance_to_ocean_mean_m)  1 0.03845 0.038448  0.3553 0.5729
    ## Residuals                      6 0.64934 0.108223

``` r
pr_model <- lm(log(pd.obs) ~ log(max_depth), data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(pr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ log(max_depth), data = phylo_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.36256 -0.21746 -0.03948  0.09788  0.53772 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     6.873279   0.776298   8.854 0.000115 ***
    ## log(max_depth) -0.008431   0.246758  -0.034 0.973851    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3385 on 6 degrees of freedom
    ## Multiple R-squared:  0.0001945,  Adjusted R-squared:  -0.1664 
    ## F-statistic: 0.001167 on 1 and 6 DF,  p-value: 0.9739

``` r
anova(pr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                Df  Sum Sq  Mean Sq F value Pr(>F)
    ## log(max_depth)  1 0.00013 0.000134  0.0012 0.9739
    ## Residuals       6 0.68765 0.114609

## PR habitat island Linear Models

``` r
PR_model <- lm(log(pd.obs) ~ volume_m3_w_chemocline + distance_to_ocean_mean_m + tidal_lag_time_minutes + max_depth + logArea, data = phylo_sespd_env[c(3,6,8,10,13:15,23),])
summary(PR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ volume_m3_w_chemocline + distance_to_ocean_mean_m + 
    ##     tidal_lag_time_minutes + max_depth + logArea, data = phylo_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##         3         6        10        13        14        15        23 
    ##  0.098067  0.004722  0.043633 -0.080041 -0.043487 -0.073022  0.050127 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)               7.222e+00  8.635e-01   8.363   0.0758 .
    ## volume_m3_w_chemocline    2.175e-06  8.007e-07   2.716   0.2246  
    ## distance_to_ocean_mean_m -2.023e-03  2.583e-03  -0.783   0.5769  
    ## tidal_lag_time_minutes    9.373e-04  2.889e-03   0.324   0.8003  
    ## max_depth                -2.432e-02  1.841e-02  -1.321   0.4124  
    ## logArea                   8.342e-02  1.145e-01   0.728   0.5993  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1664 on 1 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.9605, Adjusted R-squared:  0.7631 
    ## F-statistic: 4.865 on 5 and 1 DF,  p-value: 0.3307

``` r
p_values <- summary(PR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(PR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                          Df  Sum Sq Mean Sq F value Pr(>F)
    ## volume_m3_w_chemocline    1 0.55037 0.55037 19.8793 0.1405
    ## distance_to_ocean_mean_m  1 0.00008 0.00008  0.0028 0.9662
    ## tidal_lag_time_minutes    1 0.07420 0.07420  2.6799 0.3491
    ## max_depth                 1 0.03415 0.03415  1.2335 0.4667
    ## logArea                   1 0.01468 0.01468  0.5303 0.5993
    ## Residuals                 1 0.02769 0.02769

``` r
# Get R-squared value
r_squared <- summary(PR_model)$r.squared

# Get number of observations
n <- nrow(phylo_sespd_env[c(3,6,8,10,13:15,23),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(PR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 9.730346

``` r
PR_model <- lm(pd.obs ~ temperature_median + salinity_median + oxygen_median + pH_median, data = phylo_sespd_env[c(3,6,8,10,13:15,23),])
summary(PR_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = phylo_sespd_env[c(3, 6, 
    ##     8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##       3       6       8      10      13      14      15      23 
    ## -181.81  401.17  135.41 -693.53  694.26   60.03 -727.06  311.52 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -36037.3    38214.7  -0.943   0.4152  
    ## temperature_median   -458.3      528.3  -0.868   0.4495  
    ## salinity_median      -346.4      911.3  -0.380   0.7292  
    ## oxygen_median       -1114.3     1029.5  -1.082   0.3583  
    ## pH_median            8908.2     3312.9   2.689   0.0745 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 775.6 on 3 degrees of freedom
    ## Multiple R-squared:  0.8094, Adjusted R-squared:  0.5552 
    ## F-statistic: 3.184 on 4 and 3 DF,  p-value: 0.1843

``` r
p_values <- summary(PR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(PR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##                    Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## temperature_median  1 1091362 1091362  1.8143 0.27070  
    ## salinity_median     1 1086118 1086118  1.8056 0.27163  
    ## oxygen_median       1 1135172 1135172  1.8872 0.26318  
    ## pH_median           1 4349426 4349426  7.2307 0.07448 .
    ## Residuals           3 1804575  601525                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get R-squared value
r_squared <- summary(PR_model)$r.squared

# Get number of observations
n <- nrow(phylo_sespd_env[c(3,6,8,10,13:15,23),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(PR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 3.184438

``` r
PR_model <- lm(pd.obs ~ pH_median, data = phylo_sespd_env[c(3,6,8,10,13:15,23),])
summary(PR_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ pH_median, data = phylo_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -955.4 -507.1  170.7  483.4  731.2 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   -44106      13149  -3.354   0.0153 *
    ## pH_median       6042       1678   3.600   0.0114 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 706.6 on 6 degrees of freedom
    ## Multiple R-squared:  0.6836, Adjusted R-squared:  0.6308 
    ## F-statistic: 12.96 on 1 and 6 DF,  p-value: 0.01136

``` r
anova(PR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##           Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## pH_median  1 6471003 6471003  12.961 0.01136 *
    ## Residuals  6 2995650  499275                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(phylo_sespd_env[c(3,6,8,10,13:15,23),3], phylo_sespd_env[c(3,6,8,10,13:15,23),18])
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Island biogeography common relationships
pr_model <- lm(log(pd.obs) ~ logArea, data = phylo_sespd_env[c(3,6,8,10,13:15,23),])
summary(pr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ logArea, data = phylo_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.45385 -0.00715  0.07400  0.12000  0.14339 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  5.92809    0.53631  11.053 3.26e-05 ***
    ## logArea      0.21641    0.05497   3.937  0.00765 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2156 on 6 degrees of freedom
    ## Multiple R-squared:  0.7209, Adjusted R-squared:  0.6744 
    ## F-statistic:  15.5 on 1 and 6 DF,  p-value: 0.007654

``` r
anova(pr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##           Df  Sum Sq Mean Sq F value   Pr(>F)   
    ## logArea    1 0.72064 0.72064  15.498 0.007654 **
    ## Residuals  6 0.27900 0.04650                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
prz_model <- lm(pd.obs.z ~ logArea, data = phylo_sespd_env[c(3,6,8,10,13:15,23),])
summary(prz_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs.z ~ logArea, data = phylo_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.1360 -0.6441 -0.2743  0.2703  1.7740 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   4.0920     2.7879   1.468   0.1925  
    ## logArea      -0.7389     0.2858  -2.586   0.0415 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.121 on 6 degrees of freedom
    ## Multiple R-squared:  0.527,  Adjusted R-squared:  0.4482 
    ## F-statistic: 6.686 on 1 and 6 DF,  p-value: 0.04145

``` r
anova(prz_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs.z
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## logArea    1 8.4007  8.4007  6.6858 0.04145 *
    ## Residuals  6 7.5390  1.2565                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
pr_model <- lm(log(pd.obs) ~ log(distance_to_ocean_mean_m), data = phylo_sespd_env[c(3,6,8,10,13:15,23),])
summary(pr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ log(distance_to_ocean_mean_m), data = phylo_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5907 -0.3123  0.1144  0.2656  0.4400 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)                     6.6877     1.8025    3.71  0.00997 **
    ## log(distance_to_ocean_mean_m)   0.2762     0.3732    0.74  0.48714   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3907 on 6 degrees of freedom
    ## Multiple R-squared:  0.08367,    Adjusted R-squared:  -0.06905 
    ## F-statistic: 0.5478 on 1 and 6 DF,  p-value: 0.4871

``` r
anova(pr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                               Df  Sum Sq  Mean Sq F value Pr(>F)
    ## log(distance_to_ocean_mean_m)  1 0.08364 0.083638  0.5478 0.4871
    ## Residuals                      6 0.91600 0.152666

``` r
pr_model <- lm(log(pd.obs) ~ log(max_depth), data = phylo_sespd_env[c(3,6,8,10,13:15,23),])
summary(pr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ log(max_depth), data = phylo_sespd_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.46962 -0.13917 -0.01723  0.22159  0.39293 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      7.2935     0.3772  19.337 1.24e-06 ***
    ## log(max_depth)   0.3039     0.1511   2.011   0.0911 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3155 on 6 degrees of freedom
    ## Multiple R-squared:  0.4026, Adjusted R-squared:  0.303 
    ## F-statistic: 4.043 on 1 and 6 DF,  p-value: 0.09107

``` r
anova(pr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## log(max_depth)  1 0.40241 0.40241  4.0428 0.09107 .
    ## Residuals       6 0.59722 0.09954                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Plot phylogenetic richness for all sites with more than 1 species

``` r
pr_plot <- ggplot(data = phylo_sespd_env[-20,], mapping = aes(y = log(pd.obs), x = logArea, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = phylo_sespd_env[-20,1], size = 5, point.padding = 3) +
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
  labs(x="Log Area (m"^"2"~")", y="Log PRic", colour = "Site type", fill = "Site type")
pr_plot <- pr_plot + guides(color = guide_legend(override.aes = list(label = "")))
pr_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/pr_plot_area.png", plot = pr_plot, width = 6, height = 4, units = "in")

pr_plot <- ggplot(data = phylo_sespd_env[-20,], mapping = aes(y = log(pd.obs), x = log(distance_to_ocean_mean_m), color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = phylo_sespd_env[-20,1], size = 5, point.padding = 3) +
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
  labs(x="Log Isolation", y="Log PRic", colour = "Site type", fill = "Site type")
pr_plot <- pr_plot + guides(color = guide_legend(override.aes = list(label = "")))
pr_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/pr_plot_dist.png", plot = pr_plot, width = 6, height = 4, units = "in")

pr_plot <- ggplot(data = phylo_sespd_env[-20,], mapping = aes(y = log(pd.obs), x = log(max_depth), color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = phylo_sespd_env[-20,1], size = 5, point.padding = 3) +
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
  labs(x="Log Age", y="Log PRic", colour = "Site type", fill = "Site type")
pr_plot <- pr_plot + guides(color = guide_legend(override.aes = list(label = "")))
pr_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/pr_plot_depth.png", plot = pr_plot, width = 6, height = 4, units = "in")
```

## Phylo Z score vs p value plots

``` r
phylo_sesmpd_plot <- ggplot(data = phylo_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = phylo_sesmpd_env[,1], size = 5, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "A") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sesmpd_plot <- phylo_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sesmpd_plot
```

    ## Warning: ggrepel: 14 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/create%20PD%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sesmpd_plot_z.png", phylo_sesmpd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 17 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_phylo_sesmpd_env <- subset(phylo_sesmpd_env, Stratification == "Ocean")
ordered_subset_phylo_sesmpd_env <- subset_phylo_sesmpd_env[order(subset_phylo_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_phylo_sesmpd_env[2,7]
Q3 <- ordered_subset_phylo_sesmpd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -3.218805

``` r
Q3 + 1.5*IQR1
```

    ## [1] -0.5399805

``` r
ordered_subset_phylo_sesmpd_env$mpd.obs.z
```

    ## [1] -3.1459677 -2.2142459 -2.0637407 -1.8389517 -1.5445397 -0.9156388

``` r
subset_phylo_sesmpd_env <- subset(phylo_sesmpd_env, Stratification == "Mixed")
ordered_subset_phylo_sesmpd_env <- subset_phylo_sesmpd_env[order(subset_phylo_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_phylo_sesmpd_env[3,7]
Q3 <- ordered_subset_phylo_sesmpd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -3.182811

``` r
Q3 + 1.5*IQR1
```

    ## [1] -0.09392653

``` r
ordered_subset_phylo_sesmpd_env$mpd.obs.z
```

    ## [1] -2.801591 -2.084557 -2.024480 -1.808298 -1.434609 -1.252258 -1.212086
    ## [8] -1.091182

``` r
subset_phylo_sesmpd_env <- subset(phylo_sesmpd_env, Stratification == "Stratified")
ordered_subset_phylo_sesmpd_env <- subset_phylo_sesmpd_env[order(subset_phylo_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_phylo_sesmpd_env[3,7]
Q3 <- ordered_subset_phylo_sesmpd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -1.867856

``` r
Q3 + 1.5*IQR1
```

    ## [1] 0.3278997

``` r
ordered_subset_phylo_sesmpd_env$mpd.obs.z
```

    ## [1] -1.3003579 -1.1387666 -1.0444477 -0.7442054 -0.6698887 -0.4955087 -0.3242263
    ## [8] -0.1934505

``` r
phylo_sesmntd_plot <- ggplot(data = phylo_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = phylo_sesmntd_env[,1], size = 5, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "B") +
  ylim(c(0,1))  +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sesmntd_plot <- phylo_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sesmntd_plot
```

    ## Warning: ggrepel: 10 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/create%20PD%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sesmntd_plot_z.png", phylo_sesmntd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 16 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_phylo_sesmntd_env <- subset(phylo_sesmntd_env, Stratification == "Ocean")
ordered_subset_phylo_sesmntd_env <- subset_phylo_sesmntd_env[order(subset_phylo_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_phylo_sesmntd_env[2,7]
Q3 <- ordered_subset_phylo_sesmntd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -8.40514

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.657627

``` r
ordered_subset_phylo_sesmntd_env$mntd.obs.z
```

    ## [1] -5.9695951 -4.6316025 -4.2064054 -2.7480364 -2.1159109 -0.6965152

``` r
subset_phylo_sesmntd_env <- subset(phylo_sesmntd_env, Stratification == "Mixed")
ordered_subset_phylo_sesmntd_env <- subset_phylo_sesmntd_env[order(subset_phylo_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_phylo_sesmntd_env[3,7]
Q3 <- ordered_subset_phylo_sesmntd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -7.467293

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.252799

``` r
ordered_subset_phylo_sesmntd_env$mntd.obs.z
```

    ## [1] -5.277404 -4.278256 -4.197259 -3.962508 -2.450007 -2.017235 -1.618733
    ## [8] -1.230735

``` r
subset_phylo_sesmntd_env <- subset(phylo_sesmntd_env, Stratification == "Stratified")
ordered_subset_phylo_sesmntd_env <- subset_phylo_sesmntd_env[order(subset_phylo_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_phylo_sesmntd_env[3,7]
Q3 <- ordered_subset_phylo_sesmntd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -3.585302

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.156256

``` r
ordered_subset_phylo_sesmntd_env$mntd.obs.z
```

    ## [1] -2.7346379 -1.9627476 -1.8072177 -1.6483935 -1.0244152 -0.6218281 -0.3051538
    ## [8] -0.1578621

``` r
phylo_sespd_plot <- ggplot(data = phylo_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = phylo_sespd_env[,1], size = 5, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "C") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sespd_plot <- phylo_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sespd_plot
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_text_repel()`).

    ## Warning: ggrepel: 11 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/create%20PD%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sespd_plot_z.png", phylo_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_text_repel()`).

    ## Warning: ggrepel: 16 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_phylo_sespd_env <- subset(phylo_sespd_env, Stratification == "Ocean")
ordered_subset_phylo_sespd_env <- subset_phylo_sespd_env[order(subset_phylo_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_phylo_sespd_env[2,7]
Q3 <- ordered_subset_phylo_sespd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -7.515035

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.044698

``` r
ordered_subset_phylo_sespd_env$pd.obs.z
```

    ## [1] -6.5822163 -4.3051349 -4.1079259 -3.1473418 -2.1652015 -0.7819302

``` r
subset_phylo_sespd_env <- subset(phylo_sespd_env, Stratification == "Mixed")
ordered_subset_phylo_sespd_env <- subset_phylo_sespd_env[order(subset_phylo_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_phylo_sespd_env[3,7]
Q3 <- ordered_subset_phylo_sespd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -6.127817

``` r
Q3 + 1.5*IQR1
```

    ## [1] 0.1346758

``` r
ordered_subset_phylo_sespd_env$pd.obs.z
```

    ## [1] -5.630572 -4.043431 -3.779382 -3.730212 -2.306180 -2.213759 -1.512303
    ## [8] -1.131156

``` r
subset_phylo_sespd_env <- subset(phylo_sespd_env, Stratification == "Stratified")
ordered_subset_phylo_sespd_env <- subset_phylo_sespd_env[order(subset_phylo_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_phylo_sespd_env[3,7]
Q3 <- ordered_subset_phylo_sespd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -2.379815

``` r
Q3 + 1.5*IQR1
```

    ## [1] 0.560136

``` r
ordered_subset_phylo_sespd_env$pd.obs.z
```

    ## [1] -2.4428886 -1.8347508 -1.2773334 -1.2608703 -0.8130289 -0.5423456 -0.3549293
    ## [8] -0.2565106

## Phylo Z score vs obs value plots

``` r
phylo_sesmpd_plot <- ggplot(data = phylo_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
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
  labs(y= "Observed MPD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,300)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=300, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sesmpd_plot <- phylo_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sesmpd_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sesmpd_plot_obs+z.png", phylo_sesmpd_plot, width = 6, height = 4, units = "in")

phylo_sesmntd_plot <- ggplot(data = phylo_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
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
  labs(y= "Observed MNTD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,250))  +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=250, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sesmntd_plot <- phylo_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sesmntd_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sesmntd_plot_obs+z.png", phylo_sesmntd_plot, width = 6, height = 4, units = "in")

phylo_sespd_plot <- ggplot(data = phylo_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
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
  labs(y= "Observed PD", x = "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,6000)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=6000, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sespd_plot <- phylo_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sespd_plot
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sespd_plot_obs+z.png", phylo_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

## Sphylo Z score vs p value plots

``` r
sphylo_sesmpd_plot <- ggplot(data = sphylo_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = sphylo_sesmpd_env[,1], size = 5, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "D") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=-0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sesmpd_plot <- sphylo_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sesmpd_plot
```

    ## Warning: ggrepel: 7 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sesmpd_plot_z.png", sphylo_sesmpd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_sphylo_sesmpd_env <- subset(sphylo_sesmpd_env, Stratification == "Ocean")
ordered_subset_sphylo_sesmpd_env <- subset_sphylo_sesmpd_env[order(subset_sphylo_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_sphylo_sesmpd_env[2,7]
Q3 <- ordered_subset_sphylo_sesmpd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -1.631156

``` r
Q3 + 1.5*IQR1
```

    ## [1] -0.3731182

``` r
ordered_subset_sphylo_sesmpd_env$mpd.obs.z
```

    ## [1] -2.7509049 -1.1593920 -1.0817916 -1.0653181 -0.8448825  0.1692666

``` r
subset_sphylo_sesmpd_env <- subset(sphylo_sesmpd_env, Stratification == "Mixed")
ordered_subset_sphylo_sesmpd_env <- subset_sphylo_sesmpd_env[order(subset_sphylo_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_sphylo_sesmpd_env[3,7]
Q3 <- ordered_subset_sphylo_sesmpd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -2.003568

``` r
Q3 + 1.5*IQR1
```

    ## [1] 0.7123557

``` r
ordered_subset_sphylo_sesmpd_env$mpd.obs.z
```

    ## [1] -1.5102017 -1.0749683 -0.9850968 -0.6006071 -0.5469298 -0.3061158 -0.2497164
    ## [8]  0.9535480

``` r
subset_sphylo_sesmpd_env <- subset(sphylo_sesmpd_env, Stratification == "Stratified")
ordered_subset_sphylo_sesmpd_env <- subset_sphylo_sesmpd_env[order(subset_sphylo_sesmpd_env$mpd.obs.z), ]
Q1 <- ordered_subset_sphylo_sesmpd_env[3,7]
Q3 <- ordered_subset_sphylo_sesmpd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -1.844118

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.027974

``` r
ordered_subset_sphylo_sesmpd_env$mpd.obs.z
```

    ## [1] -1.42961216 -0.89311964 -0.76708384 -0.52321330 -0.47901191 -0.04906082
    ## [7]  0.04385659  0.62797589

``` r
sphylo_sesmntd_plot <- ggplot(data = sphylo_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = sphylo_sesmntd_env[,1], size = 5, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Site type", fill = "Site type", tag = "E") +
  ylim(c(0,1))  +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sesmntd_plot <- sphylo_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sesmntd_plot
```

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sesmntd_plot_z.png", sphylo_sesmntd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 12 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_sphylo_sesmntd_env <- subset(sphylo_sesmntd_env, Stratification == "Ocean")
ordered_subset_sphylo_sesmntd_env <- subset_sphylo_sesmntd_env[order(subset_sphylo_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_sphylo_sesmntd_env[2,7]
Q3 <- ordered_subset_sphylo_sesmntd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -5.433269

``` r
Q3 + 1.5*IQR1
```

    ## [1] 2.497295

``` r
ordered_subset_sphylo_sesmntd_env$mntd.obs.z
```

    ## [1] -3.0070823 -2.4593073 -1.3082997 -0.8703830 -0.4766662  0.7431485

``` r
subset_sphylo_sesmntd_env <- subset(sphylo_sesmntd_env, Stratification == "Mixed")
ordered_subset_sphylo_sesmntd_env <- subset_sphylo_sesmntd_env[order(subset_sphylo_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_sphylo_sesmntd_env[3,7]
Q3 <- ordered_subset_sphylo_sesmntd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -3.370706

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.258472

``` r
ordered_subset_sphylo_sesmntd_env$mntd.obs.z
```

    ## [1] -1.8123664 -1.6998078 -1.6347643 -1.5787903 -0.4807683 -0.4774698 -0.1102438
    ## [8]  0.0744964

``` r
subset_sphylo_sesmntd_env <- subset(sphylo_sesmntd_env, Stratification == "Stratified")
ordered_subset_sphylo_sesmntd_env <- subset_sphylo_sesmntd_env[order(subset_sphylo_sesmntd_env$mntd.obs.z), ]
Q1 <- ordered_subset_sphylo_sesmntd_env[3,7]
Q3 <- ordered_subset_sphylo_sesmntd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -3.793801

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.854389

``` r
ordered_subset_sphylo_sesmntd_env$mntd.obs.z
```

    ## [1] -2.2321700 -1.7175372 -1.6757300 -1.0413421 -1.0407344 -0.2636825  0.3121675
    ## [8]  0.3514270

``` r
sphylo_sespd_plot <- ggplot(data = sphylo_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = sphylo_sespd_env[,1], size = 5, point.padding = 3) +
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
  labs(y= "p-value", x = "z-score", colour = "Site type", fill = "Site type", tag = "F") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sespd_plot <- sphylo_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sespd_plot
```

    ## Warning: ggrepel: 12 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sespd_plot_z.png", sphylo_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 13 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
subset_sphylo_sespd_env <- subset(sphylo_sespd_env, Stratification == "Ocean")
ordered_subset_sphylo_sespd_env <- subset_sphylo_sespd_env[order(subset_sphylo_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_sphylo_sespd_env[2,7]
Q3 <- ordered_subset_sphylo_sespd_env[5,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -4.291437

``` r
Q3 + 1.5*IQR1
```

    ## [1] 0.8040616

``` r
ordered_subset_sphylo_sespd_env$pd.obs.z
```

    ## [1] -3.7568835 -2.3806252 -1.3011029 -1.1177769 -1.1067505  0.8390991

``` r
subset_sphylo_sespd_env <- subset(sphylo_sespd_env, Stratification == "Mixed")
ordered_subset_sphylo_sespd_env <- subset_sphylo_sespd_env[order(subset_sphylo_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_sphylo_sespd_env[3,7]
Q3 <- ordered_subset_sphylo_sespd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -2.650386

``` r
Q3 + 1.5*IQR1
```

    ## [1] 0.8824724

``` r
ordered_subset_sphylo_sespd_env$pd.obs.z
```

    ## [1] -2.176664485 -1.505530529 -1.325564365 -1.045768699 -0.795475299
    ## [6] -0.442349652 -0.002337534  0.429460560

``` r
subset_sphylo_sespd_env <- subset(sphylo_sespd_env, Stratification == "Stratified")
ordered_subset_sphylo_sespd_env <- subset_sphylo_sespd_env[order(subset_sphylo_sespd_env$pd.obs.z), ]
Q1 <- ordered_subset_sphylo_sespd_env[3,7]
Q3 <- ordered_subset_sphylo_sespd_env[6,7]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -3.038822

``` r
Q3 + 1.5*IQR1
```

    ## [1] 1.72407

``` r
ordered_subset_sphylo_sespd_env$pd.obs.z
```

    ## [1] -2.02033203 -1.66357563 -1.25273763 -1.10365093 -0.88779549 -0.06201469
    ## [7]  0.19722676  0.23472021

## Sphylo Z score vs obs value plots

``` r
sphylo_sesmpd_plot <- ggplot(data = sphylo_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
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
  labs(y= "Observed MPD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,250)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=250, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sesmpd_plot <- sphylo_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sesmpd_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sesmpd_plot_obs+z.png", sphylo_sesmpd_plot, width = 6, height = 4, units = "in")

sphylo_sesmntd_plot <- ggplot(data = sphylo_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
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
  labs(y= "Observed MNTD", x= "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,250))  +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=250, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sesmntd_plot <- sphylo_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sesmntd_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sesmntd_plot_obs+z.png", sphylo_sesmntd_plot, width = 6, height = 4, units = "in")

sphylo_sespd_plot <- ggplot(data = sphylo_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
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
  labs(y= "Observed PD", x = "z-score", colour = "Site type", fill = "Site type") +
  ylim(c(0,6000)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=6000, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sespd_plot <- sphylo_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sespd_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sespd_plot_obs+z.png", sphylo_sespd_plot, width = 6, height = 4, units = "in")
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
phylo_sesmpd_env_chem <- phylo_sesmpd_env[-c(9, 16, 18, 20),]

man_phylo_sesmpd_env_chem <- manova(cbind(temperature_median, salinity_median, oxygen_median, pH_median) ~ mpd.obs.z, data = phylo_sesmpd_env_chem)
summary.aov(man_phylo_sesmpd_env_chem)

lm_phylo_sesmpd_env_chem <- lm(mpd.obs.z ~ temperature_median * salinity_median * oxygen_median * pH_median, data = phylo_sesmpd_env_chem)
summary(lm_phylo_sesmpd_env_chem)

step( object = lm_phylo_sesmpd_env_chem,     # start at the full model
       direction = "backward"   # allow it remove predictors but not add them
)

phylo_sesmpd_env_chem_model <- glmulti("mpd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = phylo_sesmpd_env, method = "h")
```

## Plot MPD vs Oxygen Concentration

``` r
lm_phylo_sesmpd_te <- lm(mpd.obs.z ~tidal_efficiency, data = phylo_sesmpd_env_phys)
summary(lm_phylo_sesmpd_te)

phylo_sesmpd_tidal_efficiency_plot <- ggplot(data = phylo_sesmpd_env, mapping = aes(x = tidal_efficiency, y = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Tidal Efficiency", y= "sesMPD") +
  ylim(c(-6,2))
phylo_sesmpd_tidal_efficiency_plot <- phylo_sesmpd_tidal_efficiency_plot + geom_smooth(mapping = aes(x = tidal_efficiency, y = mpd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
phylo_sesmpd_tidal_efficiency_plot
ggsave("phylo_sesmpd_tidal_efficiency_plot.png", phylo_sesmpd_tidal_efficiency_plot, width = 6, height = 4, units = "in")
```

## Plot MPD vs Tidal efficiency

``` r
lm_phylo_sesmpd_te <- lm(mpd.obs.z ~tidal_efficiency, data = phylo_sesmpd_env_phys)
summary(lm_phylo_sesmpd_te)

phylo_sesmpd_tidal_efficiency_plot <- ggplot(data = phylo_sesmpd_env, mapping = aes(x = tidal_efficiency, y = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Tidal Efficiency", y= "sesMPD") +
  ylim(c(-6,2))
phylo_sesmpd_tidal_efficiency_plot <- phylo_sesmpd_tidal_efficiency_plot + geom_smooth(mapping = aes(x = tidal_efficiency, y = mpd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
phylo_sesmpd_tidal_efficiency_plot
ggsave("phylo_sesmpd_tidal_efficiency_plot.png", phylo_sesmpd_tidal_efficiency_plot, width = 6, height = 4, units = "in")
```

## Regression of MNTD

``` r
#MNTD
man_phylo_sesmntd_env_chem <- manova(cbind(temperature_median, salinity_median, oxygen_median, pH_median) ~ mntd.obs.z, data = phylo_sesmntd_env)
summary.aov(man_phylo_sesmntd_env_chem)

phylo_sesmntd_env_chem_model <- glmulti("mntd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = phylo_sesmntd_env, method = "h")

man_phylo_sesmntd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth) ~ mntd.obs.z, data = phylo_sesmntd_env)
summary.aov(man_phylo_sesmntd_env_phys)

phylo_sesmntd_env_phys_model <- glmulti("mntd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = phylo_sesmntd_env, method = "h")

#phylo_sesmntd_env_model <- glmulti("mntd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "conductivity_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = phylo_sesmntd_env, method = "h")

phylo_mntd_env_model <- lm(mntd.obs.z ~ -1 + volume_m3_w_chemocline + distance_to_ocean_min_m + tidal_efficiency + salinity_median, data = phylo_sesmntd_env)
summary(phylo_mntd_env_model)
anova(phylo_mntd_env_model)
```

## Plot MNTD vs Oxygen Concentration

``` r
phylo_sesmntd_oxygen_plot <- ggplot(data = phylo_sesmntd_env, mapping = aes(x = oxygen_median, y = mntd.obs.z, color = Environment)) + 
  geom_point(stat = 'identity',
    size = 12,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = .3, end = .95, discrete = TRUE, option = "H") +
  theme_bw() +
  theme(text = element_text(size = 36),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "sesMNTD") +
  ylim(c(-7,1))
phylo_sesmntd_oxygen_plot <- phylo_sesmntd_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = mntd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("phylo_sesmntd_oxygen_plot.png", phylo_sesmntd_oxygen_plot, width = 2, height = 2)
```

## Plot MNTD vs Tidal Efficiency

``` r
phylo_sesmntd_tidal_efficiency_plot <- ggplot(data = phylo_sesmntd_env, mapping = aes(x = tidal_efficiency, y = mntd.obs.z, color = Environment)) + 
  geom_point(stat = 'identity',
    size = 12,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = .3, end = .95, discrete = TRUE, option = "H") +
  theme_bw() +
  theme(text = element_text(size = 36),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Tidal Efficiency", y= "sesMNTD") +
  ylim(c(-7,1))
phylo_sesmntd_tidal_efficiency_plot <- phylo_sesmntd_tidal_efficiency_plot + geom_smooth(mapping = aes(x = tidal_efficiency, y = mntd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("phylo_sesmntd_tidal_efficiency_plot.png", phylo_sesmntd_tidal_efficiency_plot, width = 2, height = 2)
```

## Regression of PD

``` r
#PD
man_phylo_sespd_env_chem <- manova(cbind(temperature_median, salinity_median, oxygen_median, pH_median) ~ pd.obs.z, data = phylo_sespd_env)
summary.aov(man_phylo_sespd_env_chem)

phylo_sespd_env_chem_model <- glmulti("pd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = phylo_sespd_env, method = "h")

man_phylo_sespd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth, logArea) ~ pd.obs.z, data = phylo_sespd_env)
summary.aov(man_phylo_sespd_env_phys)

phylo_sespd_env_phys_model <- glmulti("pd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = phylo_sespd_env, method = "h")

#phylo_sespd_env_model <- glmulti("pd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "conductivity_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = phylo_sespd_env, method = "h")

phylo_sespd_env_model <- lm(pd.obs.z ~ -1 + volume_m3_w_chemocline + distance_to_ocean_min_m + tidal_efficiency + salinity_median, data = phylo_sespd_env)
summary(phylo_pd_env_model)
anova(phylo_pd_env_model)
```

## Plot PD vs Oxygen Concentration

``` r
phylo_sespd_env_plot <- ggplot(data = phylo_sespd_env, mapping = aes(x = env_variable, y = pd.obs.z, color = Environment, shape = Environment, label = X)) + 
  geom_point(stat = 'identity',
    size = 12,
    alpha = 1) + 
  geom_text_repel(min.segment.length = 2, show.legend = FALSE) +
  scale_color_viridis(alpha = 1, begin = .3, end = .95, discrete = TRUE, option = "H") +
  theme_bw() +
  theme(text = element_text(size = 36),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Env Variable", y= "sesFPD")
phylo_sespd_env_plot
ggsave("phylo_sespd_env_plot.png", phylo_sespd_env_plot, width = 2, height = 2)
```

## Plot PD vs Tidal Efficiency

``` r
phylo_sespd_env_plot <- ggplot(data = phylo_sespd_env, mapping = aes(x = env_variable, y = pd.obs.z, color = Environment, shape = Environment, label = X)) + 
  geom_point(stat = 'identity',
    size = 12,
    alpha = 1) + 
  geom_text_repel(min.segment.length = 2, show.legend = FALSE) +
  scale_color_viridis(alpha = 1, begin = .3, end = .95, discrete = TRUE, option = "H") +
  theme_bw() +
  theme(text = element_text(size = 36),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Env Variable", y= "sesFPD")
phylo_sespd_env_plot
ggsave("phylo_sespd_env_plot.png", phylo_sespd_env_plot, width = 2, height = 2)
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

## Beta Diversity of phylogenetic relatedness

``` r
# Do the following to order lakes by the circulation pattern, then delete env data
####https://pedrohbraga.github.io/CommunityPhylogenetics-Workshop/CommunityPhylogenetics-Workshop.html#between-assemblage-phylogenetic-structure####
presabs_env <- merge(pres_abs_by_lake, env_byLake, by = 'X')
presabs_env <- presabs_env[order(presabs_env$Environment),]
row.names(presabs_env) <- presabs_env$X
presabs_env_ordered <- presabs_env[,-c(1, 1692:1722)]

phylo_betap_jac <- phylo.beta.pair(presabs_env_ordered[-23,], phylo_tree, index.family = "jaccard")

# Turnover matrix:
phylo.jac.turn <- phylo_betap_jac$phylo.beta.jtu %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("Sites") + ylab("Sites") + labs(fill = "Turnover") + # edit axis and legend titles
  scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels

# Nestedness matrix:
phylo.jac.nest <- phylo_betap_jac$phylo.beta.jne %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("Sites") + ylab("Sites") + labs(fill = "Nestedness") + # edit axis and legend titles
    scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels

# Beta Diversity matrix:
phylo.jac <- phylo_betap_jac$phylo.beta.jac %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("Sites") + ylab("Sites") + labs(fill = "Beta Diversity") + # edit axis and legend titles
    scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels
phylo.jac

# plot both heatmaps next to each other
phylo_beta_div_plot <- gridExtra::grid.arrange(phylo.jac.turn, phylo.jac.nest, ncol = 2)
ggsave("phylo_beta_div_plot.png", phylo_beta_div_plot, width = 6, height = 4.5)
```
