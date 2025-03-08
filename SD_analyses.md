Species diversity analyses
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

    ##   |                  |          |   0%  |                  |          |   3%                                                                          |                  |.         |   7% [Bringing everything together load modifying files packages]             |                  |.         |  10%                                                                          |                  |.         |  13% [Bringing everything together load in modifying files]                   |                  |..        |  17%                                                                          |                  |..        |  20% [Check species names across files]                                       |                  |..        |  23%                                                                          |                  |...       |  27% [Modify environment data]                                                |                  |...       |  30%                                                                          |                  |...       |  33% [Modify incidence matrices]                                              |                  |....      |  37%                                                                          |                  |....      |  40% [Modify phylogeny]                                                       |                  |....      |  43%                                                                          |                  |.....     |  47% [Modify trait data]                                                      |                  |.....     |  50%                                                                          |                  |.....     |  53% [Modify stratification data frames]                                      |                  |......    |  57%                                                                          |                  |......    |  60% [Modify stratification trait data]                                       |                  |......    |  63%                                                                          |                  |.......   |  67% [Stratification trait data tests]                                        |                  |.......   |  70%                                                                          |                  |.......   |  73% [Modify site trait data frames]                                          |                  |........  |  77%                                                                          |                  |........  |  80% [Site trait data tests]                                                  |                  |........  |  83%                                                                          |                  |......... |  87% [unnamed-chunk-3]                                                        |                  |......... |  90%                                                                          |                  |......... |  93% [unnamed-chunk-4]                                                        |                  |..........|  97%                                                                          |                  |..........| 100% [Bringing everything together load out modified files and session info]

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
library(ggvenn)
```

    ## Loading required package: grid

``` r
library(dplyr)

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

env_cont <- "#FFB90F"

trait_cont <- "#68228B"

trait_cat <- "#DEB887"
```

### Edit input files

``` r
strat_fish <- stratified_fish[,-c(1,10)]
mix_fish <- mixed_fish[,-c(1,10)]
oc_fish <- ocean_fish[,-c(1,8)]

env <- env[order(env$X),]
presabs_lake <- presabs_lake[order(row.names(presabs_lake)),]
surveyed_sites_lake <- surveyed_sites_lake[order(row.names(surveyed_sites_lake)),]

stratification_group_ref <- env[,19]
stratification_group <- env[surveyed_sites,19]
```

# Species Diversity

## SD alpha diversity

### SD alpha bar plot

``` r
#Plot species richness with a bar plot
SR_env$Stratification <- factor(SR_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

(SR_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = reorder(X, row_sum, decreasing = T), y = row_sum, color = Stratification, fill = Stratification, )) + 
  geom_bar(stat = 'identity',
    alpha = 1) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Site", y="Species Richness", colour = "Site type", fill = "Site type"))
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20bar%20plot-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SR_plot.jpg", plot = SR_plot, width = 6, height = 6, units = "in")
```

### SD alpha venn diagram plot

``` r
site_type_pres <- strat_presabs_lake[c(24:26),]
site_type_pres <- t(site_type_pres)
site_type_pres <- site_type_pres[which(rowSums(site_type_pres) > 0),]

# Convert numerical values to logical
site_type_pres_logical <- as.data.frame(site_type_pres > 0)

# Check the structure of the transformed data
str(site_type_pres_logical)
```

    ## 'data.frame':    249 obs. of  3 variables:
    ##  $ Ocean sites     : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ Mixed lakes     : logi  TRUE FALSE TRUE FALSE TRUE FALSE ...
    ##  $ Stratified lakes: logi  FALSE FALSE FALSE FALSE FALSE FALSE ...

``` r
# Create the Venn diagram
(venn_plot <- ggvenn(site_type_pres_logical,
                    show_percentage = F,
                    fill_color = c("#EE6363", "#87CEFA", "#6E8B3D"),
                    fill_alpha = 0.7,
                    stroke_alpha = 0,
                    stroke_size = 0.5, 
                    set_name_size = 5,
                    text_size = 5))
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20venn%20diagram%20plot-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/venn_plot.jpg", plot = venn_plot, width = 4.25, height = 4, units = "in")
```

### SD alpha subsample venn diagram plots

``` r
# Randomly sample six columns from stratified_fish and mixed_fish
sampled_stratified <- sample(names(strat_fish)[-1], 6)
sampled_mixed <- sample(names(mix_fish)[-1], 6)

# Select sampled columns from stratified_fish and mixed_fish
stratified_subset <- select(stratified_fish, all_of(sampled_stratified))
mixed_subset <- select(mixed_fish, all_of(sampled_mixed))

# Get rid of species not present in the 6 remaining sites
stratified_subset <- stratified_subset[which(rowSums(stratified_subset[, -1]) > 0),]
mixed_subset <- mixed_subset[which(rowSums(mixed_subset[, -1]) > 0),]

# Make site type columns in each dataframe
stratified_subset[,7] <- rowSums(stratified_subset[,c(1:6)])
colnames(stratified_subset)[colnames(stratified_subset) == "V7"] <- "Stratified lakes"

mixed_subset[,7] <- rowSums(mixed_subset[,c(1:6)])
colnames(mixed_subset)[colnames(mixed_subset) == "V7"] <- "Mixed lakes"

oc_fish[,7] <- rowSums(oc_fish[,c(1:6)])
colnames(oc_fish)[colnames(oc_fish) == "V7"] <- "Ocean sites"

# Merge the two data frames, replacing missing values with 0
oc_mix_merge <- merge(oc_fish, mixed_subset, by = 0, all = TRUE)
row.names(oc_mix_merge) <- oc_mix_merge$Row.names
oc_mix_merge <- oc_mix_merge[,-1]
oc_mix_strat_merge <- merge(oc_mix_merge, stratified_subset, by = 0, all = TRUE)
row.names(oc_mix_strat_merge) <- oc_mix_strat_merge$Row.names
oc_mix_strat_merge <- oc_mix_strat_merge[,-1]

# Keep only site type
keep <- c("Ocean sites", "Mixed lakes", "Stratified lakes")
site_type_pres_sample <- oc_mix_strat_merge[,keep]

# Replace NA with 0
site_type_pres_sample[is.na(site_type_pres_sample)] <- 0

## Remove all species not found in any locations
# Identifies which rows are greater than 0
site_type_pres_sample <- site_type_pres_sample[which(rowSums(site_type_pres_sample) > 0),]

site_type_pres_sample_logical <- as.data.frame(site_type_pres_sample > 0)

venn_plot <- ggvenn(site_type_pres_sample_logical,
                    show_percentage = F,
                    fill_color = c("#EE6363", "#87CEFA", "#6E8B3D"),
                    fill_alpha = 0.7,
                    stroke_alpha = 0,
                    stroke_size = 0.5, 
                    set_name_size = 5,
                    text_size = 5)
venn_plot
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20venn%20diagram%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/venn_plot_s10.jpg", plot = venn_plot, width = 4.25, height = 4, units = "in")
```

### SD alpha outliers for each stratification type

``` r
outlier_SD_alpha <- SR_env[surveyed_sites,] %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(row_sum, 0.25),
    Q3 = quantile(row_sum, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = row_sum < lower_bound | row_sum > upper_bound
  )

outlier_SD_alpha <- as.data.frame(outlier_SD_alpha)
row.names(outlier_SD_alpha) <- outlier_SD_alpha$X

# Create the plot
(outlier_SD_alpha_plot <- ggplot(outlier_SD_alpha, aes(x = Stratification, y = row_sum, fill = Stratification)) +
  geom_violin(trim = FALSE) +
  geom_jitter(aes(color = is_outlier, fill = is_outlier), width = 0.2, size = 4, alpha = 0.7) +
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
  guides(fill = "none") + 
  labs(y = "Species Richness", x = "Site type", color = "Outlier", fill = "Site type", tag = "a"))
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/outlier_SD_alpha.jpg", outlier_SD_alpha_plot, width = 6.26, height = 6, units = "in")
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

mod <-  lm(row_sum ~ temperature_median + salinity_median + oxygen_median + pH_median, SR_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = SR_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -25.9542 -10.3708  -0.4769   8.9325  29.5684 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -541.813    201.886  -2.684   0.0178 *
    ## temperature_median   -4.950      3.874  -1.278   0.2221  
    ## salinity_median       0.201      1.973   0.102   0.9203  
    ## oxygen_median         6.058      5.343   1.134   0.2758  
    ## pH_median            90.737     34.922   2.598   0.0210 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 16.28 on 14 degrees of freedom
    ##   (4 observations deleted due to missingness)
    ## Multiple R-squared:  0.7996, Adjusted R-squared:  0.7424 
    ## F-statistic: 13.97 on 4 and 14 DF,  p-value: 8.549e-05

``` r
anova(mod)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                    Df Sum Sq Mean Sq F value    Pr(>F)    
    ## temperature_median  1 1618.4  1618.4  6.1049 0.0269418 *  
    ## salinity_median     1 7595.9  7595.9 28.6536 0.0001019 ***
    ## oxygen_median       1 3808.8  3808.8 14.3679 0.0019877 ** 
    ## pH_median           1 1789.6  1789.6  6.7509 0.0210487 *  
    ## Residuals          14 3711.3   265.1                      
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
mod <-  lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -31.549  -9.598  -4.173   9.619  47.290 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        -131.1738   147.7497  -0.888  0.38866   
    ## temperature_median   -0.7881     4.1484  -0.190  0.85187   
    ## salinity_median       4.1712     1.4680   2.841  0.01238 * 
    ## oxygen_median        15.2171     4.7218   3.223  0.00569 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 19.15 on 15 degrees of freedom
    ##   (4 observations deleted due to missingness)
    ## Multiple R-squared:  0.703,  Adjusted R-squared:  0.6436 
    ## F-statistic: 11.84 on 3 and 15 DF,  p-value: 0.000309

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

mod <-  lm(row_sum ~ volume_m3_w_chemocline + volume_m3 + surface_area_m2 + distance_to_ocean_min_m + distance_to_ocean_mean_m + distance_to_ocean_median_m + tidal_lag_time_minutes + tidal_efficiency + perimeter_fromSat + max_depth + logArea, SR_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ volume_m3_w_chemocline + volume_m3 + surface_area_m2 + 
    ##     distance_to_ocean_min_m + distance_to_ocean_mean_m + distance_to_ocean_median_m + 
    ##     tidal_lag_time_minutes + tidal_efficiency + perimeter_fromSat + 
    ##     max_depth + logArea, data = SR_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -27.5820  -6.9408   0.9783   8.2989  28.1100 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                -1.263e+02  1.267e+02  -0.997   0.3447  
    ## volume_m3_w_chemocline      5.099e-05  4.629e-05   1.102   0.2992  
    ## volume_m3                  -4.415e-05  4.795e-05  -0.921   0.3812  
    ## surface_area_m2            -8.853e-05  4.614e-05  -1.919   0.0872 .
    ## distance_to_ocean_min_m     3.566e-02  2.530e-01   0.141   0.8910  
    ## distance_to_ocean_mean_m    1.580e-01  1.075e+00   0.147   0.8863  
    ## distance_to_ocean_median_m -2.599e-01  9.878e-01  -0.263   0.7984  
    ## tidal_lag_time_minutes     -9.226e-02  2.965e-01  -0.311   0.7628  
    ## tidal_efficiency            8.936e+00  7.900e+01   0.113   0.9124  
    ## perimeter_fromSat          -1.692e-02  1.953e-02  -0.866   0.4087  
    ## max_depth                  -3.370e-01  1.301e+00  -0.259   0.8015  
    ## logArea                     2.046e+01  1.264e+01   1.619   0.1399  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 21.04 on 9 degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## Multiple R-squared:  0.7342, Adjusted R-squared:  0.4093 
    ## F-statistic:  2.26 on 11 and 9 DF,  p-value: 0.1156

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
mod <-  lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -42.936 -17.133   0.982  10.144  68.145 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             14.99041   36.14360   0.415   0.6832  
    ## distance_to_ocean_min_m -0.22020    0.08150  -2.702   0.0146 *
    ## max_depth                0.09856    0.59911   0.165   0.8712  
    ## logArea                  3.66034    3.59139   1.019   0.3216  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 25.95 on 18 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.3943, Adjusted R-squared:  0.2933 
    ## F-statistic: 3.906 on 3 and 18 DF,  p-value: 0.02604

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

### SD alpha SRic & stratification

``` r
### Stratification
# Run ANOVA
anova_logSRic_result <- aov(log(row_sum) ~ Stratification, data = SR_env[surveyed_sites,])
summary(anova_logSRic_result)
```

    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Stratification  2 22.346  11.173   27.41 2.51e-06 ***
    ## Residuals      19  7.745   0.408                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Perform Tukey's HSD test
tukey_logSRic_result <- TukeyHSD(anova_logSRic_result)
print(tukey_logSRic_result)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = log(row_sum) ~ Stratification, data = SR_env[surveyed_sites, ])
    ## 
    ## $Stratification
    ##                         diff        lwr        upr     p adj
    ## Mixed-Ocean      -0.03070513 -0.9066479  0.8452376 0.9956384
    ## Stratified-Ocean -2.11249529 -2.9884380 -1.2365526 0.0000197
    ## Stratified-Mixed -2.08179016 -2.8927555 -1.2708248 0.0000087

``` r
# Sites
N <- length(anova_logSRic_result$residuals)
# Categories
k <- length(anova_logSRic_result$coefficients)
# Sum of squares
SS <- sum(resid(anova_logSRic_result)^2) # from anova, the residual SS
# Mean squares
MS <- SS/(N-k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

# Get differences
Stratification <- tukey_logSRic_result$Stratification

# q-values
(OM_q <- (abs(Stratification[1,1]))/USE)
```

    ## [1] 0.1259392

``` r
(MS_q <- (abs(Stratification[3,1]))/BSE)
```

    ## [1] 9.222747

``` r
(SO_q <- (abs(Stratification[2,1]))/USE)
```

    ## [1] 8.664544

``` r
# Run ANOVA
anova_SRic_result <- aov(row_sum ~ Stratification, data = SR_env[surveyed_sites,])
summary(anova_SRic_result)
```

    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Stratification  2  10743    5371   11.01 0.000667 ***
    ## Residuals      19   9268     488                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Perform Tukey's HSD test
tukey_SRic_result <- TukeyHSD(anova_SRic_result)
print(tukey_SRic_result)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = row_sum ~ Stratification, data = SR_env[surveyed_sites, ])
    ## 
    ## $Stratification
    ##                        diff       lwr       upr     p adj
    ## Mixed-Ocean        1.041667 -29.25978  31.34312 0.9958049
    ## Stratified-Ocean -45.333333 -75.63478 -15.03188 0.0033013
    ## Stratified-Mixed -46.375000 -74.42869 -18.32131 0.0013468

``` r
# Sites
N <- length(anova_SRic_result$residuals)
# Categories
k <- length(anova_SRic_result$coefficients)
# Sum of squares
SS <- sum(resid(anova_SRic_result)^2) # from anova, the residual SS
# Mean squares
MS <- SS/(N-k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

# Get differences
Stratification <- tukey_SRic_result$Stratification

# q-values
(OM_q <- (abs(Stratification[1,1]))/USE)
```

    ## [1] 0.1235068

``` r
(MS_q <- (abs(Stratification[3,1]))/BSE)
```

    ## [1] 5.939085

``` r
(SO_q <- (abs(Stratification[2,1]))/USE)
```

    ## [1] 5.375017

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx
```

### SD alpha SRic env and geo linear models

``` r
### Environmental
# Surveyed sites
SD_logSRic_env_model <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[surveyed_sites_env,])
summary(SD_logSRic_env_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[surveyed_sites_env, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.73012 -0.20864 -0.06934  0.27642  0.88318 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -7.63922    3.81124  -2.004  0.06343 .  
    ## temperature_median  0.03262    0.10701   0.305  0.76470    
    ## salinity_median     0.26035    0.03787   6.876 5.27e-06 ***
    ## oxygen_median       0.39484    0.12180   3.242  0.00548 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.494 on 15 degrees of freedom
    ## Multiple R-squared:  0.8715, Adjusted R-squared:  0.8458 
    ## F-statistic: 33.92 on 3 and 15 DF,  p-value: 6.32e-07

``` r
anova(SD_logSRic_env_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                    Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## temperature_median  1  2.5993  2.5993  10.652  0.005234 ** 
    ## salinity_median     1 19.6670 19.6670  80.596 2.025e-07 ***
    ## oxygen_median       1  2.5644  2.5644  10.509  0.005475 ** 
    ## Residuals          15  3.6603  0.2440                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##       2.537061e-01       1.000000e+00       2.109736e-05       2.190101e-02

``` r
SD_SRic_env_model <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[surveyed_sites_env,])
summary(SD_SRic_env_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[surveyed_sites_env, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -31.549  -9.598  -4.173   9.619  47.290 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        -131.1738   147.7497  -0.888  0.38866   
    ## temperature_median   -0.7881     4.1484  -0.190  0.85187   
    ## salinity_median       4.1712     1.4680   2.841  0.01238 * 
    ## oxygen_median        15.2171     4.7218   3.223  0.00569 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 19.15 on 15 degrees of freedom
    ## Multiple R-squared:  0.703,  Adjusted R-squared:  0.6436 
    ## F-statistic: 11.84 on 3 and 15 DF,  p-value: 0.000309

``` r
anova(SD_SRic_env_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                    Df Sum Sq Mean Sq F value    Pr(>F)    
    ## temperature_median  1 1618.4  1618.4   4.413 0.0529836 .  
    ## salinity_median     1 7595.9  7595.9  20.712 0.0003824 ***
    ## oxygen_median       1 3808.8  3808.8  10.386 0.0056926 ** 
    ## Residuals          15 5500.9   366.7                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##         1.00000000         1.00000000         0.04950985         0.02277051

``` r
# Mixed and stratified lakes
SD_logSRic_env_MS_model <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_env_MS_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[mixed_stratified_lakes, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.71621 -0.27163 -0.06681  0.37624  0.91222 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -7.73929    4.51110  -1.716   0.1119    
    ## temperature_median  0.03109    0.12473   0.249   0.8074    
    ## salinity_median     0.26303    0.04265   6.168 4.81e-05 ***
    ## oxygen_median       0.41650    0.15412   2.702   0.0192 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5464 on 12 degrees of freedom
    ## Multiple R-squared:  0.8494, Adjusted R-squared:  0.8117 
    ## F-statistic: 22.56 on 3 and 12 DF,  p-value: 3.194e-05

``` r
anova(SD_logSRic_env_MS_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                    Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## temperature_median  1  2.8913  2.8913  9.6854  0.008985 ** 
    ## salinity_median     1 15.1315 15.1315 50.6883 1.214e-05 ***
    ## oxygen_median       1  2.1802  2.1802  7.3032  0.019222 *  
    ## Residuals          12  3.5822  0.2985                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_MS_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##       0.4476551535       1.0000000000       0.0001925785       0.0768885411

``` r
SD_SRic_env_MS_model <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[mixed_stratified_lakes,])
summary(SD_SRic_env_MS_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[mixed_stratified_lakes, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -31.897  -9.129  -2.596   7.932  45.355 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)         -94.826    168.390  -0.563   0.5837  
    ## temperature_median   -1.931      4.656  -0.415   0.6857  
    ## salinity_median       4.176      1.592   2.623   0.0223 *
    ## oxygen_median        15.034      5.753   2.613   0.0227 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 20.39 on 12 degrees of freedom
    ## Multiple R-squared:  0.6648, Adjusted R-squared:  0.581 
    ## F-statistic: 7.934 on 3 and 12 DF,  p-value: 0.003509

``` r
anova(SD_SRic_env_MS_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                    Df Sum Sq Mean Sq F value   Pr(>F)   
    ## temperature_median  1 2184.6  2184.6  5.2521 0.040803 * 
    ## salinity_median     1 4875.1  4875.1 11.7205 0.005046 **
    ## oxygen_median       1 2840.3  2840.3  6.8286 0.022669 * 
    ## Residuals          12 4991.4   415.9                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_MS_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##         1.00000000         1.00000000         0.08902474         0.09067730

``` r
# Ocean sites and mixed lakes
SD_logSRic_env_OM_model <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_env_OM_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[ocean_mixed_sites_env, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.52504 -0.28816 -0.07468  0.20202  0.64145 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -8.29859   18.14153  -0.457    0.661
    ## temperature_median -0.09415    0.25046  -0.376    0.718
    ## salinity_median     0.39271    0.46891   0.838    0.430
    ## oxygen_median       0.43003    0.37824   1.137    0.293
    ## 
    ## Residual standard error: 0.4802 on 7 degrees of freedom
    ## Multiple R-squared:  0.466,  Adjusted R-squared:  0.2371 
    ## F-statistic: 2.036 on 3 and 7 DF,  p-value: 0.1975

``` r
anova(SD_logSRic_env_OM_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                    Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## temperature_median  1 0.00311 0.00311  0.0135 0.91083  
    ## salinity_median     1 1.10728 1.10728  4.8018 0.06456 .
    ## oxygen_median       1 0.29806 0.29806  1.2926 0.29299  
    ## Residuals           7 1.61419 0.23060                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_OM_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
SD_SRic_env_OM_model <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[ocean_mixed_sites_env,])
summary(SD_SRic_env_OM_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[ocean_mixed_sites_env, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -24.059 -18.215  -3.151  11.182  38.375 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         -37.380    926.047  -0.040    0.969
    ## temperature_median   -8.798     12.785  -0.688    0.514
    ## salinity_median       7.098     23.936   0.297    0.775
    ## oxygen_median        26.197     19.308   1.357    0.217
    ## 
    ## Residual standard error: 24.51 on 7 degrees of freedom
    ## Multiple R-squared:  0.4012, Adjusted R-squared:  0.1446 
    ## F-statistic: 1.563 on 3 and 7 DF,  p-value: 0.2814

``` r
anova(SD_SRic_env_OM_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                    Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1   60.4   60.40  0.1005 0.7605
    ## salinity_median     1 1651.6 1651.59  2.7487 0.1413
    ## oxygen_median       1 1106.2 1106.16  1.8410 0.2170
    ## Residuals           7 4206.0  600.86

``` r
p_values <- summary(SD_SRic_env_OM_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          1.0000000          1.0000000          0.8678638

``` r
# Stratified lakes and ocean sites
SD_logSRic_env_SO_model <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_env_SO_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[ocean_stratified_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.74751 -0.14210 -0.02821  0.14682  0.84158 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -8.58582    4.18049  -2.054 0.079087 .  
    ## temperature_median  0.07539    0.11959   0.630 0.548449    
    ## salinity_median     0.25128    0.04292   5.854 0.000628 ***
    ## oxygen_median       0.34943    0.13583   2.573 0.036866 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5124 on 7 degrees of freedom
    ## Multiple R-squared:  0.8875, Adjusted R-squared:  0.8392 
    ## F-statistic:  18.4 on 3 and 7 DF,  p-value: 0.001063

``` r
anova(SD_logSRic_env_SO_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                    Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## temperature_median  1  0.3028  0.3028  1.1537 0.3184081    
    ## salinity_median     1 12.4503 12.4503 47.4289 0.0002341 ***
    ## oxygen_median       1  1.7374  1.7374  6.6184 0.0368663 *  
    ## Residuals           7  1.8375  0.2625                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_SO_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##        0.316349339        1.000000000        0.002511461        0.147465030

``` r
SD_SRic_env_SO_model <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[ocean_stratified_sites_env,])
summary(SD_SRic_env_SO_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[ocean_stratified_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -12.023  -6.847  -2.633   3.538  17.340 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        -206.9563    94.1879  -2.197  0.06399 . 
    ## temperature_median    1.8111     2.6945   0.672  0.52306   
    ## salinity_median       4.2382     0.9671   4.383  0.00322 **
    ## oxygen_median        13.2060     3.0602   4.315  0.00350 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 11.54 on 7 degrees of freedom
    ## Multiple R-squared:  0.8848, Adjusted R-squared:  0.8354 
    ## F-statistic: 17.92 on 3 and 7 DF,  p-value: 0.001153

``` r
anova(SD_SRic_env_SO_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                    Df Sum Sq Mean Sq F value    Pr(>F)    
    ## temperature_median  1   84.0    84.0  0.6302 0.4533626    
    ## salinity_median     1 4597.8  4597.8 34.5044 0.0006154 ***
    ## oxygen_median       1 2481.5  2481.5 18.6227 0.0034998 ** 
    ## Residuals           7  932.8   133.3                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_SO_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##         0.25595351         1.00000000         0.01289859         0.01399934

``` r
# Mixed lakes
SD_logSRic_env_M_model <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[mixed_lakes,])
summary(SD_logSRic_env_M_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##     FLK     HLO     LLN     MLN     NLN     NLU     OLO     ULN 
    ##  0.1488 -0.2163  0.5206 -0.3148  0.6425 -0.4922 -0.4585  0.1700 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -14.6394    27.2353  -0.538    0.619
    ## temperature_median  -0.0965     0.3930  -0.246    0.818
    ## salinity_median      0.5768     0.6332   0.911    0.414
    ## oxygen_median        0.5189     0.5122   1.013    0.368
    ## 
    ## Residual standard error: 0.5774 on 4 degrees of freedom
    ## Multiple R-squared:  0.4991, Adjusted R-squared:  0.1234 
    ## F-statistic: 1.329 on 3 and 4 DF,  p-value: 0.3825

``` r
anova(SD_logSRic_env_M_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                    Df  Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1 0.24929 0.24929  0.7479 0.4359
    ## salinity_median     1 0.73706 0.73706  2.2112 0.2112
    ## oxygen_median       1 0.34222 0.34222  1.0267 0.3683
    ## Residuals           4 1.33334 0.33334

``` r
p_values <- summary(SD_logSRic_env_M_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
SD_SRic_env_M_model <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[mixed_lakes,])
summary(SD_SRic_env_M_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##     FLK     HLO     LLN     MLN     NLN     NLU     OLO     ULN 
    ##   7.440 -13.112  28.023 -11.289  30.672 -23.921 -20.011   2.198 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)          -32.78    1304.37  -0.025    0.981
    ## temperature_median   -15.30      18.82  -0.813    0.462
    ## salinity_median       12.84      30.32   0.424    0.694
    ## oxygen_median         27.71      24.53   1.130    0.322
    ## 
    ## Residual standard error: 27.65 on 4 degrees of freedom
    ## Multiple R-squared:  0.4914, Adjusted R-squared:  0.1099 
    ## F-statistic: 1.288 on 3 and 4 DF,  p-value: 0.3928

``` r
anova(SD_SRic_env_M_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                    Df  Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1 1209.80 1209.80  1.5823 0.2769
    ## salinity_median     1  769.05  769.05  1.0058 0.3726
    ## oxygen_median       1  975.73  975.73  1.2762 0.3218
    ## Residuals           4 3058.30  764.57

``` r
p_values <- summary(SD_SRic_env_M_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
# Stratified lakes
SD_logSRic_env_S_model <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[stratified_lakes,])
summary(SD_logSRic_env_S_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##      BCM      CLM      GLK      HLM      NLK      OTM      SLN      TLN 
    ##  0.37953  0.04308  0.17921  0.24651 -0.28497 -0.87230 -0.12917  0.43811 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -3.19314    6.27810  -0.509    0.638
    ## temperature_median  0.04973    0.13905   0.358    0.739
    ## salinity_median     0.14298    0.09750   1.466    0.216
    ## oxygen_median      -0.16727    0.42987  -0.389    0.717
    ## 
    ## Residual standard error: 0.5678 on 4 degrees of freedom
    ## Multiple R-squared:  0.6596, Adjusted R-squared:  0.4042 
    ## F-statistic: 2.583 on 3 and 4 DF,  p-value: 0.1908

``` r
anova(SD_logSRic_env_S_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                    Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## temperature_median  1 0.00289 0.00289  0.0090 0.92916  
    ## salinity_median     1 2.44663 2.44663  7.5892 0.05112 .
    ## oxygen_median       1 0.04881 0.04881  0.1514 0.71702  
    ## Residuals           4 1.28953 0.32238                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_S_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          1.0000000          0.8657002          1.0000000

``` r
SD_SRic_env_S_model <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[stratified_lakes,])
summary(SD_SRic_env_S_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##     BCM     CLM     GLK     HLM     NLK     OTM     SLN     TLN 
    ##  2.3022  1.9724  0.9048  2.9002 -3.2143 -6.6194 -1.4092  3.1633 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -12.67663   50.83889  -0.249    0.815
    ## temperature_median  -0.02215    1.12601  -0.020    0.985
    ## salinity_median      1.04089    0.78954   1.318    0.258
    ## oxygen_median       -2.41766    3.48100  -0.695    0.526
    ## 
    ## Residual standard error: 4.598 on 4 degrees of freedom
    ## Multiple R-squared:  0.6936, Adjusted R-squared:  0.4638 
    ## F-statistic: 3.019 on 3 and 4 DF,  p-value: 0.1568

``` r
anova(SD_SRic_env_S_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                    Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## temperature_median  1   4.710   4.710  0.2228 0.66150  
    ## salinity_median     1 176.532 176.532  8.3506 0.04457 *
    ## oxygen_median       1  10.197  10.197  0.4824 0.52558  
    ## Residuals           4  84.560  21.140                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_S_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
### Geographical
# Surveyed sites
SD_logSRic_geo_model <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[surveyed_sites_geo,])
summary(SD_logSRic_geo_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[surveyed_sites_geo, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.70334 -0.46841  0.06636  0.39118  0.88860 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              3.909348   1.051170   3.719  0.00171 ** 
    ## distance_to_ocean_min_m -0.012065   0.002357  -5.119 8.55e-05 ***
    ## max_depth                0.003980   0.017364   0.229  0.82145    
    ## logArea                  0.001813   0.105246   0.017  0.98646    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7453 on 17 degrees of freedom
    ## Multiple R-squared:  0.6565, Adjusted R-squared:  0.5959 
    ## F-statistic: 10.83 on 3 and 17 DF,  p-value: 0.0003245

``` r
anova(SD_logSRic_geo_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                         Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## distance_to_ocean_min_m  1 18.0016 18.0016 32.4116 2.641e-05 ***
    ## max_depth                1  0.0434  0.0434  0.0781    0.7832    
    ## logArea                  1  0.0002  0.0002  0.0003    0.9865    
    ## Residuals               17  9.4419  0.5554                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_geo_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             0.006823287             0.000341928             1.000000000 
    ##                 logArea 
    ##             1.000000000

``` r
SD_SRic_geo_model <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[surveyed_sites_geo,])
summary(SD_SRic_geo_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[surveyed_sites_geo, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -32.593 -17.991   4.612  13.843  31.925 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)             30.83775   28.70356   1.074  0.29768   
    ## distance_to_ocean_min_m -0.24672    0.06436  -3.834  0.00133 **
    ## max_depth                0.32181    0.47414   0.679  0.50646   
    ## logArea                  1.65805    2.87389   0.577  0.57154   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 20.35 on 17 degrees of freedom
    ## Multiple R-squared:  0.5301, Adjusted R-squared:  0.4472 
    ## F-statistic: 6.393 on 3 and 17 DF,  p-value: 0.004253

``` r
anova(SD_SRic_geo_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                         Df Sum Sq Mean Sq F value    Pr(>F)    
    ## distance_to_ocean_min_m  1 7257.0  7257.0 17.5235 0.0006199 ***
    ## max_depth                1  548.1   548.1  1.3236 0.2658752    
    ## logArea                  1  137.8   137.8  0.3329 0.5715436    
    ## Residuals               17 7040.2   414.1                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             1.000000000             0.005321597             1.000000000 
    ##                 logArea 
    ##             1.000000000

``` r
# Mixed and stratified lakes
SD_logSRic_geo_MS_model <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[mixed_stratified_lakes_geo,])
summary(SD_logSRic_geo_MS_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[mixed_stratified_lakes_geo, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6198 -0.3072  0.2026  0.3977  0.9029 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)              1.492123   2.345867   0.636  0.53774   
    ## distance_to_ocean_min_m -0.012616   0.003058  -4.126  0.00168 **
    ## max_depth               -0.032934   0.030515  -1.079  0.30356   
    ## logArea                  0.325574   0.281517   1.156  0.27198   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.802 on 11 degrees of freedom
    ## Multiple R-squared:  0.6473, Adjusted R-squared:  0.5511 
    ## F-statistic: 6.729 on 3 and 11 DF,  p-value: 0.00766

``` r
anova(SD_logSRic_geo_MS_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                         Df  Sum Sq Mean Sq F value   Pr(>F)   
    ## distance_to_ocean_min_m  1 12.0787 12.0787 18.7781 0.001188 **
    ## max_depth                1  0.0454  0.0454  0.0705 0.795496   
    ## logArea                  1  0.8603  0.8603  1.3375 0.271977   
    ## Residuals               11  7.0756  0.6432                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_geo_MS_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             1.000000000             0.006732835             1.000000000 
    ##                 logArea 
    ##             1.000000000

``` r
SD_SRic_geo_MS_model <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[mixed_stratified_lakes_geo,])
summary(SD_SRic_geo_MS_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[mixed_stratified_lakes_geo, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -26.141 -14.940   4.736   9.903  26.910 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)             -35.90394   54.17233  -0.663  0.52112   
    ## distance_to_ocean_min_m  -0.24598    0.07061  -3.484  0.00512 **
    ## max_depth                -0.82801    0.70467  -1.175  0.26480   
    ## logArea                  10.61515    6.50098   1.633  0.13077   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 18.52 on 11 degrees of freedom
    ## Multiple R-squared:  0.5732, Adjusted R-squared:  0.4568 
    ## F-statistic: 4.925 on 3 and 11 DF,  p-value: 0.02084

``` r
anova(SD_SRic_geo_MS_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                         Df Sum Sq Mean Sq F value   Pr(>F)   
    ## distance_to_ocean_min_m  1 4140.1  4140.1 12.0695 0.005202 **
    ## max_depth                1   13.5    13.5  0.0393 0.846470   
    ## logArea                  1  914.6   914.6  2.6662 0.130769   
    ## Residuals               11 3773.2   343.0                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_MS_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              1.00000000              0.02046022              1.00000000 
    ##                 logArea 
    ##              0.52307591

``` r
# Ocean sites and mixed lakes
SD_logSRic_geo_OM_model <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_mixed_sites_geo,])
summary(SD_logSRic_geo_OM_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[ocean_mixed_sites_geo, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.58045 -0.23931 -0.02709  0.21184  0.60673 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)              3.328124   0.767088   4.339  0.00188 **
    ## distance_to_ocean_min_m -0.002071   0.003761  -0.551  0.59532   
    ## max_depth                0.031774   0.013097   2.426  0.03823 * 
    ## logArea                  0.010811   0.067896   0.159  0.87700   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4223 on 9 degrees of freedom
    ## Multiple R-squared:  0.5054, Adjusted R-squared:  0.3405 
    ## F-statistic: 3.065 on 3 and 9 DF,  p-value: 0.08386

``` r
anova(SD_logSRic_geo_OM_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                         Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## distance_to_ocean_min_m  1 0.42772 0.42772  2.3983 0.15587  
    ## max_depth                1 1.20777 1.20777  6.7723 0.02863 *
    ## logArea                  1 0.00452 0.00452  0.0254 0.87700  
    ## Residuals                9 1.60507 0.17834                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_geo_OM_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             0.007524722             1.000000000             0.152931616 
    ##                 logArea 
    ##             1.000000000

``` r
SD_SRic_geo_OM_model <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_mixed_sites_geo,])
summary(SD_SRic_geo_OM_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[ocean_mixed_sites_geo, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -27.663 -10.323  -1.235   7.805  26.814 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             17.49338   31.98707   0.547   0.5978  
    ## distance_to_ocean_min_m -0.05212    0.15682  -0.332   0.7472  
    ## max_depth                1.34837    0.54615   2.469   0.0356 *
    ## logArea                  1.52769    2.83122   0.540   0.6026  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 17.61 on 9 degrees of freedom
    ## Multiple R-squared:  0.5363, Adjusted R-squared:  0.3818 
    ## F-statistic:  3.47 on 3 and 9 DF,  p-value: 0.06405

``` r
anova(SD_SRic_geo_OM_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                         Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## distance_to_ocean_min_m  1  756.80  756.80  2.4405 0.15268  
    ## max_depth                1 2381.04 2381.04  7.6782 0.02172 *
    ## logArea                  1   90.29   90.29  0.2912 0.60257  
    ## Residuals                9 2790.95  310.11                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_OM_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               1.0000000               0.1425413 
    ##                 logArea 
    ##               1.0000000

``` r
# Stratified lakes and ocean sites
SD_logSRic_geo_SO_model <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_stratified_sites,])
summary(SD_logSRic_geo_SO_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[ocean_stratified_sites, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.5114 -0.1857  0.1292  0.4833  0.6771 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)              2.807256   1.265080   2.219  0.05078 . 
    ## distance_to_ocean_min_m -0.011234   0.002530  -4.440  0.00126 **
    ## max_depth                0.007569   0.019657   0.385  0.70826   
    ## logArea                  0.066059   0.118135   0.559  0.58834   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.717 on 10 degrees of freedom
    ## Multiple R-squared:  0.7478, Adjusted R-squared:  0.6721 
    ## F-statistic: 9.883 on 3 and 10 DF,  p-value: 0.002455

``` r
anova(SD_logSRic_geo_SO_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                         Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## distance_to_ocean_min_m  1 14.8167 14.8167 28.8208 0.0003152 ***
    ## max_depth                1  0.2646  0.2646  0.5148 0.4895003    
    ## logArea                  1  0.1608  0.1608  0.3127 0.5883374    
    ## Residuals               10  5.1410  0.5141                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_geo_SO_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             0.203110617             0.005020788             1.000000000 
    ##                 logArea 
    ##             1.000000000

``` r
SD_SRic_geo_SO_model <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_stratified_sites,])
summary(SD_SRic_geo_SO_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[ocean_stratified_sites, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -28.019 -12.884   1.303  13.644  26.058 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             10.04799   35.53443   0.283   0.7831  
    ## distance_to_ocean_min_m -0.21631    0.07107  -3.044   0.0124 *
    ## max_depth                0.29471    0.55214   0.534   0.6052  
    ## logArea                  2.81947    3.31824   0.850   0.4154  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 20.14 on 10 degrees of freedom
    ## Multiple R-squared:  0.6062, Adjusted R-squared:  0.4881 
    ## F-statistic: 5.132 on 3 and 10 DF,  p-value: 0.02099

``` r
anova(SD_SRic_geo_SO_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                         Df Sum Sq Mean Sq F value   Pr(>F)   
    ## distance_to_ocean_min_m  1 5519.5  5519.5 13.6079 0.004184 **
    ## max_depth                1  432.5   432.5  1.0662 0.326125   
    ## logArea                  1  292.8   292.8  0.7220 0.415378   
    ## Residuals               10 4056.1   405.6                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_SO_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              1.00000000              0.04954771              1.00000000 
    ##                 logArea 
    ##              1.00000000

``` r
# Mixed lakes
SD_logSRic_geo_M_model <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[mixed_lakes_geo,])
summary(SD_logSRic_geo_M_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[mixed_lakes_geo, ])
    ## 
    ## Residuals:
    ##      FLK      HLO      MLN      NLN      NLU      OLO      ULN 
    ##  0.34018 -0.20253  0.07438  0.06918 -0.04383 -0.49408  0.25669 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              1.760524   1.369824   1.285    0.289
    ## distance_to_ocean_min_m -0.005783   0.007053  -0.820    0.472
    ## max_depth                0.009129   0.028794   0.317    0.772
    ## logArea                  0.236244   0.164058   1.440    0.245
    ## 
    ## Residual standard error: 0.3996 on 3 degrees of freedom
    ## Multiple R-squared:  0.745,  Adjusted R-squared:  0.4901 
    ## F-statistic: 2.922 on 3 and 3 DF,  p-value: 0.201

``` r
anova(SD_logSRic_geo_M_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.64740 0.64740  4.0548 0.1375
    ## max_depth                1 0.42113 0.42113  2.6376 0.2028
    ## logArea                  1 0.33108 0.33108  2.0736 0.2455
    ## Residuals                3 0.47899 0.15966

``` r
p_values <- summary(SD_logSRic_geo_M_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               1.0000000               1.0000000 
    ##                 logArea 
    ##               0.9819636

``` r
SD_SRic_geo_M_model <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[mixed_lakes_geo,])
summary(SD_SRic_geo_M_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[mixed_lakes_geo, ])
    ## 
    ## Residuals:
    ##     FLK     HLO     MLN     NLN     NLU     OLO     ULN 
    ##  10.283 -10.556   5.711   1.971  -0.466 -17.170  10.227 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)             -37.07542   50.58903  -0.733    0.517
    ## distance_to_ocean_min_m  -0.28145    0.26048  -1.080    0.359
    ## max_depth                 0.04393    1.06339   0.041    0.970
    ## logArea                  10.77730    6.05885   1.779    0.173
    ## 
    ## Residual standard error: 14.76 on 3 degrees of freedom
    ## Multiple R-squared:  0.7753, Adjusted R-squared:  0.5506 
    ## F-statistic: 3.451 on 3 and 3 DF,  p-value: 0.1681

``` r
anova(SD_SRic_geo_M_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 1043.18 1043.18  4.7904 0.1164
    ## max_depth                1  522.23  522.23  2.3981 0.2192
    ## logArea                  1  689.01  689.01  3.1640 0.1733
    ## Residuals                3  653.29  217.76

``` r
p_values <- summary(SD_SRic_geo_M_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               1.0000000               1.0000000 
    ##                 logArea 
    ##               0.6933258

``` r
# Stratified lakes
SD_logSRic_geo_S_model <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[stratified_lakes,])
summary(SD_logSRic_geo_S_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##     BCM     CLM     GLK     HLM     NLK     OTM     SLN     TLN 
    ##  0.3198 -0.1201 -0.4622  0.6912  0.2227 -1.0490 -0.1905  0.5881 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              1.499758   4.710289   0.318    0.766
    ## distance_to_ocean_min_m -0.006280   0.004002  -1.569    0.192
    ## max_depth               -0.002996   0.056849  -0.053    0.961
    ## logArea                  0.126598   0.582703   0.217    0.839
    ## 
    ## Residual standard error: 0.7649 on 4 degrees of freedom
    ## Multiple R-squared:  0.3822, Adjusted R-squared:  -0.0812 
    ## F-statistic: 0.8248 on 3 and 4 DF,  p-value: 0.5448

``` r
anova(SD_logSRic_geo_S_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 1.35465 1.35465  2.3154 0.2027
    ## max_depth                1 0.06535 0.06535  0.1117 0.7550
    ## logArea                  1 0.02762 0.02762  0.0472 0.8386
    ## Residuals                4 2.34024 0.58506

``` r
p_values <- summary(SD_logSRic_geo_S_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               0.7668993               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
SD_SRic_geo_S_model <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[stratified_lakes,])
summary(SD_SRic_geo_S_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##      BCM      CLM      GLK      HLM      NLK      OTM      SLN      TLN 
    ##  1.93145  0.32193 -5.31399  6.78700  0.09536 -7.80458 -0.90951  4.89232 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              8.66169   39.40929   0.220    0.837
    ## distance_to_ocean_min_m -0.05508    0.03349  -1.645    0.175
    ## max_depth               -0.01023    0.47564  -0.022    0.984
    ## logArea                  0.74304    4.87527   0.152    0.886
    ## 
    ## Residual standard error: 6.4 on 4 degrees of freedom
    ## Multiple R-squared:  0.4065, Adjusted R-squared:  -0.03871 
    ## F-statistic: 0.9131 on 3 and 4 DF,  p-value: 0.5102

``` r
anova(SD_SRic_geo_S_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 108.235 108.235  2.6428 0.1793
    ## max_depth                1   2.995   2.995  0.0731 0.8002
    ## logArea                  1   0.951   0.951  0.0232 0.8862
    ## Residuals                4 163.819  40.955

``` r
p_values <- summary(SD_SRic_geo_S_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                 1.00000                 0.70149                 1.00000 
    ##                 logArea 
    ##                 1.00000

### SD alpha z-scores with env and geo correlated variables

``` r
# Log Area
SD_alpha_LA_plot <- ggplot(data = SR_env[surveyed_sites,], mapping = aes(y = log(row_sum), x = logArea, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = SR_env[surveyed_sites,1], size = 5, point.padding = 3) +
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
  labs(x="Log Area (m"^"2"~")", y="Log SRic", colour = "Site type", fill = "Site type", tag = "c")
(SD_alpha_LA_plot <- SD_alpha_LA_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_alpha_LA_plot.jpg", plot = SD_alpha_LA_plot, width = 6.26, height = 6, units = "in")

# Distance from the ocean mean
SD_alpha_D_plot <- ggplot(data = SR_env[surveyed_sites,], mapping = aes(y = row_sum, x = distance_to_ocean_min_m, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = SR_env[surveyed_sites,1], size = 5, point.padding = 3) +
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
  labs(x="Isolation (m)", y="SRic", colour = "Site type", fill = "Site type", tag = "b")
(SD_alpha_D_plot <- SD_alpha_D_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_alpha_D_plot.jpg", plot = SD_alpha_D_plot, width = 6.26, height = 6, units = "in")

# Max depth
SD_alpha_MD_plot <- ggplot(data = SR_env[surveyed_sites,], mapping = aes(y = row_sum, x = max_depth, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = SR_env[surveyed_sites,1], size = 5, point.padding = 3) +
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
  labs(x="Age (m)", y="SRic", colour = "Site type", fill = "Site type", tag = "a")
(SD_alpha_MD_plot <- SD_alpha_MD_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_alpha_MD_plot.jpg", plot = SD_alpha_MD_plot, width = 6.26, height = 6, units = "in")
```

## SD beta diversity

### SD beta diversity distance calculations plus dispersion and PERMANOVA

``` r
### Regular
# Reference
SD_beta_ref_dist<- BAT::beta(presabs_lake, abund = F)
(SD_beta_ref_BD <- betadisper(SD_beta_ref_dist$Btotal, stratification_group_ref))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_ref_dist$Btotal, group =
    ## stratification_group_ref)
    ## 
    ## No. of Positive Eigenvalues: 22
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##  Reference      Ocean      Mixed Stratified 
    ##     0.0000     0.5255     0.5147     0.5340 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 22 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.8509 1.0638 0.5853 0.5356 0.5082 0.4669 0.4567 0.4082

``` r
(SD_beta_ref_AOV <- anova(SD_beta_ref_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value    Pr(>F)    
    ## Groups     3 0.26484 0.088279  14.931 3.107e-05 ***
    ## Residuals 19 0.11233 0.005912                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_ref_THSD <- TukeyHSD(SD_beta_ref_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                             diff        lwr       upr     p adj
    ## Ocean-Reference       0.52552924  0.2919988 0.7590597 0.0000249
    ## Mixed-Reference       0.51472959  0.2854072 0.7440519 0.0000258
    ## Stratified-Reference  0.53404498  0.3047226 0.7633673 0.0000159
    ## Mixed-Ocean          -0.01079966 -0.1275649 0.1059656 0.9936256
    ## Stratified-Ocean      0.00851574 -0.1082495 0.1252810 0.9968402
    ## Stratified-Mixed      0.01931540 -0.0887882 0.1274190 0.9575115

``` r
(SD_beta_ref_PM <- adonis2(SD_beta_ref_dist$Btotal ~ env[,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_ref_dist$Btotal ~ env[, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     3   2.6696 0.30291 2.7521  0.001 ***
    ## Residual 19   6.1435 0.69709                  
    ## Total    22   8.8131 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_ref_PM_pair <- pairwise.adonis(SD_beta_ref_dist$Btotal, env[,19], p.adjust.m = "bonferroni", perm = 999))
```

    ## Set of permutations < 'minperm'. Generating entire set.

    ##                     pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted
    ## 1     Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.001      0.006
    ## 2     Stratified vs Ocean  1 1.3093089 3.913330 0.2459152   0.001      0.006
    ## 3 Stratified vs Reference  1 0.6261283 1.909357 0.2143091   0.112      0.672
    ## 4          Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.057      0.342
    ## 5      Mixed vs Reference  1 0.5979203 1.966332 0.2193017   0.106      0.636
    ## 6      Ocean vs Reference  1 0.5599217 1.628213 0.2456489   0.139      0.834
    ##   sig
    ## 1   *
    ## 2   *
    ## 3    
    ## 4    
    ## 5    
    ## 6

``` r
# Surveyed sites
SD_beta_dist<- BAT::beta(surveyed_sites_lake, abund = F)
(SD_beta_BD <- betadisper(SD_beta_dist$Btotal, stratification_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_dist$Btotal, group = stratification_group)
    ## 
    ## No. of Positive Eigenvalues: 21
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.5255     0.5147     0.5340 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 21 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.8509 1.0570 0.5602 0.5088 0.4687 0.4575 0.4098 0.3594

``` r
(SD_beta_AOV <- anova(SD_beta_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq  Mean Sq F value Pr(>F)
    ## Groups     2 0.001498 0.000749  0.1267 0.8818
    ## Residuals 19 0.112366 0.005914

``` r
(SD_beta_THSD <- TukeyHSD(SD_beta_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                          diff         lwr        upr     p adj
    ## Mixed-Ocean      -0.010799644 -0.11630981 0.09471052 0.9634834
    ## Stratified-Ocean  0.008515751 -0.09699441 0.11402591 0.9771176
    ## Stratified-Mixed  0.019315395 -0.07836803 0.11699882 0.8710688

``` r
(SD_beta_PM <- adonis2(SD_beta_dist$Btotal ~ env[surveyed_sites,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_dist$Btotal ~ env[surveyed_sites, 19], permutations = 999)
    ##          Df SumOfSqs      R2     F Pr(>F)    
    ## Model     2   2.1121 0.25583 3.266  0.001 ***
    ## Residual 19   6.1435 0.74417                 
    ## Total    21   8.2555 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_PM_pair <- pairwise.adonis(SD_beta_dist$Btotal, env[surveyed_sites,19], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.3093089 3.913330 0.2459152   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.057      0.171

``` r
# Without LCN
SD_beta_wo_LCN_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_LCN,], abund = F)
(SD_beta_wo_LCN_BD <- betadisper(SD_beta_wo_LCN_dist$Btotal, stratification_group[-8]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_wo_LCN_dist$Btotal, group =
    ## stratification_group[-8])
    ## 
    ## No. of Positive Eigenvalues: 20
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.4758     0.5147     0.5340 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 20 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.8508 0.9729 0.5584 0.4992 0.4632 0.4280 0.3706 0.3397

``` r
(SD_beta_wo_LCN_AOV <- anova(SD_beta_wo_LCN_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.010500 0.0052498  1.2539 0.3092
    ## Residuals 18 0.075362 0.0041868

``` r
(SD_beta_wo_LCN_THSD <- TukeyHSD(SD_beta_wo_LCN_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                        diff         lwr       upr     p adj
    ## Mixed-Ocean      0.03896759 -0.05517575 0.1331109 0.5522541
    ## Stratified-Ocean 0.05828298 -0.03586035 0.1524263 0.2794657
    ## Stratified-Mixed 0.01931539 -0.06325377 0.1018846 0.8234561

``` r
(SD_beta_wo_LCN_PM <- adonis2(SD_beta_wo_LCN_dist$Btotal ~ env[surveyed_sites_wo_LCN,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_wo_LCN_dist$Btotal ~ env[surveyed_sites_wo_LCN, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   2.2539 0.28733 3.6286  0.001 ***
    ## Residual 18   5.5904 0.71267                  
    ## Total    20   7.8444 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_wo_LCN_PM_pair <- pairwise.adonis(SD_beta_wo_LCN_dist$Btotal, env[surveyed_sites_wo_LCN,19], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.3841204 4.397984 0.2856207   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.6396174 2.135322 0.1625633   0.009      0.027   .

``` r
# Without TLN and HLM
SD_beta_wo_TLN_HLM_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_TLN_HLM,], abund = F)
(SD_beta_wo_TLN_HLM_BD <- betadisper(SD_beta_wo_TLN_HLM_dist$Btotal, stratification_group[-c(5,21)]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_wo_TLN_HLM_dist$Btotal, group =
    ## stratification_group[-c(5, 21)])
    ## 
    ## No. of Positive Eigenvalues: 19
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.5255     0.5147     0.4944 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 19 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.7803 0.9765 0.5390 0.4934 0.4611 0.4078 0.3448 0.3059

``` r
(SD_beta_wo_TLN_HLM_AOV <- anova(SD_beta_wo_TLN_HLM_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.003007 0.0015036  0.2325  0.795
    ## Residuals 17 0.109936 0.0064668

``` r
(SD_beta_wo_TLN_HLM_THSD <- TukeyHSD(SD_beta_wo_TLN_HLM_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff        lwr        upr     p adj
    ## Mixed-Ocean      -0.01079965 -0.1222128 0.10061351 0.9665532
    ## Stratified-Ocean -0.03108653 -0.1501922 0.08801915 0.7839629
    ## Stratified-Mixed -0.02028688 -0.1317000 0.09112628 0.8874562

``` r
(SD_beta_wo_TLN_HLM_PM <- adonis2(SD_beta_wo_TLN_HLM_dist$Btotal ~ env[surveyed_sites_wo_TLN_HLM,19], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_wo_TLN_HLM_dist$Btotal ~ env[surveyed_sites_wo_TLN_HLM, 19], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   2.1465 0.28727 3.4259  0.001 ***
    ## Residual 17   5.3256 0.71273                  
    ## Total    19   7.4721 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_wo_TLN_HLM_PM_pair <- pairwise.adonis(SD_beta_wo_TLN_HLM_dist$Btotal, env[surveyed_sites_wo_TLN_HLM,19], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.4219030 4.731564 0.2827927   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.3260066 4.147587 0.2931657   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.051      0.153

``` r
# Mixed and stratified lakes
SD_beta_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], abund = F)
# Ocean sites and mixed lakes
SD_beta_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites,], abund = F)
# Stratified lakes and ocean sites
SD_beta_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites,], abund = F)
# Stratified lakes
SD_beta_S_dist<- BAT::beta(surveyed_sites_lake[stratified_lakes,], abund = F)


### Environmental
# Surveyed sites
SD_beta_env_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_env,], abund = F)
# Mixed and stratified lakes
SD_beta_env_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], abund = F)
# Ocean sites and mixed lakes
SD_beta_env_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites_env,], abund = F)
# Stratified lakes and ocean sites
SD_beta_env_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites_env,], abund = F)
# Mixed lakes
SD_beta_env_M_dist<- BAT::beta(surveyed_sites_lake[mixed_lakes,], abund = F)


### Geographic
# Surveyed sites
SD_beta_geo_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_geo,], abund = F)
# Mixed and stratified lakes
SD_beta_geo_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes_geo,], abund = F)
# Ocean sites and mixed lakes
SD_beta_geo_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites_geo,], abund = F)
# Stratified lakes and ocean sites
SD_beta_geo_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites,], abund = F)
# Mixed lakes
SD_beta_geo_M_dist<- BAT::beta(surveyed_sites_lake[mixed_lakes_geo,], abund = F)
```

### SD beta q values of dispersion statistics

``` r
# Sites
N <- 22
# Categories
k <- 3
# Sum of squares
SS <- SD_beta_AOV$`Sum Sq`[2]
# Mean of squares
MS <- SD_beta_AOV$`Mean Sq`[2]
# q value
q.value <- qtukey(p = 0.95, nmeans = k, df = N - k)
# Balanced
BSE <- sqrt(MS / 8)
# Unbalanced
USE <- sqrt((MS/2)*((1/6)+(1/8)))

group <- SD_beta_THSD$group

(OM_q <- (abs(group[1,1]))/USE)
```

    ## [1] 0.3677399

``` r
(MS_q <- (abs(group[2,1]))/BSE)
```

    ## [1] 0.3132043

``` r
(SO_q <- (abs(group[3,1]))/USE)
```

    ## [1] 0.6577108

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx
```

### SD beta richness dispersions

- The dissimilarity between sites of the same site type

``` r
# Dispersion within-groups
SD_beta_BD_dist <- SD_beta_BD$distances
SD_beta_BD_dist <- as.data.frame(SD_beta_BD_dist)
SD_beta_BD_dist$X <- row.names(SD_beta_BD_dist)
SD_beta_BD_dist_env <- merge(SD_beta_BD_dist, env[surveyed_sites,], by = "X", sort = F)

SD_beta_BD_dist_env$Stratification <- factor(SD_beta_BD_dist_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

(SD_beta_BD_dist_env_plot <- ggplot(SD_beta_BD_dist_env, aes(x = Stratification, y = SD_beta_BD_dist, color = Stratification, fill = Stratification)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Stratification, fill = Stratification)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.8,
            width = 0.1) +
    geom_text_repel(data = SD_beta_BD_dist_env, label = SD_beta_BD_dist_env$X, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
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
  labs(color = "Stratification", tag = "a"))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20dispersion%20boxplot-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_beta_BD_dist_plot.jpg", SD_beta_BD_dist_env_plot, width = 4.88, height = 6, units = "in")
```

### SD beta NMDS

``` r
### Regular
# Reference
SD_beta_ref_NMDS <- metaMDS(SD_beta_ref_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1049903 
    ## Run 1 stress 0.10835 
    ## Run 2 stress 0.1138985 
    ## Run 3 stress 0.1050144 
    ## ... Procrustes: rmse 0.006660779  max resid 0.0237507 
    ## Run 4 stress 0.1097988 
    ## Run 5 stress 0.1063128 
    ## Run 6 stress 0.1063128 
    ## Run 7 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214944  max resid 0.07899447 
    ## Run 8 stress 0.1050144 
    ## ... Procrustes: rmse 0.006664169  max resid 0.02378272 
    ## Run 9 stress 0.1063128 
    ## Run 10 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214776  max resid 0.07898263 
    ## Run 11 stress 0.1160563 
    ## Run 12 stress 0.1063128 
    ## Run 13 stress 0.1051402 
    ## ... Procrustes: rmse 0.02306978  max resid 0.07847144 
    ## Run 14 stress 0.1051402 
    ## ... Procrustes: rmse 0.02306892  max resid 0.07846916 
    ## Run 15 stress 0.1084964 
    ## Run 16 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307332  max resid 0.07848453 
    ## Run 17 stress 0.1138391 
    ## Run 18 stress 0.1309381 
    ## Run 19 stress 0.1097987 
    ## Run 20 stress 0.1050144 
    ## ... Procrustes: rmse 0.006662835  max resid 0.02376326 
    ## Run 21 stress 0.1049903 
    ## ... New best solution
    ## ... Procrustes: rmse 5.041242e-05  max resid 0.0001357551 
    ## ... Similar to previous best
    ## Run 22 stress 0.1088547 
    ## Run 23 stress 0.1136646 
    ## Run 24 stress 0.1079381 
    ## Run 25 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307298  max resid 0.07847561 
    ## Run 26 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214679  max resid 0.07900182 
    ## Run 27 stress 0.1097987 
    ## Run 28 stress 0.1160549 
    ## Run 29 stress 0.1088547 
    ## Run 30 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214336  max resid 0.07897697 
    ## Run 31 stress 0.1051402 
    ## ... Procrustes: rmse 0.02306705  max resid 0.07845978 
    ## Run 32 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 0.00672745  max resid 0.02574535 
    ## Run 33 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284677  max resid 0.07725608 
    ## Run 34 stress 0.1135596 
    ## Run 35 stress 0.1097988 
    ## Run 36 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284434  max resid 0.0773192 
    ## Run 37 stress 0.1138984 
    ## Run 38 stress 0.1135598 
    ## Run 39 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 2.231828e-05  max resid 5.510739e-05 
    ## ... Similar to previous best
    ## Run 40 stress 0.1058423 
    ## Run 41 stress 0.1063128 
    ## Run 42 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283882  max resid 0.07730488 
    ## Run 43 stress 0.1136647 
    ## Run 44 stress 0.1102722 
    ## Run 45 stress 0.1084965 
    ## Run 46 stress 0.113594 
    ## Run 47 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180166  max resid 0.0768478 
    ## Run 48 stress 0.1050338 
    ## ... Procrustes: rmse 0.00627063  max resid 0.02192583 
    ## Run 49 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179955  max resid 0.07684598 
    ## Run 50 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179482  max resid 0.07682597 
    ## Run 51 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180225  max resid 0.07684902 
    ## Run 52 stress 0.1136647 
    ## Run 53 stress 0.1050338 
    ## ... Procrustes: rmse 0.00625746  max resid 0.02187497 
    ## Run 54 stress 0.1083499 
    ## Run 55 stress 0.1160542 
    ## Run 56 stress 0.1084965 
    ## Run 57 stress 0.1082722 
    ## Run 58 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283904  max resid 0.07729661 
    ## Run 59 stress 0.1082723 
    ## Run 60 stress 0.1160537 
    ## Run 61 stress 0.1097988 
    ## Run 62 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180333  max resid 0.07685975 
    ## Run 63 stress 0.1079381 
    ## Run 64 stress 0.1135598 
    ## Run 65 stress 0.1049903 
    ## ... Procrustes: rmse 0.006789844  max resid 0.02608757 
    ## Run 66 stress 0.1084964 
    ## Run 67 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180038  max resid 0.07683339 
    ## Run 68 stress 0.1050144 
    ## ... Procrustes: rmse 0.00968137  max resid 0.0261364 
    ## Run 69 stress 0.1079383 
    ## Run 70 stress 0.1079382 
    ## Run 71 stress 0.1058423 
    ## Run 72 stress 0.1136646 
    ## Run 73 stress 0.1160554 
    ## Run 74 stress 0.1063128 
    ## Run 75 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283956  max resid 0.07730252 
    ## Run 76 stress 0.113594 
    ## Run 77 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284046  max resid 0.07727938 
    ## Run 78 stress 0.1082722 
    ## Run 79 stress 0.1136646 
    ## Run 80 stress 0.1083995 
    ## Run 81 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180514  max resid 0.0768807 
    ## Run 82 stress 0.1049903 
    ## ... Procrustes: rmse 0.006721998  max resid 0.02577359 
    ## Run 83 stress 0.1165307 
    ## Run 84 stress 0.1083501 
    ## Run 85 stress 0.1082722 
    ## Run 86 stress 0.1063128 
    ## Run 87 stress 0.1138393 
    ## Run 88 stress 0.1083986 
    ## Run 89 stress 0.1051402 
    ## ... Procrustes: rmse 0.0217998  max resid 0.0768182 
    ## Run 90 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 5.166222e-06  max resid 1.079039e-05 
    ## ... Similar to previous best
    ## Run 91 stress 0.1050338 
    ## ... Procrustes: rmse 0.006251407  max resid 0.02183972 
    ## Run 92 stress 0.1097988 
    ## Run 93 stress 0.1135597 
    ## Run 94 stress 0.1082722 
    ## Run 95 stress 0.1063128 
    ## Run 96 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283926  max resid 0.07729398 
    ## Run 97 stress 0.1097988 
    ## Run 98 stress 0.1102722 
    ## Run 99 stress 0.1083499 
    ## Run 100 stress 0.1138396 
    ## Run 101 stress 0.1079382 
    ## Run 102 stress 0.1063128 
    ## Run 103 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180121  max resid 0.0768554 
    ## Run 104 stress 0.110272 
    ## Run 105 stress 0.10835 
    ## Run 106 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180008  max resid 0.0768512 
    ## Run 107 stress 0.1050144 
    ## ... Procrustes: rmse 0.009654864  max resid 0.02603846 
    ## Run 108 stress 0.1063128 
    ## Run 109 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218054  max resid 0.07686872 
    ## Run 110 stress 0.1083999 
    ## Run 111 stress 0.1049903 
    ## ... Procrustes: rmse 0.006744564  max resid 0.02587922 
    ## Run 112 stress 0.113594 
    ## Run 113 stress 0.1050338 
    ## ... Procrustes: rmse 0.006264207  max resid 0.02189879 
    ## Run 114 stress 0.1082722 
    ## Run 115 stress 0.1138397 
    ## Run 116 stress 0.113594 
    ## Run 117 stress 0.1135597 
    ## Run 118 stress 0.1136646 
    ## Run 119 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284225  max resid 0.07730177 
    ## Run 120 stress 0.113594 
    ## Run 121 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218027  max resid 0.07686205 
    ## Run 122 stress 0.1136646 
    ## Run 123 stress 0.1084965 
    ## Run 124 stress 0.1135598 
    ## Run 125 stress 0.1050144 
    ## ... Procrustes: rmse 0.00967222  max resid 0.02615904 
    ## Run 126 stress 0.107938 
    ## Run 127 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180321  max resid 0.07686458 
    ## Run 128 stress 0.1084966 
    ## Run 129 stress 0.1135599 
    ## Run 130 stress 0.1049791 
    ## ... Procrustes: rmse 1.263883e-05  max resid 3.583027e-05 
    ## ... Similar to previous best
    ## Run 131 stress 0.1097989 
    ## Run 132 stress 0.1138392 
    ## Run 133 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179782  max resid 0.07684228 
    ## Run 134 stress 0.1063128 
    ## Run 135 stress 0.1135598 
    ## Run 136 stress 0.1097987 
    ## Run 137 stress 0.107938 
    ## Run 138 stress 0.1136646 
    ## Run 139 stress 0.1050338 
    ## ... Procrustes: rmse 0.0062529  max resid 0.02184952 
    ## Run 140 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284015  max resid 0.07731121 
    ## Run 141 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283835  max resid 0.07729217 
    ## Run 142 stress 0.130938 
    ## Run 143 stress 0.1079381 
    ## Run 144 stress 0.1050338 
    ## ... Procrustes: rmse 0.006259175  max resid 0.02186992 
    ## Run 145 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283839  max resid 0.07728413 
    ## Run 146 stress 0.1063128 
    ## Run 147 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283874  max resid 0.07725267 
    ## Run 148 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283856  max resid 0.07729357 
    ## Run 149 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181299  max resid 0.07689737 
    ## Run 150 stress 0.1174543 
    ## Run 151 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228372  max resid 0.07729928 
    ## Run 152 stress 0.1063128 
    ## Run 153 stress 0.1160542 
    ## Run 154 stress 0.107938 
    ## Run 155 stress 0.1050338 
    ## ... Procrustes: rmse 0.006259205  max resid 0.02187611 
    ## Run 156 stress 0.1083501 
    ## Run 157 stress 0.1063128 
    ## Run 158 stress 0.1138985 
    ## Run 159 stress 0.1050338 
    ## ... Procrustes: rmse 0.006255816  max resid 0.02183905 
    ## Run 160 stress 0.1084964 
    ## Run 161 stress 0.1050144 
    ## ... Procrustes: rmse 0.00965405  max resid 0.02603138 
    ## Run 162 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283851  max resid 0.07729677 
    ## Run 163 stress 0.1063128 
    ## Run 164 stress 0.1051924 
    ## ... Procrustes: rmse 0.0228392  max resid 0.07730291 
    ## Run 165 stress 0.1049903 
    ## ... Procrustes: rmse 0.006723878  max resid 0.02578611 
    ## Run 166 stress 0.1082722 
    ## Run 167 stress 0.1135597 
    ## Run 168 stress 0.1083995 
    ## Run 169 stress 0.110272 
    ## Run 170 stress 0.1088548 
    ## Run 171 stress 0.1165298 
    ## Run 172 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180528  max resid 0.07686926 
    ## Run 173 stress 0.1049903 
    ## ... Procrustes: rmse 0.006730098  max resid 0.02581112 
    ## Run 174 stress 0.1097987 
    ## Run 175 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181164  max resid 0.07689075 
    ## Run 176 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180143  max resid 0.07684032 
    ## Run 177 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180284  max resid 0.07686519 
    ## Run 178 stress 0.110272 
    ## Run 179 stress 0.1049903 
    ## ... Procrustes: rmse 0.006731493  max resid 0.0258198 
    ## Run 180 stress 0.1050338 
    ## ... Procrustes: rmse 0.006253849  max resid 0.02186121 
    ## Run 181 stress 0.1084966 
    ## Run 182 stress 0.1097987 
    ## Run 183 stress 0.1082722 
    ## Run 184 stress 0.1084966 
    ## Run 185 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283896  max resid 0.07731032 
    ## Run 186 stress 0.1049903 
    ## ... Procrustes: rmse 0.00670723  max resid 0.02570827 
    ## Run 187 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284055  max resid 0.07733041 
    ## Run 188 stress 0.1063128 
    ## Run 189 stress 0.1097988 
    ## Run 190 stress 0.1079382 
    ## Run 191 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283993  max resid 0.07728379 
    ## Run 192 stress 0.1135597 
    ## Run 193 stress 0.1082722 
    ## Run 194 stress 0.1088547 
    ## Run 195 stress 0.1049903 
    ## ... Procrustes: rmse 0.006729109  max resid 0.02581637 
    ## Run 196 stress 0.1084964 
    ## Run 197 stress 0.1063128 
    ## Run 198 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284268  max resid 0.07731191 
    ## Run 199 stress 0.1084965 
    ## Run 200 stress 0.1063128 
    ## Run 201 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181354  max resid 0.07689966 
    ## Run 202 stress 0.1050144 
    ## ... Procrustes: rmse 0.009663064  max resid 0.02608245 
    ## Run 203 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180268  max resid 0.07686118 
    ## Run 204 stress 0.113594 
    ## Run 205 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283901  max resid 0.07730671 
    ## Run 206 stress 0.1135942 
    ## Run 207 stress 0.108399 
    ## Run 208 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283979  max resid 0.07729205 
    ## Run 209 stress 0.1050338 
    ## ... Procrustes: rmse 0.006264939  max resid 0.02190182 
    ## Run 210 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284003  max resid 0.07728693 
    ## Run 211 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180086  max resid 0.07685595 
    ## Run 212 stress 0.1138394 
    ## Run 213 stress 0.1063128 
    ## Run 214 stress 0.1135599 
    ## Run 215 stress 0.1050144 
    ## ... Procrustes: rmse 0.009661832  max resid 0.02609227 
    ## Run 216 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284811  max resid 0.07734023 
    ## Run 217 stress 0.1136646 
    ## Run 218 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284  max resid 0.0773097 
    ## Run 219 stress 0.1135941 
    ## Run 220 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284226  max resid 0.07732381 
    ## Run 221 stress 0.1049791 
    ## ... Procrustes: rmse 2.851055e-06  max resid 7.780921e-06 
    ## ... Similar to previous best
    ## Run 222 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283825  max resid 0.0772978 
    ## Run 223 stress 0.1136647 
    ## Run 224 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283799  max resid 0.07730444 
    ## Run 225 stress 0.1165306 
    ## Run 226 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283923  max resid 0.07730425 
    ## Run 227 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283639  max resid 0.07729144 
    ## Run 228 stress 0.107938 
    ## Run 229 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283883  max resid 0.07730958 
    ## Run 230 stress 0.1049903 
    ## ... Procrustes: rmse 0.006730109  max resid 0.02581518 
    ## Run 231 stress 0.1083503 
    ## Run 232 stress 0.1049903 
    ## ... Procrustes: rmse 0.006724151  max resid 0.02578765 
    ## Run 233 stress 0.1082722 
    ## Run 234 stress 0.1063128 
    ## Run 235 stress 0.1049791 
    ## ... Procrustes: rmse 8.166714e-06  max resid 1.738993e-05 
    ## ... Similar to previous best
    ## Run 236 stress 0.113594 
    ## Run 237 stress 0.1079382 
    ## Run 238 stress 0.1088548 
    ## Run 239 stress 0.1084965 
    ## Run 240 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181185  max resid 0.07689431 
    ## Run 241 stress 0.1049791 
    ## ... Procrustes: rmse 7.545705e-06  max resid 2.853272e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.1102723 
    ## Run 243 stress 0.1050338 
    ## ... Procrustes: rmse 0.006253411  max resid 0.02184583 
    ## Run 244 stress 0.1084965 
    ## Run 245 stress 0.1136646 
    ## Run 246 stress 0.1097988 
    ## Run 247 stress 0.1084964 
    ## Run 248 stress 0.1135597 
    ## Run 249 stress 0.1136646 
    ## Run 250 stress 0.1083499 
    ## Run 251 stress 0.10835 
    ## Run 252 stress 0.1136646 
    ## Run 253 stress 0.1083996 
    ## Run 254 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180143  max resid 0.07685685 
    ## Run 255 stress 0.1049903 
    ## ... Procrustes: rmse 0.006737128  max resid 0.02584606 
    ## Run 256 stress 0.1136646 
    ## Run 257 stress 0.1083503 
    ## Run 258 stress 0.1050338 
    ## ... Procrustes: rmse 0.006251712  max resid 0.02185114 
    ## Run 259 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283766  max resid 0.07732118 
    ## Run 260 stress 0.1097987 
    ## Run 261 stress 0.1084964 
    ## Run 262 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180239  max resid 0.07685997 
    ## Run 263 stress 0.1160551 
    ## Run 264 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180279  max resid 0.07686242 
    ## Run 265 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283467  max resid 0.07728588 
    ## Run 266 stress 0.1083501 
    ## Run 267 stress 0.11653 
    ## Run 268 stress 0.1058423 
    ## Run 269 stress 0.1084967 
    ## Run 270 stress 0.1160563 
    ## Run 271 stress 0.1309384 
    ## Run 272 stress 0.1136646 
    ## Run 273 stress 0.1049903 
    ## ... Procrustes: rmse 0.006787992  max resid 0.02608447 
    ## Run 274 stress 0.1135597 
    ## Run 275 stress 0.11356 
    ## Run 276 stress 0.1083995 
    ## Run 277 stress 0.1102722 
    ## Run 278 stress 0.1088548 
    ## Run 279 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180357  max resid 0.07686744 
    ## Run 280 stress 0.1050338 
    ## ... Procrustes: rmse 0.006253716  max resid 0.0218494 
    ## Run 281 stress 0.1050338 
    ## ... Procrustes: rmse 0.006258794  max resid 0.02187452 
    ## Run 282 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284049  max resid 0.07729866 
    ## Run 283 stress 0.1088547 
    ## Run 284 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283995  max resid 0.07730369 
    ## Run 285 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283755  max resid 0.07728033 
    ## Run 286 stress 0.1050338 
    ## ... Procrustes: rmse 0.006267361  max resid 0.0219103 
    ## Run 287 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228423  max resid 0.0772992 
    ## Run 288 stress 0.1050338 
    ## ... Procrustes: rmse 0.00626002  max resid 0.02188093 
    ## Run 289 stress 0.1138394 
    ## Run 290 stress 0.1138986 
    ## Run 291 stress 0.1135598 
    ## Run 292 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180439  max resid 0.07687998 
    ## Run 293 stress 0.1079381 
    ## Run 294 stress 0.113594 
    ## Run 295 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283913  max resid 0.07731091 
    ## Run 296 stress 0.1049791 
    ## ... Procrustes: rmse 1.949235e-05  max resid 5.898326e-05 
    ## ... Similar to previous best
    ## Run 297 stress 0.1082722 
    ## Run 298 stress 0.1082722 
    ## Run 299 stress 0.1088548 
    ## Run 300 stress 0.1082722 
    ## Run 301 stress 0.1135598 
    ## Run 302 stress 0.11356 
    ## Run 303 stress 0.1088547 
    ## Run 304 stress 0.1050338 
    ## ... Procrustes: rmse 0.006253706  max resid 0.02186234 
    ## Run 305 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283811  max resid 0.07728623 
    ## Run 306 stress 0.1049903 
    ## ... Procrustes: rmse 0.006715384  max resid 0.02574749 
    ## Run 307 stress 0.1063128 
    ## Run 308 stress 0.1050144 
    ## ... Procrustes: rmse 0.009666193  max resid 0.02611834 
    ## Run 309 stress 0.1138392 
    ## Run 310 stress 0.1063128 
    ## Run 311 stress 0.1084965 
    ## Run 312 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283894  max resid 0.07730007 
    ## Run 313 stress 0.1088547 
    ## Run 314 stress 0.1084966 
    ## Run 315 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180278  max resid 0.07686151 
    ## Run 316 stress 0.1063128 
    ## Run 317 stress 0.1083501 
    ## Run 318 stress 0.1088547 
    ## Run 319 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284346  max resid 0.07732832 
    ## Run 320 stress 0.1136646 
    ## Run 321 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180442  max resid 0.07686666 
    ## Run 322 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180115  max resid 0.07683839 
    ## Run 323 stress 0.1084965 
    ## Run 324 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180322  max resid 0.07686171 
    ## Run 325 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180235  max resid 0.07685932 
    ## Run 326 stress 0.113594 
    ## Run 327 stress 0.1160579 
    ## Run 328 stress 0.1084965 
    ## Run 329 stress 0.1079381 
    ## Run 330 stress 0.110272 
    ## Run 331 stress 0.1136646 
    ## Run 332 stress 0.1138399 
    ## Run 333 stress 0.1135597 
    ## Run 334 stress 0.1160549 
    ## Run 335 stress 0.1138391 
    ## Run 336 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180576  max resid 0.07687068 
    ## Run 337 stress 0.1050144 
    ## ... Procrustes: rmse 0.009645441  max resid 0.02601247 
    ## Run 338 stress 0.1063128 
    ## Run 339 stress 0.1050144 
    ## ... Procrustes: rmse 0.009663366  max resid 0.02610561 
    ## Run 340 stress 0.1088548 
    ## Run 341 stress 0.1138395 
    ## Run 342 stress 0.113594 
    ## Run 343 stress 0.1160546 
    ## Run 344 stress 0.1050144 
    ## ... Procrustes: rmse 0.009662149  max resid 0.02608956 
    ## Run 345 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283857  max resid 0.07730425 
    ## Run 346 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228393  max resid 0.07730386 
    ## Run 347 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218033  max resid 0.07686646 
    ## Run 348 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180241  max resid 0.07685547 
    ## Run 349 stress 0.1049903 
    ## ... Procrustes: rmse 0.006728462  max resid 0.02581981 
    ## Run 350 stress 0.1174543 
    ## Run 351 stress 0.1097988 
    ## Run 352 stress 0.1084964 
    ## Run 353 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283839  max resid 0.07729786 
    ## Run 354 stress 0.1063128 
    ## Run 355 stress 0.1135598 
    ## Run 356 stress 0.1088547 
    ## Run 357 stress 0.1049903 
    ## ... Procrustes: rmse 0.006724408  max resid 0.0257933 
    ## Run 358 stress 0.1063128 
    ## Run 359 stress 0.1063128 
    ## Run 360 stress 0.116057 
    ## Run 361 stress 0.1082722 
    ## Run 362 stress 0.1050338 
    ## ... Procrustes: rmse 0.006271173  max resid 0.021948 
    ## Run 363 stress 0.1135596 
    ## Run 364 stress 0.1050144 
    ## ... Procrustes: rmse 0.00966003  max resid 0.02607897 
    ## Run 365 stress 0.1050144 
    ## ... Procrustes: rmse 0.009661801  max resid 0.02609353 
    ## Run 366 stress 0.1138985 
    ## Run 367 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283939  max resid 0.07732149 
    ## Run 368 stress 0.1063128 
    ## Run 369 stress 0.1083995 
    ## Run 370 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228357  max resid 0.07726396 
    ## Run 371 stress 0.1049903 
    ## ... Procrustes: rmse 0.006784834  max resid 0.02607045 
    ## Run 372 stress 0.1050144 
    ## ... Procrustes: rmse 0.009659817  max resid 0.02608544 
    ## Run 373 stress 0.113594 
    ## Run 374 stress 0.1051924 
    ## ... Procrustes: rmse 0.0228365  max resid 0.07729159 
    ## Run 375 stress 0.1084965 
    ## Run 376 stress 0.1135599 
    ## Run 377 stress 0.1083499 
    ## Run 378 stress 0.1083998 
    ## Run 379 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180148  max resid 0.07685683 
    ## Run 380 stress 0.1050338 
    ## ... Procrustes: rmse 0.00625698  max resid 0.02186482 
    ## Run 381 stress 0.1058423 
    ## Run 382 stress 0.1088548 
    ## Run 383 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284167  max resid 0.07729144 
    ## Run 384 stress 0.1135598 
    ## Run 385 stress 0.1083992 
    ## Run 386 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283734  max resid 0.07729462 
    ## Run 387 stress 0.1084964 
    ## Run 388 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283882  max resid 0.07729735 
    ## Run 389 stress 0.1050144 
    ## ... Procrustes: rmse 0.009677133  max resid 0.0261864 
    ## Run 390 stress 0.1135597 
    ## Run 391 stress 0.1135598 
    ## Run 392 stress 0.1138986 
    ## Run 393 stress 0.1135597 
    ## Run 394 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180271  max resid 0.07686574 
    ## Run 395 stress 0.1084965 
    ## Run 396 stress 0.1102721 
    ## Run 397 stress 0.1136648 
    ## Run 398 stress 0.1063128 
    ## Run 399 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180006  max resid 0.07684583 
    ## Run 400 stress 0.1079381 
    ## Run 401 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283883  max resid 0.07729708 
    ## Run 402 stress 0.1136646 
    ## Run 403 stress 0.1088547 
    ## Run 404 stress 0.1136646 
    ## Run 405 stress 0.1135596 
    ## Run 406 stress 0.1084966 
    ## Run 407 stress 0.1050144 
    ## ... Procrustes: rmse 0.009661372  max resid 0.0260933 
    ## Run 408 stress 0.1136647 
    ## Run 409 stress 0.1049791 
    ## ... Procrustes: rmse 9.180088e-06  max resid 1.955639e-05 
    ## ... Similar to previous best
    ## Run 410 stress 0.1138398 
    ## Run 411 stress 0.1079381 
    ## Run 412 stress 0.1135597 
    ## Run 413 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179463  max resid 0.07683419 
    ## Run 414 stress 0.1051924 
    ## ... Procrustes: rmse 0.0228417  max resid 0.07726474 
    ## Run 415 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181056  max resid 0.07688711 
    ## Run 416 stress 0.1135598 
    ## Run 417 stress 0.1084003 
    ## Run 418 stress 0.1135597 
    ## Run 419 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180866  max resid 0.07688134 
    ## Run 420 stress 0.1136646 
    ## Run 421 stress 0.1063128 
    ## Run 422 stress 0.1083502 
    ## Run 423 stress 0.1082722 
    ## Run 424 stress 0.113594 
    ## Run 425 stress 0.1138985 
    ## Run 426 stress 0.1138392 
    ## Run 427 stress 0.1050144 
    ## ... Procrustes: rmse 0.009660713  max resid 0.02608449 
    ## Run 428 stress 0.1136646 
    ## Run 429 stress 0.1088547 
    ## Run 430 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180623  max resid 0.07687229 
    ## Run 431 stress 0.1084965 
    ## Run 432 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284005  max resid 0.07730717 
    ## Run 433 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180169  max resid 0.07685796 
    ## Run 434 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283829  max resid 0.07731227 
    ## Run 435 stress 0.1136646 
    ## Run 436 stress 0.1063128 
    ## Run 437 stress 0.1049903 
    ## ... Procrustes: rmse 0.006786484  max resid 0.02607791 
    ## Run 438 stress 0.1135941 
    ## Run 439 stress 0.1136648 
    ## Run 440 stress 0.1082722 
    ## Run 441 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180983  max resid 0.07688889 
    ## Run 442 stress 0.107938 
    ## Run 443 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218081  max resid 0.07687868 
    ## Run 444 stress 0.1138393 
    ## Run 445 stress 0.1084964 
    ## Run 446 stress 0.1136646 
    ## Run 447 stress 0.1049791 
    ## ... Procrustes: rmse 4.059476e-06  max resid 8.442318e-06 
    ## ... Similar to previous best
    ## Run 448 stress 0.1049903 
    ## ... Procrustes: rmse 0.006743849  max resid 0.02587742 
    ## Run 449 stress 0.1063128 
    ## Run 450 stress 0.1136647 
    ## Run 451 stress 0.1135597 
    ## Run 452 stress 0.113594 
    ## Run 453 stress 0.1063128 
    ## Run 454 stress 0.1082722 
    ## Run 455 stress 0.1083501 
    ## Run 456 stress 0.1160549 
    ## Run 457 stress 0.1097989 
    ## Run 458 stress 0.1084965 
    ## Run 459 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180134  max resid 0.0768568 
    ## Run 460 stress 0.1138392 
    ## Run 461 stress 0.1049791 
    ## ... Procrustes: rmse 4.990716e-06  max resid 1.459603e-05 
    ## ... Similar to previous best
    ## Run 462 stress 0.1082722 
    ## Run 463 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283954  max resid 0.07731462 
    ## Run 464 stress 0.1063128 
    ## Run 465 stress 0.1084964 
    ## Run 466 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228384  max resid 0.07730561 
    ## Run 467 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283945  max resid 0.07730607 
    ## Run 468 stress 0.1049791 
    ## ... Procrustes: rmse 2.865123e-06  max resid 7.655268e-06 
    ## ... Similar to previous best
    ## Run 469 stress 0.1079383 
    ## Run 470 stress 0.10835 
    ## Run 471 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180526  max resid 0.07686837 
    ## Run 472 stress 0.1063128 
    ## Run 473 stress 0.1082722 
    ## Run 474 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284005  max resid 0.0772822 
    ## Run 475 stress 0.1083503 
    ## Run 476 stress 0.1160582 
    ## Run 477 stress 0.1136648 
    ## Run 478 stress 0.1063128 
    ## Run 479 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180247  max resid 0.07685919 
    ## Run 480 stress 0.1050144 
    ## ... Procrustes: rmse 0.009673728  max resid 0.02614735 
    ## Run 481 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284108  max resid 0.0772482 
    ## Run 482 stress 0.1136646 
    ## Run 483 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283888  max resid 0.07728605 
    ## Run 484 stress 0.1049903 
    ## ... Procrustes: rmse 0.006785162  max resid 0.02607111 
    ## Run 485 stress 0.1138397 
    ## Run 486 stress 0.1083995 
    ## Run 487 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283738  max resid 0.0773036 
    ## Run 488 stress 0.1063128 
    ## Run 489 stress 0.1084967 
    ## Run 490 stress 0.110272 
    ## Run 491 stress 0.1058423 
    ## Run 492 stress 0.1084964 
    ## Run 493 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180281  max resid 0.07686263 
    ## Run 494 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283984  max resid 0.07729028 
    ## Run 495 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283972  max resid 0.07730656 
    ## Run 496 stress 0.1050338 
    ## ... Procrustes: rmse 0.006256759  max resid 0.02186573 
    ## Run 497 stress 0.1138393 
    ## Run 498 stress 0.10835 
    ## Run 499 stress 0.1049791 
    ## ... Procrustes: rmse 1.915058e-05  max resid 5.218063e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.113594 
    ## *** Best solution repeated 11 times

``` r
# Surveyed sites 
SD_beta_NMDS <- metaMDS(SD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.1074432 
    ## Run 2 stress 0.1092947 
    ## Run 3 stress 0.1092128 
    ## Run 4 stress 0.08926081 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04356292  max resid 0.1169066 
    ## Run 5 stress 0.08938538 
    ## ... Procrustes: rmse 0.0118372  max resid 0.03954363 
    ## Run 6 stress 0.09228502 
    ## Run 7 stress 0.0893854 
    ## ... Procrustes: rmse 0.0118613  max resid 0.03958843 
    ## Run 8 stress 0.08938547 
    ## ... Procrustes: rmse 0.01190719  max resid 0.03951706 
    ## Run 9 stress 0.1067502 
    ## Run 10 stress 0.08938968 
    ## ... Procrustes: rmse 0.03587163  max resid 0.118166 
    ## Run 11 stress 0.1056899 
    ## Run 12 stress 0.1065226 
    ## Run 13 stress 0.1092092 
    ## Run 14 stress 0.08926094 
    ## ... Procrustes: rmse 0.0006381279  max resid 0.002009707 
    ## ... Similar to previous best
    ## Run 15 stress 0.1108454 
    ## Run 16 stress 0.08926075 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002125167  max resid 0.0006852753 
    ## ... Similar to previous best
    ## Run 17 stress 0.0893898 
    ## ... Procrustes: rmse 0.03588569  max resid 0.1182016 
    ## Run 18 stress 0.09039142 
    ## Run 19 stress 0.1074431 
    ## Run 20 stress 0.105265 
    ## Run 21 stress 0.1076303 
    ## Run 22 stress 0.1065378 
    ## Run 23 stress 0.1052652 
    ## Run 24 stress 0.0892608 
    ## ... Procrustes: rmse 6.908828e-05  max resid 0.000177892 
    ## ... Similar to previous best
    ## Run 25 stress 0.109221 
    ## Run 26 stress 0.08926075 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003486935  max resid 0.000885618 
    ## ... Similar to previous best
    ## Run 27 stress 0.08926085 
    ## ... Procrustes: rmse 8.251579e-05  max resid 0.0002760663 
    ## ... Similar to previous best
    ## Run 28 stress 0.1052651 
    ## Run 29 stress 0.09503456 
    ## Run 30 stress 0.0918088 
    ## Run 31 stress 0.08938565 
    ## ... Procrustes: rmse 0.01188112  max resid 0.0396626 
    ## Run 32 stress 0.1087863 
    ## Run 33 stress 0.09130086 
    ## Run 34 stress 0.09099546 
    ## Run 35 stress 0.1065382 
    ## Run 36 stress 0.1063192 
    ## Run 37 stress 0.08926083 
    ## ... Procrustes: rmse 0.0005046571  max resid 0.001493585 
    ## ... Similar to previous best
    ## Run 38 stress 0.1096131 
    ## Run 39 stress 0.09503445 
    ## Run 40 stress 0.08938541 
    ## ... Procrustes: rmse 0.01183821  max resid 0.03971627 
    ## Run 41 stress 0.08951718 
    ## ... Procrustes: rmse 0.03510671  max resid 0.1157307 
    ## Run 42 stress 0.1071321 
    ## Run 43 stress 0.1063198 
    ## Run 44 stress 0.08926109 
    ## ... Procrustes: rmse 0.0002169146  max resid 0.0007132469 
    ## ... Similar to previous best
    ## Run 45 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002160164  max resid 0.0006857899 
    ## ... Similar to previous best
    ## Run 46 stress 0.09087359 
    ## Run 47 stress 0.09228507 
    ## Run 48 stress 0.1064189 
    ## Run 49 stress 0.08926068 
    ## ... Procrustes: rmse 0.0001145347  max resid 0.0003244718 
    ## ... Similar to previous best
    ## Run 50 stress 0.0913009 
    ## Run 51 stress 0.09130085 
    ## Run 52 stress 0.08938557 
    ## ... Procrustes: rmse 0.01178693  max resid 0.03974529 
    ## Run 53 stress 0.08946334 
    ## ... Procrustes: rmse 0.03752059  max resid 0.1177802 
    ## Run 54 stress 0.09021169 
    ## Run 55 stress 0.08938539 
    ## ... Procrustes: rmse 0.01185932  max resid 0.03967211 
    ## Run 56 stress 0.09180873 
    ## Run 57 stress 0.08946663 
    ## ... Procrustes: rmse 0.03320551  max resid 0.1163162 
    ## Run 58 stress 0.09228494 
    ## Run 59 stress 0.09087353 
    ## Run 60 stress 0.08946664 
    ## ... Procrustes: rmse 0.03318178  max resid 0.1162938 
    ## Run 61 stress 0.08951718 
    ## ... Procrustes: rmse 0.03506387  max resid 0.1156785 
    ## Run 62 stress 0.106759 
    ## Run 63 stress 0.09503419 
    ## Run 64 stress 0.1066201 
    ## Run 65 stress 0.08938548 
    ## ... Procrustes: rmse 0.01179174  max resid 0.03966886 
    ## Run 66 stress 0.08926078 
    ## ... Procrustes: rmse 0.0002381731  max resid 0.0008009997 
    ## ... Similar to previous best
    ## Run 67 stress 0.09039133 
    ## Run 68 stress 0.1071321 
    ## Run 69 stress 0.08926069 
    ## ... Procrustes: rmse 0.0001252308  max resid 0.000390432 
    ## ... Similar to previous best
    ## Run 70 stress 0.08938963 
    ## ... Procrustes: rmse 0.03596608  max resid 0.1183169 
    ## Run 71 stress 0.08938962 
    ## ... Procrustes: rmse 0.0359405  max resid 0.1182715 
    ## Run 72 stress 0.08926073 
    ## ... Procrustes: rmse 0.0001805345  max resid 0.0006049593 
    ## ... Similar to previous best
    ## Run 73 stress 0.08938971 
    ## ... Procrustes: rmse 0.03592156  max resid 0.1182392 
    ## Run 74 stress 0.1071321 
    ## Run 75 stress 0.1071322 
    ## Run 76 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595502  max resid 0.1182948 
    ## Run 77 stress 0.0961281 
    ## Run 78 stress 0.1064192 
    ## Run 79 stress 0.08926065 
    ## ... Procrustes: rmse 7.733673e-05  max resid 0.0001804148 
    ## ... Similar to previous best
    ## Run 80 stress 0.08938541 
    ## ... Procrustes: rmse 0.01186086  max resid 0.03965305 
    ## Run 81 stress 0.1092127 
    ## Run 82 stress 0.08938971 
    ## ... Procrustes: rmse 0.0359228  max resid 0.1182461 
    ## Run 83 stress 0.08926088 
    ## ... Procrustes: rmse 0.0003103638  max resid 0.001009088 
    ## ... Similar to previous best
    ## Run 84 stress 0.08938971 
    ## ... Procrustes: rmse 0.03591934  max resid 0.1182367 
    ## Run 85 stress 0.09503438 
    ## Run 86 stress 0.09503436 
    ## Run 87 stress 0.08926069 
    ## ... Procrustes: rmse 0.000146293  max resid 0.0003555153 
    ## ... Similar to previous best
    ## Run 88 stress 0.08946333 
    ## ... Procrustes: rmse 0.03754354  max resid 0.1178118 
    ## Run 89 stress 0.09503418 
    ## Run 90 stress 0.1087575 
    ## Run 91 stress 0.1079016 
    ## Run 92 stress 0.1108319 
    ## Run 93 stress 0.1076302 
    ## Run 94 stress 0.09109129 
    ## Run 95 stress 0.08938968 
    ## ... Procrustes: rmse 0.03592533  max resid 0.1182471 
    ## Run 96 stress 0.1071322 
    ## Run 97 stress 0.1091237 
    ## Run 98 stress 0.1091882 
    ## Run 99 stress 0.08946344 
    ## ... Procrustes: rmse 0.03747701  max resid 0.1177195 
    ## Run 100 stress 0.08926067 
    ## ... Procrustes: rmse 0.0001137576  max resid 0.0003884165 
    ## ... Similar to previous best
    ## Run 101 stress 0.1087863 
    ## Run 102 stress 0.1092363 
    ## Run 103 stress 0.1087583 
    ## Run 104 stress 0.08926068 
    ## ... Procrustes: rmse 0.000119075  max resid 0.000389769 
    ## ... Similar to previous best
    ## Run 105 stress 0.1092095 
    ## Run 106 stress 0.0893855 
    ## ... Procrustes: rmse 0.01179339  max resid 0.03971527 
    ## Run 107 stress 0.1056902 
    ## Run 108 stress 0.09130096 
    ## Run 109 stress 0.09503424 
    ## Run 110 stress 0.1091504 
    ## Run 111 stress 0.1087579 
    ## Run 112 stress 0.0893897 
    ## ... Procrustes: rmse 0.03592143  max resid 0.1182409 
    ## Run 113 stress 0.09021166 
    ## Run 114 stress 0.08938543 
    ## ... Procrustes: rmse 0.01178986  max resid 0.03969839 
    ## Run 115 stress 0.1091423 
    ## Run 116 stress 0.1071321 
    ## Run 117 stress 0.08946651 
    ## ... Procrustes: rmse 0.03322706  max resid 0.1163612 
    ## Run 118 stress 0.1056906 
    ## Run 119 stress 0.1071321 
    ## Run 120 stress 0.08926087 
    ## ... Procrustes: rmse 0.0003082703  max resid 0.0009697544 
    ## ... Similar to previous best
    ## Run 121 stress 0.09039132 
    ## Run 122 stress 0.1103708 
    ## Run 123 stress 0.1075421 
    ## Run 124 stress 0.09109108 
    ## Run 125 stress 0.08938539 
    ## ... Procrustes: rmse 0.01180561  max resid 0.03965774 
    ## Run 126 stress 0.09228482 
    ## Run 127 stress 0.09503427 
    ## Run 128 stress 0.08938544 
    ## ... Procrustes: rmse 0.01184001  max resid 0.03970396 
    ## Run 129 stress 0.08926087 
    ## ... Procrustes: rmse 0.0002863011  max resid 0.0009117649 
    ## ... Similar to previous best
    ## Run 130 stress 0.09228514 
    ## Run 131 stress 0.09044623 
    ## Run 132 stress 0.08926078 
    ## ... Procrustes: rmse 0.000200489  max resid 0.0006293495 
    ## ... Similar to previous best
    ## Run 133 stress 0.08926095 
    ## ... Procrustes: rmse 0.0003506454  max resid 0.001146728 
    ## ... Similar to previous best
    ## Run 134 stress 0.0977559 
    ## Run 135 stress 0.1056904 
    ## Run 136 stress 0.08938962 
    ## ... Procrustes: rmse 0.03596515  max resid 0.118312 
    ## Run 137 stress 0.10675 
    ## Run 138 stress 0.08938961 
    ## ... Procrustes: rmse 0.0359627  max resid 0.1183037 
    ## Run 139 stress 0.1052647 
    ## Run 140 stress 0.0904461 
    ## Run 141 stress 0.1075423 
    ## Run 142 stress 0.0894668 
    ## ... Procrustes: rmse 0.03317703  max resid 0.116287 
    ## Run 143 stress 0.1067379 
    ## Run 144 stress 0.09039135 
    ## Run 145 stress 0.08946331 
    ## ... Procrustes: rmse 0.03749812  max resid 0.1177446 
    ## Run 146 stress 0.09109114 
    ## Run 147 stress 0.08951721 
    ## ... Procrustes: rmse 0.03507565  max resid 0.1156984 
    ## Run 148 stress 0.09503464 
    ## Run 149 stress 0.1063195 
    ## Run 150 stress 0.09503436 
    ## Run 151 stress 0.0894633 
    ## ... Procrustes: rmse 0.03751043  max resid 0.1177634 
    ## Run 152 stress 0.1067366 
    ## Run 153 stress 0.0950342 
    ## Run 154 stress 0.09612827 
    ## Run 155 stress 0.1065235 
    ## Run 156 stress 0.09039143 
    ## Run 157 stress 0.09018709 
    ## Run 158 stress 0.09592174 
    ## Run 159 stress 0.0908861 
    ## Run 160 stress 0.1056905 
    ## Run 161 stress 0.08938552 
    ## ... Procrustes: rmse 0.01178102  max resid 0.03972579 
    ## Run 162 stress 0.1061297 
    ## Run 163 stress 0.09039134 
    ## Run 164 stress 0.1074435 
    ## Run 165 stress 0.08926076 
    ## ... Procrustes: rmse 0.0002533862  max resid 0.0007596354 
    ## ... Similar to previous best
    ## Run 166 stress 0.110832 
    ## Run 167 stress 0.08938543 
    ## ... Procrustes: rmse 0.01184709  max resid 0.03962934 
    ## Run 168 stress 0.1056894 
    ## Run 169 stress 0.1108632 
    ## Run 170 stress 0.08926067 
    ## ... Procrustes: rmse 0.0001148057  max resid 0.0003803329 
    ## ... Similar to previous best
    ## Run 171 stress 0.1103726 
    ## Run 172 stress 0.08926071 
    ## ... Procrustes: rmse 0.0001709952  max resid 0.0004262393 
    ## ... Similar to previous best
    ## Run 173 stress 0.08938539 
    ## ... Procrustes: rmse 0.01183523  max resid 0.03970324 
    ## Run 174 stress 0.1092091 
    ## Run 175 stress 0.09503446 
    ## Run 176 stress 0.1063194 
    ## Run 177 stress 0.1052649 
    ## Run 178 stress 0.0892607 
    ## ... Procrustes: rmse 0.0001664556  max resid 0.0005397584 
    ## ... Similar to previous best
    ## Run 179 stress 0.08926068 
    ## ... Procrustes: rmse 0.0001055186  max resid 0.0002796113 
    ## ... Similar to previous best
    ## Run 180 stress 0.1091824 
    ## Run 181 stress 0.1071321 
    ## Run 182 stress 0.0893897 
    ## ... Procrustes: rmse 0.03597682  max resid 0.1183416 
    ## Run 183 stress 0.09180888 
    ## Run 184 stress 0.09592185 
    ## Run 185 stress 0.1086565 
    ## Run 186 stress 0.09087339 
    ## Run 187 stress 0.1087584 
    ## Run 188 stress 0.08926072 
    ## ... Procrustes: rmse 0.0002134183  max resid 0.0006871266 
    ## ... Similar to previous best
    ## Run 189 stress 0.1087575 
    ## Run 190 stress 0.08926096 
    ## ... Procrustes: rmse 0.0003644461  max resid 0.001190095 
    ## ... Similar to previous best
    ## Run 191 stress 0.08926083 
    ## ... Procrustes: rmse 0.0002794266  max resid 0.0008901981 
    ## ... Similar to previous best
    ## Run 192 stress 0.08938547 
    ## ... Procrustes: rmse 0.01186359  max resid 0.03975705 
    ## Run 193 stress 0.09109109 
    ## Run 194 stress 0.09109109 
    ## Run 195 stress 0.1092091 
    ## Run 196 stress 0.1111389 
    ## Run 197 stress 0.08926096 
    ## ... Procrustes: rmse 0.0003622581  max resid 0.001135443 
    ## ... Similar to previous best
    ## Run 198 stress 0.1074434 
    ## Run 199 stress 0.08951719 
    ## ... Procrustes: rmse 0.03504268  max resid 0.1156219 
    ## Run 200 stress 0.1060434 
    ## Run 201 stress 0.1087573 
    ## Run 202 stress 0.09021165 
    ## Run 203 stress 0.1087866 
    ## Run 204 stress 0.1086566 
    ## Run 205 stress 0.1056902 
    ## Run 206 stress 0.1056901 
    ## Run 207 stress 0.09018714 
    ## Run 208 stress 0.08946337 
    ## ... Procrustes: rmse 0.03748969  max resid 0.117731 
    ## Run 209 stress 0.1060432 
    ## Run 210 stress 0.09018712 
    ## Run 211 stress 0.08938552 
    ## ... Procrustes: rmse 0.01177012  max resid 0.03971208 
    ## Run 212 stress 0.1071321 
    ## Run 213 stress 0.08946657 
    ## ... Procrustes: rmse 0.03319476  max resid 0.1163076 
    ## Run 214 stress 0.1052646 
    ## Run 215 stress 0.1092094 
    ## Run 216 stress 0.09039139 
    ## Run 217 stress 0.1079015 
    ## Run 218 stress 0.08926095 
    ## ... Procrustes: rmse 0.0003509572  max resid 0.001103609 
    ## ... Similar to previous best
    ## Run 219 stress 0.08926119 
    ## ... Procrustes: rmse 0.0004802194  max resid 0.00149856 
    ## ... Similar to previous best
    ## Run 220 stress 0.08926065 
    ## ... Procrustes: rmse 8.787195e-05  max resid 0.0002917444 
    ## ... Similar to previous best
    ## Run 221 stress 0.09018707 
    ## Run 222 stress 0.1092097 
    ## Run 223 stress 0.08938552 
    ## ... Procrustes: rmse 0.01191647  max resid 0.03971264 
    ## Run 224 stress 0.1052646 
    ## Run 225 stress 0.1063184 
    ## Run 226 stress 0.09044616 
    ## Run 227 stress 0.1074435 
    ## Run 228 stress 0.10675 
    ## Run 229 stress 0.08926075 
    ## ... Procrustes: rmse 0.0001942414  max resid 0.00064676 
    ## ... Similar to previous best
    ## Run 230 stress 0.08926079 
    ## ... Procrustes: rmse 0.0001997842  max resid 0.0005419935 
    ## ... Similar to previous best
    ## Run 231 stress 0.09178295 
    ## Run 232 stress 0.0903913 
    ## Run 233 stress 0.0913009 
    ## Run 234 stress 0.1087865 
    ## Run 235 stress 0.08938965 
    ## ... Procrustes: rmse 0.03593272  max resid 0.1182594 
    ## Run 236 stress 0.08926063 
    ## ... Procrustes: rmse 2.491401e-05  max resid 5.890453e-05 
    ## ... Similar to previous best
    ## Run 237 stress 0.1052651 
    ## Run 238 stress 0.09503414 
    ## Run 239 stress 0.08926093 
    ## ... Procrustes: rmse 0.0003414242  max resid 0.001111124 
    ## ... Similar to previous best
    ## Run 240 stress 0.09099526 
    ## Run 241 stress 0.08926072 
    ## ... Procrustes: rmse 0.0001598881  max resid 0.0005353202 
    ## ... Similar to previous best
    ## Run 242 stress 0.08926108 
    ## ... Procrustes: rmse 0.000425753  max resid 0.001398692 
    ## ... Similar to previous best
    ## Run 243 stress 0.1087569 
    ## Run 244 stress 0.08938542 
    ## ... Procrustes: rmse 0.01186737  max resid 0.03965138 
    ## Run 245 stress 0.08938544 
    ## ... Procrustes: rmse 0.01178423  max resid 0.03963798 
    ## Run 246 stress 0.1056902 
    ## Run 247 stress 0.08938964 
    ## ... Procrustes: rmse 0.03596999  max resid 0.1183202 
    ## Run 248 stress 0.1064187 
    ## Run 249 stress 0.0892607 
    ## ... Procrustes: rmse 0.0001449202  max resid 0.0004897291 
    ## ... Similar to previous best
    ## Run 250 stress 0.1071321 
    ## Run 251 stress 0.1066206 
    ## Run 252 stress 0.1064429 
    ## Run 253 stress 0.08938553 
    ## ... Procrustes: rmse 0.01177663  max resid 0.03972586 
    ## Run 254 stress 0.1104468 
    ## Run 255 stress 0.08951731 
    ## ... Procrustes: rmse 0.03501784  max resid 0.1155866 
    ## Run 256 stress 0.08938537 
    ## ... Procrustes: rmse 0.0118307  max resid 0.03968727 
    ## Run 257 stress 0.1074433 
    ## Run 258 stress 0.1075421 
    ## Run 259 stress 0.08946653 
    ## ... Procrustes: rmse 0.03323333  max resid 0.1163836 
    ## Run 260 stress 0.08938556 
    ## ... Procrustes: rmse 0.01184559  max resid 0.03957196 
    ## Run 261 stress 0.1056904 
    ## Run 262 stress 0.1092559 
    ## Run 263 stress 0.1092363 
    ## Run 264 stress 0.09018707 
    ## Run 265 stress 0.0892608 
    ## ... Procrustes: rmse 0.0002868751  max resid 0.0008728202 
    ## ... Similar to previous best
    ## Run 266 stress 0.0950346 
    ## Run 267 stress 0.1071321 
    ## Run 268 stress 0.08938546 
    ## ... Procrustes: rmse 0.01186541  max resid 0.03974175 
    ## Run 269 stress 0.09099544 
    ## Run 270 stress 0.1052652 
    ## Run 271 stress 0.08946651 
    ## ... Procrustes: rmse 0.03320913  max resid 0.1163294 
    ## Run 272 stress 0.1052647 
    ## Run 273 stress 0.09592132 
    ## Run 274 stress 0.09262368 
    ## Run 275 stress 0.09592176 
    ## Run 276 stress 0.08926066 
    ## ... Procrustes: rmse 0.0001175163  max resid 0.000259629 
    ## ... Similar to previous best
    ## Run 277 stress 0.08926088 
    ## ... Procrustes: rmse 0.0002976286  max resid 0.0009635915 
    ## ... Similar to previous best
    ## Run 278 stress 0.08926073 
    ## ... Procrustes: rmse 0.0001927845  max resid 0.0004801767 
    ## ... Similar to previous best
    ## Run 279 stress 0.08938542 
    ## ... Procrustes: rmse 0.01179147  max resid 0.03964383 
    ## Run 280 stress 0.1091824 
    ## Run 281 stress 0.1074432 
    ## Run 282 stress 0.1061301 
    ## Run 283 stress 0.08951724 
    ## ... Procrustes: rmse 0.03503522  max resid 0.1156204 
    ## Run 284 stress 0.1064191 
    ## Run 285 stress 0.08946655 
    ## ... Procrustes: rmse 0.03319985  max resid 0.1163127 
    ## Run 286 stress 0.09021165 
    ## Run 287 stress 0.08946328 
    ## ... Procrustes: rmse 0.03751184  max resid 0.117764 
    ## Run 288 stress 0.1071321 
    ## Run 289 stress 0.08926087 
    ## ... Procrustes: rmse 0.000307363  max resid 0.0009676253 
    ## ... Similar to previous best
    ## Run 290 stress 0.1092129 
    ## Run 291 stress 0.08926064 
    ## ... Procrustes: rmse 8.572634e-05  max resid 0.0002093765 
    ## ... Similar to previous best
    ## Run 292 stress 0.09503425 
    ## Run 293 stress 0.1075421 
    ## Run 294 stress 0.1079015 
    ## Run 295 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595492  max resid 0.1182851 
    ## Run 296 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595556  max resid 0.1182901 
    ## Run 297 stress 0.08926072 
    ## ... Procrustes: rmse 0.0002030603  max resid 0.000606089 
    ## ... Similar to previous best
    ## Run 298 stress 0.09099549 
    ## Run 299 stress 0.1067578 
    ## Run 300 stress 0.1092091 
    ## Run 301 stress 0.08946659 
    ## ... Procrustes: rmse 0.03318873  max resid 0.1162931 
    ## Run 302 stress 0.09099594 
    ## Run 303 stress 0.08946651 
    ## ... Procrustes: rmse 0.0332127  max resid 0.1163359 
    ## Run 304 stress 0.1110866 
    ## Run 305 stress 0.106043 
    ## Run 306 stress 0.09087364 
    ## Run 307 stress 0.1101453 
    ## Run 308 stress 0.0910911 
    ## Run 309 stress 0.1067374 
    ## Run 310 stress 0.1087864 
    ## Run 311 stress 0.09099553 
    ## Run 312 stress 0.09503445 
    ## Run 313 stress 0.08938966 
    ## ... Procrustes: rmse 0.03593072  max resid 0.1182573 
    ## Run 314 stress 0.08938964 
    ## ... Procrustes: rmse 0.03593768  max resid 0.1182742 
    ## Run 315 stress 0.1085134 
    ## Run 316 stress 0.08938973 
    ## ... Procrustes: rmse 0.03591626  max resid 0.1182319 
    ## Run 317 stress 0.08946332 
    ## ... Procrustes: rmse 0.0374926  max resid 0.1177428 
    ## Run 318 stress 0.09109108 
    ## Run 319 stress 0.1061301 
    ## Run 320 stress 0.08946328 
    ## ... Procrustes: rmse 0.03752395  max resid 0.1177875 
    ## Run 321 stress 0.1092127 
    ## Run 322 stress 0.09178298 
    ## Run 323 stress 0.08926083 
    ## ... Procrustes: rmse 0.0002801712  max resid 0.0008768491 
    ## ... Similar to previous best
    ## Run 324 stress 0.09018713 
    ## Run 325 stress 0.09503417 
    ## Run 326 stress 0.08938542 
    ## ... Procrustes: rmse 0.01184889  max resid 0.03972647 
    ## Run 327 stress 0.1061298 
    ## Run 328 stress 0.1074433 
    ## Run 329 stress 0.1091237 
    ## Run 330 stress 0.08938977 
    ## ... Procrustes: rmse 0.03591014  max resid 0.1182257 
    ## Run 331 stress 0.108064 
    ## Run 332 stress 0.1074431 
    ## Run 333 stress 0.1071321 
    ## Run 334 stress 0.1075422 
    ## Run 335 stress 0.1052647 
    ## Run 336 stress 0.09039138 
    ## Run 337 stress 0.1091827 
    ## Run 338 stress 0.09109109 
    ## Run 339 stress 0.08938977 
    ## ... Procrustes: rmse 0.03591907  max resid 0.1182405 
    ## Run 340 stress 0.08926081 
    ## ... Procrustes: rmse 0.0002630062  max resid 0.0008241587 
    ## ... Similar to previous best
    ## Run 341 stress 0.08938542 
    ## ... Procrustes: rmse 0.01180302  max resid 0.03970115 
    ## Run 342 stress 0.1067487 
    ## Run 343 stress 0.08926076 
    ## ... Procrustes: rmse 0.0002491873  max resid 0.0007445184 
    ## ... Similar to previous best
    ## Run 344 stress 0.08938965 
    ## ... Procrustes: rmse 0.03596697  max resid 0.1183224 
    ## Run 345 stress 0.08926072 
    ## ... Procrustes: rmse 0.0001825687  max resid 0.0004471603 
    ## ... Similar to previous best
    ## Run 346 stress 0.1074433 
    ## Run 347 stress 0.08926065 
    ## ... Procrustes: rmse 8.610295e-05  max resid 0.0002251721 
    ## ... Similar to previous best
    ## Run 348 stress 0.1060431 
    ## Run 349 stress 0.09228483 
    ## Run 350 stress 0.1087864 
    ## Run 351 stress 0.09228485 
    ## Run 352 stress 0.09088601 
    ## Run 353 stress 0.1074433 
    ## Run 354 stress 0.08938975 
    ## ... Procrustes: rmse 0.03591344  max resid 0.1182304 
    ## Run 355 stress 0.08938966 
    ## ... Procrustes: rmse 0.0359693  max resid 0.1183217 
    ## Run 356 stress 0.09087349 
    ## Run 357 stress 0.08926068 
    ## ... Procrustes: rmse 0.0001297591  max resid 0.0002768609 
    ## ... Similar to previous best
    ## Run 358 stress 0.08938539 
    ## ... Procrustes: rmse 0.0118238  max resid 0.03970784 
    ## Run 359 stress 0.09228504 
    ## Run 360 stress 0.0908735 
    ## Run 361 stress 0.09503421 
    ## Run 362 stress 0.1065356 
    ## Run 363 stress 0.1052647 
    ## Run 364 stress 0.0908734 
    ## Run 365 stress 0.1074431 
    ## Run 366 stress 0.09088601 
    ## Run 367 stress 0.08926111 
    ## ... Procrustes: rmse 0.0003366255  max resid 0.001053031 
    ## ... Similar to previous best
    ## Run 368 stress 0.08946673 
    ## ... Procrustes: rmse 0.03317097  max resid 0.1162766 
    ## Run 369 stress 0.08926074 
    ## ... Procrustes: rmse 0.0001996398  max resid 0.0006728563 
    ## ... Similar to previous best
    ## Run 370 stress 0.09592158 
    ## Run 371 stress 0.08946329 
    ## ... Procrustes: rmse 0.03753148  max resid 0.1178014 
    ## Run 372 stress 0.1101445 
    ## Run 373 stress 0.1071321 
    ## Run 374 stress 0.1075422 
    ## Run 375 stress 0.0909961 
    ## Run 376 stress 0.09503417 
    ## Run 377 stress 0.1071321 
    ## Run 378 stress 0.09130089 
    ## Run 379 stress 0.08938962 
    ## ... Procrustes: rmse 0.03596136  max resid 0.1182964 
    ## Run 380 stress 0.1089903 
    ## Run 381 stress 0.08926073 
    ## ... Procrustes: rmse 0.0002024719  max resid 0.0006359162 
    ## ... Similar to previous best
    ## Run 382 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595531  max resid 0.1182967 
    ## Run 383 stress 0.1064187 
    ## Run 384 stress 0.09503413 
    ## Run 385 stress 0.09130085 
    ## Run 386 stress 0.08926101 
    ## ... Procrustes: rmse 0.0003900543  max resid 0.001273453 
    ## ... Similar to previous best
    ## Run 387 stress 0.105265 
    ## Run 388 stress 0.08938965 
    ## ... Procrustes: rmse 0.03597468  max resid 0.1183305 
    ## Run 389 stress 0.08938541 
    ## ... Procrustes: rmse 0.01185816  max resid 0.03965022 
    ## Run 390 stress 0.08926066 
    ## ... Procrustes: rmse 9.620532e-05  max resid 0.0003011429 
    ## ... Similar to previous best
    ## Run 391 stress 0.08926075 
    ## ... Procrustes: rmse 0.0002118091  max resid 0.0006778515 
    ## ... Similar to previous best
    ## Run 392 stress 0.109348 
    ## Run 393 stress 0.09021167 
    ## Run 394 stress 0.1074433 
    ## Run 395 stress 0.09088608 
    ## Run 396 stress 0.0901871 
    ## Run 397 stress 0.106538 
    ## Run 398 stress 0.0893854 
    ## ... Procrustes: rmse 0.01184747  max resid 0.03971091 
    ## Run 399 stress 0.08926076 
    ## ... Procrustes: rmse 0.0002256893  max resid 0.0007212139 
    ## ... Similar to previous best
    ## Run 400 stress 0.08946329 
    ## ... Procrustes: rmse 0.03750302  max resid 0.1177522 
    ## Run 401 stress 0.09099535 
    ## Run 402 stress 0.09039132 
    ## Run 403 stress 0.09180876 
    ## Run 404 stress 0.08926084 
    ## ... Procrustes: rmse 0.0002887006  max resid 0.0009055635 
    ## ... Similar to previous best
    ## Run 405 stress 0.09228484 
    ## Run 406 stress 0.1060429 
    ## Run 407 stress 0.1052647 
    ## Run 408 stress 0.08926094 
    ## ... Procrustes: rmse 0.0003507355  max resid 0.00114498 
    ## ... Similar to previous best
    ## Run 409 stress 0.08926086 
    ## ... Procrustes: rmse 0.0002707554  max resid 0.0008997917 
    ## ... Similar to previous best
    ## Run 410 stress 0.08946331 
    ## ... Procrustes: rmse 0.03749837  max resid 0.1177473 
    ## Run 411 stress 0.09503431 
    ## Run 412 stress 0.1071322 
    ## Run 413 stress 0.1074432 
    ## Run 414 stress 0.08938961 
    ## ... Procrustes: rmse 0.03594869  max resid 0.1182831 
    ## Run 415 stress 0.1071322 
    ## Run 416 stress 0.1056899 
    ## Run 417 stress 0.09130088 
    ## Run 418 stress 0.09018709 
    ## Run 419 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183649  max resid 0.03969973 
    ## Run 420 stress 0.1070075 
    ## Run 421 stress 0.08946335 
    ## ... Procrustes: rmse 0.03748772  max resid 0.1177325 
    ## Run 422 stress 0.08926082 
    ## ... Procrustes: rmse 0.000256067  max resid 0.0007989944 
    ## ... Similar to previous best
    ## Run 423 stress 0.08926074 
    ## ... Procrustes: rmse 0.0001967938  max resid 0.0006164098 
    ## ... Similar to previous best
    ## Run 424 stress 0.1074434 
    ## Run 425 stress 0.1074432 
    ## Run 426 stress 0.1087862 
    ## Run 427 stress 0.09228482 
    ## Run 428 stress 0.1056902 
    ## Run 429 stress 0.08946654 
    ## ... Procrustes: rmse 0.03320783  max resid 0.1163269 
    ## Run 430 stress 0.0903913 
    ## Run 431 stress 0.1056903 
    ## Run 432 stress 0.09109111 
    ## Run 433 stress 0.08938541 
    ## ... Procrustes: rmse 0.01183235  max resid 0.0396572 
    ## Run 434 stress 0.09130087 
    ## Run 435 stress 0.1092131 
    ## Run 436 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595717  max resid 0.1182985 
    ## Run 437 stress 0.09039129 
    ## Run 438 stress 0.08938971 
    ## ... Procrustes: rmse 0.03591934  max resid 0.1182303 
    ## Run 439 stress 0.08926064 
    ## ... Procrustes: rmse 8.645075e-05  max resid 0.0002007189 
    ## ... Similar to previous best
    ## Run 440 stress 0.09503415 
    ## Run 441 stress 0.0917828 
    ## Run 442 stress 0.090886 
    ## Run 443 stress 0.08926089 
    ## ... Procrustes: rmse 0.0003432301  max resid 0.001013697 
    ## ... Similar to previous best
    ## Run 444 stress 0.08938969 
    ## ... Procrustes: rmse 0.03592419  max resid 0.1182308 
    ## Run 445 stress 0.1087863 
    ## Run 446 stress 0.1087863 
    ## Run 447 stress 0.09130085 
    ## Run 448 stress 0.0904461 
    ## Run 449 stress 0.1056894 
    ## Run 450 stress 0.09088601 
    ## Run 451 stress 0.08926115 
    ## ... Procrustes: rmse 0.0004572491  max resid 0.001507529 
    ## ... Similar to previous best
    ## Run 452 stress 0.1110409 
    ## Run 453 stress 0.08946654 
    ## ... Procrustes: rmse 0.03319954  max resid 0.1163134 
    ## Run 454 stress 0.1074433 
    ## Run 455 stress 0.1074434 
    ## Run 456 stress 0.08926076 
    ## ... Procrustes: rmse 0.0002421897  max resid 0.0007265176 
    ## ... Similar to previous best
    ## Run 457 stress 0.1079014 
    ## Run 458 stress 0.1052654 
    ## Run 459 stress 0.09503428 
    ## Run 460 stress 0.08938556 
    ## ... Procrustes: rmse 0.0118344  max resid 0.03963456 
    ## Run 461 stress 0.09018709 
    ## Run 462 stress 0.08938556 
    ## ... Procrustes: rmse 0.01175338  max resid 0.03959034 
    ## Run 463 stress 0.09180889 
    ## Run 464 stress 0.1101452 
    ## Run 465 stress 0.1065367 
    ## Run 466 stress 0.1065374 
    ## Run 467 stress 0.1065368 
    ## Run 468 stress 0.09594309 
    ## Run 469 stress 0.09044609 
    ## Run 470 stress 0.08938538 
    ## ... Procrustes: rmse 0.01185201  max resid 0.03967467 
    ## Run 471 stress 0.08938965 
    ## ... Procrustes: rmse 0.03597116  max resid 0.1183255 
    ## Run 472 stress 0.0910911 
    ## Run 473 stress 0.09262365 
    ## Run 474 stress 0.09130085 
    ## Run 475 stress 0.1112206 
    ## Run 476 stress 0.1068388 
    ## Run 477 stress 0.1087862 
    ## Run 478 stress 0.1074433 
    ## Run 479 stress 0.1092095 
    ## Run 480 stress 0.08926079 
    ## ... Procrustes: rmse 0.0002548275  max resid 0.0008012324 
    ## ... Similar to previous best
    ## Run 481 stress 0.1092094 
    ## Run 482 stress 0.08946664 
    ## ... Procrustes: rmse 0.03318191  max resid 0.1162864 
    ## Run 483 stress 0.1104478 
    ## Run 484 stress 0.1079015 
    ## Run 485 stress 0.09503414 
    ## Run 486 stress 0.09503414 
    ## Run 487 stress 0.110145 
    ## Run 488 stress 0.1062596 
    ## Run 489 stress 0.1105341 
    ## Run 490 stress 0.1074431 
    ## Run 491 stress 0.1075421 
    ## Run 492 stress 0.1071321 
    ## Run 493 stress 0.1076303 
    ## Run 494 stress 0.09021165 
    ## Run 495 stress 0.1052648 
    ## Run 496 stress 0.1095625 
    ## Run 497 stress 0.08946655 
    ## ... Procrustes: rmse 0.03323927  max resid 0.1164013 
    ## Run 498 stress 0.09228493 
    ## Run 499 stress 0.1064187 
    ## Run 500 stress 0.08938965 
    ## ... Procrustes: rmse 0.0359353  max resid 0.1182635 
    ## *** Best solution repeated 63 times

``` r
round(SD_beta_NMDS$stress, digits = 2)
```

    ## [1] 0.09

``` r
SD_beta_rep_NMDS <- metaMDS(SD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.3323292 
    ## Run 1 stress 0.3438849 
    ## Run 2 stress 0.3313826 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1875562  max resid 0.3235173 
    ## Run 3 stress 0.3438453 
    ## Run 4 stress 0.3396612 
    ## Run 5 stress 0.3345052 
    ## Run 6 stress 0.3334622 
    ## Run 7 stress 0.3399265 
    ## Run 8 stress 0.3292168 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1935698  max resid 0.3172013 
    ## Run 9 stress 0.3438632 
    ## Run 10 stress 0.3456987 
    ## Run 11 stress 0.3322101 
    ## Run 12 stress 0.3363682 
    ## Run 13 stress 0.330415 
    ## Run 14 stress 0.3429478 
    ## Run 15 stress 0.3427563 
    ## Run 16 stress 0.3309067 
    ## Run 17 stress 0.3605288 
    ## Run 18 stress 0.336504 
    ## Run 19 stress 0.3363981 
    ## Run 20 stress 0.3384437 
    ## Run 21 stress 0.3389782 
    ## Run 22 stress 0.342806 
    ## Run 23 stress 0.3398736 
    ## Run 24 stress 0.3351316 
    ## Run 25 stress 0.3410116 
    ## Run 26 stress 0.3313125 
    ## Run 27 stress 0.3448259 
    ## Run 28 stress 0.3348797 
    ## Run 29 stress 0.3342216 
    ## Run 30 stress 0.3292655 
    ## ... Procrustes: rmse 0.1802376  max resid 0.3160985 
    ## Run 31 stress 0.3343054 
    ## Run 32 stress 0.3280802 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1736106  max resid 0.386634 
    ## Run 33 stress 0.3396287 
    ## Run 34 stress 0.3454045 
    ## Run 35 stress 0.3394914 
    ## Run 36 stress 0.3356234 
    ## Run 37 stress 0.3384904 
    ## Run 38 stress 0.3353835 
    ## Run 39 stress 0.3318157 
    ## Run 40 stress 0.3363238 
    ## Run 41 stress 0.3553604 
    ## Run 42 stress 0.344804 
    ## Run 43 stress 0.3419658 
    ## Run 44 stress 0.3300061 
    ## Run 45 stress 0.3413934 
    ## Run 46 stress 0.3370891 
    ## Run 47 stress 0.3326376 
    ## Run 48 stress 0.3328753 
    ## Run 49 stress 0.3409522 
    ## Run 50 stress 0.3294603 
    ## Run 51 stress 0.3322351 
    ## Run 52 stress 0.3383577 
    ## Run 53 stress 0.3318293 
    ## Run 54 stress 0.3361588 
    ## Run 55 stress 0.3386114 
    ## Run 56 stress 0.339543 
    ## Run 57 stress 0.3351969 
    ## Run 58 stress 0.3287107 
    ## Run 59 stress 0.3425121 
    ## Run 60 stress 0.3427424 
    ## Run 61 stress 0.3465535 
    ## Run 62 stress 0.3372928 
    ## Run 63 stress 0.3362757 
    ## Run 64 stress 0.3434557 
    ## Run 65 stress 0.3469166 
    ## Run 66 stress 0.3424888 
    ## Run 67 stress 0.339164 
    ## Run 68 stress 0.325534 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1789896  max resid 0.3506382 
    ## Run 69 stress 0.3288532 
    ## Run 70 stress 0.3386997 
    ## Run 71 stress 0.332183 
    ## Run 72 stress 0.3291617 
    ## Run 73 stress 0.3415592 
    ## Run 74 stress 0.3337055 
    ## Run 75 stress 0.3313418 
    ## Run 76 stress 0.3290986 
    ## Run 77 stress 0.3404822 
    ## Run 78 stress 0.3371796 
    ## Run 79 stress 0.3808259 
    ## Run 80 stress 0.336428 
    ## Run 81 stress 0.3330478 
    ## Run 82 stress 0.3366321 
    ## Run 83 stress 0.3434565 
    ## Run 84 stress 0.3346109 
    ## Run 85 stress 0.3313328 
    ## Run 86 stress 0.3316578 
    ## Run 87 stress 0.3359086 
    ## Run 88 stress 0.3411914 
    ## Run 89 stress 0.3466392 
    ## Run 90 stress 0.3347433 
    ## Run 91 stress 0.341855 
    ## Run 92 stress 0.3367039 
    ## Run 93 stress 0.3338584 
    ## Run 94 stress 0.3298436 
    ## Run 95 stress 0.3355686 
    ## Run 96 stress 0.3372788 
    ## Run 97 stress 0.3343666 
    ## Run 98 stress 0.3404573 
    ## Run 99 stress 0.3372977 
    ## Run 100 stress 0.3335597 
    ## Run 101 stress 0.343761 
    ## Run 102 stress 0.3301475 
    ## Run 103 stress 0.3319147 
    ## Run 104 stress 0.3519345 
    ## Run 105 stress 0.3426011 
    ## Run 106 stress 0.3440509 
    ## Run 107 stress 0.3464932 
    ## Run 108 stress 0.3339012 
    ## Run 109 stress 0.3481119 
    ## Run 110 stress 0.3335243 
    ## Run 111 stress 0.3311747 
    ## Run 112 stress 0.3341786 
    ## Run 113 stress 0.3448591 
    ## Run 114 stress 0.3418538 
    ## Run 115 stress 0.3304943 
    ## Run 116 stress 0.3337315 
    ## Run 117 stress 0.3381972 
    ## Run 118 stress 0.3307728 
    ## Run 119 stress 0.3374929 
    ## Run 120 stress 0.3460875 
    ## Run 121 stress 0.335774 
    ## Run 122 stress 0.3335681 
    ## Run 123 stress 0.3396239 
    ## Run 124 stress 0.3420946 
    ## Run 125 stress 0.3324662 
    ## Run 126 stress 0.3383415 
    ## Run 127 stress 0.3330182 
    ## Run 128 stress 0.3468863 
    ## Run 129 stress 0.3442853 
    ## Run 130 stress 0.331445 
    ## Run 131 stress 0.3442807 
    ## Run 132 stress 0.3318222 
    ## Run 133 stress 0.3376553 
    ## Run 134 stress 0.3359784 
    ## Run 135 stress 0.3363374 
    ## Run 136 stress 0.3438177 
    ## Run 137 stress 0.3319157 
    ## Run 138 stress 0.336262 
    ## Run 139 stress 0.3304532 
    ## Run 140 stress 0.3398185 
    ## Run 141 stress 0.3292565 
    ## Run 142 stress 0.3350046 
    ## Run 143 stress 0.343697 
    ## Run 144 stress 0.3348975 
    ## Run 145 stress 0.3468298 
    ## Run 146 stress 0.3331619 
    ## Run 147 stress 0.3451933 
    ## Run 148 stress 0.3347245 
    ## Run 149 stress 0.3365854 
    ## Run 150 stress 0.348793 
    ## Run 151 stress 0.3377054 
    ## Run 152 stress 0.3420979 
    ## Run 153 stress 0.336511 
    ## Run 154 stress 0.3297906 
    ## Run 155 stress 0.3348004 
    ## Run 156 stress 0.3284851 
    ## Run 157 stress 0.3301419 
    ## Run 158 stress 0.3476329 
    ## Run 159 stress 0.3337126 
    ## Run 160 stress 0.3386043 
    ## Run 161 stress 0.3437201 
    ## Run 162 stress 0.3409777 
    ## Run 163 stress 0.3471249 
    ## Run 164 stress 0.3407451 
    ## Run 165 stress 0.329707 
    ## Run 166 stress 0.3385911 
    ## Run 167 stress 0.3321302 
    ## Run 168 stress 0.3317374 
    ## Run 169 stress 0.3338228 
    ## Run 170 stress 0.33707 
    ## Run 171 stress 0.3487593 
    ## Run 172 stress 0.3404889 
    ## Run 173 stress 0.3467759 
    ## Run 174 stress 0.3301916 
    ## Run 175 stress 0.3323833 
    ## Run 176 stress 0.339324 
    ## Run 177 stress 0.3373948 
    ## Run 178 stress 0.3314459 
    ## Run 179 stress 0.327329 
    ## Run 180 stress 0.327391 
    ## Run 181 stress 0.3401802 
    ## Run 182 stress 0.3405971 
    ## Run 183 stress 0.3258649 
    ## ... Procrustes: rmse 0.181246  max resid 0.3402828 
    ## Run 184 stress 0.3292904 
    ## Run 185 stress 0.3399171 
    ## Run 186 stress 0.3445407 
    ## Run 187 stress 0.3364159 
    ## Run 188 stress 0.3323584 
    ## Run 189 stress 0.332056 
    ## Run 190 stress 0.348917 
    ## Run 191 stress 0.3283292 
    ## Run 192 stress 0.3433115 
    ## Run 193 stress 0.3461102 
    ## Run 194 stress 0.3474996 
    ## Run 195 stress 0.3512773 
    ## Run 196 stress 0.3295535 
    ## Run 197 stress 0.3270134 
    ## Run 198 stress 0.3440879 
    ## Run 199 stress 0.3403343 
    ## Run 200 stress 0.3331106 
    ## Run 201 stress 0.3386676 
    ## Run 202 stress 0.3429706 
    ## Run 203 stress 0.3439385 
    ## Run 204 stress 0.3280086 
    ## Run 205 stress 0.3420105 
    ## Run 206 stress 0.3336888 
    ## Run 207 stress 0.33001 
    ## Run 208 stress 0.3242537 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1847945  max resid 0.3362757 
    ## Run 209 stress 0.331819 
    ## Run 210 stress 0.328811 
    ## Run 211 stress 0.3337002 
    ## Run 212 stress 0.3324974 
    ## Run 213 stress 0.3450134 
    ## Run 214 stress 0.328687 
    ## Run 215 stress 0.3394514 
    ## Run 216 stress 0.3438811 
    ## Run 217 stress 0.3328294 
    ## Run 218 stress 0.3475207 
    ## Run 219 stress 0.3592185 
    ## Run 220 stress 0.3374848 
    ## Run 221 stress 0.3329102 
    ## Run 222 stress 0.3419074 
    ## Run 223 stress 0.3440986 
    ## Run 224 stress 0.3483398 
    ## Run 225 stress 0.3320563 
    ## Run 226 stress 0.3357645 
    ## Run 227 stress 0.3434306 
    ## Run 228 stress 0.3336558 
    ## Run 229 stress 0.3356499 
    ## Run 230 stress 0.3491683 
    ## Run 231 stress 0.3295081 
    ## Run 232 stress 0.3351802 
    ## Run 233 stress 0.3420812 
    ## Run 234 stress 0.3297691 
    ## Run 235 stress 0.3386121 
    ## Run 236 stress 0.3344759 
    ## Run 237 stress 0.3343056 
    ## Run 238 stress 0.3483676 
    ## Run 239 stress 0.334009 
    ## Run 240 stress 0.3355841 
    ## Run 241 stress 0.3410549 
    ## Run 242 stress 0.3438514 
    ## Run 243 stress 0.3438266 
    ## Run 244 stress 0.3371121 
    ## Run 245 stress 0.3318583 
    ## Run 246 stress 0.3406586 
    ## Run 247 stress 0.3465097 
    ## Run 248 stress 0.3477493 
    ## Run 249 stress 0.340242 
    ## Run 250 stress 0.3309121 
    ## Run 251 stress 0.3344815 
    ## Run 252 stress 0.3448323 
    ## Run 253 stress 0.338528 
    ## Run 254 stress 0.341805 
    ## Run 255 stress 0.3379833 
    ## Run 256 stress 0.3443346 
    ## Run 257 stress 0.3298366 
    ## Run 258 stress 0.3295429 
    ## Run 259 stress 0.3285142 
    ## Run 260 stress 0.3339231 
    ## Run 261 stress 0.3368483 
    ## Run 262 stress 0.3322527 
    ## Run 263 stress 0.3376101 
    ## Run 264 stress 0.3409741 
    ## Run 265 stress 0.337151 
    ## Run 266 stress 0.3309346 
    ## Run 267 stress 0.3450239 
    ## Run 268 stress 0.3304357 
    ## Run 269 stress 0.3345948 
    ## Run 270 stress 0.3307988 
    ## Run 271 stress 0.3277823 
    ## Run 272 stress 0.3299923 
    ## Run 273 stress 0.3335028 
    ## Run 274 stress 0.3358454 
    ## Run 275 stress 0.3451111 
    ## Run 276 stress 0.331443 
    ## Run 277 stress 0.3467951 
    ## Run 278 stress 0.3443746 
    ## Run 279 stress 0.3390654 
    ## Run 280 stress 0.3284132 
    ## Run 281 stress 0.3398202 
    ## Run 282 stress 0.3397652 
    ## Run 283 stress 0.3320031 
    ## Run 284 stress 0.3355362 
    ## Run 285 stress 0.33795 
    ## Run 286 stress 0.3336938 
    ## Run 287 stress 0.3399643 
    ## Run 288 stress 0.3386867 
    ## Run 289 stress 0.3334185 
    ## Run 290 stress 0.336862 
    ## Run 291 stress 0.3316441 
    ## Run 292 stress 0.3344898 
    ## Run 293 stress 0.3402623 
    ## Run 294 stress 0.3286306 
    ## Run 295 stress 0.3398111 
    ## Run 296 stress 0.3293315 
    ## Run 297 stress 0.3429316 
    ## Run 298 stress 0.329996 
    ## Run 299 stress 0.3353126 
    ## Run 300 stress 0.3280773 
    ## Run 301 stress 0.3330726 
    ## Run 302 stress 0.3330743 
    ## Run 303 stress 0.3408364 
    ## Run 304 stress 0.3441878 
    ## Run 305 stress 0.3491708 
    ## Run 306 stress 0.3350697 
    ## Run 307 stress 0.3513105 
    ## Run 308 stress 0.3476291 
    ## Run 309 stress 0.329246 
    ## Run 310 stress 0.3435741 
    ## Run 311 stress 0.3283813 
    ## Run 312 stress 0.3497369 
    ## Run 313 stress 0.3376306 
    ## Run 314 stress 0.3362834 
    ## Run 315 stress 0.3350173 
    ## Run 316 stress 0.3344166 
    ## Run 317 stress 0.3425978 
    ## Run 318 stress 0.3413768 
    ## Run 319 stress 0.3300395 
    ## Run 320 stress 0.3361225 
    ## Run 321 stress 0.3311536 
    ## Run 322 stress 0.3456079 
    ## Run 323 stress 0.3295504 
    ## Run 324 stress 0.3369879 
    ## Run 325 stress 0.3372771 
    ## Run 326 stress 0.3351684 
    ## Run 327 stress 0.326703 
    ## Run 328 stress 0.3363481 
    ## Run 329 stress 0.3421487 
    ## Run 330 stress 0.3369109 
    ## Run 331 stress 0.3382731 
    ## Run 332 stress 0.3326322 
    ## Run 333 stress 0.3388145 
    ## Run 334 stress 0.3425502 
    ## Run 335 stress 0.3314287 
    ## Run 336 stress 0.3340427 
    ## Run 337 stress 0.3369072 
    ## Run 338 stress 0.341641 
    ## Run 339 stress 0.3506775 
    ## Run 340 stress 0.344877 
    ## Run 341 stress 0.3317099 
    ## Run 342 stress 0.336769 
    ## Run 343 stress 0.325952 
    ## Run 344 stress 0.3336531 
    ## Run 345 stress 0.3394699 
    ## Run 346 stress 0.3401614 
    ## Run 347 stress 0.3317083 
    ## Run 348 stress 0.3335079 
    ## Run 349 stress 0.3473721 
    ## Run 350 stress 0.3303376 
    ## Run 351 stress 0.335463 
    ## Run 352 stress 0.3450322 
    ## Run 353 stress 0.3309541 
    ## Run 354 stress 0.3352746 
    ## Run 355 stress 0.3344638 
    ## Run 356 stress 0.3302149 
    ## Run 357 stress 0.3375589 
    ## Run 358 stress 0.3402802 
    ## Run 359 stress 0.3282772 
    ## Run 360 stress 0.330535 
    ## Run 361 stress 0.3362457 
    ## Run 362 stress 0.331544 
    ## Run 363 stress 0.3409137 
    ## Run 364 stress 0.3321062 
    ## Run 365 stress 0.3338079 
    ## Run 366 stress 0.3411693 
    ## Run 367 stress 0.3316805 
    ## Run 368 stress 0.3341696 
    ## Run 369 stress 0.3279373 
    ## Run 370 stress 0.3344835 
    ## Run 371 stress 0.337042 
    ## Run 372 stress 0.3463835 
    ## Run 373 stress 0.3355207 
    ## Run 374 stress 0.344082 
    ## Run 375 stress 0.3295473 
    ## Run 376 stress 0.33597 
    ## Run 377 stress 0.3330093 
    ## Run 378 stress 0.3437215 
    ## Run 379 stress 0.3389428 
    ## Run 380 stress 0.3392227 
    ## Run 381 stress 0.3348904 
    ## Run 382 stress 0.349012 
    ## Run 383 stress 0.3360242 
    ## Run 384 stress 0.3448889 
    ## Run 385 stress 0.3453637 
    ## Run 386 stress 0.3343948 
    ## Run 387 stress 0.3280195 
    ## Run 388 stress 0.3318716 
    ## Run 389 stress 0.3277585 
    ## Run 390 stress 0.3330286 
    ## Run 391 stress 0.3324469 
    ## Run 392 stress 0.3353504 
    ## Run 393 stress 0.3355447 
    ## Run 394 stress 0.3346846 
    ## Run 395 stress 0.329276 
    ## Run 396 stress 0.3307716 
    ## Run 397 stress 0.3340668 
    ## Run 398 stress 0.3359915 
    ## Run 399 stress 0.3407872 
    ## Run 400 stress 0.333277 
    ## Run 401 stress 0.343555 
    ## Run 402 stress 0.3430576 
    ## Run 403 stress 0.3340725 
    ## Run 404 stress 0.3392725 
    ## Run 405 stress 0.3439954 
    ## Run 406 stress 0.3354198 
    ## Run 407 stress 0.3351865 
    ## Run 408 stress 0.3351945 
    ## Run 409 stress 0.3319463 
    ## Run 410 stress 0.3466648 
    ## Run 411 stress 0.3483544 
    ## Run 412 stress 0.3398834 
    ## Run 413 stress 0.3486056 
    ## Run 414 stress 0.3381068 
    ## Run 415 stress 0.3356493 
    ## Run 416 stress 0.3339222 
    ## Run 417 stress 0.3430428 
    ## Run 418 stress 0.3325122 
    ## Run 419 stress 0.3392161 
    ## Run 420 stress 0.342003 
    ## Run 421 stress 0.3368581 
    ## Run 422 stress 0.3345518 
    ## Run 423 stress 0.3285676 
    ## Run 424 stress 0.3320915 
    ## Run 425 stress 0.3359493 
    ## Run 426 stress 0.3303717 
    ## Run 427 stress 0.3289929 
    ## Run 428 stress 0.3453682 
    ## Run 429 stress 0.3368912 
    ## Run 430 stress 0.3371861 
    ## Run 431 stress 0.3458225 
    ## Run 432 stress 0.3330074 
    ## Run 433 stress 0.3346727 
    ## Run 434 stress 0.3426274 
    ## Run 435 stress 0.3627181 
    ## Run 436 stress 0.3481037 
    ## Run 437 stress 0.3419785 
    ## Run 438 stress 0.3393398 
    ## Run 439 stress 0.3397185 
    ## Run 440 stress 0.3491664 
    ## Run 441 stress 0.340655 
    ## Run 442 stress 0.33132 
    ## Run 443 stress 0.3309528 
    ## Run 444 stress 0.3461735 
    ## Run 445 stress 0.3268308 
    ## Run 446 stress 0.3291617 
    ## Run 447 stress 0.3375991 
    ## Run 448 stress 0.3359141 
    ## Run 449 stress 0.3453079 
    ## Run 450 stress 0.3386905 
    ## Run 451 stress 0.325853 
    ## Run 452 stress 0.33332 
    ## Run 453 stress 0.3272548 
    ## Run 454 stress 0.3296631 
    ## Run 455 stress 0.3313429 
    ## Run 456 stress 0.3309157 
    ## Run 457 stress 0.3343762 
    ## Run 458 stress 0.3366072 
    ## Run 459 stress 0.3380961 
    ## Run 460 stress 0.3296041 
    ## Run 461 stress 0.3394707 
    ## Run 462 stress 0.3341267 
    ## Run 463 stress 0.3360659 
    ## Run 464 stress 0.3355347 
    ## Run 465 stress 0.3480926 
    ## Run 466 stress 0.3309201 
    ## Run 467 stress 0.3473959 
    ## Run 468 stress 0.3456119 
    ## Run 469 stress 0.3428278 
    ## Run 470 stress 0.3363642 
    ## Run 471 stress 0.3453023 
    ## Run 472 stress 0.3322731 
    ## Run 473 stress 0.3428237 
    ## Run 474 stress 0.340132 
    ## Run 475 stress 0.3314935 
    ## Run 476 stress 0.3404583 
    ## Run 477 stress 0.3458538 
    ## Run 478 stress 0.347377 
    ## Run 479 stress 0.3354524 
    ## Run 480 stress 0.3397587 
    ## Run 481 stress 0.3325377 
    ## Run 482 stress 0.3408348 
    ## Run 483 stress 0.3451274 
    ## Run 484 stress 0.3377276 
    ## Run 485 stress 0.3425459 
    ## Run 486 stress 0.3332066 
    ## Run 487 stress 0.337178 
    ## Run 488 stress 0.3431712 
    ## Run 489 stress 0.3461457 
    ## Run 490 stress 0.3426486 
    ## Run 491 stress 0.3365851 
    ## Run 492 stress 0.3326357 
    ## Run 493 stress 0.3424392 
    ## Run 494 stress 0.3415257 
    ## Run 495 stress 0.3345646 
    ## Run 496 stress 0.3346753 
    ## Run 497 stress 0.3419071 
    ## Run 498 stress 0.3321605 
    ## Run 499 stress 0.3335407 
    ## Run 500 stress 0.3289537 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##    500: stress ratio > sratmax

``` r
SD_beta_ric_NMDS <- metaMDS(SD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.01879552 
    ## Run 1 stress 0.01884066 
    ## ... Procrustes: rmse 0.003386456  max resid 0.006558902 
    ## ... Similar to previous best
    ## Run 2 stress 0.01879584 
    ## ... Procrustes: rmse 0.001282623  max resid 0.002641529 
    ## ... Similar to previous best
    ## Run 3 stress 0.01880071 
    ## ... Procrustes: rmse 0.002340931  max resid 0.004828637 
    ## ... Similar to previous best
    ## Run 4 stress 0.02492194 
    ## Run 5 stress 0.02520901 
    ## Run 6 stress 0.02486139 
    ## Run 7 stress 0.0189414 
    ## ... Procrustes: rmse 0.009145635  max resid 0.01860392 
    ## Run 8 stress 0.01892328 
    ## ... Procrustes: rmse 0.009345857  max resid 0.01917718 
    ## Run 9 stress 0.02509444 
    ## Run 10 stress 0.0250945 
    ## Run 11 stress 0.02509442 
    ## Run 12 stress 0.02508855 
    ## Run 13 stress 0.02509469 
    ## Run 14 stress 0.01879573 
    ## ... Procrustes: rmse 0.0001031232  max resid 0.0002116188 
    ## ... Similar to previous best
    ## Run 15 stress 0.02492163 
    ## Run 16 stress 0.01888279 
    ## ... Procrustes: rmse 0.007090969  max resid 0.01433773 
    ## Run 17 stress 0.02486129 
    ## Run 18 stress 0.02486125 
    ## Run 19 stress 0.01879592 
    ## ... Procrustes: rmse 0.001302413  max resid 0.00268246 
    ## ... Similar to previous best
    ## Run 20 stress 0.02509446 
    ## Run 21 stress 0.01880161 
    ## ... Procrustes: rmse 0.00256221  max resid 0.005276269 
    ## ... Similar to previous best
    ## Run 22 stress 0.02520872 
    ## Run 23 stress 0.01879576 
    ## ... Procrustes: rmse 0.001245923  max resid 0.002565979 
    ## ... Similar to previous best
    ## Run 24 stress 0.01931573 
    ## Run 25 stress 0.02520893 
    ## Run 26 stress 0.01884347 
    ## ... Procrustes: rmse 0.004528189  max resid 0.008987882 
    ## ... Similar to previous best
    ## Run 27 stress 0.01879575 
    ## ... Procrustes: rmse 0.001245474  max resid 0.002565027 
    ## ... Similar to previous best
    ## Run 28 stress 0.01879609 
    ## ... Procrustes: rmse 0.001369666  max resid 0.002820436 
    ## ... Similar to previous best
    ## Run 29 stress 0.0248613 
    ## Run 30 stress 0.01879555 
    ## ... Procrustes: rmse 0.001146595  max resid 0.002361525 
    ## ... Similar to previous best
    ## Run 31 stress 0.0188753 
    ## ... Procrustes: rmse 0.006502295  max resid 0.0131099 
    ## Run 32 stress 0.01879758 
    ## ... Procrustes: rmse 0.001802097  max resid 0.003710752 
    ## ... Similar to previous best
    ## Run 33 stress 0.0249218 
    ## Run 34 stress 0.01883905 
    ## ... Procrustes: rmse 0.004085166  max resid 0.008051717 
    ## ... Similar to previous best
    ## Run 35 stress 0.01879584 
    ## ... Procrustes: rmse 0.001281423  max resid 0.002639038 
    ## ... Similar to previous best
    ## Run 36 stress 0.01879573 
    ## ... Procrustes: rmse 0.001225312  max resid 0.002523635 
    ## ... Similar to previous best
    ## Run 37 stress 0.01882447 
    ## ... Procrustes: rmse 0.004797154  max resid 0.009867612 
    ## ... Similar to previous best
    ## Run 38 stress 0.02509432 
    ## Run 39 stress 0.01884774 
    ## ... Procrustes: rmse 0.004899515  max resid 0.00976925 
    ## ... Similar to previous best
    ## Run 40 stress 0.02509442 
    ## Run 41 stress 0.01883116 
    ## ... Procrustes: rmse 0.0031373  max resid 0.00604766 
    ## ... Similar to previous best
    ## Run 42 stress 0.01883006 
    ## ... Procrustes: rmse 0.002966119  max resid 0.005694754 
    ## ... Similar to previous best
    ## Run 43 stress 0.01879569 
    ## ... Procrustes: rmse 0.001214233  max resid 0.002500788 
    ## ... Similar to previous best
    ## Run 44 stress 0.02509455 
    ## Run 45 stress 0.01879585 
    ## ... Procrustes: rmse 0.001285329  max resid 0.002647037 
    ## ... Similar to previous best
    ## Run 46 stress 0.01879585 
    ## ... Procrustes: rmse 0.001274816  max resid 0.002625266 
    ## ... Similar to previous best
    ## Run 47 stress 0.02486121 
    ## Run 48 stress 0.01879974 
    ## ... Procrustes: rmse 0.002251457  max resid 0.004636479 
    ## ... Similar to previous best
    ## Run 49 stress 0.02486119 
    ## Run 50 stress 0.01879874 
    ## ... Procrustes: rmse 0.002059818  max resid 0.004244694 
    ## ... Similar to previous best
    ## Run 51 stress 0.0188033 
    ## ... Procrustes: rmse 0.00280205  max resid 0.005770067 
    ## ... Similar to previous best
    ## Run 52 stress 0.02509459 
    ## Run 53 stress 0.02231871 
    ## Run 54 stress 0.02492193 
    ## Run 55 stress 0.01879562 
    ## ... Procrustes: rmse 4.096588e-05  max resid 8.259932e-05 
    ## ... Similar to previous best
    ## Run 56 stress 0.01894011 
    ## ... Procrustes: rmse 0.009583226  max resid 0.01950297 
    ## Run 57 stress 0.024922 
    ## Run 58 stress 0.0192859 
    ## ... Procrustes: rmse 0.01814475  max resid 0.03733788 
    ## Run 59 stress 0.02520921 
    ## Run 60 stress 0.02509462 
    ## Run 61 stress 0.01879605 
    ## ... Procrustes: rmse 0.001347856  max resid 0.002775473 
    ## ... Similar to previous best
    ## Run 62 stress 0.01880337 
    ## ... Procrustes: rmse 0.002759383  max resid 0.005680869 
    ## ... Similar to previous best
    ## Run 63 stress 0.01879577 
    ## ... Procrustes: rmse 0.001249274  max resid 0.002572891 
    ## ... Similar to previous best
    ## Run 64 stress 0.01880897 
    ## ... Procrustes: rmse 0.003471992  max resid 0.007148094 
    ## ... Similar to previous best
    ## Run 65 stress 0.02486101 
    ## Run 66 stress 0.0250943 
    ## Run 67 stress 0.01879563 
    ## ... Procrustes: rmse 0.001186301  max resid 0.002443256 
    ## ... Similar to previous best
    ## Run 68 stress 0.01881807 
    ## ... Procrustes: rmse 0.004307313  max resid 0.008862705 
    ## ... Similar to previous best
    ## Run 69 stress 0.01879567 
    ## ... Procrustes: rmse 0.001209086  max resid 0.002490171 
    ## ... Similar to previous best
    ## Run 70 stress 0.2039098 
    ## Run 71 stress 0.01888634 
    ## ... Procrustes: rmse 0.007979787  max resid 0.01639181 
    ## Run 72 stress 0.0187957 
    ## ... Procrustes: rmse 8.761384e-05  max resid 0.0001794773 
    ## ... Similar to previous best
    ## Run 73 stress 0.02492162 
    ## Run 74 stress 0.01879562 
    ## ... Procrustes: rmse 0.00118147  max resid 0.002433293 
    ## ... Similar to previous best
    ## Run 75 stress 0.01879587 
    ## ... Procrustes: rmse 0.000160331  max resid 0.0003290979 
    ## ... Similar to previous best
    ## Run 76 stress 0.01909035 
    ## ... Procrustes: rmse 0.01409279  max resid 0.02886821 
    ## Run 77 stress 0.02509433 
    ## Run 78 stress 0.01885034 
    ## ... Procrustes: rmse 0.005863203  max resid 0.0120555 
    ## Run 79 stress 0.02509479 
    ## Run 80 stress 0.01879576 
    ## ... Procrustes: rmse 0.001244808  max resid 0.002563706 
    ## ... Similar to previous best
    ## Run 81 stress 0.01879728 
    ## ... Procrustes: rmse 0.001740824  max resid 0.003584455 
    ## ... Similar to previous best
    ## Run 82 stress 0.02486129 
    ## Run 83 stress 0.02492185 
    ## Run 84 stress 0.02494809 
    ## Run 85 stress 0.2237308 
    ## Run 86 stress 0.02509439 
    ## Run 87 stress 0.02509468 
    ## Run 88 stress 0.01879594 
    ## ... Procrustes: rmse 0.0001763142  max resid 0.0003617768 
    ## ... Similar to previous best
    ## Run 89 stress 0.02509437 
    ## Run 90 stress 0.02486142 
    ## Run 91 stress 0.02509456 
    ## Run 92 stress 0.02486118 
    ## Run 93 stress 0.01879549 
    ## ... New best solution
    ## ... Procrustes: rmse 2.597712e-05  max resid 5.264656e-05 
    ## ... Similar to previous best
    ## Run 94 stress 0.02492179 
    ## Run 95 stress 0.02509451 
    ## Run 96 stress 0.01879563 
    ## ... Procrustes: rmse 0.001161038  max resid 0.002391555 
    ## ... Similar to previous best
    ## Run 97 stress 0.01884157 
    ## ... Procrustes: rmse 0.004331918  max resid 0.008572963 
    ## ... Similar to previous best
    ## Run 98 stress 0.02509456 
    ## Run 99 stress 0.01881418 
    ## ... Procrustes: rmse 0.003735448  max resid 0.007691621 
    ## ... Similar to previous best
    ## Run 100 stress 0.01879587 
    ## ... Procrustes: rmse 0.001269307  max resid 0.002614397 
    ## ... Similar to previous best
    ## Run 101 stress 0.01880811 
    ## ... Procrustes: rmse 0.002874672  max resid 0.005920682 
    ## ... Similar to previous best
    ## Run 102 stress 0.02486115 
    ## Run 103 stress 0.02492178 
    ## Run 104 stress 0.02486131 
    ## Run 105 stress 0.02509463 
    ## Run 106 stress 0.02509436 
    ## Run 107 stress 0.01881243 
    ## ... Procrustes: rmse 0.00377564  max resid 0.007771309 
    ## ... Similar to previous best
    ## Run 108 stress 0.02492184 
    ## Run 109 stress 0.01939749 
    ## Run 110 stress 0.01911528 
    ## ... Procrustes: rmse 0.01360131  max resid 0.02784848 
    ## Run 111 stress 0.02492218 
    ## Run 112 stress 0.01879579 
    ## ... Procrustes: rmse 0.0001520801  max resid 0.0003117539 
    ## ... Similar to previous best
    ## Run 113 stress 0.01879984 
    ## ... Procrustes: rmse 0.002248242  max resid 0.004629887 
    ## ... Similar to previous best
    ## Run 114 stress 0.02509454 
    ## Run 115 stress 0.01879589 
    ## ... Procrustes: rmse 0.001277828  max resid 0.002631911 
    ## ... Similar to previous best
    ## Run 116 stress 0.02509469 
    ## Run 117 stress 0.02509473 
    ## Run 118 stress 0.01888201 
    ## ... Procrustes: rmse 0.007451497  max resid 0.01529841 
    ## Run 119 stress 0.01911532 
    ## ... Procrustes: rmse 0.0146751  max resid 0.03006424 
    ## Run 120 stress 0.02509453 
    ## Run 121 stress 0.01879594 
    ## ... Procrustes: rmse 0.00129702  max resid 0.002671431 
    ## ... Similar to previous best
    ## Run 122 stress 0.01883605 
    ## ... Procrustes: rmse 0.003681879  max resid 0.007196239 
    ## ... Similar to previous best
    ## Run 123 stress 0.01879589 
    ## ... Procrustes: rmse 0.0001931565  max resid 0.0003957958 
    ## ... Similar to previous best
    ## Run 124 stress 0.01879573 
    ## ... Procrustes: rmse 0.0001289386  max resid 0.0002641534 
    ## ... Similar to previous best
    ## Run 125 stress 0.02492189 
    ## Run 126 stress 0.01891216 
    ## ... Procrustes: rmse 0.009024298  max resid 0.01852489 
    ## Run 127 stress 0.02509458 
    ## Run 128 stress 0.01968734 
    ## Run 129 stress 0.01879587 
    ## ... Procrustes: rmse 0.001200208  max resid 0.002471754 
    ## ... Similar to previous best
    ## Run 130 stress 0.02486119 
    ## Run 131 stress 0.02492187 
    ## Run 132 stress 0.02509473 
    ## Run 133 stress 0.02492178 
    ## Run 134 stress 0.02486134 
    ## Run 135 stress 0.02086657 
    ## Run 136 stress 0.01879577 
    ## ... Procrustes: rmse 0.001214895  max resid 0.002502454 
    ## ... Similar to previous best
    ## Run 137 stress 0.01879536 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001014238  max resid 0.002089378 
    ## ... Similar to previous best
    ## Run 138 stress 0.01879595 
    ## ... Procrustes: rmse 0.0002872627  max resid 0.0005918252 
    ## ... Similar to previous best
    ## Run 139 stress 0.01879556 
    ## ... Procrustes: rmse 0.0001118894  max resid 0.0002305314 
    ## ... Similar to previous best
    ## Run 140 stress 0.0188494 
    ## ... Procrustes: rmse 0.004018686  max resid 0.007891384 
    ## ... Similar to previous best
    ## Run 141 stress 0.02492169 
    ## Run 142 stress 0.01879581 
    ## ... Procrustes: rmse 0.0002289056  max resid 0.0004716184 
    ## ... Similar to previous best
    ## Run 143 stress 0.0188407 
    ## ... Procrustes: rmse 0.00480021  max resid 0.009874164 
    ## ... Similar to previous best
    ## Run 144 stress 0.0250944 
    ## Run 145 stress 0.02494716 
    ## Run 146 stress 0.02509455 
    ## Run 147 stress 0.02509448 
    ## Run 148 stress 0.01890221 
    ## ... Procrustes: rmse 0.007156181  max resid 0.01471466 
    ## Run 149 stress 0.02509481 
    ## Run 150 stress 0.01879589 
    ## ... Procrustes: rmse 0.001198894  max resid 0.002470538 
    ## ... Similar to previous best
    ## Run 151 stress 0.01883914 
    ## ... Procrustes: rmse 0.003130669  max resid 0.005982965 
    ## ... Similar to previous best
    ## Run 152 stress 0.02509461 
    ## Run 153 stress 0.02492171 
    ## Run 154 stress 0.02486102 
    ## Run 155 stress 0.01882779 
    ## ... Procrustes: rmse 0.001762197  max resid 0.004030995 
    ## ... Similar to previous best
    ## Run 156 stress 0.01879583 
    ## ... Procrustes: rmse 0.0002182485  max resid 0.0004502862 
    ## ... Similar to previous best
    ## Run 157 stress 0.01879543 
    ## ... Procrustes: rmse 0.0009878446  max resid 0.00203697 
    ## ... Similar to previous best
    ## Run 158 stress 0.01882304 
    ## ... Procrustes: rmse 0.003651258  max resid 0.007515223 
    ## ... Similar to previous best
    ## Run 159 stress 0.02509451 
    ## Run 160 stress 0.01879688 
    ## ... Procrustes: rmse 0.0005937059  max resid 0.001223122 
    ## ... Similar to previous best
    ## Run 161 stress 0.01886441 
    ## ... Procrustes: rmse 0.005057708  max resid 0.01008798 
    ## Run 162 stress 0.01899244 
    ## ... Procrustes: rmse 0.009450374  max resid 0.01936927 
    ## Run 163 stress 0.02492188 
    ## Run 164 stress 0.018796 
    ## ... Procrustes: rmse 0.001240639  max resid 0.002556512 
    ## ... Similar to previous best
    ## Run 165 stress 0.01884532 
    ## ... Procrustes: rmse 0.003697271  max resid 0.007204691 
    ## ... Similar to previous best
    ## Run 166 stress 0.01879578 
    ## ... Procrustes: rmse 0.0002162654  max resid 0.000445574 
    ## ... Similar to previous best
    ## Run 167 stress 0.0250948 
    ## Run 168 stress 0.01879615 
    ## ... Procrustes: rmse 0.000341681  max resid 0.0007037108 
    ## ... Similar to previous best
    ## Run 169 stress 0.02492198 
    ## Run 170 stress 0.02509449 
    ## Run 171 stress 0.02492477 
    ## Run 172 stress 0.02486132 
    ## Run 173 stress 0.01886425 
    ## ... Procrustes: rmse 0.006030308  max resid 0.01239962 
    ## Run 174 stress 0.01879556 
    ## ... Procrustes: rmse 0.0001039147  max resid 0.0002138668 
    ## ... Similar to previous best
    ## Run 175 stress 0.0187958 
    ## ... Procrustes: rmse 0.001171499  max resid 0.002414331 
    ## ... Similar to previous best
    ## Run 176 stress 0.01879564 
    ## ... Procrustes: rmse 0.0001541779  max resid 0.0003176656 
    ## ... Similar to previous best
    ## Run 177 stress 0.01879578 
    ## ... Procrustes: rmse 0.0002184219  max resid 0.0004500178 
    ## ... Similar to previous best
    ## Run 178 stress 0.02509462 
    ## Run 179 stress 0.02492867 
    ## Run 180 stress 0.01886694 
    ## ... Procrustes: rmse 0.00616136  max resid 0.01266714 
    ## Run 181 stress 0.01879572 
    ## ... Procrustes: rmse 0.0001870822  max resid 0.000385402 
    ## ... Similar to previous best
    ## Run 182 stress 0.02492167 
    ## Run 183 stress 0.0188288 
    ## ... Procrustes: rmse 0.00183309  max resid 0.004471591 
    ## ... Similar to previous best
    ## Run 184 stress 0.3513408 
    ## Run 185 stress 0.02509461 
    ## Run 186 stress 0.01879597 
    ## ... Procrustes: rmse 0.0002828473  max resid 0.0005825322 
    ## ... Similar to previous best
    ## Run 187 stress 0.02492178 
    ## Run 188 stress 0.02509469 
    ## Run 189 stress 0.01882102 
    ## ... Procrustes: rmse 0.003369221  max resid 0.006932919 
    ## ... Similar to previous best
    ## Run 190 stress 0.0248614 
    ## Run 191 stress 0.02486108 
    ## Run 192 stress 0.01879576 
    ## ... Procrustes: rmse 0.0002043881  max resid 0.0004210426 
    ## ... Similar to previous best
    ## Run 193 stress 0.01884804 
    ## ... Procrustes: rmse 0.0039151  max resid 0.007673097 
    ## ... Similar to previous best
    ## Run 194 stress 0.01883578 
    ## ... Procrustes: rmse 0.00277144  max resid 0.00543191 
    ## ... Similar to previous best
    ## Run 195 stress 0.0187961 
    ## ... Procrustes: rmse 0.0003343805  max resid 0.0006887907 
    ## ... Similar to previous best
    ## Run 196 stress 0.01879591 
    ## ... Procrustes: rmse 0.0002716799  max resid 0.0005597282 
    ## ... Similar to previous best
    ## Run 197 stress 0.02509452 
    ## Run 198 stress 0.02509476 
    ## Run 199 stress 0.01879598 
    ## ... Procrustes: rmse 0.0002829992  max resid 0.0005828511 
    ## ... Similar to previous best
    ## Run 200 stress 0.0187956 
    ## ... Procrustes: rmse 0.001079349  max resid 0.002224984 
    ## ... Similar to previous best
    ## Run 201 stress 0.01885766 
    ## ... Procrustes: rmse 0.004625372  max resid 0.009178353 
    ## ... Similar to previous best
    ## Run 202 stress 0.02509446 
    ## Run 203 stress 0.02509437 
    ## Run 204 stress 0.02492203 
    ## Run 205 stress 0.01879592 
    ## ... Procrustes: rmse 0.00120918  max resid 0.002491813 
    ## ... Similar to previous best
    ## Run 206 stress 0.02509479 
    ## Run 207 stress 0.01879557 
    ## ... Procrustes: rmse 0.0001191959  max resid 0.0002455775 
    ## ... Similar to previous best
    ## Run 208 stress 0.0250945 
    ## Run 209 stress 0.01879572 
    ## ... Procrustes: rmse 0.0001919816  max resid 0.0003955316 
    ## ... Similar to previous best
    ## Run 210 stress 0.01879582 
    ## ... Procrustes: rmse 0.0002345565  max resid 0.0004832531 
    ## ... Similar to previous best
    ## Run 211 stress 0.01879554 
    ## ... Procrustes: rmse 0.0001023915  max resid 0.0002108959 
    ## ... Similar to previous best
    ## Run 212 stress 0.02509446 
    ## Run 213 stress 0.02509481 
    ## Run 214 stress 0.02509467 
    ## Run 215 stress 0.01879556 
    ## ... Procrustes: rmse 0.0001154012  max resid 0.0002377671 
    ## ... Similar to previous best
    ## Run 216 stress 0.02509459 
    ## Run 217 stress 0.01880163 
    ## ... Procrustes: rmse 0.001525272  max resid 0.003143567 
    ## ... Similar to previous best
    ## Run 218 stress 0.01883033 
    ## ... Procrustes: rmse 0.002057057  max resid 0.004820795 
    ## ... Similar to previous best
    ## Run 219 stress 0.02509442 
    ## Run 220 stress 0.01890114 
    ## ... Procrustes: rmse 0.007600283  max resid 0.01561209 
    ## Run 221 stress 0.02486121 
    ## Run 222 stress 0.02509455 
    ## Run 223 stress 0.01883028 
    ## ... Procrustes: rmse 0.002041652  max resid 0.004794968 
    ## ... Similar to previous best
    ## Run 224 stress 0.01879574 
    ## ... Procrustes: rmse 0.0001981952  max resid 0.0004082915 
    ## ... Similar to previous best
    ## Run 225 stress 0.02492188 
    ## Run 226 stress 0.01879588 
    ## ... Procrustes: rmse 0.0002582801  max resid 0.0005321277 
    ## ... Similar to previous best
    ## Run 227 stress 0.02486134 
    ## Run 228 stress 0.02520891 
    ## Run 229 stress 0.01879559 
    ## ... Procrustes: rmse 0.000129024  max resid 0.0002658294 
    ## ... Similar to previous best
    ## Run 230 stress 0.0249218 
    ## Run 231 stress 0.01889181 
    ## ... Procrustes: rmse 0.007224693  max resid 0.0148446 
    ## Run 232 stress 0.0187959 
    ## ... Procrustes: rmse 0.000265849  max resid 0.0005477108 
    ## ... Similar to previous best
    ## Run 233 stress 0.01881594 
    ## ... Procrustes: rmse 0.002738268  max resid 0.005633362 
    ## ... Similar to previous best
    ## Run 234 stress 0.0249216 
    ## Run 235 stress 0.0250946 
    ## Run 236 stress 0.02492197 
    ## Run 237 stress 0.01879988 
    ## ... Procrustes: rmse 0.00124025  max resid 0.002556105 
    ## ... Similar to previous best
    ## Run 238 stress 0.01879572 
    ## ... Procrustes: rmse 0.0001909514  max resid 0.0003934096 
    ## ... Similar to previous best
    ## Run 239 stress 0.02492172 
    ## Run 240 stress 0.01879584 
    ## ... Procrustes: rmse 0.0002328468  max resid 0.0004796175 
    ## ... Similar to previous best
    ## Run 241 stress 0.01897443 
    ## ... Procrustes: rmse 0.01004171  max resid 0.02060487 
    ## Run 242 stress 0.01879567 
    ## ... Procrustes: rmse 0.000168936  max resid 0.0003480748 
    ## ... Similar to previous best
    ## Run 243 stress 0.3726999 
    ## Run 244 stress 0.02486128 
    ## Run 245 stress 0.02486098 
    ## Run 246 stress 0.02509449 
    ## Run 247 stress 0.01884233 
    ## ... Procrustes: rmse 0.003435751  max resid 0.006642674 
    ## ... Similar to previous best
    ## Run 248 stress 0.01879591 
    ## ... Procrustes: rmse 0.0002713032  max resid 0.0005589392 
    ## ... Similar to previous best
    ## Run 249 stress 0.01883591 
    ## ... Procrustes: rmse 0.004516846  max resid 0.009293733 
    ## ... Similar to previous best
    ## Run 250 stress 0.02492205 
    ## Run 251 stress 0.02087431 
    ## Run 252 stress 0.02509447 
    ## Run 253 stress 0.02520908 
    ## Run 254 stress 0.02492354 
    ## Run 255 stress 0.01879612 
    ## ... Procrustes: rmse 0.0003395177  max resid 0.0006992367 
    ## ... Similar to previous best
    ## Run 256 stress 0.02509445 
    ## Run 257 stress 0.01881306 
    ## ... Procrustes: rmse 0.002798213  max resid 0.005762505 
    ## ... Similar to previous best
    ## Run 258 stress 0.01888372 
    ## ... Procrustes: rmse 0.005939692  max resid 0.01193829 
    ## Run 259 stress 0.02509467 
    ## Run 260 stress 0.02509479 
    ## Run 261 stress 0.01879592 
    ## ... Procrustes: rmse 0.0002643518  max resid 0.0005444374 
    ## ... Similar to previous best
    ## Run 262 stress 0.01879596 
    ## ... Procrustes: rmse 0.0002852698  max resid 0.0005876516 
    ## ... Similar to previous best
    ## Run 263 stress 0.02486106 
    ## Run 264 stress 0.01880101 
    ## ... Procrustes: rmse 0.001430483  max resid 0.002948091 
    ## ... Similar to previous best
    ## Run 265 stress 0.01879556 
    ## ... Procrustes: rmse 0.001056053  max resid 0.002177189 
    ## ... Similar to previous best
    ## Run 266 stress 0.02509466 
    ## Run 267 stress 0.018833 
    ## ... Procrustes: rmse 0.004308971  max resid 0.008864657 
    ## ... Similar to previous best
    ## Run 268 stress 0.01969417 
    ## Run 269 stress 0.01879598 
    ## ... Procrustes: rmse 0.001224298  max resid 0.002522697 
    ## ... Similar to previous best
    ## Run 270 stress 0.02492173 
    ## Run 271 stress 0.01879599 
    ## ... Procrustes: rmse 0.0003009137  max resid 0.0006198802 
    ## ... Similar to previous best
    ## Run 272 stress 0.02509442 
    ## Run 273 stress 0.02509452 
    ## Run 274 stress 0.01879594 
    ## ... Procrustes: rmse 0.001219596  max resid 0.002513214 
    ## ... Similar to previous best
    ## Run 275 stress 0.01879566 
    ## ... Procrustes: rmse 0.001108134  max resid 0.002284155 
    ## ... Similar to previous best
    ## Run 276 stress 0.0248614 
    ## Run 277 stress 0.02509456 
    ## Run 278 stress 0.018796 
    ## ... Procrustes: rmse 0.0002926696  max resid 0.0006027728 
    ## ... Similar to previous best
    ## Run 279 stress 0.01879612 
    ## ... Procrustes: rmse 0.0003408575  max resid 0.0007020155 
    ## ... Similar to previous best
    ## Run 280 stress 0.02486105 
    ## Run 281 stress 0.01899845 
    ## ... Procrustes: rmse 0.010763  max resid 0.02207206 
    ## Run 282 stress 0.0195276 
    ## Run 283 stress 0.02509442 
    ## Run 284 stress 0.01879574 
    ## ... Procrustes: rmse 0.0001985682  max resid 0.0004091169 
    ## ... Similar to previous best
    ## Run 285 stress 0.01883792 
    ## ... Procrustes: rmse 0.003006275  max resid 0.005711942 
    ## ... Similar to previous best
    ## Run 286 stress 0.01879575 
    ## ... Procrustes: rmse 0.0001899887  max resid 0.0003911995 
    ## ... Similar to previous best
    ## Run 287 stress 0.02492807 
    ## Run 288 stress 0.01879565 
    ## ... Procrustes: rmse 0.0001586503  max resid 0.0003268955 
    ## ... Similar to previous best
    ## Run 289 stress 0.02509448 
    ## Run 290 stress 0.0188499 
    ## ... Procrustes: rmse 0.004062049  max resid 0.007985094 
    ## ... Similar to previous best
    ## Run 291 stress 0.0187959 
    ## ... Procrustes: rmse 0.001203071  max resid 0.002479261 
    ## ... Similar to previous best
    ## Run 292 stress 0.01879583 
    ## ... Procrustes: rmse 0.001182659  max resid 0.002437318 
    ## ... Similar to previous best
    ## Run 293 stress 0.018832 
    ## ... Procrustes: rmse 0.002305529  max resid 0.00508177 
    ## ... Similar to previous best
    ## Run 294 stress 0.01879535 
    ## ... New best solution
    ## ... Procrustes: rmse 4.70815e-06  max resid 9.281997e-06 
    ## ... Similar to previous best
    ## Run 295 stress 0.02520885 
    ## Run 296 stress 0.02509436 
    ## Run 297 stress 0.01881519 
    ## ... Procrustes: rmse 0.003027204  max resid 0.006234778 
    ## ... Similar to previous best
    ## Run 298 stress 0.01879983 
    ## ... Procrustes: rmse 0.001235013  max resid 0.002545086 
    ## ... Similar to previous best
    ## Run 299 stress 0.01896171 
    ## ... Procrustes: rmse 0.009615875  max resid 0.01973867 
    ## Run 300 stress 0.02509469 
    ## Run 301 stress 0.02509452 
    ## Run 302 stress 0.02520915 
    ## Run 303 stress 0.01879548 
    ## ... Procrustes: rmse 7.815654e-05  max resid 0.0001611102 
    ## ... Similar to previous best
    ## Run 304 stress 0.01880356 
    ## ... Procrustes: rmse 0.001800823  max resid 0.003711427 
    ## ... Similar to previous best
    ## Run 305 stress 0.01879569 
    ## ... Procrustes: rmse 0.0001812415  max resid 0.0003734972 
    ## ... Similar to previous best
    ## Run 306 stress 0.01879569 
    ## ... Procrustes: rmse 0.001115514  max resid 0.00229919 
    ## ... Similar to previous best
    ## Run 307 stress 0.02509462 
    ## Run 308 stress 0.01883676 
    ## ... Procrustes: rmse 0.002880375  max resid 0.005503255 
    ## ... Similar to previous best
    ## Run 309 stress 0.02492186 
    ## Run 310 stress 0.02509463 
    ## Run 311 stress 0.02486116 
    ## Run 312 stress 0.02509474 
    ## Run 313 stress 0.01883679 
    ## ... Procrustes: rmse 0.00287562  max resid 0.005499546 
    ## ... Similar to previous best
    ## Run 314 stress 0.224043 
    ## Run 315 stress 0.01894275 
    ## ... Procrustes: rmse 0.008665268  max resid 0.01760636 
    ## Run 316 stress 0.01879596 
    ## ... Procrustes: rmse 0.0002880997  max resid 0.0005935531 
    ## ... Similar to previous best
    ## Run 317 stress 0.01879608 
    ## ... Procrustes: rmse 0.0003219693  max resid 0.0006630066 
    ## ... Similar to previous best
    ## Run 318 stress 0.01890792 
    ## ... Procrustes: rmse 0.007258533  max resid 0.01468733 
    ## Run 319 stress 0.01879573 
    ## ... Procrustes: rmse 0.0002002127  max resid 0.0004125505 
    ## ... Similar to previous best
    ## Run 320 stress 0.01880887 
    ## ... Procrustes: rmse 0.002405239  max resid 0.004956645 
    ## ... Similar to previous best
    ## Run 321 stress 0.02509468 
    ## Run 322 stress 0.01879826 
    ## ... Procrustes: rmse 0.000904901  max resid 0.001865659 
    ## ... Similar to previous best
    ## Run 323 stress 0.01884114 
    ## ... Procrustes: rmse 0.003320915  max resid 0.006395238 
    ## ... Similar to previous best
    ## Run 324 stress 0.02509445 
    ## Run 325 stress 0.01879599 
    ## ... Procrustes: rmse 0.0003015471  max resid 0.0006211261 
    ## ... Similar to previous best
    ## Run 326 stress 0.02492197 
    ## Run 327 stress 0.01879583 
    ## ... Procrustes: rmse 0.000234151  max resid 0.0004823291 
    ## ... Similar to previous best
    ## Run 328 stress 0.02486127 
    ## Run 329 stress 0.01879568 
    ## ... Procrustes: rmse 0.0001782749  max resid 0.0003673875 
    ## ... Similar to previous best
    ## Run 330 stress 0.01879583 
    ## ... Procrustes: rmse 0.0002305779  max resid 0.0004750079 
    ## ... Similar to previous best
    ## Run 331 stress 0.01883724 
    ## ... Procrustes: rmse 0.00293895  max resid 0.005565008 
    ## ... Similar to previous best
    ## Run 332 stress 0.0250948 
    ## Run 333 stress 0.01879581 
    ## ... Procrustes: rmse 0.000233453  max resid 0.0004810139 
    ## ... Similar to previous best
    ## Run 334 stress 0.0250943 
    ## Run 335 stress 0.02492201 
    ## Run 336 stress 0.01883531 
    ## ... Procrustes: rmse 0.002723627  max resid 0.005398483 
    ## ... Similar to previous best
    ## Run 337 stress 0.0187985 
    ## ... Procrustes: rmse 0.000974679  max resid 0.002008506 
    ## ... Similar to previous best
    ## Run 338 stress 0.0250946 
    ## Run 339 stress 0.02520895 
    ## Run 340 stress 0.02509473 
    ## Run 341 stress 0.01879814 
    ## ... Procrustes: rmse 0.0009074343  max resid 0.001869644 
    ## ... Similar to previous best
    ## Run 342 stress 0.02509465 
    ## Run 343 stress 0.01879575 
    ## ... Procrustes: rmse 0.0001894191  max resid 0.0003900853 
    ## ... Similar to previous best
    ## Run 344 stress 0.01879558 
    ## ... Procrustes: rmse 0.0001258719  max resid 0.0002593933 
    ## ... Similar to previous best
    ## Run 345 stress 0.01879501 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0005264267  max resid 0.001087685 
    ## ... Similar to previous best
    ## Run 346 stress 0.02509443 
    ## Run 347 stress 0.02520908 
    ## Run 348 stress 0.02509467 
    ## Run 349 stress 0.02492187 
    ## Run 350 stress 0.024922 
    ## Run 351 stress 0.01879742 
    ## ... Procrustes: rmse 0.00123093  max resid 0.002539381 
    ## ... Similar to previous best
    ## Run 352 stress 0.02520898 
    ## Run 353 stress 0.01879578 
    ## ... Procrustes: rmse 0.0007493077  max resid 0.001546899 
    ## ... Similar to previous best
    ## Run 354 stress 0.01879568 
    ## ... Procrustes: rmse 0.0006970536  max resid 0.001439272 
    ## ... Similar to previous best
    ## Run 355 stress 0.01879592 
    ## ... Procrustes: rmse 0.0007992871  max resid 0.001649982 
    ## ... Similar to previous best
    ## Run 356 stress 0.02520908 
    ## Run 357 stress 0.0188694 
    ## ... Procrustes: rmse 0.00587692  max resid 0.01180681 
    ## Run 358 stress 0.01891287 
    ## ... Procrustes: rmse 0.008014874  max resid 0.01625733 
    ## Run 359 stress 0.0187958 
    ## ... Procrustes: rmse 0.0007561102  max resid 0.001560876 
    ## ... Similar to previous best
    ## Run 360 stress 0.02509475 
    ## Run 361 stress 0.02492175 
    ## Run 362 stress 0.01884031 
    ## ... Procrustes: rmse 0.003747924  max resid 0.007318247 
    ## ... Similar to previous best
    ## Run 363 stress 0.02492183 
    ## Run 364 stress 0.02492178 
    ## Run 365 stress 0.02509473 
    ## Run 366 stress 0.01879565 
    ## ... Procrustes: rmse 0.0006893543  max resid 0.001423366 
    ## ... Similar to previous best
    ## Run 367 stress 0.02486121 
    ## Run 368 stress 0.02520895 
    ## Run 369 stress 0.01879566 
    ## ... Procrustes: rmse 0.0005793502  max resid 0.00118634 
    ## ... Similar to previous best
    ## Run 370 stress 0.01945336 
    ## Run 371 stress 0.02486129 
    ## Run 372 stress 0.02486103 
    ## Run 373 stress 0.02509479 
    ## Run 374 stress 0.01881223 
    ## ... Procrustes: rmse 0.003280106  max resid 0.006759021 
    ## ... Similar to previous best
    ## Run 375 stress 0.01879588 
    ## ... Procrustes: rmse 0.0006744925  max resid 0.001382036 
    ## ... Similar to previous best
    ## Run 376 stress 0.02486127 
    ## Run 377 stress 0.01883732 
    ## ... Procrustes: rmse 0.003442102  max resid 0.006662741 
    ## ... Similar to previous best
    ## Run 378 stress 0.01881409 
    ## ... Procrustes: rmse 0.003455057  max resid 0.007118236 
    ## ... Similar to previous best
    ## Run 379 stress 0.01879556 
    ## ... Procrustes: rmse 0.0006423758  max resid 0.001326427 
    ## ... Similar to previous best
    ## Run 380 stress 0.01879606 
    ## ... Procrustes: rmse 0.0008435534  max resid 0.001741227 
    ## ... Similar to previous best
    ## Run 381 stress 0.01879601 
    ## ... Procrustes: rmse 0.0008394934  max resid 0.001732738 
    ## ... Similar to previous best
    ## Run 382 stress 0.01879555 
    ## ... Procrustes: rmse 0.0005275098  max resid 0.001079698 
    ## ... Similar to previous best
    ## Run 383 stress 0.02492196 
    ## Run 384 stress 0.01879887 
    ## ... Procrustes: rmse 0.001270831  max resid 0.00262116 
    ## ... Similar to previous best
    ## Run 385 stress 0.02492183 
    ## Run 386 stress 0.02486111 
    ## Run 387 stress 0.02486118 
    ## Run 388 stress 0.0250946 
    ## Run 389 stress 0.02493415 
    ## Run 390 stress 0.02492171 
    ## Run 391 stress 0.0249219 
    ## Run 392 stress 0.02509475 
    ## Run 393 stress 0.02509437 
    ## Run 394 stress 0.02509457 
    ## Run 395 stress 0.02486141 
    ## Run 396 stress 0.02509441 
    ## Run 397 stress 0.02509453 
    ## Run 398 stress 0.02492168 
    ## Run 399 stress 0.01883955 
    ## ... Procrustes: rmse 0.003633492  max resid 0.007072059 
    ## ... Similar to previous best
    ## Run 400 stress 0.02509449 
    ## Run 401 stress 0.02492201 
    ## Run 402 stress 0.02492199 
    ## Run 403 stress 0.02509461 
    ## Run 404 stress 0.02492172 
    ## Run 405 stress 0.02025599 
    ## Run 406 stress 0.01897562 
    ## ... Procrustes: rmse 0.01060327  max resid 0.02175102 
    ## Run 407 stress 0.01905173 
    ## ... Procrustes: rmse 0.01254712  max resid 0.02560478 
    ## Run 408 stress 0.02509461 
    ## Run 409 stress 0.01879592 
    ## ... Procrustes: rmse 0.0008040447  max resid 0.00165967 
    ## ... Similar to previous best
    ## Run 410 stress 0.02492181 
    ## Run 411 stress 0.01881817 
    ## ... Procrustes: rmse 0.002072851  max resid 0.004034448 
    ## ... Similar to previous best
    ## Run 412 stress 0.02486101 
    ## Run 413 stress 0.0187956 
    ## ... Procrustes: rmse 0.0005516469  max resid 0.001129373 
    ## ... Similar to previous best
    ## Run 414 stress 0.01883237 
    ## ... Procrustes: rmse 0.002838486  max resid 0.005402714 
    ## ... Similar to previous best
    ## Run 415 stress 0.01879587 
    ## ... Procrustes: rmse 0.000670373  max resid 0.001373568 
    ## ... Similar to previous best
    ## Run 416 stress 0.01879578 
    ## ... Procrustes: rmse 0.0007472237  max resid 0.001542633 
    ## ... Similar to previous best
    ## Run 417 stress 0.02509436 
    ## Run 418 stress 0.01879895 
    ## ... Procrustes: rmse 0.001577907  max resid 0.003257485 
    ## ... Similar to previous best
    ## Run 419 stress 0.02509447 
    ## Run 420 stress 0.01882428 
    ## ... Procrustes: rmse 0.004273187  max resid 0.008798626 
    ## ... Similar to previous best
    ## Run 421 stress 0.02492175 
    ## Run 422 stress 0.02509476 
    ## Run 423 stress 0.02486115 
    ## Run 424 stress 0.02486103 
    ## Run 425 stress 0.02497138 
    ## Run 426 stress 0.01879967 
    ## ... Procrustes: rmse 0.001729922  max resid 0.003567526 
    ## ... Similar to previous best
    ## Run 427 stress 0.02492184 
    ## Run 428 stress 0.01890849 
    ## ... Procrustes: rmse 0.007542655  max resid 0.01527012 
    ## Run 429 stress 0.02509468 
    ## Run 430 stress 0.02492171 
    ## Run 431 stress 0.02520904 
    ## Run 432 stress 0.01884178 
    ## ... Procrustes: rmse 0.003746563  max resid 0.007315777 
    ## ... Similar to previous best
    ## Run 433 stress 0.02492204 
    ## Run 434 stress 0.01887727 
    ## ... Procrustes: rmse 0.006321508  max resid 0.01273478 
    ## Run 435 stress 0.02509475 
    ## Run 436 stress 0.0188359 
    ## ... Procrustes: rmse 0.003282473  max resid 0.006320268 
    ## ... Similar to previous best
    ## Run 437 stress 0.01884159 
    ## ... Procrustes: rmse 0.005386378  max resid 0.01108197 
    ## Run 438 stress 0.0187955 
    ## ... Procrustes: rmse 0.0005002805  max resid 0.001023678 
    ## ... Similar to previous best
    ## Run 439 stress 0.02492192 
    ## Run 440 stress 0.02509463 
    ## Run 441 stress 0.01892334 
    ## ... Procrustes: rmse 0.008367355  max resid 0.01698803 
    ## Run 442 stress 0.01879711 
    ## ... Procrustes: rmse 0.0009542854  max resid 0.001968764 
    ## ... Similar to previous best
    ## Run 443 stress 0.02486118 
    ## Run 444 stress 0.01889146 
    ## ... Procrustes: rmse 0.007049129  max resid 0.01425075 
    ## Run 445 stress 0.0250945 
    ## Run 446 stress 0.01882039 
    ## ... Procrustes: rmse 0.003917949  max resid 0.008069954 
    ## ... Similar to previous best
    ## Run 447 stress 0.02509447 
    ## Run 448 stress 0.01883091 
    ## ... Procrustes: rmse 0.00263216  max resid 0.00522559 
    ## ... Similar to previous best
    ## Run 449 stress 0.01879581 
    ## ... Procrustes: rmse 0.0007547594  max resid 0.001557999 
    ## ... Similar to previous best
    ## Run 450 stress 0.01880407 
    ## ... Procrustes: rmse 0.00227295  max resid 0.004684219 
    ## ... Similar to previous best
    ## Run 451 stress 0.02492197 
    ## Run 452 stress 0.01879598 
    ## ... Procrustes: rmse 0.0008219137  max resid 0.001696612 
    ## ... Similar to previous best
    ## Run 453 stress 0.02492195 
    ## Run 454 stress 0.224043 
    ## Run 455 stress 0.0187958 
    ## ... Procrustes: rmse 0.0007544966  max resid 0.001557605 
    ## ... Similar to previous best
    ## Run 456 stress 0.02509467 
    ## Run 457 stress 0.02509449 
    ## Run 458 stress 0.0198496 
    ## Run 459 stress 0.01879582 
    ## ... Procrustes: rmse 0.0006517078  max resid 0.001335193 
    ## ... Similar to previous best
    ## Run 460 stress 0.01879573 
    ## ... Procrustes: rmse 0.0007217755  max resid 0.001490228 
    ## ... Similar to previous best
    ## Run 461 stress 0.01879579 
    ## ... Procrustes: rmse 0.0007518906  max resid 0.001552157 
    ## ... Similar to previous best
    ## Run 462 stress 0.0187962 
    ## ... Procrustes: rmse 0.0009080685  max resid 0.001873911 
    ## ... Similar to previous best
    ## Run 463 stress 0.01883747 
    ## ... Procrustes: rmse 0.00344382  max resid 0.006668608 
    ## ... Similar to previous best
    ## Run 464 stress 0.02509441 
    ## Run 465 stress 0.01880839 
    ## ... Procrustes: rmse 0.00289665  max resid 0.005970462 
    ## ... Similar to previous best
    ## Run 466 stress 0.02492184 
    ## Run 467 stress 0.01879558 
    ## ... Procrustes: rmse 0.0006531206  max resid 0.001348621 
    ## ... Similar to previous best
    ## Run 468 stress 0.02492193 
    ## Run 469 stress 0.02486131 
    ## Run 470 stress 0.02509444 
    ## Run 471 stress 0.02486139 
    ## Run 472 stress 0.02486121 
    ## Run 473 stress 0.01879593 
    ## ... Procrustes: rmse 0.0008049887  max resid 0.001661533 
    ## ... Similar to previous best
    ## Run 474 stress 0.02509465 
    ## Run 475 stress 0.02509461 
    ## Run 476 stress 0.01879573 
    ## ... Procrustes: rmse 0.0006069671  max resid 0.00124311 
    ## ... Similar to previous best
    ## Run 477 stress 0.02509448 
    ## Run 478 stress 0.02520897 
    ## Run 479 stress 0.01879575 
    ## ... Procrustes: rmse 0.0007356405  max resid 0.001518757 
    ## ... Similar to previous best
    ## Run 480 stress 0.02509454 
    ## Run 481 stress 0.01879606 
    ## ... Procrustes: rmse 0.0008463631  max resid 0.001746675 
    ## ... Similar to previous best
    ## Run 482 stress 0.0190477 
    ## ... Procrustes: rmse 0.01243991  max resid 0.02539899 
    ## Run 483 stress 0.02509459 
    ## Run 484 stress 0.01886464 
    ## ... Procrustes: rmse 0.005484688  max resid 0.0109822 
    ## Run 485 stress 0.02492184 
    ## Run 486 stress 0.02520916 
    ## Run 487 stress 0.02509465 
    ## Run 488 stress 0.02509469 
    ## Run 489 stress 0.01881694 
    ## ... Procrustes: rmse 0.003612986  max resid 0.007439135 
    ## ... Similar to previous best
    ## Run 490 stress 0.02486135 
    ## Run 491 stress 0.02492194 
    ## Run 492 stress 0.01879588 
    ## ... Procrustes: rmse 0.0007905706  max resid 0.001631958 
    ## ... Similar to previous best
    ## Run 493 stress 0.024922 
    ## Run 494 stress 0.02683512 
    ## Run 495 stress 0.02509458 
    ## Run 496 stress 0.01879598 
    ## ... Procrustes: rmse 0.0008195506  max resid 0.001691476 
    ## ... Similar to previous best
    ## Run 497 stress 0.02509474 
    ## Run 498 stress 0.02492184 
    ## Run 499 stress 0.01879649 
    ## ... Procrustes: rmse 0.001004197  max resid 0.002071996 
    ## ... Similar to previous best
    ## Run 500 stress 0.02486125 
    ## *** Best solution repeated 53 times

``` r
# Mixed and stratified lakes
SD_beta_MS_NMDS <- metaMDS(SD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09337244 
    ## Run 2 stress 0.08773472 
    ## Run 3 stress 0.09539207 
    ## Run 4 stress 0.09590631 
    ## Run 5 stress 0.09030393 
    ## Run 6 stress 0.09407967 
    ## Run 7 stress 0.08503631 
    ## ... Procrustes: rmse 0.0001665063  max resid 0.0003846985 
    ## ... Similar to previous best
    ## Run 8 stress 0.09168948 
    ## Run 9 stress 0.09030395 
    ## Run 10 stress 0.1026215 
    ## Run 11 stress 0.08973869 
    ## Run 12 stress 0.09417785 
    ## Run 13 stress 0.0940802 
    ## Run 14 stress 0.09286084 
    ## Run 15 stress 0.1052851 
    ## Run 16 stress 0.08773486 
    ## Run 17 stress 0.08973864 
    ## Run 18 stress 0.09168942 
    ## Run 19 stress 0.09145338 
    ## Run 20 stress 0.08973866 
    ## Run 21 stress 0.09969984 
    ## Run 22 stress 0.09721199 
    ## Run 23 stress 0.09407984 
    ## Run 24 stress 0.09908184 
    ## Run 25 stress 0.09445514 
    ## Run 26 stress 0.09539201 
    ## Run 27 stress 0.08773467 
    ## Run 28 stress 0.09337242 
    ## Run 29 stress 0.09969974 
    ## Run 30 stress 0.2950554 
    ## Run 31 stress 0.09590791 
    ## Run 32 stress 0.08503503 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001774933  max resid 0.004426589 
    ## ... Similar to previous best
    ## Run 33 stress 0.09337262 
    ## Run 34 stress 0.09407988 
    ## Run 35 stress 0.08503507 
    ## ... Procrustes: rmse 2.12791e-05  max resid 4.698781e-05 
    ## ... Similar to previous best
    ## Run 36 stress 0.08503508 
    ## ... Procrustes: rmse 3.17736e-05  max resid 8.504146e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.09407974 
    ## Run 38 stress 0.2567616 
    ## Run 39 stress 0.08773469 
    ## Run 40 stress 0.09168936 
    ## Run 41 stress 0.09030401 
    ## Run 42 stress 0.08973871 
    ## Run 43 stress 0.09030395 
    ## Run 44 stress 0.09159083 
    ## Run 45 stress 0.08503509 
    ## ... Procrustes: rmse 3.011159e-05  max resid 5.851224e-05 
    ## ... Similar to previous best
    ## Run 46 stress 0.09159085 
    ## Run 47 stress 0.09337249 
    ## Run 48 stress 0.09417798 
    ## Run 49 stress 0.09030401 
    ## Run 50 stress 0.2574363 
    ## Run 51 stress 0.08773476 
    ## Run 52 stress 0.08973869 
    ## Run 53 stress 0.09721197 
    ## Run 54 stress 0.09535556 
    ## Run 55 stress 0.09403413 
    ## Run 56 stress 0.09145329 
    ## Run 57 stress 0.0946445 
    ## Run 58 stress 0.09159085 
    ## Run 59 stress 0.245554 
    ## Run 60 stress 0.09268402 
    ## Run 61 stress 0.08973892 
    ## Run 62 stress 0.0928609 
    ## Run 63 stress 0.09030406 
    ## Run 64 stress 0.08503658 
    ## ... Procrustes: rmse 0.001997722  max resid 0.0048612 
    ## ... Similar to previous best
    ## Run 65 stress 0.09268335 
    ## Run 66 stress 0.08973886 
    ## Run 67 stress 0.09539207 
    ## Run 68 stress 0.09145321 
    ## Run 69 stress 0.09159085 
    ## Run 70 stress 0.08440271 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01419619  max resid 0.04323321 
    ## Run 71 stress 0.08503533 
    ## Run 72 stress 0.1038042 
    ## Run 73 stress 0.09464474 
    ## Run 74 stress 0.08773465 
    ## Run 75 stress 0.1026216 
    ## Run 76 stress 0.09969976 
    ## Run 77 stress 0.08503458 
    ## Run 78 stress 0.09168951 
    ## Run 79 stress 0.09030416 
    ## Run 80 stress 0.1072864 
    ## Run 81 stress 0.09168945 
    ## Run 82 stress 0.09539205 
    ## Run 83 stress 0.09374279 
    ## Run 84 stress 0.0953551 
    ## Run 85 stress 0.09669858 
    ## Run 86 stress 0.08503472 
    ## Run 87 stress 0.09374275 
    ## Run 88 stress 0.09268415 
    ## Run 89 stress 0.0996997 
    ## Run 90 stress 0.09268355 
    ## Run 91 stress 0.09407976 
    ## Run 92 stress 0.09286083 
    ## Run 93 stress 0.09407434 
    ## Run 94 stress 0.09145324 
    ## Run 95 stress 0.09407979 
    ## Run 96 stress 0.09145328 
    ## Run 97 stress 0.09403419 
    ## Run 98 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002096588  max resid 0.0004274813 
    ## ... Similar to previous best
    ## Run 99 stress 0.0946591 
    ## Run 100 stress 0.09030395 
    ## Run 101 stress 0.09381853 
    ## Run 102 stress 0.09407983 
    ## Run 103 stress 0.09721199 
    ## Run 104 stress 0.09030415 
    ## Run 105 stress 0.08440267 
    ## ... Procrustes: rmse 0.0002735801  max resid 0.0005500995 
    ## ... Similar to previous best
    ## Run 106 stress 0.08440261 
    ## ... Procrustes: rmse 0.0001634022  max resid 0.000370657 
    ## ... Similar to previous best
    ## Run 107 stress 0.09308978 
    ## Run 108 stress 0.09286089 
    ## Run 109 stress 0.08773473 
    ## Run 110 stress 0.0897387 
    ## Run 111 stress 0.1038077 
    ## Run 112 stress 0.09465905 
    ## Run 113 stress 0.09407473 
    ## Run 114 stress 0.09535498 
    ## Run 115 stress 0.09445481 
    ## Run 116 stress 0.09337242 
    ## Run 117 stress 0.09407965 
    ## Run 118 stress 0.08503556 
    ## Run 119 stress 0.0897387 
    ## Run 120 stress 0.09030405 
    ## Run 121 stress 0.08440272 
    ## ... Procrustes: rmse 0.0002445681  max resid 0.0004353752 
    ## ... Similar to previous best
    ## Run 122 stress 0.09030394 
    ## Run 123 stress 0.08973864 
    ## Run 124 stress 0.08973874 
    ## Run 125 stress 0.09539194 
    ## Run 126 stress 0.09030406 
    ## Run 127 stress 0.0877347 
    ## Run 128 stress 0.09030404 
    ## Run 129 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001083965  max resid 0.0002121338 
    ## ... Similar to previous best
    ## Run 130 stress 0.09145323 
    ## Run 131 stress 0.09286083 
    ## Run 132 stress 0.08440264 
    ## ... Procrustes: rmse 0.0001697302  max resid 0.0003050418 
    ## ... Similar to previous best
    ## Run 133 stress 0.09445465 
    ## Run 134 stress 0.08773466 
    ## Run 135 stress 0.09268361 
    ## Run 136 stress 0.09381868 
    ## Run 137 stress 0.09286085 
    ## Run 138 stress 0.09445514 
    ## Run 139 stress 0.09308979 
    ## Run 140 stress 0.09407127 
    ## Run 141 stress 0.09168945 
    ## Run 142 stress 0.09465918 
    ## Run 143 stress 0.0940799 
    ## Run 144 stress 0.08773465 
    ## Run 145 stress 0.09145322 
    ## Run 146 stress 0.09268382 
    ## Run 147 stress 0.0937423 
    ## Run 148 stress 0.09268341 
    ## Run 149 stress 0.0877347 
    ## Run 150 stress 0.09337253 
    ## Run 151 stress 0.1052851 
    ## Run 152 stress 0.09381845 
    ## Run 153 stress 0.09407965 
    ## Run 154 stress 0.08973862 
    ## Run 155 stress 0.09145322 
    ## Run 156 stress 0.09168956 
    ## Run 157 stress 0.08503479 
    ## Run 158 stress 0.08973868 
    ## Run 159 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001416539  max resid 0.000290276 
    ## ... Similar to previous best
    ## Run 160 stress 0.09760823 
    ## Run 161 stress 0.08503626 
    ## Run 162 stress 0.08503485 
    ## Run 163 stress 0.09030396 
    ## Run 164 stress 0.100461 
    ## Run 165 stress 0.09407976 
    ## Run 166 stress 0.09977543 
    ## Run 167 stress 0.09407143 
    ## Run 168 stress 0.095392 
    ## Run 169 stress 0.09374252 
    ## Run 170 stress 0.09374279 
    ## Run 171 stress 0.09969967 
    ## Run 172 stress 0.09407129 
    ## Run 173 stress 0.09539193 
    ## Run 174 stress 0.08440256 
    ## ... Procrustes: rmse 3.550379e-05  max resid 6.576815e-05 
    ## ... Similar to previous best
    ## Run 175 stress 0.08440252 
    ## ... Procrustes: rmse 4.098408e-05  max resid 8.415367e-05 
    ## ... Similar to previous best
    ## Run 176 stress 0.09416132 
    ## Run 177 stress 0.103804 
    ## Run 178 stress 0.09159083 
    ## Run 179 stress 0.09030404 
    ## Run 180 stress 0.09337304 
    ## Run 181 stress 0.09159084 
    ## Run 182 stress 0.09407996 
    ## Run 183 stress 0.09168958 
    ## Run 184 stress 0.09416134 
    ## Run 185 stress 0.1052855 
    ## Run 186 stress 0.08503507 
    ## Run 187 stress 0.09969979 
    ## Run 188 stress 0.09268325 
    ## Run 189 stress 0.09268353 
    ## Run 190 stress 0.09268371 
    ## Run 191 stress 0.08503477 
    ## Run 192 stress 0.08503638 
    ## Run 193 stress 0.08773487 
    ## Run 194 stress 0.09374313 
    ## Run 195 stress 0.09321518 
    ## Run 196 stress 0.09308987 
    ## Run 197 stress 0.09159095 
    ## Run 198 stress 0.09407992 
    ## Run 199 stress 0.09539212 
    ## Run 200 stress 0.09465914 
    ## Run 201 stress 0.09407992 
    ## Run 202 stress 0.09465909 
    ## Run 203 stress 0.08503504 
    ## Run 204 stress 0.09380577 
    ## Run 205 stress 0.08440252 
    ## ... Procrustes: rmse 4.407493e-05  max resid 8.197116e-05 
    ## ... Similar to previous best
    ## Run 206 stress 0.09286084 
    ## Run 207 stress 0.09590673 
    ## Run 208 stress 0.09407968 
    ## Run 209 stress 0.09030404 
    ## Run 210 stress 0.09445503 
    ## Run 211 stress 0.09447306 
    ## Run 212 stress 0.09400525 
    ## Run 213 stress 0.09168961 
    ## Run 214 stress 0.09145322 
    ## Run 215 stress 0.09417817 
    ## Run 216 stress 0.0926834 
    ## Run 217 stress 0.09286083 
    ## Run 218 stress 0.08503477 
    ## Run 219 stress 0.09308978 
    ## Run 220 stress 0.08773477 
    ## Run 221 stress 0.09159083 
    ## Run 222 stress 0.09159099 
    ## Run 223 stress 0.09465905 
    ## Run 224 stress 0.09590932 
    ## Run 225 stress 0.08503572 
    ## Run 226 stress 0.09337261 
    ## Run 227 stress 0.09465917 
    ## Run 228 stress 0.0877347 
    ## Run 229 stress 0.09159083 
    ## Run 230 stress 0.09407142 
    ## Run 231 stress 0.09407965 
    ## Run 232 stress 0.09416133 
    ## Run 233 stress 0.09416157 
    ## Run 234 stress 0.09445471 
    ## Run 235 stress 0.09268394 
    ## Run 236 stress 0.09030399 
    ## Run 237 stress 0.08503493 
    ## Run 238 stress 0.08503578 
    ## Run 239 stress 0.09374375 
    ## Run 240 stress 0.0941614 
    ## Run 241 stress 0.09408016 
    ## Run 242 stress 0.09464478 
    ## Run 243 stress 0.09268387 
    ## Run 244 stress 0.1026215 
    ## Run 245 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001785266  max resid 0.0003568953 
    ## ... Similar to previous best
    ## Run 246 stress 0.08503457 
    ## Run 247 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 2.904328e-05  max resid 7.189795e-05 
    ## ... Similar to previous best
    ## Run 248 stress 0.09407119 
    ## Run 249 stress 0.09465909 
    ## Run 250 stress 0.09407965 
    ## Run 251 stress 0.09030408 
    ## Run 252 stress 0.09159093 
    ## Run 253 stress 0.09030396 
    ## Run 254 stress 0.08973864 
    ## Run 255 stress 0.09404381 
    ## Run 256 stress 0.08773479 
    ## Run 257 stress 0.08973865 
    ## Run 258 stress 0.08503508 
    ## Run 259 stress 0.0850347 
    ## Run 260 stress 0.09030396 
    ## Run 261 stress 0.09712931 
    ## Run 262 stress 0.09337315 
    ## Run 263 stress 0.0940799 
    ## Run 264 stress 0.09286088 
    ## Run 265 stress 0.09268338 
    ## Run 266 stress 0.09337257 
    ## Run 267 stress 0.08773476 
    ## Run 268 stress 0.08973879 
    ## Run 269 stress 0.08973867 
    ## Run 270 stress 0.08440261 
    ## ... Procrustes: rmse 0.0001420203  max resid 0.0002396478 
    ## ... Similar to previous best
    ## Run 271 stress 0.0941777 
    ## Run 272 stress 0.08503486 
    ## Run 273 stress 0.09268342 
    ## Run 274 stress 0.09465905 
    ## Run 275 stress 0.09145325 
    ## Run 276 stress 0.0938058 
    ## Run 277 stress 0.08440267 
    ## ... Procrustes: rmse 0.0002330507  max resid 0.0003869669 
    ## ... Similar to previous best
    ## Run 278 stress 0.1052851 
    ## Run 279 stress 0.08973867 
    ## Run 280 stress 0.09407993 
    ## Run 281 stress 0.0946591 
    ## Run 282 stress 0.08503455 
    ## Run 283 stress 0.1124502 
    ## Run 284 stress 0.1038032 
    ## Run 285 stress 0.09969992 
    ## Run 286 stress 0.0877347 
    ## Run 287 stress 0.09465906 
    ## Run 288 stress 0.08503553 
    ## Run 289 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001461724  max resid 0.0003186347 
    ## ... Similar to previous best
    ## Run 290 stress 0.08503479 
    ## Run 291 stress 0.09030412 
    ## Run 292 stress 0.08773466 
    ## Run 293 stress 0.09381845 
    ## Run 294 stress 0.08503625 
    ## Run 295 stress 0.09539205 
    ## Run 296 stress 0.09465909 
    ## Run 297 stress 0.0933728 
    ## Run 298 stress 0.09465904 
    ## Run 299 stress 0.08973862 
    ## Run 300 stress 0.09447311 
    ## Run 301 stress 0.09286084 
    ## Run 302 stress 0.09030394 
    ## Run 303 stress 0.09407966 
    ## Run 304 stress 0.09030411 
    ## Run 305 stress 0.2663943 
    ## Run 306 stress 0.09159083 
    ## Run 307 stress 0.08503477 
    ## Run 308 stress 0.09268337 
    ## Run 309 stress 0.09407986 
    ## Run 310 stress 0.09286087 
    ## Run 311 stress 0.09465913 
    ## Run 312 stress 0.09308956 
    ## Run 313 stress 0.08503486 
    ## Run 314 stress 0.08503486 
    ## Run 315 stress 0.08773468 
    ## Run 316 stress 0.09407965 
    ## Run 317 stress 0.09465911 
    ## Run 318 stress 0.09969987 
    ## Run 319 stress 0.09445479 
    ## Run 320 stress 0.09030409 
    ## Run 321 stress 0.09407108 
    ## Run 322 stress 0.08503694 
    ## Run 323 stress 0.08773475 
    ## Run 324 stress 0.0850346 
    ## Run 325 stress 0.09407151 
    ## Run 326 stress 0.08773466 
    ## Run 327 stress 0.09168937 
    ## Run 328 stress 0.2456114 
    ## Run 329 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001211116  max resid 0.0002662112 
    ## ... Similar to previous best
    ## Run 330 stress 0.09145324 
    ## Run 331 stress 0.09465905 
    ## Run 332 stress 0.09030401 
    ## Run 333 stress 0.09374276 
    ## Run 334 stress 0.08973864 
    ## Run 335 stress 0.09969981 
    ## Run 336 stress 0.09030393 
    ## Run 337 stress 0.09145326 
    ## Run 338 stress 0.09669861 
    ## Run 339 stress 0.09465906 
    ## Run 340 stress 0.09445473 
    ## Run 341 stress 0.0903041 
    ## Run 342 stress 0.09464457 
    ## Run 343 stress 0.08440267 
    ## ... Procrustes: rmse 0.0002408982  max resid 0.0005135312 
    ## ... Similar to previous best
    ## Run 344 stress 0.08973882 
    ## Run 345 stress 0.08503462 
    ## Run 346 stress 0.1038039 
    ## Run 347 stress 0.09380577 
    ## Run 348 stress 0.08973887 
    ## Run 349 stress 0.09969976 
    ## Run 350 stress 0.09159091 
    ## Run 351 stress 0.09404395 
    ## Run 352 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001107112  max resid 0.0003119414 
    ## ... Similar to previous best
    ## Run 353 stress 0.0938058 
    ## Run 354 stress 0.09308976 
    ## Run 355 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001282554  max resid 0.0002367468 
    ## ... Similar to previous best
    ## Run 356 stress 0.09030411 
    ## Run 357 stress 0.09760814 
    ## Run 358 stress 0.09539202 
    ## Run 359 stress 0.0877347 
    ## Run 360 stress 0.09417777 
    ## Run 361 stress 0.09374186 
    ## Run 362 stress 0.0850348 
    ## Run 363 stress 0.08503517 
    ## Run 364 stress 0.09308974 
    ## Run 365 stress 0.09030396 
    ## Run 366 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001101179  max resid 0.0002015994 
    ## ... Similar to previous best
    ## Run 367 stress 0.08773481 
    ## Run 368 stress 0.08440274 
    ## ... Procrustes: rmse 0.0002817007  max resid 0.0004682902 
    ## ... Similar to previous best
    ## Run 369 stress 0.09030394 
    ## Run 370 stress 0.09030398 
    ## Run 371 stress 0.09159085 
    ## Run 372 stress 0.09721197 
    ## Run 373 stress 0.08973875 
    ## Run 374 stress 0.09268333 
    ## Run 375 stress 0.09407997 
    ## Run 376 stress 0.09337259 
    ## Run 377 stress 0.09408005 
    ## Run 378 stress 0.09408018 
    ## Run 379 stress 0.09465904 
    ## Run 380 stress 0.09286086 
    ## Run 381 stress 0.08503483 
    ## Run 382 stress 0.08773478 
    ## Run 383 stress 0.09465915 
    ## Run 384 stress 0.09308978 
    ## Run 385 stress 0.09407997 
    ## Run 386 stress 0.09374348 
    ## Run 387 stress 0.09445477 
    ## Run 388 stress 0.09465912 
    ## Run 389 stress 0.09286089 
    ## Run 390 stress 0.09168939 
    ## Run 391 stress 0.08503472 
    ## Run 392 stress 0.09381845 
    ## Run 393 stress 0.09145325 
    ## Run 394 stress 0.09159086 
    ## Run 395 stress 0.09145331 
    ## Run 396 stress 0.09168944 
    ## Run 397 stress 0.09380583 
    ## Run 398 stress 0.08503474 
    ## Run 399 stress 0.09407971 
    ## Run 400 stress 0.0916893 
    ## Run 401 stress 0.09539204 
    ## Run 402 stress 0.08440253 
    ## ... Procrustes: rmse 8.693336e-05  max resid 0.0001555412 
    ## ... Similar to previous best
    ## Run 403 stress 0.09465914 
    ## Run 404 stress 0.1052849 
    ## Run 405 stress 0.09416147 
    ## Run 406 stress 0.08973884 
    ## Run 407 stress 0.08973899 
    ## Run 408 stress 0.09447314 
    ## Run 409 stress 0.09030404 
    ## Run 410 stress 0.0926834 
    ## Run 411 stress 0.09308954 
    ## Run 412 stress 0.09969992 
    ## Run 413 stress 0.09465909 
    ## Run 414 stress 0.09268336 
    ## Run 415 stress 0.09159084 
    ## Run 416 stress 0.0940797 
    ## Run 417 stress 0.08503483 
    ## Run 418 stress 0.09286091 
    ## Run 419 stress 0.09030413 
    ## Run 420 stress 0.0850346 
    ## Run 421 stress 0.08973863 
    ## Run 422 stress 0.09908183 
    ## Run 423 stress 0.09407974 
    ## Run 424 stress 0.09030406 
    ## Run 425 stress 0.09380584 
    ## Run 426 stress 0.09407979 
    ## Run 427 stress 0.09535432 
    ## Run 428 stress 0.09337236 
    ## Run 429 stress 0.0850359 
    ## Run 430 stress 0.100461 
    ## Run 431 stress 0.09337294 
    ## Run 432 stress 0.1052851 
    ## Run 433 stress 0.09159085 
    ## Run 434 stress 0.09159089 
    ## Run 435 stress 0.09465912 
    ## Run 436 stress 0.09760819 
    ## Run 437 stress 0.08503578 
    ## Run 438 stress 0.09969976 
    ## Run 439 stress 0.09407988 
    ## Run 440 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001214693  max resid 0.0002173051 
    ## ... Similar to previous best
    ## Run 441 stress 0.09268345 
    ## Run 442 stress 0.09465904 
    ## Run 443 stress 0.09407999 
    ## Run 444 stress 0.09286083 
    ## Run 445 stress 0.09159084 
    ## Run 446 stress 0.09407997 
    ## Run 447 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001201023  max resid 0.0002204082 
    ## ... Similar to previous best
    ## Run 448 stress 0.08973862 
    ## Run 449 stress 0.09308971 
    ## Run 450 stress 0.09286086 
    ## Run 451 stress 0.09159085 
    ## Run 452 stress 0.09286083 
    ## Run 453 stress 0.09535597 
    ## Run 454 stress 0.09465904 
    ## Run 455 stress 0.09030394 
    ## Run 456 stress 0.09539202 
    ## Run 457 stress 0.1005995 
    ## Run 458 stress 0.09535448 
    ## Run 459 stress 0.0996997 
    ## Run 460 stress 0.08503568 
    ## Run 461 stress 0.1038041 
    ## Run 462 stress 0.09145328 
    ## Run 463 stress 0.09374273 
    ## Run 464 stress 0.09286085 
    ## Run 465 stress 0.08440252 
    ## ... Procrustes: rmse 3.200145e-05  max resid 7.646624e-05 
    ## ... Similar to previous best
    ## Run 466 stress 0.08773468 
    ## Run 467 stress 0.09465906 
    ## Run 468 stress 0.09407983 
    ## Run 469 stress 0.08773477 
    ## Run 470 stress 0.08973865 
    ## Run 471 stress 0.09416148 
    ## Run 472 stress 0.09969979 
    ## Run 473 stress 0.09286083 
    ## Run 474 stress 0.09539198 
    ## Run 475 stress 0.08503467 
    ## Run 476 stress 0.0850366 
    ## Run 477 stress 0.1004611 
    ## Run 478 stress 0.09030399 
    ## Run 479 stress 0.08503514 
    ## Run 480 stress 0.09407405 
    ## Run 481 stress 0.08973887 
    ## Run 482 stress 0.09145321 
    ## Run 483 stress 0.09168934 
    ## Run 484 stress 0.09168942 
    ## Run 485 stress 0.09417769 
    ## Run 486 stress 0.09760832 
    ## Run 487 stress 0.09908181 
    ## Run 488 stress 0.09535495 
    ## Run 489 stress 0.08973879 
    ## Run 490 stress 0.2488816 
    ## Run 491 stress 0.0937424 
    ## Run 492 stress 0.09159087 
    ## Run 493 stress 0.09445514 
    ## Run 494 stress 0.09030397 
    ## Run 495 stress 0.1038042 
    ## Run 496 stress 0.08503574 
    ## Run 497 stress 0.08973872 
    ## Run 498 stress 0.09168957 
    ## Run 499 stress 0.09407372 
    ## Run 500 stress 0.08773467 
    ## *** Best solution repeated 14 times

``` r
# Ocean sites and mixed lakes
SD_beta_OM_NMDS <- metaMDS(SD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.07365792 
    ## ... Procrustes: rmse 0.000145991  max resid 0.0003474959 
    ## ... Similar to previous best
    ## Run 2 stress 0.07365785 
    ## ... Procrustes: rmse 6.010115e-05  max resid 0.0001430031 
    ## ... Similar to previous best
    ## Run 3 stress 0.07732926 
    ## Run 4 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744815  max resid 0.0510522 
    ## Run 5 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744489  max resid 0.0510263 
    ## Run 6 stress 0.07378231 
    ## ... Procrustes: rmse 0.0174382  max resid 0.05104709 
    ## Run 7 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744131  max resid 0.05105024 
    ## Run 8 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745755  max resid 0.05107937 
    ## Run 9 stress 0.07629237 
    ## Run 10 stress 0.08288045 
    ## Run 11 stress 0.0762924 
    ## Run 12 stress 0.08000467 
    ## Run 13 stress 0.07732934 
    ## Run 14 stress 0.08288036 
    ## Run 15 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001065638  max resid 0.0002507206 
    ## ... Similar to previous best
    ## Run 16 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744517  max resid 0.05105387 
    ## Run 17 stress 0.0800047 
    ## Run 18 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744238  max resid 0.05102623 
    ## Run 19 stress 0.07365785 
    ## ... Procrustes: rmse 6.560353e-05  max resid 0.000150859 
    ## ... Similar to previous best
    ## Run 20 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745827  max resid 0.05111814 
    ## Run 21 stress 0.08000472 
    ## Run 22 stress 0.08000474 
    ## Run 23 stress 0.08000474 
    ## Run 24 stress 0.07365783 
    ## ... Procrustes: rmse 6.498558e-06  max resid 1.420445e-05 
    ## ... Similar to previous best
    ## Run 25 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743383  max resid 0.05099222 
    ## Run 26 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744502  max resid 0.05104682 
    ## Run 27 stress 0.07629237 
    ## Run 28 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744546  max resid 0.05104975 
    ## Run 29 stress 0.07629239 
    ## Run 30 stress 0.07629239 
    ## Run 31 stress 0.08288059 
    ## Run 32 stress 0.07629233 
    ## Run 33 stress 0.0800047 
    ## Run 34 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744989  max resid 0.05106654 
    ## Run 35 stress 0.07732937 
    ## Run 36 stress 0.07732928 
    ## Run 37 stress 0.07365788 
    ## ... Procrustes: rmse 0.000106878  max resid 0.0002546623 
    ## ... Similar to previous best
    ## Run 38 stress 0.07378229 
    ## ... Procrustes: rmse 0.01743992  max resid 0.05099794 
    ## Run 39 stress 0.08000469 
    ## Run 40 stress 0.0762924 
    ## Run 41 stress 0.07629234 
    ## Run 42 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743014  max resid 0.05098205 
    ## Run 43 stress 0.07629233 
    ## Run 44 stress 0.07629237 
    ## Run 45 stress 0.07629239 
    ## Run 46 stress 0.08288035 
    ## Run 47 stress 0.07629234 
    ## Run 48 stress 0.07629236 
    ## Run 49 stress 0.07365784 
    ## ... Procrustes: rmse 4.481112e-05  max resid 0.0001041842 
    ## ... Similar to previous best
    ## Run 50 stress 0.07365783 
    ## ... Procrustes: rmse 3.275451e-06  max resid 6.626655e-06 
    ## ... Similar to previous best
    ## Run 51 stress 0.07629235 
    ## Run 52 stress 0.07732948 
    ## Run 53 stress 0.07629241 
    ## Run 54 stress 0.2967991 
    ## Run 55 stress 0.07365786 
    ## ... Procrustes: rmse 7.751811e-05  max resid 0.0001763995 
    ## ... Similar to previous best
    ## Run 56 stress 0.07365784 
    ## ... Procrustes: rmse 4.176063e-05  max resid 9.776506e-05 
    ## ... Similar to previous best
    ## Run 57 stress 0.07732932 
    ## Run 58 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744416  max resid 0.05105542 
    ## Run 59 stress 0.08288055 
    ## Run 60 stress 0.08000467 
    ## Run 61 stress 0.0762924 
    ## Run 62 stress 0.08000475 
    ## Run 63 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744606  max resid 0.05104413 
    ## Run 64 stress 0.0737823 
    ## ... Procrustes: rmse 0.01746763  max resid 0.0511602 
    ## Run 65 stress 0.3356087 
    ## Run 66 stress 0.07629238 
    ## Run 67 stress 0.07629234 
    ## Run 68 stress 0.0823376 
    ## Run 69 stress 0.07365783 
    ## ... Procrustes: rmse 1.04423e-05  max resid 2.406068e-05 
    ## ... Similar to previous best
    ## Run 70 stress 0.07629237 
    ## Run 71 stress 0.07732926 
    ## Run 72 stress 0.08233757 
    ## Run 73 stress 0.07629236 
    ## Run 74 stress 0.08288052 
    ## Run 75 stress 0.07365784 
    ## ... Procrustes: rmse 3.891564e-05  max resid 9.313582e-05 
    ## ... Similar to previous best
    ## Run 76 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744705  max resid 0.05104691 
    ## Run 77 stress 0.0762924 
    ## Run 78 stress 0.0773293 
    ## Run 79 stress 0.07365785 
    ## ... Procrustes: rmse 5.929359e-05  max resid 0.0001401906 
    ## ... Similar to previous best
    ## Run 80 stress 0.0800047 
    ## Run 81 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001085708  max resid 0.0002567269 
    ## ... Similar to previous best
    ## Run 82 stress 0.0828804 
    ## Run 83 stress 0.07629235 
    ## Run 84 stress 0.07629234 
    ## Run 85 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745272  max resid 0.05109046 
    ## Run 86 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743321  max resid 0.05099538 
    ## Run 87 stress 0.08288033 
    ## Run 88 stress 0.08000474 
    ## Run 89 stress 0.07732934 
    ## Run 90 stress 0.08000467 
    ## Run 91 stress 0.07365788 
    ## ... Procrustes: rmse 0.000113432  max resid 0.0002694059 
    ## ... Similar to previous best
    ## Run 92 stress 0.07629239 
    ## Run 93 stress 0.08000474 
    ## Run 94 stress 0.0828805 
    ## Run 95 stress 0.07629236 
    ## Run 96 stress 0.08000468 
    ## Run 97 stress 0.07365785 
    ## ... Procrustes: rmse 6.285587e-05  max resid 0.0001466517 
    ## ... Similar to previous best
    ## Run 98 stress 0.07629239 
    ## Run 99 stress 0.07629233 
    ## Run 100 stress 0.07629238 
    ## Run 101 stress 0.07365783 
    ## ... Procrustes: rmse 2.821636e-05  max resid 6.644087e-05 
    ## ... Similar to previous best
    ## Run 102 stress 0.07629233 
    ## Run 103 stress 0.07629237 
    ## Run 104 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744345  max resid 0.0510576 
    ## Run 105 stress 0.07365784 
    ## ... Procrustes: rmse 2.919899e-05  max resid 6.309968e-05 
    ## ... Similar to previous best
    ## Run 106 stress 0.07629234 
    ## Run 107 stress 0.08000467 
    ## Run 108 stress 0.07629241 
    ## Run 109 stress 0.07365783 
    ## ... Procrustes: rmse 2.359721e-05  max resid 5.129847e-05 
    ## ... Similar to previous best
    ## Run 110 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001386003  max resid 0.0003277707 
    ## ... Similar to previous best
    ## Run 111 stress 0.0737823 
    ## ... Procrustes: rmse 0.0174589  max resid 0.0510742 
    ## Run 112 stress 0.07365783 
    ## ... Procrustes: rmse 1.037391e-05  max resid 2.365405e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.07732932 
    ## Run 114 stress 0.08000468 
    ## Run 115 stress 0.08233759 
    ## Run 116 stress 0.0800047 
    ## Run 117 stress 0.0762924 
    ## Run 118 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744378  max resid 0.05103763 
    ## Run 119 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745751  max resid 0.05109423 
    ## Run 120 stress 0.08000468 
    ## Run 121 stress 0.07629234 
    ## Run 122 stress 0.07629234 
    ## Run 123 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001106357  max resid 0.0002517416 
    ## ... Similar to previous best
    ## Run 124 stress 0.07629237 
    ## Run 125 stress 0.07365783 
    ## ... Procrustes: rmse 7.065412e-06  max resid 1.728443e-05 
    ## ... Similar to previous best
    ## Run 126 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744551  max resid 0.05106913 
    ## Run 127 stress 0.07365786 
    ## ... Procrustes: rmse 2.462614e-05  max resid 4.112016e-05 
    ## ... Similar to previous best
    ## Run 128 stress 0.07732927 
    ## Run 129 stress 0.08288048 
    ## Run 130 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744287  max resid 0.05102909 
    ## Run 131 stress 0.07365789 
    ## ... Procrustes: rmse 6.2963e-05  max resid 0.0001224705 
    ## ... Similar to previous best
    ## Run 132 stress 0.08233755 
    ## Run 133 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744929  max resid 0.05106627 
    ## Run 134 stress 0.07365784 
    ## ... Procrustes: rmse 4.113324e-05  max resid 9.730508e-05 
    ## ... Similar to previous best
    ## Run 135 stress 0.07629238 
    ## Run 136 stress 0.07629234 
    ## Run 137 stress 0.07365784 
    ## ... Procrustes: rmse 5.217842e-05  max resid 0.0001188794 
    ## ... Similar to previous best
    ## Run 138 stress 0.08000469 
    ## Run 139 stress 0.08288062 
    ## Run 140 stress 0.08000468 
    ## Run 141 stress 0.07629244 
    ## Run 142 stress 0.07365786 
    ## ... Procrustes: rmse 8.967558e-05  max resid 0.000209899 
    ## ... Similar to previous best
    ## Run 143 stress 0.07365783 
    ## ... Procrustes: rmse 2.412312e-05  max resid 5.494935e-05 
    ## ... Similar to previous best
    ## Run 144 stress 0.07629243 
    ## Run 145 stress 0.07629235 
    ## Run 146 stress 0.08233768 
    ## Run 147 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 8.694906e-06  max resid 1.910868e-05 
    ## ... Similar to previous best
    ## Run 148 stress 0.07365784 
    ## ... Procrustes: rmse 4.281764e-05  max resid 0.0001011682 
    ## ... Similar to previous best
    ## Run 149 stress 0.07365783 
    ## ... Procrustes: rmse 1.323881e-05  max resid 2.881838e-05 
    ## ... Similar to previous best
    ## Run 150 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744564  max resid 0.05103418 
    ## Run 151 stress 0.08288041 
    ## Run 152 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744091  max resid 0.05105289 
    ## Run 153 stress 0.07629233 
    ## Run 154 stress 0.0736579 
    ## ... Procrustes: rmse 0.000135571  max resid 0.0003139011 
    ## ... Similar to previous best
    ## Run 155 stress 0.07629255 
    ## Run 156 stress 0.07629236 
    ## Run 157 stress 0.07629243 
    ## Run 158 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745307  max resid 0.05105342 
    ## Run 159 stress 0.08288044 
    ## Run 160 stress 0.07378233 
    ## ... Procrustes: rmse 0.0174399  max resid 0.05105123 
    ## Run 161 stress 0.08000468 
    ## Run 162 stress 0.07732928 
    ## Run 163 stress 0.07365783 
    ## ... Procrustes: rmse 1.68702e-05  max resid 3.609673e-05 
    ## ... Similar to previous best
    ## Run 164 stress 0.08233757 
    ## Run 165 stress 0.07732936 
    ## Run 166 stress 0.07365785 
    ## ... Procrustes: rmse 3.627782e-05  max resid 8.25421e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744734  max resid 0.05104889 
    ## Run 168 stress 0.07629254 
    ## Run 169 stress 0.07629235 
    ## Run 170 stress 0.08288035 
    ## Run 171 stress 0.08233755 
    ## Run 172 stress 0.07732942 
    ## Run 173 stress 0.07629239 
    ## Run 174 stress 0.07629238 
    ## Run 175 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001228401  max resid 0.0002900484 
    ## ... Similar to previous best
    ## Run 176 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001016496  max resid 0.0002352799 
    ## ... Similar to previous best
    ## Run 177 stress 0.08000478 
    ## Run 178 stress 0.07365784 
    ## ... Procrustes: rmse 5.137009e-05  max resid 0.000119303 
    ## ... Similar to previous best
    ## Run 179 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001080527  max resid 0.0002578195 
    ## ... Similar to previous best
    ## Run 180 stress 0.07365785 
    ## ... Procrustes: rmse 7.583877e-05  max resid 0.0001798537 
    ## ... Similar to previous best
    ## Run 181 stress 0.07629238 
    ## Run 182 stress 0.3584397 
    ## Run 183 stress 0.08233765 
    ## Run 184 stress 0.07629233 
    ## Run 185 stress 0.07365783 
    ## ... Procrustes: rmse 3.039144e-05  max resid 6.958866e-05 
    ## ... Similar to previous best
    ## Run 186 stress 0.07732931 
    ## Run 187 stress 0.0800047 
    ## Run 188 stress 0.07629237 
    ## Run 189 stress 0.07365785 
    ## ... Procrustes: rmse 6.927946e-05  max resid 0.0001635821 
    ## ... Similar to previous best
    ## Run 190 stress 0.07629241 
    ## Run 191 stress 0.07365788 
    ## ... Procrustes: rmse 0.000116844  max resid 0.000270757 
    ## ... Similar to previous best
    ## Run 192 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001115391  max resid 0.0002612115 
    ## ... Similar to previous best
    ## Run 193 stress 0.0800047 
    ## Run 194 stress 0.08000468 
    ## Run 195 stress 0.07629233 
    ## Run 196 stress 0.08000471 
    ## Run 197 stress 0.07629233 
    ## Run 198 stress 0.08233761 
    ## Run 199 stress 0.08000467 
    ## Run 200 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174441  max resid 0.05101975 
    ## Run 201 stress 0.07365783 
    ## ... Procrustes: rmse 2.596616e-05  max resid 6.158473e-05 
    ## ... Similar to previous best
    ## Run 202 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745739  max resid 0.05108006 
    ## Run 203 stress 0.08000467 
    ## Run 204 stress 0.07629234 
    ## Run 205 stress 0.07629239 
    ## Run 206 stress 0.0762924 
    ## Run 207 stress 0.07732949 
    ## Run 208 stress 0.07732939 
    ## Run 209 stress 0.07365786 
    ## ... Procrustes: rmse 9.002124e-05  max resid 0.0002119162 
    ## ... Similar to previous best
    ## Run 210 stress 0.07365783 
    ## ... Procrustes: rmse 7.310039e-06  max resid 1.550988e-05 
    ## ... Similar to previous best
    ## Run 211 stress 0.07732938 
    ## Run 212 stress 0.07629234 
    ## Run 213 stress 0.08233753 
    ## Run 214 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745374  max resid 0.05105063 
    ## Run 215 stress 0.07365783 
    ## ... Procrustes: rmse 2.036399e-05  max resid 4.908081e-05 
    ## ... Similar to previous best
    ## Run 216 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744658  max resid 0.05105896 
    ## Run 217 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744612  max resid 0.05103849 
    ## Run 218 stress 0.07378226 
    ## ... Procrustes: rmse 0.017448  max resid 0.05105586 
    ## Run 219 stress 0.08000467 
    ## Run 220 stress 0.07378232 
    ## ... Procrustes: rmse 0.01746167  max resid 0.05113824 
    ## Run 221 stress 0.07629234 
    ## Run 222 stress 0.07365783 
    ## ... Procrustes: rmse 1.247274e-05  max resid 2.055165e-05 
    ## ... Similar to previous best
    ## Run 223 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744413  max resid 0.05103229 
    ## Run 224 stress 0.07732948 
    ## Run 225 stress 0.07365783 
    ## ... Procrustes: rmse 1.192792e-05  max resid 2.759517e-05 
    ## ... Similar to previous best
    ## Run 226 stress 0.08000476 
    ## Run 227 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001296395  max resid 0.0003071113 
    ## ... Similar to previous best
    ## Run 228 stress 0.08000468 
    ## Run 229 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744689  max resid 0.05103811 
    ## Run 230 stress 0.07629238 
    ## Run 231 stress 0.07629235 
    ## Run 232 stress 0.07629242 
    ## Run 233 stress 0.07732944 
    ## Run 234 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744904  max resid 0.05106123 
    ## Run 235 stress 0.07629243 
    ## Run 236 stress 0.07629239 
    ## Run 237 stress 0.08000475 
    ## Run 238 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744642  max resid 0.05106134 
    ## Run 239 stress 0.07365784 
    ## ... Procrustes: rmse 5.360135e-05  max resid 0.0001261246 
    ## ... Similar to previous best
    ## Run 240 stress 0.07629238 
    ## Run 241 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744115  max resid 0.05102855 
    ## Run 242 stress 0.07365784 
    ## ... Procrustes: rmse 3.759218e-05  max resid 8.893661e-05 
    ## ... Similar to previous best
    ## Run 243 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744588  max resid 0.05102419 
    ## Run 244 stress 0.07365783 
    ## ... Procrustes: rmse 2.600746e-05  max resid 6.216258e-05 
    ## ... Similar to previous best
    ## Run 245 stress 0.07629242 
    ## Run 246 stress 0.08000474 
    ## Run 247 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745812  max resid 0.05109985 
    ## Run 248 stress 0.07629248 
    ## Run 249 stress 0.07365784 
    ## ... Procrustes: rmse 5.83389e-05  max resid 0.0001357837 
    ## ... Similar to previous best
    ## Run 250 stress 0.08000468 
    ## Run 251 stress 0.07732938 
    ## Run 252 stress 0.08288044 
    ## Run 253 stress 0.07732925 
    ## Run 254 stress 0.07365783 
    ## ... Procrustes: rmse 1.35298e-05  max resid 3.208114e-05 
    ## ... Similar to previous best
    ## Run 255 stress 0.08288046 
    ## Run 256 stress 0.07629235 
    ## Run 257 stress 0.07629234 
    ## Run 258 stress 0.08233761 
    ## Run 259 stress 0.07365783 
    ## ... Procrustes: rmse 2.597037e-05  max resid 6.010332e-05 
    ## ... Similar to previous best
    ## Run 260 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744397  max resid 0.051057 
    ## Run 261 stress 0.0762924 
    ## Run 262 stress 0.08000469 
    ## Run 263 stress 0.08288046 
    ## Run 264 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744597  max resid 0.05103794 
    ## Run 265 stress 0.08233762 
    ## Run 266 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744525  max resid 0.05102043 
    ## Run 267 stress 0.07629241 
    ## Run 268 stress 0.07732939 
    ## Run 269 stress 0.07365785 
    ## ... Procrustes: rmse 2.121874e-05  max resid 4.244829e-05 
    ## ... Similar to previous best
    ## Run 270 stress 0.07365784 
    ## ... Procrustes: rmse 4.310873e-05  max resid 0.0001011715 
    ## ... Similar to previous best
    ## Run 271 stress 0.08000468 
    ## Run 272 stress 0.07365784 
    ## ... Procrustes: rmse 5.460347e-05  max resid 0.000129291 
    ## ... Similar to previous best
    ## Run 273 stress 0.08000468 
    ## Run 274 stress 0.08233757 
    ## Run 275 stress 0.08000468 
    ## Run 276 stress 0.3578478 
    ## Run 277 stress 0.07365787 
    ## ... Procrustes: rmse 9.551768e-05  max resid 0.0002267133 
    ## ... Similar to previous best
    ## Run 278 stress 0.07629234 
    ## Run 279 stress 0.07629235 
    ## Run 280 stress 0.08000468 
    ## Run 281 stress 0.07629233 
    ## Run 282 stress 0.07365783 
    ## ... Procrustes: rmse 8.739884e-06  max resid 2.059875e-05 
    ## ... Similar to previous best
    ## Run 283 stress 0.08000468 
    ## Run 284 stress 0.07365784 
    ## ... Procrustes: rmse 5.986316e-05  max resid 0.0001424526 
    ## ... Similar to previous best
    ## Run 285 stress 0.07629234 
    ## Run 286 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744658  max resid 0.051043 
    ## Run 287 stress 0.07629248 
    ## Run 288 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001090312  max resid 0.0002539208 
    ## ... Similar to previous best
    ## Run 289 stress 0.07629237 
    ## Run 290 stress 0.08000469 
    ## Run 291 stress 0.07629239 
    ## Run 292 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745358  max resid 0.05109787 
    ## Run 293 stress 0.07629236 
    ## Run 294 stress 0.08233759 
    ## Run 295 stress 0.0828806 
    ## Run 296 stress 0.07378228 
    ## ... Procrustes: rmse 0.01746356  max resid 0.05113435 
    ## Run 297 stress 0.07629236 
    ## Run 298 stress 0.07629234 
    ## Run 299 stress 0.0762924 
    ## Run 300 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744841  max resid 0.05105782 
    ## Run 301 stress 0.08233762 
    ## Run 302 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174425  max resid 0.05102052 
    ## Run 303 stress 0.07732931 
    ## Run 304 stress 0.07629235 
    ## Run 305 stress 0.07365786 
    ## ... Procrustes: rmse 6.493955e-05  max resid 0.0001440782 
    ## ... Similar to previous best
    ## Run 306 stress 0.07365784 
    ## ... Procrustes: rmse 4.304986e-05  max resid 0.000102566 
    ## ... Similar to previous best
    ## Run 307 stress 0.08000476 
    ## Run 308 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744844  max resid 0.05105838 
    ## Run 309 stress 0.07629233 
    ## Run 310 stress 0.08233755 
    ## Run 311 stress 0.07732931 
    ## Run 312 stress 0.07629238 
    ## Run 313 stress 0.07629233 
    ## Run 314 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174468  max resid 0.05105401 
    ## Run 315 stress 0.07629235 
    ## Run 316 stress 0.08288053 
    ## Run 317 stress 0.07365784 
    ## ... Procrustes: rmse 5.986793e-05  max resid 0.000142174 
    ## ... Similar to previous best
    ## Run 318 stress 0.08000473 
    ## Run 319 stress 0.0800047 
    ## Run 320 stress 0.07365783 
    ## ... Procrustes: rmse 2.157704e-05  max resid 5.042822e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.07629233 
    ## Run 322 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744917  max resid 0.05106037 
    ## Run 323 stress 0.08000472 
    ## Run 324 stress 0.3232345 
    ## Run 325 stress 0.08233751 
    ## Run 326 stress 0.07629233 
    ## Run 327 stress 0.07365786 
    ## ... Procrustes: rmse 8.879233e-05  max resid 0.0002086319 
    ## ... Similar to previous best
    ## Run 328 stress 0.07629238 
    ## Run 329 stress 0.07365783 
    ## ... Procrustes: rmse 9.54179e-06  max resid 2.259454e-05 
    ## ... Similar to previous best
    ## Run 330 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744324  max resid 0.05100933 
    ## Run 331 stress 0.07365784 
    ## ... Procrustes: rmse 2.207421e-05  max resid 4.469527e-05 
    ## ... Similar to previous best
    ## Run 332 stress 0.07365784 
    ## ... Procrustes: rmse 4.561863e-05  max resid 0.0001065622 
    ## ... Similar to previous best
    ## Run 333 stress 0.07629251 
    ## Run 334 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745422  max resid 0.05108938 
    ## Run 335 stress 0.08000469 
    ## Run 336 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744523  max resid 0.05103191 
    ## Run 337 stress 0.07365783 
    ## ... Procrustes: rmse 9.086694e-06  max resid 2.132169e-05 
    ## ... Similar to previous best
    ## Run 338 stress 0.08000469 
    ## Run 339 stress 0.07365783 
    ## ... Procrustes: rmse 3.397779e-05  max resid 7.926632e-05 
    ## ... Similar to previous best
    ## Run 340 stress 0.07365784 
    ## ... Procrustes: rmse 5.629177e-05  max resid 0.0001330473 
    ## ... Similar to previous best
    ## Run 341 stress 0.07365787 
    ## ... Procrustes: rmse 9.806776e-05  max resid 0.0002277519 
    ## ... Similar to previous best
    ## Run 342 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744333  max resid 0.05101428 
    ## Run 343 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001146742  max resid 0.0002701022 
    ## ... Similar to previous best
    ## Run 344 stress 0.07365789 
    ## ... Procrustes: rmse 0.000125367  max resid 0.0002899735 
    ## ... Similar to previous best
    ## Run 345 stress 0.07732934 
    ## Run 346 stress 0.07378229 
    ## ... Procrustes: rmse 0.0174461  max resid 0.05103295 
    ## Run 347 stress 0.08000467 
    ## Run 348 stress 0.07629248 
    ## Run 349 stress 0.08000469 
    ## Run 350 stress 0.07365783 
    ## ... Procrustes: rmse 3.561786e-05  max resid 8.386948e-05 
    ## ... Similar to previous best
    ## Run 351 stress 0.07365785 
    ## ... Procrustes: rmse 7.202333e-05  max resid 0.0001695382 
    ## ... Similar to previous best
    ## Run 352 stress 0.07629248 
    ## Run 353 stress 0.08233774 
    ## Run 354 stress 0.07365786 
    ## ... Procrustes: rmse 8.065131e-05  max resid 0.0001905321 
    ## ... Similar to previous best
    ## Run 355 stress 0.08288054 
    ## Run 356 stress 0.08288039 
    ## Run 357 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001344476  max resid 0.0003144677 
    ## ... Similar to previous best
    ## Run 358 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174473  max resid 0.05104351 
    ## Run 359 stress 0.07365784 
    ## ... Procrustes: rmse 4.968067e-05  max resid 0.0001177099 
    ## ... Similar to previous best
    ## Run 360 stress 0.08000472 
    ## Run 361 stress 0.0736579 
    ## ... Procrustes: rmse 0.000125066  max resid 0.0002975083 
    ## ... Similar to previous best
    ## Run 362 stress 0.08233753 
    ## Run 363 stress 0.08000469 
    ## Run 364 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001348893  max resid 0.0003176285 
    ## ... Similar to previous best
    ## Run 365 stress 0.07365784 
    ## ... Procrustes: rmse 4.751027e-05  max resid 0.0001122867 
    ## ... Similar to previous best
    ## Run 366 stress 0.08000468 
    ## Run 367 stress 0.07365784 
    ## ... Procrustes: rmse 4.938345e-05  max resid 0.0001178862 
    ## ... Similar to previous best
    ## Run 368 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744561  max resid 0.05106223 
    ## Run 369 stress 0.07365784 
    ## ... Procrustes: rmse 5.289085e-05  max resid 0.0001260152 
    ## ... Similar to previous best
    ## Run 370 stress 0.07365783 
    ## ... Procrustes: rmse 2.213408e-05  max resid 5.17376e-05 
    ## ... Similar to previous best
    ## Run 371 stress 0.07629245 
    ## Run 372 stress 0.07365783 
    ## ... Procrustes: rmse 2.483172e-05  max resid 5.913018e-05 
    ## ... Similar to previous best
    ## Run 373 stress 0.07365784 
    ## ... Procrustes: rmse 5.900804e-05  max resid 0.0001386865 
    ## ... Similar to previous best
    ## Run 374 stress 0.08000468 
    ## Run 375 stress 0.08000467 
    ## Run 376 stress 0.08233756 
    ## Run 377 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001016625  max resid 0.0002419685 
    ## ... Similar to previous best
    ## Run 378 stress 0.07365783 
    ## ... Procrustes: rmse 3.026164e-05  max resid 7.203424e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.07629249 
    ## Run 380 stress 0.0773295 
    ## Run 381 stress 0.08000469 
    ## Run 382 stress 0.0762924 
    ## Run 383 stress 0.07629234 
    ## Run 384 stress 0.0800047 
    ## Run 385 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744786  max resid 0.0510583 
    ## Run 386 stress 0.07365783 
    ## ... Procrustes: rmse 6.315266e-06  max resid 1.296957e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174483  max resid 0.05105542 
    ## Run 388 stress 0.07732941 
    ## Run 389 stress 0.0823376 
    ## Run 390 stress 0.08000468 
    ## Run 391 stress 0.07365785 
    ## ... Procrustes: rmse 5.503223e-05  max resid 0.0001201608 
    ## ... Similar to previous best
    ## Run 392 stress 0.07365783 
    ## ... Procrustes: rmse 1.474456e-05  max resid 3.489062e-05 
    ## ... Similar to previous best
    ## Run 393 stress 0.07365784 
    ## ... Procrustes: rmse 2.981279e-05  max resid 6.668888e-05 
    ## ... Similar to previous best
    ## Run 394 stress 0.08000478 
    ## Run 395 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744905  max resid 0.05106917 
    ## Run 396 stress 0.2829247 
    ## Run 397 stress 0.07629238 
    ## Run 398 stress 0.07365784 
    ## ... Procrustes: rmse 4.441697e-05  max resid 0.0001063559 
    ## ... Similar to previous best
    ## Run 399 stress 0.08000467 
    ## Run 400 stress 0.07365786 
    ## ... Procrustes: rmse 5.69042e-05  max resid 0.0001344992 
    ## ... Similar to previous best
    ## Run 401 stress 0.07629239 
    ## Run 402 stress 0.07629237 
    ## Run 403 stress 0.08000468 
    ## Run 404 stress 0.07365784 
    ## ... Procrustes: rmse 2.902872e-05  max resid 6.900471e-05 
    ## ... Similar to previous best
    ## Run 405 stress 0.08000467 
    ## Run 406 stress 0.07365784 
    ## ... Procrustes: rmse 4.930891e-05  max resid 0.0001170902 
    ## ... Similar to previous best
    ## Run 407 stress 0.08000469 
    ## Run 408 stress 0.08000468 
    ## Run 409 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744556  max resid 0.05103196 
    ## Run 410 stress 0.07629247 
    ## Run 411 stress 0.07378229 
    ## ... Procrustes: rmse 0.01743708  max resid 0.05097964 
    ## Run 412 stress 0.07629234 
    ## Run 413 stress 0.0828806 
    ## Run 414 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744555  max resid 0.05106757 
    ## Run 415 stress 0.07365784 
    ## ... Procrustes: rmse 3.649903e-05  max resid 8.638508e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.07365784 
    ## ... Procrustes: rmse 5.374841e-05  max resid 0.0001281702 
    ## ... Similar to previous best
    ## Run 417 stress 0.07365784 
    ## ... Procrustes: rmse 3.755967e-05  max resid 8.802415e-05 
    ## ... Similar to previous best
    ## Run 418 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745892  max resid 0.0510809 
    ## Run 419 stress 0.07732935 
    ## Run 420 stress 0.07365785 
    ## ... Procrustes: rmse 7.357579e-05  max resid 0.0001717678 
    ## ... Similar to previous best
    ## Run 421 stress 0.07629253 
    ## Run 422 stress 0.07365784 
    ## ... Procrustes: rmse 4.79463e-05  max resid 0.0001109395 
    ## ... Similar to previous best
    ## Run 423 stress 0.0800047 
    ## Run 424 stress 0.07378233 
    ## ... Procrustes: rmse 0.01744205  max resid 0.05099109 
    ## Run 425 stress 0.07365796 
    ## ... Procrustes: rmse 5.378094e-05  max resid 9.832069e-05 
    ## ... Similar to previous best
    ## Run 426 stress 0.07629233 
    ## Run 427 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744403  max resid 0.05101979 
    ## Run 428 stress 0.07629237 
    ## Run 429 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744628  max resid 0.05103988 
    ## Run 430 stress 0.07365785 
    ## ... Procrustes: rmse 8.366952e-05  max resid 0.0001979666 
    ## ... Similar to previous best
    ## Run 431 stress 0.08233752 
    ## Run 432 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745189  max resid 0.05106628 
    ## Run 433 stress 0.07629237 
    ## Run 434 stress 0.08288054 
    ## Run 435 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001456994  max resid 0.000343912 
    ## ... Similar to previous best
    ## Run 436 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745235  max resid 0.05109255 
    ## Run 437 stress 0.08000468 
    ## Run 438 stress 0.08000472 
    ## Run 439 stress 0.0762924 
    ## Run 440 stress 0.07629238 
    ## Run 441 stress 0.07365785 
    ## ... Procrustes: rmse 4.486902e-05  max resid 0.0001058479 
    ## ... Similar to previous best
    ## Run 442 stress 0.08000473 
    ## Run 443 stress 0.08233762 
    ## Run 444 stress 0.07365783 
    ## ... Procrustes: rmse 1.030003e-05  max resid 1.649573e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.07365783 
    ## ... Procrustes: rmse 2.387201e-05  max resid 5.729707e-05 
    ## ... Similar to previous best
    ## Run 446 stress 0.08000467 
    ## Run 447 stress 0.07732928 
    ## Run 448 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744268  max resid 0.05105793 
    ## Run 449 stress 0.07629238 
    ## Run 450 stress 0.07365786 
    ## ... Procrustes: rmse 5.351693e-05  max resid 0.0001268396 
    ## ... Similar to previous best
    ## Run 451 stress 0.07629238 
    ## Run 452 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744489  max resid 0.05102756 
    ## Run 453 stress 0.07365783 
    ## ... Procrustes: rmse 3.385859e-05  max resid 7.851956e-05 
    ## ... Similar to previous best
    ## Run 454 stress 0.07365785 
    ## ... Procrustes: rmse 5.561491e-05  max resid 0.0001288708 
    ## ... Similar to previous best
    ## Run 455 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001140659  max resid 0.0002676517 
    ## ... Similar to previous best
    ## Run 456 stress 0.07365783 
    ## ... Procrustes: rmse 2.030677e-05  max resid 3.610015e-05 
    ## ... Similar to previous best
    ## Run 457 stress 0.08000467 
    ## Run 458 stress 0.07365785 
    ## ... Procrustes: rmse 7.24686e-05  max resid 0.0001703192 
    ## ... Similar to previous best
    ## Run 459 stress 0.2505361 
    ## Run 460 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001487427  max resid 0.000352347 
    ## ... Similar to previous best
    ## Run 461 stress 0.07629235 
    ## Run 462 stress 0.08000467 
    ## Run 463 stress 0.08000478 
    ## Run 464 stress 0.07365783 
    ## ... Procrustes: rmse 1.441744e-05  max resid 3.393038e-05 
    ## ... Similar to previous best
    ## Run 465 stress 0.0762924 
    ## Run 466 stress 0.07365783 
    ## ... Procrustes: rmse 2.1885e-05  max resid 5.124711e-05 
    ## ... Similar to previous best
    ## Run 467 stress 0.07732945 
    ## Run 468 stress 0.07629233 
    ## Run 469 stress 0.08000472 
    ## Run 470 stress 0.07365783 
    ## ... Procrustes: rmse 5.846581e-06  max resid 1.376187e-05 
    ## ... Similar to previous best
    ## Run 471 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001376172  max resid 0.0003284067 
    ## ... Similar to previous best
    ## Run 472 stress 0.07365784 
    ## ... Procrustes: rmse 2.860341e-05  max resid 6.727608e-05 
    ## ... Similar to previous best
    ## Run 473 stress 0.0762924 
    ## Run 474 stress 0.07365785 
    ## ... Procrustes: rmse 6.949225e-05  max resid 0.0001634391 
    ## ... Similar to previous best
    ## Run 475 stress 0.3495714 
    ## Run 476 stress 0.07365784 
    ## ... Procrustes: rmse 4.729113e-05  max resid 0.0001117711 
    ## ... Similar to previous best
    ## Run 477 stress 0.08000469 
    ## Run 478 stress 0.07365784 
    ## ... Procrustes: rmse 6.178546e-05  max resid 0.0001441269 
    ## ... Similar to previous best
    ## Run 479 stress 0.08233754 
    ## Run 480 stress 0.07365783 
    ## ... Procrustes: rmse 2.168244e-05  max resid 4.5359e-05 
    ## ... Similar to previous best
    ## Run 481 stress 0.08000474 
    ## Run 482 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744665  max resid 0.05106431 
    ## Run 483 stress 0.07365786 
    ## ... Procrustes: rmse 8.509218e-05  max resid 0.0002029351 
    ## ... Similar to previous best
    ## Run 484 stress 0.08000468 
    ## Run 485 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744544  max resid 0.05105864 
    ## Run 486 stress 0.07629238 
    ## Run 487 stress 0.07732927 
    ## Run 488 stress 0.0762924 
    ## Run 489 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744712  max resid 0.05104566 
    ## Run 490 stress 0.07629237 
    ## Run 491 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745333  max resid 0.0510813 
    ## Run 492 stress 0.3495562 
    ## Run 493 stress 0.07629244 
    ## Run 494 stress 0.07365783 
    ## ... Procrustes: rmse 1.520263e-05  max resid 3.61126e-05 
    ## ... Similar to previous best
    ## Run 495 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744138  max resid 0.05102137 
    ## Run 496 stress 0.07629233 
    ## Run 497 stress 0.07629241 
    ## Run 498 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743133  max resid 0.05094625 
    ## Run 499 stress 0.07365783 
    ## ... Procrustes: rmse 1.182244e-05  max resid 2.094783e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744493  max resid 0.05104598 
    ## *** Best solution repeated 102 times

``` r
# Stratified lakes and ocean sites
SD_beta_SO_NMDS <- metaMDS(SD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.08340314 
    ## Run 2 stress 0.08448435 
    ## Run 3 stress 0.069782 
    ## ... New best solution
    ## ... Procrustes: rmse 0.09984689  max resid 0.2618872 
    ## Run 4 stress 0.07844947 
    ## Run 5 stress 0.06978191 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002416744  max resid 0.0006345176 
    ## ... Similar to previous best
    ## Run 6 stress 0.06942777 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01323484  max resid 0.03353099 
    ## Run 7 stress 0.06942777 
    ## ... Procrustes: rmse 1.057661e-05  max resid 2.73964e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.08448437 
    ## Run 9 stress 0.07844936 
    ## Run 10 stress 0.08340293 
    ## Run 11 stress 0.07428313 
    ## Run 12 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326256  max resid 0.03334686 
    ## Run 13 stress 0.08448438 
    ## Run 14 stress 0.07970526 
    ## Run 15 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327296  max resid 0.03336761 
    ## Run 16 stress 0.08340288 
    ## Run 17 stress 0.07970526 
    ## Run 18 stress 0.0742832 
    ## Run 19 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 1.014477e-05  max resid 2.811755e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.08448455 
    ## Run 21 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.079334e-05  max resid 5.34837e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.07250812 
    ## Run 23 stress 0.07250812 
    ## Run 24 stress 0.07428313 
    ## Run 25 stress 0.08340287 
    ## Run 26 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319762  max resid 0.03321046 
    ## Run 27 stress 0.07428316 
    ## Run 28 stress 0.07428315 
    ## Run 29 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 4.632387e-06  max resid 1.254352e-05 
    ## ... Similar to previous best
    ## Run 30 stress 0.08448435 
    ## Run 31 stress 0.07970525 
    ## Run 32 stress 0.0742832 
    ## Run 33 stress 0.08340287 
    ## Run 34 stress 0.0844844 
    ## Run 35 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325238  max resid 0.03332426 
    ## Run 36 stress 0.07250814 
    ## Run 37 stress 0.06942777 
    ## ... Procrustes: rmse 3.180464e-05  max resid 8.349226e-05 
    ## ... Similar to previous best
    ## Run 38 stress 0.07250812 
    ## Run 39 stress 0.06978192 
    ## ... Procrustes: rmse 0.013251  max resid 0.03332039 
    ## Run 40 stress 0.07250815 
    ## Run 41 stress 0.07250812 
    ## Run 42 stress 0.07428313 
    ## Run 43 stress 0.07250812 
    ## Run 44 stress 0.08340286 
    ## Run 45 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322795  max resid 0.0332708 
    ## Run 46 stress 0.06942776 
    ## ... Procrustes: rmse 2.298417e-05  max resid 5.925565e-05 
    ## ... Similar to previous best
    ## Run 47 stress 0.07970525 
    ## Run 48 stress 0.0742832 
    ## Run 49 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319902  max resid 0.03321134 
    ## Run 50 stress 0.07428313 
    ## Run 51 stress 0.06942777 
    ## ... Procrustes: rmse 5.266239e-05  max resid 0.000136942 
    ## ... Similar to previous best
    ## Run 52 stress 0.07250812 
    ## Run 53 stress 0.07428313 
    ## Run 54 stress 0.08340286 
    ## Run 55 stress 0.07428318 
    ## Run 56 stress 0.08340288 
    ## Run 57 stress 0.07970525 
    ## Run 58 stress 0.06978195 
    ## ... Procrustes: rmse 0.01317716  max resid 0.03316378 
    ## Run 59 stress 0.07970527 
    ## Run 60 stress 0.06942776 
    ## ... Procrustes: rmse 2.322742e-05  max resid 6.151698e-05 
    ## ... Similar to previous best
    ## Run 61 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327269  max resid 0.03336911 
    ## Run 62 stress 0.0834029 
    ## Run 63 stress 0.07250815 
    ## Run 64 stress 0.08340287 
    ## Run 65 stress 0.06978199 
    ## ... Procrustes: rmse 0.0132786  max resid 0.03337454 
    ## Run 66 stress 0.08340295 
    ## Run 67 stress 0.0742832 
    ## Run 68 stress 0.06942776 
    ## ... Procrustes: rmse 6.612964e-06  max resid 1.777931e-05 
    ## ... Similar to previous best
    ## Run 69 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326287  max resid 0.03334633 
    ## Run 70 stress 0.08340289 
    ## Run 71 stress 0.08340288 
    ## Run 72 stress 0.07970525 
    ## Run 73 stress 0.06942776 
    ## ... Procrustes: rmse 2.26498e-06  max resid 6.534731e-06 
    ## ... Similar to previous best
    ## Run 74 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326961  max resid 0.03335527 
    ## Run 75 stress 0.07428315 
    ## Run 76 stress 0.08340289 
    ## Run 77 stress 0.07250813 
    ## Run 78 stress 0.07428312 
    ## Run 79 stress 0.07428314 
    ## Run 80 stress 0.08340287 
    ## Run 81 stress 0.07428314 
    ## Run 82 stress 0.07428313 
    ## Run 83 stress 0.08340294 
    ## Run 84 stress 0.06978201 
    ## ... Procrustes: rmse 0.01328574  max resid 0.03339313 
    ## Run 85 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327168  max resid 0.03336252 
    ## Run 86 stress 0.07970526 
    ## Run 87 stress 0.07970525 
    ## Run 88 stress 0.06942777 
    ## ... Procrustes: rmse 3.79443e-05  max resid 9.876074e-05 
    ## ... Similar to previous best
    ## Run 89 stress 0.07428314 
    ## Run 90 stress 0.08448441 
    ## Run 91 stress 0.06942776 
    ## ... Procrustes: rmse 1.429829e-06  max resid 3.574335e-06 
    ## ... Similar to previous best
    ## Run 92 stress 0.08340291 
    ## Run 93 stress 0.07250815 
    ## Run 94 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322265  max resid 0.03327212 
    ## Run 95 stress 0.07844923 
    ## Run 96 stress 0.06942776 
    ## ... Procrustes: rmse 4.767808e-06  max resid 1.314352e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.07428313 
    ## Run 98 stress 0.08340286 
    ## Run 99 stress 0.07428312 
    ## Run 100 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325328  max resid 0.03332688 
    ## Run 101 stress 0.07428314 
    ## Run 102 stress 0.07250812 
    ## Run 103 stress 0.07970526 
    ## Run 104 stress 0.07428315 
    ## Run 105 stress 0.06978196 
    ## ... Procrustes: rmse 0.01317432  max resid 0.03315642 
    ## Run 106 stress 0.07428319 
    ## Run 107 stress 0.06978203 
    ## ... Procrustes: rmse 0.01328918  max resid 0.03339794 
    ## Run 108 stress 0.07970526 
    ## Run 109 stress 0.0844846 
    ## Run 110 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324042  max resid 0.03330092 
    ## Run 111 stress 0.07970525 
    ## Run 112 stress 0.07428313 
    ## Run 113 stress 0.07970528 
    ## Run 114 stress 0.06942776 
    ## ... Procrustes: rmse 4.905116e-06  max resid 8.872174e-06 
    ## ... Similar to previous best
    ## Run 115 stress 0.06942776 
    ## ... Procrustes: rmse 2.918073e-05  max resid 7.63737e-05 
    ## ... Similar to previous best
    ## Run 116 stress 0.07428313 
    ## Run 117 stress 0.07428314 
    ## Run 118 stress 0.07250814 
    ## Run 119 stress 0.08340288 
    ## Run 120 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327899  max resid 0.03337985 
    ## Run 121 stress 0.08340289 
    ## Run 122 stress 0.07250812 
    ## Run 123 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322533  max resid 0.03326776 
    ## Run 124 stress 0.07970525 
    ## Run 125 stress 0.08340289 
    ## Run 126 stress 0.07428313 
    ## Run 127 stress 0.06942776 
    ## ... Procrustes: rmse 1.142238e-05  max resid 3.090194e-05 
    ## ... Similar to previous best
    ## Run 128 stress 0.06978195 
    ## ... Procrustes: rmse 0.01317725  max resid 0.03316501 
    ## Run 129 stress 0.07250815 
    ## Run 130 stress 0.06978192 
    ## ... Procrustes: rmse 0.01320856  max resid 0.03322738 
    ## Run 131 stress 0.07428313 
    ## Run 132 stress 0.06942776 
    ## ... Procrustes: rmse 1.727359e-05  max resid 4.346585e-05 
    ## ... Similar to previous best
    ## Run 133 stress 0.07428313 
    ## Run 134 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326965  max resid 0.03335708 
    ## Run 135 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325336  max resid 0.03332706 
    ## Run 136 stress 0.08340293 
    ## Run 137 stress 0.06942776 
    ## ... Procrustes: rmse 6.603196e-06  max resid 1.78976e-05 
    ## ... Similar to previous best
    ## Run 138 stress 0.07428313 
    ## Run 139 stress 0.07970525 
    ## Run 140 stress 0.07428317 
    ## Run 141 stress 0.07250814 
    ## Run 142 stress 0.07428315 
    ## Run 143 stress 0.07250815 
    ## Run 144 stress 0.07428313 
    ## Run 145 stress 0.07250815 
    ## Run 146 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323721  max resid 0.03329236 
    ## Run 147 stress 0.07428312 
    ## Run 148 stress 0.07844945 
    ## Run 149 stress 0.06978191 
    ## ... Procrustes: rmse 0.0132431  max resid 0.03330384 
    ## Run 150 stress 0.07428313 
    ## Run 151 stress 0.07970525 
    ## Run 152 stress 0.08340289 
    ## Run 153 stress 0.07250812 
    ## Run 154 stress 0.07970526 
    ## Run 155 stress 0.08340289 
    ## Run 156 stress 0.07428313 
    ## Run 157 stress 0.06942776 
    ## ... Procrustes: rmse 2.638675e-06  max resid 5.466618e-06 
    ## ... Similar to previous best
    ## Run 158 stress 0.08340288 
    ## Run 159 stress 0.07428318 
    ## Run 160 stress 0.07428313 
    ## Run 161 stress 0.08340287 
    ## Run 162 stress 0.07250812 
    ## Run 163 stress 0.06942776 
    ## ... Procrustes: rmse 1.77397e-05  max resid 4.422465e-05 
    ## ... Similar to previous best
    ## Run 164 stress 0.07428313 
    ## Run 165 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322924  max resid 0.03327658 
    ## Run 166 stress 0.07970526 
    ## Run 167 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320379  max resid 0.03322317 
    ## Run 168 stress 0.07970525 
    ## Run 169 stress 0.06942776 
    ## ... Procrustes: rmse 5.836321e-06  max resid 1.619465e-05 
    ## ... Similar to previous best
    ## Run 170 stress 0.07428318 
    ## Run 171 stress 0.07428313 
    ## Run 172 stress 0.07250816 
    ## Run 173 stress 0.06942777 
    ## ... Procrustes: rmse 4.577085e-05  max resid 0.0001188892 
    ## ... Similar to previous best
    ## Run 174 stress 0.07428322 
    ## Run 175 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132039  max resid 0.03322172 
    ## Run 176 stress 0.07844947 
    ## Run 177 stress 0.08340291 
    ## Run 178 stress 0.07970525 
    ## Run 179 stress 0.06942776 
    ## ... Procrustes: rmse 1.967584e-05  max resid 5.212263e-05 
    ## ... Similar to previous best
    ## Run 180 stress 0.07970525 
    ## Run 181 stress 0.06942776 
    ## ... Procrustes: rmse 2.106464e-05  max resid 5.706365e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.07428313 
    ## Run 183 stress 0.07844929 
    ## Run 184 stress 0.07970525 
    ## Run 185 stress 0.2671616 
    ## Run 186 stress 0.07428315 
    ## Run 187 stress 0.08340286 
    ## Run 188 stress 0.07970525 
    ## Run 189 stress 0.08340292 
    ## Run 190 stress 0.08340286 
    ## Run 191 stress 0.07844928 
    ## Run 192 stress 0.07250812 
    ## Run 193 stress 0.08340296 
    ## Run 194 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326348  max resid 0.03334564 
    ## Run 195 stress 0.06978191 
    ## ... Procrustes: rmse 0.01322574  max resid 0.03328344 
    ## Run 196 stress 0.07970525 
    ## Run 197 stress 0.07428316 
    ## Run 198 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327075  max resid 0.03335925 
    ## Run 199 stress 0.06942776 
    ## ... Procrustes: rmse 6.935238e-06  max resid 1.909085e-05 
    ## ... Similar to previous best
    ## Run 200 stress 0.08448435 
    ## Run 201 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322374  max resid 0.03326369 
    ## Run 202 stress 0.07250813 
    ## Run 203 stress 0.07970525 
    ## Run 204 stress 0.07970526 
    ## Run 205 stress 0.07250813 
    ## Run 206 stress 0.06942777 
    ## ... Procrustes: rmse 3.070117e-05  max resid 8.081157e-05 
    ## ... Similar to previous best
    ## Run 207 stress 0.08448438 
    ## Run 208 stress 0.08340289 
    ## Run 209 stress 0.06942776 
    ## ... Procrustes: rmse 8.048764e-06  max resid 2.41649e-05 
    ## ... Similar to previous best
    ## Run 210 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323062  max resid 0.03327264 
    ## Run 211 stress 0.07970526 
    ## Run 212 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325767  max resid 0.03333662 
    ## Run 213 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324119  max resid 0.03329942 
    ## Run 214 stress 0.07970526 
    ## Run 215 stress 0.07250812 
    ## Run 216 stress 0.08448436 
    ## Run 217 stress 0.06978201 
    ## ... Procrustes: rmse 0.01315889  max resid 0.03312343 
    ## Run 218 stress 0.06942776 
    ## ... Procrustes: rmse 2.644848e-05  max resid 6.943209e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.07428313 
    ## Run 220 stress 0.06942778 
    ## ... Procrustes: rmse 6.438157e-05  max resid 0.0001668069 
    ## ... Similar to previous best
    ## Run 221 stress 0.06942776 
    ## ... Procrustes: rmse 6.582763e-06  max resid 1.571026e-05 
    ## ... Similar to previous best
    ## Run 222 stress 0.06978201 
    ## ... Procrustes: rmse 0.0131595  max resid 0.03312473 
    ## Run 223 stress 0.08340287 
    ## Run 224 stress 0.08340288 
    ## Run 225 stress 0.08448436 
    ## Run 226 stress 0.07970525 
    ## Run 227 stress 0.08340287 
    ## Run 228 stress 0.06942776 
    ## ... Procrustes: rmse 5.972338e-06  max resid 1.434991e-05 
    ## ... Similar to previous best
    ## Run 229 stress 0.08340293 
    ## Run 230 stress 0.06978196 
    ## ... Procrustes: rmse 0.01317552  max resid 0.03315918 
    ## Run 231 stress 0.08448436 
    ## Run 232 stress 0.08340287 
    ## Run 233 stress 0.07428315 
    ## Run 234 stress 0.07250813 
    ## Run 235 stress 0.06942776 
    ## ... Procrustes: rmse 4.280303e-06  max resid 1.33372e-05 
    ## ... Similar to previous best
    ## Run 236 stress 0.08448451 
    ## Run 237 stress 0.08340287 
    ## Run 238 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324715  max resid 0.03331084 
    ## Run 239 stress 0.07844924 
    ## Run 240 stress 0.07844912 
    ## Run 241 stress 0.07428321 
    ## Run 242 stress 0.06942776 
    ## ... Procrustes: rmse 1.7692e-05  max resid 4.659527e-05 
    ## ... Similar to previous best
    ## Run 243 stress 0.08340288 
    ## Run 244 stress 0.07250812 
    ## Run 245 stress 0.08340287 
    ## Run 246 stress 0.07250812 
    ## Run 247 stress 0.07970525 
    ## Run 248 stress 0.06978201 
    ## ... Procrustes: rmse 0.01328552  max resid 0.03339165 
    ## Run 249 stress 0.07428313 
    ## Run 250 stress 0.07428313 
    ## Run 251 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325023  max resid 0.0333195 
    ## Run 252 stress 0.07970525 
    ## Run 253 stress 0.07844944 
    ## Run 254 stress 0.07970525 
    ## Run 255 stress 0.07250812 
    ## Run 256 stress 0.08448447 
    ## Run 257 stress 0.07250813 
    ## Run 258 stress 0.07970525 
    ## Run 259 stress 0.06942776 
    ## ... Procrustes: rmse 1.489185e-06  max resid 4.01204e-06 
    ## ... Similar to previous best
    ## Run 260 stress 0.07250814 
    ## Run 261 stress 0.07428318 
    ## Run 262 stress 0.08448444 
    ## Run 263 stress 0.07428313 
    ## Run 264 stress 0.07250812 
    ## Run 265 stress 0.07428313 
    ## Run 266 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323149  max resid 0.03328015 
    ## Run 267 stress 0.07428313 
    ## Run 268 stress 0.07428316 
    ## Run 269 stress 0.0834029 
    ## Run 270 stress 0.06978211 
    ## ... Procrustes: rmse 0.01320949  max resid 0.03324371 
    ## Run 271 stress 0.08448451 
    ## Run 272 stress 0.06942776 
    ## ... Procrustes: rmse 1.275884e-05  max resid 3.322663e-05 
    ## ... Similar to previous best
    ## Run 273 stress 0.07844938 
    ## Run 274 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325542  max resid 0.03333054 
    ## Run 275 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321683  max resid 0.03324895 
    ## Run 276 stress 0.07250816 
    ## Run 277 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323392  max resid 0.03328674 
    ## Run 278 stress 0.06942776 
    ## ... Procrustes: rmse 7.814045e-06  max resid 1.907844e-05 
    ## ... Similar to previous best
    ## Run 279 stress 0.08340287 
    ## Run 280 stress 0.07844937 
    ## Run 281 stress 0.07970525 
    ## Run 282 stress 0.06942777 
    ## ... Procrustes: rmse 1.84793e-05  max resid 5.552474e-05 
    ## ... Similar to previous best
    ## Run 283 stress 0.08340291 
    ## Run 284 stress 0.07428314 
    ## Run 285 stress 0.07428313 
    ## Run 286 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324719  max resid 0.0333144 
    ## Run 287 stress 0.07250814 
    ## Run 288 stress 0.07250813 
    ## Run 289 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324926  max resid 0.03331641 
    ## Run 290 stress 0.07970526 
    ## Run 291 stress 0.06942776 
    ## ... Procrustes: rmse 5.88945e-06  max resid 1.428204e-05 
    ## ... Similar to previous best
    ## Run 292 stress 0.07428313 
    ## Run 293 stress 0.07250813 
    ## Run 294 stress 0.06978203 
    ## ... Procrustes: rmse 0.01328808  max resid 0.03339248 
    ## Run 295 stress 0.08340286 
    ## Run 296 stress 0.08340289 
    ## Run 297 stress 0.07428313 
    ## Run 298 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323616  max resid 0.03328735 
    ## Run 299 stress 0.07428315 
    ## Run 300 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324773  max resid 0.03331143 
    ## Run 301 stress 0.07428313 
    ## Run 302 stress 0.07428313 
    ## Run 303 stress 0.07428314 
    ## Run 304 stress 0.06942776 
    ## ... Procrustes: rmse 1.181059e-05  max resid 2.951206e-05 
    ## ... Similar to previous best
    ## Run 305 stress 0.08448435 
    ## Run 306 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323565  max resid 0.03328834 
    ## Run 307 stress 0.07970525 
    ## Run 308 stress 0.07970526 
    ## Run 309 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326778  max resid 0.03335481 
    ## Run 310 stress 0.07428316 
    ## Run 311 stress 0.07844927 
    ## Run 312 stress 0.08340287 
    ## Run 313 stress 0.06942776 
    ## ... Procrustes: rmse 2.062014e-06  max resid 5.792973e-06 
    ## ... Similar to previous best
    ## Run 314 stress 0.07250812 
    ## Run 315 stress 0.07250813 
    ## Run 316 stress 0.06942778 
    ## ... Procrustes: rmse 5.511837e-05  max resid 0.0001422516 
    ## ... Similar to previous best
    ## Run 317 stress 0.07250812 
    ## Run 318 stress 0.08448436 
    ## Run 319 stress 0.0844844 
    ## Run 320 stress 0.07970526 
    ## Run 321 stress 0.08448443 
    ## Run 322 stress 0.07250813 
    ## Run 323 stress 0.07250812 
    ## Run 324 stress 0.08340289 
    ## Run 325 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323468  max resid 0.03328794 
    ## Run 326 stress 0.07428312 
    ## Run 327 stress 0.07970526 
    ## Run 328 stress 0.06942776 
    ## ... Procrustes: rmse 3.131572e-06  max resid 5.176852e-06 
    ## ... Similar to previous best
    ## Run 329 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322458  max resid 0.03326578 
    ## Run 330 stress 0.08340287 
    ## Run 331 stress 0.07428313 
    ## Run 332 stress 0.06942776 
    ## ... Procrustes: rmse 1.679348e-05  max resid 4.442838e-05 
    ## ... Similar to previous best
    ## Run 333 stress 0.07844945 
    ## Run 334 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326084  max resid 0.03334356 
    ## Run 335 stress 0.07428318 
    ## Run 336 stress 0.07970525 
    ## Run 337 stress 0.07250821 
    ## Run 338 stress 0.08448443 
    ## Run 339 stress 0.07428314 
    ## Run 340 stress 0.07970526 
    ## Run 341 stress 0.07428314 
    ## Run 342 stress 0.07970525 
    ## Run 343 stress 0.07428317 
    ## Run 344 stress 0.06942777 
    ## ... Procrustes: rmse 3.031109e-05  max resid 7.984066e-05 
    ## ... Similar to previous best
    ## Run 345 stress 0.07250812 
    ## Run 346 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324092  max resid 0.03330425 
    ## Run 347 stress 0.07428314 
    ## Run 348 stress 0.08340286 
    ## Run 349 stress 0.06942777 
    ## ... Procrustes: rmse 4.228997e-05  max resid 0.0001099847 
    ## ... Similar to previous best
    ## Run 350 stress 0.06942776 
    ## ... Procrustes: rmse 4.527004e-06  max resid 1.251508e-05 
    ## ... Similar to previous best
    ## Run 351 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323109  max resid 0.03327866 
    ## Run 352 stress 0.0784495 
    ## Run 353 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325098  max resid 0.03332223 
    ## Run 354 stress 0.07970525 
    ## Run 355 stress 0.07970526 
    ## Run 356 stress 0.07428315 
    ## Run 357 stress 0.07428313 
    ## Run 358 stress 0.06942782 
    ## ... Procrustes: rmse 4.238439e-05  max resid 0.0001189512 
    ## ... Similar to previous best
    ## Run 359 stress 0.06942776 
    ## ... Procrustes: rmse 7.183125e-06  max resid 1.896618e-05 
    ## ... Similar to previous best
    ## Run 360 stress 0.07844913 
    ## Run 361 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327113  max resid 0.03336337 
    ## Run 362 stress 0.0834029 
    ## Run 363 stress 0.07428313 
    ## Run 364 stress 0.07250812 
    ## Run 365 stress 0.07428313 
    ## Run 366 stress 0.07844946 
    ## Run 367 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.776373e-06  max resid 5.248332e-06 
    ## ... Similar to previous best
    ## Run 368 stress 0.08340287 
    ## Run 369 stress 0.07844913 
    ## Run 370 stress 0.06942777 
    ## ... Procrustes: rmse 4.52313e-05  max resid 0.0001165416 
    ## ... Similar to previous best
    ## Run 371 stress 0.07428317 
    ## Run 372 stress 0.07250812 
    ## Run 373 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323583  max resid 0.03328801 
    ## Run 374 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319893  max resid 0.03320905 
    ## Run 375 stress 0.07428313 
    ## Run 376 stress 0.06942778 
    ## ... Procrustes: rmse 4.773985e-05  max resid 0.0001227761 
    ## ... Similar to previous best
    ## Run 377 stress 0.07250813 
    ## Run 378 stress 0.08340288 
    ## Run 379 stress 0.06942776 
    ## ... Procrustes: rmse 1.845775e-05  max resid 4.815124e-05 
    ## ... Similar to previous best
    ## Run 380 stress 0.08340291 
    ## Run 381 stress 0.07428313 
    ## Run 382 stress 0.07970526 
    ## Run 383 stress 0.06978199 
    ## ... Procrustes: rmse 0.01316706  max resid 0.03314155 
    ## Run 384 stress 0.0784493 
    ## Run 385 stress 0.07250814 
    ## Run 386 stress 0.06942777 
    ## ... Procrustes: rmse 4.304785e-05  max resid 0.000111264 
    ## ... Similar to previous best
    ## Run 387 stress 0.07250813 
    ## Run 388 stress 0.07428313 
    ## Run 389 stress 0.07250812 
    ## Run 390 stress 0.07428313 
    ## Run 391 stress 0.07970526 
    ## Run 392 stress 0.08340286 
    ## Run 393 stress 0.07428319 
    ## Run 394 stress 0.0834029 
    ## Run 395 stress 0.07428313 
    ## Run 396 stress 0.07428314 
    ## Run 397 stress 0.06942776 
    ## ... Procrustes: rmse 6.892623e-06  max resid 1.15853e-05 
    ## ... Similar to previous best
    ## Run 398 stress 0.07250814 
    ## Run 399 stress 0.07970525 
    ## Run 400 stress 0.07250814 
    ## Run 401 stress 0.08340298 
    ## Run 402 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324088  max resid 0.0332986 
    ## Run 403 stress 0.06978192 
    ## ... Procrustes: rmse 0.01318952  max resid 0.03318866 
    ## Run 404 stress 0.08448441 
    ## Run 405 stress 0.06942777 
    ## ... Procrustes: rmse 3.17392e-05  max resid 8.395733e-05 
    ## ... Similar to previous best
    ## Run 406 stress 0.06942777 
    ## ... Procrustes: rmse 3.50898e-05  max resid 9.173645e-05 
    ## ... Similar to previous best
    ## Run 407 stress 0.07428321 
    ## Run 408 stress 0.07428316 
    ## Run 409 stress 0.07250812 
    ## Run 410 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325163  max resid 0.03332351 
    ## Run 411 stress 0.07428315 
    ## Run 412 stress 0.06942777 
    ## ... Procrustes: rmse 2.726289e-05  max resid 7.189077e-05 
    ## ... Similar to previous best
    ## Run 413 stress 0.06942776 
    ## ... Procrustes: rmse 2.548848e-05  max resid 6.578845e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.06942778 
    ## ... Procrustes: rmse 6.126673e-05  max resid 0.0001581035 
    ## ... Similar to previous best
    ## Run 415 stress 0.07428314 
    ## Run 416 stress 0.07970525 
    ## Run 417 stress 0.08448452 
    ## Run 418 stress 0.07428312 
    ## Run 419 stress 0.07970526 
    ## Run 420 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132313  max resid 0.0332779 
    ## Run 421 stress 0.07970525 
    ## Run 422 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326107  max resid 0.03333991 
    ## Run 423 stress 0.07428313 
    ## Run 424 stress 0.07428313 
    ## Run 425 stress 0.07250813 
    ## Run 426 stress 0.07250812 
    ## Run 427 stress 0.06942776 
    ## ... Procrustes: rmse 4.744948e-06  max resid 1.061901e-05 
    ## ... Similar to previous best
    ## Run 428 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326109  max resid 0.03333726 
    ## Run 429 stress 0.06942776 
    ## ... Procrustes: rmse 4.538179e-06  max resid 1.140416e-05 
    ## ... Similar to previous best
    ## Run 430 stress 0.07428316 
    ## Run 431 stress 0.07250814 
    ## Run 432 stress 0.08340288 
    ## Run 433 stress 0.07250812 
    ## Run 434 stress 0.07428313 
    ## Run 435 stress 0.08448443 
    ## Run 436 stress 0.07428313 
    ## Run 437 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324925  max resid 0.03331418 
    ## Run 438 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324495  max resid 0.0333063 
    ## Run 439 stress 0.06942776 
    ## ... Procrustes: rmse 1.241303e-05  max resid 3.181727e-05 
    ## ... Similar to previous best
    ## Run 440 stress 0.07250812 
    ## Run 441 stress 0.08340287 
    ## Run 442 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324822  max resid 0.03331337 
    ## Run 443 stress 0.06942776 
    ## ... Procrustes: rmse 2.802109e-06  max resid 5.70225e-06 
    ## ... Similar to previous best
    ## Run 444 stress 0.07970526 
    ## Run 445 stress 0.08448449 
    ## Run 446 stress 0.08340287 
    ## Run 447 stress 0.06942776 
    ## ... Procrustes: rmse 8.422739e-06  max resid 2.168266e-05 
    ## ... Similar to previous best
    ## Run 448 stress 0.07250812 
    ## Run 449 stress 0.06942776 
    ## ... Procrustes: rmse 7.680043e-06  max resid 2.05796e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.07428323 
    ## Run 451 stress 0.07428313 
    ## Run 452 stress 0.06942776 
    ## ... Procrustes: rmse 1.965804e-05  max resid 5.061619e-05 
    ## ... Similar to previous best
    ## Run 453 stress 0.06942776 
    ## ... Procrustes: rmse 2.01134e-05  max resid 5.175524e-05 
    ## ... Similar to previous best
    ## Run 454 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324247  max resid 0.03330132 
    ## Run 455 stress 0.07428313 
    ## Run 456 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326864  max resid 0.03335453 
    ## Run 457 stress 0.06942776 
    ## ... Procrustes: rmse 1.4494e-05  max resid 3.725295e-05 
    ## ... Similar to previous best
    ## Run 458 stress 0.08340288 
    ## Run 459 stress 0.07428313 
    ## Run 460 stress 0.08340287 
    ## Run 461 stress 0.07428313 
    ## Run 462 stress 0.07428314 
    ## Run 463 stress 0.07428312 
    ## Run 464 stress 0.07428313 
    ## Run 465 stress 0.08340288 
    ## Run 466 stress 0.07970525 
    ## Run 467 stress 0.06942776 
    ## ... Procrustes: rmse 6.644093e-06  max resid 1.592033e-05 
    ## ... Similar to previous best
    ## Run 468 stress 0.08448439 
    ## Run 469 stress 0.07970525 
    ## Run 470 stress 0.07970526 
    ## Run 471 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324517  max resid 0.03330453 
    ## Run 472 stress 0.07844946 
    ## Run 473 stress 0.06942776 
    ## ... Procrustes: rmse 5.179237e-06  max resid 1.235051e-05 
    ## ... Similar to previous best
    ## Run 474 stress 0.06942777 
    ## ... Procrustes: rmse 2.457869e-05  max resid 6.543924e-05 
    ## ... Similar to previous best
    ## Run 475 stress 0.07428313 
    ## Run 476 stress 0.08340292 
    ## Run 477 stress 0.08448447 
    ## Run 478 stress 0.06978194 
    ## ... Procrustes: rmse 0.01318246  max resid 0.03317484 
    ## Run 479 stress 0.07428313 
    ## Run 480 stress 0.07428313 
    ## Run 481 stress 0.07970525 
    ## Run 482 stress 0.07428313 
    ## Run 483 stress 0.07428318 
    ## Run 484 stress 0.06942776 
    ## ... Procrustes: rmse 1.933429e-05  max resid 4.988539e-05 
    ## ... Similar to previous best
    ## Run 485 stress 0.06942776 
    ## ... Procrustes: rmse 1.097239e-05  max resid 2.808806e-05 
    ## ... Similar to previous best
    ## Run 486 stress 0.06942776 
    ## ... Procrustes: rmse 8.92535e-06  max resid 2.194563e-05 
    ## ... Similar to previous best
    ## Run 487 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327164  max resid 0.03336073 
    ## Run 488 stress 0.07250813 
    ## Run 489 stress 0.07428314 
    ## Run 490 stress 0.07428318 
    ## Run 491 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317156  max resid 0.03314975 
    ## Run 492 stress 0.07428315 
    ## Run 493 stress 0.07250813 
    ## Run 494 stress 0.07970525 
    ## Run 495 stress 0.07250812 
    ## Run 496 stress 0.07250814 
    ## Run 497 stress 0.07970526 
    ## Run 498 stress 0.08340288 
    ## Run 499 stress 0.08448437 
    ## Run 500 stress 0.07970526 
    ## *** Best solution repeated 26 times

``` r
# Stratified lakes
SD_beta_S_NMDS <- metaMDS(SD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1410449 
    ## Run 1 stress 0.1572668 
    ## Run 2 stress 0.1551863 
    ## Run 3 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2673133  max resid 0.5425958 
    ## Run 4 stress 0.1415299 
    ## Run 5 stress 0.1771969 
    ## Run 6 stress 0.1663544 
    ## Run 7 stress 0.1771969 
    ## Run 8 stress 0.1320258 
    ## ... Procrustes: rmse 8.73086e-07  max resid 1.774296e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.1693717 
    ## Run 10 stress 0.1410451 
    ## Run 11 stress 0.1587623 
    ## Run 12 stress 0.1929519 
    ## Run 13 stress 0.1320258 
    ## ... Procrustes: rmse 9.636628e-07  max resid 1.749835e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.1587623 
    ## Run 15 stress 0.1771969 
    ## Run 16 stress 0.1551863 
    ## Run 17 stress 0.2687485 
    ## Run 18 stress 0.1572668 
    ## Run 19 stress 0.1572668 
    ## Run 20 stress 0.141045 
    ## Run 21 stress 0.1572668 
    ## Run 22 stress 0.2501064 
    ## Run 23 stress 0.1771966 
    ## Run 24 stress 0.2514078 
    ## Run 25 stress 0.1445811 
    ## Run 26 stress 0.1415299 
    ## Run 27 stress 0.1628606 
    ## Run 28 stress 0.1410455 
    ## Run 29 stress 0.200816 
    ## Run 30 stress 0.1415299 
    ## Run 31 stress 0.1407206 
    ## Run 32 stress 0.1796865 
    ## Run 33 stress 0.1320258 
    ## ... Procrustes: rmse 2.366243e-06  max resid 4.877333e-06 
    ## ... Similar to previous best
    ## Run 34 stress 0.1693717 
    ## Run 35 stress 0.1415299 
    ## Run 36 stress 0.2496382 
    ## Run 37 stress 0.1572668 
    ## Run 38 stress 0.1407298 
    ## Run 39 stress 0.1320258 
    ## ... Procrustes: rmse 1.516765e-06  max resid 2.295735e-06 
    ## ... Similar to previous best
    ## Run 40 stress 0.1407298 
    ## Run 41 stress 0.1320258 
    ## ... Procrustes: rmse 6.494454e-06  max resid 1.290495e-05 
    ## ... Similar to previous best
    ## Run 42 stress 0.1415299 
    ## Run 43 stress 0.2711712 
    ## Run 44 stress 0.1320258 
    ## ... Procrustes: rmse 1.28658e-06  max resid 2.361999e-06 
    ## ... Similar to previous best
    ## Run 45 stress 0.1830332 
    ## Run 46 stress 0.1415299 
    ## Run 47 stress 0.1407208 
    ## Run 48 stress 0.1407298 
    ## Run 49 stress 0.1407298 
    ## Run 50 stress 0.2496385 
    ## Run 51 stress 0.1802587 
    ## Run 52 stress 0.2510079 
    ## Run 53 stress 0.1929519 
    ## Run 54 stress 0.1320258 
    ## ... Procrustes: rmse 1.128831e-06  max resid 2.2269e-06 
    ## ... Similar to previous best
    ## Run 55 stress 0.1970303 
    ## Run 56 stress 0.1693717 
    ## Run 57 stress 0.2073091 
    ## Run 58 stress 0.1415299 
    ## Run 59 stress 0.2385029 
    ## Run 60 stress 0.1598154 
    ## Run 61 stress 0.1929519 
    ## Run 62 stress 0.2085297 
    ## Run 63 stress 0.2385029 
    ## Run 64 stress 0.1407207 
    ## Run 65 stress 0.1572668 
    ## Run 66 stress 0.1628606 
    ## Run 67 stress 0.1383681 
    ## Run 68 stress 0.1320258 
    ## ... Procrustes: rmse 5.009096e-07  max resid 9.694989e-07 
    ## ... Similar to previous best
    ## Run 69 stress 0.1320258 
    ## ... Procrustes: rmse 4.982985e-07  max resid 9.135885e-07 
    ## ... Similar to previous best
    ## Run 70 stress 0.1383681 
    ## Run 71 stress 0.1628606 
    ## Run 72 stress 0.1572668 
    ## Run 73 stress 0.2514078 
    ## Run 74 stress 0.1410455 
    ## Run 75 stress 0.1320258 
    ## ... Procrustes: rmse 2.527955e-06  max resid 5.177629e-06 
    ## ... Similar to previous best
    ## Run 76 stress 0.3083098 
    ## Run 77 stress 0.1383681 
    ## Run 78 stress 0.1415299 
    ## Run 79 stress 0.141045 
    ## Run 80 stress 0.1551863 
    ## Run 81 stress 0.1434197 
    ## Run 82 stress 0.1320258 
    ## ... Procrustes: rmse 1.128978e-06  max resid 2.327319e-06 
    ## ... Similar to previous best
    ## Run 83 stress 0.1572668 
    ## Run 84 stress 0.1407207 
    ## Run 85 stress 0.1929519 
    ## Run 86 stress 0.1410452 
    ## Run 87 stress 0.1410448 
    ## Run 88 stress 0.1410451 
    ## Run 89 stress 0.1998955 
    ## Run 90 stress 0.2073091 
    ## Run 91 stress 0.1572668 
    ## Run 92 stress 0.1410452 
    ## Run 93 stress 0.1583016 
    ## Run 94 stress 0.2842805 
    ## Run 95 stress 0.1434197 
    ## Run 96 stress 0.2008162 
    ## Run 97 stress 0.2227526 
    ## Run 98 stress 0.1383681 
    ## Run 99 stress 0.257865 
    ## Run 100 stress 0.1581642 
    ## Run 101 stress 0.3083108 
    ## Run 102 stress 0.1320258 
    ## ... Procrustes: rmse 8.945042e-07  max resid 1.557005e-06 
    ## ... Similar to previous best
    ## Run 103 stress 0.1802752 
    ## Run 104 stress 0.1796865 
    ## Run 105 stress 0.1383681 
    ## Run 106 stress 0.1551863 
    ## Run 107 stress 0.1410449 
    ## Run 108 stress 0.2615852 
    ## Run 109 stress 0.1583016 
    ## Run 110 stress 0.1572668 
    ## Run 111 stress 0.2008162 
    ## Run 112 stress 0.2073091 
    ## Run 113 stress 0.1410448 
    ## Run 114 stress 0.1663544 
    ## Run 115 stress 0.1407298 
    ## Run 116 stress 0.1383682 
    ## Run 117 stress 0.1830332 
    ## Run 118 stress 0.1320258 
    ## ... Procrustes: rmse 1.000322e-06  max resid 1.725805e-06 
    ## ... Similar to previous best
    ## Run 119 stress 0.1383681 
    ## Run 120 stress 0.1410452 
    ## Run 121 stress 0.1415299 
    ## Run 122 stress 0.1693717 
    ## Run 123 stress 0.3083098 
    ## Run 124 stress 0.1383683 
    ## Run 125 stress 0.1415299 
    ## Run 126 stress 0.2501065 
    ## Run 127 stress 0.1598154 
    ## Run 128 stress 0.2422668 
    ## Run 129 stress 0.1410449 
    ## Run 130 stress 0.1415299 
    ## Run 131 stress 0.1407207 
    ## Run 132 stress 0.2600627 
    ## Run 133 stress 0.1407207 
    ## Run 134 stress 0.1407298 
    ## Run 135 stress 0.2422742 
    ## Run 136 stress 0.2008166 
    ## Run 137 stress 0.3083098 
    ## Run 138 stress 0.1415299 
    ## Run 139 stress 0.1410448 
    ## Run 140 stress 0.1663544 
    ## Run 141 stress 0.1320258 
    ## ... Procrustes: rmse 9.690482e-07  max resid 2.044611e-06 
    ## ... Similar to previous best
    ## Run 142 stress 0.1551863 
    ## Run 143 stress 0.1415299 
    ## Run 144 stress 0.1771969 
    ## Run 145 stress 0.1415299 
    ## Run 146 stress 0.1407298 
    ## Run 147 stress 0.1445811 
    ## Run 148 stress 0.1415299 
    ## Run 149 stress 0.1320258 
    ## ... Procrustes: rmse 1.890973e-06  max resid 3.828569e-06 
    ## ... Similar to previous best
    ## Run 150 stress 0.1693717 
    ## Run 151 stress 0.1410449 
    ## Run 152 stress 0.2385029 
    ## Run 153 stress 0.1572668 
    ## Run 154 stress 0.1598154 
    ## Run 155 stress 0.2073091 
    ## Run 156 stress 0.1572668 
    ## Run 157 stress 0.1383682 
    ## Run 158 stress 0.2510079 
    ## Run 159 stress 0.1383681 
    ## Run 160 stress 0.1998955 
    ## Run 161 stress 0.1320258 
    ## ... Procrustes: rmse 6.582705e-07  max resid 1.229436e-06 
    ## ... Similar to previous best
    ## Run 162 stress 0.1583016 
    ## Run 163 stress 0.1410454 
    ## Run 164 stress 0.1551863 
    ## Run 165 stress 0.1962476 
    ## Run 166 stress 0.1383681 
    ## Run 167 stress 0.1415299 
    ## Run 168 stress 0.2538018 
    ## Run 169 stress 0.1407298 
    ## Run 170 stress 0.1830332 
    ## Run 171 stress 0.2005461 
    ## Run 172 stress 0.1551863 
    ## Run 173 stress 0.1383681 
    ## Run 174 stress 0.1415299 
    ## Run 175 stress 0.1383681 
    ## Run 176 stress 0.1320258 
    ## ... Procrustes: rmse 7.873329e-07  max resid 1.362168e-06 
    ## ... Similar to previous best
    ## Run 177 stress 0.1628606 
    ## Run 178 stress 0.2615852 
    ## Run 179 stress 0.1383682 
    ## Run 180 stress 0.2848522 
    ## Run 181 stress 0.2373121 
    ## Run 182 stress 0.2456775 
    ## Run 183 stress 0.1320258 
    ## ... Procrustes: rmse 1.068632e-06  max resid 1.931057e-06 
    ## ... Similar to previous best
    ## Run 184 stress 0.1410448 
    ## Run 185 stress 0.1970303 
    ## Run 186 stress 0.1407298 
    ## Run 187 stress 0.1415299 
    ## Run 188 stress 0.1410448 
    ## Run 189 stress 0.2008158 
    ## Run 190 stress 0.1320258 
    ## ... Procrustes: rmse 1.370158e-06  max resid 2.041027e-06 
    ## ... Similar to previous best
    ## Run 191 stress 0.1663544 
    ## Run 192 stress 0.1693717 
    ## Run 193 stress 0.141045 
    ## Run 194 stress 0.1320258 
    ## ... Procrustes: rmse 1.161777e-06  max resid 2.284659e-06 
    ## ... Similar to previous best
    ## Run 195 stress 0.1407298 
    ## Run 196 stress 0.1771966 
    ## Run 197 stress 0.1415299 
    ## Run 198 stress 0.1929519 
    ## Run 199 stress 0.1551863 
    ## Run 200 stress 0.1383681 
    ## Run 201 stress 0.1663544 
    ## Run 202 stress 0.2118436 
    ## Run 203 stress 0.1415299 
    ## Run 204 stress 0.1410448 
    ## Run 205 stress 0.1410449 
    ## Run 206 stress 0.1383682 
    ## Run 207 stress 0.1383681 
    ## Run 208 stress 0.1320258 
    ## ... Procrustes: rmse 2.011657e-06  max resid 4.121472e-06 
    ## ... Similar to previous best
    ## Run 209 stress 0.2456775 
    ## Run 210 stress 0.1415299 
    ## Run 211 stress 0.1383681 
    ## Run 212 stress 0.1551863 
    ## Run 213 stress 0.300185 
    ## Run 214 stress 0.141045 
    ## Run 215 stress 0.2962402 
    ## Run 216 stress 0.1434197 
    ## Run 217 stress 0.1410448 
    ## Run 218 stress 0.2615852 
    ## Run 219 stress 0.1929519 
    ## Run 220 stress 0.1407298 
    ## Run 221 stress 0.1407298 
    ## Run 222 stress 0.2456775 
    ## Run 223 stress 0.1970303 
    ## Run 224 stress 0.1628606 
    ## Run 225 stress 0.1551863 
    ## Run 226 stress 0.2612821 
    ## Run 227 stress 0.1583016 
    ## Run 228 stress 0.1407298 
    ## Run 229 stress 0.1830332 
    ## Run 230 stress 0.2227526 
    ## Run 231 stress 0.1830332 
    ## Run 232 stress 0.1410448 
    ## Run 233 stress 0.2683391 
    ## Run 234 stress 0.1663544 
    ## Run 235 stress 0.1415299 
    ## Run 236 stress 0.2501068 
    ## Run 237 stress 0.2550046 
    ## Run 238 stress 0.2422667 
    ## Run 239 stress 0.1407298 
    ## Run 240 stress 0.1587623 
    ## Run 241 stress 0.1410451 
    ## Run 242 stress 0.1383681 
    ## Run 243 stress 0.2550043 
    ## Run 244 stress 0.1320258 
    ## ... Procrustes: rmse 8.749755e-07  max resid 1.659166e-06 
    ## ... Similar to previous best
    ## Run 245 stress 0.1383681 
    ## Run 246 stress 0.1663544 
    ## Run 247 stress 0.2514078 
    ## Run 248 stress 0.1320258 
    ## ... Procrustes: rmse 1.634753e-06  max resid 3.191281e-06 
    ## ... Similar to previous best
    ## Run 249 stress 0.3002107 
    ## Run 250 stress 0.1410449 
    ## Run 251 stress 0.1383681 
    ## Run 252 stress 0.2373121 
    ## Run 253 stress 0.1320258 
    ## ... Procrustes: rmse 2.036872e-06  max resid 4.227717e-06 
    ## ... Similar to previous best
    ## Run 254 stress 0.1693717 
    ## Run 255 stress 0.1410448 
    ## Run 256 stress 0.1383682 
    ## Run 257 stress 0.1434197 
    ## Run 258 stress 0.2615852 
    ## Run 259 stress 0.1410448 
    ## Run 260 stress 0.140721 
    ## Run 261 stress 0.1407207 
    ## Run 262 stress 0.1383682 
    ## Run 263 stress 0.1407207 
    ## Run 264 stress 0.141045 
    ## Run 265 stress 0.1410453 
    ## Run 266 stress 0.1383681 
    ## Run 267 stress 0.1410448 
    ## Run 268 stress 0.141045 
    ## Run 269 stress 0.1415299 
    ## Run 270 stress 0.1410454 
    ## Run 271 stress 0.1415299 
    ## Run 272 stress 0.1320258 
    ## ... Procrustes: rmse 2.584192e-06  max resid 5.264979e-06 
    ## ... Similar to previous best
    ## Run 273 stress 0.1410451 
    ## Run 274 stress 0.141045 
    ## Run 275 stress 0.1663544 
    ## Run 276 stress 0.1407298 
    ## Run 277 stress 0.3083098 
    ## Run 278 stress 0.1320258 
    ## ... Procrustes: rmse 8.587582e-07  max resid 1.713168e-06 
    ## ... Similar to previous best
    ## Run 279 stress 0.1407207 
    ## Run 280 stress 0.1407298 
    ## Run 281 stress 0.1572668 
    ## Run 282 stress 0.1776793 
    ## Run 283 stress 0.1796865 
    ## Run 284 stress 0.1771966 
    ## Run 285 stress 0.3083098 
    ## Run 286 stress 0.1693717 
    ## Run 287 stress 0.1693717 
    ## Run 288 stress 0.1796864 
    ## Run 289 stress 0.1320258 
    ## ... Procrustes: rmse 3.671632e-07  max resid 6.565061e-07 
    ## ... Similar to previous best
    ## Run 290 stress 0.141045 
    ## Run 291 stress 0.1796865 
    ## Run 292 stress 0.1415299 
    ## Run 293 stress 0.2422667 
    ## Run 294 stress 0.1415299 
    ## Run 295 stress 0.2501066 
    ## Run 296 stress 0.1407207 
    ## Run 297 stress 0.1320258 
    ## ... Procrustes: rmse 1.174146e-06  max resid 2.321989e-06 
    ## ... Similar to previous best
    ## Run 298 stress 0.2227526 
    ## Run 299 stress 0.1771966 
    ## Run 300 stress 0.1407298 
    ## Run 301 stress 0.1663544 
    ## Run 302 stress 0.1320258 
    ## ... Procrustes: rmse 1.10261e-06  max resid 2.201319e-06 
    ## ... Similar to previous best
    ## Run 303 stress 0.1383681 
    ## Run 304 stress 0.1551863 
    ## Run 305 stress 0.1320258 
    ## ... Procrustes: rmse 1.302169e-06  max resid 2.565257e-06 
    ## ... Similar to previous best
    ## Run 306 stress 0.1320258 
    ## ... Procrustes: rmse 9.836742e-07  max resid 1.951641e-06 
    ## ... Similar to previous best
    ## Run 307 stress 0.1551863 
    ## Run 308 stress 0.1383681 
    ## Run 309 stress 0.1410453 
    ## Run 310 stress 0.2538018 
    ## Run 311 stress 0.2612821 
    ## Run 312 stress 0.1383682 
    ## Run 313 stress 0.1415299 
    ## Run 314 stress 0.2562346 
    ## Run 315 stress 0.1320258 
    ## ... Procrustes: rmse 1.180088e-06  max resid 2.404674e-06 
    ## ... Similar to previous best
    ## Run 316 stress 0.1598154 
    ## Run 317 stress 0.1410452 
    ## Run 318 stress 0.1407298 
    ## Run 319 stress 0.1445811 
    ## Run 320 stress 0.1572668 
    ## Run 321 stress 0.2008163 
    ## Run 322 stress 0.1434197 
    ## Run 323 stress 0.1320258 
    ## ... Procrustes: rmse 5.068038e-07  max resid 6.959846e-07 
    ## ... Similar to previous best
    ## Run 324 stress 0.1551863 
    ## Run 325 stress 0.2385029 
    ## Run 326 stress 0.1434197 
    ## Run 327 stress 0.1434197 
    ## Run 328 stress 0.1693717 
    ## Run 329 stress 0.141045 
    ## Run 330 stress 0.1383681 
    ## Run 331 stress 0.1572668 
    ## Run 332 stress 0.1320258 
    ## ... Procrustes: rmse 2.05528e-06  max resid 4.043447e-06 
    ## ... Similar to previous best
    ## Run 333 stress 0.1693717 
    ## Run 334 stress 0.1320258 
    ## ... Procrustes: rmse 5.25992e-07  max resid 8.386666e-07 
    ## ... Similar to previous best
    ## Run 335 stress 0.2008162 
    ## Run 336 stress 0.1415299 
    ## Run 337 stress 0.1572668 
    ## Run 338 stress 0.1383682 
    ## Run 339 stress 0.1383682 
    ## Run 340 stress 0.1407207 
    ## Run 341 stress 0.1320258 
    ## ... Procrustes: rmse 2.72911e-06  max resid 5.642821e-06 
    ## ... Similar to previous best
    ## Run 342 stress 0.1410449 
    ## Run 343 stress 0.2848516 
    ## Run 344 stress 0.1771969 
    ## Run 345 stress 0.1407206 
    ## Run 346 stress 0.2385029 
    ## Run 347 stress 0.1551863 
    ## Run 348 stress 0.3083098 
    ## Run 349 stress 0.1383681 
    ## Run 350 stress 0.141045 
    ## Run 351 stress 0.1415299 
    ## Run 352 stress 0.1383681 
    ## Run 353 stress 0.1383682 
    ## Run 354 stress 0.1415299 
    ## Run 355 stress 0.1771966 
    ## Run 356 stress 0.1970303 
    ## Run 357 stress 0.1970303 
    ## Run 358 stress 0.1410449 
    ## Run 359 stress 0.1802752 
    ## Run 360 stress 0.1572668 
    ## Run 361 stress 0.1551863 
    ## Run 362 stress 0.1320258 
    ## ... Procrustes: rmse 1.040497e-06  max resid 2.069703e-06 
    ## ... Similar to previous best
    ## Run 363 stress 0.1771969 
    ## Run 364 stress 0.1410453 
    ## Run 365 stress 0.1572668 
    ## Run 366 stress 0.1407298 
    ## Run 367 stress 0.1771966 
    ## Run 368 stress 0.1415299 
    ## Run 369 stress 0.1407298 
    ## Run 370 stress 0.1407207 
    ## Run 371 stress 0.1998955 
    ## Run 372 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 3.279102e-07  max resid 5.207305e-07 
    ## ... Similar to previous best
    ## Run 373 stress 0.1998955 
    ## Run 374 stress 0.1415299 
    ## Run 375 stress 0.1407298 
    ## Run 376 stress 0.1572668 
    ## Run 377 stress 0.1410451 
    ## Run 378 stress 0.1830332 
    ## Run 379 stress 0.2612821 
    ## Run 380 stress 0.1320258 
    ## ... Procrustes: rmse 3.176976e-06  max resid 6.40222e-06 
    ## ... Similar to previous best
    ## Run 381 stress 0.2074935 
    ## Run 382 stress 0.1796865 
    ## Run 383 stress 0.1551863 
    ## Run 384 stress 0.1410448 
    ## Run 385 stress 0.3083098 
    ## Run 386 stress 0.1383681 
    ## Run 387 stress 0.1445811 
    ## Run 388 stress 0.2085297 
    ## Run 389 stress 0.2422742 
    ## Run 390 stress 0.1407298 
    ## Run 391 stress 0.1320258 
    ## ... Procrustes: rmse 1.854005e-06  max resid 3.68279e-06 
    ## ... Similar to previous best
    ## Run 392 stress 0.1410448 
    ## Run 393 stress 0.1551863 
    ## Run 394 stress 0.1410453 
    ## Run 395 stress 0.2385029 
    ## Run 396 stress 0.2373121 
    ## Run 397 stress 0.1407209 
    ## Run 398 stress 0.1407298 
    ## Run 399 stress 0.1383681 
    ## Run 400 stress 0.1320258 
    ## ... Procrustes: rmse 2.230768e-06  max resid 4.571131e-06 
    ## ... Similar to previous best
    ## Run 401 stress 0.1383682 
    ## Run 402 stress 0.1407207 
    ## Run 403 stress 0.1445811 
    ## Run 404 stress 0.1410448 
    ## Run 405 stress 0.1929519 
    ## Run 406 stress 0.1320258 
    ## ... Procrustes: rmse 1.320456e-06  max resid 2.626365e-06 
    ## ... Similar to previous best
    ## Run 407 stress 0.1407207 
    ## Run 408 stress 0.2385029 
    ## Run 409 stress 0.1693717 
    ## Run 410 stress 0.1587623 
    ## Run 411 stress 0.1583016 
    ## Run 412 stress 0.1383682 
    ## Run 413 stress 0.2615852 
    ## Run 414 stress 0.1551863 
    ## Run 415 stress 0.1410454 
    ## Run 416 stress 0.1407209 
    ## Run 417 stress 0.1320258 
    ## ... Procrustes: rmse 1.057471e-06  max resid 2.038786e-06 
    ## ... Similar to previous best
    ## Run 418 stress 0.1434197 
    ## Run 419 stress 0.1445811 
    ## Run 420 stress 0.1410451 
    ## Run 421 stress 0.1320258 
    ## ... Procrustes: rmse 1.100759e-06  max resid 2.243869e-06 
    ## ... Similar to previous best
    ## Run 422 stress 0.1410452 
    ## Run 423 stress 0.1551863 
    ## Run 424 stress 0.1410451 
    ## Run 425 stress 0.1410449 
    ## Run 426 stress 0.1551863 
    ## Run 427 stress 0.2385029 
    ## Run 428 stress 0.1415299 
    ## Run 429 stress 0.1628606 
    ## Run 430 stress 0.2711187 
    ## Run 431 stress 0.1551863 
    ## Run 432 stress 0.1320258 
    ## ... Procrustes: rmse 1.907109e-06  max resid 3.818657e-06 
    ## ... Similar to previous best
    ## Run 433 stress 0.1410453 
    ## Run 434 stress 0.1830332 
    ## Run 435 stress 0.1415299 
    ## Run 436 stress 0.1320258 
    ## ... Procrustes: rmse 1.295983e-06  max resid 2.681863e-06 
    ## ... Similar to previous best
    ## Run 437 stress 0.1628606 
    ## Run 438 stress 0.1970303 
    ## Run 439 stress 0.1383681 
    ## Run 440 stress 0.1771969 
    ## Run 441 stress 0.1410448 
    ## Run 442 stress 0.3072832 
    ## Run 443 stress 0.1383681 
    ## Run 444 stress 0.1802752 
    ## Run 445 stress 0.2422667 
    ## Run 446 stress 0.2385029 
    ## Run 447 stress 0.1410454 
    ## Run 448 stress 0.1320258 
    ## ... Procrustes: rmse 1.146506e-06  max resid 2.313977e-06 
    ## ... Similar to previous best
    ## Run 449 stress 0.1407298 
    ## Run 450 stress 0.1434197 
    ## Run 451 stress 0.1320258 
    ## ... Procrustes: rmse 1.878138e-06  max resid 3.685263e-06 
    ## ... Similar to previous best
    ## Run 452 stress 0.1383681 
    ## Run 453 stress 0.1383682 
    ## Run 454 stress 0.1572668 
    ## Run 455 stress 0.1320258 
    ## ... Procrustes: rmse 1.269674e-06  max resid 2.56531e-06 
    ## ... Similar to previous best
    ## Run 456 stress 0.1410453 
    ## Run 457 stress 0.1407298 
    ## Run 458 stress 0.1415299 
    ## Run 459 stress 0.1407298 
    ## Run 460 stress 0.1410449 
    ## Run 461 stress 0.2538018 
    ## Run 462 stress 0.1434197 
    ## Run 463 stress 0.1551863 
    ## Run 464 stress 0.1407298 
    ## Run 465 stress 0.1383681 
    ## Run 466 stress 0.1551863 
    ## Run 467 stress 0.1407298 
    ## Run 468 stress 0.1320258 
    ## ... Procrustes: rmse 3.209718e-06  max resid 6.602061e-06 
    ## ... Similar to previous best
    ## Run 469 stress 0.1415299 
    ## Run 470 stress 0.1415299 
    ## Run 471 stress 0.1407298 
    ## Run 472 stress 0.2422742 
    ## Run 473 stress 0.2852061 
    ## Run 474 stress 0.1445811 
    ## Run 475 stress 0.1407206 
    ## Run 476 stress 0.1551863 
    ## Run 477 stress 0.1407298 
    ## Run 478 stress 0.1962476 
    ## Run 479 stress 0.1407298 
    ## Run 480 stress 0.1970303 
    ## Run 481 stress 0.2008161 
    ## Run 482 stress 0.2842805 
    ## Run 483 stress 0.1383681 
    ## Run 484 stress 0.1410449 
    ## Run 485 stress 0.1583016 
    ## Run 486 stress 0.2422667 
    ## Run 487 stress 0.1551863 
    ## Run 488 stress 0.1445811 
    ## Run 489 stress 0.1415299 
    ## Run 490 stress 0.1445811 
    ## Run 491 stress 0.1445811 
    ## Run 492 stress 0.1407298 
    ## Run 493 stress 0.1415299 
    ## Run 494 stress 0.1693717 
    ## Run 495 stress 0.1320258 
    ## ... Procrustes: rmse 1.431879e-06  max resid 2.993339e-06 
    ## ... Similar to previous best
    ## Run 496 stress 0.1410449 
    ## Run 497 stress 0.1410453 
    ## Run 498 stress 0.1445811 
    ## Run 499 stress 0.1415299 
    ## Run 500 stress 0.1415299 
    ## *** Best solution repeated 14 times

``` r
### Environmental
# Surveyed sites 
SD_beta_env_NMDS <- metaMDS(SD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07669178 
    ## Run 1 stress 0.07871748 
    ## Run 2 stress 0.07871764 
    ## Run 3 stress 0.08471635 
    ## Run 4 stress 0.08453168 
    ## Run 5 stress 0.07669164 
    ## ... New best solution
    ## ... Procrustes: rmse 9.586928e-05  max resid 0.0001631174 
    ## ... Similar to previous best
    ## Run 6 stress 0.07731534 
    ## Run 7 stress 0.08433691 
    ## Run 8 stress 0.08682677 
    ## Run 9 stress 0.07868337 
    ## Run 10 stress 0.08433693 
    ## Run 11 stress 0.07669153 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003206114  max resid 0.0005526275 
    ## ... Similar to previous best
    ## Run 12 stress 0.08104934 
    ## Run 13 stress 0.07951498 
    ## Run 14 stress 0.0852703 
    ## Run 15 stress 0.08684079 
    ## Run 16 stress 0.08175688 
    ## Run 17 stress 0.08104911 
    ## Run 18 stress 0.07868327 
    ## Run 19 stress 0.08459622 
    ## Run 20 stress 0.07948976 
    ## Run 21 stress 0.08104913 
    ## Run 22 stress 0.0810489 
    ## Run 23 stress 0.08305314 
    ## Run 24 stress 0.0847158 
    ## Run 25 stress 0.08175692 
    ## Run 26 stress 0.08499328 
    ## Run 27 stress 0.08842148 
    ## Run 28 stress 0.08175694 
    ## Run 29 stress 0.0831737 
    ## Run 30 stress 0.08012968 
    ## Run 31 stress 0.081757 
    ## Run 32 stress 0.08360269 
    ## Run 33 stress 0.07731538 
    ## Run 34 stress 0.08375705 
    ## Run 35 stress 0.08314409 
    ## Run 36 stress 0.0794898 
    ## Run 37 stress 0.07871737 
    ## Run 38 stress 0.08375705 
    ## Run 39 stress 0.08421021 
    ## Run 40 stress 0.08097776 
    ## Run 41 stress 0.0824961 
    ## Run 42 stress 0.08265425 
    ## Run 43 stress 0.0766915 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000102889  max resid 0.0002086796 
    ## ... Similar to previous best
    ## Run 44 stress 0.08433701 
    ## Run 45 stress 0.08003186 
    ## Run 46 stress 0.08240494 
    ## Run 47 stress 0.07951483 
    ## Run 48 stress 0.08317343 
    ## Run 49 stress 0.08360266 
    ## Run 50 stress 0.08182125 
    ## Run 51 stress 0.08355807 
    ## Run 52 stress 0.08182128 
    ## Run 53 stress 0.0810491 
    ## Run 54 stress 0.08265366 
    ## Run 55 stress 0.08104919 
    ## Run 56 stress 0.08582973 
    ## Run 57 stress 0.0801294 
    ## Run 58 stress 0.07669164 
    ## ... Procrustes: rmse 0.0002377993  max resid 0.0005301038 
    ## ... Similar to previous best
    ## Run 59 stress 0.08485176 
    ## Run 60 stress 0.07960687 
    ## Run 61 stress 0.0831735 
    ## Run 62 stress 0.07960657 
    ## Run 63 stress 0.07731519 
    ## Run 64 stress 0.0794336 
    ## Run 65 stress 0.07871744 
    ## Run 66 stress 0.07871726 
    ## Run 67 stress 0.08485134 
    ## Run 68 stress 0.07731523 
    ## Run 69 stress 0.08379905 
    ## Run 70 stress 0.0794339 
    ## Run 71 stress 0.08097754 
    ## Run 72 stress 0.08464803 
    ## Run 73 stress 0.08704886 
    ## Run 74 stress 0.08097771 
    ## Run 75 stress 0.08175698 
    ## Run 76 stress 0.08360257 
    ## Run 77 stress 0.0830555 
    ## Run 78 stress 0.08530158 
    ## Run 79 stress 0.08379866 
    ## Run 80 stress 0.0831731 
    ## Run 81 stress 0.08175711 
    ## Run 82 stress 0.0826533 
    ## Run 83 stress 0.08283538 
    ## Run 84 stress 0.07868394 
    ## Run 85 stress 0.07669173 
    ## ... Procrustes: rmse 0.0002822272  max resid 0.0004978482 
    ## ... Similar to previous best
    ## Run 86 stress 0.08104886 
    ## Run 87 stress 0.08527016 
    ## Run 88 stress 0.08175719 
    ## Run 89 stress 0.08968625 
    ## Run 90 stress 0.08217343 
    ## Run 91 stress 0.07871735 
    ## Run 92 stress 0.07960663 
    ## Run 93 stress 0.08360277 
    ## Run 94 stress 0.08957495 
    ## Run 95 stress 0.08175712 
    ## Run 96 stress 0.08835165 
    ## Run 97 stress 0.0801295 
    ## Run 98 stress 0.08012986 
    ## Run 99 stress 0.08418572 
    ## Run 100 stress 0.08305496 
    ## Run 101 stress 0.0766917 
    ## ... Procrustes: rmse 0.0002742813  max resid 0.0004811274 
    ## ... Similar to previous best
    ## Run 102 stress 0.08144838 
    ## Run 103 stress 0.08360261 
    ## Run 104 stress 0.08265372 
    ## Run 105 stress 0.08355838 
    ## Run 106 stress 0.08104884 
    ## Run 107 stress 0.08097769 
    ## Run 108 stress 0.07871726 
    ## Run 109 stress 0.08453202 
    ## Run 110 stress 0.08097747 
    ## Run 111 stress 0.08175695 
    ## Run 112 stress 0.08013388 
    ## Run 113 stress 0.07960688 
    ## Run 114 stress 0.09091579 
    ## Run 115 stress 0.07949009 
    ## Run 116 stress 0.08175697 
    ## Run 117 stress 0.07669157 
    ## ... Procrustes: rmse 0.0001720628  max resid 0.0003406619 
    ## ... Similar to previous best
    ## Run 118 stress 0.0850036 
    ## Run 119 stress 0.08527016 
    ## Run 120 stress 0.08485152 
    ## Run 121 stress 0.08422698 
    ## Run 122 stress 0.07669161 
    ## ... Procrustes: rmse 0.0002089533  max resid 0.0003670586 
    ## ... Similar to previous best
    ## Run 123 stress 0.081449 
    ## Run 124 stress 0.08421002 
    ## Run 125 stress 0.08104887 
    ## Run 126 stress 0.08337151 
    ## Run 127 stress 0.08387162 
    ## Run 128 stress 0.08013552 
    ## Run 129 stress 0.07669164 
    ## ... Procrustes: rmse 0.0002281958  max resid 0.0004004287 
    ## ... Similar to previous best
    ## Run 130 stress 0.08471585 
    ## Run 131 stress 0.08684103 
    ## Run 132 stress 0.08500322 
    ## Run 133 stress 0.0786835 
    ## Run 134 stress 0.2801452 
    ## Run 135 stress 0.08672705 
    ## Run 136 stress 0.08485159 
    ## Run 137 stress 0.08241485 
    ## Run 138 stress 0.09224576 
    ## Run 139 stress 0.08835158 
    ## Run 140 stress 0.0794901 
    ## Run 141 stress 0.08337165 
    ## Run 142 stress 0.0817569 
    ## Run 143 stress 0.08471582 
    ## Run 144 stress 0.08486913 
    ## Run 145 stress 0.0800319 
    ## Run 146 stress 0.08175704 
    ## Run 147 stress 0.08104917 
    ## Run 148 stress 0.08329776 
    ## Run 149 stress 0.08835157 
    ## Run 150 stress 0.08375737 
    ## Run 151 stress 0.08842068 
    ## Run 152 stress 0.08433687 
    ## Run 153 stress 0.08217348 
    ## Run 154 stress 0.08957554 
    ## Run 155 stress 0.08317385 
    ## Run 156 stress 0.0831734 
    ## Run 157 stress 0.08283558 
    ## Run 158 stress 0.08144948 
    ## Run 159 stress 0.07669163 
    ## ... Procrustes: rmse 0.0002088407  max resid 0.0003649692 
    ## ... Similar to previous best
    ## Run 160 stress 0.08337135 
    ## Run 161 stress 0.0846478 
    ## Run 162 stress 0.08370298 
    ## Run 163 stress 0.07943366 
    ## Run 164 stress 0.08243398 
    ## Run 165 stress 0.07871744 
    ## Run 166 stress 0.08665276 
    ## Run 167 stress 0.08003198 
    ## Run 168 stress 0.08013927 
    ## Run 169 stress 0.0796068 
    ## Run 170 stress 0.08553087 
    ## Run 171 stress 0.07960685 
    ## Run 172 stress 0.08012966 
    ## Run 173 stress 0.08265472 
    ## Run 174 stress 0.08527015 
    ## Run 175 stress 0.08443327 
    ## Run 176 stress 0.08314435 
    ## Run 177 stress 0.08182108 
    ## Run 178 stress 0.0794338 
    ## Run 179 stress 0.08218533 
    ## Run 180 stress 0.08265427 
    ## Run 181 stress 0.07669168 
    ## ... Procrustes: rmse 0.0002465528  max resid 0.0004240404 
    ## ... Similar to previous best
    ## Run 182 stress 0.08835154 
    ## Run 183 stress 0.0786839 
    ## Run 184 stress 0.0827496 
    ## Run 185 stress 0.0818212 
    ## Run 186 stress 0.08355817 
    ## Run 187 stress 0.08337158 
    ## Run 188 stress 0.08243394 
    ## Run 189 stress 0.08470959 
    ## Run 190 stress 0.08476039 
    ## Run 191 stress 0.09186061 
    ## Run 192 stress 0.08003196 
    ## Run 193 stress 0.08265464 
    ## Run 194 stress 0.0824148 
    ## Run 195 stress 0.08249604 
    ## Run 196 stress 0.0843371 
    ## Run 197 stress 0.07943326 
    ## Run 198 stress 0.08175692 
    ## Run 199 stress 0.0806393 
    ## Run 200 stress 0.08283547 
    ## Run 201 stress 0.08243374 
    ## Run 202 stress 0.08477454 
    ## Run 203 stress 0.08175701 
    ## Run 204 stress 0.07960694 
    ## Run 205 stress 0.08433699 
    ## Run 206 stress 0.07948961 
    ## Run 207 stress 0.08003206 
    ## Run 208 stress 0.08485128 
    ## Run 209 stress 0.08104918 
    ## Run 210 stress 0.0883514 
    ## Run 211 stress 0.08314437 
    ## Run 212 stress 0.07669152 
    ## ... Procrustes: rmse 9.146438e-05  max resid 0.0001847884 
    ## ... Similar to previous best
    ## Run 213 stress 0.08957573 
    ## Run 214 stress 0.0801294 
    ## Run 215 stress 0.07868385 
    ## Run 216 stress 0.08355839 
    ## Run 217 stress 0.0818939 
    ## Run 218 stress 0.08250112 
    ## Run 219 stress 0.07868344 
    ## Run 220 stress 0.08175703 
    ## Run 221 stress 0.08433695 
    ## Run 222 stress 0.08329754 
    ## Run 223 stress 0.08682716 
    ## Run 224 stress 0.08233234 
    ## Run 225 stress 0.08189396 
    ## Run 226 stress 0.07731516 
    ## Run 227 stress 0.08182108 
    ## Run 228 stress 0.08189378 
    ## Run 229 stress 0.08433686 
    ## Run 230 stress 0.08175704 
    ## Run 231 stress 0.07669153 
    ## ... Procrustes: rmse 8.216081e-05  max resid 0.0001807846 
    ## ... Similar to previous best
    ## Run 232 stress 0.09196181 
    ## Run 233 stress 0.08957492 
    ## Run 234 stress 0.08252139 
    ## Run 235 stress 0.08243386 
    ## Run 236 stress 0.08003221 
    ## Run 237 stress 0.08375742 
    ## Run 238 stress 0.08317311 
    ## Run 239 stress 0.08265463 
    ## Run 240 stress 0.08337156 
    ## Run 241 stress 0.07949034 
    ## Run 242 stress 0.08250164 
    ## Run 243 stress 0.08252117 
    ## Run 244 stress 0.08527016 
    ## Run 245 stress 0.08835155 
    ## Run 246 stress 0.08175715 
    ## Run 247 stress 0.08249624 
    ## Run 248 stress 0.08012961 
    ## Run 249 stress 0.0830546 
    ## Run 250 stress 0.259577 
    ## Run 251 stress 0.08182124 
    ## Run 252 stress 0.08355803 
    ## Run 253 stress 0.08485148 
    ## Run 254 stress 0.08252116 
    ## Run 255 stress 0.08175711 
    ## Run 256 stress 0.0852703 
    ## Run 257 stress 0.08175694 
    ## Run 258 stress 0.08305348 
    ## Run 259 stress 0.08305308 
    ## Run 260 stress 0.08189395 
    ## Run 261 stress 0.08317291 
    ## Run 262 stress 0.08104887 
    ## Run 263 stress 0.08433692 
    ## Run 264 stress 0.08265423 
    ## Run 265 stress 0.0826534 
    ## Run 266 stress 0.0818212 
    ## Run 267 stress 0.08362769 
    ## Run 268 stress 0.08104897 
    ## Run 269 stress 0.08835177 
    ## Run 270 stress 0.08422695 
    ## Run 271 stress 0.08104899 
    ## Run 272 stress 0.08003208 
    ## Run 273 stress 0.081049 
    ## Run 274 stress 0.08252146 
    ## Run 275 stress 0.08687395 
    ## Run 276 stress 0.08003196 
    ## Run 277 stress 0.08265353 
    ## Run 278 stress 0.08957564 
    ## Run 279 stress 0.07960688 
    ## Run 280 stress 0.08175691 
    ## Run 281 stress 0.08104921 
    ## Run 282 stress 0.08337107 
    ## Run 283 stress 0.0817569 
    ## Run 284 stress 0.0824179 
    ## Run 285 stress 0.07669158 
    ## ... Procrustes: rmse 0.0001844393  max resid 0.0003662974 
    ## ... Similar to previous best
    ## Run 286 stress 0.07871743 
    ## Run 287 stress 0.07951472 
    ## Run 288 stress 0.08175694 
    ## Run 289 stress 0.08421013 
    ## Run 290 stress 0.08104927 
    ## Run 291 stress 0.08012981 
    ## Run 292 stress 0.08182119 
    ## Run 293 stress 0.08182124 
    ## Run 294 stress 0.0794337 
    ## Run 295 stress 0.08063952 
    ## Run 296 stress 0.08012956 
    ## Run 297 stress 0.07868325 
    ## Run 298 stress 0.08317342 
    ## Run 299 stress 0.08379853 
    ## Run 300 stress 0.08283553 
    ## Run 301 stress 0.07871765 
    ## Run 302 stress 0.07731512 
    ## Run 303 stress 0.08317328 
    ## Run 304 stress 0.08175694 
    ## Run 305 stress 0.0928285 
    ## Run 306 stress 0.08243377 
    ## Run 307 stress 0.08835182 
    ## Run 308 stress 0.08097746 
    ## Run 309 stress 0.08175688 
    ## Run 310 stress 0.08835178 
    ## Run 311 stress 0.08104885 
    ## Run 312 stress 0.0842101 
    ## Run 313 stress 0.08433697 
    ## Run 314 stress 0.08144953 
    ## Run 315 stress 0.08252131 
    ## Run 316 stress 0.08003192 
    ## Run 317 stress 0.08553085 
    ## Run 318 stress 0.07731544 
    ## Run 319 stress 0.08684096 
    ## Run 320 stress 0.08314436 
    ## Run 321 stress 0.08355842 
    ## Run 322 stress 0.0800322 
    ## Run 323 stress 0.08553066 
    ## Run 324 stress 0.08317342 
    ## Run 325 stress 0.08175689 
    ## Run 326 stress 0.08527038 
    ## Run 327 stress 0.08175697 
    ## Run 328 stress 0.08175688 
    ## Run 329 stress 0.08477411 
    ## Run 330 stress 0.08097757 
    ## Run 331 stress 0.07871754 
    ## Run 332 stress 0.08175713 
    ## Run 333 stress 0.0836026 
    ## Run 334 stress 0.07948954 
    ## Run 335 stress 0.07669159 
    ## ... Procrustes: rmse 0.0001755354  max resid 0.0002997309 
    ## ... Similar to previous best
    ## Run 336 stress 0.07948965 
    ## Run 337 stress 0.08501465 
    ## Run 338 stress 0.08104898 
    ## Run 339 stress 0.08104885 
    ## Run 340 stress 0.08929085 
    ## Run 341 stress 0.08433704 
    ## Run 342 stress 0.0796068 
    ## Run 343 stress 0.08433692 
    ## Run 344 stress 0.08370294 
    ## Run 345 stress 0.081049 
    ## Run 346 stress 0.0796068 
    ## Run 347 stress 0.08687337 
    ## Run 348 stress 0.07871738 
    ## Run 349 stress 0.08233242 
    ## Run 350 stress 0.08097768 
    ## Run 351 stress 0.0836026 
    ## Run 352 stress 0.08233221 
    ## Run 353 stress 0.08375714 
    ## Run 354 stress 0.08182127 
    ## Run 355 stress 0.08360277 
    ## Run 356 stress 0.08337139 
    ## Run 357 stress 0.08305412 
    ## Run 358 stress 0.07669164 
    ## ... Procrustes: rmse 0.0002346236  max resid 0.0004075614 
    ## ... Similar to previous best
    ## Run 359 stress 0.07868368 
    ## Run 360 stress 0.08443332 
    ## Run 361 stress 0.08252117 
    ## Run 362 stress 0.07868353 
    ## Run 363 stress 0.08317342 
    ## Run 364 stress 0.08265423 
    ## Run 365 stress 0.08013781 
    ## Run 366 stress 0.08265445 
    ## Run 367 stress 0.07943373 
    ## Run 368 stress 0.07731523 
    ## Run 369 stress 0.08374681 
    ## Run 370 stress 0.07868391 
    ## Run 371 stress 0.08012949 
    ## Run 372 stress 0.07669157 
    ## ... Procrustes: rmse 0.0001615923  max resid 0.0002836833 
    ## ... Similar to previous best
    ## Run 373 stress 0.08182111 
    ## Run 374 stress 0.08003181 
    ## Run 375 stress 0.08175691 
    ## Run 376 stress 0.08470925 
    ## Run 377 stress 0.08104886 
    ## Run 378 stress 0.08433708 
    ## Run 379 stress 0.07731528 
    ## Run 380 stress 0.08957527 
    ## Run 381 stress 0.07868331 
    ## Run 382 stress 0.08527047 
    ## Run 383 stress 0.08144711 
    ## Run 384 stress 0.07951482 
    ## Run 385 stress 0.0825213 
    ## Run 386 stress 0.0817569 
    ## Run 387 stress 0.08144928 
    ## Run 388 stress 0.08252141 
    ## Run 389 stress 0.08898769 
    ## Run 390 stress 0.07871729 
    ## Run 391 stress 0.08443251 
    ## Run 392 stress 0.08553085 
    ## Run 393 stress 0.08317375 
    ## Run 394 stress 0.0851084 
    ## Run 395 stress 0.0889878 
    ## Run 396 stress 0.08616026 
    ## Run 397 stress 0.07960675 
    ## Run 398 stress 0.08218541 
    ## Run 399 stress 0.08175692 
    ## Run 400 stress 0.08104891 
    ## Run 401 stress 0.08012991 
    ## Run 402 stress 0.08433697 
    ## Run 403 stress 0.08443323 
    ## Run 404 stress 0.08842065 
    ## Run 405 stress 0.08957607 
    ## Run 406 stress 0.0817571 
    ## Run 407 stress 0.08189421 
    ## Run 408 stress 0.08104901 
    ## Run 409 stress 0.09103935 
    ## Run 410 stress 0.08104907 
    ## Run 411 stress 0.08615993 
    ## Run 412 stress 0.08355829 
    ## Run 413 stress 0.08175689 
    ## Run 414 stress 0.08485162 
    ## Run 415 stress 0.0846479 
    ## Run 416 stress 0.07949012 
    ## Run 417 stress 0.08012969 
    ## Run 418 stress 0.08189394 
    ## Run 419 stress 0.08175692 
    ## Run 420 stress 0.0766915 
    ## ... New best solution
    ## ... Procrustes: rmse 4.891311e-05  max resid 9.145644e-05 
    ## ... Similar to previous best
    ## Run 421 stress 0.08012938 
    ## Run 422 stress 0.09034799 
    ## Run 423 stress 0.08175698 
    ## Run 424 stress 0.08527026 
    ## Run 425 stress 0.08375708 
    ## Run 426 stress 0.08104909 
    ## Run 427 stress 0.08175712 
    ## Run 428 stress 0.07669153 
    ## ... Procrustes: rmse 0.0001432061  max resid 0.000246962 
    ## ... Similar to previous best
    ## Run 429 stress 0.09012077 
    ## Run 430 stress 0.07948971 
    ## Run 431 stress 0.07951486 
    ## Run 432 stress 0.08265433 
    ## Run 433 stress 0.08241507 
    ## Run 434 stress 0.0824149 
    ## Run 435 stress 0.08317295 
    ## Run 436 stress 0.08418572 
    ## Run 437 stress 0.08012969 
    ## Run 438 stress 0.08097765 
    ## Run 439 stress 0.08470918 
    ## Run 440 stress 0.08337145 
    ## Run 441 stress 0.08217351 
    ## Run 442 stress 0.08421029 
    ## Run 443 stress 0.08265435 
    ## Run 444 stress 0.07669152 
    ## ... Procrustes: rmse 6.492641e-05  max resid 0.0001123396 
    ## ... Similar to previous best
    ## Run 445 stress 0.08252129 
    ## Run 446 stress 0.08182122 
    ## Run 447 stress 0.07868355 
    ## Run 448 stress 0.07943325 
    ## Run 449 stress 0.08433685 
    ## Run 450 stress 0.08265367 
    ## Run 451 stress 0.0809776 
    ## Run 452 stress 0.07669157 
    ## ... Procrustes: rmse 0.0001348677  max resid 0.000235758 
    ## ... Similar to previous best
    ## Run 453 stress 0.08175698 
    ## Run 454 stress 0.08375706 
    ## Run 455 stress 0.08339772 
    ## Run 456 stress 0.08175712 
    ## Run 457 stress 0.07868338 
    ## Run 458 stress 0.07948947 
    ## Run 459 stress 0.0833711 
    ## Run 460 stress 0.0766915 
    ## ... Procrustes: rmse 1.852153e-05  max resid 3.995073e-05 
    ## ... Similar to previous best
    ## Run 461 stress 0.08305405 
    ## Run 462 stress 0.08249622 
    ## Run 463 stress 0.08957601 
    ## Run 464 stress 0.08104918 
    ## Run 465 stress 0.08182134 
    ## Run 466 stress 0.09022445 
    ## Run 467 stress 0.07949025 
    ## Run 468 stress 0.08337209 
    ## Run 469 stress 0.08252118 
    ## Run 470 stress 0.08182107 
    ## Run 471 stress 0.07871729 
    ## Run 472 stress 0.08175689 
    ## Run 473 stress 0.08218526 
    ## Run 474 stress 0.08104907 
    ## Run 475 stress 0.08362759 
    ## Run 476 stress 0.08464785 
    ## Run 477 stress 0.081757 
    ## Run 478 stress 0.08337182 
    ## Run 479 stress 0.086841 
    ## Run 480 stress 0.08305347 
    ## Run 481 stress 0.08835171 
    ## Run 482 stress 0.08888012 
    ## Run 483 stress 0.08104913 
    ## Run 484 stress 0.08012961 
    ## Run 485 stress 0.08175692 
    ## Run 486 stress 0.08471589 
    ## Run 487 stress 0.08443463 
    ## Run 488 stress 0.08355826 
    ## Run 489 stress 0.08265406 
    ## Run 490 stress 0.08104894 
    ## Run 491 stress 0.09091602 
    ## Run 492 stress 0.08507061 
    ## Run 493 stress 0.08104901 
    ## Run 494 stress 0.08249616 
    ## Run 495 stress 0.08337128 
    ## Run 496 stress 0.07943348 
    ## Run 497 stress 0.08355829 
    ## Run 498 stress 0.08265436 
    ## Run 499 stress 0.08265388 
    ## Run 500 stress 0.0843369 
    ## *** Best solution repeated 5 times

``` r
# Mixed and stratified lakes
SD_beta_env_MS_NMDS <- metaMDS(SD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09030397 
    ## Run 2 stress 0.08503656 
    ## ... Procrustes: rmse 0.0002475252  max resid 0.0005238487 
    ## ... Similar to previous best
    ## Run 3 stress 0.08773467 
    ## Run 4 stress 0.09374299 
    ## Run 5 stress 0.09030394 
    ## Run 6 stress 0.08503474 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001584528  max resid 0.004067439 
    ## ... Similar to previous best
    ## Run 7 stress 0.09760822 
    ## Run 8 stress 0.09445476 
    ## Run 9 stress 0.09030403 
    ## Run 10 stress 0.0946592 
    ## Run 11 stress 0.0937431 
    ## Run 12 stress 0.09416133 
    ## Run 13 stress 0.09286088 
    ## Run 14 stress 0.09337297 
    ## Run 15 stress 0.09400524 
    ## Run 16 stress 0.09669834 
    ## Run 17 stress 0.08503476 
    ## ... Procrustes: rmse 1.680714e-05  max resid 4.444495e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.09030394 
    ## Run 19 stress 0.09374214 
    ## Run 20 stress 0.09465914 
    ## Run 21 stress 0.09145333 
    ## Run 22 stress 0.09030398 
    ## Run 23 stress 0.09268326 
    ## Run 24 stress 0.08440252 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01417374  max resid 0.04293991 
    ## Run 25 stress 0.09407967 
    ## Run 26 stress 0.09535512 
    ## Run 27 stress 0.09145324 
    ## Run 28 stress 0.09407966 
    ## Run 29 stress 0.2776025 
    ## Run 30 stress 0.2793156 
    ## Run 31 stress 0.09374332 
    ## Run 32 stress 0.09535478 
    ## Run 33 stress 0.09168951 
    ## Run 34 stress 0.08973862 
    ## Run 35 stress 0.1038033 
    ## Run 36 stress 0.08503498 
    ## Run 37 stress 0.0903041 
    ## Run 38 stress 0.09268344 
    ## Run 39 stress 0.103804 
    ## Run 40 stress 0.09286092 
    ## Run 41 stress 0.09407375 
    ## Run 42 stress 0.08440267 
    ## ... Procrustes: rmse 0.0001560514  max resid 0.0002926043 
    ## ... Similar to previous best
    ## Run 43 stress 0.09159089 
    ## Run 44 stress 0.08503514 
    ## Run 45 stress 0.09465915 
    ## Run 46 stress 0.08440266 
    ## ... Procrustes: rmse 0.0001503278  max resid 0.000258281 
    ## ... Similar to previous best
    ## Run 47 stress 0.09159095 
    ## Run 48 stress 0.0850346 
    ## Run 49 stress 0.08973868 
    ## Run 50 stress 0.09286086 
    ## Run 51 stress 0.09030403 
    ## Run 52 stress 0.100461 
    ## Run 53 stress 0.08773474 
    ## Run 54 stress 0.09403398 
    ## Run 55 stress 0.09337238 
    ## Run 56 stress 0.08440279 
    ## ... Procrustes: rmse 0.0002429418  max resid 0.0004548439 
    ## ... Similar to previous best
    ## Run 57 stress 0.08973877 
    ## Run 58 stress 0.0850349 
    ## Run 59 stress 0.09407969 
    ## Run 60 stress 0.08440269 
    ## ... Procrustes: rmse 0.0001751389  max resid 0.0003069605 
    ## ... Similar to previous best
    ## Run 61 stress 0.09030396 
    ## Run 62 stress 0.09374244 
    ## Run 63 stress 0.09030396 
    ## Run 64 stress 0.2767877 
    ## Run 65 stress 0.09308983 
    ## Run 66 stress 0.09268331 
    ## Run 67 stress 0.09712949 
    ## Run 68 stress 0.0940797 
    ## Run 69 stress 0.09337272 
    ## Run 70 stress 0.0940797 
    ## Run 71 stress 0.2510663 
    ## Run 72 stress 0.09610765 
    ## Run 73 stress 0.09380586 
    ## Run 74 stress 0.09539206 
    ## Run 75 stress 0.09168935 
    ## Run 76 stress 0.09337256 
    ## Run 77 stress 0.08773466 
    ## Run 78 stress 0.0897389 
    ## Run 79 stress 0.09286083 
    ## Run 80 stress 0.09030398 
    ## Run 81 stress 0.1084716 
    ## Run 82 stress 0.09374365 
    ## Run 83 stress 0.09308951 
    ## Run 84 stress 0.0916895 
    ## Run 85 stress 0.090304 
    ## Run 86 stress 0.09030399 
    ## Run 87 stress 0.09403444 
    ## Run 88 stress 0.09268318 
    ## Run 89 stress 0.0916893 
    ## Run 90 stress 0.09145329 
    ## Run 91 stress 0.09417813 
    ## Run 92 stress 0.08973883 
    ## Run 93 stress 0.09407967 
    ## Run 94 stress 0.09539192 
    ## Run 95 stress 0.09268397 
    ## Run 96 stress 0.09465911 
    ## Run 97 stress 0.09416145 
    ## Run 98 stress 0.09030403 
    ## Run 99 stress 0.09337306 
    ## Run 100 stress 0.09337265 
    ## Run 101 stress 0.08503489 
    ## Run 102 stress 0.08503459 
    ## Run 103 stress 0.08503508 
    ## Run 104 stress 0.09337232 
    ## Run 105 stress 0.08503465 
    ## Run 106 stress 0.09268331 
    ## Run 107 stress 0.09380587 
    ## Run 108 stress 0.09465905 
    ## Run 109 stress 0.09030394 
    ## Run 110 stress 0.09268328 
    ## Run 111 stress 0.08503482 
    ## Run 112 stress 0.08773484 
    ## Run 113 stress 0.09268359 
    ## Run 114 stress 0.09030398 
    ## Run 115 stress 0.09030401 
    ## Run 116 stress 0.09465904 
    ## Run 117 stress 0.09030397 
    ## Run 118 stress 0.08503498 
    ## Run 119 stress 0.08503556 
    ## Run 120 stress 0.1053419 
    ## Run 121 stress 0.09030399 
    ## Run 122 stress 0.09286087 
    ## Run 123 stress 0.08503655 
    ## Run 124 stress 0.09407994 
    ## Run 125 stress 0.09374223 
    ## Run 126 stress 0.09145325 
    ## Run 127 stress 0.09337285 
    ## Run 128 stress 0.09159084 
    ## Run 129 stress 0.2456114 
    ## Run 130 stress 0.0941776 
    ## Run 131 stress 0.08503511 
    ## Run 132 stress 0.09286083 
    ## Run 133 stress 0.09268372 
    ## Run 134 stress 0.09374258 
    ## Run 135 stress 0.09407999 
    ## Run 136 stress 0.09286089 
    ## Run 137 stress 0.1103049 
    ## Run 138 stress 0.09286095 
    ## Run 139 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 4.188014e-05  max resid 8.738082e-05 
    ## ... Similar to previous best
    ## Run 140 stress 0.09337286 
    ## Run 141 stress 0.09407995 
    ## Run 142 stress 0.09337254 
    ## Run 143 stress 0.09374308 
    ## Run 144 stress 0.09969969 
    ## Run 145 stress 0.09407971 
    ## Run 146 stress 0.09308958 
    ## Run 147 stress 0.09464421 
    ## Run 148 stress 0.09030408 
    ## Run 149 stress 0.08773482 
    ## Run 150 stress 0.09374248 
    ## Run 151 stress 0.089739 
    ## Run 152 stress 0.09337265 
    ## Run 153 stress 0.09465904 
    ## Run 154 stress 0.0850351 
    ## Run 155 stress 0.09407965 
    ## Run 156 stress 0.09337231 
    ## Run 157 stress 0.08503525 
    ## Run 158 stress 0.08973886 
    ## Run 159 stress 0.1038048 
    ## Run 160 stress 0.08773465 
    ## Run 161 stress 0.09407998 
    ## Run 162 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001273356  max resid 0.0002518185 
    ## ... Similar to previous best
    ## Run 163 stress 0.08973871 
    ## Run 164 stress 0.09407989 
    ## Run 165 stress 0.1038076 
    ## Run 166 stress 0.09286088 
    ## Run 167 stress 0.08973862 
    ## Run 168 stress 0.09321717 
    ## Run 169 stress 0.08503484 
    ## Run 170 stress 0.09159103 
    ## Run 171 stress 0.09159083 
    ## Run 172 stress 0.09030394 
    ## Run 173 stress 0.09465905 
    ## Run 174 stress 0.1026216 
    ## Run 175 stress 0.09380584 
    ## Run 176 stress 0.09407993 
    ## Run 177 stress 0.0976083 
    ## Run 178 stress 0.09374278 
    ## Run 179 stress 0.09030397 
    ## Run 180 stress 0.09416132 
    ## Run 181 stress 0.09408017 
    ## Run 182 stress 0.09030411 
    ## Run 183 stress 0.09381848 
    ## Run 184 stress 0.09590971 
    ## Run 185 stress 0.08773467 
    ## Run 186 stress 0.1038045 
    ## Run 187 stress 0.09268387 
    ## Run 188 stress 0.08440256 
    ## ... Procrustes: rmse 9.57535e-05  max resid 0.0001997057 
    ## ... Similar to previous best
    ## Run 189 stress 0.08440253 
    ## ... Procrustes: rmse 8.291417e-05  max resid 0.0002688933 
    ## ... Similar to previous best
    ## Run 190 stress 0.09030397 
    ## Run 191 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001550045  max resid 0.0003235093 
    ## ... Similar to previous best
    ## Run 192 stress 0.09416159 
    ## Run 193 stress 0.0940055 
    ## Run 194 stress 0.08973882 
    ## Run 195 stress 0.09464454 
    ## Run 196 stress 0.0926836 
    ## Run 197 stress 0.09407989 
    ## Run 198 stress 0.09145335 
    ## Run 199 stress 0.09380579 
    ## Run 200 stress 0.09030407 
    ## Run 201 stress 0.092684 
    ## Run 202 stress 0.09145329 
    ## Run 203 stress 0.09760826 
    ## Run 204 stress 0.2914152 
    ## Run 205 stress 0.09145328 
    ## Run 206 stress 0.09590955 
    ## Run 207 stress 0.08440269 
    ## ... Procrustes: rmse 0.0002952624  max resid 0.0005889052 
    ## ... Similar to previous best
    ## Run 208 stress 0.09030406 
    ## Run 209 stress 0.09380589 
    ## Run 210 stress 0.105285 
    ## Run 211 stress 0.08503501 
    ## Run 212 stress 0.0941614 
    ## Run 213 stress 0.09407964 
    ## Run 214 stress 0.3411713 
    ## Run 215 stress 0.08773467 
    ## Run 216 stress 0.09030395 
    ## Run 217 stress 0.0976083 
    ## Run 218 stress 0.08440252 
    ## ... Procrustes: rmse 4.456407e-05  max resid 8.446398e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.09286087 
    ## Run 220 stress 0.0959065 
    ## Run 221 stress 0.08503483 
    ## Run 222 stress 0.09030411 
    ## Run 223 stress 0.0877347 
    ## Run 224 stress 0.1052854 
    ## Run 225 stress 0.09159095 
    ## Run 226 stress 0.09168949 
    ## Run 227 stress 0.09407987 
    ## Run 228 stress 0.09145327 
    ## Run 229 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001794722  max resid 0.0003428451 
    ## ... Similar to previous best
    ## Run 230 stress 0.09465906 
    ## Run 231 stress 0.09408006 
    ## Run 232 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001016826  max resid 0.0002147035 
    ## ... Similar to previous best
    ## Run 233 stress 0.09268354 
    ## Run 234 stress 0.08503493 
    ## Run 235 stress 0.09030408 
    ## Run 236 stress 0.08503558 
    ## Run 237 stress 0.09465905 
    ## Run 238 stress 0.08503548 
    ## Run 239 stress 0.09380581 
    ## Run 240 stress 0.2574698 
    ## Run 241 stress 0.09337265 
    ## Run 242 stress 0.09408001 
    ## Run 243 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001696601  max resid 0.0003307451 
    ## ... Similar to previous best
    ## Run 244 stress 0.0996997 
    ## Run 245 stress 0.09159084 
    ## Run 246 stress 0.09268323 
    ## Run 247 stress 0.09374373 
    ## Run 248 stress 0.09159106 
    ## Run 249 stress 0.09721202 
    ## Run 250 stress 0.09447324 
    ## Run 251 stress 0.09417786 
    ## Run 252 stress 0.09417776 
    ## Run 253 stress 0.1052848 
    ## Run 254 stress 0.09308955 
    ## Run 255 stress 0.09539206 
    ## Run 256 stress 0.08440271 
    ## ... Procrustes: rmse 0.000230696  max resid 0.0004376471 
    ## ... Similar to previous best
    ## Run 257 stress 0.09145322 
    ## Run 258 stress 0.09465905 
    ## Run 259 stress 0.08440253 
    ## ... Procrustes: rmse 4.834779e-05  max resid 0.0001085914 
    ## ... Similar to previous best
    ## Run 260 stress 0.09168928 
    ## Run 261 stress 0.08503572 
    ## Run 262 stress 0.09407964 
    ## Run 263 stress 0.1054547 
    ## Run 264 stress 0.1052849 
    ## Run 265 stress 0.09464494 
    ## Run 266 stress 0.09407983 
    ## Run 267 stress 0.09168942 
    ## Run 268 stress 0.09168943 
    ## Run 269 stress 0.09030399 
    ## Run 270 stress 0.3244864 
    ## Run 271 stress 0.0897389 
    ## Run 272 stress 0.09321382 
    ## Run 273 stress 0.1026216 
    ## Run 274 stress 0.09030399 
    ## Run 275 stress 0.09168946 
    ## Run 276 stress 0.09337232 
    ## Run 277 stress 0.08503545 
    ## Run 278 stress 0.09407974 
    ## Run 279 stress 0.09407145 
    ## Run 280 stress 0.09308969 
    ## Run 281 stress 0.09030395 
    ## Run 282 stress 0.09407985 
    ## Run 283 stress 0.09030398 
    ## Run 284 stress 0.09308949 
    ## Run 285 stress 0.09407978 
    ## Run 286 stress 0.09030394 
    ## Run 287 stress 0.09374232 
    ## Run 288 stress 0.09030394 
    ## Run 289 stress 0.09447305 
    ## Run 290 stress 0.08973865 
    ## Run 291 stress 0.08503496 
    ## Run 292 stress 0.09407966 
    ## Run 293 stress 0.09286084 
    ## Run 294 stress 0.09030401 
    ## Run 295 stress 0.08773471 
    ## Run 296 stress 0.09380577 
    ## Run 297 stress 0.09374355 
    ## Run 298 stress 0.09159087 
    ## Run 299 stress 0.09268394 
    ## Run 300 stress 0.09159094 
    ## Run 301 stress 0.0844028 
    ## ... Procrustes: rmse 0.0002421809  max resid 0.0004207618 
    ## ... Similar to previous best
    ## Run 302 stress 0.09159091 
    ## Run 303 stress 0.09159099 
    ## Run 304 stress 0.09590849 
    ## Run 305 stress 0.08973867 
    ## Run 306 stress 0.09030396 
    ## Run 307 stress 0.09030394 
    ## Run 308 stress 0.08773466 
    ## Run 309 stress 0.1052851 
    ## Run 310 stress 0.08503455 
    ## Run 311 stress 0.09465907 
    ## Run 312 stress 0.09445489 
    ## Run 313 stress 0.1038101 
    ## Run 314 stress 0.09400534 
    ## Run 315 stress 0.09159088 
    ## Run 316 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 3.623603e-05  max resid 6.785446e-05 
    ## ... Similar to previous best
    ## Run 317 stress 0.09407494 
    ## Run 318 stress 0.0940799 
    ## Run 319 stress 0.09969974 
    ## Run 320 stress 0.09465913 
    ## Run 321 stress 0.09286089 
    ## Run 322 stress 0.09286083 
    ## Run 323 stress 0.09969961 
    ## Run 324 stress 0.09030407 
    ## Run 325 stress 0.09760827 
    ## Run 326 stress 0.09464462 
    ## Run 327 stress 0.08773468 
    ## Run 328 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001133709  max resid 0.0002064253 
    ## ... Similar to previous best
    ## Run 329 stress 0.09168928 
    ## Run 330 stress 0.08503509 
    ## Run 331 stress 0.09159092 
    ## Run 332 stress 0.08440262 
    ## ... Procrustes: rmse 0.0002040463  max resid 0.000342566 
    ## ... Similar to previous best
    ## Run 333 stress 0.08503462 
    ## Run 334 stress 0.09374235 
    ## Run 335 stress 0.09145325 
    ## Run 336 stress 0.09268352 
    ## Run 337 stress 0.09374334 
    ## Run 338 stress 0.09030405 
    ## Run 339 stress 0.09030397 
    ## Run 340 stress 0.09590917 
    ## Run 341 stress 0.09030396 
    ## Run 342 stress 0.09030393 
    ## Run 343 stress 0.08773465 
    ## Run 344 stress 0.09380597 
    ## Run 345 stress 0.1038038 
    ## Run 346 stress 0.09400535 
    ## Run 347 stress 0.09337288 
    ## Run 348 stress 0.08773481 
    ## Run 349 stress 0.08973862 
    ## Run 350 stress 0.090304 
    ## Run 351 stress 0.2772424 
    ## Run 352 stress 0.09407982 
    ## Run 353 stress 0.09465905 
    ## Run 354 stress 0.0930898 
    ## Run 355 stress 0.09286085 
    ## Run 356 stress 0.09969973 
    ## Run 357 stress 0.08440252 
    ## ... Procrustes: rmse 5.664472e-05  max resid 0.0001556455 
    ## ... Similar to previous best
    ## Run 358 stress 0.09168936 
    ## Run 359 stress 0.09404379 
    ## Run 360 stress 0.08503489 
    ## Run 361 stress 0.08440275 
    ## ... Procrustes: rmse 0.0002729005  max resid 0.0005674682 
    ## ... Similar to previous best
    ## Run 362 stress 0.09030404 
    ## Run 363 stress 0.09445472 
    ## Run 364 stress 0.09308958 
    ## Run 365 stress 0.09030399 
    ## Run 366 stress 0.09445489 
    ## Run 367 stress 0.09381845 
    ## Run 368 stress 0.09374313 
    ## Run 369 stress 0.09374371 
    ## Run 370 stress 0.08773469 
    ## Run 371 stress 0.08503584 
    ## Run 372 stress 0.1038049 
    ## Run 373 stress 0.09380591 
    ## Run 374 stress 0.09030393 
    ## Run 375 stress 0.08503455 
    ## Run 376 stress 0.09669842 
    ## Run 377 stress 0.09168937 
    ## Run 378 stress 0.3191846 
    ## Run 379 stress 0.09168927 
    ## Run 380 stress 0.08503485 
    ## Run 381 stress 0.08773466 
    ## Run 382 stress 0.09447306 
    ## Run 383 stress 0.09465909 
    ## Run 384 stress 0.09145321 
    ## Run 385 stress 0.09416136 
    ## Run 386 stress 0.09286089 
    ## Run 387 stress 0.09145337 
    ## Run 388 stress 0.09407118 
    ## Run 389 stress 0.09416135 
    ## Run 390 stress 0.0915909 
    ## Run 391 stress 0.09286084 
    ## Run 392 stress 0.09308967 
    ## Run 393 stress 0.09374262 
    ## Run 394 stress 0.08503471 
    ## Run 395 stress 0.08503501 
    ## Run 396 stress 0.1095796 
    ## Run 397 stress 0.09145333 
    ## Run 398 stress 0.09286085 
    ## Run 399 stress 0.09445502 
    ## Run 400 stress 0.09286091 
    ## Run 401 stress 0.09159083 
    ## Run 402 stress 0.09407307 
    ## Run 403 stress 0.08773479 
    ## Run 404 stress 0.08440273 
    ## ... Procrustes: rmse 0.0002812859  max resid 0.0005125415 
    ## ... Similar to previous best
    ## Run 405 stress 0.09374328 
    ## Run 406 stress 0.0850363 
    ## Run 407 stress 0.09465921 
    ## Run 408 stress 0.09286084 
    ## Run 409 stress 0.09159088 
    ## Run 410 stress 0.3347612 
    ## Run 411 stress 0.09030395 
    ## Run 412 stress 0.09268353 
    ## Run 413 stress 0.09969962 
    ## Run 414 stress 0.08503492 
    ## Run 415 stress 0.0940797 
    ## Run 416 stress 0.09380575 
    ## Run 417 stress 0.0897388 
    ## Run 418 stress 0.09337258 
    ## Run 419 stress 0.2692051 
    ## Run 420 stress 0.09407964 
    ## Run 421 stress 0.08503487 
    ## Run 422 stress 0.09168961 
    ## Run 423 stress 0.09308968 
    ## Run 424 stress 0.09465914 
    ## Run 425 stress 0.08440263 
    ## ... Procrustes: rmse 0.0002066062  max resid 0.0003773076 
    ## ... Similar to previous best
    ## Run 426 stress 0.09337254 
    ## Run 427 stress 0.09380576 
    ## Run 428 stress 0.09168945 
    ## Run 429 stress 0.08503619 
    ## Run 430 stress 0.08973879 
    ## Run 431 stress 0.08773476 
    ## Run 432 stress 0.08773474 
    ## Run 433 stress 0.09403415 
    ## Run 434 stress 0.09381847 
    ## Run 435 stress 0.09407976 
    ## Run 436 stress 0.09407122 
    ## Run 437 stress 0.09337291 
    ## Run 438 stress 0.09145324 
    ## Run 439 stress 0.09159098 
    ## Run 440 stress 0.09465905 
    ## Run 441 stress 0.1038042 
    ## Run 442 stress 0.09286086 
    ## Run 443 stress 0.08773468 
    ## Run 444 stress 0.09416144 
    ## Run 445 stress 0.09159097 
    ## Run 446 stress 0.09969975 
    ## Run 447 stress 0.09417764 
    ## Run 448 stress 0.09337291 
    ## Run 449 stress 0.09416131 
    ## Run 450 stress 0.09447311 
    ## Run 451 stress 0.09145322 
    ## Run 452 stress 0.09374267 
    ## Run 453 stress 0.08440274 
    ## ... Procrustes: rmse 0.0002758579  max resid 0.0005259953 
    ## ... Similar to previous best
    ## Run 454 stress 0.0944731 
    ## Run 455 stress 0.09145335 
    ## Run 456 stress 0.08503464 
    ## Run 457 stress 0.08440275 
    ## ... Procrustes: rmse 0.0002612117  max resid 0.0004758861 
    ## ... Similar to previous best
    ## Run 458 stress 0.08503582 
    ## Run 459 stress 0.09590827 
    ## Run 460 stress 0.09321809 
    ## Run 461 stress 0.09407975 
    ## Run 462 stress 0.08440261 
    ## ... Procrustes: rmse 0.0001937748  max resid 0.0003546034 
    ## ... Similar to previous best
    ## Run 463 stress 0.09159084 
    ## Run 464 stress 0.0959063 
    ## Run 465 stress 0.09286086 
    ## Run 466 stress 0.09445464 
    ## Run 467 stress 0.09400544 
    ## Run 468 stress 0.08503611 
    ## Run 469 stress 0.08973873 
    ## Run 470 stress 0.0926832 
    ## Run 471 stress 0.09159085 
    ## Run 472 stress 0.1052849 
    ## Run 473 stress 0.09590949 
    ## Run 474 stress 0.09159088 
    ## Run 475 stress 0.0941613 
    ## Run 476 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001134099  max resid 0.0002329767 
    ## ... Similar to previous best
    ## Run 477 stress 0.09407139 
    ## Run 478 stress 0.250235 
    ## Run 479 stress 0.09168961 
    ## Run 480 stress 0.09407981 
    ## Run 481 stress 0.08503503 
    ## Run 482 stress 0.09465905 
    ## Run 483 stress 0.09030405 
    ## Run 484 stress 0.09030406 
    ## Run 485 stress 0.09407976 
    ## Run 486 stress 0.08973868 
    ## Run 487 stress 0.09416131 
    ## Run 488 stress 0.09465919 
    ## Run 489 stress 0.0844026 
    ## ... Procrustes: rmse 0.0001846265  max resid 0.0003845609 
    ## ... Similar to previous best
    ## Run 490 stress 0.09308956 
    ## Run 491 stress 0.09400517 
    ## Run 492 stress 0.09407977 
    ## Run 493 stress 0.1038047 
    ## Run 494 stress 0.08773471 
    ## Run 495 stress 0.08503523 
    ## Run 496 stress 0.09465905 
    ## Run 497 stress 0.09286083 
    ## Run 498 stress 0.09030398 
    ## Run 499 stress 0.090304 
    ## Run 500 stress 0.08973868 
    ## *** Best solution repeated 12 times

``` r
# Ocean sites and mixed lakes
SD_beta_env_OM_NMDS <- metaMDS(SD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06623458 
    ## Run 1 stress 0.09468937 
    ## Run 2 stress 0.07166965 
    ## Run 3 stress 0.06623458 
    ## ... Procrustes: rmse 9.102238e-05  max resid 0.0001439618 
    ## ... Similar to previous best
    ## Run 4 stress 0.06477946 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1534115  max resid 0.2655186 
    ## Run 5 stress 0.0647796 
    ## ... Procrustes: rmse 8.404286e-05  max resid 0.0001340701 
    ## ... Similar to previous best
    ## Run 6 stress 0.08226172 
    ## Run 7 stress 0.07411955 
    ## Run 8 stress 0.06623459 
    ## Run 9 stress 0.08226159 
    ## Run 10 stress 0.06477831 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001221808  max resid 0.002050114 
    ## ... Similar to previous best
    ## Run 11 stress 0.08226168 
    ## Run 12 stress 0.06623457 
    ## Run 13 stress 0.06477878 
    ## ... Procrustes: rmse 0.0007092622  max resid 0.00118992 
    ## ... Similar to previous best
    ## Run 14 stress 0.07411946 
    ## Run 15 stress 0.07411949 
    ## Run 16 stress 0.0822616 
    ## Run 17 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 0.08936929  max resid 0.2152668 
    ## Run 18 stress 0.08226157 
    ## Run 19 stress 0.0716697 
    ## Run 20 stress 0.2394588 
    ## Run 21 stress 0.06477963 
    ## Run 22 stress 0.06477956 
    ## Run 23 stress 0.06477857 
    ## Run 24 stress 0.06124525 
    ## ... Procrustes: rmse 0.01527423  max resid 0.04227979 
    ## Run 25 stress 0.06477954 
    ## Run 26 stress 0.06477853 
    ## Run 27 stress 0.0612452 
    ## ... Procrustes: rmse 0.01524887  max resid 0.0422138 
    ## Run 28 stress 0.2362366 
    ## Run 29 stress 0.06623459 
    ## Run 30 stress 0.0716697 
    ## Run 31 stress 0.07411939 
    ## Run 32 stress 0.06623463 
    ## Run 33 stress 0.06623458 
    ## Run 34 stress 0.07166965 
    ## Run 35 stress 0.06623462 
    ## Run 36 stress 0.07411941 
    ## Run 37 stress 0.07166981 
    ## Run 38 stress 0.06623458 
    ## Run 39 stress 0.07411943 
    ## Run 40 stress 0.06623467 
    ## Run 41 stress 0.06477879 
    ## Run 42 stress 0.06124503 
    ## ... Procrustes: rmse 0.0144499  max resid 0.04007272 
    ## Run 43 stress 0.06477851 
    ## Run 44 stress 0.07166976 
    ## Run 45 stress 0.06623462 
    ## Run 46 stress 0.06623467 
    ## Run 47 stress 0.06623463 
    ## Run 48 stress 0.0741194 
    ## Run 49 stress 0.07411944 
    ## Run 50 stress 0.0647794 
    ## Run 51 stress 0.0612336 
    ## ... Procrustes: rmse 6.775862e-06  max resid 1.480571e-05 
    ## ... Similar to previous best
    ## Run 52 stress 0.07411942 
    ## Run 53 stress 0.06477898 
    ## Run 54 stress 0.06623457 
    ## Run 55 stress 0.08226163 
    ## Run 56 stress 0.07166964 
    ## Run 57 stress 0.06124511 
    ## ... Procrustes: rmse 0.01514561  max resid 0.04193978 
    ## Run 58 stress 0.3132604 
    ## Run 59 stress 0.07411939 
    ## Run 60 stress 0.06623458 
    ## Run 61 stress 0.06477862 
    ## Run 62 stress 0.07166974 
    ## Run 63 stress 0.0647796 
    ## Run 64 stress 0.06623467 
    ## Run 65 stress 0.07411956 
    ## Run 66 stress 0.061245 
    ## ... Procrustes: rmse 0.01507917  max resid 0.04175968 
    ## Run 67 stress 0.06477956 
    ## Run 68 stress 0.08226171 
    ## Run 69 stress 0.07166971 
    ## Run 70 stress 0.07411951 
    ## Run 71 stress 0.06623468 
    ## Run 72 stress 0.06123361 
    ## ... Procrustes: rmse 2.412551e-05  max resid 3.718916e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.06623459 
    ## Run 74 stress 0.0612449 
    ## ... Procrustes: rmse 0.01495062  max resid 0.04141484 
    ## Run 75 stress 0.06124512 
    ## ... Procrustes: rmse 0.01436801  max resid 0.03985236 
    ## Run 76 stress 0.07166982 
    ## Run 77 stress 0.07166965 
    ## Run 78 stress 0.07166964 
    ## Run 79 stress 0.07411943 
    ## Run 80 stress 0.0741194 
    ## Run 81 stress 0.0612336 
    ## ... Procrustes: rmse 6.144387e-06  max resid 1.403185e-05 
    ## ... Similar to previous best
    ## Run 82 stress 0.06123361 
    ## ... Procrustes: rmse 3.063537e-05  max resid 4.631179e-05 
    ## ... Similar to previous best
    ## Run 83 stress 0.06124494 
    ## ... Procrustes: rmse 0.01500239  max resid 0.04155483 
    ## Run 84 stress 0.09468857 
    ## Run 85 stress 0.07411951 
    ## Run 86 stress 0.06623467 
    ## Run 87 stress 0.06623458 
    ## Run 88 stress 0.08226161 
    ## Run 89 stress 0.06124491 
    ## ... Procrustes: rmse 0.0145899  max resid 0.04044794 
    ## Run 90 stress 0.06124525 
    ## ... Procrustes: rmse 0.0152802  max resid 0.04229672 
    ## Run 91 stress 0.06124522 
    ## ... Procrustes: rmse 0.01526009  max resid 0.04224248 
    ## Run 92 stress 0.06477932 
    ## Run 93 stress 0.06623461 
    ## Run 94 stress 0.07411947 
    ## Run 95 stress 0.09469093 
    ## Run 96 stress 0.06123361 
    ## ... Procrustes: rmse 4.497211e-05  max resid 6.910979e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.07411944 
    ## Run 98 stress 0.07166973 
    ## Run 99 stress 0.06623457 
    ## Run 100 stress 0.06477899 
    ## Run 101 stress 0.3131441 
    ## Run 102 stress 0.06477999 
    ## Run 103 stress 0.06124486 
    ## ... Procrustes: rmse 0.01484632  max resid 0.04113489 
    ## Run 104 stress 0.06124507 
    ## ... Procrustes: rmse 0.01514731  max resid 0.04194145 
    ## Run 105 stress 0.06124515 
    ## ... Procrustes: rmse 0.01520838  max resid 0.04210558 
    ## Run 106 stress 0.06477867 
    ## Run 107 stress 0.07411943 
    ## Run 108 stress 0.07411945 
    ## Run 109 stress 0.08226172 
    ## Run 110 stress 0.07166981 
    ## Run 111 stress 0.07411939 
    ## Run 112 stress 0.07166965 
    ## Run 113 stress 0.08226175 
    ## Run 114 stress 0.2362366 
    ## Run 115 stress 0.06623461 
    ## Run 116 stress 0.06477862 
    ## Run 117 stress 0.08226156 
    ## Run 118 stress 0.06477898 
    ## Run 119 stress 0.2324715 
    ## Run 120 stress 0.06477846 
    ## Run 121 stress 0.06478002 
    ## Run 122 stress 0.0612452 
    ## ... Procrustes: rmse 0.01430044  max resid 0.03967034 
    ## Run 123 stress 0.07411939 
    ## Run 124 stress 0.07166967 
    ## Run 125 stress 0.06623461 
    ## Run 126 stress 0.06623458 
    ## Run 127 stress 0.06124504 
    ## ... Procrustes: rmse 0.01510544  max resid 0.04183127 
    ## Run 128 stress 0.06623458 
    ## Run 129 stress 0.2394588 
    ## Run 130 stress 0.07166975 
    ## Run 131 stress 0.06477883 
    ## Run 132 stress 0.07166969 
    ## Run 133 stress 0.07411943 
    ## Run 134 stress 0.07411942 
    ## Run 135 stress 0.08226172 
    ## Run 136 stress 0.06477837 
    ## Run 137 stress 0.06124497 
    ## ... Procrustes: rmse 0.01503893  max resid 0.04165178 
    ## Run 138 stress 0.06477848 
    ## Run 139 stress 0.06477945 
    ## Run 140 stress 0.06623457 
    ## Run 141 stress 0.06623459 
    ## Run 142 stress 0.08226179 
    ## Run 143 stress 0.06477867 
    ## Run 144 stress 0.06477868 
    ## Run 145 stress 0.0612451 
    ## ... Procrustes: rmse 0.01517613  max resid 0.04201842 
    ## Run 146 stress 0.07411954 
    ## Run 147 stress 0.07166964 
    ## Run 148 stress 0.06124493 
    ## ... Procrustes: rmse 0.01499568  max resid 0.04153658 
    ## Run 149 stress 0.07411938 
    ## Run 150 stress 0.08226166 
    ## Run 151 stress 0.06477984 
    ## Run 152 stress 0.06623467 
    ## Run 153 stress 0.06623463 
    ## Run 154 stress 0.06477932 
    ## Run 155 stress 0.07411941 
    ## Run 156 stress 0.07411941 
    ## Run 157 stress 0.07411942 
    ## Run 158 stress 0.07411939 
    ## Run 159 stress 0.3022789 
    ## Run 160 stress 0.2548127 
    ## Run 161 stress 0.06124499 
    ## ... Procrustes: rmse 0.01501995  max resid 0.04160321 
    ## Run 162 stress 0.06623463 
    ## Run 163 stress 0.06623457 
    ## Run 164 stress 0.3266727 
    ## Run 165 stress 0.06623458 
    ## Run 166 stress 0.06124502 
    ## ... Procrustes: rmse 0.0151041  max resid 0.04182641 
    ## Run 167 stress 0.07166981 
    ## Run 168 stress 0.06477846 
    ## Run 169 stress 0.0716697 
    ## Run 170 stress 0.08226155 
    ## Run 171 stress 0.06623458 
    ## Run 172 stress 0.0662347 
    ## Run 173 stress 0.08226169 
    ## Run 174 stress 0.0612336 
    ## ... Procrustes: rmse 1.503299e-05  max resid 2.76867e-05 
    ## ... Similar to previous best
    ## Run 175 stress 0.07411941 
    ## Run 176 stress 0.07166984 
    ## Run 177 stress 0.06623464 
    ## Run 178 stress 0.06623461 
    ## Run 179 stress 0.2381016 
    ## Run 180 stress 0.06623459 
    ## Run 181 stress 0.07411948 
    ## Run 182 stress 0.2324715 
    ## Run 183 stress 0.07166974 
    ## Run 184 stress 0.0822617 
    ## Run 185 stress 0.06124525 
    ## ... Procrustes: rmse 0.01527999  max resid 0.04229551 
    ## Run 186 stress 0.07166975 
    ## Run 187 stress 0.08226173 
    ## Run 188 stress 0.06623458 
    ## Run 189 stress 0.0662346 
    ## Run 190 stress 0.06623458 
    ## Run 191 stress 0.0716697 
    ## Run 192 stress 0.06623458 
    ## Run 193 stress 0.06124511 
    ## ... Procrustes: rmse 0.01517768  max resid 0.04202197 
    ## Run 194 stress 0.07411947 
    ## Run 195 stress 0.2345379 
    ## Run 196 stress 0.07411952 
    ## Run 197 stress 0.2362366 
    ## Run 198 stress 0.07411938 
    ## Run 199 stress 0.08226158 
    ## Run 200 stress 0.06477839 
    ## Run 201 stress 0.06623464 
    ## Run 202 stress 0.06124499 
    ## ... Procrustes: rmse 0.01507658  max resid 0.04175206 
    ## Run 203 stress 0.06623458 
    ## Run 204 stress 0.08226157 
    ## Run 205 stress 0.06123361 
    ## ... Procrustes: rmse 3.757509e-05  max resid 8.564611e-05 
    ## ... Similar to previous best
    ## Run 206 stress 0.2629205 
    ## Run 207 stress 0.06623461 
    ## Run 208 stress 0.0647784 
    ## Run 209 stress 0.07166982 
    ## Run 210 stress 0.07411948 
    ## Run 211 stress 0.07411938 
    ## Run 212 stress 0.06623458 
    ## Run 213 stress 0.07411938 
    ## Run 214 stress 0.07411939 
    ## Run 215 stress 0.07411948 
    ## Run 216 stress 0.31192 
    ## Run 217 stress 0.06477877 
    ## Run 218 stress 0.08226169 
    ## Run 219 stress 0.07411938 
    ## Run 220 stress 0.07166992 
    ## Run 221 stress 0.07411943 
    ## Run 222 stress 0.06623461 
    ## Run 223 stress 0.06477915 
    ## Run 224 stress 0.08226155 
    ## Run 225 stress 0.3275212 
    ## Run 226 stress 0.06623457 
    ## Run 227 stress 0.06623457 
    ## Run 228 stress 0.06123361 
    ## ... Procrustes: rmse 4.459174e-05  max resid 5.393787e-05 
    ## ... Similar to previous best
    ## Run 229 stress 0.07411948 
    ## Run 230 stress 0.07166984 
    ## Run 231 stress 0.06124507 
    ## ... Procrustes: rmse 0.01440897  max resid 0.03996171 
    ## Run 232 stress 0.0612336 
    ## ... Procrustes: rmse 1.3514e-05  max resid 2.363323e-05 
    ## ... Similar to previous best
    ## Run 233 stress 0.07166964 
    ## Run 234 stress 0.07411952 
    ## Run 235 stress 0.0741194 
    ## Run 236 stress 0.06124495 
    ## ... Procrustes: rmse 0.01501918  max resid 0.0415987 
    ## Run 237 stress 0.0716697 
    ## Run 238 stress 0.0612336 
    ## ... Procrustes: rmse 1.516564e-05  max resid 3.5498e-05 
    ## ... Similar to previous best
    ## Run 239 stress 0.06623461 
    ## Run 240 stress 0.0741194 
    ## Run 241 stress 0.2381016 
    ## Run 242 stress 0.07411939 
    ## Run 243 stress 0.06124523 
    ## ... Procrustes: rmse 0.01438629  max resid 0.03990503 
    ## Run 244 stress 0.06124517 
    ## ... Procrustes: rmse 0.01522209  max resid 0.04214083 
    ## Run 245 stress 0.3092172 
    ## Run 246 stress 0.09623037 
    ## Run 247 stress 0.07166965 
    ## Run 248 stress 0.06623457 
    ## Run 249 stress 0.07166965 
    ## Run 250 stress 0.07166985 
    ## Run 251 stress 0.06623459 
    ## Run 252 stress 0.06477837 
    ## Run 253 stress 0.07166965 
    ## Run 254 stress 0.08226165 
    ## Run 255 stress 0.06623458 
    ## Run 256 stress 0.06123361 
    ## ... Procrustes: rmse 3.786204e-05  max resid 6.215773e-05 
    ## ... Similar to previous best
    ## Run 257 stress 0.07411949 
    ## Run 258 stress 0.2569846 
    ## Run 259 stress 0.0822617 
    ## Run 260 stress 0.07166975 
    ## Run 261 stress 0.0647792 
    ## Run 262 stress 0.08226164 
    ## Run 263 stress 0.07166986 
    ## Run 264 stress 0.0716697 
    ## Run 265 stress 0.07166965 
    ## Run 266 stress 0.07166989 
    ## Run 267 stress 0.2569851 
    ## Run 268 stress 0.0662346 
    ## Run 269 stress 0.06623462 
    ## Run 270 stress 0.08226177 
    ## Run 271 stress 0.07166972 
    ## Run 272 stress 0.06623462 
    ## Run 273 stress 0.06623458 
    ## Run 274 stress 0.06124494 
    ## ... Procrustes: rmse 0.01500764  max resid 0.04156918 
    ## Run 275 stress 0.0930082 
    ## Run 276 stress 0.0741194 
    ## Run 277 stress 0.07166967 
    ## Run 278 stress 0.06623459 
    ## Run 279 stress 0.06477922 
    ## Run 280 stress 0.06623458 
    ## Run 281 stress 0.07411939 
    ## Run 282 stress 0.0741194 
    ## Run 283 stress 0.0822616 
    ## Run 284 stress 0.06124508 
    ## ... Procrustes: rmse 0.01515545  max resid 0.04196218 
    ## Run 285 stress 0.08226156 
    ## Run 286 stress 0.0647797 
    ## Run 287 stress 0.07411938 
    ## Run 288 stress 0.07166977 
    ## Run 289 stress 0.06623461 
    ## Run 290 stress 0.07411938 
    ## Run 291 stress 0.07166985 
    ## Run 292 stress 0.07166983 
    ## Run 293 stress 0.06623462 
    ## Run 294 stress 0.06124495 
    ## ... Procrustes: rmse 0.01502405  max resid 0.04161152 
    ## Run 295 stress 0.06123362 
    ## ... Procrustes: rmse 5.888168e-05  max resid 7.963548e-05 
    ## ... Similar to previous best
    ## Run 296 stress 0.06477834 
    ## Run 297 stress 0.06477877 
    ## Run 298 stress 0.06124512 
    ## ... Procrustes: rmse 0.01519067  max resid 0.04205797 
    ## Run 299 stress 0.06623457 
    ## Run 300 stress 0.06623457 
    ## Run 301 stress 0.06623458 
    ## Run 302 stress 0.06477875 
    ## Run 303 stress 0.06623458 
    ## Run 304 stress 0.06477966 
    ## Run 305 stress 0.07411939 
    ## Run 306 stress 0.06623465 
    ## Run 307 stress 0.07166973 
    ## Run 308 stress 0.07411951 
    ## Run 309 stress 0.2548998 
    ## Run 310 stress 0.07411947 
    ## Run 311 stress 0.0647787 
    ## Run 312 stress 0.0612336 
    ## ... Procrustes: rmse 2.066436e-05  max resid 2.596851e-05 
    ## ... Similar to previous best
    ## Run 313 stress 0.08226165 
    ## Run 314 stress 0.07411939 
    ## Run 315 stress 0.0647788 
    ## Run 316 stress 0.06623461 
    ## Run 317 stress 0.07166972 
    ## Run 318 stress 0.06477905 
    ## Run 319 stress 0.07411938 
    ## Run 320 stress 0.0612336 
    ## ... Procrustes: rmse 5.104208e-06  max resid 9.334498e-06 
    ## ... Similar to previous best
    ## Run 321 stress 0.07411942 
    ## Run 322 stress 0.08226159 
    ## Run 323 stress 0.06623457 
    ## Run 324 stress 0.061245 
    ## ... Procrustes: rmse 0.01447133  max resid 0.04012905 
    ## Run 325 stress 0.07166973 
    ## Run 326 stress 0.06623468 
    ## Run 327 stress 0.06623463 
    ## Run 328 stress 0.07166977 
    ## Run 329 stress 0.08226164 
    ## Run 330 stress 0.07411942 
    ## Run 331 stress 0.07166966 
    ## Run 332 stress 0.0612452 
    ## ... Procrustes: rmse 0.01525037  max resid 0.04221713 
    ## Run 333 stress 0.06623457 
    ## Run 334 stress 0.064779 
    ## Run 335 stress 0.06477975 
    ## Run 336 stress 0.06477874 
    ## Run 337 stress 0.0822616 
    ## Run 338 stress 0.07166969 
    ## Run 339 stress 0.07411948 
    ## Run 340 stress 0.07411955 
    ## Run 341 stress 0.07411952 
    ## Run 342 stress 0.07411949 
    ## Run 343 stress 0.07411949 
    ## Run 344 stress 0.06623457 
    ## Run 345 stress 0.07411953 
    ## Run 346 stress 0.07411944 
    ## Run 347 stress 0.3407338 
    ## Run 348 stress 0.06477956 
    ## Run 349 stress 0.06123361 
    ## ... Procrustes: rmse 4.460391e-05  max resid 6.487879e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.07166975 
    ## Run 351 stress 0.07411941 
    ## Run 352 stress 0.06124504 
    ## ... Procrustes: rmse 0.01511964  max resid 0.0418683 
    ## Run 353 stress 0.0612336 
    ## ... Procrustes: rmse 4.988285e-06  max resid 9.869566e-06 
    ## ... Similar to previous best
    ## Run 354 stress 0.07411947 
    ## Run 355 stress 0.07411941 
    ## Run 356 stress 0.06623458 
    ## Run 357 stress 0.06477932 
    ## Run 358 stress 0.2381016 
    ## Run 359 stress 0.07411944 
    ## Run 360 stress 0.06623459 
    ## Run 361 stress 0.06623458 
    ## Run 362 stress 0.06623462 
    ## Run 363 stress 0.07411957 
    ## Run 364 stress 0.08226171 
    ## Run 365 stress 0.08226158 
    ## Run 366 stress 0.06623458 
    ## Run 367 stress 0.07411956 
    ## Run 368 stress 0.06123362 
    ## ... Procrustes: rmse 6.743864e-05  max resid 9.689169e-05 
    ## ... Similar to previous best
    ## Run 369 stress 0.06123363 
    ## ... Procrustes: rmse 7.064702e-05  max resid 0.0001018396 
    ## ... Similar to previous best
    ## Run 370 stress 0.07166969 
    ## Run 371 stress 0.2394588 
    ## Run 372 stress 0.06623459 
    ## Run 373 stress 0.07411939 
    ## Run 374 stress 0.07166976 
    ## Run 375 stress 0.06477869 
    ## Run 376 stress 0.06623459 
    ## Run 377 stress 0.06124516 
    ## ... Procrustes: rmse 0.01520778  max resid 0.04210166 
    ## Run 378 stress 0.06477887 
    ## Run 379 stress 0.07411947 
    ## Run 380 stress 0.06477862 
    ## Run 381 stress 0.0647798 
    ## Run 382 stress 0.06623463 
    ## Run 383 stress 0.06623463 
    ## Run 384 stress 0.06124501 
    ## ... Procrustes: rmse 0.01509928  max resid 0.04181335 
    ## Run 385 stress 0.06623461 
    ## Run 386 stress 0.07166964 
    ## Run 387 stress 0.06477907 
    ## Run 388 stress 0.07411941 
    ## Run 389 stress 0.3257213 
    ## Run 390 stress 0.06623466 
    ## Run 391 stress 0.0612449 
    ## ... Procrustes: rmse 0.01495019  max resid 0.04141386 
    ## Run 392 stress 0.07166964 
    ## Run 393 stress 0.06623462 
    ## Run 394 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 4.252324e-06  max resid 9.218351e-06 
    ## ... Similar to previous best
    ## Run 395 stress 0.07411941 
    ## Run 396 stress 0.07166984 
    ## Run 397 stress 0.06623459 
    ## Run 398 stress 0.06123364 
    ## ... Procrustes: rmse 8.382955e-05  max resid 0.000109156 
    ## ... Similar to previous best
    ## Run 399 stress 0.07411946 
    ## Run 400 stress 0.06124519 
    ## ... Procrustes: rmse 0.01520108  max resid 0.04208253 
    ## Run 401 stress 0.06124511 
    ## ... Procrustes: rmse 0.01518325  max resid 0.04203845 
    ## Run 402 stress 0.06477857 
    ## Run 403 stress 0.0741194 
    ## Run 404 stress 0.06623458 
    ## Run 405 stress 0.07411939 
    ## Run 406 stress 0.07411947 
    ## Run 407 stress 0.0647786 
    ## Run 408 stress 0.06477875 
    ## Run 409 stress 0.07411945 
    ## Run 410 stress 0.06124494 
    ## ... Procrustes: rmse 0.0150044  max resid 0.04156008 
    ## Run 411 stress 0.06477895 
    ## Run 412 stress 0.07166964 
    ## Run 413 stress 0.07411942 
    ## Run 414 stress 0.06623459 
    ## Run 415 stress 0.07411951 
    ## Run 416 stress 0.0647787 
    ## Run 417 stress 0.06623459 
    ## Run 418 stress 0.07411948 
    ## Run 419 stress 0.07166969 
    ## Run 420 stress 0.0741194 
    ## Run 421 stress 0.07166976 
    ## Run 422 stress 0.0612336 
    ## ... Procrustes: rmse 7.448502e-06  max resid 1.620068e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 1.48377e-06  max resid 3.266868e-06 
    ## ... Similar to previous best
    ## Run 424 stress 0.06623459 
    ## Run 425 stress 0.07411945 
    ## Run 426 stress 0.07166979 
    ## Run 427 stress 0.2361811 
    ## Run 428 stress 0.07411951 
    ## Run 429 stress 0.06477979 
    ## Run 430 stress 0.0612336 
    ## ... Procrustes: rmse 3.427132e-06  max resid 7.243951e-06 
    ## ... Similar to previous best
    ## Run 431 stress 0.06124497 
    ## ... Procrustes: rmse 0.0150414  max resid 0.04165901 
    ## Run 432 stress 0.06623462 
    ## Run 433 stress 0.2788939 
    ## Run 434 stress 0.06124504 
    ## ... Procrustes: rmse 0.01512475  max resid 0.04188153 
    ## Run 435 stress 0.06623461 
    ## Run 436 stress 0.06623464 
    ## Run 437 stress 0.07411939 
    ## Run 438 stress 0.07166965 
    ## Run 439 stress 0.3407336 
    ## Run 440 stress 0.06623458 
    ## Run 441 stress 0.06124522 
    ## ... Procrustes: rmse 0.01524486  max resid 0.04220093 
    ## Run 442 stress 0.06477925 
    ## Run 443 stress 0.06477841 
    ## Run 444 stress 0.06478007 
    ## Run 445 stress 0.318841 
    ## Run 446 stress 0.0741194 
    ## Run 447 stress 0.06623458 
    ## Run 448 stress 0.09300819 
    ## Run 449 stress 0.07411949 
    ## Run 450 stress 0.07411942 
    ## Run 451 stress 0.07411944 
    ## Run 452 stress 0.06124518 
    ## ... Procrustes: rmse 0.01523374  max resid 0.04217304 
    ## Run 453 stress 0.3131447 
    ## Run 454 stress 0.06623459 
    ## Run 455 stress 0.07166964 
    ## Run 456 stress 0.06477902 
    ## Run 457 stress 0.06623457 
    ## Run 458 stress 0.08226157 
    ## Run 459 stress 0.3407334 
    ## Run 460 stress 0.0741195 
    ## Run 461 stress 0.07166964 
    ## Run 462 stress 0.06477927 
    ## Run 463 stress 0.07411939 
    ## Run 464 stress 0.07166986 
    ## Run 465 stress 0.07411939 
    ## Run 466 stress 0.07166966 
    ## Run 467 stress 0.07411939 
    ## Run 468 stress 0.06477874 
    ## Run 469 stress 0.06623458 
    ## Run 470 stress 0.06124505 
    ## ... Procrustes: rmse 0.0144238  max resid 0.04000262 
    ## Run 471 stress 0.06623462 
    ## Run 472 stress 0.06477903 
    ## Run 473 stress 0.07411938 
    ## Run 474 stress 0.06477858 
    ## Run 475 stress 0.06623457 
    ## Run 476 stress 0.07166967 
    ## Run 477 stress 0.06477898 
    ## Run 478 stress 0.06124489 
    ## ... Procrustes: rmse 0.01462906  max resid 0.04055352 
    ## Run 479 stress 0.3131401 
    ## Run 480 stress 0.06623464 
    ## Run 481 stress 0.07166969 
    ## Run 482 stress 0.06124521 
    ## ... Procrustes: rmse 0.0152521  max resid 0.04222256 
    ## Run 483 stress 0.0741194 
    ## Run 484 stress 0.0612336 
    ## ... Procrustes: rmse 5.267521e-06  max resid 7.560318e-06 
    ## ... Similar to previous best
    ## Run 485 stress 0.0612336 
    ## ... Procrustes: rmse 6.389185e-06  max resid 8.878155e-06 
    ## ... Similar to previous best
    ## Run 486 stress 0.07411938 
    ## Run 487 stress 0.06623458 
    ## Run 488 stress 0.0612336 
    ## ... Procrustes: rmse 1.497401e-05  max resid 3.186498e-05 
    ## ... Similar to previous best
    ## Run 489 stress 0.08226171 
    ## Run 490 stress 0.06477847 
    ## Run 491 stress 0.07166978 
    ## Run 492 stress 0.0741194 
    ## Run 493 stress 0.06623462 
    ## Run 494 stress 0.07411941 
    ## Run 495 stress 0.06623464 
    ## Run 496 stress 0.3163823 
    ## Run 497 stress 0.06477856 
    ## Run 498 stress 0.07411939 
    ## Run 499 stress 0.06124511 
    ## ... Procrustes: rmse 0.01516177  max resid 0.04197867 
    ## Run 500 stress 0.3266728 
    ## *** Best solution repeated 5 times

``` r
# Stratified lakes and ocean sites
SD_beta_env_SO_NMDS <- metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.59483e-05 
    ## Run 1 stress 9.223633e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002312947  max resid 0.0004752459 
    ## ... Similar to previous best
    ## Run 2 stress 9.193509e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001680792  max resid 0.0003430405 
    ## ... Similar to previous best
    ## Run 3 stress 9.377116e-05 
    ## ... Procrustes: rmse 0.0001756523  max resid 0.0003243409 
    ## ... Similar to previous best
    ## Run 4 stress 6.638567e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000157136  max resid 0.000254787 
    ## ... Similar to previous best
    ## Run 5 stress 8.978155e-05 
    ## ... Procrustes: rmse 0.0001518929  max resid 0.0002564586 
    ## ... Similar to previous best
    ## Run 6 stress 9.99658e-05 
    ## ... Procrustes: rmse 0.0001569401  max resid 0.0002783499 
    ## ... Similar to previous best
    ## Run 7 stress 9.412639e-05 
    ## ... Procrustes: rmse 0.0001644417  max resid 0.0002643521 
    ## ... Similar to previous best
    ## Run 8 stress 9.881569e-05 
    ## ... Procrustes: rmse 0.0001838514  max resid 0.0002430791 
    ## ... Similar to previous best
    ## Run 9 stress 9.954425e-05 
    ## ... Procrustes: rmse 0.0001899413  max resid 0.0003399305 
    ## ... Similar to previous best
    ## Run 10 stress 9.326105e-05 
    ## ... Procrustes: rmse 0.0001676692  max resid 0.000226465 
    ## ... Similar to previous best
    ## Run 11 stress 9.945592e-05 
    ## ... Procrustes: rmse 0.0001388786  max resid 0.000248527 
    ## ... Similar to previous best
    ## Run 12 stress 9.858674e-05 
    ## ... Procrustes: rmse 0.000170175  max resid 0.0002467002 
    ## ... Similar to previous best
    ## Run 13 stress 9.643292e-05 
    ## ... Procrustes: rmse 0.0001992431  max resid 0.0003724107 
    ## ... Similar to previous best
    ## Run 14 stress 9.537787e-05 
    ## ... Procrustes: rmse 0.0001571221  max resid 0.0002462271 
    ## ... Similar to previous best
    ## Run 15 stress 8.998064e-05 
    ## ... Procrustes: rmse 0.0001666252  max resid 0.0003055136 
    ## ... Similar to previous best
    ## Run 16 stress 9.164079e-05 
    ## ... Procrustes: rmse 0.000159732  max resid 0.0003088936 
    ## ... Similar to previous best
    ## Run 17 stress 9.92111e-05 
    ## ... Procrustes: rmse 0.0001890847  max resid 0.0003301817 
    ## ... Similar to previous best
    ## Run 18 stress 9.38798e-05 
    ## ... Procrustes: rmse 0.000185303  max resid 0.000260112 
    ## ... Similar to previous best
    ## Run 19 stress 9.930602e-05 
    ## ... Procrustes: rmse 0.0002027389  max resid 0.0002904644 
    ## ... Similar to previous best
    ## Run 20 stress 9.760799e-05 
    ## ... Procrustes: rmse 0.0001636672  max resid 0.0002729452 
    ## ... Similar to previous best
    ## Run 21 stress 9.460935e-05 
    ## ... Procrustes: rmse 0.0001275426  max resid 0.0002293844 
    ## ... Similar to previous best
    ## Run 22 stress 8.155902e-05 
    ## ... Procrustes: rmse 0.0001190639  max resid 0.0002093566 
    ## ... Similar to previous best
    ## Run 23 stress 9.916047e-05 
    ## ... Procrustes: rmse 0.0001775531  max resid 0.000252558 
    ## ... Similar to previous best
    ## Run 24 stress 9.429274e-05 
    ## ... Procrustes: rmse 0.0001845724  max resid 0.0002618779 
    ## ... Similar to previous best
    ## Run 25 stress 9.514488e-05 
    ## ... Procrustes: rmse 0.0001744706  max resid 0.0002326312 
    ## ... Similar to previous best
    ## Run 26 stress 9.401297e-05 
    ## ... Procrustes: rmse 0.0001885195  max resid 0.0003587545 
    ## ... Similar to previous best
    ## Run 27 stress 9.307215e-05 
    ## ... Procrustes: rmse 0.0001407215  max resid 0.0002469108 
    ## ... Similar to previous best
    ## Run 28 stress 9.988988e-05 
    ## ... Procrustes: rmse 0.0001865364  max resid 0.0002454031 
    ## ... Similar to previous best
    ## Run 29 stress 9.514933e-05 
    ## ... Procrustes: rmse 0.000128567  max resid 0.0002281245 
    ## ... Similar to previous best
    ## Run 30 stress 9.596761e-05 
    ## ... Procrustes: rmse 0.0001940976  max resid 0.0002820187 
    ## ... Similar to previous best
    ## Run 31 stress 9.930043e-05 
    ## ... Procrustes: rmse 0.0001477245  max resid 0.0002043743 
    ## ... Similar to previous best
    ## Run 32 stress 9.942034e-05 
    ## ... Procrustes: rmse 0.0001666149  max resid 0.0003048509 
    ## ... Similar to previous best
    ## Run 33 stress 9.575804e-05 
    ## ... Procrustes: rmse 0.0001771845  max resid 0.0002305678 
    ## ... Similar to previous best
    ## Run 34 stress 9.9166e-05 
    ## ... Procrustes: rmse 0.0001917428  max resid 0.0002688857 
    ## ... Similar to previous best
    ## Run 35 stress 9.937133e-05 
    ## ... Procrustes: rmse 0.0001714072  max resid 0.0002675837 
    ## ... Similar to previous best
    ## Run 36 stress 9.164874e-05 
    ## ... Procrustes: rmse 0.0001259055  max resid 0.0001936962 
    ## ... Similar to previous best
    ## Run 37 stress 8.772659e-05 
    ## ... Procrustes: rmse 0.0001702564  max resid 0.0002464663 
    ## ... Similar to previous best
    ## Run 38 stress 9.96489e-05 
    ## ... Procrustes: rmse 0.0001595781  max resid 0.0002532308 
    ## ... Similar to previous best
    ## Run 39 stress 9.710714e-05 
    ## ... Procrustes: rmse 0.0001699903  max resid 0.0002337962 
    ## ... Similar to previous best
    ## Run 40 stress 9.358095e-05 
    ## ... Procrustes: rmse 0.0001924958  max resid 0.0003635954 
    ## ... Similar to previous best
    ## Run 41 stress 9.894538e-05 
    ## ... Procrustes: rmse 0.0002033217  max resid 0.0003833286 
    ## ... Similar to previous best
    ## Run 42 stress 9.663445e-05 
    ## ... Procrustes: rmse 0.0001636661  max resid 0.0002584808 
    ## ... Similar to previous best
    ## Run 43 stress 9.612057e-05 
    ## ... Procrustes: rmse 0.0001675715  max resid 0.000265141 
    ## ... Similar to previous best
    ## Run 44 stress 9.704634e-05 
    ## ... Procrustes: rmse 0.0001969731  max resid 0.000285072 
    ## ... Similar to previous best
    ## Run 45 stress 9.102965e-05 
    ## ... Procrustes: rmse 0.0001545269  max resid 0.000217516 
    ## ... Similar to previous best
    ## Run 46 stress 9.655077e-05 
    ## ... Procrustes: rmse 0.0001895829  max resid 0.0002674421 
    ## ... Similar to previous best
    ## Run 47 stress 9.758735e-05 
    ## ... Procrustes: rmse 0.0001279563  max resid 0.0002306545 
    ## ... Similar to previous best
    ## Run 48 stress 9.043905e-05 
    ## ... Procrustes: rmse 0.0001948701  max resid 0.0002782607 
    ## ... Similar to previous best
    ## Run 49 stress 9.31955e-05 
    ## ... Procrustes: rmse 0.0001750948  max resid 0.0003178706 
    ## ... Similar to previous best
    ## Run 50 stress 8.349232e-05 
    ## ... Procrustes: rmse 0.0001368933  max resid 0.0002236917 
    ## ... Similar to previous best
    ## Run 51 stress 9.346509e-05 
    ## ... Procrustes: rmse 0.0001748494  max resid 0.0003215763 
    ## ... Similar to previous best
    ## Run 52 stress 9.450891e-05 
    ## ... Procrustes: rmse 0.0001617865  max resid 0.0002589336 
    ## ... Similar to previous best
    ## Run 53 stress 9.595014e-05 
    ## ... Procrustes: rmse 0.0001848209  max resid 0.0002681901 
    ## ... Similar to previous best
    ## Run 54 stress 9.551536e-05 
    ## ... Procrustes: rmse 0.0001868727  max resid 0.0002636775 
    ## ... Similar to previous best
    ## Run 55 stress 9.49017e-05 
    ## ... Procrustes: rmse 0.0001940592  max resid 0.0003631648 
    ## ... Similar to previous best
    ## Run 56 stress 9.568571e-05 
    ## ... Procrustes: rmse 0.0001634151  max resid 0.0002532693 
    ## ... Similar to previous best
    ## Run 57 stress 8.844251e-05 
    ## ... Procrustes: rmse 0.0001483219  max resid 0.0002120603 
    ## ... Similar to previous best
    ## Run 58 stress 9.600089e-05 
    ## ... Procrustes: rmse 0.0001671763  max resid 0.000266939 
    ## ... Similar to previous best
    ## Run 59 stress 9.909326e-05 
    ## ... Procrustes: rmse 0.000172237  max resid 0.0002668631 
    ## ... Similar to previous best
    ## Run 60 stress 9.449876e-05 
    ## ... Procrustes: rmse 0.0001430378  max resid 0.0002267011 
    ## ... Similar to previous best
    ## Run 61 stress 8.857325e-05 
    ## ... Procrustes: rmse 0.000110559  max resid 0.0002013046 
    ## ... Similar to previous best
    ## Run 62 stress 9.047907e-05 
    ## ... Procrustes: rmse 0.0001198243  max resid 0.0002150028 
    ## ... Similar to previous best
    ## Run 63 stress 9.876061e-05 
    ## ... Procrustes: rmse 0.0001983267  max resid 0.0002820278 
    ## ... Similar to previous best
    ## Run 64 stress 9.528354e-05 
    ## ... Procrustes: rmse 0.0001787319  max resid 0.0003325396 
    ## ... Similar to previous best
    ## Run 65 stress 7.465618e-05 
    ## ... Procrustes: rmse 0.000117797  max resid 0.0002259549 
    ## ... Similar to previous best
    ## Run 66 stress 9.807591e-05 
    ## ... Procrustes: rmse 9.704719e-05  max resid 0.0001765711 
    ## ... Similar to previous best
    ## Run 67 stress 9.850381e-05 
    ## ... Procrustes: rmse 0.0001379634  max resid 0.0002436745 
    ## ... Similar to previous best
    ## Run 68 stress 9.389605e-05 
    ## ... Procrustes: rmse 0.0001727807  max resid 0.0002258043 
    ## ... Similar to previous best
    ## Run 69 stress 9.652158e-05 
    ## ... Procrustes: rmse 0.0001691563  max resid 0.0002280461 
    ## ... Similar to previous best
    ## Run 70 stress 8.822568e-05 
    ## ... Procrustes: rmse 0.000184474  max resid 0.0003213242 
    ## ... Similar to previous best
    ## Run 71 stress 9.648375e-05 
    ## ... Procrustes: rmse 0.0001545445  max resid 0.0002650903 
    ## ... Similar to previous best
    ## Run 72 stress 9.816928e-05 
    ## ... Procrustes: rmse 0.0001591467  max resid 0.0002265996 
    ## ... Similar to previous best
    ## Run 73 stress 9.804906e-05 
    ## ... Procrustes: rmse 0.0001380564  max resid 0.0002592024 
    ## ... Similar to previous best
    ## Run 74 stress 9.568106e-05 
    ## ... Procrustes: rmse 0.0001656755  max resid 0.0002296408 
    ## ... Similar to previous best
    ## Run 75 stress 9.948164e-05 
    ## ... Procrustes: rmse 0.0002032035  max resid 0.0003856257 
    ## ... Similar to previous best
    ## Run 76 stress 9.161997e-05 
    ## ... Procrustes: rmse 0.0001501027  max resid 0.000252919 
    ## ... Similar to previous best
    ## Run 77 stress 9.353838e-05 
    ## ... Procrustes: rmse 0.0001963255  max resid 0.0003728214 
    ## ... Similar to previous best
    ## Run 78 stress 8.842032e-05 
    ## ... Procrustes: rmse 0.0001765632  max resid 0.0003226583 
    ## ... Similar to previous best
    ## Run 79 stress 9.103106e-05 
    ## ... Procrustes: rmse 0.0001473346  max resid 0.000232315 
    ## ... Similar to previous best
    ## Run 80 stress 9.435541e-05 
    ## ... Procrustes: rmse 0.0001960785  max resid 0.0003390569 
    ## ... Similar to previous best
    ## Run 81 stress 0.273256 
    ## Run 82 stress 8.905977e-05 
    ## ... Procrustes: rmse 0.0001802671  max resid 0.0003403215 
    ## ... Similar to previous best
    ## Run 83 stress 9.139084e-05 
    ## ... Procrustes: rmse 0.0001775675  max resid 0.0002574818 
    ## ... Similar to previous best
    ## Run 84 stress 9.293432e-05 
    ## ... Procrustes: rmse 0.0001833706  max resid 0.0003188872 
    ## ... Similar to previous best
    ## Run 85 stress 9.968276e-05 
    ## ... Procrustes: rmse 0.0001561119  max resid 0.0002503769 
    ## ... Similar to previous best
    ## Run 86 stress 9.121638e-05 
    ## ... Procrustes: rmse 0.0001585951  max resid 0.0002581973 
    ## ... Similar to previous best
    ## Run 87 stress 0.2330486 
    ## Run 88 stress 9.503783e-05 
    ## ... Procrustes: rmse 0.0001573845  max resid 0.0002441898 
    ## ... Similar to previous best
    ## Run 89 stress 9.477843e-05 
    ## ... Procrustes: rmse 0.00017523  max resid 0.0003250007 
    ## ... Similar to previous best
    ## Run 90 stress 0.2362213 
    ## Run 91 stress 9.932971e-05 
    ## ... Procrustes: rmse 0.000176538  max resid 0.0003231115 
    ## ... Similar to previous best
    ## Run 92 stress 9.126886e-05 
    ## ... Procrustes: rmse 0.00015082  max resid 0.0003108281 
    ## ... Similar to previous best
    ## Run 93 stress 9.636063e-05 
    ## ... Procrustes: rmse 0.000113943  max resid 0.0001716555 
    ## ... Similar to previous best
    ## Run 94 stress 9.989693e-05 
    ## ... Procrustes: rmse 0.0001380296  max resid 0.0002477971 
    ## ... Similar to previous best
    ## Run 95 stress 9.764188e-05 
    ## ... Procrustes: rmse 0.0001876658  max resid 0.0003290264 
    ## ... Similar to previous best
    ## Run 96 stress 9.989268e-05 
    ## ... Procrustes: rmse 0.0001903753  max resid 0.0003501896 
    ## ... Similar to previous best
    ## Run 97 stress 9.513679e-05 
    ## ... Procrustes: rmse 0.0001825181  max resid 0.0003097739 
    ## ... Similar to previous best
    ## Run 98 stress 9.884485e-05 
    ## ... Procrustes: rmse 0.0002128482  max resid 0.0003867819 
    ## ... Similar to previous best
    ## Run 99 stress 9.166568e-05 
    ## ... Procrustes: rmse 0.000174419  max resid 0.0002533174 
    ## ... Similar to previous best
    ## Run 100 stress 9.857364e-05 
    ## ... Procrustes: rmse 0.0001685619  max resid 0.000266841 
    ## ... Similar to previous best
    ## Run 101 stress 8.987704e-05 
    ## ... Procrustes: rmse 0.0001540117  max resid 0.000288606 
    ## ... Similar to previous best
    ## Run 102 stress 9.791817e-05 
    ## ... Procrustes: rmse 0.0001448147  max resid 0.0002581497 
    ## ... Similar to previous best
    ## Run 103 stress 9.913934e-05 
    ## ... Procrustes: rmse 0.0001663648  max resid 0.0002546568 
    ## ... Similar to previous best
    ## Run 104 stress 8.670043e-05 
    ## ... Procrustes: rmse 0.0001644837  max resid 0.0002865242 
    ## ... Similar to previous best
    ## Run 105 stress 9.474979e-05 
    ## ... Procrustes: rmse 0.0001768557  max resid 0.0003024958 
    ## ... Similar to previous best
    ## Run 106 stress 9.720187e-05 
    ## ... Procrustes: rmse 0.0001701697  max resid 0.0002682529 
    ## ... Similar to previous best
    ## Run 107 stress 9.466135e-05 
    ## ... Procrustes: rmse 0.0001733907  max resid 0.0002276189 
    ## ... Similar to previous best
    ## Run 108 stress 9.190467e-05 
    ## ... Procrustes: rmse 0.0001744784  max resid 0.0003174543 
    ## ... Similar to previous best
    ## Run 109 stress 9.318479e-05 
    ## ... Procrustes: rmse 0.000118262  max resid 0.0001968305 
    ## ... Similar to previous best
    ## Run 110 stress 9.203759e-05 
    ## ... Procrustes: rmse 0.000140416  max resid 0.0002349566 
    ## ... Similar to previous best
    ## Run 111 stress 9.767133e-05 
    ## ... Procrustes: rmse 0.0001727865  max resid 0.0003311078 
    ## ... Similar to previous best
    ## Run 112 stress 9.51854e-05 
    ## ... Procrustes: rmse 0.0001781949  max resid 0.0003226936 
    ## ... Similar to previous best
    ## Run 113 stress 9.728676e-05 
    ## ... Procrustes: rmse 0.0001269036  max resid 0.0002280761 
    ## ... Similar to previous best
    ## Run 114 stress 9.670303e-05 
    ## ... Procrustes: rmse 0.0001919109  max resid 0.0002806139 
    ## ... Similar to previous best
    ## Run 115 stress 9.124539e-05 
    ## ... Procrustes: rmse 0.0001443963  max resid 0.0002083129 
    ## ... Similar to previous best
    ## Run 116 stress 9.564818e-05 
    ## ... Procrustes: rmse 0.0001920829  max resid 0.0003474543 
    ## ... Similar to previous best
    ## Run 117 stress 0.2391424 
    ## Run 118 stress 9.699357e-05 
    ## ... Procrustes: rmse 0.0001756396  max resid 0.0003236652 
    ## ... Similar to previous best
    ## Run 119 stress 9.44901e-05 
    ## ... Procrustes: rmse 0.0001935084  max resid 0.0002819235 
    ## ... Similar to previous best
    ## Run 120 stress 9.155516e-05 
    ## ... Procrustes: rmse 0.0001668553  max resid 0.0002210431 
    ## ... Similar to previous best
    ## Run 121 stress 9.871741e-05 
    ## ... Procrustes: rmse 0.0001254526  max resid 0.0002205872 
    ## ... Similar to previous best
    ## Run 122 stress 9.432557e-05 
    ## ... Procrustes: rmse 0.0001759114  max resid 0.0003104001 
    ## ... Similar to previous best
    ## Run 123 stress 9.516332e-05 
    ## ... Procrustes: rmse 0.0001650878  max resid 0.000265654 
    ## ... Similar to previous best
    ## Run 124 stress 9.976372e-05 
    ## ... Procrustes: rmse 0.0002052395  max resid 0.0002878764 
    ## ... Similar to previous best
    ## Run 125 stress 9.947001e-05 
    ## ... Procrustes: rmse 0.000139172  max resid 0.0002438813 
    ## ... Similar to previous best
    ## Run 126 stress 9.74269e-05 
    ## ... Procrustes: rmse 0.0001793917  max resid 0.0002809764 
    ## ... Similar to previous best
    ## Run 127 stress 8.847738e-05 
    ## ... Procrustes: rmse 0.0001257993  max resid 0.0002817355 
    ## ... Similar to previous best
    ## Run 128 stress 9.404538e-05 
    ## ... Procrustes: rmse 0.0001501077  max resid 0.0002127147 
    ## ... Similar to previous best
    ## Run 129 stress 9.615069e-05 
    ## ... Procrustes: rmse 0.0001326577  max resid 0.0002340454 
    ## ... Similar to previous best
    ## Run 130 stress 9.742831e-05 
    ## ... Procrustes: rmse 0.0001487813  max resid 0.0003168242 
    ## ... Similar to previous best
    ## Run 131 stress 9.985216e-05 
    ## ... Procrustes: rmse 0.0002038661  max resid 0.0002918999 
    ## ... Similar to previous best
    ## Run 132 stress 9.51555e-05 
    ## ... Procrustes: rmse 0.0001820159  max resid 0.0002613349 
    ## ... Similar to previous best
    ## Run 133 stress 9.277099e-05 
    ## ... Procrustes: rmse 0.0002018041  max resid 0.000375037 
    ## ... Similar to previous best
    ## Run 134 stress 9.713721e-05 
    ## ... Procrustes: rmse 0.0001279583  max resid 0.0002171007 
    ## ... Similar to previous best
    ## Run 135 stress 9.169185e-05 
    ## ... Procrustes: rmse 0.0001421278  max resid 0.0002326576 
    ## ... Similar to previous best
    ## Run 136 stress 9.863737e-05 
    ## ... Procrustes: rmse 0.0002026877  max resid 0.0003760588 
    ## ... Similar to previous best
    ## Run 137 stress 9.202499e-05 
    ## ... Procrustes: rmse 0.000156836  max resid 0.0002528325 
    ## ... Similar to previous best
    ## Run 138 stress 8.852613e-05 
    ## ... Procrustes: rmse 0.0001290459  max resid 0.0002838465 
    ## ... Similar to previous best
    ## Run 139 stress 9.466604e-05 
    ## ... Procrustes: rmse 0.000163074  max resid 0.0002611643 
    ## ... Similar to previous best
    ## Run 140 stress 8.745237e-05 
    ## ... Procrustes: rmse 0.0001449809  max resid 0.0002437385 
    ## ... Similar to previous best
    ## Run 141 stress 9.765717e-05 
    ## ... Procrustes: rmse 0.0001228794  max resid 0.0002521675 
    ## ... Similar to previous best
    ## Run 142 stress 9.707016e-05 
    ## ... Procrustes: rmse 0.000161234  max resid 0.0002583378 
    ## ... Similar to previous best
    ## Run 143 stress 9.046939e-05 
    ## ... Procrustes: rmse 0.0001064512  max resid 0.0002372997 
    ## ... Similar to previous best
    ## Run 144 stress 8.28243e-05 
    ## ... Procrustes: rmse 0.0001480114  max resid 0.0002384205 
    ## ... Similar to previous best
    ## Run 145 stress 9.579295e-05 
    ## ... Procrustes: rmse 0.0002003637  max resid 0.0003472856 
    ## ... Similar to previous best
    ## Run 146 stress 9.638195e-05 
    ## ... Procrustes: rmse 0.0002082736  max resid 0.0003821595 
    ## ... Similar to previous best
    ## Run 147 stress 9.768325e-05 
    ## ... Procrustes: rmse 0.0001588849  max resid 0.0003128833 
    ## ... Similar to previous best
    ## Run 148 stress 9.53109e-05 
    ## ... Procrustes: rmse 0.0001664452  max resid 0.0002684348 
    ## ... Similar to previous best
    ## Run 149 stress 9.436961e-05 
    ## ... Procrustes: rmse 0.0001676049  max resid 0.00025684 
    ## ... Similar to previous best
    ## Run 150 stress 9.693925e-05 
    ## ... Procrustes: rmse 0.0001652943  max resid 0.0002646773 
    ## ... Similar to previous best
    ## Run 151 stress 9.503265e-05 
    ## ... Procrustes: rmse 0.0001572957  max resid 0.0002437185 
    ## ... Similar to previous best
    ## Run 152 stress 9.715252e-05 
    ## ... Procrustes: rmse 0.0001661531  max resid 0.0002622498 
    ## ... Similar to previous best
    ## Run 153 stress 9.20982e-05 
    ## ... Procrustes: rmse 0.0001464932  max resid 0.0002432928 
    ## ... Similar to previous best
    ## Run 154 stress 8.053994e-05 
    ## ... Procrustes: rmse 0.0001153881  max resid 0.0001973453 
    ## ... Similar to previous best
    ## Run 155 stress 9.472493e-05 
    ## ... Procrustes: rmse 0.0001304254  max resid 0.0002326325 
    ## ... Similar to previous best
    ## Run 156 stress 9.053228e-05 
    ## ... Procrustes: rmse 0.0001535921  max resid 0.0003041192 
    ## ... Similar to previous best
    ## Run 157 stress 9.21688e-05 
    ## ... Procrustes: rmse 0.0001881277  max resid 0.0003550713 
    ## ... Similar to previous best
    ## Run 158 stress 9.361659e-05 
    ## ... Procrustes: rmse 0.0001268576  max resid 0.0002263534 
    ## ... Similar to previous best
    ## Run 159 stress 9.47612e-05 
    ## ... Procrustes: rmse 0.0001608575  max resid 0.0002569874 
    ## ... Similar to previous best
    ## Run 160 stress 9.317204e-05 
    ## ... Procrustes: rmse 0.0001670823  max resid 0.0002445854 
    ## ... Similar to previous best
    ## Run 161 stress 9.854706e-05 
    ## ... Procrustes: rmse 0.0001830254  max resid 0.0002394982 
    ## ... Similar to previous best
    ## Run 162 stress 9.721221e-05 
    ## ... Procrustes: rmse 0.0001509929  max resid 0.000229608 
    ## ... Similar to previous best
    ## Run 163 stress 8.897878e-05 
    ## ... Procrustes: rmse 0.0001114101  max resid 0.0002123775 
    ## ... Similar to previous best
    ## Run 164 stress 9.019353e-05 
    ## ... Procrustes: rmse 0.0001594519  max resid 0.0003022705 
    ## ... Similar to previous best
    ## Run 165 stress 9.941694e-05 
    ## ... Procrustes: rmse 0.0001860215  max resid 0.0003400244 
    ## ... Similar to previous best
    ## Run 166 stress 8.690887e-05 
    ## ... Procrustes: rmse 0.0001204856  max resid 0.0002377591 
    ## ... Similar to previous best
    ## Run 167 stress 9.973881e-05 
    ## ... Procrustes: rmse 0.0001525844  max resid 0.0003175537 
    ## ... Similar to previous best
    ## Run 168 stress 0.2280257 
    ## Run 169 stress 9.926431e-05 
    ## ... Procrustes: rmse 0.0001699915  max resid 0.0002694616 
    ## ... Similar to previous best
    ## Run 170 stress 9.785798e-05 
    ## ... Procrustes: rmse 0.0002129321  max resid 0.0003897873 
    ## ... Similar to previous best
    ## Run 171 stress 9.844599e-05 
    ## ... Procrustes: rmse 0.0001750825  max resid 0.0003048737 
    ## ... Similar to previous best
    ## Run 172 stress 9.756246e-05 
    ## ... Procrustes: rmse 0.0001791734  max resid 0.0003352188 
    ## ... Similar to previous best
    ## Run 173 stress 9.657865e-05 
    ## ... Procrustes: rmse 0.0002004122  max resid 0.0002823697 
    ## ... Similar to previous best
    ## Run 174 stress 9.611122e-05 
    ## ... Procrustes: rmse 0.0002005245  max resid 0.000372316 
    ## ... Similar to previous best
    ## Run 175 stress 9.887442e-05 
    ## ... Procrustes: rmse 0.0002037435  max resid 0.000375673 
    ## ... Similar to previous best
    ## Run 176 stress 9.381e-05 
    ## ... Procrustes: rmse 0.0001655642  max resid 0.0003075618 
    ## ... Similar to previous best
    ## Run 177 stress 9.54131e-05 
    ## ... Procrustes: rmse 0.0001252963  max resid 0.000233938 
    ## ... Similar to previous best
    ## Run 178 stress 9.024593e-05 
    ## ... Procrustes: rmse 0.0001273049  max resid 0.0002837606 
    ## ... Similar to previous best
    ## Run 179 stress 9.395019e-05 
    ## ... Procrustes: rmse 0.0002012256  max resid 0.0003761731 
    ## ... Similar to previous best
    ## Run 180 stress 9.660724e-05 
    ## ... Procrustes: rmse 0.0001880881  max resid 0.0003341627 
    ## ... Similar to previous best
    ## Run 181 stress 9.757437e-05 
    ## ... Procrustes: rmse 0.0002102379  max resid 0.0003849321 
    ## ... Similar to previous best
    ## Run 182 stress 9.882552e-05 
    ## ... Procrustes: rmse 0.0002000543  max resid 0.0003496368 
    ## ... Similar to previous best
    ## Run 183 stress 8.913489e-05 
    ## ... Procrustes: rmse 9.878227e-05  max resid 0.0001843246 
    ## ... Similar to previous best
    ## Run 184 stress 9.860531e-05 
    ## ... Procrustes: rmse 0.0001989849  max resid 0.0003708232 
    ## ... Similar to previous best
    ## Run 185 stress 9.238896e-05 
    ## ... Procrustes: rmse 0.0001107015  max resid 0.000203744 
    ## ... Similar to previous best
    ## Run 186 stress 8.993914e-05 
    ## ... Procrustes: rmse 0.0001840258  max resid 0.0003523352 
    ## ... Similar to previous best
    ## Run 187 stress 9.133577e-05 
    ## ... Procrustes: rmse 0.0001236647  max resid 0.0002395561 
    ## ... Similar to previous best
    ## Run 188 stress 9.822757e-05 
    ## ... Procrustes: rmse 0.0001907148  max resid 0.0002773026 
    ## ... Similar to previous best
    ## Run 189 stress 9.852968e-05 
    ## ... Procrustes: rmse 0.0001995117  max resid 0.0003688511 
    ## ... Similar to previous best
    ## Run 190 stress 9.472917e-05 
    ## ... Procrustes: rmse 0.0001724549  max resid 0.000243672 
    ## ... Similar to previous best
    ## Run 191 stress 9.641935e-05 
    ## ... Procrustes: rmse 0.0001634631  max resid 0.0002595856 
    ## ... Similar to previous best
    ## Run 192 stress 9.668836e-05 
    ## ... Procrustes: rmse 0.000157105  max resid 0.0002580225 
    ## ... Similar to previous best
    ## Run 193 stress 9.761464e-05 
    ## ... Procrustes: rmse 0.0002100691  max resid 0.0003859451 
    ## ... Similar to previous best
    ## Run 194 stress 9.606347e-05 
    ## ... Procrustes: rmse 0.0001586104  max resid 0.000273085 
    ## ... Similar to previous best
    ## Run 195 stress 9.636921e-05 
    ## ... Procrustes: rmse 0.0002067508  max resid 0.0003800469 
    ## ... Similar to previous best
    ## Run 196 stress 9.794862e-05 
    ## ... Procrustes: rmse 0.0001770085  max resid 0.000248099 
    ## ... Similar to previous best
    ## Run 197 stress 9.478307e-05 
    ## ... Procrustes: rmse 0.0001627131  max resid 0.0002598287 
    ## ... Similar to previous best
    ## Run 198 stress 0.2507977 
    ## Run 199 stress 9.541604e-05 
    ## ... Procrustes: rmse 0.0001888802  max resid 0.0002660189 
    ## ... Similar to previous best
    ## Run 200 stress 9.707888e-05 
    ## ... Procrustes: rmse 0.0001252956  max resid 0.0001979434 
    ## ... Similar to previous best
    ## Run 201 stress 9.760589e-05 
    ## ... Procrustes: rmse 0.0001554547  max resid 0.0002589395 
    ## ... Similar to previous best
    ## Run 202 stress 9.419843e-05 
    ## ... Procrustes: rmse 0.0001800663  max resid 0.000308653 
    ## ... Similar to previous best
    ## Run 203 stress 9.693898e-05 
    ## ... Procrustes: rmse 0.0001639346  max resid 0.000309701 
    ## ... Similar to previous best
    ## Run 204 stress 9.619169e-05 
    ## ... Procrustes: rmse 0.000168643  max resid 0.0003083635 
    ## ... Similar to previous best
    ## Run 205 stress 9.666654e-05 
    ## ... Procrustes: rmse 0.0001658119  max resid 0.0002669096 
    ## ... Similar to previous best
    ## Run 206 stress 9.521109e-05 
    ## ... Procrustes: rmse 0.0001689302  max resid 0.000313374 
    ## ... Similar to previous best
    ## Run 207 stress 9.849785e-05 
    ## ... Procrustes: rmse 0.0001963754  max resid 0.0003443182 
    ## ... Similar to previous best
    ## Run 208 stress 9.540466e-05 
    ## ... Procrustes: rmse 0.0001737861  max resid 0.0002280271 
    ## ... Similar to previous best
    ## Run 209 stress 8.045203e-05 
    ## ... Procrustes: rmse 0.0001213385  max resid 0.0002660776 
    ## ... Similar to previous best
    ## Run 210 stress 9.890165e-05 
    ## ... Procrustes: rmse 0.0001739867  max resid 0.0002389377 
    ## ... Similar to previous best
    ## Run 211 stress 9.689012e-05 
    ## ... Procrustes: rmse 0.0001304729  max resid 0.000234282 
    ## ... Similar to previous best
    ## Run 212 stress 9.749385e-05 
    ## ... Procrustes: rmse 0.0001631014  max resid 0.0002539036 
    ## ... Similar to previous best
    ## Run 213 stress 9.201572e-05 
    ## ... Procrustes: rmse 0.0001863897  max resid 0.0002723916 
    ## ... Similar to previous best
    ## Run 214 stress 9.804194e-05 
    ## ... Procrustes: rmse 0.0001456118  max resid 0.0002439135 
    ## ... Similar to previous best
    ## Run 215 stress 9.677403e-05 
    ## ... Procrustes: rmse 0.0001691855  max resid 0.0002339501 
    ## ... Similar to previous best
    ## Run 216 stress 9.695093e-05 
    ## ... Procrustes: rmse 0.0001605371  max resid 0.0002249435 
    ## ... Similar to previous best
    ## Run 217 stress 9.874638e-05 
    ## ... Procrustes: rmse 0.0001385135  max resid 0.0002465165 
    ## ... Similar to previous best
    ## Run 218 stress 9.877414e-05 
    ## ... Procrustes: rmse 0.000194465  max resid 0.0002785418 
    ## ... Similar to previous best
    ## Run 219 stress 9.946893e-05 
    ## ... Procrustes: rmse 0.0001676888  max resid 0.0002621806 
    ## ... Similar to previous best
    ## Run 220 stress 8.68334e-05 
    ## ... Procrustes: rmse 0.0001415407  max resid 0.0002368904 
    ## ... Similar to previous best
    ## Run 221 stress 8.822364e-05 
    ## ... Procrustes: rmse 0.0001354282  max resid 0.0002560702 
    ## ... Similar to previous best
    ## Run 222 stress 8.900966e-05 
    ## ... Procrustes: rmse 0.0001045191  max resid 0.0001990698 
    ## ... Similar to previous best
    ## Run 223 stress 9.027548e-05 
    ## ... Procrustes: rmse 0.0001683999  max resid 0.0003205483 
    ## ... Similar to previous best
    ## Run 224 stress 9.994161e-05 
    ## ... Procrustes: rmse 0.0002093534  max resid 0.0002917303 
    ## ... Similar to previous best
    ## Run 225 stress 9.554381e-05 
    ## ... Procrustes: rmse 0.000162452  max resid 0.0002853369 
    ## ... Similar to previous best
    ## Run 226 stress 9.694995e-05 
    ## ... Procrustes: rmse 0.0001624176  max resid 0.0002609297 
    ## ... Similar to previous best
    ## Run 227 stress 8.806528e-05 
    ## ... Procrustes: rmse 0.0001616566  max resid 0.0003138267 
    ## ... Similar to previous best
    ## Run 228 stress 9.683707e-05 
    ## ... Procrustes: rmse 0.0001407792  max resid 0.0002215822 
    ## ... Similar to previous best
    ## Run 229 stress 9.60336e-05 
    ## ... Procrustes: rmse 0.0001362021  max resid 0.0002477559 
    ## ... Similar to previous best
    ## Run 230 stress 9.947853e-05 
    ## ... Procrustes: rmse 0.0001350509  max resid 0.0002427468 
    ## ... Similar to previous best
    ## Run 231 stress 9.08999e-05 
    ## ... Procrustes: rmse 0.0001077851  max resid 0.0001760347 
    ## ... Similar to previous best
    ## Run 232 stress 7.920435e-05 
    ## ... Procrustes: rmse 0.0001544603  max resid 0.0002987643 
    ## ... Similar to previous best
    ## Run 233 stress 9.319996e-05 
    ## ... Procrustes: rmse 0.0001617223  max resid 0.0002579155 
    ## ... Similar to previous best
    ## Run 234 stress 9.656259e-05 
    ## ... Procrustes: rmse 0.0001636635  max resid 0.0002660066 
    ## ... Similar to previous best
    ## Run 235 stress 9.732297e-05 
    ## ... Procrustes: rmse 0.000152467  max resid 0.0002442395 
    ## ... Similar to previous best
    ## Run 236 stress 9.730338e-05 
    ## ... Procrustes: rmse 0.0001986448  max resid 0.0002834441 
    ## ... Similar to previous best
    ## Run 237 stress 9.3728e-05 
    ## ... Procrustes: rmse 0.0001436906  max resid 0.0002133773 
    ## ... Similar to previous best
    ## Run 238 stress 9.536232e-05 
    ## ... Procrustes: rmse 0.0001305428  max resid 0.0002356424 
    ## ... Similar to previous best
    ## Run 239 stress 9.243017e-05 
    ## ... Procrustes: rmse 0.0001565496  max resid 0.00025549 
    ## ... Similar to previous best
    ## Run 240 stress 9.833565e-05 
    ## ... Procrustes: rmse 0.0001837643  max resid 0.000337832 
    ## ... Similar to previous best
    ## Run 241 stress 9.010557e-05 
    ## ... Procrustes: rmse 0.000135345  max resid 0.0002630927 
    ## ... Similar to previous best
    ## Run 242 stress 9.926844e-05 
    ## ... Procrustes: rmse 0.0001607055  max resid 0.000270634 
    ## ... Similar to previous best
    ## Run 243 stress 9.291358e-05 
    ## ... Procrustes: rmse 0.0001882049  max resid 0.0003481151 
    ## ... Similar to previous best
    ## Run 244 stress 9.896175e-05 
    ## ... Procrustes: rmse 0.0001740694  max resid 0.0002728651 
    ## ... Similar to previous best
    ## Run 245 stress 9.675863e-05 
    ## ... Procrustes: rmse 0.0001207469  max resid 0.0002086071 
    ## ... Similar to previous best
    ## Run 246 stress 8.92708e-05 
    ## ... Procrustes: rmse 0.000116775  max resid 0.0002093668 
    ## ... Similar to previous best
    ## Run 247 stress 9.318638e-05 
    ## ... Procrustes: rmse 0.0001874295  max resid 0.0003294368 
    ## ... Similar to previous best
    ## Run 248 stress 9.263304e-05 
    ## ... Procrustes: rmse 0.0001496202  max resid 0.0002354438 
    ## ... Similar to previous best
    ## Run 249 stress 8.99625e-05 
    ## ... Procrustes: rmse 0.0001470308  max resid 0.0002372106 
    ## ... Similar to previous best
    ## Run 250 stress 9.677103e-05 
    ## ... Procrustes: rmse 0.0001790817  max resid 0.0002327681 
    ## ... Similar to previous best
    ## Run 251 stress 9.829873e-05 
    ## ... Procrustes: rmse 0.0001343784  max resid 0.0002268344 
    ## ... Similar to previous best
    ## Run 252 stress 9.405618e-05 
    ## ... Procrustes: rmse 0.0002038834  max resid 0.0002863515 
    ## ... Similar to previous best
    ## Run 253 stress 9.93115e-05 
    ## ... Procrustes: rmse 0.0001329026  max resid 0.0002375385 
    ## ... Similar to previous best
    ## Run 254 stress 9.047433e-05 
    ## ... Procrustes: rmse 0.0001554274  max resid 0.0002560541 
    ## ... Similar to previous best
    ## Run 255 stress 9.102898e-05 
    ## ... Procrustes: rmse 0.000145644  max resid 0.0002071436 
    ## ... Similar to previous best
    ## Run 256 stress 8.567028e-05 
    ## ... Procrustes: rmse 0.0001341635  max resid 0.0002156395 
    ## ... Similar to previous best
    ## Run 257 stress 9.164455e-05 
    ## ... Procrustes: rmse 0.0001311889  max resid 0.0001973924 
    ## ... Similar to previous best
    ## Run 258 stress 9.253546e-05 
    ## ... Procrustes: rmse 0.0001866425  max resid 0.0002755812 
    ## ... Similar to previous best
    ## Run 259 stress 9.485006e-05 
    ## ... Procrustes: rmse 0.0001916048  max resid 0.0003373832 
    ## ... Similar to previous best
    ## Run 260 stress 9.665085e-05 
    ## ... Procrustes: rmse 0.0001745698  max resid 0.0003149743 
    ## ... Similar to previous best
    ## Run 261 stress 9.951476e-05 
    ## ... Procrustes: rmse 0.0001739992  max resid 0.0002709006 
    ## ... Similar to previous best
    ## Run 262 stress 9.455726e-05 
    ## ... Procrustes: rmse 0.0001651542  max resid 0.0002650735 
    ## ... Similar to previous best
    ## Run 263 stress 9.925873e-05 
    ## ... Procrustes: rmse 0.0001686105  max resid 0.0002596608 
    ## ... Similar to previous best
    ## Run 264 stress 9.414925e-05 
    ## ... Procrustes: rmse 0.000149544  max resid 0.0002636327 
    ## ... Similar to previous best
    ## Run 265 stress 9.279863e-05 
    ## ... Procrustes: rmse 0.0002010607  max resid 0.0003742282 
    ## ... Similar to previous best
    ## Run 266 stress 9.833623e-05 
    ## ... Procrustes: rmse 0.0001980478  max resid 0.0002854035 
    ## ... Similar to previous best
    ## Run 267 stress 9.268982e-05 
    ## ... Procrustes: rmse 0.0001701231  max resid 0.0002233956 
    ## ... Similar to previous best
    ## Run 268 stress 0.2203301 
    ## Run 269 stress 9.920846e-05 
    ## ... Procrustes: rmse 0.0002038372  max resid 0.0002909314 
    ## ... Similar to previous best
    ## Run 270 stress 9.064101e-05 
    ## ... Procrustes: rmse 0.0001971869  max resid 0.0003711112 
    ## ... Similar to previous best
    ## Run 271 stress 9.977763e-05 
    ## ... Procrustes: rmse 0.0002003971  max resid 0.0002880908 
    ## ... Similar to previous best
    ## Run 272 stress 9.678059e-05 
    ## ... Procrustes: rmse 0.0001701727  max resid 0.0003217613 
    ## ... Similar to previous best
    ## Run 273 stress 9.831389e-05 
    ## ... Procrustes: rmse 0.0002048624  max resid 0.00028687 
    ## ... Similar to previous best
    ## Run 274 stress 8.657499e-05 
    ## ... Procrustes: rmse 0.0001413534  max resid 0.0002417736 
    ## ... Similar to previous best
    ## Run 275 stress 8.315228e-05 
    ## ... Procrustes: rmse 0.0001475786  max resid 0.0002876681 
    ## ... Similar to previous best
    ## Run 276 stress 9.255275e-05 
    ## ... Procrustes: rmse 0.0001911  max resid 0.0003616881 
    ## ... Similar to previous best
    ## Run 277 stress 0.206342 
    ## Run 278 stress 8.697154e-05 
    ## ... Procrustes: rmse 0.0001129338  max resid 0.0002015308 
    ## ... Similar to previous best
    ## Run 279 stress 9.83845e-05 
    ## ... Procrustes: rmse 0.0001286795  max resid 0.000221981 
    ## ... Similar to previous best
    ## Run 280 stress 9.817973e-05 
    ## ... Procrustes: rmse 0.0001659921  max resid 0.000265604 
    ## ... Similar to previous best
    ## Run 281 stress 9.759371e-05 
    ## ... Procrustes: rmse 0.0001348009  max resid 0.000239897 
    ## ... Similar to previous best
    ## Run 282 stress 9.39767e-05 
    ## ... Procrustes: rmse 0.000159041  max resid 0.0002587323 
    ## ... Similar to previous best
    ## Run 283 stress 9.665909e-05 
    ## ... Procrustes: rmse 0.0001857852  max resid 0.0003267442 
    ## ... Similar to previous best
    ## Run 284 stress 9.711954e-05 
    ## ... Procrustes: rmse 0.000160228  max resid 0.0002650471 
    ## ... Similar to previous best
    ## Run 285 stress 9.05552e-05 
    ## ... Procrustes: rmse 0.0001210946  max resid 0.0002234693 
    ## ... Similar to previous best
    ## Run 286 stress 9.269434e-05 
    ## ... Procrustes: rmse 0.0001793018  max resid 0.000322609 
    ## ... Similar to previous best
    ## Run 287 stress 8.948097e-05 
    ## ... Procrustes: rmse 0.0001707439  max resid 0.0003126971 
    ## ... Similar to previous best
    ## Run 288 stress 9.779844e-05 
    ## ... Procrustes: rmse 0.0001394119  max resid 0.0002143824 
    ## ... Similar to previous best
    ## Run 289 stress 9.622258e-05 
    ## ... Procrustes: rmse 0.0001578034  max resid 0.0002171058 
    ## ... Similar to previous best
    ## Run 290 stress 8.290468e-05 
    ## ... Procrustes: rmse 0.0001318636  max resid 0.000230553 
    ## ... Similar to previous best
    ## Run 291 stress 9.654733e-05 
    ## ... Procrustes: rmse 0.0001227142  max resid 0.0002058598 
    ## ... Similar to previous best
    ## Run 292 stress 9.895427e-05 
    ## ... Procrustes: rmse 0.0002108653  max resid 0.0003846706 
    ## ... Similar to previous best
    ## Run 293 stress 9.889646e-05 
    ## ... Procrustes: rmse 0.000138853  max resid 0.0002476143 
    ## ... Similar to previous best
    ## Run 294 stress 9.912168e-05 
    ## ... Procrustes: rmse 0.0001711969  max resid 0.0002682232 
    ## ... Similar to previous best
    ## Run 295 stress 9.980394e-05 
    ## ... Procrustes: rmse 0.0001947085  max resid 0.0003625205 
    ## ... Similar to previous best
    ## Run 296 stress 9.8192e-05 
    ## ... Procrustes: rmse 0.0001718481  max resid 0.0002691312 
    ## ... Similar to previous best
    ## Run 297 stress 9.858042e-05 
    ## ... Procrustes: rmse 0.0001701426  max resid 0.0002688706 
    ## ... Similar to previous best
    ## Run 298 stress 9.830497e-05 
    ## ... Procrustes: rmse 0.000132427  max resid 0.0002253351 
    ## ... Similar to previous best
    ## Run 299 stress 0.2413829 
    ## Run 300 stress 9.875188e-05 
    ## ... Procrustes: rmse 0.0001868926  max resid 0.0003273801 
    ## ... Similar to previous best
    ## Run 301 stress 8.892244e-05 
    ## ... Procrustes: rmse 0.0001164928  max resid 0.000209548 
    ## ... Similar to previous best
    ## Run 302 stress 9.575313e-05 
    ## ... Procrustes: rmse 0.0002078517  max resid 0.0003821912 
    ## ... Similar to previous best
    ## Run 303 stress 9.860378e-05 
    ## ... Procrustes: rmse 0.0002002774  max resid 0.0003721437 
    ## ... Similar to previous best
    ## Run 304 stress 9.723302e-05 
    ## ... Procrustes: rmse 0.0001706389  max resid 0.0002708278 
    ## ... Similar to previous best
    ## Run 305 stress 8.683752e-05 
    ## ... Procrustes: rmse 0.0001551865  max resid 0.0002979449 
    ## ... Similar to previous best
    ## Run 306 stress 9.582707e-05 
    ## ... Procrustes: rmse 0.0002069494  max resid 0.0002901624 
    ## ... Similar to previous best
    ## Run 307 stress 9.821051e-05 
    ## ... Procrustes: rmse 0.0002010854  max resid 0.0002889049 
    ## ... Similar to previous best
    ## Run 308 stress 9.925018e-05 
    ## ... Procrustes: rmse 0.0001602443  max resid 0.000235378 
    ## ... Similar to previous best
    ## Run 309 stress 9.545041e-05 
    ## ... Procrustes: rmse 0.0001308903  max resid 0.0002350079 
    ## ... Similar to previous best
    ## Run 310 stress 9.587696e-05 
    ## ... Procrustes: rmse 0.0001617533  max resid 0.0002544262 
    ## ... Similar to previous best
    ## Run 311 stress 8.331105e-05 
    ## ... Procrustes: rmse 0.0001327586  max resid 0.0002207541 
    ## ... Similar to previous best
    ## Run 312 stress 9.957691e-05 
    ## ... Procrustes: rmse 0.0001998195  max resid 0.0003749569 
    ## ... Similar to previous best
    ## Run 313 stress 9.767953e-05 
    ## ... Procrustes: rmse 9.920977e-05  max resid 0.0002233438 
    ## ... Similar to previous best
    ## Run 314 stress 9.846704e-05 
    ## ... Procrustes: rmse 0.0001705105  max resid 0.0002661562 
    ## ... Similar to previous best
    ## Run 315 stress 9.756597e-05 
    ## ... Procrustes: rmse 0.0002024284  max resid 0.0002847871 
    ## ... Similar to previous best
    ## Run 316 stress 9.838183e-05 
    ## ... Procrustes: rmse 0.0001549507  max resid 0.000233972 
    ## ... Similar to previous best
    ## Run 317 stress 9.988516e-05 
    ## ... Procrustes: rmse 0.0001760429  max resid 0.0002726174 
    ## ... Similar to previous best
    ## Run 318 stress 9.580908e-05 
    ## ... Procrustes: rmse 0.0001666801  max resid 0.0002381244 
    ## ... Similar to previous best
    ## Run 319 stress 9.042455e-05 
    ## ... Procrustes: rmse 0.0001865226  max resid 0.0003532961 
    ## ... Similar to previous best
    ## Run 320 stress 9.856628e-05 
    ## ... Procrustes: rmse 0.0001611591  max resid 0.0002620088 
    ## ... Similar to previous best
    ## Run 321 stress 9.691951e-05 
    ## ... Procrustes: rmse 0.0001314004  max resid 0.0002223519 
    ## ... Similar to previous best
    ## Run 322 stress 9.965973e-05 
    ## ... Procrustes: rmse 0.0002052476  max resid 0.0002921847 
    ## ... Similar to previous best
    ## Run 323 stress 9.733166e-05 
    ## ... Procrustes: rmse 0.0001681567  max resid 0.0002360427 
    ## ... Similar to previous best
    ## Run 324 stress 9.378296e-05 
    ## ... Procrustes: rmse 0.0001573582  max resid 0.0002532391 
    ## ... Similar to previous best
    ## Run 325 stress 9.350665e-05 
    ## ... Procrustes: rmse 0.0001359011  max resid 0.0002238785 
    ## ... Similar to previous best
    ## Run 326 stress 9.895181e-05 
    ## ... Procrustes: rmse 0.0001716321  max resid 0.000268728 
    ## ... Similar to previous best
    ## Run 327 stress 9.412485e-05 
    ## ... Procrustes: rmse 0.0001387016  max resid 0.000229412 
    ## ... Similar to previous best
    ## Run 328 stress 9.737049e-05 
    ## ... Procrustes: rmse 0.0001293103  max resid 0.000214152 
    ## ... Similar to previous best
    ## Run 329 stress 8.907798e-05 
    ## ... Procrustes: rmse 0.0001600007  max resid 0.0002376904 
    ## ... Similar to previous best
    ## Run 330 stress 9.748706e-05 
    ## ... Procrustes: rmse 0.0001925933  max resid 0.0002692137 
    ## ... Similar to previous best
    ## Run 331 stress 9.976721e-05 
    ## ... Procrustes: rmse 0.0001672987  max resid 0.0002687606 
    ## ... Similar to previous best
    ## Run 332 stress 9.72098e-05 
    ## ... Procrustes: rmse 0.0001449736  max resid 0.0002189577 
    ## ... Similar to previous best
    ## Run 333 stress 9.817505e-05 
    ## ... Procrustes: rmse 0.0001733663  max resid 0.0003163883 
    ## ... Similar to previous best
    ## Run 334 stress 9.85543e-05 
    ## ... Procrustes: rmse 0.0001520782  max resid 0.0003239142 
    ## ... Similar to previous best
    ## Run 335 stress 9.866775e-05 
    ## ... Procrustes: rmse 0.0001647128  max resid 0.000227884 
    ## ... Similar to previous best
    ## Run 336 stress 9.41349e-05 
    ## ... Procrustes: rmse 0.0001816096  max resid 0.0003394585 
    ## ... Similar to previous best
    ## Run 337 stress 9.530558e-05 
    ## ... Procrustes: rmse 0.0001657895  max resid 0.0002626424 
    ## ... Similar to previous best
    ## Run 338 stress 9.423563e-05 
    ## ... Procrustes: rmse 0.0001774618  max resid 0.0003177879 
    ## ... Similar to previous best
    ## Run 339 stress 9.822573e-05 
    ## ... Procrustes: rmse 9.144533e-05  max resid 0.0001463863 
    ## ... Similar to previous best
    ## Run 340 stress 9.392505e-05 
    ## ... Procrustes: rmse 0.0001729116  max resid 0.0003227375 
    ## ... Similar to previous best
    ## Run 341 stress 9.381955e-05 
    ## ... Procrustes: rmse 0.0001615555  max resid 0.0002222572 
    ## ... Similar to previous best
    ## Run 342 stress 8.904024e-05 
    ## ... Procrustes: rmse 0.0001661046  max resid 0.0003155753 
    ## ... Similar to previous best
    ## Run 343 stress 9.282913e-05 
    ## ... Procrustes: rmse 0.000174003  max resid 0.0003196779 
    ## ... Similar to previous best
    ## Run 344 stress 9.926711e-05 
    ## ... Procrustes: rmse 0.0001660365  max resid 0.0002336767 
    ## ... Similar to previous best
    ## Run 345 stress 9.493684e-05 
    ## ... Procrustes: rmse 0.0001564835  max resid 0.0002364375 
    ## ... Similar to previous best
    ## Run 346 stress 9.695467e-05 
    ## ... Procrustes: rmse 0.0001797892  max resid 0.0002349475 
    ## ... Similar to previous best
    ## Run 347 stress 9.739209e-05 
    ## ... Procrustes: rmse 0.0001908767  max resid 0.0002738015 
    ## ... Similar to previous best
    ## Run 348 stress 8.726958e-05 
    ## ... Procrustes: rmse 0.0001053387  max resid 0.0001885569 
    ## ... Similar to previous best
    ## Run 349 stress 9.914539e-05 
    ## ... Procrustes: rmse 0.0001893892  max resid 0.0003312401 
    ## ... Similar to previous best
    ## Run 350 stress 9.468307e-05 
    ## ... Procrustes: rmse 0.0001917606  max resid 0.0003558581 
    ## ... Similar to previous best
    ## Run 351 stress 9.04949e-05 
    ## ... Procrustes: rmse 0.0001833438  max resid 0.0002720647 
    ## ... Similar to previous best
    ## Run 352 stress 9.302407e-05 
    ## ... Procrustes: rmse 0.0001897021  max resid 0.0003597972 
    ## ... Similar to previous best
    ## Run 353 stress 8.367648e-05 
    ## ... Procrustes: rmse 0.0001049946  max resid 0.0001964054 
    ## ... Similar to previous best
    ## Run 354 stress 9.816511e-05 
    ## ... Procrustes: rmse 0.0002000882  max resid 0.0002842982 
    ## ... Similar to previous best
    ## Run 355 stress 9.834197e-05 
    ## ... Procrustes: rmse 0.0001933943  max resid 0.0002701956 
    ## ... Similar to previous best
    ## Run 356 stress 9.866024e-05 
    ## ... Procrustes: rmse 0.000200011  max resid 0.0003705295 
    ## ... Similar to previous best
    ## Run 357 stress 9.363058e-05 
    ## ... Procrustes: rmse 0.0001937121  max resid 0.0003684377 
    ## ... Similar to previous best
    ## Run 358 stress 9.361138e-05 
    ## ... Procrustes: rmse 0.0001561754  max resid 0.0002551252 
    ## ... Similar to previous best
    ## Run 359 stress 9.195608e-05 
    ## ... Procrustes: rmse 0.0001695659  max resid 0.0003417633 
    ## ... Similar to previous best
    ## Run 360 stress 9.996493e-05 
    ## ... Procrustes: rmse 0.0001405876  max resid 0.0002488634 
    ## ... Similar to previous best
    ## Run 361 stress 9.371368e-05 
    ## ... Procrustes: rmse 0.0001904467  max resid 0.0002784555 
    ## ... Similar to previous best
    ## Run 362 stress 9.715489e-05 
    ## ... Procrustes: rmse 0.0001870129  max resid 0.0003284672 
    ## ... Similar to previous best
    ## Run 363 stress 9.707038e-05 
    ## ... Procrustes: rmse 0.0001870795  max resid 0.000329583 
    ## ... Similar to previous best
    ## Run 364 stress 9.6712e-05 
    ## ... Procrustes: rmse 0.000158308  max resid 0.0002248616 
    ## ... Similar to previous best
    ## Run 365 stress 8.312804e-05 
    ## ... Procrustes: rmse 0.0001234252  max resid 0.0002282461 
    ## ... Similar to previous best
    ## Run 366 stress 9.221369e-05 
    ## ... Procrustes: rmse 0.0001398176  max resid 0.0002339126 
    ## ... Similar to previous best
    ## Run 367 stress 9.613662e-05 
    ## ... Procrustes: rmse 0.0001612894  max resid 0.0002491526 
    ## ... Similar to previous best
    ## Run 368 stress 9.617881e-05 
    ## ... Procrustes: rmse 0.0001887746  max resid 0.0002712722 
    ## ... Similar to previous best
    ## Run 369 stress 9.449571e-05 
    ## ... Procrustes: rmse 0.0001739945  max resid 0.0002279874 
    ## ... Similar to previous best
    ## Run 370 stress 9.485541e-05 
    ## ... Procrustes: rmse 0.0001710892  max resid 0.0003276026 
    ## ... Similar to previous best
    ## Run 371 stress 9.983637e-05 
    ## ... Procrustes: rmse 0.0001702354  max resid 0.0002643282 
    ## ... Similar to previous best
    ## Run 372 stress 9.610222e-05 
    ## ... Procrustes: rmse 0.0001157559  max resid 0.0001725057 
    ## ... Similar to previous best
    ## Run 373 stress 9.401659e-05 
    ## ... Procrustes: rmse 0.0001575176  max resid 0.0002488782 
    ## ... Similar to previous best
    ## Run 374 stress 9.894066e-05 
    ## ... Procrustes: rmse 0.0002063901  max resid 0.0002885629 
    ## ... Similar to previous best
    ## Run 375 stress 9.385635e-05 
    ## ... Procrustes: rmse 0.0001603652  max resid 0.0002608172 
    ## ... Similar to previous best
    ## Run 376 stress 9.561122e-05 
    ## ... Procrustes: rmse 0.0001720045  max resid 0.0002364597 
    ## ... Similar to previous best
    ## Run 377 stress 9.466608e-05 
    ## ... Procrustes: rmse 0.0001652565  max resid 0.0002418491 
    ## ... Similar to previous best
    ## Run 378 stress 8.672714e-05 
    ## ... Procrustes: rmse 0.0001277316  max resid 0.0002695029 
    ## ... Similar to previous best
    ## Run 379 stress 9.731837e-05 
    ## ... Procrustes: rmse 0.0002008478  max resid 0.0003711038 
    ## ... Similar to previous best
    ## Run 380 stress 8.283797e-05 
    ## ... Procrustes: rmse 0.0001488884  max resid 0.0002732432 
    ## ... Similar to previous best
    ## Run 381 stress 9.417448e-05 
    ## ... Procrustes: rmse 0.0001720361  max resid 0.0002281304 
    ## ... Similar to previous best
    ## Run 382 stress 9.155365e-05 
    ## ... Procrustes: rmse 0.0001235296  max resid 0.0002241673 
    ## ... Similar to previous best
    ## Run 383 stress 9.040165e-05 
    ## ... Procrustes: rmse 0.0001428014  max resid 0.0002387404 
    ## ... Similar to previous best
    ## Run 384 stress 8.756264e-05 
    ## ... Procrustes: rmse 0.0001465959  max resid 0.0002500481 
    ## ... Similar to previous best
    ## Run 385 stress 9.891256e-05 
    ## ... Procrustes: rmse 0.0001640454  max resid 0.0002771608 
    ## ... Similar to previous best
    ## Run 386 stress 9.495106e-05 
    ## ... Procrustes: rmse 0.0001923189  max resid 0.0003391866 
    ## ... Similar to previous best
    ## Run 387 stress 9.356774e-05 
    ## ... Procrustes: rmse 0.0001086112  max resid 0.0001889994 
    ## ... Similar to previous best
    ## Run 388 stress 9.593969e-05 
    ## ... Procrustes: rmse 0.0001318733  max resid 0.0002716734 
    ## ... Similar to previous best
    ## Run 389 stress 9.579496e-05 
    ## ... Procrustes: rmse 0.0001320386  max resid 0.00023533 
    ## ... Similar to previous best
    ## Run 390 stress 9.254791e-05 
    ## ... Procrustes: rmse 0.0001660248  max resid 0.0003066974 
    ## ... Similar to previous best
    ## Run 391 stress 9.619061e-05 
    ## ... Procrustes: rmse 0.0001909887  max resid 0.0003290116 
    ## ... Similar to previous best
    ## Run 392 stress 9.654579e-05 
    ## ... Procrustes: rmse 0.0001510968  max resid 0.0002130958 
    ## ... Similar to previous best
    ## Run 393 stress 9.433824e-05 
    ## ... Procrustes: rmse 0.0001953132  max resid 0.0003395202 
    ## ... Similar to previous best
    ## Run 394 stress 9.478906e-05 
    ## ... Procrustes: rmse 0.0001987791  max resid 0.0002809303 
    ## ... Similar to previous best
    ## Run 395 stress 9.325704e-05 
    ## ... Procrustes: rmse 0.0001505005  max resid 0.0002404923 
    ## ... Similar to previous best
    ## Run 396 stress 9.855896e-05 
    ## ... Procrustes: rmse 0.0002144955  max resid 0.0003919974 
    ## ... Similar to previous best
    ## Run 397 stress 9.871234e-05 
    ## ... Procrustes: rmse 0.0002003418  max resid 0.0002857502 
    ## ... Similar to previous best
    ## Run 398 stress 9.468988e-05 
    ## ... Procrustes: rmse 0.0001256156  max resid 0.0002301119 
    ## ... Similar to previous best
    ## Run 399 stress 9.997593e-05 
    ## ... Procrustes: rmse 0.0001974034  max resid 0.0003496691 
    ## ... Similar to previous best
    ## Run 400 stress 9.960995e-05 
    ## ... Procrustes: rmse 0.0002032367  max resid 0.0002917695 
    ## ... Similar to previous best
    ## Run 401 stress 8.814638e-05 
    ## ... Procrustes: rmse 0.0001594184  max resid 0.0002346676 
    ## ... Similar to previous best
    ## Run 402 stress 9.893349e-05 
    ## ... Procrustes: rmse 0.000182307  max resid 0.0003355017 
    ## ... Similar to previous best
    ## Run 403 stress 9.591108e-05 
    ## ... Procrustes: rmse 0.0001843929  max resid 0.000341195 
    ## ... Similar to previous best
    ## Run 404 stress 9.194073e-05 
    ## ... Procrustes: rmse 0.0001627455  max resid 0.0002508047 
    ## ... Similar to previous best
    ## Run 405 stress 9.654597e-05 
    ## ... Procrustes: rmse 0.0001144342  max resid 0.0002206355 
    ## ... Similar to previous best
    ## Run 406 stress 9.52889e-05 
    ## ... Procrustes: rmse 0.0001291004  max resid 0.0002319051 
    ## ... Similar to previous best
    ## Run 407 stress 9.940198e-05 
    ## ... Procrustes: rmse 0.0001870566  max resid 0.000339961 
    ## ... Similar to previous best
    ## Run 408 stress 9.724836e-05 
    ## ... Procrustes: rmse 0.0001684021  max resid 0.0002669762 
    ## ... Similar to previous best
    ## Run 409 stress 8.563868e-05 
    ## ... Procrustes: rmse 0.0001388771  max resid 0.0002345048 
    ## ... Similar to previous best
    ## Run 410 stress 9.184688e-05 
    ## ... Procrustes: rmse 0.0001718135  max resid 0.0003168248 
    ## ... Similar to previous best
    ## Run 411 stress 9.957984e-05 
    ## ... Procrustes: rmse 0.0001859018  max resid 0.0003411834 
    ## ... Similar to previous best
    ## Run 412 stress 9.736771e-05 
    ## ... Procrustes: rmse 0.000167471  max resid 0.0002300363 
    ## ... Similar to previous best
    ## Run 413 stress 9.89429e-05 
    ## ... Procrustes: rmse 0.0001560284  max resid 0.0002561536 
    ## ... Similar to previous best
    ## Run 414 stress 7.601375e-05 
    ## ... Procrustes: rmse 0.0001128684  max resid 0.0002125301 
    ## ... Similar to previous best
    ## Run 415 stress 9.960958e-05 
    ## ... Procrustes: rmse 0.000148922  max resid 0.0002446932 
    ## ... Similar to previous best
    ## Run 416 stress 0.2336195 
    ## Run 417 stress 9.843643e-05 
    ## ... Procrustes: rmse 0.0002007469  max resid 0.0002877264 
    ## ... Similar to previous best
    ## Run 418 stress 9.733631e-05 
    ## ... Procrustes: rmse 0.0002042278  max resid 0.0002857287 
    ## ... Similar to previous best
    ## Run 419 stress 9.406737e-05 
    ## ... Procrustes: rmse 0.0001781506  max resid 0.0003268338 
    ## ... Similar to previous best
    ## Run 420 stress 8.958328e-05 
    ## ... Procrustes: rmse 0.0001513448  max resid 0.0002478954 
    ## ... Similar to previous best
    ## Run 421 stress 9.833657e-05 
    ## ... Procrustes: rmse 0.0001783052  max resid 0.0003336652 
    ## ... Similar to previous best
    ## Run 422 stress 9.721933e-05 
    ## ... Procrustes: rmse 0.0001707734  max resid 0.0002294632 
    ## ... Similar to previous best
    ## Run 423 stress 9.454747e-05 
    ## ... Procrustes: rmse 0.0001845044  max resid 0.0002649394 
    ## ... Similar to previous best
    ## Run 424 stress 9.514893e-05 
    ## ... Procrustes: rmse 0.0001921402  max resid 0.000339152 
    ## ... Similar to previous best
    ## Run 425 stress 9.644372e-05 
    ## ... Procrustes: rmse 0.0001157479  max resid 0.0001905361 
    ## ... Similar to previous best
    ## Run 426 stress 9.704325e-05 
    ## ... Procrustes: rmse 0.0001366004  max resid 0.000294324 
    ## ... Similar to previous best
    ## Run 427 stress 9.560081e-05 
    ## ... Procrustes: rmse 0.0001585421  max resid 0.0002420554 
    ## ... Similar to previous best
    ## Run 428 stress 9.270638e-05 
    ## ... Procrustes: rmse 0.0001461871  max resid 0.0002372635 
    ## ... Similar to previous best
    ## Run 429 stress 9.897131e-05 
    ## ... Procrustes: rmse 0.000184021  max resid 0.0003382805 
    ## ... Similar to previous best
    ## Run 430 stress 9.546826e-05 
    ## ... Procrustes: rmse 0.0001903896  max resid 0.0003405856 
    ## ... Similar to previous best
    ## Run 431 stress 9.632434e-05 
    ## ... Procrustes: rmse 0.0001623031  max resid 0.0002753517 
    ## ... Similar to previous best
    ## Run 432 stress 9.45876e-05 
    ## ... Procrustes: rmse 0.0001535719  max resid 0.000251201 
    ## ... Similar to previous best
    ## Run 433 stress 9.544259e-05 
    ## ... Procrustes: rmse 0.0001803175  max resid 0.0003212906 
    ## ... Similar to previous best
    ## Run 434 stress 0.2332537 
    ## Run 435 stress 9.556467e-05 
    ## ... Procrustes: rmse 0.0001086594  max resid 0.0001783161 
    ## ... Similar to previous best
    ## Run 436 stress 9.229796e-05 
    ## ... Procrustes: rmse 0.000164161  max resid 0.000241926 
    ## ... Similar to previous best
    ## Run 437 stress 8.635701e-05 
    ## ... Procrustes: rmse 6.987772e-05  max resid 0.0001278772 
    ## ... Similar to previous best
    ## Run 438 stress 9.801165e-05 
    ## ... Procrustes: rmse 0.0001611484  max resid 0.0002368282 
    ## ... Similar to previous best
    ## Run 439 stress 9.55809e-05 
    ## ... Procrustes: rmse 0.0001658112  max resid 0.0002630104 
    ## ... Similar to previous best
    ## Run 440 stress 8.934585e-05 
    ## ... Procrustes: rmse 0.0001356525  max resid 0.0002010306 
    ## ... Similar to previous best
    ## Run 441 stress 9.883868e-05 
    ## ... Procrustes: rmse 0.0001384526  max resid 0.0002453608 
    ## ... Similar to previous best
    ## Run 442 stress 8.653367e-05 
    ## ... Procrustes: rmse 0.0001251737  max resid 0.0002484407 
    ## ... Similar to previous best
    ## Run 443 stress 9.619156e-05 
    ## ... Procrustes: rmse 0.0001972678  max resid 0.0003689439 
    ## ... Similar to previous best
    ## Run 444 stress 9.556538e-05 
    ## ... Procrustes: rmse 0.0001943983  max resid 0.0002797288 
    ## ... Similar to previous best
    ## Run 445 stress 9.693328e-05 
    ## ... Procrustes: rmse 0.0002011902  max resid 0.000283152 
    ## ... Similar to previous best
    ## Run 446 stress 9.660114e-05 
    ## ... Procrustes: rmse 0.0002025291  max resid 0.0002844719 
    ## ... Similar to previous best
    ## Run 447 stress 9.048737e-05 
    ## ... Procrustes: rmse 0.0001061293  max resid 0.0001805948 
    ## ... Similar to previous best
    ## Run 448 stress 9.470794e-05 
    ## ... Procrustes: rmse 0.0002007484  max resid 0.0002832627 
    ## ... Similar to previous best
    ## Run 449 stress 9.731775e-05 
    ## ... Procrustes: rmse 0.0001968931  max resid 0.0002823875 
    ## ... Similar to previous best
    ## Run 450 stress 9.788432e-05 
    ## ... Procrustes: rmse 0.00018387  max resid 0.000338552 
    ## ... Similar to previous best
    ## Run 451 stress 9.199516e-05 
    ## ... Procrustes: rmse 0.0001353823  max resid 0.0002954041 
    ## ... Similar to previous best
    ## Run 452 stress 8.670128e-05 
    ## ... Procrustes: rmse 8.359674e-05  max resid 0.0001482749 
    ## ... Similar to previous best
    ## Run 453 stress 9.746746e-05 
    ## ... Procrustes: rmse 0.0001709113  max resid 0.0002707074 
    ## ... Similar to previous best
    ## Run 454 stress 9.706076e-05 
    ## ... Procrustes: rmse 0.0001639548  max resid 0.000318925 
    ## ... Similar to previous best
    ## Run 455 stress 9.539989e-05 
    ## ... Procrustes: rmse 0.000150099  max resid 0.0002480393 
    ## ... Similar to previous best
    ## Run 456 stress 9.97976e-05 
    ## ... Procrustes: rmse 0.0001450338  max resid 0.0003057437 
    ## ... Similar to previous best
    ## Run 457 stress 9.908296e-05 
    ## ... Procrustes: rmse 0.000139002  max resid 0.0002477683 
    ## ... Similar to previous best
    ## Run 458 stress 9.80133e-05 
    ## ... Procrustes: rmse 0.0001706661  max resid 0.00026791 
    ## ... Similar to previous best
    ## Run 459 stress 9.964289e-05 
    ## ... Procrustes: rmse 0.0001861929  max resid 0.0003559155 
    ## ... Similar to previous best
    ## Run 460 stress 7.489213e-05 
    ## ... Procrustes: rmse 0.0001199095  max resid 0.0001814915 
    ## ... Similar to previous best
    ## Run 461 stress 9.815885e-05 
    ## ... Procrustes: rmse 0.0001666837  max resid 0.0002593423 
    ## ... Similar to previous best
    ## Run 462 stress 8.567316e-05 
    ## ... Procrustes: rmse 0.0001150888  max resid 0.0002448557 
    ## ... Similar to previous best
    ## Run 463 stress 9.668222e-05 
    ## ... Procrustes: rmse 0.0001814354  max resid 0.0003327563 
    ## ... Similar to previous best
    ## Run 464 stress 9.291184e-05 
    ## ... Procrustes: rmse 0.000138607  max resid 0.0002402569 
    ## ... Similar to previous best
    ## Run 465 stress 9.810035e-05 
    ## ... Procrustes: rmse 0.000213331  max resid 0.0003888524 
    ## ... Similar to previous best
    ## Run 466 stress 9.840329e-05 
    ## ... Procrustes: rmse 0.0001212614  max resid 0.0001893135 
    ## ... Similar to previous best
    ## Run 467 stress 7.729977e-05 
    ## ... Procrustes: rmse 0.0001496727  max resid 0.0002844003 
    ## ... Similar to previous best
    ## Run 468 stress 9.638382e-05 
    ## ... Procrustes: rmse 0.0001935107  max resid 0.0003400561 
    ## ... Similar to previous best
    ## Run 469 stress 9.573949e-05 
    ## ... Procrustes: rmse 0.0001962955  max resid 0.0003678693 
    ## ... Similar to previous best
    ## Run 470 stress 9.131711e-05 
    ## ... Procrustes: rmse 0.0001690103  max resid 0.0003071703 
    ## ... Similar to previous best
    ## Run 471 stress 9.746096e-05 
    ## ... Procrustes: rmse 0.0001956864  max resid 0.0003549209 
    ## ... Similar to previous best
    ## Run 472 stress 9.241141e-05 
    ## ... Procrustes: rmse 0.0001251547  max resid 0.0002238281 
    ## ... Similar to previous best
    ## Run 473 stress 9.532308e-05 
    ## ... Procrustes: rmse 0.0001743309  max resid 0.0003166881 
    ## ... Similar to previous best
    ## Run 474 stress 8.667493e-05 
    ## ... Procrustes: rmse 0.0001561216  max resid 0.0002146008 
    ## ... Similar to previous best
    ## Run 475 stress 9.964195e-05 
    ## ... Procrustes: rmse 0.0001362602  max resid 0.0002455906 
    ## ... Similar to previous best
    ## Run 476 stress 9.584715e-05 
    ## ... Procrustes: rmse 0.00016244  max resid 0.0002625209 
    ## ... Similar to previous best
    ## Run 477 stress 9.813186e-05 
    ## ... Procrustes: rmse 0.0001625932  max resid 0.0002303681 
    ## ... Similar to previous best
    ## Run 478 stress 9.823777e-05 
    ## ... Procrustes: rmse 0.0001994285  max resid 0.0003740157 
    ## ... Similar to previous best
    ## Run 479 stress 9.60227e-05 
    ## ... Procrustes: rmse 8.422263e-05  max resid 0.0001567275 
    ## ... Similar to previous best
    ## Run 480 stress 9.634659e-05 
    ## ... Procrustes: rmse 0.0001226423  max resid 0.000199872 
    ## ... Similar to previous best
    ## Run 481 stress 9.914085e-05 
    ## ... Procrustes: rmse 0.0001683222  max resid 0.0002674434 
    ## ... Similar to previous best
    ## Run 482 stress 9.402677e-05 
    ## ... Procrustes: rmse 0.0001752299  max resid 0.0003289191 
    ## ... Similar to previous best
    ## Run 483 stress 9.861011e-05 
    ## ... Procrustes: rmse 0.0001755724  max resid 0.0003194085 
    ## ... Similar to previous best
    ## Run 484 stress 9.78257e-05 
    ## ... Procrustes: rmse 0.0001880808  max resid 0.0003294329 
    ## ... Similar to previous best
    ## Run 485 stress 9.927014e-05 
    ## ... Procrustes: rmse 0.0001842661  max resid 0.0003157706 
    ## ... Similar to previous best
    ## Run 486 stress 9.838546e-05 
    ## ... Procrustes: rmse 0.0002006201  max resid 0.0002851012 
    ## ... Similar to previous best
    ## Run 487 stress 9.515723e-05 
    ## ... Procrustes: rmse 0.0001885468  max resid 0.000275883 
    ## ... Similar to previous best
    ## Run 488 stress 9.479108e-05 
    ## ... Procrustes: rmse 0.0001587301  max resid 0.0002506299 
    ## ... Similar to previous best
    ## Run 489 stress 9.524574e-05 
    ## ... Procrustes: rmse 0.0001782474  max resid 0.0003325187 
    ## ... Similar to previous best
    ## Run 490 stress 9.884834e-05 
    ## ... Procrustes: rmse 0.0001659481  max resid 0.0002638846 
    ## ... Similar to previous best
    ## Run 491 stress 9.470997e-05 
    ## ... Procrustes: rmse 0.0001877038  max resid 0.0002708746 
    ## ... Similar to previous best
    ## Run 492 stress 0.2280257 
    ## Run 493 stress 9.612624e-05 
    ## ... Procrustes: rmse 0.000155494  max resid 0.0002317165 
    ## ... Similar to previous best
    ## Run 494 stress 9.222119e-05 
    ## ... Procrustes: rmse 0.0001266149  max resid 0.0002689349 
    ## ... Similar to previous best
    ## Run 495 stress 9.396998e-05 
    ## ... Procrustes: rmse 0.0001745857  max resid 0.0002762644 
    ## ... Similar to previous best
    ## Run 496 stress 9.214312e-05 
    ## ... Procrustes: rmse 0.000188266  max resid 0.0002772514 
    ## ... Similar to previous best
    ## Run 497 stress 9.713737e-05 
    ## ... Procrustes: rmse 0.0001713199  max resid 0.0003152799 
    ## ... Similar to previous best
    ## Run 498 stress 8.879461e-05 
    ## ... Procrustes: rmse 0.0001768585  max resid 0.0002658432 
    ## ... Similar to previous best
    ## Run 499 stress 9.868518e-05 
    ## ... Procrustes: rmse 0.0001651788  max resid 0.0003255468 
    ## ... Similar to previous best
    ## Run 500 stress 9.502479e-05 
    ## ... Procrustes: rmse 0.0001953796  max resid 0.0002778712 
    ## ... Similar to previous best
    ## *** Best solution repeated 485 times

    ## Warning in metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
SD_beta_env_M_NMDS <- metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 0.1491957 
    ## Run 2 stress 0.1990774 
    ## Run 3 stress 0.1990774 
    ## Run 4 stress 0.0004814814 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04159832  max resid 0.08170322 
    ## Run 5 stress 0.1990774 
    ## Run 6 stress 9.619117e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01601418  max resid 0.02187912 
    ## Run 7 stress 9.996951e-05 
    ## ... Procrustes: rmse 0.0104663  max resid 0.01428039 
    ## Run 8 stress 0.001096824 
    ## Run 9 stress 0.001302764 
    ## Run 10 stress 8.870058e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001359688  max resid 0.0002471419 
    ## ... Similar to previous best
    ## Run 11 stress 0.0004675386 
    ## ... Procrustes: rmse 0.01583304  max resid 0.02167642 
    ## Run 12 stress 0.001253046 
    ## Run 13 stress 9.315546e-05 
    ## ... Procrustes: rmse 0.0001509121  max resid 0.0002277066 
    ## ... Similar to previous best
    ## Run 14 stress 0.0002898364 
    ## ... Procrustes: rmse 0.01245973  max resid 0.01702131 
    ## Run 15 stress 9.678494e-05 
    ## ... Procrustes: rmse 0.0001585914  max resid 0.0002492605 
    ## ... Similar to previous best
    ## Run 16 stress 0.2848515 
    ## Run 17 stress 7.522102e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001746976  max resid 0.0002826454 
    ## ... Similar to previous best
    ## Run 18 stress 0.2842805 
    ## Run 19 stress 9.124073e-05 
    ## ... Procrustes: rmse 0.0001763086  max resid 0.000274083 
    ## ... Similar to previous best
    ## Run 20 stress 9.977711e-05 
    ## ... Procrustes: rmse 0.0001873419  max resid 0.0003110329 
    ## ... Similar to previous best
    ## Run 21 stress 0.001249014 
    ## Run 22 stress 0.00033222 
    ## ... Procrustes: rmse 0.01325596  max resid 0.0182057 
    ## Run 23 stress 9.997297e-05 
    ## ... Procrustes: rmse 0.007220591  max resid 0.009877059 
    ## Run 24 stress 9.082118e-05 
    ## ... Procrustes: rmse 0.0001197595  max resid 0.0001919506 
    ## ... Similar to previous best
    ## Run 25 stress 9.25371e-05 
    ## ... Procrustes: rmse 8.489818e-05  max resid 0.0001415685 
    ## ... Similar to previous best
    ## Run 26 stress 0.001359751 
    ## Run 27 stress 0.1491957 
    ## Run 28 stress 0.0004747743 
    ## ... Procrustes: rmse 0.01587107  max resid 0.02181367 
    ## Run 29 stress 0.0005189915 
    ## ... Procrustes: rmse 0.01657658  max resid 0.02279094 
    ## Run 30 stress 0.0001804673 
    ## ... Procrustes: rmse 0.01623531  max resid 0.02218089 
    ## Run 31 stress 0.001183935 
    ## Run 32 stress 9.938693e-05 
    ## ... Procrustes: rmse 0.0001952586  max resid 0.0004177234 
    ## ... Similar to previous best
    ## Run 33 stress 0.001236667 
    ## Run 34 stress 9.961084e-05 
    ## ... Procrustes: rmse 0.007153207  max resid 0.009789945 
    ## Run 35 stress 0.001398102 
    ## Run 36 stress 8.274575e-05 
    ## ... Procrustes: rmse 0.002272745  max resid 0.003054225 
    ## ... Similar to previous best
    ## Run 37 stress 9.485043e-05 
    ## ... Procrustes: rmse 0.004979663  max resid 0.006806289 
    ## ... Similar to previous best
    ## Run 38 stress 0.001211614 
    ## Run 39 stress 0.0002954817 
    ## ... Procrustes: rmse 0.01249421  max resid 0.01715515 
    ## Run 40 stress 8.95776e-05 
    ## ... Procrustes: rmse 0.0001869028  max resid 0.0004125883 
    ## ... Similar to previous best
    ## Run 41 stress 0.001381997 
    ## Run 42 stress 0.1990774 
    ## Run 43 stress 8.780718e-05 
    ## ... Procrustes: rmse 7.430811e-05  max resid 0.000147288 
    ## ... Similar to previous best
    ## Run 44 stress 0.001199147 
    ## Run 45 stress 0.1491957 
    ## Run 46 stress 9.454193e-05 
    ## ... Procrustes: rmse 0.0002288266  max resid 0.0004500697 
    ## ... Similar to previous best
    ## Run 47 stress 0.3083098 
    ## Run 48 stress 0.001313446 
    ## Run 49 stress 0.0003297609 
    ## ... Procrustes: rmse 0.02194843  max resid 0.0300827 
    ## Run 50 stress 0.2440272 
    ## Run 51 stress 0.001242845 
    ## Run 52 stress 0.3083099 
    ## Run 53 stress 0.1990774 
    ## Run 54 stress 0.0012699 
    ## Run 55 stress 0.1491957 
    ## Run 56 stress 8.665159e-05 
    ## ... Procrustes: rmse 0.0001064801  max resid 0.0001897533 
    ## ... Similar to previous best
    ## Run 57 stress 0.001145956 
    ## Run 58 stress 8.382743e-05 
    ## ... Procrustes: rmse 0.001977656  max resid 0.002660965 
    ## ... Similar to previous best
    ## Run 59 stress 0.0004961418 
    ## ... Procrustes: rmse 0.01611504  max resid 0.02215335 
    ## Run 60 stress 0.0002907056 
    ## ... Procrustes: rmse 0.0123939  max resid 0.01701625 
    ## Run 61 stress 0.1990774 
    ## Run 62 stress 0.001178982 
    ## Run 63 stress 3.543133e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001184958  max resid 0.0001672832 
    ## ... Similar to previous best
    ## Run 64 stress 0.0008374869 
    ## Run 65 stress 0.2520602 
    ## Run 66 stress 0.2520602 
    ## Run 67 stress 0.2842805 
    ## Run 68 stress 9.105847e-05 
    ## ... Procrustes: rmse 0.0001100525  max resid 0.0001708476 
    ## ... Similar to previous best
    ## Run 69 stress 9.833488e-05 
    ## ... Procrustes: rmse 0.0001430695  max resid 0.0002231815 
    ## ... Similar to previous best
    ## Run 70 stress 0.001388916 
    ## Run 71 stress 9.414963e-05 
    ## ... Procrustes: rmse 0.00012427  max resid 0.0002005488 
    ## ... Similar to previous best
    ## Run 72 stress 9.592112e-05 
    ## ... Procrustes: rmse 0.004795425  max resid 0.006576032 
    ## ... Similar to previous best
    ## Run 73 stress 8.215304e-05 
    ## ... Procrustes: rmse 0.0001226341  max resid 0.0002021224 
    ## ... Similar to previous best
    ## Run 74 stress 0.00124782 
    ## Run 75 stress 9.940549e-05 
    ## ... Procrustes: rmse 0.0001376375  max resid 0.0002215323 
    ## ... Similar to previous best
    ## Run 76 stress 0.001017096 
    ## Run 77 stress 0.001247634 
    ## Run 78 stress 0.001202385 
    ## Run 79 stress 0.001194139 
    ## Run 80 stress 9.990282e-05 
    ## ... Procrustes: rmse 0.0001789043  max resid 0.0002872734 
    ## ... Similar to previous best
    ## Run 81 stress 0.001266852 
    ## Run 82 stress 0.00142466 
    ## Run 83 stress 0.1491957 
    ## Run 84 stress 0.1491957 
    ## Run 85 stress 0.1491957 
    ## Run 86 stress 0.0003761894 
    ## ... Procrustes: rmse 0.01420416  max resid 0.01953305 
    ## Run 87 stress 0.0001033376 
    ## ... Procrustes: rmse 0.007430748  max resid 0.01018601 
    ## Run 88 stress 9.504739e-05 
    ## ... Procrustes: rmse 0.0001350488  max resid 0.0002052793 
    ## ... Similar to previous best
    ## Run 89 stress 0.001317465 
    ## Run 90 stress 8.829034e-05 
    ## ... Procrustes: rmse 0.0001190134  max resid 0.0002207319 
    ## ... Similar to previous best
    ## Run 91 stress 9.706984e-05 
    ## ... Procrustes: rmse 0.0001902204  max resid 0.0003395985 
    ## ... Similar to previous best
    ## Run 92 stress 8.750996e-05 
    ## ... Procrustes: rmse 0.0001127385  max resid 0.0001807065 
    ## ... Similar to previous best
    ## Run 93 stress 0.0004220717 
    ## ... Procrustes: rmse 0.01504703  max resid 0.02069595 
    ## Run 94 stress 9.095737e-05 
    ## ... Procrustes: rmse 0.0001170346  max resid 0.0001872158 
    ## ... Similar to previous best
    ## Run 95 stress 0.1776022 
    ## Run 96 stress 9.32949e-05 
    ## ... Procrustes: rmse 0.0001344496  max resid 0.000180701 
    ## ... Similar to previous best
    ## Run 97 stress 0.001255227 
    ## Run 98 stress 9.854614e-05 
    ## ... Procrustes: rmse 0.0001948726  max resid 0.0003264514 
    ## ... Similar to previous best
    ## Run 99 stress 0.305249 
    ## Run 100 stress 9.753998e-05 
    ## ... Procrustes: rmse 0.0001392163  max resid 0.0001848073 
    ## ... Similar to previous best
    ## Run 101 stress 0.001059593 
    ## Run 102 stress 9.431499e-05 
    ## ... Procrustes: rmse 0.0001913311  max resid 0.0003308601 
    ## ... Similar to previous best
    ## Run 103 stress 0.0004662097 
    ## ... Procrustes: rmse 0.01579251  max resid 0.02172421 
    ## Run 104 stress 9.794789e-05 
    ## ... Procrustes: rmse 0.0001964721  max resid 0.0003339453 
    ## ... Similar to previous best
    ## Run 105 stress 0.0009586571 
    ## Run 106 stress 0.1776011 
    ## Run 107 stress 0.001231452 
    ## Run 108 stress 0.001233901 
    ## Run 109 stress 9.03918e-05 
    ## ... Procrustes: rmse 0.0001096914  max resid 0.0001884908 
    ## ... Similar to previous best
    ## Run 110 stress 0.00124016 
    ## Run 111 stress 0.001099367 
    ## Run 112 stress 7.952888e-05 
    ## ... Procrustes: rmse 6.405608e-05  max resid 0.0001097852 
    ## ... Similar to previous best
    ## Run 113 stress 0.001356294 
    ## Run 114 stress 0.0003657198 
    ## ... Procrustes: rmse 0.01400334  max resid 0.01925582 
    ## Run 115 stress 0.0004960162 
    ## ... Procrustes: rmse 0.01631404  max resid 0.02244371 
    ## Run 116 stress 0.001253019 
    ## Run 117 stress 0.2440272 
    ## Run 118 stress 9.418898e-05 
    ## ... Procrustes: rmse 0.0001677426  max resid 0.0002677386 
    ## ... Similar to previous best
    ## Run 119 stress 0.1491957 
    ## Run 120 stress 0.001290726 
    ## Run 121 stress 0.001223077 
    ## Run 122 stress 9.064193e-05 
    ## ... Procrustes: rmse 0.0001762749  max resid 0.0003161963 
    ## ... Similar to previous best
    ## Run 123 stress 0.1491957 
    ## Run 124 stress 0.1491957 
    ## Run 125 stress 9.592433e-05 
    ## ... Procrustes: rmse 0.0001661367  max resid 0.000294757 
    ## ... Similar to previous best
    ## Run 126 stress 0.001288818 
    ## Run 127 stress 9.777753e-05 
    ## ... Procrustes: rmse 0.0001454416  max resid 0.0002184342 
    ## ... Similar to previous best
    ## Run 128 stress 0.1491957 
    ## Run 129 stress 9.703286e-05 
    ## ... Procrustes: rmse 0.0001552307  max resid 0.0002598264 
    ## ... Similar to previous best
    ## Run 130 stress 0.000418264 
    ## ... Procrustes: rmse 0.0149781  max resid 0.0206009 
    ## Run 131 stress 9.877316e-05 
    ## ... Procrustes: rmse 0.0001754429  max resid 0.0002387914 
    ## ... Similar to previous best
    ## Run 132 stress 0.001275432 
    ## Run 133 stress 9.149709e-05 
    ## ... Procrustes: rmse 0.00821414  max resid 0.01129354 
    ## Run 134 stress 8.642476e-05 
    ## ... Procrustes: rmse 0.0001747974  max resid 0.0003094182 
    ## ... Similar to previous best
    ## Run 135 stress 0.3083098 
    ## Run 136 stress 0.1491957 
    ## Run 137 stress 0.00122085 
    ## Run 138 stress 0.1491957 
    ## Run 139 stress 8.860942e-05 
    ## ... Procrustes: rmse 0.0001556486  max resid 0.0002529289 
    ## ... Similar to previous best
    ## Run 140 stress 9.317209e-05 
    ## ... Procrustes: rmse 0.000179936  max resid 0.000304743 
    ## ... Similar to previous best
    ## Run 141 stress 0.001388212 
    ## Run 142 stress 0.1990774 
    ## Run 143 stress 0.2848518 
    ## Run 144 stress 0.1990774 
    ## Run 145 stress 0.1491957 
    ## Run 146 stress 8.007727e-05 
    ## ... Procrustes: rmse 0.0001419772  max resid 0.0002539906 
    ## ... Similar to previous best
    ## Run 147 stress 0.001167165 
    ## Run 148 stress 0.1776012 
    ## Run 149 stress 0.001428021 
    ## Run 150 stress 0.001091246 
    ## Run 151 stress 9.106768e-05 
    ## ... Procrustes: rmse 0.000135739  max resid 0.000222121 
    ## ... Similar to previous best
    ## Run 152 stress 0.0001157629 
    ## ... Procrustes: rmse 0.00786663  max resid 0.0107876 
    ## Run 153 stress 9.656792e-05 
    ## ... Procrustes: rmse 0.001898829  max resid 0.002571839 
    ## ... Similar to previous best
    ## Run 154 stress 0.001342269 
    ## Run 155 stress 0.001345736 
    ## Run 156 stress 9.766901e-05 
    ## ... Procrustes: rmse 0.0001979701  max resid 0.0003403754 
    ## ... Similar to previous best
    ## Run 157 stress 8.189438e-05 
    ## ... Procrustes: rmse 0.0001543396  max resid 0.0002803417 
    ## ... Similar to previous best
    ## Run 158 stress 0.3078523 
    ## Run 159 stress 0.0007800325 
    ## Run 160 stress 0.001124559 
    ## Run 161 stress 0.0004926327 
    ## ... Procrustes: rmse 0.01625589  max resid 0.02236469 
    ## Run 162 stress 0.001293458 
    ## Run 163 stress 9.079466e-05 
    ## ... Procrustes: rmse 0.0001682887  max resid 0.0003138049 
    ## ... Similar to previous best
    ## Run 164 stress 0.2696744 
    ## Run 165 stress 0.1491957 
    ## Run 166 stress 0.1990774 
    ## Run 167 stress 0.1491957 
    ## Run 168 stress 8.527128e-05 
    ## ... Procrustes: rmse 0.0001251802  max resid 0.0002165315 
    ## ... Similar to previous best
    ## Run 169 stress 0.0001492052 
    ## ... Procrustes: rmse 0.008932581  max resid 0.01225862 
    ## Run 170 stress 0.00136315 
    ## Run 171 stress 8.552246e-05 
    ## ... Procrustes: rmse 0.0001578869  max resid 0.0002698628 
    ## ... Similar to previous best
    ## Run 172 stress 0.001386966 
    ## Run 173 stress 0.001373049 
    ## Run 174 stress 0.001349202 
    ## Run 175 stress 9.268522e-05 
    ## ... Procrustes: rmse 0.004060058  max resid 0.005557782 
    ## ... Similar to previous best
    ## Run 176 stress 9.321367e-05 
    ## ... Procrustes: rmse 0.0001537918  max resid 0.0002054361 
    ## ... Similar to previous best
    ## Run 177 stress 8.727167e-05 
    ## ... Procrustes: rmse 0.0001742091  max resid 0.0003050129 
    ## ... Similar to previous best
    ## Run 178 stress 7.757878e-05 
    ## ... Procrustes: rmse 0.0004875225  max resid 0.0006754089 
    ## ... Similar to previous best
    ## Run 179 stress 9.378349e-05 
    ## ... Procrustes: rmse 0.0001899473  max resid 0.0003299507 
    ## ... Similar to previous best
    ## Run 180 stress 0.001145149 
    ## Run 181 stress 8.362426e-05 
    ## ... Procrustes: rmse 0.003223905  max resid 0.004399877 
    ## ... Similar to previous best
    ## Run 182 stress 9.913849e-05 
    ## ... Procrustes: rmse 0.0001986466  max resid 0.0003411534 
    ## ... Similar to previous best
    ## Run 183 stress 0.001068051 
    ## Run 184 stress 4.70587e-05 
    ## ... Procrustes: rmse 0.0002972121  max resid 0.0004478884 
    ## ... Similar to previous best
    ## Run 185 stress 9.085085e-05 
    ## ... Procrustes: rmse 0.001461654  max resid 0.001973557 
    ## ... Similar to previous best
    ## Run 186 stress 0.000102807 
    ## ... Procrustes: rmse 0.007409121  max resid 0.01015613 
    ## Run 187 stress 0.2797322 
    ## Run 188 stress 9.482705e-05 
    ## ... Procrustes: rmse 0.0001723375  max resid 0.0002552597 
    ## ... Similar to previous best
    ## Run 189 stress 9.256381e-05 
    ## ... Procrustes: rmse 0.0001056929  max resid 0.0001688379 
    ## ... Similar to previous best
    ## Run 190 stress 0.001303237 
    ## Run 191 stress 0.0009930427 
    ## Run 192 stress 9.973139e-05 
    ## ... Procrustes: rmse 0.0001200429  max resid 0.0001490142 
    ## ... Similar to previous best
    ## Run 193 stress 9.455472e-05 
    ## ... Procrustes: rmse 0.0001300944  max resid 0.0002142602 
    ## ... Similar to previous best
    ## Run 194 stress 9.126876e-05 
    ## ... Procrustes: rmse 0.0001728852  max resid 0.0002853965 
    ## ... Similar to previous best
    ## Run 195 stress 9.111247e-05 
    ## ... Procrustes: rmse 0.0001635545  max resid 0.0002458536 
    ## ... Similar to previous best
    ## Run 196 stress 9.409565e-05 
    ## ... Procrustes: rmse 0.0001596941  max resid 0.0002219901 
    ## ... Similar to previous best
    ## Run 197 stress 0.1491957 
    ## Run 198 stress 0.000211344 
    ## ... Procrustes: rmse 0.01755342  max resid 0.02416469 
    ## Run 199 stress 7.469304e-05 
    ## ... Procrustes: rmse 0.0001118743  max resid 0.0002281296 
    ## ... Similar to previous best
    ## Run 200 stress 8.1866e-05 
    ## ... Procrustes: rmse 0.0002493493  max resid 0.0004479974 
    ## ... Similar to previous best
    ## Run 201 stress 0.1491957 
    ## Run 202 stress 9.23575e-05 
    ## ... Procrustes: rmse 0.000183465  max resid 0.0003164616 
    ## ... Similar to previous best
    ## Run 203 stress 9.877252e-05 
    ## ... Procrustes: rmse 0.0001193421  max resid 0.0001852157 
    ## ... Similar to previous best
    ## Run 204 stress 0.177601 
    ## Run 205 stress 0.000452639 
    ## ... Procrustes: rmse 0.01558313  max resid 0.02143616 
    ## Run 206 stress 9.222989e-05 
    ## ... Procrustes: rmse 0.0001215577  max resid 0.0001952064 
    ## ... Similar to previous best
    ## Run 207 stress 0.0004561568 
    ## ... Procrustes: rmse 0.0156423  max resid 0.02151739 
    ## Run 208 stress 0.0001343812 
    ## ... Procrustes: rmse 0.01399259  max resid 0.01924332 
    ## Run 209 stress 0.001275534 
    ## Run 210 stress 0.1990774 
    ## Run 211 stress 0.1990774 
    ## Run 212 stress 8.841606e-05 
    ## ... Procrustes: rmse 0.0006536318  max resid 0.000888697 
    ## ... Similar to previous best
    ## Run 213 stress 0.1990774 
    ## Run 214 stress 8.944206e-05 
    ## ... Procrustes: rmse 0.0001703857  max resid 0.0002905396 
    ## ... Similar to previous best
    ## Run 215 stress 9.338604e-05 
    ## ... Procrustes: rmse 0.0001212404  max resid 0.0002037645 
    ## ... Similar to previous best
    ## Run 216 stress 0.0004915772 
    ## ... Procrustes: rmse 0.01623735  max resid 0.02233838 
    ## Run 217 stress 0.1491957 
    ## Run 218 stress 0.1491957 
    ## Run 219 stress 9.030291e-05 
    ## ... Procrustes: rmse 0.0001404443  max resid 0.000190342 
    ## ... Similar to previous best
    ## Run 220 stress 0.00128834 
    ## Run 221 stress 9.467403e-05 
    ## ... Procrustes: rmse 0.005995658  max resid 0.008223606 
    ## Run 222 stress 0.1990774 
    ## Run 223 stress 0.2520602 
    ## Run 224 stress 9.999235e-05 
    ## ... Procrustes: rmse 0.007309004  max resid 0.01001802 
    ## Run 225 stress 0.001284447 
    ## Run 226 stress 0.0004708033 
    ## ... Procrustes: rmse 0.01588255  max resid 0.02184885 
    ## Run 227 stress 9.442877e-05 
    ## ... Procrustes: rmse 0.0001691731  max resid 0.0002644056 
    ## ... Similar to previous best
    ## Run 228 stress 0.0005098922 
    ## ... Procrustes: rmse 0.01653952  max resid 0.02275998 
    ## Run 229 stress 0.1491957 
    ## Run 230 stress 0.001368021 
    ## Run 231 stress 0.00116326 
    ## Run 232 stress 9.432017e-05 
    ## ... Procrustes: rmse 0.0001348768  max resid 0.000195406 
    ## ... Similar to previous best
    ## Run 233 stress 0.001300987 
    ## Run 234 stress 8.225985e-05 
    ## ... Procrustes: rmse 0.0001477251  max resid 0.0002561004 
    ## ... Similar to previous best
    ## Run 235 stress 8.713239e-05 
    ## ... Procrustes: rmse 0.0001617462  max resid 0.0002854218 
    ## ... Similar to previous best
    ## Run 236 stress 0.001333616 
    ## Run 237 stress 0.0004964793 
    ## ... Procrustes: rmse 0.01631339  max resid 0.02244439 
    ## Run 238 stress 8.907733e-05 
    ## ... Procrustes: rmse 0.0001291599  max resid 0.0001737565 
    ## ... Similar to previous best
    ## Run 239 stress 0.001227963 
    ## Run 240 stress 0.1491957 
    ## Run 241 stress 7.967986e-05 
    ## ... Procrustes: rmse 0.0001193331  max resid 0.000155663 
    ## ... Similar to previous best
    ## Run 242 stress 9.598175e-05 
    ## ... Procrustes: rmse 0.0001290368  max resid 0.0002174616 
    ## ... Similar to previous best
    ## Run 243 stress 0.312013 
    ## Run 244 stress 0.00133682 
    ## Run 245 stress 0.001255661 
    ## Run 246 stress 0.001312864 
    ## Run 247 stress 0.0004531177 
    ## ... Procrustes: rmse 0.01559068  max resid 0.02144628 
    ## Run 248 stress 9.038764e-05 
    ## ... Procrustes: rmse 0.006517366  max resid 0.009014103 
    ## Run 249 stress 9.066042e-05 
    ## ... Procrustes: rmse 0.0001650268  max resid 0.000259029 
    ## ... Similar to previous best
    ## Run 250 stress 9.772262e-05 
    ## ... Procrustes: rmse 0.0001463586  max resid 0.00020773 
    ## ... Similar to previous best
    ## Run 251 stress 0.001232947 
    ## Run 252 stress 9.229606e-05 
    ## ... Procrustes: rmse 0.0001232514  max resid 0.0002445163 
    ## ... Similar to previous best
    ## Run 253 stress 8.844191e-05 
    ## ... Procrustes: rmse 0.0001332497  max resid 0.0001815975 
    ## ... Similar to previous best
    ## Run 254 stress 0.0005176059 
    ## ... Procrustes: rmse 0.01613224  max resid 0.02219212 
    ## Run 255 stress 0.2373687 
    ## Run 256 stress 0.2570111 
    ## Run 257 stress 0.001094909 
    ## Run 258 stress 8.960976e-05 
    ## ... Procrustes: rmse 0.0001134701  max resid 0.0001828187 
    ## ... Similar to previous best
    ## Run 259 stress 0.1491957 
    ## Run 260 stress 0.001273541 
    ## Run 261 stress 0.001151245 
    ## Run 262 stress 0.001308944 
    ## Run 263 stress 9.464508e-05 
    ## ... Procrustes: rmse 0.000148719  max resid 0.0002362382 
    ## ... Similar to previous best
    ## Run 264 stress 0.001272633 
    ## Run 265 stress 0.001339882 
    ## Run 266 stress 7.646179e-05 
    ## ... Procrustes: rmse 0.0005995817  max resid 0.0007882026 
    ## ... Similar to previous best
    ## Run 267 stress 7.475178e-05 
    ## ... Procrustes: rmse 0.0001328785  max resid 0.0002223199 
    ## ... Similar to previous best
    ## Run 268 stress 0.0004972514 
    ## ... Procrustes: rmse 0.01633429  max resid 0.02247196 
    ## Run 269 stress 0.001191327 
    ## Run 270 stress 8.526922e-05 
    ## ... Procrustes: rmse 0.000123899  max resid 0.0001699683 
    ## ... Similar to previous best
    ## Run 271 stress 9.510153e-05 
    ## ... Procrustes: rmse 0.0001665608  max resid 0.0002411974 
    ## ... Similar to previous best
    ## Run 272 stress 0.2848521 
    ## Run 273 stress 9.765053e-05 
    ## ... Procrustes: rmse 0.0001963824  max resid 0.0003338108 
    ## ... Similar to previous best
    ## Run 274 stress 0.001379993 
    ## Run 275 stress 9.944619e-05 
    ## ... Procrustes: rmse 0.0001473624  max resid 0.000237073 
    ## ... Similar to previous best
    ## Run 276 stress 9.402337e-05 
    ## ... Procrustes: rmse 0.0001263675  max resid 0.0002031867 
    ## ... Similar to previous best
    ## Run 277 stress 9.079213e-05 
    ## ... Procrustes: rmse 0.0001373757  max resid 0.0001880536 
    ## ... Similar to previous best
    ## Run 278 stress 0.0005046468 
    ## ... Procrustes: rmse 0.01645524  max resid 0.02263943 
    ## Run 279 stress 0.0007110958 
    ## Run 280 stress 0.001226851 
    ## Run 281 stress 0.0008971784 
    ## Run 282 stress 0.2494706 
    ## Run 283 stress 0.001384955 
    ## Run 284 stress 0.001226793 
    ## Run 285 stress 0.0009764483 
    ## Run 286 stress 0.001306147 
    ## Run 287 stress 0.001177508 
    ## Run 288 stress 0.00136151 
    ## Run 289 stress 9.854047e-05 
    ## ... Procrustes: rmse 0.0001313797  max resid 0.0002644368 
    ## ... Similar to previous best
    ## Run 290 stress 8.927909e-05 
    ## ... Procrustes: rmse 0.0001702977  max resid 0.0002888386 
    ## ... Similar to previous best
    ## Run 291 stress 0.001192989 
    ## Run 292 stress 0.001308413 
    ## Run 293 stress 0.2797319 
    ## Run 294 stress 8.501763e-05 
    ## ... Procrustes: rmse 0.000160635  max resid 0.0002833347 
    ## ... Similar to previous best
    ## Run 295 stress 0.001399334 
    ## Run 296 stress 9.832899e-05 
    ## ... Procrustes: rmse 0.0001784959  max resid 0.0002634914 
    ## ... Similar to previous best
    ## Run 297 stress 0.2852145 
    ## Run 298 stress 0.001225519 
    ## Run 299 stress 0.001172708 
    ## Run 300 stress 9.167531e-05 
    ## ... Procrustes: rmse 0.0001783329  max resid 0.0003232627 
    ## ... Similar to previous best
    ## Run 301 stress 0.001351664 
    ## Run 302 stress 0.1990774 
    ## Run 303 stress 9.811354e-05 
    ## ... Procrustes: rmse 0.0001411543  max resid 0.0001867661 
    ## ... Similar to previous best
    ## Run 304 stress 0.2528294 
    ## Run 305 stress 0.001396459 
    ## Run 306 stress 0.3083098 
    ## Run 307 stress 9.400024e-05 
    ## ... Procrustes: rmse 0.0001759639  max resid 0.0002990641 
    ## ... Similar to previous best
    ## Run 308 stress 0.001321067 
    ## Run 309 stress 0.001435246 
    ## Run 310 stress 0.0004619159 
    ## ... Procrustes: rmse 0.0157405  max resid 0.0216535 
    ## Run 311 stress 0.001371264 
    ## Run 312 stress 0.2842805 
    ## Run 313 stress 9.612158e-05 
    ## ... Procrustes: rmse 0.0001293538  max resid 0.0002203937 
    ## ... Similar to previous best
    ## Run 314 stress 8.650658e-05 
    ## ... Procrustes: rmse 0.001764936  max resid 0.002372597 
    ## ... Similar to previous best
    ## Run 315 stress 0.001452589 
    ## Run 316 stress 0.0004386308 
    ## ... Procrustes: rmse 0.0153394  max resid 0.02109923 
    ## Run 317 stress 9.723694e-05 
    ## ... Procrustes: rmse 8.901334e-05  max resid 0.0001364789 
    ## ... Similar to previous best
    ## Run 318 stress 9.218339e-05 
    ## ... Procrustes: rmse 0.0001549765  max resid 0.0002105717 
    ## ... Similar to previous best
    ## Run 319 stress 0.0004567 
    ## ... Procrustes: rmse 0.02580505  max resid 0.035583 
    ## Run 320 stress 9.688085e-05 
    ## ... Procrustes: rmse 0.0001311918  max resid 0.000210742 
    ## ... Similar to previous best
    ## Run 321 stress 8.735232e-05 
    ## ... Procrustes: rmse 0.0001676983  max resid 0.0003083135 
    ## ... Similar to previous best
    ## Run 322 stress 0.001239482 
    ## Run 323 stress 9.562743e-05 
    ## ... Procrustes: rmse 0.000469027  max resid 0.0006361466 
    ## ... Similar to previous best
    ## Run 324 stress 0.001433072 
    ## Run 325 stress 0.001194113 
    ## Run 326 stress 9.583175e-05 
    ## ... Procrustes: rmse 0.0001729658  max resid 0.0002683552 
    ## ... Similar to previous best
    ## Run 327 stress 0.001179588 
    ## Run 328 stress 9.586375e-05 
    ## ... Procrustes: rmse 0.0001486391  max resid 0.0002265506 
    ## ... Similar to previous best
    ## Run 329 stress 0.001153942 
    ## Run 330 stress 0.0005956866 
    ## Run 331 stress 0.1776018 
    ## Run 332 stress 0.1491957 
    ## Run 333 stress 0.2440272 
    ## Run 334 stress 9.044717e-05 
    ## ... Procrustes: rmse 0.0001585147  max resid 0.0002367323 
    ## ... Similar to previous best
    ## Run 335 stress 0.1776017 
    ## Run 336 stress 0.177601 
    ## Run 337 stress 0.1776018 
    ## Run 338 stress 0.001278853 
    ## Run 339 stress 0.0008073694 
    ## Run 340 stress 0.1491957 
    ## Run 341 stress 0.0001724254 
    ## ... Procrustes: rmse 0.009606359  max resid 0.0131884 
    ## Run 342 stress 0.1491957 
    ## Run 343 stress 0.001360798 
    ## Run 344 stress 7.8047e-05 
    ## ... Procrustes: rmse 0.0005316992  max resid 0.0007856923 
    ## ... Similar to previous best
    ## Run 345 stress 7.996118e-05 
    ## ... Procrustes: rmse 0.0001478972  max resid 0.0002577649 
    ## ... Similar to previous best
    ## Run 346 stress 0.0003882488 
    ## ... Procrustes: rmse 0.0144287  max resid 0.01984269 
    ## Run 347 stress 9.158583e-05 
    ## ... Procrustes: rmse 0.0001685674  max resid 0.0002805547 
    ## ... Similar to previous best
    ## Run 348 stress 0.1990774 
    ## Run 349 stress 9.050115e-05 
    ## ... Procrustes: rmse 0.0001823248  max resid 0.0003130145 
    ## ... Similar to previous best
    ## Run 350 stress 0.0007506382 
    ## Run 351 stress 0.2520602 
    ## Run 352 stress 9.52919e-05 
    ## ... Procrustes: rmse 0.0001039969  max resid 0.0001429054 
    ## ... Similar to previous best
    ## Run 353 stress 9.894884e-05 
    ## ... Procrustes: rmse 0.0001447739  max resid 0.0002146558 
    ## ... Similar to previous best
    ## Run 354 stress 0.0004115092 
    ## ... Procrustes: rmse 0.01465185  max resid 0.0201504 
    ## Run 355 stress 0.1491957 
    ## Run 356 stress 6.523343e-05 
    ## ... Procrustes: rmse 0.0009545666  max resid 0.001317122 
    ## ... Similar to previous best
    ## Run 357 stress 0.001265547 
    ## Run 358 stress 8.537889e-05 
    ## ... Procrustes: rmse 0.0001301714  max resid 0.0002428446 
    ## ... Similar to previous best
    ## Run 359 stress 0.1491957 
    ## Run 360 stress 9.354784e-05 
    ## ... Procrustes: rmse 0.0001873652  max resid 0.0003252549 
    ## ... Similar to previous best
    ## Run 361 stress 0.001172119 
    ## Run 362 stress 0.001163703 
    ## Run 363 stress 0.1491957 
    ## Run 364 stress 0.0001156163 
    ## ... Procrustes: rmse 0.00785548  max resid 0.01077222 
    ## Run 365 stress 0.001330303 
    ## Run 366 stress 9.91305e-05 
    ## ... Procrustes: rmse 0.0001792  max resid 0.0003117066 
    ## ... Similar to previous best
    ## Run 367 stress 9.896008e-05 
    ## ... Procrustes: rmse 0.0001351055  max resid 0.0002171786 
    ## ... Similar to previous best
    ## Run 368 stress 0.00125258 
    ## Run 369 stress 0.1990774 
    ## Run 370 stress 0.001424259 
    ## Run 371 stress 0.1990774 
    ## Run 372 stress 0.001389624 
    ## Run 373 stress 0.0004013553 
    ## ... Procrustes: rmse 0.01467183  max resid 0.02017865 
    ## Run 374 stress 9.134928e-05 
    ## ... Procrustes: rmse 0.0001288958  max resid 0.0001914671 
    ## ... Similar to previous best
    ## Run 375 stress 0.3008322 
    ## Run 376 stress 0.001257818 
    ## Run 377 stress 0.001391148 
    ## Run 378 stress 9.666593e-05 
    ## ... Procrustes: rmse 0.000192113  max resid 0.0003324194 
    ## ... Similar to previous best
    ## Run 379 stress 8.190525e-05 
    ## ... Procrustes: rmse 0.0008622471  max resid 0.001142595 
    ## ... Similar to previous best
    ## Run 380 stress 0.2578647 
    ## Run 381 stress 0.001272149 
    ## Run 382 stress 0.1990774 
    ## Run 383 stress 0.1776013 
    ## Run 384 stress 0.0004644557 
    ## ... Procrustes: rmse 0.0157857  max resid 0.02171602 
    ## Run 385 stress 9.474943e-05 
    ## ... Procrustes: rmse 0.0001117629  max resid 0.0001516999 
    ## ... Similar to previous best
    ## Run 386 stress 0.0012412 
    ## Run 387 stress 9.994463e-05 
    ## ... Procrustes: rmse 0.0001201167  max resid 0.0001590702 
    ## ... Similar to previous best
    ## Run 388 stress 0.1491957 
    ## Run 389 stress 0.001168578 
    ## Run 390 stress 0.001195602 
    ## Run 391 stress 9.106861e-05 
    ## ... Procrustes: rmse 0.0001821822  max resid 0.0003147067 
    ## ... Similar to previous best
    ## Run 392 stress 8.391596e-05 
    ## ... Procrustes: rmse 0.0001169482  max resid 0.0001837705 
    ## ... Similar to previous best
    ## Run 393 stress 0.1776019 
    ## Run 394 stress 9.910133e-05 
    ## ... Procrustes: rmse 0.0001371615  max resid 0.0001802686 
    ## ... Similar to previous best
    ## Run 395 stress 8.767783e-05 
    ## ... Procrustes: rmse 0.0001302135  max resid 0.0001737866 
    ## ... Similar to previous best
    ## Run 396 stress 0.0004723032 
    ## ... Procrustes: rmse 0.01586428  max resid 0.02182065 
    ## Run 397 stress 0.001206938 
    ## Run 398 stress 9.565801e-05 
    ## ... Procrustes: rmse 0.0001295544  max resid 0.000215876 
    ## ... Similar to previous best
    ## Run 399 stress 0.001271457 
    ## Run 400 stress 0.001345717 
    ## Run 401 stress 8.911356e-05 
    ## ... Procrustes: rmse 0.0001556453  max resid 0.0002752515 
    ## ... Similar to previous best
    ## Run 402 stress 0.2494706 
    ## Run 403 stress 0.001291846 
    ## Run 404 stress 8.734365e-05 
    ## ... Procrustes: rmse 0.0001335082  max resid 0.0002170429 
    ## ... Similar to previous best
    ## Run 405 stress 0.001425523 
    ## Run 406 stress 9.071942e-05 
    ## ... Procrustes: rmse 0.0001386924  max resid 0.0002186326 
    ## ... Similar to previous best
    ## Run 407 stress 0.1491957 
    ## Run 408 stress 0.001438186 
    ## Run 409 stress 0.001298415 
    ## Run 410 stress 7.861337e-05 
    ## ... Procrustes: rmse 0.0001009116  max resid 0.0001527884 
    ## ... Similar to previous best
    ## Run 411 stress 0.3083098 
    ## Run 412 stress 0.1491957 
    ## Run 413 stress 8.754358e-05 
    ## ... Procrustes: rmse 0.00016519  max resid 0.0002760873 
    ## ... Similar to previous best
    ## Run 414 stress 0.0004516121 
    ## ... Procrustes: rmse 0.01552428  max resid 0.02135282 
    ## Run 415 stress 9.699512e-05 
    ## ... Procrustes: rmse 0.00013321  max resid 0.0002241223 
    ## ... Similar to previous best
    ## Run 416 stress 8.937681e-05 
    ## ... Procrustes: rmse 0.003133418  max resid 0.004352755 
    ## ... Similar to previous best
    ## Run 417 stress 9.333438e-05 
    ## ... Procrustes: rmse 0.0001643427  max resid 0.0002858147 
    ## ... Similar to previous best
    ## Run 418 stress 9.205602e-05 
    ## ... Procrustes: rmse 0.0001683432  max resid 0.0002944476 
    ## ... Similar to previous best
    ## Run 419 stress 0.001371569 
    ## Run 420 stress 0.001219452 
    ## Run 421 stress 0.001428316 
    ## Run 422 stress 0.0004090896 
    ## ... Procrustes: rmse 0.01480148  max resid 0.02035553 
    ## Run 423 stress 9.86535e-05 
    ## ... Procrustes: rmse 0.0001857339  max resid 0.0002998831 
    ## ... Similar to previous best
    ## Run 424 stress 8.910688e-05 
    ## ... Procrustes: rmse 0.0001059816  max resid 0.0001660725 
    ## ... Similar to previous best
    ## Run 425 stress 0.0004839038 
    ## ... Procrustes: rmse 0.01611311  max resid 0.02216684 
    ## Run 426 stress 0.1776011 
    ## Run 427 stress 0.0004867898 
    ## ... Procrustes: rmse 0.01615875  max resid 0.02222576 
    ## Run 428 stress 9.918182e-05 
    ## ... Procrustes: rmse 0.006737667  max resid 0.00924315 
    ## Run 429 stress 0.001190331 
    ## Run 430 stress 0.2520602 
    ## Run 431 stress 8.2223e-05 
    ## ... Procrustes: rmse 0.002650642  max resid 0.003608044 
    ## ... Similar to previous best
    ## Run 432 stress 0.001323757 
    ## Run 433 stress 0.00126288 
    ## Run 434 stress 0.001211514 
    ## Run 435 stress 0.0007809951 
    ## Run 436 stress 8.516384e-05 
    ## ... Procrustes: rmse 0.0001258793  max resid 0.0001691659 
    ## ... Similar to previous best
    ## Run 437 stress 9.624695e-05 
    ## ... Procrustes: rmse 0.000195281  max resid 0.0003368421 
    ## ... Similar to previous best
    ## Run 438 stress 0.1491957 
    ## Run 439 stress 0.0004622417 
    ## ... Procrustes: rmse 0.01574739  max resid 0.02166293 
    ## Run 440 stress 0.00027032 
    ## ... Procrustes: rmse 0.01202313  max resid 0.01652349 
    ## Run 441 stress 0.1990774 
    ## Run 442 stress 8.158486e-05 
    ## ... Procrustes: rmse 0.0001020002  max resid 0.0001752977 
    ## ... Similar to previous best
    ## Run 443 stress 9.22646e-05 
    ## ... Procrustes: rmse 0.0001217887  max resid 0.0002053412 
    ## ... Similar to previous best
    ## Run 444 stress 0.001370091 
    ## Run 445 stress 0.001195501 
    ## Run 446 stress 0.0004280008 
    ## ... Procrustes: rmse 0.01515232  max resid 0.02084109 
    ## Run 447 stress 9.046782e-05 
    ## ... Procrustes: rmse 0.0001256244  max resid 0.0001625548 
    ## ... Similar to previous best
    ## Run 448 stress 9.061681e-05 
    ## ... Procrustes: rmse 0.0001000025  max resid 0.000158397 
    ## ... Similar to previous best
    ## Run 449 stress 0.001477704 
    ## Run 450 stress 0.1990774 
    ## Run 451 stress 0.0005482779 
    ## Run 452 stress 9.993824e-05 
    ## ... Procrustes: rmse 0.01202602  max resid 0.01653136 
    ## Run 453 stress 0.0009469194 
    ## Run 454 stress 0.0003372065 
    ## ... Procrustes: rmse 0.01344643  max resid 0.01848811 
    ## Run 455 stress 0.0004508478 
    ## ... Procrustes: rmse 0.01555178  max resid 0.02139287 
    ## Run 456 stress 9.057227e-05 
    ## ... Procrustes: rmse 0.0001707294  max resid 0.0002961386 
    ## ... Similar to previous best
    ## Run 457 stress 0.001341793 
    ## Run 458 stress 0.001382918 
    ## Run 459 stress 9.057964e-05 
    ## ... Procrustes: rmse 0.0001834588  max resid 0.0003197875 
    ## ... Similar to previous best
    ## Run 460 stress 0.0002873135 
    ## ... Procrustes: rmse 0.02047157  max resid 0.02820023 
    ## Run 461 stress 7.799312e-05 
    ## ... Procrustes: rmse 0.0003684133  max resid 0.0005595124 
    ## ... Similar to previous best
    ## Run 462 stress 0.001333074 
    ## Run 463 stress 9.873108e-05 
    ## ... Procrustes: rmse 0.0001983319  max resid 0.0003358472 
    ## ... Similar to previous best
    ## Run 464 stress 0.0004389057 
    ## ... Procrustes: rmse 0.01534366  max resid 0.02110857 
    ## Run 465 stress 9.837982e-05 
    ## ... Procrustes: rmse 0.0001794072  max resid 0.0002783074 
    ## ... Similar to previous best
    ## Run 466 stress 0.2797318 
    ## Run 467 stress 9.941912e-05 
    ## ... Procrustes: rmse 0.0001556177  max resid 0.000207864 
    ## ... Similar to previous best
    ## Run 468 stress 0.001321491 
    ## Run 469 stress 9.445276e-05 
    ## ... Procrustes: rmse 0.0001083118  max resid 0.0001737862 
    ## ... Similar to previous best
    ## Run 470 stress 0.001166818 
    ## Run 471 stress 0.2520602 
    ## Run 472 stress 8.201675e-05 
    ## ... Procrustes: rmse 0.0001543516  max resid 0.0002743737 
    ## ... Similar to previous best
    ## Run 473 stress 0.0005132575 
    ## ... Procrustes: rmse 0.01645827  max resid 0.02264367 
    ## Run 474 stress 0.001358139 
    ## Run 475 stress 0.001073918 
    ## Run 476 stress 0.001385965 
    ## Run 477 stress 0.001308309 
    ## Run 478 stress 0.0007031529 
    ## Run 479 stress 8.47256e-05 
    ## ... Procrustes: rmse 9.603496e-05  max resid 0.0001605783 
    ## ... Similar to previous best
    ## Run 480 stress 0.0004964104 
    ## ... Procrustes: rmse 0.01616089  max resid 0.02223165 
    ## Run 481 stress 0.001236342 
    ## Run 482 stress 9.67675e-05 
    ## ... Procrustes: rmse 0.0001116076  max resid 0.0002150827 
    ## ... Similar to previous best
    ## Run 483 stress 0.0004206555 
    ## ... Procrustes: rmse 0.01502041  max resid 0.02065938 
    ## Run 484 stress 0.0004929941 
    ## ... Procrustes: rmse 0.0162596  max resid 0.0223697 
    ## Run 485 stress 9.260375e-05 
    ## ... Procrustes: rmse 0.0001681711  max resid 0.0002537665 
    ## ... Similar to previous best
    ## Run 486 stress 9.692984e-05 
    ## ... Procrustes: rmse 0.0001810106  max resid 0.0003053762 
    ## ... Similar to previous best
    ## Run 487 stress 9.341659e-05 
    ## ... Procrustes: rmse 0.0001354271  max resid 0.0001819683 
    ## ... Similar to previous best
    ## Run 488 stress 0.0003521676 
    ## ... Procrustes: rmse 0.01374193  max resid 0.0188958 
    ## Run 489 stress 0.1776017 
    ## Run 490 stress 0.2373687 
    ## Run 491 stress 8.971416e-05 
    ## ... Procrustes: rmse 0.007600668  max resid 0.01044631 
    ## Run 492 stress 0.0002652559 
    ## ... Procrustes: rmse 0.01192297  max resid 0.01638563 
    ## Run 493 stress 0.1491957 
    ## Run 494 stress 0.1491957 
    ## Run 495 stress 8.725604e-05 
    ## ... Procrustes: rmse 0.0001607855  max resid 0.0002881957 
    ## ... Similar to previous best
    ## Run 496 stress 0.001327872 
    ## Run 497 stress 0.0001594025 
    ## ... Procrustes: rmse 0.009234454  max resid 0.0126752 
    ## Run 498 stress 0.1491957 
    ## Run 499 stress 9.91009e-05 
    ## ... Procrustes: rmse 0.0001778808  max resid 0.0002742074 
    ## ... Similar to previous best
    ## Run 500 stress 9.861527e-05 
    ## ... Procrustes: rmse 0.0001341455  max resid 0.0002168281 
    ## ... Similar to previous best
    ## *** Best solution repeated 155 times

    ## Warning in metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
SD_beta_geo_NMDS <- metaMDS(SD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09069997 
    ## Run 1 stress 0.1121185 
    ## Run 2 stress 0.09098484 
    ## ... Procrustes: rmse 0.009820074  max resid 0.03255582 
    ## Run 3 stress 0.0900685 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04610056  max resid 0.1186678 
    ## Run 4 stress 0.09006845 
    ## ... New best solution
    ## ... Procrustes: rmse 0.00020414  max resid 0.0004807274 
    ## ... Similar to previous best
    ## Run 5 stress 0.09552773 
    ## Run 6 stress 0.09172844 
    ## Run 7 stress 0.0902259 
    ## ... Procrustes: rmse 0.01345913  max resid 0.04270567 
    ## Run 8 stress 0.08994994 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03614325  max resid 0.1175437 
    ## Run 9 stress 0.1078371 
    ## Run 10 stress 0.1078463 
    ## Run 11 stress 0.08994993 
    ## ... New best solution
    ## ... Procrustes: rmse 9.633869e-06  max resid 2.035336e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.09069431 
    ## Run 13 stress 0.09552826 
    ## Run 14 stress 0.1082381 
    ## Run 15 stress 0.09098471 
    ## Run 16 stress 0.09022592 
    ## ... Procrustes: rmse 0.03753703  max resid 0.1187777 
    ## Run 17 stress 0.09172812 
    ## Run 18 stress 0.1078903 
    ## Run 19 stress 0.1079954 
    ## Run 20 stress 0.09022588 
    ## ... Procrustes: rmse 0.03752402  max resid 0.118736 
    ## Run 21 stress 0.09006844 
    ## ... Procrustes: rmse 0.0361524  max resid 0.1167604 
    ## Run 22 stress 0.1082381 
    ## Run 23 stress 0.1110028 
    ## Run 24 stress 0.09069999 
    ## Run 25 stress 0.1078901 
    ## Run 26 stress 0.09022589 
    ## ... Procrustes: rmse 0.03752678  max resid 0.1187439 
    ## Run 27 stress 0.1079961 
    ## Run 28 stress 0.09004936 
    ## ... Procrustes: rmse 0.01075461  max resid 0.03867101 
    ## Run 29 stress 0.1078463 
    ## Run 30 stress 0.09001003 
    ## ... Procrustes: rmse 0.007604779  max resid 0.02498384 
    ## Run 31 stress 0.09001008 
    ## ... Procrustes: rmse 0.007613215  max resid 0.0251729 
    ## Run 32 stress 0.09022599 
    ## ... Procrustes: rmse 0.03750131  max resid 0.1186635 
    ## Run 33 stress 0.1100425 
    ## Run 34 stress 0.09103364 
    ## Run 35 stress 0.09006858 
    ## ... Procrustes: rmse 0.0361846  max resid 0.1168032 
    ## Run 36 stress 0.1074569 
    ## Run 37 stress 0.09004943 
    ## ... Procrustes: rmse 0.01074919  max resid 0.03863247 
    ## Run 38 stress 0.09172815 
    ## Run 39 stress 0.09098469 
    ## Run 40 stress 0.09098472 
    ## Run 41 stress 0.1084408 
    ## Run 42 stress 0.1121179 
    ## Run 43 stress 0.09004947 
    ## ... Procrustes: rmse 0.01074704  max resid 0.03860633 
    ## Run 44 stress 0.09007306 
    ## ... Procrustes: rmse 0.01357981  max resid 0.03843225 
    ## Run 45 stress 0.1082181 
    ## Run 46 stress 0.09001018 
    ## ... Procrustes: rmse 0.007619685  max resid 0.0252708 
    ## Run 47 stress 0.09069435 
    ## Run 48 stress 0.1120472 
    ## Run 49 stress 0.1078877 
    ## Run 50 stress 0.09022592 
    ## ... Procrustes: rmse 0.03753524  max resid 0.118772 
    ## Run 51 stress 0.09006849 
    ## ... Procrustes: rmse 0.03613764  max resid 0.1167378 
    ## Run 52 stress 0.09190044 
    ## Run 53 stress 0.1082182 
    ## Run 54 stress 0.09069997 
    ## Run 55 stress 0.09001009 
    ## ... Procrustes: rmse 0.007614669  max resid 0.02519165 
    ## Run 56 stress 0.09006847 
    ## ... Procrustes: rmse 0.03615722  max resid 0.1167378 
    ## Run 57 stress 0.09004941 
    ## ... Procrustes: rmse 0.01075355  max resid 0.0386653 
    ## Run 58 stress 0.08995 
    ## ... Procrustes: rmse 0.0001028551  max resid 0.0001723136 
    ## ... Similar to previous best
    ## Run 59 stress 0.09098469 
    ## Run 60 stress 0.0902259 
    ## ... Procrustes: rmse 0.03751267  max resid 0.1187023 
    ## Run 61 stress 0.09001005 
    ## ... Procrustes: rmse 0.007611403  max resid 0.02512561 
    ## Run 62 stress 0.09190034 
    ## Run 63 stress 0.09103358 
    ## Run 64 stress 0.09006845 
    ## ... Procrustes: rmse 0.03617007  max resid 0.1167813 
    ## Run 65 stress 0.1104391 
    ## Run 66 stress 0.09172814 
    ## Run 67 stress 0.1078471 
    ## Run 68 stress 0.09639134 
    ## Run 69 stress 0.1078896 
    ## Run 70 stress 0.09172812 
    ## Run 71 stress 0.0909847 
    ## Run 72 stress 0.09001016 
    ## ... Procrustes: rmse 0.007617298  max resid 0.0252431 
    ## Run 73 stress 0.1084409 
    ## Run 74 stress 0.09098469 
    ## Run 75 stress 0.09007292 
    ## ... Procrustes: rmse 0.01358904  max resid 0.03851927 
    ## Run 76 stress 0.09103365 
    ## Run 77 stress 0.0963923 
    ## Run 78 stress 0.09103358 
    ## Run 79 stress 0.09006848 
    ## ... Procrustes: rmse 0.03618276  max resid 0.1168032 
    ## Run 80 stress 0.09006856 
    ## ... Procrustes: rmse 0.03611968  max resid 0.1167028 
    ## Run 81 stress 0.08994993 
    ## ... New best solution
    ## ... Procrustes: rmse 2.880305e-05  max resid 5.016138e-05 
    ## ... Similar to previous best
    ## Run 82 stress 0.107596 
    ## Run 83 stress 0.1082182 
    ## Run 84 stress 0.09098477 
    ## Run 85 stress 0.09007296 
    ## ... Procrustes: rmse 0.01358045  max resid 0.03845078 
    ## Run 86 stress 0.1104395 
    ## Run 87 stress 0.09006861 
    ## ... Procrustes: rmse 0.03620967  max resid 0.116845 
    ## Run 88 stress 0.1074586 
    ## Run 89 stress 0.09552801 
    ## Run 90 stress 0.09006847 
    ## ... Procrustes: rmse 0.03617892  max resid 0.1167817 
    ## Run 91 stress 0.0900685 
    ## ... Procrustes: rmse 0.03619111  max resid 0.1168147 
    ## Run 92 stress 0.09172811 
    ## Run 93 stress 0.1124333 
    ## Run 94 stress 0.09001007 
    ## ... Procrustes: rmse 0.00761517  max resid 0.02519254 
    ## Run 95 stress 0.09006852 
    ## ... Procrustes: rmse 0.03619494  max resid 0.1168241 
    ## Run 96 stress 0.09103358 
    ## Run 97 stress 0.1118971 
    ## Run 98 stress 0.1104394 
    ## Run 99 stress 0.1082183 
    ## Run 100 stress 0.1078011 
    ## Run 101 stress 0.08994997 
    ## ... Procrustes: rmse 0.0001294761  max resid 0.0002133714 
    ## ... Similar to previous best
    ## Run 102 stress 0.09006846 
    ## ... Procrustes: rmse 0.03618036  max resid 0.1167725 
    ## Run 103 stress 0.1082182 
    ## Run 104 stress 0.09552806 
    ## Run 105 stress 0.09022591 
    ## ... Procrustes: rmse 0.03752371  max resid 0.1187351 
    ## Run 106 stress 0.09022593 
    ## ... Procrustes: rmse 0.03751154  max resid 0.1186966 
    ## Run 107 stress 0.09552769 
    ## Run 108 stress 0.09004942 
    ## ... Procrustes: rmse 0.01074497  max resid 0.038604 
    ## Run 109 stress 0.1082182 
    ## Run 110 stress 0.09172838 
    ## Run 111 stress 0.1074593 
    ## Run 112 stress 0.09190036 
    ## Run 113 stress 0.1079949 
    ## Run 114 stress 0.09103363 
    ## Run 115 stress 0.09098469 
    ## Run 116 stress 0.09069996 
    ## Run 117 stress 0.1078012 
    ## Run 118 stress 0.1080016 
    ## Run 119 stress 0.09552752 
    ## Run 120 stress 0.1079967 
    ## Run 121 stress 0.09006854 
    ## ... Procrustes: rmse 0.03619916  max resid 0.1168313 
    ## Run 122 stress 0.1100422 
    ## Run 123 stress 0.09006849 
    ## ... Procrustes: rmse 0.03618981  max resid 0.1168143 
    ## Run 124 stress 0.09552778 
    ## Run 125 stress 0.09006852 
    ## ... Procrustes: rmse 0.03614075  max resid 0.1167462 
    ## Run 126 stress 0.09552806 
    ## Run 127 stress 0.09022594 
    ## ... Procrustes: rmse 0.03754435  max resid 0.1188001 
    ## Run 128 stress 0.1077265 
    ## Run 129 stress 0.113197 
    ## Run 130 stress 0.09069998 
    ## Run 131 stress 0.1078014 
    ## Run 132 stress 0.1082182 
    ## Run 133 stress 0.08995013 
    ## ... Procrustes: rmse 0.0002277016  max resid 0.0003836801 
    ## ... Similar to previous best
    ## Run 134 stress 0.09069431 
    ## Run 135 stress 0.09006869 
    ## ... Procrustes: rmse 0.03610919  max resid 0.1167027 
    ## Run 136 stress 0.1074567 
    ## Run 137 stress 0.1078011 
    ## Run 138 stress 0.1100421 
    ## Run 139 stress 0.09001005 
    ## ... Procrustes: rmse 0.007613965  max resid 0.025155 
    ## Run 140 stress 0.1075961 
    ## Run 141 stress 0.09172813 
    ## Run 142 stress 0.09004935 
    ## ... Procrustes: rmse 0.01075611  max resid 0.03868692 
    ## Run 143 stress 0.1081628 
    ## Run 144 stress 0.08995004 
    ## ... Procrustes: rmse 0.0001623984  max resid 0.0002785246 
    ## ... Similar to previous best
    ## Run 145 stress 0.09022601 
    ## ... Procrustes: rmse 0.03755364  max resid 0.1188288 
    ## Run 146 stress 0.09103358 
    ## Run 147 stress 0.09172823 
    ## Run 148 stress 0.09069999 
    ## Run 149 stress 0.08995007 
    ## ... Procrustes: rmse 0.0001915469  max resid 0.0003228694 
    ## ... Similar to previous best
    ## Run 150 stress 0.1081581 
    ## Run 151 stress 0.09172824 
    ## Run 152 stress 0.09004939 
    ## ... Procrustes: rmse 0.01076255  max resid 0.03872789 
    ## Run 153 stress 0.09006852 
    ## ... Procrustes: rmse 0.03612997  max resid 0.1167462 
    ## Run 154 stress 0.09006853 
    ## ... Procrustes: rmse 0.03612833  max resid 0.1167451 
    ## Run 155 stress 0.09103358 
    ## Run 156 stress 0.09004936 
    ## ... Procrustes: rmse 0.01075293  max resid 0.03866586 
    ## Run 157 stress 0.09006861 
    ## ... Procrustes: rmse 0.03611484  max resid 0.1167428 
    ## Run 158 stress 0.09172824 
    ## Run 159 stress 0.1107994 
    ## Run 160 stress 0.1078908 
    ## Run 161 stress 0.09172811 
    ## Run 162 stress 0.0902259 
    ## ... Procrustes: rmse 0.03753548  max resid 0.1187719 
    ## Run 163 stress 0.1078011 
    ## Run 164 stress 0.1078011 
    ## Run 165 stress 0.09172819 
    ## Run 166 stress 0.09103358 
    ## Run 167 stress 0.1082383 
    ## Run 168 stress 0.1078898 
    ## Run 169 stress 0.1107952 
    ## Run 170 stress 0.1078011 
    ## Run 171 stress 0.09022593 
    ## ... Procrustes: rmse 0.03751225  max resid 0.1186988 
    ## Run 172 stress 0.1078013 
    ## Run 173 stress 0.09069996 
    ## Run 174 stress 0.09552814 
    ## Run 175 stress 0.09006873 
    ## ... Procrustes: rmse 0.03610794  max resid 0.1167467 
    ## Run 176 stress 0.09006862 
    ## ... Procrustes: rmse 0.0362104  max resid 0.1168477 
    ## Run 177 stress 0.1100426 
    ## Run 178 stress 0.1078464 
    ## Run 179 stress 0.08995009 
    ## ... Procrustes: rmse 0.0001789389  max resid 0.0003012075 
    ## ... Similar to previous best
    ## Run 180 stress 0.09007294 
    ## ... Procrustes: rmse 0.01358712  max resid 0.03851247 
    ## Run 181 stress 0.1082381 
    ## Run 182 stress 0.1107992 
    ## Run 183 stress 0.1082381 
    ## Run 184 stress 0.1079973 
    ## Run 185 stress 0.09006854 
    ## ... Procrustes: rmse 0.03612576  max resid 0.1167442 
    ## Run 186 stress 0.09001006 
    ## ... Procrustes: rmse 0.007614591  max resid 0.02517832 
    ## Run 187 stress 0.09552763 
    ## Run 188 stress 0.09172821 
    ## Run 189 stress 0.09022595 
    ## ... Procrustes: rmse 0.03752298  max resid 0.1187335 
    ## Run 190 stress 0.09172811 
    ## Run 191 stress 0.09172814 
    ## Run 192 stress 0.09006863 
    ## ... Procrustes: rmse 0.03621191  max resid 0.1168515 
    ## Run 193 stress 0.09022591 
    ## ... Procrustes: rmse 0.03751257  max resid 0.1187031 
    ## Run 194 stress 0.08994999 
    ## ... Procrustes: rmse 9.904405e-05  max resid 0.0001669572 
    ## ... Similar to previous best
    ## Run 195 stress 0.1084408 
    ## Run 196 stress 0.09001017 
    ## ... Procrustes: rmse 0.00761979  max resid 0.02530426 
    ## Run 197 stress 0.1077274 
    ## Run 198 stress 0.09001002 
    ## ... Procrustes: rmse 0.007606434  max resid 0.02507981 
    ## Run 199 stress 0.1079962 
    ## Run 200 stress 0.0902259 
    ## ... Procrustes: rmse 0.03751625  max resid 0.1187117 
    ## Run 201 stress 0.09022588 
    ## ... Procrustes: rmse 0.03752589  max resid 0.1187436 
    ## Run 202 stress 0.1078468 
    ## Run 203 stress 0.09022591 
    ## ... Procrustes: rmse 0.03753862  max resid 0.1187823 
    ## Run 204 stress 0.1078897 
    ## Run 205 stress 0.108238 
    ## Run 206 stress 0.1082382 
    ## Run 207 stress 0.09022591 
    ## ... Procrustes: rmse 0.03750946  max resid 0.1186933 
    ## Run 208 stress 0.1109058 
    ## Run 209 stress 0.09004937 
    ## ... Procrustes: rmse 0.01075071  max resid 0.03864839 
    ## Run 210 stress 0.09006847 
    ## ... Procrustes: rmse 0.03618253  max resid 0.116799 
    ## Run 211 stress 0.09172813 
    ## Run 212 stress 0.1123445 
    ## Run 213 stress 0.1077697 
    ## Run 214 stress 0.08995007 
    ## ... Procrustes: rmse 0.0001857186  max resid 0.0003128607 
    ## ... Similar to previous best
    ## Run 215 stress 0.1120714 
    ## Run 216 stress 0.09022589 
    ## ... Procrustes: rmse 0.03753078  max resid 0.1187556 
    ## Run 217 stress 0.09001004 
    ## ... Procrustes: rmse 0.007612325  max resid 0.02514787 
    ## Run 218 stress 0.1121184 
    ## Run 219 stress 0.09001016 
    ## ... Procrustes: rmse 0.007620132  max resid 0.02529261 
    ## Run 220 stress 0.1078909 
    ## Run 221 stress 0.09022588 
    ## ... Procrustes: rmse 0.03752679  max resid 0.1187457 
    ## Run 222 stress 0.09022589 
    ## ... Procrustes: rmse 0.0375227  max resid 0.1187326 
    ## Run 223 stress 0.09172813 
    ## Run 224 stress 0.1078903 
    ## Run 225 stress 0.09007294 
    ## ... Procrustes: rmse 0.01359154  max resid 0.03853829 
    ## Run 226 stress 0.1078358 
    ## Run 227 stress 0.09006859 
    ## ... Procrustes: rmse 0.03620486  max resid 0.1168413 
    ## Run 228 stress 0.09552804 
    ## Run 229 stress 0.09639126 
    ## Run 230 stress 0.1077692 
    ## Run 231 stress 0.1079975 
    ## Run 232 stress 0.1078903 
    ## Run 233 stress 0.09006867 
    ## ... Procrustes: rmse 0.03621125  max resid 0.1168542 
    ## Run 234 stress 0.09172832 
    ## Run 235 stress 0.09172829 
    ## Run 236 stress 0.09552797 
    ## Run 237 stress 0.08994993 
    ## ... Procrustes: rmse 6.159119e-05  max resid 0.000100844 
    ## ... Similar to previous best
    ## Run 238 stress 0.1078471 
    ## Run 239 stress 0.09552788 
    ## Run 240 stress 0.1116815 
    ## Run 241 stress 0.1084408 
    ## Run 242 stress 0.08994996 
    ## ... Procrustes: rmse 8.480544e-05  max resid 0.0001435337 
    ## ... Similar to previous best
    ## Run 243 stress 0.09639257 
    ## Run 244 stress 0.0900729 
    ## ... Procrustes: rmse 0.01359107  max resid 0.03853746 
    ## Run 245 stress 0.09004935 
    ## ... Procrustes: rmse 0.01075255  max resid 0.03866143 
    ## Run 246 stress 0.0955281 
    ## Run 247 stress 0.107596 
    ## Run 248 stress 0.09006864 
    ## ... Procrustes: rmse 0.03611061  max resid 0.1167407 
    ## Run 249 stress 0.09022601 
    ## ... Procrustes: rmse 0.03753397  max resid 0.1187681 
    ## Run 250 stress 0.09006847 
    ## ... Procrustes: rmse 0.03618019  max resid 0.1167983 
    ## Run 251 stress 0.09190024 
    ## Run 252 stress 0.09552774 
    ## Run 253 stress 0.1077682 
    ## Run 254 stress 0.1078469 
    ## Run 255 stress 0.09006861 
    ## ... Procrustes: rmse 0.0362062  max resid 0.1168448 
    ## Run 256 stress 0.1078901 
    ## Run 257 stress 0.09098476 
    ## Run 258 stress 0.1077266 
    ## Run 259 stress 0.1074572 
    ## Run 260 stress 0.1078466 
    ## Run 261 stress 0.1077274 
    ## Run 262 stress 0.09006847 
    ## ... Procrustes: rmse 0.03617118  max resid 0.1167826 
    ## Run 263 stress 0.09639128 
    ## Run 264 stress 0.09069998 
    ## Run 265 stress 0.08995006 
    ## ... Procrustes: rmse 0.000182844  max resid 0.0003083436 
    ## ... Similar to previous best
    ## Run 266 stress 0.08994995 
    ## ... Procrustes: rmse 7.316845e-05  max resid 0.000124703 
    ## ... Similar to previous best
    ## Run 267 stress 0.09552753 
    ## Run 268 stress 0.09006859 
    ## ... Procrustes: rmse 0.03612037  max resid 0.116746 
    ## Run 269 stress 0.09007296 
    ## ... Procrustes: rmse 0.01358015  max resid 0.03844635 
    ## Run 270 stress 0.08995001 
    ## ... Procrustes: rmse 0.0001318518  max resid 0.0002220328 
    ## ... Similar to previous best
    ## Run 271 stress 0.1078463 
    ## Run 272 stress 0.1074563 
    ## Run 273 stress 0.09006857 
    ## ... Procrustes: rmse 0.03620152  max resid 0.1168362 
    ## Run 274 stress 0.0900685 
    ## ... Procrustes: rmse 0.03614128  max resid 0.1167539 
    ## Run 275 stress 0.1074565 
    ## Run 276 stress 0.1100423 
    ## Run 277 stress 0.1084406 
    ## Run 278 stress 0.0902259 
    ## ... Procrustes: rmse 0.03753456  max resid 0.1187695 
    ## Run 279 stress 0.1074558 
    ## Run 280 stress 0.09006862 
    ## ... Procrustes: rmse 0.03611361  max resid 0.1167405 
    ## Run 281 stress 0.09098491 
    ## Run 282 stress 0.09022588 
    ## ... Procrustes: rmse 0.0375282  max resid 0.1187501 
    ## Run 283 stress 0.09022589 
    ## ... Procrustes: rmse 0.03751669  max resid 0.1187154 
    ## Run 284 stress 0.08994999 
    ## ... Procrustes: rmse 0.0001290828  max resid 0.0002197247 
    ## ... Similar to previous best
    ## Run 285 stress 0.1078344 
    ## Run 286 stress 0.1079964 
    ## Run 287 stress 0.1075955 
    ## Run 288 stress 0.1077275 
    ## Run 289 stress 0.1121189 
    ## Run 290 stress 0.0910336 
    ## Run 291 stress 0.1104395 
    ## Run 292 stress 0.09004935 
    ## ... Procrustes: rmse 0.01075299  max resid 0.03866553 
    ## Run 293 stress 0.0900685 
    ## ... Procrustes: rmse 0.03619228  max resid 0.1168189 
    ## Run 294 stress 0.1133814 
    ## Run 295 stress 0.09007291 
    ## ... Procrustes: rmse 0.01358571  max resid 0.03849497 
    ## Run 296 stress 0.09022601 
    ## ... Procrustes: rmse 0.03755219  max resid 0.1188233 
    ## Run 297 stress 0.09007299 
    ## ... Procrustes: rmse 0.01357859  max resid 0.03843257 
    ## Run 298 stress 0.1079975 
    ## Run 299 stress 0.1078905 
    ## Run 300 stress 0.09006862 
    ## ... Procrustes: rmse 0.03611587  max resid 0.1167001 
    ## Run 301 stress 0.09639141 
    ## Run 302 stress 0.1110025 
    ## Run 303 stress 0.1104388 
    ## Run 304 stress 0.1132789 
    ## Run 305 stress 0.1079962 
    ## Run 306 stress 0.09172811 
    ## Run 307 stress 0.1081625 
    ## Run 308 stress 0.09006846 
    ## ... Procrustes: rmse 0.03616987  max resid 0.1167779 
    ## Run 309 stress 0.09006859 
    ## ... Procrustes: rmse 0.03611881  max resid 0.1167403 
    ## Run 310 stress 0.09069997 
    ## Run 311 stress 0.09552823 
    ## Run 312 stress 0.09172813 
    ## Run 313 stress 0.09172848 
    ## Run 314 stress 0.1079966 
    ## Run 315 stress 0.08995005 
    ## ... Procrustes: rmse 0.0001786148  max resid 0.0003017974 
    ## ... Similar to previous best
    ## Run 316 stress 0.08995 
    ## ... Procrustes: rmse 0.000136458  max resid 0.00023094 
    ## ... Similar to previous best
    ## Run 317 stress 0.09172827 
    ## Run 318 stress 0.08994996 
    ## ... Procrustes: rmse 8.573982e-05  max resid 0.0001455865 
    ## ... Similar to previous best
    ## Run 319 stress 0.08995006 
    ## ... Procrustes: rmse 0.0001836204  max resid 0.0003095984 
    ## ... Similar to previous best
    ## Run 320 stress 0.09001008 
    ## ... Procrustes: rmse 0.007604444  max resid 0.02493328 
    ## Run 321 stress 0.09103358 
    ## Run 322 stress 0.1082382 
    ## Run 323 stress 0.09004938 
    ## ... Procrustes: rmse 0.01075152  max resid 0.03866038 
    ## Run 324 stress 0.110439 
    ## Run 325 stress 0.09098492 
    ## Run 326 stress 0.1081563 
    ## Run 327 stress 0.09007299 
    ## ... Procrustes: rmse 0.01360066  max resid 0.03861024 
    ## Run 328 stress 0.09098486 
    ## Run 329 stress 0.0902259 
    ## ... Procrustes: rmse 0.03751793  max resid 0.1187171 
    ## Run 330 stress 0.09001002 
    ## ... Procrustes: rmse 0.007610144  max resid 0.02506809 
    ## Run 331 stress 0.09006853 
    ## ... Procrustes: rmse 0.03619751  max resid 0.116826 
    ## Run 332 stress 0.1077276 
    ## Run 333 stress 0.09001011 
    ## ... Procrustes: rmse 0.007616942  max resid 0.02523356 
    ## Run 334 stress 0.1082381 
    ## Run 335 stress 0.11251 
    ## Run 336 stress 0.1078012 
    ## Run 337 stress 0.1082384 
    ## Run 338 stress 0.1121181 
    ## Run 339 stress 0.1077275 
    ## Run 340 stress 0.09069433 
    ## Run 341 stress 0.09552774 
    ## Run 342 stress 0.09022588 
    ## ... Procrustes: rmse 0.03752472  max resid 0.1187385 
    ## Run 343 stress 0.09022589 
    ## ... Procrustes: rmse 0.03752607  max resid 0.118741 
    ## Run 344 stress 0.0963914 
    ## Run 345 stress 0.09098481 
    ## Run 346 stress 0.1077265 
    ## Run 347 stress 0.09022593 
    ## ... Procrustes: rmse 0.03751227  max resid 0.1186997 
    ## Run 348 stress 0.09001012 
    ## ... Procrustes: rmse 0.007618248  max resid 0.02525443 
    ## Run 349 stress 0.1116812 
    ## Run 350 stress 0.110042 
    ## Run 351 stress 0.08995005 
    ## ... Procrustes: rmse 0.0001810608  max resid 0.0003055464 
    ## ... Similar to previous best
    ## Run 352 stress 0.09004937 
    ## ... Procrustes: rmse 0.01075039  max resid 0.038648 
    ## Run 353 stress 0.1083748 
    ## Run 354 stress 0.1074568 
    ## Run 355 stress 0.1078907 
    ## Run 356 stress 0.09098469 
    ## Run 357 stress 0.1078876 
    ## Run 358 stress 0.09069996 
    ## Run 359 stress 0.09098469 
    ## Run 360 stress 0.09006855 
    ## ... Procrustes: rmse 0.03619865  max resid 0.1168308 
    ## Run 361 stress 0.1135133 
    ## Run 362 stress 0.09098469 
    ## Run 363 stress 0.09001006 
    ## ... Procrustes: rmse 0.007614811  max resid 0.02518749 
    ## Run 364 stress 0.09022588 
    ## ... Procrustes: rmse 0.03752596  max resid 0.1187421 
    ## Run 365 stress 0.08995003 
    ## ... Procrustes: rmse 0.0001435862  max resid 0.0002420818 
    ## ... Similar to previous best
    ## Run 366 stress 0.09022589 
    ## ... Procrustes: rmse 0.03752863  max resid 0.1187504 
    ## Run 367 stress 0.1079968 
    ## Run 368 stress 0.1078876 
    ## Run 369 stress 0.09022589 
    ## ... Procrustes: rmse 0.03752309  max resid 0.1187336 
    ## Run 370 stress 0.110795 
    ## Run 371 stress 0.1078367 
    ## Run 372 stress 0.1074543 
    ## Run 373 stress 0.1110027 
    ## Run 374 stress 0.1078904 
    ## Run 375 stress 0.0900494 
    ## ... Procrustes: rmse 0.01074572  max resid 0.03861225 
    ## Run 376 stress 0.1075953 
    ## Run 377 stress 0.1123445 
    ## Run 378 stress 0.1078011 
    ## Run 379 stress 0.1077266 
    ## Run 380 stress 0.09103372 
    ## Run 381 stress 0.09007293 
    ## ... Procrustes: rmse 0.01359504  max resid 0.03857097 
    ## Run 382 stress 0.09001017 
    ## ... Procrustes: rmse 0.007621101  max resid 0.02530337 
    ## Run 383 stress 0.09190026 
    ## Run 384 stress 0.1074578 
    ## Run 385 stress 0.09006857 
    ## ... Procrustes: rmse 0.03619875  max resid 0.116827 
    ## Run 386 stress 0.09006858 
    ## ... Procrustes: rmse 0.03620531  max resid 0.1168397 
    ## Run 387 stress 0.1078906 
    ## Run 388 stress 0.09006852 
    ## ... Procrustes: rmse 0.03613166  max resid 0.1167276 
    ## Run 389 stress 0.09172831 
    ## Run 390 stress 0.1100423 
    ## Run 391 stress 0.09006869 
    ## ... Procrustes: rmse 0.03621829  max resid 0.1168598 
    ## Run 392 stress 0.09069434 
    ## Run 393 stress 0.1136617 
    ## Run 394 stress 0.1078906 
    ## Run 395 stress 0.09639117 
    ## Run 396 stress 0.09069431 
    ## Run 397 stress 0.09069996 
    ## Run 398 stress 0.09022601 
    ## ... Procrustes: rmse 0.03750349  max resid 0.1186704 
    ## Run 399 stress 0.1107948 
    ## Run 400 stress 0.1121188 
    ## Run 401 stress 0.1100422 
    ## Run 402 stress 0.09006867 
    ## ... Procrustes: rmse 0.03610696  max resid 0.1167405 
    ## Run 403 stress 0.09639234 
    ## Run 404 stress 0.1078474 
    ## Run 405 stress 0.08994993 
    ## ... Procrustes: rmse 3.566702e-05  max resid 5.783153e-05 
    ## ... Similar to previous best
    ## Run 406 stress 0.1100424 
    ## Run 407 stress 0.09006854 
    ## ... Procrustes: rmse 0.03619748  max resid 0.1168299 
    ## Run 408 stress 0.09022593 
    ## ... Procrustes: rmse 0.03751553  max resid 0.118711 
    ## Run 409 stress 0.09552755 
    ## Run 410 stress 0.09022588 
    ## ... Procrustes: rmse 0.03752856  max resid 0.1187507 
    ## Run 411 stress 0.09639128 
    ## Run 412 stress 0.09069431 
    ## Run 413 stress 0.09098488 
    ## Run 414 stress 0.1078368 
    ## Run 415 stress 0.09007297 
    ## ... Procrustes: rmse 0.01358004  max resid 0.03844585 
    ## Run 416 stress 0.09006844 
    ## ... Procrustes: rmse 0.03616655  max resid 0.1167762 
    ## Run 417 stress 0.0899501 
    ## ... Procrustes: rmse 0.0002158897  max resid 0.0003638298 
    ## ... Similar to previous best
    ## Run 418 stress 0.1124329 
    ## Run 419 stress 0.0900494 
    ## ... Procrustes: rmse 0.01074717  max resid 0.03862458 
    ## Run 420 stress 0.09098469 
    ## Run 421 stress 0.1074567 
    ## Run 422 stress 0.0900686 
    ## ... Procrustes: rmse 0.03611909  max resid 0.1167035 
    ## Run 423 stress 0.09189032 
    ## Run 424 stress 0.09172821 
    ## Run 425 stress 0.1074554 
    ## Run 426 stress 0.09022593 
    ## ... Procrustes: rmse 0.03753  max resid 0.1187555 
    ## Run 427 stress 0.09552744 
    ## Run 428 stress 0.1078908 
    ## Run 429 stress 0.09006852 
    ## ... Procrustes: rmse 0.03614698  max resid 0.1167392 
    ## Run 430 stress 0.090226 
    ## ... Procrustes: rmse 0.03753783  max resid 0.1187789 
    ## Run 431 stress 0.09004935 
    ## ... Procrustes: rmse 0.01075571  max resid 0.03868464 
    ## Run 432 stress 0.09172829 
    ## Run 433 stress 0.1079955 
    ## Run 434 stress 0.09172811 
    ## Run 435 stress 0.1081574 
    ## Run 436 stress 0.09004936 
    ## ... Procrustes: rmse 0.01075851  max resid 0.03870138 
    ## Run 437 stress 0.1081631 
    ## Run 438 stress 0.09069431 
    ## Run 439 stress 0.09172831 
    ## Run 440 stress 0.0902259 
    ## ... Procrustes: rmse 0.03751899  max resid 0.1187203 
    ## Run 441 stress 0.09006855 
    ## ... Procrustes: rmse 0.03620031  max resid 0.1168317 
    ## Run 442 stress 0.1077266 
    ## Run 443 stress 0.1078465 
    ## Run 444 stress 0.1108 
    ## Run 445 stress 0.08995004 
    ## ... Procrustes: rmse 0.0001720371  max resid 0.0002913076 
    ## ... Similar to previous best
    ## Run 446 stress 0.08994993 
    ## ... Procrustes: rmse 5.623487e-05  max resid 9.48439e-05 
    ## ... Similar to previous best
    ## Run 447 stress 0.1082182 
    ## Run 448 stress 0.09098469 
    ## Run 449 stress 0.09007294 
    ## ... Procrustes: rmse 0.01358381  max resid 0.03848311 
    ## Run 450 stress 0.09190035 
    ## Run 451 stress 0.1124331 
    ## Run 452 stress 0.09006866 
    ## ... Procrustes: rmse 0.03621497  max resid 0.116855 
    ## Run 453 stress 0.0917284 
    ## Run 454 stress 0.107837 
    ## Run 455 stress 0.09069996 
    ## Run 456 stress 0.09190037 
    ## Run 457 stress 0.09001012 
    ## ... Procrustes: rmse 0.007616757  max resid 0.02525462 
    ## Run 458 stress 0.1078879 
    ## Run 459 stress 0.09006854 
    ## ... Procrustes: rmse 0.03619625  max resid 0.1168273 
    ## Run 460 stress 0.1078881 
    ## Run 461 stress 0.1074556 
    ## Run 462 stress 0.1074592 
    ## Run 463 stress 0.09022591 
    ## ... Procrustes: rmse 0.03753507  max resid 0.118771 
    ## Run 464 stress 0.09007302 
    ## ... Procrustes: rmse 0.01357748  max resid 0.03842392 
    ## Run 465 stress 0.09069996 
    ## Run 466 stress 0.108163 
    ## Run 467 stress 0.0902259 
    ## ... Procrustes: rmse 0.0375155  max resid 0.1187101 
    ## Run 468 stress 0.08994999 
    ## ... Procrustes: rmse 0.0001146528  max resid 0.0001937642 
    ## ... Similar to previous best
    ## Run 469 stress 0.09006855 
    ## ... Procrustes: rmse 0.03612418  max resid 0.1167462 
    ## Run 470 stress 0.1083747 
    ## Run 471 stress 0.09103361 
    ## Run 472 stress 0.1074569 
    ## Run 473 stress 0.107888 
    ## Run 474 stress 0.08994995 
    ## ... Procrustes: rmse 5.84008e-05  max resid 0.000106288 
    ## ... Similar to previous best
    ## Run 475 stress 0.09022589 
    ## ... Procrustes: rmse 0.03751985  max resid 0.1187213 
    ## Run 476 stress 0.1100419 
    ## Run 477 stress 0.1082381 
    ## Run 478 stress 0.1077672 
    ## Run 479 stress 0.09007294 
    ## ... Procrustes: rmse 0.01358211  max resid 0.03846473 
    ## Run 480 stress 0.0902259 
    ## ... Procrustes: rmse 0.03753423  max resid 0.1187679 
    ## Run 481 stress 0.09552753 
    ## Run 482 stress 0.09022594 
    ## ... Procrustes: rmse 0.03750995  max resid 0.1186922 
    ## Run 483 stress 0.09069431 
    ## Run 484 stress 0.1077685 
    ## Run 485 stress 0.09552814 
    ## Run 486 stress 0.1120702 
    ## Run 487 stress 0.1078464 
    ## Run 488 stress 0.08995001 
    ## ... Procrustes: rmse 0.0001396194  max resid 0.0002349161 
    ## ... Similar to previous best
    ## Run 489 stress 0.1078464 
    ## Run 490 stress 0.108238 
    ## Run 491 stress 0.1078907 
    ## Run 492 stress 0.09069996 
    ## Run 493 stress 0.09552785 
    ## Run 494 stress 0.1077274 
    ## Run 495 stress 0.1083747 
    ## Run 496 stress 0.1082382 
    ## Run 497 stress 0.1100419 
    ## Run 498 stress 0.1079964 
    ## Run 499 stress 0.09172835 
    ## Run 500 stress 0.09007302 
    ## ... Procrustes: rmse 0.01357705  max resid 0.03841726 
    ## *** Best solution repeated 27 times

``` r
# Mixed and stratified lakes
SD_beta_geo_MS_NMDS <- metaMDS(SD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08659673 
    ## Run 1 stress 0.09820592 
    ## Run 2 stress 0.09769937 
    ## Run 3 stress 0.09926935 
    ## Run 4 stress 0.09669175 
    ## Run 5 stress 0.09719297 
    ## Run 6 stress 0.09610472 
    ## Run 7 stress 0.1009084 
    ## Run 8 stress 0.1008329 
    ## Run 9 stress 0.09610453 
    ## Run 10 stress 0.09610473 
    ## Run 11 stress 0.1015318 
    ## Run 12 stress 0.09639168 
    ## Run 13 stress 0.09610465 
    ## Run 14 stress 0.08709981 
    ## Run 15 stress 0.08709981 
    ## Run 16 stress 0.281041 
    ## Run 17 stress 0.08659676 
    ## ... Procrustes: rmse 0.0001186313  max resid 0.000254124 
    ## ... Similar to previous best
    ## Run 18 stress 0.09659298 
    ## Run 19 stress 0.09639142 
    ## Run 20 stress 0.09610478 
    ## Run 21 stress 0.09769985 
    ## Run 22 stress 0.0992693 
    ## Run 23 stress 0.09610459 
    ## Run 24 stress 0.09820574 
    ## Run 25 stress 0.09820566 
    ## Run 26 stress 0.09610479 
    ## Run 27 stress 0.09769989 
    ## Run 28 stress 0.0963914 
    ## Run 29 stress 0.09610464 
    ## Run 30 stress 0.09820567 
    ## Run 31 stress 0.09781857 
    ## Run 32 stress 0.09933322 
    ## Run 33 stress 0.08659673 
    ## ... New best solution
    ## ... Procrustes: rmse 3.236052e-05  max resid 7.019225e-05 
    ## ... Similar to previous best
    ## Run 34 stress 0.09659308 
    ## Run 35 stress 0.09847314 
    ## Run 36 stress 0.09974298 
    ## Run 37 stress 0.09639136 
    ## Run 38 stress 0.09719305 
    ## Run 39 stress 0.09610472 
    ## Run 40 stress 0.08709978 
    ## Run 41 stress 0.08659674 
    ## ... Procrustes: rmse 4.75751e-05  max resid 0.0001038209 
    ## ... Similar to previous best
    ## Run 42 stress 0.09639142 
    ## Run 43 stress 0.09927017 
    ## Run 44 stress 0.1015318 
    ## Run 45 stress 0.09610459 
    ## Run 46 stress 0.09974298 
    ## Run 47 stress 0.08659681 
    ## ... Procrustes: rmse 0.0001779207  max resid 0.0003851491 
    ## ... Similar to previous best
    ## Run 48 stress 0.08709981 
    ## Run 49 stress 0.08659674 
    ## ... Procrustes: rmse 8.80103e-05  max resid 0.0001913449 
    ## ... Similar to previous best
    ## Run 50 stress 0.09659302 
    ## Run 51 stress 0.09610478 
    ## Run 52 stress 0.09639178 
    ## Run 53 stress 0.09610484 
    ## Run 54 stress 0.09847314 
    ## Run 55 stress 0.08659673 
    ## ... Procrustes: rmse 2.731516e-05  max resid 5.52905e-05 
    ## ... Similar to previous best
    ## Run 56 stress 0.1001417 
    ## Run 57 stress 0.0965929 
    ## Run 58 stress 0.097193 
    ## Run 59 stress 0.1084406 
    ## Run 60 stress 0.08709995 
    ## Run 61 stress 0.08709981 
    ## Run 62 stress 0.1008328 
    ## Run 63 stress 0.09926897 
    ## Run 64 stress 0.09610469 
    ## Run 65 stress 0.09974301 
    ## Run 66 stress 0.09610468 
    ## Run 67 stress 0.09610476 
    ## Run 68 stress 0.0865968 
    ## ... Procrustes: rmse 0.0001238296  max resid 0.000274711 
    ## ... Similar to previous best
    ## Run 69 stress 0.09639158 
    ## Run 70 stress 0.1008328 
    ## Run 71 stress 0.09639163 
    ## Run 72 stress 0.09610499 
    ## Run 73 stress 0.08709979 
    ## Run 74 stress 0.08659676 
    ## ... Procrustes: rmse 9.463751e-05  max resid 0.0002067671 
    ## ... Similar to previous best
    ## Run 75 stress 0.1017606 
    ## Run 76 stress 0.09639152 
    ## Run 77 stress 0.09926943 
    ## Run 78 stress 0.09769993 
    ## Run 79 stress 0.08709986 
    ## Run 80 stress 0.09769969 
    ## Run 81 stress 0.09659303 
    ## Run 82 stress 0.1088714 
    ## Run 83 stress 0.09610472 
    ## Run 84 stress 0.1001415 
    ## Run 85 stress 0.09639136 
    ## Run 86 stress 0.09610461 
    ## Run 87 stress 0.09926848 
    ## Run 88 stress 0.09899135 
    ## Run 89 stress 0.09926845 
    ## Run 90 stress 0.08709982 
    ## Run 91 stress 0.09974297 
    ## Run 92 stress 0.09820579 
    ## Run 93 stress 0.09610453 
    ## Run 94 stress 0.08709979 
    ## Run 95 stress 0.08659672 
    ## ... New best solution
    ## ... Procrustes: rmse 4.638434e-06  max resid 7.893485e-06 
    ## ... Similar to previous best
    ## Run 96 stress 0.08709978 
    ## Run 97 stress 0.09820565 
    ## Run 98 stress 0.09678141 
    ## Run 99 stress 0.09847317 
    ## Run 100 stress 0.08659679 
    ## ... Procrustes: rmse 0.0001508602  max resid 0.0003262685 
    ## ... Similar to previous best
    ## Run 101 stress 0.09719296 
    ## Run 102 stress 0.08659673 
    ## ... Procrustes: rmse 3.758931e-05  max resid 8.013369e-05 
    ## ... Similar to previous best
    ## Run 103 stress 0.09678112 
    ## Run 104 stress 0.0963915 
    ## Run 105 stress 0.09719296 
    ## Run 106 stress 0.09719296 
    ## Run 107 stress 0.0870998 
    ## Run 108 stress 0.1015318 
    ## Run 109 stress 0.08710001 
    ## Run 110 stress 0.08709981 
    ## Run 111 stress 0.09639135 
    ## Run 112 stress 0.2664092 
    ## Run 113 stress 0.3282725 
    ## Run 114 stress 0.09659313 
    ## Run 115 stress 0.1084402 
    ## Run 116 stress 0.08709979 
    ## Run 117 stress 0.09719302 
    ## Run 118 stress 0.0870998 
    ## Run 119 stress 0.09610462 
    ## Run 120 stress 0.09926953 
    ## Run 121 stress 0.1009083 
    ## Run 122 stress 0.08709983 
    ## Run 123 stress 0.08709985 
    ## Run 124 stress 0.09719302 
    ## Run 125 stress 0.08709993 
    ## Run 126 stress 0.3268701 
    ## Run 127 stress 0.1009082 
    ## Run 128 stress 0.09820571 
    ## Run 129 stress 0.09639162 
    ## Run 130 stress 0.09639156 
    ## Run 131 stress 0.09639136 
    ## Run 132 stress 0.09610486 
    ## Run 133 stress 0.09610487 
    ## Run 134 stress 0.09847314 
    ## Run 135 stress 0.0965934 
    ## Run 136 stress 0.1055395 
    ## Run 137 stress 0.09610461 
    ## Run 138 stress 0.08659673 
    ## ... Procrustes: rmse 1.829661e-05  max resid 3.592285e-05 
    ## ... Similar to previous best
    ## Run 139 stress 0.09639166 
    ## Run 140 stress 0.08709978 
    ## Run 141 stress 0.1055398 
    ## Run 142 stress 0.09610462 
    ## Run 143 stress 0.2491548 
    ## Run 144 stress 0.08709978 
    ## Run 145 stress 0.09659306 
    ## Run 146 stress 0.08709992 
    ## Run 147 stress 0.09659312 
    ## Run 148 stress 0.1008329 
    ## Run 149 stress 0.1015319 
    ## Run 150 stress 0.08709979 
    ## Run 151 stress 0.09769999 
    ## Run 152 stress 0.08659675 
    ## ... Procrustes: rmse 7.980152e-05  max resid 0.0001712809 
    ## ... Similar to previous best
    ## Run 153 stress 0.0984733 
    ## Run 154 stress 0.08709989 
    ## Run 155 stress 0.1009081 
    ## Run 156 stress 0.09639143 
    ## Run 157 stress 0.09610465 
    ## Run 158 stress 0.09659313 
    ## Run 159 stress 0.08659676 
    ## ... Procrustes: rmse 8.519075e-05  max resid 0.0001843371 
    ## ... Similar to previous best
    ## Run 160 stress 0.09610466 
    ## Run 161 stress 0.08659673 
    ## ... Procrustes: rmse 5.562694e-05  max resid 0.0001188072 
    ## ... Similar to previous best
    ## Run 162 stress 0.08709979 
    ## Run 163 stress 0.1009082 
    ## Run 164 stress 0.09610455 
    ## Run 165 stress 0.09926882 
    ## Run 166 stress 0.09820567 
    ## Run 167 stress 0.0963916 
    ## Run 168 stress 0.08659677 
    ## ... Procrustes: rmse 0.0001304788  max resid 0.000312189 
    ## ... Similar to previous best
    ## Run 169 stress 0.08709981 
    ## Run 170 stress 0.09847313 
    ## Run 171 stress 0.08709982 
    ## Run 172 stress 0.08659672 
    ## ... New best solution
    ## ... Procrustes: rmse 6.192705e-06  max resid 1.303191e-05 
    ## ... Similar to previous best
    ## Run 173 stress 0.09610462 
    ## Run 174 stress 0.08709979 
    ## Run 175 stress 0.08709979 
    ## Run 176 stress 0.09899135 
    ## Run 177 stress 0.08709987 
    ## Run 178 stress 0.0966916 
    ## Run 179 stress 0.09610455 
    ## Run 180 stress 0.0992683 
    ## Run 181 stress 0.09634629 
    ## Run 182 stress 0.09820565 
    ## Run 183 stress 0.08659673 
    ## ... Procrustes: rmse 8.873383e-06  max resid 2.09494e-05 
    ## ... Similar to previous best
    ## Run 184 stress 0.08659677 
    ## ... Procrustes: rmse 0.0001101178  max resid 0.0002365398 
    ## ... Similar to previous best
    ## Run 185 stress 0.08659682 
    ## ... Procrustes: rmse 0.0001723731  max resid 0.0003722601 
    ## ... Similar to previous best
    ## Run 186 stress 0.09769994 
    ## Run 187 stress 0.08709987 
    ## Run 188 stress 0.08709979 
    ## Run 189 stress 0.09639141 
    ## Run 190 stress 0.08709981 
    ## Run 191 stress 0.08659679 
    ## ... Procrustes: rmse 0.0001333276  max resid 0.0002898973 
    ## ... Similar to previous best
    ## Run 192 stress 0.1001415 
    ## Run 193 stress 0.09926937 
    ## Run 194 stress 0.09719296 
    ## Run 195 stress 0.08709991 
    ## Run 196 stress 0.09926919 
    ## Run 197 stress 0.08709979 
    ## Run 198 stress 0.1055396 
    ## Run 199 stress 0.09719304 
    ## Run 200 stress 0.0961046 
    ## Run 201 stress 0.0870998 
    ## Run 202 stress 0.08659675 
    ## ... Procrustes: rmse 8.870742e-05  max resid 0.00019167 
    ## ... Similar to previous best
    ## Run 203 stress 0.08659675 
    ## ... Procrustes: rmse 8.310911e-05  max resid 0.0001796948 
    ## ... Similar to previous best
    ## Run 204 stress 0.1017607 
    ## Run 205 stress 0.08709978 
    ## Run 206 stress 0.08709986 
    ## Run 207 stress 0.1055396 
    ## Run 208 stress 0.08659676 
    ## ... Procrustes: rmse 0.0001080927  max resid 0.0002340721 
    ## ... Similar to previous best
    ## Run 209 stress 0.1001415 
    ## Run 210 stress 0.1009083 
    ## Run 211 stress 0.09926993 
    ## Run 212 stress 0.09899138 
    ## Run 213 stress 0.09719298 
    ## Run 214 stress 0.08709979 
    ## Run 215 stress 0.09974298 
    ## Run 216 stress 0.08709987 
    ## Run 217 stress 0.08659682 
    ## ... Procrustes: rmse 0.0001570994  max resid 0.0003384447 
    ## ... Similar to previous best
    ## Run 218 stress 0.08659673 
    ## ... Procrustes: rmse 3.510894e-05  max resid 7.568723e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.1001417 
    ## Run 220 stress 0.1001416 
    ## Run 221 stress 0.09719296 
    ## Run 222 stress 0.1017606 
    ## Run 223 stress 0.1017606 
    ## Run 224 stress 0.09719299 
    ## Run 225 stress 0.09719297 
    ## Run 226 stress 0.09974298 
    ## Run 227 stress 0.1009083 
    ## Run 228 stress 0.08709987 
    ## Run 229 stress 0.1055802 
    ## Run 230 stress 0.09820568 
    ## Run 231 stress 0.09926975 
    ## Run 232 stress 0.09659289 
    ## Run 233 stress 0.08709984 
    ## Run 234 stress 0.09659315 
    ## Run 235 stress 0.08709978 
    ## Run 236 stress 0.09678117 
    ## Run 237 stress 0.09820589 
    ## Run 238 stress 0.1140719 
    ## Run 239 stress 0.0961045 
    ## Run 240 stress 0.09610487 
    ## Run 241 stress 0.1001416 
    ## Run 242 stress 0.09669149 
    ## Run 243 stress 0.09678117 
    ## Run 244 stress 0.3107442 
    ## Run 245 stress 0.1015317 
    ## Run 246 stress 0.1009083 
    ## Run 247 stress 0.0984732 
    ## Run 248 stress 0.0963914 
    ## Run 249 stress 0.09669177 
    ## Run 250 stress 0.09926887 
    ## Run 251 stress 0.09899138 
    ## Run 252 stress 0.08709988 
    ## Run 253 stress 0.09639151 
    ## Run 254 stress 0.09770026 
    ## Run 255 stress 0.1083837 
    ## Run 256 stress 0.09974304 
    ## Run 257 stress 0.08659675 
    ## ... Procrustes: rmse 8.914153e-05  max resid 0.0001926819 
    ## ... Similar to previous best
    ## Run 258 stress 0.09639152 
    ## Run 259 stress 0.09610463 
    ## Run 260 stress 0.09669174 
    ## Run 261 stress 0.1001415 
    ## Run 262 stress 0.08709984 
    ## Run 263 stress 0.08709982 
    ## Run 264 stress 0.0870999 
    ## Run 265 stress 0.09610464 
    ## Run 266 stress 0.09610503 
    ## Run 267 stress 0.08659672 
    ## ... New best solution
    ## ... Procrustes: rmse 4.988099e-06  max resid 1.046062e-05 
    ## ... Similar to previous best
    ## Run 268 stress 0.09669183 
    ## Run 269 stress 0.09678121 
    ## Run 270 stress 0.09770004 
    ## Run 271 stress 0.09926872 
    ## Run 272 stress 0.09639136 
    ## Run 273 stress 0.08709986 
    ## Run 274 stress 0.08659673 
    ## ... Procrustes: rmse 3.051906e-05  max resid 6.561353e-05 
    ## ... Similar to previous best
    ## Run 275 stress 0.09769948 
    ## Run 276 stress 0.09769979 
    ## Run 277 stress 0.09899121 
    ## Run 278 stress 0.1017605 
    ## Run 279 stress 0.09678133 
    ## Run 280 stress 0.09610468 
    ## Run 281 stress 0.09769965 
    ## Run 282 stress 0.09974306 
    ## Run 283 stress 0.0993334 
    ## Run 284 stress 0.09639138 
    ## Run 285 stress 0.1009079 
    ## Run 286 stress 0.09719297 
    ## Run 287 stress 0.08659673 
    ## ... Procrustes: rmse 3.713556e-05  max resid 5.834173e-05 
    ## ... Similar to previous best
    ## Run 288 stress 0.1084399 
    ## Run 289 stress 0.08659673 
    ## ... Procrustes: rmse 4.17388e-05  max resid 9.250224e-05 
    ## ... Similar to previous best
    ## Run 290 stress 0.09820565 
    ## Run 291 stress 0.1009081 
    ## Run 292 stress 0.1023117 
    ## Run 293 stress 0.1030425 
    ## Run 294 stress 0.08659676 
    ## ... Procrustes: rmse 0.0001047294  max resid 0.0002262833 
    ## ... Similar to previous best
    ## Run 295 stress 0.1001416 
    ## Run 296 stress 0.09669145 
    ## Run 297 stress 0.1084405 
    ## Run 298 stress 0.09610451 
    ## Run 299 stress 0.09719301 
    ## Run 300 stress 0.1084405 
    ## Run 301 stress 0.09639153 
    ## Run 302 stress 0.08659674 
    ## ... Procrustes: rmse 6.853426e-05  max resid 0.000148278 
    ## ... Similar to previous best
    ## Run 303 stress 0.08709978 
    ## Run 304 stress 0.08709981 
    ## Run 305 stress 0.09669152 
    ## Run 306 stress 0.1009083 
    ## Run 307 stress 0.1025583 
    ## Run 308 stress 0.09669164 
    ## Run 309 stress 0.09820567 
    ## Run 310 stress 0.09610484 
    ## Run 311 stress 0.1001416 
    ## Run 312 stress 0.09719305 
    ## Run 313 stress 0.330392 
    ## Run 314 stress 0.09610452 
    ## Run 315 stress 0.09678119 
    ## Run 316 stress 0.09610471 
    ## Run 317 stress 0.09820566 
    ## Run 318 stress 0.08709987 
    ## Run 319 stress 0.09659323 
    ## Run 320 stress 0.09770022 
    ## Run 321 stress 0.1084406 
    ## Run 322 stress 0.09719296 
    ## Run 323 stress 0.09678134 
    ## Run 324 stress 0.0967811 
    ## Run 325 stress 0.09974311 
    ## Run 326 stress 0.09974318 
    ## Run 327 stress 0.09678135 
    ## Run 328 stress 0.08709989 
    ## Run 329 stress 0.09781866 
    ## Run 330 stress 0.09610459 
    ## Run 331 stress 0.09926946 
    ## Run 332 stress 0.09769967 
    ## Run 333 stress 0.09770019 
    ## Run 334 stress 0.09769984 
    ## Run 335 stress 0.09639165 
    ## Run 336 stress 0.09610479 
    ## Run 337 stress 0.08709989 
    ## Run 338 stress 0.08659688 
    ## ... Procrustes: rmse 0.0002198988  max resid 0.0004760237 
    ## ... Similar to previous best
    ## Run 339 stress 0.08709988 
    ## Run 340 stress 0.1017604 
    ## Run 341 stress 0.09974312 
    ## Run 342 stress 0.1001415 
    ## Run 343 stress 0.08709981 
    ## Run 344 stress 0.08659672 
    ## ... New best solution
    ## ... Procrustes: rmse 1.859865e-06  max resid 3.030863e-06 
    ## ... Similar to previous best
    ## Run 345 stress 0.09610468 
    ## Run 346 stress 0.09847334 
    ## Run 347 stress 0.08709981 
    ## Run 348 stress 0.09610455 
    ## Run 349 stress 0.09820565 
    ## Run 350 stress 0.08709978 
    ## Run 351 stress 0.1015317 
    ## Run 352 stress 0.08709987 
    ## Run 353 stress 0.09639157 
    ## Run 354 stress 0.09719319 
    ## Run 355 stress 0.09610481 
    ## Run 356 stress 0.09639158 
    ## Run 357 stress 0.08659687 
    ## ... Procrustes: rmse 0.0002171432  max resid 0.0004711477 
    ## ... Similar to previous best
    ## Run 358 stress 0.09610502 
    ## Run 359 stress 0.1083834 
    ## Run 360 stress 0.09669182 
    ## Run 361 stress 0.1084405 
    ## Run 362 stress 0.0865968 
    ## ... Procrustes: rmse 0.0001555918  max resid 0.0003370609 
    ## ... Similar to previous best
    ## Run 363 stress 0.1001415 
    ## Run 364 stress 0.09610461 
    ## Run 365 stress 0.0961047 
    ## Run 366 stress 0.09610464 
    ## Run 367 stress 0.1001415 
    ## Run 368 stress 0.09926915 
    ## Run 369 stress 0.09639136 
    ## Run 370 stress 0.08659673 
    ## ... Procrustes: rmse 5.458522e-05  max resid 0.0001191171 
    ## ... Similar to previous best
    ## Run 371 stress 0.08659684 
    ## ... Procrustes: rmse 0.0001913488  max resid 0.0004163375 
    ## ... Similar to previous best
    ## Run 372 stress 0.1017603 
    ## Run 373 stress 0.09639166 
    ## Run 374 stress 0.09678116 
    ## Run 375 stress 0.08709986 
    ## Run 376 stress 0.1023119 
    ## Run 377 stress 0.09719299 
    ## Run 378 stress 0.09678116 
    ## Run 379 stress 0.08659679 
    ## ... Procrustes: rmse 0.0001185234  max resid 0.0002557821 
    ## ... Similar to previous best
    ## Run 380 stress 0.0967811 
    ## Run 381 stress 0.0963916 
    ## Run 382 stress 0.08709987 
    ## Run 383 stress 0.0870998 
    ## Run 384 stress 0.1017607 
    ## Run 385 stress 0.08709996 
    ## Run 386 stress 0.1001416 
    ## Run 387 stress 0.08709978 
    ## Run 388 stress 0.08659675 
    ## ... Procrustes: rmse 9.75001e-05  max resid 0.0002029656 
    ## ... Similar to previous best
    ## Run 389 stress 0.09639136 
    ## Run 390 stress 0.3345084 
    ## Run 391 stress 0.1023118 
    ## Run 392 stress 0.09926913 
    ## Run 393 stress 0.09634639 
    ## Run 394 stress 0.0870998 
    ## Run 395 stress 0.09639139 
    ## Run 396 stress 0.0961045 
    ## Run 397 stress 0.0870998 
    ## Run 398 stress 0.09639153 
    ## Run 399 stress 0.09610458 
    ## Run 400 stress 0.09719299 
    ## Run 401 stress 0.09926815 
    ## Run 402 stress 0.08709986 
    ## Run 403 stress 0.2519916 
    ## Run 404 stress 0.08709985 
    ## Run 405 stress 0.08709986 
    ## Run 406 stress 0.09926901 
    ## Run 407 stress 0.1001416 
    ## Run 408 stress 0.08659676 
    ## ... Procrustes: rmse 9.449617e-05  max resid 0.0002033554 
    ## ... Similar to previous best
    ## Run 409 stress 0.09669163 
    ## Run 410 stress 0.09769976 
    ## Run 411 stress 0.1084402 
    ## Run 412 stress 0.08659675 
    ## ... Procrustes: rmse 0.0001032663  max resid 0.0002246931 
    ## ... Similar to previous best
    ## Run 413 stress 0.1055397 
    ## Run 414 stress 0.09847313 
    ## Run 415 stress 0.09610476 
    ## Run 416 stress 0.1001416 
    ## Run 417 stress 0.1055803 
    ## Run 418 stress 0.1009083 
    ## Run 419 stress 0.09669152 
    ## Run 420 stress 0.09820567 
    ## Run 421 stress 0.09639143 
    ## Run 422 stress 0.08659673 
    ## ... Procrustes: rmse 5.60188e-06  max resid 1.092875e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.09634638 
    ## Run 424 stress 0.09847316 
    ## Run 425 stress 0.09719302 
    ## Run 426 stress 0.1083837 
    ## Run 427 stress 0.09610482 
    ## Run 428 stress 0.1001415 
    ## Run 429 stress 0.09899095 
    ## Run 430 stress 0.09610462 
    ## Run 431 stress 0.09926834 
    ## Run 432 stress 0.09769961 
    ## Run 433 stress 0.09639152 
    ## Run 434 stress 0.09820566 
    ## Run 435 stress 0.09974305 
    ## Run 436 stress 0.09669146 
    ## Run 437 stress 0.09669153 
    ## Run 438 stress 0.09669149 
    ## Run 439 stress 0.08709979 
    ## Run 440 stress 0.1015317 
    ## Run 441 stress 0.0870998 
    ## Run 442 stress 0.0997431 
    ## Run 443 stress 0.1001416 
    ## Run 444 stress 0.08659673 
    ## ... Procrustes: rmse 3.967269e-05  max resid 8.733763e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.09719306 
    ## Run 446 stress 0.09639146 
    ## Run 447 stress 0.0963915 
    ## Run 448 stress 0.1088716 
    ## Run 449 stress 0.08659681 
    ## ... Procrustes: rmse 0.0001365452  max resid 0.0002970131 
    ## ... Similar to previous best
    ## Run 450 stress 0.09847326 
    ## Run 451 stress 0.09974314 
    ## Run 452 stress 0.09974298 
    ## Run 453 stress 0.08709983 
    ## Run 454 stress 0.09639151 
    ## Run 455 stress 0.08709996 
    ## Run 456 stress 0.1084404 
    ## Run 457 stress 0.08709981 
    ## Run 458 stress 0.09678121 
    ## Run 459 stress 0.08659675 
    ## ... Procrustes: rmse 7.583209e-05  max resid 0.0001631205 
    ## ... Similar to previous best
    ## Run 460 stress 0.09639145 
    ## Run 461 stress 0.09719298 
    ## Run 462 stress 0.1023117 
    ## Run 463 stress 0.1001415 
    ## Run 464 stress 0.100908 
    ## Run 465 stress 0.08659674 
    ## ... Procrustes: rmse 5.961702e-05  max resid 0.0001283423 
    ## ... Similar to previous best
    ## Run 466 stress 0.3099054 
    ## Run 467 stress 0.09769987 
    ## Run 468 stress 0.2668227 
    ## Run 469 stress 0.1008329 
    ## Run 470 stress 0.09678122 
    ## Run 471 stress 0.09678119 
    ## Run 472 stress 0.09719311 
    ## Run 473 stress 0.09610472 
    ## Run 474 stress 0.100833 
    ## Run 475 stress 0.09678112 
    ## Run 476 stress 0.0870998 
    ## Run 477 stress 0.1008329 
    ## Run 478 stress 0.09610481 
    ## Run 479 stress 0.08659673 
    ## ... Procrustes: rmse 9.365084e-06  max resid 1.679689e-05 
    ## ... Similar to previous best
    ## Run 480 stress 0.09820565 
    ## Run 481 stress 0.09678131 
    ## Run 482 stress 0.09610472 
    ## Run 483 stress 0.1084405 
    ## Run 484 stress 0.09639158 
    ## Run 485 stress 0.09610481 
    ## Run 486 stress 0.0976995 
    ## Run 487 stress 0.08659676 
    ## ... Procrustes: rmse 9.63527e-05  max resid 0.0002104326 
    ## ... Similar to previous best
    ## Run 488 stress 0.09669175 
    ## Run 489 stress 0.09678112 
    ## Run 490 stress 0.09669174 
    ## Run 491 stress 0.0992691 
    ## Run 492 stress 0.09770091 
    ## Run 493 stress 0.0967811 
    ## Run 494 stress 0.1009084 
    ## Run 495 stress 0.09659305 
    ## Run 496 stress 0.08659676 
    ## ... Procrustes: rmse 0.0001157467  max resid 0.0002521622 
    ## ... Similar to previous best
    ## Run 497 stress 0.09899078 
    ## Run 498 stress 0.09610465 
    ## Run 499 stress 0.1030417 
    ## Run 500 stress 0.1055802 
    ## *** Best solution repeated 17 times

``` r
# Ocean sites and mixed lakes
SD_beta_geo_OM_NMDS <- metaMDS(SD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06908149 
    ## Run 1 stress 0.06969272 
    ## Run 2 stress 0.06908185 
    ## ... Procrustes: rmse 0.0002842703  max resid 0.0007587292 
    ## ... Similar to previous best
    ## Run 3 stress 0.06899846 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01027662  max resid 0.02912237 
    ## Run 4 stress 0.06969217 
    ## Run 5 stress 0.06969202 
    ## Run 6 stress 0.06899835 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002755458  max resid 0.0005260092 
    ## ... Similar to previous best
    ## Run 7 stress 0.0696371 
    ## Run 8 stress 0.06963712 
    ## Run 9 stress 0.06899846 
    ## ... Procrustes: rmse 0.0001296013  max resid 0.0003110094 
    ## ... Similar to previous best
    ## Run 10 stress 0.06963701 
    ## Run 11 stress 0.06899837 
    ## ... Procrustes: rmse 0.0001601371  max resid 0.0003044439 
    ## ... Similar to previous best
    ## Run 12 stress 0.0696922 
    ## Run 13 stress 0.06969249 
    ## Run 14 stress 0.06908113 
    ## ... Procrustes: rmse 0.009794495  max resid 0.02773047 
    ## Run 15 stress 0.06969272 
    ## Run 16 stress 0.06899835 
    ## ... New best solution
    ## ... Procrustes: rmse 1.883552e-05  max resid 3.504743e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.06899849 
    ## ... Procrustes: rmse 0.00017037  max resid 0.0003702665 
    ## ... Similar to previous best
    ## Run 18 stress 0.06969198 
    ## Run 19 stress 0.06899835 
    ## ... Procrustes: rmse 0.0001108227  max resid 0.000203721 
    ## ... Similar to previous best
    ## Run 20 stress 0.2613752 
    ## Run 21 stress 0.06908176 
    ## ... Procrustes: rmse 0.01033836  max resid 0.02923806 
    ## Run 22 stress 0.06899846 
    ## ... Procrustes: rmse 0.0001441809  max resid 0.0002821771 
    ## ... Similar to previous best
    ## Run 23 stress 0.06963702 
    ## Run 24 stress 0.06899837 
    ## ... Procrustes: rmse 4.379686e-05  max resid 8.676191e-05 
    ## ... Similar to previous best
    ## Run 25 stress 0.06969193 
    ## Run 26 stress 0.06908122 
    ## ... Procrustes: rmse 0.009879213  max resid 0.02794527 
    ## Run 27 stress 0.06969246 
    ## Run 28 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001014662  max resid 0.0002171534 
    ## ... Similar to previous best
    ## Run 29 stress 0.06969193 
    ## Run 30 stress 0.06899834 
    ## ... New best solution
    ## ... Procrustes: rmse 3.544201e-05  max resid 6.818423e-05 
    ## ... Similar to previous best
    ## Run 31 stress 0.06963702 
    ## Run 32 stress 0.069637 
    ## Run 33 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001514697  max resid 0.0002904689 
    ## ... Similar to previous best
    ## Run 34 stress 0.06908164 
    ## ... Procrustes: rmse 0.01027324  max resid 0.02905555 
    ## Run 35 stress 0.06969218 
    ## Run 36 stress 0.06899834 
    ## ... Procrustes: rmse 4.520983e-05  max resid 8.265528e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.06969234 
    ## Run 38 stress 0.06969225 
    ## Run 39 stress 0.06899849 
    ## ... Procrustes: rmse 0.0001819836  max resid 0.0003568277 
    ## ... Similar to previous best
    ## Run 40 stress 0.069637 
    ## Run 41 stress 0.06969309 
    ## Run 42 stress 0.06969283 
    ## Run 43 stress 0.06969283 
    ## Run 44 stress 0.06899839 
    ## ... Procrustes: rmse 0.0001528972  max resid 0.0003006996 
    ## ... Similar to previous best
    ## Run 45 stress 0.06899835 
    ## ... Procrustes: rmse 5.025107e-05  max resid 0.0001044475 
    ## ... Similar to previous best
    ## Run 46 stress 0.06963699 
    ## Run 47 stress 0.06908157 
    ## ... Procrustes: rmse 0.01022229  max resid 0.02891422 
    ## Run 48 stress 0.06899836 
    ## ... Procrustes: rmse 8.24802e-05  max resid 0.000151374 
    ## ... Similar to previous best
    ## Run 49 stress 0.06969215 
    ## Run 50 stress 0.0696371 
    ## Run 51 stress 0.06908172 
    ## ... Procrustes: rmse 0.01032961  max resid 0.02921157 
    ## Run 52 stress 0.06969224 
    ## Run 53 stress 0.0690815 
    ## ... Procrustes: rmse 0.01017171  max resid 0.02875951 
    ## Run 54 stress 0.06969244 
    ## Run 55 stress 0.06908181 
    ## ... Procrustes: rmse 0.01039019  max resid 0.02938279 
    ## Run 56 stress 0.06969284 
    ## Run 57 stress 0.06899834 
    ## ... Procrustes: rmse 1.319788e-05  max resid 2.921722e-05 
    ## ... Similar to previous best
    ## Run 58 stress 0.06963699 
    ## Run 59 stress 0.06908178 
    ## ... Procrustes: rmse 0.01037424  max resid 0.02933503 
    ## Run 60 stress 0.06969252 
    ## Run 61 stress 0.06908081 
    ## ... Procrustes: rmse 0.009500133  max resid 0.02690802 
    ## Run 62 stress 0.06963725 
    ## Run 63 stress 0.06969255 
    ## Run 64 stress 0.06908049 
    ## ... Procrustes: rmse 0.009223441  max resid 0.02614253 
    ## Run 65 stress 0.06969227 
    ## Run 66 stress 0.06908177 
    ## ... Procrustes: rmse 0.01036729  max resid 0.02931648 
    ## Run 67 stress 0.06969231 
    ## Run 68 stress 0.06908067 
    ## ... Procrustes: rmse 0.009365204  max resid 0.02651265 
    ## Run 69 stress 0.06969284 
    ## Run 70 stress 0.06908145 
    ## ... Procrustes: rmse 0.01011966  max resid 0.02864598 
    ## Run 71 stress 0.06969243 
    ## Run 72 stress 0.06963717 
    ## Run 73 stress 0.06969288 
    ## Run 74 stress 0.06969235 
    ## Run 75 stress 0.06969251 
    ## Run 76 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001245178  max resid 0.0002609064 
    ## ... Similar to previous best
    ## Run 77 stress 0.06899846 
    ## ... Procrustes: rmse 0.00019091  max resid 0.0003749665 
    ## ... Similar to previous best
    ## Run 78 stress 0.06969273 
    ## Run 79 stress 0.06908114 
    ## ... Procrustes: rmse 0.009834682  max resid 0.02783924 
    ## Run 80 stress 0.06969232 
    ## Run 81 stress 0.06969265 
    ## Run 82 stress 0.06969198 
    ## Run 83 stress 0.06969203 
    ## Run 84 stress 0.06899836 
    ## ... Procrustes: rmse 3.836129e-05  max resid 6.446142e-05 
    ## ... Similar to previous best
    ## Run 85 stress 0.06969308 
    ## Run 86 stress 0.06908166 
    ## ... Procrustes: rmse 0.01027921  max resid 0.02906921 
    ## Run 87 stress 0.06969261 
    ## Run 88 stress 0.06899835 
    ## ... Procrustes: rmse 6.057125e-05  max resid 0.0001090924 
    ## ... Similar to previous best
    ## Run 89 stress 0.2690151 
    ## Run 90 stress 0.06969269 
    ## Run 91 stress 0.06969249 
    ## Run 92 stress 0.06899838 
    ## ... Procrustes: rmse 9.670626e-05  max resid 0.0001810299 
    ## ... Similar to previous best
    ## Run 93 stress 0.06969204 
    ## Run 94 stress 0.06969269 
    ## Run 95 stress 0.06908115 
    ## ... Procrustes: rmse 0.009849103  max resid 0.02787773 
    ## Run 96 stress 0.06969213 
    ## Run 97 stress 0.06963705 
    ## Run 98 stress 0.06899839 
    ## ... Procrustes: rmse 0.0001117994  max resid 0.0002049776 
    ## ... Similar to previous best
    ## Run 99 stress 0.06969214 
    ## Run 100 stress 0.06899835 
    ## ... Procrustes: rmse 5.456742e-05  max resid 0.000103424 
    ## ... Similar to previous best
    ## Run 101 stress 0.06963719 
    ## Run 102 stress 0.06969231 
    ## Run 103 stress 0.06899838 
    ## ... Procrustes: rmse 0.0001379617  max resid 0.0002634113 
    ## ... Similar to previous best
    ## Run 104 stress 0.06908133 
    ## ... Procrustes: rmse 0.01002536  max resid 0.02835751 
    ## Run 105 stress 0.06969246 
    ## Run 106 stress 0.06969214 
    ## Run 107 stress 0.06963709 
    ## Run 108 stress 0.06969312 
    ## Run 109 stress 0.06963723 
    ## Run 110 stress 0.06899837 
    ## ... Procrustes: rmse 9.066305e-05  max resid 0.0001656049 
    ## ... Similar to previous best
    ## Run 111 stress 0.06963709 
    ## Run 112 stress 0.06969197 
    ## Run 113 stress 0.06899834 
    ## ... New best solution
    ## ... Procrustes: rmse 2.036832e-05  max resid 4.899589e-05 
    ## ... Similar to previous best
    ## Run 114 stress 0.06908169 
    ## ... Procrustes: rmse 0.01031146  max resid 0.0291552 
    ## Run 115 stress 0.06899846 
    ## ... Procrustes: rmse 0.0001896218  max resid 0.0003710204 
    ## ... Similar to previous best
    ## Run 116 stress 0.06963721 
    ## Run 117 stress 0.0690819 
    ## ... Procrustes: rmse 0.0103785  max resid 0.02934498 
    ## Run 118 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001315364  max resid 0.000261117 
    ## ... Similar to previous best
    ## Run 119 stress 0.06899847 
    ## ... Procrustes: rmse 0.0001987846  max resid 0.0003554772 
    ## ... Similar to previous best
    ## Run 120 stress 0.06969226 
    ## Run 121 stress 0.06969305 
    ## Run 122 stress 0.06969283 
    ## Run 123 stress 0.06963715 
    ## Run 124 stress 0.06899841 
    ## ... Procrustes: rmse 0.0001548425  max resid 0.0002754513 
    ## ... Similar to previous best
    ## Run 125 stress 0.06969213 
    ## Run 126 stress 0.06963705 
    ## Run 127 stress 0.069692 
    ## Run 128 stress 0.06908155 
    ## ... Procrustes: rmse 0.01021352  max resid 0.02888016 
    ## Run 129 stress 0.06969191 
    ## Run 130 stress 0.06969246 
    ## Run 131 stress 0.06969302 
    ## Run 132 stress 0.06969243 
    ## Run 133 stress 0.06908189 
    ## ... Procrustes: rmse 0.01042461  max resid 0.02946445 
    ## Run 134 stress 0.06969192 
    ## Run 135 stress 0.06963709 
    ## Run 136 stress 0.06963723 
    ## Run 137 stress 0.06969246 
    ## Run 138 stress 0.06963702 
    ## Run 139 stress 0.06969224 
    ## Run 140 stress 0.06899855 
    ## ... Procrustes: rmse 0.0002321676  max resid 0.0004879172 
    ## ... Similar to previous best
    ## Run 141 stress 0.3336558 
    ## Run 142 stress 0.06969282 
    ## Run 143 stress 0.069637 
    ## Run 144 stress 0.3432082 
    ## Run 145 stress 0.06899836 
    ## ... Procrustes: rmse 6.103995e-05  max resid 0.0001230643 
    ## ... Similar to previous best
    ## Run 146 stress 0.06969245 
    ## Run 147 stress 0.06969223 
    ## Run 148 stress 0.06908165 
    ## ... Procrustes: rmse 0.01029035  max resid 0.02909563 
    ## Run 149 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001489149  max resid 0.0002946943 
    ## ... Similar to previous best
    ## Run 150 stress 0.0690815 
    ## ... Procrustes: rmse 0.01017175  max resid 0.02876659 
    ## Run 151 stress 0.06969267 
    ## Run 152 stress 0.06963706 
    ## Run 153 stress 0.06969262 
    ## Run 154 stress 0.06969263 
    ## Run 155 stress 0.06899843 
    ## ... Procrustes: rmse 0.0001741985  max resid 0.0003476839 
    ## ... Similar to previous best
    ## Run 156 stress 0.06899839 
    ## ... Procrustes: rmse 0.000132119  max resid 0.000247848 
    ## ... Similar to previous best
    ## Run 157 stress 0.0696929 
    ## Run 158 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001220098  max resid 0.0002741588 
    ## ... Similar to previous best
    ## Run 159 stress 0.06963716 
    ## Run 160 stress 0.06969253 
    ## Run 161 stress 0.06969233 
    ## Run 162 stress 0.06908197 
    ## ... Procrustes: rmse 0.01048534  max resid 0.02963677 
    ## Run 163 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001276499  max resid 0.0002780449 
    ## ... Similar to previous best
    ## Run 164 stress 0.0690818 
    ## ... Procrustes: rmse 0.01038217  max resid 0.02935114 
    ## Run 165 stress 0.06908155 
    ## ... Procrustes: rmse 0.01021456  max resid 0.02888643 
    ## Run 166 stress 0.06969252 
    ## Run 167 stress 0.06969255 
    ## Run 168 stress 0.06908164 
    ## ... Procrustes: rmse 0.0102652  max resid 0.02902355 
    ## Run 169 stress 0.06969274 
    ## Run 170 stress 0.06969207 
    ## Run 171 stress 0.06969196 
    ## Run 172 stress 0.06969258 
    ## Run 173 stress 0.06969271 
    ## Run 174 stress 0.06908178 
    ## ... Procrustes: rmse 0.01037365  max resid 0.02932765 
    ## Run 175 stress 0.06969207 
    ## Run 176 stress 0.06908178 
    ## ... Procrustes: rmse 0.01032017  max resid 0.02918064 
    ## Run 177 stress 0.06908162 
    ## ... Procrustes: rmse 0.0102646  max resid 0.02902368 
    ## Run 178 stress 0.06908157 
    ## ... Procrustes: rmse 0.01023176  max resid 0.02893255 
    ## Run 179 stress 0.06908192 
    ## ... Procrustes: rmse 0.01045242  max resid 0.02954329 
    ## Run 180 stress 0.06899847 
    ## ... Procrustes: rmse 0.0002092293  max resid 0.0004338449 
    ## ... Similar to previous best
    ## Run 181 stress 0.06969199 
    ## Run 182 stress 0.06908196 
    ## ... Procrustes: rmse 0.01041163  max resid 0.02942546 
    ## Run 183 stress 0.06969238 
    ## Run 184 stress 0.06969218 
    ## Run 185 stress 0.06969234 
    ## Run 186 stress 0.06899837 
    ## ... Procrustes: rmse 0.0001005057  max resid 0.0001973589 
    ## ... Similar to previous best
    ## Run 187 stress 0.06899836 
    ## ... Procrustes: rmse 8.073241e-05  max resid 0.0002077812 
    ## ... Similar to previous best
    ## Run 188 stress 0.06963704 
    ## Run 189 stress 0.06899845 
    ## ... Procrustes: rmse 0.0001968012  max resid 0.0003766301 
    ## ... Similar to previous best
    ## Run 190 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001167537  max resid 0.0002527257 
    ## ... Similar to previous best
    ## Run 191 stress 0.06899847 
    ## ... Procrustes: rmse 0.0001838217  max resid 0.0003702829 
    ## ... Similar to previous best
    ## Run 192 stress 0.06963716 
    ## Run 193 stress 0.06969237 
    ## Run 194 stress 0.06908158 
    ## ... Procrustes: rmse 0.01022358  max resid 0.0289103 
    ## Run 195 stress 0.06899841 
    ## ... Procrustes: rmse 9.932128e-05  max resid 0.0001886216 
    ## ... Similar to previous best
    ## Run 196 stress 0.06969221 
    ## Run 197 stress 0.06969226 
    ## Run 198 stress 0.06908116 
    ## ... Procrustes: rmse 0.009861421  max resid 0.02790174 
    ## Run 199 stress 0.06969198 
    ## Run 200 stress 0.06969266 
    ## Run 201 stress 0.06963702 
    ## Run 202 stress 0.06969255 
    ## Run 203 stress 0.06969277 
    ## Run 204 stress 0.06908211 
    ## ... Procrustes: rmse 0.01052295  max resid 0.02973843 
    ## Run 205 stress 0.06969305 
    ## Run 206 stress 0.06969289 
    ## Run 207 stress 0.06908198 
    ## ... Procrustes: rmse 0.01043401  max resid 0.02949895 
    ## Run 208 stress 0.06899845 
    ## ... Procrustes: rmse 0.0001913642  max resid 0.000359175 
    ## ... Similar to previous best
    ## Run 209 stress 0.06969244 
    ## Run 210 stress 0.06963728 
    ## Run 211 stress 0.06969273 
    ## Run 212 stress 0.06899845 
    ## ... Procrustes: rmse 0.0001857155  max resid 0.0003569511 
    ## ... Similar to previous best
    ## Run 213 stress 0.06899843 
    ## ... Procrustes: rmse 0.0001773485  max resid 0.0003527651 
    ## ... Similar to previous best
    ## Run 214 stress 0.06969274 
    ## Run 215 stress 0.06899835 
    ## ... Procrustes: rmse 7.285617e-05  max resid 0.00015071 
    ## ... Similar to previous best
    ## Run 216 stress 0.06899843 
    ## ... Procrustes: rmse 0.0001686126  max resid 0.0003392085 
    ## ... Similar to previous best
    ## Run 217 stress 0.06908175 
    ## ... Procrustes: rmse 0.01035399  max resid 0.02927094 
    ## Run 218 stress 0.06969284 
    ## Run 219 stress 0.06969264 
    ## Run 220 stress 0.0690817 
    ## ... Procrustes: rmse 0.01032403  max resid 0.02917742 
    ## Run 221 stress 0.06969277 
    ## Run 222 stress 0.06899834 
    ## ... Procrustes: rmse 1.48493e-05  max resid 3.750473e-05 
    ## ... Similar to previous best
    ## Run 223 stress 0.06969238 
    ## Run 224 stress 0.06899839 
    ## ... Procrustes: rmse 0.0001237021  max resid 0.0002213593 
    ## ... Similar to previous best
    ## Run 225 stress 0.06963701 
    ## Run 226 stress 0.06969227 
    ## Run 227 stress 0.06969208 
    ## Run 228 stress 0.06899834 
    ## ... Procrustes: rmse 2.274703e-05  max resid 5.433257e-05 
    ## ... Similar to previous best
    ## Run 229 stress 0.06969252 
    ## Run 230 stress 0.06899842 
    ## ... Procrustes: rmse 0.000163736  max resid 0.0003070131 
    ## ... Similar to previous best
    ## Run 231 stress 0.06969231 
    ## Run 232 stress 0.06899835 
    ## ... Procrustes: rmse 6.342819e-05  max resid 0.0001262326 
    ## ... Similar to previous best
    ## Run 233 stress 0.0696371 
    ## Run 234 stress 0.06969213 
    ## Run 235 stress 0.06899839 
    ## ... Procrustes: rmse 0.0001382847  max resid 0.0002698666 
    ## ... Similar to previous best
    ## Run 236 stress 0.06963726 
    ## Run 237 stress 0.06969226 
    ## Run 238 stress 0.06963699 
    ## Run 239 stress 0.06899836 
    ## ... Procrustes: rmse 5.776913e-05  max resid 0.0001225211 
    ## ... Similar to previous best
    ## Run 240 stress 0.06899845 
    ## ... Procrustes: rmse 0.0001964136  max resid 0.0003877982 
    ## ... Similar to previous best
    ## Run 241 stress 0.06899839 
    ## ... Procrustes: rmse 0.0001326407  max resid 0.0002413464 
    ## ... Similar to previous best
    ## Run 242 stress 0.3432084 
    ## Run 243 stress 0.06908198 
    ## ... Procrustes: rmse 0.01048662  max resid 0.029646 
    ## Run 244 stress 0.06969266 
    ## Run 245 stress 0.3532058 
    ## Run 246 stress 0.06908087 
    ## ... Procrustes: rmse 0.009556228  max resid 0.0270586 
    ## Run 247 stress 0.06908129 
    ## ... Procrustes: rmse 0.009985149  max resid 0.02824801 
    ## Run 248 stress 0.06899848 
    ## ... Procrustes: rmse 0.0001900961  max resid 0.00037618 
    ## ... Similar to previous best
    ## Run 249 stress 0.06908186 
    ## ... Procrustes: rmse 0.01042543  max resid 0.02947218 
    ## Run 250 stress 0.06899835 
    ## ... Procrustes: rmse 4.919595e-05  max resid 0.000103914 
    ## ... Similar to previous best
    ## Run 251 stress 0.06969286 
    ## Run 252 stress 0.06969194 
    ## Run 253 stress 0.06908193 
    ## ... Procrustes: rmse 0.01042004  max resid 0.02946037 
    ## Run 254 stress 0.06969218 
    ## Run 255 stress 0.0696931 
    ## Run 256 stress 0.06969196 
    ## Run 257 stress 0.06963725 
    ## Run 258 stress 0.06899836 
    ## ... Procrustes: rmse 7.833593e-05  max resid 0.0001736958 
    ## ... Similar to previous best
    ## Run 259 stress 0.06899839 
    ## ... Procrustes: rmse 0.0001247117  max resid 0.000239095 
    ## ... Similar to previous best
    ## Run 260 stress 0.06899836 
    ## ... Procrustes: rmse 8.61718e-05  max resid 0.0001622671 
    ## ... Similar to previous best
    ## Run 261 stress 0.06963726 
    ## Run 262 stress 0.0696927 
    ## Run 263 stress 0.06969262 
    ## Run 264 stress 0.06969273 
    ## Run 265 stress 0.06969294 
    ## Run 266 stress 0.06908165 
    ## ... Procrustes: rmse 0.01028554  max resid 0.02908225 
    ## Run 267 stress 0.06899843 
    ## ... Procrustes: rmse 0.000171198  max resid 0.0003410037 
    ## ... Similar to previous best
    ## Run 268 stress 0.06899843 
    ## ... Procrustes: rmse 0.0001750309  max resid 0.0003370668 
    ## ... Similar to previous best
    ## Run 269 stress 0.06908044 
    ## ... Procrustes: rmse 0.009190102  max resid 0.02604403 
    ## Run 270 stress 0.0696928 
    ## Run 271 stress 0.06969299 
    ## Run 272 stress 0.06963713 
    ## Run 273 stress 0.0690815 
    ## ... Procrustes: rmse 0.01017446  max resid 0.02877366 
    ## Run 274 stress 0.06969217 
    ## Run 275 stress 0.06963708 
    ## Run 276 stress 0.06899836 
    ## ... Procrustes: rmse 5.865168e-05  max resid 0.000117128 
    ## ... Similar to previous best
    ## Run 277 stress 0.06969267 
    ## Run 278 stress 0.06908159 
    ## ... Procrustes: rmse 0.01023166  max resid 0.02893642 
    ## Run 279 stress 0.06969264 
    ## Run 280 stress 0.06969199 
    ## Run 281 stress 0.06899838 
    ## ... Procrustes: rmse 0.0001163638  max resid 0.0002051589 
    ## ... Similar to previous best
    ## Run 282 stress 0.06899835 
    ## ... Procrustes: rmse 6.97746e-05  max resid 0.0001424823 
    ## ... Similar to previous best
    ## Run 283 stress 0.3532055 
    ## Run 284 stress 0.06899838 
    ## ... Procrustes: rmse 0.0001119951  max resid 0.0002113262 
    ## ... Similar to previous best
    ## Run 285 stress 0.06969204 
    ## Run 286 stress 0.06969195 
    ## Run 287 stress 0.06899843 
    ## ... Procrustes: rmse 0.0001830585  max resid 0.0003614652 
    ## ... Similar to previous best
    ## Run 288 stress 0.06963719 
    ## Run 289 stress 0.06969284 
    ## Run 290 stress 0.06969247 
    ## Run 291 stress 0.06969255 
    ## Run 292 stress 0.06969203 
    ## Run 293 stress 0.06969253 
    ## Run 294 stress 0.06969235 
    ## Run 295 stress 0.0690818 
    ## ... Procrustes: rmse 0.01038609  max resid 0.0293595 
    ## Run 296 stress 0.06899835 
    ## ... Procrustes: rmse 6.385447e-05  max resid 0.0001336839 
    ## ... Similar to previous best
    ## Run 297 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001034988  max resid 0.0002255746 
    ## ... Similar to previous best
    ## Run 298 stress 0.06899834 
    ## ... Procrustes: rmse 3.408082e-05  max resid 7.695648e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.06899835 
    ## ... Procrustes: rmse 5.913479e-05  max resid 0.0001526526 
    ## ... Similar to previous best
    ## Run 300 stress 0.06969198 
    ## Run 301 stress 0.06969244 
    ## Run 302 stress 0.06969283 
    ## Run 303 stress 0.06969285 
    ## Run 304 stress 0.06969242 
    ## Run 305 stress 0.06899835 
    ## ... Procrustes: rmse 2.663027e-05  max resid 5.733587e-05 
    ## ... Similar to previous best
    ## Run 306 stress 0.0690817 
    ## ... Procrustes: rmse 0.0103204  max resid 0.02917815 
    ## Run 307 stress 0.06969291 
    ## Run 308 stress 0.06969237 
    ## Run 309 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001260745  max resid 0.0002599426 
    ## ... Similar to previous best
    ## Run 310 stress 0.06969208 
    ## Run 311 stress 0.06969304 
    ## Run 312 stress 0.06969214 
    ## Run 313 stress 0.06969193 
    ## Run 314 stress 0.0690814 
    ## ... Procrustes: rmse 0.01009124  max resid 0.02854086 
    ## Run 315 stress 0.06963711 
    ## Run 316 stress 0.06899836 
    ## ... Procrustes: rmse 8.765595e-05  max resid 0.0001613145 
    ## ... Similar to previous best
    ## Run 317 stress 0.06969205 
    ## Run 318 stress 0.06969278 
    ## Run 319 stress 0.06969258 
    ## Run 320 stress 0.06969251 
    ## Run 321 stress 0.06963706 
    ## Run 322 stress 0.06963723 
    ## Run 323 stress 0.0696926 
    ## Run 324 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001245861  max resid 0.0002495464 
    ## ... Similar to previous best
    ## Run 325 stress 0.06969239 
    ## Run 326 stress 0.06899843 
    ## ... Procrustes: rmse 0.0001727837  max resid 0.000327968 
    ## ... Similar to previous best
    ## Run 327 stress 0.06908187 
    ## ... Procrustes: rmse 0.0104225  max resid 0.02945982 
    ## Run 328 stress 0.06969223 
    ## Run 329 stress 0.06899841 
    ## ... Procrustes: rmse 9.949071e-05  max resid 0.0002112079 
    ## ... Similar to previous best
    ## Run 330 stress 0.06963702 
    ## Run 331 stress 0.06963715 
    ## Run 332 stress 0.06969268 
    ## Run 333 stress 0.06899844 
    ## ... Procrustes: rmse 0.0001885145  max resid 0.0004115524 
    ## ... Similar to previous best
    ## Run 334 stress 0.06899842 
    ## ... Procrustes: rmse 0.0001650027  max resid 0.0003148456 
    ## ... Similar to previous best
    ## Run 335 stress 0.069637 
    ## Run 336 stress 0.06908164 
    ## ... Procrustes: rmse 0.01027617  max resid 0.029055 
    ## Run 337 stress 0.06899834 
    ## ... Procrustes: rmse 3.42306e-05  max resid 7.222842e-05 
    ## ... Similar to previous best
    ## Run 338 stress 0.06969201 
    ## Run 339 stress 0.06899837 
    ## ... Procrustes: rmse 0.0001089686  max resid 0.0001994146 
    ## ... Similar to previous best
    ## Run 340 stress 0.06969277 
    ## Run 341 stress 0.069637 
    ## Run 342 stress 0.06908227 
    ## ... Procrustes: rmse 0.01058385  max resid 0.0299029 
    ## Run 343 stress 0.06969263 
    ## Run 344 stress 0.06969273 
    ## Run 345 stress 0.06969204 
    ## Run 346 stress 0.06969217 
    ## Run 347 stress 0.06963704 
    ## Run 348 stress 0.06969263 
    ## Run 349 stress 0.06963712 
    ## Run 350 stress 0.06899838 
    ## ... Procrustes: rmse 0.0001228026  max resid 0.0002456001 
    ## ... Similar to previous best
    ## Run 351 stress 0.06969194 
    ## Run 352 stress 0.0689984 
    ## ... Procrustes: rmse 0.0001244898  max resid 0.0002252329 
    ## ... Similar to previous best
    ## Run 353 stress 0.06899837 
    ## ... Procrustes: rmse 0.0001097326  max resid 0.0002012647 
    ## ... Similar to previous best
    ## Run 354 stress 0.06908153 
    ## ... Procrustes: rmse 0.01019204  max resid 0.02882309 
    ## Run 355 stress 0.06969215 
    ## Run 356 stress 0.06963734 
    ## Run 357 stress 0.06908162 
    ## ... Procrustes: rmse 0.0102692  max resid 0.02903663 
    ## Run 358 stress 0.06899837 
    ## ... Procrustes: rmse 0.0001035974  max resid 0.0001897636 
    ## ... Similar to previous best
    ## Run 359 stress 0.06963705 
    ## Run 360 stress 0.06899846 
    ## ... Procrustes: rmse 0.0002095576  max resid 0.0004633988 
    ## ... Similar to previous best
    ## Run 361 stress 0.2612326 
    ## Run 362 stress 0.0690819 
    ## ... Procrustes: rmse 0.0104409  max resid 0.02951312 
    ## Run 363 stress 0.06969277 
    ## Run 364 stress 0.06899835 
    ## ... Procrustes: rmse 6.2214e-05  max resid 0.0001111828 
    ## ... Similar to previous best
    ## Run 365 stress 0.0696371 
    ## Run 366 stress 0.06969216 
    ## Run 367 stress 0.06963726 
    ## Run 368 stress 0.06908195 
    ## ... Procrustes: rmse 0.01046832  max resid 0.02958788 
    ## Run 369 stress 0.06908188 
    ## ... Procrustes: rmse 0.01039488  max resid 0.02936347 
    ## Run 370 stress 0.06899848 
    ## ... Procrustes: rmse 0.0002255522  max resid 0.0004301586 
    ## ... Similar to previous best
    ## Run 371 stress 0.3014963 
    ## Run 372 stress 0.06969256 
    ## Run 373 stress 0.0696922 
    ## Run 374 stress 0.06963718 
    ## Run 375 stress 0.06908109 
    ## ... Procrustes: rmse 0.009766162  max resid 0.02763799 
    ## Run 376 stress 0.06969274 
    ## Run 377 stress 0.06908084 
    ## ... Procrustes: rmse 0.009529107  max resid 0.02698632 
    ## Run 378 stress 0.06969321 
    ## Run 379 stress 0.06969258 
    ## Run 380 stress 0.06963702 
    ## Run 381 stress 0.06969241 
    ## Run 382 stress 0.06969207 
    ## Run 383 stress 0.06963718 
    ## Run 384 stress 0.06899841 
    ## ... Procrustes: rmse 0.0001621967  max resid 0.0003251832 
    ## ... Similar to previous best
    ## Run 385 stress 0.06963709 
    ## Run 386 stress 0.06899835 
    ## ... Procrustes: rmse 4.771578e-05  max resid 9.826881e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.06908179 
    ## ... Procrustes: rmse 0.01037964  max resid 0.02934261 
    ## Run 388 stress 0.06969198 
    ## Run 389 stress 0.06899838 
    ## ... Procrustes: rmse 0.0001197703  max resid 0.0002348153 
    ## ... Similar to previous best
    ## Run 390 stress 0.06963699 
    ## Run 391 stress 0.06899849 
    ## ... Procrustes: rmse 0.0001864591  max resid 0.0003476195 
    ## ... Similar to previous best
    ## Run 392 stress 0.06969248 
    ## Run 393 stress 0.06963703 
    ## Run 394 stress 0.06969255 
    ## Run 395 stress 0.06963713 
    ## Run 396 stress 0.06969249 
    ## Run 397 stress 0.06908177 
    ## ... Procrustes: rmse 0.01035265  max resid 0.02924605 
    ## Run 398 stress 0.06969259 
    ## Run 399 stress 0.06963704 
    ## Run 400 stress 0.06963733 
    ## Run 401 stress 0.3336295 
    ## Run 402 stress 0.06899837 
    ## ... Procrustes: rmse 9.766588e-05  max resid 0.0001732957 
    ## ... Similar to previous best
    ## Run 403 stress 0.06963702 
    ## Run 404 stress 0.06969221 
    ## Run 405 stress 0.0696372 
    ## Run 406 stress 0.06899839 
    ## ... Procrustes: rmse 6.955936e-05  max resid 0.0001225392 
    ## ... Similar to previous best
    ## Run 407 stress 0.06899835 
    ## ... Procrustes: rmse 6.104069e-05  max resid 0.0001473616 
    ## ... Similar to previous best
    ## Run 408 stress 0.06899844 
    ## ... Procrustes: rmse 0.0001893144  max resid 0.0003884856 
    ## ... Similar to previous best
    ## Run 409 stress 0.06969208 
    ## Run 410 stress 0.06899836 
    ## ... Procrustes: rmse 6.837561e-05  max resid 0.000124152 
    ## ... Similar to previous best
    ## Run 411 stress 0.06969274 
    ## Run 412 stress 0.06899843 
    ## ... Procrustes: rmse 0.0001710356  max resid 0.0003477734 
    ## ... Similar to previous best
    ## Run 413 stress 0.069693 
    ## Run 414 stress 0.06908153 
    ## ... Procrustes: rmse 0.01020047  max resid 0.02884533 
    ## Run 415 stress 0.069637 
    ## Run 416 stress 0.06969239 
    ## Run 417 stress 0.06963714 
    ## Run 418 stress 0.06899836 
    ## ... Procrustes: rmse 7.641839e-05  max resid 0.0001349724 
    ## ... Similar to previous best
    ## Run 419 stress 0.06908202 
    ## ... Procrustes: rmse 0.01051113  max resid 0.02970519 
    ## Run 420 stress 0.06899838 
    ## ... Procrustes: rmse 0.0001153483  max resid 0.0002285161 
    ## ... Similar to previous best
    ## Run 421 stress 0.06963711 
    ## Run 422 stress 0.06899836 
    ## ... Procrustes: rmse 7.290903e-05  max resid 0.0001343616 
    ## ... Similar to previous best
    ## Run 423 stress 0.06969226 
    ## Run 424 stress 0.3173128 
    ## Run 425 stress 0.06899839 
    ## ... Procrustes: rmse 0.000137682  max resid 0.0002834049 
    ## ... Similar to previous best
    ## Run 426 stress 0.06969241 
    ## Run 427 stress 0.06969295 
    ## Run 428 stress 0.06899834 
    ## ... Procrustes: rmse 1.593203e-05  max resid 3.996834e-05 
    ## ... Similar to previous best
    ## Run 429 stress 0.06899838 
    ## ... Procrustes: rmse 0.0001126116  max resid 0.0002019125 
    ## ... Similar to previous best
    ## Run 430 stress 0.06963708 
    ## Run 431 stress 0.06899846 
    ## ... Procrustes: rmse 0.0002050195  max resid 0.0003896101 
    ## ... Similar to previous best
    ## Run 432 stress 0.06969238 
    ## Run 433 stress 0.069082 
    ## ... Procrustes: rmse 0.01044532  max resid 0.02952984 
    ## Run 434 stress 0.06899842 
    ## ... Procrustes: rmse 0.0001652942  max resid 0.0003161136 
    ## ... Similar to previous best
    ## Run 435 stress 0.06969231 
    ## Run 436 stress 0.06963701 
    ## Run 437 stress 0.06969195 
    ## Run 438 stress 0.06963699 
    ## Run 439 stress 0.06963714 
    ## Run 440 stress 0.3187482 
    ## Run 441 stress 0.069692 
    ## Run 442 stress 0.0690815 
    ## ... Procrustes: rmse 0.01013425  max resid 0.02866349 
    ## Run 443 stress 0.06969257 
    ## Run 444 stress 0.06908203 
    ## ... Procrustes: rmse 0.01051935  max resid 0.02972909 
    ## Run 445 stress 0.06908179 
    ## ... Procrustes: rmse 0.01036667  max resid 0.02930868 
    ## Run 446 stress 0.06899834 
    ## ... Procrustes: rmse 2.297552e-05  max resid 5.297071e-05 
    ## ... Similar to previous best
    ## Run 447 stress 0.06969306 
    ## Run 448 stress 0.06963699 
    ## Run 449 stress 0.06899835 
    ## ... Procrustes: rmse 3.433458e-05  max resid 6.311111e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.06969209 
    ## Run 451 stress 0.06899845 
    ## ... Procrustes: rmse 0.0001860343  max resid 0.0003897171 
    ## ... Similar to previous best
    ## Run 452 stress 0.0696921 
    ## Run 453 stress 0.06963699 
    ## Run 454 stress 0.06969256 
    ## Run 455 stress 0.06908156 
    ## ... Procrustes: rmse 0.01022537  max resid 0.02891379 
    ## Run 456 stress 0.06969321 
    ## Run 457 stress 0.06899845 
    ## ... Procrustes: rmse 0.0002012989  max resid 0.0004014401 
    ## ... Similar to previous best
    ## Run 458 stress 0.06899835 
    ## ... Procrustes: rmse 2.702669e-05  max resid 6.352592e-05 
    ## ... Similar to previous best
    ## Run 459 stress 0.06969192 
    ## Run 460 stress 0.06963709 
    ## Run 461 stress 0.0696928 
    ## Run 462 stress 0.06963708 
    ## Run 463 stress 0.06899838 
    ## ... Procrustes: rmse 0.0001008612  max resid 0.0001837115 
    ## ... Similar to previous best
    ## Run 464 stress 0.06908122 
    ## ... Procrustes: rmse 0.009921302  max resid 0.02807035 
    ## Run 465 stress 0.06969228 
    ## Run 466 stress 0.06969293 
    ## Run 467 stress 0.06969273 
    ## Run 468 stress 0.06969247 
    ## Run 469 stress 0.06908191 
    ## ... Procrustes: rmse 0.01045306  max resid 0.029547 
    ## Run 470 stress 0.06969251 
    ## Run 471 stress 0.06969218 
    ## Run 472 stress 0.06908155 
    ## ... Procrustes: rmse 0.01021205  max resid 0.02887792 
    ## Run 473 stress 0.06899848 
    ## ... Procrustes: rmse 0.0001101302  max resid 0.0001895376 
    ## ... Similar to previous best
    ## Run 474 stress 0.06899838 
    ## ... Procrustes: rmse 0.0001112603  max resid 0.0002004136 
    ## ... Similar to previous best
    ## Run 475 stress 0.06899836 
    ## ... Procrustes: rmse 8.305733e-05  max resid 0.0001548122 
    ## ... Similar to previous best
    ## Run 476 stress 0.06969249 
    ## Run 477 stress 0.06908152 
    ## ... Procrustes: rmse 0.01018138  max resid 0.02879425 
    ## Run 478 stress 0.06899843 
    ## ... Procrustes: rmse 0.0001821894  max resid 0.0003923373 
    ## ... Similar to previous best
    ## Run 479 stress 0.06969257 
    ## Run 480 stress 0.06969264 
    ## Run 481 stress 0.353205 
    ## Run 482 stress 0.06908176 
    ## ... Procrustes: rmse 0.01030818  max resid 0.02915039 
    ## Run 483 stress 0.06908161 
    ## ... Procrustes: rmse 0.01021012  max resid 0.02886856 
    ## Run 484 stress 0.06908143 
    ## ... Procrustes: rmse 0.01011651  max resid 0.02859829 
    ## Run 485 stress 0.06969218 
    ## Run 486 stress 0.0696929 
    ## Run 487 stress 0.06899846 
    ## ... Procrustes: rmse 0.0002052403  max resid 0.0004059703 
    ## ... Similar to previous best
    ## Run 488 stress 0.06899835 
    ## ... Procrustes: rmse 4.926099e-05  max resid 0.0001061352 
    ## ... Similar to previous best
    ## Run 489 stress 0.06899836 
    ## ... Procrustes: rmse 9.308798e-05  max resid 0.0001879785 
    ## ... Similar to previous best
    ## Run 490 stress 0.06963713 
    ## Run 491 stress 0.06969216 
    ## Run 492 stress 0.06969233 
    ## Run 493 stress 0.06908209 
    ## ... Procrustes: rmse 0.01054114  max resid 0.0297924 
    ## Run 494 stress 0.06899847 
    ## ... Procrustes: rmse 0.0002110141  max resid 0.0004141163 
    ## ... Similar to previous best
    ## Run 495 stress 0.06969269 
    ## Run 496 stress 0.06908155 
    ## ... Procrustes: rmse 0.01020683  max resid 0.02886454 
    ## Run 497 stress 0.069693 
    ## Run 498 stress 0.06963699 
    ## Run 499 stress 0.06969194 
    ## Run 500 stress 0.06908077 
    ## ... Procrustes: rmse 0.009460901  max resid 0.02677997 
    ## *** Best solution repeated 98 times

``` r
# Stratified lakes and ocean sites
SD_beta_geo_SO_NMDS <- metaMDS(SD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.06978199 
    ## ... New best solution
    ## ... Procrustes: rmse 0.09984319  max resid 0.2618872 
    ## Run 2 stress 0.08340288 
    ## Run 3 stress 0.06978189 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001677749  max resid 0.0004305126 
    ## ... Similar to previous best
    ## Run 4 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01322097  max resid 0.03349804 
    ## Run 5 stress 0.08340288 
    ## Run 6 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323681  max resid 0.0332945 
    ## Run 7 stress 0.07428319 
    ## Run 8 stress 0.07250814 
    ## Run 9 stress 0.08340295 
    ## Run 10 stress 0.06942776 
    ## ... Procrustes: rmse 1.006604e-05  max resid 2.60222e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.07428313 
    ## Run 12 stress 0.069782 
    ## ... Procrustes: rmse 0.01328149  max resid 0.0333883 
    ## Run 13 stress 0.06942777 
    ## ... Procrustes: rmse 3.68242e-05  max resid 9.415146e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.06942777 
    ## ... Procrustes: rmse 2.90952e-05  max resid 7.734092e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.08340286 
    ## Run 16 stress 0.07970525 
    ## Run 17 stress 0.06942777 
    ## ... Procrustes: rmse 1.506615e-05  max resid 4.334692e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.08340288 
    ## Run 19 stress 0.07250812 
    ## Run 20 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132382  max resid 0.03329695 
    ## Run 21 stress 0.08340286 
    ## Run 22 stress 0.08340289 
    ## Run 23 stress 0.08340287 
    ## Run 24 stress 0.07250814 
    ## Run 25 stress 0.07970526 
    ## Run 26 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325459  max resid 0.03333072 
    ## Run 27 stress 0.08340288 
    ## Run 28 stress 0.07428314 
    ## Run 29 stress 0.07844918 
    ## Run 30 stress 0.07970528 
    ## Run 31 stress 0.07970525 
    ## Run 32 stress 0.08340287 
    ## Run 33 stress 0.07428314 
    ## Run 34 stress 0.08448461 
    ## Run 35 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321347  max resid 0.03324642 
    ## Run 36 stress 0.06942777 
    ## ... Procrustes: rmse 2.520797e-05  max resid 6.674985e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.07970525 
    ## Run 38 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323151  max resid 0.03328313 
    ## Run 39 stress 0.07970525 
    ## Run 40 stress 0.08448444 
    ## Run 41 stress 0.08448439 
    ## Run 42 stress 0.06978192 
    ## ... Procrustes: rmse 0.0132517  max resid 0.03332632 
    ## Run 43 stress 0.07970526 
    ## Run 44 stress 0.07428313 
    ## Run 45 stress 0.07428313 
    ## Run 46 stress 0.0834029 
    ## Run 47 stress 0.08340295 
    ## Run 48 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324379  max resid 0.03331031 
    ## Run 49 stress 0.07970525 
    ## Run 50 stress 0.06942776 
    ## ... Procrustes: rmse 1.535006e-06  max resid 3.295422e-06 
    ## ... Similar to previous best
    ## Run 51 stress 0.08340289 
    ## Run 52 stress 0.08340298 
    ## Run 53 stress 0.08340294 
    ## Run 54 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325492  max resid 0.03333063 
    ## Run 55 stress 0.07250812 
    ## Run 56 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324323  max resid 0.03330777 
    ## Run 57 stress 0.07428318 
    ## Run 58 stress 0.0834029 
    ## Run 59 stress 0.07250812 
    ## Run 60 stress 0.06942776 
    ## ... Procrustes: rmse 1.289621e-05  max resid 3.213254e-05 
    ## ... Similar to previous best
    ## Run 61 stress 0.07428314 
    ## Run 62 stress 0.07250812 
    ## Run 63 stress 0.07428316 
    ## Run 64 stress 0.07970525 
    ## Run 65 stress 0.07970526 
    ## Run 66 stress 0.07250812 
    ## Run 67 stress 0.07250812 
    ## Run 68 stress 0.07428318 
    ## Run 69 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325706  max resid 0.03333418 
    ## Run 70 stress 0.07428318 
    ## Run 71 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325934  max resid 0.03334066 
    ## Run 72 stress 0.07250819 
    ## Run 73 stress 0.07970525 
    ## Run 74 stress 0.06942777 
    ## ... Procrustes: rmse 3.341445e-05  max resid 8.552684e-05 
    ## ... Similar to previous best
    ## Run 75 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323217  max resid 0.03328379 
    ## Run 76 stress 0.07428313 
    ## Run 77 stress 0.07250813 
    ## Run 78 stress 0.07428314 
    ## Run 79 stress 0.07250812 
    ## Run 80 stress 0.07428314 
    ## Run 81 stress 0.07970525 
    ## Run 82 stress 0.07428312 
    ## Run 83 stress 0.07428313 
    ## Run 84 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323303  max resid 0.03328821 
    ## Run 85 stress 0.07428313 
    ## Run 86 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132093  max resid 0.03323711 
    ## Run 87 stress 0.06942777 
    ## ... Procrustes: rmse 1.133874e-05  max resid 2.203287e-05 
    ## ... Similar to previous best
    ## Run 88 stress 0.07250812 
    ## Run 89 stress 0.08448448 
    ## Run 90 stress 0.08340287 
    ## Run 91 stress 0.0834029 
    ## Run 92 stress 0.06942778 
    ## ... Procrustes: rmse 4.457435e-05  max resid 0.0001174445 
    ## ... Similar to previous best
    ## Run 93 stress 0.07844949 
    ## Run 94 stress 0.07428313 
    ## Run 95 stress 0.2666427 
    ## Run 96 stress 0.06942776 
    ## ... Procrustes: rmse 2.667158e-05  max resid 6.893935e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.08340298 
    ## Run 98 stress 0.07428317 
    ## Run 99 stress 0.0834029 
    ## Run 100 stress 0.07428313 
    ## Run 101 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325719  max resid 0.03333732 
    ## Run 102 stress 0.06942776 
    ## ... Procrustes: rmse 5.122717e-06  max resid 1.410615e-05 
    ## ... Similar to previous best
    ## Run 103 stress 0.06978192 
    ## ... Procrustes: rmse 0.0132512  max resid 0.03332542 
    ## Run 104 stress 0.07428315 
    ## Run 105 stress 0.08448443 
    ## Run 106 stress 0.06942778 
    ## ... Procrustes: rmse 5.114935e-05  max resid 0.0001319363 
    ## ... Similar to previous best
    ## Run 107 stress 0.06942777 
    ## ... Procrustes: rmse 3.13175e-05  max resid 8.001607e-05 
    ## ... Similar to previous best
    ## Run 108 stress 0.07970526 
    ## Run 109 stress 0.06942776 
    ## ... Procrustes: rmse 1.10466e-05  max resid 1.633487e-05 
    ## ... Similar to previous best
    ## Run 110 stress 0.07250812 
    ## Run 111 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324839  max resid 0.03331945 
    ## Run 112 stress 0.07428314 
    ## Run 113 stress 0.06942777 
    ## ... Procrustes: rmse 5.181716e-05  max resid 0.0001341365 
    ## ... Similar to previous best
    ## Run 114 stress 0.07428313 
    ## Run 115 stress 0.07844954 
    ## Run 116 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318891  max resid 0.03319236 
    ## Run 117 stress 0.08340291 
    ## Run 118 stress 0.07844922 
    ## Run 119 stress 0.07428317 
    ## Run 120 stress 0.07250815 
    ## Run 121 stress 0.08340289 
    ## Run 122 stress 0.07970525 
    ## Run 123 stress 0.08340288 
    ## Run 124 stress 0.07428315 
    ## Run 125 stress 0.06942777 
    ## ... Procrustes: rmse 3.430246e-05  max resid 8.797808e-05 
    ## ... Similar to previous best
    ## Run 126 stress 0.06978199 
    ## ... Procrustes: rmse 0.01316541  max resid 0.03313977 
    ## Run 127 stress 0.07428317 
    ## Run 128 stress 0.08340296 
    ## Run 129 stress 0.07428316 
    ## Run 130 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325314  max resid 0.03332769 
    ## Run 131 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 5.836151e-06  max resid 1.493671e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.06978196 
    ## ... Procrustes: rmse 0.01327154  max resid 0.03336649 
    ## Run 133 stress 0.07428313 
    ## Run 134 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326716  max resid 0.0333556 
    ## Run 135 stress 0.07428316 
    ## Run 136 stress 0.07970526 
    ## Run 137 stress 0.06942777 
    ## ... Procrustes: rmse 4.895903e-05  max resid 0.0001260029 
    ## ... Similar to previous best
    ## Run 138 stress 0.07250812 
    ## Run 139 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326318  max resid 0.03334854 
    ## Run 140 stress 0.08340291 
    ## Run 141 stress 0.06942776 
    ## ... Procrustes: rmse 8.970035e-06  max resid 2.443536e-05 
    ## ... Similar to previous best
    ## Run 142 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323402  max resid 0.0332861 
    ## Run 143 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324141  max resid 0.03330411 
    ## Run 144 stress 0.07428314 
    ## Run 145 stress 0.07970525 
    ## Run 146 stress 0.06978194 
    ## ... Procrustes: rmse 0.01318231  max resid 0.03317443 
    ## Run 147 stress 0.07970525 
    ## Run 148 stress 0.07428321 
    ## Run 149 stress 0.08340291 
    ## Run 150 stress 0.08340289 
    ## Run 151 stress 0.08448448 
    ## Run 152 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325219  max resid 0.03332425 
    ## Run 153 stress 0.0784495 
    ## Run 154 stress 0.07428313 
    ## Run 155 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326108  max resid 0.03334007 
    ## Run 156 stress 0.08340287 
    ## Run 157 stress 0.06942776 
    ## ... Procrustes: rmse 1.198872e-05  max resid 3.096201e-05 
    ## ... Similar to previous best
    ## Run 158 stress 0.07250815 
    ## Run 159 stress 0.07428313 
    ## Run 160 stress 0.08340287 
    ## Run 161 stress 0.08340295 
    ## Run 162 stress 0.07250813 
    ## Run 163 stress 0.07250815 
    ## Run 164 stress 0.07250815 
    ## Run 165 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132367  max resid 0.03329208 
    ## Run 166 stress 0.07428317 
    ## Run 167 stress 0.07250812 
    ## Run 168 stress 0.06942776 
    ## ... Procrustes: rmse 1.863301e-06  max resid 3.891056e-06 
    ## ... Similar to previous best
    ## Run 169 stress 0.08448437 
    ## Run 170 stress 0.06978199 
    ## ... Procrustes: rmse 0.01328036  max resid 0.03338099 
    ## Run 171 stress 0.07428313 
    ## Run 172 stress 0.07250812 
    ## Run 173 stress 0.07428318 
    ## Run 174 stress 0.07970526 
    ## Run 175 stress 0.07250812 
    ## Run 176 stress 0.06942776 
    ## ... Procrustes: rmse 2.724589e-05  max resid 6.943934e-05 
    ## ... Similar to previous best
    ## Run 177 stress 0.06942776 
    ## ... Procrustes: rmse 6.378871e-06  max resid 1.638844e-05 
    ## ... Similar to previous best
    ## Run 178 stress 0.06942776 
    ## ... Procrustes: rmse 1.871082e-05  max resid 4.814863e-05 
    ## ... Similar to previous best
    ## Run 179 stress 0.07970526 
    ## Run 180 stress 0.08448454 
    ## Run 181 stress 0.07970525 
    ## Run 182 stress 0.07428316 
    ## Run 183 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325543  max resid 0.03333121 
    ## Run 184 stress 0.08340292 
    ## Run 185 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327921  max resid 0.03337931 
    ## Run 186 stress 0.07844957 
    ## Run 187 stress 0.08448455 
    ## Run 188 stress 0.08340291 
    ## Run 189 stress 0.07428313 
    ## Run 190 stress 0.07250812 
    ## Run 191 stress 0.07428314 
    ## Run 192 stress 0.07250812 
    ## Run 193 stress 0.07250812 
    ## Run 194 stress 0.07428313 
    ## Run 195 stress 0.07970525 
    ## Run 196 stress 0.08340291 
    ## Run 197 stress 0.06942776 
    ## ... Procrustes: rmse 9.969883e-06  max resid 2.58753e-05 
    ## ... Similar to previous best
    ## Run 198 stress 0.07428319 
    ## Run 199 stress 0.08340286 
    ## Run 200 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324886  max resid 0.03331725 
    ## Run 201 stress 0.07844938 
    ## Run 202 stress 0.06942776 
    ## ... Procrustes: rmse 1.780836e-05  max resid 4.578415e-05 
    ## ... Similar to previous best
    ## Run 203 stress 0.07428313 
    ## Run 204 stress 0.06942776 
    ## ... Procrustes: rmse 2.166511e-06  max resid 6.225289e-06 
    ## ... Similar to previous best
    ## Run 205 stress 0.08340295 
    ## Run 206 stress 0.07250812 
    ## Run 207 stress 0.07844934 
    ## Run 208 stress 0.07250813 
    ## Run 209 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321353  max resid 0.0332425 
    ## Run 210 stress 0.07970527 
    ## Run 211 stress 0.07250815 
    ## Run 212 stress 0.07428313 
    ## Run 213 stress 0.07250814 
    ## Run 214 stress 0.07250815 
    ## Run 215 stress 0.06942779 
    ## ... Procrustes: rmse 4.097955e-05  max resid 0.0001123742 
    ## ... Similar to previous best
    ## Run 216 stress 0.0742832 
    ## Run 217 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322276  max resid 0.03325998 
    ## Run 218 stress 0.06942777 
    ## ... Procrustes: rmse 5.287765e-05  max resid 0.0001359543 
    ## ... Similar to previous best
    ## Run 219 stress 0.07970525 
    ## Run 220 stress 0.07428315 
    ## Run 221 stress 0.08448446 
    ## Run 222 stress 0.0844844 
    ## Run 223 stress 0.07428313 
    ## Run 224 stress 0.07970525 
    ## Run 225 stress 0.0697819 
    ## ... Procrustes: rmse 0.01324118  max resid 0.03330198 
    ## Run 226 stress 0.07970525 
    ## Run 227 stress 0.06942776 
    ## ... Procrustes: rmse 2.950052e-05  max resid 7.602901e-05 
    ## ... Similar to previous best
    ## Run 228 stress 0.07844967 
    ## Run 229 stress 0.07970526 
    ## Run 230 stress 0.07250815 
    ## Run 231 stress 0.07428316 
    ## Run 232 stress 0.07844941 
    ## Run 233 stress 0.06942776 
    ## ... Procrustes: rmse 1.712322e-05  max resid 4.456489e-05 
    ## ... Similar to previous best
    ## Run 234 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325161  max resid 0.03332365 
    ## Run 235 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323524  max resid 0.03328863 
    ## Run 236 stress 0.07970525 
    ## Run 237 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322048  max resid 0.03325693 
    ## Run 238 stress 0.07970526 
    ## Run 239 stress 0.07428316 
    ## Run 240 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322105  max resid 0.03325865 
    ## Run 241 stress 0.06942777 
    ## ... Procrustes: rmse 4.750211e-05  max resid 0.0001221801 
    ## ... Similar to previous best
    ## Run 242 stress 0.07428313 
    ## Run 243 stress 0.07428315 
    ## Run 244 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324772  max resid 0.03331473 
    ## Run 245 stress 0.07970526 
    ## Run 246 stress 0.07250813 
    ## Run 247 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324669  max resid 0.0333123 
    ## Run 248 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324558  max resid 0.03330981 
    ## Run 249 stress 0.07250812 
    ## Run 250 stress 0.07250812 
    ## Run 251 stress 0.08340299 
    ## Run 252 stress 0.07428317 
    ## Run 253 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319987  max resid 0.03321278 
    ## Run 254 stress 0.08340286 
    ## Run 255 stress 0.069782 
    ## ... Procrustes: rmse 0.01328219  max resid 0.03338559 
    ## Run 256 stress 0.07970526 
    ## Run 257 stress 0.07844923 
    ## Run 258 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326653  max resid 0.03334921 
    ## Run 259 stress 0.07970525 
    ## Run 260 stress 0.06942776 
    ## ... Procrustes: rmse 1.416633e-05  max resid 3.683025e-05 
    ## ... Similar to previous best
    ## Run 261 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323662  max resid 0.03329168 
    ## Run 262 stress 0.07970526 
    ## Run 263 stress 0.08340295 
    ## Run 264 stress 0.07428316 
    ## Run 265 stress 0.06942778 
    ## ... Procrustes: rmse 5.925525e-05  max resid 0.000152636 
    ## ... Similar to previous best
    ## Run 266 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325545  max resid 0.03333045 
    ## Run 267 stress 0.07970525 
    ## Run 268 stress 0.07250812 
    ## Run 269 stress 0.07428314 
    ## Run 270 stress 0.06942777 
    ## ... Procrustes: rmse 3.742785e-05  max resid 9.656716e-05 
    ## ... Similar to previous best
    ## Run 271 stress 0.06942776 
    ## ... Procrustes: rmse 1.695523e-05  max resid 4.36653e-05 
    ## ... Similar to previous best
    ## Run 272 stress 0.07250813 
    ## Run 273 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326719  max resid 0.03335254 
    ## Run 274 stress 0.07970526 
    ## Run 275 stress 0.08340289 
    ## Run 276 stress 0.07970525 
    ## Run 277 stress 0.07844935 
    ## Run 278 stress 0.06942776 
    ## ... Procrustes: rmse 1.572398e-06  max resid 3.88304e-06 
    ## ... Similar to previous best
    ## Run 279 stress 0.08340287 
    ## Run 280 stress 0.08340286 
    ## Run 281 stress 0.07250812 
    ## Run 282 stress 0.06942776 
    ## ... Procrustes: rmse 7.839553e-06  max resid 2.103466e-05 
    ## ... Similar to previous best
    ## Run 283 stress 0.06942776 
    ## ... Procrustes: rmse 1.369078e-06  max resid 2.987463e-06 
    ## ... Similar to previous best
    ## Run 284 stress 0.07250815 
    ## Run 285 stress 0.07250813 
    ## Run 286 stress 0.07250812 
    ## Run 287 stress 0.07970526 
    ## Run 288 stress 0.06942776 
    ## ... Procrustes: rmse 2.695925e-05  max resid 7.047052e-05 
    ## ... Similar to previous best
    ## Run 289 stress 0.07250813 
    ## Run 290 stress 0.07428313 
    ## Run 291 stress 0.06942776 
    ## ... Procrustes: rmse 7.54591e-06  max resid 1.953864e-05 
    ## ... Similar to previous best
    ## Run 292 stress 0.08340289 
    ## Run 293 stress 0.07250813 
    ## Run 294 stress 0.07970526 
    ## Run 295 stress 0.07428313 
    ## Run 296 stress 0.07250812 
    ## Run 297 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323907  max resid 0.03330073 
    ## Run 298 stress 0.07970526 
    ## Run 299 stress 0.07844917 
    ## Run 300 stress 0.07428321 
    ## Run 301 stress 0.08340291 
    ## Run 302 stress 0.07970525 
    ## Run 303 stress 0.07250812 
    ## Run 304 stress 0.07428313 
    ## Run 305 stress 0.07970525 
    ## Run 306 stress 0.07970526 
    ## Run 307 stress 0.07428317 
    ## Run 308 stress 0.08448451 
    ## Run 309 stress 0.07428314 
    ## Run 310 stress 0.08340293 
    ## Run 311 stress 0.07970525 
    ## Run 312 stress 0.07250812 
    ## Run 313 stress 0.07428317 
    ## Run 314 stress 0.08340289 
    ## Run 315 stress 0.06942776 
    ## ... Procrustes: rmse 2.585013e-05  max resid 6.664275e-05 
    ## ... Similar to previous best
    ## Run 316 stress 0.0834029 
    ## Run 317 stress 0.07970525 
    ## Run 318 stress 0.07428313 
    ## Run 319 stress 0.07970526 
    ## Run 320 stress 0.06942776 
    ## ... Procrustes: rmse 6.139874e-06  max resid 9.165004e-06 
    ## ... Similar to previous best
    ## Run 321 stress 0.07428316 
    ## Run 322 stress 0.08340286 
    ## Run 323 stress 0.06942776 
    ## ... Procrustes: rmse 1.052954e-05  max resid 2.709369e-05 
    ## ... Similar to previous best
    ## Run 324 stress 0.07250813 
    ## Run 325 stress 0.08340287 
    ## Run 326 stress 0.07970525 
    ## Run 327 stress 0.0784493 
    ## Run 328 stress 0.08340286 
    ## Run 329 stress 0.07428312 
    ## Run 330 stress 0.06942776 
    ## ... Procrustes: rmse 5.234238e-06  max resid 1.137628e-05 
    ## ... Similar to previous best
    ## Run 331 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326807  max resid 0.03335616 
    ## Run 332 stress 0.06978194 
    ## ... Procrustes: rmse 0.01318422  max resid 0.03317925 
    ## Run 333 stress 0.08448438 
    ## Run 334 stress 0.07428322 
    ## Run 335 stress 0.08340298 
    ## Run 336 stress 0.08448452 
    ## Run 337 stress 0.06978196 
    ## ... Procrustes: rmse 0.01327064  max resid 0.03336136 
    ## Run 338 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326433  max resid 0.03335175 
    ## Run 339 stress 0.06942776 
    ## ... Procrustes: rmse 2.765419e-05  max resid 7.134086e-05 
    ## ... Similar to previous best
    ## Run 340 stress 0.07970526 
    ## Run 341 stress 0.07428313 
    ## Run 342 stress 0.07250812 
    ## Run 343 stress 0.07970526 
    ## Run 344 stress 0.06942776 
    ## ... Procrustes: rmse 2.049617e-05  max resid 5.30848e-05 
    ## ... Similar to previous best
    ## Run 345 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325767  max resid 0.03333881 
    ## Run 346 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323617  max resid 0.03329279 
    ## Run 347 stress 0.07428316 
    ## Run 348 stress 0.07428314 
    ## Run 349 stress 0.08340292 
    ## Run 350 stress 0.07250813 
    ## Run 351 stress 0.07428315 
    ## Run 352 stress 0.08448437 
    ## Run 353 stress 0.06978201 
    ## ... Procrustes: rmse 0.01327272  max resid 0.03336124 
    ## Run 354 stress 0.07970526 
    ## Run 355 stress 0.08340287 
    ## Run 356 stress 0.08340288 
    ## Run 357 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323938  max resid 0.03330229 
    ## Run 358 stress 0.08340286 
    ## Run 359 stress 0.08340287 
    ## Run 360 stress 0.07250812 
    ## Run 361 stress 0.08340294 
    ## Run 362 stress 0.07428313 
    ## Run 363 stress 0.07970525 
    ## Run 364 stress 0.07428313 
    ## Run 365 stress 0.07250813 
    ## Run 366 stress 0.06942776 
    ## ... Procrustes: rmse 2.392486e-05  max resid 6.184051e-05 
    ## ... Similar to previous best
    ## Run 367 stress 0.08340287 
    ## Run 368 stress 0.06942776 
    ## ... Procrustes: rmse 5.247701e-06  max resid 1.000592e-05 
    ## ... Similar to previous best
    ## Run 369 stress 0.07250812 
    ## Run 370 stress 0.0784496 
    ## Run 371 stress 0.08448462 
    ## Run 372 stress 0.08340289 
    ## Run 373 stress 0.06942776 
    ## ... Procrustes: rmse 6.499763e-06  max resid 1.874993e-05 
    ## ... Similar to previous best
    ## Run 374 stress 0.07428312 
    ## Run 375 stress 0.07428312 
    ## Run 376 stress 0.07428313 
    ## Run 377 stress 0.08448436 
    ## Run 378 stress 0.07970525 
    ## Run 379 stress 0.07428317 
    ## Run 380 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321488  max resid 0.03324572 
    ## Run 381 stress 0.06942778 
    ## ... Procrustes: rmse 1.790872e-05  max resid 5.741621e-05 
    ## ... Similar to previous best
    ## Run 382 stress 0.06942779 
    ## ... Procrustes: rmse 7.13013e-05  max resid 0.00018418 
    ## ... Similar to previous best
    ## Run 383 stress 0.06942776 
    ## ... Procrustes: rmse 1.996331e-06  max resid 5.181281e-06 
    ## ... Similar to previous best
    ## Run 384 stress 0.08340288 
    ## Run 385 stress 0.06978201 
    ## ... Procrustes: rmse 0.01316037  max resid 0.03312656 
    ## Run 386 stress 0.06942777 
    ## ... Procrustes: rmse 2.516736e-05  max resid 7.102825e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325612  max resid 0.03333168 
    ## Run 388 stress 0.08340286 
    ## Run 389 stress 0.08448448 
    ## Run 390 stress 0.07844949 
    ## Run 391 stress 0.07250812 
    ## Run 392 stress 0.07428313 
    ## Run 393 stress 0.07250816 
    ## Run 394 stress 0.07428313 
    ## Run 395 stress 0.07250816 
    ## Run 396 stress 0.083403 
    ## Run 397 stress 0.08340299 
    ## Run 398 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324714  max resid 0.03331586 
    ## Run 399 stress 0.08340286 
    ## Run 400 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324576  max resid 0.03331449 
    ## Run 401 stress 0.07428315 
    ## Run 402 stress 0.07428316 
    ## Run 403 stress 0.06942776 
    ## ... Procrustes: rmse 5.554057e-06  max resid 1.407695e-05 
    ## ... Similar to previous best
    ## Run 404 stress 0.07970526 
    ## Run 405 stress 0.07844915 
    ## Run 406 stress 0.07250812 
    ## Run 407 stress 0.08340286 
    ## Run 408 stress 0.07250812 
    ## Run 409 stress 0.07970526 
    ## Run 410 stress 0.08340287 
    ## Run 411 stress 0.08340291 
    ## Run 412 stress 0.08448452 
    ## Run 413 stress 0.07428313 
    ## Run 414 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323834  max resid 0.03329398 
    ## Run 415 stress 0.08448438 
    ## Run 416 stress 0.08340289 
    ## Run 417 stress 0.08340288 
    ## Run 418 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323413  max resid 0.03328641 
    ## Run 419 stress 0.08340288 
    ## Run 420 stress 0.07428315 
    ## Run 421 stress 0.08448447 
    ## Run 422 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326033  max resid 0.03334303 
    ## Run 423 stress 0.08340287 
    ## Run 424 stress 0.07970526 
    ## Run 425 stress 0.07970526 
    ## Run 426 stress 0.07250814 
    ## Run 427 stress 0.06942776 
    ## ... Procrustes: rmse 1.597155e-05  max resid 4.100845e-05 
    ## ... Similar to previous best
    ## Run 428 stress 0.07428313 
    ## Run 429 stress 0.08448436 
    ## Run 430 stress 0.08448436 
    ## Run 431 stress 0.06942776 
    ## ... Procrustes: rmse 2.784724e-05  max resid 7.174064e-05 
    ## ... Similar to previous best
    ## Run 432 stress 0.07970525 
    ## Run 433 stress 0.07428313 
    ## Run 434 stress 0.07428313 
    ## Run 435 stress 0.07428313 
    ## Run 436 stress 0.07428317 
    ## Run 437 stress 0.08340293 
    ## Run 438 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327524  max resid 0.03337016 
    ## Run 439 stress 0.06978201 
    ## ... Procrustes: rmse 0.01328377  max resid 0.03338468 
    ## Run 440 stress 0.07250812 
    ## Run 441 stress 0.07250813 
    ## Run 442 stress 0.08340292 
    ## Run 443 stress 0.06978192 
    ## ... Procrustes: rmse 0.0131944  max resid 0.03320071 
    ## Run 444 stress 0.06942776 
    ## ... Procrustes: rmse 2.452755e-05  max resid 6.307977e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.07428316 
    ## Run 446 stress 0.07428313 
    ## Run 447 stress 0.06942778 
    ## ... Procrustes: rmse 6.122699e-05  max resid 0.0001581899 
    ## ... Similar to previous best
    ## Run 448 stress 0.08340287 
    ## Run 449 stress 0.06942776 
    ## ... Procrustes: rmse 7.842895e-06  max resid 2.019703e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.0844844 
    ## Run 451 stress 0.0834029 
    ## Run 452 stress 0.06978196 
    ## ... Procrustes: rmse 0.01327089  max resid 0.03336428 
    ## Run 453 stress 0.07250813 
    ## Run 454 stress 0.07970525 
    ## Run 455 stress 0.07250812 
    ## Run 456 stress 0.08340288 
    ## Run 457 stress 0.06942776 
    ## ... Procrustes: rmse 7.870654e-06  max resid 2.023445e-05 
    ## ... Similar to previous best
    ## Run 458 stress 0.06942777 
    ## ... Procrustes: rmse 4.370328e-05  max resid 0.0001123834 
    ## ... Similar to previous best
    ## Run 459 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319117  max resid 0.0331953 
    ## Run 460 stress 0.07428313 
    ## Run 461 stress 0.07844956 
    ## Run 462 stress 0.07970525 
    ## Run 463 stress 0.07428313 
    ## Run 464 stress 0.07428319 
    ## Run 465 stress 0.07970526 
    ## Run 466 stress 0.07428315 
    ## Run 467 stress 0.07428313 
    ## Run 468 stress 0.07428316 
    ## Run 469 stress 0.06978197 
    ## ... Procrustes: rmse 0.01324814  max resid 0.03332087 
    ## Run 470 stress 0.07970525 
    ## Run 471 stress 0.07250812 
    ## Run 472 stress 0.07428314 
    ## Run 473 stress 0.07428314 
    ## Run 474 stress 0.06942776 
    ## ... Procrustes: rmse 1.820864e-05  max resid 4.686068e-05 
    ## ... Similar to previous best
    ## Run 475 stress 0.07970525 
    ## Run 476 stress 0.08340289 
    ## Run 477 stress 0.07250815 
    ## Run 478 stress 0.07970525 
    ## Run 479 stress 0.06978192 
    ## ... Procrustes: rmse 0.01323624  max resid 0.03329442 
    ## Run 480 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321889  max resid 0.03324896 
    ## Run 481 stress 0.08340297 
    ## Run 482 stress 0.07428315 
    ## Run 483 stress 0.06942776 
    ## ... Procrustes: rmse 2.176325e-05  max resid 5.663777e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.06942778 
    ## ... Procrustes: rmse 2.024627e-05  max resid 6.121299e-05 
    ## ... Similar to previous best
    ## Run 485 stress 0.06942776 
    ## ... Procrustes: rmse 1.825816e-06  max resid 4.258577e-06 
    ## ... Similar to previous best
    ## Run 486 stress 0.08340291 
    ## Run 487 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322283  max resid 0.03325862 
    ## Run 488 stress 0.07428315 
    ## Run 489 stress 0.08340286 
    ## Run 490 stress 0.06942776 
    ## ... Procrustes: rmse 3.317963e-05  max resid 8.553084e-05 
    ## ... Similar to previous best
    ## Run 491 stress 0.0742832 
    ## Run 492 stress 0.06942776 
    ## ... Procrustes: rmse 3.418102e-05  max resid 8.794079e-05 
    ## ... Similar to previous best
    ## Run 493 stress 0.06942777 
    ## ... Procrustes: rmse 2.210453e-05  max resid 6.079599e-05 
    ## ... Similar to previous best
    ## Run 494 stress 0.07970525 
    ## Run 495 stress 0.08340287 
    ## Run 496 stress 0.08340292 
    ## Run 497 stress 0.0742832 
    ## Run 498 stress 0.06942776 
    ## ... Procrustes: rmse 1.478338e-05  max resid 4.01011e-05 
    ## ... Similar to previous best
    ## Run 499 stress 0.069782 
    ## ... Procrustes: rmse 0.01323368  max resid 0.03328012 
    ## Run 500 stress 0.0844844 
    ## *** Best solution repeated 54 times

``` r
# Mixed lakes
SD_beta_geo_M_NMDS <- metaMDS(SD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.994289e-05 
    ## Run 1 stress 0.2290616 
    ## Run 2 stress 0.157905 
    ## Run 3 stress 9.101415e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01992343  max resid 0.0299692 
    ## Run 4 stress 7.801362e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1639863  max resid 0.2245242 
    ## Run 5 stress 0.0004385328 
    ## ... Procrustes: rmse 0.02946319  max resid 0.06029971 
    ## Run 6 stress 0.2264602 
    ## Run 7 stress 9.816578e-05 
    ## ... Procrustes: rmse 0.17387  max resid 0.2637629 
    ## Run 8 stress 0.293122 
    ## Run 9 stress 9.084845e-05 
    ## ... Procrustes: rmse 0.1982669  max resid 0.3074761 
    ## Run 10 stress 9.257352e-05 
    ## ... Procrustes: rmse 0.02801196  max resid 0.03814212 
    ## Run 11 stress 9.798821e-05 
    ## ... Procrustes: rmse 0.2289713  max resid 0.3657448 
    ## Run 12 stress 9.821032e-05 
    ## ... Procrustes: rmse 0.02799329  max resid 0.03805229 
    ## Run 13 stress 9.06845e-05 
    ## ... Procrustes: rmse 0.02344944  max resid 0.03198878 
    ## Run 14 stress 8.795876e-05 
    ## ... Procrustes: rmse 0.2197221  max resid 0.3478398 
    ## Run 15 stress 9.78915e-05 
    ## ... Procrustes: rmse 0.02788429  max resid 0.03802805 
    ## Run 16 stress 9.226611e-05 
    ## ... Procrustes: rmse 0.06821332  max resid 0.09571956 
    ## Run 17 stress 9.211331e-05 
    ## ... Procrustes: rmse 0.0279963  max resid 0.03807368 
    ## Run 18 stress 9.3206e-05 
    ## ... Procrustes: rmse 0.0279849  max resid 0.03811941 
    ## Run 19 stress 9.354975e-05 
    ## ... Procrustes: rmse 0.2170104  max resid 0.3426513 
    ## Run 20 stress 8.965876e-05 
    ## ... Procrustes: rmse 0.0279999  max resid 0.03808309 
    ## Run 21 stress 8.642283e-05 
    ## ... Procrustes: rmse 0.02663386  max resid 0.03680385 
    ## Run 22 stress 0.2124849 
    ## Run 23 stress 8.868705e-05 
    ## ... Procrustes: rmse 0.2014412  max resid 0.3133108 
    ## Run 24 stress 0.2569849 
    ## Run 25 stress 0.2264602 
    ## Run 26 stress 8.745759e-05 
    ## ... Procrustes: rmse 0.2258451  max resid 0.3597172 
    ## Run 27 stress 9.979137e-05 
    ## ... Procrustes: rmse 0.02788193  max resid 0.03802394 
    ## Run 28 stress 9.521645e-05 
    ## ... Procrustes: rmse 0.02799094  max resid 0.03806174 
    ## Run 29 stress 9.322177e-05 
    ## ... Procrustes: rmse 0.02808286  max resid 0.03805044 
    ## Run 30 stress 8.953261e-05 
    ## ... Procrustes: rmse 0.04446425  max resid 0.06175419 
    ## Run 31 stress 0.0003978527 
    ## ... Procrustes: rmse 0.02919998  max resid 0.05926202 
    ## Run 32 stress 9.409612e-05 
    ## ... Procrustes: rmse 0.1386342  max resid 0.2040191 
    ## Run 33 stress 0.0004431821 
    ## ... Procrustes: rmse 0.02949068  max resid 0.06041345 
    ## Run 34 stress 0.0004238716 
    ## ... Procrustes: rmse 0.02936634  max resid 0.05993118 
    ## Run 35 stress 0.2124849 
    ## Run 36 stress 0.157905 
    ## Run 37 stress 0.2124849 
    ## Run 38 stress 9.546573e-05 
    ## ... Procrustes: rmse 0.0280834  max resid 0.03804076 
    ## Run 39 stress 9.471596e-05 
    ## ... Procrustes: rmse 0.2286749  max resid 0.3651781 
    ## Run 40 stress 9.31845e-05 
    ## ... Procrustes: rmse 0.2170181  max resid 0.342658 
    ## Run 41 stress 9.822229e-05 
    ## ... Procrustes: rmse 0.02617808  max resid 0.03610784 
    ## Run 42 stress 8.661886e-05 
    ## ... Procrustes: rmse 0.223055  max resid 0.3541714 
    ## Run 43 stress 9.674671e-05 
    ## ... Procrustes: rmse 0.1093968  max resid 0.157562 
    ## Run 44 stress 9.348484e-05 
    ## ... Procrustes: rmse 0.07220948  max resid 0.1016033 
    ## Run 45 stress 9.067523e-05 
    ## ... Procrustes: rmse 0.1236499  max resid 0.1799244 
    ## Run 46 stress 9.540868e-05 
    ## ... Procrustes: rmse 0.02799154  max resid 0.03806032 
    ## Run 47 stress 9.876226e-05 
    ## ... Procrustes: rmse 0.02798979  max resid 0.03804573 
    ## Run 48 stress 0.2601501 
    ## Run 49 stress 0.0001678864 
    ## ... Procrustes: rmse 0.0278884  max resid 0.05197959 
    ## Run 50 stress 9.001043e-05 
    ## ... Procrustes: rmse 0.02790146  max resid 0.03806159 
    ## Run 51 stress 0.2544152 
    ## Run 52 stress 8.669117e-05 
    ## ... Procrustes: rmse 0.02791127  max resid 0.03807957 
    ## Run 53 stress 9.64575e-05 
    ## ... Procrustes: rmse 0.1899495  max resid 0.2923403 
    ## Run 54 stress 0.0004164256 
    ## ... Procrustes: rmse 0.02931124  max resid 0.0597021 
    ## Run 55 stress 0.2601494 
    ## Run 56 stress 9.996105e-05 
    ## ... Procrustes: rmse 0.02788316  max resid 0.0380258 
    ## Run 57 stress 9.678533e-05 
    ## ... Procrustes: rmse 0.2330485  max resid 0.3737161 
    ## Run 58 stress 9.597225e-05 
    ## ... Procrustes: rmse 0.2136288  max resid 0.3362149 
    ## Run 59 stress 0.0002327584 
    ## ... Procrustes: rmse 0.02821393  max resid 0.05437118 
    ## Run 60 stress 9.431396e-05 
    ## ... Procrustes: rmse 0.164425  max resid 0.2473608 
    ## Run 61 stress 0.2253365 
    ## Run 62 stress 0.2264602 
    ## Run 63 stress 7.152482e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.004655905  max resid 0.006349094 
    ## ... Similar to previous best
    ## Run 64 stress 0.0004248714 
    ## ... Procrustes: rmse 0.02555784  max resid 0.05365979 
    ## Run 65 stress 0.0002585873 
    ## ... Procrustes: rmse 0.02424823  max resid 0.04893444 
    ## Run 66 stress 9.915129e-05 
    ## ... Procrustes: rmse 0.1549683  max resid 0.2303389 
    ## Run 67 stress 0.0004126465 
    ## ... Procrustes: rmse 0.02545997  max resid 0.05335099 
    ## Run 68 stress 9.926189e-05 
    ## ... Procrustes: rmse 0.02324262  max resid 0.03170506 
    ## Run 69 stress 0.293122 
    ## Run 70 stress 8.826567e-05 
    ## ... Procrustes: rmse 0.2254244  max resid 0.3576532 
    ## Run 71 stress 0.000409809 
    ## ... Procrustes: rmse 0.02537606  max resid 0.05308611 
    ## Run 72 stress 9.347653e-05 
    ## ... Procrustes: rmse 0.2314884  max resid 0.3695064 
    ## Run 73 stress 8.112339e-05 
    ## ... Procrustes: rmse 0.238589  max resid 0.3835621 
    ## Run 74 stress 9.485029e-05 
    ## ... Procrustes: rmse 0.006683532  max resid 0.009122971 
    ## Run 75 stress 8.627185e-05 
    ## ... Procrustes: rmse 0.02304808  max resid 0.03756918 
    ## Run 76 stress 9.529686e-05 
    ## ... Procrustes: rmse 0.02334814  max resid 0.03174007 
    ## Run 77 stress 9.987979e-05 
    ## ... Procrustes: rmse 0.1040634  max resid 0.1487729 
    ## Run 78 stress 9.439672e-05 
    ## ... Procrustes: rmse 0.02336367  max resid 0.03176656 
    ## Run 79 stress 8.733052e-05 
    ## ... Procrustes: rmse 0.02341127  max resid 0.03213694 
    ## Run 80 stress 8.788245e-05 
    ## ... Procrustes: rmse 0.02337432  max resid 0.03183563 
    ## Run 81 stress 8.198378e-05 
    ## ... Procrustes: rmse 0.02346316  max resid 0.03213966 
    ## Run 82 stress 9.652065e-05 
    ## ... Procrustes: rmse 0.02335089  max resid 0.03173348 
    ## Run 83 stress 0.0003687329 
    ## ... Procrustes: rmse 0.02510801  max resid 0.05219332 
    ## Run 84 stress 0.157905 
    ## Run 85 stress 0.26012 
    ## Run 86 stress 0.2544152 
    ## Run 87 stress 9.758613e-05 
    ## ... Procrustes: rmse 0.2124484  max resid 0.3328235 
    ## Run 88 stress 0.1648944 
    ## Run 89 stress 0.0003998286 
    ## ... Procrustes: rmse 0.02534213  max resid 0.0529724 
    ## Run 90 stress 9.228591e-05 
    ## ... Procrustes: rmse 0.02337876  max resid 0.03182466 
    ## Run 91 stress 9.562783e-05 
    ## ... Procrustes: rmse 0.02346606  max resid 0.0321439 
    ## Run 92 stress 9.228562e-05 
    ## ... Procrustes: rmse 0.2026907  max resid 0.3145714 
    ## Run 93 stress 8.976119e-05 
    ## ... Procrustes: rmse 0.02339433  max resid 0.03213694 
    ## Run 94 stress 9.827677e-05 
    ## ... Procrustes: rmse 0.2356278  max resid 0.3776909 
    ## Run 95 stress 9.733713e-05 
    ## ... Procrustes: rmse 0.1475493  max resid 0.2179223 
    ## Run 96 stress 9.911398e-05 
    ## ... Procrustes: rmse 0.1896841  max resid 0.2908583 
    ## Run 97 stress 9.617317e-05 
    ## ... Procrustes: rmse 0.04239924  max resid 0.05867171 
    ## Run 98 stress 9.97382e-05 
    ## ... Procrustes: rmse 0.2169354  max resid 0.3413443 
    ## Run 99 stress 8.640536e-05 
    ## ... Procrustes: rmse 0.02341064  max resid 0.03213985 
    ## Run 100 stress 0.2124849 
    ## Run 101 stress 8.309955e-05 
    ## ... Procrustes: rmse 0.2385374  max resid 0.3834281 
    ## Run 102 stress 9.013162e-05 
    ## ... Procrustes: rmse 0.1334258  max resid 0.1948667 
    ## Run 103 stress 9.126502e-05 
    ## ... Procrustes: rmse 0.02337902  max resid 0.03182278 
    ## Run 104 stress 8.648807e-05 
    ## ... Procrustes: rmse 0.1953038  max resid 0.3010221 
    ## Run 105 stress 0.2570537 
    ## Run 106 stress 9.139502e-05 
    ## ... Procrustes: rmse 0.2254078  max resid 0.3576194 
    ## Run 107 stress 0.2933059 
    ## Run 108 stress 9.824595e-05 
    ## ... Procrustes: rmse 0.1433797  max resid 0.2109697 
    ## Run 109 stress 0.293306 
    ## Run 110 stress 0.2253366 
    ## Run 111 stress 7.547718e-05 
    ## ... Procrustes: rmse 0.2371788  max resid 0.380752 
    ## Run 112 stress 9.670465e-05 
    ## ... Procrustes: rmse 0.05381662  max resid 0.07490077 
    ## Run 113 stress 0.2570537 
    ## Run 114 stress 0.157905 
    ## Run 115 stress 8.363647e-05 
    ## ... Procrustes: rmse 0.2025644  max resid 0.3144161 
    ## Run 116 stress 9.763775e-05 
    ## ... Procrustes: rmse 0.01101962  max resid 0.01498761 
    ## Run 117 stress 9.674e-05 
    ## ... Procrustes: rmse 0.1523484  max resid 0.2258503 
    ## Run 118 stress 9.303626e-05 
    ## ... Procrustes: rmse 0.06611836  max resid 0.09252965 
    ## Run 119 stress 9.560465e-05 
    ## ... Procrustes: rmse 0.02324965  max resid 0.03172172 
    ## Run 120 stress 9.80692e-05 
    ## ... Procrustes: rmse 0.2001234  max resid 0.3097569 
    ## Run 121 stress 9.48748e-05 
    ## ... Procrustes: rmse 0.02334769  max resid 0.03174163 
    ## Run 122 stress 9.533635e-05 
    ## ... Procrustes: rmse 0.1948327  max resid 0.3001444 
    ## Run 123 stress 9.000967e-05 
    ## ... Procrustes: rmse 0.02337398  max resid 0.03178575 
    ## Run 124 stress 9.61198e-05 
    ## ... Procrustes: rmse 0.1736208  max resid 0.2623163 
    ## Run 125 stress 9.680494e-05 
    ## ... Procrustes: rmse 0.02344097  max resid 0.03171282 
    ## Run 126 stress 8.205698e-05 
    ## ... Procrustes: rmse 0.1679586  max resid 0.2524388 
    ## Run 127 stress 9.799368e-05 
    ## ... Procrustes: rmse 0.2014518  max resid 0.312305 
    ## Run 128 stress 9.250727e-05 
    ## ... Procrustes: rmse 0.2114068  max resid 0.330865 
    ## Run 129 stress 9.610447e-05 
    ## ... Procrustes: rmse 0.02335651  max resid 0.03174688 
    ## Run 130 stress 8.962987e-05 
    ## ... Procrustes: rmse 0.1519653  max resid 0.225362 
    ## Run 131 stress 9.966809e-05 
    ## ... Procrustes: rmse 0.01563652  max resid 0.02152708 
    ## Run 132 stress 0.2569849 
    ## Run 133 stress 0.2264602 
    ## Run 134 stress 0.0004016764 
    ## ... Procrustes: rmse 0.02537106  max resid 0.05306613 
    ## Run 135 stress 0.157905 
    ## Run 136 stress 9.444775e-05 
    ## ... Procrustes: rmse 0.1204292  max resid 0.1740689 
    ## Run 137 stress 9.688035e-05 
    ## ... Procrustes: rmse 0.02324679  max resid 0.03171493 
    ## Run 138 stress 9.770261e-05 
    ## ... Procrustes: rmse 0.0214023  max resid 0.02918103 
    ## Run 139 stress 9.804322e-05 
    ## ... Procrustes: rmse 0.2244904  max resid 0.3558585 
    ## Run 140 stress 9.175975e-05 
    ## ... Procrustes: rmse 0.211611  max resid 0.3312731 
    ## Run 141 stress 9.386079e-05 
    ## ... Procrustes: rmse 0.2198061  max resid 0.3468414 
    ## Run 142 stress 9.858377e-05 
    ## ... Procrustes: rmse 0.02339081  max resid 0.03258632 
    ## Run 143 stress 0.157905 
    ## Run 144 stress 9.499607e-05 
    ## ... Procrustes: rmse 0.2130158  max resid 0.3339557 
    ## Run 145 stress 0.1648944 
    ## Run 146 stress 9.640754e-05 
    ## ... Procrustes: rmse 0.234917  max resid 0.376303 
    ## Run 147 stress 9.785973e-05 
    ## ... Procrustes: rmse 0.02346594  max resid 0.03214329 
    ## Run 148 stress 9.141778e-05 
    ## ... Procrustes: rmse 0.02304068  max resid 0.03680602 
    ## Run 149 stress 9.579933e-05 
    ## ... Procrustes: rmse 0.1880181  max resid 0.2878092 
    ## Run 150 stress 0.157905 
    ## Run 151 stress 8.26394e-05 
    ## ... Procrustes: rmse 0.1662415  max resid 0.2495795 
    ## Run 152 stress 9.962196e-05 
    ## ... Procrustes: rmse 0.02334612  max resid 0.03172736 
    ## Run 153 stress 9.865435e-05 
    ## ... Procrustes: rmse 0.02324785  max resid 0.03171737 
    ## Run 154 stress 0.157905 
    ## Run 155 stress 8.672754e-05 
    ## ... Procrustes: rmse 0.2102988  max resid 0.3288 
    ## Run 156 stress 0.2264602 
    ## Run 157 stress 0.2569849 
    ## Run 158 stress 9.388378e-05 
    ## ... Procrustes: rmse 0.1851028  max resid 0.2824999 
    ## Run 159 stress 9.305083e-05 
    ## ... Procrustes: rmse 0.2243959  max resid 0.3557691 
    ## Run 160 stress 0.0004027276 
    ## ... Procrustes: rmse 0.02537895  max resid 0.05309279 
    ## Run 161 stress 0.0004160674 
    ## ... Procrustes: rmse 0.0254709  max resid 0.05338563 
    ## Run 162 stress 0.293133 
    ## Run 163 stress 0.2253371 
    ## Run 164 stress 9.388822e-05 
    ## ... Procrustes: rmse 0.1281045  max resid 0.186298 
    ## Run 165 stress 0.2601498 
    ## Run 166 stress 0.2253365 
    ## Run 167 stress 0.0004044925 
    ## ... Procrustes: rmse 0.02539386  max resid 0.05314005 
    ## Run 168 stress 8.957035e-05 
    ## ... Procrustes: rmse 0.1627036  max resid 0.2435311 
    ## Run 169 stress 9.439269e-05 
    ## ... Procrustes: rmse 0.2316629  max resid 0.3697809 
    ## Run 170 stress 0.2264602 
    ## Run 171 stress 0.293122 
    ## Run 172 stress 0.2601669 
    ## Run 173 stress 9.576187e-05 
    ## ... Procrustes: rmse 0.02344062  max resid 0.03171871 
    ## Run 174 stress 9.734351e-05 
    ## ... Procrustes: rmse 0.1388917  max resid 0.2037331 
    ## Run 175 stress 9.186964e-05 
    ## ... Procrustes: rmse 0.09911591  max resid 0.1412549 
    ## Run 176 stress 0.2601482 
    ## Run 177 stress 9.373567e-05 
    ## ... Procrustes: rmse 0.02325337  max resid 0.03172656 
    ## Run 178 stress 0.2570537 
    ## Run 179 stress 9.952541e-05 
    ## ... Procrustes: rmse 0.02334397  max resid 0.03172161 
    ## Run 180 stress 0.2908327 
    ## Run 181 stress 0.157905 
    ## Run 182 stress 0.2927618 
    ## Run 183 stress 9.050509e-05 
    ## ... Procrustes: rmse 0.02325172  max resid 0.03388662 
    ## Run 184 stress 0.157905 
    ## Run 185 stress 9.219574e-05 
    ## ... Procrustes: rmse 0.2135079  max resid 0.3348429 
    ## Run 186 stress 0.2570537 
    ## Run 187 stress 8.977385e-05 
    ## ... Procrustes: rmse 0.02338073  max resid 0.03183204 
    ## Run 188 stress 0.2124849 
    ## Run 189 stress 9.787983e-05 
    ## ... Procrustes: rmse 0.1934248  max resid 0.2975783 
    ## Run 190 stress 8.485347e-05 
    ## ... Procrustes: rmse 0.02336437  max resid 0.03178233 
    ## Run 191 stress 0.0003850926 
    ## ... Procrustes: rmse 0.02521308  max resid 0.05254928 
    ## Run 192 stress 0.293122 
    ## Run 193 stress 0.2290617 
    ## Run 194 stress 0.260162 
    ## Run 195 stress 9.563314e-05 
    ## ... Procrustes: rmse 0.02335256  max resid 0.0317422 
    ## Run 196 stress 9.973931e-05 
    ## ... Procrustes: rmse 0.2358024  max resid 0.3780075 
    ## Run 197 stress 9.951611e-05 
    ## ... Procrustes: rmse 0.008519689  max resid 0.01156092 
    ## Run 198 stress 8.704764e-05 
    ## ... Procrustes: rmse 0.1923499  max resid 0.295602 
    ## Run 199 stress 9.087521e-05 
    ## ... Procrustes: rmse 0.02319145  max resid 0.03551219 
    ## Run 200 stress 9.526182e-05 
    ## ... Procrustes: rmse 0.1777507  max resid 0.2695417 
    ## Run 201 stress 0.0004067887 
    ## ... Procrustes: rmse 0.02540413  max resid 0.05317247 
    ## Run 202 stress 0.0003201968 
    ## ... Procrustes: rmse 0.02472098  max resid 0.0508208 
    ## Run 203 stress 0.157905 
    ## Run 204 stress 9.877295e-05 
    ## ... Procrustes: rmse 0.2018388  max resid 0.3129112 
    ## Run 205 stress 9.351531e-05 
    ## ... Procrustes: rmse 0.0233918  max resid 0.03268512 
    ## Run 206 stress 0.2124849 
    ## Run 207 stress 9.70796e-05 
    ## ... Procrustes: rmse 0.01709782  max resid 0.02354722 
    ## Run 208 stress 9.444313e-05 
    ## ... Procrustes: rmse 0.02336068  max resid 0.03175264 
    ## Run 209 stress 9.648175e-05 
    ## ... Procrustes: rmse 0.02324518  max resid 0.03171106 
    ## Run 210 stress 9.321676e-05 
    ## ... Procrustes: rmse 0.05236381  max resid 0.0727519 
    ## Run 211 stress 9.28453e-05 
    ## ... Procrustes: rmse 0.02302485  max resid 0.03732513 
    ## Run 212 stress 8.899883e-05 
    ## ... Procrustes: rmse 0.1291078  max resid 0.1879883 
    ## Run 213 stress 8.963272e-05 
    ## ... Procrustes: rmse 0.2343206  max resid 0.3750948 
    ## Run 214 stress 0.2124849 
    ## Run 215 stress 9.515017e-05 
    ## ... Procrustes: rmse 0.2328149  max resid 0.3721023 
    ## Run 216 stress 9.396253e-05 
    ## ... Procrustes: rmse 0.2019991  max resid 0.3132937 
    ## Run 217 stress 8.520633e-05 
    ## ... Procrustes: rmse 0.009394276  max resid 0.01278159 
    ## Run 218 stress 9.333147e-05 
    ## ... Procrustes: rmse 0.02346131  max resid 0.032134 
    ## Run 219 stress 0.0003452355 
    ## ... Procrustes: rmse 0.02491968  max resid 0.05154217 
    ## Run 220 stress 9.501539e-05 
    ## ... Procrustes: rmse 0.02334999  max resid 0.0317396 
    ## Run 221 stress 0.0004075227 
    ## ... Procrustes: rmse 0.02541789  max resid 0.05321771 
    ## Run 222 stress 0.0004088343 
    ## ... Procrustes: rmse 0.02542878  max resid 0.05325205 
    ## Run 223 stress 9.61234e-05 
    ## ... Procrustes: rmse 0.02339094  max resid 0.03213654 
    ## Run 224 stress 0.2124849 
    ## Run 225 stress 8.963831e-05 
    ## ... Procrustes: rmse 0.2114678  max resid 0.3310527 
    ## Run 226 stress 9.699822e-05 
    ## ... Procrustes: rmse 0.2317561  max resid 0.3700317 
    ## Run 227 stress 9.551024e-05 
    ## ... Procrustes: rmse 0.02338946  max resid 0.03213716 
    ## Run 228 stress 8.987281e-05 
    ## ... Procrustes: rmse 0.02335827  max resid 0.03176057 
    ## Run 229 stress 8.091536e-05 
    ## ... Procrustes: rmse 0.1691633  max resid 0.2545516 
    ## Run 230 stress 0.2253365 
    ## Run 231 stress 9.731103e-05 
    ## ... Procrustes: rmse 0.2262702  max resid 0.3593014 
    ## Run 232 stress 9.557055e-05 
    ## ... Procrustes: rmse 0.02291383  max resid 0.03129224 
    ## Run 233 stress 8.949418e-05 
    ## ... Procrustes: rmse 0.2122524  max resid 0.3324794 
    ## Run 234 stress 8.741674e-05 
    ## ... Procrustes: rmse 0.2048509  max resid 0.3185636 
    ## Run 235 stress 9.683721e-05 
    ## ... Procrustes: rmse 0.2361783  max resid 0.3787533 
    ## Run 236 stress 8.881334e-05 
    ## ... Procrustes: rmse 0.1651849  max resid 0.2477321 
    ## Run 237 stress 8.806661e-05 
    ## ... Procrustes: rmse 0.2372664  max resid 0.3809315 
    ## Run 238 stress 9.409673e-05 
    ## ... Procrustes: rmse 0.2130703  max resid 0.3340025 
    ## Run 239 stress 0.157905 
    ## Run 240 stress 9.393929e-05 
    ## ... Procrustes: rmse 0.2033077  max resid 0.3156334 
    ## Run 241 stress 0.2124849 
    ## Run 242 stress 0.2927601 
    ## Run 243 stress 8.838134e-05 
    ## ... Procrustes: rmse 0.01884195  max resid 0.02580558 
    ## Run 244 stress 9.199699e-05 
    ## ... Procrustes: rmse 0.1978818  max resid 0.3057106 
    ## Run 245 stress 0.0004321979 
    ## ... Procrustes: rmse 0.02561757  max resid 0.05384532 
    ## Run 246 stress 8.073707e-05 
    ## ... Procrustes: rmse 0.1705401  max resid 0.2569995 
    ## Run 247 stress 9.150471e-05 
    ## ... Procrustes: rmse 0.08521029  max resid 0.1204255 
    ## Run 248 stress 8.69925e-05 
    ## ... Procrustes: rmse 0.2086372  max resid 0.3257307 
    ## Run 249 stress 9.631683e-05 
    ## ... Procrustes: rmse 0.02344179  max resid 0.03171409 
    ## Run 250 stress 9.678931e-05 
    ## ... Procrustes: rmse 0.2214078  max resid 0.3500097 
    ## Run 251 stress 8.007998e-05 
    ## ... Procrustes: rmse 0.2385913  max resid 0.3835236 
    ## Run 252 stress 9.617721e-05 
    ## ... Procrustes: rmse 0.02344267  max resid 0.03171563 
    ## Run 253 stress 9.482805e-05 
    ## ... Procrustes: rmse 0.02334942  max resid 0.03174064 
    ## Run 254 stress 7.805445e-05 
    ## ... Procrustes: rmse 0.229548  max resid 0.3656349 
    ## Run 255 stress 0.0001539111 
    ## ... Procrustes: rmse 0.02351059  max resid 0.04509878 
    ## Run 256 stress 9.030431e-05 
    ## ... Procrustes: rmse 0.1282844  max resid 0.1866191 
    ## Run 257 stress 7.893109e-05 
    ## ... Procrustes: rmse 0.2385999  max resid 0.3835403 
    ## Run 258 stress 9.761747e-05 
    ## ... Procrustes: rmse 0.2341443  max resid 0.3747511 
    ## Run 259 stress 0.157905 
    ## Run 260 stress 9.004246e-05 
    ## ... Procrustes: rmse 0.02336811  max resid 0.03181988 
    ## Run 261 stress 9.516614e-05 
    ## ... Procrustes: rmse 0.1377176  max resid 0.2018615 
    ## Run 262 stress 0.2124849 
    ## Run 263 stress 9.864927e-05 
    ## ... Procrustes: rmse 0.2325085  max resid 0.371518 
    ## Run 264 stress 8.579411e-05 
    ## ... Procrustes: rmse 0.1987632  max resid 0.3073508 
    ## Run 265 stress 9.08742e-05 
    ## ... Procrustes: rmse 0.1999151  max resid 0.3094391 
    ## Run 266 stress 9.827027e-05 
    ## ... Procrustes: rmse 0.1548149  max resid 0.2300774 
    ## Run 267 stress 9.539375e-05 
    ## ... Procrustes: rmse 0.02344137  max resid 0.03171858 
    ## Run 268 stress 9.720488e-05 
    ## ... Procrustes: rmse 0.2108004  max resid 0.3297151 
    ## Run 269 stress 9.388345e-05 
    ## ... Procrustes: rmse 0.02342946  max resid 0.03212067 
    ## Run 270 stress 9.809444e-05 
    ## ... Procrustes: rmse 0.2332515  max resid 0.3729814 
    ## Run 271 stress 9.656065e-05 
    ## ... Procrustes: rmse 0.2334765  max resid 0.3734286 
    ## Run 272 stress 8.884279e-05 
    ## ... Procrustes: rmse 0.02346689  max resid 0.03213014 
    ## Run 273 stress 0.2568897 
    ## Run 274 stress 9.995949e-05 
    ## ... Procrustes: rmse 0.02336807  max resid 0.03178727 
    ## Run 275 stress 9.770969e-05 
    ## ... Procrustes: rmse 0.07584605  max resid 0.1066582 
    ## Run 276 stress 9.760772e-05 
    ## ... Procrustes: rmse 0.02325648  max resid 0.03173493 
    ## Run 277 stress 9.607702e-05 
    ## ... Procrustes: rmse 0.02335916  max resid 0.03175283 
    ## Run 278 stress 9.21006e-05 
    ## ... Procrustes: rmse 0.2376931  max resid 0.381755 
    ## Run 279 stress 9.606388e-05 
    ## ... Procrustes: rmse 0.02336228  max resid 0.03175114 
    ## Run 280 stress 9.36615e-05 
    ## ... Procrustes: rmse 0.2330455  max resid 0.3725054 
    ## Run 281 stress 9.655096e-05 
    ## ... Procrustes: rmse 0.1978487  max resid 0.3056409 
    ## Run 282 stress 9.925674e-05 
    ## ... Procrustes: rmse 0.1110304  max resid 0.1595119 
    ## Run 283 stress 8.589419e-05 
    ## ... Procrustes: rmse 0.02346549  max resid 0.03214077 
    ## Run 284 stress 9.007918e-05 
    ## ... Procrustes: rmse 0.2120097  max resid 0.3320154 
    ## Run 285 stress 0.0001042112 
    ## ... Procrustes: rmse 0.02322188  max resid 0.04279864 
    ## Run 286 stress 8.817565e-05 
    ## ... Procrustes: rmse 0.02346234  max resid 0.03212959 
    ## Run 287 stress 9.567808e-05 
    ## ... Procrustes: rmse 0.1877221  max resid 0.2871732 
    ## Run 288 stress 0.2601497 
    ## Run 289 stress 8.923371e-05 
    ## ... Procrustes: rmse 0.2385747  max resid 0.3834956 
    ## Run 290 stress 9.147133e-05 
    ## ... Procrustes: rmse 0.07132423  max resid 0.1000017 
    ## Run 291 stress 9.660381e-05 
    ## ... Procrustes: rmse 0.1339753  max resid 0.1956612 
    ## Run 292 stress 0.2570537 
    ## Run 293 stress 9.473481e-05 
    ## ... Procrustes: rmse 0.01104296  max resid 0.01501644 
    ## Run 294 stress 9.41788e-05 
    ## ... Procrustes: rmse 0.2150016  max resid 0.3377506 
    ## Run 295 stress 0.0004115356 
    ## ... Procrustes: rmse 0.02544753  max resid 0.05331166 
    ## Run 296 stress 0.2570537 
    ## Run 297 stress 7.93177e-05 
    ## ... Procrustes: rmse 0.1469016  max resid 0.2169363 
    ## Run 298 stress 6.441672e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.02341581  max resid 0.0319603 
    ## Run 299 stress 9.820892e-05 
    ## ... Procrustes: rmse 0.2358174  max resid 0.3721375 
    ## Run 300 stress 9.203726e-05 
    ## ... Procrustes: rmse 0.1124657  max resid 0.1593656 
    ## Run 301 stress 9.562852e-05 
    ## ... Procrustes: rmse 0.0002225718  max resid 0.0003600281 
    ## ... Similar to previous best
    ## Run 302 stress 9.760097e-05 
    ## ... Procrustes: rmse 0.2175171  max resid 0.3366007 
    ## Run 303 stress 9.82656e-05 
    ## ... Procrustes: rmse 0.2528753  max resid 0.4063775 
    ## Run 304 stress 9.206203e-05 
    ## ... Procrustes: rmse 0.1809899  max resid 0.2705484 
    ## Run 305 stress 8.9575e-05 
    ## ... Procrustes: rmse 0.07538735  max resid 0.1047278 
    ## Run 306 stress 0.2570537 
    ## Run 307 stress 9.664606e-05 
    ## ... Procrustes: rmse 0.0001854758  max resid 0.000299472 
    ## ... Similar to previous best
    ## Run 308 stress 9.402633e-05 
    ## ... Procrustes: rmse 0.250522  max resid 0.4015592 
    ## Run 309 stress 9.530971e-05 
    ## ... Procrustes: rmse 0.2456074  max resid 0.3915812 
    ## Run 310 stress 9.327114e-05 
    ## ... Procrustes: rmse 0.0002172499  max resid 0.0003527163 
    ## ... Similar to previous best
    ## Run 311 stress 0.2601489 
    ## Run 312 stress 8.703109e-05 
    ## ... Procrustes: rmse 0.208646  max resid 0.3201278 
    ## Run 313 stress 9.135339e-05 
    ## ... Procrustes: rmse 0.162896  max resid 0.2396721 
    ## Run 314 stress 9.981646e-05 
    ## ... Procrustes: rmse 0.2522891  max resid 0.4052044 
    ## Run 315 stress 8.426994e-05 
    ## ... Procrustes: rmse 0.0001769642  max resid 0.0002537157 
    ## ... Similar to previous best
    ## Run 316 stress 0.2569849 
    ## Run 317 stress 0.2124849 
    ## Run 318 stress 9.83313e-05 
    ## ... Procrustes: rmse 0.1567734  max resid 0.2293544 
    ## Run 319 stress 9.588779e-05 
    ## ... Procrustes: rmse 0.2232309  max resid 0.347594 
    ## Run 320 stress 8.700249e-05 
    ## ... Procrustes: rmse 0.2124345  max resid 0.3272186 
    ## Run 321 stress 4.043551e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2527535  max resid 0.4061352 
    ## Run 322 stress 9.634462e-05 
    ## ... Procrustes: rmse 0.02346103  max resid 0.03710865 
    ## Run 323 stress 6.829968e-05 
    ## ... Procrustes: rmse 0.2519538  max resid 0.3443342 
    ## Run 324 stress 0.293122 
    ## Run 325 stress 9.379886e-05 
    ## ... Procrustes: rmse 0.2527105  max resid 0.3452181 
    ## Run 326 stress 9.807786e-05 
    ## ... Procrustes: rmse 0.2526401  max resid 0.3451974 
    ## Run 327 stress 9.577035e-05 
    ## ... Procrustes: rmse 0.06520242  max resid 0.1000071 
    ## Run 328 stress 0.2601501 
    ## Run 329 stress 9.131831e-05 
    ## ... Procrustes: rmse 0.1652083  max resid 0.2348124 
    ## Run 330 stress 0.2546681 
    ## Run 331 stress 9.64244e-05 
    ## ... Procrustes: rmse 0.2527672  max resid 0.3454477 
    ## Run 332 stress 0.0003938074 
    ## ... Procrustes: rmse 0.2512273  max resid 0.3575332 
    ## Run 333 stress 9.645317e-05 
    ## ... Procrustes: rmse 0.2222465  max resid 0.3055838 
    ## Run 334 stress 9.973833e-05 
    ## ... Procrustes: rmse 0.05572797  max resid 0.08607424 
    ## Run 335 stress 9.879633e-05 
    ## ... Procrustes: rmse 0.2374079  max resid 0.3248587 
    ## Run 336 stress 9.15042e-05 
    ## ... Procrustes: rmse 0.02540757  max resid 0.04011193 
    ## Run 337 stress 8.964841e-05 
    ## ... Procrustes: rmse 0.02369217  max resid 0.03745025 
    ## Run 338 stress 8.959784e-05 
    ## ... Procrustes: rmse 0.002788942  max resid 0.004498695 
    ## ... Similar to previous best
    ## Run 339 stress 9.138661e-05 
    ## ... Procrustes: rmse 0.04828358  max resid 0.07500511 
    ## Run 340 stress 9.79981e-05 
    ## ... Procrustes: rmse 0.2526405  max resid 0.3451969 
    ## Run 341 stress 9.496116e-05 
    ## ... Procrustes: rmse 0.2526472  max resid 0.3452039 
    ## Run 342 stress 8.929349e-05 
    ## ... Procrustes: rmse 0.04745132  max resid 0.07383169 
    ## Run 343 stress 9.345179e-05 
    ## ... Procrustes: rmse 0.2527098  max resid 0.3452171 
    ## Run 344 stress 0.2124849 
    ## Run 345 stress 8.890562e-05 
    ## ... Procrustes: rmse 0.006429631  max resid 0.01034128 
    ## Run 346 stress 9.318892e-05 
    ## ... Procrustes: rmse 0.2097101  max resid 0.289964 
    ## Run 347 stress 9.427639e-05 
    ## ... Procrustes: rmse 0.03911852  max resid 0.06121642 
    ## Run 348 stress 9.849265e-05 
    ## ... Procrustes: rmse 0.2527091  max resid 0.345209 
    ## Run 349 stress 9.305293e-05 
    ## ... Procrustes: rmse 0.004739471  max resid 0.007582899 
    ## ... Similar to previous best
    ## Run 350 stress 0.2933059 
    ## Run 351 stress 9.75702e-05 
    ## ... Procrustes: rmse 0.005829983  max resid 0.0093144 
    ## Run 352 stress 9.035742e-05 
    ## ... Procrustes: rmse 0.05068979  max resid 0.07862979 
    ## Run 353 stress 0.0004097294 
    ## ... Procrustes: rmse 0.2512048  max resid 0.3577378 
    ## Run 354 stress 9.418272e-05 
    ## ... Procrustes: rmse 0.25276  max resid 0.345213 
    ## Run 355 stress 8.296954e-05 
    ## ... Procrustes: rmse 0.2527212  max resid 0.3452422 
    ## Run 356 stress 0.157905 
    ## Run 357 stress 0.157905 
    ## Run 358 stress 0.00042663 
    ## ... Procrustes: rmse 0.2511729  max resid 0.3580168 
    ## Run 359 stress 0.000406353 
    ## ... Procrustes: rmse 0.251206  max resid 0.3577207 
    ## Run 360 stress 0.0003879623 
    ## ... Procrustes: rmse 0.2512383  max resid 0.357445 
    ## Run 361 stress 8.777199e-05 
    ## ... Procrustes: rmse 0.2527395  max resid 0.3452559 
    ## Run 362 stress 9.281286e-05 
    ## ... Procrustes: rmse 0.1440053  max resid 0.2078293 
    ## Run 363 stress 9.285153e-05 
    ## ... Procrustes: rmse 0.2527096  max resid 0.3452632 
    ## Run 364 stress 0.0002499986 
    ## ... Procrustes: rmse 0.2515114  max resid 0.3551022 
    ## Run 365 stress 9.192027e-05 
    ## ... Procrustes: rmse 0.0479762  max resid 0.07456368 
    ## Run 366 stress 9.888674e-05 
    ## ... Procrustes: rmse 0.08147705  max resid 0.1234337 
    ## Run 367 stress 9.408358e-05 
    ## ... Procrustes: rmse 0.2527108  max resid 0.3452141 
    ## Run 368 stress 0.293122 
    ## Run 369 stress 0.2601496 
    ## Run 370 stress 0.293122 
    ## Run 371 stress 9.662485e-05 
    ## ... Procrustes: rmse 0.09096751  max resid 0.1368943 
    ## Run 372 stress 9.770555e-05 
    ## ... Procrustes: rmse 0.006167315  max resid 0.009837928 
    ## Run 373 stress 9.469413e-05 
    ## ... Procrustes: rmse 0.02017404  max resid 0.03198833 
    ## Run 374 stress 9.250519e-05 
    ## ... Procrustes: rmse 0.1306832  max resid 0.1905434 
    ## Run 375 stress 9.254102e-05 
    ## ... Procrustes: rmse 0.04419149  max resid 0.06884675 
    ## Run 376 stress 0.293306 
    ## Run 377 stress 9.070449e-05 
    ## ... Procrustes: rmse 0.2527229  max resid 0.3452621 
    ## Run 378 stress 0.157905 
    ## Run 379 stress 9.316598e-05 
    ## ... Procrustes: rmse 0.2435212  max resid 0.3328821 
    ## Run 380 stress 9.13339e-05 
    ## ... Procrustes: rmse 0.2484719  max resid 0.3394678 
    ## Run 381 stress 6.847673e-05 
    ## ... Procrustes: rmse 0.00502974  max resid 0.008055649 
    ## Run 382 stress 0.2570537 
    ## Run 383 stress 9.596155e-05 
    ## ... Procrustes: rmse 0.2310053  max resid 0.3166528 
    ## Run 384 stress 9.594387e-05 
    ## ... Procrustes: rmse 0.2527058  max resid 0.3452638 
    ## Run 385 stress 8.65988e-05 
    ## ... Procrustes: rmse 0.2413449  max resid 0.3300235 
    ## Run 386 stress 9.003086e-05 
    ## ... Procrustes: rmse 0.02025407  max resid 0.03213608 
    ## Run 387 stress 8.979317e-05 
    ## ... Procrustes: rmse 0.05262298  max resid 0.0815605 
    ## Run 388 stress 9.486866e-05 
    ## ... Procrustes: rmse 0.2526463  max resid 0.345203 
    ## Run 389 stress 0.0003383392 
    ## ... Procrustes: rmse 0.2513499  max resid 0.3565379 
    ## Run 390 stress 0.2124849 
    ## Run 391 stress 8.204783e-05 
    ## ... Procrustes: rmse 0.2527568  max resid 0.3454262 
    ## Run 392 stress 6.945103e-05 
    ## ... Procrustes: rmse 0.1295625  max resid 0.1892127 
    ## Run 393 stress 9.763361e-05 
    ## ... Procrustes: rmse 0.008609705  max resid 0.01375849 
    ## Run 394 stress 8.120672e-05 
    ## ... Procrustes: rmse 0.09532299  max resid 0.1428777 
    ## Run 395 stress 0.2264602 
    ## Run 396 stress 0.2253365 
    ## Run 397 stress 9.315396e-05 
    ## ... Procrustes: rmse 0.1797908  max resid 0.2529907 
    ## Run 398 stress 8.362246e-05 
    ## ... Procrustes: rmse 0.252722  max resid 0.3452391 
    ## Run 399 stress 9.269293e-05 
    ## ... Procrustes: rmse 0.02523932  max resid 0.03992786 
    ## Run 400 stress 7.996183e-05 
    ## ... Procrustes: rmse 0.2527327  max resid 0.3452828 
    ## Run 401 stress 8.836013e-05 
    ## ... Procrustes: rmse 0.06152871  max resid 0.09465845 
    ## Run 402 stress 9.88227e-05 
    ## ... Procrustes: rmse 0.2526396  max resid 0.3451936 
    ## Run 403 stress 9.101986e-05 
    ## ... Procrustes: rmse 0.01520499  max resid 0.02416112 
    ## Run 404 stress 8.882408e-05 
    ## ... Procrustes: rmse 0.1052757  max resid 0.1566035 
    ## Run 405 stress 8.91172e-05 
    ## ... Procrustes: rmse 0.1785024  max resid 0.2513932 
    ## Run 406 stress 8.730903e-05 
    ## ... Procrustes: rmse 0.1145559  max resid 0.1691827 
    ## Run 407 stress 0.2601487 
    ## Run 408 stress 8.855642e-05 
    ## ... Procrustes: rmse 0.03902471  max resid 0.06110079 
    ## Run 409 stress 9.829788e-05 
    ## ... Procrustes: rmse 0.2527682  max resid 0.3454544 
    ## Run 410 stress 9.72717e-05 
    ## ... Procrustes: rmse 0.2527047  max resid 0.345211 
    ## Run 411 stress 9.10058e-05 
    ## ... Procrustes: rmse 0.1176578  max resid 0.1733277 
    ## Run 412 stress 9.608028e-05 
    ## ... Procrustes: rmse 0.1328004  max resid 0.1933123 
    ## Run 413 stress 0.2570419 
    ## Run 414 stress 9.413412e-05 
    ## ... Procrustes: rmse 0.05403867  max resid 0.08360181 
    ## Run 415 stress 0.0004133358 
    ## ... Procrustes: rmse 0.2512112  max resid 0.3576854 
    ## Run 416 stress 0.2570537 
    ## Run 417 stress 9.351766e-05 
    ## ... Procrustes: rmse 0.252709  max resid 0.3452593 
    ## Run 418 stress 8.999259e-05 
    ## ... Procrustes: rmse 0.1352623  max resid 0.1965768 
    ## Run 419 stress 0.2546681 
    ## Run 420 stress 9.109698e-05 
    ## ... Procrustes: rmse 0.1025491  max resid 0.1528786 
    ## Run 421 stress 0.0001603994 
    ## ... Procrustes: rmse 0.2517458  max resid 0.3531943 
    ## Run 422 stress 9.637884e-05 
    ## ... Procrustes: rmse 0.01344762  max resid 0.02140566 
    ## Run 423 stress 9.390159e-05 
    ## ... Procrustes: rmse 0.0007249493  max resid 0.001132753 
    ## ... Similar to previous best
    ## Run 424 stress 8.997977e-05 
    ## ... Procrustes: rmse 0.1828396  max resid 0.2567707 
    ## Run 425 stress 0.2253371 
    ## Run 426 stress 9.930516e-05 
    ## ... Procrustes: rmse 0.0297975  max resid 0.04697896 
    ## Run 427 stress 0.2264602 
    ## Run 428 stress 7.677055e-05 
    ## ... Procrustes: rmse 0.0001216646  max resid 0.0002261099 
    ## ... Similar to previous best
    ## Run 429 stress 9.931136e-05 
    ## ... Procrustes: rmse 0.2527028  max resid 0.3452624 
    ## Run 430 stress 9.498609e-05 
    ## ... Procrustes: rmse 0.2527089  max resid 0.3452143 
    ## Run 431 stress 9.488027e-05 
    ## ... Procrustes: rmse 0.0376939  max resid 0.05904017 
    ## Run 432 stress 9.597313e-05 
    ## ... Procrustes: rmse 0.006357523  max resid 0.01015535 
    ## Run 433 stress 8.830503e-05 
    ## ... Procrustes: rmse 0.2344276  max resid 0.3210406 
    ## Run 434 stress 9.731414e-05 
    ## ... Procrustes: rmse 0.1482996  max resid 0.213342 
    ## Run 435 stress 0.0004226261 
    ## ... Procrustes: rmse 0.2511802  max resid 0.3579572 
    ## Run 436 stress 8.72618e-05 
    ## ... Procrustes: rmse 0.009221818  max resid 0.01478289 
    ## Run 437 stress 9.159032e-05 
    ## ... Procrustes: rmse 0.008633094  max resid 0.01377593 
    ## Run 438 stress 0.2917768 
    ## Run 439 stress 9.793003e-05 
    ## ... Procrustes: rmse 0.252642  max resid 0.3451953 
    ## Run 440 stress 9.939016e-05 
    ## ... Procrustes: rmse 0.2527668  max resid 0.3454397 
    ## Run 441 stress 9.315931e-05 
    ## ... Procrustes: rmse 0.05476574  max resid 0.0846797 
    ## Run 442 stress 9.477582e-05 
    ## ... Procrustes: rmse 0.04059398  max resid 0.06340915 
    ## Run 443 stress 9.762102e-05 
    ## ... Procrustes: rmse 0.23877  max resid 0.3267163 
    ## Run 444 stress 9.776418e-05 
    ## ... Procrustes: rmse 0.2526464  max resid 0.3451988 
    ## Run 445 stress 0.0004062375 
    ## ... Procrustes: rmse 0.2512114  max resid 0.3577149 
    ## Run 446 stress 9.870901e-05 
    ## ... Procrustes: rmse 0.03491064  max resid 0.0548297 
    ## Run 447 stress 9.726481e-05 
    ## ... Procrustes: rmse 0.2527052  max resid 0.3452109 
    ## Run 448 stress 9.634115e-05 
    ## ... Procrustes: rmse 0.2526441  max resid 0.3452046 
    ## Run 449 stress 9.290386e-05 
    ## ... Procrustes: rmse 0.1066132  max resid 0.1584251 
    ## Run 450 stress 9.435265e-05 
    ## ... Procrustes: rmse 0.01025234  max resid 0.0163478 
    ## Run 451 stress 0.0004239304 
    ## ... Procrustes: rmse 0.2511773  max resid 0.3579766 
    ## Run 452 stress 9.393369e-05 
    ## ... Procrustes: rmse 0.02131544  max resid 0.03374859 
    ## Run 453 stress 9.591116e-05 
    ## ... Procrustes: rmse 0.2527591  max resid 0.3452003 
    ## Run 454 stress 0.157905 
    ## Run 455 stress 9.892011e-05 
    ## ... Procrustes: rmse 0.1675262  max resid 0.237723 
    ## Run 456 stress 8.593533e-05 
    ## ... Procrustes: rmse 0.05783104  max resid 0.08922316 
    ## Run 457 stress 8.539301e-05 
    ## ... Procrustes: rmse 0.03733429  max resid 0.05849753 
    ## Run 458 stress 9.434862e-05 
    ## ... Procrustes: rmse 0.2527271  max resid 0.3452508 
    ## Run 459 stress 0.2601482 
    ## Run 460 stress 8.521725e-05 
    ## ... Procrustes: rmse 0.2526779  max resid 0.345238 
    ## Run 461 stress 0.0004030517 
    ## ... Procrustes: rmse 0.2512115  max resid 0.3576725 
    ## Run 462 stress 0.2253365 
    ## Run 463 stress 0.2264602 
    ## Run 464 stress 0.2546681 
    ## Run 465 stress 0.157905 
    ## Run 466 stress 9.788845e-05 
    ## ... Procrustes: rmse 0.02657366  max resid 0.04199738 
    ## Run 467 stress 0.2124849 
    ## Run 468 stress 0.000424667 
    ## ... Procrustes: rmse 0.2511758  max resid 0.3579916 
    ## Run 469 stress 9.726151e-05 
    ## ... Procrustes: rmse 0.2527587  max resid 0.3451974 
    ## Run 470 stress 9.802549e-05 
    ## ... Procrustes: rmse 0.2526403  max resid 0.3451949 
    ## Run 471 stress 0.0002133087 
    ## ... Procrustes: rmse 0.2516005  max resid 0.3543703 
    ## Run 472 stress 0.0004244925 
    ## ... Procrustes: rmse 0.2511775  max resid 0.3579901 
    ## Run 473 stress 0.2927619 
    ## Run 474 stress 8.269194e-05 
    ## ... Procrustes: rmse 0.0148055  max resid 0.0236054 
    ## Run 475 stress 7.566233e-05 
    ## ... Procrustes: rmse 0.07431426  max resid 0.1132333 
    ## Run 476 stress 9.476609e-05 
    ## ... Procrustes: rmse 0.0005408612  max resid 0.0008682862 
    ## ... Similar to previous best
    ## Run 477 stress 9.471501e-05 
    ## ... Procrustes: rmse 0.05564902  max resid 0.08600005 
    ## Run 478 stress 9.348232e-05 
    ## ... Procrustes: rmse 0.1491345  max resid 0.214469 
    ## Run 479 stress 0.0002039265 
    ## ... Procrustes: rmse 0.251624  max resid 0.3541749 
    ## Run 480 stress 8.845079e-05 
    ## ... Procrustes: rmse 0.002270576  max resid 0.003634676 
    ## ... Similar to previous best
    ## Run 481 stress 9.228658e-05 
    ## ... Procrustes: rmse 0.003144695  max resid 0.005097319 
    ## ... Similar to previous best
    ## Run 482 stress 8.969359e-05 
    ## ... Procrustes: rmse 0.2527128  max resid 0.3452258 
    ## Run 483 stress 0.2601498 
    ## Run 484 stress 9.790445e-05 
    ## ... Procrustes: rmse 0.007136824  max resid 0.01140534 
    ## Run 485 stress 9.524245e-05 
    ## ... Procrustes: rmse 0.2474539  max resid 0.3381089 
    ## Run 486 stress 0.2933059 
    ## Run 487 stress 0.2601512 
    ## Run 488 stress 9.587177e-05 
    ## ... Procrustes: rmse 0.01041334  max resid 0.01659868 
    ## Run 489 stress 9.450739e-05 
    ## ... Procrustes: rmse 0.2527071  max resid 0.3452149 
    ## Run 490 stress 8.874493e-05 
    ## ... Procrustes: rmse 0.1582941  max resid 0.2261059 
    ## Run 491 stress 9.62752e-05 
    ## ... Procrustes: rmse 0.02357572  max resid 0.0373144 
    ## Run 492 stress 9.615578e-05 
    ## ... Procrustes: rmse 0.2107946  max resid 0.2913106 
    ## Run 493 stress 9.252492e-05 
    ## ... Procrustes: rmse 0.04538731  max resid 0.07065355 
    ## Run 494 stress 9.320111e-05 
    ## ... Procrustes: rmse 0.05662703  max resid 0.08741804 
    ## Run 495 stress 9.320785e-05 
    ## ... Procrustes: rmse 0.01487895  max resid 0.02373401 
    ## Run 496 stress 9.121298e-05 
    ## ... Procrustes: rmse 0.2527373  max resid 0.345205 
    ## Run 497 stress 0.0002800608 
    ## ... Procrustes: rmse 0.2514458  max resid 0.3556614 
    ## Run 498 stress 0.2264602 
    ## Run 499 stress 8.560302e-05 
    ## ... Procrustes: rmse 0.2524173  max resid 0.3448653 
    ## Run 500 stress 8.372196e-05 
    ## ... Procrustes: rmse 0.0006518015  max resid 0.0009822583 
    ## ... Similar to previous best
    ## *** Best solution repeated 8 times

    ## Warning in metaMDS(SD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, :
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

### SD beta env and geo correlated variables using envfit

``` r
### Environmental
# Surveyed sites 
# For figure
(SD_beta_env_ef <- envfit(SD_beta_env_NMDS, env[surveyed_sites_env,c(34)], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                                   NMDS1    NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)]  0.99781 -0.06607 0.7388  0.026 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_env_efp <- p.adjust.envfit(SD_beta_env_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                                   NMDS1    NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)]  0.99781 -0.06607 0.7388  0.026 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by stratification
(SD_beta_env_A_ef <- envfit(SD_beta_env_NMDS, env[surveyed_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.48790 -0.87290 0.0586  0.837  
    ## salinity_median     0.99781 -0.06607 0.7388  0.031 *
    ## oxygen_median       0.79303  0.60919 0.5362  0.816  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_env_A_efp <- p.adjust.envfit(SD_beta_env_A_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.48790 -0.87290 0.0586  1.000  
    ## salinity_median     0.99781 -0.06607 0.7388  0.093 .
    ## oxygen_median       0.79303  0.60919 0.5362  1.000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
(SD_beta_env_MS_ef <- envfit(SD_beta_env_MS_NMDS, env[mixed_stratified_lakes,environment], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.56538 -0.82483 0.0778  0.787  
    ## salinity_median     0.99173  0.12832 0.7458  0.039 *
    ## oxygen_median       0.94279  0.33338 0.3552  0.905  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_env_MS_efp <- p.adjust.envfit(SD_beta_env_MS_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median -0.56538 -0.82483 0.0778  1.000
    ## salinity_median     0.99173  0.12832 0.7458  0.117
    ## oxygen_median       0.94279  0.33338 0.3552  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
(SD_beta_env_OM_ef <- envfit(SD_beta_env_OM_NMDS, env[ocean_mixed_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)   
    ## temperature_median -0.10614  0.99435 0.0705  0.715   
    ## salinity_median    -0.78732 -0.61655 0.7842  0.009 **
    ## oxygen_median      -0.96712  0.25433 0.7188  0.052 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_env_OM_efp <- p.adjust.envfit(SD_beta_env_OM_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.10614  0.99435 0.0705  1.000  
    ## salinity_median    -0.78732 -0.61655 0.7842  0.027 *
    ## oxygen_median      -0.96712  0.25433 0.7188  0.156  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
(SD_beta_env_SO_ef <- envfit(SD_beta_env_SO_NMDS, env[ocean_stratified_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)  
    ## temperature_median  0.00019429  1.00000000 0.0911  0.510  
    ## salinity_median    -0.00044050 -1.00000000 0.6762  0.057 .
    ## oxygen_median      -0.00105747  1.00000000 0.7608  0.234  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_env_SO_efp <- p.adjust.envfit(SD_beta_env_SO_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median  0.00019429  1.00000000 0.0911  1.000
    ## salinity_median    -0.00044050 -1.00000000 0.6762  0.171
    ## oxygen_median      -0.00105747  1.00000000 0.7608  0.702
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed lakes
(SD_beta_env_M_ef <- envfit(SD_beta_env_M_NMDS, env[mixed_lakes,environment], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median -0.00010131  1.00000000 0.0841  0.808
    ## salinity_median     0.00075654  1.00000000 0.4002  0.277
    ## oxygen_median       0.00308549 -1.00000000 0.5830  0.122
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_env_M_efp <- p.adjust.envfit(SD_beta_env_M_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                          NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median -0.00010131  1.00000000 0.0841  1.000
    ## salinity_median     0.00075654  1.00000000 0.4002  0.831
    ## oxygen_median       0.00308549 -1.00000000 0.5830  0.366
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes
(SD_beta_S_ef <- envfit(SD_beta_S_NMDS, env[stratified_lakes,c(2,6,8,26,31:32)], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median       -0.35881 -0.93341 0.1965  0.614
    ## salinity_median          -0.70036  0.71379 0.5991  0.129
    ## oxygen_median             0.98543 -0.17006 0.5185  0.185
    ## distance_to_ocean_mean_m  0.01326 -0.99991 0.1726  0.574
    ## max_depth                -0.49997 -0.86604 0.0383  0.884
    ## logArea                  -0.26513 -0.96421 0.2144  0.558
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_S_efp <- p.adjust.envfit(SD_beta_S_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                             NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median       -0.35881 -0.93341 0.1965  1.000
    ## salinity_median          -0.70036  0.71379 0.5991  0.774
    ## oxygen_median             0.98543 -0.17006 0.5185  1.000
    ## distance_to_ocean_mean_m  0.01326 -0.99991 0.1726  1.000
    ## max_depth                -0.49997 -0.86604 0.0383  1.000
    ## logArea                  -0.26513 -0.96421 0.2144  1.000
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
# For figure
(SD_beta_geo_ef <- envfit(SD_beta_geo_NMDS, env[surveyed_sites_geo,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.98644 -0.16411 0.6448  0.124  
    ## max_depth               -0.15363 -0.98813 0.1504  0.628  
    ## logArea                  0.21786 -0.97598 0.2870  0.067 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_geo_efp <- p.adjust.envfit(SD_beta_geo_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.98644 -0.16411 0.6448  0.372
    ## max_depth               -0.15363 -0.98813 0.1504  1.000
    ## logArea                  0.21786 -0.97598 0.2870  0.201
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by stratification
(SD_beta_geo_A_ef <- envfit(SD_beta_geo_NMDS, env[surveyed_sites_geo,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.98644 -0.16411 0.6448  0.122  
    ## max_depth               -0.15363 -0.98813 0.1504  0.652  
    ## logArea                  0.21786 -0.97598 0.2870  0.089 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_geo_A_efp <- p.adjust.envfit(SD_beta_geo_A_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.98644 -0.16411 0.6448  0.366
    ## max_depth               -0.15363 -0.98813 0.1504  1.000
    ## logArea                  0.21786 -0.97598 0.2870  0.267
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
(SD_beta_geo_MS_ef <- envfit(SD_beta_geo_MS_NMDS, env[mixed_stratified_lakes_geo,geography], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.94998  0.31231 0.5323  0.149
    ## max_depth               -0.24666 -0.96910 0.1830  0.692
    ## logArea                 -0.88116  0.47282 0.0201  0.981
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_geo_MS_efp <- p.adjust.envfit(SD_beta_geo_MS_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.94998  0.31231 0.5323  0.447
    ## max_depth               -0.24666 -0.96910 0.1830  1.000
    ## logArea                 -0.88116  0.47282 0.0201  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
(SD_beta_geo_OM_ef <- envfit(SD_beta_geo_OM_NMDS, env[ocean_mixed_sites_geo,geography], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m  0.85078 -0.52552 0.4208  0.033 *
    ## max_depth               -0.70976 -0.70444 0.5550  0.028 *
    ## logArea                 -0.77935  0.62659 0.3227  0.241  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_geo_OM_efp <- p.adjust.envfit(SD_beta_geo_OM_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m  0.85078 -0.52552 0.4208  0.099 .
    ## max_depth               -0.70976 -0.70444 0.5550  0.084 .
    ## logArea                 -0.77935  0.62659 0.3227  0.723  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
(SD_beta_geo_SO_ef <- envfit(SD_beta_geo_SO_NMDS, env[ocean_stratified_sites,geography], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m  0.86963  0.49371 0.6793  0.192
    ## max_depth                0.68534  0.72823 0.0510  0.962
    ## logArea                 -0.52042  0.85391 0.2187  0.545
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_geo_SO_efp <- p.adjust.envfit(SD_beta_geo_SO_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m  0.86963  0.49371 0.6793  0.576
    ## max_depth                0.68534  0.72823 0.0510  1.000
    ## logArea                 -0.52042  0.85391 0.2187  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed lakes
(SD_beta_geo_M_ef <- envfit(SD_beta_geo_M_NMDS, env[mixed_lakes_geo,geography], permutations = 999, na.rm = TRUE))
```

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m -0.96071 -0.27756 0.7754  0.041 * 
    ## max_depth                0.85615  0.51672 0.8457  0.088 . 
    ## logArea                  0.19585  0.98063 0.8267  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 5039

``` r
(SD_beta_geo_M_efp <- p.adjust.envfit(SD_beta_geo_M_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m -0.96071 -0.27756 0.7754  0.123   
    ## max_depth                0.85615  0.51672 0.8457  0.264   
    ## logArea                  0.19585  0.98063 0.8267  0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 5039

### SD beta Mantel correlation tests

``` r
### Environmental
# Surveyed sites
env_dist_t <- dist(scaled_env[surveyed_sites_env,c(1)], method = "euclidean")
(SD_beta_env_mant_t <- mantel(SD_beta_env_dist$Btotal, env_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_dist$Btotal, ydis = env_dist_t, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.241 
    ##       Significance: 0.231 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.300 0.329 0.347 0.373 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_dist_s <- dist(scaled_env[surveyed_sites_env,c(2)], method = "euclidean")
(SD_beta_env_mant_s <- mantel(SD_beta_env_dist$Btotal, env_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_dist$Btotal, ydis = env_dist_s, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6697 
    ##       Significance: 0.006 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.566 0.597 0.631 0.654 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_dist_o <- dist(scaled_env[surveyed_sites_env,c(3)], method = "euclidean")
(SD_beta_env_mant_o <- mantel(SD_beta_env_dist$Btotal, env_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_dist$Btotal, ydis = env_dist_o, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r:  0.47 
    ##       Significance: 0.491 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.598 0.633 0.658 0.705 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_mant_pv <- rbind(SD_beta_env_mant_t$signif, SD_beta_env_mant_s$signif, SD_beta_env_mant_o$signif)
SD_beta_env_mant_pv <- SD_beta_env_mant_pv[,1]
(SD_beta_env_mant_pv <- p.adjust(SD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.693 0.018 1.000

``` r
# Mixed and stratified lakes
env_MS_dist_t <- dist(scaled_env[mixed_stratified_lakes,c(1)], method = "euclidean")
(SD_beta_env_MS_mant_t <- mantel(SD_beta_env_MS_dist$Btotal, env_MS_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_t,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2368 
    ##       Significance: 0.345 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.334 0.367 0.390 0.433 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_MS_dist_s <- dist(scaled_env[mixed_stratified_lakes,c(2)], method = "euclidean")
(SD_beta_env_MS_mant_s <- mantel(SD_beta_env_MS_dist$Btotal, env_MS_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_s,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6676 
    ##       Significance: 0.014 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.554 0.603 0.645 0.670 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_MS_dist_o <- dist(scaled_env[mixed_stratified_lakes,c(3)], method = "euclidean")
(SD_beta_env_MS_mant_o <- mantel(SD_beta_env_MS_dist$Btotal, env_MS_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_o,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2979 
    ##       Significance: 0.601 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.489 0.528 0.558 0.593 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_MS_mant_pv <- rbind(SD_beta_env_MS_mant_t$signif, SD_beta_env_MS_mant_s$signif, SD_beta_env_MS_mant_o$signif)
SD_beta_env_MS_mant_pv <- SD_beta_env_MS_mant_pv[,1]
(SD_beta_env_MS_mant_pv <- p.adjust(SD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.042 1.000

``` r
# Ocean sites and mixed lakes
env_OM_dist_t <- dist(scaled_env[ocean_mixed_sites_env,c(1)], method = "euclidean")
(SD_beta_env_OM_mant_t <- mantel(SD_beta_env_OM_dist$Btotal, env_OM_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_t,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1214 
    ##       Significance: 0.747 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.128 0.177 0.226 0.403 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_OM_dist_s <- dist(scaled_env[ocean_mixed_sites_env,c(2)], method = "euclidean")
(SD_beta_env_OM_mant_s <- mantel(SD_beta_env_OM_dist$Btotal, env_OM_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_s,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4258 
    ##       Significance: 0.055 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.359 0.435 0.485 0.545 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_OM_dist_o <- dist(scaled_env[ocean_mixed_sites_env,c(3)], method = "euclidean")
(SD_beta_env_OM_mant_o <- mantel(SD_beta_env_OM_dist$Btotal, env_OM_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_o,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5319 
    ##       Significance: 0.033 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.355 0.471 0.549 0.644 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_OM_mant_pv <- rbind(SD_beta_env_OM_mant_t$signif, SD_beta_env_OM_mant_s$signif, SD_beta_env_OM_mant_o$signif)
SD_beta_env_OM_mant_pv <- SD_beta_env_OM_mant_pv[,1]
(SD_beta_env_OM_mant_pv <- p.adjust(SD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.165 0.099

``` r
# Stratified lakes and ocean sites
env_SO_dist_t <- dist(scaled_env[ocean_stratified_sites_env,c(1)], method = "euclidean")
(SD_beta_env_SO_mant_t <- mantel(SD_beta_env_SO_dist$Btotal, env_SO_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_t,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.008188 
    ##       Significance: 0.222 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0271 0.0483 0.0648 0.0762 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_SO_dist_s <- dist(scaled_env[ocean_stratified_sites_env,c(2)], method = "euclidean")
(SD_beta_env_SO_mant_s <- mantel(SD_beta_env_SO_dist$Btotal, env_SO_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_s,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5528 
    ##       Significance: 0.015 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.460 0.501 0.537 0.564 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_SO_dist_o <- dist(scaled_env[ocean_stratified_sites_env,c(3)], method = "euclidean")
(SD_beta_env_SO_mant_o <- mantel(SD_beta_env_SO_dist$Btotal, env_SO_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_o,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6567 
    ##       Significance: 0.386 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.708 0.730 0.747 0.767 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_SO_mant_pv <- rbind(SD_beta_env_SO_mant_t$signif, SD_beta_env_SO_mant_s$signif, SD_beta_env_SO_mant_o$signif)
SD_beta_env_SO_mant_pv <- SD_beta_env_SO_mant_pv[,1]
(SD_beta_env_SO_mant_pv <- p.adjust(SD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 0.666 0.045 1.000

``` r
# Mixed lakes
env_M_dist_t <- dist(scaled_env[mixed_lakes,c(1)], method = "euclidean")
(SD_beta_env_M_mant_t <- mantel(SD_beta_env_M_dist$Btotal, env_M_dist_t, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_M_dist$Btotal, ydis = env_M_dist_t,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1352 
    ##       Significance: 0.738 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.251 0.353 0.450 0.677 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_M_dist_s <- dist(scaled_env[mixed_lakes,c(2)], method = "euclidean")
(SD_beta_env_M_mant_s <- mantel(SD_beta_env_M_dist$Btotal, env_M_dist_s, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_M_dist$Btotal, ydis = env_M_dist_s,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1877 
    ##       Significance: 0.152 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.249 0.312 0.357 0.396 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_M_dist_o <- dist(scaled_env[mixed_lakes,c(3)], method = "euclidean")
(SD_beta_env_M_mant_o <- mantel(SD_beta_env_M_dist$Btotal, env_M_dist_o, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_M_dist$Btotal, ydis = env_M_dist_o,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2781 
    ##       Significance: 0.106 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.296 0.457 0.543 0.668 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_M_mant_pv <- rbind(SD_beta_env_M_mant_t$signif, SD_beta_env_M_mant_s$signif, SD_beta_env_M_mant_o$signif)
SD_beta_env_M_mant_pv <- SD_beta_env_M_mant_pv[,1]
(SD_beta_env_M_mant_pv <- p.adjust(SD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.456 0.318

``` r
# Stratified lakes
env_geo_S_dist <- dist(scaled_env[stratified_lakes,c(1:3,9,14:15)], method = "euclidean")
(SD_beta_S_mant <- mantel(SD_beta_S_dist$Btotal, env_geo_S_dist, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_S_dist$Btotal, ydis = env_geo_S_dist, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.112 
    ##       Significance: 0.261 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.240 0.293 0.331 0.381 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
geo_dist_dmean <- dist(scaled_env[surveyed_sites_geo,c(9)], method = "euclidean")
(SD_beta_geo_mant_dmean <- mantel(SD_beta_geo_dist$Btotal, geo_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_dist$Btotal, ydis = geo_dist_dmean,      method = "spearman", permutations = 999, strata = env[surveyed_sites_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4572 
    ##       Significance: 0.406 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.554 0.591 0.620 0.642 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_md <- dist(scaled_env[surveyed_sites_geo,c(14)], method = "euclidean")
(SD_beta_geo_mant_md <- mantel(SD_beta_geo_dist$Btotal, geo_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_dist$Btotal, ydis = geo_dist_md, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.04407 
    ##       Significance: 0.683 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.193 0.223 0.255 0.284 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_la <- dist(scaled_env[surveyed_sites_geo,c(15)], method = "euclidean")
(SD_beta_geo_mant_la <- mantel(SD_beta_geo_dist$Btotal, geo_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_dist$Btotal, ydis = geo_dist_la, method = "spearman",      permutations = 999, strata = env[surveyed_sites_geo, 19],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.07044 
    ##       Significance: 0.469 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.124 0.138 0.148 0.156 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_mant_pv <- rbind(SD_beta_geo_mant_dmean$signif, SD_beta_geo_mant_md$signif, SD_beta_geo_mant_la$signif)
SD_beta_geo_mant_pv <- SD_beta_geo_mant_pv[,1]
(SD_beta_geo_mant_pv <- p.adjust(SD_beta_geo_mant_pv, method = "bonferroni"))
```

    ## [1] 1 1 1

``` r
# Mixed and stratified lakes 
geo_MS_dist_dmean <- dist(scaled_env[mixed_stratified_lakes_geo,c(9)], method = "euclidean")
(SD_beta_geo_MS_mant_dmean <- mantel(SD_beta_geo_MS_dist$Btotal, geo_MS_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_dmean,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2753 
    ##       Significance: 0.393 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.426 0.469 0.511 0.532 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_md <- dist(scaled_env[mixed_stratified_lakes_geo,c(14)], method = "euclidean")
(SD_beta_geo_MS_mant_md <- mantel(SD_beta_geo_MS_dist$Btotal, geo_MS_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_md,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02831 
    ##       Significance: 0.823 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.201 0.256 0.304 0.346 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_la <- dist(scaled_env[mixed_stratified_lakes_geo,c(15)], method = "euclidean")
(SD_beta_geo_MS_mant_la <- mantel(SD_beta_geo_MS_dist$Btotal, geo_MS_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_la,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02183 
    ##       Significance: 0.802 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.139 0.169 0.190 0.223 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_MS_mant_pv <- rbind(SD_beta_geo_MS_mant_dmean$signif, SD_beta_geo_MS_mant_md$signif, SD_beta_geo_MS_mant_la$signif)
SD_beta_geo_MS_mant_pv <- SD_beta_geo_MS_mant_pv[,1]
(SD_beta_geo_MS_mant_pv <- p.adjust(SD_beta_geo_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 1 1 1

``` r
# Ocean sites and mixed lakes
geo_OM_dist_dmean <- dist(scaled_env[ocean_mixed_sites_geo,c(9)], method = "euclidean")
(SD_beta_geo_OM_mant_dmean <- mantel(SD_beta_geo_OM_dist$Btotal, geo_OM_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_dmean,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.09197 
    ##       Significance: 0.525 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.257 0.284 0.315 0.345 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_md <- dist(scaled_env[ocean_mixed_sites_geo,c(14)], method = "euclidean")
(SD_beta_geo_OM_mant_md <- mantel(SD_beta_geo_OM_dist$Btotal, geo_OM_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_md,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2229 
    ##       Significance: 0.08 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.191 0.276 0.340 0.423 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_la <- dist(scaled_env[ocean_mixed_sites_geo,c(15)], method = "euclidean")
(SD_beta_geo_OM_mant_la <- mantel(SD_beta_geo_OM_dist$Btotal, geo_OM_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_geo,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_la,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_geo,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1894 
    ##       Significance: 0.122 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.201 0.235 0.263 0.293 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_OM_mant_pv <- rbind( SD_beta_geo_OM_mant_dmean$signif, SD_beta_geo_OM_mant_md$signif, SD_beta_geo_OM_mant_la$signif)
SD_beta_geo_OM_mant_pv <- SD_beta_geo_OM_mant_pv[,1]
(SD_beta_geo_OM_mant_pv <- p.adjust(SD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.240 0.366

``` r
# Stratified lakes and ocean sites
geo_SO_dist_dmean <- dist(scaled_env[ocean_stratified_sites,c(9)], method = "euclidean")
(SD_beta_geo_SO_mant_dmean <- mantel(SD_beta_geo_SO_dist$Btotal, geo_SO_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_dmean,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5858 
    ##       Significance: 0.429 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.678 0.701 0.718 0.745 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_SO_dist_md <- dist(scaled_env[ocean_stratified_sites,c(14)], method = "euclidean")
(SD_beta_geo_SO_mant_md <- mantel(SD_beta_geo_SO_dist$Btotal, geo_SO_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_md,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.007031 
    ##       Significance: 0.62 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.145 0.175 0.207 0.243 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_SO_dist_la <- dist(scaled_env[ocean_stratified_sites,c(15)], method = "euclidean")
(SD_beta_geo_SO_mant_la <- mantel(SD_beta_geo_SO_dist$Btotal, geo_SO_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_la,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2461 
    ##       Significance: 0.25 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.285 0.302 0.317 0.330 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_SO_mant_pv <- rbind( SD_beta_geo_SO_mant_dmean$signif, SD_beta_geo_SO_mant_md$signif, SD_beta_geo_SO_mant_la$signif)
SD_beta_geo_SO_mant_pv <- SD_beta_geo_SO_mant_pv[,1]
(SD_beta_geo_SO_mant_pv <- p.adjust(SD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.00 1.00 0.75

``` r
# Mixed lakes
geo_M_dist_dmean <- dist(scaled_env[mixed_lakes_geo,c(9)], method = "euclidean")
(SD_beta_geo_M_mant_dmean <- mantel(SD_beta_geo_M_dist$Btotal, geo_M_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_dmean,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.09227 
    ##       Significance: 0.614 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.318 0.423 0.550 0.684 
    ## Permutation: free
    ## Number of permutations: 5039

``` r
geo_M_dist_md <- dist(scaled_env[mixed_lakes_geo,c(14)], method = "euclidean")
(SD_beta_geo_M_mant_md <- mantel(SD_beta_geo_M_dist$Btotal, geo_M_dist_md, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_md,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.653 
    ##       Significance: 0.014 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.311 0.430 0.529 0.714 
    ## Permutation: free
    ## Number of permutations: 5039

``` r
geo_M_dist_la <- dist(scaled_env[mixed_lakes_geo,c(15)], method = "euclidean")
(SD_beta_geo_M_mant_la <- mantel(SD_beta_geo_M_dist$Btotal, geo_M_dist_la, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_la,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.387 
    ##       Significance: 0.038 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.247 0.343 0.447 0.544 
    ## Permutation: free
    ## Number of permutations: 5039

``` r
# Adjust p-values
SD_beta_geo_M_mant_pv <- rbind( SD_beta_geo_M_mant_dmean$signif, SD_beta_geo_M_mant_md$signif, SD_beta_geo_M_mant_la$signif)
SD_beta_geo_M_mant_pv <- SD_beta_geo_M_mant_pv[,1]
(SD_beta_geo_M_mant_pv <- p.adjust(SD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.042 0.114

### SD beta NMDS ordination plots

``` r
# SD beta total NMDS scores
SD_beta_NMDS_data.scores <- as.data.frame(scores(SD_beta_NMDS))
SD_beta_NMDS_data.scores$Stratification <- env[surveyed_sites,19]
SD_beta_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
SD_beta_NMDS_data.scores$Stratification <- factor(SD_beta_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# Get significantly correlated environmental variables
SD_beta_env_ef_coord_cont <- as.data.frame(scores(SD_beta_env_ef, "vectors")) * ordiArrowMul(SD_beta_env_ef)
SD_beta_env_row_names <- c("S")
# Assign the new row names to the data frame
SD_beta_env_ef_coord_cont <- data.frame(row.names = SD_beta_env_row_names, SD_beta_env_ef_coord_cont)

# SD_beta_geo_ef_coord_cont <- as.data.frame(scores(SD_beta_geo_ef, "vectors")) * ordiArrowMul(SD_beta_geo_ef)
# SD_beta_geo_row_names <- c("minD")
# # Assign the new row names to the data frame
# SD_beta_geo_ef_coord_cont <- data.frame(row.names = SD_beta_geo_row_names, SD_beta_geo_ef_coord_cont)

# Plot NMDS ordination of SD beta total with CI = 0.95
SD_beta_ef_plot <- ggplot(data = SD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.95) +
  geom_point(data = SD_beta_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = SD_beta_env_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = SD_beta_env_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
            color = env_cont, label = row.names(SD_beta_env_ef_coord_cont), size = 5) +
  geom_text_repel(data = SD_beta_NMDS_data.scores, label = SD_beta_NMDS_data.scores$Lakes, 
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
           label = paste("Stress: ", round(SD_beta_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81))
(SD_beta_ef_plot <- SD_beta_ef_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_beta_NMDS.jpg", SD_beta_ef_plot, width = 6, height = 6, units = "in")


# Plot NMDS ordination of SD beta total with CI = 0.90
SD_beta_ef_plot_CI90 <- ggplot(data = SD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.90) +
  geom_point(data = SD_beta_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = SD_beta_env_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = SD_beta_env_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
            color = env_cont, label = row.names(SD_beta_env_ef_coord_cont), size = 5) +
  geom_text_repel(data = SD_beta_NMDS_data.scores, label = SD_beta_NMDS_data.scores$Lakes, 
                  size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 16), 
        legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_rect(fill = NA), 
        axis.text = element_text(color = "black"), 
        legend.key = element_blank()) +
  annotate("text", x = -0.6, y = 0.5, size = 5,
           label = paste("Stress: ", round(SD_beta_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81))
(SD_beta_ef_plot_CI90 <- SD_beta_ef_plot_CI90 + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_beta_NMDS_CI90.jpg", SD_beta_ef_plot_CI90, width = 6, height = 6, units = "in")


# SD beta replacement NMDS scores
SD_beta_rep_NMDS_data.scores <- as.data.frame(scores(SD_beta_rep_NMDS))
SD_beta_rep_NMDS_data.scores$Stratification <- env[surveyed_sites,19]
SD_beta_rep_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
SD_beta_rep_NMDS_data.scores$Stratification <- factor(SD_beta_rep_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta replacement with CI = 0.95
SD_beta_rep_plot <- ggplot(data = SD_beta_rep_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.95) +
  geom_point(data = SD_beta_rep_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = SD_beta_rep_NMDS_data.scores, label = SD_beta_rep_NMDS_data.scores$Lakes, 
                  size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 14), 
        legend.position = "bottom", 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_rect(fill = NA), 
        axis.text = element_text(color = "black"), 
        legend.key = element_blank()) +
  annotate("text", x = -0.7, y = 0.8, size = 5,
           label = paste("Stress: ", round(SD_beta_rep_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "a") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(SD_beta_rep_plot <- SD_beta_rep_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_beta_rep_NMDS.jpg", SD_beta_rep_plot, width = 4.88, height = 6, units = "in")


# SD beta richness NMDS scores
SD_beta_ric_NMDS_data.scores <- as.data.frame(scores(SD_beta_ric_NMDS))
SD_beta_ric_NMDS_data.scores$Stratification <- env[surveyed_sites,19]
SD_beta_ric_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
SD_beta_ric_NMDS_data.scores$Stratification <- factor(SD_beta_ric_NMDS_data.scores$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta richness with CI = 0.95
SD_beta_ric_plot <- ggplot(data = SD_beta_ric_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.95) +
  geom_point(data = SD_beta_ric_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = SD_beta_ric_NMDS_data.scores, label = SD_beta_ric_NMDS_data.scores$Lakes, 
                  size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 14), 
        legend.position = "bottom", 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_rect(fill = NA), 
        axis.text = element_text(color = "black"), 
        legend.key = element_blank()) +
  annotate("text", x = -0.7, y = 0.8, size = 5,
           label = paste("Stress: ", round(SD_beta_ric_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "b") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(SD_beta_ric_plot <- SD_beta_ric_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-4.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_beta_ric_NMDS.jpg", SD_beta_ric_plot, width = 4.88, height = 6, units = "in")


# SD beta ref NMDS scores
SD_beta_ref_NMDS_data.scores <- as.data.frame(scores(SD_beta_ref_NMDS))
SD_beta_ref_NMDS_data.scores$Stratification <- env[,19]
SD_beta_ref_NMDS_data.scores$Lakes <- env[,1]
SD_beta_ref_NMDS_data.scores$Stratification <- factor(SD_beta_ref_NMDS_data.scores$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta ref with CI = 0.95
SD_beta_ref_plot <- ggplot(data = SD_beta_ref_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Stratification)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Stratification), type = "t", level = 0.95) +
  geom_point(data = SD_beta_ref_NMDS_data.scores, aes(color = Stratification, fill = Stratification), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = SD_beta_ref_NMDS_data.scores, label = SD_beta_ref_NMDS_data.scores$Lakes, 
                  size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 14), 
        legend.position = "bottom", 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.border = element_rect(fill = NA), 
        axis.text = element_text(color = "black"), 
        legend.key = element_blank()) +
  annotate("text", x = -0.7, y = 0.8, size = 5,
           label = paste("Stress: ", round(SD_beta_ref_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "a") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(SD_beta_ref_plot <- SD_beta_ref_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_path()`).

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-5.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_beta_ref_NMDS.jpg", SD_beta_ref_plot, width = 4.88, height = 6, units = "in")
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_path()`).

### SD beta dendrogram

``` r
# Cluster communities using average-linkage algorithm
SD_beta_dist_clust <- hclust(SD_beta_dist$Btotal, method = "average")

# The following will create a dendrogram of SD_beta putting things in one dimensional space
# Open a jpg device
png("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_beta_dendrogram.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

# Create your plot using the plot() function
plot(SD_beta_dist_clust, 
     xlab = "Sites",
     ylab = "SD beta",
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
png("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_beta_dg_silhouette.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

Si <- numeric(nrow(presabs_lake[surveyed_sites,]))
for (k in 2:(nrow(presabs_lake[surveyed_sites,])-1))
{
  sil<-silhouette(cutree(SD_beta_dist_clust, k=k), SD_beta_dist$Btotal)
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
png("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_beta_dg_partitioning.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

SD_beta_pam <- pam(SD_beta_dist$Btotal,2)
plot(SD_beta_pam)

# Close the jpg device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### SD beta outliers

``` r
# Determine NMDS1 outliers
ordered_SD_beta_NMDS1_data.scores <- SD_beta_NMDS_data.scores[order(SD_beta_NMDS_data.scores$NMDS1), ]

outlier_SD_beta_NMDS1 <- ordered_SD_beta_NMDS1_data.scores %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(NMDS1, 0.25),
    Q3 = quantile(NMDS1, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = NMDS1 < lower_bound | NMDS1 > upper_bound
  )

outlier_SD_beta_NMDS1 <- as.data.frame(outlier_SD_beta_NMDS1)

row.names(outlier_SD_beta_NMDS1) <- outlier_SD_beta_NMDS1$X

# Create the plot
(outlier_SD_beta_NMDS1_plot <- ggplot(outlier_SD_beta_NMDS1, aes(x = Stratification, y = NMDS1, fill = Stratification)) +
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

![](SD_analyses_files/figure-gfm/SD%20beta%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/outlier_SD_beta_NMDS1.jpg", outlier_SD_beta_NMDS1_plot, width = 6.26, height = 6, units = "in")


# Determine NMDS2 outliers
ordered_SD_beta_NMDS2_data.scores <- SD_beta_NMDS_data.scores[order(SD_beta_NMDS_data.scores$NMDS2), ]

outlier_SD_beta_NMDS2 <- ordered_SD_beta_NMDS2_data.scores %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(NMDS2, 0.25),
    Q3 = quantile(NMDS2, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = NMDS2 < lower_bound | NMDS2 > upper_bound
  )

outlier_SD_beta_NMDS2 <- as.data.frame(outlier_SD_beta_NMDS2)

row.names(outlier_SD_beta_NMDS2) <- outlier_SD_beta_NMDS2$X

# Create the plot
(outlier_SD_beta_NMDS2_plot <- ggplot(outlier_SD_beta_NMDS2, aes(x = Stratification, y = NMDS2, fill = Stratification)) +
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

![](SD_analyses_files/figure-gfm/SD%20beta%20outliers-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/outlier_SD_beta_NMDS2.jpg", outlier_SD_beta_NMDS2_plot, width = 6.26, height = 6, units = "in")


# Make SD_beta distances into a matrix
SD_beta_dist_matrix <- as.matrix(SD_beta_dist$Btotal)

# Set the diagonal elements to NA
diag(SD_beta_dist_matrix) <- NA

# Get the labels of the sites
sites <- attr(SD_beta_dist_matrix, "Labels")

# Calculate the mean dist for each group
SD_beta_mean_dist <- aggregate(SD_beta_dist_matrix, by = list(stratification_group), FUN = mean, na.rm = TRUE)

# Rename rows, delete column, and transpose matrix and make it a data frame
row.names(SD_beta_mean_dist) <- SD_beta_mean_dist$Group.1
SD_beta_mean_dist <- SD_beta_mean_dist[,-1] 
SD_beta_mean_dist <- as.data.frame(t(SD_beta_mean_dist))

# Add Stratification column
SD_beta_mean_dist$Stratification <- env[surveyed_sites,19]
SD_beta_mean_dist$Stratification <- factor(SD_beta_mean_dist$Stratification, levels = c("Ocean", "Mixed", "Stratified"))


# Determine ocean sites outliers based on distance from other ocean sites
outlier_SD_beta_mean_dist_ocean <- SD_beta_mean_dist[ocean_sites,] %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(Ocean, 0.25),
    Q3 = quantile(Ocean, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = Ocean < lower_bound | Ocean > upper_bound
  )

outlier_SD_beta_mean_dist_ocean <- as.data.frame(outlier_SD_beta_mean_dist_ocean)
row.names(outlier_SD_beta_mean_dist_ocean) <- outlier_SD_beta_mean_dist_ocean$X
outlier_SD_beta_mean_dist_ocean$Distances <- outlier_SD_beta_mean_dist_ocean$Ocean
outlier_SD_beta_mean_dist_ocean <- outlier_SD_beta_mean_dist_ocean[,-c(1:3)]


# Determine mixed lake outliers based on distance from other mixed lakes
outlier_SD_beta_mean_dist_mixed <- SD_beta_mean_dist[mixed_lakes,] %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(Mixed, 0.25),
    Q3 = quantile(Mixed, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = Mixed < lower_bound | Mixed > upper_bound
  )

outlier_SD_beta_mean_dist_mixed <- as.data.frame(outlier_SD_beta_mean_dist_mixed)
row.names(outlier_SD_beta_mean_dist_mixed) <- outlier_SD_beta_mean_dist_mixed$X
outlier_SD_beta_mean_dist_mixed$Distances <- outlier_SD_beta_mean_dist_mixed$Mixed
outlier_SD_beta_mean_dist_mixed <- outlier_SD_beta_mean_dist_mixed[,-c(1:3)]


# Determine stratified lake outliers based on distance from other stratified lakes
outlier_SD_beta_mean_dist_stratified <- SD_beta_mean_dist[stratified_lakes,] %>%
  group_by(Stratification) %>%
  mutate(
    Q1 = quantile(Stratified, 0.25),
    Q3 = quantile(Stratified, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = Stratified < lower_bound | Stratified > upper_bound
  )

outlier_SD_beta_mean_dist_stratified <- as.data.frame(outlier_SD_beta_mean_dist_stratified)
row.names(outlier_SD_beta_mean_dist_stratified) <- outlier_SD_beta_mean_dist_stratified$X
outlier_SD_beta_mean_dist_stratified$Distances <- outlier_SD_beta_mean_dist_stratified$Stratified
outlier_SD_beta_mean_dist_stratified <- outlier_SD_beta_mean_dist_stratified[,-c(1:3)]

outlier_SD_beta_mean_dist <- rbind(outlier_SD_beta_mean_dist_ocean, outlier_SD_beta_mean_dist_mixed, outlier_SD_beta_mean_dist_stratified)

(outlier_SD_beta_mean_dist_plot <- ggplot(outlier_SD_beta_mean_dist, aes(x = Stratification, y = Distances, fill = Stratification)) +
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

![](SD_analyses_files/figure-gfm/SD%20beta%20outliers-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/outlier_SD_beta_mean_dist.jpg", outlier_SD_beta_mean_dist_plot, width = 6.26, height = 6, units = "in")
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
    ## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggvenn_0.1.10        pairwiseAdonis_0.4.1 cluster_2.1.8       
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
    ##  [52] quantreg_6.00           MASS_7.3-64             lava_1.8.1             
    ##  [55] scatterplot3d_0.3-44    ModelMetrics_1.2.2.2    tools_4.4.3            
    ##  [58] foreign_0.8-88          future.apply_1.11.3     nnet_7.3-20            
    ##  [61] glue_1.8.0              quadprog_1.5-8          checkmate_2.3.2        
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

#### Determine what families of fish are most present in the lakes

``` r
library(rfishbase)
species_sums_tax <- site_sums
species_sums_tax$X <- gsub("_", " ", species_sums_tax$X)
species_sums_tax <- species_sums_tax$X
species_sums_tax <- rfishbase::load_taxa() %>% 
  filter(Species %in% species_sums_tax) %>%
  collect()
species_sums_tax$Species <- gsub(" ", "_", species_sums_tax$Species)
species_sums_tax <- merge(species_sums_tax, species_sums, by.x = "Species", by.y = "X")
species_sums_tax <- species_sums_tax[, -c(2:4, 7:8)]
species_sums_tax <- species_sums_tax[order(species_sums_tax$row_sum),]

species_sums_tax_summary <- species_sums_tax %>%
  group_by(Family) %>%
  summarise(across(.cols = row_sum, .fns = sum, na.rm = TRUE))

family <- ggplot(data = species_sums_tax_summary) + 
  geom_bar(aes(x = reorder(Family, row_sum), y = row_sum),
    stat = 'identity',
    col="green", 
    fill="green4", 
    alpha = .5) + 
  theme_bw() +
  theme(plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.line = element_line(color = "black")) +
  labs(x="Family", y="Family Richness")
family
ggsave("family.jpg", family, width = 5, height = 4.5)

order <- ggplot(data = species_sums_tax) + 
  geom_bar(aes(x = reorder(Order, row_sum), y = row_sum),
    stat = 'identity',
    col="green", 
    fill="green4", 
    alpha = .5) + 
  theme_bw() +
  theme(plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.line = element_line(color = "black")) +
  labs(x="Order", y="Richness in Lakess")
order
#ggsave("order.jpg", order)

detach("package:rfishbase")
```
