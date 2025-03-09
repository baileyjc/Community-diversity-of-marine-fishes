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
SD_logSRic_env_am <- aov(log(row_sum) ~ (temperature_median + salinity_median + oxygen_median) * Stratification, data = SR_env[surveyed_sites_env,])
summary(SD_logSRic_env_am)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median                 1  2.599   2.599   7.928   0.0226 *  
    ## salinity_median                    1 19.667  19.667  59.986 5.51e-05 ***
    ## oxygen_median                      1  2.564   2.564   7.822   0.0233 *  
    ## Stratification                     2  0.144   0.072   0.219   0.8079    
    ## temperature_median:Stratification  2  0.230   0.115   0.350   0.7148    
    ## salinity_median:Stratification     2  0.318   0.159   0.485   0.6328    
    ## oxygen_median:Stratification       1  0.346   0.346   1.056   0.3342    
    ## Residuals                          8  2.623   0.328                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.1585043341 0.0003856326 0.1632139530 1.0000000000 1.0000000000
    ## [6] 1.0000000000 1.0000000000           NA

``` r
SD_logSRic_env_lm <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[surveyed_sites_env,])
summary(SD_logSRic_env_lm)
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
p_values <- summary(SD_logSRic_env_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##       2.537061e-01       1.000000e+00       2.109736e-05       2.190101e-02

``` r
SD_SRic_env_am <- aov(row_sum ~ (temperature_median + salinity_median + oxygen_median) * Stratification, data = SR_env[surveyed_sites_env,])
summary(SD_SRic_env_am)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)   
    ## temperature_median                 1   1618    1618   4.120 0.0769 . 
    ## salinity_median                    1   7596    7596  19.335 0.0023 **
    ## oxygen_median                      1   3809    3809   9.695 0.0144 * 
    ## Stratification                     2     71      36   0.091 0.9144   
    ## temperature_median:Stratification  2   1276     638   1.624 0.2559   
    ## salinity_median:Stratification     2    344     172   0.438 0.6601   
    ## oxygen_median:Stratification       1    667     667   1.699 0.2287   
    ## Residuals                          8   3143     393                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.53828754 0.01606838 0.10055228 1.00000000 1.00000000 1.00000000 1.00000000
    ## [8]         NA

``` r
SD_SRic_env_lm <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[surveyed_sites_env,])
summary(SD_SRic_env_lm)
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
p_values <- summary(SD_SRic_env_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##         1.00000000         1.00000000         0.04950985         0.02277051

``` r
# Mixed and stratified lakes
SD_logSRic_env_MS_am <- aov(log(row_sum) ~ (temperature_median + salinity_median + oxygen_median) * Stratification, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_env_MS_am)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median                 1  2.891   2.891   8.819 0.017881 *  
    ## salinity_median                    1 15.131  15.131  46.152 0.000139 ***
    ## oxygen_median                      1  2.180   2.180   6.650 0.032682 *  
    ## Stratification                     1  0.119   0.119   0.364 0.563202    
    ## temperature_median:Stratification  1  0.190   0.190   0.580 0.468197    
    ## salinity_median:Stratification     1  0.304   0.304   0.927 0.363923    
    ## oxygen_median:Stratification       1  0.346   0.346   1.056 0.334200    
    ## Residuals                          8  2.623   0.328                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_MS_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.1251675438 0.0009711235 0.2287765520 1.0000000000 1.0000000000
    ## [6] 1.0000000000 1.0000000000           NA

``` r
SD_logSRic_env_MS_lm <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_env_MS_lm)
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
anova(SD_logSRic_env_MS_lm)
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
p_values <- summary(SD_logSRic_env_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##       0.4476551535       1.0000000000       0.0001925785       0.0768885411

``` r
SD_SRic_env_MS_am <- aov(row_sum ~ (temperature_median + salinity_median + oxygen_median) * Stratification, data = SR_env[mixed_stratified_lakes,])
summary(SD_SRic_env_MS_am)
```

    ##                                   Df Sum Sq Mean Sq F value  Pr(>F)   
    ## temperature_median                 1   2185    2185   5.561 0.04609 * 
    ## salinity_median                    1   4875    4875  12.409 0.00782 **
    ## oxygen_median                      1   2840    2840   7.230 0.02755 * 
    ## Stratification                     1      1       1   0.002 0.96894   
    ## temperature_median:Stratification  1    875     875   2.228 0.17387   
    ## salinity_median:Stratification     1    305     305   0.777 0.40377   
    ## oxygen_median:Stratification       1    667     667   1.699 0.22872   
    ## Residuals                          8   3143     393                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_MS_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.3226481 0.0547052 0.1928173 1.0000000 1.0000000 1.0000000 1.0000000
    ## [8]        NA

``` r
SD_SRic_env_MS_lm <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[mixed_stratified_lakes,])
summary(SD_SRic_env_MS_lm)
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
anova(SD_SRic_env_MS_lm)
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
p_values <- summary(SD_SRic_env_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##         1.00000000         1.00000000         0.08902474         0.09067730

``` r
# Ocean sites and mixed lakes
SD_logSRic_env_OM_am <- aov(log(row_sum) ~ (temperature_median + salinity_median + oxygen_median) * Stratification, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_env_OM_am)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median                 1 0.0031  0.0031   0.009  0.928
    ## salinity_median                    1 1.1073  1.1073   3.322  0.142
    ## oxygen_median                      1 0.2981  0.2981   0.894  0.398
    ## Stratification                     1 0.2353  0.2353   0.706  0.448
    ## temperature_median:Stratification  1 0.0279  0.0279   0.084  0.787
    ## salinity_median:Stratification     1 0.0177  0.0177   0.053  0.829
    ## Residuals                          4 1.3333  0.3333

``` r
p_values <- summary(SD_logSRic_env_OM_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 0.8546885 1.0000000 1.0000000 1.0000000 1.0000000        NA

``` r
SD_logSRic_env_OM_lm <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_env_OM_lm)
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
p_values <- summary(SD_logSRic_env_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
SD_SRic_env_OM_am <- aov(row_sum ~ (temperature_median + salinity_median + oxygen_median) * Stratification, data = SR_env[ocean_mixed_sites_env,])
summary(SD_SRic_env_OM_am)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median                 1   60.4    60.4   0.079  0.793
    ## salinity_median                    1 1651.6  1651.6   2.160  0.216
    ## oxygen_median                      1 1106.2  1106.2   1.447  0.295
    ## Stratification                     1  569.3   569.3   0.745  0.437
    ## temperature_median:Stratification  1  548.7   548.7   0.718  0.445
    ## salinity_median:Stratification     1   29.7    29.7   0.039  0.853
    ## Residuals                          4 3058.3   764.6

``` r
p_values <- summary(SD_SRic_env_OM_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1]  1  1  1  1  1  1 NA

``` r
SD_SRic_env_OM_lm <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[ocean_mixed_sites_env,])
summary(SD_SRic_env_OM_lm)
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
p_values <- summary(SD_SRic_env_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          1.0000000          1.0000000          0.8678638

``` r
# Stratified lakes and ocean sites
SD_logSRic_env_SO_am <- aov(log(row_sum) ~ (temperature_median + salinity_median + oxygen_median) * Stratification, data = SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_env_SO_am)
```

    ##                                   Df Sum Sq Mean Sq F value  Pr(>F)   
    ## temperature_median                 1  0.303   0.303   0.939 0.38732   
    ## salinity_median                    1 12.450  12.450  38.620 0.00341 **
    ## oxygen_median                      1  1.737   1.737   5.389 0.08101 . 
    ## Stratification                     1  0.398   0.398   1.235 0.32870   
    ## temperature_median:Stratification  1  0.122   0.122   0.378 0.57184   
    ## salinity_median:Stratification     1  0.028   0.028   0.086 0.78337   
    ## Residuals                          4  1.290   0.322                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_SO_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.00000000 0.02047409 0.48604559 1.00000000 1.00000000 1.00000000         NA

``` r
SD_logSRic_env_SO_lm <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_env_SO_lm)
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
p_values <- summary(SD_logSRic_env_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##        0.316349339        1.000000000        0.002511461        0.147465030

``` r
SD_SRic_env_SO_am <- aov(row_sum ~ (temperature_median + salinity_median + oxygen_median) * Stratification, data = SR_env[ocean_stratified_sites_env,])
summary(SD_SRic_env_SO_am)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median                 1     84      84   3.972 0.117049    
    ## salinity_median                    1   4598    4598 217.490 0.000123 ***
    ## oxygen_median                      1   2482    2482 117.384 0.000412 ***
    ## Stratification                     1    171     171   8.088 0.046684 *  
    ## temperature_median:Stratification  1    610     610  28.877 0.005793 ** 
    ## salinity_median:Stratification     1     67      67   3.159 0.150163    
    ## Residuals                          4     85      21                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_SO_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.7022941073 0.0007382902 0.0024706703 0.2801063330 0.0347580162
    ## [6] 0.9009791825           NA

``` r
SD_SRic_env_SO_lm <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[ocean_stratified_sites_env,])
summary(SD_SRic_env_SO_lm)
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
p_values <- summary(SD_SRic_env_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##         0.25595351         1.00000000         0.01289859         0.01399934

``` r
# Ocean sites
SD_logSRic_env_O_lm <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[ocean_sites,])
summary(SD_logSRic_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         15.9560        NaN     NaN      NaN
    ## temperature_median   0.3599        NaN     NaN      NaN
    ## salinity_median     -0.6833        NaN     NaN      NaN
    ## oxygen_median            NA         NA      NA       NA
    ## 
    ## Residual standard error: NaN on 0 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
anova(SD_logSRic_env_O_lm)
```

    ## Warning in anova.lm(SD_logSRic_env_O_lm): ANOVA F-tests on an essentially
    ## perfect fit are unreliable

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                    Df   Sum Sq  Mean Sq F value Pr(>F)
    ## temperature_median  1 0.150753 0.150753     NaN    NaN
    ## salinity_median     1 0.013834 0.013834     NaN    NaN
    ## Residuals           0 0.000000      NaN

``` r
p_values <- summary(SD_logSRic_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median 
    ##                NaN                NaN                NaN

``` r
SD_SRic_env_O_lm <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[ocean_sites,])
summary(SD_SRic_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median, data = SR_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)          798.05        NaN     NaN      NaN
    ## temperature_median    24.62        NaN     NaN      NaN
    ## salinity_median      -44.50        NaN     NaN      NaN
    ## oxygen_median            NA         NA      NA       NA
    ## 
    ## Residual standard error: NaN on 0 degrees of freedom
    ##   (3 observations deleted due to missingness)
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
anova(SD_SRic_env_O_lm)
```

    ## Warning in anova.lm(SD_SRic_env_O_lm): ANOVA F-tests on an essentially perfect
    ## fit are unreliable

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                    Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1 713.98  713.98     NaN    NaN
    ## salinity_median     1  58.69   58.69     NaN    NaN
    ## Residuals           0   0.00     NaN

``` r
p_values <- summary(SD_SRic_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median 
    ##                NaN                NaN                NaN

``` r
# Mixed lakes
SD_logSRic_env_M_lm <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[mixed_lakes,])
summary(SD_logSRic_env_M_lm)
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
anova(SD_logSRic_env_M_lm)
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
p_values <- summary(SD_logSRic_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
SD_SRic_env_M_lm <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[mixed_lakes,])
summary(SD_SRic_env_M_lm)
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
anova(SD_SRic_env_M_lm)
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
p_values <- summary(SD_SRic_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
# Stratified lakes
SD_logSRic_env_S_lm <- lm(log(row_sum) ~ temperature_median + salinity_median + oxygen_median, SR_env[stratified_lakes,])
summary(SD_logSRic_env_S_lm)
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
anova(SD_logSRic_env_S_lm)
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
p_values <- summary(SD_logSRic_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##          1.0000000          1.0000000          0.8657002          1.0000000

``` r
SD_SRic_env_S_lm <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median, SR_env[stratified_lakes,])
summary(SD_SRic_env_S_lm)
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
anova(SD_SRic_env_S_lm)
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
p_values <- summary(SD_SRic_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept) temperature_median    salinity_median      oxygen_median 
    ##                  1                  1                  1                  1

``` r
### Geographical
# Surveyed sites
SD_logSRic_geo_am <- aov(log(row_sum) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = SR_env[surveyed_sites,])
summary(SD_logSRic_geo_am)
```

    ##                                        Df Sum Sq Mean Sq F value  Pr(>F)    
    ## distance_to_ocean_min_m                 1 17.419  17.419  46.808 4.5e-05 ***
    ## max_depth                               1  0.021   0.021   0.055 0.81908    
    ## logArea                                 1  0.136   0.136   0.366 0.55882    
    ## Stratification                          2  7.445   3.723  10.003 0.00411 ** 
    ## distance_to_ocean_min_m:Stratification  2  0.218   0.109   0.292 0.75273    
    ## max_depth:Stratification                2  0.505   0.253   0.679 0.52907    
    ## logArea:Stratification                  2  0.625   0.313   0.840 0.45998    
    ## Residuals                              10  3.721   0.372                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_geo_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
SD_logSRic_geo_lm <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[surveyed_sites,])
summary(SD_logSRic_geo_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[surveyed_sites, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7491 -0.4274  0.0483  0.3671  1.6760 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              3.519588   1.161407   3.030 0.007193 ** 
    ## distance_to_ocean_min_m -0.011413   0.002619  -4.358 0.000379 ***
    ## max_depth               -0.001511   0.019251  -0.078 0.938314    
    ## logArea                  0.051058   0.115403   0.442 0.663444    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8338 on 18 degrees of freedom
    ## Multiple R-squared:  0.5841, Adjusted R-squared:  0.5148 
    ## F-statistic: 8.426 on 3 and 18 DF,  p-value: 0.001038

``` r
p_values <- summary(SD_logSRic_geo_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             0.028770087             0.001515849             1.000000000 
    ##                 logArea 
    ##             1.000000000

``` r
SD_SRic_geo_am <- aov(row_sum ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = SR_env[surveyed_sites,])
summary(SD_SRic_geo_am)
```

    ##                                        Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m                 1   6762    6762  20.745 0.00105 **
    ## max_depth                               1    428     428   1.314 0.27838   
    ## logArea                                 1    699     699   2.146 0.17367   
    ## Stratification                          2   5631    2816   8.638 0.00662 **
    ## distance_to_ocean_min_m:Stratification  2    149      74   0.229 0.79972   
    ## max_depth:Stratification                2   1554     777   2.383 0.14244   
    ## logArea:Stratification                  2   1527     764   2.343 0.14638   
    ## Residuals                              10   3260     326                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
SD_SRic_geo_lm <- lm(row_sum ~  max_depth + logArea * Stratification, SR_env[surveyed_sites,])
summary(SD_SRic_geo_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ max_depth + logArea * Stratification, 
    ##     data = SR_env[surveyed_sites, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -27.6617  -7.6425   0.1206   9.7574  23.2474 
    ## 
    ## Coefficients:
    ##                                   Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                        24.5607    33.3088   0.737   0.4723  
    ## max_depth                           0.8677     0.4735   1.833   0.0868 .
    ## logArea                             1.3683     2.8789   0.475   0.6414  
    ## StratificationMixed              -111.7171    52.8812  -2.113   0.0518 .
    ## StratificationStratified           46.0307    69.7863   0.660   0.5195  
    ## logArea:StratificationMixed        12.0791     5.0741   2.381   0.0310 *
    ## logArea:StratificationStratified   -9.5152     6.8430  -1.391   0.1847  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 16.04 on 15 degrees of freedom
    ## Multiple R-squared:  0.8071, Adjusted R-squared:   0.73 
    ## F-statistic: 10.46 on 6 and 15 DF,  p-value: 0.0001214

``` r
p_values <- summary(SD_SRic_geo_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##                      (Intercept)                        max_depth 
    ##                        1.0000000                        0.6075203 
    ##                          logArea              StratificationMixed 
    ##                        1.0000000                        0.3626745 
    ##         StratificationStratified      logArea:StratificationMixed 
    ##                        1.0000000                        0.2168482 
    ## logArea:StratificationStratified 
    ##                        1.0000000

``` r
# Mixed and stratified lakes
SD_logSRic_geo_MS_am <- aov(log(row_sum) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_geo_MS_am)
```

    ##                                        Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m                 1 12.814  12.814  34.386 0.000377 ***
    ## max_depth                               1  0.075   0.075   0.202 0.664772    
    ## logArea                                 1  2.999   2.999   8.048 0.021915 *  
    ## Stratification                          1  4.624   4.624  12.410 0.007814 ** 
    ## distance_to_ocean_min_m:Stratification  1  0.000   0.000   0.000 0.991222    
    ## max_depth:Stratification                1  0.246   0.246   0.659 0.440430    
    ## logArea:Stratification                  1  0.046   0.046   0.123 0.734712    
    ## Residuals                               8  2.981   0.373                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_geo_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
SD_logSRic_geo_MS_lm <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_geo_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[mixed_stratified_lakes, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6649 -0.3645  0.1885  0.4100  0.9013 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)              0.109909   2.024761   0.054  0.95760   
    ## distance_to_ocean_min_m -0.012883   0.003083  -4.178  0.00128 **
    ## max_depth               -0.048875   0.027370  -1.786  0.09942 . 
    ## logArea                  0.503683   0.235942   2.135  0.05409 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8112 on 12 degrees of freedom
    ## Multiple R-squared:  0.668,  Adjusted R-squared:  0.585 
    ## F-statistic: 8.048 on 3 and 12 DF,  p-value: 0.003321

``` r
p_values <- summary(SD_logSRic_geo_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             1.000000000             0.005123546             0.397682450 
    ##                 logArea 
    ##             0.216372794

``` r
SD_SRic_geo_MS_am <- aov(row_sum ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = SR_env[mixed_stratified_lakes,])
summary(SD_SRic_geo_MS_am)
```

    ##                                        Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m                 1   4715    4715  21.003 0.00180 **
    ## max_depth                               1      1       1   0.005 0.94705   
    ## logArea                                 1   4426    4426  19.718 0.00217 **
    ## Stratification                          1   2800    2800  12.474 0.00771 **
    ## distance_to_ocean_min_m:Stratification  1    110     110   0.491 0.50343   
    ## max_depth:Stratification                1    703     703   3.131 0.11476   
    ## logArea:Stratification                  1    340     340   1.514 0.25342   
    ## Residuals                               8   1796     224                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
SD_SRic_geo_MS_lm <- lm(row_sum ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, SR_env[mixed_stratified_lakes,])
summary(SD_SRic_geo_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ (distance_to_ocean_min_m + max_depth + 
    ##     logArea) * Stratification, data = SR_env[mixed_stratified_lakes, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -28.0861  -5.3118   0.2086   4.9038  22.0988 
    ## 
    ## Coefficients:
    ##                                                   Estimate Std. Error t value
    ## (Intercept)                                      -99.47259   41.77624  -2.381
    ## distance_to_ocean_min_m                           -0.08445    0.24707  -0.342
    ## max_depth                                          0.09009    1.07945   0.083
    ## logArea                                           16.35700    5.54093   2.952
    ## StratificationStratified                         108.13429  101.28289   1.068
    ## distance_to_ocean_min_m:StratificationStratified   0.02938    0.25922   0.113
    ## max_depth:StratificationStratified                -0.10032    1.55089  -0.065
    ## logArea:StratificationStratified                 -15.61396   12.68790  -1.231
    ##                                                  Pr(>|t|)  
    ## (Intercept)                                        0.0445 *
    ## distance_to_ocean_min_m                            0.7413  
    ## max_depth                                          0.9355  
    ## logArea                                            0.0184 *
    ## StratificationStratified                           0.3168  
    ## distance_to_ocean_min_m:StratificationStratified   0.9126  
    ## max_depth:StratificationStratified                 0.9500  
    ## logArea:StratificationStratified                   0.2534  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 14.98 on 8 degrees of freedom
    ## Multiple R-squared:  0.8794, Adjusted R-squared:  0.7739 
    ## F-statistic: 8.334 on 7 and 8 DF,  p-value: 0.003854

``` r
p_values <- summary(SD_SRic_geo_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##                                      (Intercept) 
    ##                                        0.3557648 
    ##                          distance_to_ocean_min_m 
    ##                                        1.0000000 
    ##                                        max_depth 
    ##                                        1.0000000 
    ##                                          logArea 
    ##                                        0.1469491 
    ##                         StratificationStratified 
    ##                                        1.0000000 
    ## distance_to_ocean_min_m:StratificationStratified 
    ##                                        1.0000000 
    ##               max_depth:StratificationStratified 
    ##                                        1.0000000 
    ##                 logArea:StratificationStratified 
    ##                                        1.0000000

``` r
# Ocean sites and mixed lakes
SD_logSRic_geo_OM_am <- aov(log(row_sum) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_geo_OM_am)
```

    ##                                        Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1 0.0712  0.0712   0.309 0.5981  
    ## max_depth                               1 1.6259  1.6259   7.063 0.0376 *
    ## logArea                                 1 0.1182  0.1182   0.514 0.5005  
    ## Stratification                          1 0.0121  0.0121   0.052 0.8266  
    ## distance_to_ocean_min_m:Stratification  1 0.0101  0.0101   0.044 0.8410  
    ## max_depth:Stratification                1 0.1230  0.1230   0.534 0.4923  
    ## logArea:Stratification                  1 0.6182  0.6182   2.686 0.1524  
    ## Residuals                               6 1.3812  0.2302                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_geo_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
SD_logSRic_geo_OM_lm <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_mixed_sites,])
summary(SD_logSRic_geo_OM_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[ocean_mixed_sites, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.65264 -0.22920 -0.08159  0.33157  0.59872 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)             2.803905   0.773535   3.625  0.00465 **
    ## distance_to_ocean_min_m 0.001322   0.003526   0.375  0.71548   
    ## max_depth               0.033177   0.014335   2.314  0.04318 * 
    ## logArea                 0.051836   0.069817   0.742  0.47489   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4631 on 10 degrees of freedom
    ## Multiple R-squared:  0.4584, Adjusted R-squared:  0.296 
    ## F-statistic: 2.822 on 3 and 10 DF,  p-value: 0.09321

``` r
p_values <- summary(SD_logSRic_geo_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              0.01861009              1.00000000              0.17272621 
    ##                 logArea 
    ##              1.00000000

``` r
SD_SRic_geo_OM_am <- aov(row_sum ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = SR_env[ocean_mixed_sites,])
summary(SD_SRic_geo_OM_am)
```

    ##                                        Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m                 1     14      14   0.027 0.8751  
    ## max_depth                               1   3536    3536   6.853 0.0397 *
    ## logArea                                 1    702     702   1.361 0.2876  
    ## Stratification                          1     16      16   0.032 0.8641  
    ## distance_to_ocean_min_m:Stratification  1      2       2   0.005 0.9476  
    ## max_depth:Stratification                1    120     120   0.232 0.6468  
    ## logArea:Stratification                  1   1509    1509   2.925 0.1381  
    ## Residuals                               6   3096     516                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_OM_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
SD_SRic_geo_OM_lm <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_mixed_sites,])
summary(SD_SRic_geo_OM_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[ocean_mixed_sites, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -32.006  -9.850  -4.391  14.406  36.018 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -14.0425    36.3794  -0.386   0.7076  
    ## distance_to_ocean_min_m   0.1520     0.1658   0.917   0.3809  
    ## max_depth                 1.4328     0.6742   2.125   0.0595 .
    ## logArea                   3.9957     3.2835   1.217   0.2516  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 21.78 on 10 degrees of freedom
    ## Multiple R-squared:  0.4727, Adjusted R-squared:  0.3145 
    ## F-statistic: 2.988 on 3 and 10 DF,  p-value: 0.08246

``` r
p_values <- summary(SD_SRic_geo_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               1.0000000               0.2379787 
    ##                 logArea 
    ##               1.0000000

``` r
# Stratified lakes and ocean sites
SD_logSRic_geo_SO_am <- aov(log(row_sum) ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_geo_SO_am)
```

    ##                                        Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m                 1 14.817  14.817  28.858 0.00171 **
    ## max_depth                               1  0.265   0.265   0.515 0.49978   
    ## logArea                                 1  0.161   0.161   0.313 0.59604   
    ## Stratification                          1  1.849   1.849   3.601 0.10653   
    ## distance_to_ocean_min_m:Stratification  1  0.009   0.009   0.017 0.90111   
    ## max_depth:Stratification                1  0.162   0.162   0.316 0.59413   
    ## logArea:Stratification                  1  0.041   0.041   0.079 0.78809   
    ## Residuals                               6  3.081   0.513                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_geo_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
SD_logSRic_geo_SO_lm <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_stratified_sites,])
summary(SD_logSRic_geo_SO_lm)
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
p_values <- summary(SD_logSRic_geo_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##             0.203110617             0.005020788             1.000000000 
    ##                 logArea 
    ##             1.000000000

``` r
SD_SRic_geo_SO_am <- aov(row_sum ~ (distance_to_ocean_min_m + max_depth + logArea) * Stratification, data = SR_env[ocean_stratified_sites,])
summary(SD_SRic_geo_SO_am)
```

    ##                                        Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m                 1   5520    5520  20.349 0.00406 **
    ## max_depth                               1    432     432   1.594 0.25355   
    ## logArea                                 1    293     293   1.080 0.33883   
    ## Stratification                          1   1572    1572   5.794 0.05279 . 
    ## distance_to_ocean_min_m:Stratification  1      1       1   0.005 0.94809   
    ## max_depth:Stratification                1    853     853   3.146 0.12647   
    ## logArea:Stratification                  1      3       3   0.009 0.92647   
    ## Residuals                               6   1627     271                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_SO_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
SD_SRic_geo_SO_lm <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_stratified_sites,])
summary(SD_SRic_geo_SO_lm)
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
p_values <- summary(SD_SRic_geo_SO_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##              1.00000000              0.04954771              1.00000000 
    ##                 logArea 
    ##              1.00000000

``` r
# Ocean sites
SD_logSRic_geo_M_lm <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_sites,])
summary(SD_logSRic_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[ocean_sites, ])
    ## 
    ## Residuals:
    ##        IBK        LCN        NCN        OOO        OOM        RCA 
    ## -1.430e-03 -2.889e-01  8.726e-02 -4.591e-01  6.622e-01  2.767e-17 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              3.835905   1.558146   2.462    0.133
    ## distance_to_ocean_min_m -0.007523   0.037011  -0.203    0.858
    ## max_depth                0.031838   0.027959   1.139    0.373
    ## logArea                 -0.033553   0.138658  -0.242    0.831
    ## 
    ## Residual standard error: 0.6084 on 2 degrees of freedom
    ## Multiple R-squared:  0.4282, Adjusted R-squared:  -0.4294 
    ## F-statistic: 0.4993 on 3 and 2 DF,  p-value: 0.7198

``` r
anova(SD_logSRic_geo_M_lm)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1 0.01808 0.01808  0.0488 0.8456
    ## max_depth                1 0.51470 0.51470  1.3905 0.3596
    ## logArea                  1 0.02168 0.02168  0.0586 0.8313
    ## Residuals                2 0.74033 0.37016

``` r
p_values <- summary(SD_logSRic_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.5315567               1.0000000               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
SD_SRic_geo_M_lm <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_sites,])
summary(SD_SRic_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[ocean_sites, ])
    ## 
    ## Residuals:
    ##        IBK        LCN        NCN        OOO        OOM        RCA 
    ##  5.607e-01 -1.064e+01  2.356e+00 -2.178e+01  2.950e+01  1.451e-15 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              38.2137    69.2809   0.552    0.637
    ## distance_to_ocean_min_m  -0.4885     1.6456  -0.297    0.795
    ## max_depth                 1.5906     1.2432   1.279    0.329
    ## logArea                  -0.5173     6.1653  -0.084    0.941
    ## 
    ## Residual standard error: 27.05 on 2 degrees of freedom
    ## Multiple R-squared:  0.5087, Adjusted R-squared:  -0.2284 
    ## F-statistic: 0.6901 on 3 and 2 DF,  p-value: 0.6372

``` r
anova(SD_SRic_geo_M_lm)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                         Df  Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1    1.63    1.63  0.0022 0.9666
    ## max_depth                1 1508.40 1508.40  2.0612 0.2876
    ## logArea                  1    5.15    5.15  0.0070 0.9408
    ## Residuals                2 1463.65  731.82

``` r
p_values <- summary(SD_SRic_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                       1                       1                       1 
    ##                 logArea 
    ##                       1

``` r
# Mixed lakes
SD_logSRic_geo_M_lm <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[mixed_lakes,])
summary(SD_logSRic_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##       FLK       HLO       LLN       MLN       NLN       NLU       OLO       ULN 
    ##  0.251959 -0.196322  0.284201 -0.008349 -0.024467  0.138767 -0.634468  0.188679 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)              0.958065   1.116068   0.858    0.439
    ## distance_to_ocean_min_m -0.003250   0.006601  -0.492    0.648
    ## max_depth                0.009723   0.028838   0.337    0.753
    ## logArea                  0.308001   0.148028   2.081    0.106
    ## 
    ## Residual standard error: 0.4003 on 4 degrees of freedom
    ## Multiple R-squared:  0.7592, Adjusted R-squared:  0.5787 
    ## F-statistic: 4.205 on 3 and 4 DF,  p-value: 0.09952

``` r
anova(SD_logSRic_geo_M_lm)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                         Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## distance_to_ocean_min_m  1 0.18157 0.18157  1.1332 0.34708  
    ## max_depth                1 1.14585 1.14585  7.1518 0.05556 .
    ## logArea                  1 0.69363 0.69363  4.3293 0.10594  
    ## Residuals                4 0.64087 0.16022                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               1.0000000               1.0000000 
    ##                 logArea 
    ##               0.4237593

``` r
SD_SRic_geo_M_lm <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[mixed_lakes,])
summary(SD_SRic_geo_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ distance_to_ocean_min_m + max_depth + 
    ##     logArea, data = SR_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##      FLK      HLO      LLN      MLN      NLN      NLU      OLO      ULN 
    ##   3.4234 -10.0735  22.0988  -0.7218  -5.3111  13.7320 -28.0861   4.9383 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -99.47259   56.32146  -1.766   0.1521  
    ## distance_to_ocean_min_m  -0.08445    0.33310  -0.254   0.8123  
    ## max_depth                 0.09009    1.45528   0.062   0.9536  
    ## logArea                  16.35700    7.47012   2.190   0.0937 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 20.2 on 4 degrees of freedom
    ## Multiple R-squared:  0.7286, Adjusted R-squared:  0.525 
    ## F-statistic: 3.579 on 3 and 4 DF,  p-value: 0.1249

``` r
anova(SD_SRic_geo_M_lm)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                         Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## distance_to_ocean_min_m  1   94.74   94.74  0.2322 0.65507  
    ## max_depth                1 2329.80 2329.80  5.7101 0.07520 .
    ## logArea                  1 1956.27 1956.27  4.7946 0.09373 .
    ## Residuals                4 1632.06  408.02                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.6084792               1.0000000               1.0000000 
    ##                 logArea 
    ##               0.3749120

``` r
# Stratified lakes
SD_logSRic_geo_S_lm <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[stratified_lakes,])
summary(SD_logSRic_geo_S_lm)
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
anova(SD_logSRic_geo_S_lm)
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
p_values <- summary(SD_logSRic_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               0.7668993               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
SD_SRic_geo_S_lm <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[stratified_lakes,])
summary(SD_SRic_geo_S_lm)
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
anova(SD_SRic_geo_S_lm)
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
p_values <- summary(SD_SRic_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
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
SD_beta_ref_dist <- BAT::beta(presabs_lake, abund = F)
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
    ## 3 Stratified vs Reference  1 0.6261283 1.909357 0.2143091   0.113      0.678
    ## 4          Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.047      0.282
    ## 5      Mixed vs Reference  1 0.5979203 1.966332 0.2193017   0.110      0.660
    ## 6      Ocean vs Reference  1 0.5599217 1.628213 0.2456489   0.146      0.876
    ##   sig
    ## 1   *
    ## 2   *
    ## 3    
    ## 4    
    ## 5    
    ## 6

``` r
# Surveyed sites
SD_beta_dist <- BAT::beta(surveyed_sites_lake, abund = F)
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
(SD_beta_PM <- adonis2(SD_beta_dist$Btotal ~ env[surveyed_sites,19], permutations = 999, method = "euclidean"))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_dist$Btotal ~ env[surveyed_sites, 19], permutations = 999, method = "euclidean")
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
    ## 3      Mixed vs Ocean  1 0.6396174 2.135322 0.1625633   0.007      0.021   .

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
    ## 2 Stratified vs Ocean  1 1.3260066 4.147587 0.2931657   0.003      0.009   *
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.043      0.129

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
SD_beta_geo_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites,], abund = F)
# Mixed and stratified lakes
SD_beta_geo_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], abund = F)
# Ocean sites and mixed lakes
SD_beta_geo_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites,], abund = F)
# Stratified lakes and ocean sites
SD_beta_geo_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites,], abund = F)
# Mixed lakes
SD_beta_geo_M_dist<- BAT::beta(surveyed_sites_lake[mixed_lakes,], abund = F)
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
    ## Run 1 stress 0.1050144 
    ## ... Procrustes: rmse 0.006661692  max resid 0.02376129 
    ## Run 2 stress 0.1050144 
    ## ... Procrustes: rmse 0.006664834  max resid 0.02374565 
    ## Run 3 stress 0.113594 
    ## Run 4 stress 0.1049903 
    ## ... New best solution
    ## ... Procrustes: rmse 8.923336e-05  max resid 0.0002414561 
    ## ... Similar to previous best
    ## Run 5 stress 0.1079381 
    ## Run 6 stress 0.1135941 
    ## Run 7 stress 0.1079381 
    ## Run 8 stress 0.1097987 
    ## Run 9 stress 0.1174536 
    ## Run 10 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214578  max resid 0.07901451 
    ## Run 11 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214119  max resid 0.07900233 
    ## Run 12 stress 0.1084965 
    ## Run 13 stress 0.107938 
    ## Run 14 stress 0.1051402 
    ## ... Procrustes: rmse 0.02308045  max resid 0.07849322 
    ## Run 15 stress 0.1088547 
    ## Run 16 stress 0.1138396 
    ## Run 17 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214228  max resid 0.07900419 
    ## Run 18 stress 0.1102722 
    ## Run 19 stress 0.1051402 
    ## ... Procrustes: rmse 0.02308358  max resid 0.07849688 
    ## Run 20 stress 0.1082722 
    ## Run 21 stress 0.1051402 
    ## ... Procrustes: rmse 0.02308706  max resid 0.07851602 
    ## Run 22 stress 0.1079384 
    ## Run 23 stress 0.1136646 
    ## Run 24 stress 0.1049903 
    ## ... Procrustes: rmse 8.526299e-05  max resid 0.0002406072 
    ## ... Similar to previous best
    ## Run 25 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 0.006739336  max resid 0.0258252 
    ## Run 26 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284419  max resid 0.0773373 
    ## Run 27 stress 0.1102721 
    ## Run 28 stress 0.1138986 
    ## Run 29 stress 0.1138392 
    ## Run 30 stress 0.1135596 
    ## Run 31 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284535  max resid 0.07735169 
    ## Run 32 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181619  max resid 0.07691643 
    ## Run 33 stress 0.1136646 
    ## Run 34 stress 0.1088547 
    ## Run 35 stress 0.1050144 
    ## ... Procrustes: rmse 0.009674972  max resid 0.0261667 
    ## Run 36 stress 0.1063128 
    ## Run 37 stress 0.1084964 
    ## Run 38 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284838  max resid 0.07731969 
    ## Run 39 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 2.665064e-05  max resid 8.02093e-05 
    ## ... Similar to previous best
    ## Run 40 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181082  max resid 0.07689681 
    ## Run 41 stress 0.113594 
    ## Run 42 stress 0.1063128 
    ## Run 43 stress 0.1082722 
    ## Run 44 stress 0.1049903 
    ## ... Procrustes: rmse 0.006731431  max resid 0.02582758 
    ## Run 45 stress 0.1063128 
    ## Run 46 stress 0.1135598 
    ## Run 47 stress 0.1063128 
    ## Run 48 stress 0.1050144 
    ## ... Procrustes: rmse 0.009652554  max resid 0.02604809 
    ## Run 49 stress 0.1135597 
    ## Run 50 stress 0.1136647 
    ## Run 51 stress 0.113594 
    ## Run 52 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180424  max resid 0.07685902 
    ## Run 53 stress 0.113594 
    ## Run 54 stress 0.110272 
    ## Run 55 stress 0.117454 
    ## Run 56 stress 0.1082722 
    ## Run 57 stress 0.1309381 
    ## Run 58 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180195  max resid 0.07686563 
    ## Run 59 stress 0.1135596 
    ## Run 60 stress 0.1063128 
    ## Run 61 stress 0.1063128 
    ## Run 62 stress 0.1136646 
    ## Run 63 stress 0.1136646 
    ## Run 64 stress 0.1135598 
    ## Run 65 stress 0.1083994 
    ## Run 66 stress 0.1050144 
    ## ... Procrustes: rmse 0.009679991  max resid 0.02620041 
    ## Run 67 stress 0.1084964 
    ## Run 68 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180616  max resid 0.07688032 
    ## Run 69 stress 0.113594 
    ## Run 70 stress 0.1160549 
    ## Run 71 stress 0.1135597 
    ## Run 72 stress 0.1165296 
    ## Run 73 stress 0.1084965 
    ## Run 74 stress 0.1097987 
    ## Run 75 stress 0.1083995 
    ## Run 76 stress 0.1082722 
    ## Run 77 stress 0.1083995 
    ## Run 78 stress 0.1063128 
    ## Run 79 stress 0.1165307 
    ## Run 80 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.159295e-05  max resid 2.838021e-05 
    ## ... Similar to previous best
    ## Run 81 stress 0.1063128 
    ## Run 82 stress 0.1063128 
    ## Run 83 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180307  max resid 0.07685489 
    ## Run 84 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180355  max resid 0.07686777 
    ## Run 85 stress 0.1136646 
    ## Run 86 stress 0.1050338 
    ## ... Procrustes: rmse 0.00627445  max resid 0.02196159 
    ## Run 87 stress 0.1051402 
    ## ... Procrustes: rmse 0.021808  max resid 0.07687803 
    ## Run 88 stress 0.1138983 
    ## Run 89 stress 0.1136647 
    ## Run 90 stress 0.1084965 
    ## Run 91 stress 0.1063128 
    ## Run 92 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180186  max resid 0.07684831 
    ## Run 93 stress 0.1174532 
    ## Run 94 stress 0.1174532 
    ## Run 95 stress 0.1058423 
    ## Run 96 stress 0.1063128 
    ## Run 97 stress 0.1135596 
    ## Run 98 stress 0.1084964 
    ## Run 99 stress 0.1084964 
    ## Run 100 stress 0.1136648 
    ## Run 101 stress 0.1050144 
    ## ... Procrustes: rmse 0.009660312  max resid 0.02609132 
    ## Run 102 stress 0.1082722 
    ## Run 103 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179965  max resid 0.07684972 
    ## Run 104 stress 0.1136646 
    ## Run 105 stress 0.117454 
    ## Run 106 stress 0.1135598 
    ## Run 107 stress 0.1058423 
    ## Run 108 stress 0.1049903 
    ## ... Procrustes: rmse 0.00673531  max resid 0.02584488 
    ## Run 109 stress 0.113594 
    ## Run 110 stress 0.1079381 
    ## Run 111 stress 0.1082722 
    ## Run 112 stress 0.1138987 
    ## Run 113 stress 0.107938 
    ## Run 114 stress 0.1138395 
    ## Run 115 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283738  max resid 0.07730028 
    ## Run 116 stress 0.10835 
    ## Run 117 stress 0.1088547 
    ## Run 118 stress 0.1138986 
    ## Run 119 stress 0.1079384 
    ## Run 120 stress 0.1097989 
    ## Run 121 stress 0.1079382 
    ## Run 122 stress 0.1097988 
    ## Run 123 stress 0.107938 
    ## Run 124 stress 0.1160548 
    ## Run 125 stress 0.1135596 
    ## Run 126 stress 0.1049903 
    ## ... Procrustes: rmse 0.006678855  max resid 0.02558739 
    ## Run 127 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283848  max resid 0.07730799 
    ## Run 128 stress 0.1135599 
    ## Run 129 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 9.908662e-06  max resid 2.835957e-05 
    ## ... Similar to previous best
    ## Run 130 stress 0.1049792 
    ## ... Procrustes: rmse 4.291626e-05  max resid 0.0001118639 
    ## ... Similar to previous best
    ## Run 131 stress 0.1135598 
    ## Run 132 stress 0.1135598 
    ## Run 133 stress 0.1083501 
    ## Run 134 stress 0.1136647 
    ## Run 135 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283851  max resid 0.07731024 
    ## Run 136 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181396  max resid 0.07689829 
    ## Run 137 stress 0.1082722 
    ## Run 138 stress 0.1049903 
    ## ... Procrustes: rmse 0.006689627  max resid 0.0256294 
    ## Run 139 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180419  max resid 0.07686982 
    ## Run 140 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283874  max resid 0.07728984 
    ## Run 141 stress 0.1160541 
    ## Run 142 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284019  max resid 0.07729793 
    ## Run 143 stress 0.1136647 
    ## Run 144 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180274  max resid 0.07685698 
    ## Run 145 stress 0.1174543 
    ## Run 146 stress 0.1050338 
    ## ... Procrustes: rmse 0.006256143  max resid 0.02185867 
    ## Run 147 stress 0.1063128 
    ## Run 148 stress 0.1135597 
    ## Run 149 stress 0.1049791 
    ## ... Procrustes: rmse 2.670435e-05  max resid 6.728505e-05 
    ## ... Similar to previous best
    ## Run 150 stress 0.1084965 
    ## Run 151 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283789  max resid 0.07731424 
    ## Run 152 stress 0.1102723 
    ## Run 153 stress 0.1084965 
    ## Run 154 stress 0.1049791 
    ## ... Procrustes: rmse 2.023148e-05  max resid 5.697942e-05 
    ## ... Similar to previous best
    ## Run 155 stress 0.1082722 
    ## Run 156 stress 0.1063128 
    ## Run 157 stress 0.1138391 
    ## Run 158 stress 0.1135599 
    ## Run 159 stress 0.1138985 
    ## Run 160 stress 0.1083499 
    ## Run 161 stress 0.1136646 
    ## Run 162 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180284  max resid 0.07685731 
    ## Run 163 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180227  max resid 0.0768565 
    ## Run 164 stress 0.1138989 
    ## Run 165 stress 0.1063128 
    ## Run 166 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180279  max resid 0.07686803 
    ## Run 167 stress 0.1082722 
    ## Run 168 stress 0.1082722 
    ## Run 169 stress 0.1309398 
    ## Run 170 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283877  max resid 0.07730708 
    ## Run 171 stress 0.1097988 
    ## Run 172 stress 0.10835 
    ## Run 173 stress 0.1136646 
    ## Run 174 stress 0.1083998 
    ## Run 175 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180882  max resid 0.07687838 
    ## Run 176 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180414  max resid 0.07686309 
    ## Run 177 stress 0.1049792 
    ## ... Procrustes: rmse 5.122861e-05  max resid 0.0001581903 
    ## ... Similar to previous best
    ## Run 178 stress 0.1138393 
    ## Run 179 stress 0.1088547 
    ## Run 180 stress 0.1084965 
    ## Run 181 stress 0.1135597 
    ## Run 182 stress 0.1083499 
    ## Run 183 stress 0.1309391 
    ## Run 184 stress 0.1088547 
    ## Run 185 stress 0.1083993 
    ## Run 186 stress 0.1058423 
    ## Run 187 stress 0.1084965 
    ## Run 188 stress 0.1063128 
    ## Run 189 stress 0.1049791 
    ## ... Procrustes: rmse 2.620453e-05  max resid 5.678416e-05 
    ## ... Similar to previous best
    ## Run 190 stress 0.1084964 
    ## Run 191 stress 0.1136646 
    ## Run 192 stress 0.1050144 
    ## ... Procrustes: rmse 0.009661471  max resid 0.02608227 
    ## Run 193 stress 0.1049791 
    ## ... Procrustes: rmse 3.662614e-06  max resid 8.87618e-06 
    ## ... Similar to previous best
    ## Run 194 stress 0.10835 
    ## Run 195 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283804  max resid 0.07727682 
    ## Run 196 stress 0.1050338 
    ## ... Procrustes: rmse 0.006264045  max resid 0.02191251 
    ## Run 197 stress 0.1063128 
    ## Run 198 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283771  max resid 0.07728562 
    ## Run 199 stress 0.1135596 
    ## Run 200 stress 0.1049791 
    ## ... Procrustes: rmse 3.071702e-05  max resid 9.283852e-05 
    ## ... Similar to previous best
    ## Run 201 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283878  max resid 0.07728905 
    ## Run 202 stress 0.1084965 
    ## Run 203 stress 0.1063128 
    ## Run 204 stress 0.1136647 
    ## Run 205 stress 0.1063128 
    ## Run 206 stress 0.1174538 
    ## Run 207 stress 0.1049792 
    ## ... Procrustes: rmse 5.593824e-05  max resid 0.000163418 
    ## ... Similar to previous best
    ## Run 208 stress 0.1138393 
    ## Run 209 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180071  max resid 0.07683012 
    ## Run 210 stress 0.1082722 
    ## Run 211 stress 0.1084965 
    ## Run 212 stress 0.1079381 
    ## Run 213 stress 0.1136646 
    ## Run 214 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284165  max resid 0.07731566 
    ## Run 215 stress 0.1079381 
    ## Run 216 stress 0.1088548 
    ## Run 217 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283672  max resid 0.07726979 
    ## Run 218 stress 0.1049903 
    ## ... Procrustes: rmse 0.006757291  max resid 0.02593565 
    ## Run 219 stress 0.1050144 
    ## ... Procrustes: rmse 0.009688232  max resid 0.02622279 
    ## Run 220 stress 0.1049791 
    ## ... Procrustes: rmse 2.532077e-05  max resid 7.053014e-05 
    ## ... Similar to previous best
    ## Run 221 stress 0.1063128 
    ## Run 222 stress 0.1082722 
    ## Run 223 stress 0.1049903 
    ## ... Procrustes: rmse 0.006689056  max resid 0.02562057 
    ## Run 224 stress 0.1165308 
    ## Run 225 stress 0.1063128 
    ## Run 226 stress 0.1135596 
    ## Run 227 stress 0.1049792 
    ## ... Procrustes: rmse 5.868495e-05  max resid 0.0001900296 
    ## ... Similar to previous best
    ## Run 228 stress 0.1135597 
    ## Run 229 stress 0.1165308 
    ## Run 230 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284289  max resid 0.07724457 
    ## Run 231 stress 0.1102722 
    ## Run 232 stress 0.1063128 
    ## Run 233 stress 0.1050144 
    ## ... Procrustes: rmse 0.009661645  max resid 0.0260875 
    ## Run 234 stress 0.1049791 
    ## ... Procrustes: rmse 2.193644e-05  max resid 4.842802e-05 
    ## ... Similar to previous best
    ## Run 235 stress 0.11653 
    ## Run 236 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283946  max resid 0.07730152 
    ## Run 237 stress 0.1079383 
    ## Run 238 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180646  max resid 0.07687137 
    ## Run 239 stress 0.1050144 
    ## ... Procrustes: rmse 0.009664766  max resid 0.02610494 
    ## Run 240 stress 0.1049903 
    ## ... Procrustes: rmse 0.006690623  max resid 0.02563384 
    ## Run 241 stress 0.1079381 
    ## Run 242 stress 0.1135598 
    ## Run 243 stress 0.1083502 
    ## Run 244 stress 0.1058708 
    ## Run 245 stress 0.1102723 
    ## Run 246 stress 0.1135596 
    ## Run 247 stress 0.1138399 
    ## Run 248 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180526  max resid 0.0768671 
    ## Run 249 stress 0.1136647 
    ## Run 250 stress 0.1136647 
    ## Run 251 stress 0.1050338 
    ## ... Procrustes: rmse 0.006255351  max resid 0.02185568 
    ## Run 252 stress 0.10835 
    ## Run 253 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180283  max resid 0.07686077 
    ## Run 254 stress 0.1063128 
    ## Run 255 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283923  max resid 0.07731451 
    ## Run 256 stress 0.1082722 
    ## Run 257 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228388  max resid 0.07730944 
    ## Run 258 stress 0.1088547 
    ## Run 259 stress 0.1050144 
    ## ... Procrustes: rmse 0.009661693  max resid 0.02606755 
    ## Run 260 stress 0.113594 
    ## Run 261 stress 0.1084965 
    ## Run 262 stress 0.1050144 
    ## ... Procrustes: rmse 0.009662901  max resid 0.02609029 
    ## Run 263 stress 0.1174528 
    ## Run 264 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180082  max resid 0.07682799 
    ## Run 265 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283979  max resid 0.07730813 
    ## Run 266 stress 0.1050144 
    ## ... Procrustes: rmse 0.009693937  max resid 0.0262814 
    ## Run 267 stress 0.107938 
    ## Run 268 stress 0.1088547 
    ## Run 269 stress 0.1082722 
    ## Run 270 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283944  max resid 0.07732392 
    ## Run 271 stress 0.1063128 
    ## Run 272 stress 0.113594 
    ## Run 273 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283837  max resid 0.07728116 
    ## Run 274 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180344  max resid 0.07686181 
    ## Run 275 stress 0.1084965 
    ## Run 276 stress 0.1050144 
    ## ... Procrustes: rmse 0.00966462  max resid 0.02608931 
    ## Run 277 stress 0.1088547 
    ## Run 278 stress 0.1082722 
    ## Run 279 stress 0.107938 
    ## Run 280 stress 0.1079382 
    ## Run 281 stress 0.1049903 
    ## ... Procrustes: rmse 0.00670788  max resid 0.02570992 
    ## Run 282 stress 0.1102722 
    ## Run 283 stress 0.113594 
    ## Run 284 stress 0.1135597 
    ## Run 285 stress 0.1050144 
    ## ... Procrustes: rmse 0.009659809  max resid 0.02608467 
    ## Run 286 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180897  max resid 0.07687859 
    ## Run 287 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180049  max resid 0.07682522 
    ## Run 288 stress 0.1058423 
    ## Run 289 stress 0.1088547 
    ## Run 290 stress 0.1084965 
    ## Run 291 stress 0.1050144 
    ## ... Procrustes: rmse 0.009648854  max resid 0.02603158 
    ## Run 292 stress 0.1049903 
    ## ... Procrustes: rmse 0.006753847  max resid 0.02591958 
    ## Run 293 stress 0.1135598 
    ## Run 294 stress 0.1102723 
    ## Run 295 stress 0.1058705 
    ## Run 296 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284089  max resid 0.07731041 
    ## Run 297 stress 0.1050338 
    ## ... Procrustes: rmse 0.006255793  max resid 0.02185943 
    ## Run 298 stress 0.1083501 
    ## Run 299 stress 0.1138395 
    ## Run 300 stress 0.1136646 
    ## Run 301 stress 0.1083499 
    ## Run 302 stress 0.1160545 
    ## Run 303 stress 0.1160558 
    ## Run 304 stress 0.1063129 
    ## Run 305 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179658  max resid 0.07683999 
    ## Run 306 stress 0.1136646 
    ## Run 307 stress 0.1049791 
    ## ... Procrustes: rmse 2.420696e-05  max resid 6.285405e-05 
    ## ... Similar to previous best
    ## Run 308 stress 0.1138392 
    ## Run 309 stress 0.1050144 
    ## ... Procrustes: rmse 0.009661291  max resid 0.02608801 
    ## Run 310 stress 0.1136646 
    ## Run 311 stress 0.1084965 
    ## Run 312 stress 0.1135596 
    ## Run 313 stress 0.1084964 
    ## Run 314 stress 0.1174545 
    ## Run 315 stress 0.1082722 
    ## Run 316 stress 0.1083994 
    ## Run 317 stress 0.1084964 
    ## Run 318 stress 0.1058423 
    ## Run 319 stress 0.1049903 
    ## ... Procrustes: rmse 0.006796604  max resid 0.02612165 
    ## Run 320 stress 0.1136646 
    ## Run 321 stress 0.1135597 
    ## Run 322 stress 0.1063128 
    ## Run 323 stress 0.1136646 
    ## Run 324 stress 0.1136646 
    ## Run 325 stress 0.1063128 
    ## Run 326 stress 0.1088548 
    ## Run 327 stress 0.1138393 
    ## Run 328 stress 0.1138392 
    ## Run 329 stress 0.1083988 
    ## Run 330 stress 0.1160544 
    ## Run 331 stress 0.1079381 
    ## Run 332 stress 0.1135599 
    ## Run 333 stress 0.1135596 
    ## Run 334 stress 0.113594 
    ## Run 335 stress 0.1063128 
    ## Run 336 stress 0.1050338 
    ## ... Procrustes: rmse 0.006268475  max resid 0.02191242 
    ## Run 337 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283959  max resid 0.07730548 
    ## Run 338 stress 0.1063128 
    ## Run 339 stress 0.1135598 
    ## Run 340 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283945  max resid 0.07731455 
    ## Run 341 stress 0.1102722 
    ## Run 342 stress 0.1135598 
    ## Run 343 stress 0.1049791 
    ## ... Procrustes: rmse 1.262201e-05  max resid 3.588133e-05 
    ## ... Similar to previous best
    ## Run 344 stress 0.1084965 
    ## Run 345 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180921  max resid 0.07687988 
    ## Run 346 stress 0.1165295 
    ## Run 347 stress 0.1079382 
    ## Run 348 stress 0.113594 
    ## Run 349 stress 0.1135598 
    ## Run 350 stress 0.1049903 
    ## ... Procrustes: rmse 0.006747635  max resid 0.02589085 
    ## Run 351 stress 0.1138398 
    ## Run 352 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180223  max resid 0.07685041 
    ## Run 353 stress 0.1063128 
    ## Run 354 stress 0.1050338 
    ## ... Procrustes: rmse 0.006254059  max resid 0.02185653 
    ## Run 355 stress 0.1049791 
    ## ... Procrustes: rmse 8.230155e-06  max resid 1.960895e-05 
    ## ... Similar to previous best
    ## Run 356 stress 0.107938 
    ## Run 357 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218063  max resid 0.07687089 
    ## Run 358 stress 0.1051403 
    ## ... Procrustes: rmse 0.02180338  max resid 0.0768622 
    ## Run 359 stress 0.1135596 
    ## Run 360 stress 0.1063128 
    ## Run 361 stress 0.1058423 
    ## Run 362 stress 0.116055 
    ## Run 363 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180292  max resid 0.07686195 
    ## Run 364 stress 0.1160546 
    ## Run 365 stress 0.1050144 
    ## ... Procrustes: rmse 0.009669367  max resid 0.02607736 
    ## Run 366 stress 0.1063128 
    ## Run 367 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283775  max resid 0.07728465 
    ## Run 368 stress 0.1097987 
    ## Run 369 stress 0.1135597 
    ## Run 370 stress 0.1058708 
    ## Run 371 stress 0.1138392 
    ## Run 372 stress 0.1063128 
    ## Run 373 stress 0.113594 
    ## Run 374 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284251  max resid 0.07725027 
    ## Run 375 stress 0.1088547 
    ## Run 376 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180348  max resid 0.07686157 
    ## Run 377 stress 0.1079382 
    ## Run 378 stress 0.1083501 
    ## Run 379 stress 0.1084965 
    ## Run 380 stress 0.1049903 
    ## ... Procrustes: rmse 0.006727322  max resid 0.02579854 
    ## Run 381 stress 0.1136646 
    ## Run 382 stress 0.1050338 
    ## ... Procrustes: rmse 0.006257203  max resid 0.02187428 
    ## Run 383 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180289  max resid 0.0768594 
    ## Run 384 stress 0.1088547 
    ## Run 385 stress 0.1063128 
    ## Run 386 stress 0.1079383 
    ## Run 387 stress 0.1102723 
    ## Run 388 stress 0.1063128 
    ## Run 389 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283893  max resid 0.07730396 
    ## Run 390 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284181  max resid 0.07730786 
    ## Run 391 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283923  max resid 0.07728822 
    ## Run 392 stress 0.1083994 
    ## Run 393 stress 0.1063128 
    ## Run 394 stress 0.1050144 
    ## ... Procrustes: rmse 0.009675381  max resid 0.02617191 
    ## Run 395 stress 0.1136646 
    ## Run 396 stress 0.1050144 
    ## ... Procrustes: rmse 0.009663134  max resid 0.02609498 
    ## Run 397 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180392  max resid 0.07686482 
    ## Run 398 stress 0.1050144 
    ## ... Procrustes: rmse 0.009658012  max resid 0.02606388 
    ## Run 399 stress 0.1063128 
    ## Run 400 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180144  max resid 0.0768395 
    ## Run 401 stress 0.1160538 
    ## Run 402 stress 0.1051402 
    ## ... Procrustes: rmse 0.02178877  max resid 0.07681068 
    ## Run 403 stress 0.1135598 
    ## Run 404 stress 0.1049903 
    ## ... Procrustes: rmse 0.006736351  max resid 0.02583873 
    ## Run 405 stress 0.1160543 
    ## Run 406 stress 0.113594 
    ## Run 407 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283935  max resid 0.07730014 
    ## Run 408 stress 0.1083997 
    ## Run 409 stress 0.1088548 
    ## Run 410 stress 0.1174548 
    ## Run 411 stress 0.1084966 
    ## Run 412 stress 0.109799 
    ## Run 413 stress 0.1136647 
    ## Run 414 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180211  max resid 0.07684351 
    ## Run 415 stress 0.1049903 
    ## ... Procrustes: rmse 0.00671077  max resid 0.0257183 
    ## Run 416 stress 0.113594 
    ## Run 417 stress 0.107938 
    ## Run 418 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284261  max resid 0.07727199 
    ## Run 419 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180786  max resid 0.07687553 
    ## Run 420 stress 0.1160539 
    ## Run 421 stress 0.1063128 
    ## Run 422 stress 0.1136646 
    ## Run 423 stress 0.1050338 
    ## ... Procrustes: rmse 0.006259759  max resid 0.02189196 
    ## Run 424 stress 0.1063128 
    ## Run 425 stress 0.1088547 
    ## Run 426 stress 0.1083499 
    ## Run 427 stress 0.1058423 
    ## Run 428 stress 0.1088548 
    ## Run 429 stress 0.1079381 
    ## Run 430 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180266  max resid 0.0768583 
    ## Run 431 stress 0.1063128 
    ## Run 432 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180323  max resid 0.07686709 
    ## Run 433 stress 0.1063128 
    ## Run 434 stress 0.107938 
    ## Run 435 stress 0.1084964 
    ## Run 436 stress 0.10835 
    ## Run 437 stress 0.1088547 
    ## Run 438 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181063  max resid 0.07688639 
    ## Run 439 stress 0.1097987 
    ## Run 440 stress 0.1063128 
    ## Run 441 stress 0.107938 
    ## Run 442 stress 0.1083501 
    ## Run 443 stress 0.1063129 
    ## Run 444 stress 0.1049791 
    ## ... Procrustes: rmse 4.720104e-06  max resid 1.388798e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.113594 
    ## Run 446 stress 0.1160551 
    ## Run 447 stress 0.1049903 
    ## ... Procrustes: rmse 0.006744727  max resid 0.02587563 
    ## Run 448 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180458  max resid 0.07686443 
    ## Run 449 stress 0.1084965 
    ## Run 450 stress 0.1050144 
    ## ... Procrustes: rmse 0.009657749  max resid 0.02605725 
    ## Run 451 stress 0.1160544 
    ## Run 452 stress 0.1049791 
    ## ... Procrustes: rmse 5.446421e-06  max resid 1.661182e-05 
    ## ... Similar to previous best
    ## Run 453 stress 0.1138398 
    ## Run 454 stress 0.1088548 
    ## Run 455 stress 0.1079382 
    ## Run 456 stress 0.1049903 
    ## ... Procrustes: rmse 0.006739585  max resid 0.02585556 
    ## Run 457 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180323  max resid 0.07686073 
    ## Run 458 stress 0.1082722 
    ## Run 459 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180264  max resid 0.07685115 
    ## Run 460 stress 0.1082722 
    ## Run 461 stress 0.1082722 
    ## Run 462 stress 0.10835 
    ## Run 463 stress 0.1084966 
    ## Run 464 stress 0.1135598 
    ## Run 465 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284092  max resid 0.07734428 
    ## Run 466 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283883  max resid 0.07727517 
    ## Run 467 stress 0.1083994 
    ## Run 468 stress 0.1160545 
    ## Run 469 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218009  max resid 0.07682971 
    ## Run 470 stress 0.1050338 
    ## ... Procrustes: rmse 0.006254744  max resid 0.02186264 
    ## Run 471 stress 0.1160542 
    ## Run 472 stress 0.1136647 
    ## Run 473 stress 0.10835 
    ## Run 474 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283946  max resid 0.07729949 
    ## Run 475 stress 0.1049903 
    ## ... Procrustes: rmse 0.006712782  max resid 0.02571947 
    ## Run 476 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180285  max resid 0.07686086 
    ## Run 477 stress 0.1084964 
    ## Run 478 stress 0.1088547 
    ## Run 479 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181296  max resid 0.07689228 
    ## Run 480 stress 0.1135599 
    ## Run 481 stress 0.1058423 
    ## Run 482 stress 0.1050144 
    ## ... Procrustes: rmse 0.009683232  max resid 0.02620066 
    ## Run 483 stress 0.1084965 
    ## Run 484 stress 0.1138392 
    ## Run 485 stress 0.1138394 
    ## Run 486 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180372  max resid 0.07686203 
    ## Run 487 stress 0.113594 
    ## Run 488 stress 0.1160548 
    ## Run 489 stress 0.1083996 
    ## Run 490 stress 0.1088547 
    ## Run 491 stress 0.1138391 
    ## Run 492 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283467  max resid 0.07728282 
    ## Run 493 stress 0.10835 
    ## Run 494 stress 0.1079381 
    ## Run 495 stress 0.1097988 
    ## Run 496 stress 0.1160545 
    ## Run 497 stress 0.1049903 
    ## ... Procrustes: rmse 0.006777043  max resid 0.0260287 
    ## Run 498 stress 0.1136646 
    ## Run 499 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180792  max resid 0.07687598 
    ## Run 500 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283979  max resid 0.07730427 
    ## *** Best solution repeated 17 times

``` r
# Surveyed sites 
SD_beta_NMDS <- metaMDS(SD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.1056896 
    ## Run 2 stress 0.08946651 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0255846  max resid 0.0777671 
    ## Run 3 stress 0.1052647 
    ## Run 4 stress 0.0903913 
    ## Run 5 stress 0.0894633 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01224812  max resid 0.03851522 
    ## Run 6 stress 0.09039149 
    ## Run 7 stress 0.1128766 
    ## Run 8 stress 0.1111016 
    ## Run 9 stress 0.08938554 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03576317  max resid 0.1192182 
    ## Run 10 stress 0.08938538 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002101881  max resid 0.000569513 
    ## ... Similar to previous best
    ## Run 11 stress 0.08926082 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01181321  max resid 0.04004315 
    ## Run 12 stress 0.09262344 
    ## Run 13 stress 0.08926082 
    ## ... New best solution
    ## ... Procrustes: rmse 3.693238e-05  max resid 8.360272e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.1092211 
    ## Run 15 stress 0.08951717 
    ## ... Procrustes: rmse 0.03511727  max resid 0.1157079 
    ## Run 16 stress 0.08926082 
    ## ... Procrustes: rmse 2.110522e-05  max resid 4.710193e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.1071321 
    ## Run 18 stress 0.08951738 
    ## ... Procrustes: rmse 0.03506199  max resid 0.1156456 
    ## Run 19 stress 0.1071321 
    ## Run 20 stress 0.08938538 
    ## ... Procrustes: rmse 0.01182108  max resid 0.03980232 
    ## Run 21 stress 0.08938538 
    ## ... Procrustes: rmse 0.01184006  max resid 0.03980784 
    ## Run 22 stress 0.1091827 
    ## Run 23 stress 0.09503445 
    ## Run 24 stress 0.09039141 
    ## Run 25 stress 0.109213 
    ## Run 26 stress 0.09775535 
    ## Run 27 stress 0.109188 
    ## Run 28 stress 0.09503424 
    ## Run 29 stress 0.08951723 
    ## ... Procrustes: rmse 0.03508763  max resid 0.1156756 
    ## Run 30 stress 0.08946651 
    ## ... Procrustes: rmse 0.0332735  max resid 0.1164182 
    ## Run 31 stress 0.1067893 
    ## Run 32 stress 0.08926095 
    ## ... Procrustes: rmse 9.082001e-05  max resid 0.0002648038 
    ## ... Similar to previous best
    ## Run 33 stress 0.1087578 
    ## Run 34 stress 0.08938548 
    ## ... Procrustes: rmse 0.01178243  max resid 0.03983955 
    ## Run 35 stress 0.08946328 
    ## ... Procrustes: rmse 0.03756339  max resid 0.1178354 
    ## Run 36 stress 0.090886 
    ## Run 37 stress 0.1052653 
    ## Run 38 stress 0.1104181 
    ## Run 39 stress 0.08938961 
    ## ... Procrustes: rmse 0.0360013  max resid 0.1183586 
    ## Run 40 stress 0.08926076 
    ## ... New best solution
    ## ... Procrustes: rmse 5.155112e-05  max resid 0.0001174567 
    ## ... Similar to previous best
    ## Run 41 stress 0.1052647 
    ## Run 42 stress 0.08938538 
    ## ... Procrustes: rmse 0.01182114  max resid 0.03975061 
    ## Run 43 stress 0.08938969 
    ## ... Procrustes: rmse 0.03596071  max resid 0.1182965 
    ## Run 44 stress 0.09180895 
    ## Run 45 stress 0.08938547 
    ## ... Procrustes: rmse 0.0118368  max resid 0.03979455 
    ## Run 46 stress 0.1061298 
    ## Run 47 stress 0.1086562 
    ## Run 48 stress 0.08938538 
    ## ... Procrustes: rmse 0.01180244  max resid 0.03974087 
    ## Run 49 stress 0.09021165 
    ## Run 50 stress 0.1075422 
    ## Run 51 stress 0.09088601 
    ## Run 52 stress 0.1089907 
    ## Run 53 stress 0.1099391 
    ## Run 54 stress 0.1056906 
    ## Run 55 stress 0.08938546 
    ## ... Procrustes: rmse 0.01186772  max resid 0.03974686 
    ## Run 56 stress 0.1091448 
    ## Run 57 stress 0.09088605 
    ## Run 58 stress 0.08946657 
    ## ... Procrustes: rmse 0.03324462  max resid 0.116371 
    ## Run 59 stress 0.107902 
    ## Run 60 stress 0.09109127 
    ## Run 61 stress 0.09099542 
    ## Run 62 stress 0.08938547 
    ## ... Procrustes: rmse 0.01175156  max resid 0.03974204 
    ## Run 63 stress 0.1103726 
    ## Run 64 stress 0.1087556 
    ## Run 65 stress 0.09592141 
    ## Run 66 stress 0.1092362 
    ## Run 67 stress 0.09088613 
    ## Run 68 stress 0.1092096 
    ## Run 69 stress 0.1064193 
    ## Run 70 stress 0.09228494 
    ## Run 71 stress 0.09262349 
    ## Run 72 stress 0.1076305 
    ## Run 73 stress 0.1071321 
    ## Run 74 stress 0.09039143 
    ## Run 75 stress 0.08946331 
    ## ... Procrustes: rmse 0.0375355  max resid 0.1178038 
    ## Run 76 stress 0.1056901 
    ## Run 77 stress 0.1071322 
    ## Run 78 stress 0.08926106 
    ## ... Procrustes: rmse 0.0001871246  max resid 0.0005694961 
    ## ... Similar to previous best
    ## Run 79 stress 0.09044628 
    ## Run 80 stress 0.1093313 
    ## Run 81 stress 0.105265 
    ## Run 82 stress 0.1052649 
    ## Run 83 stress 0.1075422 
    ## Run 84 stress 0.109213 
    ## Run 85 stress 0.08938961 
    ## ... Procrustes: rmse 0.03599655  max resid 0.1183537 
    ## Run 86 stress 0.08938981 
    ## ... Procrustes: rmse 0.03594243  max resid 0.1182734 
    ## Run 87 stress 0.08938962 
    ## ... Procrustes: rmse 0.03599715  max resid 0.1183471 
    ## Run 88 stress 0.09039129 
    ## Run 89 stress 0.1056904 
    ## Run 90 stress 0.09099526 
    ## Run 91 stress 0.08951719 
    ## ... Procrustes: rmse 0.03511936  max resid 0.1157122 
    ## Run 92 stress 0.0893854 
    ## ... Procrustes: rmse 0.01181279  max resid 0.03975168 
    ## Run 93 stress 0.08938557 
    ## ... Procrustes: rmse 0.01186633  max resid 0.03967331 
    ## Run 94 stress 0.08926082 
    ## ... Procrustes: rmse 4.991033e-05  max resid 0.0001671994 
    ## ... Similar to previous best
    ## Run 95 stress 0.09130087 
    ## Run 96 stress 0.1074433 
    ## Run 97 stress 0.106257 
    ## Run 98 stress 0.1092093 
    ## Run 99 stress 0.1079016 
    ## Run 100 stress 0.08926064 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002601774  max resid 0.0008196697 
    ## ... Similar to previous best
    ## Run 101 stress 0.09099545 
    ## Run 102 stress 0.1076302 
    ## Run 103 stress 0.0977554 
    ## Run 104 stress 0.08926071 
    ## ... Procrustes: rmse 0.0002073209  max resid 0.0006029231 
    ## ... Similar to previous best
    ## Run 105 stress 0.09109108 
    ## Run 106 stress 0.09228489 
    ## Run 107 stress 0.09503429 
    ## Run 108 stress 0.08946669 
    ## ... Procrustes: rmse 0.03316956  max resid 0.1162626 
    ## Run 109 stress 0.110373 
    ## Run 110 stress 0.1063194 
    ## Run 111 stress 0.1091877 
    ## Run 112 stress 0.08938542 
    ## ... Procrustes: rmse 0.01186723  max resid 0.03963609 
    ## Run 113 stress 0.09503429 
    ## Run 114 stress 0.1101455 
    ## Run 115 stress 0.08938977 
    ## ... Procrustes: rmse 0.03590646  max resid 0.1182185 
    ## Run 116 stress 0.08946335 
    ## ... Procrustes: rmse 0.03748228  max resid 0.1177235 
    ## Run 117 stress 0.109213 
    ## Run 118 stress 0.08926083 
    ## ... Procrustes: rmse 0.0002699907  max resid 0.0007552884 
    ## ... Similar to previous best
    ## Run 119 stress 0.1091239 
    ## Run 120 stress 0.1075421 
    ## Run 121 stress 0.1091829 
    ## Run 122 stress 0.09109109 
    ## Run 123 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 4.517706e-05  max resid 0.0001127842 
    ## ... Similar to previous best
    ## Run 124 stress 0.0894666 
    ## ... Procrustes: rmse 0.03319062  max resid 0.1162954 
    ## Run 125 stress 0.1089906 
    ## Run 126 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595859  max resid 0.1182955 
    ## Run 127 stress 0.0893854 
    ## ... Procrustes: rmse 0.01181941  max resid 0.03969909 
    ## Run 128 stress 0.09503416 
    ## Run 129 stress 0.1071321 
    ## Run 130 stress 0.09592164 
    ## Run 131 stress 0.105265 
    ## Run 132 stress 0.08938547 
    ## ... Procrustes: rmse 0.01191466  max resid 0.03970776 
    ## Run 133 stress 0.09503461 
    ## Run 134 stress 0.0903913 
    ## Run 135 stress 0.1063199 
    ## Run 136 stress 0.09503447 
    ## Run 137 stress 0.09109108 
    ## Run 138 stress 0.1065223 
    ## Run 139 stress 0.0908736 
    ## Run 140 stress 0.08946661 
    ## ... Procrustes: rmse 0.03318921  max resid 0.1162895 
    ## Run 141 stress 0.0894634 
    ## ... Procrustes: rmse 0.03748057  max resid 0.1177095 
    ## Run 142 stress 0.1052648 
    ## Run 143 stress 0.08926084 
    ## ... Procrustes: rmse 0.0002765644  max resid 0.0007876147 
    ## ... Similar to previous best
    ## Run 144 stress 0.08938551 
    ## ... Procrustes: rmse 0.0118473  max resid 0.0397788 
    ## Run 145 stress 0.08938969 
    ## ... Procrustes: rmse 0.03592159  max resid 0.1182366 
    ## Run 146 stress 0.1062561 
    ## Run 147 stress 0.09087339 
    ## Run 148 stress 0.09021166 
    ## Run 149 stress 0.09039134 
    ## Run 150 stress 0.1056899 
    ## Run 151 stress 0.09503429 
    ## Run 152 stress 0.0892609 
    ## ... Procrustes: rmse 0.0003374285  max resid 0.001123617 
    ## ... Similar to previous best
    ## Run 153 stress 0.08926115 
    ## ... Procrustes: rmse 0.000468608  max resid 0.001316792 
    ## ... Similar to previous best
    ## Run 154 stress 0.1091449 
    ## Run 155 stress 0.1074434 
    ## Run 156 stress 0.08938962 
    ## ... Procrustes: rmse 0.03595106  max resid 0.1182811 
    ## Run 157 stress 0.09039144 
    ## Run 158 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003435094  max resid 0.001142035 
    ## ... Similar to previous best
    ## Run 159 stress 0.1071322 
    ## Run 160 stress 0.1064188 
    ## Run 161 stress 0.0894665 
    ## ... Procrustes: rmse 0.03322067  max resid 0.1163409 
    ## Run 162 stress 0.1075421 
    ## Run 163 stress 0.1092092 
    ## Run 164 stress 0.0903913 
    ## Run 165 stress 0.09178305 
    ## Run 166 stress 0.08946665 
    ## ... Procrustes: rmse 0.03318234  max resid 0.1162821 
    ## Run 167 stress 0.1092092 
    ## Run 168 stress 0.08938964 
    ## ... Procrustes: rmse 0.03593447  max resid 0.1182579 
    ## Run 169 stress 0.09087344 
    ## Run 170 stress 0.08938558 
    ## ... Procrustes: rmse 0.01179437  max resid 0.03978054 
    ## Run 171 stress 0.1075423 
    ## Run 172 stress 0.1076306 
    ## Run 173 stress 0.09039129 
    ## Run 174 stress 0.08926067 
    ## ... Procrustes: rmse 0.0001334083  max resid 0.0003773907 
    ## ... Similar to previous best
    ## Run 175 stress 0.08926105 
    ## ... Procrustes: rmse 0.0004247327  max resid 0.001406962 
    ## ... Similar to previous best
    ## Run 176 stress 0.08938963 
    ## ... Procrustes: rmse 0.03596834  max resid 0.1183092 
    ## Run 177 stress 0.1092092 
    ## Run 178 stress 0.1087866 
    ## Run 179 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003114899  max resid 0.001031788 
    ## ... Similar to previous best
    ## Run 180 stress 0.1071321 
    ## Run 181 stress 0.08926069 
    ## ... Procrustes: rmse 0.0001658658  max resid 0.0004852333 
    ## ... Similar to previous best
    ## Run 182 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595704  max resid 0.1182917 
    ## Run 183 stress 0.09503444 
    ## Run 184 stress 0.08926081 
    ## ... Procrustes: rmse 0.0002790409  max resid 0.0009200021 
    ## ... Similar to previous best
    ## Run 185 stress 0.09178297 
    ## Run 186 stress 0.1075422 
    ## Run 187 stress 0.1091882 
    ## Run 188 stress 0.08926098 
    ## ... Procrustes: rmse 0.0003868935  max resid 0.001294153 
    ## ... Similar to previous best
    ## Run 189 stress 0.09087355 
    ## Run 190 stress 0.1121881 
    ## Run 191 stress 0.1074432 
    ## Run 192 stress 0.08951723 
    ## ... Procrustes: rmse 0.0350351  max resid 0.1156041 
    ## Run 193 stress 0.1092127 
    ## Run 194 stress 0.1071323 
    ## Run 195 stress 0.09503466 
    ## Run 196 stress 0.09039129 
    ## Run 197 stress 0.09088607 
    ## Run 198 stress 0.1110594 
    ## Run 199 stress 0.1067516 
    ## Run 200 stress 0.1075422 
    ## Run 201 stress 0.0894667 
    ## ... Procrustes: rmse 0.03318578  max resid 0.1162985 
    ## Run 202 stress 0.1067584 
    ## Run 203 stress 0.1067477 
    ## Run 204 stress 0.09039131 
    ## Run 205 stress 0.09503466 
    ## Run 206 stress 0.1101453 
    ## Run 207 stress 0.1087863 
    ## Run 208 stress 0.08938541 
    ## ... Procrustes: rmse 0.01187569  max resid 0.03977336 
    ## Run 209 stress 0.08938968 
    ## ... Procrustes: rmse 0.03597701  max resid 0.1183359 
    ## Run 210 stress 0.08946653 
    ## ... Procrustes: rmse 0.03320472  max resid 0.1163169 
    ## Run 211 stress 0.1092096 
    ## Run 212 stress 0.08946654 
    ## ... Procrustes: rmse 0.03319994  max resid 0.116307 
    ## Run 213 stress 0.09039129 
    ## Run 214 stress 0.1105341 
    ## Run 215 stress 0.1092128 
    ## Run 216 stress 0.08946652 
    ## ... Procrustes: rmse 0.03323212  max resid 0.1163699 
    ## Run 217 stress 0.1074434 
    ## Run 218 stress 0.1062559 
    ## Run 219 stress 0.1079018 
    ## Run 220 stress 0.09044615 
    ## Run 221 stress 0.0909957 
    ## Run 222 stress 0.1086568 
    ## Run 223 stress 0.1071322 
    ## Run 224 stress 0.1071322 
    ## Run 225 stress 0.08946662 
    ## ... Procrustes: rmse 0.03318825  max resid 0.1162878 
    ## Run 226 stress 0.1079015 
    ## Run 227 stress 0.1079017 
    ## Run 228 stress 0.09087341 
    ## Run 229 stress 0.1060429 
    ## Run 230 stress 0.1072678 
    ## Run 231 stress 0.08926101 
    ## ... Procrustes: rmse 0.0004101908  max resid 0.001158032 
    ## ... Similar to previous best
    ## Run 232 stress 0.1092131 
    ## Run 233 stress 0.08938965 
    ## ... Procrustes: rmse 0.03593722  max resid 0.1182643 
    ## Run 234 stress 0.1076304 
    ## Run 235 stress 0.08938559 
    ## ... Procrustes: rmse 0.01175544  max resid 0.03961936 
    ## Run 236 stress 0.08938555 
    ## ... Procrustes: rmse 0.01175051  max resid 0.03970397 
    ## Run 237 stress 0.08938979 
    ## ... Procrustes: rmse 0.03591337  max resid 0.1182294 
    ## Run 238 stress 0.08926078 
    ## ... Procrustes: rmse 0.0002649005  max resid 0.0007433967 
    ## ... Similar to previous best
    ## Run 239 stress 0.08926099 
    ## ... Procrustes: rmse 0.0004022045  max resid 0.001131037 
    ## ... Similar to previous best
    ## Run 240 stress 0.0908861 
    ## Run 241 stress 0.1056902 
    ## Run 242 stress 0.1062598 
    ## Run 243 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 5.529206e-05  max resid 0.0001453224 
    ## ... Similar to previous best
    ## Run 244 stress 0.08938538 
    ## ... Procrustes: rmse 0.01180372  max resid 0.03963239 
    ## Run 245 stress 0.09130096 
    ## Run 246 stress 0.08938962 
    ## ... Procrustes: rmse 0.03595468  max resid 0.1182992 
    ## Run 247 stress 0.109213 
    ## Run 248 stress 0.1052649 
    ## Run 249 stress 0.1096133 
    ## Run 250 stress 0.09087348 
    ## Run 251 stress 0.08938541 
    ## ... Procrustes: rmse 0.01183826  max resid 0.03967018 
    ## Run 252 stress 0.1110408 
    ## Run 253 stress 0.1105345 
    ## Run 254 stress 0.09099602 
    ## Run 255 stress 0.1052652 
    ## Run 256 stress 0.1063194 
    ## Run 257 stress 0.09018712 
    ## Run 258 stress 0.1074433 
    ## Run 259 stress 0.1052653 
    ## Run 260 stress 0.08926102 
    ## ... Procrustes: rmse 0.0004240408  max resid 0.001348319 
    ## ... Similar to previous best
    ## Run 261 stress 0.09503443 
    ## Run 262 stress 0.09503451 
    ## Run 263 stress 0.0909956 
    ## Run 264 stress 0.08946328 
    ## ... Procrustes: rmse 0.03750809  max resid 0.1177616 
    ## Run 265 stress 0.09018707 
    ## Run 266 stress 0.08926064 
    ## ... Procrustes: rmse 7.086604e-05  max resid 0.0002399916 
    ## ... Similar to previous best
    ## Run 267 stress 0.1079017 
    ## Run 268 stress 0.1103723 
    ## Run 269 stress 0.089261 
    ## ... Procrustes: rmse 0.0004120524  max resid 0.00131672 
    ## ... Similar to previous best
    ## Run 270 stress 0.09503435 
    ## Run 271 stress 0.09018708 
    ## Run 272 stress 0.09503414 
    ## Run 273 stress 0.09503416 
    ## Run 274 stress 0.09503414 
    ## Run 275 stress 0.08926071 
    ## ... Procrustes: rmse 0.0002018797  max resid 0.000604547 
    ## ... Similar to previous best
    ## Run 276 stress 0.08926086 
    ## ... Procrustes: rmse 0.0003339247  max resid 0.0009853941 
    ## ... Similar to previous best
    ## Run 277 stress 0.09228484 
    ## Run 278 stress 0.1074432 
    ## Run 279 stress 0.09130096 
    ## Run 280 stress 0.09088606 
    ## Run 281 stress 0.08946653 
    ## ... Procrustes: rmse 0.03322738  max resid 0.1163814 
    ## Run 282 stress 0.08938974 
    ## ... Procrustes: rmse 0.03591285  max resid 0.1182337 
    ## Run 283 stress 0.08926082 
    ## ... Procrustes: rmse 0.0002713678  max resid 0.0009192772 
    ## ... Similar to previous best
    ## Run 284 stress 0.1061308 
    ## Run 285 stress 0.09018707 
    ## Run 286 stress 0.1092127 
    ## Run 287 stress 0.0908734 
    ## Run 288 stress 0.1087866 
    ## Run 289 stress 0.0901871 
    ## Run 290 stress 0.1056907 
    ## Run 291 stress 0.08938539 
    ## ... Procrustes: rmse 0.011836  max resid 0.03967173 
    ## Run 292 stress 0.08946653 
    ## ... Procrustes: rmse 0.0331959  max resid 0.1163143 
    ## Run 293 stress 0.09088603 
    ## Run 294 stress 0.09503418 
    ## Run 295 stress 0.08938557 
    ## ... Procrustes: rmse 0.01177074  max resid 0.03970469 
    ## Run 296 stress 0.1086563 
    ## Run 297 stress 0.08938538 
    ## ... Procrustes: rmse 0.0118393  max resid 0.03963038 
    ## Run 298 stress 0.09087342 
    ## Run 299 stress 0.09130085 
    ## Run 300 stress 0.0901871 
    ## Run 301 stress 0.0894665 
    ## ... Procrustes: rmse 0.03320977  max resid 0.1163354 
    ## Run 302 stress 0.08938961 
    ## ... Procrustes: rmse 0.03593881  max resid 0.1182703 
    ## Run 303 stress 0.08926073 
    ## ... Procrustes: rmse 0.0002255617  max resid 0.0006700424 
    ## ... Similar to previous best
    ## Run 304 stress 0.1089908 
    ## Run 305 stress 0.08946652 
    ## ... Procrustes: rmse 0.03319882  max resid 0.1163131 
    ## Run 306 stress 0.1089903 
    ## Run 307 stress 0.1074432 
    ## Run 308 stress 0.08938552 
    ## ... Procrustes: rmse 0.01188063  max resid 0.03958886 
    ## Run 309 stress 0.1071322 
    ## Run 310 stress 0.1074432 
    ## Run 311 stress 0.1079016 
    ## Run 312 stress 0.08926093 
    ## ... Procrustes: rmse 0.0003350983  max resid 0.001023947 
    ## ... Similar to previous best
    ## Run 313 stress 0.08926063 
    ## ... Procrustes: rmse 6.94574e-05  max resid 0.0001636776 
    ## ... Similar to previous best
    ## Run 314 stress 0.08946337 
    ## ... Procrustes: rmse 0.03747316  max resid 0.1177132 
    ## Run 315 stress 0.09228508 
    ## Run 316 stress 0.08938966 
    ## ... Procrustes: rmse 0.03596663  max resid 0.1183221 
    ## Run 317 stress 0.09044612 
    ## Run 318 stress 0.08938541 
    ## ... Procrustes: rmse 0.01185877  max resid 0.03967639 
    ## Run 319 stress 0.106256 
    ## Run 320 stress 0.08938968 
    ## ... Procrustes: rmse 0.03592001  max resid 0.1182347 
    ## Run 321 stress 0.1052649 
    ## Run 322 stress 0.09503413 
    ## Run 323 stress 0.1052652 
    ## Run 324 stress 0.09503434 
    ## Run 325 stress 0.1062567 
    ## Run 326 stress 0.08946333 
    ## ... Procrustes: rmse 0.03747963  max resid 0.1177262 
    ## Run 327 stress 0.09503427 
    ## Run 328 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 3.812119e-05  max resid 0.000100502 
    ## ... Similar to previous best
    ## Run 329 stress 0.08938553 
    ## ... Procrustes: rmse 0.01184515  max resid 0.03958961 
    ## Run 330 stress 0.08946652 
    ## ... Procrustes: rmse 0.03321009  max resid 0.1163275 
    ## Run 331 stress 0.09503428 
    ## Run 332 stress 0.09503424 
    ## Run 333 stress 0.09088627 
    ## Run 334 stress 0.0950342 
    ## Run 335 stress 0.1103697 
    ## Run 336 stress 0.1075421 
    ## Run 337 stress 0.08926077 
    ## ... Procrustes: rmse 0.0002329213  max resid 0.000717916 
    ## ... Similar to previous best
    ## Run 338 stress 0.106758 
    ## Run 339 stress 0.08926065 
    ## ... Procrustes: rmse 9.271713e-05  max resid 0.0002912259 
    ## ... Similar to previous best
    ## Run 340 stress 0.109145 
    ## Run 341 stress 0.1065384 
    ## Run 342 stress 0.09109108 
    ## Run 343 stress 0.08938537 
    ## ... Procrustes: rmse 0.01184333  max resid 0.03969203 
    ## Run 344 stress 0.08926089 
    ## ... Procrustes: rmse 0.0003303319  max resid 0.001022852 
    ## ... Similar to previous best
    ## Run 345 stress 0.08938555 
    ## ... Procrustes: rmse 0.01190812  max resid 0.0396396 
    ## Run 346 stress 0.08946331 
    ## ... Procrustes: rmse 0.03749737  max resid 0.1177382 
    ## Run 347 stress 0.1052648 
    ## Run 348 stress 0.09039148 
    ## Run 349 stress 0.1074433 
    ## Run 350 stress 0.08938973 
    ## ... Procrustes: rmse 0.03598627  max resid 0.1183363 
    ## Run 351 stress 0.1056898 
    ## Run 352 stress 0.08946329 
    ## ... Procrustes: rmse 0.03752809  max resid 0.1177821 
    ## Run 353 stress 0.09021168 
    ## Run 354 stress 0.08946335 
    ## ... Procrustes: rmse 0.03749015  max resid 0.1177272 
    ## Run 355 stress 0.08951716 
    ## ... Procrustes: rmse 0.03506158  max resid 0.1156399 
    ## Run 356 stress 0.09039129 
    ## Run 357 stress 0.1087864 
    ## Run 358 stress 0.1091238 
    ## Run 359 stress 0.1061298 
    ## Run 360 stress 0.09592153 
    ## Run 361 stress 0.1092126 
    ## Run 362 stress 0.1064426 
    ## Run 363 stress 0.09087358 
    ## Run 364 stress 0.08946328 
    ## ... Procrustes: rmse 0.03752654  max resid 0.1177833 
    ## Run 365 stress 0.08926103 
    ## ... Procrustes: rmse 0.0004107689  max resid 0.001306173 
    ## ... Similar to previous best
    ## Run 366 stress 0.08951737 
    ## ... Procrustes: rmse 0.03502279  max resid 0.1155969 
    ## Run 367 stress 0.1052651 
    ## Run 368 stress 0.1071321 
    ## Run 369 stress 0.1099389 
    ## Run 370 stress 0.09088605 
    ## Run 371 stress 0.1056892 
    ## Run 372 stress 0.08926098 
    ## ... Procrustes: rmse 0.0003261901  max resid 0.000936182 
    ## ... Similar to previous best
    ## Run 373 stress 0.0894665 
    ## ... Procrustes: rmse 0.03322358  max resid 0.1163491 
    ## Run 374 stress 0.09503442 
    ## Run 375 stress 0.09018712 
    ## Run 376 stress 0.1064188 
    ## Run 377 stress 0.1056894 
    ## Run 378 stress 0.1096134 
    ## Run 379 stress 0.09039133 
    ## Run 380 stress 0.0894634 
    ## ... Procrustes: rmse 0.03747937  max resid 0.1177083 
    ## Run 381 stress 0.09592143 
    ## Run 382 stress 0.1066192 
    ## Run 383 stress 0.1052651 
    ## Run 384 stress 0.1087572 
    ## Run 385 stress 0.09044608 
    ## Run 386 stress 0.09109109 
    ## Run 387 stress 0.1067885 
    ## Run 388 stress 0.08938542 
    ## ... Procrustes: rmse 0.01186477  max resid 0.03973622 
    ## Run 389 stress 0.1086563 
    ## Run 390 stress 0.1092126 
    ## Run 391 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003335173  max resid 0.001067117 
    ## ... Similar to previous best
    ## Run 392 stress 0.09039134 
    ## Run 393 stress 0.08938539 
    ## ... Procrustes: rmse 0.0118495  max resid 0.03972377 
    ## Run 394 stress 0.106861 
    ## Run 395 stress 0.1085128 
    ## Run 396 stress 0.1093314 
    ## Run 397 stress 0.1074433 
    ## Run 398 stress 0.1085128 
    ## Run 399 stress 0.0977554 
    ## Run 400 stress 0.1091502 
    ## Run 401 stress 0.09503413 
    ## Run 402 stress 0.09130086 
    ## Run 403 stress 0.1076303 
    ## Run 404 stress 0.1075421 
    ## Run 405 stress 0.1122893 
    ## Run 406 stress 0.08926102 
    ## ... Procrustes: rmse 0.0004005859  max resid 0.001295037 
    ## ... Similar to previous best
    ## Run 407 stress 0.1076305 
    ## Run 408 stress 0.08938961 
    ## ... Procrustes: rmse 0.0359627  max resid 0.1183022 
    ## Run 409 stress 0.08938962 
    ## ... Procrustes: rmse 0.03594363  max resid 0.1182713 
    ## Run 410 stress 0.08926072 
    ## ... Procrustes: rmse 0.000181471  max resid 0.00055957 
    ## ... Similar to previous best
    ## Run 411 stress 0.09039144 
    ## Run 412 stress 0.08946667 
    ## ... Procrustes: rmse 0.03318695  max resid 0.1162979 
    ## Run 413 stress 0.09039144 
    ## Run 414 stress 0.1080635 
    ## Run 415 stress 0.1071321 
    ## Run 416 stress 0.1052646 
    ## Run 417 stress 0.08938542 
    ## ... Procrustes: rmse 0.01183462  max resid 0.03969879 
    ## Run 418 stress 0.08926065 
    ## ... Procrustes: rmse 9.490957e-05  max resid 0.0002546857 
    ## ... Similar to previous best
    ## Run 419 stress 0.08938964 
    ## ... Procrustes: rmse 0.03596789  max resid 0.1183096 
    ## Run 420 stress 0.1063202 
    ## Run 421 stress 0.1092093 
    ## Run 422 stress 0.09775574 
    ## Run 423 stress 0.1105345 
    ## Run 424 stress 0.09592168 
    ## Run 425 stress 0.09775564 
    ## Run 426 stress 0.08946652 
    ## ... Procrustes: rmse 0.03323269  max resid 0.1163723 
    ## Run 427 stress 0.09021167 
    ## Run 428 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595237  max resid 0.1182809 
    ## Run 429 stress 0.09044612 
    ## Run 430 stress 0.09088607 
    ## Run 431 stress 0.1067511 
    ## Run 432 stress 0.09039142 
    ## Run 433 stress 0.08926086 
    ## ... Procrustes: rmse 0.0002844738  max resid 0.0008982199 
    ## ... Similar to previous best
    ## Run 434 stress 0.09503417 
    ## Run 435 stress 0.08938966 
    ## ... Procrustes: rmse 0.03592991  max resid 0.1182475 
    ## Run 436 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595728  max resid 0.1182906 
    ## Run 437 stress 0.10569 
    ## Run 438 stress 0.08938983 
    ## ... Procrustes: rmse 0.03590487  max resid 0.118214 
    ## Run 439 stress 0.106537 
    ## Run 440 stress 0.08938962 
    ## ... Procrustes: rmse 0.035964  max resid 0.1182979 
    ## Run 441 stress 0.08938977 
    ## ... Procrustes: rmse 0.0359141  max resid 0.1182255 
    ## Run 442 stress 0.0908735 
    ## Run 443 stress 0.1074431 
    ## Run 444 stress 0.1076303 
    ## Run 445 stress 0.0908735 
    ## Run 446 stress 0.1087554 
    ## Run 447 stress 0.1076308 
    ## Run 448 stress 0.08926093 
    ## ... Procrustes: rmse 0.0003474733  max resid 0.001110621 
    ## ... Similar to previous best
    ## Run 449 stress 0.08926065 
    ## ... Procrustes: rmse 8.903322e-05  max resid 0.0002667182 
    ## ... Similar to previous best
    ## Run 450 stress 0.09130105 
    ## Run 451 stress 0.08938967 
    ## ... Procrustes: rmse 0.03592739  max resid 0.1182463 
    ## Run 452 stress 0.0903913 
    ## Run 453 stress 0.08946653 
    ## ... Procrustes: rmse 0.03323641  max resid 0.1163924 
    ## Run 454 stress 0.09503423 
    ## Run 455 stress 0.08926076 
    ## ... Procrustes: rmse 0.0002284908  max resid 0.0007023854 
    ## ... Similar to previous best
    ## Run 456 stress 0.09503416 
    ## Run 457 stress 0.1061297 
    ## Run 458 stress 0.1060091 
    ## Run 459 stress 0.1060434 
    ## Run 460 stress 0.1066198 
    ## Run 461 stress 0.08926102 
    ## ... Procrustes: rmse 0.0004020654  max resid 0.001304023 
    ## ... Similar to previous best
    ## Run 462 stress 0.08938971 
    ## ... Procrustes: rmse 0.03592007  max resid 0.1182337 
    ## Run 463 stress 0.09228499 
    ## Run 464 stress 0.1052653 
    ## Run 465 stress 0.1067576 
    ## Run 466 stress 0.08938547 
    ## ... Procrustes: rmse 0.0118848  max resid 0.03965127 
    ## Run 467 stress 0.08926066 
    ## ... Procrustes: rmse 0.0001114664  max resid 0.0003355212 
    ## ... Similar to previous best
    ## Run 468 stress 0.09087361 
    ## Run 469 stress 0.1076305 
    ## Run 470 stress 0.08926064 
    ## ... Procrustes: rmse 8.109847e-05  max resid 0.0002105521 
    ## ... Similar to previous best
    ## Run 471 stress 0.1067507 
    ## Run 472 stress 0.08946654 
    ## ... Procrustes: rmse 0.03323855  max resid 0.116394 
    ## Run 473 stress 0.1071323 
    ## Run 474 stress 0.09503455 
    ## Run 475 stress 0.1079014 
    ## Run 476 stress 0.1108231 
    ## Run 477 stress 0.08926081 
    ## ... Procrustes: rmse 0.0002710364  max resid 0.0007754122 
    ## ... Similar to previous best
    ## Run 478 stress 0.1092363 
    ## Run 479 stress 0.08946332 
    ## ... Procrustes: rmse 0.03749549  max resid 0.1177338 
    ## Run 480 stress 0.1056908 
    ## Run 481 stress 0.0903913 
    ## Run 482 stress 0.08938538 
    ## ... Procrustes: rmse 0.01182998  max resid 0.03970307 
    ## Run 483 stress 0.08938962 
    ## ... Procrustes: rmse 0.03596505  max resid 0.118303 
    ## Run 484 stress 0.09044608 
    ## Run 485 stress 0.1108316 
    ## Run 486 stress 0.1067378 
    ## Run 487 stress 0.08938538 
    ## ... Procrustes: rmse 0.0118512  max resid 0.03970183 
    ## Run 488 stress 0.09087344 
    ## Run 489 stress 0.1096134 
    ## Run 490 stress 0.106419 
    ## Run 491 stress 0.1091827 
    ## Run 492 stress 0.1065387 
    ## Run 493 stress 0.08926076 
    ## ... Procrustes: rmse 0.0002085322  max resid 0.0005491421 
    ## ... Similar to previous best
    ## Run 494 stress 0.09775564 
    ## Run 495 stress 0.09109108 
    ## Run 496 stress 0.08926075 
    ## ... Procrustes: rmse 0.0002096744  max resid 0.0006893673 
    ## ... Similar to previous best
    ## Run 497 stress 0.09594317 
    ## Run 498 stress 0.09087344 
    ## Run 499 stress 0.09109115 
    ## Run 500 stress 0.0903913 
    ## *** Best solution repeated 20 times

``` r
round(SD_beta_NMDS$stress, digits = 2)
```

    ## [1] 0.09

``` r
SD_beta_rep_NMDS <- metaMDS(SD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.3323292 
    ## Run 1 stress 0.3306184 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1916454  max resid 0.3089509 
    ## Run 2 stress 0.3448005 
    ## Run 3 stress 0.3293859 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1693854  max resid 0.3210874 
    ## Run 4 stress 0.336057 
    ## Run 5 stress 0.3272826 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1718214  max resid 0.3740805 
    ## Run 6 stress 0.3396427 
    ## Run 7 stress 0.334652 
    ## Run 8 stress 0.3403966 
    ## Run 9 stress 0.3371535 
    ## Run 10 stress 0.337087 
    ## Run 11 stress 0.3336313 
    ## Run 12 stress 0.3399159 
    ## Run 13 stress 0.3466634 
    ## Run 14 stress 0.3399437 
    ## Run 15 stress 0.3488721 
    ## Run 16 stress 0.336652 
    ## Run 17 stress 0.3314219 
    ## Run 18 stress 0.3365103 
    ## Run 19 stress 0.336224 
    ## Run 20 stress 0.3543134 
    ## Run 21 stress 0.34004 
    ## Run 22 stress 0.3256439 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1732788  max resid 0.3655753 
    ## Run 23 stress 0.3375896 
    ## Run 24 stress 0.3482841 
    ## Run 25 stress 0.3383079 
    ## Run 26 stress 0.3483446 
    ## Run 27 stress 0.3300171 
    ## Run 28 stress 0.336239 
    ## Run 29 stress 0.3336752 
    ## Run 30 stress 0.3478879 
    ## Run 31 stress 0.3355325 
    ## Run 32 stress 0.3271787 
    ## Run 33 stress 0.3250324 
    ## ... New best solution
    ## ... Procrustes: rmse 0.09806118  max resid 0.299236 
    ## Run 34 stress 0.3296517 
    ## Run 35 stress 0.3508565 
    ## Run 36 stress 0.3417951 
    ## Run 37 stress 0.3295395 
    ## Run 38 stress 0.3306041 
    ## Run 39 stress 0.3345168 
    ## Run 40 stress 0.3452914 
    ## Run 41 stress 0.3382304 
    ## Run 42 stress 0.3340149 
    ## Run 43 stress 0.3363579 
    ## Run 44 stress 0.3336488 
    ## Run 45 stress 0.3378255 
    ## Run 46 stress 0.3343624 
    ## Run 47 stress 0.3451413 
    ## Run 48 stress 0.3407802 
    ## Run 49 stress 0.3344265 
    ## Run 50 stress 0.342994 
    ## Run 51 stress 0.3424141 
    ## Run 52 stress 0.3415506 
    ## Run 53 stress 0.3451086 
    ## Run 54 stress 0.3392158 
    ## Run 55 stress 0.3325979 
    ## Run 56 stress 0.3477402 
    ## Run 57 stress 0.3391733 
    ## Run 58 stress 0.3452808 
    ## Run 59 stress 0.3268028 
    ## Run 60 stress 0.3430295 
    ## Run 61 stress 0.336196 
    ## Run 62 stress 0.3496362 
    ## Run 63 stress 0.3338524 
    ## Run 64 stress 0.3683505 
    ## Run 65 stress 0.3330676 
    ## Run 66 stress 0.3314143 
    ## Run 67 stress 0.3310428 
    ## Run 68 stress 0.3359952 
    ## Run 69 stress 0.3371309 
    ## Run 70 stress 0.3391221 
    ## Run 71 stress 0.3379961 
    ## Run 72 stress 0.3390901 
    ## Run 73 stress 0.3404841 
    ## Run 74 stress 0.3415231 
    ## Run 75 stress 0.3328859 
    ## Run 76 stress 0.3375168 
    ## Run 77 stress 0.3372153 
    ## Run 78 stress 0.332224 
    ## Run 79 stress 0.3455194 
    ## Run 80 stress 0.3308744 
    ## Run 81 stress 0.3448311 
    ## Run 82 stress 0.3352453 
    ## Run 83 stress 0.3331709 
    ## Run 84 stress 0.3396623 
    ## Run 85 stress 0.3404924 
    ## Run 86 stress 0.3264052 
    ## Run 87 stress 0.3367305 
    ## Run 88 stress 0.3344131 
    ## Run 89 stress 0.3510925 
    ## Run 90 stress 0.3473632 
    ## Run 91 stress 0.3255598 
    ## Run 92 stress 0.338366 
    ## Run 93 stress 0.3366753 
    ## Run 94 stress 0.3335184 
    ## Run 95 stress 0.3391456 
    ## Run 96 stress 0.3253591 
    ## ... Procrustes: rmse 0.1615848  max resid 0.2545902 
    ## Run 97 stress 0.3313523 
    ## Run 98 stress 0.3456956 
    ## Run 99 stress 0.3293005 
    ## Run 100 stress 0.3299482 
    ## Run 101 stress 0.330111 
    ## Run 102 stress 0.3320738 
    ## Run 103 stress 0.3457612 
    ## Run 104 stress 0.330774 
    ## Run 105 stress 0.345568 
    ## Run 106 stress 0.3355662 
    ## Run 107 stress 0.3336955 
    ## Run 108 stress 0.3476024 
    ## Run 109 stress 0.330582 
    ## Run 110 stress 0.3461369 
    ## Run 111 stress 0.3334147 
    ## Run 112 stress 0.3336039 
    ## Run 113 stress 0.336643 
    ## Run 114 stress 0.3318333 
    ## Run 115 stress 0.3348642 
    ## Run 116 stress 0.3300022 
    ## Run 117 stress 0.3260298 
    ## Run 118 stress 0.339881 
    ## Run 119 stress 0.3398536 
    ## Run 120 stress 0.3446897 
    ## Run 121 stress 0.3349779 
    ## Run 122 stress 0.3465134 
    ## Run 123 stress 0.3274758 
    ## Run 124 stress 0.3343949 
    ## Run 125 stress 0.3343408 
    ## Run 126 stress 0.3322712 
    ## Run 127 stress 0.3378431 
    ## Run 128 stress 0.3503773 
    ## Run 129 stress 0.3443223 
    ## Run 130 stress 0.3452183 
    ## Run 131 stress 0.3452301 
    ## Run 132 stress 0.3289151 
    ## Run 133 stress 0.3345582 
    ## Run 134 stress 0.3454015 
    ## Run 135 stress 0.3358072 
    ## Run 136 stress 0.3369794 
    ## Run 137 stress 0.3346412 
    ## Run 138 stress 0.3257587 
    ## Run 139 stress 0.3424557 
    ## Run 140 stress 0.3362222 
    ## Run 141 stress 0.3327009 
    ## Run 142 stress 0.3393849 
    ## Run 143 stress 0.3358381 
    ## Run 144 stress 0.3408749 
    ## Run 145 stress 0.3482661 
    ## Run 146 stress 0.3373506 
    ## Run 147 stress 0.3357906 
    ## Run 148 stress 0.3364683 
    ## Run 149 stress 0.3312232 
    ## Run 150 stress 0.3387899 
    ## Run 151 stress 0.3468862 
    ## Run 152 stress 0.3328475 
    ## Run 153 stress 0.3300402 
    ## Run 154 stress 0.3300827 
    ## Run 155 stress 0.3293276 
    ## Run 156 stress 0.3372 
    ## Run 157 stress 0.3458751 
    ## Run 158 stress 0.3327084 
    ## Run 159 stress 0.3428003 
    ## Run 160 stress 0.3473125 
    ## Run 161 stress 0.3375503 
    ## Run 162 stress 0.3391988 
    ## Run 163 stress 0.3344368 
    ## Run 164 stress 0.3278636 
    ## Run 165 stress 0.333854 
    ## Run 166 stress 0.3361694 
    ## Run 167 stress 0.3280158 
    ## Run 168 stress 0.3412723 
    ## Run 169 stress 0.3436201 
    ## Run 170 stress 0.3352257 
    ## Run 171 stress 0.3332616 
    ## Run 172 stress 0.3426341 
    ## Run 173 stress 0.3415438 
    ## Run 174 stress 0.3424015 
    ## Run 175 stress 0.3365648 
    ## Run 176 stress 0.3434689 
    ## Run 177 stress 0.3444093 
    ## Run 178 stress 0.3293808 
    ## Run 179 stress 0.3297878 
    ## Run 180 stress 0.3293072 
    ## Run 181 stress 0.3397353 
    ## Run 182 stress 0.3458577 
    ## Run 183 stress 0.3306373 
    ## Run 184 stress 0.3305166 
    ## Run 185 stress 0.3451842 
    ## Run 186 stress 0.3343526 
    ## Run 187 stress 0.3367599 
    ## Run 188 stress 0.3306301 
    ## Run 189 stress 0.331652 
    ## Run 190 stress 0.346437 
    ## Run 191 stress 0.3447681 
    ## Run 192 stress 0.3425145 
    ## Run 193 stress 0.3350113 
    ## Run 194 stress 0.3330558 
    ## Run 195 stress 0.3362395 
    ## Run 196 stress 0.3380044 
    ## Run 197 stress 0.3431637 
    ## Run 198 stress 0.335574 
    ## Run 199 stress 0.3418125 
    ## Run 200 stress 0.3363125 
    ## Run 201 stress 0.3339536 
    ## Run 202 stress 0.3332009 
    ## Run 203 stress 0.3391468 
    ## Run 204 stress 0.3371942 
    ## Run 205 stress 0.3365802 
    ## Run 206 stress 0.3398538 
    ## Run 207 stress 0.3372515 
    ## Run 208 stress 0.3479114 
    ## Run 209 stress 0.339255 
    ## Run 210 stress 0.3367238 
    ## Run 211 stress 0.3423283 
    ## Run 212 stress 0.3458068 
    ## Run 213 stress 0.340338 
    ## Run 214 stress 0.3433331 
    ## Run 215 stress 0.3301571 
    ## Run 216 stress 0.335421 
    ## Run 217 stress 0.3326543 
    ## Run 218 stress 0.338818 
    ## Run 219 stress 0.331696 
    ## Run 220 stress 0.3409211 
    ## Run 221 stress 0.3276307 
    ## Run 222 stress 0.34217 
    ## Run 223 stress 0.3344409 
    ## Run 224 stress 0.3441779 
    ## Run 225 stress 0.3497043 
    ## Run 226 stress 0.3452665 
    ## Run 227 stress 0.338802 
    ## Run 228 stress 0.3335768 
    ## Run 229 stress 0.3284577 
    ## Run 230 stress 0.3317051 
    ## Run 231 stress 0.3301637 
    ## Run 232 stress 0.3351514 
    ## Run 233 stress 0.3372208 
    ## Run 234 stress 0.3464213 
    ## Run 235 stress 0.3415771 
    ## Run 236 stress 0.3376187 
    ## Run 237 stress 0.3461666 
    ## Run 238 stress 0.3327579 
    ## Run 239 stress 0.3353412 
    ## Run 240 stress 0.337401 
    ## Run 241 stress 0.3487918 
    ## Run 242 stress 0.3341217 
    ## Run 243 stress 0.3459297 
    ## Run 244 stress 0.3388976 
    ## Run 245 stress 0.3324412 
    ## Run 246 stress 0.3391332 
    ## Run 247 stress 0.336486 
    ## Run 248 stress 0.3336161 
    ## Run 249 stress 0.330099 
    ## Run 250 stress 0.3387445 
    ## Run 251 stress 0.3333039 
    ## Run 252 stress 0.337562 
    ## Run 253 stress 0.3276275 
    ## Run 254 stress 0.3426449 
    ## Run 255 stress 0.3300677 
    ## Run 256 stress 0.3390506 
    ## Run 257 stress 0.3347102 
    ## Run 258 stress 0.3360186 
    ## Run 259 stress 0.3402413 
    ## Run 260 stress 0.3315791 
    ## Run 261 stress 0.3313651 
    ## Run 262 stress 0.3282535 
    ## Run 263 stress 0.3349381 
    ## Run 264 stress 0.3429691 
    ## Run 265 stress 0.3319365 
    ## Run 266 stress 0.3329796 
    ## Run 267 stress 0.3287811 
    ## Run 268 stress 0.3312572 
    ## Run 269 stress 0.3477099 
    ## Run 270 stress 0.3415574 
    ## Run 271 stress 0.3367802 
    ## Run 272 stress 0.3463022 
    ## Run 273 stress 0.3378215 
    ## Run 274 stress 0.3450312 
    ## Run 275 stress 0.3274177 
    ## Run 276 stress 0.3450135 
    ## Run 277 stress 0.3384852 
    ## Run 278 stress 0.3310594 
    ## Run 279 stress 0.326647 
    ## Run 280 stress 0.3423183 
    ## Run 281 stress 0.342361 
    ## Run 282 stress 0.328326 
    ## Run 283 stress 0.3325008 
    ## Run 284 stress 0.3270742 
    ## Run 285 stress 0.3453253 
    ## Run 286 stress 0.3424762 
    ## Run 287 stress 0.3416774 
    ## Run 288 stress 0.3271887 
    ## Run 289 stress 0.3312615 
    ## Run 290 stress 0.3306418 
    ## Run 291 stress 0.3450883 
    ## Run 292 stress 0.3328289 
    ## Run 293 stress 0.3363632 
    ## Run 294 stress 0.3360192 
    ## Run 295 stress 0.3307591 
    ## Run 296 stress 0.3420011 
    ## Run 297 stress 0.3383462 
    ## Run 298 stress 0.3423872 
    ## Run 299 stress 0.3290949 
    ## Run 300 stress 0.3288042 
    ## Run 301 stress 0.3298677 
    ## Run 302 stress 0.3325427 
    ## Run 303 stress 0.3313848 
    ## Run 304 stress 0.3339915 
    ## Run 305 stress 0.3408928 
    ## Run 306 stress 0.3308685 
    ## Run 307 stress 0.3483708 
    ## Run 308 stress 0.3360886 
    ## Run 309 stress 0.3451515 
    ## Run 310 stress 0.3406811 
    ## Run 311 stress 0.3382779 
    ## Run 312 stress 0.3341657 
    ## Run 313 stress 0.3325174 
    ## Run 314 stress 0.3353834 
    ## Run 315 stress 0.3323063 
    ## Run 316 stress 0.3305454 
    ## Run 317 stress 0.3384583 
    ## Run 318 stress 0.3346242 
    ## Run 319 stress 0.3462414 
    ## Run 320 stress 0.3283292 
    ## Run 321 stress 0.3435624 
    ## Run 322 stress 0.3332996 
    ## Run 323 stress 0.3393289 
    ## Run 324 stress 0.3312976 
    ## Run 325 stress 0.3316343 
    ## Run 326 stress 0.3272862 
    ## Run 327 stress 0.3391592 
    ## Run 328 stress 0.3386346 
    ## Run 329 stress 0.3408712 
    ## Run 330 stress 0.3338393 
    ## Run 331 stress 0.3516617 
    ## Run 332 stress 0.3487116 
    ## Run 333 stress 0.3364274 
    ## Run 334 stress 0.3381845 
    ## Run 335 stress 0.3362271 
    ## Run 336 stress 0.3436912 
    ## Run 337 stress 0.3463537 
    ## Run 338 stress 0.3353951 
    ## Run 339 stress 0.335266 
    ## Run 340 stress 0.344893 
    ## Run 341 stress 0.329421 
    ## Run 342 stress 0.3510964 
    ## Run 343 stress 0.3371341 
    ## Run 344 stress 0.3322303 
    ## Run 345 stress 0.3368047 
    ## Run 346 stress 0.3320175 
    ## Run 347 stress 0.3351794 
    ## Run 348 stress 0.3469221 
    ## Run 349 stress 0.3300805 
    ## Run 350 stress 0.3495855 
    ## Run 351 stress 0.3586376 
    ## Run 352 stress 0.3321288 
    ## Run 353 stress 0.3298286 
    ## Run 354 stress 0.3387824 
    ## Run 355 stress 0.3306316 
    ## Run 356 stress 0.3464725 
    ## Run 357 stress 0.3438852 
    ## Run 358 stress 0.3434475 
    ## Run 359 stress 0.3350762 
    ## Run 360 stress 0.3347775 
    ## Run 361 stress 0.3383303 
    ## Run 362 stress 0.3436871 
    ## Run 363 stress 0.3322377 
    ## Run 364 stress 0.3322084 
    ## Run 365 stress 0.3396524 
    ## Run 366 stress 0.3343758 
    ## Run 367 stress 0.3285498 
    ## Run 368 stress 0.3494227 
    ## Run 369 stress 0.3354337 
    ## Run 370 stress 0.3296573 
    ## Run 371 stress 0.3367396 
    ## Run 372 stress 0.3451383 
    ## Run 373 stress 0.3302975 
    ## Run 374 stress 0.3325745 
    ## Run 375 stress 0.335965 
    ## Run 376 stress 0.3348493 
    ## Run 377 stress 0.3406791 
    ## Run 378 stress 0.3340341 
    ## Run 379 stress 0.3369316 
    ## Run 380 stress 0.3311781 
    ## Run 381 stress 0.3324124 
    ## Run 382 stress 0.3332508 
    ## Run 383 stress 0.3468928 
    ## Run 384 stress 0.3372599 
    ## Run 385 stress 0.3329757 
    ## Run 386 stress 0.3389478 
    ## Run 387 stress 0.3302331 
    ## Run 388 stress 0.3339792 
    ## Run 389 stress 0.3469041 
    ## Run 390 stress 0.3288077 
    ## Run 391 stress 0.3357564 
    ## Run 392 stress 0.3392322 
    ## Run 393 stress 0.3471547 
    ## Run 394 stress 0.3310347 
    ## Run 395 stress 0.3387824 
    ## Run 396 stress 0.3332794 
    ## Run 397 stress 0.3392693 
    ## Run 398 stress 0.3415964 
    ## Run 399 stress 0.3354021 
    ## Run 400 stress 0.3404188 
    ## Run 401 stress 0.3434067 
    ## Run 402 stress 0.3480428 
    ## Run 403 stress 0.3371717 
    ## Run 404 stress 0.3433285 
    ## Run 405 stress 0.3376022 
    ## Run 406 stress 0.3332544 
    ## Run 407 stress 0.341873 
    ## Run 408 stress 0.3345908 
    ## Run 409 stress 0.3351288 
    ## Run 410 stress 0.3323179 
    ## Run 411 stress 0.3327278 
    ## Run 412 stress 0.3337442 
    ## Run 413 stress 0.3389151 
    ## Run 414 stress 0.3363576 
    ## Run 415 stress 0.3315407 
    ## Run 416 stress 0.3451787 
    ## Run 417 stress 0.333958 
    ## Run 418 stress 0.3387611 
    ## Run 419 stress 0.3421256 
    ## Run 420 stress 0.336987 
    ## Run 421 stress 0.3397518 
    ## Run 422 stress 0.343032 
    ## Run 423 stress 0.3412077 
    ## Run 424 stress 0.3397638 
    ## Run 425 stress 0.3347666 
    ## Run 426 stress 0.34301 
    ## Run 427 stress 0.3366934 
    ## Run 428 stress 0.3430105 
    ## Run 429 stress 0.345173 
    ## Run 430 stress 0.3312092 
    ## Run 431 stress 0.3332438 
    ## Run 432 stress 0.3396772 
    ## Run 433 stress 0.3411422 
    ## Run 434 stress 0.3396067 
    ## Run 435 stress 0.3326692 
    ## Run 436 stress 0.3417766 
    ## Run 437 stress 0.3313559 
    ## Run 438 stress 0.3312047 
    ## Run 439 stress 0.3318501 
    ## Run 440 stress 0.3301771 
    ## Run 441 stress 0.3272557 
    ## Run 442 stress 0.3306222 
    ## Run 443 stress 0.3418833 
    ## Run 444 stress 0.338653 
    ## Run 445 stress 0.3330074 
    ## Run 446 stress 0.3323304 
    ## Run 447 stress 0.3315667 
    ## Run 448 stress 0.3425225 
    ## Run 449 stress 0.3403069 
    ## Run 450 stress 0.3263427 
    ## Run 451 stress 0.3475483 
    ## Run 452 stress 0.335549 
    ## Run 453 stress 0.3265491 
    ## Run 454 stress 0.3328157 
    ## Run 455 stress 0.3269641 
    ## Run 456 stress 0.3367019 
    ## Run 457 stress 0.3385347 
    ## Run 458 stress 0.3308772 
    ## Run 459 stress 0.3344515 
    ## Run 460 stress 0.3328896 
    ## Run 461 stress 0.3472763 
    ## Run 462 stress 0.3310567 
    ## Run 463 stress 0.3345094 
    ## Run 464 stress 0.3461267 
    ## Run 465 stress 0.3360306 
    ## Run 466 stress 0.3491614 
    ## Run 467 stress 0.3277834 
    ## Run 468 stress 0.3439696 
    ## Run 469 stress 0.3314294 
    ## Run 470 stress 0.3356512 
    ## Run 471 stress 0.3340342 
    ## Run 472 stress 0.3333477 
    ## Run 473 stress 0.3374337 
    ## Run 474 stress 0.3484569 
    ## Run 475 stress 0.3426884 
    ## Run 476 stress 0.3332416 
    ## Run 477 stress 0.3339756 
    ## Run 478 stress 0.3323285 
    ## Run 479 stress 0.3401479 
    ## Run 480 stress 0.3334673 
    ## Run 481 stress 0.3432618 
    ## Run 482 stress 0.333268 
    ## Run 483 stress 0.3330933 
    ## Run 484 stress 0.3436084 
    ## Run 485 stress 0.342681 
    ## Run 486 stress 0.3337038 
    ## Run 487 stress 0.3341259 
    ## Run 488 stress 0.3379648 
    ## Run 489 stress 0.3327116 
    ## Run 490 stress 0.3346649 
    ## Run 491 stress 0.3450202 
    ## Run 492 stress 0.3261307 
    ## Run 493 stress 0.3344665 
    ## Run 494 stress 0.3336232 
    ## Run 495 stress 0.3346445 
    ## Run 496 stress 0.3456115 
    ## Run 497 stress 0.3370892 
    ## Run 498 stress 0.3284497 
    ## Run 499 stress 0.3311372 
    ## Run 500 stress 0.3253116 
    ## ... Procrustes: rmse 0.1435542  max resid 0.30756 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##    500: stress ratio > sratmax

``` r
SD_beta_ric_NMDS <- metaMDS(SD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.01879552 
    ## Run 1 stress 0.02492266 
    ## Run 2 stress 0.02486147 
    ## Run 3 stress 0.01879548 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001108455  max resid 0.002283015 
    ## ... Similar to previous best
    ## Run 4 stress 0.01879573 
    ## ... Procrustes: rmse 0.0001106299  max resid 0.0002274745 
    ## ... Similar to previous best
    ## Run 5 stress 0.02486138 
    ## Run 6 stress 0.02509463 
    ## Run 7 stress 0.01898896 
    ## ... Procrustes: rmse 0.01002382  max resid 0.02040751 
    ## Run 8 stress 0.01879601 
    ## ... Procrustes: rmse 0.0002311334  max resid 0.0004759161 
    ## ... Similar to previous best
    ## Run 9 stress 0.0189045 
    ## ... Procrustes: rmse 0.007659828  max resid 0.0157347 
    ## Run 10 stress 0.02486136 
    ## Run 11 stress 0.01879585 
    ## ... Procrustes: rmse 0.0001787881  max resid 0.0003683367 
    ## ... Similar to previous best
    ## Run 12 stress 0.01883273 
    ## ... Procrustes: rmse 0.002340837  max resid 0.005131153 
    ## ... Similar to previous best
    ## Run 13 stress 0.02509443 
    ## Run 14 stress 0.02486134 
    ## Run 15 stress 0.02509489 
    ## Run 16 stress 0.02492193 
    ## Run 17 stress 0.0189141 
    ## ... Procrustes: rmse 0.007311311  max resid 0.01479974 
    ## Run 18 stress 0.02509467 
    ## Run 19 stress 0.01879783 
    ## ... Procrustes: rmse 0.000765146  max resid 0.001576405 
    ## ... Similar to previous best
    ## Run 20 stress 0.02509445 
    ## Run 21 stress 0.01879594 
    ## ... Procrustes: rmse 0.001288178  max resid 0.002654646 
    ## ... Similar to previous best
    ## Run 22 stress 0.02509469 
    ## Run 23 stress 0.02520888 
    ## Run 24 stress 0.01879586 
    ## ... Procrustes: rmse 0.0001798143  max resid 0.0003704188 
    ## ... Similar to previous best
    ## Run 25 stress 0.0250945 
    ## Run 26 stress 0.02509447 
    ## Run 27 stress 0.02520909 
    ## Run 28 stress 0.02509454 
    ## Run 29 stress 0.01879565 
    ## ... Procrustes: rmse 8.943821e-05  max resid 0.0001842607 
    ## ... Similar to previous best
    ## Run 30 stress 0.02509477 
    ## Run 31 stress 0.01879888 
    ## ... Procrustes: rmse 0.0009846865  max resid 0.002032522 
    ## ... Similar to previous best
    ## Run 32 stress 0.01879537 
    ## ... New best solution
    ## ... Procrustes: rmse 5.911404e-05  max resid 0.0001217605 
    ## ... Similar to previous best
    ## Run 33 stress 0.01879559 
    ## ... Procrustes: rmse 0.0001102848  max resid 0.0002271453 
    ## ... Similar to previous best
    ## Run 34 stress 0.02509465 
    ## Run 35 stress 0.01893221 
    ## ... Procrustes: rmse 0.008707751  max resid 0.01787846 
    ## Run 36 stress 0.02509451 
    ## Run 37 stress 0.01884405 
    ## ... Procrustes: rmse 0.003580162  max resid 0.00695386 
    ## ... Similar to previous best
    ## Run 38 stress 0.02486136 
    ## Run 39 stress 0.02492191 
    ## Run 40 stress 0.02520878 
    ## Run 41 stress 0.02486098 
    ## Run 42 stress 0.01883191 
    ## ... Procrustes: rmse 0.002277922  max resid 0.00505952 
    ## ... Similar to previous best
    ## Run 43 stress 0.01883098 
    ## ... Procrustes: rmse 0.002152009  max resid 0.004939366 
    ## ... Similar to previous best
    ## Run 44 stress 0.0187956 
    ## ... Procrustes: rmse 0.0001181516  max resid 0.0002432966 
    ## ... Similar to previous best
    ## Run 45 stress 0.02486131 
    ## Run 46 stress 0.02486141 
    ## Run 47 stress 0.01889558 
    ## ... Procrustes: rmse 0.006706848  max resid 0.01353925 
    ## Run 48 stress 0.02509462 
    ## Run 49 stress 0.01879582 
    ## ... Procrustes: rmse 0.00119052  max resid 0.002453591 
    ## ... Similar to previous best
    ## Run 50 stress 0.01879549 
    ## ... Procrustes: rmse 0.00103342  max resid 0.002130806 
    ## ... Similar to previous best
    ## Run 51 stress 0.01879577 
    ## ... Procrustes: rmse 0.0001844639  max resid 0.0003795179 
    ## ... Similar to previous best
    ## Run 52 stress 0.2716757 
    ## Run 53 stress 0.02492199 
    ## Run 54 stress 0.01879585 
    ## ... Procrustes: rmse 0.0002338353  max resid 0.0004816405 
    ## ... Similar to previous best
    ## Run 55 stress 0.0187956 
    ## ... Procrustes: rmse 0.0001211146  max resid 0.0002494433 
    ## ... Similar to previous best
    ## Run 56 stress 0.01879573 
    ## ... Procrustes: rmse 0.0001833384  max resid 0.0003777252 
    ## ... Similar to previous best
    ## Run 57 stress 0.01979497 
    ## Run 58 stress 0.01887497 
    ## ... Procrustes: rmse 0.005666015  max resid 0.01136342 
    ## Run 59 stress 0.02520914 
    ## Run 60 stress 0.01879587 
    ## ... Procrustes: rmse 0.0002434569  max resid 0.0005015661 
    ## ... Similar to previous best
    ## Run 61 stress 0.02509455 
    ## Run 62 stress 0.01879928 
    ## ... Procrustes: rmse 0.001057112  max resid 0.002182721 
    ## ... Similar to previous best
    ## Run 63 stress 0.01879602 
    ## ... Procrustes: rmse 0.0002922925  max resid 0.0006018288 
    ## ... Similar to previous best
    ## Run 64 stress 0.01885189 
    ## ... Procrustes: rmse 0.00533887  max resid 0.01097665 
    ## Run 65 stress 0.02509455 
    ## Run 66 stress 0.02509467 
    ## Run 67 stress 0.01879584 
    ## ... Procrustes: rmse 0.000231245  max resid 0.0004764383 
    ## ... Similar to previous best
    ## Run 68 stress 0.01879596 
    ## ... Procrustes: rmse 0.001236388  max resid 0.002547829 
    ## ... Similar to previous best
    ## Run 69 stress 0.01882908 
    ## ... Procrustes: rmse 0.001866581  max resid 0.004459229 
    ## ... Similar to previous best
    ## Run 70 stress 0.02486143 
    ## Run 71 stress 0.02492183 
    ## Run 72 stress 0.01879882 
    ## ... Procrustes: rmse 0.001034342  max resid 0.002131459 
    ## ... Similar to previous best
    ## Run 73 stress 0.01879559 
    ## ... Procrustes: rmse 0.0001199436  max resid 0.0002469078 
    ## ... Similar to previous best
    ## Run 74 stress 0.02509444 
    ## Run 75 stress 0.02509446 
    ## Run 76 stress 0.02492332 
    ## Run 77 stress 0.01879568 
    ## ... Procrustes: rmse 0.001127624  max resid 0.002324354 
    ## ... Similar to previous best
    ## Run 78 stress 0.02486142 
    ## Run 79 stress 0.01885816 
    ## ... Procrustes: rmse 0.00458825  max resid 0.009100662 
    ## ... Similar to previous best
    ## Run 80 stress 0.02492193 
    ## Run 81 stress 0.01879563 
    ## ... Procrustes: rmse 0.0001376328  max resid 0.000283503 
    ## ... Similar to previous best
    ## Run 82 stress 0.01879571 
    ## ... Procrustes: rmse 0.0001767482  max resid 0.0003640522 
    ## ... Similar to previous best
    ## Run 83 stress 0.3097771 
    ## Run 84 stress 0.01879578 
    ## ... Procrustes: rmse 0.0002055151  max resid 0.0004234236 
    ## ... Similar to previous best
    ## Run 85 stress 0.01880138 
    ## ... Procrustes: rmse 0.001417569  max resid 0.002921948 
    ## ... Similar to previous best
    ## Run 86 stress 0.01879573 
    ## ... Procrustes: rmse 0.0001855121  max resid 0.0003821823 
    ## ... Similar to previous best
    ## Run 87 stress 0.01884245 
    ## ... Procrustes: rmse 0.003426926  max resid 0.006623223 
    ## ... Similar to previous best
    ## Run 88 stress 0.0189128 
    ## ... Procrustes: rmse 0.007477089  max resid 0.01514023 
    ## Run 89 stress 0.01879573 
    ## ... Procrustes: rmse 0.0001848346  max resid 0.0003808001 
    ## ... Similar to previous best
    ## Run 90 stress 0.01945072 
    ## Run 91 stress 0.01879554 
    ## ... Procrustes: rmse 9.424212e-05  max resid 0.0001941952 
    ## ... Similar to previous best
    ## Run 92 stress 0.02509449 
    ## Run 93 stress 0.02520918 
    ## Run 94 stress 0.01879887 
    ## ... Procrustes: rmse 0.0010426  max resid 0.002152327 
    ## ... Similar to previous best
    ## Run 95 stress 0.02509462 
    ## Run 96 stress 0.02509457 
    ## Run 97 stress 0.01939636 
    ## Run 98 stress 0.02486127 
    ## Run 99 stress 0.02509433 
    ## Run 100 stress 0.01879826 
    ## ... Procrustes: rmse 0.0009023071  max resid 0.001859252 
    ## ... Similar to previous best
    ## Run 101 stress 0.02492164 
    ## Run 102 stress 0.01879716 
    ## ... Procrustes: rmse 0.0006600814  max resid 0.001360015 
    ## ... Similar to previous best
    ## Run 103 stress 0.01879856 
    ## ... Procrustes: rmse 0.0009803095  max resid 0.002020059 
    ## ... Similar to previous best
    ## Run 104 stress 0.01883457 
    ## ... Procrustes: rmse 0.002590715  max resid 0.005304474 
    ## ... Similar to previous best
    ## Run 105 stress 0.02492194 
    ## Run 106 stress 0.02492174 
    ## Run 107 stress 0.02509458 
    ## Run 108 stress 0.01880518 
    ## ... Procrustes: rmse 0.001463087  max resid 0.002912287 
    ## ... Similar to previous best
    ## Run 109 stress 0.01879606 
    ## ... Procrustes: rmse 0.0003022753  max resid 0.000622691 
    ## ... Similar to previous best
    ## Run 110 stress 0.02509446 
    ## Run 111 stress 0.01886145 
    ## ... Procrustes: rmse 0.004858352  max resid 0.009667675 
    ## ... Similar to previous best
    ## Run 112 stress 0.01879576 
    ## ... Procrustes: rmse 0.0001962984  max resid 0.0004044568 
    ## ... Similar to previous best
    ## Run 113 stress 0.0187953 
    ## ... New best solution
    ## ... Procrustes: rmse 5.521606e-05  max resid 0.0001134708 
    ## ... Similar to previous best
    ## Run 114 stress 0.02520892 
    ## Run 115 stress 0.02509439 
    ## Run 116 stress 0.01926015 
    ## ... Procrustes: rmse 0.01659314  max resid 0.03410865 
    ## Run 117 stress 0.02517303 
    ## Run 118 stress 0.01883741 
    ## ... Procrustes: rmse 0.002933747  max resid 0.005556701 
    ## ... Similar to previous best
    ## Run 119 stress 0.02492181 
    ## Run 120 stress 0.01879584 
    ## ... Procrustes: rmse 0.0002869503  max resid 0.0005911953 
    ## ... Similar to previous best
    ## Run 121 stress 0.02492171 
    ## Run 122 stress 0.01879595 
    ## ... Procrustes: rmse 0.0003329415  max resid 0.0006858547 
    ## ... Similar to previous best
    ## Run 123 stress 0.01879585 
    ## ... Procrustes: rmse 0.0002913019  max resid 0.0006000335 
    ## ... Similar to previous best
    ## Run 124 stress 0.01879591 
    ## ... Procrustes: rmse 0.00031612  max resid 0.0006510752 
    ## ... Similar to previous best
    ## Run 125 stress 0.02509464 
    ## Run 126 stress 0.01879607 
    ## ... Procrustes: rmse 0.0003641592  max resid 0.0007502893 
    ## ... Similar to previous best
    ## Run 127 stress 0.02492185 
    ## Run 128 stress 0.02520865 
    ## Run 129 stress 0.02492188 
    ## Run 130 stress 0.02520909 
    ## Run 131 stress 0.01879544 
    ## ... Procrustes: rmse 8.67046e-05  max resid 0.0001786916 
    ## ... Similar to previous best
    ## Run 132 stress 0.01879593 
    ## ... Procrustes: rmse 0.0003247087  max resid 0.0006688409 
    ## ... Similar to previous best
    ## Run 133 stress 0.01879587 
    ## ... Procrustes: rmse 0.0003021016  max resid 0.0006222308 
    ## ... Similar to previous best
    ## Run 134 stress 0.01885003 
    ## ... Procrustes: rmse 0.003976044  max resid 0.007799653 
    ## ... Similar to previous best
    ## Run 135 stress 0.01879584 
    ## ... Procrustes: rmse 0.0002755463  max resid 0.0005670342 
    ## ... Similar to previous best
    ## Run 136 stress 0.02486134 
    ## Run 137 stress 0.01883539 
    ## ... Procrustes: rmse 0.002768525  max resid 0.005418476 
    ## ... Similar to previous best
    ## Run 138 stress 0.01887336 
    ## ... Procrustes: rmse 0.006500782  max resid 0.01336346 
    ## Run 139 stress 0.01879602 
    ## ... Procrustes: rmse 0.001200012  max resid 0.002472654 
    ## ... Similar to previous best
    ## Run 140 stress 0.01953257 
    ## Run 141 stress 0.02509448 
    ## Run 142 stress 0.01879575 
    ## ... Procrustes: rmse 0.001094529  max resid 0.002255882 
    ## ... Similar to previous best
    ## Run 143 stress 0.0187958 
    ## ... Procrustes: rmse 0.0002519075  max resid 0.000518042 
    ## ... Similar to previous best
    ## Run 144 stress 0.02492157 
    ## Run 145 stress 0.01880235 
    ## ... Procrustes: rmse 0.001629841  max resid 0.003357627 
    ## ... Similar to previous best
    ## Run 146 stress 0.02486125 
    ## Run 147 stress 0.01879555 
    ## ... Procrustes: rmse 0.0001537117  max resid 0.0003166293 
    ## ... Similar to previous best
    ## Run 148 stress 0.01898715 
    ## ... Procrustes: rmse 0.01041992  max resid 0.02137854 
    ## Run 149 stress 0.01882481 
    ## ... Procrustes: rmse 0.003811006  max resid 0.007844836 
    ## ... Similar to previous best
    ## Run 150 stress 0.01883207 
    ## ... Procrustes: rmse 0.002354372  max resid 0.005103906 
    ## ... Similar to previous best
    ## Run 151 stress 0.01883358 
    ## ... Procrustes: rmse 0.002320901  max resid 0.00509357 
    ## ... Similar to previous best
    ## Run 152 stress 0.01901717 
    ## ... Procrustes: rmse 0.01110639  max resid 0.02264243 
    ## Run 153 stress 0.01884084 
    ## ... Procrustes: rmse 0.003339475  max resid 0.006436009 
    ## ... Similar to previous best
    ## Run 154 stress 0.02509466 
    ## Run 155 stress 0.01903752 
    ## ... Procrustes: rmse 0.01184776  max resid 0.0242836 
    ## Run 156 stress 0.01883271 
    ## ... Procrustes: rmse 0.004357219  max resid 0.008967471 
    ## ... Similar to previous best
    ## Run 157 stress 0.02492266 
    ## Run 158 stress 0.01883219 
    ## ... Procrustes: rmse 0.002333193  max resid 0.005087492 
    ## ... Similar to previous best
    ## Run 159 stress 0.01879568 
    ## ... Procrustes: rmse 0.000217033  max resid 0.0004470911 
    ## ... Similar to previous best
    ## Run 160 stress 0.01879577 
    ## ... Procrustes: rmse 0.0002507565  max resid 0.0005160598 
    ## ... Similar to previous best
    ## Run 161 stress 0.01885144 
    ## ... Procrustes: rmse 0.00415523  max resid 0.008182557 
    ## ... Similar to previous best
    ## Run 162 stress 0.02486136 
    ## Run 163 stress 0.01887947 
    ## ... Procrustes: rmse 0.005715317  max resid 0.01146245 
    ## Run 164 stress 0.02509471 
    ## Run 165 stress 0.01881827 
    ## ... Procrustes: rmse 0.001595883  max resid 0.003039196 
    ## ... Similar to previous best
    ## Run 166 stress 0.2038604 
    ## Run 167 stress 0.01879573 
    ## ... Procrustes: rmse 0.0002414669  max resid 0.0004973903 
    ## ... Similar to previous best
    ## Run 168 stress 0.01882983 
    ## ... Procrustes: rmse 0.004173901  max resid 0.008589823 
    ## ... Similar to previous best
    ## Run 169 stress 0.02509478 
    ## Run 170 stress 0.02509438 
    ## Run 171 stress 0.02492171 
    ## Run 172 stress 0.0187995 
    ## ... Procrustes: rmse 0.001209216  max resid 0.00249138 
    ## ... Similar to previous best
    ## Run 173 stress 0.02509447 
    ## Run 174 stress 0.01882975 
    ## ... Procrustes: rmse 0.002014562  max resid 0.004699878 
    ## ... Similar to previous best
    ## Run 175 stress 0.02509488 
    ## Run 176 stress 0.01879584 
    ## ... Procrustes: rmse 0.0002860784  max resid 0.0005894022 
    ## ... Similar to previous best
    ## Run 177 stress 0.01879596 
    ## ... Procrustes: rmse 0.0003249763  max resid 0.0006695866 
    ## ... Similar to previous best
    ## Run 178 stress 0.02492204 
    ## Run 179 stress 0.01879568 
    ## ... Procrustes: rmse 0.0002104872  max resid 0.0004337241 
    ## ... Similar to previous best
    ## Run 180 stress 0.02509441 
    ## Run 181 stress 0.02509456 
    ## Run 182 stress 0.01879585 
    ## ... Procrustes: rmse 0.000290261  max resid 0.0005979517 
    ## ... Similar to previous best
    ## Run 183 stress 0.01879597 
    ## ... Procrustes: rmse 0.0003321864  max resid 0.0006838472 
    ## ... Similar to previous best
    ## Run 184 stress 0.01883885 
    ## ... Procrustes: rmse 0.003143775  max resid 0.006012734 
    ## ... Similar to previous best
    ## Run 185 stress 0.02486127 
    ## Run 186 stress 0.01881195 
    ## ... Procrustes: rmse 0.002772888  max resid 0.005711888 
    ## ... Similar to previous best
    ## Run 187 stress 0.01881216 
    ## ... Procrustes: rmse 0.002791652  max resid 0.005749873 
    ## ... Similar to previous best
    ## Run 188 stress 0.01895644 
    ## ... Procrustes: rmse 0.009185194  max resid 0.01868246 
    ## Run 189 stress 0.02520912 
    ## Run 190 stress 0.01905668 
    ## ... Procrustes: rmse 0.01233962  max resid 0.02528876 
    ## Run 191 stress 0.02509451 
    ## Run 192 stress 0.01879519 
    ## ... New best solution
    ## ... Procrustes: rmse 8.244056e-05  max resid 0.0001695376 
    ## ... Similar to previous best
    ## Run 193 stress 0.01879968 
    ## ... Procrustes: rmse 0.001331481  max resid 0.002743746 
    ## ... Similar to previous best
    ## Run 194 stress 0.02492198 
    ## Run 195 stress 0.02492188 
    ## Run 196 stress 0.02104334 
    ## Run 197 stress 0.0187956 
    ## ... Procrustes: rmse 0.0002623743  max resid 0.0005406484 
    ## ... Similar to previous best
    ## Run 198 stress 0.0188103 
    ## ... Procrustes: rmse 0.002684916  max resid 0.005530741 
    ## ... Similar to previous best
    ## Run 199 stress 0.02486143 
    ## Run 200 stress 0.01879557 
    ## ... Procrustes: rmse 0.0002424366  max resid 0.0004995966 
    ## ... Similar to previous best
    ## Run 201 stress 0.02509422 
    ## Run 202 stress 0.01882988 
    ## ... Procrustes: rmse 0.002034518  max resid 0.004675105 
    ## ... Similar to previous best
    ## Run 203 stress 0.02509442 
    ## Run 204 stress 0.02492183 
    ## Run 205 stress 0.01885884 
    ## ... Procrustes: rmse 0.004829835  max resid 0.009610243 
    ## ... Similar to previous best
    ## Run 206 stress 0.01883776 
    ## ... Procrustes: rmse 0.003107111  max resid 0.00593615 
    ## ... Similar to previous best
    ## Run 207 stress 0.02509465 
    ## Run 208 stress 0.02492459 
    ## Run 209 stress 0.01886656 
    ## ... Procrustes: rmse 0.00530736  max resid 0.01061156 
    ## Run 210 stress 0.01879773 
    ## ... Procrustes: rmse 0.0006979813  max resid 0.00143549 
    ## ... Similar to previous best
    ## Run 211 stress 0.02509478 
    ## Run 212 stress 0.01879543 
    ## ... Procrustes: rmse 0.0008640268  max resid 0.001781576 
    ## ... Similar to previous best
    ## Run 213 stress 0.02486105 
    ## Run 214 stress 0.01879567 
    ## ... Procrustes: rmse 0.000985378  max resid 0.002030881 
    ## ... Similar to previous best
    ## Run 215 stress 0.01879526 
    ## ... Procrustes: rmse 6.259176e-05  max resid 0.0001289383 
    ## ... Similar to previous best
    ## Run 216 stress 0.0188128 
    ## ... Procrustes: rmse 0.002930814  max resid 0.006035933 
    ## ... Similar to previous best
    ## Run 217 stress 0.02486118 
    ## Run 218 stress 0.01882894 
    ## ... Procrustes: rmse 0.001971513  max resid 0.004464975 
    ## ... Similar to previous best
    ## Run 219 stress 0.0187956 
    ## ... Procrustes: rmse 0.0002601871  max resid 0.0005361203 
    ## ... Similar to previous best
    ## Run 220 stress 0.0248613 
    ## Run 221 stress 0.01880373 
    ## ... Procrustes: rmse 0.001948057  max resid 0.004014068 
    ## ... Similar to previous best
    ## Run 222 stress 0.01879621 
    ## ... Procrustes: rmse 0.0004671903  max resid 0.0009615939 
    ## ... Similar to previous best
    ## Run 223 stress 0.02492181 
    ## Run 224 stress 0.01882628 
    ## ... Procrustes: rmse 0.001779547  max resid 0.003623445 
    ## ... Similar to previous best
    ## Run 225 stress 0.02509463 
    ## Run 226 stress 0.01884601 
    ## ... Procrustes: rmse 0.00387448  max resid 0.007586156 
    ## ... Similar to previous best
    ## Run 227 stress 0.01885173 
    ## ... Procrustes: rmse 0.004328576  max resid 0.008553272 
    ## ... Similar to previous best
    ## Run 228 stress 0.02493016 
    ## Run 229 stress 0.01883225 
    ## ... Procrustes: rmse 0.001693353  max resid 0.003239125 
    ## ... Similar to previous best
    ## Run 230 stress 0.01879591 
    ## ... Procrustes: rmse 0.001072295  max resid 0.00220947 
    ## ... Similar to previous best
    ## Run 231 stress 0.02486135 
    ## Run 232 stress 0.01880249 
    ## ... Procrustes: rmse 0.001779688  max resid 0.003667274 
    ## ... Similar to previous best
    ## Run 233 stress 0.01882034 
    ## ... Procrustes: rmse 0.003481594  max resid 0.007163879 
    ## ... Similar to previous best
    ## Run 234 stress 0.01879567 
    ## ... Procrustes: rmse 0.0002936436  max resid 0.000605058 
    ## ... Similar to previous best
    ## Run 235 stress 0.01879547 
    ## ... Procrustes: rmse 0.0001927828  max resid 0.0003972783 
    ## ... Similar to previous best
    ## Run 236 stress 0.01879584 
    ## ... Procrustes: rmse 0.0003662194  max resid 0.0007545023 
    ## ... Similar to previous best
    ## Run 237 stress 0.01879553 
    ## ... Procrustes: rmse 0.0002269151  max resid 0.0004675991 
    ## ... Similar to previous best
    ## Run 238 stress 0.02492192 
    ## Run 239 stress 0.02492196 
    ## Run 240 stress 0.02492176 
    ## Run 241 stress 0.02486122 
    ## Run 242 stress 0.3155852 
    ## Run 243 stress 0.0250945 
    ## Run 244 stress 0.01883548 
    ## ... Procrustes: rmse 0.002849477  max resid 0.005454209 
    ## ... Similar to previous best
    ## Run 245 stress 0.02509471 
    ## Run 246 stress 0.0187957 
    ## ... Procrustes: rmse 0.0003053624  max resid 0.000629178 
    ## ... Similar to previous best
    ## Run 247 stress 0.01879606 
    ## ... Procrustes: rmse 0.0004345401  max resid 0.0008949604 
    ## ... Similar to previous best
    ## Run 248 stress 0.02492191 
    ## Run 249 stress 0.02509463 
    ## Run 250 stress 0.01879571 
    ## ... Procrustes: rmse 0.0003125939  max resid 0.0006440852 
    ## ... Similar to previous best
    ## Run 251 stress 0.02509472 
    ## Run 252 stress 0.02509462 
    ## Run 253 stress 0.02492189 
    ## Run 254 stress 0.02492174 
    ## Run 255 stress 0.01885318 
    ## ... Procrustes: rmse 0.005615548  max resid 0.01154848 
    ## Run 256 stress 0.01879572 
    ## ... Procrustes: rmse 0.00100542  max resid 0.002072069 
    ## ... Similar to previous best
    ## Run 257 stress 0.02486125 
    ## Run 258 stress 0.01884714 
    ## ... Procrustes: rmse 0.003969301  max resid 0.00778847 
    ## ... Similar to previous best
    ## Run 259 stress 0.0188996 
    ## ... Procrustes: rmse 0.005387621  max resid 0.01077428 
    ## Run 260 stress 0.02509441 
    ## Run 261 stress 0.02509465 
    ## Run 262 stress 0.01879764 
    ## ... Procrustes: rmse 0.0009155659  max resid 0.001885971 
    ## ... Similar to previous best
    ## Run 263 stress 0.01880221 
    ## ... Procrustes: rmse 0.001604934  max resid 0.003304444 
    ## ... Similar to previous best
    ## Run 264 stress 0.01879792 
    ## ... Procrustes: rmse 0.0009727846  max resid 0.002004927 
    ## ... Similar to previous best
    ## Run 265 stress 0.02486147 
    ## Run 266 stress 0.0249217 
    ## Run 267 stress 0.02486138 
    ## Run 268 stress 0.02520905 
    ## Run 269 stress 0.02509439 
    ## Run 270 stress 0.02509464 
    ## Run 271 stress 0.01879565 
    ## ... Procrustes: rmse 0.0002842108  max resid 0.0005856315 
    ## ... Similar to previous best
    ## Run 272 stress 0.02509459 
    ## Run 273 stress 0.01879841 
    ## ... Procrustes: rmse 0.001087339  max resid 0.002240239 
    ## ... Similar to previous best
    ## Run 274 stress 0.01882916 
    ## ... Procrustes: rmse 0.001990084  max resid 0.00451738 
    ## ... Similar to previous best
    ## Run 275 stress 0.01885553 
    ## ... Procrustes: rmse 0.004567229  max resid 0.009056018 
    ## ... Similar to previous best
    ## Run 276 stress 0.02520897 
    ## Run 277 stress 0.0188772 
    ## ... Procrustes: rmse 0.005424167  max resid 0.01085375 
    ## Run 278 stress 0.02509442 
    ## Run 279 stress 0.01879573 
    ## ... Procrustes: rmse 0.0003241965  max resid 0.0006679948 
    ## ... Similar to previous best
    ## Run 280 stress 0.02509434 
    ## Run 281 stress 0.02520889 
    ## Run 282 stress 0.02520903 
    ## Run 283 stress 0.024922 
    ## Run 284 stress 0.02509464 
    ## Run 285 stress 0.01926894 
    ## ... Procrustes: rmse 0.01678488  max resid 0.03448838 
    ## Run 286 stress 0.02492196 
    ## Run 287 stress 0.02509466 
    ## Run 288 stress 0.01879575 
    ## ... Procrustes: rmse 0.0003323291  max resid 0.0006847223 
    ## ... Similar to previous best
    ## Run 289 stress 0.02509454 
    ## Run 290 stress 0.02486115 
    ## Run 291 stress 0.01879568 
    ## ... Procrustes: rmse 0.0002984271  max resid 0.0006149243 
    ## ... Similar to previous best
    ## Run 292 stress 0.01879561 
    ## ... Procrustes: rmse 0.0002662607  max resid 0.000548727 
    ## ... Similar to previous best
    ## Run 293 stress 0.01879554 
    ## ... Procrustes: rmse 0.0002304031  max resid 0.0004747907 
    ## ... Similar to previous best
    ## Run 294 stress 0.01879584 
    ## ... Procrustes: rmse 0.001052267  max resid 0.002168322 
    ## ... Similar to previous best
    ## Run 295 stress 0.01880467 
    ## ... Procrustes: rmse 0.002026031  max resid 0.004172803 
    ## ... Similar to previous best
    ## Run 296 stress 0.01879613 
    ## ... Procrustes: rmse 0.000468797  max resid 0.0009656083 
    ## ... Similar to previous best
    ## Run 297 stress 0.01881164 
    ## ... Procrustes: rmse 0.002776582  max resid 0.005720382 
    ## ... Similar to previous best
    ## Run 298 stress 0.01879587 
    ## ... Procrustes: rmse 0.001064964  max resid 0.002194766 
    ## ... Similar to previous best
    ## Run 299 stress 0.01936436 
    ## Run 300 stress 0.01884473 
    ## ... Procrustes: rmse 0.00376977  max resid 0.007362596 
    ## ... Similar to previous best
    ## Run 301 stress 0.01879589 
    ## ... Procrustes: rmse 0.0003904656  max resid 0.0008044661 
    ## ... Similar to previous best
    ## Run 302 stress 0.01879576 
    ## ... Procrustes: rmse 0.0003317577  max resid 0.0006835407 
    ## ... Similar to previous best
    ## Run 303 stress 0.0249216 
    ## Run 304 stress 0.02509464 
    ## Run 305 stress 0.01879567 
    ## ... Procrustes: rmse 0.0002930729  max resid 0.0006038744 
    ## ... Similar to previous best
    ## Run 306 stress 0.01879609 
    ## ... Procrustes: rmse 0.0004521403  max resid 0.0009312774 
    ## ... Similar to previous best
    ## Run 307 stress 0.01879577 
    ## ... Procrustes: rmse 0.0003386195  max resid 0.0006976987 
    ## ... Similar to previous best
    ## Run 308 stress 0.02492272 
    ## Run 309 stress 0.02509455 
    ## Run 310 stress 0.02509472 
    ## Run 311 stress 0.02492193 
    ## Run 312 stress 0.02492463 
    ## Run 313 stress 0.01883456 
    ## ... Procrustes: rmse 0.002436589  max resid 0.00515412 
    ## ... Similar to previous best
    ## Run 314 stress 0.01879579 
    ## ... Procrustes: rmse 0.0003490831  max resid 0.0007192476 
    ## ... Similar to previous best
    ## Run 315 stress 0.02486107 
    ## Run 316 stress 0.02486127 
    ## Run 317 stress 0.01887021 
    ## ... Procrustes: rmse 0.005534827  max resid 0.01108924 
    ## Run 318 stress 0.01879588 
    ## ... Procrustes: rmse 0.0003873779  max resid 0.0007981352 
    ## ... Similar to previous best
    ## Run 319 stress 0.01899849 
    ## ... Procrustes: rmse 0.01061474  max resid 0.0216317 
    ## Run 320 stress 0.02509457 
    ## Run 321 stress 0.02509462 
    ## Run 322 stress 0.01879562 
    ## ... Procrustes: rmse 0.0009618392  max resid 0.001982538 
    ## ... Similar to previous best
    ## Run 323 stress 0.01879558 
    ## ... Procrustes: rmse 0.0002491516  max resid 0.000513485 
    ## ... Similar to previous best
    ## Run 324 stress 0.02492178 
    ## Run 325 stress 0.02016619 
    ## Run 326 stress 0.0187959 
    ## ... Procrustes: rmse 0.0003628993  max resid 0.0007474742 
    ## ... Similar to previous best
    ## Run 327 stress 0.01879586 
    ## ... Procrustes: rmse 0.000367101  max resid 0.0007562945 
    ## ... Similar to previous best
    ## Run 328 stress 0.01879883 
    ## ... Procrustes: rmse 0.001173078  max resid 0.002417056 
    ## ... Similar to previous best
    ## Run 329 stress 0.01879578 
    ## ... Procrustes: rmse 0.0003390989  max resid 0.0006985928 
    ## ... Similar to previous best
    ## Run 330 stress 0.01954926 
    ## Run 331 stress 0.01879567 
    ## ... Procrustes: rmse 0.0002903037  max resid 0.0005980951 
    ## ... Similar to previous best
    ## Run 332 stress 0.02492168 
    ## Run 333 stress 0.01879555 
    ## ... Procrustes: rmse 0.000235831  max resid 0.0004859922 
    ## ... Similar to previous best
    ## Run 334 stress 0.02492207 
    ## Run 335 stress 0.01892174 
    ## ... Procrustes: rmse 0.007979302  max resid 0.01618339 
    ## Run 336 stress 0.01879576 
    ## ... Procrustes: rmse 0.0003156579  max resid 0.0006501461 
    ## ... Similar to previous best
    ## Run 337 stress 0.02509434 
    ## Run 338 stress 0.02492189 
    ## Run 339 stress 0.02509465 
    ## Run 340 stress 0.01894225 
    ## ... Procrustes: rmse 0.008766941  max resid 0.01781584 
    ## Run 341 stress 0.01879587 
    ## ... Procrustes: rmse 0.0003674313  max resid 0.0007569669 
    ## ... Similar to previous best
    ## Run 342 stress 0.2038604 
    ## Run 343 stress 0.01884145 
    ## ... Procrustes: rmse 0.003471991  max resid 0.006724156 
    ## ... Similar to previous best
    ## Run 344 stress 0.02492183 
    ## Run 345 stress 0.02509473 
    ## Run 346 stress 0.01879592 
    ## ... Procrustes: rmse 0.0003886143  max resid 0.0008004831 
    ## ... Similar to previous best
    ## Run 347 stress 0.01885268 
    ## ... Procrustes: rmse 0.005577697  max resid 0.01147009 
    ## Run 348 stress 0.01880338 
    ## ... Procrustes: rmse 0.001870496  max resid 0.003853519 
    ## ... Similar to previous best
    ## Run 349 stress 0.02492176 
    ## Run 350 stress 0.01903977 
    ## ... Procrustes: rmse 0.01180997  max resid 0.02409352 
    ## Run 351 stress 0.0187958 
    ## ... Procrustes: rmse 0.0003550264  max resid 0.000731492 
    ## ... Similar to previous best
    ## Run 352 stress 0.01879579 
    ## ... Procrustes: rmse 0.001040504  max resid 0.002144125 
    ## ... Similar to previous best
    ## Run 353 stress 0.01879589 
    ## ... Procrustes: rmse 0.0003790574  max resid 0.0007810563 
    ## ... Similar to previous best
    ## Run 354 stress 0.02492195 
    ## Run 355 stress 0.01885011 
    ## ... Procrustes: rmse 0.005454221  max resid 0.01121743 
    ## Run 356 stress 0.01879564 
    ## ... Procrustes: rmse 0.0002785159  max resid 0.000573897 
    ## ... Similar to previous best
    ## Run 357 stress 0.02492197 
    ## Run 358 stress 0.01879597 
    ## ... Procrustes: rmse 0.000411178  max resid 0.0008469876 
    ## ... Similar to previous best
    ## Run 359 stress 0.0249244 
    ## Run 360 stress 0.01885751 
    ## ... Procrustes: rmse 0.004690059  max resid 0.009317855 
    ## ... Similar to previous best
    ## Run 361 stress 0.02520914 
    ## Run 362 stress 0.01887092 
    ## ... Procrustes: rmse 0.006444933  max resid 0.01325067 
    ## Run 363 stress 0.02509443 
    ## Run 364 stress 0.0250946 
    ## Run 365 stress 0.2038604 
    ## Run 366 stress 0.01879554 
    ## ... Procrustes: rmse 0.000924634  max resid 0.001906076 
    ## ... Similar to previous best
    ## Run 367 stress 0.0187956 
    ## ... Procrustes: rmse 0.000262091  max resid 0.0005400662 
    ## ... Similar to previous best
    ## Run 368 stress 0.01879545 
    ## ... Procrustes: rmse 0.0001816721  max resid 0.0003743771 
    ## ... Similar to previous best
    ## Run 369 stress 0.02509467 
    ## Run 370 stress 0.02509451 
    ## Run 371 stress 0.01879552 
    ## ... Procrustes: rmse 0.0002203034  max resid 0.0004539551 
    ## ... Similar to previous best
    ## Run 372 stress 0.0250943 
    ## Run 373 stress 0.02492194 
    ## Run 374 stress 0.01879595 
    ## ... Procrustes: rmse 0.0004140178  max resid 0.0008529953 
    ## ... Similar to previous best
    ## Run 375 stress 0.02509454 
    ## Run 376 stress 0.01879747 
    ## ... Procrustes: rmse 0.0008755543  max resid 0.001804096 
    ## ... Similar to previous best
    ## Run 377 stress 0.01879588 
    ## ... Procrustes: rmse 0.0003732776  max resid 0.0007689975 
    ## ... Similar to previous best
    ## Run 378 stress 0.02520878 
    ## Run 379 stress 0.02509466 
    ## Run 380 stress 0.02509445 
    ## Run 381 stress 0.0250948 
    ## Run 382 stress 0.01879739 
    ## ... Procrustes: rmse 0.0008489377  max resid 0.001748931 
    ## ... Similar to previous best
    ## Run 383 stress 0.01879594 
    ## ... Procrustes: rmse 0.0004101049  max resid 0.0008448847 
    ## ... Similar to previous best
    ## Run 384 stress 0.02509479 
    ## Run 385 stress 0.01879901 
    ## ... Procrustes: rmse 0.001196184  max resid 0.002464347 
    ## ... Similar to previous best
    ## Run 386 stress 0.02520889 
    ## Run 387 stress 0.02520887 
    ## Run 388 stress 0.01879577 
    ## ... Procrustes: rmse 0.0003396673  max resid 0.0006998312 
    ## ... Similar to previous best
    ## Run 389 stress 0.0248611 
    ## Run 390 stress 0.02509462 
    ## Run 391 stress 0.02486151 
    ## Run 392 stress 0.02047069 
    ## Run 393 stress 0.01879561 
    ## ... Procrustes: rmse 0.0002674077  max resid 0.0005510203 
    ## ... Similar to previous best
    ## Run 394 stress 0.02492184 
    ## Run 395 stress 0.01879585 
    ## ... Procrustes: rmse 0.0003759082  max resid 0.0007745047 
    ## ... Similar to previous best
    ## Run 396 stress 0.01879566 
    ## ... Procrustes: rmse 0.0009605516  max resid 0.001982817 
    ## ... Similar to previous best
    ## Run 397 stress 0.02492183 
    ## Run 398 stress 0.01879517 
    ## ... New best solution
    ## ... Procrustes: rmse 1.233297e-05  max resid 2.519583e-05 
    ## ... Similar to previous best
    ## Run 399 stress 0.01879574 
    ## ... Procrustes: rmse 0.0003336868  max resid 0.0006874689 
    ## ... Similar to previous best
    ## Run 400 stress 0.02509449 
    ## Run 401 stress 0.0194364 
    ## Run 402 stress 0.01907533 
    ## ... Procrustes: rmse 0.01230364  max resid 0.02521408 
    ## Run 403 stress 0.01879579 
    ## ... Procrustes: rmse 0.0003553307  max resid 0.0007320546 
    ## ... Similar to previous best
    ## Run 404 stress 0.01906761 
    ## ... Procrustes: rmse 0.01159442  max resid 0.02363791 
    ## Run 405 stress 0.01882592 
    ## ... Procrustes: rmse 0.003992276  max resid 0.008214473 
    ## ... Similar to previous best
    ## Run 406 stress 0.01879573 
    ## ... Procrustes: rmse 0.0003337569  max resid 0.0006876161 
    ## ... Similar to previous best
    ## Run 407 stress 0.01892101 
    ## ... Procrustes: rmse 0.007963583  max resid 0.01615211 
    ## Run 408 stress 0.02492203 
    ## Run 409 stress 0.02486133 
    ## Run 410 stress 0.01879566 
    ## ... Procrustes: rmse 0.0003003476  max resid 0.0006187822 
    ## ... Similar to previous best
    ## Run 411 stress 0.01879577 
    ## ... Procrustes: rmse 0.0003506974  max resid 0.0007225079 
    ## ... Similar to previous best
    ## Run 412 stress 0.02509466 
    ## Run 413 stress 0.01879554 
    ## ... Procrustes: rmse 0.0009121179  max resid 0.001880334 
    ## ... Similar to previous best
    ## Run 414 stress 0.01879598 
    ## ... Procrustes: rmse 0.000427664  max resid 0.0008810344 
    ## ... Similar to previous best
    ## Run 415 stress 0.02463292 
    ## Run 416 stress 0.0249219 
    ## Run 417 stress 0.01888838 
    ## ... Procrustes: rmse 0.006512448  max resid 0.01313315 
    ## Run 418 stress 0.02486142 
    ## Run 419 stress 0.02492179 
    ## Run 420 stress 0.02492233 
    ## Run 421 stress 0.02492177 
    ## Run 422 stress 0.01880473 
    ## ... Procrustes: rmse 0.002078345  max resid 0.004282887 
    ## ... Similar to previous best
    ## Run 423 stress 0.02509472 
    ## Run 424 stress 0.02492188 
    ## Run 425 stress 0.02509462 
    ## Run 426 stress 0.02492852 
    ## Run 427 stress 0.01879584 
    ## ... Procrustes: rmse 0.001050836  max resid 0.002165327 
    ## ... Similar to previous best
    ## Run 428 stress 0.02492204 
    ## Run 429 stress 0.01887251 
    ## ... Procrustes: rmse 0.006548608  max resid 0.01346222 
    ## Run 430 stress 0.0187982 
    ## ... Procrustes: rmse 0.001052506  max resid 0.002168048 
    ## ... Similar to previous best
    ## Run 431 stress 0.01883065 
    ## ... Procrustes: rmse 0.002217732  max resid 0.004915072 
    ## ... Similar to previous best
    ## Run 432 stress 0.01884554 
    ## ... Procrustes: rmse 0.003841255  max resid 0.0075164 
    ## ... Similar to previous best
    ## Run 433 stress 0.02520913 
    ## Run 434 stress 0.02509473 
    ## Run 435 stress 0.02492194 
    ## Run 436 stress 0.01883439 
    ## ... Procrustes: rmse 0.002712033  max resid 0.005353776 
    ## ... Similar to previous best
    ## Run 437 stress 0.02509435 
    ## Run 438 stress 0.02509453 
    ## Run 439 stress 0.01879592 
    ## ... Procrustes: rmse 0.0004023598  max resid 0.0008287096 
    ## ... Similar to previous best
    ## Run 440 stress 0.02509453 
    ## Run 441 stress 0.01882961 
    ## ... Procrustes: rmse 0.00207449  max resid 0.004700971 
    ## ... Similar to previous best
    ## Run 442 stress 0.02492189 
    ## Run 443 stress 0.01879588 
    ## ... Procrustes: rmse 0.0003975567  max resid 0.0008190264 
    ## ... Similar to previous best
    ## Run 444 stress 0.02492192 
    ## Run 445 stress 0.01879569 
    ## ... Procrustes: rmse 0.0003173412  max resid 0.0006537364 
    ## ... Similar to previous best
    ## Run 446 stress 0.01883834 
    ## ... Procrustes: rmse 0.004802804  max resid 0.009880396 
    ## ... Similar to previous best
    ## Run 447 stress 0.02486108 
    ## Run 448 stress 0.01881494 
    ## ... Procrustes: rmse 0.003140953  max resid 0.006468604 
    ## ... Similar to previous best
    ## Run 449 stress 0.01883435 
    ## ... Procrustes: rmse 0.004561227  max resid 0.009384206 
    ## ... Similar to previous best
    ## Run 450 stress 0.0252087 
    ## Run 451 stress 0.01879599 
    ## ... Procrustes: rmse 0.001097047  max resid 0.00226018 
    ## ... Similar to previous best
    ## Run 452 stress 0.02509471 
    ## Run 453 stress 0.02492483 
    ## Run 454 stress 0.02509445 
    ## Run 455 stress 0.02509476 
    ## Run 456 stress 0.01886939 
    ## ... Procrustes: rmse 0.005376954  max resid 0.01076109 
    ## Run 457 stress 0.02509485 
    ## Run 458 stress 0.01879574 
    ## ... Procrustes: rmse 0.0003346337  max resid 0.0006894294 
    ## ... Similar to previous best
    ## Run 459 stress 0.01883588 
    ## ... Procrustes: rmse 0.004250294  max resid 0.00874733 
    ## ... Similar to previous best
    ## Run 460 stress 0.02509476 
    ## Run 461 stress 0.02486143 
    ## Run 462 stress 0.02520909 
    ## Run 463 stress 0.01882906 
    ## ... Procrustes: rmse 0.001997119  max resid 0.00451071 
    ## ... Similar to previous best
    ## Run 464 stress 0.01879578 
    ## ... Procrustes: rmse 0.001024138  max resid 0.002110483 
    ## ... Similar to previous best
    ## Run 465 stress 0.01892139 
    ## ... Procrustes: rmse 0.007979593  max resid 0.01618497 
    ## Run 466 stress 0.0187958 
    ## ... Procrustes: rmse 0.0003600395  max resid 0.0007416647 
    ## ... Similar to previous best
    ## Run 467 stress 0.018851 
    ## ... Procrustes: rmse 0.005513945  max resid 0.01133947 
    ## Run 468 stress 0.02509461 
    ## Run 469 stress 0.02509464 
    ## Run 470 stress 0.02509456 
    ## Run 471 stress 0.01879565 
    ## ... Procrustes: rmse 0.0002991602  max resid 0.0006163381 
    ## ... Similar to previous best
    ## Run 472 stress 0.01886242 
    ## ... Procrustes: rmse 0.004943823  max resid 0.00985766 
    ## ... Similar to previous best
    ## Run 473 stress 0.01879558 
    ## ... Procrustes: rmse 0.0002643144  max resid 0.0005445257 
    ## ... Similar to previous best
    ## Run 474 stress 0.01879556 
    ## ... Procrustes: rmse 0.0002481466  max resid 0.0005111467 
    ## ... Similar to previous best
    ## Run 475 stress 0.01879584 
    ## ... Procrustes: rmse 0.0003827431  max resid 0.0007885137 
    ## ... Similar to previous best
    ## Run 476 stress 0.02509451 
    ## Run 477 stress 0.02486109 
    ## Run 478 stress 0.02492167 
    ## Run 479 stress 0.02509447 
    ## Run 480 stress 0.02492169 
    ## Run 481 stress 0.01879593 
    ## ... Procrustes: rmse 0.0004195543  max resid 0.0008643108 
    ## ... Similar to previous best
    ## Run 482 stress 0.02492179 
    ## Run 483 stress 0.2038604 
    ## Run 484 stress 0.0250948 
    ## Run 485 stress 0.3652817 
    ## Run 486 stress 0.02509458 
    ## Run 487 stress 0.02492355 
    ## Run 488 stress 0.02509492 
    ## Run 489 stress 0.01892475 
    ## ... Procrustes: rmse 0.006904548  max resid 0.01415627 
    ## Run 490 stress 0.01887843 
    ## ... Procrustes: rmse 0.006546693  max resid 0.0134625 
    ## Run 491 stress 0.02509439 
    ## Run 492 stress 0.02509461 
    ## Run 493 stress 0.01879583 
    ## ... Procrustes: rmse 0.000350488  max resid 0.0007219749 
    ## ... Similar to previous best
    ## Run 494 stress 0.01885562 
    ## ... Procrustes: rmse 0.004601674  max resid 0.009130989 
    ## ... Similar to previous best
    ## Run 495 stress 0.0250946 
    ## Run 496 stress 0.0190279 
    ## ... Procrustes: rmse 0.01150471  max resid 0.02346154 
    ## Run 497 stress 0.01879589 
    ## ... Procrustes: rmse 0.0003976291  max resid 0.0008190142 
    ## ... Similar to previous best
    ## Run 498 stress 0.02509455 
    ## Run 499 stress 0.3654894 
    ## Run 500 stress 0.02486135 
    ## *** Best solution repeated 37 times

``` r
# Mixed and stratified lakes
SD_beta_MS_NMDS <- metaMDS(SD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09286085 
    ## Run 2 stress 0.2758293 
    ## Run 3 stress 0.09407982 
    ## Run 4 stress 0.09908183 
    ## Run 5 stress 0.08973869 
    ## Run 6 stress 0.09030395 
    ## Run 7 stress 0.09337266 
    ## Run 8 stress 0.09286086 
    ## Run 9 stress 0.08440254 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01442928  max resid 0.04294496 
    ## Run 10 stress 0.1052852 
    ## Run 11 stress 0.09030397 
    ## Run 12 stress 0.08503458 
    ## Run 13 stress 0.09030403 
    ## Run 14 stress 0.09168955 
    ## Run 15 stress 0.0997 
    ## Run 16 stress 0.09465905 
    ## Run 17 stress 0.08503551 
    ## Run 18 stress 0.09030395 
    ## Run 19 stress 0.0850374 
    ## Run 20 stress 0.09030393 
    ## Run 21 stress 0.1038043 
    ## Run 22 stress 0.09407978 
    ## Run 23 stress 0.09145329 
    ## Run 24 stress 0.09464463 
    ## Run 25 stress 0.09374249 
    ## Run 26 stress 0.2693957 
    ## Run 27 stress 0.09145344 
    ## Run 28 stress 0.09286094 
    ## Run 29 stress 0.08973898 
    ## Run 30 stress 0.09408006 
    ## Run 31 stress 0.09337258 
    ## Run 32 stress 0.08503509 
    ## Run 33 stress 0.09030397 
    ## Run 34 stress 0.08440256 
    ## ... Procrustes: rmse 4.347531e-05  max resid 7.600404e-05 
    ## ... Similar to previous best
    ## Run 35 stress 0.09416145 
    ## Run 36 stress 0.08773466 
    ## Run 37 stress 0.1038044 
    ## Run 38 stress 0.08973865 
    ## Run 39 stress 0.09464433 
    ## Run 40 stress 0.09030404 
    ## Run 41 stress 0.09030412 
    ## Run 42 stress 0.09030396 
    ## Run 43 stress 0.08503506 
    ## Run 44 stress 0.09969979 
    ## Run 45 stress 0.09030404 
    ## Run 46 stress 0.09445464 
    ## Run 47 stress 0.0916894 
    ## Run 48 stress 0.09268408 
    ## Run 49 stress 0.09030395 
    ## Run 50 stress 0.08440261 
    ## ... Procrustes: rmse 9.067834e-05  max resid 0.0001650576 
    ## ... Similar to previous best
    ## Run 51 stress 0.09030395 
    ## Run 52 stress 0.1038044 
    ## Run 53 stress 0.2567615 
    ## Run 54 stress 0.09030404 
    ## Run 55 stress 0.09030402 
    ## Run 56 stress 0.0915909 
    ## Run 57 stress 0.08503482 
    ## Run 58 stress 0.09168933 
    ## Run 59 stress 0.1052849 
    ## Run 60 stress 0.09159083 
    ## Run 61 stress 0.09159094 
    ## Run 62 stress 0.09969979 
    ## Run 63 stress 0.09374308 
    ## Run 64 stress 0.09145337 
    ## Run 65 stress 0.09381845 
    ## Run 66 stress 0.08773488 
    ## Run 67 stress 0.09380579 
    ## Run 68 stress 0.09159084 
    ## Run 69 stress 0.09168947 
    ## Run 70 stress 0.09465906 
    ## Run 71 stress 0.09286083 
    ## Run 72 stress 0.09465914 
    ## Run 73 stress 0.08973891 
    ## Run 74 stress 0.09539192 
    ## Run 75 stress 0.09145321 
    ## Run 76 stress 0.0926842 
    ## Run 77 stress 0.09535615 
    ## Run 78 stress 0.08503517 
    ## Run 79 stress 0.09465909 
    ## Run 80 stress 0.09374282 
    ## Run 81 stress 0.09286084 
    ## Run 82 stress 0.09969965 
    ## Run 83 stress 0.08503459 
    ## Run 84 stress 0.08503479 
    ## Run 85 stress 0.1038052 
    ## Run 86 stress 0.09407966 
    ## Run 87 stress 0.09403435 
    ## Run 88 stress 0.08440255 
    ## ... Procrustes: rmse 4.215701e-05  max resid 0.0001025503 
    ## ... Similar to previous best
    ## Run 89 stress 0.09030393 
    ## Run 90 stress 0.0933726 
    ## Run 91 stress 0.09286091 
    ## Run 92 stress 0.09400516 
    ## Run 93 stress 0.09407981 
    ## Run 94 stress 0.08503459 
    ## Run 95 stress 0.09030394 
    ## Run 96 stress 0.3053468 
    ## Run 97 stress 0.09969965 
    ## Run 98 stress 0.09145323 
    ## Run 99 stress 0.08503495 
    ## Run 100 stress 0.09268358 
    ## Run 101 stress 0.09030407 
    ## Run 102 stress 0.0844026 
    ## ... Procrustes: rmse 8.820229e-05  max resid 0.000184348 
    ## ... Similar to previous best
    ## Run 103 stress 0.1038046 
    ## Run 104 stress 0.09337263 
    ## Run 105 stress 0.09407968 
    ## Run 106 stress 0.09381846 
    ## Run 107 stress 0.09030404 
    ## Run 108 stress 0.09464467 
    ## Run 109 stress 0.1052851 
    ## Run 110 stress 0.09407995 
    ## Run 111 stress 0.08973865 
    ## Run 112 stress 0.1005995 
    ## Run 113 stress 0.09760818 
    ## Run 114 stress 0.09445476 
    ## Run 115 stress 0.08440271 
    ## ... Procrustes: rmse 0.000169778  max resid 0.0003557725 
    ## ... Similar to previous best
    ## Run 116 stress 0.1026215 
    ## Run 117 stress 0.09535456 
    ## Run 118 stress 0.105285 
    ## Run 119 stress 0.09447311 
    ## Run 120 stress 0.08973892 
    ## Run 121 stress 0.08503614 
    ## Run 122 stress 0.09669846 
    ## Run 123 stress 0.09408 
    ## Run 124 stress 0.1052849 
    ## Run 125 stress 0.09407105 
    ## Run 126 stress 0.09030396 
    ## Run 127 stress 0.0877347 
    ## Run 128 stress 0.09535584 
    ## Run 129 stress 0.09400535 
    ## Run 130 stress 0.09030397 
    ## Run 131 stress 0.09030414 
    ## Run 132 stress 0.09286084 
    ## Run 133 stress 0.09374346 
    ## Run 134 stress 0.09286091 
    ## Run 135 stress 0.09407976 
    ## Run 136 stress 0.1052852 
    ## Run 137 stress 0.09337234 
    ## Run 138 stress 0.09286083 
    ## Run 139 stress 0.08440255 
    ## ... Procrustes: rmse 4.698414e-05  max resid 0.0001057374 
    ## ... Similar to previous best
    ## Run 140 stress 0.08973864 
    ## Run 141 stress 0.08440252 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001605586  max resid 0.0003014966 
    ## ... Similar to previous best
    ## Run 142 stress 0.09381851 
    ## Run 143 stress 0.08973878 
    ## Run 144 stress 0.08503562 
    ## Run 145 stress 0.08503492 
    ## Run 146 stress 0.09465908 
    ## Run 147 stress 0.0966988 
    ## Run 148 stress 0.09407972 
    ## Run 149 stress 0.09374256 
    ## Run 150 stress 0.09030395 
    ## Run 151 stress 0.09407965 
    ## Run 152 stress 0.09145322 
    ## Run 153 stress 0.09168926 
    ## Run 154 stress 0.1052849 
    ## Run 155 stress 0.09417768 
    ## Run 156 stress 0.08773478 
    ## Run 157 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 7.364682e-05  max resid 0.000134539 
    ## ... Similar to previous best
    ## Run 158 stress 0.09030395 
    ## Run 159 stress 0.09969977 
    ## Run 160 stress 0.09760814 
    ## Run 161 stress 0.09030423 
    ## Run 162 stress 0.0926838 
    ## Run 163 stress 0.0976082 
    ## Run 164 stress 0.09535439 
    ## Run 165 stress 0.092684 
    ## Run 166 stress 0.09969965 
    ## Run 167 stress 0.09145336 
    ## Run 168 stress 0.100461 
    ## Run 169 stress 0.09159084 
    ## Run 170 stress 0.09374177 
    ## Run 171 stress 0.09381864 
    ## Run 172 stress 0.09337238 
    ## Run 173 stress 0.08773475 
    ## Run 174 stress 0.09407991 
    ## Run 175 stress 0.09268391 
    ## Run 176 stress 0.09374223 
    ## Run 177 stress 0.08503467 
    ## Run 178 stress 0.09286083 
    ## Run 179 stress 0.09268355 
    ## Run 180 stress 0.08503551 
    ## Run 181 stress 0.09030408 
    ## Run 182 stress 0.09465905 
    ## Run 183 stress 0.09030403 
    ## Run 184 stress 0.09381845 
    ## Run 185 stress 0.1038044 
    ## Run 186 stress 0.08440265 
    ## ... Procrustes: rmse 0.0001460882  max resid 0.00030532 
    ## ... Similar to previous best
    ## Run 187 stress 0.09268342 
    ## Run 188 stress 0.09308968 
    ## Run 189 stress 0.09308968 
    ## Run 190 stress 0.09465916 
    ## Run 191 stress 0.09908182 
    ## Run 192 stress 0.08773465 
    ## Run 193 stress 0.08503468 
    ## Run 194 stress 0.09268379 
    ## Run 195 stress 0.08773469 
    ## Run 196 stress 0.09168926 
    ## Run 197 stress 0.09030393 
    ## Run 198 stress 0.09416134 
    ## Run 199 stress 0.08503473 
    ## Run 200 stress 0.09407965 
    ## Run 201 stress 0.09535474 
    ## Run 202 stress 0.08440251 
    ## ... Procrustes: rmse 1.048306e-05  max resid 2.936114e-05 
    ## ... Similar to previous best
    ## Run 203 stress 0.08973878 
    ## Run 204 stress 0.09268371 
    ## Run 205 stress 0.09445499 
    ## Run 206 stress 0.09286085 
    ## Run 207 stress 0.09416139 
    ## Run 208 stress 0.09407124 
    ## Run 209 stress 0.09465904 
    ## Run 210 stress 0.09403426 
    ## Run 211 stress 0.09337249 
    ## Run 212 stress 0.09445463 
    ## Run 213 stress 0.08503654 
    ## Run 214 stress 0.09416132 
    ## Run 215 stress 0.0940799 
    ## Run 216 stress 0.09407969 
    ## Run 217 stress 0.09168947 
    ## Run 218 stress 0.1026216 
    ## Run 219 stress 0.09168935 
    ## Run 220 stress 0.0959065 
    ## Run 221 stress 0.09416137 
    ## Run 222 stress 0.1103049 
    ## Run 223 stress 0.09030404 
    ## Run 224 stress 0.09159087 
    ## Run 225 stress 0.0915909 
    ## Run 226 stress 0.08973871 
    ## Run 227 stress 0.09381852 
    ## Run 228 stress 0.2465446 
    ## Run 229 stress 0.09465905 
    ## Run 230 stress 0.0937424 
    ## Run 231 stress 0.09760825 
    ## Run 232 stress 0.09337271 
    ## Run 233 stress 0.09337272 
    ## Run 234 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001613527  max resid 0.0003101183 
    ## ... Similar to previous best
    ## Run 235 stress 0.09407979 
    ## Run 236 stress 0.09381844 
    ## Run 237 stress 0.0850348 
    ## Run 238 stress 0.09969997 
    ## Run 239 stress 0.09030396 
    ## Run 240 stress 0.08503574 
    ## Run 241 stress 0.09268331 
    ## Run 242 stress 0.09408016 
    ## Run 243 stress 0.09970005 
    ## Run 244 stress 0.08503628 
    ## Run 245 stress 0.09465905 
    ## Run 246 stress 0.090304 
    ## Run 247 stress 0.08503493 
    ## Run 248 stress 0.08440277 
    ## ... Procrustes: rmse 0.000287373  max resid 0.0005205855 
    ## ... Similar to previous best
    ## Run 249 stress 0.09030399 
    ## Run 250 stress 0.09030396 
    ## Run 251 stress 0.08773472 
    ## Run 252 stress 0.09321865 
    ## Run 253 stress 0.2758242 
    ## Run 254 stress 0.08503469 
    ## Run 255 stress 0.09286084 
    ## Run 256 stress 0.09407984 
    ## Run 257 stress 0.09407998 
    ## Run 258 stress 0.0915909 
    ## Run 259 stress 0.09268409 
    ## Run 260 stress 0.08973866 
    ## Run 261 stress 0.08503518 
    ## Run 262 stress 0.09380578 
    ## Run 263 stress 0.09464464 
    ## Run 264 stress 0.09168941 
    ## Run 265 stress 0.09610782 
    ## Run 266 stress 0.09464407 
    ## Run 267 stress 0.08973862 
    ## Run 268 stress 0.0850349 
    ## Run 269 stress 0.09407983 
    ## Run 270 stress 0.09407994 
    ## Run 271 stress 0.08773469 
    ## Run 272 stress 0.09380588 
    ## Run 273 stress 0.09337281 
    ## Run 274 stress 0.09159084 
    ## Run 275 stress 0.09381849 
    ## Run 276 stress 0.08503522 
    ## Run 277 stress 0.09286085 
    ## Run 278 stress 0.09445465 
    ## Run 279 stress 0.09030395 
    ## Run 280 stress 0.09268389 
    ## Run 281 stress 0.09969978 
    ## Run 282 stress 0.09381852 
    ## Run 283 stress 0.09590881 
    ## Run 284 stress 0.08503484 
    ## Run 285 stress 0.09030397 
    ## Run 286 stress 0.09374365 
    ## Run 287 stress 0.09337282 
    ## Run 288 stress 0.08773479 
    ## Run 289 stress 0.09159083 
    ## Run 290 stress 0.08503602 
    ## Run 291 stress 0.1072862 
    ## Run 292 stress 0.09159092 
    ## Run 293 stress 0.08973885 
    ## Run 294 stress 0.08503498 
    ## Run 295 stress 0.08503607 
    ## Run 296 stress 0.08773484 
    ## Run 297 stress 0.08440251 
    ## ... Procrustes: rmse 2.650668e-05  max resid 5.268337e-05 
    ## ... Similar to previous best
    ## Run 298 stress 0.09030395 
    ## Run 299 stress 0.097212 
    ## Run 300 stress 0.09445467 
    ## Run 301 stress 0.09286087 
    ## Run 302 stress 0.08503456 
    ## Run 303 stress 0.09159095 
    ## Run 304 stress 0.09539199 
    ## Run 305 stress 0.09464471 
    ## Run 306 stress 0.1026217 
    ## Run 307 stress 0.09539209 
    ## Run 308 stress 0.09408013 
    ## Run 309 stress 0.09969984 
    ## Run 310 stress 0.0937436 
    ## Run 311 stress 0.09145325 
    ## Run 312 stress 0.09321766 
    ## Run 313 stress 0.0844026 
    ## ... Procrustes: rmse 0.0001626222  max resid 0.0002838775 
    ## ... Similar to previous best
    ## Run 314 stress 0.08440274 
    ## ... Procrustes: rmse 0.0002486255  max resid 0.0004539454 
    ## ... Similar to previous best
    ## Run 315 stress 0.09416131 
    ## Run 316 stress 0.09407988 
    ## Run 317 stress 0.08503456 
    ## Run 318 stress 0.0844026 
    ## ... Procrustes: rmse 0.000192804  max resid 0.0003837267 
    ## ... Similar to previous best
    ## Run 319 stress 0.275825 
    ## Run 320 stress 0.09030401 
    ## Run 321 stress 0.09145322 
    ## Run 322 stress 0.09030404 
    ## Run 323 stress 0.09030395 
    ## Run 324 stress 0.08773467 
    ## Run 325 stress 0.08440253 
    ## ... Procrustes: rmse 6.993842e-05  max resid 0.0001357005 
    ## ... Similar to previous best
    ## Run 326 stress 0.08503503 
    ## Run 327 stress 0.09286084 
    ## Run 328 stress 0.09145322 
    ## Run 329 stress 0.09337229 
    ## Run 330 stress 0.09159085 
    ## Run 331 stress 0.08503456 
    ## Run 332 stress 0.09447305 
    ## Run 333 stress 0.0938058 
    ## Run 334 stress 0.09407405 
    ## Run 335 stress 0.09268359 
    ## Run 336 stress 0.08440251 
    ## ... Procrustes: rmse 2.325016e-05  max resid 6.00075e-05 
    ## ... Similar to previous best
    ## Run 337 stress 0.0941615 
    ## Run 338 stress 0.09030398 
    ## Run 339 stress 0.09268364 
    ## Run 340 stress 0.08973862 
    ## Run 341 stress 0.09268426 
    ## Run 342 stress 0.09145322 
    ## Run 343 stress 0.08503612 
    ## Run 344 stress 0.09407968 
    ## Run 345 stress 0.09145326 
    ## Run 346 stress 0.09407144 
    ## Run 347 stress 0.0976084 
    ## Run 348 stress 0.09407985 
    ## Run 349 stress 0.09268337 
    ## Run 350 stress 0.09407972 
    ## Run 351 stress 0.0996997 
    ## Run 352 stress 0.09465908 
    ## Run 353 stress 0.09407987 
    ## Run 354 stress 0.09416153 
    ## Run 355 stress 0.09465905 
    ## Run 356 stress 0.09268375 
    ## Run 357 stress 0.08503579 
    ## Run 358 stress 0.09465905 
    ## Run 359 stress 0.08503579 
    ## Run 360 stress 0.09407974 
    ## Run 361 stress 0.09286087 
    ## Run 362 stress 0.08503534 
    ## Run 363 stress 0.2693675 
    ## Run 364 stress 0.09465911 
    ## Run 365 stress 0.09286084 
    ## Run 366 stress 0.09416148 
    ## Run 367 stress 0.09465904 
    ## Run 368 stress 0.09145321 
    ## Run 369 stress 0.09408007 
    ## Run 370 stress 0.09408005 
    ## Run 371 stress 0.09407996 
    ## Run 372 stress 0.09407983 
    ## Run 373 stress 0.09286093 
    ## Run 374 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001174398  max resid 0.0002197512 
    ## ... Similar to previous best
    ## Run 375 stress 0.08440263 
    ## ... Procrustes: rmse 0.0001881364  max resid 0.0003461793 
    ## ... Similar to previous best
    ## Run 376 stress 0.09268353 
    ## Run 377 stress 0.09159087 
    ## Run 378 stress 0.08440268 
    ## ... Procrustes: rmse 0.0002369098  max resid 0.00043459 
    ## ... Similar to previous best
    ## Run 379 stress 0.09168942 
    ## Run 380 stress 0.0877347 
    ## Run 381 stress 0.09159084 
    ## Run 382 stress 0.0933723 
    ## Run 383 stress 0.09159093 
    ## Run 384 stress 0.09268318 
    ## Run 385 stress 0.0946591 
    ## Run 386 stress 0.09465906 
    ## Run 387 stress 0.0916894 
    ## Run 388 stress 0.3251303 
    ## Run 389 stress 0.09286085 
    ## Run 390 stress 0.09145322 
    ## Run 391 stress 0.09159083 
    ## Run 392 stress 0.09969975 
    ## Run 393 stress 0.0897387 
    ## Run 394 stress 0.09030393 
    ## Run 395 stress 0.09590634 
    ## Run 396 stress 0.09407979 
    ## Run 397 stress 0.09535484 
    ## Run 398 stress 0.09030394 
    ## Run 399 stress 0.09416131 
    ## Run 400 stress 0.09969985 
    ## Run 401 stress 0.09030393 
    ## Run 402 stress 0.09030395 
    ## Run 403 stress 0.09416145 
    ## Run 404 stress 0.08503486 
    ## Run 405 stress 0.09030396 
    ## Run 406 stress 0.09159098 
    ## Run 407 stress 0.09030406 
    ## Run 408 stress 0.08973886 
    ## Run 409 stress 0.09407996 
    ## Run 410 stress 0.09381847 
    ## Run 411 stress 0.08440251 
    ## ... Procrustes: rmse 2.476219e-05  max resid 5.143501e-05 
    ## ... Similar to previous best
    ## Run 412 stress 0.09286088 
    ## Run 413 stress 0.1026215 
    ## Run 414 stress 0.09030397 
    ## Run 415 stress 0.09969977 
    ## Run 416 stress 0.09159087 
    ## Run 417 stress 0.08503595 
    ## Run 418 stress 0.09535416 
    ## Run 419 stress 0.09168947 
    ## Run 420 stress 0.09321829 
    ## Run 421 stress 0.09268322 
    ## Run 422 stress 0.09464429 
    ## Run 423 stress 0.09145322 
    ## Run 424 stress 0.09337262 
    ## Run 425 stress 0.09400517 
    ## Run 426 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001031956  max resid 0.0001979214 
    ## ... Similar to previous best
    ## Run 427 stress 0.09145327 
    ## Run 428 stress 0.0914533 
    ## Run 429 stress 0.08440255 
    ## ... Procrustes: rmse 0.000112341  max resid 0.0002156147 
    ## ... Similar to previous best
    ## Run 430 stress 0.09030394 
    ## Run 431 stress 0.08440267 
    ## ... Procrustes: rmse 0.0002310236  max resid 0.0004244096 
    ## ... Similar to previous best
    ## Run 432 stress 0.09145324 
    ## Run 433 stress 0.09407119 
    ## Run 434 stress 0.09159086 
    ## Run 435 stress 0.09407975 
    ## Run 436 stress 0.09374313 
    ## Run 437 stress 0.1026215 
    ## Run 438 stress 0.08440273 
    ## ... Procrustes: rmse 0.0001961489  max resid 0.0003669491 
    ## ... Similar to previous best
    ## Run 439 stress 0.09374232 
    ## Run 440 stress 0.09400517 
    ## Run 441 stress 0.09535469 
    ## Run 442 stress 0.09159083 
    ## Run 443 stress 0.09380591 
    ## Run 444 stress 0.2547731 
    ## Run 445 stress 0.09145325 
    ## Run 446 stress 0.08773465 
    ## Run 447 stress 0.09465907 
    ## Run 448 stress 0.09381862 
    ## Run 449 stress 0.09030403 
    ## Run 450 stress 0.08503511 
    ## Run 451 stress 0.09407309 
    ## Run 452 stress 0.09145326 
    ## Run 453 stress 0.09286096 
    ## Run 454 stress 0.09374317 
    ## Run 455 stress 0.09400537 
    ## Run 456 stress 0.09337236 
    ## Run 457 stress 0.100461 
    ## Run 458 stress 0.100461 
    ## Run 459 stress 0.100461 
    ## Run 460 stress 0.1038047 
    ## Run 461 stress 0.09337255 
    ## Run 462 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001524321  max resid 0.0003155118 
    ## ... Similar to previous best
    ## Run 463 stress 0.08973897 
    ## Run 464 stress 0.08503569 
    ## Run 465 stress 0.08973888 
    ## Run 466 stress 0.1005995 
    ## Run 467 stress 0.08973871 
    ## Run 468 stress 0.09159084 
    ## Run 469 stress 0.09970006 
    ## Run 470 stress 0.08503564 
    ## Run 471 stress 0.09969968 
    ## Run 472 stress 0.08503638 
    ## Run 473 stress 0.08773486 
    ## Run 474 stress 0.08503488 
    ## Run 475 stress 0.08503505 
    ## Run 476 stress 0.09374349 
    ## Run 477 stress 0.08973886 
    ## Run 478 stress 0.09380581 
    ## Run 479 stress 0.09445481 
    ## Run 480 stress 0.09407107 
    ## Run 481 stress 0.09445497 
    ## Run 482 stress 0.08503491 
    ## Run 483 stress 0.09721197 
    ## Run 484 stress 0.09407121 
    ## Run 485 stress 0.09374278 
    ## Run 486 stress 0.08773467 
    ## Run 487 stress 0.08503492 
    ## Run 488 stress 0.08503541 
    ## Run 489 stress 0.0903041 
    ## Run 490 stress 0.09268369 
    ## Run 491 stress 0.08503461 
    ## Run 492 stress 0.09407405 
    ## Run 493 stress 0.09590857 
    ## Run 494 stress 0.09969994 
    ## Run 495 stress 0.08503462 
    ## Run 496 stress 0.08440264 
    ## ... Procrustes: rmse 0.000184617  max resid 0.0003895914 
    ## ... Similar to previous best
    ## Run 497 stress 0.0897387 
    ## Run 498 stress 0.09268326 
    ## Run 499 stress 0.09308953 
    ## Run 500 stress 0.09286086 
    ## *** Best solution repeated 21 times

``` r
# Ocean sites and mixed lakes
SD_beta_OM_NMDS <- metaMDS(SD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743114  max resid 0.05098619 
    ## Run 2 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744206  max resid 0.05102372 
    ## Run 3 stress 0.08233762 
    ## Run 4 stress 0.07732931 
    ## Run 5 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745075  max resid 0.05108304 
    ## Run 6 stress 0.07732929 
    ## Run 7 stress 0.07629249 
    ## Run 8 stress 0.08000474 
    ## Run 9 stress 0.0800047 
    ## Run 10 stress 0.08000476 
    ## Run 11 stress 0.07629235 
    ## Run 12 stress 0.07629235 
    ## Run 13 stress 0.07629234 
    ## Run 14 stress 0.349573 
    ## Run 15 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744385  max resid 0.05103864 
    ## Run 16 stress 0.08233754 
    ## Run 17 stress 0.08000467 
    ## Run 18 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744483  max resid 0.05101679 
    ## Run 19 stress 0.07629234 
    ## Run 20 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744078  max resid 0.05105154 
    ## Run 21 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744714  max resid 0.05106567 
    ## Run 22 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744158  max resid 0.05100054 
    ## Run 23 stress 0.07732933 
    ## Run 24 stress 0.07378228 
    ## ... Procrustes: rmse 0.01743993  max resid 0.05100373 
    ## Run 25 stress 0.08000474 
    ## Run 26 stress 0.07365789 
    ## ... Procrustes: rmse 9.354631e-05  max resid 0.0002236235 
    ## ... Similar to previous best
    ## Run 27 stress 0.07378229 
    ## ... Procrustes: rmse 0.01746868  max resid 0.05115694 
    ## Run 28 stress 0.07365785 
    ## ... Procrustes: rmse 8.399723e-05  max resid 0.0001968643 
    ## ... Similar to previous best
    ## Run 29 stress 0.0736579 
    ## ... Procrustes: rmse 0.000121612  max resid 0.0002908531 
    ## ... Similar to previous best
    ## Run 30 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174476  max resid 0.05103826 
    ## Run 31 stress 0.07365785 
    ## ... Procrustes: rmse 8.05607e-05  max resid 0.0001883294 
    ## ... Similar to previous best
    ## Run 32 stress 0.07629239 
    ## Run 33 stress 0.07365784 
    ## ... Procrustes: rmse 4.748744e-05  max resid 0.0001117472 
    ## ... Similar to previous best
    ## Run 34 stress 0.07365785 
    ## ... Procrustes: rmse 7.085728e-05  max resid 0.0001657926 
    ## ... Similar to previous best
    ## Run 35 stress 0.07365784 
    ## ... Procrustes: rmse 2.614674e-05  max resid 6.152397e-05 
    ## ... Similar to previous best
    ## Run 36 stress 0.08233767 
    ## Run 37 stress 0.07365784 
    ## ... Procrustes: rmse 6.602508e-05  max resid 0.0001557514 
    ## ... Similar to previous best
    ## Run 38 stress 0.07732934 
    ## Run 39 stress 0.07629237 
    ## Run 40 stress 0.08000468 
    ## Run 41 stress 0.07365784 
    ## ... Procrustes: rmse 4.062358e-05  max resid 9.356882e-05 
    ## ... Similar to previous best
    ## Run 42 stress 0.08233758 
    ## Run 43 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743425  max resid 0.05099933 
    ## Run 44 stress 0.07629235 
    ## Run 45 stress 0.07629242 
    ## Run 46 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174558  max resid 0.05107198 
    ## Run 47 stress 0.07629233 
    ## Run 48 stress 0.07629235 
    ## Run 49 stress 0.08000469 
    ## Run 50 stress 0.07365786 
    ## ... Procrustes: rmse 8.24342e-05  max resid 0.0001956162 
    ## ... Similar to previous best
    ## Run 51 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744686  max resid 0.05105099 
    ## Run 52 stress 0.07365784 
    ## ... Procrustes: rmse 4.570873e-05  max resid 0.0001046064 
    ## ... Similar to previous best
    ## Run 53 stress 0.07629245 
    ## Run 54 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744267  max resid 0.05102328 
    ## Run 55 stress 0.07629234 
    ## Run 56 stress 0.08233767 
    ## Run 57 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743922  max resid 0.05101926 
    ## Run 58 stress 0.07365787 
    ## ... Procrustes: rmse 2.688051e-05  max resid 5.008701e-05 
    ## ... Similar to previous best
    ## Run 59 stress 0.08233767 
    ## Run 60 stress 0.07365784 
    ## ... Procrustes: rmse 6.198014e-05  max resid 0.0001444725 
    ## ... Similar to previous best
    ## Run 61 stress 0.08288057 
    ## Run 62 stress 0.07365784 
    ## ... Procrustes: rmse 2.481607e-05  max resid 5.57673e-05 
    ## ... Similar to previous best
    ## Run 63 stress 0.07732933 
    ## Run 64 stress 0.07629245 
    ## Run 65 stress 0.07629238 
    ## Run 66 stress 0.07365784 
    ## ... Procrustes: rmse 5.15912e-05  max resid 0.0001195973 
    ## ... Similar to previous best
    ## Run 67 stress 0.07629242 
    ## Run 68 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744459  max resid 0.05104956 
    ## Run 69 stress 0.0737823 
    ## ... Procrustes: rmse 0.01743462  max resid 0.05097135 
    ## Run 70 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001391153  max resid 0.0003246924 
    ## ... Similar to previous best
    ## Run 71 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743884  max resid 0.05102448 
    ## Run 72 stress 0.0823376 
    ## Run 73 stress 0.08000469 
    ## Run 74 stress 0.07629235 
    ## Run 75 stress 0.07732931 
    ## Run 76 stress 0.2754058 
    ## Run 77 stress 0.07365784 
    ## ... Procrustes: rmse 3.691563e-05  max resid 8.859721e-05 
    ## ... Similar to previous best
    ## Run 78 stress 0.08288047 
    ## Run 79 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744442  max resid 0.05105382 
    ## Run 80 stress 0.07365794 
    ## ... Procrustes: rmse 0.0001679254  max resid 0.0003973572 
    ## ... Similar to previous best
    ## Run 81 stress 0.07629233 
    ## Run 82 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744368  max resid 0.0510403 
    ## Run 83 stress 0.07365794 
    ## ... Procrustes: rmse 0.0001150213  max resid 0.0002739468 
    ## ... Similar to previous best
    ## Run 84 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001137489  max resid 0.0002645219 
    ## ... Similar to previous best
    ## Run 85 stress 0.08000476 
    ## Run 86 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744889  max resid 0.05105737 
    ## Run 87 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001599893  max resid 0.000378033 
    ## ... Similar to previous best
    ## Run 88 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744502  max resid 0.05104812 
    ## Run 89 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001233429  max resid 0.0002897722 
    ## ... Similar to previous best
    ## Run 90 stress 0.2560005 
    ## Run 91 stress 0.08000475 
    ## Run 92 stress 0.07365786 
    ## ... Procrustes: rmse 8.348028e-05  max resid 0.0001962254 
    ## ... Similar to previous best
    ## Run 93 stress 0.07732937 
    ## Run 94 stress 0.07732932 
    ## Run 95 stress 0.08000478 
    ## Run 96 stress 0.07378229 
    ## ... Procrustes: rmse 0.01747384  max resid 0.05115771 
    ## Run 97 stress 0.08288059 
    ## Run 98 stress 0.07629248 
    ## Run 99 stress 0.07629237 
    ## Run 100 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744802  max resid 0.05107336 
    ## Run 101 stress 0.07732929 
    ## Run 102 stress 0.08233756 
    ## Run 103 stress 0.07629235 
    ## Run 104 stress 0.08000467 
    ## Run 105 stress 0.07365786 
    ## ... Procrustes: rmse 9.829462e-05  max resid 0.0002321142 
    ## ... Similar to previous best
    ## Run 106 stress 0.07365783 
    ## ... Procrustes: rmse 1.196179e-05  max resid 2.298857e-05 
    ## ... Similar to previous best
    ## Run 107 stress 0.07732938 
    ## Run 108 stress 0.08233752 
    ## Run 109 stress 0.07365784 
    ## ... Procrustes: rmse 4.550303e-05  max resid 0.0001077889 
    ## ... Similar to previous best
    ## Run 110 stress 0.07365785 
    ## ... Procrustes: rmse 7.100802e-05  max resid 0.0001679897 
    ## ... Similar to previous best
    ## Run 111 stress 0.07629236 
    ## Run 112 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745266  max resid 0.05108734 
    ## Run 113 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744326  max resid 0.05103389 
    ## Run 114 stress 0.07629235 
    ## Run 115 stress 0.07365792 
    ## ... Procrustes: rmse 0.000162159  max resid 0.0003772498 
    ## ... Similar to previous best
    ## Run 116 stress 0.08288034 
    ## Run 117 stress 0.08000468 
    ## Run 118 stress 0.07629235 
    ## Run 119 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744634  max resid 0.0510526 
    ## Run 120 stress 0.07378229 
    ## ... Procrustes: rmse 0.01747311  max resid 0.05117154 
    ## Run 121 stress 0.07629237 
    ## Run 122 stress 0.07732937 
    ## Run 123 stress 0.0737823 
    ## ... Procrustes: rmse 0.01743676  max resid 0.05104109 
    ## Run 124 stress 0.07365787 
    ## ... Procrustes: rmse 6.268315e-05  max resid 0.0001493224 
    ## ... Similar to previous best
    ## Run 125 stress 0.07365784 
    ## ... Procrustes: rmse 3.761201e-05  max resid 8.097459e-05 
    ## ... Similar to previous best
    ## Run 126 stress 0.08000468 
    ## Run 127 stress 0.08000469 
    ## Run 128 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001088476  max resid 0.000250682 
    ## ... Similar to previous best
    ## Run 129 stress 0.07629238 
    ## Run 130 stress 0.08233776 
    ## Run 131 stress 0.08288049 
    ## Run 132 stress 0.08000468 
    ## Run 133 stress 0.07629234 
    ## Run 134 stress 0.08000477 
    ## Run 135 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744502  max resid 0.05105109 
    ## Run 136 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744261  max resid 0.05102934 
    ## Run 137 stress 0.3495699 
    ## Run 138 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001304235  max resid 0.0003084238 
    ## ... Similar to previous best
    ## Run 139 stress 0.07629233 
    ## Run 140 stress 0.08000472 
    ## Run 141 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744461  max resid 0.0510458 
    ## Run 142 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744147  max resid 0.0510219 
    ## Run 143 stress 0.07365783 
    ## ... Procrustes: rmse 9.844158e-06  max resid 2.353012e-05 
    ## ... Similar to previous best
    ## Run 144 stress 0.0800047 
    ## Run 145 stress 0.07629238 
    ## Run 146 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744623  max resid 0.05105268 
    ## Run 147 stress 0.07365785 
    ## ... Procrustes: rmse 6.931168e-05  max resid 0.0001605887 
    ## ... Similar to previous best
    ## Run 148 stress 0.07629234 
    ## Run 149 stress 0.0800047 
    ## Run 150 stress 0.07732935 
    ## Run 151 stress 0.07732933 
    ## Run 152 stress 0.07365785 
    ## ... Procrustes: rmse 7.326411e-05  max resid 0.0001718585 
    ## ... Similar to previous best
    ## Run 153 stress 0.07365783 
    ## ... Procrustes: rmse 8.753953e-06  max resid 2.011783e-05 
    ## ... Similar to previous best
    ## Run 154 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744123  max resid 0.0510521 
    ## Run 155 stress 0.08233756 
    ## Run 156 stress 0.08000468 
    ## Run 157 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174443  max resid 0.05103945 
    ## Run 158 stress 0.07629236 
    ## Run 159 stress 0.07629235 
    ## Run 160 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744461  max resid 0.05104331 
    ## Run 161 stress 0.0762924 
    ## Run 162 stress 0.07629236 
    ## Run 163 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174486  max resid 0.05106313 
    ## Run 164 stress 0.07629234 
    ## Run 165 stress 0.07629244 
    ## Run 166 stress 0.07378232 
    ## ... Procrustes: rmse 0.01744429  max resid 0.05107494 
    ## Run 167 stress 0.07629236 
    ## Run 168 stress 0.07365783 
    ## ... Procrustes: rmse 2.627841e-05  max resid 6.246357e-05 
    ## ... Similar to previous best
    ## Run 169 stress 0.07365787 
    ## ... Procrustes: rmse 8.272767e-05  max resid 0.0001979281 
    ## ... Similar to previous best
    ## Run 170 stress 0.07365788 
    ## ... Procrustes: rmse 0.000103836  max resid 0.0002393831 
    ## ... Similar to previous best
    ## Run 171 stress 0.07629233 
    ## Run 172 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744145  max resid 0.0510471 
    ## Run 173 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744795  max resid 0.05105862 
    ## Run 174 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744117  max resid 0.05101676 
    ## Run 175 stress 0.07365784 
    ## ... Procrustes: rmse 3.073812e-05  max resid 7.308292e-05 
    ## ... Similar to previous best
    ## Run 176 stress 0.0800047 
    ## Run 177 stress 0.08000469 
    ## Run 178 stress 0.07629239 
    ## Run 179 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744526  max resid 0.05104762 
    ## Run 180 stress 0.07365783 
    ## ... Procrustes: rmse 1.0822e-05  max resid 2.490801e-05 
    ## ... Similar to previous best
    ## Run 181 stress 0.07732927 
    ## Run 182 stress 0.07629234 
    ## Run 183 stress 0.08000472 
    ## Run 184 stress 0.07365783 
    ## ... Procrustes: rmse 2.081189e-05  max resid 4.877966e-05 
    ## ... Similar to previous best
    ## Run 185 stress 0.08233756 
    ## Run 186 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001317483  max resid 0.0003146076 
    ## ... Similar to previous best
    ## Run 187 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745814  max resid 0.05109747 
    ## Run 188 stress 0.07629246 
    ## Run 189 stress 0.07365783 
    ## ... Procrustes: rmse 2.955517e-05  max resid 6.902474e-05 
    ## ... Similar to previous best
    ## Run 190 stress 0.07629233 
    ## Run 191 stress 0.0828805 
    ## Run 192 stress 0.07629234 
    ## Run 193 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745192  max resid 0.05108977 
    ## Run 194 stress 0.07629234 
    ## Run 195 stress 0.07365785 
    ## ... Procrustes: rmse 3.016233e-05  max resid 5.55925e-05 
    ## ... Similar to previous best
    ## Run 196 stress 0.07365784 
    ## ... Procrustes: rmse 5.114894e-05  max resid 0.0001193947 
    ## ... Similar to previous best
    ## Run 197 stress 0.0773293 
    ## Run 198 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744652  max resid 0.05102421 
    ## Run 199 stress 0.07378227 
    ## ... Procrustes: rmse 0.01742778  max resid 0.05097301 
    ## Run 200 stress 0.08288047 
    ## Run 201 stress 0.07365785 
    ## ... Procrustes: rmse 6.399685e-05  max resid 0.0001511885 
    ## ... Similar to previous best
    ## Run 202 stress 0.07378233 
    ## ... Procrustes: rmse 0.0174368  max resid 0.05104721 
    ## Run 203 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744453  max resid 0.05105344 
    ## Run 204 stress 0.07365783 
    ## ... Procrustes: rmse 3.714807e-05  max resid 8.814352e-05 
    ## ... Similar to previous best
    ## Run 205 stress 0.07365791 
    ## ... Procrustes: rmse 0.00013614  max resid 0.0003241156 
    ## ... Similar to previous best
    ## Run 206 stress 0.07365783 
    ## ... Procrustes: rmse 4.817075e-06  max resid 1.098544e-05 
    ## ... Similar to previous best
    ## Run 207 stress 0.07365783 
    ## ... Procrustes: rmse 7.320628e-06  max resid 1.602934e-05 
    ## ... Similar to previous best
    ## Run 208 stress 0.08000468 
    ## Run 209 stress 0.08233762 
    ## Run 210 stress 0.07365783 
    ## ... Procrustes: rmse 1.604685e-05  max resid 3.716423e-05 
    ## ... Similar to previous best
    ## Run 211 stress 0.07629241 
    ## Run 212 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744605  max resid 0.05105696 
    ## Run 213 stress 0.08233755 
    ## Run 214 stress 0.07378227 
    ## ... Procrustes: rmse 0.01746525  max resid 0.05112849 
    ## Run 215 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744601  max resid 0.05105381 
    ## Run 216 stress 0.07629234 
    ## Run 217 stress 0.07365783 
    ## ... Procrustes: rmse 3.421832e-06  max resid 6.612473e-06 
    ## ... Similar to previous best
    ## Run 218 stress 0.0737823 
    ## ... Procrustes: rmse 0.01743635  max resid 0.05103751 
    ## Run 219 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744495  max resid 0.05103133 
    ## Run 220 stress 0.07365787 
    ## ... Procrustes: rmse 0.00010771  max resid 0.0002540672 
    ## ... Similar to previous best
    ## Run 221 stress 0.07629234 
    ## Run 222 stress 0.07365785 
    ## ... Procrustes: rmse 7.158136e-05  max resid 0.0001678689 
    ## ... Similar to previous best
    ## Run 223 stress 0.07629237 
    ## Run 224 stress 0.0773293 
    ## Run 225 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743083  max resid 0.05099353 
    ## Run 226 stress 0.07629236 
    ## Run 227 stress 0.07629237 
    ## Run 228 stress 0.07365788 
    ## ... Procrustes: rmse 0.000102816  max resid 0.0002372925 
    ## ... Similar to previous best
    ## Run 229 stress 0.07629234 
    ## Run 230 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001376329  max resid 0.0003229887 
    ## ... Similar to previous best
    ## Run 231 stress 0.07365783 
    ## ... Procrustes: rmse 1.67321e-05  max resid 3.635399e-05 
    ## ... Similar to previous best
    ## Run 232 stress 0.07629238 
    ## Run 233 stress 0.07629239 
    ## Run 234 stress 0.08233752 
    ## Run 235 stress 0.07629238 
    ## Run 236 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744513  max resid 0.05103067 
    ## Run 237 stress 0.07629243 
    ## Run 238 stress 0.07629235 
    ## Run 239 stress 0.08233768 
    ## Run 240 stress 0.07365791 
    ## ... Procrustes: rmse 0.000141029  max resid 0.0003285236 
    ## ... Similar to previous best
    ## Run 241 stress 0.07365784 
    ## ... Procrustes: rmse 4.917218e-05  max resid 0.0001155403 
    ## ... Similar to previous best
    ## Run 242 stress 0.07629243 
    ## Run 243 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745219  max resid 0.05108853 
    ## Run 244 stress 0.07629243 
    ## Run 245 stress 0.07365787 
    ## ... Procrustes: rmse 9.540267e-05  max resid 0.0002250661 
    ## ... Similar to previous best
    ## Run 246 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001255107  max resid 0.0003007531 
    ## ... Similar to previous best
    ## Run 247 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001151293  max resid 0.000272117 
    ## ... Similar to previous best
    ## Run 248 stress 0.08288052 
    ## Run 249 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744194  max resid 0.05105417 
    ## Run 250 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001093014  max resid 0.0002540622 
    ## ... Similar to previous best
    ## Run 251 stress 0.0773294 
    ## Run 252 stress 0.07629237 
    ## Run 253 stress 0.07365783 
    ## ... Procrustes: rmse 4.178458e-06  max resid 6.361297e-06 
    ## ... Similar to previous best
    ## Run 254 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744403  max resid 0.05103589 
    ## Run 255 stress 0.07732937 
    ## Run 256 stress 0.08000467 
    ## Run 257 stress 0.07629246 
    ## Run 258 stress 0.07629246 
    ## Run 259 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001007558  max resid 0.0002350385 
    ## ... Similar to previous best
    ## Run 260 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744559  max resid 0.05105282 
    ## Run 261 stress 0.08000469 
    ## Run 262 stress 0.07629233 
    ## Run 263 stress 0.07629234 
    ## Run 264 stress 0.08000477 
    ## Run 265 stress 0.08000467 
    ## Run 266 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744327  max resid 0.05104395 
    ## Run 267 stress 0.08000467 
    ## Run 268 stress 0.08233759 
    ## Run 269 stress 0.07629235 
    ## Run 270 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744405  max resid 0.05103939 
    ## Run 271 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744552  max resid 0.05105375 
    ## Run 272 stress 0.0800047 
    ## Run 273 stress 0.07629235 
    ## Run 274 stress 0.3582286 
    ## Run 275 stress 0.07629238 
    ## Run 276 stress 0.07365784 
    ## ... Procrustes: rmse 5.126635e-05  max resid 0.0001192673 
    ## ... Similar to previous best
    ## Run 277 stress 0.07365784 
    ## ... Procrustes: rmse 5.458441e-05  max resid 0.0001264828 
    ## ... Similar to previous best
    ## Run 278 stress 0.07365785 
    ## ... Procrustes: rmse 5.929454e-05  max resid 0.0001416554 
    ## ... Similar to previous best
    ## Run 279 stress 0.08000467 
    ## Run 280 stress 0.07629235 
    ## Run 281 stress 0.07629235 
    ## Run 282 stress 0.08000467 
    ## Run 283 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744234  max resid 0.0510064 
    ## Run 284 stress 0.07732938 
    ## Run 285 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744294  max resid 0.05102243 
    ## Run 286 stress 0.0800048 
    ## Run 287 stress 0.07629234 
    ## Run 288 stress 0.326292 
    ## Run 289 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744199  max resid 0.05105211 
    ## Run 290 stress 0.08233762 
    ## Run 291 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744483  max resid 0.05104559 
    ## Run 292 stress 0.08288035 
    ## Run 293 stress 0.08000467 
    ## Run 294 stress 0.07378232 
    ## ... Procrustes: rmse 0.0174362  max resid 0.05104358 
    ## Run 295 stress 0.08233762 
    ## Run 296 stress 0.08000477 
    ## Run 297 stress 0.07629235 
    ## Run 298 stress 0.07365786 
    ## ... Procrustes: rmse 4.388076e-05  max resid 9.485404e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744557  max resid 0.05105411 
    ## Run 300 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174495  max resid 0.05105966 
    ## Run 301 stress 0.07365783 
    ## ... Procrustes: rmse 1.334353e-05  max resid 3.127969e-05 
    ## ... Similar to previous best
    ## Run 302 stress 0.07629249 
    ## Run 303 stress 0.0762924 
    ## Run 304 stress 0.08000477 
    ## Run 305 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744739  max resid 0.05104374 
    ## Run 306 stress 0.08233767 
    ## Run 307 stress 0.0800047 
    ## Run 308 stress 0.07365783 
    ## ... Procrustes: rmse 2.257795e-05  max resid 5.288172e-05 
    ## ... Similar to previous best
    ## Run 309 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744436  max resid 0.05101963 
    ## Run 310 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001124255  max resid 0.0002640302 
    ## ... Similar to previous best
    ## Run 311 stress 0.0800047 
    ## Run 312 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745201  max resid 0.05108958 
    ## Run 313 stress 0.07732931 
    ## Run 314 stress 0.07629253 
    ## Run 315 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744699  max resid 0.05104677 
    ## Run 316 stress 0.0773293 
    ## Run 317 stress 0.07629234 
    ## Run 318 stress 0.07629238 
    ## Run 319 stress 0.07629233 
    ## Run 320 stress 0.08233753 
    ## Run 321 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 1.082061e-05  max resid 2.532412e-05 
    ## ... Similar to previous best
    ## Run 322 stress 0.08000468 
    ## Run 323 stress 0.07629235 
    ## Run 324 stress 0.07365783 
    ## ... Procrustes: rmse 4.290666e-06  max resid 9.5297e-06 
    ## ... Similar to previous best
    ## Run 325 stress 0.08233753 
    ## Run 326 stress 0.07629234 
    ## Run 327 stress 0.07732927 
    ## Run 328 stress 0.07629236 
    ## Run 329 stress 0.0800047 
    ## Run 330 stress 0.07629234 
    ## Run 331 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744771  max resid 0.05105623 
    ## Run 332 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744653  max resid 0.05104977 
    ## Run 333 stress 0.07629236 
    ## Run 334 stress 0.07365795 
    ## ... Procrustes: rmse 0.0001892181  max resid 0.0004474562 
    ## ... Similar to previous best
    ## Run 335 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744852  max resid 0.05105322 
    ## Run 336 stress 0.08233752 
    ## Run 337 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744172  max resid 0.05102107 
    ## Run 338 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744492  max resid 0.05102922 
    ## Run 339 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001557067  max resid 0.0003634109 
    ## ... Similar to previous best
    ## Run 340 stress 0.0800047 
    ## Run 341 stress 0.07365783 
    ## ... Procrustes: rmse 1.620126e-05  max resid 3.554202e-05 
    ## ... Similar to previous best
    ## Run 342 stress 0.07365783 
    ## ... Procrustes: rmse 1.075681e-05  max resid 2.544952e-05 
    ## ... Similar to previous best
    ## Run 343 stress 0.07629244 
    ## Run 344 stress 0.0828804 
    ## Run 345 stress 0.07629236 
    ## Run 346 stress 0.08000467 
    ## Run 347 stress 0.07629237 
    ## Run 348 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744782  max resid 0.05105493 
    ## Run 349 stress 0.07629234 
    ## Run 350 stress 0.08000469 
    ## Run 351 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745383  max resid 0.05106768 
    ## Run 352 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745051  max resid 0.05105499 
    ## Run 353 stress 0.07629237 
    ## Run 354 stress 0.08000468 
    ## Run 355 stress 0.07629242 
    ## Run 356 stress 0.07732933 
    ## Run 357 stress 0.0800047 
    ## Run 358 stress 0.07629234 
    ## Run 359 stress 0.08233768 
    ## Run 360 stress 0.07365783 
    ## ... Procrustes: rmse 3.968991e-05  max resid 9.435377e-05 
    ## ... Similar to previous best
    ## Run 361 stress 0.08000467 
    ## Run 362 stress 0.08288048 
    ## Run 363 stress 0.07365783 
    ## ... Procrustes: rmse 1.626352e-05  max resid 4.024417e-05 
    ## ... Similar to previous best
    ## Run 364 stress 0.07365786 
    ## ... Procrustes: rmse 2.252857e-05  max resid 4.05535e-05 
    ## ... Similar to previous best
    ## Run 365 stress 0.07365783 
    ## ... Procrustes: rmse 2.957601e-05  max resid 6.932923e-05 
    ## ... Similar to previous best
    ## Run 366 stress 0.07629237 
    ## Run 367 stress 0.08000477 
    ## Run 368 stress 0.08233758 
    ## Run 369 stress 0.08233763 
    ## Run 370 stress 0.0773294 
    ## Run 371 stress 0.3582286 
    ## Run 372 stress 0.08000467 
    ## Run 373 stress 0.08000471 
    ## Run 374 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744762  max resid 0.05105763 
    ## Run 375 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 3.783696e-06  max resid 8.414755e-06 
    ## ... Similar to previous best
    ## Run 376 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744624  max resid 0.05104609 
    ## Run 377 stress 0.07365784 
    ## ... Procrustes: rmse 5.648263e-05  max resid 0.0001335349 
    ## ... Similar to previous best
    ## Run 378 stress 0.07732929 
    ## Run 379 stress 0.07732939 
    ## Run 380 stress 0.07629237 
    ## Run 381 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744894  max resid 0.05106651 
    ## Run 382 stress 0.07378228 
    ## ... Procrustes: rmse 0.01746858  max resid 0.05114793 
    ## Run 383 stress 0.08233755 
    ## Run 384 stress 0.07629241 
    ## Run 385 stress 0.07365784 
    ## ... Procrustes: rmse 5.102944e-05  max resid 0.0001201418 
    ## ... Similar to previous best
    ## Run 386 stress 0.07629234 
    ## Run 387 stress 0.07365783 
    ## ... Procrustes: rmse 1.337687e-05  max resid 2.869369e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.07629242 
    ## Run 389 stress 0.0800047 
    ## Run 390 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745363  max resid 0.05107422 
    ## Run 391 stress 0.07629248 
    ## Run 392 stress 0.08000484 
    ## Run 393 stress 0.07365785 
    ## ... Procrustes: rmse 7.590111e-05  max resid 0.0001788488 
    ## ... Similar to previous best
    ## Run 394 stress 0.08233763 
    ## Run 395 stress 0.07732931 
    ## Run 396 stress 0.07378231 
    ## ... Procrustes: rmse 0.01744627  max resid 0.05101694 
    ## Run 397 stress 0.07629239 
    ## Run 398 stress 0.07378231 
    ## ... Procrustes: rmse 0.01744093  max resid 0.0510547 
    ## Run 399 stress 0.08000469 
    ## Run 400 stress 0.07629236 
    ## Run 401 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744382  max resid 0.05105543 
    ## Run 402 stress 0.08000467 
    ## Run 403 stress 0.07629235 
    ## Run 404 stress 0.07732937 
    ## Run 405 stress 0.2754058 
    ## Run 406 stress 0.07732934 
    ## Run 407 stress 0.07732943 
    ## Run 408 stress 0.07365786 
    ## ... Procrustes: rmse 9.110713e-05  max resid 0.0002146077 
    ## ... Similar to previous best
    ## Run 409 stress 0.07365784 
    ## ... Procrustes: rmse 5.501094e-05  max resid 0.0001303719 
    ## ... Similar to previous best
    ## Run 410 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744678  max resid 0.05104935 
    ## Run 411 stress 0.08000468 
    ## Run 412 stress 0.07365783 
    ## ... Procrustes: rmse 8.997249e-06  max resid 1.372949e-05 
    ## ... Similar to previous best
    ## Run 413 stress 0.07629247 
    ## Run 414 stress 0.07629234 
    ## Run 415 stress 0.07365784 
    ## ... Procrustes: rmse 3.275415e-05  max resid 7.632428e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174461  max resid 0.05103558 
    ## Run 417 stress 0.07629246 
    ## Run 418 stress 0.08000467 
    ## Run 419 stress 0.07365784 
    ## ... Procrustes: rmse 6.446638e-05  max resid 0.0001516002 
    ## ... Similar to previous best
    ## Run 420 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744256  max resid 0.05105556 
    ## Run 421 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001399431  max resid 0.0003322173 
    ## ... Similar to previous best
    ## Run 422 stress 0.07732931 
    ## Run 423 stress 0.08000469 
    ## Run 424 stress 0.07378235 
    ## ... Procrustes: rmse 0.01743717  max resid 0.05096253 
    ## Run 425 stress 0.07629236 
    ## Run 426 stress 0.0800047 
    ## Run 427 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001283219  max resid 0.0003029843 
    ## ... Similar to previous best
    ## Run 428 stress 0.07629238 
    ## Run 429 stress 0.07629248 
    ## Run 430 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744288  max resid 0.05103931 
    ## Run 431 stress 0.07365784 
    ## ... Procrustes: rmse 2.303485e-05  max resid 5.429308e-05 
    ## ... Similar to previous best
    ## Run 432 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745972  max resid 0.05112658 
    ## Run 433 stress 0.08288058 
    ## Run 434 stress 0.07365784 
    ## ... Procrustes: rmse 6.161471e-05  max resid 0.0001435415 
    ## ... Similar to previous best
    ## Run 435 stress 0.07365786 
    ## ... Procrustes: rmse 8.566016e-05  max resid 0.0002016644 
    ## ... Similar to previous best
    ## Run 436 stress 0.07365784 
    ## ... Procrustes: rmse 5.549308e-05  max resid 0.0001327881 
    ## ... Similar to previous best
    ## Run 437 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743964  max resid 0.05101698 
    ## Run 438 stress 0.0762924 
    ## Run 439 stress 0.07378229 
    ## ... Procrustes: rmse 0.01743977  max resid 0.05101721 
    ## Run 440 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174313  max resid 0.05097214 
    ## Run 441 stress 0.07629246 
    ## Run 442 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743808  max resid 0.05101097 
    ## Run 443 stress 0.07378226 
    ## ... Procrustes: rmse 0.017454  max resid 0.05107507 
    ## Run 444 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744986  max resid 0.05105179 
    ## Run 445 stress 0.08000467 
    ## Run 446 stress 0.07378233 
    ## ... Procrustes: rmse 0.0174363  max resid 0.05104438 
    ## Run 447 stress 0.07732937 
    ## Run 448 stress 0.07629235 
    ## Run 449 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745157  max resid 0.05106274 
    ## Run 450 stress 0.07629245 
    ## Run 451 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744653  max resid 0.05104494 
    ## Run 452 stress 0.07732926 
    ## Run 453 stress 0.08233769 
    ## Run 454 stress 0.07629242 
    ## Run 455 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744411  max resid 0.05106026 
    ## Run 456 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001037333  max resid 0.0002435506 
    ## ... Similar to previous best
    ## Run 457 stress 0.07365795 
    ## ... Procrustes: rmse 0.0001795808  max resid 0.0004173671 
    ## ... Similar to previous best
    ## Run 458 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001035233  max resid 0.0002465027 
    ## ... Similar to previous best
    ## Run 459 stress 0.08000467 
    ## Run 460 stress 0.08000469 
    ## Run 461 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001332607  max resid 0.0003159677 
    ## ... Similar to previous best
    ## Run 462 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744311  max resid 0.0510331 
    ## Run 463 stress 0.07365783 
    ## ... Procrustes: rmse 7.917613e-06  max resid 1.763518e-05 
    ## ... Similar to previous best
    ## Run 464 stress 0.0800047 
    ## Run 465 stress 0.08288046 
    ## Run 466 stress 0.07629235 
    ## Run 467 stress 0.0828805 
    ## Run 468 stress 0.08000471 
    ## Run 469 stress 0.07365784 
    ## ... Procrustes: rmse 5.457107e-05  max resid 0.0001280543 
    ## ... Similar to previous best
    ## Run 470 stress 0.08000467 
    ## Run 471 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744733  max resid 0.0510526 
    ## Run 472 stress 0.08000474 
    ## Run 473 stress 0.07629237 
    ## Run 474 stress 0.08000469 
    ## Run 475 stress 0.07365784 
    ## ... Procrustes: rmse 5.245479e-05  max resid 0.0001233545 
    ## ... Similar to previous best
    ## Run 476 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745007  max resid 0.05106575 
    ## Run 477 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001095092  max resid 0.000253535 
    ## ... Similar to previous best
    ## Run 478 stress 0.08000468 
    ## Run 479 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744386  max resid 0.05103648 
    ## Run 480 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744619  max resid 0.05104727 
    ## Run 481 stress 0.07365785 
    ## ... Procrustes: rmse 6.851163e-05  max resid 0.0001618656 
    ## ... Similar to previous best
    ## Run 482 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744612  max resid 0.05105911 
    ## Run 483 stress 0.08000467 
    ## Run 484 stress 0.08000469 
    ## Run 485 stress 0.07365783 
    ## ... Procrustes: rmse 2.021781e-05  max resid 4.434953e-05 
    ## ... Similar to previous best
    ## Run 486 stress 0.08233766 
    ## Run 487 stress 0.3169451 
    ## Run 488 stress 0.07629247 
    ## Run 489 stress 0.07732925 
    ## Run 490 stress 0.08000468 
    ## Run 491 stress 0.07629234 
    ## Run 492 stress 0.08288053 
    ## Run 493 stress 0.07365784 
    ## ... Procrustes: rmse 3.590998e-05  max resid 8.603076e-05 
    ## ... Similar to previous best
    ## Run 494 stress 0.07629236 
    ## Run 495 stress 0.07378227 
    ## ... Procrustes: rmse 0.01746211  max resid 0.05109991 
    ## Run 496 stress 0.07629247 
    ## Run 497 stress 0.07365783 
    ## ... Procrustes: rmse 1.503822e-05  max resid 3.540938e-05 
    ## ... Similar to previous best
    ## Run 498 stress 0.07629235 
    ## Run 499 stress 0.07629234 
    ## Run 500 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744991  max resid 0.05106766 
    ## *** Best solution repeated 28 times

``` r
# Stratified lakes and ocean sites
SD_beta_SO_NMDS <- metaMDS(SD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.07970526 
    ## ... New best solution
    ## ... Procrustes: rmse 1.590264e-05  max resid 2.865443e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.08340293 
    ## Run 3 stress 0.07428314 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1089855  max resid 0.2637097 
    ## Run 4 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04475015  max resid 0.1235811 
    ## Run 5 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 5.214059e-05  max resid 0.0001348546 
    ## ... Similar to previous best
    ## Run 6 stress 0.07844914 
    ## Run 7 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323751  max resid 0.0332957 
    ## Run 8 stress 0.07970525 
    ## Run 9 stress 0.07250812 
    ## Run 10 stress 0.07844947 
    ## Run 11 stress 0.08340291 
    ## Run 12 stress 0.08448436 
    ## Run 13 stress 0.06942778 
    ## ... Procrustes: rmse 3.544909e-05  max resid 9.220694e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.07428319 
    ## Run 15 stress 0.08340286 
    ## Run 16 stress 0.08340286 
    ## Run 17 stress 0.07844944 
    ## Run 18 stress 0.07250814 
    ## Run 19 stress 0.0834029 
    ## Run 20 stress 0.06942776 
    ## ... Procrustes: rmse 5.333978e-05  max resid 0.0001377119 
    ## ... Similar to previous best
    ## Run 21 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 3.240393e-05  max resid 7.695054e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.07250812 
    ## Run 23 stress 0.08340292 
    ## Run 24 stress 0.07428314 
    ## Run 25 stress 0.06978195 
    ## ... Procrustes: rmse 0.0131743  max resid 0.03315126 
    ## Run 26 stress 0.07250812 
    ## Run 27 stress 0.08340287 
    ## Run 28 stress 0.08340289 
    ## Run 29 stress 0.07844933 
    ## Run 30 stress 0.07250812 
    ## Run 31 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 1.294976e-05  max resid 2.476207e-05 
    ## ... Similar to previous best
    ## Run 32 stress 0.0834029 
    ## Run 33 stress 0.07250812 
    ## Run 34 stress 0.07970525 
    ## Run 35 stress 0.06942777 
    ## ... Procrustes: rmse 4.35343e-05  max resid 0.0001121614 
    ## ... Similar to previous best
    ## Run 36 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326217  max resid 0.03334263 
    ## Run 37 stress 0.07428322 
    ## Run 38 stress 0.08448443 
    ## Run 39 stress 0.06942776 
    ## ... Procrustes: rmse 2.415254e-05  max resid 6.266023e-05 
    ## ... Similar to previous best
    ## Run 40 stress 0.06942776 
    ## ... Procrustes: rmse 2.756138e-05  max resid 7.101818e-05 
    ## ... Similar to previous best
    ## Run 41 stress 0.06942776 
    ## ... Procrustes: rmse 1.562015e-05  max resid 4.037807e-05 
    ## ... Similar to previous best
    ## Run 42 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324278  max resid 0.03330297 
    ## Run 43 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327588  max resid 0.03336946 
    ## Run 44 stress 0.08340289 
    ## Run 45 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323892  max resid 0.03329449 
    ## Run 46 stress 0.06942776 
    ## ... Procrustes: rmse 8.096473e-06  max resid 2.097533e-05 
    ## ... Similar to previous best
    ## Run 47 stress 0.07428313 
    ## Run 48 stress 0.0742832 
    ## Run 49 stress 0.06942776 
    ## ... Procrustes: rmse 2.159209e-05  max resid 5.565921e-05 
    ## ... Similar to previous best
    ## Run 50 stress 0.06978198 
    ## ... Procrustes: rmse 0.01326764  max resid 0.03335649 
    ## Run 51 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322927  max resid 0.03327433 
    ## Run 52 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326764  max resid 0.0333575 
    ## Run 53 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.226796e-06  max resid 4.998292e-06 
    ## ... Similar to previous best
    ## Run 54 stress 0.07250812 
    ## Run 55 stress 0.08448455 
    ## Run 56 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317185  max resid 0.03315179 
    ## Run 57 stress 0.06942776 
    ## ... Procrustes: rmse 4.443552e-06  max resid 1.066425e-05 
    ## ... Similar to previous best
    ## Run 58 stress 0.07428316 
    ## Run 59 stress 0.07428314 
    ## Run 60 stress 0.0834029 
    ## Run 61 stress 0.07970525 
    ## Run 62 stress 0.07250813 
    ## Run 63 stress 0.08340289 
    ## Run 64 stress 0.08340289 
    ## Run 65 stress 0.07970525 
    ## Run 66 stress 0.06942778 
    ## ... Procrustes: rmse 6.004431e-05  max resid 0.0001557261 
    ## ... Similar to previous best
    ## Run 67 stress 0.07970526 
    ## Run 68 stress 0.07428318 
    ## Run 69 stress 0.06942776 
    ## ... Procrustes: rmse 2.710328e-05  max resid 7.070704e-05 
    ## ... Similar to previous best
    ## Run 70 stress 0.07428321 
    ## Run 71 stress 0.07250812 
    ## Run 72 stress 0.07428315 
    ## Run 73 stress 0.08340289 
    ## Run 74 stress 0.06942776 
    ## ... Procrustes: rmse 5.416819e-06  max resid 1.495535e-05 
    ## ... Similar to previous best
    ## Run 75 stress 0.07428316 
    ## Run 76 stress 0.07428323 
    ## Run 77 stress 0.07428315 
    ## Run 78 stress 0.07250812 
    ## Run 79 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325678  max resid 0.03333181 
    ## Run 80 stress 0.08340287 
    ## Run 81 stress 0.07844913 
    ## Run 82 stress 0.07970527 
    ## Run 83 stress 0.07970525 
    ## Run 84 stress 0.07970525 
    ## Run 85 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324539  max resid 0.03331017 
    ## Run 86 stress 0.0742832 
    ## Run 87 stress 0.07428316 
    ## Run 88 stress 0.08340298 
    ## Run 89 stress 0.06942776 
    ## ... Procrustes: rmse 2.657588e-05  max resid 6.915566e-05 
    ## ... Similar to previous best
    ## Run 90 stress 0.07970525 
    ## Run 91 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320791  max resid 0.03323043 
    ## Run 92 stress 0.06942776 
    ## ... Procrustes: rmse 5.299861e-06  max resid 1.020649e-05 
    ## ... Similar to previous best
    ## Run 93 stress 0.0844844 
    ## Run 94 stress 0.07250812 
    ## Run 95 stress 0.07428314 
    ## Run 96 stress 0.07970526 
    ## Run 97 stress 0.06942776 
    ## ... Procrustes: rmse 6.805422e-06  max resid 1.801975e-05 
    ## ... Similar to previous best
    ## Run 98 stress 0.07428314 
    ## Run 99 stress 0.08448442 
    ## Run 100 stress 0.07428314 
    ## Run 101 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324953  max resid 0.03331751 
    ## Run 102 stress 0.07428317 
    ## Run 103 stress 0.06942776 
    ## ... Procrustes: rmse 2.589855e-05  max resid 6.882937e-05 
    ## ... Similar to previous best
    ## Run 104 stress 0.08340287 
    ## Run 105 stress 0.07250812 
    ## Run 106 stress 0.07428313 
    ## Run 107 stress 0.06942777 
    ## ... Procrustes: rmse 5.223555e-05  max resid 0.0001350156 
    ## ... Similar to previous best
    ## Run 108 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319317  max resid 0.03319783 
    ## Run 109 stress 0.07970527 
    ## Run 110 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324457  max resid 0.03330711 
    ## Run 111 stress 0.07428319 
    ## Run 112 stress 0.06942776 
    ## ... Procrustes: rmse 3.871654e-06  max resid 1.086147e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.07428313 
    ## Run 114 stress 0.07250813 
    ## Run 115 stress 0.0742832 
    ## Run 116 stress 0.07970525 
    ## Run 117 stress 0.07428312 
    ## Run 118 stress 0.08340287 
    ## Run 119 stress 0.07250812 
    ## Run 120 stress 0.08340286 
    ## Run 121 stress 0.07970525 
    ## Run 122 stress 0.06942776 
    ## ... Procrustes: rmse 3.866522e-06  max resid 7.609217e-06 
    ## ... Similar to previous best
    ## Run 123 stress 0.07428313 
    ## Run 124 stress 0.07970526 
    ## Run 125 stress 0.07428319 
    ## Run 126 stress 0.0742832 
    ## Run 127 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327614  max resid 0.03337299 
    ## Run 128 stress 0.07428315 
    ## Run 129 stress 0.06942776 
    ## ... Procrustes: rmse 4.460136e-06  max resid 1.223867e-05 
    ## ... Similar to previous best
    ## Run 130 stress 0.08340292 
    ## Run 131 stress 0.08340289 
    ## Run 132 stress 0.07250815 
    ## Run 133 stress 0.06942776 
    ## ... Procrustes: rmse 2.897048e-05  max resid 7.514605e-05 
    ## ... Similar to previous best
    ## Run 134 stress 0.08340288 
    ## Run 135 stress 0.07844912 
    ## Run 136 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323858  max resid 0.03329311 
    ## Run 137 stress 0.06942777 
    ## ... Procrustes: rmse 3.289809e-05  max resid 8.670702e-05 
    ## ... Similar to previous best
    ## Run 138 stress 0.07428313 
    ## Run 139 stress 0.07970525 
    ## Run 140 stress 0.06942778 
    ## ... Procrustes: rmse 6.897669e-05  max resid 0.0001779039 
    ## ... Similar to previous best
    ## Run 141 stress 0.07970525 
    ## Run 142 stress 0.07428312 
    ## Run 143 stress 0.08340287 
    ## Run 144 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323263  max resid 0.03328219 
    ## Run 145 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323932  max resid 0.03329465 
    ## Run 146 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320496  max resid 0.03322623 
    ## Run 147 stress 0.08340297 
    ## Run 148 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.560718e-06  max resid 5.718922e-06 
    ## ... Similar to previous best
    ## Run 149 stress 0.07428314 
    ## Run 150 stress 0.06942777 
    ## ... Procrustes: rmse 1.323242e-05  max resid 4.178059e-05 
    ## ... Similar to previous best
    ## Run 151 stress 0.07428313 
    ## Run 152 stress 0.08340293 
    ## Run 153 stress 0.07428317 
    ## Run 154 stress 0.07970525 
    ## Run 155 stress 0.07428313 
    ## Run 156 stress 0.07250812 
    ## Run 157 stress 0.07250812 
    ## Run 158 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320748  max resid 0.0332291 
    ## Run 159 stress 0.07428316 
    ## Run 160 stress 0.07428313 
    ## Run 161 stress 0.07428313 
    ## Run 162 stress 0.069782 
    ## ... Procrustes: rmse 0.01327928  max resid 0.03338407 
    ## Run 163 stress 0.06942776 
    ## ... Procrustes: rmse 7.696701e-06  max resid 1.939757e-05 
    ## ... Similar to previous best
    ## Run 164 stress 0.07970525 
    ## Run 165 stress 0.07428313 
    ## Run 166 stress 0.07970525 
    ## Run 167 stress 0.07250812 
    ## Run 168 stress 0.07250814 
    ## Run 169 stress 0.08448443 
    ## Run 170 stress 0.07428312 
    ## Run 171 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326334  max resid 0.0333505 
    ## Run 172 stress 0.07428315 
    ## Run 173 stress 0.0844844 
    ## Run 174 stress 0.07970527 
    ## Run 175 stress 0.07250812 
    ## Run 176 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327374  max resid 0.03337123 
    ## Run 177 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324275  max resid 0.03330444 
    ## Run 178 stress 0.06978196 
    ## ... Procrustes: rmse 0.01327055  max resid 0.03336449 
    ## Run 179 stress 0.07844948 
    ## Run 180 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323487  max resid 0.0332896 
    ## Run 181 stress 0.0784495 
    ## Run 182 stress 0.07428313 
    ## Run 183 stress 0.07970525 
    ## Run 184 stress 0.07428313 
    ## Run 185 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323659  max resid 0.03329283 
    ## Run 186 stress 0.07428321 
    ## Run 187 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323644  max resid 0.03329454 
    ## Run 188 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324361  max resid 0.0333055 
    ## Run 189 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322791  max resid 0.03327555 
    ## Run 190 stress 0.07250813 
    ## Run 191 stress 0.07970525 
    ## Run 192 stress 0.07428313 
    ## Run 193 stress 0.06978197 
    ## ... Procrustes: rmse 0.0132732  max resid 0.03336858 
    ## Run 194 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322963  max resid 0.0332817 
    ## Run 195 stress 0.07428314 
    ## Run 196 stress 0.07250812 
    ## Run 197 stress 0.08340288 
    ## Run 198 stress 0.07428313 
    ## Run 199 stress 0.06942776 
    ## ... Procrustes: rmse 9.975425e-06  max resid 2.986367e-05 
    ## ... Similar to previous best
    ## Run 200 stress 0.07970526 
    ## Run 201 stress 0.07428315 
    ## Run 202 stress 0.07428313 
    ## Run 203 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325314  max resid 0.03332816 
    ## Run 204 stress 0.07970525 
    ## Run 205 stress 0.07428317 
    ## Run 206 stress 0.07428313 
    ## Run 207 stress 0.07970525 
    ## Run 208 stress 0.07428319 
    ## Run 209 stress 0.07844937 
    ## Run 210 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325766  max resid 0.03333763 
    ## Run 211 stress 0.08340296 
    ## Run 212 stress 0.07844966 
    ## Run 213 stress 0.07250812 
    ## Run 214 stress 0.07970525 
    ## Run 215 stress 0.06942776 
    ## ... Procrustes: rmse 1.154338e-05  max resid 2.906104e-05 
    ## ... Similar to previous best
    ## Run 216 stress 0.07844937 
    ## Run 217 stress 0.07428315 
    ## Run 218 stress 0.07428313 
    ## Run 219 stress 0.07428313 
    ## Run 220 stress 0.07428313 
    ## Run 221 stress 0.07428313 
    ## Run 222 stress 0.07970525 
    ## Run 223 stress 0.07970525 
    ## Run 224 stress 0.07250812 
    ## Run 225 stress 0.08448453 
    ## Run 226 stress 0.07970525 
    ## Run 227 stress 0.07428316 
    ## Run 228 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325741  max resid 0.03333739 
    ## Run 229 stress 0.07250814 
    ## Run 230 stress 0.07250812 
    ## Run 231 stress 0.08448441 
    ## Run 232 stress 0.06942776 
    ## ... Procrustes: rmse 5.172112e-06  max resid 1.57092e-05 
    ## ... Similar to previous best
    ## Run 233 stress 0.07428315 
    ## Run 234 stress 0.07428313 
    ## Run 235 stress 0.07250812 
    ## Run 236 stress 0.07428314 
    ## Run 237 stress 0.07428316 
    ## Run 238 stress 0.07844913 
    ## Run 239 stress 0.07250814 
    ## Run 240 stress 0.07428313 
    ## Run 241 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320816  max resid 0.03323363 
    ## Run 242 stress 0.07970526 
    ## Run 243 stress 0.07970526 
    ## Run 244 stress 0.08340287 
    ## Run 245 stress 0.07428313 
    ## Run 246 stress 0.06942776 
    ## ... Procrustes: rmse 2.207528e-05  max resid 5.636722e-05 
    ## ... Similar to previous best
    ## Run 247 stress 0.07428315 
    ## Run 248 stress 0.08340287 
    ## Run 249 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325671  max resid 0.03333518 
    ## Run 250 stress 0.07428312 
    ## Run 251 stress 0.0784494 
    ## Run 252 stress 0.08340291 
    ## Run 253 stress 0.06942776 
    ## ... Procrustes: rmse 1.437859e-05  max resid 3.653478e-05 
    ## ... Similar to previous best
    ## Run 254 stress 0.07428319 
    ## Run 255 stress 0.07970525 
    ## Run 256 stress 0.06942777 
    ## ... Procrustes: rmse 4.209342e-05  max resid 0.0001078518 
    ## ... Similar to previous best
    ## Run 257 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318949  max resid 0.03319312 
    ## Run 258 stress 0.07428313 
    ## Run 259 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323948  max resid 0.03330138 
    ## Run 260 stress 0.07428314 
    ## Run 261 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132351  max resid 0.03328968 
    ## Run 262 stress 0.07250814 
    ## Run 263 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324397  max resid 0.03330949 
    ## Run 264 stress 0.07970525 
    ## Run 265 stress 0.08340287 
    ## Run 266 stress 0.07428318 
    ## Run 267 stress 0.08340289 
    ## Run 268 stress 0.07428313 
    ## Run 269 stress 0.07428316 
    ## Run 270 stress 0.07428313 
    ## Run 271 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321497  max resid 0.03324775 
    ## Run 272 stress 0.069782 
    ## ... Procrustes: rmse 0.01327071  max resid 0.03336926 
    ## Run 273 stress 0.07428313 
    ## Run 274 stress 0.07428317 
    ## Run 275 stress 0.07970526 
    ## Run 276 stress 0.07970525 
    ## Run 277 stress 0.07970526 
    ## Run 278 stress 0.06942776 
    ## ... Procrustes: rmse 1.723375e-05  max resid 4.449854e-05 
    ## ... Similar to previous best
    ## Run 279 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323247  max resid 0.03328633 
    ## Run 280 stress 0.07970525 
    ## Run 281 stress 0.06942776 
    ## ... Procrustes: rmse 1.724268e-05  max resid 4.469122e-05 
    ## ... Similar to previous best
    ## Run 282 stress 0.07250815 
    ## Run 283 stress 0.07250813 
    ## Run 284 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321674  max resid 0.03325292 
    ## Run 285 stress 0.07428315 
    ## Run 286 stress 0.0784492 
    ## Run 287 stress 0.08340288 
    ## Run 288 stress 0.07428313 
    ## Run 289 stress 0.07428319 
    ## Run 290 stress 0.08340287 
    ## Run 291 stress 0.06942776 
    ## ... Procrustes: rmse 4.863711e-06  max resid 1.063988e-05 
    ## ... Similar to previous best
    ## Run 292 stress 0.08340287 
    ## Run 293 stress 0.07250815 
    ## Run 294 stress 0.07970526 
    ## Run 295 stress 0.07428316 
    ## Run 296 stress 0.08340287 
    ## Run 297 stress 0.07970525 
    ## Run 298 stress 0.06942776 
    ## ... Procrustes: rmse 1.034783e-05  max resid 2.757289e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.07428313 
    ## Run 300 stress 0.07428314 
    ## Run 301 stress 0.06942776 
    ## ... Procrustes: rmse 2.345508e-05  max resid 4.7156e-05 
    ## ... Similar to previous best
    ## Run 302 stress 0.07428314 
    ## Run 303 stress 0.06942776 
    ## ... Procrustes: rmse 8.847386e-06  max resid 2.696032e-05 
    ## ... Similar to previous best
    ## Run 304 stress 0.07428313 
    ## Run 305 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326531  max resid 0.03335161 
    ## Run 306 stress 0.08340288 
    ## Run 307 stress 0.08340292 
    ## Run 308 stress 0.08448438 
    ## Run 309 stress 0.08448451 
    ## Run 310 stress 0.08340293 
    ## Run 311 stress 0.07970526 
    ## Run 312 stress 0.07428313 
    ## Run 313 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320116  max resid 0.03321682 
    ## Run 314 stress 0.06942777 
    ## ... Procrustes: rmse 3.333401e-05  max resid 8.583215e-05 
    ## ... Similar to previous best
    ## Run 315 stress 0.07428318 
    ## Run 316 stress 0.06942776 
    ## ... Procrustes: rmse 7.266268e-06  max resid 1.976254e-05 
    ## ... Similar to previous best
    ## Run 317 stress 0.07428314 
    ## Run 318 stress 0.08340287 
    ## Run 319 stress 0.07250812 
    ## Run 320 stress 0.07844968 
    ## Run 321 stress 0.07250813 
    ## Run 322 stress 0.07250812 
    ## Run 323 stress 0.07428315 
    ## Run 324 stress 0.07250812 
    ## Run 325 stress 0.06942777 
    ## ... Procrustes: rmse 4.378414e-05  max resid 0.0001127727 
    ## ... Similar to previous best
    ## Run 326 stress 0.08448445 
    ## Run 327 stress 0.07970525 
    ## Run 328 stress 0.07250812 
    ## Run 329 stress 0.0834029 
    ## Run 330 stress 0.08448452 
    ## Run 331 stress 0.06942776 
    ## ... Procrustes: rmse 5.228884e-06  max resid 1.224599e-05 
    ## ... Similar to previous best
    ## Run 332 stress 0.08340287 
    ## Run 333 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324616  max resid 0.03331528 
    ## Run 334 stress 0.07428318 
    ## Run 335 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320915  max resid 0.03323449 
    ## Run 336 stress 0.07970525 
    ## Run 337 stress 0.07250813 
    ## Run 338 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325913  max resid 0.03334159 
    ## Run 339 stress 0.07250812 
    ## Run 340 stress 0.07844951 
    ## Run 341 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325481  max resid 0.03333561 
    ## Run 342 stress 0.07428313 
    ## Run 343 stress 0.06942776 
    ## ... Procrustes: rmse 1.412447e-05  max resid 3.276454e-05 
    ## ... Similar to previous best
    ## Run 344 stress 0.07428313 
    ## Run 345 stress 0.06942777 
    ## ... Procrustes: rmse 5.408066e-05  max resid 0.0001385695 
    ## ... Similar to previous best
    ## Run 346 stress 0.07250812 
    ## Run 347 stress 0.07844931 
    ## Run 348 stress 0.08340287 
    ## Run 349 stress 0.08340287 
    ## Run 350 stress 0.07250812 
    ## Run 351 stress 0.06942776 
    ## ... Procrustes: rmse 1.676375e-05  max resid 4.271635e-05 
    ## ... Similar to previous best
    ## Run 352 stress 0.08340286 
    ## Run 353 stress 0.06942776 
    ## ... Procrustes: rmse 6.234046e-06  max resid 1.62562e-05 
    ## ... Similar to previous best
    ## Run 354 stress 0.07428315 
    ## Run 355 stress 0.08340287 
    ## Run 356 stress 0.07250816 
    ## Run 357 stress 0.07428317 
    ## Run 358 stress 0.08340295 
    ## Run 359 stress 0.07428312 
    ## Run 360 stress 0.07250813 
    ## Run 361 stress 0.07250812 
    ## Run 362 stress 0.06942776 
    ## ... Procrustes: rmse 1.201157e-05  max resid 3.168059e-05 
    ## ... Similar to previous best
    ## Run 363 stress 0.08340287 
    ## Run 364 stress 0.08340295 
    ## Run 365 stress 0.08340296 
    ## Run 366 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326369  max resid 0.03335027 
    ## Run 367 stress 0.06942777 
    ## ... Procrustes: rmse 4.480591e-05  max resid 0.0001150312 
    ## ... Similar to previous best
    ## Run 368 stress 0.08340289 
    ## Run 369 stress 0.08448435 
    ## Run 370 stress 0.06978191 
    ## ... Procrustes: rmse 0.01321123  max resid 0.03323802 
    ## Run 371 stress 0.08340291 
    ## Run 372 stress 0.07970525 
    ## Run 373 stress 0.07428313 
    ## Run 374 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323217  max resid 0.03328438 
    ## Run 375 stress 0.08340286 
    ## Run 376 stress 0.06942777 
    ## ... Procrustes: rmse 3.87069e-05  max resid 0.000100408 
    ## ... Similar to previous best
    ## Run 377 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323518  max resid 0.03329209 
    ## Run 378 stress 0.07970525 
    ## Run 379 stress 0.07428316 
    ## Run 380 stress 0.07428316 
    ## Run 381 stress 0.07250812 
    ## Run 382 stress 0.08340286 
    ## Run 383 stress 0.06978193 
    ## ... Procrustes: rmse 0.01319107  max resid 0.0331947 
    ## Run 384 stress 0.08448442 
    ## Run 385 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323958  max resid 0.03329929 
    ## Run 386 stress 0.08340287 
    ## Run 387 stress 0.06942776 
    ## ... Procrustes: rmse 1.842408e-05  max resid 4.808212e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.07428313 
    ## Run 389 stress 0.07250812 
    ## Run 390 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319958  max resid 0.03321501 
    ## Run 391 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326291  max resid 0.03335057 
    ## Run 392 stress 0.07250812 
    ## Run 393 stress 0.06942777 
    ## ... Procrustes: rmse 4.546386e-05  max resid 0.0001167582 
    ## ... Similar to previous best
    ## Run 394 stress 0.07250812 
    ## Run 395 stress 0.07428313 
    ## Run 396 stress 0.07250812 
    ## Run 397 stress 0.07428318 
    ## Run 398 stress 0.08448435 
    ## Run 399 stress 0.0844844 
    ## Run 400 stress 0.07970525 
    ## Run 401 stress 0.07428314 
    ## Run 402 stress 0.06942776 
    ## ... Procrustes: rmse 2.813928e-06  max resid 7.287753e-06 
    ## ... Similar to previous best
    ## Run 403 stress 0.07250813 
    ## Run 404 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323849  max resid 0.03329756 
    ## Run 405 stress 0.08340293 
    ## Run 406 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324696  max resid 0.03331567 
    ## Run 407 stress 0.07970525 
    ## Run 408 stress 0.07428313 
    ## Run 409 stress 0.08340287 
    ## Run 410 stress 0.07970525 
    ## Run 411 stress 0.07250812 
    ## Run 412 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322596  max resid 0.0332702 
    ## Run 413 stress 0.07428313 
    ## Run 414 stress 0.06942776 
    ## ... Procrustes: rmse 5.954697e-06  max resid 1.447681e-05 
    ## ... Similar to previous best
    ## Run 415 stress 0.07970526 
    ## Run 416 stress 0.07428314 
    ## Run 417 stress 0.07970525 
    ## Run 418 stress 0.07250814 
    ## Run 419 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132332  max resid 0.03328759 
    ## Run 420 stress 0.07428313 
    ## Run 421 stress 0.06942776 
    ## ... Procrustes: rmse 5.433515e-06  max resid 1.203202e-05 
    ## ... Similar to previous best
    ## Run 422 stress 0.08340292 
    ## Run 423 stress 0.07428321 
    ## Run 424 stress 0.06942776 
    ## ... Procrustes: rmse 6.498993e-06  max resid 1.487954e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.07970525 
    ## Run 426 stress 0.08340287 
    ## Run 427 stress 0.08340286 
    ## Run 428 stress 0.06942777 
    ## ... Procrustes: rmse 3.672446e-05  max resid 9.53087e-05 
    ## ... Similar to previous best
    ## Run 429 stress 0.0742832 
    ## Run 430 stress 0.08448439 
    ## Run 431 stress 0.07250813 
    ## Run 432 stress 0.07428318 
    ## Run 433 stress 0.07428318 
    ## Run 434 stress 0.07970526 
    ## Run 435 stress 0.07428313 
    ## Run 436 stress 0.07970526 
    ## Run 437 stress 0.07250812 
    ## Run 438 stress 0.06942776 
    ## ... Procrustes: rmse 7.329557e-06  max resid 1.903862e-05 
    ## ... Similar to previous best
    ## Run 439 stress 0.07844954 
    ## Run 440 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325612  max resid 0.03333305 
    ## Run 441 stress 0.08340287 
    ## Run 442 stress 0.08448441 
    ## Run 443 stress 0.07970526 
    ## Run 444 stress 0.08448443 
    ## Run 445 stress 0.07844915 
    ## Run 446 stress 0.07970525 
    ## Run 447 stress 0.07970526 
    ## Run 448 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323449  max resid 0.03328836 
    ## Run 449 stress 0.06942777 
    ## ... Procrustes: rmse 3.683298e-05  max resid 9.530413e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.07250812 
    ## Run 451 stress 0.07250813 
    ## Run 452 stress 0.08340298 
    ## Run 453 stress 0.07970525 
    ## Run 454 stress 0.07970527 
    ## Run 455 stress 0.07844946 
    ## Run 456 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323034  max resid 0.03327797 
    ## Run 457 stress 0.07970525 
    ## Run 458 stress 0.07250813 
    ## Run 459 stress 0.07844945 
    ## Run 460 stress 0.07250812 
    ## Run 461 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326068  max resid 0.03334504 
    ## Run 462 stress 0.07428314 
    ## Run 463 stress 0.07428317 
    ## Run 464 stress 0.07844928 
    ## Run 465 stress 0.07250813 
    ## Run 466 stress 0.07428314 
    ## Run 467 stress 0.08340288 
    ## Run 468 stress 0.0742832 
    ## Run 469 stress 0.07428314 
    ## Run 470 stress 0.07428317 
    ## Run 471 stress 0.08340286 
    ## Run 472 stress 0.07428317 
    ## Run 473 stress 0.07428314 
    ## Run 474 stress 0.08448457 
    ## Run 475 stress 0.08340287 
    ## Run 476 stress 0.08448437 
    ## Run 477 stress 0.08340286 
    ## Run 478 stress 0.07970525 
    ## Run 479 stress 0.07428316 
    ## Run 480 stress 0.08340287 
    ## Run 481 stress 0.07428312 
    ## Run 482 stress 0.07970525 
    ## Run 483 stress 0.08340291 
    ## Run 484 stress 0.08340293 
    ## Run 485 stress 0.07250812 
    ## Run 486 stress 0.08340296 
    ## Run 487 stress 0.07428313 
    ## Run 488 stress 0.06942776 
    ## ... Procrustes: rmse 8.768351e-06  max resid 2.742415e-05 
    ## ... Similar to previous best
    ## Run 489 stress 0.07970526 
    ## Run 490 stress 0.0742832 
    ## Run 491 stress 0.08448444 
    ## Run 492 stress 0.07428316 
    ## Run 493 stress 0.07428313 
    ## Run 494 stress 0.07428322 
    ## Run 495 stress 0.0742832 
    ## Run 496 stress 0.07428313 
    ## Run 497 stress 0.07428313 
    ## Run 498 stress 0.06942776 
    ## ... Procrustes: rmse 6.398391e-06  max resid 1.582371e-05 
    ## ... Similar to previous best
    ## Run 499 stress 0.07970527 
    ## Run 500 stress 0.07970525 
    ## *** Best solution repeated 37 times

``` r
# Stratified lakes
SD_beta_S_NMDS <- metaMDS(SD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1410449 
    ## Run 1 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2673132  max resid 0.5425958 
    ## Run 2 stress 0.1407206 
    ## Run 3 stress 0.1410453 
    ## Run 4 stress 0.2550044 
    ## Run 5 stress 0.1434197 
    ## Run 6 stress 0.1693717 
    ## Run 7 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 3.850228e-07  max resid 8.447166e-07 
    ## ... Similar to previous best
    ## Run 8 stress 0.1407298 
    ## Run 9 stress 0.1410451 
    ## Run 10 stress 0.1320258 
    ## ... Procrustes: rmse 2.498045e-06  max resid 4.928453e-06 
    ## ... Similar to previous best
    ## Run 11 stress 0.1410448 
    ## Run 12 stress 0.1410449 
    ## Run 13 stress 0.1551863 
    ## Run 14 stress 0.1583016 
    ## Run 15 stress 0.1410452 
    ## Run 16 stress 0.1572668 
    ## Run 17 stress 0.1407207 
    ## Run 18 stress 0.1383682 
    ## Run 19 stress 0.1572668 
    ## Run 20 stress 0.1410448 
    ## Run 21 stress 0.1802587 
    ## Run 22 stress 0.1587623 
    ## Run 23 stress 0.1383682 
    ## Run 24 stress 0.1410452 
    ## Run 25 stress 0.1320258 
    ## ... Procrustes: rmse 2.799069e-06  max resid 5.824704e-06 
    ## ... Similar to previous best
    ## Run 26 stress 0.1551863 
    ## Run 27 stress 0.1410453 
    ## Run 28 stress 0.1410453 
    ## Run 29 stress 0.1415299 
    ## Run 30 stress 0.2008158 
    ## Run 31 stress 0.1410449 
    ## Run 32 stress 0.1383682 
    ## Run 33 stress 0.1962476 
    ## Run 34 stress 0.141045 
    ## Run 35 stress 0.1320258 
    ## ... Procrustes: rmse 2.073987e-06  max resid 4.390476e-06 
    ## ... Similar to previous best
    ## Run 36 stress 0.2227526 
    ## Run 37 stress 0.1572668 
    ## Run 38 stress 0.1415299 
    ## Run 39 stress 0.1410448 
    ## Run 40 stress 0.1434197 
    ## Run 41 stress 0.2074935 
    ## Run 42 stress 0.1320258 
    ## ... Procrustes: rmse 2.514794e-06  max resid 5.30412e-06 
    ## ... Similar to previous best
    ## Run 43 stress 0.1551863 
    ## Run 44 stress 0.1572668 
    ## Run 45 stress 0.1830332 
    ## Run 46 stress 0.1320258 
    ## ... Procrustes: rmse 1.52116e-06  max resid 3.019656e-06 
    ## ... Similar to previous best
    ## Run 47 stress 0.1598154 
    ## Run 48 stress 0.2799115 
    ## Run 49 stress 0.1445811 
    ## Run 50 stress 0.1572668 
    ## Run 51 stress 0.1830332 
    ## Run 52 stress 0.1410452 
    ## Run 53 stress 0.1407207 
    ## Run 54 stress 0.1693717 
    ## Run 55 stress 0.1407298 
    ## Run 56 stress 0.1320258 
    ## ... Procrustes: rmse 9.562175e-07  max resid 1.978439e-06 
    ## ... Similar to previous best
    ## Run 57 stress 0.1320258 
    ## ... Procrustes: rmse 3.002653e-06  max resid 6.249861e-06 
    ## ... Similar to previous best
    ## Run 58 stress 0.1407298 
    ## Run 59 stress 0.1693717 
    ## Run 60 stress 0.1551863 
    ## Run 61 stress 0.1407207 
    ## Run 62 stress 0.1407298 
    ## Run 63 stress 0.1693717 
    ## Run 64 stress 0.1407207 
    ## Run 65 stress 0.2550046 
    ## Run 66 stress 0.1663544 
    ## Run 67 stress 0.1320258 
    ## ... Procrustes: rmse 6.514467e-07  max resid 1.068175e-06 
    ## ... Similar to previous best
    ## Run 68 stress 0.1970303 
    ## Run 69 stress 0.1415299 
    ## Run 70 stress 0.1407208 
    ## Run 71 stress 0.1587623 
    ## Run 72 stress 0.200816 
    ## Run 73 stress 0.1771969 
    ## Run 74 stress 0.1551863 
    ## Run 75 stress 0.1551863 
    ## Run 76 stress 0.1830332 
    ## Run 77 stress 0.1572668 
    ## Run 78 stress 0.1383681 
    ## Run 79 stress 0.1771969 
    ## Run 80 stress 0.1320258 
    ## ... Procrustes: rmse 1.093625e-06  max resid 2.36577e-06 
    ## ... Similar to previous best
    ## Run 81 stress 0.1551863 
    ## Run 82 stress 0.1415299 
    ## Run 83 stress 0.1628606 
    ## Run 84 stress 0.1410452 
    ## Run 85 stress 0.1410451 
    ## Run 86 stress 0.1410449 
    ## Run 87 stress 0.1407208 
    ## Run 88 stress 0.1410449 
    ## Run 89 stress 0.1551863 
    ## Run 90 stress 0.1383681 
    ## Run 91 stress 0.1320258 
    ## ... Procrustes: rmse 2.576603e-06  max resid 5.472237e-06 
    ## ... Similar to previous best
    ## Run 92 stress 0.1445811 
    ## Run 93 stress 0.1320258 
    ## ... Procrustes: rmse 4.499744e-07  max resid 8.67052e-07 
    ## ... Similar to previous best
    ## Run 94 stress 0.1415299 
    ## Run 95 stress 0.1320258 
    ## ... Procrustes: rmse 1.417031e-06  max resid 3.041168e-06 
    ## ... Similar to previous best
    ## Run 96 stress 0.1407206 
    ## Run 97 stress 0.2456775 
    ## Run 98 stress 0.1383682 
    ## Run 99 stress 0.2384387 
    ## Run 100 stress 0.1415299 
    ## Run 101 stress 0.1970303 
    ## Run 102 stress 0.3120129 
    ## Run 103 stress 0.1693717 
    ## Run 104 stress 0.3083098 
    ## Run 105 stress 0.1320258 
    ## ... Procrustes: rmse 9.289054e-07  max resid 2.052894e-06 
    ## ... Similar to previous best
    ## Run 106 stress 0.1415299 
    ## Run 107 stress 0.1415299 
    ## Run 108 stress 0.1663544 
    ## Run 109 stress 0.1415299 
    ## Run 110 stress 0.1445811 
    ## Run 111 stress 0.1410448 
    ## Run 112 stress 0.1407298 
    ## Run 113 stress 0.2714945 
    ## Run 114 stress 0.1383681 
    ## Run 115 stress 0.1771969 
    ## Run 116 stress 0.2085297 
    ## Run 117 stress 0.1771969 
    ## Run 118 stress 0.1407298 
    ## Run 119 stress 0.2520602 
    ## Run 120 stress 0.1929519 
    ## Run 121 stress 0.2385029 
    ## Run 122 stress 0.1551863 
    ## Run 123 stress 0.1407208 
    ## Run 124 stress 0.1320258 
    ## ... Procrustes: rmse 5.771195e-07  max resid 1.131809e-06 
    ## ... Similar to previous best
    ## Run 125 stress 0.1445811 
    ## Run 126 stress 0.1407207 
    ## Run 127 stress 0.1572668 
    ## Run 128 stress 0.1407298 
    ## Run 129 stress 0.1445811 
    ## Run 130 stress 0.1415299 
    ## Run 131 stress 0.1802752 
    ## Run 132 stress 0.2015531 
    ## Run 133 stress 0.2520602 
    ## Run 134 stress 0.1410451 
    ## Run 135 stress 0.1598154 
    ## Run 136 stress 0.2015531 
    ## Run 137 stress 0.1551863 
    ## Run 138 stress 0.1551863 
    ## Run 139 stress 0.1410449 
    ## Run 140 stress 0.1693717 
    ## Run 141 stress 0.1410448 
    ## Run 142 stress 0.1410452 
    ## Run 143 stress 0.1410449 
    ## Run 144 stress 0.2085297 
    ## Run 145 stress 0.2385029 
    ## Run 146 stress 0.1407298 
    ## Run 147 stress 0.1320258 
    ## ... Procrustes: rmse 8.420639e-07  max resid 1.798076e-06 
    ## ... Similar to previous best
    ## Run 148 stress 0.2538018 
    ## Run 149 stress 0.1415299 
    ## Run 150 stress 0.2842805 
    ## Run 151 stress 0.1410448 
    ## Run 152 stress 0.1320258 
    ## ... Procrustes: rmse 1.213393e-06  max resid 2.318277e-06 
    ## ... Similar to previous best
    ## Run 153 stress 0.1693717 
    ## Run 154 stress 0.1693717 
    ## Run 155 stress 0.2224298 
    ## Run 156 stress 0.1407298 
    ## Run 157 stress 0.1410453 
    ## Run 158 stress 0.1598154 
    ## Run 159 stress 0.1410452 
    ## Run 160 stress 0.1628606 
    ## Run 161 stress 0.1551863 
    ## Run 162 stress 0.2015531 
    ## Run 163 stress 0.2823279 
    ## Run 164 stress 0.1572668 
    ## Run 165 stress 0.2085297 
    ## Run 166 stress 0.1410452 
    ## Run 167 stress 0.200816 
    ## Run 168 stress 0.1407298 
    ## Run 169 stress 0.1572668 
    ## Run 170 stress 0.1383682 
    ## Run 171 stress 0.1383681 
    ## Run 172 stress 0.1802752 
    ## Run 173 stress 0.1383681 
    ## Run 174 stress 0.3083098 
    ## Run 175 stress 0.1415299 
    ## Run 176 stress 0.1320258 
    ## ... Procrustes: rmse 1.600202e-06  max resid 3.086979e-06 
    ## ... Similar to previous best
    ## Run 177 stress 0.141045 
    ## Run 178 stress 0.1551863 
    ## Run 179 stress 0.2510079 
    ## Run 180 stress 0.1415299 
    ## Run 181 stress 0.1802587 
    ## Run 182 stress 0.1383682 
    ## Run 183 stress 0.1320258 
    ## ... Procrustes: rmse 3.44507e-06  max resid 6.961282e-06 
    ## ... Similar to previous best
    ## Run 184 stress 0.1383681 
    ## Run 185 stress 0.1320258 
    ## ... Procrustes: rmse 1.77395e-06  max resid 3.765974e-06 
    ## ... Similar to previous best
    ## Run 186 stress 0.1572668 
    ## Run 187 stress 0.1415299 
    ## Run 188 stress 0.1415299 
    ## Run 189 stress 0.2227526 
    ## Run 190 stress 0.1415299 
    ## Run 191 stress 0.1410452 
    ## Run 192 stress 0.2842805 
    ## Run 193 stress 0.1320258 
    ## ... Procrustes: rmse 2.306547e-06  max resid 4.586992e-06 
    ## ... Similar to previous best
    ## Run 194 stress 0.1583016 
    ## Run 195 stress 0.1771969 
    ## Run 196 stress 0.1320258 
    ## ... Procrustes: rmse 1.436919e-06  max resid 2.681773e-06 
    ## ... Similar to previous best
    ## Run 197 stress 0.1551863 
    ## Run 198 stress 0.1415299 
    ## Run 199 stress 0.1320258 
    ## ... Procrustes: rmse 1.542882e-06  max resid 3.005992e-06 
    ## ... Similar to previous best
    ## Run 200 stress 0.2015531 
    ## Run 201 stress 0.1434197 
    ## Run 202 stress 0.1410448 
    ## Run 203 stress 0.1407208 
    ## Run 204 stress 0.1970303 
    ## Run 205 stress 0.1415299 
    ## Run 206 stress 0.1445811 
    ## Run 207 stress 0.1383681 
    ## Run 208 stress 0.1771966 
    ## Run 209 stress 0.1383681 
    ## Run 210 stress 0.1407206 
    ## Run 211 stress 0.2385029 
    ## Run 212 stress 0.1407298 
    ## Run 213 stress 0.1407298 
    ## Run 214 stress 0.1445811 
    ## Run 215 stress 0.1572668 
    ## Run 216 stress 0.1320258 
    ## ... Procrustes: rmse 7.511981e-07  max resid 1.688886e-06 
    ## ... Similar to previous best
    ## Run 217 stress 0.141045 
    ## Run 218 stress 0.1445811 
    ## Run 219 stress 0.1445811 
    ## Run 220 stress 0.2514078 
    ## Run 221 stress 0.3037154 
    ## Run 222 stress 0.1551863 
    ## Run 223 stress 0.1572668 
    ## Run 224 stress 0.1572668 
    ## Run 225 stress 0.141045 
    ## Run 226 stress 0.1407298 
    ## Run 227 stress 0.1693717 
    ## Run 228 stress 0.1551863 
    ## Run 229 stress 0.1383681 
    ## Run 230 stress 0.1383681 
    ## Run 231 stress 0.1551863 
    ## Run 232 stress 0.1407298 
    ## Run 233 stress 0.1830332 
    ## Run 234 stress 0.1572668 
    ## Run 235 stress 0.1929519 
    ## Run 236 stress 0.200816 
    ## Run 237 stress 0.1551863 
    ## Run 238 stress 0.1583016 
    ## Run 239 stress 0.1551863 
    ## Run 240 stress 0.1998955 
    ## Run 241 stress 0.1415299 
    ## Run 242 stress 0.1572668 
    ## Run 243 stress 0.1415299 
    ## Run 244 stress 0.1551863 
    ## Run 245 stress 0.2496387 
    ## Run 246 stress 0.1415299 
    ## Run 247 stress 0.1771969 
    ## Run 248 stress 0.1320258 
    ## ... Procrustes: rmse 1.134313e-06  max resid 2.452017e-06 
    ## ... Similar to previous best
    ## Run 249 stress 0.1929519 
    ## Run 250 stress 0.1572668 
    ## Run 251 stress 0.1445811 
    ## Run 252 stress 0.1572668 
    ## Run 253 stress 0.1320258 
    ## ... Procrustes: rmse 7.683677e-07  max resid 1.224558e-06 
    ## ... Similar to previous best
    ## Run 254 stress 0.1693717 
    ## Run 255 stress 0.1407298 
    ## Run 256 stress 0.1583016 
    ## Run 257 stress 0.2332576 
    ## Run 258 stress 0.1410449 
    ## Run 259 stress 0.1551863 
    ## Run 260 stress 0.1407298 
    ## Run 261 stress 0.1415299 
    ## Run 262 stress 0.1320258 
    ## ... Procrustes: rmse 8.425253e-07  max resid 1.527239e-06 
    ## ... Similar to previous best
    ## Run 263 stress 0.1802752 
    ## Run 264 stress 0.1383681 
    ## Run 265 stress 0.1320258 
    ## ... Procrustes: rmse 8.494969e-07  max resid 1.908295e-06 
    ## ... Similar to previous best
    ## Run 266 stress 0.1581642 
    ## Run 267 stress 0.1320258 
    ## ... Procrustes: rmse 7.044247e-07  max resid 1.088792e-06 
    ## ... Similar to previous best
    ## Run 268 stress 0.1320258 
    ## ... Procrustes: rmse 7.203405e-07  max resid 1.621875e-06 
    ## ... Similar to previous best
    ## Run 269 stress 0.1434197 
    ## Run 270 stress 0.1383681 
    ## Run 271 stress 0.2885554 
    ## Run 272 stress 0.2501065 
    ## Run 273 stress 0.2852164 
    ## Run 274 stress 0.2385029 
    ## Run 275 stress 0.2085297 
    ## Run 276 stress 0.1771969 
    ## Run 277 stress 0.1551863 
    ## Run 278 stress 0.1320258 
    ## ... Procrustes: rmse 1.33279e-06  max resid 2.517146e-06 
    ## ... Similar to previous best
    ## Run 279 stress 0.1410452 
    ## Run 280 stress 0.2008163 
    ## Run 281 stress 0.1410452 
    ## Run 282 stress 0.1830332 
    ## Run 283 stress 0.1320258 
    ## ... Procrustes: rmse 9.706365e-07  max resid 2.067323e-06 
    ## ... Similar to previous best
    ## Run 284 stress 0.1415299 
    ## Run 285 stress 0.1410449 
    ## Run 286 stress 0.1771969 
    ## Run 287 stress 0.1383682 
    ## Run 288 stress 0.1796864 
    ## Run 289 stress 0.2008162 
    ## Run 290 stress 0.1572668 
    ## Run 291 stress 0.1970303 
    ## Run 292 stress 0.1407298 
    ## Run 293 stress 0.1551863 
    ## Run 294 stress 0.1572668 
    ## Run 295 stress 0.2703379 
    ## Run 296 stress 0.2385029 
    ## Run 297 stress 0.200816 
    ## Run 298 stress 0.2385029 
    ## Run 299 stress 0.1407298 
    ## Run 300 stress 0.2578647 
    ## Run 301 stress 0.1320258 
    ## ... Procrustes: rmse 2.029514e-06  max resid 4.269323e-06 
    ## ... Similar to previous best
    ## Run 302 stress 0.1407298 
    ## Run 303 stress 0.1970303 
    ## Run 304 stress 0.1415299 
    ## Run 305 stress 0.1383681 
    ## Run 306 stress 0.2015531 
    ## Run 307 stress 0.1776793 
    ## Run 308 stress 0.1415299 
    ## Run 309 stress 0.1598154 
    ## Run 310 stress 0.1415299 
    ## Run 311 stress 0.1796865 
    ## Run 312 stress 0.1551863 
    ## Run 313 stress 0.1383681 
    ## Run 314 stress 0.1320258 
    ## ... Procrustes: rmse 1.848157e-06  max resid 3.69102e-06 
    ## ... Similar to previous best
    ## Run 315 stress 0.1693717 
    ## Run 316 stress 0.1415299 
    ## Run 317 stress 0.1320258 
    ## ... Procrustes: rmse 3.326582e-07  max resid 7.258581e-07 
    ## ... Similar to previous best
    ## Run 318 stress 0.2074935 
    ## Run 319 stress 0.1407208 
    ## Run 320 stress 0.1320258 
    ## ... Procrustes: rmse 1.928995e-06  max resid 3.811097e-06 
    ## ... Similar to previous best
    ## Run 321 stress 0.1410454 
    ## Run 322 stress 0.1320258 
    ## ... Procrustes: rmse 4.885388e-07  max resid 1.006627e-06 
    ## ... Similar to previous best
    ## Run 323 stress 0.1383682 
    ## Run 324 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 2.047216e-07  max resid 3.053363e-07 
    ## ... Similar to previous best
    ## Run 325 stress 0.1581642 
    ## Run 326 stress 0.1802587 
    ## Run 327 stress 0.1572668 
    ## Run 328 stress 0.2422742 
    ## Run 329 stress 0.1693717 
    ## Run 330 stress 0.2538018 
    ## Run 331 stress 0.1693717 
    ## Run 332 stress 0.1415299 
    ## Run 333 stress 0.2854381 
    ## Run 334 stress 0.2550043 
    ## Run 335 stress 0.1581642 
    ## Run 336 stress 0.1320258 
    ## ... Procrustes: rmse 2.614798e-06  max resid 5.186456e-06 
    ## ... Similar to previous best
    ## Run 337 stress 0.1771969 
    ## Run 338 stress 0.2848518 
    ## Run 339 stress 0.1572668 
    ## Run 340 stress 0.2562346 
    ## Run 341 stress 0.1415299 
    ## Run 342 stress 0.1407298 
    ## Run 343 stress 0.1551863 
    ## Run 344 stress 0.1410448 
    ## Run 345 stress 0.1415299 
    ## Run 346 stress 0.1415299 
    ## Run 347 stress 0.2456775 
    ## Run 348 stress 0.1320258 
    ## ... Procrustes: rmse 7.570135e-07  max resid 1.562995e-06 
    ## ... Similar to previous best
    ## Run 349 stress 0.1410449 
    ## Run 350 stress 0.1320258 
    ## ... Procrustes: rmse 1.925907e-06  max resid 3.865978e-06 
    ## ... Similar to previous best
    ## Run 351 stress 0.1320258 
    ## ... Procrustes: rmse 1.240729e-06  max resid 2.61039e-06 
    ## ... Similar to previous best
    ## Run 352 stress 0.1320258 
    ## ... Procrustes: rmse 2.990074e-06  max resid 6.168111e-06 
    ## ... Similar to previous best
    ## Run 353 stress 0.1407298 
    ## Run 354 stress 0.1598154 
    ## Run 355 stress 0.1929519 
    ## Run 356 stress 0.1929519 
    ## Run 357 stress 0.1410448 
    ## Run 358 stress 0.1771969 
    ## Run 359 stress 0.1410451 
    ## Run 360 stress 0.1383681 
    ## Run 361 stress 0.1693717 
    ## Run 362 stress 0.1693718 
    ## Run 363 stress 0.1383682 
    ## Run 364 stress 0.1407298 
    ## Run 365 stress 0.2578648 
    ## Run 366 stress 0.1410453 
    ## Run 367 stress 0.1551863 
    ## Run 368 stress 0.1320258 
    ## ... Procrustes: rmse 5.506629e-07  max resid 1.14684e-06 
    ## ... Similar to previous best
    ## Run 369 stress 0.1383681 
    ## Run 370 stress 0.1415299 
    ## Run 371 stress 0.1796864 
    ## Run 372 stress 0.1445811 
    ## Run 373 stress 0.1551863 
    ## Run 374 stress 0.2008158 
    ## Run 375 stress 0.1383681 
    ## Run 376 stress 0.1415299 
    ## Run 377 stress 0.1434197 
    ## Run 378 stress 0.1407208 
    ## Run 379 stress 0.1663544 
    ## Run 380 stress 0.1383682 
    ## Run 381 stress 0.1830332 
    ## Run 382 stress 0.1598154 
    ## Run 383 stress 0.1320258 
    ## ... Procrustes: rmse 4.897871e-07  max resid 9.54228e-07 
    ## ... Similar to previous best
    ## Run 384 stress 0.257864 
    ## Run 385 stress 0.3083098 
    ## Run 386 stress 0.1407298 
    ## Run 387 stress 0.1551863 
    ## Run 388 stress 0.1796865 
    ## Run 389 stress 0.1551863 
    ## Run 390 stress 0.1410449 
    ## Run 391 stress 0.1410449 
    ## Run 392 stress 0.1771966 
    ## Run 393 stress 0.1410448 
    ## Run 394 stress 0.1663544 
    ## Run 395 stress 0.1320258 
    ## ... Procrustes: rmse 2.710501e-06  max resid 4.861409e-06 
    ## ... Similar to previous best
    ## Run 396 stress 0.1320258 
    ## ... Procrustes: rmse 5.912072e-07  max resid 8.407643e-07 
    ## ... Similar to previous best
    ## Run 397 stress 0.1320258 
    ## ... Procrustes: rmse 2.029282e-06  max resid 4.203489e-06 
    ## ... Similar to previous best
    ## Run 398 stress 0.1551863 
    ## Run 399 stress 0.1415299 
    ## Run 400 stress 0.1410448 
    ## Run 401 stress 0.1929519 
    ## Run 402 stress 0.2008158 
    ## Run 403 stress 0.1320258 
    ## ... Procrustes: rmse 4.016063e-07  max resid 6.851829e-07 
    ## ... Similar to previous best
    ## Run 404 stress 0.1598154 
    ## Run 405 stress 0.1415299 
    ## Run 406 stress 0.1407206 
    ## Run 407 stress 0.1415299 
    ## Run 408 stress 0.1929519 
    ## Run 409 stress 0.1598154 
    ## Run 410 stress 0.1587623 
    ## Run 411 stress 0.1415299 
    ## Run 412 stress 0.1320258 
    ## ... Procrustes: rmse 1.397247e-06  max resid 2.815913e-06 
    ## ... Similar to previous best
    ## Run 413 stress 0.1551863 
    ## Run 414 stress 0.1583016 
    ## Run 415 stress 0.2562346 
    ## Run 416 stress 0.1383682 
    ## Run 417 stress 0.2227526 
    ## Run 418 stress 0.1383682 
    ## Run 419 stress 0.1962476 
    ## Run 420 stress 0.1320258 
    ## ... Procrustes: rmse 1.027127e-06  max resid 2.051036e-06 
    ## ... Similar to previous best
    ## Run 421 stress 0.1410451 
    ## Run 422 stress 0.1572668 
    ## Run 423 stress 0.2501067 
    ## Run 424 stress 0.1415299 
    ## Run 425 stress 0.1407298 
    ## Run 426 stress 0.1551863 
    ## Run 427 stress 0.3083098 
    ## Run 428 stress 0.1551863 
    ## Run 429 stress 0.1320258 
    ## ... Procrustes: rmse 2.492389e-06  max resid 5.178212e-06 
    ## ... Similar to previous best
    ## Run 430 stress 0.1383682 
    ## Run 431 stress 0.1628606 
    ## Run 432 stress 0.2074935 
    ## Run 433 stress 0.1407298 
    ## Run 434 stress 0.1407209 
    ## Run 435 stress 0.1802587 
    ## Run 436 stress 0.1407206 
    ## Run 437 stress 0.1407207 
    ## Run 438 stress 0.2422742 
    ## Run 439 stress 0.1410452 
    ## Run 440 stress 0.1572668 
    ## Run 441 stress 0.1587623 
    ## Run 442 stress 0.1320258 
    ## ... Procrustes: rmse 6.782075e-07  max resid 1.24709e-06 
    ## ... Similar to previous best
    ## Run 443 stress 0.1383682 
    ## Run 444 stress 0.1572668 
    ## Run 445 stress 0.1410448 
    ## Run 446 stress 0.1407298 
    ## Run 447 stress 0.1572668 
    ## Run 448 stress 0.1572668 
    ## Run 449 stress 0.2711187 
    ## Run 450 stress 0.2501066 
    ## Run 451 stress 0.1320258 
    ## ... Procrustes: rmse 1.677283e-06  max resid 3.449009e-06 
    ## ... Similar to previous best
    ## Run 452 stress 0.1415299 
    ## Run 453 stress 0.1572668 
    ## Run 454 stress 0.2842805 
    ## Run 455 stress 0.1583016 
    ## Run 456 stress 0.1415299 
    ## Run 457 stress 0.1407298 
    ## Run 458 stress 0.1407298 
    ## Run 459 stress 0.1445811 
    ## Run 460 stress 0.1434197 
    ## Run 461 stress 0.1320258 
    ## ... Procrustes: rmse 2.269653e-06  max resid 4.592319e-06 
    ## ... Similar to previous best
    ## Run 462 stress 0.1415299 
    ## Run 463 stress 0.1962476 
    ## Run 464 stress 0.1383682 
    ## Run 465 stress 0.1830332 
    ## Run 466 stress 0.2085297 
    ## Run 467 stress 0.1445811 
    ## Run 468 stress 0.1320258 
    ## ... Procrustes: rmse 1.626097e-06  max resid 3.4127e-06 
    ## ... Similar to previous best
    ## Run 469 stress 0.1407298 
    ## Run 470 stress 0.1693717 
    ## Run 471 stress 0.1415299 
    ## Run 472 stress 0.1383681 
    ## Run 473 stress 0.1415299 
    ## Run 474 stress 0.1551863 
    ## Run 475 stress 0.1320258 
    ## ... Procrustes: rmse 2.885116e-07  max resid 5.384423e-07 
    ## ... Similar to previous best
    ## Run 476 stress 0.1383681 
    ## Run 477 stress 0.1410449 
    ## Run 478 stress 0.1320258 
    ## ... Procrustes: rmse 1.1176e-06  max resid 2.327715e-06 
    ## ... Similar to previous best
    ## Run 479 stress 0.1572668 
    ## Run 480 stress 0.2600627 
    ## Run 481 stress 0.1383682 
    ## Run 482 stress 0.1970303 
    ## Run 483 stress 0.1410448 
    ## Run 484 stress 0.1970303 
    ## Run 485 stress 0.1572668 
    ## Run 486 stress 0.1410454 
    ## Run 487 stress 0.1572668 
    ## Run 488 stress 0.1693717 
    ## Run 489 stress 0.1551863 
    ## Run 490 stress 0.1410448 
    ## Run 491 stress 0.1551863 
    ## Run 492 stress 0.1693717 
    ## Run 493 stress 0.1415299 
    ## Run 494 stress 0.1383681 
    ## Run 495 stress 0.1587623 
    ## Run 496 stress 0.1410448 
    ## Run 497 stress 0.1415299 
    ## Run 498 stress 0.1796865 
    ## Run 499 stress 0.1693717 
    ## Run 500 stress 0.1320258 
    ## ... Procrustes: rmse 1.599334e-06  max resid 3.297736e-06 
    ## ... Similar to previous best
    ## *** Best solution repeated 22 times

``` r
### Environmental
# Surveyed sites 
SD_beta_env_NMDS <- metaMDS(SD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07669178 
    ## Run 1 stress 0.08317359 
    ## Run 2 stress 0.08252135 
    ## Run 3 stress 0.08063924 
    ## Run 4 stress 0.07669157 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003119394  max resid 0.0006650224 
    ## ... Similar to previous best
    ## Run 5 stress 0.07871738 
    ## Run 6 stress 0.08835171 
    ## Run 7 stress 0.08003203 
    ## Run 8 stress 0.08252144 
    ## Run 9 stress 0.08175702 
    ## Run 10 stress 0.07731523 
    ## Run 11 stress 0.0766916 
    ## ... Procrustes: rmse 0.0001622025  max resid 0.0003536783 
    ## ... Similar to previous best
    ## Run 12 stress 0.08012943 
    ## Run 13 stress 0.0766915 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001472288  max resid 0.0002615547 
    ## ... Similar to previous best
    ## Run 14 stress 0.0837571 
    ## Run 15 stress 0.0831733 
    ## Run 16 stress 0.08553054 
    ## Run 17 stress 0.08527022 
    ## Run 18 stress 0.08486852 
    ## Run 19 stress 0.0800321 
    ## Run 20 stress 0.08433703 
    ## Run 21 stress 0.07669151 
    ## ... Procrustes: rmse 3.832017e-05  max resid 6.892472e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.07669157 
    ## ... Procrustes: rmse 0.000196546  max resid 0.0003399814 
    ## ... Similar to previous best
    ## Run 23 stress 0.08218516 
    ## Run 24 stress 0.08252117 
    ## Run 25 stress 0.0887651 
    ## Run 26 stress 0.08012966 
    ## Run 27 stress 0.08577448 
    ## Run 28 stress 0.08104907 
    ## Run 29 stress 0.08003196 
    ## Run 30 stress 0.08182126 
    ## Run 31 stress 0.07669184 
    ## ... Procrustes: rmse 0.0003152522  max resid 0.0005411205 
    ## ... Similar to previous best
    ## Run 32 stress 0.07943326 
    ## Run 33 stress 0.0810492 
    ## Run 34 stress 0.08012972 
    ## Run 35 stress 0.08175695 
    ## Run 36 stress 0.08464769 
    ## Run 37 stress 0.08835189 
    ## Run 38 stress 0.08471583 
    ## Run 39 stress 0.07669171 
    ## ... Procrustes: rmse 0.0001649285  max resid 0.0002806262 
    ## ... Similar to previous best
    ## Run 40 stress 0.08433696 
    ## Run 41 stress 0.08329746 
    ## Run 42 stress 0.0766918 
    ## ... Procrustes: rmse 0.0003626531  max resid 0.0006232961 
    ## ... Similar to previous best
    ## Run 43 stress 0.08464766 
    ## Run 44 stress 0.08329766 
    ## Run 45 stress 0.08003184 
    ## Run 46 stress 0.07669178 
    ## ... Procrustes: rmse 0.000349883  max resid 0.0006135346 
    ## ... Similar to previous best
    ## Run 47 stress 0.08305317 
    ## Run 48 stress 0.07731538 
    ## Run 49 stress 0.08684077 
    ## Run 50 stress 0.07868337 
    ## Run 51 stress 0.09160885 
    ## Run 52 stress 0.08459621 
    ## Run 53 stress 0.08337145 
    ## Run 54 stress 0.08175688 
    ## Run 55 stress 0.08252761 
    ## Run 56 stress 0.08265373 
    ## Run 57 stress 0.08433703 
    ## Run 58 stress 0.08684122 
    ## Run 59 stress 0.08375752 
    ## Run 60 stress 0.08175693 
    ## Run 61 stress 0.07669172 
    ## ... Procrustes: rmse 0.0002316399  max resid 0.0003965573 
    ## ... Similar to previous best
    ## Run 62 stress 0.08274939 
    ## Run 63 stress 0.08433712 
    ## Run 64 stress 0.08250108 
    ## Run 65 stress 0.08337175 
    ## Run 66 stress 0.08433684 
    ## Run 67 stress 0.08305364 
    ## Run 68 stress 0.08684136 
    ## Run 69 stress 0.08252134 
    ## Run 70 stress 0.07731514 
    ## Run 71 stress 0.08175691 
    ## Run 72 stress 0.08175692 
    ## Run 73 stress 0.08104921 
    ## Run 74 stress 0.08527034 
    ## Run 75 stress 0.08283547 
    ## Run 76 stress 0.08252126 
    ## Run 77 stress 0.08233203 
    ## Run 78 stress 0.07960662 
    ## Run 79 stress 0.08929098 
    ## Run 80 stress 0.08459587 
    ## Run 81 stress 0.0766916 
    ## ... Procrustes: rmse 0.0002025148  max resid 0.0003505464 
    ## ... Similar to previous best
    ## Run 82 stress 0.08104884 
    ## Run 83 stress 0.08305328 
    ## Run 84 stress 0.08012988 
    ## Run 85 stress 0.07868364 
    ## Run 86 stress 0.08339545 
    ## Run 87 stress 0.2358287 
    ## Run 88 stress 0.07731512 
    ## Run 89 stress 0.08362753 
    ## Run 90 stress 0.07948967 
    ## Run 91 stress 0.0809776 
    ## Run 92 stress 0.08175715 
    ## Run 93 stress 0.0809774 
    ## Run 94 stress 0.08433685 
    ## Run 95 stress 0.07871741 
    ## Run 96 stress 0.08144915 
    ## Run 97 stress 0.07960674 
    ## Run 98 stress 0.08360262 
    ## Run 99 stress 0.08175715 
    ## Run 100 stress 0.08317311 
    ## Run 101 stress 0.08835168 
    ## Run 102 stress 0.08175706 
    ## Run 103 stress 0.08433703 
    ## Run 104 stress 0.08175701 
    ## Run 105 stress 0.0773152 
    ## Run 106 stress 0.08182133 
    ## Run 107 stress 0.08274889 
    ## Run 108 stress 0.08097747 
    ## Run 109 stress 0.07731544 
    ## Run 110 stress 0.08013768 
    ## Run 111 stress 0.08249603 
    ## Run 112 stress 0.08003212 
    ## Run 113 stress 0.07669157 
    ## ... Procrustes: rmse 0.0001270029  max resid 0.0002208032 
    ## ... Similar to previous best
    ## Run 114 stress 0.0868275 
    ## Run 115 stress 0.08433691 
    ## Run 116 stress 0.08097745 
    ## Run 117 stress 0.08835135 
    ## Run 118 stress 0.07669165 
    ## ... Procrustes: rmse 0.0001988873  max resid 0.0004030401 
    ## ... Similar to previous best
    ## Run 119 stress 0.08375708 
    ## Run 120 stress 0.08012975 
    ## Run 121 stress 0.07871759 
    ## Run 122 stress 0.08012973 
    ## Run 123 stress 0.08144918 
    ## Run 124 stress 0.07669154 
    ## ... Procrustes: rmse 0.0001450294  max resid 0.0002495002 
    ## ... Similar to previous best
    ## Run 125 stress 0.08144956 
    ## Run 126 stress 0.0837033 
    ## Run 127 stress 0.07669178 
    ## ... Procrustes: rmse 0.0002878611  max resid 0.0005001415 
    ## ... Similar to previous best
    ## Run 128 stress 0.07868344 
    ## Run 129 stress 0.08012951 
    ## Run 130 stress 0.08305385 
    ## Run 131 stress 0.07868391 
    ## Run 132 stress 0.08241498 
    ## Run 133 stress 0.07669163 
    ## ... Procrustes: rmse 0.0002544012  max resid 0.000443432 
    ## ... Similar to previous best
    ## Run 134 stress 0.08876495 
    ## Run 135 stress 0.0809774 
    ## Run 136 stress 0.08104906 
    ## Run 137 stress 0.0821734 
    ## Run 138 stress 0.08104894 
    ## Run 139 stress 0.08175718 
    ## Run 140 stress 0.08144961 
    ## Run 141 stress 0.08175691 
    ## Run 142 stress 0.07669168 
    ## ... Procrustes: rmse 0.0002863095  max resid 0.0005004094 
    ## ... Similar to previous best
    ## Run 143 stress 0.07731539 
    ## Run 144 stress 0.08012976 
    ## Run 145 stress 0.08583941 
    ## Run 146 stress 0.08265431 
    ## Run 147 stress 0.07731528 
    ## Run 148 stress 0.09012168 
    ## Run 149 stress 0.08684117 
    ## Run 150 stress 0.08189387 
    ## Run 151 stress 0.07669154 
    ## ... Procrustes: rmse 0.0001626511  max resid 0.0002817551 
    ## ... Similar to previous best
    ## Run 152 stress 0.08104914 
    ## Run 153 stress 0.08459576 
    ## Run 154 stress 0.07871728 
    ## Run 155 stress 0.0818942 
    ## Run 156 stress 0.08175711 
    ## Run 157 stress 0.07949 
    ## Run 158 stress 0.08097742 
    ## Run 159 stress 0.07871726 
    ## Run 160 stress 0.08485109 
    ## Run 161 stress 0.07731519 
    ## Run 162 stress 0.08329767 
    ## Run 163 stress 0.08957587 
    ## Run 164 stress 0.08485121 
    ## Run 165 stress 0.08471018 
    ## Run 166 stress 0.07949028 
    ## Run 167 stress 0.08182123 
    ## Run 168 stress 0.08003201 
    ## Run 169 stress 0.08252116 
    ## Run 170 stress 0.08243379 
    ## Run 171 stress 0.08684363 
    ## Run 172 stress 0.08104909 
    ## Run 173 stress 0.08175697 
    ## Run 174 stress 0.08360275 
    ## Run 175 stress 0.07731518 
    ## Run 176 stress 0.09063327 
    ## Run 177 stress 0.08252141 
    ## Run 178 stress 0.08003216 
    ## Run 179 stress 0.07731533 
    ## Run 180 stress 0.0766915 
    ## ... New best solution
    ## ... Procrustes: rmse 1.481726e-05  max resid 2.663697e-05 
    ## ... Similar to previous best
    ## Run 181 stress 0.08510898 
    ## Run 182 stress 0.0931204 
    ## Run 183 stress 0.08249612 
    ## Run 184 stress 0.08500374 
    ## Run 185 stress 0.0831734 
    ## Run 186 stress 0.08097749 
    ## Run 187 stress 0.0844324 
    ## Run 188 stress 0.08305541 
    ## Run 189 stress 0.08421033 
    ## Run 190 stress 0.08506981 
    ## Run 191 stress 0.07951485 
    ## Run 192 stress 0.08355818 
    ## Run 193 stress 0.08233228 
    ## Run 194 stress 0.0845964 
    ## Run 195 stress 0.08233235 
    ## Run 196 stress 0.08218526 
    ## Run 197 stress 0.08360269 
    ## Run 198 stress 0.08355838 
    ## Run 199 stress 0.07669153 
    ## ... Procrustes: rmse 9.327998e-05  max resid 0.0001618292 
    ## ... Similar to previous best
    ## Run 200 stress 0.07871763 
    ## Run 201 stress 0.08175688 
    ## Run 202 stress 0.08443288 
    ## Run 203 stress 0.08265369 
    ## Run 204 stress 0.08252119 
    ## Run 205 stress 0.07669152 
    ## ... Procrustes: rmse 8.555821e-05  max resid 0.00018104 
    ## ... Similar to previous best
    ## Run 206 stress 0.0814489 
    ## Run 207 stress 0.08507056 
    ## Run 208 stress 0.0836026 
    ## Run 209 stress 0.08175705 
    ## Run 210 stress 0.07669174 
    ## ... Procrustes: rmse 0.000274852  max resid 0.0004757507 
    ## ... Similar to previous best
    ## Run 211 stress 0.08097755 
    ## Run 212 stress 0.07871759 
    ## Run 213 stress 0.0824337 
    ## Run 214 stress 0.08012962 
    ## Run 215 stress 0.07960693 
    ## Run 216 stress 0.08243369 
    ## Run 217 stress 0.08433687 
    ## Run 218 stress 0.08241486 
    ## Run 219 stress 0.08957523 
    ## Run 220 stress 0.0818211 
    ## Run 221 stress 0.09022524 
    ## Run 222 stress 0.08379855 
    ## Run 223 stress 0.07669154 
    ## ... Procrustes: rmse 9.094095e-05  max resid 0.0001573687 
    ## ... Similar to previous best
    ## Run 224 stress 0.08684121 
    ## Run 225 stress 0.08375709 
    ## Run 226 stress 0.08229486 
    ## Run 227 stress 0.08527046 
    ## Run 228 stress 0.08317321 
    ## Run 229 stress 0.08527021 
    ## Run 230 stress 0.08717117 
    ## Run 231 stress 0.08421038 
    ## Run 232 stress 0.08175701 
    ## Run 233 stress 0.07669175 
    ## ... Procrustes: rmse 0.0003219003  max resid 0.00055425 
    ## ... Similar to previous best
    ## Run 234 stress 0.08337123 
    ## Run 235 stress 0.07871743 
    ## Run 236 stress 0.08527016 
    ## Run 237 stress 0.08097759 
    ## Run 238 stress 0.08527025 
    ## Run 239 stress 0.08530144 
    ## Run 240 stress 0.07943375 
    ## Run 241 stress 0.07948946 
    ## Run 242 stress 0.08337178 
    ## Run 243 stress 0.08360268 
    ## Run 244 stress 0.0841857 
    ## Run 245 stress 0.08556508 
    ## Run 246 stress 0.08274972 
    ## Run 247 stress 0.08175704 
    ## Run 248 stress 0.08459691 
    ## Run 249 stress 0.07669158 
    ## ... Procrustes: rmse 0.0001147921  max resid 0.0001716341 
    ## ... Similar to previous best
    ## Run 250 stress 0.0830535 
    ## Run 251 stress 0.0848514 
    ## Run 252 stress 0.08104912 
    ## Run 253 stress 0.08493084 
    ## Run 254 stress 0.08422754 
    ## Run 255 stress 0.08012949 
    ## Run 256 stress 0.08104923 
    ## Run 257 stress 0.08265435 
    ## Run 258 stress 0.08241489 
    ## Run 259 stress 0.08097753 
    ## Run 260 stress 0.08104892 
    ## Run 261 stress 0.07669162 
    ## ... Procrustes: rmse 0.0001802921  max resid 0.0003114834 
    ## ... Similar to previous best
    ## Run 262 stress 0.07871753 
    ## Run 263 stress 0.08104919 
    ## Run 264 stress 0.08175694 
    ## Run 265 stress 0.0826535 
    ## Run 266 stress 0.07871747 
    ## Run 267 stress 0.08493102 
    ## Run 268 stress 0.08329773 
    ## Run 269 stress 0.08097742 
    ## Run 270 stress 0.08360268 
    ## Run 271 stress 0.08337155 
    ## Run 272 stress 0.08144922 
    ## Run 273 stress 0.08182113 
    ## Run 274 stress 0.08305444 
    ## Run 275 stress 0.08355803 
    ## Run 276 stress 0.08003176 
    ## Run 277 stress 0.0830553 
    ## Run 278 stress 0.08459588 
    ## Run 279 stress 0.08243363 
    ## Run 280 stress 0.08464787 
    ## Run 281 stress 0.0766917 
    ## ... Procrustes: rmse 0.0001388603  max resid 0.0002435332 
    ## ... Similar to previous best
    ## Run 282 stress 0.08337165 
    ## Run 283 stress 0.0843369 
    ## Run 284 stress 0.08175695 
    ## Run 285 stress 0.08182121 
    ## Run 286 stress 0.08144943 
    ## Run 287 stress 0.081757 
    ## Run 288 stress 0.07871743 
    ## Run 289 stress 0.08265429 
    ## Run 290 stress 0.08485117 
    ## Run 291 stress 0.08175707 
    ## Run 292 stress 0.08233218 
    ## Run 293 stress 0.08252762 
    ## Run 294 stress 0.07669159 
    ## ... Procrustes: rmse 0.0001523032  max resid 0.0002647501 
    ## ... Similar to previous best
    ## Run 295 stress 0.08337193 
    ## Run 296 stress 0.08265416 
    ## Run 297 stress 0.08443441 
    ## Run 298 stress 0.07948983 
    ## Run 299 stress 0.08375704 
    ## Run 300 stress 0.08265378 
    ## Run 301 stress 0.08182138 
    ## Run 302 stress 0.08487575 
    ## Run 303 stress 0.0801299 
    ## Run 304 stress 0.0787176 
    ## Run 305 stress 0.08957616 
    ## Run 306 stress 0.08265437 
    ## Run 307 stress 0.08553056 
    ## Run 308 stress 0.07669172 
    ## ... Procrustes: rmse 0.0003087239  max resid 0.0005329573 
    ## ... Similar to previous best
    ## Run 309 stress 0.07669155 
    ## ... Procrustes: rmse 0.0001638961  max resid 0.0003355308 
    ## ... Similar to previous best
    ## Run 310 stress 0.08443237 
    ## Run 311 stress 0.0818939 
    ## Run 312 stress 0.08104922 
    ## Run 313 stress 0.0824961 
    ## Run 314 stress 0.08433722 
    ## Run 315 stress 0.08265461 
    ## Run 316 stress 0.08144897 
    ## Run 317 stress 0.08375764 
    ## Run 318 stress 0.0786833 
    ## Run 319 stress 0.07669173 
    ## ... Procrustes: rmse 0.0003091489  max resid 0.000623477 
    ## ... Similar to previous best
    ## Run 320 stress 0.08305423 
    ## Run 321 stress 0.08250159 
    ## Run 322 stress 0.08684143 
    ## Run 323 stress 0.08337192 
    ## Run 324 stress 0.08013612 
    ## Run 325 stress 0.07960666 
    ## Run 326 stress 0.08283591 
    ## Run 327 stress 0.07868337 
    ## Run 328 stress 0.08097743 
    ## Run 329 stress 0.08252118 
    ## Run 330 stress 0.08433693 
    ## Run 331 stress 0.08097747 
    ## Run 332 stress 0.08144802 
    ## Run 333 stress 0.08013356 
    ## Run 334 stress 0.08144925 
    ## Run 335 stress 0.08337147 
    ## Run 336 stress 0.08013582 
    ## Run 337 stress 0.08175695 
    ## Run 338 stress 0.07669163 
    ## ... Procrustes: rmse 0.0001982094  max resid 0.0003469549 
    ## ... Similar to previous best
    ## Run 339 stress 0.08355803 
    ## Run 340 stress 0.0845964 
    ## Run 341 stress 0.08422739 
    ## Run 342 stress 0.08433685 
    ## Run 343 stress 0.0833719 
    ## Run 344 stress 0.08360273 
    ## Run 345 stress 0.08104884 
    ## Run 346 stress 0.07943323 
    ## Run 347 stress 0.08097742 
    ## Run 348 stress 0.0795149 
    ## Run 349 stress 0.08104893 
    ## Run 350 stress 0.07669172 
    ## ... Procrustes: rmse 0.0002262628  max resid 0.0004361037 
    ## ... Similar to previous best
    ## Run 351 stress 0.07871728 
    ## Run 352 stress 0.08283573 
    ## Run 353 stress 0.08500361 
    ## Run 354 stress 0.08182116 
    ## Run 355 stress 0.08218514 
    ## Run 356 stress 0.07943377 
    ## Run 357 stress 0.07669191 
    ## ... Procrustes: rmse 0.0003243699  max resid 0.0005626297 
    ## ... Similar to previous best
    ## Run 358 stress 0.07731517 
    ## Run 359 stress 0.08252131 
    ## Run 360 stress 0.08615948 
    ## Run 361 stress 0.0830553 
    ## Run 362 stress 0.07960666 
    ## Run 363 stress 0.08175705 
    ## Run 364 stress 0.08104902 
    ## Run 365 stress 0.08682654 
    ## Run 366 stress 0.08003185 
    ## Run 367 stress 0.07948984 
    ## Run 368 stress 0.08665279 
    ## Run 369 stress 0.07871759 
    ## Run 370 stress 0.08329764 
    ## Run 371 stress 0.08144967 
    ## Run 372 stress 0.08175698 
    ## Run 373 stress 0.08243366 
    ## Run 374 stress 0.08433703 
    ## Run 375 stress 0.07669168 
    ## ... Procrustes: rmse 0.0002131448  max resid 0.000366856 
    ## ... Similar to previous best
    ## Run 376 stress 0.08097763 
    ## Run 377 stress 0.08182107 
    ## Run 378 stress 0.08684082 
    ## Run 379 stress 0.08418576 
    ## Run 380 stress 0.08337185 
    ## Run 381 stress 0.08553063 
    ## Run 382 stress 0.08305568 
    ## Run 383 stress 0.08337154 
    ## Run 384 stress 0.08421034 
    ## Run 385 stress 0.08175695 
    ## Run 386 stress 0.08003211 
    ## Run 387 stress 0.07669161 
    ## ... Procrustes: rmse 0.0001823672  max resid 0.0003169827 
    ## ... Similar to previous best
    ## Run 388 stress 0.08360258 
    ## Run 389 stress 0.08274931 
    ## Run 390 stress 0.08360264 
    ## Run 391 stress 0.08375752 
    ## Run 392 stress 0.08249623 
    ## Run 393 stress 0.07871727 
    ## Run 394 stress 0.07960682 
    ## Run 395 stress 0.0821734 
    ## Run 396 stress 0.08104884 
    ## Run 397 stress 0.08012988 
    ## Run 398 stress 0.07948982 
    ## Run 399 stress 0.08265398 
    ## Run 400 stress 0.08243393 
    ## Run 401 stress 0.08317309 
    ## Run 402 stress 0.08950023 
    ## Run 403 stress 0.08097742 
    ## Run 404 stress 0.08249602 
    ## Run 405 stress 0.08305358 
    ## Run 406 stress 0.08418564 
    ## Run 407 stress 0.08104884 
    ## Run 408 stress 0.07871733 
    ## Run 409 stress 0.08485134 
    ## Run 410 stress 0.08175694 
    ## Run 411 stress 0.07731518 
    ## Run 412 stress 0.08217351 
    ## Run 413 stress 0.08305527 
    ## Run 414 stress 0.0868407 
    ## Run 415 stress 0.07669163 
    ## ... Procrustes: rmse 0.000197776  max resid 0.000342889 
    ## ... Similar to previous best
    ## Run 416 stress 0.07669159 
    ## ... Procrustes: rmse 0.0002004637  max resid 0.0003739402 
    ## ... Similar to previous best
    ## Run 417 stress 0.08252142 
    ## Run 418 stress 0.08362755 
    ## Run 419 stress 0.08250197 
    ## Run 420 stress 0.07669164 
    ## ... Procrustes: rmse 0.0002332753  max resid 0.0004088819 
    ## ... Similar to previous best
    ## Run 421 stress 0.08144924 
    ## Run 422 stress 0.07669151 
    ## ... Procrustes: rmse 5.121433e-05  max resid 9.587754e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.08530182 
    ## Run 424 stress 0.08104922 
    ## Run 425 stress 0.08189404 
    ## Run 426 stress 0.08317348 
    ## Run 427 stress 0.07669151 
    ## ... Procrustes: rmse 0.0001051372  max resid 0.0001834106 
    ## ... Similar to previous best
    ## Run 428 stress 0.08144922 
    ## Run 429 stress 0.0794905 
    ## Run 430 stress 0.07871736 
    ## Run 431 stress 0.08527024 
    ## Run 432 stress 0.08218532 
    ## Run 433 stress 0.08433693 
    ## Run 434 stress 0.081049 
    ## Run 435 stress 0.08182106 
    ## Run 436 stress 0.08360257 
    ## Run 437 stress 0.07951489 
    ## Run 438 stress 0.08360272 
    ## Run 439 stress 0.08360271 
    ## Run 440 stress 0.07868326 
    ## Run 441 stress 0.08584015 
    ## Run 442 stress 0.09070247 
    ## Run 443 stress 0.08097768 
    ## Run 444 stress 0.07669167 
    ## ... Procrustes: rmse 0.00026643  max resid 0.0004589983 
    ## ... Similar to previous best
    ## Run 445 stress 0.0826532 
    ## Run 446 stress 0.07871726 
    ## Run 447 stress 0.0854727 
    ## Run 448 stress 0.08553065 
    ## Run 449 stress 0.08382979 
    ## Run 450 stress 0.08443298 
    ## Run 451 stress 0.07960659 
    ## Run 452 stress 0.08097762 
    ## Run 453 stress 0.08242512 
    ## Run 454 stress 0.09028691 
    ## Run 455 stress 0.0794337 
    ## Run 456 stress 0.08684109 
    ## Run 457 stress 0.08957603 
    ## Run 458 stress 0.08684121 
    ## Run 459 stress 0.07951477 
    ## Run 460 stress 0.08175691 
    ## Run 461 stress 0.08370297 
    ## Run 462 stress 0.0830541 
    ## Run 463 stress 0.0817569 
    ## Run 464 stress 0.08003173 
    ## Run 465 stress 0.08927527 
    ## Run 466 stress 0.082415 
    ## Run 467 stress 0.08428762 
    ## Run 468 stress 0.08283536 
    ## Run 469 stress 0.08175702 
    ## Run 470 stress 0.08443413 
    ## Run 471 stress 0.08835132 
    ## Run 472 stress 0.07948968 
    ## Run 473 stress 0.08514309 
    ## Run 474 stress 0.08477002 
    ## Run 475 stress 0.08317331 
    ## Run 476 stress 0.08182123 
    ## Run 477 stress 0.08314416 
    ## Run 478 stress 0.08684185 
    ## Run 479 stress 0.08317355 
    ## Run 480 stress 0.08249604 
    ## Run 481 stress 0.08175689 
    ## Run 482 stress 0.0809777 
    ## Run 483 stress 0.08104914 
    ## Run 484 stress 0.07948996 
    ## Run 485 stress 0.07951508 
    ## Run 486 stress 0.08003206 
    ## Run 487 stress 0.07868384 
    ## Run 488 stress 0.08233221 
    ## Run 489 stress 0.08421001 
    ## Run 490 stress 0.08583933 
    ## Run 491 stress 0.08317315 
    ## Run 492 stress 0.07669166 
    ## ... Procrustes: rmse 0.0001525027  max resid 0.0002905296 
    ## ... Similar to previous best
    ## Run 493 stress 0.08360258 
    ## Run 494 stress 0.082496 
    ## Run 495 stress 0.08013015 
    ## Run 496 stress 0.07943369 
    ## Run 497 stress 0.0827493 
    ## Run 498 stress 0.07731513 
    ## Run 499 stress 0.08337113 
    ## Run 500 stress 0.08968588 
    ## *** Best solution repeated 25 times

``` r
# Mixed and stratified lakes
SD_beta_env_MS_NMDS <- metaMDS(SD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09445477 
    ## Run 2 stress 0.09286086 
    ## Run 3 stress 0.0976084 
    ## Run 4 stress 0.2776223 
    ## Run 5 stress 0.0940711 
    ## Run 6 stress 0.08773471 
    ## Run 7 stress 0.09337276 
    ## Run 8 stress 0.08503506 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001768205  max resid 0.004411504 
    ## ... Similar to previous best
    ## Run 9 stress 0.09159083 
    ## Run 10 stress 0.08773484 
    ## Run 11 stress 0.09400534 
    ## Run 12 stress 0.09030395 
    ## Run 13 stress 0.08440276 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01420246  max resid 0.04326642 
    ## Run 14 stress 0.09760826 
    ## Run 15 stress 0.09374291 
    ## Run 16 stress 0.090304 
    ## Run 17 stress 0.1038083 
    ## Run 18 stress 0.09408012 
    ## Run 19 stress 0.09030394 
    ## Run 20 stress 0.09374232 
    ## Run 21 stress 0.0971296 
    ## Run 22 stress 0.09407995 
    ## Run 23 stress 0.1038049 
    ## Run 24 stress 0.08440254 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001873013  max resid 0.0003240579 
    ## ... Similar to previous best
    ## Run 25 stress 0.09030412 
    ## Run 26 stress 0.09408001 
    ## Run 27 stress 0.09969983 
    ## Run 28 stress 0.1006348 
    ## Run 29 stress 0.08503512 
    ## Run 30 stress 0.09590489 
    ## Run 31 stress 0.09159084 
    ## Run 32 stress 0.0850365 
    ## Run 33 stress 0.09168934 
    ## Run 34 stress 0.09268376 
    ## Run 35 stress 0.09030399 
    ## Run 36 stress 0.09286085 
    ## Run 37 stress 0.09030393 
    ## Run 38 stress 0.08440254 
    ## ... New best solution
    ## ... Procrustes: rmse 2.72139e-05  max resid 5.319013e-05 
    ## ... Similar to previous best
    ## Run 39 stress 0.09380587 
    ## Run 40 stress 0.09970011 
    ## Run 41 stress 0.3159104 
    ## Run 42 stress 0.0940797 
    ## Run 43 stress 0.09030394 
    ## Run 44 stress 0.08773471 
    ## Run 45 stress 0.09447337 
    ## Run 46 stress 0.09590952 
    ## Run 47 stress 0.09308975 
    ## Run 48 stress 0.09337292 
    ## Run 49 stress 0.09168944 
    ## Run 50 stress 0.08973865 
    ## Run 51 stress 0.09030397 
    ## Run 52 stress 0.0850347 
    ## Run 53 stress 0.2709182 
    ## Run 54 stress 0.09969982 
    ## Run 55 stress 0.09712994 
    ## Run 56 stress 0.09380588 
    ## Run 57 stress 0.09337229 
    ## Run 58 stress 0.09308968 
    ## Run 59 stress 0.09407977 
    ## Run 60 stress 0.08503477 
    ## Run 61 stress 0.09337232 
    ## Run 62 stress 0.1038039 
    ## Run 63 stress 0.09321496 
    ## Run 64 stress 0.09380584 
    ## Run 65 stress 0.08773479 
    ## Run 66 stress 0.09712966 
    ## Run 67 stress 0.08773474 
    ## Run 68 stress 0.08773475 
    ## Run 69 stress 0.09030403 
    ## Run 70 stress 0.08973871 
    ## Run 71 stress 0.09308973 
    ## Run 72 stress 0.09030412 
    ## Run 73 stress 0.08503464 
    ## Run 74 stress 0.08503656 
    ## Run 75 stress 0.09374285 
    ## Run 76 stress 0.08773467 
    ## Run 77 stress 0.09030404 
    ## Run 78 stress 0.09407975 
    ## Run 79 stress 0.09286087 
    ## Run 80 stress 0.0903041 
    ## Run 81 stress 0.08773482 
    ## Run 82 stress 0.09030397 
    ## Run 83 stress 0.09400515 
    ## Run 84 stress 0.08440254 
    ## ... New best solution
    ## ... Procrustes: rmse 1.309205e-05  max resid 3.286076e-05 
    ## ... Similar to previous best
    ## Run 85 stress 0.08973887 
    ## Run 86 stress 0.09168944 
    ## Run 87 stress 0.09159095 
    ## Run 88 stress 0.09168932 
    ## Run 89 stress 0.09464391 
    ## Run 90 stress 0.09337235 
    ## Run 91 stress 0.09145335 
    ## Run 92 stress 0.09447296 
    ## Run 93 stress 0.08440255 
    ## ... Procrustes: rmse 2.673996e-05  max resid 5.464605e-05 
    ## ... Similar to previous best
    ## Run 94 stress 0.09407975 
    ## Run 95 stress 0.09337316 
    ## Run 96 stress 0.09407986 
    ## Run 97 stress 0.08440263 
    ## ... Procrustes: rmse 0.0002719865  max resid 0.0005382955 
    ## ... Similar to previous best
    ## Run 98 stress 0.08503478 
    ## Run 99 stress 0.08773479 
    ## Run 100 stress 0.105285 
    ## Run 101 stress 0.09465908 
    ## Run 102 stress 0.08503509 
    ## Run 103 stress 0.08773468 
    ## Run 104 stress 0.09030397 
    ## Run 105 stress 0.09416131 
    ## Run 106 stress 0.08773465 
    ## Run 107 stress 0.1006332 
    ## Run 108 stress 0.09417779 
    ## Run 109 stress 0.0933729 
    ## Run 110 stress 0.08973881 
    ## Run 111 stress 0.09374253 
    ## Run 112 stress 0.09465904 
    ## Run 113 stress 0.2461313 
    ## Run 114 stress 0.09445485 
    ## Run 115 stress 0.09404382 
    ## Run 116 stress 0.09030394 
    ## Run 117 stress 0.09416132 
    ## Run 118 stress 0.09403397 
    ## Run 119 stress 0.09159086 
    ## Run 120 stress 0.08503461 
    ## Run 121 stress 0.08503541 
    ## Run 122 stress 0.08503483 
    ## Run 123 stress 0.09268366 
    ## Run 124 stress 0.09030415 
    ## Run 125 stress 0.09760834 
    ## Run 126 stress 0.09159087 
    ## Run 127 stress 0.09374258 
    ## Run 128 stress 0.2466508 
    ## Run 129 stress 0.2510692 
    ## Run 130 stress 0.09268369 
    ## Run 131 stress 0.09286103 
    ## Run 132 stress 0.09464509 
    ## Run 133 stress 0.08503525 
    ## Run 134 stress 0.09268319 
    ## Run 135 stress 0.09168938 
    ## Run 136 stress 0.09030397 
    ## Run 137 stress 0.08773471 
    ## Run 138 stress 0.08440268 
    ## ... Procrustes: rmse 0.000346989  max resid 0.0006878923 
    ## ... Similar to previous best
    ## Run 139 stress 0.09145322 
    ## Run 140 stress 0.09539199 
    ## Run 141 stress 0.08503569 
    ## Run 142 stress 0.08773486 
    ## Run 143 stress 0.0940798 
    ## Run 144 stress 0.09286084 
    ## Run 145 stress 0.09407108 
    ## Run 146 stress 0.08440263 
    ## ... Procrustes: rmse 0.0001131374  max resid 0.0001900718 
    ## ... Similar to previous best
    ## Run 147 stress 0.09030398 
    ## Run 148 stress 0.09380576 
    ## Run 149 stress 0.093743 
    ## Run 150 stress 0.08973862 
    ## Run 151 stress 0.09407983 
    ## Run 152 stress 0.08973879 
    ## Run 153 stress 0.09407984 
    ## Run 154 stress 0.09465905 
    ## Run 155 stress 0.08503479 
    ## Run 156 stress 0.09407994 
    ## Run 157 stress 0.09286095 
    ## Run 158 stress 0.1038057 
    ## Run 159 stress 0.09030395 
    ## Run 160 stress 0.0897387 
    ## Run 161 stress 0.08503605 
    ## Run 162 stress 0.09400551 
    ## Run 163 stress 0.09030397 
    ## Run 164 stress 0.0996997 
    ## Run 165 stress 0.08440262 
    ## ... Procrustes: rmse 9.205571e-05  max resid 0.0001664642 
    ## ... Similar to previous best
    ## Run 166 stress 0.08503493 
    ## Run 167 stress 0.08973862 
    ## Run 168 stress 0.09465904 
    ## Run 169 stress 0.09465909 
    ## Run 170 stress 0.1054526 
    ## Run 171 stress 0.09030396 
    ## Run 172 stress 0.08440253 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001908463  max resid 0.0003737926 
    ## ... Similar to previous best
    ## Run 173 stress 0.08973863 
    ## Run 174 stress 0.08503492 
    ## Run 175 stress 0.100461 
    ## Run 176 stress 0.1038047 
    ## Run 177 stress 0.09268392 
    ## Run 178 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000120922  max resid 0.000222387 
    ## ... Similar to previous best
    ## Run 179 stress 0.09416134 
    ## Run 180 stress 0.08773482 
    ## Run 181 stress 0.1038039 
    ## Run 182 stress 0.09969966 
    ## Run 183 stress 0.09159085 
    ## Run 184 stress 0.09407969 
    ## Run 185 stress 0.09374219 
    ## Run 186 stress 0.09308982 
    ## Run 187 stress 0.08503627 
    ## Run 188 stress 0.09145324 
    ## Run 189 stress 0.09381849 
    ## Run 190 stress 0.100461 
    ## Run 191 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001956266  max resid 0.0003756378 
    ## ... Similar to previous best
    ## Run 192 stress 0.09445509 
    ## Run 193 stress 0.09337274 
    ## Run 194 stress 0.09168926 
    ## Run 195 stress 0.09465905 
    ## Run 196 stress 0.08973863 
    ## Run 197 stress 0.09268398 
    ## Run 198 stress 0.09539192 
    ## Run 199 stress 0.09159084 
    ## Run 200 stress 0.08440265 
    ## ... Procrustes: rmse 0.0002444149  max resid 0.0004784581 
    ## ... Similar to previous best
    ## Run 201 stress 0.08440264 
    ## ... Procrustes: rmse 0.000189033  max resid 0.0003949592 
    ## ... Similar to previous best
    ## Run 202 stress 0.08503643 
    ## Run 203 stress 0.09159084 
    ## Run 204 stress 0.09416131 
    ## Run 205 stress 0.08503533 
    ## Run 206 stress 0.09712978 
    ## Run 207 stress 0.09969969 
    ## Run 208 stress 0.09030402 
    ## Run 209 stress 0.08973881 
    ## Run 210 stress 0.0850347 
    ## Run 211 stress 0.08440255 
    ## ... Procrustes: rmse 3.445424e-05  max resid 7.982048e-05 
    ## ... Similar to previous best
    ## Run 212 stress 0.09590897 
    ## Run 213 stress 0.0903041 
    ## Run 214 stress 0.09465914 
    ## Run 215 stress 0.09286083 
    ## Run 216 stress 0.09145324 
    ## Run 217 stress 0.08440261 
    ## ... Procrustes: rmse 0.0001526176  max resid 0.0003159987 
    ## ... Similar to previous best
    ## Run 218 stress 0.09416131 
    ## Run 219 stress 0.09465907 
    ## Run 220 stress 0.09464461 
    ## Run 221 stress 0.09030403 
    ## Run 222 stress 0.09408008 
    ## Run 223 stress 0.09308958 
    ## Run 224 stress 0.09321767 
    ## Run 225 stress 0.09308968 
    ## Run 226 stress 0.09760836 
    ## Run 227 stress 0.09465908 
    ## Run 228 stress 0.09308955 
    ## Run 229 stress 0.09969982 
    ## Run 230 stress 0.09408004 
    ## Run 231 stress 0.08440262 
    ## ... Procrustes: rmse 0.0001532131  max resid 0.0002863219 
    ## ... Similar to previous best
    ## Run 232 stress 0.08440267 
    ## ... Procrustes: rmse 0.0001764349  max resid 0.0003331399 
    ## ... Similar to previous best
    ## Run 233 stress 0.09447321 
    ## Run 234 stress 0.1026217 
    ## Run 235 stress 0.09381862 
    ## Run 236 stress 0.09590604 
    ## Run 237 stress 0.09308974 
    ## Run 238 stress 0.09030408 
    ## Run 239 stress 0.09535401 
    ## Run 240 stress 0.0850351 
    ## Run 241 stress 0.09337281 
    ## Run 242 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001001018  max resid 0.0001993081 
    ## ... Similar to previous best
    ## Run 243 stress 0.08973866 
    ## Run 244 stress 0.08503605 
    ## Run 245 stress 0.08503492 
    ## Run 246 stress 0.09374246 
    ## Run 247 stress 0.09407974 
    ## Run 248 stress 0.09030408 
    ## Run 249 stress 0.09145323 
    ## Run 250 stress 0.0940742 
    ## Run 251 stress 0.094034 
    ## Run 252 stress 0.09407983 
    ## Run 253 stress 0.09407996 
    ## Run 254 stress 0.0937422 
    ## Run 255 stress 0.08503501 
    ## Run 256 stress 0.09447322 
    ## Run 257 stress 0.09337231 
    ## Run 258 stress 0.09286084 
    ## Run 259 stress 0.08973874 
    ## Run 260 stress 0.09908182 
    ## Run 261 stress 0.08973894 
    ## Run 262 stress 0.09447332 
    ## Run 263 stress 0.09168939 
    ## Run 264 stress 0.09268357 
    ## Run 265 stress 0.0937428 
    ## Run 266 stress 0.09407984 
    ## Run 267 stress 0.09321839 
    ## Run 268 stress 0.08503466 
    ## Run 269 stress 0.1038051 
    ## Run 270 stress 0.09380591 
    ## Run 271 stress 0.2531479 
    ## Run 272 stress 0.09465909 
    ## Run 273 stress 0.0850351 
    ## Run 274 stress 0.09465905 
    ## Run 275 stress 0.1054541 
    ## Run 276 stress 0.09445512 
    ## Run 277 stress 0.09465905 
    ## Run 278 stress 0.0933726 
    ## Run 279 stress 0.09407967 
    ## Run 280 stress 0.09286092 
    ## Run 281 stress 0.08973872 
    ## Run 282 stress 0.09374261 
    ## Run 283 stress 0.08503459 
    ## Run 284 stress 0.08973873 
    ## Run 285 stress 0.09465906 
    ## Run 286 stress 0.09030395 
    ## Run 287 stress 0.09308956 
    ## Run 288 stress 0.09159085 
    ## Run 289 stress 0.09030404 
    ## Run 290 stress 0.09168959 
    ## Run 291 stress 0.09268368 
    ## Run 292 stress 0.09407981 
    ## Run 293 stress 0.09407999 
    ## Run 294 stress 0.09464465 
    ## Run 295 stress 0.09030409 
    ## Run 296 stress 0.1084713 
    ## Run 297 stress 0.08503698 
    ## Run 298 stress 0.0877347 
    ## Run 299 stress 0.09465918 
    ## Run 300 stress 0.09286083 
    ## Run 301 stress 0.08973887 
    ## Run 302 stress 0.09969977 
    ## Run 303 stress 0.090304 
    ## Run 304 stress 0.09407986 
    ## Run 305 stress 0.09268385 
    ## Run 306 stress 0.09145323 
    ## Run 307 stress 0.08440254 
    ## ... Procrustes: rmse 7.948336e-05  max resid 0.0001554421 
    ## ... Similar to previous best
    ## Run 308 stress 0.0937426 
    ## Run 309 stress 0.09407991 
    ## Run 310 stress 0.09159096 
    ## Run 311 stress 0.09408007 
    ## Run 312 stress 0.08440266 
    ## ... Procrustes: rmse 0.0002430328  max resid 0.0005445842 
    ## ... Similar to previous best
    ## Run 313 stress 0.09407969 
    ## Run 314 stress 0.09465905 
    ## Run 315 stress 0.09159093 
    ## Run 316 stress 0.1038044 
    ## Run 317 stress 0.09407992 
    ## Run 318 stress 0.09760818 
    ## Run 319 stress 0.09337255 
    ## Run 320 stress 0.08973879 
    ## Run 321 stress 0.09539201 
    ## Run 322 stress 0.09030401 
    ## Run 323 stress 0.09407105 
    ## Run 324 stress 0.08973866 
    ## Run 325 stress 0.09465905 
    ## Run 326 stress 0.09380596 
    ## Run 327 stress 0.09337276 
    ## Run 328 stress 0.09535403 
    ## Run 329 stress 0.09030394 
    ## Run 330 stress 0.08973867 
    ## Run 331 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 3.109152e-05  max resid 8.721627e-05 
    ## ... Similar to previous best
    ## Run 332 stress 0.09407987 
    ## Run 333 stress 0.0877348 
    ## Run 334 stress 0.09381846 
    ## Run 335 stress 0.08440263 
    ## ... Procrustes: rmse 0.0001842141  max resid 0.0003260029 
    ## ... Similar to previous best
    ## Run 336 stress 0.09030395 
    ## Run 337 stress 0.09374359 
    ## Run 338 stress 0.09465905 
    ## Run 339 stress 0.09712953 
    ## Run 340 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001738108  max resid 0.0003605752 
    ## ... Similar to previous best
    ## Run 341 stress 0.09308964 
    ## Run 342 stress 0.09030397 
    ## Run 343 stress 0.0940743 
    ## Run 344 stress 0.1054528 
    ## Run 345 stress 0.0850349 
    ## Run 346 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001465006  max resid 0.000259495 
    ## ... Similar to previous best
    ## Run 347 stress 0.0915909 
    ## Run 348 stress 0.09286093 
    ## Run 349 stress 0.08440251 
    ## ... Procrustes: rmse 1.896554e-05  max resid 4.448451e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.09407993 
    ## Run 351 stress 0.09268385 
    ## Run 352 stress 0.09407985 
    ## Run 353 stress 0.09308947 
    ## Run 354 stress 0.09381858 
    ## Run 355 stress 0.08973873 
    ## Run 356 stress 0.09030393 
    ## Run 357 stress 0.09407998 
    ## Run 358 stress 0.09168953 
    ## Run 359 stress 0.09168949 
    ## Run 360 stress 0.09416131 
    ## Run 361 stress 0.08773472 
    ## Run 362 stress 0.0850355 
    ## Run 363 stress 0.08503503 
    ## Run 364 stress 0.09374328 
    ## Run 365 stress 0.09159084 
    ## Run 366 stress 0.09590685 
    ## Run 367 stress 0.09030398 
    ## Run 368 stress 0.09145321 
    ## Run 369 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001347092  max resid 0.000277533 
    ## ... Similar to previous best
    ## Run 370 stress 0.0926838 
    ## Run 371 stress 0.08503659 
    ## Run 372 stress 0.09380577 
    ## Run 373 stress 0.09465905 
    ## Run 374 stress 0.08503469 
    ## Run 375 stress 0.09380579 
    ## Run 376 stress 0.09159084 
    ## Run 377 stress 0.09168944 
    ## Run 378 stress 0.09145321 
    ## Run 379 stress 0.08973885 
    ## Run 380 stress 0.09168953 
    ## Run 381 stress 0.0915909 
    ## Run 382 stress 0.09337294 
    ## Run 383 stress 0.09030402 
    ## Run 384 stress 0.09590644 
    ## Run 385 stress 0.0946443 
    ## Run 386 stress 0.1038055 
    ## Run 387 stress 0.0930895 
    ## Run 388 stress 0.08773467 
    ## Run 389 stress 0.08440262 
    ## ... Procrustes: rmse 0.0001935104  max resid 0.0003609724 
    ## ... Similar to previous best
    ## Run 390 stress 0.09969973 
    ## Run 391 stress 0.09030395 
    ## Run 392 stress 0.09308956 
    ## Run 393 stress 0.08440269 
    ## ... Procrustes: rmse 0.0002707378  max resid 0.0005613854 
    ## ... Similar to previous best
    ## Run 394 stress 0.09308957 
    ## Run 395 stress 0.2918071 
    ## Run 396 stress 0.090304 
    ## Run 397 stress 0.09286089 
    ## Run 398 stress 0.08503677 
    ## Run 399 stress 0.090304 
    ## Run 400 stress 0.09337262 
    ## Run 401 stress 0.1026215 
    ## Run 402 stress 0.08440261 
    ## ... Procrustes: rmse 0.0001494601  max resid 0.0002607431 
    ## ... Similar to previous best
    ## Run 403 stress 0.08503502 
    ## Run 404 stress 0.09416146 
    ## Run 405 stress 0.09268402 
    ## Run 406 stress 0.09465904 
    ## Run 407 stress 0.09416133 
    ## Run 408 stress 0.09760814 
    ## Run 409 stress 0.09447315 
    ## Run 410 stress 0.0953549 
    ## Run 411 stress 0.09465918 
    ## Run 412 stress 0.08773468 
    ## Run 413 stress 0.08773471 
    ## Run 414 stress 0.09590923 
    ## Run 415 stress 0.09969987 
    ## Run 416 stress 0.09590948 
    ## Run 417 stress 0.0940343 
    ## Run 418 stress 0.09337238 
    ## Run 419 stress 0.09407106 
    ## Run 420 stress 0.1026215 
    ## Run 421 stress 0.09408005 
    ## Run 422 stress 0.09407975 
    ## Run 423 stress 0.09337237 
    ## Run 424 stress 0.08440261 
    ## ... Procrustes: rmse 0.0002035873  max resid 0.000417329 
    ## ... Similar to previous best
    ## Run 425 stress 0.105285 
    ## Run 426 stress 0.2500205 
    ## Run 427 stress 0.08503615 
    ## Run 428 stress 0.09465907 
    ## Run 429 stress 0.09465906 
    ## Run 430 stress 0.08973872 
    ## Run 431 stress 0.09159098 
    ## Run 432 stress 0.09337263 
    ## Run 433 stress 0.08440254 
    ## ... Procrustes: rmse 5.203896e-05  max resid 9.710184e-05 
    ## ... Similar to previous best
    ## Run 434 stress 0.09407969 
    ## Run 435 stress 0.09030397 
    ## Run 436 stress 0.08503485 
    ## Run 437 stress 0.09407988 
    ## Run 438 stress 0.09407414 
    ## Run 439 stress 0.09168933 
    ## Run 440 stress 0.08773481 
    ## Run 441 stress 0.08503575 
    ## Run 442 stress 0.08440267 
    ## ... Procrustes: rmse 0.0002537913  max resid 0.0005298024 
    ## ... Similar to previous best
    ## Run 443 stress 0.09337297 
    ## Run 444 stress 0.08440253 
    ## ... Procrustes: rmse 4.29376e-05  max resid 9.312742e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.09159084 
    ## Run 446 stress 0.09145328 
    ## Run 447 stress 0.09539194 
    ## Run 448 stress 0.0946591 
    ## Run 449 stress 0.09030402 
    ## Run 450 stress 0.1054545 
    ## Run 451 stress 0.09416135 
    ## Run 452 stress 0.08973868 
    ## Run 453 stress 0.08440251 
    ## ... Procrustes: rmse 3.616115e-05  max resid 6.779314e-05 
    ## ... Similar to previous best
    ## Run 454 stress 0.09337271 
    ## Run 455 stress 0.09145324 
    ## Run 456 stress 0.09159089 
    ## Run 457 stress 0.08973883 
    ## Run 458 stress 0.09374323 
    ## Run 459 stress 0.09407131 
    ## Run 460 stress 0.09286095 
    ## Run 461 stress 0.09268361 
    ## Run 462 stress 0.08503632 
    ## Run 463 stress 0.0844026 
    ## ... Procrustes: rmse 0.0001675473  max resid 0.0002970864 
    ## ... Similar to previous best
    ## Run 464 stress 0.1052849 
    ## Run 465 stress 0.09159101 
    ## Run 466 stress 0.09268322 
    ## Run 467 stress 0.0850363 
    ## Run 468 stress 0.08503597 
    ## Run 469 stress 0.09374323 
    ## Run 470 stress 0.091591 
    ## Run 471 stress 0.09337279 
    ## Run 472 stress 0.09374293 
    ## Run 473 stress 0.09407973 
    ## Run 474 stress 0.08503475 
    ## Run 475 stress 0.09539195 
    ## Run 476 stress 0.09159086 
    ## Run 477 stress 0.08440251 
    ## ... Procrustes: rmse 4.154655e-05  max resid 8.740569e-05 
    ## ... Similar to previous best
    ## Run 478 stress 0.08973863 
    ## Run 479 stress 0.09445486 
    ## Run 480 stress 0.09308946 
    ## Run 481 stress 0.08503571 
    ## Run 482 stress 0.08973877 
    ## Run 483 stress 0.09308973 
    ## Run 484 stress 0.255569 
    ## Run 485 stress 0.08503491 
    ## Run 486 stress 0.09465905 
    ## Run 487 stress 0.09168968 
    ## Run 488 stress 0.09286085 
    ## Run 489 stress 0.09159085 
    ## Run 490 stress 0.1026215 
    ## Run 491 stress 0.09465917 
    ## Run 492 stress 0.08503512 
    ## Run 493 stress 0.09407981 
    ## Run 494 stress 0.09416145 
    ## Run 495 stress 0.09407986 
    ## Run 496 stress 0.09465909 
    ## Run 497 stress 0.09969974 
    ## Run 498 stress 0.1038047 
    ## Run 499 stress 0.1026217 
    ## Run 500 stress 0.08773476 
    ## *** Best solution repeated 16 times

``` r
# Ocean sites and mixed lakes
SD_beta_env_OM_NMDS <- metaMDS(SD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06623458 
    ## Run 1 stress 0.07411943 
    ## Run 2 stress 0.07411944 
    ## Run 3 stress 0.06123361 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1442997  max resid 0.2619055 
    ## Run 4 stress 0.3121118 
    ## Run 5 stress 0.07411945 
    ## Run 6 stress 0.06623462 
    ## Run 7 stress 0.06123361 
    ## ... Procrustes: rmse 2.933377e-05  max resid 3.933845e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.0741194 
    ## Run 9 stress 0.06623459 
    ## Run 10 stress 0.08226174 
    ## Run 11 stress 0.07411954 
    ## Run 12 stress 0.06623457 
    ## Run 13 stress 0.06124502 
    ## ... Procrustes: rmse 0.01510959  max resid 0.04184091 
    ## Run 14 stress 0.07166964 
    ## Run 15 stress 0.07411941 
    ## Run 16 stress 0.07166981 
    ## Run 17 stress 0.06623466 
    ## Run 18 stress 0.07411938 
    ## Run 19 stress 0.06623461 
    ## Run 20 stress 0.0716697 
    ## Run 21 stress 0.07166974 
    ## Run 22 stress 0.07411938 
    ## Run 23 stress 0.06623457 
    ## Run 24 stress 0.06124511 
    ## ... Procrustes: rmse 0.01518118  max resid 0.04203366 
    ## Run 25 stress 0.06477866 
    ## Run 26 stress 0.07411938 
    ## Run 27 stress 0.08226175 
    ## Run 28 stress 0.06124517 
    ## ... Procrustes: rmse 0.0152171  max resid 0.04213033 
    ## Run 29 stress 0.06477837 
    ## Run 30 stress 0.06123361 
    ## ... Procrustes: rmse 2.637192e-06  max resid 3.320319e-06 
    ## ... Similar to previous best
    ## Run 31 stress 0.0741194 
    ## Run 32 stress 0.07166971 
    ## Run 33 stress 0.06477856 
    ## Run 34 stress 0.06124493 
    ## ... Procrustes: rmse 0.0150067  max resid 0.04156569 
    ## Run 35 stress 0.06477827 
    ## Run 36 stress 0.07411942 
    ## Run 37 stress 0.06123361 
    ## ... Procrustes: rmse 5.29321e-05  max resid 6.403247e-05 
    ## ... Similar to previous best
    ## Run 38 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 1.639413e-05  max resid 2.235421e-05 
    ## ... Similar to previous best
    ## Run 39 stress 0.07166983 
    ## Run 40 stress 0.06124516 
    ## ... Procrustes: rmse 0.01521858  max resid 0.04213308 
    ## Run 41 stress 0.0662346 
    ## Run 42 stress 0.06477888 
    ## Run 43 stress 0.06477884 
    ## Run 44 stress 0.07166987 
    ## Run 45 stress 0.06124486 
    ## ... Procrustes: rmse 0.0147679  max resid 0.04092585 
    ## Run 46 stress 0.07411951 
    ## Run 47 stress 0.06477861 
    ## Run 48 stress 0.07411946 
    ## Run 49 stress 0.3155894 
    ## Run 50 stress 0.07411951 
    ## Run 51 stress 0.06477874 
    ## Run 52 stress 0.0741194 
    ## Run 53 stress 0.07411939 
    ## Run 54 stress 0.06477847 
    ## Run 55 stress 0.06623462 
    ## Run 56 stress 0.06124506 
    ## ... Procrustes: rmse 0.01513473  max resid 0.04190904 
    ## Run 57 stress 0.06623459 
    ## Run 58 stress 0.06124511 
    ## ... Procrustes: rmse 0.01436732  max resid 0.03985041 
    ## Run 59 stress 0.0741194 
    ## Run 60 stress 0.06477832 
    ## Run 61 stress 0.06623457 
    ## Run 62 stress 0.2324715 
    ## Run 63 stress 0.07166969 
    ## Run 64 stress 0.06623459 
    ## Run 65 stress 0.0647792 
    ## Run 66 stress 0.0662346 
    ## Run 67 stress 0.3005077 
    ## Run 68 stress 0.06623458 
    ## Run 69 stress 0.3113447 
    ## Run 70 stress 0.0662346 
    ## Run 71 stress 0.06623462 
    ## Run 72 stress 0.3117921 
    ## Run 73 stress 0.06623458 
    ## Run 74 stress 0.06623466 
    ## Run 75 stress 0.3157939 
    ## Run 76 stress 0.08226156 
    ## Run 77 stress 0.07411941 
    ## Run 78 stress 0.07411941 
    ## Run 79 stress 0.07411951 
    ## Run 80 stress 0.06477911 
    ## Run 81 stress 0.07411939 
    ## Run 82 stress 0.09300862 
    ## Run 83 stress 0.06123361 
    ## ... Procrustes: rmse 1.184891e-05  max resid 1.891698e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.07411938 
    ## Run 85 stress 0.06623459 
    ## Run 86 stress 0.06477955 
    ## Run 87 stress 0.06623458 
    ## Run 88 stress 0.07411943 
    ## Run 89 stress 0.0741195 
    ## Run 90 stress 0.0662346 
    ## Run 91 stress 0.06477946 
    ## Run 92 stress 0.3041507 
    ## Run 93 stress 0.08226164 
    ## Run 94 stress 0.06477876 
    ## Run 95 stress 0.06124519 
    ## ... Procrustes: rmse 0.01523152  max resid 0.04216532 
    ## Run 96 stress 0.06124512 
    ## ... Procrustes: rmse 0.01509557  max resid 0.04179936 
    ## Run 97 stress 0.06623462 
    ## Run 98 stress 0.06124511 
    ## ... Procrustes: rmse 0.01514849  max resid 0.04194734 
    ## Run 99 stress 0.06623461 
    ## Run 100 stress 0.06623466 
    ## Run 101 stress 0.2361811 
    ## Run 102 stress 0.06477968 
    ## Run 103 stress 0.06623459 
    ## Run 104 stress 0.07166975 
    ## Run 105 stress 0.07166964 
    ## Run 106 stress 0.06123361 
    ## ... Procrustes: rmse 1.791761e-05  max resid 2.31082e-05 
    ## ... Similar to previous best
    ## Run 107 stress 0.06623459 
    ## Run 108 stress 0.06623458 
    ## Run 109 stress 0.07411942 
    ## Run 110 stress 0.0741194 
    ## Run 111 stress 0.06477916 
    ## Run 112 stress 0.06623457 
    ## Run 113 stress 0.06623457 
    ## Run 114 stress 0.06477954 
    ## Run 115 stress 0.0647787 
    ## Run 116 stress 0.06124504 
    ## ... Procrustes: rmse 0.01511855  max resid 0.04186504 
    ## Run 117 stress 0.08226167 
    ## Run 118 stress 0.07166966 
    ## Run 119 stress 0.06124499 
    ## ... Procrustes: rmse 0.01507417  max resid 0.04174583 
    ## Run 120 stress 0.2381016 
    ## Run 121 stress 0.07411946 
    ## Run 122 stress 0.07166965 
    ## Run 123 stress 0.06477899 
    ## Run 124 stress 0.07411943 
    ## Run 125 stress 0.07411941 
    ## Run 126 stress 0.07411939 
    ## Run 127 stress 0.06124521 
    ## ... Procrustes: rmse 0.01523269  max resid 0.04216776 
    ## Run 128 stress 0.06623458 
    ## Run 129 stress 0.08226156 
    ## Run 130 stress 0.07166976 
    ## Run 131 stress 0.07411939 
    ## Run 132 stress 0.0612451 
    ## ... Procrustes: rmse 0.01516018  max resid 0.04197485 
    ## Run 133 stress 0.06124505 
    ## ... Procrustes: rmse 0.01513206  max resid 0.04190049 
    ## Run 134 stress 0.06124515 
    ## ... Procrustes: rmse 0.01520913  max resid 0.04210689 
    ## Run 135 stress 0.07411939 
    ## Run 136 stress 0.07166969 
    ## Run 137 stress 0.07411941 
    ## Run 138 stress 0.061245 
    ## ... Procrustes: rmse 0.01507768  max resid 0.04175574 
    ## Run 139 stress 0.07411949 
    ## Run 140 stress 0.06623458 
    ## Run 141 stress 0.06623458 
    ## Run 142 stress 0.07411938 
    ## Run 143 stress 0.2345379 
    ## Run 144 stress 0.07411944 
    ## Run 145 stress 0.07166985 
    ## Run 146 stress 0.06477865 
    ## Run 147 stress 0.06123361 
    ## ... Procrustes: rmse 3.620275e-05  max resid 8.500745e-05 
    ## ... Similar to previous best
    ## Run 148 stress 0.06124508 
    ## ... Procrustes: rmse 0.01515665  max resid 0.0419666 
    ## Run 149 stress 0.2324715 
    ## Run 150 stress 0.07166978 
    ## Run 151 stress 0.06477948 
    ## Run 152 stress 0.07411951 
    ## Run 153 stress 0.06477866 
    ## Run 154 stress 0.06477829 
    ## Run 155 stress 0.06477827 
    ## Run 156 stress 0.07411946 
    ## Run 157 stress 0.06623458 
    ## Run 158 stress 0.06623469 
    ## Run 159 stress 0.07411944 
    ## Run 160 stress 0.06123361 
    ## ... Procrustes: rmse 2.764685e-05  max resid 4.886814e-05 
    ## ... Similar to previous best
    ## Run 161 stress 0.07166966 
    ## Run 162 stress 0.07166985 
    ## Run 163 stress 0.07411942 
    ## Run 164 stress 0.06124487 
    ## ... Procrustes: rmse 0.01486878  max resid 0.04119608 
    ## Run 165 stress 0.07411944 
    ## Run 166 stress 0.07411938 
    ## Run 167 stress 0.07411946 
    ## Run 168 stress 0.06623459 
    ## Run 169 stress 0.06124496 
    ## ... Procrustes: rmse 0.01503964  max resid 0.04165291 
    ## Run 170 stress 0.07411947 
    ## Run 171 stress 0.06477862 
    ## Run 172 stress 0.07411941 
    ## Run 173 stress 0.2859295 
    ## Run 174 stress 0.07411942 
    ## Run 175 stress 0.06477864 
    ## Run 176 stress 0.07166974 
    ## Run 177 stress 0.07411951 
    ## Run 178 stress 0.0612336 
    ## ... Procrustes: rmse 1.357944e-05  max resid 1.674376e-05 
    ## ... Similar to previous best
    ## Run 179 stress 0.2361811 
    ## Run 180 stress 0.2394588 
    ## Run 181 stress 0.06124512 
    ## ... Procrustes: rmse 0.0151819  max resid 0.04203308 
    ## Run 182 stress 0.0612336 
    ## ... Procrustes: rmse 9.636823e-06  max resid 1.668597e-05 
    ## ... Similar to previous best
    ## Run 183 stress 0.07166977 
    ## Run 184 stress 0.06623461 
    ## Run 185 stress 0.06477932 
    ## Run 186 stress 0.06623463 
    ## Run 187 stress 0.06123361 
    ## ... Procrustes: rmse 2.901981e-05  max resid 3.679273e-05 
    ## ... Similar to previous best
    ## Run 188 stress 0.06124499 
    ## ... Procrustes: rmse 0.01447972  max resid 0.0401522 
    ## Run 189 stress 0.06623458 
    ## Run 190 stress 0.07411949 
    ## Run 191 stress 0.06477864 
    ## Run 192 stress 0.08226161 
    ## Run 193 stress 0.06623457 
    ## Run 194 stress 0.06623457 
    ## Run 195 stress 0.06124513 
    ## ... Procrustes: rmse 0.01519168  max resid 0.0420617 
    ## Run 196 stress 0.06477831 
    ## Run 197 stress 0.06477883 
    ## Run 198 stress 0.06477872 
    ## Run 199 stress 0.07166986 
    ## Run 200 stress 0.06477919 
    ## Run 201 stress 0.0647785 
    ## Run 202 stress 0.06623458 
    ## Run 203 stress 0.07166966 
    ## Run 204 stress 0.06623458 
    ## Run 205 stress 0.07411942 
    ## Run 206 stress 0.06623457 
    ## Run 207 stress 0.07411941 
    ## Run 208 stress 0.07166969 
    ## Run 209 stress 0.06623458 
    ## Run 210 stress 0.07166974 
    ## Run 211 stress 0.3275207 
    ## Run 212 stress 0.0612336 
    ## ... Procrustes: rmse 6.056613e-06  max resid 9.437179e-06 
    ## ... Similar to previous best
    ## Run 213 stress 0.06477896 
    ## Run 214 stress 0.06623457 
    ## Run 215 stress 0.06623463 
    ## Run 216 stress 0.3114814 
    ## Run 217 stress 0.0662346 
    ## Run 218 stress 0.06477843 
    ## Run 219 stress 0.07166969 
    ## Run 220 stress 0.07166973 
    ## Run 221 stress 0.061245 
    ## ... Procrustes: rmse 0.01508041  max resid 0.04176339 
    ## Run 222 stress 0.06623458 
    ## Run 223 stress 0.08226163 
    ## Run 224 stress 0.06124509 
    ## ... Procrustes: rmse 0.01516128  max resid 0.04197873 
    ## Run 225 stress 0.06477919 
    ## Run 226 stress 0.3407336 
    ## Run 227 stress 0.08226156 
    ## Run 228 stress 0.06623458 
    ## Run 229 stress 0.07166972 
    ## Run 230 stress 0.06124487 
    ## ... Procrustes: rmse 0.0146821  max resid 0.04069461 
    ## Run 231 stress 0.0612336 
    ## ... Procrustes: rmse 1.023798e-05  max resid 1.298367e-05 
    ## ... Similar to previous best
    ## Run 232 stress 0.06623463 
    ## Run 233 stress 0.06477905 
    ## Run 234 stress 0.06477827 
    ## Run 235 stress 0.06123361 
    ## ... Procrustes: rmse 3.582073e-05  max resid 4.944384e-05 
    ## ... Similar to previous best
    ## Run 236 stress 0.0662346 
    ## Run 237 stress 0.0741194 
    ## Run 238 stress 0.07411947 
    ## Run 239 stress 0.2548998 
    ## Run 240 stress 0.06477826 
    ## Run 241 stress 0.06477827 
    ## Run 242 stress 0.07411939 
    ## Run 243 stress 0.07411938 
    ## Run 244 stress 0.07166974 
    ## Run 245 stress 0.07411941 
    ## Run 246 stress 0.06623458 
    ## Run 247 stress 0.06124496 
    ## ... Procrustes: rmse 0.0150376  max resid 0.04164834 
    ## Run 248 stress 0.06124525 
    ## ... Procrustes: rmse 0.01428728  max resid 0.03963685 
    ## Run 249 stress 0.07411938 
    ## Run 250 stress 0.0741194 
    ## Run 251 stress 0.06477866 
    ## Run 252 stress 0.07411952 
    ## Run 253 stress 0.06623461 
    ## Run 254 stress 0.07411947 
    ## Run 255 stress 0.06477839 
    ## Run 256 stress 0.08226172 
    ## Run 257 stress 0.064779 
    ## Run 258 stress 0.06623458 
    ## Run 259 stress 0.07411945 
    ## Run 260 stress 0.07411951 
    ## Run 261 stress 0.0662346 
    ## Run 262 stress 0.06623457 
    ## Run 263 stress 0.06623463 
    ## Run 264 stress 0.06623467 
    ## Run 265 stress 0.07411952 
    ## Run 266 stress 0.07411951 
    ## Run 267 stress 0.06623457 
    ## Run 268 stress 0.06124512 
    ## ... Procrustes: rmse 0.01518717  max resid 0.04204858 
    ## Run 269 stress 0.07166966 
    ## Run 270 stress 0.07411945 
    ## Run 271 stress 0.07166965 
    ## Run 272 stress 0.0612336 
    ## ... Procrustes: rmse 1.108329e-05  max resid 1.785852e-05 
    ## ... Similar to previous best
    ## Run 273 stress 0.06477896 
    ## Run 274 stress 0.3192438 
    ## Run 275 stress 0.06124511 
    ## ... Procrustes: rmse 0.01517076  max resid 0.04200563 
    ## Run 276 stress 0.0822617 
    ## Run 277 stress 0.06477854 
    ## Run 278 stress 0.0612451 
    ## ... Procrustes: rmse 0.01517655  max resid 0.04202037 
    ## Run 279 stress 0.06123361 
    ## ... Procrustes: rmse 2.673782e-05  max resid 4.022377e-05 
    ## ... Similar to previous best
    ## Run 280 stress 0.06623461 
    ## Run 281 stress 0.07411938 
    ## Run 282 stress 0.07166965 
    ## Run 283 stress 0.07166965 
    ## Run 284 stress 0.07166965 
    ## Run 285 stress 0.0612336 
    ## ... Procrustes: rmse 2.855557e-06  max resid 4.91195e-06 
    ## ... Similar to previous best
    ## Run 286 stress 0.08226174 
    ## Run 287 stress 0.06623467 
    ## Run 288 stress 0.3270095 
    ## Run 289 stress 0.06477919 
    ## Run 290 stress 0.2345379 
    ## Run 291 stress 0.06124507 
    ## ... Procrustes: rmse 0.01514217  max resid 0.0419272 
    ## Run 292 stress 0.06124493 
    ## ... Procrustes: rmse 0.01498652  max resid 0.04151035 
    ## Run 293 stress 0.0662346 
    ## Run 294 stress 0.07166986 
    ## Run 295 stress 0.06124512 
    ## ... Procrustes: rmse 0.01518627  max resid 0.04204633 
    ## Run 296 stress 0.06123361 
    ## ... Procrustes: rmse 2.264174e-05  max resid 3.17006e-05 
    ## ... Similar to previous best
    ## Run 297 stress 0.06623458 
    ## Run 298 stress 0.07166973 
    ## Run 299 stress 0.07166974 
    ## Run 300 stress 0.06124504 
    ## ... Procrustes: rmse 0.01512201  max resid 0.04187393 
    ## Run 301 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 3.720058e-06  max resid 7.390619e-06 
    ## ... Similar to previous best
    ## Run 302 stress 0.0741194 
    ## Run 303 stress 0.06123361 
    ## ... Procrustes: rmse 3.876224e-05  max resid 7.130677e-05 
    ## ... Similar to previous best
    ## Run 304 stress 0.07411949 
    ## Run 305 stress 0.08226168 
    ## Run 306 stress 0.07166979 
    ## Run 307 stress 0.07411941 
    ## Run 308 stress 0.06124516 
    ## ... Procrustes: rmse 0.01433884  max resid 0.03977428 
    ## Run 309 stress 0.06623458 
    ## Run 310 stress 0.06623458 
    ## Run 311 stress 0.08226165 
    ## Run 312 stress 0.0741194 
    ## Run 313 stress 0.06623458 
    ## Run 314 stress 0.06477826 
    ## Run 315 stress 0.06623465 
    ## Run 316 stress 0.06623464 
    ## Run 317 stress 0.07411938 
    ## Run 318 stress 0.06623468 
    ## Run 319 stress 0.06623457 
    ## Run 320 stress 0.06623468 
    ## Run 321 stress 0.06124521 
    ## ... Procrustes: rmse 0.01525205  max resid 0.04222184 
    ## Run 322 stress 0.07166971 
    ## Run 323 stress 0.07166977 
    ## Run 324 stress 0.08226164 
    ## Run 325 stress 0.07166969 
    ## Run 326 stress 0.07166965 
    ## Run 327 stress 0.0716697 
    ## Run 328 stress 0.06123361 
    ## ... Procrustes: rmse 1.620356e-05  max resid 2.508088e-05 
    ## ... Similar to previous best
    ## Run 329 stress 0.0662346 
    ## Run 330 stress 0.06123363 
    ## ... Procrustes: rmse 8.30817e-05  max resid 0.0001011278 
    ## ... Similar to previous best
    ## Run 331 stress 0.08226165 
    ## Run 332 stress 0.06477936 
    ## Run 333 stress 0.07411942 
    ## Run 334 stress 0.06477844 
    ## Run 335 stress 0.2362366 
    ## Run 336 stress 0.07411956 
    ## Run 337 stress 0.06477928 
    ## Run 338 stress 0.3008523 
    ## Run 339 stress 0.0612336 
    ## ... Procrustes: rmse 4.088918e-06  max resid 6.625399e-06 
    ## ... Similar to previous best
    ## Run 340 stress 0.06124508 
    ## ... Procrustes: rmse 0.01514723  max resid 0.04194282 
    ## Run 341 stress 0.06623458 
    ## Run 342 stress 0.06477931 
    ## Run 343 stress 0.07411938 
    ## Run 344 stress 0.06623459 
    ## Run 345 stress 0.0612336 
    ## ... Procrustes: rmse 9.815401e-06  max resid 2.324479e-05 
    ## ... Similar to previous best
    ## Run 346 stress 0.3311511 
    ## Run 347 stress 0.06124503 
    ## ... Procrustes: rmse 0.01443847  max resid 0.04004089 
    ## Run 348 stress 0.06477838 
    ## Run 349 stress 0.06623458 
    ## Run 350 stress 0.07411941 
    ## Run 351 stress 0.3252888 
    ## Run 352 stress 0.07411944 
    ## Run 353 stress 0.06124522 
    ## ... Procrustes: rmse 0.01523624  max resid 0.04217703 
    ## Run 354 stress 0.0716697 
    ## Run 355 stress 0.06623463 
    ## Run 356 stress 0.07166977 
    ## Run 357 stress 0.07411946 
    ## Run 358 stress 0.0612336 
    ## ... Procrustes: rmse 9.911791e-06  max resid 1.616144e-05 
    ## ... Similar to previous best
    ## Run 359 stress 0.07411941 
    ## Run 360 stress 0.2518589 
    ## Run 361 stress 0.3389606 
    ## Run 362 stress 0.06124503 
    ## ... Procrustes: rmse 0.0151152  max resid 0.04185539 
    ## Run 363 stress 0.06124498 
    ## ... Procrustes: rmse 0.01505835  max resid 0.04170365 
    ## Run 364 stress 0.3257199 
    ## Run 365 stress 0.06477864 
    ## Run 366 stress 0.06124512 
    ## ... Procrustes: rmse 0.01518606  max resid 0.04204443 
    ## Run 367 stress 0.06623458 
    ## Run 368 stress 0.06477902 
    ## Run 369 stress 0.06623458 
    ## Run 370 stress 0.06623462 
    ## Run 371 stress 0.07411952 
    ## Run 372 stress 0.06124486 
    ## ... Procrustes: rmse 0.01476133  max resid 0.04090769 
    ## Run 373 stress 0.0741194 
    ## Run 374 stress 0.0716697 
    ## Run 375 stress 0.07166994 
    ## Run 376 stress 0.07411938 
    ## Run 377 stress 0.06477863 
    ## Run 378 stress 0.06623459 
    ## Run 379 stress 0.3407334 
    ## Run 380 stress 0.06477835 
    ## Run 381 stress 0.0741194 
    ## Run 382 stress 0.06478001 
    ## Run 383 stress 0.06477962 
    ## Run 384 stress 0.06623457 
    ## Run 385 stress 0.06124512 
    ## ... Procrustes: rmse 0.01518897  max resid 0.04205329 
    ## Run 386 stress 0.06623462 
    ## Run 387 stress 0.06623467 
    ## Run 388 stress 0.0647787 
    ## Run 389 stress 0.06623458 
    ## Run 390 stress 0.08226172 
    ## Run 391 stress 0.07411939 
    ## Run 392 stress 0.06623459 
    ## Run 393 stress 0.07166967 
    ## Run 394 stress 0.07166986 
    ## Run 395 stress 0.0612451 
    ## ... Procrustes: rmse 0.01516305  max resid 0.04198214 
    ## Run 396 stress 0.0741195 
    ## Run 397 stress 0.06623459 
    ## Run 398 stress 0.06477867 
    ## Run 399 stress 0.06124512 
    ## ... Procrustes: rmse 0.01518935  max resid 0.0420545 
    ## Run 400 stress 0.08226162 
    ## Run 401 stress 0.06124516 
    ## ... Procrustes: rmse 0.01433886  max resid 0.03977424 
    ## Run 402 stress 0.06623457 
    ## Run 403 stress 0.06623457 
    ## Run 404 stress 0.3156069 
    ## Run 405 stress 0.06623459 
    ## Run 406 stress 0.0741195 
    ## Run 407 stress 0.07411942 
    ## Run 408 stress 0.07411938 
    ## Run 409 stress 0.06477964 
    ## Run 410 stress 0.07411948 
    ## Run 411 stress 0.07166969 
    ## Run 412 stress 0.07411939 
    ## Run 413 stress 0.3407342 
    ## Run 414 stress 0.06124517 
    ## ... Procrustes: rmse 0.01522314  max resid 0.04214427 
    ## Run 415 stress 0.07411943 
    ## Run 416 stress 0.06477837 
    ## Run 417 stress 0.07411939 
    ## Run 418 stress 0.06477901 
    ## Run 419 stress 0.08226158 
    ## Run 420 stress 0.07411946 
    ## Run 421 stress 0.06477994 
    ## Run 422 stress 0.0741195 
    ## Run 423 stress 0.06124487 
    ## ... Procrustes: rmse 0.01486546  max resid 0.04118726 
    ## Run 424 stress 0.09623049 
    ## Run 425 stress 0.08226163 
    ## Run 426 stress 0.2345379 
    ## Run 427 stress 0.06124488 
    ## ... Procrustes: rmse 0.0149065  max resid 0.04129707 
    ## Run 428 stress 0.07411948 
    ## Run 429 stress 0.07166976 
    ## Run 430 stress 0.0822618 
    ## Run 431 stress 0.0612336 
    ## ... Procrustes: rmse 7.095431e-06  max resid 1.372765e-05 
    ## ... Similar to previous best
    ## Run 432 stress 0.07411939 
    ## Run 433 stress 0.07166975 
    ## Run 434 stress 0.07411944 
    ## Run 435 stress 0.07166986 
    ## Run 436 stress 0.07411943 
    ## Run 437 stress 0.06124516 
    ## ... Procrustes: rmse 0.01520221  max resid 0.04208644 
    ## Run 438 stress 0.06124487 
    ## ... Procrustes: rmse 0.0148273  max resid 0.04108263 
    ## Run 439 stress 0.07411938 
    ## Run 440 stress 0.07411938 
    ## Run 441 stress 0.06124517 
    ## ... Procrustes: rmse 0.01522284  max resid 0.04214456 
    ## Run 442 stress 0.06623461 
    ## Run 443 stress 0.07411951 
    ## Run 444 stress 0.06623464 
    ## Run 445 stress 0.0741194 
    ## Run 446 stress 0.06124526 
    ## ... Procrustes: rmse 0.01528336  max resid 0.0423058 
    ## Run 447 stress 0.3407337 
    ## Run 448 stress 0.06124512 
    ## ... Procrustes: rmse 0.01518156  max resid 0.04203403 
    ## Run 449 stress 0.07411942 
    ## Run 450 stress 0.06124488 
    ## ... Procrustes: rmse 0.01488971  max resid 0.04125179 
    ## Run 451 stress 0.07411946 
    ## Run 452 stress 0.06623461 
    ## Run 453 stress 0.06123362 
    ## ... Procrustes: rmse 5.687549e-05  max resid 9.821215e-05 
    ## ... Similar to previous best
    ## Run 454 stress 0.07166974 
    ## Run 455 stress 0.06123361 
    ## ... Procrustes: rmse 4.195969e-05  max resid 9.698413e-05 
    ## ... Similar to previous best
    ## Run 456 stress 0.07411941 
    ## Run 457 stress 0.06123362 
    ## ... Procrustes: rmse 6.232847e-05  max resid 7.892266e-05 
    ## ... Similar to previous best
    ## Run 458 stress 0.06623463 
    ## Run 459 stress 0.06623457 
    ## Run 460 stress 0.3114121 
    ## Run 461 stress 0.06123361 
    ## ... Procrustes: rmse 2.737034e-05  max resid 3.893202e-05 
    ## ... Similar to previous best
    ## Run 462 stress 0.06477889 
    ## Run 463 stress 0.0612336 
    ## ... Procrustes: rmse 4.367887e-06  max resid 7.500591e-06 
    ## ... Similar to previous best
    ## Run 464 stress 0.06477831 
    ## Run 465 stress 0.06623458 
    ## Run 466 stress 0.07411944 
    ## Run 467 stress 0.0612451 
    ## ... Procrustes: rmse 0.01516199  max resid 0.04197948 
    ## Run 468 stress 0.07166986 
    ## Run 469 stress 0.0741194 
    ## Run 470 stress 0.0612336 
    ## ... Procrustes: rmse 3.482698e-06  max resid 6.497366e-06 
    ## ... Similar to previous best
    ## Run 471 stress 0.06124511 
    ## ... Procrustes: rmse 0.01515857  max resid 0.04196947 
    ## Run 472 stress 0.06623457 
    ## Run 473 stress 0.06477947 
    ## Run 474 stress 0.08226168 
    ## Run 475 stress 0.0612336 
    ## ... Procrustes: rmse 1.349991e-05  max resid 1.870392e-05 
    ## ... Similar to previous best
    ## Run 476 stress 0.07411949 
    ## Run 477 stress 0.08226161 
    ## Run 478 stress 0.07411938 
    ## Run 479 stress 0.06124514 
    ## ... Procrustes: rmse 0.01520072  max resid 0.04208492 
    ## Run 480 stress 0.07166969 
    ## Run 481 stress 0.07411945 
    ## Run 482 stress 0.07166988 
    ## Run 483 stress 0.07411943 
    ## Run 484 stress 0.07411939 
    ## Run 485 stress 0.06623459 
    ## Run 486 stress 0.06477878 
    ## Run 487 stress 0.07411946 
    ## Run 488 stress 0.06477908 
    ## Run 489 stress 0.07166979 
    ## Run 490 stress 0.2381016 
    ## Run 491 stress 0.06623466 
    ## Run 492 stress 0.07411948 
    ## Run 493 stress 0.06477894 
    ## Run 494 stress 0.06477964 
    ## Run 495 stress 0.07411944 
    ## Run 496 stress 0.06623457 
    ## Run 497 stress 0.06124493 
    ## ... Procrustes: rmse 0.01498572  max resid 0.04150882 
    ## Run 498 stress 0.07166978 
    ## Run 499 stress 0.08226164 
    ## Run 500 stress 0.3116178 
    ## *** Best solution repeated 15 times

``` r
# Stratified lakes and ocean sites
SD_beta_env_SO_NMDS <- metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.59483e-05 
    ## Run 1 stress 9.488499e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 4.383661e-05  max resid 8.649917e-05 
    ## ... Similar to previous best
    ## Run 2 stress 9.947058e-05 
    ## ... Procrustes: rmse 9.541678e-05  max resid 0.0002271953 
    ## ... Similar to previous best
    ## Run 3 stress 9.629336e-05 
    ## ... Procrustes: rmse 7.568999e-05  max resid 0.0001001005 
    ## ... Similar to previous best
    ## Run 4 stress 9.658084e-05 
    ## ... Procrustes: rmse 0.0001651279  max resid 0.0003621121 
    ## ... Similar to previous best
    ## Run 5 stress 9.764563e-05 
    ## ... Procrustes: rmse 0.0002256791  max resid 0.0003858341 
    ## ... Similar to previous best
    ## Run 6 stress 0.2693706 
    ## Run 7 stress 9.246052e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001629857  max resid 0.0003721912 
    ## ... Similar to previous best
    ## Run 8 stress 9.715297e-05 
    ## ... Procrustes: rmse 0.0001129189  max resid 0.0001804974 
    ## ... Similar to previous best
    ## Run 9 stress 9.36039e-05 
    ## ... Procrustes: rmse 0.0001087168  max resid 0.0001806611 
    ## ... Similar to previous best
    ## Run 10 stress 8.834176e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 8.226231e-05  max resid 0.000163308 
    ## ... Similar to previous best
    ## Run 11 stress 9.633394e-05 
    ## ... Procrustes: rmse 0.0001540728  max resid 0.0002571639 
    ## ... Similar to previous best
    ## Run 12 stress 9.332205e-05 
    ## ... Procrustes: rmse 0.0001616308  max resid 0.0002366367 
    ## ... Similar to previous best
    ## Run 13 stress 9.912122e-05 
    ## ... Procrustes: rmse 0.0001944507  max resid 0.0003624985 
    ## ... Similar to previous best
    ## Run 14 stress 9.62954e-05 
    ## ... Procrustes: rmse 0.0002021587  max resid 0.0003667254 
    ## ... Similar to previous best
    ## Run 15 stress 8.927238e-05 
    ## ... Procrustes: rmse 0.00017025  max resid 0.0003654976 
    ## ... Similar to previous best
    ## Run 16 stress 8.635319e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002114553  max resid 0.0003110103 
    ## ... Similar to previous best
    ## Run 17 stress 9.566856e-05 
    ## ... Procrustes: rmse 0.0002208114  max resid 0.0003478102 
    ## ... Similar to previous best
    ## Run 18 stress 9.160745e-05 
    ## ... Procrustes: rmse 0.0002066585  max resid 0.0004229575 
    ## ... Similar to previous best
    ## Run 19 stress 9.798047e-05 
    ## ... Procrustes: rmse 0.0002247235  max resid 0.0003285442 
    ## ... Similar to previous best
    ## Run 20 stress 8.613716e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002328939  max resid 0.00042288 
    ## ... Similar to previous best
    ## Run 21 stress 9.566545e-05 
    ## ... Procrustes: rmse 0.0001683474  max resid 0.000341583 
    ## ... Similar to previous best
    ## Run 22 stress 9.836223e-05 
    ## ... Procrustes: rmse 2.532657e-05  max resid 4.522387e-05 
    ## ... Similar to previous best
    ## Run 23 stress 8.799395e-05 
    ## ... Procrustes: rmse 0.000220966  max resid 0.0003345683 
    ## ... Similar to previous best
    ## Run 24 stress 9.992686e-05 
    ## ... Procrustes: rmse 0.0002381672  max resid 0.0003861881 
    ## ... Similar to previous best
    ## Run 25 stress 9.943758e-05 
    ## ... Procrustes: rmse 0.0001751914  max resid 0.0003276125 
    ## ... Similar to previous best
    ## Run 26 stress 9.586844e-05 
    ## ... Procrustes: rmse 0.0002007249  max resid 0.0004127402 
    ## ... Similar to previous best
    ## Run 27 stress 9.973932e-05 
    ## ... Procrustes: rmse 0.0001722779  max resid 0.0003644504 
    ## ... Similar to previous best
    ## Run 28 stress 9.623545e-05 
    ## ... Procrustes: rmse 0.0001985797  max resid 0.0004317061 
    ## ... Similar to previous best
    ## Run 29 stress 9.463544e-05 
    ## ... Procrustes: rmse 0.0001945797  max resid 0.0004114994 
    ## ... Similar to previous best
    ## Run 30 stress 9.358027e-05 
    ## ... Procrustes: rmse 0.0002048842  max resid 0.0003858862 
    ## ... Similar to previous best
    ## Run 31 stress 8.826064e-05 
    ## ... Procrustes: rmse 0.0001732605  max resid 0.0003477318 
    ## ... Similar to previous best
    ## Run 32 stress 9.380794e-05 
    ## ... Procrustes: rmse 0.000164013  max resid 0.0003347423 
    ## ... Similar to previous best
    ## Run 33 stress 9.071053e-05 
    ## ... Procrustes: rmse 8.200697e-05  max resid 0.0001952915 
    ## ... Similar to previous best
    ## Run 34 stress 9.625039e-05 
    ## ... Procrustes: rmse 0.0001911326  max resid 0.0004661209 
    ## ... Similar to previous best
    ## Run 35 stress 7.031983e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001262805  max resid 0.0002651503 
    ## ... Similar to previous best
    ## Run 36 stress 9.886018e-05 
    ## ... Procrustes: rmse 0.0002153287  max resid 0.0004082079 
    ## ... Similar to previous best
    ## Run 37 stress 9.567212e-05 
    ## ... Procrustes: rmse 0.0002232786  max resid 0.000312571 
    ## ... Similar to previous best
    ## Run 38 stress 9.947841e-05 
    ## ... Procrustes: rmse 0.0001065972  max resid 0.0001588088 
    ## ... Similar to previous best
    ## Run 39 stress 9.569433e-05 
    ## ... Procrustes: rmse 0.0002070889  max resid 0.0003879838 
    ## ... Similar to previous best
    ## Run 40 stress 9.467065e-05 
    ## ... Procrustes: rmse 0.0002183815  max resid 0.0003088697 
    ## ... Similar to previous best
    ## Run 41 stress 9.762293e-05 
    ## ... Procrustes: rmse 0.0001984605  max resid 0.0003660933 
    ## ... Similar to previous best
    ## Run 42 stress 9.406556e-05 
    ## ... Procrustes: rmse 0.0002002221  max resid 0.0003652513 
    ## ... Similar to previous best
    ## Run 43 stress 9.94107e-05 
    ## ... Procrustes: rmse 0.0001751072  max resid 0.0002887309 
    ## ... Similar to previous best
    ## Run 44 stress 9.999236e-05 
    ## ... Procrustes: rmse 0.0002191078  max resid 0.0004097768 
    ## ... Similar to previous best
    ## Run 45 stress 9.422517e-05 
    ## ... Procrustes: rmse 0.0002081823  max resid 0.0003893413 
    ## ... Similar to previous best
    ## Run 46 stress 9.351667e-05 
    ## ... Procrustes: rmse 0.0001610622  max resid 0.0002748263 
    ## ... Similar to previous best
    ## Run 47 stress 9.848895e-05 
    ## ... Procrustes: rmse 0.0001842348  max resid 0.0003622965 
    ## ... Similar to previous best
    ## Run 48 stress 9.316986e-05 
    ## ... Procrustes: rmse 0.0001175745  max resid 0.0001944806 
    ## ... Similar to previous best
    ## Run 49 stress 9.557409e-05 
    ## ... Procrustes: rmse 0.0001789364  max resid 0.0003562862 
    ## ... Similar to previous best
    ## Run 50 stress 9.042818e-05 
    ## ... Procrustes: rmse 0.0001738468  max resid 0.0003311991 
    ## ... Similar to previous best
    ## Run 51 stress 9.704997e-05 
    ## ... Procrustes: rmse 0.0001749166  max resid 0.0003497618 
    ## ... Similar to previous best
    ## Run 52 stress 9.983435e-05 
    ## ... Procrustes: rmse 0.0002151133  max resid 0.0004103087 
    ## ... Similar to previous best
    ## Run 53 stress 8.673282e-05 
    ## ... Procrustes: rmse 0.0001621916  max resid 0.0003380364 
    ## ... Similar to previous best
    ## Run 54 stress 9.288801e-05 
    ## ... Procrustes: rmse 0.0001845324  max resid 0.0003491943 
    ## ... Similar to previous best
    ## Run 55 stress 9.790196e-05 
    ## ... Procrustes: rmse 0.0001299274  max resid 0.0002344983 
    ## ... Similar to previous best
    ## Run 56 stress 9.950677e-05 
    ## ... Procrustes: rmse 0.0002079806  max resid 0.0003893054 
    ## ... Similar to previous best
    ## Run 57 stress 8.806025e-05 
    ## ... Procrustes: rmse 0.0001846816  max resid 0.0003453745 
    ## ... Similar to previous best
    ## Run 58 stress 9.948631e-05 
    ## ... Procrustes: rmse 9.5104e-05  max resid 0.0001307971 
    ## ... Similar to previous best
    ## Run 59 stress 9.361566e-05 
    ## ... Procrustes: rmse 0.0001549189  max resid 0.0002656168 
    ## ... Similar to previous best
    ## Run 60 stress 9.279654e-05 
    ## ... Procrustes: rmse 0.0001718157  max resid 0.0003552409 
    ## ... Similar to previous best
    ## Run 61 stress 8.960391e-05 
    ## ... Procrustes: rmse 0.0001237328  max resid 0.0002349809 
    ## ... Similar to previous best
    ## Run 62 stress 9.937511e-05 
    ## ... Procrustes: rmse 0.0002160568  max resid 0.0004063021 
    ## ... Similar to previous best
    ## Run 63 stress 9.916689e-05 
    ## ... Procrustes: rmse 0.0001742446  max resid 0.0002887948 
    ## ... Similar to previous best
    ## Run 64 stress 9.976547e-05 
    ## ... Procrustes: rmse 0.0001455547  max resid 0.0002691354 
    ## ... Similar to previous best
    ## Run 65 stress 9.737231e-05 
    ## ... Procrustes: rmse 0.0002146462  max resid 0.0003940839 
    ## ... Similar to previous best
    ## Run 66 stress 9.990789e-05 
    ## ... Procrustes: rmse 0.0001953323  max resid 0.0003261697 
    ## ... Similar to previous best
    ## Run 67 stress 9.594324e-05 
    ## ... Procrustes: rmse 0.0001330355  max resid 0.0002640186 
    ## ... Similar to previous best
    ## Run 68 stress 9.445498e-05 
    ## ... Procrustes: rmse 0.0001724328  max resid 0.0003311036 
    ## ... Similar to previous best
    ## Run 69 stress 9.919836e-05 
    ## ... Procrustes: rmse 0.000134582  max resid 0.0002134543 
    ## ... Similar to previous best
    ## Run 70 stress 9.897937e-05 
    ## ... Procrustes: rmse 0.0001731106  max resid 0.0002871645 
    ## ... Similar to previous best
    ## Run 71 stress 9.393533e-05 
    ## ... Procrustes: rmse 0.0001875389  max resid 0.0003350131 
    ## ... Similar to previous best
    ## Run 72 stress 9.740618e-05 
    ## ... Procrustes: rmse 0.0002062307  max resid 0.000392615 
    ## ... Similar to previous best
    ## Run 73 stress 9.83397e-05 
    ## ... Procrustes: rmse 0.0001778548  max resid 0.0003033885 
    ## ... Similar to previous best
    ## Run 74 stress 9.081202e-05 
    ## ... Procrustes: rmse 0.0002008772  max resid 0.0003826036 
    ## ... Similar to previous best
    ## Run 75 stress 9.697319e-05 
    ## ... Procrustes: rmse 0.0001842566  max resid 0.0003192516 
    ## ... Similar to previous best
    ## Run 76 stress 8.858614e-05 
    ## ... Procrustes: rmse 0.0001939302  max resid 0.0003760784 
    ## ... Similar to previous best
    ## Run 77 stress 8.382709e-05 
    ## ... Procrustes: rmse 0.000178617  max resid 0.0003154145 
    ## ... Similar to previous best
    ## Run 78 stress 9.533831e-05 
    ## ... Procrustes: rmse 0.0001614678  max resid 0.0002636017 
    ## ... Similar to previous best
    ## Run 79 stress 9.559955e-05 
    ## ... Procrustes: rmse 0.0001826934  max resid 0.0003024427 
    ## ... Similar to previous best
    ## Run 80 stress 9.7357e-05 
    ## ... Procrustes: rmse 0.0001731092  max resid 0.0002725981 
    ## ... Similar to previous best
    ## Run 81 stress 9.498442e-05 
    ## ... Procrustes: rmse 0.0002178554  max resid 0.0003060967 
    ## ... Similar to previous best
    ## Run 82 stress 9.668358e-05 
    ## ... Procrustes: rmse 0.0002027636  max resid 0.0003697579 
    ## ... Similar to previous best
    ## Run 83 stress 9.819363e-05 
    ## ... Procrustes: rmse 0.0002002799  max resid 0.0003877973 
    ## ... Similar to previous best
    ## Run 84 stress 9.614252e-05 
    ## ... Procrustes: rmse 0.0002019137  max resid 0.000367664 
    ## ... Similar to previous best
    ## Run 85 stress 9.583674e-05 
    ## ... Procrustes: rmse 0.0001628522  max resid 0.0002744181 
    ## ... Similar to previous best
    ## Run 86 stress 9.506525e-05 
    ## ... Procrustes: rmse 0.0001627726  max resid 0.0002758641 
    ## ... Similar to previous best
    ## Run 87 stress 9.454926e-05 
    ## ... Procrustes: rmse 0.00011046  max resid 0.0001944478 
    ## ... Similar to previous best
    ## Run 88 stress 9.838008e-05 
    ## ... Procrustes: rmse 0.0001844143  max resid 0.0002606379 
    ## ... Similar to previous best
    ## Run 89 stress 8.93065e-05 
    ## ... Procrustes: rmse 9.912254e-05  max resid 0.0002295921 
    ## ... Similar to previous best
    ## Run 90 stress 9.74341e-05 
    ## ... Procrustes: rmse 0.0002178182  max resid 0.0003406568 
    ## ... Similar to previous best
    ## Run 91 stress 9.846621e-05 
    ## ... Procrustes: rmse 0.0002148739  max resid 0.000405367 
    ## ... Similar to previous best
    ## Run 92 stress 9.675013e-05 
    ## ... Procrustes: rmse 0.0001192772  max resid 0.0001979202 
    ## ... Similar to previous best
    ## Run 93 stress 9.937302e-05 
    ## ... Procrustes: rmse 0.0002130239  max resid 0.000414675 
    ## ... Similar to previous best
    ## Run 94 stress 9.386116e-05 
    ## ... Procrustes: rmse 0.0001367403  max resid 0.0002747671 
    ## ... Similar to previous best
    ## Run 95 stress 9.808084e-05 
    ## ... Procrustes: rmse 0.0001992833  max resid 0.0003783933 
    ## ... Similar to previous best
    ## Run 96 stress 9.706163e-05 
    ## ... Procrustes: rmse 0.0002118015  max resid 0.0004002716 
    ## ... Similar to previous best
    ## Run 97 stress 9.94746e-05 
    ## ... Procrustes: rmse 0.0001351875  max resid 0.0002138811 
    ## ... Similar to previous best
    ## Run 98 stress 9.382219e-05 
    ## ... Procrustes: rmse 0.00018703  max resid 0.0003593094 
    ## ... Similar to previous best
    ## Run 99 stress 9.880717e-05 
    ## ... Procrustes: rmse 0.0001733692  max resid 0.0002872877 
    ## ... Similar to previous best
    ## Run 100 stress 9.333509e-05 
    ## ... Procrustes: rmse 0.0001873299  max resid 0.000355343 
    ## ... Similar to previous best
    ## Run 101 stress 9.247786e-05 
    ## ... Procrustes: rmse 0.0001958638  max resid 0.00038544 
    ## ... Similar to previous best
    ## Run 102 stress 9.942476e-05 
    ## ... Procrustes: rmse 0.0002218832  max resid 0.0003483524 
    ## ... Similar to previous best
    ## Run 103 stress 9.633711e-05 
    ## ... Procrustes: rmse 0.0002176005  max resid 0.0003032151 
    ## ... Similar to previous best
    ## Run 104 stress 9.222963e-05 
    ## ... Procrustes: rmse 0.0002000271  max resid 0.0003827686 
    ## ... Similar to previous best
    ## Run 105 stress 9.229943e-05 
    ## ... Procrustes: rmse 0.000120345  max resid 0.0002005538 
    ## ... Similar to previous best
    ## Run 106 stress 9.649241e-05 
    ## ... Procrustes: rmse 0.0002118707  max resid 0.0003986147 
    ## ... Similar to previous best
    ## Run 107 stress 9.591656e-05 
    ## ... Procrustes: rmse 0.0001642686  max resid 0.0003432355 
    ## ... Similar to previous best
    ## Run 108 stress 9.043076e-05 
    ## ... Procrustes: rmse 0.0001404683  max resid 0.000311727 
    ## ... Similar to previous best
    ## Run 109 stress 9.337721e-05 
    ## ... Procrustes: rmse 0.0001747451  max resid 0.0003382742 
    ## ... Similar to previous best
    ## Run 110 stress 9.737691e-05 
    ## ... Procrustes: rmse 0.0002127156  max resid 0.0004047524 
    ## ... Similar to previous best
    ## Run 111 stress 9.892567e-05 
    ## ... Procrustes: rmse 0.0002137058  max resid 0.0004017099 
    ## ... Similar to previous best
    ## Run 112 stress 9.969404e-05 
    ## ... Procrustes: rmse 0.0002208425  max resid 0.0004036729 
    ## ... Similar to previous best
    ## Run 113 stress 9.762985e-05 
    ## ... Procrustes: rmse 0.0001878443  max resid 0.0003464863 
    ## ... Similar to previous best
    ## Run 114 stress 9.830817e-05 
    ## ... Procrustes: rmse 0.0002353306  max resid 0.0004168767 
    ## ... Similar to previous best
    ## Run 115 stress 8.490237e-05 
    ## ... Procrustes: rmse 0.0001425301  max resid 0.0002854242 
    ## ... Similar to previous best
    ## Run 116 stress 8.781482e-05 
    ## ... Procrustes: rmse 0.0001219739  max resid 0.0002565479 
    ## ... Similar to previous best
    ## Run 117 stress 9.323223e-05 
    ## ... Procrustes: rmse 0.0002064027  max resid 0.0003935134 
    ## ... Similar to previous best
    ## Run 118 stress 9.112183e-05 
    ## ... Procrustes: rmse 0.000112748  max resid 0.0001935066 
    ## ... Similar to previous best
    ## Run 119 stress 9.384495e-05 
    ## ... Procrustes: rmse 0.0001218772  max resid 0.0002099551 
    ## ... Similar to previous best
    ## Run 120 stress 9.576999e-05 
    ## ... Procrustes: rmse 0.0001977736  max resid 0.0003697606 
    ## ... Similar to previous best
    ## Run 121 stress 9.500419e-05 
    ## ... Procrustes: rmse 0.0001651197  max resid 0.0002749935 
    ## ... Similar to previous best
    ## Run 122 stress 7.664898e-05 
    ## ... Procrustes: rmse 0.0001184457  max resid 0.0002270956 
    ## ... Similar to previous best
    ## Run 123 stress 9.496954e-05 
    ## ... Procrustes: rmse 0.0001243648  max resid 0.0002025825 
    ## ... Similar to previous best
    ## Run 124 stress 9.967042e-05 
    ## ... Procrustes: rmse 0.0002093107  max resid 0.0003770136 
    ## ... Similar to previous best
    ## Run 125 stress 9.711313e-05 
    ## ... Procrustes: rmse 0.0002123547  max resid 0.0003995034 
    ## ... Similar to previous best
    ## Run 126 stress 9.899342e-05 
    ## ... Procrustes: rmse 0.0001739891  max resid 0.0002866252 
    ## ... Similar to previous best
    ## Run 127 stress 8.335312e-05 
    ## ... Procrustes: rmse 0.0001179693  max resid 0.0002600548 
    ## ... Similar to previous best
    ## Run 128 stress 9.435193e-05 
    ## ... Procrustes: rmse 0.0001357539  max resid 0.000272489 
    ## ... Similar to previous best
    ## Run 129 stress 9.747886e-05 
    ## ... Procrustes: rmse 0.0001291676  max resid 0.0002079292 
    ## ... Similar to previous best
    ## Run 130 stress 8.561055e-05 
    ## ... Procrustes: rmse 0.0001898796  max resid 0.0003655454 
    ## ... Similar to previous best
    ## Run 131 stress 9.164227e-05 
    ## ... Procrustes: rmse 0.0002084488  max resid 0.0003914351 
    ## ... Similar to previous best
    ## Run 132 stress 9.914541e-05 
    ## ... Procrustes: rmse 0.000182696  max resid 0.0003643134 
    ## ... Similar to previous best
    ## Run 133 stress 9.75783e-05 
    ## ... Procrustes: rmse 0.0002208931  max resid 0.0003093844 
    ## ... Similar to previous best
    ## Run 134 stress 9.853337e-05 
    ## ... Procrustes: rmse 0.0001917725  max resid 0.0003213863 
    ## ... Similar to previous best
    ## Run 135 stress 9.565398e-05 
    ## ... Procrustes: rmse 0.0001412995  max resid 0.0002785092 
    ## ... Similar to previous best
    ## Run 136 stress 9.591548e-05 
    ## ... Procrustes: rmse 0.000185795  max resid 0.0003302009 
    ## ... Similar to previous best
    ## Run 137 stress 8.188412e-05 
    ## ... Procrustes: rmse 0.0001817183  max resid 0.0003546784 
    ## ... Similar to previous best
    ## Run 138 stress 9.678449e-05 
    ## ... Procrustes: rmse 0.0001284493  max resid 0.0002086785 
    ## ... Similar to previous best
    ## Run 139 stress 9.824285e-05 
    ## ... Procrustes: rmse 0.0001665122  max resid 0.0002800701 
    ## ... Similar to previous best
    ## Run 140 stress 9.745248e-05 
    ## ... Procrustes: rmse 0.0001528794  max resid 0.0002634662 
    ## ... Similar to previous best
    ## Run 141 stress 8.75463e-05 
    ## ... Procrustes: rmse 0.0001683824  max resid 0.0003265867 
    ## ... Similar to previous best
    ## Run 142 stress 8.563023e-05 
    ## ... Procrustes: rmse 0.0001793723  max resid 0.0003549343 
    ## ... Similar to previous best
    ## Run 143 stress 9.962565e-05 
    ## ... Procrustes: rmse 0.0002164325  max resid 0.0004067235 
    ## ... Similar to previous best
    ## Run 144 stress 8.713137e-05 
    ## ... Procrustes: rmse 0.0001633113  max resid 0.0002360811 
    ## ... Similar to previous best
    ## Run 145 stress 9.673616e-05 
    ## ... Procrustes: rmse 0.0001857341  max resid 0.0003163262 
    ## ... Similar to previous best
    ## Run 146 stress 9.256307e-05 
    ## ... Procrustes: rmse 0.000163773  max resid 0.0002578751 
    ## ... Similar to previous best
    ## Run 147 stress 9.265389e-05 
    ## ... Procrustes: rmse 0.0001214358  max resid 0.0002010955 
    ## ... Similar to previous best
    ## Run 148 stress 9.949793e-05 
    ## ... Procrustes: rmse 9.047631e-05  max resid 0.0001170243 
    ## ... Similar to previous best
    ## Run 149 stress 9.784592e-05 
    ## ... Procrustes: rmse 0.0001682911  max resid 0.0002860843 
    ## ... Similar to previous best
    ## Run 150 stress 9.15259e-05 
    ## ... Procrustes: rmse 0.0002184013  max resid 0.0003890768 
    ## ... Similar to previous best
    ## Run 151 stress 9.298684e-05 
    ## ... Procrustes: rmse 7.784365e-05  max resid 0.0001199768 
    ## ... Similar to previous best
    ## Run 152 stress 8.989773e-05 
    ## ... Procrustes: rmse 0.0001866568  max resid 0.0003541573 
    ## ... Similar to previous best
    ## Run 153 stress 9.534413e-05 
    ## ... Procrustes: rmse 0.00016679  max resid 0.0002584203 
    ## ... Similar to previous best
    ## Run 154 stress 9.757175e-05 
    ## ... Procrustes: rmse 0.0001116253  max resid 0.0002267284 
    ## ... Similar to previous best
    ## Run 155 stress 9.905754e-05 
    ## ... Procrustes: rmse 0.0001332746  max resid 0.000211947 
    ## ... Similar to previous best
    ## Run 156 stress 8.811124e-05 
    ## ... Procrustes: rmse 0.0001488369  max resid 0.0002627579 
    ## ... Similar to previous best
    ## Run 157 stress 9.588546e-05 
    ## ... Procrustes: rmse 0.0001240588  max resid 0.0002039154 
    ## ... Similar to previous best
    ## Run 158 stress 9.896748e-05 
    ## ... Procrustes: rmse 0.0001338192  max resid 0.0002145742 
    ## ... Similar to previous best
    ## Run 159 stress 8.912177e-05 
    ## ... Procrustes: rmse 0.0001121651  max resid 0.0001922076 
    ## ... Similar to previous best
    ## Run 160 stress 9.96542e-05 
    ## ... Procrustes: rmse 0.00021077  max resid 0.000388511 
    ## ... Similar to previous best
    ## Run 161 stress 9.990475e-05 
    ## ... Procrustes: rmse 0.0002126107  max resid 0.0004144363 
    ## ... Similar to previous best
    ## Run 162 stress 9.348451e-05 
    ## ... Procrustes: rmse 0.0001633817  max resid 0.000278915 
    ## ... Similar to previous best
    ## Run 163 stress 9.009278e-05 
    ## ... Procrustes: rmse 0.0001113669  max resid 0.0001709155 
    ## ... Similar to previous best
    ## Run 164 stress 9.260401e-05 
    ## ... Procrustes: rmse 0.0001790151  max resid 0.0003012624 
    ## ... Similar to previous best
    ## Run 165 stress 9.51924e-05 
    ## ... Procrustes: rmse 0.0002045201  max resid 0.0003719862 
    ## ... Similar to previous best
    ## Run 166 stress 8.918381e-05 
    ## ... Procrustes: rmse 0.0001892829  max resid 0.0003490351 
    ## ... Similar to previous best
    ## Run 167 stress 9.523799e-05 
    ## ... Procrustes: rmse 0.0002092355  max resid 0.0003929672 
    ## ... Similar to previous best
    ## Run 168 stress 9.69711e-05 
    ## ... Procrustes: rmse 0.0002027075  max resid 0.0003808849 
    ## ... Similar to previous best
    ## Run 169 stress 8.205271e-05 
    ## ... Procrustes: rmse 9.036478e-05  max resid 0.0001726947 
    ## ... Similar to previous best
    ## Run 170 stress 9.639714e-05 
    ## ... Procrustes: rmse 0.0002088819  max resid 0.0003384634 
    ## ... Similar to previous best
    ## Run 171 stress 9.887834e-05 
    ## ... Procrustes: rmse 0.0002261021  max resid 0.0003535292 
    ## ... Similar to previous best
    ## Run 172 stress 9.293013e-05 
    ## ... Procrustes: rmse 0.0001899121  max resid 0.0003681051 
    ## ... Similar to previous best
    ## Run 173 stress 9.444966e-05 
    ## ... Procrustes: rmse 0.0001717843  max resid 0.0003486684 
    ## ... Similar to previous best
    ## Run 174 stress 9.395551e-05 
    ## ... Procrustes: rmse 0.0002047968  max resid 0.0003904177 
    ## ... Similar to previous best
    ## Run 175 stress 9.752069e-05 
    ## ... Procrustes: rmse 0.0002052401  max resid 0.0003796681 
    ## ... Similar to previous best
    ## Run 176 stress 9.183324e-05 
    ## ... Procrustes: rmse 0.0001713366  max resid 0.0003507943 
    ## ... Similar to previous best
    ## Run 177 stress 9.189126e-05 
    ## ... Procrustes: rmse 8.218522e-05  max resid 0.000125168 
    ## ... Similar to previous best
    ## Run 178 stress 9.852533e-05 
    ## ... Procrustes: rmse 0.0001976048  max resid 0.0003783611 
    ## ... Similar to previous best
    ## Run 179 stress 9.658868e-05 
    ## ... Procrustes: rmse 0.0001148954  max resid 0.000249209 
    ## ... Similar to previous best
    ## Run 180 stress 9.997575e-05 
    ## ... Procrustes: rmse 0.0001790249  max resid 0.0003516607 
    ## ... Similar to previous best
    ## Run 181 stress 9.956778e-05 
    ## ... Procrustes: rmse 0.0002084161  max resid 0.0003911075 
    ## ... Similar to previous best
    ## Run 182 stress 9.810619e-05 
    ## ... Procrustes: rmse 0.0002070838  max resid 0.000384347 
    ## ... Similar to previous best
    ## Run 183 stress 9.568234e-05 
    ## ... Procrustes: rmse 0.0002013909  max resid 0.000367691 
    ## ... Similar to previous best
    ## Run 184 stress 8.966778e-05 
    ## ... Procrustes: rmse 7.509594e-05  max resid 0.0001215624 
    ## ... Similar to previous best
    ## Run 185 stress 9.827709e-05 
    ## ... Procrustes: rmse 0.0002321378  max resid 0.0004076111 
    ## ... Similar to previous best
    ## Run 186 stress 9.880101e-05 
    ## ... Procrustes: rmse 0.0001730346  max resid 0.0002870184 
    ## ... Similar to previous best
    ## Run 187 stress 9.50811e-05 
    ## ... Procrustes: rmse 0.0002072105  max resid 0.0003922265 
    ## ... Similar to previous best
    ## Run 188 stress 9.624241e-05 
    ## ... Procrustes: rmse 7.07893e-05  max resid 0.0001434708 
    ## ... Similar to previous best
    ## Run 189 stress 9.360498e-05 
    ## ... Procrustes: rmse 0.000198208  max resid 0.0003891429 
    ## ... Similar to previous best
    ## Run 190 stress 9.671869e-05 
    ## ... Procrustes: rmse 0.0001829617  max resid 0.0003628433 
    ## ... Similar to previous best
    ## Run 191 stress 9.918203e-05 
    ## ... Procrustes: rmse 0.0002298526  max resid 0.0003177536 
    ## ... Similar to previous best
    ## Run 192 stress 9.761182e-05 
    ## ... Procrustes: rmse 0.0002176866  max resid 0.0003439053 
    ## ... Similar to previous best
    ## Run 193 stress 8.89841e-05 
    ## ... Procrustes: rmse 0.0001553646  max resid 0.000269516 
    ## ... Similar to previous best
    ## Run 194 stress 9.02278e-05 
    ## ... Procrustes: rmse 0.0001951465  max resid 0.0003916414 
    ## ... Similar to previous best
    ## Run 195 stress 6.209174e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000102173  max resid 0.0002222895 
    ## ... Similar to previous best
    ## Run 196 stress 9.333863e-05 
    ## ... Procrustes: rmse 0.0001664867  max resid 0.0002428166 
    ## ... Similar to previous best
    ## Run 197 stress 9.90853e-05 
    ## ... Procrustes: rmse 0.0001924498  max resid 0.0002739374 
    ## ... Similar to previous best
    ## Run 198 stress 9.040519e-05 
    ## ... Procrustes: rmse 0.0001180571  max resid 0.0001920411 
    ## ... Similar to previous best
    ## Run 199 stress 9.80151e-05 
    ## ... Procrustes: rmse 0.0001920565  max resid 0.000279861 
    ## ... Similar to previous best
    ## Run 200 stress 9.487424e-05 
    ## ... Procrustes: rmse 0.0001498165  max resid 0.0002138783 
    ## ... Similar to previous best
    ## Run 201 stress 9.637099e-05 
    ## ... Procrustes: rmse 0.0001512378  max resid 0.0002097871 
    ## ... Similar to previous best
    ## Run 202 stress 9.370028e-05 
    ## ... Procrustes: rmse 0.0001683077  max resid 0.0002672617 
    ## ... Similar to previous best
    ## Run 203 stress 8.29006e-05 
    ## ... Procrustes: rmse 0.0001248669  max resid 0.000207438 
    ## ... Similar to previous best
    ## Run 204 stress 9.045757e-05 
    ## ... Procrustes: rmse 0.00017718  max resid 0.000249968 
    ## ... Similar to previous best
    ## Run 205 stress 9.948205e-05 
    ## ... Procrustes: rmse 0.0001256419  max resid 0.0001861852 
    ## ... Similar to previous best
    ## Run 206 stress 9.579645e-05 
    ## ... Procrustes: rmse 0.0001580882  max resid 0.0002501961 
    ## ... Similar to previous best
    ## Run 207 stress 9.362527e-05 
    ## ... Procrustes: rmse 0.0001700616  max resid 0.0002445165 
    ## ... Similar to previous best
    ## Run 208 stress 9.995815e-05 
    ## ... Procrustes: rmse 0.00014932  max resid 0.0002221031 
    ## ... Similar to previous best
    ## Run 209 stress 9.740085e-05 
    ## ... Procrustes: rmse 0.0001701411  max resid 0.000257695 
    ## ... Similar to previous best
    ## Run 210 stress 9.494169e-05 
    ## ... Procrustes: rmse 0.00012853  max resid 0.0001746174 
    ## ... Similar to previous best
    ## Run 211 stress 9.845604e-05 
    ## ... Procrustes: rmse 0.0001925748  max resid 0.000262587 
    ## ... Similar to previous best
    ## Run 212 stress 8.537268e-05 
    ## ... Procrustes: rmse 0.0001232896  max resid 0.0002150391 
    ## ... Similar to previous best
    ## Run 213 stress 8.801137e-05 
    ## ... Procrustes: rmse 0.0001072602  max resid 0.0001716196 
    ## ... Similar to previous best
    ## Run 214 stress 9.364616e-05 
    ## ... Procrustes: rmse 0.000111853  max resid 0.0002330158 
    ## ... Similar to previous best
    ## Run 215 stress 9.684476e-05 
    ## ... Procrustes: rmse 0.0001661562  max resid 0.0002580504 
    ## ... Similar to previous best
    ## Run 216 stress 9.601287e-05 
    ## ... Procrustes: rmse 0.0001850325  max resid 0.0002681436 
    ## ... Similar to previous best
    ## Run 217 stress 9.610282e-05 
    ## ... Procrustes: rmse 0.0001769623  max resid 0.0002583502 
    ## ... Similar to previous best
    ## Run 218 stress 9.358799e-05 
    ## ... Procrustes: rmse 0.000141855  max resid 0.0001982323 
    ## ... Similar to previous best
    ## Run 219 stress 6.79201e-05 
    ## ... Procrustes: rmse 6.75114e-05  max resid 0.0001015529 
    ## ... Similar to previous best
    ## Run 220 stress 9.682519e-05 
    ## ... Procrustes: rmse 0.0001885597  max resid 0.000277559 
    ## ... Similar to previous best
    ## Run 221 stress 9.739609e-05 
    ## ... Procrustes: rmse 0.0001887059  max resid 0.0002758465 
    ## ... Similar to previous best
    ## Run 222 stress 9.914225e-05 
    ## ... Procrustes: rmse 0.0001535035  max resid 0.000235261 
    ## ... Similar to previous best
    ## Run 223 stress 0.2559352 
    ## Run 224 stress 9.840341e-05 
    ## ... Procrustes: rmse 0.0001397949  max resid 0.0002152112 
    ## ... Similar to previous best
    ## Run 225 stress 8.988902e-05 
    ## ... Procrustes: rmse 0.0001403682  max resid 0.0002477651 
    ## ... Similar to previous best
    ## Run 226 stress 9.301916e-05 
    ## ... Procrustes: rmse 0.0001778192  max resid 0.0002626132 
    ## ... Similar to previous best
    ## Run 227 stress 9.880545e-05 
    ## ... Procrustes: rmse 0.0001675906  max resid 0.0002379329 
    ## ... Similar to previous best
    ## Run 228 stress 9.867643e-05 
    ## ... Procrustes: rmse 0.0001594771  max resid 0.0002403955 
    ## ... Similar to previous best
    ## Run 229 stress 9.373478e-05 
    ## ... Procrustes: rmse 0.0001703529  max resid 0.0002583401 
    ## ... Similar to previous best
    ## Run 230 stress 9.843183e-05 
    ## ... Procrustes: rmse 0.0001944179  max resid 0.0002640151 
    ## ... Similar to previous best
    ## Run 231 stress 9.666449e-05 
    ## ... Procrustes: rmse 0.0001528095  max resid 0.0002360101 
    ## ... Similar to previous best
    ## Run 232 stress 9.916714e-05 
    ## ... Procrustes: rmse 0.0001636211  max resid 0.0002344593 
    ## ... Similar to previous best
    ## Run 233 stress 9.931641e-05 
    ## ... Procrustes: rmse 0.000194814  max resid 0.0002826033 
    ## ... Similar to previous best
    ## Run 234 stress 9.83293e-05 
    ## ... Procrustes: rmse 0.0001596196  max resid 0.0002353995 
    ## ... Similar to previous best
    ## Run 235 stress 9.91264e-05 
    ## ... Procrustes: rmse 0.0001652895  max resid 0.0002555792 
    ## ... Similar to previous best
    ## Run 236 stress 9.912312e-05 
    ## ... Procrustes: rmse 0.0001570165  max resid 0.0002407903 
    ## ... Similar to previous best
    ## Run 237 stress 9.876953e-05 
    ## ... Procrustes: rmse 0.0001840426  max resid 0.0002542901 
    ## ... Similar to previous best
    ## Run 238 stress 9.919035e-05 
    ## ... Procrustes: rmse 0.0001775717  max resid 0.0002571646 
    ## ... Similar to previous best
    ## Run 239 stress 8.551691e-05 
    ## ... Procrustes: rmse 0.0001152849  max resid 0.0001725689 
    ## ... Similar to previous best
    ## Run 240 stress 9.5566e-05 
    ## ... Procrustes: rmse 0.0001576815  max resid 0.0002481363 
    ## ... Similar to previous best
    ## Run 241 stress 8.902405e-05 
    ## ... Procrustes: rmse 0.0001126018  max resid 0.0001538766 
    ## ... Similar to previous best
    ## Run 242 stress 9.136014e-05 
    ## ... Procrustes: rmse 0.0001308536  max resid 0.0002368619 
    ## ... Similar to previous best
    ## Run 243 stress 9.900326e-05 
    ## ... Procrustes: rmse 0.0001599215  max resid 0.0002379851 
    ## ... Similar to previous best
    ## Run 244 stress 9.912833e-05 
    ## ... Procrustes: rmse 0.0001698713  max resid 0.0002414963 
    ## ... Similar to previous best
    ## Run 245 stress 9.120006e-05 
    ## ... Procrustes: rmse 0.0001351334  max resid 0.0001924155 
    ## ... Similar to previous best
    ## Run 246 stress 9.764019e-05 
    ## ... Procrustes: rmse 0.0001624124  max resid 0.0002491815 
    ## ... Similar to previous best
    ## Run 247 stress 9.918422e-05 
    ## ... Procrustes: rmse 0.0001928241  max resid 0.0002754701 
    ## ... Similar to previous best
    ## Run 248 stress 9.448305e-05 
    ## ... Procrustes: rmse 0.0001495483  max resid 0.0002320234 
    ## ... Similar to previous best
    ## Run 249 stress 9.669042e-05 
    ## ... Procrustes: rmse 0.0001718877  max resid 0.0002544301 
    ## ... Similar to previous best
    ## Run 250 stress 9.552105e-05 
    ## ... Procrustes: rmse 0.0001594385  max resid 0.0002317804 
    ## ... Similar to previous best
    ## Run 251 stress 9.819628e-05 
    ## ... Procrustes: rmse 0.0001911318  max resid 0.0003007673 
    ## ... Similar to previous best
    ## Run 252 stress 9.528891e-05 
    ## ... Procrustes: rmse 0.0001626572  max resid 0.0002419881 
    ## ... Similar to previous best
    ## Run 253 stress 9.950057e-05 
    ## ... Procrustes: rmse 0.0001745172  max resid 0.0002555941 
    ## ... Similar to previous best
    ## Run 254 stress 9.545149e-05 
    ## ... Procrustes: rmse 0.0001766618  max resid 0.0002493164 
    ## ... Similar to previous best
    ## Run 255 stress 0.249999 
    ## Run 256 stress 9.325133e-05 
    ## ... Procrustes: rmse 0.0001494355  max resid 0.000279149 
    ## ... Similar to previous best
    ## Run 257 stress 9.834251e-05 
    ## ... Procrustes: rmse 0.0001929106  max resid 0.0002814564 
    ## ... Similar to previous best
    ## Run 258 stress 9.65726e-05 
    ## ... Procrustes: rmse 0.0001736305  max resid 0.0002492085 
    ## ... Similar to previous best
    ## Run 259 stress 9.295154e-05 
    ## ... Procrustes: rmse 0.0001721164  max resid 0.0002474637 
    ## ... Similar to previous best
    ## Run 260 stress 7.461827e-05 
    ## ... Procrustes: rmse 0.0001503064  max resid 0.0002451853 
    ## ... Similar to previous best
    ## Run 261 stress 8.856623e-05 
    ## ... Procrustes: rmse 0.0001310173  max resid 0.0002137671 
    ## ... Similar to previous best
    ## Run 262 stress 9.945886e-05 
    ## ... Procrustes: rmse 0.0001791474  max resid 0.0002575814 
    ## ... Similar to previous best
    ## Run 263 stress 9.829258e-05 
    ## ... Procrustes: rmse 0.0001916164  max resid 0.0002649633 
    ## ... Similar to previous best
    ## Run 264 stress 9.200427e-05 
    ## ... Procrustes: rmse 0.0001346055  max resid 0.0002073225 
    ## ... Similar to previous best
    ## Run 265 stress 9.883084e-05 
    ## ... Procrustes: rmse 0.0001835315  max resid 0.0002688863 
    ## ... Similar to previous best
    ## Run 266 stress 9.75117e-05 
    ## ... Procrustes: rmse 0.00014832  max resid 0.0002315921 
    ## ... Similar to previous best
    ## Run 267 stress 9.435175e-05 
    ## ... Procrustes: rmse 0.0001478563  max resid 0.0002316134 
    ## ... Similar to previous best
    ## Run 268 stress 9.257741e-05 
    ## ... Procrustes: rmse 0.000144107  max resid 0.000226866 
    ## ... Similar to previous best
    ## Run 269 stress 9.784367e-05 
    ## ... Procrustes: rmse 0.0001865313  max resid 0.0002755127 
    ## ... Similar to previous best
    ## Run 270 stress 8.855568e-05 
    ## ... Procrustes: rmse 0.0001356373  max resid 0.0001962036 
    ## ... Similar to previous best
    ## Run 271 stress 9.592342e-05 
    ## ... Procrustes: rmse 0.000167456  max resid 0.0003023624 
    ## ... Similar to previous best
    ## Run 272 stress 9.483917e-05 
    ## ... Procrustes: rmse 0.0001890396  max resid 0.0002954726 
    ## ... Similar to previous best
    ## Run 273 stress 9.938919e-05 
    ## ... Procrustes: rmse 0.0001498109  max resid 0.0002466007 
    ## ... Similar to previous best
    ## Run 274 stress 0.2732559 
    ## Run 275 stress 9.316667e-05 
    ## ... Procrustes: rmse 0.0001429039  max resid 0.0002281316 
    ## ... Similar to previous best
    ## Run 276 stress 9.657235e-05 
    ## ... Procrustes: rmse 0.0001733392  max resid 0.0002451701 
    ## ... Similar to previous best
    ## Run 277 stress 9.457869e-05 
    ## ... Procrustes: rmse 0.0001511372  max resid 0.0002067306 
    ## ... Similar to previous best
    ## Run 278 stress 9.823767e-05 
    ## ... Procrustes: rmse 0.0001876233  max resid 0.0002702466 
    ## ... Similar to previous best
    ## Run 279 stress 9.51981e-05 
    ## ... Procrustes: rmse 0.0001656526  max resid 0.0002999646 
    ## ... Similar to previous best
    ## Run 280 stress 9.923455e-05 
    ## ... Procrustes: rmse 0.0001949358  max resid 0.0002633708 
    ## ... Similar to previous best
    ## Run 281 stress 9.389887e-05 
    ## ... Procrustes: rmse 0.0001728651  max resid 0.0002439441 
    ## ... Similar to previous best
    ## Run 282 stress 9.789544e-05 
    ## ... Procrustes: rmse 0.0001832059  max resid 0.0002633609 
    ## ... Similar to previous best
    ## Run 283 stress 9.412217e-05 
    ## ... Procrustes: rmse 0.0001274426  max resid 0.0002169687 
    ## ... Similar to previous best
    ## Run 284 stress 0.2280257 
    ## Run 285 stress 9.373129e-05 
    ## ... Procrustes: rmse 0.0001721662  max resid 0.0002480783 
    ## ... Similar to previous best
    ## Run 286 stress 9.390376e-05 
    ## ... Procrustes: rmse 0.0001515855  max resid 0.000226936 
    ## ... Similar to previous best
    ## Run 287 stress 9.193485e-05 
    ## ... Procrustes: rmse 0.0001422889  max resid 0.0002234529 
    ## ... Similar to previous best
    ## Run 288 stress 9.403191e-05 
    ## ... Procrustes: rmse 0.0001275993  max resid 0.000223348 
    ## ... Similar to previous best
    ## Run 289 stress 9.807766e-05 
    ## ... Procrustes: rmse 0.0001010977  max resid 0.0001597552 
    ## ... Similar to previous best
    ## Run 290 stress 9.146814e-05 
    ## ... Procrustes: rmse 0.0001469817  max resid 0.0002087116 
    ## ... Similar to previous best
    ## Run 291 stress 8.615179e-05 
    ## ... Procrustes: rmse 0.0001144132  max resid 0.0001793127 
    ## ... Similar to previous best
    ## Run 292 stress 9.158234e-05 
    ## ... Procrustes: rmse 0.0001235662  max resid 0.0002195469 
    ## ... Similar to previous best
    ## Run 293 stress 9.483893e-05 
    ## ... Procrustes: rmse 0.0001754786  max resid 0.0002488492 
    ## ... Similar to previous best
    ## Run 294 stress 9.889514e-05 
    ## ... Procrustes: rmse 0.0001326353  max resid 0.0002350577 
    ## ... Similar to previous best
    ## Run 295 stress 9.993429e-05 
    ## ... Procrustes: rmse 0.000166225  max resid 0.0002450508 
    ## ... Similar to previous best
    ## Run 296 stress 9.605874e-05 
    ## ... Procrustes: rmse 0.0001526612  max resid 0.0002305621 
    ## ... Similar to previous best
    ## Run 297 stress 9.510355e-05 
    ## ... Procrustes: rmse 0.0001662005  max resid 0.0003004054 
    ## ... Similar to previous best
    ## Run 298 stress 8.685287e-05 
    ## ... Procrustes: rmse 0.0001484775  max resid 0.0002360068 
    ## ... Similar to previous best
    ## Run 299 stress 9.555487e-05 
    ## ... Procrustes: rmse 0.0001546411  max resid 0.0002080986 
    ## ... Similar to previous best
    ## Run 300 stress 9.391144e-05 
    ## ... Procrustes: rmse 0.0001470374  max resid 0.0002331851 
    ## ... Similar to previous best
    ## Run 301 stress 9.358043e-05 
    ## ... Procrustes: rmse 0.0001634586  max resid 0.0002229077 
    ## ... Similar to previous best
    ## Run 302 stress 9.771406e-05 
    ## ... Procrustes: rmse 0.0001498751  max resid 0.0002834057 
    ## ... Similar to previous best
    ## Run 303 stress 9.559127e-05 
    ## ... Procrustes: rmse 0.0001595812  max resid 0.0002332031 
    ## ... Similar to previous best
    ## Run 304 stress 9.412554e-05 
    ## ... Procrustes: rmse 0.0001510667  max resid 0.0002262503 
    ## ... Similar to previous best
    ## Run 305 stress 9.586945e-05 
    ## ... Procrustes: rmse 0.000134279  max resid 0.0001945728 
    ## ... Similar to previous best
    ## Run 306 stress 9.20335e-05 
    ## ... Procrustes: rmse 0.0001518542  max resid 0.0002289843 
    ## ... Similar to previous best
    ## Run 307 stress 9.979242e-05 
    ## ... Procrustes: rmse 0.0001817017  max resid 0.0002501169 
    ## ... Similar to previous best
    ## Run 308 stress 9.652242e-05 
    ## ... Procrustes: rmse 0.0001553321  max resid 0.0002200118 
    ## ... Similar to previous best
    ## Run 309 stress 8.358358e-05 
    ## ... Procrustes: rmse 0.0001035286  max resid 0.0001520967 
    ## ... Similar to previous best
    ## Run 310 stress 9.666422e-05 
    ## ... Procrustes: rmse 0.0001482369  max resid 0.0002321137 
    ## ... Similar to previous best
    ## Run 311 stress 9.384706e-05 
    ## ... Procrustes: rmse 0.0001369077  max resid 0.0002724381 
    ## ... Similar to previous best
    ## Run 312 stress 9.999656e-05 
    ## ... Procrustes: rmse 0.0001865189  max resid 0.0002603158 
    ## ... Similar to previous best
    ## Run 313 stress 9.547632e-05 
    ## ... Procrustes: rmse 0.000187731  max resid 0.0002584932 
    ## ... Similar to previous best
    ## Run 314 stress 9.468684e-05 
    ## ... Procrustes: rmse 0.000159976  max resid 0.0002911346 
    ## ... Similar to previous best
    ## Run 315 stress 9.156821e-05 
    ## ... Procrustes: rmse 0.0001694138  max resid 0.0002453811 
    ## ... Similar to previous best
    ## Run 316 stress 6.21026e-05 
    ## ... Procrustes: rmse 8.334402e-05  max resid 0.0001691899 
    ## ... Similar to previous best
    ## Run 317 stress 9.692587e-05 
    ## ... Procrustes: rmse 0.0001335209  max resid 0.0002271918 
    ## ... Similar to previous best
    ## Run 318 stress 9.870386e-05 
    ## ... Procrustes: rmse 0.0001560622  max resid 0.0002458382 
    ## ... Similar to previous best
    ## Run 319 stress 9.699638e-05 
    ## ... Procrustes: rmse 0.0001573813  max resid 0.0002234173 
    ## ... Similar to previous best
    ## Run 320 stress 8.756354e-05 
    ## ... Procrustes: rmse 0.0001519746  max resid 0.000282197 
    ## ... Similar to previous best
    ## Run 321 stress 9.843205e-05 
    ## ... Procrustes: rmse 0.0001743118  max resid 0.0002466102 
    ## ... Similar to previous best
    ## Run 322 stress 9.847091e-05 
    ## ... Procrustes: rmse 0.0001706616  max resid 0.0003055026 
    ## ... Similar to previous best
    ## Run 323 stress 8.115706e-05 
    ## ... Procrustes: rmse 0.0001151154  max resid 0.0001761646 
    ## ... Similar to previous best
    ## Run 324 stress 8.801747e-05 
    ## ... Procrustes: rmse 0.0001398741  max resid 0.000215902 
    ## ... Similar to previous best
    ## Run 325 stress 8.203677e-05 
    ## ... Procrustes: rmse 0.000123237  max resid 0.0001753631 
    ## ... Similar to previous best
    ## Run 326 stress 8.699284e-05 
    ## ... Procrustes: rmse 0.0001322055  max resid 0.0002084152 
    ## ... Similar to previous best
    ## Run 327 stress 9.554855e-05 
    ## ... Procrustes: rmse 0.0001634472  max resid 0.0002353337 
    ## ... Similar to previous best
    ## Run 328 stress 8.799569e-05 
    ## ... Procrustes: rmse 0.0001294642  max resid 0.0002213081 
    ## ... Similar to previous best
    ## Run 329 stress 9.946158e-05 
    ## ... Procrustes: rmse 0.0001842502  max resid 0.0002705617 
    ## ... Similar to previous best
    ## Run 330 stress 9.385527e-05 
    ## ... Procrustes: rmse 0.0001092651  max resid 0.0001994064 
    ## ... Similar to previous best
    ## Run 331 stress 9.891723e-05 
    ## ... Procrustes: rmse 0.0001917272  max resid 0.0002746043 
    ## ... Similar to previous best
    ## Run 332 stress 9.260143e-05 
    ## ... Procrustes: rmse 0.0001438605  max resid 0.0002293207 
    ## ... Similar to previous best
    ## Run 333 stress 9.616713e-05 
    ## ... Procrustes: rmse 0.0001884893  max resid 0.0002586849 
    ## ... Similar to previous best
    ## Run 334 stress 9.27521e-05 
    ## ... Procrustes: rmse 0.0001561269  max resid 0.0002273161 
    ## ... Similar to previous best
    ## Run 335 stress 0.2280257 
    ## Run 336 stress 9.912847e-05 
    ## ... Procrustes: rmse 0.0001587516  max resid 0.0002432318 
    ## ... Similar to previous best
    ## Run 337 stress 9.985271e-05 
    ## ... Procrustes: rmse 0.0001837794  max resid 0.0002556565 
    ## ... Similar to previous best
    ## Run 338 stress 9.349104e-05 
    ## ... Procrustes: rmse 0.0001444775  max resid 0.0002024354 
    ## ... Similar to previous best
    ## Run 339 stress 9.882305e-05 
    ## ... Procrustes: rmse 0.0001626403  max resid 0.0002363148 
    ## ... Similar to previous best
    ## Run 340 stress 9.939281e-05 
    ## ... Procrustes: rmse 0.0001819446  max resid 0.0002542181 
    ## ... Similar to previous best
    ## Run 341 stress 8.943893e-05 
    ## ... Procrustes: rmse 0.0001192597  max resid 0.0002658903 
    ## ... Similar to previous best
    ## Run 342 stress 9.398672e-05 
    ## ... Procrustes: rmse 0.0001634373  max resid 0.0002982005 
    ## ... Similar to previous best
    ## Run 343 stress 9.899312e-05 
    ## ... Procrustes: rmse 0.0001670497  max resid 0.0002319478 
    ## ... Similar to previous best
    ## Run 344 stress 9.570401e-05 
    ## ... Procrustes: rmse 0.0001164972  max resid 0.0001937535 
    ## ... Similar to previous best
    ## Run 345 stress 9.545354e-05 
    ## ... Procrustes: rmse 0.0001597137  max resid 0.0002303253 
    ## ... Similar to previous best
    ## Run 346 stress 9.806368e-05 
    ## ... Procrustes: rmse 0.0001471368  max resid 0.0002232328 
    ## ... Similar to previous best
    ## Run 347 stress 9.29148e-05 
    ## ... Procrustes: rmse 0.0001457686  max resid 0.000194656 
    ## ... Similar to previous best
    ## Run 348 stress 9.992808e-05 
    ## ... Procrustes: rmse 0.000160974  max resid 0.000229011 
    ## ... Similar to previous best
    ## Run 349 stress 9.107458e-05 
    ## ... Procrustes: rmse 0.0001583303  max resid 0.0002917028 
    ## ... Similar to previous best
    ## Run 350 stress 8.627994e-05 
    ## ... Procrustes: rmse 0.0001554673  max resid 0.0002335827 
    ## ... Similar to previous best
    ## Run 351 stress 9.693938e-05 
    ## ... Procrustes: rmse 0.0001883612  max resid 0.0002973852 
    ## ... Similar to previous best
    ## Run 352 stress 9.083373e-05 
    ## ... Procrustes: rmse 0.0001529148  max resid 0.0002233642 
    ## ... Similar to previous best
    ## Run 353 stress 9.366947e-05 
    ## ... Procrustes: rmse 0.0001374172  max resid 0.0002431469 
    ## ... Similar to previous best
    ## Run 354 stress 9.614336e-05 
    ## ... Procrustes: rmse 0.0001875193  max resid 0.0002750508 
    ## ... Similar to previous best
    ## Run 355 stress 9.888211e-05 
    ## ... Procrustes: rmse 0.0001400486  max resid 0.0002410549 
    ## ... Similar to previous best
    ## Run 356 stress 9.719044e-05 
    ## ... Procrustes: rmse 0.0001883789  max resid 0.0002715435 
    ## ... Similar to previous best
    ## Run 357 stress 9.623798e-05 
    ## ... Procrustes: rmse 0.0001649095  max resid 0.0002992894 
    ## ... Similar to previous best
    ## Run 358 stress 9.885e-05 
    ## ... Procrustes: rmse 0.0001878328  max resid 0.0002965997 
    ## ... Similar to previous best
    ## Run 359 stress 9.368251e-05 
    ## ... Procrustes: rmse 0.0001431421  max resid 0.0002148312 
    ## ... Similar to previous best
    ## Run 360 stress 9.45334e-05 
    ## ... Procrustes: rmse 0.000149646  max resid 0.0002329655 
    ## ... Similar to previous best
    ## Run 361 stress 9.488527e-05 
    ## ... Procrustes: rmse 0.0001534323  max resid 0.0002210495 
    ## ... Similar to previous best
    ## Run 362 stress 9.224345e-05 
    ## ... Procrustes: rmse 0.0001240217  max resid 0.0002118022 
    ## ... Similar to previous best
    ## Run 363 stress 9.866105e-05 
    ## ... Procrustes: rmse 0.0001885585  max resid 0.0002608886 
    ## ... Similar to previous best
    ## Run 364 stress 9.989551e-05 
    ## ... Procrustes: rmse 0.0001327656  max resid 0.0001885234 
    ## ... Similar to previous best
    ## Run 365 stress 9.157051e-05 
    ## ... Procrustes: rmse 0.0001415492  max resid 0.0002262571 
    ## ... Similar to previous best
    ## Run 366 stress 9.386109e-05 
    ## ... Procrustes: rmse 0.0001738327  max resid 0.000247395 
    ## ... Similar to previous best
    ## Run 367 stress 9.745507e-05 
    ## ... Procrustes: rmse 0.0001810035  max resid 0.0002664196 
    ## ... Similar to previous best
    ## Run 368 stress 9.696943e-05 
    ## ... Procrustes: rmse 0.0001479573  max resid 0.0001976701 
    ## ... Similar to previous best
    ## Run 369 stress 9.875978e-05 
    ## ... Procrustes: rmse 0.0001608864  max resid 0.0002373922 
    ## ... Similar to previous best
    ## Run 370 stress 9.926015e-05 
    ## ... Procrustes: rmse 0.0001574292  max resid 0.0002346074 
    ## ... Similar to previous best
    ## Run 371 stress 9.363098e-05 
    ## ... Procrustes: rmse 0.000167551  max resid 0.0002434883 
    ## ... Similar to previous best
    ## Run 372 stress 9.314393e-05 
    ## ... Procrustes: rmse 0.0001711637  max resid 0.0002589514 
    ## ... Similar to previous best
    ## Run 373 stress 9.862346e-05 
    ## ... Procrustes: rmse 0.000169995  max resid 0.0002461297 
    ## ... Similar to previous best
    ## Run 374 stress 9.669487e-05 
    ## ... Procrustes: rmse 0.0001845944  max resid 0.0002927098 
    ## ... Similar to previous best
    ## Run 375 stress 9.972517e-05 
    ## ... Procrustes: rmse 0.0001402481  max resid 0.0002215946 
    ## ... Similar to previous best
    ## Run 376 stress 9.495823e-05 
    ## ... Procrustes: rmse 0.0001568423  max resid 0.0002304359 
    ## ... Similar to previous best
    ## Run 377 stress 9.35597e-05 
    ## ... Procrustes: rmse 0.0001592498  max resid 0.0002299393 
    ## ... Similar to previous best
    ## Run 378 stress 9.558458e-05 
    ## ... Procrustes: rmse 0.0001431893  max resid 0.0002139597 
    ## ... Similar to previous best
    ## Run 379 stress 9.900018e-05 
    ## ... Procrustes: rmse 0.0001720117  max resid 0.0002598605 
    ## ... Similar to previous best
    ## Run 380 stress 8.146896e-05 
    ## ... Procrustes: rmse 0.000132142  max resid 0.0002068222 
    ## ... Similar to previous best
    ## Run 381 stress 9.920172e-05 
    ## ... Procrustes: rmse 0.0001866301  max resid 0.0002920143 
    ## ... Similar to previous best
    ## Run 382 stress 9.659936e-05 
    ## ... Procrustes: rmse 0.0001850193  max resid 0.0002671326 
    ## ... Similar to previous best
    ## Run 383 stress 9.526679e-05 
    ## ... Procrustes: rmse 0.0001810552  max resid 0.000254561 
    ## ... Similar to previous best
    ## Run 384 stress 9.395223e-05 
    ## ... Procrustes: rmse 0.0001458805  max resid 0.0002073013 
    ## ... Similar to previous best
    ## Run 385 stress 9.683914e-05 
    ## ... Procrustes: rmse 0.0001561309  max resid 0.0002224754 
    ## ... Similar to previous best
    ## Run 386 stress 9.719225e-05 
    ## ... Procrustes: rmse 0.0001614544  max resid 0.0002479917 
    ## ... Similar to previous best
    ## Run 387 stress 9.732165e-05 
    ## ... Procrustes: rmse 0.0001794017  max resid 0.000251113 
    ## ... Similar to previous best
    ## Run 388 stress 9.515093e-05 
    ## ... Procrustes: rmse 0.0001189594  max resid 0.0001747081 
    ## ... Similar to previous best
    ## Run 389 stress 9.962365e-05 
    ## ... Procrustes: rmse 0.0001735344  max resid 0.0002412936 
    ## ... Similar to previous best
    ## Run 390 stress 9.250548e-05 
    ## ... Procrustes: rmse 0.0001288333  max resid 0.000179589 
    ## ... Similar to previous best
    ## Run 391 stress 9.623066e-05 
    ## ... Procrustes: rmse 0.0001861378  max resid 0.0002693958 
    ## ... Similar to previous best
    ## Run 392 stress 9.768717e-05 
    ## ... Procrustes: rmse 0.0001820008  max resid 0.000253255 
    ## ... Similar to previous best
    ## Run 393 stress 9.403049e-05 
    ## ... Procrustes: rmse 0.0001480708  max resid 0.0002307946 
    ## ... Similar to previous best
    ## Run 394 stress 9.967815e-05 
    ## ... Procrustes: rmse 0.0001587612  max resid 0.0002092732 
    ## ... Similar to previous best
    ## Run 395 stress 9.610987e-05 
    ## ... Procrustes: rmse 0.0001226309  max resid 0.000251956 
    ## ... Similar to previous best
    ## Run 396 stress 9.44961e-05 
    ## ... Procrustes: rmse 0.0001467399  max resid 0.0002032391 
    ## ... Similar to previous best
    ## Run 397 stress 9.935191e-05 
    ## ... Procrustes: rmse 0.0001508367  max resid 0.0002112219 
    ## ... Similar to previous best
    ## Run 398 stress 9.422025e-05 
    ## ... Procrustes: rmse 0.0001608171  max resid 0.0002302384 
    ## ... Similar to previous best
    ## Run 399 stress 9.654331e-05 
    ## ... Procrustes: rmse 0.000189702  max resid 0.0002993505 
    ## ... Similar to previous best
    ## Run 400 stress 9.361596e-05 
    ## ... Procrustes: rmse 0.0001414515  max resid 0.000219429 
    ## ... Similar to previous best
    ## Run 401 stress 9.089788e-05 
    ## ... Procrustes: rmse 0.0001468958  max resid 0.0002200312 
    ## ... Similar to previous best
    ## Run 402 stress 9.558786e-05 
    ## ... Procrustes: rmse 0.0001535173  max resid 0.0002176762 
    ## ... Similar to previous best
    ## Run 403 stress 9.569052e-05 
    ## ... Procrustes: rmse 0.0001575316  max resid 0.0002337286 
    ## ... Similar to previous best
    ## Run 404 stress 9.422938e-05 
    ## ... Procrustes: rmse 0.0001159484  max resid 0.0001690338 
    ## ... Similar to previous best
    ## Run 405 stress 9.537554e-05 
    ## ... Procrustes: rmse 0.0001681546  max resid 0.0002404831 
    ## ... Similar to previous best
    ## Run 406 stress 9.851461e-05 
    ## ... Procrustes: rmse 0.0001398471  max resid 0.0002091736 
    ## ... Similar to previous best
    ## Run 407 stress 8.685486e-05 
    ## ... Procrustes: rmse 0.0001393578  max resid 0.0002126701 
    ## ... Similar to previous best
    ## Run 408 stress 8.976204e-05 
    ## ... Procrustes: rmse 0.0001576454  max resid 0.0002351217 
    ## ... Similar to previous best
    ## Run 409 stress 9.416257e-05 
    ## ... Procrustes: rmse 0.0001679081  max resid 0.0002437563 
    ## ... Similar to previous best
    ## Run 410 stress 9.291837e-05 
    ## ... Procrustes: rmse 0.000156579  max resid 0.0002247906 
    ## ... Similar to previous best
    ## Run 411 stress 9.809897e-05 
    ## ... Procrustes: rmse 0.0001727254  max resid 0.0002593537 
    ## ... Similar to previous best
    ## Run 412 stress 9.997993e-05 
    ## ... Procrustes: rmse 0.0001948768  max resid 0.0002771108 
    ## ... Similar to previous best
    ## Run 413 stress 9.039281e-05 
    ## ... Procrustes: rmse 0.0001353183  max resid 0.0002100092 
    ## ... Similar to previous best
    ## Run 414 stress 9.882683e-05 
    ## ... Procrustes: rmse 0.00016287  max resid 0.0002155109 
    ## ... Similar to previous best
    ## Run 415 stress 9.976616e-05 
    ## ... Procrustes: rmse 0.0001667494  max resid 0.0002544009 
    ## ... Similar to previous best
    ## Run 416 stress 9.710254e-05 
    ## ... Procrustes: rmse 0.0001929478  max resid 0.0003030748 
    ## ... Similar to previous best
    ## Run 417 stress 9.187237e-05 
    ## ... Procrustes: rmse 0.0001706871  max resid 0.0002444145 
    ## ... Similar to previous best
    ## Run 418 stress 9.607931e-05 
    ## ... Procrustes: rmse 0.000181447  max resid 0.0002549588 
    ## ... Similar to previous best
    ## Run 419 stress 9.844573e-05 
    ## ... Procrustes: rmse 0.0001565135  max resid 0.000238442 
    ## ... Similar to previous best
    ## Run 420 stress 9.329475e-05 
    ## ... Procrustes: rmse 0.0001759976  max resid 0.0002510734 
    ## ... Similar to previous best
    ## Run 421 stress 8.453141e-05 
    ## ... Procrustes: rmse 0.0001384035  max resid 0.0002094699 
    ## ... Similar to previous best
    ## Run 422 stress 9.802469e-05 
    ## ... Procrustes: rmse 0.0001939284  max resid 0.0003048814 
    ## ... Similar to previous best
    ## Run 423 stress 9.703298e-05 
    ## ... Procrustes: rmse 0.0001499144  max resid 0.0002306897 
    ## ... Similar to previous best
    ## Run 424 stress 8.621611e-05 
    ## ... Procrustes: rmse 0.0001564412  max resid 0.0002321026 
    ## ... Similar to previous best
    ## Run 425 stress 9.517517e-05 
    ## ... Procrustes: rmse 0.0001343839  max resid 0.0002015209 
    ## ... Similar to previous best
    ## Run 426 stress 8.747067e-05 
    ## ... Procrustes: rmse 0.0001536228  max resid 0.0002322707 
    ## ... Similar to previous best
    ## Run 427 stress 9.270617e-05 
    ## ... Procrustes: rmse 0.0001450109  max resid 0.0002288008 
    ## ... Similar to previous best
    ## Run 428 stress 9.971781e-05 
    ## ... Procrustes: rmse 0.0001479245  max resid 0.0002021313 
    ## ... Similar to previous best
    ## Run 429 stress 9.515586e-05 
    ## ... Procrustes: rmse 0.0001698587  max resid 0.0002504567 
    ## ... Similar to previous best
    ## Run 430 stress 7.748567e-05 
    ## ... Procrustes: rmse 0.0001319145  max resid 0.0002113876 
    ## ... Similar to previous best
    ## Run 431 stress 8.634467e-05 
    ## ... Procrustes: rmse 0.0001494661  max resid 0.000226226 
    ## ... Similar to previous best
    ## Run 432 stress 9.430796e-05 
    ## ... Procrustes: rmse 0.0001471301  max resid 0.000224904 
    ## ... Similar to previous best
    ## Run 433 stress 0.2857987 
    ## Run 434 stress 9.371003e-05 
    ## ... Procrustes: rmse 0.0001268873  max resid 0.0002226562 
    ## ... Similar to previous best
    ## Run 435 stress 9.992765e-05 
    ## ... Procrustes: rmse 0.0001918831  max resid 0.0002719277 
    ## ... Similar to previous best
    ## Run 436 stress 9.774892e-05 
    ## ... Procrustes: rmse 0.0001823045  max resid 0.0002537768 
    ## ... Similar to previous best
    ## Run 437 stress 9.296262e-05 
    ## ... Procrustes: rmse 0.0001343308  max resid 0.0002367179 
    ## ... Similar to previous best
    ## Run 438 stress 9.986502e-05 
    ## ... Procrustes: rmse 0.000186747  max resid 0.0002705715 
    ## ... Similar to previous best
    ## Run 439 stress 9.835095e-05 
    ## ... Procrustes: rmse 0.0001903191  max resid 0.0002734821 
    ## ... Similar to previous best
    ## Run 440 stress 8.636974e-05 
    ## ... Procrustes: rmse 0.0001047333  max resid 0.0001771473 
    ## ... Similar to previous best
    ## Run 441 stress 9.860925e-05 
    ## ... Procrustes: rmse 0.0001836607  max resid 0.0002570444 
    ## ... Similar to previous best
    ## Run 442 stress 9.037631e-05 
    ## ... Procrustes: rmse 0.0001643687  max resid 0.0002383471 
    ## ... Similar to previous best
    ## Run 443 stress 9.083799e-05 
    ## ... Procrustes: rmse 0.0001504216  max resid 0.0002282693 
    ## ... Similar to previous best
    ## Run 444 stress 9.589654e-05 
    ## ... Procrustes: rmse 0.00015439  max resid 0.000208287 
    ## ... Similar to previous best
    ## Run 445 stress 9.603594e-05 
    ## ... Procrustes: rmse 0.0001489205  max resid 0.0002126809 
    ## ... Similar to previous best
    ## Run 446 stress 9.25288e-05 
    ## ... Procrustes: rmse 0.0001215939  max resid 0.0001901816 
    ## ... Similar to previous best
    ## Run 447 stress 9.69901e-05 
    ## ... Procrustes: rmse 0.0001587601  max resid 0.0002341636 
    ## ... Similar to previous best
    ## Run 448 stress 8.446225e-05 
    ## ... Procrustes: rmse 0.0001613765  max resid 0.000236945 
    ## ... Similar to previous best
    ## Run 449 stress 9.821129e-05 
    ## ... Procrustes: rmse 0.0001672831  max resid 0.0002363632 
    ## ... Similar to previous best
    ## Run 450 stress 9.741417e-05 
    ## ... Procrustes: rmse 0.0001551235  max resid 0.0002372825 
    ## ... Similar to previous best
    ## Run 451 stress 9.393054e-05 
    ## ... Procrustes: rmse 0.0001614066  max resid 0.0002493159 
    ## ... Similar to previous best
    ## Run 452 stress 9.158089e-05 
    ## ... Procrustes: rmse 0.0001597353  max resid 0.0002209601 
    ## ... Similar to previous best
    ## Run 453 stress 9.407491e-05 
    ## ... Procrustes: rmse 0.0001496015  max resid 0.0002305774 
    ## ... Similar to previous best
    ## Run 454 stress 8.681943e-05 
    ## ... Procrustes: rmse 0.0001159845  max resid 0.000189142 
    ## ... Similar to previous best
    ## Run 455 stress 9.255288e-05 
    ## ... Procrustes: rmse 0.0001239471  max resid 0.0002156721 
    ## ... Similar to previous best
    ## Run 456 stress 9.990632e-05 
    ## ... Procrustes: rmse 0.0001962134  max resid 0.0003082166 
    ## ... Similar to previous best
    ## Run 457 stress 9.931701e-05 
    ## ... Procrustes: rmse 0.0001651664  max resid 0.000239539 
    ## ... Similar to previous best
    ## Run 458 stress 9.495815e-05 
    ## ... Procrustes: rmse 0.0001466962  max resid 0.000235417 
    ## ... Similar to previous best
    ## Run 459 stress 9.707065e-05 
    ## ... Procrustes: rmse 0.0001588125  max resid 0.0002403654 
    ## ... Similar to previous best
    ## Run 460 stress 9.471927e-05 
    ## ... Procrustes: rmse 0.0001556143  max resid 0.0002357374 
    ## ... Similar to previous best
    ## Run 461 stress 8.049281e-05 
    ## ... Procrustes: rmse 0.0001279633  max resid 0.0002062089 
    ## ... Similar to previous best
    ## Run 462 stress 9.550511e-05 
    ## ... Procrustes: rmse 0.0001767226  max resid 0.0002653824 
    ## ... Similar to previous best
    ## Run 463 stress 7.259243e-05 
    ## ... Procrustes: rmse 9.546768e-05  max resid 0.0001743721 
    ## ... Similar to previous best
    ## Run 464 stress 9.590468e-05 
    ## ... Procrustes: rmse 0.0001416655  max resid 0.0002201079 
    ## ... Similar to previous best
    ## Run 465 stress 9.201612e-05 
    ## ... Procrustes: rmse 0.0001695472  max resid 0.0002464448 
    ## ... Similar to previous best
    ## Run 466 stress 9.663522e-05 
    ## ... Procrustes: rmse 0.0001556743  max resid 0.0002188827 
    ## ... Similar to previous best
    ## Run 467 stress 9.70993e-05 
    ## ... Procrustes: rmse 0.0001513993  max resid 0.0001997059 
    ## ... Similar to previous best
    ## Run 468 stress 9.20503e-05 
    ## ... Procrustes: rmse 0.0001678896  max resid 0.0002455309 
    ## ... Similar to previous best
    ## Run 469 stress 8.908408e-05 
    ## ... Procrustes: rmse 0.0001503037  max resid 0.0002383243 
    ## ... Similar to previous best
    ## Run 470 stress 9.941466e-05 
    ## ... Procrustes: rmse 0.0001377229  max resid 0.0002186626 
    ## ... Similar to previous best
    ## Run 471 stress 9.926986e-05 
    ## ... Procrustes: rmse 0.0001880439  max resid 0.000267796 
    ## ... Similar to previous best
    ## Run 472 stress 8.753481e-05 
    ## ... Procrustes: rmse 0.0001585633  max resid 0.0002325834 
    ## ... Similar to previous best
    ## Run 473 stress 7.067529e-05 
    ## ... Procrustes: rmse 9.066996e-05  max resid 0.0001762451 
    ## ... Similar to previous best
    ## Run 474 stress 9.435215e-05 
    ## ... Procrustes: rmse 0.0001542615  max resid 0.0002841129 
    ## ... Similar to previous best
    ## Run 475 stress 8.748809e-05 
    ## ... Procrustes: rmse 0.0001199418  max resid 0.0001620972 
    ## ... Similar to previous best
    ## Run 476 stress 9.488481e-05 
    ## ... Procrustes: rmse 0.0001896586  max resid 0.0002989586 
    ## ... Similar to previous best
    ## Run 477 stress 9.117495e-05 
    ## ... Procrustes: rmse 0.0001349516  max resid 0.0001848788 
    ## ... Similar to previous best
    ## Run 478 stress 9.933717e-05 
    ## ... Procrustes: rmse 0.00018563  max resid 0.0002673267 
    ## ... Similar to previous best
    ## Run 479 stress 9.472498e-05 
    ## ... Procrustes: rmse 0.0001796312  max resid 0.000253554 
    ## ... Similar to previous best
    ## Run 480 stress 7.72422e-05 
    ## ... Procrustes: rmse 0.0001045435  max resid 0.0001773966 
    ## ... Similar to previous best
    ## Run 481 stress 9.989135e-05 
    ## ... Procrustes: rmse 0.0002001873  max resid 0.0003123119 
    ## ... Similar to previous best
    ## Run 482 stress 9.776928e-05 
    ## ... Procrustes: rmse 0.0001497109  max resid 0.0002112631 
    ## ... Similar to previous best
    ## Run 483 stress 9.63944e-05 
    ## ... Procrustes: rmse 0.000165389  max resid 0.0002335695 
    ## ... Similar to previous best
    ## Run 484 stress 8.857011e-05 
    ## ... Procrustes: rmse 0.0001262849  max resid 0.0002083165 
    ## ... Similar to previous best
    ## Run 485 stress 9.993265e-05 
    ## ... Procrustes: rmse 0.0001862421  max resid 0.0002588066 
    ## ... Similar to previous best
    ## Run 486 stress 9.373725e-05 
    ## ... Procrustes: rmse 0.0001590416  max resid 0.0002284848 
    ## ... Similar to previous best
    ## Run 487 stress 9.500105e-05 
    ## ... Procrustes: rmse 0.0001443194  max resid 0.0002030236 
    ## ... Similar to previous best
    ## Run 488 stress 9.167993e-05 
    ## ... Procrustes: rmse 0.0001796895  max resid 0.0002519412 
    ## ... Similar to previous best
    ## Run 489 stress 9.949385e-05 
    ## ... Procrustes: rmse 0.0001836495  max resid 0.0002585993 
    ## ... Similar to previous best
    ## Run 490 stress 9.772262e-05 
    ## ... Procrustes: rmse 0.0001497932  max resid 0.000213349 
    ## ... Similar to previous best
    ## Run 491 stress 9.499319e-05 
    ## ... Procrustes: rmse 0.0001554728  max resid 0.0002409077 
    ## ... Similar to previous best
    ## Run 492 stress 9.128643e-05 
    ## ... Procrustes: rmse 0.0001825914  max resid 0.0002869096 
    ## ... Similar to previous best
    ## Run 493 stress 9.04093e-05 
    ## ... Procrustes: rmse 0.000103039  max resid 0.0001755775 
    ## ... Similar to previous best
    ## Run 494 stress 9.280804e-05 
    ## ... Procrustes: rmse 0.0001438691  max resid 0.0002011523 
    ## ... Similar to previous best
    ## Run 495 stress 9.882841e-05 
    ## ... Procrustes: rmse 0.000187493  max resid 0.0002695995 
    ## ... Similar to previous best
    ## Run 496 stress 8.954012e-05 
    ## ... Procrustes: rmse 0.0001292586  max resid 0.0001834355 
    ## ... Similar to previous best
    ## Run 497 stress 9.773938e-05 
    ## ... Procrustes: rmse 0.0001954174  max resid 0.0003062476 
    ## ... Similar to previous best
    ## Run 498 stress 9.819601e-05 
    ## ... Procrustes: rmse 0.0001766139  max resid 0.0002490161 
    ## ... Similar to previous best
    ## Run 499 stress 9.65106e-05 
    ## ... Procrustes: rmse 0.0001482077  max resid 0.0001928119 
    ## ... Similar to previous best
    ## Run 500 stress 9.911217e-05 
    ## ... Procrustes: rmse 0.0001672165  max resid 0.0002944883 
    ## ... Similar to previous best
    ## *** Best solution repeated 300 times

    ## Warning in metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
SD_beta_env_M_NMDS <- metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 0.00122774 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0008914  max resid 0.001236069 
    ## ... Similar to previous best
    ## Run 2 stress 0.001252834 
    ## ... Procrustes: rmse 0.001780581  max resid 0.00290898 
    ## ... Similar to previous best
    ## Run 3 stress 9.638293e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04231483  max resid 0.05830261 
    ## Run 4 stress 0.0004892091 
    ## ... Procrustes: rmse 0.01618147  max resid 0.02221873 
    ## Run 5 stress 0.1990774 
    ## Run 6 stress 0.000374026 
    ## ... Procrustes: rmse 0.01407394  max resid 0.01931149 
    ## Run 7 stress 0.001389718 
    ## Run 8 stress 9.857057e-05 
    ## ... Procrustes: rmse 0.0001747454  max resid 0.0003002355 
    ## ... Similar to previous best
    ## Run 9 stress 0.001375242 
    ## Run 10 stress 0.001086718 
    ## Run 11 stress 9.152607e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001536589  max resid 0.0002846703 
    ## ... Similar to previous best
    ## Run 12 stress 0.1491957 
    ## Run 13 stress 0.001119895 
    ## Run 14 stress 9.748143e-05 
    ## ... Procrustes: rmse 0.0001731329  max resid 0.0003267018 
    ## ... Similar to previous best
    ## Run 15 stress 9.840994e-05 
    ## ... Procrustes: rmse 0.0001183101  max resid 0.0001828553 
    ## ... Similar to previous best
    ## Run 16 stress 0.2795549 
    ## Run 17 stress 0.001287459 
    ## Run 18 stress 9.351897e-05 
    ## ... Procrustes: rmse 0.0001988445  max resid 0.0003411386 
    ## ... Similar to previous best
    ## Run 19 stress 0.001158352 
    ## Run 20 stress 9.565334e-05 
    ## ... Procrustes: rmse 0.0001649673  max resid 0.0002789865 
    ## ... Similar to previous best
    ## Run 21 stress 0.001221775 
    ## Run 22 stress 8.711467e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001612593  max resid 0.0002975446 
    ## ... Similar to previous best
    ## Run 23 stress 0.001280058 
    ## Run 24 stress 0.2528294 
    ## Run 25 stress 8.778084e-05 
    ## ... Procrustes: rmse 0.0001838488  max resid 0.0003634357 
    ## ... Similar to previous best
    ## Run 26 stress 0.2528294 
    ## Run 27 stress 0.0004933566 
    ## ... Procrustes: rmse 0.0161316  max resid 0.02208213 
    ## Run 28 stress 0.001265016 
    ## Run 29 stress 0.001371104 
    ## Run 30 stress 0.2851274 
    ## Run 31 stress 0.1776018 
    ## Run 32 stress 0.0004954065 
    ## ... Procrustes: rmse 0.0162583  max resid 0.0222494 
    ## Run 33 stress 0.0004961697 
    ## ... Procrustes: rmse 0.01627035  max resid 0.02226445 
    ## Run 34 stress 0.00110943 
    ## Run 35 stress 0.1491957 
    ## Run 36 stress 0.001146149 
    ## Run 37 stress 0.001325038 
    ## Run 38 stress 0.1776011 
    ## Run 39 stress 0.0001912205 
    ## ... Procrustes: rmse 0.01007297  max resid 0.01371444 
    ## Run 40 stress 9.192897e-05 
    ## ... Procrustes: rmse 0.0001585698  max resid 0.0002814574 
    ## ... Similar to previous best
    ## Run 41 stress 8.328652e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001801869  max resid 0.0003685552 
    ## ... Similar to previous best
    ## Run 42 stress 0.001126125 
    ## Run 43 stress 0.0004646047 
    ## ... Procrustes: rmse 0.01580171  max resid 0.02158928 
    ## Run 44 stress 8.708324e-05 
    ## ... Procrustes: rmse 0.0001708679  max resid 0.0002573854 
    ## ... Similar to previous best
    ## Run 45 stress 0.0003434169 
    ## ... Procrustes: rmse 0.01358545  max resid 0.0185299 
    ## Run 46 stress 9.155469e-05 
    ## ... Procrustes: rmse 0.0001817648  max resid 0.0003398307 
    ## ... Similar to previous best
    ## Run 47 stress 9.680761e-05 
    ## ... Procrustes: rmse 9.474841e-05  max resid 0.0001374337 
    ## ... Similar to previous best
    ## Run 48 stress 0.001359683 
    ## Run 49 stress 0.00133854 
    ## Run 50 stress 0.1491957 
    ## Run 51 stress 9.15066e-05 
    ## ... Procrustes: rmse 0.0001891115  max resid 0.0003517398 
    ## ... Similar to previous best
    ## Run 52 stress 9.94258e-05 
    ## ... Procrustes: rmse 0.0001260411  max resid 0.0002103493 
    ## ... Similar to previous best
    ## Run 53 stress 0.001395003 
    ## Run 54 stress 0.1776012 
    ## Run 55 stress 0.0006877413 
    ## Run 56 stress 0.001310299 
    ## Run 57 stress 0.001058101 
    ## Run 58 stress 0.001417144 
    ## Run 59 stress 0.1491957 
    ## Run 60 stress 0.0001324972 
    ## ... Procrustes: rmse 0.01380233  max resid 0.01896468 
    ## Run 61 stress 0.2440272 
    ## Run 62 stress 0.001334222 
    ## Run 63 stress 0.001328109 
    ## Run 64 stress 0.001413261 
    ## Run 65 stress 0.2570108 
    ## Run 66 stress 0.00130648 
    ## Run 67 stress 0.0004663691 
    ## ... Procrustes: rmse 0.01568564  max resid 0.02142723 
    ## Run 68 stress 8.498412e-05 
    ## ... Procrustes: rmse 0.0001730207  max resid 0.0002868255 
    ## ... Similar to previous best
    ## Run 69 stress 8.677998e-05 
    ## ... Procrustes: rmse 2.466496e-05  max resid 3.472606e-05 
    ## ... Similar to previous best
    ## Run 70 stress 9.22722e-05 
    ## ... Procrustes: rmse 0.0001246819  max resid 0.0002213862 
    ## ... Similar to previous best
    ## Run 71 stress 9.984177e-05 
    ## ... Procrustes: rmse 9.459464e-05  max resid 0.0001525659 
    ## ... Similar to previous best
    ## Run 72 stress 0.1491957 
    ## Run 73 stress 6.929495e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001655  max resid 0.0002411352 
    ## ... Similar to previous best
    ## Run 74 stress 0.001179409 
    ## Run 75 stress 0.1776012 
    ## Run 76 stress 8.946816e-05 
    ## ... Procrustes: rmse 0.0001744706  max resid 0.0002982088 
    ## ... Similar to previous best
    ## Run 77 stress 9.500569e-05 
    ## ... Procrustes: rmse 0.0002045198  max resid 0.0003162171 
    ## ... Similar to previous best
    ## Run 78 stress 9.439972e-05 
    ## ... Procrustes: rmse 0.0001834166  max resid 0.0003240263 
    ## ... Similar to previous best
    ## Run 79 stress 0.2842805 
    ## Run 80 stress 0.1990774 
    ## Run 81 stress 9.807309e-05 
    ## ... Procrustes: rmse 0.0001890423  max resid 0.000334327 
    ## ... Similar to previous best
    ## Run 82 stress 0.001230345 
    ## Run 83 stress 0.001340815 
    ## Run 84 stress 0.1491957 
    ## Run 85 stress 0.1990774 
    ## Run 86 stress 0.0004834935 
    ## ... Procrustes: rmse 0.01613092  max resid 0.02221958 
    ## Run 87 stress 9.063282e-05 
    ## ... Procrustes: rmse 0.000145783  max resid 0.0002204516 
    ## ... Similar to previous best
    ## Run 88 stress 0.001368688 
    ## Run 89 stress 0.0006543723 
    ## Run 90 stress 9.179492e-05 
    ## ... Procrustes: rmse 0.0001314796  max resid 0.0001832433 
    ## ... Similar to previous best
    ## Run 91 stress 0.001322711 
    ## Run 92 stress 8.96682e-05 
    ## ... Procrustes: rmse 0.0001787197  max resid 0.0003273867 
    ## ... Similar to previous best
    ## Run 93 stress 0.001304877 
    ## Run 94 stress 0.0004041158 
    ## ... Procrustes: rmse 0.01475008  max resid 0.02031129 
    ## Run 95 stress 0.1776015 
    ## Run 96 stress 0.001177468 
    ## Run 97 stress 8.786935e-05 
    ## ... Procrustes: rmse 0.008968777  max resid 0.01244639 
    ## Run 98 stress 0.0004722332 
    ## ... Procrustes: rmse 0.01591381  max resid 0.02191995 
    ## Run 99 stress 0.1990774 
    ## Run 100 stress 9.779324e-05 
    ## ... Procrustes: rmse 0.0001827197  max resid 0.0003260257 
    ## ... Similar to previous best
    ## Run 101 stress 0.001239816 
    ## Run 102 stress 6.649401e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001313259  max resid 0.0002430071 
    ## ... Similar to previous best
    ## Run 103 stress 9.33653e-05 
    ## ... Procrustes: rmse 0.000100389  max resid 0.0001821404 
    ## ... Similar to previous best
    ## Run 104 stress 0.001234201 
    ## Run 105 stress 0.1990774 
    ## Run 106 stress 9.645742e-05 
    ## ... Procrustes: rmse 0.0001031976  max resid 0.0001916334 
    ## ... Similar to previous best
    ## Run 107 stress 9.584316e-05 
    ## ... Procrustes: rmse 9.536833e-05  max resid 0.0001835062 
    ## ... Similar to previous best
    ## Run 108 stress 8.383819e-05 
    ## ... Procrustes: rmse 0.002961958  max resid 0.004012073 
    ## ... Similar to previous best
    ## Run 109 stress 9.881204e-05 
    ## ... Procrustes: rmse 0.0001007678  max resid 0.0001660484 
    ## ... Similar to previous best
    ## Run 110 stress 9.589269e-05 
    ## ... Procrustes: rmse 0.006141683  max resid 0.008395509 
    ## Run 111 stress 9.196027e-05 
    ## ... Procrustes: rmse 7.087646e-05  max resid 0.0001317206 
    ## ... Similar to previous best
    ## Run 112 stress 0.1491957 
    ## Run 113 stress 9.484953e-05 
    ## ... Procrustes: rmse 5.178129e-05  max resid 7.769181e-05 
    ## ... Similar to previous best
    ## Run 114 stress 9.347662e-05 
    ## ... Procrustes: rmse 0.0001862166  max resid 0.0002931669 
    ## ... Similar to previous best
    ## Run 115 stress 9.133616e-05 
    ## ... Procrustes: rmse 0.0001376124  max resid 0.0003009512 
    ## ... Similar to previous best
    ## Run 116 stress 8.242401e-05 
    ## ... Procrustes: rmse 0.0001542516  max resid 0.0002998301 
    ## ... Similar to previous best
    ## Run 117 stress 9.012172e-05 
    ## ... Procrustes: rmse 0.0001640771  max resid 0.0002947551 
    ## ... Similar to previous best
    ## Run 118 stress 0.001258267 
    ## Run 119 stress 0.0004212003 
    ## ... Procrustes: rmse 0.01499086  max resid 0.02058995 
    ## Run 120 stress 9.964124e-05 
    ## ... Procrustes: rmse 0.0001903332  max resid 0.0003403265 
    ## ... Similar to previous best
    ## Run 121 stress 9.397435e-05 
    ## ... Procrustes: rmse 0.0001416234  max resid 0.0002047295 
    ## ... Similar to previous best
    ## Run 122 stress 0.1990774 
    ## Run 123 stress 8.170685e-05 
    ## ... Procrustes: rmse 0.0001197581  max resid 0.0001784926 
    ## ... Similar to previous best
    ## Run 124 stress 0.2848523 
    ## Run 125 stress 9.387338e-05 
    ## ... Procrustes: rmse 7.382294e-05  max resid 0.0001384785 
    ## ... Similar to previous best
    ## Run 126 stress 9.524086e-05 
    ## ... Procrustes: rmse 0.0001738054  max resid 0.0003072054 
    ## ... Similar to previous best
    ## Run 127 stress 7.415674e-05 
    ## ... Procrustes: rmse 0.002108861  max resid 0.002769129 
    ## ... Similar to previous best
    ## Run 128 stress 9.943006e-05 
    ## ... Procrustes: rmse 0.0001676972  max resid 0.0003323332 
    ## ... Similar to previous best
    ## Run 129 stress 0.001092341 
    ## Run 130 stress 0.1990774 
    ## Run 131 stress 0.0001085017 
    ## ... Procrustes: rmse 0.007574496  max resid 0.0103556 
    ## Run 132 stress 0.0002714787 
    ## ... Procrustes: rmse 0.01202178  max resid 0.01649274 
    ## Run 133 stress 0.001333374 
    ## Run 134 stress 0.001301554 
    ## Run 135 stress 0.1491957 
    ## Run 136 stress 9.864e-05 
    ## ... Procrustes: rmse 0.0001453617  max resid 0.0003130356 
    ## ... Similar to previous best
    ## Run 137 stress 9.182696e-05 
    ## ... Procrustes: rmse 0.0001211148  max resid 0.0002120994 
    ## ... Similar to previous best
    ## Run 138 stress 9.573505e-05 
    ## ... Procrustes: rmse 0.0001580811  max resid 0.0003002107 
    ## ... Similar to previous best
    ## Run 139 stress 0.001297806 
    ## Run 140 stress 0.1990774 
    ## Run 141 stress 0.001375338 
    ## Run 142 stress 9.040981e-05 
    ## ... Procrustes: rmse 0.006558315  max resid 0.008927759 
    ## Run 143 stress 0.001355676 
    ## Run 144 stress 0.001003501 
    ## Run 145 stress 8.978982e-05 
    ## ... Procrustes: rmse 0.0001600548  max resid 0.0002908142 
    ## ... Similar to previous best
    ## Run 146 stress 0.1776011 
    ## Run 147 stress 0.001382228 
    ## Run 148 stress 0.001216918 
    ## Run 149 stress 0.1491957 
    ## Run 150 stress 8.113346e-05 
    ## ... Procrustes: rmse 9.240516e-05  max resid 0.0001872875 
    ## ... Similar to previous best
    ## Run 151 stress 0.2509018 
    ## Run 152 stress 0.1491957 
    ## Run 153 stress 0.001275584 
    ## Run 154 stress 9.429854e-05 
    ## ... Procrustes: rmse 0.0001038705  max resid 0.0001878508 
    ## ... Similar to previous best
    ## Run 155 stress 0.1776011 
    ## Run 156 stress 0.001347928 
    ## Run 157 stress 0.001357628 
    ## Run 158 stress 0.001382534 
    ## Run 159 stress 8.639605e-05 
    ## ... Procrustes: rmse 0.000115456  max resid 0.0001958896 
    ## ... Similar to previous best
    ## Run 160 stress 0.001346444 
    ## Run 161 stress 9.420789e-05 
    ## ... Procrustes: rmse 0.0001274782  max resid 0.00021629 
    ## ... Similar to previous best
    ## Run 162 stress 0.1491957 
    ## Run 163 stress 0.001252933 
    ## Run 164 stress 0.1776016 
    ## Run 165 stress 7.43299e-05 
    ## ... Procrustes: rmse 0.0001089245  max resid 0.0001570492 
    ## ... Similar to previous best
    ## Run 166 stress 0.001464798 
    ## Run 167 stress 0.001377248 
    ## Run 168 stress 0.001167021 
    ## Run 169 stress 0.001380295 
    ## Run 170 stress 0.000463226 
    ## ... Procrustes: rmse 0.01572074  max resid 0.02159759 
    ## Run 171 stress 9.02221e-05 
    ## ... Procrustes: rmse 0.0001306253  max resid 0.0001839038 
    ## ... Similar to previous best
    ## Run 172 stress 8.360454e-05 
    ## ... Procrustes: rmse 5.54103e-05  max resid 0.0001095337 
    ## ... Similar to previous best
    ## Run 173 stress 0.1491957 
    ## Run 174 stress 9.664786e-05 
    ## ... Procrustes: rmse 0.0001911897  max resid 0.0003151866 
    ## ... Similar to previous best
    ## Run 175 stress 0.1776011 
    ## Run 176 stress 9.480707e-05 
    ## ... Procrustes: rmse 0.0001770058  max resid 0.0003314401 
    ## ... Similar to previous best
    ## Run 177 stress 0.1491957 
    ## Run 178 stress 0.1491957 
    ## Run 179 stress 0.1491957 
    ## Run 180 stress 0.1491957 
    ## Run 181 stress 0.00034844 
    ## ... Procrustes: rmse 0.02254335  max resid 0.03098996 
    ## Run 182 stress 9.340335e-05 
    ## ... Procrustes: rmse 0.0001726946  max resid 0.0003084986 
    ## ... Similar to previous best
    ## Run 183 stress 0.0001717417 
    ## ... Procrustes: rmse 0.01581595  max resid 0.02168465 
    ## Run 184 stress 0.00118003 
    ## Run 185 stress 0.2578645 
    ## Run 186 stress 0.001283839 
    ## Run 187 stress 0.1491957 
    ## Run 188 stress 9.992161e-05 
    ## ... Procrustes: rmse 0.007262115  max resid 0.009924502 
    ## Run 189 stress 8.671843e-05 
    ## ... Procrustes: rmse 0.0001630437  max resid 0.0003110707 
    ## ... Similar to previous best
    ## Run 190 stress 0.001118226 
    ## Run 191 stress 0.0004594232 
    ## ... Procrustes: rmse 0.01562261  max resid 0.02146169 
    ## Run 192 stress 9.796932e-05 
    ## ... Procrustes: rmse 0.0001727795  max resid 0.0002938094 
    ## ... Similar to previous best
    ## Run 193 stress 8.886476e-05 
    ## ... Procrustes: rmse 0.0001199412  max resid 0.0002011812 
    ## ... Similar to previous best
    ## Run 194 stress 0.2832379 
    ## Run 195 stress 0.001391191 
    ## Run 196 stress 0.001423196 
    ## Run 197 stress 0.1776013 
    ## Run 198 stress 0.001262521 
    ## Run 199 stress 9.050155e-05 
    ## ... Procrustes: rmse 9.824154e-05  max resid 0.0001739124 
    ## ... Similar to previous best
    ## Run 200 stress 0.001177874 
    ## Run 201 stress 0.1491957 
    ## Run 202 stress 0.0005082877 
    ## ... Procrustes: rmse 0.01647325  max resid 0.02263535 
    ## Run 203 stress 8.969957e-05 
    ## ... Procrustes: rmse 0.0001207682  max resid 0.000202312 
    ## ... Similar to previous best
    ## Run 204 stress 0.001152659 
    ## Run 205 stress 0.001276267 
    ## Run 206 stress 0.1990774 
    ## Run 207 stress 0.001378705 
    ## Run 208 stress 9.460987e-05 
    ## ... Procrustes: rmse 0.0001024203  max resid 0.0001747463 
    ## ... Similar to previous best
    ## Run 209 stress 9.488712e-05 
    ## ... Procrustes: rmse 0.0001384414  max resid 0.0002077426 
    ## ... Similar to previous best
    ## Run 210 stress 0.0009880706 
    ## Run 211 stress 0.279732 
    ## Run 212 stress 0.001436901 
    ## Run 213 stress 0.1990774 
    ## Run 214 stress 9.43629e-05 
    ## ... Procrustes: rmse 0.0001654501  max resid 0.000265091 
    ## ... Similar to previous best
    ## Run 215 stress 0.3078457 
    ## Run 216 stress 9.103636e-05 
    ## ... Procrustes: rmse 9.637748e-05  max resid 0.0001681418 
    ## ... Similar to previous best
    ## Run 217 stress 0.2373687 
    ## Run 218 stress 0.0007080171 
    ## Run 219 stress 8.485487e-05 
    ## ... Procrustes: rmse 0.0001544438  max resid 0.0003002335 
    ## ... Similar to previous best
    ## Run 220 stress 0.001250122 
    ## Run 221 stress 0.0004849539 
    ## ... Procrustes: rmse 0.01608955  max resid 0.02210583 
    ## Run 222 stress 0.0004066847 
    ## ... Procrustes: rmse 0.01469975  max resid 0.02019008 
    ## Run 223 stress 0.2494706 
    ## Run 224 stress 0.001257817 
    ## Run 225 stress 0.001166407 
    ## Run 226 stress 9.156092e-05 
    ## ... Procrustes: rmse 9.385651e-05  max resid 0.0001580676 
    ## ... Similar to previous best
    ## Run 227 stress 0.0004663935 
    ## ... Procrustes: rmse 0.01510591  max resid 0.02075085 
    ## Run 228 stress 9.881691e-05 
    ## ... Procrustes: rmse 9.119966e-05  max resid 0.0001511024 
    ## ... Similar to previous best
    ## Run 229 stress 0.177602 
    ## Run 230 stress 0.1491957 
    ## Run 231 stress 0.0004668326 
    ## ... Procrustes: rmse 0.01578347  max resid 0.02168479 
    ## Run 232 stress 9.549408e-05 
    ## ... Procrustes: rmse 0.0001695028  max resid 0.0003019798 
    ## ... Similar to previous best
    ## Run 233 stress 0.001172007 
    ## Run 234 stress 0.001350475 
    ## Run 235 stress 9.834607e-05 
    ## ... Procrustes: rmse 0.0001532073  max resid 0.000216581 
    ## ... Similar to previous best
    ## Run 236 stress 9.922895e-05 
    ## ... Procrustes: rmse 7.628338e-05  max resid 0.0001423102 
    ## ... Similar to previous best
    ## Run 237 stress 0.0004364579 
    ## ... Procrustes: rmse 0.01526072  max resid 0.02096066 
    ## Run 238 stress 9.15273e-05 
    ## ... Procrustes: rmse 0.000169125  max resid 0.0003180213 
    ## ... Similar to previous best
    ## Run 239 stress 0.001128899 
    ## Run 240 stress 9.952985e-05 
    ## ... Procrustes: rmse 0.003550553  max resid 0.004827633 
    ## ... Similar to previous best
    ## Run 241 stress 8.771599e-05 
    ## ... Procrustes: rmse 8.993655e-05  max resid 0.0001606389 
    ## ... Similar to previous best
    ## Run 242 stress 0.0003446015 
    ## ... Procrustes: rmse 0.02241601  max resid 0.03081364 
    ## Run 243 stress 0.0004393279 
    ## ... Procrustes: rmse 0.01531123  max resid 0.02103136 
    ## Run 244 stress 5.840364e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.003418943  max resid 0.00462807 
    ## ... Similar to previous best
    ## Run 245 stress 0.2520602 
    ## Run 246 stress 0.001199635 
    ## Run 247 stress 9.768546e-05 
    ## ... Procrustes: rmse 0.003400442  max resid 0.004506479 
    ## ... Similar to previous best
    ## Run 248 stress 0.001121148 
    ## Run 249 stress 0.0004864001 
    ## ... Procrustes: rmse 0.01268778  max resid 0.0175048 
    ## Run 250 stress 0.2842805 
    ## Run 251 stress 0.001328013 
    ## Run 252 stress 9.992985e-05 
    ## ... Procrustes: rmse 0.003401269  max resid 0.00450704 
    ## ... Similar to previous best
    ## Run 253 stress 0.001158336 
    ## Run 254 stress 0.1491957 
    ## Run 255 stress 0.2528294 
    ## Run 256 stress 0.2854321 
    ## Run 257 stress 0.1491957 
    ## Run 258 stress 0.0004040971 
    ## ... Procrustes: rmse 0.01126491  max resid 0.01553811 
    ## Run 259 stress 9.696135e-05 
    ## ... Procrustes: rmse 0.003383341  max resid 0.004515251 
    ## ... Similar to previous best
    ## Run 260 stress 8.906106e-05 
    ## ... Procrustes: rmse 0.003414717  max resid 0.004530509 
    ## ... Similar to previous best
    ## Run 261 stress 0.001403267 
    ## Run 262 stress 0.2440272 
    ## Run 263 stress 9.952274e-05 
    ## ... Procrustes: rmse 0.003474968  max resid 0.004527818 
    ## ... Similar to previous best
    ## Run 264 stress 8.078684e-05 
    ## ... Procrustes: rmse 0.003491448  max resid 0.004849259 
    ## ... Similar to previous best
    ## Run 265 stress 9.244912e-05 
    ## ... Procrustes: rmse 0.003392093  max resid 0.004545665 
    ## ... Similar to previous best
    ## Run 266 stress 8.457577e-05 
    ## ... Procrustes: rmse 0.003414585  max resid 0.004545494 
    ## ... Similar to previous best
    ## Run 267 stress 0.001294636 
    ## Run 268 stress 0.001336642 
    ## Run 269 stress 0.001201551 
    ## Run 270 stress 0.001227891 
    ## Run 271 stress 0.2832508 
    ## Run 272 stress 0.001097792 
    ## Run 273 stress 0.2570108 
    ## Run 274 stress 0.0004859638 
    ## ... Procrustes: rmse 0.01269144  max resid 0.01750419 
    ## Run 275 stress 7.925252e-05 
    ## ... Procrustes: rmse 0.003427601  max resid 0.004578765 
    ## ... Similar to previous best
    ## Run 276 stress 9.278376e-05 
    ## ... Procrustes: rmse 0.003419814  max resid 0.004552865 
    ## ... Similar to previous best
    ## Run 277 stress 0.001157994 
    ## Run 278 stress 9.703638e-05 
    ## ... Procrustes: rmse 0.003382024  max resid 0.004510377 
    ## ... Similar to previous best
    ## Run 279 stress 0.1491957 
    ## Run 280 stress 0.0001425403 
    ## ... Procrustes: rmse 0.005272206  max resid 0.007267696 
    ## Run 281 stress 0.0008755862 
    ## Run 282 stress 9.182368e-05 
    ## ... Procrustes: rmse 0.003473266  max resid 0.004811277 
    ## ... Similar to previous best
    ## Run 283 stress 9.124555e-05 
    ## ... Procrustes: rmse 0.003415052  max resid 0.00452695 
    ## ... Similar to previous best
    ## Run 284 stress 0.001277906 
    ## Run 285 stress 0.0003956299 
    ## ... Procrustes: rmse 0.01106036  max resid 0.01525424 
    ## Run 286 stress 0.001329371 
    ## Run 287 stress 0.0004625656 
    ## ... Procrustes: rmse 0.01229718  max resid 0.01696032 
    ## Run 288 stress 0.2848517 
    ## Run 289 stress 8.340703e-05 
    ## ... Procrustes: rmse 0.003464187  max resid 0.004563804 
    ## ... Similar to previous best
    ## Run 290 stress 0.001389863 
    ## Run 291 stress 0.001226472 
    ## Run 292 stress 0.1491957 
    ## Run 293 stress 9.007904e-05 
    ## ... Procrustes: rmse 0.003481748  max resid 0.004566361 
    ## ... Similar to previous best
    ## Run 294 stress 9.781806e-05 
    ## ... Procrustes: rmse 0.003383461  max resid 0.004514103 
    ## ... Similar to previous best
    ## Run 295 stress 0.0004765913 
    ## ... Procrustes: rmse 0.01252984  max resid 0.01728134 
    ## Run 296 stress 0.001327021 
    ## Run 297 stress 0.1491957 
    ## Run 298 stress 0.001253087 
    ## Run 299 stress 0.001407867 
    ## Run 300 stress 8.579255e-05 
    ## ... Procrustes: rmse 0.003423808  max resid 0.004561905 
    ## ... Similar to previous best
    ## Run 301 stress 0.0004886381 
    ## ... Procrustes: rmse 0.01271061  max resid 0.0175308 
    ## Run 302 stress 8.969574e-05 
    ## ... Procrustes: rmse 0.003485033  max resid 0.004561534 
    ## ... Similar to previous best
    ## Run 303 stress 0.1776016 
    ## Run 304 stress 9.651082e-05 
    ## ... Procrustes: rmse 0.003404521  max resid 0.004524661 
    ## ... Similar to previous best
    ## Run 305 stress 9.776969e-05 
    ## ... Procrustes: rmse 0.003480853  max resid 0.004589624 
    ## ... Similar to previous best
    ## Run 306 stress 0.001348691 
    ## Run 307 stress 6.74526e-05 
    ## ... Procrustes: rmse 0.003174947  max resid 0.00433381 
    ## ... Similar to previous best
    ## Run 308 stress 0.1491957 
    ## Run 309 stress 9.033006e-05 
    ## ... Procrustes: rmse 0.003443077  max resid 0.004534518 
    ## ... Similar to previous best
    ## Run 310 stress 0.1990774 
    ## Run 311 stress 0.001419843 
    ## Run 312 stress 9.409469e-05 
    ## ... Procrustes: rmse 0.003460664  max resid 0.004577618 
    ## ... Similar to previous best
    ## Run 313 stress 9.678229e-05 
    ## ... Procrustes: rmse 0.003456988  max resid 0.004563035 
    ## ... Similar to previous best
    ## Run 314 stress 8.962518e-05 
    ## ... Procrustes: rmse 0.003381064  max resid 0.004539372 
    ## ... Similar to previous best
    ## Run 315 stress 0.001248478 
    ## Run 316 stress 0.001360572 
    ## Run 317 stress 0.0001121717 
    ## ... Procrustes: rmse 0.004279294  max resid 0.005897786 
    ## ... Similar to previous best
    ## Run 318 stress 9.885972e-05 
    ## ... Procrustes: rmse 0.003404959  max resid 0.004555833 
    ## ... Similar to previous best
    ## Run 319 stress 9.26782e-05 
    ## ... Procrustes: rmse 0.003416183  max resid 0.00454323 
    ## ... Similar to previous best
    ## Run 320 stress 0.001164458 
    ## Run 321 stress 9.119529e-05 
    ## ... Procrustes: rmse 0.003411507  max resid 0.004531624 
    ## ... Similar to previous best
    ## Run 322 stress 0.1990774 
    ## Run 323 stress 0.0004743917 
    ## ... Procrustes: rmse 0.012497  max resid 0.01723646 
    ## Run 324 stress 9.175603e-05 
    ## ... Procrustes: rmse 0.003404521  max resid 0.004520714 
    ## ... Similar to previous best
    ## Run 325 stress 0.0004629244 
    ## ... Procrustes: rmse 0.01228372  max resid 0.01694319 
    ## Run 326 stress 9.286646e-05 
    ## ... Procrustes: rmse 0.003448515  max resid 0.004528709 
    ## ... Similar to previous best
    ## Run 327 stress 0.3083098 
    ## Run 328 stress 0.1491957 
    ## Run 329 stress 0.001335258 
    ## Run 330 stress 0.0004893843 
    ## ... Procrustes: rmse 0.01274612  max resid 0.01757981 
    ## Run 331 stress 0.0009553928 
    ## Run 332 stress 0.001125783 
    ## Run 333 stress 0.001341184 
    ## Run 334 stress 9.151568e-05 
    ## ... Procrustes: rmse 0.001272752  max resid 0.001768215 
    ## ... Similar to previous best
    ## Run 335 stress 9.414688e-05 
    ## ... Procrustes: rmse 0.003473102  max resid 0.004807301 
    ## ... Similar to previous best
    ## Run 336 stress 9.364141e-05 
    ## ... Procrustes: rmse 0.003481449  max resid 0.00456443 
    ## ... Similar to previous best
    ## Run 337 stress 0.1776011 
    ## Run 338 stress 9.54936e-05 
    ## ... Procrustes: rmse 0.003455663  max resid 0.004570737 
    ## ... Similar to previous best
    ## Run 339 stress 9.592547e-05 
    ## ... Procrustes: rmse 0.003481698  max resid 0.004583063 
    ## ... Similar to previous best
    ## Run 340 stress 0.001310871 
    ## Run 341 stress 0.001268069 
    ## Run 342 stress 0.1491957 
    ## Run 343 stress 0.001077433 
    ## Run 344 stress 0.001174573 
    ## Run 345 stress 9.912593e-05 
    ## ... Procrustes: rmse 0.003404848  max resid 0.004560263 
    ## ... Similar to previous best
    ## Run 346 stress 0.00100894 
    ## Run 347 stress 0.001351895 
    ## Run 348 stress 0.0004752652 
    ## ... Procrustes: rmse 0.0125117  max resid 0.01725634 
    ## Run 349 stress 0.0005035827 
    ## ... Procrustes: rmse 0.01297998  max resid 0.01790235 
    ## Run 350 stress 9.680243e-05 
    ## ... Procrustes: rmse 0.0034667  max resid 0.004799646 
    ## ... Similar to previous best
    ## Run 351 stress 0.1491957 
    ## Run 352 stress 8.939472e-05 
    ## ... Procrustes: rmse 0.003423658  max resid 0.004555561 
    ## ... Similar to previous best
    ## Run 353 stress 0.000494181 
    ## ... Procrustes: rmse 0.0128269  max resid 0.0176925 
    ## Run 354 stress 0.2570107 
    ## Run 355 stress 0.2494706 
    ## Run 356 stress 0.001364449 
    ## Run 357 stress 0.001180684 
    ## Run 358 stress 0.1990774 
    ## Run 359 stress 0.1776015 
    ## Run 360 stress 0.1776011 
    ## Run 361 stress 0.001231967 
    ## Run 362 stress 0.1776011 
    ## Run 363 stress 0.1491957 
    ## Run 364 stress 0.001343166 
    ## Run 365 stress 9.696406e-05 
    ## ... Procrustes: rmse 0.003369879  max resid 0.004506879 
    ## ... Similar to previous best
    ## Run 366 stress 0.0005997646 
    ## Run 367 stress 9.387683e-05 
    ## ... Procrustes: rmse 0.00343623  max resid 0.004526195 
    ## ... Similar to previous best
    ## Run 368 stress 0.00131442 
    ## Run 369 stress 0.001367169 
    ## Run 370 stress 0.001245359 
    ## Run 371 stress 0.001382202 
    ## Run 372 stress 8.520424e-05 
    ## ... Procrustes: rmse 0.0034811  max resid 0.004556129 
    ## ... Similar to previous best
    ## Run 373 stress 0.0004809649 
    ## ... Procrustes: rmse 0.01260306  max resid 0.01738254 
    ## Run 374 stress 8.176958e-05 
    ## ... Procrustes: rmse 0.003425041  max resid 0.004624094 
    ## ... Similar to previous best
    ## Run 375 stress 0.2578638 
    ## Run 376 stress 9.449507e-05 
    ## ... Procrustes: rmse 0.003396083  max resid 0.00452824 
    ## ... Similar to previous best
    ## Run 377 stress 0.001202962 
    ## Run 378 stress 8.830909e-05 
    ## ... Procrustes: rmse 0.00348718  max resid 0.004845781 
    ## ... Similar to previous best
    ## Run 379 stress 0.001115112 
    ## Run 380 stress 9.219679e-05 
    ## ... Procrustes: rmse 0.003413288  max resid 0.00453398 
    ## ... Similar to previous best
    ## Run 381 stress 0.2852153 
    ## Run 382 stress 0.0006067781 
    ## Run 383 stress 9.906385e-05 
    ## ... Procrustes: rmse 0.0034296  max resid 0.004737051 
    ## ... Similar to previous best
    ## Run 384 stress 0.000313672 
    ## ... Procrustes: rmse 0.009504694  max resid 0.0131076 
    ## Run 385 stress 0.0004707437 
    ## ... Procrustes: rmse 0.0124321  max resid 0.01714663 
    ## Run 386 stress 0.0004534543 
    ## ... Procrustes: rmse 0.01214101  max resid 0.01674322 
    ## Run 387 stress 9.258802e-05 
    ## ... Procrustes: rmse 0.003478697  max resid 0.004806255 
    ## ... Similar to previous best
    ## Run 388 stress 9.846305e-05 
    ## ... Procrustes: rmse 0.003421242  max resid 0.004563784 
    ## ... Similar to previous best
    ## Run 389 stress 0.001316233 
    ## Run 390 stress 0.1491957 
    ## Run 391 stress 9.871803e-05 
    ## ... Procrustes: rmse 0.003458576  max resid 0.004582866 
    ## ... Similar to previous best
    ## Run 392 stress 0.001362894 
    ## Run 393 stress 9.383779e-05 
    ## ... Procrustes: rmse 0.003395583  max resid 0.004556399 
    ## ... Similar to previous best
    ## Run 394 stress 9.715648e-05 
    ## ... Procrustes: rmse 0.003447353  max resid 0.004538189 
    ## ... Similar to previous best
    ## Run 395 stress 0.001287087 
    ## Run 396 stress 9.327996e-05 
    ## ... Procrustes: rmse 0.003407776  max resid 0.004533022 
    ## ... Similar to previous best
    ## Run 397 stress 9.361182e-05 
    ## ... Procrustes: rmse 0.003481057  max resid 0.004575638 
    ## ... Similar to previous best
    ## Run 398 stress 0.0004814238 
    ## ... Procrustes: rmse 0.01261573  max resid 0.01739986 
    ## Run 399 stress 0.1491957 
    ## Run 400 stress 8.458759e-05 
    ## ... Procrustes: rmse 0.003469269  max resid 0.004805433 
    ## ... Similar to previous best
    ## Run 401 stress 9.359293e-05 
    ## ... Procrustes: rmse 0.003420912  max resid 0.00454718 
    ## ... Similar to previous best
    ## Run 402 stress 9.604538e-05 
    ## ... Procrustes: rmse 0.002183149  max resid 0.003025949 
    ## ... Similar to previous best
    ## Run 403 stress 9.531427e-05 
    ## ... Procrustes: rmse 0.00335416  max resid 0.004531805 
    ## ... Similar to previous best
    ## Run 404 stress 8.406035e-05 
    ## ... Procrustes: rmse 0.003430372  max resid 0.004556493 
    ## ... Similar to previous best
    ## Run 405 stress 0.0005156612 
    ## ... Procrustes: rmse 0.01264333  max resid 0.01744814 
    ## Run 406 stress 0.0004623202 
    ## ... Procrustes: rmse 0.01229211  max resid 0.01695327 
    ## Run 407 stress 9.60305e-05 
    ## ... Procrustes: rmse 0.003403976  max resid 0.004533351 
    ## ... Similar to previous best
    ## Run 408 stress 9.25782e-05 
    ## ... Procrustes: rmse 0.003412192  max resid 0.00452788 
    ## ... Similar to previous best
    ## Run 409 stress 0.1990774 
    ## Run 410 stress 9.772447e-05 
    ## ... Procrustes: rmse 0.003402836  max resid 0.004541337 
    ## ... Similar to previous best
    ## Run 411 stress 9.445385e-05 
    ## ... Procrustes: rmse 0.003405722  max resid 0.004509566 
    ## ... Similar to previous best
    ## Run 412 stress 0.1491957 
    ## Run 413 stress 0.0002605912 
    ## ... Procrustes: rmse 0.01882804  max resid 0.03177227 
    ## Run 414 stress 0.1776013 
    ## Run 415 stress 7.92618e-05 
    ## ... Procrustes: rmse 0.00342861  max resid 0.005502367 
    ## ... Similar to previous best
    ## Run 416 stress 8.76919e-05 
    ## ... Procrustes: rmse 0.003398173  max resid 0.006456762 
    ## ... Similar to previous best
    ## Run 417 stress 9.249667e-05 
    ## ... Procrustes: rmse 0.003406755  max resid 0.00452637 
    ## ... Similar to previous best
    ## Run 418 stress 9.326975e-05 
    ## ... Procrustes: rmse 0.005798057  max resid 0.01274596 
    ## Run 419 stress 0.001099393 
    ## Run 420 stress 0.0005084117 
    ## ... Procrustes: rmse 0.01305752  max resid 0.0180094 
    ## Run 421 stress 0.0004305885 
    ## ... Procrustes: rmse 0.01173953  max resid 0.01619075 
    ## Run 422 stress 0.1776013 
    ## Run 423 stress 0.0004496901 
    ## ... Procrustes: rmse 0.01206649  max resid 0.01664462 
    ## Run 424 stress 9.59778e-05 
    ## ... Procrustes: rmse 0.003454793  max resid 0.005278501 
    ## ... Similar to previous best
    ## Run 425 stress 0.1776019 
    ## Run 426 stress 0.1491957 
    ## Run 427 stress 0.000232136 
    ## ... Procrustes: rmse 0.01774994  max resid 0.03025725 
    ## Run 428 stress 0.001350937 
    ## Run 429 stress 9.152048e-05 
    ## ... Procrustes: rmse 0.00354123  max resid 0.007641206 
    ## ... Similar to previous best
    ## Run 430 stress 0.1990774 
    ## Run 431 stress 9.191622e-05 
    ## ... Procrustes: rmse 0.003453834  max resid 0.004601011 
    ## ... Similar to previous best
    ## Run 432 stress 7.979362e-05 
    ## ... Procrustes: rmse 0.003490359  max resid 0.004834748 
    ## ... Similar to previous best
    ## Run 433 stress 0.2354286 
    ## Run 434 stress 0.0004764311 
    ## ... Procrustes: rmse 0.01251795  max resid 0.017266 
    ## Run 435 stress 0.001459034 
    ## Run 436 stress 9.712643e-05 
    ## ... Procrustes: rmse 0.003444013  max resid 0.00452351 
    ## ... Similar to previous best
    ## Run 437 stress 6.29307e-05 
    ## ... Procrustes: rmse 0.003408595  max resid 0.005656717 
    ## ... Similar to previous best
    ## Run 438 stress 0.0001531957 
    ## ... Procrustes: rmse 0.01436953  max resid 0.02548263 
    ## Run 439 stress 0.001349364 
    ## Run 440 stress 0.001411292 
    ## Run 441 stress 0.1491957 
    ## Run 442 stress 0.0001045353 
    ## ... Procrustes: rmse 0.00401551  max resid 0.005533804 
    ## ... Similar to previous best
    ## Run 443 stress 9.859907e-05 
    ## ... Procrustes: rmse 0.003410684  max resid 0.004520364 
    ## ... Similar to previous best
    ## Run 444 stress 0.001347441 
    ## Run 445 stress 0.001184171 
    ## Run 446 stress 0.001285716 
    ## Run 447 stress 0.001415532 
    ## Run 448 stress 0.00114869 
    ## Run 449 stress 0.2570118 
    ## Run 450 stress 0.1491957 
    ## Run 451 stress 9.154147e-05 
    ## ... Procrustes: rmse 0.003482119  max resid 0.004564771 
    ## ... Similar to previous best
    ## Run 452 stress 0.00113685 
    ## Run 453 stress 0.2842805 
    ## Run 454 stress 9.421925e-05 
    ## ... Procrustes: rmse 0.00309969  max resid 0.004231327 
    ## ... Similar to previous best
    ## Run 455 stress 0.0004989715 
    ## ... Procrustes: rmse 0.01290676  max resid 0.01780127 
    ## Run 456 stress 0.0003703675 
    ## ... Procrustes: rmse 0.01063703  max resid 0.01467013 
    ## Run 457 stress 0.0004877096 
    ## ... Procrustes: rmse 0.01261414  max resid 0.01739998 
    ## Run 458 stress 0.0009131903 
    ## Run 459 stress 8.276448e-05 
    ## ... Procrustes: rmse 0.003482728  max resid 0.004549304 
    ## ... Similar to previous best
    ## Run 460 stress 0.001296095 
    ## Run 461 stress 8.246072e-05 
    ## ... Procrustes: rmse 0.001369807  max resid 0.001823326 
    ## ... Similar to previous best
    ## Run 462 stress 0.00111481 
    ## Run 463 stress 8.538748e-05 
    ## ... Procrustes: rmse 0.003438712  max resid 0.007070971 
    ## ... Similar to previous best
    ## Run 464 stress 8.00355e-05 
    ## ... Procrustes: rmse 0.00307334  max resid 0.004238671 
    ## ... Similar to previous best
    ## Run 465 stress 9.997021e-05 
    ## ... Procrustes: rmse 0.003849746  max resid 0.005305086 
    ## ... Similar to previous best
    ## Run 466 stress 7.875987e-05 
    ## ... Procrustes: rmse 0.003470239  max resid 0.004816115 
    ## ... Similar to previous best
    ## Run 467 stress 9.560011e-05 
    ## ... Procrustes: rmse 0.003410992  max resid 0.005795266 
    ## ... Similar to previous best
    ## Run 468 stress 0.001180027 
    ## Run 469 stress 0.0004692891 
    ## ... Procrustes: rmse 0.01241146  max resid 0.01711819 
    ## Run 470 stress 0.1776012 
    ## Run 471 stress 0.001277149 
    ## Run 472 stress 0.001479224 
    ## Run 473 stress 0.0005105901 
    ## ... Procrustes: rmse 0.01299198  max resid 0.01791645 
    ## Run 474 stress 0.001519924 
    ## Run 475 stress 0.001366012 
    ## Run 476 stress 0.001389493 
    ## Run 477 stress 0.1491957 
    ## Run 478 stress 0.001235923 
    ## Run 479 stress 9.601122e-05 
    ## ... Procrustes: rmse 0.003444177  max resid 0.005373667 
    ## ... Similar to previous best
    ## Run 480 stress 0.0004170529 
    ## ... Procrustes: rmse 0.01149866  max resid 0.01585942 
    ## Run 481 stress 7.485755e-05 
    ## ... Procrustes: rmse 0.003437073  max resid 0.005442215 
    ## ... Similar to previous best
    ## Run 482 stress 0.001282705 
    ## Run 483 stress 9.366745e-05 
    ## ... Procrustes: rmse 0.003458089  max resid 0.004574669 
    ## ... Similar to previous best
    ## Run 484 stress 0.2696744 
    ## Run 485 stress 9.404163e-05 
    ## ... Procrustes: rmse 0.003374618  max resid 0.004520749 
    ## ... Similar to previous best
    ## Run 486 stress 0.001163848 
    ## Run 487 stress 0.001292274 
    ## Run 488 stress 0.0004832686 
    ## ... Procrustes: rmse 0.01264642  max resid 0.01744219 
    ## Run 489 stress 0.001292837 
    ## Run 490 stress 0.001239097 
    ## Run 491 stress 5.122197e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.003483842  max resid 0.004701384 
    ## ... Similar to previous best
    ## Run 492 stress 0.001271576 
    ## Run 493 stress 0.1776022 
    ## Run 494 stress 5.80808e-05 
    ## ... Procrustes: rmse 2.477649e-05  max resid 4.345275e-05 
    ## ... Similar to previous best
    ## Run 495 stress 0.001288466 
    ## Run 496 stress 0.001384324 
    ## Run 497 stress 0.0001103221 
    ## ... Procrustes: rmse 0.007702737  max resid 0.01051706 
    ## Run 498 stress 8.602824e-05 
    ## ... Procrustes: rmse 0.007040121  max resid 0.009756202 
    ## Run 499 stress 9.692333e-05 
    ## ... Procrustes: rmse 0.0001439195  max resid 0.0002068751 
    ## ... Similar to previous best
    ## Run 500 stress 0.001315705 
    ## *** Best solution repeated 3 times

    ## Warning in metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
SD_beta_geo_NMDS <- metaMDS(SD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.1067587 
    ## Run 2 stress 0.09099556 
    ## Run 3 stress 0.09039143 
    ## ... Procrustes: rmse 0.01220161  max resid 0.0329574 
    ## Run 4 stress 0.1071321 
    ## Run 5 stress 0.09088605 
    ## Run 6 stress 0.1052648 
    ## Run 7 stress 0.08938961 
    ## ... New best solution
    ## ... Procrustes: rmse 0.02649993  max resid 0.07738992 
    ## Run 8 stress 0.08938544 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0372362  max resid 0.1196088 
    ## Run 9 stress 0.107902 
    ## Run 10 stress 0.09503435 
    ## Run 11 stress 0.09775551 
    ## Run 12 stress 0.08946338 
    ## ... Procrustes: rmse 0.03575926  max resid 0.1200205 
    ## Run 13 stress 0.1079022 
    ## Run 14 stress 0.1071321 
    ## Run 15 stress 0.08938544 
    ## ... Procrustes: rmse 0.0002179929  max resid 0.0006762695 
    ## ... Similar to previous best
    ## Run 16 stress 0.0893854 
    ## ... New best solution
    ## ... Procrustes: rmse 5.324299e-05  max resid 0.0001789261 
    ## ... Similar to previous best
    ## Run 17 stress 0.1108635 
    ## Run 18 stress 0.09503437 
    ## Run 19 stress 0.1052649 
    ## Run 20 stress 0.0892608 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01182969  max resid 0.04010147 
    ## Run 21 stress 0.1110886 
    ## Run 22 stress 0.08938541 
    ## ... Procrustes: rmse 0.01183798  max resid 0.03982764 
    ## Run 23 stress 0.1076304 
    ## Run 24 stress 0.08926068 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000115916  max resid 0.0003399188 
    ## ... Similar to previous best
    ## Run 25 stress 0.1091826 
    ## Run 26 stress 0.08938538 
    ## ... Procrustes: rmse 0.01181358  max resid 0.03973732 
    ## Run 27 stress 0.1075421 
    ## Run 28 stress 0.1074432 
    ## Run 29 stress 0.09612797 
    ## Run 30 stress 0.1064188 
    ## Run 31 stress 0.08938962 
    ## ... Procrustes: rmse 0.03598964  max resid 0.1183458 
    ## Run 32 stress 0.09087354 
    ## Run 33 stress 0.09044614 
    ## Run 34 stress 0.1067585 
    ## Run 35 stress 0.08946656 
    ## ... Procrustes: rmse 0.03322298  max resid 0.1163433 
    ## Run 36 stress 0.08926072 
    ## ... Procrustes: rmse 3.985556e-05  max resid 0.0001287386 
    ## ... Similar to previous best
    ## Run 37 stress 0.1063194 
    ## Run 38 stress 0.1111367 
    ## Run 39 stress 0.09262334 
    ## Run 40 stress 0.09018708 
    ## Run 41 stress 0.09088615 
    ## Run 42 stress 0.108757 
    ## Run 43 stress 0.1091452 
    ## Run 44 stress 0.08938542 
    ## ... Procrustes: rmse 0.01179064  max resid 0.03974609 
    ## Run 45 stress 0.08938961 
    ## ... Procrustes: rmse 0.03597324  max resid 0.1183184 
    ## Run 46 stress 0.106535 
    ## Run 47 stress 0.09021166 
    ## Run 48 stress 0.1056908 
    ## Run 49 stress 0.1105341 
    ## Run 50 stress 0.08938541 
    ## ... Procrustes: rmse 0.01184223  max resid 0.03976891 
    ## Run 51 stress 0.10626 
    ## Run 52 stress 0.08926072 
    ## ... Procrustes: rmse 3.762975e-05  max resid 0.0001124333 
    ## ... Similar to previous best
    ## Run 53 stress 0.1079014 
    ## Run 54 stress 0.1101449 
    ## Run 55 stress 0.08946652 
    ## ... Procrustes: rmse 0.03323482  max resid 0.1163634 
    ## Run 56 stress 0.08946328 
    ## ... Procrustes: rmse 0.03754385  max resid 0.1178204 
    ## Run 57 stress 0.08946651 
    ## ... Procrustes: rmse 0.03325105  max resid 0.1163731 
    ## Run 58 stress 0.08938972 
    ## ... Procrustes: rmse 0.03594188  max resid 0.1182706 
    ## Run 59 stress 0.08938973 
    ## ... Procrustes: rmse 0.03593996  max resid 0.1182675 
    ## Run 60 stress 0.08926082 
    ## ... Procrustes: rmse 0.0001286467  max resid 0.000419274 
    ## ... Similar to previous best
    ## Run 61 stress 0.08938541 
    ## ... Procrustes: rmse 0.0118512  max resid 0.03973974 
    ## Run 62 stress 0.0903913 
    ## Run 63 stress 0.09087342 
    ## Run 64 stress 0.09503416 
    ## Run 65 stress 0.08938544 
    ## ... Procrustes: rmse 0.01180542  max resid 0.03975504 
    ## Run 66 stress 0.08926064 
    ## ... New best solution
    ## ... Procrustes: rmse 7.641294e-05  max resid 0.0002382811 
    ## ... Similar to previous best
    ## Run 67 stress 0.08926066 
    ## ... Procrustes: rmse 0.0001039359  max resid 0.0003291034 
    ## ... Similar to previous best
    ## Run 68 stress 0.08926066 
    ## ... Procrustes: rmse 3.829983e-05  max resid 9.313366e-05 
    ## ... Similar to previous best
    ## Run 69 stress 0.09044623 
    ## Run 70 stress 0.09180889 
    ## Run 71 stress 0.1060095 
    ## Run 72 stress 0.1052647 
    ## Run 73 stress 0.08938546 
    ## ... Procrustes: rmse 0.01187976  max resid 0.03972565 
    ## Run 74 stress 0.0895172 
    ## ... Procrustes: rmse 0.03505314  max resid 0.1156364 
    ## Run 75 stress 0.08926072 
    ## ... Procrustes: rmse 0.0001140444  max resid 0.0003556047 
    ## ... Similar to previous best
    ## Run 76 stress 0.09178294 
    ## Run 77 stress 0.1092098 
    ## Run 78 stress 0.08938966 
    ## ... Procrustes: rmse 0.03594036  max resid 0.1182704 
    ## Run 79 stress 0.1075423 
    ## Run 80 stress 0.09044617 
    ## Run 81 stress 0.1067588 
    ## Run 82 stress 0.1064428 
    ## Run 83 stress 0.09503414 
    ## Run 84 stress 0.1063188 
    ## Run 85 stress 0.1091241 
    ## Run 86 stress 0.08938965 
    ## ... Procrustes: rmse 0.0359421  max resid 0.1182709 
    ## Run 87 stress 0.09021171 
    ## Run 88 stress 0.1076304 
    ## Run 89 stress 0.1056904 
    ## Run 90 stress 0.08946654 
    ## ... Procrustes: rmse 0.03321308  max resid 0.1163331 
    ## Run 91 stress 0.1052647 
    ## Run 92 stress 0.08926104 
    ## ... Procrustes: rmse 0.0003332463  max resid 0.001092367 
    ## ... Similar to previous best
    ## Run 93 stress 0.08938538 
    ## ... Procrustes: rmse 0.01180353  max resid 0.03969748 
    ## Run 94 stress 0.0892609 
    ## ... Procrustes: rmse 0.0002310231  max resid 0.0007470169 
    ## ... Similar to previous best
    ## Run 95 stress 0.1101447 
    ## Run 96 stress 0.09018711 
    ## Run 97 stress 0.09087344 
    ## Run 98 stress 0.09018712 
    ## Run 99 stress 0.1075421 
    ## Run 100 stress 0.09612777 
    ## Run 101 stress 0.09592153 
    ## Run 102 stress 0.08926101 
    ## ... Procrustes: rmse 0.0003041835  max resid 0.0009462921 
    ## ... Similar to previous best
    ## Run 103 stress 0.08938964 
    ## ... Procrustes: rmse 0.03595149  max resid 0.1182843 
    ## Run 104 stress 0.08938966 
    ## ... Procrustes: rmse 0.03598506  max resid 0.1183469 
    ## Run 105 stress 0.09503451 
    ## Run 106 stress 0.1061302 
    ## Run 107 stress 0.1087863 
    ## Run 108 stress 0.09044608 
    ## Run 109 stress 0.08926095 
    ## ... Procrustes: rmse 0.0004385433  max resid 0.001311982 
    ## ... Similar to previous best
    ## Run 110 stress 0.1092128 
    ## Run 111 stress 0.1067375 
    ## Run 112 stress 0.0894665 
    ## ... Procrustes: rmse 0.03323488  max resid 0.1163696 
    ## Run 113 stress 0.08951719 
    ## ... Procrustes: rmse 0.03506189  max resid 0.1156446 
    ## Run 114 stress 0.09039129 
    ## Run 115 stress 0.1074432 
    ## Run 116 stress 0.08938962 
    ## ... Procrustes: rmse 0.03597564  max resid 0.1183259 
    ## Run 117 stress 0.09088603 
    ## Run 118 stress 0.0893854 
    ## ... Procrustes: rmse 0.01179948  max resid 0.0397136 
    ## Run 119 stress 0.08938966 
    ## ... Procrustes: rmse 0.03593926  max resid 0.118267 
    ## Run 120 stress 0.08926067 
    ## ... Procrustes: rmse 6.436444e-05  max resid 0.0002186741 
    ## ... Similar to previous best
    ## Run 121 stress 0.09178297 
    ## Run 122 stress 0.09775539 
    ## Run 123 stress 0.09087339 
    ## Run 124 stress 0.08938967 
    ## ... Procrustes: rmse 0.03594016  max resid 0.1182716 
    ## Run 125 stress 0.08926104 
    ## ... Procrustes: rmse 0.0003436486  max resid 0.001047044 
    ## ... Similar to previous best
    ## Run 126 stress 0.1087579 
    ## Run 127 stress 0.1074433 
    ## Run 128 stress 0.089261 
    ## ... Procrustes: rmse 0.0003171663  max resid 0.001000028 
    ## ... Similar to previous best
    ## Run 129 stress 0.10569 
    ## Run 130 stress 0.1061302 
    ## Run 131 stress 0.08938962 
    ## ... Procrustes: rmse 0.03595291  max resid 0.1182885 
    ## Run 132 stress 0.08946336 
    ## ... Procrustes: rmse 0.03749402  max resid 0.1177364 
    ## Run 133 stress 0.1056899 
    ## Run 134 stress 0.09503434 
    ## Run 135 stress 0.1075421 
    ## Run 136 stress 0.1101451 
    ## Run 137 stress 0.09087369 
    ## Run 138 stress 0.1105343 
    ## Run 139 stress 0.08938547 
    ## ... Procrustes: rmse 0.01186082  max resid 0.03976362 
    ## Run 140 stress 0.1071322 
    ## Run 141 stress 0.1071322 
    ## Run 142 stress 0.09039133 
    ## Run 143 stress 0.1060431 
    ## Run 144 stress 0.09503451 
    ## Run 145 stress 0.09018707 
    ## Run 146 stress 0.08946343 
    ## ... Procrustes: rmse 0.03748307  max resid 0.1177194 
    ## Run 147 stress 0.1056895 
    ## Run 148 stress 0.1062569 
    ## Run 149 stress 0.1092127 
    ## Run 150 stress 0.106257 
    ## Run 151 stress 0.08946663 
    ## ... Procrustes: rmse 0.03319671  max resid 0.1163069 
    ## Run 152 stress 0.09044611 
    ## Run 153 stress 0.1056909 
    ## Run 154 stress 0.105265 
    ## Run 155 stress 0.09178285 
    ## Run 156 stress 0.1071322 
    ## Run 157 stress 0.08926095 
    ## ... Procrustes: rmse 0.0002854301  max resid 0.0008943582 
    ## ... Similar to previous best
    ## Run 158 stress 0.1091831 
    ## Run 159 stress 0.08926069 
    ## ... Procrustes: rmse 0.0002357854  max resid 0.000767551 
    ## ... Similar to previous best
    ## Run 160 stress 0.09109109 
    ## Run 161 stress 0.08951721 
    ## ... Procrustes: rmse 0.03505375  max resid 0.115637 
    ## Run 162 stress 0.09594333 
    ## Run 163 stress 0.09039135 
    ## Run 164 stress 0.09021172 
    ## Run 165 stress 0.1074432 
    ## Run 166 stress 0.1074432 
    ## Run 167 stress 0.08926081 
    ## ... Procrustes: rmse 0.000204876  max resid 0.000638851 
    ## ... Similar to previous best
    ## Run 168 stress 0.106751 
    ## Run 169 stress 0.1056897 
    ## Run 170 stress 0.1052651 
    ## Run 171 stress 0.09775566 
    ## Run 172 stress 0.09088603 
    ## Run 173 stress 0.1067381 
    ## Run 174 stress 0.09087341 
    ## Run 175 stress 0.08938544 
    ## ... Procrustes: rmse 0.01185551  max resid 0.03974793 
    ## Run 176 stress 0.09044612 
    ## Run 177 stress 0.08938548 
    ## ... Procrustes: rmse 0.0117594  max resid 0.03962561 
    ## Run 178 stress 0.0893856 
    ## ... Procrustes: rmse 0.01187842  max resid 0.03963235 
    ## Run 179 stress 0.08938967 
    ## ... Procrustes: rmse 0.0359838  max resid 0.1183435 
    ## Run 180 stress 0.08926072 
    ## ... Procrustes: rmse 0.0001185248  max resid 0.0003820297 
    ## ... Similar to previous best
    ## Run 181 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183465  max resid 0.03969542 
    ## Run 182 stress 0.1076304 
    ## Run 183 stress 0.08938961 
    ## ... Procrustes: rmse 0.03596538  max resid 0.1183077 
    ## Run 184 stress 0.08926089 
    ## ... Procrustes: rmse 0.0002568574  max resid 0.0008266076 
    ## ... Similar to previous best
    ## Run 185 stress 0.08926079 
    ## ... Procrustes: rmse 0.0001823095  max resid 0.0005887006 
    ## ... Similar to previous best
    ## Run 186 stress 0.1092126 
    ## Run 187 stress 0.1105343 
    ## Run 188 stress 0.108839 
    ## Run 189 stress 0.08938969 
    ## ... Procrustes: rmse 0.03593622  max resid 0.1182641 
    ## Run 190 stress 0.109295 
    ## Run 191 stress 0.08926083 
    ## ... Procrustes: rmse 0.0001609044  max resid 0.0004983123 
    ## ... Similar to previous best
    ## Run 192 stress 0.09039135 
    ## Run 193 stress 0.08946668 
    ## ... Procrustes: rmse 0.03320962  max resid 0.1163334 
    ## Run 194 stress 0.08946654 
    ## ... Procrustes: rmse 0.03321293  max resid 0.1163295 
    ## Run 195 stress 0.1087572 
    ## Run 196 stress 0.1061296 
    ## Run 197 stress 0.08946653 
    ## ... Procrustes: rmse 0.03324983  max resid 0.1164125 
    ## Run 198 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 3.288293e-05  max resid 0.000101151 
    ## ... Similar to previous best
    ## Run 199 stress 0.08938967 
    ## ... Procrustes: rmse 0.03593948  max resid 0.1182681 
    ## Run 200 stress 0.08938967 
    ## ... Procrustes: rmse 0.03593258  max resid 0.1182587 
    ## Run 201 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183915  max resid 0.03968365 
    ## Run 202 stress 0.0894666 
    ## ... Procrustes: rmse 0.03319461  max resid 0.1163063 
    ## Run 203 stress 0.09503428 
    ## Run 204 stress 0.1071321 
    ## Run 205 stress 0.08926073 
    ## ... Procrustes: rmse 0.0001575529  max resid 0.0004862798 
    ## ... Similar to previous best
    ## Run 206 stress 0.1068389 
    ## Run 207 stress 0.1092363 
    ## Run 208 stress 0.08938978 
    ## ... Procrustes: rmse 0.03591482  max resid 0.1182274 
    ## Run 209 stress 0.08938971 
    ## ... Procrustes: rmse 0.03592407  max resid 0.1182443 
    ## Run 210 stress 0.08926066 
    ## ... Procrustes: rmse 0.0001056718  max resid 0.0002919293 
    ## ... Similar to previous best
    ## Run 211 stress 0.1074433 
    ## Run 212 stress 0.1086562 
    ## Run 213 stress 0.1112199 
    ## Run 214 stress 0.109145 
    ## Run 215 stress 0.1071322 
    ## Run 216 stress 0.1079014 
    ## Run 217 stress 0.09503445 
    ## Run 218 stress 0.08926141 
    ## ... Procrustes: rmse 0.0004723989  max resid 0.001509901 
    ## ... Similar to previous best
    ## Run 219 stress 0.0902117 
    ## Run 220 stress 0.08938962 
    ## ... Procrustes: rmse 0.03594829  max resid 0.1182807 
    ## Run 221 stress 0.1071321 
    ## Run 222 stress 0.09503417 
    ## Run 223 stress 0.08938962 
    ## ... Procrustes: rmse 0.03597037  max resid 0.1183179 
    ## Run 224 stress 0.08938963 
    ## ... Procrustes: rmse 0.03594493  max resid 0.118278 
    ## Run 225 stress 0.1108317 
    ## Run 226 stress 0.0908734 
    ## Run 227 stress 0.0910911 
    ## Run 228 stress 0.08946651 
    ## ... Procrustes: rmse 0.03323183  max resid 0.1163675 
    ## Run 229 stress 0.08926064 
    ## ... Procrustes: rmse 4.445566e-05  max resid 0.0001520668 
    ## ... Similar to previous best
    ## Run 230 stress 0.08926095 
    ## ... Procrustes: rmse 0.0003168341  max resid 0.001065323 
    ## ... Similar to previous best
    ## Run 231 stress 0.09039131 
    ## Run 232 stress 0.08946331 
    ## ... Procrustes: rmse 0.03750235  max resid 0.1177507 
    ## Run 233 stress 0.1056901 
    ## Run 234 stress 0.1071322 
    ## Run 235 stress 0.0893898 
    ## ... Procrustes: rmse 0.03591437  max resid 0.1182174 
    ## Run 236 stress 0.1071323 
    ## Run 237 stress 0.09503426 
    ## Run 238 stress 0.1052647 
    ## Run 239 stress 0.0904461 
    ## Run 240 stress 0.0893854 
    ## ... Procrustes: rmse 0.01184917  max resid 0.03972867 
    ## Run 241 stress 0.1087866 
    ## Run 242 stress 0.1087866 
    ## Run 243 stress 0.1056893 
    ## Run 244 stress 0.08938545 
    ## ... Procrustes: rmse 0.01180164  max resid 0.0396935 
    ## Run 245 stress 0.1093315 
    ## Run 246 stress 0.08926108 
    ## ... Procrustes: rmse 0.0004820147  max resid 0.001427959 
    ## ... Similar to previous best
    ## Run 247 stress 0.1056904 
    ## Run 248 stress 0.1079019 
    ## Run 249 stress 0.1079014 
    ## Run 250 stress 0.09039146 
    ## Run 251 stress 0.09503437 
    ## Run 252 stress 0.1071322 
    ## Run 253 stress 0.1056899 
    ## Run 254 stress 0.1063187 
    ## Run 255 stress 0.08938967 
    ## ... Procrustes: rmse 0.03593315  max resid 0.1182598 
    ## Run 256 stress 0.08946656 
    ## ... Procrustes: rmse 0.03320227  max resid 0.1163156 
    ## Run 257 stress 0.09021165 
    ## Run 258 stress 0.08938981 
    ## ... Procrustes: rmse 0.03591091  max resid 0.1182269 
    ## Run 259 stress 0.1071033 
    ## Run 260 stress 0.08926068 
    ## ... Procrustes: rmse 0.0001054112  max resid 0.000319306 
    ## ... Similar to previous best
    ## Run 261 stress 0.09099524 
    ## Run 262 stress 0.1052648 
    ## Run 263 stress 0.1056901 
    ## Run 264 stress 0.08926103 
    ## ... Procrustes: rmse 0.0003678403  max resid 0.001204281 
    ## ... Similar to previous best
    ## Run 265 stress 0.09044615 
    ## Run 266 stress 0.08946651 
    ## ... Procrustes: rmse 0.03323228  max resid 0.1163613 
    ## Run 267 stress 0.09503444 
    ## Run 268 stress 0.08946661 
    ## ... Procrustes: rmse 0.03319625  max resid 0.1163096 
    ## Run 269 stress 0.08938543 
    ## ... Procrustes: rmse 0.01179529  max resid 0.03970914 
    ## Run 270 stress 0.09018707 
    ## Run 271 stress 0.08938961 
    ## ... Procrustes: rmse 0.03596788  max resid 0.1183122 
    ## Run 272 stress 0.09109128 
    ## Run 273 stress 0.109188 
    ## Run 274 stress 0.0893856 
    ## ... Procrustes: rmse 0.01176778  max resid 0.03975413 
    ## Run 275 stress 0.10569 
    ## Run 276 stress 0.1056902 
    ## Run 277 stress 0.09503431 
    ## Run 278 stress 0.09039141 
    ## Run 279 stress 0.08926083 
    ## ... Procrustes: rmse 0.0002124156  max resid 0.0006958944 
    ## ... Similar to previous best
    ## Run 280 stress 0.0959431 
    ## Run 281 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003972241  max resid 0.001181646 
    ## ... Similar to previous best
    ## Run 282 stress 0.08938539 
    ## ... Procrustes: rmse 0.01185149  max resid 0.03967948 
    ## Run 283 stress 0.1056904 
    ## Run 284 stress 0.08926088 
    ## ... Procrustes: rmse 0.0002831882  max resid 0.0008880909 
    ## ... Similar to previous best
    ## Run 285 stress 0.08938977 
    ## ... Procrustes: rmse 0.03591676  max resid 0.1182317 
    ## Run 286 stress 0.1111367 
    ## Run 287 stress 0.0977556 
    ## Run 288 stress 0.09178286 
    ## Run 289 stress 0.1075422 
    ## Run 290 stress 0.08938555 
    ## ... Procrustes: rmse 0.01187909  max resid 0.03961767 
    ## Run 291 stress 0.09018708 
    ## Run 292 stress 0.08938566 
    ## ... Procrustes: rmse 0.01191108  max resid 0.03961956 
    ## Run 293 stress 0.09180873 
    ## Run 294 stress 0.08938543 
    ## ... Procrustes: rmse 0.01184335  max resid 0.03964302 
    ## Run 295 stress 0.09503413 
    ## Run 296 stress 0.091783 
    ## Run 297 stress 0.08938545 
    ## ... Procrustes: rmse 0.01186083  max resid 0.03974912 
    ## Run 298 stress 0.1108452 
    ## Run 299 stress 0.10569 
    ## Run 300 stress 0.1085127 
    ## Run 301 stress 0.0893855 
    ## ... Procrustes: rmse 0.01179158  max resid 0.0397235 
    ## Run 302 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595422  max resid 0.1182849 
    ## Run 303 stress 0.08946329 
    ## ... Procrustes: rmse 0.0375093  max resid 0.1177601 
    ## Run 304 stress 0.08938967 
    ## ... Procrustes: rmse 0.03597833  max resid 0.1183416 
    ## Run 305 stress 0.08946334 
    ## ... Procrustes: rmse 0.03749496  max resid 0.1177403 
    ## Run 306 stress 0.1067385 
    ## Run 307 stress 0.1084539 
    ## Run 308 stress 0.1060428 
    ## Run 309 stress 0.08938557 
    ## ... Procrustes: rmse 0.01177025  max resid 0.03974451 
    ## Run 310 stress 0.0893854 
    ## ... Procrustes: rmse 0.01184356  max resid 0.03972443 
    ## Run 311 stress 0.108656 
    ## Run 312 stress 0.0894633 
    ## ... Procrustes: rmse 0.03750604  max resid 0.1177576 
    ## Run 313 stress 0.1052649 
    ## Run 314 stress 0.08926114 
    ## ... Procrustes: rmse 0.0004233875  max resid 0.001309267 
    ## ... Similar to previous best
    ## Run 315 stress 0.09021165 
    ## Run 316 stress 0.0893855 
    ## ... Procrustes: rmse 0.0118839  max resid 0.03964648 
    ## Run 317 stress 0.09503423 
    ## Run 318 stress 0.08938981 
    ## ... Procrustes: rmse 0.03591361  max resid 0.1182265 
    ## Run 319 stress 0.08938966 
    ## ... Procrustes: rmse 0.0359351  max resid 0.1182603 
    ## Run 320 stress 0.0895172 
    ## ... Procrustes: rmse 0.0350431  max resid 0.1156211 
    ## Run 321 stress 0.08926071 
    ## ... Procrustes: rmse 0.0001816649  max resid 0.0004232455 
    ## ... Similar to previous best
    ## Run 322 stress 0.08938538 
    ## ... Procrustes: rmse 0.0118426  max resid 0.03968756 
    ## Run 323 stress 0.08938961 
    ## ... Procrustes: rmse 0.03596018  max resid 0.1182979 
    ## Run 324 stress 0.08946653 
    ## ... Procrustes: rmse 0.03321093  max resid 0.1163303 
    ## Run 325 stress 0.1075421 
    ## Run 326 stress 0.0908734 
    ## Run 327 stress 0.1056898 
    ## Run 328 stress 0.09099574 
    ## Run 329 stress 0.09087339 
    ## Run 330 stress 0.08926097 
    ## ... Procrustes: rmse 0.000323006  max resid 0.001079524 
    ## ... Similar to previous best
    ## Run 331 stress 0.08938539 
    ## ... Procrustes: rmse 0.01182805  max resid 0.03971713 
    ## Run 332 stress 0.08938962 
    ## ... Procrustes: rmse 0.03594617  max resid 0.1182784 
    ## Run 333 stress 0.08938541 
    ## ... Procrustes: rmse 0.01185902  max resid 0.03966188 
    ## Run 334 stress 0.109213 
    ## Run 335 stress 0.1061303 
    ## Run 336 stress 0.09228482 
    ## Run 337 stress 0.1056906 
    ## Run 338 stress 0.08938537 
    ## ... Procrustes: rmse 0.01182884  max resid 0.0396824 
    ## Run 339 stress 0.10613 
    ## Run 340 stress 0.09503413 
    ## Run 341 stress 0.1074434 
    ## Run 342 stress 0.09594329 
    ## Run 343 stress 0.08926064 
    ## ... Procrustes: rmse 8.565391e-05  max resid 0.0003021812 
    ## ... Similar to previous best
    ## Run 344 stress 0.09592152 
    ## Run 345 stress 0.1067809 
    ## Run 346 stress 0.09775595 
    ## Run 347 stress 0.08938981 
    ## ... Procrustes: rmse 0.03599805  max resid 0.1183553 
    ## Run 348 stress 0.09503436 
    ## Run 349 stress 0.08926092 
    ## ... Procrustes: rmse 0.0003038673  max resid 0.0009622088 
    ## ... Similar to previous best
    ## Run 350 stress 0.1074433 
    ## Run 351 stress 0.09503428 
    ## Run 352 stress 0.08946651 
    ## ... Procrustes: rmse 0.03322911  max resid 0.116358 
    ## Run 353 stress 0.08926077 
    ## ... Procrustes: rmse 0.000194948  max resid 0.0006170096 
    ## ... Similar to previous best
    ## Run 354 stress 0.09130085 
    ## Run 355 stress 0.0913011 
    ## Run 356 stress 0.08951717 
    ## ... Procrustes: rmse 0.03506176  max resid 0.1156498 
    ## Run 357 stress 0.08938545 
    ## ... Procrustes: rmse 0.01185633  max resid 0.03974194 
    ## Run 358 stress 0.09109115 
    ## Run 359 stress 0.1067491 
    ## Run 360 stress 0.0894666 
    ## ... Procrustes: rmse 0.03325999  max resid 0.116416 
    ## Run 361 stress 0.1071321 
    ## Run 362 stress 0.1075421 
    ## Run 363 stress 0.08938544 
    ## ... Procrustes: rmse 0.0118467  max resid 0.03966923 
    ## Run 364 stress 0.09087369 
    ## Run 365 stress 0.08926108 
    ## ... Procrustes: rmse 0.0003419634  max resid 0.001071775 
    ## ... Similar to previous best
    ## Run 366 stress 0.08938542 
    ## ... Procrustes: rmse 0.01184644  max resid 0.03972183 
    ## Run 367 stress 0.1096137 
    ## Run 368 stress 0.1060435 
    ## Run 369 stress 0.08926066 
    ## ... Procrustes: rmse 0.0001187479  max resid 0.0002882481 
    ## ... Similar to previous best
    ## Run 370 stress 0.1067377 
    ## Run 371 stress 0.110145 
    ## Run 372 stress 0.09021166 
    ## Run 373 stress 0.1052649 
    ## Run 374 stress 0.1052652 
    ## Run 375 stress 0.0904461 
    ## Run 376 stress 0.1085125 
    ## Run 377 stress 0.08938549 
    ## ... Procrustes: rmse 0.01186751  max resid 0.03976651 
    ## Run 378 stress 0.1074433 
    ## Run 379 stress 0.09018707 
    ## Run 380 stress 0.1060428 
    ## Run 381 stress 0.08946328 
    ## ... Procrustes: rmse 0.03751922  max resid 0.117781 
    ## Run 382 stress 0.1101445 
    ## Run 383 stress 0.08926064 
    ## ... Procrustes: rmse 0.0001271689  max resid 0.0003578053 
    ## ... Similar to previous best
    ## Run 384 stress 0.08926063 
    ## ... Procrustes: rmse 8.265326e-05  max resid 0.0001789717 
    ## ... Similar to previous best
    ## Run 385 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595998  max resid 0.1183086 
    ## Run 386 stress 0.08938964 
    ## ... Procrustes: rmse 0.03597596  max resid 0.1183301 
    ## Run 387 stress 0.08926113 
    ## ... Procrustes: rmse 0.0004152515  max resid 0.001300012 
    ## ... Similar to previous best
    ## Run 388 stress 0.1075423 
    ## Run 389 stress 0.08926113 
    ## ... Procrustes: rmse 0.0004131624  max resid 0.001369627 
    ## ... Similar to previous best
    ## Run 390 stress 0.0893854 
    ## ... Procrustes: rmse 0.01185659  max resid 0.03967475 
    ## Run 391 stress 0.1074431 
    ## Run 392 stress 0.1056892 
    ## Run 393 stress 0.1086566 
    ## Run 394 stress 0.09109108 
    ## Run 395 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595208  max resid 0.1182889 
    ## Run 396 stress 0.09612808 
    ## Run 397 stress 0.09018708 
    ## Run 398 stress 0.1052647 
    ## Run 399 stress 0.08946654 
    ## ... Procrustes: rmse 0.03320758  max resid 0.1163403 
    ## Run 400 stress 0.0902117 
    ## Run 401 stress 0.08938976 
    ## ... Procrustes: rmse 0.03591698  max resid 0.118234 
    ## Run 402 stress 0.09503418 
    ## Run 403 stress 0.1086564 
    ## Run 404 stress 0.09109128 
    ## Run 405 stress 0.1075421 
    ## Run 406 stress 0.08938556 
    ## ... Procrustes: rmse 0.01188597  max resid 0.03963825 
    ## Run 407 stress 0.0893854 
    ## ... Procrustes: rmse 0.01184959  max resid 0.03973036 
    ## Run 408 stress 0.08938964 
    ## ... Procrustes: rmse 0.03594068  max resid 0.1182703 
    ## Run 409 stress 0.08946331 
    ## ... Procrustes: rmse 0.03750246  max resid 0.1177764 
    ## Run 410 stress 0.09087346 
    ## Run 411 stress 0.1076302 
    ## Run 412 stress 0.08951717 
    ## ... Procrustes: rmse 0.03507144  max resid 0.1156572 
    ## Run 413 stress 0.1087567 
    ## Run 414 stress 0.08938961 
    ## ... Procrustes: rmse 0.0359681  max resid 0.1183126 
    ## Run 415 stress 0.1092128 
    ## Run 416 stress 0.1076305 
    ## Run 417 stress 0.08926084 
    ## ... Procrustes: rmse 0.0002263235  max resid 0.000701866 
    ## ... Similar to previous best
    ## Run 418 stress 0.08926098 
    ## ... Procrustes: rmse 0.0003307137  max resid 0.001084595 
    ## ... Similar to previous best
    ## Run 419 stress 0.08926069 
    ## ... Procrustes: rmse 0.0001038876  max resid 0.0002564159 
    ## ... Similar to previous best
    ## Run 420 stress 0.08926093 
    ## ... Procrustes: rmse 0.0004012339  max resid 0.001196222 
    ## ... Similar to previous best
    ## Run 421 stress 0.1063193 
    ## Run 422 stress 0.1071321 
    ## Run 423 stress 0.08938966 
    ## ... Procrustes: rmse 0.03593424  max resid 0.1182613 
    ## Run 424 stress 0.0892609 
    ## ... Procrustes: rmse 0.0003870278  max resid 0.001153135 
    ## ... Similar to previous best
    ## Run 425 stress 0.09099598 
    ## Run 426 stress 0.1087863 
    ## Run 427 stress 0.1071322 
    ## Run 428 stress 0.1074433 
    ## Run 429 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183699  max resid 0.03970476 
    ## Run 430 stress 0.1101442 
    ## Run 431 stress 0.1068603 
    ## Run 432 stress 0.1074433 
    ## Run 433 stress 0.09130086 
    ## Run 434 stress 0.08946335 
    ## ... Procrustes: rmse 0.03749183  max resid 0.1177319 
    ## Run 435 stress 0.1103716 
    ## Run 436 stress 0.08946336 
    ## ... Procrustes: rmse 0.03749148  max resid 0.1177324 
    ## Run 437 stress 0.09039135 
    ## Run 438 stress 0.1063196 
    ## Run 439 stress 0.1075421 
    ## Run 440 stress 0.0894665 
    ## ... Procrustes: rmse 0.03323093  max resid 0.116369 
    ## Run 441 stress 0.09021166 
    ## Run 442 stress 0.1074432 
    ## Run 443 stress 0.09039138 
    ## Run 444 stress 0.1056902 
    ## Run 445 stress 0.08938537 
    ## ... Procrustes: rmse 0.01183106  max resid 0.03969875 
    ## Run 446 stress 0.1060431 
    ## Run 447 stress 0.1108319 
    ## Run 448 stress 0.09109111 
    ## Run 449 stress 0.1071018 
    ## Run 450 stress 0.1075421 
    ## Run 451 stress 0.08951718 
    ## ... Procrustes: rmse 0.03508427  max resid 0.1156947 
    ## Run 452 stress 0.08946655 
    ## ... Procrustes: rmse 0.03320527  max resid 0.1163213 
    ## Run 453 stress 0.09021171 
    ## Run 454 stress 0.1075421 
    ## Run 455 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003158884  max resid 0.0009755345 
    ## ... Similar to previous best
    ## Run 456 stress 0.1091876 
    ## Run 457 stress 0.1071322 
    ## Run 458 stress 0.105691 
    ## Run 459 stress 0.08938961 
    ## ... Procrustes: rmse 0.03596871  max resid 0.1183155 
    ## Run 460 stress 0.09018711 
    ## Run 461 stress 0.1079016 
    ## Run 462 stress 0.1079014 
    ## Run 463 stress 0.1108635 
    ## Run 464 stress 0.08926086 
    ## ... Procrustes: rmse 0.0002714197  max resid 0.0008755974 
    ## ... Similar to previous best
    ## Run 465 stress 0.09018707 
    ## Run 466 stress 0.09021166 
    ## Run 467 stress 0.1085129 
    ## Run 468 stress 0.1052648 
    ## Run 469 stress 0.1075422 
    ## Run 470 stress 0.08926114 
    ## ... Procrustes: rmse 0.0004203599  max resid 0.001369203 
    ## ... Similar to previous best
    ## Run 471 stress 0.089261 
    ## ... Procrustes: rmse 0.0003478017  max resid 0.001156362 
    ## ... Similar to previous best
    ## Run 472 stress 0.109213 
    ## Run 473 stress 0.1108231 
    ## Run 474 stress 0.09503414 
    ## Run 475 stress 0.08946653 
    ## ... Procrustes: rmse 0.03324342  max resid 0.1164011 
    ## Run 476 stress 0.1071321 
    ## Run 477 stress 0.09178285 
    ## Run 478 stress 0.09503437 
    ## Run 479 stress 0.1052651 
    ## Run 480 stress 0.09178301 
    ## Run 481 stress 0.1087573 
    ## Run 482 stress 0.08946661 
    ## ... Procrustes: rmse 0.03319602  max resid 0.1163074 
    ## Run 483 stress 0.10569 
    ## Run 484 stress 0.09612779 
    ## Run 485 stress 0.1056899 
    ## Run 486 stress 0.08938962 
    ## ... Procrustes: rmse 0.03597087  max resid 0.1183187 
    ## Run 487 stress 0.1066197 
    ## Run 488 stress 0.08951716 
    ## ... Procrustes: rmse 0.03506166  max resid 0.1156448 
    ## Run 489 stress 0.09018717 
    ## Run 490 stress 0.1064198 
    ## Run 491 stress 0.09039135 
    ## Run 492 stress 0.08926081 
    ## ... Procrustes: rmse 0.0003247869  max resid 0.000970927 
    ## ... Similar to previous best
    ## Run 493 stress 0.1092094 
    ## Run 494 stress 0.09018708 
    ## Run 495 stress 0.08938546 
    ## ... Procrustes: rmse 0.01187496  max resid 0.03965324 
    ## Run 496 stress 0.09039135 
    ## Run 497 stress 0.09130091 
    ## Run 498 stress 0.08938539 
    ## ... Procrustes: rmse 0.01185809  max resid 0.03969068 
    ## Run 499 stress 0.1071321 
    ## Run 500 stress 0.09018707 
    ## *** Best solution repeated 34 times

``` r
# Mixed and stratified lakes
SD_beta_geo_MS_NMDS <- metaMDS(SD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09030406 
    ## Run 2 stress 0.09969964 
    ## Run 3 stress 0.09590892 
    ## Run 4 stress 0.3368708 
    ## Run 5 stress 0.08973873 
    ## Run 6 stress 0.0940799 
    ## Run 7 stress 0.09159088 
    ## Run 8 stress 0.09408001 
    ## Run 9 stress 0.09030403 
    ## Run 10 stress 0.09407965 
    ## Run 11 stress 0.100461 
    ## Run 12 stress 0.09030397 
    ## Run 13 stress 0.3131048 
    ## Run 14 stress 0.09465912 
    ## Run 15 stress 0.0937436 
    ## Run 16 stress 0.09159093 
    ## Run 17 stress 0.08440255 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01442908  max resid 0.04293066 
    ## Run 18 stress 0.0850348 
    ## Run 19 stress 0.09030397 
    ## Run 20 stress 0.08503493 
    ## Run 21 stress 0.09380598 
    ## Run 22 stress 0.09380592 
    ## Run 23 stress 0.09590885 
    ## Run 24 stress 0.09416138 
    ## Run 25 stress 0.09977535 
    ## Run 26 stress 0.09535554 
    ## Run 27 stress 0.0933728 
    ## Run 28 stress 0.09286083 
    ## Run 29 stress 0.09145323 
    ## Run 30 stress 0.09286088 
    ## Run 31 stress 0.08503458 
    ## Run 32 stress 0.0850349 
    ## Run 33 stress 0.09407975 
    ## Run 34 stress 0.09268366 
    ## Run 35 stress 0.1038045 
    ## Run 36 stress 0.09268411 
    ## Run 37 stress 0.09308967 
    ## Run 38 stress 0.09380582 
    ## Run 39 stress 0.09465906 
    ## Run 40 stress 0.09159091 
    ## Run 41 stress 0.09159083 
    ## Run 42 stress 0.09030399 
    ## Run 43 stress 0.09407982 
    ## Run 44 stress 0.08503467 
    ## Run 45 stress 0.09407997 
    ## Run 46 stress 0.09374242 
    ## Run 47 stress 0.09159088 
    ## Run 48 stress 0.09268388 
    ## Run 49 stress 0.09407983 
    ## Run 50 stress 0.09407965 
    ## Run 51 stress 0.08773469 
    ## Run 52 stress 0.09030397 
    ## Run 53 stress 0.09374253 
    ## Run 54 stress 0.08503468 
    ## Run 55 stress 0.105342 
    ## Run 56 stress 0.09374195 
    ## Run 57 stress 0.09337266 
    ## Run 58 stress 0.09407973 
    ## Run 59 stress 0.08773483 
    ## Run 60 stress 0.09268406 
    ## Run 61 stress 0.08773478 
    ## Run 62 stress 0.08973883 
    ## Run 63 stress 0.09030407 
    ## Run 64 stress 0.09380579 
    ## Run 65 stress 0.09268388 
    ## Run 66 stress 0.09030394 
    ## Run 67 stress 0.08503456 
    ## Run 68 stress 0.1038047 
    ## Run 69 stress 0.08503485 
    ## Run 70 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 3.900798e-05  max resid 7.483793e-05 
    ## ... Similar to previous best
    ## Run 71 stress 0.105285 
    ## Run 72 stress 0.0844027 
    ## ... Procrustes: rmse 0.0002103333  max resid 0.0003896178 
    ## ... Similar to previous best
    ## Run 73 stress 0.09416131 
    ## Run 74 stress 0.09268375 
    ## Run 75 stress 0.09308964 
    ## Run 76 stress 0.08440261 
    ## ... Procrustes: rmse 0.0002466466  max resid 0.000488812 
    ## ... Similar to previous best
    ## Run 77 stress 0.09374286 
    ## Run 78 stress 0.09407969 
    ## Run 79 stress 0.09760832 
    ## Run 80 stress 0.08440259 
    ## ... Procrustes: rmse 0.0002179749  max resid 0.0004495898 
    ## ... Similar to previous best
    ## Run 81 stress 0.1026215 
    ## Run 82 stress 0.08503465 
    ## Run 83 stress 0.08773468 
    ## Run 84 stress 0.09286083 
    ## Run 85 stress 0.08973867 
    ## Run 86 stress 0.0926837 
    ## Run 87 stress 0.09590895 
    ## Run 88 stress 0.09159083 
    ## Run 89 stress 0.09407965 
    ## Run 90 stress 0.08503635 
    ## Run 91 stress 0.09145326 
    ## Run 92 stress 0.09145331 
    ## Run 93 stress 0.09268377 
    ## Run 94 stress 0.09407984 
    ## Run 95 stress 0.09407989 
    ## Run 96 stress 0.09407969 
    ## Run 97 stress 0.09159085 
    ## Run 98 stress 0.09030398 
    ## Run 99 stress 0.0953546 
    ## Run 100 stress 0.09030396 
    ## Run 101 stress 0.0933729 
    ## Run 102 stress 0.09381874 
    ## Run 103 stress 0.09030404 
    ## Run 104 stress 0.09403415 
    ## Run 105 stress 0.08973867 
    ## Run 106 stress 0.09286083 
    ## Run 107 stress 0.09969964 
    ## Run 108 stress 0.09465907 
    ## Run 109 stress 0.09408006 
    ## Run 110 stress 0.09407968 
    ## Run 111 stress 0.09408024 
    ## Run 112 stress 0.1006352 
    ## Run 113 stress 0.09286086 
    ## Run 114 stress 0.08503566 
    ## Run 115 stress 0.09030403 
    ## Run 116 stress 0.09321727 
    ## Run 117 stress 0.08503514 
    ## Run 118 stress 0.09381852 
    ## Run 119 stress 0.09145322 
    ## Run 120 stress 0.09286088 
    ## Run 121 stress 0.09030399 
    ## Run 122 stress 0.09030395 
    ## Run 123 stress 0.1038048 
    ## Run 124 stress 0.09417915 
    ## Run 125 stress 0.08773479 
    ## Run 126 stress 0.08973895 
    ## Run 127 stress 0.08440264 
    ## ... Procrustes: rmse 0.000171353  max resid 0.0002991885 
    ## ... Similar to previous best
    ## Run 128 stress 0.08973885 
    ## Run 129 stress 0.09286092 
    ## Run 130 stress 0.09286085 
    ## Run 131 stress 0.09408006 
    ## Run 132 stress 0.08773484 
    ## Run 133 stress 0.09159101 
    ## Run 134 stress 0.08773474 
    ## Run 135 stress 0.09416133 
    ## Run 136 stress 0.09590926 
    ## Run 137 stress 0.09535493 
    ## Run 138 stress 0.08973863 
    ## Run 139 stress 0.08503475 
    ## Run 140 stress 0.09416136 
    ## Run 141 stress 0.08773473 
    ## Run 142 stress 0.09030401 
    ## Run 143 stress 0.09374219 
    ## Run 144 stress 0.1019601 
    ## Run 145 stress 0.09308974 
    ## Run 146 stress 0.09030394 
    ## Run 147 stress 0.09268356 
    ## Run 148 stress 0.09159094 
    ## Run 149 stress 0.09721198 
    ## Run 150 stress 0.09400542 
    ## Run 151 stress 0.0996998 
    ## Run 152 stress 0.09145331 
    ## Run 153 stress 0.09417793 
    ## Run 154 stress 0.09030401 
    ## Run 155 stress 0.09308974 
    ## Run 156 stress 0.09407978 
    ## Run 157 stress 0.08440261 
    ## ... Procrustes: rmse 0.0001271429  max resid 0.0002311091 
    ## ... Similar to previous best
    ## Run 158 stress 0.09337257 
    ## Run 159 stress 0.09168928 
    ## Run 160 stress 0.08503485 
    ## Run 161 stress 0.2914222 
    ## Run 162 stress 0.09308955 
    ## Run 163 stress 0.09416136 
    ## Run 164 stress 0.09408003 
    ## Run 165 stress 0.0996998 
    ## Run 166 stress 0.1038031 
    ## Run 167 stress 0.0940713 
    ## Run 168 stress 0.09465905 
    ## Run 169 stress 0.09721199 
    ## Run 170 stress 0.09407974 
    ## Run 171 stress 0.09286095 
    ## Run 172 stress 0.09030412 
    ## Run 173 stress 0.0903042 
    ## Run 174 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001940337  max resid 0.0004182207 
    ## ... Similar to previous best
    ## Run 175 stress 0.09407978 
    ## Run 176 stress 0.08503468 
    ## Run 177 stress 0.09030401 
    ## Run 178 stress 0.09030397 
    ## Run 179 stress 0.09760818 
    ## Run 180 stress 0.0877348 
    ## Run 181 stress 0.09308984 
    ## Run 182 stress 0.09308966 
    ## Run 183 stress 0.09159086 
    ## Run 184 stress 0.09268327 
    ## Run 185 stress 0.09268359 
    ## Run 186 stress 0.08973876 
    ## Run 187 stress 0.08503462 
    ## Run 188 stress 0.08773466 
    ## Run 189 stress 0.09969967 
    ## Run 190 stress 0.09464486 
    ## Run 191 stress 0.09535415 
    ## Run 192 stress 0.09168938 
    ## Run 193 stress 0.09030393 
    ## Run 194 stress 0.09407974 
    ## Run 195 stress 0.08973878 
    ## Run 196 stress 0.09403431 
    ## Run 197 stress 0.09145333 
    ## Run 198 stress 0.08773497 
    ## Run 199 stress 0.09030394 
    ## Run 200 stress 0.09308973 
    ## Run 201 stress 0.09168931 
    ## Run 202 stress 0.103804 
    ## Run 203 stress 0.08773467 
    ## Run 204 stress 0.09030397 
    ## Run 205 stress 0.09268378 
    ## Run 206 stress 0.09465907 
    ## Run 207 stress 0.08503603 
    ## Run 208 stress 0.09374278 
    ## Run 209 stress 0.09381844 
    ## Run 210 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001029808  max resid 0.0001824764 
    ## ... Similar to previous best
    ## Run 211 stress 0.0897387 
    ## Run 212 stress 0.08440267 
    ## ... Procrustes: rmse 0.0002673438  max resid 0.0005354047 
    ## ... Similar to previous best
    ## Run 213 stress 0.09030402 
    ## Run 214 stress 0.08973869 
    ## Run 215 stress 0.09337293 
    ## Run 216 stress 0.09268404 
    ## Run 217 stress 0.337224 
    ## Run 218 stress 0.09400516 
    ## Run 219 stress 0.09337238 
    ## Run 220 stress 0.09030398 
    ## Run 221 stress 0.08773473 
    ## Run 222 stress 0.09381844 
    ## Run 223 stress 0.1038045 
    ## Run 224 stress 0.0914533 
    ## Run 225 stress 0.09030397 
    ## Run 226 stress 0.09030401 
    ## Run 227 stress 0.08773467 
    ## Run 228 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 9.415591e-05  max resid 0.0001913188 
    ## ... Similar to previous best
    ## Run 229 stress 0.08773466 
    ## Run 230 stress 0.09407967 
    ## Run 231 stress 0.08773479 
    ## Run 232 stress 0.09286083 
    ## Run 233 stress 0.08773478 
    ## Run 234 stress 0.09407138 
    ## Run 235 stress 0.09030401 
    ## Run 236 stress 0.09159085 
    ## Run 237 stress 0.0877348 
    ## Run 238 stress 0.08973869 
    ## Run 239 stress 0.09145321 
    ## Run 240 stress 0.0937429 
    ## Run 241 stress 0.09030413 
    ## Run 242 stress 0.09286083 
    ## Run 243 stress 0.340943 
    ## Run 244 stress 0.08973864 
    ## Run 245 stress 0.09286095 
    ## Run 246 stress 0.08973881 
    ## Run 247 stress 0.08503456 
    ## Run 248 stress 0.09465906 
    ## Run 249 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001582749  max resid 0.0003315048 
    ## ... Similar to previous best
    ## Run 250 stress 0.09268349 
    ## Run 251 stress 0.09159084 
    ## Run 252 stress 0.09407131 
    ## Run 253 stress 0.09159088 
    ## Run 254 stress 0.09030395 
    ## Run 255 stress 0.09374305 
    ## Run 256 stress 0.09308961 
    ## Run 257 stress 0.0976084 
    ## Run 258 stress 0.09465904 
    ## Run 259 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001720534  max resid 0.0003321848 
    ## ... Similar to previous best
    ## Run 260 stress 0.09407976 
    ## Run 261 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001590761  max resid 0.0003202713 
    ## ... Similar to previous best
    ## Run 262 stress 0.09159096 
    ## Run 263 stress 0.09407993 
    ## Run 264 stress 0.09408009 
    ## Run 265 stress 0.09407985 
    ## Run 266 stress 0.09374236 
    ## Run 267 stress 0.09465904 
    ## Run 268 stress 0.09381856 
    ## Run 269 stress 0.09268384 
    ## Run 270 stress 0.1052849 
    ## Run 271 stress 0.09407973 
    ## Run 272 stress 0.09145321 
    ## Run 273 stress 0.08773482 
    ## Run 274 stress 0.09669905 
    ## Run 275 stress 0.09407973 
    ## Run 276 stress 0.09407979 
    ## Run 277 stress 0.103804 
    ## Run 278 stress 0.09159097 
    ## Run 279 stress 0.0850367 
    ## Run 280 stress 0.08973864 
    ## Run 281 stress 0.09168949 
    ## Run 282 stress 0.09374332 
    ## Run 283 stress 0.09969976 
    ## Run 284 stress 0.09407966 
    ## Run 285 stress 0.1052848 
    ## Run 286 stress 0.09030398 
    ## Run 287 stress 0.08773477 
    ## Run 288 stress 0.08973862 
    ## Run 289 stress 0.09168957 
    ## Run 290 stress 0.09416133 
    ## Run 291 stress 0.09407452 
    ## Run 292 stress 0.09465911 
    ## Run 293 stress 0.09380577 
    ## Run 294 stress 0.09321657 
    ## Run 295 stress 0.09400548 
    ## Run 296 stress 0.09268341 
    ## Run 297 stress 0.09465905 
    ## Run 298 stress 0.090304 
    ## Run 299 stress 0.09465907 
    ## Run 300 stress 0.09145327 
    ## Run 301 stress 0.0850351 
    ## Run 302 stress 0.09337239 
    ## Run 303 stress 0.08440269 
    ## ... Procrustes: rmse 0.000297817  max resid 0.0005487075 
    ## ... Similar to previous best
    ## Run 304 stress 0.09407986 
    ## Run 305 stress 0.09145324 
    ## Run 306 stress 0.08503518 
    ## Run 307 stress 0.08503494 
    ## Run 308 stress 0.08973875 
    ## Run 309 stress 0.09286096 
    ## Run 310 stress 0.09465904 
    ## Run 311 stress 0.09168955 
    ## Run 312 stress 0.09286083 
    ## Run 313 stress 0.08773484 
    ## Run 314 stress 0.09145321 
    ## Run 315 stress 0.09308957 
    ## Run 316 stress 0.09969965 
    ## Run 317 stress 0.09407994 
    ## Run 318 stress 0.09374275 
    ## Run 319 stress 0.09407971 
    ## Run 320 stress 0.08973871 
    ## Run 321 stress 0.08773473 
    ## Run 322 stress 0.09286085 
    ## Run 323 stress 0.09380583 
    ## Run 324 stress 0.09407999 
    ## Run 325 stress 0.09030403 
    ## Run 326 stress 0.09268394 
    ## Run 327 stress 0.08503613 
    ## Run 328 stress 0.09308946 
    ## Run 329 stress 0.0946591 
    ## Run 330 stress 0.08503465 
    ## Run 331 stress 0.09159087 
    ## Run 332 stress 0.09407981 
    ## Run 333 stress 0.09416137 
    ## Run 334 stress 0.09168946 
    ## Run 335 stress 0.085035 
    ## Run 336 stress 0.09407991 
    ## Run 337 stress 0.09145321 
    ## Run 338 stress 0.2804667 
    ## Run 339 stress 0.1038048 
    ## Run 340 stress 0.08503504 
    ## Run 341 stress 0.09407966 
    ## Run 342 stress 0.1006347 
    ## Run 343 stress 0.08440256 
    ## ... Procrustes: rmse 6.686357e-05  max resid 0.0001309008 
    ## ... Similar to previous best
    ## Run 344 stress 0.09308972 
    ## Run 345 stress 0.09721199 
    ## Run 346 stress 0.1038037 
    ## Run 347 stress 0.09308979 
    ## Run 348 stress 0.09374271 
    ## Run 349 stress 0.0897387 
    ## Run 350 stress 0.09969961 
    ## Run 351 stress 0.08503472 
    ## Run 352 stress 0.090304 
    ## Run 353 stress 0.09268334 
    ## Run 354 stress 0.08503573 
    ## Run 355 stress 0.09590607 
    ## Run 356 stress 0.09465914 
    ## Run 357 stress 0.103804 
    ## Run 358 stress 0.08973891 
    ## Run 359 stress 0.08503459 
    ## Run 360 stress 0.09145333 
    ## Run 361 stress 0.08973882 
    ## Run 362 stress 0.09308973 
    ## Run 363 stress 0.09321743 
    ## Run 364 stress 0.1104214 
    ## Run 365 stress 0.09030414 
    ## Run 366 stress 0.09669905 
    ## Run 367 stress 0.09145326 
    ## Run 368 stress 0.08503557 
    ## Run 369 stress 0.0877347 
    ## Run 370 stress 0.09030418 
    ## Run 371 stress 0.09417767 
    ## Run 372 stress 0.09145335 
    ## Run 373 stress 0.09381858 
    ## Run 374 stress 0.09168935 
    ## Run 375 stress 0.09168947 
    ## Run 376 stress 0.09407975 
    ## Run 377 stress 0.1052848 
    ## Run 378 stress 0.09145322 
    ## Run 379 stress 0.08773483 
    ## Run 380 stress 0.09407974 
    ## Run 381 stress 0.09268371 
    ## Run 382 stress 0.09721204 
    ## Run 383 stress 0.09030397 
    ## Run 384 stress 0.08503565 
    ## Run 385 stress 0.09030398 
    ## Run 386 stress 0.09286085 
    ## Run 387 stress 0.08973869 
    ## Run 388 stress 0.08440267 
    ## ... Procrustes: rmse 0.0001736626  max resid 0.0003912158 
    ## ... Similar to previous best
    ## Run 389 stress 0.09403418 
    ## Run 390 stress 0.09337259 
    ## Run 391 stress 0.09145334 
    ## Run 392 stress 0.09969972 
    ## Run 393 stress 0.09539212 
    ## Run 394 stress 0.09168966 
    ## Run 395 stress 0.09590945 
    ## Run 396 stress 0.1038039 
    ## Run 397 stress 0.08503502 
    ## Run 398 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001581501  max resid 0.0003011043 
    ## ... Similar to previous best
    ## Run 399 stress 0.0940797 
    ## Run 400 stress 0.09464471 
    ## Run 401 stress 0.08503539 
    ## Run 402 stress 0.09268359 
    ## Run 403 stress 0.08773475 
    ## Run 404 stress 0.08503469 
    ## Run 405 stress 0.09159103 
    ## Run 406 stress 0.09407967 
    ## Run 407 stress 0.09416139 
    ## Run 408 stress 0.09407981 
    ## Run 409 stress 0.09268381 
    ## Run 410 stress 0.09407968 
    ## Run 411 stress 0.09030393 
    ## Run 412 stress 0.09321863 
    ## Run 413 stress 0.08503477 
    ## Run 414 stress 0.1026215 
    ## Run 415 stress 0.09286091 
    ## Run 416 stress 0.08973865 
    ## Run 417 stress 0.09030393 
    ## Run 418 stress 0.0926837 
    ## Run 419 stress 0.08440273 
    ## ... Procrustes: rmse 0.0002985945  max resid 0.0005468971 
    ## ... Similar to previous best
    ## Run 420 stress 0.09286083 
    ## Run 421 stress 0.09969973 
    ## Run 422 stress 0.09969973 
    ## Run 423 stress 0.09465906 
    ## Run 424 stress 0.09268422 
    ## Run 425 stress 0.09380575 
    ## Run 426 stress 0.09159084 
    ## Run 427 stress 0.246648 
    ## Run 428 stress 0.08973882 
    ## Run 429 stress 0.09374175 
    ## Run 430 stress 0.09337272 
    ## Run 431 stress 0.09030397 
    ## Run 432 stress 0.09268363 
    ## Run 433 stress 0.09760813 
    ## Run 434 stress 0.09408005 
    ## Run 435 stress 0.09400533 
    ## Run 436 stress 0.09168949 
    ## Run 437 stress 0.09159093 
    ## Run 438 stress 0.09535545 
    ## Run 439 stress 0.09407992 
    ## Run 440 stress 0.08973867 
    ## Run 441 stress 0.09721203 
    ## Run 442 stress 0.09268362 
    ## Run 443 stress 0.08973872 
    ## Run 444 stress 0.08503467 
    ## Run 445 stress 0.09374314 
    ## Run 446 stress 0.09416133 
    ## Run 447 stress 0.08503461 
    ## Run 448 stress 0.09760816 
    ## Run 449 stress 0.09407982 
    ## Run 450 stress 0.09447305 
    ## Run 451 stress 0.09969974 
    ## Run 452 stress 0.0940798 
    ## Run 453 stress 0.08973868 
    ## Run 454 stress 0.08973868 
    ## Run 455 stress 0.09030399 
    ## Run 456 stress 0.3309726 
    ## Run 457 stress 0.08503458 
    ## Run 458 stress 0.09337256 
    ## Run 459 stress 0.09590938 
    ## Run 460 stress 0.08773482 
    ## Run 461 stress 0.09669886 
    ## Run 462 stress 0.09908195 
    ## Run 463 stress 0.09374387 
    ## Run 464 stress 0.08773468 
    ## Run 465 stress 0.08773472 
    ## Run 466 stress 0.09969979 
    ## Run 467 stress 0.08440269 
    ## ... Procrustes: rmse 0.0002883389  max resid 0.0005355342 
    ## ... Similar to previous best
    ## Run 468 stress 0.09407965 
    ## Run 469 stress 0.09374264 
    ## Run 470 stress 0.08503639 
    ## Run 471 stress 0.09145325 
    ## Run 472 stress 0.08973873 
    ## Run 473 stress 0.09539215 
    ## Run 474 stress 0.09535511 
    ## Run 475 stress 0.09145322 
    ## Run 476 stress 0.09030394 
    ## Run 477 stress 0.08973884 
    ## Run 478 stress 0.1026216 
    ## Run 479 stress 0.08973896 
    ## Run 480 stress 0.09145327 
    ## Run 481 stress 0.08503472 
    ## Run 482 stress 0.09712929 
    ## Run 483 stress 0.2914089 
    ## Run 484 stress 0.0937431 
    ## Run 485 stress 0.08973868 
    ## Run 486 stress 0.08773474 
    ## Run 487 stress 0.09969992 
    ## Run 488 stress 0.09403426 
    ## Run 489 stress 0.09145331 
    ## Run 490 stress 0.09308966 
    ## Run 491 stress 0.09159083 
    ## Run 492 stress 0.09969975 
    ## Run 493 stress 0.09969973 
    ## Run 494 stress 0.0940798 
    ## Run 495 stress 0.09321445 
    ## Run 496 stress 0.09030394 
    ## Run 497 stress 0.0850349 
    ## Run 498 stress 0.0903042 
    ## Run 499 stress 0.0961073 
    ## Run 500 stress 0.09268401 
    ## *** Best solution repeated 10 times

``` r
# Ocean sites and mixed lakes
SD_beta_geo_OM_NMDS <- metaMDS(SD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.07365789 
    ## ... Procrustes: rmse 9.778151e-05  max resid 0.0002341448 
    ## ... Similar to previous best
    ## Run 2 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743715  max resid 0.05104384 
    ## Run 3 stress 0.07378233 
    ## ... Procrustes: rmse 0.01743601  max resid 0.0510454 
    ## Run 4 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743654  max resid 0.05104323 
    ## Run 5 stress 0.07629236 
    ## Run 6 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001107722  max resid 0.0002634582 
    ## ... Similar to previous best
    ## Run 7 stress 0.07365784 
    ## ... Procrustes: rmse 3.055253e-05  max resid 7.18721e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.07365783 
    ## ... Procrustes: rmse 1.80354e-05  max resid 4.28204e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.07378233 
    ## ... Procrustes: rmse 0.01743708  max resid 0.05096773 
    ## Run 10 stress 0.07365787 
    ## ... Procrustes: rmse 9.50394e-05  max resid 0.0002250483 
    ## ... Similar to previous best
    ## Run 11 stress 0.07629233 
    ## Run 12 stress 0.08000468 
    ## Run 13 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744388  max resid 0.05105194 
    ## Run 14 stress 0.07629236 
    ## Run 15 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744518  max resid 0.05104835 
    ## Run 16 stress 0.07629234 
    ## Run 17 stress 0.07365784 
    ## ... Procrustes: rmse 4.016299e-05  max resid 9.474154e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744308  max resid 0.05104221 
    ## Run 19 stress 0.07378232 
    ## ... Procrustes: rmse 0.01742703  max resid 0.05100505 
    ## Run 20 stress 0.07629233 
    ## Run 21 stress 0.08000473 
    ## Run 22 stress 0.07629235 
    ## Run 23 stress 0.07378229 
    ## ... Procrustes: rmse 0.01746069  max resid 0.05108928 
    ## Run 24 stress 0.07629252 
    ## Run 25 stress 0.08000468 
    ## Run 26 stress 0.08000467 
    ## Run 27 stress 0.08000468 
    ## Run 28 stress 0.08000469 
    ## Run 29 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743124  max resid 0.05098186 
    ## Run 30 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744339  max resid 0.05102609 
    ## Run 31 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744383  max resid 0.05104193 
    ## Run 32 stress 0.07365784 
    ## ... Procrustes: rmse 1.681741e-05  max resid 3.737185e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.07365785 
    ## ... Procrustes: rmse 7.230507e-05  max resid 0.0001657678 
    ## ... Similar to previous best
    ## Run 34 stress 0.07629234 
    ## Run 35 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744166  max resid 0.05101537 
    ## Run 36 stress 0.07365783 
    ## ... Procrustes: rmse 9.031818e-06  max resid 1.593063e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.07629244 
    ## Run 38 stress 0.08288037 
    ## Run 39 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744995  max resid 0.05107198 
    ## Run 40 stress 0.07629237 
    ## Run 41 stress 0.07365784 
    ## ... Procrustes: rmse 5.436211e-05  max resid 0.0001279511 
    ## ... Similar to previous best
    ## Run 42 stress 0.08000475 
    ## Run 43 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744252  max resid 0.05103519 
    ## Run 44 stress 0.07629237 
    ## Run 45 stress 0.07365784 
    ## ... Procrustes: rmse 3.506226e-05  max resid 8.242658e-05 
    ## ... Similar to previous best
    ## Run 46 stress 0.07365789 
    ## ... Procrustes: rmse 9.798387e-05  max resid 0.0002188052 
    ## ... Similar to previous best
    ## Run 47 stress 0.07629234 
    ## Run 48 stress 0.07629238 
    ## Run 49 stress 0.07732944 
    ## Run 50 stress 0.07629241 
    ## Run 51 stress 0.08288038 
    ## Run 52 stress 0.08000467 
    ## Run 53 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174522  max resid 0.05108658 
    ## Run 54 stress 0.07365783 
    ## ... Procrustes: rmse 7.465216e-06  max resid 1.469343e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.0762924 
    ## Run 56 stress 0.07732937 
    ## Run 57 stress 0.07365785 
    ## ... Procrustes: rmse 5.620199e-05  max resid 0.0001290042 
    ## ... Similar to previous best
    ## Run 58 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744237  max resid 0.05103019 
    ## Run 59 stress 0.07629235 
    ## Run 60 stress 0.07629245 
    ## Run 61 stress 0.07629245 
    ## Run 62 stress 0.08288045 
    ## Run 63 stress 0.07629233 
    ## Run 64 stress 0.07629234 
    ## Run 65 stress 0.07629234 
    ## Run 66 stress 0.0823376 
    ## Run 67 stress 0.07629234 
    ## Run 68 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744251  max resid 0.0510302 
    ## Run 69 stress 0.07365785 
    ## ... Procrustes: rmse 6.675655e-05  max resid 0.0001540244 
    ## ... Similar to previous best
    ## Run 70 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001153322  max resid 0.0002701782 
    ## ... Similar to previous best
    ## Run 71 stress 0.08233756 
    ## Run 72 stress 0.08233774 
    ## Run 73 stress 0.07365784 
    ## ... Procrustes: rmse 4.967592e-05  max resid 0.0001185606 
    ## ... Similar to previous best
    ## Run 74 stress 0.07378229 
    ## ... Procrustes: rmse 0.01743957  max resid 0.05099277 
    ## Run 75 stress 0.08233775 
    ## Run 76 stress 0.07629243 
    ## Run 77 stress 0.07365785 
    ## ... Procrustes: rmse 8.596528e-05  max resid 0.0002019016 
    ## ... Similar to previous best
    ## Run 78 stress 0.07365785 
    ## ... Procrustes: rmse 7.141511e-05  max resid 0.000166923 
    ## ... Similar to previous best
    ## Run 79 stress 0.0762924 
    ## Run 80 stress 0.07732928 
    ## Run 81 stress 0.0800047 
    ## Run 82 stress 0.07732926 
    ## Run 83 stress 0.08000471 
    ## Run 84 stress 0.3495733 
    ## Run 85 stress 0.08000467 
    ## Run 86 stress 0.08233758 
    ## Run 87 stress 0.07629238 
    ## Run 88 stress 0.07365785 
    ## ... Procrustes: rmse 6.058623e-05  max resid 0.0001432932 
    ## ... Similar to previous best
    ## Run 89 stress 0.07629237 
    ## Run 90 stress 0.07365784 
    ## ... Procrustes: rmse 5.615601e-05  max resid 0.0001310854 
    ## ... Similar to previous best
    ## Run 91 stress 0.07629234 
    ## Run 92 stress 0.07365784 
    ## ... Procrustes: rmse 3.941945e-05  max resid 9.011483e-05 
    ## ... Similar to previous best
    ## Run 93 stress 0.07365796 
    ## ... Procrustes: rmse 0.0001877627  max resid 0.0004409634 
    ## ... Similar to previous best
    ## Run 94 stress 0.07365784 
    ## ... Procrustes: rmse 2.37954e-05  max resid 4.882513e-05 
    ## ... Similar to previous best
    ## Run 95 stress 0.07629235 
    ## Run 96 stress 0.08000468 
    ## Run 97 stress 0.07365786 
    ## ... Procrustes: rmse 8.646797e-05  max resid 0.0002018622 
    ## ... Similar to previous best
    ## Run 98 stress 0.07732926 
    ## Run 99 stress 0.07365785 
    ## ... Procrustes: rmse 6.729375e-05  max resid 0.0001596299 
    ## ... Similar to previous best
    ## Run 100 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001114424  max resid 0.0002593598 
    ## ... Similar to previous best
    ## Run 101 stress 0.08288044 
    ## Run 102 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744158  max resid 0.05101786 
    ## Run 103 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743436  max resid 0.05097275 
    ## Run 104 stress 0.08000471 
    ## Run 105 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174187  max resid 0.05095428 
    ## Run 106 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743822  max resid 0.05104767 
    ## Run 107 stress 0.08000471 
    ## Run 108 stress 0.07732929 
    ## Run 109 stress 0.07629246 
    ## Run 110 stress 0.07629233 
    ## Run 111 stress 0.07365784 
    ## ... Procrustes: rmse 4.524159e-05  max resid 0.0001046394 
    ## ... Similar to previous best
    ## Run 112 stress 0.07365784 
    ## ... Procrustes: rmse 1.411117e-05  max resid 2.196264e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.3352968 
    ## Run 114 stress 0.07365785 
    ## ... Procrustes: rmse 5.369776e-05  max resid 0.0001194222 
    ## ... Similar to previous best
    ## Run 115 stress 0.07629246 
    ## Run 116 stress 0.07365784 
    ## ... Procrustes: rmse 3.104695e-05  max resid 7.262791e-05 
    ## ... Similar to previous best
    ## Run 117 stress 0.07365784 
    ## ... Procrustes: rmse 6.115662e-05  max resid 0.000143521 
    ## ... Similar to previous best
    ## Run 118 stress 0.08000468 
    ## Run 119 stress 0.07365785 
    ## ... Procrustes: rmse 6.426136e-05  max resid 0.0001524064 
    ## ... Similar to previous best
    ## Run 120 stress 0.07629234 
    ## Run 121 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744448  max resid 0.05105578 
    ## Run 122 stress 0.07629234 
    ## Run 123 stress 0.08233767 
    ## Run 124 stress 0.08000468 
    ## Run 125 stress 0.07365784 
    ## ... Procrustes: rmse 5.353228e-05  max resid 0.0001251709 
    ## ... Similar to previous best
    ## Run 126 stress 0.07629236 
    ## Run 127 stress 0.07732933 
    ## Run 128 stress 0.08000471 
    ## Run 129 stress 0.07365785 
    ## ... Procrustes: rmse 6.247956e-05  max resid 0.0001471133 
    ## ... Similar to previous best
    ## Run 130 stress 0.07378235 
    ## ... Procrustes: rmse 0.01745855  max resid 0.05105861 
    ## Run 131 stress 0.07629236 
    ## Run 132 stress 0.07365796 
    ## ... Procrustes: rmse 0.0001777851  max resid 0.0004234639 
    ## ... Similar to previous best
    ## Run 133 stress 0.08000472 
    ## Run 134 stress 0.07629246 
    ## Run 135 stress 0.08233769 
    ## Run 136 stress 0.07629234 
    ## Run 137 stress 0.07629236 
    ## Run 138 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745554  max resid 0.05108403 
    ## Run 139 stress 0.07629233 
    ## Run 140 stress 0.08288034 
    ## Run 141 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744381  max resid 0.05105754 
    ## Run 142 stress 0.08000468 
    ## Run 143 stress 0.07629233 
    ## Run 144 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745275  max resid 0.0510691 
    ## Run 145 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 8.268065e-06  max resid 1.856655e-05 
    ## ... Similar to previous best
    ## Run 146 stress 0.07629239 
    ## Run 147 stress 0.07629244 
    ## Run 148 stress 0.07365783 
    ## ... Procrustes: rmse 3.38048e-05  max resid 7.92531e-05 
    ## ... Similar to previous best
    ## Run 149 stress 0.07629236 
    ## Run 150 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001573913  max resid 0.0003646453 
    ## ... Similar to previous best
    ## Run 151 stress 0.0736579 
    ## ... Procrustes: rmse 0.000143192  max resid 0.0003389423 
    ## ... Similar to previous best
    ## Run 152 stress 0.07629235 
    ## Run 153 stress 0.07365784 
    ## ... Procrustes: rmse 5.788147e-05  max resid 0.0001360831 
    ## ... Similar to previous best
    ## Run 154 stress 0.08000468 
    ## Run 155 stress 0.07378228 
    ## ... Procrustes: rmse 0.01743825  max resid 0.05099678 
    ## Run 156 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744582  max resid 0.05101991 
    ## Run 157 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174419  max resid 0.05104067 
    ## Run 158 stress 0.07365785 
    ## ... Procrustes: rmse 8.023787e-05  max resid 0.000188925 
    ## ... Similar to previous best
    ## Run 159 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001432345  max resid 0.000332344 
    ## ... Similar to previous best
    ## Run 160 stress 0.07365785 
    ## ... Procrustes: rmse 6.918426e-05  max resid 0.0001629268 
    ## ... Similar to previous best
    ## Run 161 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744372  max resid 0.05105557 
    ## Run 162 stress 0.07629236 
    ## Run 163 stress 0.0823376 
    ## Run 164 stress 0.07629247 
    ## Run 165 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744413  max resid 0.05104914 
    ## Run 166 stress 0.07365791 
    ## ... Procrustes: rmse 0.000152462  max resid 0.0003562909 
    ## ... Similar to previous best
    ## Run 167 stress 0.2574761 
    ## Run 168 stress 0.07365784 
    ## ... Procrustes: rmse 5.232658e-05  max resid 0.0001216692 
    ## ... Similar to previous best
    ## Run 169 stress 0.07732928 
    ## Run 170 stress 0.0828805 
    ## Run 171 stress 0.07732929 
    ## Run 172 stress 0.08000467 
    ## Run 173 stress 0.07732937 
    ## Run 174 stress 0.07629234 
    ## Run 175 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001496182  max resid 0.0003486595 
    ## ... Similar to previous best
    ## Run 176 stress 0.08000469 
    ## Run 177 stress 0.07365783 
    ## ... Procrustes: rmse 3.277953e-05  max resid 7.654784e-05 
    ## ... Similar to previous best
    ## Run 178 stress 0.07629236 
    ## Run 179 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174433  max resid 0.05103536 
    ## Run 180 stress 0.07365783 
    ## ... Procrustes: rmse 9.045994e-06  max resid 1.947154e-05 
    ## ... Similar to previous best
    ## Run 181 stress 0.07629233 
    ## Run 182 stress 0.07732933 
    ## Run 183 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744403  max resid 0.05104827 
    ## Run 184 stress 0.07629248 
    ## Run 185 stress 0.07629238 
    ## Run 186 stress 0.07629239 
    ## Run 187 stress 0.07629234 
    ## Run 188 stress 0.08000469 
    ## Run 189 stress 0.2559407 
    ## Run 190 stress 0.07629236 
    ## Run 191 stress 0.07365786 
    ## ... Procrustes: rmse 2.764379e-05  max resid 5.559269e-05 
    ## ... Similar to previous best
    ## Run 192 stress 0.07629241 
    ## Run 193 stress 0.07629234 
    ## Run 194 stress 0.07629233 
    ## Run 195 stress 0.07629242 
    ## Run 196 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744266  max resid 0.05103419 
    ## Run 197 stress 0.08288052 
    ## Run 198 stress 0.08288063 
    ## Run 199 stress 0.07365784 
    ## ... Procrustes: rmse 1.370168e-05  max resid 2.426574e-05 
    ## ... Similar to previous best
    ## Run 200 stress 0.07629247 
    ## Run 201 stress 0.3577963 
    ## Run 202 stress 0.0773294 
    ## Run 203 stress 0.08233754 
    ## Run 204 stress 0.08000468 
    ## Run 205 stress 0.07365785 
    ## ... Procrustes: rmse 7.215865e-05  max resid 0.0001716275 
    ## ... Similar to previous best
    ## Run 206 stress 0.07365785 
    ## ... Procrustes: rmse 8.071421e-05  max resid 0.0001915218 
    ## ... Similar to previous best
    ## Run 207 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744271  max resid 0.05102303 
    ## Run 208 stress 0.07629242 
    ## Run 209 stress 0.08000469 
    ## Run 210 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744489  max resid 0.05104834 
    ## Run 211 stress 0.07365783 
    ## ... Procrustes: rmse 2.713581e-05  max resid 6.318472e-05 
    ## ... Similar to previous best
    ## Run 212 stress 0.08000468 
    ## Run 213 stress 0.07629234 
    ## Run 214 stress 0.08000476 
    ## Run 215 stress 0.07732938 
    ## Run 216 stress 0.08000467 
    ## Run 217 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744444  max resid 0.05105755 
    ## Run 218 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001277503  max resid 0.0002975519 
    ## ... Similar to previous best
    ## Run 219 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001188746  max resid 0.0002828154 
    ## ... Similar to previous best
    ## Run 220 stress 0.07365784 
    ## ... Procrustes: rmse 5.095885e-05  max resid 0.0001204786 
    ## ... Similar to previous best
    ## Run 221 stress 0.07629253 
    ## Run 222 stress 0.08000467 
    ## Run 223 stress 0.07629234 
    ## Run 224 stress 0.08000468 
    ## Run 225 stress 0.0800047 
    ## Run 226 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743711  max resid 0.05104951 
    ## Run 227 stress 0.07365783 
    ## ... Procrustes: rmse 1.183266e-05  max resid 2.61661e-05 
    ## ... Similar to previous best
    ## Run 228 stress 0.08000471 
    ## Run 229 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744373  max resid 0.05104619 
    ## Run 230 stress 0.08000468 
    ## Run 231 stress 0.07378229 
    ## ... Procrustes: rmse 0.01743245  max resid 0.05097182 
    ## Run 232 stress 0.07629234 
    ## Run 233 stress 0.08000467 
    ## Run 234 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744363  max resid 0.05103324 
    ## Run 235 stress 0.08000468 
    ## Run 236 stress 0.07629234 
    ## Run 237 stress 0.08000467 
    ## Run 238 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744349  max resid 0.05103785 
    ## Run 239 stress 0.0800047 
    ## Run 240 stress 0.07365794 
    ## ... Procrustes: rmse 0.0001784096  max resid 0.0004197603 
    ## ... Similar to previous best
    ## Run 241 stress 0.07365783 
    ## ... Procrustes: rmse 1.304822e-05  max resid 2.885083e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744547  max resid 0.05105859 
    ## Run 243 stress 0.08288047 
    ## Run 244 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744359  max resid 0.0510285 
    ## Run 245 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001005261  max resid 0.0002380282 
    ## ... Similar to previous best
    ## Run 246 stress 0.08000468 
    ## Run 247 stress 0.0773294 
    ## Run 248 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744239  max resid 0.05105614 
    ## Run 249 stress 0.07365783 
    ## ... Procrustes: rmse 8.758871e-06  max resid 1.726614e-05 
    ## ... Similar to previous best
    ## Run 250 stress 0.07629235 
    ## Run 251 stress 0.07629233 
    ## Run 252 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744467  max resid 0.05104635 
    ## Run 253 stress 0.08233758 
    ## Run 254 stress 0.07365783 
    ## ... Procrustes: rmse 2.897218e-06  max resid 6.715731e-06 
    ## ... Similar to previous best
    ## Run 255 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743921  max resid 0.0510259 
    ## Run 256 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744657  max resid 0.05105965 
    ## Run 257 stress 0.261224 
    ## Run 258 stress 0.07365784 
    ## ... Procrustes: rmse 4.25826e-05  max resid 9.958792e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.08000468 
    ## Run 260 stress 0.07629244 
    ## Run 261 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174318  max resid 0.05099018 
    ## Run 262 stress 0.07365786 
    ## ... Procrustes: rmse 6.735332e-05  max resid 0.0001607433 
    ## ... Similar to previous best
    ## Run 263 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745072  max resid 0.05106671 
    ## Run 264 stress 0.07629238 
    ## Run 265 stress 0.07629242 
    ## Run 266 stress 0.07629243 
    ## Run 267 stress 0.07629236 
    ## Run 268 stress 0.07629233 
    ## Run 269 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744337  max resid 0.05103076 
    ## Run 270 stress 0.07629233 
    ## Run 271 stress 0.07629236 
    ## Run 272 stress 0.08000473 
    ## Run 273 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744018  max resid 0.05100779 
    ## Run 274 stress 0.07365783 
    ## ... Procrustes: rmse 1.702325e-05  max resid 3.345562e-05 
    ## ... Similar to previous best
    ## Run 275 stress 0.08000467 
    ## Run 276 stress 0.08288038 
    ## Run 277 stress 0.08000471 
    ## Run 278 stress 0.08000469 
    ## Run 279 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744195  max resid 0.05101123 
    ## Run 280 stress 0.0762924 
    ## Run 281 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744691  max resid 0.05106069 
    ## Run 282 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744449  max resid 0.05105751 
    ## Run 283 stress 0.07365787 
    ## ... Procrustes: rmse 8.934835e-05  max resid 0.0002036355 
    ## ... Similar to previous best
    ## Run 284 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744602  max resid 0.05104555 
    ## Run 285 stress 0.08000468 
    ## Run 286 stress 0.07629242 
    ## Run 287 stress 0.07365785 
    ## ... Procrustes: rmse 5.63776e-05  max resid 0.0001253733 
    ## ... Similar to previous best
    ## Run 288 stress 0.2840605 
    ## Run 289 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744788  max resid 0.0510645 
    ## Run 290 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744302  max resid 0.05103229 
    ## Run 291 stress 0.07629247 
    ## Run 292 stress 0.07629237 
    ## Run 293 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744234  max resid 0.05102415 
    ## Run 294 stress 0.07732947 
    ## Run 295 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744674  max resid 0.05105914 
    ## Run 296 stress 0.07629233 
    ## Run 297 stress 0.08233766 
    ## Run 298 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745196  max resid 0.05107862 
    ## Run 299 stress 0.07378232 
    ## ... Procrustes: rmse 0.01746635  max resid 0.05115924 
    ## Run 300 stress 0.07732926 
    ## Run 301 stress 0.07732943 
    ## Run 302 stress 0.07365783 
    ## ... Procrustes: rmse 2.985506e-05  max resid 7.015019e-05 
    ## ... Similar to previous best
    ## Run 303 stress 0.07629244 
    ## Run 304 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001097972  max resid 0.0002557698 
    ## ... Similar to previous best
    ## Run 305 stress 0.07365784 
    ## ... Procrustes: rmse 3.78061e-05  max resid 8.883341e-05 
    ## ... Similar to previous best
    ## Run 306 stress 0.07629234 
    ## Run 307 stress 0.07365786 
    ## ... Procrustes: rmse 8.692727e-05  max resid 0.0002009745 
    ## ... Similar to previous best
    ## Run 308 stress 0.07378244 
    ## ... Procrustes: rmse 0.01743484  max resid 0.05101411 
    ## Run 309 stress 0.07365784 
    ## ... Procrustes: rmse 3.632685e-05  max resid 8.550958e-05 
    ## ... Similar to previous best
    ## Run 310 stress 0.07629233 
    ## Run 311 stress 0.07365784 
    ## ... Procrustes: rmse 2.724754e-05  max resid 5.695577e-05 
    ## ... Similar to previous best
    ## Run 312 stress 0.08000468 
    ## Run 313 stress 0.315272 
    ## Run 314 stress 0.0737823 
    ## ... Procrustes: rmse 0.01746889  max resid 0.05117056 
    ## Run 315 stress 0.0737823 
    ## ... Procrustes: rmse 0.01746126  max resid 0.05114157 
    ## Run 316 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744108  max resid 0.05103191 
    ## Run 317 stress 0.07629236 
    ## Run 318 stress 0.07365785 
    ## ... Procrustes: rmse 6.834921e-05  max resid 0.0001589471 
    ## ... Similar to previous best
    ## Run 319 stress 0.07365785 
    ## ... Procrustes: rmse 7.551337e-05  max resid 0.0001770391 
    ## ... Similar to previous best
    ## Run 320 stress 0.07629233 
    ## Run 321 stress 0.07629235 
    ## Run 322 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745926  max resid 0.05111268 
    ## Run 323 stress 0.07629235 
    ## Run 324 stress 0.07629243 
    ## Run 325 stress 0.07629235 
    ## Run 326 stress 0.07365795 
    ## ... Procrustes: rmse 0.0001825972  max resid 0.0004321064 
    ## ... Similar to previous best
    ## Run 327 stress 0.3225161 
    ## Run 328 stress 0.08000468 
    ## Run 329 stress 0.07629235 
    ## Run 330 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744525  max resid 0.05103699 
    ## Run 331 stress 0.07629242 
    ## Run 332 stress 0.07629237 
    ## Run 333 stress 0.07365785 
    ## ... Procrustes: rmse 2.692137e-05  max resid 5.883135e-05 
    ## ... Similar to previous best
    ## Run 334 stress 0.07365784 
    ## ... Procrustes: rmse 5.561883e-05  max resid 0.000128027 
    ## ... Similar to previous best
    ## Run 335 stress 0.07365786 
    ## ... Procrustes: rmse 8.045377e-05  max resid 0.0001917999 
    ## ... Similar to previous best
    ## Run 336 stress 0.08000468 
    ## Run 337 stress 0.07365783 
    ## ... Procrustes: rmse 1.580963e-05  max resid 3.687398e-05 
    ## ... Similar to previous best
    ## Run 338 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744438  max resid 0.05103972 
    ## Run 339 stress 0.08000468 
    ## Run 340 stress 0.07629253 
    ## Run 341 stress 0.0762924 
    ## Run 342 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744401  max resid 0.05104554 
    ## Run 343 stress 0.07629234 
    ## Run 344 stress 0.07365783 
    ## ... Procrustes: rmse 3.425202e-05  max resid 8.053523e-05 
    ## ... Similar to previous best
    ## Run 345 stress 0.07365783 
    ## ... Procrustes: rmse 1.722389e-05  max resid 3.733253e-05 
    ## ... Similar to previous best
    ## Run 346 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744701  max resid 0.05103836 
    ## Run 347 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744254  max resid 0.05105698 
    ## Run 348 stress 0.08233758 
    ## Run 349 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744309  max resid 0.05104352 
    ## Run 350 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744654  max resid 0.05106403 
    ## Run 351 stress 0.08000468 
    ## Run 352 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744637  max resid 0.0510684 
    ## Run 353 stress 0.07365783 
    ## ... Procrustes: rmse 8.801219e-06  max resid 1.84996e-05 
    ## ... Similar to previous best
    ## Run 354 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744504  max resid 0.05105543 
    ## Run 355 stress 0.07629244 
    ## Run 356 stress 0.07365786 
    ## ... Procrustes: rmse 8.564071e-05  max resid 0.0002045109 
    ## ... Similar to previous best
    ## Run 357 stress 0.0773293 
    ## Run 358 stress 0.07629237 
    ## Run 359 stress 0.07365784 
    ## ... Procrustes: rmse 4.11379e-05  max resid 9.686183e-05 
    ## ... Similar to previous best
    ## Run 360 stress 0.07732948 
    ## Run 361 stress 0.07629237 
    ## Run 362 stress 0.07732932 
    ## Run 363 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744293  max resid 0.05103291 
    ## Run 364 stress 0.07365784 
    ## ... Procrustes: rmse 4.159289e-05  max resid 0.0001008 
    ## ... Similar to previous best
    ## Run 365 stress 0.07365784 
    ## ... Procrustes: rmse 3.207223e-05  max resid 6.705189e-05 
    ## ... Similar to previous best
    ## Run 366 stress 0.07629234 
    ## Run 367 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001473718  max resid 0.0003475594 
    ## ... Similar to previous best
    ## Run 368 stress 0.08000468 
    ## Run 369 stress 0.07629237 
    ## Run 370 stress 0.07629249 
    ## Run 371 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745129  max resid 0.05104907 
    ## Run 372 stress 0.07365784 
    ## ... Procrustes: rmse 6.085029e-05  max resid 0.0001398339 
    ## ... Similar to previous best
    ## Run 373 stress 0.08000467 
    ## Run 374 stress 0.08288034 
    ## Run 375 stress 0.07365785 
    ## ... Procrustes: rmse 6.832273e-05  max resid 0.000161171 
    ## ... Similar to previous best
    ## Run 376 stress 0.07365784 
    ## ... Procrustes: rmse 4.015846e-05  max resid 9.451094e-05 
    ## ... Similar to previous best
    ## Run 377 stress 0.07629235 
    ## Run 378 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743872  max resid 0.05101663 
    ## Run 379 stress 0.08233761 
    ## Run 380 stress 0.07365789 
    ## ... Procrustes: rmse 6.037902e-05  max resid 0.0001150277 
    ## ... Similar to previous best
    ## Run 381 stress 0.07629234 
    ## Run 382 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 3.791131e-06  max resid 7.869445e-06 
    ## ... Similar to previous best
    ## Run 383 stress 0.07365786 
    ## ... Procrustes: rmse 8.906505e-05  max resid 0.000209789 
    ## ... Similar to previous best
    ## Run 384 stress 0.08288042 
    ## Run 385 stress 0.08000467 
    ## Run 386 stress 0.08000467 
    ## Run 387 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744715  max resid 0.05104297 
    ## Run 388 stress 0.08000468 
    ## Run 389 stress 0.08233764 
    ## Run 390 stress 0.0762924 
    ## Run 391 stress 0.07629245 
    ## Run 392 stress 0.08000469 
    ## Run 393 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001193604  max resid 0.0002830731 
    ## ... Similar to previous best
    ## Run 394 stress 0.07365783 
    ## ... Procrustes: rmse 2.641186e-05  max resid 6.217644e-05 
    ## ... Similar to previous best
    ## Run 395 stress 0.07365794 
    ## ... Procrustes: rmse 0.0001709771  max resid 0.0004056055 
    ## ... Similar to previous best
    ## Run 396 stress 0.08000469 
    ## Run 397 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001451689  max resid 0.0003395526 
    ## ... Similar to previous best
    ## Run 398 stress 0.07365783 
    ## ... Procrustes: rmse 2.718651e-05  max resid 6.277024e-05 
    ## ... Similar to previous best
    ## Run 399 stress 0.08233757 
    ## Run 400 stress 0.07629242 
    ## Run 401 stress 0.08233771 
    ## Run 402 stress 0.07365798 
    ## ... Procrustes: rmse 0.0002085719  max resid 0.0004874446 
    ## ... Similar to previous best
    ## Run 403 stress 0.07378231 
    ## ... Procrustes: rmse 0.01746417  max resid 0.05113975 
    ## Run 404 stress 0.0762925 
    ## Run 405 stress 0.07365784 
    ## ... Procrustes: rmse 3.03249e-05  max resid 7.283079e-05 
    ## ... Similar to previous best
    ## Run 406 stress 0.0800047 
    ## Run 407 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744488  max resid 0.0510456 
    ## Run 408 stress 0.08000477 
    ## Run 409 stress 0.07365786 
    ## ... Procrustes: rmse 9.166e-05  max resid 0.0002146571 
    ## ... Similar to previous best
    ## Run 410 stress 0.07732943 
    ## Run 411 stress 0.08000468 
    ## Run 412 stress 0.07378227 
    ## ... Procrustes: rmse 0.017447  max resid 0.05105311 
    ## Run 413 stress 0.07365788 
    ## ... Procrustes: rmse 9.074001e-05  max resid 0.0002165846 
    ## ... Similar to previous best
    ## Run 414 stress 0.07629241 
    ## Run 415 stress 0.07629236 
    ## Run 416 stress 0.08000472 
    ## Run 417 stress 0.07629237 
    ## Run 418 stress 0.07629241 
    ## Run 419 stress 0.0800047 
    ## Run 420 stress 0.07365785 
    ## ... Procrustes: rmse 5.422939e-05  max resid 0.0001299356 
    ## ... Similar to previous best
    ## Run 421 stress 0.07732936 
    ## Run 422 stress 0.07365784 
    ## ... Procrustes: rmse 2.743658e-05  max resid 6.437667e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.07365784 
    ## ... Procrustes: rmse 3.925023e-05  max resid 8.789175e-05 
    ## ... Similar to previous best
    ## Run 424 stress 0.08288039 
    ## Run 425 stress 0.07629243 
    ## Run 426 stress 0.08288046 
    ## Run 427 stress 0.07629237 
    ## Run 428 stress 0.07629246 
    ## Run 429 stress 0.07365783 
    ## ... Procrustes: rmse 1.894243e-05  max resid 4.405142e-05 
    ## ... Similar to previous best
    ## Run 430 stress 0.07365789 
    ## ... Procrustes: rmse 7.624895e-05  max resid 0.0001608773 
    ## ... Similar to previous best
    ## Run 431 stress 0.07629239 
    ## Run 432 stress 0.07732938 
    ## Run 433 stress 0.07365784 
    ## ... Procrustes: rmse 3.914244e-05  max resid 9.191951e-05 
    ## ... Similar to previous best
    ## Run 434 stress 0.07629236 
    ## Run 435 stress 0.07365786 
    ## ... Procrustes: rmse 8.316238e-05  max resid 0.0001978822 
    ## ... Similar to previous best
    ## Run 436 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745589  max resid 0.05106007 
    ## Run 437 stress 0.07732936 
    ## Run 438 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744616  max resid 0.05106222 
    ## Run 439 stress 0.0800047 
    ## Run 440 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744385  max resid 0.05106537 
    ## Run 441 stress 0.07629243 
    ## Run 442 stress 0.07365783 
    ## ... Procrustes: rmse 3.562264e-06  max resid 8.779013e-06 
    ## ... Similar to previous best
    ## Run 443 stress 0.07365784 
    ## ... Procrustes: rmse 1.96734e-05  max resid 4.217154e-05 
    ## ... Similar to previous best
    ## Run 444 stress 0.07365783 
    ## ... Procrustes: rmse 1.56863e-05  max resid 3.581641e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.3411846 
    ## Run 446 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744738  max resid 0.05104692 
    ## Run 447 stress 0.07365783 
    ## ... Procrustes: rmse 1.969245e-05  max resid 4.589598e-05 
    ## ... Similar to previous best
    ## Run 448 stress 0.07629233 
    ## Run 449 stress 0.07629248 
    ## Run 450 stress 0.08000473 
    ## Run 451 stress 0.08000474 
    ## Run 452 stress 0.07365783 
    ## ... Procrustes: rmse 1.823819e-05  max resid 4.264133e-05 
    ## ... Similar to previous best
    ## Run 453 stress 0.07629245 
    ## Run 454 stress 0.07365783 
    ## ... Procrustes: rmse 2.090776e-05  max resid 4.883297e-05 
    ## ... Similar to previous best
    ## Run 455 stress 0.0773294 
    ## Run 456 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001260478  max resid 0.0002987281 
    ## ... Similar to previous best
    ## Run 457 stress 0.07629245 
    ## Run 458 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743839  max resid 0.05101806 
    ## Run 459 stress 0.07629238 
    ## Run 460 stress 0.07365786 
    ## ... Procrustes: rmse 4.300349e-05  max resid 9.842988e-05 
    ## ... Similar to previous best
    ## Run 461 stress 0.07365783 
    ## ... Procrustes: rmse 4.961834e-06  max resid 1.048533e-05 
    ## ... Similar to previous best
    ## Run 462 stress 0.07378226 
    ## ... Procrustes: rmse 0.017445  max resid 0.05103349 
    ## Run 463 stress 0.08000468 
    ## Run 464 stress 0.07629235 
    ## Run 465 stress 0.08288061 
    ## Run 466 stress 0.07629235 
    ## Run 467 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744723  max resid 0.05105192 
    ## Run 468 stress 0.07629246 
    ## Run 469 stress 0.07365783 
    ## ... Procrustes: rmse 2.062947e-05  max resid 4.663759e-05 
    ## ... Similar to previous best
    ## Run 470 stress 0.07732927 
    ## Run 471 stress 0.08000476 
    ## Run 472 stress 0.07629244 
    ## Run 473 stress 0.07365787 
    ## ... Procrustes: rmse 8.563378e-05  max resid 0.0001924506 
    ## ... Similar to previous best
    ## Run 474 stress 0.07365785 
    ## ... Procrustes: rmse 5.980127e-05  max resid 0.0001424952 
    ## ... Similar to previous best
    ## Run 475 stress 0.08000471 
    ## Run 476 stress 0.07629241 
    ## Run 477 stress 0.07629237 
    ## Run 478 stress 0.07629235 
    ## Run 479 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744664  max resid 0.05106292 
    ## Run 480 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744869  max resid 0.05106123 
    ## Run 481 stress 0.3495698 
    ## Run 482 stress 0.08000469 
    ## Run 483 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174468  max resid 0.05106063 
    ## Run 484 stress 0.07629245 
    ## Run 485 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001225434  max resid 0.000286903 
    ## ... Similar to previous best
    ## Run 486 stress 0.07629235 
    ## Run 487 stress 0.0762924 
    ## Run 488 stress 0.07629233 
    ## Run 489 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174483  max resid 0.05106435 
    ## Run 490 stress 0.07732936 
    ## Run 491 stress 0.07629236 
    ## Run 492 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745012  max resid 0.05106999 
    ## Run 493 stress 0.08233759 
    ## Run 494 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745834  max resid 0.05110928 
    ## Run 495 stress 0.08000472 
    ## Run 496 stress 0.07365784 
    ## ... Procrustes: rmse 6.729203e-05  max resid 0.0001578056 
    ## ... Similar to previous best
    ## Run 497 stress 0.08000474 
    ## Run 498 stress 0.08233761 
    ## Run 499 stress 0.07629235 
    ## Run 500 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744851  max resid 0.05106746 
    ## *** Best solution repeated 32 times

``` r
# Stratified lakes and ocean sites
SD_beta_geo_SO_NMDS <- metaMDS(SD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.07428316 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1089999  max resid 0.2637185 
    ## Run 2 stress 0.07428314 
    ## ... New best solution
    ## ... Procrustes: rmse 9.329515e-05  max resid 0.0001884952 
    ## ... Similar to previous best
    ## Run 3 stress 0.08340288 
    ## Run 4 stress 0.07970525 
    ## Run 5 stress 0.08340287 
    ## Run 6 stress 0.07250812 
    ## ... New best solution
    ## ... Procrustes: rmse 0.06048847  max resid 0.1271576 
    ## Run 7 stress 0.06978189 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04136793  max resid 0.1124467 
    ## Run 8 stress 0.06942777 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01321109  max resid 0.03347847 
    ## Run 9 stress 0.07428313 
    ## Run 10 stress 0.07250814 
    ## Run 11 stress 0.07970525 
    ## Run 12 stress 0.06942777 
    ## ... Procrustes: rmse 4.669623e-06  max resid 1.19201e-05 
    ## ... Similar to previous best
    ## Run 13 stress 0.08448438 
    ## Run 14 stress 0.07428313 
    ## Run 15 stress 0.07844926 
    ## Run 16 stress 0.07428314 
    ## Run 17 stress 0.07970526 
    ## Run 18 stress 0.069782 
    ## ... Procrustes: rmse 0.01327263  max resid 0.03336632 
    ## Run 19 stress 0.0844844 
    ## Run 20 stress 0.06978194 
    ## ... Procrustes: rmse 0.01317359  max resid 0.03315722 
    ## Run 21 stress 0.07250812 
    ## Run 22 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 4.275765e-05  max resid 0.0001104165 
    ## ... Similar to previous best
    ## Run 23 stress 0.07428313 
    ## Run 24 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325755  max resid 0.03333333 
    ## Run 25 stress 0.08340287 
    ## Run 26 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 1.049233e-05  max resid 2.53673e-05 
    ## ... Similar to previous best
    ## Run 27 stress 0.08340287 
    ## Run 28 stress 0.07970525 
    ## Run 29 stress 0.07250812 
    ## Run 30 stress 0.07970526 
    ## Run 31 stress 0.07250812 
    ## Run 32 stress 0.07428313 
    ## Run 33 stress 0.06978191 
    ## ... Procrustes: rmse 0.0131994  max resid 0.03321309 
    ## Run 34 stress 0.07970525 
    ## Run 35 stress 0.07428314 
    ## Run 36 stress 0.07428313 
    ## Run 37 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 3.280413e-06  max resid 7.772922e-06 
    ## ... Similar to previous best
    ## Run 38 stress 0.08340287 
    ## Run 39 stress 0.07250812 
    ## Run 40 stress 0.07428313 
    ## Run 41 stress 0.07428316 
    ## Run 42 stress 0.07428318 
    ## Run 43 stress 0.07428314 
    ## Run 44 stress 0.07428313 
    ## Run 45 stress 0.08340287 
    ## Run 46 stress 0.07428318 
    ## Run 47 stress 0.07970525 
    ## Run 48 stress 0.07970525 
    ## Run 49 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323472  max resid 0.03328794 
    ## Run 50 stress 0.07970526 
    ## Run 51 stress 0.07428313 
    ## Run 52 stress 0.07844949 
    ## Run 53 stress 0.07250812 
    ## Run 54 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326481  max resid 0.03334779 
    ## Run 55 stress 0.08448436 
    ## Run 56 stress 0.07970526 
    ## Run 57 stress 0.07970525 
    ## Run 58 stress 0.07250813 
    ## Run 59 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322631  max resid 0.03327055 
    ## Run 60 stress 0.08340287 
    ## Run 61 stress 0.07428313 
    ## Run 62 stress 0.07428319 
    ## Run 63 stress 0.07970526 
    ## Run 64 stress 0.07970526 
    ## Run 65 stress 0.06942776 
    ## ... Procrustes: rmse 1.391927e-05  max resid 3.601293e-05 
    ## ... Similar to previous best
    ## Run 66 stress 0.08340288 
    ## Run 67 stress 0.07428319 
    ## Run 68 stress 0.06942776 
    ## ... Procrustes: rmse 1.131323e-05  max resid 1.94564e-05 
    ## ... Similar to previous best
    ## Run 69 stress 0.07428316 
    ## Run 70 stress 0.07428313 
    ## Run 71 stress 0.07844929 
    ## Run 72 stress 0.07428313 
    ## Run 73 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326669  max resid 0.03335458 
    ## Run 74 stress 0.07970525 
    ## Run 75 stress 0.08340295 
    ## Run 76 stress 0.06978198 
    ## ... Procrustes: rmse 0.0131698  max resid 0.03314987 
    ## Run 77 stress 0.06942777 
    ## ... Procrustes: rmse 4.548993e-05  max resid 0.0001170903 
    ## ... Similar to previous best
    ## Run 78 stress 0.07428313 
    ## Run 79 stress 0.07970527 
    ## Run 80 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325381  max resid 0.03332897 
    ## Run 81 stress 0.06942777 
    ## ... Procrustes: rmse 5.671335e-05  max resid 0.0001463086 
    ## ... Similar to previous best
    ## Run 82 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132378  max resid 0.03329499 
    ## Run 83 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326406  max resid 0.0333514 
    ## Run 84 stress 0.06978205 
    ## ... Procrustes: rmse 0.01328548  max resid 0.03338672 
    ## Run 85 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326505  max resid 0.03335416 
    ## Run 86 stress 0.08340292 
    ## Run 87 stress 0.07428314 
    ## Run 88 stress 0.07970525 
    ## Run 89 stress 0.06942777 
    ## ... Procrustes: rmse 3.934007e-05  max resid 0.0001008681 
    ## ... Similar to previous best
    ## Run 90 stress 0.07844917 
    ## Run 91 stress 0.07428316 
    ## Run 92 stress 0.08340289 
    ## Run 93 stress 0.08448436 
    ## Run 94 stress 0.06942776 
    ## ... Procrustes: rmse 8.418771e-06  max resid 2.206009e-05 
    ## ... Similar to previous best
    ## Run 95 stress 0.06942776 
    ## ... Procrustes: rmse 1.304829e-05  max resid 3.392417e-05 
    ## ... Similar to previous best
    ## Run 96 stress 0.06942776 
    ## ... Procrustes: rmse 1.544504e-05  max resid 4.001136e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324285  max resid 0.03330555 
    ## Run 98 stress 0.069782 
    ## ... Procrustes: rmse 0.01328159  max resid 0.03338481 
    ## Run 99 stress 0.07250815 
    ## Run 100 stress 0.07250812 
    ## Run 101 stress 0.08340288 
    ## Run 102 stress 0.06942776 
    ## ... Procrustes: rmse 1.106498e-06  max resid 2.744924e-06 
    ## ... Similar to previous best
    ## Run 103 stress 0.07970526 
    ## Run 104 stress 0.07970526 
    ## Run 105 stress 0.07970525 
    ## Run 106 stress 0.06942776 
    ## ... Procrustes: rmse 1.674956e-05  max resid 4.352869e-05 
    ## ... Similar to previous best
    ## Run 107 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326952  max resid 0.03336009 
    ## Run 108 stress 0.07428313 
    ## Run 109 stress 0.07970526 
    ## Run 110 stress 0.08340289 
    ## Run 111 stress 0.07970526 
    ## Run 112 stress 0.07250812 
    ## Run 113 stress 0.07250813 
    ## Run 114 stress 0.07428313 
    ## Run 115 stress 0.06942776 
    ## ... Procrustes: rmse 2.590177e-05  max resid 6.76671e-05 
    ## ... Similar to previous best
    ## Run 116 stress 0.07250816 
    ## Run 117 stress 0.07428313 
    ## Run 118 stress 0.08340287 
    ## Run 119 stress 0.06942776 
    ## ... Procrustes: rmse 1.301652e-05  max resid 3.579747e-05 
    ## ... Similar to previous best
    ## Run 120 stress 0.07428313 
    ## Run 121 stress 0.07250812 
    ## Run 122 stress 0.07970525 
    ## Run 123 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322831  max resid 0.03327386 
    ## Run 124 stress 0.07250813 
    ## Run 125 stress 0.06978196 
    ## ... Procrustes: rmse 0.01317518  max resid 0.03315925 
    ## Run 126 stress 0.08448435 
    ## Run 127 stress 0.07970525 
    ## Run 128 stress 0.06942776 
    ## ... Procrustes: rmse 1.003012e-05  max resid 2.636383e-05 
    ## ... Similar to previous best
    ## Run 129 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323131  max resid 0.03328152 
    ## Run 130 stress 0.07428315 
    ## Run 131 stress 0.06978206 
    ## ... Procrustes: rmse 0.01328792  max resid 0.03339123 
    ## Run 132 stress 0.07970525 
    ## Run 133 stress 0.07250812 
    ## Run 134 stress 0.06942776 
    ## ... Procrustes: rmse 8.10019e-06  max resid 2.074245e-05 
    ## ... Similar to previous best
    ## Run 135 stress 0.07250812 
    ## Run 136 stress 0.07428313 
    ## Run 137 stress 0.07428313 
    ## Run 138 stress 0.07250812 
    ## Run 139 stress 0.07250812 
    ## Run 140 stress 0.07844936 
    ## Run 141 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326719  max resid 0.03335783 
    ## Run 142 stress 0.06942776 
    ## ... Procrustes: rmse 1.31372e-05  max resid 3.434272e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.07250812 
    ## Run 144 stress 0.06942778 
    ## ... Procrustes: rmse 2.399292e-05  max resid 6.847315e-05 
    ## ... Similar to previous best
    ## Run 145 stress 0.08340297 
    ## Run 146 stress 0.07428322 
    ## Run 147 stress 0.07428314 
    ## Run 148 stress 0.08340296 
    ## Run 149 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324567  max resid 0.03330977 
    ## Run 150 stress 0.08340287 
    ## Run 151 stress 0.07250812 
    ## Run 152 stress 0.07970525 
    ## Run 153 stress 0.07428314 
    ## Run 154 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327816  max resid 0.03337901 
    ## Run 155 stress 0.07250812 
    ## Run 156 stress 0.07428314 
    ## Run 157 stress 0.07428313 
    ## Run 158 stress 0.07970525 
    ## Run 159 stress 0.08448437 
    ## Run 160 stress 0.07970526 
    ## Run 161 stress 0.0844846 
    ## Run 162 stress 0.07428314 
    ## Run 163 stress 0.07428313 
    ## Run 164 stress 0.07970525 
    ## Run 165 stress 0.07970526 
    ## Run 166 stress 0.07428314 
    ## Run 167 stress 0.07428313 
    ## Run 168 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 4.022397e-06  max resid 9.103205e-06 
    ## ... Similar to previous best
    ## Run 169 stress 0.07970525 
    ## Run 170 stress 0.07428313 
    ## Run 171 stress 0.06978196 
    ## ... Procrustes: rmse 0.01327009  max resid 0.03336269 
    ## Run 172 stress 0.07970525 
    ## Run 173 stress 0.06978193 
    ## ... Procrustes: rmse 0.0132593  max resid 0.03333951 
    ## Run 174 stress 0.07428314 
    ## Run 175 stress 0.07428312 
    ## Run 176 stress 0.07970525 
    ## Run 177 stress 0.07970525 
    ## Run 178 stress 0.07844957 
    ## Run 179 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 4.561308e-06  max resid 1.056014e-05 
    ## ... Similar to previous best
    ## Run 180 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319207  max resid 0.03319427 
    ## Run 181 stress 0.07250812 
    ## Run 182 stress 0.07970525 
    ## Run 183 stress 0.06942776 
    ## ... Procrustes: rmse 2.323878e-05  max resid 6.012832e-05 
    ## ... Similar to previous best
    ## Run 184 stress 0.07250813 
    ## Run 185 stress 0.07970525 
    ## Run 186 stress 0.06978193 
    ## ... Procrustes: rmse 0.0132574  max resid 0.03333142 
    ## Run 187 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318824  max resid 0.03318589 
    ## Run 188 stress 0.07250812 
    ## Run 189 stress 0.08448437 
    ## Run 190 stress 0.0844845 
    ## Run 191 stress 0.07428312 
    ## Run 192 stress 0.07428319 
    ## Run 193 stress 0.07428313 
    ## Run 194 stress 0.06942777 
    ## ... Procrustes: rmse 5.350526e-05  max resid 0.0001386241 
    ## ... Similar to previous best
    ## Run 195 stress 0.07970525 
    ## Run 196 stress 0.08340287 
    ## Run 197 stress 0.069782 
    ## ... Procrustes: rmse 0.01328191  max resid 0.03338099 
    ## Run 198 stress 0.07970525 
    ## Run 199 stress 0.07250812 
    ## Run 200 stress 0.07428316 
    ## Run 201 stress 0.06978201 
    ## ... Procrustes: rmse 0.01327807  max resid 0.03336893 
    ## Run 202 stress 0.07250812 
    ## Run 203 stress 0.07428317 
    ## Run 204 stress 0.0834029 
    ## Run 205 stress 0.07428313 
    ## Run 206 stress 0.07970526 
    ## Run 207 stress 0.07428313 
    ## Run 208 stress 0.07250812 
    ## Run 209 stress 0.06942777 
    ## ... Procrustes: rmse 5.734506e-05  max resid 0.0001478451 
    ## ... Similar to previous best
    ## Run 210 stress 0.07970525 
    ## Run 211 stress 0.07970525 
    ## Run 212 stress 0.07428313 
    ## Run 213 stress 0.06942778 
    ## ... Procrustes: rmse 6.069004e-05  max resid 0.0001564412 
    ## ... Similar to previous best
    ## Run 214 stress 0.08340287 
    ## Run 215 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327619  max resid 0.0333694 
    ## Run 216 stress 0.07970525 
    ## Run 217 stress 0.07250814 
    ## Run 218 stress 0.0834029 
    ## Run 219 stress 0.06942777 
    ## ... Procrustes: rmse 3.817915e-05  max resid 9.865026e-05 
    ## ... Similar to previous best
    ## Run 220 stress 0.07428314 
    ## Run 221 stress 0.07250812 
    ## Run 222 stress 0.07428318 
    ## Run 223 stress 0.06942778 
    ## ... Procrustes: rmse 5.523197e-05  max resid 0.0001420035 
    ## ... Similar to previous best
    ## Run 224 stress 0.06978195 
    ## ... Procrustes: rmse 0.01318246  max resid 0.03316995 
    ## Run 225 stress 0.07250817 
    ## Run 226 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324468  max resid 0.03330523 
    ## Run 227 stress 0.08340291 
    ## Run 228 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319214  max resid 0.03319275 
    ## Run 229 stress 0.07844922 
    ## Run 230 stress 0.07970525 
    ## Run 231 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326744  max resid 0.03335396 
    ## Run 232 stress 0.07428312 
    ## Run 233 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325668  max resid 0.03333008 
    ## Run 234 stress 0.07970525 
    ## Run 235 stress 0.08340292 
    ## Run 236 stress 0.07428323 
    ## Run 237 stress 0.07844945 
    ## Run 238 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323188  max resid 0.03327615 
    ## Run 239 stress 0.06942776 
    ## ... Procrustes: rmse 2.083518e-05  max resid 5.263107e-05 
    ## ... Similar to previous best
    ## Run 240 stress 0.07428313 
    ## Run 241 stress 0.07970526 
    ## Run 242 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321384  max resid 0.03323982 
    ## Run 243 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324523  max resid 0.03330703 
    ## Run 244 stress 0.07428317 
    ## Run 245 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325311  max resid 0.03332346 
    ## Run 246 stress 0.06942776 
    ## ... Procrustes: rmse 1.608738e-05  max resid 4.111919e-05 
    ## ... Similar to previous best
    ## Run 247 stress 0.07428313 
    ## Run 248 stress 0.06942776 
    ## ... Procrustes: rmse 6.152302e-06  max resid 1.671823e-05 
    ## ... Similar to previous best
    ## Run 249 stress 0.06978197 
    ## ... Procrustes: rmse 0.01324545  max resid 0.03331273 
    ## Run 250 stress 0.07250812 
    ## Run 251 stress 0.08448454 
    ## Run 252 stress 0.06942776 
    ## ... Procrustes: rmse 2.362366e-05  max resid 6.12636e-05 
    ## ... Similar to previous best
    ## Run 253 stress 0.08340292 
    ## Run 254 stress 0.06942778 
    ## ... Procrustes: rmse 5.658804e-05  max resid 0.0001464579 
    ## ... Similar to previous best
    ## Run 255 stress 0.07428317 
    ## Run 256 stress 0.0834029 
    ## Run 257 stress 0.07970525 
    ## Run 258 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325717  max resid 0.03333094 
    ## Run 259 stress 0.07970526 
    ## Run 260 stress 0.07250812 
    ## Run 261 stress 0.07428313 
    ## Run 262 stress 0.07428315 
    ## Run 263 stress 0.07250812 
    ## Run 264 stress 0.08340287 
    ## Run 265 stress 0.07970526 
    ## Run 266 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321873  max resid 0.03325054 
    ## Run 267 stress 0.069782 
    ## ... Procrustes: rmse 0.01327814  max resid 0.03337189 
    ## Run 268 stress 0.07428317 
    ## Run 269 stress 0.07428315 
    ## Run 270 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319939  max resid 0.03320779 
    ## Run 271 stress 0.06942777 
    ## ... Procrustes: rmse 2.762285e-05  max resid 7.418893e-05 
    ## ... Similar to previous best
    ## Run 272 stress 0.08340287 
    ## Run 273 stress 0.07970525 
    ## Run 274 stress 0.08340286 
    ## Run 275 stress 0.07428317 
    ## Run 276 stress 0.07428313 
    ## Run 277 stress 0.07970525 
    ## Run 278 stress 0.07428312 
    ## Run 279 stress 0.07250813 
    ## Run 280 stress 0.07428313 
    ## Run 281 stress 0.07428321 
    ## Run 282 stress 0.07428313 
    ## Run 283 stress 0.07970525 
    ## Run 284 stress 0.07428314 
    ## Run 285 stress 0.07428317 
    ## Run 286 stress 0.07970526 
    ## Run 287 stress 0.07428313 
    ## Run 288 stress 0.07428313 
    ## Run 289 stress 0.07428313 
    ## Run 290 stress 0.07970525 
    ## Run 291 stress 0.07428315 
    ## Run 292 stress 0.07428319 
    ## Run 293 stress 0.08340287 
    ## Run 294 stress 0.06942777 
    ## ... Procrustes: rmse 4.062995e-05  max resid 0.0001060898 
    ## ... Similar to previous best
    ## Run 295 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323566  max resid 0.03328679 
    ## Run 296 stress 0.08448441 
    ## Run 297 stress 0.08340287 
    ## Run 298 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325514  max resid 0.03332264 
    ## Run 299 stress 0.07970525 
    ## Run 300 stress 0.06942776 
    ## ... Procrustes: rmse 7.054468e-06  max resid 1.742835e-05 
    ## ... Similar to previous best
    ## Run 301 stress 0.08340286 
    ## Run 302 stress 0.07250816 
    ## Run 303 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319505  max resid 0.03319791 
    ## Run 304 stress 0.08340289 
    ## Run 305 stress 0.07428313 
    ## Run 306 stress 0.06942776 
    ## ... Procrustes: rmse 5.075539e-06  max resid 1.20423e-05 
    ## ... Similar to previous best
    ## Run 307 stress 0.08340291 
    ## Run 308 stress 0.07250814 
    ## Run 309 stress 0.07428321 
    ## Run 310 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326219  max resid 0.03334123 
    ## Run 311 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322642  max resid 0.03326699 
    ## Run 312 stress 0.06942776 
    ## ... Procrustes: rmse 6.125965e-06  max resid 1.510512e-05 
    ## ... Similar to previous best
    ## Run 313 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321945  max resid 0.03325157 
    ## Run 314 stress 0.07428314 
    ## Run 315 stress 0.07250813 
    ## Run 316 stress 0.08340289 
    ## Run 317 stress 0.06942777 
    ## ... Procrustes: rmse 5.553508e-05  max resid 0.0001431398 
    ## ... Similar to previous best
    ## Run 318 stress 0.08340297 
    ## Run 319 stress 0.08340296 
    ## Run 320 stress 0.06942776 
    ## ... Procrustes: rmse 5.449535e-06  max resid 1.344347e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324543  max resid 0.03330734 
    ## Run 322 stress 0.07428314 
    ## Run 323 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321281  max resid 0.03323745 
    ## Run 324 stress 0.0834029 
    ## Run 325 stress 0.08340288 
    ## Run 326 stress 0.07428316 
    ## Run 327 stress 0.07428314 
    ## Run 328 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326295  max resid 0.03334322 
    ## Run 329 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324068  max resid 0.03329591 
    ## Run 330 stress 0.07250814 
    ## Run 331 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325876  max resid 0.03333549 
    ## Run 332 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324457  max resid 0.03330272 
    ## Run 333 stress 0.06942778 
    ## ... Procrustes: rmse 5.90335e-05  max resid 0.0001520192 
    ## ... Similar to previous best
    ## Run 334 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323238  max resid 0.03327849 
    ## Run 335 stress 0.08340294 
    ## Run 336 stress 0.07844919 
    ## Run 337 stress 0.06942778 
    ## ... Procrustes: rmse 6.178452e-05  max resid 0.0001590622 
    ## ... Similar to previous best
    ## Run 338 stress 0.08340291 
    ## Run 339 stress 0.07970526 
    ## Run 340 stress 0.0844845 
    ## Run 341 stress 0.07250813 
    ## Run 342 stress 0.07970525 
    ## Run 343 stress 0.06942778 
    ## ... Procrustes: rmse 5.941345e-05  max resid 0.0001530194 
    ## ... Similar to previous best
    ## Run 344 stress 0.06942777 
    ## ... Procrustes: rmse 5.281051e-05  max resid 0.0001360458 
    ## ... Similar to previous best
    ## Run 345 stress 0.08340287 
    ## Run 346 stress 0.06942776 
    ## ... Procrustes: rmse 6.290705e-06  max resid 1.584503e-05 
    ## ... Similar to previous best
    ## Run 347 stress 0.07428314 
    ## Run 348 stress 0.07428315 
    ## Run 349 stress 0.08340287 
    ## Run 350 stress 0.08340287 
    ## Run 351 stress 0.07428314 
    ## Run 352 stress 0.07250813 
    ## Run 353 stress 0.07428315 
    ## Run 354 stress 0.06942776 
    ## ... Procrustes: rmse 1.490412e-05  max resid 3.846394e-05 
    ## ... Similar to previous best
    ## Run 355 stress 0.07250812 
    ## Run 356 stress 0.07428315 
    ## Run 357 stress 0.07428315 
    ## Run 358 stress 0.07970526 
    ## Run 359 stress 0.07428313 
    ## Run 360 stress 0.08340287 
    ## Run 361 stress 0.07250812 
    ## Run 362 stress 0.07970525 
    ## Run 363 stress 0.07970525 
    ## Run 364 stress 0.08448436 
    ## Run 365 stress 0.07250812 
    ## Run 366 stress 0.0844844 
    ## Run 367 stress 0.08340288 
    ## Run 368 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325168  max resid 0.03332135 
    ## Run 369 stress 0.07250812 
    ## Run 370 stress 0.07250813 
    ## Run 371 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321907  max resid 0.03324373 
    ## Run 372 stress 0.07844943 
    ## Run 373 stress 0.07428313 
    ## Run 374 stress 0.07970526 
    ## Run 375 stress 0.06942777 
    ## ... Procrustes: rmse 1.713502e-05  max resid 4.788023e-05 
    ## ... Similar to previous best
    ## Run 376 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324895  max resid 0.03331664 
    ## Run 377 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322395  max resid 0.03326186 
    ## Run 378 stress 0.07428313 
    ## Run 379 stress 0.08340296 
    ## Run 380 stress 0.08340291 
    ## Run 381 stress 0.07970525 
    ## Run 382 stress 0.08340287 
    ## Run 383 stress 0.08340286 
    ## Run 384 stress 0.08340292 
    ## Run 385 stress 0.07250814 
    ## Run 386 stress 0.07250812 
    ## Run 387 stress 0.08340297 
    ## Run 388 stress 0.0834029 
    ## Run 389 stress 0.08448447 
    ## Run 390 stress 0.08340287 
    ## Run 391 stress 0.07428313 
    ## Run 392 stress 0.08340287 
    ## Run 393 stress 0.0834029 
    ## Run 394 stress 0.07844914 
    ## Run 395 stress 0.08340295 
    ## Run 396 stress 0.07428313 
    ## Run 397 stress 0.07428316 
    ## Run 398 stress 0.08448441 
    ## Run 399 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325392  max resid 0.03332575 
    ## Run 400 stress 0.07428313 
    ## Run 401 stress 0.06942776 
    ## ... Procrustes: rmse 7.205435e-06  max resid 1.518821e-05 
    ## ... Similar to previous best
    ## Run 402 stress 0.07250812 
    ## Run 403 stress 0.07428313 
    ## Run 404 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327368  max resid 0.03336194 
    ## Run 405 stress 0.07428313 
    ## Run 406 stress 0.07970525 
    ## Run 407 stress 0.07970525 
    ## Run 408 stress 0.07970525 
    ## Run 409 stress 0.07428313 
    ## Run 410 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322782  max resid 0.03326799 
    ## Run 411 stress 0.06942776 
    ## ... Procrustes: rmse 5.958854e-06  max resid 1.489072e-05 
    ## ... Similar to previous best
    ## Run 412 stress 0.07428314 
    ## Run 413 stress 0.06942776 
    ## ... Procrustes: rmse 6.964478e-06  max resid 1.76369e-05 
    ## ... Similar to previous best
    ## Run 414 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319558  max resid 0.03319696 
    ## Run 415 stress 0.07970525 
    ## Run 416 stress 0.07428314 
    ## Run 417 stress 0.07970525 
    ## Run 418 stress 0.07428313 
    ## Run 419 stress 0.06942776 
    ## ... Procrustes: rmse 2.242549e-05  max resid 5.980338e-05 
    ## ... Similar to previous best
    ## Run 420 stress 0.06978194 
    ## ... Procrustes: rmse 0.01318659  max resid 0.03317961 
    ## Run 421 stress 0.06942777 
    ## ... Procrustes: rmse 4.659838e-05  max resid 0.0001202412 
    ## ... Similar to previous best
    ## Run 422 stress 0.06942776 
    ## ... Procrustes: rmse 1.303093e-05  max resid 3.338773e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.06942776 
    ## ... Procrustes: rmse 1.474698e-05  max resid 3.675806e-05 
    ## ... Similar to previous best
    ## Run 424 stress 0.07428313 
    ## Run 425 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323645  max resid 0.03328908 
    ## Run 426 stress 0.07250813 
    ## Run 427 stress 0.08340287 
    ## Run 428 stress 0.07250814 
    ## Run 429 stress 0.07428313 
    ## Run 430 stress 0.07428313 
    ## Run 431 stress 0.07428313 
    ## Run 432 stress 0.08340289 
    ## Run 433 stress 0.07250814 
    ## Run 434 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323644  max resid 0.03328809 
    ## Run 435 stress 0.08340287 
    ## Run 436 stress 0.08340295 
    ## Run 437 stress 0.06942776 
    ## ... Procrustes: rmse 1.332827e-05  max resid 3.016889e-05 
    ## ... Similar to previous best
    ## Run 438 stress 0.08340288 
    ## Run 439 stress 0.07428314 
    ## Run 440 stress 0.06942777 
    ## ... Procrustes: rmse 4.306707e-05  max resid 0.0001105304 
    ## ... Similar to previous best
    ## Run 441 stress 0.07250812 
    ## Run 442 stress 0.07970525 
    ## Run 443 stress 0.06978198 
    ## ... Procrustes: rmse 0.0132762  max resid 0.03336847 
    ## Run 444 stress 0.07428314 
    ## Run 445 stress 0.07250812 
    ## Run 446 stress 0.07428313 
    ## Run 447 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324173  max resid 0.03330072 
    ## Run 448 stress 0.08340286 
    ## Run 449 stress 0.07970525 
    ## Run 450 stress 0.06942778 
    ## ... Procrustes: rmse 3.797576e-05  max resid 9.522516e-05 
    ## ... Similar to previous best
    ## Run 451 stress 0.07428313 
    ## Run 452 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321611  max resid 0.03324521 
    ## Run 453 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319803  max resid 0.03320818 
    ## Run 454 stress 0.0834029 
    ## Run 455 stress 0.06978194 
    ## ... Procrustes: rmse 0.01319414  max resid 0.03319703 
    ## Run 456 stress 0.07250814 
    ## Run 457 stress 0.069782 
    ## ... Procrustes: rmse 0.01328245  max resid 0.03338213 
    ## Run 458 stress 0.07428315 
    ## Run 459 stress 0.07970525 
    ## Run 460 stress 0.06942776 
    ## ... Procrustes: rmse 1.780829e-05  max resid 4.600842e-05 
    ## ... Similar to previous best
    ## Run 461 stress 0.0844844 
    ## Run 462 stress 0.06978198 
    ## ... Procrustes: rmse 0.01317546  max resid 0.03316086 
    ## Run 463 stress 0.08340288 
    ## Run 464 stress 0.07428313 
    ## Run 465 stress 0.07428313 
    ## Run 466 stress 0.08340289 
    ## Run 467 stress 0.07250812 
    ## Run 468 stress 0.07428313 
    ## Run 469 stress 0.07428319 
    ## Run 470 stress 0.07970525 
    ## Run 471 stress 0.06942777 
    ## ... Procrustes: rmse 3.602083e-05  max resid 9.291313e-05 
    ## ... Similar to previous best
    ## Run 472 stress 0.07250814 
    ## Run 473 stress 0.07970526 
    ## Run 474 stress 0.07970525 
    ## Run 475 stress 0.08340287 
    ## Run 476 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322645  max resid 0.03327101 
    ## Run 477 stress 0.07428313 
    ## Run 478 stress 0.08340286 
    ## Run 479 stress 0.07250814 
    ## Run 480 stress 0.07428316 
    ## Run 481 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.393226e-06  max resid 4.45052e-06 
    ## ... Similar to previous best
    ## Run 482 stress 0.08340287 
    ## Run 483 stress 0.08340287 
    ## Run 484 stress 0.07428313 
    ## Run 485 stress 0.07428315 
    ## Run 486 stress 0.07970525 
    ## Run 487 stress 0.08340288 
    ## Run 488 stress 0.07250812 
    ## Run 489 stress 0.07250812 
    ## Run 490 stress 0.08448449 
    ## Run 491 stress 0.06942776 
    ## ... Procrustes: rmse 1.171123e-05  max resid 3.081524e-05 
    ## ... Similar to previous best
    ## Run 492 stress 0.06942776 
    ## ... Procrustes: rmse 2.47585e-05  max resid 6.4223e-05 
    ## ... Similar to previous best
    ## Run 493 stress 0.08340292 
    ## Run 494 stress 0.07970525 
    ## Run 495 stress 0.08340287 
    ## Run 496 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325543  max resid 0.03332836 
    ## Run 497 stress 0.08340288 
    ## Run 498 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323832  max resid 0.03329446 
    ## Run 499 stress 0.07970526 
    ## Run 500 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323009  max resid 0.03327486 
    ## *** Best solution repeated 3 times

``` r
# Mixed lakes
SD_beta_geo_M_NMDS <- metaMDS(SD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 9.385805e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04251051  max resid 0.06411063 
    ## Run 2 stress 0.1491957 
    ## Run 3 stress 0.0003319497 
    ## ... Procrustes: rmse 0.01020165  max resid 0.01397215 
    ## Run 4 stress 9.49567e-05 
    ## ... Procrustes: rmse 0.003120758  max resid 0.00509169 
    ## ... Similar to previous best
    ## Run 5 stress 0.001373121 
    ## Run 6 stress 9.690307e-05 
    ## ... Procrustes: rmse 0.003082603  max resid 0.00542227 
    ## ... Similar to previous best
    ## Run 7 stress 9.223362e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.003073836  max resid 0.004162484 
    ## ... Similar to previous best
    ## Run 8 stress 9.669258e-05 
    ## ... Procrustes: rmse 0.0001246702  max resid 0.0002099595 
    ## ... Similar to previous best
    ## Run 9 stress 0.1491957 
    ## Run 10 stress 0.00132002 
    ## Run 11 stress 9.176299e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002088168  max resid 0.0004523973 
    ## ... Similar to previous best
    ## Run 12 stress 0.001208037 
    ## Run 13 stress 0.0003931159 
    ## ... Procrustes: rmse 0.01451227  max resid 0.01977848 
    ## Run 14 stress 9.447733e-05 
    ## ... Procrustes: rmse 0.0002116071  max resid 0.0004550833 
    ## ... Similar to previous best
    ## Run 15 stress 0.0004698312 
    ## ... Procrustes: rmse 0.02609407  max resid 0.03592343 
    ## Run 16 stress 0.2568281 
    ## Run 17 stress 0.1990774 
    ## Run 18 stress 0.001306122 
    ## Run 19 stress 0.2528294 
    ## Run 20 stress 9.961436e-05 
    ## ... Procrustes: rmse 0.0002060555  max resid 0.0004250465 
    ## ... Similar to previous best
    ## Run 21 stress 9.490936e-05 
    ## ... Procrustes: rmse 8.97548e-05  max resid 0.0002046831 
    ## ... Similar to previous best
    ## Run 22 stress 9.821137e-05 
    ## ... Procrustes: rmse 0.0002133024  max resid 0.0004505246 
    ## ... Similar to previous best
    ## Run 23 stress 0.001358855 
    ## Run 24 stress 0.0006302693 
    ## Run 25 stress 0.001332189 
    ## Run 26 stress 0.0004372863 
    ## ... Procrustes: rmse 0.01530967  max resid 0.02087891 
    ## Run 27 stress 0.1491957 
    ## Run 28 stress 0.001349012 
    ## Run 29 stress 0.001324817 
    ## Run 30 stress 9.909813e-05 
    ## ... Procrustes: rmse 9.383043e-05  max resid 0.0002188224 
    ## ... Similar to previous best
    ## Run 31 stress 0.001308593 
    ## Run 32 stress 0.0004908225 
    ## ... Procrustes: rmse 0.01622062  max resid 0.02213646 
    ## Run 33 stress 9.489152e-05 
    ## ... Procrustes: rmse 0.0002104259  max resid 0.0004535149 
    ## ... Similar to previous best
    ## Run 34 stress 9.863235e-05 
    ## ... Procrustes: rmse 0.0003170766  max resid 0.000672264 
    ## ... Similar to previous best
    ## Run 35 stress 9.574837e-05 
    ## ... Procrustes: rmse 9.548698e-05  max resid 0.0001502589 
    ## ... Similar to previous best
    ## Run 36 stress 0.2842598 
    ## Run 37 stress 0.001210504 
    ## Run 38 stress 9.860197e-05 
    ## ... Procrustes: rmse 9.724156e-05  max resid 0.0002200655 
    ## ... Similar to previous best
    ## Run 39 stress 9.875848e-05 
    ## ... Procrustes: rmse 1.891681e-05  max resid 3.240171e-05 
    ## ... Similar to previous best
    ## Run 40 stress 9.488212e-05 
    ## ... Procrustes: rmse 0.0002299037  max resid 0.0004953395 
    ## ... Similar to previous best
    ## Run 41 stress 0.000460053 
    ## ... Procrustes: rmse 0.01570378  max resid 0.02142279 
    ## Run 42 stress 0.00132613 
    ## Run 43 stress 9.086565e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001986677  max resid 0.0003167602 
    ## ... Similar to previous best
    ## Run 44 stress 7.357026e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001748988  max resid 0.000346801 
    ## ... Similar to previous best
    ## Run 45 stress 0.0001179727 
    ## ... Procrustes: rmse 0.007913707  max resid 0.01075708 
    ## Run 46 stress 0.1990774 
    ## Run 47 stress 9.51992e-05 
    ## ... Procrustes: rmse 0.0001803362  max resid 0.0002716879 
    ## ... Similar to previous best
    ## Run 48 stress 0.1990774 
    ## Run 49 stress 9.717845e-05 
    ## ... Procrustes: rmse 0.000181527  max resid 0.0003630912 
    ## ... Similar to previous best
    ## Run 50 stress 0.0001369948 
    ## ... Procrustes: rmse 0.01408809  max resid 0.01929771 
    ## Run 51 stress 9.927645e-05 
    ## ... Procrustes: rmse 0.0001017712  max resid 0.0002239208 
    ## ... Similar to previous best
    ## Run 52 stress 0.0005422299 
    ## ... Procrustes: rmse 0.02809098  max resid 0.03867536 
    ## Run 53 stress 0.2528294 
    ## Run 54 stress 9.932256e-05 
    ## ... Procrustes: rmse 0.0001815127  max resid 0.0004004804 
    ## ... Similar to previous best
    ## Run 55 stress 0.0007360964 
    ## Run 56 stress 0.001260072 
    ## Run 57 stress 0.001244756 
    ## Run 58 stress 0.2842805 
    ## Run 59 stress 9.222623e-05 
    ## ... Procrustes: rmse 3.662693e-05  max resid 6.008639e-05 
    ## ... Similar to previous best
    ## Run 60 stress 0.0008516668 
    ## Run 61 stress 0.0009295168 
    ## Run 62 stress 9.479536e-05 
    ## ... Procrustes: rmse 0.0001851181  max resid 0.0003732493 
    ## ... Similar to previous best
    ## Run 63 stress 0.001241933 
    ## Run 64 stress 0.0004634857 
    ## ... Procrustes: rmse 0.01573997  max resid 0.02155788 
    ## Run 65 stress 0.0002478628 
    ## ... Procrustes: rmse 0.01148327  max resid 0.01568336 
    ## Run 66 stress 0.1491957 
    ## Run 67 stress 5.672007e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002267248  max resid 0.0004763996 
    ## ... Similar to previous best
    ## Run 68 stress 0.0012102 
    ## Run 69 stress 9.7829e-05 
    ## ... Procrustes: rmse 0.0002610526  max resid 0.0005146319 
    ## ... Similar to previous best
    ## Run 70 stress 0.00133906 
    ## Run 71 stress 8.353678e-05 
    ## ... Procrustes: rmse 0.00103516  max resid 0.001874441 
    ## ... Similar to previous best
    ## Run 72 stress 9.990649e-05 
    ## ... Procrustes: rmse 0.0002772386  max resid 0.0005000099 
    ## ... Similar to previous best
    ## Run 73 stress 0.001303574 
    ## Run 74 stress 9.340008e-05 
    ## ... Procrustes: rmse 0.0002394539  max resid 0.0005169394 
    ## ... Similar to previous best
    ## Run 75 stress 9.730308e-05 
    ## ... Procrustes: rmse 0.0002322217  max resid 0.000486131 
    ## ... Similar to previous best
    ## Run 76 stress 0.001146706 
    ## Run 77 stress 0.00122739 
    ## Run 78 stress 8.832155e-05 
    ## ... Procrustes: rmse 0.000193774  max resid 0.0003844805 
    ## ... Similar to previous best
    ## Run 79 stress 0.0001412545 
    ## ... Procrustes: rmse 0.008676199  max resid 0.01245234 
    ## Run 80 stress 9.329317e-05 
    ## ... Procrustes: rmse 0.0002293829  max resid 0.0004959591 
    ## ... Similar to previous best
    ## Run 81 stress 0.2440272 
    ## Run 82 stress 0.001216468 
    ## Run 83 stress 0.2848435 
    ## Run 84 stress 0.0004929951 
    ## ... Procrustes: rmse 0.01624903  max resid 0.02290665 
    ## Run 85 stress 9.705635e-05 
    ## ... Procrustes: rmse 0.001496918  max resid 0.001976322 
    ## ... Similar to previous best
    ## Run 86 stress 6.895432e-05 
    ## ... Procrustes: rmse 0.0004474156  max resid 0.0009603524 
    ## ... Similar to previous best
    ## Run 87 stress 9.967737e-05 
    ## ... Procrustes: rmse 0.0001164749  max resid 0.0001799022 
    ## ... Similar to previous best
    ## Run 88 stress 0.000481163 
    ## ... Procrustes: rmse 0.01580127  max resid 0.02228823 
    ## Run 89 stress 8.446277e-05 
    ## ... Procrustes: rmse 0.002213521  max resid 0.003547885 
    ## ... Similar to previous best
    ## Run 90 stress 0.2854408 
    ## Run 91 stress 9.665577e-05 
    ## ... Procrustes: rmse 0.0002690381  max resid 0.0005099252 
    ## ... Similar to previous best
    ## Run 92 stress 0.2797323 
    ## Run 93 stress 0.0005135237 
    ## ... Procrustes: rmse 0.01657024  max resid 0.02334935 
    ## Run 94 stress 0.0012988 
    ## Run 95 stress 9.928761e-05 
    ## ... Procrustes: rmse 0.0002587931  max resid 0.000500383 
    ## ... Similar to previous best
    ## Run 96 stress 9.929414e-05 
    ## ... Procrustes: rmse 0.000283957  max resid 0.0005297436 
    ## ... Similar to previous best
    ## Run 97 stress 0.0005102782 
    ## ... Procrustes: rmse 0.01600291  max resid 0.02256728 
    ## Run 98 stress 0.001293281 
    ## Run 99 stress 0.001284777 
    ## Run 100 stress 0.177602 
    ## Run 101 stress 9.484744e-05 
    ## ... Procrustes: rmse 0.0002362001  max resid 0.0005089974 
    ## ... Similar to previous best
    ## Run 102 stress 9.359927e-05 
    ## ... Procrustes: rmse 0.0002171198  max resid 0.0004671866 
    ## ... Similar to previous best
    ## Run 103 stress 0.00118842 
    ## Run 104 stress 0.0004099574 
    ## ... Procrustes: rmse 0.01481321  max resid 0.0209246 
    ## Run 105 stress 0.001329496 
    ## Run 106 stress 0.1776022 
    ## Run 107 stress 0.2842805 
    ## Run 108 stress 9.957072e-05 
    ## ... Procrustes: rmse 0.0002418145  max resid 0.0004956105 
    ## ... Similar to previous best
    ## Run 109 stress 0.00106385 
    ## Run 110 stress 0.1491957 
    ## Run 111 stress 0.1491957 
    ## Run 112 stress 0.0005162695 
    ## ... Procrustes: rmse 0.01662962  max resid 0.02343103 
    ## Run 113 stress 9.738378e-05 
    ## ... Procrustes: rmse 0.0001354722  max resid 0.0002289578 
    ## ... Similar to previous best
    ## Run 114 stress 0.001200747 
    ## Run 115 stress 8.798062e-05 
    ## ... Procrustes: rmse 0.000159313  max resid 0.000310821 
    ## ... Similar to previous best
    ## Run 116 stress 0.1491957 
    ## Run 117 stress 0.001300159 
    ## Run 118 stress 0.1990774 
    ## Run 119 stress 0.001321767 
    ## Run 120 stress 0.001454811 
    ## Run 121 stress 8.67312e-05 
    ## ... Procrustes: rmse 0.0002858627  max resid 0.0004522097 
    ## ... Similar to previous best
    ## Run 122 stress 0.000315294 
    ## ... Procrustes: rmse 0.01296259  max resid 0.01837021 
    ## Run 123 stress 0.0003583224 
    ## ... Procrustes: rmse 0.01384846  max resid 0.01959276 
    ## Run 124 stress 8.599527e-05 
    ## ... Procrustes: rmse 0.0002317947  max resid 0.0005116045 
    ## ... Similar to previous best
    ## Run 125 stress 9.191119e-05 
    ## ... Procrustes: rmse 0.0002377482  max resid 0.0004813713 
    ## ... Similar to previous best
    ## Run 126 stress 5.617703e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003335025  max resid 0.000538198 
    ## ... Similar to previous best
    ## Run 127 stress 0.001268484 
    ## Run 128 stress 0.2848513 
    ## Run 129 stress 0.001387587 
    ## Run 130 stress 9.297539e-05 
    ## ... Procrustes: rmse 0.0001779253  max resid 0.0004342849 
    ## ... Similar to previous best
    ## Run 131 stress 0.00113688 
    ## Run 132 stress 0.001310392 
    ## Run 133 stress 6.104403e-05 
    ## ... Procrustes: rmse 0.0001787254  max resid 0.0002494663 
    ## ... Similar to previous best
    ## Run 134 stress 9.396359e-05 
    ## ... Procrustes: rmse 0.0001587533  max resid 0.0003593734 
    ## ... Similar to previous best
    ## Run 135 stress 0.0007070386 
    ## Run 136 stress 8.766957e-05 
    ## ... Procrustes: rmse 0.0001902469  max resid 0.0004550837 
    ## ... Similar to previous best
    ## Run 137 stress 0.0004678998 
    ## ... Procrustes: rmse 0.01569212  max resid 0.02161062 
    ## Run 138 stress 0.0007730686 
    ## Run 139 stress 0.0004662034 
    ## ... Procrustes: rmse 0.0156616  max resid 0.02156774 
    ## Run 140 stress 0.1491957 
    ## Run 141 stress 9.387368e-05 
    ## ... Procrustes: rmse 0.0002389433  max resid 0.0004902641 
    ## ... Similar to previous best
    ## Run 142 stress 0.001355298 
    ## Run 143 stress 9.609298e-05 
    ## ... Procrustes: rmse 0.000226626  max resid 0.0004929444 
    ## ... Similar to previous best
    ## Run 144 stress 0.2842805 
    ## Run 145 stress 0.0004565369 
    ## ... Procrustes: rmse 0.01524111  max resid 0.02098899 
    ## Run 146 stress 0.1776022 
    ## Run 147 stress 0.000259455 
    ## ... Procrustes: rmse 0.01944957  max resid 0.02719816 
    ## Run 148 stress 0.001210614 
    ## Run 149 stress 0.001205478 
    ## Run 150 stress 9.587371e-05 
    ## ... Procrustes: rmse 0.0002085588  max resid 0.0003814083 
    ## ... Similar to previous best
    ## Run 151 stress 0.001155166 
    ## Run 152 stress 0.0004994073 
    ## ... Procrustes: rmse 0.0269893  max resid 0.03763522 
    ## Run 153 stress 9.565599e-05 
    ## ... Procrustes: rmse 0.0001415239  max resid 0.0003326147 
    ## ... Similar to previous best
    ## Run 154 stress 8.326276e-05 
    ## ... Procrustes: rmse 0.002782694  max resid 0.004195673 
    ## ... Similar to previous best
    ## Run 155 stress 9.525942e-05 
    ## ... Procrustes: rmse 0.0002018869  max resid 0.0003283859 
    ## ... Similar to previous best
    ## Run 156 stress 0.001188131 
    ## Run 157 stress 0.001297982 
    ## Run 158 stress 0.001290382 
    ## Run 159 stress 0.0004712962 
    ## ... Procrustes: rmse 0.0156423  max resid 0.02154244 
    ## Run 160 stress 0.1491957 
    ## Run 161 stress 7.819835e-05 
    ## ... Procrustes: rmse 0.0001302598  max resid 0.0003025291 
    ## ... Similar to previous best
    ## Run 162 stress 0.001406065 
    ## Run 163 stress 9.434943e-05 
    ## ... Procrustes: rmse 0.000237267  max resid 0.0004846774 
    ## ... Similar to previous best
    ## Run 164 stress 7.88243e-05 
    ## ... Procrustes: rmse 0.0002919958  max resid 0.0004822814 
    ## ... Similar to previous best
    ## Run 165 stress 0.1491957 
    ## Run 166 stress 0.00110035 
    ## Run 167 stress 9.520794e-05 
    ## ... Procrustes: rmse 0.0001970993  max resid 0.0004776718 
    ## ... Similar to previous best
    ## Run 168 stress 0.001377163 
    ## Run 169 stress 8.803107e-05 
    ## ... Procrustes: rmse 0.0002232461  max resid 0.000462269 
    ## ... Similar to previous best
    ## Run 170 stress 0.001150793 
    ## Run 171 stress 0.1990774 
    ## Run 172 stress 0.1491957 
    ## Run 173 stress 0.001247728 
    ## Run 174 stress 0.1776017 
    ## Run 175 stress 0.1491957 
    ## Run 176 stress 0.1776012 
    ## Run 177 stress 9.545463e-05 
    ## ... Procrustes: rmse 0.0001541544  max resid 0.0002914637 
    ## ... Similar to previous best
    ## Run 178 stress 9.307071e-05 
    ## ... Procrustes: rmse 0.0001549695  max resid 0.0003058871 
    ## ... Similar to previous best
    ## Run 179 stress 0.001034401 
    ## Run 180 stress 9.257008e-05 
    ## ... Procrustes: rmse 0.0001459174  max resid 0.0002808805 
    ## ... Similar to previous best
    ## Run 181 stress 0.001185073 
    ## Run 182 stress 0.00046088 
    ## ... Procrustes: rmse 0.01557202  max resid 0.02144518 
    ## Run 183 stress 0.1491957 
    ## Run 184 stress 8.84517e-05 
    ## ... Procrustes: rmse 0.0001579739  max resid 0.0003183443 
    ## ... Similar to previous best
    ## Run 185 stress 9.281612e-05 
    ## ... Procrustes: rmse 0.0006523585  max resid 0.001248136 
    ## ... Similar to previous best
    ## Run 186 stress 0.0007958109 
    ## Run 187 stress 8.227875e-05 
    ## ... Procrustes: rmse 7.250418e-05  max resid 0.0001126807 
    ## ... Similar to previous best
    ## Run 188 stress 0.0003977921 
    ## ... Procrustes: rmse 0.01444092  max resid 0.019885 
    ## Run 189 stress 9.193983e-05 
    ## ... Procrustes: rmse 0.0002432543  max resid 0.0004794982 
    ## ... Similar to previous best
    ## Run 190 stress 0.00133066 
    ## Run 191 stress 0.0005103161 
    ## ... Procrustes: rmse 0.01638466  max resid 0.02256554 
    ## Run 192 stress 9.258637e-05 
    ## ... Procrustes: rmse 0.0002200177  max resid 0.0004938054 
    ## ... Similar to previous best
    ## Run 193 stress 0.0009902948 
    ## Run 194 stress 0.001370978 
    ## Run 195 stress 0.0005179828 
    ## ... Procrustes: rmse 0.02749171  max resid 0.03833201 
    ## Run 196 stress 0.000536931 
    ## ... Procrustes: rmse 0.01682076  max resid 0.02316961 
    ## Run 197 stress 0.001413316 
    ## Run 198 stress 9.775499e-05 
    ## ... Procrustes: rmse 0.0002454857  max resid 0.000513506 
    ## ... Similar to previous best
    ## Run 199 stress 0.1776019 
    ## Run 200 stress 9.9843e-05 
    ## ... Procrustes: rmse 0.0002324531  max resid 0.0005155131 
    ## ... Similar to previous best
    ## Run 201 stress 0.1990774 
    ## Run 202 stress 0.0005062258 
    ## ... Procrustes: rmse 0.01632367  max resid 0.0224842 
    ## Run 203 stress 9.762662e-05 
    ## ... Procrustes: rmse 0.0002456915  max resid 0.0004841834 
    ## ... Similar to previous best
    ## Run 204 stress 0.1990774 
    ## Run 205 stress 9.318025e-05 
    ## ... Procrustes: rmse 0.0002157804  max resid 0.0004250406 
    ## ... Similar to previous best
    ## Run 206 stress 0.0007058303 
    ## Run 207 stress 9.04057e-05 
    ## ... Procrustes: rmse 0.0002039285  max resid 0.000362601 
    ## ... Similar to previous best
    ## Run 208 stress 0.1990774 
    ## Run 209 stress 0.2440272 
    ## Run 210 stress 0.0005132578 
    ## ... Procrustes: rmse 0.01644214  max resid 0.02264526 
    ## Run 211 stress 9.825066e-05 
    ## ... Procrustes: rmse 0.0002308682  max resid 0.0005052311 
    ## ... Similar to previous best
    ## Run 212 stress 8.839938e-05 
    ## ... Procrustes: rmse 0.000558378  max resid 0.001101945 
    ## ... Similar to previous best
    ## Run 213 stress 0.001342882 
    ## Run 214 stress 9.310801e-05 
    ## ... Procrustes: rmse 0.0001543441  max resid 0.0003465228 
    ## ... Similar to previous best
    ## Run 215 stress 7.907285e-05 
    ## ... Procrustes: rmse 0.0001264791  max resid 0.0002989473 
    ## ... Similar to previous best
    ## Run 216 stress 9.910659e-05 
    ## ... Procrustes: rmse 0.0001707095  max resid 0.0003668809 
    ## ... Similar to previous best
    ## Run 217 stress 0.0004504926 
    ## ... Procrustes: rmse 0.01538474  max resid 0.02118678 
    ## Run 218 stress 0.0003009655 
    ## ... Procrustes: rmse 0.01255054  max resid 0.01727625 
    ## Run 219 stress 9.644249e-05 
    ## ... Procrustes: rmse 0.0001580537  max resid 0.0002995749 
    ## ... Similar to previous best
    ## Run 220 stress 0.0004561793 
    ## ... Procrustes: rmse 0.01547355  max resid 0.02131784 
    ## Run 221 stress 0.001321273 
    ## Run 222 stress 0.001229059 
    ## Run 223 stress 0.1491957 
    ## Run 224 stress 0.1491957 
    ## Run 225 stress 9.783281e-05 
    ## ... Procrustes: rmse 0.0001997954  max resid 0.0004870332 
    ## ... Similar to previous best
    ## Run 226 stress 8.318519e-05 
    ## ... Procrustes: rmse 0.0001509386  max resid 0.0002992045 
    ## ... Similar to previous best
    ## Run 227 stress 0.001317969 
    ## Run 228 stress 0.1990774 
    ## Run 229 stress 0.001409137 
    ## Run 230 stress 0.0004684557 
    ## ... Procrustes: rmse 0.01569412  max resid 0.02162052 
    ## Run 231 stress 9.526035e-05 
    ## ... Procrustes: rmse 0.0001767787  max resid 0.0003253968 
    ## ... Similar to previous best
    ## Run 232 stress 0.0007546476 
    ## Run 233 stress 0.00124377 
    ## Run 234 stress 0.1990774 
    ## Run 235 stress 9.085238e-05 
    ## ... Procrustes: rmse 9.622295e-05  max resid 0.0001914132 
    ## ... Similar to previous best
    ## Run 236 stress 0.1491957 
    ## Run 237 stress 9.278987e-05 
    ## ... Procrustes: rmse 0.0001829564  max resid 0.0004358366 
    ## ... Similar to previous best
    ## Run 238 stress 8.333966e-05 
    ## ... Procrustes: rmse 0.0002019347  max resid 0.0003598442 
    ## ... Similar to previous best
    ## Run 239 stress 0.0004610716 
    ## ... Procrustes: rmse 0.01554958  max resid 0.02141806 
    ## Run 240 stress 9.936419e-05 
    ## ... Procrustes: rmse 0.0001766913  max resid 0.0004107094 
    ## ... Similar to previous best
    ## Run 241 stress 8.475299e-05 
    ## ... Procrustes: rmse 0.0002024861  max resid 0.0003619216 
    ## ... Similar to previous best
    ## Run 242 stress 9.706017e-05 
    ## ... Procrustes: rmse 0.000384689  max resid 0.0005237094 
    ## ... Similar to previous best
    ## Run 243 stress 0.2509018 
    ## Run 244 stress 0.1776011 
    ## Run 245 stress 0.0004688072 
    ## ... Procrustes: rmse 0.01570393  max resid 0.02162315 
    ## Run 246 stress 0.001330022 
    ## Run 247 stress 9.292092e-05 
    ## ... Procrustes: rmse 0.0002389197  max resid 0.000491335 
    ## ... Similar to previous best
    ## Run 248 stress 0.001323314 
    ## Run 249 stress 9.221074e-05 
    ## ... Procrustes: rmse 0.0001948  max resid 0.0004089978 
    ## ... Similar to previous best
    ## Run 250 stress 0.0003981724 
    ## ... Procrustes: rmse 0.014458  max resid 0.01990795 
    ## Run 251 stress 0.00103874 
    ## Run 252 stress 9.154268e-05 
    ## ... Procrustes: rmse 0.0001696568  max resid 0.0003805926 
    ## ... Similar to previous best
    ## Run 253 stress 0.0006193287 
    ## Run 254 stress 0.2373687 
    ## Run 255 stress 0.2646973 
    ## Run 256 stress 0.1990774 
    ## Run 257 stress 0.001271635 
    ## Run 258 stress 0.001285742 
    ## Run 259 stress 9.266313e-05 
    ## ... Procrustes: rmse 0.0001630256  max resid 0.0003302164 
    ## ... Similar to previous best
    ## Run 260 stress 9.833637e-05 
    ## ... Procrustes: rmse 0.0002144298  max resid 0.0004474408 
    ## ... Similar to previous best
    ## Run 261 stress 9.76956e-05 
    ## ... Procrustes: rmse 0.0001988634  max resid 0.0004838656 
    ## ... Similar to previous best
    ## Run 262 stress 0.0004913634 
    ## ... Procrustes: rmse 0.01608459  max resid 0.02215324 
    ## Run 263 stress 8.310213e-05 
    ## ... Procrustes: rmse 0.001295216  max resid 0.001773316 
    ## ... Similar to previous best
    ## Run 264 stress 0.001229637 
    ## Run 265 stress 0.0009940021 
    ## Run 266 stress 9.602012e-05 
    ## ... Procrustes: rmse 0.006209236  max resid 0.00854124 
    ## Run 267 stress 0.1776019 
    ## Run 268 stress 0.3083098 
    ## Run 269 stress 0.0004791678 
    ## ... Procrustes: rmse 0.0264409  max resid 0.03687585 
    ## Run 270 stress 0.0004647744 
    ## ... Procrustes: rmse 0.01563229  max resid 0.02152212 
    ## Run 271 stress 0.1491957 
    ## Run 272 stress 0.0004771074 
    ## ... Procrustes: rmse 0.01584034  max resid 0.02181536 
    ## Run 273 stress 0.2848522 
    ## Run 274 stress 0.3083102 
    ## Run 275 stress 0.001255895 
    ## Run 276 stress 0.1491957 
    ## Run 277 stress 0.0008588891 
    ## Run 278 stress 0.1990774 
    ## Run 279 stress 0.0001273524 
    ## ... Procrustes: rmse 0.008100926  max resid 0.01113565 
    ## Run 280 stress 0.2520602 
    ## Run 281 stress 0.001365903 
    ## Run 282 stress 0.2440272 
    ## Run 283 stress 0.001267662 
    ## Run 284 stress 0.0001205133 
    ## ... Procrustes: rmse 0.007872597  max resid 0.0108204 
    ## Run 285 stress 0.000565941 
    ## Run 286 stress 7.966478e-05 
    ## ... Procrustes: rmse 0.0001874042  max resid 0.0002982055 
    ## ... Similar to previous best
    ## Run 287 stress 9.825532e-05 
    ## ... Procrustes: rmse 7.189568e-05  max resid 0.0001063916 
    ## ... Similar to previous best
    ## Run 288 stress 0.2440272 
    ## Run 289 stress 0.0004670386 
    ## ... Procrustes: rmse 0.01567753  max resid 0.02159057 
    ## Run 290 stress 8.122709e-05 
    ## ... Procrustes: rmse 0.0001398853  max resid 0.0002745759 
    ## ... Similar to previous best
    ## Run 291 stress 0.1990774 
    ## Run 292 stress 8.501623e-05 
    ## ... Procrustes: rmse 0.0001549721  max resid 0.0003300159 
    ## ... Similar to previous best
    ## Run 293 stress 0.001356128 
    ## Run 294 stress 0.2797329 
    ## Run 295 stress 9.160697e-05 
    ## ... Procrustes: rmse 0.001516116  max resid 0.002426345 
    ## ... Similar to previous best
    ## Run 296 stress 9.242087e-05 
    ## ... Procrustes: rmse 0.0002233646  max resid 0.0004900084 
    ## ... Similar to previous best
    ## Run 297 stress 9.858916e-05 
    ## ... Procrustes: rmse 0.0001754087  max resid 0.0003420704 
    ## ... Similar to previous best
    ## Run 298 stress 9.739389e-05 
    ## ... Procrustes: rmse 0.0001881885  max resid 0.0004592763 
    ## ... Similar to previous best
    ## Run 299 stress 9.174858e-05 
    ## ... Procrustes: rmse 0.0001591834  max resid 0.0003315136 
    ## ... Similar to previous best
    ## Run 300 stress 9.260252e-05 
    ## ... Procrustes: rmse 0.0002167069  max resid 0.0004618302 
    ## ... Similar to previous best
    ## Run 301 stress 9.562286e-05 
    ## ... Procrustes: rmse 0.0002415146  max resid 0.0004953612 
    ## ... Similar to previous best
    ## Run 302 stress 0.001170793 
    ## Run 303 stress 0.0004541976 
    ## ... Procrustes: rmse 0.0154578  max resid 0.02128731 
    ## Run 304 stress 9.296182e-05 
    ## ... Procrustes: rmse 0.0001510303  max resid 0.0002972009 
    ## ... Similar to previous best
    ## Run 305 stress 0.001275413 
    ## Run 306 stress 0.1776015 
    ## Run 307 stress 0.0004983557 
    ## ... Procrustes: rmse 0.01615979  max resid 0.02225651 
    ## Run 308 stress 0.001079911 
    ## Run 309 stress 8.589062e-05 
    ## ... Procrustes: rmse 0.004157651  max resid 0.006032907 
    ## ... Similar to previous best
    ## Run 310 stress 0.001301099 
    ## Run 311 stress 9.814459e-05 
    ## ... Procrustes: rmse 0.0001420928  max resid 0.000333541 
    ## ... Similar to previous best
    ## Run 312 stress 0.001205732 
    ## Run 313 stress 9.352097e-05 
    ## ... Procrustes: rmse 0.0001402786  max resid 0.0003280102 
    ## ... Similar to previous best
    ## Run 314 stress 0.1990774 
    ## Run 315 stress 0.001428633 
    ## Run 316 stress 0.001400111 
    ## Run 317 stress 8.923985e-05 
    ## ... Procrustes: rmse 0.0002182878  max resid 0.0004755526 
    ## ... Similar to previous best
    ## Run 318 stress 0.001277335 
    ## Run 319 stress 9.634759e-05 
    ## ... Procrustes: rmse 0.0002267371  max resid 0.0004927555 
    ## ... Similar to previous best
    ## Run 320 stress 9.061899e-05 
    ## ... Procrustes: rmse 0.000221018  max resid 0.0003868754 
    ## ... Similar to previous best
    ## Run 321 stress 0.0003294598 
    ## ... Procrustes: rmse 0.01313908  max resid 0.01808845 
    ## Run 322 stress 0.0003475584 
    ## ... Procrustes: rmse 0.01349967  max resid 0.01858572 
    ## Run 323 stress 9.464954e-05 
    ## ... Procrustes: rmse 0.0001775624  max resid 0.0003402596 
    ## ... Similar to previous best
    ## Run 324 stress 0.1491957 
    ## Run 325 stress 0.2696744 
    ## Run 326 stress 0.0009800226 
    ## Run 327 stress 0.001355693 
    ## Run 328 stress 9.759566e-05 
    ## ... Procrustes: rmse 0.0001709696  max resid 0.0004161353 
    ## ... Similar to previous best
    ## Run 329 stress 9.155794e-05 
    ## ... Procrustes: rmse 5.797756e-05  max resid 9.0459e-05 
    ## ... Similar to previous best
    ## Run 330 stress 9.680669e-05 
    ## ... Procrustes: rmse 0.005741183  max resid 0.00789958 
    ## Run 331 stress 9.446984e-05 
    ## ... Procrustes: rmse 0.0001636841  max resid 0.000369076 
    ## ... Similar to previous best
    ## Run 332 stress 0.3083098 
    ## Run 333 stress 9.766898e-05 
    ## ... Procrustes: rmse 0.0002292105  max resid 0.0004951541 
    ## ... Similar to previous best
    ## Run 334 stress 9.385319e-05 
    ## ... Procrustes: rmse 0.0001938281  max resid 0.0004694571 
    ## ... Similar to previous best
    ## Run 335 stress 9.554387e-05 
    ## ... Procrustes: rmse 0.0001524131  max resid 0.0002844557 
    ## ... Similar to previous best
    ## Run 336 stress 0.177601 
    ## Run 337 stress 0.2520602 
    ## Run 338 stress 0.0001722544 
    ## ... Procrustes: rmse 0.01584313  max resid 0.02221186 
    ## Run 339 stress 7.783913e-05 
    ## ... Procrustes: rmse 0.0002188439  max resid 0.0004538879 
    ## ... Similar to previous best
    ## Run 340 stress 9.731859e-05 
    ## ... Procrustes: rmse 0.0001894989  max resid 0.0004626246 
    ## ... Similar to previous best
    ## Run 341 stress 0.1990774 
    ## Run 342 stress 0.2509018 
    ## Run 343 stress 0.1491957 
    ## Run 344 stress 0.001329491 
    ## Run 345 stress 0.0003575833 
    ## ... Procrustes: rmse 0.01369401  max resid 0.0188544 
    ## Run 346 stress 0.001205357 
    ## Run 347 stress 0.0001372529 
    ## ... Procrustes: rmse 0.008416517  max resid 0.0115712 
    ## Run 348 stress 0.001265649 
    ## Run 349 stress 5.282491e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000750312  max resid 0.00135895 
    ## ... Similar to previous best
    ## Run 350 stress 0.001239827 
    ## Run 351 stress 9.602856e-05 
    ## ... Procrustes: rmse 0.0009193068  max resid 0.002108654 
    ## ... Similar to previous best
    ## Run 352 stress 0.00139978 
    ## Run 353 stress 0.001232695 
    ## Run 354 stress 9.181048e-05 
    ## ... Procrustes: rmse 0.0006923627  max resid 0.001086802 
    ## ... Similar to previous best
    ## Run 355 stress 9.104515e-05 
    ## ... Procrustes: rmse 0.0006904646  max resid 0.001084045 
    ## ... Similar to previous best
    ## Run 356 stress 9.338105e-05 
    ## ... Procrustes: rmse 0.0007192188  max resid 0.001104799 
    ## ... Similar to previous best
    ## Run 357 stress 0.001199435 
    ## Run 358 stress 0.1990774 
    ## Run 359 stress 0.1491957 
    ## Run 360 stress 9.431948e-05 
    ## ... Procrustes: rmse 0.0007337682  max resid 0.001142486 
    ## ... Similar to previous best
    ## Run 361 stress 9.119339e-05 
    ## ... Procrustes: rmse 0.0006910012  max resid 0.001145683 
    ## ... Similar to previous best
    ## Run 362 stress 8.680639e-05 
    ## ... Procrustes: rmse 0.0006881805  max resid 0.001099643 
    ## ... Similar to previous best
    ## Run 363 stress 9.589128e-05 
    ## ... Procrustes: rmse 0.0007313821  max resid 0.001097407 
    ## ... Similar to previous best
    ## Run 364 stress 0.177601 
    ## Run 365 stress 8.884422e-05 
    ## ... Procrustes: rmse 0.002830711  max resid 0.005132457 
    ## ... Similar to previous best
    ## Run 366 stress 0.0004868159 
    ## ... Procrustes: rmse 0.01600921  max resid 0.02340345 
    ## Run 367 stress 9.871016e-05 
    ## ... Procrustes: rmse 0.0007512537  max resid 0.001079574 
    ## ... Similar to previous best
    ## Run 368 stress 0.001137683 
    ## Run 369 stress 8.351607e-05 
    ## ... Procrustes: rmse 0.0005453074  max resid 0.000707513 
    ## ... Similar to previous best
    ## Run 370 stress 0.3083098 
    ## Run 371 stress 0.001207546 
    ## Run 372 stress 0.1491957 
    ## Run 373 stress 0.001270476 
    ## Run 374 stress 0.001304929 
    ## Run 375 stress 9.882248e-05 
    ## ... Procrustes: rmse 0.0008173863  max resid 0.001802093 
    ## ... Similar to previous best
    ## Run 376 stress 0.0004664245 
    ## ... Procrustes: rmse 0.01562817  max resid 0.02286425 
    ## Run 377 stress 0.1491957 
    ## Run 378 stress 9.22979e-05 
    ## ... Procrustes: rmse 0.0007218095  max resid 0.001121835 
    ## ... Similar to previous best
    ## Run 379 stress 0.001234498 
    ## Run 380 stress 9.211688e-05 
    ## ... Procrustes: rmse 0.0006159577  max resid 0.000919823 
    ## ... Similar to previous best
    ## Run 381 stress 0.0002954514 
    ## ... Procrustes: rmse 0.01243667  max resid 0.01846874 
    ## Run 382 stress 9.166033e-05 
    ## ... Procrustes: rmse 0.0006682483  max resid 0.0008891132 
    ## ... Similar to previous best
    ## Run 383 stress 9.09773e-05 
    ## ... Procrustes: rmse 0.0006950149  max resid 0.001095379 
    ## ... Similar to previous best
    ## Run 384 stress 0.0008210986 
    ## Run 385 stress 0.001108998 
    ## Run 386 stress 9.141978e-05 
    ## ... Procrustes: rmse 0.0007107418  max resid 0.0009908665 
    ## ... Similar to previous best
    ## Run 387 stress 9.904596e-05 
    ## ... Procrustes: rmse 0.0007368677  max resid 0.001079668 
    ## ... Similar to previous best
    ## Run 388 stress 8.978735e-05 
    ## ... Procrustes: rmse 0.0006990779  max resid 0.001141551 
    ## ... Similar to previous best
    ## Run 389 stress 0.1491957 
    ## Run 390 stress 0.0004632374 
    ## ... Procrustes: rmse 0.01561277  max resid 0.02285559 
    ## Run 391 stress 0.1491957 
    ## Run 392 stress 0.001496206 
    ## Run 393 stress 0.0004949009 
    ## ... Procrustes: rmse 0.01610095  max resid 0.0235409 
    ## Run 394 stress 0.000967622 
    ## Run 395 stress 0.001179721 
    ## Run 396 stress 0.001206091 
    ## Run 397 stress 9.998901e-05 
    ## ... Procrustes: rmse 0.007176975  max resid 0.01119047 
    ## Run 398 stress 9.694928e-05 
    ## ... Procrustes: rmse 0.0006980216  max resid 0.001095922 
    ## ... Similar to previous best
    ## Run 399 stress 0.177601 
    ## Run 400 stress 0.1491957 
    ## Run 401 stress 9.941083e-05 
    ## ... Procrustes: rmse 0.0006922525  max resid 0.001141839 
    ## ... Similar to previous best
    ## Run 402 stress 9.942192e-05 
    ## ... Procrustes: rmse 0.0006889713  max resid 0.001138937 
    ## ... Similar to previous best
    ## Run 403 stress 0.0003226528 
    ## ... Procrustes: rmse 0.01299989  max resid 0.01924718 
    ## Run 404 stress 9.365594e-05 
    ## ... Procrustes: rmse 0.000717332  max resid 0.001056904 
    ## ... Similar to previous best
    ## Run 405 stress 9.854464e-05 
    ## ... Procrustes: rmse 0.0006099641  max resid 0.0009189153 
    ## ... Similar to previous best
    ## Run 406 stress 0.000503438 
    ## ... Procrustes: rmse 0.01627798  max resid 0.02377542 
    ## Run 407 stress 0.001155437 
    ## Run 408 stress 0.1491957 
    ## Run 409 stress 0.001074357 
    ## Run 410 stress 9.419034e-05 
    ## ... Procrustes: rmse 0.000683704  max resid 0.001106414 
    ## ... Similar to previous best
    ## Run 411 stress 0.1491957 
    ## Run 412 stress 0.00119059 
    ## Run 413 stress 7.256668e-05 
    ## ... Procrustes: rmse 0.0008387939  max resid 0.001903833 
    ## ... Similar to previous best
    ## Run 414 stress 0.001363547 
    ## Run 415 stress 0.001398084 
    ## Run 416 stress 0.00138534 
    ## Run 417 stress 0.000343348 
    ## ... Procrustes: rmse 0.02165839  max resid 0.02989167 
    ## Run 418 stress 0.001251147 
    ## Run 419 stress 9.555348e-05 
    ## ... Procrustes: rmse 0.0007081708  max resid 0.0009956797 
    ## ... Similar to previous best
    ## Run 420 stress 0.001433226 
    ## Run 421 stress 8.872132e-05 
    ## ... Procrustes: rmse 0.0006945196  max resid 0.001082485 
    ## ... Similar to previous best
    ## Run 422 stress 9.105006e-05 
    ## ... Procrustes: rmse 0.0006686728  max resid 0.0009406325 
    ## ... Similar to previous best
    ## Run 423 stress 9.505108e-05 
    ## ... Procrustes: rmse 0.004177081  max resid 0.005777095 
    ## ... Similar to previous best
    ## Run 424 stress 0.1776012 
    ## Run 425 stress 8.261643e-05 
    ## ... Procrustes: rmse 0.0007524396  max resid 0.001053791 
    ## ... Similar to previous best
    ## Run 426 stress 9.897496e-05 
    ## ... Procrustes: rmse 0.0007071459  max resid 0.001131326 
    ## ... Similar to previous best
    ## Run 427 stress 0.1776012 
    ## Run 428 stress 0.0002204341 
    ## ... Procrustes: rmse 0.01072215  max resid 0.01609847 
    ## Run 429 stress 0.0004900598 
    ## ... Procrustes: rmse 0.01605254  max resid 0.02346344 
    ## Run 430 stress 0.001195026 
    ## Run 431 stress 9.720611e-05 
    ## ... Procrustes: rmse 0.0006103585  max resid 0.0009236427 
    ## ... Similar to previous best
    ## Run 432 stress 9.672166e-05 
    ## ... Procrustes: rmse 0.002251144  max resid 0.004296237 
    ## ... Similar to previous best
    ## Run 433 stress 0.001316236 
    ## Run 434 stress 0.001471075 
    ## Run 435 stress 0.001310231 
    ## Run 436 stress 0.2354286 
    ## Run 437 stress 0.2852177 
    ## Run 438 stress 0.2848516 
    ## Run 439 stress 0.1491957 
    ## Run 440 stress 0.0004923454 
    ## ... Procrustes: rmse 0.01610113  max resid 0.02353249 
    ## Run 441 stress 0.0004741785 
    ## ... Procrustes: rmse 0.01579726  max resid 0.02311092 
    ## Run 442 stress 0.0004948733 
    ## ... Procrustes: rmse 0.01614302  max resid 0.02358798 
    ## Run 443 stress 0.001301796 
    ## Run 444 stress 0.1990774 
    ## Run 445 stress 0.001288201 
    ## Run 446 stress 0.1491957 
    ## Run 447 stress 0.001118917 
    ## Run 448 stress 7.861207e-05 
    ## ... Procrustes: rmse 0.0007301112  max resid 0.001092216 
    ## ... Similar to previous best
    ## Run 449 stress 0.0009768964 
    ## Run 450 stress 0.001299494 
    ## Run 451 stress 0.001437495 
    ## Run 452 stress 7.812569e-05 
    ## ... Procrustes: rmse 0.00069903  max resid 0.001062372 
    ## ... Similar to previous best
    ## Run 453 stress 0.1491957 
    ## Run 454 stress 8.638612e-05 
    ## ... Procrustes: rmse 0.0007171371  max resid 0.001068443 
    ## ... Similar to previous best
    ## Run 455 stress 8.861774e-05 
    ## ... Procrustes: rmse 0.003704065  max resid 0.005112893 
    ## ... Similar to previous best
    ## Run 456 stress 9.131849e-05 
    ## ... Procrustes: rmse 0.00069046  max resid 0.001119137 
    ## ... Similar to previous best
    ## Run 457 stress 0.001269286 
    ## Run 458 stress 0.0004128666 
    ## ... Procrustes: rmse 0.02382341  max resid 0.03288879 
    ## Run 459 stress 0.001350351 
    ## Run 460 stress 0.001283238 
    ## Run 461 stress 9.242834e-05 
    ## ... Procrustes: rmse 0.0006815968  max resid 0.001129745 
    ## ... Similar to previous best
    ## Run 462 stress 0.001152708 
    ## Run 463 stress 0.00120163 
    ## Run 464 stress 0.001238082 
    ## Run 465 stress 0.1491957 
    ## Run 466 stress 0.001342126 
    ## Run 467 stress 0.001360495 
    ## Run 468 stress 0.1990774 
    ## Run 469 stress 9.896179e-05 
    ## ... Procrustes: rmse 0.0007132652  max resid 0.0009904231 
    ## ... Similar to previous best
    ## Run 470 stress 0.2520602 
    ## Run 471 stress 0.3082973 
    ## Run 472 stress 0.001358938 
    ## Run 473 stress 8.517406e-05 
    ## ... Procrustes: rmse 0.0007480894  max resid 0.001083031 
    ## ... Similar to previous best
    ## Run 474 stress 0.001183785 
    ## Run 475 stress 0.1491957 
    ## Run 476 stress 9.76973e-05 
    ## ... Procrustes: rmse 0.0006452184  max resid 0.000959437 
    ## ... Similar to previous best
    ## Run 477 stress 9.347255e-05 
    ## ... Procrustes: rmse 0.0006511098  max resid 0.001017009 
    ## ... Similar to previous best
    ## Run 478 stress 9.399936e-05 
    ## ... Procrustes: rmse 0.000613933  max resid 0.0009194022 
    ## ... Similar to previous best
    ## Run 479 stress 9.868626e-05 
    ## ... Procrustes: rmse 0.0006800745  max resid 0.001119308 
    ## ... Similar to previous best
    ## Run 480 stress 9.215531e-05 
    ## ... Procrustes: rmse 0.0006552598  max resid 0.001014231 
    ## ... Similar to previous best
    ## Run 481 stress 0.001373986 
    ## Run 482 stress 8.953672e-05 
    ## ... Procrustes: rmse 0.0006917637  max resid 0.001100177 
    ## ... Similar to previous best
    ## Run 483 stress 0.3083098 
    ## Run 484 stress 0.0006398816 
    ## Run 485 stress 9.279897e-05 
    ## ... Procrustes: rmse 0.000697336  max resid 0.001104539 
    ## ... Similar to previous best
    ## Run 486 stress 0.001232884 
    ## Run 487 stress 9.354037e-05 
    ## ... Procrustes: rmse 0.000697343  max resid 0.001103198 
    ## ... Similar to previous best
    ## Run 488 stress 0.0004929687 
    ## ... Procrustes: rmse 0.01611174  max resid 0.0235466 
    ## Run 489 stress 7.302428e-05 
    ## ... Procrustes: rmse 0.0004245703  max resid 0.0005492416 
    ## ... Similar to previous best
    ## Run 490 stress 0.001218454 
    ## Run 491 stress 0.001197503 
    ## Run 492 stress 0.0001609757 
    ## ... Procrustes: rmse 0.01459459  max resid 0.0201233 
    ## Run 493 stress 0.001324079 
    ## Run 494 stress 6.76976e-05 
    ## ... Procrustes: rmse 0.00289575  max resid 0.005226854 
    ## ... Similar to previous best
    ## Run 495 stress 9.421513e-05 
    ## ... Procrustes: rmse 0.0007333526  max resid 0.001136844 
    ## ... Similar to previous best
    ## Run 496 stress 0.1491957 
    ## Run 497 stress 0.001314259 
    ## Run 498 stress 0.001279088 
    ## Run 499 stress 9.07211e-05 
    ## ... Procrustes: rmse 0.0006853701  max resid 0.001108064 
    ## ... Similar to previous best
    ## Run 500 stress 0.1776019 
    ## *** Best solution repeated 55 times

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
    ##                                    NMDS1     NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)]  0.997860 -0.065421 0.7388  0.038 *
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
    ##                                    NMDS1     NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)]  0.997860 -0.065421 0.7388  0.038 *
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
    ## temperature_median -0.48784 -0.87293 0.0586  0.816  
    ## salinity_median     0.99786 -0.06542 0.7388  0.031 *
    ## oxygen_median       0.79330  0.60883 0.5361  0.848  
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
    ## temperature_median -0.48784 -0.87293 0.0586  1.000  
    ## salinity_median     0.99786 -0.06542 0.7388  0.093 .
    ## oxygen_median       0.79330  0.60883 0.5361  1.000  
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
    ## temperature_median -0.56523 -0.82493 0.0778  0.789  
    ## salinity_median     0.99174  0.12828 0.7458  0.033 *
    ## oxygen_median       0.94263  0.33383 0.3552  0.899  
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
    ## temperature_median -0.56523 -0.82493 0.0778  1.000  
    ## salinity_median     0.99174  0.12828 0.7458  0.099 .
    ## oxygen_median       0.94263  0.33383 0.3552  1.000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ## temperature_median -0.10614  0.99435 0.0705  0.670  
    ## salinity_median    -0.78731 -0.61655 0.7842  0.013 *
    ## oxygen_median      -0.96711  0.25435 0.7188  0.043 *
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
    ## salinity_median    -0.78731 -0.61655 0.7842  0.039 *
    ## oxygen_median      -0.96711  0.25435 0.7188  0.129  
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
    ## temperature_median  0.00010397  1.00000000 0.1172  0.384
    ## salinity_median    -0.00061311 -1.00000000 0.5232  0.445
    ## oxygen_median      -0.00125668  1.00000000 0.7255  0.510
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
    ## temperature_median  0.00010397  1.00000000 0.1172      1
    ## salinity_median    -0.00061311 -1.00000000 0.5232      1
    ## oxygen_median      -0.00125668  1.00000000 0.7255      1
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
    ## temperature_median -0.00006954 -1.00000000 0.2313  0.482  
    ## salinity_median     0.00048733  1.00000000 0.4599  0.211  
    ## oxygen_median       0.00072584  1.00000000 0.6373  0.083 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ## temperature_median -0.00006954 -1.00000000 0.2313  1.000
    ## salinity_median     0.00048733  1.00000000 0.4599  0.633
    ## oxygen_median       0.00072584  1.00000000 0.6373  0.249
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
    ## temperature_median       -0.35881 -0.93341 0.1964  0.597
    ## salinity_median          -0.70036  0.71379 0.5991  0.122
    ## oxygen_median             0.98543 -0.17006 0.5185  0.168
    ## distance_to_ocean_mean_m  0.01326 -0.99991 0.1726  0.564
    ## max_depth                -0.49998 -0.86604 0.0383  0.889
    ## logArea                  -0.26513 -0.96421 0.2144  0.542
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
    ## temperature_median       -0.35881 -0.93341 0.1964  1.000
    ## salinity_median          -0.70036  0.71379 0.5991  0.732
    ## oxygen_median             0.98543 -0.17006 0.5185  1.000
    ## distance_to_ocean_mean_m  0.01326 -0.99991 0.1726  1.000
    ## max_depth                -0.49998 -0.86604 0.0383  1.000
    ## logArea                  -0.26513 -0.96421 0.2144  1.000
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
# For figure
(SD_beta_geo_ef <- envfit(SD_beta_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.99979  0.02049 0.6139  0.152
    ## max_depth               -0.17727 -0.98416 0.1185  0.772
    ## logArea                  0.26565 -0.96407 0.2080  0.117
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
    ## distance_to_ocean_min_m -0.99979  0.02049 0.6139  0.456
    ## max_depth               -0.17727 -0.98416 0.1185  1.000
    ## logArea                  0.26565 -0.96407 0.2080  0.351
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by stratification
(SD_beta_geo_A_ef <- envfit(SD_beta_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.99979  0.02049 0.6139  0.143
    ## max_depth               -0.17727 -0.98416 0.1185  0.786
    ## logArea                  0.26565 -0.96407 0.2080  0.117
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
    ## distance_to_ocean_min_m -0.99979  0.02049 0.6139  0.429
    ## max_depth               -0.17727 -0.98416 0.1185  1.000
    ## logArea                  0.26565 -0.96407 0.2080  0.351
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
(SD_beta_geo_MS_ef <- envfit(SD_beta_geo_MS_NMDS, env[mixed_stratified_lakes,geography], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.94035  0.34020 0.5170  0.137
    ## max_depth               -0.21480 -0.97666 0.1851  0.636
    ## logArea                 -0.06544  0.99786 0.0118  0.960
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
    ## distance_to_ocean_min_m -0.94035  0.34020 0.5170  0.411
    ## max_depth               -0.21480 -0.97666 0.1851  1.000
    ## logArea                 -0.06544  0.99786 0.0118  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
(SD_beta_geo_OM_ef <- envfit(SD_beta_geo_OM_NMDS, env[ocean_mixed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m  0.86178 -0.50728 0.3198  0.075 .
    ## max_depth               -0.68923 -0.72454 0.5689  0.015 *
    ## logArea                 -0.84909  0.52824 0.2993  0.195  
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
    ## distance_to_ocean_min_m  0.86178 -0.50728 0.3198  0.225  
    ## max_depth               -0.68923 -0.72454 0.5689  0.045 *
    ## logArea                 -0.84909  0.52824 0.2993  0.585  
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
    ## distance_to_ocean_min_m  0.86963  0.49370 0.6793  0.179
    ## max_depth                0.68537  0.72820 0.0510  0.962
    ## logArea                 -0.52042  0.85391 0.2187  0.510
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
    ## distance_to_ocean_min_m  0.86963  0.49370 0.6793  0.537
    ## max_depth                0.68537  0.72820 0.0510  1.000
    ## logArea                 -0.52042  0.85391 0.2187  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed lakes
(SD_beta_geo_M_ef <- envfit(SD_beta_geo_M_NMDS, env[mixed_lakes,geography], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                              NMDS1      NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m -0.0019666  1.0000000 0.6727  0.059 . 
    ## max_depth                0.0047648 -0.9999900 0.8424  0.012 * 
    ## logArea                  0.0038670  0.9999900 0.8464  0.003 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_geo_M_efp <- p.adjust.envfit(SD_beta_geo_M_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                              NMDS1      NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m -0.0019666  1.0000000 0.6727  0.177   
    ## max_depth                0.0047648 -0.9999900 0.8424  0.036 * 
    ## logArea                  0.0038670  0.9999900 0.8464  0.009 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

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
    ##       Significance: 0.248 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.296 0.323 0.352 0.365 
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
    ## 0.574 0.605 0.628 0.649 
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
    ##       Significance: 0.468 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.587 0.626 0.665 0.701 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_mant_pv <- rbind(SD_beta_env_mant_t$signif, SD_beta_env_mant_s$signif, SD_beta_env_mant_o$signif)
SD_beta_env_mant_pv <- SD_beta_env_mant_pv[,1]
(SD_beta_env_mant_pv <- p.adjust(SD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.744 0.018 1.000

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
    ##       Significance: 0.321 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.328 0.361 0.381 0.402 
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
    ##       Significance: 0.007 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.546 0.588 0.630 0.666 
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
    ##       Significance: 0.589 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.475 0.525 0.564 0.615 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_MS_mant_pv <- rbind(SD_beta_env_MS_mant_t$signif, SD_beta_env_MS_mant_s$signif, SD_beta_env_MS_mant_o$signif)
SD_beta_env_MS_mant_pv <- SD_beta_env_MS_mant_pv[,1]
(SD_beta_env_MS_mant_pv <- p.adjust(SD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.963 0.021 1.000

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
    ##       Significance: 0.755 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.118 0.156 0.195 0.372 
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
    ##       Significance: 0.062 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.367 0.441 0.493 0.549 
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
    ##       Significance: 0.038 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.342 0.474 0.598 0.671 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_OM_mant_pv <- rbind(SD_beta_env_OM_mant_t$signif, SD_beta_env_OM_mant_s$signif, SD_beta_env_OM_mant_o$signif)
SD_beta_env_OM_mant_pv <- SD_beta_env_OM_mant_pv[,1]
(SD_beta_env_OM_mant_pv <- p.adjust(SD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.186 0.114

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
    ##       Significance: 0.218 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0282 0.0532 0.0699 0.0828 
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
    ##       Significance: 0.008 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.454 0.482 0.517 0.547 
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
    ##       Significance: 0.382 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.713 0.735 0.755 0.769 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_SO_mant_pv <- rbind(SD_beta_env_SO_mant_t$signif, SD_beta_env_SO_mant_s$signif, SD_beta_env_SO_mant_o$signif)
SD_beta_env_SO_mant_pv <- SD_beta_env_SO_mant_pv[,1]
(SD_beta_env_SO_mant_pv <- p.adjust(SD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 0.654 0.024 1.000

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
    ##       Significance: 0.756 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.259 0.337 0.459 0.716 
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
    ##       Significance: 0.153 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.262 0.320 0.373 0.426 
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
    ##       Significance: 0.084 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.227 0.376 0.472 0.545 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_M_mant_pv <- rbind(SD_beta_env_M_mant_t$signif, SD_beta_env_M_mant_s$signif, SD_beta_env_M_mant_o$signif)
SD_beta_env_M_mant_pv <- SD_beta_env_M_mant_pv[,1]
(SD_beta_env_M_mant_pv <- p.adjust(SD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.459 0.252

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
    ##       Significance: 0.242 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.217 0.304 0.368 0.420 
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
geo_dist_dmean <- dist(scaled_env[surveyed_sites,c(9)], method = "euclidean")
(SD_beta_geo_mant_dmean <- mantel(SD_beta_geo_dist$Btotal, geo_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_dist$Btotal, ydis = geo_dist_dmean,      method = "spearman", permutations = 999, strata = env[surveyed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4392 
    ##       Significance: 0.448 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.541 0.569 0.590 0.613 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_md <- dist(scaled_env[surveyed_sites,c(14)], method = "euclidean")
(SD_beta_geo_mant_md <- mantel(SD_beta_geo_dist$Btotal, geo_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_dist$Btotal, ydis = geo_dist_md, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.06903 
    ##       Significance: 0.628 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.191 0.225 0.247 0.280 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_la <- dist(scaled_env[surveyed_sites,c(15)], method = "euclidean")
(SD_beta_geo_mant_la <- mantel(SD_beta_geo_dist$Btotal, geo_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_dist$Btotal, ydis = geo_dist_la, method = "spearman",      permutations = 999, strata = env[surveyed_sites, 19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.04299 
    ##       Significance: 0.412 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0818 0.0959 0.1053 0.1140 
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
geo_MS_dist_dmean <- dist(scaled_env[mixed_stratified_lakes,c(9)], method = "euclidean")
(SD_beta_geo_MS_mant_dmean <- mantel(SD_beta_geo_MS_dist$Btotal, geo_MS_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_dmean,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.299 
    ##       Significance: 0.418 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.435 0.471 0.508 0.546 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_md <- dist(scaled_env[mixed_stratified_lakes,c(14)], method = "euclidean")
(SD_beta_geo_MS_mant_md <- mantel(SD_beta_geo_MS_dist$Btotal, geo_MS_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_md,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.007317 
    ##       Significance: 0.732 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.209 0.256 0.295 0.315 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_la <- dist(scaled_env[mixed_stratified_lakes,c(15)], method = "euclidean")
(SD_beta_geo_MS_mant_la <- mantel(SD_beta_geo_MS_dist$Btotal, geo_MS_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_la,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.04736 
    ##       Significance: 0.76 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0665 0.0897 0.1054 0.1286 
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
geo_OM_dist_dmean <- dist(scaled_env[ocean_mixed_sites,c(9)], method = "euclidean")
(SD_beta_geo_OM_mant_dmean <- mantel(SD_beta_geo_OM_dist$Btotal, geo_OM_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_dmean,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.08705 
    ##       Significance: 0.583 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.238 0.272 0.292 0.334 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_md <- dist(scaled_env[ocean_mixed_sites,c(14)], method = "euclidean")
(SD_beta_geo_OM_mant_md <- mantel(SD_beta_geo_OM_dist$Btotal, geo_OM_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_md,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2349 
    ##       Significance: 0.063 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.183 0.260 0.311 0.374 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_la <- dist(scaled_env[ocean_mixed_sites,c(15)], method = "euclidean")
(SD_beta_geo_OM_mant_la <- mantel(SD_beta_geo_OM_dist$Btotal, geo_OM_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,19]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_la,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          19], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2007 
    ##       Significance: 0.107 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.205 0.242 0.274 0.302 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_OM_mant_pv <- rbind( SD_beta_geo_OM_mant_dmean$signif, SD_beta_geo_OM_mant_md$signif, SD_beta_geo_OM_mant_la$signif)
SD_beta_geo_OM_mant_pv <- SD_beta_geo_OM_mant_pv[,1]
(SD_beta_geo_OM_mant_pv <- p.adjust(SD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.189 0.321

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
    ##       Significance: 0.482 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.679 0.708 0.729 0.751 
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
    ##       Significance: 0.642 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.135 0.168 0.192 0.244 
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
    ##       Significance: 0.231 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.277 0.293 0.306 0.323 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_SO_mant_pv <- rbind( SD_beta_geo_SO_mant_dmean$signif, SD_beta_geo_SO_mant_md$signif, SD_beta_geo_SO_mant_la$signif)
SD_beta_geo_SO_mant_pv <- SD_beta_geo_SO_mant_pv[,1]
(SD_beta_geo_SO_mant_pv <- p.adjust(SD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.693

``` r
# Mixed lakes
geo_M_dist_dmean <- dist(scaled_env[mixed_lakes,c(9)], method = "euclidean")
(SD_beta_geo_M_mant_dmean <- mantel(SD_beta_geo_M_dist$Btotal, geo_M_dist_dmean, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_dmean,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.02135 
    ##       Significance: 0.479 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.256 0.383 0.462 0.707 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_M_dist_md <- dist(scaled_env[mixed_lakes,c(14)], method = "euclidean")
(SD_beta_geo_M_mant_md <- mantel(SD_beta_geo_M_dist$Btotal, geo_M_dist_md, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_md,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6829 
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.275 0.414 0.524 0.579 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_M_dist_la <- dist(scaled_env[mixed_lakes,c(15)], method = "euclidean")
(SD_beta_geo_M_mant_la <- mantel(SD_beta_geo_M_dist$Btotal, geo_M_dist_la, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_M_dist$Btotal, ydis = geo_M_dist_la,      method = "spearman", permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.416 
    ##       Significance: 0.027 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.227 0.313 0.419 0.502 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_M_mant_pv <- rbind( SD_beta_geo_M_mant_dmean$signif, SD_beta_geo_M_mant_md$signif, SD_beta_geo_M_mant_la$signif)
SD_beta_geo_M_mant_pv <- SD_beta_geo_M_mant_pv[,1]
(SD_beta_geo_M_mant_pv <- p.adjust(SD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.006 0.081

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
